#!/usr/bin/env python3
"""fuzz_small.py - randomized differential tests on small adversarial inputs.

Each round builds a random multi-record genome from the shapes that broke optik before, namely reverse-complement palindromes (a k-mer and its reverse complement in one super-k-mer), exact and reverse-complement repeats, homopolymers, tandem repeats, N runs, lowercase and records shorter than k, at a random odd or even k from 3 to 32.
Each check compares with a pure-Python oracle or an invariance, never with optik itself.

  fm / spss       query counts (positive, negative, invalid) equal the truth
                  at -c 1 (build-fm, build-spss) and at -c 2 (build-fm -c,
                  build-spss --compact -c)
  spss --ids      one id per distinct canonical k-mer, the same on both strands, and no id shared by two k-mers
  dbruijn         the GFA segments partition the solid k-mer set exactly (-c 1 and -c 2), KC:i equals the summed abundances, and the GFA is byte-identical at -t 1 and -t 4
  spss-update     (base | added) - removed equals the truth, with one id per k-mer afterwards
  fastq.gz        a gzipped FASTQ copy of the genome gives the same answers
  map             one primary record per read (short reads and reads of 500-3000 bp with 5% edits), with identical SAM at -t 1 and -t 4

Usage: python3 bench/fuzz_small.py [--optik target/release/optik] [--rounds 80] [--seed 1]
"""

import argparse
import enum
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from diagnostics import DETAIL_CHARS
from dna import revcomp
from sam_spec import (
    DBG_MAX_K,
    SAM_TAG_PREFIX_LEN,
    GfaSegmentCol,
    SamFlag,
)

HERE = os.path.dirname(os.path.abspath(__file__))


class Shape(enum.IntEnum):
    """Adversarial sequence shapes a fuzz genome is built from, each drawn with equal chance.

    The values are the draws of rng.randrange(len(Shape)), so reordering them changes every seeded genome.
    """

    RANDOM = 0
    PALINDROME = 1
    HOMOPOLYMER = 2
    TANDEM = 3
    N_RUN = 4
    REPEAT = 5
    """A copy of an earlier part of the record."""
    RC_REPEAT = 6
    """A reverse-complemented copy of an earlier part of the record."""
    LOWERCASE = 7


class Edit(enum.IntEnum):
    """Edits applied to a long read, each drawn with equal chance.

    As with Shape, the values are the rng draws themselves.
    """

    SUBSTITUTION = 0
    INSERTION = 1
    DELETION = 2


# A genome is at most MAX_RECORDS records of at most MAX_PARTS parts, each part one Shape.
MAX_RECORDS = 6
MAX_PARTS = 12
MAX_RANDOM_LEN = 400
MAX_PALINDROME_HALF = 30
MAX_HOMOPOLYMER_LEN = 200
MAX_TANDEM_UNIT = 7
MIN_TANDEM_COPIES = 2
MAX_TANDEM_COPIES = 60
MAX_N_RUN = 20
MAX_LOWERCASE_LEN = 60
# k values a round draws from: odd and even, around the 2-bit word boundaries
KMER_CHOICES = [3, 4, 5, 8, 11, 15, 16, 20, 21, 25, 31, 32]
# a round's seed is seed x ROUND_SEED_MULT + round, a prime so seeds and rounds never collide
ROUND_SEED_MULT = 100003
# Each record contributes QUERIES_PER_RECORD query windows of up to QUERY_EXTRA bases past k.
# A window may be reverse-complemented, and with QUERY_MUTATE_PROB it also gets one random base.
# RANDOM_QUERIES fully random queries, up to RANDOM_QUERY_EXTRA bases past k, are added on top.
QUERIES_PER_RECORD = 20
QUERY_EXTRA = 40
QUERY_MUTATE_PROB = 0.3
RANDOM_QUERIES = 50
RANDOM_QUERY_EXTRA = 30
# spss-update removes every REMOVE_STRIDE-th query
REMOVE_STRIDE = 3
# short reads: records of at least SHORT_READ_MIN_REC bases give SHORT_READS_PER_REC reads of SHORT_READ_MIN_LEN..SHORT_READ_MAX_LEN bases with up to MAX_READ_SUBS random bases
SHORT_READ_MIN_REC = 60
SHORT_READS_PER_REC = 15
SHORT_READ_MIN_LEN = 30
SHORT_READ_MAX_LEN = 300
MAX_READ_SUBS = 6
REVCOMP_PROB = 0.5
# long reads (the chaining path): records of at least LONG_READ_MIN_REC bases give LONG_READS_PER_REC reads of LONG_READ_MIN_LEN..LONG_READ_MAX_LEN bases, one edit per LONG_READ_EDIT_SPACING bases
LONG_READ_MIN_REC = 600
LONG_READS_PER_REC = 3
LONG_READ_MIN_LEN = 500
LONG_READ_MAX_LEN = 3000
LONG_READ_EDIT_SPACING = 20
# harness
DEFAULT_ROUNDS = 80
DEFAULT_SEED = 1
RUN_TIMEOUT_S = 600


def kmers(seq, kmer_len):
    """Each k-mer window of s, upper-cased, or None for a window holding a non-ACGT base."""
    seq = seq.upper()
    for i in range(len(seq) - kmer_len + 1):
        window = seq[i : i + kmer_len]
        yield window if all(base in "ACGT" for base in window) else None


def genome(rng):
    """Random genome of one to six records built from the adversarial shapes."""
    recs = []
    for i in range(rng.randint(1, MAX_RECORDS)):
        parts = []
        for _ in range(rng.randint(1, MAX_PARTS)):
            kind = Shape(rng.randrange(len(Shape)))
            if kind == Shape.RANDOM:
                parts.append("".join(rng.choice("ACGT") for _ in range(rng.randint(1, MAX_RANDOM_LEN))))
            elif kind == Shape.PALINDROME:
                half = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, MAX_PALINDROME_HALF)))
                parts.append(half + revcomp(half))
            elif kind == Shape.HOMOPOLYMER:
                parts.append(rng.choice("ACGT") * rng.randint(1, MAX_HOMOPOLYMER_LEN))
            elif kind == Shape.TANDEM:
                unit = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, MAX_TANDEM_UNIT)))
                parts.append(unit * rng.randint(MIN_TANDEM_COPIES, MAX_TANDEM_COPIES))
            elif kind == Shape.N_RUN:
                parts.append("N" * rng.randint(1, MAX_N_RUN))
            elif kind == Shape.REPEAT and parts:
                parts.append(parts[rng.randrange(len(parts))])
            elif kind == Shape.RC_REPEAT and parts:
                parts.append(revcomp(parts[rng.randrange(len(parts))]))
            else:  # Shape.LOWERCASE, or a repeat drawn before the record has a part to copy
                parts.append(
                    "".join(rng.choice("ACGT") for _ in range(rng.randint(1, MAX_LOWERCASE_LEN))).lower()
                )
        recs.append((f"r{i}", "".join(parts)))
    return recs


def main():
    """Run every fuzz round and return 1 if any check failed, else 0."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    parser.add_argument("--rounds", type=int, default=DEFAULT_ROUNDS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    optik = args.optik
    fails = []
    checks = 0

    def run(cmd):
        return subprocess.run(
            cmd, capture_output=True, text=True, timeout=RUN_TIMEOUT_S, check=False
        )

    def report_value(out, key):
        """Integer value of the report line starting with `key`."""
        for line in out.splitlines():
            if line.startswith(key):
                return int(line.split("=")[1].split()[0])
        return None

    with tempfile.TemporaryDirectory(prefix="optik_fuzz_") as tmp_dir:

        def in_dir(name):
            """Path of f in this round's directory."""
            return os.path.join(tmp_dir, name)

        def write(path, recs):
            """Write records as single-line FASTA."""
            with open(path, "w") as handle:
                for name, seq in recs:
                    handle.write(f">{name}\n{seq}\n")

        for rnd in range(args.rounds):
            rng = random.Random(args.seed * ROUND_SEED_MULT + rnd)
            kmer_len = rng.choice(KMER_CHOICES)
            recs = genome(rng)
            write(in_dir("g.fa"), recs)
            tag = f"seed {args.seed} round {rnd} k={kmer_len}"
            cnt = {}
            for _, rec_seq in recs:
                for window in kmers(rec_seq, kmer_len):
                    if window:
                        canon = min(window, revcomp(window))
                        cnt[canon] = cnt.get(canon, 0) + 1
            queries = []
            for _, rec_seq in recs:
                if len(rec_seq) >= kmer_len:
                    for _ in range(QUERIES_PER_RECORD):
                        start = rng.randrange(len(rec_seq) - kmer_len + 1)
                        query_seq = rec_seq[start : start + rng.randint(kmer_len, kmer_len + QUERY_EXTRA)]
                        if rng.random() < QUERY_MUTATE_PROB:
                            query_seq = revcomp(query_seq)
                        if rng.random() < QUERY_MUTATE_PROB:
                            query_seq = list(query_seq)
                            query_seq[rng.randrange(len(query_seq))] = rng.choice("ACGTN")
                            query_seq = "".join(query_seq)
                        queries.append(query_seq)
            queries += [
                "".join(rng.choice("ACGT") for _ in range(rng.randint(kmer_len, kmer_len + RANDOM_QUERY_EXTRA)))
                for _ in range(RANDOM_QUERIES)
            ]
            write(in_dir("q.fa"), [(f"q{i}", query_seq) for i, query_seq in enumerate(queries)])

            def want(solid, queries=queries, kmer_len=kmer_len):
                """Expected (positive, negative, invalid) query counts against the solid set."""
                pos = neg = inv = 0
                for query_seq in queries:
                    for kmer in kmers(query_seq, kmer_len):
                        if kmer is None:
                            inv += 1
                        elif min(kmer, revcomp(kmer)) in solid:
                            pos += 1
                        else:
                            neg += 1
                return pos, neg, inv

            # Membership indexes at -c 1 and -c 2
            for min_count in (1, 2):
                solid = {kmer for kmer, count in cnt.items() if count >= min_count}
                builds = (
                    [["build-fm"], ["build-spss"]]
                    if min_count == 1
                    else [
                        ["build-fm", "-c", "2"],
                        ["build-spss", "--compact", "-c", "2"],
                    ]
                )
                for build_args in builds:
                    proc = run(
                        [
                            optik,
                            *build_args,
                            "-i",
                            in_dir("g.fa"),
                            "-k",
                            str(kmer_len),
                            "-t",
                            "4",
                            "-o",
                            in_dir("i.idx"),
                        ]
                    )
                    checks += 1
                    # A genome with no solid k-mer must be rejected, so a successful build is itself a failure.
                    if not solid:
                        if proc.returncode == 0:
                            fails.append(
                                f"{tag} {' '.join(build_args)}: no solid k-mers, yet rc 0"
                            )
                        continue
                    if proc.returncode != 0:
                        fails.append(
                            f"{tag} {' '.join(build_args)}: rc {proc.returncode}: {proc.stderr.strip()[-DETAIL_CHARS:]}"
                        )
                        continue
                    # Only the SPSS index can report per-k-mer ids, so --ids is checked for that build alone.
                    ids = build_args[0] == "build-spss"
                    proc = run(
                        [
                            optik,
                            "query",
                            "-i",
                            in_dir("i.idx"),
                            "-q",
                            in_dir("q.fa"),
                            "-k",
                            str(kmer_len),
                            "-t",
                            "3",
                            *(["--ids", in_dir("ids.txt")] if ids else []),
                        ]
                    )
                    got = tuple(
                        report_value(proc.stdout, f"num_{kind}_kmers")
                        for kind in ("positive", "negative", "invalid")
                    )
                    checks += 1
                    if got != want(solid):
                        fails.append(
                            f"{tag} {' '.join(build_args)} query: {got} != truth {want(solid)}"
                        )
                    # The ids must be a bijection on canonical k-mers, with each k-mer getting one id and no id shared by two k-mers.
                    if ids and proc.returncode == 0:
                        ids_by_kmer, inv = {}, {}
                        for query_seq, line in zip(
                            queries,
                            Path(in_dir("ids.txt")).read_text().split("\n"),
                            strict=False,
                        ):
                            for kmer, kmer_id in zip(kmers(query_seq, kmer_len), line.split(), strict=False):
                                if kmer is not None and kmer_id != "-1":
                                    ids_by_kmer.setdefault(min(kmer, revcomp(kmer)), set()).add(kmer_id)
                        for kmer, kmer_ids in ids_by_kmer.items():
                            for kmer_id in kmer_ids:
                                inv.setdefault(kmer_id, set()).add(kmer)
                        checks += 1
                        multi = sum(len(ids) > 1 for ids in ids_by_kmer.values())
                        shared = sum(len(kmers_of_id) > 1 for kmers_of_id in inv.values())
                        if multi or shared:
                            fails.append(
                                f"{tag} {' '.join(build_args)} --ids: {multi} k-mers with several ids, "
                                f"{shared} ids shared by several k-mers"
                            )

            # spss-update, (base | add) - remove
            def canon_set(seqs, kmer_len=kmer_len):
                """Canonical k-mers of every sequence in seqs."""
                out = set()
                for seq in seqs:
                    for kmer in kmers(seq, kmer_len):
                        if kmer:
                            out.add(min(kmer, revcomp(kmer)))
                return out

            # The update test builds from half the records, adds the other half and removes some query k-mers, then compares against the same edit done on Python sets.
            half = max(1, len(recs) // 2)
            base, addr = recs[:half], recs[half:]
            remr = [(f"x{i}", query_seq) for i, query_seq in enumerate(queries[::REMOVE_STRIDE])]
            write(in_dir("base.fa"), base)
            write(in_dir("add.fa"), addr or [("e", "")])
            write(in_dir("rem.fa"), remr)
            upd = (
                canon_set(rec_seq for _, rec_seq in base) | canon_set(rec_seq for _, rec_seq in addr)
            ) - canon_set(query_seq for _, query_seq in remr)
            if canon_set(rec_seq for _, rec_seq in base):
                proc = run(
                    [
                        optik,
                        "build-spss",
                        "-i",
                        in_dir("base.fa"),
                        "-k",
                        str(kmer_len),
                        "-o",
                        in_dir("u.spss"),
                    ]
                )
                proc = (
                    run(
                        [
                            optik,
                            "spss-update",
                            "-i",
                            in_dir("u.spss"),
                            "--add",
                            in_dir("add.fa"),
                            "--remove",
                            in_dir("rem.fa"),
                            "-o",
                            in_dir("u2.spss"),
                        ]
                    )
                    if proc.returncode == 0
                    else proc
                )
                checks += 1
                # An update that leaves the set empty may fail, while any other failure is reported.
                if proc.returncode != 0:
                    if upd:
                        fails.append(
                            f"{tag} spss-update: rc {proc.returncode}: {proc.stderr.strip()[-DETAIL_CHARS:]}"
                        )
                else:
                    proc = run(
                        [
                            optik,
                            "query",
                            "-i",
                            in_dir("u2.spss"),
                            "-q",
                            in_dir("q.fa"),
                            "-k",
                            str(kmer_len),
                            "-t",
                            "2",
                            "--ids",
                            in_dir("uids.txt"),
                        ]
                    )
                    got = tuple(
                        report_value(proc.stdout, f"num_{kind}_kmers")
                        for kind in ("positive", "negative", "invalid")
                    )
                    if got != want(upd):
                        fails.append(
                            f"{tag} spss-update query: {got} != truth {want(upd)}"
                        )
                    ids_by_kmer = {}
                    for query_seq, line in zip(
                        queries,
                        Path(in_dir("uids.txt")).read_text().split("\n"),
                        strict=False,
                    ):
                        for kmer, kmer_id in zip(kmers(query_seq, kmer_len), line.split(), strict=False):
                            if kmer is not None and kmer_id != "-1":
                                ids_by_kmer.setdefault(min(kmer, revcomp(kmer)), set()).add(kmer_id)
                    checks += 1
                    if any(len(ids) > 1 for ids in ids_by_kmer.values()) or len(
                        {next(iter(ids)) for ids in ids_by_kmer.values()}
                    ) != len(ids_by_kmer):
                        fails.append(
                            f"{tag} spss-update --ids: not one distinct id per k-mer"
                        )

            # FASTQ.gz input gives the same index answers.
            import gzip

            with gzip.open(in_dir("g.fq.gz"), "wt") as handle:
                for name, rec_seq in recs:
                    handle.write(f"@{name}\n{rec_seq}\n+\n{'I' * len(rec_seq)}\n")
            if cnt:
                for build_args in (["build-fm"], ["build-spss"]):
                    build_proc = run(
                        [
                            optik,
                            *build_args,
                            "-i",
                            in_dir("g.fq.gz"),
                            "-k",
                            str(kmer_len),
                            "-o",
                            in_dir("z.idx"),
                        ]
                    )
                    query_proc = (
                        run(
                            [
                                optik,
                                "query",
                                "-i",
                                in_dir("z.idx"),
                                "-q",
                                in_dir("q.fa"),
                                "-k",
                                str(kmer_len),
                            ]
                        )
                        if build_proc.returncode == 0
                        else build_proc
                    )
                    got = tuple(
                        report_value(query_proc.stdout, f"num_{kind}_kmers")
                        for kind in ("positive", "negative", "invalid")
                    )
                    checks += 1
                    if got != want(set(cnt)):
                        fails.append(
                            f"{tag} {build_args[0]} from fastq.gz: {got} != truth {want(set(cnt))}"
                        )

            # dbruijn segments partition the solid set, with exact KC
            if kmer_len <= DBG_MAX_K:
                for min_count in (1, 2):
                    solid = {kmer for kmer, count in cnt.items() if count >= min_count}
                    gfas = []
                    for threads in ("1", "4"):
                        proc = run(
                            [
                                optik,
                                "dbruijn",
                                "-i",
                                in_dir("g.fa"),
                                "-k",
                                str(kmer_len),
                                "-c",
                                str(min_count),
                                "-t",
                                threads,
                                "-o",
                                in_dir(f"d{threads}.gfa"),
                            ]
                        )
                        gfas.append(
                            (
                                proc.returncode,
                                Path(in_dir(f"d{threads}.gfa")).read_text()
                                if proc.returncode == 0
                                else "",
                            )
                        )
                    checks += 2
                    if gfas[0] != gfas[1]:
                        fails.append(
                            f"{tag} dbruijn -c {min_count}: GFA differs between -t 1 and -t 4"
                        )
                    rc0, text = gfas[0]
                    # No k-mer at all is an error, while k-mers of which none is solid give an empty graph, exiting 0 with a warning.
                    if not cnt:
                        if rc0 == 0:
                            fails.append(
                                f"{tag} dbruijn -c {min_count}: input has no k-mers, yet rc 0"
                            )
                        continue
                    if rc0 != 0:
                        fails.append(f"{tag} dbruijn -c {min_count}: rc {rc0}")
                        continue
                    seen, dup, kc_bad = set(), 0, 0
                    for line in text.splitlines():
                        if not line.startswith("S\t"):
                            continue
                        fields = line.split("\t")
                        seq = fields[GfaSegmentCol.SEQ]
                        kc_tag = int(next(field for field in fields if field.startswith("KC:i:"))[SAM_TAG_PREFIX_LEN:])
                        tot = 0
                        for window in kmers(seq, kmer_len):
                            kmer = min(window, revcomp(window))
                            if kmer in seen:
                                dup += 1
                            seen.add(kmer)
                            tot += cnt.get(kmer, 0)
                        kc_bad += tot != kc_tag
                    if seen != solid or dup or kc_bad:
                        fails.append(
                            f"{tag} dbruijn -c {min_count}: {len(seen ^ solid)} k-mers off the solid set, "
                            f"{dup} duplicated, {kc_bad} wrong KC"
                        )

            # map gives one primary per read and the same output at every thread count.
            reads = []
            for _, rec_seq in recs:
                if len(rec_seq) >= SHORT_READ_MIN_REC:
                    for _ in range(SHORT_READS_PER_REC):
                        length = rng.randint(SHORT_READ_MIN_LEN, min(len(rec_seq), SHORT_READ_MAX_LEN))
                        start = rng.randrange(len(rec_seq) - length + 1)
                        read_seq = list(rec_seq[start : start + length])
                        for _ in range(rng.randint(0, MAX_READ_SUBS)):
                            read_seq[rng.randrange(len(read_seq))] = rng.choice("ACGTN")
                        read_seq = "".join(read_seq)
                        reads.append(revcomp(read_seq) if rng.random() < REVCOMP_PROB else read_seq)
                if len(rec_seq) >= LONG_READ_MIN_REC:  # the long-read (chaining) path
                    for _ in range(LONG_READS_PER_REC):
                        length = rng.randint(LONG_READ_MIN_LEN, min(len(rec_seq), LONG_READ_MAX_LEN))
                        start = rng.randrange(len(rec_seq) - length + 1)
                        read_seq = list(rec_seq[start : start + length])
                        # Edits land at random positions about once per LONG_READ_EDIT_SPACING bases, mixing all three kinds.
                        for _ in range(length // LONG_READ_EDIT_SPACING):
                            pos = rng.randrange(len(read_seq))
                            edit = Edit(rng.randrange(len(Edit)))
                            if edit == Edit.SUBSTITUTION:
                                read_seq[pos] = rng.choice("ACGT")
                            elif edit == Edit.INSERTION:
                                read_seq.insert(pos, rng.choice("ACGT"))
                            elif len(read_seq) > 1:
                                del read_seq[pos]
                        read_seq = "".join(read_seq)
                        reads.append(revcomp(read_seq) if rng.random() < REVCOMP_PROB else read_seq)
            if reads:
                write(in_dir("r.fa"), [(f"m{i}", read_seq) for i, read_seq in enumerate(reads)])
                sams = []
                # The same reads at -t 1 and -t 4 must give identical records.
                for threads in ("1", "4"):
                    proc = run(
                        [
                            optik,
                            "map",
                            "-r",
                            in_dir("g.fa"),
                            "-q",
                            in_dir("r.fa"),
                            "-t",
                            threads,
                            "-o",
                            in_dir(f"m{threads}.sam"),
                        ]
                    )
                    # @PG CL records the command line (here -t), so it differs by design.
                    sam_text = ""
                    if proc.returncode == 0:
                        with open(in_dir(f"m{threads}.sam")) as sam_handle:
                            sam_text = "".join(line for line in sam_handle if not line.startswith("@PG"))
                    sams.append((proc.returncode, sam_text))
                checks += 2
                if sams[0] != sams[1]:
                    fails.append(f"{tag} map: SAM differs between -t 1 and -t 4")
                prim = {}
                for line in sams[0][1].splitlines():
                    if line.startswith("@"):
                        continue
                    fields = line.split("\t")
                    if not int(fields[1]) & SamFlag.SECONDARY_SUPPLEMENTARY:
                        prim[fields[0]] = prim.get(fields[0], 0) + 1
                if (
                    sams[0][0] != 0
                    or len(prim) != len(reads)
                    or any(count != 1 for count in prim.values())
                ):
                    fails.append(
                        f"{tag} map: {len(prim)} reads with a primary of {len(reads)} "
                        f"(rc {sams[0][0]})"
                    )

    for failure in fails:
        print("  [FAIL]", failure)
    print(f"fuzz_small: {args.rounds} rounds, {checks} checks, {len(fails)} failures")
    print("fuzz_small: ALL PASS" if not fails else "fuzz_small: FAILED")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
