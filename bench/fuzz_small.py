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
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from sam_spec import DBG_MAX_K, SAM_FLAGS_SECONDARY_SUPPLEMENTARY

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    """Reverse complement of s."""
    return s.translate(COMP)[::-1]


def kmers(s, k):
    """Each k-mer window of s, upper-cased, or None for a window holding a non-ACGT base."""
    s = s.upper()
    for i in range(len(s) - k + 1):
        w = s[i : i + k]
        yield w if all(c in "ACGT" for c in w) else None


def genome(r):
    """Random genome of one to six records built from the adversarial shapes."""
    recs = []
    for i in range(r.randint(1, 6)):
        s = []
        for _ in range(r.randint(1, 12)):
            kind = r.randrange(8)
            if kind == 0:
                s.append("".join(r.choice("ACGT") for _ in range(r.randint(1, 400))))
            elif kind == 1:
                y = "".join(r.choice("ACGT") for _ in range(r.randint(1, 30)))
                s.append(y + rc(y))
            elif kind == 2:
                s.append(r.choice("ACGT") * r.randint(1, 200))
            elif kind == 3:
                u = "".join(r.choice("ACGT") for _ in range(r.randint(1, 7)))
                s.append(u * r.randint(2, 60))
            elif kind == 4:
                s.append("N" * r.randint(1, 20))
            elif kind == 5 and s:
                s.append(s[r.randrange(len(s))])
            elif kind == 6 and s:
                s.append(rc(s[r.randrange(len(s))]))
            else:
                s.append(
                    "".join(r.choice("ACGT") for _ in range(r.randint(1, 60))).lower()
                )
        recs.append((f"r{i}", "".join(s)))
    return recs


def main():
    """Run every fuzz round and return 1 if any check failed, else 0."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    ap.add_argument("--rounds", type=int, default=80)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    optik = a.optik
    fails = []
    checks = 0

    def run(cmd):
        return subprocess.run(
            cmd, capture_output=True, text=True, timeout=600, check=False
        )

    def kv(out, key):
        """Integer value of the report line starting with `key`."""
        for line in out.splitlines():
            if line.startswith(key):
                return int(line.split("=")[1].split()[0])
        return None

    with tempfile.TemporaryDirectory(prefix="optik_fuzz_") as td:

        def in_dir(f):
            """Path of f in this round's directory."""
            return os.path.join(td, f)

        def write(path, recs):
            """Write records as single-line FASTA."""
            with open(path, "w") as f:
                for n, s in recs:
                    f.write(f">{n}\n{s}\n")

        for rnd in range(a.rounds):
            r = random.Random(a.seed * 100003 + rnd)
            k = r.choice([3, 4, 5, 8, 11, 15, 16, 20, 21, 25, 31, 32])
            recs = genome(r)
            write(in_dir("g.fa"), recs)
            tag = f"seed {a.seed} round {rnd} k={k}"
            cnt = {}
            for _, s in recs:
                for w in kmers(s, k):
                    if w:
                        c = min(w, rc(w))
                        cnt[c] = cnt.get(c, 0) + 1
            q = []
            for _, s in recs:
                if len(s) >= k:
                    for _ in range(20):
                        st = r.randrange(len(s) - k + 1)
                        w = s[st : st + r.randint(k, k + 40)]
                        if r.random() < 0.3:
                            w = rc(w)
                        if r.random() < 0.3:
                            w = list(w)
                            w[r.randrange(len(w))] = r.choice("ACGTN")
                            w = "".join(w)
                        q.append(w)
            q += [
                "".join(r.choice("ACGT") for _ in range(r.randint(k, k + 30)))
                for _ in range(50)
            ]
            write(in_dir("q.fa"), [(f"q{i}", w) for i, w in enumerate(q)])

            def want(solid, q=q, k=k):
                """Expected (positive, negative, invalid) query counts against the solid set."""
                pos = neg = inv = 0
                for w in q:
                    for x in kmers(w, k):
                        if x is None:
                            inv += 1
                        elif min(x, rc(x)) in solid:
                            pos += 1
                        else:
                            neg += 1
                return pos, neg, inv

            # Membership indexes at -c 1 and -c 2
            for c in (1, 2):
                solid = {x for x, v in cnt.items() if v >= c}
                builds = (
                    [["build-fm"], ["build-spss"]]
                    if c == 1
                    else [
                        ["build-fm", "-c", "2"],
                        ["build-spss", "--compact", "-c", "2"],
                    ]
                )
                for b in builds:
                    p = run(
                        [
                            optik,
                            *b,
                            "-i",
                            in_dir("g.fa"),
                            "-k",
                            str(k),
                            "-t",
                            "4",
                            "-o",
                            in_dir("i.idx"),
                        ]
                    )
                    checks += 1
                    if not solid:
                        if p.returncode == 0:
                            fails.append(
                                f"{tag} {' '.join(b)}: no solid k-mers, yet rc 0"
                            )
                        continue
                    if p.returncode != 0:
                        fails.append(
                            f"{tag} {' '.join(b)}: rc {p.returncode}: {p.stderr.strip()[-160:]}"
                        )
                        continue
                    ids = b[0] == "build-spss"
                    p = run(
                        [
                            optik,
                            "query",
                            "-i",
                            in_dir("i.idx"),
                            "-q",
                            in_dir("q.fa"),
                            "-k",
                            str(k),
                            "-t",
                            "3",
                            *(["--ids", in_dir("ids.txt")] if ids else []),
                        ]
                    )
                    got = tuple(
                        kv(p.stdout, f"num_{x}_kmers")
                        for x in ("positive", "negative", "invalid")
                    )
                    checks += 1
                    if got != want(solid):
                        fails.append(
                            f"{tag} {' '.join(b)} query: {got} != truth {want(solid)}"
                        )
                    if ids and p.returncode == 0:
                        by, inv = {}, {}
                        for w, line in zip(
                            q,
                            Path(in_dir("ids.txt")).read_text().split("\n"),
                            strict=False,
                        ):
                            for x, v in zip(kmers(w, k), line.split(), strict=False):
                                if x is not None and v != "-1":
                                    by.setdefault(min(x, rc(x)), set()).add(v)
                        for x, vs in by.items():
                            for v in vs:
                                inv.setdefault(v, set()).add(x)
                        checks += 1
                        multi = sum(len(v) > 1 for v in by.values())
                        shared = sum(len(v) > 1 for v in inv.values())
                        if multi or shared:
                            fails.append(
                                f"{tag} {' '.join(b)} --ids: {multi} k-mers with several ids, "
                                f"{shared} ids shared by several k-mers"
                            )

            # spss-update, (base | add) - remove
            def canon_set(seqs, k=k):
                """Canonical k-mers of every sequence in seqs."""
                out = set()
                for x in seqs:
                    for w in kmers(x, k):
                        if w:
                            out.add(min(w, rc(w)))
                return out

            half = max(1, len(recs) // 2)
            base, addr = recs[:half], recs[half:]
            remr = [(f"x{i}", w) for i, w in enumerate(q[::3])]
            write(in_dir("base.fa"), base)
            write(in_dir("add.fa"), addr or [("e", "")])
            write(in_dir("rem.fa"), remr)
            upd = (
                canon_set(s for _, s in base) | canon_set(s for _, s in addr)
            ) - canon_set(w for _, w in remr)
            if canon_set(s for _, s in base):
                p = run(
                    [
                        optik,
                        "build-spss",
                        "-i",
                        in_dir("base.fa"),
                        "-k",
                        str(k),
                        "-o",
                        in_dir("u.spss"),
                    ]
                )
                p = (
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
                    if p.returncode == 0
                    else p
                )
                checks += 1
                if p.returncode != 0:
                    if upd:
                        fails.append(
                            f"{tag} spss-update: rc {p.returncode}: {p.stderr.strip()[-160:]}"
                        )
                else:
                    p = run(
                        [
                            optik,
                            "query",
                            "-i",
                            in_dir("u2.spss"),
                            "-q",
                            in_dir("q.fa"),
                            "-k",
                            str(k),
                            "-t",
                            "2",
                            "--ids",
                            in_dir("uids.txt"),
                        ]
                    )
                    got = tuple(
                        kv(p.stdout, f"num_{x}_kmers")
                        for x in ("positive", "negative", "invalid")
                    )
                    if got != want(upd):
                        fails.append(
                            f"{tag} spss-update query: {got} != truth {want(upd)}"
                        )
                    by = {}
                    for w, line in zip(
                        q,
                        Path(in_dir("uids.txt")).read_text().split("\n"),
                        strict=False,
                    ):
                        for x, v in zip(kmers(w, k), line.split(), strict=False):
                            if x is not None and v != "-1":
                                by.setdefault(min(x, rc(x)), set()).add(v)
                    checks += 1
                    if any(len(v) > 1 for v in by.values()) or len(
                        {next(iter(v)) for v in by.values()}
                    ) != len(by):
                        fails.append(
                            f"{tag} spss-update --ids: not one distinct id per k-mer"
                        )

            # FASTQ.gz input gives the same index answers.
            import gzip

            with gzip.open(in_dir("g.fq.gz"), "wt") as f:
                for n, x in recs:
                    f.write(f"@{n}\n{x}\n+\n{'I' * len(x)}\n")
            if cnt:
                for b in (["build-fm"], ["build-spss"]):
                    p1 = run(
                        [
                            optik,
                            *b,
                            "-i",
                            in_dir("g.fq.gz"),
                            "-k",
                            str(k),
                            "-o",
                            in_dir("z.idx"),
                        ]
                    )
                    p2 = (
                        run(
                            [
                                optik,
                                "query",
                                "-i",
                                in_dir("z.idx"),
                                "-q",
                                in_dir("q.fa"),
                                "-k",
                                str(k),
                            ]
                        )
                        if p1.returncode == 0
                        else p1
                    )
                    got = tuple(
                        kv(p2.stdout, f"num_{x}_kmers")
                        for x in ("positive", "negative", "invalid")
                    )
                    checks += 1
                    if got != want(set(cnt)):
                        fails.append(
                            f"{tag} {b[0]} from fastq.gz: {got} != truth {want(set(cnt))}"
                        )

            # dbruijn segments partition the solid set, with exact KC
            if k <= DBG_MAX_K:
                for c in (1, 2):
                    solid = {x for x, v in cnt.items() if v >= c}
                    gfas = []
                    for t in ("1", "4"):
                        p = run(
                            [
                                optik,
                                "dbruijn",
                                "-i",
                                in_dir("g.fa"),
                                "-k",
                                str(k),
                                "-c",
                                str(c),
                                "-t",
                                t,
                                "-o",
                                in_dir(f"d{t}.gfa"),
                            ]
                        )
                        gfas.append(
                            (
                                p.returncode,
                                Path(in_dir(f"d{t}.gfa")).read_text()
                                if p.returncode == 0
                                else "",
                            )
                        )
                    checks += 2
                    if gfas[0] != gfas[1]:
                        fails.append(
                            f"{tag} dbruijn -c {c}: GFA differs between -t 1 and -t 4"
                        )
                    rc0, text = gfas[0]
                    # No k-mer at all is an error, while k-mers of which none is solid give an empty graph, exiting 0 with a warning.
                    if not cnt:
                        if rc0 == 0:
                            fails.append(
                                f"{tag} dbruijn -c {c}: input has no k-mers, yet rc 0"
                            )
                        continue
                    if rc0 != 0:
                        fails.append(f"{tag} dbruijn -c {c}: rc {rc0}")
                        continue
                    seen, dup, kc_bad = set(), 0, 0
                    for line in text.splitlines():
                        if not line.startswith("S\t"):
                            continue
                        f = line.split("\t")
                        seq = f[2]
                        kc = int(next(x for x in f if x.startswith("KC:i:"))[5:])
                        tot = 0
                        for w in kmers(seq, k):
                            x = min(w, rc(w))
                            if x in seen:
                                dup += 1
                            seen.add(x)
                            tot += cnt.get(x, 0)
                        kc_bad += tot != kc
                    if seen != solid or dup or kc_bad:
                        fails.append(
                            f"{tag} dbruijn -c {c}: {len(seen ^ solid)} k-mers off the solid set, "
                            f"{dup} duplicated, {kc_bad} wrong KC"
                        )

            # map gives one primary per read and the same output at every thread count.
            reads = []
            for _, s in recs:
                if len(s) >= 60:
                    for _ in range(15):
                        length = r.randint(30, min(len(s), 300))
                        st = r.randrange(len(s) - length + 1)
                        w = list(s[st : st + length])
                        for _ in range(r.randint(0, 6)):
                            w[r.randrange(len(w))] = r.choice("ACGTN")
                        w = "".join(w)
                        reads.append(rc(w) if r.random() < 0.5 else w)
                if len(s) >= 600:  # the long-read (chaining) path
                    for _ in range(3):
                        length = r.randint(500, min(len(s), 3000))
                        st = r.randrange(len(s) - length + 1)
                        w = list(s[st : st + length])
                        for _ in range(length // 20):
                            j = r.randrange(len(w))
                            op = r.randrange(3)
                            if op == 0:
                                w[j] = r.choice("ACGT")
                            elif op == 1:
                                w.insert(j, r.choice("ACGT"))
                            elif len(w) > 1:
                                del w[j]
                        w = "".join(w)
                        reads.append(rc(w) if r.random() < 0.5 else w)
            if reads:
                write(in_dir("r.fa"), [(f"m{i}", w) for i, w in enumerate(reads)])
                sams = []
                for t in ("1", "4"):
                    p = run(
                        [
                            optik,
                            "map",
                            "-r",
                            in_dir("g.fa"),
                            "-q",
                            in_dir("r.fa"),
                            "-t",
                            t,
                            "-o",
                            in_dir(f"m{t}.sam"),
                        ]
                    )
                    # @PG CL records the command line (here -t), so it differs by design.
                    sam_text = ""
                    if p.returncode == 0:
                        with open(in_dir(f"m{t}.sam")) as fh:
                            sam_text = "".join(x for x in fh if not x.startswith("@PG"))
                    sams.append((p.returncode, sam_text))
                checks += 2
                if sams[0] != sams[1]:
                    fails.append(f"{tag} map: SAM differs between -t 1 and -t 4")
                prim = {}
                for line in sams[0][1].splitlines():
                    if line.startswith("@"):
                        continue
                    f = line.split("\t")
                    if not int(f[1]) & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                        prim[f[0]] = prim.get(f[0], 0) + 1
                if (
                    sams[0][0] != 0
                    or len(prim) != len(reads)
                    or any(v != 1 for v in prim.values())
                ):
                    fails.append(
                        f"{tag} map: {len(prim)} reads with a primary of {len(reads)} "
                        f"(rc {sams[0][0]})"
                    )

    for f in fails:
        print("  [FAIL]", f)
    print(f"fuzz_small: {a.rounds} rounds, {checks} checks, {len(fails)} failures")
    print("fuzz_small: ALL PASS" if not fails else "fuzz_small: FAILED")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
