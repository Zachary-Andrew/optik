#!/usr/bin/env python3
"""edge_cases.py - adversarial differential tests of optik's active commands.

Each check compares optik's output with an independent pure-Python oracle on inputs built to hit parser and kernel edge cases.
These include lowercase and soft-masked bases, N runs, IUPAC codes, CRLF line endings, wrapped FASTA, FASTQ, records shorter than k, empty records, even k with its reverse-complement palindromes, k = 32, small k with dense graphs and cycles, reads hanging off sequence ends, reads with N and multiple references.

  query   num_kmers, positive, negative and invalid counts for build, build-spss and build-fm indexes at several k, single- and multi-threaded
  dbg     dbruijn's GFA against the validate_unitigs oracle (k-mer set, duplicates, paths, compaction, maximality, exact links) at several k and -c
  map     primaries at the simulated strand and origin, POS minus the leading clip, SAM validity of SEQ and CIGAR length and a recomputed NM, well-formed unmapped records and one primary per read
  memory  long reads whose reference-derived core is flanked by tens of kb of foreign sequence, or interrupted by an unrelated 4.9 kb block, must map at 4 threads with a peak RSS of at most 1 GiB, a primary on the core and no record aligning unrelated sequence.
          These guard against earlier regressions, O(s^2) flank and gap WFA memory of GBs per thread and unrelated gaps written as one M run.
  membound peak RSS of the query and graph-build paths on random k-mers, a long single record and a low-complexity genome
  corner  gzip (multi-member, recognized by magic), stdin, multi-line and CRLF FASTQ, '@' quality lines, a missing final newline, blank lines, trailing whitespace, RNA U, long and duplicate read names, reads across two contigs, a circular origin, an N gap or a contig end, the 499/500/501 bp routing boundary, the limits of map -e, -w and -k, abundance above 2^24, --ids with duplicate input k-mers and spss-update corner cases

Usage: python3 bench/edge_cases.py [--optik target/release/optik]
       [--only query,dbg,map,memory,membound,corner]
Any failed check makes the exit status non-zero.
"""

import argparse
import gzip
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import validate_unitigs
from checklog import CheckLog
from diagnostics import DETAIL_CHARS
from dna import revcomp
from sam_spec import (
    DBG_MAX_K,
    SAM_MAX_QNAME_LEN,
    SAM_MIN_COLUMNS,
    SAM_TAG_PREFIX_LEN,
    SamCol,
    SamFlag,
)

# seconds before one optik run counts as hung
RUN_TIMEOUT_S = 3600
# examples listed with a failed check: wrong placements and memory-test records, and missing or duplicated reads
MAX_EXAMPLES = 3
MAX_LISTED_READS = 5
# characters of a CIGAR shown in an example
CIGAR_EXAMPLE_CHARS = 30
LONG_CIGAR_EXAMPLE_CHARS = 40
# a read is placed correctly within max(MIN_ORIGIN_TOL, length / ORIGIN_TOL_DIV) bp of its origin
MIN_ORIGIN_TOL = 10
ORIGIN_TOL_DIV = 10
# memory tests: a primary is on the core within CORE_POS_TOL bp, and no record may have NM above MAX_NM_FRAC of its aligned bases
CORE_POS_TOL = 100
MAX_NM_FRAC = 0.2
DEFAULT_SEED = 11
# placeholder for a k-mer position the --ids output did not reach, the absent-k-mer id
MISSING_ID = -1
check = CheckLog()
FAILS = check.failed


def canon(seq):
    """Canonical form of k-mer string s, upper-cased."""
    seq = seq.upper()
    rev = revcomp(seq)
    return min(seq, rev)


def write_fasta(path, recs, width=60, crlf=False, empty_extra=False):
    """Write records as FASTA wrapped at `width`, optionally with CRLF and a trailing empty record."""
    newline = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as handle:
        for name, seq in recs:
            handle.write(f">{name}{newline}")
            handle.writelines(seq[i : i + width] + newline for i in range(0, len(seq), width))
        if empty_extra:
            handle.write(f">empty{newline}")


def write_fastq(path, recs):
    """Write records as FASTQ with constant quality 'I'."""
    with open(path, "w") as handle:
        handle.writelines(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in recs)


def run(cmd, **kwargs):
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=RUN_TIMEOUT_S, **kwargs, check=False
    )


def kmers_of(seq, kmer_len):
    """Valid canonical k-mers (ACGT only, case-insensitive) and the total number of positions."""
    total = max(0, len(seq) - kmer_len + 1)
    out = []
    upper_seq = seq.upper()
    for i in range(total):
        window = upper_seq[i : i + kmer_len]
        if all(base in "ACGT" for base in window):
            out.append(canon(window))
    return total, out


def make_reference(rng, n_rec=6, length=6000):
    """Random reference records seeded with lowercase stretches, N runs, IUPAC codes and repeats."""
    recs = []
    for r in range(n_rec):
        seq = list("".join(rng.choice("ACGT") for _ in range(length)))
        # A soft-masked (lowercase) stretch, an N run, IUPAC codes and an inverted repeat
        mask_start = rng.randrange(0, length - 400)
        for i in range(mask_start, mask_start + 300):
            seq[i] = seq[i].lower()
        n_start = rng.randrange(0, length - 50)
        for i in range(n_start, n_start + rng.randrange(1, 40)):
            seq[i] = "N"
        for _ in range(5):
            seq[rng.randrange(length)] = rng.choice("RYKMSWn")
        seq = "".join(seq)
        inv_start = rng.randrange(0, length - 500)
        seq = seq[: inv_start + 200] + revcomp(seq[inv_start : inv_start + 200]) + seq[inv_start + 200 :]
        recs.append((f"ref{r} description text", seq))
    recs.append(("short", "ACGTACGTAC"))  # shorter than k
    recs.append(("repeat", recs[0][1][1000:2500]))  # duplicate content
    return recs


def make_queries(rng, refs, kmer_len):
    """Query reads cut from the references with random strand, case and single-base changes, plus edge records."""
    queries = []
    seqs = [seq for _, seq in refs if len(seq) > 400]
    for i in range(400):
        src = rng.choice(seqs)
        read_len = rng.randrange(1, 300)
        start = rng.randrange(0, len(src) - read_len)
        query = src[start : start + read_len]
        if rng.random() < 0.5:
            query = revcomp(query)
        if rng.random() < 0.3:
            query = query.lower()
        if rng.random() < 0.3 and len(query) > 2:
            pos = rng.randrange(len(query))
            query = query[:pos] + rng.choice("ACGTN") + query[pos + 1 :]
        queries.append((f"q{i}", query))
    queries += [
        ("rand", "".join(rng.choice("ACGT") for _ in range(500))),
        ("allN", "N" * 80),
        ("tiny", "ACG"),
        ("exactk", seqs[0][100 : 100 + kmer_len]),
        ("iupac", "ACGTRYACGT" * 10),
    ]
    return queries


# Query suite
def suite_query(optik, tmp_dir, rng):
    """Compare query counts of every index type with the oracle at several k and thread counts."""
    print("== query")
    refs = make_reference(rng)
    ref_fa = os.path.join(tmp_dir, "ref.fa")
    write_fasta(ref_fa, refs, crlf=True, empty_extra=True)
    for kmer_len in (31, 32, 21, 20, 15):
        index = set()
        for _, seq in refs:
            index.update(kmers_of(seq, kmer_len)[1])
        queries = make_queries(rng, refs, kmer_len)
        # The same queries go in as FASTA, FASTQ and wrapped FASTA, and the expected counts come from a Python set of the reference k-mers.
        for fmt in ("fa", "fq", "fa_wrapped"):
            qpath = os.path.join(tmp_dir, f"q_{kmer_len}.{fmt}")
            if fmt == "fq":
                write_fastq(qpath, queries)
            else:
                write_fasta(qpath, queries, width=(37 if fmt == "fa_wrapped" else 10**9))
            exp_total = exp_valid = exp_pos = 0
            for _, query in queries:
                num_windows, kmers = kmers_of(query, kmer_len)
                exp_total += num_windows
                exp_valid += len(kmers)
                exp_pos += sum(1 for kmer in kmers if kmer in index)
            for sub in ("build-spss", "build-fm"):
                if sub == "build-fm" and kmer_len > DBG_MAX_K:
                    continue
                idx = os.path.join(tmp_dir, f"idx_{sub}_{kmer_len}")
                build_proc = run([optik, sub, "-i", ref_fa, "-k", str(kmer_len), "-o", idx])
                if build_proc.returncode != 0:
                    check(f"query k={kmer_len} {sub} build", False, build_proc.stderr.strip()[-DETAIL_CHARS:])
                    continue
                for threads in (1, 4):
                    proc = run(
                        [
                            optik,
                            "query",
                            "-i",
                            idx,
                            "-q",
                            qpath,
                            "-k",
                            str(kmer_len),
                            "-t",
                            str(threads),
                        ]
                    )
                    got = {
                        match.group(1): int(match.group(2))
                        for match in re.finditer(
                            r"^(num_\w+) = (\d+)", proc.stdout, re.MULTILINE
                        )
                    }
                    want = {
                        "num_kmers": exp_total,
                        "num_positive_kmers": exp_pos,
                        "num_negative_kmers": exp_valid - exp_pos,
                        "num_invalid_kmers": exp_total - exp_valid,
                    }
                    bad = {
                        key: (got.get(key), value)
                        for key, value in want.items()
                        if got.get(key) != value
                    }
                    check(
                        f"query k={kmer_len} {sub} {fmt} t={threads}",
                        proc.returncode == 0 and not bad,
                        str(bad) if bad else "",
                    )


# de Bruijn graph suite
def suite_dbg(optik, tmp_dir, rng):
    """Check dbruijn's GFA with the validate_unitigs oracle on reads with cycles, palindromes and hairpins."""
    print("== dbg")
    base = "".join(rng.choice("ACGT") for _ in range(3000))
    reads = []
    for i in range(1500):
        read_len = rng.randrange(20, 160)
        start = rng.randrange(0, len(base) - read_len)
        read_seq = base[start : start + read_len]
        if rng.random() < 0.5:
            read_seq = revcomp(read_seq)
        if rng.random() < 0.2:
            pos = rng.randrange(len(read_seq))
            read_seq = read_seq[:pos] + rng.choice("ACGT") + read_seq[pos + 1 :]
        if rng.random() < 0.1:
            read_seq = read_seq.lower()
        if rng.random() < 0.05:
            pos = rng.randrange(len(read_seq))
            read_seq = read_seq[:pos] + "N" + read_seq[pos + 1 :]
        reads.append((f"r{i}", read_seq))
    # A tandem repeat forming a cycle, a palindrome-rich stretch and a hairpin
    unit = "ACGTTGCAAC"
    reads += [
        ("cycle", unit * 12),
        ("pal", "ACGT" * 30),
        ("hairpin", base[:60] + revcomp(base[:60])),
    ]
    fq_path = os.path.join(tmp_dir, "dbg.fq")
    write_fastq(fq_path, reads)
    fa_path = os.path.join(tmp_dir, "dbg.fa")
    write_fasta(fa_path, reads, width=50, crlf=True)
    for kmer_len in (31, 21, 16, 9, 5):
        for min_count in (1, 2, 3):
            for inp in (fq_path, fa_path):
                gfa = os.path.join(tmp_dir, f"g_{kmer_len}_{min_count}.gfa")
                proc = run(
                    [
                        optik,
                        "dbruijn",
                        "-i",
                        inp,
                        "-k",
                        str(kmer_len),
                        "-c",
                        str(min_count),
                        "-o",
                        gfa,
                        "-t",
                        "4",
                    ]
                )
                tag = f"dbg k={kmer_len} c={min_count} {os.path.basename(inp)}"
                if proc.returncode != 0:
                    check(tag, False, proc.stderr.strip()[-DETAIL_CHARS:])
                    continue
                res = validate_unitigs.validate(inp, gfa, kmer_len, min_count)
                bad = {name: check for name, check in res.items() if not check["pass"]}
                check(tag, not bad, str(bad)[:DETAIL_CHARS] if bad else "")


# Mapping suite
def nm_from(seq, cigar, ref, pos1):
    """Edit distance of a SAM alignment recomputed from SEQ, CIGAR and the reference."""
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    read_pos, ref_pos, nm_count = 0, pos1 - 1, 0
    for op_len, cigar_op in ops:
        op_len = int(op_len)
        if cigar_op == "M":
            nm_count += sum(
                1
                for j in range(op_len)
                if seq[read_pos + j].upper() != ref[ref_pos + j].upper()
                or ref[ref_pos + j].upper() not in "ACGT"
            )
            read_pos += op_len
            ref_pos += op_len
        elif cigar_op in "I":
            nm_count += op_len
            read_pos += op_len
        elif cigar_op == "D":
            nm_count += op_len
            ref_pos += op_len
        elif cigar_op == "S":
            read_pos += op_len
    return nm_count


def suite_map(optik, tmp_dir, rng):
    """Map reads at sequence ends, with N and across references, checking origin, SAM validity and parity."""
    print("== map")
    refs = [
        (f"chr{i}", "".join(rng.choice("ACGT") for _ in range(200000)))
        for i in range(3)
    ]
    refs.append(("tiny", "ACGTACGTACGTTTGACCA"))
    ref_fa = os.path.join(tmp_dir, "mref.fa")
    write_fasta(ref_fa, [(name, seq[:5000].lower() + seq[5000:]) for name, seq in refs], crlf=True)
    ref = {name: seq.upper() for name, seq in refs}
    truth, reads = {}, []
    for i in range(3000):
        name, src_seq = rng.choice(refs[:3])
        read_len = rng.choice([30, 60, 100, 150, 151, 300, 800, 2500])
        # Include reads at both ends of a reference.
        start = rng.choice([0, len(src_seq) - read_len, rng.randrange(0, len(src_seq) - read_len)])
        read_seq = list(src_seq[start : start + read_len])
        for j in range(len(read_seq)):
            if rng.random() < 0.01:
                read_seq[j] = rng.choice("ACGT")
        if rng.random() < 0.1:
            read_seq[rng.randrange(read_len)] = "N"
        read_seq = "".join(read_seq)
        strand = "-" if rng.random() < 0.5 else "+"
        if strand == "-":
            read_seq = revcomp(read_seq)
        if rng.random() < 0.1:
            read_seq = read_seq.lower()
        truth[f"m{i}"] = (name, start, strand, read_len)
        reads.append((f"m{i}", read_seq))
    # Degenerate reads (too short, all N, empty, random) must still give one unmapped record each.
    reads += [
        ("short", "ACGTAC"),
        ("allN", "N" * 200),
        ("empty", ""),
        ("rand", "".join(rng.choice("ACGT") for _ in range(150))),
    ]
    for fmt in ("fq", "fa"):
        qpath = os.path.join(tmp_dir, f"mreads.{fmt}")
        (write_fastq if fmt == "fq" else write_fasta)(qpath, reads)
        sam = os.path.join(tmp_dir, f"m_{fmt}.sam")
        map_proc = run([optik, "map", "-r", ref_fa, "-q", qpath, "-o", sam, "-t", "4"])
        if map_proc.returncode != 0:
            check(f"map {fmt} exit", False, map_proc.stderr.strip()[-DETAIL_CHARS:])
            continue
        prim = Counter()
        unmapped_names = []
        bad_len = bad_nm = bad_seq = wrong = unmapped = bad_unmapped = 0
        examples = []
        readseq = dict(reads)
        with open(sam) as handle:
            for line in handle:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < SAM_MIN_COLUMNS:
                    bad_unmapped += 1
                    continue
                flag = int(fields[SamCol.FLAG])
                if not flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                    prim[fields[SamCol.QNAME]] += 1
                if flag & SamFlag.UNMAPPED:
                    if fields[SamCol.RNAME] != "*" or fields[SamCol.CIGAR] != "*":
                        bad_unmapped += 1
                    if fields[SamCol.QNAME] in truth and not flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                        unmapped += 1
                        unmapped_names.append(fields[SamCol.QNAME])
                    continue
                seq = fields[SamCol.SEQ]
                want = readseq[fields[SamCol.QNAME]]
                # A mapped record must include SEQ and NM:i, so a writer that drops them cannot pass.
                if (
                    seq == "*"
                    or seq.upper()
                    != (revcomp(want) if flag & SamFlag.REVERSE else want).upper()
                ):
                    bad_seq += 1
                ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
                cig_ok = seq != "*" and sum(
                    int(op_len) for op_len, cigar_op in ops if cigar_op in "MIS=X"
                ) == len(seq)
                if not cig_ok:
                    bad_len += 1
                nm_tag = next((int(tag[SAM_TAG_PREFIX_LEN:]) for tag in fields[SamCol.TAGS:] if tag.startswith("NM:i:")), None)
                # NM can be recomputed only when the CIGAR consumes SEQ.
                # NM must equal the edits the CIGAR implies against the reference, recomputed here from the bases.
                if (
                    nm_tag is None
                    or not cig_ok
                    or nm_tag != nm_from(seq, fields[SamCol.CIGAR], ref[fields[SamCol.RNAME]], int(fields[SamCol.POS]))
                ):
                    bad_nm += 1
                if fields[SamCol.QNAME] in truth and not flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                    name, true_pos, strand, read_len = truth[fields[SamCol.QNAME]]
                    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                    correct = (
                        fields[SamCol.RNAME] == name
                        and ((flag & SamFlag.REVERSE) > 0) == (strand == "-")
                        and abs(int(fields[SamCol.POS]) - 1 - lead - true_pos) <= max(MIN_ORIGIN_TOL, read_len // ORIGIN_TOL_DIV)
                    )
                    if not correct:
                        win = ref[fields[SamCol.RNAME]][int(fields[SamCol.POS]) - 1 : int(fields[SamCol.POS]) - 1 + read_len]
                        correct = win == ref[name][true_pos : true_pos + read_len] and nm_tag == 0
                    if not correct:
                        wrong += 1
                        if len(examples) < MAX_EXAMPLES:
                            examples.append(
                                (fields[SamCol.QNAME], truth[fields[SamCol.QNAME]], fields[SamCol.FLAG], fields[SamCol.RNAME], fields[SamCol.POS], fields[SamCol.CIGAR][:CIGAR_EXAMPLE_CHARS])
                            )
        missing = [read_name for read_name, _ in reads if prim[read_name] == 0]
        dup = [read_name for read_name, count in prim.items() if count > 1]
        check(
            f"map {fmt}: one primary record per read",
            not missing and not dup,
            f"missing={missing[:MAX_LISTED_READS]} dup={dup[:MAX_LISTED_READS]}",
        )
        check(f"map {fmt}: SEQ orientation", bad_seq == 0, f"bad={bad_seq}")
        check(f"map {fmt}: CIGAR consumes SEQ", bad_len == 0, f"bad={bad_len}")
        check(f"map {fmt}: NM == CIGAR vs reference", bad_nm == 0, f"bad={bad_nm}")
        check(
            f"map {fmt}: unmapped records well-formed",
            bad_unmapped == 0,
            f"bad={bad_unmapped}",
        )
        n_long = sum(1 for truth_rec in truth.values() if truth_rec[3] >= 60)
        check(
            f"map {fmt}: placements correct",
            wrong == 0,
            f"wrong={wrong} unmapped={unmapped}/{len(truth)} (>=60bp: {n_long}) e.g. {examples}",
        )
        by_len = Counter(truth[read_name][3] for read_name in unmapped_names)
        mm2 = os.path.join(HERE, "..", "temp/minimap2/minimap2")
        if fmt == "fq" and os.path.isfile(mm2):
            # Sensitivity parity with the same reads and truth rule as minimap2 -x sr/map-ont
            o_ok = len(truth) - wrong - unmapped
            m_ok = 0
            for preset, min_len, max_len in (("sr", 0, 1000), ("map-ont", 1000, 10**9)):
                sub = os.path.join(tmp_dir, f"mm_{preset}.fq")
                write_fastq(
                    sub, [(name, seq) for name, seq in reads if name in truth and min_len <= len(seq) < max_len]
                )
                mm_proc = run([mm2, "-a", "-x", preset, ref_fa, sub])
                seen = set()
                for line in mm_proc.stdout.splitlines():
                    if line.startswith("@"):
                        continue
                    fields = line.split("\t")
                    flag = int(fields[SamCol.FLAG])
                    if flag & SamFlag.NOT_PRIMARY or fields[SamCol.QNAME] in seen:
                        continue
                    seen.add(fields[SamCol.QNAME])
                    name, true_pos, strand, read_len = truth[fields[SamCol.QNAME]]
                    ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
                    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                    m_ok += (
                        fields[SamCol.RNAME] == name
                        and ((flag & SamFlag.REVERSE) > 0) == (strand == "-")
                        and abs(int(fields[SamCol.POS]) - 1 - lead - true_pos) <= max(MIN_ORIGIN_TOL, read_len // ORIGIN_TOL_DIV)
                    )
            check(
                f"map {fmt}: correct placements >= minimap2",
                o_ok >= m_ok,
                f"optik {o_ok} vs minimap2 {m_ok} of {len(truth)}; "
                f"optik unmapped by length={dict(by_len)}",
            )


def suite_memory(optik, tmp_dir, rng):
    """Reads with off-reference flanks.
    Before flank extension was chunked, one 15-40 kb flank cost GBs of wavefront memory per thread, which exhausted a 32 GB host."""
    print("== memory")
    ref = "".join(rng.choice("ACGT") for _ in range(300000))
    ref_fa = os.path.join(tmp_dir, "memref.fa")
    write_fasta(ref_fa, [("ref", ref)])
    truth, reads = {}, []
    for i in range(48):
        read_len = rng.randrange(3000, 5000)
        true_pos = rng.randrange(0, len(ref) - read_len)
        core = list(ref[true_pos : true_pos + read_len])
        for j in range(read_len):
            if rng.random() < 0.05:
                core[j] = rng.choice("ACGT")
        left = "".join(rng.choice("ACGT") for _ in range(rng.randrange(0, 20000)))
        right = "".join(rng.choice("ACGT") for _ in range(rng.randrange(10000, 35000)))
        read_seq = left + "".join(core) + right
        strand = "-" if rng.random() < 0.5 else "+"
        reads.append((f"m{i}", revcomp(read_seq) if strand == "-" else read_seq))
        truth[f"m{i}"] = (true_pos, strand)
    # An unrelated internal 4.9 kb block between two matching parts of 5-6 kb must yield a primary and a supplementary record, never one record aligning the block.
    for i in range(16):
        true_pos = rng.randrange(0, len(ref) - 20000)
        read_seq = (
            ref[true_pos : true_pos + 6000]
            + "".join(rng.choice("ACGT") for _ in range(4900))
            + ref[true_pos + 10900 : true_pos + 16000]
        )
        reads.append((f"u{i}", read_seq))
        truth[f"u{i}"] = (true_pos, "+")
    rd_fa = os.path.join(tmp_dir, "memreads.fa")
    write_fasta(rd_fa, reads)
    # The gate is the child's measured peak RSS.
    # The address-space limit only protects the host, because under it the aligner's allocation failure is caught and the flank soft-clipped, so the limit alone cannot reveal the bug.
    cap = 6 << 30
    sam = os.path.join(tmp_dir, "mem.sam")
    with open(sam, "w") as out, open(os.path.join(tmp_dir, "mem.err"), "w") as err:
        proc = subprocess.Popen(
            # The shell sets the limit in the child before exec, since preexec_fn is unsafe in a threaded parent.
            [
                "bash",
                "-c",
                f'ulimit -v {cap >> 10} && exec "$@"',
                "bash",
                optik,
                "map",
                "-r",
                ref_fa,
                "-q",
                rd_fa,
                "-t",
                "4",
            ],
            stdout=out,
            stderr=err,
        )
        _, status, rusage = os.wait4(proc.pid, 0)
    peak_mb = rusage.ru_maxrss / 1024
    check(
        "memory: map exits 0",
        os.waitstatus_to_exitcode(status) == 0,
        f"exit {os.waitstatus_to_exitcode(status)}",
    )
    check(
        "memory: peak RSS <= 1 GiB at 4 threads", peak_mb <= 1024, f"{peak_mb:.0f} MB"
    )
    sam_out = subprocess.CompletedProcess([], 0, Path(sam).read_text(), "")
    bad = []
    correct = 0
    for line in sam_out.stdout.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        flag = int(fields[SamCol.FLAG])
        if flag & SamFlag.NOT_PRIMARY or fields[SamCol.QNAME] not in truth:
            continue
        true_pos, strand = truth[fields[SamCol.QNAME]]
        good = ((flag & SamFlag.REVERSE) > 0) == (strand == "-") and abs(
            int(fields[SamCol.POS]) - 1 - true_pos
        ) <= CORE_POS_TOL
        correct += good
        if not good:
            bad.append(
                f"{fields[SamCol.QNAME]} flag={flag} pos={fields[SamCol.POS]} want={true_pos + 1}{strand} cigar={fields[SamCol.CIGAR][:LONG_CIGAR_EXAMPLE_CHARS]}"
            )
    check(
        "memory: primaries on the reference-derived core",
        correct == len(truth),
        f"{correct}/{len(truth)} {bad[:MAX_EXAMPLES]}",
    )
    # No garbage, meaning every aligned record has NM of at most 20% of its aligned bases
    worst = []
    for line in sam_out.stdout.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
            continue
        aligned = sum(int(op_len) for op_len, cigar_op in re.findall(r"(\d+)([MI])", fields[SamCol.CIGAR]))
        nm_tag = int(next(tag for tag in fields[SamCol.TAGS:] if tag.startswith("NM:i:"))[SAM_TAG_PREFIX_LEN:])
        if nm_tag > MAX_NM_FRAC * aligned:
            worst.append(f"{fields[SamCol.QNAME]} {fields[SamCol.CIGAR][:CIGAR_EXAMPLE_CHARS]} NM={nm_tag}")
    check(
        "memory: no record aligns unrelated sequence (NM <= 20% of aligned)",
        not worst,
        f"{len(worst)} records: {worst[:MAX_EXAMPLES]}",
    )
    supp = sum(
        1
        for line in sam_out.stdout.splitlines()
        if not line.startswith("@")
        and line.split("\t")[0].startswith("u")
        and int(line.split("\t")[1]) & SamFlag.SUPPLEMENTARY
    )
    check(
        "memory: unrelated internal block splits into primary + supplementary",
        supp == 16,
        f"{supp}/16 supplementaries",
    )


def suite_memory_bounds(optik, tmp_dir, rng):
    """Peak RSS of the query and graph-build paths, guarding the fixes of 2026-10-04.

    A query used to keep the whole query file mapped plus a slice per record (about 88 MB for 1M random k-mers, against about 13 MB now).
    The count, graph and walk pipeline held a rank table with 8-byte keys, the counter's freed but unreleased output and buffers glibc retained (about 100 bytes per k-mer, against about 38 now)."""
    import gzip as _gz

    print("== memory: query and graph-build peaks")

    gnu_time = shutil.which("time") or "/usr/bin/time"

    def peak(cmd):
        """Peak RSS of cmd alone, measured through GNU time.

        A child forked from this Python process starts with Python's resident pages counted in its ru_maxrss (Linux folds the address space before exec into it), which here is about 100 MB of test data.
        GNU time forks the command from its own small process."""
        proc = subprocess.run(
            [gnu_time, "-f", "%M", *cmd],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        last = (proc.stderr.strip().splitlines() or ["0"])[-1]
        return proc.returncode, (int(last) / 1024 if last.isdigit() else float("inf"))

    bound_dir = os.path.join(tmp_dir, "membound")
    os.makedirs(bound_dir, exist_ok=True)
    genome = "".join(rng.choice("ACGT") for _ in range(3_000_000))
    gfa = os.path.join(bound_dir, "g.fa")
    write_fasta(gfa, [("g", genome)], width=10**9)  # a single long line per record
    # Query memory must stay flat with the number of records, and index builds must stay within a fixed number of bytes per k-mer.
    idx = os.path.join(bound_dir, "g.spss")
    run([optik, "build-spss", "-i", gfa, "-k", "31", "-o", idx])
    query_fa = os.path.join(bound_dir, "rand.fa")
    with open(query_fa, "w") as handle:
        handle.writelines(
            f">{i}\n{''.join(rng.choice('ACGT') for _ in range(31))}\n"
            for i in range(1_000_000)
        )
    with open(query_fa, "rb") as in_handle, _gz.open(query_fa + ".gz", "wb", compresslevel=1) as out_handle:
        out_handle.write(in_handle.read())
    for path in (query_fa, query_fa + ".gz"):
        for threads in ("1", "8"):
            rc_, peak_mb = peak([optik, "query", "-i", idx, "-q", path, "-k", "31", "-t", threads])
            check(
                f"memory: query 1M records ({os.path.basename(path)}, -t {threads}) peak <= 40 MB",
                rc_ == 0 and peak_mb <= 40,
                f"{peak_mb:.0f} MB",
            )
    n_kmers = len(genome) - 30
    for name, cmd in (
        (
            "build-fm",
            [
                optik,
                "build-fm",
                "-i",
                gfa,
                "-k",
                "31",
                "-t",
                "8",
                "-o",
                os.path.join(bound_dir, "g.fm"),
            ],
        ),
        (
            "dbruijn",
            [
                optik,
                "dbruijn",
                "-i",
                gfa,
                "-k",
                "31",
                "-c",
                "1",
                "-t",
                "8",
                "-o",
                os.path.join(bound_dir, "g.gfa"),
            ],
        ),
        (
            "build-spss --compact",
            [
                optik,
                "build-spss",
                "-i",
                gfa,
                "-k",
                "31",
                "--compact",
                "-t",
                "8",
                "-o",
                os.path.join(bound_dir, "gc.spss"),
            ],
        ),
    ):
        rc_, peak_mb = peak(cmd)
        per = (peak_mb - 10) * 1048576 / n_kmers
        check(
            f"memory: {name} peak <= 60 bytes per k-mer (single 3 Mbp record)",
            rc_ == 0 and per <= 60,
            f"{peak_mb:.0f} MB = {per:.0f} B/k-mer",
        )
    # Low complexity, 8 Mbp holding 3 distinct k-mers.
    # The counter's per-partition table was once sized by occurrences, about 270 MB per pass here, instead of growing with the distinct k-mers it holds.
    # Megabases of poly-A and CG repeat hold almost no distinct k-mers, which must keep memory low instead of blowing up a bucket.
    lowc_fa = os.path.join(bound_dir, "lowc.fa")
    write_fasta(lowc_fa, [("a", "A" * 4_000_000), ("c", "CG" * 2_000_000)], width=10**9)
    for name, cmd in (
        (
            "dbruijn",
            [
                optik,
                "dbruijn",
                "-i",
                lowc_fa,
                "-k",
                "31",
                "-t",
                "8",
                "-o",
                os.path.join(bound_dir, "lc.gfa"),
            ],
        ),
        (
            "build-fm",
            [
                optik,
                "build-fm",
                "-i",
                lowc_fa,
                "-k",
                "31",
                "-t",
                "8",
                "-o",
                os.path.join(bound_dir, "lc.fm"),
            ],
        ),
    ):
        rc_, peak_mb = peak(cmd)
        check(
            f"memory: {name} on 8 Mbp of poly-A / CG repeat (3 distinct k-mers) peak <= 100 MB",
            rc_ == 0 and peak_mb <= 100,
            f"{peak_mb:.0f} MB",
        )


@dataclass
class CornerFixtures:
    """Carries the binary, genomes and reads that the corner sections share.

    It also holds the suite's path and query helpers, so each section can build its own inputs.
    """

    optik_bin: str
    genome: str
    genome2: str
    reads: list
    want: int
    in_dir: Callable
    random_dna: Callable
    counts_of: Callable
    all_found: Callable
    sam: Callable


def corner_input_formats(fx):
    """Feeds unusual input formats, such as multi-member gzip, stdin, CRLF or multi-line FASTQ and '@' quality lines."""
    optik_bin, genome, genome2, reads, in_dir, random_dna, counts_of, all_found, sam = fx.optik_bin, fx.genome, fx.genome2, fx.reads, fx.in_dir, fx.random_dna, fx.counts_of, fx.all_found, fx.sam
    print("== corner: input formats")
    half = len(reads) // 2
    with open(in_dir("mm.fa.gz"), "wb") as handle:  # multi-member gzip, as bgzip writes
        # A multi-member gzip, two compressed parts concatenated, must read as one file.
        handle.writelines(
            gzip.compress("".join(f">{name}\n{seq}\n" for name, seq in part).encode())
            for part in (reads[:half], reads[half:])
        )
    all_found("multi-member gzip", in_dir("mm.fa.gz"))
    shutil.copy(in_dir("mm.fa.gz"), in_dir("mm_noext.fa"))
    all_found("gzip detected by magic bytes (no .gz suffix)", in_dir("mm_noext.fa"))
    all_found("query from stdin '-'", "-", input=Path(in_dir("reads.fa")).read_text())
    with open(in_dir("ml.fq"), "w") as handle:
        handle.writelines(
            f"@{name}\n{seq[:70]}\n{seq[70:]}\n+\n{'I' * 70}\n{'I' * 80}\n" for name, seq in reads
        )
    all_found("multi-line FASTQ", in_dir("ml.fq"))
    with open(in_dir("crlf.fq"), "w", newline="") as handle:
        handle.writelines(f"@{name}\r\n{seq}\r\n+\r\n{'I' * 150}\r\n" for name, seq in reads)
    all_found("CRLF FASTQ", in_dir("crlf.fq"))
    with open(in_dir("qat.fq"), "w") as handle:
        handle.writelines(f"@{name}\n{seq}\n+\n@{'I' * 149}\n" for name, seq in reads)
    all_found("FASTQ quality line starting with '@'", in_dir("qat.fq"))
    Path(in_dir("noeol.fa")).write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in reads).rstrip("\n")
    )
    all_found("FASTA without a final newline", in_dir("noeol.fa"))
    Path(in_dir("blank.fa")).write_text(
        "".join(f">{name}\n{seq[:75]}\n\n{seq[75:]}\n" for name, seq in reads)
    )
    all_found("blank line inside a record", in_dir("blank.fa"))
    Path(in_dir("ws.fa")).write_text(
        "".join(f">{name} desc\t\n{seq[:75]} \t\n{seq[75:]}\t \n" for name, seq in reads)
    )
    all_found("trailing spaces/tabs on wrapped sequence lines", in_dir("ws.fa"))
    Path(in_dir("rna.fa")).write_text(
        "".join(f">{name}\n{seq.replace('T', 'U')}\n" for name, seq in reads)
    )
    all_found("RNA (U read as T)", in_dir("rna.fa"))

    proc = run(
        [
            optik_bin,
            "map",
            "-r",
            in_dir("ref.fa"),
            "-q",
            in_dir("ws.fa"),
            "-o",
            in_dir("ws.sam"),
        ]
    )
    recs = sam(in_dir("ws.sam"))
    clean = all(re.fullmatch(r"[ACGTN*]+", fields[SamCol.SEQ]) for fields in recs)
    check(
        "corner: map SEQ carries no whitespace from the input lines",
        proc.returncode == 0 and clean and len(recs) == len(reads),
    )
    proc = run(
        [
            optik_bin,
            "map",
            "-r",
            in_dir("ref.fa"),
            "-q",
            in_dir("rna.fa"),
            "-o",
            in_dir("rna.sam"),
        ]
    )
    mapped = sum(1 for fields in sam(in_dir("rna.sam")) if not int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED)
    check(
        "corner: map RNA reads (U as T)",
        proc.returncode == 0 and mapped == len(reads),
        f"{mapped}/{len(reads)} mapped",
    )

    # One record of several hundred kb is cut into window-exact slices when -t > 1, and every count must equal the oracle at every thread count and index type.
    # One record far longer than a work chunk is split across threads, and the counts must still match the oracle at every -t.
    big = (
        genome[: 300000 % len(genome)]
        + "N" * 50
        + revcomp(genome)
        + "NNNN"
        + genome2
        + random_dna(200000)
        + genome[:150000]
        + "N"
        + genome2[::-1]
    )[:700000]
    index = set(kmers_of(genome, 31)[1]) | set(kmers_of(genome2, 31)[1])
    recs = [
        ("big", big),
        ("short", genome[500:540]),
        ("tiny", "ACGT"),
        ("big2", revcomp(big[:200000])),
    ]
    write_fasta(in_dir("bigq.fa"), recs, width=10**9)
    want_big = {
        "num_kmers": 0,
        "num_positive_kmers": 0,
        "num_negative_kmers": 0,
        "num_invalid_kmers": 0,
    }
    for _, seq in recs:
        tot, kmers = kmers_of(seq, 31)
        pos = sum(1 for kmer in kmers if kmer in index)
        want_big["num_kmers"] += tot
        want_big["num_positive_kmers"] += pos
        want_big["num_negative_kmers"] += len(kmers) - pos
        want_big["num_invalid_kmers"] += tot - len(kmers)
    for sub in ("build-spss", "build-fm"):
        idx = in_dir(f"bigidx_{sub}")
        run([optik_bin, sub, "-i", in_dir("ref.fa"), "-k", "31", "-o", idx])
        for threads in (1, 2, 8, 16):
            counts = counts_of(
                run(
                    [
                        optik_bin,
                        "query",
                        "-i",
                        idx,
                        "-q",
                        in_dir("bigq.fa"),
                        "-k",
                        "31",
                        "-t",
                        str(threads),
                    ]
                )
            )
            bad = {key: (counts.get(key), value) for key, value in want_big.items() if counts.get(key) != value}
            check(
                f"corner: long record split across threads, {sub} -t {threads}: counts == oracle",
                not bad,
                str(bad),
            )


def corner_read_names(fx):
    """Checks that long, duplicate and unusual read names survive into SAM."""
    optik_bin, reads, in_dir, sam = fx.optik_bin, fx.reads, fx.in_dir, fx.sam
    print("== corner: read names")
    # A 300-character name exceeds SAM's QNAME limit, so map truncates it with a warning and samtools must still read the file.
    long_name = "x" * 300
    write_fasta(
        in_dir("names.fa"),
        [
            (long_name, reads[0][1]),
            ("dup extra words", reads[1][1]),
            ("dup", reads[2][1]),
        ],
        width=10**9,
    )
    proc = run(
        [
            optik_bin,
            "map",
            "-r",
            in_dir("ref.fa"),
            "-q",
            in_dir("names.fa"),
            "-o",
            in_dir("names.sam"),
        ]
    )
    names = [fields[SamCol.QNAME] for fields in sam(in_dir("names.sam"))]
    check(
        "corner: QNAME > 254 truncated to the SAM limit (with a warning)",
        proc.returncode == 0 and names[0] == "x" * SAM_MAX_QNAME_LEN and "truncated" in proc.stderr,
        f"{len(names[0])} chars",
    )
    check(
        "corner: header words after the name dropped; duplicate names kept",
        names[1:] == ["dup", "dup"],
        str(names[1:]),
    )
    if shutil.which("samtools"):
        view_proc = run(["samtools", "view", "-c", in_dir("names.sam")])
        check(
            "corner: samtools reads the SAM with a long read name",
            view_proc.returncode == 0,
            view_proc.stderr.strip()[-DETAIL_CHARS:],
        )


def corner_mapping_geometry(fx):
    """Places reads where geometry gets awkward.

    Some cross contigs, N gaps or a circular origin, others sit at a contig end or right on the long-read routing boundary.
    """
    optik_bin, genome, genome2, in_dir, random_dna, sam = fx.optik_bin, fx.genome, fx.genome2, fx.in_dir, fx.random_dna, fx.sam
    print("== corner: mapping geometry")
    span = genome[-75:] + genome2[:75]
    circ = genome[-3000:] + genome[:3000]
    write_fasta(
        in_dir("geo.fa"),
        [
            ("span", span),
            ("circ", circ),
            ("longer", genome2 + random_dna(500)),
            ("b499", genome[1000:1499]),
            ("b500", genome[2000:2500]),
            ("b501", genome[3000:3501]),
        ],
        width=10**9,
    )
    run(
        [
            optik_bin,
            "map",
            "-r",
            in_dir("ref.fa"),
            "-q",
            in_dir("geo.fa"),
            "-o",
            in_dir("geo.sam"),
        ]
    )
    recs_by_name = {}
    for fields in sam(in_dir("geo.sam")):
        recs_by_name.setdefault(fields[SamCol.QNAME], []).append(fields)

    def recs_of(name):
        """(flag, ref, pos, CIGAR) of every record of read `name`."""
        return [(int(fields[SamCol.FLAG]), fields[SamCol.RNAME], int(fields[SamCol.POS]), fields[SamCol.CIGAR]) for fields in recs_by_name.get(name, [])]

    def split_ok(name, parts):
        """Check for one primary and one supplementary with these (ref, pos, CIGAR) parts.

        Each must have an SA tag, and which half is primary is a tie.
        Clips compare as S, while the primary must soft-clip and the supplementary hard-clip, the convention of minimap2 and BWA-MEM, whose supplementary SEQ holds only the aligned bases."""
        recs = recs_of(name)
        flags = sorted(flag for flag, *_ in recs)
        has_sa = all(any(tag.startswith("SA:Z:") for tag in fields[SamCol.TAGS:]) for fields in recs_by_name.get(name, []))
        clips = all(
            ("H" not in cigar)
            if not flag & SamFlag.SUPPLEMENTARY
            else ("S" not in cigar and "H" in cigar)
            for flag, _, _, cigar in recs
        )
        # Once the clip styles have passed, H and S compare as one, so the parts match whichever half the mapper made primary.
        norm = sorted((rec[1], rec[2], rec[3].replace("H", "S")) for rec in recs)
        return flags == [0, 2048] and norm == sorted(parts) and has_sa and clips

    check(
        "corner: short read across two contigs -> primary + supplementary",
        split_ok("span", [("chr2", 1, "75S75M"), ("chr1", 99926, "75M75S")]),
        str(recs_of("span")),
    )
    check(
        "corner: read across a circular origin -> primary + supplementary",
        split_ok("circ", [("chr1", 97001, "3000M3000S"), ("chr1", 1, "3000S3000M")]),
        str(recs_of("circ")),
    )
    check(
        "corner: read longer than its contig: every contig base aligned to the end",
        recs_of("longer") == [(0, "chr2", 1, "30000M500S")],
        str(recs_of("longer")),
    )
    # These three lengths sit on either side of the switch from the short-read path to chaining, and each must still map end to end.
    for read_name, pos, read_len in (("b499", 1001, 499), ("b500", 2001, 500), ("b501", 3001, 501)):
        check(
            f"corner: {read_len} bp read at the short/long routing boundary",
            recs_of(read_name) == [(0, "chr1", pos, f"{read_len}M")],
            str(recs_of(read_name)),
        )
    # The gapped reference puts 50 kb of N between two parts of the genome, with one read spanning the gap and one made only of N.
    write_fasta(
        in_dir("nref.fa"),
        [("gapped", genome[:40000] + "N" * 50000 + genome[40000:80000])],
        width=10**9,
    )
    write_fasta(
        in_dir("gapread.fa"),
        [("across_gap", genome[39000:41000]), ("in_N", "N" * 200)],
        width=10**9,
    )
    run(
        [
            optik_bin,
            "map",
            "-r",
            in_dir("nref.fa"),
            "-q",
            in_dir("gapread.fa"),
            "-o",
            in_dir("gap.sam"),
        ]
    )
    recs_by_name = {}
    for fields in sam(in_dir("gap.sam")):
        recs_by_name.setdefault(fields[SamCol.QNAME], []).append(fields)
    check(
        "corner: read across a 50 kb N gap: both sides aligned up to the gap",
        split_ok(
            "across_gap",
            [("gapped", 90001, "1000S1000M"), ("gapped", 39001, "1000M1000S")],
        ),
        str(recs_of("across_gap")),
    )
    check(
        "corner: all-N read is one unmapped record",
        [rec[0] & SamFlag.UNMAPPED for rec in recs_of("in_N")] == [SamFlag.UNMAPPED],
        str(recs_of("in_N")),
    )


def corner_map_parameter_limits(fx):
    """Probes the limits of map -e, -w and -k."""
    optik_bin, genome, reads, in_dir, sam = fx.optik_bin, fx.genome, fx.reads, fx.in_dir, fx.sam
    print("== corner: map parameter limits")
    for extra in (
        ["-e", "0"],
        ["-w", "1"],
        ["-k", "32"],
        ["-k", "8"],
        ["-k", "32", "-w", "1"],
    ):
        proc = run(
            [
                optik_bin,
                "map",
                "-r",
                in_dir("ref.fa"),
                "-q",
                in_dir("reads.fa"),
                "-o",
                in_dir("par.sam"),
                *extra,
            ]
        )
        mapped_recs = (
            [fields for fields in sam(in_dir("par.sam")) if not int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED]
            if proc.returncode == 0
            else []
        )
        at_origin = sum(1 for fields in mapped_recs if genome.find(fields[SamCol.SEQ]) + 1 == int(fields[SamCol.POS]))
        check(
            f"corner: map {' '.join(extra)}: every exact read at its origin",
            proc.returncode == 0 and at_origin == len(reads),
            f"{at_origin}/{len(reads)}",
        )


def corner_abundance_and_ids(fx):
    """Tests abundances above the in-slot range, hairpin k-mers and --ids on input with duplicate k-mers."""
    optik_bin, genome, in_dir, random_dna = fx.optik_bin, fx.genome, fx.in_dir, fx.random_dna
    print("== corner: abundance and k-mer ids")
    # Two k-mers above the graph's in-slot count range (2^22 - 1), the poly-A k-mer and the CG-repeat k-mer (whose two strands share one canonical k-mer)
    # Abundances above 2^22 and 2^24 must stay exact instead of saturating at a packed field's limit.
    Path(in_dir("homo.fa")).write_text(
        ">h\n" + "A" * 17_000_000 + "\n>g\n" + "CG" * 2_200_000 + "\n"
    )
    proc = run(
        [
            optik_bin,
            "dbruijn",
            "-i",
            in_dir("homo.fa"),
            "-k",
            "31",
            "-c",
            "1",
            "-t",
            "4",
            "-o",
            in_dir("homo.gfa"),
        ]
    )
    kc_values = (
        re.findall(r"KC:i:(\d+)", Path(in_dir("homo.gfa")).read_text())
        if proc.returncode == 0
        else []
    )
    check(
        "corner: dbruijn abundances above 2^22 and 2^24 exact (no saturation)",
        sorted(map(int, kc_values)) == [4_400_000 - 30, 17_000_000 - 30],
        f"KC {kc_values}",
    )
    write_fasta(
        in_dir("dupk.fa"),
        [("a", genome[:5000]), ("b", genome[:5000]), ("c", revcomp(genome[:5000]))],
        width=10**9,
    )
    proc = run(
        [
            optik_bin,
            "build-spss",
            "-i",
            in_dir("dupk.fa"),
            "-k",
            "31",
            "-o",
            in_dir("dupk.spss"),
        ]
    )
    distinct_match = re.search(r"distinct_kmers\t(\d+)", proc.stdout)
    check(
        "corner: build-spss distinct_kmers counts distinct k-mers, not positions",
        distinct_match is not None and int(distinct_match.group(1)) == 5000 - 30,
        f"{distinct_match.group(1) if distinct_match else None}",
    )
    check(
        "corner: build-spss warns that duplicate k-mers make --ids a non-bijection",
        "not a bijection" in proc.stderr,
    )
    write_fasta(
        in_dir("idq.fa"), [("f", genome[:5000]), ("r", revcomp(genome[:5000]))], width=10**9
    )
    proc = run(
        [
            optik_bin,
            "query",
            "-i",
            in_dir("dupk.spss"),
            "-q",
            in_dir("idq.fa"),
            "-k",
            "31",
            "--ids",
            in_dir("ids.txt"),
        ]
    )
    lines = [[], []]
    if proc.returncode == 0:
        with open(in_dir("ids.txt")) as handle:
            lines = [list(map(int, line.split())) for line in handle]
    check(
        "corner: --ids with duplicate input k-mers: one id per distinct k-mer, strand-independent",
        len(lines) == 2
        and len(set(lines[0])) == 5000 - 30
        and min(lines[0]) >= 0
        and lines[0] == lines[1][::-1],
    )

    # Hairpins.
    # A 36 bp reverse-complement palindrome y + rc(y) holds pairs of k-mers (at offsets i and 5 - i) that are each other's reverse complement and share one minimizer occurrence, i.e. one super-k-mer record.
    # Both strands must get the same single id.
    def palindrome(length):
        """A random n-base sequence followed by its reverse complement."""
        half = random_dna(length)
        return half + revcomp(half)

    # Palindromes are their own reverse complements, which tests that --ids still gives one id per distinct k-mer on either strand.
    hairpin_seq = "".join(random_dna(40) + palindrome(18) for _ in range(200)) + random_dna(40)
    write_fasta(in_dir("hairpin.fa"), [("h", hairpin_seq)], width=10**9)
    proc = run(
        [
            optik_bin,
            "build-spss",
            "-i",
            in_dir("hairpin.fa"),
            "-k",
            "31",
            "-o",
            in_dir("hairpin.spss"),
        ]
    )
    ids_proc = (
        run(
            [
                optik_bin,
                "query",
                "-i",
                in_dir("hairpin.spss"),
                "-q",
                in_dir("hairpin.fa"),
                "-k",
                "31",
                "--ids",
                in_dir("hp.txt"),
            ]
        )
        if proc.returncode == 0
        else proc
    )
    ids = (
        list(map(int, Path(in_dir("hp.txt")).read_text().split()))
        if ids_proc.returncode == 0
        else []
    )
    ids_by_kmer = {}
    for i, kmer in enumerate(
        kmers_of(hairpin_seq, 31)[1]
    ):  # canonical, at every position of an N-free sequence
        ids_by_kmer.setdefault(kmer, set()).add(ids[i] if i < len(ids) else MISSING_ID)
    distinct_match = re.search(r"distinct_kmers\t(\d+)", proc.stdout) if proc.returncode == 0 else None
    check(
        "corner: hairpin k-mers: one id per distinct k-mer on both strands, distinct count exact",
        len(ids) == len(hairpin_seq) - 30
        and all(len(id_set) == 1 for id_set in ids_by_kmer.values())
        and len({next(iter(id_set)) for id_set in ids_by_kmer.values()}) == len(ids_by_kmer)
        and distinct_match is not None
        and int(distinct_match.group(1)) == len(ids_by_kmer),
        f"{sum(len(id_set) > 1 for id_set in ids_by_kmer.values())} multi-id, distinct {distinct_match.group(1) if distinct_match else None} vs {len(ids_by_kmer)}",
    )


def corner_spss_update(fx):
    """Runs spss-update on edits that change nothing, on removing every k-mer and on in-place updates."""
    optik_bin, genome, genome2, want, in_dir, random_dna, counts_of = fx.optik_bin, fx.genome, fx.genome2, fx.want, fx.in_dir, fx.random_dna, fx.counts_of
    print("== corner: spss-update")
    # spss-update gets k-mers already present, absent ones and the same k-mers added and removed, and the counts after each must match.
    write_fasta(in_dir("present.fa"), [("p", genome[:2000])])
    write_fasta(in_dir("absent.fa"), [("x", random_dna(2000))])
    write_fasta(in_dir("everything.fa"), [("e1", genome), ("e2", genome2)])
    before = len(set(kmers_of(genome, 31)[1]) | set(kmers_of(genome2, 31)[1]))
    for name, extra, after, pos in (
        ("add already-present k-mers", ["--add", in_dir("present.fa")], before, want),
        ("remove absent k-mers", ["--remove", in_dir("absent.fa")], before, want),
        (
            "add and remove the same k-mers",
            ["--add", in_dir("absent.fa"), "--remove", in_dir("absent.fa")],
            before,
            want,
        ),
        ("remove every k-mer", ["--remove", in_dir("everything.fa")], 0, 0),
    ):
        proc = run(
            [
                optik_bin,
                "spss-update",
                "-i",
                in_dir("ref.spss"),
                "-o",
                in_dir("u.spss"),
                *extra,
            ]
        )
        after_match = re.search(r"kmers_after\t(\d+)", proc.stdout)
        counts = counts_of(
            run(
                [
                    optik_bin,
                    "query",
                    "-i",
                    in_dir("u.spss"),
                    "-q",
                    in_dir("reads.fa"),
                    "-k",
                    "31",
                ]
            )
        )
        check(
            f"corner: spss-update {name}",
            proc.returncode == 0
            and after_match is not None
            and int(after_match.group(1)) == after
            and counts.get("num_positive_kmers") == pos,
            f"kmers_after {after_match.group(1) if after_match else None} (want {after}), positives {counts.get('num_positive_kmers')} (want {pos})",
        )
    shutil.copy(in_dir("ref.spss"), in_dir("inplace.spss"))
    write_fasta(in_dir("add.fa"), [("a", random_dna(5000))])
    proc = run(
        [
            optik_bin,
            "spss-update",
            "-i",
            in_dir("inplace.spss"),
            "--add",
            in_dir("add.fa"),
            "-o",
            in_dir("inplace.spss"),
        ]
    )
    counts = counts_of(
        run(
            [
                optik_bin,
                "query",
                "-i",
                in_dir("inplace.spss"),
                "-q",
                in_dir("add.fa"),
                "-k",
                "31",
            ]
        )
    )
    check(
        "corner: spss-update in place (-i X -o X) keeps a valid index with the additions",
        proc.returncode == 0 and counts.get("num_positive_kmers") == 5000 - 30,
        str(counts),
    )


def suite_corner(optik, tmp_dir, rng):
    """Input-format, geometry and parameter corner cases.

    Each was found by a probe that showed wrong or silently partial output before its fix."""

    optik_bin = optik
    corner_dir = os.path.join(tmp_dir, "corner")
    os.makedirs(corner_dir, exist_ok=True)

    def in_dir(name):
        """Path of `name` in this suite's directory."""
        return os.path.join(corner_dir, name)

    def random_dna(length):
        """Random ACGT string of length n."""
        return "".join(rng.choice("ACGT") for _ in range(length))

    # The corner sections run on this two-contig genome, on 200 exact reads drawn from it (every other one reverse-complemented) and on an SPSS index of it.
    genome, genome2 = random_dna(100000), random_dna(30000)
    write_fasta(in_dir("ref.fa"), [("chr1", genome), ("chr2", genome2)], width=10**9)
    reads = [
        (f"r{i}", genome[start : start + 150])
        for i, start in enumerate(rng.sample(range(99000), 200))
    ]
    reads = [(name, revcomp(seq) if i % 2 else seq) for i, (name, seq) in enumerate(reads)]
    write_fasta(in_dir("reads.fa"), reads, width=10**9)
    run(
        [
            optik_bin,
            "build-spss",
            "-i",
            in_dir("ref.fa"),
            "-k",
            "31",
            "-o",
            in_dir("ref.spss"),
        ]
    )
    # Every 150 bp read holds 120 k-mers of length 31, all present in the index.
    want = 200 * (150 - 31 + 1)

    def counts_of(proc):
        """The num_* counters printed by an optik query run."""
        return {
            match.group(1): int(match.group(2))
            for match in re.finditer(r"^(num_\w+) = (\d+)", proc.stdout, re.MULTILINE)
        }

    def query_counts(path, **kwargs):
        """Query the SPSS index with `path`, returning (exit code, num_kmers, num_positive_kmers)."""
        proc = run(
            [optik_bin, "query", "-i", in_dir("ref.spss"), "-q", path, "-k", "31"], **kwargs
        )
        counters = counts_of(proc)
        return proc.returncode, counters.get("num_kmers"), counters.get("num_positive_kmers")

    def all_found(name, path, **kwargs):
        """Check that querying `path` finds all `want` k-mers."""
        rcode, num_kmers, pos = query_counts(path, **kwargs)
        check(
            f"corner: {name}",
            rcode == 0 and num_kmers == want and pos == want,
            f"exit {rcode}, {pos}/{num_kmers} positive (want {want}/{want})",
        )

    def sam(path):
        """SAM records of a file as field lists, header excluded."""
        with open(path) as handle:
            return [
                line.rstrip("\n").split("\t") for line in handle if not line.startswith("@")
            ]
    fx = CornerFixtures(optik_bin=optik_bin, genome=genome, genome2=genome2, reads=reads, want=want, in_dir=in_dir, random_dna=random_dna, counts_of=counts_of, all_found=all_found, sam=sam)
    for section in (corner_input_formats, corner_read_names, corner_mapping_geometry, corner_map_parameter_limits, corner_abundance_and_ids, corner_spss_update):
        section(fx)


def main():
    """Run the suites chosen with --only."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    parser.add_argument("--only", default="query,dbg,map,memory,membound,corner")
    parser.add_argument("-s", "--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    with tempfile.TemporaryDirectory(prefix="optik_edge_") as tmp_dir:
        for name in args.only.split(","):
            {
                "query": suite_query,
                "dbg": suite_dbg,
                "map": suite_map,
                "memory": suite_memory,
                "membound": suite_memory_bounds,
                "corner": suite_corner,
            }[name](args.optik, tmp_dir, rng)
    print(f"edge_cases: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
