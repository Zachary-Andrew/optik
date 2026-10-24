#!/usr/bin/env python3
"""edge_cases.py - adversarial differential tests of optik's active commands.

Each check compares optik's output with an independent pure-Python oracle on inputs built to hit parser and kernel edge cases.
These include lowercase and soft-masked bases, N runs, IUPAC codes, CRLF line endings, wrapped FASTA, FASTQ, records shorter than k, empty records, even k with its reverse-complement palindromes, k = 32, small k with dense graphs and cycles, reads hanging off sequence ends, reads with N and multiple references.

  query   num_kmers, positive, negative and invalid counts for build, build-spss and build-fm indexes at several k, single- and multi-threaded
  dbg     dbruijn's GFA against the validate_unitigs oracle (k-mer set, duplicates, paths, compaction, maximality, exact links) at several k and -c
  map     primaries at the simulated strand and origin, POS minus the leading clip, SAM validity of SEQ and CIGAR length and a recomputed NM, well-formed unmapped records and one primary per read
  memory  long reads whose reference-derived core is flanked by tens of kb of foreign sequence, or interrupted by an unrelated 4.9 kb block, must map at 4 threads with a peak RSS of at most 1 GiB, a primary on the core and no record aligning unrelated sequence.
          These guard against earlier regressions, O(s^2) flank and gap WFA memory of GBs per thread and unrelated gaps written as one M run.
  membound peak RSS of the query and graph-build paths on random k-mers,
          a long single record and a low-complexity genome
  corner  gzip (multi-member, recognized by magic), stdin, multi-line and CRLF FASTQ, '@' quality lines, a missing final newline, blank lines, trailing whitespace, RNA U, long and duplicate read names, reads across two contigs, a circular origin, an N gap or a contig end, the 499/500/501 bp routing boundary, the limits of map -e, -w and -k, abundance above 2^24, --ids with duplicate input k-mers and spss-update corner cases

Usage: python3 bench/edge_cases.py [--optik target/release/optik]
       [--only query,dbg,map,memory,membound,corner]
Any failed check makes the exit status non-zero.
"""

import argparse
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import validate_unitigs as vu
from sam_spec import (
    DBG_MAX_K,
    SAM_FLAG_REVERSE,
    SAM_FLAG_SUPPLEMENTARY,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_NOT_PRIMARY,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
    SAM_MIN_COLUMNS,
)

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
FAILS = []


def rc(s):
    """Reverse complement of s."""
    return s.translate(COMP)[::-1]


def canon(s):
    """Canonical form of k-mer string s, upper-cased."""
    s = s.upper()
    r = rc(s)
    return min(s, r)


def check(name, ok, detail=""):
    """Report a check as PASS or FAIL, adding failures to FAILS."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def write_fasta(path, recs, width=60, crlf=False, empty_extra=False):
    """Write records as FASTA wrapped at `width`, optionally with CRLF and a trailing empty record."""
    nl = "\r\n" if crlf else "\n"
    with open(path, "w", newline="") as fh:
        for name, s in recs:
            fh.write(f">{name}{nl}")
            fh.writelines(s[i : i + width] + nl for i in range(0, len(s), width))
        if empty_extra:
            fh.write(f">empty{nl}")


def write_fastq(path, recs):
    """Write records as FASTQ with constant quality 'I'."""
    with open(path, "w") as fh:
        fh.writelines(f"@{name}\n{s}\n+\n{'I' * len(s)}\n" for name, s in recs)


def run(cmd, **kw):
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=3600, **kw, check=False
    )


def kmers_of(s, k):
    """Valid canonical k-mers (ACGT only, case-insensitive) and the total number of positions."""
    total = max(0, len(s) - k + 1)
    out = []
    u = s.upper()
    for i in range(total):
        w = u[i : i + k]
        if all(c in "ACGT" for c in w):
            out.append(canon(w))
    return total, out


def make_reference(rng, n_rec=6, length=6000):
    """Random reference records seeded with lowercase stretches, N runs, IUPAC codes and repeats."""
    recs = []
    for r in range(n_rec):
        s = list("".join(rng.choice("ACGT") for _ in range(length)))
        # A soft-masked (lowercase) stretch, an N run, IUPAC codes and an inverted repeat
        a = rng.randrange(0, length - 400)
        for i in range(a, a + 300):
            s[i] = s[i].lower()
        b = rng.randrange(0, length - 50)
        for i in range(b, b + rng.randrange(1, 40)):
            s[i] = "N"
        for _ in range(5):
            s[rng.randrange(length)] = rng.choice("RYKMSWn")
        s = "".join(s)
        c = rng.randrange(0, length - 500)
        s = s[: c + 200] + rc(s[c : c + 200]) + s[c + 200 :]
        recs.append((f"ref{r} description text", s))
    recs.append(("short", "ACGTACGTAC"))  # shorter than k
    recs.append(("repeat", recs[0][1][1000:2500]))  # duplicate content
    return recs


def make_queries(rng, refs, k):
    """Query reads cut from the references with random strand, case and single-base changes, plus edge records."""
    qs = []
    seqs = [s for _, s in refs if len(s) > 400]
    for i in range(400):
        s = rng.choice(seqs)
        read_len = rng.randrange(1, 300)
        p = rng.randrange(0, len(s) - read_len)
        q = s[p : p + read_len]
        if rng.random() < 0.5:
            q = rc(q)
        if rng.random() < 0.3:
            q = q.lower()
        if rng.random() < 0.3 and len(q) > 2:
            j = rng.randrange(len(q))
            q = q[:j] + rng.choice("ACGTN") + q[j + 1 :]
        qs.append((f"q{i}", q))
    qs += [
        ("rand", "".join(rng.choice("ACGT") for _ in range(500))),
        ("allN", "N" * 80),
        ("tiny", "ACG"),
        ("exactk", seqs[0][100 : 100 + k]),
        ("iupac", "ACGTRYACGT" * 10),
    ]
    return qs


# Query suite
def suite_query(optik, td, rng):
    """Compare query counts of every index type with the oracle at several k and thread counts."""
    print("== query")
    refs = make_reference(rng)
    ref_fa = os.path.join(td, "ref.fa")
    write_fasta(ref_fa, refs, crlf=True, empty_extra=True)
    for k in (31, 32, 21, 20, 15):
        index = set()
        for _, s in refs:
            index.update(kmers_of(s, k)[1])
        qs = make_queries(rng, refs, k)
        for fmt in ("fa", "fq", "fa_wrapped"):
            qpath = os.path.join(td, f"q_{k}.{fmt}")
            if fmt == "fq":
                write_fastq(qpath, qs)
            else:
                write_fasta(qpath, qs, width=(37 if fmt == "fa_wrapped" else 10**9))
            exp_total = exp_valid = exp_pos = 0
            for _, q in qs:
                t, ks = kmers_of(q, k)
                exp_total += t
                exp_valid += len(ks)
                exp_pos += sum(1 for x in ks if x in index)
            for sub in ("build-spss", "build-fm"):
                if sub == "build-fm" and k > DBG_MAX_K:
                    continue
                idx = os.path.join(td, f"idx_{sub}_{k}")
                b = run([optik, sub, "-i", ref_fa, "-k", str(k), "-o", idx])
                if b.returncode != 0:
                    check(f"query k={k} {sub} build", False, b.stderr.strip()[-200:])
                    continue
                for t in (1, 4):
                    p = run(
                        [
                            optik,
                            "query",
                            "-i",
                            idx,
                            "-q",
                            qpath,
                            "-k",
                            str(k),
                            "-t",
                            str(t),
                        ]
                    )
                    got = {
                        m.group(1): int(m.group(2))
                        for m in re.finditer(
                            r"^(num_\w+) = (\d+)", p.stdout, re.MULTILINE
                        )
                    }
                    want = {
                        "num_kmers": exp_total,
                        "num_positive_kmers": exp_pos,
                        "num_negative_kmers": exp_valid - exp_pos,
                        "num_invalid_kmers": exp_total - exp_valid,
                    }
                    bad = {
                        key: (got.get(key), v)
                        for key, v in want.items()
                        if got.get(key) != v
                    }
                    check(
                        f"query k={k} {sub} {fmt} t={t}",
                        p.returncode == 0 and not bad,
                        str(bad) if bad else "",
                    )


# de Bruijn graph suite
def suite_dbg(optik, td, rng):
    """Check dbruijn's GFA with the validate_unitigs oracle on reads with cycles, palindromes and hairpins."""
    print("== dbg")
    base = "".join(rng.choice("ACGT") for _ in range(3000))
    reads = []
    for i in range(1500):
        read_len = rng.randrange(20, 160)
        p = rng.randrange(0, len(base) - read_len)
        r = base[p : p + read_len]
        if rng.random() < 0.5:
            r = rc(r)
        if rng.random() < 0.2:
            j = rng.randrange(len(r))
            r = r[:j] + rng.choice("ACGT") + r[j + 1 :]
        if rng.random() < 0.1:
            r = r.lower()
        if rng.random() < 0.05:
            j = rng.randrange(len(r))
            r = r[:j] + "N" + r[j + 1 :]
        reads.append((f"r{i}", r))
    # A tandem repeat forming a cycle, a palindrome-rich stretch and a hairpin
    unit = "ACGTTGCAAC"
    reads += [
        ("cycle", unit * 12),
        ("pal", "ACGT" * 30),
        ("hairpin", base[:60] + rc(base[:60])),
    ]
    fq = os.path.join(td, "dbg.fq")
    write_fastq(fq, reads)
    fa = os.path.join(td, "dbg.fa")
    write_fasta(fa, reads, width=50, crlf=True)
    for k in (31, 21, 16, 9, 5):
        for c in (1, 2, 3):
            for inp in (fq, fa):
                gfa = os.path.join(td, f"g_{k}_{c}.gfa")
                p = run(
                    [
                        optik,
                        "dbruijn",
                        "-i",
                        inp,
                        "-k",
                        str(k),
                        "-c",
                        str(c),
                        "-o",
                        gfa,
                        "-t",
                        "4",
                    ]
                )
                tag = f"dbg k={k} c={c} {os.path.basename(inp)}"
                if p.returncode != 0:
                    check(tag, False, p.stderr.strip()[-200:])
                    continue
                res = vu.validate(inp, gfa, k, c)
                bad = {n: v for n, v in res.items() if not v["pass"]}
                check(tag, not bad, str(bad)[:300] if bad else "")


# Mapping suite
def nm_from(seq, cigar, ref, pos1):
    """Edit distance of a SAM alignment recomputed from SEQ, CIGAR and the reference."""
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    qi, ri, nm = 0, pos1 - 1, 0
    for a, o in ops:
        a = int(a)
        if o == "M":
            nm += sum(
                1
                for j in range(a)
                if seq[qi + j].upper() != ref[ri + j].upper()
                or ref[ri + j].upper() not in "ACGT"
            )
            qi += a
            ri += a
        elif o in "I":
            nm += a
            qi += a
        elif o == "D":
            nm += a
            ri += a
        elif o == "S":
            qi += a
    return nm


def suite_map(optik, td, rng):
    """Map reads at sequence ends, with N and across references, checking origin, SAM validity and parity."""
    print("== map")
    refs = [
        (f"chr{i}", "".join(rng.choice("ACGT") for _ in range(200000)))
        for i in range(3)
    ]
    refs.append(("tiny", "ACGTACGTACGTTTGACCA"))
    ref_fa = os.path.join(td, "mref.fa")
    write_fasta(ref_fa, [(n, s[:5000].lower() + s[5000:]) for n, s in refs], crlf=True)
    ref = {n: s.upper() for n, s in refs}
    truth, reads = {}, []
    for i in range(3000):
        name, s = rng.choice(refs[:3])
        read_len = rng.choice([30, 60, 100, 150, 151, 300, 800, 2500])
        # Include reads at both ends of a reference.
        p = rng.choice([0, len(s) - read_len, rng.randrange(0, len(s) - read_len)])
        r = list(s[p : p + read_len])
        for j in range(len(r)):
            if rng.random() < 0.01:
                r[j] = rng.choice("ACGT")
        if rng.random() < 0.1:
            r[rng.randrange(read_len)] = "N"
        r = "".join(r)
        strand = "-" if rng.random() < 0.5 else "+"
        if strand == "-":
            r = rc(r)
        if rng.random() < 0.1:
            r = r.lower()
        truth[f"m{i}"] = (name, p, strand, read_len)
        reads.append((f"m{i}", r))
    reads += [
        ("short", "ACGTAC"),
        ("allN", "N" * 200),
        ("empty", ""),
        ("rand", "".join(rng.choice("ACGT") for _ in range(150))),
    ]
    for fmt in ("fq", "fa"):
        qpath = os.path.join(td, f"mreads.{fmt}")
        (write_fastq if fmt == "fq" else write_fasta)(qpath, reads)
        sam = os.path.join(td, f"m_{fmt}.sam")
        p = run([optik, "map", "-r", ref_fa, "-q", qpath, "-o", sam, "-t", "4"])
        if p.returncode != 0:
            check(f"map {fmt} exit", False, p.stderr.strip()[-300:])
            continue
        prim = Counter()
        unmapped_names = []
        bad_len = bad_nm = bad_seq = wrong = unmapped = bad_unmapped = 0
        examples = []
        readseq = dict(reads)
        with open(sam) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                t = line.rstrip("\n").split("\t")
                if len(t) < SAM_MIN_COLUMNS:
                    bad_unmapped += 1
                    continue
                f = int(t[1])
                if not f & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                    prim[t[0]] += 1
                if f & SAM_FLAG_UNMAPPED:
                    if t[2] != "*" or t[5] != "*":
                        bad_unmapped += 1
                    if t[0] in truth and not f & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                        unmapped += 1
                        unmapped_names.append(t[0])
                    continue
                seq = t[9]
                want = readseq[t[0]]
                # A mapped record must carry SEQ and NM:i, so a writer that drops them cannot pass.
                if (
                    seq == "*"
                    or seq.upper()
                    != (rc(want) if f & SAM_FLAG_REVERSE else want).upper()
                ):
                    bad_seq += 1
                ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
                cig_ok = seq != "*" and sum(
                    int(a) for a, o in ops if o in "MIS=X"
                ) == len(seq)
                if not cig_ok:
                    bad_len += 1
                nm = next((int(x[5:]) for x in t[11:] if x.startswith("NM:i:")), None)
                # NM can be recomputed only when the CIGAR consumes SEQ.
                if (
                    nm is None
                    or not cig_ok
                    or nm != nm_from(seq, t[5], ref[t[2]], int(t[3]))
                ):
                    bad_nm += 1
                if t[0] in truth and not f & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                    name, p0, strand, read_len = truth[t[0]]
                    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                    ok = (
                        t[2] == name
                        and ((f & SAM_FLAG_REVERSE) > 0) == (strand == "-")
                        and abs(int(t[3]) - 1 - lead - p0) <= max(10, read_len // 10)
                    )
                    if not ok:
                        win = ref[t[2]][int(t[3]) - 1 : int(t[3]) - 1 + read_len]
                        ok = win == ref[name][p0 : p0 + read_len] and nm == 0
                    if not ok:
                        wrong += 1
                        if len(examples) < 3:
                            examples.append(
                                (t[0], truth[t[0]], t[1], t[2], t[3], t[5][:30])
                            )
        missing = [r for r, _ in reads if prim[r] == 0]
        dup = [r for r, n in prim.items() if n > 1]
        check(
            f"map {fmt}: one primary record per read",
            not missing and not dup,
            f"missing={missing[:5]} dup={dup[:5]}",
        )
        check(f"map {fmt}: SEQ orientation", bad_seq == 0, f"bad={bad_seq}")
        check(f"map {fmt}: CIGAR consumes SEQ", bad_len == 0, f"bad={bad_len}")
        check(f"map {fmt}: NM == CIGAR vs reference", bad_nm == 0, f"bad={bad_nm}")
        check(
            f"map {fmt}: unmapped records well-formed",
            bad_unmapped == 0,
            f"bad={bad_unmapped}",
        )
        n_long = sum(1 for v in truth.values() if v[3] >= 60)
        check(
            f"map {fmt}: placements correct",
            wrong == 0,
            f"wrong={wrong} unmapped={unmapped}/{len(truth)} (>=60bp: {n_long}) e.g. {examples}",
        )
        by_len = Counter(truth[r][3] for r in unmapped_names)
        mm2 = os.path.join(HERE, "..", "temp/minimap2/minimap2")
        if fmt == "fq" and os.path.isfile(mm2):
            # Sensitivity parity with the same reads and truth rule as minimap2 -x sr/map-ont
            o_ok = len(truth) - wrong - unmapped
            m_ok = 0
            for preset, lo, hi in (("sr", 0, 1000), ("map-ont", 1000, 10**9)):
                sub = os.path.join(td, f"mm_{preset}.fq")
                write_fastq(
                    sub, [(n, s) for n, s in reads if n in truth and lo <= len(s) < hi]
                )
                pm = run([mm2, "-a", "-x", preset, ref_fa, sub])
                seen = set()
                for line in pm.stdout.splitlines():
                    if line.startswith("@"):
                        continue
                    t = line.split("\t")
                    f = int(t[1])
                    if f & SAM_FLAGS_NOT_PRIMARY or t[0] in seen:
                        continue
                    seen.add(t[0])
                    name, p0, strand, read_len = truth[t[0]]
                    ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
                    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                    m_ok += (
                        t[2] == name
                        and ((f & SAM_FLAG_REVERSE) > 0) == (strand == "-")
                        and abs(int(t[3]) - 1 - lead - p0) <= max(10, read_len // 10)
                    )
            check(
                f"map {fmt}: correct placements >= minimap2",
                o_ok >= m_ok,
                f"optik {o_ok} vs minimap2 {m_ok} of {len(truth)}; "
                f"optik unmapped by length={dict(by_len)}",
            )


def suite_memory(optik, td, rng):
    """Reads with off-reference flanks.
    Before flank extension was chunked, one 15-40 kb flank cost GBs of wavefront memory per thread, which exhausted a 32 GB host."""
    print("== memory")
    ref = "".join(rng.choice("ACGT") for _ in range(300000))
    ref_fa = os.path.join(td, "memref.fa")
    write_fasta(ref_fa, [("ref", ref)])
    truth, reads = {}, []
    for i in range(48):
        read_len = rng.randrange(3000, 5000)
        p0 = rng.randrange(0, len(ref) - read_len)
        core = list(ref[p0 : p0 + read_len])
        for j in range(read_len):
            if rng.random() < 0.05:
                core[j] = rng.choice("ACGT")
        left = "".join(rng.choice("ACGT") for _ in range(rng.randrange(0, 20000)))
        right = "".join(rng.choice("ACGT") for _ in range(rng.randrange(10000, 35000)))
        r = left + "".join(core) + right
        strand = "-" if rng.random() < 0.5 else "+"
        reads.append((f"m{i}", rc(r) if strand == "-" else r))
        truth[f"m{i}"] = (p0, strand)
    # An unrelated internal 4.9 kb block between two genuine parts of 5-6 kb must yield a primary and a supplementary record, never one record aligning the block.
    for i in range(16):
        p0 = rng.randrange(0, len(ref) - 20000)
        r = (
            ref[p0 : p0 + 6000]
            + "".join(rng.choice("ACGT") for _ in range(4900))
            + ref[p0 + 10900 : p0 + 16000]
        )
        reads.append((f"u{i}", r))
        truth[f"u{i}"] = (p0, "+")
    rd_fa = os.path.join(td, "memreads.fa")
    write_fasta(rd_fa, reads)
    # The gate is the child's measured peak RSS.
    # The address-space limit only protects the host, because under it the aligner's allocation failure is caught and the flank soft-clipped, so the limit alone cannot reveal the bug.
    cap = 6 << 30
    sam = os.path.join(td, "mem.sam")
    with open(sam, "w") as out, open(os.path.join(td, "mem.err"), "w") as err:
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
        _, status, ru = os.wait4(proc.pid, 0)
    peak_mb = ru.ru_maxrss / 1024
    check(
        "memory: map exits 0",
        os.waitstatus_to_exitcode(status) == 0,
        f"exit {os.waitstatus_to_exitcode(status)}",
    )
    check(
        "memory: peak RSS <= 1 GiB at 4 threads", peak_mb <= 1024, f"{peak_mb:.0f} MB"
    )
    p = subprocess.CompletedProcess([], 0, Path(sam).read_text(), "")
    bad = []
    ok = 0
    for line in p.stdout.splitlines():
        if line.startswith("@"):
            continue
        t = line.split("\t")
        f = int(t[1])
        if f & SAM_FLAGS_NOT_PRIMARY or t[0] not in truth:
            continue
        p0, strand = truth[t[0]]
        good = ((f & SAM_FLAG_REVERSE) > 0) == (strand == "-") and abs(
            int(t[3]) - 1 - p0
        ) <= 100
        ok += good
        if not good:
            bad.append(
                f"{t[0]} flag={f} pos={t[3]} want={p0 + 1}{strand} cigar={t[5][:40]}"
            )
    check(
        "memory: primaries on the reference-derived core",
        ok == len(truth),
        f"{ok}/{len(truth)} {bad[:3]}",
    )
    # No garbage, meaning every aligned record has NM of at most 20% of its aligned bases
    worst = []
    for line in p.stdout.splitlines():
        if line.startswith("@"):
            continue
        t = line.split("\t")
        if int(t[1]) & SAM_FLAG_UNMAPPED:
            continue
        aligned = sum(int(n) for n, o in re.findall(r"(\d+)([MI])", t[5]))
        nm = int(next(x for x in t[11:] if x.startswith("NM:i:"))[5:])
        if nm > 0.2 * aligned:
            worst.append(f"{t[0]} {t[5][:30]} NM={nm}")
    check(
        "memory: no record aligns unrelated sequence (NM <= 20% of aligned)",
        not worst,
        f"{len(worst)} records: {worst[:3]}",
    )
    supp = sum(
        1
        for line in p.stdout.splitlines()
        if not line.startswith("@")
        and line.split("\t")[0].startswith("u")
        and int(line.split("\t")[1]) & SAM_FLAG_SUPPLEMENTARY
    )
    check(
        "memory: unrelated internal block splits into primary + supplementary",
        supp == 16,
        f"{supp}/16 supplementaries",
    )


def suite_memory_bounds(optik, td, rng):
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
        p = subprocess.run(
            [gnu_time, "-f", "%M", *cmd],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        last = (p.stderr.strip().splitlines() or ["0"])[-1]
        return p.returncode, (int(last) / 1024 if last.isdigit() else float("inf"))

    d = os.path.join(td, "membound")
    os.makedirs(d, exist_ok=True)
    genome = "".join(rng.choice("ACGT") for _ in range(3_000_000))
    gfa = os.path.join(d, "g.fa")
    write_fasta(gfa, [("g", genome)], width=10**9)  # a single long line per record
    idx = os.path.join(d, "g.spss")
    run([optik, "build-spss", "-i", gfa, "-k", "31", "-o", idx])
    q = os.path.join(d, "rand.fa")
    with open(q, "w") as fh:
        fh.writelines(
            f">{i}\n{''.join(rng.choice('ACGT') for _ in range(31))}\n"
            for i in range(1_000_000)
        )
    with open(q, "rb") as fi, _gz.open(q + ".gz", "wb", compresslevel=1) as fo:
        fo.write(fi.read())
    for path in (q, q + ".gz"):
        for t in ("1", "8"):
            rc_, mb = peak([optik, "query", "-i", idx, "-q", path, "-k", "31", "-t", t])
            check(
                f"memory: query 1M records ({os.path.basename(path)}, -t {t}) peak <= 40 MB",
                rc_ == 0 and mb <= 40,
                f"{mb:.0f} MB",
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
                os.path.join(d, "g.fm"),
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
                os.path.join(d, "g.gfa"),
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
                os.path.join(d, "gc.spss"),
            ],
        ),
    ):
        rc_, mb = peak(cmd)
        per = (mb - 10) * 1048576 / n_kmers
        check(
            f"memory: {name} peak <= 60 bytes per k-mer (single 3 Mbp record)",
            rc_ == 0 and per <= 60,
            f"{mb:.0f} MB = {per:.0f} B/k-mer",
        )
    # Low complexity, 8 Mbp holding 3 distinct k-mers.
    # The counter's per-partition table was once sized by occurrences, about 270 MB per pass here, instead of growing with the distinct k-mers it holds.
    lc = os.path.join(d, "lowc.fa")
    write_fasta(lc, [("a", "A" * 4_000_000), ("c", "CG" * 2_000_000)], width=10**9)
    for name, cmd in (
        (
            "dbruijn",
            [
                optik,
                "dbruijn",
                "-i",
                lc,
                "-k",
                "31",
                "-t",
                "8",
                "-o",
                os.path.join(d, "lc.gfa"),
            ],
        ),
        (
            "build-fm",
            [
                optik,
                "build-fm",
                "-i",
                lc,
                "-k",
                "31",
                "-t",
                "8",
                "-o",
                os.path.join(d, "lc.fm"),
            ],
        ),
    ):
        rc_, mb = peak(cmd)
        check(
            f"memory: {name} on 8 Mbp of poly-A / CG repeat (3 distinct k-mers) peak <= 100 MB",
            rc_ == 0 and mb <= 100,
            f"{mb:.0f} MB",
        )


def suite_corner(optik, td, rng):
    """Input-format, geometry and parameter corner cases.

    Each was found by a probe that showed wrong or silently partial output before its fix."""
    import gzip
    import shutil

    optik_bin = optik
    d = os.path.join(td, "corner")
    os.makedirs(d, exist_ok=True)

    def in_dir(name):
        """Path of `name` in this suite's directory."""
        return os.path.join(d, name)

    def random_dna(n):
        """Random ACGT string of length n."""
        return "".join(rng.choice("ACGT") for _ in range(n))

    genome, genome2 = random_dna(100000), random_dna(30000)
    write_fasta(in_dir("ref.fa"), [("chr1", genome), ("chr2", genome2)], width=10**9)
    reads = [
        (f"r{i}", genome[p : p + 150])
        for i, p in enumerate(rng.sample(range(99000), 200))
    ]
    reads = [(n, rc(s) if i % 2 else s) for i, (n, s) in enumerate(reads)]
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
    want = 200 * (150 - 31 + 1)

    def counts_of(p):
        """The num_* counters printed by an optik query run."""
        return {
            m.group(1): int(m.group(2))
            for m in re.finditer(r"^(num_\w+) = (\d+)", p.stdout, re.MULTILINE)
        }

    def q(path, **kw):
        """Query the SPSS index with `path`, returning (exit code, num_kmers, num_positive_kmers)."""
        p = run(
            [optik_bin, "query", "-i", in_dir("ref.spss"), "-q", path, "-k", "31"], **kw
        )
        c = counts_of(p)
        return p.returncode, c.get("num_kmers"), c.get("num_positive_kmers")

    def all_found(name, path, **kw):
        """Check that querying `path` finds all `want` k-mers."""
        rcode, nk, pos = q(path, **kw)
        check(
            f"corner: {name}",
            rcode == 0 and nk == want and pos == want,
            f"exit {rcode}, {pos}/{nk} positive (want {want}/{want})",
        )

    def sam(path):
        """SAM records of a file as field lists, header excluded."""
        with open(path) as fh:
            return [
                line.rstrip("\n").split("\t") for line in fh if not line.startswith("@")
            ]

    print("== corner: input formats")
    half = len(reads) // 2
    with open(in_dir("mm.fa.gz"), "wb") as fh:  # multi-member gzip, as bgzip writes
        fh.writelines(
            gzip.compress("".join(f">{n}\n{s}\n" for n, s in part).encode())
            for part in (reads[:half], reads[half:])
        )
    all_found("multi-member gzip", in_dir("mm.fa.gz"))
    shutil.copy(in_dir("mm.fa.gz"), in_dir("mm_noext.fa"))
    all_found("gzip detected by magic bytes (no .gz suffix)", in_dir("mm_noext.fa"))
    all_found("query from stdin '-'", "-", input=Path(in_dir("reads.fa")).read_text())
    with open(in_dir("ml.fq"), "w") as fh:
        fh.writelines(
            f"@{n}\n{s[:70]}\n{s[70:]}\n+\n{'I' * 70}\n{'I' * 80}\n" for n, s in reads
        )
    all_found("multi-line FASTQ", in_dir("ml.fq"))
    with open(in_dir("crlf.fq"), "w", newline="") as fh:
        fh.writelines(f"@{n}\r\n{s}\r\n+\r\n{'I' * 150}\r\n" for n, s in reads)
    all_found("CRLF FASTQ", in_dir("crlf.fq"))
    with open(in_dir("qat.fq"), "w") as fh:
        fh.writelines(f"@{n}\n{s}\n+\n@{'I' * 149}\n" for n, s in reads)
    all_found("FASTQ quality line starting with '@'", in_dir("qat.fq"))
    Path(in_dir("noeol.fa")).write_text(
        "".join(f">{n}\n{s}\n" for n, s in reads).rstrip("\n")
    )
    all_found("FASTA without a final newline", in_dir("noeol.fa"))
    Path(in_dir("blank.fa")).write_text(
        "".join(f">{n}\n{s[:75]}\n\n{s[75:]}\n" for n, s in reads)
    )
    all_found("blank line inside a record", in_dir("blank.fa"))
    Path(in_dir("ws.fa")).write_text(
        "".join(f">{n} desc\t\n{s[:75]} \t\n{s[75:]}\t \n" for n, s in reads)
    )
    all_found("trailing spaces/tabs on wrapped sequence lines", in_dir("ws.fa"))
    Path(in_dir("rna.fa")).write_text(
        "".join(f">{n}\n{s.replace('T', 'U')}\n" for n, s in reads)
    )
    all_found("RNA (U read as T)", in_dir("rna.fa"))

    p = run(
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
    clean = all(re.fullmatch(r"[ACGTN*]+", t[9]) for t in recs)
    check(
        "corner: map SEQ carries no whitespace from the input lines",
        p.returncode == 0 and clean and len(recs) == len(reads),
    )
    p = run(
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
    mapped = sum(1 for t in sam(in_dir("rna.sam")) if not int(t[1]) & SAM_FLAG_UNMAPPED)
    check(
        "corner: map RNA reads (U as T)",
        p.returncode == 0 and mapped == len(reads),
        f"{mapped}/{len(reads)} mapped",
    )

    # One record of several hundred kb is cut into window-exact slices when -t > 1, and every count must equal the oracle at every thread count and index type.
    big = (
        genome[: 300000 % len(genome)]
        + "N" * 50
        + rc(genome)
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
        ("big2", rc(big[:200000])),
    ]
    write_fasta(in_dir("bigq.fa"), recs, width=10**9)
    want_big = {
        "num_kmers": 0,
        "num_positive_kmers": 0,
        "num_negative_kmers": 0,
        "num_invalid_kmers": 0,
    }
    for _, sq in recs:
        tot, ks = kmers_of(sq, 31)
        pos = sum(1 for x in ks if x in index)
        want_big["num_kmers"] += tot
        want_big["num_positive_kmers"] += pos
        want_big["num_negative_kmers"] += len(ks) - pos
        want_big["num_invalid_kmers"] += tot - len(ks)
    for sub in ("build-spss", "build-fm"):
        idx = in_dir(f"bigidx_{sub}")
        run([optik_bin, sub, "-i", in_dir("ref.fa"), "-k", "31", "-o", idx])
        for t in (1, 2, 8, 16):
            c = counts_of(
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
                        str(t),
                    ]
                )
            )
            bad = {x: (c.get(x), v) for x, v in want_big.items() if c.get(x) != v}
            check(
                f"corner: long record split across threads, {sub} -t {t}: counts == oracle",
                not bad,
                str(bad),
            )

    print("== corner: read names")
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
    p = run(
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
    names = [t[0] for t in sam(in_dir("names.sam"))]
    check(
        "corner: QNAME > 254 truncated to the SAM limit (with a warning)",
        p.returncode == 0 and names[0] == "x" * 254 and "truncated" in p.stderr,
        f"{len(names[0])} chars",
    )
    check(
        "corner: header words after the name dropped; duplicate names kept",
        names[1:] == ["dup", "dup"],
        str(names[1:]),
    )
    if shutil.which("samtools"):
        sv = run(["samtools", "view", "-c", in_dir("names.sam")])
        check(
            "corner: samtools reads the SAM with a long read name",
            sv.returncode == 0,
            sv.stderr.strip()[-120:],
        )

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
    p = run(
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
    by = {}
    for t in sam(in_dir("geo.sam")):
        by.setdefault(t[0], []).append(t)

    def recs_of(name):
        """(flag, ref, pos, CIGAR) of every record of read `name`."""
        return [(int(t[1]), t[2], int(t[3]), t[5]) for t in by.get(name, [])]

    def split_ok(name, parts):
        """Check for one primary and one supplementary with these (ref, pos, CIGAR) parts.

        Each must carry an SA tag, and which half is primary is a tie.
        Clips compare as S, while the primary must soft-clip and the supplementary hard-clip, the convention of minimap2 and BWA-MEM, whose supplementary SEQ holds only the aligned bases."""
        r = recs_of(name)
        flags = sorted(f for f, *_ in r)
        sa = all(any(x.startswith("SA:Z:") for x in t[11:]) for t in by.get(name, []))
        clips = all(
            ("H" not in c)
            if not f & SAM_FLAG_SUPPLEMENTARY
            else ("S" not in c and "H" in c)
            for f, _, _, c in r
        )
        norm = sorted((x[1], x[2], x[3].replace("H", "S")) for x in r)
        return flags == [0, 2048] and norm == sorted(parts) and sa and clips

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
    for nm, pos, ln in (("b499", 1001, 499), ("b500", 2001, 500), ("b501", 3001, 501)):
        check(
            f"corner: {ln} bp read at the short/long routing boundary",
            recs_of(nm) == [(0, "chr1", pos, f"{ln}M")],
            str(recs_of(nm)),
        )
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
    p = run(
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
    by = {}
    for t in sam(in_dir("gap.sam")):
        by.setdefault(t[0], []).append(t)
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
        [r[0] & SAM_FLAG_UNMAPPED for r in recs_of("in_N")] == [SAM_FLAG_UNMAPPED],
        str(recs_of("in_N")),
    )

    print("== corner: map parameter limits")
    for extra in (
        ["-e", "0"],
        ["-w", "1"],
        ["-k", "32"],
        ["-k", "8"],
        ["-k", "32", "-w", "1"],
    ):
        p = run(
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
        ok = (
            [t for t in sam(in_dir("par.sam")) if not int(t[1]) & SAM_FLAG_UNMAPPED]
            if p.returncode == 0
            else []
        )
        at_origin = sum(1 for t in ok if genome.find(t[9]) + 1 == int(t[3]))
        check(
            f"corner: map {' '.join(extra)}: every exact read at its origin",
            p.returncode == 0 and at_origin == len(reads),
            f"{at_origin}/{len(reads)}",
        )

    print("== corner: abundance and k-mer ids")
    # Two k-mers above the graph's in-slot count range (2^22 - 1), the poly-A k-mer and the CG-repeat k-mer (whose two strands share one canonical k-mer)
    Path(in_dir("homo.fa")).write_text(
        ">h\n" + "A" * 17_000_000 + "\n>g\n" + "CG" * 2_200_000 + "\n"
    )
    p = run(
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
    kc = (
        re.findall(r"KC:i:(\d+)", Path(in_dir("homo.gfa")).read_text())
        if p.returncode == 0
        else []
    )
    check(
        "corner: dbruijn abundances above 2^22 and 2^24 exact (no saturation)",
        sorted(map(int, kc)) == [4_400_000 - 30, 17_000_000 - 30],
        f"KC {kc}",
    )
    write_fasta(
        in_dir("dupk.fa"),
        [("a", genome[:5000]), ("b", genome[:5000]), ("c", rc(genome[:5000]))],
        width=10**9,
    )
    p = run(
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
    dk = re.search(r"distinct_kmers\t(\d+)", p.stdout)
    check(
        "corner: build-spss distinct_kmers counts distinct k-mers, not positions",
        dk is not None and int(dk.group(1)) == 5000 - 30,
        f"{dk.group(1) if dk else None}",
    )
    check(
        "corner: build-spss warns that duplicate k-mers make --ids a non-bijection",
        "not a bijection" in p.stderr,
    )
    write_fasta(
        in_dir("idq.fa"), [("f", genome[:5000]), ("r", rc(genome[:5000]))], width=10**9
    )
    p = run(
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
    if p.returncode == 0:
        with open(in_dir("ids.txt")) as fh:
            lines = [list(map(int, line.split())) for line in fh]
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
    def palindrome(n):
        """A random n-base sequence followed by its reverse complement."""
        y = random_dna(n)
        return y + rc(y)

    hp = "".join(random_dna(40) + palindrome(18) for _ in range(200)) + random_dna(40)
    write_fasta(in_dir("hairpin.fa"), [("h", hp)], width=10**9)
    p = run(
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
    p2 = (
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
        if p.returncode == 0
        else p
    )
    ids = (
        list(map(int, Path(in_dir("hp.txt")).read_text().split()))
        if p2.returncode == 0
        else []
    )
    by = {}
    for i, km in enumerate(
        kmers_of(hp, 31)[1]
    ):  # canonical, at every position of an N-free sequence
        by.setdefault(km, set()).add(ids[i] if i < len(ids) else -1)
    dk = re.search(r"distinct_kmers\t(\d+)", p.stdout) if p.returncode == 0 else None
    check(
        "corner: hairpin k-mers: one id per distinct k-mer on both strands, distinct count exact",
        len(ids) == len(hp) - 30
        and all(len(v) == 1 for v in by.values())
        and len({next(iter(v)) for v in by.values()}) == len(by)
        and dk is not None
        and int(dk.group(1)) == len(by),
        f"{sum(len(v) > 1 for v in by.values())} multi-id, distinct {dk.group(1) if dk else None} vs {len(by)}",
    )

    print("== corner: spss-update")
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
        p = run(
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
        ka = re.search(r"kmers_after\t(\d+)", p.stdout)
        c = counts_of(
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
            p.returncode == 0
            and ka is not None
            and int(ka.group(1)) == after
            and c.get("num_positive_kmers") == pos,
            f"kmers_after {ka.group(1) if ka else None} (want {after}), positives {c.get('num_positive_kmers')} (want {pos})",
        )
    shutil.copy(in_dir("ref.spss"), in_dir("inplace.spss"))
    write_fasta(in_dir("add.fa"), [("a", random_dna(5000))])
    p = run(
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
    c = counts_of(
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
        p.returncode == 0 and c.get("num_positive_kmers") == 5000 - 30,
        str(c),
    )


def main():
    """Run the suites chosen with --only."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    ap.add_argument("--only", default="query,dbg,map,memory,membound,corner")
    ap.add_argument("-s", "--seed", type=int, default=11)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    with tempfile.TemporaryDirectory(prefix="optik_edge_") as td:
        for name in args.only.split(","):
            {
                "query": suite_query,
                "dbg": suite_dbg,
                "map": suite_map,
                "memory": suite_memory,
                "membound": suite_memory_bounds,
                "corner": suite_corner,
            }[name](args.optik, td, rng)
    print(f"edge_cases: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
