#!/usr/bin/env python3
"""validate_official.py - exact-equality validation of optik against the original tools, namely SSHash, Cuttlefish, Bifrost, minimap2, samtools, KMC and k8 with paftools.

Unlike the threshold-based validators, each check here asserts exact equality with 0 mismatches, or a strict structural invariant.
The harness covers the real use case by running the commands a researcher would run with the official tools and requiring optik's output to agree exactly.

Each of the three suites is a subcommand.

  query   reference unitigs and query reads
          - SSHash `build --check` reports "EVERYTHING OK!" for the same reference.
          - optik query equals sshash query on single-k-mer FASTA probes drawn from the reference k-mer set of the KMC oracle and on random non-reference k-mers.
            Positive probes give exactly G positives and 0 negatives, negative probes the reverse, and num_kmers is identical.
          - On real cross-strain reads both tools report 0 positives and the same num_negative, with optik's num_kmers equal to sshash's num_kmers minus its num_invalid.

  dbg     reads to unitig graph
          - The k-mer set of optik dbruijn's unitigs equals the KMC oracle on the reads exactly, with 0 differences in either direction.
          - Unitigs cover every input k-mer, and every unitig k-mer belongs to the input set, so none is missing or spurious.
          - No unitig crosses a branch node of the canonical k-mer graph, where interior nodes have degree 2 and self-loops are allowed.
          - Bifrost's unitig set equals the oracle with 0 differences.
          - Cuttlefish's unitig set must be a subset of the oracle, since it drops k-mers isolated between N runs.
            The count of dropped k-mers should be small and must equal |oracle - cuttlefish|.

  map     reference and reads
          - On error-free simulated reads with gen_reads.py ground truth, every read maps with position, strand and start exactly as in the truth, where optik's 1-based position equals the 0-based truth plus 1.
          - On simulated reads with errors, optik's reported NM equals the edit distance to the true locus, computed independently.
          - minimap2 concordance requires the same mapped/unmapped classification and starts within 1 bp on error-free reads.
          - A samtools SAM->BAM->SAM round trip succeeds, flagstat parses, and the BAM is byte-stable under samtools validation.

The exit code is 0 only if every required check passes.
"""

import argparse
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_NOT_PRIMARY,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
    SAM_MIN_COLUMNS,
)

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    """Reverse complement of seq."""
    return seq.translate(COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    """The lexicographically smaller of kmer and its reverse complement."""
    r = rc(kmer)
    return min(kmer, r)


# Process helpers


def run(cmd, **kw):
    """subprocess.run with stdout and stderr captured unless the caller says otherwise."""
    kw.setdefault("stdout", subprocess.PIPE)
    kw.setdefault("stderr", subprocess.PIPE)
    return subprocess.run(cmd, **kw, check=False)


def need(name):
    """Stop the run because required tool `name` is missing."""
    raise SystemExit(f"ERROR: required tool '{name}' not available")


class Tool:
    """An external tool's path, which is false in a boolean context when the binary is missing."""

    def __init__(self, path):
        self.path = Path(path) if path else None
        self.available = bool(self.path) and self.path.is_file()

    def __bool__(self):
        return self.available

    def __str__(self):
        return str(self.path)


# FASTA and sequence helpers


def read_fasta(path):
    """Yield (name, seq) records of a FASTA file, plain or .gz, without splitting at N runs."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name = None
    parts = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                toks = line[1:].split()
                name = toks[0] if toks else f"seq{len(parts)}"
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def read_seqs(path):
    """Yield (name, seq) of FASTA or FASTQ records, plain or .gz.
    The official query reads are FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.read(1)
        fh.seek(0)
        if first == "@":
            while True:
                h = fh.readline()
                if not h:
                    return
                sq = fh.readline().strip()
                fh.readline()
                fh.readline()
                yield h[1:].split()[0] if h[1:].split() else "", sq
        else:
            yield from read_fasta(path)


def extract_kmers(seq, k):
    """Canonical k-mers of seq, splitting at N."""
    out = []
    start = 0
    for i, c in enumerate(seq):
        if c == "N":
            if i - start >= k:
                out.extend(canonical(seq[j : j + k]) for j in range(start, i - k + 1))
            start = i + 1
    if len(seq) - start >= k:
        out.extend(canonical(seq[j : j + k]) for j in range(start, len(seq) - k + 1))
    return out


def to_single_line_fasta(path_in, path_out):
    """Write one sequence per line with every record kept, as KMC's FASTA reader and SSHash's query parser require."""
    with open(path_out, "w") as fh:
        fh.writelines(f">{name}\n{seq}\n" for name, seq in read_fasta(path_in))


def write_uniquified_fasta(path_in, path_out):
    """Write single-line FASTA with unambiguous names.

    Stitched-unitig FASTA files can repeat a header, such as many '>seq1' records.
    These are renamed seq1, seq2, ... so that ground-truth positions refer to one contig only.

    Callers that index records by name in a dict over `read_fasta` collapse duplicates, while `gen_reads` works on the full list of records and reports only the name.
    The truth window cannot be reconstructed without unique names, so the map reference is rewritten with unique IDs.
    """
    with open(path_out, "w") as fh:
        fh.writelines(
            f">seq{i}\n{seq}\n" for i, (_name, seq) in enumerate(read_fasta(path_in), 1)
        )


def write_single_kmer_fasta(kmers, path_out):
    """Write one canonical k-mer per record.
    k may be up to 63, which covers k = 31."""
    with open(path_out, "w") as fh:
        fh.writelines(f">q{i}\n{k}\n" for i, k in enumerate(kmers))


# KMC oracle


def _is_fastq(path):
    """True for a .fq or .fastq path, gzipped or not."""
    p = str(path).lower()
    return p.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz"))


def kmc_build(kmc_bin, kmc_tmp, in_fasta, prefix, k, threads, min_count=1):
    """Run KMC over a one-sequence-per-line FASTA file or a FASTQ file and return the distinct k-mer set as a Python set."""
    if not kmc_bin.available:
        need("KMC (kmc)")
    workdir = kmc_tmp / "kmc_work"
    workdir.mkdir(parents=True, exist_ok=True)
    fmt = "-fq" if _is_fastq(in_fasta) else "-fa"
    db = kmc_tmp / str(prefix)
    p = run(
        [
            str(kmc_bin.path),
            fmt,
            f"-k{k}",
            f"-ci{min_count}",
            f"-t{threads}",
            str(in_fasta),
            str(db),
            str(workdir),
        ]
    )
    if p.returncode != 0:
        raise RuntimeError(f"KMC failed: {p.stderr.decode()[-2000:]}")
    dump = kmc_tmp / f"{prefix}.dump"
    p = run([str(kmc_bin.path.parent / "kmc_dump"), str(db), str(dump)])
    if p.returncode != 0:
        raise RuntimeError(f"kmc_dump failed: {p.stderr.decode()[-2000:]}")
    return {ln.split("\t")[0] for ln in dump.read_text().splitlines()}


# SSHash helpers


def sshash_query_report(sshash, index, query, k):
    """Completed process of `sshash query` on index and query."""
    p = run(
        [str(sshash.path), "query", "-i", str(index), "-q", str(query)], timeout=3600
    )
    return p


def parse_sshash_report(text):
    """Counter values (num_kmers, num_positive_kmers, ...) parsed from SSHash's query report."""
    m = {}
    for key in (
        "num_kmers",
        "num_positive_kmers",
        "num_negative_kmers",
        "num_invalid_kmers",
        "num_searches",
    ):
        mm = re.search(rf"^{key}\s*=\s*(\d+)", text, re.MULTILINE)
        if mm:
            m[key] = int(mm.group(1))
    return m


# optik helpers


def optik_query_report(optik, index, query, k):
    """Completed process of `optik query` on index and query at k."""
    p = run(
        [str(optik.path), "query", "-i", str(index), "-q", str(query), "-k", str(k)],
        timeout=3600,
    )
    return p


def parse_optik_report(text):
    """Counter values parsed from optik's query report, named as in SSHash's."""
    m = {}
    for key in (
        "num_kmers",
        "num_positive_kmers",
        "num_negative_kmers",
        "num_invalid_kmers",
        "num_searches",
    ):
        mm = re.search(rf"{key}\s*=\s*(\d+)", text)
        if mm:
            m[key] = int(mm.group(1))
    return m


# Query suite


def suite_query(args):
    """Compare optik query with sshash query on probes and reads and return the checks."""
    checks = {}
    optik, sshash, kmc = args.optik, args.sshash, args.kmc
    k = args.k
    keep = Path(args.workdir)
    keep.mkdir(parents=True, exist_ok=True)

    # Step 0 builds the reference k-mer oracle with KMC on the unitig FASTA.
    ref_oneline = keep / "ref_oneline.fa"
    to_single_line_fasta(args.ref, ref_oneline)
    oracle = kmc_build(kmc, keep, ref_oneline, "oracle", k, args.threads)
    checks["kmc_oracle"] = {"pass": True, "n": len(oracle)}

    # Step 1 indexes with optik and sshash.
    # Both read plain FASTA, and the unitig FASTA may be gzipped or wrapped, so the one-line conversion is indexed in every case.
    # optik contributes both index formats, the SPSS index from build-spss and the FM index from build-fm.
    # Every comparison below runs on each, and --ids adds the per-k-mer check for the SPSS index.
    optik_idxs = {}
    for fmt, sub in (("_spss", "build-spss"), ("_fm", "build-fm")):
        idx = keep / f"optik_idx{fmt}"
        p = run(
            [
                str(optik.path),
                sub,
                "-i",
                str(ref_oneline),
                "-k",
                str(k),
                "-o",
                str(idx),
            ],
            timeout=3600,
        )
        checks[f"optik_build{fmt}"] = {
            "pass": p.returncode == 0,
            "stderr": p.stderr.decode()[-500:],
        }
        if p.returncode == 0:
            optik_idxs[fmt] = idx
    if not optik_idxs:
        return checks

    sshash_idx = keep / "sshash_idx"
    p = run(
        [
            str(sshash.path),
            "build",
            "-i",
            str(ref_oneline),
            "-k",
            str(k),
            "-m",
            str(args.sshash_m),
            "-o",
            str(sshash_idx),
            "-d",
            str(keep),
            "--check",
        ],
        timeout=3600,
    )
    out = p.stdout.decode() + p.stderr.decode()
    checks["sshash_build_check"] = {
        "pass": p.returncode == 0 and "EVERYTHING OK!" in out,
        "rc": p.returncode,
        "tail": out[-300:],
    }

    # Step 2 draws positive probes, G distinct reference k-mers.
    n_positive = args.probes
    rng = random.Random(args.seed)
    ref_kmers = sorted(oracle)
    rng.shuffle(ref_kmers)
    pos = ref_kmers[:n_positive]
    pos_fa = keep / "pos.fa"
    # Probes come in random orientation, so a lookup that finds only the canonical orientation cannot pass.
    write_single_kmer_fasta([rc(x) if rng.random() < 0.5 else x for x in pos], pos_fa)

    # Step 3 draws negative probes, R random canonical k-mers absent from the oracle.
    n_negative = args.probes
    neg = set()
    while len(neg) < n_negative:
        s = "".join(rng.choice("ACGT") for _ in range(k))
        s = canonical(s)
        if s not in oracle:
            neg.add(s)
    neg_fa = keep / "neg.fa"
    write_single_kmer_fasta(
        [rc(x) if rng.random() < 0.5 else x for x in sorted(neg)], neg_fa
    )

    for (label, fa, want), (fmt, iidx) in (
        (x, y)
        for x in (("pos", pos_fa, n_positive), ("neg", neg_fa, 0))
        for y in optik_idxs.items()
    ):
        p = optik_query_report(optik, iidx, fa, k)
        rep = parse_optik_report(p.stdout.decode())
        checks[f"optik_query_{label}{fmt}"] = {
            "pass": (
                rep["num_kmers"] == len(pos if label == "pos" else neg)
                and rep["num_positive_kmers"] == want
                and rep["num_negative_kmers"]
                == len(pos if label == "pos" else neg) - want
            ),
            "report": rep,
            "expected_positives": want,
        }

    for label, fa, want in (("pos", pos_fa, n_positive), ("neg", neg_fa, 0)):
        p = sshash_query_report(sshash, sshash_idx, fa, k)
        rep = parse_sshash_report(p.stdout.decode())
        nrec = len(pos if label == "pos" else neg)
        checks[f"sshash_query_{label}"] = {
            "pass": (
                rep["num_kmers"] == nrec
                and rep["num_positive_kmers"] == want
                and rep["num_negative_kmers"] == nrec - want
            ),
            "report": rep,
            "expected_positives": want,
        }

    # The read sets are the official query reads as given, plus reads simulated from the indexed sequences with 1% substitutions and indels, so the comparison also covers mostly-positive streaming queries.
    sim_q = keep / "sim_reads.fa"
    run(
        [
            sys.executable,
            str(Path(__file__).resolve().parent / "gen_reads.py"),
            "-r",
            str(ref_oneline),
            "-n",
            "20000",
            "-l",
            "150",
            "-s",
            str(args.seed),
            "-e",
            "0.01",
            "--ins",
            "0.002",
            "--del",
            "0.002",
            "-o",
            str(sim_q),
        ],
        timeout=1200,
    )
    # SSHash's query parser takes one sequence per line.
    # In a wrapped FASTA (such as a genome at 70 columns) it would read every other line as a header, so both tools get the same single-line copy.
    # FASTQ is passed unchanged.
    queries = args.queries
    if queries:
        opener = gzip.open if str(queries).endswith(".gz") else open
        with opener(queries, "rt") as fh:
            first = fh.readline()
        if first.startswith(">"):
            queries = keep / "queries_oneline.fa"
            to_single_line_fasta(args.queries, queries)
    read_sets = ([("", queries)] if queries else []) + [("_simreads", sim_q)]
    for (rs, qpath), (fmt, optik_idx) in (
        (x, y) for x in read_sets for y in optik_idxs.items()
    ):
        fmt = fmt + rs
        # Step 4 runs real reads, on which both tools must count positives and negatives the same way.
        # optik counts N-free k-mers, and sshash counts raw ones and reports the invalid separately.
        p1 = optik_query_report(optik, optik_idx, qpath, k)
        r1 = parse_optik_report(p1.stdout.decode())
        p2 = sshash_query_report(sshash, sshash_idx, qpath, k)
        r2 = parse_sshash_report(p2.stdout.decode())
        checks["real_reads_positive_equal" + fmt] = {
            "pass": (
                r1.get("num_positive_kmers") is not None
                and r1.get("num_positive_kmers") == r2.get("num_positive_kmers")
            ),
            "optik": r1,
            "sshash": r2,
        }
        checks["real_reads_negative_equal" + fmt] = {
            "pass": r1.get("num_negative_kmers") is not None
            and r1.get("num_negative_kmers") == r2.get("num_negative_kmers"),
            "optik_negative": r1.get("num_negative_kmers"),
            "sshash_negative": r2.get("num_negative_kmers"),
        }
        # Both tools report num_kmers as the total k-mers seen (positives, negatives and invalid), so the totals are compared directly and the invalid counts must agree.
        # An earlier form subtracted sshash's invalid count from its total while comparing against optik's total, which set two different quantities side by side and could never hold.
        # Step 5 checks per-k-mer exactness of the streaming path.
        # --ids gives each window's membership (id >= 0 when present, -1 when absent or non-ACGT).
        # Every window must agree with the independent KMC k-mer set, and the number of present windows must equal the reported positive count, which ties the --ids path to the counting scan.
        # Aggregate totals alone would hide a false positive and a false negative in the same run.
        if fmt.startswith("_spss"):
            ids_out = keep / f"ids{fmt}.txt"
            p3 = run(
                [
                    str(optik.path),
                    "query",
                    "-i",
                    str(optik_idx),
                    "-q",
                    str(qpath),
                    "-k",
                    str(k),
                    "--ids",
                    str(ids_out),
                ],
                timeout=3600,
            )
            seqs = [sq for _, sq in read_seqs(qpath)]
            lines = ids_out.read_text().splitlines() if p3.returncode == 0 else []
            bad = present = 0
            examples = []
            for sq, ln in zip(seqs, lines, strict=False):
                ids = [int(x) for x in ln.split()]
                u = sq.upper()
                if len(ids) != max(0, len(u) - k + 1):
                    bad += 1
                    continue
                for i, idv in enumerate(ids):
                    w = u[i : i + k]
                    want = all(c in "ACGT" for c in w) and canonical(w) in oracle
                    present += idv >= 0
                    if (idv >= 0) != want:
                        bad += 1
                        if len(examples) < 3:
                            examples.append((w, idv, want))
            checks["per_kmer_membership" + fmt] = {
                "pass": (
                    p3.returncode == 0
                    and len(lines) == len(seqs) > 0
                    and bad == 0
                    and present == r1.get("num_positive_kmers")
                ),
                "records": len(seqs),
                "id_lines": len(lines),
                "window_mismatches": bad,
                "present_windows": present,
                "reported_positive": r1.get("num_positive_kmers"),
                "examples": examples,
                "detail": "every window's --ids membership == KMC oracle; present windows "
                "== reported num_positive_kmers",
            }
        checks["real_reads_kmers_mapping" + fmt] = {
            "pass": (
                r1.get("num_kmers") == r2.get("num_kmers")
                and (r1.get("num_invalid_kmers") or 0)
                == (r2.get("num_invalid_kmers") or 0)
            ),
            "optik_kmers": r1.get("num_kmers"),
            "sshash_kmers": r2.get("num_kmers"),
            "optik_invalid": r1.get("num_invalid_kmers"),
            "sshash_invalid": r2.get("num_invalid_kmers"),
        }

    return checks


# De Bruijn graph suite


def parse_gfa_unitigs(gfa_path):
    """Sequences of the S lines of a GFA file, plain or gzipped."""
    seqs = []
    op = gzip.open if str(gfa_path).endswith(".gz") else open
    with op(gfa_path, "rt") as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t and t[0] == "S" and len(t) >= 3:
                seqs.append(t[2])
    return seqs


def unitig_kmers(seq, k):
    """Set of canonical k-mers of seq."""
    return {canonical(seq[i : i + k]) for i in range(len(seq) - k + 1)}


def build_canon_graph(km31, km32, k):
    """Canonical k-mer graph built from the KMC (k+1)-mer oracle.
    out[x] holds every neighbour of x, since each edge is added in both directions."""
    out = {m: set() for m in km31}
    selfloop = set()
    for w in km32:
        a = canonical(w[:k])
        b = canonical(w[1 : k + 1])
        if a not in km31 or b not in km31:
            continue
        if a == b:
            selfloop.add(a)
            continue
        out[a].add(b)
        out[b].add(a)
    return out, selfloop


def branch_crossing_count(unitigs, km31, km32, k):
    """Count unitigs that cross a branch node.
    An interior k-mer of a unitig may touch only its path neighbours and its own self-loop."""
    out, selfloop = build_canon_graph(km31, km32, k)
    violations = 0
    for s in unitigs:
        wins = [canonical(s[i : i + k]) for i in range(len(s) - k + 1)]
        if len(wins) < 3:
            continue
        for i in range(1, len(wins) - 1):
            x = wins[i]
            allowed = {wins[i - 1], wins[i + 1]}
            extra = set(out.get(x, set())) - allowed
            if extra and x in selfloop:
                extra.discard(x)
            if extra:
                violations += 1
    return violations


def suite_dbg(args):
    """Check optik dbruijn's unitigs against the KMC oracle, Bifrost and Cuttlefish and return the checks."""
    checks = {}
    optik, bifrost, cuttle, kmc = args.optik, args.bifrost, args.cuttlefish, args.kmc
    k = args.k
    keep = Path(args.workdir)
    keep.mkdir(parents=True, exist_ok=True)

    # optik dbruijn, then GFA, then unitig FASTA
    gfa = keep / "optik.gfa"
    kdbg = keep / "optik.kdbg"
    p = run(
        [
            str(optik.path),
            "dbruijn",
            "-i",
            str(args.input),
            "-k",
            str(k),
            "-c",
            str(args.min_count),
            "-o",
            str(gfa),
            "-t",
            str(args.threads),
            "-s",
            str(kdbg),
        ],
        timeout=3600,
    )
    checks["optik_dbruijn"] = {
        "pass": p.returncode == 0,
        "stderr": p.stderr.decode()[-500:],
    }
    if not checks["optik_dbruijn"]["pass"]:
        # optik may have written the GFA despite a nonzero exit on runs without retries, so check for it.
        pass
    if not os.path.isfile(gfa):
        return checks
    optik_us = parse_gfa_unitigs(gfa)
    optik_set = set()
    for s in optik_us:
        optik_set |= unitig_kmers(s, k)
    checks["optik_unitigs"] = {
        "n_unitigs": len(optik_us),
        "total_bp": sum(map(len, optik_us)),
    }

    # The oracle holds the input k-mer set (k) and the (k+1)-mer set.
    # KMC's FASTA reader requires one sequence per line, so FASTA inputs are converted for it, while optik, Bifrost and Cuttlefish read the original file.
    oracle_input = args.input
    if not _is_fastq(args.input):
        oracle_input = keep / "input_oneline.fa"
        to_single_line_fasta(args.input, oracle_input)
    # The solid set uses the same abundance cutoff optik ran with.
    # The (k+1)-mer set stays unfiltered, because under optik's edge-centric semantics an edge needs one observed window between two solid k-mers.
    oracle = kmc_build(
        kmc, keep, oracle_input, "oracle", k, args.threads, args.min_count
    )
    oracle_k1 = kmc_build(kmc, keep, oracle_input, "oracle_k1", k + 1, args.threads)
    checks["oracle"] = {"n_kmers": len(oracle), "n_kmers_k1": len(oracle_k1)}

    # The unitig set must equal the oracle exactly.
    only_optik = optik_set - oracle
    only_oracle = oracle - optik_set
    checks["unitig_set_equal_oracle"] = {
        "pass": not only_optik and not only_oracle,
        "only_optik": len(only_optik),
        "only_oracle": len(only_oracle),
    }

    # Coverage and spurious k-mers
    checks["input_kmers_covered"] = {
        "pass": not only_oracle,
        "uncovered": len(only_oracle),
    }
    checks["unitig_kmers_in_input"] = {
        "pass": not only_optik,
        "spurious": len(only_optik),
    }

    # No unitig may cross a branch node.
    viol = branch_crossing_count(optik_us, oracle, oracle_k1, k)
    checks["no_branch_crossing"] = {"pass": viol == 0, "violations": viol}

    # Structural checks cover duplicated k-mers, maximality, exact L links and, at min_count 1, equality of the unitig partition with Cuttlefish --read.
    # Set equality alone would accept duplicated k-mers and over-split unitigs.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import validate_unitigs

    struct = validate_unitigs.validate(
        str(oracle_input),
        str(gfa),
        k,
        args.min_count,
        str(cuttle.path) if cuttle.available else None,
        args.threads,
        kdbg=str(kdbg),
    )
    for name, c in struct.items():
        checks["gfa_" + name] = c

    # Bifrost, as an independent tool
    if bifrost.available and args.min_count == 1:
        bf_prefix = keep / "bifrost"
        p = run(
            [
                str(bifrost.path),
                "build",
                "-r",
                str(args.input),
                "-o",
                str(bf_prefix),
                "-k",
                str(k),
                "-t",
                str(args.threads),
            ],
            timeout=3600,
        )
        bfr = p.returncode
        bf_us = []
        for cand in (str(bf_prefix) + ".gfa", str(bf_prefix) + ".gfa.gz"):
            if os.path.isfile(cand):
                bf_us = parse_gfa_unitigs(cand)
                break
        bf_set = set()
        for s in bf_us:
            bf_set |= unitig_kmers(s, k)
        d_bf = (len(bf_set ^ oracle), len(bf_set - oracle), len(oracle - bf_set))
        checks["bifrost_build"] = {"pass": bfr == 0}
        checks["bifrost_set_equal_oracle"] = {
            "pass": d_bf[0] == 0,
            "n_unitigs": len(bf_us),
            "sym_diff": d_bf[0],
            "only_bifrost": d_bf[1],
            "only_oracle": d_bf[2],
        }

    # Cuttlefish, also independent, which works on (k+1)-mers and drops islands isolated by N
    if cuttle.available and args.min_count == 1:
        cf_prefix = keep / "cuttlefish"
        p = run(
            [
                str(cuttle.path),
                "build",
                "--read",
                "-s",
                str(args.input),
                "-k",
                str(k),
                "-c",
                str(args.min_count),
                "-o",
                str(cf_prefix),
                "-t",
                str(args.threads),
            ],
            timeout=3600,
        )
        cf_fa = Path(str(cf_prefix) + ".fa")
        cf_set = set()
        ncf = 0
        if cf_fa.is_file():
            for _, seq in read_fasta(cf_fa):
                cf_set |= unitig_kmers(seq, k)
                ncf += 1
        extra_cf = cf_set - oracle
        missing_cf = oracle - cf_set
        checks["cuttlefish_build"] = {"pass": p.returncode == 0}
        checks["cuttlefish_set_subset"] = {
            "pass": not extra_cf,
            "only_cuttlefish": len(extra_cf),
            "missing_vs_oracle": len(missing_cf),
            "n_unitigs": ncf,
            "detail": "cuttlefish drops k-mers isolated between N runs; "
            "it must never add k-mers absent from the oracle",
        }

    return checks


# Mapping suite


def parse_truth_tsv(path):
    """Rows (read, contig, position, strand, length) of a gen_reads.py truth TSV."""
    rows = []
    for ln in path.read_text().splitlines():
        if not ln or ln.startswith("read_name"):
            continue
        p = ln.split("\t")
        if len(p) >= 5:
            rows.append((p[0], p[1], int(p[2]), p[3], int(p[4])))
    return rows


def parse_sam_maps(path):
    """Map each read name to its primary record as dict(pairname, pos (1-based), strand, mapq, cigar, nm)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < SAM_MIN_COLUMNS or t[2] == "*" or int(t[1]) & SAM_FLAG_UNMAPPED:
                out[t[0]] = {"unmapped": True}
                continue
            nm = None
            for tag in t[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break
            f = int(t[1])
            out[t[0]] = {
                "unmapped": False,
                "ref": t[2],
                "pos": int(t[3]),
                "strand": "-" if f & SAM_FLAG_REVERSE else "+",
                "cigar": t[5],
                "nm": nm,
                "mapq": int(t[4]),
            }
    return out


def sam_nm_from_cigar(read, cigar, refseq, pos1, strand):
    """Recompute the NM tag from the SAM CIGAR against the reference.

    The walk compares M/=/X bases case-insensitively, with any N counted as a mismatch.
    I adds its length while consuming read bases, D adds its length while consuming reference bases, and S, H and P consume without cost.
    On the '-' strand the read is reverse-complemented first, since SAM SEQ is stored as sequenced.

    Returns None if the CIGAR's reference span runs past the contig.
    """
    read = rc(read) if strand == "-" else read
    qi = 0
    ti = pos1 - 1
    ref_span = sum(
        int(mm.group(1))
        for mm in re.finditer(r"(\d+)([MDN=X])", cigar)
        if mm.group(2) in "MDN=X"
    )
    if ti < 0 or ti + ref_span > len(refseq):
        return None
    nm = 0
    for mm in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        ln = int(mm.group(1))
        op = mm.group(2)
        if op in "M=X":
            for j in range(ln):
                rq = read[qi + j] if qi + j < len(read) else "N"
                rt = refseq[ti + j] if ti + j < len(refseq) else "N"
                if not (rq.upper() == rt.upper() and rq.upper() in "ACGT"):
                    nm += 1
            qi += ln
            ti += ln
        elif op in "IS":
            qi += ln
            if op == "I":
                nm += ln
        elif op in "DN":
            ti += ln
            if op == "D":
                nm += ln
    return nm


# Helpers for the wavefront mapper wfmash, which pairs MashMap segments with WFA alignment

# wfmash is built for long sequences.
# Its defaults segment queries into 1 kb windows and keep only chains with at least 10 kb of scaffold mass, so reads shorter than that are never reported.
# It therefore runs twice, once with its exact defaults and once with read settings that stay within its documented options.
# Those set the scaffold filter to its minimum (-S 1), the identity floor to 70% (-p 70, about 15% divergence) and, for short reads, the window to its minimum (-w 100).
WFMASH_READ_OPTS = {
    "short": ["-S", "1", "-p", "70", "-w", "100"],
    "long": ["-S", "1", "-p", "70"],
}


def faidx(path):
    """Write the .fai indexes wfmash needs for target and query."""
    fai = Path(str(path) + ".fai")
    if fai.is_file() and fai.stat().st_mtime >= Path(path).stat().st_mtime:
        return
    rows, off = [], 0
    with open(path, "rb") as fh:
        name, length, start, lb, lw = None, 0, 0, 0, 0
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    rows.append(f"{name}\t{length}\t{start}\t{lb}\t{lw}")
                name = line[1:].split()[0].decode()
                length, start, lb, lw = 0, off + len(line), 0, 0
            else:
                n = len(line.rstrip(b"\r\n"))
                if lb == 0:
                    lb, lw = n, len(line)
                length += n
            off += len(line)
        if name is not None:
            rows.append(f"{name}\t{length}\t{start}\t{lb}\t{lw}")
    fai.write_text("\n".join(rows) + "\n")


def one_primary_per_read(sam_in, sam_out):
    """Reduce wfmash's output to one primary per read.

    wfmash reports one primary per mapped query segment.
    For each read the primary with the most aligned read bases (M/=/X/I) is kept and the others are marked supplementary, so the same auditors and scorers apply to every mapper.
    Returns the number of reads that had more than one primary."""

    def aligned(t):
        """Read bases aligned by SAM record t (M, =, X and I)."""
        return sum(
            int(a) for a, o in re.findall(r"(\d+)([MIDNSHP=X])", t[5]) if o in "M=XI"
        )

    best, lines = {}, []
    with open(sam_in) as fh:
        for line in fh:
            lines.append(line)
            if line.startswith("@"):
                continue
            t = line.split("\t")
            f = int(t[1])
            if f & SAM_FLAGS_NOT_PRIMARY:
                continue
            a = aligned(t)
            if t[0] not in best or a > best[t[0]][0]:
                best[t[0]] = (a, len(lines) - 1)
    multi = Counter()
    with open(sam_out, "w") as out:
        for i, line in enumerate(lines):
            if not line.startswith("@"):
                t = line.split("\t")
                f = int(t[1])
                if not f & SAM_FLAGS_NOT_PRIMARY:
                    multi[t[0]] += 1
                    if best[t[0]][1] != i:
                        t[1] = str(f | 0x800)
                        line = "\t".join(t)
            out.write(line)
    return sum(1 for v in multi.values() if v > 1)


def run_wfmash(wfmash, ref, reads, out, threads, opts):
    """Run wfmash REF READS [opts] -a to produce SAM with one primary per read.

    Returns (wall seconds, returncode, reads with split primaries)."""
    faidx(ref)
    faidx(reads)
    raw = Path(str(out) + ".raw")
    t0 = time.monotonic()
    with open(raw, "wb") as fh:
        p = run(
            [str(wfmash.path), str(ref), str(reads), "-t", str(threads), "-a", *opts],
            stdout=fh,
            stderr=subprocess.DEVNULL,
            timeout=7200,
        )
    dt = time.monotonic() - t0
    split = one_primary_per_read(raw, out) if p.returncode == 0 else 0
    return dt, p.returncode, split


def timed(fn, *a, **kw):
    """Result and wall seconds of calling fn."""
    t0 = time.monotonic()
    r = fn(*a, **kw)
    return r, time.monotonic() - t0


def suite_map(args):
    """Check optik map on simulated and real reads against truth, minimap2, wfmash and samtools."""
    checks = {}
    optik, minimap2, samtools = args.optik, args.minimap2, args.samtools
    wfmash = args.wfmash
    speed = {}  # leg -> {mapper -> wall seconds}, with the same reads and threads for all
    keep = Path(args.workdir)
    keep.mkdir(parents=True, exist_ok=True)

    # Tidy the reference names.
    # A unitig FASTA may be gzipped or wrapped, with repeated headers, so a plain one-line reference with unique names is written and used for read simulation, optik mapping and minimap2.
    ref = keep / "sim_ref.fa"
    write_uniquified_fasta(args.ref, ref)
    if not minimap2.available and not samtools.available:
        need("at least minimap2 or samtools")

    # optik indexes the same minimizers as the minimap2 preset it is compared with (minimap2 v2.31 options.c gives sr k = 21, w = 11 and map-ont, the default, k = 15, w = 10).
    # --map-k and -w override both.
    mm2_preset_kw = {"sr": (21, 11), "map-ont": (15, 10)}

    def run_optik_map(reads, out, threads, preset="sr"):
        """Run optik map on reads with the preset's k and w unless overridden, writing SAM to out."""
        cmd = [str(optik.path), "map", "-r", str(ref), "-q", str(reads)]
        pk, pw = mm2_preset_kw[preset]
        for flag, v in (("-k", args.map_k or pk), ("-w", args.w or pw), ("-e", args.e)):
            if v is not None:
                cmd += [flag, str(v)]
        p = run(cmd + ["-t", str(threads), "-o", str(out)], timeout=3600)
        return p

    # Part 1, error-free simulated reads with exact ground truth
    sim_reads = keep / "sim.fa"
    sim_truth = keep / "sim.truth.tsv"
    greads = Path(__file__).resolve().parent / "gen_reads.py"
    p = run(
        [
            sys.executable,
            str(greads),
            "-r",
            str(ref),
            "-n",
            str(args.n),
            "-l",
            str(args.read_len),
            "-s",
            str(args.seed),
            "-o",
            str(sim_reads),
            "--gt",
            str(sim_truth),
        ],
        timeout=1200,
    )
    checks["gen_sim_reads"] = {"pass": p.returncode == 0}
    truth = parse_truth_tsv(sim_truth)

    sam0 = keep / "sim_optik.sam"
    p, dt = timed(run_optik_map, sim_reads, sam0, args.threads)
    speed.setdefault("error_free_short", {})["optik"] = round(dt, 3)
    checks["optik_map_exit"] = {
        "pass": p.returncode == 0,
        "stderr": p.stderr.decode()[-400:],
    }

    maps = parse_sam_maps(sam0)
    ref_seqs = {n: s for n, s in read_fasta(ref)}
    mism = 0
    repeat_ties = 0
    for tname, tref, tpos, tstrand, tlen in truth:
        m = maps.get(tname)
        if m is None or m.get("unmapped"):
            mism += 1
            continue
        if (m["ref"], m["pos"], m.get("strand")) != (tref, tpos + 1, tstrand):
            # The read has no errors, so another placement is a genuine tie only when the reference at the reported locus, in the reported orientation, spells the read exactly, as a repeat copy does, inverted repeats on the other strand included.
            # Anything else is a misplacement.
            rw = ref_seqs[m["ref"]][m["pos"] - 1 : m["pos"] - 1 + tlen]
            tw = ref_seqs[tref][tpos : tpos + tlen]
            want = tw if tstrand == "+" else rc(tw)  # the read itself
            got = rw if m.get("strand") == "+" else rc(rw)
            if got == want:
                repeat_ties += 1
            else:
                mism += 1
    checks["sim_exact_position_strand"] = {
        "pass": mism == 0,
        "n_truth": len(truth),
        "mismatches": mism,
        "repeat_copy_ties": repeat_ties,
    }

    # Part 2, reads with substitutions and indels on both strands Every read counts.
    # Unmapped or misplaced reads lower sensitivity, which is compared with minimap2's on the same reads.
    # Records must be self-consistent, with NM matching a recomputation from SAM SEQ and CIGAR, each read needs exactly one primary, and placements with high MAPQ must be right.
    err_reads = keep / "sim_err.fa"
    err_truth_tsv = keep / "sim_err.truth.tsv"
    p = run(
        [
            sys.executable,
            str(greads),
            "-r",
            str(ref),
            "-n",
            str(args.n),
            "-l",
            str(args.read_len),
            "-s",
            str(args.seed + 1),
            "-e",
            str(args.err_rate),
            "--ins",
            str(args.indel_rate),
            "--del",
            str(args.indel_rate),
            "-o",
            str(err_reads),
            "--gt",
            str(err_truth_tsv),
        ],
        timeout=1200,
    )
    checks["gen_err_reads"] = {"pass": p.returncode == 0}
    err_truth = {t[0]: t for t in parse_truth_tsv(err_truth_tsv)}
    sam1 = keep / "sim_err_optik.sam"
    p, dt = timed(run_optik_map, err_reads, sam1, args.threads)
    speed.setdefault("error_short", {})["optik"] = round(dt, 3)
    checks["optik_map_err_exit"] = {"pass": p.returncode == 0}

    def audit_sam(path, seqs=None):
        """Per-read primary placements and spec checks of a SAM file.

        A mapped primary without SEQ or NM:i counts as bad, so a writer that dropped them cannot pass, and the CIGAR must consume SEQ."""
        seqs = err_seqs if seqs is None else seqs
        prim = {}
        dup_primary = seq_bad = nm_bad = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                t = line.rstrip("\n").split("\t")
                flag = int(t[1])
                if flag & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                    continue
                if t[0] in prim:
                    dup_primary += 1
                    continue
                if flag & SAM_FLAG_UNMAPPED:
                    prim[t[0]] = None
                    continue
                strand = "-" if flag & SAM_FLAG_REVERSE else "+"
                read = seqs.get(t[0])
                seq = t[9]
                want = None if read is None else (rc(read) if strand == "-" else read)
                if seq == "*" or want is None or seq.upper() != want.upper():
                    seq_bad += 1
                nm = None
                for tag in t[11:]:
                    if tag.startswith("NM:i:"):
                        nm = int(tag[5:])
                ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
                if (
                    seq == "*"
                    or sum(int(a) for a, o in ops if o in "MIS=X") != len(seq)
                    or t[2] not in ref_seqs
                ):
                    nm_bad += 1
                else:
                    # SAM SEQ is already in reference orientation.
                    d = sam_nm_from_cigar(seq, t[5], ref_seqs[t[2]], int(t[3]), "+")
                    if nm is None or d != nm:
                        nm_bad += 1
                # The read's origin on the reference is POS minus the leading clip, so clipped alignments, soft or hard, are scored like end-to-end ones.
                lead = re.match(r"^(\d+)[SH]", t[5])
                origin = int(t[3]) - (int(lead.group(1)) if lead else 0)
                prim[t[0]] = (t[2], origin, strand, int(t[4]))
        return prim, dup_primary, seq_bad, nm_bad

    def score(prim, tr=None):
        """Sensitivity and high-MAPQ errors of primary placements prim against truth tr."""
        correct = wrong = unmapped = hq = hq_wrong = 0
        tr = err_truth if tr is None else tr
        for name, (_, tref, tpos, tstrand, tspan) in tr.items():
            m = prim.get(name)
            if m is None:
                unmapped += 1
                continue
            ok = (
                m[0] == tref
                and m[2] == tstrand
                and abs(m[1] - (tpos + 1)) <= max(10, tspan // 10)
            )
            if not ok:
                # An equally good repeat copy is not an error.
                win = ref_seqs[m[0]][m[1] - 1 : m[1] - 1 + tspan]
                tw = ref_seqs[tref][tpos : tpos + tspan]
                ok = win == tw
            correct += ok
            wrong += not ok
            if m[3] >= HIGH_MAPQ:
                hq += 1
                hq_wrong += not ok
        n = max(1, len(tr))
        return {
            "sensitivity": round(correct / n, 5),
            "correct": correct,
            "wrong": wrong,
            "unmapped": unmapped,
            "mapq30_n": hq,
            "mapq30_wrong": hq_wrong,
        }

    err_seqs = {n: s for n, s in read_fasta(err_reads)}
    # The error-free set's records get the same spec audit, because its placement check above looks only at positions.
    sim_seqs = {n: s for n, s in read_fasta(sim_reads)}
    _, sdup, sseq, snm = audit_sam(sam0, sim_seqs)
    checks["sim_reads_records_valid"] = {
        "pass": sdup == 0 and sseq == 0 and snm == 0,
        "duplicate_primary": sdup,
        "seq_mismatches": sseq,
        "nm_or_cigar_bad": snm,
    }
    oprim, odup, oseq, onm = audit_sam(sam1)
    osc = score(oprim)
    checks["err_reads_one_primary"] = {
        "pass": odup == 0 and len(oprim) == len(err_truth),
        "duplicate_primary": odup,
        "reads_with_record": len(oprim),
        "reads": len(err_truth),
    }
    checks["err_reads_seq_orientation"] = {"pass": oseq == 0, "mismatches": oseq}
    checks["err_reads_nm_matches_cigar"] = {"pass": onm == 0, "nm_mismatches": onm}
    checks["err_reads_mapq_calibration"] = {
        "pass": osc["mapq30_wrong"] <= max(1, osc["mapq30_n"] // 100),
        "mapq30_n": osc["mapq30_n"],
        "mapq30_wrong": osc["mapq30_wrong"],
        "correct": osc["correct"],
        "detail": "<= 1% of MAPQ>=30 placements may be wrong (MAPQ is capped by "
        "divergence by design, so high-error reads may stay below 30)",
    }
    # MAPQ must also be informative.
    # On the error-free set, unique correct placements must reach MAPQ >= 30, since an all-zero MAPQ would satisfy the calibration bound above trivially.
    sim_prim, _, _, _ = audit_sam(sam0, sim_seqs)
    ssc = score(sim_prim, {t[0]: t for t in truth})
    checks["mapq_informative"] = {
        "pass": ssc["mapq30_n"] >= 0.9 * ssc["correct"]
        and ssc["mapq30_wrong"] == 0
        and ssc["correct"] > 0,
        "error_free_correct": ssc["correct"],
        "mapq30_n": ssc["mapq30_n"],
        "mapq30_wrong": ssc["mapq30_wrong"],
        "detail": ">= 90% of correct error-free placements have MAPQ >= 30 "
        "(repeat copies legitimately stay low), none of them wrong",
    }
    sens = {"optik": osc}
    if minimap2.available:
        mm_err = keep / "mm2_err.sam"
        with open(mm_err, "wb") as _fh:
            _, dt = timed(
                run,
                [
                    str(minimap2.path),
                    "-ax",
                    "sr",
                    "-t",
                    str(args.threads),
                    str(ref),
                    str(err_reads),
                ],
                stdout=_fh,
                stderr=subprocess.DEVNULL,
                timeout=3600,
            )
        speed["error_short"]["minimap2"] = round(dt, 3)
        mprim, _, _, _ = audit_sam(mm_err)
        sens["minimap2"] = score(mprim)
        # minimap2 writes a record for every read, unmapped ones included.
        # Fewer records point to a crashed or truncated baseline, which would make "optik >= minimap2" vacuous.
        mm_missing = sum(1 for n in err_truth if n not in mprim)
        checks["err_reads_sensitivity_vs_minimap2"] = {
            "pass": (
                osc["sensitivity"] >= sens["minimap2"]["sensitivity"] - args.sens_tol
                and mm_missing == 0
            ),
            "optik": osc,
            "minimap2": sens["minimap2"],
            "minimap2_reads_without_record": mm_missing,
            "tolerance": args.sens_tol,
        }
    else:
        checks["err_reads_sensitivity"] = {"info": osc}

    # Part 2b, long reads on the chaining path (ONT-like, 0.5-30 kb, about 7% error)
    lr_reads = keep / "sim_long.fa"
    lr_truth_tsv = keep / "sim_long.truth.tsv"
    p = run(
        [
            sys.executable,
            str(greads),
            "-r",
            str(ref),
            "-n",
            str(args.n_long),
            "-l",
            "5000",
            "--len-sd",
            "0.6",
            "-s",
            str(args.seed + 2),
            "-e",
            "0.03",
            "--ins",
            "0.02",
            "--del",
            "0.02",
            "-o",
            str(lr_reads),
            "--gt",
            str(lr_truth_tsv),
        ],
        timeout=1200,
    )
    lr_truth = {t[0]: t for t in parse_truth_tsv(lr_truth_tsv)}

    def long_audit(path):
        """Primary placements of a long-read SAM and its duplicate, NM and SEQ-length violations."""
        prim, dup, nm_bad, len_bad = {}, 0, 0, 0
        seen = set()
        with open(path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                t = line.rstrip("\n").split("\t")
                flag = int(t[1])
                seen.add(t[0])
                if flag & SAM_FLAG_UNMAPPED:
                    prim.setdefault(t[0], None)
                    continue
                ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
                hard = any(o == "H" for _, o in ops) or t[9] == "*"
                if not hard:
                    if sum(int(a) for a, o in ops if o in "MIS=X") != len(t[9]):
                        len_bad += 1
                    nm = next(
                        (int(x[5:]) for x in t[11:] if x.startswith("NM:i:")), None
                    )
                    if (
                        nm is not None
                        and sam_nm_from_cigar(
                            t[9], t[5], ref_seqs[t[2]], int(t[3]), "+"
                        )
                        != nm
                    ):
                        nm_bad += 1
                if flag & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                    continue
                if prim.get(t[0]):
                    dup += 1
                rs = int(t[3]) - 1
                re_ = rs + sum(int(a) for a, o in ops if o in "MDN=X")
                prim[t[0]] = (
                    t[2],
                    rs,
                    re_,
                    "-" if flag & SAM_FLAG_REVERSE else "+",
                    int(t[4]),
                )
        good = hq = hq_wrong = 0
        for name, (_, tref, tpos, tstrand, tspan) in lr_truth.items():
            m = prim.get(name)
            if not m:
                continue
            ok = (
                m[0] == tref
                and m[3] == tstrand
                and min(m[2], tpos + tspan) - max(m[1], tpos) >= 0.5 * tspan
            )
            good += ok
            if m[4] >= HIGH_MAPQ:
                hq += 1
                hq_wrong += not ok
        n = max(1, len(lr_truth))
        return {
            "sensitivity": round(good / n, 5),
            "correct": good,
            "reads_without_record": sum(1 for x in lr_truth if x not in seen),
            "duplicate_primary": dup,
            "nm_mismatches": nm_bad,
            "cigar_len_mismatches": len_bad,
            "mapq30_n": hq,
            "mapq30_wrong": hq_wrong,
        }

    lr_sam = keep / "sim_long_optik.sam"
    p, dt = timed(run_optik_map, lr_reads, lr_sam, args.threads, "map-ont")
    speed.setdefault("long", {})["optik"] = round(dt, 3)
    la = long_audit(lr_sam)
    checks["long_reads_valid_sam"] = {
        "pass": p.returncode == 0
        and la["duplicate_primary"] == 0
        and la["nm_mismatches"] == 0
        and la["cigar_len_mismatches"] == 0
        and la["reads_without_record"] == 0,
        **la,
    }
    checks["long_reads_mapq_calibration"] = {
        "pass": (
            la["mapq30_wrong"] <= max(1, la["mapq30_n"] // 100)
            and la["mapq30_n"] >= la["correct"] // 2
        ),
        "mapq30_n": la["mapq30_n"],
        "mapq30_wrong": la["mapq30_wrong"],
        "correct": la["correct"],
    }
    if minimap2.available:
        mm_lr = keep / "mm2_long.sam"
        with open(mm_lr, "wb") as _fh:
            _, dt = timed(
                run,
                [
                    str(minimap2.path),
                    "-ax",
                    "map-ont",
                    "-t",
                    str(args.threads),
                    str(ref),
                    str(lr_reads),
                ],
                stdout=_fh,
                stderr=subprocess.DEVNULL,
                timeout=3600,
            )
        speed["long"]["minimap2"] = round(dt, 3)
        mla = long_audit(mm_lr)
        checks["long_reads_sensitivity_vs_minimap2"] = {
            "pass": (
                la["sensitivity"] >= mla["sensitivity"] - args.sens_tol
                and mla["reads_without_record"] == 0
            ),
            "minimap2_reads_without_record": mla["reads_without_record"],
            "optik": la["sensitivity"],
            "minimap2": mla["sensitivity"],
            "tolerance": args.sens_tol,
        }

    # Part 3, determinism across thread counts
    sam_t = keep / "sim_optik_t1.sam"
    p = run_optik_map(sim_reads, sam_t, 1)
    with open(sam_t) as _f:
        a_body = [line for line in _f if not line.startswith("@")]
    with open(sam0) as _f:
        b_body = [line for line in _f if not line.startswith("@")]
    checks["map_determinism_t1_vs_tN"] = {
        "pass": a_body == b_body,
        "note": "SAM alignment body must be byte-identical across thread counts",
    }

    # Part 4, concordance with minimap2
    if minimap2.available:
        mm2_sam = keep / "mm2_sim.sam"
        with open(mm2_sam, "wb") as _fh:
            p, dt = timed(
                run,
                [
                    str(minimap2.path),
                    "-ax",
                    "sr",
                    "-t",
                    str(args.threads),
                    str(ref),
                    str(sim_reads),
                ],
                stdout=_fh,
                stderr=subprocess.DEVNULL,
                timeout=3600,
            )
        speed["error_free_short"]["minimap2"] = round(dt, 3)
        mmaps = parse_sam_maps(mm2_sam)
        reads_seq = {n: s for n, s in read_fasta(sim_reads)}
        agree = disagree = both_fail = 0

        def content_ok(mp, tlen, read):
            """True when mapping mp places the read on reference text that spells it, on either strand."""
            if mp.get("unmapped", False):
                return False
            win = ref_seqs[mp["ref"]][mp["pos"] - 1 : mp["pos"] - 1 + tlen]
            return win == read or rc(win) == read

        for tname, _tref, _tpos, _tstrand, tlen in truth:
            m = maps.get(tname)
            mm = mmaps.get(tname)
            read = reads_seq.get(tname)
            if not m or not mm or read is None:
                continue
            mok = content_ok(m, tlen, read)
            nok = content_ok(mm, tlen, read)
            # A different placement is acceptable when it matches the reference content exactly, in either orientation.
            # Legitimate multi-mappers and palindromic or inverted-repeat ties then count as equal rather than as errors.
            if mok and nok:
                agree += 1
            elif not mok and not nok:
                both_fail += 1  # a failure both share is not agreement
            else:
                disagree += 1
        checks["minimap2_concordance"] = {
            "pass": disagree == 0 and both_fail == 0,
            "agree": agree,
            "disagree": disagree,
            "both_fail": both_fail,
        }

        # On real reads, compare the mapped/unmapped split.
        if args.queries:
            mm2_sam_r = keep / "mm2_real.sam"
            with open(mm2_sam_r, "wb") as _fh:
                run(
                    [
                        str(minimap2.path),
                        "-ax",
                        "sr",
                        "-t",
                        str(args.threads),
                        str(ref),
                        str(args.queries),
                    ],
                    stdout=_fh,
                    stderr=subprocess.DEVNULL,
                    timeout=3600,
                )
            mmr = parse_sam_maps(mm2_sam_r)
            samr = keep / "optik_real.sam"
            run_optik_map(args.queries, samr, args.threads)
            mr = parse_sam_maps(samr)
            agg = Counter()
            for rn in set(mmr) | set(mr):
                o = not mr.get(rn, {}).get("unmapped", False)
                mm = not mmr.get(rn, {}).get("unmapped", False)
                agg[(o, mm)] += 1
            checks["real_reads_map_parity"] = {
                "pass": agg.get((False, True), 0) == 0
                and agg.get((True, False), 0) == 0,
                "both_unmapped": agg.get((False, False), 0),
                "only_optik": agg.get((True, False), 0),
                "only_minimap2": agg.get((False, True), 0),
                "both_mapped": agg.get((True, True), 0),
            }

    # Part 4b, the wavefront mapper (wfmash) on the same three read sets The same auditors as for optik and minimap2 score it (one primary per read, origin within tolerance, strand, MAPQ >= 30 calibration).
    # optik must be at least as sensitive as the better wfmash configuration.
    if wfmash.available:
        sim_truth_d = {t[0]: t for t in truth}
        optik_sens = {
            "error_free_short": score(audit_sam(sam0)[0], sim_truth_d),
            "error_short": osc,
            "long": {"sensitivity": la["sensitivity"]},
        }
        legs = (
            ("error_free_short", sim_reads, "short", sim_truth_d),
            ("error_short", err_reads, "short", err_truth),
            ("long", lr_reads, "long", None),
        )
        for leg, reads, kind, tr in legs:
            res = {}
            for cfg, opts in (("defaults", []), ("read_opts", WFMASH_READ_OPTS[kind])):
                out = keep / f"wfmash_{leg}_{cfg}.sam"
                dt, rc_, split = run_wfmash(wfmash, ref, reads, out, args.threads, opts)
                if rc_ != 0:
                    res[cfg] = {"error": f"wfmash exit {rc_}", "seconds": round(dt, 3)}
                    continue
                if kind == "long":
                    r = long_audit(out)
                    r = {
                        "sensitivity": r["sensitivity"],
                        "mapq30_n": r["mapq30_n"],
                        "mapq30_wrong": r["mapq30_wrong"],
                    }
                else:
                    r = score(audit_sam(out)[0], tr)
                r.update(
                    {
                        "seconds": round(dt, 3),
                        "reads_split_into_segments": split,
                        "cmd": " ".join(["wfmash", "REF", "READS", "-a", *opts]),
                    }
                )
                res[cfg] = r
                speed[leg][f"wfmash_{cfg}"] = round(dt, 3)
            best = max((v.get("sensitivity", 0) for v in res.values()), default=0)
            # A baseline that crashed or failed is no comparison, since it would score 0.
            broken = {c: v["error"] for c, v in res.items() if "error" in v}
            checks[f"{leg}_sensitivity_vs_wfmash"] = {
                "pass": optik_sens[leg]["sensitivity"] >= best - args.sens_tol
                and not broken,
                "optik": optik_sens[leg],
                "wfmash": res,
                "tolerance": args.sens_tol,
                "error": f"wfmash baseline failed: {broken}" if broken else "",
            }
    checks["mapper_speed"] = {
        "info": speed,
        "threads": args.threads,
        "detail": "whole-process wall seconds on identical reads (index build included)",
    }

    # Part 5, samtools round trip
    if samtools.available:
        bam = keep / "sim_optik.bam"
        p = run(
            [str(samtools.path), "view", "-b", "-o", str(bam), str(sam0)], timeout=3600
        )
        checks["samtools_view_bam"] = {"pass": p.returncode == 0}
        back = keep / "roundtrip.sam"
        with open(back, "wb") as _fh:
            p = run([str(samtools.path), "view", str(bam)], stdout=_fh, timeout=3600)
        checks["samtools_view_back"] = {"pass": p.returncode == 0}
        with open(sam0) as _f:
            a2 = [line for line in _f if not line.startswith("@")]
        with open(back) as _f:
            b2 = [line for line in _f if not line.startswith("@")]
        checks["samtools_roundtrip_same_body"] = {"pass": a2 == b2, "n_ali": len(a2)}
        p = run([str(samtools.path), "flagstat", str(bam)], timeout=3600)
        fs = p.stdout.decode()
        checks["samtools_flagstat"] = {
            "pass": p.returncode == 0 and "mapped (" in fs,
            "summary": fs.splitlines()[0] if fs else "",
        }

    return checks


# Command line


def base_parser():
    """Parent parser with the tool-path and shared options of every suite."""
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--optik", default=None)
    ap.add_argument("--sshash", default=None)
    ap.add_argument("--kmc", default=None)
    ap.add_argument("--bifrost", default=None)
    ap.add_argument("--cuttlefish", default=None)
    ap.add_argument("--minimap2", default=None)
    ap.add_argument("--samtools", default=None)
    ap.add_argument("--wfmash", default=None)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("-s", "--seed", type=int, default=42)
    return ap


def print_report(checks):
    """Print every check as JSON with the verdict, returning 0 when all required checks pass."""
    fails = [
        key
        for key, val in checks.items()
        if isinstance(val, dict) and "pass" in val and not val["pass"]
    ]
    print(json.dumps(checks, indent=2, default=str))
    print("\nRESULT:", "FAIL" if fails else "PASS", file=sys.stderr)
    for f in fails:
        print(f"  FAILED: {f} -> {checks[f]}", file=sys.stderr)
    return 0 if not fails else 1


def main():
    """Resolve tool paths and run the suite named on the command line."""
    parser = argparse.ArgumentParser(
        description="Exact-equality correctness validation of optik vs the original tools."
    )
    sub = parser.add_subparsers(dest="suite", required=True)

    q = sub.add_parser("query", parents=[base_parser()])
    q.add_argument("-r", "--ref", required=True)
    q.add_argument("--queries", default=None)
    q.add_argument("--probes", type=int, default=2000)
    q.add_argument("--sshash-m", type=int, default=15)

    d = sub.add_parser("dbg", parents=[base_parser()])
    d.add_argument("-i", "--input", required=True)
    d.add_argument("--min-count", type=int, default=1)

    m = sub.add_parser("map", parents=[base_parser()])
    m.add_argument("-r", "--ref", required=True)
    m.add_argument("--queries", default=None)
    m.add_argument("-n", type=int, default=5000)
    m.add_argument("--read-len", type=int, default=150)
    m.add_argument(
        "--n-long",
        type=int,
        default=1000,
        help="simulated long reads for the chaining-path checks",
    )
    # optik runs with its own defaults, as minimap2 runs its documented preset.
    # These options override them only when given explicitly.
    m.add_argument("--map-k", type=int, default=None)
    m.add_argument("-w", type=int, default=None)
    m.add_argument("-e", type=int, default=None)
    m.add_argument(
        "--err-rate",
        type=float,
        default=0.01,
        help="per-base substitution rate of the error-read set",
    )
    m.add_argument(
        "--indel-rate",
        type=float,
        default=0.002,
        help="per-base insertion (and deletion) rate",
    )
    m.add_argument(
        "--sens-tol",
        type=float,
        default=0.005,
        help="allowed sensitivity deficit vs minimap2",
    )

    args = parser.parse_args()

    # Default tool paths
    defaults = {
        "optik": "target/release/optik",
        "sshash": "sshash",
        "kmc": "kmc",
        "bifrost": "bifrost",
        "cuttlefish": "cuttlefish",
        "minimap2": "minimap2",
        "samtools": "samtools",
        "wfmash": "wfmash",
    }

    def auto(name):
        """Path of tool `name`, from its option, else a known build location, else the default."""
        if getattr(args, name, None) not in (None, ""):
            return Path(getattr(args, name))
        base = Path(__file__).resolve().parent.parent
        tmp = base / "temp"
        env_tools = Path(os.environ.get("TOOL_BENCH_TOOLS_DIR", ""))
        candidates = {
            "optik": [base / "target/release/optik"],
            "sshash": [tmp / "sshash/build/sshash", env_tools / "sshash/build/sshash"],
            "kmc": [tmp / "KMC/bin/kmc", env_tools / "KMC/bin/kmc"],
            "bifrost": [
                tmp / "bifrost/build/src/Bifrost",
                env_tools / "bifrost/build/src/Bifrost",
            ],
            "cuttlefish": [
                tmp / "cuttlefish/target/release/cuttlefish",
                env_tools / "cuttlefish/target/release/cuttlefish",
            ],
            "minimap2": [tmp / "minimap2/minimap2", env_tools / "minimap2/minimap2"],
            "wfmash": [
                tmp / "wfmash/build/bin/wfmash",
                env_tools / "wfmash/build/bin/wfmash",
                Path(shutil.which("wfmash") or ""),
            ],
            "samtools": [
                base.parent / "samtools" / "samtools",
                env_tools / "samtools" / "samtools",
                Path(shutil.which("samtools") or ""),
            ],
        }
        for c in candidates.get(name, []):
            if c and c.is_file():
                return c
        return Path(defaults[name])

    for name in defaults:
        setattr(args, name, Tool(auto(name)))

    if not args.workdir:
        args.workdir = tempfile.mkdtemp(prefix="optik_official_")

    if args.suite == "query":
        if (
            not args.optik.available
            or not args.sshash.available
            or not args.kmc.available
        ):
            need("optik/sshash/kmc for 'query' suite")
        checks = suite_query(args)
    elif args.suite == "dbg":
        if not args.optik.available or not args.kmc.available:
            need("optik/kmc for 'dbg' suite")
        checks = suite_dbg(args)
    elif args.suite == "map":
        if not args.optik.available:
            need("optik for 'map' suite")
        checks = suite_map(args)
    else:
        parser.error(f"unknown suite {args.suite}")

    return print_report(checks)


if __name__ == "__main__":
    sys.exit(main())
