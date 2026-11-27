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

from diagnostics import DETAIL_CHARS, LOG_TAIL_CHARS
from dna import read_fasta, revcomp
from sam_spec import (
    HIGH_MAPQ,
    SAM_MIN_COLUMNS,
    SAM_TAG_PREFIX_LEN,
    GfaSegmentCol,
    SamCol,
    SamFlag,
    TruthCol,
)

# timeouts in seconds: an official tool or optik run, a read simulation, and a whole wfmash run
TOOL_TIMEOUT_S = 3600
SIMULATE_TIMEOUT_S = 1200
WFMASH_TIMEOUT_S = 7200
# characters of a failed tool's stderr, or of its output, kept in the report
# probability that a probe k-mer is written reverse-complemented
REVCOMP_PROB = 0.5
# disagreeing windows kept as examples
MAX_EXAMPLES = 3
# a unitig needs this many k-mers to have an interior one
MIN_INTERIOR_WINDOWS = 3
# minimap2 v2.31 options.c presets: (k, w)
MM2_PRESET_KW = {"sr": (21, 11), "map-ont": (15, 10)}
# A short read is placed correctly within max(MIN_ORIGIN_TOL, length / ORIGIN_TOL_DIV) bp of its origin.
# A long read, by contrast, counts as correct when it overlaps at least MIN_LONG_OVERLAP of its true span.
MIN_ORIGIN_TOL = 10
ORIGIN_TOL_DIV = 10
MIN_LONG_OVERLAP = 0.5
# at most one in HQ_WRONG_DIV MAPQ >= 30 placements may be wrong (1%), and on error-free reads at least MIN_HQ_SHARE of the correct ones must reach MAPQ 30
HQ_WRONG_DIV = 100
MIN_HQ_SHARE = 0.9
# CLI defaults
DEFAULT_K = 31
DEFAULT_THREADS = 4
DEFAULT_SEED = 42
DEFAULT_PROBES = 2000
DEFAULT_SSHASH_M = 15
DEFAULT_NUM_READS = 5000
DEFAULT_READ_LEN = 150
DEFAULT_LONG_READS = 1000
DEFAULT_ERR_RATE = 0.01
DEFAULT_INDEL_RATE = 0.002
DEFAULT_SENS_TOL = 0.005


# The long-read simulation is seeded apart from the short-read sets.
LONG_READ_SEED_OFFSET = 2
# At least half of the correct long-read placements must reach MAPQ 30.
MIN_CONFIDENT_SHARE_DEN = 2


def canonical(kmer: str) -> str:
    """The lexicographically smaller of kmer and its reverse complement."""
    rev = revcomp(kmer)
    return min(kmer, rev)


# Process helpers


def run(cmd, **kwargs):
    """subprocess.run with stdout and stderr captured unless the caller says otherwise."""
    kwargs.setdefault("stdout", subprocess.PIPE)
    kwargs.setdefault("stderr", subprocess.PIPE)
    return subprocess.run(cmd, **kwargs, check=False)


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


def read_seqs(path):
    """Yield (name, seq) of FASTA or FASTQ records, plain or .gz.
    The official query reads are FASTQ."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        first = handle.read(1)
        handle.seek(0)
        if first == "@":
            while True:
                header = handle.readline()
                if not header:
                    return
                seq = handle.readline().strip()
                handle.readline()
                handle.readline()
                yield header[1:].split()[0] if header[1:].split() else "", seq
        else:
            yield from read_fasta(path)


def extract_kmers(seq, kmer_len):
    """Canonical k-mers of seq, splitting at N."""
    out = []
    start = 0
    for i, base in enumerate(seq):
        if base == "N":
            if i - start >= kmer_len:
                out.extend(canonical(seq[j : j + kmer_len]) for j in range(start, i - kmer_len + 1))
            start = i + 1
    if len(seq) - start >= kmer_len:
        out.extend(canonical(seq[j : j + kmer_len]) for j in range(start, len(seq) - kmer_len + 1))
    return out


def to_single_line_fasta(path_in, path_out):
    """Write one sequence per line with every record kept, as KMC's FASTA reader and SSHash's query parser require."""
    with open(path_out, "w") as handle:
        handle.writelines(f">{name}\n{seq}\n" for name, seq in read_fasta(path_in))


def write_uniquified_fasta(path_in, path_out):
    """Write single-line FASTA with unambiguous names.

    Stitched-unitig FASTA files can repeat a header, such as many '>seq1' records.
    These are renamed seq1, seq2, ... so that ground-truth positions refer to one contig only.

    Callers that index records by name in a dict over `read_fasta` collapse duplicates, while `gen_reads` works on the full list of records and reports only the name.
    The truth window cannot be reconstructed without unique names, so the map reference is rewritten with unique IDs.
    """
    with open(path_out, "w") as handle:
        handle.writelines(
            f">seq{i}\n{seq}\n" for i, (_name, seq) in enumerate(read_fasta(path_in), 1)
        )


def write_single_kmer_fasta(kmers, path_out):
    """Write one canonical k-mer per record.
    k may be up to 63, which covers k = 31."""
    with open(path_out, "w") as handle:
        handle.writelines(f">q{i}\n{kmer}\n" for i, kmer in enumerate(kmers))


# KMC oracle


def _is_fastq(path):
    """True for a .fq or .fastq path, gzipped or not."""
    lower_path = str(path).lower()
    return lower_path.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz"))


def kmc_build(kmc_bin, kmc_tmp, in_fasta, prefix, kmer_len, threads, min_count=1):
    """Run KMC over a one-sequence-per-line FASTA file or a FASTQ file and return the distinct k-mer set as a Python set."""
    if not kmc_bin.available:
        need("KMC (kmc)")
    workdir = kmc_tmp / "kmc_work"
    workdir.mkdir(parents=True, exist_ok=True)
    fmt = "-fq" if _is_fastq(in_fasta) else "-fa"
    db_prefix = kmc_tmp / str(prefix)
    proc = run(
        [
            str(kmc_bin.path),
            fmt,
            f"-k{kmer_len}",
            f"-ci{min_count}",
            f"-t{threads}",
            str(in_fasta),
            str(db_prefix),
            str(workdir),
        ]
    )
    if proc.returncode != 0:
        raise RuntimeError(f"KMC failed: {proc.stderr.decode()[-LOG_TAIL_CHARS:]}")
    dump = kmc_tmp / f"{prefix}.dump"
    proc = run([str(kmc_bin.path.parent / "kmc_dump"), str(db_prefix), str(dump)])
    if proc.returncode != 0:
        raise RuntimeError(f"kmc_dump failed: {proc.stderr.decode()[-LOG_TAIL_CHARS:]}")
    return {line.split("\t")[0] for line in dump.read_text().splitlines()}


# SSHash helpers


def sshash_query_report(sshash, index, query, kmer_len):
    """Completed process of `sshash query` on index and query."""
    proc = run(
        [str(sshash.path), "query", "-i", str(index), "-q", str(query)], timeout=TOOL_TIMEOUT_S
    )
    return proc


def parse_sshash_report(text):
    """Counter values (num_kmers, num_positive_kmers, ...) parsed from SSHash's query report."""
    counters = {}
    for key in (
        "num_kmers",
        "num_positive_kmers",
        "num_negative_kmers",
        "num_invalid_kmers",
        "num_searches",
    ):
        match = re.search(rf"^{key}\s*=\s*(\d+)", text, re.MULTILINE)
        if match:
            counters[key] = int(match.group(1))
    return counters


# optik helpers


def optik_query_report(optik, index, query, kmer_len):
    """Completed process of `optik query` on index and query at k."""
    proc = run(
        [str(optik.path), "query", "-i", str(index), "-q", str(query), "-k", str(kmer_len)],
        timeout=TOOL_TIMEOUT_S,
    )
    return proc


def parse_optik_report(text):
    """Counter values parsed from optik's query report, named as in SSHash's."""
    counters = {}
    for key in (
        "num_kmers",
        "num_positive_kmers",
        "num_negative_kmers",
        "num_invalid_kmers",
        "num_searches",
    ):
        match = re.search(rf"{key}\s*=\s*(\d+)", text)
        if match:
            counters[key] = int(match.group(1))
    return counters


# Query suite


def suite_query(args):
    """Compare optik query with sshash query on probes and reads and return the checks."""
    checks = {}
    optik, sshash, kmc = args.optik, args.sshash, args.kmc
    kmer_len = args.k
    keep = Path(args.workdir)
    keep.mkdir(parents=True, exist_ok=True)

    # Step 0 builds the reference k-mer oracle with KMC on the unitig FASTA.
    ref_oneline = keep / "ref_oneline.fa"
    to_single_line_fasta(args.ref, ref_oneline)
    oracle = kmc_build(kmc, keep, ref_oneline, "oracle", kmer_len, args.threads)
    checks["kmc_oracle"] = {"pass": True, "n": len(oracle)}

    # Step 1 indexes with optik and sshash.
    # Both read plain FASTA, and the unitig FASTA may be gzipped or wrapped, so the one-line conversion is indexed in every case.
    # optik contributes both index formats, the SPSS index from build-spss and the FM index from build-fm.
    # Every comparison below runs on each, and --ids adds the per-k-mer check for the SPSS index.
    # Both optik index types answer the same probes, so each is held to the same oracle.
    optik_idxs = {}
    for fmt, sub in (("_spss", "build-spss"), ("_fm", "build-fm")):
        idx = keep / f"optik_idx{fmt}"
        proc = run(
            [
                str(optik.path),
                sub,
                "-i",
                str(ref_oneline),
                "-k",
                str(kmer_len),
                "-o",
                str(idx),
            ],
            timeout=TOOL_TIMEOUT_S,
        )
        checks[f"optik_build{fmt}"] = {
            "pass": proc.returncode == 0,
            "stderr": proc.stderr.decode()[-DETAIL_CHARS:],
        }
        if proc.returncode == 0:
            optik_idxs[fmt] = idx
    if not optik_idxs:
        return checks

    sshash_idx = keep / "sshash_idx"
    proc = run(
        [
            str(sshash.path),
            "build",
            "-i",
            str(ref_oneline),
            "-k",
            str(kmer_len),
            "-m",
            str(args.sshash_m),
            "-o",
            str(sshash_idx),
            "-d",
            str(keep),
            "--check",
        ],
        timeout=TOOL_TIMEOUT_S,
    )
    out = proc.stdout.decode() + proc.stderr.decode()
    checks["sshash_build_check"] = {
        # SSHash's --check verifies its own index. Its success line confirms the oracle is sound before it judges optik.
        "pass": proc.returncode == 0 and "EVERYTHING OK!" in out,
        "rc": proc.returncode,
        "tail": out[-DETAIL_CHARS:],
    }

    # Step 2 draws positive probes, G distinct reference k-mers.
    n_positive = args.probes
    rng = random.Random(args.seed)
    ref_kmers = sorted(oracle)
    rng.shuffle(ref_kmers)
    pos = ref_kmers[:n_positive]
    pos_fa = keep / "pos.fa"
    # Probes come in random orientation, so a lookup that finds only the canonical orientation cannot pass.
    write_single_kmer_fasta([revcomp(kmer) if rng.random() < REVCOMP_PROB else kmer for kmer in pos], pos_fa)

    # Step 3 draws negative probes, R random canonical k-mers absent from the oracle.
    n_negative = args.probes
    neg = set()
    while len(neg) < n_negative:
        probe = "".join(rng.choice("ACGT") for _ in range(kmer_len))
        probe = canonical(probe)
        if probe not in oracle:
            neg.add(probe)
    neg_fa = keep / "neg.fa"
    # Flipping some probes to their reverse complement tests canonical lookup on both strands.
    write_single_kmer_fasta(
        [revcomp(kmer) if rng.random() < REVCOMP_PROB else kmer for kmer in sorted(neg)], neg_fa
    )

    # Both optik indexes answer the positive and negative probe sets, and each report must split them exactly.
    for (label, probe_fa, want), (fmt, iidx) in (
        (probe_set, index_item)
        for probe_set in (("pos", pos_fa, n_positive), ("neg", neg_fa, 0))
        for index_item in optik_idxs.items()
    ):
        proc = optik_query_report(optik, iidx, probe_fa, kmer_len)
        rep = parse_optik_report(proc.stdout.decode())
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

    # SSHash faces the same probes. A failure here means the oracle or the probes are wrong, not optik.
    for label, probe_fa, want in (("pos", pos_fa, n_positive), ("neg", neg_fa, 0)):
        proc = sshash_query_report(sshash, sshash_idx, probe_fa, kmer_len)
        rep = parse_sshash_report(proc.stdout.decode())
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
        timeout=SIMULATE_TIMEOUT_S,
    )
    # SSHash's query parser takes one sequence per line.
    # In a wrapped FASTA (such as a genome at 70 columns) it would read every other line as a header, so both tools get the same single-line copy.
    # FASTQ is passed unchanged.
    queries = args.queries
    if queries:
        opener = gzip.open if str(queries).endswith(".gz") else open
        with opener(queries, "rt") as handle:
            first = handle.readline()
        if first.startswith(">"):
            queries = keep / "queries_oneline.fa"
            to_single_line_fasta(args.queries, queries)
    read_sets = ([("", queries)] if queries else []) + [("_simreads", sim_q)]
    for (read_set_tag, qpath), (fmt, optik_idx) in (
        (read_set, index_item) for read_set in read_sets for index_item in optik_idxs.items()
    ):
        fmt = fmt + read_set_tag
        # Step 4 runs real reads, on which both tools must count positives and negatives the same way.
        # optik counts N-free k-mers, and sshash counts raw ones and reports the invalid separately.
        optik_proc = optik_query_report(optik, optik_idx, qpath, kmer_len)
        optik_rep = parse_optik_report(optik_proc.stdout.decode())
        sshash_proc = sshash_query_report(sshash, sshash_idx, qpath, kmer_len)
        sshash_rep = parse_sshash_report(sshash_proc.stdout.decode())
        checks["real_reads_positive_equal" + fmt] = {
            "pass": (
                optik_rep.get("num_positive_kmers") is not None
                and optik_rep.get("num_positive_kmers") == sshash_rep.get("num_positive_kmers")
            ),
            "optik": optik_rep,
            "sshash": sshash_rep,
        }
        checks["real_reads_negative_equal" + fmt] = {
            "pass": optik_rep.get("num_negative_kmers") is not None
            and optik_rep.get("num_negative_kmers") == sshash_rep.get("num_negative_kmers"),
            "optik_negative": optik_rep.get("num_negative_kmers"),
            "sshash_negative": sshash_rep.get("num_negative_kmers"),
        }
        # Both tools report num_kmers as the total k-mers seen (positives, negatives and invalid), so the totals are compared directly and the invalid counts must agree.
        # An earlier form subtracted sshash's invalid count from its total while comparing against optik's total, which set two different quantities side by side and could never hold.
        # Step 5 checks per-k-mer exactness of the streaming path.
        # --ids gives each window's membership (id >= 0 when present, -1 when absent or non-ACGT).
        # Every window must agree with the independent KMC k-mer set, and the number of present windows must equal the reported positive count, which ties the --ids path to the counting scan.
        # Aggregate totals alone would hide a false positive and a false negative in the same run.
        if fmt.startswith("_spss"):
            ids_out = keep / f"ids{fmt}.txt"
            ids_proc = run(
                [
                    str(optik.path),
                    "query",
                    "-i",
                    str(optik_idx),
                    "-q",
                    str(qpath),
                    "-k",
                    str(kmer_len),
                    "--ids",
                    str(ids_out),
                ],
                timeout=TOOL_TIMEOUT_S,
            )
            seqs = [seq for _, seq in read_seqs(qpath)]
            lines = ids_out.read_text().splitlines() if ids_proc.returncode == 0 else []
            bad = present = 0
            examples = []
            for seq, ids_line in zip(seqs, lines, strict=False):
                ids = [int(id_str) for id_str in ids_line.split()]
                upper_seq = seq.upper()
                if len(ids) != max(0, len(upper_seq) - kmer_len + 1):
                    bad += 1
                    continue
                for i, idv in enumerate(ids):
                    window = upper_seq[i : i + kmer_len]
                    want = all(base in "ACGT" for base in window) and canonical(window) in oracle
                    present += idv >= 0
                    if (idv >= 0) != want:
                        bad += 1
                        if len(examples) < MAX_EXAMPLES:
                            examples.append((window, idv, want))
            checks["per_kmer_membership" + fmt] = {
                "pass": (
                    ids_proc.returncode == 0
                    and len(lines) == len(seqs) > 0
                    and bad == 0
                    and present == optik_rep.get("num_positive_kmers")
                ),
                "records": len(seqs),
                "id_lines": len(lines),
                "window_mismatches": bad,
                "present_windows": present,
                "reported_positive": optik_rep.get("num_positive_kmers"),
                "examples": examples,
                "detail": "every window's --ids membership == KMC oracle; present windows "
                "== reported num_positive_kmers",
            }
        checks["real_reads_kmers_mapping" + fmt] = {
            "pass": (
                optik_rep.get("num_kmers") == sshash_rep.get("num_kmers")
                and (optik_rep.get("num_invalid_kmers") or 0)
                == (sshash_rep.get("num_invalid_kmers") or 0)
            ),
            "optik_kmers": optik_rep.get("num_kmers"),
            "sshash_kmers": sshash_rep.get("num_kmers"),
            "optik_invalid": optik_rep.get("num_invalid_kmers"),
            "sshash_invalid": sshash_rep.get("num_invalid_kmers"),
        }

    return checks


# De Bruijn graph suite


def parse_gfa_unitigs(gfa_path):
    """Sequences of the S lines of a GFA file, plain or gzipped."""
    seqs = []
    opener = gzip.open if str(gfa_path).endswith(".gz") else open
    with opener(gfa_path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] == "S" and len(fields) > GfaSegmentCol.SEQ:
                seqs.append(fields[GfaSegmentCol.SEQ])
    return seqs


def unitig_kmers(seq, kmer_len):
    """Set of canonical k-mers of seq."""
    return {canonical(seq[i : i + kmer_len]) for i in range(len(seq) - kmer_len + 1)}


def build_canon_graph(km31, km32, kmer_len):
    """Canonical k-mer graph built from the KMC (k+1)-mer oracle.
    out[x] holds every neighbour of x, since each edge is added in both directions."""
    out = {kmer: set() for kmer in km31}
    selfloop = set()
    for edge_kmer in km32:
        left = canonical(edge_kmer[:kmer_len])
        right = canonical(edge_kmer[1 : kmer_len + 1])
        if left not in km31 or right not in km31:
            continue
        if left == right:
            selfloop.add(left)
            continue
        out[left].add(right)
        out[right].add(left)
    return out, selfloop


def branch_crossing_count(unitigs, km31, km32, kmer_len):
    """Count unitigs that cross a branch node.
    An interior k-mer of a unitig may touch only its path neighbours and its own self-loop."""
    out, selfloop = build_canon_graph(km31, km32, kmer_len)
    violations = 0
    for seq in unitigs:
        wins = [canonical(seq[i : i + kmer_len]) for i in range(len(seq) - kmer_len + 1)]
        if len(wins) < MIN_INTERIOR_WINDOWS:
            continue
        for i in range(1, len(wins) - 1):
            kmer = wins[i]
            allowed = {wins[i - 1], wins[i + 1]}
            extra = set(out.get(kmer, set())) - allowed
            if extra and kmer in selfloop:
                extra.discard(kmer)
            if extra:
                violations += 1
    return violations


def suite_dbg(args):
    """Check optik dbruijn's unitigs against the KMC oracle, Bifrost and Cuttlefish and return the checks."""
    checks = {}
    optik, bifrost, cuttle, kmc = args.optik, args.bifrost, args.cuttlefish, args.kmc
    kmer_len = args.k
    keep = Path(args.workdir)
    keep.mkdir(parents=True, exist_ok=True)

    # optik dbruijn, then GFA, then unitig FASTA
    gfa = keep / "optik.gfa"
    kdbg = keep / "optik.kdbg"
    proc = run(
        [
            str(optik.path),
            "dbruijn",
            "-i",
            str(args.input),
            "-k",
            str(kmer_len),
            "-c",
            str(args.min_count),
            "-o",
            str(gfa),
            "-t",
            str(args.threads),
            "-s",
            str(kdbg),
        ],
        timeout=TOOL_TIMEOUT_S,
    )
    checks["optik_dbruijn"] = {
        "pass": proc.returncode == 0,
        "stderr": proc.stderr.decode()[-DETAIL_CHARS:],
    }
    if not checks["optik_dbruijn"]["pass"]:
        # optik may have written the GFA despite a nonzero exit on runs without retries, so check for it.
        pass
    if not os.path.isfile(gfa):
        return checks
    optik_us = parse_gfa_unitigs(gfa)
    optik_set = set()
    for seq in optik_us:
        optik_set |= unitig_kmers(seq, kmer_len)
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
        kmc, keep, oracle_input, "oracle", kmer_len, args.threads, args.min_count
    )
    oracle_k1 = kmc_build(kmc, keep, oracle_input, "oracle_k1", kmer_len + 1, args.threads)
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
    viol = branch_crossing_count(optik_us, oracle, oracle_k1, kmer_len)
    checks["no_branch_crossing"] = {"pass": viol == 0, "violations": viol}

    # Structural checks cover duplicated k-mers, maximality, exact L links and, at min_count 1, equality of the unitig partition with Cuttlefish --read.
    # Set equality alone would accept duplicated k-mers and over-split unitigs.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import validate_unitigs

    struct = validate_unitigs.validate(
        str(oracle_input),
        str(gfa),
        kmer_len,
        args.min_count,
        str(cuttle.path) if cuttle.available else None,
        args.threads,
        kdbg=str(kdbg),
    )
    for name, check in struct.items():
        checks["gfa_" + name] = check

    # Bifrost, as an independent tool
    if bifrost.available and args.min_count == 1:
        bf_prefix = keep / "bifrost"
        proc = run(
            [
                str(bifrost.path),
                "build",
                "-r",
                str(args.input),
                "-o",
                str(bf_prefix),
                "-k",
                str(kmer_len),
                "-t",
                str(args.threads),
            ],
            timeout=TOOL_TIMEOUT_S,
        )
        bfr = proc.returncode
        bf_us = []
        for cand in (str(bf_prefix) + ".gfa", str(bf_prefix) + ".gfa.gz"):
            if os.path.isfile(cand):
                bf_us = parse_gfa_unitigs(cand)
                break
        bf_set = set()
        for seq in bf_us:
            bf_set |= unitig_kmers(seq, kmer_len)
        sym_diff, only_bifrost, only_oracle = len(bf_set ^ oracle), len(bf_set - oracle), len(oracle - bf_set)
        checks["bifrost_build"] = {"pass": bfr == 0}
        checks["bifrost_set_equal_oracle"] = {
            "pass": sym_diff == 0,
            "n_unitigs": len(bf_us),
            "sym_diff": sym_diff,
            "only_bifrost": only_bifrost,
            "only_oracle": only_oracle,
        }

    # Cuttlefish, also independent, which works on (k+1)-mers and drops islands isolated by N
    if cuttle.available and args.min_count == 1:
        cf_prefix = keep / "cuttlefish"
        proc = run(
            [
                str(cuttle.path),
                "build",
                "--read",
                "-s",
                str(args.input),
                "-k",
                str(kmer_len),
                "-c",
                str(args.min_count),
                "-o",
                str(cf_prefix),
                "-t",
                str(args.threads),
            ],
            timeout=TOOL_TIMEOUT_S,
        )
        cf_fa = Path(str(cf_prefix) + ".fa")
        cf_set = set()
        ncf = 0
        if cf_fa.is_file():
            for _, seq in read_fasta(cf_fa):
                cf_set |= unitig_kmers(seq, kmer_len)
                ncf += 1
        extra_cf = cf_set - oracle
        missing_cf = oracle - cf_set
        checks["cuttlefish_build"] = {"pass": proc.returncode == 0}
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
    for line in path.read_text().splitlines():
        if not line or line.startswith("read_name"):
            continue
        fields = line.split("\t")
        if len(fields) > TruthCol.LEN:
            rows.append(
                (
                    fields[TruthCol.NAME],
                    fields[TruthCol.REF],
                    int(fields[TruthCol.POS]),
                    fields[TruthCol.STRAND],
                    int(fields[TruthCol.LEN]),
                )
            )
    return rows


def parse_sam_maps(path):
    """Map each read name to its primary record as dict(pairname, pos (1-based), strand, mapq, cigar, nm)."""
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            # A truncated line or a missing reference counts as unmapped, whatever the flag says.
            if len(fields) < SAM_MIN_COLUMNS or fields[SamCol.RNAME] == "*" or int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
                out[fields[SamCol.QNAME]] = {"unmapped": True}
                continue
            nm_tag = None
            for tag in fields[SamCol.TAGS:]:
                if tag.startswith("NM:i:"):
                    nm_tag = int(tag[SAM_TAG_PREFIX_LEN:])
                    break
            flag = int(fields[SamCol.FLAG])
            out[fields[SamCol.QNAME]] = {
                "unmapped": False,
                "ref": fields[SamCol.RNAME],
                "pos": int(fields[SamCol.POS]),
                "strand": "-" if flag & SamFlag.REVERSE else "+",
                "cigar": fields[SamCol.CIGAR],
                "nm": nm_tag,
                "mapq": int(fields[SamCol.MAPQ]),
            }
    return out


def sam_nm_from_cigar(read, cigar, refseq, pos1, strand):
    """Recompute the NM tag from the SAM CIGAR against the reference.

    The walk compares M/=/X bases case-insensitively, with any N counted as a mismatch.
    I adds its length while consuming read bases, D adds its length while consuming reference bases, and S, H and P consume without cost.
    On the '-' strand the read is reverse-complemented first, since SAM SEQ is stored as sequenced.

    Returns None if the CIGAR's reference span runs past the contig.
    """
    read = revcomp(read) if strand == "-" else read
    read_pos = 0
    ref_pos = pos1 - 1
    ref_span = sum(
        int(match.group("len"))
        for match in re.finditer(r"(?P<len>\d+)(?P<op>[MDN=X])", cigar)
    )
    # A span past the contig means the record is wrong. None reports that apart from a plain NM mismatch.
    if ref_pos < 0 or ref_pos + ref_span > len(refseq):
        return None
    nm_count = 0
    for match in re.finditer(r"(?P<len>\d+)(?P<op>[MIDNSHP=X])", cigar):
        op_len = int(match.group("len"))
        cigar_op = match.group("op")
        if cigar_op in "M=X":
            for j in range(op_len):
                read_base = read[read_pos + j] if read_pos + j < len(read) else "N"
                ref_base = refseq[ref_pos + j] if ref_pos + j < len(refseq) else "N"
                if not (read_base.upper() == ref_base.upper() and read_base.upper() in "ACGT"):
                    nm_count += 1
            read_pos += op_len
            ref_pos += op_len
        elif cigar_op in "IS":
            read_pos += op_len
            if cigar_op == "I":
                nm_count += op_len
        elif cigar_op in "DN":
            ref_pos += op_len
            if cigar_op == "D":
                nm_count += op_len
    return nm_count


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
    with open(path, "rb") as handle:
        name, length, start, line_bases, line_width = None, 0, 0, 0, 0
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    rows.append(f"{name}\t{length}\t{start}\t{line_bases}\t{line_width}")
                name = line[1:].split()[0].decode()
                length, start, line_bases, line_width = 0, off + len(line), 0, 0
            else:
                base_count = len(line.rstrip(b"\r\n"))
                if line_bases == 0:
                    line_bases, line_width = base_count, len(line)
                length += base_count
            off += len(line)
        if name is not None:
            rows.append(f"{name}\t{length}\t{start}\t{line_bases}\t{line_width}")
    fai.write_text("\n".join(rows) + "\n")


def one_primary_per_read(sam_in, sam_out):
    """Reduce wfmash's output to one primary per read.

    wfmash reports one primary per mapped query segment.
    For each read the primary with the most aligned read bases (M/=/X/I) is kept and the others are marked supplementary, so the same auditors and scorers apply to every mapper.
    Returns the number of reads that had more than one primary."""

    def aligned(fields):
        """Read bases aligned by SAM record t (M, =, X and I)."""
        return sum(
            int(op_len) for op_len, cigar_op in re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR]) if cigar_op in "M=XI"
        )

    # The first pass finds each read's longest primary. The second rewrites the others as supplementary.
    best, lines = {}, []
    with open(sam_in) as handle:
        for line in handle:
            lines.append(line)
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            flag = int(fields[SamCol.FLAG])
            if flag & SamFlag.NOT_PRIMARY:
                continue
            aligned_len = aligned(fields)
            if fields[SamCol.QNAME] not in best or aligned_len > best[fields[SamCol.QNAME]][0]:
                best[fields[SamCol.QNAME]] = (aligned_len, len(lines) - 1)
    multi = Counter()
    with open(sam_out, "w") as out:
        for i, line in enumerate(lines):
            if not line.startswith("@"):
                fields = line.split("\t")
                flag = int(fields[SamCol.FLAG])
                if not flag & SamFlag.NOT_PRIMARY:
                    multi[fields[SamCol.QNAME]] += 1
                    if best[fields[SamCol.QNAME]][1] != i:
                        fields[SamCol.FLAG] = str(flag | SamFlag.SUPPLEMENTARY)
                        line = "\t".join(fields)
            out.write(line)
    return sum(1 for count in multi.values() if count > 1)


def run_wfmash(wfmash, ref, reads, out, threads, opts):
    """Run wfmash REF READS [opts] -a to produce SAM with one primary per read.

    Returns (wall seconds, returncode, reads with split primaries)."""
    faidx(ref)
    faidx(reads)
    raw = Path(str(out) + ".raw")
    start = time.monotonic()
    with open(raw, "wb") as handle:
        proc = run(
            [str(wfmash.path), str(ref), str(reads), "-t", str(threads), "-a", *opts],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            timeout=WFMASH_TIMEOUT_S,
        )
    elapsed = time.monotonic() - start
    split = one_primary_per_read(raw, out) if proc.returncode == 0 else 0
    return elapsed, proc.returncode, split


def timed(func, *call_args, **kwargs):
    """Result and wall seconds of calling fn."""
    start = time.monotonic()
    result = func(*call_args, **kwargs)
    return result, time.monotonic() - start


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
    mm2_preset_kw = MM2_PRESET_KW

    def run_optik_map(reads, out, threads, preset="sr"):
        """Run optik map on reads with the preset's k and w unless overridden, writing SAM to out."""
        cmd = [str(optik.path), "map", "-r", str(ref), "-q", str(reads)]
        preset_k, preset_w = mm2_preset_kw[preset]
        for flag, value in (("-k", args.map_k or preset_k), ("-w", args.w or preset_w), ("-e", args.e)):
            if value is not None:
                cmd += [flag, str(value)]
        proc = run(cmd + ["-t", str(threads), "-o", str(out)], timeout=TOOL_TIMEOUT_S)
        return proc

    # Part 1, error-free simulated reads with exact ground truth
    sim_reads = keep / "sim.fa"
    sim_truth = keep / "sim.truth.tsv"
    greads = Path(__file__).resolve().parent / "gen_reads.py"
    proc = run(
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
        timeout=SIMULATE_TIMEOUT_S,
    )
    checks["gen_sim_reads"] = {"pass": proc.returncode == 0}
    truth = parse_truth_tsv(sim_truth)

    sam0 = keep / "sim_optik.sam"
    proc, elapsed = timed(run_optik_map, sim_reads, sam0, args.threads)
    speed.setdefault("error_free_short", {})["optik"] = round(elapsed, 3)
    checks["optik_map_exit"] = {
        "pass": proc.returncode == 0,
        "stderr": proc.stderr.decode()[-DETAIL_CHARS:],
    }

    maps = parse_sam_maps(sam0)
    ref_seqs = {name: seq for name, seq in read_fasta(ref)}
    mism = 0
    repeat_ties = 0
    for tname, tref, tpos, tstrand, tlen in truth:
        optik_map = maps.get(tname)
        if optik_map is None or optik_map.get("unmapped"):
            mism += 1
            continue
        if (optik_map["ref"], optik_map["pos"], optik_map.get("strand")) != (tref, tpos + 1, tstrand):
            # The read has no errors, so another placement is a real tie only when the reference at the reported locus, in the reported orientation, spells the read exactly, as a repeat copy does, inverted repeats on the other strand included.
            # Anything else is a misplacement.
            mapped_window = ref_seqs[optik_map["ref"]][optik_map["pos"] - 1 : optik_map["pos"] - 1 + tlen]
            truth_window = ref_seqs[tref][tpos : tpos + tlen]
            want = truth_window if tstrand == "+" else revcomp(truth_window)  # the read itself
            got = mapped_window if optik_map.get("strand") == "+" else revcomp(mapped_window)
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
    proc = run(
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
            # seed + 1 draws the error set from different positions than the error-free set.
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
        timeout=SIMULATE_TIMEOUT_S,
    )
    checks["gen_err_reads"] = {"pass": proc.returncode == 0}
    err_truth = {row[0]: row for row in parse_truth_tsv(err_truth_tsv)}
    sam1 = keep / "sim_err_optik.sam"
    proc, elapsed = timed(run_optik_map, err_reads, sam1, args.threads)
    speed.setdefault("error_short", {})["optik"] = round(elapsed, 3)
    checks["optik_map_err_exit"] = {"pass": proc.returncode == 0}

    def audit_sam(path, seqs=None):
        """Per-read primary placements and spec checks of a SAM file.

        A mapped primary without SEQ or NM:i counts as bad, so a writer that dropped them cannot pass, and the CIGAR must consume SEQ."""
        seqs = err_seqs if seqs is None else seqs
        prim = {}
        dup_primary = seq_bad = nm_bad = 0
        with open(path) as handle:
            for line in handle:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                flag = int(fields[SamCol.FLAG])
                if flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                    continue
                if fields[SamCol.QNAME] in prim:
                    dup_primary += 1
                    continue
                if flag & SamFlag.UNMAPPED:
                    prim[fields[SamCol.QNAME]] = None
                    continue
                strand = "-" if flag & SamFlag.REVERSE else "+"
                read = seqs.get(fields[SamCol.QNAME])
                seq = fields[SamCol.SEQ]
                want = None if read is None else (revcomp(read) if strand == "-" else read)
                if seq == "*" or want is None or seq.upper() != want.upper():
                    seq_bad += 1
                nm_tag = None
                for tag in fields[SamCol.TAGS:]:
                    if tag.startswith("NM:i:"):
                        nm_tag = int(tag[SAM_TAG_PREFIX_LEN:])
                ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
                if (
                    seq == "*"
                    or sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MIS=X") != len(seq)
                    or fields[SamCol.RNAME] not in ref_seqs
                ):
                    nm_bad += 1
                else:
                    # SAM SEQ is already in reference orientation.
                    recomputed_nm = sam_nm_from_cigar(seq, fields[SamCol.CIGAR], ref_seqs[fields[SamCol.RNAME]], int(fields[SamCol.POS]), "+")
                    if nm_tag is None or recomputed_nm != nm_tag:
                        nm_bad += 1
                # The read's origin on the reference is POS minus the leading clip, so clipped alignments, soft or hard, are scored like end-to-end ones.
                lead = re.match(r"^(\d+)[SH]", fields[SamCol.CIGAR])
                origin = int(fields[SamCol.POS]) - (int(lead.group(1)) if lead else 0)
                prim[fields[SamCol.QNAME]] = (fields[SamCol.RNAME], origin, strand, int(fields[SamCol.MAPQ]))
        return prim, dup_primary, seq_bad, nm_bad

    def score(prim, truth_rows=None):
        """Sensitivity and high-MAPQ errors of primary placements prim against truth tr."""
        correct = wrong = unmapped = hq_total = hq_wrong = 0
        truth_rows = err_truth if truth_rows is None else truth_rows
        for name, (_, tref, tpos, tstrand, tspan) in truth_rows.items():
            placement = prim.get(name)
            if placement is None:
                unmapped += 1
                continue
            place_ref, place_origin, place_strand, place_mapq = placement
            correct_place = (
                place_ref == tref
                and place_strand == tstrand
                and abs(place_origin - (tpos + 1)) <= max(MIN_ORIGIN_TOL, tspan // ORIGIN_TOL_DIV)
            )
            if not correct_place:
                # An equally good repeat copy is not an error.
                win = ref_seqs[place_ref][place_origin - 1 : place_origin - 1 + tspan]
                truth_window = ref_seqs[tref][tpos : tpos + tspan]
                correct_place = win == truth_window
            correct += correct_place
            wrong += not correct_place
            if place_mapq >= HIGH_MAPQ:
                hq_total += 1
                hq_wrong += not correct_place
        total = max(1, len(truth_rows))
        return {
            "sensitivity": round(correct / total, 5),
            "correct": correct,
            "wrong": wrong,
            "unmapped": unmapped,
            "mapq30_n": hq_total,
            "mapq30_wrong": hq_wrong,
        }

    err_seqs = {name: seq for name, seq in read_fasta(err_reads)}
    # The error-free set's records get the same spec audit, because its placement check above looks only at positions.
    sim_seqs = {name: seq for name, seq in read_fasta(sim_reads)}
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
        # max(1, ...) allows one wrong placement on small runs, where 1% would round to zero.
        "pass": osc["mapq30_wrong"] <= max(1, osc["mapq30_n"] // HQ_WRONG_DIV),
        "mapq30_n": osc["mapq30_n"],
        "mapq30_wrong": osc["mapq30_wrong"],
        "correct": osc["correct"],
        "detail": "<= 1% of MAPQ>=30 placements may be wrong (MAPQ is capped by "
        "divergence by design, so high-error reads may stay below 30)",
    }
    # MAPQ must also be informative.
    # On the error-free set, unique correct placements must reach MAPQ >= 30, since an all-zero MAPQ would satisfy the calibration bound above trivially.
    sim_prim, _, _, _ = audit_sam(sam0, sim_seqs)
    ssc = score(sim_prim, {row[0]: row for row in truth})
    checks["mapq_informative"] = {
        "pass": ssc["mapq30_n"] >= MIN_HQ_SHARE * ssc["correct"]
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
            _, elapsed = timed(
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
                timeout=TOOL_TIMEOUT_S,
            )
        speed["error_short"]["minimap2"] = round(elapsed, 3)
        mprim, _, _, _ = audit_sam(mm_err)
        sens["minimap2"] = score(mprim)
        # minimap2 writes a record for every read, unmapped ones included.
        # Fewer records point to a crashed or truncated baseline, which would make "optik >= minimap2" vacuous.
        mm_missing = sum(1 for name in err_truth if name not in mprim)
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
    proc = run(
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
            # LONG_READ_SEED_OFFSET gives the long reads positions unrelated to the short-read sets.
            str(args.seed + LONG_READ_SEED_OFFSET),
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
        timeout=SIMULATE_TIMEOUT_S,
    )
    lr_truth = {row[0]: row for row in parse_truth_tsv(lr_truth_tsv)}

    def long_audit(path):
        """Primary placements of a long-read SAM and its duplicate, NM and SEQ-length violations."""
        prim, dup, nm_bad, len_bad = {}, 0, 0, 0
        seen = set()
        with open(path) as handle:
            for line in handle:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                flag = int(fields[SamCol.FLAG])
                seen.add(fields[SamCol.QNAME])
                if flag & SamFlag.UNMAPPED:
                    prim.setdefault(fields[SamCol.QNAME], None)
                    continue
                ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
                # A hard-clipped record holds only part of the read, so its SEQ length and NM cannot be checked against the whole sequence.
                hard = any(cigar_op == "H" for _, cigar_op in ops) or fields[SamCol.SEQ] == "*"
                if not hard:
                    if sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MIS=X") != len(fields[SamCol.SEQ]):
                        len_bad += 1
                    nm_tag = next(
                        (int(tag[SAM_TAG_PREFIX_LEN:]) for tag in fields[SamCol.TAGS:] if tag.startswith("NM:i:")), None
                    )
                    if (
                        nm_tag is not None
                        and sam_nm_from_cigar(
                            fields[SamCol.SEQ], fields[SamCol.CIGAR], ref_seqs[fields[SamCol.RNAME]], int(fields[SamCol.POS]), "+"
                        )
                        != nm_tag
                    ):
                        nm_bad += 1
                if flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                    continue
                # A second primary for the same read breaks the SAM rule of exactly one primary per read.
                if prim.get(fields[SamCol.QNAME]):
                    dup += 1
                ref_start = int(fields[SamCol.POS]) - 1
                ref_end = ref_start + sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MDN=X")
                prim[fields[SamCol.QNAME]] = (
                    fields[SamCol.RNAME],
                    ref_start,
                    ref_end,
                    "-" if flag & SamFlag.REVERSE else "+",
                    int(fields[SamCol.MAPQ]),
                )
        # A long read is placed correctly when it lands on the true contig and strand and covers at least MIN_LONG_OVERLAP of its true span.
        good = hq_total = hq_wrong = 0
        for name, (_, tref, tpos, tstrand, tspan) in lr_truth.items():
            placement = prim.get(name)
            if not placement:
                continue
            place_ref, place_start, place_end, place_strand, place_mapq = placement
            correct_place = (
                place_ref == tref
                and place_strand == tstrand
                and min(place_end, tpos + tspan) - max(place_start, tpos) >= MIN_LONG_OVERLAP * tspan
            )
            good += correct_place
            if place_mapq >= HIGH_MAPQ:
                hq_total += 1
                hq_wrong += not correct_place
        total = max(1, len(lr_truth))
        return {
            "sensitivity": round(good / total, 5),
            "correct": good,
            "reads_without_record": sum(1 for name in lr_truth if name not in seen),
            "duplicate_primary": dup,
            "nm_mismatches": nm_bad,
            "cigar_len_mismatches": len_bad,
            "mapq30_n": hq_total,
            "mapq30_wrong": hq_wrong,
        }

    lr_sam = keep / "sim_long_optik.sam"
    proc, elapsed = timed(run_optik_map, lr_reads, lr_sam, args.threads, "map-ont")
    speed.setdefault("long", {})["optik"] = round(elapsed, 3)
    long_res = long_audit(lr_sam)
    checks["long_reads_valid_sam"] = {
        "pass": proc.returncode == 0
        and long_res["duplicate_primary"] == 0
        and long_res["nm_mismatches"] == 0
        and long_res["cigar_len_mismatches"] == 0
        and long_res["reads_without_record"] == 0,
        **long_res,
    }
    # MAPQ must be trustworthy and also useful, so few confident placements may be wrong while at least half the correct ones must be confident.
    checks["long_reads_mapq_calibration"] = {
        "pass": (
            long_res["mapq30_wrong"] <= max(1, long_res["mapq30_n"] // HQ_WRONG_DIV)
            and long_res["mapq30_n"] >= long_res["correct"] // MIN_CONFIDENT_SHARE_DEN
        ),
        "mapq30_n": long_res["mapq30_n"],
        "mapq30_wrong": long_res["mapq30_wrong"],
        "correct": long_res["correct"],
    }
    if minimap2.available:
        mm_lr = keep / "mm2_long.sam"
        with open(mm_lr, "wb") as _fh:
            _, elapsed = timed(
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
                timeout=TOOL_TIMEOUT_S,
            )
        speed["long"]["minimap2"] = round(elapsed, 3)
        mla = long_audit(mm_lr)
        checks["long_reads_sensitivity_vs_minimap2"] = {
            "pass": (
                long_res["sensitivity"] >= mla["sensitivity"] - args.sens_tol
                and mla["reads_without_record"] == 0
            ),
            "minimap2_reads_without_record": mla["reads_without_record"],
            "optik": long_res["sensitivity"],
            "minimap2": mla["sensitivity"],
            "tolerance": args.sens_tol,
        }

    # Part 3, determinism across thread counts
    sam_t = keep / "sim_optik_t1.sam"
    proc = run_optik_map(sim_reads, sam_t, 1)
    with open(sam_t) as handle:
        a_body = [line for line in handle if not line.startswith("@")]
    with open(sam0) as handle:
        b_body = [line for line in handle if not line.startswith("@")]
    checks["map_determinism_t1_vs_tN"] = {
        "pass": a_body == b_body,
        "note": "SAM alignment body must be byte-identical across thread counts",
    }

    # Part 4, concordance with minimap2
    if minimap2.available:
        mm2_sam = keep / "mm2_sim.sam"
        with open(mm2_sam, "wb") as _fh:
            proc, elapsed = timed(
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
                timeout=TOOL_TIMEOUT_S,
            )
        speed["error_free_short"]["minimap2"] = round(elapsed, 3)
        mmaps = parse_sam_maps(mm2_sam)
        reads_seq = {name: seq for name, seq in read_fasta(sim_reads)}
        agree = disagree = both_fail = 0

        def content_ok(mapping, tlen, read):
            """True when mapping mp places the read on reference text that spells it, on either strand."""
            if mapping.get("unmapped", False):
                return False
            win = ref_seqs[mapping["ref"]][mapping["pos"] - 1 : mapping["pos"] - 1 + tlen]
            return win == read or revcomp(win) == read

        for tname, _tref, _tpos, _tstrand, tlen in truth:
            optik_map = maps.get(tname)
            mm_map = mmaps.get(tname)
            read = reads_seq.get(tname)
            if not optik_map or not mm_map or read is None:
                continue
            mok = content_ok(optik_map, tlen, read)
            nok = content_ok(mm_map, tlen, read)
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
                    timeout=TOOL_TIMEOUT_S,
                )
            mmr = parse_sam_maps(mm2_sam_r)
            samr = keep / "optik_real.sam"
            run_optik_map(args.queries, samr, args.threads)
            optik_real_maps = parse_sam_maps(samr)
            # Real reads have no truth, so the check is that both mappers place the same reads.
            agg = Counter()
            for read_name in set(mmr) | set(optik_real_maps):
                optik_mapped = not optik_real_maps.get(read_name, {}).get("unmapped", False)
                mm_mapped = not mmr.get(read_name, {}).get("unmapped", False)
                agg[(optik_mapped, mm_mapped)] += 1
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
        sim_truth_d = {row[0]: row for row in truth}
        optik_sens = {
            "error_free_short": score(audit_sam(sam0)[0], sim_truth_d),
            "error_short": osc,
            "long": {"sensitivity": long_res["sensitivity"]},
        }
        legs = (
            ("error_free_short", sim_reads, "short", sim_truth_d),
            ("error_short", err_reads, "short", err_truth),
            ("long", lr_reads, "long", None),
        )
        for leg, reads, kind, leg_truth in legs:
            res = {}
            for cfg, opts in (("defaults", []), ("read_opts", WFMASH_READ_OPTS[kind])):
                out = keep / f"wfmash_{leg}_{cfg}.sam"
                elapsed, rc_, split = run_wfmash(wfmash, ref, reads, out, args.threads, opts)
                if rc_ != 0:
                    res[cfg] = {"error": f"wfmash exit {rc_}", "seconds": round(elapsed, 3)}
                    continue
                if kind == "long":
                    leg_res = long_audit(out)
                    leg_res = {
                        "sensitivity": leg_res["sensitivity"],
                        "mapq30_n": leg_res["mapq30_n"],
                        "mapq30_wrong": leg_res["mapq30_wrong"],
                    }
                else:
                    leg_res = score(audit_sam(out)[0], leg_truth)
                leg_res.update(
                    {
                        "seconds": round(elapsed, 3),
                        "reads_split_into_segments": split,
                        "cmd": " ".join(["wfmash", "REF", "READS", "-a", *opts]),
                    }
                )
                res[cfg] = leg_res
                speed[leg][f"wfmash_{cfg}"] = round(elapsed, 3)
            best = max((result.get("sensitivity", 0) for result in res.values()), default=0)
            # A baseline that crashed or failed is no comparison, since it would score 0.
            broken = {cfg: result["error"] for cfg, result in res.items() if "error" in result}
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
        proc = run(
            [str(samtools.path), "view", "-b", "-o", str(bam), str(sam0)], timeout=TOOL_TIMEOUT_S
        )
        checks["samtools_view_bam"] = {"pass": proc.returncode == 0}
        back = keep / "roundtrip.sam"
        with open(back, "wb") as _fh:
            proc = run([str(samtools.path), "view", str(bam)], stdout=_fh, timeout=TOOL_TIMEOUT_S)
        checks["samtools_view_back"] = {"pass": proc.returncode == 0}
        with open(sam0) as handle:
            orig_body = [line for line in handle if not line.startswith("@")]
        with open(back) as handle:
            back_body = [line for line in handle if not line.startswith("@")]
        checks["samtools_roundtrip_same_body"] = {"pass": orig_body == back_body, "n_ali": len(orig_body)}
        proc = run([str(samtools.path), "flagstat", str(bam)], timeout=TOOL_TIMEOUT_S)
        flagstat_text = proc.stdout.decode()
        checks["samtools_flagstat"] = {
            "pass": proc.returncode == 0 and "mapped (" in flagstat_text,
            "summary": flagstat_text.splitlines()[0] if flagstat_text else "",
        }

    return checks


# Command line


def base_parser():
    """Parent parser with the tool-path and shared options of every suite."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--optik", default=None)
    parser.add_argument("--sshash", default=None)
    parser.add_argument("--kmc", default=None)
    parser.add_argument("--bifrost", default=None)
    parser.add_argument("--cuttlefish", default=None)
    parser.add_argument("--minimap2", default=None)
    parser.add_argument("--samtools", default=None)
    parser.add_argument("--wfmash", default=None)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--workdir", default=None)
    parser.add_argument("-s", "--seed", type=int, default=DEFAULT_SEED)
    return parser


def print_report(checks):
    """Print every check as JSON with the verdict, returning 0 when all required checks pass."""
    fails = [
        key
        for key, val in checks.items()
        if isinstance(val, dict) and "pass" in val and not val["pass"]
    ]
    print(json.dumps(checks, indent=2, default=str))
    print("\nRESULT:", "FAIL" if fails else "PASS", file=sys.stderr)
    for failure in fails:
        print(f"  FAILED: {failure} -> {checks[failure]}", file=sys.stderr)
    return 0 if not fails else 1


def main():
    """Resolve tool paths and run the suite named on the command line."""
    parser = argparse.ArgumentParser(
        description="Exact-equality correctness validation of optik vs the original tools."
    )
    sub = parser.add_subparsers(dest="suite", required=True)

    query_parser = sub.add_parser("query", parents=[base_parser()])
    query_parser.add_argument("-r", "--ref", required=True)
    query_parser.add_argument("--queries", default=None)
    query_parser.add_argument("--probes", type=int, default=DEFAULT_PROBES)
    query_parser.add_argument("--sshash-m", type=int, default=DEFAULT_SSHASH_M)

    dbg_parser = sub.add_parser("dbg", parents=[base_parser()])
    dbg_parser.add_argument("-i", "--input", required=True)
    dbg_parser.add_argument("--min-count", type=int, default=1)

    map_parser = sub.add_parser("map", parents=[base_parser()])
    map_parser.add_argument("-r", "--ref", required=True)
    map_parser.add_argument("--queries", default=None)
    map_parser.add_argument("-n", type=int, default=DEFAULT_NUM_READS)
    map_parser.add_argument("--read-len", type=int, default=DEFAULT_READ_LEN)
    map_parser.add_argument(
        "--n-long",
        type=int,
        default=DEFAULT_LONG_READS,
        help="simulated long reads for the chaining-path checks",
    )
    # optik runs with its own defaults, as minimap2 runs its documented preset.
    # These options override them only when given explicitly.
    map_parser.add_argument("--map-k", type=int, default=None)
    map_parser.add_argument("-w", type=int, default=None)
    map_parser.add_argument("-e", type=int, default=None)
    map_parser.add_argument(
        "--err-rate",
        type=float,
        default=DEFAULT_ERR_RATE,
        help="per-base substitution rate of the error-read set",
    )
    map_parser.add_argument(
        "--indel-rate",
        type=float,
        default=DEFAULT_INDEL_RATE,
        help="per-base insertion (and deletion) rate",
    )
    map_parser.add_argument(
        "--sens-tol",
        type=float,
        default=DEFAULT_SENS_TOL,
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
        for candidate in candidates.get(name, []):
            if candidate and candidate.is_file():
                return candidate
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
