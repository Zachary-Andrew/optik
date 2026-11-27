#!/usr/bin/env python3
"""Validate optik's de Bruijn graph against an independent oracle, Bifrost and Cuttlefish.

Every check compares with the oracle that validate_unitigs.oracle builds from the raw input, and the official tools run in the mode that matches optik's cutoff, with their actual output files compared k-mer for k-mer.

Checks
  optik_build          optik dbruijn runs and writes a GFA
  gfa_compliance       the GFA has a VN:Z:1.0 header, well-formed S lines with ACGT sequences, and only known record types
  optik_*              validate_unitigs' structural checks (k-mer set, duplicates, paths, compaction, maximality, exact links, tags) and, at -c 1, the exact Cuttlefish partition
  bifrost_kmer_set_equal  Bifrost's k-mer set equals its filter model's
                       oracle (-r for -c 1, -s for -c 2)
  cuttlefish_kmer_set  Cuttlefish's k-mer set equals the endpoints of the (k+1)-mer edges seen at least c times

Requirements
    - the optik binary, which runs `optik dbruijn`
    - optionally the Bifrost binary, for the cross-tool comparison
    - optionally the Cuttlefish binary, for the cross-tool comparison

Usage
    python validate_dbg.py -i reads.fa -k 31 [--bifrost /path/to/bifrost]
                                       [--cuttlefish /path/to/cuttlefish]

A JSON report goes to stdout and a human-readable summary to stderr.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_unitigs
from diagnostics import DETAIL_CHARS
from dna import read_fasta
from sam_spec import (
    DBG_MAX_K,
    GfaSegmentCol,
)

# seconds allowed per graph build, and GFA issues listed
BUILD_TIMEOUT_S = 3600
ISSUES_SHOWN = 5
LINE_PREVIEW_CHARS = 80
# build threads by default
DEFAULT_THREADS = 4
GFA_SUFFIX = ".gfa"


# Bifrost's two filters, -r keeping every k-mer and -s dropping those seen once, by the cutoff each gives.
BIFROST_CUTOFF_FLAGS = {1: "-r", 2: "-s"}


def parse_gfa_unitigs(gfa_path):
    """List the (name, sequence, tags) of every S line of a GFA file."""
    unitigs = []
    gfa_file = Path(gfa_path)
    if not gfa_file.exists() or gfa_file.stat().st_size == 0:
        return unitigs

    opener = gzip.open if str(gfa_file).endswith(".gz") else open
    mode = "rt" if str(gfa_file).endswith(".gz") else "r"
    with opener(gfa_file, mode) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith(("H", "C", "L", "P")):
                continue
            if line.startswith("S\t"):
                parts = line.split("\t")
                if len(parts) >= GfaSegmentCol.SEQ + 1:
                    name = parts[GfaSegmentCol.NAME]
                    seq = parts[GfaSegmentCol.SEQ]
                    tags = {}
                    for field in parts[GfaSegmentCol.TAGS:]:
                        if ":" in field:
                            key, _, val = field.partition(":")
                            tags[key] = val
                    unitigs.append((name, seq, tags))
    return unitigs


def validate_gfa_compliance(gfa_path):
    """Check basic GFA 1.0 compliance, returning (valid, issues)."""
    issues = []
    gfa_file = Path(gfa_path)
    if not gfa_file.exists():
        return False, ["file does not exist"]

    has_header = False
    has_segments = False
    opener = gzip.open if str(gfa_file).endswith(".gz") else open
    mode = "rt" if str(gfa_file).endswith(".gz") else "r"
    with opener(gfa_file, mode) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # The header must declare GFA 1.0, and every segment must hold one sequence of nucleotide letters.
            if fields[0] == "H":
                has_header = True
                if "VN:Z:1.0" not in line:
                    issues.append("missing VN:Z:1.0 in header")
            elif fields[0] == "S":
                has_segments = True
                if len(fields) < GfaSegmentCol.SEQ + 1:
                    issues.append(f"segment line too short: {line[:LINE_PREVIEW_CHARS]}")
                else:
                    seq = fields[GfaSegmentCol.SEQ]
                    if any(char in seq for char in "\t "):
                        issues.append(f"sequence contains whitespace: {fields[1]}")
                    valid_bases = set("ACGTNacgtn")
                    invalid = set(seq) - valid_bases
                    if invalid:
                        issues.append(f"invalid bases in {fields[1]}: {invalid}")
            # Containments and paths are allowed by GFA 1.0, but any other record type is an error.
            elif fields[0] not in ("H", "S", "C", "L", "P"):
                issues.append(f"unknown record type: {fields[0]}")

    if not has_header:
        issues.append("no header line found")
    if not has_segments:
        issues.append("no segment lines found")

    return len(issues) == 0, issues


def run_optik_dbruijn(input_fa, kmer_len, min_count, out_dir, optik_bin, threads=DEFAULT_THREADS):
    """Run optik dbruijn and return (gfa_path, stdout_text)."""
    gfa_path = os.path.join(out_dir, "optik.gfa")
    cmd = [
        optik_bin,
        "dbruijn",
        "-i",
        input_fa,
        "-k",
        str(kmer_len),
        "-c",
        str(min_count),
        "-o",
        gfa_path,
        "-s",
        gfa_path[: -len(GFA_SUFFIX)] + ".kdbg",
        "-t",
        str(threads),
    ]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=BUILD_TIMEOUT_S, check=False
    )
    if proc.returncode != 0:
        return None, None, f"optik dbruijn failed: {proc.stderr[:DETAIL_CHARS]}"
    return gfa_path, proc.stdout, None


def run_bifrost(input_fa, kmer_len, min_count, out_dir, bifrost_bin, threads=DEFAULT_THREADS):
    """Run Bifrost in the mode matching min_count, returning (gfa_path, cmd, err).

    Bifrost documents two filters.
    -r keeps every k-mer, i.e. min_count 1, and -s drops k-mers seen exactly once, i.e. min_count 2.
    Other cutoffs have no Bifrost equivalent.
    """
    mode = BIFROST_CUTOFF_FLAGS.get(min_count)
    if mode is None:
        return None, None, f"no Bifrost mode equals min_count={min_count}"
    prefix = os.path.join(out_dir, "bifrost", "graph")
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    cmd = [
        bifrost_bin,
        "build",
        mode,
        input_fa,
        "-k",
        str(kmer_len),
        "-o",
        prefix,
        "-t",
        str(threads),
    ]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=BUILD_TIMEOUT_S, check=False
    )
    gfa = prefix + ".gfa.gz"
    if not os.path.isfile(gfa) and os.path.isfile(prefix + ".gfa"):
        gfa = prefix + ".gfa"
    if proc.returncode != 0 or not os.path.isfile(gfa):
        return None, cmd, f"bifrost build failed: {proc.stderr[-DETAIL_CHARS:]}"
    return gfa, cmd, None


def canon_kmer_ints(seqs, kmer_len):
    """Canonical packed k-mers of sequences, in the oracle's encoding."""
    out = set()
    for seq in seqs:
        for kmer in validate_unitigs.windows(seq, kmer_len):
            out.add(validate_unitigs.canon(kmer, kmer_len))
    return out


def validate(args):
    """Run every check and return the report dict.

    Each check compares with the independent oracle built from the raw input (validate_unitigs.oracle).
    The official tools run in their matching mode, and their actual output files are compared k-mer for k-mer."""
    checks = {}
    kmer_len = args.k
    if kmer_len > DBG_MAX_K:
        return {"k_supported": {"pass": False, "error": "dbruijn needs k <= 31"}}
    tmpdir = tempfile.mkdtemp(prefix="vd_")
    try:
        gfa_path, _, err = run_optik_dbruijn(
            args.input, kmer_len, args.min_count, tmpdir, args.optik, args.threads
        )
        if err:
            checks["optik_build"] = {"pass": False, "error": err}
            return checks
        checks["optik_build"] = {"pass": True}
        valid, issues = validate_gfa_compliance(gfa_path)
        checks["gfa_compliance"] = {"pass": valid, "issues": issues[:ISSUES_SHOWN]}

        orc = validate_unitigs.oracle(args.input, kmer_len, args.min_count)
        solid = orc[0]
        checks["oracle"] = {
            "solid_kmers": len(solid),
            "detail": f"canonical k-mers with count >= {args.min_count}",
        }
        # Structural checks, plus the exact Cuttlefish partition settled against the oracle at c = 1
        cf_bin = (
            args.cuttlefish
            if args.cuttlefish and os.path.isfile(args.cuttlefish)
            else None
        )
        for name, check in validate_unitigs.validate(
            args.input,
            gfa_path,
            kmer_len,
            args.min_count,
            cf_bin,
            args.threads,
            orc=orc,
            kdbg=gfa_path[: -len(GFA_SUFFIX)] + ".kdbg",
        ).items():
            checks["optik_" + name if name != "cuttlefish_partition" else name] = check

        if args.bifrost and os.path.isfile(args.bifrost):
            bf_gfa, cmd, err = run_bifrost(
                args.input, kmer_len, args.min_count, tmpdir, args.bifrost, args.threads
            )
            if err and cmd is None:
                checks["bifrost_kmer_set_equal"] = {"skipped": err}
            elif err:
                checks["bifrost_kmer_set_equal"] = {
                    "pass": False,
                    "cmd": " ".join(cmd),
                    "error": err,
                }
            else:
                bifrost_kmers = canon_kmer_ints((seq for _, seq, _ in parse_gfa_unitigs(bf_gfa)), kmer_len)
                checks["bifrost_kmer_set_equal"] = {
                    "pass": bifrost_kmers == solid,
                    "cmd": " ".join(cmd),
                    "bifrost_kmers": len(bifrost_kmers),
                    "only_bifrost": len(bifrost_kmers - solid),
                    "only_oracle": len(solid - bifrost_kmers),
                    "detail": "Bifrost's unitig k-mers == oracle solid set (optik's "
                    "equality is optik_kmer_set_equal)",
                }

        if cf_bin:
            pref = os.path.join(tmpdir, "cuttlefish", "cf")
            os.makedirs(os.path.dirname(pref), exist_ok=True)
            cmd = [
                cf_bin,
                "build",
                "--read",
                "-s",
                args.input,
                "-k",
                str(kmer_len),
                "-c",
                str(args.min_count),
                "-o",
                pref,
                "-t",
                str(args.threads),
                "-w",
                os.path.dirname(pref),
            ]
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=BUILD_TIMEOUT_S, check=False
            )
            cf_fa = pref + ".fa"
            if proc.returncode != 0 or not os.path.isfile(cf_fa):
                checks["cuttlefish_kmer_set"] = {
                    "pass": False,
                    "cmd": " ".join(cmd),
                    "error": proc.stderr[-DETAIL_CHARS:],
                }
            else:
                # Cuttlefish applies the cutoff to (k+1)-mer edges, whereas optik and Bifrost apply it to k-mers, so its exact expected set is the endpoints of edges seen at least c times (dbg_oracles.filter_oracles).
                from dbg_oracles import filter_oracles

                _, edge_filter = filter_oracles(args.input, kmer_len, args.min_count)
                cf_seqs = {seq for _, raw_seq in read_fasta(cf_fa) for seq in [raw_seq.upper()]}
                tbl = str.maketrans("ACGT", "TGCA")
                cfk = set()
                for seq in cf_seqs:
                    rev = seq.translate(tbl)[::-1]
                    for i in range(len(seq) - kmer_len + 1):
                        fwd_kmer, rev_kmer = seq[i : i + kmer_len], rev[len(seq) - kmer_len - i : len(seq) - i]
                        cfk.add(min(rev_kmer, fwd_kmer))
                checks["cuttlefish_kmer_set"] = {
                    "pass": cfk == edge_filter,
                    "cmd": " ".join(cmd),
                    "cuttlefish_kmers": len(cfk),
                    "edge_filter_oracle": len(edge_filter),
                    "only_cuttlefish": len(cfk - edge_filter),
                    "only_oracle": len(edge_filter - cfk),
                    "detail": "Cuttlefish output == canonical endpoints of (k+1)-mers "
                    "seen >= c times (its documented edge filter), exactly",
                }
    finally:
        import shutil

        shutil.rmtree(tmpdir, ignore_errors=True)
    return checks


def main():
    """Validate one input and print the JSON report and summary."""
    parser = argparse.ArgumentParser(
        description="Validate optik de Bruijn graph correctness"
    )
    parser.add_argument("-i", dest="input", required=True, help="Input reads/FASTA")
    parser.add_argument("-k", type=int, required=True, help="K-mer size")
    parser.add_argument(
        "--optik", default="target/release/optik", help="Path to optik binary"
    )
    parser.add_argument(
        "--bifrost", default=None, help="Path to Bifrost binary (optional)"
    )
    parser.add_argument(
        "--cuttlefish", default=None, help="Path to Cuttlefish binary (optional)"
    )
    parser.add_argument(
        "--min-count", type=int, default=1, help="Minimum k-mer count [default: 1]"
    )
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="Thread count")
    args = parser.parse_args()

    if not os.path.isfile(args.optik):
        print(f"ERROR: optik binary not found at {args.optik}", file=sys.stderr)
        sys.exit(1)

    checks = validate(args)

    # Entries without a pass field are measurements rather than checks, so they do not count toward the verdict.
    all_pass = all(
        check.get("pass", True)
        for check in checks.values()
        if isinstance(check, dict) and "pass" in check
    )
    overall = "PASS" if all_pass else "FAIL"

    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  DBG Validation: {overall}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    # Each check prints on one line with its coverage, Jaccard index, detail or error when it has one.
    for name, check in checks.items():
        if not isinstance(check, dict) or "pass" not in check:
            continue
        status = "PASS" if check.get("pass", False) else "FAIL"
        line = f"  [{status}] {name}"
        if "coverage_pct" in check:
            line += f"  ({check['coverage_pct']}%)"
        if "jaccard_pct" in check:
            line += f"  (jaccard={check['jaccard_pct']}%)"
        if "detail" in check:
            line += f"  -- {check['detail']}"
        if check.get("error"):
            line += f"  ERROR: {check['error']}"
        print(line, file=sys.stderr)
    print(file=sys.stderr)

    report = {
        "overall": overall,
        "checks": checks,
    }
    print(json.dumps(report, indent=2))

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
