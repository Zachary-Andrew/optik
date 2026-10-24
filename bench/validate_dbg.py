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
import validate_unitigs as vu
from sam_spec import DBG_MAX_K


def read_fasta(path):
    """Read FASTA or FASTQ records, gzipped or not.

    FASTA records start with '>' and FASTQ records with '@'.
    In FASTQ the sequence section is followed by a '+' line and a quality section with exactly as many lines, which is skipped.
    The format is detected from the first header line, and empty headers get sequential names seq_0, seq_1 and so on.
    """
    sequences = []
    name = None
    parts = []
    seq_idx = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt" if str(path).endswith(".gz") else "r"
    have = False
    fastq = False
    seq_lines = 0  # sequence lines in the current section
    qual_left = 0  # quality lines still to skip
    with opener(path, mode) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if qual_left > 0:  # skip a quality line verbatim
                qual_left -= 1
                continue
            if not line:
                continue
            first = line[0]
            if first in (">", "@"):
                if have:
                    sequences.append((name, "".join(parts)))
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                fastq = first == "@"
                have = True
                seq_lines = 0
                parts = []
            elif first == "+" and have and fastq:
                qual_left = seq_lines  # mirror the sequence line count
            elif have:
                parts.append(line)
                if fastq:
                    seq_lines += 1
    if have:
        sequences.append((name, "".join(parts)))
    return sequences


def parse_gfa_unitigs(gfa_path):
    """List the (name, sequence, tags) of every S line of a GFA file."""
    unitigs = []
    p = Path(gfa_path)
    if not p.exists() or p.stat().st_size == 0:
        return unitigs

    opener = gzip.open if str(p).endswith(".gz") else open
    mode = "rt" if str(p).endswith(".gz") else "r"
    with opener(p, mode) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith(("H", "C", "L", "P")):
                continue
            if line.startswith("S\t"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    name = parts[1]
                    seq = parts[2]
                    tags = {}
                    for field in parts[3:]:
                        if ":" in field:
                            key, _, val = field.partition(":")
                            tags[key] = val
                    unitigs.append((name, seq, tags))
    return unitigs


def validate_gfa_compliance(gfa_path):
    """Check basic GFA 1.0 compliance, returning (valid, issues)."""
    issues = []
    p = Path(gfa_path)
    if not p.exists():
        return False, ["file does not exist"]

    has_header = False
    has_segments = False
    opener = gzip.open if str(p).endswith(".gz") else open
    mode = "rt" if str(p).endswith(".gz") else "r"
    with opener(p, mode) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] == "H":
                has_header = True
                if "VN:Z:1.0" not in line:
                    issues.append("missing VN:Z:1.0 in header")
            elif fields[0] == "S":
                has_segments = True
                if len(fields) < 3:
                    issues.append(f"segment line too short: {line[:80]}")
                else:
                    seq = fields[2]
                    if any(c in seq for c in "\t "):
                        issues.append(f"sequence contains whitespace: {fields[1]}")
                    valid_bases = set("ACGTNacgtn")
                    invalid = set(seq) - valid_bases
                    if invalid:
                        issues.append(f"invalid bases in {fields[1]}: {invalid}")
            elif fields[0] not in ("H", "S", "C", "L", "P"):
                issues.append(f"unknown record type: {fields[0]}")

    if not has_header:
        issues.append("no header line found")
    if not has_segments:
        issues.append("no segment lines found")

    return len(issues) == 0, issues


def run_optik_dbruijn(input_fa, k, min_count, out_dir, optik_bin, threads=4):
    """Run optik dbruijn and return (gfa_path, stdout_text)."""
    gfa_path = os.path.join(out_dir, "optik.gfa")
    cmd = [
        optik_bin,
        "dbruijn",
        "-i",
        input_fa,
        "-k",
        str(k),
        "-c",
        str(min_count),
        "-o",
        gfa_path,
        "-s",
        gfa_path[:-4] + ".kdbg",
        "-t",
        str(threads),
    ]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=3600, check=False
    )
    if proc.returncode != 0:
        return None, None, f"optik dbruijn failed: {proc.stderr[:500]}"
    return gfa_path, proc.stdout, None


def run_bifrost(input_fa, k, min_count, out_dir, bifrost_bin, threads=4):
    """Run Bifrost in the mode matching min_count, returning (gfa_path, cmd, err).

    Bifrost documents two filters.
    -r keeps every k-mer, i.e. min_count 1, and -s drops k-mers seen exactly once, i.e. min_count 2.
    Other cutoffs have no Bifrost equivalent.
    """
    mode = {1: "-r", 2: "-s"}.get(min_count)
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
        str(k),
        "-o",
        prefix,
        "-t",
        str(threads),
    ]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=3600, check=False
    )
    gfa = prefix + ".gfa.gz"
    if not os.path.isfile(gfa) and os.path.isfile(prefix + ".gfa"):
        gfa = prefix + ".gfa"
    if proc.returncode != 0 or not os.path.isfile(gfa):
        return None, cmd, f"bifrost build failed: {proc.stderr[-500:]}"
    return gfa, cmd, None


def canon_kmer_ints(seqs, k):
    """Canonical packed k-mers of sequences, in the oracle's encoding."""
    out = set()
    for s in seqs:
        for v in vu.windows(s, k):
            out.add(vu.canon(v, k))
    return out


def validate(args):
    """Run every check and return the report dict.

    Each check compares with the independent oracle built from the raw input (validate_unitigs.oracle).
    The official tools run in their matching mode, and their actual output files are compared k-mer for k-mer."""
    checks = {}
    k = args.k
    if k > DBG_MAX_K:
        return {"k_supported": {"pass": False, "error": "dbruijn needs k <= 31"}}
    tmpdir = tempfile.mkdtemp(prefix="vd_")
    try:
        gfa_path, _, err = run_optik_dbruijn(
            args.input, k, args.min_count, tmpdir, args.optik, args.threads
        )
        if err:
            checks["optik_build"] = {"pass": False, "error": err}
            return checks
        checks["optik_build"] = {"pass": True}
        valid, issues = validate_gfa_compliance(gfa_path)
        checks["gfa_compliance"] = {"pass": valid, "issues": issues[:5]}

        orc = vu.oracle(args.input, k, args.min_count)
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
        for name, c in vu.validate(
            args.input,
            gfa_path,
            k,
            args.min_count,
            cf_bin,
            args.threads,
            orc=orc,
            kdbg=gfa_path[:-4] + ".kdbg",
        ).items():
            checks["optik_" + name if name != "cuttlefish_partition" else name] = c

        if args.bifrost and os.path.isfile(args.bifrost):
            bf_gfa, cmd, err = run_bifrost(
                args.input, k, args.min_count, tmpdir, args.bifrost, args.threads
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
                bf = canon_kmer_ints((s for _, s, _ in parse_gfa_unitigs(bf_gfa)), k)
                checks["bifrost_kmer_set_equal"] = {
                    "pass": bf == solid,
                    "cmd": " ".join(cmd),
                    "bifrost_kmers": len(bf),
                    "only_bifrost": len(bf - solid),
                    "only_oracle": len(solid - bf),
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
                str(k),
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
                cmd, capture_output=True, text=True, timeout=3600, check=False
            )
            cf_fa = pref + ".fa"
            if proc.returncode != 0 or not os.path.isfile(cf_fa):
                checks["cuttlefish_kmer_set"] = {
                    "pass": False,
                    "cmd": " ".join(cmd),
                    "error": proc.stderr[-500:],
                }
            else:
                # Cuttlefish applies the cutoff to (k+1)-mer edges, whereas optik and Bifrost apply it to k-mers, so its exact expected set is the endpoints of edges seen at least c times (dbg_oracles.filter_oracles).
                from dbg_oracles import filter_oracles

                _, ef = filter_oracles(args.input, k, args.min_count)
                cf = {s for _, s0 in read_fasta(cf_fa) for s in [s0.upper()]}
                tbl = str.maketrans("ACGT", "TGCA")
                cfk = set()
                for s in cf:
                    r = s.translate(tbl)[::-1]
                    for i in range(len(s) - k + 1):
                        a, b = s[i : i + k], r[len(s) - k - i : len(s) - i]
                        cfk.add(min(b, a))
                checks["cuttlefish_kmer_set"] = {
                    "pass": cfk == ef,
                    "cmd": " ".join(cmd),
                    "cuttlefish_kmers": len(cfk),
                    "edge_filter_oracle": len(ef),
                    "only_cuttlefish": len(cfk - ef),
                    "only_oracle": len(ef - cfk),
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
    parser.add_argument("-t", "--threads", type=int, default=4, help="Thread count")
    args = parser.parse_args()

    if not os.path.isfile(args.optik):
        print(f"ERROR: optik binary not found at {args.optik}", file=sys.stderr)
        sys.exit(1)

    checks = validate(args)

    all_pass = all(
        c.get("pass", True)
        for c in checks.values()
        if isinstance(c, dict) and "pass" in c
    )
    overall = "PASS" if all_pass else "FAIL"

    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  DBG Validation: {overall}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
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
