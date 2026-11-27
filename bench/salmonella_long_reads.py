#!/usr/bin/env python3
"""Salmonella long-read suite of 50 kb reads, configuration B with injected errors.

optik map, on its chaining path with its own defaults, and the official wfmash tool at its defaults map the same simulated reads, drawn from both strands of the reference.
One rule scores both tools' placements against the simulation truth (same reference and strand, read origin within 1% of the read length), and optik must place every read and at least as many as wfmash.
optik's SAM also goes through wgatools validate and pafcov, an independent check of CIGARs and positions.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import hw
import validate_mapping
from diagnostics import DETAIL_CHARS
from sam_spec import PafCol
from units import PERCENT

# Seconds allowed for a tool run and for the read simulation.
RUN_TIMEOUT_S = 3600
SIMULATE_TIMEOUT_S = 300
# command-line defaults: reads, read length, error rate, seed
DEFAULT_READS = 30
DEFAULT_READ_LEN = 50000
DEFAULT_ERROR = 0.001
DEFAULT_SEED = 42
# a placement is correct within max(MIN_POS_TOL, read length / POS_TOL_DIV) bases
MIN_POS_TOL = 10
POS_TOL_DIV = 100
# best_hit tuples: (aligned query length, reference, strand, read origin)
BEST_REF, BEST_STRAND, BEST_ORIGIN = 1, 2, 3
SECONDS_DECIMALS = 3


def run(cmd, timeout=RUN_TIMEOUT_S, check=True):
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout, check=False
    )
    if check and proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed rc={proc.returncode}: {proc.stderr[:DETAIL_CHARS]}")
    return proc


def ensure_faidx(fasta_path):
    """Create a .fai for `fasta_path`, since wfmash needs indexes for both target and query.

    samtools is used when available.
    Otherwise the plain FASTA index format is written directly, with name, length, sequence offset, bases per line and line width including the newline."""
    if os.path.exists(fasta_path + ".fai"):
        return
    try:
        proc = subprocess.run(
            ["samtools", "faidx", fasta_path],
            capture_output=True,
            timeout=SIMULATE_TIMEOUT_S,
            check=False,
        )
        if proc.returncode == 0 and os.path.exists(fasta_path + ".fai"):
            return
    except (OSError, subprocess.TimeoutExpired):
        pass
    rows = []
    with open(fasta_path, "rb") as handle:
        data = handle.read()
    pos = 0
    data_len = len(data)
    while pos < data_len:
        if data[pos : pos + 1] != b">":
            pos += 1
            continue
        head_end = data.index(b"\n", pos)
        name = data[pos + 1 : head_end].split(b" ")[0].decode()
        seq_start = head_end + 1
        seq_end = seq_start
        while seq_end < data_len and data[seq_end : seq_end + 1] != b">":
            seq_end += 1
        block = data[seq_start:seq_end]
        lines = block.split(b"\n")
        lines = [line for line in lines if line]
        length = sum(len(line) for line in lines)
        line_bases = len(lines[0]) if lines else 0
        line_width = line_bases + 1 if lines else 0
        rows.append(f"{name}\t{length}\t{seq_start}\t{line_bases}\t{line_width}")
        pos = seq_end
    if rows:
        with open(fasta_path + ".fai", "w") as handle:
            handle.write("\n".join(rows) + "\n")


def main():
    """Simulate the reads, map them with optik and wfmash, and score both against the truth."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--ref", required=True)
    parser.add_argument("--optik", required=True)
    parser.add_argument("--wfmash", required=True)
    parser.add_argument("--wgatools", required=True)
    parser.add_argument("-n", type=int, default=DEFAULT_READS)
    parser.add_argument("-l", type=int, default=DEFAULT_READ_LEN)
    parser.add_argument("-e", type=float, default=DEFAULT_ERROR)
    parser.add_argument("-s", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of both mappers (default: the usable CPUs)",
    )
    args = parser.parse_args()

    tmp = tempfile.mkdtemp(prefix="sal50k_b_")
    reads = os.path.join(tmp, "reads.fa")
    truth_tsv = os.path.join(tmp, "truth.tsv")
    gen = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_reads.py")
    run(
        [
            sys.executable,
            gen,
            "-r",
            args.ref,
            "-n",
            str(args.n),
            "-l",
            str(args.l),
            "-e",
            str(args.e),
            "-s",
            str(args.s),
            "-o",
            reads,
            "--gt",
            truth_tsv,
        ]
    )
    truth = validate_mapping.parse_ground_truth(truth_tsv)
    tol = max(MIN_POS_TOL, args.l // POS_TOL_DIV)

    # optik map with its own defaults, where long reads take the chaining path
    sam = os.path.join(tmp, "optik.sam")
    optik_cmd = [
        args.optik,
        "map",
        "-r",
        args.ref,
        "-q",
        reads,
        "-o",
        sam,
        "-t",
        str(args.threads),
    ]
    start = time.monotonic()
    run(optik_cmd)
    optik_s = round(time.monotonic() - start, SECONDS_DECIMALS)
    recs = validate_mapping.parse_sam(Path(sam).read_text())
    optik_ok = 0
    for name, truth_rec in truth.items():
        rec = recs.get(name)
        if (
            rec
            and rec["is_mapped"]
            and rec["ref"] == truth_rec["ref_name"]
            and rec["is_reverse"] == (truth_rec["strand"] == "-")
            and abs(rec["pos"] - 1 - validate_mapping.leading_s(rec["cigar"]) - truth_rec["pos"]) <= tol
        ):
            optik_ok += 1
    mapped = sum(1 for rec in recs.values() if rec["is_mapped"])

    # wgatools round trip on optik's SAM
    ref_seqs = validate_mapping.load_ref_seqs(args.ref)
    ref_lens = {name: len(seq) for name, seq in ref_seqs.items()}
    sam_text = Path(sam).read_text()
    dropped = []
    paf_lines = validate_mapping.sam_to_paf(sam_text, ref_lens, dropped)
    paf = os.path.join(tmp, "optik.paf")
    Path(paf).write_text("\n".join(paf_lines) + "\n")
    wg_result, werr = validate_mapping.run_wgatools_paf_checks(paf, args.wgatools, tmp)
    # Every mapped record has SEQ equal to the read, a CIGAR that consumes it and an NM that recomputes.
    audit = validate_mapping.audit_records(sam_text, dict(validate_mapping.read_fasta(reads)), ref_seqs)

    # wfmash with its official defaults, writing PAF to stdout
    ensure_faidx(args.ref)
    ensure_faidx(reads)
    wfm_paf = os.path.join(tmp, "wfmash.paf")
    start = time.monotonic()
    with open(wfm_paf, "w") as handle:
        wfm_proc = subprocess.run(
            [args.wfmash, args.ref, reads, "-t", str(args.threads)],
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
            timeout=RUN_TIMEOUT_S,
            check=False,
        )
    wfmash_s = round(time.monotonic() - start, SECONDS_DECIMALS)
    wfm_lines = (
        [line for line in Path(wfm_paf).read_text().splitlines() if line]
        if wfm_proc.returncode == 0
        else []
    )
    # Per read the record with the longest query block counts, and its read origin is tstart - qstart on + and tstart - (qlen - qend) on -.
    best = {}
    for line in wfm_lines:
        fields = line.split("\t")
        qlen, q_start, q_end, strand, tname, t_start = (
            int(fields[PafCol.QLEN]),
            int(fields[PafCol.QSTART]),
            int(fields[PafCol.QEND]),
            fields[PafCol.STRAND],
            fields[PafCol.TNAME],
            int(fields[PafCol.TSTART]),
        )
        origin = t_start - q_start if strand == "+" else t_start - (qlen - q_end)
        if fields[PafCol.QNAME] not in best or q_end - q_start > best[fields[PafCol.QNAME]][0]:
            best[fields[PafCol.QNAME]] = (q_end - q_start, tname, strand, origin)
    # A wfmash read counts as correct with the right contig, strand and origin within tol.
    wfm_ok = 0
    for name, truth_rec in truth.items():
        best_hit = best.get(name)
        if (
            best_hit
            and best_hit[BEST_REF] == truth_rec["ref_name"]
            and best_hit[BEST_STRAND] == truth_rec["strand"]
            and abs(best_hit[BEST_ORIGIN] - truth_rec["pos"]) <= tol
        ):
            wfm_ok += 1

    report = {
        "config": f"{args.n}x{args.l}bp @ {args.e * PERCENT:.1f}% errors, seed {args.s}",
        "optik_cmd": " ".join(
            os.path.basename(arg) if i == 0 else arg for i, arg in enumerate(optik_cmd)
        ),
        "wfmash_cmd": "wfmash REF READS -t 8",
        "optik_mapped": mapped,
        "optik_correct": optik_ok,
        "wfmash_mapped": len(best),
        "wfmash_correct": wfm_ok,
        "optik_seconds": optik_s,
        "wfmash_seconds": wfmash_s,
        "wfmash_rc": wfm_proc.returncode,
        "wfmash_paf_records": len(wfm_lines),
        "reads": args.n,
        "position_tolerance_bp": tol,
        "wgatools_validate": wg_result,
        "wgatools_error": werr,
        "unconvertible_records": len(dropped),
        "record_audit": audit,
    }
    # optik must place every read and do at least as well as wfmash. Its records must also pass the audit and wgatools.
    passed = (
        len(truth) == args.n
        and optik_ok == args.n
        and optik_ok >= wfm_ok
        and wfm_proc.returncode == 0
        and wfm_ok > 0
        and not dropped
        and audit["records"] >= args.n
        and not (
            audit["seq_bad"]
            or audit["cigar_bad"]
            or audit["span_bad"]
            or audit["nm_bad"]
        )
        and (wg_result or {}).get("validate_total_records") == len(paf_lines)
        and (wg_result or {}).get("validate_query_invalid") == 0
        and (wg_result or {}).get("validate_target_invalid") == 0
        and (wg_result or {}).get("pafcov_covered_ref_bases", 0) > 0
    )
    report["overall"] = "PASS" if passed else "FAIL"
    print(json.dumps(report, indent=1, default=str))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
