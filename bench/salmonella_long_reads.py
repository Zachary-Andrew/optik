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
import validate_mapping as vm


def run(cmd, timeout=3600, check=True):
    p = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout, check=False
    )
    if check and p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed rc={p.returncode}: {p.stderr[:500]}")
    return p


def ensure_faidx(fasta_path):
    """Create a .fai for `fasta_path`, since wfmash needs indexes for both target and query.

    samtools is used when available.
    Otherwise the plain FASTA index format is written directly, with name, length, sequence offset, bases per line and line width including the newline."""
    if os.path.exists(fasta_path + ".fai"):
        return
    try:
        p = subprocess.run(
            ["samtools", "faidx", fasta_path],
            capture_output=True,
            timeout=300,
            check=False,
        )
        if p.returncode == 0 and os.path.exists(fasta_path + ".fai"):
            return
    except (OSError, subprocess.TimeoutExpired):
        pass
    rows = []
    with open(fasta_path, "rb") as f:
        data = f.read()
    pos = 0
    n = len(data)
    while pos < n:
        if data[pos : pos + 1] != b">":
            pos += 1
            continue
        head_end = data.index(b"\n", pos)
        name = data[pos + 1 : head_end].split(b" ")[0].decode()
        seq_start = head_end + 1
        seq_end = seq_start
        while seq_end < n and data[seq_end : seq_end + 1] != b">":
            seq_end += 1
        block = data[seq_start:seq_end]
        lines = block.split(b"\n")
        lines = [ln for ln in lines if ln]
        length = sum(len(ln) for ln in lines)
        line_bases = len(lines[0]) if lines else 0
        line_width = line_bases + 1 if lines else 0
        rows.append(f"{name}\t{length}\t{seq_start}\t{line_bases}\t{line_width}")
        pos = seq_end
    if rows:
        with open(fasta_path + ".fai", "w") as f:
            f.write("\n".join(rows) + "\n")


def main():
    """Simulate the reads, map them with optik and wfmash, and score both against the truth."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ref", required=True)
    ap.add_argument("--optik", required=True)
    ap.add_argument("--wfmash", required=True)
    ap.add_argument("--wgatools", required=True)
    ap.add_argument("-n", type=int, default=30)
    ap.add_argument("-l", type=int, default=50000)
    ap.add_argument("-e", type=float, default=0.001)
    ap.add_argument("-s", type=int, default=42)
    ap.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of both mappers (default: the usable CPUs)",
    )
    args = ap.parse_args()

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
    truth = vm.parse_ground_truth(truth_tsv)
    tol = max(10, args.l // 100)

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
    t0 = time.monotonic()
    run(optik_cmd)
    optik_s = round(time.monotonic() - t0, 3)
    recs = vm.parse_sam(Path(sam).read_text())
    optik_ok = 0
    for name, g in truth.items():
        r = recs.get(name)
        if (
            r
            and r["is_mapped"]
            and r["ref"] == g["ref_name"]
            and r["is_reverse"] == (g["strand"] == "-")
            and abs(r["pos"] - 1 - vm.leading_s(r["cigar"]) - g["pos"]) <= tol
        ):
            optik_ok += 1
    mapped = sum(1 for r in recs.values() if r["is_mapped"])

    # wgatools round trip on optik's SAM
    ref_seqs = vm.load_ref_seqs(args.ref)
    ref_lens = {n: len(s) for n, s in ref_seqs.items()}
    sam_text = Path(sam).read_text()
    dropped = []
    paf_lines = vm.sam_to_paf(sam_text, ref_lens, dropped)
    paf = os.path.join(tmp, "optik.paf")
    Path(paf).write_text("\n".join(paf_lines) + "\n")
    wg, werr = vm.run_wgatools_paf_checks(paf, args.wgatools, tmp)
    # Every mapped record has SEQ equal to the read, a CIGAR that consumes it and an NM that recomputes.
    audit = vm.audit_records(sam_text, dict(vm.read_fasta(reads)), ref_seqs)

    # wfmash with its official defaults, writing PAF to stdout
    ensure_faidx(args.ref)
    ensure_faidx(reads)
    wfm_paf = os.path.join(tmp, "wfmash.paf")
    t0 = time.monotonic()
    with open(wfm_paf, "w") as fh:
        p = subprocess.run(
            [args.wfmash, args.ref, reads, "-t", str(args.threads)],
            stdout=fh,
            stderr=subprocess.PIPE,
            text=True,
            timeout=3600,
            check=False,
        )
    wfmash_s = round(time.monotonic() - t0, 3)
    wfm_lines = (
        [line for line in Path(wfm_paf).read_text().splitlines() if line]
        if p.returncode == 0
        else []
    )
    # Per read the record with the longest query block counts, and its read origin is tstart - qstart on + and tstart - (qlen - qend) on -.
    best = {}
    for line in wfm_lines:
        t = line.split("\t")
        qlen, qs, qe, strand, tname, ts = (
            int(t[1]),
            int(t[2]),
            int(t[3]),
            t[4],
            t[5],
            int(t[7]),
        )
        origin = ts - qs if strand == "+" else ts - (qlen - qe)
        if t[0] not in best or qe - qs > best[t[0]][0]:
            best[t[0]] = (qe - qs, tname, strand, origin)
    wfm_ok = 0
    for name, g in truth.items():
        b = best.get(name)
        if (
            b
            and b[1] == g["ref_name"]
            and b[2] == g["strand"]
            and abs(b[3] - g["pos"]) <= tol
        ):
            wfm_ok += 1

    report = {
        "config": f"{args.n}x{args.l}bp @ {args.e * 100:.1f}% errors, seed {args.s}",
        "optik_cmd": " ".join(
            os.path.basename(x) if i == 0 else x for i, x in enumerate(optik_cmd)
        ),
        "wfmash_cmd": "wfmash REF READS -t 8",
        "optik_mapped": mapped,
        "optik_correct": optik_ok,
        "wfmash_mapped": len(best),
        "wfmash_correct": wfm_ok,
        "optik_seconds": optik_s,
        "wfmash_seconds": wfmash_s,
        "wfmash_rc": p.returncode,
        "wfmash_paf_records": len(wfm_lines),
        "reads": args.n,
        "position_tolerance_bp": tol,
        "wgatools_validate": wg,
        "wgatools_error": werr,
        "unconvertible_records": len(dropped),
        "record_audit": audit,
    }
    ok = (
        len(truth) == args.n
        and optik_ok == args.n
        and optik_ok >= wfm_ok
        and p.returncode == 0
        and wfm_ok > 0
        and not dropped
        and audit["records"] >= args.n
        and not (
            audit["seq_bad"]
            or audit["cigar_bad"]
            or audit["span_bad"]
            or audit["nm_bad"]
        )
        and (wg or {}).get("validate_total_records") == len(paf_lines)
        and (wg or {}).get("validate_query_invalid") == 0
        and (wg or {}).get("validate_target_invalid") == 0
        and (wg or {}).get("pafcov_covered_ref_bases", 0) > 0
    )
    report["overall"] = "PASS" if ok else "FAIL"
    print(json.dumps(report, indent=1, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
