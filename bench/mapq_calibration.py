#!/usr/bin/env python3
"""mapq_calibration.py - MAPQ calibration of optik map against minimap2.

The script builds a deterministic repeat-rich genome, simulates reads with bench/gen_reads.py, which writes a truth TSV, and maps them with both tools.
For each MAPQ bin it reports the number of primaries and of wrong ones next to the count their MAPQs imply, the sum of 10^(-MAPQ/10).
A calibrated mapper's wrong count stays at or below the implied count, and a placement with MAPQ >= 30 should be wrong less than 0.1% of the time.

The genome
  A random backbone with segmental duplications at 0-3% divergence, two or
  three copies each, plus three structures that defeat the uniqueness
  evidence of a seed-based mapper.
    - Tandem arrays (16-200 bp units, 4-40 copies, 0-2% divergence).
      A read inside one aligns equally well one period over, often within the same candidate window.
    - A high-copy interspersed family (300 bp, 2500 copies, 2-15% divergence, Alu-like).
      Its minimizers exceed the seed index's occurrence cap, so most copies never become candidates.
    - A satellite block (171 bp unit, 2000 copies, 1-5% divergence, alpha-satellite-like).
  Real chromosomes have all three.  On hg38 chr1 and chr21, reads were placed
  wrongly at MAPQ 30-60 through exactly these until they were modelled here.

--gate exits non-zero if any optik MAPQ bin (MAPQ > 0) has more than
2 x implied + 3 wrong reads with a Poisson tail below 1e-3, or if optik's error rate at MAPQ >= 30 exceeds 0.1%.

--multi is the generalization check.  The constants of the MAPQ model were
fitted on simulated reads of 14 genomes (E. coli, Salmonella, four synthetic repeat genomes, hg38 chr1, chr2, chr7, chr17, chr19, chr21, chrX and the whole assembly).
--multi maps fresh reads, from seeds never used for fitting, of the repeat genome of seed 303, itself one of the fitted genomes, and of every --refs genome.
Short reads are 100, 150 and 250 bp at 1, 4 and 8% error, and long reads are HiFi-like (15 kb, 0.5%) and ONT-like (10 kb, 7%).
Each bin's wrong count is compared with the sum of 10^(-MAPQ/10) over its reads, and a bin fails when wrong > 2 x implied + 3 and the Poisson tail P(X >= wrong) < 1e-3.
At MAPQ >= 30 at most 0.1% may be wrong.
Long reads' primary records must not carry an indel of 50 bp or more that the simulation never made, and more than 0.5% of records doing so fails, since a chain that jumps between tandem-repeat copies writes false indels the size of the unit.

Usage: python3 bench/mapq_calibration.py [--reads N] [--error E] [--threads T]
           [--genome-seed S | --ref REF.fa] [--gate]
           [--out DIR] [--optik BIN] [--minimap2 BIN]
       python3 bench/mapq_calibration.py --multi [--refs A.fa,B.fa] [--reads N] [--gate]
"""

import argparse
import json
import math
import os
import random
import re
import subprocess
import sys

import hw
import report as results_file
from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_NOT_PRIMARY,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 255)]


def make_genome(path, seed=17, backbone=4_000_000, dups=600):
    """Write the deterministic repeat-rich genome of the given seed to path."""
    rng = random.Random(seed)
    g = [rng.choice("ACGT") for _ in range(backbone)]
    for _ in range(dups):
        length = rng.randrange(500, 5000)
        src = rng.randrange(0, backbone - length)
        dst = rng.randrange(0, backbone - length)
        div = rng.choice([0.0, 0.002, 0.005, 0.01, 0.02, 0.03])
        seg = g[src : src + length]
        seg = [rng.choice("ACGT") if rng.random() < div else c for c in seg]
        g[dst : dst + length] = seg

    def mut(seq, div):
        """Copy of seq with each base replaced at random with probability div."""
        return [rng.choice("ACGT") if rng.random() < div else c for c in seq]

    # Tandem arrays
    for _ in range(300):
        unit = [rng.choice("ACGT") for _ in range(rng.randrange(16, 201))]
        arr = []
        div = rng.choice([0.0, 0.005, 0.01, 0.02])
        for _ in range(rng.randrange(4, 41)):
            arr += mut(unit, div)
        dst = rng.randrange(0, backbone - len(arr))
        g[dst : dst + len(arr)] = arr
    # High-copy interspersed family
    fam = [rng.choice("ACGT") for _ in range(300)]
    for _ in range(2500):
        cp = mut(fam, rng.choice([0.02, 0.05, 0.10, 0.15]))
        dst = rng.randrange(0, backbone - len(cp))
        g[dst : dst + len(cp)] = cp
    # Satellite block
    sat = [rng.choice("ACGT") for _ in range(171)]
    block = []
    for _ in range(2000):
        block += mut(sat, rng.choice([0.01, 0.02, 0.03, 0.05]))
    dst = rng.randrange(0, backbone - len(block))
    g[dst : dst + len(block)] = block
    with open(path, "w") as f:
        f.write(">rep\n")
        s = "".join(g)
        f.writelines(s[i : i + 80] + "\n" for i in range(0, len(s), 80))


def score(sam, tsv):
    """Per-bin (correct, wrong, implied wrong) counts of a SAM against its truth TSV.

    Returns the counts, the unmapped reads and the total reads."""
    gt = {}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != "read_name":
                gt[f[0]] = (f[1], int(f[2]), f[3], int(f[4]))
    counts = {b: [0, 0, 0.0] for b in BINS}  # correct, wrong, implied wrong
    seen, unmapped = set(), 0
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            t = line.split("\t")
            f = int(t[1])
            if f & SAM_FLAGS_SECONDARY_SUPPLEMENTARY or t[0] in seen:
                continue
            seen.add(t[0])
            if f & SAM_FLAG_UNMAPPED:
                unmapped += 1
                continue
            ref, p0, strand, length = gt[t[0]]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                t[2] == ref
                and ((f & SAM_FLAG_REVERSE) > 0) == (strand == "-")
                and abs(int(t[3]) - 1 - lead - p0) <= max(20, length // 10)
            )
            q = int(t[4])
            for b in BINS:
                if b[0] <= q <= b[1]:
                    counts[b][0 if good else 1] += 1
                    counts[b][2] += 10 ** (-q / 10)
    return counts, unmapped, len(gt)


def poisson_tail(w, lam):
    """P(X >= w) for X ~ Poisson(lam)."""
    if w <= 0:
        return 1.0
    term, cdf = math.exp(-lam), 0.0
    for i in range(w):
        cdf += term
        term *= lam / (i + 1)
    return max(0.0, 1.0 - cdf)


def judge(label, rows, failed, long_reads=False):
    """Print wrong and implied counts per MAPQ bin for (mapq, correct) rows and record failures (see --multi).

    Returns the line as a row of the results file."""
    bins = {b: [0, 0, 0.0] for b in BINS}
    for q, ok in rows:
        for b in BINS:
            if b[0] <= q <= b[1]:
                x = bins[b]
                x[0] += 1
                x[1] += not ok
                x[2] += 10 ** (-q / 10)
                break
    hq = sum(x[0] for b, x in bins.items() if b[0] >= HIGH_MAPQ)
    hqw = sum(x[1] for b, x in bins.items() if b[0] >= HIGH_MAPQ)
    cells = []
    for (a, b), (n, w, lam) in bins.items():
        if not n:
            continue
        cells.append(f"{a}-{b}:{w}/{lam:.1f}" if a else f"0:{n}")
        if a > 0 and w > 2 * lam + 3 and poisson_tail(w, lam) < 1e-3:
            failed.append(f"{label}: MAPQ {a}-{b} {w} wrong vs {lam:.1f} implied")
    print(
        f"  {label:34} MAPQ>=30 {hqw}/{hq} wrong | wrong/implied: " + "  ".join(cells)
    )
    bad = hqw > 0.001 * max(hq, 1)
    if bad:
        failed.append(f"{label}: MAPQ>=30 {hqw}/{hq} wrong (> 0.1%)")
    return {
        "set": label,
        "primaries": len(rows),
        "mapq30": hq,
        "mapq30_wrong": hqw,
        "wrong/implied by MAPQ bin": "  ".join(cells),
    }


def sam_rows(sam, tsv, tol=15, big_indel=None):
    """(mapq, correct) of each primary record.
    With big_indel set, also the number of primaries carrying an indel of at least big_indel bases."""
    gt = {}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != "read_name":
                gt[f[0]] = (f[1], int(f[2]), f[3], int(f[4]))
    rows, big, seen = [], 0, set()
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            t = line.split("\t", 7)
            f = int(t[1])
            if f & SAM_FLAGS_NOT_PRIMARY or t[0] in seen:
                continue
            seen.add(t[0])
            ref, p0, strand, length = gt[t[0]]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                t[2] == ref
                and ((f & SAM_FLAG_REVERSE) > 0) == (strand == "-")
                and abs(int(t[3]) - 1 - lead - p0)
                <= max(tol, length // 10 if big_indel else tol)
            )
            rows.append((int(t[4]), good))
            if big_indel and any(int(n) >= big_indel for n, o in ops if o in "ID"):
                big += 1
    return rows, big


def multi(args):
    """The --multi generalization check over several genomes, read lengths, error rates and long reads."""
    refs = []
    syn = os.path.join(args.out, "repeats_s303.fa")
    if not os.path.exists(syn):
        make_genome(syn, seed=303)
    refs.append(("synthetic_s303", syn))
    for r in args.refs.split(",") if args.refs else []:
        if os.path.exists(r):
            refs.append((os.path.splitext(os.path.basename(r))[0], os.path.abspath(r)))
        else:
            print(f"  [SKIP] {r}: missing")
    failed = []
    table = []
    gen = os.path.join(HERE, "gen_reads.py")
    for name, ref in refs:
        print(f"== {name}")
        for length in (100, 150, 250):
            reads = os.path.join(args.out, f"m_{name}_L{length}.fa")
            tsv = reads[:-3] + ".tsv"
            if not os.path.exists(tsv):
                with open(reads, "w") as fo, open(tsv, "w") as ft:
                    ft.write("read_name\tref_name\tpos\tstrand\tread_length\n")
                    for e in (0.01, 0.04, 0.08):
                        tf, tt = reads + ".part", tsv + ".part"
                        subprocess.run(
                            [
                                sys.executable,
                                gen,
                                "-r",
                                ref,
                                "-n",
                                str(args.reads),
                                "-l",
                                str(length),
                                "-e",
                                str(e),
                                "--ins",
                                str(e / 10),
                                "--del",
                                str(e / 10),
                                "-s",
                                str(length + int(e * 1000)),
                                "-o",
                                tf,
                                "--gt",
                                tt,
                            ],
                            check=True,
                            stdout=subprocess.DEVNULL,
                        )
                        with open(tf) as fh:
                            fo.writelines(
                                f">e{e}_{line[1:]}" if line[0] == ">" else line
                                for line in fh
                            )
                        with open(tt) as fh:
                            ft.writelines(f"e{e}_{line}" for line in list(fh)[1:])
                        os.remove(tf)
                        os.remove(tt)
            sam = reads[:-3] + ".sam"
            with open(sam, "w") as out:
                subprocess.run(
                    [
                        args.optik,
                        "map",
                        "-r",
                        ref,
                        "-q",
                        reads,
                        "-t",
                        str(args.threads),
                    ],
                    stdout=out,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            rows, _ = sam_rows(sam, tsv)
            table.append(judge(f"{name} {length} bp", rows, failed))
        for kind, opt in (
            (
                "hifi",
                [
                    "-l",
                    "15000",
                    "--len-sd",
                    "0.3",
                    "-e",
                    "0.003",
                    "--ins",
                    "0.001",
                    "--del",
                    "0.001",
                ],
            ),
            (
                "ont",
                [
                    "-l",
                    "10000",
                    "--len-sd",
                    "0.6",
                    "-e",
                    "0.03",
                    "--ins",
                    "0.02",
                    "--del",
                    "0.02",
                ],
            ),
        ):
            reads = os.path.join(args.out, f"m_{name}_{kind}.fa")
            tsv = reads[:-3] + ".tsv"
            if not os.path.exists(tsv):
                subprocess.run(
                    [
                        sys.executable,
                        gen,
                        "-r",
                        ref,
                        "-n",
                        str(max(args.reads // 10, 500)),
                    ]
                    + opt
                    + ["-s", "77", "-o", reads, "--gt", tsv],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
            sam = reads[:-3] + ".sam"
            with open(sam, "w") as out:
                subprocess.run(
                    [
                        args.optik,
                        "map",
                        "-r",
                        ref,
                        "-q",
                        reads,
                        "-t",
                        str(args.threads),
                    ],
                    stdout=out,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            rows, big = sam_rows(sam, tsv, big_indel=50)
            table.append(judge(f"{name} {kind}", rows, failed, long_reads=True))
            table[-1]["false_50bp_indels"] = big
            print(f"  {'':34} primaries with a false >= 50 bp indel: {big}/{len(rows)}")
            if big > 0.005 * max(len(rows), 1):
                failed.append(
                    f"{name} {kind}: {big}/{len(rows)} primaries carry a false >= 50 bp indel"
                )
    for x in failed:
        print("  [FAIL]", x)
    print("mapq_calibration --multi: " + ("FAILED" if failed else "PASS"))
    results_file.write_section(
        args.report,
        "mapq_multi",
        f"optik MAPQ calibration on fresh simulated reads ({args.reads} per error rate) of "
        "several genomes, read lengths and long-read kinds; a bin fails when wrong > 2 x implied + 3",
        [("", table)],
        optik=args.optik,
        notes=[f"- FAIL {x}" for x in failed]
        or ["All bins and MAPQ>=30 error rates PASS"],
    )
    return 1 if failed and args.gate else 0


def main():
    """Run one calibration, or the --multi check."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--reads", type=int, default=1_000_000)
    ap.add_argument(
        "--threads",
        type=int,
        default=max(1, hw.cpus() // 2),
        help="threads of every mapper (default: half the usable CPUs, optik's own default)",
    )
    ap.add_argument("--out", default=os.path.join(ROOT, "temp", "mapq_calibration"))
    ap.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    ap.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    ap.add_argument("--error", type=float, default=0.01, help="read substitution rate")
    ap.add_argument(
        "--ref", help="map to this reference instead of the synthetic genome"
    )
    ap.add_argument(
        "--genome-seed",
        type=int,
        default=17,
        help="seed of the synthetic repeat-rich genome (17 was used to fit the MAPQ cap)",
    )
    ap.add_argument(
        "--gate",
        action="store_true",
        help="exit 1 if optik has an overconfident MAPQ bin (wrong > 2 x implied + 3, "
        "Poisson tail < 1e-3) or > 0.1%% MAPQ>=30 errors",
    )
    ap.add_argument(
        "--multi",
        action="store_true",
        help="fresh reads of several genomes x read lengths x error rates, short and long reads (see above)",
    )
    ap.add_argument(
        "--refs",
        help="--multi: comma-separated real genomes to add to the synthetic one",
    )
    results_file.add_argument(ap)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.multi:
        if args.reads == 1_000_000:
            args.reads = 20000
        return multi(args)
    if args.ref:
        ref = os.path.abspath(args.ref)
        gtag = os.path.splitext(os.path.basename(ref))[0]
    else:
        ref = os.path.join(args.out, f"repeats_s{args.genome_seed}.fa")
        if not os.path.exists(ref):
            make_genome(ref, seed=args.genome_seed)
        gtag = f"rep{args.genome_seed}"
    tag = f"{gtag}_e{args.error:g}_n{args.reads}"
    reads, tsv = (
        os.path.join(args.out, f"reads_{tag}.fa"),
        os.path.join(args.out, f"reads_{tag}.tsv"),
    )
    if not os.path.exists(tsv):
        subprocess.run(
            [
                sys.executable,
                os.path.join(HERE, "gen_reads.py"),
                "-r",
                ref,
                "-n",
                str(args.reads),
                "-l",
                "150",
                "-e",
                str(args.error),
                "--ins",
                str(args.error / 10),
                "--del",
                str(args.error / 10),
                "--rc",
                "0.5",
                "-s",
                "23",
                "-o",
                reads,
                "--gt",
                tsv,
            ],
            check=True,
        )
    runs = {
        "optik": [args.optik, "map", "-r", ref, "-q", reads, "-t", str(args.threads)],
        "minimap2": [args.minimap2, "-ax", "sr", "-t", str(args.threads), ref, reads],
    }
    report = {}
    failed = []
    table = []
    for name, cmd in runs.items():
        sam = os.path.join(args.out, f"{name}_{tag}.sam")
        with open(sam, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        counts, unmapped, n = score(sam, tsv)
        report[name] = {
            "unmapped": unmapped,
            "reads": n,
            "bins": {f"{a}-{b}": c for (a, b), c in counts.items()},
        }
        print(f"{name}: {n} reads, {unmapped} unmapped")
        # A bin is over-confident when its wrong reads exceed the expectation from their MAPQs (the sum of 10^(-MAPQ/10)) by more than 2x + 3 with a Poisson tail below 1e-3, the --multi test.
        # Comparing with the bin's lowest MAPQ and no significance test once flagged a bin over one wrong read in 40k.
        print(
            f"  {'MAPQ':>7} {'reads':>8} {'wrong':>7} {'implied':>9} {'empirical':>10}"
        )
        for (a, b), (ok, bad, lam) in counts.items():
            tot = ok + bad
            emp = bad / tot if tot else 0.0
            flag = (
                "  OVERCONFIDENT"
                if a > 0 and bad > 2 * lam + 3 and poisson_tail(bad, lam) < 1e-3
                else ""
            )
            if flag and name == "optik":
                failed.append(
                    f"optik MAPQ {a}-{b}: {bad}/{tot} wrong vs {lam:.1f} implied"
                )
            print(
                f"  {a:>3}-{b:<3} {tot:8d} {bad:7d} {lam:9.1f} {100 * emp:9.3f}%{flag}"
            )
            table.append(
                {
                    "mapper": name,
                    "MAPQ": f"{a}-{b}",
                    "reads": tot,
                    "wrong": bad,
                    "implied": round(lam, 1),
                    "empirical_%": round(100 * emp, 3),
                    "overconfident": "yes" if flag else "",
                }
            )
        hq = [(ok, bad) for (a, _), (ok, bad, _) in counts.items() if a >= HIGH_MAPQ]
        hq_wrong = sum(b for _, b in hq)
        hq_n = sum(o + b for o, b in hq)
        print(
            f"  MAPQ>=30: {hq_n} reads, {hq_wrong} wrong ({100 * hq_wrong / max(hq_n, 1):.4f}%)"
        )
        report[name]["mapq30"] = {"reads": hq_n, "wrong": hq_wrong}
        if name == "optik" and hq_wrong > 0.001 * max(hq_n, 1):
            failed.append(f"optik MAPQ>=30: {hq_wrong}/{hq_n} wrong (> 0.1%)")
    with open(os.path.join(args.out, f"report_{tag}.json"), "w") as f:
        json.dump(report, f, indent=1)
    results_file.write_section(
        args.report,
        f"mapq_{tag}",
        f"MAPQ calibration of optik and minimap2 on {args.reads} simulated 150 bp reads "
        f"({gtag}, error {args.error:g}): wrong reads per MAPQ bin against the number its "
        "MAPQs imply",
        [
            ("bins", table),
            (
                "MAPQ>=30",
                [
                    {"mapper": m, **r["mapq30"], "unmapped": r["unmapped"]}
                    for m, r in report.items()
                ],
            ),
        ],
        optik=args.optik,
        notes=[f"- FAIL {x}" for x in failed],
    )
    if args.gate:
        for x in failed:
            print("  [FAIL]", x)
        print("mapq_calibration: " + ("FAILED" if failed else "PASS"))
        return 1 if failed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
