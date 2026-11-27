#!/usr/bin/env python3
"""mapq_calibration.py - MAPQ calibration of optik map against minimap2.

The script builds a deterministic repeat-rich genome, simulates reads with bench/gen_reads.py, which writes a truth TSV, and maps them with both tools.
For each MAPQ bin it reports the number of primaries and of wrong ones next to the count their MAPQs imply, the sum of 10^(-MAPQ/10).
A calibrated mapper's wrong count stays at or below the implied count, and a placement with MAPQ >= 30 should be wrong less than 0.1% of the time.

The genome
  A random backbone with segmental duplications at 0-3% divergence, two or three copies each, plus three structures that defeat the uniqueness evidence of a seed-based mapper.
    - Tandem arrays (16-200 bp units, 4-40 copies, 0-2% divergence).
      A read inside one aligns equally well one period over, often within the same candidate window.
    - A high-copy interspersed family (300 bp, 2500 copies, 2-15% divergence, Alu-like).
      Its minimizers exceed the seed index's occurrence cap, so most copies never become candidates.
    - A satellite block (171 bp unit, 2000 copies, 1-5% divergence, alpha-satellite-like).
  Real chromosomes have all three.
  On hg38 chr1 and chr21, reads were placed wrongly at MAPQ 30-60 through exactly these until they were modelled here.

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
    PHRED_SCALE,
    SamCol,
    SamFlag,
    TruthCol,
)
from units import PERCENT

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 255)]
# A MAPQ bin is over-confident when its wrong count beats twice the implied count plus a small slack.
# This script and real_reads.py also require the Poisson tail P(X >= wrong) to be below OVERCONFIDENT_P.
OVERCONFIDENT_FACTOR = 2
OVERCONFIDENT_SLACK = 3
OVERCONFIDENT_P = 1e-3


# A MAPQ bin counts its correct (or all) placements, its wrong ones and the errors its MAPQs imply, in this order.
BIN_CORRECT = BIN_TOTAL = 0
BIN_WRONG = 1
BIN_IMPLIED = 2
# optik's own default thread count, one per two logical CPUs.
LOGICAL_CPUS_PER_THREAD = 2


def overconfident(wrong, implied):
    """Whether `wrong` placements exceed OVERCONFIDENT_FACTOR times the `implied` count of a bin's MAPQ plus OVERCONFIDENT_SLACK."""
    return wrong > OVERCONFIDENT_FACTOR * implied + OVERCONFIDENT_SLACK
# Largest share of MAPQ >= 30 placements that may be wrong (0.1%)
HQ_MAX_WRONG_RATE = 0.001
# In long reads, an indel this long that the simulation never made is false, and at most this share of primaries may contain one
BIG_INDEL_LEN = 50
MAX_FALSE_INDEL_RATE = 0.005
# A placement is correct within max(tolerance, read length / ORIGIN_TOL_DIV) bp of the truth.
# score() passes SCORE_ORIGIN_TOL as the tolerance, while sam_rows() uses DEFAULT_ORIGIN_TOL.
SCORE_ORIGIN_TOL = 20
DEFAULT_ORIGIN_TOL = 15
ORIGIN_TOL_DIV = 10
# Lines are split only through the columns read here (QNAME..CIGAR), so the rest stays one field.
SAM_SPLIT_MAX = SamCol.PNEXT

# The synthetic repeat-rich genome (see the module docstring)
DEFAULT_GENOME_SEED = 17  # the seed the MAPQ cap was fitted on
MULTI_GENOME_SEED = 303
BACKBONE_LEN = 4_000_000
NUM_DUPS = 600
DUP_LEN_RANGE = (500, 5000)  # randrange bounds of a segmental duplication's length
DUP_DIVERGENCES = [0.0, 0.002, 0.005, 0.01, 0.02, 0.03]
NUM_TANDEM_ARRAYS = 300
TANDEM_UNIT_LEN_RANGE = (16, 201)  # randrange bounds: 16-200 bp units
TANDEM_COPIES_RANGE = (4, 41)  # randrange bounds: 4-40 copies
TANDEM_DIVERGENCES = [0.0, 0.005, 0.01, 0.02]
FAMILY_UNIT_LEN = 300  # the Alu-like interspersed family
FAMILY_COPIES = 2500
FAMILY_DIVERGENCES = [0.02, 0.05, 0.10, 0.15]
SATELLITE_UNIT_LEN = 171  # the alpha-satellite-like block
SATELLITE_COPIES = 2000
SATELLITE_DIVERGENCES = [0.01, 0.02, 0.03, 0.05]
FASTA_LINE_WIDTH = 80

# --multi short-read sets: read lengths and substitution rates
MULTI_READ_LENGTHS = (100, 150, 250)
MULTI_ERROR_RATES = (0.01, 0.04, 0.08)
# insertion and deletion rates are each the substitution rate divided by this
INDEL_RATE_DIV = 10
# a short-read set's seed is its read length plus its error rate times this
SEED_ERR_SCALE = 1000
# long-read sets have a tenth of the short-read count, but at least MIN_LONG_READS
LONG_READS_DIV = 10
MIN_LONG_READS = 500
# reads per run: DEFAULT_READS, or MULTI_DEFAULT_READS per error rate under --multi when --reads is left at its default
DEFAULT_READS = 1_000_000
MULTI_DEFAULT_READS = 20000
DEFAULT_ERROR = 0.01
FASTA_EXT = ".fa"


def make_genome(path, seed=DEFAULT_GENOME_SEED, backbone=BACKBONE_LEN, dups=NUM_DUPS):
    """Write the deterministic repeat-rich genome of the given seed to path."""
    rng = random.Random(seed)
    genome = [rng.choice("ACGT") for _ in range(backbone)]
    for _ in range(dups):
        length = rng.randrange(*DUP_LEN_RANGE)
        src = rng.randrange(0, backbone - length)
        dst = rng.randrange(0, backbone - length)
        div = rng.choice(DUP_DIVERGENCES)
        seg = genome[src : src + length]
        seg = [rng.choice("ACGT") if rng.random() < div else base for base in seg]
        genome[dst : dst + length] = seg

    def mut(seq, div):
        """Copy of seq with each base replaced at random with probability div."""
        return [rng.choice("ACGT") if rng.random() < div else base for base in seq]

    # Tandem arrays
    for _ in range(NUM_TANDEM_ARRAYS):
        unit = [rng.choice("ACGT") for _ in range(rng.randrange(*TANDEM_UNIT_LEN_RANGE))]
        arr = []
        div = rng.choice(TANDEM_DIVERGENCES)
        for _ in range(rng.randrange(*TANDEM_COPIES_RANGE)):
            arr += mut(unit, div)
        dst = rng.randrange(0, backbone - len(arr))
        genome[dst : dst + len(arr)] = arr
    # High-copy interspersed family
    fam = [rng.choice("ACGT") for _ in range(FAMILY_UNIT_LEN)]
    for _ in range(FAMILY_COPIES):
        copy_seq = mut(fam, rng.choice(FAMILY_DIVERGENCES))
        dst = rng.randrange(0, backbone - len(copy_seq))
        genome[dst : dst + len(copy_seq)] = copy_seq
    # Satellite block
    sat = [rng.choice("ACGT") for _ in range(SATELLITE_UNIT_LEN)]
    block = []
    for _ in range(SATELLITE_COPIES):
        block += mut(sat, rng.choice(SATELLITE_DIVERGENCES))
    dst = rng.randrange(0, backbone - len(block))
    genome[dst : dst + len(block)] = block
    with open(path, "w") as handle:
        handle.write(">rep\n")
        seq = "".join(genome)
        handle.writelines(seq[i : i + FASTA_LINE_WIDTH] + "\n" for i in range(0, len(seq), FASTA_LINE_WIDTH))


def score(sam, tsv):
    """Per-bin (correct, wrong, implied wrong) counts of a SAM against its truth TSV.

    Returns the counts, the unmapped reads and the total reads."""
    truth = {}
    with open(tsv) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[TruthCol.NAME] != "read_name":
                truth[fields[TruthCol.NAME]] = (fields[TruthCol.REF], int(fields[TruthCol.POS]), fields[TruthCol.STRAND], int(fields[TruthCol.LEN]))
    counts = {bin_range: [0, 0, 0.0] for bin_range in BINS}  # correct, wrong, implied wrong
    seen, unmapped = set(), 0
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            rec_fields = line.split("\t")
            flag = int(rec_fields[SamCol.FLAG])
            # Only the first primary of each read is scored.
            if flag & SamFlag.SECONDARY_SUPPLEMENTARY or rec_fields[SamCol.QNAME] in seen:
                continue
            seen.add(rec_fields[SamCol.QNAME])
            if flag & SamFlag.UNMAPPED:
                unmapped += 1
                continue
            ref, true_pos, strand, length = truth[rec_fields[SamCol.QNAME]]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", rec_fields[SamCol.CIGAR])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                rec_fields[SamCol.RNAME] == ref
                and ((flag & SamFlag.REVERSE) > 0) == (strand == "-")
                and abs(int(rec_fields[SamCol.POS]) - 1 - lead - true_pos) <= max(SCORE_ORIGIN_TOL, length // ORIGIN_TOL_DIV)
            )
            mapq = int(rec_fields[SamCol.MAPQ])
            for bin_range in BINS:
                if bin_range[0] <= mapq <= bin_range[1]:
                    counts[bin_range][BIN_CORRECT if good else BIN_WRONG] += 1
                    # Each read adds its MAPQ's error probability, so the bin's implied wrong count is their sum.
                    counts[bin_range][BIN_IMPLIED] += PHRED_SCALE ** (-mapq / PHRED_SCALE)
    return counts, unmapped, len(truth)


def poisson_tail(threshold, lam):
    """P(X >= w) for X ~ Poisson(lam)."""
    if threshold <= 0:
        return 1.0
    term, cdf = math.exp(-lam), 0.0
    for i in range(threshold):
        cdf += term
        term *= lam / (i + 1)
    return max(0.0, 1.0 - cdf)


def judge(label, rows, failed, long_reads=False):
    """Print wrong and implied counts per MAPQ bin for (mapq, correct) rows and record failures (see --multi).

    Returns the line as a row of the results file."""
    bins = {bin_range: [0, 0, 0.0] for bin_range in BINS}
    for mapq, correct in rows:
        for bin_range in BINS:
            if bin_range[0] <= mapq <= bin_range[1]:
                bin_counts = bins[bin_range]
                bin_counts[BIN_TOTAL] += 1
                bin_counts[BIN_WRONG] += not correct
                bin_counts[BIN_IMPLIED] += PHRED_SCALE ** (-mapq / PHRED_SCALE)
                break
    hq_total = sum(bin_counts[BIN_TOTAL] for bin_range, bin_counts in bins.items() if bin_range[0] >= HIGH_MAPQ)
    hqw = sum(bin_counts[BIN_WRONG] for bin_range, bin_counts in bins.items() if bin_range[0] >= HIGH_MAPQ)
    cells = []
    for (bin_lo, bin_hi), (count, wrong, lam) in bins.items():
        if not count:
            continue
        cells.append(f"{bin_lo}-{bin_hi}:{wrong}/{lam:.1f}" if bin_lo else f"0:{count}")
        if bin_lo > 0 and overconfident(wrong, lam) and poisson_tail(wrong, lam) < OVERCONFIDENT_P:
            failed.append(f"{label}: MAPQ {bin_lo}-{bin_hi} {wrong} wrong vs {lam:.1f} implied")
    print(
        f"  {label:34} MAPQ>=30 {hqw}/{hq_total} wrong | wrong/implied: " + "  ".join(cells)
    )
    bad = hqw > HQ_MAX_WRONG_RATE * max(hq_total, 1)
    if bad:
        failed.append(f"{label}: MAPQ>=30 {hqw}/{hq_total} wrong (> 0.1%)")
    return {
        "set": label,
        "primaries": len(rows),
        "mapq30": hq_total,
        "mapq30_wrong": hqw,
        "wrong/implied by MAPQ bin": "  ".join(cells),
    }


def sam_rows(sam, tsv, tol=DEFAULT_ORIGIN_TOL, big_indel=None):
    """(mapq, correct) of each primary record.
    With big_indel set, also the number of primaries carrying an indel of at least big_indel bases."""
    truth = {}
    with open(tsv) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[TruthCol.NAME] != "read_name":
                truth[fields[TruthCol.NAME]] = (fields[TruthCol.REF], int(fields[TruthCol.POS]), fields[TruthCol.STRAND], int(fields[TruthCol.LEN]))
    rows, big, seen = [], 0, set()
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            rec_fields = line.split("\t", SAM_SPLIT_MAX)
            flag = int(rec_fields[SamCol.FLAG])
            if flag & SamFlag.NOT_PRIMARY or rec_fields[SamCol.QNAME] in seen:
                continue
            seen.add(rec_fields[SamCol.QNAME])
            ref, true_pos, strand, length = truth[rec_fields[SamCol.QNAME]]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", rec_fields[SamCol.CIGAR])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                rec_fields[SamCol.RNAME] == ref
                and ((flag & SamFlag.REVERSE) > 0) == (strand == "-")
                and abs(int(rec_fields[SamCol.POS]) - 1 - lead - true_pos)
                <= max(tol, length // ORIGIN_TOL_DIV if big_indel else tol)
            )
            rows.append((int(rec_fields[SamCol.MAPQ]), good))
            if big_indel and any(int(op_len) >= big_indel for op_len, cigar_op in ops if cigar_op in "ID"):
                big += 1
    return rows, big


def multi(args):
    """The --multi generalization check over several genomes, read lengths, error rates and long reads."""
    refs = []
    syn = os.path.join(args.out, "repeats_s303.fa")
    if not os.path.exists(syn):
        make_genome(syn, seed=MULTI_GENOME_SEED)
    refs.append(("synthetic_s303", syn))
    for ref_path in args.refs.split(",") if args.refs else []:
        if os.path.exists(ref_path):
            refs.append((os.path.splitext(os.path.basename(ref_path))[0], os.path.abspath(ref_path)))
        else:
            print(f"  [SKIP] {ref_path}: missing")
    failed = []
    table = []
    gen = os.path.join(HERE, "gen_reads.py")
    for name, ref in refs:
        print(f"== {name}")
        for length in MULTI_READ_LENGTHS:
            reads = os.path.join(args.out, f"m_{name}_L{length}.fa")
            tsv = reads[: -len(FASTA_EXT)] + ".tsv"
            # Each read length gets one file mixing every error rate, and the names are prefixed with their rate so the truth rows stay unique.
            # The seed comes from the length and rate, so every set differs while a rerun repeats it exactly.
            if not os.path.exists(tsv):
                with open(reads, "w") as reads_out, open(tsv, "w") as truth_out:
                    truth_out.write("read_name\tref_name\tpos\tstrand\tread_length\n")
                    for err_rate in MULTI_ERROR_RATES:
                        tmp_reads, tmp_truth = reads + ".part", tsv + ".part"
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
                                str(err_rate),
                                "--ins",
                                str(err_rate / INDEL_RATE_DIV),
                                "--del",
                                str(err_rate / INDEL_RATE_DIV),
                                "-s",
                                str(length + int(err_rate * SEED_ERR_SCALE)),
                                "-o",
                                tmp_reads,
                                "--gt",
                                tmp_truth,
                            ],
                            check=True,
                            stdout=subprocess.DEVNULL,
                        )
                        with open(tmp_reads) as handle:
                            reads_out.writelines(
                                f">e{err_rate}_{line[1:]}" if line[0] == ">" else line
                                for line in handle
                            )
                        with open(tmp_truth) as handle:
                            truth_out.writelines(f"e{err_rate}_{line}" for line in list(handle)[1:])
                        os.remove(tmp_reads)
                        os.remove(tmp_truth)
            # Every set is mapped with optik's defaults, so the calibration judges the MAPQ that users get.
            sam = reads[: -len(FASTA_EXT)] + ".sam"
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
        # The long-read sets imitate HiFi (15 kb, 0.3% substitutions) and ONT (10 kb, 3% substitutions with 2% insertions and 2% deletions).
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
            # An existing truth TSV means the cached read set is complete, since gen_reads.py writes it after the reads.
            tsv = reads[: -len(FASTA_EXT)] + ".tsv"
            if not os.path.exists(tsv):
                subprocess.run(
                    [
                        sys.executable,
                        gen,
                        "-r",
                        ref,
                        "-n",
                        # Long reads are costlier to map, so fewer are simulated.
                        str(max(args.reads // LONG_READS_DIV, MIN_LONG_READS)),
                    ]
                    + opt
                    + ["-s", "77", "-o", reads, "--gt", tsv],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
            sam = reads[: -len(FASTA_EXT)] + ".sam"
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
            # A false large indel is a long-read failure that MAPQ bins cannot show.
            rows, big = sam_rows(sam, tsv, big_indel=BIG_INDEL_LEN)
            table.append(judge(f"{name} {kind}", rows, failed, long_reads=True))
            table[-1]["false_50bp_indels"] = big
            print(f"  {'':34} primaries with a false >= 50 bp indel: {big}/{len(rows)}")
            if big > MAX_FALSE_INDEL_RATE * max(len(rows), 1):
                failed.append(
                    f"{name} {kind}: {big}/{len(rows)} primaries carry a false >= 50 bp indel"
                )
    for failure in failed:
        print("  [FAIL]", failure)
    print("mapq_calibration --multi: " + ("FAILED" if failed else "PASS"))
    results_file.write_section(
        args.report,
        "mapq_multi",
        f"optik MAPQ calibration on fresh simulated reads ({args.reads} per error rate) of "
        "several genomes, read lengths and long-read kinds; a bin fails when wrong > 2 x implied + 3",
        [("", table)],
        optik=args.optik,
        notes=[f"- FAIL {failure}" for failure in failed]
        or ["All bins and MAPQ>=30 error rates PASS"],
    )
    return 1 if failed and args.gate else 0


def main():
    """Run one calibration, or the --multi check."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reads", type=int, default=DEFAULT_READS)
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, hw.cpus() // LOGICAL_CPUS_PER_THREAD),
        help="threads of every mapper (default: half the usable CPUs, optik's own default)",
    )
    parser.add_argument("--out", default=os.path.join(ROOT, "temp", "mapq_calibration"))
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    parser.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    parser.add_argument("--error", type=float, default=DEFAULT_ERROR, help="read substitution rate")
    parser.add_argument(
        "--ref", help="map to this reference instead of the synthetic genome"
    )
    parser.add_argument(
        "--genome-seed",
        type=int,
        default=DEFAULT_GENOME_SEED,
        help="seed of the synthetic repeat-rich genome (17 was used to fit the MAPQ cap)",
    )
    parser.add_argument(
        "--gate",
        action="store_true",
        help="exit 1 if optik has an overconfident MAPQ bin (wrong > 2 x implied + 3, "
        "Poisson tail < 1e-3) or > 0.1%% MAPQ>=30 errors",
    )
    parser.add_argument(
        "--multi",
        action="store_true",
        help="fresh reads of several genomes x read lengths x error rates, short and long reads (see above)",
    )
    parser.add_argument(
        "--refs",
        help="--multi: comma-separated real genomes to add to the synthetic one",
    )
    results_file.add_argument(parser)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.multi:
        if args.reads == DEFAULT_READS:
            args.reads = MULTI_DEFAULT_READS
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
                str(args.error / INDEL_RATE_DIV),
                "--del",
                str(args.error / INDEL_RATE_DIV),
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
        counts, unmapped, num_reads = score(sam, tsv)
        report[name] = {
            "unmapped": unmapped,
            "reads": num_reads,
            "bins": {f"{bin_lo}-{bin_hi}": bin_counts for (bin_lo, bin_hi), bin_counts in counts.items()},
        }
        print(f"{name}: {num_reads} reads, {unmapped} unmapped")
        # A bin is over-confident when its wrong reads exceed the expectation from their MAPQs (the sum of 10^(-MAPQ/10)) by more than 2x + 3 with a Poisson tail below 1e-3, the --multi test.
        # Comparing with the bin's lowest MAPQ and no significance test once flagged a bin over one wrong read in 40k.
        print(
            f"  {'MAPQ':>7} {'reads':>8} {'wrong':>7} {'implied':>9} {'empirical':>10}"
        )
        for (bin_lo, bin_hi), (num_ok, bad, lam) in counts.items():
            tot = num_ok + bad
            emp = bad / tot if tot else 0.0
            flag = (
                "  OVERCONFIDENT"
                if bin_lo > 0 and overconfident(bad, lam) and poisson_tail(bad, lam) < OVERCONFIDENT_P
                else ""
            )
            if flag and name == "optik":
                failed.append(
                    f"optik MAPQ {bin_lo}-{bin_hi}: {bad}/{tot} wrong vs {lam:.1f} implied"
                )
            print(
                f"  {bin_lo:>3}-{bin_hi:<3} {tot:8d} {bad:7d} {lam:9.1f} {PERCENT * emp:9.3f}%{flag}"
            )
            table.append(
                {
                    "mapper": name,
                    "MAPQ": f"{bin_lo}-{bin_hi}",
                    "reads": tot,
                    "wrong": bad,
                    "implied": round(lam, 1),
                    "empirical_%": round(PERCENT * emp, 3),
                    "overconfident": "yes" if flag else "",
                }
            )
        hq_counts = [(num_ok, bad) for (bin_lo, _), (num_ok, bad, _) in counts.items() if bin_lo >= HIGH_MAPQ]
        hq_wrong = sum(num_bad for _, num_bad in hq_counts)
        hq_n = sum(num_ok + num_bad for num_ok, num_bad in hq_counts)
        print(
            f"  MAPQ>=30: {hq_n} reads, {hq_wrong} wrong ({PERCENT * hq_wrong / max(hq_n, 1):.4f}%)"
        )
        report[name]["mapq30"] = {"reads": hq_n, "wrong": hq_wrong}
        if name == "optik" and hq_wrong > HQ_MAX_WRONG_RATE * max(hq_n, 1):
            failed.append(f"optik MAPQ>=30: {hq_wrong}/{hq_n} wrong (> 0.1%)")
    # The JSON keeps every tool's full report, while the results file gets the summary table.
    with open(os.path.join(args.out, f"report_{tag}.json"), "w") as handle:
        json.dump(report, handle, indent=1)
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
                    {"mapper": mapper, **mapper_report["mapq30"], "unmapped": mapper_report["unmapped"]}
                    for mapper, mapper_report in report.items()
                ],
            ),
        ],
        optik=args.optik,
        notes=[f"- FAIL {failure}" for failure in failed],
    )
    if args.gate:
        for failure in failed:
            print("  [FAIL]", failure)
        print("mapq_calibration: " + ("FAILED" if failed else "PASS"))
        return 1 if failed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
