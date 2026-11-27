#!/usr/bin/env python3
"""wgsim_check.py - placement and MAPQ checked with samtools' own simulator and evaluator, so that neither the reads nor the scoring are optik's code.

Reads come from wgsim (samtools), with uniform substitution errors plus sample mutations at -r, of which a share -R are indels.
A mapper's primary records are scored by `wgsim_eval.pl alneval`, Heng Li's evaluator shipped with samtools.
A record is correct when its leftmost (forward) or rightmost (reverse) aligned base, with clipping undone, lies within 5 bp of the simulated fragment end it came from.
For each MAPQ decile the script prints the wrong count (alneval) next to the expected one (the sum of 10^(-MAPQ/10) over the bin's primaries), and the reads and errors at MAPQ >= 30.

The mappers are optik map, minimap2 -ax sr and, when --bwa-index is given, bwa mem, each run paired on R1 R2 and single-end on R1.

The evaluator has one caveat.
A read with a wgsim indel within about 5 bp of its end is soft-clipped by every mapper, and its clipped end can fall just outside the 5 bp window although the locus is right.

Usage:
  wgsim_check.py --ref REF.fa --out DIR [--pairs 500000] [--len 150]
                 [--error 0.01] [--bwa-index PREFIX] [--threads <usable CPUs>] [--warmup 1] [--runs 3]
                 [--mem SIZE] [--seed 7]
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import hw
import measure
import memcap
import report as results_file
from diagnostics import DETAIL_CHARS
from sam_spec import (
    HIGH_MAPQ,
    PHRED_SCALE,
    SamCol,
)

# MAPQ bins are deciles.
# Each bin holds [primaries, wrong (alneval), wrong implied by the MAPQs].
MAPQ_DECILE = 10
BIN_PRIMARIES = 0
BIN_WRONG = 1
BIN_IMPLIED = 2
# alneval rows read "<decile>x  wrong  total  ...", at least this many fields
ALNEVAL_MIN_FIELDS = 4
# command-line defaults: pairs, read length, error, mutation rate, indel share, outer distance and sd, alneval window, seed, measured runs
DEFAULT_PAIRS = 500000
DEFAULT_READ_LEN = 150
DEFAULT_ERROR = 0.01
DEFAULT_MUT = 0.001
DEFAULT_INDEL_FRAC = 0.15
DEFAULT_INSERT = 400
DEFAULT_INSERT_SD = 50
DEFAULT_WINDOW = 5
DEFAULT_SEED = 7
DEFAULT_RUNS = 3


# Mates are numbered as in wgsim's .r1 and .r2 output files.
MATE1, MATE2 = 1, 2


def score(sam, samtools, evaluator, window=DEFAULT_WINDOW):
    """{decile: [primaries, wrong, implied]} and (MAPQ>=30 reads, wrong) from alneval's output."""
    bins = {}
    view = subprocess.Popen(
        [samtools, "view", "-F", "0x904", sam], stdout=subprocess.PIPE, text=True
    )
    for line in view.stdout:
        mapq = int(line.split("\t", SamCol.CIGAR)[SamCol.MAPQ])
        counts = bins.setdefault(mapq // MAPQ_DECILE, [0, 0, 0.0])
        counts[BIN_PRIMARIES] += 1
        counts[BIN_IMPLIED] += PHRED_SCALE ** (-mapq / PHRED_SCALE)
    view.wait()
    view = subprocess.Popen(
        [samtools, "view", "-F", "0x900", sam], stdout=subprocess.PIPE
    )
    eval_proc = subprocess.run(
        ["perl", evaluator, "alneval", "-g", str(window), "-"],
        stdin=view.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    view.wait()
    for row in eval_proc.stdout.splitlines():
        fields = row.split()  # "06x  wrong / total  cumulative  ratio"
        if len(fields) >= ALNEVAL_MIN_FIELDS and fields[0].endswith("x"):
            bins.setdefault(int(fields[0][:-1]), [0, 0, 0.0])[BIN_WRONG] = int(fields[1])
    high_bins = [counts for decile, counts in bins.items() if decile >= HIGH_MAPQ // MAPQ_DECILE]
    return bins, (sum(counts[BIN_PRIMARIES] for counts in high_bins), sum(counts[BIN_WRONG] for counts in high_bins))


def main():
    """Simulate reads with wgsim, map them with each mapper and score them with alneval."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--ref", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--pairs", type=int, default=DEFAULT_PAIRS)
    parser.add_argument("--len", type=int, default=DEFAULT_READ_LEN)
    parser.add_argument("--error", type=float, default=DEFAULT_ERROR)
    parser.add_argument(
        "--mut", type=float, default=DEFAULT_MUT, help="wgsim -r: sample mutation rate"
    )
    parser.add_argument(
        "--indel-frac",
        type=float,
        default=DEFAULT_INDEL_FRAC,
        help="wgsim -R: share of mutations that are indels",
    )
    parser.add_argument("--insert", type=int, default=DEFAULT_INSERT, help="wgsim -d: outer distance")
    parser.add_argument("--insert-sd", type=int, default=DEFAULT_INSERT_SD, help="wgsim -s")
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help="wgsim_eval.pl alneval -g: correct within INT bp (its default 5; BWA-MEM's paper 20)",
    )
    parser.add_argument(
        "--se-both",
        action="store_true",
        help="single-end runs map R1 and R2 together as 2N reads (BWA-MEM Fig. 1, top)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--seed-kw",
        default="21,11",
        help="optik's minimizer k,w: minimap2 -ax sr's (v2.31 options.c), so both "
        "index the same seeds",
    )
    parser.add_argument(
        "--threads",
        default=str(hw.cpus()),
        help="threads of every mapper (default: the usable CPUs)",
    )
    parser.add_argument(
        "--mem",
        default=hw.default_cap_str(),
        help="resident-memory cap of each mapper ('' = uncapped; default: usable memory less 2 GiB)",
    )
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    parser.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    parser.add_argument(
        "--bwa", default=shutil.which("bwa") or os.path.join(ROOT, "temp/bwa/bwa")
    )
    parser.add_argument("--bwa-index", help="BWA index prefix (adds bwa mem)")
    parser.add_argument("--samtools", default=shutil.which("samtools"))
    parser.add_argument("--wgsim", default=shutil.which("wgsim"))
    parser.add_argument("--evaluator", default=shutil.which("wgsim_eval.pl"))
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="map again even where a SAM exists (the simulated reads are kept)",
    )
    parser.add_argument(
        "--warmup", type=int, default=1, help="warmup runs of every mapper, discarded"
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=DEFAULT_RUNS,
        help="measured runs of every mapper: the mean wall time and the peak RSS over them",
    )
    results_file.add_argument(parser)
    args = parser.parse_args()
    for name in ("samtools", "wgsim", "evaluator"):
        if not getattr(args, name):
            sys.exit(
                f"error: {name} not found (samtools ships wgsim and wgsim_eval.pl)"
            )
    os.makedirs(args.out, exist_ok=True)

    # The file tag names the simulation parameters, and non-default sample settings extend it, so cached reads are never reused under different settings.
    tag = f"l{args.len}e{args.error}"
    if (args.mut, args.indel_frac, args.insert, args.insert_sd) != (DEFAULT_MUT, DEFAULT_INDEL_FRAC, DEFAULT_INSERT, DEFAULT_INSERT_SD):
        tag += f"r{args.mut}R{args.indel_frac}d{args.insert}s{args.insert_sd}"
    reads1, reads2 = (os.path.join(args.out, f"{tag}.r{mate}.fq") for mate in (MATE1, MATE2))
    # An existing second mate file is taken to be a finished simulation from an earlier run.
    if not os.path.exists(reads2):
        with open(os.devnull, "w") as devnull:
            subprocess.run(
                [
                    args.wgsim,
                    "-S",
                    str(args.seed),
                    "-N",
                    str(args.pairs),
                    "-1",
                    str(args.len),
                    "-2",
                    str(args.len),
                    "-d",
                    str(args.insert),
                    "-s",
                    str(args.insert_sd),
                    "-e",
                    str(args.error),
                    "-r",
                    str(args.mut),
                    "-R",
                    str(args.indel_frac),
                    args.ref,
                    reads1,
                    reads2,
                ],
                stdout=devnull,
                stderr=devnull,
                check=True,
            )

    seed_k, seed_w = args.seed_kw.split(",")
    tools = {
        "optik": lambda reads: (
            [args.optik, "map", "-r", args.ref, "-k", seed_k, "-w", seed_w, "-t", args.threads]
            + [arg for read_path in reads for arg in ("-q", read_path)]
        ),
        "minimap2": lambda reads: (
            [args.minimap2, "-ax", "sr", "-t", args.threads, args.ref] + reads
        ),
    }
    if args.bwa_index:
        tools["bwa"] = lambda reads: (
            [args.bwa, "mem", "-t", args.threads, args.bwa_index] + reads
        )

    se_reads = [reads1]
    # The joined file is written under a .part name and renamed when complete, so an interrupted copy is never taken for a finished one.
    if args.se_both:
        se_reads = [os.path.join(args.out, f"{tag}.r12.fq")]
        if not os.path.exists(se_reads[0]):
            with open(se_reads[0] + ".part", "w") as out_file:
                for read_path in (reads1, reads2):
                    with open(read_path) as handle:
                        shutil.copyfileobj(handle, out_file)
            os.replace(se_reads[0] + ".part", se_reads[0])
    print(
        f"wgsim {args.pairs} pairs x {args.len} bp, error {args.error}, mutations {args.mut} "
        f"({args.indel_frac:.0%} indels), insert {args.insert}+-{args.insert_sd}; scored by wgsim_eval.pl "
        f"alneval -g {args.window}; wrong / implied per MAPQ decile"
    )
    # The reference, the reads and bwa's index are pre-read before every run (measure.py), so every mapper starts with warm inputs.
    idx = (
        [args.bwa_index + ext for ext in (".bwt", ".sa", ".pac", ".ann", ".amb")]
        if args.bwa_index
        else []
    )
    rows = []
    for tool, cmd in tools.items():
        for mode, reads in (("pe", [reads1, reads2]), ("se", se_reads)):
            sam = os.path.join(args.out, f"{tag}.{tool}.{mode}.sam")
            usage = ""
            row = {"mapper": tool, "mode": mode}
            if args.fresh or not os.path.exists(sam) or os.path.getsize(sam) == 0:
                # Every mapper runs under the same memcap resident-memory cap, on any Linux, and measure.py discards the --warmup runs and measures --runs.
                # The SAM scored is the last run's.
                measurement = measure.measure(
                    [str(part) for part in cmd(reads)],
                    args.warmup,
                    args.runs,
                    memcap.parse_size(args.mem),
                    [args.ref, *reads, *idx],
                    scratch=os.path.join(args.out, ".measure"),
                    tag=f"{tool}.{mode}",
                )
                Path(sam + ".log").write_text(measurement.err)
                if not measurement.succeeded:
                    # One mapper failing is a result to report, not a reason to stop the others.
                    why = "CAPPED" if measurement.capped else "FAILED"
                    print(f"  {tool:9s}{mode}  {why} exit {measurement.code}: {measurement.err[-DETAIL_CHARS:]!r}")
                    rows.append({**row, "status": f"{why} exit {measurement.code}"})
                    continue
                shutil.copyfile(measurement.out, sam)
                wall_sd = f"{measurement.wall['sd']:.3g}" if measurement.wall.get("sd") is not None else "-"
                usage = (
                    f" wall {measurement.wall['mean']:.2f} s rss {measurement.rss():.0f} MB sd {wall_sd} runs {measurement.runs} "
                    f"warmup {measurement.warmup}"
                )
                row.update(wall_s=measurement.wall, rss_MB=round(measurement.rss()))
            bins, (hi_n, hi_w) = score(sam, args.samtools, args.evaluator, args.window)
            cells = "  ".join(
                f"{decile * MAPQ_DECILE}-{decile * MAPQ_DECILE + MAPQ_DECILE - 1}:{counts[BIN_WRONG]}/{counts[BIN_IMPLIED]:.1f}"
                for decile, counts in sorted(bins.items())
                if counts[BIN_PRIMARIES]
            )
            print(
                f"  {tool:9s}{mode}  MAPQ>=30 {hi_n:>8} wrong {hi_w:>4} | {cells} |{usage}",
                flush=True,
            )
            rows.append(
                {
                    **row,
                    "mapq30": hi_n,
                    "mapq30_wrong": hi_w,
                    "wrong/implied by MAPQ decile": cells,
                }
            )
    results_file.write_section(
        args.report,
        f"wgsim_{tag}",
        f"Short-read accuracy on wgsim reads ({args.pairs} pairs x {args.len} bp, error {args.error}, "
        f"mutations {args.mut}, insert {args.insert}±{args.insert_sd}), scored by wgsim_eval.pl alneval "
        f"-g {args.window} (wall s mean ± sd, peak RSS MB; SAMs reused unless --fresh show no time)",
        [("", rows)],
        optik=args.optik,
        notes=[f"SAMs and logs: `{os.path.relpath(args.out, ROOT)}`"],
    )


if __name__ == "__main__":
    main()
