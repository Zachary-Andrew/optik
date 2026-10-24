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


def score(sam, samtools, evaluator, window=5):
    """{decile: [primaries, wrong, implied]} and (MAPQ>=30 reads, wrong) from alneval's output."""
    bins = {}
    view = subprocess.Popen(
        [samtools, "view", "-F", "0x904", sam], stdout=subprocess.PIPE, text=True
    )
    for line in view.stdout:
        q = int(line.split("\t", 5)[4])
        b = bins.setdefault(q // 10, [0, 0, 0.0])
        b[0] += 1
        b[2] += 10 ** (-q / 10)
    view.wait()
    view = subprocess.Popen(
        [samtools, "view", "-F", "0x900", sam], stdout=subprocess.PIPE
    )
    ev = subprocess.run(
        ["perl", evaluator, "alneval", "-g", str(window), "-"],
        stdin=view.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    view.wait()
    for row in ev.stdout.splitlines():
        t = row.split()  # "06x  wrong / total  cumulative  ratio"
        if len(t) >= 4 and t[0].endswith("x"):
            bins.setdefault(int(t[0][:-1]), [0, 0, 0.0])[1] = int(t[1])
    hi = [b for d, b in bins.items() if d >= 3]
    return bins, (sum(b[0] for b in hi), sum(b[1] for b in hi))


def main():
    """Simulate reads with wgsim, map them with each mapper and score them with alneval."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pairs", type=int, default=500000)
    ap.add_argument("--len", type=int, default=150)
    ap.add_argument("--error", type=float, default=0.01)
    ap.add_argument(
        "--mut", type=float, default=0.001, help="wgsim -r: sample mutation rate"
    )
    ap.add_argument(
        "--indel-frac",
        type=float,
        default=0.15,
        help="wgsim -R: share of mutations that are indels",
    )
    ap.add_argument("--insert", type=int, default=400, help="wgsim -d: outer distance")
    ap.add_argument("--insert-sd", type=int, default=50, help="wgsim -s")
    ap.add_argument(
        "--window",
        type=int,
        default=5,
        help="wgsim_eval.pl alneval -g: correct within INT bp (its default 5; BWA-MEM's paper 20)",
    )
    ap.add_argument(
        "--se-both",
        action="store_true",
        help="single-end runs map R1 and R2 together as 2N reads (BWA-MEM Fig. 1, top)",
    )
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument(
        "--seed-kw",
        default="21,11",
        help="optik's minimizer k,w: minimap2 -ax sr's (v2.31 options.c), so both "
        "index the same seeds",
    )
    ap.add_argument(
        "--threads",
        default=str(hw.cpus()),
        help="threads of every mapper (default: the usable CPUs)",
    )
    ap.add_argument(
        "--mem",
        default=hw.default_cap_str(),
        help="resident-memory cap of each mapper ('' = uncapped; default: usable memory less 2 GiB)",
    )
    ap.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    ap.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    ap.add_argument(
        "--bwa", default=shutil.which("bwa") or os.path.join(ROOT, "temp/bwa/bwa")
    )
    ap.add_argument("--bwa-index", help="BWA index prefix (adds bwa mem)")
    ap.add_argument("--samtools", default=shutil.which("samtools"))
    ap.add_argument("--wgsim", default=shutil.which("wgsim"))
    ap.add_argument("--evaluator", default=shutil.which("wgsim_eval.pl"))
    ap.add_argument(
        "--fresh",
        action="store_true",
        help="map again even where a SAM exists (the simulated reads are kept)",
    )
    ap.add_argument(
        "--warmup", type=int, default=1, help="warmup runs of every mapper, discarded"
    )
    ap.add_argument(
        "--runs",
        type=int,
        default=3,
        help="measured runs of every mapper: the mean wall time and the peak RSS over them",
    )
    results_file.add_argument(ap)
    a = ap.parse_args()
    for name in ("samtools", "wgsim", "evaluator"):
        if not getattr(a, name):
            sys.exit(
                f"error: {name} not found (samtools ships wgsim and wgsim_eval.pl)"
            )
    os.makedirs(a.out, exist_ok=True)

    tag = f"l{a.len}e{a.error}"
    if (a.mut, a.indel_frac, a.insert, a.insert_sd) != (0.001, 0.15, 400, 50):
        tag += f"r{a.mut}R{a.indel_frac}d{a.insert}s{a.insert_sd}"
    r1, r2 = (os.path.join(a.out, f"{tag}.r{i}.fq") for i in (1, 2))
    if not os.path.exists(r2):
        with open(os.devnull, "w") as dn:
            subprocess.run(
                [
                    a.wgsim,
                    "-S",
                    str(a.seed),
                    "-N",
                    str(a.pairs),
                    "-1",
                    str(a.len),
                    "-2",
                    str(a.len),
                    "-d",
                    str(a.insert),
                    "-s",
                    str(a.insert_sd),
                    "-e",
                    str(a.error),
                    "-r",
                    str(a.mut),
                    "-R",
                    str(a.indel_frac),
                    a.ref,
                    r1,
                    r2,
                ],
                stdout=dn,
                stderr=dn,
                check=True,
            )

    sk, sw = a.seed_kw.split(",")
    tools = {
        "optik": lambda reads: (
            [a.optik, "map", "-r", a.ref, "-k", sk, "-w", sw, "-t", a.threads]
            + [a for r in reads for a in ("-q", r)]
        ),
        "minimap2": lambda reads: (
            [a.minimap2, "-ax", "sr", "-t", a.threads, a.ref] + reads
        ),
    }
    if a.bwa_index:
        tools["bwa"] = lambda reads: (
            [a.bwa, "mem", "-t", a.threads, a.bwa_index] + reads
        )

    se = [r1]
    if a.se_both:
        se = [os.path.join(a.out, f"{tag}.r12.fq")]
        if not os.path.exists(se[0]):
            with open(se[0] + ".part", "w") as o:
                for r in (r1, r2):
                    with open(r) as f:
                        shutil.copyfileobj(f, o)
            os.replace(se[0] + ".part", se[0])
    print(
        f"wgsim {a.pairs} pairs x {a.len} bp, error {a.error}, mutations {a.mut} "
        f"({a.indel_frac:.0%} indels), insert {a.insert}+-{a.insert_sd}; scored by wgsim_eval.pl "
        f"alneval -g {a.window}; wrong / implied per MAPQ decile"
    )
    # The reference, the reads and bwa's index are pre-read before every run (measure.py), so every mapper starts with warm inputs.
    idx = (
        [a.bwa_index + e for e in (".bwt", ".sa", ".pac", ".ann", ".amb")]
        if a.bwa_index
        else []
    )
    rows = []
    for t, cmd in tools.items():
        for mode, reads in (("pe", [r1, r2]), ("se", se)):
            sam = os.path.join(a.out, f"{tag}.{t}.{mode}.sam")
            usage = ""
            row = {"mapper": t, "mode": mode}
            if a.fresh or not os.path.exists(sam) or os.path.getsize(sam) == 0:
                # Every mapper runs under the same memcap resident-memory cap, on any Linux, and measure.py discards the --warmup runs and measures --runs.
                # The SAM scored is the last run's.
                r = measure.measure(
                    [str(c) for c in cmd(reads)],
                    a.warmup,
                    a.runs,
                    memcap.parse_size(a.mem),
                    [a.ref, *reads, *idx],
                    scratch=os.path.join(a.out, ".measure"),
                    tag=f"{t}.{mode}",
                )
                Path(sam + ".log").write_text(r.err)
                if not r.ok:
                    # One mapper failing is a result to report, not a reason to stop the others.
                    why = "CAPPED" if r.capped else "FAILED"
                    print(f"  {t:9s}{mode}  {why} exit {r.code}: {r.err[-200:]!r}")
                    rows.append({**row, "status": f"{why} exit {r.code}"})
                    continue
                shutil.copyfile(r.out, sam)
                sd = f"{r.wall['sd']:.3g}" if r.wall.get("sd") is not None else "-"
                usage = (
                    f" wall {r.wall['mean']:.2f} s rss {r.rss():.0f} MB sd {sd} runs {r.runs} "
                    f"warmup {r.warmup}"
                )
                row.update(wall_s=r.wall, rss_MB=round(r.rss()))
            bins, (hi_n, hi_w) = score(sam, a.samtools, a.evaluator, a.window)
            cells = "  ".join(
                f"{d * 10}-{d * 10 + 9}:{b[1]}/{b[2]:.1f}"
                for d, b in sorted(bins.items())
                if b[0]
            )
            print(
                f"  {t:9s}{mode}  MAPQ>=30 {hi_n:>8} wrong {hi_w:>4} | {cells} |{usage}",
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
        a.report,
        f"wgsim_{tag}",
        f"Short-read accuracy on wgsim reads ({a.pairs} pairs x {a.len} bp, error {a.error}, "
        f"mutations {a.mut}, insert {a.insert}±{a.insert_sd}), scored by wgsim_eval.pl alneval "
        f"-g {a.window} (wall s mean ± sd, peak RSS MB; SAMs reused unless --fresh show no time)",
        [("", rows)],
        optik=a.optik,
        notes=[f"SAMs and logs: `{os.path.relpath(a.out, ROOT)}`"],
    )


if __name__ == "__main__":
    main()
