#!/usr/bin/env python3
"""micro_bench.py - per-operation hardware counters of optik's hot kernels (bench/cpp/micro_kernels.cpp, with the harness in bench/cpp/nanobench.hpp), without root.
It applies nanoBench's method (the difference of n and 2n operations, a warm-up and the median of repetitions) to unprivileged perf_event counters.

The calibration kernels are checked before any kernel result is reported.
  nop / add_dependent       1 instruction per op (plus 2/64 for the loop), and add at 1 cycle
  add_independent           4 instructions per op at about 1 cycle
  chase_l1                  the L1 load-to-use latency (3-8 cycles) with no L1D misses
  chase_dram                over 100 cycles with about 1 L1D miss per op
If the harness cannot open counters (perf_event_paranoid > 2, or a VM without a PMU) or a calibration fails, the run is reported as unusable instead of printing numbers that mean nothing.

Usage:
  python3 bench/micro_bench.py [--quick] [--json OUT] [--warmup 5] [--runs 31]
  (tool_bench.py --micro runs it and stores the result with the timings)
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import report as results_file
from diagnostics import LOG_TAIL_CHARS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# discarded warm-ups and measured repetitions of every kernel by default
DEFAULT_WARMUP = 5
DEFAULT_RUNS = 31
# lines of a failed build's or run's stderr kept
STDERR_TAIL_LINES = 6

CALIBRATION = {
    # kernel -> [(metric, lo, hi), ...]
    "nop": [("instructions", 1.0, 1.06)],
    "add_dependent": [("instructions", 1.0, 1.06), ("cycles", 0.95, 1.10)],
    "add_independent": [("instructions", 4.0, 4.06), ("cycles", 0.95, 1.50)],
    "chase_l1": [
        ("instructions", 1.0, 1.06),
        ("cycles", 3.0, 8.0),
        ("l1d_misses", -0.05, 0.05),
    ],
    "chase_dram": [("cycles", 100.0, 5000.0), ("l1d_misses", 0.8, 1.2)],
}


def build(out):
    """Compile the micro-kernel harness to `out`."""
    cxx = os.environ.get("CXX", "g++")
    cmd = [
        cxx,
        "-std=c++17",
        "-O3",
        "-march=native",
        "-funroll-loops",
        "-fopenmp",
        "-pthread",
        "-I",
        str(ROOT / "include"),
        "-I",
        str(ROOT / "src"),
        "-I",
        str(HERE / "cpp"),
        str(HERE / "cpp" / "micro_kernels.cpp"),
        "-o",
        str(out),
    ]
    # Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0.
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        env=dict(os.environ, NIX_ENFORCE_NO_NATIVE="0"),
    )
    if proc.returncode != 0:
        raise RuntimeError("micro_kernels build failed:\n" + proc.stderr[-LOG_TAIL_CHARS:])


def run(quick=False, workdir=None, warmup=DEFAULT_WARMUP, runs=DEFAULT_RUNS):
    work_dir = Path(workdir or "/tmp/optik_micro")
    work_dir.mkdir(parents=True, exist_ok=True)
    exe = work_dir / "micro_kernels"
    build(exe)
    # For every kernel `warmup` measurements are discarded, then the median and quartiles of `runs` are kept.
    env = dict(os.environ, NB_WARMUP=str(warmup), NB_RUNS=str(runs))
    proc = subprocess.run(
        [str(exe)] + (["--quick"] if quick else []),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    rows = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
    return evaluate(rows, proc.returncode, proc.stderr, quick)


def evaluate(rows, returncode=0, stderr="", quick=False):
    """Usable only if the counters opened and every calibration kernel lies inside its expected range.
    The function is pure and unit-tested in test_harness.py."""
    res = {
        "kernels": rows,
        "stderr": stderr.strip().splitlines()[-STDERR_TAIL_LINES:],
        "calibration": {},
        "usable": True,
    }
    by_kernel = {row["kernel"]: row for row in rows}
    if returncode != 0 or not rows:
        res["usable"] = False
        res["reason"] = f"micro_kernels exited {returncode}"
        return res
    # Hardware counters may be locked down, and the reason then names perf_event_paranoid so the fix is clear.
    errs = [row for row in rows if "error" in row]
    if errs:
        res["usable"] = False
        res["reason"] = (
            f"counters unavailable: {errs[0]['error']} (perf_event_paranoid "
        )
        try:
            res["reason"] += (
                Path("/proc/sys/kernel/perf_event_paranoid").read_text().strip() + ")"
            )
        except OSError:
            res["reason"] += "?)"
        return res
    # Each calibration kernel has known counter ranges, and a value outside them means the counters cannot be trusted on this machine.
    for kernel, rules in CALIBRATION.items():
        if kernel not in by_kernel:
            if kernel == "chase_dram" and quick:
                continue
            res["calibration"][kernel] = {"pass": False, "missing": True}
            res["usable"] = False
            continue
        bad = [
            (metric, by_kernel[kernel].get(metric), low, high)
            for metric, low, high in rules
            if not (low <= by_kernel[kernel].get(metric, float("nan")) <= high)
        ]
        res["calibration"][kernel] = {"pass": not bad, "out_of_range": bad}
        if bad:
            res["usable"] = False
    if not res["usable"]:
        res["reason"] = "calibration failed: " + ", ".join(
            kernel for kernel, verdict in res["calibration"].items() if not verdict["pass"]
        )
    return res


def table(res):
    """Text table of the kernel results, or the reason the run is unusable."""
    lines = []
    if not res["usable"]:
        return f"  micro-benchmarks unusable: {res.get('reason')}"
    lines.append(
        f"  {'kernel':24} {'cycles':>9} {'instr':>9} {'IPC':>5} {'L1D':>7} {'LLC':>7} {'dTLB':>6} {'br-miss':>7}"
    )
    for row in res["kernels"]:
        cycles, instrs = row["cycles"], row["instructions"]
        lines.append(
            f"  {row['kernel']:24} {cycles:9.1f} {instrs:9.1f} {instrs / cycles if cycles > 0 else 0:5.2f} {row['l1d_misses']:7.2f} "
            f"{row['llc_misses']:7.2f} {row['dtlb_misses']:6.2f} {row['branch_misses']:7.2f}"
        )
    lines.append(
        f"  (per operation; median of {res.get('runs', DEFAULT_RUNS)} repetitions of n vs 2n after "
        f"{res.get('warmup', DEFAULT_WARMUP)} discarded warm-ups; calibration: "
        + ", ".join(
            f"{kernel} {'ok' if verdict['pass'] else 'FAIL'}"
            for kernel, verdict in res["calibration"].items()
        )
        + ")"
    )
    return "\n".join(lines)


def main():
    """Build and run the micro benchmarks and print or save the results."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--quick", action="store_true", help="cache-resident sizes only, no DRAM chase"
    )
    parser.add_argument("--json", help="write the result here")
    parser.add_argument("--workdir", default="/tmp/optik_micro")
    parser.add_argument(
        "--warmup",
        type=int,
        default=DEFAULT_WARMUP,
        help="warm-up measurements per kernel, discarded",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=DEFAULT_RUNS,
        help="kept measurements per kernel (median and quartiles)",
    )
    results_file.add_argument(parser)
    args = parser.parse_args()
    res = run(args.quick, args.workdir, args.warmup, args.runs)
    res["warmup"], res["runs"] = args.warmup, args.runs
    print(table(res))
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1))
    # Only trustworthy counters reach the results file, while an unusable run records why instead.
    if res["usable"]:
        rows = [
            {
                "kernel": row["kernel"],
                "cycles": row["cycles"],
                "instructions": row["instructions"],
                "IPC": row["instructions"] / row["cycles"] if row["cycles"] > 0 else 0.0,
                "L1D_misses": row["l1d_misses"],
                "LLC_misses": row["llc_misses"],
                "dTLB_misses": row["dtlb_misses"],
                "branch_misses": row["branch_misses"],
            }
            for row in res["kernels"]
        ]
        calib = [
            {"kernel": kernel, "gate": verdict["pass"]} for kernel, verdict in res["calibration"].items()
        ]
        notes = []
    else:
        rows, calib, notes = [], [], [f"Unusable on this machine: {res.get('reason')}"]
    results_file.write_section(
        args.report,
        "micro",
        f"Hardware counters per operation of optik's hot kernels (median of {args.runs} repetitions "
        f"of n against 2n after {args.warmup} warm-ups)",
        [("kernels", rows), ("calibration", calib)],
        notes=notes,
    )
    sys.exit(0 if res["usable"] else 1)


if __name__ == "__main__":
    main()
