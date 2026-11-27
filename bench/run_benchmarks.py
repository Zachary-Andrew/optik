#!/usr/bin/env python3
"""run_benchmarks.py - one entry point for optik's validations and benchmarks.

  python3 bench/run_benchmarks.py list
      lists every bench, what it checks or measures, what it needs, and the section it writes in the results file
  python3 bench/run_benchmarks.py run --only small_sets,scaling [--warmup 1 --runs 3]
  python3 bench/run_benchmarks.py run --all
      runs the named benches one after another, each with its usual defaults
      and collects all of their sections in one results file (bench/results/RESULTS.md, or --report PATH)
  python3 bench/run_benchmarks.py run --only wgsim -- --pairs 100000
      passes everything after -- to each selected bench

A bench whose inputs or tools are missing is skipped with the command that provides them.
The 'runs' section of the results file lists every bench this script ran, its command, exit status and wall time, so one file holds both the metrics and how they were made.
`python3 bench/report.py` lists the sections, and `python3 bench/report.py --json` prints them as data.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import genomes
import hw
import report as results_file
from units import SECONDS_PER_MINUTE

# Fetched into temp/genomes/ by bench/genomes.py before a bench that reads it runs
ECOLI = os.path.join(genomes.DEST, genomes.GENOMES["ecoli"][0])
OPTIK = os.path.join(ROOT, "target", "release", "optik")
WORK = os.path.join(ROOT, "temp", "run_benchmarks")
MINIMAP2 = os.path.join(ROOT, "temp", "minimap2", "minimap2")
SCALING = os.path.join(ROOT, "temp", "scaling")
REAL = os.path.join(ROOT, "temp", "real")


def exe(name):
    """True when `name` is on PATH."""
    return shutil.which(name) is not None


def path(target):
    """True when the file or directory `p` exists."""
    return os.path.exists(target)


# Each bench names its script and arguments, the shared flags it takes, what it needs as (check, what, how to get it), the results-file sections it writes, and the wall time of its last full run on a 16-CPU, 31 GiB machine (temp/audit/chain.log and temp/run_benchmarks/validations.log, 2026-10-23), where one was timed.
BENCHES = {
    "validations": {
        "cmd": ["bash", "bench/run_validations.sh"],
        "flags": (),
        "what": "Correctness suite: exact equality with SSHash, Cuttlefish, Bifrost, KMC and "
        "minimap2 oracles, kernel stress tests, SAM conventions and error paths",
        "needs": [(lambda: path(OPTIK), "optik", "cargo build --release")],
        "sections": "runs (pass or fail; details in its log)",
        "time": "28 min",
    },
    "small_sets": {
        "cmd": ["python3", "bench/compare_tools_small_sets.py", "--mode", "both"],
        "flags": ("--warmup", "--runs", "--threads"),
        "what": "optik against SSHash, SBWT, FMSI, Cuttlefish, Bifrost, minimap2 and wfmash on "
        "small genomes, every output verified",
        "needs": [
            (lambda: path(OPTIK), "optik", "cargo build --release"),
            (lambda: exe("hyperfine"), "hyperfine", "see README, Install"),
            (lambda: path(MINIMAP2), "the compared tools", "bench/build_tools.sh"),
        ],
        "sections": "small_sets",
        "time": "75 min",
    },
    "published_smoke": {
        "cmd": ["python3", "bench/compare_tools_published_sets.py", "smoke"],
        "flags": ("--threads",),
        "what": "Every published dataset on a prefix of its input: each tool builds and "
        "queries, and the query tools agree",
        "needs": [(lambda: path(OPTIK), "optik", "cargo build --release")],
        "sections": "published_smoke",
        "time": "6 min",
    },
    "published": {
        "cmd": [
            "python3",
            "bench/compare_tools_published_sets.py",
            "bench",
            "--mode",
            "both",
        ],
        "flags": ("--warmup", "--runs", "--threads"),
        "what": "The papers' own datasets and protocols for SSHash, SBWT, FMSI, Cuttlefish, "
        "GGCAT, minimap2, BWA-MEM, Winnowmap and BiWFA",
        "needs": [
            (
                lambda: path(os.path.join(ROOT, "temp", "published_sets")),
                "the published datasets",
                "python3 bench/compare_tools_published_sets.py fetch",
            )
        ],
        "sections": "published_query, published_dbg, published_map, published_align",
        "time": "about a day",
    },
    "scaling": {
        "cmd": ["python3", "bench/scaling_bench.py", os.path.join(WORK, "scaling")],
        "flags": ("--warmup", "--runs", "--threads"),
        "what": "Index, graph and mapping cost of every tool on E. coli, chr21 and chr1",
        "needs": [
            (
                lambda: path(os.path.join(SCALING, "chr21.fa")),
                "chr21 and chr1",
                "bench/build_tools.sh scaling",
            )
        ],
        "sections": "scaling",
        "time": "6.5 min (E. coli and chr21)",
    },
    "dbg": {
        "cmd": ["python3", "bench/dbg_bench.py", os.path.join(WORK, "dbg")],
        "flags": ("--warmup", "--runs", "--threads"),
        "what": "optik dbruijn against Cuttlefish, GGCAT and Bifrost on read sets, GRCh38 and "
        "a gut metagenome",
        "needs": [(lambda: path(OPTIK), "optik", "cargo build --release")],
        "sections": "dbg",
        "time": "",
    },
    "wgsim": {
        "cmd": [
            "python3",
            "bench/wgsim_check.py",
            "--ref",
            ECOLI,
            "--out",
            os.path.join(WORK, "wgsim"),
        ],
        "flags": ("--warmup", "--runs", "--threads"),
        "what": "Short-read accuracy per MAPQ decile on wgsim reads, scored by wgsim_eval.pl",
        "needs": [
            (
                lambda: exe("wgsim") and exe("wgsim_eval.pl"),
                "wgsim",
                "install samtools",
            ),
            (lambda: path(MINIMAP2), "minimap2", "bench/build_tools.sh minimap2"),
        ],
        "sections": "wgsim_<read settings>",
        "time": "",
    },
    "paired": {
        "cmd": [
            "python3",
            "bench/paired_bench.py",
            "--ref",
            ECOLI,
            "--out",
            os.path.join(WORK, "paired"),
            "--gate",
        ],
        "flags": ("--threads",),
        "what": "Paired-end mapping: placement, proper-pair flags, TLEN, mate fields and "
        "samtools/Picard downstream checks",
        "needs": [(lambda: exe("samtools"), "samtools", "see README, Install")],
        "sections": "paired",
        "time": "",
    },
    "mapq": {
        "cmd": [
            "python3",
            "bench/mapq_calibration.py",
            "--multi",
            "--gate",
            "--out",
            os.path.join(WORK, "mapq"),
        ],
        "flags": ("--threads",),
        "what": "MAPQ calibration of optik on fresh simulated short and long reads of a "
        "repeat-rich genome",
        "needs": [(lambda: path(OPTIK), "optik", "cargo build --release")],
        "sections": "mapq_multi",
        "time": "",
    },
    "real_anchored": {
        "cmd": [
            "python3",
            "bench/real_reads.py",
            "anchored",
            "--out",
            os.path.join(WORK, "real_anchored"),
        ],
        "flags": ("--threads",),
        "what": "MAPQ of real HG002 read mates anchored by a uniquely placed partner",
        "needs": [
            (
                lambda: False,
                "--ref, --r1 and --r2 after --",
                "see bench/real_reads.py --help",
            )
        ],
        "sections": "real_anchored",
        "time": "4 min (HG002 chr20 slice)",
    },
    "real_variants": {
        "cmd": [
            "python3",
            "bench/real_reads.py",
            "variants",
            "--out",
            os.path.join(WORK, "real_variants"),
        ],
        "flags": ("--threads",),
        "what": "Germline variants called from each mapper's SAM against the GIAB benchmark",
        "needs": [
            (
                lambda: False,
                "--ref, --r1, --r2, --truth-vcf and --truth-bed after --",
                "see bench/real_reads.py --help",
            )
        ],
        "sections": "real_variants",
        "time": "2.5 min (HG002 chr20 slice)",
    },
    "micro": {
        "cmd": ["python3", "bench/micro_bench.py"],
        "flags": (),
        "what": "Hardware counters per operation of optik's hot kernels, with calibration",
        "needs": [
            (lambda: exe("g++") or exe("c++"), "a C++ compiler", "see README, Install")
        ],
        "sections": "micro",
        "time": "",
    },
}


def missing(bench, extra):
    """The needs of bench `b` that are not met, as 'what (how)' strings.

    A need whose check cannot be met by files on disk, such as real reads given after --, counts as met when extra arguments were given."""
    out = []
    for check, what, how in bench["needs"]:
        if check():
            continue
        if extra and what.endswith("after --"):
            continue
        out.append(f"{what} ({how})")
    return out


def cmd_list():
    """Prints every bench with what it does, what it needs and the section it writes."""
    for name, bench in BENCHES.items():
        print(f"{name}")
        print(f"  {bench['what']}")
        print(
            f"  command  {' '.join(os.path.relpath(arg, ROOT) if arg.startswith(ROOT) else arg for arg in bench['cmd'])}"
        )
        print(f"  writes   {bench['sections']}")
        if bench["time"]:
            print(f"  time     {bench['time']} on 16 CPUs")
        gaps = missing(bench, [])
        print("  ready" if not gaps else "  needs    " + "; ".join(gaps))
        print()


def cmd_run(args, extra):
    """Runs the selected benches in order and writes the 'runs' section."""
    names = list(BENCHES) if args.all else [name for name in args.only.split(",") if name]
    bad = [name for name in names if name not in BENCHES]
    if bad or not names:
        sys.exit(
            f"unknown bench {bad}; `run_benchmarks.py list` shows them"
            if bad
            else "name benches with --only or use --all"
        )
    os.makedirs(WORK, exist_ok=True)
    # Benches of earlier invocations keep their rows, and a bench run again replaces its own.
    earlier = []
    if args.report and os.path.exists(args.report):
        old = results_file.read_report(args.report).get("runs", {}).get("tables", {})
        earlier = [row for row in old.get("", []) if row.get("bench") not in names]
    rows = []

    def record():
        """Writes the runs section, after every bench so an interrupted campaign keeps what ran."""
        results_file.write_section(
            args.report,
            "runs",
            "The benches run by bench/run_benchmarks.py, with their exit status and wall time",
            [("", earlier + rows)],
        )

    for name in names:
        bench = BENCHES[name]
        gaps = missing(bench, extra)
        if gaps:
            print(f"[skip] {name}: needs {'; '.join(gaps)}", flush=True)
            rows.append(
                {
                    "bench": name,
                    "status": "SKIP",
                    "wall_min": None,
                    "command": "",
                    "log": "needs " + "; ".join(gaps),
                }
            )
            record()
            continue
        cmd = list(bench["cmd"])
        for flag in bench["flags"]:
            value = getattr(args, flag.removeprefix("--"))
            if value is not None:
                cmd += [flag, str(value)]
        if name != "validations":
            cmd += ["--report", args.report]
        cmd += extra
        if ECOLI in cmd:
            genomes.path("ecoli")
        log = os.path.join(WORK, f"{name}.log")
        print(
            f"[run] {name}: {' '.join(cmd)}\n      log {os.path.relpath(log, ROOT)}",
            flush=True,
        )
        start = time.time()
        with open(log, "w") as log_file:
            exit_code = subprocess.run(
                cmd, cwd=ROOT, stdout=log_file, stderr=subprocess.STDOUT, check=False
            ).returncode
        mins = (time.time() - start) / SECONDS_PER_MINUTE
        print(
            f"[{'ok' if exit_code == 0 else 'FAIL'}] {name}: exit {exit_code} after {mins:.1f} min",
            flush=True,
        )
        rows.append(
            {
                "bench": name,
                "status": "PASS" if exit_code == 0 else f"FAIL exit {exit_code}",
                "wall_min": round(mins, 1),
                "command": " ".join(
                    os.path.relpath(arg, ROOT) if arg.startswith(ROOT) else arg for arg in cmd
                ),
                "log": os.path.relpath(log, ROOT),
            }
        )
        record()
    count = {
        status: sum(row["status"].startswith(status) for row in rows)
        for status in ("PASS", "FAIL", "SKIP")
    }
    print(
        f"{count['PASS']} passed, {count['FAIL']} failed, {count['SKIP']} skipped -> {args.report}"
    )
    return 1 if count["FAIL"] else 0


def main():
    """Lists or runs the benches."""
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        sep_idx = argv.index("--")
        argv, extra = argv[:sep_idx], argv[sep_idx + 1 :]
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="show every bench, what it needs and what it writes")
    run_parser = sub.add_parser("run", help="run benches into one results file")
    group = run_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--only", help="comma list of bench names (see list)")
    group.add_argument(
        "--all", action="store_true", help="every bench, in the order of list"
    )
    run_parser.add_argument(
        "--warmup",
        type=int,
        help="warmup runs of every timed command (each bench's default when omitted)",
    )
    run_parser.add_argument(
        "--runs",
        type=int,
        help="measured runs of every timed command (each bench's default when omitted)",
    )
    run_parser.add_argument(
        "--threads",
        type=int,
        help=f"threads of every tool (default: the usable CPUs, {hw.cpus()} here)",
    )
    results_file.add_argument(run_parser)
    args = parser.parse_args(argv)
    if args.cmd == "list":
        cmd_list()
        return 0
    return cmd_run(args, extra)


if __name__ == "__main__":
    sys.exit(main())
