#!/usr/bin/env python3
"""compare_tools_small_sets.py - optik against every published tool it implements, on small inputs.

The inputs are SSHash's E. coli, 100-Salmonella and Penicillium unitig sets for the k-mer indexes, the Salmonella LT2 genome and 30x simulated reads for the graphs, simulated short, ONT-like, HiFi-like, chimeric and random reads for mapping, random sequence pairs for alignment, and a repeat-rich synthetic genome for MAPQ calibration. bench/compare_tools_published_sets.py runs the same comparisons on the datasets of the tools' own papers.

Two modes, run alone or together

  parity   identical inputs and matched parameters for every tool in a section, i.e. k = 31, minimizer m = 15 where a tool takes one, and the same alignment penalties and thread counts, with whole-process wall time and single-threaded queries.
           This table compares algorithms.
  optimal  every tool at its own documented best settings, meaning its presets, its recommended minimizer length, variant or query mode, and all physical cores where it parallelises.
           This is the table a user would see.

Sections
  kmer   exact k-mer membership and dictionary lookups with optik build-spss, optik build-fm, SSHash, SBWT and FMSI with KmerCamel on SSHash's unitig sets
  dbg    compacted de Bruijn graphs from optik dbruijn, Cuttlefish and Bifrost on a genome and on 30x simulated reads
  map    read mapping with optik map, minimap2 and, for long reads, wfmash on short, ONT-like, HiFi-like, chimeric and random reads with truth
  align  exact gap-affine pairwise alignment, optik WFA and BiWFA against WFA2-lib's WFA and BiWFA
  mapq   MAPQ calibration of optik against minimap2 by bench/mapq_calibration.py

Each row carries a correctness gate.
A row whose output is wrong is marked FAIL in the tables and is never silently timed.
The gates are
  kmer   num_kmers and num_positive equal an independent brute-force truth
  dbg    the unitig k-mer set equals KMC's solid k-mers at the same cutoff, and no k-mer appears in two unitigs
  map    every SAM record consumes SEQ exactly, NM recomputed from the reference matches NM:i, each read has one primary, and accuracy is scored against the truth
  align  all exact methods report the same optimal scores, and each CIGAR re-scores to its own score

Timing goes through bench/measure.py, which runs hyperfine on the whole process, drops the --warmup runs and computes statistics over the --runs measured runs, with each run's inputs cached beforehand.
Peak RSS and block I/O of every measured run come from wait4 in rusage_exec.
The correctness gates read the measured runs' own output.
Per-k-mer answers are the exception, since one untimed run writes them and the timed runs send them to /dev/null, as counting tools write none.
Every command line, tool revision, compiler and machine detail is saved in the run directory.

The tool paths and sequence helpers are in small_common.py, and the external correctness oracles built on KMC, Picard and samtools in external_oracles.py.

Usage:
  python3 bench/compare_tools_small_sets.py [--mode parity|optimal|both]
      [--sections kmer,dbg,map,align,mapq] [--runs 3] [--warmup 1]
      [--threads <usable CPUs>] [--out bench/results/compare_tools_small_sets] [--quick]"""

import argparse
import datetime
import json
import os
import shutil
import sys
from pathlib import Path

import hw
import measure
import report as results_file
from small_common import BIN, TEMP_DIR
from small_sets.align import (
    parse_wfa_bench_summary,
    section_align,
)
from small_sets.common import (
    DEFAULT_RUNS,
    HERE,
    ROOT,
    Bench,
    _fmt,
    build_wfab,
    check_optik_fresh,
    environment,
)
from small_sets.dbg import (
    gfa_or_fasta_seqs,
    section_dbg,
)
from small_sets.kmer import fmsi_counts, parse_counts, sbwt_counts, section_kmer
from small_sets.mapping import (
    section_map,
    write_fai,
)
from small_sets.mapq import section_mapq
from units import BYTES_PER_GIB

# Other bench scripts import these through this module, so they stay part of its interface.
__all__ = [
    "Bench",
    "_fmt",
    "build_wfab",
    "check_optik_fresh",
    "environment",
    "fmsi_counts",
    "gfa_or_fasta_seqs",
    "parse_counts",
    "parse_wfa_bench_summary",
    "sbwt_counts",
    "write_fai",
]


def report(bench, out, env, dest=""):
    """Write results.json and the Markdown tables of a run to `out`, and the run's section of the sectioned results file `dest`."""
    md_lines = [
        f"# optik vs published tools, small sets ({env['date']})",
        "",
        (
            f"- CPU: {env['cpu']} ({env['physical_cores']} cores, {env['cores_threads']}), {env['memory']}, "
            f"kernel {env['kernel']}, THP {env['thp']}, governor {env['governor']}"
        ),
        f"- {env['gcc']}; {env['rustc']}; optik {env['optik_git']}",
        f"- Builds: {env['build_flags']}",
        "- Tool revisions: " + ", ".join(f"{tool} {version}" for tool, version in env["tools"].items()),
        (
            f"- Times: hyperfine whole-process mean ± sd over {bench.runs} measured runs after {bench.warmup} discarded "
            "warmup run(s), inputs read into the page cache before each (bench/measure.py); RSS: peak over the "
            "measured runs (wait4). Alignment: the harness's own timer, same warmup/runs."
        ),
        (
            "- Gate: PASS = output verified (see compare_tools_small_sets.py docstring); FAIL rows are wrong "
            "output, not just slow. Each tool is judged under its own documented model (Cuttlefish: cutoff on "
            "(k+1)-mer edges). `sam_spec` (SAM conventions such as wfmash's one primary per segment) and "
            "`calibrated` (MAPQ over-confidence) are measured results reported for every tool, not gates."
        ),
        "",
    ]
    # This table is written by hand from each tool's documentation, so it has to be updated when a tool gains or drops a feature.
    md_lines += [
        "## Use-case capabilities (documented features, not measured)",
        "",
        "| k-mer index | membership + streaming | dense ids (Lookup) | Access(id) | weights | add/remove after build |",
        "|---|---|---|---|---|---|",
        "| optik build-spss | yes | yes (`query --ids`) | yes (library) | no | yes (`spss-update`) |",
        "| SSHash | yes | yes (library / `bench`) | yes | yes (`--weighted`) | no |",
        "| SBWT | yes | colex ranks (`search` output) | yes | no | no |",
        "| FMSI | yes | no | no | no | set operations on masked superstrings |",
        "| optik build-fm | yes | no | no | no | no |",
        "",
        "| dBG | graph model | min-count filter | output |",
        "|---|---|---|---|",
        "| optik dbruijn | edge-centric (observed (k+1)-mers) | exact k-mer counts | GFA1 (S + L) |",
        "| Cuttlefish 2/3 | edge-centric | exact, on (k+1)-mer edges (`--read`) | FASTA unitigs |",
        "| Bifrost | node-centric | singletons only, Bloom filter (approximate) | GFA1 |",
        "",
        "| mapper | output | short reads | long reads |",
        "|---|---|---|---|",
        "| optik map | SAM, primary + supplementary (SA), MAPQ | yes | yes (chaining + WFA/BiWFA) |",
        "| minimap2 | SAM/PAF, primary + supplementary | yes (`-x sr`) | yes (`map-ont`, `map-hifi`) |",
        "| wfmash | SAM/PAF, one primary per segment | not designed for | yes (WFA) |",
        "",
    ]
    rows_by_section = {}
    for row in bench.rows:
        rows_by_section.setdefault((row["section"], row["mode"]), []).append(row)
    for (sec, mode), rows in rows_by_section.items():
        md_lines.append(f"## {sec} -- {mode}")
        md_lines.append("")
        # The columns follow the first row's order, and keys that only later rows report are appended, since tools report different metrics.
        cols = [
            col for col in rows[0] if col not in ("section", "mode", "gate_detail", "cmd")
        ]
        for row in rows[1:]:
            for col in row:
                if col not in cols and col not in ("section", "mode", "gate_detail", "cmd"):
                    cols.append(col)
        md_lines.append("| " + " | ".join(cols) + " |")
        md_lines.append("|" + "---|" * len(cols))
        for row in rows:
            cells = []
            for col in cols:
                cell = row.get(col, "")
                if col == "gate":
                    cell = "PASS" if cell else "FAIL"
                # Timed cells hold hyperfine's statistics and print as mean ± sd.
                elif isinstance(cell, dict) and "mean" in cell:
                    cell = f"{cell['mean']:.4g}" + (
                        f" ± {cell['sd']:.2g}" if cell.get("sd") is not None else ""
                    )
                elif isinstance(cell, float):
                    cell = f"{cell:.4g}"
                cells.append(str(cell))
            md_lines.append("| " + " | ".join(cells) + " |")
        # Gate failures follow their table, so a reader sees why a row failed without opening the logs.
        fails = [row for row in rows if not row.get("gate")]
        if fails:
            md_lines.append("")
            md_lines.append("Gate failures:")
            for row in fails:
                md_lines.append(
                    f"- {row.get('tool') or row.get('method')} / {row.get('dataset') or row.get('length')}: "
                    f"{row.get('gate_detail')}"
                )
        md_lines.append("")
    Path(os.path.join(out, "tables.md")).write_text("\n".join(md_lines) + "\n")
    Path(os.path.join(out, "results.json")).write_text(
        json.dumps(bench.rows, indent=1, default=str)
    )
    # Commands and gate details are dropped from the results-file tables but remain in results.json and tables.md.
    hidden = ("section", "mode", "gate_detail", "cmd")
    tables = [
        (
            f"{sec} {mode}",
            [{key: value for key, value in row.items() if key not in hidden} for row in rows],
        )
        for (sec, mode), rows in rows_by_section.items()
    ]
    fails = [
        f"- FAIL {row['section']} {row['mode']} {row.get('tool') or row.get('method')} / "
        f"{row.get('dataset') or row.get('length')}: {row.get('gate_detail')}"
        for row in bench.rows
        if not row.get("gate")
    ]
    results_file.write_section(
        dest,
        "small_sets",
        "optik against SSHash, SBWT, FMSI, Cuttlefish, Bifrost, minimap2 and wfmash on small sets"
        f" (wall s mean ± sd over {bench.runs} runs after {bench.warmup} warmup, peak RSS MB,"
        " gate PASS when the output is verified)",
        tables,
        notes=fails
        + [f"Raw rows, logs and the full tables: `{os.path.relpath(out, ROOT)}`"],
    )


def main():
    """Check the binaries, then run the sections and modes the command line selects."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--mode", default="both", choices=["parity", "optimal", "both"])
    parser.add_argument("--sections", default="kmer,dbg,map,align,mapq")
    parser.add_argument(
        "--runs",
        type=int,
        default=DEFAULT_RUNS,
        help="measured runs per command (statistics over these)",
    )
    parser.add_argument(
        "--warmup", type=int, default=1, help="warmup runs per command, discarded"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every parallel tool (default: the usable CPUs)",
    )
    parser.add_argument(
        "--cap-gb",
        type=int,
        default=hw.default_cap_bytes() // BYTES_PER_GIB,
        help="per-command resident-memory cap in GiB (memcap; default: usable memory less 2 GiB)",
    )
    parser.add_argument(
        "--out", default=os.path.join(HERE, "results", "compare_tools_small_sets")
    )
    parser.add_argument(
        "--work",
        default=os.path.join(TEMP_DIR, "compare_tools_small_sets"),
        help="inputs and indexes",
    )
    parser.add_argument("--quick", action="store_true", help="small subset (smoke test)")
    results_file.add_argument(parser)
    args = parser.parse_args()
    # Tools run inside per-dataset directories, where relative paths would resolve wrongly.
    args.out, args.work = os.path.abspath(args.out), os.path.abspath(args.work)
    check_optik_fresh()
    stamp = datetime.datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")
    out = os.path.join(args.out, stamp + ("-quick" if args.quick else ""))
    os.makedirs(out, exist_ok=True)
    os.makedirs(args.work, exist_ok=True)
    for binary in (
        "optik",
        "sshash",
        "sbwt",
        "fmsi",
        "kmercamel",
        "cuttlefish",
        "bifrost",
        "minimap2",
        "wfmash",
        "kmc",
        "kmc_dump",
    ):
        if not os.access(BIN[binary], os.X_OK):
            sys.exit(f"missing {binary}: {BIN[binary]} (bench/build_tools.sh)")
    if not shutil.which("hyperfine"):
        sys.exit("hyperfine not on PATH")
    measure.rusage_exec()  # built up front, or rebuilt when its source is newer
    env = environment(out)
    bench = Bench(out, args.runs, args.warmup, args.cap_gb, os.path.join(args.work, "scratch"))
    modes = ["parity", "optimal"] if args.mode == "both" else [args.mode]
    secs = args.sections.split(",")
    # The report is rewritten after every mode, so an interrupted run keeps finished tables.
    for mode in modes:
        if "kmer" in secs:
            section_kmer(bench, mode, args.work, args.threads, args.quick)
        if "dbg" in secs:
            section_dbg(bench, mode, args.work, args.threads, args.quick)
        if "map" in secs:
            section_map(bench, mode, args.work, args.threads, args.quick)
        if "align" in secs:
            section_align(bench, mode, args.work, args.quick)
        report(bench, out, env, args.report)
    if "mapq" in secs:
        section_mapq(bench, "calibration", args.work, args.threads, args.quick)
    report(bench, out, env, args.report)
    nfail = sum(1 for row in bench.rows if not row.get("gate"))
    print(f"\n{len(bench.rows)} rows, {nfail} gate failures -> {out}/tables.md")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
