#!/usr/bin/env python3
"""compare_tools_published_sets.py - optik against the published tools on the datasets their own most recent papers used.

Each dataset runs in the section it came from, query data in query, graph data in dbg and so on, and the way that paper used it.
The same comparisons on small inputs are in bench/compare_tools_small_sets.py.

Every dataset is taken from the paper or from the data or code deposit the paper cites, and its `source` field quotes where.
Nothing is substituted.
A dataset that cannot be obtained or run as the paper describes is listed with the reason and never fetched.
The catalogue itself is in published_catalog.py.
Downloads, input formatting, tool command lines and the mapping evaluators are in published_fetch.py, published_inputs.py, published_tools.py and published_mapping.py.

  python3 bench/compare_tools_published_sets.py list [--md]
  python3 bench/compare_tools_published_sets.py tools          # GGCAT, pbsim2, k8, Winnowmap, mason2 (pinned)
  python3 bench/compare_tools_published_sets.py fetch [--tier auto|S|M|L] [--section query] [--only ID,..]
  python3 bench/compare_tools_published_sets.py smoke [--only ID,..]   # formats on a slice, every tool
  python3 bench/compare_tools_published_sets.py status [--watch 10]    # download progress, rate, ETA
  python3 bench/compare_tools_published_sets.py bench [--tier auto|S|M|L] [--section query]
                                                      [--mode parity|optimal|both] [--plan paper]
                                                      [--warmup 1] [--runs 3] [--threads <usable CPUs>]
                                                      [--k-list 15,19,21,23,27,31] [--t-list 1,2,4,8,16]

Matched sweeps (--k-list, --t-list) run every tool at the same k and the same thread count, so k = 15 on 1 thread is compared with k = 15 on 1 thread.
A dataset is rebuilt at a listed k only where its input allows it, which rules out the SSHash 2026 deposits that exist only at k = 31.
The sweep compares tools on sweep-chr1, sweep-pbsim2-chm13 and sweep-mason2-hg38.
Steps without a thread option, such as FMSI's build and every compared tool's queries, are measured once at the first thread count, while optik's queries are swept as its own scaling.
Mappers get the same seed k and w = 10, within minimap2's limit of k <= 28, while bwa mem has no comparable k.
On read sets the graph tools also run at cutoff 2, the only cutoff Bifrost takes.

What runs depends on the machine rather than on a fixed list.

  disk     Before fetching, a dataset needs room for its download plus the estimated decompressed input, indexes and temporary files, as disk_need computes them.
  memory   A tool is SKIPped on a dataset when its memory need exceeds what this machine allows, i.e. MemTotal or a cgroup limit the bench runs under.
           The need is the figure its paper reports for that dataset where there is one, otherwise the tool's bytes per k-mer from bench/scaling_bench.py runs times the dataset's k-mers.
           Tools that keep to a memory budget are not gated.
           SSHash's -g and KMC's -m get their authors' 16 GiB or less, Cuttlefish's -m half of this machine's memory, and optik dbruijn its own default of half the available memory with cgroup limits included, or in parity mode --max-mem equal to Cuttlefish's -m.
           SBWT's -m is "not strictly enforced" according to its help, so SBWT is gated like the rest.
  cap      A watchdog holds every run's whole process group under the usable memory minus 2 GiB, which --cap-gb lowers.
           bench/memcap.py works on any Linux without systemd or root.
           A run that still exceeds it is killed and reported as CAP, never as a wrong answer.

`list` shows what would run on this machine and why the rest would not, and the same command on a larger machine runs more of the same protocol.
Tiers S, M and L limit a run by size class, and --force ignores the preflight.
cf2-white-spruce (1.1 TB) is opt-in and runs only when --only names it.
A lower cap leaves room for other work, and the memory gate follows it.

  python3 bench/compare_tools_published_sets.py bench --cap-gb 24

Inputs are formatted for every tool with the rules the other benches use, applied while streaming so that multi-GB sets never sit in memory.

  - One upper-case sequence per line, with non-ACGT bases splitting the record as scaling_bench.prep and categories/_base.unwrap_fasta do.
    SSHash reads one sequence per line, and a k-mer across a removed base exists in no sequence.
  - Queries are split into ACGT runs of at least k by categories/_base.write_acgt_segments, one file shared by all tools, because FMSI answers per base rather than per k-mer for a read with N.
  - validate_mapping.sanitize_ref cuts reference headers to their first token, so minimap2 and optik name each sequence alike.
  - The k-mer indexes all read one SPSS, the deposit's eulertigs where the paper supplies them, else those of GGCAT, the papers' tool, else Cuttlefish's unitigs, the common k-mer set of scaling_bench.

`smoke` runs every tool on the first records of the formatted datasets, downloading only a prefix of files not fetched yet, so a format problem shows up before hours of benchmarking.

The sections below follow the measurements of their papers, which is what runs here.

  query  The SSHash 2026 protocol (Pibiri & Patro, Bioinformatics btag264), taken from the paper and its benchmark code at github.com/jermp/kmer_sets_benchmark @ c73a606.
         Eulertigs from GGCAT are the input.
         10^6 positive k-mers are sampled from the set with half of them reverse-complemented, 10^6 more are uniform random, and streaming queries use the paper's high-hit FASTQ per collection.
         Space is measured on disk.
         The FMSI datasets (Sladky et al., Bioinf. Adv. 2025) come from its Zenodo deposit, with the SPSS recomputed at each k as the paper did, at k = 23, the paper's query comparison, and 31, plus the deposited 10%-subsampled sets.
         SBWT (Alanko et al. 2023) contributes the 3682 E. coli genomes of its Zenodo deposit and its gut reads, at k = 31.
  dbg    Cuttlefish 2 (Khan et al., Genome Biol 2022) provides the C. elegans VC2010 reference at f0 = 1, the C. elegans reads DRR008444 at f0 = 12, GRCh38 at f0 = 1 and its Table 1 read sets, all at k = 27 as in its tables.
  map    BWA-MEM (Li 2013, arXiv:1303.3997) used wgsim with 10^6 pairs of 101 bp from the human reference, 1.5% substitution errors, 0.2% indel variants and insert N(500, 50^2), paired and single-end, and aligned E. coli K-12 MG1655, NC_000913, against strain 536, NC_008253.
         minimap2 2018 (Bioinformatics 34:3094) used 10 M mason2 pairs.
         minimap2 2021 (btab705) scored pbsim2 R94 reads from CHM13 v1.1, GCA_009914755.3, with paftools.js mapeval at MAPQ >= 10, and judged CHM13 HiFi reads against Winnowmap2 with paftools.js pafcmp.
  align  BiWFA (Marco-Sola et al., Bioinformatics 2023) aligned the 48 Nanopore MinION pairs of at least 500 kbp shipped in github.com/smarco/BiWFA-paper, and pairs from the paper's generator, WFA2-lib's generate_dataset, at 100 kbp, 1 Mbp and 2 Mbp with 10% and 20% error, with penalties x=4, o=6, e=2, which is the paper's --affine-penalties 0,4,6,2.
         optik's BiWFA is compared only with WFA2-lib's BiWFA in MemoryUltralow mode.
         Full-memory WFA cannot hold these pairs, its memory being quadratic in the score, and is compared on compare_tools_small_sets.py's 150 bp to 100 kbp pairs.

The correctness gates come from external tools, not from this script.

  query  KMC's distinct k-mers of the formatted input equal the count the paper publishes where it does, and any tool reporting a count agrees with KMC.
         On each workload, a tool's positives equal KMC's count of query windows in the set from kmc_tools intersect, and all tools queried an equal number of windows.
         On the positive workload no tool misses a k-mer, and on the negative and streaming workloads all tools report the same positives.
  dbg    Each tool's unitig k-mers, counted by KMC, are compared with KMC's k-mers of the input under the tool's own model, by kmc_tools subtract in both directions, so nothing is missing and nothing is extra.
         For Cuttlefish on reads the model is the k-mers of KMC's (k+1)-mers at f0.
         For Bifrost on reads, extras may only be input k-mers its Bloom filter kept.
         No k-mer may sit in two unitigs, and since maximal unitigs are unique, optik and Cuttlefish, which share the edge-centric model, emit the same unitig and base counts on references.
  map    The papers' own evaluators judge the mappers, namely samtools' wgsim_eval.pl and minimap2's paftools.js mapeval, pafcmp and call.
         All SAMs also go through Picard ValidateSamFile and samtools calmd.
  align  Both BiWFAs must report identical optimal scores, and each CIGAR is re-scored.

Timing goes through bench/measure.py.
Every timed command runs under hyperfine, which discards --warmup runs and measures --runs, by default 1 and 3 for every command alike.
Rows report the mean, standard deviation, median, minimum and maximum over the measured runs, with peak RSS and block I/O per run, and every run starts with the command's inputs already cached in memory.
Correctness is judged on the measured runs' own output.
Per-k-mer query answers come from one untimed run, and the timed runs send them to /dev/null, since counting tools write none.
Builds and graph tools use --threads, by default every logical CPU.
Queries run on one thread, as every query paper measures them, and optik's multi-threaded query is an extra, labelled row.
optik and minimap2 index inside the run, so they also report when their index was ready as index_ready_s and the mapping time after it as mapping_s, while bwa reads a prebuilt index whose build is timed once as bwa_index_s.
The BiWFA rows use the harness's own alignment timer with the same warmup and runs."""

import argparse
import concurrent.futures
import datetime
import logging
import os
import shutil
import subprocess
import sys

import compare_tools_small_sets
import hw
import measure
import report as results_file
from diagnostics import DETAIL_CHARS
from published.align import bench_align
from published.common import (
    DEFAULT_FETCH_JOBS,
    DEFAULT_RUNS,
    FETCH_ERRORS,
    HERE,
    PARITY_T,
    ROUND_FIELDS,
    TOOL_BUDGET_SHARE,
    add_status,
)
from published.dbg import bench_dbg
from published.fetch import build_tools, fetch_one
from published.machine import disk_need, machine, preflight, selected
from published.mapping import bench_map
from published.query import bench_query
from published.smoke import smoke
from published.status import cmd_list, cmd_status, report
from published_catalog import CATALOG, DATA_ROOT, TEMP_DIR
from units import BYTES_PER_GB, MIB_PER_GIB

LOG = logging.getLogger(__name__)


# The memory cap leaves this much of the machine to the system and the bench itself.
CAP_HEADROOM_GIB = 2


def sweep_ks(ds_id, dataset, k_list):
    """The k values a dataset runs at.

    Without --k-list this is the paper's k for the dataset.
    With it, the listed k apply wherever the dataset can be rebuilt at them, which excludes the SSHash 2026 deposits, whose eulertigs exist only at k = 31, and FMSI's 10% subsample at k = 23.
    For map datasets it is the seed k of the as_ sweep sets, and None leaves each mapper at its own default for the paper datasets."""
    if dataset["section"] == "map":
        # Only the matched-k sweep sets.
        # A paper's own map dataset runs its preset once.
        return (k_list or dataset.get("k")) if dataset.get("as_") else None
    if dataset["section"] not in ("query", "dbg"):
        return None
    if not k_list:
        return dataset["k"]
    if ds_id.startswith("sshash26"):
        return [kmer_len for kmer_len in k_list if kmer_len in dataset["k"]] or dataset["k"]
    return k_list


def plan_rounds(ds_id, dataset, base_t):
    """List the rounds --plan paper runs for one dataset.

    Each round is (mode, k list, threads, first_t, sweep_t, label, extra dataset fields).
    The papers' rows run each tool at its optimal settings on all threads at the paper's k (queries add parity rows, shared where the commands coincide).
    k scaling (d["kscale"]) repeats the tools over the k list with PARITY_T threads and parity settings.
    Thread scaling covers hw.thread_points(d["tscale"], base_t) at the paper's k, or at the largest k of a sweep set.
    Sweep and scaling sets have no paper round of their own."""
    sec = dataset["section"]
    rounds = []
    kmer_list = dataset.get("k") if sec in ("query", "dbg") else None
    if not ds_id.startswith(("sweep-", "scale-", "kscale-")):
        for mode in ["parity", "optimal"] if sec == "query" else ["optimal"]:
            rounds.append(
                (mode, kmer_list, 1 if sec == "align" else base_t, True, False, "paper")
            )
    if dataset.get("kscale"):
        # The streaming query, a whole read run on one thread for every tool, happens only at the paper's k, the largest of the sweep.
        # Rounds at other k measure builds, sizes and lookups.
        rounds.append(
            (
                "parity",
                dataset["kscale"],
                PARITY_T,
                True,
                False,
                "kscale",
                {"stream_k": [max(dataset["kscale"])]} if sec == "query" else {},
            )
        )
    if dataset.get("tscale"):
        max_k = [max(kmer_list)] if kmer_list else None
        # Every query paper measures queries on one thread, so thread scaling re-measures only the builds.
        # If the k-scaling round already ran at k0, it holds every tool's queries and the builds of the tools without a thread option.
        # These rounds then rebuild the multi-threaded builders and rerun only optik's queries at t > 1, an extra row for optik's own multi-threaded query.
        covered = sec == "query" and dataset.get("kscale") and max_k and max_k[0] in dataset["kscale"]
        for j, round_threads in enumerate(hw.thread_points(dataset["tscale"], base_t)):
            rounds.append(
                ("parity", max_k, round_threads, j == 0 and not covered, True, "tscale")
                + (({"no_queries": True} if round_threads == 1 else {}) if covered else {},)
            )
    return [planned if len(planned) == ROUND_FIELDS else planned + ({},) for planned in rounds]


def main():
    """Dispatch to list, tools, status, fetch, bench or smoke from the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub_parser = sub.add_parser("list")
    sub_parser.add_argument("--md", action="store_true")
    sub.add_parser("tools")
    sub_parser = sub.add_parser(
        "status", help="download progress of every dataset (fetch running or not)"
    )
    sub_parser.add_argument("--watch", type=float, default=0, help="refresh every N seconds")
    for name in ("fetch", "bench", "smoke"):
        sub_parser = sub.add_parser(name)
        sub_parser.add_argument(
            "--tier",
            default="auto",
            choices=["auto", "S", "M", "L"],
            help="auto (default): every dataset this machine has the disk and tools for",
        )
        sub_parser.add_argument(
            "--force", action="store_true", help="ignore the disk/tool preflight"
        )
        sub_parser.add_argument("--section", default=None, help="query,dbg,map,align")
        sub_parser.add_argument(
            "--only", default=None, help="comma-separated dataset ids (ignores --tier)"
        )
        sub_parser.add_argument(
            "--skip", default=None, help="comma-separated dataset ids to leave out"
        )
        sub_parser.add_argument(
            "--threads",
            type=int,
            default=hw.cpus(),
            help="threads of every parallel tool (default: the usable CPUs)",
        )
        if name in ("bench", "smoke"):
            results_file.add_argument(sub_parser)
        if name == "smoke":
            sub_parser.add_argument(
                "--no-download",
                action="store_true",
                help="only datasets already fetched",
            )
        if name == "fetch":
            sub_parser.add_argument(
                "--jobs", type=int, default=DEFAULT_FETCH_JOBS, help="datasets downloaded at once"
            )
        if name == "bench":
            sub_parser.add_argument(
                "--mode", default="both", choices=["parity", "optimal", "both"]
            )
            sub_parser.add_argument(
                "--runs",
                type=int,
                default=DEFAULT_RUNS,
                help="measured runs of every timed command (statistics over these)",
            )
            sub_parser.add_argument(
                "--warmup",
                type=int,
                default=1,
                help="warmup runs of every timed command, discarded",
            )
            sub_parser.add_argument(
                "--cap-gb",
                type=int,
                default=0,
                help="per-process memory cap (default: usable memory - 2 GiB)",
            )
            sub_parser.add_argument("--work", default=os.path.join(TEMP_DIR, "published_sets_work"))
            sub_parser.add_argument(
                "--keep-work",
                action="store_true",
                help="keep each dataset's working files (default: deleted once its rows are recorded)",
            )
            sub_parser.add_argument(
                "--k-list",
                default=None,
                help="sweep k: every tool at each listed k (e.g. 15,19,21,23,27,31)",
            )
            sub_parser.add_argument(
                "--t-list",
                default=None,
                help="sweep threads: every tool at each listed thread count (e.g. 1,2,4,8,16); "
                "the first should be 1",
            )
            sub_parser.add_argument(
                "--plan",
                default=None,
                choices=["paper"],
                help="paper: the papers' protocol instead of the full matrix -- every dataset once "
                "at each tool's optimal settings on all threads (query: parity rows too, "
                "shared where commands match), dbg parity only where the memory budget binds, "
                "k scaling (kscale) at 8 threads, or every usable CPU when fewer, and thread "
                "scaling (tscale) from 1 up to --threads where a dataset "
                "defines them; queries timed over 5 runs (SSHash 2026), longer runs once",
            )
            sub_parser.add_argument(
                "--skip-rounds",
                default="",
                help="--plan paper: comma-separated DATASET:LABEL rounds to leave out (resuming "
                "a run whose rounds already finished, e.g. sweep-chr1:kscale; their rows "
                "stay in the earlier run's output)",
            )
            sub_parser.add_argument(
                "--dry-run", action="store_true", help="print the rounds and exit"
            )
            sub_parser.add_argument(
                "--out",
                default=os.path.join(HERE, "results", "compare_tools_published_sets"),
            )
    args = parser.parse_args()
    # Tools run inside per-dataset directories, where relative paths would resolve wrongly.
    for opt in ("out", "work"):
        if getattr(args, opt, None):
            setattr(args, opt, os.path.abspath(getattr(args, opt)))
    if args.cmd == "list":
        return cmd_list(args)
    if args.cmd == "tools":
        return build_tools()
    if args.cmd == "status":
        return cmd_status(args)
    if args.cmd in ("bench", "smoke") and not getattr(args, "dry_run", False):
        compare_tools_small_sets.check_optik_fresh()
    mach = machine()
    if getattr(args, "cap_gb", 0):
        # A lower cap leaves the remaining memory to other work, and the memory gate plans within it.
        mach["mem_gib"] = min(mach["mem_gib"], args.cap_gb + CAP_HEADROOM_GIB)
    print(
        f"this machine: {mach['mem_gib']} GiB memory, {mach['disk_free_gb']} GB free disk"
    )
    measure.rusage_exec()  # built up front, or rebuilt when its source is newer
    if args.cmd == "smoke":
        return smoke(args, mach)
    skipped = {}
    for ds_id, dataset in CATALOG.items():
        if "excluded" not in dataset and (
            not args.section or dataset["section"] in args.section.split(",")
        ):
            runnable, why = preflight(ds_id, dataset, mach)
            if not runnable:
                skipped[ds_id] = "; ".join(why)
    todo = selected(args, mach)
    if args.cmd == "fetch":
        # Sources throttle each connection (Zenodo about 250 kB/s) rather than the total, so datasets download side by side.
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            futs = {pool.submit(fetch_one, ds_id, dataset): ds_id for ds_id, dataset in todo}
            bad = []
            for future in concurrent.futures.as_completed(futs):
                try:
                    future.result()
                # Every failed dataset is reported, not just the first.
                except FETCH_ERRORS as fetch_err:
                    bad.append(futs[future])
                    print(f"  {futs[future]}: FAILED ({fetch_err})", flush=True)
        if bad:
            raise SystemExit("fetch failed for: " + ", ".join(bad))
        return
    modes = ["parity", "optimal"] if args.mode == "both" else [args.mode]
    k_list = [int(value) for value in args.k_list.split(",")] if args.k_list else None
    t_list = [int(value) for value in args.t_list.split(",")] if args.t_list else [args.threads]
    if args.t_list and t_list[0] != 1:
        # The compared tools' single-threaded queries run in the first round only.
        raise SystemExit(
            "--t-list must start at 1 (the round that also measures the single-threaded tools)"
        )
    # Unless --threads is given, the plan's paper rounds use every hardware thread.
    plan_t = args.threads
    # A dry run walks the same round plan without running anything, so its listing matches what bench would do.
    if args.dry_run:
        num_rounds = 0
        skip_r = set(filter(None, args.skip_rounds.split(",")))
        for ds_id, dataset in todo:
            rounds = (
                plan_rounds(ds_id, dataset, plan_t)
                if args.plan == "paper"
                else [
                    (mode, sweep_ks(ds_id, dataset, k_list), round_threads, False, False, "", {})
                    for mode in modes
                    for round_threads in t_list
                ]
            )
            for mode, ks_r, round_threads, first_t, _, label, pool in rounds:
                if f"{ds_id}:{label}" in skip_r:
                    print(
                        f"  {dataset['section']:5s} {ds_id:28s} {mode:8s} {label:7s} skipped (--skip-rounds)"
                    )
                    continue
                num_rounds += 1
                what = (
                    f" streaming at k = {pool['stream_k'][0]} only"
                    if pool.get("stream_k") and dataset.get("stream")
                    else " builds only"
                    if pool.get("no_queries")
                    else " builds + optik queries"
                    if label == "tscale" and dataset["section"] == "query" and not first_t
                    else ""
                )
                print(
                    f"  {dataset['section']:5s} {ds_id:28s} {mode:8s} {label:7s} k={ks_r or 'default'} t={round_threads}{what}"
                )
        print(f"{num_rounds} rounds")
        return
    out = os.path.join(
        args.out, datetime.datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")
    )
    os.makedirs(out, exist_ok=True)
    compare_tools_small_sets.environment(out)
    cap = args.cap_gb or max(1, int(mach["mem_gib"] - CAP_HEADROOM_GIB))
    bench = compare_tools_small_sets.Bench(
        out, args.runs, args.warmup, cap, os.path.join(args.work, "_scratch")
    )
    for ds_id, dataset in todo:
        if (
            not os.path.exists(os.path.join(DATA_ROOT.path, ds_id, "PROVENANCE.json"))
            and not dataset.get("local")
            and dataset["files"]
        ):
            print(
                f"  {ds_id}: not fetched (run: compare_tools_published_sets.py fetch --only {ds_id})"
            )
            continue
        # Check free disk again here.
        # The start-up check saw the disk before the earlier datasets wrote and removed their files, and a dataset that cannot fit becomes a SKIP row instead of a tool failing on a full disk partway through.
        need, how = disk_need(ds_id, dataset)
        free = shutil.disk_usage(DATA_ROOT.path).free / BYTES_PER_GB
        if need > free and not args.force:
            print(
                f"\n=== {dataset['section']} {ds_id}: skipped -- needs ~{need:.0f} GB free disk ({how}), {free:.0f} GB free"
            )
            add_status(
                bench,
                "skip",
                section=dataset["section"],
                mode="-",
                dataset=ds_id,
                tool="(all)",
                phase="disk",
                gate_detail=f"needs ~{need:.0f} GB free disk ({how}); {free:.0f} GB free at its turn",
            )
            report(bench, out, mach, skipped, args.report)
            continue
        kmer_list = sweep_ks(ds_id, dataset, k_list)
        if not args.keep_work:
            # Leftovers of an earlier, interrupted run are never reused.
            for entry in os.listdir(args.work) if os.path.isdir(args.work) else []:
                if entry == ds_id or entry.startswith(ds_id + "_"):
                    shutil.rmtree(os.path.join(args.work, entry), ignore_errors=True)

        def run_round(
            mode,
            kmer_list,
            round_threads,
            first_t,
            sweep_t,
            label="",
            only_tools=None,
            extra=None,
            *,
            ds_id=ds_id,
            dataset=dataset,
        ):
            """Run one planned round of dataset `i` and record a failed round as FAIL rows."""
            print(
                f"\n=== {dataset['section']} {ds_id}: {dataset.get('title')} ({mode}{', ' + label if label else ''}, "
                f"k {kmer_list or 'default'}, {round_threads} threads)",
                flush=True,
            )
            round_dataset = dict(dataset, **(extra or {}))
            try:
                if dataset["section"] == "query":
                    bench_query(
                        bench,
                        ds_id,
                        dict(round_dataset, k=kmer_list),
                        mode,
                        args.work,
                        round_threads,
                        mach,
                        first_t=first_t,
                        sweep_t=sweep_t,
                    )
                elif dataset["section"] == "dbg":
                    bench_dbg(
                        bench,
                        ds_id,
                        dict(round_dataset, k=kmer_list),
                        mode,
                        args.work,
                        round_threads,
                        mach,
                        only_tools=only_tools,
                    )
                elif dataset["section"] == "map":
                    for map_k in kmer_list or [None]:
                        bench_map(
                            bench, ds_id, dict(round_dataset, map_k=map_k), mode, args.work, round_threads, mach
                        )
                elif dataset["section"] == "align":
                    bench_align(bench, ds_id, dataset, mode, args.work)
            except Exception as err:
                # A failure while preparing the dataset's input (a crashed tool, an unreadable file) or a bench bug fails only this round.
                # Its rows so far are kept and the run continues.
                # A bench bug prints its traceback so nobody mistakes it for a tool's failure.
                det = f"{type(err).__name__}: {err}" + (
                    f" ({err.output})" if getattr(err, "output", None) else ""
                )
                if not isinstance(
                    err, (subprocess.CalledProcessError, EOFError, OSError, ValueError)
                ):
                    LOG.exception("bench error in round %s of %s", label or mode, ds_id)
                    det = "bench error, " + det
                bench.add(
                    section=dataset["section"],
                    mode=mode,
                    dataset=ds_id,
                    threads=round_threads,
                    tool="(input)",
                    phase="prepare",
                    gate=False,
                    gate_detail=det[-DETAIL_CHARS:],
                )
            report(bench, out, mach, skipped, args.report)

        if args.plan == "paper":
            skip_r = set(filter(None, args.skip_rounds.split(",")))
            for mode, ks_r, round_threads, first_t, sweep_t, label, pool in plan_rounds(ds_id, dataset, plan_t):
                if f"{ds_id}:{label}" in skip_r:
                    print(
                        f"\n=== {dataset['section']} {ds_id}: {label} round skipped (--skip-rounds)",
                        flush=True,
                    )
                    continue
                if dataset["section"] == "map" and label == "kscale":
                    pool = dict(pool, no_bwa=True)
                run_round(mode, ks_r, round_threads, first_t, sweep_t, label, extra=pool)
                if dataset["section"] == "dbg" and label == "paper":
                    # Parity runs only where it changes a command and binds.
                    # optik's --max-mem is set to Cuttlefish's -m budget, run when a tool's optimal peak exceeded it.
                    budget_mb = max(1, int(mach["mem_gib"] * TOOL_BUDGET_SHARE)) * MIB_PER_GIB
                    peaks = [
                        row.get("rss_mb") or 0
                        for row in bench.rows
                        if row.get("dataset") == ds_id
                        and row.get("mode") == "optimal"
                        and row.get("threads") == round_threads
                        and row.get("tool") in ("optik", "cuttlefish")
                    ]
                    # --max-mem only changes optik when a measured optimal peak exceeds budget_mb. Otherwise the parity row is recorded as not run.
                    if any(peak > budget_mb for peak in peaks):
                        run_round(
                            "parity",
                            ks_r,
                            round_threads,
                            True,
                            False,
                            "budget binds: optik --max-mem",
                            only_tools=["optik"],
                        )
                    else:
                        add_status(
                            bench,
                            "not run",
                            section="dbg",
                            mode="parity",
                            dataset=ds_id,
                            threads=round_threads,
                            tool="optik",
                            gate_detail=f"every optimal peak is within the {budget_mb // MIB_PER_GIB} GB budget, "
                            "so --max-mem would not bind (identical to the optimal row)",
                        )
        else:
            for mode in modes:
                if (
                    dataset["section"] in ("dbg", "map", "align")
                    and mode == "parity"
                    and len(modes) > 1
                ):
                    # With --mode both, each tool runs one configuration, the papers' own settings.
                    # For dbg that means the paper's k and f0, with every tool on all threads at its documented best memory setting.
                    # The sweep, --mode parity alone, instead runs dbg with optik's --max-mem equal to Cuttlefish's -m.
                    continue
                for ti, round_threads in enumerate(t_list if dataset["section"] != "align" else [1]):
                    run_round(mode, kmer_list, round_threads, ti == 0, args.t_list is not None)
        if not args.keep_work:
            # With its rows recorded, the dataset's working files (formatted inputs, indexes, SAMs, KMC databases) are removed, so only one dataset's working space is needed at a time.
            for work_entry in os.listdir(args.work):
                if work_entry == ds_id or work_entry.startswith(ds_id + "_"):
                    shutil.rmtree(os.path.join(args.work, work_entry), ignore_errors=True)
    report(bench, out, mach, skipped, args.report)
    print(f"\nwrote {out}/SUMMARY.md")
    # A failed correctness gate fails the run.
    # SKIP and CAP rows have no gate.
    nfail = sum(1 for row in bench.rows if row.get("gate") is False)
    if nfail:
        print(f"{nfail} gate failures")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
