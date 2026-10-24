#!/usr/bin/env python3
"""tool_bench.py - comparative benchmark suite for k-mer indexing tools.

Runs the tools of each category, grouped by underlying mechanism, and compares each with the matching optik variant at that tool's best configuration.

    search-dict   dictionary k-mer lookups (sshash, sbwt) against optik-spss built by build-spss over super-k-mer strings with a minimizer MPHF
    search-fm     FM-index lookups by fmsi with kmercamel against optik-fm, the k-order BWT over the unitigs that build-fm builds
    debruijn      de Bruijn graph construction (bifrost, cuttlefish) against optik dbruijn
    map           sequence mapping (minimap2) against optik map, which needs reads from the genome the unitigs came from

Each tool of the selected categories runs its build and query phases under hyperfine.
Timing and resource statistics are written to JSON and CSV.

Usage
  ./tool_bench.py --list-datasets
  ./tool_bench.py --list-tools
  ./tool_bench.py --dataset ecoli_k31 --tools search-dict
  ./tool_bench.py --dataset salmonella_k31 --all --runs 5 --warmup 2

Results go to bench/results/ by default.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from categories import align, dbg, fmsi, sbwt, ssh
from config import BENCHMARK, DATASETS

# Errors of a round or tool that the run survives, reported with their traceback
LOG = logging.getLogger(__name__)

CSV_COLUMN_NAMES = [
    "category",
    "tool",
    "build_s",
    "build_s_stddev",
    "build_s_min",
    "build_s_max",
    "build_runs",
    "idx_mb",
    "bits_per_kmer",
    "build_rss_mb",
    "query_s",
    "query_s_stddev",
    "query_s_min",
    "query_s_max",
    "query_runs",
    "ns_per_query",
    "queries_s",
    "rand_ns",
    "rand_bench_ns",
    "rand_neg_bench_ns",
    "access_bench_ns",
    "query_rss_mb",
    "query_kmers",
    "query_positive",
    "num_searches",
    # High-hit workload of reads sampled from the indexed sequences, and the tool's own in-memory size accounting, where idx_mb is the bytes on disk
    "pos_query_s",
    "pos_ns_per_kmer",
    "pos_positive_frac",
    "pos_query_kmers",
    "pos_rss_mb",
    "idx_mb_inmem",
    # Secondary parity row on all cores (same subprocess basis, -t threads_query)
    "query_s_allT",
    "ns_per_query_allT",
    "queries_s_allT",
    "rand_ns_allT",
    # Warm in-process metrics, labelled separately from the parity basis
    "warm_ns_per_query",
    "warm_queries_s",
    "warm_rand_ns",
    "warm_query_kmers",
    "warm_query_positive",
    "warm_rand_positive",
    "warm_rand_negative",
    "warm_bench_s",
    "warm_rand_bench_ns",
    "warm_rand_neg_bench_ns",
    "warm_access_bench_ns",
    "warm_query_rss_mb",
    "warm_rand_rss_mb",
    "kmercamel_s",
    # de Bruijn graph metrics
    "num_unitigs",
    "total_length",
    "n50",
    "max_unitig",
    "avg_unitig",
    "min_unitig",
    "valid_gfa",
    "input_kmers",
    "high_confidence_length",
    "distinct_kmers",
    # Alignment metrics
    "mapped_reads",
    "total_reads",
    "reads_per_second",
    "mapping_rate",
    "index_s",
    "index_s_stddev",
    "align_s",
    "align_s_stddev",
    "avg_edit_dist",
    "median_edit_dist",
    "edit_dist_min",
    "edit_dist_max",
    # Concordance metrics
    "concordance_pct",
    "position_agree_pct",
]
"""Ordered list of columns for the CSV output.  Each runner may produce
a subset of these, and missing values become empty cells."""


# Metric keys of each category, in display order.
# Every tool of a category is normalized to exactly this key set and order before printing, so all tools show the same metrics in the same order.
# A metric that does not apply to a tool is shown as n/a, not left out.
CATEGORY_METRIC_ORDER = {
    "search-dict": [
        "build_s",
        "build_s_stddev",
        "build_s_min",
        "build_s_max",
        "build_runs",
        "build_rss_mb",
        "idx_mb",
        "bits_per_kmer",
        "distinct_kmers",
        "query_s",
        "query_s_stddev",
        "query_s_min",
        "query_s_max",
        "query_runs",
        "query_rss_mb",
        "query_kmers",
        "query_positive",
        "query_negative",
        "num_searches",
        "ns_per_query",
        "queries_s",
        "pos_query_s",
        "pos_ns_per_kmer",
        "pos_positive_frac",
        "idx_mb_inmem",
        "rand_ns",
        "rand_rss_mb",
        "rand_positive",
        "rand_bench_ns",
        "rand_neg_bench_ns",
        "access_bench_ns",
        "bench_s",
    ],
    "search-fm": [
        "build_s",
        "build_s_stddev",
        "build_s_min",
        "build_s_max",
        "build_runs",
        "build_rss_mb",
        "idx_mb",
        "bits_per_kmer",
        "distinct_kmers",
        "query_s",
        "query_s_stddev",
        "query_s_min",
        "query_s_max",
        "query_runs",
        "query_rss_mb",
        "query_kmers",
        "query_positive",
        "query_negative",
        "num_searches",
        "ns_per_query",
        "queries_s",
        "pos_query_s",
        "pos_ns_per_kmer",
        "pos_positive_frac",
        "idx_mb_inmem",
        "rand_ns",
        "rand_rss_mb",
        "rand_positive",
        "rand_negative",
        "pos_positive",
        "kmercamel_s",
    ],
    "debruijn": [
        "build_s",
        "build_s_stddev",
        "build_s_min",
        "build_s_max",
        "build_runs",
        "build_rss_mb",
        "idx_mb",
        "num_unitigs",
        "total_length",
        "max_unitig",
        "avg_unitig",
        "min_unitig",
        "n50",
        "valid_gfa",
        "distinct_kmers",
        "input_kmers",
        "high_confidence_length",
    ],
    "map": [
        "build_s",
        "build_s_stddev",
        "build_s_min",
        "build_s_max",
        "build_runs",
        "build_rss_mb",
        "query_rss_mb",
        "idx_mb",
        "index_s",
        "index_s_stddev",
        "align_s",
        "align_s_stddev",
        "mapped_reads",
        "total_reads",
        "reads_per_second",
        "mapping_rate",
        "correct_reads",
        "wrong_reads",
        "unmapped_reads",
        "sensitivity",
        "mapq30_reads",
        "mapq30_wrong",
        "avg_edit_dist",
        "median_edit_dist",
        "edit_dist_min",
        "edit_dist_max",
    ],
}


def normalize_results(results, cat_key):
    """Reorder and complete a runner's result dict to its category's schema.

    Each tool is normalized to the category's metric set (CATEGORY_METRIC_ORDER), so every tool prints the same metrics in the same order.
    A missing metric becomes None and prints as n/a, and keys specific to one tool are kept after the shared ones.

    Error dicts of the form {"error": ...} are returned unchanged.
    """
    if "error" in results:
        return results
    order = CATEGORY_METRIC_ORDER.get(cat_key)
    if not order:
        return results
    normalized = {}
    for key in order:
        normalized[key] = results.get(key)
    for key, val in results.items():
        if key not in normalized:
            normalized[key] = val
    return normalized


PARITY_FIELDS = (
    "query_kmers",
    "query_positive",
    "rand_positive",
    "rand_negative",
    "pos_query_kmers",
    "pos_positive",
)


def query_parity(tools):
    """Cross-tool equality of query answers on identical workloads."""
    fields = {}
    for key in PARITY_FIELDS:
        vals = {
            t: r.get(key)
            for t, r in tools.items()
            if isinstance(r, dict) and r.get(key) is not None and "error" not in r
        }
        if len(vals) >= 2:
            fields[key] = vals
    ok = all(len(set(v.values())) == 1 for v in fields.values())
    # A tool that errored answered nothing, so parity is not established, and the report names the tool and the reason instead of a bare FAIL.
    errored = {
        t: str(r.get("error"))[:200]
        for t, r in tools.items()
        if isinstance(r, dict) and "error" in r
    }
    return {
        "pass": ok and bool(fields) and not errored,
        "fields": fields,
        "errored": errored,
    }


def graph_parity(tools, ds=None):
    """Check each dBG tool's output k-mer set, as dumped by categories/dbg.py.

    It must equal exactly the oracle of the tool's documented filter model, computed from the raw input, and every k-mer must occupy exactly one unitig position."""
    report, ok = [], True
    any_tool = next(
        (r for r in tools.values() if isinstance(r, dict) and "graph_kmers_file" in r),
        None,
    )
    if ds is None or any_tool is None:
        return {"pass": False, "report": ["no dBG outputs"]}
    k, c = int(ds["k"]), int(any_tool.get("min_count", ds.get("dbg_min_count", 2)))
    from dbg_oracles import filter_oracles

    kf, ef = filter_oracles(ds.get("dbg_input", ds["queries"]), k, c)
    report.append(f"oracle (c={c}): k-mer filter {len(kf)}, edge filter {len(ef)}")
    model = {
        "optik": ("k-mer filter", kf),
        "bifrost": ("k-mer filter", kf),
        "cuttlefish": ("edge filter", ef),
    }
    for t, r in tools.items():
        if not isinstance(r, dict) or "graph_kmers_file" not in r or t not in model:
            continue
        with open(r["graph_kmers_file"]) as fh:
            got = set(fh.read().split())
        name, want = model[t]
        dup = r.get("graph_kmer_occurrences", 0) - r.get("graph_kmers", 0)
        good = got == want and dup == 0
        ok &= good
        report.append(
            f"[{'ok' if good else 'MISMATCH'}] {t} ({r.get('mode', '')}) == "
            f"{name}: {len(got)} k-mers, only_tool={len(got - want)} "
            f"only_oracle={len(want - got)}, duplicated positions {dup}"
        )
    return {"pass": ok, "report": report}


def _map_runners():
    """Run minimap2 and optik with their defaults on the same truth-annotated reads.

    wfmash, built on MashMap segments of 1 kb by default, is a long-sequence aligner that short-read pipelines do not use.
    It joins only when the simulated reads are 1 kb or longer (OPTIK_BENCH_MAP_LEN), with read settings inside its documented options.
    On 150 bp reads it was no fair baseline, since an earlier se_k31 run gave it 93.6% correct at about 4k reads/s."""

    runners = [align.Minimap2Runner, align.OptikMapRunner]
    if int(os.environ.get("OPTIK_BENCH_MAP_LEN", "150")) >= 1000:
        runners.append(align.WfmashReadsRunner)
    return runners


TOOL_CATEGORIES = {
    "search-dict": {
        "label": "Dict k-mer lookup",
        "desc": "k-mer dictionaries (sshash, sbwt) vs optik-spss (build-spss)",
        "runners": [
            ssh.SshashRunner,
            sbwt.SbwtFastRunner,
            sbwt.SbwtSmallRunner,
            ssh.OptikSpssRunner,
        ],
        "enabled": lambda: any(
            [
                ssh.SshashRunner.available(),
                sbwt.SbwtFastRunner.available(),
                sbwt.SbwtSmallRunner.available(),
                ssh.OptikSpssRunner.available(),
            ]
        ),
    },
    "search-fm": {
        "label": "FM-index k-mer lookup",
        "desc": "FM-index (kmercamel -S + fmsi) vs optik-fm (k-order BWT over the unitigs, build-fm)",
        "runners": [
            fmsi.FmsiRunner,
            ssh.OptikFMRunner,
        ],
        "enabled": lambda: any(
            [
                fmsi.FmsiRunner.available(),
                ssh.OptikFMRunner.available(),
            ]
        ),
    },
    "debruijn": {
        "label": "De Bruijn graph construction",
        "runners": [
            dbg.BifrostRunner,
            dbg.CuttlefishRunner,
            dbg.OptikDbruijnRunner,
        ],
        "enabled": lambda: any(
            [
                dbg.BifrostRunner.available(),
                dbg.CuttlefishRunner.available(),
                dbg.OptikDbruijnRunner.available(),
            ]
        ),
    },
    "map": {
        "label": "Sequence mapping",
        # Reads are simulated with a fixed seed and a truth file from the dataset's genome (align._map_reference), so every mapper sees the same input.
        "desc": "Read mapping: minimap2 -x sr and optik map (wfmash too for reads >= 1 kb) on identical simulated reads from the genome, scored against truth",
        "runners": _map_runners(),
        "enabled": lambda: any(
            [
                align.Minimap2Runner.available(),
                align.OptikMapRunner.available(),
                align.WfmashRunner.available(),
            ]
        ),
    },
}


def cpu_freq_info():
    """Current cpufreq governor and frequency in MHz of CPU 0.

    The sysfs cpufreq files are read so that results carry the frequency state they were measured under.
    Returns (governor, mhz), with None for any value the machine does not expose, as without cpufreq.
    """
    base = Path("/sys/devices/system/cpu/cpu0/cpufreq")
    governor = mhz = None
    try:
        gov = (base / "scaling_governor").read_text().strip()
        governor = gov if gov else None
    except OSError:
        pass
    for fname in ("scaling_cur_freq", "cpuinfo_cur_freq"):
        try:
            cur = (base / fname).read_text().strip()
            if cur:
                mhz = round(int(cur) / 1000)
                break
        except OSError:
            continue
    return governor, mhz


def list_datasets():
    """Print a table of all configured datasets."""
    print(f"{'Name':<20} {'Label':<25} {'k':>4}  Description")
    print("-" * 80)
    for name, ds in DATASETS.items():
        print(f"{name:<20} {ds['label']:<25} {ds['k']:>4}  {ds.get('desc', '')}")


def list_tools():
    """Print a table showing every registered tool and its availability."""
    print(f"{'Category':<15} {'Tool':<20} {'Available':<12}  Binary")
    print("-" * 75)
    for cat_key, cat in TOOL_CATEGORIES.items():
        for runner_cls in cat["runners"]:
            runner = runner_cls(None, None)
            avail = "yes" if runner.available() else "no"
            binary = runner_cls.binary_path() or "n/a"
            print(f"{cat_key:<15} {runner_cls.tool_name():<20} {avail:<12}  {binary}")


def validate_dataset(dataset_name):
    """Check that the dataset exists and that its files are on disk."""
    if dataset_name not in DATASETS:
        print(f"error: unknown dataset '{dataset_name}'")
        print(f"available datasets: {', '.join(DATASETS.keys())}")
        sys.exit(1)
    ds = DATASETS[dataset_name]
    missing = []
    for key in ["unitigs", "queries"]:
        p = Path(ds[key])
        if not p.exists():
            missing.append(str(p))
    if "dbg_input" in ds:
        p = Path(ds["dbg_input"])
        if not p.exists():
            missing.append(str(p))
    if missing:
        print(f"error: dataset '{dataset_name}' has missing files:")
        for m in missing:
            print(f"  {m}")
        print("edit config.py or download the required files")
        sys.exit(1)
    return ds


def resolve_categories(names):
    """Map a comma-separated category list to the runner configuration of each category."""
    selected = {}
    for name in names.split(","):
        name = name.strip()
        if name not in TOOL_CATEGORIES:
            print(f"error: unknown category '{name}'")
            print(f"available: {', '.join(TOOL_CATEGORIES.keys())}")
            sys.exit(1)
        selected[name] = TOOL_CATEGORIES[name]
    return selected


def format_metric(value, fmt=".2f", placeholder="n/a"):
    """Format a metric value, or return a placeholder for None."""
    if value is None:
        return placeholder
    return f"{value:{fmt}}"


def print_metrics(results):
    """Print each metric of a runner result dict on its own line.

    None prints as n/a, so every tool of a category shows the same metric lines in the same order.
    """
    for key, val in results.items():
        if val is None:
            print(f"    {key}: n/a")
        elif isinstance(val, float):
            print(f"    {key}: {val:.4f}")
        elif isinstance(val, (int, str)):
            print(f"    {key}: {val}")
        else:
            print(f"    {key}: {val}")


def print_summary(results):
    """Print a summary table across all tools and categories."""
    print()

    has_dbg = any(
        "num_unitigs" in tool_res
        for cat in results.get("categories", {}).values()
        for tool_res in cat.get("tools", {}).values()
    )
    has_map = any(
        "mapped_reads" in tool_res
        for cat in results.get("categories", {}).values()
        for tool_res in cat.get("tools", {}).values()
    )
    has_thread_sweep = any(
        isinstance(tool_res.get("build_sweep_rows"), list)
        or isinstance(tool_res.get("query_sweep_rows"), list)
        or isinstance(tool_res.get("warm_sweep_rows"), list)
        or isinstance(tool_res.get("align_sweep_rows"), list)
        for cat in results.get("categories", {}).values()
        for tool_res in cat.get("tools", {}).values()
    )

    # Column widths follow the longest names, where labels used to be cut at 20 characters.
    cw = max(
        [len("Category")]
        + [len(c.get("label", "")) for c in results.get("categories", {}).values()]
    )
    tw = max(
        [len("Tool")]
        + [
            len(t) + 1
            for c in results.get("categories", {}).values()
            for t in c.get("tools", {})
        ]
    )
    header_parts = [
        f"{'Category':<{cw}}",
        f"{'Tool':<{tw}}",
        f"{'Build(s)':>14}",
        f"{'Idx(MB)':>7}",
        f"{'b/kmer':>6}",
        f"{'RSS(MB)':>7}",
        f"{'Seq(ns)':>7}",
        f"{'Rand(ns)':>7}",
        f"{'Query/s':>10}",
    ]

    if has_thread_sweep:
        header_parts.append(f"{'Spd@T':>7}")
        header_parts.append(f"{'RSSxT':>7}")

    if has_dbg:
        header_parts.extend(
            [
                f"{'Unitigs':>8}",
                f"{'TotalLen':>10}",
                f"{'MaxU':>7}",
                f"{'AvgU':>7}",
                f"{'N50':>6}",
            ]
        )

    if has_map:
        header_parts.extend(
            [
                f"{'Mapped':>7}",
                f"{'Rate%':>6}",
                f"{'R/s':>10}",
                f"{'AvgED':>6}",
                f"{'MedED':>6}",
                f"{'Corr%':>6}",
                f"{'MQ30':>7}",
                f"{'MQ30x':>6}",
            ]
        )

    header = "  ".join(header_parts)
    print(header)
    print("-" * len(header))

    for cat in results.get("categories", {}).values():
        for tool_name, tool_res in cat.get("tools", {}).items():
            if "error" in tool_res:
                print(
                    "  ".join(
                        [
                            f"{cat['label']:<{cw}}",
                            f"{tool_name:<{tw}}",
                            f"{'ERROR':>14}",
                        ]
                    )
                )
                continue

            bs = tool_res.get("build_s")
            bs_sd = tool_res.get("build_s_stddev")
            im = tool_res.get("idx_mb")
            bpk = tool_res.get("bits_per_kmer")
            rss = tool_res.get("build_rss_mb")
            seq_ns = tool_res.get("ns_per_query")
            rand_ns = tool_res.get("rand_ns")
            qps = tool_res.get("queries_s")

            if bs is not None and bs_sd is not None and bs_sd > 0:
                bs_str = f"{bs:.2f} +/- {bs_sd:.2f}"
            elif bs is not None:
                bs_str = f"{bs:.2f}"
            else:
                bs_str = "n/a"

            im_str = format_metric(im, ".2f")
            bpk_str = format_metric(bpk, ".2f")
            rss_str = format_metric(rss, ".1f")
            seq_str = format_metric(seq_ns, ".0f")
            rand_str = format_metric(rand_ns, ".0f")
            qps_str = format_metric(qps, ".0f")

            tool_disp = tool_name

            row_parts = [
                f"{cat['label']:<{cw}}",
                f"{tool_disp:<{tw}}",
                f"{bs_str:>14}",
                f"{im_str:>7}",
                f"{bpk_str:>6}",
                f"{rss_str:>7}",
                f"{seq_str:>7}",
                f"{rand_str:>7}",
                f"{qps_str:>10}",
            ]

            if has_thread_sweep:
                rows = (
                    tool_res.get("build_sweep_rows")
                    or tool_res.get("query_sweep_rows")
                    or tool_res.get("align_sweep_rows")
                    or tool_res.get("warm_sweep_rows")
                    or []
                )
                spd = rows[-1].get("speedup") if rows else None
                rss_scale = (
                    tool_res.get("build_rss_scale")
                    or tool_res.get("query_rss_scale")
                    or tool_res.get("align_rss_scale")
                    or tool_res.get("warm_rss_scale")
                )
                row_parts.append(f"{format_metric(spd, '.2f'):>7}")
                row_parts.append(f"{format_metric(rss_scale, '.2f'):>7}")

            if has_dbg:
                unitigs = tool_res.get("num_unitigs")
                total_len = tool_res.get("total_length")
                max_u = tool_res.get("max_unitig")
                avg_u = tool_res.get("avg_unitig")
                n50 = tool_res.get("n50")
                row_parts.extend(
                    [
                        f"{format_metric(unitigs, '.0f'):>8}",
                        f"{format_metric(total_len, '.0f'):>10}",
                        f"{format_metric(max_u, '.0f'):>7}",
                        f"{format_metric(avg_u, '.1f'):>7}",
                        f"{format_metric(n50, '.0f'):>6}",
                    ]
                )

            if has_map:
                mapped = tool_res.get("mapped_reads")
                rate = tool_res.get("mapping_rate")
                rps = tool_res.get("reads_per_second")
                avg_ed = tool_res.get("avg_edit_dist")
                med_ed = tool_res.get("median_edit_dist")
                row_parts.extend(
                    [
                        f"{format_metric(mapped, '.0f'):>7}",
                        f"{format_metric(rate, '.1f'):>6}",
                        f"{format_metric(rps, '.0f'):>10}",
                        f"{format_metric(avg_ed, '.1f'):>6}",
                        f"{format_metric(med_ed, '.1f'):>6}",
                        # Placement accuracy against the simulation truth, since a mapper that places every read wrongly still shows 100% in Rate%
                        f"{format_metric(None if tool_res.get('sensitivity') is None else 100 * tool_res['sensitivity'], '.2f'):>6}",
                        # Confident reads (MAPQ >= 30), the denominator of MQ30x
                        f"{format_metric(tool_res.get('mapq30_reads'), '.0f'):>7}",
                        f"{format_metric(tool_res.get('mapq30_wrong'), '.0f'):>6}",
                    ]
                )

            print("  ".join(row_parts))


def git_revision():
    """HEAD of the repository, with -dirty appended when tracked files have uncommitted changes, so results are never credited to code that did not produce them."""
    root = Path(__file__).resolve().parent.parent
    try:
        rev = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        ).stdout.strip()
        return (rev + "-dirty") if dirty else rev
    except (OSError, subprocess.SubprocessError):
        return None


def main():
    """Parse the arguments, run the benchmarks and write the results."""
    parser = argparse.ArgumentParser(
        description="tool_bench.py: comparative benchmark for k-mer indexing tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  ./tool_bench.py --list-datasets
  ./tool_bench.py --list-tools
  ./tool_bench.py --dataset ecoli_k31 --all --runs 3
  ./tool_bench.py --dataset ecoli_k31 --tools search-dict --runs 5 --warmup 2
  ./tool_bench.py --dataset ecoli_k31 --tools search-dict,search-fm,debruijn
  ./tool_bench.py --reference unitigs.fa.gz --query reads.fq.gz -k 31 --all
        """,
    )
    parser.add_argument(
        "--list-datasets", action="store_true", help="List available datasets"
    )
    parser.add_argument(
        "--list-tools", action="store_true", help="List available tools"
    )
    parser.add_argument("--dataset", "-d", help="Dataset name to benchmark")
    parser.add_argument(
        "--reference", "-r", help="Reference/unitig FASTA file (plain or .gz)"
    )
    parser.add_argument("--query", "-q", help="Query FASTA/FASTQ file (plain or .gz)")
    parser.add_argument(
        "-k",
        "--kmer",
        type=int,
        help="k-mer length (required when using --reference/--query)",
    )
    parser.add_argument("--tools", "-t", help="Comma-separated category names to run")
    parser.add_argument(
        "--all", "-a", action="store_true", help="Run all available tool categories"
    )
    parser.add_argument(
        "--runs", type=int, default=BENCHMARK["runs"], help="Measurement runs per phase"
    )
    parser.add_argument(
        "--warmup", type=int, default=BENCHMARK["warmup"], help="Warmup runs per phase"
    )
    parser.add_argument("--output", "-o", default=None, help="Output JSON file path")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument(
        "--skip-map",
        action="store_true",
        help="Skip the mapping category (minimap2 + optik map)",
    )
    parser.add_argument(
        "--micro",
        action="store_true",
        help="Also run the kernel micro-benchmarks (bench/micro_bench.py): "
        "per-operation cycles, instructions and cache/TLB/branch misses "
        "from unprivileged perf counters, nanoBench method, no root",
    )
    parser.add_argument(
        "--micro-quick",
        action="store_true",
        help="With --micro: cache-resident sizes only",
    )
    parser.add_argument(
        "--pin",
        action="store_true",
        help="Pin all benchmarked commands to one thread per physical core (taskset); "
        "stabilises wall-clock and cycle measurements against SMT/frequency jumps",
    )
    parser.add_argument(
        "--threads-sweep",
        default=BENCHMARK.get("threads", {}).get("query_sweep", "auto"),
        help="Thread-count sweep for every benchmarked phase whose tool accepts a "
        "thread count, including the in-process (warm) phases. 'auto' (default) expands to powers of "
        "two up to the hardware threads plus the physical-core count; pass a "
        "comma list, e.g. '1,2,4,8', to pin the sweep points; 'none' disables "
        "the sweep.  Per-thread timings, throughput, and peak-RSS scaling are "
        "recorded.",
    )

    args = parser.parse_args()

    if args.list_datasets:
        list_datasets()
        return

    if args.list_tools:
        list_tools()
        return

    has_direct_files = args.reference and args.query

    if not args.dataset and not has_direct_files:
        parser.print_help()
        print("\nerror: provide --dataset or --reference + --query + --k")
        sys.exit(1)

    if args.dataset and has_direct_files:
        print("error: cannot use --dataset together with --reference/--query")
        sys.exit(1)

    if has_direct_files:
        if not args.kmer:
            print("error: -k/--kmer is required when using --reference/--query")
            sys.exit(1)
        ref_path = Path(args.reference)
        qry_path = Path(args.query)
        missing = []
        if not ref_path.exists():
            missing.append(str(ref_path))
        if not qry_path.exists():
            missing.append(str(qry_path))
        if missing:
            print("error: file(s) not found:")
            for m in missing:
                print(f"  {m}")
            sys.exit(1)
        ds = {
            "label": f"{ref_path.stem} / {qry_path.stem}",
            "unitigs": ref_path.resolve(),
            "queries": qry_path.resolve(),
            "dbg_input": qry_path.resolve(),
            "k": args.kmer,
            "desc": f"custom reference: {ref_path.name}, query: {qry_path.name}",
        }
        dataset_label = ds["label"]
    else:
        ds = validate_dataset(args.dataset)
        dataset_label = args.dataset

    if args.all:
        selected_cats = {k: v for k, v in TOOL_CATEGORIES.items() if v["enabled"]()}
    elif args.tools:
        selected_cats = resolve_categories(args.tools)
    else:
        print("error: specify --tools or --all")
        sys.exit(1)

    if args.skip_map and "map" in selected_cats:
        del selected_cats["map"]
        if not args.all:
            print("  [skip] map category (--skip-map)")
    if args.all and args.skip_map:
        print("  [skip] map category (--skip-map)")

    if not selected_cats:
        print("no tools available in the selected categories")
        print("run './tool_bench.py --list-tools' to see what is installed")
        sys.exit(1)

    safe_label = dataset_label.replace("/", "_").replace(" ", "_")

    ctx = {
        "runs": args.runs,
        "warmup": args.warmup,
        "verbose": args.verbose,
        "tmp_dir": BENCHMARK["tmp_dir"],
        "threads_build": BENCHMARK["threads"]["build"],
        "threads_query": BENCHMARK["threads"]["query"],
        "threads_sweep": args.threads_sweep,
    }

    if args.pin:
        from categories._base import physical_core_list

        phys = physical_core_list()
        if phys:
            ctx["pin_cpus"] = phys
            print(f"  [info] pinning to physical cores: {','.join(map(str, phys))}")
        else:
            print("  [warn] could not enumerate physical cores; --pin ignored")

    results = {
        "metadata": {
            "tool_bench_version": "2.2",
            "dataset": dataset_label,
            "dataset_info": ds,
            "runs": args.runs,
            "warmup": args.warmup,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pin_cpus": ctx.get("pin_cpus"),
            "threads_sweep": args.threads_sweep,
            "git_revision": git_revision(),
            "cpu_governor": cpu_freq_info()[0],
            "cpu_mhz": cpu_freq_info()[1],
        },
        "categories": {},
    }

    benchmark_separator = "=" * 60

    for cat_key, cat in selected_cats.items():
        cat_results = {"label": cat["label"], "tools": {}}
        print(f"\n{benchmark_separator}")
        print(f"  Category: {cat['label']}")
        if "desc" in cat:
            print(f"  {cat['desc']}")
        print(f"{benchmark_separator}")

        for runner_cls in cat["runners"]:
            tool_name = runner_cls.tool_name()
            runner = runner_cls(ds, ctx)
            if not runner.available():
                print(f"  [skip] {tool_name}: binary not found")
                continue
            print(f"\n  tool: {tool_name}")
            try:
                tool_results = runner.run()
                tool_results = normalize_results(tool_results, cat_key)
                cat_results["tools"][tool_name] = tool_results
                print_metrics(tool_results)
            except Exception as exc:
                # One tool's failure is recorded with its traceback, and the other tools still run.
                LOG.exception("  [error] %s", tool_name)
                cat_results["tools"][tool_name] = {"error": str(exc)}

        # Exactness parity requires every tool that answered the same workload to report the same k-mer and positive counts, because timings of a tool that answers differently are not comparable.
        if cat_key in ("search-dict", "search-fm"):
            cat_results["parity"] = query_parity(cat_results["tools"])
            par = cat_results["parity"]
            print(f"\n  parity: {'PASS' if par['pass'] else 'FAIL'}")
            for t, err in par["errored"].items():
                print(f"    [ERROR] {t} produced no answers: {err}")
            if not par["fields"]:
                print("    [n/a] fewer than two tools answered: nothing to compare")
            for key, vals in par["fields"].items():
                mark = "ok " if len(set(vals.values())) == 1 else "MISMATCH"
                print(f"    [{mark}] {key}: {vals}")

        if cat_key == "debruijn":
            cat_results["graph_parity"] = graph_parity(cat_results["tools"], ds)
            gp = cat_results["graph_parity"]
            print(f"\n  graph parity: {'PASS' if gp['pass'] else 'FAIL'}")
            for line in gp["report"]:
                print(f"    {line}")

        # Graph comparison of the debruijn category
        if cat_key == "debruijn":
            tools = cat_results["tools"]
            optik_dbg = tools.get("optik", {})
            if "error" not in optik_dbg and "num_unitigs" in optik_dbg:
                print(
                    f"\n  optik graph: {optik_dbg.get('num_unitigs', 'n/a')} unitigs, "
                    f"{optik_dbg.get('total_length', 'n/a')} bp total, "
                    f"N50={optik_dbg.get('n50', 'n/a')}"
                )
                for ref_tool in ["bifrost", "cuttlefish"]:
                    ref_res = tools.get(ref_tool, {})
                    if "error" not in ref_res and "num_unitigs" in ref_res:
                        unitig_ratio = (
                            round(optik_dbg["num_unitigs"] / ref_res["num_unitigs"], 2)
                            if ref_res["num_unitigs"]
                            else None
                        )
                        length_ratio = (
                            round(
                                optik_dbg["total_length"] / ref_res["total_length"], 2
                            )
                            if ref_res["total_length"]
                            else None
                        )
                        print(
                            f"  vs {ref_tool}: unitig_ratio={unitig_ratio}, length_ratio={length_ratio}"
                        )

        results["categories"][cat_key] = cat_results

    if args.micro:
        # Kernel-level counters sit next to hyperfine's end-to-end timings, from the same revision and machine, in one results file.
        import micro_bench

        print(f"\n{benchmark_separator}")
        print(
            "  Kernel micro-benchmarks (unprivileged perf counters, nanoBench method)"
        )
        print(f"{benchmark_separator}")
        mres = micro_bench.run(quick=args.micro_quick)
        results["micro"] = mres
        print(micro_bench.table(mres))

    print(f"\n{benchmark_separator}")
    print("  Summary")
    print(f"{benchmark_separator}")
    print_summary(results)

    results_dir = Path(__file__).resolve().parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = (
        Path(args.output)
        if args.output
        else results_dir / f"tool_bench_{safe_label}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nresults written to {out_path}")

    csv_path = out_path.with_suffix(".csv")

    # Scalar keys generated per run (per-thread sweep columns, speedups, RSS scales and so on) are collected across all tools, so the CSV header stays complete whichever sweep points and tools were active.
    extra_cols = []
    for cat in results.get("categories", {}).values():
        for tool_res in cat.get("tools", {}).values():
            for k, v in tool_res.items():
                if k in CSV_COLUMN_NAMES:
                    continue
                if isinstance(v, (list, dict)):
                    continue
                if k not in extra_cols:
                    extra_cols.append(k)
    header_cols = list(CSV_COLUMN_NAMES) + sorted(extra_cols)

    with open(csv_path, "w") as csv:
        csv.write(",".join(header_cols) + "\n")
        for cat_key, cat in results.get("categories", {}).items():
            for tool_name, tool_res in cat.get("tools", {}).items():
                if "error" in tool_res:
                    empty_row = [cat_key, tool_name, "error"] + [""] * (
                        len(header_cols) - 2
                    )
                    csv.write(",".join(empty_row) + "\n")
                    continue
                row = [cat_key, tool_name]
                for col in header_cols[2:]:
                    val = tool_res.get(col, "")
                    row.append(str(val) if val is not None else "")
                csv.write(",".join(row) + "\n")
    print(f"csv summary written to {csv_path}")


if __name__ == "__main__":
    main()
