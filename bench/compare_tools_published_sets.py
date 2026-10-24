#!/usr/bin/env python3
"""compare_tools_published_sets.py - optik against the published tools on the datasets their own most recent papers used.

Each dataset runs in the section it came from, query data in query, graph data in dbg and so on, and the way that paper used it.
The same comparisons on small inputs live in bench/compare_tools_small_sets.py.

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
The BiWFA rows use the harness's own alignment timer with the same warmup and runs.
"""

import argparse
import concurrent.futures
import datetime
import glob
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import threading
import time
from pathlib import Path

# Errors of a round or tool that the run survives, reported with their traceback
LOG = logging.getLogger(__name__)
# The errors a download or a source lookup raises, rangefetch's SystemExit among them
FETCH_ERRORS = (OSError, ValueError, KeyError, subprocess.SubprocessError, SystemExit)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import compare_tools_small_sets as pb
import hw
import measure
import rangefetch
import report as results_file
from external_oracles import (
    SAMTOOLS,
    external_sam_checks,
    judge_unitigs,
    kmc_count,
    kmc_db,
    kmc_diff,
    kmc_hits,
    picard_cmd,
    same_unitigs,
)
from published_catalog import (
    BWA,
    DBG_TOOLS,
    K8_URL,
    QUERY_TOOLS,
    SEQAN_TAG,
    TIERS,
    TOOLS,
    D,
    T,
)
from published_fetch import (
    CURL,
    RANGED_MIN,
    _ranged,
    checksum,
    ena_files,
    ena_prefetch,
    file_urls,
    unzip_ncbi,
    verify_by_source,
    zenodo_file,
)
from published_inputs import (
    _records,
    join_mates,
    query_segments,
    random_kmers,
    reference_input,
    sample_kmers,
    spss_input,
    sshash_m,
    windows_in,
)
from published_mapping import (
    MM2_PRESET_KW,
    PAFTOOLS,
    WGSIM_FAIL,
    WGSIM_ROW,
    bwa_index_hg38,
    bwa_index_s,
    count_records,
    index_ready_s,
    mapeval_calibration,
    paftools_call_snvs,
    primaries,
    seq_bytes,
)
from published_tools import (
    NO_BUILD_THREADS,
    OWN_OUTPUT,
    count_of,
    dev_env,
    ggcat_bin,
    ggcat_version,
    have_tool,
    idx_files,
    query_tools,
    tool_path,
)
from small_common import BIN

# OPTIK_BENCH_DATA names another data root (bench/fetch_check.py fetches into a scratch one)
DATA = os.environ.get("OPTIK_BENCH_DATA") or os.path.join(T, "published_sets")


# Machine resources and gates


def machine():
    """Usable memory in GiB, hw.mem_bytes(), i.e. MemTotal or the tightest cgroup memory.max on this process's path such as a systemd-run -p MemoryMax= scope, and the free disk in GB where the datasets go."""
    mem = hw.mem_bytes() / 2**30
    os.makedirs(DATA, exist_ok=True)
    return {
        "mem_gib": round(mem, 1),
        "disk_free_gb": round(shutil.disk_usage(DATA).free / 1e9, 1),
    }


# Peak memory per distinct k-mer for each tool, by k, as bench/scaling_bench.py measured it on human chr1.
# The scaling logs in SCALING_LOGS are read when they exist.
# Otherwise the fallback below applies, taken from the chr1 run of 2026-10-12 at optik e6202c3, with each tool's peak RSS in MB and the distinct k-mer count _kmers behind it.
SCALING_LOGS = [
    os.path.join(T, "scaling_rerun.log"),
    os.path.join(T, "scaling_rerun_k1521.log"),
]
_PHASE_TOOL = {
    ("optik-spss", "build (default m)"): "optik-spss",
    ("optik-fm", "build"): "optik-fm",
    ("fmsi", "build (kmercamel + index)"): "fmsi",
    ("sbwt-fast", "build"): "sbwt",
    ("cuttlefish", "dbg build (--ref)"): "cuttlefish",
    ("bifrost", "dbg build (-r)"): "bifrost",
    ("optik", "dbg build (-c 1, GFA + links)"): "optik",
}
_FALLBACK_BPK = {
    31: {
        "optik-spss": 1358,
        "optik-fm": 959,
        "fmsi": 3415,
        "sbwt": 2333,
        "cuttlefish": 2215,
        "bifrost": 2433,
        "optik": 1225,
        "_kmers": 204_155_258,
    },
    21: {
        "optik-spss": 2077,
        "cuttlefish": 2545,
        "bifrost": 3130,
        "optik": 1437,
        "_kmers": 189_802_223,
    },
    15: {
        "optik-spss": 3463,
        "cuttlefish": 2796,
        "bifrost": 4343,
        "optik": 5191,
        "_kmers": 111_838_235,
    },
}


def bytes_per_kmer():
    """{k: {tool: bytes per distinct k-mer}} from the chr1 rows of the scaling logs, or from the fallback run above."""
    peaks = {}
    for log in SCALING_LOGS:
        if not os.path.exists(log):
            continue
        with open(log) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) < 6 or not t[0].startswith("chr1") or t[5] != "ok":
                    continue
                k = int(t[0].split("k=")[1]) if "k=" in t[0] else 31
                row = peaks.setdefault(k, {})
                if t[1] == "optik-spss" and t[2] == "build (default m)":
                    m = re.search(r"kmers=(\d+)", " ".join(t[6:]))
                    if m:
                        row["_kmers"] = int(m.group(1))
                tool = _PHASE_TOOL.get((t[1], t[2]))
                if tool and t[4]:
                    row[tool] = max(row.get(tool, 0), float(t[4]))
    src = peaks if 31 in peaks and "_kmers" in peaks[31] else _FALLBACK_BPK
    out = {}
    for k, row in src.items():
        if "_kmers" in row:
            out[k] = {
                tool: mb * 2**20 / row["_kmers"]
                for tool, mb in row.items()
                if tool != "_kmers"
            }
    return out, (
        "scaling logs "
        + ", ".join(os.path.basename(x) for x in SCALING_LOGS if os.path.exists(x))
        if src is peaks
        else "chr1 scaling run of 2026-10-12 (fallback table)"
    )


_BPK = None


def est_mem_gb(tool, k, kmers):
    """Estimated peak in GB of `tool` building over `kmers` distinct k-mers.

    The tool's measured bytes per k-mer are taken at the nearest measured k at or below k, since shorter k costs more per k-mer, and at k = 31 when nothing smaller was measured.
    None if the tool was never measured."""
    global _BPK
    if _BPK is None:
        _BPK = bytes_per_kmer()
    table, _ = _BPK
    ks = sorted(table)
    for kk in [x for x in ks if x <= k][::-1] + [31]:
        if tool in table.get(kk, {}):
            return table[kk][tool] * kmers / 1e9 + 0.2
    return None


def tool_skip(dd, tool, mach, k=None, kmers=None):
    """Why `tool` cannot run on this dataset here, or None.

    A memory need the paper reports for this dataset wins over the measured estimate.
    Tools that keep to a memory budget (BUDGETED) are given this machine's and are never skipped for it.
    SBWT's -m is "not strictly enforced", so SBWT is estimated like the others."""
    if tool in dd.get("skip_tools", {}):
        return dd["skip_tools"][tool]
    allow_gb = mach["mem_gib"] * 2**30 / 1e9
    pub = dd.get("mem", {}).get(tool)
    if pub:
        if pub[0] > allow_gb:
            return (
                f"needs {pub[0]} GB ({pub[1]}); this machine allows {allow_gb:.1f} GB"
            )
        return None
    if tool in BUDGETED or not kmers:
        return None
    est = est_mem_gb(tool, k or 31, kmers)
    if est and est > allow_gb:
        return (
            f"estimated {est:.1f} GB ({kmers / 1e6:.0f} M k-mers x measured bytes per k-mer, "
            f"{_BPK[1]}); this machine allows {allow_gb:.1f} GB"
        )
    return None


# Tools given this machine's memory budget that keep to it, namely SSHash (-g), optik dbruijn (--max-mem, where a graph over it is a named error) and Cuttlefish (-m).
BUDGETED = ("sshash", "optik", "cuttlefish")

_SIZES = os.path.join(DATA, "_download_sizes.json")


def download_gb(i, d):
    """Bytes the sources publish for the dataset's files (Zenodo, ENA or the HTTP Content-Length), cached.
    None where a source gives no size."""
    cache = json.loads(Path(_SIZES).read_text()) if os.path.exists(_SIZES) else {}
    if i in cache:
        return cache[i]
    tot = 0
    try:
        ena_prefetch([f["ena"] for f in d["files"] if "ena" in f])
        for f in d["files"]:
            if "zenodo" in f:
                tot += zenodo_file(f["zenodo"], f["key"])[3]
            elif "ena" in f:
                tot += sum(x[3] for x in ena_files(f["ena"]))
            elif "url" in f:
                hdr = subprocess.run(
                    ["curl", "-fsSLI", "--retry", "4", f["url"]],
                    capture_output=True,
                    text=True,
                    check=False,
                ).stdout
                sizes = re.findall(r"(?im)^content-length:\s*(\d+)", hdr)
                tot += int(sizes[-1]) if sizes else 0
            elif "ncbi" in f:
                tot += 1e9  # a ~3 Gbp assembly, zipped
    except FETCH_ERRORS as e:
        print(f"  {i}: size lookup failed ({e})")
        return None
    cache[i] = round(tot / 1e9, 2)
    Path(_SIZES).write_text(json.dumps(cache, indent=1))
    return cache[i]


def present_gb(i):
    """GB already on disk under dataset `i`'s directory."""
    p = os.path.join(DATA, i)
    if not os.path.isdir(p):
        return 0.0
    return (
        sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(p) for f in fs)
        / 1e9
    )


# Index bits per k-mer of each tool (chr1, k = 31, the same scaling run), which tells what one dataset's indexes take on disk next to their input.
INDEX_BITS = {
    "optik-spss": 8.37,
    "sshash": 7.91,
    "sbwt": 10.54,
    "fmsi": 3.28,
    "optik-fm": 5.43,
}


def disk_parts(i, d):
    """(download GB still to fetch, working GB while the dataset runs, how).

    Compressed FASTA expands about 4.5x, because gzip and xz both store DNA near 2 bits per base, i.e. 8 / 1.8 bytes.
    A k-mer set's SPSS holds about one base per k-mer, and SSHash and KMC write temporary files about the size of their input.
    A mapping run holds its reads, one SAM per mapper at about 2 kB per simulated 150 bp read pair end over all files, and a copy of the reference."""
    gb = download_gb(i, d) or 0.0
    left = max(0.0, gb - present_gb(i))
    if d.get("as_") and not d["files"]:
        # a scaling set reads its source's files, or a share of them, so it downloads nothing
        gb = (download_gb(d["as_"], D[d["as_"]]) or 0.0) * d.get("runs_frac", 1.0)
        left = 0.0
    if d["section"] == "query":
        kmers = d.get("kmers31") or gb * 4.5e9
        spss = kmers / 1e9
        idx = sum(INDEX_BITS.values()) * kmers / 8e9
        work = (spss + idx) * len(d.get("k", [31])) + spss
        how = f"download {gb:.1f} + SPSS {spss:.1f} + indexes {idx:.1f} + temp {spss:.1f} GB per k"
    elif d["section"] == "dbg":
        # Read sets, as measured on cf2-gut with 48 GB of FASTQ.gz.
        # The temporary files of optik and Cuttlefish 3 peaked at 103-106 GB, 2.2x the input, and the KMC databases the rows are judged by stay on disk beside them, about another 1x.
        work = gb * 3.0 if d.get("kind") == "reads" else gb * 4.5 * 2
        how = (
            f"input {gb:.1f} GB + graph temporaries and KMC databases (3x)"
            if d.get("kind") == "reads"
            else f"download {gb:.1f} GB + decompressed genome and unitigs"
        )
    else:
        local = (
            os.path.getsize(d["local"]) / 1e9
            if d.get("local") and os.path.exists(d["local"])
            else 0.0
        )
        work = gb * 3 + local * 2 + d.get("pairs", 0) * 2 * 2000 / 1e9
        how = f"download {gb:.1f} GB + reads, SAMs and reference copies"
    return left, work, how


def disk_need(i, d):
    """(GB, how) one dataset needs free, counting its remaining download and its working files."""
    left, work, how = disk_parts(i, d)
    return left + work, how


def preflight(i, d, mach):
    """(runnable, reasons) of dataset `i` on this machine."""
    if "excluded" in d:
        return False, ["excluded: " + d["excluded"]]
    why = []
    miss = [t for t in d.get("needs", []) if not have_tool(t)]
    if d.get("spss") and not (
        have_tool("ggcat/bin/ggcat") or os.path.exists(BIN["cuttlefish"])
    ):
        miss.append("GGCAT or Cuttlefish (the SPSS step)")
    if d["section"] == "map":
        # the external checks every mapper's SAM goes through
        if not SAMTOOLS:
            miss.append("samtools (PATH, or bench/build_tools.sh samtools)")
        if not picard_cmd():
            miss.append("Picard (PATH, or bench/build_tools.sh picard; needs java 17+)")
    if d.get("local") and not os.path.exists(d["local"]):
        why.append(f"needs {d['local']} (bench/build_tools.sh hg38 / scaling)")
    if (
        d.get("as_")
        and D[d["as_"]]["files"]
        and not os.path.exists(os.path.join(DATA, d["as_"], "PROVENANCE.json"))
    ):
        why.append(f"needs {d['as_']} fetched (fetch --only {d['as_']})")
    if d.get("stream_from") and not os.path.exists(
        os.path.join(DATA, d["stream_from"], d["stream"])
    ):
        why.append(f"needs {d['stream']} (fetch --only {d['stream_from']})")
    if miss:
        why.append(
            "not installed: "
            + ", ".join(miss)
            + " (on PATH, or `compare_tools_published_sets.py tools` "
            "for GGCAT, pbsim2, k8, Winnowmap/meryl and mason2)"
        )
    if d["files"] and download_gb(i, d) is None:
        why.append(
            "download size could not be determined (source unreachable?); --force to try anyway"
        )
    need, how = disk_need(i, d)
    if need > mach["disk_free_gb"]:
        why.append(
            f"needs ~{need:.0f} GB free disk ({how}), {mach['disk_free_gb']:.0f} GB free"
        )
    return not why, why


def selected(args, mach):
    """The (id, dataset) pairs this invocation runs.

    Applies --only, --skip, --section and --tier, then each dataset's preflight, and finally checks that the selection fits the free disk as a whole (unless --force)."""
    out, dl, peak = [], 0.0, 0.0
    only = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if getattr(args, "skip", None) else set()
    unknown = (skip | (only or set())) - set(D)
    if unknown:
        raise SystemExit("unknown dataset ids: " + ", ".join(sorted(unknown)))
    for i, d in D.items():
        if (only is not None and i not in only) or i in skip:
            continue
        if args.section and d["section"] not in args.section.split(","):
            continue
        if "excluded" in d:
            continue
        # An opt-in dataset runs only when --only names it, never through --tier or --force.
        if "opt_in" in d and only is None:
            continue
        if only is None and args.tier != "auto" and TIERS[d["tier"]] > TIERS[args.tier]:
            continue
        ok, why = preflight(i, d, mach)
        if not ok and not args.force:
            print(f"  {i}: skipped on this machine -- {'; '.join(why)}")
            continue
        # The selection must fit as a whole.
        # Downloads stay and add up, while working files go after each dataset, so only the largest of those is needed at once.
        left, work, how = disk_parts(i, d)
        if dl + left + max(peak, work) > mach["disk_free_gb"] and not args.force:
            print(
                f"  {i}: skipped -- needs ~{left + work:.0f} GB ({how}); "
                f"{mach['disk_free_gb'] - dl - peak:.0f} GB left after the datasets selected before it"
            )
            continue
        dl += left
        peak = max(peak, work)
        out.append((i, d))
    return out


# Allocation failures only.
# A log line with "killed" or "oom" in it (a tool killed by another signal, "killed 12 tips") is not a memory cap, which memcap's own flag decides.
MEM_SIG = re.compile(
    r"std::bad_alloc|cannot allocate memory|out of memory|memory allocation failed",
    re.IGNORECASE,
)


def failure(bench, rc_, rss, err):
    """The status row of a failed run.

    A run that died against the memory cap hit a resource limit of this machine and is reported as CAP.
    Any other non-zero exit, an abort without an allocation failure included, is a FAIL.
    The cap applies to the resident memory of the run's whole process group (memcap's watchdog), so CAP means the tool needed more than this machine allows."""
    if bench.last_capped:
        return (
            "cap",
            f"resident memory reached the {bench.cap >> 30} GiB cap (exit {rc_}, peak {rss:.0f} MB)",
        )
    m = MEM_SIG.search(err[-2000:])
    if m:
        return (
            "cap",
            f"memory allocation failed ({m.group(0)!r}, exit {rc_}, peak {rss:.0f} MB)",
        )
    return "fail", f"exit {rc_}: {err[-300:]}"


def add_status(bench, status, **row):
    """Add a SKIP, CAP or NOT RUN row, which is neither a pass nor a wrong answer."""
    row["status"] = status
    row["gate"] = None
    bench.rows.append(row)
    print(
        f"  [{status.upper()}] "
        + "  ".join(
            f"{k}={pb._fmt(v)}"
            for k, v in row.items()
            if k not in ("section", "mode", "gate", "status")
        ),
        flush=True,
    )


# Fetching


def _verified_copy(url, name, dest):
    """A complete, verified copy of url's file under another dataset's directory whose PROVENANCE.json records the same URL as verified, or None."""
    for prov in glob.glob(os.path.join(DATA, "*", "PROVENANCE.json")):
        try:
            with open(prov) as f:
                recs = json.load(f).get("files", [])
        except (OSError, ValueError):
            continue
        for r in recs:
            p = r.get("path", "")
            if (
                r.get("url") == url
                and r.get("verified")
                and os.path.basename(p) == name
                and os.path.isfile(p)
                and not os.path.samefile(os.path.dirname(p), os.path.dirname(dest))
            ):
                return p
    return None


def download(url, dest, algo=None, digest=None, size=None):
    """Fetch url to dest with curl, resuming partial files.

    The file is verified against the source's checksum when the source publishes one (Zenodo, ENA).
    Otherwise its MD5 is recorded."""
    if (
        os.path.exists(dest)
        and (size is None or os.path.getsize(dest) == size)
        and (digest is None or checksum(dest, algo) == digest)
    ):
        return {
            "url": url,
            "path": dest,
            "verified": digest is not None,
            "cached": True,
        }
    same = _verified_copy(url, os.path.basename(dest), dest)
    if same:
        # The same file of the same source was already verified for another dataset, as when CHM13's NCBI zip serves two, so it is linked instead of fetched again.
        for stale in (dest, dest + ".part"):
            if os.path.exists(stale):
                os.remove(stale)
        os.link(same, dest)
        print(
            f"    {os.path.basename(dest)}: linked from {os.path.relpath(same, DATA)}",
            flush=True,
        )
        return {
            "url": url,
            "path": dest,
            "verified": True,
            "cached": True,
            "linked_from": same,
        }
    part = dest + ".part"
    # -sS because several downloads run at once (fetch --jobs) and progress bars would interleave
    print(
        f"    downloading {os.path.basename(dest)}"
        + (f" ({size / 1e9:.2f} GB)" if size else ""),
        flush=True,
    )
    total = size if size is not None else rangefetch.content_length(url)
    if total is not None and total >= RANGED_MIN and _ranged(url):
        # Many paths limit throughput per connection (2026-10-16 through a VPN, about 1.5 MB/s each and 10 MB/s over 8 or more).
        # Large files are therefore split into byte ranges spread over connections that every download shares (bench/rangefetch.py).
        rangefetch.scheduler(
            log=print, progress=os.path.join(DATA, "_fetch_progress.json")
        ).fetch(url, part, total)
    else:
        # A stalled transfer resumes with -C - where the server serves ranges.
        # Where it does not, as NCBI Datasets builds its zips per request, the transfer restarts from zero, because curl refuses to resume there and each retry would fail the same way.
        resume = total is not None and _ranged(url)
        for attempt in range(8):
            if not resume and os.path.exists(part):
                os.remove(part)
            if (
                subprocess.run(
                    CURL + (["-C", "-"] if resume else []) + ["-o", part, url],
                    check=False,
                ).returncode
                == 0
            ):
                break
            time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit(f"{dest}: download failed after 8 attempts")
    if size is not None and os.path.getsize(part) != size:
        raise SystemExit(f"{dest}: size {os.path.getsize(part)} != published {size}")
    got = checksum(part, algo or "md5")
    if digest is not None and got != digest:
        os.remove(part)
        raise SystemExit(f"{dest}: {algo} {got} != published {digest}")
    os.replace(part, dest)
    return {
        "url": url,
        "path": dest,
        algo or "md5": got,
        "verified": digest is not None,
        "date": datetime.datetime.now().astimezone().date().isoformat(),
    }


def fetch_one(i, d):
    """Download, verify and unpack every file of dataset `i`, then write its PROVENANCE.json."""
    if "excluded" in d:
        print(f"  {i}: excluded -- {d['excluded']}")
        return
    out = os.path.join(DATA, i)
    os.makedirs(out, exist_ok=True)
    prov = {
        "id": i,
        "title": d.get("title"),
        "paper": d["paper"],
        "source": d["source"],
        "use": d.get("use"),
        "files": [],
    }
    files = file_urls(d)
    # All files go to the range scheduler together, so a connection moves on to the next file's ranges while the last ranges of a file finish.
    with concurrent.futures.ThreadPoolExecutor(max(1, len(files))) as ex:
        recs = list(
            ex.map(lambda f: download(f[0], os.path.join(out, f[1]), *f[2:]), files)
        )
    prov["files"] += recs
    for (_url, name, *_), rec in zip(files, recs, strict=False):
        if name.endswith(".zip"):
            fa = os.path.join(out, name[:-4] + ".fa")
            src_fa = rec.get("linked_from", "")[:-4] + ".fa"
            if (
                not os.path.exists(fa)
                and rec.get("linked_from")
                and os.path.isfile(src_fa)
            ):
                os.link(src_fa, fa)  # the linked zip's own extraction
            unzip_ncbi(os.path.join(out, name), fa)
    if d.get("local"):
        if not os.path.exists(d["local"]):
            raise SystemExit(
                f"{i}: local input {d['local']} missing (bench/build_tools.sh hg38)"
            )
        prov["files"].append({"path": d["local"], "local": True})
    for rec in prov["files"]:
        if not rec.get("local") and not rec.get("verified"):
            rec["verified_by"] = verify_by_source(rec["url"], rec["path"])
            rec["verified"] = True
    Path(os.path.join(out, "PROVENANCE.json")).write_text(json.dumps(prov, indent=1))
    print(f"  {i}: ok ({len(prov['files'])} files)", flush=True)


# Tool builds


def build_tools():
    """Build the tools this bench needs beyond bench/build_tools.sh.

    These are GGCAT, PBSIM2, k8 to run paftools.js, Winnowmap with meryl, and Mason's simulator, all under TOOLS and linked there by name."""
    os.makedirs(TOOLS, exist_ok=True)
    if not os.path.exists(ggcat_bin()):
        # Latest default branch.
        # The deposit's eulertigs came from 14b2853, and the k-mer set is the same for every version, which the KMC gate checks.
        subprocess.run(
            [
                "cargo",
                "install",
                "--locked",
                "--git",
                "https://github.com/algbio/ggcat",
                "--root",
                os.path.join(TOOLS, "ggcat"),
                "ggcat_cmdline",
            ],
            check=True,
        )
    pb2 = os.path.join(TOOLS, "pbsim2")
    if not os.path.exists(os.path.join(pb2, "src/pbsim")):
        if not os.path.exists(pb2):
            subprocess.run(
                ["git", "clone", "-q", "https://github.com/yukiteruono/pbsim2", pb2],
                check=True,
            )
        # A fresh clone's timestamps make `make` rerun autoconf, which then lacks config.guess.
        subprocess.run(
            "autoreconf -fi && ./configure && make -j8", shell=True, cwd=pb2, check=True
        )
    k8 = os.path.join(TOOLS, "k8")
    if not os.path.exists(k8):
        tb = os.path.join(TOOLS, "k8.tar.bz2")
        download(K8_URL, tb)
        with tarfile.open(tb) as tf:
            m = next(x for x in tf.getmembers() if x.name.endswith("k8-x86_64-Linux"))
            m.name = "k8"
            tf.extract(m, TOOLS, filter="data")
        os.chmod(k8, 0o755)
    wm = os.path.join(TOOLS, "winnowmap_src")
    if not os.path.exists(os.path.join(wm, "bin/winnowmap")):
        if not os.path.exists(wm):
            subprocess.run(
                [
                    "git",
                    "clone",
                    "-q",
                    "--depth",
                    "1",
                    "https://github.com/marbl/Winnowmap",
                    wm,
                ],
                check=True,
            )
        # Older meryl copies include <sys/sysctl.h>, which glibc 2.32 and later no longer ship.
        # meryl uses it only on BSD and macOS, so each such include gets a guard.
        for r, _, fs in os.walk(os.path.join(wm, "ext")):
            for f in fs:
                fp = os.path.join(r, f)
                if not f.endswith((".C", ".H", ".c", ".h", ".cc", ".hpp", ".cpp")):
                    continue
                txt = Path(fp).read_text(errors="replace")
                if (
                    "#include <sys/sysctl.h>\n" in txt
                    and "__has_include(<sys/sysctl.h>)" not in txt
                ):
                    Path(fp).write_text(
                        txt.replace(
                            "#include <sys/sysctl.h>\n",
                            "#if __has_include(<sys/sysctl.h>)\n#include <sys/sysctl.h>\n#endif\n",
                        )
                    )
        subprocess.run(
            ["make", "-j8"], cwd=wm, check=True, env=dev_env(["curl", "openssl"])
        )
    sq = os.path.join(TOOLS, "seqan_src")
    mb = os.path.join(TOOLS, "mason_build")
    if not os.path.exists(os.path.join(mb, "bin/mason_simulator")):
        if not os.path.exists(sq):
            subprocess.run(
                [
                    "git",
                    "clone",
                    "-q",
                    "--depth",
                    "1",
                    "--branch",
                    SEQAN_TAG,
                    "https://github.com/seqan/seqan",
                    sq,
                ],
                check=True,
            )
        os.makedirs(mb, exist_ok=True)
        subprocess.run(
            [
                "cmake",
                sq,
                "-DCMAKE_BUILD_TYPE=Release",
                "-DSEQAN_BUILD_SYSTEM=APP:mason2",
            ],
            cwd=mb,
            check=True,
        )
        subprocess.run(["make", "-j8", "mason_simulator"], cwd=mb, check=True)
    for name, target in (
        ("winnowmap", os.path.join(wm, "bin/winnowmap")),
        ("meryl", os.path.join(wm, "bin/meryl")),
        ("mason_simulator", os.path.join(mb, "bin/mason_simulator")),
    ):
        lnk = os.path.join(TOOLS, name)
        if not os.path.lexists(lnk):
            os.symlink(target, lnk)
    print(
        f"ggcat {ggcat_bin()}\npbsim2 {pb2}/src/pbsim\nk8 {k8}\nwinnowmap/meryl {wm}/bin\n"
        f"mason_simulator {mb}/bin/mason_simulator"
    )


# SPSS formatting
# Streams the rules the other benches use, as the module docstring lists them.


SPSS_TRIES = 3  # GGCAT eulertig attempts before its unitigs stand in (GGCAT issue #65)


def make_spss(inputs, k, out, threads, unitigs=None):
    """Build the SPSS of a k-mer set that every index is built from.

    GGCAT eulertigs from the papers' tool at its latest version are used when GGCAT is built, and otherwise Cuttlefish's unitigs, the common k-mer set of scaling_bench.
    `unitigs` gives the reason GGCAT's plain unitigs stand in for its eulertigs.
    An SPSS is accepted only when its k-mer set equals the input's, meaning KMC finds the same distinct count and nothing outside the input.
    GGCAT's eulertig step aborts on some sets and silently drops k-mers on others, as reported in issue #65, and chm13 at k = 23 lost 2 of 2,337,245,516 in one run.
    A failed attempt is therefore retried, and after SPSS_TRIES failures the unitigs are used.
    Returns the path and the tool that built it."""
    if os.path.exists(out):
        return out, Path(out + ".by").read_text().strip() if os.path.exists(
            out + ".by"
        ) else "?"
    work = os.path.dirname(out)
    kmc_gb = max(1, min(16, int(machine()["mem_gib"] * 0.6)))
    raw_db = kmc_db(inputs, k, 1, os.path.join(work, "kmc_raw"), threads, True, kmc_gb)
    raw_n = kmc_count(raw_db)
    tries = [] if unitigs or not os.path.exists(ggcat_bin()) else [False] * SPSS_TRIES
    notes = []
    for ut in tries + [True]:
        try:
            by = _make_spss(
                inputs,
                k,
                out + ".try",
                threads,
                work,
                ut
                and (
                    unitigs
                    or f"eulertigs: {len(tries)} GGCAT attempts failed, issue #65"
                ),
            )
        except subprocess.CalledProcessError as e:
            tail = (
                Path(out + ".try.log")
                .read_text(errors="replace")
                .strip()
                .splitlines()[-2:]
                if os.path.exists(out + ".try.log")
                else []
            )
            notes.append(
                "aborted"
                + (f" ({tail[-2] if len(tail) > 1 else tail[-1]})" if tail else "")
            )
            if ut:
                raise subprocess.CalledProcessError(
                    e.returncode, e.cmd[:2], output="; ".join(notes)
                ) from None
            continue
        spss_db = kmc_db(
            [out + ".try"], k, 1, os.path.join(work, "kmc_try"), threads, True, kmc_gb
        )
        n, extra = (
            kmc_count(spss_db),
            kmc_diff(spss_db, raw_db, os.path.join(work, "kmc_try_extra")),
        )
        if n == raw_n and extra == 0:
            break
        notes.append(
            f"k-mer set differs from the input's ({n} vs {raw_n} distinct, {extra} not in the input)"
        )
        if ut:
            raise subprocess.CalledProcessError(
                1, ["ggcat", "build"], output="; ".join(notes)
            )
    for ext in ("", ".log"):
        if os.path.exists(out + ".try" + ext):
            os.replace(out + ".try" + ext, out + ext)
    if notes:
        by += (
            f" (accepted on attempt {len(notes) + 1}; earlier: "
            + "; ".join(notes)
            + ")"
        )
    Path(out + ".by").write_text(by + "\n")
    return out, by


def _make_spss(inputs, k, out, threads, work, unitigs=None):
    """Run GGCAT once, or Cuttlefish when GGCAT is absent, and format its output as `out`.
    Returns a description of the tool and mode used."""
    if os.path.exists(ggcat_bin()):
        lst = out + ".list"
        Path(lst).write_text("\n".join(inputs) + "\n")
        with open(out + ".log", "w") as lg:
            subprocess.run(
                [
                    ggcat_bin(),
                    "build",
                    "-k",
                    str(k),
                    "-j",
                    str(threads),
                    "-s",
                    "1",
                    *([] if unitigs else ["--eulertigs"]),
                    "-l",
                    lst,
                    "-o",
                    out + ".raw.fa",
                    "--temp-dir",
                    work,
                ],
                check=True,
                stdout=lg,
                stderr=subprocess.STDOUT,
            )
        by = f"GGCAT {ggcat_version()} " + (
            f"unitigs ({unitigs})" if unitigs else "--eulertigs"
        )
    else:
        with open(out + ".log", "w") as lg:
            subprocess.run(
                [
                    BIN["cuttlefish"],
                    "build",
                    "--ref",
                    *[a for x in inputs for a in ("-s", x)],
                    "-k",
                    str(k),
                    "-t",
                    str(threads),
                    "-o",
                    out + ".cf",
                    "-w",
                    work,
                ],
                check=True,
                stdout=lg,
                stderr=subprocess.STDOUT,
            )
        os.replace(out + ".cf.fa", out + ".raw.fa")
        by = "Cuttlefish unitigs (GGCAT not built)"
    for f in (out, out + ".part"):
        if os.path.exists(f):
            os.remove(f)
    spss_input([out + ".raw.fa"], out, k)
    os.remove(out + ".raw.fa")
    return by


# Query benchmark


def query_inputs(i, dd, k, variant, work, threads, limit_bases=None, nq=10**6):
    """Prepare the formatted SPSS and query files of one (dataset, k, variant).

    Both query modes share them.
    `limit_bases` cuts a slice for `smoke`."""
    src = os.path.join(DATA, i)
    base = os.path.join(work, f"{i}_k{k}_{variant}" + ("_smoke" if limit_bases else ""))
    os.makedirs(base, exist_ok=True)
    unitigs_fa = os.path.join(base, "input.fa")
    by = "deposit (as published)"
    if i.startswith("sshash26"):
        spss_input([os.path.join(src, dd["input"])], unitigs_fa, k, limit_bases)
    elif variant == "sub10":
        spss_input([os.path.join(src, dd["subsampled"])], unitigs_fa, k, limit_bases)
    else:
        if i == "sbwt23-ecoli3682":
            xd = os.path.join(src, "genomes")
            if not os.path.exists(xd):
                with tarfile.open(os.path.join(src, "coli3682_dataset.tar.gz")) as tf:
                    tf.extractall(xd, filter="data")
            raw = sorted(
                os.path.join(r, f)
                for r, _, fs in os.walk(xd)
                for f in fs
                if f.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz"))
            )
        elif dd.get("reads"):
            raw = sorted(
                os.path.join(src, f) for f in os.listdir(src) if f.endswith(".fastq.gz")
            )
            if dd.get("join_mates"):
                joined = os.path.join(base, "spots.fa")
                join_mates(raw, joined, limit_bases)
                raw = [joined]
        elif dd.get("local"):
            raw = [dd["local"]]
        else:
            raw = [os.path.join(src, dd["input"])]
        plain = os.path.join(base, "plain.fa")
        spss_input(raw, plain, k, limit_bases, dd.get("delete_nonacgt", False))
        unitigs_fa, by = make_spss(
            [plain], k, os.path.join(base, "spss.fa"), threads, dd.get("spss_unitigs")
        )
    n = min(nq, 10**4) if limit_bases else nq
    wl = {
        "pos_halfrc": sample_kmers(unitigs_fa, k, n, os.path.join(base, "q_pos.fa")),
        "neg_random": random_kmers(k, n, os.path.join(base, "q_neg.fa")),
        "_load": random_kmers(k, 1, os.path.join(base, "q_load.fa"), seed=3),
    }
    # stream_k (k-scaling rounds) restricts the streaming workload to these k.
    # FMSI's k sweep measures size and lookups, whereas SSHash's streaming numbers are a k = 31 measurement.
    if (
        dd.get("stream")
        and (k == 31 or dd.get("stream_from"))
        and k in dd.get("stream_k", [k])
    ):
        q = os.path.join(base, "q_stream.fa")
        query_segments(
            os.path.join(DATA, dd.get("stream_from", i), dd["stream"]),
            q,
            k,
            2 * 10**6 if limit_bases else None,
        )
        wl["stream_" + dd["stream"].split(".")[0]] = q
    kmers = spss_input([], unitigs_fa, k)[2]
    return unitigs_fa, by, wl, kmers


def bench_query(bench, i, dd, mode, work, threads, mach, first_t=True, sweep_t=False):
    """Run one thread count `threads` of a dataset.

    first_t marks the sweep's first, single-thread round, which also runs the tools with no thread option for a step.
    Later rounds rebuild only the multi-threaded builders and rerun only optik's queries (sweep_t runs optik queries at `threads`)."""
    for k in dd["k"]:
        variants = ["full"] + (["sub10"] if dd.get("subsampled") and k == 23 else [])
        for variant in variants:
            unitigs_fa, by, wl, _windows = query_inputs(
                i, dd, k, variant, work, threads
            )
            load = wl.pop("_load")
            # The gates below judge only the rows this round adds.
            # A k-scaling and a thread-scaling round can share (dataset, k, variant, mode, threads), and the other round's rows were already judged against its own answers.
            row0 = len(bench.rows)
            # Each mode and thread count gets a fresh directory for its indexes and answers, so a tool that exits 0 without writing is never scored on an older file.
            d = os.path.join(work, f"{i}_k{k}_{variant}", f"{mode}_t{threads}")
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            pub = dd.get("kmers31") if k == 31 and variant == "full" else None
            # KMC counts the indexed set's distinct k-mers independently of every compared tool.
            # When the paper publishes a count, the two must agree, which also shows that formatting kept the paper's set.
            kmc_gb = max(1, min(16, int(mach["mem_gib"] * 0.6)))
            kdb = os.path.join(os.path.dirname(d), "kmc_set")
            kmc_n = kmc_count(kmc_db([unitigs_fa], k, 1, kdb, threads, True, kmc_gb))
            kmers = pub or kmc_n
            input_ok = pub is None or kmc_n == pub
            if pub is not None and dd.get("kmers31_forward"):
                # The paper counts k-mers as written, so KMC runs without canonical form on the reads.
                fwd = kmc_count(
                    kmc_db(
                        [os.path.join(os.path.dirname(d), "plain.fa")],
                        k,
                        1,
                        os.path.join(os.path.dirname(d), "kmc_forward"),
                        threads,
                        True,
                        kmc_gb,
                        canonical=False,
                    )
                )
                input_ok = abs(fwd - pub) <= dd.get("kmers31_tol", 0) * pub
                kmers = kmc_n
            if first_t:
                fw = dd.get("kmers31_forward") and pub is not None
                bench.add(
                    section="query",
                    mode=mode,
                    dataset=i,
                    variant=variant,
                    k=k,
                    tool="kmc (input)",
                    phase="input",
                    kmers_reported=kmc_n,
                    **({"kmers_forward": fwd} if fw else {}),
                    kmers_published=pub,
                    gate=input_ok,
                    gate_detail=(
                        f"forward k-mers {fwd} vs the paper's {pub} ({fwd / pub - 1:+.2%}, "
                        f"allowed +-{dd.get('kmers31_tol', 0):.1%}: see the dataset's note)"
                        if fw
                        else ""
                        if input_ok
                        else f"KMC counts {kmc_n} distinct k-mers in "
                        f"the formatted input, the paper publishes {pub}"
                    ),
                    params=f"SPSS: {by}"
                    + ("; reads = spots, mates joined" if dd.get("join_mates") else "")
                    + (
                        "; non-ACGT deleted within records"
                        if dd.get("delete_nonacgt")
                        else ""
                    ),
                )
            expected = {w: windows_in(q, k) for w, q in wl.items()}
            # KMC's answer for every workload is the set of query k-mer windows it finds.
            kmc_truth = {
                w: kmc_hits(
                    kdb,
                    q,
                    k,
                    os.path.join(os.path.dirname(d), f"kmc_{w}"),
                    threads,
                    kmc_gb,
                )
                for w, q in wl.items()
                if w != "_load"
            }
            # SSHash requires -m and has no default, so the paper's m is used where given and k/2 otherwise.
            m = sshash_m(dd, k, kmc_n)
            tools = query_tools(
                mode,
                d,
                unitigs_fa,
                k,
                m,
                threads,
                mach,
                bench.cap >> 30,
                threads if sweep_t else 1,
            )
            res, nk = {}, {}
            for tool, t in tools.items():
                if not first_t and tool in NO_BUILD_THREADS:
                    continue  # without a thread option, its t = 1 row stands
                tag = f"{mode}/{i}/k{k}/{variant}/t{threads}/{tool}"
                # FMSI's build takes no thread option, so its row records what it ran with.
                base = {
                    "section": "query",
                    "mode": mode,
                    "dataset": i,
                    "variant": variant,
                    "k": k,
                    "threads": 1 if tool in NO_BUILD_THREADS else threads,
                    "tool": tool,
                    "answers": t["answers"],
                    "params": t["params"] + f"; SPSS: {by}",
                }
                # A command that is the same in both modes (every compared tool, and optik-fm) is measured once in the parity round, and its rows also stand in the optimal table.
                # optik-spss differs (-m), and its --ids view queries that index.
                same = [
                    r
                    for r in bench.rows
                    if mode == "optimal"
                    and not tool.startswith("optik-spss")
                    and r.get("mode") == "parity"
                    and r.get("section") == "query"
                    and (
                        r.get("dataset"),
                        r.get("k"),
                        r.get("variant"),
                        r.get("threads"),
                        r.get("tool"),
                    )
                    == (i, k, variant, base["threads"], tool)
                    and r.get("params") == base["params"]
                ]
                if same:
                    for r in same:
                        nr = dict(
                            r,
                            mode="optimal",
                            params=r["params"]
                            + "; the parity run's command: measured once",
                        )
                        bench.rows.append(nr)
                        if r.get("phase") in wl and r.get("kmers") is not None:
                            res.setdefault(r["phase"], {})[tool] = (
                                r["kmers"],
                                r["positive"],
                            )
                        if r.get("phase") == "build":
                            nk[tool] = r.get("kmers_reported")
                    continue
                if t["build"] is None:
                    # This is a query mode of another tool's index, with no build of its own.
                    if t["reuse"] not in nk:
                        continue  # that index was skipped or failed
                    lcmd, _ = t["query"](load)
                    idx = idx_files(tools[t["reuse"]], d)
                    tl = (
                        bench.measure(
                            tag + "/load",
                            t.get("timed", lambda q, lcmd=lcmd: lcmd)(load),
                            cwd=d,
                            inputs=[load, *idx],
                        ).wall_row()
                        or {}
                    )
                if t["build"] is not None:
                    why = tool_skip(dd, tool, mach, k, kmers)
                    if why:
                        add_status(
                            bench, "skip", **base, phase="build", gate_detail=why
                        )
                        continue
                    bout = os.path.join(d, f"{tool}.build.out")
                    # The SPSS is read in before each run (measure.py), so all runs see one cache state.
                    rb = bench.measure(
                        tag + "/build",
                        t["build"],
                        cwd=d,
                        stdout=bout,
                        inputs=[unitigs_fa],
                    )
                    err = (
                        Path(bout).read_text() if os.path.exists(bout) else ""
                    ) + rb.err
                    if not rb.ok:
                        kind, det = failure(
                            bench, rb.code, rb.rss() or bench.cap / 2**20, err
                        )
                        if kind == "cap":
                            add_status(
                                bench, "cap", **base, phase="build", gate_detail=det
                            )
                        else:
                            bench.add(
                                **base,
                                phase="build",
                                gate=False,
                                gate_detail="build " + det,
                            )
                        continue
                    idx = idx_files(t, d)
                    size = sum(os.path.getsize(f) for f in idx)
                    n = re.search(r"distinct_kmers\s+(\d+)", err) or re.search(
                        r"num_kmers\s*=?\s*(\d+)", err
                    )
                    nk[tool] = int(n.group(1)) if n else None
                    # A load-only run, which the per-k-mer time below subtracts so that index loading, seconds at this scale, does not count as query cost.
                    lcmd, _ = t["query"](load)
                    tl = (
                        bench.measure(
                            tag + "/load",
                            t["timed"](load) if "timed" in t else lcmd,
                            cwd=d,
                            inputs=[load, *idx],
                        ).wall_row()
                        or {}
                    )
                    extra = {}
                    if t.get("extra"):
                        ecmd, efile = t["extra"]
                        if (
                            bench.once(
                                tag + "/build-norc",
                                ecmd,
                                cwd=d,
                                stdout=os.path.join(d, "norc.out"),
                            )[0]
                            == 0
                        ):
                            extra["index_mb_without_rc"] = round(
                                os.path.getsize(efile) / 2**20, 2
                            )
                    bench.add(
                        **base,
                        phase="build",
                        build_s=rb.wall_row(),
                        build_rss_mb=round(rb.rss(), 1),
                        index_mb=round(size / 2**20, 2),
                        **extra,
                        kmers_reported=nk[tool],
                        load_s=tl,
                        gate=True,
                        gate_detail="",
                    )
                for w, q in wl.items():
                    if dd.get("no_queries"):
                        break  # a build-only round leaves the queries to another round
                    if not first_t and not tool.startswith("optik"):
                        break  # single-threaded queries belong to the t = 1 round
                    cmd, outf = t["query"](q)
                    capt = (
                        os.path.join(d, "q.out")
                        if tool in OWN_OUTPUT
                        else (outf or os.path.join(d, "q.out"))
                    )
                    qidx = idx_files(tools[t["reuse"]] if t["build"] is None else t, d)
                    if t["answers"] == "per-k-mer":
                        # One untimed run writes the per-k-mer answers that the gate reads.
                        # The timed runs discard them in /dev/null, since the counting tools write none.
                        rc_, _, rss, err = bench.once(
                            f"{tag}/{w} (answers)",
                            cmd,
                            cwd=d,
                            stdout=None if tool == "sbwt" else capt,
                        )
                        rq = (
                            bench.measure(
                                f"{tag}/{w}",
                                t["timed"](q) if "timed" in t else cmd,
                                cwd=d,
                                inputs=[q, *qidx],
                            )
                            if rc_ == 0
                            else None
                        )
                    else:
                        rq = bench.measure(
                            f"{tag}/{w}", cmd, cwd=d, stdout=capt, inputs=[q, *qidx]
                        )
                        rc_, rss, err = (
                            (0 if rq.ok else rq.code),
                            rq.rss() or bench.cap / 2**20,
                            rq.err,
                        )
                    if rc_ == 0 and not rq.ok:
                        rc_, rss, err = rq.code, bench.cap / 2**20, rq.err
                    if rc_ != 0:
                        kind, det = failure(bench, rc_, rss, err)
                        if kind == "cap":
                            add_status(bench, "cap", **base, phase=w, gate_detail=det)
                        else:
                            bench.add(
                                **base, phase=w, gate=False, gate_detail="query " + det
                            )
                        continue
                    got = count_of(tool, outf, capt, err)
                    res.setdefault(w, {})[tool] = got
                    bench.add(
                        **base,
                        phase=w,
                        query_s=rq.wall_row(),
                        query_rss_mb=round(rq.rss(), 1),
                        kmers=got[0],
                        positive=got[1],
                        gate=True,
                        gate_detail="(cross-tool gate below)",
                        _load=tl.get("mean"),
                    )
            # Gates on published or cross-tool distinct counts and on per-workload agreement
            ref_n = kmc_n
            for r in bench.rows[row0:]:
                # Only this round's rows count.
                # FMSI's (threads = 1, no thread option) belong to the first round.
                rt = (
                    threads
                    if r.get("tool") not in NO_BUILD_THREADS
                    else (1 if first_t else None)
                )
                if (
                    r.get("dataset"),
                    r.get("k"),
                    r.get("variant"),
                    r.get("mode"),
                    r.get("threads"),
                ) != (i, k, variant, mode, rt) or r.get("section") != "query":
                    continue
                if r.get("status"):
                    continue
                if r.get("phase") == "build" and "index_mb" in r and ref_n:
                    r["bits_per_kmer"] = round(8 * r["index_mb"] * 2**20 / ref_n, 3)
                    if (
                        r.get("kmers_reported") is not None
                        and r["kmers_reported"] != ref_n
                    ):
                        r["gate"] = False
                        r["gate_detail"] = (
                            f"distinct k-mers {r['kmers_reported']} != KMC's {ref_n}"
                        )
                w = r.get("phase")
                if w in res:
                    # A tool's positives must equal KMC's, and it must have queried all windows of the query file.
                    # Matching the others' total is not enough, since one tool alone, or a parser that read nothing, agrees with itself.
                    mine, kq = res[w].get(r["tool"], (None, None)), kmc_truth.get(w)
                    ok = mine[1] == kq and mine[0] == expected[w]
                    r["kmc_positives"] = kq
                    r["gate"] = ok
                    r["gate_detail"] = (
                        ""
                        if ok
                        else (
                            f"positives {mine[1]} != KMC {kq}"
                            if mine[1] != kq
                            else f"queried {mine[0]} windows of the file's {expected[w]}"
                        )
                    )
                    if (
                        "mean" in r.get("query_s", {})
                        and r.get("kmers")
                        and "_load" in r
                    ):
                        net = r["query_s"]["mean"] - (r.pop("_load") or 0.0)
                        r["ns_per_kmer_net"] = round(
                            max(net, 0.0) * 1e9 / r["kmers"], 1
                        )
                r.pop("_load", None)
                if not input_ok and r.get("tool") != "kmc (input)":
                    # All tools agreed with KMC, but the input is not the paper's set.
                    r["gate"] = False
                    r["gate_detail"] = (
                        (r.get("gate_detail") or "")
                        + f"; the formatted input has {kmc_n} "
                        f"distinct k-mers, the paper publishes {pub}"
                    ).lstrip("; ")


# De Bruijn graph benchmark

# The KMC oracles of external_oracles.py, shared with compare_tools_small_sets
KMC = BIN["kmc"]


def dbg_inputs(i, dd, k, work, limit_bases=None):
    """Return the inputs of the graph tools.

    A reference is formatted with one record per line, upper case and first-token headers.
    Read files are passed as fetched, since every tool reads FASTQ.gz, and `smoke` takes a slice of them."""
    src = os.path.join(DATA, dd.get("as_", i))
    base = os.path.join(work, f"{i}_k{k}" + ("_smoke" if limit_bases else ""))
    os.makedirs(base, exist_ok=True)
    if dd["kind"] == "reads":
        fq = sorted(
            os.path.join(src, f) for f in os.listdir(src) if f.endswith(".fastq.gz")
        )
        if dd.get("runs_frac"):
            # A fixed subsample of the first runs, both mates, that hold runs_frac of the compressed bytes
            runs = sorted(
                {re.sub(r"_[12]\.fastq\.gz$", "", os.path.basename(f)) for f in fq}
            )
            total = sum(os.path.getsize(f) for f in fq)
            keep, acc, want = set(), 0, dd["runs_frac"] * total
            for r in runs:
                sz = sum(
                    os.path.getsize(f)
                    for f in fq
                    if os.path.basename(f).startswith(r + "_")
                )
                # The run that crosses the target is kept only if keeping it lands nearer the target.
                if acc + sz > want and acc + sz - want >= want - acc:
                    break
                keep.add(r)
                acc += sz
                if acc >= want:
                    break
            if not keep and runs:
                # A prefix-only set, as in smoke, can be smaller than one run's share, so its first run stands in.
                first = runs[0]
                keep.add(first)
                acc = sum(
                    os.path.getsize(f)
                    for f in fq
                    if os.path.basename(f).startswith(first + "_")
                )
            fq = [
                f
                for f in fq
                if re.sub(r"_[12]\.fastq\.gz$", "", os.path.basename(f)) in keep
            ]
            print(
                f"  {i}: the first {len(keep)} runs of {len(runs)}, {acc / total:.1%} of the compressed bytes "
                f"(target {dd['runs_frac']:.1%}; whole runs, both mates)",
                flush=True,
            )
        if not limit_bases:
            return fq, True
        out = []
        for j, f in enumerate(fq[:2]):
            o = os.path.join(base, f"slice{j}.fq")
            if not os.path.exists(o):
                with open(o, "w") as w:
                    w.writelines(
                        f"@{n}\n{sq}\n+\n{'I' * len(sq)}\n"
                        for n, sq in _records(f, limit_bases // 2)
                    )
            out.append(o)
        return out, True
    srcf = dd["local"] if dd.get("local") else os.path.join(src, dd["input"])
    return [reference_input(srcf, os.path.join(base, "ref.fa"), limit_bases)], False


class DiskPeak:
    """Peak disk space a run takes on the filesystem of `path`.

    The filesystem's used bytes, as statvfs and df report them, are sampled every second and the starting value is subtracted.
    Counting files in the work directory would miss optik's spill files, which O_TMPFILE never links, and Cuttlefish's deleted buckets.
    Because the bench runs one tool at a time, the filesystem's growth belongs to that tool, output file included.
    Cuttlefish 2 reports its intermediate disk use the same way in its Tables S1 and S3."""

    def __init__(self, path):
        self.path, self.peak, self._stop = path, 0, threading.Event()

    def _used(self):
        """Bytes in use on the filesystem."""
        st = os.statvfs(self.path)
        return (st.f_blocks - st.f_bfree) * st.f_frsize

    def __enter__(self):
        self.base = self._used()
        self._t = threading.Thread(target=self._run, daemon=True)
        self._t.start()
        return self

    def _run(self):
        """Sampler thread body, polling once a second until __exit__."""
        while not self._stop.wait(1.0):
            self.peak = max(self.peak, self._used() - self.base)

    def __exit__(self, *a):
        self._stop.set()
        self._t.join()
        self.peak = max(self.peak, self._used() - self.base)


def bench_dbg(
    bench, i, dd, mode, work, threads, mach, limit_bases=None, only_tools=None
):
    """Run every graph tool at the same k, cutoff and thread count.

    On reads the paper's f0 compares only two tools, because optik and Cuttlefish take a cutoff while Bifrost's is fixed at 2.
    A matched round at cutoff 2 therefore adds Bifrost."""
    kmc_gb = max(1, min(16, int(mach["mem_gib"] * 0.6)))
    for k, c in [
        (k, c)
        for k in dd["k"]
        for c in dict.fromkeys(
            [dd["cutoff"]]
            + (
                [2]
                if dd.get("kind") == "reads"
                and not limit_bases
                and not dd.get("no_c2_round")
                else []
            )
        )
    ]:
        inputs, fq = dbg_inputs(i, dd, k, work, limit_bases)
        # Outputs go to a directory per thread count, while the input-only KMC databases are shared.
        shared = os.path.join(
            work, f"{i}_k{k}" + ("_smoke" if limit_bases else ""), f"kmc_c{c}"
        )
        d = os.path.join(
            work,
            f"{i}_k{k}" + ("_smoke" if limit_bases else ""),
            f"{mode}_c{c}_t{threads}",
        )
        os.makedirs(d, exist_ok=True)
        os.makedirs(shared, exist_ok=True)
        truth = kmc_db(
            inputs, k, c, os.path.join(shared, "truth"), threads, not fq, kmc_gb
        )
        # A reference's distinct k-mers feed the memory gate.
        # Read sets use published figures only.
        n_truth = None if fq else kmc_count(truth)
        t = threads
        # Tools that take a memory budget share one, half of this machine's (or --cap-gb's) memory.
        # It sets Cuttlefish's -m and, in parity mode, optik's --max-mem.
        # In optimal mode optik keeps its documented default of half of MemAvailable at launch, cgroup limits included.
        # A larger --max-mem only keeps more counting buffers resident at the same speed.
        budget = max(1, int(mach["mem_gib"] * 0.5))
        tools = {
            "optik": (
                [
                    BIN["optik"],
                    "dbruijn",
                    *[a for x in inputs for a in ("-i", x)],
                    "-k",
                    str(k),
                    "-c",
                    str(c),
                    "-t",
                    str(t),
                    "-o",
                    f"{d}/optik.fa",
                ]
                + (["--max-mem", f"{budget}G"] if mode == "parity" else []),
                f"{d}/optik.fa",
                "node",
            ),
            "cuttlefish": (
                [
                    BIN["cuttlefish"],
                    "build",
                    *[a for x in inputs for a in ("-s", x)],
                    "-k",
                    str(k),
                    "-t",
                    str(t),
                    "-m",
                    str(budget),
                    "-o",
                    f"{d}/cf",
                    "-w",
                    d,
                    "--read" if fq else "--ref",
                ]
                + (["-c", str(c)] if fq else []),
                f"{d}/cf.fa",
                "edge",
            ),
            "bifrost": (
                [
                    BIN["bifrost"],
                    "build",
                    *[a for x in inputs for a in (["-s"] if fq else ["-r"]) + [x]],
                    "-k",
                    str(k),
                    "-t",
                    str(t),
                    "-o",
                    f"{d}/bf",
                    "-f",
                    "-n",
                    "-N",
                ],
                f"{d}/bf.fasta",
                "node",
            ),
        }
        stats = {}
        row0 = len(
            bench.rows
        )  # the unitig cross-check below marks only this round's rows
        for tool, (cmd, outp, model) in tools.items():
            if (only_tools and tool not in only_tools) or (
                dd.get("tools") and tool not in dd["tools"]
            ):
                continue
            tag = f"{mode}/{i}/k{k}/c{c}/t{threads}/{tool}"
            params = " ".join(map(str, cmd[1:]))
            if tool == "bifrost" and fq and c != 2:
                continue  # Bifrost's cutoff is fixed at 2, so the c = 2 round compares it
            why = None if limit_bases else tool_skip(dd, tool, mach, k, n_truth)
            if why:
                add_status(
                    bench,
                    "skip",
                    section="dbg",
                    mode=mode,
                    dataset=i,
                    k=k,
                    cutoff=c,
                    threads=threads,
                    tool=tool,
                    gate_detail=why,
                )
                continue
            # Every run, warmups included, starts with the inputs in the page cache (measure.py), since the KMC judging evicts them.
            # The disk peak spans all runs, with each run's output replacing the previous one.
            with DiskPeak(d) as dp:
                r = bench.measure(tag, cmd, cwd=d, inputs=inputs)
            rc_, rss, err = (0 if r.ok else r.code), r.rss() or bench.cap / 2**20, r.err
            wrote, readb = r.written, r.read  # median block I/O of a measured run
            if rc_ != 0 and failure(bench, rc_, rss, err)[0] == "cap":
                add_status(
                    bench,
                    "cap",
                    section="dbg",
                    mode=mode,
                    dataset=i,
                    k=k,
                    cutoff=c,
                    threads=threads,
                    tool=tool,
                    gate_detail=failure(bench, rc_, rss, err)[1],
                    params=params,
                )
                continue
            if rc_ != 0 or not os.path.exists(outp):
                alt = [
                    p
                    for p in (
                        outp + ".gz",
                        outp.replace(".fasta", ".fa"),
                        outp.replace(".fasta", ".gfa"),
                    )
                    if os.path.exists(p)
                ]
                if rc_ != 0 or not alt:
                    bench.add(
                        section="dbg",
                        mode=mode,
                        dataset=i,
                        k=k,
                        cutoff=c,
                        threads=threads,
                        tool=tool,
                        gate=False,
                        gate_detail=f"exit {rc_}, no output: {err[-300:]}",
                        params=params,
                    )
                    continue
                outp = alt[0]
            fa = os.path.join(d, f"{tool}.unitigs.fa")
            nu = bp = 0
            with open(fa, "w") as o:
                for s in pb.gfa_or_fasta_seqs(outp):
                    o.write(f">{nu}\n{s}\n")
                    nu += 1
                    bp += len(s)
            ok, det, n_tool, n_model = judge_unitigs(
                tool, fa, inputs, k, c, fq, d, threads, kmc_gb, truth, shared
            )
            stats[tool] = fa
            bench.add(
                section="dbg",
                mode=mode,
                dataset=i,
                k=k,
                tool=tool,
                cutoff=c,
                threads=threads,
                wall_s=r.wall_row(),
                rss_mb=round(rss, 1),
                written_gb=None if wrote is None else round(wrote / 1e9, 2),
                read_gb=None if readb is None else round(readb / 1e9, 2),
                peak_disk_gb=round(dp.peak / 1e9, 2),
                unitigs=nu,
                total_bp=bp,
                kmers=n_tool,
                truth_kmers=n_model,
                graph=f"{model}-centric",
                gate=ok,
                gate_detail=det,
                params=params,
            )
        # Maximal unitigs are unique, so on a genome optik's must match Cuttlefish's string for string.
        if (
            "optik" in stats
            and "cuttlefish" in stats
            and not fq
            and not same_unitigs(stats["optik"], stats["cuttlefish"])
        ):
            for r in bench.rows[row0:]:
                if (
                    r.get("dataset"),
                    r.get("k"),
                    r.get("cutoff"),
                    r.get("threads"),
                    r.get("mode"),
                    r.get("tool"),
                ) == (i, k, c, threads, mode, "optik"):
                    r["gate"] = False
                    r["gate_detail"] += (
                        "; unitigs differ from Cuttlefish's (maximal unitigs are unique)"
                    )


# Mapping benchmark


def bench_map(bench, i, dd, mode, work, threads, mach):
    """Run one round of a map dataset at a given mode, k and thread count.

    The reference copy, the simulated reads and a truth mapping are made once in <work>/<id>_inputs and shared by every mode and round, so all mappers see the same reads.
    Each round writes to a fresh directory, which keeps an older round's output from being scored as this one's."""
    kk = dd.get("map_k")
    inp = os.path.join(work, f"{i}_inputs")
    d = os.path.join(work, f"{i}_{mode}_k{kk or 'default'}_t{threads}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    os.makedirs(inp, exist_ok=True)
    src = os.path.join(DATA, i)
    mm2, opt = BIN["minimap2"], BIN["optik"]
    paftools = PAFTOOLS
    bi = dd.get("as_", i)  # a sweep dataset runs its source's protocol
    base_row = {
        "section": "map",
        "mode": mode,
        "dataset": i,
        "k": kk or "default",
        "threads": threads,
    }
    if bi == "bwamem-wgsim-hg38":
        # BWA-MEM Fig. 1, run through the repo's wgsim harness at the paper's parameters with this bench's binaries and memory cap.
        # --fresh makes every mapper map again in this round.
        bwa_prefix = bwa_index_hg38() if os.path.exists(BWA) else None
        base_row["k"] = MM2_PRESET_KW["sr"][
            0
        ]  # optik takes the k and w of minimap2 -ax sr
        rep_f = os.path.join(d, "report.txt")
        with open(rep_f, "w") as rep:
            r = subprocess.run(
                [
                    sys.executable,
                    os.path.join(HERE, "wgsim_check.py"),
                    # Its rows join this bench's map section, so it writes no section of its own.
                    "--report",
                    "",
                    "--ref",
                    dd["local"],
                    "--out",
                    inp,
                    "--pairs",
                    str(dd["pairs"]),
                    "--len",
                    str(dd["length"]),
                    "--error",
                    "0.015",
                    "--mut",
                    "0.002",
                    "--indel-frac",
                    "1.0",
                    "--insert",
                    "500",
                    "--insert-sd",
                    "50",
                    "--window",
                    "20",
                    "--se-both",
                    "--threads",
                    str(threads),
                    "--fresh",
                    "--optik",
                    opt,
                    "--minimap2",
                    mm2,
                    "--bwa",
                    BWA,
                    "--seed-kw",
                    ",".join(map(str, MM2_PRESET_KW["sr"])),
                    "--mem",
                    f"{bench.cap >> 30}G",
                    "--warmup",
                    str(bench.warmup),
                    "--runs",
                    str(bench.runs),
                    *(["--bwa-index", bwa_prefix] if bwa_prefix else []),
                ],
                stdout=rep,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        if r.returncode != 0:
            bench.add(
                **base_row,
                tool="wgsim_check",
                gate=False,
                gate_detail=f"wgsim_check exit {r.returncode}: {r.stderr[-300:]}",
            )
            return
        r1 = next(f for f in os.listdir(inp) if f.endswith(".r1.fq"))
        # Every simulated read is mapped twice, as a pair and again single-end (--se-both).
        # wgsim may write a pair or so fewer than requested.
        with open(os.path.join(inp, r1)) as fh:
            want = 2 * sum(1 for _ in fh) // 4
        with open(rep_f) as fh:
            for line in fh:
                f = WGSIM_FAIL.match(line)
                if f:
                    kind = "cap" if f.group(3) == "CAPPED" else None
                    if kind:
                        add_status(
                            bench,
                            "cap",
                            **base_row,
                            tool=f.group(1),
                            reads=f.group(2),
                            gate_detail=f"memory cap ({bench.cap >> 30} GiB): {f.group(5)}",
                        )
                    else:
                        bench.add(
                            **base_row,
                            tool=f.group(1),
                            reads=f.group(2),
                            gate=False,
                            gate_detail=f"exit {f.group(4)}: {f.group(5)}",
                        )
                    continue
                m = WGSIM_ROW.match(line)
                if not m:
                    continue
                tool, ends = m.group(1), m.group(2)
                sam = [f for f in os.listdir(inp) if f.endswith(f".{tool}.{ends}.sam")]
                got = primaries(os.path.join(inp, sam[0])) if len(sam) == 1 else -1
                ext_ok, ext = (
                    external_sam_checks(os.path.join(inp, sam[0]), dd["local"])
                    if len(sam) == 1
                    else (False, "no SAM")
                )
                # Index-ready time and mapping time come from the last measured run's log, as in the other map rows.
                logf = os.path.join(inp, sam[0]) + ".log" if len(sam) == 1 else ""
                ready = (
                    index_ready_s(tool, Path(logf).read_text())
                    if logf and os.path.exists(logf)
                    else None
                )
                split = (
                    {}
                    if ready is None or not m.group(6)
                    else {
                        "index_ready_s": round(ready, 2),
                        "mapping_s": round(max(0.0, float(m.group(6)) - ready), 2),
                    }
                )
                bench.add(
                    **base_row,
                    tool=tool,
                    reads=ends,
                    wall_s=(
                        {
                            "mean": float(m.group(6)),
                            "sd": None
                            if m.group(8) in (None, "-")
                            else float(m.group(8)),
                            "runs": int(m.group(9)) if m.group(9) else None,
                            "warmup": int(m.group(10)) if m.group(10) else None,
                        }
                        if m.group(6)
                        else None
                    ),
                    rss_mb=int(m.group(7)) if m.group(7) else None,
                    bwa_index_s=bwa_index_s(bwa_prefix)
                    if tool == "bwa" and bwa_prefix
                    else None,
                    mapq30=int(m.group(3)),
                    mapq30_wrong=int(m.group(4)),
                    deciles=m.group(5).strip(),
                    **split,
                    primaries=got,
                    external="pass" if ext_ok else "fail",
                    gate=got == want and ext_ok,
                    params=(
                        "bwa mem on a prebuilt index (bwa_index_s: its build)"
                        if tool == "bwa"
                        else "index built from the FASTA inside the run"
                    )
                    + "; wgsim_check.py, BWA-MEM Fig. 1",
                    gate_detail="; ".join(
                        x
                        for x in (
                            ""
                            if got == want
                            else f"{got} primary records for {want} reads",
                            ext,
                        )
                        if x
                    ),
                )
        return
    tools, want, bwa_prefix, fail = {}, None, None, []
    if bi == "mm2-pbsim2-chm13":
        ref = reference_input(
            os.path.join(DATA, bi, "GCA_009914755.3.fa"), os.path.join(inp, "ref.fa")
        )
        if not os.path.exists(ref + ".fai"):
            subprocess.run([SAMTOOLS or "samtools", "faidx", ref], check=True)
        fq = os.path.join(inp, "reads.fa")  # paftools.js pbsim2fq writes FASTA
        if not os.path.exists(fq):
            # The paper gives no --seed.
            # A fixed one makes a rerun map the same reads.
            subprocess.run(
                [
                    os.path.join(TOOLS, "pbsim2/src/pbsim"),
                    "--seed",
                    "7",
                    "--depth",
                    "1",
                    "--hmm_model",
                    os.path.join(TOOLS, "pbsim2/data/R94.model"),
                    "--length-min",
                    "5000",
                    "--length-mean",
                    "20000",
                    "--accuracy-mean",
                    "0.95",
                    "--prefix",
                    f"{inp}/sd",
                    ref,
                ],
                check=True,
                cwd=inp,
            )
            mafs = sorted(
                os.path.join(inp, f)
                for f in os.listdir(inp)
                if f.startswith("sd_") and f.endswith(".maf")
            )
            with open(fq + ".part", "w") as o:
                subprocess.run(
                    paftools + ["pbsim2fq", ref + ".fai", *mafs], stdout=o, check=True
                )
            os.replace(fq + ".part", fq)
        if dd.get("subset"):
            # Keep every n-th read.
            # pbsim2 writes reads chromosome by chromosome, so the first reads would all come from chr1, while a stride keeps each chromosome's share.
            sub = os.path.join(inp, f"reads_stride{dd['subset']}.fa")
            if not os.path.exists(sub):
                step = max(1, count_records(fq) // dd["subset"])
                with open(sub + ".part", "w") as o:
                    n = 0
                    for j, (nm, sq) in enumerate(_records(fq)):
                        if j % step == 0 and n < dd["subset"]:
                            o.write(f">{nm}\n{sq}\n")
                            n += 1
                os.replace(sub + ".part", sub)
            fq = sub
        want = count_records(fq)
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                fq,
                "-t",
                str(threads),
                "-o",
                f"{d}/optik.sam",
            ],
            "minimap2": [
                mm2,
                "-ax",
                "map-ont",
                "-t",
                str(threads),
                "-o",
                f"{d}/minimap2.sam",
                ref,
                fq,
            ],
        }
        evaluate = "mapeval"
    elif bi == "mm2-mason2":
        ref = reference_input(dd["local"], os.path.join(inp, "ref.fa"))
        r1, r2 = os.path.join(inp, "r1.fq"), os.path.join(inp, "r2.fq")
        if not os.path.exists(r2):
            truth = os.path.join(inp, "mason.sam")
            subprocess.run(
                [
                    tool_path("mason_simulator"),
                    "-ir",
                    ref,
                    "-n",
                    str(dd["pairs"]),
                    "--illumina-read-length",
                    "150",
                    "--illumina-prob-mismatch-scale",
                    "2.5",
                    "-o",
                    os.path.join(inp, "m1.fq"),
                    "-or",
                    os.path.join(inp, "m2.fq"),
                    "-oa",
                    truth,
                    "--seed",
                    "7",
                    "--num-threads",
                    str(threads),
                ],
                check=True,
                cwd=inp,
            )

            # Mates come out as consecutive records, NAME/1 then NAME/2, with both truths in NAME.
            # They are streamed one record at a time.
            # The converted FASTQ of 10 M pairs is 7 GB, and holding it here would take the page cache the timed mappers read from.
            def base(rec):
                """Read name without its /1 or /2 mate suffix."""
                return re.sub(r"/[12]$", "", rec[0].split()[0])

            conv = subprocess.Popen(
                paftools + ["mason2fq", truth], stdout=subprocess.PIPE, text=True
            )
            with open(r1 + ".part", "w") as a, open(r2 + ".part", "w") as b:
                prev = None
                while True:
                    rec = [conv.stdout.readline().rstrip("\n") for _ in range(4)]
                    if not rec[3] and not rec[0]:
                        break
                    if (
                        prev
                        and base(prev) == base(rec)
                        and prev[0].split()[0].endswith("/1")
                    ):
                        a.write("\n".join(prev) + "\n")
                        b.write("\n".join(rec) + "\n")
                        prev = None
                    else:
                        prev = rec
            if conv.wait() != 0:
                raise subprocess.CalledProcessError(
                    conv.returncode, paftools + ["mason2fq", truth]
                )
            os.replace(r1 + ".part", r1)
            os.replace(r2 + ".part", r2)
        with open(r1) as fh:
            pairs = sum(1 for _ in fh) // 4
        if pairs != dd["pairs"]:
            fail.append(f"{pairs} mason2 pairs paired up of {dd['pairs']} simulated")
        want = 2 * pairs
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                r1,
                "-q",
                r2,
                "-t",
                str(threads),
                "-o",
                f"{d}/optik.sam",
            ],
            # -F 800 is minimap2's default and matches the paper's "maximum fragment length to 800 bp".
            "minimap2": [
                mm2,
                "-ax",
                "sr",
                "-t",
                str(threads),
                "-o",
                f"{d}/minimap2.sam",
                ref,
                r1,
                r2,
            ],
        }
        if os.path.exists(BWA):
            bwa_prefix = bwa_index_hg38()
            tools["bwa"] = [
                "bash",
                "-c",
                f"{BWA} mem -t {threads} {bwa_prefix} {r1} {r2} > {d}/bwa.sam",
            ]
        evaluate = "mapeval"
    elif bi == "mm2-winnowmap-truth":
        ref = reference_input(
            os.path.join(src, "GCA_009914755.3.fa"), os.path.join(inp, "ref.fa")
        )
        fq = min(
            os.path.join(src, f) for f in os.listdir(src) if f.endswith(".fastq.gz")
        )
        base = os.path.join(inp, "winnowmap.paf")
        if not os.path.exists(base):
            # Winnowmap's README recipe counts repetitive 15-mers with meryl, then maps HiFi with -x map-pb.
            subprocess.run(
                [tool_path("meryl"), "count", "k=15", "output", f"{inp}/merylDB", ref],
                check=True,
            )
            with open(f"{inp}/repetitive_k15.txt", "w") as o:
                subprocess.run(
                    [
                        tool_path("meryl"),
                        "print",
                        "greater-than",
                        "distinct=0.9998",
                        f"{inp}/merylDB",
                    ],
                    stdout=o,
                    check=True,
                )
            with open(base + ".part", "w") as o:
                subprocess.run(
                    [
                        tool_path("winnowmap"),
                        "-W",
                        f"{inp}/repetitive_k15.txt",
                        "-x",
                        "map-pb",
                        "-t",
                        str(threads),
                        ref,
                        fq,
                    ],
                    stdout=o,
                    check=True,
                )
            os.replace(base + ".part", base)
        want = count_records(fq)
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                fq,
                "-t",
                str(threads),
                "-o",
                f"{d}/optik.sam",
            ],
            # SAM with base-level alignments, as optik writes.
            # PAF without -c would skip minimap2's alignment step and time a different, cheaper product.
            "minimap2": [
                mm2,
                "-ax",
                "map-hifi",
                "-t",
                str(threads),
                "-o",
                f"{d}/minimap2.sam",
                ref,
                fq,
            ],
        }
        evaluate = "pafcmp"
    elif bi == "bwamem-ecoli-asm":
        q = reference_input(
            os.path.join(src, "NC_000913.3.fa"), os.path.join(inp, "query.fa")
        )
        r = ref = reference_input(
            os.path.join(src, "NC_008253.1.fa"), os.path.join(inp, "ref.fa")
        )
        want = 1  # a single query sequence yields one primary record
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                r,
                "-q",
                q,
                "-t",
                str(threads),
                "-o",
                f"{d}/optik.sam",
            ],
            "minimap2": [
                mm2,
                "-ax",
                "asm20",
                "--eqx",
                "-t",
                str(threads),
                "-o",
                f"{d}/minimap2.sam",
                r,
                q,
            ],
        }
        if os.path.exists(BWA):
            # The paper's own run uses BWA-MEM with defaults on the whole MG1655 sequence.
            if not os.path.exists(r + ".bwt"):
                t0 = time.time()
                subprocess.run([BWA, "index", r], check=True, capture_output=True)
                Path(r + ".index_s").write_text(f"{time.time() - t0:.2f}")
            tools["bwa"] = [
                "bash",
                "-c",
                f"{BWA} mem -t {threads} {r} {q} > {d}/bwa.sam",
            ]
            bwa_prefix = r
        evaluate = "call"
    loc = locals()
    warm_set = [
        loc[x] for x in ("ref", "fq", "r1", "r2", "q") if isinstance(loc.get(x), str)
    ] + (
        [bwa_prefix + e for e in (".bwt", ".sa", ".pac", ".ann", ".amb")]
        if bwa_prefix
        else []
    )
    if dd.get("no_bwa") and "bwa" in tools:
        # bwa mem has no comparable k, so in k-scaling rounds its row at this thread count comes from the thread-scaling round.
        del tools["bwa"]
    if kk:
        # Matched seeds give optik and minimap2 the same k and w = 10, with their presets setting the rest.
        # bwa mem keeps its default, since its -k is a minimum exact-match length.
        tools["optik"] = (
            tools["optik"][:2] + ["-k", str(kk), "-w", "10"] + tools["optik"][2:]
        )
        if kk > 28:
            del tools["minimap2"]
            add_status(
                bench,
                "skip",
                **base_row,
                tool="minimap2",
                gate_detail="minimap2 supports k <= 28",
            )
        else:
            tools["minimap2"] = (
                tools["minimap2"][:3]
                + ["-k", str(kk), "-w", "10"]
                + tools["minimap2"][3:]
            )
    elif "minimap2" in tools and "optik" in tools:
        # Without a sweep, optik takes the k and w of the minimap2 preset the paper runs, so both index the same minimizers.
        # optik's own choice of k is not compared here.
        pk, pw = MM2_PRESET_KW[tools["minimap2"][tools["minimap2"].index("-ax") + 1]]
        tools["optik"] = (
            tools["optik"][:2] + ["-k", str(pk), "-w", str(pw)] + tools["optik"][2:]
        )
        base_row["k"] = pk
    for tool, cmd in tools.items():
        tag = f"{mode}/{i}/k{kk or 'default'}/t{threads}/{tool}"
        # Every mapper run starts with identical cached inputs (measure.py).
        # Every run writes its SAM, and the gate below reads the one from the last measured run.
        r = bench.measure(tag, cmd, cwd=d, inputs=warm_set)
        rc_, rss, err = (0 if r.ok else r.code), r.rss() or bench.cap / 2**20, r.err
        if rc_ != 0:
            kind, det = failure(bench, rc_, rss, err)
            if kind == "cap":
                add_status(bench, "cap", **base_row, tool=tool, gate_detail=det)
            else:
                bench.add(**base_row, tool=tool, gate=False, gate_detail=det)
            continue
        # optik and minimap2 build their index from the FASTA inside the run, whereas bwa mem reads a prebuilt one whose build time is reported beside it (bwa_index_s).
        shown = cmd[-1] if cmd[:2] == ["bash", "-c"] else " ".join(map(str, cmd))
        row = dict(
            **base_row,
            tool=tool,
            wall_s=r.wall_row(),
            rss_mb=round(rss, 1),
            params=shown
            + (
                "  [prebuilt bwa index; its build in bwa_index_s]"
                if tool == "bwa"
                else "  [index built from the FASTA inside this run]"
            ),
        )
        if tool == "bwa":
            row["bwa_index_s"] = bwa_index_s(bwa_prefix)
        # optik and minimap2 index the reference inside every run and log when the index was ready, counted from process start.
        # The remainder of the run is mapping alone, so the three mappers can be compared with and without the index build.
        ready = index_ready_s(tool, err)
        if ready is not None and r.ok:
            row["index_ready_s"] = round(ready, 2)
            row["mapping_s"] = round(max(0.0, r.wall["mean"] - ready), 2)
        sam_out = f"{d}/{tool}.sam"
        if not os.path.exists(sam_out) or os.path.getsize(sam_out) == 0:
            bench.add(**row, gate=False, gate_detail="exit 0 but no SAM written")
            continue
        ext_ok, ext = external_sam_checks(sam_out, ref)
        n_prim = primaries(sam_out)
        row.update(
            external="pass" if ext_ok else "fail",
            primaries=n_prim,
            seq_mb=round(seq_bytes(sam_out) / 1e6, 2),
        )
        det = [
            x
            for x in fail
            + [
                ext,
                "" if n_prim == want else f"{n_prim} primary records for {want} reads",
            ]
            if x
        ]
        ok = ext_ok and n_prim == want and not fail
        try:
            if evaluate == "mapeval":
                ev = subprocess.run(
                    paftools + ["mapeval", sam_out],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                Path(f"{d}/{tool}.mapeval").write_text(ev)
                # Rows read "Q mapq n_in_bin n_wrong_in_bin cum_err_rate cum_mapped" in descending mapq.
                # The papers plot the cumulative error rate at MAPQ >= 10.
                q10 = [
                    x.split()
                    for x in ev.splitlines()
                    if x.startswith("Q\t") and int(x.split()[1]) >= 10
                ]
                r10 = q10[-1] if q10 else None
                if r10 is None:
                    det.append("mapeval produced no MAPQ >= 10 row")
                bench.add(
                    **row,
                    mapq10_mapped=int(r10[5]) if r10 else None,
                    mapq10_err_rate=float(r10[4]) if r10 else None,
                    calibrated=mapeval_calibration(ev),
                    gate=ok and r10 is not None,
                    gate_detail="; ".join(det),
                )
            elif evaluate == "pafcmp":
                # Both mappers' SAM go through the same conversion.
                test = f"{d}/{tool}.paf"
                with open(test, "w") as o:
                    subprocess.run(
                        paftools + ["sam2paf", sam_out], stdout=o, check=True
                    )
                cmp_ = subprocess.run(
                    paftools + ["pafcmp", base, test],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                Path(f"{d}/{tool}.pafcmp").write_text(cmp_)
                if not cmp_.strip():
                    det.append("pafcmp produced nothing")
                bench.add(
                    **row,
                    pafcmp=" | ".join(
                        x.strip() for x in cmp_.strip().splitlines()[-4:]
                    ),
                    gate=ok and bool(cmp_.strip()),
                    gate_detail="; ".join(det),
                )
            else:
                # paftools.js call, minimap2's assembly-variant caller, counts substitutions over the bases one alignment covers.
                # Picard's check above covers CIGAR/SEQ consistency.
                call_snv, cov1 = paftools_call_snvs(sam_out, ref, d, tool)
                if not cov1:
                    det.append("paftools.js call: no bases covered by one alignment")
                bench.add(
                    **row,
                    substitutions=call_snv,
                    aligned_bases_cov1=cov1,
                    gate=ok and cov1 > 0,
                    gate_detail="; ".join(det),
                )
        except subprocess.CalledProcessError as e:
            bench.add(
                **row,
                gate=False,
                gate_detail="; ".join(
                    det
                    + [
                        f"{evaluate} failed (exit {e.returncode}): {(e.stderr or '')[-200:]}"
                    ]
                ),
            )


# Pairwise alignment benchmark


def bench_align(bench, i, dd, mode, work):
    """Time optik align and BiWFA on the pairs of dataset `i` and check that their scores agree.

    biwfa-ont-ul48 unpacks the BiWFA paper's 48 ultra-long ONT pairs.
    The other sets are made with WFA2-lib's generate_dataset over the dataset's (length, error) grid."""
    d = os.path.join(work, f"{i}")
    os.makedirs(d, exist_ok=True)
    exe = pb.build_wfab(d)
    sets = []
    if i == "biwfa-ont-ul48":
        tsv = os.path.join(d, "pairs.tsv")
        if not os.path.exists(tsv):
            with (
                tarfile.open(
                    os.path.join(DATA, i, "ONT_MinION_UL.500kbps.48.tar.gz")
                ) as tf,
                open(tsv + ".part", "w") as o,
            ):
                for m in sorted(tf.getmembers(), key=lambda x: x.name):
                    if m.isfile() and m.name.endswith(".seq"):
                        lines = tf.extractfile(m).read().decode().split("\n")
                        p = (
                            next(x for x in lines if x.startswith(">"))[1:]
                            .strip()
                            .upper()
                        )
                        t = (
                            next(x for x in lines if x.startswith("<"))[1:]
                            .strip()
                            .upper()
                        )
                        o.write(f"{p}\t{t}\n")
            os.replace(tsv + ".part", tsv)
        sets.append(("ont_ul_500k_48", tsv, 48))
    else:
        gen = os.path.join(T, "WFA2-lib/bin/generate_dataset")
        for read_len, e in dd["grid"]:
            tsv = os.path.join(d, f"sim_l{read_len}_e{e}.tsv")
            if not os.path.exists(tsv):
                seq = os.path.join(d, f"sim_l{read_len}_e{e}.seq")
                subprocess.run(
                    [gen, "-o", seq, "-n", "1", "-l", str(read_len), "-e", str(e)],
                    check=True,
                )
                lines = Path(seq).read_text().split("\n")
                with open(tsv + ".part", "w") as o:
                    for a, b in zip(lines[0::2], lines[1::2], strict=False):
                        if a.startswith(">") and b.startswith("<"):
                            o.write(f"{a[1:]}\t{b[1:]}\n")
                os.replace(tsv + ".part", tsv)
            sets.append((f"sim_l{read_len}_e{e}", tsv, 1))
    for name, tsv, want in sets:
        with open(tsv) as fh:
            npairs = sum(1 for _ in fh)
        scores = {}
        # Only BiWFA runs here.
        # These pairs of 100 kbp to 2 Mbp are the BiWFA paper's long-sequence case, where full-memory WFA needs memory quadratic in the score, about 1 TB at 1 Mbp and 10% error, and cannot run on any machine.
        # Full-memory WFA is compared where it is a real choice, on the pairs up to 100 kbp in compare_tools_small_sets.py and on optik map's gaps of 1 kb or less.
        for m in ["optik-biwfa", "wfa2-ultralow"]:
            out = os.path.join(d, f"scores_{m}_{name}.txt")
            res = os.path.join(d, "res.txt")
            # The harness times its own alignments, which is BiWFA's measure and leaves out parsing.
            # Warmups are dropped, and the statistics come from the --runs.
            runs_ = bench.repeat(
                f"{mode}/{i}/{name}/{m}", [exe, m, tsv, out], stdout=res
            )
            rc_, rss, e = runs_[-1] if runs_ else (-1, 0.0, "no run")
            rss = max(x[1] for x in runs_) if runs_ else 0.0
            if rc_ != 0:
                kind, det = failure(bench, rc_, rss, e)
                if kind == "cap":
                    add_status(
                        bench,
                        "cap",
                        section="align",
                        mode=mode,
                        dataset=i,
                        set=name,
                        method=m,
                        gate_detail=det,
                    )
                else:
                    bench.add(
                        section="align",
                        mode=mode,
                        dataset=i,
                        set=name,
                        method=m,
                        gate=False,
                        gate_detail=det,
                    )
                continue
            # One "align_seconds S invalid N pairs P" line per kept run
            lines = [x[2].split() for x in runs_]
            secs_all = [float(x[1]) for x in lines]
            invalid, aligned = (
                max(int(x[3]) for x in lines),
                min(int(x[5]) for x in lines),
            )
            secs = sum(secs_all) / len(secs_all)
            secs_sd = (
                (sum((x - secs) ** 2 for x in secs_all) / (len(secs_all) - 1)) ** 0.5
                if len(secs_all) > 1
                else None
            )
            scores[m] = Path(out).read_text().split()
            det = [
                x
                for x in (
                    "" if invalid == 0 else f"{invalid} CIGARs do not re-score",
                    ""
                    if npairs == want == aligned == len(scores[m])
                    else f"{aligned} pairs aligned, {len(scores[m])} scores, {npairs} in the set, "
                    f"{want} expected",
                )
                if x
            ]
            bench.add(
                section="align",
                mode=mode,
                dataset=i,
                set=name,
                method=m,
                pairs=npairs,
                align_s=round(secs, 3),
                align_s_sd=None if secs_sd is None else round(secs_sd, 3),
                runs=len(secs_all),
                warmup=bench.warmup,
                rss_mb=round(rss, 1),
                gate=not det,
                gate_detail="; ".join(det),
            )
        # The two exact methods must agree pair for pair.
        # If either is missing there is nothing to cross-check, and that fails as well.
        ref = next(iter(scores.values()), None)
        why = (
            "scores differ across exact methods"
            if any(v != ref for v in scores.values())
            else ""
            if len(scores) == 2
            else "only one method ran: scores not cross-checked"
        )
        if why:
            for r in bench.rows:
                if (
                    r.get("dataset") == i
                    and r.get("set") == name
                    and not r.get("status")
                ):
                    r["gate"] = False
                    r["gate_detail"] = (
                        (r.get("gate_detail") or "") + "; " + why
                    ).lstrip("; ")


# Smoke test

SMOKE_BASES = 2_000_000


def fetch_prefix(i, d, nbytes=24 << 20):
    """Fetch the first nbytes of each file not yet fetched, through HTTP Range.

    That is enough for a smoke slice, and a compressed prefix decompresses up to its cut."""
    src = os.path.join(DATA, i)
    os.makedirs(src, exist_ok=True)
    if d.get("kind") == "reads" or d.get("reads"):
        urls = file_urls(dict(d, files=d["files"][:1]))[
            :2
        ]  # two read files are enough to show the format
    else:
        urls = file_urls(d)
    for url, name, _, _, _ in urls:
        dest = os.path.join(src, name)
        if os.path.exists(dest) or os.path.exists(dest + ".prefix"):
            continue
        if name.endswith((".zip", ".tar.gz")):
            continue  # an archive is only readable whole
        subprocess.run(
            ["curl", "-fsL", "-r", f"0-{nbytes - 1}", "-o", dest + ".prefix", url],
            check=True,
        )
    for f in os.listdir(src):  # prefixes keep the real names, in a side directory
        if f.endswith(".prefix"):
            os.makedirs(os.path.join(src, "_prefix"), exist_ok=True)
            os.replace(os.path.join(src, f), os.path.join(src, "_prefix", f[:-7]))


def smoke(args, mach):
    """Run every tool of the dataset's section on the first 2 Mbp or so of the formatted inputs.

    The formats are right when every tool runs and the gates pass on that slice.
    Nothing is written to the results directory.
    Returns 1 when any gate fails, else 0."""
    global DATA
    out = os.path.join(T, "published_sets_smoke")
    os.makedirs(out, exist_ok=True)
    bench = pb.Bench(
        out, 1, 0, max(1, int(mach["mem_gib"] - 2)), os.path.join(out, "_scratch")
    )
    only = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if args.skip else set()
    nfail = 0
    checks = []  # (dataset, check, passed, detail) of the query rounds
    for i, d in D.items():
        if (
            "excluded" in d
            or ("opt_in" in d and not only)
            or (only and i not in only)
            or i in skip
            or (args.section and d["section"] not in args.section.split(","))
        ):
            continue
        if d["section"] not in ("query", "dbg") or d.get("local"):
            continue
        # A scaling set with no files of its own runs on its source set's data (`as_`).
        src_id = d["as_"] if d.get("as_") and not d["files"] else i
        src = os.path.join(DATA, src_id)
        fetched = os.path.exists(os.path.join(src, "PROVENANCE.json"))
        if not fetched:
            if args.no_download:
                continue
            if any(
                f.get("key", f.get("url", "")).endswith((".tar.gz", ".zip"))
                or "ncbi" in f
                for f in D[src_id]["files"]
            ):
                print(
                    f"\n=== smoke {d['section']} {i}: an archive, which needs the full download "
                    f"(fetch --only {i}, then smoke)"
                )
                continue
            try:
                fetch_prefix(src_id, D[src_id])
            # A source outage skips this set, not the whole run.
            except FETCH_ERRORS as e:
                print(
                    f"\n=== smoke {d['section']} {i}: prefix download failed ({str(e)[:200]})"
                )
                continue
        saved = DATA
        if not fetched:
            # The prefixes stand in for the files under the files' own names.
            DATA = os.path.join(T, "published_sets_smoke_data")
            os.makedirs(os.path.join(DATA, src_id), exist_ok=True)
            for f in (
                os.listdir(os.path.join(src, "_prefix"))
                if os.path.isdir(os.path.join(src, "_prefix"))
                else []
            ):
                lnk = os.path.join(DATA, src_id, f)
                if not os.path.exists(lnk):
                    os.symlink(os.path.join(src, "_prefix", f), lnk)
        print(
            f"\n=== smoke {d['section']} {i} ({'fetched' if fetched else 'prefix'})",
            flush=True,
        )
        try:
            if d["section"] == "query":
                k = d["k"][-1]
                unitigs_fa, by, wl, windows = query_inputs(
                    i,
                    d,
                    k,
                    "full",
                    os.path.join(out, "work"),
                    args.threads,
                    limit_bases=SMOKE_BASES,
                )
                wl.pop("_load")
                dq = os.path.join(out, "work", f"{i}_k{k}_full_smoke", "optimal")
                os.makedirs(dq, exist_ok=True)
                tools = query_tools(
                    "optimal",
                    dq,
                    unitigs_fa,
                    k,
                    sshash_m(d, k, windows),
                    args.threads,
                    mach,
                    bench.cap >> 30,
                )
                counts = {}
                for tool, t in tools.items():
                    if (
                        t["build"] is not None
                    ):  # None marks a query mode of another tool's index
                        rc_, _wall, _rss, err = bench.once(
                            f"smoke/{i}/{tool}/build", t["build"], cwd=dq
                        )
                        if rc_ != 0:
                            print(f"  [FAIL] {tool} build: {err.strip()[-300:]}")
                            checks.append(
                                (i, f"{tool} build", False, err.strip()[-200:])
                            )
                            nfail += 1
                            continue
                    for w, q in wl.items():
                        cmd, outf = t["query"](q)
                        capt = (
                            os.path.join(dq, "q.out")
                            if tool in OWN_OUTPUT
                            else (outf or os.path.join(dq, "q.out"))
                        )
                        rc_, _wall, _rss, err = bench.once(
                            f"smoke/{i}/{tool}/{w}",
                            cmd,
                            cwd=dq,
                            stdout=None if tool == "sbwt" else capt,
                        )
                        if rc_ != 0:
                            print(f"  [FAIL] {tool} {w}: {err.strip()[-300:]}")
                            checks.append((i, f"{tool} {w}", False, err.strip()[-200:]))
                            nfail += 1
                            continue
                        counts.setdefault(w, {})[tool] = count_of(
                            tool, outf, capt, err
                        )[1]
                for w, c in counts.items():
                    same = len(set(c.values())) == 1
                    print(
                        f"  [{'PASS' if same else 'FAIL'}] {w}: positives {c}  (SPSS: {by}, "
                        f"{windows} k-mer windows)"
                    )
                    checks.append((i, f"{w} positives", same, str(c)))
                    nfail += not same
            else:
                bench_dbg(
                    bench,
                    i,
                    d,
                    "optimal",
                    os.path.join(out, "work"),
                    args.threads,
                    mach,
                    limit_bases=SMOKE_BASES,
                )
        except Exception as e:
            # Any failure of one dataset, a bench bug included, is reported with its traceback and the smoke run goes on.
            LOG.exception("  [FAIL] %s", i)
            checks.append((i, "run", False, f"{type(e).__name__}: {str(e)[-200:]}"))
            nfail += 1
        finally:
            DATA = saved
    # The graph rounds record their gates as rows.
    # SKIP and CAP rows have no gate.
    nfail += sum(1 for r in bench.rows if r.get("gate") is False)
    print(f"\nsmoke: {nfail} gate failures")
    hidden = ("section", "params", "gate_detail", "status", "gate")
    results_file.write_section(
        args.report,
        "published_smoke",
        "Functional check of every published dataset on a prefix of its input: each tool builds "
        "and queries, query tools agree on positives, and graph outputs pass their gates",
        [
            (
                "query",
                [
                    {
                        "dataset": i,
                        "check": c,
                        "gate": "PASS" if ok else "FAIL",
                        "detail": det,
                    }
                    for i, c, ok, det in checks
                ],
            ),
            (
                "dbg",
                [
                    {
                        **{k: v for k, v in r.items() if k not in hidden},
                        "gate": gate_text(r).strip(),
                    }
                    for r in bench.rows
                ],
            ),
        ],
        notes=[f"Logs: `{os.path.relpath(out, ROOT)}`"],
    )
    return 1 if nfail else 0


# Report


# Headings of the paper sections in the sectioned results file
SECTION_TITLES = {
    "query": "k-mer membership queries against SSHash, SBWT and FMSI on their papers' datasets",
    "dbg": "compacted de Bruijn graphs against Cuttlefish and GGCAT on their papers' datasets",
    "map": "read mapping against minimap2, BWA-MEM and Winnowmap on their papers' simulations",
    "align": "pairwise alignment against WFA2-lib BiWFA on its paper's sequence pairs",
}


def gate_text(r):
    """PASS, FAIL or the status (SKIP, CAP) of one row, with its detail when it did not pass."""
    if r.get("status"):
        return r["status"].upper() + " " + str(r.get("gate_detail", ""))
    return "PASS" if r.get("gate") else "FAIL " + str(r.get("gate_detail", ""))


def report(bench, out, mach, skipped, dest=""):
    """Write results.json and SUMMARY.md of a bench run to `out`, and one section per paper section of the sectioned results file `dest`."""
    Path(os.path.join(out, "results.json")).write_text(
        json.dumps(
            {"machine": mach, "skipped": skipped, "rows": bench.rows},
            indent=1,
            default=str,
        )
    )
    md_lines = [
        "# optik vs published tools, on the tools' own published datasets",
        "",
        (
            "bench/compare_tools_published_sets.py; provenance of each dataset in its PROVENANCE.json under "
            f"{os.path.relpath(DATA, ROOT)}/."
        ),
        "",
        (
            f"Machine: {mach['mem_gib']} GiB usable memory (per-process cap {bench.cap >> 30} GiB), "
            f"{mach['disk_free_gb']} GB free disk. PASS / FAIL are correctness gates; SKIP is a tool whose "
            "reported or estimated memory need exceeds this machine; CAP is a run that died against the cap."
        ),
        "",
    ]
    for sec in ("query", "dbg", "map", "align"):
        rows = [r for r in bench.rows if r.get("section") == sec]
        if not rows:
            continue
        keys = []
        for r in rows:
            keys += [
                k
                for k in r
                if k not in keys
                and k not in ("section", "params", "gate_detail", "status")
            ]
        md_lines += [
            f"## {sec}",
            "",
            "| " + " | ".join(keys) + " | gate |",
            "|" + "---|" * (len(keys) + 1),
        ]
        for r in rows:
            st = gate_text(r)
            md_lines.append(
                "| "
                + " | ".join(pb._fmt(r.get(k, "")) for k in keys)
                + " | "
                + st
                + " |"
            )
        md_lines.append("")
    if skipped:
        md_lines += ["## Not run on this machine", ""]
        md_lines += [
            f"- `{i}` ({D[i].get('title')}): {why}" for i, why in skipped.items()
        ]
        md_lines.append("")
    excl = [(i, d) for i, d in D.items() if "excluded" in d]
    md_lines += ["## Not runnable anywhere (from the papers, with the reason)", ""]
    md_lines += [
        f"- `{i}` ({d.get('title')}; {d['paper']}): {d['excluded']}" for i, d in excl
    ]
    Path(os.path.join(out, "SUMMARY.md")).write_text("\n".join(md_lines) + "\n")
    hidden = ("section", "params", "gate_detail", "status", "gate")
    for sec, title in SECTION_TITLES.items():
        rows = [r for r in bench.rows if r.get("section") == sec]
        if not rows:
            continue
        flat = [
            {
                **{k: v for k, v in r.items() if k not in hidden},
                "gate": gate_text(r).strip(),
            }
            for r in rows
        ]
        notes = [
            f"- `{i}` not run here: {why}"
            for i, why in skipped.items()
            if D[i].get("section") == sec
        ]
        results_file.write_section(
            dest,
            f"published_{sec}",
            title + f" (per-process cap {bench.cap >> 30} GiB)",
            [("", flat)],
            notes=notes
            + [f"Raw rows, parameters and logs: `{os.path.relpath(out, ROOT)}`"],
        )


def pid_alive(pid):
    """True when process `pid` exists, so a fetch killed seconds ago is not shown as running.

    A progress file without a pid, written before pids were recorded, counts as alive."""
    if pid is None:
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # it exists, under another user
    return True


def part_done(part):
    """Bytes of an in-progress download that are already on disk.

    These are the completed ranges the range scheduler recorded in
    <part>.done or, failing that, the size of a plain curl .part, which curl
    writes from the start."""
    st = part + ".done"
    if os.path.exists(st):
        try:
            return sum(e - s for s, e in json.loads(Path(st).read_text()))
        except (OSError, ValueError):
            pass
    return os.path.getsize(part) if os.path.exists(part) else 0


def cmd_status(args):
    """Print every dataset's download state.

    For each dataset this shows the bytes its sources publish (from the size cache that fetch fills) and the bytes on disk, counting finished files and the written part of files in progress.
    It also shows whether the dataset is verified (PROVENANCE.json appears once all its files have passed their checksums) and the fetch's measured rate and ETA, taken from DATA/_fetch_progress.json, which the range scheduler rewrites every few seconds during a fetch."""
    while True:
        sizes = json.loads(Path(_SIZES).read_text()) if os.path.exists(_SIZES) else {}
        prog = os.path.join(DATA, "_fetch_progress.json")
        pg = None
        if os.path.exists(prog):
            try:
                pg = json.loads(Path(prog).read_text())
            except (OSError, ValueError):
                pg = None
        live = (
            pg is not None
            and time.time() - pg.get("time", 0) < 60
            and pid_alive(pg.get("pid"))
        )
        rows, left_all = [], 0.0
        for i, d in D.items():
            if "excluded" in d or not d.get("files"):
                continue
            want = sizes.get(i)
            out = os.path.join(DATA, i)
            have = 0
            parts = 0
            if os.path.isdir(out):
                for f in os.listdir(out):
                    pth = os.path.join(out, f)
                    if f.endswith(".part"):
                        have += part_done(pth)
                        parts += 1
                    elif (
                        os.path.isfile(pth)
                        and not f.endswith((".done", ".json", ".tmp"))
                        and not f.startswith("_")
                    ):
                        have += os.path.getsize(pth)
            verified = os.path.exists(os.path.join(out, "PROVENANCE.json"))
            if verified:
                # Count the downloaded files themselves, not what was unpacked or derived from them.
                recs = json.loads(
                    Path(os.path.join(out, "PROVENANCE.json")).read_text()
                ).get("files", [])
                have = sum(
                    os.path.getsize(r["path"])
                    for r in recs
                    if not r.get("local")
                    and r.get("path")
                    and os.path.exists(r["path"])
                )
                state = "verified"
            elif parts:
                state = "downloading" if live else "partial (fetch not running)"
            elif have:
                state = "partial"
            else:
                state = "-"
            pct = (
                "100%"
                if verified
                else (f"{min(100.0, 100 * have / (want * 1e9)):.1f}%" if want else "?")
            )
            if "opt_in" in d and not verified:
                state = "opt-in, not fetched by default" + (
                    " (partial)" if have else ""
                )
            elif want and not verified:
                left_all += max(0.0, want * 1e9 - have)
            rows.append((i, d["section"], want, have / 1e9, pct, state))
        print(
            f"{'dataset':28s} {'section':7s} {'source GB':>9s} {'on disk GB':>10s} {'done':>7s}  state"
        )
        for i, sec, want, have, pct, state in rows:
            print(
                f"{i:28s} {sec:7s} {('?' if want is None else f'{want:.2f}'):>9s} {have:10.2f} {pct:>7s}  {state}"
            )
        if live:
            rate = pg.get("rate_bytes_s", 0.0)
            infl = sum(f["left"] for f in pg.get("files", []))
            eta = (
                f"{infl / rate / 60:.0f} min for the files in flight"
                if rate > 0
                else "rate not measured yet"
            )
            print(
                f"\nfetch running: {rate / 1e6:.1f} MB/s over {pg.get('connections')} connections "
                f"(target {pg.get('target')}); {infl / 1e9:.1f} GB left in flight; ETA {eta}"
            )
        else:
            print(
                "\nno fetch running"
                + (f" (last progress {time.ctime(pg['time'])})" if pg else "")
            )
        print(
            f"not yet fetched and verified: {left_all / 1e9:.1f} GB over all listed datasets "
            "(fetch takes the ones this machine has disk for: `fetch --tier auto`)"
        )
        if not args.watch:
            return
        time.sleep(args.watch)
        print()


# What a section measures, following the latest paper of the compared tools, and how they run.
# A tool's row names the regime its authors built it for, so it is compared where it is meant to be strong.
PROTOCOL = {
    "query": {
        "papers": "SSHash 2026 (btag264): bits/k-mer, build time and memory, per-k-mer time of positive "
        "lookups (half reverse complements), negative lookups (random k-mers) and streaming lookups "
        "over a read set, queries on one thread; FMSI 2025: index size across k (15/23/31 on the big "
        "sets) and query time; SBWT 2023: per-k-mer search over the same sets",
        "metrics": "build_s (mean±sd), build_rss_mb, index_mb, bits_per_kmer, load_s, query_s, "
        "ns_per_kmer_net (load subtracted), positives against KMC",
        "tools": "SSHash: its paper's m (FMSI's rule on FMSI sets, else k/2), -g budget, single-threaded "
        "queries -- built for streaming lookups over reads; SBWT: plain-matrix with reverse "
        "complements, -m 16 (its authors' benchmark build) -- built for streaming search; FMSI: "
        "KmerCamel -S + fmsi index, query -S -O (its paper) -- built for the smallest index; "
        "optik-spss: default m chosen from a sample (and SSHash's m in parity rounds); optik-fm: "
        "k-order BWT over the unitigs",
    },
    "dbg": {
        "papers": "Cuttlefish 2 (Genome Biol 2022) and Cuttlefish 3 (2025): wall time, peak RSS, intermediate "
        "disk (bytes written, peak), thread scaling (HG004 20x, k = 27) and k scaling",
        "metrics": "wall_s (mean±sd), rss_mb, written_gb, read_gb, peak_disk_gb, unitigs, k-mers against KMC",
        "tools": "Cuttlefish 3: --ref on genomes, --read -c f0 on reads, -m half this machine's memory -- "
        "built for read sets at a small memory budget; Bifrost: -r / -s (cutoff 2 only) where its "
        "memory fits -- built for dynamic, colored graphs (not measured here); optik dbruijn: its "
        "default memory budget, and --max-mem = Cuttlefish's -m where a budget binds",
    },
    "map": {
        "papers": "minimap2 (2018) mapeval: error rate against mapped fraction by MAPQ; BWA-MEM (2013): "
        "wgsim reads, wrong vs MAPQ; Winnowmap2 (2022): agreement with its mappings (pafcmp)",
        "metrics": "wall_s (mean±sd), rss_mb, index_ready_s / mapping_s (optik, minimap2), bwa_index_s, "
        "mapq10_err_rate, mapq30 wrong, MAPQ calibration, Picard and calmd on every SAM",
        "tools": "minimap2: the preset of the read type (sr, map-ont, map-hifi, asm20); bwa mem: defaults on "
        "a prebuilt index -- built for short-read accuracy; optik: the k and w of the minimap2 "
        "preset (sweeps: the same k and w = 10 for both)",
    },
    "align": {
        "papers": "BiWFA (Bioinformatics 2023): time and memory per alignment of long pairs, one thread",
        "metrics": "align_s (mean±sd over the harness's own timer), rss_mb, identical optimal scores",
        "tools": "WFA2-lib BiWFA (MemoryUltralow), the paper's implementation; optik's BiWFA; full-memory WFA "
        "only on compare_tools_small_sets.py's shorter pairs, where it can run",
    },
}


def cmd_list(args):
    """Print the dataset catalogue by section, with each dataset's protocol and whether it runs here."""
    mach = machine()
    _, src = bytes_per_kmer()
    print(
        f"this machine: {mach['mem_gib']} GiB memory, {mach['disk_free_gb']} GB free disk; "
        f"memory estimates from {src}\n"
    )
    if args.md:
        print(
            "| id | what | section | paper | source | use | here |\n|---|---|---|---|---|---|---|"
        )
    for sec in ("query", "dbg", "map", "align", "count"):
        if not args.md:
            print(f"== {sec}")
            pr = PROTOCOL.get(sec)
            if pr:
                print(
                    f"  protocol: {pr['papers']}\n  metrics:  {pr['metrics']}\n  tools:    {pr['tools']}\n"
                )
        for i, d in D.items():
            if d["section"] != sec:
                continue
            use = d.get("excluded") and "EXCLUDED: " + d["excluded"] or d.get("use", "")
            if d.get("note"):
                use += " NOTE: " + d["note"]
            if "excluded" in d:
                here = "excluded"
            elif "opt_in" in d:
                here = d["opt_in"]
            else:
                ok, why = preflight(i, d, mach)
                here = "runnable" if ok else "; ".join(why)
                gb = download_gb(i, d)
                if gb:
                    here += f" ({gb:.1f} GB download)"
                kmers = d.get("kmers31")
                tl = [
                    f"{t}: {w}"
                    for t in (
                        QUERY_TOOLS
                        if sec == "query"
                        else DBG_TOOLS
                        if sec == "dbg"
                        else ()
                    )
                    if (
                        w := tool_skip(
                            d,
                            t,
                            mach,
                            31 if sec == "query" else d.get("k", [31])[0],
                            kmers,
                        )
                    )
                ]
                if tl:
                    here += "; skips " + "; ".join(tl)
            if args.md:
                print(
                    f"| `{i}` | {d.get('title')} | {sec} | {d['paper']} | {d['source']} | {use} | {here} |"
                )
            else:
                print(
                    f"  {i:24} {d.get('title')}\n  {'':24} {d['paper']}\n  {'':24} {d['source']}\n"
                    f"  {'':24} {use}\n  {'':24} here: {here}\n"
                )


def sweep_ks(i, d, k_list):
    """The k values a dataset runs at.

    Without --k-list this is the paper's k for the dataset.
    With it, the listed k apply wherever the dataset can be rebuilt at them, which excludes the SSHash 2026 deposits, whose eulertigs exist only at k = 31, and FMSI's 10% subsample at k = 23.
    For map datasets it is the seed k of the as_ sweep sets, and None leaves each mapper at its own default for the paper datasets."""
    if d["section"] == "map":
        # Only the matched-k sweep sets.
        # A paper's own map dataset runs its preset once.
        return (k_list or d.get("k")) if d.get("as_") else None
    if d["section"] not in ("query", "dbg"):
        return None
    if not k_list:
        return d["k"]
    if i.startswith("sshash26"):
        return [k for k in k_list if k in d["k"]] or d["k"]
    return k_list


# Matched thread count of the k-scaling plan.
# Cuttlefish 2 Table 4 uses 8 threads, and a machine with fewer usable CPUs runs at all of them, which the round's label records.
PAPER_PARITY_T = 8
PARITY_T = min(PAPER_PARITY_T, hw.cpus())


def plan_rounds(i, d, base_t):
    """List the rounds --plan paper runs for one dataset.

    Each round is (mode, k list, threads, first_t, sweep_t, label, extra dataset fields).
    The papers' rows run each tool at its optimal settings on all threads at the paper's k (queries add parity rows, shared where the commands coincide).
    k scaling (d["kscale"]) repeats the tools over the k list with PARITY_T threads and parity settings.
    Thread scaling covers hw.thread_points(d["tscale"], base_t) at the paper's k, or at the largest k of a sweep set.
    Sweep and scaling sets have no paper round of their own."""
    sec = d["section"]
    rounds = []
    ks = d.get("k") if sec in ("query", "dbg") else None
    if not i.startswith(("sweep-", "scale-", "kscale-")):
        for m in ["parity", "optimal"] if sec == "query" else ["optimal"]:
            rounds.append(
                (m, ks, 1 if sec == "align" else base_t, True, False, "paper")
            )
    if d.get("kscale"):
        # The streaming query, a whole read run on one thread for every tool, happens only at the paper's k, the largest of the sweep.
        # Rounds at other k measure builds, sizes and lookups.
        rounds.append(
            (
                "parity",
                d["kscale"],
                PARITY_T,
                True,
                False,
                "kscale",
                {"stream_k": [max(d["kscale"])]} if sec == "query" else {},
            )
        )
    if d.get("tscale"):
        k0 = [max(ks)] if ks else None
        # Every query paper measures queries on one thread, so thread scaling re-measures only the builds.
        # If the k-scaling round already ran at k0, it holds every tool's queries and the builds of the tools without a thread option.
        # These rounds then rebuild the multi-threaded builders and rerun only optik's queries at t > 1, an extra row for optik's own multi-threaded query.
        covered = sec == "query" and d.get("kscale") and k0 and k0[0] in d["kscale"]
        for j, t in enumerate(hw.thread_points(d["tscale"], base_t)):
            rounds.append(
                ("parity", k0, t, j == 0 and not covered, True, "tscale")
                + (({"no_queries": True} if t == 1 else {}) if covered else {},)
            )
    return [r if len(r) == 7 else r + ({},) for r in rounds]


def main():
    """Dispatch to list, tools, status, fetch, bench or smoke from the command line."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list")
    p.add_argument("--md", action="store_true")
    sub.add_parser("tools")
    p = sub.add_parser(
        "status", help="download progress of every dataset (fetch running or not)"
    )
    p.add_argument("--watch", type=float, default=0, help="refresh every N seconds")
    for name in ("fetch", "bench", "smoke"):
        p = sub.add_parser(name)
        p.add_argument(
            "--tier",
            default="auto",
            choices=["auto", "S", "M", "L"],
            help="auto (default): every dataset this machine has the disk and tools for",
        )
        p.add_argument(
            "--force", action="store_true", help="ignore the disk/tool preflight"
        )
        p.add_argument("--section", default=None, help="query,dbg,map,align")
        p.add_argument(
            "--only", default=None, help="comma-separated dataset ids (ignores --tier)"
        )
        p.add_argument(
            "--skip", default=None, help="comma-separated dataset ids to leave out"
        )
        p.add_argument(
            "--threads",
            type=int,
            default=hw.cpus(),
            help="threads of every parallel tool (default: the usable CPUs)",
        )
        if name in ("bench", "smoke"):
            results_file.add_argument(p)
        if name == "smoke":
            p.add_argument(
                "--no-download",
                action="store_true",
                help="only datasets already fetched",
            )
        if name == "fetch":
            p.add_argument(
                "--jobs", type=int, default=6, help="datasets downloaded at once"
            )
        if name == "bench":
            p.add_argument(
                "--mode", default="both", choices=["parity", "optimal", "both"]
            )
            p.add_argument(
                "--runs",
                type=int,
                default=3,
                help="measured runs of every timed command (statistics over these)",
            )
            p.add_argument(
                "--warmup",
                type=int,
                default=1,
                help="warmup runs of every timed command, discarded",
            )
            p.add_argument(
                "--cap-gb",
                type=int,
                default=0,
                help="per-process memory cap (default: usable memory - 2 GiB)",
            )
            p.add_argument("--work", default=os.path.join(T, "published_sets_work"))
            p.add_argument(
                "--keep-work",
                action="store_true",
                help="keep each dataset's working files (default: deleted once its rows are recorded)",
            )
            p.add_argument(
                "--k-list",
                default=None,
                help="sweep k: every tool at each listed k (e.g. 15,19,21,23,27,31)",
            )
            p.add_argument(
                "--t-list",
                default=None,
                help="sweep threads: every tool at each listed thread count (e.g. 1,2,4,8,16); "
                "the first should be 1",
            )
            p.add_argument(
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
            p.add_argument(
                "--skip-rounds",
                default="",
                help="--plan paper: comma-separated DATASET:LABEL rounds to leave out (resuming "
                "a run whose rounds already finished, e.g. sweep-chr1:kscale; their rows "
                "stay in the earlier run's output)",
            )
            p.add_argument(
                "--dry-run", action="store_true", help="print the rounds and exit"
            )
            p.add_argument(
                "--out",
                default=os.path.join(HERE, "results", "compare_tools_published_sets"),
            )
    args = ap.parse_args()
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
        pb.check_optik_fresh()
    mach = machine()
    if getattr(args, "cap_gb", 0):
        # A lower cap leaves the remaining memory to other work, and the memory gate plans within it.
        mach["mem_gib"] = min(mach["mem_gib"], args.cap_gb + 2)
    print(
        f"this machine: {mach['mem_gib']} GiB memory, {mach['disk_free_gb']} GB free disk"
    )
    measure.rusage_exec()  # built up front, or rebuilt when its source is newer
    if args.cmd == "smoke":
        return smoke(args, mach)
    skipped = {}
    for i, d in D.items():
        if "excluded" not in d and (
            not args.section or d["section"] in args.section.split(",")
        ):
            ok, why = preflight(i, d, mach)
            if not ok:
                skipped[i] = "; ".join(why)
    todo = selected(args, mach)
    if args.cmd == "fetch":
        # Sources throttle each connection (Zenodo about 250 kB/s) rather than the total, so datasets download side by side.
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as ex:
            futs = {ex.submit(fetch_one, i, d): i for i, d in todo}
            bad = []
            for f in concurrent.futures.as_completed(futs):
                try:
                    f.result()
                # Every failed dataset is reported, not just the first.
                except FETCH_ERRORS as e:
                    bad.append(futs[f])
                    print(f"  {futs[f]}: FAILED ({e})", flush=True)
        if bad:
            raise SystemExit("fetch failed for: " + ", ".join(bad))
        return
    modes = ["parity", "optimal"] if args.mode == "both" else [args.mode]
    k_list = [int(x) for x in args.k_list.split(",")] if args.k_list else None
    t_list = [int(x) for x in args.t_list.split(",")] if args.t_list else [args.threads]
    if args.t_list and t_list[0] != 1:
        # The compared tools' single-threaded queries run in the first round only.
        raise SystemExit(
            "--t-list must start at 1 (the round that also measures the single-threaded tools)"
        )
    # Unless --threads is given, the plan's paper rounds use every hardware thread.
    plan_t = args.threads
    if args.dry_run:
        n = 0
        skip_r = set(filter(None, args.skip_rounds.split(",")))
        for i, d in todo:
            rounds = (
                plan_rounds(i, d, plan_t)
                if args.plan == "paper"
                else [
                    (m, sweep_ks(i, d, k_list), t, False, False, "", {})
                    for m in modes
                    for t in t_list
                ]
            )
            for mode, ks_r, t, first_t, _, label, ex in rounds:
                if f"{i}:{label}" in skip_r:
                    print(
                        f"  {d['section']:5s} {i:28s} {mode:8s} {label:7s} skipped (--skip-rounds)"
                    )
                    continue
                n += 1
                what = (
                    f" streaming at k = {ex['stream_k'][0]} only"
                    if ex.get("stream_k") and d.get("stream")
                    else " builds only"
                    if ex.get("no_queries")
                    else " builds + optik queries"
                    if label == "tscale" and d["section"] == "query" and not first_t
                    else ""
                )
                print(
                    f"  {d['section']:5s} {i:28s} {mode:8s} {label:7s} k={ks_r or 'default'} t={t}{what}"
                )
        print(f"{n} rounds")
        return
    out = os.path.join(
        args.out, datetime.datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")
    )
    os.makedirs(out, exist_ok=True)
    pb.environment(out)
    cap = args.cap_gb or max(1, int(mach["mem_gib"] - 2))
    bench = pb.Bench(
        out, args.runs, args.warmup, cap, os.path.join(args.work, "_scratch")
    )
    for i, d in todo:
        if (
            not os.path.exists(os.path.join(DATA, i, "PROVENANCE.json"))
            and not d.get("local")
            and d["files"]
        ):
            print(
                f"  {i}: not fetched (run: compare_tools_published_sets.py fetch --only {i})"
            )
            continue
        # Check free disk again here.
        # The start-up check saw the disk before the earlier datasets wrote and removed their files, and a dataset that cannot fit becomes a SKIP row instead of a tool failing on a full disk partway through.
        need, how = disk_need(i, d)
        free = shutil.disk_usage(DATA).free / 1e9
        if need > free and not args.force:
            print(
                f"\n=== {d['section']} {i}: skipped -- needs ~{need:.0f} GB free disk ({how}), {free:.0f} GB free"
            )
            add_status(
                bench,
                "skip",
                section=d["section"],
                mode="-",
                dataset=i,
                tool="(all)",
                phase="disk",
                gate_detail=f"needs ~{need:.0f} GB free disk ({how}); {free:.0f} GB free at its turn",
            )
            report(bench, out, mach, skipped, args.report)
            continue
        ks = sweep_ks(i, d, k_list)
        if not args.keep_work:
            # Leftovers of an earlier, interrupted run are never reused.
            for entry in os.listdir(args.work) if os.path.isdir(args.work) else []:
                if entry == i or entry.startswith(i + "_"):
                    shutil.rmtree(os.path.join(args.work, entry), ignore_errors=True)

        def run_round(
            mode,
            ks,
            t,
            first_t,
            sweep_t,
            label="",
            only_tools=None,
            extra=None,
            *,
            i=i,
            d=d,
        ):
            """Run one planned round of dataset `i` and record a failed round as FAIL rows."""
            print(
                f"\n=== {d['section']} {i}: {d.get('title')} ({mode}{', ' + label if label else ''}, "
                f"k {ks or 'default'}, {t} threads)",
                flush=True,
            )
            dd = dict(d, **(extra or {}))
            try:
                if d["section"] == "query":
                    bench_query(
                        bench,
                        i,
                        dict(dd, k=ks),
                        mode,
                        args.work,
                        t,
                        mach,
                        first_t=first_t,
                        sweep_t=sweep_t,
                    )
                elif d["section"] == "dbg":
                    bench_dbg(
                        bench,
                        i,
                        dict(dd, k=ks),
                        mode,
                        args.work,
                        t,
                        mach,
                        only_tools=only_tools,
                    )
                elif d["section"] == "map":
                    for kk in ks or [None]:
                        bench_map(
                            bench, i, dict(dd, map_k=kk), mode, args.work, t, mach
                        )
                elif d["section"] == "align":
                    bench_align(bench, i, d, mode, args.work)
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
                    LOG.exception("bench error in round %s of %s", label or mode, i)
                    det = "bench error, " + det
                bench.add(
                    section=d["section"],
                    mode=mode,
                    dataset=i,
                    threads=t,
                    tool="(input)",
                    phase="prepare",
                    gate=False,
                    gate_detail=det[-400:],
                )
            report(bench, out, mach, skipped, args.report)

        if args.plan == "paper":
            skip_r = set(filter(None, args.skip_rounds.split(",")))
            for mode, ks_r, t, first_t, sweep_t, label, ex in plan_rounds(i, d, plan_t):
                if f"{i}:{label}" in skip_r:
                    print(
                        f"\n=== {d['section']} {i}: {label} round skipped (--skip-rounds)",
                        flush=True,
                    )
                    continue
                if d["section"] == "map" and label == "kscale":
                    ex = dict(ex, no_bwa=True)
                run_round(mode, ks_r, t, first_t, sweep_t, label, extra=ex)
                if d["section"] == "dbg" and label == "paper":
                    # Parity runs only where it changes a command and binds.
                    # optik's --max-mem is set to Cuttlefish's -m budget, run when a tool's optimal peak exceeded it.
                    budget_mb = max(1, int(mach["mem_gib"] * 0.5)) * 1024
                    peaks = [
                        r.get("rss_mb") or 0
                        for r in bench.rows
                        if r.get("dataset") == i
                        and r.get("mode") == "optimal"
                        and r.get("threads") == t
                        and r.get("tool") in ("optik", "cuttlefish")
                    ]
                    if any(p > budget_mb for p in peaks):
                        run_round(
                            "parity",
                            ks_r,
                            t,
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
                            dataset=i,
                            threads=t,
                            tool="optik",
                            gate_detail=f"every optimal peak is within the {budget_mb // 1024} GB budget, "
                            "so --max-mem would not bind (identical to the optimal row)",
                        )
        else:
            for mode in modes:
                if (
                    d["section"] in ("dbg", "map", "align")
                    and mode == "parity"
                    and len(modes) == 2
                ):
                    # --mode both runs one configuration each, the papers' own settings.  For dbg that
                    # is the paper's k and f0 with every tool on all threads at its documented best memory setting.
                    # --mode parity alone, the sweep, runs dbg with optik's --max-mem equal to Cuttlefish's -m.
                    continue
                for ti, t in enumerate(t_list if d["section"] != "align" else [1]):
                    run_round(mode, ks, t, ti == 0, args.t_list is not None)
        if not args.keep_work:
            # With its rows recorded, the dataset's working files (formatted inputs, indexes, SAMs, KMC databases) are removed, so only one dataset's working space is needed at a time.
            for e in os.listdir(args.work):
                if e == i or e.startswith(i + "_"):
                    shutil.rmtree(os.path.join(args.work, e), ignore_errors=True)
    report(bench, out, mach, skipped, args.report)
    print(f"\nwrote {out}/SUMMARY.md")
    # A failed correctness gate fails the run.
    # SKIP and CAP rows have no gate.
    nfail = sum(1 for r in bench.rows if r.get("gate") is False)
    if nfail:
        print(f"{nfail} gate failures")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
