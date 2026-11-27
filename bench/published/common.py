"""Holds the constants and helpers that every phase of the published-set comparison uses."""

import logging
import os
import re
import subprocess
import sys

import compare_tools_small_sets
import hw
from diagnostics import DETAIL_CHARS, LOG_TAIL_CHARS
from published_catalog import DATA_ROOT, TEMP_DIR
from small_common import BIN
from units import BYTES_PER_GIB, BYTES_PER_MIB

# Errors of a round or tool that do not stop the run, reported with their traceback
LOG = logging.getLogger(__name__)


# The errors a download or a source lookup raises, rangefetch's SystemExit among them
FETCH_ERRORS = (OSError, ValueError, KeyError, subprocess.SubprocessError, SystemExit)


# This package is one level below bench/, which holds the sibling modules and scripts.
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


ROOT = os.path.dirname(HERE)


sys.path.insert(0, HERE)


# k of the published k-mer counts (kmers31), of the scaling run behind the memory model, and of every query set's streaming rows
REFERENCE_K = 31


# k at which the subsampled (sub10) variant of a query set is published
SUBSAMPLED_K = 23


# GB added to every peak estimate for a tool's fixed overhead
EST_MEM_OVERHEAD_GB = 0.2


# For the disk model, compressed FASTA expands DECOMPRESSED_PER_GB times and holds about KMERS_PER_COMPRESSED_GB k-mers per compressed GB.
# A read set's working files take READS_WORK_FACTOR times its input, and a genome GENOME_WORK_FACTOR times its decompressed size.
# A mapping set keeps MAP_DOWNLOAD_FACTOR times its download and MAP_LOCAL_FACTOR times its local reference, plus about SAM_BYTES_PER_READ of SAM per read.
DECOMPRESSED_PER_GB = 4.5


KMERS_PER_COMPRESSED_GB = 4.5e9


READS_WORK_FACTOR = 3.0


GENOME_WORK_FACTOR = 2


MAP_DOWNLOAD_FACTOR = 3


MAP_LOCAL_FACTOR = 2


SAM_BYTES_PER_READ = 2000


MATES_PER_PAIR = 2


# size assumed for an NCBI assembly zip, whose size the source does not publish
NCBI_ASSEMBLY_BYTES = 1e9


# whole-file downloads: attempts, and a pause growing by DOWNLOAD_RETRY_STEP_S per attempt
DOWNLOAD_ATTEMPTS = 8


DOWNLOAD_RETRY_STEP_S = 5


PREFIX_EXT = ".prefix"


EXEC_MODE = 0o755


# KMC gets KMC_MEM_SHARE of the usable memory, capped at KMC_MAX_MEM_GB.
# Memory-budgeted tools, meanwhile, share TOOL_BUDGET_SHARE of it.
KMC_MAX_MEM_GB = 16


KMC_MEM_SHARE = 0.6


TOOL_BUDGET_SHARE = 0.5


# query workloads: k-mers per workload, a smoke run's cap and its streaming bases, and the seed of the one-k-mer load probe
MAX_QUERIES = 10**6


SMOKE_MAX_QUERIES = 10**4


SMOKE_STREAM_BASES = 2 * 10**6


LOAD_PROBE_SEED = 3


# minimap2 accepts k up to this
MM2_MAX_K = 28


# argv positions where a k/w override goes: after "optik map" and after "minimap2 -ax PRESET"
OPTIK_ARGV_PREFIX = 2


MM2_ARGV_PREFIX = 3


# mapeval rows are compared at MAPQ >= MAPEVAL_MAPQ_CUT, as the papers plot them
MAPEVAL_MAPQ_CUT = 10


# summary lines of pafcmp's output kept in a row
PAFCMP_SUMMARY_LINES = 4


# the BiWFA paper's ONT set holds this many pairs
BIWFA_ONT_PAIRS = 48


# smoke runs fetch this prefix of each file
SMOKE_PREFIX_BYTES = 24 * BYTES_PER_MIB


# fetch progress newer than this many seconds means a fetch is running
PROGRESS_FRESH_S = 60


# planned rounds: (mode, k list, threads, first round, sweep, label, extras)
ROUND_FIELDS = 7


# CLI defaults
DEFAULT_FETCH_JOBS = 6


DEFAULT_RUNS = 3


# exit status reported when a command produced no run at all, and the primaries count without a SAM
NO_RUN_EXIT = -1


NO_SAM = -1


# Peak memory per distinct k-mer for each tool, by k, as bench/scaling_bench.py measured it on human chr1.
# The scaling logs in SCALING_LOGS are read when they exist.
# Otherwise the fallback below applies, taken from the chr1 run of 2026-10-12 at optik e6202c3, with each tool's peak RSS in MB and the distinct k-mer count _kmers behind it.
SCALING_LOGS = [
    os.path.join(TEMP_DIR, "scaling_rerun.log"),
    os.path.join(TEMP_DIR, "scaling_rerun_k1521.log"),
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


# SSHash (-g), optik dbruijn (--max-mem) and Cuttlefish (-m) take this machine's memory budget and keep to it.
# A graph over the budget makes optik dbruijn fail with a named error.
BUDGETED = ("sshash", "optik", "cuttlefish")


_SIZES = os.path.join(DATA_ROOT.path, "_download_sizes.json")


# Index bits per k-mer of each tool (chr1, k = 31, the same scaling run), which shows how much disk one dataset's indexes take next to their input.
INDEX_BITS = {
    "optik-spss": 8.37,
    "sshash": 7.91,
    "sbwt": 10.54,
    "fmsi": 3.28,
    "optik-fm": 5.43,
}


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
            f"resident memory reached the {bench.cap // BYTES_PER_GIB} GiB cap (exit {rc_}, peak {rss:.0f} MB)",
        )
    # measure.Result keeps LOG_TAIL_CHARS of stderr, and all of it is searched for an allocation failure.
    match = MEM_SIG.search(err[-LOG_TAIL_CHARS:])
    if match:
        return (
            "cap",
            f"memory allocation failed ({match.group(0)!r}, exit {rc_}, peak {rss:.0f} MB)",
        )
    return "fail", f"exit {rc_}: {err[-DETAIL_CHARS:]}"


def add_status(bench, status, **row):
    """Add a SKIP, CAP or NOT RUN row, which is neither a pass nor a wrong answer."""
    row["status"] = status
    row["gate"] = None
    bench.rows.append(row)
    print(
        f"  [{status.upper()}] "
        + "  ".join(
            f"{key}={compare_tools_small_sets._fmt(value)}"
            for key, value in row.items()
            if key not in ("section", "mode", "gate", "status")
        ),
        flush=True,
    )


SPSS_TRIES = 3  # GGCAT eulertig attempts before its unitigs stand in (GGCAT issue #65)


# The KMC oracles of external_oracles.py, shared with compare_tools_small_sets
KMC = BIN["kmc"]


SMOKE_BASES = 2_000_000


# Headings of the paper sections in the sectioned results file
SECTION_TITLES = {
    "query": "k-mer membership queries against SSHash, SBWT and FMSI on their papers' datasets",
    "dbg": "compacted de Bruijn graphs against Cuttlefish and GGCAT on their papers' datasets",
    "map": "read mapping against minimap2, BWA-MEM and Winnowmap on their papers' simulations",
    "align": "pairwise alignment against WFA2-lib BiWFA on its paper's sequence pairs",
}


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


# Matched thread count of the k-scaling plan.
# Cuttlefish 2 Table 4 uses 8 threads, and a machine with fewer usable CPUs runs at all of them, which the round's label records.
PAPER_PARITY_T = 8


PARITY_T = min(PAPER_PARITY_T, hw.cpus())
