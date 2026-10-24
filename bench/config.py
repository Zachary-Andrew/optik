"""Tool paths, dataset definitions and benchmark parameters.

This module is the single point of configuration of the benchmark suite.
PATHS points at the local tool installations and DATASETS holds the benchmark inputs.

Each dataset entry provides
    label     short display name for summary tables
    unitigs   stitched unitigs in FASTA format, .gz accepted
    queries   query sequences in FASTA or FASTQ format, .gz accepted
    k         k-mer length for indexing and querying
    desc      longer description shown by --list-datasets
"""

import os
from pathlib import Path

import hw

# Root directory of the optik project, one level above bench/
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# With TOOL_BENCH_TOOLS_DIR set, every external tool is resolved inside that directory, the output of build_tools.sh.
# Otherwise the older layout applies, where each tool lives in a sibling directory of optik/ under the bio/ parent (../sshash/, ../bifrost/ and so on).
#
# optik itself is always built locally with cargo.
#
# A tool that is not installed has None or a missing file as its entry, and the runner checks availability with check_binary() and skips it.

_TOOLS_DIR = os.environ.get("TOOL_BENCH_TOOLS_DIR")

# A temp/ directory next to bench/ is used when it exists.
if not _TOOLS_DIR:
    _auto = PROJECT_ROOT / "temp"
    if _auto.is_dir():
        _TOOLS_DIR = str(_auto)

if _TOOLS_DIR:
    _td = Path(_TOOLS_DIR)
    PATHS = {
        "optik": PROJECT_ROOT / "target/release/optik",
        "sshash": _td / "sshash/build/sshash",
        "fmsi": _td / "fmsi/fmsi",
        "kmercamel": _td / "kmercamel/kmercamel",
        "bifrost": _td / "bifrost/build/src/Bifrost",
        "cuttlefish": _td / "cuttlefish/target/release/cuttlefish",
        "minimap2": _td / "minimap2/minimap2",
        "sbwt": _td / "SBWT/build/bin/sbwt",
        "wfmash": _td / "wfmash/build/bin/wfmash",
        "wgatools": _td / "wgatools/target/release/wgatools",
    }
else:
    PATHS = {
        "optik": PROJECT_ROOT / "target/release/optik",
        "sshash": PROJECT_ROOT.parent / "sshash/build/sshash",
        "fmsi": PROJECT_ROOT.parent / "fmsi/fmsi",
        "kmercamel": PROJECT_ROOT.parent / "kmercamel/kmercamel",
        "bifrost": PROJECT_ROOT.parent / "bifrost/build/src/Bifrost",
        "cuttlefish": PROJECT_ROOT.parent / "cuttlefish/target/release/cuttlefish",
        "minimap2": PROJECT_ROOT.parent / "minimap2/minimap2",
        "sbwt": PROJECT_ROOT.parent / "SBWT/build/bin/sbwt",
        "wfmash": PROJECT_ROOT.parent / "wfmash/build/bin/wfmash",
        "wgatools": PROJECT_ROOT.parent / "wgatools/target/release/wgatools",
    }

# The unitig files come from the SSHash data repository and hold stitched unitigs of the reference genomes.
# Query files are sequencing reads (FASTQ) or whole-genome sequences (FASTA).
#
# With TOOL_BENCH_TOOLS_DIR set, data files are resolved inside that directory (in the cloned sshash/data/ subtree), and otherwise in the older sibling-directory layout.


def _sshash_data(relpath):
    """Resolve a path inside the sshash data directory.

    When build_tools.sh builds the tools into $TOOLS_DIR, the sshash repository is cloned there with its data/ subtree.
    Without TOOL_BENCH_TOOLS_DIR the older sibling path is used.
    """
    if _TOOLS_DIR:
        td = Path(_TOOLS_DIR)
        return td / "sshash" / "data" / relpath
    return PROJECT_ROOT.parent / "sshash" / "data" / relpath


DATASETS = {
    "ecoli_k31": {
        "label": "E. coli k31",
        "unitigs": _sshash_data("unitigs_stitched/ecoli1_k31_ust.fa.gz"),
        # ecoli1.fasta is the genome these unitigs were built from (99.87% of its k-mers are positive and the rest hold N).
        # SSHash's demo read set SRR5833294.10K comes from another organism, and 0 of its 460,000 k-mers occur in this index, so as the query it measured only misses and as the graph input it built an 11 kb graph of 300 unitigs.
        "queries": _sshash_data("queries/ecoli1.fasta.gz"),
        "dbg_input": _sshash_data("queries/ecoli1.fasta.gz"),
        "k": 31,
        "dbg_min_count": 1,
        "desc": "E. coli (ecoli1), 31-mer unitigs; queried with its own genome (all-positive stream)",
    },
    "salmonella_k31": {
        "label": "Salmonella k31",
        "unitigs": _sshash_data("unitigs_stitched/salmonella_enterica_k31_ust.fa.gz"),
        "queries": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "dbg_input": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "k": 31,
        "dbg_min_count": 1,
        "desc": "S. enterica, 31-mer unitigs",
    },
    "se_k31": {
        "label": "SE k31",
        "unitigs": _sshash_data("unitigs_stitched/se.ust.k31.fa.gz"),
        "queries": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "dbg_input": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "k": 31,
        "dbg_min_count": 1,
        # SE is S. enterica, whose genome hits 97.5% of these unitigs' 31-mers (E. coli 1.8%).
        # The query, graph and mapping inputs used to be ecoli1, so every non-index row of an "SE" run repeated the E. coli run.
        "desc": "SE dataset (S. enterica), 31-mer unitigs; queried with the S. enterica genome",
    },
    "se_k47": {
        "label": "SE k47",
        "unitigs": _sshash_data("unitigs_stitched/se.ust.k47.fa.gz"),
        "queries": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "dbg_input": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "k": 47,
        "dbg_min_count": 1,
        "desc": "SE dataset (S. enterica), 47-mer unitigs",
    },
    "se_k63": {
        "label": "SE k63",
        "unitigs": _sshash_data("unitigs_stitched/se.ust.k63.fa.gz"),
        "queries": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "dbg_input": _sshash_data("queries/salmonella_enterica.fasta.gz"),
        "k": 63,
        "dbg_min_count": 1,
        "desc": "SE dataset (S. enterica), 63-mer unitigs",
    },
}

# OPTIK_BIN points every harness at another optik binary, for A/B runs or for a candidate build while a benchmark keeps using target/release.
if os.environ.get("OPTIK_BIN"):
    PATHS["optik"] = Path(os.environ["OPTIK_BIN"]).resolve()

# --runs, --warmup and the other options override these on the command line.

BENCHMARK = {
    "runs": 5,
    "warmup": 2,
    "timeout": 3600,
    # builds use every logical CPU, and queries run single-threaded, as every query paper measures them
    "threads": {"build": hw.cpus(), "query": 1, "query_sweep": "auto"},
    "max_ram_gb": 16,
    "tmp_dir": Path("/tmp/tool_bench"),
}
