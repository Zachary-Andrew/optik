"""Holds the constants, the Bench harness and the helpers that every small-set section uses."""

import contextlib
import datetime
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

import measure
import memcap
from diagnostics import DETAIL_CHARS
from small_common import BIN, TEMP_DIR
from units import BYTES_PER_GIB, BYTES_PER_MIB, KIB_PER_MIB

# This package is one level below bench/, which holds the sibling modules and scripts.
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


ROOT = os.path.dirname(HERE)


sys.path.insert(0, HERE)


# seconds a machine-description shell command may take
SHELL_TIMEOUT_S = 30


# k of SSHash's published unitig sets, which every k-mer index is built and queried at
SET_K = 31


# KMC's memory for the truth databases, in GB
KMC_MEM_GB = 4


# SSHash minimizer length of the parity runs (k/2 for the sets without a published m)
PARITY_M = 15


# The k-mer workload is KMER_READS reads of KMER_READ_LEN bases, each inside a single unitig and half of them reverse-complemented.
# Reads get KMER_ERR_RATE errors (SUB_SHARE of them substitutions) plus KMER_INDEL_RATE insertions and deletions.
# RANDOM_KMERS random k-mers are added to the queries as well.
KMER_SEED = 7


KMER_READS = 20000


KMER_READ_LEN = 150


REVCOMP_PROB = 0.5


KMER_ERR_RATE = 0.01


KMER_INDEL_RATE = 0.002


RANDOM_KMERS = 1_000_000


# id that SBWT and optik --ids write for an absent k-mer, and ids_check's result when the run failed
ABSENT_KMER_ID = -1


IDS_CHECK_FAILED = -1


# A short read is placed correctly within max(MIN_ORIGIN_TOL, length / ORIGIN_TOL_DIV) bp.
# A long read, on the other hand, needs MIN_LONG_OVERLAP of its alignment to overlap its true span.
MIN_ORIGIN_TOL = 20


ORIGIN_TOL_DIV = 10


MIN_LONG_OVERLAP = 0.5


MAX_EXAMPLES = 3


# synthetic 100 Mbp genome: SYN_CHROMS chromosomes of SYN_CHROM_LEN bases
SYN_SEED = 100


SYN_CHROMS = 4


SYN_CHROM_LEN = 25_000_000


# chimeric ONT reads get CHIMERA_TAIL_RANGE (randrange bounds) bases of foreign sequence appended
CHIMERA_SEED = 11


CHIMERA_TAIL_RANGE = (5000, 20000)


# negative control: random reads that must not map, at most MAX_RANDOM_MAPPED_PCT percent of them
NEG_SEED = 3


NEG_READS = 50000


NEG_READ_LEN = 150


MAX_RANDOM_MAPPED_PCT = 0.1


# Align grid cells are (length, error rate, pairs), and --quick keeps only the first QUICK_ALIGN_CELLS.
ALIGN_GRID = [
    (150, 0.01, 20000),
    (150, 0.05, 20000),
    (1000, 0.05, 2000),
    (1000, 0.15, 500),
    (10000, 0.05, 40),
    (10000, 0.15, 10),
    (100000, 0.05, 2),
]


QUICK_ALIGN_CELLS = 3


# a cell's seed is length x ALIGN_SEED_LEN_MULT + error rate x ALIGN_SEED_ERR_SCALE
ALIGN_SEED_LEN_MULT = 13


ALIGN_SEED_ERR_SCALE = 1000


# an error is a substitution, insertion or deletion in these shares
SUB_SHARE = 0.6


INS_SHARE = 0.2


DEL_SHARE = 0.2


DEFAULT_RUNS = 3


# exit status reported when a command produced no run at all
NO_RUN_EXIT = -1


class Bench:
    """Runs and times the bench's commands and collects its result rows."""

    def __init__(self, out, runs, warmup, cap_gb, scratch):
        self.out, self.runs, self.warmup = out, runs, warmup
        # Hyperfine JSON, stderr and rusage files go to one directory per process, so two benches sharing a work directory never read each other's results.
        self.scratch = os.path.join(scratch, str(os.getpid()))
        os.makedirs(self.scratch, exist_ok=True)
        self.cap = cap_gb * BYTES_PER_GIB
        self.last_capped = False
        self.last_io = (None, None)
        self.cmdlog = measure.AppendLog(os.path.join(out, "commands.log"))
        self.rows = []

    def log(self, tag, cmd):
        """Append a tagged command line to the run's command log."""
        self.cmdlog.write(f"[{tag}] {' '.join(map(str, cmd))}\n")

    def once(self, tag, cmd, cwd=None, stdout=None, cap=None):
        """Run a command once, untimed, and return (returncode, wall_s, peak_rss_mb, stderr).

        Used for an index the timed runs read, or for an extra build whose size is reported.
        The command is launched via rusage_exec, a 1 MB process, because a child forked from this large Python process would report the parent's RSS as its own peak.
        `cap` gives the resident bytes the whole run may use, enforced by memcap's watchdog over its process group on any Linux.
        It defaults to self.cap, and 0 means no cap.
        Every number a results row reports comes from measure(), not from here."""
        self.log(tag, cmd)
        errf = os.path.join(self.scratch, "stderr.tmp")
        resf = os.path.join(self.scratch, "rusage.tmp")
        if os.path.exists(resf):
            os.remove(resf)
        with contextlib.ExitStack() as stack:
            err_handle = stack.enter_context(open(errf, "w"))
            out = (
                stack.enter_context(open(stdout, "w")) if stdout else subprocess.DEVNULL
            )
            result = memcap.run(
                [measure.rusage_exec(), "0", resf] + [str(arg) for arg in cmd],
                self.cap if cap is None else cap,
                cwd=cwd,
                stdout=out,
                stderr=err_handle,
            )
        self.last_capped = result.capped
        if not os.path.exists(resf):  # the group was killed before rusage_exec wrote it
            return (
                result.returncode,
                0.0,
                (self.cap if cap is None else cap) / BYTES_PER_MIB,
                Path(errf).read_text(),
            )
        run = measure.RusageLine.parse(Path(resf).read_text())
        return run.exit_code, run.wall_s, run.peak_rss_kib / KIB_PER_MIB, Path(errf).read_text()

    def repeat(self, tag, cmd, stdout, cwd=None):
        """Run a program that times itself, such as a microbenchmark harness printing its own seconds.

        The first --warmup runs are thrown away and the next --runs kept.
        Returns [(returncode, peak_rss_mb, stdout text)] of the kept runs, stopping at the first failure, whose entry comes last."""
        kept = []
        for i in range(self.warmup + self.runs):
            rc_, _, rss, err = self.once(
                f"{tag} ({'warmup' if i < self.warmup else 'run'} {i + 1})",
                cmd,
                cwd=cwd,
                stdout=stdout,
            )
            if rc_ != 0:
                return kept + [(rc_, rss, err)]
            if i >= self.warmup:
                kept.append((rc_, rss, Path(stdout).read_text()))
        return kept

    def measure(self, tag, cmd, cwd=None, stdout=None, inputs=(), cap=None):
        """Time a command through bench/measure.py.

        Hyperfine runs with this bench's --warmup and --runs, drops the warmups and computes statistics over the measured runs.
        Each run reports peak RSS and block I/O, and the inputs are cached before every run.
        `stdout` names where the last measured run's standard output is copied for the correctness gate.
        Returns a measure.Result, whose .wall_row() is the {mean, sd, ...} a row prints."""
        result = measure.measure(
            cmd,
            self.warmup,
            self.runs,
            self.cap if cap is None else cap,
            inputs,
            scratch=self.scratch,
            cwd=cwd,
            tag=tag,
            log=self.cmdlog,
        )
        if stdout and os.path.exists(result.out):
            shutil.copyfile(result.out, stdout)
        self.last_capped = result.capped
        self.last_io = (result.written, result.read)
        return result

    def add(self, **row):
        """Append a result row and print it, failing the row if any measurement failed."""
        # A failed measurement makes a failed row, never a silent PASS.
        for key, value in row.items():
            if isinstance(value, dict) and "error" in value:
                row["gate"] = False
                row["gate_detail"] = (
                    row.get("gate_detail") or ""
                ) + f" timing of {key} failed: {value['error'][-DETAIL_CHARS:]}"
        self.rows.append(row)
        flag = "PASS" if row.get("gate") else "FAIL"
        keys = [
            key for key in row if key not in ("section", "mode", "gate", "gate_detail", "cmd")
        ]
        print(
            f"  [{flag}] " + "  ".join(f"{key}={_fmt(row[key])}" for key in keys), flush=True
        )
        if not row.get("gate"):
            print(f"         gate: {row.get('gate_detail')}", flush=True)


def _shell_quote(arg):
    """c quoted for a POSIX shell when it holds anything beyond plain path characters."""
    arg = str(arg)
    return (
        arg if re.fullmatch(r"[\w./=:,+-]+", arg) else "'" + arg.replace("'", "'\\''") + "'"
    )


def _fmt(value):
    """Format a row value for the tables, printing timings as mean±sd seconds."""
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, dict) and "mean" in value:
        # Hyperfine reports no standard deviation for a single run.
        return f"{value['mean']:.4g}s" + (
            f"±{value['sd']:.2g}" if value.get("sd") is not None else ""
        )
    return str(value)


def environment(out):
    """Record the machine, compiler and tool revisions of this run in `out`."""

    def shell_output(cmd):
        """Output of shell command cmd, or "n/a (reason)" if it fails."""
        try:
            return subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=SHELL_TIMEOUT_S, check=False
            ).stdout.strip()
        except (subprocess.SubprocessError, OSError) as err:
            return f"n/a ({err})"

    # Everything that can move a timing is recorded, from the CPU model and governor to transparent huge pages and the build flags.
    env = {
        "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "cpu": shell_output("lscpu | grep 'Model name' | sed 's/.*: *//'"),
        "cores_threads": shell_output("nproc") + " logical",
        "physical_cores": shell_output("lscpu -p=core | grep -v '^#' | sort -u | wc -l"),
        "memory": shell_output("free -g | awk '/Mem/{print $2\" GiB\"}'"),
        "kernel": platform.release(),
        "thp": shell_output("cat /sys/kernel/mm/transparent_hugepage/enabled"),
        "governor": shell_output("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "gcc": shell_output("g++ --version | head -1"),
        "rustc": shell_output("rustc --version"),
        "optik_git": shell_output(f"git -C {ROOT} rev-parse --short HEAD")
        + (" (dirty)" if shell_output(f"git -C {ROOT} status --porcelain -uno") else ""),
        "build_flags": "optik: -O3 -march=native (build.rs) + target-cpu=native; competitors: "
        "bench/build_tools.sh NATIVE=1 (-O3 -march=native; SBWT -O3 -mbmi2 -mpopcnt)",
        "tools": {},
    }
    # Each competitor's source revision pins exactly which build was measured.
    for name in (
        "sshash",
        "SBWT",
        "fmsi",
        "kmercamel",
        "cuttlefish",
        "bifrost",
        "minimap2",
        "wfmash",
        "WFA2-lib",
        "KMC",
    ):
        env["tools"][name] = shell_output(
            f"git -C {TEMP_DIR}/{name} rev-parse --short HEAD 2>/dev/null || echo n/a"
        )
    env["tools"]["minimap2_version"] = shell_output(f"{BIN['minimap2']} --version")
    env["tools"]["cuttlefish_version"] = shell_output(f"{BIN['cuttlefish']} --version")
    Path(os.path.join(out, "environment.json")).write_text(json.dumps(env, indent=1))
    return env


def check_optik_fresh():
    """Refuse to benchmark an optik binary older than any of its sources, since a stale build would be measured and reported as the current code."""
    exe = BIN["optik"]
    srcs = [
        os.path.join(root, fname)
        for top in ("src", "include", "tools/src")
        for root, _, files in os.walk(os.path.join(ROOT, top))
        for fname in files
        if fname.endswith((".rs", ".cpp", ".hpp", ".h"))
    ] + [os.path.join(ROOT, name) for name in ("build.rs", "Cargo.toml")]
    newest = max(srcs, key=os.path.getmtime)
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(newest):
        sys.exit(
            f"error: {exe} is older than {os.path.relpath(newest, ROOT)}: run `cargo build --release` first"
        )


def build_wfab(out_dir):
    """Path of the align bench binary (optik's WFA/BiWFA and WFA2-lib).

    It is rebuilt whenever a source is newer, since a binary kept across code changes would time an old aligner."""
    exe = os.path.join(out_dir, "wfab")
    # the compiler's inputs, plus the header they include, which only counts for staleness
    compiled = [
        os.path.join(HERE, "cpp/wfa_affine_bench.cpp"),
        os.path.join(TEMP_DIR, "WFA2-lib/bindings/cpp/WFAligner.cpp"),
        os.path.join(TEMP_DIR, "WFA2-lib/lib/libwfa.a"),
    ]
    srcs = [os.path.join(ROOT, "include/wfa_affine.hpp"), *compiled]
    if not os.path.exists(exe) or os.path.getmtime(exe) < max(
        os.path.getmtime(path) for path in srcs
    ):
        subprocess.run(
            [
                "g++",
                "-std=c++17",
                "-O3",
                "-march=native",
                "-I",
                os.path.join(ROOT, "include"),
                "-I",
                os.path.join(TEMP_DIR, "WFA2-lib"),
                *compiled,
                "-fopenmp",
                "-o",
                exe + ".part",
            ],
            check=True,
            # Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0.
            env=dict(os.environ, NIX_ENFORCE_NO_NATIVE="0"),
        )
        os.replace(exe + ".part", exe)
    return exe


KMER_SETS = {
    "ecoli": ("ecoli1_k31_ust.fa.gz", 13),  # SSHash README gives -m 13 for one genome
    "salmonella_100": (
        "salmonella_100_k31_ust.fa.gz",
        15,
    ),  # -m 15 for the 100-strain set (README)
    "penicillium": (
        "penicillium_chrysogenum_k31_ust.fa.gz",
        15,
    ),  # no published m, so k/2
}
