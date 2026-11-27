"""Base class and helper functions shared by every tool benchmark runner.

The ToolRunner base class defines the interface each runner implements.
The module-level functions run commands, wrap hyperfine, decompress inputs and generate random k-mers.

Nix library paths
    Some tools, kmercamel among them, need shared libraries from the Nix store.
    Those directories are detected and prepended to LD_LIBRARY_PATH so that subprocesses find the right .so files.

RSS measurement
    Hyperfine records each run's VmHWM from /proc when available, and the bytes are reported in megabytes.
"""

import json
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import hw
import measure
from diagnostics import DETAIL_CHARS
from fastq import (
    FASTQ_LINES_PER_RECORD,
    FASTQ_SEQ_LINE,
)
from units import BYTES_PER_MIB

# Nix store paths of shared libraries that some tools need in the development environment.
#
# Only list a library here when a tool cannot find it on its own, such as glpk, which kmercamel and fmsi link in a non-standard way.
# glibc, libgcc and libstdc++ store paths must not be hardcoded.
# Every dynamically linked binary already refers to the glibc it was built against, and Nix store paths are content-hashed per generation.
# A hardcoded path goes stale after the next system rebuild while the old generation stays on disk until garbage collection, so it would be injected into every subprocess's LD_LIBRARY_PATH ahead of the correct one.
# That mismatch produces errors such as
#   undefined symbol: __nptl_change_stack_perm, version GLIBC_PRIVATE
# for every tool tool_bench launches, however correctly each one was built.
NIX_GLPK_LIB = "/nix/store/halm4fcvdny1842yaww22w160zjzc88k-glpk-5.0/lib"

NIX_LIB_PATHS = [NIX_GLPK_LIB]


# Build threads when the runtime context names none
DEFAULT_BUILD_THREADS = 8
# Hyperfine measurement runs when the runtime context names none
DEFAULT_HYPERFINE_RUNS = 3
# last_bad value before any non-ACGT base was seen
NO_BAD_BASE = -1


# RSS scaling compares the first and last rows of a thread sweep, so it needs two of them.
MIN_SCALING_ROWS = 2
# Packed k-mers hold two bits per base, so BASE_MASK - code complements a base.
BITS_PER_BASE = 2
BASE_MASK = (1 << BITS_PER_BASE) - 1
# Lemire's method draws 64-bit words and keeps the high half of the product.
WORD_BITS = 64


def _build_ld_path():
    """Prepend Nix store library directories to the existing LD_LIBRARY_PATH.

    Some tools in this suite (kmercamel, fmsi) link against GLPK and other libraries installed through Nix rather than in system paths, and their subprocesses fail with missing .so errors without these directories.
    """
    current = os.environ.get("LD_LIBRARY_PATH", "")
    for lib_dir in NIX_LIB_PATHS:
        if os.path.isdir(lib_dir) and lib_dir not in current:
            current = f"{lib_dir}:{current}"
    return current


# The library path is built once at import, so every subprocess inherits it.
DEFAULT_ENV_VARS = {"LD_LIBRARY_PATH": _build_ld_path()} if _build_ld_path() else {}


def check_binary(path):
    """True if `path` is an executable file on disk.

    Runners check this before benchmarking a tool.
    Skipping missing tools quietly is left to the CLI layer above.
    """
    file_path = Path(path)
    return file_path.is_file() and os.access(file_path, os.X_OK)


def run_cmd(cmd, timeout=None, verbose=False, env=None, cwd=None):
    """Run a command once and return the CompletedProcess and wall time.

    This is used for non-timed operations like parsing tool output.
    Benchmark timing is handled by run_hyperfine instead.

    Parameters
    cmd : list[str]
        The command and arguments to execute.
    timeout : int or None
        Kill the process after this many seconds.
    verbose : bool
        Print the command and its output to stderr.
    env : dict or None
        Extra environment variables merged over the defaults.
    cwd : str or None
        Working directory for the subprocess.

    Returns
    tuple[subprocess.CompletedProcess, float]
        The process result object and elapsed wall time in seconds.
    """
    if verbose:
        print(f"  [cmd] {' '.join(str(arg) for arg in cmd)}", file=sys.stderr)

    # Environment = the OS environment, then the Nix library paths, then any extras.
    merged_env = {**os.environ, **DEFAULT_ENV_VARS, **(env or {})}
    start = time.time()
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=merged_env,
        cwd=cwd,
        check=False,
    )
    elapsed = time.time() - start

    if verbose and proc.stdout:
        for line in proc.stdout.strip().split("\n"):
            print(f"    [out] {line}", file=sys.stderr)
    if verbose and proc.stderr:
        for line in proc.stderr.strip().split("\n"):
            print(f"    [err] {line}", file=sys.stderr)

    return proc, elapsed


def _locate_hyperfine():
    """Search the usual locations for the hyperfine binary.

    Hyperfine gives calibrated wall-clock timing with warmup runs, a standard deviation and peak RSS.
    When it is missing, timing falls back to a single run through run_cmd.
    """
    candidates = [
        "/run/current-system/sw/bin/hyperfine",
        "/usr/bin/hyperfine",
        "hyperfine",
    ]
    for candidate in candidates:
        if shutil.which(candidate) or (
            Path(candidate).is_file() and os.access(candidate, os.X_OK)
        ):
            return candidate
    return None


# hyperfine is located once, at import.
HYPERFINE_BIN = _locate_hyperfine()


def physical_core_list():
    """One logical CPU id per physical core, in ascending physical-core order.

    Reads /proc/cpuinfo (physical id + core id) so benchmarks pinned with taskset run one thread per physical core instead of spreading across SMT siblings.
    Returns None when the parsing fails (callers then skip pinning).
    """
    first_by_core = {}
    try:
        for block in Path("/proc/cpuinfo").read_text().split("\n\n"):
            proc = phys = core = None
            for line in block.splitlines():
                if line.startswith("processor"):
                    proc = line.split(":", 1)[1].strip()
                elif line.startswith("physical id"):
                    phys = line.split(":", 1)[1].strip()
                elif line.startswith("core id"):
                    core = line.split(":", 1)[1].strip()
            if proc is not None and phys is not None and core is not None:
                first_by_core.setdefault((int(phys), int(core)), int(proc))
    except OSError:
        return None
    if not first_by_core:
        return None
    return [first_by_core[core] for core in sorted(first_by_core)]


def taskset_prefix(cpu_list):
    """taskset prefix to pin a command to `cpu_list`, or [] if unavailable."""
    if not cpu_list:
        return []
    taskset = shutil.which("taskset")
    if not taskset:
        return []
    cpus = ",".join(str(cpu) for cpu in cpu_list)
    return [taskset, "-c", cpus]


def run_hyperfine(
    cmd,
    runs=DEFAULT_HYPERFINE_RUNS,
    warmup=1,
    timeout=None,
    verbose=False,
    env=None,
    cwd=None,
    ignore_failure=False,
    pin_cpus=None,
):
    """Run a command under hyperfine and return timing statistics.

    Hyperfine runs the command several times after a warmup phase and reports the mean, standard deviation, minimum, maximum and per-run RSS.
    The command is wrapped so that the stdout and stderr of the last measured run are kept for parsing.

    Parameters
    cmd : list[str] or str
        The command to benchmark.
        Lists are shell-quoted and joined.
    runs : int
        Number of measurement runs after warmup.
    warmup : int
        Number of warmup runs discarded from statistics.
    timeout : int or None
        Total time limit for hyperfine.
    verbose : bool
        Print the wrapped command to stderr.
    env : dict or None
        Extra environment variables.
    cwd : str or None
        Working directory for the subprocess.
    ignore_failure : bool
        Pass --ignore-failure to hyperfine so it continues even when the benchmarked command returns a non-zero exit code.
        The exit codes are recorded in the JSON results.
    pin_cpus : list[int] or None
        When given, the command is run under `taskset -c <cpus>` so all hyperfine runs and their threads stay on the listed CPUs.

    Returns
    dict
        Keys include mean_s, stddev_s, min_s, max_s, the times_s list, mean_rss_mb, max_rss_mb, stdout, stderr, and the raw hyperfine result dict.
        On error the dict has a single key "error".
    """
    # Without hyperfine, a single run through run_cmd stands in.
    if HYPERFINE_BIN is None:
        if verbose:
            print(
                "  [warn] hyperfine not found; falling back to single run",
                file=sys.stderr,
            )
        proc, elapsed = run_cmd(cmd, timeout=timeout, verbose=verbose, env=env)
        return {
            "mean_s": elapsed,
            "stddev_s": 0.0,
            "min_s": elapsed,
            "max_s": elapsed,
            "times_s": [elapsed],
            "mean_rss_mb": None,
            "max_rss_mb": None,
            "exit_codes": [proc.returncode],
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "proc": proc,
        }

    # Hyperfine's output files go to a temporary directory.
    work_dir = Path(tempfile.mkdtemp(prefix="hf_"))
    stdout_path = work_dir / "stdout.txt"
    stderr_path = work_dir / "stderr.txt"
    json_path = work_dir / "results.json"

    # Hyperfine discards command output by default, so the command string itself redirects stdout and stderr to files that are read once hyperfine finishes.
    if isinstance(cmd, (list, tuple)):
        cmd_str = " ".join(shlex.quote(str(arg)) for arg in cmd)
    else:
        cmd_str = str(cmd)
    prefix = taskset_prefix(pin_cpus)
    if prefix:
        cmd_str = " ".join(shlex.quote(str(arg)) for arg in prefix) + " " + cmd_str
    # With the redirect inside the timed shell, hyperfine never buffers the tool's output, and the last run's files remain for the correctness check.
    wrapped = f"{{ {cmd_str}; }} >{stdout_path} 2>{stderr_path}"

    hf_command = [
        HYPERFINE_BIN,
        "--runs",
        str(runs),
        "--warmup",
        str(warmup),
        "--export-json",
        str(json_path),
        "--style",
        "none",
        "--shell",
        "bash",
    ]
    if ignore_failure:
        hf_command.append("--ignore-failure")
    hf_command.append(wrapped)

    if verbose:
        print(f"  [hyperfine] {wrapped}", file=sys.stderr)

    merged_env = {**os.environ, **DEFAULT_ENV_VARS, **(env or {})}
    try:
        proc = subprocess.run(
            hf_command,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=merged_env,
            cwd=cwd,
            check=False,
        )
    except subprocess.TimeoutExpired:
        shutil.rmtree(work_dir, ignore_errors=True)
        return {"error": "hyperfine timed out"}

    if proc.returncode != 0:
        err = (proc.stderr or "")[:DETAIL_CHARS]
        shutil.rmtree(work_dir, ignore_errors=True)
        return {"error": f"hyperfine failed: {err}"}

    # Read the JSON results hyperfine wrote.
    try:
        hf_data = json.loads(json_path.read_text())
    except (json.JSONDecodeError, FileNotFoundError) as exc:
        shutil.rmtree(work_dir, ignore_errors=True)
        return {"error": f"hyperfine JSON parse: {exc}"}

    first_result = measure.hyperfine_result(hf_data)
    stdout_text = stdout_path.read_text() if stdout_path.exists() else ""
    stderr_text = stderr_path.read_text() if stderr_path.exists() else ""

    shutil.rmtree(work_dir, ignore_errors=True)

    # Per-run VmHWM is in bytes and is reported in megabytes.
    memory_usage_bytes = first_result.get("memory_usage_byte") or []
    if memory_usage_bytes:
        mean_rss = round(
            sum(memory_usage_bytes) / len(memory_usage_bytes) / BYTES_PER_MIB, 1
        )
        max_rss = round(max(memory_usage_bytes) / BYTES_PER_MIB, 1)
    else:
        mean_rss = None
        max_rss = None

    return {
        "mean_s": first_result["mean"],
        "stddev_s": first_result["stddev"]
        if first_result["stddev"] is not None
        else 0.0,
        "median_s": first_result.get("median", first_result["mean"]),
        "min_s": first_result["min"],
        "max_s": first_result["max"],
        "times_s": first_result.get("times", [first_result["mean"]]),
        "mean_rss_mb": mean_rss,
        "max_rss_mb": max_rss,
        "exit_codes": first_result.get("exit_codes", [0]),
        "stdout": stdout_text,
        "stderr": stderr_text,
        "raw": first_result,
    }


def sweep_thread_counts(ctx):
    """Resolve the configured thread sweep to a concrete list of thread counts.

    The sweep spec is read from ``ctx["threads_sweep"]`` and takes one of three forms.
        "auto"   hw.thread_points(), powers of two up to the usable CPUs, plus the exact CPU count when it is not a power of two, so 12 CPUs give 1, 2, 4, 8, 12.
        "none"   or any false value disables the sweep and gives [].
        "1,2,4"  explicit comma-separated thread counts.

    Returns a de-duplicated list of positive integers in ascending order.
    """
    spec = str(ctx.get("threads_sweep", "auto")).strip().lower()
    if not spec or spec == "none":
        return []
    if spec == "auto":
        return hw.thread_points()
    counts = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            counts.append(int(part))
        except ValueError:
            continue  # skip non-numeric fields
    seen = set()
    return [count for count in counts if not (count in seen or seen.add(count))]


class ThreadSweepMixin:
    """Thread-count sweep with per-thread RSS scaling for threaded tools.

    A runner mixes this in next to ToolRunner and calls run_sweep() for each benchmark phase whose command accepts a thread count through ``-t``, ``--threads`` or similar.
    The phase's command runs under hyperfine at every count in the sweep, which records mean wall time and peak RSS per count along with the derived speedup and RSS-scaling ratios.

    Metric keys emitted for a phase ``phase`` such as "build" or "query"
        {phase}_s_t{N}          mean seconds at N threads
        {phase}_rss_t{N}        peak RSS in MB at N threads
        {phase}_speedup_t{N}    s(t1) / s(tN) for N above the first count
        {phase}_rss_scale       rss(tN) / rss(t1) across the sweep
        {phase}_sweep_rows      ordered list of {threads, s, rss, speedup}
    """

    def sweep_threads(self):
        """The concrete thread counts to sweep, resolved from ctx."""
        return sweep_thread_counts(self.ctx)

    def sweep_default(self, phase):
        """Thread count used for that phase's headline (parity) metrics."""
        key = {
            "build": "threads_build",
            "index": "threads_build",
            "align": "threads_build",
        }.get(phase, "threads_query")
        return int(self.ctx.get(key, 1))

    def run_sweep(
        self,
        phase,
        cmd_factory,
        gate=None,
        timeout=None,
        ignore_failure=False,
        **hf_kwargs,
    ):
        """Run ``cmd_factory(threads)`` under hyperfine at each swept count.

        Parameters
        phase : str
            Metric prefix, "build" or "query".
        cmd_factory : callable(threads:int) -> list[str] | None
            Builds the benchmark command for a thread count.
            Returning None skips that count, as for a tool that clamps thread counts.
        gate : callable(CompletedProcess) -> bool | None
            Optional correctness gate, run once and untimed at the first swept count before any timing.
            A False result aborts the sweep with a RuntimeError, and None runs no gate.
        timeout, ignore_failure
            Forwarded to hyperfine.
        **hf_kwargs
            Further keyword arguments for run_hyperfine, such as env and pin_cpus.
            ``env`` also applies to the gate run.

        Returns
        (flat, rows, headline_hf, headline_threads)
            flat              per-thread and scaling metrics, as the class doc lists
            rows              ordered [{threads, s, rss, speedup}, ...]
            headline_hf       hyperfine dict of the sweep point at the runner's default thread count, or None when the sweep is disabled and the caller falls back to a normal run at the default count
            headline_threads  thread count of headline_hf, or None
        """
        sweep = self.sweep_threads()
        flat = {}
        rows = []
        headline_hf = None
        headline_threads = None
        if not sweep:
            return flat, rows, headline_hf, headline_threads

        default = self.sweep_default(phase)
        base_time = None
        verbose = self.ctx.get("verbose", False)
        for i, threads in enumerate(sweep):
            cmd = cmd_factory(threads)
            if cmd is None:
                continue
            # The first thread count runs once outside hyperfine to check its output, so a sweep never times a broken command.
            if i == 0 and gate is not None:
                gproc, _ = run_cmd(
                    cmd,
                    timeout=timeout,
                    verbose=verbose,
                    env=hf_kwargs.get("env"),
                    cwd=hf_kwargs.get("cwd"),
                )
                if gproc.returncode != 0 or gate(gproc) is False:
                    err = (gproc.stderr or gproc.stdout or "")[:DETAIL_CHARS]
                    raise RuntimeError(
                        f"{self.tool_name()} {phase} thread-sweep gate failed at t={threads}: {err}"
                    )
            hyperfine_res = self.run_hyperfine(
                cmd, timeout=timeout, ignore_failure=ignore_failure, **hf_kwargs
            )
            if "error" in hyperfine_res or hyperfine_res.get("mean_s") is None:
                continue
            secs = hyperfine_res["mean_s"]
            rss = hyperfine_res.get("max_rss_mb")
            # Speedups are relative to the first thread count that ran, which is normally one thread.
            if base_time is None:
                base_time = secs
                speedup = 1.0
            else:
                speedup = round(base_time / secs, 3) if secs else None
            rows.append(
                {
                    "threads": threads,
                    "s": round(secs, 4),
                    "rss": rss,
                    "speedup": speedup,
                }
            )
            flat[f"{phase}_s_t{threads}"] = round(secs, 4)
            flat[f"{phase}_rss_t{threads}"] = rss
            if threads > 1 and secs:
                flat[f"{phase}_speedup_t{threads}"] = round(base_time / secs, 3)
            # The measurement at the default thread count also becomes the tool's headline figure.
            if threads == default:
                headline_hf = hyperfine_res
                headline_threads = threads

        if len(rows) >= MIN_SCALING_ROWS:
            first_rss = rows[0].get("rss")
            last_rss = rows[-1].get("rss")
            if first_rss and last_rss:
                flat[f"{phase}_rss_scale"] = round(last_rss / first_rss, 3)
        flat[f"{phase}_sweep_rows"] = rows
        return flat, rows, headline_hf, headline_threads


class ToolRunner:
    """Base class for all tool benchmark runners.

    Every runner subclass implements
        tool_name()    human-readable name for output tables
        binary_path()  filesystem path to the tool binary, or None
        run()          execute the benchmark phases and return a result dict

    Attributes
    ds : dict
        Dataset configuration with the unitigs path, queries path and k.
    ctx : dict
        Runtime context with verbosity, thread counts and hyperfine settings.
    """

    def __init__(self, dataset, ctx):
        self.dataset = dataset or {}
        self.ctx = ctx or {}

    @classmethod
    def tool_name(cls):
        """Human-readable tool name used in output tables."""
        raise NotImplementedError

    @classmethod
    def binary_path(cls):
        """Filesystem path of the tool binary, or None."""
        raise NotImplementedError

    @classmethod
    def available(cls):
        """Return True if the tool binary exists and is executable."""
        return bool(cls.binary_path())

    def tmpdir(self, prefix):
        """Create a temporary directory for this runner's files.

        The directory is cleaned up after run() completes via the cleanup() method in the finally block.
        """
        parent = Path(self.ctx.get("tmp_dir", "/tmp/tool_bench"))
        parent.mkdir(parents=True, exist_ok=True)
        return tempfile.mkdtemp(prefix=f"{prefix}_", dir=parent)

    def cleanup(self, directory):
        """Remove a temporary directory and all its contents."""
        shutil.rmtree(directory, ignore_errors=True)

    @property
    def _hf_runs(self):
        """Number of hyperfine measurement runs from the runtime context."""
        return self.ctx.get("runs", DEFAULT_HYPERFINE_RUNS)

    @property
    def _hf_warmup(self):
        """Number of hyperfine warmup runs from the runtime context."""
        return self.ctx.get("warmup", 1)

    def run_hyperfine(self, cmd, **kwargs):
        """Wrap run_hyperfine() with this runner's default settings.

        Passes through the configured run count, warmup count, verbosity and CPU-pinning automatically so subclasses do not need to repeat them.
        """
        kwargs.setdefault("runs", self._hf_runs)
        kwargs.setdefault("warmup", self._hf_warmup)
        kwargs.setdefault("verbose", self.ctx.get("verbose", False))
        kwargs.setdefault("pin_cpus", self.ctx.get("pin_cpus"))
        return run_hyperfine(cmd, **kwargs)

    def run(self):
        """Execute all benchmark phases for this tool.

        Returns a flat dict of metric name to value that is merged into the output JSON and CSV.
        On error the dict should contain an "error" key with a descriptive message.
        """
        raise NotImplementedError


def ensure_uncompressed(path, temp_dir, strip_blank_lines=False):
    """Decompress a .gz file and return the path of the plain copy.

    A file that is not gzipped is returned unchanged.
    Decompressed copies go to temp_dir, so they disappear when the runner finishes.

    Decompression is kept apart from tool execution because only some tools, such as sshash build, minimap2 and kmercamel, read .gz natively, and build times should not include gzip decompression.
    """
    src_path = Path(path)
    if not str(src_path).endswith(".gz"):
        return str(src_path)
    out = Path(temp_dir) / src_path.stem
    if not out.exists():
        import gzip

        with gzip.open(src_path, "rt") as f_in, open(out, "w") as f_out:
            if strip_blank_lines:
                last = ""
                for line in f_in:
                    if line.strip():
                        f_out.write(line)
                        last = line
                # Several tools (sbwt, minimap2) need a trailing newline, which the source archive may omit, so the rewritten file always ends with one.
                if last and not last.endswith("\n"):
                    f_out.write("\n")
            else:
                f_out.write(f_in.read())
    return str(out)


def _canonical_kmer_int(kmer):
    """Integer encoding of the canonical k-mer, the smaller of forward and reverse complement.

    Each base takes 2 bits (A=0, C=1, G=2, T=3).
    The forward and reverse-complement encodings are both high-aligned, with the first base in the most significant bits.
    """
    base_code = {"A": 0, "C": 1, "G": 2, "T": 3}
    fwd = 0
    rev = 0
    for base_char in kmer:
        code = base_code.get(base_char.upper(), 0)
        fwd = (fwd << BITS_PER_BASE) | code
    for base_char in reversed(kmer):
        code = BASE_MASK - base_code.get(base_char.upper(), 0)
        rev = (rev << BITS_PER_BASE) | code
    return min(fwd, rev)


def _iter_fasta_seqs(path):
    """Yield sequence strings from a FASTA file, gzipped or not."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    seq = []
    with opener(str(path), "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    yield "".join(seq)
                    seq.clear()
            else:
                seq.append(line)
        if seq:
            yield "".join(seq)


def count_kmers(path, kmer_len, distinct=False):
    """Count the k-mers of a FASTA file.

    Parameters
    path : str
        FASTA file, optionally gzipped.
    k : int
        k-mer length.
    distinct : bool
        True counts distinct *canonical* k-mers, and False the total count (the sum of L - k + 1).

    k-mers holding a base other than A/C/G/T are skipped, as the indexed tools skip them.

    The counts fill the ``distinct_kmers`` and ``input_kmers`` metrics, which some tools do not report themselves but which matter when comparing tools of one category.
    """
    total = 0
    seen = set()
    valid = {"A", "C", "G", "T"}
    for seq in _iter_fasta_seqs(path):
        num_windows = len(seq) - kmer_len + 1
        if num_windows <= 0:
            continue
        if distinct:
            for i in range(num_windows):
                kmer = seq[i : i + kmer_len]
                if kmer and all(base in valid for base in kmer):
                    seen.add(_canonical_kmer_int(kmer))
        else:
            total += num_windows
    if distinct:
        return len(seen) if seen else None
    return total if total > 0 else None


def count_invalid_windows(path, kmer_len):
    """Count the k-mer windows of a FASTA/FASTQ query file, plain or .gz, that hold a non-ACGT base.

    SSHash and optik report these as num_invalid_kmers.
    Tools that print only present or absent per window, SBWT and FMSI, count them as absent, so their negative count becomes total - positive - invalid to keep one meaning across the table."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    seqs = []
    with opener(str(path), "rt") as handle:
        first = handle.read(1)
        handle.seek(0)
        if first == "@":
            for i, line in enumerate(handle):
                if i % FASTQ_LINES_PER_RECORD == FASTQ_SEQ_LINE:
                    seqs.append(line.strip())
        else:
            seqs = list(_iter_fasta_seqs(path))
    bad = 0
    for seq in seqs:
        upper_seq = seq.upper()
        last_bad = NO_BAD_BASE  # index of the most recent non-ACGT base
        for i, base in enumerate(upper_seq):
            if base not in "ACGT":
                last_bad = i
            if i >= kmer_len - 1 and last_bad > i - kmer_len:
                bad += 1
    return bad


def write_acgt_segments(path, kmer_len, out_path):
    """Write each maximal ACGT run of at least k bases in a query file as its own FASTA record.

    The input is FASTA or FASTQ, plain or .gz.
    Returns (valid windows, invalid windows), so every tool is asked about the same k-mers and windows containing N are accounted separately.

    This serves tools whose per-record answer stops being one value per k-mer window once a record contains a non-ACGT base.
    FMSI, for example, prints one character per base for such records, which cannot be attributed to windows."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt") as handle:
        first = handle.read(1)
        handle.seek(0)
        if first == "@":
            seqs = [line.strip() for i, line in enumerate(handle) if i % FASTQ_LINES_PER_RECORD == FASTQ_SEQ_LINE]
        else:
            seqs = list(_iter_fasta_seqs(path))
    valid = total = 0
    with open(out_path, "w") as out:
        rec_idx = 0
        for seq in seqs:
            total += max(0, len(seq) - kmer_len + 1)
            for run in re.findall(r"[ACGT]+", seq.upper()):
                if len(run) >= kmer_len:
                    out.write(f">s{rec_idx}\n{run}\n")
                    rec_idx += 1
                    valid += len(run) - kmer_len + 1
    return valid, total - valid


POSITIVE_READ_COUNT = 20000
# length and seed of the positive query reads
POSITIVE_READ_LEN = 150
POSITIVE_READ_SEED = 7
# probability that a positive read is reverse-complemented
REVCOMP_PROB = 0.5
# seed of the random k-mer workload
RANDOM_KMER_SEED = 42
# bases a random k-mer draws from
NUM_BASES = 4


def generate_positive_reads(
    fasta_path, output_path, count=POSITIVE_READ_COUNT, read_len=POSITIVE_READ_LEN, seed=POSITIVE_READ_SEED
):
    """Write a high-hit query FASTA made of substrings of the indexed sequences.

    Reads are sampled uniformly from sequences at least `read_len` long, with shorter ones taken whole, and half are reverse-complemented.
    Nearly every query k-mer is then positive, and consecutive k-mers share locality like real reads from the indexed genome.
    Query sets dominated by negatives would only time cheap rejections, not lookups.
    """
    import random as random_mod

    rng = random_mod.Random(seed)
    seqs = []
    parts = []
    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if parts:
                    seqs.append("".join(parts))
                parts = []
            elif line:
                parts.append(line.upper())
    if parts:
        seqs.append("".join(parts))
    seqs = [seq for seq in seqs if seq]
    if not seqs:
        raise ValueError(f"no sequences in {fasta_path}")
    weights = [len(seq) for seq in seqs]
    comp = str.maketrans("ACGT", "TGCA")
    with open(output_path, "w") as out:
        for i in range(count):
            read = rng.choices(seqs, weights=weights)[0]
            if len(read) > read_len:
                start = rng.randrange(len(read) - read_len + 1)
                read = read[start : start + read_len]
            if rng.random() < REVCOMP_PROB:
                read = read.translate(comp)[::-1]
            out.write(f">p{i}\n{read}\n")
    return output_path


def generate_random_kmers(output_path, kmer_len, count, seed=RANDOM_KMER_SEED):
    """Write a FASTA file of uniformly random, independent k-mers.

    Each k-mer is drawn uniformly from the 4^k possible sequences.
    Being independent, consecutive k-mers share no cache locality, which is how k-mer membership papers usually measure worst-case lookup throughput.

    Bases come from Lemire's unbiased range reduction (Lemire 2019), so the draws carry no modulo bias.

    Parameters
    output_path : str
        Path of the FASTA file to write.
    k : int
        Length of each k-mer.
    count : int
        Number of k-mers to generate.
    seed : int
        Random seed, for reproducibility.
    """
    rng = random.Random(seed)
    bases = ["A", "C", "G", "T"]
    with open(output_path, "w") as handle:
        for i in range(count):
            kmer = "".join(bases[lemire_randbelow(rng, NUM_BASES)] for _ in range(kmer_len))
            handle.write(f">{i}\n{kmer}\n")


def lemire_randbelow(rng, bound):
    """Unbiased random integer in [0, n) by Lemire's method (Lemire 2019).

    Plain ``x % n`` carries modulo bias whenever n does not divide the PRNG range.
    Lemire's scheme maps the 64-bit draw into [0, n) uniformly and redraws only on a rare threshold collision.
    A power-of-two n reduces to a masked draw.
    """
    if bound <= 0:
        raise ValueError("n must be >= 1")
    if bound & (bound - 1) == 0:
        return rng.getrandbits(bound.bit_length() - 1)
    threshold = (-bound) % bound
    mask64 = (1 << WORD_BITS) - 1
    while True:
        draw = rng.getrandbits(WORD_BITS)
        product = draw * bound
        if (product & mask64) >= threshold:
            return product >> WORD_BITS


def gate_query_report(text, tool, phase):
    """Correctness gate on a tool's query report.

    Before any timing phase may run, the tool's own accounting must add up, meaning reported queries equal positives plus negatives.
    Otherwise a descriptive ValueError is raised, so a runner fails fast with a clear message instead of timing a broken query path.
    Returns the positive ratio when the report parses, and None when the tool prints no count report, in which case only the exit status is checked upstream.
    """

    def num(key):
        """Integer value printed after `key` in the report, or None."""
        for line in text.split("\n"):
            line = line.strip()
            if line.startswith((key + " ", key + "\t")):
                after = line[len(key) :].lstrip("= \t:") or "0"
                tok = after.split()[0]
                return int(tok.split("/")[0])
        return None

    total = num("num_kmers")
    pos = num("num_positive_kmers")
    neg = num("num_negative_kmers")
    invalid = num("num_invalid_kmers")
    if total is None:
        return None
    if pos is None and neg is None:
        return None
    pos = pos if pos is not None else 0
    neg = neg if neg is not None else 0
    invalid = invalid if invalid is not None else 0
    if pos + neg + invalid != total:
        raise ValueError(
            f"{tool} {phase}: query accounting mismatch "
            f"(num_kmers={total} != positive {pos} + negative {neg}"
            f"{f' + invalid {invalid}' if invalid else ''}); "
            "refusing to time an inconsistent query report"
        )
    return pos / total if total else None


def unwrap_fasta(path, out_dir=None):
    """Join multi-line FASTA sequences into single-line records.

    Query parsers that expect one sequence per line need this.
    Without --multiline, SSHash's query would treat alternate lines of a wrapped record as headers.
    The output goes to out_dir (by default next to the input) with a _uw_ prefix, and gzipped input is read transparently.
    """
    in_path = Path(path)
    if not in_path.exists() and in_path.with_suffix(in_path.suffix + ".gz").exists():
        return unwrap_fasta(str(in_path.with_suffix(in_path.suffix + ".gz")), out_dir)

    opener = __import__("gzip").open if str(path).endswith(".gz") else open
    records = []
    with opener(path, "rt") as handle:
        buffer = []
        for line in handle:
            if line.startswith(">"):
                if buffer:
                    records.append("".join(buffer))
                buffer = [line]
            else:
                buffer.append(line.strip())
        if buffer:
            records.append("".join(buffer))

    out_path = str(Path(out_dir or Path(path).parent) / ("_uw_" + Path(path).name))
    with open(out_path, "w") as handle:
        for rec in records:
            if rec.startswith(">"):
                # Each record becomes a header line followed by its joined sequence.
                handle.write(rec + "\n")
    return out_path
