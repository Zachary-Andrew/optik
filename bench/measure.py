#!/usr/bin/env python3
"""The single timing path of every optik benchmark.

measure(cmd, warmup, runs, ...) runs `cmd` under hyperfine (-N, no shell) with the --warmup W and --runs R of the bench's command line, all under a resident-memory cap (memcap.py, applied to the whole process group on any Linux without root).
Each run goes through rusage_exec (bench/cpp/rusage_exec.c) in append mode, so every run, warmups first, leaves one line of wall time, peak RSS and block I/O.
The first W lines are dropped and exactly R must remain.
Statistics cover the R measured runs only.

  wall      hyperfine's mean, standard deviation, median, min and max
  user/sys  hyperfine's mean user and system CPU seconds
  rss_mb    peak resident set, median and max over the runs
  written / read   block I/O bytes per run (median)

Before every run, warmups included, hyperfine's --prepare reads the command's input files into the page cache (memcap.py --warm), so no run and no tool is timed on a cold read that another skips.
Inputs larger than half the available memory cannot stay cached and are streamed by every run alike.

The command's stdout and stderr go to files (rusage_exec -o/-e) truncated on each run, so afterwards they hold the last measured run's output, which the bench's correctness gate reads.
A run that fails (non-zero exit, or killed at the cap) ends the measurement, with Result.succeeded False, .capped telling whether the cap killed it and .err holding its stderr.
"""

import json
import os
import shlex
import subprocess
import sys
from dataclasses import (
    dataclass,
    field,
)
from typing import NamedTuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import statistics

import memcap
from diagnostics import LOG_TAIL_CHARS
from units import KIB_PER_MIB

# exit code of a Result whose failure is the bench's own (a count mismatch), not the command's
NO_EXIT_CODE = -1


class RusageLine(NamedTuple):
    """One line of rusage_exec output, describing a single run of the command."""

    exit_code: int
    wall_s: float
    peak_rss_kib: int
    bytes_written: int | None
    bytes_read: int | None

    @classmethod
    def parse(cls, line):
        """Parse "code wall_s rss_kib [written read]".

        The two block-I/O fields are None when rusage_exec omits them.
        """
        code, wall_s, rss_kib, *io_fields = line.split()
        # Newer rusage_exec builds append bytes written and bytes read, and older lines stop after RSS.
        written, read = (int(io_fields[0]), int(io_fields[1])) if len(io_fields) >= 2 else (None, None)
        return cls(int(code), float(wall_s), int(rss_kib), written, read)


RUSAGE_SRC = os.path.join(HERE, "cpp", "rusage_exec.c")
RUSAGE_EXEC = os.path.join(HERE, "..", "temp", "measure", "rusage_exec")


def rusage_exec():
    """The rusage_exec binary, rebuilt when its source is newer."""
    exe = os.path.abspath(RUSAGE_EXEC)
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(RUSAGE_SRC):
        os.makedirs(os.path.dirname(exe), exist_ok=True)
        subprocess.run(["gcc", "-O2", "-o", exe + ".part", RUSAGE_SRC], check=True)
        os.replace(exe + ".part", exe)
    return exe


@dataclass
class Result:
    """Outcome and statistics of one measure() call."""

    succeeded: bool
    capped: bool = False
    code: int = 0
    warmup: int = 0
    runs: int = 0
    wall: dict = field(default_factory=dict)  # mean, sd, median, min, max (s)
    user_s: float = None
    sys_s: float = None
    rss_mb: dict = field(default_factory=dict)  # median, max
    written: int = None  # bytes per run (median)
    read: int = None
    out: str = ""  # stdout file of the last run
    err: str = ""  # stderr text of the last run
    per_run: list = field(
        default_factory=list
    )  # [(wall s, rss kb, written, read)] of the kept runs

    def wall_row(self):
        """The wall-time cell of a results row, {mean, sd, median, min, max, runs, warmup}, which the benches print as mean±sd."""
        return dict(self.wall, runs=self.runs, warmup=self.warmup) if self.succeeded else None

    def rss(self):
        """Peak resident MB over the measured runs."""
        return self.rss_mb.get("max")

    def row(self, prefix=""):
        """Flat fields for a results row."""
        col_prefix = prefix
        fields = {
            f"{col_prefix}wall_s": self.wall.get("mean"),
            f"{col_prefix}wall_sd": self.wall.get("sd"),
            f"{col_prefix}wall_median": self.wall.get("median"),
            f"{col_prefix}wall_min": self.wall.get("min"),
            f"{col_prefix}wall_max": self.wall.get("max"),
            f"{col_prefix}user_s": self.user_s,
            f"{col_prefix}sys_s": self.sys_s,
            f"{col_prefix}rss_mb": self.rss_mb.get("max"),
            f"{col_prefix}rss_mb_median": self.rss_mb.get("median"),
            f"{col_prefix}warmup": self.warmup,
            f"{col_prefix}runs": self.runs,
        }
        if self.written is not None:
            fields[f"{col_prefix}written_bytes"] = self.written
            fields[f"{col_prefix}read_bytes"] = self.read
        return fields


def _read(path):
    """Text of the file at path, or an empty string when it is missing."""
    if not os.path.exists(path):
        return ""
    with open(path, errors="replace") as handle:
        return handle.read()


# Seconds per unit of a hyperfine 2 time value
_SECONDS = {"second": 1.0, "millisecond": 1e-3, "microsecond": 1e-6, "nanosecond": 1e-9}


class AppendLog:
    """A command log that appends each write to its file and closes it again, so no handle stays open."""

    def __init__(self, path):
        self.path = path

    def write(self, text):
        """Append `text` to the log file."""
        with open(self.path, "a") as handle:
            handle.write(text)

    def flush(self):
        """Nothing to flush, since every write closes the file."""


def hyperfine_result(doc):
    """The first result of a hyperfine JSON export, with the field names of hyperfine 1.x.

    hyperfine 1.x lists `times`, `mean`, `stddev`, `median`, `min`, `max`, `user`, `system`, `exit_codes` and `memory_usage_byte` per result.
    hyperfine 2 instead gives each run in `measurements` and the statistics in `summary`, so those are converted here and both versions read alike.
    """
    first_result = doc["results"][0]
    if "times" in first_result:
        return first_result
    runs = first_result.get("measurements", [])
    summ = first_result.get("summary", {})

    def sec(time_val):
        """A hyperfine 2 time value in seconds."""
        return time_val["value"] * _SECONDS.get(time_val.get("unit", "second"), 1.0)

    wall = summ.get("time_wall_clock", {})
    scale = _SECONDS.get(wall.get("unit", "second"), 1.0)
    out = {
        "command": first_result.get("command"),
        "times": [sec(measurement["time_wall_clock"]) for measurement in runs],
        "exit_codes": [measurement.get("exit_code", 0) for measurement in runs],
        "memory_usage_byte": [
            int(measurement["memory_peak_resident"]["value"])
            for measurement in runs
            if "memory_peak_resident" in measurement
        ],
    }
    for stat in ("mean", "median", "min", "max"):
        out[stat] = wall[stat] * scale if stat in wall else None
    out["stddev"] = (
        wall["stddev"] * scale if len(runs) > 1 and "stddev" in wall else None
    )
    for stat, src in (("user", "time_user"), ("system", "time_system")):
        stats = summ.get(src, {})
        out[stat] = (
            stats["mean"] * _SECONDS.get(stats.get("unit", "second"), 1.0)
            if "mean" in stats
            else None
        )
    return out


def _median(values):
    """Median of xs, the mean of the two middle values for an even count."""
    return statistics.median(values)


def measure(
    cmd, warmup, runs, cap=0, inputs=(), scratch=None, cwd=None, tag="cmd", log=None
):
    """Time the argument list `cmd` with W warmups and R measured runs, as the module docstring describes.

    `cap` gives the resident bytes for the whole run, or 0 for none, `inputs` the files pre-read before every run, and `scratch` a directory for the per-run files (default <cwd or .>/.measure).
    `log` is a file object that receives the full hyperfine command line."""
    if runs < 1:
        raise ValueError("runs must be >= 1")
    scratch = os.path.abspath(scratch or os.path.join(cwd or ".", ".measure"))
    os.makedirs(scratch, exist_ok=True)
    res_f, out_f, err_f, json_path = (
        os.path.join(scratch, name)
        for name in ("rusage.txt", "stdout.txt", "stderr.txt", "hyperfine.json")
    )
    # Results left by an earlier call would be read as this one's, so they are removed first.
    for path in (res_f, json_path):
        if os.path.exists(path):
            os.remove(path)
    # rusage_exec appends one line of exit code, wall time, peak RSS and I/O per run, while hyperfine times the wrapper as a whole.
    inner = [rusage_exec(), "-a", "-o", out_f, "-e", err_f, "0", res_f] + [
        str(arg) for arg in cmd
    ]
    hyperfine_cmd = [
        "hyperfine",
        "-N",
        "--warmup",
        str(warmup),
        "--runs",
        str(runs),
        "--export-json",
        json_path,
        "--command-name",
        tag,
        "--style",
        "none",
    ]
    ins = [os.path.abspath(os.path.join(cwd or ".", path)) for path in inputs if path]
    ins = [path for path in ins if os.path.isfile(path)]
    # The inputs are read back into the page cache before every run, so each run starts from the same cache state.
    if ins:
        hyperfine_cmd += [
            # Before each run, --prepare reads the inputs into the page cache so every run starts equally warm.
            "--prepare",
            shlex.join(
                [sys.executable, os.path.join(HERE, "memcap.py"), "--warm", *ins]
            ),
        ]
    hyperfine_cmd.append(shlex.join(inner))
    if log:
        log.write(f"[{tag}] {shlex.join(hyperfine_cmd)}\n")
        log.flush()
    proc = memcap.run(
        hyperfine_cmd, cap, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    err = _read(err_f)
    lines = [RusageLine.parse(line) for line in _read(res_f).splitlines() if line.strip()]
    if proc.returncode != 0 or proc.capped:
        # The rusage wrapper knows the tool's own exit code. hyperfine's code only says that some run failed.
        code = lines[-1].exit_code if lines else proc.returncode
        return Result(
            succeeded=False,
            capped=bool(proc.capped),
            code=code,
            warmup=warmup,
            runs=runs,
            err=err or (proc.stderr or "")[-LOG_TAIL_CHARS:],
            out=out_f,
        )
    # One rusage line per run is expected, so any other count means a run was lost or the file was shared.
    if len(lines) != warmup + runs:
        return Result(
            succeeded=False,
            code=NO_EXIT_CODE,
            warmup=warmup,
            runs=runs,
            err=f"measure: {len(lines)} rusage lines for {warmup} warmups + {runs} runs",
        )
    # Warmup runs are timed by hyperfine but left out of the statistics.
    kept = lines[warmup:]
    hyperfine = hyperfine_result(json.loads(_read(json_path)))
    times = hyperfine["times"]
    if len(times) != runs:
        return Result(
            succeeded=False,
            code=NO_EXIT_CODE,
            warmup=warmup,
            runs=runs,
            err=f"measure: hyperfine reported {len(times)} runs, {runs} asked",
        )
    # Wall time comes from hyperfine and memory from the rusage wrapper, both from the same runs.
    rss = [run.peak_rss_kib / KIB_PER_MIB for run in kept]
    written = [run.bytes_written for run in kept if run.bytes_written is not None]
    read_bytes = [run.bytes_read for run in kept if run.bytes_read is not None]
    return Result(
        succeeded=True,
        code=0,
        warmup=warmup,
        runs=runs,
        wall={
            "mean": hyperfine["mean"],
            "sd": hyperfine["stddev"] if runs > 1 else None,
            "median": hyperfine["median"],
            "min": hyperfine["min"],
            "max": hyperfine["max"],
        },
        user_s=hyperfine.get("user"),
        sys_s=hyperfine.get("system"),
        rss_mb={"median": _median(rss), "max": max(rss)},
        written=int(_median(written)) if written else None,
        read=int(_median(read_bytes)) if read_bytes else None,
        out=out_f,
        err=err,
        per_run=[
            (run.wall_s, run.peak_rss_kib, run.bytes_written, run.bytes_read)
            for run in kept
        ],
    )


if __name__ == "__main__":
    # python3 bench/measure.py W R cmd... gives a quick check of the path.
    warmup, runs, *cmd = sys.argv[1:]
    result = measure(cmd, int(warmup), int(runs))
    print(json.dumps(result.__dict__, indent=1, default=str))
