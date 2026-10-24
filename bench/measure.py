#!/usr/bin/env python3
"""measure.py - the single timing path of every optik benchmark.

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
A run that fails (non-zero exit, or killed at the cap) ends the measurement, with Result.ok False, .capped telling whether the cap killed it and .err holding its stderr.
"""

import json
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import memcap

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

    ok: bool
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
        return dict(self.wall, runs=self.runs, warmup=self.warmup) if self.ok else None

    def rss(self):
        """Peak resident MB over the measured runs."""
        return self.rss_mb.get("max")

    def row(self, prefix=""):
        """Flat fields for a results row."""
        p = prefix
        d = {
            f"{p}wall_s": self.wall.get("mean"),
            f"{p}wall_sd": self.wall.get("sd"),
            f"{p}wall_median": self.wall.get("median"),
            f"{p}wall_min": self.wall.get("min"),
            f"{p}wall_max": self.wall.get("max"),
            f"{p}user_s": self.user_s,
            f"{p}sys_s": self.sys_s,
            f"{p}rss_mb": self.rss_mb.get("max"),
            f"{p}rss_mb_median": self.rss_mb.get("median"),
            f"{p}warmup": self.warmup,
            f"{p}runs": self.runs,
        }
        if self.written is not None:
            d[f"{p}written_bytes"] = self.written
            d[f"{p}read_bytes"] = self.read
        return d


def _read(path):
    """Text of the file at path, or an empty string when it is missing."""
    if not os.path.exists(path):
        return ""
    with open(path, errors="replace") as f:
        return f.read()


# Seconds per unit of a hyperfine 2 time value
_SECONDS = {"second": 1.0, "millisecond": 1e-3, "microsecond": 1e-6, "nanosecond": 1e-9}


class AppendLog:
    """A command log that appends each write to its file and closes it again, so no handle stays open."""

    def __init__(self, path):
        self.path = path

    def write(self, text):
        """Append `text` to the log file."""
        with open(self.path, "a") as fh:
            fh.write(text)

    def flush(self):
        """Nothing to flush, since every write closes the file."""


def hyperfine_result(doc):
    """The first result of a hyperfine JSON export, with the field names of hyperfine 1.x.

    hyperfine 1.x lists `times`, `mean`, `stddev`, `median`, `min`, `max`, `user`, `system`, `exit_codes` and `memory_usage_byte` per result.
    hyperfine 2 instead gives each run in `measurements` and the statistics in `summary`, so those are converted here and both versions read alike.
    """
    r = doc["results"][0]
    if "times" in r:
        return r
    runs = r.get("measurements", [])
    summ = r.get("summary", {})

    def sec(v):
        """A hyperfine 2 time value in seconds."""
        return v["value"] * _SECONDS.get(v.get("unit", "second"), 1.0)

    wall = summ.get("time_wall_clock", {})
    scale = _SECONDS.get(wall.get("unit", "second"), 1.0)
    out = {
        "command": r.get("command"),
        "times": [sec(m["time_wall_clock"]) for m in runs],
        "exit_codes": [m.get("exit_code", 0) for m in runs],
        "memory_usage_byte": [
            int(m["memory_peak_resident"]["value"])
            for m in runs
            if "memory_peak_resident" in m
        ],
    }
    for k in ("mean", "median", "min", "max"):
        out[k] = wall[k] * scale if k in wall else None
    out["stddev"] = (
        wall["stddev"] * scale if len(runs) > 1 and "stddev" in wall else None
    )
    for k, src in (("user", "time_user"), ("system", "time_system")):
        s = summ.get(src, {})
        out[k] = (
            s["mean"] * _SECONDS.get(s.get("unit", "second"), 1.0)
            if "mean" in s
            else None
        )
    return out


def _median(xs):
    """Median of xs, the mean of the two middle values for an even count."""
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


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
    res_f, out_f, err_f, js = (
        os.path.join(scratch, x)
        for x in ("rusage.txt", "stdout.txt", "stderr.txt", "hyperfine.json")
    )
    for f in (res_f, js):
        if os.path.exists(f):
            os.remove(f)
    inner = [rusage_exec(), "-a", "-o", out_f, "-e", err_f, "0", res_f] + [
        str(c) for c in cmd
    ]
    hf = [
        "hyperfine",
        "-N",
        "--warmup",
        str(warmup),
        "--runs",
        str(runs),
        "--export-json",
        js,
        "--command-name",
        tag,
        "--style",
        "none",
    ]
    ins = [os.path.abspath(os.path.join(cwd or ".", x)) for x in inputs if x]
    ins = [x for x in ins if os.path.isfile(x)]
    if ins:
        hf += [
            "--prepare",
            shlex.join(
                [sys.executable, os.path.join(HERE, "memcap.py"), "--warm", *ins]
            ),
        ]
    hf.append(shlex.join(inner))
    if log:
        log.write(f"[{tag}] {shlex.join(hf)}\n")
        log.flush()
    p = memcap.run(
        hf, cap, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    err = _read(err_f)
    lines = [x.split() for x in _read(res_f).splitlines()]
    if p.returncode != 0 or p.capped:
        code = int(lines[-1][0]) if lines else p.returncode
        return Result(
            ok=False,
            capped=bool(p.capped),
            code=code,
            warmup=warmup,
            runs=runs,
            err=err or (p.stderr or "")[-2000:],
            out=out_f,
        )
    if len(lines) != warmup + runs:
        return Result(
            ok=False,
            code=-1,
            warmup=warmup,
            runs=runs,
            err=f"measure: {len(lines)} rusage lines for {warmup} warmups + {runs} runs",
        )
    kept = lines[warmup:]
    h = hyperfine_result(json.loads(_read(js)))
    times = h["times"]
    if len(times) != runs:
        return Result(
            ok=False,
            code=-1,
            warmup=warmup,
            runs=runs,
            err=f"measure: hyperfine reported {len(times)} runs, {runs} asked",
        )
    rss = [int(x[2]) / 1024 for x in kept]
    wr = [int(x[3]) for x in kept if len(x) >= 5]
    rd = [int(x[4]) for x in kept if len(x) >= 5]
    return Result(
        ok=True,
        code=0,
        warmup=warmup,
        runs=runs,
        wall={
            "mean": h["mean"],
            "sd": h["stddev"] if runs > 1 else None,
            "median": h["median"],
            "min": h["min"],
            "max": h["max"],
        },
        user_s=h.get("user"),
        sys_s=h.get("system"),
        rss_mb={"median": _median(rss), "max": max(rss)},
        written=int(_median(wr)) if wr else None,
        read=int(_median(rd)) if rd else None,
        out=out_f,
        err=err,
        per_run=[
            (
                float(x[1]),
                int(x[2]),
                int(x[3]) if len(x) >= 5 else None,
                int(x[4]) if len(x) >= 5 else None,
            )
            for x in kept
        ],
    )


if __name__ == "__main__":
    # python3 bench/measure.py W R cmd... gives a quick check of the path.
    w, r, *c = sys.argv[1:]
    m = measure(c, int(w), int(r))
    print(json.dumps(m.__dict__, indent=1, default=str))
