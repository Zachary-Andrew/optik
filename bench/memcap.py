"""Run a command under a resident-memory cap on any Linux.

The command starts in its own process group, and a watchdog thread sums the resident memory (/proc/PID/statm) of every process in that group every 0.2 s, killing the whole group once the sum passes the cap.
It needs no systemd, cgroup delegation or root, and it caps resident memory rather than address space (RLIMIT_AS would kill tools that reserve far more than they touch, SSHash among them).
A run killed this way has `.capped = True` and returncode -9.
Each entry point returns the following::

  memcap.run(cmd, cap_bytes, **subprocess_kwargs) -> CompletedProcess (+ .capped)
  memcap.popen(cmd, cap_bytes, **popen_kwargs)    -> Popen (+ .capped after exit)
  memcap.parse_size("26G")                        -> bytes
  memcap.warm(paths)                              -> bytes read into the page cache
  python3 memcap.py --warm FILE...                (the same as a command, for hyperfine --prepare)
"""

import os
import signal
import subprocess
import threading
import time

from units import BYTES_PER_GIB, BYTES_PER_KIB, BYTES_PER_MIB, BYTES_PER_TIB

# Inputs are pre-read in chunks of this many bytes, and only when they take at most 1/WARM_MAX_SHARE of the available memory.
READ_CHUNK_BYTES = 1 << 24
WARM_MAX_SHARE = 2

PAGE = os.sysconf("SC_PAGE_SIZE")
POLL_S = 0.2


# /proc/<pid>/stat puts the command name in parentheses, and pgrp is the third field after it (proc(5)).
STAT_COMM_END = ") "
STAT_PGRP = 2


def parse_size(size):
    """Bytes of a size such as '26G', '512M' or '1048576', with '' or None giving 0 for no cap."""
    if not size:
        return 0
    size = str(size).strip().upper()
    mult = {"K": BYTES_PER_KIB, "M": BYTES_PER_MIB, "G": BYTES_PER_GIB, "T": BYTES_PER_TIB}
    return int(float(size[:-1]) * mult[size[-1]]) if size[-1] in mult else int(size)


def group_rss(pgid):
    """Resident bytes of every process whose process group is pgid."""
    tot = 0
    for pid_dir in os.listdir("/proc"):
        if not pid_dir.isdigit():
            continue
        try:
            with open(f"/proc/{pid_dir}/stat") as handle:
                stat_line = handle.read()
            # Fields after the parenthesised command name are state, ppid, pgrp and so on.
            if int(stat_line[stat_line.rindex(STAT_COMM_END) + len(STAT_COMM_END) :].split()[STAT_PGRP]) != pgid:
                continue
            with open(f"/proc/{pid_dir}/statm") as handle:
                tot += int(handle.read().split()[1]) * PAGE
        except (OSError, ValueError, IndexError):
            continue  # the process exited while being read
    return tot


class _Watch(threading.Thread):
    """Watchdog thread that kills a process group once its resident memory passes the cap."""

    def __init__(self, proc, cap):
        super().__init__(daemon=True)
        self.proc, self.cap, self.capped = proc, cap, False

    def run(self):
        while self.proc.poll() is None:
            if group_rss(self.proc.pid) > self.cap:
                self.capped = True
                try:
                    os.killpg(self.proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass  # the group already exited
                return
            time.sleep(POLL_S)


def popen(cmd, cap=0, **kwargs):
    """subprocess.Popen in its own process group, watched against `cap` resident bytes, 0 meaning no cap.
    After the process exits, `.capped` tells whether the watchdog killed it."""
    kwargs.setdefault("start_new_session", True)
    proc = subprocess.Popen(cmd, **kwargs)
    proc._watch = _Watch(proc, cap) if cap else None
    if proc._watch:
        proc._watch.start()
    return proc


def capped(proc):
    """True when the watchdog of process p killed it at the cap."""
    return bool(getattr(proc, "_watch", None) and proc._watch.capped)


def run(cmd, cap=0, timeout=None, check=False, input=None, **kwargs):
    """subprocess.run under the cap, returning a CompletedProcess that carries `.capped`.

    Ctrl-C or a timeout kills the whole group, not only its leader."""
    if kwargs.pop("capture_output", False):  # as subprocess.run handles it
        kwargs["stdout"] = kwargs["stderr"] = subprocess.PIPE
    if input is not None and "stdin" not in kwargs:
        kwargs["stdin"] = subprocess.PIPE
    proc = popen(cmd, cap, **kwargs)
    try:
        out, err = proc.communicate(input=input, timeout=timeout)
    except BaseException:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass  # the group already exited
        proc.wait()
        raise
    if proc._watch:
        proc._watch.join()
    result = subprocess.CompletedProcess(proc.args, proc.returncode, out, err)
    result.capped = capped(proc)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, proc.args, out, err)
    return result


def warm(paths):
    """Read input files into the page cache before the tools run, so no tool is timed reading from disk while the others read from memory.

    Inputs larger than half the available memory cannot stay cached and are left alone, so every tool streams them from disk alike.
    Returns the bytes read."""
    paths = [path for path in dict.fromkeys(paths) if path and os.path.isfile(path)]
    tot = sum(os.path.getsize(path) for path in paths)
    with open("/proc/meminfo") as handle:
        avail = (
            next(int(line.split()[1]) for line in handle if line.startswith("MemAvailable:")) * BYTES_PER_KIB
        )
    if not paths or tot > avail // WARM_MAX_SHARE:
        return 0
    for path in paths:
        with open(path, "rb") as handle:
            while handle.read(READ_CHUNK_BYTES):
                pass
    return tot


if __name__ == "__main__":
    # `python3 bench/memcap.py --warm FILE...` runs warm() as a command for hyperfine --prepare, so every timed run begins with identical cached inputs.
    import sys

    match sys.argv[1:]:
        case ["--warm", *paths]:
            warm(paths)
        case _:
            sys.exit("usage: memcap.py --warm FILE...")
