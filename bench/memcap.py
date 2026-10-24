"""memcap.py - run a command under a resident-memory cap on any Linux.

The command starts in its own process group, and a watchdog thread sums the resident memory (/proc/PID/statm) of every process in that group every 0.2 s, killing the whole group once the sum passes the cap.
It needs no systemd, cgroup delegation or root, and it caps memory actually used rather than address space (RLIMIT_AS would kill tools that reserve far more than they touch, SSHash among them).
A run killed this way has `.capped = True` and returncode -9.

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

PAGE = os.sysconf("SC_PAGE_SIZE")
POLL_S = 0.2


def parse_size(s):
    """Bytes of a size such as '26G', '512M' or '1048576', with '' or None giving 0 for no cap."""
    if not s:
        return 0
    s = str(s).strip().upper()
    mult = {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}
    return int(float(s[:-1]) * mult[s[-1]]) if s[-1] in mult else int(s)


def group_rss(pgid):
    """Resident bytes of every process whose process group is pgid."""
    tot = 0
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                st = f.read()
            # Fields after the parenthesised command name are state, ppid, pgrp and so on.
            if int(st[st.rindex(")") + 2 :].split()[2]) != pgid:
                continue
            with open(f"/proc/{d}/statm") as f:
                tot += int(f.read().split()[1]) * PAGE
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
                    pass
                return
            time.sleep(POLL_S)


def popen(cmd, cap=0, **kw):
    """subprocess.Popen in its own process group, watched against `cap` resident bytes, 0 meaning no cap.
    After the process exits, `.capped` tells whether the watchdog killed it."""
    kw.setdefault("start_new_session", True)
    p = subprocess.Popen(cmd, **kw)
    p._watch = _Watch(p, cap) if cap else None
    if p._watch:
        p._watch.start()
    return p


def capped(p):
    """True when the watchdog of process p killed it at the cap."""
    return bool(getattr(p, "_watch", None) and p._watch.capped)


def run(cmd, cap=0, timeout=None, check=False, input=None, **kw):
    """subprocess.run under the cap, returning a CompletedProcess that carries `.capped`.

    Ctrl-C or a timeout kills the whole group, not only its leader."""
    if kw.pop("capture_output", False):  # as subprocess.run handles it
        kw["stdout"] = kw["stderr"] = subprocess.PIPE
    if input is not None and "stdin" not in kw:
        kw["stdin"] = subprocess.PIPE
    p = popen(cmd, cap, **kw)
    try:
        out, err = p.communicate(input=input, timeout=timeout)
    except BaseException:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        p.wait()
        raise
    if p._watch:
        p._watch.join()
    r = subprocess.CompletedProcess(p.args, p.returncode, out, err)
    r.capped = capped(p)
    if check and r.returncode != 0:
        raise subprocess.CalledProcessError(r.returncode, p.args, out, err)
    return r


def warm(paths):
    """Read input files into the page cache before the tools run, so no tool is timed reading from disk while the others read from memory.

    Inputs larger than half the available memory cannot stay cached and are left alone, so every tool streams them from disk alike.
    Returns the bytes read."""
    paths = [p for p in dict.fromkeys(paths) if p and os.path.isfile(p)]
    tot = sum(os.path.getsize(p) for p in paths)
    with open("/proc/meminfo") as fh:
        avail = (
            next(int(x.split()[1]) for x in fh if x.startswith("MemAvailable:")) * 1024
        )
    if not paths or tot > avail // 2:
        return 0
    for p in paths:
        with open(p, "rb") as f:
            while f.read(1 << 24):
                pass
    return tot


if __name__ == "__main__":
    # `python3 bench/memcap.py --warm FILE...` runs warm() as a command for hyperfine --prepare, so every timed run begins with identical cached inputs.
    import sys

    if len(sys.argv) >= 2 and sys.argv[1] == "--warm":
        warm(sys.argv[2:])
    else:
        sys.exit("usage: memcap.py --warm FILE...")
