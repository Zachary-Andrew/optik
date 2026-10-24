"""hw.py - the CPUs and memory this process may use, so every bench scales with its machine.

optik resolves its default thread count from Rust's available_parallelism(), which honors the CPU affinity mask and a cgroup v2 CPU quota, and its default counting budget from MemAvailable and the cgroup v2 memory limits on its path (include/spill.hpp).
The benches read the same limits here, so a bench under taskset, in a container or in a systemd-run scope sizes its threads and memory caps to what the tools themselves will see.

  cpus()               CPUs usable by this process
  mem_bytes()          MemTotal, bounded by every cgroup v2 memory.max on the path
  default_cap_bytes()  the default resident-memory cap of a bench run
  thread_points(lo)    1, 2, 4, ... up to cpus(), plus cpus() itself
"""

import math
import os

CGROUP_ROOT = "/sys/fs/cgroup"
# Room left outside a bench's memory cap for the bench itself, the page cache of its inputs and the rest of the system.
CAP_HEADROOM_BYTES = 2 << 30


def _cgroup_dirs():
    """The cgroup v2 directories from this process's own cgroup up to the root.

    Empty under cgroup v1 or without /proc/self/cgroup."""
    try:
        with open("/proc/self/cgroup") as f:
            rel = next((x[3:].strip() for x in f if x.startswith("0::")), None)
    except OSError:
        return []
    if rel is None:
        return []
    dirs, path = [], CGROUP_ROOT + (rel if rel != "/" else "")
    while path.startswith(CGROUP_ROOT):
        dirs.append(path)
        if path == CGROUP_ROOT:
            break
        path = os.path.dirname(path)
    return dirs


def _read(path):
    """Stripped text of a small file, or None when it cannot be read."""
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return None


def cpus():
    """CPUs this process may use.

    The affinity mask, bounded by the tightest cgroup v2 cpu.max quota (quota over period, rounded up) on the process's path.
    This is the count Rust's available_parallelism() gives optik."""
    n = len(os.sched_getaffinity(0))
    for d in _cgroup_dirs():
        v = _read(os.path.join(d, "cpu.max"))
        if not v:
            continue
        quota, _, period = v.partition(" ")
        if quota != "max" and period:
            n = min(n, max(1, math.ceil(int(quota) / int(period))))
    return max(1, n)


def mem_bytes():
    """MemTotal, bounded by the tightest cgroup v2 memory.max on this process's path."""
    mem = None
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemTotal:"):
                mem = int(line.split()[1]) * 1024
                break
    for d in _cgroup_dirs():
        v = _read(os.path.join(d, "memory.max"))
        if v and v != "max":
            mem = min(mem, int(v))
    return mem


def default_cap_bytes():
    """The default resident-memory cap of a bench run, mem_bytes() less CAP_HEADROOM_BYTES.

    Never below half of mem_bytes(), so a small machine still gets a usable cap."""
    mem = mem_bytes()
    return max(mem // 2, mem - CAP_HEADROOM_BYTES)


def default_cap_str():
    """default_cap_bytes() as a size string memcap.parse_size() reads, in whole MiB."""
    return f"{default_cap_bytes() >> 20}M"


def thread_points(lo=1, n=None):
    """Thread counts of a scaling sweep, i.e. powers of two from @lo up to @n, plus @n itself when it is not a power of two.
    @n defaults to cpus().

    With 12 CPUs this gives 1, 2, 4, 8, 12, and with 32 it gives 1 to 32."""
    n = cpus() if n is None else n
    lo = max(1, min(lo, n))
    pts, t = [], 1
    while t <= n:
        if t >= lo:
            pts.append(t)
        t *= 2
    if not pts or pts[-1] != n:
        pts.append(n)
    if pts[0] != lo:
        pts.insert(0, lo)
    return pts


if __name__ == "__main__":
    print(f"cpus {cpus()}")
    print(f"mem_gib {mem_bytes() / 2**30:.1f}")
    print(f"default_cap {default_cap_str()}")
    print(f"thread_points {thread_points()}")
