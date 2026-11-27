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

from units import BYTES_PER_GIB, BYTES_PER_KIB, BYTES_PER_MIB

CGROUP_ROOT = "/sys/fs/cgroup"
# the cgroup v2 line of /proc/self/cgroup starts with this prefix
CGROUP_V2_PREFIX = "0::"
# Room left outside a bench's memory cap for the bench itself, the page cache of its inputs and the rest of the system.
CAP_HEADROOM_BYTES = 2 << 30


# The default cap never drops below half of the memory.
MIN_CAP_SHARE_DEN = 2
# Thread sweeps double the count at each step.
SWEEP_GROWTH = 2


def _cgroup_dirs():
    """The cgroup v2 directories from this process's own cgroup up to the root.

    Empty under cgroup v1 or without /proc/self/cgroup."""
    try:
        with open("/proc/self/cgroup") as handle:
            rel = next((line[len(CGROUP_V2_PREFIX):].strip() for line in handle if line.startswith(CGROUP_V2_PREFIX)), None)
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
        with open(path) as handle:
            return handle.read().strip()
    except OSError:
        return None


def cpus():
    """CPUs this process may use.

    The affinity mask, bounded by the tightest cgroup v2 cpu.max quota (quota over period, rounded up) on the process's path.
    This is the count Rust's available_parallelism() gives optik."""
    count = len(os.sched_getaffinity(0))
    for cgroup_dir in _cgroup_dirs():
        cpu_max = _read(os.path.join(cgroup_dir, "cpu.max"))
        if not cpu_max:
            continue
        quota, _, period = cpu_max.partition(" ")
        if quota != "max" and period:
            count = min(count, max(1, math.ceil(int(quota) / int(period))))
    return max(1, count)


def mem_bytes():
    """MemTotal, bounded by the tightest cgroup v2 memory.max on this process's path."""
    mem = None
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemTotal:"):
                mem = int(line.split()[1]) * BYTES_PER_KIB
                break
    for cgroup_dir in _cgroup_dirs():
        mem_max = _read(os.path.join(cgroup_dir, "memory.max"))
        if mem_max and mem_max != "max":
            mem = min(mem, int(mem_max))
    return mem


def default_cap_bytes():
    """The default resident-memory cap of a bench run, mem_bytes() less CAP_HEADROOM_BYTES.

    Never below half of mem_bytes(), so a small machine still gets a usable cap."""
    mem = mem_bytes()
    return max(mem // MIN_CAP_SHARE_DEN, mem - CAP_HEADROOM_BYTES)


def default_cap_str():
    """default_cap_bytes() as a size string memcap.parse_size() reads, in whole MiB."""
    return f"{default_cap_bytes() // BYTES_PER_MIB}M"


def thread_points(min_threads=1, max_threads=None):
    """Thread counts of a scaling sweep, i.e. powers of two from @lo up to @n, plus @n itself when it is not a power of two.
    @n defaults to cpus().

    With 12 CPUs this gives 1, 2, 4, 8, 12, and with 32 it gives 1 to 32."""
    max_threads = cpus() if max_threads is None else max_threads
    min_threads = max(1, min(min_threads, max_threads))
    pts, threads = [], 1
    while threads <= max_threads:
        if threads >= min_threads:
            pts.append(threads)
        threads *= SWEEP_GROWTH
    if not pts or pts[-1] != max_threads:
        pts.append(max_threads)
    if pts[0] != min_threads:
        pts.insert(0, min_threads)
    return pts


if __name__ == "__main__":
    print(f"cpus {cpus()}")
    print(f"mem_gib {mem_bytes() / BYTES_PER_GIB:.1f}")
    print(f"default_cap {default_cap_str()}")
    print(f"thread_points {thread_points()}")
