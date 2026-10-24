#!/usr/bin/env python3
"""thread_check.py - checks that optik keeps no more threads busy than -t.

-t is a budget for everything a command keeps running at once, its worker pools, its input readers, each decoding a .gz inside its own thread, and the OpenMP regions of the C++ core.
The script runs all subcommands on small generated inputs at -t 1, 2 and 4 and without -t, and watches the whole process tree through /proc, with neither root nor perf.
Every 5 ms it reads each task's time on a CPU (/proc/PID/task/TID/schedstat) and counts the tasks in state R, and at the end it takes the CPU time the tree used from wait4.

Gates per run
  busy      over every window of at least BUSY_WINDOW_S, the tree's time on a CPU per wall second is at most t + BUSY_SLACK.
            A whole extra thread at work for more than a few milliseconds fails it, while a thread that is only runnable, waiting for a core or yielding while a pool winds down, uses no CPU and does not count.
            The report still shows the peak number of tasks in state R and the share of samples above t, for reference.
  cpu       (user + system) / wall <= t + 0.1
  default   without -t the budget is half the usable CPUs, hw.cpus() with affinity and cgroup quota applied, so the same two gates apply at t = cpus / 2, and dbruijn reports "with <t> thread(s)"

The inputs, generated with fixed seeds in --work, are a 2 Mbp random genome, 30x of 150 bp reads from it as FASTA and as gzip (two files, for parallel readers), and paired FASTQ.gz for the mapper.

Usage: python3 bench/thread_check.py [--optik target/release/optik] [--work DIR] The exit status is 0 when every gate passes.
"""

import argparse
import gzip
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

import hw

R = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# The CPUs optik's available_parallelism() sees, with affinity and cgroup quota applied.
CPUS = hw.cpus()
DEFAULT_T = max(1, CPUS // 2)
# Sampling period, the shortest span the busy gate judges and the CPU it allows above t.
# A whole extra thread working for a quarter of a window, 12.5 ms, exceeds the slack.
SAMPLE_S = 0.005
BUSY_WINDOW_S = 0.05
BUSY_SLACK = 0.25


def tree(pid):
    """pid and every descendant, found through the parent links in /proc/*/stat."""
    kids = {}
    for e in os.listdir("/proc"):
        if not e.isdigit():
            continue
        try:
            with open(f"/proc/{e}/stat") as f:
                s = f.read()
        except OSError:
            continue
        ppid = int(s[s.rindex(")") + 2 :].split()[1])
        kids.setdefault(ppid, []).append(int(e))
    out, todo = [], [pid]
    while todo:
        p = todo.pop()
        out.append(p)
        todo += kids.get(p, [])
    return out


def running(pids):
    """Tasks in state R across the given processes, and each task's time on a CPU in ns, keyed by thread id, from the first field of its schedstat."""
    n, on_cpu = 0, {}
    for p in pids:
        try:
            tids = os.listdir(f"/proc/{p}/task")
        except OSError:
            continue
        for t in tids:
            try:
                with open(f"/proc/{p}/task/{t}/stat") as f:
                    s = f.read()
                with open(f"/proc/{p}/task/{t}/schedstat") as f:
                    on_cpu[int(t)] = int(f.read().split()[0])
            except (OSError, ValueError, IndexError):
                continue
            if s[s.rindex(")") + 2] == "R":
                n += 1
    return n, on_cpu


def peak_busy(samples, window):
    """Highest CPU time per wall second of the tree over any span of at least
    @window seconds, from (time, on-CPU ns per thread) samples.

    A thread's time counts from the sample it first appears in, and the time a thread used between its last sample and its exit is lost, which can only lower the result."""
    times, cum, prev, total = [], [], {}, 0
    for when, on_cpu in samples:
        for tid, ns in on_cpu.items():
            total += max(0, ns - prev.get(tid, ns))
        prev = on_cpu
        times.append(when)
        cum.append(total)
    best, j = 0.0, 0
    for i in range(len(times)):
        while j < len(times) and times[j] - times[i] < window:
            j += 1
        if j == len(times):
            break
        best = max(best, (cum[j] - cum[i]) / 1e9 / (times[j] - times[i]))
    return best


def cpu_of(cmd, errf):
    """Run cmd and return (exit code, wall s, user + system s of the whole tree, samples of (time, tasks in state R, on-CPU ns per thread), stderr).
    wait4's rusage includes the children the command reaped."""
    t0 = time.monotonic()
    pid = os.fork()
    if pid == 0:
        try:
            fd = os.open(os.devnull, os.O_WRONLY)
            os.dup2(fd, 1)
            err_fd = os.open(errf, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(err_fd, 2)
            os.execvp(cmd[0], cmd)
        finally:
            os._exit(127)
    samples = []
    while True:
        done, status, ru = os.wait4(pid, os.WNOHANG)
        if done:
            break
        n_running, on_cpu = running(tree(pid))
        samples.append((time.monotonic(), n_running, on_cpu))
        time.sleep(SAMPLE_S)
    wall = time.monotonic() - t0
    code = os.waitstatus_to_exitcode(status)
    return code, wall, ru.ru_utime + ru.ru_stime, samples, Path(errf).read_text()


def make_inputs(d):
    """Generate the genome, reads and paired FASTQ.gz in d once, returning their paths."""
    os.makedirs(d, exist_ok=True)
    g = os.path.join(d, "genome.fa")
    if not os.path.exists(g):
        r = random.Random(5)
        seq = "".join(r.choice("ACGT") for _ in range(2_000_000))
        with open(g, "w") as f:
            f.write(
                ">g\n"
                + "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80))
                + "\n"
            )
    else:
        with open(g) as fh:
            seq = "".join(x.strip() for x in fh if not x.startswith(">"))
    comp = str.maketrans("ACGT", "TGCA")

    def reads(path, n, seed, fq=False, opener=open):
        """Write n random 150 bp reads of the genome as FASTA or FASTQ."""
        r = random.Random(seed)
        with opener(path, "wt") as f:
            for i in range(n):
                p = r.randrange(len(seq) - 150)
                s = list(seq[p : p + 150])
                for j in range(150):
                    if r.random() < 0.005:
                        s[j] = r.choice("ACGT")
                s = "".join(s)
                if r.random() < 0.5:
                    s = s.translate(comp)[::-1]
                f.write(f"@r{i}\n{s}\n+\n{'I' * 150}\n" if fq else f">r{i}\n{s}\n")

    out = {"genome": g}
    for name, n, seed, fq, op in (
        ("reads1.fa.gz", 200_000, 1, False, gzip.open),
        ("reads2.fa.gz", 200_000, 2, False, gzip.open),
        ("reads.fa", 200_000, 3, False, open),
    ):
        p = os.path.join(d, name)
        if not os.path.exists(p):
            reads(p, n, seed, fq, op)
        out[name] = p
    # Paired mates come from the same fragments, with R2 the reverse complement of the far end.
    r1, r2 = os.path.join(d, "pe_1.fq.gz"), os.path.join(d, "pe_2.fq.gz")
    if not os.path.exists(r2):
        r = random.Random(7)
        with gzip.open(r1, "wt") as a, gzip.open(r2, "wt") as b:
            for i in range(100_000):
                p = r.randrange(len(seq) - 500)
                m1 = seq[p : p + 150]
                m2 = seq[p + 350 : p + 500].translate(comp)[::-1]
                a.write(f"@p{i}/1\n{m1}\n+\n{'I' * 150}\n")
                b.write(f"@p{i}/2\n{m2}\n+\n{'I' * 150}\n")
    out["r1"], out["r2"] = r1, r2
    return out


def main():
    """Run every subcommand at each thread budget and apply the gates."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--optik", default=os.path.join(R, "target/release/optik"))
    ap.add_argument("--work", default=os.path.join(R, "temp/thread_check"))
    ap.add_argument(
        "--threads", default="1,2,4", help="explicit -t values (plus the default)"
    )
    a = ap.parse_args()
    optik = a.optik
    inputs = make_inputs(a.work)

    def in_work(f):
        """Path of f in the work directory."""
        return os.path.join(a.work, f)

    # The indexes the query rows read are built once and not judged here.
    subprocess.run(
        [
            optik,
            "build-spss",
            "-i",
            inputs["genome"],
            "-k",
            "31",
            "-t",
            "2",
            "-o",
            in_work("g.spss"),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            optik,
            "build-fm",
            "-i",
            inputs["genome"],
            "-k",
            "31",
            "-t",
            "2",
            "-o",
            in_work("g.fm"),
        ],
        check=True,
        capture_output=True,
    )
    jobs = {
        "dbruijn genome": [
            optik,
            "dbruijn",
            "-i",
            inputs["genome"],
            "-k",
            "31",
            "-c",
            "1",
            "-o",
            in_work("o.gfa"),
        ],
        "dbruijn reads .gz x2 -c 2": [
            optik,
            "dbruijn",
            "-i",
            inputs["reads1.fa.gz"],
            "-i",
            inputs["reads2.fa.gz"],
            "-k",
            "31",
            "-c",
            "2",
            "-o",
            in_work("o.fa"),
        ],
        "dbruijn reads --graph table": [
            optik,
            "dbruijn",
            "-i",
            inputs["reads.fa"],
            "-k",
            "31",
            "-c",
            "2",
            "--graph",
            "table",
            "-o",
            in_work("o2.fa"),
        ],
        "dbruijn reads -k 15 (table)": [
            optik,
            "dbruijn",
            "-i",
            inputs["reads.fa"],
            "-k",
            "15",
            "-c",
            "2",
            "-o",
            in_work("o3.fa"),
        ],
        "build-spss genome": [
            optik,
            "build-spss",
            "-i",
            inputs["genome"],
            "-k",
            "31",
            "-o",
            in_work("o.spss"),
        ],
        "build-spss --compact -c 2": [
            optik,
            "build-spss",
            "--compact",
            "-i",
            inputs["reads.fa"],
            "-k",
            "31",
            "-c",
            "2",
            "-o",
            in_work("o2.spss"),
        ],
        "build-fm": [
            optik,
            "build-fm",
            "-i",
            inputs["genome"],
            "-k",
            "31",
            "-o",
            in_work("o.fm"),
        ],
        "query spss": [
            optik,
            "query",
            "-i",
            in_work("g.spss"),
            "-q",
            inputs["reads1.fa.gz"],
        ],
        "query fm": [
            optik,
            "query",
            "-i",
            in_work("g.fm"),
            "-q",
            inputs["reads1.fa.gz"],
        ],
        "map single-end .gz": [
            optik,
            "map",
            "-r",
            inputs["genome"],
            "-q",
            inputs["reads1.fa.gz"],
            "-o",
            in_work("o.sam"),
        ],
        "map paired .gz": [
            optik,
            "map",
            "-r",
            inputs["genome"],
            "-q",
            inputs["r1"],
            "-q",
            inputs["r2"],
            "-o",
            in_work("o2.sam"),
        ],
    }
    ts = [int(x) for x in a.threads.split(",")]
    print(f"{CPUS} usable CPUs; default budget {DEFAULT_T}")
    print(
        f"{'command':30s} {'-t':>4s} {'budget':>6s} {'wall s':>7s} {'cpu/wall':>8s} {'busy':>5s} "
        f"{'max R':>5s} {'R>t %':>6s}  result"
    )
    bad = 0
    for name, cmd in jobs.items():
        for t in ts + [None]:
            full = cmd + (["-t", str(t)] if t is not None else [])
            budget = t if t is not None else DEFAULT_T
            code, wall, cpu, samples, err = cpu_of(full, in_work("stderr.txt"))
            n_r = [n for _, n, _ in samples]
            over = sum(1 for n in n_r if n > budget)
            frac = over / len(samples) if samples else 0.0
            util = cpu / wall if wall > 0 else 0.0
            busy = peak_busy([(w, c) for w, _, c in samples], BUSY_WINDOW_S)
            why = []
            if code != 0:
                why.append(
                    f"exit {code}: {err.strip().splitlines()[-1] if err.strip() else ''}"
                )
            if busy > budget + BUSY_SLACK:
                why.append(
                    f"{busy:.2f} CPUs busy over {BUSY_WINDOW_S} s > {budget} + {BUSY_SLACK}"
                )
            if util > budget + 0.1:
                why.append(f"cpu/wall {util:.2f} > {budget} + 0.1")
            if t is None and name.startswith("dbruijn"):
                m = re.search(r"with (\d+) thread\(s\)", err)
                if not m or int(m.group(1)) != DEFAULT_T:
                    why.append(
                        f"reports {m.group(1) if m else '?'} threads, default should be {DEFAULT_T}"
                    )
            bad += bool(why)
            print(
                f"{name:30s} {('-' if t is None else str(t)):>4s} {budget:6d} {wall:7.2f} {util:8.2f} "
                f"{busy:5.2f} {max(n_r, default=0):5d} {100 * frac:6.2f}  {'PASS' if not why else 'FAIL ' + '; '.join(why)}",
                flush=True,
            )
    print("THREAD CHECK: " + ("PASS" if not bad else f"FAIL ({bad} runs)"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
