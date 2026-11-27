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
from units import NS_PER_SECOND, PERCENT

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# The CPU count optik's available_parallelism() reports, with affinity and cgroup quota applied.
CPUS = hw.cpus()
DEFAULT_T = max(1, CPUS // 2)
# Sampling period, the shortest span the busy gate judges and the CPU it allows above t.
# A whole extra thread working for a quarter of a window, 12.5 ms, exceeds the slack.
SAMPLE_S = 0.005
BUSY_WINDOW_S = 0.05
BUSY_SLACK = 0.25
# whole-run CPU/wall may exceed the thread budget by this many CPUs
UTIL_SLACK = 0.1
# a child whose exec failed exits with the shell's "command not found" status
EXEC_FAILED_EXIT = 127
# permissions of a captured stderr file
ERR_FILE_MODE = 0o644

# The inputs are synthetic, made of a random genome, single-end read sets and a paired library.
GENOME_SEED = 5
GENOME_LEN = 2_000_000
FASTA_LINE_WIDTH = 80
READ_LEN = 150
READ_SUB_RATE = 0.005  # per-base substitution probability
REVCOMP_PROB = 0.5  # probability that a read is taken from the reverse strand
FAKE_QUAL = "I"  # FASTQ quality character of every base
# (file name, reads, seed, FASTQ?, opener) of each single-end set
READS_PER_SET = 200_000
READ_SETS = (
    ("reads1.fa.gz", READS_PER_SET, 1, False, gzip.open),
    ("reads2.fa.gz", READS_PER_SET, 2, False, gzip.open),
    ("reads.fa", READS_PER_SET, 3, False, open),
)
PAIRS_SEED = 7
NUM_PAIRS = 100_000
FRAGMENT_LEN = 500  # R1 is the first READ_LEN bases of the fragment, R2 the reverse complement of the last READ_LEN


# /proc/<pid>/stat puts the command name in parentheses, then state and ppid (proc(5)).
STAT_COMM_END = ") "
STAT_PPID = 1
# The child's stderr descriptor, which execvp passes on.
STDERR_FILENO = 2


def tree(pid):
    """pid and every descendant, found through the parent links in /proc/*/stat."""
    kids = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat") as handle:
                stat_line = handle.read()
        except OSError:
            continue  # the process exited while /proc was being scanned
        ppid = int(stat_line[stat_line.rindex(STAT_COMM_END) + len(STAT_COMM_END) :].split()[STAT_PPID])
        kids.setdefault(ppid, []).append(int(entry))
    out, todo = [], [pid]
    while todo:
        proc_pid = todo.pop()
        out.append(proc_pid)
        todo += kids.get(proc_pid, [])
    return out


def running(pids):
    """Tasks in state R across the given processes, and each task's time on a CPU in ns, keyed by thread id, from the first field of its schedstat."""
    num_running, on_cpu = 0, {}
    for proc_pid in pids:
        try:
            tids = os.listdir(f"/proc/{proc_pid}/task")
        except OSError:
            continue  # the process exited while /proc was being scanned
        for tid in tids:
            try:
                with open(f"/proc/{proc_pid}/task/{tid}/stat") as handle:
                    stat_line = handle.read()
                with open(f"/proc/{proc_pid}/task/{tid}/schedstat") as handle:
                    on_cpu[int(tid)] = int(handle.read().split()[0])
            except (OSError, ValueError, IndexError):
                continue
            if stat_line[stat_line.rindex(STAT_COMM_END) + len(STAT_COMM_END)] == "R":
                num_running += 1
    return num_running, on_cpu


def peak_busy(samples, window):
    """Highest CPU time per wall second of the tree over any span of at least
    @window seconds, from (time, on-CPU ns per thread) samples.

    A thread's time counts from the sample it first appears in, and the time a thread used between its last sample and its exit is lost, which can only lower the result."""
    times, cum, prev, total = [], [], {}, 0
    for when, on_cpu in samples:
        for tid, cpu_ns in on_cpu.items():
            total += max(0, cpu_ns - prev.get(tid, cpu_ns))
        prev = on_cpu
        times.append(when)
        cum.append(total)
    best, j = 0.0, 0
    for i in range(len(times)):
        while j < len(times) and times[j] - times[i] < window:
            j += 1
        if j == len(times):
            break
        best = max(best, (cum[j] - cum[i]) / NS_PER_SECOND / (times[j] - times[i]))
    return best


def cpu_of(cmd, errf):
    """Run cmd and return (exit code, wall s, user + system s of the whole tree, samples of (time, tasks in state R, on-CPU ns per thread), stderr).
    wait4's rusage includes the children the command reaped."""
    start = time.monotonic()
    pid = os.fork()
    if pid == 0:
        try:
            null_fd = os.open(os.devnull, os.O_WRONLY)
            os.dup2(null_fd, 1)
            err_fd = os.open(errf, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, ERR_FILE_MODE)
            os.dup2(err_fd, STDERR_FILENO)
            os.execvp(cmd[0], cmd)
        finally:
            os._exit(EXEC_FAILED_EXIT)
    samples = []
    while True:
        done, status, rusage = os.wait4(pid, os.WNOHANG)
        if done:
            break
        n_running, on_cpu = running(tree(pid))
        samples.append((time.monotonic(), n_running, on_cpu))
        time.sleep(SAMPLE_S)
    wall = time.monotonic() - start
    code = os.waitstatus_to_exitcode(status)
    return code, wall, rusage.ru_utime + rusage.ru_stime, samples, Path(errf).read_text()


def make_inputs(work_dir):
    """Generate the genome, reads and paired FASTQ.gz in d once, returning their paths."""
    os.makedirs(work_dir, exist_ok=True)
    genome_path = os.path.join(work_dir, "genome.fa")
    if not os.path.exists(genome_path):
        rng = random.Random(GENOME_SEED)
        seq = "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))
        with open(genome_path, "w") as handle:
            handle.write(
                ">g\n"
                + "\n".join(seq[i : i + FASTA_LINE_WIDTH] for i in range(0, len(seq), FASTA_LINE_WIDTH))
                + "\n"
            )
    else:
        with open(genome_path) as handle:
            seq = "".join(line.strip() for line in handle if not line.startswith(">"))
    comp = str.maketrans("ACGT", "TGCA")

    def reads(path, num_reads, seed, fastq=False, opener=open):
        """Write n random 150 bp reads of the genome as FASTA or FASTQ."""
        rng = random.Random(seed)
        with opener(path, "wt") as handle:
            for i in range(num_reads):
                pos = rng.randrange(len(seq) - READ_LEN)
                read_seq = list(seq[pos : pos + READ_LEN])
                for j in range(READ_LEN):
                    if rng.random() < READ_SUB_RATE:
                        read_seq[j] = rng.choice("ACGT")
                read_seq = "".join(read_seq)
                if rng.random() < REVCOMP_PROB:
                    read_seq = read_seq.translate(comp)[::-1]
                handle.write(f"@r{i}\n{read_seq}\n+\n{FAKE_QUAL * READ_LEN}\n" if fastq else f">r{i}\n{read_seq}\n")

    out = {"genome": genome_path}
    for name, num_reads, seed, fastq, opener in READ_SETS:
        path = os.path.join(work_dir, name)
        if not os.path.exists(path):
            reads(path, num_reads, seed, fastq, opener)
        out[name] = path
    # Paired mates come from the same fragments, with R2 the reverse complement of the far end.
    mate1_path, mate2_path = os.path.join(work_dir, "pe_1.fq.gz"), os.path.join(work_dir, "pe_2.fq.gz")
    if not os.path.exists(mate2_path):
        rng = random.Random(PAIRS_SEED)
        with gzip.open(mate1_path, "wt") as mate1_handle, gzip.open(mate2_path, "wt") as mate2_handle:
            for i in range(NUM_PAIRS):
                path = rng.randrange(len(seq) - FRAGMENT_LEN)
                mate1 = seq[path : path + READ_LEN]
                mate2 = seq[path + FRAGMENT_LEN - READ_LEN : path + FRAGMENT_LEN].translate(comp)[::-1]
                mate1_handle.write(f"@p{i}/1\n{mate1}\n+\n{FAKE_QUAL * READ_LEN}\n")
                mate2_handle.write(f"@p{i}/2\n{mate2}\n+\n{FAKE_QUAL * READ_LEN}\n")
    out["r1"], out["r2"] = mate1_path, mate2_path
    return out


def main():
    """Run every subcommand at each thread budget and apply the gates."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    parser.add_argument("--work", default=os.path.join(ROOT, "temp/thread_check"))
    parser.add_argument(
        "--threads", default="1,2,4", help="explicit -t values (plus the default)"
    )
    args = parser.parse_args()
    optik = args.optik
    inputs = make_inputs(args.work)

    def in_work(name):
        """Path of f in the work directory."""
        return os.path.join(args.work, name)

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
    # Each job is one command line, run once per -t value and once without -t to check the default budget.
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
    thread_counts = [int(value) for value in args.threads.split(",")]
    print(f"{CPUS} usable CPUs; default budget {DEFAULT_T}")
    print(
        f"{'command':30s} {'-t':>4s} {'budget':>6s} {'wall s':>7s} {'cpu/wall':>8s} {'busy':>5s} "
        f"{'max R':>5s} {'R>t %':>6s}  result"
    )
    bad = 0
    for name, cmd in jobs.items():
        for threads in thread_counts + [None]:
            full = cmd + (["-t", str(threads)] if threads is not None else [])
            budget = threads if threads is not None else DEFAULT_T
            code, wall, cpu, samples, err = cpu_of(full, in_work("stderr.txt"))
            n_r = [num_running for _, num_running, _ in samples]
            over = sum(1 for num_running in n_r if num_running > budget)
            frac = over / len(samples) if samples else 0.0
            util = cpu / wall if wall > 0 else 0.0
            # busy is the most CPU time per wall second over any BUSY_WINDOW_S span, which catches bursts that the whole-run average hides.
            busy = peak_busy([(when, on_cpu) for when, _, on_cpu in samples], BUSY_WINDOW_S)
            why = []
            if code != 0:
                why.append(
                    f"exit {code}: {err.strip().splitlines()[-1] if err.strip() else ''}"
                )
            if busy > budget + BUSY_SLACK:
                why.append(
                    f"{busy:.2f} CPUs busy over {BUSY_WINDOW_S} s > {budget} + {BUSY_SLACK}"
                )
            if util > budget + UTIL_SLACK:
                why.append(f"cpu/wall {util:.2f} > {budget} + {UTIL_SLACK}")
            # Without -t, dbruijn must also report that it used the default budget.
            if threads is None and name.startswith("dbruijn"):
                match = re.search(r"with (\d+) thread\(s\)", err)
                if not match or int(match.group(1)) != DEFAULT_T:
                    why.append(
                        f"reports {match.group(1) if match else '?'} threads, default should be {DEFAULT_T}"
                    )
            bad += bool(why)
            print(
                f"{name:30s} {('-' if threads is None else str(threads)):>4s} {budget:6d} {wall:7.2f} {util:8.2f} "
                f"{busy:5.2f} {max(n_r, default=0):5d} {PERCENT * frac:6.2f}  {'PASS' if not why else 'FAIL ' + '; '.join(why)}",
                flush=True,
            )
    print("THREAD CHECK: " + ("PASS" if not bad else f"FAIL ({bad} runs)"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
