#!/usr/bin/env python3
"""compare_tools_small_sets.py - optik against every published tool it implements, on small inputs.

The inputs are SSHash's E. coli, 100-Salmonella and Penicillium unitig sets for the k-mer indexes, the Salmonella LT2 genome and 30x simulated reads for the graphs, simulated short, ONT-like, HiFi-like, chimeric and random reads for mapping, random sequence pairs for alignment, and a repeat-rich synthetic genome for MAPQ calibration.
bench/compare_tools_published_sets.py runs the same comparisons on the datasets of the tools' own papers.

Two modes, run alone or together

  parity   identical inputs and matched parameters for every tool in a section, i.e. k = 31, minimizer m = 15 where a tool takes one, and the same alignment penalties and thread counts, with whole-process wall time and single-threaded queries.
           This table compares algorithms.
  optimal  every tool at its own documented best settings, meaning its presets, its recommended minimizer length, variant or query mode, and all physical cores where it parallelises.
           This is the table a user would see.

Sections
  kmer   exact k-mer membership and dictionary lookups with optik build-spss, optik build-fm, SSHash, SBWT and FMSI with KmerCamel on SSHash's unitig sets
  dbg    compacted de Bruijn graphs from optik dbruijn, Cuttlefish and Bifrost on a genome and on 30x simulated reads
  map    read mapping with optik map, minimap2 and, for long reads, wfmash on short, ONT-like, HiFi-like, chimeric and random reads with truth
  align  exact gap-affine pairwise alignment, optik WFA and BiWFA against WFA2-lib's WFA and BiWFA
  mapq   MAPQ calibration of optik against minimap2 by bench/mapq_calibration.py

Each row carries a correctness gate.
A row whose output is wrong is marked FAIL in the tables and is never silently timed.
The gates are
  kmer   num_kmers and num_positive equal an independent brute-force truth
  dbg    the unitig k-mer set equals KMC's solid k-mers at the same cutoff, and no k-mer appears in two unitigs
  map    every SAM record consumes SEQ exactly, NM recomputed from the reference matches NM:i, each read has one primary, and accuracy is scored against the truth
  align  all exact methods report the same optimal scores, and each CIGAR re-scores to its own score

Timing goes through bench/measure.py, which runs hyperfine on the whole process, drops the --warmup runs and computes statistics over the --runs measured runs, with each run's inputs cached beforehand.
Peak RSS and block I/O of every measured run come from wait4 in rusage_exec.
The correctness gates read the measured runs' own output.
Per-k-mer answers are the exception, since one untimed run writes them and the timed runs send them to /dev/null, as counting tools write none.
Every command line, tool revision, compiler and machine detail is saved in the run directory.

The tool paths and sequence helpers live in small_common.py, and the external correctness oracles built on KMC, Picard and samtools in external_oracles.py.

Usage:
  python3 bench/compare_tools_small_sets.py [--mode parity|optimal|both]
      [--sections kmer,dbg,map,align,mapq] [--runs 3] [--warmup 1]
      [--threads <usable CPUs>] [--out bench/results/compare_tools_small_sets] [--quick]
"""

import argparse
import contextlib
import datetime
import gzip
import json
import os
import platform
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import hw
import measure
import memcap
import report as results_file
from external_oracles import (
    external_sam_checks,
    judge_unitigs,
    kmc_count,
    kmc_db,
    kmc_hits,
    mem_budget_gb,
    same_unitigs,
)
from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
)
from small_common import (
    BIN,
    DATA,
    T,
    kmer_codes,
    mutate,
    rc,
    read_fasta,
    write_fasta,
)

# Execution and timing


class Bench:
    """Runs and times the bench's commands and collects its result rows."""

    def __init__(self, out, runs, warmup, cap_gb, scratch):
        self.out, self.runs, self.warmup = out, runs, warmup
        # Hyperfine JSON, stderr and rusage files go to one directory per process, so two benches sharing a work directory never read each other's results.
        self.scratch = os.path.join(scratch, str(os.getpid()))
        os.makedirs(self.scratch, exist_ok=True)
        self.cap = cap_gb << 30
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
            ef = stack.enter_context(open(errf, "w"))
            out = (
                stack.enter_context(open(stdout, "w")) if stdout else subprocess.DEVNULL
            )
            r = memcap.run(
                [measure.rusage_exec(), "0", resf] + [str(c) for c in cmd],
                self.cap if cap is None else cap,
                cwd=cwd,
                stdout=out,
                stderr=ef,
            )
        self.last_capped = r.capped
        if not os.path.exists(resf):  # the group was killed before rusage_exec wrote it
            return (
                r.returncode,
                0.0,
                (self.cap if cap is None else cap) / 2**20,
                Path(errf).read_text(),
            )
        code, wall, kb = Path(resf).read_text().split()[:3]
        return int(code), float(wall), int(kb) / 1024.0, Path(errf).read_text()

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
        r = measure.measure(
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
        if stdout and os.path.exists(r.out):
            shutil.copyfile(r.out, stdout)
        self.last_capped = r.capped
        self.last_io = (r.written, r.read)
        return r

    def add(self, **row):
        """Append a result row and print it, failing the row if any measurement failed."""
        # A failed measurement makes a failed row, never a silent PASS.
        for k, v in row.items():
            if isinstance(v, dict) and "error" in v:
                row["gate"] = False
                row["gate_detail"] = (
                    row.get("gate_detail") or ""
                ) + f" timing of {k} failed: {v['error'][-200:]}"
        self.rows.append(row)
        flag = "PASS" if row.get("gate") else "FAIL"
        keys = [
            k for k in row if k not in ("section", "mode", "gate", "gate_detail", "cmd")
        ]
        print(
            f"  [{flag}] " + "  ".join(f"{k}={_fmt(row[k])}" for k in keys), flush=True
        )
        if not row.get("gate"):
            print(f"         gate: {row.get('gate_detail')}", flush=True)


def _q(c):
    """c quoted for a POSIX shell when it holds anything beyond plain path characters."""
    c = str(c)
    return (
        c if re.fullmatch(r"[\w./=:,+-]+", c) else "'" + c.replace("'", "'\\''") + "'"
    )


def _fmt(v):
    """Format a row value for the tables, printing timings as mean±sd seconds."""
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, dict) and "mean" in v:
        # Hyperfine reports no standard deviation for a single run.
        return f"{v['mean']:.4g}s" + (
            f"±{v['sd']:.2g}" if v.get("sd") is not None else ""
        )
    return str(v)


def environment(out):
    """Record the machine, compiler and tool revisions of this run in `out`."""

    def sh(cmd):
        """Output of shell command cmd, or "n/a (reason)" if it fails."""
        try:
            return subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=30, check=False
            ).stdout.strip()
        except (subprocess.SubprocessError, OSError) as e:
            return f"n/a ({e})"

    env = {
        "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "cpu": sh("lscpu | grep 'Model name' | sed 's/.*: *//'"),
        "cores_threads": sh("nproc") + " logical",
        "physical_cores": sh("lscpu -p=core | grep -v '^#' | sort -u | wc -l"),
        "memory": sh("free -g | awk '/Mem/{print $2\" GiB\"}'"),
        "kernel": platform.release(),
        "thp": sh("cat /sys/kernel/mm/transparent_hugepage/enabled"),
        "governor": sh("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "gcc": sh("g++ --version | head -1"),
        "rustc": sh("rustc --version"),
        "optik_git": sh(f"git -C {ROOT} rev-parse --short HEAD")
        + (" (dirty)" if sh(f"git -C {ROOT} status --porcelain -uno") else ""),
        "build_flags": "optik: -O3 -march=native (build.rs) + target-cpu=native; competitors: "
        "bench/build_tools.sh NATIVE=1 (-O3 -march=native; SBWT -O3 -mbmi2 -mpopcnt)",
        "tools": {},
    }
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
        env["tools"][name] = sh(
            f"git -C {T}/{name} rev-parse --short HEAD 2>/dev/null || echo n/a"
        )
    env["tools"]["minimap2_version"] = sh(f"{BIN['minimap2']} --version")
    env["tools"]["cuttlefish_version"] = sh(f"{BIN['cuttlefish']} --version")
    Path(os.path.join(out, "environment.json")).write_text(json.dumps(env, indent=1))
    return env


# Sequence helpers


def check_optik_fresh():
    """Refuse to benchmark an optik binary older than any of its sources, since a stale build would be measured and reported as the current code."""
    exe = BIN["optik"]
    srcs = [
        os.path.join(r, f)
        for top in ("src", "include", "tools/src")
        for r, _, fs in os.walk(os.path.join(ROOT, top))
        for f in fs
        if f.endswith((".rs", ".cpp", ".hpp", ".h"))
    ] + [os.path.join(ROOT, x) for x in ("build.rs", "Cargo.toml")]
    newest = max(srcs, key=os.path.getmtime)
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(newest):
        sys.exit(
            f"error: {exe} is older than {os.path.relpath(newest, ROOT)}: run `cargo build --release` first"
        )


def build_wfab(d):
    """Path of the align bench binary (optik's WFA/BiWFA and WFA2-lib).

    It is rebuilt whenever a source is newer, since a binary kept across code changes would time an old aligner."""
    exe = os.path.join(d, "wfab")
    srcs = [
        os.path.join(ROOT, "include/wfa_affine.hpp"),
        os.path.join(HERE, "cpp/wfa_affine_bench.cpp"),
        os.path.join(T, "WFA2-lib/bindings/cpp/WFAligner.cpp"),
        os.path.join(T, "WFA2-lib/lib/libwfa.a"),
    ]
    if not os.path.exists(exe) or os.path.getmtime(exe) < max(
        os.path.getmtime(x) for x in srcs
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
                os.path.join(T, "WFA2-lib"),
                srcs[1],
                srcs[2],
                srcs[3],
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


# k-mer index section

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


def kmer_workloads(work, unitigs, rng):
    """Write the query files that every tool reads.

    They hold reads taken inside single unitigs (100% positive), the same reads with 1% substitutions and 0.4% indels, 1M random 31-mers, and the real SRR5833294 reads from another strain."""
    long_u = [u for _, u in unitigs if len(u) >= 150]
    pos, err = [], []
    for i in range(20000):
        u = rng.choice(long_u)
        p = rng.randrange(0, len(u) - 150 + 1)
        r = u[p : p + 150]
        r = rc(r) if rng.random() < 0.5 else r
        pos.append((f"p{i}", r))
        err.append((f"e{i}", mutate(r, rng, 0.01 * 0.6, 0.002, 0.002)))
    rand = [
        (f"x{i}", "".join(rng.choice("ACGT") for _ in range(31)))
        for i in range(1_000_000)
    ]
    srr = list(read_fasta(os.path.join(DATA, "queries/SRR5833294.10K.fastq.gz")))
    files = {}
    for name, recs in (
        ("reads_from_index", pos),
        ("reads_1pct_err", err),
        ("random_31mers", rand),
        ("SRR5833294", srr),
    ):
        files[name] = os.path.join(work, f"q_{name}.fa")
        write_fasta(files[name], recs)
    return files


def kmer_truth(set_db, path, d, threads, k=31):
    """(k-mer windows, positives) expected for a query file.

    Every tool reports all len - k + 1 windows of a read as num_kmers, with windows containing N counted negative.
    The positives are KMC's count of query windows whose k-mer is in the set (kmc_hits)."""
    total = sum(max(0, len(s) - k + 1) for _, s in read_fasta(path))
    return total, kmc_hits(set_db, path, k, d, threads)


def parse_counts(text):
    """(num_kmers, num_positive_kmers) from a query tool's report, None where absent."""
    k = re.search(r"num_kmers = (\d+)", text)
    p = re.search(r"num_positive_kmers = (\d+)", text)
    return (int(k.group(1)) if k else None, int(p.group(1)) if p else None)


def sbwt_counts(path):
    """(windows, positives) of SBWT's search output, one id per window with -1 for absent."""
    total = pos = 0
    with open(path) as f:
        for line in f:
            ids = line.split()
            total += len(ids)
            pos += sum(1 for x in ids if x != "-1")
    return total, pos


def fmsi_counts(path):
    """(windows, positives) of FMSI's query output, one 0/1 symbol per window."""
    total = pos = 0
    with open(path) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) == 2:
                total += len(t[1])
                pos += t[1].count("1")
    return total, pos


def section_kmer(bench, mode, work, threads, quick):
    """Build and query every k-mer index on SSHash's unitig sets and gate each row against KMC."""
    print(f"\n=== kmer ({mode})", flush=True)
    rng = random.Random(7)
    sets = ["ecoli"] if quick else list(KMER_SETS)
    for ds in sets:
        fname, sshash_m = KMER_SETS[ds]
        d = os.path.join(work, f"kmer_{ds}")
        os.makedirs(d, exist_ok=True)
        unitigs_fa = os.path.join(d, "unitigs.fa")
        if not os.path.exists(unitigs_fa):
            write_fasta(
                unitigs_fa, read_fasta(os.path.join(DATA, "unitigs_stitched", fname))
            )
        unitigs = list(read_fasta(unitigs_fa))
        # KMC gives the truth independently of every compared tool, namely the set's distinct k-mers and, per workload, its k-mer windows and how many of them KMC finds in the set.
        set_db = kmc_db(
            [unitigs_fa], 31, 1, os.path.join(d, "kmc_set"), threads, True, 4
        )
        n = kmc_count(set_db)
        wl = kmer_workloads(d, unitigs, rng)
        truth = {
            w: kmer_truth(set_db, p, os.path.join(d, f"kmc_{w}"), threads)
            for w, p in wl.items()
        }
        print(f"  {ds}: {n} distinct k-mers; truth {truth}", flush=True)
        m_par = 15
        bt = threads if mode == "optimal" else 1
        # Queries run on one thread in both modes.
        # SSHash, SBWT and FMSI have no query threads, so a multi-threaded optik row would not compare like with like.
        # optik's own -t scaling gets a separate row below.
        qt = 1
        tools = {
            "optik-spss": {
                "build": [
                    BIN["optik"],
                    "build-spss",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-m",
                    str(m_par if mode == "parity" else 0),
                    "-t",
                    str(bt),
                    "-o",
                    f"{d}/optik.spss",
                ],
                "files": [f"{d}/optik.spss"],
                "query": lambda q, d=d, qt=qt: (
                    [
                        BIN["optik"],
                        "query",
                        "-i",
                        f"{d}/optik.spss",
                        "-q",
                        q,
                        "-k",
                        "31",
                        "-t",
                        str(qt),
                    ],
                    None,
                ),
                "params": f"k=31 m={m_par if mode == 'parity' else 'auto'} build -t {bt}, query -t {qt}",
            },
            "optik-fm": {
                "build": [
                    BIN["optik"],
                    "build-fm",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-c",
                    "1",
                    "-t",
                    str(bt),
                    "-o",
                    f"{d}/optik.fm",
                ],
                "files": [f"{d}/optik.fm"],
                "query": lambda q, d=d, qt=qt: (
                    [
                        BIN["optik"],
                        "query",
                        "-i",
                        f"{d}/optik.fm",
                        "-q",
                        q,
                        "-k",
                        "31",
                        "-t",
                        str(qt),
                    ],
                    None,
                ),
                "params": f"k=31 build -t {bt}, query -t {qt}",
            },
            "sshash": {
                # The SSHash 2026 authors' command uses -g 16 (or 80% of this machine's memory).
                "build": [
                    BIN["sshash"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-m",
                    str(m_par if mode == "parity" else sshash_m),
                    "-g",
                    str(mem_budget_gb()),
                    "-t",
                    str(bt),
                    "-d",
                    d,
                    "-o",
                    f"{d}/idx.sshash",
                ],
                "files": [f"{d}/idx.sshash"],
                "query": lambda q, d=d: (
                    [BIN["sshash"], "query", "-i", f"{d}/idx.sshash", "-q", q],
                    None,
                ),
                "params": f"k=31 m={m_par if mode == 'parity' else sshash_m} build -t {bt}, query single-threaded",
            },
            "sbwt": {
                # The SSHash 2026 authors run SBWT with -m 16 (or 80% of this machine) and plain-matrix.
                "build": [
                    BIN["sbwt"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-o",
                    f"{d}/idx.sbwt",
                    "-k",
                    "31",
                    "--variant",
                    "plain-matrix",
                    "--add-reverse-complements",
                    "-m",
                    str(mem_budget_gb()),
                    "-t",
                    str(bt),
                    "--temp-dir",
                    d,
                ],
                "files": [f"{d}/idx.sbwt"],
                "query": lambda q, d=d: (
                    [
                        BIN["sbwt"],
                        "search",
                        "-i",
                        f"{d}/idx.sbwt",
                        "-q",
                        q,
                        "-o",
                        f"{d}/sbwt.out",
                    ],
                    f"{d}/sbwt.out",
                ),
                # Timed runs send the per-k-mer ids to /dev/null, where hyperfine also sends FMSI's stdout.
                # Neither tool then spends time on file writes the counters skip.
                "timed": lambda q, d=d: [
                    BIN["sbwt"],
                    "search",
                    "-i",
                    f"{d}/idx.sbwt",
                    "-q",
                    q,
                    "-o",
                    "/dev/null",
                ],
                "answers": "per-k-mer",
                "params": f"k=31 plain-matrix (+ reverse complements) -m {mem_budget_gb()} build -t {bt}, "
                "search single-threaded",
            },
            "fmsi": {
                # FMSI indexes a masked superstring, so its build is KmerCamel plus fmsi index.
                # Both modes use the FMSI settings of the SSHash 2026 authors' kmer_sets_benchmark, kmercamel compute -S for SPSS input, fmsi index, and queries with the kLCP array and max-one masks through -S -O.
                # Those are FMSI's own optimizations, not knobs to match.
                "build": [
                    "bash",
                    "-c",
                    (
                        f"{BIN['kmercamel']} compute -S -k 31 "
                        f"-M {d}/ms.fa -o /dev/null {unitigs_fa} && {BIN['fmsi']} index -k 31 {d}/ms.fa"
                    ),
                ],
                "files_glob": f"{d}/ms.fa.fmsi.",
                "query": lambda q, d=d: (
                    [
                        BIN["fmsi"],
                        "query",
                        "-S",
                        "-O",
                        "-q",
                        q,
                        "-k",
                        "31",
                        f"{d}/ms.fa",
                    ],
                    f"{d}/fmsi.out",
                ),
                "params": "k=31 KmerCamel greedy max-one superstring (-S unitig input) + fmsi query -S -O "
                "(kLCP streaming, max-one; the SSHash 2026 benchmark's settings)",
            },
        }

        # Dictionary workload mapping k-mers to dense ids with optik query --ids.
        # SBWT's search output is also per-k-mer ids, its colex ranks, whereas SSHash has no CLI for per-k-mer ids, since its Lookup exists only in the library and its benchmarks.
        def ids_check(path, qfile, d=d, n=n):
            """Count invalid entries of an optik --ids output for query file qfile.

            Each window's id must be present exactly where SBWT finds the k-mer, lie below n, and map one-to-one onto k-mers."""
            id_of, kmer_of, bad = {}, {}, 0
            lines = Path(path).read_text().split("\n")
            # Presence per window comes from SBWT's own search output on the same query file, an external tool writing one colex id per window and -1 when absent.
            sb = os.path.join(d, "sbwt_ids.out")
            subprocess.run(
                [BIN["sbwt"], "search", "-i", f"{d}/idx.sbwt", "-q", qfile, "-o", sb],
                check=True,
                capture_output=True,
            )
            present = [
                [int(x) != -1 for x in text_line.split()]
                for text_line in Path(sb).read_text().split("\n")
            ]
            for (_, sq), line, pres in zip(
                read_fasta(qfile), lines, present, strict=False
            ):
                vals = list(map(int, line.split())) if line else []
                if len(vals) != max(0, len(sq) - 30) or len(pres) != len(vals):
                    bad += 1
                    continue
                for i, v in enumerate(vals):
                    cc = kmer_codes(sq[i : i + 31], 31)
                    c = cc[0] if cc else None
                    if (
                        pres[i] != (v >= 0)
                        or v >= 0
                        and (
                            v >= n
                            or id_of.setdefault(c, v) != v
                            or kmer_of.setdefault(v, c) != c
                        )
                    ):
                        bad += 1
            return bad

        for tool, t in tools.items():
            # Inputs are pre-read before each timed run (measure.py).
            r = bench.measure(
                f"{mode}/{ds}/{tool}/build", t["build"], cwd=d, inputs=[unitigs_fa]
            )
            if not r.ok:
                bench.add(
                    section="kmer",
                    mode=mode,
                    dataset=ds,
                    tool=tool,
                    phase="build",
                    gate=False,
                    gate_detail=f"build exit {r.code}: {r.err[-300:]}",
                    params=t["params"],
                )
                continue
            if tool == "optik-spss" and mode != "parity":
                # With -m 0 optik picks m from the input, and the row names the m it picked.
                mm = re.search(
                    r"minimizer_len\t(\d+)",
                    Path(r.out).read_text() if os.path.exists(r.out) else "",
                )
                t["params"] = t["params"].replace(
                    "m=auto", f"m=auto ({mm.group(1) if mm else '?'})"
                )
            idx = (
                list(t["files"])
                if "files" in t
                else [
                    os.path.join(d, f)
                    for f in os.listdir(d)
                    if os.path.join(d, f).startswith(t["files_glob"])
                ]
            )
            size = sum(os.path.getsize(f) for f in idx)
            bench.add(
                section="kmer",
                mode=mode,
                dataset=ds,
                tool=tool,
                phase="build",
                build_s=r.wall_row(),
                build_rss_mb=round(r.rss(), 1),
                index_mb=round(size / 2**20, 3),
                bits_per_kmer=round(8 * size / n, 3),
                gate=True,
                gate_detail="",
                params=t["params"],
            )
            per_kmer = (
                t.get("answers") or ("per-k-mer" if tool == "fmsi" else "counts")
            ) == "per-k-mer"
            for w, q in wl.items():
                cmd, outf = t["query"](q)
                capt = outf or os.path.join(d, "q.out")
                if per_kmer:
                    # One untimed run writes the per-k-mer answers to a file for the gate.
                    # Timed runs write them to /dev/null, so their timings carry no output writes.
                    rc_, _, _, err = bench.once(
                        f"{mode}/{ds}/{tool}/{w} (answers)",
                        cmd,
                        cwd=d,
                        stdout=None if tool in ("sbwt",) else capt,
                    )
                    r = (
                        bench.measure(
                            f"{mode}/{ds}/{tool}/{w}",
                            t["timed"](q) if "timed" in t else cmd,
                            cwd=d,
                            inputs=[q, *idx],
                        )
                        if rc_ == 0
                        else None
                    )
                else:
                    r = bench.measure(
                        f"{mode}/{ds}/{tool}/{w}",
                        cmd,
                        cwd=d,
                        stdout=capt,
                        inputs=[q, *idx],
                    )
                    rc_, err = r.code, r.err
                if rc_ != 0 or not r.ok:
                    bench.add(
                        section="kmer",
                        mode=mode,
                        dataset=ds,
                        tool=tool,
                        phase=w,
                        gate=False,
                        gate_detail=f"query exit {rc_ or r.code}: {(err or r.err)[-300:]}",
                        params=t["params"],
                    )
                    continue
                if tool == "sbwt":
                    got = sbwt_counts(outf)
                elif tool == "fmsi":
                    got = fmsi_counts(capt)
                else:
                    got = parse_counts(Path(capt).read_text() + err)
                if tool == "fmsi":
                    # FMSI prints one symbol per base rather than per k-mer for a read containing N, so only its positives can be compared.
                    ok = got[1] == truth[w][1]
                else:
                    ok = got == truth[w]
                nk = truth[w][0]
                bench.add(
                    section="kmer",
                    mode=mode,
                    dataset=ds,
                    tool=tool,
                    phase=w,
                    answers=t.get("answers")
                    or ("per-k-mer" if tool == "fmsi" else "counts"),
                    query_s=r.wall_row(),
                    ns_per_kmer=round(r.wall["mean"] * 1e9 / nk, 1),
                    query_rss_mb=round(r.rss(), 1),
                    kmers=got[0],
                    positive=got[1],
                    gate=ok,
                    gate_detail="" if ok else f"counts {got} != truth {truth[w]}",
                    params=t["params"],
                )
                if (
                    mode == "optimal"
                    and threads > 1
                    and tool in ("optik-spss", "optik-fm")
                ):
                    # optik's own query scaling runs the same query on all threads and must give the same answers.
                    cmd_t = [
                        x if x != str(qt) or cmd[i - 1] != "-t" else str(threads)
                        for i, x in enumerate(cmd)
                    ]
                    r = bench.measure(
                        f"{mode}/{ds}/{tool}/{w}/t{threads}",
                        cmd_t,
                        cwd=d,
                        stdout=capt,
                        inputs=[q, *idx],
                    )
                    got_t = (
                        parse_counts(Path(capt).read_text() + r.err) if r.ok else None
                    )
                    ok_t = got_t is not None and got_t[1] == truth[w][1]
                    bench.add(
                        section="kmer",
                        mode=mode,
                        dataset=ds,
                        tool=f"{tool} -t {threads} (own scaling)",
                        phase=w,
                        answers="counts",
                        query_s=r.wall_row(),
                        ns_per_kmer=round(r.wall["mean"] * 1e9 / nk, 1)
                        if r.ok
                        else None,
                        query_rss_mb=round(r.rss(), 1) if r.ok else None,
                        kmers=got_t[0] if got_t else None,
                        positive=got_t[1] if got_t else None,
                        gate=ok_t,
                        gate_detail=""
                        if ok_t
                        else f"counts {got_t} != truth {truth[w]}",
                        params=t["params"]
                        + f"; query -t {threads}: not a comparison row (others cannot thread)",
                    )
        # optik's dictionary ids are checked once every index exists, with presence judged by SBWT's search.
        if os.path.exists(f"{d}/optik.spss") and os.path.exists(f"{d}/idx.sbwt"):
            t = tools["optik-spss"]
            q = wl["reads_1pct_err"]
            cmd = [
                BIN["optik"],
                "query",
                "-i",
                f"{d}/optik.spss",
                "-q",
                q,
                "-k",
                "31",
                "-t",
                str(qt),
                "--ids",
                f"{d}/ids.txt",
            ]
            rc_, _, _, err = bench.once(
                f"{mode}/{ds}/optik-spss/ids (answers)",
                cmd,
                cwd=d,
                stdout=os.path.join(d, "q.out"),
            )
            bad = ids_check(f"{d}/ids.txt", q) if rc_ == 0 else -1
            # Timed with the ids sent to /dev/null, like the other per-k-mer answers.
            r = bench.measure(
                f"{mode}/{ds}/optik-spss/ids",
                cmd[:-1] + ["/dev/null"],
                cwd=d,
                inputs=[q, f"{d}/optik.spss"],
            )
            bench.add(
                section="kmer",
                mode=mode,
                dataset=ds,
                tool="optik-spss --ids",
                phase="reads_1pct_err",
                answers="per-k-mer",
                query_s=r.wall_row(),
                ns_per_kmer=round(r.wall["mean"] * 1e9 / truth["reads_1pct_err"][0], 1)
                if r.ok
                else None,
                query_rss_mb=round(r.rss(), 1) if r.ok else None,
                gate=bad == 0 and r.ok,
                gate_detail=""
                if bad == 0
                else f"{bad} id inconsistencies (presence vs SBWT, range, bijection)",
                params=t["params"]
                + " --ids (dense k-mer ids, SSHash Lookup semantics)",
            )


# De Bruijn graph section


def gfa_or_fasta_seqs(path):
    """Segment sequences of a GFA file, or the records of a FASTA file."""
    op = gzip.open if path.endswith(".gz") else open
    if ".gfa" in path:
        with op(path, "rt") as f:
            return [line.split("\t")[2] for line in f if line.startswith("S\t")]
    return [s for _, s in read_fasta(path)]


def section_dbg(bench, mode, work, threads, quick):
    """Build compacted de Bruijn graphs with every graph tool and judge each against KMC."""
    print(f"\n=== dbg ({mode})", flush=True)
    d = os.path.join(work, "dbg")
    os.makedirs(d, exist_ok=True)
    genome = os.path.join(T, "salmonella_LT2.fa")
    reads = os.path.join(d, "reads30x.fa")
    if not os.path.exists(reads):
        subprocess.run(
            [
                sys.executable,
                os.path.join(HERE, "gen_reads.py"),
                "-r",
                genome,
                "-c",
                "30",
                "-l",
                "150",
                "-e",
                "0.005",
                "--ins",
                "0.0005",
                "--del",
                "0.0005",
                "--rc",
                "0.5",
                "-s",
                "31",
                "-o",
                reads,
            ],
            check=True,
        )
    cases = [("genome", genome, 1)] + ([] if quick else [("reads_30x", reads, 2)])
    for name, inp, c in cases:
        reads_in = name.startswith("reads")
        o = os.path.join(d, f"{name}_{mode}")
        os.makedirs(o, exist_ok=True)
        truth = kmc_db([inp], 31, c, os.path.join(o, "truth"), threads, True, 4)
        print(f"  {name}: KMC solid k-mers (c >= {c}) = {kmc_count(truth)}", flush=True)
        cf_mode = "--read" if reads_in else "--ref"
        tools = {
            # All tools write the same artifact, uncompressed FASTA unitigs without GFA links or gzip, so the timings compare graph construction rather than output formats.
            "optik": (
                [
                    BIN["optik"],
                    "dbruijn",
                    "-i",
                    inp,
                    "-k",
                    "31",
                    "-c",
                    str(c),
                    "-t",
                    str(threads),
                    "-o",
                    f"{o}/optik.fa",
                ],
                f"{o}/optik.fa",
                f"k=31 -c {c} -t {threads}; FASTA unitigs",
            ),
            "cuttlefish": (
                [
                    BIN["cuttlefish"],
                    "build",
                    "-s",
                    inp,
                    "-k",
                    "31",
                    "-t",
                    str(threads),
                    "-o",
                    f"{o}/cf",
                    cf_mode,
                ]
                + (["-c", str(c)] if cf_mode == "--read" else []),
                f"{o}/cf.fa",
                f"k=31 {cf_mode}"
                + (f" -c {c} (edge cutoff)" if cf_mode == "--read" else "")
                + f" -t {threads}",
            ),
            # For Bifrost, -r keeps every k-mer and -s drops k-mers seen once through its Bloom filter.
            "bifrost": (
                [
                    BIN["bifrost"],
                    "build",
                    "-r" if c == 1 else "-s",
                    inp,
                    "-k",
                    "31",
                    "-t",
                    str(threads),
                    "-o",
                    f"{o}/bf",
                    "-f",
                    "-n",
                    "-N",
                ],
                f"{o}/bf.fasta",
                (
                    f"k=31 {'-r (all k-mers)' if c == 1 else '-s (drop singletons)'} -t {threads}; "
                    "-f -n -N (FASTA, uncompressed, no index file)"
                ),
            ),
        }
        written = {}
        for tool, (cmd, outp, params) in tools.items():
            # The KMC judging evicts the input from the cache, so it is read in again before every run.
            r = bench.measure(f"{mode}/dbg/{name}/{tool}", cmd, cwd=o, inputs=[inp])
            rc_, err = (0 if r.ok else r.code), r.err
            if rc_ != 0 or not os.path.exists(outp):
                alt = outp.replace(".gfa.gz", ".gfa")
                if rc_ == 0 and os.path.exists(alt):
                    outp = alt
                else:
                    bench.add(
                        section="dbg",
                        mode=mode,
                        dataset=name,
                        tool=tool,
                        gate=False,
                        gate_detail=f"exit {rc_}, output {outp} missing: {err[-300:]}",
                        params=params,
                    )
                    continue
            seqs = gfa_or_fasta_seqs(outp)
            fa = os.path.join(o, f"{tool}.unitigs.fa")
            write_fasta(fa, ((str(i), x) for i, x in enumerate(seqs)))
            written[tool] = fa
            # KMC judges each tool under that tool's own documented graph model.
            ok, det, nd, n_model = judge_unitigs(
                tool, fa, [inp], 31, c, reads_in, o, threads, 4, truth
            )
            lens = sorted((len(x) for x in seqs), reverse=True)
            half, acc, n50 = sum(lens) / 2, 0, 0
            for read_len in lens:
                acc += read_len
                if acc >= half:
                    n50 = read_len
                    break
            bench.add(
                section="dbg",
                mode=mode,
                dataset=name,
                tool=tool,
                wall_s=r.wall_row(),
                rss_mb=round(r.rss(), 1),
                unitigs=len(lens),
                total_bp=sum(lens),
                n50=n50,
                kmers=nd,
                truth_kmers=n_model,
                graph="edge-centric" if tool == "cuttlefish" else "node-centric",
                gate=ok,
                gate_detail=det,
                params=params,
            )
        # Maximal unitigs are unique, so on the genome optik's must equal Cuttlefish's string for string.
        if (
            not reads_in
            and "optik" in written
            and "cuttlefish" in written
            and not same_unitigs(written["optik"], written["cuttlefish"])
        ):
            for r in bench.rows:
                if (
                    r.get("section"),
                    r.get("dataset"),
                    r.get("mode"),
                    r.get("tool"),
                ) == (
                    "dbg",
                    name,
                    mode,
                    "optik",
                ):
                    r["gate"] = False
                    r["gate_detail"] = (
                        r["gate_detail"] + "; " if r["gate_detail"] else ""
                    ) + "unitigs differ from Cuttlefish's (maximal unitigs are unique)"


# Mapping section


def _ref_span(t):
    """Reference bases spanned by the CIGAR of SAM fields t."""
    return sum(
        int(x) for x, o in re.findall(r"(\d+)([MIDNSHP=X])", t[5]) if o in "MDN=X"
    )


def sam_check(sam, ref, truth, interval=False, reads=None):
    """Validity and accuracy of a SAM against the reference and the truth.

    A read is placed correctly on the right contig and strand when, for short reads, POS - 1 - leading clip lies within max(20, L/10) of the origin.
    With interval=True (long or chimeric reads, whose clips may be foreign sequence) the aligned reference interval must instead overlap the truth interval by at least half its length.

    `reads` is given for optik's rows and selects strict mode.
    optik always writes SEQ and NM:i and never hard-clips, so a mapped record lacking SEQ or NM, holding H, or whose SEQ is not the read (reverse-complemented for FLAG 16) is invalid.
    The lenient rules below exist only for the documented formats of other tools.

    Two SAM conventions are not output errors for the other tools.
    wfmash writes every mapped segment of a read as a primary record, none flagged supplementary, and omits unmapped reads.
    For those tools both are counted in spec_notes, and the read is judged by its longest primary alignment, the rule of minimap2's paftools.js mapeval (-m 0), which does not depend on output order.
    For optik (strict) they stay invalid, since optik writes one primary per read and every read."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    seqs = {}
    for n, s in read_fasta(ref):
        seqs[n] = s.upper()
    prim = {}
    bad = []
    spec = []
    strict = reads is not None
    nm_tot = aln_tot = 0
    with open(sam) as f:
        for line in f:
            if line[0] == "@":
                continue
            t = line.rstrip("\n").split("\t")
            fl = int(t[1])
            if not fl & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                if t[0] in prim:
                    (bad if strict else spec).append(f"{t[0]}: two primaries")
                    if not strict and _ref_span(t) <= _ref_span(prim[t[0]]):
                        continue  # judge the read by its longest primary
                prim[t[0]] = t
            if fl & SAM_FLAG_UNMAPPED:
                continue
            ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
            qlen = sum(int(x) for x, o in ops if o in "MIS=X")
            if reads is not None:
                rd = reads.get(t[0])
                want = (
                    None
                    if rd is None
                    else (rd.translate(comp)[::-1] if fl & SAM_FLAG_REVERSE else rd)
                )
                if (
                    t[9] == "*"
                    or any(o == "H" for _, o in ops)
                    or not any(y.startswith("NM:i:") for y in t[11:])
                    or want is None
                    or t[9].upper() != want.upper()
                ):
                    bad.append(f"{t[0]}: missing SEQ/NM, hard clip, or SEQ != read")
                    continue
            if t[9] != "*" and qlen != len(t[9]):
                bad.append(f"{t[0]}: CIGAR consumes {qlen} != SEQ {len(t[9])}")
                continue
            if t[9] == "*" or any(o == "H" for _, o in ops):
                continue  # hard-clipped or without SEQ, so NM cannot be recomputed
            r = seqs[t[2]]
            rp, qp, nm, al = int(t[3]) - 1, 0, 0, 0
            for x, o in ops:
                x = int(x)
                if o in "M=X":
                    seg = t[9][qp : qp + x].upper()
                    nm += sum(
                        1 for i in range(x) if seg[i] != r[rp + i] or seg[i] == "N"
                    )
                    qp += x
                    rp += x
                    al += x
                elif o == "I":
                    nm += x
                    qp += x
                    al += x
                elif o == "D":
                    nm += x
                    rp += x
                elif o == "S":
                    qp += x
            tag = [y for y in t[11:] if y.startswith("NM:i:")]
            if tag and int(tag[0][5:]) != nm:
                bad.append(f"{t[0]}: NM:i {tag[0][5:]} != recomputed {nm}")
            nm_tot += nm
            aln_tot += al
    missing = [n for n in truth if n not in prim]
    if missing:
        (bad if strict else spec).append(
            f"{len(missing)} reads without a primary record"
        )
    ok_reads = mapped = wrong = hq_wrong = 0
    for n, g in truth.items():
        t = prim.get(n)
        if not t or int(t[1]) & SAM_FLAG_UNMAPPED:
            continue
        mapped += 1
        if g is None:
            continue
        rname, p0, strand, read_len = g
        ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
        same = t[2] == rname and ((int(t[1]) & SAM_FLAG_REVERSE) > 0) == (strand == "-")
        if interval:
            a = int(t[3]) - 1
            b = a + sum(int(x) for x, o in ops if o in "MD=X")
            good = (
                same and b > a and min(b, p0 + read_len) - max(a, p0) >= 0.5 * (b - a)
            )
        else:
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = same and abs(int(t[3]) - 1 - lead - p0) <= max(20, read_len // 10)
        ok_reads += good
        if not good:
            wrong += 1
            hq_wrong += int(t[4]) >= HIGH_MAPQ
    n = len(truth)
    return {
        "valid": not bad,
        "invalid_examples": bad[:3],
        "invalid_count": len(bad),
        "spec_notes": spec[:3],
        "spec_count": len(spec),
        "mapped_pct": round(100 * mapped / n, 3),
        "correct_pct": round(100 * ok_reads / n, 3),
        "wrong": wrong,
        "wrong_mapq30": hq_wrong,
        "nm_per_aligned_base": round(nm_tot / aln_tot, 4) if aln_tot else None,
    }


def write_fai(path):
    """Write a samtools-compatible .fai with name, length, offset, bases per line and bytes per line."""
    out = []
    name = None
    with open(path, "rb") as f:
        data = f.read()
    pos = 0
    n = len(data)
    while pos < n:
        nl = data.find(b"\n", pos)
        nl = n if nl < 0 else nl
        line = data[pos:nl]
        if line.startswith(b">"):
            tok = line[1:].split()
            name = tok[0].decode() if tok else f"seq{len(out)}"
            out.append([name, 0, nl + 1, 0, 0])
        elif name is not None:
            if out[-1][3] == 0:
                out[-1][3] = len(line)
                out[-1][4] = len(line) + 1
            out[-1][1] += len(line)
        pos = nl + 1
    with open(path + ".fai", "w") as f:
        f.writelines("\t".join(map(str, r)) + "\n" for r in out)


def load_truth(tsv):
    """Map each read name to (contig, position, strand, length) from a truth TSV."""
    gt = {}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != "read_name":
                gt[f[0]] = (f[1], int(f[2]), f[3], int(f[4]))
    return gt


def section_map(bench, mode, work, threads, quick):
    """Map simulated read sets with every mapper and score validity and accuracy against truth."""
    print(f"\n=== map ({mode})", flush=True)
    d = os.path.join(work, "map")
    os.makedirs(d, exist_ok=True)
    sal = os.path.join(T, "salmonella_LT2.fa")
    gen = os.path.join(HERE, "gen_reads.py")
    sets = {
        # name -> (reference, reads, truth, kind)
        "salmonella_150bp": (
            sal,
            {
                "n": 50000,
                "length": 150,
                "e": 0.006,
                "ins": 0.002,
                "dele": 0.002,
                "s": 2,
            },
            "short",
        ),
        "salmonella_ONT_20kb": (
            sal,
            {
                "n": 400,
                "length": 20000,
                "sd": 0.5,
                "e": 0.04,
                "ins": 0.03,
                "dele": 0.03,
                "s": 7,
            },
            "ont",
        ),
        "salmonella_HiFi_15kb": (
            sal,
            {
                "n": 400,
                "length": 15000,
                "sd": 0.3,
                "e": 0.002,
                "ins": 0.001,
                "dele": 0.001,
                "s": 8,
            },
            "hifi",
        ),
    }
    if not quick:
        syn = os.path.join(work, "syn100m.fa")
        if not os.path.exists(syn):
            r = random.Random(100)
            write_fasta(
                syn,
                [
                    (f"chr{i}", "".join(r.choice("ACGT") for _ in range(25_000_000)))
                    for i in range(4)
                ],
            )
        sets["synthetic100M_150bp"] = (
            syn,
            {"n": 300000, "length": 150, "e": 0.01, "s": 1},
            "short",
        )
    files = {}
    for name, (ref, p, kind) in sets.items():
        fa, tsv = os.path.join(d, f"{name}.fa"), os.path.join(d, f"{name}.tsv")
        if not os.path.exists(tsv):
            cmd = [
                sys.executable,
                gen,
                "-r",
                ref,
                "-n",
                str(p["n"]),
                "-l",
                str(p["length"]),
                "-e",
                str(p["e"]),
                "--rc",
                "0.5",
                "-s",
                str(p["s"]),
                "-o",
                fa,
                "--gt",
                tsv,
            ]
            if "ins" in p:
                cmd += ["--ins", str(p["ins"]), "--del", str(p["dele"])]
            if "sd" in p:
                cmd += ["--len-sd", str(p["sd"])]
            subprocess.run(cmd, check=True)
        files[name] = (ref, fa, load_truth(tsv), kind)
    # Chimeric long reads are ONT reads with 5-20 kb of foreign sequence appended.
    ch_fa = os.path.join(d, "salmonella_ONT_chimeric.fa")
    ont_ref, ont_fa, ont_truth, _ = files["salmonella_ONT_20kb"]
    if not os.path.exists(ch_fa):
        r = random.Random(11)
        write_fasta(
            ch_fa,
            [
                (
                    n,
                    s
                    + "".join(
                        r.choice("ACGT") for _ in range(r.randrange(5000, 20000))
                    ),
                )
                for n, s in read_fasta(ont_fa)
            ],
        )
    files["salmonella_ONT_chimeric"] = (ont_ref, ch_fa, ont_truth, "ont")
    # Negative control, random reads that must not map
    neg = os.path.join(d, "random_150bp.fa")
    if not os.path.exists(neg):
        r = random.Random(3)
        write_fasta(
            neg,
            [
                (f"x{i}", "".join(r.choice("ACGT") for _ in range(150)))
                for i in range(50000)
            ],
        )
    files["random_150bp"] = (sal, neg, {f"x{i}": None for i in range(50000)}, "short")
    for ref, fa, _truth, _kind in files.values():
        for f in (ref, fa):  # wfmash needs .fai indexes
            if not os.path.exists(f + ".fai"):
                write_fai(f)
    preset = {"short": "sr", "ont": "map-ont", "hifi": "map-hifi"}
    # k and w of those presets, from minimap2 v2.31 options.c.
    # mm_idxopt_init sets k = 15, w = 10 for map-ont, the default, map-hifi uses k = 19, w = 19, and sr k = 21, w = 11.
    mm2_preset_kw = {"sr": (21, 11), "map-ont": (15, 10), "map-hifi": (19, 19)}
    thread_set = [1, threads] if mode == "parity" else [threads]
    for name, (ref, fa, truth, kind) in files.items():
        for t in thread_set:
            # In both modes optik indexes the same minimizers as minimap2's preset, with its k and w, and each mapper keeps its own other settings.
            pk, pw = mm2_preset_kw[preset[kind]]
            tools = {
                "optik": (
                    [
                        BIN["optik"],
                        "map",
                        "-r",
                        ref,
                        "-q",
                        fa,
                        "-t",
                        str(t),
                        "-k",
                        str(pk),
                        "-w",
                        str(pw),
                    ],
                    f"-k {pk} -w {pw} (minimap2 -ax {preset[kind]}'s seeds)",
                ),
            }
            mm = [BIN["minimap2"], "-ax", preset[kind], "-t", str(t)]
            tools["minimap2"] = (
                mm + [ref, fa],
                f"-ax {preset[kind]} (preset: k {pk}, w {pw})",
            )
            if kind != "short":
                # Wavefront-based mapper with WFA alignment.
                # Parity uses the read mode -S 1 -p 70 that run_validations validates, and optimal its defaults.
                # -k sets the same k-mer length, while its -w is a sketch window in bases rather than a minimizer window and keeps its default.
                ro = ["-S", "1", "-p", "70"] if mode == "parity" else []
                tools["wfmash"] = (
                    [BIN["wfmash"], ref, fa, "-a", "-k", str(pk), "-t", str(t)] + ro,
                    f"-a -k {pk} "
                    + ("-S 1 -p 70 (read mode)" if ro else "(otherwise defaults)"),
                )
            for tool, (cmd, params) in tools.items():
                sam = os.path.join(d, f"{name}_{tool}_{mode}_t{t}.sam")
                # Each mapper writes its SAM in every timed run, so all do the same work, and the gate reads the SAM of the last measured run.
                r = bench.measure(
                    f"{mode}/map/{name}/{tool}/t{t}",
                    cmd,
                    cwd=d,
                    stdout=sam,
                    inputs=[ref, fa],
                )
                rc_, err = (0 if r.ok else r.code), r.err
                if rc_ != 0:
                    bench.add(
                        section="map",
                        mode=mode,
                        dataset=name,
                        tool=tool,
                        threads=t,
                        gate=False,
                        gate_detail=f"exit {rc_}: {err[-300:]}",
                        params=params,
                    )
                    continue
                chk = sam_check(
                    sam,
                    ref,
                    truth,
                    interval=kind != "short",
                    reads=dict(read_fasta(fa)) if tool == "optik" else None,
                )
                gate = chk["valid"]
                det = (
                    ""
                    if gate
                    else f"{chk['invalid_count']} invalid: {chk['invalid_examples']}"
                )
                spec_note = ""
                if chk["spec_count"]:
                    multi = sum(1 for x in chk["spec_notes"] if "two primaries" in x)
                    spec_note = (
                        f"{chk['spec_count']} SAM-convention notes, e.g. {chk['spec_notes']}"
                        + (
                            " (one primary per mapped segment, none flagged supplementary; "
                            "the longest is judged)"
                            if multi
                            else ""
                        )
                    )
                # The same SAM also goes through Picard ValidateSamFile and samtools calmd.
                ext_ok, ext = external_sam_checks(sam, ref)
                if not ext_ok:
                    gate = False
                    det = (det + "; " if det else "") + ext
                if name == "random_150bp" and chk["mapped_pct"] > 0.1:
                    # Random sequence must not map, so this row FAILS.
                    # It used to only append a note and still pass.
                    gate = False
                    det += f" negative control: {chk['mapped_pct']}% of random reads mapped"
                bench.add(
                    section="map",
                    mode=mode,
                    dataset=name,
                    tool=tool,
                    threads=t,
                    wall_s=r.wall_row(),
                    rss_mb=round(r.rss(), 1),
                    mapped_pct=chk["mapped_pct"],
                    correct_pct=chk["correct_pct"],
                    wrong=chk["wrong"],
                    wrong_mapq30=chk["wrong_mapq30"],
                    nm_rate=chk["nm_per_aligned_base"],
                    sam_spec=spec_note or "conforms",
                    gate=gate,
                    gate_detail=det,
                    params=params,
                )


# Pairwise alignment section


def section_align(bench, mode, work, quick):
    """Align random pairs over a length and error grid with optik and WFA2-lib and compare scores."""
    print(f"\n=== align ({mode})", flush=True)
    d = os.path.join(work, "align")
    os.makedirs(d, exist_ok=True)
    exe = build_wfab(d)
    grid = [
        (150, 0.01, 20000),
        (150, 0.05, 20000),
        (1000, 0.05, 2000),
        (1000, 0.15, 500),
        (10000, 0.05, 40),
        (10000, 0.15, 10),
        (100000, 0.05, 2),
    ]
    if quick:
        grid = grid[:3]
    methods = ["optik-wfa", "optik-biwfa", "wfa2-high", "wfa2-ultralow"]
    for read_len, err, npairs in grid:
        rng = random.Random(read_len * 13 + int(err * 1000))
        pf = os.path.join(d, f"pairs_{read_len}_{err}.tsv")
        if not os.path.exists(pf):
            with open(pf, "w") as f:
                for _ in range(npairs):
                    p = "".join(rng.choice("ACGT") for _ in range(read_len))
                    f.write(
                        p
                        + "\t"
                        + mutate(p, rng, err * 0.6, err * 0.2, err * 0.2)
                        + "\n"
                    )
        scores = {}
        for m in methods:
            out = os.path.join(d, f"scores_{m}_{read_len}_{err}.txt")
            res = os.path.join(d, "res.txt")
            # The harness times its own alignments, BiWFA's measure, which leaves out input parsing.
            # Warmup runs are dropped from the statistics.
            runs_ = bench.repeat(
                f"{mode}/align/{read_len}/{err}/{m}", [exe, m, pf, out], stdout=res
            )
            if not runs_ or runs_[-1][0] != 0:
                rc_, _, e = runs_[-1] if runs_ else (-1, 0, "no run")
                bench.add(
                    section="align",
                    mode=mode,
                    length=read_len,
                    error=err,
                    method=m,
                    gate=False,
                    gate_detail=f"exit {rc_} (memory cap {bench.cap >> 30} GiB?): {e[-200:]}",
                )
                continue
            secs = [float(x[2].split()[1]) for x in runs_]
            invalid = max(int(x[2].split()[3]) for x in runs_)
            scores[m] = Path(out).read_text().split()
            us = [x * 1e6 / npairs for x in secs]
            mean = sum(us) / len(us)
            sd = (
                (sum((x - mean) ** 2 for x in us) / (len(us) - 1)) ** 0.5
                if len(us) > 1
                else None
            )
            bench.add(
                section="align",
                mode=mode,
                length=read_len,
                error=err,
                method=m,
                pairs=npairs,
                us_per_pair=round(mean, 2),
                us_per_pair_sd=round(sd, 2) if sd is not None else None,
                runs=len(us),
                warmup=bench.warmup,
                rss_mb=round(max(x[1] for x in runs_), 1),
                gate=invalid == 0,
                gate_detail="" if invalid == 0 else f"{invalid} CIGARs do not re-score",
            )
        ref = next(iter(scores.values()), None)
        agree = all(v == ref for v in scores.values())
        for r in bench.rows:
            if (
                r.get("section") == "align"
                and r.get("length") == read_len
                and r.get("error") == err
                and r.get("mode") == mode
            ):
                r["scores_agree"] = agree
                if not agree:
                    r["gate"] = False
                    r["gate_detail"] = (
                        r.get("gate_detail") or ""
                    ) + " optimal scores differ between methods"


# MAPQ calibration section


def section_mapq(bench, mode, work, threads, quick):
    """Run bench/mapq_calibration.py at several error rates and report over-confident MAPQ bins."""
    print(f"\n=== mapq ({mode})", flush=True)
    d = os.path.join(work, "mapq")
    for args in (
        ["--error", "0.01"],
        ["--error", "0.04"],
        ["--error", "0.08"],
        ["--ref", os.path.join(T, "salmonella_LT2.fa"), "--error", "0.04"],
    ):
        if quick and args[1] != "0.04":
            continue
        n = "200000" if quick else "1000000"
        p = subprocess.run(
            [
                sys.executable,
                os.path.join(HERE, "mapq_calibration.py"),
                # Its rows join this bench's section, so it writes no section of its own.
                "--report",
                "",
                "--reads",
                n,
                "--threads",
                str(threads),
                "--out",
                d,
            ]
            + args,
            capture_output=True,
            text=True,
            check=False,
        )
        bench.log(f"{mode}/mapq", ["mapq_calibration.py", "--reads", n] + args)
        txt = p.stdout
        for tool in ("optik", "minimap2"):
            m = re.search(
                rf"^{tool}:.*?MAPQ>=30: (\d+) reads, (\d+) wrong",
                txt,
                re.DOTALL | re.MULTILINE,
            )
            block = (
                txt.split(f"{tool}:")[1].split("minimap2:")[0]
                if f"{tool}:" in txt
                else ""
            )
            over = block.count("OVERCONFIDENT")
            # A missing report, whether the calibration crashed or its format changed, is a failure, not "0 overconfident bins".
            parsed = p.returncode == 0 and m is not None and bool(block)
            bench.add(
                section="mapq",
                mode=mode,
                setting=" ".join(a.replace(T + "/", "") for a in args),
                tool=tool,
                mapq30_reads=int(m.group(1)) if m else None,
                mapq30_wrong=int(m.group(2)) if m else None,
                overconfident_bins=over,
                calibrated=(
                    "yes"
                    if over == 0
                    else f"no: {over} bins over-confident (wrong > 2x implied + 3, Poisson tail < 1e-3)"
                )
                if parsed
                else "n/a",
                # The gate only requires that the calibration ran and was parsed.
                # Whether a mapper is calibrated is the result, reported the same way for every tool.
                gate=parsed,
                gate_detail=""
                if parsed
                else (
                    f"mapq_calibration.py produced no report (rc {p.returncode}): "
                    f"{p.stderr.strip()[-200:]}"
                ),
            )


# Report


def report(bench, out, env, dest=""):
    """Write results.json and the Markdown tables of a run to `out`, and the run's section of the sectioned results file `dest`."""
    md = [
        f"# optik vs published tools, small sets ({env['date']})",
        "",
        (
            f"- CPU: {env['cpu']} ({env['physical_cores']} cores, {env['cores_threads']}), {env['memory']}, "
            f"kernel {env['kernel']}, THP {env['thp']}, governor {env['governor']}"
        ),
        f"- {env['gcc']}; {env['rustc']}; optik {env['optik_git']}",
        f"- Builds: {env['build_flags']}",
        "- Tool revisions: " + ", ".join(f"{k} {v}" for k, v in env["tools"].items()),
        (
            f"- Times: hyperfine whole-process mean ± sd over {bench.runs} measured runs after {bench.warmup} discarded "
            "warmup run(s), inputs read into the page cache before each (bench/measure.py); RSS: peak over the "
            "measured runs (wait4). Alignment: the harness's own timer, same warmup/runs."
        ),
        (
            "- Gate: PASS = output verified (see compare_tools_small_sets.py docstring); FAIL rows are wrong "
            "output, not just slow. Each tool is judged under its own documented model (Cuttlefish: cutoff on "
            "(k+1)-mer edges). `sam_spec` (SAM conventions such as wfmash's one primary per segment) and "
            "`calibrated` (MAPQ over-confidence) are measured results reported for every tool, not gates."
        ),
        "",
    ]
    md += [
        "## Use-case capabilities (documented features, not measured)",
        "",
        "| k-mer index | membership + streaming | dense ids (Lookup) | Access(id) | weights | add/remove after build |",
        "|---|---|---|---|---|---|",
        "| optik build-spss | yes | yes (`query --ids`) | yes (library) | no | yes (`spss-update`) |",
        "| SSHash | yes | yes (library / `bench`) | yes | yes (`--weighted`) | no |",
        "| SBWT | yes | colex ranks (`search` output) | yes | no | no |",
        "| FMSI | yes | no | no | no | set operations on masked superstrings |",
        "| optik build-fm | yes | no | no | no | no |",
        "",
        "| dBG | graph model | min-count filter | output |",
        "|---|---|---|---|",
        "| optik dbruijn | edge-centric (observed (k+1)-mers) | exact k-mer counts | GFA1 (S + L) |",
        "| Cuttlefish 2/3 | edge-centric | exact, on (k+1)-mer edges (`--read`) | FASTA unitigs |",
        "| Bifrost | node-centric | singletons only, Bloom filter (approximate) | GFA1 |",
        "",
        "| mapper | output | short reads | long reads |",
        "|---|---|---|---|",
        "| optik map | SAM, primary + supplementary (SA), MAPQ | yes | yes (chaining + WFA/BiWFA) |",
        "| minimap2 | SAM/PAF, primary + supplementary | yes (`-x sr`) | yes (`map-ont`, `map-hifi`) |",
        "| wfmash | SAM/PAF, one primary per segment | not designed for | yes (WFA) |",
        "",
    ]
    by = {}
    for r in bench.rows:
        by.setdefault((r["section"], r["mode"]), []).append(r)
    for (sec, mode), rows in by.items():
        md.append(f"## {sec} -- {mode}")
        md.append("")
        cols = [
            k for k in rows[0] if k not in ("section", "mode", "gate_detail", "cmd")
        ]
        for r in rows[1:]:
            for k in r:
                if k not in cols and k not in ("section", "mode", "gate_detail", "cmd"):
                    cols.append(k)
        md.append("| " + " | ".join(cols) + " |")
        md.append("|" + "---|" * len(cols))
        for r in rows:
            cells = []
            for k in cols:
                v = r.get(k, "")
                if k == "gate":
                    v = "PASS" if v else "FAIL"
                elif isinstance(v, dict) and "mean" in v:
                    v = f"{v['mean']:.4g}" + (
                        f" ± {v['sd']:.2g}" if v.get("sd") is not None else ""
                    )
                elif isinstance(v, float):
                    v = f"{v:.4g}"
                cells.append(str(v))
            md.append("| " + " | ".join(cells) + " |")
        fails = [r for r in rows if not r.get("gate")]
        if fails:
            md.append("")
            md.append("Gate failures:")
            for r in fails:
                md.append(
                    f"- {r.get('tool') or r.get('method')} / {r.get('dataset') or r.get('length')}: "
                    f"{r.get('gate_detail')}"
                )
        md.append("")
    Path(os.path.join(out, "tables.md")).write_text("\n".join(md) + "\n")
    Path(os.path.join(out, "results.json")).write_text(
        json.dumps(bench.rows, indent=1, default=str)
    )
    hidden = ("section", "mode", "gate_detail", "cmd")
    tables = [
        (
            f"{sec} {mode}",
            [{k: v for k, v in r.items() if k not in hidden} for r in rows],
        )
        for (sec, mode), rows in by.items()
    ]
    fails = [
        f"- FAIL {r['section']} {r['mode']} {r.get('tool') or r.get('method')} / "
        f"{r.get('dataset') or r.get('length')}: {r.get('gate_detail')}"
        for r in bench.rows
        if not r.get("gate")
    ]
    results_file.write_section(
        dest,
        "small_sets",
        "optik against SSHash, SBWT, FMSI, Cuttlefish, Bifrost, minimap2 and wfmash on small sets"
        f" (wall s mean ± sd over {bench.runs} runs after {bench.warmup} warmup, peak RSS MB,"
        " gate PASS when the output is verified)",
        tables,
        notes=fails
        + [f"Raw rows, logs and the full tables: `{os.path.relpath(out, ROOT)}`"],
    )


def main():
    """Check the binaries, then run the sections and modes the command line selects."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mode", default="both", choices=["parity", "optimal", "both"])
    ap.add_argument("--sections", default="kmer,dbg,map,align,mapq")
    ap.add_argument(
        "--runs",
        type=int,
        default=3,
        help="measured runs per command (statistics over these)",
    )
    ap.add_argument(
        "--warmup", type=int, default=1, help="warmup runs per command, discarded"
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every parallel tool (default: the usable CPUs)",
    )
    ap.add_argument(
        "--cap-gb",
        type=int,
        default=hw.default_cap_bytes() >> 30,
        help="per-command resident-memory cap in GiB (memcap; default: usable memory less 2 GiB)",
    )
    ap.add_argument(
        "--out", default=os.path.join(HERE, "results", "compare_tools_small_sets")
    )
    ap.add_argument(
        "--work",
        default=os.path.join(T, "compare_tools_small_sets"),
        help="inputs and indexes",
    )
    ap.add_argument("--quick", action="store_true", help="small subset (smoke test)")
    results_file.add_argument(ap)
    a = ap.parse_args()
    # Tools run inside per-dataset directories, where relative paths would resolve wrongly.
    a.out, a.work = os.path.abspath(a.out), os.path.abspath(a.work)
    check_optik_fresh()
    stamp = datetime.datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")
    out = os.path.join(a.out, stamp + ("-quick" if a.quick else ""))
    os.makedirs(out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    for b in (
        "optik",
        "sshash",
        "sbwt",
        "fmsi",
        "kmercamel",
        "cuttlefish",
        "bifrost",
        "minimap2",
        "wfmash",
        "kmc",
        "kmc_dump",
    ):
        if not os.access(BIN[b], os.X_OK):
            sys.exit(f"missing {b}: {BIN[b]} (bench/build_tools.sh)")
    if not shutil.which("hyperfine"):
        sys.exit("hyperfine not on PATH")
    measure.rusage_exec()  # built up front, or rebuilt when its source is newer
    env = environment(out)
    bench = Bench(out, a.runs, a.warmup, a.cap_gb, os.path.join(a.work, "scratch"))
    modes = ["parity", "optimal"] if a.mode == "both" else [a.mode]
    secs = a.sections.split(",")
    for mode in modes:
        if "kmer" in secs:
            section_kmer(bench, mode, a.work, a.threads, a.quick)
        if "dbg" in secs:
            section_dbg(bench, mode, a.work, a.threads, a.quick)
        if "map" in secs:
            section_map(bench, mode, a.work, a.threads, a.quick)
        if "align" in secs:
            section_align(bench, mode, a.work, a.quick)
        report(bench, out, env, a.report)
    if "mapq" in secs:
        section_mapq(bench, "calibration", a.work, a.threads, a.quick)
    report(bench, out, env, a.report)
    nfail = sum(1 for r in bench.rows if not r.get("gate"))
    print(f"\n{len(bench.rows)} rows, {nfail} gate failures -> {out}/tables.md")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
