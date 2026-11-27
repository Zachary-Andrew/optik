#!/usr/bin/env python3
"""dbg_bench.py - the measurements behind the README's de Bruijn graph section.

Each table there is regenerated here from fixed inputs and seeds.
Every command runs alone under a resident-memory cap (bench/memcap.py, any Linux, no systemd) and is timed by bench/measure.py, which means hyperfine discards the --warmup runs and reports the mean and standard deviation over the --runs measured runs, each run's input files are cached first, and every run records its peak RSS and disk traffic.
Every command runs on --threads, by default the usable CPUs.

Sections (--only)
  reads      simulated reads.
             Salmonella LT2 at 30x (150 bp, 0.5% substitutions, 0.05% insertions and 0.05% deletions) with -k 31 -c 2, and hg38 chr1 at 10x (150 bp, 0.5% substitutions) with -k 31 -c 2, comparing partitioned compaction with --graph table, without a limit and with --max-mem 4G, plus build-fm and build-spss --compact on the same reads
  modes      hg38 chr1 as GFA with links, partitioned compaction against --graph table at k = 31, 21 and 15 (-c 1), and the 10x reads at k = 21, -c 3
  grch38     GRCh38 primary (chr1-22, X, Y, M) at -k 31 -c 1, with optik GFA and FASTA, Cuttlefish --ref, and optik's compact graph (--graph table --max-mem 23G, which falls back to it when the table does not fit, run only under a cap of at least 24 GiB).
             Unitig sets are compared by bench/cpp/unitig_digest.
  index      build-fm and build-spss --compact on GRCh38 (k = 31) with bits per k-mer, and 20,000 random 150 bp windows of the genome queried exactly and with one substitution each.
             Both indexes must report every exact window k-mer present and agree on the substituted ones.
  gut        Cuttlefish 2's gut reads, 45 GB of FASTQ.gz in 18 files, at -k 27 -c 2, with optik at its default, optik --max-mem 13G and Cuttlefish 3 -m 13 --read -c 2.
             Reported are wall time, peak RSS, the bytes written and read from getrusage in bench/cpp/rusage_exec.c, and peak disk from statvfs sampled every second with output included.
             optik's two outputs must be byte-identical.

Scaling with genome size and the other categories (queries, FM, map) are in bench/scaling_bench.py.

Output goes to <out>/rows.tsv (section, input, tool, setting, wall s, peak MB, status, extra), <out>/meta.txt (git revision, SHA-256 of the optik binary, command line) and the per-command logs in <out>/logs/.

Usage: python3 bench/dbg_bench.py <out> [--only reads,modes,grch38,index,gut] [--cap SIZE]
                                       [--threads <usable CPUs>] [--warmup 1] [--runs 3]
"""

import argparse
import glob
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))

import hw
import memcap
import report as results_file
import scaling_bench
from config import PATHS
from units import BITS_PER_BYTE, BYTES_PER_GB, BYTES_PER_GIB

OPTIK = os.path.join(ROOT, "target/release/optik")
SAL = os.path.join(ROOT, "temp/salmonella_LT2.fa")
CHR1 = os.path.join(ROOT, "temp/scaling/chr1.fa")
HG38 = os.path.join(ROOT, "temp/scaling/hg38_primary.fa")
GUT = os.path.join(ROOT, "temp/published_sets/cf2-gut")
DIGEST_SRC = os.path.join(ROOT, "bench/cpp/unitig_digest.cpp")
SECTIONS = ["reads", "modes", "grch38", "index", "gut"]
# The compact-graph row of grch38 fixes the counting limit so the human node table cannot fit, and needs this much cap above it for the graph itself.
EF_MAX_MEM = "23G"
EF_MIN_CAP = 24 << 30
THREADS = str(hw.cpus())  # thread count of every command, set from --threads

# read coverage of the Salmonella and chr1 read sets
SAL_COVERAGE = 30
CHR1_COVERAGE = 10
# index section: windows of WINDOW_LEN bases (seeded by WINDOWS_SEED), each holding WINDOW_LEN - INDEX_K + 1 k-mers
WINDOWS_SEED = 7
NUM_WINDOWS = 20000
WINDOW_LEN = 150
INDEX_K = 31
DEFAULT_RUNS = 3


# Both k-mer indexes are built, each with the flags its build needs.
INDEX_BUILDS = (("build-fm", []), ("build-spss", ["--compact"]))


class Bench:
    """Runs and times commands through scaling_bench.Runner and writes rows.tsv."""

    def __init__(self, out, cap, work, warmup, runs):
        self.out = out
        self.cap = cap
        self.logs = os.path.join(out, "logs")
        os.makedirs(self.logs, exist_ok=True)
        self.rows = os.path.join(out, "rows.tsv")
        # Only the Runner's measure.py path is used.
        # Its row() is not, since rows go to rows.tsv.
        self.run = scaling_bench.Runner(
            os.devnull, cap, THREADS, warmup, runs, scratch=os.path.join(work, ".measure")
        )

    def row(self, section, inp, tool, setting, wall, peak_mb, status, extra=""):
        """Append one row to rows.tsv, adding the spread of the measurement it reports."""
        # The spread of the measurement this row reports, whose mean is `wall`
        last = self.run.last
        if last is not None and last.succeeded and wall is not None and wall == last.wall["mean"]:
            wall_sd = f"{last.wall['sd']:.3g}" if last.wall.get("sd") is not None else "-"
            extra = (f"sd={wall_sd} runs={last.runs} warmup={last.warmup} " + extra).rstrip()
        line = "\t".join(
            [
                section,
                inp,
                tool,
                setting,
                "" if wall is None else f"{wall:.2f}",
                "" if peak_mb is None else f"{peak_mb:.0f}",
                status,
                extra,
            ]
        )
        print(line, flush=True)
        with open(self.rows, "a") as handle:
            handle.write(line + "\n")

    def cmd(self, name, cmd):
        """Measure `cmd` through bench/measure.py, returning mean wall s, peak MB, status and output.

        The mean covers the --runs that follow the discarded warmups.
        The last run's output goes to logs/<name>.log."""
        wall, peak_mb, status, output = self.run.run(cmd)
        with open(os.path.join(self.logs, name + ".log"), "w") as handle:
            handle.write(" ".join(map(str, cmd)) + "\n" + (output or ""))
        return wall, peak_mb, status, output


def stat(out, key):
    """Value of a 'key<TAB>value' line in optik's report, or an empty string."""
    match = re.search(rf"^{key}\t(\S+)", out or "", re.MULTILINE)
    return match.group(1) if match else ""


def digest_tool(work):
    """Path of the unitig_digest binary, compiled on first use."""
    exe = os.path.join(work, "unitig_digest")
    if not os.path.exists(exe):
        subprocess.run(["g++", "-O2", "-o", exe, DIGEST_SRC], check=True)
    return exe


def digest(exe, path):
    """unitig_digest's order-independent digest of the unitig set in path."""
    return subprocess.run(
        [exe, path], capture_output=True, text=True, check=True
    ).stdout.strip()


def reads(work, ref, cov, name, extra_err):
    """Simulated reads from bench/gen_reads.py with seed 42, made once."""
    out = os.path.join(work, name)
    if not os.path.exists(out):
        subprocess.run(
            [
                sys.executable,
                os.path.join(ROOT, "bench/gen_reads.py"),
                "-r",
                ref,
                "-c",
                str(cov),
                "-l",
                "150",
                "-e",
                "0.005",
                *extra_err,
                "-s",
                "42",
                "-o",
                out + ".tmp",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        os.replace(out + ".tmp", out)
    return out


def sec_reads(bench, work):
    """The reads section.
    Returns the chr1 10x read file that the modes section reuses."""
    sal = reads(work, SAL, SAL_COVERAGE, "sal30x.fa", ["--ins", "0.0005", "--del", "0.0005"])
    wall, peak_mb, status, output = bench.cmd(
        "sal30x",
        [
            OPTIK,
            "dbruijn",
            "-i",
            sal,
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            THREADS,
            "-o",
            os.path.join(work, "sal30x.gfa"),
        ],
    )
    bench.row(
        "reads",
        "Salmonella 30x",
        "optik",
        f"-k 31 -c 2 -t {THREADS}, GFA",
        wall,
        peak_mb,
        status,
        f"unitigs={stat(output, 'unitigs')}",
    )
    # The 10x chr1 read set is built with both unitig builders, which must write byte-identical graphs.
    c10 = reads(work, CHR1, CHR1_COVERAGE, "chr1_10x.fa", [])
    for mode in ("part", "table"):
        for lim in ([], ["--max-mem", "4G"]):
            setting = f"--graph {mode}" + (" --max-mem 4G" if lim else "")
            wall, peak_mb, status, output = bench.cmd(
                f"chr1_10x_{mode}{'_4G' if lim else ''}",
                [
                    OPTIK,
                    "dbruijn",
                    "-i",
                    c10,
                    "-k",
                    "31",
                    "-c",
                    "2",
                    "-t",
                    THREADS,
                    "--graph",
                    mode,
                    *lim,
                    "-o",
                    os.path.join(work, f"c10_{mode}.gfa"),
                ],
            )
            bench.row(
                "reads",
                "chr1 10x",
                "optik dbruijn",
                setting + f", -k 31 -c 2 -t {THREADS}",
                wall,
                peak_mb,
                status,
                f"unitigs={stat(output, 'unitigs')}",
            )
    part_gfa, table_gfa = (os.path.join(work, f"c10_{mode_name}.gfa") for mode_name in ("part", "table"))
    same = subprocess.run(["cmp", "-s", part_gfa, table_gfa], check=False).returncode == 0
    bench.row(
        "reads",
        "chr1 10x",
        "optik dbruijn",
        "part vs table GFA",
        None,
        None,
        "same" if same else "DIFFERENT",
    )
    os.remove(part_gfa)
    os.remove(table_gfa)
    for sub, extra in INDEX_BUILDS:
        index_path = os.path.join(work, f"c10.{sub}")
        wall, peak_mb, status, _ = bench.cmd(
            f"chr1_10x_{sub}",
            [OPTIK, sub, *extra, "-i", c10, "-k", "31", "-c", "2", "-t", THREADS, "-o", index_path],
        )
        bench.row(
            "reads",
            "chr1 10x",
            f"optik {sub}",
            " ".join(extra + [f"-k 31 -c 2 -t {THREADS}"]),
            wall,
            peak_mb,
            status,
        )
        if os.path.exists(index_path):
            os.remove(index_path)
    return c10


def sec_modes(bench, work, c10):
    """The modes section, partitioned compaction against --graph table at several k and cutoffs."""
    # Each k runs both unitig builders on chr1, which must agree byte for byte, across the range of k users pick.
    for kmer_arg in ("31", "21", "15"):
        outs = []
        for mode in ("part", "table"):
            gfa_path = os.path.join(work, f"chr1_k{kmer_arg}_{mode}.gfa")
            wall, peak_mb, status, output = bench.cmd(
                f"chr1_k{kmer_arg}_{mode}",
                [
                    OPTIK,
                    "dbruijn",
                    "-i",
                    CHR1,
                    "-k",
                    kmer_arg,
                    "-c",
                    "1",
                    "-t",
                    THREADS,
                    "--graph",
                    mode,
                    "-o",
                    gfa_path,
                ],
            )
            bench.row(
                "modes",
                "chr1",
                "optik dbruijn",
                f"--graph {mode}, -k {kmer_arg} -c 1 -t {THREADS}, GFA",
                wall,
                peak_mb,
                status,
                f"unitigs={stat(output, 'unitigs')}",
            )
            outs.append(gfa_path)
        same = subprocess.run(["cmp", "-s", *outs], check=False).returncode == 0
        bench.row(
            "modes",
            "chr1",
            "optik dbruijn",
            f"part vs table GFA, k = {kmer_arg}",
            None,
            None,
            "same" if same else "DIFFERENT",
        )
        # The chr1 graphs are large, so each pair is deleted once compared.
        for gfa_path in outs:
            os.remove(gfa_path)
    # The 10x reads at -c 3 add the solid-k-mer filter, under which the two builders must still agree.
    outs = []
    for mode in ("part", "table"):
        gfa_path = os.path.join(work, f"c10_k21_{mode}.gfa")
        wall, peak_mb, status, output = bench.cmd(
            f"chr1_10x_k21_c3_{mode}",
            [
                OPTIK,
                "dbruijn",
                "-i",
                c10,
                "-k",
                "21",
                "-c",
                "3",
                "-t",
                THREADS,
                "--graph",
                mode,
                "-o",
                gfa_path,
            ],
        )
        bench.row(
            "modes",
            "chr1 10x",
            "optik dbruijn",
            f"--graph {mode}, -k 21 -c 3 -t {THREADS}, GFA",
            wall,
            peak_mb,
            status,
            f"unitigs={stat(output, 'unitigs')}",
        )
        outs.append(gfa_path)
    same = subprocess.run(["cmp", "-s", *outs], check=False).returncode == 0
    bench.row(
        "modes",
        "chr1 10x",
        "optik dbruijn",
        "part vs table GFA, k = 21 -c 3",
        None,
        None,
        "same" if same else "DIFFERENT",
    )


def sec_grch38(bench, work):
    """The grch38 section, optik and Cuttlefish on GRCh38 with the unitig sets compared by digest."""
    if not os.path.exists(HG38):
        print(
            f"[skip] grch38: {HG38} missing (bench/build_tools.sh hg38)",
            file=sys.stderr,
        )
        return
    digest_bin = digest_tool(work)
    runs = [
        (
            "optik",
            "partitioned, GFA with links",
            [
                OPTIK,
                "dbruijn",
                "-i",
                HG38,
                "-k",
                "31",
                "-c",
                "1",
                "-t",
                THREADS,
                "-o",
                os.path.join(work, "hg.gfa"),
            ],
            "hg.gfa",
        ),
        (
            "optik",
            "partitioned, unitig FASTA",
            [
                OPTIK,
                "dbruijn",
                "-i",
                HG38,
                "-k",
                "31",
                "-c",
                "1",
                "-t",
                THREADS,
                "-o",
                os.path.join(work, "hg.fa"),
            ],
            "hg.fa",
        ),
        (
            "cuttlefish",
            "--ref, unitig FASTA",
            [
                str(PATHS["cuttlefish"]),
                "build",
                "--ref",
                "-s",
                HG38,
                "-k",
                "31",
                "-t",
                THREADS,
                "-o",
                os.path.join(work, "hgcf"),
                "-w",
                work,
            ],
            "hgcf.fa",
        ),
        # The node table of 2.5e9 k-mers, about 40 GB, does not fit, so the counter builds its compact graph, 19 GB plus 20% working room, hence EF_MAX_MEM under EF_MIN_CAP.
        (
            "optik",
            f"--graph table --max-mem {EF_MAX_MEM} (compact graph)",
            [
                OPTIK,
                "dbruijn",
                "-i",
                HG38,
                "-k",
                "31",
                "-c",
                "1",
                "-t",
                THREADS,
                "--graph",
                "table",
                "--max-mem",
                EF_MAX_MEM,
                "-o",
                os.path.join(work, "hg_ef.fa"),
            ],
            "hg_ef.fa",
        ),
    ]
    # Every GRCh38 run is reduced to an order-independent digest of its unitigs, and all digests must equal the first one.
    ref = None
    for tool, setting, cmd, outf in runs:
        if outf == "hg_ef.fa" and memcap.parse_size(bench.cap) < EF_MIN_CAP:
            reason = f"needs a cap of at least {EF_MIN_CAP // BYTES_PER_GIB} GiB, got {bench.cap}"
            bench.row("grch38", "GRCh38", tool, setting, None, None, "skipped", reason)
            continue
        wall, peak_mb, status, output = bench.cmd(f"grch38_{outf}", cmd)
        out_path = os.path.join(work, outf)
        dig = digest(digest_bin, out_path) if status == "ok" and os.path.exists(out_path) else ""
        ref = ref or dig
        bench.row(
            "grch38",
            "GRCh38",
            tool,
            setting + f", -k 31 -t {THREADS}",
            wall,
            peak_mb,
            status,
            f"digest={dig}"
            + ("" if not dig or dig == ref else " DIFFERENT")
            + (f" links={stat(output, 'links')}" if outf.endswith(".gfa") else ""),
        )
        if os.path.exists(out_path):
            os.remove(out_path)


def windows(work):
    """Write 20,000 random N-free 150 bp windows of GRCh38 with seed 7, and the same windows with one substitution each at a random position."""
    exact_path, sb_ = os.path.join(work, "win_exact.fa"), os.path.join(work, "win_sub.fa")
    if os.path.exists(sb_):
        return exact_path, sb_
    seqs, cur = [], []
    with open(HG38) as handle:
        for line in handle:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    lens = [len(seq) for seq in seqs]
    tot = sum(lens)
    rng = random.Random(WINDOWS_SEED)
    with open(exact_path, "w") as exact_handle, open(sb_, "w") as sub_handle:
        count = 0
        while count < NUM_WINDOWS:
            offset = rng.randrange(tot)
            seq_idx = 0
            while offset >= lens[seq_idx]:
                offset -= lens[seq_idx]
                seq_idx += 1
            if offset + WINDOW_LEN > lens[seq_idx]:
                continue
            window = seqs[seq_idx][offset : offset + WINDOW_LEN]
            if set(window) - set("ACGT"):
                continue
            sub_pos = rng.randrange(WINDOW_LEN)
            mutated = window[:sub_pos] + rng.choice([base for base in "ACGT" if base != window[sub_pos]]) + window[sub_pos + 1 :]
            exact_handle.write(f">{count}\n{window}\n")
            sub_handle.write(f">{count}\n{mutated}\n")
            count += 1
    return exact_path, sb_


def sec_index(bench, work):
    """The index section, which builds both indexes on GRCh38 and queries exact and substituted windows."""
    idx, res = {}, []
    for sub, extra in INDEX_BUILDS:
        index_path = os.path.join(work, f"hg.{sub}")
        wall, peak_mb, status, output = bench.cmd(
            f"grch38_{sub}",
            [OPTIK, sub, *extra, "-i", HG38, "-k", "31", "-t", THREADS, "-o", index_path],
        )
        num_match = re.search(r"distinct_kmers\t(\d+)", output or "")
        res.append((sub, extra, wall, peak_mb, status, int(num_match.group(1)) if num_match else None))
        if status == "ok":
            idx[sub] = index_path
    # Bits per distinct k-mer, both over the count the builds report (build-spss prints it)
    # Both builders index the same k-mer set, so whichever reported a count sizes bits per k-mer for both.
    num_kmers = next((distinct for *_, distinct in res if distinct), None)
    for sub, extra, wall, peak_mb, status, distinct in res:
        bits = f"{os.path.getsize(idx[sub]) * BITS_PER_BYTE / num_kmers:.2f}" if sub in idx and num_kmers else ""
        bench.row(
            "index",
            "GRCh38",
            f"optik {sub}",
            " ".join(extra + [f"-k 31 -t {THREADS}"]),
            wall,
            peak_mb,
            status,
            f"distinct_kmers={distinct} bits_per_kmer={bits}",
        )
    exact_path, sub_path = windows(work)
    got = {}
    for sub, index_path in idx.items():
        for workload, query_path in (("exact", exact_path), ("one substitution", sub_path)):
            wall, peak_mb, status, output = bench.cmd(
                f"grch38_{sub}_query_{workload.split()[0]}",
                [OPTIK, "query", "-i", index_path, "-q", query_path, "-k", "31", "-t", THREADS],
            )
            pos, tot = scaling_bench.report_value(output, "num_positive_kmers"), scaling_bench.report_value(output, "num_kmers")
            got[(sub, workload)] = pos
            # Exact windows come from the genome and must hit fully. Substituted ones have no fixed answer.
            passed = status == "ok" and tot == NUM_WINDOWS * (WINDOW_LEN - INDEX_K + 1) and (workload != "exact" or pos == tot)
            bench.row(
                "index",
                "GRCh38 windows",
                f"optik {sub}",
                f"query {workload}",
                wall,
                peak_mb,
                status if passed else "FAIL",
                f"kmers={tot} present={pos}",
            )
    # Substituted windows have no oracle, but both index types must still give the same answer.
    if len(idx) == len(INDEX_BUILDS):
        agree = (
            got[("build-fm", "one substitution")]
            == got[("build-spss", "one substitution")]
        )
        bench.row(
            "index",
            "GRCh38 windows",
            "optik",
            "build-fm vs build-spss, substituted windows",
            None,
            None,
            "same" if agree else "DIFFERENT",
        )
    for index_path in idx.values():
        os.remove(index_path)


class DiskPeak(threading.Thread):
    """Filesystem growth under `path`, sampled every second through statvfs."""

    def __init__(self, path):
        super().__init__(daemon=True)
        self.path, self.peak, self.stop = path, 0, threading.Event()
        self.base = self.used()

    def used(self):
        """Bytes in use on the filesystem."""
        fs_stat = os.statvfs(self.path)
        return (fs_stat.f_blocks - fs_stat.f_bfree) * fs_stat.f_frsize

    def run(self):
        while not self.stop.wait(1.0):
            self.peak = max(self.peak, self.used() - self.base)


def sec_gut(bench, work):
    """The gut section, optik and Cuttlefish 3 on Cuttlefish 2's gut read set."""
    # The files the fetch recorded, since PROVENANCE.json lists nine paired individuals in 18 files
    want = (
        sorted(
            entry["path"]
            for entry in json.loads(Path(os.path.join(GUT, "PROVENANCE.json")).read_text())[
                "files"
            ]
        )
        if os.path.exists(os.path.join(GUT, "PROVENANCE.json"))
        else []
    )
    ins = [path for path in want if os.path.exists(path)]
    if not want or len(ins) != len(want):
        bench.row(
            "gut",
            "cf2-gut",
            "-",
            "-",
            None,
            None,
            "skipped",
            f"{len(ins)} of {len(want)} files in {GUT}",
        )
        return
    # On the gut read set, optik runs at its default and again with --max-mem 13G, and Cuttlefish 3 at its documented setting.
    runs = [
        (
            "optik",
            "default",
            [
                OPTIK,
                "dbruijn",
                *[arg for path in ins for arg in ("-i", path)],
                "-k",
                "27",
                "-c",
                "2",
                "-t",
                THREADS,
                "-o",
                os.path.join(work, "gut.fa"),
            ],
        ),
        (
            "optik",
            "--max-mem 13G",
            [
                OPTIK,
                "dbruijn",
                *[arg for path in ins for arg in ("-i", path)],
                "-k",
                "27",
                "-c",
                "2",
                "-t",
                THREADS,
                "--max-mem",
                "13G",
                "-o",
                os.path.join(work, "gut13.fa"),
            ],
        ),
        (
            "cuttlefish",
            "-m 13",
            [
                str(PATHS["cuttlefish"]),
                "build",
                *[arg for path in ins for arg in ("-s", path)],
                "-k",
                "27",
                "-t",
                THREADS,
                "-m",
                "13",
                "--read",
                "-c",
                "2",
                "-o",
                os.path.join(work, "gutcf"),
                "-w",
                work,
            ],
        ),
    ]
    for tool, setting, cmd in runs:
        # The disk peak spans every run, warmups included, with each run's output replacing the last.
        disk_peak = DiskPeak(work)
        disk_peak.start()
        wall, peak_mb, status, output = bench.cmd(f"gut_{tool}_{setting.split()[0].strip('-')}", cmd)
        disk_peak.stop.set()
        disk_peak.join()
        last = bench.run.last
        io_stats = (
            f"written_gb={last.written / BYTES_PER_GB:.1f} read_gb={last.read / BYTES_PER_GB:.1f}"
            if last is not None and last.succeeded and last.written is not None
            else ""
        )
        bench.row(
            "gut",
            "cf2-gut",
            tool,
            setting + f", -k 27 -c 2 -t {THREADS}",
            wall,
            peak_mb,
            status,
            f"{io_stats} peak_disk_gb={disk_peak.peak / BYTES_PER_GB:.0f} unitigs={stat(output, 'unitigs')}",
        )
    gut_fa, gut13_fa = os.path.join(work, "gut.fa"), os.path.join(work, "gut13.fa")
    if os.path.exists(gut_fa) and os.path.exists(gut13_fa):
        same = subprocess.run(["cmp", "-s", gut_fa, gut13_fa], check=False).returncode == 0
        bench.row(
            "gut",
            "cf2-gut",
            "optik",
            "default vs --max-mem 13G output",
            None,
            None,
            "same" if same else "DIFFERENT",
        )
    for path in glob.glob(os.path.join(work, "gut*.fa")):
        os.remove(path)


def main():
    """Entry point that runs the sections named by --only."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("out")
    parser.add_argument("--only", default=",".join(SECTIONS))
    parser.add_argument(
        "--cap",
        default=hw.default_cap_str(),
        help="per-command resident-memory cap (memcap; default: usable memory less 2 GiB)",
    )
    parser.add_argument("--work", default=os.path.join(ROOT, "temp/dbg_bench_work"))
    parser.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every command (default: the usable CPUs)",
    )
    parser.add_argument(
        "--warmup", type=int, default=1, help="warmup runs of every command, discarded"
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=DEFAULT_RUNS,
        help="measured runs of every command (statistics over these)",
    )
    results_file.add_argument(parser)
    args = parser.parse_args()
    global THREADS
    THREADS = str(args.threads)
    only = args.only.split(",")
    bad = [section for section in only if section not in SECTIONS]
    if bad:
        sys.exit(f"unknown section(s) {bad}; known: {SECTIONS}")
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.work, exist_ok=True)
    # meta.txt records the git revision and the binary's hash, so a result can be traced to the build that made it.
    rev = subprocess.run(
        ["git", "-C", ROOT, "describe", "--always", "--dirty"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    sha = hashlib.sha256(Path(OPTIK).read_bytes()).hexdigest()
    with open(os.path.join(args.out, "meta.txt"), "w") as handle:
        handle.write(
            f"git {rev}\noptik sha256 {sha}\ncap {args.cap}\nthreads {THREADS}\nwarmup {args.warmup}\nruns {args.runs}\n"
            f"started {time.ctime()}\n"
            f"cmd {' '.join(sys.argv)}\n"
        )
    bench = Bench(args.out, args.cap, args.work, args.warmup, args.runs)
    # Rows of earlier runs into the same directory stay in rows.tsv but not in this run's report.
    first_row = 0
    if os.path.exists(bench.rows):
        with open(bench.rows) as handle:
            first_row = sum(1 for _ in handle)
    c10 = None
    if "reads" in only or "modes" in only:
        c10 = (
            sec_reads(bench, args.work)
            if "reads" in only
            else reads(args.work, CHR1, CHR1_COVERAGE, "chr1_10x.fa", [])
        )
    if "modes" in only:
        sec_modes(bench, args.work, c10)
    if "grch38" in only:
        sec_grch38(bench, args.work)
    if "index" in only:
        sec_index(bench, args.work)
    if "gut" in only:
        sec_gut(bench, args.work)
    with open(os.path.join(args.out, "meta.txt"), "a") as handle:
        handle.write(f"finished {time.ctime()}\n")
    rows_by_section = {}
    with open(bench.rows) as handle:
        for line in list(handle)[first_row:]:
            row = dict(zip(TSV_COLS, line.rstrip("\n").split("\t"), strict=False))
            rows_by_section.setdefault(row.pop("section"), []).append(row)
    results_file.write_section(
        args.report,
        "dbg",
        "optik dbruijn against Cuttlefish, GGCAT and Bifrost on read sets, modes, GRCh38, "
        "index building and a gut metagenome (wall s mean over the measured runs, peak RSS MB)",
        list(rows_by_section.items()),
        notes=[
            f"Raw TSV, metadata and per-command logs: `{os.path.relpath(args.out, ROOT)}`"
        ],
    )


# Columns of rows.tsv, in order, as the results file names them
TSV_COLS = (
    "section",
    "input",
    "tool",
    "setting",
    "wall_s",
    "rss_MB",
    "status",
    "detail",
)


if __name__ == "__main__":
    main()
