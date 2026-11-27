#!/usr/bin/env python3
"""scaling_bench.py - optik against the reference tools on growing genomes.

tool_bench.py's datasets are a few Mbp, where process start-up and I/O dominate and nothing scales.
This harness runs every category on the 4.6 Mbp E. coli genome, the 47 Mbp human chr21 and the 249 Mbp chr1.
Commands run alone under the resident-memory cap of bench/memcap.py, on any Linux, and is timed by bench/measure.py.
That means hyperfine discards the --warmup runs and reports the mean and standard deviation of the measured runs, each run's input files are cached first, and the peak RSS of each run comes from wait4 through rusage_exec, a small process, so the Python parent's resident pages are not folded in.

Inputs are the E. coli K-12 genome, which bench/genomes.py downloads into temp/genomes/ on first use, and hg38 chr21 and chr1, which `bench/build_tools.sh scaling` fetches into temp/scaling/.
Simulated reads come from bench/gen_reads.py with fixed seeds, and random 31-mers from a fixed RNG.

Per genome
  dbg        Cuttlefish (--ref), optik dbruijn -c 1 and Bifrost (-r)
  dict       SSHash build and query against optik build-spss and query, on Cuttlefish's unitigs (the common k-mer set), with reads, random k-mers (negative) and random genome k-mers on random strands (positive).
             Positive counts must agree with SSHash.
  fm         KmerCamel with FMSI against optik build-fm, on random negative and positive k-mers
  sbwt       SBWT plain-matrix build and random negative and positive queries
  map        minimap2 -ax sr against optik map on 200k reads of 150 bp (1% substitutions, 0.1% indels), scored against the simulation truth (correct means right strand and origin within 15 bp, with MAPQ >= 30 errors counted)

SSHash is built at m = k/2, since it has no default and no paper gives one for these inputs.
optik-spss runs untuned at its default m, which it picks per input from a sampled estimate of size plus lookup cost.
Every query tool runs on one thread, because SSHash, SBWT and FMSI have no query threads, and optik's queries are also timed at --threads as an optik-only row.
Query times cover whole processes, index loading included.

Output is one TSV row per (genome, tool, phase) in <workdir>/scaling.tsv, holding the mean wall time with sd, runs and warmup in its last column, plus <workdir>/scaling.meta (git revision, CPU governor, warmup and runs, command line).

Usage: python3 bench/scaling_bench.py <workdir> [--genomes ecoli,chr21,chr1] [--k 15,21,25,31] [--skip-map]
                                      [--threads <usable CPUs>] [--cap SIZE] [--warmup 1] [--runs 3]
"""

import argparse
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "bench"))
import math

import genomes
import hw
import measure
import memcap
import report as results_file
from config import PATHS
from sam_spec import (
    HIGH_MAPQ,
    SamCol,
    SamFlag,
    TruthCol,
)
from units import BITS_PER_BYTE, BYTES_PER_MIB

# measured runs of every command by default
DEFAULT_RUNS = 3
# Negative probes are KMER_SET_SIZE uniform random k-mers (seed NEGATIVE_SEED).
# Positive probes are as many genome k-mers (seed POSITIVE_SEED), each reverse-complemented with probability RC_PROB.
KMER_SET_SIZE = 1_000_000
NEGATIVE_SEED = 9
POSITIVE_SEED = 11
RC_PROB = 0.5
# a mapped read is correct when its origin lies within this many bases of the truth
ORIGIN_TOL = 15

GENOMES = {
    "ecoli": os.path.join(genomes.DEST, genomes.GENOMES["ecoli"][0]),
    "chr21": os.path.join(ROOT, "temp/scaling/chr21.fa"),
    "chr1": os.path.join(ROOT, "temp/scaling/chr1.fa"),
}


# An FMSI query result line holds the k-mer and its presence bit string.
# SSHash needs -m, and no paper gives one for these inputs, so m = k / 2.
K_PER_MINIMIZER_LEN = 2


class Runner:
    """Runs every command through bench/measure.py.

    Hyperfine ignores the warmups and computes statistics over the --runs measured runs.
    Each run's input files are read into the page cache first, each run's peak RSS is recorded, and the whole run stays under the resident-memory cap (memcap)."""

    def __init__(self, out, cap, threads, warmup=1, runs=DEFAULT_RUNS, scratch=None):
        self.out, self.cap, self.threads = out, cap, str(threads)
        self.warmup, self.runs = warmup, runs
        self.scratch = scratch or os.path.join(
            os.path.dirname(os.path.abspath(out)), ".measure"
        )
        self.last = None

    def run(self, cmd, stdout=None):
        """Return (mean wall s, peak MB, status, stdout + stderr text) of the last measured run.

        `stdout` names where that run's standard output is copied, in which case it is left out of the returned text."""
        ins = [arg for arg in map(str, cmd[1:]) if os.path.isfile(arg)]
        result = measure.measure(
            cmd,
            self.warmup,
            self.runs,
            memcap.parse_size(self.cap),
            ins,
            scratch=self.scratch,
            tag=os.path.basename(str(cmd[0])),
        )
        self.last = result
        text = result.err
        if os.path.exists(result.out):
            if stdout:
                shutil.copyfile(result.out, stdout)
            else:
                text = Path(result.out).read_text(errors="replace") + text
        if not result.succeeded:
            return None, None, "capped" if result.capped else f"rc{result.code}", text
        return result.wall["mean"], result.rss(), "ok", text

    run_once = run  # every command follows one policy, the same warmup and runs

    def row(self, genome, tool, phase, wall, peak_mb, status, extra=""):
        """Append one TSV row, adding the spread of the measurement it reports."""
        # The spread of the measurement this row reports, whose mean is `wall`
        result = self.last
        if result is not None and result.succeeded and wall is not None and wall == result.wall["mean"]:
            wall_sd = f"{result.wall['sd']:.3g}" if result.wall.get("sd") is not None else "-"
            extra = (f"sd={wall_sd} runs={result.runs} warmup={result.warmup} " + extra).rstrip()
        line = "\t".join(
            [
                genome,
                tool,
                phase,
                "" if wall is None else f"{wall:.2f}",
                "" if peak_mb is None else f"{peak_mb:.0f}",
                status,
                extra,
            ]
        )
        print(line, flush=True)
        with open(self.out, "a") as handle:
            handle.write(line + "\n")


def report_value(text, key):
    """Integer value of a 'key = N' line in text, or None."""
    match = re.search(rf"^{key}\s*=\s*(\d+)", text, re.MULTILINE)
    return int(match.group(1)) if match else None


def prep(work, name, src):
    """Write one genome's single-record FASTA and simulated reads, returning (dir, genome, reads)."""
    genome_dir = os.path.join(work, name)
    os.makedirs(genome_dir, exist_ok=True)
    genome_path = os.path.join(genome_dir, "genome.fa")
    if not os.path.exists(genome_path):
        # One record per genome, because joining several records would create k-mers that span their boundary and exist in no sequence.
        with open(src) as handle:
            n_hdr = sum(1 for line in handle if line.startswith(">"))
        if n_hdr != 1:
            raise SystemExit(
                f"{src}: {n_hdr} FASTA records; scaling_bench expects one chromosome per file"
            )
        with open(src) as in_handle, open(genome_path, "w") as out_handle:
            out_handle.write(">g\n")
            for line in in_handle:
                if not line.startswith(">"):
                    out_handle.write(line.strip().upper())
            out_handle.write("\n")
    reads = os.path.join(genome_dir, "reads.fa")
    if not os.path.exists(reads):
        subprocess.run(
            [
                sys.executable,
                os.path.join(ROOT, "bench/gen_reads.py"),
                "-r",
                genome_path,
                "-n",
                "200000",
                "-l",
                "150",
                "-e",
                "0.01",
                "--ins",
                "0.001",
                "--del",
                "0.001",
                "-s",
                "5",
                "-o",
                reads,
                "--gt",
                reads + ".tsv",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    return genome_dir, genome_path, reads


def kmer_sets(genome_dir, genome_path, kmer_len):
    """Write 1M uniform random k-mers for negative lookups, since they almost never occur in a genome, and 1M random genome k-mers on random strands for positive lookups."""
    rnd = os.path.join(genome_dir, f"rand_k{kmer_len}.fa")
    if not os.path.exists(rnd):
        rng = random.Random(NEGATIVE_SEED)
        with open(rnd, "w") as out_handle:
            out_handle.writelines(
                f">{i}\n{''.join(rng.choice('ACGT') for _ in range(kmer_len))}\n"
                for i in range(KMER_SET_SIZE)
            )
    rpos = os.path.join(genome_dir, f"rpos_k{kmer_len}.fa")
    if not os.path.exists(rpos):
        genome_seq = Path(genome_path).read_text().split("\n", 1)[1].replace("\n", "")
        rng = random.Random(POSITIVE_SEED)
        comp = str.maketrans("ACGT", "TGCA")
        with open(rpos, "w") as out_handle:
            count = 0
            while count < KMER_SET_SIZE:
                pos = rng.randrange(len(genome_seq) - kmer_len)
                kmer = genome_seq[pos : pos + kmer_len]
                if "N" in kmer:
                    continue
                if rng.random() < RC_PROB:
                    kmer = kmer.translate(comp)[::-1]
                out_handle.write(f">{count}\n{kmer}\n")
                count += 1
    return rnd, rpos


def fmsi_positives(path):
    """Count present k-mers in `fmsi query` output, one bitstring per record."""
    count = 0
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            # An FMSI result line is the query k-mer, a tab and its presence bit string.
            if len(parts) >= 2:
                count += parts[-1].count("1")
    return count


def sbwt_positives(path):
    """Count present k-mers in `sbwt search` output, one id per k-mer with -1 when absent."""
    count = 0
    with open(path) as handle:
        for line in handle:
            count += sum(1 for kmer_id in line.split() if kmer_id != "-1")
    return count


def map_accuracy(sam, truth_tsv):
    """Return (correct, total, MAPQ>=30 wrong, MAPQ>=30) from the first primary of each read.

    A read is correct when its strand matches and its origin (POS minus a leading clip) lies within 15.
    The MAPQ>=30 count is reported as the denominator of the wrong count, since a mapper that rates fewer reads confident trivially makes fewer confident errors."""
    truth = {}
    with open(truth_tsv) as handle:
        for line in handle:
            fields = line.split()
            if fields and fields[TruthCol.NAME] != "read_name":
                truth[fields[TruthCol.NAME]] = (int(fields[TruthCol.POS]), fields[TruthCol.STRAND])
    correct = hq_wrong = hq_total = 0
    done = set()
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            fields = line.split("\t", SamCol.RNEXT)
            flag = int(fields[SamCol.FLAG])
            if flag & SamFlag.NOT_PRIMARY or fields[SamCol.QNAME] in done:
                continue
            done.add(fields[SamCol.QNAME])
            clip_match = re.match(r"^(\d+)[SH]", fields[SamCol.CIGAR])
            origin = int(fields[SamCol.POS]) - 1 - (int(clip_match.group(1)) if clip_match else 0)
            truth_rec = truth.get(fields[SamCol.QNAME])
            if truth_rec:
                good = abs(origin - truth_rec[0]) <= ORIGIN_TOL and (
                    ("-" if flag & SamFlag.REVERSE else "+") == truth_rec[1]
                )
                correct += good
                hq_total += int(fields[SamCol.MAPQ]) >= HIGH_MAPQ
                hq_wrong += (not good) and int(fields[SamCol.MAPQ]) >= HIGH_MAPQ
    return correct, len(truth), hq_wrong, hq_total


def run_k(run, optik, bins, threads, genome, genome_dir, genome_path, reads, kmer_len, tag_k):
    """Run every k-dependent category at one k, in <genome dir>/k<k>."""
    k_str = str(kmer_len)
    name = f"{genome} k={kmer_len}" if tag_k else genome
    k_dir = os.path.join(genome_dir, f"k{kmer_len}")
    os.makedirs(k_dir, exist_ok=True)

    def in_dir(name):
        """Path of f in this k's directory."""
        return os.path.join(k_dir, name)

    rnd, rpos = kmer_sets(genome_dir, genome_path, kmer_len)
    num_kmers = None
    optik_pos = {}  # workload -> optik-spss positive count
    # De Bruijn graphs, whose Cuttlefish unitigs become the common k-mer set below
    wall, peak_mb, status, _ = run.run(
        [
            bins["cuttlefish"],
            "build",
            "--ref",
            "-s",
            genome_path,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("cf"),
            "-w",
            k_dir,
        ]
    )
    run.row(name, "cuttlefish", "dbg build (--ref)", wall, peak_mb, status)
    unitigs_fa = in_dir("cf.fa")
    # Unitig FASTA first, like for like with Cuttlefish, then the GFA with links
    wall, peak_mb, status, output = run.run(
        [
            optik,
            "dbruijn",
            "-i",
            genome_path,
            "-k",
            k_str,
            "-c",
            "1",
            "-t",
            threads,
            "-o",
            in_dir("optik.fa"),
        ]
    )
    run.row(
        name,
        "optik",
        "dbg build (-c 1, unitig FASTA)",
        wall,
        peak_mb,
        status,
        " ".join(re.findall(r"^unitigs\t\d+", output, re.MULTILINE)),
    )
    wall, peak_mb, status, output = run.run(
        [
            optik,
            "dbruijn",
            "-i",
            genome_path,
            "-k",
            k_str,
            "-c",
            "1",
            "-t",
            threads,
            "-o",
            in_dir("optik.gfa"),
        ]
    )
    run.row(
        name,
        "optik",
        "dbg build (-c 1, GFA + links)",
        wall,
        peak_mb,
        status,
        " ".join(re.findall(r"^unitigs\t\d+", output, re.MULTILINE)),
    )
    wall, peak_mb, status, _ = run.run(
        [
            bins["bifrost"],
            "build",
            "-r",
            genome_path,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("bf"),
        ]
    )
    run.row(name, "bifrost", "dbg build (-r)", wall, peak_mb, status)
    # Exact k-mer dictionaries.
    # SSHash requires -m and no paper gives one for these inputs, so m = k/2.
    minimizer_len = str(kmer_len // K_PER_MINIMIZER_LEN)
    wall, peak_mb, status, _ = run.run(
        [
            bins["sshash"],
            "build",
            "-i",
            unitigs_fa,
            "-k",
            k_str,
            "-m",
            minimizer_len,
            "-t",
            threads,
            "-d",
            k_dir,
            "-o",
            in_dir("s.idx"),
        ]
    )
    run.row(
        name,
        "sshash",
        f"build (m={minimizer_len})",
        wall,
        peak_mb,
        status,
        f"bytes={os.path.getsize(in_dir('s.idx')) if status == 'ok' else None}",
    )
    best = (minimizer_len, None) if status == "ok" else None
    wall, peak_mb, status, output = run.run(
        [
            optik,
            "build-spss",
            "-i",
            unitigs_fa,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("o.spss"),
        ]
    )
    # optik's own distinct k-mer count is the denominator for both tools' bits per k-mer.
    num_kmers = re.search(r"distinct_kmers\t(\d+)", output)
    num_kmers = int(num_kmers.group(1)) if num_kmers else None
    minimizer_match = re.search(r"minimizer_len\t(\d+)", output)
    run.row(
        name,
        "optik-spss",
        "build (default m)",
        wall,
        peak_mb,
        status,
        f"kmers={num_kmers} m={minimizer_match.group(1) if minimizer_match else '?'}",
    )
    for tool, path in (
        (f"sshash (m={best[0] if best else '?'})", "s.idx"),
        ("optik-spss", "o.spss"),
    ):
        if os.path.exists(in_dir(path)) and num_kmers:
            run.row(
                name,
                tool,
                "index bits/k-mer",
                None,
                None,
                "ok",
                f"{os.path.getsize(in_dir(path)) * BITS_PER_BYTE / num_kmers:.2f} ({os.path.getsize(in_dir(path)) / BYTES_PER_MIB:.0f} MB)",
            )
    for workload, query_path in (("reads", reads), ("random", rnd), ("random positive", rpos)):
        wall, peak_mb, status, output = run.run(
            [bins["sshash"], "query", "-i", in_dir("s.idx"), "-q", query_path]
        )
        ref_pos = report_value(output, "num_positive_kmers")
        run.row(
            name,
            "sshash",
            f"query {workload} (1 thread)",
            wall,
            peak_mb,
            status,
            f"kmers={report_value(output, 'num_kmers')} pos={ref_pos}",
        )
        # SSHash queries run on one thread, so optik runs at -t 1 for a like-for-like time and again at full threads.
        # Its positive count must equal SSHash's in both runs.
        for threads_arg in ("1", threads):
            wall, peak_mb, status, output = run.run(
                [optik, "query", "-i", in_dir("o.spss"), "-q", query_path, "-k", k_str, "-t", threads_arg]
            )
            pos = report_value(output, "num_positive_kmers")
            if threads_arg == "1":
                optik_pos[workload] = pos
            agree = "" if ref_pos is None or pos == ref_pos else " MISMATCH-vs-sshash"
            run.row(
                name,
                "optik-spss",
                f"query {workload} (-t {threads_arg})",
                wall,
                peak_mb,
                status,
                f"kmers={report_value(output, 'num_kmers')} pos={pos}{agree}",
            )
    # FM indexes.
    # -S is KmerCamel's fast path for unitig input, as in FMSI's own workflow.
    wall1, mb1, st1, _ = run.run(
        [
            bins["kmercamel"],
            "compute",
            "-k",
            k_str,
            "-S",
            "-M",
            in_dir("ms.fa"),
            "-o",
            "/dev/null",
            unitigs_fa,
        ]
    )
    sd1 = run.last.wall.get("sd") if st1 == "ok" else None
    wall2, mb2, st2, _ = (
        run.run([bins["fmsi"], "index", "-k", k_str, in_dir("ms.fa")])
        if st1 == "ok"
        else (None, None, "skipped", "")
    )
    # Two measured steps, so the row is the sum of their means and its sd combines both spreads.
    sd2 = run.last.wall.get("sd") if st2 == "ok" else None
    wall_sd = (
        # The two stages vary independently, so their deviations add in quadrature.
        f"{math.hypot(sd1, sd2):.3g}"
        if sd1 is not None and sd2 is not None
        else "-"
    )
    run.row(
        name,
        "fmsi",
        "build (kmercamel + index)",
        None if wall2 is None else wall1 + wall2,
        None if mb2 is None else max(mb1, mb2),
        st1 if st1 != "ok" else st2,
        (
            f"sd={wall_sd} runs={run.runs} warmup={run.warmup} kmercamel {wall1:.2f}s/{mb1:.0f}MB"
            if wall2 is not None
            else ""
        ),
    )
    wall, peak_mb, status, _ = run.run(
        [
            optik,
            "build-fm",
            "-i",
            unitigs_fa,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("o.fm"),
        ]
    )
    run.row(name, "optik-fm", "build", wall, peak_mb, status)
    # FMSI queries single-threaded, so optik runs at -t 1 for a fair comparison.
    for workload, query_path in (("random", rnd), ("random positive", rpos)):
        wall, peak_mb, status, output = run.run(
            [optik, "query", "-i", in_dir("o.fm"), "-q", query_path, "-k", k_str, "-t", "1"]
        )
        opos = report_value(output, "num_positive_kmers")
        run.row(name, "optik-fm", f"query {workload} (-t 1)", wall, peak_mb, status, f"pos={opos}")
        if st2 == "ok":
            wall, peak_mb, status, _ = run.run(
                [bins["fmsi"], "query", "-q", query_path, "-k", k_str, in_dir("ms.fa")],
                stdout=in_dir("fmsi.out"),
            )
            fpos = fmsi_positives(in_dir("fmsi.out")) if status == "ok" else None
            run.row(
                name,
                "fmsi",
                f"query {workload}",
                wall,
                peak_mb,
                status,
                # Both index the same unitigs, so their positive counts must agree.
                f"pos={fpos}" + ("" if fpos == opos else " MISMATCH-vs-optik"),
            )
    # SBWT, plain-matrix with reverse complements added
    wall, peak_mb, status, _ = run.run(
        [
            bins["sbwt"],
            "build",
            "-i",
            unitigs_fa,
            "-o",
            in_dir("x.sbwt"),
            "-k",
            k_str,
            "--variant",
            "plain-matrix",
            "--add-reverse-complements",
            "-t",
            threads,
            "--temp-dir",
            k_dir,
        ]
    )
    run.row(name, "sbwt-fast", "build", wall, peak_mb, status)
    if status == "ok":
        for workload, query_path in (("random", rnd), ("random positive", rpos)):
            wall, peak_mb, status, _ = run.run(
                [
                    bins["sbwt"],
                    "search",
                    "-i",
                    in_dir("x.sbwt"),
                    "-q",
                    query_path,
                    "-o",
                    in_dir("sbwt.out"),
                ]
            )
            spos = sbwt_positives(in_dir("sbwt.out")) if status == "ok" else None
            # Random k-mers do occur in a genome at small k (about 1% at k = 15), so the reference is optik's answer on the same file rather than 0.
            want = optik_pos.get(workload)
            run.row(
                name,
                "sbwt-fast",
                f"query {workload}",
                wall,
                peak_mb,
                status,
                f"pos={spos}" + ("" if spos == want else f" MISMATCH-vs-optik-{want}"),
            )
    # Index sizes per distinct k-mer, with every tool on the same k-mer set
    if num_kmers:
        sizes = {
            "optik-fm": [in_dir("o.fm")],
            "fmsi": [str(path) for path in Path(k_dir).glob("ms.fa.fmsi.*")],
            "sbwt-fast": [in_dir("x.sbwt")],
        }
        for tool, files in sizes.items():
            tot = sum(os.path.getsize(path) for path in files if os.path.exists(path))
            if tot:
                run.row(
                    name,
                    tool,
                    "index bits/k-mer",
                    None,
                    None,
                    "ok",
                    f"{tot * BITS_PER_BYTE / num_kmers:.2f} ({tot / BYTES_PER_MIB:.0f} MB)",
                )


def main():
    """Run every genome and k, then the mapping comparison."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("workdir")
    parser.add_argument("--genomes", default="ecoli,chr21,chr1")
    parser.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every parallel tool (default: the usable CPUs)",
    )
    parser.add_argument(
        "--cap",
        default=hw.default_cap_str(),
        help="per-command resident-memory cap (memcap; default: usable memory less 2 GiB)",
    )
    parser.add_argument(
        "--k",
        default="31",
        help="comma list of odd k (Cuttlefish needs odd k; SSHash and Bifrost k <= 31)",
    )
    parser.add_argument(
        "--skip-map", action="store_true", help="skip the (k-independent) mapping rows"
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
    os.makedirs(args.workdir, exist_ok=True)
    out = os.path.join(args.workdir, "scaling.tsv")
    # Rows of earlier runs in the same workdir stay in the TSV but not in this run's report.
    first_row = 0
    if os.path.exists(out):
        with open(out) as handle:
            first_row = sum(1 for _ in handle)
    run = Runner(out, args.cap, args.threads, args.warmup, args.runs)
    optik = str(PATHS["optik"])
    bins = {
        tool: str(PATHS[tool])
        for tool in (
            "sshash",
            "sbwt",
            "kmercamel",
            "cuttlefish",
            "bifrost",
            "minimap2",
            "fmsi",
        )
    }
    threads = run.threads
    with open(os.path.join(args.workdir, "scaling.meta"), "w") as handle:
        rev = subprocess.run(
            ["git", "-C", ROOT, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
        # A bare revision would claim results for code that was not committed.
        if subprocess.run(
            ["git", "-C", ROOT, "status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip():
            rev += "-dirty"
        gov = ""
        try:
            gov = (
                Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
                .read_text()
                .strip()
            )
        except OSError:
            pass  # no cpufreq here, so the governor is recorded empty
        handle.write(
            f"git\t{rev}\ngovernor\t{gov}\nthreads\t{threads}\ncap\t{args.cap}\nwarmup\t{args.warmup}\nruns\t{args.runs}\n"
            f"argv\t{' '.join(sys.argv)}\n"
        )

    for name in args.genomes.split(","):
        # E. coli is fetched on first use, and the human chromosomes come from build_tools.sh scaling.
        src = genomes.path("ecoli") if name == "ecoli" else GENOMES[name]
        if not os.path.exists(src):
            print(
                f"[skip] {name}: {src} missing (bench/build_tools.sh scaling)",
                file=sys.stderr,
            )
            continue
        genome_dir, genome_path, reads = prep(args.workdir, name, src)
        for kmer_len in [int(k_value) for k_value in args.k.split(",")]:
            run_k(
                run,
                optik,
                bins,
                threads,
                name,
                genome_dir,
                genome_path,
                reads,
                kmer_len,
                len(args.k.split(",")) > 1,
            )
        if args.skip_map:
            continue

        def in_dir(name, genome_dir=genome_dir):
            """Path of f in this genome's directory."""
            return os.path.join(genome_dir, name)

        # Mapping
        for tool, cmd in (
            ("minimap2", [bins["minimap2"], "-ax", "sr", "-t", threads, genome_path, reads]),
            # The same minimizers as minimap2 -ax sr (v2.31 options.c, k = 21, w = 11)
            (
                "optik",
                [
                    optik,
                    "map",
                    "-r",
                    genome_path,
                    "-q",
                    reads,
                    "-k",
                    "21",
                    "-w",
                    "11",
                    "-t",
                    threads,
                ],
            ),
        ):
            sam = in_dir(f"{tool}.sam")
            wall, peak_mb, status, _ = run.run(cmd, stdout=sam)
            correct, total, hqw, hq_total = (
                map_accuracy(sam, reads + ".tsv") if status == "ok" else (0, 0, 0, 0)
            )
            run.row(
                name,
                tool,
                f"map 200k reads (-t {threads}, index incl.)",
                wall,
                peak_mb,
                status,
                f"correct={correct}/{total} mapq30={hq_total} mapq30_wrong={hqw}",
            )
    report_rows(out, first_row, args.report, args.workdir)


# Columns of scaling.tsv, in order, as the results file names them
TSV_COLS = ("genome", "tool", "phase", "wall_s", "rss_MB", "status", "detail")


def report_rows(tsv, first_row, dest, workdir):
    """Writes the rows this run appended to scaling.tsv as the 'scaling' section, one table per genome."""
    rows_by_genome = {}
    with open(tsv) as handle:
        for line in list(handle)[first_row:]:
            row = dict(zip(TSV_COLS, line.rstrip("\n").split("\t"), strict=False))
            rows_by_genome.setdefault(row.pop("genome"), []).append(row)
    results_file.write_section(
        dest,
        "scaling",
        "Index, graph and mapping cost of every tool as genome size grows "
        "(wall s mean over the measured runs, peak RSS MB)",
        list(rows_by_genome.items()),
        notes=[f"Raw TSV and run metadata: `{os.path.relpath(workdir, ROOT)}`"],
    )


if __name__ == "__main__":
    main()
