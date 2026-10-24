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

Output is one TSV row per (genome, tool, phase) in <workdir>/scaling.tsv, holding the mean wall time with sd, runs and warmup in its last column, plus
<workdir>/scaling.meta (git revision, CPU governor, warmup and runs, command
line).

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

R = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(R, "bench"))
import genomes
import hw
import measure
import memcap
import report as results_file
from config import PATHS
from sam_spec import HIGH_MAPQ, SAM_FLAG_REVERSE, SAM_FLAGS_NOT_PRIMARY

GENOMES = {
    "ecoli": os.path.join(genomes.DEST, genomes.GENOMES["ecoli"][0]),
    "chr21": os.path.join(R, "temp/scaling/chr21.fa"),
    "chr1": os.path.join(R, "temp/scaling/chr1.fa"),
}


class Runner:
    """Runs every command through bench/measure.py.

    Hyperfine ignores the warmups and computes statistics over the --runs measured runs.
    Each run's input files are read into the page cache first, each run's peak RSS is recorded, and the whole run stays under the resident-memory cap (memcap)."""

    def __init__(self, out, cap, threads, warmup=1, runs=3, scratch=None):
        self.out, self.cap, self.threads = out, cap, str(threads)
        self.warmup, self.runs = warmup, runs
        self.scratch = scratch or os.path.join(
            os.path.dirname(os.path.abspath(out)), ".measure"
        )
        self.last = None

    def run(self, cmd, stdout=None):
        """Return (mean wall s, peak MB, status, stdout + stderr text) of the last measured run.

        `stdout` names where that run's standard output is copied, in which case it is left out of the returned text."""
        ins = [c for c in map(str, cmd[1:]) if os.path.isfile(c)]
        r = measure.measure(
            cmd,
            self.warmup,
            self.runs,
            memcap.parse_size(self.cap),
            ins,
            scratch=self.scratch,
            tag=os.path.basename(str(cmd[0])),
        )
        self.last = r
        text = r.err
        if os.path.exists(r.out):
            if stdout:
                shutil.copyfile(r.out, stdout)
            else:
                text = Path(r.out).read_text(errors="replace") + text
        if not r.ok:
            return None, None, "capped" if r.capped else f"rc{r.code}", text
        return r.wall["mean"], r.rss(), "ok", text

    run_once = run  # every command follows one policy, the same warmup and runs

    def row(self, genome, tool, phase, wall, mb, status, extra=""):
        """Append one TSV row, adding the spread of the measurement it reports."""
        # The spread of the measurement this row reports, whose mean is `wall`
        r = self.last
        if r is not None and r.ok and wall is not None and wall == r.wall["mean"]:
            sd = f"{r.wall['sd']:.3g}" if r.wall.get("sd") is not None else "-"
            extra = (f"sd={sd} runs={r.runs} warmup={r.warmup} " + extra).rstrip()
        line = "\t".join(
            [
                genome,
                tool,
                phase,
                "" if wall is None else f"{wall:.2f}",
                "" if mb is None else f"{mb:.0f}",
                status,
                extra,
            ]
        )
        print(line, flush=True)
        with open(self.out, "a") as f:
            f.write(line + "\n")


def kv(text, key):
    """Integer value of a 'key = N' line in text, or None."""
    m = re.search(rf"^{key}\s*=\s*(\d+)", text, re.MULTILINE)
    return int(m.group(1)) if m else None


def prep(work, name, src):
    """Write one genome's single-record FASTA and simulated reads, returning (dir, genome, reads)."""
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    up = os.path.join(d, "genome.fa")
    if not os.path.exists(up):
        # One record per genome, because joining several records would create k-mers that span their boundary and exist in no sequence.
        with open(src) as fh:
            n_hdr = sum(1 for line in fh if line.startswith(">"))
        if n_hdr != 1:
            raise SystemExit(
                f"{src}: {n_hdr} FASTA records; scaling_bench expects one chromosome per file"
            )
        with open(src) as f, open(up, "w") as o:
            o.write(">g\n")
            for line in f:
                if not line.startswith(">"):
                    o.write(line.strip().upper())
            o.write("\n")
    reads = os.path.join(d, "reads.fa")
    if not os.path.exists(reads):
        subprocess.run(
            [
                sys.executable,
                os.path.join(R, "bench/gen_reads.py"),
                "-r",
                up,
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
    return d, up, reads


def kmer_sets(d, up, k):
    """Write 1M uniform random k-mers for negative lookups, since they almost never occur in a genome, and 1M random genome k-mers on random strands for positive lookups."""
    rnd = os.path.join(d, f"rand_k{k}.fa")
    if not os.path.exists(rnd):
        r = random.Random(9)
        with open(rnd, "w") as o:
            o.writelines(
                f">{i}\n{''.join(r.choice('ACGT') for _ in range(k))}\n"
                for i in range(1_000_000)
            )
    rpos = os.path.join(d, f"rpos_k{k}.fa")
    if not os.path.exists(rpos):
        g = Path(up).read_text().split("\n", 1)[1].replace("\n", "")
        r = random.Random(11)
        comp = str.maketrans("ACGT", "TGCA")
        with open(rpos, "w") as o:
            n = 0
            while n < 1_000_000:
                i = r.randrange(len(g) - k)
                w = g[i : i + k]
                if "N" in w:
                    continue
                if r.random() < 0.5:
                    w = w.translate(comp)[::-1]
                o.write(f">{n}\n{w}\n")
                n += 1
    return rnd, rpos


def fmsi_positives(path):
    """Count present k-mers in `fmsi query` output, one bitstring per record."""
    n = 0
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                n += parts[-1].count("1")
    return n


def sbwt_positives(path):
    """Count present k-mers in `sbwt search` output, one id per k-mer with -1 when absent."""
    n = 0
    with open(path) as fh:
        for line in fh:
            n += sum(1 for x in line.split() if x != "-1")
    return n


def map_accuracy(sam, truth_tsv):
    """Return (correct, total, MAPQ>=30 wrong, MAPQ>=30) from the first primary of each read.

    A read is correct when its strand matches and its origin (POS minus a leading clip) lies within 15.
    The MAPQ>=30 count is reported as the denominator of the wrong count, since a mapper that rates fewer reads confident trivially makes fewer confident errors."""
    truth = {}
    with open(truth_tsv) as fh:
        for line in fh:
            t = line.split()
            if t and t[0] != "read_name":
                truth[t[0]] = (int(t[2]), t[3])
    ok = hq_wrong = hq = 0
    done = set()
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            t = line.split("\t", 6)
            f = int(t[1])
            if f & SAM_FLAGS_NOT_PRIMARY or t[0] in done:
                continue
            done.add(t[0])
            m = re.match(r"^(\d+)[SH]", t[5])
            origin = int(t[3]) - 1 - (int(m.group(1)) if m else 0)
            g = truth.get(t[0])
            if g:
                good = abs(origin - g[0]) <= 15 and (
                    ("-" if f & SAM_FLAG_REVERSE else "+") == g[1]
                )
                ok += good
                hq += int(t[4]) >= HIGH_MAPQ
                hq_wrong += (not good) and int(t[4]) >= HIGH_MAPQ
    return ok, len(truth), hq_wrong, hq


def run_k(run, optik, bins, threads, genome, gd, up, reads, k, tag_k):
    """Run every k-dependent category at one k, in <genome dir>/k<k>."""
    k_str = str(k)
    name = f"{genome} k={k}" if tag_k else genome
    d = os.path.join(gd, f"k{k}")
    os.makedirs(d, exist_ok=True)

    def in_dir(f):
        """Path of f in this k's directory."""
        return os.path.join(d, f)

    rnd, rpos = kmer_sets(gd, up, k)
    nk = None
    optik_pos = {}  # workload -> optik-spss positive count
    # De Bruijn graphs, whose Cuttlefish unitigs become the common k-mer set below
    w, mb, st, _ = run.run(
        [
            bins["cuttlefish"],
            "build",
            "--ref",
            "-s",
            up,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("cf"),
            "-w",
            d,
        ]
    )
    run.row(name, "cuttlefish", "dbg build (--ref)", w, mb, st)
    unitigs_fa = in_dir("cf.fa")
    # Unitig FASTA first, like for like with Cuttlefish, then the GFA with links
    w, mb, st, o = run.run(
        [
            optik,
            "dbruijn",
            "-i",
            up,
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
        w,
        mb,
        st,
        " ".join(re.findall(r"^unitigs\t\d+", o, re.MULTILINE)),
    )
    w, mb, st, o = run.run(
        [
            optik,
            "dbruijn",
            "-i",
            up,
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
        w,
        mb,
        st,
        " ".join(re.findall(r"^unitigs\t\d+", o, re.MULTILINE)),
    )
    w, mb, st, _ = run.run(
        [
            bins["bifrost"],
            "build",
            "-r",
            up,
            "-k",
            k_str,
            "-t",
            threads,
            "-o",
            in_dir("bf"),
        ]
    )
    run.row(name, "bifrost", "dbg build (-r)", w, mb, st)
    # Exact k-mer dictionaries.
    # SSHash requires -m and no paper gives one for these inputs, so m = k/2.
    m = str(k // 2)
    w, mb, st, _ = run.run(
        [
            bins["sshash"],
            "build",
            "-i",
            unitigs_fa,
            "-k",
            k_str,
            "-m",
            m,
            "-t",
            threads,
            "-d",
            d,
            "-o",
            in_dir("s.idx"),
        ]
    )
    run.row(
        name,
        "sshash",
        f"build (m={m})",
        w,
        mb,
        st,
        f"bytes={os.path.getsize(in_dir('s.idx')) if st == 'ok' else None}",
    )
    best = (m, None) if st == "ok" else None
    w, mb, st, o = run.run(
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
    nk = re.search(r"distinct_kmers\t(\d+)", o)
    nk = int(nk.group(1)) if nk else None
    om = re.search(r"minimizer_len\t(\d+)", o)
    run.row(
        name,
        "optik-spss",
        "build (default m)",
        w,
        mb,
        st,
        f"kmers={nk} m={om.group(1) if om else '?'}",
    )
    for tool, path in (
        (f"sshash (m={best[0] if best else '?'})", "s.idx"),
        ("optik-spss", "o.spss"),
    ):
        if os.path.exists(in_dir(path)) and nk:
            run.row(
                name,
                tool,
                "index bits/k-mer",
                None,
                None,
                "ok",
                f"{os.path.getsize(in_dir(path)) * 8 / nk:.2f} ({os.path.getsize(in_dir(path)) / 1048576:.0f} MB)",
            )
    for wl, q in (("reads", reads), ("random", rnd), ("random positive", rpos)):
        w, mb, st, o = run.run(
            [bins["sshash"], "query", "-i", in_dir("s.idx"), "-q", q]
        )
        ref_pos = kv(o, "num_positive_kmers")
        run.row(
            name,
            "sshash",
            f"query {wl} (1 thread)",
            w,
            mb,
            st,
            f"kmers={kv(o, 'num_kmers')} pos={ref_pos}",
        )
        for t in ("1", threads):
            w, mb, st, o = run.run(
                [optik, "query", "-i", in_dir("o.spss"), "-q", q, "-k", k_str, "-t", t]
            )
            pos = kv(o, "num_positive_kmers")
            if t == "1":
                optik_pos[wl] = pos
            agree = "" if ref_pos is None or pos == ref_pos else " MISMATCH-vs-sshash"
            run.row(
                name,
                "optik-spss",
                f"query {wl} (-t {t})",
                w,
                mb,
                st,
                f"kmers={kv(o, 'num_kmers')} pos={pos}{agree}",
            )
    # FM indexes.
    # -S is KmerCamel's fast path for unitig input, as in FMSI's own workflow.
    w1, mb1, st1, _ = run.run(
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
    w2, mb2, st2, _ = (
        run.run([bins["fmsi"], "index", "-k", k_str, in_dir("ms.fa")])
        if st1 == "ok"
        else (None, None, "skipped", "")
    )
    # Two measured steps, so the row is the sum of their means and its sd combines both spreads.
    sd2 = run.last.wall.get("sd") if st2 == "ok" else None
    sd = (
        f"{(sd1**2 + sd2**2) ** 0.5:.3g}"
        if sd1 is not None and sd2 is not None
        else "-"
    )
    run.row(
        name,
        "fmsi",
        "build (kmercamel + index)",
        None if w2 is None else w1 + w2,
        None if mb2 is None else max(mb1, mb2),
        st1 if st1 != "ok" else st2,
        (
            f"sd={sd} runs={run.runs} warmup={run.warmup} kmercamel {w1:.2f}s/{mb1:.0f}MB"
            if w2 is not None
            else ""
        ),
    )
    w, mb, st, _ = run.run(
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
    run.row(name, "optik-fm", "build", w, mb, st)
    for wl, q in (("random", rnd), ("random positive", rpos)):
        w, mb, st, o = run.run(
            [optik, "query", "-i", in_dir("o.fm"), "-q", q, "-k", k_str, "-t", "1"]
        )
        opos = kv(o, "num_positive_kmers")
        run.row(name, "optik-fm", f"query {wl} (-t 1)", w, mb, st, f"pos={opos}")
        if st2 == "ok":
            w, mb, st, _ = run.run(
                [bins["fmsi"], "query", "-q", q, "-k", k_str, in_dir("ms.fa")],
                stdout=in_dir("fmsi.out"),
            )
            fpos = fmsi_positives(in_dir("fmsi.out")) if st == "ok" else None
            run.row(
                name,
                "fmsi",
                f"query {wl}",
                w,
                mb,
                st,
                f"pos={fpos}" + ("" if fpos == opos else " MISMATCH-vs-optik"),
            )
    # SBWT, plain-matrix with reverse complements added
    w, mb, st, _ = run.run(
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
            d,
        ]
    )
    run.row(name, "sbwt-fast", "build", w, mb, st)
    if st == "ok":
        for wl, q in (("random", rnd), ("random positive", rpos)):
            w, mb, st, _ = run.run(
                [
                    bins["sbwt"],
                    "search",
                    "-i",
                    in_dir("x.sbwt"),
                    "-q",
                    q,
                    "-o",
                    in_dir("sbwt.out"),
                ]
            )
            spos = sbwt_positives(in_dir("sbwt.out")) if st == "ok" else None
            # Random k-mers do occur in a genome at small k (about 1% at k = 15), so the reference is optik's answer on the same file rather than 0.
            want = optik_pos.get(wl)
            run.row(
                name,
                "sbwt-fast",
                f"query {wl}",
                w,
                mb,
                st,
                f"pos={spos}" + ("" if spos == want else f" MISMATCH-vs-optik-{want}"),
            )
    # Index sizes per distinct k-mer, with every tool on the same k-mer set
    if nk:
        sizes = {
            "optik-fm": [in_dir("o.fm")],
            "fmsi": [str(x) for x in Path(d).glob("ms.fa.fmsi.*")],
            "sbwt-fast": [in_dir("x.sbwt")],
        }
        for tool, files in sizes.items():
            tot = sum(os.path.getsize(f) for f in files if os.path.exists(f))
            if tot:
                run.row(
                    name,
                    tool,
                    "index bits/k-mer",
                    None,
                    None,
                    "ok",
                    f"{tot * 8 / nk:.2f} ({tot / 1048576:.0f} MB)",
                )


def main():
    """Run every genome and k, then the mapping comparison."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("workdir")
    ap.add_argument("--genomes", default="ecoli,chr21,chr1")
    ap.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every parallel tool (default: the usable CPUs)",
    )
    ap.add_argument(
        "--cap",
        default=hw.default_cap_str(),
        help="per-command resident-memory cap (memcap; default: usable memory less 2 GiB)",
    )
    ap.add_argument(
        "--k",
        default="31",
        help="comma list of odd k (Cuttlefish needs odd k; SSHash and Bifrost k <= 31)",
    )
    ap.add_argument(
        "--skip-map", action="store_true", help="skip the (k-independent) mapping rows"
    )
    ap.add_argument(
        "--warmup", type=int, default=1, help="warmup runs of every command, discarded"
    )
    ap.add_argument(
        "--runs",
        type=int,
        default=3,
        help="measured runs of every command (statistics over these)",
    )
    results_file.add_argument(ap)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    out = os.path.join(a.workdir, "scaling.tsv")
    # Rows of earlier runs in the same workdir stay in the TSV but not in this run's report.
    first_row = 0
    if os.path.exists(out):
        with open(out) as fh:
            first_row = sum(1 for _ in fh)
    run = Runner(out, a.cap, a.threads, a.warmup, a.runs)
    optik = str(PATHS["optik"])
    bins = {
        t: str(PATHS[t])
        for t in (
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
    with open(os.path.join(a.workdir, "scaling.meta"), "w") as f:
        rev = subprocess.run(
            ["git", "-C", R, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
        # A bare revision would claim results for code that was not committed.
        if subprocess.run(
            ["git", "-C", R, "status", "--porcelain", "--untracked-files=no"],
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
            pass
        f.write(
            f"git\t{rev}\ngovernor\t{gov}\nthreads\t{threads}\ncap\t{a.cap}\nwarmup\t{a.warmup}\nruns\t{a.runs}\n"
            f"argv\t{' '.join(sys.argv)}\n"
        )

    for name in a.genomes.split(","):
        # E. coli is fetched on first use, and the human chromosomes come from build_tools.sh scaling.
        src = genomes.path("ecoli") if name == "ecoli" else GENOMES[name]
        if not os.path.exists(src):
            print(
                f"[skip] {name}: {src} missing (bench/build_tools.sh scaling)",
                file=sys.stderr,
            )
            continue
        d, up, reads = prep(a.workdir, name, src)
        for k in [int(x) for x in a.k.split(",")]:
            run_k(
                run,
                optik,
                bins,
                threads,
                name,
                d,
                up,
                reads,
                k,
                len(a.k.split(",")) > 1,
            )
        if a.skip_map:
            continue

        def in_dir(f, d=d):
            """Path of f in this genome's directory."""
            return os.path.join(d, f)

        # Mapping
        for tool, cmd in (
            ("minimap2", [bins["minimap2"], "-ax", "sr", "-t", threads, up, reads]),
            # The same minimizers as minimap2 -ax sr (v2.31 options.c, k = 21, w = 11)
            (
                "optik",
                [
                    optik,
                    "map",
                    "-r",
                    up,
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
            w, mb, st, _ = run.run(cmd, stdout=sam)
            ok, n, hqw, hq = (
                map_accuracy(sam, reads + ".tsv") if st == "ok" else (0, 0, 0, 0)
            )
            run.row(
                name,
                tool,
                f"map 200k reads (-t {threads}, index incl.)",
                w,
                mb,
                st,
                f"correct={ok}/{n} mapq30={hq} mapq30_wrong={hqw}",
            )
    report_rows(out, first_row, a.report, a.workdir)


# Columns of scaling.tsv, in order, as the results file names them
TSV_COLS = ("genome", "tool", "phase", "wall_s", "rss_MB", "status", "detail")


def report_rows(tsv, first_row, dest, workdir):
    """Writes the rows this run appended to scaling.tsv as the 'scaling' section, one table per genome."""
    by = {}
    with open(tsv) as f:
        for line in list(f)[first_row:]:
            r = dict(zip(TSV_COLS, line.rstrip("\n").split("\t"), strict=False))
            by.setdefault(r.pop("genome"), []).append(r)
    results_file.write_section(
        dest,
        "scaling",
        "Index, graph and mapping cost of every tool as genome size grows "
        "(wall s mean over the measured runs, peak RSS MB)",
        list(by.items()),
        notes=[f"Raw TSV and run metadata: `{os.path.relpath(workdir, R)}`"],
    )


if __name__ == "__main__":
    main()
