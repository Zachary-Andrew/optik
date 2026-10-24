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

R = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(R, "bench"))

import hw
import memcap
import report as results_file
import scaling_bench as sb
from config import PATHS

OPTIK = os.path.join(R, "target/release/optik")
SAL = os.path.join(R, "temp/salmonella_LT2.fa")
CHR1 = os.path.join(R, "temp/scaling/chr1.fa")
HG38 = os.path.join(R, "temp/scaling/hg38_primary.fa")
GUT = os.path.join(R, "temp/published_sets/cf2-gut")
DIGEST_SRC = os.path.join(R, "bench/cpp/unitig_digest.cpp")
SECTIONS = ["reads", "modes", "grch38", "index", "gut"]
# The compact-graph row of grch38 fixes the counting limit so the human node table cannot fit, and needs this much cap above it for the graph itself.
EF_MAX_MEM = "23G"
EF_MIN_CAP = 24 << 30
T = str(hw.cpus())  # thread count of every command, set from --threads


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
        self.run = sb.Runner(
            os.devnull, cap, T, warmup, runs, scratch=os.path.join(work, ".measure")
        )

    def row(self, section, inp, tool, setting, wall, mb, status, extra=""):
        """Append one row to rows.tsv, adding the spread of the measurement it reports."""
        # The spread of the measurement this row reports, whose mean is `wall`
        r = self.run.last
        if r is not None and r.ok and wall is not None and wall == r.wall["mean"]:
            sd = f"{r.wall['sd']:.3g}" if r.wall.get("sd") is not None else "-"
            extra = (f"sd={sd} runs={r.runs} warmup={r.warmup} " + extra).rstrip()
        line = "\t".join(
            [
                section,
                inp,
                tool,
                setting,
                "" if wall is None else f"{wall:.2f}",
                "" if mb is None else f"{mb:.0f}",
                status,
                extra,
            ]
        )
        print(line, flush=True)
        with open(self.rows, "a") as f:
            f.write(line + "\n")

    def cmd(self, name, cmd):
        """Measure `cmd` through bench/measure.py, returning mean wall s, peak MB, status and output.

        The mean covers the --runs that follow the discarded warmups.
        The last run's output goes to logs/<name>.log."""
        r = self.run.run(cmd)
        with open(os.path.join(self.logs, name + ".log"), "w") as f:
            f.write(" ".join(map(str, cmd)) + "\n" + (r[3] or ""))
        return r


def stat(out, key):
    """Value of a 'key<TAB>value' line in optik's report, or an empty string."""
    m = re.search(rf"^{key}\t(\S+)", out or "", re.MULTILINE)
    return m.group(1) if m else ""


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
                os.path.join(R, "bench/gen_reads.py"),
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
    sal = reads(work, SAL, 30, "sal30x.fa", ["--ins", "0.0005", "--del", "0.0005"])
    w, mb, st, o = bench.cmd(
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
            T,
            "-o",
            os.path.join(work, "sal30x.gfa"),
        ],
    )
    bench.row(
        "reads",
        "Salmonella 30x",
        "optik",
        f"-k 31 -c 2 -t {T}, GFA",
        w,
        mb,
        st,
        f"unitigs={stat(o, 'unitigs')}",
    )
    c10 = reads(work, CHR1, 10, "chr1_10x.fa", [])
    for mode in ("part", "table"):
        for lim in ([], ["--max-mem", "4G"]):
            setting = f"--graph {mode}" + (" --max-mem 4G" if lim else "")
            w, mb, st, o = bench.cmd(
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
                    T,
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
                setting + f", -k 31 -c 2 -t {T}",
                w,
                mb,
                st,
                f"unitigs={stat(o, 'unitigs')}",
            )
    a, b = (os.path.join(work, f"c10_{m}.gfa") for m in ("part", "table"))
    same = subprocess.run(["cmp", "-s", a, b], check=False).returncode == 0
    bench.row(
        "reads",
        "chr1 10x",
        "optik dbruijn",
        "part vs table GFA",
        None,
        None,
        "same" if same else "DIFFERENT",
    )
    os.remove(a)
    os.remove(b)
    for sub, extra in (("build-fm", []), ("build-spss", ["--compact"])):
        p = os.path.join(work, f"c10.{sub}")
        w, mb, st, _ = bench.cmd(
            f"chr1_10x_{sub}",
            [OPTIK, sub, *extra, "-i", c10, "-k", "31", "-c", "2", "-t", T, "-o", p],
        )
        bench.row(
            "reads",
            "chr1 10x",
            f"optik {sub}",
            " ".join(extra + [f"-k 31 -c 2 -t {T}"]),
            w,
            mb,
            st,
        )
        if os.path.exists(p):
            os.remove(p)
    return c10


def sec_modes(bench, work, c10):
    """The modes section, partitioned compaction against --graph table at several k and cutoffs."""
    for k in ("31", "21", "15"):
        outs = []
        for mode in ("part", "table"):
            o_ = os.path.join(work, f"chr1_k{k}_{mode}.gfa")
            w, mb, st, o = bench.cmd(
                f"chr1_k{k}_{mode}",
                [
                    OPTIK,
                    "dbruijn",
                    "-i",
                    CHR1,
                    "-k",
                    k,
                    "-c",
                    "1",
                    "-t",
                    T,
                    "--graph",
                    mode,
                    "-o",
                    o_,
                ],
            )
            bench.row(
                "modes",
                "chr1",
                "optik dbruijn",
                f"--graph {mode}, -k {k} -c 1 -t {T}, GFA",
                w,
                mb,
                st,
                f"unitigs={stat(o, 'unitigs')}",
            )
            outs.append(o_)
        same = subprocess.run(["cmp", "-s", *outs], check=False).returncode == 0
        bench.row(
            "modes",
            "chr1",
            "optik dbruijn",
            f"part vs table GFA, k = {k}",
            None,
            None,
            "same" if same else "DIFFERENT",
        )
        for o_ in outs:
            os.remove(o_)
    outs = []
    for mode in ("part", "table"):
        o_ = os.path.join(work, f"c10_k21_{mode}.gfa")
        w, mb, st, o = bench.cmd(
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
                T,
                "--graph",
                mode,
                "-o",
                o_,
            ],
        )
        bench.row(
            "modes",
            "chr1 10x",
            "optik dbruijn",
            f"--graph {mode}, -k 21 -c 3 -t {T}, GFA",
            w,
            mb,
            st,
            f"unitigs={stat(o, 'unitigs')}",
        )
        outs.append(o_)
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
    dg = digest_tool(work)
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
                T,
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
                T,
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
                T,
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
                T,
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
    ref = None
    for tool, setting, cmd, outf in runs:
        if outf == "hg_ef.fa" and memcap.parse_size(bench.cap) < EF_MIN_CAP:
            reason = f"needs a cap of at least {EF_MIN_CAP >> 30} GiB, got {bench.cap}"
            bench.row("grch38", "GRCh38", tool, setting, None, None, "skipped", reason)
            continue
        w, mb, st, o = bench.cmd(f"grch38_{outf}", cmd)
        p = os.path.join(work, outf)
        dig = digest(dg, p) if st == "ok" and os.path.exists(p) else ""
        ref = ref or dig
        bench.row(
            "grch38",
            "GRCh38",
            tool,
            setting + f", -k 31 -t {T}",
            w,
            mb,
            st,
            f"digest={dig}"
            + ("" if not dig or dig == ref else " DIFFERENT")
            + (f" links={stat(o, 'links')}" if outf.endswith(".gfa") else ""),
        )
        if os.path.exists(p):
            os.remove(p)


def windows(work):
    """Write 20,000 random N-free 150 bp windows of GRCh38 with seed 7, and the same windows with one substitution each at a random position."""
    ex, sb_ = os.path.join(work, "win_exact.fa"), os.path.join(work, "win_sub.fa")
    if os.path.exists(sb_):
        return ex, sb_
    seqs, cur = [], []
    with open(HG38) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    lens = [len(s) for s in seqs]
    tot = sum(lens)
    r = random.Random(7)
    with open(ex, "w") as e, open(sb_, "w") as s:
        n = 0
        while n < 20000:
            x = r.randrange(tot)
            i = 0
            while x >= lens[i]:
                x -= lens[i]
                i += 1
            if x + 150 > lens[i]:
                continue
            w = seqs[i][x : x + 150]
            if set(w) - set("ACGT"):
                continue
            p = r.randrange(150)
            m = w[:p] + r.choice([b for b in "ACGT" if b != w[p]]) + w[p + 1 :]
            e.write(f">{n}\n{w}\n")
            s.write(f">{n}\n{m}\n")
            n += 1
    return ex, sb_


def sec_index(bench, work):
    """The index section, which builds both indexes on GRCh38 and queries exact and substituted windows."""
    idx, res = {}, []
    for sub, extra in (("build-fm", []), ("build-spss", ["--compact"])):
        p = os.path.join(work, f"hg.{sub}")
        w, mb, st, o = bench.cmd(
            f"grch38_{sub}",
            [OPTIK, sub, *extra, "-i", HG38, "-k", "31", "-t", T, "-o", p],
        )
        n = re.search(r"distinct_kmers\t(\d+)", o or "")
        res.append((sub, extra, w, mb, st, int(n.group(1)) if n else None))
        if st == "ok":
            idx[sub] = p
    # Bits per distinct k-mer, both over the count the builds report (build-spss prints it)
    nk = next((r[5] for r in res if r[5]), None)
    for sub, extra, w, mb, st, n in res:
        bits = f"{os.path.getsize(idx[sub]) * 8 / nk:.2f}" if sub in idx and nk else ""
        bench.row(
            "index",
            "GRCh38",
            f"optik {sub}",
            " ".join(extra + [f"-k 31 -t {T}"]),
            w,
            mb,
            st,
            f"distinct_kmers={n} bits_per_kmer={bits}",
        )
    ex, su = windows(work)
    got = {}
    for sub, p in idx.items():
        for wl, q in (("exact", ex), ("one substitution", su)):
            w, mb, st, o = bench.cmd(
                f"grch38_{sub}_query_{wl.split()[0]}",
                [OPTIK, "query", "-i", p, "-q", q, "-k", "31", "-t", T],
            )
            pos, tot = sb.kv(o, "num_positive_kmers"), sb.kv(o, "num_kmers")
            got[(sub, wl)] = pos
            ok = st == "ok" and tot == 20000 * 120 and (wl != "exact" or pos == tot)
            bench.row(
                "index",
                "GRCh38 windows",
                f"optik {sub}",
                f"query {wl}",
                w,
                mb,
                st if ok else "FAIL",
                f"kmers={tot} present={pos}",
            )
    if len(idx) == 2:
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
    for p in idx.values():
        os.remove(p)


class DiskPeak(threading.Thread):
    """Filesystem growth under `path`, sampled every second through statvfs."""

    def __init__(self, path):
        super().__init__(daemon=True)
        self.path, self.peak, self.stop = path, 0, threading.Event()
        self.base = self.used()

    def used(self):
        """Bytes in use on the filesystem."""
        st = os.statvfs(self.path)
        return (st.f_blocks - st.f_bfree) * st.f_frsize

    def run(self):
        while not self.stop.wait(1.0):
            self.peak = max(self.peak, self.used() - self.base)


def sec_gut(bench, work):
    """The gut section, optik and Cuttlefish 3 on Cuttlefish 2's gut read set."""
    # The files the fetch recorded, since PROVENANCE.json lists nine paired individuals in 18 files
    want = (
        sorted(
            f["path"]
            for f in json.loads(Path(os.path.join(GUT, "PROVENANCE.json")).read_text())[
                "files"
            ]
        )
        if os.path.exists(os.path.join(GUT, "PROVENANCE.json"))
        else []
    )
    ins = [f for f in want if os.path.exists(f)]
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
    runs = [
        (
            "optik",
            "default",
            [
                OPTIK,
                "dbruijn",
                *[a for f in ins for a in ("-i", f)],
                "-k",
                "27",
                "-c",
                "2",
                "-t",
                T,
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
                *[a for f in ins for a in ("-i", f)],
                "-k",
                "27",
                "-c",
                "2",
                "-t",
                T,
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
                *[a for f in ins for a in ("-s", f)],
                "-k",
                "27",
                "-t",
                T,
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
        dp = DiskPeak(work)
        dp.start()
        w, mb, st, o = bench.cmd(f"gut_{tool}_{setting.split()[0].strip('-')}", cmd)
        dp.stop.set()
        dp.join()
        r = bench.run.last
        io = (
            f"written_gb={r.written / 1e9:.1f} read_gb={r.read / 1e9:.1f}"
            if r is not None and r.ok and r.written is not None
            else ""
        )
        bench.row(
            "gut",
            "cf2-gut",
            tool,
            setting + f", -k 27 -c 2 -t {T}",
            w,
            mb,
            st,
            f"{io} peak_disk_gb={dp.peak / 1e9:.0f} unitigs={stat(o, 'unitigs')}",
        )
    a, b = os.path.join(work, "gut.fa"), os.path.join(work, "gut13.fa")
    if os.path.exists(a) and os.path.exists(b):
        same = subprocess.run(["cmp", "-s", a, b], check=False).returncode == 0
        bench.row(
            "gut",
            "cf2-gut",
            "optik",
            "default vs --max-mem 13G output",
            None,
            None,
            "same" if same else "DIFFERENT",
        )
    for x in glob.glob(os.path.join(work, "gut*.fa")):
        os.remove(x)


def main():
    """Entry point that runs the sections named by --only."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("out")
    ap.add_argument("--only", default=",".join(SECTIONS))
    ap.add_argument(
        "--cap",
        default=hw.default_cap_str(),
        help="per-command resident-memory cap (memcap; default: usable memory less 2 GiB)",
    )
    ap.add_argument("--work", default=os.path.join(R, "temp/dbg_bench_work"))
    ap.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every command (default: the usable CPUs)",
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
    global T
    T = str(a.threads)
    only = a.only.split(",")
    bad = [s for s in only if s not in SECTIONS]
    if bad:
        sys.exit(f"unknown section(s) {bad}; known: {SECTIONS}")
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    rev = subprocess.run(
        ["git", "-C", R, "describe", "--always", "--dirty"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    sha = hashlib.sha256(Path(OPTIK).read_bytes()).hexdigest()
    with open(os.path.join(a.out, "meta.txt"), "w") as f:
        f.write(
            f"git {rev}\noptik sha256 {sha}\ncap {a.cap}\nthreads {T}\nwarmup {a.warmup}\nruns {a.runs}\n"
            f"started {time.ctime()}\n"
            f"cmd {' '.join(sys.argv)}\n"
        )
    bench = Bench(a.out, a.cap, a.work, a.warmup, a.runs)
    # Rows of earlier runs into the same directory stay in rows.tsv but not in this run's report.
    first_row = 0
    if os.path.exists(bench.rows):
        with open(bench.rows) as fh:
            first_row = sum(1 for _ in fh)
    c10 = None
    if "reads" in only or "modes" in only:
        c10 = (
            sec_reads(bench, a.work)
            if "reads" in only
            else reads(a.work, CHR1, 10, "chr1_10x.fa", [])
        )
    if "modes" in only:
        sec_modes(bench, a.work, c10)
    if "grch38" in only:
        sec_grch38(bench, a.work)
    if "index" in only:
        sec_index(bench, a.work)
    if "gut" in only:
        sec_gut(bench, a.work)
    with open(os.path.join(a.out, "meta.txt"), "a") as f:
        f.write(f"finished {time.ctime()}\n")
    by = {}
    with open(bench.rows) as f:
        for line in list(f)[first_row:]:
            r = dict(zip(TSV_COLS, line.rstrip("\n").split("\t"), strict=False))
            by.setdefault(r.pop("section"), []).append(r)
    results_file.write_section(
        a.report,
        "dbg",
        "optik dbruijn against Cuttlefish, GGCAT and Bifrost on read sets, modes, GRCh38, "
        "index building and a gut metagenome (wall s mean over the measured runs, peak RSS MB)",
        list(by.items()),
        notes=[
            f"Raw TSV, metadata and per-command logs: `{os.path.relpath(a.out, R)}`"
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
