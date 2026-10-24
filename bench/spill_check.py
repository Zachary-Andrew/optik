#!/usr/bin/env python3
"""spill_check.py - checks that the k-mer counter's disk spill (--max-mem, --tmp-dir) changes where the partition buffers live but never the result.

Every command that counts k-mers (dbruijn -c 1 and -c 2, build-fm, build-spss --compact) runs at 1, 4 and 16 threads with its buffers forced to disk (OPTIK_SPILL_TEST, a tiny per-thread share so every worker spills many times).
The output must be byte-identical to the in-memory run, the spill must have happened, and no file may be left in the spill directory.
The real paths follow.
They are --max-mem on a read set, spilling in pass 1 and before the graph is allocated, a graph larger than --max-mem, which must give a clear error rather than an OOM kill, an unusable spill directory, which must give an error rather than a crash, and malformed flags.

Usage: python3 bench/spill_check.py [--optik BIN] [--genome FASTA] [--workdir DIR]
"""

import argparse
import filecmp
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from config import DATASETS, PATHS

FAILS = []


def check(name, ok, detail=""):
    """Log one check, recording it in FAILS when it fails."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def run(cmd, env=None):
    e = dict(os.environ)
    for v in ("OPTIK_SPILL_TEST", "OPTIK_DBG_MODE", "OPTIK_DBG_COMPACT"):
        e.pop(v, None)
    if env:
        e.update(env)
    return subprocess.run(
        [str(c) for c in cmd], capture_output=True, text=True, env=e, check=False
    )


def main():
    """Run every spill check and exit non-zero on any failure."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--optik", default=str(PATHS["optik"]))
    ap.add_argument("--genome", default=str(DATASETS["ecoli_k31"].get("dbg_input")))
    ap.add_argument("--workdir", default="/tmp/optik_spill_check")
    a = ap.parse_args()
    optik = a.optik
    wd = Path(a.workdir)
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    g = wd / "genome.fa"
    opener = gzip.open if a.genome.endswith(".gz") else open
    with opener(a.genome, "rt") as fi, open(g, "w") as fo:
        shutil.copyfileobj(fi, fo)
    reads = wd / "reads.fa"
    subprocess.run(
        [
            sys.executable,
            str(HERE / "gen_reads.py"),
            "-r",
            str(g),
            "-c",
            "12",
            "-l",
            "150",
            "-e",
            "0.005",
            "-s",
            "5",
            "-o",
            str(reads),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    sp = wd / "spill"
    sp.mkdir()

    cmds = {
        "dbruijn -c 1": lambda out, t: [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "1",
            "-t",
            t,
            "-o",
            out,
        ],
        "dbruijn -c 2 (k=21)": lambda out, t: [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "21",
            "-c",
            "2",
            "-t",
            t,
            "-o",
            out,
        ],
        "build-fm -c 2": lambda out, t: [
            optik,
            "build-fm",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            t,
            "-o",
            out,
        ],
        "build-spss --compact -c 2": lambda out, t: [
            optik,
            "build-spss",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "--compact",
            "-t",
            t,
            "-o",
            out,
        ],
    }
    for name, mk in cmds.items():
        for t in ("1", "4", "16"):
            mem_out, sp_out = wd / "mem.out", wd / "spill.out"
            p1 = run(mk(mem_out, t))
            p2 = run(mk(sp_out, t) + ["--tmp-dir", sp], {"OPTIK_SPILL_TEST": "65536"})
            same = (
                p1.returncode == 0
                and p2.returncode == 0
                and filecmp.cmp(mem_out, sp_out, shallow=False)
            )
            check(
                f"{name} -t {t}: spilled output byte-identical to in-memory",
                same and "spilling" in p2.stderr and not any(sp.iterdir()),
                f"rc={p1.returncode}/{p2.returncode} spilled={'spilling' in p2.stderr} "
                f"left={len(list(sp.iterdir()))} {p2.stderr.strip()[-160:] if p2.returncode else ''}",
            )

    # The compact graph (Elias-Fano keys), chosen when the table does not fit and forced here.
    # Every output must equal the table's, with its solid-k-mer buckets both in memory and spilled.
    for name, mk in cmds.items():
        for env in (
            {"OPTIK_DBG_MODE": "table", "OPTIK_DBG_COMPACT": "1"},
            {
                "OPTIK_DBG_MODE": "table",
                "OPTIK_DBG_COMPACT": "1",
                "OPTIK_SPILL_TEST": "65536",
            },
        ):
            mem_out, cp_out = wd / "mem.out", wd / "compact.out"
            p1 = run(mk(mem_out, "8"), {"OPTIK_DBG_MODE": "table"})
            p2 = run(mk(cp_out, "8") + ["--tmp-dir", sp], env)
            same = (
                p1.returncode == 0
                and p2.returncode == 0
                and filecmp.cmp(mem_out, cp_out, shallow=False)
            )
            check(
                f"{name}: compact graph{' + spill' if 'OPTIK_SPILL_TEST' in env else ''} output "
                f"byte-identical to the table's",
                same and "built compact" in p2.stderr and not any(sp.iterdir()),
                f"rc={p1.returncode}/{p2.returncode} {p2.stderr.strip()[-160:] if p2.returncode else ''}",
            )
    k1, k2 = wd / "t.kdbg", wd / "c.kdbg"
    p1 = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-o",
            wd / "t.gfa",
            "-s",
            k1,
        ]
    )
    p2 = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-o",
            wd / "c.gfa",
            "-s",
            k2,
        ],
        {"OPTIK_DBG_MODE": "table", "OPTIK_DBG_COMPACT": "1"},
    )
    check(
        "dbruijn -s: the compact graph's archive equals the table's",
        p1.returncode == p2.returncode == 0 and filecmp.cmp(k1, k2, shallow=False),
    )

    # Partitioned compaction (dbg_part.hpp, the default from k = 19).
    # Every output must equal the table's at 1, 4 and 16 threads, with its fragment and end-key streams in memory and spilled, at odd and even k.
    part_cmds = dict(cmds)
    part_cmds["dbruijn -c 1 (k=16, palindromes)"] = lambda out, t: [
        optik,
        "dbruijn",
        "-i",
        reads,
        "-k",
        "16",
        "-c",
        "1",
        "-t",
        t,
        "-o",
        out,
    ]
    part_cmds["dbruijn -c 3 (k=25)"] = lambda out, t: [
        optik,
        "dbruijn",
        "-i",
        reads,
        "-k",
        "25",
        "-c",
        "3",
        "-t",
        t,
        "-o",
        out,
    ]
    for name, mk in part_cmds.items():
        for t in ("1", "4", "16"):
            for extra in ({}, {"OPTIK_SPILL_TEST": "65536"}):
                a_out, b_out = wd / "table.out", wd / "part.out"
                p1 = run(mk(a_out, t), {"OPTIK_DBG_MODE": "table"})
                p2 = run(
                    mk(b_out, t) + ["--tmp-dir", sp], dict(extra, OPTIK_DBG_MODE="part")
                )
                same = (
                    p1.returncode == 0
                    and p2.returncode == 0
                    and filecmp.cmp(a_out, b_out, shallow=False)
                )
                check(
                    f"{name} -t {t}: partitioned{' + spill' if extra else ''} == table",
                    same and not any(sp.iterdir()),
                    f"rc={p1.returncode}/{p2.returncode} {p2.stderr.strip()[-160:] if p2.returncode else ''}",
                )

    # The real budget path, with buffers past --max-mem and then the spill before insertion
    base = wd / "base.gfa"
    run([optik, "dbruijn", "-i", reads, "-k", "31", "-c", "2", "-t", "8", "-o", base])
    lim = wd / "lim.gfa"
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            "8",
            "-o",
            lim,
            "--max-mem",
            "100M",
            "--tmp-dir",
            sp,
        ]
    )
    check(
        "--max-mem 100M on a read set: spills, same GFA",
        p.returncode == 0
        and "spilling" in p.stderr
        and filecmp.cmp(base, lim, shallow=False),
        f"rc={p.returncode} {p.stderr.strip()[-160:]}",
    )
    # The default spill directory is the output's.
    od = wd / "outdir"
    od.mkdir()
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            "8",
            "-o",
            od / "g.gfa",
        ],
        {"OPTIK_SPILL_TEST": "65536"},
    )
    check(
        "default spill directory is the output's",
        p.returncode == 0 and f"spilling them to {od}" in p.stderr,
        p.stderr.strip()[-160:],
    )
    # A table larger than --max-mem switches to the compact graph, with the same GFA.
    t1, t2 = wd / "tab.gfa", wd / "lim_compact.gfa"
    run([optik, "dbruijn", "-i", g, "-k", "31", "-c", "1", "-o", t1])
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            g,
            "-k",
            "31",
            "-c",
            "1",
            "-o",
            t2,
            "--max-mem",
            "64M",
            "--tmp-dir",
            sp,
        ],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "table larger than --max-mem: built compact, same GFA",
        p.returncode == 0
        and "built compact" in p.stderr
        and filecmp.cmp(t1, t2, shallow=False),
        f"rc={p.returncode} {p.stderr.strip()[-200:]}",
    )
    # Even the compact graph larger than the limit gives a clear error rather than an OOM kill.
    import random as _r

    rg = _r.Random(4)
    big = wd / "big.fa"
    with open(big, "w") as fh:
        fh.write(">r\n" + "".join(rg.choice("ACGT") for _ in range(12_000_000)) + "\n")
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            big,
            "-k",
            "31",
            "-c",
            "1",
            "-o",
            wd / "x.gfa",
            "--max-mem",
            "64M",
            "--tmp-dir",
            sp,
        ],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "compact graph larger than --max-mem: error naming the need, no OOM",
        p.returncode == 1
        and "more than the memory limit" in p.stderr
        and "--max-mem" in p.stderr,
        f"rc={p.returncode} {p.stderr.strip()[-200:]}",
    )
    # Partitioned compaction keeps that genome under 64 MB of counting state.
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            big,
            "-k",
            "31",
            "-c",
            "1",
            "-o",
            wd / "x_part.gfa",
            "--max-mem",
            "64M",
            "--tmp-dir",
            sp,
        ],
        {"OPTIK_DBG_MODE": "part"},
    )
    p_ref = run(
        [optik, "dbruijn", "-i", big, "-k", "31", "-c", "1", "-o", wd / "x_ref.gfa"],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "partitioned compaction under --max-mem 64M: same GFA as the unlimited table",
        p.returncode == 0
        and p_ref.returncode == 0
        and filecmp.cmp(wd / "x_part.gfa", wd / "x_ref.gfa", shallow=False),
        f"rc={p.returncode} {p.stderr.strip()[-160:]}",
    )
    # Unusable spill directories
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-o",
            wd / "x.gfa",
            "--tmp-dir",
            wd / "missing",
        ]
    )
    check(
        "--tmp-dir that does not exist: error",
        p.returncode == 1 and "does not exist" in p.stderr,
        p.stderr.strip()[-160:],
    )
    ro = wd / "readonly"
    ro.mkdir()
    ro.chmod(0o500)
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-o",
            wd / "x.gfa",
            "--tmp-dir",
            ro,
        ],
        {"OPTIK_SPILL_TEST": "65536"},
    )
    ro.chmod(0o700)
    writable_anyway = os.access(ro, os.W_OK) and os.geteuid() == 0
    check(
        "unwritable --tmp-dir when a spill is needed: error, not a crash",
        writable_anyway
        or (p.returncode == 1 and "cannot create a spill file" in p.stderr),
        f"rc={p.returncode} {p.stderr.strip()[-160:]}",
    )
    # A spill write that fails mid-run (a file size limit gives EFBIG, standing in for a full disk) must give an error naming it, with no crash and no partial output.
    import resource
    import signal

    def small_fsize():
        """Child pre-exec hook that limits file size to 200 kB and ignores SIGXFSZ."""
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (200_000, 200_000))

    e = dict(os.environ, OPTIK_SPILL_TEST="65536")
    for cmd in (
        ["dbruijn", "-o", str(wd / "f.gfa")],
        ["build-fm", "-o", str(wd / "f.fm")],
    ):
        p = subprocess.run(
            [optik]
            + cmd
            + ["-i", str(reads), "-k", "31", "-t", "8", "--tmp-dir", str(sp)],
            env=e,
            preexec_fn=small_fsize,
            capture_output=True,
            text=True,
            check=False,
        )
        check(
            f"{cmd[0]}: a failing spill write is an error (not a crash or a short result)",
            p.returncode == 1
            and "spill write failed" in p.stderr
            and not Path(cmd[2]).exists(),
            f"rc={p.returncode} {p.stderr.strip()[-160:]}",
        )
    # Several inputs spill like one.
    half = wd / "h1.fa", wd / "h2.fa"
    lines = Path(reads).read_text().splitlines(True)
    cut = len(lines) // 4 * 2  # on a record boundary, with 2 lines per read
    Path(half[0]).write_text("".join(lines[:cut]))
    Path(half[1]).write_text("".join(lines[cut:]))
    p = run(
        [
            optik,
            "dbruijn",
            "-i",
            half[0],
            "-i",
            half[1],
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            "8",
            "-o",
            wd / "two.gfa",
            "--tmp-dir",
            sp,
        ],
        {"OPTIK_SPILL_TEST": "65536"},
    )
    check(
        "-i R1 -i R2 with a spill == one input in memory",
        p.returncode == 0
        and "spilling" in p.stderr
        and filecmp.cmp(base, wd / "two.gfa", shallow=False),
        p.stderr.strip()[-160:],
    )
    for bad in ("0", "-5G", "abc", "10K"):
        p = run([optik, "dbruijn", "-i", reads, "-o", wd / "x.gfa", "--max-mem", bad])
        check(
            f"--max-mem {bad}: rejected",
            p.returncode == 1 and "--max-mem" in p.stderr,
            p.stderr.strip()[-120:],
        )
    p = run([optik, "build-spss", "-i", g, "-o", wd / "x.spss", "--max-mem", "1G"])
    check(
        "build-spss --max-mem without --compact: rejected",
        p.returncode == 1 and "--compact" in p.stderr,
        p.stderr.strip()[-120:],
    )

    print(f"\nspill_check: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}")
    if not FAILS:
        shutil.rmtree(wd, ignore_errors=True)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
