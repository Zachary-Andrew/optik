#!/usr/bin/env python3
"""spill_check.py - checks that the k-mer counter's disk spill (--max-mem, --tmp-dir) changes where the partition buffers are kept but never the result.

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
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from checklog import CheckLog
from config import (
    DATASETS,
    PATHS,
)
from diagnostics import DETAIL_CHARS

# optik's exit status for a rejected flag or a clean runtime error
OPTIK_ERROR_EXIT = 1
# a random single-record genome whose graph exceeds the tiny --max-mem of the out-of-memory check
BIG_GENOME_SEED = 4
BIG_GENOME_LEN = 12_000_000
# permissions of the spill directory made unwritable, and restored afterwards so it can be removed
READONLY_DIR_MODE = 0o500
OWNER_RWX_MODE = 0o700
# file size limit (bytes) that makes the spill run out of disk space
FSIZE_LIMIT_BYTES = 200_000
# lines of a FASTA read record
FASTA_LINES_PER_READ = 2

check = CheckLog()
FAILS = check.failed


def run(cmd, env=None):
    child_env = dict(os.environ)
    for var in ("OPTIK_SPILL_TEST", "OPTIK_DBG_MODE", "OPTIK_DBG_COMPACT"):
        child_env.pop(var, None)
    if env:
        child_env.update(env)
    return subprocess.run(
        [str(arg) for arg in cmd], capture_output=True, text=True, env=child_env, check=False
    )


@dataclass
class SpillInputs:
    """Holds the binary, the work directory and the inputs that every spill section runs on."""

    optik: str
    work_dir: Path
    genome_fa: Path
    reads: Path
    spill_dir: Path
    cmds: dict
    base: Path


def check_forced_spill(fx):
    """Forces every counting command to spill at 1, 4 and 16 threads and compares the output with the in-memory run."""
    work_dir, spill_dir, cmds = fx.work_dir, fx.spill_dir, fx.cmds
    for name, make_cmd in cmds.items():
        for threads in ("1", "4", "16"):
            mem_out, sp_out = work_dir / "mem.out", work_dir / "spill.out"
            mem_proc = run(make_cmd(mem_out, threads))
            spill_proc = run(make_cmd(sp_out, threads) + ["--tmp-dir", spill_dir], {"OPTIK_SPILL_TEST": "65536"})
            same = (
                mem_proc.returncode == 0
                and spill_proc.returncode == 0
                and filecmp.cmp(mem_out, sp_out, shallow=False)
            )
            check(
                f"{name} -t {threads}: spilled output byte-identical to in-memory",
                same and "spilling" in spill_proc.stderr and not any(spill_dir.iterdir()),
                f"rc={mem_proc.returncode}/{spill_proc.returncode} spilled={'spilling' in spill_proc.stderr} "
                f"left={len(list(spill_dir.iterdir()))} {spill_proc.stderr.strip()[-DETAIL_CHARS:] if spill_proc.returncode else ''}",
            )


def check_compact_graph(fx):
    """Forces the compact Elias-Fano graph and checks it against the table, with its buckets both in memory and spilled."""
    optik, work_dir, reads, spill_dir, cmds = fx.optik, fx.work_dir, fx.reads, fx.spill_dir, fx.cmds
    # The compact graph (Elias-Fano keys), chosen when the table does not fit and forced here.
    # Every output must equal the table's, with its solid-k-mer buckets both in memory and spilled.
    for name, make_cmd in cmds.items():
        for env in (
            {"OPTIK_DBG_MODE": "table", "OPTIK_DBG_COMPACT": "1"},
            {
                "OPTIK_DBG_MODE": "table",
                "OPTIK_DBG_COMPACT": "1",
                "OPTIK_SPILL_TEST": "65536",
            },
        ):
            mem_out, cp_out = work_dir / "mem.out", work_dir / "compact.out"
            mem_proc = run(make_cmd(mem_out, "8"), {"OPTIK_DBG_MODE": "table"})
            spill_proc = run(make_cmd(cp_out, "8") + ["--tmp-dir", spill_dir], env)
            same = (
                mem_proc.returncode == 0
                and spill_proc.returncode == 0
                and filecmp.cmp(mem_out, cp_out, shallow=False)
            )
            check(
                f"{name}: compact graph{' + spill' if 'OPTIK_SPILL_TEST' in env else ''} output "
                f"byte-identical to the table's",
                same and "built compact" in spill_proc.stderr and not any(spill_dir.iterdir()),
                f"rc={mem_proc.returncode}/{spill_proc.returncode} {spill_proc.stderr.strip()[-DETAIL_CHARS:] if spill_proc.returncode else ''}",
            )
    table_kdbg, compact_kdbg = work_dir / "t.kdbg", work_dir / "c.kdbg"
    mem_proc = run(
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
            work_dir / "t.gfa",
            "-s",
            table_kdbg,
        ]
    )
    spill_proc = run(
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
            work_dir / "c.gfa",
            "-s",
            compact_kdbg,
        ],
        {"OPTIK_DBG_MODE": "table", "OPTIK_DBG_COMPACT": "1"},
    )
    check(
        "dbruijn -s: the compact graph's archive equals the table's",
        mem_proc.returncode == spill_proc.returncode == 0 and filecmp.cmp(table_kdbg, compact_kdbg, shallow=False),
    )


def check_partitioned_compaction(fx):
    """Runs partitioned compaction at 1, 4 and 16 threads, in memory and spilled, and checks it against the table."""
    optik, work_dir, reads, spill_dir, cmds = fx.optik, fx.work_dir, fx.reads, fx.spill_dir, fx.cmds
    # Partitioned compaction (dbg_part.hpp, the default from k = 19).
    # Every output must equal the table's at 1, 4 and 16 threads, with its fragment and end-key streams in memory and spilled, at odd and even k.
    part_cmds = dict(cmds)
    part_cmds["dbruijn -c 1 (k=16, palindromes)"] = lambda out, threads: [
        optik,
        "dbruijn",
        "-i",
        reads,
        "-k",
        "16",
        "-c",
        "1",
        "-t",
        threads,
        "-o",
        out,
    ]
    part_cmds["dbruijn -c 3 (k=25)"] = lambda out, threads: [
        optik,
        "dbruijn",
        "-i",
        reads,
        "-k",
        "25",
        "-c",
        "3",
        "-t",
        threads,
        "-o",
        out,
    ]
    for name, make_cmd in part_cmds.items():
        for threads in ("1", "4", "16"):
            for extra in ({}, {"OPTIK_SPILL_TEST": "65536"}):
                a_out, b_out = work_dir / "table.out", work_dir / "part.out"
                mem_proc = run(make_cmd(a_out, threads), {"OPTIK_DBG_MODE": "table"})
                spill_proc = run(
                    make_cmd(b_out, threads) + ["--tmp-dir", spill_dir], dict(extra, OPTIK_DBG_MODE="part")
                )
                same = (
                    mem_proc.returncode == 0
                    and spill_proc.returncode == 0
                    and filecmp.cmp(a_out, b_out, shallow=False)
                )
                check(
                    f"{name} -t {threads}: partitioned{' + spill' if extra else ''} == table",
                    same and not any(spill_dir.iterdir()),
                    f"rc={mem_proc.returncode}/{spill_proc.returncode} {spill_proc.stderr.strip()[-DETAIL_CHARS:] if spill_proc.returncode else ''}",
                )


def check_budget_path(fx):
    """Pushes a read set past --max-mem so it spills in pass 1 and again before the graph is allocated."""
    optik, work_dir, reads, spill_dir, base = fx.optik, fx.work_dir, fx.reads, fx.spill_dir, fx.base
    # The real budget path, with buffers past --max-mem and then the spill before insertion
    run([optik, "dbruijn", "-i", reads, "-k", "31", "-c", "2", "-t", "8", "-o", base])
    lim = work_dir / "lim.gfa"
    proc = run(
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
            spill_dir,
        ]
    )
    check(
        "--max-mem 100M on a read set: spills, same GFA",
        proc.returncode == 0
        and "spilling" in proc.stderr
        and filecmp.cmp(base, lim, shallow=False),
        f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )


def check_default_spill_dir(fx):
    """Confirms the spill directory defaults to the output's directory."""
    optik, work_dir, reads = fx.optik, fx.work_dir, fx.reads
    # The default spill directory is the output's.
    out_dir = work_dir / "outdir"
    out_dir.mkdir()
    proc = run(
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
            out_dir / "g.gfa",
        ],
        {"OPTIK_SPILL_TEST": "65536"},
    )
    check(
        "default spill directory is the output's",
        proc.returncode == 0 and f"spilling them to {out_dir}" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )


def check_table_over_limit(fx):
    """Lets a table larger than --max-mem switch to the compact graph, which must give the same GFA."""
    optik, work_dir, genome_fa, spill_dir = fx.optik, fx.work_dir, fx.genome_fa, fx.spill_dir
    # A table larger than --max-mem switches to the compact graph, with the same GFA.
    table_gfa, compact_gfa = work_dir / "tab.gfa", work_dir / "lim_compact.gfa"
    run([optik, "dbruijn", "-i", genome_fa, "-k", "31", "-c", "1", "-o", table_gfa])
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            genome_fa,
            "-k",
            "31",
            "-c",
            "1",
            "-o",
            compact_gfa,
            "--max-mem",
            "64M",
            "--tmp-dir",
            spill_dir,
        ],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "table larger than --max-mem: built compact, same GFA",
        proc.returncode == 0
        and "built compact" in proc.stderr
        and filecmp.cmp(table_gfa, compact_gfa, shallow=False),
        f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )


def check_graph_over_limit(fx):
    """Expects a clear error when even the compact graph exceeds the limit, and checks that partitioned compaction still fits."""
    optik, work_dir, spill_dir = fx.optik, fx.work_dir, fx.spill_dir
    # Even the compact graph larger than the limit gives a clear error rather than an OOM kill.
    import random as random_mod

    rng = random_mod.Random(BIG_GENOME_SEED)
    big = work_dir / "big.fa"
    with open(big, "w") as handle:
        handle.write(">r\n" + "".join(rng.choice("ACGT") for _ in range(BIG_GENOME_LEN)) + "\n")
    proc = run(
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
            work_dir / "x.gfa",
            "--max-mem",
            "64M",
            "--tmp-dir",
            spill_dir,
        ],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "compact graph larger than --max-mem: error naming the need, no OOM",
        proc.returncode == OPTIK_ERROR_EXIT
        and "more than the memory limit" in proc.stderr
        and "--max-mem" in proc.stderr,
        f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )
    # Partitioned compaction keeps that genome under 64 MB of counting state.
    proc = run(
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
            work_dir / "x_part.gfa",
            "--max-mem",
            "64M",
            "--tmp-dir",
            spill_dir,
        ],
        {"OPTIK_DBG_MODE": "part"},
    )
    p_ref = run(
        [optik, "dbruijn", "-i", big, "-k", "31", "-c", "1", "-o", work_dir / "x_ref.gfa"],
        {"OPTIK_DBG_MODE": "table"},
    )
    check(
        "partitioned compaction under --max-mem 64M: same GFA as the unlimited table",
        proc.returncode == 0
        and p_ref.returncode == 0
        and filecmp.cmp(work_dir / "x_part.gfa", work_dir / "x_ref.gfa", shallow=False),
        f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )


def check_unusable_spill_dirs(fx):
    """Points the spill at directories that can't be used, which must fail with an error."""
    optik, work_dir, reads = fx.optik, fx.work_dir, fx.reads
    # Unusable spill directories
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-o",
            work_dir / "x.gfa",
            "--tmp-dir",
            work_dir / "missing",
        ]
    )
    check(
        "--tmp-dir that does not exist: error",
        proc.returncode == OPTIK_ERROR_EXIT and "does not exist" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    readonly_dir = work_dir / "readonly"
    readonly_dir.mkdir()
    readonly_dir.chmod(READONLY_DIR_MODE)
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-o",
            work_dir / "x.gfa",
            "--tmp-dir",
            readonly_dir,
        ],
        {"OPTIK_SPILL_TEST": "65536"},
    )
    readonly_dir.chmod(OWNER_RWX_MODE)
    writable_anyway = os.access(readonly_dir, os.W_OK) and os.geteuid() == 0
    check(
        "unwritable --tmp-dir when a spill is needed: error, not a crash",
        writable_anyway
        or (proc.returncode == OPTIK_ERROR_EXIT and "cannot create a spill file" in proc.stderr),
        f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )


def check_failing_spill_write(fx):
    """Makes a spill write fail mid-run, standing in for a full disk, and expects an error."""
    optik, work_dir, reads, spill_dir = fx.optik, fx.work_dir, fx.reads, fx.spill_dir
    # A spill write that fails mid-run (a file size limit gives EFBIG, standing in for a full disk) must give an error naming it, with no crash and no partial output.
    import resource
    import signal

    def small_fsize():
        """Child pre-exec hook that limits file size to 200 kB and ignores SIGXFSZ."""
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (FSIZE_LIMIT_BYTES, FSIZE_LIMIT_BYTES))

    spill_env = dict(os.environ, OPTIK_SPILL_TEST="65536")
    for cmd in (
        ["dbruijn", "-o", str(work_dir / "f.gfa")],
        ["build-fm", "-o", str(work_dir / "f.fm")],
    ):
        proc = subprocess.run(
            [optik]
            + cmd
            + ["-i", str(reads), "-k", "31", "-t", "8", "--tmp-dir", str(spill_dir)],
            env=spill_env,
            preexec_fn=small_fsize,
            capture_output=True,
            text=True,
            check=False,
        )
        check(
            f"{cmd[0]}: a failing spill write is an error (not a crash or a short result)",
            proc.returncode == OPTIK_ERROR_EXIT
            and "spill write failed" in proc.stderr
            and not Path(cmd[cmd.index("-o") + 1]).exists(),
            f"rc={proc.returncode} {proc.stderr.strip()[-DETAIL_CHARS:]}",
        )


def check_several_inputs(fx):
    """Confirms that several inputs spill exactly like one joined file."""
    optik, work_dir, genome_fa, reads, spill_dir, base = fx.optik, fx.work_dir, fx.genome_fa, fx.reads, fx.spill_dir, fx.base
    # Several inputs spill like one.
    half = work_dir / "h1.fa", work_dir / "h2.fa"
    lines = Path(reads).read_text().splitlines(True)
    cut = len(lines) // (len(half) * FASTA_LINES_PER_READ) * FASTA_LINES_PER_READ  # on a record boundary, halfway
    Path(half[0]).write_text("".join(lines[:cut]))
    Path(half[1]).write_text("".join(lines[cut:]))
    proc = run(
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
            work_dir / "two.gfa",
            "--tmp-dir",
            spill_dir,
        ],
        # A tiny test budget forces spilling on a small input.
        {"OPTIK_SPILL_TEST": "65536"},
    )
    check(
        "-i R1 -i R2 with a spill == one input in memory",
        proc.returncode == 0
        and "spilling" in proc.stderr
        and filecmp.cmp(base, work_dir / "two.gfa", shallow=False),
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    # 10K parses but falls below the 64M floor.
    for bad in ("0", "-5G", "abc", "10K"):
        proc = run([optik, "dbruijn", "-i", reads, "-o", work_dir / "x.gfa", "--max-mem", bad])
        check(
            f"--max-mem {bad}: rejected",
            proc.returncode == OPTIK_ERROR_EXIT and "--max-mem" in proc.stderr,
            proc.stderr.strip()[-DETAIL_CHARS:],
        )
    proc = run([optik, "build-spss", "-i", genome_fa, "-o", work_dir / "x.spss", "--max-mem", "1G"])
    check(
        "build-spss --max-mem without --compact: rejected",
        proc.returncode == OPTIK_ERROR_EXIT and "--compact" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )


def main():
    """Run every spill check and exit non-zero on any failure."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--optik", default=str(PATHS["optik"]))
    parser.add_argument("--genome", default=str(DATASETS["ecoli_k31"].get("dbg_input")))
    parser.add_argument("--workdir", default="/tmp/optik_spill_check")
    args = parser.parse_args()
    optik = args.optik
    work_dir = Path(args.workdir)
    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True)
    genome_fa = work_dir / "genome.fa"
    opener = gzip.open if args.genome.endswith(".gz") else open
    with opener(args.genome, "rt") as src_file, open(genome_fa, "w") as dst_file:
        shutil.copyfileobj(src_file, dst_file)
    reads = work_dir / "reads.fa"
    subprocess.run(
        [
            sys.executable,
            str(HERE / "gen_reads.py"),
            "-r",
            str(genome_fa),
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
    spill_dir = work_dir / "spill"
    spill_dir.mkdir()

    cmds = {
        "dbruijn -c 1": lambda out, threads: [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "1",
            "-t",
            threads,
            "-o",
            out,
        ],
        "dbruijn -c 2 (k=21)": lambda out, threads: [
            optik,
            "dbruijn",
            "-i",
            reads,
            "-k",
            "21",
            "-c",
            "2",
            "-t",
            threads,
            "-o",
            out,
        ],
        "build-fm -c 2": lambda out, threads: [
            optik,
            "build-fm",
            "-i",
            reads,
            "-k",
            "31",
            "-c",
            "2",
            "-t",
            threads,
            "-o",
            out,
        ],
        "build-spss --compact -c 2": lambda out, threads: [
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
            threads,
            "-o",
            out,
        ],
    }
    # check_budget_path builds this in-memory graph, and check_several_inputs compares against it.
    base = work_dir / "base.gfa"
    fx = SpillInputs(optik=optik, work_dir=work_dir, genome_fa=genome_fa, reads=reads, spill_dir=spill_dir, cmds=cmds, base=base)
    for section in (check_forced_spill, check_compact_graph, check_partitioned_compaction, check_budget_path, check_default_spill_dir, check_table_over_limit, check_graph_over_limit, check_unusable_spill_dirs, check_failing_spill_write, check_several_inputs):
        section(fx)
    print(f"\nspill_check: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}")
    if not FAILS:
        shutil.rmtree(work_dir, ignore_errors=True)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
