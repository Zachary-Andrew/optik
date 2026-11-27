#!/usr/bin/env python3
"""boundary_cases.py - explicit and implicit boundary values of the optik CLI.

For every active command, each parameter is driven to and past its limits and each input to its degenerate shapes.
The expectations are explicit.
  ok    exit 0, with exact output where an oracle exists
  err   exit 1 with an "error:" message, never a crash by signal or abort, a Rust panic or a silent fall-back to a default
A crash, meaning exit >= 128 or a kill by a signal, always fails.

Usage: python3 bench/boundary_cases.py [--optik target/release/optik]
"""

import argparse
import gzip
import os
import random
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from checklog import CheckLog
from diagnostics import DETAIL_CHARS
from edge_cases import kmers_of
from sam_spec import SamFlag

# seconds before one optik run counts as hung
RUN_TIMEOUT_S = 600
# a shell reports death by signal N as exit status SIGNAL_EXIT_BASE + N
SIGNAL_EXIT_BASE = 128

check = CheckLog()
FAILS = check.failed


def run(cmd, cwd):
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=RUN_TIMEOUT_S, cwd=cwd, check=False
    )
    return proc.returncode, proc.stdout, proc.stderr


def expect(name, cmd, want, cwd):
    """Run cmd and check that it ends as `want`, "ok" or "err", without crashing."""
    exit_code, out, err = run(cmd, cwd)
    crashed = exit_code < 0 or exit_code >= SIGNAL_EXIT_BASE or "panicked" in err or "terminate called" in err
    if want == "ok":
        good = exit_code == 0 and not crashed
    else:
        good = exit_code == 1 and "error" in err.lower() and not crashed
    check(name, good, "" if good else f"rc={exit_code} err={err.strip()[-DETAIL_CHARS:]!r}")
    return exit_code, out, err


def counts(out):
    """The num_* counters of an optik query report."""
    return {
        match.group(1): int(match.group(2))
        for match in re.finditer(r"^(num_\w+) = (\d+)", out, re.MULTILINE)
    }


@dataclass
class BoundaryInputs:
    """Holds the binary, the scratch directory and the small reference, reads and queries that every section checks against."""

    optik: str
    tmp_dir: str
    ref: str
    reads: list
    queries: list
    query_fa: str
    idx31: str


def check_k_boundaries(fx):
    """Rejects k values outside 1..32 and malformed ones on every command that takes -k."""
    optik, tmp_dir = fx.optik, fx.tmp_dir
    print("== k boundaries")
    # Every command that takes -k must reject values past 1..32 and anything that is not an integer.
    for kmer_arg in ("0", "33", "-1", "x", "1.5", ""):
        for sub in ("build-spss", "build-fm", "dbruijn"):
            expect(
                f"{sub} -k {kmer_arg!r} rejected",
                [optik, sub, "-i", "ref.fa", "-k", kmer_arg, "-o", "o"],
                "err",
                tmp_dir,
            )
        expect(
            f"map -k {kmer_arg!r} rejected",
            [
                optik,
                "map",
                "-r",
                "ref.fa",
                "-q",
                "reads.fa",
                "-k",
                kmer_arg,
                "-o",
                "o.sam",
            ],
            "err",
            tmp_dir,
        )
    # Two limits are tighter than the k range itself, because the graph stores (k+1)-mers and SPSS needs m < k.
    expect(
        "dbruijn -k 32 rejected (edge set needs k+1 <= 32)",
        [optik, "dbruijn", "-i", "ref.fa", "-k", "32", "-o", "o"],
        "err",
        tmp_dir,
    )
    expect(
        "build-spss -k 1 rejected (needs m < k)",
        [optik, "build-spss", "-i", "ref.fa", "-k", "1", "-o", "o"],
        "err",
        tmp_dir,
    )


def check_query_at_k_limits(fx):
    """Builds indexes at the smallest and largest k and checks that queries stay exact there."""
    optik, tmp_dir, ref, queries, query_fa = fx.optik, fx.tmp_dir, fx.ref, fx.queries, fx.query_fa
    print("== query exactness at the k limits")
    for kmer_arg in (1, 2, 3, 31, 32):
        index = set(kmers_of(ref, kmer_arg)[1])
        tot = val = pos = 0
        for seq in queries:
            total_kmers, valid_kmers = kmers_of(seq, kmer_arg)
            tot += total_kmers
            val += len(valid_kmers)
            pos += sum(1 for kmer in valid_kmers if kmer in index)
        want = {
            "num_kmers": tot,
            "num_positive_kmers": pos,
            "num_negative_kmers": val - pos,
            "num_invalid_kmers": tot - val,
        }
        for sub in ("build-spss", "build-fm"):
            if sub == "build-spss" and kmer_arg < 2:
                continue
            idx = os.path.join(tmp_dir, f"i_{sub}_{kmer_arg}")
            exit_code, _, err = run(
                [optik, sub, "-i", "ref.fa", "-k", str(kmer_arg), "-o", idx], tmp_dir
            )
            if exit_code != 0:
                check(f"query k={kmer_arg} {sub} build", False, err.strip()[-DETAIL_CHARS:])
                continue
            exit_code, out, err = run(
                [optik, "query", "-i", idx, "-q", query_fa, "-k", str(kmer_arg)], tmp_dir
            )
            got = counts(out)
            bad = {
                key: (got.get(key), expected)
                for key, expected in want.items()
                if got.get(key) != expected
            }
            check(
                f"query k={kmer_arg} {sub}", exit_code == 0 and not bad, str(bad) if bad else ""
            )


def check_map_parameters(fx):
    """Probes the accepted and rejected values of map's -w and -e."""
    optik, tmp_dir, reads = fx.optik, fx.tmp_dir, fx.reads
    print("== map parameters")
    # Each map parameter is tried at the edges of its range, which must work, and just past them, which must fail.
    for window_arg in ("0", "1", "256"):
        expect(
            f"map -w {window_arg}",
            [
                optik,
                "map",
                "-r",
                "ref.fa",
                "-q",
                "reads.fa",
                "-w",
                window_arg,
                "-o",
                "o.sam",
            ],
            "ok",
            tmp_dir,
        )
    for window_arg in ("257", "-1", "ten"):
        expect(
            f"map -w {window_arg} rejected",
            [
                optik,
                "map",
                "-r",
                "ref.fa",
                "-q",
                "reads.fa",
                "-w",
                window_arg,
                "-o",
                "o.sam",
            ],
            "err",
            tmp_dir,
        )
    for max_edits_arg in ("auto", "0", "5", "1000"):
        expect(
            f"map -e {max_edits_arg}",
            [
                optik,
                "map",
                "-r",
                "ref.fa",
                "-q",
                "reads.fa",
                "-e",
                max_edits_arg,
                "-o",
                "o.sam",
            ],
            "ok",
            tmp_dir,
        )
    for max_edits_arg in ("-1", "x"):
        expect(
            f"map -e {max_edits_arg} rejected",
            [
                optik,
                "map",
                "-r",
                "ref.fa",
                "-q",
                "reads.fa",
                "-e",
                max_edits_arg,
                "-o",
                "o.sam",
            ],
            "err",
            tmp_dir,
        )
    # The smallest window keeps every minimizer, so every exact read must still map.
    exit_code, _, _ = run(
        [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-w", "1", "-o", "w1.sam"],
        tmp_dir,
    )
    with open(os.path.join(tmp_dir, "w1.sam")) as handle:
        mapped = sum(
            1
            for line in handle
            if not line.startswith("@")
            and not int(line.split("\t")[1]) & SamFlag.UNMAPPED
        )
    check(
        "map -w 1 maps every exact read",
        exit_code == 0 and mapped == len(reads),
        f"{mapped}/{len(reads)}",
    )


def check_thread_count_and_m(fx):
    """Probes -t, -c and -m at and past their limits."""
    optik, tmp_dir, query_fa = fx.optik, fx.tmp_dir, fx.query_fa
    print("== thread / count / m parameters")
    # The thread count accepts 0 for automatic and up to 4096, and rejects anything else.
    for total_kmers in ("0", "1", "4096"):
        expect(
            f"query -t {total_kmers}",
            [
                optik,
                "query",
                "-i",
                os.path.join(tmp_dir, "i_build-spss_31"),
                "-q",
                query_fa,
                "-k",
                "31",
                "-t",
                total_kmers,
            ],
            "ok",
            tmp_dir,
        )
    for total_kmers in ("-1", "abc", "4097", "1.5"):
        expect(
            f"query -t {total_kmers} rejected",
            [
                optik,
                "query",
                "-i",
                os.path.join(tmp_dir, "i_build-spss_31"),
                "-q",
                query_fa,
                "-k",
                "31",
                "-t",
                total_kmers,
            ],
            "err",
            tmp_dir,
        )
    for count_arg in ("4294967296", "-2", "two"):
        expect(
            f"dbruijn -c {count_arg} rejected",
            [optik, "dbruijn", "-i", "ref.fa", "-k", "21", "-c", count_arg, "-o", "o.gfa"],
            "err",
            tmp_dir,
        )
    # A cutoff no k-mer reaches is not an error, so the result is a GFA holding only its header.
    exit_code, _out, _err = expect(
        "dbruijn -c 1000000 (nothing solid) -> empty graph",
        [
            optik,
            "dbruijn",
            "-i",
            "ref.fa",
            "-k",
            "21",
            "-c",
            "1000000",
            "-o",
            "e.gfa",
        ],
        "ok",
        tmp_dir,
    )
    if exit_code == 0:
        check(
            "empty graph GFA is header-only",
            Path(os.path.join(tmp_dir, "e.gfa")).read_text() == "H\tVN:Z:1.0\n",
        )
    for minimizer_arg in ("0", "1", "20"):
        expect(
            f"build-spss -k 21 -m {minimizer_arg}",
            [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-m", minimizer_arg, "-o", "o"],
            "ok",
            tmp_dir,
        )
    for minimizer_arg in ("21", "22", "-3"):
        expect(
            f"build-spss -k 21 -m {minimizer_arg} rejected",
            [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-m", minimizer_arg, "-o", "o"],
            "err",
            tmp_dir,
        )


def check_degenerate_inputs(fx):
    """Feeds missing files, directories, empty files, bare headers and all-N records."""
    optik, tmp_dir, query_fa = fx.optik, fx.tmp_dir, fx.query_fa
    print("== degenerate inputs")
    # Missing files and directories must give a clean error from every command, whichever flag names them.
    for input_name in ("nonexist.fa", "adir"):
        for sub in ("build-spss", "build-fm", "dbruijn"):
            expect(
                f"{sub} -i {input_name} rejected",
                [optik, sub, "-i", input_name, "-k", "21", "-o", "o"],
                "err",
                tmp_dir,
            )
        expect(
            f"map -r {input_name} rejected",
            [optik, "map", "-r", input_name, "-q", "reads.fa", "-o", "o.sam"],
            "err",
            tmp_dir,
        )
        expect(
            f"map -q {input_name} rejected",
            [optik, "map", "-r", "ref.fa", "-q", input_name, "-o", "o.sam"],
            "err",
            tmp_dir,
        )
        expect(
            f"query -q {input_name} rejected",
            [
                optik,
                "query",
                "-i",
                os.path.join(tmp_dir, "i_build-spss_31"),
                "-q",
                input_name,
                "-k",
                "31",
            ],
            "err",
            tmp_dir,
        )
    expect(
        "query -i <FASTA, not an index> rejected",
        [optik, "query", "-i", "ref.fa", "-q", query_fa, "-k", "31"],
        "err",
        tmp_dir,
    )
    # Inputs with no usable k-mer are an error for the index builders, while map still writes one record per read.
    for input_name in ("empty.fa", "hdr.fa", "alln.fa"):
        # An input without any k-mer is a user error for every index and graph builder alike, while k-mers that only miss the -c cutoff still give an empty graph (see "nothing solid" above).
        for sub in ("build-spss", "build-fm"):
            expect(
                f"{sub} -i {input_name} rejected (no k-mers)",
                [optik, sub, "-i", input_name, "-k", "21", "-o", "o"],
                "err",
                tmp_dir,
            )
        expect(
            f"dbruijn -i {input_name} rejected (no k-mers)",
            [optik, "dbruijn", "-i", input_name, "-k", "21", "-o", "o.gfa"],
            "err",
            tmp_dir,
        )
        expect(
            f"map -r {input_name} (nothing to map to)",
            [optik, "map", "-r", input_name, "-q", "reads.fa", "-o", "o.sam"],
            "ok",
            tmp_dir,
        )
        exit_code, _, _ = run([optik, "map", "-r", "ref.fa", "-q", input_name, "-o", "e.sam"], tmp_dir)
        with open(os.path.join(tmp_dir, "e.sam")) as handle:
            num_records = sum(1 for line in handle if not line.startswith("@"))
        want_n = {"empty.fa": 0, "hdr.fa": 2, "alln.fa": 1}[input_name]
        check(
            f"map -q {input_name}: one record per input record",
            exit_code == 0 and num_records == want_n,
            f"{num_records} vs {want_n}",
        )
    # Leaving out a required input flag must fail with a message instead of reading something by default.
    for sub, flagset in (
        ("map", ["-q", "reads.fa", "-o", "o.sam"]),
        ("query", ["-q", query_fa]),
        ("build-spss", ["-o", "o"]),
        ("dbruijn", ["-o", "o.gfa"]),
    ):
        expect(
            f"{sub} missing required input flag", [optik, sub] + flagset, "err", tmp_dir
        )
    expect(
        "build-spss -o <unwritable>",
        [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-o", "/proc/nope/o"],
        "err",
        tmp_dir,
    )
    # A record name that is not valid UTF-8 is refused by name instead of being altered, since SAM allows only printable ASCII.
    latin1_name = b"caf\xe9"
    # The error escapes non-ASCII bytes as \xNN, which backslashreplace reproduces.
    escaped_name = latin1_name.decode("ascii", "backslashreplace")
    first_read = Path(os.path.join(tmp_dir, "reads.fa")).read_bytes().split(b"\n")[1]
    Path(os.path.join(tmp_dir, "latin1.fa")).write_bytes(b">" + latin1_name + b"\n" + first_read + b"\n")

    for role, flagset in (("-r", ["-r", "latin1.fa", "-q", "reads.fa"]), ("-q", ["-r", "ref.fa", "-q", "latin1.fa"])):
        _exit_code, _out, err = expect(
            f"map {role} <non-UTF-8 record name> rejected",
            [optik, "map", *flagset, "-o", "o.sam"],
            "err",
            tmp_dir,
        )
        check("  ... names the record's exact bytes", escaped_name in err, err.strip()[-DETAIL_CHARS:])


def check_gzip_matches_plain(fx):
    """Confirms that gzip input gives exactly the same results as the plain file."""
    optik, tmp_dir = fx.optik, fx.tmp_dir
    print("== gzip input equals plain input")
    for sub in ("build-spss", "dbruijn"):
        outs = []
        for input_name in ("ref.fa", "ref.fa.gz"):
            out_path = os.path.join(tmp_dir, f"{sub}_{input_name}.out")
            run(
                [optik, sub, "-i", input_name, "-k", "21", "-c", "1", "-o", out_path]
                if sub == "dbruijn"
                else [optik, sub, "-i", input_name, "-k", "21", "-o", out_path],
                tmp_dir,
            )
            outs.append(Path(out_path).read_bytes() if os.path.exists(out_path) else None)
        check(
            f"{sub}: .gz output identical to plain",
            outs[0] is not None and outs[0] == outs[1],
        )
    run([optik, "map", "-r", "ref.fa.gz", "-q", "reads.fa", "-o", "gz.sam"], tmp_dir)
    run([optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-o", "pl.sam"], tmp_dir)

    def body(sam_name):
        """SAM lines of file p without its @ header lines."""
        with open(os.path.join(tmp_dir, sam_name)) as handle:
            return [line for line in handle if not line.startswith("@")]

    check(
        "map: gz reference gives identical SAM body",
        body("gz.sam") == body("pl.sam"),
    )


def check_truncated_and_corrupt_inputs(fx):
    """Makes sure truncated or corrupt inputs fail outright instead of giving partial results."""
    optik, tmp_dir, ref, reads, query_fa, idx31 = fx.optik, fx.tmp_dir, fx.ref, fx.reads, fx.query_fa, fx.idx31
    print("== truncated / corrupt inputs are errors, never partial results")
    # A gzip file cut at a third of its length must fail every reader, rather than give results for the part that survived.
    recs = "".join(f">r{i}\n{read}\n" for i, read in enumerate(reads * 20))
    gzdata = gzip.compress(recs.encode())
    Path(os.path.join(tmp_dir, "trunc.fa.gz")).write_bytes(gzdata[: len(gzdata) // 3])
    expect(
        "query -q <truncated .gz> rejected",
        [optik, "query", "-i", idx31, "-q", "trunc.fa.gz", "-k", "31"],
        "err",
        tmp_dir,
    )
    _exit_code, _out, _ = expect(
        "map -q <truncated .gz> rejected",
        [optik, "map", "-r", "ref.fa", "-q", "trunc.fa.gz"],
        "err",
        tmp_dir,
    )
    for sub in ("build-spss", "build-fm"):
        expect(
            f"{sub} -i <truncated .gz> rejected",
            [optik, sub, "-i", "trunc.fa.gz", "-k", "21", "-o", "o"],
            "err",
            tmp_dir,
        )
    expect(
        "dbruijn -i <truncated .gz> rejected",
        [optik, "dbruijn", "-i", "trunc.fa.gz", "-k", "21", "-o", "o.gfa"],
        "err",
        tmp_dir,
    )
    expect(
        "map -r <truncated .gz> rejected",
        [optik, "map", "-r", "trunc.fa.gz", "-q", "reads.fa", "-o", "o.sam"],
        "err",
        tmp_dir,
    )
    Path(os.path.join(tmp_dir, "it's.fa.gz")).write_bytes(
        gzip.compress(f">a\n{ref}\n".encode())
    )
    expect(
        "gzip path with a single quote",
        [optik, "build-spss", "-i", "it's.fa.gz", "-k", "21", "-o", "o"],
        "ok",
        tmp_dir,
    )
    # A truncated or bit-flipped index must fail its checksum, unless --no-verify skips it on purpose.
    blob = Path(idx31).read_bytes()
    Path(os.path.join(tmp_dir, "half.spss")).write_bytes(blob[: len(blob) // 2])
    expect(
        "query -i <truncated SPSS index> rejected",
        [optik, "query", "-i", "half.spss", "-q", query_fa, "-k", "31"],
        "err",
        tmp_dir,
    )
    for frac in (0.3, 0.6, 0.9):
        corrupted = bytearray(blob)
        corrupted[int(len(corrupted) * frac)] ^= 0x5A
        Path(os.path.join(tmp_dir, "flip.spss")).write_bytes(bytes(corrupted))
        _exit_code, _, err = expect(
            f"query -i <SPSS index with a flipped byte at {frac:.0%}> rejected",
            [optik, "query", "-i", "flip.spss", "-q", query_fa, "-k", "31"],
            "err",
            tmp_dir,
        )
        check(f"  ... reports a checksum mismatch ({frac:.0%})", "checksum" in err)
    expect(
        "query --no-verify skips the checksum (intact index)",
        [optik, "query", "-i", idx31, "-q", query_fa, "-k", "31", "--no-verify"],
        "ok",
        tmp_dir,
    )


def check_outputs_never_overwrite_inputs(fx):
    """Refuses outputs that would overwrite an input, while in-place index updates stay safe."""
    optik, tmp_dir, query_fa, idx31 = fx.optik, fx.tmp_dir, fx.query_fa, fx.idx31
    print("== outputs never overwrite inputs")
    import shutil

    # Each command is pointed at an output path that is also one of its inputs, which must be refused before anything is truncated.
    for sub, cmd in (
        (
            "map -o <reference>",
            [optik, "map", "-r", "cp.fa", "-q", "reads.fa", "-o", "cp.fa"],
        ),
        (
            "map -o <reads>",
            [optik, "map", "-r", "ref.fa", "-q", "cp.fa", "-o", "cp.fa"],
        ),
        (
            "build-spss -o <input>",
            [optik, "build-spss", "-i", "cp.fa", "-k", "21", "-o", "cp.fa"],
        ),
        (
            "dbruijn -o <input>",
            [optik, "dbruijn", "-i", "cp.fa", "-k", "21", "-o", "cp.fa"],
        ),
        (
            "query --ids <queries>",
            [
                optik,
                "query",
                "-i",
                idx31,
                "-q",
                "cp.fa",
                "-k",
                "31",
                "--ids",
                "cp.fa",
            ],
        ),
    ):
        shutil.copy(os.path.join(tmp_dir, "ref.fa"), os.path.join(tmp_dir, "cp.fa"))
        expect(f"{sub} rejected", cmd, "err", tmp_dir)
        check(
            f"  ... input left intact ({sub})",
            Path(os.path.join(tmp_dir, "cp.fa")).read_text()
            == Path(os.path.join(tmp_dir, "ref.fa")).read_text(),
        )
    # spss-update is the one exception, since it loads the whole index before rewriting it in place.
    shutil.copy(idx31, os.path.join(tmp_dir, "inpl.spss"))
    expect(
        "spss-update -i X -o X (in place) allowed",
        [
            optik,
            "spss-update",
            "-i",
            "inpl.spss",
            "--add",
            "reads.fa",
            "-o",
            "inpl.spss",
        ],
        "ok",
        tmp_dir,
    )
    expect(
        "  ... and the updated index loads",
        [optik, "query", "-i", "inpl.spss", "-q", query_fa, "-k", "31"],
        "ok",
        tmp_dir,
    )


def check_malformed_command_lines(fx):
    """Rejects malformed command lines with a clear error."""
    optik, tmp_dir, ref, query_fa, idx31 = fx.optik, fx.tmp_dir, fx.ref, fx.query_fa, fx.idx31
    print("== malformed command lines")
    expect(
        "map unknown flag (-K) rejected",
        [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-K", "25"],
        "err",
        tmp_dir,
    )
    expect(
        "query unknown flag (--thread) rejected",
        [optik, "query", "-i", idx31, "-q", query_fa, "--thread", "4"],
        "err",
        tmp_dir,
    )
    expect(
        "map flag missing its value rejected",
        [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-k"],
        "err",
        tmp_dir,
    )
    expect(
        "spss-update --add <missing file> rejected",
        [optik, "spss-update", "-i", idx31, "--add", "nonexist.fa", "-o", "u.spss"],
        "err",
        tmp_dir,
    )
    Path(os.path.join(tmp_dir, "dupref.fa")).write_text(
        f">c\n{ref[:2000]}\n>c\n{ref[2000:]}\n"
    )
    expect(
        "map -r <duplicate sequence names> rejected (SAM @SQ must be unique)",
        [optik, "map", "-r", "dupref.fa", "-q", "reads.fa", "-o", "o.sam"],
        "err",
        tmp_dir,
    )


def main():
    """Run every boundary case and exit non-zero on any failure."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    args = parser.parse_args()
    optik = os.path.abspath(args.optik)
    rng = random.Random(5)
    with tempfile.TemporaryDirectory(prefix="optik_bnd_") as tmp_dir:
        ref = "".join(rng.choice("ACGT") for _ in range(4000))
        Path(os.path.join(tmp_dir, "ref.fa")).write_text(f">a\n{ref}\n")
        with gzip.open(os.path.join(tmp_dir, "ref.fa.gz"), "wt") as handle:
            handle.write(f">a\n{ref}\n")
        reads = [ref[i : i + 150] for i in range(0, 3800, 97)]
        Path(os.path.join(tmp_dir, "reads.fa")).write_text(
            "".join(f">r{i}\n{read}\n" for i, read in enumerate(reads))
        )
        Path(os.path.join(tmp_dir, "empty.fa")).write_text("")
        Path(os.path.join(tmp_dir, "hdr.fa")).write_text(">x\n>y\n")
        Path(os.path.join(tmp_dir, "alln.fa")).write_text(">n\n" + "N" * 500 + "\n")
        os.mkdir(os.path.join(tmp_dir, "adir"))
        query_fa = os.path.join(tmp_dir, "q.fa")
        queries = [ref[100:400], ref[1000:1100].lower(), "ACGTN" * 20, "AC", ""]
        Path(query_fa).write_text("".join(f">q{i}\n{seq}\n" for i, seq in enumerate(queries)))

        # check_query_at_k_limits builds this k = 31 SPSS index, and later sections reuse it.
        idx31 = os.path.join(tmp_dir, "i_build-spss_31")
        fx = BoundaryInputs(optik=optik, tmp_dir=tmp_dir, ref=ref, reads=reads, queries=queries, query_fa=query_fa, idx31=idx31)
        for section in (check_k_boundaries, check_query_at_k_limits, check_map_parameters, check_thread_count_and_m, check_degenerate_inputs, check_gzip_matches_plain, check_truncated_and_corrupt_inputs, check_outputs_never_overwrite_inputs, check_malformed_command_lines):
            section(fx)

    print(f"boundary_cases: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
