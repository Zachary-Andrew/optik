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
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from edge_cases import kmers_of
from sam_spec import SAM_FLAG_UNMAPPED

FAILS = []


def check(name, ok, detail=""):
    """Print PASS or FAIL for one check and collect its name on failure."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def run(cmd, cwd):
    p = subprocess.run(
        cmd, capture_output=True, text=True, timeout=600, cwd=cwd, check=False
    )
    return p.returncode, p.stdout, p.stderr


def expect(name, cmd, want, cwd):
    """Run cmd and check that it ends as `want`, "ok" or "err", without crashing."""
    rc, out, err = run(cmd, cwd)
    crashed = rc < 0 or rc >= 128 or "panicked" in err or "terminate called" in err
    if want == "ok":
        good = rc == 0 and not crashed
    else:
        good = rc == 1 and "error" in err.lower() and not crashed
    check(name, good, "" if good else f"rc={rc} err={err.strip()[-160:]!r}")
    return rc, out, err


def counts(out):
    """The num_* counters of an optik query report."""
    return {
        m.group(1): int(m.group(2))
        for m in re.finditer(r"^(num_\w+) = (\d+)", out, re.MULTILINE)
    }


def main():
    """Run every boundary case and exit non-zero on any failure."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--optik", default=os.path.join(HERE, "..", "target/release/optik"))
    args = ap.parse_args()
    optik = os.path.abspath(args.optik)
    rng = random.Random(5)
    with tempfile.TemporaryDirectory(prefix="optik_bnd_") as td:
        ref = "".join(rng.choice("ACGT") for _ in range(4000))
        Path(os.path.join(td, "ref.fa")).write_text(f">a\n{ref}\n")
        with gzip.open(os.path.join(td, "ref.fa.gz"), "wt") as fh:
            fh.write(f">a\n{ref}\n")
        reads = [ref[i : i + 150] for i in range(0, 3800, 97)]
        Path(os.path.join(td, "reads.fa")).write_text(
            "".join(f">r{i}\n{r}\n" for i, r in enumerate(reads))
        )
        Path(os.path.join(td, "empty.fa")).write_text("")
        Path(os.path.join(td, "hdr.fa")).write_text(">x\n>y\n")
        Path(os.path.join(td, "alln.fa")).write_text(">n\n" + "N" * 500 + "\n")
        os.mkdir(os.path.join(td, "adir"))
        q = os.path.join(td, "q.fa")
        qs = [ref[100:400], ref[1000:1100].lower(), "ACGTN" * 20, "AC", ""]
        Path(q).write_text("".join(f">q{i}\n{s}\n" for i, s in enumerate(qs)))

        print("== k boundaries")
        for k in ("0", "33", "-1", "x", "1.5", ""):
            for sub in ("build-spss", "build-fm", "dbruijn"):
                expect(
                    f"{sub} -k {k!r} rejected",
                    [optik, sub, "-i", "ref.fa", "-k", k, "-o", "o"],
                    "err",
                    td,
                )
            expect(
                f"map -k {k!r} rejected",
                [
                    optik,
                    "map",
                    "-r",
                    "ref.fa",
                    "-q",
                    "reads.fa",
                    "-k",
                    k,
                    "-o",
                    "o.sam",
                ],
                "err",
                td,
            )
        expect(
            "dbruijn -k 32 rejected (edge set needs k+1 <= 32)",
            [optik, "dbruijn", "-i", "ref.fa", "-k", "32", "-o", "o"],
            "err",
            td,
        )
        expect(
            "build-spss -k 1 rejected (needs m < k)",
            [optik, "build-spss", "-i", "ref.fa", "-k", "1", "-o", "o"],
            "err",
            td,
        )

        print("== query exactness at the k limits")
        for k in (1, 2, 3, 31, 32):
            index = set(kmers_of(ref, k)[1])
            tot = val = pos = 0
            for s in qs:
                t, ks = kmers_of(s, k)
                tot += t
                val += len(ks)
                pos += sum(1 for x in ks if x in index)
            want = {
                "num_kmers": tot,
                "num_positive_kmers": pos,
                "num_negative_kmers": val - pos,
                "num_invalid_kmers": tot - val,
            }
            for sub in ("build-spss", "build-fm"):
                if sub == "build-spss" and k < 2:
                    continue
                idx = os.path.join(td, f"i_{sub}_{k}")
                rc, _, err = run(
                    [optik, sub, "-i", "ref.fa", "-k", str(k), "-o", idx], td
                )
                if rc != 0:
                    check(f"query k={k} {sub} build", False, err.strip()[-160:])
                    continue
                rc, out, err = run(
                    [optik, "query", "-i", idx, "-q", q, "-k", str(k)], td
                )
                got = counts(out)
                bad = {
                    key: (got.get(key), v)
                    for key, v in want.items()
                    if got.get(key) != v
                }
                check(
                    f"query k={k} {sub}", rc == 0 and not bad, str(bad) if bad else ""
                )

        print("== map parameters")
        for w in ("0", "1", "256"):
            expect(
                f"map -w {w}",
                [
                    optik,
                    "map",
                    "-r",
                    "ref.fa",
                    "-q",
                    "reads.fa",
                    "-w",
                    w,
                    "-o",
                    "o.sam",
                ],
                "ok",
                td,
            )
        for w in ("257", "-1", "ten"):
            expect(
                f"map -w {w} rejected",
                [
                    optik,
                    "map",
                    "-r",
                    "ref.fa",
                    "-q",
                    "reads.fa",
                    "-w",
                    w,
                    "-o",
                    "o.sam",
                ],
                "err",
                td,
            )
        for e in ("auto", "0", "5", "1000"):
            expect(
                f"map -e {e}",
                [
                    optik,
                    "map",
                    "-r",
                    "ref.fa",
                    "-q",
                    "reads.fa",
                    "-e",
                    e,
                    "-o",
                    "o.sam",
                ],
                "ok",
                td,
            )
        for e in ("-1", "x"):
            expect(
                f"map -e {e} rejected",
                [
                    optik,
                    "map",
                    "-r",
                    "ref.fa",
                    "-q",
                    "reads.fa",
                    "-e",
                    e,
                    "-o",
                    "o.sam",
                ],
                "err",
                td,
            )
        rc, _, _ = run(
            [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-w", "1", "-o", "w1.sam"],
            td,
        )
        with open(os.path.join(td, "w1.sam")) as fh:
            mapped = sum(
                1
                for line in fh
                if not line.startswith("@")
                and not int(line.split("\t")[1]) & SAM_FLAG_UNMAPPED
            )
        check(
            "map -w 1 maps every exact read",
            rc == 0 and mapped == len(reads),
            f"{mapped}/{len(reads)}",
        )

        print("== thread / count / m parameters")
        for t in ("0", "1", "4096"):
            expect(
                f"query -t {t}",
                [
                    optik,
                    "query",
                    "-i",
                    os.path.join(td, "i_build-spss_31"),
                    "-q",
                    q,
                    "-k",
                    "31",
                    "-t",
                    t,
                ],
                "ok",
                td,
            )
        for t in ("-1", "abc", "4097", "1.5"):
            expect(
                f"query -t {t} rejected",
                [
                    optik,
                    "query",
                    "-i",
                    os.path.join(td, "i_build-spss_31"),
                    "-q",
                    q,
                    "-k",
                    "31",
                    "-t",
                    t,
                ],
                "err",
                td,
            )
        for c in ("4294967296", "-2", "two"):
            expect(
                f"dbruijn -c {c} rejected",
                [optik, "dbruijn", "-i", "ref.fa", "-k", "21", "-c", c, "-o", "o.gfa"],
                "err",
                td,
            )
        rc, out, err = expect(
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
            td,
        )
        if rc == 0:
            check(
                "empty graph GFA is header-only",
                Path(os.path.join(td, "e.gfa")).read_text() == "H\tVN:Z:1.0\n",
            )
        for m in ("0", "1", "20"):
            expect(
                f"build-spss -k 21 -m {m}",
                [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-m", m, "-o", "o"],
                "ok",
                td,
            )
        for m in ("21", "22", "-3"):
            expect(
                f"build-spss -k 21 -m {m} rejected",
                [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-m", m, "-o", "o"],
                "err",
                td,
            )

        print("== degenerate inputs")
        for f in ("nonexist.fa", "adir"):
            for sub in ("build-spss", "build-fm", "dbruijn"):
                expect(
                    f"{sub} -i {f} rejected",
                    [optik, sub, "-i", f, "-k", "21", "-o", "o"],
                    "err",
                    td,
                )
            expect(
                f"map -r {f} rejected",
                [optik, "map", "-r", f, "-q", "reads.fa", "-o", "o.sam"],
                "err",
                td,
            )
            expect(
                f"map -q {f} rejected",
                [optik, "map", "-r", "ref.fa", "-q", f, "-o", "o.sam"],
                "err",
                td,
            )
            expect(
                f"query -q {f} rejected",
                [
                    optik,
                    "query",
                    "-i",
                    os.path.join(td, "i_build-spss_31"),
                    "-q",
                    f,
                    "-k",
                    "31",
                ],
                "err",
                td,
            )
        expect(
            "query -i <FASTA, not an index> rejected",
            [optik, "query", "-i", "ref.fa", "-q", q, "-k", "31"],
            "err",
            td,
        )
        for f in ("empty.fa", "hdr.fa", "alln.fa"):
            # An input without any k-mer is a user error for every index and graph builder alike, while k-mers that only miss the -c cutoff still give an empty graph (see "nothing solid" above).
            for sub in ("build-spss", "build-fm"):
                expect(
                    f"{sub} -i {f} rejected (no k-mers)",
                    [optik, sub, "-i", f, "-k", "21", "-o", "o"],
                    "err",
                    td,
                )
            expect(
                f"dbruijn -i {f} rejected (no k-mers)",
                [optik, "dbruijn", "-i", f, "-k", "21", "-o", "o.gfa"],
                "err",
                td,
            )
            expect(
                f"map -r {f} (nothing to map to)",
                [optik, "map", "-r", f, "-q", "reads.fa", "-o", "o.sam"],
                "ok",
                td,
            )
            rc, _, _ = run([optik, "map", "-r", "ref.fa", "-q", f, "-o", "e.sam"], td)
            with open(os.path.join(td, "e.sam")) as fh:
                n = sum(1 for line in fh if not line.startswith("@"))
            want_n = {"empty.fa": 0, "hdr.fa": 2, "alln.fa": 1}[f]
            check(
                f"map -q {f}: one record per input record",
                rc == 0 and n == want_n,
                f"{n} vs {want_n}",
            )
        for sub, flagset in (
            ("map", ["-q", "reads.fa", "-o", "o.sam"]),
            ("query", ["-q", q]),
            ("build-spss", ["-o", "o"]),
            ("dbruijn", ["-o", "o.gfa"]),
        ):
            expect(
                f"{sub} missing required input flag", [optik, sub] + flagset, "err", td
            )
        expect(
            "build-spss -o <unwritable>",
            [optik, "build-spss", "-i", "ref.fa", "-k", "21", "-o", "/proc/nope/o"],
            "err",
            td,
        )

        print("== gzip input equals plain input")
        for sub in ("build-spss", "dbruijn"):
            outs = []
            for f in ("ref.fa", "ref.fa.gz"):
                o = os.path.join(td, f"{sub}_{f}.out")
                run(
                    [optik, sub, "-i", f, "-k", "21", "-c", "1", "-o", o]
                    if sub == "dbruijn"
                    else [optik, sub, "-i", f, "-k", "21", "-o", o],
                    td,
                )
                outs.append(Path(o).read_bytes() if os.path.exists(o) else None)
            check(
                f"{sub}: .gz output identical to plain",
                outs[0] is not None and outs[0] == outs[1],
            )
        run([optik, "map", "-r", "ref.fa.gz", "-q", "reads.fa", "-o", "gz.sam"], td)
        run([optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-o", "pl.sam"], td)

        def body(p):
            """SAM lines of file p without its @ header lines."""
            with open(os.path.join(td, p)) as fh:
                return [line for line in fh if not line.startswith("@")]

        check(
            "map: gz reference gives identical SAM body",
            body("gz.sam") == body("pl.sam"),
        )

        print("== truncated / corrupt inputs are errors, never partial results")
        recs = "".join(f">r{i}\n{r}\n" for i, r in enumerate(reads * 20))
        gzdata = gzip.compress(recs.encode())
        Path(os.path.join(td, "trunc.fa.gz")).write_bytes(gzdata[: len(gzdata) // 3])
        idx31 = os.path.join(td, "i_build-spss_31")
        expect(
            "query -q <truncated .gz> rejected",
            [optik, "query", "-i", idx31, "-q", "trunc.fa.gz", "-k", "31"],
            "err",
            td,
        )
        rc, out, _ = expect(
            "map -q <truncated .gz> rejected",
            [optik, "map", "-r", "ref.fa", "-q", "trunc.fa.gz"],
            "err",
            td,
        )
        for sub in ("build-spss", "build-fm"):
            expect(
                f"{sub} -i <truncated .gz> rejected",
                [optik, sub, "-i", "trunc.fa.gz", "-k", "21", "-o", "o"],
                "err",
                td,
            )
        expect(
            "dbruijn -i <truncated .gz> rejected",
            [optik, "dbruijn", "-i", "trunc.fa.gz", "-k", "21", "-o", "o.gfa"],
            "err",
            td,
        )
        expect(
            "map -r <truncated .gz> rejected",
            [optik, "map", "-r", "trunc.fa.gz", "-q", "reads.fa", "-o", "o.sam"],
            "err",
            td,
        )
        Path(os.path.join(td, "it's.fa.gz")).write_bytes(
            gzip.compress(f">a\n{ref}\n".encode())
        )
        expect(
            "gzip path with a single quote",
            [optik, "build-spss", "-i", "it's.fa.gz", "-k", "21", "-o", "o"],
            "ok",
            td,
        )
        blob = Path(idx31).read_bytes()
        Path(os.path.join(td, "half.spss")).write_bytes(blob[: len(blob) // 2])
        expect(
            "query -i <truncated SPSS index> rejected",
            [optik, "query", "-i", "half.spss", "-q", q, "-k", "31"],
            "err",
            td,
        )
        for frac in (0.3, 0.6, 0.9):
            b = bytearray(blob)
            b[int(len(b) * frac)] ^= 0x5A
            Path(os.path.join(td, "flip.spss")).write_bytes(bytes(b))
            rc, _, err = expect(
                f"query -i <SPSS index with a flipped byte at {frac:.0%}> rejected",
                [optik, "query", "-i", "flip.spss", "-q", q, "-k", "31"],
                "err",
                td,
            )
            check(f"  ... reports a checksum mismatch ({frac:.0%})", "checksum" in err)
        expect(
            "query --no-verify skips the checksum (intact index)",
            [optik, "query", "-i", idx31, "-q", q, "-k", "31", "--no-verify"],
            "ok",
            td,
        )

        print("== outputs never overwrite inputs")
        import shutil

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
            shutil.copy(os.path.join(td, "ref.fa"), os.path.join(td, "cp.fa"))
            expect(f"{sub} rejected", cmd, "err", td)
            check(
                f"  ... input left intact ({sub})",
                Path(os.path.join(td, "cp.fa")).read_text()
                == Path(os.path.join(td, "ref.fa")).read_text(),
            )
        shutil.copy(idx31, os.path.join(td, "inpl.spss"))
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
            td,
        )
        expect(
            "  ... and the updated index loads",
            [optik, "query", "-i", "inpl.spss", "-q", q, "-k", "31"],
            "ok",
            td,
        )

        print("== malformed command lines")
        expect(
            "map unknown flag (-K) rejected",
            [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-K", "25"],
            "err",
            td,
        )
        expect(
            "query unknown flag (--thread) rejected",
            [optik, "query", "-i", idx31, "-q", q, "--thread", "4"],
            "err",
            td,
        )
        expect(
            "map flag missing its value rejected",
            [optik, "map", "-r", "ref.fa", "-q", "reads.fa", "-k"],
            "err",
            td,
        )
        expect(
            "spss-update --add <missing file> rejected",
            [optik, "spss-update", "-i", idx31, "--add", "nonexist.fa", "-o", "u.spss"],
            "err",
            td,
        )
        Path(os.path.join(td, "dupref.fa")).write_text(
            f">c\n{ref[:2000]}\n>c\n{ref[2000:]}\n"
        )
        expect(
            "map -r <duplicate sequence names> rejected (SAM @SQ must be unique)",
            [optik, "map", "-r", "dupref.fa", "-q", "reads.fa", "-o", "o.sam"],
            "err",
            td,
        )

    print(f"boundary_cases: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
