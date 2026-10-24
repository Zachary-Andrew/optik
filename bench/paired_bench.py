#!/usr/bin/env python3
"""paired_bench.py - paired-end mapping, optik map -q R1 -q R2 against minimap2 -ax sr R1 R2 and against optik on each mate alone.

Reads come from bench/gen_pairs.py, an FR library with a truth table.
For each mapper the script reports
  - per mate, the mates placed correctly, meaning right reference and strand with the origin within 20 bp, the MAPQ >= 30 count and how many of those are wrong
  - MAPQ calibration, i.e. per MAPQ band the reads placed wrongly next to the expected number, the sum of 10^(-MAPQ/10)
  - pairs, i.e. the properly-paired flag 0x2 on pairs whose mates are both right, TLEN equal to the true fragment length minus soft-clipped outer read ends, within 2 for indels, and mate fields RNEXT, PNEXT, 0x20 and 0x8 consistent with the mate's record
  - downstream use, i.e. samtools fixmate -m, sort and markdup run, and, if Picard is available, ValidateSamFile STRICT reports no errors.

--gate fails on a mate-field inconsistency, a Picard error, an optik error
rate above 0.1% at MAPQ >= 30, an over-confident MAPQ band (wrong > 2x implied + 3), or paired optik placing fewer mates correctly than single-end optik.

Usage: python3 bench/paired_bench.py --ref REF.fa [--pairs N] [--error E]
           [--threads T] [--out DIR] [--optik BIN] [--minimap2 BIN] [--picard BIN] [--gate]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

import hw
import report as results_file
from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_MATE_REVERSE,
    SAM_FLAG_MATE_UNMAPPED,
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_READ1,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FAILS = []
CHECKS = []  # every check run, as rows of the results file


def check(name, ok, detail=""):
    """Print a check's verdict and keep the names of failed ones in FAILS."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    CHECKS.append({"check": name, "gate": bool(ok), "detail": detail})
    if not ok:
        FAILS.append(name)


def load_truth(tsv):
    """Map each pair name to (ref, pos1, strand1, pos2, strand2, fragment) from gen_pairs' truth TSV."""
    t = {}
    with open(tsv) as fh:
        for line in fh:
            f = line.split()
            if f and f[0] != "name":
                t[f[0]] = (f[1], int(f[2]), f[3], int(f[4]), f[5], int(f[6]))
    return t


def score(sam, truth, paired):
    """Per-mate and per-pair statistics of a SAM file."""
    prim = {}  # (name, mate) -> fields
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                continue
            name = f[0]
            if paired:
                mate = 1 if flag & SAM_FLAG_READ1 else 2
            else:
                name, mate = name.rsplit("/", 1)
                mate = int(mate)
            prim[(name, mate)] = f
    st = {
        "mates": 0,
        "correct": 0,
        "hq": 0,
        "hq_wrong": 0,
        "unmapped": 0,
        "bands": {},
        "proper_right": 0,
        "both_right": 0,
        "tlen_ok": 0,
        "mate_field_bad": 0,
    }
    for name, (ref, p1, s1, p2, s2, frag) in truth.items():
        ok_pair = []
        for mate, (p, s) in ((1, (p1, s1)), (2, (p2, s2))):
            f = prim.get((name, mate))
            st["mates"] += 1
            if f is None or int(f[1]) & SAM_FLAG_UNMAPPED:
                st["unmapped"] += 1
                ok_pair.append(False)
                continue
            flag, q = int(f[1]), int(f[4])
            ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                f[2] == ref
                and (("-" if flag & SAM_FLAG_REVERSE else "+") == s)
                and abs(int(f[3]) - 1 - lead - p) <= 20
            )
            st["correct"] += good
            ok_pair.append(good)
            if q >= HIGH_MAPQ:
                st["hq"] += 1
                st["hq_wrong"] += not good
            b = min(q // 10 * 10, 60)
            band = st["bands"].setdefault(b, [0, 0, 0.0])
            band[0] += 1
            band[1] += not good
            band[2] += 10 ** (-q / 10)
        if paired and all(ok_pair):
            st["both_right"] += 1
            f1, f2 = prim[(name, 1)], prim[(name, 2)]
            st["proper_right"] += bool(
                int(f1[1]) & SAM_FLAG_PROPER_PAIR and int(f2[1]) & SAM_FLAG_PROPER_PAIR
            )

            # The SAM spec has TLEN span the outermost aligned bases.
            # A soft-clipped outer read end, common in high-error reads, shortens it by the clip, so the expected value is the fragment minus the outer clips.
            def clips(f):
                """(leading, trailing) clip lengths of SAM record f."""
                ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
                lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                trail = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] in "SH" else 0
                return lead, trail

            fwd, rev = (f1, f2) if not int(f1[1]) & SAM_FLAG_REVERSE else (f2, f1)
            exp_tlen = frag - clips(fwd)[0] - clips(rev)[1]
            st["tlen_ok"] += abs(abs(int(f1[8])) - exp_tlen) <= 2 and int(
                f1[8]
            ) == -int(f2[8])
        if paired and (name, 1) in prim and (name, 2) in prim:
            for a, b in (((name, 1), (name, 2)), ((name, 2), (name, 1))):
                fa, fb = prim[a], prim[b]
                fla, flb = int(fa[1]), int(fb[1])
                bad = False
                if not flb & SAM_FLAG_UNMAPPED:
                    rnext = "=" if fb[2] == fa[2] else fb[2]
                    bad |= fa[6] != rnext or fa[7] != fb[3]
                    bad |= bool(fla & SAM_FLAG_MATE_REVERSE) != bool(
                        flb & SAM_FLAG_REVERSE
                    )
                bad |= bool(fla & SAM_FLAG_MATE_UNMAPPED) != bool(
                    flb & SAM_FLAG_UNMAPPED
                )
                st["mate_field_bad"] += bad
    return st


def run(cmd, out=None):
    with open(out, "w") if out else open(os.devnull, "w") as fh:
        p = subprocess.run(
            [str(c) for c in cmd],
            stdout=fh if out else subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if p.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{p.stderr[-2000:]}")
    return p


def main():
    """Simulate pairs, map them paired and single-end, and report or gate the statistics."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ref", required=True)
    ap.add_argument("--pairs", type=int, default=100000)
    ap.add_argument("--error", type=float, default=0.005)
    ap.add_argument("--threads", type=int, default=hw.cpus())
    ap.add_argument("--out", default=os.path.join(ROOT, "temp", "paired_bench"))
    ap.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    ap.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    ap.add_argument(
        "--picard", default=os.environ.get("PICARD_BIN") or shutil.which("picard")
    )
    ap.add_argument("--samtools", default=shutil.which("samtools"))
    ap.add_argument("--gate", action="store_true")
    results_file.add_argument(ap)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pre = os.path.join(a.out, f"pairs_e{a.error}_n{a.pairs}")
    if not os.path.exists(pre + ".tsv"):
        run(
            [
                sys.executable,
                os.path.join(HERE, "gen_pairs.py"),
                "-r",
                a.ref,
                "-n",
                a.pairs,
                "-e",
                a.error,
                "-o",
                pre,
            ]
        )
    r1, r2 = pre + "_1.fq", pre + "_2.fq"
    truth = load_truth(pre + ".tsv")
    threads = str(a.threads)
    sams = {}
    sams["optik paired"] = os.path.join(a.out, "optik_pe.sam")
    run(
        [
            a.optik,
            "map",
            "-r",
            a.ref,
            "-q",
            r1,
            "-q",
            r2,
            "-t",
            threads,
            "-R",
            "@RG\\tID:s1\\tSM:s1\\tPL:ILLUMINA",
            "-o",
            sams["optik paired"],
        ]
    )
    # Single-end runs map each mate file alone, with records renamed back to p<i>/<mate>.
    se = os.path.join(a.out, "optik_se.sam")
    with open(se, "w") as out:
        for m, rf in ((1, r1), (2, r2)):
            tmp = os.path.join(a.out, f"se{m}.sam")
            run([a.optik, "map", "-r", a.ref, "-q", rf, "-t", threads, "-o", tmp])
            with open(tmp) as fh:
                for line in fh:
                    if line[0] != "@":
                        f = line.split("\t", 1)
                        out.write(
                            f"{f[0]}/{m}\t{f[1]}"
                            if not f[0].endswith(f"/{m}")
                            else line
                        )
    sams["optik single-end"] = se
    if os.path.exists(a.minimap2):
        sams["minimap2 paired"] = os.path.join(a.out, "minimap2_pe.sam")
        run(
            [a.minimap2, "-ax", "sr", "-t", threads, a.ref, r1, r2],
            out=sams["minimap2 paired"],
        )
    res = {}
    for name, sam in sams.items():
        st = score(sam, truth, paired="paired" in name)
        res[name] = st
        bands = "  ".join(
            f"{b}:{w}/{e:.1f}"
            for b, (n, w, e) in sorted(st["bands"].items())
            if b >= 10
        )
        line = (
            f"{name:18} correct {st['correct']}/{st['mates']}  MAPQ>=30 {st['hq']} "
            f"wrong {st['hq_wrong']}  unmapped {st['unmapped']}"
        )
        if "paired" in name:
            line += (
                f"  proper|both-right {st['proper_right']}/{st['both_right']}  "
                f"TLEN right {st['tlen_ok']}  mate-field errors {st['mate_field_bad']}"
            )
        print(line)
        print(f"{'':18} wrong/implied by MAPQ band: {bands}")
    o, s1 = res["optik paired"], res["optik single-end"]
    print()
    check(
        "mate fields consistent (RNEXT, PNEXT, 0x20, 0x8 agree with the mate's record)",
        o["mate_field_bad"] == 0,
        f"{o['mate_field_bad']} inconsistent",
    )
    check(
        "pairing places at least as many mates correctly as single-end",
        o["correct"] >= s1["correct"],
        f"{o['correct']} vs {s1['correct']}",
    )
    check(
        "MAPQ >= 30 wrong <= 0.1%",
        o["hq_wrong"] <= max(1, o["hq"] // 1000),
        f"{o['hq_wrong']}/{o['hq']}",
    )
    over = [
        (b, w, round(e, 1))
        for b, (n, w, e) in sorted(o["bands"].items())
        if b >= 10 and w > 2 * e + 3
    ]
    check("no overconfident MAPQ band (wrong <= 2x implied + 3)", not over, str(over))
    check(
        "proper-pair flag on >= 99% of pairs with both mates right",
        o["proper_right"] >= 0.99 * o["both_right"],
        f"{o['proper_right']}/{o['both_right']}",
    )
    check(
        "TLEN = true fragment length on >= 99% of properly paired right pairs",
        o["tlen_ok"] >= 0.99 * o["proper_right"],
        f"{o['tlen_ok']}/{o['proper_right']}",
    )
    if a.samtools:
        fx = os.path.join(a.out, "fixmate.bam")
        st_ = os.path.join(a.out, "sorted.bam")
        md = os.path.join(a.out, "markdup.bam")
        p1 = subprocess.run(
            [a.samtools, "fixmate", "-m", sams["optik paired"], fx],
            capture_output=True,
            text=True,
            check=False,
        )
        p2 = subprocess.run(
            [a.samtools, "sort", "-o", st_, fx],
            capture_output=True,
            text=True,
            check=False,
        )
        p3 = subprocess.run(
            [a.samtools, "markdup", st_, md],
            capture_output=True,
            text=True,
            check=False,
        )
        check(
            "samtools fixmate -m | sort | markdup run on optik's paired SAM",
            p1.returncode == p2.returncode == p3.returncode == 0,
            (p1.stderr + p2.stderr + p3.stderr)[-200:],
        )
        if a.picard and not (os.path.isfile(a.picard) and os.access(a.picard, os.X_OK)):
            print(
                f"  [SKIP] Picard ValidateSamFile: {a.picard} is not an executable (PICARD_BIN stale?)"
            )
        elif a.picard:
            pv = subprocess.run(
                [
                    a.picard,
                    "ValidateSamFile",
                    "-I",
                    st_,
                    "-R",
                    a.ref,
                    "-MODE",
                    "SUMMARY",
                    "-IGNORE_WARNINGS",
                    "true",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            txt = pv.stdout + pv.stderr
            check(
                "Picard ValidateSamFile (STRICT): no errors on optik's paired BAM",
                "No errors found" in txt,
                "\n".join(
                    text_line for text_line in txt.splitlines() if "ERROR" in text_line
                )[:400],
            )
    print(f"\npaired_bench: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}")
    mappers = []
    for name, st in res.items():
        bands = "  ".join(
            f"{b}:{w}/{e:.1f}"
            for b, (n, w, e) in sorted(st["bands"].items())
            if b >= 10
        )
        row = {
            "mapper": name,
            "mates": st["mates"],
            "correct": st["correct"],
            "mapq30": st["hq"],
            "mapq30_wrong": st["hq_wrong"],
            "unmapped": st["unmapped"],
        }
        if "paired" in name:
            row.update(
                proper_of_both_right=f"{st['proper_right']}/{st['both_right']}",
                tlen_right=st["tlen_ok"],
                mate_field_errors=st["mate_field_bad"],
            )
        mappers.append({**row, "wrong/implied by MAPQ band": bands})
    results_file.write_section(
        a.report,
        "paired",
        f"Paired-end mapping on {a.pairs} gen_pairs.py pairs at error {a.error}: optik paired, "
        "optik per mate and minimap2 -ax sr, scored against the simulation truth",
        [("mappers", mappers), ("checks", CHECKS)],
        optik=a.optik,
        notes=[f"SAMs: `{os.path.relpath(a.out, ROOT)}`"],
    )
    if a.gate:
        sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
