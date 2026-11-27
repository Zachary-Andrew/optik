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
from diagnostics import DETAIL_CHARS, LOG_TAIL_CHARS
from mapq_calibration import overconfident
from sam_spec import (
    HIGH_MAPQ,
    PHRED_SCALE,
    SAM_MAX_MAPQ,
    SamCol,
    SamFlag,
)

# a mate is placed correctly when its origin is within this many bp of the truth
ORIGIN_TOL = 20
# TLEN may differ from the expected span by this many bp, for indels
TLEN_TOL = 2
# width of a MAPQ calibration band
MAPQ_BAND_WIDTH = 10
# at most one in this many MAPQ >= 30 mates may be wrong (0.1%)
HQ_WRONG_RATIO = 1000
# share of right pairs that must have the proper-pair flag set, and of those the right TLEN
MIN_PROPER_SHARE = 0.99
# argparse defaults
DEFAULT_PAIRS = 100000
DEFAULT_ERROR = 0.005

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FAILS = []
CHECKS = []  # every check run, as rows of the results file


# Mates are numbered as in their /1 and /2 name suffixes.
MATE1, MATE2 = 1, 2
# A MAPQ band counts its placements, its wrong ones and the errors its MAPQs imply, in this order.
BAND_COUNT, BAND_WRONG, BAND_IMPLIED = 0, 1, 2


def check(name, passed, detail=""):
    """Print a check's verdict and keep the names of failed ones in FAILS."""
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    CHECKS.append({"check": name, "gate": bool(passed), "detail": detail})
    if not passed:
        FAILS.append(name)


def load_truth(tsv):
    """Map each pair name to (ref, pos1, strand1, pos2, strand2, fragment) from gen_pairs' truth TSV."""
    truth = {}
    with open(tsv) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0] == "name":  # blank line or header
                continue
            name, ref, pos1, strand1, pos2, strand2, fragment, *_ = fields
            truth[name] = (ref, int(pos1), strand1, int(pos2), strand2, int(fragment))
    return truth


def score(sam, truth, paired):
    """Per-mate and per-pair statistics of a SAM file."""
    prim = {}  # (name, mate) -> fields
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            rec = line.rstrip("\n").split("\t")
            flag = int(rec[SamCol.FLAG])
            if flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                continue
            name = rec[SamCol.QNAME]
            if paired:
                mate = MATE1 if flag & SamFlag.READ1 else MATE2
            else:
                name, mate = name.rsplit("/", 1)
                mate = int(mate)
            prim[(name, mate)] = rec
    # Most counts are per mate, while proper_right and tlen_ok only count pairs whose two mates are both placed correctly.
    stats = {
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
    for name, (ref, pos1, strand1, pos2, strand2, frag) in truth.items():
        ok_pair = []
        for mate, (true_pos, true_strand) in ((MATE1, (pos1, strand1)), (MATE2, (pos2, strand2))):
            rec = prim.get((name, mate))
            stats["mates"] += 1
            if rec is None or int(rec[SamCol.FLAG]) & SamFlag.UNMAPPED:
                stats["unmapped"] += 1
                ok_pair.append(False)
                continue
            flag, mapq = int(rec[SamCol.FLAG]), int(rec[SamCol.MAPQ])
            ops = re.findall(r"(\d+)([MIDNSHP=X])", rec[SamCol.CIGAR])
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = (
                rec[SamCol.RNAME] == ref
                and (("-" if flag & SamFlag.REVERSE else "+") == true_strand)
                and abs(int(rec[SamCol.POS]) - 1 - lead - true_pos) <= ORIGIN_TOL
            )
            stats["correct"] += good
            ok_pair.append(good)
            if mapq >= HIGH_MAPQ:
                stats["hq"] += 1
                stats["hq_wrong"] += not good
            band_key = min(mapq // MAPQ_BAND_WIDTH * MAPQ_BAND_WIDTH, SAM_MAX_MAPQ)
            band = stats["bands"].setdefault(band_key, [0, 0, 0.0])
            band[BAND_COUNT] += 1
            band[BAND_WRONG] += not good
            band[BAND_IMPLIED] += PHRED_SCALE ** (-mapq / PHRED_SCALE)
        if paired and all(ok_pair):
            stats["both_right"] += 1
            rec1, rec2 = prim[(name, MATE1)], prim[(name, MATE2)]
            stats["proper_right"] += bool(
                int(rec1[SamCol.FLAG]) & SamFlag.PROPER_PAIR and int(rec2[SamCol.FLAG]) & SamFlag.PROPER_PAIR
            )

            # The SAM spec has TLEN span the outermost aligned bases.
            # A soft-clipped outer read end, common in high-error reads, shortens it by the clip, so the expected value is the fragment minus the outer clips.
            def clips(rec):
                """(leading, trailing) clip lengths of SAM record f."""
                ops = re.findall(r"(\d+)([MIDNSHP=X])", rec[SamCol.CIGAR])
                lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
                trail = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] in "SH" else 0
                return lead, trail

            fwd, rev = (rec1, rec2) if not int(rec1[SamCol.FLAG]) & SamFlag.REVERSE else (rec2, rec1)
            # Clipped outer bases shorten the aligned template, so TLEN is judged against the fragment minus those clips.
            exp_tlen = frag - clips(fwd)[0] - clips(rev)[1]
            stats["tlen_ok"] += abs(abs(int(rec1[SamCol.TLEN])) - exp_tlen) <= TLEN_TOL and int(
                rec1[SamCol.TLEN]
            ) == -int(rec2[SamCol.TLEN])
        if paired and (name, MATE1) in prim and (name, MATE2) in prim:
            for self_key, mate_key in (((name, MATE1), (name, MATE2)), ((name, MATE2), (name, MATE1))):
                self_rec, mate_rec = prim[self_key], prim[mate_key]
                fla, flb = int(self_rec[SamCol.FLAG]), int(mate_rec[SamCol.FLAG])
                bad = False
                # Each record's mate fields must mirror its mate's own placement and strand.
                if not flb & SamFlag.UNMAPPED:
                    rnext = "=" if mate_rec[SamCol.RNAME] == self_rec[SamCol.RNAME] else mate_rec[SamCol.RNAME]
                    bad |= self_rec[SamCol.RNEXT] != rnext or self_rec[SamCol.PNEXT] != mate_rec[SamCol.POS]
                    bad |= bool(fla & SamFlag.MATE_REVERSE) != bool(
                        flb & SamFlag.REVERSE
                    )
                bad |= bool(fla & SamFlag.MATE_UNMAPPED) != bool(
                    flb & SamFlag.UNMAPPED
                )
                stats["mate_field_bad"] += bad
    return stats


def run(cmd, out=None):
    with open(out, "w") if out else open(os.devnull, "w") as handle:
        proc = subprocess.run(
            [str(arg) for arg in cmd],
            stdout=handle if out else subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if proc.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{proc.stderr[-LOG_TAIL_CHARS:]}")
    return proc


def main():
    """Simulate pairs, map them paired and single-end, and report or gate the statistics."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--ref", required=True)
    parser.add_argument("--pairs", type=int, default=DEFAULT_PAIRS)
    parser.add_argument("--error", type=float, default=DEFAULT_ERROR)
    parser.add_argument("--threads", type=int, default=hw.cpus())
    parser.add_argument("--out", default=os.path.join(ROOT, "temp", "paired_bench"))
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    parser.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    parser.add_argument(
        "--picard", default=os.environ.get("PICARD_BIN") or shutil.which("picard")
    )
    parser.add_argument("--samtools", default=shutil.which("samtools"))
    parser.add_argument("--gate", action="store_true")
    results_file.add_argument(parser)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    # Simulated pairs are cached by error rate and pair count, so a rerun maps the same pairs.
    pre = os.path.join(args.out, f"pairs_e{args.error}_n{args.pairs}")
    if not os.path.exists(pre + ".tsv"):
        run(
            [
                sys.executable,
                os.path.join(HERE, "gen_pairs.py"),
                "-r",
                args.ref,
                "-n",
                args.pairs,
                "-e",
                args.error,
                "-o",
                pre,
            ]
        )
    mate1_path, mate2_path = pre + "_1.fq", pre + "_2.fq"
    truth = load_truth(pre + ".tsv")
    threads = str(args.threads)
    sams = {}
    sams["optik paired"] = os.path.join(args.out, "optik_pe.sam")
    run(
        [
            args.optik,
            "map",
            "-r",
            args.ref,
            "-q",
            mate1_path,
            "-q",
            mate2_path,
            "-t",
            threads,
            "-R",
            "@RG\\tID:s1\\tSM:s1\\tPL:ILLUMINA",
            "-o",
            sams["optik paired"],
        ]
    )
    # Single-end runs map each mate file alone, with records renamed back to p<i>/<mate>.
    se_sam = os.path.join(args.out, "optik_se.sam")
    with open(se_sam, "w") as out:
        for mate, mate_path in ((MATE1, mate1_path), (MATE2, mate2_path)):
            tmp = os.path.join(args.out, f"se{mate}.sam")
            run([args.optik, "map", "-r", args.ref, "-q", mate_path, "-t", threads, "-o", tmp])
            with open(tmp) as handle:
                for line in handle:
                    if line[0] != "@":
                        fields = line.split("\t", 1)
                        out.write(
                            f"{fields[0]}/{mate}\t{fields[1]}"
                            if not fields[0].endswith(f"/{mate}")
                            else line
                        )
    sams["optik single-end"] = se_sam
    # minimap2 is an optional reference row, so its absence only drops that row from the table.
    if os.path.exists(args.minimap2):
        sams["minimap2 paired"] = os.path.join(args.out, "minimap2_pe.sam")
        run(
            [args.minimap2, "-ax", "sr", "-t", threads, args.ref, mate1_path, mate2_path],
            out=sams["minimap2 paired"],
        )
    res = {}
    for name, sam in sams.items():
        stats = score(sam, truth, paired="paired" in name)
        res[name] = stats
        # Each MAPQ band shows its wrong placements next to the count its MAPQs imply, and the MAPQ 0 band is left out because it makes no claim.
        bands = "  ".join(
            f"{band_key}:{wrong}/{implied:.1f}"
            for band_key, (count, wrong, implied) in sorted(stats["bands"].items())
            if band_key >= MAPQ_BAND_WIDTH
        )
        line = (
            f"{name:18} correct {stats['correct']}/{stats['mates']}  MAPQ>=30 {stats['hq']} "
            f"wrong {stats['hq_wrong']}  unmapped {stats['unmapped']}"
        )
        if "paired" in name:
            line += (
                f"  proper|both-right {stats['proper_right']}/{stats['both_right']}  "
                f"TLEN right {stats['tlen_ok']}  mate-field errors {stats['mate_field_bad']}"
            )
        print(line)
        print(f"{'':18} wrong/implied by MAPQ band: {bands}")
    # Pairing has to place at least as many mates as mapping each file alone, with calibrated MAPQ and correct pair fields.
    optik_pe, optik_se = res["optik paired"], res["optik single-end"]
    print()
    check(
        "mate fields consistent (RNEXT, PNEXT, 0x20, 0x8 agree with the mate's record)",
        optik_pe["mate_field_bad"] == 0,
        f"{optik_pe['mate_field_bad']} inconsistent",
    )
    check(
        "pairing places at least as many mates correctly as single-end",
        optik_pe["correct"] >= optik_se["correct"],
        f"{optik_pe['correct']} vs {optik_se['correct']}",
    )
    check(
        "MAPQ >= 30 wrong <= 0.1%",
        optik_pe["hq_wrong"] <= max(1, optik_pe["hq"] // HQ_WRONG_RATIO),
        f"{optik_pe['hq_wrong']}/{optik_pe['hq']}",
    )
    over = [
        (band_key, wrong, round(implied, 1))
        for band_key, (count, wrong, implied) in sorted(optik_pe["bands"].items())
        if band_key >= MAPQ_BAND_WIDTH and overconfident(wrong, implied)
    ]
    check("no overconfident MAPQ band (wrong <= 2x implied + 3)", not over, str(over))
    check(
        "proper-pair flag on >= 99% of pairs with both mates right",
        optik_pe["proper_right"] >= MIN_PROPER_SHARE * optik_pe["both_right"],
        f"{optik_pe['proper_right']}/{optik_pe['both_right']}",
    )
    check(
        "TLEN = true fragment length on >= 99% of properly paired right pairs",
        optik_pe["tlen_ok"] >= MIN_PROPER_SHARE * optik_pe["proper_right"],
        f"{optik_pe['tlen_ok']}/{optik_pe['proper_right']}",
    )
    # The paired SAM must also survive the standard duplicate-marking chain, which reads the mate fields that a mapper most often gets wrong.
    if args.samtools:
        fixmate_bam = os.path.join(args.out, "fixmate.bam")
        st_ = os.path.join(args.out, "sorted.bam")
        markdup_bam = os.path.join(args.out, "markdup.bam")
        fixmate_proc = subprocess.run(
            [args.samtools, "fixmate", "-m", sams["optik paired"], fixmate_bam],
            capture_output=True,
            text=True,
            check=False,
        )
        sort_proc = subprocess.run(
            [args.samtools, "sort", "-o", st_, fixmate_bam],
            capture_output=True,
            text=True,
            check=False,
        )
        markdup_proc = subprocess.run(
            [args.samtools, "markdup", st_, markdup_bam],
            capture_output=True,
            text=True,
            check=False,
        )
        check(
            "samtools fixmate -m | sort | markdup run on optik's paired SAM",
            fixmate_proc.returncode == sort_proc.returncode == markdup_proc.returncode == 0,
            (fixmate_proc.stderr + sort_proc.stderr + markdup_proc.stderr)[-DETAIL_CHARS:],
        )
        if args.picard and not (os.path.isfile(args.picard) and os.access(args.picard, os.X_OK)):
            print(
                f"  [SKIP] Picard ValidateSamFile: {args.picard} is not an executable (PICARD_BIN stale?)"
            )
        elif args.picard:
            picard_proc = subprocess.run(
                [
                    args.picard,
                    "ValidateSamFile",
                    "-I",
                    st_,
                    "-R",
                    args.ref,
                    "-MODE",
                    "SUMMARY",
                    "-IGNORE_WARNINGS",
                    "true",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            txt = picard_proc.stdout + picard_proc.stderr
            check(
                "Picard ValidateSamFile (STRICT): no errors on optik's paired BAM",
                "No errors found" in txt,
                "\n".join(
                    text_line for text_line in txt.splitlines() if "ERROR" in text_line
                )[:DETAIL_CHARS],
            )
    print(f"\npaired_bench: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}")
    # The results file gets one row per mapper, with the same figures as the printout.
    mappers = []
    for name, stats in res.items():
        bands = "  ".join(
            f"{band_key}:{wrong}/{implied:.1f}"
            for band_key, (count, wrong, implied) in sorted(stats["bands"].items())
            if band_key >= MAPQ_BAND_WIDTH
        )
        row = {
            "mapper": name,
            "mates": stats["mates"],
            "correct": stats["correct"],
            "mapq30": stats["hq"],
            "mapq30_wrong": stats["hq_wrong"],
            "unmapped": stats["unmapped"],
        }
        if "paired" in name:
            row.update(
                proper_of_both_right=f"{stats['proper_right']}/{stats['both_right']}",
                tlen_right=stats["tlen_ok"],
                mate_field_errors=stats["mate_field_bad"],
            )
        mappers.append({**row, "wrong/implied by MAPQ band": bands})
    results_file.write_section(
        args.report,
        "paired",
        f"Paired-end mapping on {args.pairs} gen_pairs.py pairs at error {args.error}: optik paired, "
        "optik per mate and minimap2 -ax sr, scored against the simulation truth",
        [("mappers", mappers), ("checks", CHECKS)],
        optik=args.optik,
        notes=[f"SAMs: `{os.path.relpath(args.out, ROOT)}`"],
    )
    if args.gate:
        sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
