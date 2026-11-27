#!/usr/bin/env python3
"""real_reads.py - optik map on real sequencing reads, next to minimap2 (and BWA-MEM when given), where simulations cannot speak.
Real reads bring real error profiles, real repeat content and real library artefacts.

Real reads carry no placement truth, so the truth comes from the reads' own structure and from a published benchmark.

  anchored   MAPQ calibration by pair anchoring.
             The tools map the pairs as pairs.
             A mate is anchored when, in every tool's paired output, its partner has MAPQ >= 50, the pair is proper, and all tools place both mates at the same origin (within 20 bp, same strand).
             Its true locus is then fixed by a uniquely placed partner and the library's insert size.
             Each such mate is then mapped single-end by all tools, and its MAPQ is judged against the anchored locus, which covers the case that matters, a mate inside a repeat whose partner is unique.
             Mates the tools place differently when paired are left out and counted, so the truth set leans easy.
             The per-bin test is the Poisson test of mapq_calibration.py.

  variants   Germline small-variant calling against a GIAB benchmark (such as HG002 v4.2.1).
             Each tool's paired SAM goes through samtools sort, bcftools mpileup | call -mv and a QUAL >= 20, DP >= 5 filter, is normalized (bcftools norm -m -both -f REF) and is compared with the normalized truth inside the benchmark BED by RTG vcfeval (--rtg). vcfeval is haplotype- and genotype-aware and is the engine of GIAB's hap.py, so an indel written differently but giving the same haplotype counts as found.
             Without RTG, calls are matched allele by allele, which is stricter, since one haplotype written two ways counts as an FP and an FN.
             Precision and recall are reported for SNPs and indels.

  longreads  Real HiFi or ONT long reads, optik against minimap2.
             It reports primary placement agreement, disagreements at MAPQ >= 30 on either side, and indels of 50 bp or more in primary CIGARs.
             A real sample has true SVs, so indels one tool reports at a locus the other aligns without them are listed per tool.

Usage
  real_reads.py anchored  --ref REF --r1 R1.fq --r2 R2.fq --out DIR [--pairs N] [--bwa-index PREFIX]
  real_reads.py variants  --ref REF --r1 R1.fq --r2 R2.fq --truth-vcf V --truth-bed B --region chr20 --out DIR [--bwa-index PREFIX]
  real_reads.py longreads --ref REF --reads R.fq --preset map-hifi|map-ont --out DIR
Common options --optik BIN --minimap2 BIN --bwa BIN --samtools BIN --bcftools BIN --threads T
"""

import argparse
import collections
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from typing import NamedTuple

import hw
import memcap
import report as results_file
from diagnostics import LOG_TAIL_CHARS
from fastq import FASTQ_LINES_PER_RECORD
from mapq_calibration import (
    BINS,
    OVERCONFIDENT_P,
    overconfident,
    poisson_tail,
)
from sam_spec import (
    HIGH_MAPQ,
    PHRED_SCALE,
    SamCol,
    SamFlag,
    VcfCol,
)

# variants mode results per tool and type: (TP, FP, FN, precision, recall, F1)
VARIANT_F1 = 5
# SAM lines are split only as far as needed: primaries() reads QNAME..CIGAR and longreads() QNAME..RNEXT.
# Everything after that stays one field.
SAM_SPLIT_PRIMARIES = SamCol.RNEXT
SAM_SPLIT_LONGREADS = SamCol.PNEXT
# mates are numbered 1 and 2, so a mate's partner is MATE_NUMBER_SUM minus its number
MATE_NUMBER_SUM = 3
# anchored mode: a partner anchors its mate when it is a proper pair at this MAPQ or more in every tool, and placements agree within ANCHOR_TOL bp
ANCHOR_MIN_MAPQ = 50
ANCHOR_TOL = 20
# variants mode: the denominator floor of F1, and how far optik's F1 may trail another mapper's
F1_EPS = 1e-9
F1_TOLERANCE = 0.005
# In longreads mode, indels at least BIG_INDEL_LEN long are compared.
# Two primaries agree when they lie within LOCUS_TOL bp of each other.
BIG_INDEL_LEN = 50
LOCUS_TOL = 50
# two large indels are the same when within INDEL_POS_TOL bp and their lengths within max(INDEL_MIN_LEN_TOL, INDEL_LEN_TOL_FRAC x length)
INDEL_POS_TOL = 100
INDEL_MIN_LEN_TOL = 10
INDEL_LEN_TOL_FRAC = 0.2

MEM_CAP = (
    None  # set from --mem-cap, so every tool runs under its own resident-memory cap
)


# Mates are numbered as in their /1 and /2 name suffixes.
MATE1, MATE2 = 1, 2
# A MAPQ bin counts its placements, its wrong ones and the errors its MAPQs imply, in this order.
BIN_TOTAL, BIN_WRONG, BIN_IMPLIED = 0, 1, 2
# F1 is the harmonic mean of precision and recall, 2PR / (P + R).
F1_HARMONIC_FACTOR = 2


def cap_bytes():
    """Each tool's resident-memory cap, applied by memcap on any Linux.

    A tool that outgrows the machine is then killed alone rather than taking the session with it, as on a 32 GB machine running an hg38 mapper next to another job."""
    return memcap.parse_size(MEM_CAP)


def run_shell(cmd, out=None, env=None):
    """Run cmd under the memory cap, writing stdout to `out` when given, and stop on failure."""
    with open(out, "w") if out else open(os.devnull, "w") as handle:
        proc = memcap.run(
            [str(arg) for arg in cmd],
            cap_bytes(),
            stdout=handle if out else subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    if proc.returncode != 0:
        sys.exit(
            f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}\n{proc.stderr[-LOG_TAIL_CHARS:]}"
        )
    return proc


def head_fastq(src, dst, num_reads):
    """Link or copy the first n records of a FASTQ.
    n = 0 means all of them, given as a symlink."""
    if os.path.exists(dst):
        return dst
    if num_reads <= 0:
        os.symlink(os.path.abspath(src), dst)
        return dst
    with open(src) as in_handle, open(dst, "w") as out_handle:
        for i, line in enumerate(in_handle):
            if i >= FASTQ_LINES_PER_RECORD * num_reads:
                break
            out_handle.write(line)
    return dst


def origin(fields):
    """(ref, strand, origin) of a record, where origin is the leftmost aligned base minus a leading clip."""
    ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
    return fields[SamCol.RNAME], "-" if int(fields[SamCol.FLAG]) & SamFlag.REVERSE else "+", int(fields[SamCol.POS]) - 1 - lead


class Placement(NamedTuple):
    """Records where a read lies, by reference, strand and origin."""

    ref: str
    strand: str
    origin: int


class Primary(NamedTuple):
    """A read's primary alignment, as primaries() reports it."""

    ref: str
    strand: str
    origin: int
    mapq: int
    flag: int

    @property
    def placement(self):
        return Placement(self.ref, self.strand, self.origin)


def primaries(sam, paired):
    """Map name (or name/mate) to (ref, strand, origin, mapq, flag) of each primary."""
    out = {}
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            fields = line.split("\t", SAM_SPLIT_PRIMARIES)
            flag = int(fields[SamCol.FLAG])
            if flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                continue
            name = fields[SamCol.QNAME]
            if paired:
                name = f"{name}/{MATE1 if flag & SamFlag.READ1 else MATE2}"
            elif name.endswith(("/1", "/2")):
                pass
            if flag & SamFlag.UNMAPPED:
                out[name] = None
            else:
                ref_name, strand, origin_pos = origin(fields)
                out[name] = Primary(ref_name, strand, origin_pos, int(fields[SamCol.MAPQ]), flag)
    return out


def map_cmds(args, ref, reads1, reads2=None, se_tag=None):
    """{tool: argv} for paired mapping when r2 is given, else for single-end mapping."""
    reads = [reads1] + ([reads2] if reads2 else [])
    cmds = {
        "optik": [args.optik, "map", "-r", ref, "-t", args.threads]
        + [arg for read_path in reads for arg in ("-q", read_path)],
        "minimap2": [args.minimap2, "-ax", "sr", "-t", args.threads, ref] + reads,
    }
    if args.bwa_index:
        cmds["bwa"] = [args.bwa, "mem", "-t", args.threads, args.bwa_index] + reads
    return cmds


def run_maps(args, tag, ref, reads1, reads2=None):
    """Map with every tool, reusing a non-empty SAM from an earlier run, and return {tool: SAM path}."""
    sams = {}
    for tool, cmd in map_cmds(args, ref, reads1, reads2).items():
        sam = os.path.join(args.out, f"{tool}_{tag}.sam")
        if not os.path.exists(sam) or os.path.getsize(sam) == 0:
            print(f"  mapping {tag} with {tool}", flush=True)
            run_shell(cmd, out=sam + ".part")
            os.replace(sam + ".part", sam)
        sams[tool] = sam
    return sams


def calib_line(label, rows):
    """Print a MAPQ calibration line for (mapq, correct) rows and return (over-confident bins, MAPQ>=30, wrong)."""
    bins = {bin_range: [0, 0, 0.0] for bin_range in BINS}
    for mapq, correct in rows:
        for bin_range in BINS:
            if bin_range[0] <= mapq <= bin_range[1]:
                bin_counts = bins[bin_range]
                bin_counts[BIN_TOTAL] += 1
                bin_counts[BIN_WRONG] += not correct
                bin_counts[BIN_IMPLIED] += PHRED_SCALE ** (-mapq / PHRED_SCALE)
                break
    hq_total = sum(bin_counts[BIN_TOTAL] for bin_range, bin_counts in bins.items() if bin_range[0] >= HIGH_MAPQ)
    hqw = sum(bin_counts[BIN_WRONG] for bin_range, bin_counts in bins.items() if bin_range[0] >= HIGH_MAPQ)
    cells, over = [], []
    for (bin_lo, bin_hi), (count, wrong, lam) in bins.items():
        if not count:
            continue
        cells.append(f"{bin_lo}-{bin_hi}:{wrong}/{lam:.1f}" if bin_lo else f"0:{wrong}/{count}")
        # A bin is flagged only when its excess errors are also unlikely by chance under a Poisson count.
        if bin_lo > 0 and overconfident(wrong, lam) and poisson_tail(wrong, lam) < OVERCONFIDENT_P:
            over.append(f"MAPQ {bin_lo}-{bin_hi}: {wrong} wrong vs {lam:.1f} implied")
    print(
        f"  {label:10} mapped {len(rows):8}  MAPQ>=30 {hq_total:8} wrong {hqw:4} | wrong/implied: "
        + "  ".join(cells)
    )
    ROWS.append(
        {
            "tool": label,
            "mapped": len(rows),
            "mapq30": hq_total,
            "mapq30_wrong": hqw,
            "wrong/implied by MAPQ bin": "  ".join(cells),
        }
    )
    return over, hq_total, hqw


def anchored(args):
    """Anchored mode, which judges single-end MAPQ against loci fixed by uniquely placed partners."""
    reads1 = head_fastq(args.r1, os.path.join(args.out, "a_R1.fq"), args.pairs)
    reads2 = head_fastq(args.r2, os.path.join(args.out, "a_R2.fq"), args.pairs)
    paired_recs = {tool: primaries(sam, True) for tool, sam in run_maps(args, "pe", args.ref, reads1, reads2).items()}
    tools = list(paired_recs)
    truth, left_out = {}, collections.Counter()
    names = set().union(*[set(recs_by_name) for recs_by_name in paired_recs.values()])
    for read_name in names:
        base, mate_str = read_name.rsplit("/", 1)
        partner = f"{base}/{MATE_NUMBER_SUM - int(mate_str)}"
        recs = [paired_recs[tool].get(read_name) for tool in tools]
        precs = [paired_recs[tool].get(partner) for tool in tools]
        if any(rec is None for rec in recs + precs):
            left_out["a mate unmapped by some tool"] += 1
            continue
        if not all(partner_rec.mapq >= ANCHOR_MIN_MAPQ and partner_rec.flag & SamFlag.PROPER_PAIR for partner_rec in precs):
            left_out["partner not MAPQ>=50 proper in every tool"] += 1
            continue

        def same(placements):
            """True when all placements share reference and strand and lie within 20 bp."""
            return all(
                place.ref == placements[0].ref
                and place.strand == placements[0].strand
                and abs(place.origin - placements[0].origin) <= ANCHOR_TOL
                for place in placements
            )

        if not same(precs) or not same(recs):
            left_out["tools place the pair differently"] += 1
            continue
        truth[read_name] = recs[0].placement
    print(f"anchored mates: {len(truth)} of {len(names)} ({dict(left_out)})")
    with open(
        os.path.join(args.out, "anchored_truth.tsv"), "w"
    ) as handle:  # kept for offline MAPQ studies
        for read_name, (ref_name, strand, origin_pos) in truth.items():
            handle.write(f"{read_name}\t{ref_name}\t{strand}\t{origin_pos}\n")
    rep = {}
    failed = []
    for mate, fq_path in ((MATE1, reads1), (MATE2, reads2)):
        for tool, sam in run_maps(args, f"se{mate}", args.ref, fq_path).items():
            for read_name, rec in primaries(sam, False).items():
                # Mappers differ in whether they keep the /1 /2 suffix, so names are normalised to include it.
                key = read_name if read_name.endswith(("/1", "/2")) else f"{read_name}/{mate}"
                if key in truth:
                    rep.setdefault(tool, []).append(
                        rec
                        and (
                            rec.mapq,
                            rec.ref == truth[key].ref
                            and rec.strand == truth[key].strand
                            and abs(rec.origin - truth[key].origin) <= ANCHOR_TOL,
                        )
                    )
    print("single-end MAPQ of anchored mates (wrong / implied by MAPQ, per bin):")
    summary = {}
    for tool in tools:
        # Unmapped mates are recorded as None and counted separately from the MAPQ bins.
        rows = [row for row in rep.get(tool, []) if row]
        unm = sum(1 for row in rep.get(tool, []) if not row)
        over, hq_total, hqw = calib_line(tool, rows)
        print(f"  {'':10} unmapped {unm}")
        ROWS[-1]["unmapped"] = unm
        summary[tool] = (hq_total, hqw, over)
        if tool == "optik":
            failed += over
    return failed


def norm_vcf(args, vcf, out, region, bed=None):
    """Normalize a VCF to `region` and return its {snp, indel} allele keys, inside `bed` when given."""
    # -r needs an index, which is always rebuilt, because the filtered VCF is rewritten on each run and a stale index from an earlier run makes norm return nothing.
    run_shell([args.bcftools, "index", "-f", vcf])
    cmd = [
        args.bcftools,
        "norm",
        "-m",
        "-both",
        "-f",
        args.ref,
        "-r",
        region,
        "-Oz",
        "-o",
        out,
        vcf,
    ]
    run_shell(cmd)
    run_shell([args.bcftools, "index", "-f", out])
    keys = {"snp": set(), "indel": set()}
    view = [args.bcftools, "view", "-H"] + (["-T", bed] if bed else []) + [out]
    proc = subprocess.run(
        [str(arg) for arg in view], capture_output=True, text=True, check=True
    )
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        ref, alt = fields[VcfCol.REF], fields[VcfCol.ALT]
        if alt in ("*", ".") or "<" in alt:
            continue
        kind = "snp" if len(ref) == len(alt) == 1 else "indel"
        keys[kind].add((fields[VcfCol.CHROM], int(fields[VcfCol.POS]), ref, alt))
    return keys


def vcfeval(args, truth_vcf, calls_vcf, bed, out):
    """Compare calls with the truth by RTG vcfeval, haplotype- and genotype-aware.

    vcfeval is the engine of GIAB's hap.py, so an indel or SNP written differently but giving the same haplotypes counts as a match.
    Returns (TP, FP, FN) per type, with recall taken from the baseline side and precision from the call side."""
    # vcfeval reads the reference as an RTG SDF, which is built once per region and then reused.
    sdf = os.path.join(args.out, f"{args.region}.sdf")
    if not os.path.isdir(sdf):
        ref_fa = os.path.join(args.out, f"{args.region}.fa")
        run_shell([args.samtools, "faidx", args.ref, args.region], out=ref_fa)
        run_shell([args.rtg, "format", "-o", sdf, ref_fa])
    for vcf in (truth_vcf, calls_vcf):
        run_shell([args.bcftools, "index", "-t", "-f", vcf])
    # vcfeval refuses to write into an existing directory, so the previous run's output is removed first.
    shutil.rmtree(out, ignore_errors=True)
    env = dict(os.environ, RTG_MEM=os.environ.get("RTG_MEM", "6g"))
    run_shell(
        [
            args.rtg,
            "vcfeval",
            "-b",
            truth_vcf,
            "-c",
            calls_vcf,
            "-e",
            bed,
            "-t",
            sdf,
            "-o",
            out,
            "--threads",
            args.threads,
        ],
        env=env,
    )

    def count(name):
        """SNP and indel counts of one vcfeval output VCF."""
        counts = {"snp": 0, "indel": 0}
        proc = subprocess.run(
            [str(args.bcftools), "view", "-H", os.path.join(out, name)],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in proc.stdout.splitlines():
            fields = line.split("\t")
            ref, alts = fields[VcfCol.REF], [alt for alt in fields[VcfCol.ALT].split(",") if alt not in ("*", ".")]
            counts[
                "snp"
                if alts and all(len(ref) == len(alt) == 1 for alt in alts)
                else "indel"
            ] += 1
        return counts

    # True positives come in the baseline representation for recall and in the call representation for precision.
    tpb, tpc, false_pos, false_neg = (
        count(vcf_name) for vcf_name in ("tp-baseline.vcf.gz", "tp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz")
    )
    return {kind: (tpb[kind], tpc[kind], false_pos[kind], false_neg[kind]) for kind in ("snp", "indel")}


def variants(args):
    """Variants mode, which calls germline variants from each tool's mapping and scores them against GIAB."""
    bed = os.path.join(args.out, "truth_region.bed")
    # The truth regions are cut to the chromosome being called, so every tool is scored on the same intervals.
    with open(args.truth_bed) as in_handle, open(bed, "w") as out_handle:
        for line in in_handle:
            if line.split("\t")[0] == args.region:
                out_handle.write(line)
    truth = norm_vcf(
        args, args.truth_vcf, os.path.join(args.out, "truth.norm.vcf.gz"), args.region, bed
    )
    print(
        f"truth in {args.region} benchmark regions: {len(truth['snp'])} SNPs, {len(truth['indel'])} indels; "
        f"scored by {'rtg vcfeval (haplotype- and genotype-aware)' if args.rtg else 'exact allele match'}"
    )
    failed, res = [], {}
    for tool, sam in run_maps(args, "var", args.ref, args.r1, args.r2).items():
        bam = os.path.join(args.out, f"{tool}.bam")
        # Sorted BAMs and calls left by an earlier run are reused, because a rerun only needs to score them again.
        if not os.path.exists(bam + ".bai"):
            run_shell([args.samtools, "sort", "-@", args.threads, "-m", "1G", "-o", bam, sam])
            run_shell([args.samtools, "index", bam])
        raw = os.path.join(args.out, f"{tool}.calls.vcf.gz")
        if not os.path.exists(raw):
            mpileup_proc = memcap.popen(
                [
                    str(arg)
                    for arg in (
                        args.bcftools,
                        "mpileup",
                        "-f",
                        args.ref,
                        "-r",
                        args.region,
                        "-a",
                        "AD,DP",
                        "-Ou",
                        bam,
                    )
                ],
                cap_bytes(),
                stdout=subprocess.PIPE,
            )
            call_proc = subprocess.run(
                [str(args.bcftools), "call", "-mv", "-Oz", "-o", raw],
                stdin=mpileup_proc.stdout,
                check=False,
            )
            mpileup_proc.wait()
            if call_proc.returncode or mpileup_proc.returncode:
                sys.exit(f"bcftools failed for {tool}")
        filt = os.path.join(args.out, f"{tool}.filt.vcf.gz")
        # Calls below QUAL 20 or depth 5 are dropped before scoring.
        run_shell([args.bcftools, "view", "-i", "QUAL>=20 && INFO/DP>=5", "-Oz", "-o", filt, raw])
        calls = norm_vcf(
            args, filt, os.path.join(args.out, f"{tool}.norm.vcf.gz"), args.region, bed
        )
        # A tool with no SNP calls at all points to a broken pipeline, so the run stops instead of scoring it as zero recall.
        if not calls["snp"]:
            sys.exit(
                f"{tool}: no calls in {args.region} after filtering -- pipeline error, not a result"
            )
        res[tool] = {}
        eval_res = (
            vcfeval(
                args,
                os.path.join(args.out, "truth.norm.vcf.gz"),
                os.path.join(args.out, f"{tool}.norm.vcf.gz"),
                bed,
                os.path.join(args.out, f"vcfeval_{tool}"),
            )
            if args.rtg
            else None
        )
        for kind in ("snp", "indel"):
            if eval_res:
                true_pos, tpc, false_pos, false_neg = eval_res[kind]
                prec = tpc / max(1, tpc + false_pos)
            else:
                # Without RTG, exact allele match after bcftools norm, which is stricter, because equivalent representations of one haplotype count as an FP and an FN.
                true_pos = len(calls[kind] & truth[kind])
                false_pos = len(calls[kind] - truth[kind])
                prec = true_pos / max(1, true_pos + false_pos)
                false_neg = len(truth[kind] - calls[kind])
            # max(1, ...) guards against an empty call or truth set.
            rec = true_pos / max(1, true_pos + false_neg)
            f1_score = F1_HARMONIC_FACTOR * prec * rec / max(F1_EPS, prec + rec)
            res[tool][kind] = (true_pos, false_pos, false_neg, prec, rec, f1_score)
            ROWS.append(
                {
                    "tool": tool,
                    "type": kind,
                    "TP": true_pos,
                    "FP": false_pos,
                    "FN": false_neg,
                    "precision": round(prec, 4),
                    "recall": round(rec, 4),
                    "F1": round(f1_score, 4),
                }
            )
            print(
                f"  {tool:9} {kind:5} TP {true_pos:6} FP {false_pos:5} FN {false_neg:5}  precision {prec:.4f}  recall {rec:.4f}  F1 {f1_score:.4f}"
            )
    if "optik" in res:
        for other in res:
            if other == "optik":
                continue
            for kind in ("snp", "indel"):
                # optik fails only when it trails another mapper by more than the tolerance. Small F1 noise is accepted.
                if res["optik"][kind][VARIANT_F1] < res[other][kind][VARIANT_F1] - F1_TOLERANCE:
                    failed.append(
                        f"{kind} F1: optik {res['optik'][kind][VARIANT_F1]:.4f} < {other} {res[other][kind][VARIANT_F1]:.4f} - {F1_TOLERANCE}"
                    )
    return failed


def longreads(args):
    """Longreads mode, which compares optik's and minimap2's primaries and large indels on real long reads."""
    sams = {}
    for tool, cmd in (
        ("optik", [args.optik, "map", "-r", args.ref, "-q", args.reads, "-t", args.threads]),
        ("minimap2", [args.minimap2, "-ax", args.preset, "-t", args.threads, args.ref, args.reads]),
    ):
        sam = os.path.join(args.out, f"{tool}_long.sam")
        if not os.path.exists(sam) or os.path.getsize(sam) == 0:
            print(f"  mapping long reads with {tool}", flush=True)
            # Writing to .part and renaming means an interrupted run never leaves a SAM that looks complete.
            run_shell(cmd, out=sam + ".part")
            os.replace(sam + ".part", sam)
        sams[tool] = sam
    recs, big = {}, {}
    for tool, sam in sams.items():
        by_name, indels = {}, []
        with open(sam) as handle:
            for line in handle:
                if line[0] == "@":
                    continue
                fields = line.split("\t", SAM_SPLIT_LONGREADS)
                flag = int(fields[SamCol.FLAG])
                if flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                    continue
                if flag & SamFlag.UNMAPPED:
                    by_name[fields[SamCol.QNAME]] = None
                    continue
                ref_name, strand, origin_pos = origin(fields)
                by_name[fields[SamCol.QNAME]] = Primary(ref_name, strand, origin_pos, int(fields[SamCol.MAPQ]), flag)
                pos = int(fields[SamCol.POS]) - 1
                for op_len, cigar_op in re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR]):
                    op_len = int(op_len)
                    if cigar_op in "ID" and op_len >= BIG_INDEL_LEN:
                        indels.append((ref_name, pos, cigar_op, op_len, fields[SamCol.QNAME]))
                    # Only reference-consuming ops move pos, so each indel gets its reference coordinate.
                    if cigar_op in "MDN=X":
                        pos += op_len
        recs[tool], big[tool] = by_name, indels
    # A disagreement where either side reports MAPQ >= 30 means one mapper is confidently wrong, which is what dis_hq_o and dis_hq_m count.
    names = set(recs["optik"]) | set(recs["minimap2"])
    agree = dis = dis_hq_o = dis_hq_m = un_o = un_m = 0
    for read_name in names:
        optik_rec, mm_rec = recs["optik"].get(read_name), recs["minimap2"].get(read_name)
        un_o += optik_rec is None
        un_m += mm_rec is None
        if optik_rec is None or mm_rec is None:
            continue
        if (
            optik_rec.ref == mm_rec.ref
            and optik_rec.strand == mm_rec.strand
            and abs(optik_rec.origin - mm_rec.origin) <= max(LOCUS_TOL, 0)
        ):
            agree += 1
        else:
            dis += 1
            dis_hq_o += optik_rec.mapq >= HIGH_MAPQ
            dis_hq_m += mm_rec.mapq >= HIGH_MAPQ
    print(
        f"  reads {len(names)}: unmapped optik {un_o}, minimap2 {un_m}; same primary locus {agree}, "
        f"different {dis} (MAPQ>=30 there: optik {dis_hq_o}, minimap2 {dis_hq_m})"
    )

    # Large indels are either shared (same read, position within 100 bp, same operation, length within 20%) or tool-only.

    by_read = {tool: collections.defaultdict(list) for tool in big}
    for tool, indels in big.items():
        for ref_name, indel_pos, cigar_op, op_len, read_name in indels:
            by_read[tool][read_name].append((ref_name, indel_pos, cigar_op, op_len))
    only = {}
    for tool, other in (("optik", "minimap2"), ("minimap2", "optik")):
        only_count = 0
        for read_name, evs in by_read[tool].items():
            for ref_name, indel_pos, cigar_op, op_len in evs:
                if not any(
                    ref_name == other_ref
                    and cigar_op == op2
                    and abs(indel_pos - other_pos) <= INDEL_POS_TOL
                    and abs(op_len - other_len) <= max(INDEL_MIN_LEN_TOL, INDEL_LEN_TOL_FRAC * op_len)
                    for other_ref, other_pos, op2, other_len in by_read[other].get(read_name, [])
                ):
                    only_count += 1
        only[tool] = only_count
    print(
        f"  >= 50 bp indels in primaries: optik {len(big['optik'])} ({only['optik']} not in minimap2's record "
        f"of the same read), minimap2 {len(big['minimap2'])} ({only['minimap2']} not in optik's)"
    )
    for tool, other in (("optik", "minimap2"), ("minimap2", "optik")):
        ROWS.append(
            {
                "tool": tool,
                "reads": len(names),
                "unmapped": un_o if tool == "optik" else un_m,
                "same_primary_locus": agree,
                "different_locus": dis,
                "mapq30_at_different": dis_hq_o if tool == "optik" else dis_hq_m,
                "indels_50bp": len(big[tool]),
                f"indels_50bp_not_in_{other}": only[tool],
            }
        )
    return []


# Rows of the results file, appended by the mode that runs
ROWS = []
# Heading of each mode's section in the results file
MODE_TITLES = {
    "anchored": "MAPQ of real read mates whose locus is fixed by a uniquely placed partner",
    "variants": "Germline small variants called from each mapper's SAM against the GIAB benchmark",
    "longreads": "Real long reads, optik against minimap2: primary agreement and >= 50 bp indels",
}


def main():
    """Run the anchored, variants or longreads mode."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Each mode reads its own inputs, so --r1/--r2 serve anchored and variants while --reads serves longreads.
    parser.add_argument("mode", choices=("anchored", "variants", "longreads"))
    parser.add_argument("--ref", required=True)
    parser.add_argument("--r1")
    parser.add_argument("--r2")
    parser.add_argument("--reads")
    parser.add_argument("--preset", default="map-hifi")
    parser.add_argument(
        "--pairs", type=int, default=0, help="anchored: use the first N pairs (0 = all)"
    )
    parser.add_argument("--truth-vcf")
    parser.add_argument("--truth-bed")
    parser.add_argument("--region", default="chr20")
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--threads",
        default=str(hw.cpus()),
        help="threads of every tool (default: the usable CPUs)",
    )
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    parser.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    # Tool paths fall back from PATH to the copies bench/build_tools.sh puts under temp/, then to the Nix GC roots.
    parser.add_argument(
        "--bwa",
        default=shutil.which("bwa")
        or next(
            (
                path
                for path in (
                    os.path.join(ROOT, "temp/bwa/bwa"),
                    os.path.join(ROOT, "temp/gcroots/bwa/bin/bwa"),
                )
                if os.path.exists(path)
            ),
            os.path.join(ROOT, "temp/bwa/bwa"),
        ),
    )
    parser.add_argument("--bwa-index", help="BWA index prefix (adds BWA-MEM as a baseline)")
    parser.add_argument("--samtools", default=shutil.which("samtools"))
    parser.add_argument(
        "--bcftools",
        default=os.environ.get("BCFTOOLS_BIN")
        or shutil.which("bcftools")
        or (
            os.path.join(ROOT, "temp/bcftools/bcftools")
            if os.path.exists(os.path.join(ROOT, "temp/bcftools/bcftools"))
            else os.path.join(ROOT, "temp/gcroots/bcftools/bin/bcftools")
        ),
    )
    parser.add_argument(
        "--rtg",
        default=shutil.which("rtg")
        or next(
            (
                path
                for path in (
                    os.path.join(ROOT, "temp/rtg/rtg"),
                    os.path.join(ROOT, "temp/gcroots/rtg/bin/rtg"),
                )
                if os.path.exists(path)
            ),
            None,
        ),
        help="RTG Tools rtg (vcfeval: haplotype-aware comparison, as GIAB hap.py); "
        "without it calls are matched allele by allele",
    )
    parser.add_argument("--gate", action="store_true")
    parser.add_argument(
        "--mem-cap",
        default=hw.default_cap_str(),
        help="resident-memory cap of each tool ('' = uncapped; default: usable memory less 2 GiB)",
    )
    results_file.add_argument(parser)
    args = parser.parse_args()
    global MEM_CAP
    # An empty --mem-cap means no cap, which the runners take as None.
    MEM_CAP = args.mem_cap or None
    os.makedirs(args.out, exist_ok=True)
    # Each mode returns its list of failures, and an empty list is a pass.
    failed = {"anchored": anchored, "variants": variants, "longreads": longreads}[
        args.mode
    ](args)
    for failure in failed:
        print("  [FAIL]", failure)
    print(f"real_reads {args.mode}: " + ("FAILED" if failed else "PASS"))
    results_file.write_section(
        args.report,
        f"real_{args.mode}",
        MODE_TITLES[args.mode] + f" (reference {os.path.basename(args.ref)})",
        [("", ROWS)],
        optik=args.optik,
        notes=[f"- FAIL {failure}" for failure in failed]
        + [f"SAMs and logs: `{os.path.relpath(args.out, ROOT)}`"],
    )
    return 1 if failed and args.gate else 0


if __name__ == "__main__":
    sys.exit(main())
