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
             Each tool's paired SAM goes through samtools sort, bcftools mpileup | call -mv and a QUAL >= 20, DP >= 5 filter, is normalized (bcftools norm -m -both -f REF) and is compared with the normalized truth inside the benchmark BED by RTG vcfeval (--rtg).
             vcfeval is haplotype- and genotype-aware and is the engine of GIAB's hap.py, so an indel written differently but giving the same haplotype counts as found.
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

import hw
import memcap
import report as results_file
from mapq_calibration import BINS, poisson_tail
from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_READ1,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
)

MEM_CAP = (
    None  # set from --mem-cap, so every tool runs under its own resident-memory cap
)


def cap_bytes():
    """Each tool's resident-memory cap, applied by memcap on any Linux.

    A tool that outgrows the machine is then killed alone rather than taking the session with it, as on a 32 GB machine running an hg38 mapper next to another job."""
    return memcap.parse_size(MEM_CAP)


def sh(cmd, out=None, env=None):
    """Run cmd under the memory cap, writing stdout to `out` when given, and stop on failure."""
    with open(out, "w") if out else open(os.devnull, "w") as fh:
        p = memcap.run(
            [str(c) for c in cmd],
            cap_bytes(),
            stdout=fh if out else subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    if p.returncode != 0:
        sys.exit(
            f"command failed ({p.returncode}): {' '.join(map(str, cmd))}\n{p.stderr[-1500:]}"
        )
    return p


def head_fastq(src, dst, n):
    """Link or copy the first n records of a FASTQ.
    n = 0 means all of them, given as a symlink."""
    if os.path.exists(dst):
        return dst
    if n <= 0:
        os.symlink(os.path.abspath(src), dst)
        return dst
    with open(src) as f, open(dst, "w") as o:
        for i, line in enumerate(f):
            if i >= 4 * n:
                break
            o.write(line)
    return dst


def origin(f):
    """(ref, strand, origin) of a record, where origin is the leftmost aligned base minus a leading clip."""
    ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
    lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
    return f[2], "-" if int(f[1]) & SAM_FLAG_REVERSE else "+", int(f[3]) - 1 - lead


def primaries(sam, paired):
    """Map name (or name/mate) to (ref, strand, origin, mapq, flag) of each primary."""
    out = {}
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)
            flag = int(f[1])
            if flag & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                continue
            name = f[0]
            if paired:
                name = f"{name}/{1 if flag & SAM_FLAG_READ1 else 2}"
            elif name.endswith(("/1", "/2")):
                pass
            if flag & SAM_FLAG_UNMAPPED:
                out[name] = None
            else:
                r, s, o = origin(f)
                out[name] = (r, s, o, int(f[4]), flag)
    return out


def map_cmds(a, ref, r1, r2=None, se_tag=None):
    """{tool: argv} for paired mapping when r2 is given, else for single-end mapping."""
    reads = [r1] + ([r2] if r2 else [])
    cmds = {
        "optik": [a.optik, "map", "-r", ref, "-t", a.threads]
        + [a for r in reads for a in ("-q", r)],
        "minimap2": [a.minimap2, "-ax", "sr", "-t", a.threads, ref] + reads,
    }
    if a.bwa_index:
        cmds["bwa"] = [a.bwa, "mem", "-t", a.threads, a.bwa_index] + reads
    return cmds


def run_maps(a, tag, ref, r1, r2=None):
    """Map with every tool, reusing a non-empty SAM from an earlier run, and return {tool: SAM path}."""
    sams = {}
    for tool, cmd in map_cmds(a, ref, r1, r2).items():
        sam = os.path.join(a.out, f"{tool}_{tag}.sam")
        if not os.path.exists(sam) or os.path.getsize(sam) == 0:
            print(f"  mapping {tag} with {tool}", flush=True)
            sh(cmd, out=sam + ".part")
            os.replace(sam + ".part", sam)
        sams[tool] = sam
    return sams


def calib_line(label, rows):
    """Print a MAPQ calibration line for (mapq, correct) rows and return (over-confident bins, MAPQ>=30, wrong)."""
    bins = {b: [0, 0, 0.0] for b in BINS}
    for q, ok in rows:
        for b in BINS:
            if b[0] <= q <= b[1]:
                x = bins[b]
                x[0] += 1
                x[1] += not ok
                x[2] += 10 ** (-q / 10)
                break
    hq = sum(x[0] for b, x in bins.items() if b[0] >= HIGH_MAPQ)
    hqw = sum(x[1] for b, x in bins.items() if b[0] >= HIGH_MAPQ)
    cells, over = [], []
    for (lo, hi), (n, w, lam) in bins.items():
        if not n:
            continue
        cells.append(f"{lo}-{hi}:{w}/{lam:.1f}" if lo else f"0:{w}/{n}")
        if lo > 0 and w > 2 * lam + 3 and poisson_tail(w, lam) < 1e-3:
            over.append(f"MAPQ {lo}-{hi}: {w} wrong vs {lam:.1f} implied")
    print(
        f"  {label:10} mapped {len(rows):8}  MAPQ>=30 {hq:8} wrong {hqw:4} | wrong/implied: "
        + "  ".join(cells)
    )
    ROWS.append(
        {
            "tool": label,
            "mapped": len(rows),
            "mapq30": hq,
            "mapq30_wrong": hqw,
            "wrong/implied by MAPQ bin": "  ".join(cells),
        }
    )
    return over, hq, hqw


def anchored(a):
    """Anchored mode, which judges single-end MAPQ against loci fixed by uniquely placed partners."""
    r1 = head_fastq(a.r1, os.path.join(a.out, "a_R1.fq"), a.pairs)
    r2 = head_fastq(a.r2, os.path.join(a.out, "a_R2.fq"), a.pairs)
    pe = {t: primaries(s, True) for t, s in run_maps(a, "pe", a.ref, r1, r2).items()}
    tools = list(pe)
    truth, left_out = {}, collections.Counter()
    names = set().union(*[set(d) for d in pe.values()])
    for nm in names:
        base, m = nm.rsplit("/", 1)
        partner = f"{base}/{3 - int(m)}"
        recs = [pe[t].get(nm) for t in tools]
        precs = [pe[t].get(partner) for t in tools]
        if any(x is None for x in recs + precs):
            left_out["a mate unmapped by some tool"] += 1
            continue
        if not all(p[3] >= 50 and p[4] & SAM_FLAG_PROPER_PAIR for p in precs):
            left_out["partner not MAPQ>=50 proper in every tool"] += 1
            continue

        def same(xs):
            """True when all placements share reference and strand and lie within 20 bp."""
            return all(
                x[0] == xs[0][0] and x[1] == xs[0][1] and abs(x[2] - xs[0][2]) <= 20
                for x in xs
            )

        if not same(precs) or not same(recs):
            left_out["tools place the pair differently"] += 1
            continue
        truth[nm] = recs[0][:3]
    print(f"anchored mates: {len(truth)} of {len(names)} ({dict(left_out)})")
    with open(
        os.path.join(a.out, "anchored_truth.tsv"), "w"
    ) as fo:  # kept for offline MAPQ studies
        for nm, (r, s, o) in truth.items():
            fo.write(f"{nm}\t{r}\t{s}\t{o}\n")
    rep = {}
    failed = []
    for mate, fq in ((1, r1), (2, r2)):
        for t, s in run_maps(a, f"se{mate}", a.ref, fq).items():
            for nm, rec in primaries(s, False).items():
                key = nm if nm.endswith(("/1", "/2")) else f"{nm}/{mate}"
                if key in truth:
                    rep.setdefault(t, []).append(
                        rec
                        and (
                            rec[3],
                            rec[0] == truth[key][0]
                            and rec[1] == truth[key][1]
                            and abs(rec[2] - truth[key][2]) <= 20,
                        )
                    )
    print("single-end MAPQ of anchored mates (wrong / implied by MAPQ, per bin):")
    summary = {}
    for t in tools:
        rows = [x for x in rep.get(t, []) if x]
        unm = sum(1 for x in rep.get(t, []) if not x)
        over, hq, hqw = calib_line(t, rows)
        print(f"  {'':10} unmapped {unm}")
        ROWS[-1]["unmapped"] = unm
        summary[t] = (hq, hqw, over)
        if t == "optik":
            failed += over
    return failed


def norm_vcf(a, vcf, out, region, bed=None):
    """Normalize a VCF to `region` and return its {snp, indel} allele keys, inside `bed` when given."""
    # -r needs an index, which is always rebuilt, because the filtered VCF is rewritten on each run and a stale index from an earlier run makes norm return nothing.
    sh([a.bcftools, "index", "-f", vcf])
    cmd = [
        a.bcftools,
        "norm",
        "-m",
        "-both",
        "-f",
        a.ref,
        "-r",
        region,
        "-Oz",
        "-o",
        out,
        vcf,
    ]
    sh(cmd)
    sh([a.bcftools, "index", "-f", out])
    keys = {"snp": set(), "indel": set()}
    view = [a.bcftools, "view", "-H"] + (["-T", bed] if bed else []) + [out]
    p = subprocess.run(
        [str(c) for c in view], capture_output=True, text=True, check=True
    )
    for line in p.stdout.splitlines():
        f = line.split("\t")
        ref, alt = f[3], f[4]
        if alt in ("*", ".") or "<" in alt:
            continue
        kind = "snp" if len(ref) == len(alt) == 1 else "indel"
        keys[kind].add((f[0], int(f[1]), ref, alt))
    return keys


def vcfeval(a, truth_vcf, calls_vcf, bed, out):
    """Compare calls with the truth by RTG vcfeval, haplotype- and genotype-aware.

    vcfeval is the engine of GIAB's hap.py, so an indel or SNP written differently but giving the same haplotypes counts as a match.
    Returns (TP, FP, FN) per type, with recall taken from the baseline side and precision from the call side."""
    sdf = os.path.join(a.out, f"{a.region}.sdf")
    if not os.path.isdir(sdf):
        fa = os.path.join(a.out, f"{a.region}.fa")
        sh([a.samtools, "faidx", a.ref, a.region], out=fa)
        sh([a.rtg, "format", "-o", sdf, fa])
    for v in (truth_vcf, calls_vcf):
        sh([a.bcftools, "index", "-t", "-f", v])
    shutil.rmtree(out, ignore_errors=True)
    env = dict(os.environ, RTG_MEM=os.environ.get("RTG_MEM", "6g"))
    sh(
        [
            a.rtg,
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
            a.threads,
        ],
        env=env,
    )

    def count(name):
        """SNP and indel counts of one vcfeval output VCF."""
        c = {"snp": 0, "indel": 0}
        p = subprocess.run(
            [str(a.bcftools), "view", "-H", os.path.join(out, name)],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in p.stdout.splitlines():
            f = line.split("\t")
            ref, alts = f[3], [x for x in f[4].split(",") if x not in ("*", ".")]
            c[
                "snp"
                if alts and all(len(ref) == len(x) == 1 for x in alts)
                else "indel"
            ] += 1
        return c

    tpb, tpc, fp, fn = (
        count(n) for n in ("tp-baseline.vcf.gz", "tp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz")
    )
    return {k: (tpb[k], tpc[k], fp[k], fn[k]) for k in ("snp", "indel")}


def variants(a):
    """Variants mode, which calls germline variants from each tool's mapping and scores them against GIAB."""
    bed = os.path.join(a.out, "truth_region.bed")
    with open(a.truth_bed) as f, open(bed, "w") as o:
        for line in f:
            if line.split("\t")[0] == a.region:
                o.write(line)
    truth = norm_vcf(
        a, a.truth_vcf, os.path.join(a.out, "truth.norm.vcf.gz"), a.region, bed
    )
    print(
        f"truth in {a.region} benchmark regions: {len(truth['snp'])} SNPs, {len(truth['indel'])} indels; "
        f"scored by {'rtg vcfeval (haplotype- and genotype-aware)' if a.rtg else 'exact allele match'}"
    )
    failed, res = [], {}
    for t, sam in run_maps(a, "var", a.ref, a.r1, a.r2).items():
        bam = os.path.join(a.out, f"{t}.bam")
        if not os.path.exists(bam + ".bai"):
            sh([a.samtools, "sort", "-@", a.threads, "-m", "1G", "-o", bam, sam])
            sh([a.samtools, "index", bam])
        raw = os.path.join(a.out, f"{t}.calls.vcf.gz")
        if not os.path.exists(raw):
            mp = memcap.popen(
                [
                    str(x)
                    for x in (
                        a.bcftools,
                        "mpileup",
                        "-f",
                        a.ref,
                        "-r",
                        a.region,
                        "-a",
                        "AD,DP",
                        "-Ou",
                        bam,
                    )
                ],
                cap_bytes(),
                stdout=subprocess.PIPE,
            )
            cl = subprocess.run(
                [str(a.bcftools), "call", "-mv", "-Oz", "-o", raw],
                stdin=mp.stdout,
                check=False,
            )
            mp.wait()
            if cl.returncode or mp.returncode:
                sys.exit(f"bcftools failed for {t}")
        filt = os.path.join(a.out, f"{t}.filt.vcf.gz")
        sh([a.bcftools, "view", "-i", "QUAL>=20 && INFO/DP>=5", "-Oz", "-o", filt, raw])
        calls = norm_vcf(
            a, filt, os.path.join(a.out, f"{t}.norm.vcf.gz"), a.region, bed
        )
        if not calls["snp"]:
            sys.exit(
                f"{t}: no calls in {a.region} after filtering -- pipeline error, not a result"
            )
        res[t] = {}
        ve = (
            vcfeval(
                a,
                os.path.join(a.out, "truth.norm.vcf.gz"),
                os.path.join(a.out, f"{t}.norm.vcf.gz"),
                bed,
                os.path.join(a.out, f"vcfeval_{t}"),
            )
            if a.rtg
            else None
        )
        for kind in ("snp", "indel"):
            if ve:
                tp, tpc, fp, fn = ve[kind]
                prec = tpc / max(1, tpc + fp)
            else:
                # Without RTG, exact allele match after bcftools norm, which is stricter, because equivalent representations of one haplotype count as an FP and an FN.
                tp = len(calls[kind] & truth[kind])
                fp = len(calls[kind] - truth[kind])
                prec = tp / max(1, tp + fp)
                fn = len(truth[kind] - calls[kind])
            rec = tp / max(1, tp + fn)
            f1 = 2 * prec * rec / max(1e-9, prec + rec)
            res[t][kind] = (tp, fp, fn, prec, rec, f1)
            ROWS.append(
                {
                    "tool": t,
                    "type": kind,
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "precision": round(prec, 4),
                    "recall": round(rec, 4),
                    "F1": round(f1, 4),
                }
            )
            print(
                f"  {t:9} {kind:5} TP {tp:6} FP {fp:5} FN {fn:5}  precision {prec:.4f}  recall {rec:.4f}  F1 {f1:.4f}"
            )
    if "optik" in res:
        for other in res:
            if other == "optik":
                continue
            for kind in ("snp", "indel"):
                if res["optik"][kind][5] < res[other][kind][5] - 0.005:
                    failed.append(
                        f"{kind} F1: optik {res['optik'][kind][5]:.4f} < {other} {res[other][kind][5]:.4f} - 0.005"
                    )
    return failed


def longreads(a):
    """Longreads mode, which compares optik's and minimap2's primaries and large indels on real long reads."""
    sams = {}
    for t, cmd in (
        ("optik", [a.optik, "map", "-r", a.ref, "-q", a.reads, "-t", a.threads]),
        ("minimap2", [a.minimap2, "-ax", a.preset, "-t", a.threads, a.ref, a.reads]),
    ):
        sam = os.path.join(a.out, f"{t}_long.sam")
        if not os.path.exists(sam) or os.path.getsize(sam) == 0:
            print(f"  mapping long reads with {t}", flush=True)
            sh(cmd, out=sam + ".part")
            os.replace(sam + ".part", sam)
        sams[t] = sam
    recs, big = {}, {}
    for t, sam in sams.items():
        d, b = {}, []
        with open(sam) as fh:
            for line in fh:
                if line[0] == "@":
                    continue
                f = line.split("\t", 7)
                flag = int(f[1])
                if flag & SAM_FLAGS_SECONDARY_SUPPLEMENTARY:
                    continue
                if flag & SAM_FLAG_UNMAPPED:
                    d[f[0]] = None
                    continue
                r, s, o = origin(f)
                d[f[0]] = (r, s, o, int(f[4]))
                pos = int(f[3]) - 1
                for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
                    n = int(n)
                    if op in "ID" and n >= 50:
                        b.append((r, pos, op, n, f[0]))
                    if op in "MDN=X":
                        pos += n
        recs[t], big[t] = d, b
    names = set(recs["optik"]) | set(recs["minimap2"])
    agree = dis = dis_hq_o = dis_hq_m = un_o = un_m = 0
    for nm in names:
        o, m = recs["optik"].get(nm), recs["minimap2"].get(nm)
        un_o += o is None
        un_m += m is None
        if o is None or m is None:
            continue
        if o[0] == m[0] and o[1] == m[1] and abs(o[2] - m[2]) <= max(50, 0):
            agree += 1
        else:
            dis += 1
            dis_hq_o += o[3] >= HIGH_MAPQ
            dis_hq_m += m[3] >= HIGH_MAPQ
    print(
        f"  reads {len(names)}: unmapped optik {un_o}, minimap2 {un_m}; same primary locus {agree}, "
        f"different {dis} (MAPQ>=30 there: optik {dis_hq_o}, minimap2 {dis_hq_m})"
    )

    # Large indels are either shared (same read, position within 100 bp, same operation, length within 20%) or tool-only.

    by_read = {t: collections.defaultdict(list) for t in big}
    for t, b in big.items():
        for r, p, op, n, nm in b:
            by_read[t][nm].append((r, p, op, n))
    only = {}
    for t, u in (("optik", "minimap2"), ("minimap2", "optik")):
        c = 0
        for nm, evs in by_read[t].items():
            for r, p, op, n in evs:
                if not any(
                    r == r2
                    and op == op2
                    and abs(p - p2) <= 100
                    and abs(n - n2) <= max(10, 0.2 * n)
                    for r2, p2, op2, n2 in by_read[u].get(nm, [])
                ):
                    c += 1
        only[t] = c
    print(
        f"  >= 50 bp indels in primaries: optik {len(big['optik'])} ({only['optik']} not in minimap2's record "
        f"of the same read), minimap2 {len(big['minimap2'])} ({only['minimap2']} not in optik's)"
    )
    for t, u in (("optik", "minimap2"), ("minimap2", "optik")):
        ROWS.append(
            {
                "tool": t,
                "reads": len(names),
                "unmapped": un_o if t == "optik" else un_m,
                "same_primary_locus": agree,
                "different_locus": dis,
                "mapq30_at_different": dis_hq_o if t == "optik" else dis_hq_m,
                "indels_50bp": len(big[t]),
                f"indels_50bp_not_in_{u}": only[t],
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
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("mode", choices=("anchored", "variants", "longreads"))
    ap.add_argument("--ref", required=True)
    ap.add_argument("--r1")
    ap.add_argument("--r2")
    ap.add_argument("--reads")
    ap.add_argument("--preset", default="map-hifi")
    ap.add_argument(
        "--pairs", type=int, default=0, help="anchored: use the first N pairs (0 = all)"
    )
    ap.add_argument("--truth-vcf")
    ap.add_argument("--truth-bed")
    ap.add_argument("--region", default="chr20")
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--threads",
        default=str(hw.cpus()),
        help="threads of every tool (default: the usable CPUs)",
    )
    ap.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    ap.add_argument("--minimap2", default=os.path.join(ROOT, "temp/minimap2/minimap2"))
    ap.add_argument(
        "--bwa",
        default=shutil.which("bwa")
        or next(
            (
                p
                for p in (
                    os.path.join(ROOT, "temp/bwa/bwa"),
                    os.path.join(ROOT, "temp/gcroots/bwa/bin/bwa"),
                )
                if os.path.exists(p)
            ),
            os.path.join(ROOT, "temp/bwa/bwa"),
        ),
    )
    ap.add_argument("--bwa-index", help="BWA index prefix (adds BWA-MEM as a baseline)")
    ap.add_argument("--samtools", default=shutil.which("samtools"))
    ap.add_argument(
        "--bcftools",
        default=os.environ.get("BCFTOOLS_BIN")
        or shutil.which("bcftools")
        or (
            os.path.join(ROOT, "temp/bcftools/bcftools")
            if os.path.exists(os.path.join(ROOT, "temp/bcftools/bcftools"))
            else os.path.join(ROOT, "temp/gcroots/bcftools/bin/bcftools")
        ),
    )
    ap.add_argument(
        "--rtg",
        default=shutil.which("rtg")
        or next(
            (
                p
                for p in (
                    os.path.join(ROOT, "temp/rtg/rtg"),
                    os.path.join(ROOT, "temp/gcroots/rtg/bin/rtg"),
                )
                if os.path.exists(p)
            ),
            None,
        ),
        help="RTG Tools rtg (vcfeval: haplotype-aware comparison, as GIAB hap.py); "
        "without it calls are matched allele by allele",
    )
    ap.add_argument("--gate", action="store_true")
    ap.add_argument(
        "--mem-cap",
        default=hw.default_cap_str(),
        help="resident-memory cap of each tool ('' = uncapped; default: usable memory less 2 GiB)",
    )
    results_file.add_argument(ap)
    a = ap.parse_args()
    global MEM_CAP
    MEM_CAP = a.mem_cap or None
    os.makedirs(a.out, exist_ok=True)
    failed = {"anchored": anchored, "variants": variants, "longreads": longreads}[
        a.mode
    ](a)
    for x in failed:
        print("  [FAIL]", x)
    print(f"real_reads {a.mode}: " + ("FAILED" if failed else "PASS"))
    results_file.write_section(
        a.report,
        f"real_{a.mode}",
        MODE_TITLES[a.mode] + f" (reference {os.path.basename(a.ref)})",
        [("", ROWS)],
        optik=a.optik,
        notes=[f"- FAIL {x}" for x in failed]
        + [f"SAMs and logs: `{os.path.relpath(a.out, ROOT)}`"],
    )
    return 1 if failed and a.gate else 0


if __name__ == "__main__":
    sys.exit(main())
