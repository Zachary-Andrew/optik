#!/usr/bin/env python3
"""pipeline_parity.py - checks that optik's outputs serve the pipelines the reference tools serve, with the same results.

Totals and k-mer sets are checked elsewhere (tool_bench, validate_*).
Here each output is consumed the way a real analysis consumes the reference tool's output, and the downstream results are compared.

  1. k-mer dictionary, per k-mer, as used for read screening or classification by k-mer hits.
     build-spss, SBWT (--add-reverse-complements) and FMSI (KmerCamel -S masked superstring) index one set of Cuttlefish unitigs.
     Each window of every read must get the oracle's present or absent answer from all three, so per-read hit counts agree.
     optik --ids must also behave as a dictionary, with one id per canonical k-mer that is the same for both strands.
  2. Compacted dBG as a graph.
     gfatools parses the GFA of optik and of Bifrost and reports its segments and links.
     The input genome must be spelled by a walk that steps only along the graph's L lines, the property graph tools (read-to-graph mappers, bubble finders, Bandage) rely on.
     optik's unitigs must also be interchangeable with Cuttlefish's as the input of a downstream index, i.e. SSHash built on either answers the same reads identically.
  3. Mapping followed by variant calling, the main consumer of a mapper's SAM.
     A sample genome with planted SNPs and indels yields 30x reads, which optik map and minimap2 -ax sr map.
     samtools sort, index, flagstat and depth then run, followed by bcftools mpileup | call -mv --ploidy 1 | norm, and each VCF is compared with the truth set.
     The two call sets are also compared with each other.

The oracles are the external tools themselves and independent Python that uses no optik code.
A single failed check fails the run.

Usage
  python3 bench/pipeline_parity.py --genome temp/sshash/data/queries/ecoli1.fasta.gz \
      --workdir /tmp/pp [--bcftools PATH] [--threads 8]
"""

import argparse
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import validate_unitigs
from config import PATHS
from diagnostics import DETAIL_CHARS
from dna import read_fasta, write_fasta
from sam_spec import (
    HIGH_MAPQ,
    GfaLinkCol,
    GfaSegmentCol,
    SamCol,
    SamFlag,
    VcfCol,
)

# optik query --ids value of an invalid (non-ACGT) window
OPTIK_ABSENT_ID = -1
# value of a gfatools stat line that is missing, which never matches
MISSING_STAT = -1
# exact-k runs record: random flanks of EXACT_RUN_FLANK bases around the runs, and a tail of EXACT_RUN_TAIL
EXACT_RUN_SEED = 12
EXACT_RUN_FLANK = 300
EXACT_RUN_TAIL = 200
# SSHash uses minimizer length SSHASH_LONG_K_M from k = SSHASH_LONG_K_MIN on, and otherwise k / 2 + 2 but at least SSHASH_MIN_M.
SSHASH_LONG_K_MIN = 25
SSHASH_LONG_K_M = 17
SSHASH_MIN_M = 5
# Below SSHASH_LONG_K_MIN, m is k / 2 plus two.
SSHASH_K_PER_M = 2
SSHASH_M_EXTRA = 2
# Variants are planted about once every SNP_EVERY bases, and one in INDEL_EVERY is an indel of up to MAX_INDEL_LEN bases.
# Each keeps VARIANT_EDGE_MARGIN bases away from the record ends and VARIANT_CLEAR_FLANK N-free bases on either side.
PLANT_SEED = 7
SNP_EVERY = 2000
INDEL_EVERY = 15000
VARIANT_EDGE_MARGIN = 200
VARIANT_CLEAR_FLANK = 60
MAX_INDEL_LEN = 10
INSERTION_PROB = 0.5
# simulated pairs: DEFAULT_COVERAGE x the sample genome in PAIR_READ_LEN bp mates
DEFAULT_COVERAGE = 30
PAIR_READ_LEN = 150
MATES_PER_PAIR = 2
# breadth is reported as the share of bases at least this deep
MIN_BREADTH_DEPTH = 10
# Picard errors listed for optik (gated) and for minimap2 (informational)
MAX_ERRORS_SHOWN = 5
INFO_ERRORS_SHOWN = 3
# leading records of a BAM whose QUAL is checked
QUAL_CHECK_RECORDS = 200
# the standard hard filter on calls
FILTER_MIN_QUAL = 20
FILTER_MIN_DP = 5
DP_PREFIX = "DP="
# indel sites whose reads' MAPQs are collected
MAX_INDEL_SITES = 200
# parity tolerances against minimap2: recall and precision, shared true variants, MAPQ >= 30 share of indel reads, mean depth (relative), breadth, properly paired (share of reads)
CALL_TOL = 0.005
MIN_SHARED_TRUE = 0.995
INDEL_MAPQ_TOL = 0.01
DEPTH_REL_TOL = 0.02
BREADTH_TOL = 0.005
PROPER_PAIR_TOL = 0.005
# foreign reads appended to the shared read set
FOREIGN_SEED = 9
FOREIGN_READS = 2000
FOREIGN_READ_LEN = 150
# CLI defaults
DEFAULT_K = 31
DEFAULT_THREADS = 8

RESULTS = {}
FAILS = []


# A GFA link stands for an arc in each direction, which gfatools counts separately.
ARCS_PER_LINK = 2
# optik's unitig FASTA puts each record on a header line and one sequence line.
FASTA_LINES_PER_RECORD = 2


def is_snp(variant):
    """Whether a (chrom, pos, ref, alt) variant replaces one base with another."""
    _chrom, _pos, ref, alt = variant
    return len(ref) == len(alt) == 1


def is_indel(variant):
    """Whether a (chrom, pos, ref, alt) variant inserts or deletes bases."""
    _chrom, _pos, ref, alt = variant
    return len(ref) != len(alt)


def check(name, passed, **detail):
    """Record and print one check result, remembering failures in FAILS."""
    RESULTS[name] = {"pass": bool(passed), **detail}
    print(
        f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""),
        flush=True,
    )
    if not passed:
        FAILS.append(name)


def run(cmd, out=None, **kwargs):
    cmd = [str(arg) for arg in cmd]
    if out is not None:
        with open(out, "wb") as handle:
            proc = subprocess.run(
                cmd, stdout=handle, stderr=subprocess.PIPE, **kwargs, check=False
            )
    else:
        proc = subprocess.run(cmd, capture_output=True, **kwargs, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr.decode()[-DETAIL_CHARS:]}"
        )
    return proc


# Section 1, the k-mer dictionary
def kmer_parity(work_dir, genome, unitigs, reads, kmer_len, threads):
    """Section 1, which compares per-k-mer answers of optik, SBWT, FMSI and an oracle on the same reads."""
    print("== 1. k-mer dictionary, per k-mer (optik SPSS vs SBWT vs FMSI vs oracle)")
    optik, sbwt_bin, kmercamel, fmsi_bin = (
        PATHS["optik"],
        PATHS["sbwt"],
        PATHS["kmercamel"],
        PATHS["fmsi"],
    )
    recs = list(read_fasta(reads))
    # The oracle holds the canonical k-mers of the unitigs, computed in independent Python.
    kset = set()
    for _, seq in read_fasta(unitigs):
        for packed in validate_unitigs.windows(seq, kmer_len):
            kset.add(validate_unitigs.canon(packed, kmer_len))

    def oracle_vec(seq):
        """Presence of each k-mer window of s in the oracle set, None for windows with non-ACGT bases."""
        out = []
        for i in range(len(seq) - kmer_len + 1):
            window = seq[i : i + kmer_len].upper()
            if set(window) - set("ACGT"):
                out.append(None)  # invalid window
            else:
                out.append(validate_unitigs.canon(validate_unitigs.pack(window), kmer_len) in kset)
        return out

    truth = [oracle_vec(seq) for _, seq in recs]
    # optik --ids
    idx = work_dir / "o.spss"
    run([optik, "build-spss", "-i", unitigs, "-k", kmer_len, "-t", threads, "-o", idx])
    ids_path = work_dir / "o.ids"
    run([optik, "query", "-i", idx, "-q", reads, "--ids", ids_path, "-t", threads])
    optik_ids = [
        list(map(int, text_line.split()))
        for text_line in Path(ids_path).read_text().split("\n")[: len(recs)]
    ]
    # SBWT indexes both strands, as its documentation prescribes for canonical lookups.
    sbi = work_dir / "x.sbwt"
    run(
        [
            sbwt_bin,
            "build",
            "-i",
            unitigs,
            "-o",
            sbi,
            "-k",
            kmer_len,
            "--add-reverse-complements",
            "--temp-dir",
            work_dir,
            "-t",
            threads,
        ]
    )
    sbo = work_dir / "sbwt.out"
    run([sbwt_bin, "search", "-i", sbi, "-q", reads, "-o", sbo])
    sbwt = [
        list(map(int, text_line.split()))
        for text_line in Path(sbo).read_text().split("\n")[: len(recs)]
    ]
    # FMSI over KmerCamel's masked superstring (-S, since the input is unitigs)
    ms_fa = work_dir / "ms.fa"
    run([kmercamel, "compute", "-k", kmer_len, "-S", "-M", ms_fa, "-o", "/dev/null", unitigs])
    run([fmsi_bin, "index", ms_fa])
    fmo = work_dir / "fmsi.out"
    run([fmsi_bin, "query", "-q", reads, ms_fa], out=fmo)
    fmsi = [
        text_line.split("\t")[1] if "\t" in text_line else ""
        for text_line in Path(fmo).read_text().split("\n")[: len(recs)]
    ]

    n_win = n_valid = 0
    bad = Counter()
    hits = {"oracle": [], "optik": [], "sbwt": [], "fmsi": []}
    for r, truth_vec in enumerate(truth):
        hits_optik = hits_sbwt = hits_fmsi = hits_oracle = 0
        for i, want in enumerate(truth_vec):
            n_win += 1
            if want is None:
                # An invalid window is -1 in optik and absent in SBWT and FMSI.
                if optik_ids[r][i] != OPTIK_ABSENT_ID:
                    bad["optik invalid window"] += 1
                continue
            n_valid += 1
            in_optik = optik_ids[r][i] >= 0
            in_sbwt = sbwt[r][i] >= 0
            # FMSI writes one character per window. A short row leaves None, which counts as a mismatch.
            in_fmsi = fmsi[r][i] == "1" if i < len(fmsi[r]) else None
            if in_optik != want:
                bad["optik"] += 1
            if in_sbwt != want:
                bad["sbwt"] += 1
            if in_fmsi != want:
                bad["fmsi"] += 1
            hits_oracle += want
            hits_optik += in_optik
            hits_sbwt += in_sbwt
            hits_fmsi += bool(in_fmsi)
        hits["oracle"].append(hits_oracle)
        hits["optik"].append(hits_optik)
        hits["sbwt"].append(hits_sbwt)
        hits["fmsi"].append(hits_fmsi)
    check(
        "every window: optik == oracle",
        bad["optik"] == 0 and bad["optik invalid window"] == 0,
        windows=n_win,
        valid=n_valid,
        mismatches=bad["optik"],
    )
    check(
        "every window: SBWT == oracle (reference tool agrees)",
        bad["sbwt"] == 0,
        mismatches=bad["sbwt"],
    )
    check(
        "every window: FMSI == oracle (reference tool agrees)",
        bad["fmsi"] == 0,
        mismatches=bad["fmsi"],
    )
    # Read screening only looks at per-read hit counts, so those must match across all three tools as well.
    same = sum(
        optik_hits == sbwt_hits == fmsi_hits
        for optik_hits, sbwt_hits, fmsi_hits in zip(hits["optik"], hits["sbwt"], hits["fmsi"], strict=False)
    )
    check(
        "per-read k-mer hit counts identical (optik, SBWT, FMSI) -- what read screening uses",
        same == len(recs),
        reads=len(recs),
        identical=same,
    )
    # Dictionary semantics, one id per canonical k-mer, independent of strand
    by_kmer, passed = {}, True
    for (_name, seq), row in zip(recs, optik_ids, strict=False):
        for i, kmer_id in enumerate(row):
            if kmer_id < 0:
                continue
            canon_kmer = validate_unitigs.canon(validate_unitigs.pack(seq[i : i + kmer_len].upper()), kmer_len)
            passed &= by_kmer.setdefault(canon_kmer, kmer_id) == kmer_id
    ids = set(by_kmer.values())
    check(
        "optik --ids is a dictionary: one id per canonical k-mer (both strands), distinct k-mers distinct ids",
        passed and len(ids) == len(by_kmer),
        kmers=len(by_kmer),
    )


# Section 2, the compacted dBG as a graph
def mirror_edge(edge):
    """The same bidirected edge (from, orient, to, orient) read from its other end."""
    flip = {"+": "-", "-": "+"}
    from_seg, from_orient, to_seg, to_orient = edge
    return (to_seg, flip[to_orient], from_seg, flip[from_orient])


def parse_gfa(path):
    """Segments {name: sequence} and links {(from, orient, to, orient)} of a GFA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    segs, links = {}, set()
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                segs[fields[GfaSegmentCol.NAME]] = fields[GfaSegmentCol.SEQ].upper()
            elif fields[0] == "L":
                links.add(
                    (
                        fields[GfaLinkCol.FROM],
                        fields[GfaLinkCol.FROM_ORIENT],
                        fields[GfaLinkCol.TO],
                        fields[GfaLinkCol.TO_ORIENT],
                    )
                )
    return segs, links


def genome_walk(segs, links, genome, kmer_len):
    """Walk each N-free run of the genome k-mer by k-mer through the graph.

    Every k-mer must lie in a segment, and consecutive k-mers must either be consecutive in one segment or cross a segment boundary along an L line.
    Returns (k-mers walked, missing k-mers, steps without an edge)."""
    loc = {}  # canonical k-mer -> segment, offset, forward flag
    for sid, seq in segs.items():
        for i, packed in enumerate(validate_unitigs.windows(seq, kmer_len)):
            rev = validate_unitigs.revcomp_packed(packed, kmer_len)
            loc[min(packed, rev)] = (sid, i, packed <= rev)
    flip = {"+": "-", "-": "+"}
    link_set = set(links) | {(to_seg, flip[to_orient], from_seg, flip[from_orient]) for from_seg, from_orient, to_seg, to_orient in links}
    walked = missing = noedge = 0
    for _, genome_seq in read_fasta(genome):
        prev = None
        for run_ in genome_seq.upper().replace("N", " ").split():
            prev = None
            for packed in validate_unitigs.windows(run_, kmer_len):
                walked += 1
                rev = validate_unitigs.revcomp_packed(packed, kmer_len)
                canon_kmer = min(packed, rev)
                if canon_kmer not in loc:
                    missing += 1
                    prev = None
                    continue
                sid, off, fwd = loc[canon_kmer]
                # Orientation of the segment as this k-mer reads it
                orient = "+" if (packed == canon_kmer) == fwd else "-"
                n_k = len(segs[sid]) - kmer_len + 1
                pos = (
                    off if orient == "+" else n_k - 1 - off
                )  # position along the walking direction
                cur = (sid, orient, pos, n_k)
                if prev is not None:
                    prev_sid, prev_orient, prev_pos, prev_len = prev
                    inside = prev_sid == sid and prev_orient == orient and pos == prev_pos + 1
                    across = prev_pos == prev_len - 1 and pos == 0 and (prev_sid, prev_orient, sid, orient) in link_set
                    if not (inside or across):
                        noedge += 1
                prev = cur
    return walked, missing, noedge


def graph_usability(work_dir, genome, cf_unitigs, reads, kmer_len, threads, gfatools):
    """Section 2, which checks optik's GFA as graph tools read it and its unitigs as index input."""
    print(
        "== 2. compacted dBG as a graph (optik GFA vs Bifrost GFA; unitigs vs Cuttlefish)"
    )
    optik, bifrost, sshash = PATHS["optik"], PATHS["bifrost"], PATHS["sshash"]
    optik_gfa = work_dir / "optik.gfa"
    run([optik, "dbruijn", "-i", genome, "-k", kmer_len, "-c", "1", "-t", threads, "-o", optik_gfa])
    bfp = work_dir / "bf"
    run([bifrost, "build", "-r", genome, "-k", kmer_len, "-t", threads, "-o", bfp])
    bifrost_gfa = next(
        path for path in (Path(str(bfp) + ".gfa.gz"), Path(str(bfp) + ".gfa")) if path.exists()
    )
    for name, gfa_file in (("optik", optik_gfa), ("bifrost", bifrost_gfa)):
        src = gfa_file
        if str(gfa_file).endswith(".gz"):
            src = work_dir / f"{name}_plain.gfa"
            with gzip.open(gfa_file, "rb") as in_handle, open(src, "wb") as out_handle:
                shutil.copyfileobj(in_handle, out_handle)
        proc = subprocess.run(
            [gfatools, "stat", src], capture_output=True, text=True, check=False
        )
        stat = dict(
            text_line.split(": ", 1)
            for text_line in proc.stdout.splitlines()
            if ": " in text_line
        )
        segs, links = parse_gfa(gfa_file)
        # GFA writers follow different conventions.
        # optik writes each bidirected edge once and Bifrost writes both orientations.
        # gfatools merges the two orientations and counts links as distinct edges minus the self-mirrored ones that join a segment to its own reverse, and arcs as twice the edges minus the self-mirrored ones.
        edges = {min(edge, mirror_edge(edge)) for edge in links}
        selfm = sum(1 for edge in edges if edge == mirror_edge(edge))
        check(
            f"gfatools parses the {name} GFA (segments, edges, arcs and length as written)",
            proc.returncode == 0
            and int(stat.get("Number of segments", MISSING_STAT)) == len(segs)
            and int(stat.get("Number of links", MISSING_STAT)) == len(edges) - selfm
            and int(stat.get("Number of arcs", MISSING_STAT)) == ARCS_PER_LINK * len(edges) - selfm
            and int(stat.get("Total segment length", MISSING_STAT))
            == sum(map(len, segs.values())),
            segments=stat.get("Number of segments"),
            links=stat.get("Number of links"),
            l_lines=len(links),
            edges=len(edges),
            total_len=stat.get("Total segment length"),
        )
        walked, missing, noedge = genome_walk(segs, links, genome, kmer_len)
        check(
            f"{name}: the genome is a walk along the graph's L lines (every k-mer present, every step an edge)",
            missing == 0 and noedge == 0,
            kmers=walked,
            missing=missing,
            steps_without_edge=noedge,
        )
        RESULTS.setdefault("graph_shape", {})[name] = {
            "segments": len(segs),
            "links": len(links),
        }
    # The data itself is checked against independent builders.
    # The unitig sets must equal Cuttlefish's and Bifrost's as sets of orientation-free sequences, the edges encoded by optik's L lines must equal Bifrost's (each link taken as its canonical (k+1)-mer), and every KC tag must be the true summed k-mer abundance in the genome.
    comp = str.maketrans("ACGT", "TGCA")

    def rcs(seq):
        """Reverse complement of x."""
        return seq.translate(comp)[::-1]

    def canon(seq):
        """Canonical form of k-mer string x."""
        return min(seq, rcs(seq))

    oseg, olinks = parse_gfa(optik_gfa)
    bsrc = bifrost_gfa if not str(bifrost_gfa).endswith(".gz") else work_dir / "bifrost_plain.gfa"
    bseg, blinks = parse_gfa(bsrc)
    cfs = {canon(seq) for _, seq in read_fasta(cf_unitigs)}
    os_ = {canon(seq) for seq in oseg.values()}
    check(
        "unitigs identical to Cuttlefish's (orientation-free sequence sets; both edge-centric)",
        os_ == cfs and len(os_) == len(oseg),
        optik=len(os_),
        cuttlefish=len(cfs),
        only_optik=len(os_ - cfs),
        only_cf=len(cfs - os_),
    )

    def edge_set(segs, links):
        """The canonical (k+1)-mer of each GFA link."""
        out = set()
        for from_seg, from_orient, to_seg, to_orient in links:
            from_seq = segs[from_seg] if from_orient == "+" else rcs(segs[from_seg])
            to_seq = segs[to_seg] if to_orient == "+" else rcs(segs[to_seg])
            out.add(canon(from_seq[-kmer_len:] + to_seq[kmer_len - 1]))
        return out

    gk1 = set()  # the genome's (k+1)-mers, i.e. its edges
    for _, genome_seq in read_fasta(genome):
        genome_seq = genome_seq.upper()
        for i in range(len(genome_seq) - kmer_len):
            kmer = genome_seq[i : i + kmer_len + 1]
            if "N" not in kmer:
                gk1.add(canon(kmer))
    internal = set()
    for seq in oseg.values():
        for i in range(len(seq) - kmer_len):
            internal.add(canon(seq[i : i + kmer_len + 1]))
    optik_edges, bifrost_edges = edge_set(oseg, olinks), edge_set(bseg, blinks)
    check(
        "links: every L line is a genome (k+1)-mer, and links + unitig-internal edges = all of them",
        optik_edges <= gk1 and (optik_edges | internal) == gk1 and not (optik_edges & internal),
        links=len(optik_edges),
        not_in_genome=len(optik_edges - gk1),
        genome_edges_missing=len(gk1 - optik_edges - internal),
    )
    # Bifrost is node-centric by design.
    # It links two k-mers that overlap by k-1 even when that (k+1)-mer never occurs, and such branches split unitigs, so the edges at the splits become links as well.
    extra = bifrost_edges - optik_edges
    check(
        "Bifrost's extra links: unobserved k-1 overlaps, or edges inside optik's unitigs that its branches split",
        optik_edges <= bifrost_edges and extra <= ((extra - gk1) | internal),
        bifrost_only=len(extra),
        unobserved=len(extra - gk1),
        split_internal=len(extra & internal),
        unexplained=len(extra - (extra - gk1) - internal),
        optik_only=len(optik_edges - bifrost_edges),
    )
    del gk1, internal
    cnt = {}
    for _, genome_seq in read_fasta(genome):
        genome_seq = genome_seq.upper()
        for i in range(len(genome_seq) - kmer_len + 1):
            kmer = genome_seq[i : i + kmer_len]
            if "N" in kmer:
                continue
            canon_kmer = canon(kmer)
            cnt[canon_kmer] = cnt.get(canon_kmer, 0) + 1
    bad_kc = bad_km = 0
    with open(optik_gfa) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != "S":
                continue
            seq = fields[GfaSegmentCol.SEQ]
            want = sum(cnt.get(canon(seq[i : i + kmer_len]), 0) for i in range(len(seq) - kmer_len + 1))
            tags = dict(tag.split(":", 1) for tag in fields[GfaSegmentCol.TAGS:])
            kc_value = int(tags["KC"].split(":")[1])
            bad_kc += kc_value != want
            bad_km += tags["km"].split(":")[1] != f"{want / (len(seq) - kmer_len + 1):.2f}"
    check(
        "KC = the summed genome abundance of each segment's k-mers, km = KC / k-mers",
        bad_kc == 0 and bad_km == 0,
        wrong_kc=bad_kc,
        wrong_km=bad_km,
    )
    del cnt
    # A run of exactly k bases, between N runs or as a whole record, holds one k-mer.
    # Cuttlefish 3.0.3 and 3.1.0 drop such k-mers, and optik must keep them.
    rng = random.Random(EXACT_RUN_SEED)

    def rnd(length):
        """Random ACGT string of length n."""
        return "".join(rng.choice("ACGT") for _ in range(length))

    exact_runs_fa = work_dir / "exact_k_runs.fa"
    write_fasta(
        exact_runs_fa,
        [
            (
                "x",
                rnd(EXACT_RUN_FLANK)
                + "NNN"
                + rnd(kmer_len)
                + "NN"
                + rnd(EXACT_RUN_FLANK)
                + "N"
                + rnd(kmer_len + 1)
                + "NN"
                + rnd(EXACT_RUN_TAIL),
            ),
            ("y", rnd(kmer_len)),
        ],
    )
    truth_k = set()
    for _, seq in read_fasta(exact_runs_fa):
        for runs in seq.split("N"):
            for i in range(len(runs) - kmer_len + 1):
                truth_k.add(canon(runs[i : i + kmer_len]))

    def kset(path):
        """Canonical k-mers of every record in a FASTA file."""
        out = set()
        for _, seq in read_fasta(path):
            for i in range(len(seq) - kmer_len + 1):
                out.add(canon(seq[i : i + kmer_len]))
        return out

    for mode in ("part", "table"):
        nro = work_dir / f"exact_k_{mode}.fa"
        run(
            [optik, "dbruijn", "-i", exact_runs_fa, "-k", kmer_len, "-c", "1", "-o", nro, "--graph", mode]
        )
        check(
            f"k-mers of runs exactly k long are kept (--graph {mode})",
            kset(nro) == truth_k,
            truth=len(truth_k),
            got=len(kset(nro)),
        )
    ncf = work_dir / "exact_k_cf"
    run(
        [
            PATHS["cuttlefish"],
            "build",
            "--ref",
            "-s",
            exact_runs_fa,
            "-k",
            kmer_len,
            "-t",
            2,
            "-o",
            ncf,
            "-w",
            work_dir,
        ]
    )
    # Cuttlefish is only reported here, not checked. Its handling of exact-k runs is its own choice.
    print(
        f"  [info] Cuttlefish on the same input: {len(truth_k - kset(str(ncf) + '.fa'))} of "
        f"{len(truth_k)} k-mers missing (runs of exactly k bases)"
    )
    # Partitioned compaction, the default here, and the table walker must agree.
    table_gfa = work_dir / "optik_table.gfa"
    run(
        [
            optik,
            "dbruijn",
            "-i",
            genome,
            "-k",
            kmer_len,
            "-c",
            "1",
            "-t",
            threads,
            "-o",
            table_gfa,
            "--graph",
            "table",
        ]
    )
    check(
        "dbruijn --graph table writes the same GFA as the default (partitioned) build",
        Path(table_gfa).read_bytes() == Path(optik_gfa).read_bytes(),
    )
    # optik's unitigs serve as the input of a downstream index in place of Cuttlefish's, taken from dbruijn's own FASTA output (-o *.fa) as a pipeline would use it.
    # These are the same sequences as the GFA's segments.
    segs, _ = parse_gfa(optik_gfa)
    optik_unitigs = work_dir / "optik_unitigs.fa"
    run([optik, "dbruijn", "-i", genome, "-k", kmer_len, "-c", "1", "-t", threads, "-o", optik_unitigs])
    lines = Path(optik_unitigs).read_text().splitlines()
    heads, seqs = lines[0::FASTA_LINES_PER_RECORD], lines[1::FASTA_LINES_PER_RECORD]
    gfa_tags = {}
    with open(optik_gfa) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                gfa_tags[fields[GfaSegmentCol.NAME]] = fields[GfaSegmentCol.TAGS:]
    check(
        "dbruijn -o *.fa: unitig FASTA with the GFA's segments (ids, sequences, LN/KC/km)",
        len(heads) == len(segs)
        and all(head.startswith(">") for head in heads)
        and [(head[1:].split()[0], seq) for head, seq in zip(heads, seqs, strict=False)]
        == list(segs.items())
        and all(head[1:].split()[1:] == gfa_tags[head[1:].split()[0]] for head in heads),
        records=len(heads),
        segments=len(segs),
    )
    res = {}
    for name, unitigs_path in (("optik", optik_unitigs), ("cuttlefish", cf_unitigs)):
        sshash_idx = work_dir / f"{name}.sshash"
        run(
            [
                sshash,
                "build",
                "-i",
                unitigs_path,
                "-k",
                kmer_len,
                "-m",
                str(SSHASH_LONG_K_M) if kmer_len >= SSHASH_LONG_K_MIN else str(max(SSHASH_MIN_M, kmer_len // SSHASH_K_PER_M + SSHASH_M_EXTRA)),
                "-t",
                threads,
                "-o",
                sshash_idx,
            ]
        )
        proc = run([sshash, "query", "-i", sshash_idx, "-q", reads])
        txt = proc.stdout.decode() + proc.stderr.decode()
        res[name] = {
            key: value
            for key, value in (
                text_line.split(" = ")
                for text_line in txt.splitlines()
                if " = " in text_line
            )
            if key in ("num_kmers", "num_positive_kmers", "num_negative_kmers")
        }
    check(
        "SSHash built on optik's unitigs answers the reads exactly as on Cuttlefish's",
        res["optik"] == res["cuttlefish"] and bool(res["optik"]),
        optik=res["optik"],
        cuttlefish=res["cuttlefish"],
    )


# Section 3, mapping followed by variant calling
def plant_variants(ref_recs, seed, snp_every=SNP_EVERY, indel_every=INDEL_EVERY):
    """Make a sample genome from the reference with SNPs and 1-10 bp indels at spaced, N-free positions.

    Returns (sample records, truth rows of (chrom, pos1, ref, alt))."""
    rnd = random.Random(seed)
    sample, truth = [], []
    for name, seq in ref_recs:
        out, last, pos = [], 0, VARIANT_EDGE_MARGIN
        while pos < len(seq) - VARIANT_EDGE_MARGIN:
            step = rnd.randint(snp_every // 2, snp_every + snp_every // 2)
            pos += step
            if pos >= len(seq) - VARIANT_EDGE_MARGIN:
                break
            win = seq[pos - VARIANT_CLEAR_FLANK : pos + VARIANT_CLEAR_FLANK]
            if "N" in win.upper():
                continue
            out.append(seq[last:pos])
            ref_base = seq[pos].upper()
            if rnd.random() < snp_every / indel_every:  # an indel
                indel_len = rnd.randint(1, MAX_INDEL_LEN)
                if rnd.random() < INSERTION_PROB:  # an insertion after base i
                    ins = "".join(rnd.choice("ACGT") for _ in range(indel_len))
                    out.append(seq[pos] + ins)
                    truth.append((name, pos + 1, ref_base, ref_base + ins))
                    last = pos + 1
                else:  # a deletion of ln bases after base i
                    out.append(seq[pos])
                    truth.append((name, pos + 1, seq[pos : pos + 1 + indel_len].upper(), ref_base))
                    last = pos + 1 + indel_len
            else:
                alt = rnd.choice([base for base in "ACGT" if base != ref_base])
                out.append(alt)
                truth.append((name, pos + 1, ref_base, alt))
                last = pos + 1
            pos = last + VARIANT_CLEAR_FLANK
        out.append(seq[last:])
        sample.append((name, "".join(out)))
    return sample, truth


def write_vcf(path, truth, ref_recs):
    """Write the truth variants as a minimal VCF 4.2 with contig lines for ref_recs."""
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.writelines(f"##contig=<ID={name},length={len(seq)}>\n" for name, seq in ref_recs)
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.writelines(f"{chrom}\t{pos}\t.\t{ref_allele}\t{alt_allele}\t.\tPASS\t.\n" for chrom, pos, ref_allele, alt_allele in truth)


def vcf_set(path, bcftools, ref, min_qual=None, min_dp=None):
    """(chrom, pos, ref, alt) of every ALT allele after left-normalization.

    Optionally only records with QUAL >= min_qual and INFO/DP >= min_dp are kept, the hard filter that variant-calling pipelines apply."""
    proc = run([bcftools, "norm", "-f", ref, "-m", "-any", "-Ov", path])
    out = set()
    for text_line in proc.stdout.decode().splitlines():
        if text_line.startswith("#"):
            continue
        fields = text_line.split("\t")
        if min_qual is not None and fields[VcfCol.QUAL] != "." and float(fields[VcfCol.QUAL]) < min_qual:
            continue
        if min_dp is not None:
            depth = next(
                (int(field[len(DP_PREFIX) :]) for field in fields[VcfCol.INFO].split(";") if field.startswith(DP_PREFIX)),
                0,
            )
            if depth < min_dp:
                continue
        for alt in fields[VcfCol.ALT].split(","):
            out.add((fields[VcfCol.CHROM], int(fields[VcfCol.POS]), fields[VcfCol.REF], alt))
    return out


def mapping_to_variants(work_dir, genome, threads, bcftools, samtools, picard=None, cov=DEFAULT_COVERAGE):
    """Section 3, which calls variants from optik's and minimap2's mappings and scores both against the truth."""
    print(
        "== 3. mapping -> variant calling (optik map vs minimap2 -ax sr; samtools + bcftools)"
    )
    optik, minimap2 = PATHS["optik"], PATHS["minimap2"]
    ref_recs = [(name, seq) for name, seq in read_fasta(genome)]
    ref = work_dir / "ref.fa"
    write_fasta(ref, ref_recs)
    run([samtools, "faidx", ref])
    sample, truth = plant_variants(ref_recs, seed=PLANT_SEED)
    sfa = work_dir / "sample.fa"
    write_fasta(sfa, sample)
    truth_vcf = work_dir / "truth.vcf"
    write_vcf(truth_vcf, truth, ref_recs)
    # Paired-end FASTQ, as germline sequencing delivers it (FR, 400 +- 50 bp fragments, with quality falling from Q37 to Q20 over the last 20 bases, which mpileup weighs).
    # Both mappers get R1 and R2 and must pass the qualities through.
    pre = work_dir / "sample_reads"
    glen = sum(len(seq) for _, seq in sample)
    run(
        [
            sys.executable,
            HERE / "gen_pairs.py",
            "-r",
            sfa,
            "-n",
            cov * glen // (MATES_PER_PAIR * PAIR_READ_LEN),
            "-l",
            PAIR_READ_LEN,
            "-e",
            0.005,
            "--ins",
            0.0005,
            "--del",
            0.0005,
            "-s",
            11,
            "-o",
            pre,
        ]
    )
    reads1, reads2 = f"{pre}_1.fq", f"{pre}_2.fq"
    truth_variants = vcf_set(truth_vcf, bcftools, ref)
    calls, stats = {}, {}
    # A read group, as germline pipelines attach one, since GATK and Picard refuse BAMs without it.
    # Both mappers take bwa's -R form.
    read_group = "@RG\\tID:s1\\tSM:sample1\\tPL:ILLUMINA"
    for name, cmd in (
        (
            "optik",
            [
                optik,
                "map",
                "-r",
                ref,
                "-q",
                reads1,
                "-q",
                reads2,
                "-t",
                threads,
                "-R",
                read_group,
                "-o",
                work_dir / "optik.sam",
            ],
        ),
        ("minimap2", None),
    ):
        sam = work_dir / f"{name}.sam"
        if name == "minimap2":
            run(
                [minimap2, "-ax", "sr", "-t", threads, "-R", read_group, ref, reads1, reads2],
                out=sam,
            )
        else:
            run(cmd)
        # The usual germline preprocessing of fixmate -m, coordinate sort and markdup
        fixmate_bam, srt, bam = (
            work_dir / f"{name}.fixmate.bam",
            work_dir / f"{name}.sorted.bam",
            work_dir / f"{name}.bam",
        )
        run([samtools, "fixmate", "-m", sam, fixmate_bam])
        run([samtools, "sort", "-@", threads, "-o", srt, fixmate_bam])
        run([samtools, "markdup", srt, bam])
        run([samtools, "index", bam])
        quickcheck_proc = subprocess.run(
            [samtools, "quickcheck", "-v", bam], capture_output=True, check=False
        )
        flagstat_json = run([samtools, "flagstat", "-O", "json", bam]).stdout.decode()
        fsj = json.loads(flagstat_json)["QC-passed reads"]
        depth_lines = run([samtools, "depth", "-a", bam]).stdout.decode().split("\n")
        # samtools depth -a prints "contig<TAB>position<TAB>depth" for every reference position.
        depths = [int(text_line.rsplit("\t", 1)[1]) for text_line in depth_lines if text_line]
        mean = sum(depths) / max(1, len(depths))
        frac10 = sum(depth >= MIN_BREADTH_DEPTH for depth in depths) / max(1, len(depths))
        stats[name] = {
            "primary_mapped": fsj.get("primary mapped"),
            "total": fsj.get("primary"),
            "properly_paired": fsj.get("properly paired"),
            "mean_depth": round(mean, 2),
            "frac_ge10x": round(frac10, 4),
        }
        check(
            f"{name}: BAM sorts, indexes and passes samtools quickcheck",
            quickcheck_proc.returncode == 0,
        )
        n_all = int(run([samtools, "view", "-c", bam]).stdout.decode())
        n_rg = int(run([samtools, "view", "-c", "-r", "s1", bam]).stdout.decode())
        check(
            f"{name}: read group s1 in the header and on every record",
            n_all > 0 and n_rg == n_all,
            records=n_all,
            with_rg=n_rg,
        )
        if picard:
            picard_proc = subprocess.run(
                [
                    str(picard),
                    "ValidateSamFile",
                    "-I",
                    str(bam),
                    "-R",
                    str(ref),
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
            errs = [line for line in txt.splitlines() if line.startswith("ERROR")]
            if name == "optik":
                check(
                    f"{name}: Picard ValidateSamFile (STRICT) finds no errors",
                    picard_proc.returncode == 0 and "No errors found" in txt,
                    errors=errs[:MAX_ERRORS_SHOWN],
                )
            else:
                # This is the reference tool's own output and is reported, not gated (minimap2 with fixmate leaves an MC tag on a read whose mate is unmapped).
                print(
                    f"  [info] {name}: Picard ValidateSamFile: "
                    f"{'no errors' if 'No errors found' in txt else '; '.join(errs[:INFO_ERRORS_SHOWN])}"
                )
        head = run([samtools, "view", bam]).stdout.decode().splitlines()[:QUAL_CHECK_RECORDS]
        check(
            f"{name}: SAM QUAL carries the FASTQ base qualities",
            bool(head)
            and all(
                text_line.split("\t")[SamCol.QUAL] != "*"
                for text_line in head
                if not int(text_line.split("\t")[SamCol.FLAG]) & SamFlag.UNMAPPED
            ),
        )
        vcf = work_dir / f"{name}.vcf.gz"
        mpileup_proc = subprocess.Popen(
            [str(bcftools), "mpileup", "-f", str(ref), "-a", "AD,DP", "-Ou", str(bam)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        call_proc = subprocess.run(
            [str(bcftools), "call", "-mv", "--ploidy", "1", "-Oz", "-o", str(vcf)],
            stdin=mpileup_proc.stdout,
            capture_output=True,
            check=False,
        )
        mpileup_proc.wait()
        if call_proc.returncode != 0:
            raise RuntimeError(call_proc.stderr.decode()[-DETAIL_CHARS:])
        calls[name] = {
            "raw": vcf_set(vcf, bcftools, ref),
            "filtered": vcf_set(vcf, bcftools, ref, min_qual=FILTER_MIN_QUAL, min_dp=FILTER_MIN_DP),
        }
        # Reads that contain a true indel must keep a usable MAPQ, since pipelines filter at MAPQ 20 to 30 and must not lose the indel evidence.
        indel_sites = [(chrom, pos) for chrom, pos, ref_allele, alt_allele in truth if len(ref_allele) != len(alt_allele)]
        indel_mapqs = []
        for chrom, pos in indel_sites[:MAX_INDEL_SITES]:
            out = run([samtools, "view", bam, f"{chrom}:{pos}-{pos + 1}"]).stdout.decode()
            for text_line in out.splitlines():
                fields = text_line.split("\t")
                if "I" in fields[SamCol.CIGAR] or "D" in fields[SamCol.CIGAR]:
                    indel_mapqs.append(int(fields[SamCol.MAPQ]))
        stats[name]["indel_reads"] = len(indel_mapqs)
        stats[name]["indel_reads_mapq_ge30"] = round(
            sum(mapq >= HIGH_MAPQ for mapq in indel_mapqs) / max(1, len(indel_mapqs)), 4
        )
    print(f"  alignment summary: {stats}")
    RESULTS["alignment_summary"] = stats
    # Both the raw and the filtered calls are scored, because a pipeline's filter can hide a mapper's errors or expose them.
    score = {}
    for name, name_calls in calls.items():
        score[name] = {}
        for view in ("raw", "filtered"):
            row = {}
            for kind, pred in (("snp", is_snp), ("indel", is_indel)):
                truth_set = {variant for variant in truth_variants if pred(variant)}
                call_set = {variant for variant in name_calls[view] if pred(variant)}
                true_pos = len(truth_set & call_set)
                row[kind] = {
                    "truth": len(truth_set),
                    "tp": true_pos,
                    "fp": len(call_set - truth_set),
                    "fn": len(truth_set - call_set),
                    "precision": round(true_pos / max(1, len(call_set)), 4),
                    "recall": round(true_pos / max(1, len(truth_set)), 4),
                }
            score[name][view] = row
    RESULTS["variant_calls"] = score
    print(f"  variant calls vs truth (raw, and QUAL>=20 & DP>=5): {json.dumps(score)}")
    optik_score, mm_score = score["optik"]["filtered"], score["minimap2"]["filtered"]
    for kind in ("snp", "indel"):
        check(
            f"{kind}s after the standard filter: optik recall and precision >= minimap2's - 0.5%",
            optik_score[kind]["recall"] >= mm_score[kind]["recall"] - CALL_TOL
            and optik_score[kind]["precision"] >= mm_score[kind]["precision"] - CALL_TOL,
            optik=optik_score[kind],
            minimap2=mm_score[kind],
        )
    # Matching totals is not enough, so the two mappers must also find nearly the same true variants.
    optik_true = calls["optik"]["filtered"] & truth_variants
    mm_true = calls["minimap2"]["filtered"] & truth_variants
    check(
        "the same true variants are found (optik vs minimap2, filtered; >= 99.5% shared)",
        len(optik_true & mm_true) >= MIN_SHARED_TRUE * max(len(optik_true), len(mm_true)),
        both=len(optik_true & mm_true),
        only_optik=len(optik_true - mm_true),
        only_minimap2=len(mm_true - optik_true),
    )
    check(
        "reads carrying a true indel keep MAPQ >= 30 (survive pipeline MAPQ filters) as often as minimap2's",
        stats["optik"]["indel_reads_mapq_ge30"]
        >= stats["minimap2"]["indel_reads_mapq_ge30"] - INDEL_MAPQ_TOL,
        optik=stats["optik"]["indel_reads_mapq_ge30"],
        minimap2=stats["minimap2"]["indel_reads_mapq_ge30"],
    )
    optik_stats, mm_stats = stats["optik"], stats["minimap2"]
    check(
        "coverage usable like minimap2's (mean depth within 2%, >=10x breadth within 0.5%)",
        abs(optik_stats["mean_depth"] - mm_stats["mean_depth"]) <= DEPTH_REL_TOL * mm_stats["mean_depth"]
        and abs(optik_stats["frac_ge10x"] - mm_stats["frac_ge10x"]) <= BREADTH_TOL,
        optik=optik_stats,
        minimap2=mm_stats,
    )
    check(
        "properly paired like minimap2's (samtools flagstat, within 0.5% of the reads)",
        optik_stats["properly_paired"] >= mm_stats["properly_paired"] - PROPER_PAIR_TOL * mm_stats["total"],
        optik=optik_stats["properly_paired"],
        minimap2=mm_stats["properly_paired"],
    )


def main():
    """Prepare the shared inputs and run the three sections."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--genome", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--bcftools", default=shutil.which("bcftools"))
    parser.add_argument("--samtools", default=shutil.which("samtools"))
    parser.add_argument(
        "--gfatools", default=str(HERE.parent / "temp" / "gfatools" / "gfatools")
    )
    parser.add_argument(
        "--picard",
        default=os.environ.get("PICARD_BIN") or shutil.which("picard"),
        help="Picard launcher; if found, every BAM is checked with ValidateSamFile",
    )
    parser.add_argument(
        "--skip", default="", help="comma list of sections to skip: kmer,graph,map"
    )
    args = parser.parse_args()
    work_dir = Path(args.workdir)
    work_dir.mkdir(parents=True, exist_ok=True)
    skip = set(filter(None, args.skip.split(",")))
    genome = work_dir / "genome.fa"
    write_fasta(genome, read_fasta(args.genome))
    kmer_len, threads = args.k, args.threads
    # The shared inputs are Cuttlefish unitigs, which SSHash, SBWT and FMSI pipelines index, and 20k simulated reads (1% errors) plus 2k foreign ones.
    cfp = work_dir / "cf"
    run(
        [
            PATHS["cuttlefish"],
            "build",
            "--ref",
            "-s",
            genome,
            "-k",
            kmer_len,
            "-t",
            threads,
            "-o",
            cfp,
            "-w",
            work_dir,
        ]
    )
    cf_fa = Path(str(cfp) + ".fa")
    reads = work_dir / "reads.fa"
    run(
        [
            sys.executable,
            HERE / "gen_reads.py",
            "-r",
            genome,
            "-n",
            20000,
            "-l",
            150,
            "-e",
            0.01,
            "-s",
            5,
            "-o",
            reads,
        ]
    )
    rnd = random.Random(FOREIGN_SEED)
    with open(reads, "a") as handle:
        handle.writelines(
            f">foreign{i}\n{''.join(rnd.choice('ACGT') for _ in range(FOREIGN_READ_LEN))}\n"
            for i in range(FOREIGN_READS)
        )
        handle.write(">withN\n" + "ACGT" * 20 + "N" + "TTGCA" * 20 + "\n")
    if "kmer" not in skip:
        kmer_parity(work_dir, genome, cf_fa, reads, kmer_len, threads)
    if "graph" not in skip:
        if not Path(args.gfatools).exists():
            check("gfatools available", False, path=args.gfatools)
        else:
            graph_usability(work_dir, genome, cf_fa, reads, kmer_len, threads, args.gfatools)
    if "map" not in skip:
        if not args.bcftools or not args.samtools:
            check(
                "bcftools and samtools available",
                False,
                bcftools=args.bcftools,
                samtools=args.samtools,
            )
        else:
            mapping_to_variants(work_dir, genome, threads, args.bcftools, args.samtools, args.picard)
    (work_dir / "pipeline_parity.json").write_text(json.dumps(RESULTS, indent=1, default=str))
    print(
        f"\npipeline_parity: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED: ' + ', '.join(FAILS)}"
    )
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
