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
import validate_unitigs as vu
from config import PATHS
from sam_spec import HIGH_MAPQ, SAM_FLAG_UNMAPPED

RESULTS = {}
FAILS = []


def check(name, ok, **detail):
    """Record and print one check result, remembering failures in FAILS."""
    RESULTS[name] = {"pass": bool(ok), **detail}
    print(
        f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""),
        flush=True,
    )
    if not ok:
        FAILS.append(name)


def run(cmd, out=None, **kw):
    cmd = [str(c) for c in cmd]
    if out is not None:
        with open(out, "wb") as fh:
            p = subprocess.run(
                cmd, stdout=fh, stderr=subprocess.PIPE, **kw, check=False
            )
    else:
        p = subprocess.run(cmd, capture_output=True, **kw, check=False)
    if p.returncode != 0:
        raise RuntimeError(
            f"failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr.decode()[-800:]}"
        )
    return p


def read_fasta(path):
    """FASTA records, plain or gzipped, as (name, sequence) pairs."""
    op = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] in ">@":
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def write_fasta(path, recs):
    """Write (name, sequence) records as single-line FASTA."""
    with open(path, "w") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in recs)


# Section 1, the k-mer dictionary
def kmer_parity(wd, genome, unitigs, reads, k, threads):
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
    for _, s in read_fasta(unitigs):
        for v in vu.windows(s, k):
            kset.add(vu.canon(v, k))

    def oracle_vec(s):
        """Presence of each k-mer window of s in the oracle set, None for windows with non-ACGT bases."""
        out = []
        for i in range(len(s) - k + 1):
            w = s[i : i + k].upper()
            if set(w) - set("ACGT"):
                out.append(None)  # invalid window
            else:
                out.append(vu.canon(vu.pack(w), k) in kset)
        return out

    truth = [oracle_vec(s) for _, s in recs]
    # optik --ids
    idx = wd / "o.spss"
    run([optik, "build-spss", "-i", unitigs, "-k", k, "-t", threads, "-o", idx])
    ids_path = wd / "o.ids"
    run([optik, "query", "-i", idx, "-q", reads, "--ids", ids_path, "-t", threads])
    optik_ids = [
        list(map(int, text_line.split()))
        for text_line in Path(ids_path).read_text().split("\n")[: len(recs)]
    ]
    # SBWT indexes both strands, as its documentation prescribes for canonical lookups.
    sbi = wd / "x.sbwt"
    run(
        [
            sbwt_bin,
            "build",
            "-i",
            unitigs,
            "-o",
            sbi,
            "-k",
            k,
            "--add-reverse-complements",
            "--temp-dir",
            wd,
            "-t",
            threads,
        ]
    )
    sbo = wd / "sbwt.out"
    run([sbwt_bin, "search", "-i", sbi, "-q", reads, "-o", sbo])
    sbwt = [
        list(map(int, text_line.split()))
        for text_line in Path(sbo).read_text().split("\n")[: len(recs)]
    ]
    # FMSI over KmerCamel's masked superstring (-S, since the input is unitigs)
    ms = wd / "ms.fa"
    run([kmercamel, "compute", "-k", k, "-S", "-M", ms, "-o", "/dev/null", unitigs])
    run([fmsi_bin, "index", ms])
    fmo = wd / "fmsi.out"
    run([fmsi_bin, "query", "-q", reads, ms], out=fmo)
    fmsi = [
        text_line.split("\t")[1] if "\t" in text_line else ""
        for text_line in Path(fmo).read_text().split("\n")[: len(recs)]
    ]

    n_win = n_valid = 0
    bad = Counter()
    hits = {"oracle": [], "optik": [], "sbwt": [], "fmsi": []}
    for r, t in enumerate(truth):
        ho = hs = hf = ht = 0
        for i, want in enumerate(t):
            n_win += 1
            if want is None:
                # An invalid window is -1 in optik and absent in SBWT and FMSI.
                if optik_ids[r][i] != -1:
                    bad["optik invalid window"] += 1
                continue
            n_valid += 1
            o = optik_ids[r][i] >= 0
            s_ = sbwt[r][i] >= 0
            f_ = fmsi[r][i] == "1" if i < len(fmsi[r]) else None
            if o != want:
                bad["optik"] += 1
            if s_ != want:
                bad["sbwt"] += 1
            if f_ != want:
                bad["fmsi"] += 1
            ht += want
            ho += o
            hs += s_
            hf += bool(f_)
        hits["oracle"].append(ht)
        hits["optik"].append(ho)
        hits["sbwt"].append(hs)
        hits["fmsi"].append(hf)
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
    same = sum(
        a == b == c
        for a, b, c in zip(hits["optik"], hits["sbwt"], hits["fmsi"], strict=False)
    )
    check(
        "per-read k-mer hit counts identical (optik, SBWT, FMSI) -- what read screening uses",
        same == len(recs),
        reads=len(recs),
        identical=same,
    )
    # Dictionary semantics, one id per canonical k-mer, independent of strand
    by_kmer, ok = {}, True
    for (_n, s), row in zip(recs, optik_ids, strict=False):
        for i, v in enumerate(row):
            if v < 0:
                continue
            c = vu.canon(vu.pack(s[i : i + k].upper()), k)
            ok &= by_kmer.setdefault(c, v) == v
    ids = set(by_kmer.values())
    check(
        "optik --ids is a dictionary: one id per canonical k-mer (both strands), distinct k-mers distinct ids",
        ok and len(ids) == len(by_kmer),
        kmers=len(by_kmer),
    )


# Section 2, the compacted dBG as a graph
def parse_gfa(path):
    """Segments {name: sequence} and links {(from, orient, to, orient)} of a GFA file."""
    op = gzip.open if str(path).endswith(".gz") else open
    segs, links = {}, set()
    with op(path, "rt") as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t[0] == "S":
                segs[t[1]] = t[2].upper()
            elif t[0] == "L":
                links.add((t[1], t[2], t[3], t[4]))
    return segs, links


def genome_walk(segs, links, genome, k):
    """Walk each N-free run of the genome k-mer by k-mer through the graph.

    Every k-mer must lie in a segment, and consecutive k-mers must either be consecutive in one segment or cross a segment boundary along an L line.
    Returns (k-mers walked, missing k-mers, steps without an edge)."""
    loc = {}  # canonical k-mer -> segment, offset, forward flag
    for sid, s in segs.items():
        for i, v in enumerate(vu.windows(s, k)):
            r = vu.rc(v, k)
            loc[min(v, r)] = (sid, i, v <= r)
    flip = {"+": "-", "-": "+"}
    lk = set(links) | {(b, flip[ob], a, flip[oa]) for a, oa, b, ob in links}
    walked = missing = noedge = 0
    for _, g in read_fasta(genome):
        prev = None
        for run_ in g.upper().replace("N", " ").split():
            prev = None
            for v in vu.windows(run_, k):
                walked += 1
                r = vu.rc(v, k)
                c = min(v, r)
                if c not in loc:
                    missing += 1
                    prev = None
                    continue
                sid, off, fwd = loc[c]
                # Orientation of the segment as this k-mer reads it
                o = "+" if (v == c) == fwd else "-"
                n_k = len(segs[sid]) - k + 1
                pos = (
                    off if o == "+" else n_k - 1 - off
                )  # position along the walking direction
                cur = (sid, o, pos, n_k)
                if prev is not None:
                    ps, po, pp, pn = prev
                    inside = ps == sid and po == o and pos == pp + 1
                    across = pp == pn - 1 and pos == 0 and (ps, po, sid, o) in lk
                    if not (inside or across):
                        noedge += 1
                prev = cur
    return walked, missing, noedge


def graph_usability(wd, genome, cf_unitigs, reads, k, threads, gfatools):
    """Section 2, which checks optik's GFA as graph tools read it and its unitigs as index input."""
    print(
        "== 2. compacted dBG as a graph (optik GFA vs Bifrost GFA; unitigs vs Cuttlefish)"
    )
    optik, bifrost, sshash = PATHS["optik"], PATHS["bifrost"], PATHS["sshash"]
    og = wd / "optik.gfa"
    run([optik, "dbruijn", "-i", genome, "-k", k, "-c", "1", "-t", threads, "-o", og])
    bfp = wd / "bf"
    run([bifrost, "build", "-r", genome, "-k", k, "-t", threads, "-o", bfp])
    bg = next(
        p for p in (Path(str(bfp) + ".gfa.gz"), Path(str(bfp) + ".gfa")) if p.exists()
    )
    for name, g in (("optik", og), ("bifrost", bg)):
        src = g
        if str(g).endswith(".gz"):
            src = wd / f"{name}_plain.gfa"
            with gzip.open(g, "rb") as fi, open(src, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        p = subprocess.run(
            [gfatools, "stat", src], capture_output=True, text=True, check=False
        )
        stat = dict(
            text_line.split(": ", 1)
            for text_line in p.stdout.splitlines()
            if ": " in text_line
        )
        segs, links = parse_gfa(g)
        # GFA writers follow different conventions.
        # optik writes each bidirected edge once and Bifrost writes both orientations.
        # gfatools merges the two orientations and counts links as distinct edges minus the self-mirrored ones that join a segment to its own reverse, and arcs as twice the edges minus the self-mirrored ones.
        flip = {"+": "-", "-": "+"}
        edges = {min(e, (e[2], flip[e[3]], e[0], flip[e[1]])) for e in links}
        selfm = sum(1 for e in edges if e == (e[2], flip[e[3]], e[0], flip[e[1]]))
        check(
            f"gfatools parses the {name} GFA (segments, edges, arcs and length as written)",
            p.returncode == 0
            and int(stat.get("Number of segments", -1)) == len(segs)
            and int(stat.get("Number of links", -1)) == len(edges) - selfm
            and int(stat.get("Number of arcs", -1)) == 2 * len(edges) - selfm
            and int(stat.get("Total segment length", -1))
            == sum(map(len, segs.values())),
            segments=stat.get("Number of segments"),
            links=stat.get("Number of links"),
            l_lines=len(links),
            edges=len(edges),
            total_len=stat.get("Total segment length"),
        )
        walked, missing, noedge = genome_walk(segs, links, genome, k)
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

    def rcs(x):
        """Reverse complement of x."""
        return x.translate(comp)[::-1]

    def canon(x):
        """Canonical form of k-mer string x."""
        return min(x, rcs(x))

    oseg, olinks = parse_gfa(og)
    bsrc = bg if not str(bg).endswith(".gz") else wd / "bifrost_plain.gfa"
    bseg, blinks = parse_gfa(bsrc)
    cfs = {canon(q) for _, q in read_fasta(cf_unitigs)}
    os_ = {canon(q) for q in oseg.values()}
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
        for a, oa, b, ob in links:
            sa = segs[a] if oa == "+" else rcs(segs[a])
            sb = segs[b] if ob == "+" else rcs(segs[b])
            out.add(canon(sa[-k:] + sb[k - 1]))
        return out

    gk1 = set()  # the genome's (k+1)-mers, i.e. its edges
    for _, g in read_fasta(genome):
        g = g.upper()
        for i in range(len(g) - k):
            x = g[i : i + k + 1]
            if "N" not in x:
                gk1.add(canon(x))
    internal = set()
    for q in oseg.values():
        for i in range(len(q) - k):
            internal.add(canon(q[i : i + k + 1]))
    oe, be = edge_set(oseg, olinks), edge_set(bseg, blinks)
    check(
        "links: every L line is a genome (k+1)-mer, and links + unitig-internal edges = all of them",
        oe <= gk1 and (oe | internal) == gk1 and not (oe & internal),
        links=len(oe),
        not_in_genome=len(oe - gk1),
        genome_edges_missing=len(gk1 - oe - internal),
    )
    # Bifrost is node-centric by design.
    # It links two k-mers that overlap by k-1 even when that (k+1)-mer never occurs, and such branches split unitigs, so the edges at the splits become links as well.
    extra = be - oe
    check(
        "Bifrost's extra links: unobserved k-1 overlaps, or edges inside optik's unitigs that its branches split",
        oe <= be and extra <= ((extra - gk1) | internal),
        bifrost_only=len(extra),
        unobserved=len(extra - gk1),
        split_internal=len(extra & internal),
        unexplained=len(extra - (extra - gk1) - internal),
        optik_only=len(oe - be),
    )
    del gk1, internal
    cnt = {}
    for _, g in read_fasta(genome):
        g = g.upper()
        for i in range(len(g) - k + 1):
            x = g[i : i + k]
            if "N" in x:
                continue
            c = canon(x)
            cnt[c] = cnt.get(c, 0) + 1
    bad_kc = bad_km = 0
    with open(og) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t[0] != "S":
                continue
            q = t[2]
            want = sum(cnt.get(canon(q[i : i + k]), 0) for i in range(len(q) - k + 1))
            tags = dict(x.split(":", 1) for x in t[3:])
            kc = int(tags["KC"].split(":")[1])
            bad_kc += kc != want
            bad_km += tags["km"].split(":")[1] != f"{want / (len(q) - k + 1):.2f}"
    check(
        "KC = the summed genome abundance of each segment's k-mers, km = KC / k-mers",
        bad_kc == 0 and bad_km == 0,
        wrong_kc=bad_kc,
        wrong_km=bad_km,
    )
    del cnt
    # A run of exactly k bases, between N runs or as a whole record, holds one k-mer.
    # Cuttlefish 3.0.3 and 3.1.0 drop such k-mers, and optik must keep them.
    rr = random.Random(12)

    def rnd(n):
        """Random ACGT string of length n."""
        return "".join(rr.choice("ACGT") for _ in range(n))

    nr = wd / "exact_k_runs.fa"
    write_fasta(
        nr,
        [
            (
                "x",
                rnd(300)
                + "NNN"
                + rnd(k)
                + "NN"
                + rnd(300)
                + "N"
                + rnd(k + 1)
                + "NN"
                + rnd(200),
            ),
            ("y", rnd(k)),
        ],
    )
    truth_k = set()
    for _, q in read_fasta(nr):
        for runs in q.split("N"):
            for i in range(len(runs) - k + 1):
                truth_k.add(canon(runs[i : i + k]))

    def kset(path):
        """Canonical k-mers of every record in a FASTA file."""
        out = set()
        for _, q in read_fasta(path):
            for i in range(len(q) - k + 1):
                out.add(canon(q[i : i + k]))
        return out

    for mode in ("part", "table"):
        nro = wd / f"exact_k_{mode}.fa"
        run(
            [optik, "dbruijn", "-i", nr, "-k", k, "-c", "1", "-o", nro, "--graph", mode]
        )
        check(
            f"k-mers of runs exactly k long are kept (--graph {mode})",
            kset(nro) == truth_k,
            truth=len(truth_k),
            got=len(kset(nro)),
        )
    ncf = wd / "exact_k_cf"
    run(
        [
            PATHS["cuttlefish"],
            "build",
            "--ref",
            "-s",
            nr,
            "-k",
            k,
            "-t",
            2,
            "-o",
            ncf,
            "-w",
            wd,
        ]
    )
    print(
        f"  [info] Cuttlefish on the same input: {len(truth_k - kset(str(ncf) + '.fa'))} of "
        f"{len(truth_k)} k-mers missing (runs of exactly k bases)"
    )
    # Partitioned compaction, the default here, and the table walker must agree.
    tg = wd / "optik_table.gfa"
    run(
        [
            optik,
            "dbruijn",
            "-i",
            genome,
            "-k",
            k,
            "-c",
            "1",
            "-t",
            threads,
            "-o",
            tg,
            "--graph",
            "table",
        ]
    )
    check(
        "dbruijn --graph table writes the same GFA as the default (partitioned) build",
        Path(tg).read_bytes() == Path(og).read_bytes(),
    )
    # optik's unitigs serve as the input of a downstream index in place of Cuttlefish's, taken from dbruijn's own FASTA output (-o *.fa) as a pipeline would use it.
    # These are the same sequences as the GFA's segments.
    segs, _ = parse_gfa(og)
    ou = wd / "optik_unitigs.fa"
    run([optik, "dbruijn", "-i", genome, "-k", k, "-c", "1", "-t", threads, "-o", ou])
    lines = Path(ou).read_text().splitlines()
    heads, seqs = lines[0::2], lines[1::2]
    gfa_tags = {}
    with open(og) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t[0] == "S":
                gfa_tags[t[1]] = t[3:]
    check(
        "dbruijn -o *.fa: unitig FASTA with the GFA's segments (ids, sequences, LN/KC/km)",
        len(heads) == len(segs)
        and all(h.startswith(">") for h in heads)
        and [(h[1:].split()[0], q) for h, q in zip(heads, seqs, strict=False)]
        == list(segs.items())
        and all(h[1:].split()[1:] == gfa_tags[h[1:].split()[0]] for h in heads),
        records=len(heads),
        segments=len(segs),
    )
    res = {}
    for name, u in (("optik", ou), ("cuttlefish", cf_unitigs)):
        si = wd / f"{name}.sshash"
        run(
            [
                sshash,
                "build",
                "-i",
                u,
                "-k",
                k,
                "-m",
                "17" if k >= 25 else str(max(5, k // 2 + 2)),
                "-t",
                threads,
                "-o",
                si,
            ]
        )
        p = run([sshash, "query", "-i", si, "-q", reads])
        txt = p.stdout.decode() + p.stderr.decode()
        res[name] = {
            kk: v
            for kk, v in (
                text_line.split(" = ")
                for text_line in txt.splitlines()
                if " = " in text_line
            )
            if kk in ("num_kmers", "num_positive_kmers", "num_negative_kmers")
        }
    check(
        "SSHash built on optik's unitigs answers the reads exactly as on Cuttlefish's",
        res["optik"] == res["cuttlefish"] and bool(res["optik"]),
        optik=res["optik"],
        cuttlefish=res["cuttlefish"],
    )


# Section 3, mapping followed by variant calling
def plant_variants(ref_recs, seed, snp_every=2000, indel_every=15000):
    """Make a sample genome from the reference with SNPs and 1-10 bp indels at spaced, N-free positions.

    Returns (sample records, truth rows of (chrom, pos1, ref, alt))."""
    rnd = random.Random(seed)
    sample, truth = [], []
    for name, s in ref_recs:
        out, last, i = [], 0, 200
        while i < len(s) - 200:
            step = rnd.randint(snp_every // 2, snp_every * 3 // 2)
            i += step
            if i >= len(s) - 200:
                break
            win = s[i - 60 : i + 60]
            if "N" in win.upper():
                continue
            out.append(s[last:i])
            r = s[i].upper()
            if rnd.random() < snp_every / indel_every:  # an indel
                ln = rnd.randint(1, 10)
                if rnd.random() < 0.5:  # an insertion after base i
                    ins = "".join(rnd.choice("ACGT") for _ in range(ln))
                    out.append(s[i] + ins)
                    truth.append((name, i + 1, r, r + ins))
                    last = i + 1
                else:  # a deletion of ln bases after base i
                    out.append(s[i])
                    truth.append((name, i + 1, s[i : i + 1 + ln].upper(), r))
                    last = i + 1 + ln
            else:
                alt = rnd.choice([b for b in "ACGT" if b != r])
                out.append(alt)
                truth.append((name, i + 1, r, alt))
                last = i + 1
            i = last + 60
        out.append(s[last:])
        sample.append((name, "".join(out)))
    return sample, truth


def write_vcf(path, truth, ref_recs):
    """Write the truth variants as a minimal VCF 4.2 with contig lines for ref_recs."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.writelines(f"##contig=<ID={n},length={len(s)}>\n" for n, s in ref_recs)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.writelines(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n" for c, p, r, a in truth)


def vcf_set(path, bcftools, ref, min_qual=None, min_dp=None):
    """(chrom, pos, ref, alt) of every ALT allele after left-normalization.

    Optionally only records with QUAL >= min_qual and INFO/DP >= min_dp are kept, the hard filter that variant-calling pipelines apply."""
    p = run([bcftools, "norm", "-f", ref, "-m", "-any", "-Ov", path])
    out = set()
    for text_line in p.stdout.decode().splitlines():
        if text_line.startswith("#"):
            continue
        t = text_line.split("\t")
        if min_qual is not None and t[5] != "." and float(t[5]) < min_qual:
            continue
        if min_dp is not None:
            dp = next((int(x[3:]) for x in t[7].split(";") if x.startswith("DP=")), 0)
            if dp < min_dp:
                continue
        for a in t[4].split(","):
            out.add((t[0], int(t[1]), t[3], a))
    return out


def mapping_to_variants(wd, genome, threads, bcftools, samtools, picard=None, cov=30):
    """Section 3, which calls variants from optik's and minimap2's mappings and scores both against the truth."""
    print(
        "== 3. mapping -> variant calling (optik map vs minimap2 -ax sr; samtools + bcftools)"
    )
    optik, minimap2 = PATHS["optik"], PATHS["minimap2"]
    ref_recs = [(n, s) for n, s in read_fasta(genome)]
    ref = wd / "ref.fa"
    write_fasta(ref, ref_recs)
    run([samtools, "faidx", ref])
    sample, truth = plant_variants(ref_recs, seed=7)
    sfa = wd / "sample.fa"
    write_fasta(sfa, sample)
    tv = wd / "truth.vcf"
    write_vcf(tv, truth, ref_recs)
    # Paired-end FASTQ, as germline sequencing delivers it (FR, 400 +- 50 bp fragments, with quality falling from Q37 to Q20 over the last 20 bases, which mpileup weighs).
    # Both mappers get R1 and R2 and must carry the qualities through.
    pre = wd / "sample_reads"
    glen = sum(len(x) for _, x in sample)
    run(
        [
            sys.executable,
            HERE / "gen_pairs.py",
            "-r",
            sfa,
            "-n",
            cov * glen // 300,
            "-l",
            150,
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
    r1, r2 = f"{pre}_1.fq", f"{pre}_2.fq"
    truth_variants = vcf_set(tv, bcftools, ref)
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
                r1,
                "-q",
                r2,
                "-t",
                threads,
                "-R",
                read_group,
                "-o",
                wd / "optik.sam",
            ],
        ),
        ("minimap2", None),
    ):
        sam = wd / f"{name}.sam"
        if name == "minimap2":
            run(
                [minimap2, "-ax", "sr", "-t", threads, "-R", read_group, ref, r1, r2],
                out=sam,
            )
        else:
            run(cmd)
        # The usual germline preprocessing of fixmate -m, coordinate sort and markdup
        fx, srt, bam = (
            wd / f"{name}.fixmate.bam",
            wd / f"{name}.sorted.bam",
            wd / f"{name}.bam",
        )
        run([samtools, "fixmate", "-m", sam, fx])
        run([samtools, "sort", "-@", threads, "-o", srt, fx])
        run([samtools, "markdup", srt, bam])
        run([samtools, "index", bam])
        qc = subprocess.run(
            [samtools, "quickcheck", "-v", bam], capture_output=True, check=False
        )
        fs = run([samtools, "flagstat", "-O", "json", bam]).stdout.decode()
        fsj = json.loads(fs)["QC-passed reads"]
        dp = run([samtools, "depth", "-a", bam]).stdout.decode().split("\n")
        depths = [int(text_line.split("\t")[2]) for text_line in dp if text_line]
        mean = sum(depths) / max(1, len(depths))
        frac10 = sum(d >= 10 for d in depths) / max(1, len(depths))
        stats[name] = {
            "primary_mapped": fsj.get("primary mapped"),
            "total": fsj.get("primary"),
            "properly_paired": fsj.get("properly paired"),
            "mean_depth": round(mean, 2),
            "frac_ge10x": round(frac10, 4),
        }
        check(
            f"{name}: BAM sorts, indexes and passes samtools quickcheck",
            qc.returncode == 0,
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
            pv = subprocess.run(
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
            txt = pv.stdout + pv.stderr
            errs = [x for x in txt.splitlines() if x.startswith("ERROR")]
            if name == "optik":
                check(
                    f"{name}: Picard ValidateSamFile (STRICT) finds no errors",
                    pv.returncode == 0 and "No errors found" in txt,
                    errors=errs[:5],
                )
            else:
                # This is the reference tool's own output and is reported, not gated (minimap2 with fixmate leaves an MC tag on a read whose mate is unmapped).
                print(
                    f"  [info] {name}: Picard ValidateSamFile: "
                    f"{'no errors' if 'No errors found' in txt else '; '.join(errs[:3])}"
                )
        head = run([samtools, "view", bam]).stdout.decode().splitlines()[:200]
        check(
            f"{name}: SAM QUAL carries the FASTQ base qualities",
            bool(head)
            and all(
                text_line.split("\t")[10] != "*"
                for text_line in head
                if not int(text_line.split("\t")[1]) & SAM_FLAG_UNMAPPED
            ),
        )
        vcf = wd / f"{name}.vcf.gz"
        mp = subprocess.Popen(
            [str(bcftools), "mpileup", "-f", str(ref), "-a", "AD,DP", "-Ou", str(bam)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        cp = subprocess.run(
            [str(bcftools), "call", "-mv", "--ploidy", "1", "-Oz", "-o", str(vcf)],
            stdin=mp.stdout,
            capture_output=True,
            check=False,
        )
        mp.wait()
        if cp.returncode != 0:
            raise RuntimeError(cp.stderr.decode()[-500:])
        calls[name] = {
            "raw": vcf_set(vcf, bcftools, ref),
            "filtered": vcf_set(vcf, bcftools, ref, min_qual=20, min_dp=5),
        }
        # Reads carrying a true indel must keep a usable MAPQ, since pipelines filter at MAPQ 20 to 30 and must not lose the indel evidence.
        indel_sites = [(c, p) for c, p, r, a in truth if len(r) != len(a)]
        mq = []
        for c, p in indel_sites[:200]:
            out = run([samtools, "view", bam, f"{c}:{p}-{p + 1}"]).stdout.decode()
            for text_line in out.splitlines():
                f = text_line.split("\t")
                if "I" in f[5] or "D" in f[5]:
                    mq.append(int(f[4]))
        stats[name]["indel_reads"] = len(mq)
        stats[name]["indel_reads_mapq_ge30"] = round(
            sum(m >= HIGH_MAPQ for m in mq) / max(1, len(mq)), 4
        )
    print(f"  alignment summary: {stats}")
    RESULTS["alignment_summary"] = stats
    score = {}
    for name, cs in calls.items():
        score[name] = {}
        for view in ("raw", "filtered"):
            row = {}
            for kind, pred in (
                ("snp", lambda v: len(v[2]) == len(v[3]) == 1),
                ("indel", lambda v: len(v[2]) != len(v[3])),
            ):
                t_ = {v for v in truth_variants if pred(v)}
                c_ = {v for v in cs[view] if pred(v)}
                tp = len(t_ & c_)
                row[kind] = {
                    "truth": len(t_),
                    "tp": tp,
                    "fp": len(c_ - t_),
                    "fn": len(t_ - c_),
                    "precision": round(tp / max(1, len(c_)), 4),
                    "recall": round(tp / max(1, len(t_)), 4),
                }
            score[name][view] = row
    RESULTS["variant_calls"] = score
    print(f"  variant calls vs truth (raw, and QUAL>=20 & DP>=5): {json.dumps(score)}")
    o, m = score["optik"]["filtered"], score["minimap2"]["filtered"]
    for kind in ("snp", "indel"):
        check(
            f"{kind}s after the standard filter: optik recall and precision >= minimap2's - 0.5%",
            o[kind]["recall"] >= m[kind]["recall"] - 0.005
            and o[kind]["precision"] >= m[kind]["precision"] - 0.005,
            optik=o[kind],
            minimap2=m[kind],
        )
    to = calls["optik"]["filtered"] & truth_variants
    tm = calls["minimap2"]["filtered"] & truth_variants
    check(
        "the same true variants are found (optik vs minimap2, filtered; >= 99.5% shared)",
        len(to & tm) >= 0.995 * max(len(to), len(tm)),
        both=len(to & tm),
        only_optik=len(to - tm),
        only_minimap2=len(tm - to),
    )
    check(
        "reads carrying a true indel keep MAPQ >= 30 (survive pipeline MAPQ filters) as often as minimap2's",
        stats["optik"]["indel_reads_mapq_ge30"]
        >= stats["minimap2"]["indel_reads_mapq_ge30"] - 0.01,
        optik=stats["optik"]["indel_reads_mapq_ge30"],
        minimap2=stats["minimap2"]["indel_reads_mapq_ge30"],
    )
    so, sm = stats["optik"], stats["minimap2"]
    check(
        "coverage usable like minimap2's (mean depth within 2%, >=10x breadth within 0.5%)",
        abs(so["mean_depth"] - sm["mean_depth"]) <= 0.02 * sm["mean_depth"]
        and abs(so["frac_ge10x"] - sm["frac_ge10x"]) <= 0.005,
        optik=so,
        minimap2=sm,
    )
    check(
        "properly paired like minimap2's (samtools flagstat, within 0.5% of the reads)",
        so["properly_paired"] >= sm["properly_paired"] - 0.005 * sm["total"],
        optik=so["properly_paired"],
        minimap2=sm["properly_paired"],
    )


def main():
    """Prepare the shared inputs and run the three sections."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--genome", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bcftools", default=shutil.which("bcftools"))
    ap.add_argument("--samtools", default=shutil.which("samtools"))
    ap.add_argument(
        "--gfatools", default=str(HERE.parent / "temp" / "gfatools" / "gfatools")
    )
    ap.add_argument(
        "--picard",
        default=os.environ.get("PICARD_BIN") or shutil.which("picard"),
        help="Picard launcher; if found, every BAM is checked with ValidateSamFile",
    )
    ap.add_argument(
        "--skip", default="", help="comma list of sections to skip: kmer,graph,map"
    )
    a = ap.parse_args()
    wd = Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    skip = set(filter(None, a.skip.split(",")))
    genome = wd / "genome.fa"
    write_fasta(genome, read_fasta(a.genome))
    k, threads = a.k, a.threads
    # The shared inputs are Cuttlefish unitigs, which SSHash, SBWT and FMSI pipelines index, and 20k simulated reads (1% errors) plus 2k foreign ones.
    cfp = wd / "cf"
    run(
        [
            PATHS["cuttlefish"],
            "build",
            "--ref",
            "-s",
            genome,
            "-k",
            k,
            "-t",
            threads,
            "-o",
            cfp,
            "-w",
            wd,
        ]
    )
    cf = Path(str(cfp) + ".fa")
    reads = wd / "reads.fa"
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
    rnd = random.Random(9)
    with open(reads, "a") as fh:
        fh.writelines(
            f">foreign{i}\n{''.join(rnd.choice('ACGT') for _ in range(150))}\n"
            for i in range(2000)
        )
        fh.write(">withN\n" + "ACGT" * 20 + "N" + "TTGCA" * 20 + "\n")
    if "kmer" not in skip:
        kmer_parity(wd, genome, cf, reads, k, threads)
    if "graph" not in skip:
        if not Path(a.gfatools).exists():
            check("gfatools available", False, path=a.gfatools)
        else:
            graph_usability(wd, genome, cf, reads, k, threads, a.gfatools)
    if "map" not in skip:
        if not a.bcftools or not a.samtools:
            check(
                "bcftools and samtools available",
                False,
                bcftools=a.bcftools,
                samtools=a.samtools,
            )
        else:
            mapping_to_variants(wd, genome, threads, a.bcftools, a.samtools, a.picard)
    (wd / "pipeline_parity.json").write_text(json.dumps(RESULTS, indent=1, default=str))
    print(
        f"\npipeline_parity: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED: ' + ', '.join(FAILS)}"
    )
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
