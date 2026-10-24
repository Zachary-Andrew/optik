#!/usr/bin/env python3
"""validate_unitigs.py - structural validation of a compacted dBG in GFA.

An independent oracle in pure Python that uses no optik code.
From the input FASTA it derives the solid canonical k-mer set and the edge-centric bidirected edge set, in which an edge exists iff its (k+1)-mer window occurs in either orientation and both flanking k-mers are solid.
It then checks the GFA.

  kmer_set_equal     the union of segment k-mers equals the solid set
  no_duplicates      every k-mer appears at exactly one segment position (comparing sets alone hides duplicated k-mers)
  paths_valid        consecutive k-mers of a segment are joined by an edge
  interior_compact   every internal edge E->X has outdeg(E) = 1 and indeg(X) = 1
  maximal            no segment end can be compacted further, which catches over-splitting that set equality cannot
  links_exact        the L lines equal the edges between segment ends (normalized, overlap k-1)
  segments_wellformed  segment names are unique and every sequence has at
                     least k ACGT bases (a repeated name would let a
                     dict-based reader drop a segment silently)
  tags_exact         LN:i is the length, KC:i the summed abundance of the segment's k-mers over true input counts and km:f is KC divided by the number of k-mers
  kdbg_exact         (with --kdbg) the -s node table has the right k, sorted canonical keys equal to the solid set, true abundances and in/out edge masks equal to the edge-centric oracle
  cuttlefish_partition (optional) the unitig partition compared exactly with
                     Cuttlefish (--read, edge-centric), where every
                     difference must be settled against this oracle

Usage
  python3 bench/validate_unitigs.py -i in.fa -g out.gfa -k 31 [-c 1]
         [--cuttlefish temp/cuttlefish/target/release/cuttlefish]

A non-zero exit status reports at least one failed check.
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_fasta(path):
    """Yield FASTA or FASTQ records, plain or .gz, with the reader validate_dbg shares."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from validate_dbg import read_fasta as _read

    return _read(path)


def windows(seq, w):
    """Yield the packed forward w-mers of seq, skipping any that contain a non-ACGT base."""
    mask = (1 << (2 * w)) - 1
    v = 0
    valid = 0
    for ch in seq.upper():
        b = ENC.get(ch)
        if b is None:
            valid = 0
            v = 0
            continue
        v = ((v << 2) | b) & mask
        valid += 1
        if valid >= w:
            yield v


# _RC_BYTE[b] is byte b (four 2-bit bases) complemented, with its bases in reverse order.
# rc() maps the 8 bytes of a packed w-mer (w <= 32) through it and reads them back in reverse byte order, which yields the whole reverse complement at C speed.
# A loop per base made the oracle take minutes per genome.
_RC_BYTE = bytes(
    sum(((3 - ((b >> (2 * i)) & 3)) << (2 * (3 - i))) for i in range(4))
    for b in range(256)
)


def rc(v, w):
    """Reverse complement of the packed w-mer v (2 bits per base, w <= 32).

    The 32 - w unused high bases complement into the low end and are shifted out."""
    return int.from_bytes(v.to_bytes(8, "little").translate(_RC_BYTE), "big") >> (
        2 * (32 - w)
    )


def _rc_slow(v, w):
    """Base-by-base reverse complement, the reference test_harness checks rc() against."""
    r = 0
    for _ in range(w):
        r = (r << 2) | (3 - (v & 3))
        v >>= 2
    return r


def canon(v, w):
    """Canonical form of the packed w-mer v, the smaller of v and its reverse complement."""
    r = rc(v, w)
    return min(r, v)


_DIGITS = str.maketrans("ACGT", "0123")


def pack(s):
    """Packed value of an ACGT string (A=0 to T=3), i.e. its base-4 number.

    Callers pass only ACGT strings, since segments are checked before packing."""
    return int(s.translate(_DIGITS), 4) if s else 0


def _pack_slow(s):
    """Base-by-base pack(), the reference test_harness checks pack() against."""
    v = 0
    for ch in s:
        v = (v << 2) | ENC[ch]
    return v


def parse_gfa(path):
    """Parse a GFA into segments {name: SEQ}, links, S-line tags {name: {tag: value}} and the names that occur on more than one S line."""
    segs, links, tags, dup_names = {}, [], {}, []
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t[0] == "S":
                if t[1] in segs:
                    dup_names.append(t[1])
                segs[t[1]] = t[2].upper()
                tags[t[1]] = {f[:2]: f[5:] for f in t[3:] if len(f) > 5 and f[2] == ":"}
            elif t[0] == "L":
                links.append((t[1], t[2], t[3], t[4], t[5]))
    return segs, links, tags, dup_names


def oracle(input_path, k, min_count=1):
    """Independent edge-centric oracle built from the raw input.

    Returns (solid, kp1, linked, counts).
    `solid` holds the canonical k-mers with count >= min_count and `kp1` the observed forward (k+1)-mer windows.
    `linked` holds the canonical k-mers that occur at least once inside a run longer than k bases, while the rest occur only as isolated ACGT runs of length k between non-ACGT characters.
    `counts` gives the occurrence count of every canonical k-mer.
    """
    counts = Counter()
    kp1 = set()
    linked = set()
    for _, seq in read_fasta(input_path):
        for v in windows(seq, k):
            counts[canon(v, k)] += 1
        for m in windows(seq, k + 1):
            kp1.add(m)
            linked.add(canon(m >> 2, k))
            linked.add(canon(m & ((1 << (2 * k)) - 1), k))
    solid = {c for c, n in counts.items() if n >= min_count}
    return solid, kp1, linked, counts


def validate(
    input_path,
    gfa_path,
    k,
    min_count=1,
    cuttlefish=None,
    threads=4,
    orc=None,
    kdbg=None,
):
    """Run every structural check and return {check_name: {"pass": bool, ...}}."""
    args = argparse.Namespace(
        input=input_path,
        gfa=gfa_path,
        k=k,
        min_count=min_count,
        cuttlefish=cuttlefish,
        threads=threads,
    )
    kmask = (1 << (2 * k)) - 1

    solid, kp1, linked, counts = (
        orc if orc is not None else oracle(args.input, k, args.min_count)
    )

    def has_edge(cur, nxt):
        """True when oriented edge cur -> nxt, with nxt = cur[1:] + b, is in the edge-centric graph."""
        m = (cur << 2) | (nxt & 3)
        return (m in kp1 or rc(m, k + 1) in kp1) and canon(nxt, k) in solid

    def succs(cur):
        """Successor k-mers of cur along edges of the graph."""
        out = []
        for b in range(4):
            nxt = ((cur << 2) | b) & kmask
            if has_edge(cur, nxt):
                out.append(nxt)
        return out

    def indeg(nxt):
        """In-degree of the k-mer nxt, counted as the successors of its reverse complement."""
        return len(succs(rc(nxt, k)))

    def compactable(cur):
        """The k-mer cur compacts into (unique successor with in-degree 1), or None."""
        s = succs(cur)
        if len(s) != 1:
            return None
        nxt = s[0]
        if canon(nxt, k) == canon(cur, k):
            return None
        if cur == rc(cur, k) or nxt == rc(
            nxt, k
        ):  # palindromes, possible at even k, stay single
            return None
        return nxt if indeg(nxt) == 1 else None

    segs, links, tags, dup_names = parse_gfa(args.gfa)
    checks = {}

    # Segment well-formedness, checked before anything indexes by name or k-mer
    bad_seq = [sid for sid, s in segs.items() if len(s) < k or set(s) - set("ACGT")]
    checks["segments_wellformed"] = {
        "pass": not dup_names and not bad_seq and bool(segs) == bool(solid),
        "duplicate_names": dup_names[:5],
        "short_or_non_acgt": bad_seq[:5],
        "segments": len(segs),
    }
    for sid in bad_seq:  # unparseable segments are excluded, having already failed
        del segs[sid]

    # Partition
    occ = Counter()
    seg_kmers = {}
    for sid, s in segs.items():
        ks = [pack(s[i : i + k]) for i in range(len(s) - k + 1)]
        seg_kmers[sid] = ks
        for v in ks:
            occ[canon(v, k)] += 1

    # Tags, where LN is the length, KC the summed true abundance and km is KC per k-mer
    tag_bad = []
    for sid, ks in seg_kmers.items():
        tg = tags.get(sid, {})
        want_kc = sum(counts.get(canon(v, k), 0) for v in ks)
        try:
            ln_ok = int(tg.get("LN", -1)) == len(segs[sid])
            kc_ok = int(tg.get("KC", -1)) == want_kc
            km_ok = abs(float(tg.get("km", "nan")) - want_kc / len(ks)) <= 0.005 + 1e-9
        except ValueError:
            ln_ok = kc_ok = km_ok = False
        if not (ln_ok and kc_ok and km_ok):
            tag_bad.append((sid, tg, {"LN": len(segs[sid]), "KC": want_kc}))
    checks["tags_exact"] = {
        "pass": not tag_bad,
        "bad_segments": len(tag_bad),
        "examples": [str(x)[:160] for x in tag_bad[:3]],
    }
    got = set(occ)
    checks["kmer_set_equal"] = {
        "pass": got == solid,
        "only_gfa": len(got - solid),
        "only_oracle": len(solid - got),
    }
    dups = sum(1 for n in occ.values() if n > 1)
    checks["no_duplicates"] = {"pass": dups == 0, "duplicated_kmers": dups}

    # Path validity, interior compaction and maximality
    bad_path = bad_interior = non_max = 0
    for ks in seg_kmers.values():
        for a, b in itertools.pairwise(ks):
            if not has_edge(a, b):
                bad_path += 1
            elif compactable(a) != b:
                bad_interior += 1
        first, last = ks[0], ks[-1]
        # The right end is last -> nxt when compactable, unless nxt closes a cycle onto first.
        nxt = compactable(last)
        if nxt is not None and canon(nxt, k) != canon(first, k):
            non_max += 1
        # The left end is rc(first) -> prv when compactable, unless it closes onto last.
        prv = compactable(rc(first, k))
        if prv is not None and canon(prv, k) != canon(last, k):
            non_max += 1
    checks["paths_valid"] = {"pass": bad_path == 0, "bad_edges": bad_path}
    checks["interior_compact"] = {"pass": bad_interior == 0, "violations": bad_interior}
    checks["maximal"] = {"pass": non_max == 0, "extendable_ends": non_max}

    # Links
    def norm(a, oa, b, ob):
        """Orientation-independent key of a link, the smaller of it and its reverse."""
        fwd = ((a, oa), (b, ob))
        flip = {"+": "-", "-": "+"}
        rev = ((b, flip[ob]), (a, flip[oa]))
        return min(fwd, rev)

    entry = {}
    for sid, ks in seg_kmers.items():
        entry.setdefault(ks[0], []).append((sid, "+"))
        entry.setdefault(rc(ks[-1], k), []).append((sid, "-"))
    want = set()
    for sid, ks in seg_kmers.items():
        for o, cur in (("+", ks[-1]), ("-", rc(ks[0], k))):
            for nxt in succs(cur):
                for tid, to in entry.get(nxt, []):
                    want.add(norm(sid, o, tid, to))
    have = set()
    bad_ovl = 0
    for a, oa, b, ob, ovl in links:
        if ovl != f"{k - 1}M":
            bad_ovl += 1
        have.add(norm(a, oa, b, ob))
    checks["links_exact"] = {
        "pass": have == want and bad_ovl == 0,
        "n_links": len(have),
        "missing": len(want - have),
        "spurious": len(have - want),
        "bad_overlap": bad_ovl,
    }

    # Optional -s node table (OKKDBG)
    if kdbg:
        checks["kdbg_exact"] = validate_kdbg(kdbg, k, solid, counts, has_edge)

    # Optional partition parity with Cuttlefish (--read).
    # Cuttlefish filters (k+1)-mer edges by count while optik filters k-mers, so the graphs coincide only at min_count 1.
    if args.cuttlefish and os.path.isfile(args.cuttlefish) and args.min_count == 1:
        with tempfile.TemporaryDirectory() as td:
            pref = os.path.join(td, "cf")
            p = subprocess.run(
                [
                    args.cuttlefish,
                    "build",
                    "--read",
                    "-s",
                    args.input,
                    "-k",
                    str(k),
                    "-c",
                    str(args.min_count),
                    "-o",
                    pref,
                    "-t",
                    str(args.threads),
                    "-w",
                    td,
                ],
                capture_output=True,
                check=False,
            )
            cf = pref + ".fa"
            if p.returncode != 0 or not os.path.isfile(cf):
                checks["cuttlefish_partition"] = {
                    "pass": False,
                    "error": p.stderr.decode()[-400:],
                }
            else:
                # Exact partition, each unitig taken as the set of its canonical k-mers
                def sig(s):
                    """Set of canonical k-mers of sequence s."""
                    return frozenset(
                        canon(pack(s[i : i + k]), k) for i in range(len(s) - k + 1)
                    )

                ours = Counter(sig(s) for s in segs.values())
                cf_seqs = [s.upper() for _, s in read_fasta(cf)]
                theirs = Counter(sig(s) for s in cf_seqs)
                only_o, only_c = ours - theirs, theirs - ours
                # Every disagreement is settled with the independent oracle.
                # A Cuttlefish-only unitig is Cuttlefish's error if one of its internal joins is not compactable in the read edge set, which makes it a merge across a branch.
                # An optik-only unitig is explained if its k-mers are input k-mers Cuttlefish's output dropped, or if they lie inside such an invalid Cuttlefish unitig.
                cf_bad_kmers, cf_branch_cross = set(), 0
                cf_all = set()
                for sq in cf_seqs:
                    ks = [pack(sq[i : i + k]) for i in range(len(sq) - k + 1)]
                    cks = frozenset(canon(v, k) for v in ks)
                    cf_all |= cks
                    if cks in only_c and any(
                        compactable(a) != b for a, b in itertools.pairwise(ks)
                    ):
                        cf_branch_cross += 1
                        cf_bad_kmers |= cks
                unexplained_c = sum(
                    n for u, n in only_c.items() if not (u <= cf_bad_kmers)
                )
                dropped = solid - cf_all
                # Only k-mers never linked to a neighbour (isolated runs of length k between non-ACGT characters) are an accepted Cuttlefish omission.
                dropped_isolated = {x for x in dropped if x not in linked}
                unexplained_o = sum(
                    n
                    for u, n in only_o.items()
                    if not all(x in dropped_isolated or x in cf_bad_kmers for x in u)
                )
                optik_ok = all(
                    checks[c]["pass"]
                    for c in (
                        "kmer_set_equal",
                        "no_duplicates",
                        "paths_valid",
                        "interior_compact",
                        "maximal",
                    )
                )
                checks["cuttlefish_partition"] = {
                    "pass": optik_ok and unexplained_c == 0 and unexplained_o == 0,
                    "cmd": f"cuttlefish build --read -s {args.input} -k {k} -c {args.min_count}",
                    "optik_unitigs": sum(ours.values()),
                    "cuttlefish_unitigs": sum(theirs.values()),
                    "identical_unitigs": sum((ours & theirs).values()),
                    "only_optik": sum(only_o.values()),
                    "only_cuttlefish": sum(only_c.values()),
                    "cuttlefish_branch_crossings": cf_branch_cross,
                    "cuttlefish_dropped_kmers": len(dropped),
                    "cuttlefish_dropped_isolated": len(dropped_isolated),
                    "unexplained_optik": unexplained_o,
                    "unexplained_cuttlefish": unexplained_c,
                    "detail": "partitions compared exactly; each difference must be a "
                    "Cuttlefish unitig crossing a branch of the read edge set or "
                    "an isolated (N-bounded, k-length) k-mer Cuttlefish omits",
                }
    return checks


def validate_kdbg(path, k, solid, counts, has_edge):
    """Check the dbruijn -s node table against the oracle.

    The OKKDBG layout is little-endian, "OKKDBG" | u8 k | u64 n | n x u64 sorted canonical k-mers | n x u32 abundances | n x u8 masks (in<<4 | out)."""
    import struct

    data = Path(path).read_bytes()
    if data[:6] != b"OKKDBG" or len(data) < 15:
        return {"pass": False, "error": f"bad header {data[:7]!r}"}
    kk, n = data[6], struct.unpack_from("<Q", data, 7)[0]
    if kk != k or len(data) != 15 + 13 * n:
        return {
            "pass": False,
            "error": f"k={kk} n={n} size={len(data)} (want k={k}, {15 + 13 * n} bytes)",
        }
    keys = struct.unpack_from(f"<{n}Q", data, 15)
    ab = struct.unpack_from(f"<{n}I", data, 15 + 8 * n)
    masks = data[15 + 12 * n :]
    kmask = (1 << (2 * k)) - 1
    sorted_ok = all(a < b for a, b in itertools.pairwise(keys))
    canon_ok = all(canon(c, k) == c for c in keys)
    set_ok = set(keys) == solid
    ab_bad = sum(1 for c, a in zip(keys, ab, strict=False) if counts.get(c, 0) != a)
    mask_bad = 0
    for c, m in zip(keys, masks, strict=False):
        out = sum(1 << b for b in range(4) if has_edge(c, ((c << 2) | b) & kmask))
        # The in-edge P -> c is the out-edge rc(c) -> rc(P).
        # has_edge checks the target's solidity, so this also requires the predecessor P to be solid (with -c > 1, predecessors below the cutoff are not nodes).
        inn = sum(
            1 << b
            for b in range(4)
            if has_edge(rc(c, k), rc((b << (2 * (k - 1))) | (c >> 2), k))
        )
        mask_bad += m != ((inn << 4) | out)
    return {
        "pass": sorted_ok and canon_ok and set_ok and ab_bad == 0 and mask_bad == 0,
        "nodes": n,
        "sorted": sorted_ok,
        "canonical": canon_ok,
        "set_equal": set_ok,
        "abundance_mismatches": ab_bad,
        "mask_mismatches": mask_bad,
    }


def main():
    """Check one GFA against its input and print the result."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-g", "--gfa", required=True)
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("-c", "--min-count", type=int, default=1)
    ap.add_argument("--cuttlefish", default=None)
    ap.add_argument("--kdbg", default=None, help="also validate the -s node table")
    ap.add_argument("-t", "--threads", type=int, default=4)
    args = ap.parse_args()
    checks = validate(
        args.input,
        args.gfa,
        args.k,
        args.min_count,
        args.cuttlefish,
        args.threads,
        kdbg=args.kdbg,
    )
    print(json.dumps(checks, indent=2))
    fails = [n for n, c in checks.items() if not c["pass"]]
    print("RESULT:", "FAIL " + ",".join(fails) if fails else "PASS", file=sys.stderr)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
