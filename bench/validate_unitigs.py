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

from diagnostics import DETAIL_CHARS
from dna import read_fasta
from sam_spec import (
    KDBG_ABUNDANCE_BYTES,
    KDBG_COUNT_OFFSET,
    KDBG_HEADER_BYTES,
    KDBG_IN_MASK_SHIFT,
    KDBG_K_OFFSET,
    KDBG_KEY_BYTES,
    KDBG_MAGIC,
    KDBG_NODE_BYTES,
    SAM_TAG_NAME_LEN,
    SAM_TAG_PREFIX_LEN,
    GfaLinkCol,
    GfaSegmentCol,
)

# defaults of the CLI and of validate()
DEFAULT_K = 31
DEFAULT_THREADS = 4
# offending segment names listed per well-formedness failure, and tag mismatches shown (each cut to EXAMPLE_CHARS)
MAX_LISTED = 5
MAX_EXAMPLES = 3
EXAMPLE_CHARS = 160
# value of a missing LN or KC tag, which never matches
MISSING_TAG = -1
# km is written with 2 decimals, so it may differ from KC / k-mers by up to half the last digit (plus float slack)
KM_TOLERANCE = 0.005
FLOAT_SLACK = 1e-9

ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
# Packed k-mers hold two bits per base, so BASE_MASK selects one base and BASE_MASK - code complements it.
BITS_PER_BASE = 2
BASE_MASK = (1 << BITS_PER_BASE) - 1
NUM_BASES = 1 << BITS_PER_BASE
# revcomp_packed() works on the 8 bytes of a packed word, four bases to a byte.
BYTE_VALUES = 256
KMER_WORD_BYTES = 8
BASES_PER_BYTE = 8 // BITS_PER_BASE
MAX_PACKED_K = KMER_WORD_BYTES * BASES_PER_BYTE


def windows(seq, kmer_len):
    """Yield the packed forward k-mers of length `kmer_len` of seq, skipping any that contain a non-ACGT base."""
    mask = (1 << (BITS_PER_BASE * kmer_len)) - 1
    packed = 0
    valid = 0
    for char in seq.upper():
        code = ENC.get(char)
        if code is None:
            valid = 0
            packed = 0
            continue
        packed = ((packed << BITS_PER_BASE) | code) & mask
        valid += 1
        if valid >= kmer_len:
            yield packed


# _RC_BYTE[b] is byte b (four 2-bit bases) complemented, with its bases in reverse order.
# revcomp_packed() maps the 8 bytes of a packed w-mer (w <= 32) through it and reads them back in reverse byte order, which yields the whole reverse complement at C speed.
# A loop per base made the oracle take minutes per genome.
_RC_BYTE = bytes(
    sum(((BASE_MASK - ((b >> (BITS_PER_BASE * i)) & BASE_MASK)) << (BITS_PER_BASE * (BASES_PER_BYTE - 1 - i))) for i in range(BASES_PER_BYTE))
    for b in range(BYTE_VALUES)
)


def revcomp_packed(packed, kmer_len):
    """Reverse complement of the packed k-mer `packed` (2 bits per base, kmer_len <= 32).

    The 32 - kmer_len unused high bases complement into the low end and are shifted out."""
    return int.from_bytes(packed.to_bytes(KMER_WORD_BYTES, "little").translate(_RC_BYTE), "big") >> (
        BITS_PER_BASE * (MAX_PACKED_K - kmer_len)
    )


def _rc_slow(packed, kmer_len):
    """Base-by-base reverse complement, the reference test_harness checks revcomp_packed() against."""
    result = 0
    for _ in range(kmer_len):
        result = (result << BITS_PER_BASE) | (BASE_MASK - (packed & BASE_MASK))
        packed >>= BITS_PER_BASE
    return result


def canon(packed, kmer_len):
    """Canonical form of the packed k-mer `packed`, the smaller of it and its reverse complement."""
    rev = revcomp_packed(packed, kmer_len)
    return min(rev, packed)


_DIGITS = str.maketrans("ACGT", "0123")


def pack(seq):
    """Packed value of an ACGT string (A=0 to T=3), i.e. its base-4 number.

    Callers pass only ACGT strings, since segments are checked before packing."""
    return int(seq.translate(_DIGITS), NUM_BASES) if seq else 0


def _pack_slow(seq):
    """Base-by-base pack(), the reference test_harness checks pack() against."""
    packed = 0
    for char in seq:
        packed = (packed << BITS_PER_BASE) | ENC[char]
    return packed


def parse_gfa(path):
    """Parse a GFA into segments {name: SEQ}, links, S-line tags {name: {tag: value}} and the names that occur on more than one S line."""
    segs, links, tags, dup_names = {}, [], {}, []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                if fields[GfaSegmentCol.NAME] in segs:
                    dup_names.append(fields[GfaSegmentCol.NAME])
                segs[fields[GfaSegmentCol.NAME]] = fields[GfaSegmentCol.SEQ].upper()
                tags[fields[GfaSegmentCol.NAME]] = {
                    field[:SAM_TAG_NAME_LEN]: field[SAM_TAG_PREFIX_LEN:]
                    for field in fields[GfaSegmentCol.TAGS:]
                    if len(field) > SAM_TAG_PREFIX_LEN and field[SAM_TAG_NAME_LEN] == ":"
                }
            elif fields[0] == "L":
                links.append(
                    (
                        fields[GfaLinkCol.FROM],
                        fields[GfaLinkCol.FROM_ORIENT],
                        fields[GfaLinkCol.TO],
                        fields[GfaLinkCol.TO_ORIENT],
                        fields[GfaLinkCol.OVERLAP],
                    )
                )
    return segs, links, tags, dup_names


def oracle(input_path, kmer_len, min_count=1):
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
        for packed in windows(seq, kmer_len):
            counts[canon(packed, kmer_len)] += 1
        for edge in windows(seq, kmer_len + 1):
            kp1.add(edge)
            linked.add(canon(edge >> BITS_PER_BASE, kmer_len))
            linked.add(canon(edge & ((1 << (BITS_PER_BASE * kmer_len)) - 1), kmer_len))
    solid = {kmer for kmer, count in counts.items() if count >= min_count}
    return solid, kp1, linked, counts


def validate(
    input_path,
    gfa_path,
    kmer_len,
    min_count=1,
    cuttlefish=None,
    threads=DEFAULT_THREADS,
    orc=None,
    kdbg=None,
):
    """Run every structural check and return {check_name: {"pass": bool, ...}}."""
    args = argparse.Namespace(
        input=input_path,
        gfa=gfa_path,
        k=kmer_len,
        min_count=min_count,
        cuttlefish=cuttlefish,
        threads=threads,
    )
    kmask = (1 << (BITS_PER_BASE * kmer_len)) - 1

    solid, kp1, linked, counts = (
        orc if orc is not None else oracle(args.input, kmer_len, args.min_count)
    )

    def has_edge(cur, nxt):
        """True when oriented edge cur -> nxt, with nxt = cur[1:] + b, is in the edge-centric graph."""
        edge = (cur << BITS_PER_BASE) | (nxt & BASE_MASK)
        return (edge in kp1 or revcomp_packed(edge, kmer_len + 1) in kp1) and canon(nxt, kmer_len) in solid

    def succs(cur):
        """Successor k-mers of cur along edges of the graph."""
        out = []
        for b in range(NUM_BASES):
            nxt = ((cur << BITS_PER_BASE) | b) & kmask
            if has_edge(cur, nxt):
                out.append(nxt)
        return out

    def indeg(nxt):
        """In-degree of the k-mer nxt, counted as the successors of its reverse complement."""
        return len(succs(revcomp_packed(nxt, kmer_len)))

    def compactable(cur):
        """The k-mer cur compacts into (unique successor with in-degree 1), or None."""
        successors = succs(cur)
        if len(successors) != 1:
            return None
        nxt = successors[0]
        if canon(nxt, kmer_len) == canon(cur, kmer_len):
            return None
        if cur == revcomp_packed(cur, kmer_len) or nxt == revcomp_packed(
            nxt, kmer_len
        ):  # palindromes, possible at even k, stay single
            return None
        return nxt if indeg(nxt) == 1 else None

    segs, links, tags, dup_names = parse_gfa(args.gfa)
    checks = {}

    # Segment well-formedness, checked before anything indexes by name or k-mer
    bad_seq = [sid for sid, seq in segs.items() if len(seq) < kmer_len or set(seq) - set("ACGT")]
    checks["segments_wellformed"] = {
        "pass": not dup_names and not bad_seq and bool(segs) == bool(solid),
        "duplicate_names": dup_names[:MAX_LISTED],
        "short_or_non_acgt": bad_seq[:MAX_LISTED],
        "segments": len(segs),
    }
    for sid in bad_seq:  # unparseable segments are excluded, having already failed
        del segs[sid]

    # Partition
    occ = Counter()
    seg_kmers = {}
    for sid, seq in segs.items():
        seg_packed = [pack(seq[i : i + kmer_len]) for i in range(len(seq) - kmer_len + 1)]
        seg_kmers[sid] = seg_packed
        for packed in seg_packed:
            occ[canon(packed, kmer_len)] += 1

    # Tags, where LN is the length, KC the summed true abundance and km is KC per k-mer
    tag_bad = []
    for sid, seg_packed in seg_kmers.items():
        seg_tags = tags.get(sid, {})
        want_kc = sum(counts.get(canon(packed, kmer_len), 0) for packed in seg_packed)
        try:
            ln_ok = int(seg_tags.get("LN", MISSING_TAG)) == len(segs[sid])
            kc_ok = int(seg_tags.get("KC", MISSING_TAG)) == want_kc
            km_ok = abs(float(seg_tags.get("km", "nan")) - want_kc / len(seg_packed)) <= KM_TOLERANCE + FLOAT_SLACK
        except ValueError:
            ln_ok = kc_ok = km_ok = False
        if not (ln_ok and kc_ok and km_ok):
            tag_bad.append((sid, seg_tags, {"LN": len(segs[sid]), "KC": want_kc}))
    checks["tags_exact"] = {
        "pass": not tag_bad,
        "bad_segments": len(tag_bad),
        "examples": [str(bad)[:EXAMPLE_CHARS] for bad in tag_bad[:MAX_EXAMPLES]],
    }
    got = set(occ)
    checks["kmer_set_equal"] = {
        "pass": got == solid,
        "only_gfa": len(got - solid),
        "only_oracle": len(solid - got),
    }
    dups = sum(1 for count in occ.values() if count > 1)
    checks["no_duplicates"] = {"pass": dups == 0, "duplicated_kmers": dups}

    # Path validity, interior compaction and maximality
    bad_path = bad_interior = non_max = 0
    for seg_packed in seg_kmers.values():
        for left, right in itertools.pairwise(seg_packed):
            if not has_edge(left, right):
                bad_path += 1
            elif compactable(left) != right:
                bad_interior += 1
        first, last = seg_packed[0], seg_packed[-1]
        # The right end is last -> nxt when compactable, unless nxt closes a cycle onto first.
        nxt = compactable(last)
        if nxt is not None and canon(nxt, kmer_len) != canon(first, kmer_len):
            non_max += 1
        # The left end is rc(first) -> prv when compactable, unless it closes onto last.
        prv = compactable(revcomp_packed(first, kmer_len))
        if prv is not None and canon(prv, kmer_len) != canon(last, kmer_len):
            non_max += 1
    checks["paths_valid"] = {"pass": bad_path == 0, "bad_edges": bad_path}
    checks["interior_compact"] = {"pass": bad_interior == 0, "violations": bad_interior}
    checks["maximal"] = {"pass": non_max == 0, "extendable_ends": non_max}

    # Links
    def norm(from_seg, from_orient, to_seg, to_orient):
        """Orientation-independent key of a link, the smaller of it and its reverse."""
        fwd = ((from_seg, from_orient), (to_seg, to_orient))
        flip = {"+": "-", "-": "+"}
        rev = ((to_seg, flip[to_orient]), (from_seg, flip[from_orient]))
        return min(fwd, rev)

    entry = {}
    for sid, seg_packed in seg_kmers.items():
        entry.setdefault(seg_packed[0], []).append((sid, "+"))
        entry.setdefault(revcomp_packed(seg_packed[-1], kmer_len), []).append((sid, "-"))
    want = set()
    for sid, seg_packed in seg_kmers.items():
        for orient, cur in (("+", seg_packed[-1]), ("-", revcomp_packed(seg_packed[0], kmer_len))):
            for nxt in succs(cur):
                for tid, to_orient in entry.get(nxt, []):
                    want.add(norm(sid, orient, tid, to_orient))
    have = set()
    bad_ovl = 0
    for from_seg, from_orient, to_seg, to_orient, ovl in links:
        if ovl != f"{kmer_len - 1}M":
            bad_ovl += 1
        have.add(norm(from_seg, from_orient, to_seg, to_orient))
    checks["links_exact"] = {
        "pass": have == want and bad_ovl == 0,
        "n_links": len(have),
        "missing": len(want - have),
        "spurious": len(have - want),
        "bad_overlap": bad_ovl,
    }

    # Optional -s node table (OKKDBG)
    if kdbg:
        checks["kdbg_exact"] = validate_kdbg(kdbg, kmer_len, solid, counts, has_edge)

    # Optional partition parity with Cuttlefish (--read).
    # Cuttlefish filters (k+1)-mer edges by count while optik filters k-mers, so the graphs coincide only at min_count 1.
    if args.cuttlefish and os.path.isfile(args.cuttlefish) and args.min_count == 1:
        with tempfile.TemporaryDirectory() as tmp_dir:
            pref = os.path.join(tmp_dir, "cf")
            proc = subprocess.run(
                [
                    args.cuttlefish,
                    "build",
                    "--read",
                    "-s",
                    args.input,
                    "-k",
                    str(kmer_len),
                    "-c",
                    str(args.min_count),
                    "-o",
                    pref,
                    "-t",
                    str(args.threads),
                    "-w",
                    tmp_dir,
                ],
                capture_output=True,
                check=False,
            )
            cf_fa = pref + ".fa"
            if proc.returncode != 0 or not os.path.isfile(cf_fa):
                checks["cuttlefish_partition"] = {
                    "pass": False,
                    "error": proc.stderr.decode()[-DETAIL_CHARS:],
                }
            else:
                # Exact partition, each unitig taken as the set of its canonical k-mers
                def sig(seq):
                    """Set of canonical k-mers of `seq`."""
                    return frozenset(
                        canon(pack(seq[i : i + kmer_len]), kmer_len) for i in range(len(seq) - kmer_len + 1)
                    )

                ours = Counter(sig(seq) for seq in segs.values())
                cf_seqs = [seq.upper() for _, seq in read_fasta(cf_fa)]
                theirs = Counter(sig(seq) for seq in cf_seqs)
                only_o, only_c = ours - theirs, theirs - ours
                # Every disagreement is settled with the independent oracle.
                # A Cuttlefish-only unitig is Cuttlefish's error if one of its internal joins is not compactable in the read edge set, which makes it a merge across a branch.
                # An optik-only unitig is explained if its k-mers are input k-mers Cuttlefish's output dropped, or if they lie inside such an invalid Cuttlefish unitig.
                cf_bad_kmers, cf_branch_cross = set(), 0
                cf_all = set()
                for cf_seq in cf_seqs:
                    seg_packed = [pack(cf_seq[i : i + kmer_len]) for i in range(len(cf_seq) - kmer_len + 1)]
                    cks = frozenset(canon(packed, kmer_len) for packed in seg_packed)
                    cf_all |= cks
                    if cks in only_c and any(
                        compactable(left) != right for left, right in itertools.pairwise(seg_packed)
                    ):
                        cf_branch_cross += 1
                        cf_bad_kmers |= cks
                unexplained_c = sum(
                    count for unitig, count in only_c.items() if not (unitig <= cf_bad_kmers)
                )
                dropped = solid - cf_all
                # Only k-mers never linked to a neighbour (isolated runs of length k between non-ACGT characters) are an accepted Cuttlefish omission.
                dropped_isolated = {kmer for kmer in dropped if kmer not in linked}
                unexplained_o = sum(
                    count
                    for unitig, count in only_o.items()
                    if not all(kmer in dropped_isolated or kmer in cf_bad_kmers for kmer in unitig)
                )
                # Explaining differences only means something once optik's own unitigs pass every structural check.
                optik_ok = all(
                    checks[name]["pass"]
                    for name in (
                        "kmer_set_equal",
                        "no_duplicates",
                        "paths_valid",
                        "interior_compact",
                        "maximal",
                    )
                )
                checks["cuttlefish_partition"] = {
                    "pass": optik_ok and unexplained_c == 0 and unexplained_o == 0,
                    "cmd": f"cuttlefish build --read -s {args.input} -k {kmer_len} -c {args.min_count}",
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


def validate_kdbg(path, kmer_len, solid, counts, has_edge):
    """Check the dbruijn -s node table against the oracle.

    The OKKDBG layout is little-endian, "OKKDBG" | u8 k | u64 n | n x u64 sorted canonical k-mers | n x u32 abundances | n x u8 masks (in<<4 | out)."""
    import struct

    data = Path(path).read_bytes()
    if data[:KDBG_K_OFFSET] != KDBG_MAGIC or len(data) < KDBG_HEADER_BYTES:
        return {"pass": False, "error": f"bad header {data[:KDBG_COUNT_OFFSET]!r}"}
    file_k, num_nodes = data[KDBG_K_OFFSET], struct.unpack_from("<Q", data, KDBG_COUNT_OFFSET)[0]
    if file_k != kmer_len or len(data) != KDBG_HEADER_BYTES + KDBG_NODE_BYTES * num_nodes:
        return {
            "pass": False,
            "error": f"k={file_k} n={num_nodes} size={len(data)} (want k={kmer_len}, {KDBG_HEADER_BYTES + KDBG_NODE_BYTES * num_nodes} bytes)",
        }
    keys = struct.unpack_from(f"<{num_nodes}Q", data, KDBG_HEADER_BYTES)
    abundances = struct.unpack_from(f"<{num_nodes}I", data, KDBG_HEADER_BYTES + KDBG_KEY_BYTES * num_nodes)
    masks = data[KDBG_HEADER_BYTES + (KDBG_KEY_BYTES + KDBG_ABUNDANCE_BYTES) * num_nodes :]
    kmask = (1 << (BITS_PER_BASE * kmer_len)) - 1
    sorted_ok = all(prev_key < next_key for prev_key, next_key in itertools.pairwise(keys))
    canon_ok = all(canon(key, kmer_len) == key for key in keys)
    set_ok = set(keys) == solid
    ab_bad = sum(1 for key, abundance in zip(keys, abundances, strict=False) if counts.get(key, 0) != abundance)
    mask_bad = 0
    for key, mask in zip(keys, masks, strict=False):
        out = sum(1 << b for b in range(NUM_BASES) if has_edge(key, ((key << BITS_PER_BASE) | b) & kmask))
        # The in-edge P -> c is the out-edge rc(c) -> rc(P).
        # has_edge checks the target's solidity, so this also requires the predecessor P to be solid (with -c > 1, predecessors below the cutoff are not nodes).
        inn = sum(
            1 << b
            for b in range(NUM_BASES)
            if has_edge(revcomp_packed(key, kmer_len), revcomp_packed((b << (BITS_PER_BASE * (kmer_len - 1))) | (key >> BITS_PER_BASE), kmer_len))
        )
        mask_bad += mask != ((inn << KDBG_IN_MASK_SHIFT) | out)
    return {
        "pass": sorted_ok and canon_ok and set_ok and ab_bad == 0 and mask_bad == 0,
        "nodes": num_nodes,
        "sorted": sorted_ok,
        "canonical": canon_ok,
        "set_equal": set_ok,
        "abundance_mismatches": ab_bad,
        "mask_mismatches": mask_bad,
    }


def main():
    """Check one GFA against its input and print the result."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-g", "--gfa", required=True)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument("-c", "--min-count", type=int, default=1)
    parser.add_argument("--cuttlefish", default=None)
    parser.add_argument("--kdbg", default=None, help="also validate the -s node table")
    parser.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    args = parser.parse_args()
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
    fails = [name for name, check in checks.items() if not check["pass"]]
    print("RESULT:", "FAIL " + ",".join(fails) if fails else "PASS", file=sys.stderr)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
