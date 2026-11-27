#!/usr/bin/env python3
"""Validate optik's mapping against ground truth and minimap2.

Checks
  1.  Ground-truth position.
      The primary lands at the simulated position (same reference, read origin within 1 bp, NM == 0 for error-free reads), or on an exact repeat copy of it.
  2.  Ground-truth strand.
      FLAG 16 matches a simulated '-' strand.
  2b. Record structure.
      Each read has one primary, with consistent SA tags.
  2c. Record audit.
      Every mapped record's SEQ, CIGAR, span and NM agree with the read and the reference.
  3.  CIGAR consistency.
      The CIGAR consumes the read length, and NM covers the edits it implies.
  4.  Parity with minimap2.
      Scored by the same truth rule, optik places at least as many reads correctly.
  4a. Parity with wfmash.
      wfmash runs at its defaults and read settings.
  4b. PAF round trip.
      wgatools validate finds no invalid record, and pafcov covers every mapped interval.
  5.  Reproducibility.
      Runs at two thread counts give byte-identical SAM bodies.
  6.  Mapping rate check.

Requirements
    - the optik binary
    - the minimap2 binary (optional, for the cross-tool checks)
    - gen_reads.py with --gt support, which generates the ground truth

Usage
    python validate_mapping.py -r ref.fa -k 31 [--minimap2 /path/to/minimap2]

A JSON report goes to stdout and a human-readable summary to stderr.
"""

import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import hw
from diagnostics import DETAIL_CHARS
from dna import read_fasta
from sam_spec import (
    SAM_MIN_COLUMNS,
    SAM_TAG_PREFIX_LEN,
    PafCol,
    SamCol,
    SamFlag,
)
from units import PERCENT

# Reads of CHAIN_MIN_LEN or more take optik's chaining path (tools/src/map.rs).
# That path writes one primary record, whose flanks may be soft-clipped, plus read-disjoint supplementary records (flag 2048) carrying SA tags.
# The primary is then validated anchored at its origin, pos - 1 == gt + leading_S.
LONG_QUERY_LEN = 500
# Reads at least this long are mapped as long reads (minimap2 map-ont, the "long" speed row).
# Shorter reads use the short-read preset (sr) instead.
LONG_READ_MIN_LEN = 1000
# timeouts in seconds: read simulation, one optik map, one minimap2 run, each wgatools step, and a speed run
SIMULATE_TIMEOUT_S = 120
MAP_TIMEOUT_S = 600
MINIMAP2_TIMEOUT_S = 1800
WGATOOLS_TIMEOUT_S = 300
SPEED_TIMEOUT_S = 7200
# examples kept per check, and characters of a CIGAR shown in one
MAX_EXAMPLES = 5
MAX_MISMATCH_EXAMPLES = 10
AUDIT_CIGAR_CHARS = 30
MISMATCH_CIGAR_CHARS = 40
# end of the previous interval before the first one of a contig, below every coordinate
NO_INTERVAL_END = -1
# CLI defaults
DEFAULT_NUM_READS = 10000
DEFAULT_READ_LEN = 150
DEFAULT_SEED = 42


# The determinism rerun uses one thread, or this many when the main run already used one.
ALT_THREADS = 2


def sanitize_ref(ref_path, out_dir):
    """Rewrite reference headers with explicit names.

    SSHash's stitched unitigs have bare '>' headers.
    optik and gen_reads.py both name such records seq_N, but minimap2 names them differently, so every read would mismatch in cross-tool comparisons.
    Each header becomes its first whitespace token, or seq_N when empty, so all tools, minimap2 included, see the same sequence names.
    """
    opener = gzip.open if str(ref_path).endswith(".gz") else open
    mode = "rt" if str(ref_path).endswith(".gz") else "r"
    out_path = os.path.join(out_dir, "ref_clean.fa")
    seq_idx = 0
    with opener(ref_path, mode) as in_handle, open(out_path, "w") as out_handle:
        for raw in in_handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                out_handle.write(f">{name}\n")
            else:
                out_handle.write(line + "\n")
    return out_path


def parse_sam(sam_text):
    """Parse SAM text into one record per read.

    Returns read_name -> {flag, ref, pos (1-based), mapq, cigar, seq, nm, is_mapped, is_reverse}.
    """
    records = {}
    for line in sam_text.split("\n"):
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < SAM_MIN_COLUMNS:
            continue
        read_name = parts[SamCol.QNAME]
        flag = int(parts[SamCol.FLAG])
        is_unmapped = bool(flag & SamFlag.UNMAPPED)
        is_reverse = bool(flag & SamFlag.REVERSE)

        rec = {
            "flag": flag,
            "is_mapped": not is_unmapped,
            "is_reverse": is_reverse,
            "ref": parts[SamCol.RNAME] if not is_unmapped else None,
            "pos": int(parts[SamCol.POS]) if not is_unmapped else 0,
            "mapq": int(parts[SamCol.MAPQ]),
            "cigar": parts[SamCol.CIGAR],
            "seq": parts[SamCol.SEQ],
        }

        # NM:i tag
        nm_tag = None
        for field in parts[SamCol.TAGS:]:
            if field.startswith("NM:i:"):
                try:
                    nm_tag = int(field[SAM_TAG_PREFIX_LEN:])
                except ValueError:
                    pass  # a malformed NM tag leaves nm_tag unset
                break
        rec["nm"] = nm_tag

        # Only the primary record of a read is kept, skipping secondary and supplementary ones.
        if flag & (SamFlag.SECONDARY | SamFlag.SUPPLEMENTARY):
            continue
        if read_name not in records:
            records[read_name] = rec
    return records


def leading_s(cigar):
    """Length of the leading soft-clip ('S') run in a CIGAR, else 0."""
    match = re.match(r"^(\d+)S", cigar)
    return int(match.group(1)) if match else 0


def parse_ground_truth(tsv_path):
    """Parse the ground-truth TSV that gen_reads.py --gt writes.

    Returns read_name -> {ref_name, pos (0-based), strand, length}.
    """
    truth = {}
    with open(tsv_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            truth[row["read_name"]] = {
                "ref_name": row["ref_name"],
                "pos": int(row["pos"]),
                "strand": row["strand"],
                "length": int(row["read_length"]),
            }
    return truth


def validate_cigar(cigar_str, read_len, nm_tag):
    """Check that a CIGAR fits the read length and the edit distance.

    In SAM, M means "alignment match" and covers both matches and mismatches.
    Only I and D are structural edits visible in the CIGAR.
    NM also counts mismatches inside M operations, which the operation string shows only when = and X are used.

    Returns (pass, detail_string).
    """
    if cigar_str == "*":
        return True, "unmapped (CIGAR *)"

    # CIGAR as (count, op) pairs
    import re

    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)
    if not ops:
        return False, f"invalid CIGAR: {cigar_str}"

    # Query length consumed and structural edits implied
    query_consumed = 0
    structural_edits = 0
    for count_str, cigar_op in ops:
        count = int(count_str)
        if cigar_op in "M=X":
            query_consumed += count
        elif cigar_op == "I":
            query_consumed += count
            structural_edits += count
        elif cigar_op == "D":
            structural_edits += count
        elif cigar_op == "S":
            query_consumed += count
        elif cigar_op == "H":
            pass  # a hard clip consumes no query bases

    if query_consumed != read_len:
        return False, (
            f"CIGAR consumes {query_consumed} query bases, read has {read_len}"
        )

    # NM:i must be at least the structural edits, since it also counts mismatches inside M.
    if nm_tag is not None and nm_tag < structural_edits:
        return False, (
            f"CIGAR implies {structural_edits} structural edits, "
            f"NM:i={nm_tag} (must be >= structural edits)"
        )

    return True, "CIGAR consistent"


def generate_reads_with_gt(ref_path, kmer_len, num_reads, read_len, seed, out_dir):
    """Generate reads and their ground truth with gen_reads.py."""
    gen_script = Path(__file__).resolve().parent / "gen_reads.py"
    reads_fa = os.path.join(out_dir, "reads.fa")
    gt_tsv = os.path.join(out_dir, "truth.tsv")

    cmd = [
        sys.executable,
        str(gen_script),
        "-r",
        ref_path,
        "-n",
        str(num_reads),
        "-l",
        str(read_len),
        "-s",
        str(seed),
        "-o",
        reads_fa,
        "--gt",
        gt_tsv,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=SIMULATE_TIMEOUT_S, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"gen_reads.py failed: {proc.stderr[:DETAIL_CHARS]}")
    return reads_fa, gt_tsv


def run_optik_map(ref_path, reads_fa, kmer_len, optik_bin, out_dir, threads=None):
    """Run optik map and return SAM text."""
    sam_path = os.path.join(out_dir, "optik.sam")
    # optik runs with its own defaults for every read type, as minimap2 runs its preset.
    cmd = [optik_bin, "map", "-r", ref_path, "-q", reads_fa, "-o", sam_path]
    if threads:
        cmd.extend(["-t", str(threads)])
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=MAP_TIMEOUT_S, check=False)
    if proc.returncode != 0:
        return None, f"optik map failed: {proc.stderr[:DETAIL_CHARS]}"
    with open(sam_path) as handle:
        return handle.read(), None


def run_minimap2(
    ref_path, reads_fa, kmer_len, minimap2_bin, out_dir, preset="sr", threads=None
):
    """Run minimap2 with the given preset on `threads` threads, default hw.cpus(), and return (SAM text, error).

    The caller picks `sr` for short reads and `map-ont` for long reads.
    `kmer_len` is not passed to minimap2.
    """
    # The documented preset alone, with no parameter overrides
    cmd = [
        minimap2_bin,
        "-a",
        "-x",
        preset,
        "-t",
        str(threads or hw.cpus()),
        ref_path,
        reads_fa,
    ]
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=MINIMAP2_TIMEOUT_S, check=False
    )
    if proc.returncode != 0:
        return None, f"minimap2 failed: {proc.stderr[:DETAIL_CHARS]}"
    return proc.stdout, None


def load_ref_seqs(ref_path):
    """Reference sequences as a dict of name -> sequence."""
    return {name: seq for name, seq in read_fasta(ref_path)}


# PAF round trip through wgatools
#
# wgatools 1.0 and later dropped the `sam2paf` subcommand (the toolkit is now oriented to MAF and chain files), so optik's SAM is converted to PAF here.
# The 12 core PAF columns are a direct projection of SAM.
# The result then goes through two independent wgatools checks.
#
#   * `wgatools validate` recomputes query and target coordinates from the cg:Z (CIGAR) tag and flags every record whose position disagrees with its CIGAR.
#     Zero invalid records is a strong end-to-end check of optik's alignment bookkeeping.
#
#   * `wgatools pafcov` gives the per-base depth of the reference, summarized here as covered reference bases and mean depth.
#     Every read is drawn from the reference, so all mapped intervals must be covered, at a mean depth of at least 1, for the PAF to be consistent.


def sam_to_paf(sam_text, ref_lens, dropped=None):
    """Convert primary mapped SAM records to PAF lines carrying cg:Z and NM:i.

    Records that cannot be converted (a CIGAR that does not consume the read, or an empty alignment) are appended to `dropped` when a list is given.
    Skipping them silently would hide exactly the bookkeeping errors the wgatools checks are there to catch."""
    import re as _re

    paf_lines = []
    if dropped is None:
        dropped = []
    for line in sam_text.split("\n"):
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < SAM_MIN_COLUMNS:
            continue
        flag = int(parts[SamCol.FLAG])
        # Unmapped, secondary and supplementary records are skipped.
        # PAF cannot represent them, and the read's primary is the record under validation.
        if flag & (SamFlag.UNMAPPED | SamFlag.SECONDARY | SamFlag.SUPPLEMENTARY):
            continue
        ref = parts[SamCol.RNAME]
        qlen = len(parts[SamCol.SEQ]) if parts[SamCol.SEQ] != "*" else 0
        if qlen == 0:
            dropped.append(parts[SamCol.QNAME])
            continue
        tlen = ref_lens.get(ref, 0)
        tstart = int(parts[SamCol.POS]) - 1
        ops = _re.findall(r"(\d+)([MIDNSHP=X])", parts[SamCol.CIGAR])
        if not ops:
            dropped.append(parts[SamCol.QNAME])
            continue
        q_consumed = 0
        t_consumed = 0
        nmatch = 0
        # The SAM spec allows soft clips only at the ends, and PAF's cg:Z has no S operation.
        # They are stripped and the aligned query interval is reported instead of 0..qlen, with the leading S bases giving the alignment's offset in the read.
        q_lead = 0
        q_tail = 0
        aligned_ops = []
        for count_str, cigar_op in ops:
            count = int(count_str)
            if cigar_op in "M=X":
                q_consumed += count
                t_consumed += count
                nmatch += count
                aligned_ops.append((count_str, cigar_op))
            elif cigar_op == "I":
                q_consumed += count
                aligned_ops.append((count_str, cigar_op))
            elif cigar_op == "D":
                t_consumed += count
                aligned_ops.append((count_str, cigar_op))
            elif cigar_op == "S":
                q_consumed += count
        if q_consumed != qlen:
            # The CIGAR does not consume the read, which PAF cannot represent and which is a bookkeeping error, so it is reported rather than silently skipped.
            dropped.append(parts[SamCol.QNAME])
            continue
        # Leading and trailing S come from the original operation list (end-anchored only).
        first_pad = 0
        for count_str, cigar_op in ops:
            if cigar_op == "S":
                first_pad += int(count_str)
            else:
                break
        last_pad = 0
        for count_str, cigar_op in reversed(ops):
            if cigar_op == "S":
                last_pad += int(count_str)
            else:
                break
        q_lead, q_tail = first_pad, last_pad
        if not aligned_ops:
            dropped.append(parts[SamCol.QNAME])
            continue
        # PAF query coordinates refer to the read as given, while SAM SEQ and so its clips are reverse-complemented for flag 16.
        if flag & SamFlag.REVERSE:
            q_start, q_end = q_tail, qlen - q_lead
        else:
            q_start, q_end = q_lead, qlen - q_tail
        if q_end < q_start or t_consumed == 0:
            dropped.append(parts[SamCol.QNAME])
            continue
        cigar = "".join(f"{op_len}{cigar_op}" for op_len, cigar_op in aligned_ops)
        tend = tstart + t_consumed
        strand = "-" if (flag & SamFlag.REVERSE) else "+"
        blocklen = q_consumed + t_consumed - nmatch - (q_lead + q_tail)
        nm_tag = 0
        for field in parts[SamCol.TAGS:]:
            if field.startswith("NM:i:"):
                try:
                    nm_tag = int(field[SAM_TAG_PREFIX_LEN:])
                except ValueError:
                    pass  # a malformed NM tag leaves nm_tag unset
                break
        paf_lines.append(
            f"{parts[SamCol.QNAME]}\t{qlen}\t{q_start}\t{q_end}\t{strand}\t{ref}\t{tlen}\t"
            f"{tstart}\t{tend}\t{nmatch}\t{blocklen}\t{parts[SamCol.MAPQ]}\t"
            f"cg:Z:{cigar}\tnm:i:{nm_tag}"
        )
    return paf_lines


def audit_records(sam_text, reads, ref_seqs):
    """Audit every mapped record, primary and supplementary, against the input reads and the reference.

      seq      SEQ equals the read, reverse-complemented for FLAG 16.
               "*" or a missing read counts as bad.
      cigar    the CIGAR is well formed and consumes exactly len(SEQ) query bases.
      span     the reference span lies inside the contig.
      nm       NM:i is present and equals the edit count recomputed from SEQ, CIGAR and reference, counting mismatches, N, inserted and deleted bases.

    Returns a dict of counts and a few examples."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    out = {
        "records": 0,
        "seq_bad": 0,
        "cigar_bad": 0,
        "span_bad": 0,
        "nm_bad": 0,
        "examples": [],
    }
    for line in sam_text.split("\n"):
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < SAM_MIN_COLUMNS or int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
            continue
        out["records"] += 1
        flag, seq, cig = int(fields[SamCol.FLAG]), fields[SamCol.SEQ], fields[SamCol.CIGAR]
        # SEQ must be the input read, reverse-complemented on reverse-strand records.
        read = reads.get(fields[SamCol.QNAME])
        bad = []
        if (
            read is None
            or seq == "*"
            or seq.upper()
            != (read.translate(comp)[::-1] if flag & SamFlag.REVERSE else read).upper()
        ):
            out["seq_bad"] += 1
            bad.append("seq")
        # The CIGAR must parse completely and consume exactly as many bases as SEQ holds.
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
        if (
            not ops
            or "".join(op_len + cigar_op for op_len, cigar_op in ops) != cig
            or sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MIS=X") != len(seq)
        ):
            out["cigar_bad"] += 1
            bad.append("cigar")
        # The alignment must lie inside its reference sequence.
        ref = ref_seqs.get(fields[SamCol.RNAME])
        pos0 = int(fields[SamCol.POS]) - 1
        span = sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MDN=X")
        if ref is None or pos0 < 0 or pos0 + span > len(ref):
            out["span_bad"] += 1
            bad.append("span")
        nm_tag = next((field[SAM_TAG_PREFIX_LEN:] for field in fields[SamCol.TAGS:] if field.startswith("NM:i:")), None)
        # NM is recomputed from the bases only when the CIGAR and span are sound, counting N on the reference as a mismatch.
        if "cigar" not in bad and "span" not in bad:
            read_pos, ref_pos, edits = 0, pos0, 0
            for op_len, cigar_op in ops:
                op_len = int(op_len)
                if cigar_op in "M=X":
                    edits += sum(
                        1
                        for j in range(op_len)
                        if seq[read_pos + j].upper() != ref[ref_pos + j].upper()
                        or ref[ref_pos + j].upper() not in "ACGT"
                    )
                    read_pos += op_len
                    ref_pos += op_len
                elif cigar_op in "IS":
                    edits += op_len if cigar_op == "I" else 0
                    read_pos += op_len
                elif cigar_op in "DN":
                    edits += op_len if cigar_op == "D" else 0
                    ref_pos += op_len
            if nm_tag is None or not nm_tag.isdigit() or int(nm_tag) != edits:
                out["nm_bad"] += 1
                bad.append(f"nm {nm_tag} != {edits}")
        elif nm_tag is None:
            out["nm_bad"] += 1
        if bad and len(out["examples"]) < MAX_EXAMPLES:
            out["examples"].append(f"{fields[SamCol.QNAME]} flag={flag} {fields[SamCol.RNAME]}:{fields[SamCol.POS]} {cig[:AUDIT_CIGAR_CHARS]} {bad}")
    return out


def run_wgatools_paf_checks(paf_path, wgatools, tmpdir):
    """Run wgatools validate and pafcov on a PAF file.

    Returns (report_dict, error_str).
    """
    import re as _re

    # wgatools validate checks position and CIGAR consistency.
    proc = subprocess.run(
        [wgatools, "validate", paf_path],
        capture_output=True,
        text=True,
        timeout=WGATOOLS_TIMEOUT_S,
        check=False,
    )
    if proc.returncode != 0 and (proc.stderr or proc.stdout):
        # v1.1.0 panics with rc 101 on some inputs.
        # That counts as a tool failure, so it is reported instead of passing silently.
        return None, f"wgatools validate failed: {(proc.stderr or proc.stdout)[:DETAIL_CHARS]}"
    query_invalid = 0
    target_invalid = 0
    total_records = 0
    for line in proc.stdout.splitlines():
        match = _re.match(r"Total records:\s+(\d+)", line.strip())
        if match:
            total_records = int(match.group(1))
            continue
        match = _re.match(r"Query invalid records:\s+(\d+)", line.strip())
        if match:
            query_invalid = int(match.group(1))
            continue
        match = _re.match(r"Target invalid records:\s+(\d+)", line.strip())
        if match:
            target_invalid = int(match.group(1))
            continue

    # wgatools pafcov gives the per-base reference depth.
    proc2 = subprocess.run(
        [wgatools, "pafcov", paf_path, "-o", "-"],
        capture_output=True,
        text=True,
        timeout=WGATOOLS_TIMEOUT_S,
        check=False,
    )
    if proc2.returncode != 0 and (proc2.stderr or proc2.stdout):
        return None, f"wgatools pafcov failed: {(proc2.stderr or proc2.stdout)[:DETAIL_CHARS]}"
    covered_bases = 0
    depth_sum = 0
    total_ref_bases = 0
    # wgatools pafcov prints one "contig<TAB>start<TAB>end<TAB>depth" line per reference base, so each line counts as one base.
    for line in proc2.stdout.splitlines():
        try:
            _contig, _start, _end, depth_text, *_ = line.split("\t")
        except ValueError:  # fewer than four columns
            continue
        total_ref_bases += 1
        try:
            depth = int(depth_text)
        except ValueError:
            continue
        if depth >= 1:
            covered_bases += 1
            depth_sum += depth
    mean_depth = round(depth_sum / covered_bases, 3) if covered_bases else 0.0

    report = {
        "validate_total_records": total_records,
        "validate_query_invalid": query_invalid,
        "validate_target_invalid": target_invalid,
        "pafcov_covered_ref_bases": covered_bases,
        "pafcov_total_ref_bases": total_ref_bases,
        "pafcov_mean_depth": mean_depth,
    }
    return report, None


def _exact_repeat_copy(rec, ref_seqs):
    """True if the read sits on an exact duplicate copy of its locus.

    A read mapped with NM == 0 away from its simulated position is a legitimate alternate-copy placement when the reference at the reported 1-based POS spells the read in either orientation.
    The sequence then occurs identically at the true and the reported locus, so the two are a real repeat.
    """
    ref_seq = ref_seqs.get(rec["ref"])
    if not ref_seq:
        return False
    seq = rec.get("seq") or ""
    pos0 = rec["pos"] - 1
    if pos0 < 0 or pos0 + len(seq) > len(ref_seq):
        return False
    if ref_seq[pos0 : pos0 + len(seq)] == seq:
        return True
    rev_comp = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return ref_seq[pos0 : pos0 + len(seq)] == rev_comp


def validate(args):
    """Run every mapping check and return the report dict."""
    checks = {}
    kmer_len = args.k
    optik = args.optik
    tmpdir = tempfile.mkdtemp(prefix="vm_")

    try:
        # Reference headers are normalized so that optik, minimap2 and gen_reads see the same sequence names.
        ref_clean = sanitize_ref(args.ref, tmpdir)

        # Reads and ground truth
        reads_fa, gt_tsv = generate_reads_with_gt(
            ref_clean, kmer_len, args.n, args.read_len, args.seed, tmpdir
        )
        ground_truth = parse_ground_truth(gt_tsv)
        # Every requested read must exist, since with zero reads each count-based check below would pass as 0 == 0.
        checks["reads_generated"] = {
            "pass": len(ground_truth) == args.n > 0,
            "reads": len(ground_truth),
            "requested": args.n,
        }
        if not checks["reads_generated"]["pass"]:
            return checks

        # optik map
        start = time.monotonic()
        sam_text, err = run_optik_map(
            ref_clean, reads_fa, kmer_len, optik, tmpdir, threads=args.threads
        )
        optik_seconds = round(time.monotonic() - start, 3)
        if err:
            checks["optik_map"] = {"pass": False, "error": err}
            return checks
        optik_recs = parse_sam(sam_text)

        # Check 1, ground-truth position
        gt_match = 0
        gt_mismatch = 0
        gt_unmapped = 0
        gt_repeat_copy = 0
        mismatches = []
        ref_seqs = load_ref_seqs(ref_clean)
        long_mode = args.read_len >= LONG_QUERY_LEN
        for read_name, gt_info in ground_truth.items():
            rec = optik_recs.get(read_name)
            if not rec or not rec["is_mapped"]:
                gt_unmapped += 1
                continue
            # The primary record is anchored at its origin.
            # The aligned core starts on the read's simulated diagonal (long reads may soft-clip flanks), and error-free reads must align with NM == 0.
            if (
                rec["ref"] == gt_info["ref_name"]
                and rec["pos"] - 1 == gt_info["pos"] + leading_s(rec["cigar"])
                and rec["is_reverse"] == (gt_info["strand"] == "-")
                and rec["nm"] == 0
            ):
                gt_match += 1
            elif rec["nm"] == 0 and _exact_repeat_copy(rec, ref_seqs):
                # An identical repeat copy is an equally valid placement.
                gt_repeat_copy += 1
            else:
                gt_mismatch += 1
                if len(mismatches) < MAX_MISMATCH_EXAMPLES:
                    mismatches.append(
                        {
                            "read": read_name,
                            "gt_ref": gt_info["ref_name"],
                            "gt_pos": gt_info["pos"],
                            "gt_strand": gt_info["strand"],
                            "optik_ref": rec["ref"],
                            "optik_pos": rec["pos"],
                            "optik_reverse": rec["is_reverse"],
                            "cigar": rec["cigar"][:MISMATCH_CIGAR_CHARS],
                            "nm": rec["nm"],
                        }
                    )

        total_gt = len(ground_truth)
        # rate_pct is reported for reading only, since the pass condition already demands zero mismatches and zero unmapped reads.
        gt_rate = (gt_match + gt_repeat_copy) / total_gt * PERCENT if total_gt else 0
        checks["gt_position_match"] = {
            "pass": gt_mismatch == 0 and gt_unmapped == 0,
            "match": gt_match,
            "repeat_copy_accepted": gt_repeat_copy,
            "mismatch": gt_mismatch,
            "unmapped": gt_unmapped,
            "total": total_gt,
            "rate_pct": round(gt_rate, 2),
            "mismatches_sample": mismatches[:MAX_EXAMPLES],
            "detail": (
                "error-free reads: primary at pos-1 == gt + leading_S, "
                "truth strand, NM 0 (identical repeat copies accepted)"
            ),
            "long_mode": long_mode,
        }

        # Check 2, ground-truth strand.
        # gen_reads.py samples both strands, and FLAG 16 must coincide with a truth strand of '-'.
        wrong_strand = 0
        for read_name, gt_info in ground_truth.items():
            rec = optik_recs.get(read_name)
            if (
                rec
                and rec["is_mapped"]
                and rec["is_reverse"] != (gt_info["strand"] == "-")
            ):
                wrong_strand += 1
        checks["gt_strand_match"] = {
            "pass": wrong_strand == 0,
            "wrong_strand": wrong_strand,
            "total_mapped": sum(1 for rec in optik_recs.values() if rec["is_mapped"]),
            "reverse_in_truth": sum(
                1 for truth_rec in ground_truth.values() if truth_rec["strand"] == "-"
            ),
            "detail": "primary strand (FLAG 16) must equal the simulated strand",
        }

        # Check 2b, SAM record structure with one primary and consistent SA tags
        n_prim, n_supp, sa_bad = {}, 0, 0
        for line in sam_text.split("\n"):
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            flag = int(fields[SamCol.FLAG])
            if flag & SamFlag.SECONDARY:
                continue
            if flag & SamFlag.SUPPLEMENTARY:
                n_supp += 1
                if not any(field.startswith("SA:Z:") for field in fields[SamCol.TAGS:]):
                    sa_bad += 1
                continue
            n_prim[fields[SamCol.QNAME]] = n_prim.get(fields[SamCol.QNAME], 0) + 1
        dup = sum(1 for count in n_prim.values() if count > 1)
        checks["one_primary_per_read"] = {
            "pass": dup == 0 and len(n_prim) == len(ground_truth) and sa_bad == 0,
            "reads_with_primary": len(n_prim),
            "reads": len(ground_truth),
            "duplicate_primary": dup,
            "supplementary": n_supp,
            "supplementary_without_SA": sa_bad,
        }

        # Check 2c, an independent audit of every mapped record.
        # SEQ must be the read, the CIGAR must consume it, the span must fit the contig, and NM must equal the edits recomputed against the reference.
        reads_seq = {name: seq for name, seq in read_fasta(reads_fa)}
        audit = audit_records(sam_text, reads_seq, ref_seqs)
        checks["record_audit"] = {
            "pass": audit["records"] > 0
            and not (
                audit["seq_bad"] or audit["cigar_bad"] or audit["span_bad"] or audit["nm_bad"]
            ),
            **audit,
            "detail": "every mapped record: SEQ == read, CIGAR consumes SEQ, span "
            "inside the contig, NM == edits re-derived vs the reference",
        }

        # Check 3, CIGAR consistency
        cigar_pass = 0
        cigar_fail = 0
        cigar_details = []
        for read_name, rec in optik_recs.items():
            if not rec["is_mapped"]:
                continue
            read_seq = ground_truth.get(read_name, {}).get("length", len(rec["seq"]))
            passed, detail = validate_cigar(rec["cigar"], read_seq, rec["nm"])
            if passed:
                cigar_pass += 1
            else:
                cigar_fail += 1
                if len(cigar_details) < MAX_EXAMPLES:
                    cigar_details.append(
                        {
                            "read": read_name,
                            "cigar": rec["cigar"],
                            "nm": rec["nm"],
                            "detail": detail,
                        }
                    )
        checks["cigar_consistency"] = {
            "pass": cigar_fail == 0,
            "consistent": cigar_pass,
            "inconsistent": cigar_fail,
            "details_sample": cigar_details,
            "detail": "CIGAR must consume correct query length and match NM:i",
        }

        # Check 4, parity with minimap2 under the simulation truth.
        # The two tools map identical reads with their own documented defaults, and one function scores both primaries against the truth (reference, strand, and read origin pos - leading_S within 1 bp, or an exact repeat copy).
        # optik must place at least as many reads correctly.
        def truth_ok(rec, truth):
            """True when record rec places the read at its truth gt or on an exact repeat copy."""
            if not rec or not rec["is_mapped"]:
                return False
            if (
                rec["ref"] == truth["ref_name"]
                and rec["is_reverse"] == (truth["strand"] == "-")
                and abs(rec["pos"] - 1 - leading_s(rec["cigar"]) - truth["pos"]) <= 1
            ):
                return True
            return rec["nm"] == 0 and _exact_repeat_copy(rec, ref_seqs)

        speed = {"optik": optik_seconds}
        if args.minimap2 and os.path.isfile(args.minimap2):
            mm2_preset = "map-ont" if args.read_len >= LONG_READ_MIN_LEN else "sr"
            start = time.monotonic()
            mm2_sam, err = run_minimap2(
                ref_clean,
                reads_fa,
                kmer_len,
                args.minimap2,
                tmpdir,
                preset=mm2_preset,
                threads=args.threads,
            )
            speed["minimap2"] = round(time.monotonic() - start, 3)
            if err:
                checks["minimap2_truth_parity"] = {"pass": False, "error": err}
            else:
                mm2_recs = parse_sam(mm2_sam)
                # A baseline that did not report every read (a crash, or empty or truncated output) is no comparison, and "optik >= 0" would be vacuous.
                mm2_missing = sum(1 for name in ground_truth if name not in mm2_recs)
                o_ok = sum(
                    truth_ok(optik_recs.get(name), truth_rec) for name, truth_rec in ground_truth.items()
                )
                m_ok = sum(
                    truth_ok(mm2_recs.get(name), truth_rec) for name, truth_rec in ground_truth.items()
                )
                agree = sum(
                    1
                    for name in ground_truth
                    if optik_recs.get(name)
                    and mm2_recs.get(name)
                    and optik_recs[name]["is_mapped"]
                    and mm2_recs[name]["is_mapped"]
                    and optik_recs[name]["ref"] == mm2_recs[name]["ref"]
                    and abs(optik_recs[name]["pos"] - mm2_recs[name]["pos"]) <= 1
                )
                # These reads have no errors, so any optik primary with a nonzero NM is misaligned.
                o_nm = sum(
                    1
                    for rec in optik_recs.values()
                    if rec["is_mapped"] and rec["nm"] not in (0, None)
                )
                checks["minimap2_truth_parity"] = {
                    "pass": o_ok >= m_ok
                    and o_nm == 0
                    and mm2_missing == 0
                    and m_ok > 0,
                    "minimap2_reads_without_record": mm2_missing,
                    "cmd": f"minimap2 -a -x {mm2_preset} REF READS",
                    "optik_correct": o_ok,
                    "minimap2_correct": m_ok,
                    "reads": len(ground_truth),
                    "same_primary_pos": agree,
                    "optik_primary_nm_nonzero": o_nm,
                    "detail": "error-free reads: optik correct >= minimap2 correct "
                    "(same truth scoring) and every optik primary NM == 0",
                }

        # Check 4a, parity with the wavefront mapper wfmash scored against the truth.
        # It runs with defaults and with read settings (validate_official.WFMASH_READ_OPTS), and optik must place at least as many reads correctly as the better of the two.
        if getattr(args, "wfmash", None) and os.path.isfile(args.wfmash):
            from validate_official import WFMASH_READ_OPTS, faidx, one_primary_per_read

            faidx(ref_clean)
            faidx(reads_fa)
            kind = "short" if args.read_len < LONG_READ_MIN_LEN else "long"
            wres = {}
            for cfg, opts in (("defaults", []), ("read_opts", WFMASH_READ_OPTS[kind])):
                raw = os.path.join(tmpdir, f"wfmash_{cfg}.raw.sam")
                one = os.path.join(tmpdir, f"wfmash_{cfg}.sam")
                start = time.monotonic()
                with open(raw, "w") as handle:
                    wfmash_proc = subprocess.run(
                        [
                            args.wfmash,
                            ref_clean,
                            reads_fa,
                            "-a",
                            "-t",
                            str(args.threads),
                            *opts,
                        ],
                        stdout=handle,
                        stderr=subprocess.DEVNULL,
                        timeout=SPEED_TIMEOUT_S,
                        check=False,
                    )
                elapsed = round(time.monotonic() - start, 3)
                speed[f"wfmash_{cfg}"] = elapsed
                if wfmash_proc.returncode != 0:
                    wres[cfg] = {"error": f"exit {wfmash_proc.returncode}", "seconds": elapsed}
                    continue
                one_primary_per_read(raw, one)
                with open(one) as handle:
                    wrecs = parse_sam(handle.read())
                wres[cfg] = {
                    "correct": sum(
                        truth_ok(wrecs.get(name), truth_rec) for name, truth_rec in ground_truth.items()
                    ),
                    "seconds": elapsed,
                    "cmd": " ".join(["wfmash", "REF", "READS", "-a", *opts]),
                }
            o_ok = sum(truth_ok(optik_recs.get(name), truth_rec) for name, truth_rec in ground_truth.items())
            best = max((result.get("correct", 0) for result in wres.values()), default=0)
            # A baseline that crashed or failed is no comparison, so it is reported as a failure rather than scored as 0 correct.
            broken = {cfg_name: result["error"] for cfg_name, result in wres.items() if "error" in result}
            wfmash_check = {
                "optik_correct": o_ok,
                "reads": len(ground_truth),
                "wfmash": wres,
                "detail": "optik correct >= wavefront mapper (wfmash) correct, same truth scoring",
            }
            if broken:
                wfmash_check.update(
                    pass_=False,
                    error=f"wfmash baseline failed: {broken} (a negative "
                    "exit is a signal, e.g. -11 = SIGSEGV)",
                )
            elif best == 0:
                # wfmash ran but placed no read in any configuration, since it is not built for reads this short.
                # "optik >= 0" proves nothing, so this is reported rather than counted as a pass.
                wfmash_check.update(info="vacuous: wfmash placed 0 reads; not a comparison")
            else:
                wfmash_check.update(pass_=o_ok >= best)
            if "pass_" in wfmash_check:
                wfmash_check["pass"] = wfmash_check.pop("pass_")
            checks["wfmash_truth_parity"] = wfmash_check
        checks["mapper_speed"] = {
            "info": speed,
            "detail": "whole-process wall seconds, same reads",
        }

        # Check 4b, PAF round trip through wgatools.
        # optik's SAM becomes PAF with cg:Z and NM:i tags, as the helpers above write it, and goes through two independent wgatools checks.
        # `validate` recomputes query and target positions from the CIGAR, where any invalid record is a bookkeeping bug, and `pafcov` gives per-base reference depth.
        # Every read is drawn from the reference, so all mapped intervals must be covered and the mean per-base depth must be at least 1.
        if getattr(args, "wgatools", None):
            ref_lens = {name: len(seq) for name, seq in ref_seqs.items()}
            dropped = []
            paf_lines = sam_to_paf(sam_text, ref_lens, dropped)
            paf_path = os.path.join(tmpdir, "optik.paf")
            with open(paf_path, "w") as handle:
                handle.write("\n".join(paf_lines) + ("\n" if paf_lines else ""))
            wg_report, wg_err = run_wgatools_paf_checks(paf_path, args.wgatools, tmpdir)
            if wg_err:
                checks["wgatools_paf_roundtrip"] = {"pass": False, "error": wg_err}
            else:
                # pafcov must cover exactly the union of the PAF target intervals, computed here independently.
                intervals = {}
                for text_line in paf_lines:
                    fields = text_line.split("\t")
                    intervals.setdefault(fields[PafCol.TNAME], []).append((int(fields[PafCol.TSTART]), int(fields[PafCol.TEND])))
                union = 0
                for spans in intervals.values():
                    end = NO_INTERVAL_END
                    for span_start, span_end in sorted(spans):
                        span_start = max(span_start, end)
                        if span_end > span_start:
                            union += span_end - span_start
                            end = span_end
                n_prim = sum(1 for rec in optik_recs.values() if rec["is_mapped"])
                passed = (
                    not dropped
                    and len(paf_lines) == n_prim
                    and wg_report["validate_total_records"] == len(paf_lines)
                    and wg_report["validate_query_invalid"] == 0
                    and wg_report["validate_target_invalid"] == 0
                    # pafcov must count exactly the union computed here, proving the PAF coordinates are right.
                    and wg_report["pafcov_covered_ref_bases"] == union
                )
                checks["wgatools_paf_roundtrip"] = {
                    "pass": passed,
                    "paf_records": len(paf_lines),
                    "mapped_primaries": n_prim,
                    "unconvertible_records": len(dropped),
                    "expected_covered_bases": union,
                    **wg_report,
                    "detail": (
                        "every mapped primary converts to PAF; wgatools validate "
                        "parses all of them with 0 invalid; pafcov covers exactly "
                        "the union of the aligned target intervals"
                    ),
                }

        # Check 5, reproducibility.
        # The second run uses another thread count, 1 or else 2, so the SAM must not depend on it.
        # ref_clean is reused, because the main mapping above sanitized the headers (SSHash's stitched unitigs have bare '>' headers) and the second run must see the same reference byte for byte, or the SAM bodies would legitimately differ.
        other_t = 1 if args.threads != 1 else ALT_THREADS
        sam_text_2, err = run_optik_map(
            ref_clean, reads_fa, kmer_len, optik, tmpdir, threads=other_t
        )
        if err:
            checks["reproducibility"] = {"pass": False, "error": err}
        else:
            # The SAM header is dropped, since its thread-dependent timing lines may differ.
            def strip_header(text):
                """SAM text without its @ header lines."""
                lines = text.split("\n")
                body = [
                    text_line for text_line in lines if not text_line.startswith("@")
                ]
                return "\n".join(body)

            body1 = strip_header(sam_text)
            body2 = strip_header(sam_text_2)
            identical = body1 == body2
            checks["reproducibility"] = {
                "pass": identical,
                "detail": "SAM body should be byte-identical across runs",
            }

        # Check 6, mapping rate
        mapped = sum(1 for rec in optik_recs.values() if rec["is_mapped"])
        total = len(optik_recs)
        rate = mapped / total * PERCENT if total else 0
        checks["mapping_rate"] = {
            "pass": mapped == len(ground_truth),
            "mapped": mapped,
            "total": total,
            "rate_pct": round(rate, 2),
            "detail": "every error-free read drawn from the reference must map",
        }

    finally:
        import shutil

        shutil.rmtree(tmpdir, ignore_errors=True)

    return checks


def main():
    """Entry point that runs every mapping check and reports the result."""
    parser = argparse.ArgumentParser(description="Validate optik mapping correctness")
    parser.add_argument("-r", dest="ref", required=True, help="Reference FASTA")
    parser.add_argument("-k", type=int, required=True, help="K-mer size")
    parser.add_argument(
        "--optik", default="target/release/optik", help="Path to optik binary"
    )
    parser.add_argument(
        "--minimap2", default=None, help="Path to minimap2 binary (optional)"
    )
    parser.add_argument(
        "--wfmash",
        default=None,
        help="Path to wfmash (optional; wavefront mapper parity + speed)",
    )
    parser.add_argument(
        "--wgatools",
        default=None,
        help="Path to wgatools binary (optional; enables PAF round-trip checks)",
    )
    parser.add_argument(
        "-n", type=int, default=DEFAULT_NUM_READS, help="Number of reads to generate"
    )
    parser.add_argument("--read-len", type=int, default=DEFAULT_READ_LEN, help="Read length")
    parser.add_argument("-s", "--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument(
        "--threads",
        type=int,
        default=hw.cpus(),
        help="threads of every mapper (default: the usable CPUs)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.optik):
        print(f"ERROR: optik binary not found at {args.optik}", file=sys.stderr)
        sys.exit(1)

    checks = validate(args)

    all_pass = all(
        check.get("pass", False)
        for check in checks.values()
        if isinstance(check, dict) and "pass" in check
    )
    overall = "PASS" if all_pass else "FAIL"

    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  Mapping Validation: {overall}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    for name, check in checks.items():
        if "pass" not in check and "info" in check:
            status = "INFO"  # measurement, not a gate
        else:
            status = "PASS" if check.get("pass", False) else "FAIL"
        line = f"  [{status}] {name}"
        if "rate_pct" in check:
            line += f"  ({check['rate_pct']}%)"
        if "detail" in check:
            line += f"  -- {check['detail']}"
        if check.get("error"):
            line += f"  ERROR: {check['error']}"
        print(line, file=sys.stderr)
    print(file=sys.stderr)

    report = {
        "overall": overall,
        "checks": checks,
    }
    print(json.dumps(report, indent=2))

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
