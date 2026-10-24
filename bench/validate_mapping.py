#!/usr/bin/env python3
"""Validate optik's mapping against ground truth and minimap2.

Checks
  1.  Ground-truth position.
      The primary lands at the simulated position (same reference, read origin within 1 bp, NM == 0 for error-free reads), or on an exact repeat copy of it.
  2.  Ground-truth strand.
      FLAG 16 matches a simulated '-' strand.
  2b. Record structure.  One primary per read, with consistent SA tags.
  2c. Record audit.  Every mapped record's SEQ, CIGAR, span and NM agree
      with the read and the reference.
  3.  CIGAR consistency.
      The CIGAR consumes the read length, and NM covers the edits it implies.
  4.  Parity with minimap2.
      Scored by the same truth rule, optik places at least as many reads correctly.
  4a. Parity with wfmash, at its defaults and read settings.
  4b. PAF round trip.  wgatools validate finds no invalid record and
      pafcov covers every mapped interval.
  5.  Reproducibility.
      Runs at two thread counts give byte-identical SAM bodies.
  6.  Mapping rate sanity.

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
from sam_spec import (
    SAM_FLAG_REVERSE,
    SAM_FLAG_SECONDARY,
    SAM_FLAG_SUPPLEMENTARY,
    SAM_FLAG_UNMAPPED,
    SAM_MIN_COLUMNS,
)

# Reads of CHAIN_MIN_LEN or more take optik's chaining path (tools/src/map.rs).
# That path writes one primary record, whose flanks may be soft-clipped, plus read-disjoint supplementary records (flag 2048) carrying SA tags.
# The primary is then validated anchored at its origin, pos - 1 == gt + leading_S.
LONG_QUERY_LEN = 500


def read_fasta(path):
    """List the (name, sequence) records of a FASTA file, naming headerless records seq_N."""
    sequences = []
    name = None
    parts = []
    seq_idx = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(parts)))
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                parts = []
            else:
                parts.append(line)
    if name is not None:
        sequences.append((name, "".join(parts)))
    return sequences


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
    with opener(ref_path, mode) as f, open(out_path, "w") as w:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                w.write(f">{name}\n")
            else:
                w.write(line + "\n")
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
        read_name = parts[0]
        flag = int(parts[1])
        is_unmapped = bool(flag & SAM_FLAG_UNMAPPED)
        is_reverse = bool(flag & SAM_FLAG_REVERSE)

        rec = {
            "flag": flag,
            "is_mapped": not is_unmapped,
            "is_reverse": is_reverse,
            "ref": parts[2] if not is_unmapped else None,
            "pos": int(parts[3]) if not is_unmapped else 0,
            "mapq": int(parts[4]),
            "cigar": parts[5],
            "seq": parts[9],
        }

        # NM:i tag
        nm = None
        for field in parts[11:]:
            if field.startswith("NM:i:"):
                try:
                    nm = int(field[5:])
                except ValueError:
                    pass
                break
        rec["nm"] = nm

        # Only the primary record of a read is kept, skipping secondary and supplementary ones.
        if flag & (SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY):
            continue
        if read_name not in records:
            records[read_name] = rec
    return records


def leading_s(cigar):
    """Length of the leading soft-clip ('S') run in a CIGAR, else 0."""
    m = re.match(r"^(\d+)S", cigar)
    return int(m.group(1)) if m else 0


def parse_ground_truth(tsv_path):
    """Parse the ground-truth TSV that gen_reads.py --gt writes.

    Returns read_name -> {ref_name, pos (0-based), strand, length}.
    """
    gt = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gt[row["read_name"]] = {
                "ref_name": row["ref_name"],
                "pos": int(row["pos"]),
                "strand": row["strand"],
                "length": int(row["read_length"]),
            }
    return gt


def validate_cigar(cigar_str, read_len, nm):
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
    for count_str, op in ops:
        count = int(count_str)
        if op in "M=X":
            query_consumed += count
        elif op == "I":
            query_consumed += count
            structural_edits += count
        elif op == "D":
            structural_edits += count
        elif op == "S":
            query_consumed += count
        elif op == "H":
            pass  # a hard clip consumes no query bases

    if query_consumed != read_len:
        return False, (
            f"CIGAR consumes {query_consumed} query bases, read has {read_len}"
        )

    # NM:i must be at least the structural edits, since it also counts mismatches inside M.
    if nm is not None and nm < structural_edits:
        return False, (
            f"CIGAR implies {structural_edits} structural edits, "
            f"NM:i={nm} (must be >= structural edits)"
        )

    return True, "CIGAR consistent"


def generate_reads_with_gt(ref_path, k, n, read_len, seed, out_dir):
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
        str(n),
        "-l",
        str(read_len),
        "-s",
        str(seed),
        "-o",
        reads_fa,
        "--gt",
        gt_tsv,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"gen_reads.py failed: {proc.stderr[:500]}")
    return reads_fa, gt_tsv


def run_optik_map(ref_path, reads_fa, k, optik_bin, out_dir, threads=None):
    """Run optik map and return SAM text."""
    sam_path = os.path.join(out_dir, "optik.sam")
    # optik runs with its own defaults for every read type, as minimap2 runs its preset.
    cmd = [optik_bin, "map", "-r", ref_path, "-q", reads_fa, "-o", sam_path]
    if threads:
        cmd.extend(["-t", str(threads)])
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=False)
    if proc.returncode != 0:
        return None, f"optik map failed: {proc.stderr[:500]}"
    with open(sam_path) as f:
        return f.read(), None


def run_minimap2(
    ref_path, reads_fa, k, minimap2_bin, out_dir, preset="sr", threads=None
):
    """Run minimap2 with the given preset on `threads` threads, default hw.cpus(), and return (SAM text, error).

    The caller picks `sr` for short reads and `map-ont` for long reads.
    `k` is not passed to minimap2.
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
        cmd, capture_output=True, text=True, timeout=1800, check=False
    )
    if proc.returncode != 0:
        return None, f"minimap2 failed: {proc.stderr[:500]}"
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
#     Every read is drawn from the reference, so all mapped intervals must be covered, at a mean depth of at least 1, for the PAF to be sane.


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
        flag = int(parts[1])
        # Unmapped, secondary and supplementary records are skipped.
        # PAF cannot represent them, and the read's primary is the record under validation.
        if flag & (SAM_FLAG_UNMAPPED | SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY):
            continue
        ref = parts[2]
        qlen = len(parts[9]) if parts[9] != "*" else 0
        if qlen == 0:
            dropped.append(parts[0])
            continue
        tlen = ref_lens.get(ref, 0)
        tstart = int(parts[3]) - 1
        ops = _re.findall(r"(\d+)([MIDNSHP=X])", parts[5])
        if not ops:
            dropped.append(parts[0])
            continue
        q_consumed = 0
        t_consumed = 0
        nmatch = 0
        # The SAM spec allows soft clips only at the ends, and PAF's cg:Z has no S operation.
        # They are stripped and the aligned query interval is reported instead of 0..qlen, with the leading S bases giving the alignment's offset in the read.
        q_lead = 0
        q_tail = 0
        aligned_ops = []
        for count_str, op in ops:
            count = int(count_str)
            if op in "M=X":
                q_consumed += count
                t_consumed += count
                nmatch += count
                aligned_ops.append((count_str, op))
            elif op == "I":
                q_consumed += count
                aligned_ops.append((count_str, op))
            elif op == "D":
                t_consumed += count
                aligned_ops.append((count_str, op))
            elif op == "S":
                q_consumed += count
        if q_consumed != qlen:
            # The CIGAR does not consume the read, which PAF cannot represent and which is a bookkeeping error, so it is reported rather than silently skipped.
            dropped.append(parts[0])
            continue
        # Leading and trailing S come from the original operation list (end-anchored only).
        first_pad = 0
        for count_str, op in ops:
            if op == "S":
                first_pad += int(count_str)
            else:
                break
        last_pad = 0
        for count_str, op in reversed(ops):
            if op == "S":
                last_pad += int(count_str)
            else:
                break
        q_lead, q_tail = first_pad, last_pad
        if not aligned_ops:
            dropped.append(parts[0])
            continue
        # PAF query coordinates refer to the read as given, while SAM SEQ and so its clips are reverse-complemented for flag 16.
        if flag & SAM_FLAG_REVERSE:
            q_start, q_end = q_tail, qlen - q_lead
        else:
            q_start, q_end = q_lead, qlen - q_tail
        if q_end < q_start or t_consumed == 0:
            dropped.append(parts[0])
            continue
        cg = "".join(f"{n}{op}" for n, op in aligned_ops)
        tend = tstart + t_consumed
        strand = "-" if (flag & SAM_FLAG_REVERSE) else "+"
        blocklen = q_consumed + t_consumed - nmatch - (q_lead + q_tail)
        nm = 0
        for field in parts[11:]:
            if field.startswith("NM:i:"):
                try:
                    nm = int(field[5:])
                except ValueError:
                    pass
                break
        paf_lines.append(
            f"{parts[0]}\t{qlen}\t{q_start}\t{q_end}\t{strand}\t{ref}\t{tlen}\t"
            f"{tstart}\t{tend}\t{nmatch}\t{blocklen}\t{parts[4]}\t"
            f"cg:Z:{cg}\tnm:i:{nm}"
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
        t = line.split("\t")
        if len(t) < SAM_MIN_COLUMNS or int(t[1]) & SAM_FLAG_UNMAPPED:
            continue
        out["records"] += 1
        flag, seq, cig = int(t[1]), t[9], t[5]
        read = reads.get(t[0])
        bad = []
        if (
            read is None
            or seq == "*"
            or seq.upper()
            != (read.translate(comp)[::-1] if flag & SAM_FLAG_REVERSE else read).upper()
        ):
            out["seq_bad"] += 1
            bad.append("seq")
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
        if (
            not ops
            or "".join(n + o for n, o in ops) != cig
            or sum(int(n) for n, o in ops if o in "MIS=X") != len(seq)
        ):
            out["cigar_bad"] += 1
            bad.append("cigar")
        ref = ref_seqs.get(t[2])
        pos0 = int(t[3]) - 1
        span = sum(int(n) for n, o in ops if o in "MDN=X")
        if ref is None or pos0 < 0 or pos0 + span > len(ref):
            out["span_bad"] += 1
            bad.append("span")
        nm = next((x[5:] for x in t[11:] if x.startswith("NM:i:")), None)
        if "cigar" not in bad and "span" not in bad:
            qi, ri, d = 0, pos0, 0
            for n, o in ops:
                n = int(n)
                if o in "M=X":
                    d += sum(
                        1
                        for j in range(n)
                        if seq[qi + j].upper() != ref[ri + j].upper()
                        or ref[ri + j].upper() not in "ACGT"
                    )
                    qi += n
                    ri += n
                elif o in "IS":
                    d += n if o == "I" else 0
                    qi += n
                elif o in "DN":
                    d += n if o == "D" else 0
                    ri += n
            if nm is None or not nm.isdigit() or int(nm) != d:
                out["nm_bad"] += 1
                bad.append(f"nm {nm} != {d}")
        elif nm is None:
            out["nm_bad"] += 1
        if bad and len(out["examples"]) < 5:
            out["examples"].append(f"{t[0]} flag={flag} {t[2]}:{t[3]} {cig[:30]} {bad}")
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
        timeout=300,
        check=False,
    )
    if proc.returncode != 0 and (proc.stderr or proc.stdout):
        # v1.1.0 panics with rc 101 on some inputs.
        # That counts as a tool failure, so it is reported instead of passing silently.
        return None, f"wgatools validate failed: {(proc.stderr or proc.stdout)[:300]}"
    query_invalid = 0
    target_invalid = 0
    total_records = 0
    for line in proc.stdout.splitlines():
        m = _re.match(r"Total records:\s+(\d+)", line.strip())
        if m:
            total_records = int(m.group(1))
            continue
        m = _re.match(r"Query invalid records:\s+(\d+)", line.strip())
        if m:
            query_invalid = int(m.group(1))
            continue
        m = _re.match(r"Target invalid records:\s+(\d+)", line.strip())
        if m:
            target_invalid = int(m.group(1))
            continue

    # wgatools pafcov gives the per-base reference depth.
    proc2 = subprocess.run(
        [wgatools, "pafcov", paf_path, "-o", "-"],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if proc2.returncode != 0 and (proc2.stderr or proc2.stdout):
        return None, f"wgatools pafcov failed: {(proc2.stderr or proc2.stdout)[:300]}"
    covered_bases = 0
    depth_sum = 0
    total_ref_bases = 0
    for line in proc2.stdout.splitlines():
        col = line.split("\t")
        if len(col) < 4:
            continue
        total_ref_bases += 1
        try:
            depth = int(col[3])
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
    The sequence then occurs identically at the true and the reported locus, a genuine repeat.
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
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return ref_seq[pos0 : pos0 + len(seq)] == rc


def validate(args):
    """Run every mapping check and return the report dict."""
    checks = {}
    k = args.k
    optik = args.optik
    tmpdir = tempfile.mkdtemp(prefix="vm_")

    try:
        # Reference headers are normalized so that optik, minimap2 and gen_reads see the same sequence names.
        ref_clean = sanitize_ref(args.ref, tmpdir)

        # Reads and ground truth
        reads_fa, gt_tsv = generate_reads_with_gt(
            ref_clean, k, args.n, args.read_len, args.seed, tmpdir
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
        t0 = time.monotonic()
        sam_text, err = run_optik_map(
            ref_clean, reads_fa, k, optik, tmpdir, threads=args.threads
        )
        optik_seconds = round(time.monotonic() - t0, 3)
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
                if len(mismatches) < 10:
                    mismatches.append(
                        {
                            "read": read_name,
                            "gt_ref": gt_info["ref_name"],
                            "gt_pos": gt_info["pos"],
                            "gt_strand": gt_info["strand"],
                            "optik_ref": rec["ref"],
                            "optik_pos": rec["pos"],
                            "optik_reverse": rec["is_reverse"],
                            "cigar": rec["cigar"][:40],
                            "nm": rec["nm"],
                        }
                    )

        total_gt = len(ground_truth)
        gt_rate = (gt_match + gt_repeat_copy) / total_gt * 100 if total_gt else 0
        checks["gt_position_match"] = {
            "pass": gt_mismatch == 0 and gt_unmapped == 0,
            "match": gt_match,
            "repeat_copy_accepted": gt_repeat_copy,
            "mismatch": gt_mismatch,
            "unmapped": gt_unmapped,
            "total": total_gt,
            "rate_pct": round(gt_rate, 2),
            "mismatches_sample": mismatches[:5],
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
            "total_mapped": sum(1 for r in optik_recs.values() if r["is_mapped"]),
            "reverse_in_truth": sum(
                1 for g in ground_truth.values() if g["strand"] == "-"
            ),
            "detail": "primary strand (FLAG 16) must equal the simulated strand",
        }

        # Check 2b, SAM record structure with one primary and consistent SA tags
        n_prim, n_supp, sa_bad = {}, 0, 0
        for line in sam_text.split("\n"):
            if not line or line.startswith("@"):
                continue
            t = line.split("\t")
            f = int(t[1])
            if f & SAM_FLAG_SECONDARY:
                continue
            if f & SAM_FLAG_SUPPLEMENTARY:
                n_supp += 1
                if not any(x.startswith("SA:Z:") for x in t[11:]):
                    sa_bad += 1
                continue
            n_prim[t[0]] = n_prim.get(t[0], 0) + 1
        dup = sum(1 for v in n_prim.values() if v > 1)
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
        reads_seq = {n: s for n, s in read_fasta(reads_fa)}
        au = audit_records(sam_text, reads_seq, ref_seqs)
        checks["record_audit"] = {
            "pass": au["records"] > 0
            and not (
                au["seq_bad"] or au["cigar_bad"] or au["span_bad"] or au["nm_bad"]
            ),
            **au,
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
            ok, detail = validate_cigar(rec["cigar"], read_seq, rec["nm"])
            if ok:
                cigar_pass += 1
            else:
                cigar_fail += 1
                if len(cigar_details) < 5:
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
        def truth_ok(rec, gt):
            """True when record rec places the read at its truth gt or on an exact repeat copy."""
            if not rec or not rec["is_mapped"]:
                return False
            if (
                rec["ref"] == gt["ref_name"]
                and rec["is_reverse"] == (gt["strand"] == "-")
                and abs(rec["pos"] - 1 - leading_s(rec["cigar"]) - gt["pos"]) <= 1
            ):
                return True
            return rec["nm"] == 0 and _exact_repeat_copy(rec, ref_seqs)

        speed = {"optik": optik_seconds}
        if args.minimap2 and os.path.isfile(args.minimap2):
            mm2_preset = "map-ont" if args.read_len >= 1000 else "sr"
            t0 = time.monotonic()
            mm2_sam, err = run_minimap2(
                ref_clean,
                reads_fa,
                k,
                args.minimap2,
                tmpdir,
                preset=mm2_preset,
                threads=args.threads,
            )
            speed["minimap2"] = round(time.monotonic() - t0, 3)
            if err:
                checks["minimap2_truth_parity"] = {"pass": False, "error": err}
            else:
                mm2_recs = parse_sam(mm2_sam)
                # A baseline that did not report every read (a crash, or empty or truncated output) is no comparison, and "optik >= 0" would be vacuous.
                mm2_missing = sum(1 for n in ground_truth if n not in mm2_recs)
                o_ok = sum(
                    truth_ok(optik_recs.get(n), g) for n, g in ground_truth.items()
                )
                m_ok = sum(
                    truth_ok(mm2_recs.get(n), g) for n, g in ground_truth.items()
                )
                agree = sum(
                    1
                    for n in ground_truth
                    if optik_recs.get(n)
                    and mm2_recs.get(n)
                    and optik_recs[n]["is_mapped"]
                    and mm2_recs[n]["is_mapped"]
                    and optik_recs[n]["ref"] == mm2_recs[n]["ref"]
                    and abs(optik_recs[n]["pos"] - mm2_recs[n]["pos"]) <= 1
                )
                o_nm = sum(
                    1
                    for r in optik_recs.values()
                    if r["is_mapped"] and r["nm"] not in (0, None)
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
            kind = "short" if args.read_len < 1000 else "long"
            wres = {}
            for cfg, opts in (("defaults", []), ("read_opts", WFMASH_READ_OPTS[kind])):
                raw = os.path.join(tmpdir, f"wfmash_{cfg}.raw.sam")
                one = os.path.join(tmpdir, f"wfmash_{cfg}.sam")
                t0 = time.monotonic()
                with open(raw, "w") as fh:
                    pw = subprocess.run(
                        [
                            args.wfmash,
                            ref_clean,
                            reads_fa,
                            "-a",
                            "-t",
                            str(args.threads),
                            *opts,
                        ],
                        stdout=fh,
                        stderr=subprocess.DEVNULL,
                        timeout=7200,
                        check=False,
                    )
                dt = round(time.monotonic() - t0, 3)
                speed[f"wfmash_{cfg}"] = dt
                if pw.returncode != 0:
                    wres[cfg] = {"error": f"exit {pw.returncode}", "seconds": dt}
                    continue
                one_primary_per_read(raw, one)
                with open(one) as fh:
                    wrecs = parse_sam(fh.read())
                wres[cfg] = {
                    "correct": sum(
                        truth_ok(wrecs.get(n), g) for n, g in ground_truth.items()
                    ),
                    "seconds": dt,
                    "cmd": " ".join(["wfmash", "REF", "READS", "-a", *opts]),
                }
            o_ok = sum(truth_ok(optik_recs.get(n), g) for n, g in ground_truth.items())
            best = max((v.get("correct", 0) for v in wres.values()), default=0)
            # A baseline that crashed or failed is no comparison, so it is reported as a failure rather than scored as 0 correct.
            broken = {c: v["error"] for c, v in wres.items() if "error" in v}
            c = {
                "optik_correct": o_ok,
                "reads": len(ground_truth),
                "wfmash": wres,
                "detail": "optik correct >= wavefront mapper (wfmash) correct, same truth scoring",
            }
            if broken:
                c.update(
                    pass_=False,
                    error=f"wfmash baseline failed: {broken} (a negative "
                    "exit is a signal, e.g. -11 = SIGSEGV)",
                )
            elif best == 0:
                # wfmash ran but placed no read in any configuration, since it is not built for reads this short.
                # "optik >= 0" proves nothing, so this is reported rather than counted as a pass.
                c.update(info="vacuous: wfmash placed 0 reads; not a comparison")
            else:
                c.update(pass_=o_ok >= best)
            if "pass_" in c:
                c["pass"] = c.pop("pass_")
            checks["wfmash_truth_parity"] = c
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
            with open(paf_path, "w") as f:
                f.write("\n".join(paf_lines) + ("\n" if paf_lines else ""))
            wg_report, wg_err = run_wgatools_paf_checks(paf_path, args.wgatools, tmpdir)
            if wg_err:
                checks["wgatools_paf_roundtrip"] = {"pass": False, "error": wg_err}
            else:
                # pafcov must cover exactly the union of the PAF target intervals, computed here independently.
                iv = {}
                for text_line in paf_lines:
                    t = text_line.split("\t")
                    iv.setdefault(t[5], []).append((int(t[7]), int(t[8])))
                union = 0
                for spans in iv.values():
                    end = -1
                    for a, b in sorted(spans):
                        a = max(a, end)
                        if b > a:
                            union += b - a
                            end = b
                n_prim = sum(1 for r in optik_recs.values() if r["is_mapped"])
                ok = (
                    not dropped
                    and len(paf_lines) == n_prim
                    and wg_report["validate_total_records"] == len(paf_lines)
                    and wg_report["validate_query_invalid"] == 0
                    and wg_report["validate_target_invalid"] == 0
                    and wg_report["pafcov_covered_ref_bases"] == union
                )
                checks["wgatools_paf_roundtrip"] = {
                    "pass": ok,
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
        other_t = 1 if args.threads != 1 else 2
        sam_text_2, err = run_optik_map(
            ref_clean, reads_fa, k, optik, tmpdir, threads=other_t
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

        # Check 6, mapping rate sanity
        mapped = sum(1 for r in optik_recs.values() if r["is_mapped"])
        total = len(optik_recs)
        rate = mapped / total * 100 if total else 0
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
        "-n", type=int, default=10000, help="Number of reads to generate"
    )
    parser.add_argument("--read-len", type=int, default=150, help="Read length")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed")
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
        c.get("pass", False)
        for c in checks.values()
        if isinstance(c, dict) and "pass" in c
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
