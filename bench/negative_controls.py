#!/usr/bin/env python3
"""negative_controls.py - shows that every validation check can fail.

A check that cannot fail is not a check.
Each control below wraps the real optik binary and corrupts exactly one aspect of its real output (or, for a broken comparison tool, replaces the baseline with a broken fake).
It then runs the validator responsible for that aspect and requires the named checks to fail.

Each validator configuration first runs on the unmodified binary and must pass every check, since otherwise a control could "fail as required" for an unrelated reason.
The controls therefore show both directions.
  clean optik            -> all checks PASS, so no false negative here
  optik with defect X    -> the checks for X FAIL, so no false positive for X

Usage: bench/negative_controls.py <reference.fa> <reads.fq[.gz]> [--only NAME,...] It exits 0 when every clean run passes and every control fails as required.
"""

import argparse
import gzip
import json
import os
import random
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

from diagnostics import DETAIL_CHARS
from fastq import FASTQ_LINES_PER_RECORD
from sam_spec import (
    KDBG_ABUNDANCE_BYTES,
    KDBG_COUNT_OFFSET,
    KDBG_HEADER_BYTES,
    KDBG_KEY_BYTES,
    SAM_MIN_COLUMNS,
    GfaSegmentCol,
    SamCol,
    SamFlag,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OPTIK = os.environ.get("OPTIK_BIN", os.path.join(ROOT, "target/release/optik"))
TEMP_DIR = os.path.join(ROOT, "temp")
TOOLS = {
    "minimap2": os.environ.get("MINIMAP2_BIN", f"{TEMP_DIR}/minimap2/minimap2"),
    "sshash": os.environ.get("SSHASH_BIN", f"{TEMP_DIR}/sshash/build/sshash"),
    "kmc": os.environ.get("KMC_BIN", f"{TEMP_DIR}/KMC/bin/kmc"),
    "cuttlefish": os.environ.get(
        "CUTTLEFISH_BIN", f"{TEMP_DIR}/cuttlefish/target/release/cuttlefish"
    ),
    "wgatools": os.environ.get("WGATOOLS_BIN", f"{TEMP_DIR}/wgatools/target/release/wgatools"),
    "samtools": shutil.which("samtools") or "",
}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# Corruption sizes
CORRUPT_EVERY = 100  # most SAM mutants touch one matching record in this many
POS_SHIFT = 3  # pos_shift moves POS this many bases right
MISPLACE_EVERY = 33  # misplace_confident moves one primary in this many (3%) ...
MISPLACE_SHIFT = 5000  # ... this many bases away
IDS_FLIP_EVERY = 1000  # ids_flip flips one --ids entry in this many
SPLIT_MIN_LEN = 80  # split_unitig splits the first unitig at least this long ...
SPLIT_HEAD_LEN = 50  # ... into its first SPLIT_HEAD_LEN bases ...
SPLIT_TAIL_START = 20  # ... and the overlapping rest from this offset
CONCAT_SPACING_DIV = 3  # concat_unitigs joins the segments a third and two thirds into the output
# Validator runs
DEFAULT_K = 31  # k of a wrapped query when -k is not given
EXEC_MODE = 0o755  # permissions of the wrapper and fake-tool scripts
VALIDATOR_TIMEOUT_S = 7200
REF_SLICE_BASES = 600_000  # the reference slice keeps records until it holds this many bases ...
NUM_READS = 3000  # ... and the read set this many FASTQ records
EXTRA_FAILS_SHOWN = 4  # unexpected extra failed checks listed per control
LABEL_ABBREV_CHARS = 3  # a clean run's work directory is named from this many leading and trailing characters of its label


# Output corruptions
# Each takes the list of lines of a text output, or the raw bytes of a binary one, and returns the corrupted version.
# `every` picks 1 record in 100, so a check must notice a small fraction rather than a wholesale breakage.
# A mutator that also takes the argument list has this many parameters, (lines, argv).
ARGV_MUTATOR_ARGS = 2
# A value flag spans two arguments, the flag and its value.
FLAG_AND_VALUE = 2


def sam_records(lines, pred=lambda flag: not flag & SamFlag.NOT_PRIMARY, every=CORRUPT_EVERY):
    """Yield (index, fields) of every `every`-th SAM record matching pred."""
    count = 0
    for i, line in enumerate(lines):
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) >= SAM_MIN_COLUMNS and pred(int(fields[SamCol.FLAG])):
            count += 1
            if count % every == 0:
                yield i, fields


def m_strand_flip(lines):
    """Flip the strand bit of one primary in 100."""
    for line_idx, fields in sam_records(lines):
        fields[SamCol.FLAG] = str(int(fields[SamCol.FLAG]) ^ SamFlag.REVERSE)
        lines[line_idx] = "\t".join(fields)


def m_drop_records(lines):
    """Delete every 100th SAM record, primary or not."""
    drop = {line_idx for line_idx, _ in sam_records(lines, pred=lambda flag: True)}
    lines[:] = [line for i, line in enumerate(lines) if i not in drop]


def m_dup_primary(lines):
    """Duplicate one primary in 100, so those reads have two primaries."""
    for _line_idx, fields in list(sam_records(lines)):
        lines.append("\t".join(fields))


def m_pos_shift(lines):
    """Move POS 3 bases right on 1% of primaries."""
    for line_idx, fields in sam_records(lines):
        fields[SamCol.POS] = str(int(fields[SamCol.POS]) + POS_SHIFT)
        lines[line_idx] = "\t".join(fields)


def m_nm_wrong(lines):
    """Overstate NM by one on 1% of primaries."""
    for line_idx, fields in sam_records(lines):
        lines[line_idx] = "\t".join(
            re.sub(r"^NM:i:(\d+)$", lambda match: f"NM:i:{int(match.group(1)) + 1}", field)
            for field in fields
        )


def m_nm_drop(lines):
    """Remove the NM tag from every record."""
    for i, line in enumerate(lines):
        if line and not line.startswith("@"):
            lines[i] = "\t".join(
                field for field in line.split("\t") if not field.startswith("NM:i:")
            )


def m_seq_star(lines):
    """Blank SEQ to '*' on 1% of primaries."""
    for line_idx, fields in sam_records(lines):
        fields[SamCol.SEQ] = "*"
        lines[line_idx] = "\t".join(fields)


def m_seq_base(lines):
    """Swap the middle SEQ base of one primary in 100."""
    for line_idx, fields in sam_records(lines):
        seq = fields[SamCol.SEQ]
        mid = len(seq) // 2
        fields[SamCol.SEQ] = seq[:mid] + ("A" if seq[mid] != "A" else "C") + seq[mid + 1 :]
        lines[line_idx] = "\t".join(fields)


def m_cigar_short(lines):
    """Trim one base off the last CIGAR operation of 1% of primaries."""
    for line_idx, fields in sam_records(lines):
        ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
        last_len, last_op = ops[-1]
        ops[-1] = (str(int(last_len) - 1), last_op)
        fields[SamCol.CIGAR] = "".join(op_len + cigar_op for op_len, cigar_op in ops)
        lines[line_idx] = "\t".join(fields)


def m_mapq_zero(lines):
    """Set MAPQ to 0 on every mapped record."""
    for i, line in enumerate(lines):
        fields = line.split("\t")
        if (
            line
            and not line.startswith("@")
            and len(fields) >= SAM_MIN_COLUMNS
            and not int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED
        ):
            fields[SamCol.MAPQ] = "0"
            lines[i] = "\t".join(fields)


def m_misplace_confident(lines):
    """Move 3% of primaries 5 kb away with MAPQ 60, making them confidently wrong."""
    for line_idx, fields in sam_records(lines, every=MISPLACE_EVERY):
        fields[SamCol.POS] = str(int(fields[SamCol.POS]) + MISPLACE_SHIFT)
        fields[SamCol.MAPQ] = "60"
        lines[line_idx] = "\t".join(fields)


def m_shuffle_if_threads(lines, args):
    """Make the record order depend on the thread count, i.e. nondeterministic output.

    Every count above 1 shuffles with its own seed, so two runs at different counts disagree whichever counts the validator picks."""
    threads = args[args.index("-t") + 1] if "-t" in args else "0"
    if threads != "1":
        head = [line for line in lines if line.startswith("@")]
        body = [line for line in lines if line and not line.startswith("@")]
        random.Random(int(threads)).shuffle(body)
        lines[:] = head + body + [""]


def m_count_off(lines):
    """Report one positive k-mer too many in the query report."""
    lines[:] = [
        re.sub(
            r"^num_positive_kmers = (\d+)",
            lambda match: f"num_positive_kmers = {int(match.group(1)) + 1}",
            line,
        )
        for line in lines
    ]


def m_ids_flip(lines):
    """Report one present window in 1000 as absent in --ids, and vice versa."""
    count = 0
    for i, line in enumerate(lines):
        ids = line.split()
        for j, id_value in enumerate(ids):
            count += 1
            if count % IDS_FLIP_EVERY == 0:
                ids[j] = "-1" if id_value != "-1" else "0"
        lines[i] = " ".join(ids)


def m_ids_flip_first(lines):
    """Flip the first window's answer in --ids, for outputs too small for ids_flip's one in IDS_FLIP_EVERY."""
    for i, line in enumerate(lines):
        ids = line.split()
        if ids:
            ids[0] = "-1" if ids[0] != "-1" else "0"
            lines[i] = " ".join(ids)
            return


def gfa_segments(lines):
    """Indexes of the S lines of a GFA."""
    return [i for i, line in enumerate(lines) if line.startswith("S\t")]


def m_split_unitig(lines):
    """Split the first unitig of 80 bp or more into two overlapping segments."""
    for line_idx in gfa_segments(lines):
        fields = lines[line_idx].split("\t")
        if len(fields[GfaSegmentCol.SEQ]) >= SPLIT_MIN_LEN:
            seq = fields[GfaSegmentCol.SEQ]
            lines[line_idx] = "\t".join(["S", fields[GfaSegmentCol.NAME], seq[:SPLIT_HEAD_LEN]])
            lines.insert(line_idx + 1, "\t".join(["S", fields[GfaSegmentCol.NAME] + "x", seq[SPLIT_TAIL_START:]]))
            return
    raise SystemExit("split_unitig: no unitig >= 80 bp")


def m_split_short_unitig(lines, argv):
    """Split the first unitig FASTA record holding two k-mers or more into two records overlapping by k - 1 bases.
    The k-mer set stays the same, so only a check of the unitig partition can notice."""
    kmer_len = int(argv[argv.index("-k") + 1])
    for i, line in enumerate(lines):
        if line and not line.startswith(">") and len(line) > kmer_len:
            lines[i : i + 1] = [line[:kmer_len], ">split", line[1:]]
            return


def m_dup_segment_name(lines):
    """Give the second segment the first segment's name."""
    seg = gfa_segments(lines)
    first_fields, second_fields = lines[seg[0]].split("\t"), lines[seg[1]].split("\t")
    second_fields[GfaSegmentCol.NAME] = first_fields[GfaSegmentCol.NAME]
    lines[seg[1]] = "\t".join(second_fields)


def m_kc_off(lines):
    """Increase the KC tag of the middle segment by one."""
    line_idx = gfa_segments(lines)[len(gfa_segments(lines)) // 2]
    lines[line_idx] = re.sub(r"KC:i:(\d+)", lambda match: f"KC:i:{int(match.group(1)) + 1}", lines[line_idx])


def m_drop_link(lines):
    """Delete the first L line."""
    for i, line in enumerate(lines):
        if line.startswith("L\t"):
            del lines[i]
            return
    raise SystemExit("drop_link: no L line")


def m_drop_gfa_link(lines, argv):
    """Delete the first L line of a GFA output, and leave unitig FASTA alone, as one lost edge bit would."""
    if argv[argv.index("-o") + 1].endswith(".gfa"):
        for i, line in enumerate(lines):
            if line.startswith("L\t"):
                del lines[i]
                return


def m_extra_base(lines):
    """Append one base to the first segment."""
    line_idx = gfa_segments(lines)[0]
    fields = lines[line_idx].split("\t")
    fields[GfaSegmentCol.SEQ] += "A"
    lines[line_idx] = "\t".join(fields)


def m_concat_unitigs(lines):
    """Join two unrelated unitigs into one segment, which makes a path through a non-edge."""
    seg = gfa_segments(lines)
    # They are far apart in the output, which is sorted by head k-mer, so they are not graph neighbours.
    first_idx, second_idx = seg[len(seg) // CONCAT_SPACING_DIV], seg[2 * len(seg) // CONCAT_SPACING_DIV]
    first_fields, second_fields = lines[first_idx].split("\t"), lines[second_idx].split("\t")
    first_fields[GfaSegmentCol.SEQ] = first_fields[GfaSegmentCol.SEQ] + second_fields[GfaSegmentCol.SEQ]
    lines[first_idx] = "\t".join(first_fields)
    del lines[second_idx]


def m_gfa_no_header(lines):
    """Remove the GFA header line."""
    lines[:] = [line for line in lines if not line.startswith("H\t")]


def m_kdbg_abundance(data):
    """Increase one abundance in the -s node table by one.

    The OKKDBG layout is the 6-byte magic and a k byte, 7 bytes, then the u64 count n, 15 bytes so far, then n u64 keys and n u32 abundances, of which the middle one is changed."""
    num_nodes = struct.unpack_from("<Q", data, KDBG_COUNT_OFFSET)[0]
    off = KDBG_HEADER_BYTES + KDBG_KEY_BYTES * num_nodes + KDBG_ABUNDANCE_BYTES * (num_nodes // 2)
    abundance = struct.unpack_from("<I", data, off)[0]
    return data[:off] + struct.pack("<I", abundance + 1) + data[off + KDBG_ABUNDANCE_BYTES :]


def m_rc_blind_input(query_fa, kmer_len):
    """Input-side mutant emulating a lookup that finds only the canonical orientation.

    Single-k-mer records written non-canonically are replaced by an absent k-mer before the real optik reads them."""
    out = query_fa + ".rcblind"
    with open(query_fa) as in_handle, open(out, "w") as out_handle:
        for line in in_handle:
            seq = line.strip()
            if not line.startswith(">") and len(seq) == kmer_len and seq > seq.translate(COMP)[::-1]:
                line = "A" * kmer_len + "\n"
            out_handle.write(line)
    return out


# subcommand -> the output each mutant rewrites
#   ("o", fn)      text file after -o      ("stdout", fn)  captured stdout
#   ("ids", fn)    text file after --ids   ("s_bin", fn)   binary file after -s
#   ("in_q", fn)   rewrites the -q input   fn(lines, args) when it takes args
MUTANTS = {
    "strand_flip": ("map", "o", m_strand_flip),
    "drop_records": ("map", "o", m_drop_records),
    "dup_primary": ("map", "o", m_dup_primary),
    "pos_shift": ("map", "o", m_pos_shift),
    "nm_wrong": ("map", "o", m_nm_wrong),
    "nm_drop": ("map", "o", m_nm_drop),
    "seq_star": ("map", "o", m_seq_star),
    "seq_base": ("map", "o", m_seq_base),
    "cigar_short": ("map", "o", m_cigar_short),
    "mapq_zero": ("map", "o", m_mapq_zero),
    "misplace_confident": ("map", "o", m_misplace_confident),
    "nondeterministic": ("map", "o", m_shuffle_if_threads),
    "count_off": ("query", "stdout", m_count_off),
    "rc_blind": ("query", "in_q", m_rc_blind_input),
    "ids_flip": ("query", "ids", m_ids_flip),
    "ids_flip_first": ("query", "ids", m_ids_flip_first),
    "split_unitig": ("dbruijn", "o", m_split_unitig),
    "split_short_unitig": ("dbruijn", "o", m_split_short_unitig),
    "dup_segment_name": ("dbruijn", "o", m_dup_segment_name),
    "kc_off": ("dbruijn", "o", m_kc_off),
    "drop_link": ("dbruijn", "o", m_drop_link),
    "drop_gfa_link": ("dbruijn", "o", m_drop_gfa_link),
    "extra_base": ("dbruijn", "o", m_extra_base),
    "kdbg_abundance": ("dbruijn", "s_bin", m_kdbg_abundance),
    "concat_unitigs": ("dbruijn", "o", m_concat_unitigs),
    "gfa_no_header": ("dbruijn", "o", m_gfa_no_header),
}


def wrap_main(name, argv):
    """Run as the wrapper, i.e. the real optik plus one corruption."""
    sub, where, mutate = MUTANTS[name]
    if not argv or argv[0] != sub:
        os.execv(OPTIK, [OPTIK] + argv)

    def after(flag):
        """The argument following `flag` in argv, or None."""
        return argv[argv.index(flag) + 1] if flag in argv else None

    if where == "in_q":
        kmer_len = int(after("-k") or DEFAULT_K)
        argv = list(argv)
        argv[argv.index("-q") + 1] = mutate(after("-q"), kmer_len)
        os.execv(OPTIK, [OPTIK] + argv)
    if where == "stdout":
        proc = subprocess.run([OPTIK] + argv, stdout=subprocess.PIPE, check=False)
        lines = proc.stdout.decode().split("\n")
        mutate(lines)
        sys.stdout.write("\n".join(lines))
        return proc.returncode
    proc = subprocess.run([OPTIK] + argv, check=False)
    if proc.returncode != 0:
        return proc.returncode
    target = after({"o": "-o", "ids": "--ids", "s_bin": "-s"}[where])
    if target is None:
        return 0
    if where == "s_bin":
        data = Path(target).read_bytes()
        Path(target).write_bytes(mutate(data))
        return 0
    lines = Path(target).read_text().split("\n")
    if mutate.__code__.co_argcount == ARGV_MUTATOR_ARGS:
        mutate(lines, argv)
    else:
        mutate(lines)
    Path(target).write_text("\n".join(lines))
    return 0


# Validator runs
def fake_tool(tmp_dir, name, script):
    """Write an executable shell script standing in for a broken baseline tool."""
    path = os.path.join(tmp_dir, f"fake_{name}")
    with open(path, "w") as handle:
        handle.write("#!/bin/sh\n" + script + "\n")
    os.chmod(path, EXEC_MODE)
    return path


def failed_checks(report):
    """Names of the checks a validator report marks as failed."""
    checks = report.get("checks", report)
    return {
        name for name, check in checks.items() if isinstance(check, dict) and check.get("pass") is False
    }


def run_validator(cmd, tmp_dir):
    """Run a validator and return the set of failed checks and all check names."""
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=tmp_dir,
        timeout=VALIDATOR_TIMEOUT_S,
        check=False,
    )
    try:
        rep = json.loads(proc.stdout[proc.stdout.index("{") :])
    except ValueError:
        return {
            "<no report: " + (proc.stderr.strip().splitlines() or ["?"])[-1][:DETAIL_CHARS] + ">"
        }, set()
    checks = rep.get("checks", rep)
    return failed_checks(rep), {
        name for name, check in checks.items() if isinstance(check, dict) and "pass" in check
    }


def main():
    """Run every validator clean, then every control, and report which checks failed as required."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("ref")
    parser.add_argument("reads")
    parser.add_argument("--only", default=None)
    parser.add_argument("--wrap", default=None, help=argparse.SUPPRESS)
    args, _rest = parser.parse_known_args() if "--wrap" not in sys.argv else (None, None)
    if "--wrap" in sys.argv:
        wrap_idx = sys.argv.index("--wrap")
        return wrap_main(sys.argv[wrap_idx + 1], sys.argv[wrap_idx + FLAG_AND_VALUE :])

    tmp_dir = tempfile.mkdtemp(prefix="optik_negctl_")
    python = sys.executable
    # A 600 kb reference slice and 3000 reads keep each run short.
    refslice = os.path.join(tmp_dir, "ref.fa")
    with (gzip.open if args.ref.endswith(".gz") else open)(args.ref, "rt") as handle:
        seq, name, recs = [], None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = (
                    line[1:].split()[0] if line[1:].split() else f"s{len(recs)}",
                    [],
                )
            else:
                seq.append(line.strip())
        recs.append((name, "".join(seq)))
    with open(refslice, "w") as out_handle:
        total = 0
        for rec_name, rec_seq in recs:
            if total >= REF_SLICE_BASES:
                break
            out_handle.write(f">{rec_name}\n{rec_seq}\n")
            total += len(rec_seq)
    reads = os.path.join(tmp_dir, "reads.fq")
    with (
        (gzip.open if args.reads.endswith(".gz") else open)(args.reads, "rt") as handle,
        open(reads, "w") as out_handle,
    ):
        for wrap_idx, line in enumerate(handle):
            if wrap_idx >= FASTQ_LINES_PER_RECORD * NUM_READS:
                break
            out_handle.write(line)

    def wrapper(name):
        """Write the wrapper script that runs optik with mutant `name`, returning its path."""
        path = os.path.join(tmp_dir, f"optik_{name}")
        with open(path, "w") as handle:
            handle.write(
                f"#!/bin/sh\nexec {shlex.quote(python)} {shlex.quote(os.path.abspath(__file__))} "
                # The wrapper re-enters this script, which applies the mutation around the real binary.
                f'--wrap {name} "$@"\n'
            )
        os.chmod(path, EXEC_MODE)
        return path

    # Each validator runs small and seeded, so a mutant's failure reproduces quickly and identically.
    vm_ = [
        python,
        f"{HERE}/validate_mapping.py",
        "-r",
        refslice,
        "-k",
        "31",
        "-n",
        "2000",
        "-s",
        "42",
        "--minimap2",
        TOOLS["minimap2"],
        "--wgatools",
        TOOLS["wgatools"],
    ]
    vq_ = [
        python,
        f"{HERE}/validate_query.py",
        "-r",
        refslice,
        "-k",
        "31",
        "-s",
        "42",
        "--max-query-kmers",
        "5000",
        "--sshash",
        TOOLS["sshash"],
    ]
    vd_ = [
        python,
        f"{HERE}/validate_dbg.py",
        "-i",
        reads,
        "-k",
        "31",
        "--min-count",
        "1",
        "-t",
        "4",
        "--cuttlefish",
        TOOLS["cuttlefish"],
    ]
    vo_map = [
        python,
        f"{HERE}/validate_official.py",
        "map",
        "-r",
        refslice,
        "-k",
        "31",
        "-n",
        "1500",
        "--n-long",
        "100",
        "-t",
        "4",
        "--read-len",
        "150",
        "--err-rate",
        "0.02",
        "-s",
        "42",
        "--minimap2",
        TOOLS["minimap2"],
        "--samtools",
        TOOLS["samtools"],
    ]
    vo_q = [
        python,
        f"{HERE}/validate_official.py",
        "query",
        "-r",
        refslice,
        "-k",
        "31",
        "-t",
        "4",
        "--sshash",
        TOOLS["sshash"],
        "--kmc",
        TOOLS["kmc"],
        "--probes",
        "500",
        "-s",
        "42",
    ]

    # The checked-in Lean vectors, so this control needs neither pixi nor a Lean build.
    vl_ = [python, f"{HERE}/validate_lean.py", "--vectors", os.path.join(ROOT, "lean", "vectors"), "--json"]

    def with_optik(cmd, optik, work_dir):
        """cmd pointed at the optik binary, plus a work directory for validate_official."""
        extra = ["--workdir", work_dir] if "validate_official.py" in cmd[1] else []
        return cmd + ["--optik", optik] + extra

    # name -> validator, wrapper or None, required failing checks.
    # The requirement is a set, all of which must fail, or a pair of an all-of set and an any-of set.
    fake_mm2 = fake_tool(tmp_dir, "minimap2", "echo '@HD\tVN:1.6'")
    fake_wga = fake_tool(tmp_dir, "wgatools", "exit 0")
    fake_wfm = fake_tool(tmp_dir, "wfmash", "kill -SEGV $$")
    controls = {
        "strand_flip": (vm_, "strand_flip", {"gt_strand_match"}),
        "drop_records": (vm_, "drop_records", {"one_primary_per_read"}),
        "dup_primary": (vm_, "dup_primary", {"one_primary_per_read"}),
        "pos_shift": (vm_, "pos_shift", {"gt_position_match", "record_audit"}),
        "nm_wrong": (vm_, "nm_wrong", {"record_audit"}),
        "nm_drop": (vm_, "nm_drop", {"record_audit"}),
        "seq_star": (vm_, "seq_star", {"record_audit"}),
        "seq_base": (vm_, "seq_base", {"record_audit"}),
        "cigar_short": (
            vm_,
            "cigar_short",
            {"record_audit", "cigar_consistency", "wgatools_paf_roundtrip"},
        ),
        "nondeterministic": (vm_, "nondeterministic", {"reproducibility"}),
        "mapq_zero": (vo_map, "mapq_zero", {"mapq_informative"}),
        "misplace_confident": (
            vo_map,
            "misplace_confident",
            {"err_reads_mapq_calibration", "sim_exact_position_strand"},
        ),
        "count_off": (vq_, "count_off", {"spss_positive_query"}),
        "rc_blind": (vq_, "rc_blind", {"spss_positive_query", "fm_positive_query"}),
        "ids_flip": (vo_q, "ids_flip", {"per_kmer_membership_spss_simreads"}),
        "split_unitig": (
            vd_,
            "split_unitig",
            {"optik_maximal", "cuttlefish_partition"},
        ),
        "dup_segment_name": (vd_, "dup_segment_name", {"optik_segments_wellformed"}),
        "kc_off": (vd_, "kc_off", {"optik_tags_exact"}),
        "drop_link": (vd_, "drop_link", {"optik_links_exact"}),
        # The appended base may create a k-mer that exists elsewhere, leaving the set unchanged, so any structural exactness check may be the one that fires.
        "extra_base": (
            vd_,
            "extra_base",
            (
                {"optik_tags_exact"},
                {
                    "optik_kmer_set_equal",
                    "optik_no_duplicates",
                    "optik_paths_valid",
                    "optik_interior_compact",
                },
            ),
        ),
        "kdbg_abundance": (vd_, "kdbg_abundance", {"optik_kdbg_exact"}),
        "concat_unitigs": (vd_, "concat_unitigs", ({"optik_paths_valid"}, set())),
        "gfa_no_header": (vd_, "gfa_no_header", {"gfa_compliance"}),
        # The Lean vectors must catch a finer unitig partition, a GFA link lost at -c 2 and flipped presence answers
        "lean_split_unitig": (
            vl_,
            "split_short_unitig",
            {
                "lean_dbg_table",
                "lean_dbg_part",
                "lean_dbg_ef",
                "lean_dbg_solid_table",
                "lean_dbg_solid_part",
                "lean_dbg_solid_ef",
                "lean_dbg_solid_spill",
            },
        ),
        "lean_drop_link": (
            vl_,
            "drop_gfa_link",
            {"lean_dbg_edges_table", "lean_dbg_edges_part", "lean_dbg_edges_ef", "lean_dbg_edges_spill"},
        ),
        "lean_ids_flip": (vl_, "ids_flip_first", {"lean_spss_presence"}),
        # The official map suite's own record, concordance and determinism checks
        "official_nm_wrong": (
            vo_map,
            "nm_wrong",
            {
                "err_reads_nm_matches_cigar",
                "sim_reads_records_valid",
                "long_reads_valid_sam",
            },
        ),
        "official_seq_base": (
            vo_map,
            "seq_base",
            {"err_reads_seq_orientation", "sim_reads_records_valid"},
        ),
        "official_pos_shift": (
            vo_map,
            "pos_shift",
            {"minimap2_concordance", "sim_exact_position_strand"},
        ),
        "official_nondeterministic": (
            vo_map,
            "nondeterministic",
            {"map_determinism_t1_vs_tN"},
        ),
        "baseline_minimap2_empty": (
            vm_[: vm_.index("--minimap2") + 1]
            + [fake_mm2]
            + vm_[vm_.index("--minimap2") + FLAG_AND_VALUE :],
            None,
            {"minimap2_truth_parity"},
        ),
        "baseline_wgatools_silent": (
            vm_[: vm_.index("--wgatools") + 1]
            + [fake_wga]
            + vm_[vm_.index("--wgatools") + FLAG_AND_VALUE :],
            None,
            {"wgatools_paf_roundtrip"},
        ),
        "baseline_wfmash_crash": (
            vm_ + ["--wfmash", fake_wfm],
            None,
            {"wfmash_truth_parity"},
        ),
        "baseline_official_minimap2_empty": (
            vo_map[: vo_map.index("--minimap2") + 1]
            + [fake_mm2]
            + vo_map[vo_map.index("--minimap2") + FLAG_AND_VALUE :],
            None,
            {"err_reads_sensitivity_vs_minimap2", "long_reads_sensitivity_vs_minimap2"},
        ),
    }
    only = set(args.only.split(",")) if args.only else None
    bad = 0
    real_tools = {tool_path for tool_path in TOOLS.values() if tool_path}

    def absent_tools(cmd):
        """External tools a validator command needs that are absent."""
        return [arg for arg in cmd if arg in real_tools and not os.path.isfile(arg)] + [
            arg for arg in cmd if arg == "" and cmd[cmd.index(arg) - 1] == "--samtools"
        ]

    skipped = []

    # Clean runs, where every validator configuration must pass on the real binary
    print("== clean optik: every check must PASS")
    clean_names = {}
    for label, cmd in (
        ("validate_mapping", vm_),
        ("validate_query", vq_),
        ("validate_dbg", vd_),
        ("official map", vo_map),
        ("official query", vo_q),
        ("validate_lean", vl_),
    ):
        if only and not any(controls[control][0] is cmd for control in only):
            continue
        if absent_tools(cmd):
            skipped.append(f"{label} (missing {absent_tools(cmd)})")
            print(f"  [SKIP] {label}: missing {absent_tools(cmd)}")
            continue
        failed, names = run_validator(
            with_optik(cmd, OPTIK, os.path.join(tmp_dir, "w_" + label[:LABEL_ABBREV_CHARS] + label[-LABEL_ABBREV_CHARS:])), tmp_dir
        )
        clean_names[id(cmd)] = names
        passed = not failed and names
        bad += not passed
        print(
            f"  [{'PASS' if passed else 'FAIL'}] {label}: {len(names)} checks"
            + ("" if passed else f", failed {sorted(failed)}")
        )

    print("== corrupted optik / broken baselines: the named checks must FAIL")
    for name, (cmd, mutant, want) in controls.items():
        if only and name not in only:
            continue
        if absent_tools(cmd):
            skipped.append(f"{name} (missing {absent_tools(cmd)})")
            print(f"  [SKIP] {name}: missing {absent_tools(cmd)}")
            continue
        optik = wrapper(mutant) if mutant else OPTIK
        failed, names = run_validator(
            with_optik(cmd, optik, os.path.join(tmp_dir, "w_" + name)), tmp_dir
        )
        all_of, any_of = want if isinstance(want, tuple) else (want, set())
        missing = (all_of - failed) | (
            any_of if any_of and not any_of & failed else set()
        )
        want = all_of | (any_of & failed)
        passed = not missing
        bad += not passed
        extra = sorted(failed - want)
        print(
            f"  [{'PASS' if passed else 'FAIL'}] {name}: "
            + (f"failed {sorted(want)}" if passed else f"did NOT fail {sorted(missing)}")
            + (f" (also: {extra[:EXTRA_FAILS_SHOWN]})" if extra else "")
        )
    shutil.rmtree(tmp_dir, ignore_errors=True)
    print(
        f"negative controls: {'ALL PASS' if not bad else f'{bad} FAILURE(S)'}"
        + (f", {len(skipped)} SKIPPED (missing tools)" if skipped else "")
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
