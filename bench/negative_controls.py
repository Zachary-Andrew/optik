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

from sam_spec import SAM_FLAG_UNMAPPED, SAM_FLAGS_NOT_PRIMARY, SAM_MIN_COLUMNS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OPTIK = os.environ.get("OPTIK_BIN", os.path.join(ROOT, "target/release/optik"))
T = os.path.join(ROOT, "temp")
TOOLS = {
    "minimap2": os.environ.get("MINIMAP2_BIN", f"{T}/minimap2/minimap2"),
    "sshash": os.environ.get("SSHASH_BIN", f"{T}/sshash/build/sshash"),
    "kmc": os.environ.get("KMC_BIN", f"{T}/KMC/bin/kmc"),
    "cuttlefish": os.environ.get(
        "CUTTLEFISH_BIN", f"{T}/cuttlefish/target/release/cuttlefish"
    ),
    "wgatools": os.environ.get("WGATOOLS_BIN", f"{T}/wgatools/target/release/wgatools"),
    "samtools": shutil.which("samtools") or "",
}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


# Output corruptions
# Each takes the list of lines of a text output, or the raw bytes of a binary one, and returns the corrupted version.
# `every` picks 1 record in 100, so a check must notice a small fraction rather than a wholesale breakage.
def sam_records(lines, pred=lambda f: not f & SAM_FLAGS_NOT_PRIMARY, every=100):
    """Yield (index, fields) of every `every`-th SAM record matching pred."""
    n = 0
    for i, line in enumerate(lines):
        if not line or line.startswith("@"):
            continue
        t = line.split("\t")
        if len(t) >= SAM_MIN_COLUMNS and pred(int(t[1])):
            n += 1
            if n % every == 0:
                yield i, t


def m_strand_flip(lines):
    """Flip the strand bit of one primary in 100."""
    for i, t in sam_records(lines):
        t[1] = str(int(t[1]) ^ 16)
        lines[i] = "\t".join(t)


def m_drop_records(lines):
    """Delete every 100th SAM record, primary or not."""
    drop = {i for i, _ in sam_records(lines, pred=lambda f: True)}
    lines[:] = [line for i, line in enumerate(lines) if i not in drop]


def m_dup_primary(lines):
    """Duplicate one primary in 100, so those reads have two primaries."""
    for _i, t in list(sam_records(lines)):
        lines.append("\t".join(t))


def m_pos_shift(lines):
    """Move POS 3 bases right on 1% of primaries."""
    for i, t in sam_records(lines):
        t[3] = str(int(t[3]) + 3)
        lines[i] = "\t".join(t)


def m_nm_wrong(lines):
    """Overstate NM by one on 1% of primaries."""
    for i, t in sam_records(lines):
        lines[i] = "\t".join(
            re.sub(r"^NM:i:(\d+)$", lambda m: f"NM:i:{int(m.group(1)) + 1}", x)
            for x in t
        )


def m_nm_drop(lines):
    """Remove the NM tag from every record."""
    for i, line in enumerate(lines):
        if line and not line.startswith("@"):
            lines[i] = "\t".join(
                x for x in line.split("\t") if not x.startswith("NM:i:")
            )


def m_seq_star(lines):
    """Blank SEQ to '*' on 1% of primaries."""
    for i, t in sam_records(lines):
        t[9] = "*"
        lines[i] = "\t".join(t)


def m_seq_base(lines):
    """Swap the middle SEQ base of one primary in 100."""
    for i, t in sam_records(lines):
        s = t[9]
        j = len(s) // 2
        t[9] = s[:j] + ("A" if s[j] != "A" else "C") + s[j + 1 :]
        lines[i] = "\t".join(t)


def m_cigar_short(lines):
    """Trim one base off the last CIGAR operation of 1% of primaries."""
    for i, t in sam_records(lines):
        ops = re.findall(r"(\d+)([MIDNSHP=X])", t[5])
        n, o = ops[-1]
        ops[-1] = (str(int(n) - 1), o)
        t[5] = "".join(a + b for a, b in ops)
        lines[i] = "\t".join(t)


def m_mapq_zero(lines):
    """Set MAPQ to 0 on every mapped record."""
    for i, line in enumerate(lines):
        t = line.split("\t")
        if (
            line
            and not line.startswith("@")
            and len(t) >= SAM_MIN_COLUMNS
            and not int(t[1]) & SAM_FLAG_UNMAPPED
        ):
            t[4] = "0"
            lines[i] = "\t".join(t)


def m_misplace_confident(lines):
    """Move 3% of primaries 5 kb away with MAPQ 60, making them confidently wrong."""
    for i, t in sam_records(lines, every=33):
        t[3] = str(int(t[3]) + 5000)
        t[4] = "60"
        lines[i] = "\t".join(t)


def m_shuffle_if_threads(lines, args):
    """Make the record order depend on the thread count, i.e. nondeterministic output.

    Every count above 1 shuffles with its own seed, so two runs at different counts disagree whichever counts the validator picks."""
    t = args[args.index("-t") + 1] if "-t" in args else "0"
    if t != "1":
        head = [line for line in lines if line.startswith("@")]
        body = [line for line in lines if line and not line.startswith("@")]
        random.Random(int(t)).shuffle(body)
        lines[:] = head + body + [""]


def m_count_off(lines):
    """Report one positive k-mer too many in the query report."""
    lines[:] = [
        re.sub(
            r"^num_positive_kmers = (\d+)",
            lambda m: f"num_positive_kmers = {int(m.group(1)) + 1}",
            line,
        )
        for line in lines
    ]


def m_ids_flip(lines):
    """Report one present window in 1000 as absent in --ids, and vice versa."""
    n = 0
    for i, line in enumerate(lines):
        ids = line.split()
        for j, x in enumerate(ids):
            n += 1
            if n % 1000 == 0:
                ids[j] = "-1" if x != "-1" else "0"
        lines[i] = " ".join(ids)


def gfa_segments(lines):
    """Indexes of the S lines of a GFA."""
    return [i for i, line in enumerate(lines) if line.startswith("S\t")]


def m_split_unitig(lines):
    """Split the first unitig of 80 bp or more into two overlapping segments."""
    for i in gfa_segments(lines):
        t = lines[i].split("\t")
        if len(t[2]) >= 80:
            s = t[2]
            lines[i] = "\t".join(["S", t[1], s[:50]])
            lines.insert(i + 1, "\t".join(["S", t[1] + "x", s[20:]]))
            return
    raise SystemExit("split_unitig: no unitig >= 80 bp")


def m_dup_segment_name(lines):
    """Give the second segment the first segment's name."""
    seg = gfa_segments(lines)
    a, b = lines[seg[0]].split("\t"), lines[seg[1]].split("\t")
    b[1] = a[1]
    lines[seg[1]] = "\t".join(b)


def m_kc_off(lines):
    """Increase the KC tag of the middle segment by one."""
    i = gfa_segments(lines)[len(gfa_segments(lines)) // 2]
    lines[i] = re.sub(r"KC:i:(\d+)", lambda m: f"KC:i:{int(m.group(1)) + 1}", lines[i])


def m_drop_link(lines):
    """Delete the first L line."""
    for i, line in enumerate(lines):
        if line.startswith("L\t"):
            del lines[i]
            return
    raise SystemExit("drop_link: no L line")


def m_extra_base(lines):
    """Append one base to the first segment."""
    i = gfa_segments(lines)[0]
    t = lines[i].split("\t")
    t[2] += "A"
    lines[i] = "\t".join(t)


def m_concat_unitigs(lines):
    """Join two unrelated unitigs into one segment, which makes a path through a non-edge."""
    seg = gfa_segments(lines)
    # They are far apart in the output, which is sorted by head k-mer, so they are not graph neighbours.
    i, j = seg[len(seg) // 3], seg[2 * len(seg) // 3]
    a, b = lines[i].split("\t"), lines[j].split("\t")
    a[2] = a[2] + b[2]
    lines[i] = "\t".join(a)
    del lines[j]


def m_gfa_no_header(lines):
    """Remove the GFA header line."""
    lines[:] = [line for line in lines if not line.startswith("H\t")]


def m_kdbg_abundance(data):
    """Increase one abundance in the -s node table by one.

    The OKKDBG layout is the 6-byte magic and a k byte, 7 bytes, then the u64 count n, 15 bytes so far, then n u64 keys and n u32 abundances, of which the middle one is changed."""
    n = struct.unpack_from("<Q", data, 7)[0]
    off = 15 + 8 * n + 4 * (n // 2)
    v = struct.unpack_from("<I", data, off)[0]
    return data[:off] + struct.pack("<I", v + 1) + data[off + 4 :]


def m_rc_blind_input(query_fa, k):
    """Input-side mutant emulating a lookup that finds only the canonical orientation.

    Single-k-mer records written non-canonically are replaced by an absent k-mer before the real optik sees them."""
    out = query_fa + ".rcblind"
    with open(query_fa) as fh, open(out, "w") as o:
        for line in fh:
            s = line.strip()
            if not line.startswith(">") and len(s) == k and s > s.translate(COMP)[::-1]:
                line = "A" * k + "\n"
            o.write(line)
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
    "split_unitig": ("dbruijn", "o", m_split_unitig),
    "dup_segment_name": ("dbruijn", "o", m_dup_segment_name),
    "kc_off": ("dbruijn", "o", m_kc_off),
    "drop_link": ("dbruijn", "o", m_drop_link),
    "extra_base": ("dbruijn", "o", m_extra_base),
    "kdbg_abundance": ("dbruijn", "s_bin", m_kdbg_abundance),
    "concat_unitigs": ("dbruijn", "o", m_concat_unitigs),
    "gfa_no_header": ("dbruijn", "o", m_gfa_no_header),
}


def wrap_main(name, argv):
    """Run as the wrapper, i.e. the real optik plus one corruption."""
    sub, where, fn = MUTANTS[name]
    if not argv or argv[0] != sub:
        os.execv(OPTIK, [OPTIK] + argv)

    def after(flag):
        """The argument following `flag` in argv, or None."""
        return argv[argv.index(flag) + 1] if flag in argv else None

    if where == "in_q":
        k = int(after("-k") or 31)
        argv = list(argv)
        argv[argv.index("-q") + 1] = fn(after("-q"), k)
        os.execv(OPTIK, [OPTIK] + argv)
    if where == "stdout":
        p = subprocess.run([OPTIK] + argv, stdout=subprocess.PIPE, check=False)
        lines = p.stdout.decode().split("\n")
        fn(lines)
        sys.stdout.write("\n".join(lines))
        return p.returncode
    p = subprocess.run([OPTIK] + argv, check=False)
    if p.returncode != 0:
        return p.returncode
    target = after({"o": "-o", "ids": "--ids", "s_bin": "-s"}[where])
    if target is None:
        return 0
    if where == "s_bin":
        data = Path(target).read_bytes()
        Path(target).write_bytes(fn(data))
        return 0
    lines = Path(target).read_text().split("\n")
    if fn.__code__.co_argcount == 2:
        fn(lines, argv)
    else:
        fn(lines)
    Path(target).write_text("\n".join(lines))
    return 0


# Validator runs
def fake_tool(td, name, script):
    """Write an executable shell script standing in for a broken baseline tool."""
    path = os.path.join(td, f"fake_{name}")
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\n" + script + "\n")
    os.chmod(path, 0o755)
    return path


def failed_checks(report):
    """Names of the checks a validator report marks as failed."""
    checks = report.get("checks", report)
    return {
        n for n, c in checks.items() if isinstance(c, dict) and c.get("pass") is False
    }


def run_validator(cmd, td):
    """Run a validator and return the set of failed checks and all check names."""
    p = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=td,
        timeout=7200,
        check=False,
    )
    try:
        rep = json.loads(p.stdout[p.stdout.index("{") :])
    except ValueError:
        return {
            "<no report: " + (p.stderr.strip().splitlines() or ["?"])[-1][:150] + ">"
        }, set()
    checks = rep.get("checks", rep)
    return failed_checks(rep), {
        n for n, c in checks.items() if isinstance(c, dict) and "pass" in c
    }


def main():
    """Run every validator clean, then every control, and report which checks failed as required."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("ref")
    ap.add_argument("reads")
    ap.add_argument("--only", default=None)
    ap.add_argument("--wrap", default=None, help=argparse.SUPPRESS)
    args, _rest = ap.parse_known_args() if "--wrap" not in sys.argv else (None, None)
    if "--wrap" in sys.argv:
        i = sys.argv.index("--wrap")
        return wrap_main(sys.argv[i + 1], sys.argv[i + 2 :])

    td = tempfile.mkdtemp(prefix="optik_negctl_")
    py = sys.executable
    # A 600 kb reference slice and 3000 reads keep each run short.
    refslice = os.path.join(td, "ref.fa")
    with (gzip.open if args.ref.endswith(".gz") else open)(args.ref, "rt") as fh:
        seq, name, recs = [], None, []
        for line in fh:
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
    with open(refslice, "w") as o:
        total = 0
        for n, s in recs:
            if total >= 600_000:
                break
            o.write(f">{n}\n{s}\n")
            total += len(s)
    reads = os.path.join(td, "reads.fq")
    with (
        (gzip.open if args.reads.endswith(".gz") else open)(args.reads, "rt") as fh,
        open(reads, "w") as o,
    ):
        for i, line in enumerate(fh):
            if i >= 4 * 3000:
                break
            o.write(line)

    def wrapper(name):
        """Write the wrapper script that runs optik with mutant `name`, returning its path."""
        path = os.path.join(td, f"optik_{name}")
        with open(path, "w") as fh:
            fh.write(
                f"#!/bin/sh\nexec {shlex.quote(py)} {shlex.quote(os.path.abspath(__file__))} "
                f'--wrap {name} "$@"\n'
            )
        os.chmod(path, 0o755)
        return path

    vm_ = [
        py,
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
        py,
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
        py,
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
        py,
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
        py,
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

    def with_optik(cmd, optik, wd):
        """cmd pointed at the optik binary, plus a work directory for validate_official."""
        extra = ["--workdir", wd] if "validate_official.py" in cmd[1] else []
        return cmd + ["--optik", optik] + extra

    # name -> validator, wrapper or None, required failing checks.
    # The requirement is a set, all of which must fail, or a pair of an all-of set and an any-of set.
    fake_mm2 = fake_tool(td, "minimap2", "echo '@HD\tVN:1.6'")
    fake_wga = fake_tool(td, "wgatools", "exit 0")
    fake_wfm = fake_tool(td, "wfmash", "kill -SEGV $$")
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
            + vm_[vm_.index("--minimap2") + 2 :],
            None,
            {"minimap2_truth_parity"},
        ),
        "baseline_wgatools_silent": (
            vm_[: vm_.index("--wgatools") + 1]
            + [fake_wga]
            + vm_[vm_.index("--wgatools") + 2 :],
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
            + vo_map[vo_map.index("--minimap2") + 2 :],
            None,
            {"err_reads_sensitivity_vs_minimap2", "long_reads_sensitivity_vs_minimap2"},
        ),
    }
    only = set(args.only.split(",")) if args.only else None
    bad = 0
    real_tools = {v for v in TOOLS.values() if v}

    def absent_tools(cmd):
        """External tools a validator command needs that are absent."""
        return [a for a in cmd if a in real_tools and not os.path.isfile(a)] + [
            a for a in cmd if a == "" and cmd[cmd.index(a) - 1] == "--samtools"
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
    ):
        if only and not any(controls[c][0] is cmd for c in only):
            continue
        if absent_tools(cmd):
            skipped.append(f"{label} (missing {absent_tools(cmd)})")
            print(f"  [SKIP] {label}: missing {absent_tools(cmd)}")
            continue
        failed, names = run_validator(
            with_optik(cmd, OPTIK, os.path.join(td, "w_" + label[:3] + label[-3:])), td
        )
        clean_names[id(cmd)] = names
        ok = not failed and names
        bad += not ok
        print(
            f"  [{'PASS' if ok else 'FAIL'}] {label}: {len(names)} checks"
            + ("" if ok else f", failed {sorted(failed)}")
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
            with_optik(cmd, optik, os.path.join(td, "w_" + name)), td
        )
        all_of, any_of = want if isinstance(want, tuple) else (want, set())
        missing = (all_of - failed) | (
            any_of if any_of and not any_of & failed else set()
        )
        want = all_of | (any_of & failed)
        ok = not missing
        bad += not ok
        extra = sorted(failed - want)
        print(
            f"  [{'PASS' if ok else 'FAIL'}] {name}: "
            + (f"failed {sorted(want)}" if ok else f"did NOT fail {sorted(missing)}")
            + (f" (also: {extra[:4]})" if extra else "")
        )
    shutil.rmtree(td, ignore_errors=True)
    print(
        f"negative controls: {'ALL PASS' if not bad else f'{bad} FAILURE(S)'}"
        + (f", {len(skipped)} SKIPPED (missing tools)" if skipped else "")
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
