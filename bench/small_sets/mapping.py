"""Runs the mapping section, which scores each mapper's SAM against the simulated truth.

The .fai writer and truth loader it depends on are defined here too.
"""

import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass

from diagnostics import DETAIL_CHARS
from dna import read_fasta, write_fasta
from external_oracles import external_sam_checks
from published_mapping import MM2_PRESET_KW
from sam_spec import HIGH_MAPQ, SAM_TAG_PREFIX_LEN, SamCol, SamFlag, TruthCol
from small_common import BIN, TEMP_DIR
from units import PERCENT

from small_sets.common import (
    CHIMERA_SEED,
    CHIMERA_TAIL_RANGE,
    HERE,
    MAX_EXAMPLES,
    MAX_RANDOM_MAPPED_PCT,
    MIN_LONG_OVERLAP,
    MIN_ORIGIN_TOL,
    NEG_READ_LEN,
    NEG_READS,
    NEG_SEED,
    ORIGIN_TOL_DIV,
    SYN_CHROM_LEN,
    SYN_CHROMS,
    SYN_SEED,
)


def _ref_span(fields):
    """Reference bases spanned by the CIGAR of SAM fields t."""
    return sum(
        int(op_len) for op_len, cigar_op in re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR]) if cigar_op in "MDN=X"
    )


def sam_check(sam, ref, truth, interval=False, reads=None):
    """Validity and accuracy of a SAM against the reference and the truth.

    A read is placed correctly on the right contig and strand when, for short reads, POS - 1 - leading clip lies within max(20, L/10) of the origin.
    With interval=True (long or chimeric reads, whose clips may be foreign sequence) the aligned reference interval must instead overlap the truth interval by at least half its length.

    `reads` is given for optik's rows and selects strict mode.
    optik always writes SEQ and NM:i and never hard-clips, so a mapped record lacking SEQ or NM, holding H, or whose SEQ is not the read (reverse-complemented for FLAG 16) is invalid.
    The lenient rules below exist only for the documented formats of other tools.

    Two SAM conventions are not output errors for the other tools.
    wfmash writes every mapped segment of a read as a primary record, none flagged supplementary, and omits unmapped reads.
    For those tools both are counted in spec_notes, and the read is judged by its longest primary alignment, the rule of minimap2's paftools.js mapeval (-m 0), which does not depend on output order.
    For optik (strict) they stay invalid, since optik writes one primary per read and every read."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    seqs = {}
    for name, seq in read_fasta(ref):
        seqs[name] = seq.upper()
    prim = {}
    bad = []
    spec = []
    strict = reads is not None
    nm_tot = aln_tot = 0
    with open(sam) as handle:
        for line in handle:
            if line[0] == "@":
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[SamCol.FLAG])
            if not flag & SamFlag.SECONDARY_SUPPLEMENTARY:
                if fields[SamCol.QNAME] in prim:
                    (bad if strict else spec).append(f"{fields[SamCol.QNAME]}: two primaries")
                    if not strict and _ref_span(fields) <= _ref_span(prim[fields[SamCol.QNAME]]):
                        continue  # judge the read by its longest primary
                prim[fields[SamCol.QNAME]] = fields
            if flag & SamFlag.UNMAPPED:
                continue
            ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
            qlen = sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MIS=X")
            if reads is not None:
                read = reads.get(fields[SamCol.QNAME])
                want = (
                    None
                    if read is None
                    else (read.translate(comp)[::-1] if flag & SamFlag.REVERSE else read)
                )
                if (
                    fields[SamCol.SEQ] == "*"
                    or any(cigar_op == "H" for _, cigar_op in ops)
                    or not any(tag_field.startswith("NM:i:") for tag_field in fields[SamCol.TAGS:])
                    or want is None
                    or fields[SamCol.SEQ].upper() != want.upper()
                ):
                    bad.append(f"{fields[SamCol.QNAME]}: missing SEQ/NM, hard clip, or SEQ != read")
                    continue
            if fields[SamCol.SEQ] != "*" and qlen != len(fields[SamCol.SEQ]):
                bad.append(f"{fields[SamCol.QNAME]}: CIGAR consumes {qlen} != SEQ {len(fields[SamCol.SEQ])}")
                continue
            if fields[SamCol.SEQ] == "*" or any(cigar_op == "H" for _, cigar_op in ops):
                continue  # hard-clipped or without SEQ, so NM cannot be recomputed
            ref_seq = seqs[fields[SamCol.RNAME]]
            ref_pos, read_pos, nm_count, aligned = int(fields[SamCol.POS]) - 1, 0, 0, 0
            for op_len, cigar_op in ops:
                op_len = int(op_len)
                if cigar_op in "M=X":
                    seg = fields[SamCol.SEQ][read_pos : read_pos + op_len].upper()
                    nm_count += sum(
                        1 for i in range(op_len) if seg[i] != ref_seq[ref_pos + i] or seg[i] == "N"
                    )
                    read_pos += op_len
                    ref_pos += op_len
                    aligned += op_len
                elif cigar_op == "I":
                    nm_count += op_len
                    read_pos += op_len
                    aligned += op_len
                elif cigar_op == "D":
                    nm_count += op_len
                    ref_pos += op_len
                elif cigar_op == "S":
                    read_pos += op_len
            tag = [tag_field for tag_field in fields[SamCol.TAGS:] if tag_field.startswith("NM:i:")]
            if tag and int(tag[0][SAM_TAG_PREFIX_LEN:]) != nm_count:
                bad.append(f"{fields[SamCol.QNAME]}: NM:i {tag[0][SAM_TAG_PREFIX_LEN:]} != recomputed {nm_count}")
            nm_tot += nm_count
            aln_tot += aligned
    missing = [name for name in truth if name not in prim]
    if missing:
        (bad if strict else spec).append(
            f"{len(missing)} reads without a primary record"
        )
    ok_reads = mapped = wrong = hq_wrong = 0
    for name, truth_rec in truth.items():
        fields = prim.get(name)
        if not fields or int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
            continue
        mapped += 1
        if truth_rec is None:
            continue
        rname, true_pos, strand, read_len = truth_rec
        ops = re.findall(r"(\d+)([MIDNSHP=X])", fields[SamCol.CIGAR])
        same = fields[SamCol.RNAME] == rname and ((int(fields[SamCol.FLAG]) & SamFlag.REVERSE) > 0) == (strand == "-")
        if interval:
            span_start = int(fields[SamCol.POS]) - 1
            span_end = span_start + sum(int(op_len) for op_len, cigar_op in ops if cigar_op in "MD=X")
            good = (
                same and span_end > span_start and min(span_end, true_pos + read_len) - max(span_start, true_pos) >= MIN_LONG_OVERLAP * (span_end - span_start)
            )
        else:
            lead = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
            good = same and abs(int(fields[SamCol.POS]) - 1 - lead - true_pos) <= max(MIN_ORIGIN_TOL, read_len // ORIGIN_TOL_DIV)
        ok_reads += good
        if not good:
            wrong += 1
            hq_wrong += int(fields[SamCol.MAPQ]) >= HIGH_MAPQ
    total = len(truth)
    return {
        "valid": not bad,
        "invalid_examples": bad[:MAX_EXAMPLES],
        "invalid_count": len(bad),
        "spec_notes": spec[:MAX_EXAMPLES],
        "spec_count": len(spec),
        "mapped_pct": round(PERCENT * mapped / total, 3),
        "correct_pct": round(PERCENT * ok_reads / total, 3),
        "wrong": wrong,
        "wrong_mapq30": hq_wrong,
        "nm_per_aligned_base": round(nm_tot / aln_tot, 4) if aln_tot else None,
    }


@dataclass
class FaiEntry:
    """One row of a samtools .fai index."""

    name: str
    length: int = 0
    offset: int = 0
    """Byte offset of the sequence's first base."""
    line_bases: int = 0
    line_bytes: int = 0
    """Bytes of a full sequence line, the newline included."""

    def row(self):
        return f"{self.name}\t{self.length}\t{self.offset}\t{self.line_bases}\t{self.line_bytes}\n"


def write_fai(path):
    """Write a samtools-compatible .fai with name, length, offset, bases per line and bytes per line."""
    out = []
    name = None
    with open(path, "rb") as handle:
        data = handle.read()
    pos = 0
    data_len = len(data)
    while pos < data_len:
        line_end = data.find(b"\n", pos)
        line_end = data_len if line_end < 0 else line_end
        line = data[pos:line_end]
        if line.startswith(b">"):
            tok = line[1:].split()
            name = tok[0].decode() if tok else f"seq{len(out)}"
            out.append(FaiEntry(name, offset=line_end + 1))
        elif name is not None:
            entry = out[-1]
            if entry.line_bases == 0:
                entry.line_bases = len(line)
                entry.line_bytes = len(line) + 1
            entry.length += len(line)
        pos = line_end + 1
    with open(path + ".fai", "w") as handle:
        handle.writelines(entry.row() for entry in out)


def load_truth(tsv):
    """Map each read name to (contig, position, strand, length) from a truth TSV."""
    truth = {}
    with open(tsv) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[TruthCol.NAME] != "read_name":
                truth[fields[TruthCol.NAME]] = (
                    fields[TruthCol.REF],
                    int(fields[TruthCol.POS]),
                    fields[TruthCol.STRAND],
                    int(fields[TruthCol.LEN]),
                )
    return truth


def section_map(bench, mode, work, threads, quick):
    """Map simulated read sets with every mapper and score validity and accuracy against truth."""
    print(f"\n=== map ({mode})", flush=True)
    map_dir = os.path.join(work, "map")
    os.makedirs(map_dir, exist_ok=True)
    sal = os.path.join(TEMP_DIR, "salmonella_LT2.fa")
    gen = os.path.join(HERE, "gen_reads.py")
    sets = {
        # name -> (reference, reads, truth, kind)
        "salmonella_150bp": (
            sal,
            {
                "n": 50000,
                "length": 150,
                "e": 0.006,
                "ins": 0.002,
                "dele": 0.002,
                "s": 2,
            },
            "short",
        ),
        "salmonella_ONT_20kb": (
            sal,
            {
                "n": 400,
                "length": 20000,
                "sd": 0.5,
                "e": 0.04,
                "ins": 0.03,
                "dele": 0.03,
                "s": 7,
            },
            "ont",
        ),
        "salmonella_HiFi_15kb": (
            sal,
            {
                "n": 400,
                "length": 15000,
                "sd": 0.3,
                "e": 0.002,
                "ins": 0.001,
                "dele": 0.001,
                "s": 8,
            },
            "hifi",
        ),
    }
    if not quick:
        syn = os.path.join(work, "syn100m.fa")
        if not os.path.exists(syn):
            rng = random.Random(SYN_SEED)
            write_fasta(
                syn,
                [
                    (f"chr{i}", "".join(rng.choice("ACGT") for _ in range(SYN_CHROM_LEN)))
                    for i in range(SYN_CHROMS)
                ],
            )
        sets["synthetic100M_150bp"] = (
            syn,
            {"n": 300000, "length": 150, "e": 0.01, "s": 1},
            "short",
        )
    files = {}
    for name, (ref, gen_params, kind) in sets.items():
        reads_fa, tsv = os.path.join(map_dir, f"{name}.fa"), os.path.join(map_dir, f"{name}.tsv")
        if not os.path.exists(tsv):
            cmd = [
                sys.executable,
                gen,
                "-r",
                ref,
                "-n",
                str(gen_params["n"]),
                "-l",
                str(gen_params["length"]),
                "-e",
                str(gen_params["e"]),
                "--rc",
                "0.5",
                "-s",
                str(gen_params["s"]),
                "-o",
                reads_fa,
                "--gt",
                tsv,
            ]
            if "ins" in gen_params:
                cmd += ["--ins", str(gen_params["ins"]), "--del", str(gen_params["dele"])]
            if "sd" in gen_params:
                cmd += ["--len-sd", str(gen_params["sd"])]
            subprocess.run(cmd, check=True)
        files[name] = (ref, reads_fa, load_truth(tsv), kind)
    # Chimeric long reads are ONT reads with 5-20 kb of foreign sequence appended.
    ch_fa = os.path.join(map_dir, "salmonella_ONT_chimeric.fa")
    ont_ref, ont_fa, ont_truth, _ = files["salmonella_ONT_20kb"]
    if not os.path.exists(ch_fa):
        rng = random.Random(CHIMERA_SEED)
        write_fasta(
            ch_fa,
            [
                (
                    name,
                    seq
                    + "".join(
                        rng.choice("ACGT") for _ in range(rng.randrange(*CHIMERA_TAIL_RANGE))
                    ),
                )
                for name, seq in read_fasta(ont_fa)
            ],
        )
    files["salmonella_ONT_chimeric"] = (ont_ref, ch_fa, ont_truth, "ont")
    # Negative control, random reads that must not map
    neg = os.path.join(map_dir, "random_150bp.fa")
    if not os.path.exists(neg):
        rng = random.Random(NEG_SEED)
        write_fasta(
            neg,
            [
                (f"x{i}", "".join(rng.choice("ACGT") for _ in range(NEG_READ_LEN)))
                for i in range(NEG_READS)
            ],
        )
    files["random_150bp"] = (sal, neg, {f"x{i}": None for i in range(NEG_READS)}, "short")
    for ref, reads_fa, _truth, _kind in files.values():
        for path in (ref, reads_fa):  # wfmash needs .fai indexes
            if not os.path.exists(path + ".fai"):
                write_fai(path)
    preset = {"short": "sr", "ont": "map-ont", "hifi": "map-hifi"}
    # k and w of those presets, from minimap2 v2.31 options.c.
    thread_set = [1, threads] if mode == "parity" else [threads]
    for name, (ref, reads_fa, truth, kind) in files.items():
        for thread_count in thread_set:
            # In both modes optik indexes the same minimizers as minimap2's preset, with its k and w, and each mapper keeps its own other settings.
            preset_k, preset_w = MM2_PRESET_KW[preset[kind]]
            tools = {
                "optik": (
                    [
                        BIN["optik"],
                        "map",
                        "-r",
                        ref,
                        "-q",
                        reads_fa,
                        "-t",
                        str(thread_count),
                        "-k",
                        str(preset_k),
                        "-w",
                        str(preset_w),
                    ],
                    f"-k {preset_k} -w {preset_w} (minimap2 -ax {preset[kind]}'s seeds)",
                ),
            }
            mm_cmd = [BIN["minimap2"], "-ax", preset[kind], "-t", str(thread_count)]
            tools["minimap2"] = (
                mm_cmd + [ref, reads_fa],
                f"-ax {preset[kind]} (preset: k {preset_k}, w {preset_w})",
            )
            if kind != "short":
                # Wavefront-based mapper with WFA alignment.
                # Parity uses the read mode -S 1 -p 70 that run_validations validates, and optimal its defaults.
                # -k sets the same k-mer length, while its -w is a sketch window in bases rather than a minimizer window and keeps its default.
                read_opts = ["-S", "1", "-p", "70"] if mode == "parity" else []
                tools["wfmash"] = (
                    [BIN["wfmash"], ref, reads_fa, "-a", "-k", str(preset_k), "-t", str(thread_count)] + read_opts,
                    f"-a -k {preset_k} "
                    + ("-S 1 -p 70 (read mode)" if read_opts else "(otherwise defaults)"),
                )
            for tool, (cmd, params) in tools.items():
                sam = os.path.join(map_dir, f"{name}_{tool}_{mode}_t{thread_count}.sam")
                # Each mapper writes its SAM in every timed run, so all do the same work, and the gate reads the SAM of the last measured run.
                result = bench.measure(
                    f"{mode}/map/{name}/{tool}/t{thread_count}",
                    cmd,
                    cwd=map_dir,
                    stdout=sam,
                    inputs=[ref, reads_fa],
                )
                rc_, err = (0 if result.succeeded else result.code), result.err
                if rc_ != 0:
                    bench.add(
                        section="map",
                        mode=mode,
                        dataset=name,
                        tool=tool,
                        threads=thread_count,
                        gate=False,
                        gate_detail=f"exit {rc_}: {err[-DETAIL_CHARS:]}",
                        params=params,
                    )
                    continue
                chk = sam_check(
                    sam,
                    ref,
                    truth,
                    interval=kind != "short",
                    reads=dict(read_fasta(reads_fa)) if tool == "optik" else None,
                )
                gate = chk["valid"]
                det = (
                    ""
                    if gate
                    else f"{chk['invalid_count']} invalid: {chk['invalid_examples']}"
                )
                spec_note = ""
                if chk["spec_count"]:
                    multi = sum(1 for note in chk["spec_notes"] if "two primaries" in note)
                    spec_note = (
                        f"{chk['spec_count']} SAM-convention notes, e.g. {chk['spec_notes']}"
                        + (
                            " (one primary per mapped segment, none flagged supplementary; "
                            "the longest is judged)"
                            if multi
                            else ""
                        )
                    )
                # The same SAM also goes through Picard ValidateSamFile and samtools calmd.
                ext_ok, ext = external_sam_checks(sam, ref)
                if not ext_ok:
                    gate = False
                    det = (det + "; " if det else "") + ext
                if name == "random_150bp" and chk["mapped_pct"] > MAX_RANDOM_MAPPED_PCT:
                    # Random sequence must not map, so this row FAILS.
                    # It used to only append a note and still pass.
                    gate = False
                    det += f" negative control: {chk['mapped_pct']}% of random reads mapped"
                bench.add(
                    section="map",
                    mode=mode,
                    dataset=name,
                    tool=tool,
                    threads=thread_count,
                    wall_s=result.wall_row(),
                    rss_mb=round(result.rss(), 1),
                    mapped_pct=chk["mapped_pct"],
                    correct_pct=chk["correct_pct"],
                    wrong=chk["wrong"],
                    wrong_mapq30=chk["wrong_mapq30"],
                    nm_rate=chk["nm_per_aligned_base"],
                    sam_spec=spec_note or "conforms",
                    gate=gate,
                    gate_detail=det,
                    params=params,
                )
