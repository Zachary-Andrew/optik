#!/usr/bin/env python3
"""validate_tags.py - the MD and cs tags and QUAL of `optik map` checked against independent oracles.

  1. samtools calmd recomputes MD and NM from the reference for every record and warns when a record's own tag differs, and there must be no warnings.
  2. Every record's MD and both forms of cs are replayed against the reference and SEQ by an independent parser.
     Mismatches must name the reference bases, and for cs the read bases, at their positions, deletions must spell the reference and insertions the read, runs must count exactly the identical bases, and each tag must consume exactly the CIGAR's aligned query and reference.
  3. Wherever minimap2 (--MD --cs) reports the same position and CIGAR for a read, optik's MD and cs strings must equal minimap2's.

The inputs cover the short-read path (150 bp, 1% substitutions and indels), the chain path (long reads with 5% errors, soft clips and supplementary records) and reads carrying N.

Usage: python3 bench/validate_tags.py --ref ref.fa --workdir DIR [--samtools S] [--minimap2 M] Exit status 0 means all checks passed.
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from typing import NamedTuple

from checklog import CheckLog
from config import PATHS
from diagnostics import DETAIL_CHARS
from dna import read_fasta
from sam_spec import (
    SAM_TAG_NAME_LEN,
    SAM_TAG_PREFIX_LEN,
    SamCol,
    SamFlag,
)

DEFAULT_THREADS = 4
# problem records listed with a failed replay check, and QUAL characters shown per read in the bad-quality check
MAX_EXAMPLES = 3
# Short failure details quote this many examples.
SHOWN_EXAMPLES = 2
QUAL_PREVIEW_CHARS = 3
# chimeric long reads: CHIMERAS reads joining CHIMERA_HALF_LEN bases from each half of the genome, starts kept CHIMERA_MARGIN bases from a half's end
CHIMERA_SEED = 3
CHIMERAS = 40
CHIMERA_HALF_LEN = 3000
CHIMERA_MARGIN = 4000
# one short read in N_READ_EVERY gets an N at N_POS_MIN + (index mod N_POS_SPREAD)
N_READ_EVERY = 50
N_POS_MIN = 20
N_POS_SPREAD = 100
# Among the QUAL check reads, the first QUAL_SHORT_READS are QUAL_SHORT_LEN bp, then plain QUAL_LONG_LEN bp reads up to QUAL_PLAIN_READS, then chimeras of two CHIMERA_HALF_LEN halves
QUAL_SEED = 6
QUAL_READS = 2000
QUAL_SHORT_READS = 1960
QUAL_PLAIN_READS = 1980
QUAL_SHORT_LEN = 150
QUAL_LONG_LEN = 6000
REVCOMP_PROB = 0.5
# synthetic qualities cycle through QUAL_LEVELS Phred+33 characters with stride QUAL_STRIDE
PHRED_OFFSET = 33
QUAL_STRIDE = 7
QUAL_LEVELS = 41

check = CheckLog()
FAILS = check.failed


class Column(NamedTuple):
    """One aligned column, with its CIGAR op and the read and reference bases (None on the gapped side)."""

    op: str
    read: str | None
    ref: str | None


def norm_ref(seq):
    """Upper-case s with U read as T and every other non-ACGT base as N."""
    return "".join(base if base in "ACGT" else ("T" if base == "U" else "N") for base in seq.upper())


CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def replay(rec, ref):
    """Problems found when replaying one record's CIGAR, MD and cs against the reference and SEQ.
    An empty list means the record is consistent."""
    fields = rec.split("\t")
    if int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
        return []
    cig, seq, pos = fields[SamCol.CIGAR], fields[SamCol.SEQ], int(fields[SamCol.POS]) - 1
    tags = {tag_field[:SAM_TAG_NAME_LEN]: tag_field[SAM_TAG_PREFIX_LEN:] for tag_field in fields[SamCol.TAGS:]}
    ref_seq = ref[fields[SamCol.RNAME]]
    errs = []
    # Aligned columns from the CIGAR
    cols = []
    read_pos, ref_pos = 0, pos
    for op_len, cigar_op in CIG.findall(cig):
        op_len = int(op_len)
        if cigar_op == "S":
            read_pos += op_len
        elif cigar_op in "M=X":
            for _ in range(op_len):
                cols.append(Column("M", seq[read_pos], ref_seq[ref_pos]))
                read_pos += 1
                ref_pos += 1
        elif cigar_op == "I":
            for _ in range(op_len):
                cols.append(Column("I", seq[read_pos], None))
                read_pos += 1
        elif cigar_op == "D":
            for _ in range(op_len):
                cols.append(Column("D", None, ref_seq[ref_pos]))
                ref_pos += 1

    def same(read_base, ref_base):
        """True for a match of two ACGT bases."""
        return read_base == ref_base and ref_base in "ACGT"

    # NM
    nm_count = sum(1 for col_op, read_base, ref_base in cols if col_op != "M" or not same(read_base, ref_base))
    if "NM" in tags and int(tags["NM"]) != nm_count:
        errs.append(f"NM {tags['NM']} != {nm_count}")
    # MD, the expected string in SAM spec form built from the columns
    want, run = [], 0
    for i, (col_op, read_base, ref_base) in enumerate(cols):
        if col_op == "M":
            if same(read_base, ref_base):
                run += 1
            else:
                want.append(str(run))
                want.append(ref_base)
                run = 0
        elif col_op == "D":
            if i == 0 or cols[i - 1].op != "D":
                want.append(str(run))
                want.append("^")
                run = 0
            want.append(ref_base)
    want.append(str(run))
    want = "".join(want)
    if tags.get("MD") != want:
        errs.append(f"MD {tags.get('MD')} != {want}")
    # cs, replayed against the columns
    if "cs" in tags:
        cs_tag = tags["cs"]
        # The tag must split into tokens with nothing left over, so a malformed cs fails here instead of being half-checked.
        toks = re.findall(r"(:[0-9]+|=[A-Za-z]+|\*[a-z][a-z]|[-+][a-z]+)", cs_tag)
        if "".join(toks) != cs_tag:
            errs.append(f"cs unparseable: {cs_tag}")
        else:
            # Each token then consumes alignment columns and has to match what those columns hold.
            col_idx = 0
            for token in toks:
                kind, body = token[0], token[1:]
                # Short-form :N gives only a length, and each of those N columns has to be an identical M column.
                if kind == ":":
                    for _ in range(int(body)):
                        if (
                            col_idx >= len(cols)
                            or cols[col_idx].op != "M"
                            or not same(cols[col_idx].read, cols[col_idx].ref)
                        ):
                            errs.append(f"cs {token} not identical at col {col_idx}")
                            break
                        col_idx += 1
                elif kind == "=":
                    for base in body:
                        if (
                            col_idx >= len(cols)
                            or cols[col_idx].op != "M"
                            or not same(cols[col_idx].read, cols[col_idx].ref)
                            or cols[col_idx].read != base
                        ):
                            errs.append(f"cs {token} mismatch at col {col_idx}")
                            break
                        col_idx += 1
                # A substitution is one column, given as reference base then read base.
                elif kind == "*":
                    if (
                        col_idx >= len(cols)
                        or cols[col_idx].op != "M"
                        or cols[col_idx].ref.lower() != body[0]
                        or cols[col_idx].read.lower() != body[1]
                    ):
                        errs.append(f"cs {token} wrong at col {col_idx}")
                    col_idx += 1
                elif kind == "+":
                    for base in body:
                        if (
                            col_idx >= len(cols)
                            or cols[col_idx].op != "I"
                            or cols[col_idx].read.lower() != base
                        ):
                            errs.append(f"cs {token} wrong at col {col_idx}")
                            break
                        col_idx += 1
                elif kind == "-":
                    for base in body:
                        if (
                            col_idx >= len(cols)
                            or cols[col_idx].op != "D"
                            or cols[col_idx].ref.lower() != base
                        ):
                            errs.append(f"cs {token} wrong at col {col_idx}")
                            break
                        col_idx += 1
            # Tokens must account for every column, so a truncated cs tag fails too.
            if col_idx != len(cols):
                errs.append(f"cs covers {col_idx} of {len(cols)} columns")
    return errs


def run(cmd, **kwargs):
    proc = subprocess.run([str(arg) for arg in cmd], capture_output=True, **kwargs, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"failed: {' '.join(map(str, cmd))}\n{proc.stderr.decode()[-DETAIL_CHARS:]}"
        )
    return proc


def _m_columns(cigar, seq):
    """Read index of every aligned (M, = or X) column, in order."""
    cols, read_pos = [], 0
    for op_len, cigar_op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        op_len = int(op_len)
        if cigar_op in "M=X":
            cols.extend(range(read_pos, read_pos + op_len))
            read_pos += op_len
        elif cigar_op in "IS":
            read_pos += op_len
    return cols


def heal_nn(tag, value, cigar, seq):
    """The tag with every column of a reference N or read N written as a match, minimap2's convention, as by_key explains."""
    cols = _m_columns(cigar, seq)

    def isn(col):
        """True when aligned column j holds an N in the read."""
        return col < len(cols) and cols[col] < len(seq) and seq[cols[col]] in "Nn"

    if tag == "MD":
        out, run, col_idx = [], 0, 0
        for tok in re.findall(r"\d+|\^[A-Za-z]+|[A-Za-z]", value):
            if tok.isdigit():
                run += int(tok)
                col_idx += int(tok)
            elif tok[0] == "^":
                out.append(f"{run}{tok}")
                run = 0
            elif tok.upper() == "N" and isn(col_idx):
                run += 1
                col_idx += 1
            else:
                out.append(f"{run}{tok}")
                run = 0
                col_idx += 1
        return "".join(out) + str(run)
    # cs in short (":n") or long ("=ACGT") form
    long_form = "=" in value
    ops = re.findall(r"(:\d+|=[A-Za-z]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)", value)
    out, col_idx = [], 0

    def push_match(base):
        """Append one matching base to the cs string being built, merging runs."""
        if long_form:
            if out and out[-1][0] == "=":
                out[-1] += base
            else:
                out.append("=" + base)
        else:
            if out and out[-1][0] == ":":
                out[-1] = ":" + str(int(out[-1][1:]) + 1)
            else:
                out.append(":1")

    for cs_op in ops:
        if cs_op[0] == ":":
            run_len = int(cs_op[1:])
            if out and out[-1][0] == ":":
                out[-1] = ":" + str(int(out[-1][1:]) + run_len)
            else:
                out.append(cs_op)
            col_idx += run_len
        elif cs_op[0] == "=":
            for base in cs_op[1:]:
                push_match(base)
            col_idx += len(cs_op) - 1
        elif cs_op[0] == "*":
            _, ref_char, read_char = cs_op
            if ref_char == "n" and read_char == "n" and isn(col_idx):
                push_match("N")
            else:
                out.append(cs_op)
            col_idx += 1
        else:
            out.append(cs_op)
    return "".join(out)


def main():
    """Map the read sets and run every tag and QUAL check."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--ref", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--samtools", default=shutil.which("samtools"))
    parser.add_argument("--minimap2", default=str(PATHS["minimap2"]))
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    opts = parser.parse_args()
    work_dir = Path(opts.workdir)
    work_dir.mkdir(parents=True, exist_ok=True)
    optik = PATHS["optik"]
    # The reference is normalised and rewritten, so optik and the replay check see identical bases.
    ref = {name: norm_ref(seq) for name, seq in read_fasta(opts.ref)}
    refp = work_dir / "ref.fa"
    with open(refp, "w") as handle:
        handle.writelines(f">{name}\n{seq}\n" for name, seq in ref.items())
    gen = HERE / "gen_reads.py"
    # The --ins and --del rates put gap operations into cs and MD, not only substitutions.
    sets = {
        "short": [
            "-n",
            30000,
            "-l",
            150,
            "-e",
            0.01,
            "--ins",
            0.002,
            "--del",
            0.002,
            "-s",
            21,
        ],
        "long": [
            "-n",
            400,
            "-l",
            6000,
            "--len-sd",
            0.4,
            "-e",
            0.03,
            "--ins",
            0.01,
            "--del",
            0.01,
            "-s",
            22,
        ],
    }
    for name, args in sets.items():
        reads = work_dir / f"{name}.fa"
        run([sys.executable, gen, "-r", refp, "-o", reads] + args)
        if name == "long":  # chimeras, which produce supplementary records
            genome_seq = next(iter(ref.values()))
            import random as random_mod

            rng = random_mod.Random(CHIMERA_SEED)
            with open(reads, "a") as handle:
                for i in range(CHIMERAS):
                    left_start = rng.randrange(0, len(genome_seq) // 2 - CHIMERA_MARGIN)
                    right_start = rng.randrange(len(genome_seq) // 2, len(genome_seq) - CHIMERA_MARGIN)
                    handle.write(f">chim{i}\n{genome_seq[left_start : left_start + CHIMERA_HALF_LEN]}{genome_seq[right_start : right_start + CHIMERA_HALF_LEN]}\n")
        if name == "short":  # reads with N
            recs = list(read_fasta(reads))
            with open(reads, "w") as handle:
                for i, (rec_name, rec_seq) in enumerate(recs):
                    # Every N_READ_EVERY-th short read gets one N, which MD and cs must report as a mismatch.
                    if i % N_READ_EVERY == 0:
                        n_pos = N_POS_MIN + i % N_POS_SPREAD
                        rec_seq = rec_seq[:n_pos] + "N" + rec_seq[n_pos + 1 :]
                    handle.write(f">{rec_name}\n{rec_seq}\n")
        # --cs writes identical runs as :N while --cs=long spells the bases out with =, and both forms are replayed.
        for mode in ("--cs", "--cs=long"):
            sam = work_dir / f"{name}{mode.replace('=', '_')}.sam"
            run(
                [
                    optik,
                    "map",
                    "-r",
                    refp,
                    "-q",
                    reads,
                    "-t",
                    opts.threads,
                    "-o",
                    sam,
                    mode,
                ]
            )
            with open(sam) as handle:
                recs = [line.rstrip("\n") for line in handle if not line.startswith("@")]
            mapped = [
                line
                for line in recs
                if not int(line.split("\t")[SamCol.FLAG]) & SamFlag.UNMAPPED
            ]
            # Replaying rebuilds each record's alignment from its tags and compares it with the reference and SEQ.
            bad = [
                (line.split("\t")[SamCol.QNAME], problem) for line in mapped for problem in replay(line, ref)
            ]
            have = sum(1 for line in mapped if "\tMD:Z:" in line and "\tcs:Z:" in line)
            check(
                f"{name} {mode}: every mapped record carries MD and cs",
                have == len(mapped) and mapped,
                f"{have}/{len(mapped)}",
            )
            check(
                f"{name} {mode}: MD, cs and NM replay exactly against reference and SEQ",
                not bad,
                f"records={len(mapped)} problems={len(bad)} {bad[:MAX_EXAMPLES]}",
            )
            # samtools calmd recomputes MD and NM from the reference, giving a second opinion on optik's tags.
            if mode == "--cs":
                proc = subprocess.run(
                    [opts.samtools, "calmd", "-b", sam, refp],
                    capture_output=True,
                    check=False,
                )
                # calmd warns about every record whose MD or NM differs from its own recomputation.
                warn = [
                    line
                    for line in proc.stderr.decode().splitlines()
                    if "different" in line.lower()
                ]
                check(
                    f"{name}: samtools calmd recomputes the same MD and NM for every record",
                    proc.returncode == 0 and not warn,
                    f"warnings={len(warn)} {warn[:SHOWN_EXAMPLES]}",
                )
                supp = sum(
                    1
                    for line in mapped
                    if int(line.split("\t")[SamCol.FLAG]) & SamFlag.SUPPLEMENTARY
                )
                clip = sum(1 for line in mapped if "S" in line.split("\t")[SamCol.CIGAR])
                print(
                    f"    ({len(mapped)} mapped records; {supp} supplementary, {clip} with soft clips)"
                )
        # minimap2 on the same reads, where identical alignments must give identical tags.
        # In SAM mode minimap2 writes only one of MD and cs per run, so each tag form is compared with its own minimap2 run.
        preset = "sr" if name == "short" else "map-ont"

        def by_key(path, tag, heal=False):
            """Map (name, ref, pos, CIGAR) to the tag, healing N against N to a match first when asked.

            optik scores a reference N against a read N as a mismatch, as samtools calmd does when it recomputes minimap2's MD.
            minimap2 counts it as a match.
            Only that kind of column is healed, so any other difference still fails."""
            out = {}
            with open(path) as handle:
                for line in handle:
                    if line.startswith("@"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if int(fields[SamCol.FLAG]) & (SamFlag.UNMAPPED | SamFlag.SECONDARY):
                        continue
                    tags = {tag_field[:SAM_TAG_NAME_LEN]: tag_field[SAM_TAG_PREFIX_LEN:] for tag_field in fields[SamCol.TAGS:]}
                    value = tags.get(tag)
                    if heal and value is not None:
                        value = heal_nn(tag, value, fields[SamCol.CIGAR], fields[SamCol.SEQ])
                    out[(fields[SamCol.QNAME], fields[SamCol.RNAME], fields[SamCol.POS], fields[SamCol.CIGAR])] = value
            return out

        for tag, mmopt, ours in (
            ("MD", "--MD", "--cs"),
            ("cs", "--cs", "--cs"),
            ("cs", "--cs=long", "--cs_long"),
        ):
            mm_sam = work_dir / f"{name}.mm{mmopt.replace('=', '_')}.sam"
            with open(mm_sam, "wb") as handle:
                subprocess.run(
                    [
                        str(opts.minimap2),
                        "-ax",
                        preset,
                        mmopt,
                        "-t",
                        str(opts.threads),
                        str(refp),
                        str(reads),
                    ],
                    stdout=handle,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            ours_tags, mm_tags = by_key(work_dir / f"{name}{ours}.sam", tag, heal=True), by_key(mm_sam, tag)
            common = [key for key in set(ours_tags) & set(mm_tags) if mm_tags[key] is not None]
            diff = [key for key in common if ours_tags[key] != mm_tags[key]]
            check(
                f"{name}: {tag} ({mmopt}) equals minimap2's wherever position and CIGAR agree",
                len(common) > 0 and not diff,
                f"same alignments={len(common)} differing={len(diff)} {[(key, ours_tags[key], mm_tags[key]) for key in diff[:SHOWN_EXAMPLES]]}",
            )
    # QUAL copies the FASTQ base qualities into SAM, reversed on reverse-strand records and cut by the hard clips on supplementary ones, while FASTA reads keep '*'.
    fq_path = work_dir / "qual.fq"
    comp = str.maketrans("ACGT", "TGCA")
    import random as random_mod

    rng = random_mod.Random(QUAL_SEED)
    genome_seq = next(iter(ref.values()))
    quals = {}
    with open(fq_path, "w") as handle:
        for i in range(QUAL_READS):
            read_len = QUAL_SHORT_LEN if i < QUAL_SHORT_READS else QUAL_LONG_LEN
            start = rng.randrange(len(genome_seq) - 2 * read_len)
            seq_ = (
                genome_seq[start : start + read_len]
                if i < QUAL_PLAIN_READS
                else genome_seq[start : start + CHIMERA_HALF_LEN]
                + genome_seq[(start + len(genome_seq) // 2) % (len(genome_seq) - CHIMERA_MARGIN) :][:CHIMERA_HALF_LEN]
            )
            if rng.random() < REVCOMP_PROB:
                seq_ = seq_.translate(comp)[::-1]
            qual = "".join(chr(PHRED_OFFSET + ((j * QUAL_STRIDE + i) % QUAL_LEVELS)) for j in range(len(seq_)))
            quals[f"q{i}"] = qual
            handle.write(f"@q{i}\n{seq_}\n+\n{qual}\n")
    osam = work_dir / "qual.sam"
    run([optik, "map", "-r", refp, "-q", fq_path, "-t", opts.threads, "-o", osam])
    bad = rec_name = 0
    with open(osam) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if int(fields[SamCol.FLAG]) & SamFlag.UNMAPPED:
                continue
            rec_name += 1
            want = quals[fields[SamCol.QNAME]][::-1] if int(fields[SamCol.FLAG]) & SamFlag.REVERSE else quals[fields[SamCol.QNAME]]
            # Hard-clipped (supplementary) records cut QUAL as SEQ is cut, by the CIGAR's H.
            clip_match = re.match(r"^(\d+)H", fields[SamCol.CIGAR])
            lead = int(clip_match.group(1)) if clip_match else 0
            clip_match = re.search(r"(\d+)H$", fields[SamCol.CIGAR])
            trail = int(clip_match.group(1)) if clip_match else 0
            want = want[lead : len(want) - trail]
            bad += fields[SamCol.QUAL] != want or len(fields[SamCol.QUAL]) != len(fields[SamCol.SEQ])
    check(
        "FASTQ qualities are written to QUAL (reversed on reverse-strand, cut by hard clips on supplementary records)",
        rec_name > 0 and bad == 0,
        f"records={rec_name} wrong={bad}",
    )
    mm_sam = work_dir / "qual.mm.sam"
    with open(mm_sam, "wb") as handle:
        subprocess.run(
            [str(opts.minimap2), "-ax", "sr", "-t", str(opts.threads), str(refp), str(fq_path)],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    def prim(path):
        """Map (read, strand) of each primary record to its QUAL."""
        out = {}
        with open(path) as handle:
            for line in handle:
                if line.startswith("@"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if not int(fields[SamCol.FLAG]) & SamFlag.NOT_PRIMARY:
                    out[(fields[SamCol.QNAME], int(fields[SamCol.FLAG]) & SamFlag.REVERSE)] = fields[SamCol.QUAL]
        return out

    optik_quals, mm_quals = prim(osam), prim(mm_sam)
    # minimap2 reverses QUAL on reverse records too. Equal strings for the same read and strand prove the same orientation rule.
    common = set(optik_quals) & set(mm_quals)
    check(
        "QUAL equals minimap2's for the same read and strand",
        len(common) > 0 and all(optik_quals[key] == mm_quals[key] for key in common),
        f"compared={len(common)}",
    )
    # Unusable qualities (wrong length, or a byte outside Phred+33 '!'..'~') give QUAL '*' and a warning, never a SAM record that samtools rejects.
    badq = work_dir / "badqual.fq"
    base_seq = genome_seq[1000:1150]
    with open(badq, "w") as handle:
        handle.write(
            f"@len\n{base_seq}\n+\n{'I' * 149}\n@byte\n{base_seq}\n+\n{'I' * 70}\x7f{'I' * 79}\n"
            f"@good\n{base_seq}\n+\n{'I' * 150}\n"
        )
    bsam = work_dir / "badqual.sam"
    proc = subprocess.run(
        [str(optik), "map", "-r", str(refp), "-q", str(badq), "-o", str(bsam)],
        capture_output=True,
        text=True,
        check=False,
    )
    with open(bsam) as handle:
        got = {
            fields[SamCol.QNAME]: fields[SamCol.QUAL]
            for fields in (line.rstrip("\n").split("\t") for line in handle)
            if fields[SamCol.QNAME][0] != "@"
        }
    # samtools must count all three records, so no QUAL broke the SAM.
    view_proc = subprocess.run(
        [str(opts.samtools), "view", "-c", str(bsam)],
        capture_output=True,
        text=True,
        check=False,
    )
    check(
        "unusable FASTQ qualities give QUAL '*' with a warning; samtools reads the SAM",
        proc.returncode == 0
        and got == {"len": "*", "byte": "*", "good": "I" * 150}
        and "2 FASTQ records had unusable qualities" in proc.stderr
        and view_proc.returncode == 0
        and view_proc.stdout.strip() == "3",
        f"qual={ {name: qual[:QUAL_PREVIEW_CHARS] for name, qual in got.items()} } samtools={view_proc.returncode}",
    )
    print(
        f"\nvalidate_tags: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}"
    )
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
