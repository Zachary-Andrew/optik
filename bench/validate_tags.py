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
from config import PATHS
from sam_spec import (
    SAM_FLAG_REVERSE,
    SAM_FLAG_SUPPLEMENTARY,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_NOT_PRIMARY,
)

FAILS = []


def check(name, ok, detail=""):
    """Report one check, appending its name to FAILS on failure."""
    print(
        f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""),
        flush=True,
    )
    if not ok:
        FAILS.append(name)


def read_fasta(path):
    """(name, sequence) of each FASTA record."""
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def norm_ref(s):
    """Upper-case s with U read as T and every other non-ACGT base as N."""
    return "".join(c if c in "ACGT" else ("T" if c == "U" else "N") for c in s.upper())


CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def replay(rec, ref):
    """Problems found when replaying one record's CIGAR, MD and cs against the reference and SEQ.
    An empty list means the record is consistent."""
    f = rec.split("\t")
    if int(f[1]) & SAM_FLAG_UNMAPPED:
        return []
    cig, seq, pos = f[5], f[9], int(f[3]) - 1
    tags = {t[:2]: t[5:] for t in f[11:]}
    r = ref[f[2]]
    errs = []
    # Aligned columns from the CIGAR
    cols = []  # op, query base or None, reference base or None
    qi, ti = 0, pos
    for n, op in CIG.findall(cig):
        n = int(n)
        if op == "S":
            qi += n
        elif op in "M=X":
            for _ in range(n):
                cols.append(("M", seq[qi], r[ti]))
                qi += 1
                ti += 1
        elif op == "I":
            for _ in range(n):
                cols.append(("I", seq[qi], None))
                qi += 1
        elif op == "D":
            for _ in range(n):
                cols.append(("D", None, r[ti]))
                ti += 1

    def same(q, t):
        """True for a match of two ACGT bases."""
        return q == t and t in "ACGT"

    # NM
    nm = sum(1 for o, q, t in cols if o != "M" or not same(q, t))
    if "NM" in tags and int(tags["NM"]) != nm:
        errs.append(f"NM {tags['NM']} != {nm}")
    # MD, the expected string in SAM spec form built from the columns
    want, run = [], 0
    for i, (o, q, t) in enumerate(cols):
        if o == "M":
            if same(q, t):
                run += 1
            else:
                want.append(str(run))
                want.append(t)
                run = 0
        elif o == "D":
            if i == 0 or cols[i - 1][0] != "D":
                want.append(str(run))
                want.append("^")
                run = 0
            want.append(t)
    want.append(str(run))
    want = "".join(want)
    if tags.get("MD") != want:
        errs.append(f"MD {tags.get('MD')} != {want}")
    # cs, replayed against the columns
    if "cs" in tags:
        cs = tags["cs"]
        toks = re.findall(r"(:[0-9]+|=[A-Za-z]+|\*[a-z][a-z]|[-+][a-z]+)", cs)
        if "".join(toks) != cs:
            errs.append(f"cs unparseable: {cs}")
        else:
            c = 0
            for tk in toks:
                kind, body = tk[0], tk[1:]
                if kind == ":":
                    for _ in range(int(body)):
                        if (
                            c >= len(cols)
                            or cols[c][0] != "M"
                            or not same(cols[c][1], cols[c][2])
                        ):
                            errs.append(f"cs {tk} not identical at col {c}")
                            break
                        c += 1
                elif kind == "=":
                    for b in body:
                        if (
                            c >= len(cols)
                            or cols[c][0] != "M"
                            or not same(cols[c][1], cols[c][2])
                            or cols[c][1] != b
                        ):
                            errs.append(f"cs {tk} mismatch at col {c}")
                            break
                        c += 1
                elif kind == "*":
                    if (
                        c >= len(cols)
                        or cols[c][0] != "M"
                        or cols[c][2].lower() != body[0]
                        or cols[c][1].lower() != body[1]
                    ):
                        errs.append(f"cs {tk} wrong at col {c}")
                    c += 1
                elif kind == "+":
                    for b in body:
                        if (
                            c >= len(cols)
                            or cols[c][0] != "I"
                            or cols[c][1].lower() != b
                        ):
                            errs.append(f"cs {tk} wrong at col {c}")
                            break
                        c += 1
                elif kind == "-":
                    for b in body:
                        if (
                            c >= len(cols)
                            or cols[c][0] != "D"
                            or cols[c][2].lower() != b
                        ):
                            errs.append(f"cs {tk} wrong at col {c}")
                            break
                        c += 1
            if c != len(cols):
                errs.append(f"cs covers {c} of {len(cols)} columns")
    return errs


def run(cmd, **kw):
    p = subprocess.run([str(c) for c in cmd], capture_output=True, **kw, check=False)
    if p.returncode != 0:
        raise RuntimeError(
            f"failed: {' '.join(map(str, cmd))}\n{p.stderr.decode()[-500:]}"
        )
    return p


def _m_columns(cigar, seq):
    """Read index of every aligned (M, = or X) column, in order."""
    cols, q = [], 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            cols.extend(range(q, q + n))
            q += n
        elif op in "IS":
            q += n
    return cols


def heal_nn(tag, v, cigar, seq):
    """The tag with every column of a reference N or read N written as a match, minimap2's convention, as by_key explains."""
    cols = _m_columns(cigar, seq)

    def isn(j):
        """True when aligned column j holds an N in the read."""
        return j < len(cols) and cols[j] < len(seq) and seq[cols[j]] in "Nn"

    if tag == "MD":
        out, run, j = [], 0, 0
        for tok in re.findall(r"\d+|\^[A-Za-z]+|[A-Za-z]", v):
            if tok.isdigit():
                run += int(tok)
                j += int(tok)
            elif tok[0] == "^":
                out.append(f"{run}{tok}")
                run = 0
            elif tok.upper() == "N" and isn(j):
                run += 1
                j += 1
            else:
                out.append(f"{run}{tok}")
                run = 0
                j += 1
        return "".join(out) + str(run)
    # cs in short (":n") or long ("=ACGT") form
    long_form = "=" in v
    ops = re.findall(r"(:\d+|=[A-Za-z]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)", v)
    out, j = [], 0

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

    for op in ops:
        if op[0] == ":":
            n = int(op[1:])
            if out and out[-1][0] == ":":
                out[-1] = ":" + str(int(out[-1][1:]) + n)
            else:
                out.append(op)
            j += n
        elif op[0] == "=":
            for b in op[1:]:
                push_match(b)
            j += len(op) - 1
        elif op[0] == "*":
            if op[1] == "n" and op[2] == "n" and isn(j):
                push_match("N")
            else:
                out.append(op)
            j += 1
        else:
            out.append(op)
    return "".join(out)


def main():
    """Map the read sets and run every tag and QUAL check."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ref", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--samtools", default=shutil.which("samtools"))
    ap.add_argument("--minimap2", default=str(PATHS["minimap2"]))
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    wd = Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    optik = PATHS["optik"]
    ref = {n: norm_ref(s) for n, s in read_fasta(a.ref)}
    refp = wd / "ref.fa"
    with open(refp, "w") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in ref.items())
    gen = HERE / "gen_reads.py"
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
        reads = wd / f"{name}.fa"
        run([sys.executable, gen, "-r", refp, "-o", reads] + args)
        if name == "long":  # chimeras, which produce supplementary records
            g = next(iter(ref.values()))
            import random as _r

            rr = _r.Random(3)
            with open(reads, "a") as fh:
                for i in range(40):
                    a0 = rr.randrange(0, len(g) // 2 - 4000)
                    b0 = rr.randrange(len(g) // 2, len(g) - 4000)
                    fh.write(f">chim{i}\n{g[a0 : a0 + 3000]}{g[b0 : b0 + 3000]}\n")
        if name == "short":  # reads with N
            recs = list(read_fasta(reads))
            with open(reads, "w") as fh:
                for i, (n, s) in enumerate(recs):
                    if i % 50 == 0:
                        p = 20 + i % 100
                        s = s[:p] + "N" + s[p + 1 :]
                    fh.write(f">{n}\n{s}\n")
        for mode in ("--cs", "--cs=long"):
            sam = wd / f"{name}{mode.replace('=', '_')}.sam"
            run(
                [
                    optik,
                    "map",
                    "-r",
                    refp,
                    "-q",
                    reads,
                    "-t",
                    a.threads,
                    "-o",
                    sam,
                    mode,
                ]
            )
            with open(sam) as fh:
                recs = [line.rstrip("\n") for line in fh if not line.startswith("@")]
            mapped = [
                line
                for line in recs
                if not int(line.split("\t")[1]) & SAM_FLAG_UNMAPPED
            ]
            bad = [
                (line.split("\t")[0], e) for line in mapped for e in replay(line, ref)
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
                f"records={len(mapped)} problems={len(bad)} {bad[:3]}",
            )
            if mode == "--cs":
                p = subprocess.run(
                    [a.samtools, "calmd", "-b", sam, refp],
                    capture_output=True,
                    check=False,
                )
                warn = [
                    line
                    for line in p.stderr.decode().splitlines()
                    if "different" in line.lower()
                ]
                check(
                    f"{name}: samtools calmd recomputes the same MD and NM for every record",
                    p.returncode == 0 and not warn,
                    f"warnings={len(warn)} {warn[:2]}",
                )
                supp = sum(
                    1
                    for line in mapped
                    if int(line.split("\t")[1]) & SAM_FLAG_SUPPLEMENTARY
                )
                clip = sum(1 for line in mapped if "S" in line.split("\t")[5])
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
            with open(path) as fh:
                for line in fh:
                    if line.startswith("@"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if int(f[1]) & (4 | 256):
                        continue
                    t = {x[:2]: x[5:] for x in f[11:]}
                    v = t.get(tag)
                    if heal and v is not None:
                        v = heal_nn(tag, v, f[5], f[9])
                    out[(f[0], f[2], f[3], f[5])] = v
            return out

        for tag, mmopt, ours in (
            ("MD", "--MD", "--cs"),
            ("cs", "--cs", "--cs"),
            ("cs", "--cs=long", "--cs_long"),
        ):
            mm = wd / f"{name}.mm{mmopt.replace('=', '_')}.sam"
            with open(mm, "wb") as fh:
                subprocess.run(
                    [
                        str(a.minimap2),
                        "-ax",
                        preset,
                        mmopt,
                        "-t",
                        str(a.threads),
                        str(refp),
                        str(reads),
                    ],
                    stdout=fh,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            ko, km = by_key(wd / f"{name}{ours}.sam", tag, heal=True), by_key(mm, tag)
            common = [k for k in set(ko) & set(km) if km[k] is not None]
            diff = [k for k in common if ko[k] != km[k]]
            check(
                f"{name}: {tag} ({mmopt}) equals minimap2's wherever position and CIGAR agree",
                len(common) > 0 and not diff,
                f"same alignments={len(common)} differing={len(diff)} {[(k, ko[k], km[k]) for k in diff[:2]]}",
            )
    # QUAL carries the FASTQ base qualities into SAM, reversed on reverse-strand records and cut by the hard clips on supplementary ones, while FASTA reads keep '*'.
    fq = wd / "qual.fq"
    comp = str.maketrans("ACGT", "TGCA")
    import random as _r

    rr = _r.Random(6)
    g = next(iter(ref.values()))
    quals = {}
    with open(fq, "w") as fh:
        for i in range(2000):
            read_len = 150 if i < 1960 else 6000
            st = rr.randrange(len(g) - 2 * read_len)
            seq_ = (
                g[st : st + read_len]
                if i < 1980
                else g[st : st + 3000]
                + g[(st + len(g) // 2) % (len(g) - 4000) :][:3000]
            )
            if rr.random() < 0.5:
                seq_ = seq_.translate(comp)[::-1]
            q = "".join(chr(33 + ((j * 7 + i) % 41)) for j in range(len(seq_)))
            quals[f"q{i}"] = q
            fh.write(f"@q{i}\n{seq_}\n+\n{q}\n")
    osam = wd / "qual.sam"
    run([optik, "map", "-r", refp, "-q", fq, "-t", a.threads, "-o", osam])
    bad = n = 0
    with open(osam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if int(f[1]) & SAM_FLAG_UNMAPPED:
                continue
            n += 1
            want = quals[f[0]][::-1] if int(f[1]) & SAM_FLAG_REVERSE else quals[f[0]]
            # Hard-clipped (supplementary) records cut QUAL as SEQ is cut, by the CIGAR's H.
            m = re.match(r"^(\d+)H", f[5])
            lead = int(m.group(1)) if m else 0
            m = re.search(r"(\d+)H$", f[5])
            trail = int(m.group(1)) if m else 0
            want = want[lead : len(want) - trail]
            bad += f[10] != want or len(f[10]) != len(f[9])
    check(
        "FASTQ qualities are written to QUAL (reversed on reverse-strand, cut by hard clips on supplementary records)",
        n > 0 and bad == 0,
        f"records={n} wrong={bad}",
    )
    mm = wd / "qual.mm.sam"
    with open(mm, "wb") as fh:
        subprocess.run(
            [str(a.minimap2), "-ax", "sr", "-t", str(a.threads), str(refp), str(fq)],
            stdout=fh,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    def prim(path):
        """Map (read, strand) of each primary record to its QUAL."""
        out = {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                f = line.rstrip("\n").split("\t")
                if not int(f[1]) & SAM_FLAGS_NOT_PRIMARY:
                    out[(f[0], int(f[1]) & SAM_FLAG_REVERSE)] = f[10]
        return out

    po, pm = prim(osam), prim(mm)
    common = set(po) & set(pm)
    check(
        "QUAL equals minimap2's for the same read and strand",
        len(common) > 0 and all(po[k] == pm[k] for k in common),
        f"compared={len(common)}",
    )
    # Unusable qualities (wrong length, or a byte outside Phred+33 '!'..'~') give QUAL '*' and a warning, never a SAM record that samtools rejects.
    badq = wd / "badqual.fq"
    s0 = g[1000:1150]
    with open(badq, "w") as fh:
        fh.write(
            f"@len\n{s0}\n+\n{'I' * 149}\n@byte\n{s0}\n+\n{'I' * 70}\x7f{'I' * 79}\n"
            f"@good\n{s0}\n+\n{'I' * 150}\n"
        )
    bsam = wd / "badqual.sam"
    p = subprocess.run(
        [str(optik), "map", "-r", str(refp), "-q", str(badq), "-o", str(bsam)],
        capture_output=True,
        text=True,
        check=False,
    )
    with open(bsam) as fh:
        got = {
            f[0]: f[10]
            for f in (line.rstrip("\n").split("\t") for line in fh)
            if f[0][0] != "@"
        }
    v = subprocess.run(
        [str(a.samtools), "view", "-c", str(bsam)],
        capture_output=True,
        text=True,
        check=False,
    )
    check(
        "unusable FASTQ qualities give QUAL '*' with a warning; samtools reads the SAM",
        p.returncode == 0
        and got == {"len": "*", "byte": "*", "good": "I" * 150}
        and "2 FASTQ records had unusable qualities" in p.stderr
        and v.returncode == 0
        and v.stdout.strip() == "3",
        f"qual={ {k: x[:3] for k, x in got.items()} } samtools={v.returncode}",
    )
    print(
        f"\nvalidate_tags: {'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILED'}"
    )
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
