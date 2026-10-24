#!/usr/bin/env python3
"""cli_smoke.py - optik as a user runs it, with every command run end to end on real input forms, outputs checked by standard tools, and the error paths.

  query    build-spss, build-fm and build followed by query on FASTA, multi-line FASTA, FASTQ, gzip and lowercase or N input.
           Counts must equal brute force, --ids must be consistent, and spss-update add or remove followed by query must equal brute force.
  dbruijn  the GFA is validated by bench/validate_unitigs.py (k-mer set, paths, maximality, links), and -s writes the node table
  map      samtools accepts the SAM (quickcheck, flagstat, sort, index, view of a region), FASTA and gzipped FASTQ reads give the same alignments, and -o and stdout output are identical
  errors   missing files, bad k, m or -t, empty input, the wrong -k for an index and an unknown command all give a non-zero exit with a message naming the problem, never a panic or a silent success

Usage: python3 bench/cli_smoke.py [--optik target/release/optik]
"""

import argparse
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from sam_spec import (
    SAM_FLAG_PROPER_PAIR,
    SAM_FLAG_READ1,
    SAM_FLAGS_SECONDARY_SUPPLEMENTARY,
)

FAILS = []
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def check(name, ok, detail=""):
    """Print one check result, with its detail on failure, and remember failures in FAILS."""
    print(
        f"  [{'PASS' if ok else 'FAIL'}] {name}"
        + (f"  {detail}" if detail and not ok else "")
    )
    if not ok:
        FAILS.append(name)


def run(cmd, **kw):
    return subprocess.run(
        [str(c) for c in cmd],
        capture_output=True,
        text=True,
        timeout=600,
        **kw,
        check=False,
    )


def codes(s, k=31):
    """Canonical 2-bit codes of every ACGT k-mer of s, where any other base breaks the run."""
    mask, sh = (1 << 2 * k) - 1, 2 * (k - 1)
    f = r = v = 0
    out = []
    for ch in s.upper():
        c = CODE.get(ch)
        if c is None:
            f = r = v = 0
            continue
        f = ((f << 2) | c) & mask
        r = (r >> 2) | ((3 - c) << sh)
        v += 1
        if v >= k:
            out.append(min(f, r))
    return out


def counts(path, kset, k=31):
    """Brute-force (k-mer windows, positives) of a query file against the k-mer set kset."""
    op = gzip.open if path.endswith(".gz") else open
    total = pos = 0
    with op(path, "rt") as f:
        lines = f.read().split("\n")
    seqs, cur, fq = [], [], lines and lines[0].startswith("@")
    if fq:
        seqs = [
            lines[i + 1]
            for i in range(0, len(lines) - 1, 4)
            if lines[i].startswith("@")
        ]
    else:
        for line in lines:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            elif line:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
    for s in seqs:
        total += max(0, len(s) - k + 1)
        pos += sum(1 for c in codes(s, k) if c in kset)
    return total, pos


def parse_q(text):
    """(num_kmers, num_positive_kmers) of an optik query report, None where absent."""
    k = re.search(r"num_kmers = (\d+)", text)
    p = re.search(r"num_positive_kmers = (\d+)", text)
    return (int(k.group(1)) if k else None, int(p.group(1)) if p else None)


def main():
    """Run every smoke check against the optik binary and exit non-zero on any failure."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    a = ap.parse_args()
    optik = a.optik
    rng = random.Random(5)
    td = tempfile.mkdtemp(prefix="optik_cli_")
    try:
        genome = "".join(rng.choice("ACGT") for _ in range(200000))
        genome2 = "".join(rng.choice("ACGT") for _ in range(50000))
        ref = os.path.join(td, "ref.fa")
        with open(ref, "w") as f:  # multi-line, soft-masked, with two contigs
            f.write(">chrA some description\n")
            for i in range(0, len(genome), 60):
                chunk = genome[i : i + 60]
                f.write((chunk.lower() if 5000 <= i < 8000 else chunk) + "\n")
            f.write(">chrB\n" + genome2 + "\n")
        kset = set(codes(genome)) | set(codes(genome2))
        reads = []
        for i in range(2000):
            src = genome if i % 4 else genome2
            p = rng.randrange(0, len(src) - 150)
            r = list(src[p : p + 150])
            for j in range(150):
                if rng.random() < 0.01:
                    r[j] = rng.choice("ACGT")
            if i % 50 == 0:
                r[rng.randrange(150)] = "N"
            r = "".join(r)
            reads.append((f"r{i}", r.translate(COMP)[::-1] if i % 2 else r))
        rfa, rfq, rgz = (
            os.path.join(td, x) for x in ("reads.fa", "reads.fq", "reads.fq.gz")
        )
        with open(rfa, "w") as f:
            f.writelines(f">{n}\n{s}\n" for n, s in reads)
        with open(rfq, "w") as f:
            f.writelines(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in reads)
        with open(rfq, "rb") as fi, gzip.open(rgz, "wb") as fo:
            fo.write(fi.read())
        truth = counts(rfa, kset)

        print("== query")
        for cmd, idx in (("build-spss", "i.spss"), ("build-fm", "i.fm")):
            ip = os.path.join(td, idx)
            extra = ["-c", "1"] if cmd == "build-fm" else []
            p = run([optik, cmd, "-i", ref, "-k", "31", "-o", ip] + extra)
            check(f"{cmd} exits 0", p.returncode == 0, p.stderr[-200:])
            for q in (rfa, rfq, rgz):
                p = run([optik, "query", "-i", ip, "-q", q, "-k", "31"])
                got = parse_q(p.stdout + p.stderr)
                check(
                    f"{cmd} query {os.path.basename(q)} == brute force",
                    got == truth,
                    f"{got} vs {truth}",
                )
        ip = os.path.join(td, "i.spss")
        # `build` is an alias of build-spss that writes the same index file, and the removed index's -f flag is an error rather than being ignored.
        p = run(
            [
                optik,
                "build",
                "-i",
                ref,
                "-k",
                "31",
                "-o",
                os.path.join(td, "alias.spss"),
            ]
        )
        same = (
            p.returncode == 0
            and Path(os.path.join(td, "alias.spss")).read_bytes()
            == Path(ip).read_bytes()
        )
        check(
            "build is an alias of build-spss (identical index)", same, p.stderr[-200:]
        )
        p = run(
            [
                optik,
                "build-spss",
                "-i",
                ref,
                "-k",
                "31",
                "-c",
                "2",
                "-o",
                os.path.join(td, "c.spss"),
            ]
        )
        check(
            "build-spss -c without --compact is an error (it would be ignored silently)",
            p.returncode != 0 and "--compact" in p.stderr,
            p.stderr[-200:],
        )
        p = run(
            [
                optik,
                "build",
                "-i",
                ref,
                "-k",
                "31",
                "-f",
                "0.5",
                "-o",
                os.path.join(td, "f.spss"),
            ]
        )
        check(
            "build -f, a flag build-spss does not have, is an error",
            p.returncode != 0 and "-f" in p.stderr,
            p.stderr[-200:],
        )
        # The minimizer length the input picks by default, and exactness at every m.
        # m only trades size against lookup work, so the default index and indexes at other m give identical answers.
        g100 = os.path.join(td, "m100k.fa")
        Path(g100).write_text(
            ">s\n" + "".join(rng.choice("ACGT") for _ in range(100_000)) + "\n"
        )
        q100 = os.path.join(td, "m100k_q.fa")
        with open(q100, "w") as f:
            src = Path(g100).read_text().split("\n")[1]
            for qi in range(300):
                a = rng.randrange(len(src) - 60)
                f.write(
                    f">q{qi}\n{src[a : a + 60] if qi % 3 else ''.join(rng.choice('ACGT') for _ in range(60))}\n"
                )
        answers = {}
        for m in ("0", "8", "10", "12", "14"):
            mp = os.path.join(td, f"m{m}.spss")
            p = run([optik, "build-spss", "-i", g100, "-k", "21", "-m", m, "-o", mp])
            got_m = next(
                (
                    text_line.split("\t")[1]
                    for text_line in p.stdout.splitlines()
                    if text_line.startswith("minimizer_len\t")
                ),
                None,
            )
            if m == "0":
                check(
                    "build-spss default m at k=21 within [6, 19]",
                    got_m is not None and 6 <= int(got_m) <= 19,
                    f"got {got_m}",
                )
            else:
                check(f"build-spss -m {m} is used as given", got_m == m, f"got {got_m}")
            pq = run([optik, "query", "-i", mp, "-q", q100, "-k", "21"])
            answers[m] = parse_q(pq.stdout + pq.stderr)
        check(
            "query answers identical at every m (default, 8, 10, 12, 14)",
            len({str(v) for v in answers.values()}) == 1 and answers["0"],
            str(answers),
        )
        ids = os.path.join(td, "ids.txt")
        p = run([optik, "query", "-i", ip, "-q", rfa, "-k", "31", "--ids", ids])
        ok = p.returncode == 0
        seen = {}
        if ok:
            for (_n, s), line in zip(
                reads, Path(ids).read_text().split("\n"), strict=False
            ):
                vals = list(map(int, line.split()))
                ok &= len(vals) == len(s) - 30
                for i, v in enumerate(vals):
                    cc = codes(s[i : i + 31])
                    present = bool(cc) and cc[0] in kset
                    ok &= present == (v >= 0) and (
                        v < 0 or seen.setdefault(cc[0], v) == v
                    )
        check("query --ids: -1 exactly for absent k-mers, one id per k-mer", ok)
        # --ids writes one line per record, and empty records keep their empty line on every
        # input path, so plain and gzipped queries line up.
        qe = os.path.join(td, "ids_empty.fa")
        with open(qe, "w") as f:
            for i, (n, s_) in enumerate(reads[:20]):
                f.write(f">{n}\n{s_}\n")
                if i % 5 == 2:
                    f.write(f">empty{i}\n\n")
        qez = qe + ".gz"
        with open(qe, "rb") as fi, gzip.open(qez, "wb") as fo:
            fo.write(fi.read())
        outs = []
        for q in (qe, qez):
            o_ = os.path.join(td, os.path.basename(q) + ".ids")
            pp = run([optik, "query", "-i", ip, "-q", q, "--ids", o_])
            outs.append(Path(o_).read_text() if pp.returncode == 0 else None)
        check(
            "query --ids: plain and gzipped query with empty records give identical lines",
            outs[0] is not None and outs[0] == outs[1] and outs[0].count("\n") == 24,
            f"lines {[o.count(chr(10)) if o else None for o in outs]}",
        )
        # FM index integrity.
        # A flipped payload bit, a truncated file, a trailing byte and a file with another magic, such as the header of an older optik FM format, are load errors.
        fm = os.path.join(td, "i.fm")
        raw = Path(fm).read_bytes()
        base = parse_q(run([optik, "query", "-i", fm, "-q", rfa]).stdout)
        bad_cases = {
            "flipped bit": raw[: len(raw) // 2]
            + bytes([raw[len(raw) // 2] ^ 0x10])
            + raw[len(raw) // 2 + 1 :],
            "truncated": raw[:-20],
            "trailing byte": raw + b"\0",
            "different magic": b"OKFX" + raw[4:],
            "version byte after the magic (older format)": raw[:4] + b"\x02" + raw[4:],
        }
        for what, blob in bad_cases.items():
            bp = os.path.join(td, "bad.fm")
            Path(bp).write_bytes(blob)
            pp = run([optik, "query", "-i", bp, "-q", rfa])
            check(
                f"FM index with a {what} is rejected",
                pp.returncode != 0,
                pp.stdout[-120:],
            )
        check("FM index answers queries", bool(base), "")
        # An SPSS file under any other magic, an older optik header included, is refused with the advice to rebuild.
        old = os.path.join(td, "old.spss")
        with open(ip, "rb") as fh:
            Path(old).write_bytes(b"OSP7\0\0\0\0" + fh.read()[8:])
        pp = run([optik, "query", "-i", old, "-q", rfa])
        check(
            "SPSS index with another magic is rejected with a rebuild hint",
            pp.returncode != 0 and "rebuild" in pp.stderr,
            pp.stderr[-150:],
        )
        # spss-update removes chrB and adds a new sequence, and query must equal brute force.
        extra_seq = "".join(rng.choice("ACGT") for _ in range(30000))
        addf, remf = os.path.join(td, "add.fa"), os.path.join(td, "rem.fa")
        Path(addf).write_text(">new\n" + extra_seq + "\n")
        Path(remf).write_text(">chrB\n" + genome2 + "\n")
        up = os.path.join(td, "u.spss")
        p = run(
            [optik, "spss-update", "-i", ip, "--add", addf, "--remove", remf, "-o", up]
        )
        check("spss-update exits 0", p.returncode == 0, p.stderr[-200:])
        # Removal deletes every k-mer of chrB, including any shared with chrA.
        kset2 = (set(codes(genome)) - set(codes(genome2))) | set(codes(extra_seq))
        p = run([optik, "query", "-i", up, "-q", rfa, "-k", "31"])
        got, want = parse_q(p.stdout + p.stderr), counts(rfa, kset2)
        check(
            "spss-update result queries == brute force", got == want, f"{got} vs {want}"
        )
        p = run([optik, "query", "-i", ip, "-q", rfa, "-k", "21"])
        check(
            "query with the wrong -k is an error",
            p.returncode != 0 and "k=" in p.stderr,
            p.stderr[-150:],
        )
        # Without -k the index's own k is used (31 here, with the -k 31 run above as the reference).
        p_nok = run([optik, "query", "-i", ip, "-q", rfa])
        p_k = run([optik, "query", "-i", ip, "-q", rfa, "-k", "31"])
        check(
            "query without -k uses the index's k",
            p_nok.returncode == 0
            and parse_q(p_nok.stdout + p_nok.stderr)
            == parse_q(p_k.stdout + p_k.stderr),
            p_nok.stderr[-150:],
        )

        print("== dbruijn")
        gfa = os.path.join(td, "g.gfa")
        p = run(
            [
                optik,
                "dbruijn",
                "-i",
                ref,
                "-k",
                "31",
                "-c",
                "1",
                "-o",
                gfa,
                "-s",
                os.path.join(td, "g.kdbg"),
            ]
        )
        check("dbruijn exits 0", p.returncode == 0, p.stderr[-200:])
        v = run(
            [
                sys.executable,
                os.path.join(HERE, "validate_unitigs.py"),
                "--gfa",
                gfa,
                "--input",
                ref,
                "-k",
                "31",
                "-c",
                "1",
                "--kdbg",
                os.path.join(td, "g.kdbg"),
            ]
        )
        res = json.loads(v.stdout) if v.stdout.strip().startswith("{") else {}
        check(
            "GFA passes validate_unitigs (k-mers, paths, maximality, links, tags)",
            v.returncode == 0 and bool(res),
            (v.stdout + v.stderr)[-300:],
        )
        check(
            "-s node table == oracle (keys, abundances, edge masks)",
            res.get("kdbg_exact", {}).get("pass") is True,
            str(res.get("kdbg_exact"))[:200],
        )
        # Link lookup paths, segment ids stored in the graph's node words or in a separate entry index.
        # The GFA must be byte-identical, and the in-graph one must pass the oracle.
        # Inputs with many segments per node, from dense graphs at small k or short reads, are where "graph" is the default, and even k adds palindromic k-mers.
        import random as _r

        rr = _r.Random(5)
        dense = os.path.join(td, "dense.fa")
        with open(dense, "w") as f:
            f.write(">g\n" + "".join(rr.choice("ACGT") for _ in range(60000)) + "\n")
            for i in range(400):
                f.write(
                    f">r{i}\n" + "".join(rr.choice("ACGT") for _ in range(36)) + "\n"
                )
        for kk in ("9", "10", "15", "31"):
            outs = {}
            for mode in ("graph", "index"):
                gp = os.path.join(td, f"links_{mode}_{kk}.gfa")
                # the table's two link paths, which partitioned builds use neither of
                pm = run(
                    [
                        optik,
                        "dbruijn",
                        "-i",
                        dense,
                        "-k",
                        kk,
                        "-c",
                        "1",
                        "-t",
                        "3",
                        "-o",
                        gp,
                        "--graph",
                        "table",
                    ],
                    env=dict(os.environ, OPTIK_DBG_LINKS=mode),
                )
                outs[mode] = Path(gp).read_bytes() if pm.returncode == 0 else None
            check(
                f"dbruijn k={kk}: in-graph segment ids == entry index (byte-identical GFA)",
                outs["graph"] is not None and outs["graph"] == outs["index"],
            )
            v = run(
                [
                    sys.executable,
                    os.path.join(HERE, "validate_unitigs.py"),
                    "--gfa",
                    os.path.join(td, f"links_graph_{kk}.gfa"),
                    "--input",
                    dense,
                    "-k",
                    kk,
                    "-c",
                    "1",
                ]
            )
            check(
                f"dbruijn k={kk}: in-graph links pass validate_unitigs",
                v.returncode == 0,
                (v.stdout + v.stderr)[-200:],
            )

        print("== map")
        samtools = shutil.which("samtools")
        sams = {}
        for name, q in (("fa", rfa), ("fqgz", rgz)):
            out = os.path.join(td, f"m_{name}.sam")
            p = run([optik, "map", "-r", ref, "-q", q, "-t", "4", "-o", out])
            check(f"map {name} exits 0", p.returncode == 0, p.stderr[-200:])
            sams[name] = out

        def body(f):
            """Lines of a SAM file without its @PG header lines."""
            with open(f) as fh:
                return [
                    text_line for text_line in fh if not text_line.startswith("@PG")
                ]

        check(
            "FASTA and gzipped FASTQ reads give the same alignments",
            [
                text_line.split("\t")[:9]
                for text_line in body(sams["fa"])
                if not text_line.startswith("@")
            ]
            == [
                text_line.split("\t")[:9]
                for text_line in body(sams["fqgz"])
                if not text_line.startswith("@")
            ],
        )
        p = run([optik, "map", "-r", ref, "-q", rfa, "-t", "4"])
        check(
            "map to stdout == map -o",
            [
                text_line
                for text_line in p.stdout.splitlines(True)
                if not text_line.startswith("@PG")
            ]
            == body(sams["fa"]),
        )
        if samtools:
            s = sams["fa"]
            check(
                "samtools quickcheck", run([samtools, "quickcheck", s]).returncode == 0
            )
            fs = run([samtools, "flagstat", s])
            m = re.search(r"(\d+) \+ \d+ primary mapped", fs.stdout)
            check(
                "samtools flagstat: every primary mapped",
                m is not None and int(m.group(1)) == len(reads),
                fs.stdout[:200],
            )
            bam = os.path.join(td, "m.bam")
            ok = (
                run([samtools, "sort", "-o", bam, s]).returncode == 0
                and run([samtools, "index", bam]).returncode == 0
            )
            check("samtools sort + index", ok)
            reg = run([samtools, "view", "-c", bam, "chrB"])
            check(
                "samtools view of a region",
                reg.returncode == 0 and int(reg.stdout.strip() or 0) > 0,
                reg.stderr[-150:],
            )
            calmd = run(
                [samtools, "calmd", "-e", s, ref]
            )  # recomputes NM and MD from the reference
            bad = 0
            with open(s) as fh:
                orig = {
                    (text_line.split("\t")[0], text_line.split("\t")[1]): text_line
                    for text_line in fh
                    if not text_line.startswith("@")
                }
            for text_line in calmd.stdout.splitlines():
                if text_line.startswith("@"):
                    continue
                t = text_line.split("\t")
                nm = [x for x in t[11:] if x.startswith("NM:i:")]
                o = orig.get((t[0], t[1]))
                onm = (
                    [
                        x
                        for x in o.rstrip("\n").split("\t")[11:]
                        if x.startswith("NM:i:")
                    ]
                    if o
                    else []
                )
                bad += bool(nm and onm and nm[0] != onm[0])
            check(
                "samtools calmd agrees with every NM:i",
                calmd.returncode == 0 and bad == 0,
                f"{bad} differ",
            )
        else:
            check("samtools available for SAM checks", False, "samtools not on PATH")

        print("== several inputs, read groups")
        # Paired-end style, where -i R1 -i R2 of FASTA and gzipped FASTQ equals one joined file
        half = len(reads) // 2
        r1, r2gz, cat = (os.path.join(td, x) for x in ("r1.fa", "r2.fq.gz", "r12.fa"))
        with open(r1, "w") as f:
            for n_, s_ in reads[:half]:
                f.write(f">{n_}\n{s_}\n")
        with gzip.open(r2gz, "wt") as f:
            for n_, s_ in reads[half:]:
                f.write(f"@{n_}\n{s_}\n+\n{'I' * len(s_)}\n")
        shutil.copy(rfa, cat)
        for cmd in (
            ["dbruijn", "-c", "1"],
            ["build-fm"],
            ["build-spss", "--compact"],
            ["build-spss"],
        ):
            one, two = os.path.join(td, "one.out"), os.path.join(td, "two.out")
            p1 = run([optik] + cmd + ["-i", cat, "-k", "31", "-o", one])
            p2 = run([optik] + cmd + ["-i", r1, "-i", r2gz, "-k", "31", "-o", two])
            same = (
                p1.returncode == 0
                and p2.returncode == 0
                and Path(one).read_bytes() == Path(two).read_bytes()
            )
            check(
                f"{' '.join(cmd)} -i R1 -i R2.gz == the concatenated input",
                same,
                p2.stderr.strip()[-160:],
            )
        rgsam = os.path.join(td, "rg.sam")
        p = run(
            [
                optik,
                "map",
                "-r",
                ref,
                "-q",
                rfq,
                "-R",
                "@RG\\tID:s1\\tSM:x\\tPL:ILLUMINA",
                "-o",
                rgsam,
            ]
        )
        lines = Path(rgsam).read_text().splitlines() if p.returncode == 0 else []
        body = [text_line for text_line in lines if not text_line.startswith("@")]
        check(
            "map -R: @RG header line and RG:Z:s1 on every record",
            "@RG\tID:s1\tSM:x\tPL:ILLUMINA" in lines
            and body
            and all(text_line.endswith("\tRG:Z:s1") for text_line in body),
            p.stderr.strip()[-160:],
        )
        check(
            "map @PG records the command line (CL:)",
            any(
                text_line.startswith("@PG\tID:optik")
                and "\tCL:" in text_line
                and " map " in text_line
                for text_line in lines
            ),
        )

        print("== paired-end")
        pre = os.path.join(td, "pe")
        subprocess.run(
            [
                sys.executable,
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "gen_pairs.py"
                ),
                "-r",
                ref,
                "-n",
                "2000",
                "--mean",
                "300",
                "--sd",
                "30",
                "-s",
                "3",
                "-o",
                pre,
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        pe1, pe2, pesam = pre + "_1.fq", pre + "_2.fq", os.path.join(td, "pe.sam")
        p = run([optik, "map", "-r", ref, "-q", pe1, "-q", pe2, "-o", pesam])
        body = []
        if p.returncode == 0:
            with open(pesam) as fh:
                body = [x.split("\t") for x in fh if not x.startswith("@")]
        prim = [f for f in body if not int(f[1]) & SAM_FLAGS_SECONDARY_SUPPLEMENTARY]
        n_prop = sum(1 for f in prim if int(f[1]) & SAM_FLAG_PROPER_PAIR)
        check(
            "map -q R1 -q R2: every record paired (0x1, one 0x40 + one 0x80 per name), >= 95% proper",
            len(prim) == 4000
            and all(int(f[1]) & 1 for f in body)
            and sum(1 for f in prim if int(f[1]) & SAM_FLAG_READ1) == 2000
            and n_prop >= 0.95 * 4000,
            f"records={len(prim)} proper={n_prop} {p.stderr.strip()[-120:]}",
        )
        sw = shutil.which("samtools")
        if sw:
            q = subprocess.run(
                [sw, "quickcheck", pesam], capture_output=True, check=False
            )
            fs = subprocess.run(
                [sw, "flagstat", pesam], capture_output=True, text=True, check=False
            ).stdout
            check(
                "paired SAM: samtools quickcheck, flagstat counts read1/read2",
                q.returncode == 0 and "2000 + 0 read1" in fs and "2000 + 0 read2" in fs,
                fs.replace("\n", " | ")[:200],
            )
        # R2 with a record missing or out of order, or a third -q
        lines2 = Path(pe2).read_text().splitlines(True)
        short2, swap2 = os.path.join(td, "short_2.fq"), os.path.join(td, "swap_2.fq")
        Path(short2).write_text("".join(lines2[:-4]))
        Path(swap2).write_text("".join(lines2[4:8] + lines2[0:4] + lines2[8:]))
        for nm, args, want in (
            ("R2 shorter than R1", ["-q", pe1, "-q", short2], "fewer records"),
            ("mates out of order", ["-q", pe1, "-q", swap2], "mate names differ"),
            ("three read files", ["-q", pe1, "-q", pe2, "-q", pe1], "one -q"),
        ):
            p = run(
                [optik, "map", "-r", ref] + args + ["-o", os.path.join(td, "x.sam")]
            )
            check(
                f"paired error: {nm}",
                p.returncode == 1 and want in p.stderr,
                p.stderr.strip()[-160:],
            )
        p = run(
            [
                optik,
                "map",
                "-r",
                ref,
                "-q",
                short2,
                "-q",
                pe2,
                "-o",
                os.path.join(td, "x.sam"),
            ]
        )
        check(
            "paired error: R2 longer than R1",
            p.returncode == 1 and "more records" in p.stderr,
            p.stderr.strip()[-160:],
        )

        print("== --graph")
        outs = {}
        for mode in ("auto", "part", "table"):
            for kk in ("15", "31"):
                o_ = os.path.join(td, f"g_{mode}_{kk}.gfa")
                p = run(
                    [
                        optik,
                        "dbruijn",
                        "-i",
                        rfa,
                        "-k",
                        kk,
                        "-c",
                        "1",
                        "--graph",
                        mode,
                        "-o",
                        o_,
                    ]
                )
                outs[(mode, kk)] = Path(o_).read_bytes() if p.returncode == 0 else None
        check(
            "dbruijn --graph auto|part|table: the same GFA at k = 15 and 31 (part below the default 19 too)",
            all(
                outs[(m, kk)] is not None and outs[(m, kk)] == outs[("auto", kk)]
                for m in ("part", "table")
                for kk in ("15", "31")
            ),
        )
        p = run(
            [
                optik,
                "dbruijn",
                "-i",
                rfa,
                "--graph",
                "part",
                "-s",
                os.path.join(td, "x.kdbg"),
                "-o",
                os.path.join(td, "x.gfa"),
            ]
        )
        check(
            "--graph part with -s (no node table to archive): error",
            p.returncode == 1 and "-s" in p.stderr,
            p.stderr.strip()[-120:],
        )
        p = run(
            [
                optik,
                "dbruijn",
                "-i",
                rfa,
                "--graph",
                "fast",
                "-o",
                os.path.join(td, "x.gfa"),
            ]
        )
        check(
            "--graph with an unknown mode: error",
            p.returncode == 1 and "--graph" in p.stderr,
            p.stderr.strip()[-120:],
        )

        print("== full disk, closed pipes, stdin")
        spss_idx = os.path.join(td, "i.spss")
        # Every writer must fail loudly on a full disk (/dev/full, ENOSPC), index saves included, since the FM save used to return success.
        for cmd in (
            [optik, "build-fm", "-i", ref, "-o", "/dev/full"],
            [optik, "build-spss", "-i", ref, "-o", "/dev/full"],
            [optik, "dbruijn", "-i", ref, "-c", "1", "-o", "/dev/full"],
            [
                optik,
                "dbruijn",
                "-i",
                ref,
                "-c",
                "1",
                "-o",
                os.path.join(td, "x.gfa"),
                "-s",
                "/dev/full",
            ],
            [optik, "map", "-r", ref, "-q", rfa, "-o", "/dev/full"],
            [optik, "query", "-i", spss_idx, "-q", rfa, "--ids", "/dev/full"],
        ):
            p = run(cmd)
            check(
                f"full disk: {cmd[1]} {' '.join(cmd[-2:])} -> error exit 1",
                p.returncode == 1
                and "error" in p.stderr
                and "panicked" not in p.stderr,
                f"exit {p.returncode}: {p.stderr.strip()[-120:]}",
            )
        with open("/dev/full", "w") as full:
            p = subprocess.run(
                [optik, "query", "-i", spss_idx, "-q", rfa],
                stdout=full,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        check(
            "stdout on a full disk: error exit 1, no core dump",
            p.returncode == 1 and "writing to stdout failed" in p.stderr,
            f"exit {p.returncode}",
        )
        # A pipe closed by | head ends quietly, as C tools do.
        p = subprocess.run(
            f"'{optik}' map -r '{ref}' -q '{rfa}' | head -1",
            shell=True,
            capture_output=True,
            text=True,
            check=False,
        )
        check(
            "map | head -1: no panic, no core dump",
            p.stdout.startswith("@HD")
            and "panicked" not in p.stderr
            and "core" not in p.stderr,
            p.stderr.strip()[-120:],
        )
        # --ids with the queries on stdin reads them in one pass.  A second read of stdin used to
        # be empty, and the report then said 0 k-mers.
        with open(rfa) as fi:
            p = subprocess.run(
                [
                    optik,
                    "query",
                    "-i",
                    spss_idx,
                    "-q",
                    "-",
                    "--ids",
                    os.path.join(td, "ids.txt"),
                ],
                stdin=fi,
                capture_output=True,
                text=True,
                check=False,
            )
        check(
            "query -q - --ids: counts from the same pass",
            parse_q(p.stdout + p.stderr) == truth,
            f"{parse_q(p.stdout + p.stderr)} vs {truth}",
        )
        p = run([optik, "dbruijn", "-i", ref, "-o", os.path.join(td, "g2.gfa")])
        check(
            "dbruijn on a genome with the reads default -c 2: warns to use -c 1",
            p.returncode == 0 and "use -c 1" in p.stderr,
            p.stderr.strip()[-160:],
        )
        deep = os.path.join(
            td, "deep.fa"
        )  # the read set four times over, about 5x coverage
        with open(deep, "w") as f:
            f.write(Path(rfa).read_text() * 4)
        p = run(
            [optik, "dbruijn", "-i", deep, "-k", "21", "-o", os.path.join(td, "g3.gfa")]
        )
        check(
            "dbruijn on 5x reads with -c 2: no such warning",
            p.returncode == 0 and "use -c 1" not in p.stderr,
            p.stderr.strip()[-160:],
        )
        p = run(
            [
                optik,
                "dbruijn",
                "-i",
                rfa,
                "-o",
                os.path.join(td, "same.gfa"),
                "-s",
                os.path.join(td, "same.gfa"),
            ]
        )
        check(
            "two outputs on one file: error",
            p.returncode == 1 and "same file" in p.stderr,
            p.stderr[-120:],
        )

        print("== errors")
        empty = os.path.join(td, "empty.fa")
        Path(empty).write_text(">e\n\n")
        cases = [
            (
                "missing reference",
                [optik, "map", "-r", "/nonexistent.fa", "-q", rfa],
                "cannot open",
            ),
            (
                "missing reads",
                [optik, "map", "-r", ref, "-q", "/nonexistent.fa"],
                "cannot open",
            ),
            (
                "missing index",
                [optik, "query", "-i", "/nonexistent.idx", "-q", rfa],
                "cannot open",
            ),
            (
                "k out of range",
                [
                    optik,
                    "build-spss",
                    "-i",
                    ref,
                    "-k",
                    "40",
                    "-o",
                    os.path.join(td, "x"),
                ],
                "k=40",
            ),
            (
                "m >= k",
                [
                    optik,
                    "build-spss",
                    "-i",
                    ref,
                    "-k",
                    "31",
                    "-m",
                    "31",
                    "-o",
                    os.path.join(td, "x"),
                ],
                "m",
            ),
            ("bad -t", [optik, "map", "-r", ref, "-q", rfa, "-t", "abc"], "-t"),
            (
                "empty input index",
                [optik, "build-spss", "-i", empty, "-o", os.path.join(td, "x")],
                "no k-mers",
            ),
            (
                "empty dbruijn input",
                [optik, "dbruijn", "-i", empty, "-o", os.path.join(td, "x.gfa")],
                "k-mers",
            ),
            ("unknown command", [optik, "frobnicate"], "USAGE"),
            (
                "flag given twice",
                [optik, "map", "-r", ref, "-r", ref, "-q", rfa],
                "more than once",
            ),
            (
                "-k given twice",
                [optik, "dbruijn", "-i", rfa, "-k", "21", "-k", "31"],
                "more than once",
            ),
            ("stdin twice", [optik, "dbruijn", "-i", "-", "-i", "-"], "only once"),
            (
                "-R without ID",
                [optik, "map", "-r", ref, "-q", rfa, "-R", "@RG\\tSM:x"],
                "ID",
            ),
            (
                "-R not a read group",
                [optik, "map", "-r", ref, "-q", rfa, "-R", "ID:x"],
                "@RG",
            ),
            ("no arguments without a terminal", [optik], "USAGE"),
        ]
        for name, cmd, want in cases:
            p = run(cmd, stdin=subprocess.DEVNULL)
            msg = p.stderr + p.stdout
            check(
                f"error: {name}",
                p.returncode != 0 and want in msg and "panicked" not in msg,
                f"exit {p.returncode}: {msg.strip()[-160:]}",
            )
    finally:
        shutil.rmtree(td, ignore_errors=True)
    print(f"cli_smoke: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
