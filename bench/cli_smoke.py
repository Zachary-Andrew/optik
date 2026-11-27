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
from dataclasses import dataclass
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from checklog import CheckLog
from diagnostics import DETAIL_CHARS
from dna import revcomp
from fastq import FASTQ_LINES_PER_RECORD
from sam_spec import (
    SamCol,
    SamFlag,
)

# seconds before one optik run counts as hung
RUN_TIMEOUT_S = 600
# k of the brute-force k-mer oracles and of the smoke indexes
SMOKE_K = 31
# share of pairs that must be flagged proper (0x2)
MIN_PROPER_SHARE = 0.95

check = CheckLog(detail_on_pass=False)
FAILS = check.failed
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def run(cmd, **kwargs):
    return subprocess.run(
        [str(arg) for arg in cmd],
        capture_output=True,
        text=True,
        timeout=RUN_TIMEOUT_S,
        **kwargs,
        check=False,
    )


def codes(seq, kmer_len=SMOKE_K):
    """Canonical 2-bit codes of every ACGT k-mer of s, where any other base breaks the run."""
    mask, rc_shift = (1 << 2 * kmer_len) - 1, 2 * (kmer_len - 1)
    fwd = rev = run_len = 0
    out = []
    for char in seq.upper():
        code = CODE.get(char)
        if code is None:
            fwd = rev = run_len = 0
            continue
        fwd = ((fwd << 2) | code) & mask
        rev = (rev >> 2) | ((3 - code) << rc_shift)
        run_len += 1
        if run_len >= kmer_len:
            out.append(min(fwd, rev))
    return out


def counts(path, kset, kmer_len=SMOKE_K):
    """Brute-force (k-mer windows, positives) of a query file against the k-mer set kset."""
    opener = gzip.open if path.endswith(".gz") else open
    total = pos = 0
    with opener(path, "rt") as handle:
        lines = handle.read().split("\n")
    seqs, cur, is_fastq = [], [], lines and lines[0].startswith("@")
    if is_fastq:
        seqs = [
            lines[i + 1]
            for i in range(0, len(lines) - 1, FASTQ_LINES_PER_RECORD)
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
    for seq in seqs:
        total += max(0, len(seq) - kmer_len + 1)
        pos += sum(1 for code in codes(seq, kmer_len) if code in kset)
    return total, pos


def parse_q(text):
    """(num_kmers, num_positive_kmers) of an optik query report, None where absent."""
    kmers_match = re.search(r"num_kmers = (\d+)", text)
    pos_match = re.search(r"num_positive_kmers = (\d+)", text)
    return (int(kmers_match.group(1)) if kmers_match else None, int(pos_match.group(1)) if pos_match else None)


@dataclass
class Fixtures:
    """Holds what the sections share.

    Besides the binary and a scratch directory, that means the seeded RNG, both genomes, the reads and their k-mer truth.
    """

    optik: str
    tmp_dir: str
    rng: random.Random
    genome: str
    genome2: str
    ref: str
    rfa: str
    rfq: str
    rgz: str
    reads: list
    kset: set
    truth: tuple   # (k-mer windows, positives) of the reads, from counts()


def make_fixtures(optik, tmp_dir, rng):
    """Build the random genomes, reads and k-mer truth that the sections check against."""
    genome = "".join(rng.choice("ACGT") for _ in range(200000))
    genome2 = "".join(rng.choice("ACGT") for _ in range(50000))
    ref = os.path.join(tmp_dir, "ref.fa")
    with open(ref, "w") as handle:  # multi-line, soft-masked, with two contigs
        handle.write(">chrA some description\n")
        for i in range(0, len(genome), 60):
            chunk = genome[i : i + 60]
            # A lowercase stretch checks that soft-masked bases count like uppercase ones.
            handle.write((chunk.lower() if 5000 <= i < 8000 else chunk) + "\n")
        handle.write(">chrB\n" + genome2 + "\n")
    kset = set(codes(genome)) | set(codes(genome2))
    reads = []
    for i in range(2000):
        # Reads come from chrB one time in four. Odd-numbered reads are reverse-complemented and every fiftieth read gets an N.
        src = genome if i % 4 else genome2
        start = rng.randrange(0, len(src) - 150)
        read_chars = list(src[start : start + 150])
        for j in range(150):
            if rng.random() < 0.01:
                read_chars[j] = rng.choice("ACGT")
        if i % 50 == 0:
            read_chars[rng.randrange(150)] = "N"
        read_chars = "".join(read_chars)
        reads.append((f"r{i}", revcomp(read_chars) if i % 2 else read_chars))
    rfa, rfq, rgz = (
        os.path.join(tmp_dir, file_name) for file_name in ("reads.fa", "reads.fq", "reads.fq.gz")
    )
    with open(rfa, "w") as handle:
        handle.writelines(f">{rec_name}\n{rec_seq}\n" for rec_name, rec_seq in reads)
    with open(rfq, "w") as handle:
        handle.writelines(f"@{rec_name}\n{rec_seq}\n+\n{'I' * len(rec_seq)}\n" for rec_name, rec_seq in reads)
    with open(rfq, "rb") as in_handle, gzip.open(rgz, "wb") as out_handle:
        out_handle.write(in_handle.read())
    # The truth counts come from Python sets, independent of optik.
    truth = counts(rfa, kset)

    return Fixtures(optik=optik, tmp_dir=tmp_dir, rng=rng, genome=genome, genome2=genome2, ref=ref, rfa=rfa, rfq=rfq, rgz=rgz, reads=reads, kset=kset, truth=truth)


def check_query(fx):
    """Checks query, --ids, index integrity, spss-update and how -k is resolved."""
    optik, tmp_dir, rng, genome, genome2, ref, rfa, rfq, rgz, reads, kset, truth = fx.optik, fx.tmp_dir, fx.rng, fx.genome, fx.genome2, fx.ref, fx.rfa, fx.rfq, fx.rgz, fx.reads, fx.kset, fx.truth
    print("== query")
    for cmd, idx in (("build-spss", "i.spss"), ("build-fm", "i.fm")):
        index_path = os.path.join(tmp_dir, idx)
        extra = ["-c", "1"] if cmd == "build-fm" else []
        proc = run([optik, cmd, "-i", ref, "-k", "31", "-o", index_path] + extra)
        check(f"{cmd} exits 0", proc.returncode == 0, proc.stderr[-DETAIL_CHARS:])
        for query_path in (rfa, rfq, rgz):
            proc = run([optik, "query", "-i", index_path, "-q", query_path, "-k", "31"])
            got = parse_q(proc.stdout + proc.stderr)
            check(
                f"{cmd} query {os.path.basename(query_path)} == brute force",
                got == truth,
                f"{got} vs {truth}",
            )
    index_path = os.path.join(tmp_dir, "i.spss")
    # `build` is an alias of build-spss that writes the same index file, and the removed index's -f flag is an error rather than being ignored.
    proc = run(
        [
            optik,
            "build",
            "-i",
            ref,
            "-k",
            "31",
            "-o",
            os.path.join(tmp_dir, "alias.spss"),
        ]
    )
    # The alias has to write a byte-identical index, not just a working one.
    same = (
        proc.returncode == 0
        and Path(os.path.join(tmp_dir, "alias.spss")).read_bytes()
        == Path(index_path).read_bytes()
    )
    check(
        "build is an alias of build-spss (identical index)", same, proc.stderr[-DETAIL_CHARS:]
    )
    proc = run(
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
            os.path.join(tmp_dir, "c.spss"),
        ]
    )
    check(
        "build-spss -c without --compact is an error (it would be ignored silently)",
        proc.returncode != 0 and "--compact" in proc.stderr,
        proc.stderr[-DETAIL_CHARS:],
    )
    proc = run(
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
            os.path.join(tmp_dir, "f.spss"),
        ]
    )
    check(
        "build -f, a flag build-spss does not have, is an error",
        proc.returncode != 0 and "-f" in proc.stderr,
        proc.stderr[-DETAIL_CHARS:],
    )
    # The minimizer length the input picks by default, and exactness at every m.
    # m only trades size against lookup work, so the default index and indexes at other m give identical answers.
    g100 = os.path.join(tmp_dir, "m100k.fa")
    Path(g100).write_text(
        ">s\n" + "".join(rng.choice("ACGT") for _ in range(100_000)) + "\n"
    )
    q100 = os.path.join(tmp_dir, "m100k_q.fa")
    with open(q100, "w") as handle:
        src = Path(g100).read_text().split("\n")[1]
        for qi in range(300):
            offset = rng.randrange(len(src) - 60)
            handle.write(
                f">q{qi}\n{src[offset : offset + 60] if qi % 3 else ''.join(rng.choice('ACGT') for _ in range(60))}\n"
            )
    # m changes the index layout but never its contents, so every m must give the same query answers.
    answers = {}
    for min_len in ("0", "8", "10", "12", "14"):
        min_index = os.path.join(tmp_dir, f"m{min_len}.spss")
        proc = run([optik, "build-spss", "-i", g100, "-k", "21", "-m", min_len, "-o", min_index])
        got_m = next(
            (
                text_line.removeprefix("minimizer_len\t")
                for text_line in proc.stdout.splitlines()
                if text_line.startswith("minimizer_len\t")
            ),
            None,
        )
        if min_len == "0":
            check(
                "build-spss default m at k=21 within [6, 19]",
                got_m is not None and 6 <= int(got_m) <= 19,
                f"got {got_m}",
            )
        else:
            check(f"build-spss -m {min_len} is used as given", got_m == min_len, f"got {got_m}")
        min_query_proc = run([optik, "query", "-i", min_index, "-q", q100, "-k", "21"])
        answers[min_len] = parse_q(min_query_proc.stdout + min_query_proc.stderr)
    check(
        "query answers identical at every m (default, 8, 10, 12, 14)",
        len({str(answer) for answer in answers.values()}) == 1 and answers["0"],
        str(answers),
    )
    # --ids must give -1 for every absent k-mer and exactly one id to each present one.
    ids = os.path.join(tmp_dir, "ids.txt")
    proc = run([optik, "query", "-i", index_path, "-q", rfa, "-k", "31", "--ids", ids])
    passed = proc.returncode == 0
    seen = {}
    if passed:
        for (_name, read_seq), line in zip(
            reads, Path(ids).read_text().split("\n"), strict=False
        ):
            vals = list(map(int, line.split()))
            passed &= len(vals) == len(read_seq) - SMOKE_K + 1
            for i, kmer_id in enumerate(vals):
                codes_at = codes(read_seq[i : i + SMOKE_K])
                present = bool(codes_at) and codes_at[0] in kset
                passed &= present == (kmer_id >= 0) and (
                    kmer_id < 0 or seen.setdefault(codes_at[0], kmer_id) == kmer_id
                )
    check("query --ids: -1 exactly for absent k-mers, one id per k-mer", passed)
    # --ids writes one line per record, and empty records keep their empty line on every
    # input path, so plain and gzipped queries line up.
    empty_query = os.path.join(tmp_dir, "ids_empty.fa")
    with open(empty_query, "w") as handle:
        for i, (rec_name, rec_seq) in enumerate(reads[:20]):
            handle.write(f">{rec_name}\n{rec_seq}\n")
            if i % 5 == 2:
                handle.write(f">empty{i}\n\n")
    qez = empty_query + ".gz"
    with open(empty_query, "rb") as in_handle, gzip.open(qez, "wb") as out_handle:
        out_handle.write(in_handle.read())
    outs = []
    for query_path in (empty_query, qez):
        ids_out = os.path.join(tmp_dir, os.path.basename(query_path) + ".ids")
        query_proc = run([optik, "query", "-i", index_path, "-q", query_path, "--ids", ids_out])
        outs.append(Path(ids_out).read_text() if query_proc.returncode == 0 else None)
    check(
        "query --ids: plain and gzipped query with empty records give identical lines",
        outs[0] is not None and outs[0] == outs[1] and outs[0].count("\n") == 24,
        f"lines {[out.count(chr(10)) if out else None for out in outs]}",
    )
    # FM index integrity.
    # A flipped payload bit, a truncated file, a trailing byte and a file with another magic, such as the header of an older optik FM format, are load errors.
    fm_path = os.path.join(tmp_dir, "i.fm")
    raw = Path(fm_path).read_bytes()
    base = parse_q(run([optik, "query", "-i", fm_path, "-q", rfa]).stdout)
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
        bad_path = os.path.join(tmp_dir, "bad.fm")
        Path(bad_path).write_bytes(blob)
        query_proc = run([optik, "query", "-i", bad_path, "-q", rfa])
        check(
            f"FM index with a {what} is rejected",
            query_proc.returncode != 0,
            query_proc.stdout[-DETAIL_CHARS:],
        )
    check("FM index answers queries", bool(base), "")
    # An SPSS file under any other magic, an older optik header included, is refused with the advice to rebuild.
    old = os.path.join(tmp_dir, "old.spss")
    with open(index_path, "rb") as read_handle:
        Path(old).write_bytes(b"OSP7\0\0\0\0" + read_handle.read()[8:])
    query_proc = run([optik, "query", "-i", old, "-q", rfa])
    check(
        "SPSS index with another magic is rejected with a rebuild hint",
        query_proc.returncode != 0 and "rebuild" in query_proc.stderr,
        query_proc.stderr[-DETAIL_CHARS:],
    )
    # spss-update removes chrB and adds a new sequence, and query must equal brute force.
    extra_seq = "".join(rng.choice("ACGT") for _ in range(30000))
    addf, remf = os.path.join(tmp_dir, "add.fa"), os.path.join(tmp_dir, "rem.fa")
    Path(addf).write_text(">new\n" + extra_seq + "\n")
    Path(remf).write_text(">chrB\n" + genome2 + "\n")
    updated_idx = os.path.join(tmp_dir, "u.spss")
    proc = run(
        [optik, "spss-update", "-i", index_path, "--add", addf, "--remove", remf, "-o", updated_idx]
    )
    check("spss-update exits 0", proc.returncode == 0, proc.stderr[-DETAIL_CHARS:])
    # Removal deletes every k-mer of chrB, including any shared with chrA.
    kset2 = (set(codes(genome)) - set(codes(genome2))) | set(codes(extra_seq))
    proc = run([optik, "query", "-i", updated_idx, "-q", rfa, "-k", "31"])
    got, want = parse_q(proc.stdout + proc.stderr), counts(rfa, kset2)
    check(
        "spss-update result queries == brute force", got == want, f"{got} vs {want}"
    )
    proc = run([optik, "query", "-i", index_path, "-q", rfa, "-k", "21"])
    check(
        "query with the wrong -k is an error",
        proc.returncode != 0 and "k=" in proc.stderr,
        proc.stderr[-DETAIL_CHARS:],
    )
    # Without -k the index's own k is used (31 here, with the -k 31 run above as the reference).
    p_nok = run([optik, "query", "-i", index_path, "-q", rfa])
    p_k = run([optik, "query", "-i", index_path, "-q", rfa, "-k", "31"])
    check(
        "query without -k uses the index's k",
        p_nok.returncode == 0
        and parse_q(p_nok.stdout + p_nok.stderr)
        == parse_q(p_k.stdout + p_k.stderr),
        p_nok.stderr[-DETAIL_CHARS:],
    )


def check_dbruijn(fx):
    """Compares dbruijn output with the oracle and confirms both link-lookup paths agree."""
    optik, tmp_dir, ref = fx.optik, fx.tmp_dir, fx.ref
    print("== dbruijn")
    gfa = os.path.join(tmp_dir, "g.gfa")
    proc = run(
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
            os.path.join(tmp_dir, "g.kdbg"),
        ]
    )
    # The independent validator checks the GFA's k-mers, paths, maximality, links and tags, and the -s archive against its oracle.
    check("dbruijn exits 0", proc.returncode == 0, proc.stderr[-DETAIL_CHARS:])
    validate_proc = run(
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
            os.path.join(tmp_dir, "g.kdbg"),
        ]
    )
    res = json.loads(validate_proc.stdout) if validate_proc.stdout.strip().startswith("{") else {}
    check(
        "GFA passes validate_unitigs (k-mers, paths, maximality, links, tags)",
        validate_proc.returncode == 0 and bool(res),
        (validate_proc.stdout + validate_proc.stderr)[-DETAIL_CHARS:],
    )
    check(
        "-s node table == oracle (keys, abundances, edge masks)",
        res.get("kdbg_exact", {}).get("pass") is True,
        str(res.get("kdbg_exact"))[:DETAIL_CHARS],
    )
    # Link lookup paths, segment ids stored in the graph's node words or in a separate entry index.
    # The GFA must be byte-identical, and the in-graph one must pass the oracle.
    # Inputs with many segments per node, from dense graphs at small k or short reads, are where "graph" is the default, and even k adds palindromic k-mers.
    import random as random_mod

    link_rng = random_mod.Random(5)
    dense = os.path.join(tmp_dir, "dense.fa")
    with open(dense, "w") as handle:
        handle.write(">g\n" + "".join(link_rng.choice("ACGT") for _ in range(60000)) + "\n")
        handle.writelines(f">r{i}\n" + "".join(link_rng.choice("ACGT") for _ in range(36)) + "\n" for i in range(400))
    # The table builds its links in two ways, from the graph or from an index of the ends, and both must give the same GFA at each k.
    for kmer_arg in ("9", "10", "15", "31"):
        outs = {}
        for mode in ("graph", "index"):
            gfa_path = os.path.join(tmp_dir, f"links_{mode}_{kmer_arg}.gfa")
            # the table's two link paths, which partitioned builds use neither of
            build_proc = run(
                [
                    optik,
                    "dbruijn",
                    "-i",
                    dense,
                    "-k",
                    kmer_arg,
                    "-c",
                    "1",
                    "-t",
                    "3",
                    "-o",
                    gfa_path,
                    "--graph",
                    "table",
                ],
                env=dict(os.environ, OPTIK_DBG_LINKS=mode),
            )
            outs[mode] = Path(gfa_path).read_bytes() if build_proc.returncode == 0 else None
        check(
            f"dbruijn k={kmer_arg}: in-graph segment ids == entry index (byte-identical GFA)",
            outs["graph"] is not None and outs["graph"] == outs["index"],
        )
        validate_proc = run(
            [
                sys.executable,
                os.path.join(HERE, "validate_unitigs.py"),
                "--gfa",
                os.path.join(tmp_dir, f"links_graph_{kmer_arg}.gfa"),
                "--input",
                dense,
                "-k",
                kmer_arg,
                "-c",
                "1",
            ]
        )
        check(
            f"dbruijn k={kmer_arg}: in-graph links pass validate_unitigs",
            validate_proc.returncode == 0,
            (validate_proc.stdout + validate_proc.stderr)[-DETAIL_CHARS:],
        )


def check_map(fx):
    """Checks that map writes valid SAM, gives the same output at any thread count and accepts every read-file form."""
    optik, tmp_dir, ref, rfa, rgz, reads = fx.optik, fx.tmp_dir, fx.ref, fx.rfa, fx.rgz, fx.reads
    print("== map")
    samtools = shutil.which("samtools")
    sams = {}
    for name, query_path in (("fa", rfa), ("fqgz", rgz)):
        out = os.path.join(tmp_dir, f"m_{name}.sam")
        proc = run([optik, "map", "-r", ref, "-q", query_path, "-t", "4", "-o", out])
        check(f"map {name} exits 0", proc.returncode == 0, proc.stderr[-DETAIL_CHARS:])
        sams[name] = out

    def body(path):
        """Lines of a SAM file without its @PG header lines."""
        with open(path) as handle:
            return [
                text_line for text_line in handle if not text_line.startswith("@PG")
            ]

    check(
        "FASTA and gzipped FASTQ reads give the same alignments",
        [
            text_line.split("\t")[:SamCol.SEQ]
            for text_line in body(sams["fa"])
            if not text_line.startswith("@")
        ]
        == [
            text_line.split("\t")[:SamCol.SEQ]
            for text_line in body(sams["fqgz"])
            if not text_line.startswith("@")
        ],
    )
    # Writing to stdout must give the same records as -o.
    proc = run([optik, "map", "-r", ref, "-q", rfa, "-t", "4"])
    check(
        "map to stdout == map -o",
        [
            text_line
            for text_line in proc.stdout.splitlines(True)
            if not text_line.startswith("@PG")
        ]
        == body(sams["fa"]),
    )
    # samtools then reads the output as real pipelines do, sorting, indexing and recomputing MD.
    if samtools:
        sam_path = sams["fa"]
        check(
            "samtools quickcheck", run([samtools, "quickcheck", sam_path]).returncode == 0
        )
        flagstat_proc = run([samtools, "flagstat", sam_path])
        mapped_match = re.search(r"(\d+) \+ \d+ primary mapped", flagstat_proc.stdout)
        check(
            "samtools flagstat: every primary mapped",
            mapped_match is not None and int(mapped_match.group(1)) == len(reads),
            flagstat_proc.stdout[:DETAIL_CHARS],
        )
        bam = os.path.join(tmp_dir, "m.bam")
        passed = (
            run([samtools, "sort", "-o", bam, sam_path]).returncode == 0
            and run([samtools, "index", bam]).returncode == 0
        )
        check("samtools sort + index", passed)
        reg = run([samtools, "view", "-c", bam, "chrB"])
        check(
            "samtools view of a region",
            reg.returncode == 0 and int(reg.stdout.strip() or 0) > 0,
            reg.stderr[-DETAIL_CHARS:],
        )
        calmd = run(
            [samtools, "calmd", "-e", sam_path, ref]
        )  # recomputes NM and MD from the reference
        bad = 0
        with open(sam_path) as read_handle:
            # Name and flag together identify a record, since a read can have primary and supplementary lines.
            orig = {
                (text_line.split("\t")[SamCol.QNAME], text_line.split("\t")[SamCol.FLAG]): text_line
                for text_line in read_handle
                if not text_line.startswith("@")
            }
        for text_line in calmd.stdout.splitlines():
            if text_line.startswith("@"):
                continue
            fields = text_line.split("\t")
            nm_tags = [tag for tag in fields[SamCol.TAGS:] if tag.startswith("NM:i:")]
            orig_rec = orig.get((fields[SamCol.QNAME], fields[SamCol.FLAG]))
            onm = (
                [
                    tag
                    for tag in orig_rec.rstrip("\n").split("\t")[SamCol.TAGS:]
                    if tag.startswith("NM:i:")
                ]
                if orig_rec
                else []
            )
            bad += bool(nm_tags and onm and nm_tags[0] != onm[0])
        check(
            "samtools calmd agrees with every NM:i",
            calmd.returncode == 0 and bad == 0,
            f"{bad} differ",
        )
    else:
        check("samtools available for SAM checks", False, "samtools not on PATH")


def check_inputs_and_read_groups(fx):
    """Confirms that several -i or -q inputs behave like one joined file, then checks read groups."""
    optik, tmp_dir, ref, rfa, rfq, reads = fx.optik, fx.tmp_dir, fx.ref, fx.rfa, fx.rfq, fx.reads
    print("== several inputs, read groups")
    # Paired-end style, where -i R1 -i R2 of FASTA and gzipped FASTQ equals one joined file
    # Several -i inputs, plain and gzipped, must count exactly like one file holding all the reads.
    half = len(reads) // 2
    r1_path, r2gz, cat = (os.path.join(tmp_dir, file_name) for file_name in ("r1.fa", "r2.fq.gz", "r12.fa"))
    with open(r1_path, "w") as handle:
        handle.writelines(f">{mate_name}\n{rec_seq}\n" for mate_name, rec_seq in reads[:half])
    with gzip.open(r2gz, "wt") as handle:
        for mate_name, rec_seq in reads[half:]:
            handle.write(f"@{mate_name}\n{rec_seq}\n+\n{'I' * len(rec_seq)}\n")
    shutil.copy(rfa, cat)
    for cmd in (
        ["dbruijn", "-c", "1"],
        ["build-fm"],
        ["build-spss", "--compact"],
        ["build-spss"],
    ):
        one, two = os.path.join(tmp_dir, "one.out"), os.path.join(tmp_dir, "two.out")
        joined_proc = run([optik] + cmd + ["-i", cat, "-k", "31", "-o", one])
        split_proc = run([optik] + cmd + ["-i", r1_path, "-i", r2gz, "-k", "31", "-o", two])
        same = (
            joined_proc.returncode == 0
            and split_proc.returncode == 0
            and Path(one).read_bytes() == Path(two).read_bytes()
        )
        check(
            f"{' '.join(cmd)} -i R1 -i R2.gz == the concatenated input",
            same,
            split_proc.stderr.strip()[-DETAIL_CHARS:],
        )
    # -R must add the @RG line to the header and an RG:Z tag to every record, and @PG must record the command line.
    rgsam = os.path.join(tmp_dir, "rg.sam")
    proc = run(
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
    lines = Path(rgsam).read_text().splitlines() if proc.returncode == 0 else []
    body = [text_line for text_line in lines if not text_line.startswith("@")]
    check(
        "map -R: @RG header line and RG:Z:s1 on every record",
        "@RG\tID:s1\tSM:x\tPL:ILLUMINA" in lines
        and body
        and all(text_line.endswith("\tRG:Z:s1") for text_line in body),
        proc.stderr.strip()[-DETAIL_CHARS:],
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


def check_paired_end(fx):
    """Maps paired-end reads and checks proper pairs and mate fields.

    R1 and R2 files that don't match must be rejected.
    """
    optik, tmp_dir, ref = fx.optik, fx.tmp_dir, fx.ref
    print("== paired-end")
    # Paired mode must pair every record, set the read1 and read2 flags once per name and pass samtools' checks.
    pre = os.path.join(tmp_dir, "pe")
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
    pe1, pe2, pesam = pre + "_1.fq", pre + "_2.fq", os.path.join(tmp_dir, "pe.sam")
    proc = run([optik, "map", "-r", ref, "-q", pe1, "-q", pe2, "-o", pesam])
    body = []
    if proc.returncode == 0:
        with open(pesam) as read_handle:
            body = [line.split("\t") for line in read_handle if not line.startswith("@")]
    prim = [fields for fields in body if not int(fields[SamCol.FLAG]) & SamFlag.SECONDARY_SUPPLEMENTARY]
    n_prop = sum(1 for fields in prim if int(fields[SamCol.FLAG]) & SamFlag.PROPER_PAIR)
    check(
        "map -q R1 -q R2: every record paired (0x1, one 0x40 + one 0x80 per name), >= 95% proper",
        len(prim) == 4000
        and all(int(fields[SamCol.FLAG]) & SamFlag.PAIRED for fields in body)
        and sum(1 for fields in prim if int(fields[SamCol.FLAG]) & SamFlag.READ1) == 2000
        and n_prop >= MIN_PROPER_SHARE * 4000,
        f"records={len(prim)} proper={n_prop} {proc.stderr.strip()[-DETAIL_CHARS:]}",
    )
    samtools_bin = shutil.which("samtools")
    if samtools_bin:
        quickcheck_proc = subprocess.run(
            [samtools_bin, "quickcheck", pesam], capture_output=True, check=False
        )
        flagstat_out = subprocess.run(
            [samtools_bin, "flagstat", pesam], capture_output=True, text=True, check=False
        ).stdout
        check(
            "paired SAM: samtools quickcheck, flagstat counts read1/read2",
            quickcheck_proc.returncode == 0 and "2000 + 0 read1" in flagstat_out and "2000 + 0 read2" in flagstat_out,
            flagstat_out.replace("\n", " | ")[:DETAIL_CHARS],
        )
    # R2 with a record missing or out of order, or a third -q
    # R2 files that are one record short or have two records swapped must be refused, since mates are paired by position.
    lines2 = Path(pe2).read_text().splitlines(True)
    short2, swap2 = os.path.join(tmp_dir, "short_2.fq"), os.path.join(tmp_dir, "swap_2.fq")
    Path(short2).write_text("".join(lines2[:-FASTQ_LINES_PER_RECORD]))
    Path(swap2).write_text(
        "".join(
            lines2[FASTQ_LINES_PER_RECORD : 2 * FASTQ_LINES_PER_RECORD]
            + lines2[0:FASTQ_LINES_PER_RECORD]
            + lines2[2 * FASTQ_LINES_PER_RECORD :]
        )
    )
    for nm_tags, args, want in (
        ("R2 shorter than R1", ["-q", pe1, "-q", short2], "fewer records"),
        ("mates out of order", ["-q", pe1, "-q", swap2], "mate names differ"),
        ("three read files", ["-q", pe1, "-q", pe2, "-q", pe1], "one -q"),
    ):
        proc = run(
            [optik, "map", "-r", ref] + args + ["-o", os.path.join(tmp_dir, "x.sam")]
        )
        check(
            f"paired error: {nm_tags}",
            proc.returncode == 1 and want in proc.stderr,
            proc.stderr.strip()[-DETAIL_CHARS:],
        )
    proc = run(
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
            os.path.join(tmp_dir, "x.sam"),
        ]
    )
    check(
        "paired error: R2 longer than R1",
        proc.returncode == 1 and "more records" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )


def check_graph_modes(fx):
    """Makes sure --graph auto, part and table build the same graph, and that part refuses -s."""
    optik, tmp_dir, rfa = fx.optik, fx.tmp_dir, fx.rfa
    print("== --graph")
    # All three graph modes must write the same GFA, including part below the k where auto would pick it.
    outs = {}
    for mode in ("auto", "part", "table"):
        for kmer_arg in ("15", "31"):
            gfa_out = os.path.join(tmp_dir, f"g_{mode}_{kmer_arg}.gfa")
            proc = run(
                [
                    optik,
                    "dbruijn",
                    "-i",
                    rfa,
                    "-k",
                    kmer_arg,
                    "-c",
                    "1",
                    "--graph",
                    mode,
                    "-o",
                    gfa_out,
                ]
            )
            outs[(mode, kmer_arg)] = Path(gfa_out).read_bytes() if proc.returncode == 0 else None
    check(
        "dbruijn --graph auto|part|table: the same GFA at k = 15 and 31 (part below the default 19 too)",
        all(
            outs[(graph_mode, kmer_arg)] is not None and outs[(graph_mode, kmer_arg)] == outs[("auto", kmer_arg)]
            for graph_mode in ("part", "table")
            for kmer_arg in ("15", "31")
        ),
    )
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            rfa,
            "--graph",
            "part",
            "-s",
            os.path.join(tmp_dir, "x.kdbg"),
            "-o",
            os.path.join(tmp_dir, "x.gfa"),
        ]
    )
    check(
        "--graph part with -s (no node table to archive): error",
        proc.returncode == 1 and "-s" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            rfa,
            "--graph",
            "fast",
            "-o",
            os.path.join(tmp_dir, "x.gfa"),
        ]
    )
    check(
        "--graph with an unknown mode: error",
        proc.returncode == 1 and "--graph" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )


def check_io_edge_cases(fx):
    """Covers writing to a full disk, output into a closed pipe and input from stdin."""
    optik, tmp_dir, ref, rfa, truth = fx.optik, fx.tmp_dir, fx.ref, fx.rfa, fx.truth
    print("== full disk, closed pipes, stdin")
    spss_idx = os.path.join(tmp_dir, "i.spss")
    # Every writer must fail loudly on a full disk (/dev/full, ENOSPC), index saves included, since the FM save used to return success.
    # Writing to /dev/full fails on the first flush, which must give exit 1 and an error rather than a crash or a core dump.
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
            os.path.join(tmp_dir, "x.gfa"),
            "-s",
            "/dev/full",
        ],
        [optik, "map", "-r", ref, "-q", rfa, "-o", "/dev/full"],
        [optik, "query", "-i", spss_idx, "-q", rfa, "--ids", "/dev/full"],
    ):
        proc = run(cmd)
        check(
            f"full disk: {cmd[1]} {' '.join(cmd[-2:])} -> error exit 1",
            proc.returncode == 1
            and "error" in proc.stderr
            and "panicked" not in proc.stderr,
            f"exit {proc.returncode}: {proc.stderr.strip()[-DETAIL_CHARS:]}",
        )
    with open("/dev/full", "w") as full:
        proc = subprocess.run(
            [optik, "query", "-i", spss_idx, "-q", rfa],
            stdout=full,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    check(
        "stdout on a full disk: error exit 1, no core dump",
        proc.returncode == 1 and "writing to stdout failed" in proc.stderr,
        f"exit {proc.returncode}",
    )
    # A pipe closed by | head ends quietly, as C tools do.
    proc = subprocess.run(
        f"'{optik}' map -r '{ref}' -q '{rfa}' | head -1",
        shell=True,
        capture_output=True,
        text=True,
        check=False,
    )
    check(
        "map | head -1: no panic, no core dump",
        proc.stdout.startswith("@HD")
        and "panicked" not in proc.stderr
        and "core" not in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    # --ids with the queries on stdin reads them in one pass.  A second read of stdin used to
    # be empty, and the report then said 0 k-mers.
    with open(rfa) as in_handle:
        proc = subprocess.run(
            [
                optik,
                "query",
                "-i",
                spss_idx,
                "-q",
                "-",
                "--ids",
                os.path.join(tmp_dir, "ids.txt"),
            ],
            stdin=in_handle,
            capture_output=True,
            text=True,
            check=False,
        )
    check(
        "query -q - --ids: counts from the same pass",
        parse_q(proc.stdout + proc.stderr) == truth,
        f"{parse_q(proc.stdout + proc.stderr)} vs {truth}",
    )
    proc = run([optik, "dbruijn", "-i", ref, "-o", os.path.join(tmp_dir, "g2.gfa")])
    check(
        "dbruijn on a genome with the reads default -c 2: warns to use -c 1",
        proc.returncode == 0 and "use -c 1" in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    deep = os.path.join(
        tmp_dir, "deep.fa"
    )  # the read set four times over, about 5x coverage
    with open(deep, "w") as handle:
        handle.write(Path(rfa).read_text() * 4)
    proc = run(
        [optik, "dbruijn", "-i", deep, "-k", "21", "-o", os.path.join(tmp_dir, "g3.gfa")]
    )
    check(
        "dbruijn on 5x reads with -c 2: no such warning",
        proc.returncode == 0 and "use -c 1" not in proc.stderr,
        proc.stderr.strip()[-DETAIL_CHARS:],
    )
    # The archive would overwrite the GFA it shares a path with.
    proc = run(
        [
            optik,
            "dbruijn",
            "-i",
            rfa,
            "-o",
            os.path.join(tmp_dir, "same.gfa"),
            "-s",
            os.path.join(tmp_dir, "same.gfa"),
        ]
    )
    check(
        "two outputs on one file: error",
        proc.returncode == 1 and "same file" in proc.stderr,
        proc.stderr[-DETAIL_CHARS:],
    )


def check_errors(fx):
    """Runs bad invocations, each of which has to fail with a useful message."""
    optik, tmp_dir, ref, rfa = fx.optik, fx.tmp_dir, fx.ref, fx.rfa
    print("== errors")
    # Each case names a bad invocation, gives its command and the text its error message must contain.
    empty = os.path.join(tmp_dir, "empty.fa")
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
                os.path.join(tmp_dir, "x"),
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
                os.path.join(tmp_dir, "x"),
            ],
            "m",
        ),
        ("bad -t", [optik, "map", "-r", ref, "-q", rfa, "-t", "abc"], "-t"),
        (
            "empty input index",
            [optik, "build-spss", "-i", empty, "-o", os.path.join(tmp_dir, "x")],
            "no k-mers",
        ),
        (
            "empty dbruijn input",
            [optik, "dbruijn", "-i", empty, "-o", os.path.join(tmp_dir, "x.gfa")],
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
        proc = run(cmd, stdin=subprocess.DEVNULL)
        msg = proc.stderr + proc.stdout
        check(
            f"error: {name}",
            proc.returncode != 0 and want in msg and "panicked" not in msg,
            f"exit {proc.returncode}: {msg.strip()[-DETAIL_CHARS:]}",
        )


def main():
    """Run every smoke check against the optik binary and exit non-zero on any failure."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--optik", default=os.path.join(ROOT, "target/release/optik"))
    opts = parser.parse_args()
    optik = opts.optik
    rng = random.Random(5)
    tmp_dir = tempfile.mkdtemp(prefix="optik_cli_")
    try:
        fx = make_fixtures(optik, tmp_dir, rng)
        for section in (check_query, check_dbruijn, check_map, check_inputs_and_read_groups, check_paired_end, check_graph_modes, check_io_edge_cases, check_errors):
            section(fx)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    print(f"cli_smoke: {'ALL PASS' if not FAILS else f'{len(FAILS)} FAILURE(S)'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
