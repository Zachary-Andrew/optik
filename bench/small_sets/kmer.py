"""Runs the k-mer section.

Query workloads are generated with brute-force truth, and every index is then built and queried against it.
"""

import os
import random
import re
import subprocess
from pathlib import Path

from diagnostics import DETAIL_CHARS
from dna import kmer_codes, mutate, read_fasta, revcomp, write_fasta
from external_oracles import kmc_count, kmc_db, kmc_hits, mem_budget_gb
from small_common import BIN, DATA
from units import BITS_PER_BYTE, BYTES_PER_MIB, NS_PER_SECOND

from small_sets.common import (
    ABSENT_KMER_ID,
    IDS_CHECK_FAILED,
    KMC_MEM_GB,
    KMER_ERR_RATE,
    KMER_INDEL_RATE,
    KMER_READ_LEN,
    KMER_READS,
    KMER_SEED,
    KMER_SETS,
    PARITY_M,
    RANDOM_KMERS,
    REVCOMP_PROB,
    SET_K,
    SUB_SHARE,
)


def kmer_workloads(work, unitigs, rng):
    """Write the query files that every tool reads.

    They hold reads taken inside single unitigs (100% positive), the same reads with 1% substitutions and 0.4% indels, 1M random 31-mers, and the real SRR5833294 reads from another strain."""
    long_u = [unitig for _, unitig in unitigs if len(unitig) >= KMER_READ_LEN]
    pos, err = [], []
    for i in range(KMER_READS):
        unitig = rng.choice(long_u)
        start = rng.randrange(0, len(unitig) - KMER_READ_LEN + 1)
        read_seq = unitig[start : start + KMER_READ_LEN]
        read_seq = revcomp(read_seq) if rng.random() < REVCOMP_PROB else read_seq
        pos.append((f"p{i}", read_seq))
        err.append((f"e{i}", mutate(read_seq, rng, KMER_ERR_RATE * SUB_SHARE, KMER_INDEL_RATE, KMER_INDEL_RATE)))
    rand = [
        (f"x{i}", "".join(rng.choice("ACGT") for _ in range(SET_K)))
        for i in range(RANDOM_KMERS)
    ]
    srr = list(read_fasta(os.path.join(DATA, "queries/SRR5833294.10K.fastq.gz")))
    files = {}
    for name, recs in (
        ("reads_from_index", pos),
        ("reads_1pct_err", err),
        ("random_31mers", rand),
        ("SRR5833294", srr),
    ):
        files[name] = os.path.join(work, f"q_{name}.fa")
        write_fasta(files[name], recs)
    return files


def kmer_truth(set_db, path, work_dir, threads, kmer_len=SET_K):
    """(k-mer windows, positives) expected for a query file.

    Every tool reports all len - k + 1 windows of a read as num_kmers, with windows containing N counted negative.
    The positives are KMC's count of query windows whose k-mer is in the set (kmc_hits)."""
    total = sum(max(0, len(seq) - kmer_len + 1) for _, seq in read_fasta(path))
    return total, kmc_hits(set_db, path, kmer_len, work_dir, threads)


def parse_counts(text):
    """(num_kmers, num_positive_kmers) from a query tool's report, None where absent."""
    kmers_match = re.search(r"num_kmers = (\d+)", text)
    pos_match = re.search(r"num_positive_kmers = (\d+)", text)
    return (int(kmers_match.group(1)) if kmers_match else None, int(pos_match.group(1)) if pos_match else None)


def sbwt_counts(path):
    """(windows, positives) of SBWT's search output, one id per window with -1 for absent."""
    total = pos = 0
    with open(path) as handle:
        for line in handle:
            ids = line.split()
            total += len(ids)
            pos += sum(1 for kmer_id in ids if kmer_id != "-1")
    return total, pos


def fmsi_counts(path):
    """(windows, positives) of FMSI's query output, one 0/1 symbol per window."""
    total = pos = 0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            # A query result line is the k-mer, a tab and one presence bit per window.
            if len(fields) == 2:
                total += len(fields[1])
                pos += fields[1].count("1")
    return total, pos


def section_kmer(bench, mode, work, threads, quick):
    """Build and query every k-mer index on SSHash's unitig sets and gate each row against KMC."""
    print(f"\n=== kmer ({mode})", flush=True)
    rng = random.Random(KMER_SEED)
    sets = ["ecoli"] if quick else list(KMER_SETS)
    for dataset in sets:
        fname, sshash_m = KMER_SETS[dataset]
        set_dir = os.path.join(work, f"kmer_{dataset}")
        os.makedirs(set_dir, exist_ok=True)
        unitigs_fa = os.path.join(set_dir, "unitigs.fa")
        if not os.path.exists(unitigs_fa):
            write_fasta(
                unitigs_fa, read_fasta(os.path.join(DATA, "unitigs_stitched", fname))
            )
        unitigs = list(read_fasta(unitigs_fa))
        # KMC gives the truth independently of every compared tool.
        # It counts the set's distinct k-mers and, per workload, the k-mer windows and how many of them it finds in the set.
        set_db = kmc_db(
            [unitigs_fa], SET_K, 1, os.path.join(set_dir, "kmc_set"), threads, True, KMC_MEM_GB
        )
        num_kmers = kmc_count(set_db)
        workloads = kmer_workloads(set_dir, unitigs, rng)
        truth = {
            workload: kmer_truth(set_db, path, os.path.join(set_dir, f"kmc_{workload}"), threads)
            for workload, path in workloads.items()
        }
        print(f"  {dataset}: {num_kmers} distinct k-mers; truth {truth}", flush=True)
        m_par = PARITY_M
        build_threads = threads if mode == "optimal" else 1
        # Queries run on one thread in both modes.
        # SSHash, SBWT and FMSI have no query threads, so a multi-threaded optik row would not compare like with like.
        # optik's own -t scaling gets a separate row below.
        query_threads = 1
        tools = {
            "optik-spss": {
                "build": [
                    BIN["optik"],
                    "build-spss",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-m",
                    str(m_par if mode == "parity" else 0),
                    "-t",
                    str(build_threads),
                    "-o",
                    f"{set_dir}/optik.spss",
                ],
                "files": [f"{set_dir}/optik.spss"],
                "query": lambda query_path, set_dir=set_dir, query_threads=query_threads: (
                    [
                        BIN["optik"],
                        "query",
                        "-i",
                        f"{set_dir}/optik.spss",
                        "-q",
                        query_path,
                        "-k",
                        "31",
                        "-t",
                        str(query_threads),
                    ],
                    None,
                ),
                "params": f"k=31 m={m_par if mode == 'parity' else 'auto'} build -t {build_threads}, query -t {query_threads}",
            },
            "optik-fm": {
                "build": [
                    BIN["optik"],
                    "build-fm",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-c",
                    "1",
                    "-t",
                    str(build_threads),
                    "-o",
                    f"{set_dir}/optik.fm",
                ],
                "files": [f"{set_dir}/optik.fm"],
                "query": lambda query_path, set_dir=set_dir, query_threads=query_threads: (
                    [
                        BIN["optik"],
                        "query",
                        "-i",
                        f"{set_dir}/optik.fm",
                        "-q",
                        query_path,
                        "-k",
                        "31",
                        "-t",
                        str(query_threads),
                    ],
                    None,
                ),
                "params": f"k=31 build -t {build_threads}, query -t {query_threads}",
            },
            "sshash": {
                # The SSHash 2026 authors' command uses -g 16 (or 80% of this machine's memory).
                "build": [
                    BIN["sshash"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-k",
                    "31",
                    "-m",
                    str(m_par if mode == "parity" else sshash_m),
                    "-g",
                    str(mem_budget_gb()),
                    "-t",
                    str(build_threads),
                    "-d",
                    set_dir,
                    "-o",
                    f"{set_dir}/idx.sshash",
                ],
                "files": [f"{set_dir}/idx.sshash"],
                "query": lambda query_path, set_dir=set_dir: (
                    [BIN["sshash"], "query", "-i", f"{set_dir}/idx.sshash", "-q", query_path],
                    None,
                ),
                "params": f"k=31 m={m_par if mode == 'parity' else sshash_m} build -t {build_threads}, query single-threaded",
            },
            "sbwt": {
                # The SSHash 2026 authors run SBWT with -m 16 (or 80% of this machine) and plain-matrix.
                "build": [
                    BIN["sbwt"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-o",
                    f"{set_dir}/idx.sbwt",
                    "-k",
                    "31",
                    "--variant",
                    "plain-matrix",
                    "--add-reverse-complements",
                    "-m",
                    str(mem_budget_gb()),
                    "-t",
                    str(build_threads),
                    "--temp-dir",
                    set_dir,
                ],
                "files": [f"{set_dir}/idx.sbwt"],
                "query": lambda query_path, set_dir=set_dir: (
                    [
                        BIN["sbwt"],
                        "search",
                        "-i",
                        f"{set_dir}/idx.sbwt",
                        "-q",
                        query_path,
                        "-o",
                        f"{set_dir}/sbwt.out",
                    ],
                    f"{set_dir}/sbwt.out",
                ),
                # Timed runs send the per-k-mer ids to /dev/null, where hyperfine also sends FMSI's stdout.
                # Neither tool then spends time on file writes the counters skip.
                "timed": lambda query_path, set_dir=set_dir: [
                    BIN["sbwt"],
                    "search",
                    "-i",
                    f"{set_dir}/idx.sbwt",
                    "-q",
                    query_path,
                    "-o",
                    "/dev/null",
                ],
                "answers": "per-k-mer",
                "params": f"k=31 plain-matrix (+ reverse complements) -m {mem_budget_gb()} build -t {build_threads}, "
                "search single-threaded",
            },
            "fmsi": {
                # FMSI indexes a masked superstring, so its build is KmerCamel plus fmsi index.
                # Both modes use the FMSI settings of the SSHash 2026 authors' kmer_sets_benchmark, kmercamel compute -S for SPSS input, fmsi index, and queries with the kLCP array and max-one masks through -S -O.
                # Those are FMSI's own optimizations, not knobs to match.
                "build": [
                    "bash",
                    "-c",
                    (
                        f"{BIN['kmercamel']} compute -S -k 31 "
                        f"-M {set_dir}/ms.fa -o /dev/null {unitigs_fa} && {BIN['fmsi']} index -k 31 {set_dir}/ms.fa"
                    ),
                ],
                "files_glob": f"{set_dir}/ms.fa.fmsi.",
                "query": lambda query_path, set_dir=set_dir: (
                    [
                        BIN["fmsi"],
                        "query",
                        "-S",
                        "-O",
                        "-q",
                        query_path,
                        "-k",
                        "31",
                        f"{set_dir}/ms.fa",
                    ],
                    f"{set_dir}/fmsi.out",
                ),
                "params": "k=31 KmerCamel greedy max-one superstring (-S unitig input) + fmsi query -S -O "
                "(kLCP streaming, max-one; the SSHash 2026 benchmark's settings)",
            },
        }

        # Dictionary workload mapping k-mers to dense ids with optik query --ids.
        # SBWT's search output is also per-k-mer ids, its colex ranks, whereas SSHash has no CLI for per-k-mer ids, since its Lookup exists only in the library and its benchmarks.
        def ids_check(path, qfile, set_dir=set_dir, num_kmers=num_kmers):
            """Count invalid entries of an optik --ids output for query file qfile.

            Each window's id must be present exactly where SBWT finds the k-mer, lie below n, and map one-to-one onto k-mers."""
            id_of, kmer_of, bad = {}, {}, 0
            lines = Path(path).read_text().split("\n")
            # Presence per window comes from SBWT's own search output on the same query file, an external tool writing one colex id per window and -1 when absent.
            sbwt_out = os.path.join(set_dir, "sbwt_ids.out")
            subprocess.run(
                [BIN["sbwt"], "search", "-i", f"{set_dir}/idx.sbwt", "-q", qfile, "-o", sbwt_out],
                check=True,
                capture_output=True,
            )
            present = [
                [int(id_str) != ABSENT_KMER_ID for id_str in text_line.split()]
                for text_line in Path(sbwt_out).read_text().split("\n")
            ]
            for (_, seq), line, pres in zip(
                read_fasta(qfile), lines, present, strict=False
            ):
                vals = list(map(int, line.split())) if line else []
                if len(vals) != max(0, len(seq) - SET_K + 1) or len(pres) != len(vals):
                    bad += 1
                    continue
                for i, kmer_id in enumerate(vals):
                    codes = kmer_codes(seq[i : i + SET_K], SET_K)
                    code = codes[0] if codes else None
                    if (
                        pres[i] != (kmer_id >= 0)
                        or kmer_id >= 0
                        and (
                            kmer_id >= num_kmers
                            or id_of.setdefault(code, kmer_id) != kmer_id
                            or kmer_of.setdefault(kmer_id, code) != code
                        )
                    ):
                        bad += 1
            return bad

        for tool, spec in tools.items():
            # Inputs are pre-read before each timed run (measure.py).
            result = bench.measure(
                f"{mode}/{dataset}/{tool}/build", spec["build"], cwd=set_dir, inputs=[unitigs_fa]
            )
            if not result.succeeded:
                bench.add(
                    section="kmer",
                    mode=mode,
                    dataset=dataset,
                    tool=tool,
                    phase="build",
                    gate=False,
                    gate_detail=f"build exit {result.code}: {result.err[-DETAIL_CHARS:]}",
                    params=spec["params"],
                )
                continue
            if tool == "optik-spss" and mode != "parity":
                # With -m 0 optik picks m from the input, and the row names the m it picked.
                minimizer_match = re.search(
                    r"minimizer_len\t(\d+)",
                    Path(result.out).read_text() if os.path.exists(result.out) else "",
                )
                spec["params"] = spec["params"].replace(
                    "m=auto", f"m=auto ({minimizer_match.group(1) if minimizer_match else '?'})"
                )
            idx = (
                list(spec["files"])
                if "files" in spec
                else [
                    os.path.join(set_dir, fname)
                    for fname in os.listdir(set_dir)
                    if os.path.join(set_dir, fname).startswith(spec["files_glob"])
                ]
            )
            size = sum(os.path.getsize(path) for path in idx)
            bench.add(
                section="kmer",
                mode=mode,
                dataset=dataset,
                tool=tool,
                phase="build",
                build_s=result.wall_row(),
                build_rss_mb=round(result.rss(), 1),
                index_mb=round(size / BYTES_PER_MIB, 3),
                bits_per_kmer=round(BITS_PER_BYTE * size / num_kmers, 3),
                gate=True,
                gate_detail="",
                params=spec["params"],
            )
            per_kmer = (
                spec.get("answers") or ("per-k-mer" if tool == "fmsi" else "counts")
            ) == "per-k-mer"
            for workload, query_path in workloads.items():
                cmd, outf = spec["query"](query_path)
                capt = outf or os.path.join(set_dir, "q.out")
                if per_kmer:
                    # One untimed run writes the per-k-mer answers to a file for the gate.
                    # Timed runs write them to /dev/null, so their timings include no output writes.
                    rc_, _, _, err = bench.once(
                        f"{mode}/{dataset}/{tool}/{workload} (answers)",
                        cmd,
                        cwd=set_dir,
                        stdout=None if tool in ("sbwt",) else capt,
                    )
                    result = (
                        bench.measure(
                            f"{mode}/{dataset}/{tool}/{workload}",
                            spec["timed"](query_path) if "timed" in spec else cmd,
                            cwd=set_dir,
                            inputs=[query_path, *idx],
                        )
                        if rc_ == 0
                        else None
                    )
                else:
                    result = bench.measure(
                        f"{mode}/{dataset}/{tool}/{workload}",
                        cmd,
                        cwd=set_dir,
                        stdout=capt,
                        inputs=[query_path, *idx],
                    )
                    rc_, err = result.code, result.err
                if rc_ != 0 or not result.succeeded:
                    bench.add(
                        section="kmer",
                        mode=mode,
                        dataset=dataset,
                        tool=tool,
                        phase=workload,
                        gate=False,
                        gate_detail=f"query exit {rc_ or result.code}: {(err or result.err)[-DETAIL_CHARS:]}",
                        params=spec["params"],
                    )
                    continue
                if tool == "sbwt":
                    got = sbwt_counts(outf)
                elif tool == "fmsi":
                    got = fmsi_counts(capt)
                else:
                    got = parse_counts(Path(capt).read_text() + err)
                if tool == "fmsi":
                    # FMSI prints one symbol per base rather than per k-mer for a read containing N, so only its positives can be compared.
                    passed = got[1] == truth[workload][1]
                else:
                    passed = got == truth[workload]
                num_windows = truth[workload][0]
                bench.add(
                    section="kmer",
                    mode=mode,
                    dataset=dataset,
                    tool=tool,
                    phase=workload,
                    answers=spec.get("answers")
                    or ("per-k-mer" if tool == "fmsi" else "counts"),
                    query_s=result.wall_row(),
                    ns_per_kmer=round(result.wall["mean"] * NS_PER_SECOND / num_windows, 1),
                    query_rss_mb=round(result.rss(), 1),
                    kmers=got[0],
                    positive=got[1],
                    gate=passed,
                    gate_detail="" if passed else f"counts {got} != truth {truth[workload]}",
                    params=spec["params"],
                )
                if (
                    mode == "optimal"
                    and threads > 1
                    and tool in ("optik-spss", "optik-fm")
                ):
                    # optik's own query scaling runs the same query on all threads and must give the same answers.
                    cmd_t = [
                        arg if arg != str(query_threads) or cmd[i - 1] != "-t" else str(threads)
                        for i, arg in enumerate(cmd)
                    ]
                    result = bench.measure(
                        f"{mode}/{dataset}/{tool}/{workload}/t{threads}",
                        cmd_t,
                        cwd=set_dir,
                        stdout=capt,
                        inputs=[query_path, *idx],
                    )
                    got_t = (
                        parse_counts(Path(capt).read_text() + result.err) if result.succeeded else None
                    )
                    ok_t = got_t is not None and got_t[1] == truth[workload][1]
                    bench.add(
                        section="kmer",
                        mode=mode,
                        dataset=dataset,
                        tool=f"{tool} -t {threads} (own scaling)",
                        phase=workload,
                        answers="counts",
                        query_s=result.wall_row(),
                        ns_per_kmer=round(result.wall["mean"] * NS_PER_SECOND / num_windows, 1)
                        if result.succeeded
                        else None,
                        query_rss_mb=round(result.rss(), 1) if result.succeeded else None,
                        kmers=got_t[0] if got_t else None,
                        positive=got_t[1] if got_t else None,
                        gate=ok_t,
                        gate_detail=""
                        if ok_t
                        else f"counts {got_t} != truth {truth[workload]}",
                        params=spec["params"]
                        + f"; query -t {threads}: not a comparison row (others cannot thread)",
                    )
        # optik's dictionary ids are checked once every index exists, with presence judged by SBWT's search.
        if os.path.exists(f"{set_dir}/optik.spss") and os.path.exists(f"{set_dir}/idx.sbwt"):
            spec = tools["optik-spss"]
            query_path = workloads["reads_1pct_err"]
            cmd = [
                BIN["optik"],
                "query",
                "-i",
                f"{set_dir}/optik.spss",
                "-q",
                query_path,
                "-k",
                "31",
                "-t",
                str(query_threads),
                "--ids",
                f"{set_dir}/ids.txt",
            ]
            rc_, _, _, err = bench.once(
                f"{mode}/{dataset}/optik-spss/ids (answers)",
                cmd,
                cwd=set_dir,
                stdout=os.path.join(set_dir, "q.out"),
            )
            bad = ids_check(f"{set_dir}/ids.txt", query_path) if rc_ == 0 else IDS_CHECK_FAILED
            # Timed with the ids sent to /dev/null, like the other per-k-mer answers.
            result = bench.measure(
                f"{mode}/{dataset}/optik-spss/ids",
                cmd[:-1] + ["/dev/null"],
                cwd=set_dir,
                inputs=[query_path, f"{set_dir}/optik.spss"],
            )
            bench.add(
                section="kmer",
                mode=mode,
                dataset=dataset,
                tool="optik-spss --ids",
                phase="reads_1pct_err",
                answers="per-k-mer",
                query_s=result.wall_row(),
                ns_per_kmer=round(result.wall["mean"] * NS_PER_SECOND / truth["reads_1pct_err"][0], 1)
                if result.succeeded
                else None,
                query_rss_mb=round(result.rss(), 1) if result.succeeded else None,
                gate=bad == 0 and result.succeeded,
                gate_detail=""
                if bad == 0
                else f"{bad} id inconsistencies (presence vs SBWT, range, bijection)",
                params=spec["params"]
                + " --ids (dense k-mer ids, SSHash Lookup semantics)",
            )
