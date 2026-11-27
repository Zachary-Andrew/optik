"""Runs the k-mer query section.

After preparing the SPSS inputs and query sets, it times every tool's builds and lookups.
"""

import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

from external_oracles import kmc_count, kmc_db, kmc_diff, kmc_hits
from published_catalog import DATA_ROOT
from published_inputs import (
    join_mates,
    query_segments,
    random_kmers,
    sample_kmers,
    spss_input,
    sshash_m,
    windows_in,
)
from published_tools import (
    NO_BUILD_THREADS,
    OWN_OUTPUT,
    count_of,
    ggcat_bin,
    ggcat_version,
    idx_files,
    query_tools,
)
from small_common import BIN
from units import BITS_PER_BYTE, BYTES_PER_GIB, BYTES_PER_MIB, NS_PER_SECOND

from published.common import (
    KMC_MAX_MEM_GB,
    KMC_MEM_SHARE,
    LOAD_PROBE_SEED,
    MAX_QUERIES,
    REFERENCE_K,
    SMOKE_MAX_QUERIES,
    SMOKE_STREAM_BASES,
    SPSS_TRIES,
    SUBSAMPLED_K,
    add_status,
    failure,
)
from published.machine import machine, tool_skip

# A failed GGCAT try is summarised by the earlier of the last two lines of its log.
LOG_TAIL_LINES = 2
# The re-raised error names only the program and its first argument.
SHOWN_CMD_WORDS = 2


def make_spss(inputs, kmer_len, out, threads, unitigs=None):
    """Build the SPSS of a k-mer set that every index is built from.

    GGCAT eulertigs from the papers' tool at its latest version are used when GGCAT is built, and otherwise Cuttlefish's unitigs, the common k-mer set of scaling_bench.
    `unitigs` gives the reason GGCAT's plain unitigs stand in for its eulertigs.
    An SPSS is accepted only when its k-mer set equals the input's, meaning KMC finds the same distinct count and nothing outside the input.
    GGCAT's eulertig step aborts on some sets and silently drops k-mers on others, as reported in issue #65, and chm13 at k = 23 lost 2 of 2,337,245,516 in one run.
    A failed attempt is therefore retried, and after SPSS_TRIES failures the unitigs are used.
    Returns the path and the tool that built it."""
    # A finished SPSS is reused, together with the note of how it was built.
    if os.path.exists(out):
        return out, Path(out + ".by").read_text().strip() if os.path.exists(
            out + ".by"
        ) else "?"
    work = os.path.dirname(out)
    kmc_gb = max(1, min(KMC_MAX_MEM_GB, int(machine()["mem_gib"] * KMC_MEM_SHARE)))
    # KMC's k-mer set of the raw input is the yardstick every attempt is checked against.
    raw_db = kmc_db(inputs, kmer_len, 1, os.path.join(work, "kmc_raw"), threads, True, kmc_gb)
    raw_n = kmc_count(raw_db)
    # GGCAT eulertigs get SPSS_TRIES attempts, since it sometimes aborts or drops k-mers (issue #65), and maximal unitigs are the fallback.
    tries = [] if unitigs or not os.path.exists(ggcat_bin()) else [False] * SPSS_TRIES
    notes = []
    for use_unitigs in tries + [True]:
        try:
            built_by = _make_spss(
                inputs,
                kmer_len,
                out + ".try",
                threads,
                work,
                use_unitigs
                and (
                    unitigs
                    or f"eulertigs: {len(tries)} GGCAT attempts failed, issue #65"
                ),
            )
        except subprocess.CalledProcessError as err:
            tail = (
                Path(out + ".try.log")
                .read_text(errors="replace")
                .strip()
                .splitlines()[-LOG_TAIL_LINES:]
                if os.path.exists(out + ".try.log")
                else []
            )
            notes.append(
                "aborted"
                + (f" ({tail[0]})" if tail else "")
            )
            if use_unitigs:
                raise subprocess.CalledProcessError(
                    err.returncode, err.cmd[:SHOWN_CMD_WORDS], output="; ".join(notes)
                ) from None
            continue
        # An attempt is accepted only when its k-mer set equals the input's exactly, with nothing missing and nothing extra.
        spss_db = kmc_db(
            [out + ".try"], kmer_len, 1, os.path.join(work, "kmc_try"), threads, True, kmc_gb
        )
        spss_distinct, extra = (
            kmc_count(spss_db),
            kmc_diff(spss_db, raw_db, os.path.join(work, "kmc_try_extra")),
        )
        if spss_distinct == raw_n and extra == 0:
            break
        notes.append(
            f"k-mer set differs from the input's ({spss_distinct} vs {raw_n} distinct, {extra} not in the input)"
        )
        if use_unitigs:
            raise subprocess.CalledProcessError(
                1, ["ggcat", "build"], output="; ".join(notes)
            )
    # The accepted attempt is renamed into place, so a half-written SPSS is never taken for a finished one.
    for ext in ("", ".log"):
        if os.path.exists(out + ".try" + ext):
            os.replace(out + ".try" + ext, out + ext)
    if notes:
        built_by += (
            f" (accepted on attempt {len(notes) + 1}; earlier: "
            + "; ".join(notes)
            + ")"
        )
    Path(out + ".by").write_text(built_by + "\n")
    return out, built_by


def _make_spss(inputs, kmer_len, out, threads, work, unitigs=None):
    """Run GGCAT once, or Cuttlefish when GGCAT is absent, and format its output as `out`.
    Returns a description of the tool and mode used."""
    if os.path.exists(ggcat_bin()):
        # GGCAT takes its inputs as a list file, and -s 1 keeps every k-mer so the SPSS covers the whole input.
        lst = out + ".list"
        Path(lst).write_text("\n".join(inputs) + "\n")
        with open(out + ".log", "w") as log_handle:
            subprocess.run(
                [
                    ggcat_bin(),
                    "build",
                    "-k",
                    str(kmer_len),
                    "-j",
                    str(threads),
                    "-s",
                    "1",
                    *([] if unitigs else ["--eulertigs"]),
                    "-l",
                    lst,
                    "-o",
                    out + ".raw.fa",
                    "--temp-dir",
                    work,
                ],
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
        # The note records the tool, its version and the mode, because the SPSS a query runs against changes its timings.
        built_by = f"GGCAT {ggcat_version()} " + (
            f"unitigs ({unitigs})" if unitigs else "--eulertigs"
        )
    else:
        with open(out + ".log", "w") as log_handle:
            subprocess.run(
                [
                    BIN["cuttlefish"],
                    "build",
                    "--ref",
                    *[arg for input_path in inputs for arg in ("-s", input_path)],
                    "-k",
                    str(kmer_len),
                    "-t",
                    str(threads),
                    "-o",
                    out + ".cf",
                    "-w",
                    work,
                ],
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
        os.replace(out + ".cf.fa", out + ".raw.fa")
        built_by = "Cuttlefish unitigs (GGCAT not built)"
    for path in (out, out + ".part"):
        if os.path.exists(path):
            os.remove(path)
    spss_input([out + ".raw.fa"], out, kmer_len)
    os.remove(out + ".raw.fa")
    return built_by


def query_inputs(ds_id, dataset, kmer_len, variant, work, threads, limit_bases=None, max_queries=MAX_QUERIES):
    """Prepare the formatted SPSS and query files of one (dataset, k, variant).

    Both query modes share them.
    `limit_bases` cuts a slice for `smoke`."""
    src = os.path.join(DATA_ROOT.path, ds_id)
    base = os.path.join(work, f"{ds_id}_k{kmer_len}_{variant}" + ("_smoke" if limit_bases else ""))
    os.makedirs(base, exist_ok=True)
    unitigs_fa = os.path.join(base, "input.fa")
    built_by = "deposit (as published)"
    if ds_id.startswith("sshash26"):
        spss_input([os.path.join(src, dataset["input"])], unitigs_fa, kmer_len, limit_bases)
    elif variant == "sub10":
        spss_input([os.path.join(src, dataset["subsampled"])], unitigs_fa, kmer_len, limit_bases)
    else:
        if ds_id == "sbwt23-ecoli3682":
            extract_dir = os.path.join(src, "genomes")
            if not os.path.exists(extract_dir):
                with tarfile.open(os.path.join(src, "coli3682_dataset.tar.gz")) as tar:
                    tar.extractall(extract_dir, filter="data")
            raw = sorted(
                os.path.join(root, fname)
                for root, _, files in os.walk(extract_dir)
                for fname in files
                if fname.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz"))
            )
        elif dataset.get("reads"):
            raw = sorted(
                os.path.join(src, fname) for fname in os.listdir(src) if fname.endswith(".fastq.gz")
            )
            if dataset.get("join_mates"):
                joined = os.path.join(base, "spots.fa")
                join_mates(raw, joined, limit_bases)
                raw = [joined]
        elif dataset.get("local"):
            raw = [dataset["local"]]
        else:
            raw = [os.path.join(src, dataset["input"])]
        plain = os.path.join(base, "plain.fa")
        spss_input(raw, plain, kmer_len, limit_bases, dataset.get("delete_nonacgt", False))
        unitigs_fa, built_by = make_spss(
            [plain], kmer_len, os.path.join(base, "spss.fa"), threads, dataset.get("spss_unitigs")
        )
    num_queries = min(max_queries, SMOKE_MAX_QUERIES) if limit_bases else max_queries
    workloads = {
        "pos_halfrc": sample_kmers(unitigs_fa, kmer_len, num_queries, os.path.join(base, "q_pos.fa")),
        "neg_random": random_kmers(kmer_len, num_queries, os.path.join(base, "q_neg.fa")),
        "_load": random_kmers(kmer_len, 1, os.path.join(base, "q_load.fa"), seed=LOAD_PROBE_SEED),
    }
    # stream_k (k-scaling rounds) restricts the streaming workload to these k.
    # FMSI's k sweep measures size and lookups, whereas SSHash's streaming numbers are a k = 31 measurement.
    if (
        dataset.get("stream")
        and (kmer_len == REFERENCE_K or dataset.get("stream_from"))
        and kmer_len in dataset.get("stream_k", [kmer_len])
    ):
        stream_fa = os.path.join(base, "q_stream.fa")
        query_segments(
            os.path.join(DATA_ROOT.path, dataset.get("stream_from", ds_id), dataset["stream"]),
            stream_fa,
            kmer_len,
            SMOKE_STREAM_BASES if limit_bases else None,
        )
        workloads["stream_" + dataset["stream"].split(".")[0]] = stream_fa
    _, _, kmers, *_ = spss_input([], unitigs_fa, kmer_len)
    return unitigs_fa, built_by, workloads, kmers


def bench_query(bench, ds_id, dataset, mode, work, threads, mach, first_t=True, sweep_t=False):
    """Run one thread count `threads` of a dataset.

    first_t marks the sweep's first, single-thread round, which also runs the tools with no thread option for a step.
    Later rounds rebuild only the multi-threaded builders and rerun only optik's queries (sweep_t runs optik queries at `threads`)."""
    for kmer_len in dataset["k"]:
        variants = ["full"] + (["sub10"] if dataset.get("subsampled") and kmer_len == SUBSAMPLED_K else [])
        for variant in variants:
            unitigs_fa, built_by, workloads, _windows = query_inputs(
                ds_id, dataset, kmer_len, variant, work, threads
            )
            load = workloads.pop("_load")
            # The gates below judge only the rows this round adds.
            # A k-scaling and a thread-scaling round can share (dataset, k, variant, mode, threads), and the other round's rows were already judged against its own answers.
            row0 = len(bench.rows)
            # Each mode and thread count gets a fresh directory for its indexes and answers, so a tool that exits 0 without writing is never scored on an older file.
            run_dir = os.path.join(work, f"{ds_id}_k{kmer_len}_{variant}", f"{mode}_t{threads}")
            shutil.rmtree(run_dir, ignore_errors=True)
            os.makedirs(run_dir)
            pub = dataset.get("kmers31") if kmer_len == REFERENCE_K and variant == "full" else None
            # KMC counts the indexed set's distinct k-mers independently of every compared tool.
            # When the paper publishes a count, the two must agree, which also shows that formatting kept the paper's set.
            kmc_gb = max(1, min(KMC_MAX_MEM_GB, int(mach["mem_gib"] * KMC_MEM_SHARE)))
            kdb = os.path.join(os.path.dirname(run_dir), "kmc_set")
            kmc_n = kmc_count(kmc_db([unitigs_fa], kmer_len, 1, kdb, threads, True, kmc_gb))
            kmers = pub or kmc_n
            input_ok = pub is None or kmc_n == pub
            if pub is not None and dataset.get("kmers31_forward"):
                # The paper counts k-mers as written, so KMC runs without canonical form on the reads.
                fwd = kmc_count(
                    kmc_db(
                        [os.path.join(os.path.dirname(run_dir), "plain.fa")],
                        kmer_len,
                        1,
                        os.path.join(os.path.dirname(run_dir), "kmc_forward"),
                        threads,
                        True,
                        kmc_gb,
                        canonical=False,
                    )
                )
                input_ok = abs(fwd - pub) <= dataset.get("kmers31_tol", 0) * pub
                kmers = kmc_n
            if first_t:
                forward_checked = dataset.get("kmers31_forward") and pub is not None
                bench.add(
                    section="query",
                    mode=mode,
                    dataset=ds_id,
                    variant=variant,
                    k=kmer_len,
                    tool="kmc (input)",
                    phase="input",
                    kmers_reported=kmc_n,
                    **({"kmers_forward": fwd} if forward_checked else {}),
                    kmers_published=pub,
                    gate=input_ok,
                    gate_detail=(
                        f"forward k-mers {fwd} vs the paper's {pub} ({fwd / pub - 1:+.2%}, "
                        f"allowed +-{dataset.get('kmers31_tol', 0):.1%}: see the dataset's note)"
                        if forward_checked
                        else ""
                        if input_ok
                        else f"KMC counts {kmc_n} distinct k-mers in "
                        f"the formatted input, the paper publishes {pub}"
                    ),
                    params=f"SPSS: {built_by}"
                    + ("; reads = spots, mates joined" if dataset.get("join_mates") else "")
                    + (
                        "; non-ACGT deleted within records"
                        if dataset.get("delete_nonacgt")
                        else ""
                    ),
                )
            expected = {workload: windows_in(query_path, kmer_len) for workload, query_path in workloads.items()}
            # KMC's answer for every workload is the set of query k-mer windows it finds.
            kmc_truth = {
                workload: kmc_hits(
                    kdb,
                    query_path,
                    kmer_len,
                    os.path.join(os.path.dirname(run_dir), f"kmc_{workload}"),
                    threads,
                    kmc_gb,
                )
                for workload, query_path in workloads.items()
                if workload != "_load"
            }
            # SSHash requires -m and has no default, so the paper's m is used where given and k/2 otherwise.
            minimizer_len = sshash_m(dataset, kmer_len, kmc_n)
            tools = query_tools(
                mode,
                run_dir,
                unitigs_fa,
                kmer_len,
                minimizer_len,
                threads,
                mach,
                bench.cap // BYTES_PER_GIB,
                threads if sweep_t else 1,
            )
            res, kmers_reported = {}, {}
            for tool, spec in tools.items():
                if not first_t and tool in NO_BUILD_THREADS:
                    continue  # without a thread option, its t = 1 row stands
                tag = f"{mode}/{ds_id}/k{kmer_len}/{variant}/t{threads}/{tool}"
                # FMSI's build takes no thread option, so its row records what it ran with.
                base = {
                    "section": "query",
                    "mode": mode,
                    "dataset": ds_id,
                    "variant": variant,
                    "k": kmer_len,
                    "threads": 1 if tool in NO_BUILD_THREADS else threads,
                    "tool": tool,
                    "answers": spec["answers"],
                    "params": spec["params"] + f"; SPSS: {built_by}",
                }
                # A command that is the same in both modes (every compared tool, and optik-fm) is measured once in the parity round, and its rows also stand in the optimal table.
                # optik-spss differs (-m), and its --ids view queries that index.
                same = [
                    row
                    for row in bench.rows
                    if mode == "optimal"
                    and not tool.startswith("optik-spss")
                    and row.get("mode") == "parity"
                    and row.get("section") == "query"
                    and (
                        row.get("dataset"),
                        row.get("k"),
                        row.get("variant"),
                        row.get("threads"),
                        row.get("tool"),
                    )
                    == (ds_id, kmer_len, variant, base["threads"], tool)
                    and row.get("params") == base["params"]
                ]
                if same:
                    for row in same:
                        copied_row = dict(
                            row,
                            mode="optimal",
                            params=row["params"]
                            + "; the parity run's command: measured once",
                        )
                        bench.rows.append(copied_row)
                        if row.get("phase") in workloads and row.get("kmers") is not None:
                            res.setdefault(row["phase"], {})[tool] = (
                                row["kmers"],
                                row["positive"],
                            )
                        if row.get("phase") == "build":
                            kmers_reported[tool] = row.get("kmers_reported")
                    continue
                if spec["build"] is None:
                    # This is a query mode of another tool's index, with no build of its own.
                    if spec["reuse"] not in kmers_reported:
                        continue  # that index was skipped or failed
                    lcmd, _ = spec["query"](load)
                    idx = idx_files(tools[spec["reuse"]], run_dir)
                    load_timing = (
                        bench.measure(
                            tag + "/load",
                            spec.get("timed", lambda query_path, lcmd=lcmd: lcmd)(load),
                            cwd=run_dir,
                            inputs=[load, *idx],
                        ).wall_row()
                        or {}
                    )
                if spec["build"] is not None:
                    why = tool_skip(dataset, tool, mach, kmer_len, kmers)
                    if why:
                        add_status(
                            bench, "skip", **base, phase="build", gate_detail=why
                        )
                        continue
                    bout = os.path.join(run_dir, f"{tool}.build.out")
                    # The SPSS is read in before each run (measure.py), so all runs see one cache state.
                    build_res = bench.measure(
                        tag + "/build",
                        spec["build"],
                        cwd=run_dir,
                        stdout=bout,
                        inputs=[unitigs_fa],
                    )
                    err = (
                        Path(bout).read_text() if os.path.exists(bout) else ""
                    ) + build_res.err
                    if not build_res.succeeded:
                        kind, det = failure(
                            bench, build_res.code, build_res.rss() or bench.cap / BYTES_PER_MIB, err
                        )
                        if kind == "cap":
                            add_status(
                                bench, "cap", **base, phase="build", gate_detail=det
                            )
                        else:
                            bench.add(
                                **base,
                                phase="build",
                                gate=False,
                                gate_detail="build " + det,
                            )
                        continue
                    idx = idx_files(spec, run_dir)
                    size = sum(os.path.getsize(path) for path in idx)
                    distinct_match = re.search(r"distinct_kmers\s+(\d+)", err) or re.search(
                        r"num_kmers\s*=?\s*(\d+)", err
                    )
                    kmers_reported[tool] = int(distinct_match.group(1)) if distinct_match else None
                    # A load-only run, which the per-k-mer time below subtracts so that index loading, seconds at this scale, does not count as query cost.
                    lcmd, _ = spec["query"](load)
                    load_timing = (
                        bench.measure(
                            tag + "/load",
                            spec["timed"](load) if "timed" in spec else lcmd,
                            cwd=run_dir,
                            inputs=[load, *idx],
                        ).wall_row()
                        or {}
                    )
                    extra = {}
                    if spec.get("extra"):
                        ecmd, efile = spec["extra"]
                        if (
                            bench.once(
                                tag + "/build-norc",
                                ecmd,
                                cwd=run_dir,
                                stdout=os.path.join(run_dir, "norc.out"),
                            )[0]
                            == 0
                        ):
                            extra["index_mb_without_rc"] = round(
                                os.path.getsize(efile) / BYTES_PER_MIB, 2
                            )
                    bench.add(
                        **base,
                        phase="build",
                        build_s=build_res.wall_row(),
                        build_rss_mb=round(build_res.rss(), 1),
                        index_mb=round(size / BYTES_PER_MIB, 2),
                        **extra,
                        kmers_reported=kmers_reported[tool],
                        load_s=load_timing,
                        gate=True,
                        gate_detail="",
                    )
                for workload, query_path in workloads.items():
                    if dataset.get("no_queries"):
                        break  # a build-only round leaves the queries to another round
                    if not first_t and not tool.startswith("optik"):
                        break  # single-threaded queries belong to the t = 1 round
                    cmd, outf = spec["query"](query_path)
                    capt = (
                        os.path.join(run_dir, "q.out")
                        if tool in OWN_OUTPUT
                        else (outf or os.path.join(run_dir, "q.out"))
                    )
                    qidx = idx_files(tools[spec["reuse"]] if spec["build"] is None else spec, run_dir)
                    if spec["answers"] == "per-k-mer":
                        # One untimed run writes the per-k-mer answers that the gate reads.
                        # The timed runs discard them in /dev/null, since the counting tools write none.
                        rc_, _, rss, err = bench.once(
                            f"{tag}/{workload} (answers)",
                            cmd,
                            cwd=run_dir,
                            stdout=None if tool == "sbwt" else capt,
                        )
                        query_res = (
                            bench.measure(
                                f"{tag}/{workload}",
                                spec["timed"](query_path) if "timed" in spec else cmd,
                                cwd=run_dir,
                                inputs=[query_path, *qidx],
                            )
                            if rc_ == 0
                            else None
                        )
                    else:
                        query_res = bench.measure(
                            f"{tag}/{workload}", cmd, cwd=run_dir, stdout=capt, inputs=[query_path, *qidx]
                        )
                        rc_, rss, err = (
                            (0 if query_res.succeeded else query_res.code),
                            query_res.rss() or bench.cap / BYTES_PER_MIB,
                            query_res.err,
                        )
                    if rc_ == 0 and not query_res.succeeded:
                        rc_, rss, err = query_res.code, bench.cap / BYTES_PER_MIB, query_res.err
                    if rc_ != 0:
                        kind, det = failure(bench, rc_, rss, err)
                        if kind == "cap":
                            add_status(bench, "cap", **base, phase=workload, gate_detail=det)
                        else:
                            bench.add(
                                **base, phase=workload, gate=False, gate_detail="query " + det
                            )
                        continue
                    got = count_of(tool, outf, capt, err)
                    res.setdefault(workload, {})[tool] = got
                    bench.add(
                        **base,
                        phase=workload,
                        query_s=query_res.wall_row(),
                        query_rss_mb=round(query_res.rss(), 1),
                        kmers=got[0],
                        positive=got[1],
                        gate=True,
                        gate_detail="(cross-tool gate below)",
                        _load=load_timing.get("mean"),
                    )
            # Gates on published or cross-tool distinct counts and on per-workload agreement
            ref_n = kmc_n
            for row in bench.rows[row0:]:
                # Only this round's rows count.
                # FMSI's (threads = 1, no thread option) belong to the first round.
                row_threads = (
                    threads
                    if row.get("tool") not in NO_BUILD_THREADS
                    else (1 if first_t else None)
                )
                if (
                    row.get("dataset"),
                    row.get("k"),
                    row.get("variant"),
                    row.get("mode"),
                    row.get("threads"),
                ) != (ds_id, kmer_len, variant, mode, row_threads) or row.get("section") != "query":
                    continue
                if row.get("status"):
                    continue
                if row.get("phase") == "build" and "index_mb" in row and ref_n:
                    row["bits_per_kmer"] = round(BITS_PER_BYTE * row["index_mb"] * BYTES_PER_MIB / ref_n, 3)
                    if (
                        row.get("kmers_reported") is not None
                        and row["kmers_reported"] != ref_n
                    ):
                        row["gate"] = False
                        row["gate_detail"] = (
                            f"distinct k-mers {row['kmers_reported']} != KMC's {ref_n}"
                        )
                workload = row.get("phase")
                if workload in res:
                    # A tool's positives must equal KMC's, and it must have queried all windows of the query file.
                    # Matching the others' total is not enough, since one tool alone, or a parser that read nothing, agrees with itself.
                    mine, kmc_pos = res[workload].get(row["tool"], (None, None)), kmc_truth.get(workload)
                    passed = mine[1] == kmc_pos and mine[0] == expected[workload]
                    row["kmc_positives"] = kmc_pos
                    row["gate"] = passed
                    row["gate_detail"] = (
                        ""
                        if passed
                        else (
                            f"positives {mine[1]} != KMC {kmc_pos}"
                            if mine[1] != kmc_pos
                            else f"queried {mine[0]} windows of the file's {expected[workload]}"
                        )
                    )
                    if (
                        "mean" in row.get("query_s", {})
                        and row.get("kmers")
                        and "_load" in row
                    ):
                        net = row["query_s"]["mean"] - (row.pop("_load") or 0.0)
                        row["ns_per_kmer_net"] = round(
                            max(net, 0.0) * NS_PER_SECOND / row["kmers"], 1
                        )
                row.pop("_load", None)
                if not input_ok and row.get("tool") != "kmc (input)":
                    # All tools agreed with KMC, but the input is not the paper's set.
                    row["gate"] = False
                    row["gate_detail"] = (
                        (row.get("gate_detail") or "")
                        + f"; the formatted input has {kmc_n} "
                        f"distinct k-mers, the paper publishes {pub}"
                    ).lstrip("; ")
