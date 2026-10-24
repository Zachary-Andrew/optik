"""optik and SSHash k-mer index benchmarks, grouped by mechanism.

optik has one runner per index.

    OptikSpssRunner
        The SPSS index built by ``optik build-spss``, with super-k-mer strings, a
        minimizer MPHF, streaming queries and dense k-mer ids.  It is
        compared with SSHash and SBWT.

    OptikFMRunner
        The FM-index mode, a k-order BWT over the unitigs, built with
        ``optik build-fm`` at several thread counts and queried with
        ``optik query``.  It is a near-minimal exact membership index of a few bits per k-mer, compared with the
        FMSI and KmerCamel FM-index pipeline.

    SshashRunner
        The reference SSHash implementation at m = k/2 (15 for k = 31, the
        value optik-spss uses).  SSHash has no --canonical flag, since a
        k-mer and its reverse complement always share a bucket.

Methodology
    * The headline parity basis is whole-process timing under hyperfine for the build and query phases, the same measurement for optik, sshash and sbwt.
    * A warm in-process basis, labelled warm_*, covers tools with an in-process query harness such as sshash bench.
      There the index is built once and queried inside the measuring process, leaving out subprocess and FASTA parsing overhead.
      These metrics are reported separately and never mixed into the headline parity tables.
    * Correctness gate.
      Each query-timing phase first runs once untimed, and its report must add up (total == positive + negative) before hyperfine timing may run.
"""

import json
import os
import time
from pathlib import Path

from config import PATHS

from ._base import (
    ThreadSweepMixin,
    ToolRunner,
    check_binary,
    ensure_uncompressed,
    gate_query_report,
    generate_positive_reads,
    generate_random_kmers,
    run_cmd,
    unwrap_fasta,
)

# Timeouts of the benchmark phases, in seconds
BUILD_TIMEOUT = 3600
QUERY_TIMEOUT = 300
BENCH_TIMEOUT = 300

# Random independent k-mers for the random-lookup throughput benchmark.
# Each is a uniformly random sequence, so consecutive lookups share no cache locality.
RANDOM_KMER_COUNT = 1_000_000


def _parse_report_field(text, key):
    """Integer field of an sshash-style 'key = value' report, or None when missing or unparseable.

    optik query and sshash both print lines such as 'num_kmers = 12345' or 'num_positive_kmers = 456 (12.34%)'.
    sshash also prints ratios such as 'num_searches = 12/34 (...)', for which the numerator is returned.
    """
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(key):
            try:
                val = line.split("=")[1].strip().split()[0]
                return int(val.split("/")[0])
            except (ValueError, IndexError):
                pass
    return None


def _parse_tsv_report(text):
    """Parse optik's TSV metrics, ``key<TAB>value`` lines, into a dict.

    ``optik build-spss`` and ``optik build-fm`` both print their index metrics this way (build_s, rss_mb, idx_mb, bits_per, distinct_kmers and more), and numeric values become int or float.
    """
    out = {}
    for line in text.split("\n"):
        line = line.strip()
        if "\t" not in line:
            continue
        key, _, val = line.partition("\t")
        val = val.strip()
        try:
            out[key] = int(val) if val.isdigit() else float(val)
        except ValueError:
            out[key] = val
    return out


def _positive_metrics(hf):
    """Nanoseconds per k-mer and positive fraction of a high-hit query run."""
    if "error" in hf:
        return {"pos_error": hf["error"]}
    out = hf["stdout"] + hf["stderr"]
    n = _parse_report_field(out, "num_kmers")
    pos = _parse_report_field(out, "num_positive_kmers")
    res = {"pos_query_s": round(hf["mean_s"], 4), "pos_rss_mb": hf["max_rss_mb"]}
    if n:
        res["pos_ns_per_kmer"] = round(hf["mean_s"] * 1e9 / n, 2)
        res["pos_query_kmers"] = n
        if pos is not None:
            res["pos_positive_frac"] = round(pos / n, 4)
            res["pos_positive"] = pos
    return res


def _run_optik_build_query(ds, ctx, subcommand, label="optik", sweep_build=False):
    """Benchmark one of optik's serialized indexes.

    The index is built with ``optik <subcommand>`` (``build-spss`` or ``build-fm``) and queried with ``optik query``, which recognizes the index format from its magic number.

    With ``sweep_build`` True the build phase runs the auto thread-count sweep with per-thread RSS scaling, and otherwise it runs once at the configured thread count.
    """

    optik_bin = PATHS.get("optik")
    if not optik_bin or not check_binary(optik_bin):
        return {"error": f"optik binary not found at {optik_bin}"}

    class _Helper(ThreadSweepMixin, ToolRunner):
        """Minimal runner that gives this function the thread-sweep machinery."""

        def tool_name(self):
            return label

        def binary_path(self):
            return str(optik_bin)

        def run(self):
            return {}

    _h = _Helper(ds, ctx)
    td = _h.tmpdir(label)
    try:
        unitigs = ensure_uncompressed(ds["unitigs"], td)
        queries = ensure_uncompressed(ds["queries"], td)
        k = int(ds["k"])
        verbose = ctx.get("verbose", False)

        # optik query runs records in parallel on -t threads, and the index format is recognized from the file magic, so one command serves both FM and SPSS indexes.
        # The headline parity number is the query at the default thread count (threads_query, matching the single-threaded sshash query), and a thread sweep reports the scaling.
        idx_path = f"{td}/optik.idx"
        q_default = int(ctx.get("threads_query", 1))

        def query_cmd(threads):
            """optik query command, with -t when threads is positive."""
            cmd = [optik_bin, "query", "-i", idx_path, "-q", queries, "-k", str(k)]
            if threads is not None and threads > 0:
                cmd += ["-t", str(threads)]
            return cmd

        def build_cmd(threads):
            """optik build command for the subcommand, with -t when threads is positive."""
            cmd = [optik_bin, subcommand, "-i", unitigs, "-k", str(k), "-o", idx_path]
            if subcommand == "build-fm":
                cmd += ["-c", "1"]
            # the same -t as SSHash's build, since both build in parallel
            if (
                subcommand in ("build-fm", "build-spss")
                and threads is not None
                and threads > 0
            ):
                cmd += ["-t", str(threads)]
            return cmd

        # One untimed build captures the TSV index metrics and leaves an index on disk for the query phase.
        # Its output is the same at any thread count.
        ts = int(ctx.get("threads_build", 8)) if sweep_build else 1
        proc, _ = run_cmd(build_cmd(ts), timeout=BUILD_TIMEOUT, verbose=verbose, cwd=td)
        if proc.returncode != 0:
            return {"error": f"optik {subcommand} failed: {proc.stderr[:500]}"}
        build_tsv = _parse_tsv_report(proc.stdout + proc.stderr)

        gate_proc, _ = run_cmd(
            query_cmd(q_default), timeout=QUERY_TIMEOUT, verbose=verbose, cwd=td
        )
        if gate_proc.returncode != 0:
            return {"error": f"optik query gate failed: {gate_proc.stderr[:500]}"}
        try:
            gate_query_report(
                gate_proc.stdout + gate_proc.stderr, "optik", f"{subcommand}+query"
            )
        except ValueError as exc:
            return {"error": str(exc)}

        # The build phase is a thread sweep (build-fm) or a single run (build-spss).
        sweep_flat = {}
        if sweep_build:
            sweep_flat, _, hf_head, _hf_head_t = _h.run_sweep(
                "build",
                build_cmd,
                gate=lambda p: p.returncode == 0,
                timeout=BUILD_TIMEOUT,
            )
            if hf_head is not None:
                hf_build = hf_head
            else:
                hf_build = _h.run_hyperfine(
                    build_cmd(int(ctx.get("threads_build", 8))), timeout=BUILD_TIMEOUT
                )
        else:
            hf_build = _h.run_hyperfine(build_cmd(None), timeout=BUILD_TIMEOUT)
        if "error" in hf_build:
            return {"error": f"optik {subcommand} hyperfine: {hf_build['error']}"}

        result = {
            "build_s": round(hf_build["mean_s"], 4),
            "build_s_stddev": round(hf_build["stddev_s"], 4),
            "build_s_min": round(hf_build["min_s"], 4),
            "build_s_max": round(hf_build["max_s"], 4),
            "build_runs": len(hf_build["times_s"]),
            "build_rss_mb": hf_build["max_rss_mb"],
            # Index bytes on disk, the same basis as SSHash's file size.
            # The tool's own in-memory accounting is kept separately.
            "idx_mb": (
                round(os.path.getsize(idx_path) / (1024 * 1024), 4)
                if os.path.exists(idx_path)
                else None
            ),
            "bits_per_kmer": (
                round(os.path.getsize(idx_path) * 8 / build_tsv["distinct_kmers"], 2)
                if os.path.exists(idx_path) and build_tsv.get("distinct_kmers")
                else None
            ),
            "idx_mb_inmem": build_tsv.get("idx_mb"),
            "distinct_kmers": build_tsv.get("distinct_kmers"),
        }
        for extra_key in ("unitigs", "super_len", "bwt_len"):
            if extra_key in build_tsv:
                result[extra_key] = build_tsv[extra_key]
        result.update(sweep_flat)

        # The streaming query runs a -t thread sweep like the build phase, after the correctness gate above.
        # The sweep's headline is the hyperfine run at the default thread count (threads_query), which gives the parity metric next to the single-threaded sshash query.
        query_flat, _query_rows, hf_query, hf_query_t = _h.run_sweep(
            "query", query_cmd, timeout=QUERY_TIMEOUT
        )
        if hf_query is None:
            hf_query = _h.run_hyperfine(query_cmd(q_default), timeout=QUERY_TIMEOUT)
            hf_query_t = q_default
        if "error" in hf_query:
            result["error_query"] = hf_query["error"]
            return result
        result.update(query_flat)
        result["query_threads"] = hf_query_t or q_default

        result["query_s"] = round(hf_query["mean_s"], 4)
        result["query_s_stddev"] = round(hf_query["stddev_s"], 4)
        result["query_s_min"] = round(hf_query["min_s"], 4)
        result["query_s_max"] = round(hf_query["max_s"], 4)
        result["query_runs"] = len(hf_query["times_s"])
        result["query_rss_mb"] = hf_query["max_rss_mb"]

        q_out = hf_query["stdout"] + hf_query["stderr"]
        num_kmers = _parse_report_field(q_out, "num_kmers")
        num_positive = _parse_report_field(q_out, "num_positive_kmers")
        num_negative = _parse_report_field(q_out, "num_negative_kmers")
        num_searches = _parse_report_field(q_out, "num_searches")

        if num_kmers:
            result["ns_per_query"] = round((hf_query["mean_s"] * 1e9) / num_kmers, 2)
            result["queries_s"] = round(num_kmers / hf_query["mean_s"], 2)
            result["query_kmers"] = num_kmers
        result["query_positive"] = num_positive
        result["query_negative"] = num_negative
        result["num_searches"] = num_searches

        # Throughput of random independent k-mer lookups, timed by hyperfine
        rand_path = f"{td}/random_kmers.fa"
        generate_random_kmers(rand_path, k, RANDOM_KMER_COUNT)

        def query_cmd_rand(threads):
            """optik query command over the random k-mer file."""
            return [
                optik_bin,
                "query",
                "-i",
                idx_path,
                "-q",
                rand_path,
                "-k",
                str(k),
            ] + (["-t", str(threads)] if threads and threads > 0 else [])

        hf_rand = _h.run_hyperfine(query_cmd_rand(q_default), timeout=QUERY_TIMEOUT)
        if "error" not in hf_rand:
            result["rand_ns"] = round((hf_rand["mean_s"] * 1e9) / RANDOM_KMER_COUNT, 2)
            result["rand_rss_mb"] = hf_rand["max_rss_mb"]
            rand_out = hf_rand["stdout"] + hf_rand["stderr"]
            result["rand_positive"] = _parse_report_field(
                rand_out, "num_positive_kmers"
            )
            result["rand_negative"] = _parse_report_field(
                rand_out, "num_negative_kmers"
            )

        # High-hit workload, made of reads sampled from the indexed sequences
        pos_path = generate_positive_reads(unitigs, f"{td}/positive_reads.fa")
        hf_pos = _h.run_hyperfine(
            [optik_bin, "query", "-i", idx_path, "-q", pos_path, "-k", str(k)]
            + (["-t", str(q_default)] if q_default and q_default > 0 else []),
            timeout=QUERY_TIMEOUT,
        )
        result.update(_positive_metrics(hf_pos))

        return result

    finally:
        _h.cleanup(td)


class SshashRunner(ThreadSweepMixin, ToolRunner):
    """Reference SSHash implementation at m = k/2, optik-spss's default.

    SSHash builds a full dictionary supporting lookup, access and membership queries, and is compared with optik-spss.
    Its build accepts ``-t`` and is timed at each thread count like optik's (ThreadSweepMixin), while sshash query is single-threaded.
    """

    @classmethod
    def tool_name(cls):
        return "sshash"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("sshash")
        if p and check_binary(p):
            return str(p)
        return None

    @property
    def _sshash(self):
        """Path of the SSHash binary."""
        return self.binary_path()

    def run(self):
        td = self.tmpdir("sshash")
        try:
            unitigs = ensure_uncompressed(self.ds["unitigs"], td)
            queries = ensure_uncompressed(self.ds["queries"], td)
            # Without --multiline SSHash reads one sequence per line, so wrapped FASTA queries are unwrapped to one line per record.
            with open(queries) as _fh:
                if _fh.read(1) == ">":
                    queries = unwrap_fasta(queries, td)
            k = int(self.ds["k"])
            m = max(k // 2, 1)
            verbose = self.ctx.get("verbose", False)

            idx_path = f"{td}/index.sshash"
            default_t = int(self.ctx.get("threads_build", 8))

            def build_cmd(t):
                """argv of sshash build with t threads."""
                return [
                    self._sshash,
                    "build",
                    "-i",
                    unitigs,
                    "-k",
                    str(k),
                    "-m",
                    str(m),
                    "-o",
                    idx_path,
                    "-t",
                    str(t if t else default_t),
                    "-d",
                    td,  # scratch files go to the run directory, not the cwd
                ]

            # The build phase is the auto thread sweep.
            # Its first point runs once untimed as the gate, and that run is captured to parse the distinct k-mer count and to leave an index on disk for the queries.
            gate_out = {}

            def gate(p):
                """Keep the gate run's output and pass when it exited 0."""
                if p.returncode == 0:
                    gate_out["text"] = p.stdout + p.stderr
                return p.returncode == 0

            sweep_flat, _, hf_build, _hf_build_t = self.run_sweep(
                "build", build_cmd, gate=gate, timeout=BUILD_TIMEOUT
            )
            if hf_build is None:
                hf_build = self.run_hyperfine(build_cmd(default_t), cwd=td)
            if "error" in hf_build:
                return {"error": f"sshash build hyperfine: {hf_build['error']}"}

            # The distinct k-mer count comes from one untimed --verbose build, kept out of the timed runs so that SSHash's timings carry no logging overhead.
            pv, _ = run_cmd(
                build_cmd(default_t) + ["--verbose"],
                timeout=BUILD_TIMEOUT,
                verbose=verbose,
                cwd=td,
            )
            build_out = pv.stdout + pv.stderr
            distinct_kmers = None
            for line in build_out.split("\n"):
                parts = line.split()
                # The count is on the 'num_kmers N' line of a verbose build, while lines such as 'read S sequences, B bases, K kmers' lead with other counts.
                if len(parts) == 2 and parts[0] == "num_kmers" and parts[1].isdigit():
                    distinct_kmers = int(parts[1])
                    break
            idx_bytes = (
                Path(idx_path).stat().st_size if Path(idx_path).exists() else None
            )

            result = {
                "build_s": round(hf_build["mean_s"], 4),
                "build_s_stddev": round(hf_build["stddev_s"], 4),
                "build_s_min": round(hf_build["min_s"], 4),
                "build_s_max": round(hf_build["max_s"], 4),
                "build_runs": len(hf_build["times_s"]),
                "build_rss_mb": hf_build["max_rss_mb"],
                "idx_mb": round(idx_bytes / (1024 * 1024), 4) if idx_bytes else None,
                "bits_per_kmer": round((idx_bytes * 8) / distinct_kmers, 2)
                if (idx_bytes and distinct_kmers)
                else None,
                "distinct_kmers": distinct_kmers,
            }
            result.update(sweep_flat)

            bench_cmd = [self._sshash, "bench", "-i", idx_path]
            t0 = time.time()
            proc_b, _ = run_cmd(bench_cmd, timeout=BENCH_TIMEOUT, verbose=verbose)
            result["warm_bench_s"] = round(time.time() - t0, 4)

            bench_out = proc_b.stdout + proc_b.stderr
            for line in reversed(bench_out.split("\n")):
                line_s = line.strip()
                if line_s.startswith("{") and line_s.endswith("}"):
                    try:
                        bench_json = json.loads(line_s)
                        for key, val in bench_json.items():
                            kl = key.lower()
                            if "positive lookup" in kl:
                                result["warm_rand_bench_ns"] = round(float(val), 2)
                            elif "negative lookup" in kl:
                                result["warm_rand_neg_bench_ns"] = round(float(val), 2)
                            elif "access" in kl:
                                result["warm_access_bench_ns"] = round(float(val), 2)
                    except (json.JSONDecodeError, ValueError, TypeError):
                        pass
                    break

            query_cmd = [
                self._sshash,
                "query",
                "-i",
                idx_path,
                "-q",
                queries,
                # Not --multiline, which concatenates records and so counts k-mers spanning two reads, which are invalid, inflating num_kmers by about (k-1)/read_len and understating SSHash's ns/k-mer.
            ]
            # Correctness gate, one untimed run whose report must add up (total == positive + negative) before hyperfine timing may run.
            gate_proc, _ = run_cmd(
                query_cmd, timeout=QUERY_TIMEOUT, verbose=verbose, cwd=td
            )
            if gate_proc.returncode != 0:
                return {"error": f"sshash query gate failed: {gate_proc.stderr[:500]}"}
            try:
                gate_query_report(
                    gate_proc.stdout + gate_proc.stderr, "sshash", "query"
                )
            except ValueError as exc:
                return {"error": str(exc)}
            hf_query = self.run_hyperfine(query_cmd, timeout=QUERY_TIMEOUT)
            if "error" in hf_query:
                result["error_query"] = hf_query["error"]
                return result

            result["query_s"] = round(hf_query["mean_s"], 4)
            result["query_s_stddev"] = round(hf_query["stddev_s"], 4)
            result["query_s_min"] = round(hf_query["min_s"], 4)
            result["query_s_max"] = round(hf_query["max_s"], 4)
            result["query_runs"] = len(hf_query["times_s"])
            result["query_rss_mb"] = hf_query["max_rss_mb"]

            q_out = hf_query["stdout"] + hf_query["stderr"]
            num_kmers = _parse_report_field(q_out, "num_kmers")
            num_positive = _parse_report_field(q_out, "num_positive_kmers")
            num_negative = _parse_report_field(q_out, "num_negative_kmers")
            num_searches = _parse_report_field(q_out, "num_searches")

            if num_kmers:
                result["ns_per_query"] = round(
                    (hf_query["mean_s"] * 1e9) / num_kmers, 2
                )
                result["queries_s"] = round(num_kmers / hf_query["mean_s"], 2)
                result["query_kmers"] = num_kmers
                result["query_positive"] = num_positive
                result["query_negative"] = num_negative
                if num_searches is not None:
                    result["num_searches"] = num_searches

            rand_path = f"{td}/random_kmers.fa"
            generate_random_kmers(rand_path, k, RANDOM_KMER_COUNT)

            rand_cmd = [
                self._sshash,
                "query",
                "-i",
                idx_path,
                "-q",
                rand_path,
            ]
            hf_rand = self.run_hyperfine(rand_cmd, timeout=QUERY_TIMEOUT)
            if "error" not in hf_rand:
                result["rand_ns"] = round(
                    (hf_rand["mean_s"] * 1e9) / RANDOM_KMER_COUNT, 2
                )
                result["rand_rss_mb"] = hf_rand["max_rss_mb"]
                rand_out = hf_rand["stdout"] + hf_rand["stderr"]
                result["rand_positive"] = _parse_report_field(
                    rand_out, "num_positive_kmers"
                )
                result["rand_negative"] = _parse_report_field(
                    rand_out, "num_negative_kmers"
                )

            # High-hit workload, reads sampled from the indexed sequences
            pos_path = generate_positive_reads(unitigs, f"{td}/positive_reads.fa")
            hf_pos = self.run_hyperfine(
                [self._sshash, "query", "-i", idx_path, "-q", pos_path],
                timeout=QUERY_TIMEOUT,
            )
            result.update(_positive_metrics(hf_pos))

            return result

        finally:
            self.cleanup(td)


class OptikSpssRunner(ToolRunner):
    """optik's static SPSS core (build-spss), with super-k-mer strings, a minimizer MPHF, one record per bucket and streaming queries.

    This is the direct counterpart of SSHash, with the same input (unitigs), the same m = k/2 and a single-threaded query.
    """

    @classmethod
    def tool_name(cls):
        return "optik-spss"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("optik")
        if p and check_binary(p):
            return str(p)
        return None

    def run(self):
        return _run_optik_build_query(
            self.ds,
            self.ctx,
            subcommand="build-spss",
            label="optik_spss",
            sweep_build=True,
        )


class OptikFMRunner(ToolRunner):
    """optik's FM membership index (build-fm, a k-order BWT over the unitigs).

    ``optik build-fm`` builds a near-minimal exact FM index (0 false positives, 0 false negatives) of a few bits per k-mer, taking ``-t`` for the unitig-walk phase, and ``optik query`` queries it after recognizing the FM magic.
    Its build also runs at every sweep point, and the index is compared with the FMSI and KmerCamel FM-index pipeline.
    """

    @classmethod
    def tool_name(cls):
        return "optik-fm"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("optik")
        if p and check_binary(p):
            return str(p)
        return None

    def run(self):
        return _run_optik_build_query(
            self.ds, self.ctx, subcommand="build-fm", label="optik_fm", sweep_build=True
        )
