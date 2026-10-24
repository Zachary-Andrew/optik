"""SBWT (Spectral Burrows-Wheeler Transform) k-mer membership benchmarks.

SBWT is a succinct data structure that represents a k-mer set through subset rank queries.
Two variants are benchmarked.

    sbwt-fast   the plain-matrix variant, with the faster queries
    sbwt-small  the rrr-split variant, with the smaller index

Both run one binary with different --variant flags and share one benchmark implementation (_SbwtBase), so their metrics and column order match.

Metrics collected
    Build wall time (hyperfine), index size, bits per k-mer, peak RSS,
    sequential query throughput (ns per k-mer, queries per second),
    positive and negative counts, and random k-mer lookup throughput.
"""

from pathlib import Path

from config import PATHS

from ._base import (
    ThreadSweepMixin,
    ToolRunner,
    check_binary,
    count_invalid_windows,
    ensure_uncompressed,
    generate_random_kmers,
    run_cmd,
)

BUILD_TIMEOUT = 3600
QUERY_TIMEOUT = 300
RANDOM_KMER_COUNT = 1_000_000


class _SbwtBase(ThreadSweepMixin, ToolRunner):
    """Implementation shared by both SBWT benchmark runners.

    `sbwt build` takes ``-t/--n-threads`` and is timed per thread count (ThreadSweepMixin), while `sbwt search` is single-threaded and gets no sweep.
    """

    @property
    def _sbwt(self):
        """Path of the SBWT binary."""
        return self.binary_path()

    def _build_cmd(self, unitigs, k, variant, t):
        """sbwt build command for the variant on t threads, with reverse complements added."""
        cmd = [
            self._sbwt,
            "build",
            "-i",
            unitigs,
            "-o",
            f"{self._td}/index.sbwt",
            "-k",
            str(k),
            "--variant",
            variant,
            "-t",
            str(t if t else int(self.ctx.get("threads_build", 8))),
        ]
        # Both variants also index reverse complements, because the queries ask about canonical membership (every other tool answers a k-mer and its reverse complement alike).
        # Without the flag the index holds only forward k-mers, and reverse-strand query k-mers come out absent.
        cmd += ["--add-reverse-complements", "--temp-dir", self._td]
        return cmd

    def _parse_sbwt_output(self, results_path):
        """Parse SBWT search output into (total, positive) k-mer counts.

        SBWT writes one line per input read, a space-separated list of k-mer ids in which -1 marks an absent k-mer and any other value a present one.
        Returns (None, None) when the file is missing or empty.
        """
        p = Path(results_path)
        if not p.exists():
            return None, None
        total = 0
        positive = 0
        for line in p.read_text().split("\n"):
            line = line.strip()
            if not line:
                continue
            for tok in line.split():
                total += 1
                if tok != "-1":
                    positive += 1
        if total == 0:
            return None, None
        return total, positive

    def _count_unitig_kmers(self, unitigs_path, k):
        """k-mer windows of the input unitigs, used for bits per k-mer when the tool reports no count of its own.
        Unitigs hold each k-mer once, so this is the distinct k-mer count."""
        try:
            import gzip

            opener = gzip.open if str(unitigs_path).endswith(".gz") else open
            total = 0
            with opener(str(unitigs_path), "rt") as f:
                seq_len = 0
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if seq_len > 0:
                            total += max(0, seq_len - k + 1)
                        seq_len = 0
                    elif line:
                        seq_len += len(line)
                if seq_len > 0:
                    total += max(0, seq_len - k + 1)
            return total if total > 0 else None
        except (OSError, ValueError, EOFError):
            # An unreadable or truncated unitig file leaves the count unknown.
            return None

    def _run_variant(self, td, unitigs, queries, k, verbose, variant):
        """Build, gate and time one SBWT variant and return its metrics."""
        self._td = td
        idx_path = f"{td}/index.sbwt"

        # Each thread count of the sweep times the build.
        # The gate, untimed before the first point, also leaves the index on disk.
        # `sbwt build` takes -t/--n-threads for parallel construction.
        def build_cmd(t):
            """Command line of sbwt build using t threads."""
            return self._build_cmd(unitigs, k, variant, t)

        sweep_flat, _, hf_build, _ = self.run_sweep(
            "build",
            build_cmd,
            gate=lambda p: p.returncode == 0,
            timeout=BUILD_TIMEOUT,
            cwd=td,
        )
        if hf_build is None:
            hf_build = self.run_hyperfine(
                build_cmd(None), timeout=BUILD_TIMEOUT, cwd=td
            )
        if "error" in hf_build:
            return {"error": f"sbwt build hyperfine: {hf_build['error']}"}

        idx_bytes = sum(
            f.stat().st_size
            for f in Path(td).iterdir()
            if f.is_file() and "index" in f.name
        )
        distinct_kmers = self._count_unitig_kmers(unitigs, k)

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

        query_cmd = [
            self._sbwt,
            "search",
            "-i",
            idx_path,
            "-q",
            queries,
            "-o",
            f"{td}/results.txt",
        ]
        # Correctness gate, one untimed run whose search must succeed and write a parseable result file before hyperfine timing starts.
        gate_proc, _ = run_cmd(
            query_cmd, timeout=QUERY_TIMEOUT, verbose=verbose, cwd=td
        )
        if gate_proc.returncode != 0:
            return {"error": f"sbwt search gate failed: {gate_proc.stderr[:500]}"}
        g_total, _ = self._parse_sbwt_output(f"{td}/results.txt")
        if not g_total:
            return {
                "error": "sbwt search gate failed: search produced no k-mer results"
            }
        hf_query = self.run_hyperfine(
            query_cmd, timeout=QUERY_TIMEOUT, cwd=td, ignore_failure=True
        )
        if "error" in hf_query:
            result["error_query"] = hf_query["error"]
            return result

        result["query_s"] = round(hf_query["mean_s"], 4)
        result["query_s_stddev"] = round(hf_query["stddev_s"], 4)
        result["query_s_min"] = round(hf_query["min_s"], 4)
        result["query_s_max"] = round(hf_query["max_s"], 4)
        result["query_runs"] = len(hf_query["times_s"])
        result["query_rss_mb"] = hf_query["max_rss_mb"]

        total, positive = self._parse_sbwt_output(f"{td}/results.txt")
        if total:
            invalid = count_invalid_windows(queries, k)
            result["query_kmers"] = total
            result["query_positive"] = positive
            result["query_invalid"] = invalid
            result["query_negative"] = total - positive - invalid
            result["ns_per_query"] = round((hf_query["mean_s"] * 1e9) / total, 2)
            result["queries_s"] = round(total / hf_query["mean_s"], 2)

        rand_path = f"{td}/random_kmers.fa"
        generate_random_kmers(rand_path, k, RANDOM_KMER_COUNT)
        rand_cmd = [
            self._sbwt,
            "search",
            "-i",
            idx_path,
            "-q",
            rand_path,
            "-o",
            f"{td}/rand_results.txt",
        ]
        hf_rand = self.run_hyperfine(
            rand_cmd, timeout=QUERY_TIMEOUT, cwd=td, ignore_failure=True
        )
        if "error" not in hf_rand:
            result["rand_ns"] = round((hf_rand["mean_s"] * 1e9) / RANDOM_KMER_COUNT, 2)
            result["rand_rss_mb"] = hf_rand["max_rss_mb"]
            rand_total, rand_positive = self._parse_sbwt_output(
                f"{td}/rand_results.txt"
            )
            if rand_total:
                result["rand_positive"] = rand_positive
                result["rand_negative"] = rand_total - rand_positive

        return result


class SbwtFastRunner(_SbwtBase):
    """Benchmark runner for SBWT's plain-matrix variant, built for fast queries.

    --variant plain-matrix gives the fastest subset rank queries at the cost
    of a larger index.
    """

    @classmethod
    def tool_name(cls):
        return "sbwt-fast"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("sbwt")
        if p and check_binary(p):
            return str(p)
        return None

    def run(self):
        td = self.tmpdir("sbwt_fast")
        try:
            unitigs = ensure_uncompressed(self.ds["unitigs"], td)
            queries = ensure_uncompressed(
                self.ds["queries"], td, strip_blank_lines=True
            )
            k = int(self.ds["k"])
            verbose = self.ctx.get("verbose", False)
            return self._run_variant(td, unitigs, queries, k, verbose, "plain-matrix")
        finally:
            self.cleanup(td)


class SbwtSmallRunner(_SbwtBase):
    """Benchmark runner for SBWT's rrr-split variant, built for a small index.

    --variant rrr-split compresses with RRR for a smaller index at the cost of
    slower queries.
    """

    @classmethod
    def tool_name(cls):
        return "sbwt-small"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("sbwt")
        if p and check_binary(p):
            return str(p)
        return None

    def run(self):
        td = self.tmpdir("sbwt_small")
        try:
            unitigs = ensure_uncompressed(self.ds["unitigs"], td)
            queries = ensure_uncompressed(
                self.ds["queries"], td, strip_blank_lines=True
            )
            k = int(self.ds["k"])
            verbose = self.ctx.get("verbose", False)
            return self._run_variant(td, unitigs, queries, k, verbose, "rrr-split")
        finally:
            self.cleanup(td)
