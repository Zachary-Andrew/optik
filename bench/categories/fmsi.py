"""FM-index k-mer membership benchmarks of the FMSI pipeline.

Each runner implements the pipeline's three phases.
    1. kmercamel compute   builds a masked superstring from the unitigs
    2. fmsi index          builds an FM-index of the masked superstring
    3. fmsi query          queries k-mers against the FM-index

Phase 1 (kmercamel) runs once, since it is deterministic, and its wall time is added to the hyperfine-timed `fmsi index`.
build_s therefore covers the same input-to-index path as every other tool, with kmercamel_s reporting the superstring's share.
Queries use the read file the other runners use.
"""

from pathlib import Path

from config import PATHS

from ._base import (
    ToolRunner,
    check_binary,
    count_kmers,
    ensure_uncompressed,
    generate_random_kmers,
    run_cmd,
    write_acgt_segments,
)

# Timeouts of the benchmark phases, in seconds
BUILD_TIMEOUT = 3600
QUERY_TIMEOUT = 300

RANDOM_KMERS = 1_000_000  # the random lookup workload, the same size as for the other k-mer runners


class FmsiRunner(ToolRunner):
    """Benchmark runner for the FMSI pipeline of kmercamel then fmsi.

    KmerCamel computes a masked superstring from the input unitigs, and fmsi builds an FM-index on it and answers membership queries.
    The build metric is the timed fmsi index step plus KmerCamel's time, which is also reported on its own as kmercamel_s.

    FMSI is described by Sladký, Veselý and Břinda (2025, Bioinformatics Advances).
    """

    @classmethod
    def tool_name(cls):
        return "fmsi"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("fmsi")
        if p and check_binary(p):
            return str(p)
        return None

    def _kmercamel_bin(self):
        """Path of the kmercamel binary that computes the masked superstring."""
        p = PATHS.get("kmercamel")
        if p and check_binary(p):
            return str(p)
        return None

    @property
    def _fmsi(self):
        """Path of the fmsi binary."""
        return self.binary_path()

    def run(self):
        kmc_bin = self._kmercamel_bin()
        if not kmc_bin:
            return {"error": "kmercamel not found (required for FMSI pipeline)"}

        td = self.tmpdir("fmsi")
        try:
            unitigs = ensure_uncompressed(self.ds["unitigs"], td)
            queries = ensure_uncompressed(self.ds["queries"], td)
            k = int(self.ds["k"])
            verbose = self.ctx.get("verbose", False)

            # Path of the masked superstring
            ms_path = f"{td}/ms.fa"

            # Phase 1 computes the masked superstring with kmercamel.
            # It runs once outside hyperfine, because this preprocessing step produces the same output on every run, and only fmsi index is timed.
            # -S declares the input to be unitigs, KmerCamel's documented fast path for precomputed unitigs or simplitigs and FMSI's own workflow.
            # Without it the superstring step took 1.8x longer on chr21 for the same answers.
            cmd_kmc = [
                kmc_bin,
                "compute",
                "-k",
                str(k),
                "-S",
                "-M",
                ms_path,
                "-o",
                "/dev/null",
                unitigs,
            ]
            proc_kmc, kmc_s = run_cmd(cmd_kmc, timeout=BUILD_TIMEOUT, verbose=verbose)
            if proc_kmc.returncode != 0:
                return {"error": f"kmercamel failed: {proc_kmc.stderr[:500]}"}
            if not Path(ms_path).exists():
                return {"error": "kmercamel did not produce masked superstring"}

            # Distinct k-mers of the indexed set, counted from the unitigs.
            # KmerCamel's log reports the count only without -S, so parsing it silently dropped FMSI's bits/k-mer.
            distinct_kmers = count_kmers(unitigs, k, distinct=True)

            # Phase 2 runs fmsi index once to collect the index metrics.
            cmd_idx = [self._fmsi, "index", "-k", str(k), ms_path]
            proc_idx, _ = run_cmd(cmd_idx, timeout=BUILD_TIMEOUT, verbose=verbose)
            if proc_idx.returncode != 0:
                return {"error": f"fmsi index failed: {proc_idx.stderr[:500]}"}

            # The index size sums every fmsi component file.
            fmsi_components = list(Path(td).glob("*.fmsi.*"))
            idx_bytes = sum(f.stat().st_size for f in fmsi_components)

            # Hyperfine times fmsi index, and kmercamel's time is added below.
            hf_build = self.run_hyperfine(cmd_idx, timeout=BUILD_TIMEOUT)
            if "error" in hf_build:
                return {"error": f"fmsi build hyperfine: {hf_build['error']}"}

            # Phase 3 queries the read file every other tool queries.
            # FMSI streams whole records, whereas it once got one record per k-mer, cut to 50,000, a different workload.
            # For a record containing a non-ACGT base FMSI answers one character per base rather than per k-mer window, so its bits cannot be attributed to windows.
            # It is therefore given each maximal ACGT run as a record, which are exactly the valid windows of the read file.
            seg_q = f"{td}/query_segments.fa"
            n_valid, n_invalid = write_acgt_segments(queries, k, seg_q)
            cmd_q = [
                self._fmsi,
                "query",
                "-q",
                seg_q,
                "-k",
                str(k),
                ms_path,
            ]
            hf_query = self.run_hyperfine(cmd_q, timeout=QUERY_TIMEOUT)
            if "error" in hf_query:
                return {
                    "build_s": round(hf_build["mean_s"], 4),
                    "build_s_stddev": round(hf_build["stddev_s"], 4),
                    "idx_mb": round(idx_bytes / (1024 * 1024), 4)
                    if idx_bytes
                    else None,
                    "bits_per_kmer": None,
                    "distinct_kmers": distinct_kmers,
                    "error_query": hf_query["error"],
                }

            # Each output line holds a bitstring of present and absent flags for its record, so the 1s count as positive and the length as total.
            q_lines = hf_query["stdout"].strip().split("\n")
            total_present = 0
            total_kmers = 0
            for line in q_lines:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    bitstring = parts[-1]
                    total_present += bitstring.count("1")
                    total_kmers += len(bitstring)

            # The build RSS is the larger of the two pipeline phases.
            rss_vals = [
                v
                for v in [hf_build["max_rss_mb"], hf_query["max_rss_mb"]]
                if v is not None
            ]
            pipeline_rss = max(rss_vals) if rss_vals else None

            # FMSI indexes a masked superstring, so its build from the unitigs every other tool starts from is KmerCamel plus fmsi index.
            result = {
                "build_s": round(hf_build["mean_s"] + kmc_s, 4),
                "build_s_stddev": round(hf_build["stddev_s"], 4),
                "build_s_min": round(hf_build["min_s"] + kmc_s, 4),
                "build_s_max": round(hf_build["max_s"] + kmc_s, 4),
                "build_runs": len(hf_build["times_s"]),
                "build_rss_mb": pipeline_rss,
                "kmercamel_s": round(kmc_s, 2),
                "idx_mb": round(idx_bytes / (1024 * 1024), 4) if idx_bytes else None,
                "distinct_kmers": distinct_kmers,
                "query_s": round(hf_query["mean_s"], 4),
                "query_s_stddev": round(hf_query["stddev_s"], 4),
                "query_s_min": round(hf_query["min_s"], 4),
                "query_s_max": round(hf_query["max_s"], 4),
                "query_runs": len(hf_query["times_s"]),
                "query_rss_mb": hf_query["max_rss_mb"],
            }

            if total_kmers:
                result["ns_per_query"] = round(
                    (hf_query["mean_s"] * 1e9) / total_kmers, 2
                )
                result["queries_s"] = round(total_kmers / hf_query["mean_s"], 2)
                result["queries_s_max"] = (
                    round(total_kmers / hf_query["min_s"], 2)
                    if hf_query["min_s"] > 0
                    else None
                )
                if total_kmers != n_valid:
                    return {
                        **result,
                        "error": f"fmsi answered {total_kmers} k-mers for "
                        f"{n_valid} valid query windows",
                    }
                result["query_kmers"] = n_valid + n_invalid
                result["query_positive"] = total_present
                result["query_invalid"] = n_invalid
                result["query_negative"] = n_valid - total_present

            if idx_bytes and distinct_kmers:
                result["bits_per_kmer"] = round((idx_bytes * 8) / distinct_kmers, 2)

            # Phase 4 measures random k-mer lookup throughput with hyperfine.
            rand_path = f"{td}/random_kmers.fa"
            generate_random_kmers(rand_path, k, RANDOM_KMERS)
            cmd_rand = [
                self._fmsi,
                "query",
                "-q",
                rand_path,
                "-k",
                str(k),
                ms_path,
            ]
            hf_rand = self.run_hyperfine(cmd_rand, timeout=QUERY_TIMEOUT)
            if "error" not in hf_rand:
                rand_kmers = 0
                rand_present = 0
                for line in hf_rand["stdout"].strip().split("\n"):
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        bitstring = parts[-1]
                        rand_present += bitstring.count("1")
                        rand_kmers += len(bitstring)
                result["rand_ns"] = (
                    round((hf_rand["mean_s"] * 1e9) / rand_kmers, 2)
                    if rand_kmers
                    else None
                )
                result["rand_rss_mb"] = hf_rand["max_rss_mb"]
                result["rand_positive"] = rand_present
                result["rand_negative"] = (
                    rand_kmers - rand_present if rand_kmers else None
                )

            return result

        finally:
            self.cleanup(td)
