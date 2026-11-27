"""De Bruijn graph construction benchmarks.

Tools in this category

    bifrost
        Colored compacted de Bruijn graph builder.
        It builds from reads or FASTA in its -s sequence-file mode and writes a gzipped GFA of unitigs, from the same input optik gets.

    cuttlefish
        Graph compactor that writes the maximal unitigs of the graph, run in --read mode on FASTQ or FASTA, again from optik's input.

    optik
        optik's dbruijn subcommand counts canonical k-mers through minimizer-partitioned super-k-mers into one open-addressing node table and compacts maximal unitigs.
        It writes GFA 1.0 with the summed abundance KC:i, the mean abundance km:f and one L line per edge.

All three tools build from the dataset's dbg_input, so construction speed and output can be compared directly.

Metrics collected
    Build wall time from hyperfine, graph size in MB, peak RSS, the unitig count, total length and N50, and a cross-tool correctness comparison.

Correctness checks
    Each tool's graph output is parsed for its number of unitigs, its total sequence length and, where it applies, GFA validity.
    Across tools the unitig count ratio and the total length ratio between optik and the reference tools are reported, which show whether optik's graph is consistent with theirs.
"""

import json
from pathlib import Path

from config import PATHS
from diagnostics import DETAIL_CHARS
from dna import n50
from sam_spec import GfaSegmentCol
from units import BYTES_PER_MIB

from ._base import (
    DEFAULT_BUILD_THREADS,
    ThreadSweepMixin,
    ToolRunner,
    check_binary,
    count_kmers,
    ensure_uncompressed,
    run_cmd,
)

# Timeout of graph construction, in seconds
BUILD_TIMEOUT = 3600


# Every dBG tool gets the same uncompressed input, abundance cutoff and thread count, in the mode that matches optik's edge-centric graph where the tool has one.
# Each tool's output k-mer set is dumped so that tool_bench.py can require identical sets (graph_parity).


# Graph builds keep k-mers seen at least twice unless the dataset sets its own cutoff.
DEFAULT_DBG_MIN_COUNT = 2
# Bifrost's two filters, -r keeping every k-mer and -s dropping those seen once, by the cutoff each gives.
BIFROST_CUTOFF_FLAGS = {1: "-r", 2: "-s"}


def dbg_setup(runner, tmp_dir):
    """Input path, k, cutoff and build threads that every graph runner of a dataset shares."""
    dataset, ctx = runner.dataset, runner.ctx
    src = dataset.get("dbg_input", dataset["queries"])
    inp = ensure_uncompressed(src, tmp_dir)  # decompression is never timed
    kmer_len = int(dataset["k"])
    min_count = int(dataset.get("dbg_min_count", DEFAULT_DBG_MIN_COUNT))
    threads = int(ctx.get("threads_build", DEFAULT_BUILD_THREADS))
    return inp, kmer_len, min_count, threads


_CANON_TBL = str.maketrans("ACGTacgt", "TGCATGCA")


def _iter_seqs(path):
    """Yield the sequences of a FASTA or FASTQ file, or the S lines of a GFA, plain or gzipped."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("S\t"):
                yield line.split("\t")[GfaSegmentCol.SEQ].strip()
            elif line and line[0] not in ">@HLP+#" and not line.startswith("\t"):
                yield line.strip()


def dump_graph_kmers(runner, tool, paths, kmer_len):
    """Count the canonical k-mers of a tool's FASTA or GFA unitig output, possibly gzipped.

    Returns (distinct, occurrences, dump path).
    occurrences above distinct means some k-mer appears at more than one unitig position."""

    seen = set()
    occ = 0
    for path in paths:
        for seq in _iter_seqs(path):
            seq = seq.upper()
            seq_rc = seq.translate(_CANON_TBL)[::-1]
            seq_len = len(seq)
            for i in range(seq_len - kmer_len + 1):
                fwd_kmer = seq[i : i + kmer_len]
                rev_kmer = seq_rc[seq_len - kmer_len - i : seq_len - i]
                if "N" in fwd_kmer:
                    continue
                seen.add(min(rev_kmer, fwd_kmer))
                occ += 1
    out_dir = Path(runner.ctx.get("tmp_dir", "/tmp/tool_bench")) / "graph_parity"
    out_dir.mkdir(parents=True, exist_ok=True)
    dump = out_dir / f"{tool}.kmers"
    with open(dump, "w") as handle:
        handle.write("\n".join(sorted(seen)))
    return len(seen), occ, str(dump)


def _parse_gfa_stats(gfa_path):
    """Graph statistics of a GFA file, or None when it is missing or empty.

    The dict holds
        num_unitigs   number of S (segment) lines
        total_length  sum of unitig sequence lengths
        max_unitig    length of the longest unitig
        avg_unitig    mean unitig length
        n50           N50 unitig length
        valid_gfa     whether the file starts with an H VN:Z:1.0 header
    """
    path = Path(gfa_path)
    if not path.exists() or path.stat().st_size == 0:
        return None

    lengths = []
    valid_gfa = False
    try:
        with open(path, "r") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("H\t"):
                    # GFA 1.0 header
                    if "VN:Z:1.0" in line:
                        valid_gfa = True
                elif line.startswith("S\t"):
                    parts = line.split("\t")
                    if len(parts) > GfaSegmentCol.SEQ:
                        seq = parts[GfaSegmentCol.SEQ]
                        lengths.append(len(seq))
    except (OSError, UnicodeDecodeError):
        return None

    if not lengths:
        return None

    total = sum(lengths)
    unitig_n50 = n50(lengths)

    return {
        "num_unitigs": len(lengths),
        "total_length": total,
        "max_unitig": max(lengths),
        "avg_unitig": round(total / len(lengths), 1),
        "min_unitig": min(lengths),
        "n50": unitig_n50,
        "valid_gfa": valid_gfa,
    }


def _parse_gfa_gz_stats(gfa_gz_path):
    """Graph statistics of a gzipped GFA file (.gfa.gz), the form Bifrost writes."""
    import gzip

    lengths = []
    valid_gfa = False
    try:
        with gzip.open(gfa_gz_path, "rt") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                # valid_gfa needs the VN:Z:1.0 header tag, and each S line contributes one unitig length.
                if line.startswith("H\t"):
                    if "VN:Z:1.0" in line:
                        valid_gfa = True
                elif line.startswith("S\t"):
                    parts = line.split("\t")
                    if len(parts) > GfaSegmentCol.SEQ:
                        seq = parts[GfaSegmentCol.SEQ]
                        lengths.append(len(seq))
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile):
        return None

    if not lengths:
        return None

    total = sum(lengths)
    unitig_n50 = n50(lengths)

    return {
        "num_unitigs": len(lengths),
        "total_length": total,
        "max_unitig": max(lengths),
        "avg_unitig": round(total / len(lengths), 1),
        "min_unitig": min(lengths),
        "n50": unitig_n50,
        "valid_gfa": valid_gfa,
    }


def _parse_cuttlefish_stats(tmp_dir):
    """Graph statistics of Cuttlefish's output.

    Cuttlefish writes a JSON stats file with structural information, and GFA or FASTA graph files serve as fallback sources.
    """
    # The JSON stats file is the most reliable source and is read first.
    for file_path in Path(tmp_dir).iterdir():
        if file_path.suffix == ".json" and file_path.name.startswith("cdbg"):
            stats = _parse_cuttlefish_json(file_path)
            if stats:
                # The JSON has no N50, which comes from the FASTA when present.
                if "n50" not in stats:
                    for fasta_path in Path(tmp_dir).iterdir():
                        if fasta_path.suffix in (".fa", ".fasta"):
                            fa_stats = _parse_fasta_stats(fasta_path)
                            if fa_stats and "n50" in fa_stats:
                                stats["n50"] = fa_stats["n50"]
                            break
                return stats

    # Next come GFA files, .gfa1 for format 1 and .gfa2 for format 2.
    for ext in [".gfa", ".gfa1", ".gfa2"]:
        for file_path in Path(tmp_dir).iterdir():
            if file_path.suffix == ext or file_path.name.endswith(ext):
                stats = _parse_gfa_stats(file_path)
                if stats:
                    return stats

    # FASTA files are the last resort.
    for file_path in Path(tmp_dir).iterdir():
        if file_path.suffix in (".fa", ".fasta"):
            return _parse_fasta_stats(file_path)
        if file_path.suffix == ".gz" and file_path.name.endswith(".fasta.gz"):
            return _parse_fasta_gz_stats(file_path)

    return None


def _parse_fasta_stats(fasta_path):
    """Sequence statistics of a FASTA file, as the dict of _parse_gfa_stats with valid_gfa False."""
    lengths = []
    try:
        with open(fasta_path, "r") as handle:
            seq_len = 0
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if seq_len > 0:
                        lengths.append(seq_len)
                    seq_len = 0
                elif line:
                    seq_len += len(line)
            if seq_len > 0:
                lengths.append(seq_len)
    except (OSError, UnicodeDecodeError):
        return None

    if not lengths:
        return None

    total = sum(lengths)
    unitig_n50 = n50(lengths)

    return {
        "num_unitigs": len(lengths),
        "total_length": total,
        "max_unitig": max(lengths),
        "min_unitig": min(lengths),
        "avg_unitig": round(total / len(lengths), 1),
        "n50": unitig_n50,
        "valid_gfa": False,
    }


def _parse_fasta_gz_stats(fasta_gz_path):
    """Sequence statistics of a gzipped FASTA file (.fasta.gz)."""
    import gzip

    lengths = []
    try:
        with gzip.open(fasta_gz_path, "rt") as handle:
            seq_len = 0
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if seq_len > 0:
                        lengths.append(seq_len)
                    seq_len = 0
                elif line:
                    seq_len += len(line)
            if seq_len > 0:
                lengths.append(seq_len)
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile):
        return None

    if not lengths:
        return None

    total = sum(lengths)
    unitig_n50 = n50(lengths)

    return {
        "num_unitigs": len(lengths),
        "total_length": total,
        "max_unitig": max(lengths),
        "min_unitig": min(lengths),
        "avg_unitig": round(total / len(lengths), 1),
        "n50": unitig_n50,
        "valid_gfa": False,
    }


def _parse_cuttlefish_json(json_path):
    """Unitig counts, lengths and coverage from Cuttlefish's JSON stats output."""
    try:
        data = json.loads(Path(json_path).read_text())
    except (json.JSONDecodeError, FileNotFoundError):
        return None

    result = {}

    basic = data.get("basic info", {})
    if "vertex count" in basic:
        result["distinct_kmers"] = basic["vertex count"]

    contigs = data.get("contigs info", {})
    if "maximal unitig count" in contigs:
        result["num_unitigs"] = contigs["maximal unitig count"]
    if "sum maximal unitig length" in contigs:
        result["total_length"] = contigs["sum maximal unitig length"]
    if "longest maximal unitig length" in contigs:
        result["max_unitig"] = contigs["longest maximal unitig length"]
    if "avg. maximal unitig length" in contigs:
        result["avg_unitig"] = contigs["avg. maximal unitig length"]
    if "shortest maximal unitig length" in contigs:
        result["min_unitig"] = contigs["shortest maximal unitig length"]

    return result if result else None


class BifrostRunner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for Bifrost.

    Bifrost builds a colored compacted de Bruijn graph.
    Its -s sequence-file mode builds from reads or FASTA, the input mode optik uses.
    The output is a .bfi index and a gzipped GFA 1.0 graph of unitigs, .gfa.gz.

    Because the build accepts ``-t``, ThreadSweepMixin times it at each thread count of the sweep and records RSS per count.
    Measured are build time, the disk footprint of the output files, the unitig count, total sequence length and N50, and peak RSS during construction.
    """

    @classmethod
    def tool_name(cls):
        return "bifrost"

    @classmethod
    def binary_path(cls):
        path = PATHS.get("bifrost")
        if path and check_binary(path):
            return str(path)
        return None

    @property
    def _bifrost(self):
        """Path of the Bifrost binary."""
        return self.binary_path()

    def run(self):
        tmp_dir = self.tmpdir("bifrost")
        try:
            dbg_input, kmer_len, min_count, default_t = dbg_setup(self, tmp_dir)
            # Bifrost documents two filters, -r, which keeps every k-mer for cutoff 1, and -s, which drops k-mers seen once for cutoff 2.
            # It offers no other cutoff.
            input_flag = BIFROST_CUTOFF_FLAGS.get(min_count)
            if input_flag is None:
                return {"error": f"Bifrost has no mode for min_count={min_count}"}

            prefix = f"{tmp_dir}/graph"

            def build_cmd(threads):
                """Bifrost build argv for t threads."""
                return [
                    self._bifrost,
                    "build",
                    input_flag,
                    dbg_input,
                    "-k",
                    str(kmer_len),
                    "-o",
                    prefix,
                    "-t",
                    str(threads if threads else default_t),
                ]

            # The build runs the auto thread sweep.
            # Its gate runs once, untimed, at the first swept count, and a successful gate leaves the graph files on disk for the size and statistics parse below.
            sweep_flat, _, hyperfine_res, _ = self.run_sweep(
                "build",
                build_cmd,
                gate=lambda proc: proc.returncode == 0,
                timeout=BUILD_TIMEOUT,
            )
            if hyperfine_res is None:
                hyperfine_res = self.run_hyperfine(build_cmd(default_t), timeout=BUILD_TIMEOUT)
            if "error" in hyperfine_res:
                return {"error": f"bifrost build hyperfine: {hyperfine_res['error']}"}

            # The index size sums the graph output files, leaving out the decompressed input copy in td.
            idx_bytes = sum(
                file_path.stat().st_size
                for file_path in Path(tmp_dir).iterdir()
                if file_path.is_file() and file_path.name.startswith("graph")
            )

            # Bifrost writes a .bfi binary index and the .gfa.gz graph, whose GFA gives the graph statistics.
            graph_stats = None
            for file_path in Path(tmp_dir).iterdir():
                if file_path.suffix == ".gz" and file_path.name.endswith(".gfa.gz"):
                    graph_stats = _parse_gfa_gz_stats(file_path)
                    break
                elif file_path.suffix in (".gfa", ".gfa2"):
                    graph_stats = _parse_gfa_stats(file_path)
                    break

            result = {
                "build_s": round(hyperfine_res["mean_s"], 4),
                "build_s_stddev": round(hyperfine_res["stddev_s"], 4),
                "build_s_min": round(hyperfine_res["min_s"], 4),
                "build_s_max": round(hyperfine_res["max_s"], 4),
                "build_runs": len(hyperfine_res["times_s"]),
                "build_rss_mb": hyperfine_res["max_rss_mb"],
                "idx_mb": round(idx_bytes / BYTES_PER_MIB, 4) if idx_bytes else None,
            }
            result.update(sweep_flat)

            # Graph statistics, when they parsed
            if graph_stats:
                result["num_unitigs"] = graph_stats["num_unitigs"]
                result["total_length"] = graph_stats["total_length"]
                result["max_unitig"] = graph_stats["max_unitig"]
                result["avg_unitig"] = graph_stats["avg_unitig"]
                result["min_unitig"] = graph_stats.get("min_unitig")
                result["n50"] = graph_stats["n50"]
                result["valid_gfa"] = graph_stats.get("valid_gfa", False)

            result["min_count"] = min_count
            result["mode"] = input_flag
            gfa = [
                file_path
                for file_path in Path(tmp_dir).iterdir()
                if file_path.name.startswith("graph") and ".gfa" in file_path.name
            ]
            (
                result["graph_kmers"],
                result["graph_kmer_occurrences"],
                result["graph_kmers_file"],
            ) = dump_graph_kmers(self, "bifrost", gfa, kmer_len)
            result["input_kmers"] = count_kmers(dbg_input, kmer_len, distinct=False)
            result["distinct_kmers"] = count_kmers(dbg_input, kmer_len, distinct=True)

            return result

        finally:
            self.cleanup(tmp_dir)


class CuttlefishRunner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for Cuttlefish.

    Cuttlefish writes the maximal unitigs of the de Bruijn graph.
    Its --read mode builds from reads or FASTA, the input mode optik uses.
    The output is GFA 1.0 with a JSON stats file.

    As with Bifrost, the build is timed at every thread count of the sweep, recording RSS for each.
    Measured are build time, the graph's disk footprint, its unitig count, total sequence length and N50, and peak RSS during construction.
    """

    @classmethod
    def tool_name(cls):
        return "cuttlefish"

    @classmethod
    def binary_path(cls):
        path = PATHS.get("cuttlefish")
        if path and check_binary(path):
            return str(path)
        return None

    @property
    def _cuttlefish(self):
        """Path of the Cuttlefish binary."""
        return self.binary_path()

    def run(self):
        if not self.available():
            return {"error": "cuttlefish binary not found"}

        tmp_dir = self.tmpdir("cuttlefish")
        try:
            dbg_input, kmer_len, min_count, default_t = dbg_setup(self, tmp_dir)
            # --read builds the edge-centric graph optik builds, where --ref would build the
            # node-centric one.
            # The cutoff is passed explicitly.
            mode = "--read"

            def build_cmd(threads):
                """Cuttlefish build argv, t threads."""
                return [
                    self._cuttlefish,
                    "build",
                    "-s",
                    dbg_input,
                    "-k",
                    str(kmer_len),
                    "-t",
                    str(threads if threads else default_t),
                    "-c",
                    str(min_count),
                    "-o",
                    f"{tmp_dir}/cdbg",
                    mode,
                ]

            # Building sweeps the thread counts too.
            # The untimed gate on the first point also leaves the graph output for the parse below.
            sweep_flat, _, hyperfine_res, _ = self.run_sweep(
                "build",
                build_cmd,
                gate=lambda proc: proc.returncode == 0,
                timeout=BUILD_TIMEOUT,
                ignore_failure=True,
            )
            if hyperfine_res is None:
                hyperfine_res = self.run_hyperfine(
                    build_cmd(default_t), timeout=BUILD_TIMEOUT, ignore_failure=True
                )
            if "error" in hyperfine_res:
                return {"error": f"cuttlefish build hyperfine: {hyperfine_res['error']}"}

            # The index size sums the graph output files, without KMC's temporary bins.
            graph_size = sum(
                file_path.stat().st_size
                for file_path in Path(tmp_dir).iterdir()
                if file_path.is_file() and file_path.name.startswith("cdbg")
            )

            # Graph statistics come from the JSON stats file or the GFA/FASTA output.
            gfa_stats = _parse_cuttlefish_stats(tmp_dir)

            result = {
                "build_s": round(hyperfine_res["mean_s"], 4),
                "build_s_stddev": round(hyperfine_res["stddev_s"], 4),
                "build_s_min": round(hyperfine_res["min_s"], 4),
                "build_s_max": round(hyperfine_res["max_s"], 4),
                "build_runs": len(hyperfine_res["times_s"]),
                "build_rss_mb": hyperfine_res["max_rss_mb"],
                "idx_mb": round(graph_size / BYTES_PER_MIB, 4) if graph_size else None,
            }
            result.update(sweep_flat)

            # Graph statistics, when they parsed
            if gfa_stats:
                for key in (
                    "num_unitigs",
                    "total_length",
                    "max_unitig",
                    "avg_unitig",
                    "n50",
                    "min_unitig",
                    "distinct_kmers",
                ):
                    if key in gfa_stats:
                        result[key] = gfa_stats[key]

            result["min_count"] = min_count
            result["mode"] = f"{mode} -c {min_count}"
            outs = [
                file_path
                for file_path in Path(tmp_dir).iterdir()
                if file_path.name.startswith("cdbg")
                and (file_path.name.endswith(".fa") or ".gfa" in file_path.name)
            ]
            (
                result["graph_kmers"],
                result["graph_kmer_occurrences"],
                result["graph_kmers_file"],
            ) = dump_graph_kmers(self, "cuttlefish", outs, kmer_len)
            result["input_kmers"] = count_kmers(dbg_input, kmer_len, distinct=False)
            # The JSON of --read mode has no vertex count or shortest unitig, so both come from the FASTA unitig file itself, which has no GFA to check.
            if "distinct_kmers" not in result:  # defined as in the other rows
                result["distinct_kmers"] = count_kmers(dbg_input, kmer_len, distinct=True)
            fa_out = [file_path for file_path in outs if file_path.name.endswith(".fa")]
            if "min_unitig" not in result and fa_out:
                fasta_stats = _parse_fasta_stats(str(fa_out[0]))
                if fasta_stats:
                    result["min_unitig"] = fasta_stats["min_unitig"]
            result["valid_gfa"] = "n/a (FASTA output)"

            return result

        finally:
            self.cleanup(tmp_dir)


class OptikDbruijnRunner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for optik's compacted de Bruijn graph (``dbruijn``).

    optik counts the canonical k-mers of the input the reference tools get through minimizer-partitioned super-k-mers, keeps those with abundance of at least -c, and compacts maximal unitigs with a lock-free bidirected walk over one open-addressing node table.
    It writes GFA 1.0 with per-unitig abundance tags and one L line per edge, whereas Cuttlefish writes unitig FASTA only and Bifrost a GFA without abundances.

    The build phase takes ``-t`` and runs the same auto thread-count sweep as Bifrost and Cuttlefish, with an untimed correctness gate at the first swept count whose output feeds the graph statistics below.

    Correctness rests on the GFA header, the graph statistics from the TSV report and the GFA, and the k-mer set checked by the category's graph-parity oracle.
    """

    @classmethod
    def tool_name(cls):
        return "optik"

    @classmethod
    def binary_path(cls):
        path = PATHS.get("optik")
        if path and check_binary(path):
            return str(path)
        return None

    @property
    def _optik(self):
        return self.binary_path()

    def run(self):
        tmp_dir = self.tmpdir("optik_dbruijn")
        try:
            reads, kmer_len, min_count, threads = dbg_setup(self, tmp_dir)
            dbg_input = reads
            verbose = self.ctx.get("verbose", False)

            gfa_path = f"{tmp_dir}/graph.gfa"

            def build_cmd(thread_arg):
                """optik dbruijn command on t threads."""
                return [
                    self._optik,
                    "dbruijn",
                    "-i",
                    reads,
                    "-k",
                    str(kmer_len),
                    "-c",
                    str(min_count),
                    "-t",
                    str(thread_arg if thread_arg else threads),
                    "-o",
                    gfa_path,
                ]

            cmd = build_cmd(threads)

            # One untimed build gives the graph metrics of the TSV report.
            proc, _ = run_cmd(cmd, timeout=BUILD_TIMEOUT, verbose=verbose)
            if proc.returncode != 0:
                return {"error": f"optik dbruijn failed: {proc.stderr[:DETAIL_CHARS]}"}

            # optik prints the TSV report on stdout as key-value pairs such as
            #   input_kmers  4787534
            #   unitigs  1884865
            graph_stats = {}
            for line in proc.stdout.split("\n"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                # Statistics lines are 'key<TAB>value'.
                if len(parts) == 2:
                    try:
                        if parts[1].isdigit():
                            graph_stats[parts[0]] = int(parts[1])
                        else:
                            graph_stats[parts[0]] = float(parts[1])
                    except (ValueError, TypeError):
                        graph_stats[parts[0]] = parts[1]

            # GFA file size
            gfa_bytes = Path(gfa_path).stat().st_size if Path(gfa_path).exists() else 0

            # The GFA adds further graph statistics.
            gfa_stats = _parse_gfa_stats(gfa_path)

            # Thread counts are swept for the build, as for Bifrost and Cuttlefish.
            sweep_flat, _, hyperfine_res, _ = self.run_sweep(
                "build",
                build_cmd,
                gate=lambda proc: proc.returncode == 0,
                timeout=BUILD_TIMEOUT,
            )
            if hyperfine_res is None:
                hyperfine_res = self.run_hyperfine(cmd, timeout=BUILD_TIMEOUT)
            if "error" in hyperfine_res:
                return {"error": f"optik dbruijn hyperfine: {hyperfine_res['error']}"}

            result = {
                "min_count": min_count,
                "build_s": round(hyperfine_res["mean_s"], 4),
                "build_s_stddev": round(hyperfine_res["stddev_s"], 4),
                "build_s_min": round(hyperfine_res["min_s"], 4),
                "build_s_max": round(hyperfine_res["max_s"], 4),
                "build_runs": len(hyperfine_res["times_s"]),
                "build_rss_mb": hyperfine_res["max_rss_mb"],
                "idx_mb": round(gfa_bytes / BYTES_PER_MIB, 4) if gfa_bytes else None,
            }
            result.update(sweep_flat)

            # Graph statistics from the TSV report
            if "input_kmers" in graph_stats:
                result["input_kmers"] = graph_stats["input_kmers"]
            if "unitigs" in graph_stats:
                result["num_unitigs"] = graph_stats["unitigs"]
            if "total_length" in graph_stats:
                result["total_length"] = graph_stats["total_length"]
            if "n50" in graph_stats:
                result["n50"] = graph_stats["n50"]
            if "max_unitig" in graph_stats:
                result["max_unitig"] = graph_stats["max_unitig"]
            if "avg_unitig" in graph_stats:
                result["avg_unitig"] = graph_stats["avg_unitig"]
            if "high_confidence_length" in graph_stats:
                result["high_confidence_length"] = graph_stats["high_confidence_length"]

            # GFA validity, plus any graph statistics the TSV report lacked, derived from the GFA itself.
            if gfa_stats:
                result["valid_gfa"] = gfa_stats["valid_gfa"]
                for key, fallback_key in (
                    ("num_unitigs", "num_unitigs"),
                    ("total_length", "total_length"),
                    ("max_unitig", "max_unitig"),
                    ("avg_unitig", "avg_unitig"),
                    ("n50", "n50"),
                    ("min_unitig", "min_unitig"),
                ):
                    if key not in result and fallback_key in gfa_stats:
                        result[key] = gfa_stats[fallback_key]

            # Distinct k-mers derived from the input, defined the same way for every tool.
            if "distinct_kmers" not in result:
                result["distinct_kmers"] = count_kmers(dbg_input, kmer_len, distinct=True)
            result["mode"] = f"edge-centric -c {min_count}"
            (
                result["graph_kmers"],
                result["graph_kmer_occurrences"],
                result["graph_kmers_file"],
            ) = dump_graph_kmers(self, "optik", [gfa_path], kmer_len)
            result["input_kmers"] = count_kmers(dbg_input, kmer_len, distinct=False)

            return result

        finally:
            self.cleanup(tmp_dir)
