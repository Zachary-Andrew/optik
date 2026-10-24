"""Sequence alignment and read mapping benchmarks.

Tools in this category

    minimap2
        Read-to-reference and all-versus-all aligner with minimizer seeding
        and chaining (Li 2018).  Index construction and alignment
        throughput are measured separately.

    optik
        optik's map subcommand, with minimizer MPHF seeds, locus voting or
        anchor chaining, and its own gap-affine wavefront aligner.  It builds
        the index and aligns in one step.

    wfmash, the wavefront mapper
        MashMap segment mapping followed by WFA base-level alignment, built
        for long sequences with 1 kb segments by default.  tool_bench runs it
        only on reads of 1 kb or more, as real pipelines use it, with the read
        settings -S 1 -p 70 inside its documented options, the
        "wfmash-reads" row of validate_official.WFMASH_READ_OPTS.

Each mapper runs with its own documented defaults for the read type, with no seed parameters forced, on the same simulated reads with a truth file.
All are scored by the rule of bench/validate_official.py, which takes the primary per read and requires the right reference and strand with the read origin, POS minus the leading clip, within max(10, len/10), and counts MAPQ >= 30 errors.
The reference is a genome with its repeats, the dataset's map_reference or its query file when that is a genome.
_map_reference avoids unitigs, since they hold every k-mer once and never test repeat resolution or MAPQ.
The summary reports the correct fraction Corr% and the MAPQ >= 30 errors MQ30x next to Rate%.

Metrics collected
    Index wall time (hyperfine), index size, alignment throughput in reads
    per second, peak RSS of both phases, edit distance statistics, and
    mapping concordance between tools.

What the metrics mean
    Index time is how fast the reference index is built.
    Alignment time is how fast reads are mapped against the index.
    Throughput is reads per second during alignment.
    Memory is the peak RSS during index construction and alignment.
    Accuracy covers correct placements and MAPQ >= 30 errors against the
    simulation truth, the edit distance distribution and the mapping rate.

Correctness checks
    Both tools write SAM, and the comparison covers
    mapping rates, the percentage of reads each tool maps,
    concordance, i.e. agreement on which reads are mapped or unmapped,
    edit distances, the distribution of alignment quality in the NM:i tag, and
    mapping positions for reads that both tools map.

How optik differs from minimap2
    minimap2 seeds with minimizers and chains heuristically.  optik also
    seeds with minimizers but scores each candidate locus exactly with its
    own wavefront aligner (unit cost for short reads, and chaining plus
    gap-affine WFA between anchors for long reads).  minimap2 can separate
    indexing from alignment, whereas optik does both in one step.  Agreement
    in mapping rate is a sign that both implementations are correct.
"""

import sys as _sys
from pathlib import Path

from config import PATHS
from sam_spec import (
    HIGH_MAPQ,
    SAM_FLAG_REVERSE,
    SAM_FLAG_UNMAPPED,
    SAM_FLAGS_NOT_PRIMARY,
    SAM_MIN_COLUMNS,
)

from ._base import (
    ThreadSweepMixin,
    ToolRunner,
    check_binary,
    ensure_uncompressed,
    run_cmd,
)

_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# Timeouts in seconds for the index and alignment phases.
# Alignment of a large query set can be slow, so the limit is generous.
# OPTIK_BENCH_ALIGN_TIMEOUT and OPTIK_BENCH_INDEX_TIMEOUT override them.
import os as _os

from validate_official import (
    WFMASH_READ_OPTS,
)
from validate_official import (
    faidx as _faidx,
)
from validate_official import (
    one_primary_per_read as _one_primary,
)

INDEX_TIMEOUT = int(_os.environ.get("OPTIK_BENCH_INDEX_TIMEOUT", "600"))
ALIGN_TIMEOUT = int(_os.environ.get("OPTIK_BENCH_ALIGN_TIMEOUT", "600"))


# Read synthesis for genome-scale query files
#
# Some bundled datasets ship whole-genome or assembly records where the map category expects reads (se_k31's query, for example, is one wrapped multi-Mb FASTA record).
# Mapping those would make the two tools' metrics incomparable.
# Each tool splits or filters such queries under its own policy, so "mapped_reads" would count alignment fragments rather than reads, and edit distances would reflect segmentation rather than error.
# Genome-scale queries are therefore replaced with reads sampled from the reference itself, the usual setup of mapper benchmarks, so both tools see the same read-shaped input byte for byte.
#
# The sampler is bench/gen_reads.py with a fixed seed.
# Its parameters can be set from the environment.
#   OPTIK_BENCH_MAP_READS   number of reads   [50000]
#   OPTIK_BENCH_MAP_LEN     read length       [150]
#   OPTIK_BENCH_MAP_ERR     per-base error    [0.0]
#   OPTIK_BENCH_MAP_SEED    RNG seed          [42]

# : A record this long is a genome/assembly entry, not a read.
GENOME_RECORD_LEN = 5_000


def _looks_like_fastq(path):
    """True if the file starts with a FASTQ '@' header line."""
    with open(path, "rb") as f:
        return f.readline(1) == b"@"


def _fasta_max_record_len(path):
    """Length of the longest FASTA record, streamed without loading the file."""
    longest = 0
    cur = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                longest = max(longest, cur)
                cur = 0
            else:
                cur += len(line.strip())
    return max(longest, cur)


def _sanitize_fasta_headers(path):
    """Give empty FASTA headers placeholder names, returning the path unchanged when there are none.

    SSHash's stitched unitigs have empty headers (a bare '>').
    minimap2 segfaults writing SAM reference names for them, and gen_reads.py cannot split such a record name, so each empty header becomes '>seq_N'.
    """
    lines = Path(path).read_text().split("\n")
    out_lines = []
    seq_i = 0
    for line in lines:
        if line.startswith(">") and len(line.rstrip()) == 1:
            seq_i += 1
            out_lines.append(f">seq_{seq_i}")
        else:
            out_lines.append(line)
    if seq_i > 0:
        out_path = Path(path).parent / f"_sanitized_{Path(path).name}"
        out_path.write_text("\n".join(out_lines))
        return str(out_path)
    return path


def _ensure_map_reads(queries, unitigs, temp_dir):
    """Return a read-set FASTA for the mapping benchmarks.

    By default error-free reads are synthesized from ``unitigs`` with gen_reads.py and a fixed seed, so both mappers see the same input, all of which is known to map.
    Setting OPTIK_BENCH_MAP_REAL_READS=1 passes the dataset's own reads through untouched instead.
    That is off by default because a dataset's reads need not come from its reference.
    SSHash's SRR5833294 demo reads, for example, map 0 of 1000 against the ecoli1 unitigs, which would zero every mapper in the category.
    """
    import os

    real = os.environ.get("OPTIK_BENCH_MAP_REAL_READS", "").lower() in (
        "1",
        "yes",
        "on",
        "true",
    )
    if real:
        try:
            # FASTQ always holds reads.
            # The FASTA record scanner would misread it, since without '>' lines the whole file looks like one record.
            if _looks_like_fastq(queries):
                return queries
            if _fasta_max_record_len(queries) <= GENOME_RECORD_LEN:
                return queries
        except OSError:
            return queries

    n = int(os.environ.get("OPTIK_BENCH_MAP_READS", "50000"))
    length = int(os.environ.get("OPTIK_BENCH_MAP_LEN", "150"))
    # Illumina-like default of 1% substitutions and 0.1% insertions and deletions
    err = os.environ.get("OPTIK_BENCH_MAP_ERR", "0.01")
    indel = os.environ.get("OPTIK_BENCH_MAP_INDEL", "0.001")
    seed = os.environ.get("OPTIK_BENCH_MAP_SEED", "42")

    import subprocess
    import sys

    gen = Path(__file__).resolve().parents[1] / "gen_reads.py"
    src = _sanitize_fasta_headers(unitigs)
    out = Path(temp_dir) / f"map_reads_n{n}_l{length}_e{err}_i{indel}_s{seed}.fasta"
    if not out.exists():
        cmd = [
            sys.executable,
            str(gen),
            "-r",
            str(src),
            "-n",
            str(n),
            "-l",
            str(length),
            "-e",
            str(err),
            "--ins",
            str(indel),
            "--del",
            str(indel),
            "-s",
            str(seed),
            "-o",
            str(out),
            "--gt",
            str(_truth_path(out)),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0 or not out.exists():
            raise RuntimeError(
                f"read synthesis failed ({proc.returncode}): {proc.stderr[:300]}"
            )
    return str(out)


def _map_reference(ds, temp_dir):
    """Reference FASTA for the mapping category, with sanitized headers.

    This is the dataset's ``map_reference`` when set, else its query file when that is a genome (records longer than GENOME_RECORD_LEN, as in se_k31's ecoli1.fasta), else the stitched unitigs.
    Unitigs hold every k-mer once and have no repeats, so reads simulated from them and mapped back never exercise repeat resolution or MAPQ, which is where tools differ.
    OPTIK_BENCH_MAP_REF=unitigs selects the unitig reference again.
    """
    import os

    if os.environ.get("OPTIK_BENCH_MAP_REF", "").lower() != "unitigs":
        src = ds.get("map_reference")
        if src:
            return _sanitize_fasta_headers(ensure_uncompressed(src, temp_dir))
        q = ensure_uncompressed(ds["queries"], temp_dir)
        try:
            if (
                not _looks_like_fastq(q)
                and _fasta_max_record_len(q) > GENOME_RECORD_LEN
            ):
                return _sanitize_fasta_headers(q)
        except OSError:
            pass
    return _sanitize_fasta_headers(ensure_uncompressed(ds["unitigs"], temp_dir))


def _truth_path(reads):
    """Path of the truth TSV that gen_reads.py writes beside `reads`."""
    return Path(str(reads) + ".truth.tsv")


def _truth_metrics(sam_text, reads):
    """Score a SAM against the simulation truth with the rule of validate_official.score.

    The first primary record of each read is correct when its reference and strand match and the read origin (POS minus the leading S/H clip) lies within max(10, len/10) of the truth.
    Returns {} for reads without truth."""
    import re

    tp = _truth_path(reads)
    if not tp.exists():
        return {}
    truth = {}
    with open(tp) as fh:
        for line in fh:
            t = line.split()
            if t and t[0] != "read_name":
                truth[t[0]] = (t[1], int(t[2]), t[3], int(t[4]))
    prim = {}
    for line in sam_text.split("\n"):
        if not line or line.startswith("@"):
            continue
        t = line.split("\t")
        if len(t) < SAM_MIN_COLUMNS:
            continue
        f = int(t[1])
        if f & SAM_FLAGS_NOT_PRIMARY or t[0] in prim:
            continue
        lead = re.match(r"^(\d+)[SH]", t[5])
        origin = int(t[3]) - 1 - (int(lead.group(1)) if lead else 0)
        prim[t[0]] = (t[2], origin, "-" if f & SAM_FLAG_REVERSE else "+", int(t[4]))
    correct = wrong = hq_wrong = hq = 0
    for name, (ref, pos, strand, span) in truth.items():
        m = prim.get(name)
        if m is None:
            continue
        ok = m[0] == ref and m[2] == strand and abs(m[1] - pos) <= max(10, span // 10)
        correct += ok
        wrong += not ok
        hq += m[3] >= HIGH_MAPQ
        hq_wrong += (not ok) and m[3] >= HIGH_MAPQ
    n = max(1, len(truth))
    return {
        "correct_reads": correct,
        "wrong_reads": wrong,
        "unmapped_reads": len(truth) - correct - wrong,
        "sensitivity": round(correct / n, 5),
        "mapq30_reads": hq,
        "mapq30_wrong": hq_wrong,
    }


def _input_read_count(path):
    """Count the reads in the FASTA or FASTQ query file, plain or gzipped.

    This count, rather than the read names seen in a mapper's SAM, is the denominator of every rate.
    Mappers that omit unmapped records (wfmash) would otherwise get a denominator equal to the reads they mapped."""
    import gzip

    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        first = fh.read(1)
        fh.seek(0)
        if first == "@":
            return sum(1 for _ in fh) // 4
        return sum(1 for line in fh if line.startswith(">"))


def _set_rates(result, sam_stats, queries, seconds):
    """total_reads, mapping_rate and reads_per_second computed from the input reads."""
    total = _input_read_count(queries)
    result["total_reads"] = total
    if total:
        result["mapping_rate"] = round(sam_stats["mapped"] / total * 100, 2)
        if seconds:
            result["reads_per_second"] = round(total / seconds, 2)


def _parse_sam_stats(sam_text):
    """Parse SAM text into mapping statistics.

    Reads are counted by unique QNAME rather than by SAM line, because an aligner may write several supplementary or secondary records for one read.

    The returned dict holds
        mapped             reads with at least one mapped record
        total              distinct query reads seen
        edit_distances     NM:i values of mapped reads
        mapping_positions  {read_name: (ref_name, pos)} of mapped reads
        avg_edit_dist      mean edit distance, if any
        median_edit_dist   median edit distance, if any
    """
    all_names = set()
    mapped_names = set()
    edit_distances = []
    mapping_positions = {}

    for line in sam_text.split("\n"):
        # Header and empty lines carry no record.
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < SAM_MIN_COLUMNS:
            continue

        flag = int(parts[1])
        read_name = parts[0]
        all_names.add(read_name)

        # Unmapped records are skipped.
        if flag & SAM_FLAG_UNMAPPED:
            continue

        mapped_names.add(read_name)
        ref_name = parts[2]
        pos = int(parts[3])

        # The first record of a read sets its position.
        if read_name not in mapping_positions:
            mapping_positions[read_name] = (ref_name, pos)

        # The edit distance comes from the optional NM:i tag.
        for field in parts[11:]:
            if field.startswith("NM:i:"):
                try:
                    edit_distances.append(int(field[5:]))
                except (ValueError, IndexError):
                    pass

                break

    mapped = len(mapped_names)
    total = len(all_names)
    result = {
        "mapped": mapped,
        "total": total,
        "edit_distances": edit_distances,
        "mapping_positions": mapping_positions,
    }

    # Summary statistics of the edit distances
    if edit_distances:
        result["avg_edit_dist"] = round(sum(edit_distances) / len(edit_distances), 2)
        sorted_ed = sorted(edit_distances)
        mid = len(sorted_ed) // 2
        result["median_edit_dist"] = (
            sorted_ed[mid]
            if len(sorted_ed) % 2
            else round(
                (sorted_ed[mid - 1] + sorted_ed[min(mid, len(sorted_ed) - 1)]) / 2, 1
            )
        )

    return result


class Minimap2Runner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for minimap2 mapping.

    The index is built from the reference first and the query reads are then aligned against it.
    Both phases are reported separately and also summed into the total build_s of the summary table.

    Alignment accepts ``-t``, so ThreadSweepMixin times it at every swept thread count and records how RSS grows per thread.
    The index phase, `minimap2 -d`, is single-threaded, so no sweep applies to it.

    Measured are index time to build the .mmi, align time, throughput in mapped reads per second, peak RSS of both phases and the distribution of edit distances.
    Against optik, the interesting figures are index size (the .mmi against optik's in-memory index), align time and throughput, and the edit distance (NM) of the emitted alignments.
    """

    @classmethod
    def tool_name(cls):
        return "minimap2"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("minimap2")
        if p and check_binary(p):
            return str(p)
        return None

    @property
    def _minimap2(self):
        return self.binary_path()

    @staticmethod
    def _sanitize_fasta_headers(path):
        """Module-level _sanitize_fasta_headers, kept as a method for callers that use it."""
        return _sanitize_fasta_headers(path)

    def run(self):
        td = self.tmpdir("minimap2")
        try:
            unitigs = _map_reference(self.ds, td)
            queries = _ensure_map_reads(
                ensure_uncompressed(self.ds["queries"], td), unitigs, td
            )
            verbose = self.ctx.get("verbose", False)
            default_t = int(self.ctx.get("threads_build", 8))

            # Phase 1 builds the minimap2 index (.mmi) with the documented preset for the read type, sr below 1 kb and map-ont above.
            # A prebuilt .mmi fixes k and w, so an index built without -x sr would silently map with k15/w10 instead of the preset's k21/w11.
            preset = "sr" if _fasta_max_record_len(queries) < 1000 else "map-ont"
            idx_path = f"{td}/ref.mmi"
            cmd_idx = [
                self._minimap2,
                "-x",
                preset,
                "-d",
                idx_path,
                unitigs,
            ]

            proc_idx, _ = run_cmd(cmd_idx, timeout=INDEX_TIMEOUT, verbose=verbose)
            if proc_idx.returncode != 0:
                return {"error": f"minimap2 index failed: {proc_idx.stderr[:500]}"}
            idx_bytes = Path(idx_path).stat().st_size if Path(idx_path).exists() else 0

            # Index build timed by hyperfine
            hf_idx = self.run_hyperfine(cmd_idx, timeout=INDEX_TIMEOUT)
            if "error" in hf_idx:
                return {"error": f"minimap2 index hyperfine: {hf_idx['error']}"}

            # Phase 2 aligns the query reads against the index with the auto thread sweep.
            # The gate's untimed first run doubles as the probe whose SAM output feeds the mapping statistics.
            def align_cmd(t):
                """minimap2 command aligning the queries on t threads."""
                return [
                    self._minimap2,
                    "-a",
                    "-t",
                    str(t if t else default_t),
                    "-x",
                    preset,
                    idx_path,
                    queries,
                ]

            probe = {}

            def gate(p):
                """Keep the probe run's SAM and pass when it exited 0."""
                if p.returncode == 0:
                    probe["sam"] = p.stdout
                return p.returncode == 0

            sweep_flat, _, hf_align, _ = self.run_sweep(
                "align", align_cmd, gate=gate, timeout=ALIGN_TIMEOUT
            )
            if hf_align is None:
                hf_align = self.run_hyperfine(
                    align_cmd(default_t), timeout=ALIGN_TIMEOUT
                )
            if "error" in hf_align:
                return {
                    "index_s": round(hf_idx["mean_s"], 4),
                    "idx_mb": round(idx_bytes / (1024 * 1024), 4)
                    if idx_bytes
                    else None,
                    "error_align": hf_align["error"],
                }

            sam_text_probe = probe.get("sam", "")
            sam_stats_probe = _parse_sam_stats(sam_text_probe)

            # Statistics come from the last hyperfine run's SAM, or from the probe when that is empty.
            sam_text_hf = hf_align.get("stdout", "").strip()
            if sam_text_hf:
                sam_stats = _parse_sam_stats(sam_text_hf)
            else:
                sam_stats = sam_stats_probe

            result = {
                "variant": f"-x {preset}",
                "index_s": round(hf_idx["mean_s"], 4),
                "index_s_stddev": round(hf_idx["stddev_s"], 4),
                "align_s": round(hf_align["mean_s"], 4),
                "align_s_stddev": round(hf_align["stddev_s"], 4),
                "build_s": round(hf_idx["mean_s"] + hf_align["mean_s"], 4),
                "build_s_stddev": round(
                    (hf_idx["stddev_s"] ** 2 + hf_align["stddev_s"] ** 2) ** 0.5, 4
                ),
                "build_s_min": round(hf_idx["min_s"] + hf_align["min_s"], 4),
                "build_s_max": round(hf_idx["max_s"] + hf_align["max_s"], 4),
                "build_runs": len(hf_align["times_s"]),
                "idx_mb": round(idx_bytes / (1024 * 1024), 4) if idx_bytes else None,
                "mapped_reads": sam_stats["mapped"],
                "total_reads": sam_stats["total"],
                "build_rss_mb": max(
                    v
                    for v in [hf_idx["max_rss_mb"], hf_align["max_rss_mb"]]
                    if v is not None
                )
                if hf_idx["max_rss_mb"] or hf_align["max_rss_mb"]
                else None,
                "query_rss_mb": hf_align["max_rss_mb"],
            }
            result.update(sweep_flat)
            result.update(_truth_metrics(sam_text_hf or sam_text_probe, queries))

            _set_rates(result, sam_stats, queries, hf_align["mean_s"])

            # Edit distance statistics, when there are any
            if sam_stats["edit_distances"]:
                result["avg_edit_dist"] = sam_stats["avg_edit_dist"]
                result["median_edit_dist"] = sam_stats["median_edit_dist"]
                result["edit_dist_min"] = min(sam_stats["edit_distances"])
                result["edit_dist_max"] = max(sam_stats["edit_distances"])

            return result

        finally:
            self.cleanup(td)


class OptikMapRunner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for optik's read mapping.

    optik maps reads against a reference genome with minimizer seeding and its own wavefront aligner.
    It builds the reference index in memory before mapping, so index construction and alignment form a single step, unlike minimap2's two.

    Mapping gets the same per-thread-count timing (ThreadSweepMixin).

    Measured are total time for index construction plus alignment, throughput in reads mapped per second, peak RSS of the whole mapping run and the NM of the emitted alignments.
    Its mapping rate should be comparable to minimap2's, and reads both tools map should agree.
    """

    @classmethod
    def tool_name(cls):
        return "optik"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("optik")
        if p and check_binary(p):
            return str(p)
        return None

    @property
    def _optik(self):
        return self.binary_path()

    def run(self):
        td = self.tmpdir("optik_map")
        try:
            # The same sanitized reference, and so the same reference names, as the other mappers and the read truth use
            unitigs = _map_reference(self.ds, td)
            queries = _ensure_map_reads(
                ensure_uncompressed(self.ds["queries"], td), unitigs, td
            )
            # optik map builds the index and aligns in one step, with its defaults.
            sam_path = f"{td}/output.sam"
            default_t = int(self.ctx.get("threads_build", 8))

            def map_cmd(t):
                """optik map command on t threads."""
                return [
                    self._optik,
                    "map",
                    "-r",
                    unitigs,
                    "-q",
                    queries,
                    "-o",
                    sam_path,
                    "-t",
                    str(t if t else default_t),
                ]

            # Auto thread sweep, where the gate's untimed first run writes the SAM that the statistics read.
            probe = {}

            def gate(p):
                """Pass when the probe run exited 0."""
                if p.returncode == 0:
                    probe["ok"] = True
                return p.returncode == 0

            sweep_flat, _, hf, _ = self.run_sweep(
                "align", map_cmd, gate=gate, timeout=ALIGN_TIMEOUT
            )
            if hf is None:
                hf = self.run_hyperfine(map_cmd(default_t), timeout=ALIGN_TIMEOUT)
            if "error" in hf:
                return {"error": f"optik map hyperfine: {hf['error']}"}

            # Statistics come from the last hyperfine run's SAM, falling back to the file the gate run wrote.
            sam_text = Path(sam_path).read_text() if Path(sam_path).exists() else ""
            if not sam_text.strip():
                sam_text = hf.get("stdout", "")
            sam_stats = _parse_sam_stats(sam_text)

            result = {
                "variant": "defaults",
                "build_s": round(hf["mean_s"], 4),
                "build_s_stddev": round(hf["stddev_s"], 4),
                "build_s_min": round(hf["min_s"], 4),
                "build_s_max": round(hf["max_s"], 4),
                "build_runs": len(hf["times_s"]),
                "build_rss_mb": hf["max_rss_mb"],
                "query_rss_mb": hf["max_rss_mb"],
                "mapped_reads": sam_stats["mapped"],
                "total_reads": sam_stats["total"],
            }
            result.update(sweep_flat)
            result.update(_truth_metrics(sam_text, queries))

            _set_rates(result, sam_stats, queries, hf["mean_s"])

            # Edit distance statistics over the NM of the emitted alignments
            if sam_stats["edit_distances"]:
                result["avg_edit_dist"] = sam_stats["avg_edit_dist"]
                result["median_edit_dist"] = sam_stats["median_edit_dist"]
                result["edit_dist_min"] = min(sam_stats["edit_distances"])
                result["edit_dist_max"] = max(sam_stats["edit_distances"])

            return result

        finally:
            self.cleanup(td)


class WfmashRunner(ThreadSweepMixin, ToolRunner):
    """Benchmark runner for wfmash, which aligns MashMap segments with WFA.

    Like optik, the whole run of index and map is one command.
    Its -a SAM is reduced to one primary per read, the longest-aligned segment, before scoring, as in the validations.
    With `read_opts` False wfmash runs with its exact defaults."""

    read_opts = False

    @classmethod
    def tool_name(cls):
        return "wfmash-reads" if cls.read_opts else "wfmash"

    @classmethod
    def binary_path(cls):
        p = PATHS.get("wfmash")
        if p and check_binary(p):
            return str(p)
        return None

    def run(self):
        td = self.tmpdir(self.tool_name())
        try:
            unitigs = _map_reference(self.ds, td)
            queries = _ensure_map_reads(
                ensure_uncompressed(self.ds["queries"], td), unitigs, td
            )
            _faidx(unitigs)
            _faidx(queries)
            default_t = int(self.ctx.get("threads_build", 8))
            length = _fasta_max_record_len(queries)
            opts = []
            if self.read_opts:
                opts = WFMASH_READ_OPTS["short" if length < 1000 else "long"]
            raw = f"{td}/wfmash.sam"

            def map_cmd(t):
                """wfmash command on t threads."""
                # hyperfine captures stdout, so the gate run is the one that writes the SAM.
                return [
                    self.binary_path(),
                    unitigs,
                    queries,
                    "-a",
                    "-t",
                    str(t if t else default_t),
                    *opts,
                ]

            probe = {}

            def gate(p):
                """Save the probe run's SAM and pass when it exited 0."""
                if p.returncode == 0:
                    Path(raw).write_text(p.stdout)
                    probe["ok"] = True
                return p.returncode == 0

            sweep_flat, _, hf, _ = self.run_sweep(
                "align", map_cmd, gate=gate, timeout=ALIGN_TIMEOUT
            )
            if hf is None:
                p, _ = run_cmd(map_cmd(default_t), timeout=ALIGN_TIMEOUT)
                gate(p)
                hf = self.run_hyperfine(map_cmd(default_t), timeout=ALIGN_TIMEOUT)
            if "error" in hf:
                return {"error": f"wfmash hyperfine: {hf['error']}"}
            one = f"{td}/wfmash.primary.sam"
            split = _one_primary(raw, one) if Path(raw).exists() else 0
            sam_text = Path(one).read_text() if Path(one).exists() else ""
            sam_stats = _parse_sam_stats(sam_text)
            result = {
                "variant": " ".join(opts) or "defaults",
                "build_s": round(hf["mean_s"], 4),
                "build_s_stddev": round(hf["stddev_s"], 4),
                "build_s_min": round(hf["min_s"], 4),
                "build_s_max": round(hf["max_s"], 4),
                "build_runs": len(hf["times_s"]),
                "build_rss_mb": hf["max_rss_mb"],
                "query_rss_mb": hf["max_rss_mb"],
                "mapped_reads": sam_stats["mapped"],
                "total_reads": sam_stats["total"],
                "reads_split_into_segments": split,
            }
            result.update(sweep_flat)
            result.update(_truth_metrics(sam_text, queries))
            _set_rates(result, sam_stats, queries, hf["mean_s"])
            if sam_stats["edit_distances"]:
                result["avg_edit_dist"] = sam_stats["avg_edit_dist"]
                result["median_edit_dist"] = sam_stats["median_edit_dist"]
            return result
        finally:
            self.cleanup(td)


class WfmashReadsRunner(WfmashRunner):
    """wfmash with the read settings of validate_official (-S 1 -p 70)."""

    read_opts = True
