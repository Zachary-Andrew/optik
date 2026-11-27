"""Runs the de Bruijn graph section.

Inputs are prepared first, and then every graph builder runs its rounds while the disk peak is tracked.
"""

import os
import re
import threading

import compare_tools_small_sets
from diagnostics import DETAIL_CHARS
from external_oracles import judge_unitigs, kmc_count, kmc_db, same_unitigs
from published_catalog import DATA_ROOT
from published_inputs import _records, reference_input
from small_common import BIN
from units import BYTES_PER_GB, BYTES_PER_MIB

from published.common import (
    KMC_MAX_MEM_GB,
    KMC_MEM_SHARE,
    TOOL_BUDGET_SHARE,
    add_status,
    failure,
)
from published.machine import tool_skip

# A read set's slice takes the two mate files, each with half of the base budget.
MATE_FILES = 2
# Read sets also run at cutoff 2, the only cutoff Bifrost takes.
READS_CUTOFF = 2


def dbg_inputs(ds_id, dataset, kmer_len, work, limit_bases=None):
    """Return the inputs of the graph tools.

    A reference is formatted with one record per line, upper case and first-token headers.
    Read files are passed as fetched, since every tool reads FASTQ.gz, and `smoke` takes a slice of them."""
    src = os.path.join(DATA_ROOT.path, dataset.get("as_", ds_id))
    base = os.path.join(work, f"{ds_id}_k{kmer_len}" + ("_smoke" if limit_bases else ""))
    os.makedirs(base, exist_ok=True)
    if dataset["kind"] == "reads":
        fastq_files = sorted(
            os.path.join(src, path) for path in os.listdir(src) if path.endswith(".fastq.gz")
        )
        if dataset.get("runs_frac"):
            # A fixed subsample of the first runs, both mates, that hold runs_frac of the compressed bytes
            runs = sorted(
                {re.sub(r"_[12]\.fastq\.gz$", "", os.path.basename(path)) for path in fastq_files}
            )
            total = sum(os.path.getsize(path) for path in fastq_files)
            keep, acc, want = set(), 0, dataset["runs_frac"] * total
            for run_name in runs:
                run_bytes = sum(
                    os.path.getsize(path)
                    for path in fastq_files
                    if os.path.basename(path).startswith(run_name + "_")
                )
                # The run that crosses the target is kept only if keeping it lands nearer the target.
                if acc + run_bytes > want and acc + run_bytes - want >= want - acc:
                    break
                keep.add(run_name)
                acc += run_bytes
                if acc >= want:
                    break
            if not keep and runs:
                # A prefix-only set, as in smoke, can be smaller than one run's share, so its first run stands in.
                first = runs[0]
                keep.add(first)
                acc = sum(
                    os.path.getsize(path)
                    for path in fastq_files
                    if os.path.basename(path).startswith(first + "_")
                )
            fastq_files = [
                path
                for path in fastq_files
                if re.sub(r"_[12]\.fastq\.gz$", "", os.path.basename(path)) in keep
            ]
            print(
                f"  {ds_id}: the first {len(keep)} runs of {len(runs)}, {acc / total:.1%} of the compressed bytes "
                f"(target {dataset['runs_frac']:.1%}; whole runs, both mates)",
                flush=True,
            )
        if not limit_bases:
            return fastq_files, True
        out = []
        for j, fastq in enumerate(fastq_files[:MATE_FILES]):
            slice_path = os.path.join(base, f"slice{j}.fq")
            if not os.path.exists(slice_path):
                with open(slice_path, "w") as handle:
                    handle.writelines(
                        f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"
                        for name, seq in _records(fastq, limit_bases // MATE_FILES)
                    )
            out.append(slice_path)
        return out, True
    srcf = dataset["local"] if dataset.get("local") else os.path.join(src, dataset["input"])
    return [reference_input(srcf, os.path.join(base, "ref.fa"), limit_bases)], False


class DiskPeak:
    """Peak disk space a run takes on the filesystem of `path`.

    The filesystem's used bytes, as statvfs and df report them, are sampled every second and the starting value is subtracted.
    Counting files in the work directory would miss optik's spill files, which O_TMPFILE never links, and Cuttlefish's deleted buckets.
    Because the bench runs one tool at a time, the filesystem's growth belongs to that tool, output file included.
    Cuttlefish 2 reports its intermediate disk use the same way in its Tables S1 and S3."""

    def __init__(self, path):
        self.path, self.peak, self._stop = path, 0, threading.Event()

    def _used(self):
        """Bytes in use on the filesystem."""
        fs_stat = os.statvfs(self.path)
        return (fs_stat.f_blocks - fs_stat.f_bfree) * fs_stat.f_frsize

    def __enter__(self):
        self.base = self._used()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def _run(self):
        """Sampler thread body, polling once a second until __exit__."""
        while not self._stop.wait(1.0):
            self.peak = max(self.peak, self._used() - self.base)

    def __exit__(self, *exc_info):
        self._stop.set()
        self._thread.join()
        self.peak = max(self.peak, self._used() - self.base)


def bench_dbg(
    bench, ds_id, dataset, mode, work, threads, mach, limit_bases=None, only_tools=None
):
    """Run every graph tool at the same k, cutoff and thread count.

    On reads the paper's f0 compares only two tools, because optik and Cuttlefish take a cutoff while Bifrost's is fixed at 2.
    A matched round at cutoff 2 therefore adds Bifrost."""
    kmc_gb = max(1, min(KMC_MAX_MEM_GB, int(mach["mem_gib"] * KMC_MEM_SHARE)))
    for kmer_len, min_count in [
        (kmer_len, min_count)
        for kmer_len in dataset["k"]
        for min_count in dict.fromkeys(
            [dataset["cutoff"]]
            + (
                [READS_CUTOFF]
                if dataset.get("kind") == "reads"
                and not limit_bases
                and not dataset.get("no_c2_round")
                else []
            )
        )
    ]:
        inputs, reads_input = dbg_inputs(ds_id, dataset, kmer_len, work, limit_bases)
        # Outputs go to a directory per thread count, while the input-only KMC databases are shared.
        shared = os.path.join(
            work, f"{ds_id}_k{kmer_len}" + ("_smoke" if limit_bases else ""), f"kmc_c{min_count}"
        )
        run_dir = os.path.join(
            work,
            f"{ds_id}_k{kmer_len}" + ("_smoke" if limit_bases else ""),
            f"{mode}_c{min_count}_t{threads}",
        )
        os.makedirs(run_dir, exist_ok=True)
        os.makedirs(shared, exist_ok=True)
        truth = kmc_db(
            inputs, kmer_len, min_count, os.path.join(shared, "truth"), threads, not reads_input, kmc_gb
        )
        # A reference's distinct k-mers feed the memory gate.
        # Read sets use published figures only.
        n_truth = None if reads_input else kmc_count(truth)
        tool_threads = threads
        # Tools that take a memory budget share one, half of this machine's (or --cap-gb's) memory.
        # It sets Cuttlefish's -m and, in parity mode, optik's --max-mem.
        # In optimal mode optik keeps its documented default of half of MemAvailable at launch, cgroup limits included.
        # A larger --max-mem only keeps more counting buffers resident at the same speed.
        budget = max(1, int(mach["mem_gib"] * TOOL_BUDGET_SHARE))
        tools = {
            "optik": (
                [
                    BIN["optik"],
                    "dbruijn",
                    *[arg for input_path in inputs for arg in ("-i", input_path)],
                    "-k",
                    str(kmer_len),
                    "-c",
                    str(min_count),
                    "-t",
                    str(tool_threads),
                    "-o",
                    f"{run_dir}/optik.fa",
                ]
                + (["--max-mem", f"{budget}G"] if mode == "parity" else []),
                f"{run_dir}/optik.fa",
                "node",
            ),
            "cuttlefish": (
                [
                    BIN["cuttlefish"],
                    "build",
                    *[arg for input_path in inputs for arg in ("-s", input_path)],
                    "-k",
                    str(kmer_len),
                    "-t",
                    str(tool_threads),
                    "-m",
                    str(budget),
                    "-o",
                    f"{run_dir}/cf",
                    "-w",
                    run_dir,
                    "--read" if reads_input else "--ref",
                ]
                + (["-c", str(min_count)] if reads_input else []),
                f"{run_dir}/cf.fa",
                "edge",
            ),
            "bifrost": (
                [
                    BIN["bifrost"],
                    "build",
                    *[arg for input_path in inputs for arg in (["-s"] if reads_input else ["-r"]) + [input_path]],
                    "-k",
                    str(kmer_len),
                    "-t",
                    str(tool_threads),
                    "-o",
                    f"{run_dir}/bf",
                    "-f",
                    "-n",
                    "-N",
                ],
                f"{run_dir}/bf.fasta",
                "node",
            ),
        }
        stats = {}
        row0 = len(
            bench.rows
        )  # the unitig cross-check below marks only this round's rows
        for tool, (cmd, outp, model) in tools.items():
            if (only_tools and tool not in only_tools) or (
                dataset.get("tools") and tool not in dataset["tools"]
            ):
                continue
            tag = f"{mode}/{ds_id}/k{kmer_len}/c{min_count}/t{threads}/{tool}"
            params = " ".join(map(str, cmd[1:]))
            if tool == "bifrost" and reads_input and min_count != READS_CUTOFF:
                continue  # Bifrost's cutoff is fixed at 2, so the c = 2 round compares it
            why = None if limit_bases else tool_skip(dataset, tool, mach, kmer_len, n_truth)
            if why:
                add_status(
                    bench,
                    "skip",
                    section="dbg",
                    mode=mode,
                    dataset=ds_id,
                    k=kmer_len,
                    cutoff=min_count,
                    threads=threads,
                    tool=tool,
                    gate_detail=why,
                )
                continue
            # Every run, warmups included, starts with the inputs in the page cache (measure.py), since the KMC judging evicts them.
            # The disk peak spans all runs, with each run's output replacing the previous one.
            with DiskPeak(run_dir) as disk_peak:
                result = bench.measure(tag, cmd, cwd=run_dir, inputs=inputs)
            rc_, rss, err = (0 if result.succeeded else result.code), result.rss() or bench.cap / BYTES_PER_MIB, result.err
            wrote, readb = result.written, result.read  # median block I/O of a measured run
            if rc_ != 0 and failure(bench, rc_, rss, err)[0] == "cap":
                add_status(
                    bench,
                    "cap",
                    section="dbg",
                    mode=mode,
                    dataset=ds_id,
                    k=kmer_len,
                    cutoff=min_count,
                    threads=threads,
                    tool=tool,
                    gate_detail=failure(bench, rc_, rss, err)[1],
                    params=params,
                )
                continue
            # Some tools compress their output or pick their own extension, so those paths are tried before the run counts as failed.
            if rc_ != 0 or not os.path.exists(outp):
                alt = [
                    path
                    for path in (
                        outp + ".gz",
                        outp.replace(".fasta", ".fa"),
                        outp.replace(".fasta", ".gfa"),
                    )
                    if os.path.exists(path)
                ]
                if rc_ != 0 or not alt:
                    bench.add(
                        section="dbg",
                        mode=mode,
                        dataset=ds_id,
                        k=kmer_len,
                        cutoff=min_count,
                        threads=threads,
                        tool=tool,
                        gate=False,
                        gate_detail=f"exit {rc_}, no output: {err[-DETAIL_CHARS:]}",
                        params=params,
                    )
                    continue
                outp = alt[0]
            # Every tool's output, GFA or FASTA, is rewritten as one plain FASTA, so the judge reads a single format.
            unitigs_fa = os.path.join(run_dir, f"{tool}.unitigs.fa")
            num_unitigs = total_bp = 0
            with open(unitigs_fa, "w") as handle:
                for seq in compare_tools_small_sets.gfa_or_fasta_seqs(outp):
                    handle.write(f">{num_unitigs}\n{seq}\n")
                    num_unitigs += 1
                    total_bp += len(seq)
            passed, det, n_tool, n_model = judge_unitigs(
                tool, unitigs_fa, inputs, kmer_len, min_count, reads_input, run_dir, threads, kmc_gb, truth, shared
            )
            stats[tool] = unitigs_fa
            bench.add(
                section="dbg",
                mode=mode,
                dataset=ds_id,
                k=kmer_len,
                tool=tool,
                cutoff=min_count,
                threads=threads,
                wall_s=result.wall_row(),
                rss_mb=round(rss, 1),
                written_gb=None if wrote is None else round(wrote / BYTES_PER_GB, 2),
                read_gb=None if readb is None else round(readb / BYTES_PER_GB, 2),
                peak_disk_gb=round(disk_peak.peak / BYTES_PER_GB, 2),
                unitigs=num_unitigs,
                total_bp=total_bp,
                kmers=n_tool,
                truth_kmers=n_model,
                graph=f"{model}-centric",
                gate=passed,
                gate_detail=det,
                params=params,
            )
        # Maximal unitigs are unique, so on a genome optik's must match Cuttlefish's string for string.
        if (
            "optik" in stats
            and "cuttlefish" in stats
            and not reads_input
            and not same_unitigs(stats["optik"], stats["cuttlefish"])
        ):
            for row in bench.rows[row0:]:
                if (
                    row.get("dataset"),
                    row.get("k"),
                    row.get("cutoff"),
                    row.get("threads"),
                    row.get("mode"),
                    row.get("tool"),
                ) == (ds_id, kmer_len, min_count, threads, mode, "optik"):
                    row["gate"] = False
                    row["gate_detail"] += (
                        "; unitigs differ from Cuttlefish's (maximal unitigs are unique)"
                    )
