"""Runs the read-mapping section, comparing optik with the published mappers on simulated and real reads."""

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from diagnostics import DETAIL_CHARS
from external_oracles import SAMTOOLS, external_sam_checks
from fastq import FASTQ_LINES_PER_RECORD, FASTQ_QUAL_LINE, count_fastq_records
from published_catalog import BWA, DATA_ROOT, TOOLS
from published_inputs import _records, reference_input
from published_mapping import (
    MM2_PRESET_KW,
    PAFTOOLS,
    WGSIM_FAIL,
    WGSIM_FAIL_ENDS,
    WGSIM_FAIL_EXIT,
    WGSIM_FAIL_KIND,
    WGSIM_FAIL_MESSAGE,
    WGSIM_FAIL_TOOL,
    WGSIM_ROW,
    WGSIM_ROW_DECILES,
    WGSIM_ROW_ENDS,
    WGSIM_ROW_MAPQ30,
    WGSIM_ROW_MAPQ30_WRONG,
    WGSIM_ROW_RSS,
    WGSIM_ROW_RUNS,
    WGSIM_ROW_SD,
    WGSIM_ROW_TOOL,
    WGSIM_ROW_WALL,
    WGSIM_ROW_WARMUP,
    bwa_index_hg38,
    bwa_index_s,
    count_records,
    index_ready_s,
    mapeval_bins,
    mapeval_calibration,
    paftools_call_snvs,
    primaries,
    seq_bytes,
)
from published_tools import tool_path
from small_common import BIN
from units import BYTES_PER_GIB, BYTES_PER_MB, BYTES_PER_MIB

from published.common import (
    HERE,
    MAPEVAL_MAPQ_CUT,
    MATES_PER_PAIR,
    MM2_ARGV_PREFIX,
    MM2_MAX_K,
    NO_SAM,
    OPTIK_ARGV_PREFIX,
    PAFCMP_SUMMARY_LINES,
    add_status,
    failure,
)

# A shell-wrapped command shows only the script it runs.
SHELL_PREFIX = ["bash", "-c"]


def bench_map(bench, ds_id, dataset, mode, work, threads, mach):
    """Run one round of a map dataset at a given mode, k and thread count.

    The reference copy, the simulated reads and a truth mapping are made once in <work>/<id>_inputs and shared by every mode and round, so all mappers see the same reads.
    Each round writes to a fresh directory, which keeps an older round's output from being scored as this one's."""
    map_k = dataset.get("map_k")
    inp = os.path.join(work, f"{ds_id}_inputs")
    run_dir = os.path.join(work, f"{ds_id}_{mode}_k{map_k or 'default'}_t{threads}")
    shutil.rmtree(run_dir, ignore_errors=True)
    os.makedirs(run_dir)
    os.makedirs(inp, exist_ok=True)
    src = os.path.join(DATA_ROOT.path, ds_id)
    mm2, opt = BIN["minimap2"], BIN["optik"]
    paftools = PAFTOOLS
    protocol = dataset.get("as_", ds_id)  # a sweep dataset runs its source's protocol
    base_row = {
        "section": "map",
        "mode": mode,
        "dataset": ds_id,
        "k": map_k or "default",
        "threads": threads,
    }
    # Each protocol sets whichever of these input paths it uses, and the rest stay None.
    ref = reads_path = reads1 = reads2 = query_fa = None
    if protocol == "bwamem-wgsim-hg38":
        # BWA-MEM Fig. 1, run through the repo's wgsim harness at the paper's parameters with this bench's binaries and memory cap.
        # --fresh makes every mapper map again in this round.
        bwa_prefix = bwa_index_hg38() if os.path.exists(BWA) else None
        base_row["k"] = MM2_PRESET_KW["sr"][
            0
        ]  # optik takes the k and w of minimap2 -ax sr
        rep_f = os.path.join(run_dir, "report.txt")
        with open(rep_f, "w") as rep:
            wgsim_proc = subprocess.run(
                [
                    sys.executable,
                    os.path.join(HERE, "wgsim_check.py"),
                    # Its rows join this bench's map section, so it writes no section of its own.
                    "--report",
                    "",
                    "--ref",
                    dataset["local"],
                    "--out",
                    inp,
                    "--pairs",
                    str(dataset["pairs"]),
                    "--len",
                    str(dataset["length"]),
                    "--error",
                    "0.015",
                    "--mut",
                    "0.002",
                    "--indel-frac",
                    "1.0",
                    "--insert",
                    "500",
                    "--insert-sd",
                    "50",
                    "--window",
                    "20",
                    "--se-both",
                    "--threads",
                    str(threads),
                    "--fresh",
                    "--optik",
                    opt,
                    "--minimap2",
                    mm2,
                    "--bwa",
                    BWA,
                    "--seed-kw",
                    ",".join(map(str, MM2_PRESET_KW["sr"])),
                    "--mem",
                    f"{bench.cap // BYTES_PER_GIB}G",
                    "--warmup",
                    str(bench.warmup),
                    "--runs",
                    str(bench.runs),
                    *(["--bwa-index", bwa_prefix] if bwa_prefix else []),
                ],
                stdout=rep,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        if wgsim_proc.returncode != 0:
            bench.add(
                **base_row,
                tool="wgsim_check",
                gate=False,
                gate_detail=f"wgsim_check exit {wgsim_proc.returncode}: {wgsim_proc.stderr[-DETAIL_CHARS:]}",
            )
            return
        reads1 = next(fname for fname in os.listdir(inp) if fname.endswith(".r1.fq"))
        # Every simulated read is mapped twice, as a pair and again single-end (--se-both).
        # wgsim may write a pair or so fewer than requested.
        with open(os.path.join(inp, reads1)) as handle:
            want = MATES_PER_PAIR * count_fastq_records(handle)
        with open(rep_f) as handle:
            for line in handle:
                fail_match = WGSIM_FAIL.match(line)
                if fail_match:
                    kind = "cap" if fail_match.group(WGSIM_FAIL_KIND) == "CAPPED" else None
                    if kind:
                        add_status(
                            bench,
                            "cap",
                            **base_row,
                            tool=fail_match.group(WGSIM_FAIL_TOOL),
                            reads=fail_match.group(WGSIM_FAIL_ENDS),
                            gate_detail=f"memory cap ({bench.cap // BYTES_PER_GIB} GiB): {fail_match.group(WGSIM_FAIL_MESSAGE)}",
                        )
                    else:
                        bench.add(
                            **base_row,
                            tool=fail_match.group(WGSIM_FAIL_TOOL),
                            reads=fail_match.group(WGSIM_FAIL_ENDS),
                            gate=False,
                            gate_detail=f"exit {fail_match.group(WGSIM_FAIL_EXIT)}: {fail_match.group(WGSIM_FAIL_MESSAGE)}",
                        )
                    continue
                row_match = WGSIM_ROW.match(line)
                if not row_match:
                    continue
                tool, ends = row_match.group(WGSIM_ROW_TOOL), row_match.group(WGSIM_ROW_ENDS)
                sam = [fname for fname in os.listdir(inp) if fname.endswith(f".{tool}.{ends}.sam")]
                got = primaries(os.path.join(inp, sam[0])) if len(sam) == 1 else NO_SAM
                ext_ok, ext = (
                    external_sam_checks(os.path.join(inp, sam[0]), dataset["local"])
                    if len(sam) == 1
                    else (False, "no SAM")
                )
                # Index-ready time and mapping time come from the last measured run's log, as in the other map rows.
                logf = os.path.join(inp, sam[0]) + ".log" if len(sam) == 1 else ""
                ready = (
                    index_ready_s(tool, Path(logf).read_text())
                    if logf and os.path.exists(logf)
                    else None
                )
                split = (
                    {}
                    if ready is None or not row_match.group(WGSIM_ROW_WALL)
                    else {
                        "index_ready_s": round(ready, 2),
                        "mapping_s": round(max(0.0, float(row_match.group(WGSIM_ROW_WALL)) - ready), 2),
                    }
                )
                bench.add(
                    **base_row,
                    tool=tool,
                    reads=ends,
                    wall_s=(
                        {
                            "mean": float(row_match.group(WGSIM_ROW_WALL)),
                            "sd": None
                            if row_match.group(WGSIM_ROW_SD) in (None, "-")
                            else float(row_match.group(WGSIM_ROW_SD)),
                            "runs": int(row_match.group(WGSIM_ROW_RUNS)) if row_match.group(WGSIM_ROW_RUNS) else None,
                            "warmup": int(row_match.group(WGSIM_ROW_WARMUP)) if row_match.group(WGSIM_ROW_WARMUP) else None,
                        }
                        if row_match.group(WGSIM_ROW_WALL)
                        else None
                    ),
                    rss_mb=int(row_match.group(WGSIM_ROW_RSS)) if row_match.group(WGSIM_ROW_RSS) else None,
                    bwa_index_s=bwa_index_s(bwa_prefix)
                    if tool == "bwa" and bwa_prefix
                    else None,
                    mapq30=int(row_match.group(WGSIM_ROW_MAPQ30)),
                    mapq30_wrong=int(row_match.group(WGSIM_ROW_MAPQ30_WRONG)),
                    deciles=row_match.group(WGSIM_ROW_DECILES).strip(),
                    **split,
                    primaries=got,
                    external="pass" if ext_ok else "fail",
                    gate=got == want and ext_ok,
                    params=(
                        "bwa mem on a prebuilt index (bwa_index_s: its build)"
                        if tool == "bwa"
                        else "index built from the FASTA inside the run"
                    )
                    + "; wgsim_check.py, BWA-MEM Fig. 1",
                    gate_detail="; ".join(
                        reason
                        for reason in (
                            ""
                            if got == want
                            else f"{got} primary records for {want} reads",
                            ext,
                        )
                        if reason
                    ),
                )
        return
    tools, want, bwa_prefix, fail = {}, None, None, []
    if protocol == "mm2-pbsim2-chm13":
        ref = reference_input(
            os.path.join(DATA_ROOT.path, protocol, "GCA_009914755.3.fa"), os.path.join(inp, "ref.fa")
        )
        if not os.path.exists(ref + ".fai"):
            subprocess.run([SAMTOOLS or "samtools", "faidx", ref], check=True)
        reads_path = os.path.join(inp, "reads.fa")  # paftools.js pbsim2fq writes FASTA
        if not os.path.exists(reads_path):
            # The paper gives no --seed.
            # A fixed one makes a rerun map the same reads.
            subprocess.run(
                [
                    os.path.join(TOOLS, "pbsim2/src/pbsim"),
                    "--seed",
                    "7",
                    "--depth",
                    "1",
                    "--hmm_model",
                    os.path.join(TOOLS, "pbsim2/data/R94.model"),
                    "--length-min",
                    "5000",
                    "--length-mean",
                    "20000",
                    "--accuracy-mean",
                    "0.95",
                    "--prefix",
                    f"{inp}/sd",
                    ref,
                ],
                check=True,
                cwd=inp,
            )
            mafs = sorted(
                os.path.join(inp, fname)
                for fname in os.listdir(inp)
                if fname.startswith("sd_") and fname.endswith(".maf")
            )
            with open(reads_path + ".part", "w") as handle:
                subprocess.run(
                    paftools + ["pbsim2fq", ref + ".fai", *mafs], stdout=handle, check=True
                )
            os.replace(reads_path + ".part", reads_path)
        if dataset.get("subset"):
            # Keep every n-th read.
            # pbsim2 writes reads chromosome by chromosome, so the first reads would all come from chr1, while a stride keeps each chromosome's share.
            sub = os.path.join(inp, f"reads_stride{dataset['subset']}.fa")
            if not os.path.exists(sub):
                step = max(1, count_records(reads_path) // dataset["subset"])
                with open(sub + ".part", "w") as handle:
                    count = 0
                    for j, (name, seq) in enumerate(_records(reads_path)):
                        if j % step == 0 and count < dataset["subset"]:
                            handle.write(f">{name}\n{seq}\n")
                            count += 1
                os.replace(sub + ".part", sub)
            reads_path = sub
        want = count_records(reads_path)
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                reads_path,
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/optik.sam",
            ],
            "minimap2": [
                mm2,
                "-ax",
                "map-ont",
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/minimap2.sam",
                ref,
                reads_path,
            ],
        }
        evaluate = "mapeval"
    elif protocol == "mm2-mason2":
        ref = reference_input(dataset["local"], os.path.join(inp, "ref.fa"))
        reads1, reads2 = os.path.join(inp, "r1.fq"), os.path.join(inp, "r2.fq")
        if not os.path.exists(reads2):
            truth = os.path.join(inp, "mason.sam")
            subprocess.run(
                [
                    tool_path("mason_simulator"),
                    "-ir",
                    ref,
                    "-n",
                    str(dataset["pairs"]),
                    "--illumina-read-length",
                    "150",
                    "--illumina-prob-mismatch-scale",
                    "2.5",
                    "-o",
                    os.path.join(inp, "m1.fq"),
                    "-or",
                    os.path.join(inp, "m2.fq"),
                    "-oa",
                    truth,
                    "--seed",
                    "7",
                    "--num-threads",
                    str(threads),
                ],
                check=True,
                cwd=inp,
            )

            # Mates come out as consecutive records, NAME/1 then NAME/2, with both truths in NAME.
            # They are streamed one record at a time.
            # The converted FASTQ of 10 M pairs is 7 GB, and holding it here would take the page cache the timed mappers read from.
            def base(rec):
                """Read name without its /1 or /2 mate suffix."""
                return re.sub(r"/[12]$", "", rec[0].split()[0])

            conv = subprocess.Popen(
                paftools + ["mason2fq", truth], stdout=subprocess.PIPE, text=True
            )
            with open(reads1 + ".part", "w") as mate1_handle, open(reads2 + ".part", "w") as mate2_handle:
                prev = None
                while True:
                    rec = [conv.stdout.readline().rstrip("\n") for _ in range(FASTQ_LINES_PER_RECORD)]
                    if not rec[FASTQ_QUAL_LINE] and not rec[0]:
                        break
                    if (
                        prev
                        and base(prev) == base(rec)
                        and prev[0].split()[0].endswith("/1")
                    ):
                        mate1_handle.write("\n".join(prev) + "\n")
                        mate2_handle.write("\n".join(rec) + "\n")
                        prev = None
                    else:
                        prev = rec
            if conv.wait() != 0:
                raise subprocess.CalledProcessError(
                    conv.returncode, paftools + ["mason2fq", truth]
                )
            os.replace(reads1 + ".part", reads1)
            os.replace(reads2 + ".part", reads2)
        with open(reads1) as handle:
            pairs = count_fastq_records(handle)
        if pairs != dataset["pairs"]:
            fail.append(f"{pairs} mason2 pairs paired up of {dataset['pairs']} simulated")
        want = MATES_PER_PAIR * pairs
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                reads1,
                "-q",
                reads2,
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/optik.sam",
            ],
            # -F 800 is minimap2's default and matches the paper's "maximum fragment length to 800 bp".
            "minimap2": [
                mm2,
                "-ax",
                "sr",
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/minimap2.sam",
                ref,
                reads1,
                reads2,
            ],
        }
        if os.path.exists(BWA):
            bwa_prefix = bwa_index_hg38()
            tools["bwa"] = [
                "bash",
                "-c",
                f"{BWA} mem -t {threads} {bwa_prefix} {reads1} {reads2} > {run_dir}/bwa.sam",
            ]
        evaluate = "mapeval"
    elif protocol == "mm2-winnowmap-truth":
        ref = reference_input(
            os.path.join(src, "GCA_009914755.3.fa"), os.path.join(inp, "ref.fa")
        )
        reads_path = min(
            os.path.join(src, fname) for fname in os.listdir(src) if fname.endswith(".fastq.gz")
        )
        base = os.path.join(inp, "winnowmap.paf")
        if not os.path.exists(base):
            # Winnowmap's README recipe counts repetitive 15-mers with meryl, then maps HiFi with -x map-pb.
            subprocess.run(
                [tool_path("meryl"), "count", "k=15", "output", f"{inp}/merylDB", ref],
                check=True,
            )
            with open(f"{inp}/repetitive_k15.txt", "w") as handle:
                subprocess.run(
                    [
                        tool_path("meryl"),
                        "print",
                        "greater-than",
                        "distinct=0.9998",
                        f"{inp}/merylDB",
                    ],
                    stdout=handle,
                    check=True,
                )
            with open(base + ".part", "w") as handle:
                subprocess.run(
                    [
                        tool_path("winnowmap"),
                        "-W",
                        f"{inp}/repetitive_k15.txt",
                        "-x",
                        "map-pb",
                        "-t",
                        str(threads),
                        ref,
                        reads_path,
                    ],
                    stdout=handle,
                    check=True,
                )
            os.replace(base + ".part", base)
        want = count_records(reads_path)
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref,
                "-q",
                reads_path,
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/optik.sam",
            ],
            # SAM with base-level alignments, as optik writes.
            # PAF without -c would skip minimap2's alignment step and time a different, cheaper product.
            "minimap2": [
                mm2,
                "-ax",
                "map-hifi",
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/minimap2.sam",
                ref,
                reads_path,
            ],
        }
        evaluate = "pafcmp"
    elif protocol == "bwamem-ecoli-asm":
        query_fa = reference_input(
            os.path.join(src, "NC_000913.3.fa"), os.path.join(inp, "query.fa")
        )
        ref_path = ref = reference_input(
            os.path.join(src, "NC_008253.1.fa"), os.path.join(inp, "ref.fa")
        )
        want = 1  # a single query sequence yields one primary record
        tools = {
            "optik": [
                opt,
                "map",
                "-r",
                ref_path,
                "-q",
                query_fa,
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/optik.sam",
            ],
            "minimap2": [
                mm2,
                "-ax",
                "asm20",
                "--eqx",
                "-t",
                str(threads),
                "-o",
                f"{run_dir}/minimap2.sam",
                ref_path,
                query_fa,
            ],
        }
        if os.path.exists(BWA):
            # The paper's own run uses BWA-MEM with defaults on the whole MG1655 sequence.
            if not os.path.exists(ref_path + ".bwt"):
                start = time.time()
                subprocess.run([BWA, "index", ref_path], check=True, capture_output=True)
                Path(ref_path + ".index_s").write_text(f"{time.time() - start:.2f}")
            tools["bwa"] = [
                "bash",
                "-c",
                f"{BWA} mem -t {threads} {ref_path} {query_fa} > {run_dir}/bwa.sam",
            ]
            bwa_prefix = ref_path
        evaluate = "call"
    warm_set = [path for path in (ref, reads_path, reads1, reads2, query_fa) if isinstance(path, str)] + (
        [bwa_prefix + ext for ext in (".bwt", ".sa", ".pac", ".ann", ".amb")]
        if bwa_prefix
        else []
    )
    if dataset.get("no_bwa") and "bwa" in tools:
        # bwa mem has no comparable k, so in k-scaling rounds its row at this thread count comes from the thread-scaling round.
        del tools["bwa"]
    if map_k:
        # Matched seeds give optik and minimap2 the same k and w = 10, with their presets setting the rest.
        # bwa mem keeps its default, since its -k is a minimum exact-match length.
        tools["optik"] = (
            tools["optik"][:OPTIK_ARGV_PREFIX] + ["-k", str(map_k), "-w", "10"] + tools["optik"][OPTIK_ARGV_PREFIX:]
        )
        if map_k > MM2_MAX_K:
            del tools["minimap2"]
            add_status(
                bench,
                "skip",
                **base_row,
                tool="minimap2",
                gate_detail="minimap2 supports k <= 28",
            )
        else:
            tools["minimap2"] = (
                tools["minimap2"][:MM2_ARGV_PREFIX]
                + ["-k", str(map_k), "-w", "10"]
                + tools["minimap2"][MM2_ARGV_PREFIX:]
            )
    elif "minimap2" in tools and "optik" in tools:
        # Without a sweep, optik takes the k and w of the minimap2 preset the paper runs, so both index the same minimizers.
        # optik's own choice of k is not compared here.
        preset_k, preset_w = MM2_PRESET_KW[tools["minimap2"][tools["minimap2"].index("-ax") + 1]]
        tools["optik"] = (
            tools["optik"][:OPTIK_ARGV_PREFIX] + ["-k", str(preset_k), "-w", str(preset_w)] + tools["optik"][OPTIK_ARGV_PREFIX:]
        )
        base_row["k"] = preset_k
    for tool, cmd in tools.items():
        tag = f"{mode}/{ds_id}/k{map_k or 'default'}/t{threads}/{tool}"
        # Every mapper run starts with identical cached inputs (measure.py).
        # Every run writes its SAM, and the gate below reads the one from the last measured run.
        result = bench.measure(tag, cmd, cwd=run_dir, inputs=warm_set)
        rc_, rss, err = (0 if result.succeeded else result.code), result.rss() or bench.cap / BYTES_PER_MIB, result.err
        if rc_ != 0:
            kind, det = failure(bench, rc_, rss, err)
            if kind == "cap":
                add_status(bench, "cap", **base_row, tool=tool, gate_detail=det)
            else:
                bench.add(**base_row, tool=tool, gate=False, gate_detail=det)
            continue
        # optik and minimap2 build their index from the FASTA inside the run, whereas bwa mem reads a prebuilt one whose build time is reported beside it (bwa_index_s).
        shown = cmd[-1] if cmd[: len(SHELL_PREFIX)] == SHELL_PREFIX else " ".join(map(str, cmd))
        row = dict(
            **base_row,
            tool=tool,
            wall_s=result.wall_row(),
            rss_mb=round(rss, 1),
            params=shown
            + (
                "  [prebuilt bwa index; its build in bwa_index_s]"
                if tool == "bwa"
                else "  [index built from the FASTA inside this run]"
            ),
        )
        if tool == "bwa":
            row["bwa_index_s"] = bwa_index_s(bwa_prefix)
        # optik and minimap2 index the reference inside every run and log when the index was ready, counted from process start.
        # The remainder of the run is mapping alone, so the three mappers can be compared with and without the index build.
        ready = index_ready_s(tool, err)
        if ready is not None and result.succeeded:
            row["index_ready_s"] = round(ready, 2)
            row["mapping_s"] = round(max(0.0, result.wall["mean"] - ready), 2)
        sam_out = f"{run_dir}/{tool}.sam"
        if not os.path.exists(sam_out) or os.path.getsize(sam_out) == 0:
            bench.add(**row, gate=False, gate_detail="exit 0 but no SAM written")
            continue
        ext_ok, ext = external_sam_checks(sam_out, ref)
        n_prim = primaries(sam_out)
        row.update(
            external="pass" if ext_ok else "fail",
            primaries=n_prim,
            seq_mb=round(seq_bytes(sam_out) / BYTES_PER_MB, 2),
        )
        det = [
            reason
            for reason in fail
            + [
                ext,
                "" if n_prim == want else f"{n_prim} primary records for {want} reads",
            ]
            if reason
        ]
        passed = ext_ok and n_prim == want and not fail
        try:
            if evaluate == "mapeval":
                mapeval_out = subprocess.run(
                    paftools + ["mapeval", sam_out],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                Path(f"{run_dir}/{tool}.mapeval").write_text(mapeval_out)
                # Bins come in descending MAPQ, so the last one at or above the cut holds the cumulative figures the papers plot.
                q10 = [mapq_bin for mapq_bin in mapeval_bins(mapeval_out) if mapq_bin.mapq >= MAPEVAL_MAPQ_CUT]
                r10 = q10[-1] if q10 else None
                if r10 is None:
                    det.append("mapeval produced no MAPQ >= 10 row")
                bench.add(
                    **row,
                    mapq10_mapped=r10.cumulative_mapped if r10 else None,
                    mapq10_err_rate=r10.cumulative_error_rate if r10 else None,
                    calibrated=mapeval_calibration(mapeval_out),
                    gate=passed and r10 is not None,
                    gate_detail="; ".join(det),
                )
            elif evaluate == "pafcmp":
                # Both mappers' SAM go through the same conversion.
                test = f"{run_dir}/{tool}.paf"
                with open(test, "w") as handle:
                    subprocess.run(
                        paftools + ["sam2paf", sam_out], stdout=handle, check=True
                    )
                cmp_ = subprocess.run(
                    paftools + ["pafcmp", base, test],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                Path(f"{run_dir}/{tool}.pafcmp").write_text(cmp_)
                if not cmp_.strip():
                    det.append("pafcmp produced nothing")
                bench.add(
                    **row,
                    pafcmp=" | ".join(
                        line.strip() for line in cmp_.strip().splitlines()[-PAFCMP_SUMMARY_LINES:]
                    ),
                    gate=passed and bool(cmp_.strip()),
                    gate_detail="; ".join(det),
                )
            else:
                # paftools.js call, minimap2's assembly-variant caller, counts substitutions over the bases one alignment covers.
                # Picard's check above covers CIGAR/SEQ consistency.
                call_snv, cov1 = paftools_call_snvs(sam_out, ref, run_dir, tool)
                if not cov1:
                    det.append("paftools.js call: no bases covered by one alignment")
                bench.add(
                    **row,
                    substitutions=call_snv,
                    aligned_bases_cov1=cov1,
                    gate=passed and cov1 > 0,
                    gate_detail="; ".join(det),
                )
        except subprocess.CalledProcessError as call_err:
            bench.add(
                **row,
                gate=False,
                gate_detail="; ".join(
                    det
                    + [
                        f"{evaluate} failed (exit {call_err.returncode}): {(call_err.stderr or '')[-DETAIL_CHARS:]}"
                    ]
                ),
            )
