"""Mapping helpers of compare_tools_published_sets.py.

The module wraps the papers' own evaluators (paftools.js mapeval and call, wgsim_eval.pl), the prebuilt bwa index of GRCh38 and the small SAM summaries that the map section reports.
"""

import os
import re
import subprocess
import time
from pathlib import Path
from typing import NamedTuple

from external_oracles import SAMTOOLS
from mapq_calibration import overconfident
from published_catalog import (
    BWA,
    BWA_IDX,
    HG38,
    TEMP_DIR,
    TOOLS,
)
from published_inputs import _records
from sam_spec import (
    PHRED_SCALE,
    PafCol,
    SamCol,
)

PAFTOOLS = [os.path.join(TOOLS, "k8"), os.path.join(TEMP_DIR, "minimap2/misc/paftools.js")]

# over-confident bins named in a failure message
OVERCONFIDENT_BINS_SHOWN = 3
# primary_count() result when samtools fails
UNKNOWN_COUNT = -1


# Columns of a paftools.js call V line that hold the reference and alternative alleles.
V_REF_ALLELE_COL = 6
V_ALT_ALLELE_COL = 7


class MapevalBin(NamedTuple):
    """One "Q" row of paftools.js mapeval output, for the reads mapped at one MAPQ."""

    mapq: int
    reads: int
    wrong: int
    cumulative_error_rate: float
    """Wrong share of all reads mapped at this MAPQ or above."""
    cumulative_mapped: int
    """Reads mapped at this MAPQ or above."""


def mapeval_bins(mapeval_out):
    """The MAPQ bins of mapeval output in its own order, which is descending MAPQ.

    Its rows read "Q mapq reads wrong cum_err_rate cum_mapped".
    """
    bins = []
    for row in mapeval_out.splitlines():
        kind, *values = row.split("\t")
        if kind == "Q":
            mapq, reads, wrong, cumulative_error_rate, cumulative_mapped, *_ = values
            bins.append(MapevalBin(int(mapq), int(reads), int(wrong), float(cumulative_error_rate), int(cumulative_mapped)))
    return bins


def mapeval_calibration(mapeval_out):
    """Over-confident MAPQ bins in paftools.js mapeval output, whose rows read 'Q mapq n wrong ...'.

    implied = n x 10^(-Q/10), and mapq_calibration.overconfident() decides, so every bench applies the same rule."""
    over = []
    for mapq_bin in mapeval_bins(mapeval_out):
        implied = mapq_bin.reads * PHRED_SCALE ** (-mapq_bin.mapq / PHRED_SCALE)
        if overconfident(mapq_bin.wrong, implied):
            over.append(f"Q{mapq_bin.mapq}: {mapq_bin.wrong} wrong vs {implied:.1f} implied")
    return (
        "yes"
        if not over
        else f"no: {len(over)} bins over-confident ({'; '.join(over[:OVERCONFIDENT_BINS_SHOWN])})"
    )


def paftools_call_snvs(sam, ref, work_dir, tag):
    """Count substitutions called by minimap2's assembly-variant caller in the regions one alignment covers.

    paftools.js call runs on the cs-tagged PAF of the primary and supplementary records.
    Returns (SNVs, bases covered once)."""
    md_sam = os.path.join(work_dir, f"{tag}.md.sam")
    with open(md_sam, "w") as out_file:
        subprocess.run(
            [SAMTOOLS or "samtools", "calmd", sam, ref],
            stdout=out_file,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    paf = subprocess.run(
        PAFTOOLS + ["sam2paf", "-p", md_sam], capture_output=True, text=True, check=True
    ).stdout
    rows = sorted(
        (row for row in paf.splitlines() if row),
        key=lambda row: (row.split("\t")[PafCol.TNAME], int(row.split("\t")[PafCol.TSTART])),
    )
    out = subprocess.run(
        PAFTOOLS + ["call", "-l", "1000", "-L", "1000", "-"],
        input="\n".join(rows) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    snv = cov = 0
    # An "R" row reads "R contig start end" for a covered region [start, end).
    # A "V" row reads "V contig start end depth mapq ref alt" for a variant.
    for row in out.splitlines():
        fields = row.split("\t")
        if fields[0] == "R":
            _, _contig, start, end, *_ = fields
            cov += int(end) - int(start)
        elif fields[0] == "V":
            ref_allele, alt_allele = fields[V_REF_ALLELE_COL], fields[V_ALT_ALLELE_COL]
            if len(ref_allele) == 1 and len(alt_allele) == 1:
                snv += 1
    os.remove(md_sam)
    return snv, cov


def index_ready_s(tool, err):
    """Seconds from process start until the mapper's index was ready, read from its own log.

    The log is the last measured run's stderr.
    optik writes "[map] ... minimizers indexed (...) in S s", which covers reading the reference, the first read batch and the index build.
    minimap2 writes "[M::main::S*...] loaded/built the index".
    A mapper that reads a prebuilt index, as bwa does per bwa_index_s, gives None."""
    if tool == "optik":
        times = re.findall(r"minimizers indexed .*? in ([\d.]+)s", err)
    elif tool == "minimap2":
        times = re.findall(r"\[M::main::([\d.]+)\*[\d.]+\] loaded/built the index", err)
    else:
        return None
    return float(times[-1]) if times else None


def bwa_index_hg38():
    """BWA-MEM's index of HG38, or None when bwa is absent.

    BWA-MEM is the papers' short-read mapper.
    The index is built once, which takes about an hour, and the build's wall time is stored beside it (<prefix>.index_s) to be reported with bwa's rows."""
    if not os.path.exists(BWA):
        return None
    if not os.path.exists(BWA_IDX + ".bwt"):
        os.makedirs(os.path.dirname(BWA_IDX), exist_ok=True)
        if not os.path.lexists(BWA_IDX):
            os.symlink(HG38, BWA_IDX)
        print(f"  building the bwa index of {HG38} (once)", flush=True)
        start = time.time()
        subprocess.run(
            [BWA, "index", "-a", "bwtsw", BWA_IDX], check=True, capture_output=True
        )
        Path(BWA_IDX + ".index_s").write_text(f"{time.time() - start:.1f}")
    return BWA_IDX


def bwa_index_s(prefix):
    """Wall seconds of the bwa index build stored beside `prefix`, or None."""
    index_s_path = prefix + ".index_s"
    return float(Path(index_s_path).read_text()) if os.path.exists(index_s_path) else None


def primaries(sam):
    """Primary records of a SAM, those neither secondary nor supplementary."""
    result = subprocess.run(
        [SAMTOOLS or "samtools", "view", "-c", "-F", "0x900", sam],
        capture_output=True,
        text=True,
        check=False,
    )
    return int(result.stdout) if result.returncode == 0 and result.stdout.strip() else UNKNOWN_COUNT


def seq_bytes(sam):
    """Total bases stored in SEQ over all records.
    A hard-clipped supplementary record stores only its aligned part."""
    total = 0
    with open(sam) as handle:
        for line in handle:
            if line[0] != "@":
                seq_field = line.split("\t", SamCol.QUAL)[SamCol.SEQ]
                total += 0 if seq_field == "*" else len(seq_field)
    return total


def count_records(path):
    """Number of FASTA or FASTQ records in path."""
    return sum(1 for _ in _records(path))


WGSIM_ROW = re.compile(
    r"\s+(\S+)\s+(pe|se)\s+MAPQ>=30\s+(\d+) wrong\s+(\d+) \| (.*?)\s*\|"
    r"(?: wall ([\d.]+) s rss (\d+) MB(?: sd (\S+) runs (\d+) warmup (\d+))?)?\s*$"
)

WGSIM_FAIL = re.compile(r"\s+(\S+)\s+(pe|se)\s+(FAILED|CAPPED) exit (-?\d+): (.*)$")
# capture groups of WGSIM_ROW
WGSIM_ROW_TOOL = 1
WGSIM_ROW_ENDS = 2
WGSIM_ROW_MAPQ30 = 3
WGSIM_ROW_MAPQ30_WRONG = 4
WGSIM_ROW_DECILES = 5
WGSIM_ROW_WALL = 6
WGSIM_ROW_RSS = 7
WGSIM_ROW_SD = 8
WGSIM_ROW_RUNS = 9
WGSIM_ROW_WARMUP = 10
# capture groups of WGSIM_FAIL
WGSIM_FAIL_TOOL = 1
WGSIM_FAIL_ENDS = 2
WGSIM_FAIL_KIND = 3
WGSIM_FAIL_EXIT = 4
WGSIM_FAIL_MESSAGE = 5

# k and w of the minimap2 presets the map datasets run, from minimap2 v2.31 options.c.
# mm_idxopt_init sets k = 15, w = 10 for map-ont, the default.
# map-hifi and the asm* presets use k = 19, w = 19 (asm20 has w = 10), and sr uses k = 21, w = 11.
MM2_PRESET_KW = {
    "map-ont": (15, 10),
    "map-hifi": (19, 19),
    "asm20": (19, 10),
    "sr": (21, 11),
}
