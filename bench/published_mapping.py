"""Mapping helpers of compare_tools_published_sets.py.

The module wraps the papers' own evaluators (paftools.js mapeval and call, wgsim_eval.pl), the prebuilt bwa index of GRCh38 and the small SAM summaries that the map section reports.
"""

import os
import re
import subprocess
import time
from pathlib import Path

from external_oracles import SAMTOOLS
from published_catalog import BWA, BWA_IDX, HG38, TOOLS, T
from published_inputs import _records

PAFTOOLS = [os.path.join(TOOLS, "k8"), os.path.join(T, "minimap2/misc/paftools.js")]


def mapeval_calibration(ev):
    """Over-confident MAPQ bins in paftools.js mapeval output, whose rows read 'Q mapq n wrong ...'.

    A bin is over-confident when wrong > 2 x implied + 3, where implied = n x 10^(-Q/10).
    compare_tools_small_sets and mapq_calibration.py use the same rule."""
    over = []
    for x in ev.splitlines():
        t = x.split("\t")
        if t[0] == "Q":
            q, n, w = int(t[1]), int(t[2]), int(t[3])
            implied = n * 10 ** (-q / 10)
            if w > 2 * implied + 3:
                over.append(f"Q{q}: {w} wrong vs {implied:.1f} implied")
    return (
        "yes"
        if not over
        else f"no: {len(over)} bins over-confident ({'; '.join(over[:3])})"
    )


def paftools_call_snvs(sam, ref, d, tag):
    """Count substitutions called by minimap2's assembly-variant caller in the regions one alignment covers.

    paftools.js call runs on the cs-tagged PAF of the primary and supplementary records.
    Returns (SNVs, bases covered once)."""
    md = os.path.join(d, f"{tag}.md.sam")
    with open(md, "w") as o:
        subprocess.run(
            [SAMTOOLS or "samtools", "calmd", sam, ref],
            stdout=o,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    paf = subprocess.run(
        PAFTOOLS + ["sam2paf", "-p", md], capture_output=True, text=True, check=True
    ).stdout
    rows = sorted(
        (x for x in paf.splitlines() if x),
        key=lambda x: (x.split("\t")[5], int(x.split("\t")[7])),
    )
    out = subprocess.run(
        PAFTOOLS + ["call", "-l", "1000", "-L", "1000", "-"],
        input="\n".join(rows) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    snv = cov = 0
    for x in out.splitlines():
        t = x.split("\t")
        if t[0] == "R":
            cov += int(t[3]) - int(t[2])
        elif t[0] == "V" and len(t[6]) == 1 and len(t[7]) == 1:
            snv += 1
    os.remove(md)
    return snv, cov


def index_ready_s(tool, err):
    """Seconds from process start until the mapper's index was ready, read from its own log.

    The log is the last measured run's stderr.
    optik writes "[map] ...  minimizers indexed (...) in S s", which covers reading the reference, the first read batch and the index build.
    minimap2 writes "[M::main::S*...] loaded/built the index".
    A mapper that reads a prebuilt index, as bwa does per bwa_index_s, gives None."""
    if tool == "optik":
        m = re.findall(r"minimizers indexed .*? in ([\d.]+)s", err)
    elif tool == "minimap2":
        m = re.findall(r"\[M::main::([\d.]+)\*[\d.]+\] loaded/built the index", err)
    else:
        return None
    return float(m[-1]) if m else None


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
        t0 = time.time()
        subprocess.run(
            [BWA, "index", "-a", "bwtsw", BWA_IDX], check=True, capture_output=True
        )
        Path(BWA_IDX + ".index_s").write_text(f"{time.time() - t0:.1f}")
    return BWA_IDX


def bwa_index_s(prefix):
    """Wall seconds of the bwa index build stored beside `prefix`, or None."""
    f = prefix + ".index_s"
    return float(Path(f).read_text()) if os.path.exists(f) else None


def primaries(sam):
    """Primary records of a SAM, those neither secondary nor supplementary."""
    r = subprocess.run(
        [SAMTOOLS or "samtools", "view", "-c", "-F", "0x900", sam],
        capture_output=True,
        text=True,
        check=False,
    )
    return int(r.stdout) if r.returncode == 0 and r.stdout.strip() else -1


def seq_bytes(sam):
    """Total bases stored in SEQ over all records.
    A hard-clipped supplementary record stores only its aligned part."""
    n = 0
    with open(sam) as f:
        for line in f:
            if line[0] != "@":
                s = line.split("\t", 10)[9]
                n += 0 if s == "*" else len(s)
    return n


def count_records(path):
    """Number of FASTA or FASTQ records in path."""
    return sum(1 for _ in _records(path))


WGSIM_ROW = re.compile(
    r"\s+(\S+)\s+(pe|se)\s+MAPQ>=30\s+(\d+) wrong\s+(\d+) \| (.*?)\s*\|"
    r"(?: wall ([\d.]+) s rss (\d+) MB(?: sd (\S+) runs (\d+) warmup (\d+))?)?\s*$"
)

WGSIM_FAIL = re.compile(r"\s+(\S+)\s+(pe|se)\s+(FAILED|CAPPED) exit (-?\d+): (.*)$")

# k and w of the minimap2 presets the map datasets run, from minimap2 v2.31 options.c.
# mm_idxopt_init sets k = 15, w = 10 for map-ont, the default.
# map-hifi and the asm* presets use k = 19, w = 19 (asm20 has w = 10), and sr uses k = 21, w = 11.
MM2_PRESET_KW = {
    "map-ont": (15, 10),
    "map-hifi": (19, 19),
    "asm20": (19, 10),
    "sr": (21, 11),
}
