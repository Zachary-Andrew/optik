"""Runs the MAPQ calibration section."""

import os
import re
import subprocess
import sys

from diagnostics import DETAIL_CHARS
from small_common import TEMP_DIR

from small_sets.common import HERE


def section_mapq(bench, mode, work, threads, quick):
    """Run bench/mapq_calibration.py at several error rates and report over-confident MAPQ bins."""
    print(f"\n=== mapq ({mode})", flush=True)
    mapq_dir = os.path.join(work, "mapq")
    for args in (
        ["--error", "0.01"],
        ["--error", "0.04"],
        ["--error", "0.08"],
        ["--ref", os.path.join(TEMP_DIR, "salmonella_LT2.fa"), "--error", "0.04"],
    ):
        if quick and args[1] != "0.04":
            continue
        num_reads = "200000" if quick else "1000000"
        proc = subprocess.run(
            [
                sys.executable,
                os.path.join(HERE, "mapq_calibration.py"),
                # Its rows join this bench's section, so it writes no section of its own.
                "--report",
                "",
                "--reads",
                num_reads,
                "--threads",
                str(threads),
                "--out",
                mapq_dir,
            ]
            + args,
            capture_output=True,
            text=True,
            check=False,
        )
        bench.log(f"{mode}/mapq", ["mapq_calibration.py", "--reads", num_reads] + args)
        txt = proc.stdout
        for tool in ("optik", "minimap2"):
            match = re.search(
                rf"^{tool}:.*?MAPQ>=30: (?P<reads>\d+) reads, (?P<wrong>\d+) wrong",
                txt,
                re.DOTALL | re.MULTILINE,
            )
            block = (
                txt.split(f"{tool}:")[1].split("minimap2:")[0]
                if f"{tool}:" in txt
                else ""
            )
            over = block.count("OVERCONFIDENT")
            # A missing report, whether the calibration crashed or its format changed, is a failure, not "0 overconfident bins".
            parsed = proc.returncode == 0 and match is not None and bool(block)
            bench.add(
                section="mapq",
                mode=mode,
                setting=" ".join(arg.replace(TEMP_DIR + "/", "") for arg in args),
                tool=tool,
                mapq30_reads=int(match.group("reads")) if match else None,
                mapq30_wrong=int(match.group("wrong")) if match else None,
                overconfident_bins=over,
                calibrated=(
                    "yes"
                    if over == 0
                    else f"no: {over} bins over-confident (wrong > 2x implied + 3, Poisson tail < 1e-3)"
                )
                if parsed
                else "n/a",
                # The gate only requires that the calibration ran and was parsed.
                # Whether a mapper is calibrated is the result, reported the same way for every tool.
                gate=parsed,
                gate_detail=""
                if parsed
                else (
                    f"mapq_calibration.py produced no report (rc {proc.returncode}): "
                    f"{proc.stderr.strip()[-DETAIL_CHARS:]}"
                ),
            )
