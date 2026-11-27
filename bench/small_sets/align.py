"""Runs the pairwise-alignment section, timing optik's WFA against WFA2-lib's wfa_affine_bench."""

import os
import random
import statistics
from pathlib import Path

from diagnostics import DETAIL_CHARS
from dna import mutate
from units import BYTES_PER_GIB, US_PER_SECOND

from small_sets.common import (
    ALIGN_GRID,
    ALIGN_SEED_ERR_SCALE,
    ALIGN_SEED_LEN_MULT,
    DEL_SHARE,
    INS_SHARE,
    NO_RUN_EXIT,
    QUICK_ALIGN_CELLS,
    SUB_SHARE,
    build_wfab,
)

# The summary line alternates labels and values.
LABEL_AND_VALUE = 2


def parse_wfa_bench_summary(line):
    """The figures of wfa_affine_bench's summary line, "align_seconds S invalid N pairs P", keyed by their labels."""
    tokens = line.split()
    labels, values = tokens[0::LABEL_AND_VALUE], tokens[1::LABEL_AND_VALUE]
    return {label: float(value) if label == "align_seconds" else int(value) for label, value in zip(labels, values)}


def section_align(bench, mode, work, quick):
    """Align random pairs over a length and error grid with optik and WFA2-lib and compare scores."""
    print(f"\n=== align ({mode})", flush=True)
    align_dir = os.path.join(work, "align")
    os.makedirs(align_dir, exist_ok=True)
    exe = build_wfab(align_dir)
    grid = ALIGN_GRID
    if quick:
        grid = grid[:QUICK_ALIGN_CELLS]
    methods = ["optik-wfa", "optik-biwfa", "wfa2-high", "wfa2-ultralow"]
    for read_len, err, npairs in grid:
        rng = random.Random(read_len * ALIGN_SEED_LEN_MULT + int(err * ALIGN_SEED_ERR_SCALE))
        pairs_tsv = os.path.join(align_dir, f"pairs_{read_len}_{err}.tsv")
        if not os.path.exists(pairs_tsv):
            with open(pairs_tsv, "w") as handle:
                for _ in range(npairs):
                    seq = "".join(rng.choice("ACGT") for _ in range(read_len))
                    handle.write(
                        seq
                        + "\t"
                        + mutate(seq, rng, err * SUB_SHARE, err * INS_SHARE, err * DEL_SHARE)
                        + "\n"
                    )
        scores = {}
        for method in methods:
            out = os.path.join(align_dir, f"scores_{method}_{read_len}_{err}.txt")
            res = os.path.join(align_dir, "res.txt")
            # The harness times its own alignments, BiWFA's measure, which leaves out input parsing.
            # Warmup runs are dropped from the statistics.
            runs_ = bench.repeat(
                f"{mode}/align/{read_len}/{err}/{method}", [exe, method, pairs_tsv, out], stdout=res
            )
            if not runs_ or runs_[-1][0] != 0:
                rc_, _, err_out = runs_[-1] if runs_ else (NO_RUN_EXIT, 0, "no run")
                bench.add(
                    section="align",
                    mode=mode,
                    length=read_len,
                    error=err,
                    method=method,
                    gate=False,
                    gate_detail=f"exit {rc_} (memory cap {bench.cap // BYTES_PER_GIB} GiB?): {err_out[-DETAIL_CHARS:]}",
                )
                continue
            summaries = [parse_wfa_bench_summary(output) for _rc, _rss, output in runs_]
            secs = [summary["align_seconds"] for summary in summaries]
            invalid = max(summary["invalid"] for summary in summaries)
            scores[method] = Path(out).read_text().split()
            us_per_pair = [secs_val * US_PER_SECOND / npairs for secs_val in secs]
            mean = sum(us_per_pair) / len(us_per_pair)
            us_sd = statistics.stdev(us_per_pair) if len(us_per_pair) > 1 else None
            bench.add(
                section="align",
                mode=mode,
                length=read_len,
                error=err,
                method=method,
                pairs=npairs,
                us_per_pair=round(mean, 2),
                us_per_pair_sd=round(us_sd, 2) if us_sd is not None else None,
                runs=len(us_per_pair),
                warmup=bench.warmup,
                rss_mb=round(max(run[1] for run in runs_), 1),
                gate=invalid == 0,
                gate_detail="" if invalid == 0 else f"{invalid} CIGARs do not re-score",
            )
        ref = next(iter(scores.values()), None)
        agree = all(score == ref for score in scores.values())
        for row in bench.rows:
            if (
                row.get("section") == "align"
                and row.get("length") == read_len
                and row.get("error") == err
                and row.get("mode") == mode
            ):
                row["scores_agree"] = agree
                if not agree:
                    row["gate"] = False
                    row["gate_detail"] = (
                        row.get("gate_detail") or ""
                    ) + " optimal scores differ between methods"
