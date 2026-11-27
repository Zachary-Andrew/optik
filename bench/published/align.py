"""Runs the pairwise-alignment section, which pits optik's WFA against WFA2-lib and the other aligners."""

import os
import statistics
import subprocess
import tarfile
from pathlib import Path

import compare_tools_small_sets
from compare_tools_small_sets import parse_wfa_bench_summary
from published_catalog import DATA_ROOT, TEMP_DIR

from published.common import (
    BIWFA_ONT_PAIRS,
    NO_RUN_EXIT,
    add_status,
    failure,
)

# The two exact methods whose scores must agree pair for pair.
EXACT_METHODS = ("optik-biwfa", "wfa2-ultralow")
# A pair file holds a pattern line, then its text line.
LINES_PER_PAIR = 2


def bench_align(bench, ds_id, dataset, mode, work):
    """Time optik align and BiWFA on the pairs of dataset `i` and check that their scores agree.

    biwfa-ont-ul48 unpacks the BiWFA paper's 48 ultra-long ONT pairs.
    The other sets are made with WFA2-lib's generate_dataset over the dataset's (length, error) grid."""
    align_dir = os.path.join(work, f"{ds_id}")
    os.makedirs(align_dir, exist_ok=True)
    exe = compare_tools_small_sets.build_wfab(align_dir)
    sets = []
    if ds_id == "biwfa-ont-ul48":
        tsv = os.path.join(align_dir, "pairs.tsv")
        if not os.path.exists(tsv):
            with (
                tarfile.open(
                    os.path.join(DATA_ROOT.path, ds_id, "ONT_MinION_UL.500kbps.48.tar.gz")
                ) as tar,
                open(tsv + ".part", "w") as handle,
            ):
                for member in sorted(tar.getmembers(), key=lambda member_info: member_info.name):
                    if member.isfile() and member.name.endswith(".seq"):
                        lines = tar.extractfile(member).read().decode().split("\n")
                        pattern = (
                            next(line for line in lines if line.startswith(">"))[1:]
                            .strip()
                            .upper()
                        )
                        text = (
                            next(line for line in lines if line.startswith("<"))[1:]
                            .strip()
                            .upper()
                        )
                        handle.write(f"{pattern}\t{text}\n")
            os.replace(tsv + ".part", tsv)
        sets.append(("ont_ul_500k_48", tsv, BIWFA_ONT_PAIRS))
    else:
        gen = os.path.join(TEMP_DIR, "WFA2-lib/bin/generate_dataset")
        for read_len, err_rate in dataset["grid"]:
            tsv = os.path.join(align_dir, f"sim_l{read_len}_e{err_rate}.tsv")
            if not os.path.exists(tsv):
                seq = os.path.join(align_dir, f"sim_l{read_len}_e{err_rate}.seq")
                subprocess.run(
                    [gen, "-o", seq, "-n", "1", "-l", str(read_len), "-e", str(err_rate)],
                    check=True,
                )
                lines = Path(seq).read_text().split("\n")
                with open(tsv + ".part", "w") as handle:
                    for pat_line, txt_line in zip(lines[0::LINES_PER_PAIR], lines[1::LINES_PER_PAIR], strict=False):
                        if pat_line.startswith(">") and txt_line.startswith("<"):
                            handle.write(f"{pat_line[1:]}\t{txt_line[1:]}\n")
                os.replace(tsv + ".part", tsv)
            sets.append((f"sim_l{read_len}_e{err_rate}", tsv, 1))
    for name, tsv, want in sets:
        with open(tsv) as handle:
            npairs = sum(1 for _ in handle)
        scores = {}
        # Only BiWFA runs here.
        # These pairs of 100 kbp to 2 Mbp are the BiWFA paper's long-sequence case, where full-memory WFA needs memory quadratic in the score, about 1 TB at 1 Mbp and 10% error, and cannot run on any machine.
        # Full-memory WFA is compared where it is a real choice, on the pairs up to 100 kbp in compare_tools_small_sets.py and on optik map's gaps of 1 kb or less.
        for method in EXACT_METHODS:
            out = os.path.join(align_dir, f"scores_{method}_{name}.txt")
            res = os.path.join(align_dir, "res.txt")
            # The harness times its own alignments, which is BiWFA's measure and leaves out parsing.
            # Warmups are dropped, and the statistics come from the --runs.
            runs_ = bench.repeat(
                f"{mode}/{ds_id}/{name}/{method}", [exe, method, tsv, out], stdout=res
            )
            rc_, rss, err_rate = runs_[-1] if runs_ else (NO_RUN_EXIT, 0.0, "no run")
            rss = max(run[1] for run in runs_) if runs_ else 0.0
            if rc_ != 0:
                kind, det = failure(bench, rc_, rss, err_rate)
                if kind == "cap":
                    add_status(
                        bench,
                        "cap",
                        section="align",
                        mode=mode,
                        dataset=ds_id,
                        set=name,
                        method=method,
                        gate_detail=det,
                    )
                else:
                    bench.add(
                        section="align",
                        mode=mode,
                        dataset=ds_id,
                        set=name,
                        method=method,
                        gate=False,
                        gate_detail=det,
                    )
                continue
            summaries = [parse_wfa_bench_summary(output) for _rc, _rss, output in runs_]
            secs_all = [summary["align_seconds"] for summary in summaries]
            invalid, aligned = (
                max(summary["invalid"] for summary in summaries),
                min(summary["pairs"] for summary in summaries),
            )
            secs = sum(secs_all) / len(secs_all)
            secs_sd = statistics.stdev(secs_all) if len(secs_all) > 1 else None
            scores[method] = Path(out).read_text().split()
            det = [
                reason
                for reason in (
                    "" if invalid == 0 else f"{invalid} CIGARs do not re-score",
                    ""
                    if npairs == want == aligned == len(scores[method])
                    else f"{aligned} pairs aligned, {len(scores[member])} scores, {npairs} in the set, "
                    f"{want} expected",
                )
                if reason
            ]
            bench.add(
                section="align",
                mode=mode,
                dataset=ds_id,
                set=name,
                method=method,
                pairs=npairs,
                align_s=round(secs, 3),
                align_s_sd=None if secs_sd is None else round(secs_sd, 3),
                runs=len(secs_all),
                warmup=bench.warmup,
                rss_mb=round(rss, 1),
                gate=not det,
                gate_detail="; ".join(det),
            )
        # The two exact methods must agree pair for pair.
        # If either is missing there is nothing to cross-check, and that fails as well.
        ref = next(iter(scores.values()), None)
        why = (
            "scores differ across exact methods"
            if any(score != ref for score in scores.values())
            else ""
            if len(scores) == len(EXACT_METHODS)
            else "only one method ran: scores not cross-checked"
        )
        if why:
            for row in bench.rows:
                if (
                    row.get("dataset") == ds_id
                    and row.get("set") == name
                    and not row.get("status")
                ):
                    row["gate"] = False
                    row["gate_detail"] = (
                        (row.get("gate_detail") or "") + "; " + why
                    ).lstrip("; ")
