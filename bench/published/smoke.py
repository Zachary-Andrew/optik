"""Runs every tool on the first few megabases of each dataset, so formats and gates can be checked quickly."""

import logging
import os

import compare_tools_small_sets
import report as results_file
from diagnostics import DETAIL_CHARS
from published_catalog import CATALOG, DATA_ROOT, TEMP_DIR
from published_inputs import sshash_m
from published_tools import OWN_OUTPUT, count_of, query_tools
from units import BYTES_PER_GIB

from published.common import (
    FETCH_ERRORS,
    ROOT,
    SMOKE_BASES,
)
from published.dbg import bench_dbg
from published.fetch import fetch_prefix
from published.query import query_inputs
from published.status import gate_text

LOG = logging.getLogger(__name__)


# The memory cap leaves this much of the machine to the system and the bench itself.
CAP_HEADROOM_GIB = 2


def smoke(args, mach):
    """Run every tool of the dataset's section on the first 2 Mbp or so of the formatted inputs.

    The formats are right when every tool runs and the gates pass on that slice.
    Nothing is written to the results directory.
    Returns 1 when any gate fails, else 0."""
    out = os.path.join(TEMP_DIR, "published_sets_smoke")
    os.makedirs(out, exist_ok=True)
    bench = compare_tools_small_sets.Bench(
        out, 1, 0, max(1, int(mach["mem_gib"] - CAP_HEADROOM_GIB)), os.path.join(out, "_scratch")
    )
    only = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if args.skip else set()
    nfail = 0
    checks = []  # (dataset, check, passed, detail) of the query rounds
    for ds_id, dataset in CATALOG.items():
        if (
            "excluded" in dataset
            or ("opt_in" in dataset and not only)
            or (only and ds_id not in only)
            or ds_id in skip
            or (args.section and dataset["section"] not in args.section.split(","))
        ):
            continue
        if dataset["section"] not in ("query", "dbg") or dataset.get("local"):
            continue
        # A scaling set with no files of its own runs on its source set's data (`as_`).
        src_id = dataset["as_"] if dataset.get("as_") and not dataset["files"] else ds_id
        src = os.path.join(DATA_ROOT.path, src_id)
        fetched = os.path.exists(os.path.join(src, "PROVENANCE.json"))
        if not fetched:
            if args.no_download:
                continue
            if any(
                file_spec.get("key", file_spec.get("url", "")).endswith((".tar.gz", ".zip"))
                or "ncbi" in file_spec
                for file_spec in CATALOG[src_id]["files"]
            ):
                print(
                    f"\n=== smoke {dataset['section']} {ds_id}: an archive, which needs the full download "
                    f"(fetch --only {ds_id}, then smoke)"
                )
                continue
            try:
                fetch_prefix(src_id, CATALOG[src_id])
            # A source outage skips this set, not the whole run.
            except FETCH_ERRORS as fetch_err:
                print(
                    f"\n=== smoke {dataset['section']} {ds_id}: prefix download failed ({str(fetch_err)[:DETAIL_CHARS]})"
                )
                continue
        saved = DATA_ROOT.path
        if not fetched:
            # The prefixes stand in for the files under the files' own names.
            DATA_ROOT.path = os.path.join(TEMP_DIR, "published_sets_smoke_data")
            os.makedirs(os.path.join(DATA_ROOT.path, src_id), exist_ok=True)
            for fname in (
                os.listdir(os.path.join(src, "_prefix"))
                if os.path.isdir(os.path.join(src, "_prefix"))
                else []
            ):
                lnk = os.path.join(DATA_ROOT.path, src_id, fname)
                if not os.path.exists(lnk):
                    os.symlink(os.path.join(src, "_prefix", fname), lnk)
        print(
            f"\n=== smoke {dataset['section']} {ds_id} ({'fetched' if fetched else 'prefix'})",
            flush=True,
        )
        try:
            if dataset["section"] == "query":
                kmer_len = dataset["k"][-1]
                unitigs_fa, built_by, workloads, windows = query_inputs(
                    ds_id,
                    dataset,
                    kmer_len,
                    "full",
                    os.path.join(out, "work"),
                    args.threads,
                    limit_bases=SMOKE_BASES,
                )
                workloads.pop("_load")
                query_dir = os.path.join(out, "work", f"{ds_id}_k{kmer_len}_full_smoke", "optimal")
                os.makedirs(query_dir, exist_ok=True)
                tools = query_tools(
                    "optimal",
                    query_dir,
                    unitigs_fa,
                    kmer_len,
                    sshash_m(dataset, kmer_len, windows),
                    args.threads,
                    mach,
                    bench.cap // BYTES_PER_GIB,
                )
                counts = {}
                for tool, spec in tools.items():
                    if (
                        spec["build"] is not None
                    ):  # None marks a query mode of another tool's index
                        rc_, _wall, _rss, err = bench.once(
                            f"smoke/{ds_id}/{tool}/build", spec["build"], cwd=query_dir
                        )
                        if rc_ != 0:
                            print(f"  [FAIL] {tool} build: {err.strip()[-DETAIL_CHARS:]}")
                            checks.append(
                                (ds_id, f"{tool} build", False, err.strip()[-DETAIL_CHARS:])
                            )
                            nfail += 1
                            continue
                    for workload, query_path in workloads.items():
                        cmd, outf = spec["query"](query_path)
                        capt = (
                            os.path.join(query_dir, "q.out")
                            if tool in OWN_OUTPUT
                            else (outf or os.path.join(query_dir, "q.out"))
                        )
                        rc_, _wall, _rss, err = bench.once(
                            f"smoke/{ds_id}/{tool}/{workload}",
                            cmd,
                            cwd=query_dir,
                            stdout=None if tool == "sbwt" else capt,
                        )
                        if rc_ != 0:
                            print(f"  [FAIL] {tool} {workload}: {err.strip()[-DETAIL_CHARS:]}")
                            checks.append((ds_id, f"{tool} {workload}", False, err.strip()[-DETAIL_CHARS:]))
                            nfail += 1
                            continue
                        counts.setdefault(workload, {})[tool] = count_of(
                            tool, outf, capt, err
                        )[1]
                for workload, tool_counts in counts.items():
                    same = len(set(tool_counts.values())) == 1
                    print(
                        f"  [{'PASS' if same else 'FAIL'}] {workload}: positives {tool_counts}  (SPSS: {built_by}, "
                        f"{windows} k-mer windows)"
                    )
                    checks.append((ds_id, f"{workload} positives", same, str(tool_counts)))
                    nfail += not same
            else:
                bench_dbg(
                    bench,
                    ds_id,
                    dataset,
                    "optimal",
                    os.path.join(out, "work"),
                    args.threads,
                    mach,
                    limit_bases=SMOKE_BASES,
                )
        except Exception as fetch_err:
            # Any failure of one dataset, a bench bug included, is reported with its traceback and the smoke run goes on.
            LOG.exception("  [FAIL] %s", ds_id)
            checks.append((ds_id, "run", False, f"{type(fetch_err).__name__}: {str(fetch_err)[-DETAIL_CHARS:]}"))
            nfail += 1
        finally:
            DATA_ROOT.path = saved
    # The graph rounds record their gates as rows.
    # SKIP and CAP rows have no gate.
    nfail += sum(1 for row in bench.rows if row.get("gate") is False)
    print(f"\nsmoke: {nfail} gate failures")
    hidden = ("section", "params", "gate_detail", "status", "gate")
    results_file.write_section(
        args.report,
        "published_smoke",
        "Functional check of every published dataset on a prefix of its input: each tool builds "
        "and queries, query tools agree on positives, and graph outputs pass their gates",
        [
            (
                "query",
                [
                    {
                        "dataset": ds_id,
                        "check": check_name,
                        "gate": "PASS" if passed else "FAIL",
                        "detail": det,
                    }
                    for ds_id, check_name, passed, det in checks
                ],
            ),
            (
                "dbg",
                [
                    {
                        **{key: value for key, value in row.items() if key not in hidden},
                        "gate": gate_text(row).strip(),
                    }
                    for row in bench.rows
                ],
            ),
        ],
        notes=[f"Logs: `{os.path.relpath(out, ROOT)}`"],
    )
    return 1 if nfail else 0
