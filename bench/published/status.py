"""Writes the results report and gate summaries, and implements the status and list commands."""

import json
import os
import time
from pathlib import Path

import compare_tools_small_sets
import report as results_file
from published_catalog import CATALOG, DATA_ROOT, DBG_TOOLS, QUERY_TOOLS
from units import BYTES_PER_GB, BYTES_PER_GIB, BYTES_PER_MB, PERCENT, SECONDS_PER_MINUTE

from published.common import (
    _SIZES,
    PROGRESS_FRESH_S,
    PROTOCOL,
    REFERENCE_K,
    ROOT,
    SECTION_TITLES,
)
from published.machine import bytes_per_kmer, download_gb, machine, preflight, tool_skip


def gate_text(row):
    """PASS, FAIL or the status (SKIP, CAP) of one row, with its detail when it did not pass."""
    if row.get("status"):
        return row["status"].upper() + " " + str(row.get("gate_detail", ""))
    return "PASS" if row.get("gate") else "FAIL " + str(row.get("gate_detail", ""))


def report(bench, out, mach, skipped, dest=""):
    """Write results.json and SUMMARY.md of a bench run to `out`, and one section per paper section of the sectioned results file `dest`."""
    Path(os.path.join(out, "results.json")).write_text(
        json.dumps(
            {"machine": mach, "skipped": skipped, "rows": bench.rows},
            indent=1,
            default=str,
        )
    )
    md_lines = [
        "# optik vs published tools, on the tools' own published datasets",
        "",
        (
            "bench/compare_tools_published_sets.py; provenance of each dataset in its PROVENANCE.json under "
            f"{os.path.relpath(DATA_ROOT.path, ROOT)}/."
        ),
        "",
        (
            f"Machine: {mach['mem_gib']} GiB usable memory (per-process cap {bench.cap // BYTES_PER_GIB} GiB), "
            f"{mach['disk_free_gb']} GB free disk. PASS / FAIL are correctness gates; SKIP is a tool whose "
            "reported or estimated memory need exceeds this machine; CAP is a run that died against the cap."
        ),
        "",
    ]
    for sec in ("query", "dbg", "map", "align"):
        rows = [row for row in bench.rows if row.get("section") == sec]
        if not rows:
            continue
        # A section's columns are the union of its rows' keys in first-seen order, because different tools report different metrics.
        keys = []
        for row in rows:
            keys += [
                col
                for col in row
                if col not in keys
                and col not in ("section", "params", "gate_detail", "status")
            ]
        md_lines += [
            f"## {sec}",
            "",
            "| " + " | ".join(keys) + " | gate |",
            "|" + "---|" * (len(keys) + 1),
        ]
        for row in rows:
            status_text = gate_text(row)
            md_lines.append(
                "| "
                + " | ".join(compare_tools_small_sets._fmt(row.get(col, "")) for col in keys)
                + " | "
                + status_text
                + " |"
            )
        md_lines.append("")
    if skipped:
        md_lines += ["## Not run on this machine", ""]
        md_lines += [
            f"- `{ds_id}` ({CATALOG[ds_id].get('title')}): {why}" for ds_id, why in skipped.items()
        ]
        md_lines.append("")
    excl = [(ds_id, dataset) for ds_id, dataset in CATALOG.items() if "excluded" in dataset]
    md_lines += ["## Not runnable anywhere (from the papers, with the reason)", ""]
    md_lines += [
        f"- `{ds_id}` ({dataset.get('title')}; {dataset['paper']}): {dataset['excluded']}" for ds_id, dataset in excl
    ]
    Path(os.path.join(out, "SUMMARY.md")).write_text("\n".join(md_lines) + "\n")
    # The sectioned results file gets the same rows with the bookkeeping fields dropped and the gate reduced to its text.
    hidden = ("section", "params", "gate_detail", "status", "gate")
    for sec, title in SECTION_TITLES.items():
        rows = [row for row in bench.rows if row.get("section") == sec]
        if not rows:
            continue
        flat = [
            {
                **{key: value for key, value in row.items() if key not in hidden},
                "gate": gate_text(row).strip(),
            }
            for row in rows
        ]
        notes = [
            f"- `{ds_id}` not run here: {why}"
            for ds_id, why in skipped.items()
            if CATALOG[ds_id].get("section") == sec
        ]
        results_file.write_section(
            dest,
            f"published_{sec}",
            title + f" (per-process cap {bench.cap // BYTES_PER_GIB} GiB)",
            [("", flat)],
            notes=notes
            + [f"Raw rows, parameters and logs: `{os.path.relpath(out, ROOT)}`"],
        )


def pid_alive(pid):
    """True when process `pid` exists, so a fetch killed seconds ago is not shown as running.

    A progress file without a pid, written before pids were recorded, counts as alive."""
    if pid is None:
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # it exists, under another user
    return True


def part_done(part):
    """Bytes of an in-progress download that are already on disk.

    These are the completed ranges the range scheduler recorded in <part>.done.
    Failing that, it is the size of a plain curl .part, since curl writes from the start."""
    done_file = part + ".done"
    if os.path.exists(done_file):
        try:
            return sum(end - start for start, end in json.loads(Path(done_file).read_text()))
        except (OSError, ValueError):
            pass
    return os.path.getsize(part) if os.path.exists(part) else 0


def cmd_status(args):
    """Print every dataset's download state.

    For each dataset this shows the bytes its sources publish (from the size cache that fetch fills) and the bytes on disk, counting finished files and the written part of files in progress.
    It also shows whether the dataset is verified (PROVENANCE.json appears once all its files have passed their checksums) and the fetch's measured rate and ETA, taken from DATA_ROOT.path/_fetch_progress.json, which the range scheduler rewrites every few seconds during a fetch."""
    while True:
        sizes = json.loads(Path(_SIZES).read_text()) if os.path.exists(_SIZES) else {}
        prog = os.path.join(DATA_ROOT.path, "_fetch_progress.json")
        progress = None
        if os.path.exists(prog):
            try:
                progress = json.loads(Path(prog).read_text())
            except (OSError, ValueError):
                progress = None
        live = (
            progress is not None
            and time.time() - progress.get("time", 0) < PROGRESS_FRESH_S
            and pid_alive(progress.get("pid"))
        )
        rows, left_all = [], 0.0
        for ds_id, dataset in CATALOG.items():
            if "excluded" in dataset or not dataset.get("files"):
                continue
            want = sizes.get(ds_id)
            out = os.path.join(DATA_ROOT.path, ds_id)
            have = 0
            parts = 0
            if os.path.isdir(out):
                for fname in os.listdir(out):
                    pth = os.path.join(out, fname)
                    if fname.endswith(".part"):
                        have += part_done(pth)
                        parts += 1
                    elif (
                        os.path.isfile(pth)
                        and not fname.endswith((".done", ".json", ".tmp"))
                        and not fname.startswith("_")
                    ):
                        have += os.path.getsize(pth)
            verified = os.path.exists(os.path.join(out, "PROVENANCE.json"))
            if verified:
                # Count the downloaded files themselves, not what was unpacked or derived from them.
                recs = json.loads(
                    Path(os.path.join(out, "PROVENANCE.json")).read_text()
                ).get("files", [])
                have = sum(
                    os.path.getsize(rec["path"])
                    for rec in recs
                    if not rec.get("local")
                    and rec.get("path")
                    and os.path.exists(rec["path"])
                )
                state = "verified"
            elif parts:
                state = "downloading" if live else "partial (fetch not running)"
            elif have:
                state = "partial"
            else:
                state = "-"
            pct = (
                "100%"
                if verified
                else (f"{min(float(PERCENT), PERCENT * have / (want * BYTES_PER_GB)):.1f}%" if want else "?")
            )
            if "opt_in" in dataset and not verified:
                state = "opt-in, not fetched by default" + (
                    " (partial)" if have else ""
                )
            elif want and not verified:
                left_all += max(0.0, want * BYTES_PER_GB - have)
            rows.append((ds_id, dataset["section"], want, have / BYTES_PER_GB, pct, state))
        print(
            f"{'dataset':28s} {'section':7s} {'source GB':>9s} {'on disk GB':>10s} {'done':>7s}  state"
        )
        for ds_id, sec, want, have, pct, state in rows:
            print(
                f"{ds_id:28s} {sec:7s} {('?' if want is None else f'{want:.2f}'):>9s} {have:10.2f} {pct:>7s}  {state}"
            )
        if live:
            rate = progress.get("rate_bytes_s", 0.0)
            infl = sum(file_state["left"] for file_state in progress.get("files", []))
            eta = (
                f"{infl / rate / SECONDS_PER_MINUTE:.0f} min for the files in flight"
                if rate > 0
                else "rate not measured yet"
            )
            print(
                f"\nfetch running: {rate / BYTES_PER_MB:.1f} MB/s over {progress.get('connections')} connections "
                f"(target {progress.get('target')}); {infl / BYTES_PER_GB:.1f} GB left in flight; ETA {eta}"
            )
        else:
            print(
                "\nno fetch running"
                + (f" (last progress {time.ctime(progress['time'])})" if progress else "")
            )
        print(
            f"not yet fetched and verified: {left_all / BYTES_PER_GB:.1f} GB over all listed datasets "
            "(fetch takes the ones this machine has disk for: `fetch --tier auto`)"
        )
        if not args.watch:
            return
        time.sleep(args.watch)
        print()


def cmd_list(args):
    """Print the dataset catalogue by section, with each dataset's protocol and whether it runs here."""
    mach = machine()
    _, src = bytes_per_kmer()
    print(
        f"this machine: {mach['mem_gib']} GiB memory, {mach['disk_free_gb']} GB free disk; "
        f"memory estimates from {src}\n"
    )
    if args.md:
        print(
            "| id | what | section | paper | source | use | here |\n|---|---|---|---|---|---|---|"
        )
    # Each section lists its datasets with the paper protocol, and says for each whether it can run on this machine and why not.
    for sec in ("query", "dbg", "map", "align", "count"):
        if not args.md:
            print(f"== {sec}")
            protocol = PROTOCOL.get(sec)
            if protocol:
                print(
                    f"  protocol: {protocol['papers']}\n  metrics:  {protocol['metrics']}\n  tools:    {protocol['tools']}\n"
                )
        for ds_id, dataset in CATALOG.items():
            if dataset["section"] != sec:
                continue
            use = dataset.get("excluded") and "EXCLUDED: " + dataset["excluded"] or dataset.get("use", "")
            if dataset.get("note"):
                use += " NOTE: " + dataset["note"]
            # Excluded and opt-in datasets are reported as such, and every other one goes through the same preflight a bench run would use.
            if "excluded" in dataset:
                here = "excluded"
            elif "opt_in" in dataset:
                here = dataset["opt_in"]
            else:
                runnable, why = preflight(ds_id, dataset, mach)
                here = "runnable" if runnable else "; ".join(why)
                download_size = download_gb(ds_id, dataset)
                if download_size:
                    here += f" ({download_size:.1f} GB download)"
                kmers = dataset.get("kmers31")
                # Tools whose memory need exceeds this machine are named, even when the dataset itself can run.
                tool_skips = [
                    f"{tool}: {reason}"
                    for tool in (
                        QUERY_TOOLS
                        if sec == "query"
                        else DBG_TOOLS
                        if sec == "dbg"
                        else ()
                    )
                    if (
                        reason := tool_skip(
                            dataset,
                            tool,
                            mach,
                            REFERENCE_K if sec == "query" else dataset.get("k", [REFERENCE_K])[0],
                            kmers,
                        )
                    )
                ]
                if tool_skips:
                    here += "; skips " + "; ".join(tool_skips)
            if args.md:
                print(
                    f"| `{ds_id}` | {dataset.get('title')} | {sec} | {dataset['paper']} | {dataset['source']} | {use} | {here} |"
                )
            else:
                print(
                    f"  {ds_id:24} {dataset.get('title')}\n  {'':24} {dataset['paper']}\n  {'':24} {dataset['source']}\n"
                    f"  {'':24} {use}\n  {'':24} here: {here}\n"
                )
