#!/usr/bin/env python3
"""fetch_check.py - checks that the published-sets fetch works, resumes and reports.

All of it happens in a scratch data root (OPTIK_BENCH_DATA), never the real one.
  1. Fetch two small datasets (cf2-celegans-ref, fmsi25-spneumo).
     Both must end verified against their source's checksum, with PROVENANCE.json written.
  2. Start fetching a larger one (sshash26-ncbi-virus, about 140 MB), stop it once some byte ranges are on disk, and check that `status` reports it partial with those bytes.
  3. Fetch it again.
     The scheduler must resume from exactly the bytes on disk, as its "resuming ... N of M bytes" line shows, and the file must then verify against the published checksum.
  4. `status` reports all three verified.

Usage: python3 bench/fetch_check.py [--work DIR] The exit status is 0 when every step passes.
The network (curl) is needed.
"""

import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from diagnostics import LOG_TAIL_CHARS

# the partial fetch is stopped after this many seconds or once this many bytes are done
PARTIAL_FETCH_MAX_S = 600
PARTIAL_FETCH_BYTES = 20e6

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "compare_tools_published_sets.py")
SMALL = ["cf2-celegans-ref", "fmsi25-spneumo"]
BIG = "sshash26-ncbi-virus"


def run(env, *args, **kwargs):
    return subprocess.run(
        [sys.executable, "-u", BENCH, *args],
        env=env,
        capture_output=True,
        text=True,
        **kwargs,
        check=False,
    )


def done_bytes(data, dataset_id):
    """Bytes of dataset i in the completed ranges its .part.done files record."""
    dataset_dir = os.path.join(data, dataset_id)
    total = 0
    for name in os.listdir(dataset_dir) if os.path.isdir(dataset_dir) else []:
        if name.endswith(".part.done"):
            total += sum(range_end - range_start for range_start, range_end in json.loads(Path(os.path.join(dataset_dir, name)).read_text()))
    return total


def main():
    """Run the four fetch steps against a scratch data root."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--work", default=os.path.join(HERE, "..", "temp", "fetch_check"))
    args = parser.parse_args()
    data = os.path.abspath(args.work)
    shutil.rmtree(data, ignore_errors=True)
    os.makedirs(data)
    env = dict(os.environ, OPTIK_BENCH_DATA=data)
    bad = []

    def check(passed, what):
        """Print one step result and remember a failure."""
        print(f"  [{'PASS' if passed else 'FAIL'}] {what}", flush=True)
        if not passed:
            bad.append(what)

    print(f"scratch data root: {data}")
    # Step 1, small datasets
    result = run(env, "fetch", "--only", ",".join(SMALL), "--force")
    for dataset_id in SMALL:
        check(
            os.path.exists(os.path.join(data, dataset_id, "PROVENANCE.json")),
            f"{dataset_id} fetched and verified",
        )
    if result.returncode != 0:
        print(result.stdout[-LOG_TAIL_CHARS:] + result.stderr[-LOG_TAIL_CHARS:])
    # Step 2, a fetch stopped partway
    with open(os.path.join(data, "big1.log"), "w") as log:
        proc = subprocess.Popen(
            [sys.executable, "-u", BENCH, "fetch", "--only", BIG, "--force"],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        start = time.time()
        # The scheduler records completed ranges every 5 s, so wait for some but stop well before the end.
        while (
            time.time() - start < PARTIAL_FETCH_MAX_S and proc.poll() is None and done_bytes(data, BIG) < PARTIAL_FETCH_BYTES
        ):
            time.sleep(1)
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
    before = done_bytes(data, BIG)
    check(
        before > 0 and not os.path.exists(os.path.join(data, BIG, "PROVENANCE.json")),
        f"{BIG} stopped partway ({before} bytes recorded on disk)",
    )
    status_out = run(env, "status").stdout
    line = next((line for line in status_out.splitlines() if line.startswith(BIG)), "")
    check("partial" in line, f"status reports the stopped fetch: {line.strip()!r}")
    # Step 3, resume
    result = run(env, "fetch", "--only", BIG, "--force")
    resumes = re.findall(r"resuming \S+: (\d+) of (\d+) bytes already on disk", result.stdout)
    resumed = sum(int(done_str) for done_str, _ in resumes)
    check(
        resumed == before and before > 0,
        f"resumed from the {before} bytes on disk (scheduler: {resumed})",
    )
    check(
        result.returncode == 0
        and os.path.exists(os.path.join(data, BIG, "PROVENANCE.json")),
        f"{BIG} completed and verified after the resume",
    )
    if result.returncode != 0:
        print(result.stdout[-LOG_TAIL_CHARS:] + result.stderr[-LOG_TAIL_CHARS:])
    # Step 4, status
    status_out = run(env, "status").stdout
    for dataset_id in SMALL + [BIG]:
        line = next((line for line in status_out.splitlines() if line.startswith(dataset_id)), "")
        check("verified" in line and "100%" in line, f"status: {line.strip()!r}")
    print("FETCH CHECK: " + ("PASS" if not bad else f"FAIL ({len(bad)})"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
