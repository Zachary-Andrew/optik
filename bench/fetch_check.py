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

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.join(HERE, "compare_tools_published_sets.py")
SMALL = ["cf2-celegans-ref", "fmsi25-spneumo"]
BIG = "sshash26-ncbi-virus"


def run(env, *args, **kw):
    return subprocess.run(
        [sys.executable, "-u", BENCH, *args],
        env=env,
        capture_output=True,
        text=True,
        **kw,
        check=False,
    )


def done_bytes(data, i):
    """Bytes of dataset i in the completed ranges its .part.done files record."""
    d = os.path.join(data, i)
    n = 0
    for f in os.listdir(d) if os.path.isdir(d) else []:
        if f.endswith(".part.done"):
            n += sum(e - s for s, e in json.loads(Path(os.path.join(d, f)).read_text()))
    return n


def main():
    """Run the four fetch steps against a scratch data root."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", default=os.path.join(HERE, "..", "temp", "fetch_check"))
    a = ap.parse_args()
    data = os.path.abspath(a.work)
    shutil.rmtree(data, ignore_errors=True)
    os.makedirs(data)
    env = dict(os.environ, OPTIK_BENCH_DATA=data)
    bad = []

    def check(ok, what):
        """Print one step result and remember a failure."""
        print(f"  [{'PASS' if ok else 'FAIL'}] {what}", flush=True)
        if not ok:
            bad.append(what)

    print(f"scratch data root: {data}")
    # Step 1, small datasets
    r = run(env, "fetch", "--only", ",".join(SMALL), "--force")
    for i in SMALL:
        check(
            os.path.exists(os.path.join(data, i, "PROVENANCE.json")),
            f"{i} fetched and verified",
        )
    if r.returncode != 0:
        print(r.stdout[-1500:] + r.stderr[-1500:])
    # Step 2, a fetch stopped partway
    with open(os.path.join(data, "big1.log"), "w") as log:
        p = subprocess.Popen(
            [sys.executable, "-u", BENCH, "fetch", "--only", BIG, "--force"],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        t0 = time.time()
        # The scheduler records completed ranges every 5 s, so wait for some but stop well before the end.
        while (
            time.time() - t0 < 600 and p.poll() is None and done_bytes(data, BIG) < 20e6
        ):
            time.sleep(1)
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()
    before = done_bytes(data, BIG)
    check(
        before > 0 and not os.path.exists(os.path.join(data, BIG, "PROVENANCE.json")),
        f"{BIG} stopped partway ({before} bytes recorded on disk)",
    )
    st = run(env, "status").stdout
    line = next((x for x in st.splitlines() if x.startswith(BIG)), "")
    check("partial" in line, f"status reports the stopped fetch: {line.strip()!r}")
    # Step 3, resume
    r = run(env, "fetch", "--only", BIG, "--force")
    m = re.findall(r"resuming \S+: (\d+) of (\d+) bytes already on disk", r.stdout)
    resumed = sum(int(x) for x, _ in m)
    check(
        resumed == before and before > 0,
        f"resumed from the {before} bytes on disk (scheduler: {resumed})",
    )
    check(
        r.returncode == 0
        and os.path.exists(os.path.join(data, BIG, "PROVENANCE.json")),
        f"{BIG} completed and verified after the resume",
    )
    if r.returncode != 0:
        print(r.stdout[-1500:] + r.stderr[-1500:])
    # Step 4, status
    st = run(env, "status").stdout
    for i in SMALL + [BIG]:
        line = next((x for x in st.splitlines() if x.startswith(i)), "")
        check("verified" in line and "100%" in line, f"status: {line.strip()!r}")
    print("FETCH CHECK: " + ("PASS" if not bad else f"FAIL ({len(bad)})"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
