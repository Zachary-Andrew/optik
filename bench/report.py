#!/usr/bin/env python3
"""report.py - one sectioned Markdown file that holds the results of every bench.

Each bench writes its results as one section of bench/results/RESULTS.md, or of the file given by its --report flag.
A section looks like this:

  <!-- section: scaling -->
  ## scaling
  Size scaling of every tool on E. coli, chr21 and chr1

  Run 2026-10-23 14:02 UTC, git 1f4e417-dirty, optik sha256 3f9a0c..., 16 CPUs, 31 GiB,
  `python3 bench/scaling_bench.py temp/scaling --runs 3`

  ### ecoli
  | tool | phase | wall_s | rss_MB | status |
  |---|---|---|---|---|
  | optik | build k=31 | 1.52 | 210 | ok |
  <!-- end: scaling -->

Running a bench again replaces only its own section, so a campaign of benches fills one file that can be read, diffed and committed.
The tables are plain pipe tables, and read_report() parses the file back into rows.

  python3 bench/report.py [FILE]               lists the sections and their tables
  python3 bench/report.py [FILE] --json        prints every section as JSON
"""

import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_PATH = os.path.join(ROOT, "bench", "results", "RESULTS.md")
HEADER = """# optik benchmark results

Each section below comes from the bench script named in its heading, and running that bench again replaces it in place.
A section's provenance line records the date, commit, optik binary and machine of its run.
`python3 bench/run_benchmarks.py list` describes every bench, and `python3 bench/report.py --json` reads the tables back as data.
"""
OPEN = "<!-- section: {} -->"
CLOSE = "<!-- end: {} -->"
SECTION_RE = re.compile(
    r"^<!-- section: (\S+) -->\n(.*?)^<!-- end: \1 -->\n?", re.MULTILINE | re.DOTALL
)


def add_argument(ap):
    """Adds --report PATH to a bench's argument parser."""
    ap.add_argument(
        "--report",
        default=DEFAULT_PATH,
        help="sectioned results file this bench writes its section into "
        f"(default: {os.path.relpath(DEFAULT_PATH, ROOT)}; '' to skip)",
    )


def cell(v):
    """Text of one table cell.

    A measurement dict prints as its mean with the sd, a float with four significant digits, a bool as PASS or FAIL, and None as empty."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "PASS" if v else "FAIL"
    if isinstance(v, dict) and "mean" in v:
        sd = v.get("sd")
        return f"{v['mean']:.4g}" + (f" ± {sd:.2g}" if sd is not None else "")
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v).replace("|", "\\|").replace("\n", " ")


def table(rows, cols=None, drop=()):
    """Markdown lines of a pipe table of row dicts.

    Columns are `cols` when given, else every key in first-seen order less those in `drop`."""
    if not rows:
        return []
    if cols is None:
        cols = []
        for r in rows:
            cols += [k for k in r if k not in cols and k not in drop]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(cell(r.get(k)) for k in cols) + " |" for r in rows]
    return lines


def _git_rev():
    """Short commit of the checkout, with -dirty when tracked files differ from it."""
    try:
        return subprocess.run(
            ["git", "-C", ROOT, "describe", "--always", "--dirty"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
    except OSError:
        return ""


def _mem_gib():
    """MemTotal in GiB, or None when /proc/meminfo is unreadable."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) >> 20
    except OSError:
        pass
    return None


def provenance(optik=None, extra=""):
    """One line naming when, from which commit and binary, and on what machine a run was made."""
    when = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    parts = [f"Run {when}", f"git {_git_rev() or 'unknown'}"]
    optik = optik or os.path.join(ROOT, "target", "release", "optik")
    if os.path.exists(optik):
        with open(optik, "rb") as f:
            parts.append(f"optik sha256 {hashlib.sha256(f.read()).hexdigest()[:12]}")
    parts.append(f"{len(os.sched_getaffinity(0))} CPUs")
    mem = _mem_gib()
    if mem is not None:
        parts.append(f"{mem} GiB")
    parts.append(platform.machine())
    if extra:
        parts.append(extra)
    argv = " ".join(
        os.path.relpath(x, ROOT) if x.startswith(ROOT + os.sep) else x for x in sys.argv
    )
    return ", ".join(parts) + f",\n`python3 {argv}`"


def write_section(path, key, title, tables, notes=(), optik=None, extra=""):
    """Writes one bench's section into the results file at `path`.

    `tables` is a list of (subtitle, rows) pairs, or of (subtitle, rows, cols) to fix the column order, and `notes` are Markdown lines placed after the tables.
    An existing section with the same key is replaced in place, and a new one is appended.
    An empty path writes nothing."""
    if not path:
        return
    body = [OPEN.format(key), f"## {key}", title, "", provenance(optik, extra), ""]
    for t in tables:
        sub, rows = t[0], t[1]
        cols = t[2] if len(t) > 2 else None
        if not rows:
            continue
        if sub:
            body += [f"### {sub}", ""]
        body += table(rows, cols) + [""]
    for i, note in enumerate(notes):
        # A plain line after a list item would join that item, so a blank line ends the list.
        if i and notes[i - 1].startswith("- ") and not note.startswith("- "):
            body.append("")
        body.append(note)
    if notes:
        body.append("")
    body.append(CLOSE.format(key))
    text = "\n".join(body) + "\n"
    old = Path(path).read_text() if os.path.exists(path) else HEADER
    m = re.search(
        rf"^{re.escape(OPEN.format(key))}\n.*?^{re.escape(CLOSE.format(key))}\n?",
        old,
        re.MULTILINE | re.DOTALL,
    )
    if m:
        new = old[: m.start()] + text + old[m.end() :]
    else:
        new = old.rstrip("\n") + "\n\n" + text
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "w") as f:
        f.write(new)
    os.replace(tmp, path)
    print(f"report: section '{key}' written to {path}", file=sys.stderr)


def _split_row(line):
    """Cells of one pipe-table row, with escaped pipes restored."""
    cells = re.split(r"(?<!\\)\|", line.strip())[1:-1]
    return [c.strip().replace("\\|", "|") for c in cells]


def read_report(path=DEFAULT_PATH):
    """Parses a results file back into {section: {title, provenance, tables}}.

    `tables` maps each table's subtitle, '' for an untitled one, to its rows as dicts of column name to cell text."""
    out = {}
    for m in SECTION_RE.finditer(Path(path).read_text()):
        key, lines = m.group(1), m.group(2).splitlines()
        sec = {"title": "", "provenance": "", "tables": {}}
        sub, i = "", 0
        head = [x for x in lines[1:] if x.strip()][:3]
        sec["title"] = head[0] if head else ""
        sec["provenance"] = " ".join(
            x for x in head[1:3] if not x.startswith(("|", "###"))
        )
        while i < len(lines):
            line = lines[i]
            if line.startswith("### "):
                sub = line[4:].strip()
            elif (
                line.startswith("|")
                and i + 1 < len(lines)
                and lines[i + 1].startswith("|---")
            ):
                cols = _split_row(line)
                rows, i = [], i + 2
                while i < len(lines) and lines[i].startswith("|"):
                    rows.append(dict(zip(cols, _split_row(lines[i]), strict=False)))
                    i += 1
                name, n = sub, 2
                while name in sec["tables"]:
                    name, n = f"{sub} ({n})", n + 1
                sec["tables"][name] = rows
                continue
            i += 1
        out[key] = sec
    return out


def main():
    """Lists the sections of a results file, or prints them as JSON."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("file", nargs="?", default=DEFAULT_PATH)
    ap.add_argument("--json", action="store_true", help="print every section as JSON")
    a = ap.parse_args()
    if not os.path.exists(a.file):
        sys.exit(
            f"{a.file}: no results yet (run a bench, see bench/run_benchmarks.py list)"
        )
    rep = read_report(a.file)
    if a.json:
        json.dump(rep, sys.stdout, indent=1)
        print()
        return
    for key, sec in rep.items():
        print(f"{key}: {sec['title']}")
        for sub, rows in sec["tables"].items():
            print(f"  {sub or '(table)'}: {len(rows)} rows")


if __name__ == "__main__":
    main()
