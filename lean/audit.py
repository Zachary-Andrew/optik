#!/usr/bin/env python3
"""audit.py - checks that every theorem in lean/OptikProofs is proved, and proved from the standard axioms only.

It fails on any `sorry`, `admit` or project `axiom` in the sources.
It then has Lean print the axioms behind every theorem.
Allowed are propext, Classical.choice and Quot.sound, plus the `_native.bv_decide` axioms registered by the theorems in NATIVE_SOURCES.
Those record that Lean's compiled checker verified each SAT certificate, because the kernel cannot evaluate 64-bit words directly.
A theorem that uses one of those theorems inherits its native axioms, which is allowed, but a native axiom from any other theorem is not.

Run it from lean/ as python3 audit.py, with lake on PATH from elan or pixi.
"""

import pathlib
import re
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
STANDARD = {"propext", "Classical.choice", "Quot.sound"}
# The theorems that call bv_decide.
# A new bv_decide proof fails the audit until its theorem is added here.
NATIVE_SOURCES = {
    "OptikProofs.Kmer.rcWord_eq_rcSpec",
    "OptikProofs.Hash.xs33_inj",
    "OptikProofs.Hash.xs24_inj",
    "OptikProofs.Hash.xs14_inj",
    "OptikProofs.Hash.xs28_inj",
    "OptikProofs.Rolling.or_low_eq_add",
    "OptikProofs.Varint.last_byte_payload",
    "OptikProofs.Varint.last_byte_stops",
    "OptikProofs.Varint.more_byte_payload",
    "OptikProofs.Varint.more_byte_continues",
    "OptikProofs.Varint.split_low7",
    "OptikProofs.Varint.setWidth_roundtrip",
    "OptikProofs.Bpm.carry",
}
NATIVE_MARK = "._native.bv_decide."


def theorem_names():
    """Fully qualified names of every theorem declared in OptikProofs."""
    names = []
    for path in sorted((HERE / "OptikProofs").glob("*.lean")):
        names += [f"OptikProofs.{path.stem}.{name}" for name in re.findall(r"(?m)^theorem ([\w'.]+)", path.read_text())]
    return names


def main():
    """Print each problem found and return 1 if there was any."""
    problems = []
    for path in sorted((HERE / "OptikProofs").glob("*.lean")):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"\bsorry\b|\badmit\b|^axiom\b", line):
                problems.append(f"{path.name}:{number}: {line.strip()}")
    names = theorem_names()
    with tempfile.NamedTemporaryFile("w", suffix=".lean", dir=HERE, delete=False) as handle:
        handle.write("import OptikProofs\n" + "".join(f"#print axioms {name}\n" for name in names))
        script = pathlib.Path(handle.name)
    try:
        out = subprocess.run(["lake", "env", "lean", str(script)], cwd=HERE, capture_output=True, text=True, check=False)
    finally:
        script.unlink()
    report = out.stdout + out.stderr
    # Each answer names the theorem and lists its axioms in brackets, possibly over several lines.
    answers = {m.group(1): m.group(2) for m in re.finditer(r"'([\w'.]+)' depends on axioms: \[([^\]]*)\]", report, re.S)}
    answers.update({m.group(1): "" for m in re.finditer(r"'([\w'.]+)' does not depend on any axioms", report)})
    for name in names:
        if name not in answers:
            problems.append(f"{name}: no axiom report from Lean")
            continue
        axioms = {a.strip() for a in answers[name].split(",") if a.strip()}
        extra = {a for a in axioms - STANDARD if not (NATIVE_MARK in a and a.split(NATIVE_MARK)[0] in NATIVE_SOURCES)}
        if extra:
            problems.append(f"{name}: unexpected axioms {sorted(extra)}")
    for problem in problems:
        print(f"  [FAIL] {problem}")
    print(f"lean audit: {len(names)} theorems, {'ALL PROVED' if not problems else f'{len(problems)} PROBLEM(S)'}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
