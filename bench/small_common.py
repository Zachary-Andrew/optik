"""Paths, tool binaries and sequence helpers shared by compare_tools_small_sets.py, external_oracles.py and compare_tools_published_sets.py."""

import gzip
import os

from config import PATHS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


T = os.path.join(ROOT, "temp")

BIN = {
    "optik": os.path.join(ROOT, "target/release/optik"),
    "sshash": str(PATHS["sshash"]),
    "sbwt": str(PATHS["sbwt"]),
    "fmsi": str(PATHS["fmsi"]),
    "kmercamel": str(PATHS["kmercamel"]),
    "cuttlefish": str(PATHS["cuttlefish"]),
    "bifrost": str(PATHS["bifrost"]),
    "minimap2": str(PATHS["minimap2"]),
    "wfmash": str(PATHS["wfmash"]),
    "kmc": os.path.join(T, "KMC/bin/kmc"),
    "kmc_dump": os.path.join(T, "KMC/bin/kmc_dump"),
    "kmc_tools": os.path.join(T, "KMC/bin/kmc_tools"),
}

DATA = os.path.join(T, "sshash/data")

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def rc(s):
    """Reverse complement of DNA string s."""
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """Iterate over the records of a FASTA or FASTQ file (gzip allowed) as (name, sequence)."""
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith((">", "@")):
                if name is not None:
                    yield name, "".join(seq)
                tok = line[1:].split()
                name, seq = (tok[0] if tok else ""), []
                if line.startswith("@"):
                    s = next(f).strip()
                    next(f)
                    next(f)
                    yield name, s
                    name = None
            elif name is not None:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path, recs):
    """Write records atomically, so an interrupted run never leaves a partial input for the next."""
    tmp = path + ".part"
    with open(tmp, "w") as f:
        f.writelines(f">{n}\n{s}\n" for n, s in recs)
    os.replace(tmp, path)


def kmer_codes(s, k):
    """Canonical 2-bit codes of every ACGT k-mer of s, where N breaks the run."""
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    f = r = 0
    valid = 0
    out = []
    for ch in s.upper():
        c = CODE.get(ch)
        if c is None:
            valid = f = r = 0
            continue
        f = ((f << 2) | c) & mask
        r = (r >> 2) | ((3 - c) << shift)
        valid += 1
        if valid >= k:
            out.append(min(r, f))
    return out


def mutate(s, rng, sub, ins=0.0, dele=0.0):
    """Copy of s with per-base substitution, insertion and deletion rates sub, ins and dele."""
    out = []
    for c in s:
        x = rng.random()
        if x < dele:
            continue
        if x < dele + ins:
            out.append(rng.choice("ACGT"))
        out.append(rng.choice("ACGT") if rng.random() < sub else c)
    return "".join(out)
