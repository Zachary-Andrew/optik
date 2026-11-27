"""FASTA and FASTQ reading and writing, reverse complements and 2-bit k-mer codes, shared by the bench scripts.

Every validator reads its inputs through read_fasta() here, so all of them parse a file the same way.
It is still pure Python and independent of optik's C++ reader, which keeps it usable as an oracle.
"""

import gzip
import os

COMP = str.maketrans("ACGTacgt", "TGCAtgca")

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
# Packed k-mers hold two bits per base, so BASE_MASK - code complements a base.
BITS_PER_BASE = 2
BASE_MASK = (1 << BITS_PER_BASE) - 1
# N50 is the length at which the longest pieces first cover this share of the total.
N50_SHARE = 0.5


def revcomp(seq):
    """Reverse complement of `seq`, lowercase kept lowercase and anything other than ACGT left as is."""
    return seq.translate(COMP)[::-1]


def read_fasta(path):
    """Yield (name, sequence) for each record of a FASTA or FASTQ file, gzipped or not.

    The name is the header's first word, or seq_<n> for a bare '>' or '@', which is how optik names such records too.
    FASTQ may be multi-line: after the '+' line exactly as many quality lines are skipped as there were sequence lines, so a quality line starting with '@' is never taken for a header.
    Text mode turns CRLF into LF, and surrounding whitespace is stripped from every line.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    record_index = 0
    fastq = False
    seq_lines = qual_left = 0
    with opener(path, "rt") as handle:
        for raw in handle:
            if qual_left:
                qual_left -= 1
                continue
            line = raw.strip()
            if not line:
                continue
            if line[0] in ">@":
                if name is not None:
                    yield name, "".join(parts)
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{record_index}"
                record_index += 1
                fastq = line[0] == "@"
                parts, seq_lines = [], 0
            elif line[0] == "+" and fastq and name is not None:
                qual_left = seq_lines
            elif name is not None:
                parts.append(line)
                seq_lines += 1
    if name is not None:
        yield name, "".join(parts)


def write_fasta(path, recs):
    """Write (name, sequence) records as single-line FASTA through a .part file, so an interrupted run never leaves a partial input behind."""
    tmp = f"{path}.part"
    with open(tmp, "w") as handle:
        handle.writelines(f">{name}\n{seq}\n" for name, seq in recs)
    os.replace(tmp, path)


def n50(lengths):
    """Length of the piece at which the longest pieces first cover N50_SHARE of the total, or 0 for no pieces."""
    half = sum(lengths) * N50_SHARE
    covered = 0
    for length in sorted(lengths, reverse=True):
        covered += length
        if covered >= half:
            return length
    return 0


def kmer_codes(seq, kmer_len):
    """Canonical 2-bit codes of every ACGT k-mer of `seq`, where any other character breaks the run."""
    mask = (1 << (BITS_PER_BASE * kmer_len)) - 1
    shift = BITS_PER_BASE * (kmer_len - 1)
    fwd = rev = 0
    valid = 0
    out = []
    for base_char in seq.upper():
        code = CODE.get(base_char)
        if code is None:
            valid = fwd = rev = 0
            continue
        fwd = ((fwd << BITS_PER_BASE) | code) & mask
        rev = (rev >> BITS_PER_BASE) | ((BASE_MASK - code) << shift)
        valid += 1
        if valid >= kmer_len:
            out.append(min(rev, fwd))
    return out


def mutate(seq, rng, sub, ins=0.0, dele=0.0):
    """Copy of `seq` with per-base substitution, insertion and deletion rates `sub`, `ins` and `dele`."""
    out = []
    for base in seq:
        roll = rng.random()
        if roll < dele:
            continue
        if roll < dele + ins:
            out.append(rng.choice("ACGT"))
        out.append(rng.choice("ACGT") if rng.random() < sub else base)
    return "".join(out)
