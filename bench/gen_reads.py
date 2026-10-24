#!/usr/bin/env python3
"""Generate simulated sequencing reads with ground truth.

Reads are sampled from uniformly random N-free reference windows, taken from the forward strand or, with probability --rc that defaults to 0.5, the reverse strand.
They are mutated with independent per-base substitutions, always to a different base, insertions and deletions.
Read lengths are fixed at -l, or log-normal around -l when --len-sd > 0, as for long reads.

--gt writes a TSV with the columns read_name, ref_name, pos, strand and
read_length.
pos is the 0-based reference start of the sampled window, strand is + or -, and read_length is the reference span of the window (the read itself may be longer or shorter after indels).

Usage
    python gen_reads.py -r ref.fa -n 10000 -l 150 -o reads.fa --gt truth.tsv
    python gen_reads.py -r ref.fa -n 10000 -l 150 -e 0.01 --ins 0.002 --del 0.002 ...
    python gen_reads.py -r ref.fa -n 2000 -l 10000 --len-sd 0.5 -e 0.03 --ins 0.02 --del 0.02 ...

Options
    -r REF      Reference genome in FASTA format
    -n NUM      Number of reads to generate [default 10000]
    -l LEN      Read length (the mean when --len-sd > 0) [default 150]
    -e ERR      Per-base substitution rate [default 0.0]
    --ins R     Per-base insertion rate [default 0.0]
    --del R     Per-base deletion rate [default 0.0]
    --rc F      Fraction of reads from the reverse strand [default 0.5]
    --len-sd S  Log-normal sigma of the read length (0 = fixed) [default 0]
    -o OUT      Output file in FASTA format [default stdout]
    -s SEED     Random seed [default 42]
    -c COV      Target coverage, instead of -n [default 0]
    --gt GT     Ground truth TSV output
"""

import argparse
import contextlib
import random
import sys


def read_fasta(path):
    """Read a FASTA file into a list of (name, sequence) tuples.

    Sequence lines are joined up to the next header, and only the first whitespace-delimited token of a header is kept as the name.
    Empty headers (a bare '>') get sequential names seq_0, seq_1 and so on.
    """
    sequences = []
    name = None
    seq_parts = []
    seq_idx = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # A new record header, so the previous record is saved first.
                if name is not None:
                    sequences.append((name, "".join(seq_parts)))
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                seq_parts = []
            else:
                seq_parts.append(line)
    # The last record
    if name is not None:
        sequences.append((name, "".join(seq_parts)))
    return sequences


def lemire_randbelow(rng, n):
    """Unbiased random integer in [0, n) by Lemire's method (Lemire 2019).

    ``x % n`` carries modulo bias whenever n does not divide the PRNG range.
    Lemire's scheme maps the 64-bit draw into [0, n) uniformly and redraws only on a rare threshold collision.
    A power-of-two n reduces to a masked draw.
    """
    if n <= 0:
        raise ValueError("n must be >= 1")
    if n & (n - 1) == 0:
        return rng.getrandbits(n.bit_length() - 1)
    t = (-n) % n
    mask64 = (1 << 64) - 1
    while True:
        x = rng.getrandbits(64)
        m = x * n
        if (m & mask64) >= t:
            return m >> 64


COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    """Reverse complement of seq."""
    return seq.translate(COMP)[::-1]


def mutate(read, rng, sub, ins, dele):
    """Apply independent per-base substitutions, insertions and deletions."""
    if sub <= 0 and ins <= 0 and dele <= 0:
        return read
    out = []
    for c in read:
        x = rng.random()
        if x < dele:
            continue
        if x < dele + sub:
            alt = "ACGT".replace(c.upper(), "")
            out.append(
                alt[lemire_randbelow(rng, len(alt))]
                if len(alt) == 3
                else "ACGT"[lemire_randbelow(rng, 4)]
            )
        else:
            out.append(c)
        if rng.random() < ins:
            out.append("ACGT"[lemire_randbelow(rng, 4)])
    return "".join(out)


def generate_reads(
    sequences,
    n,
    read_len,
    error_rate,
    seed,
    ins_rate=0.0,
    del_rate=0.0,
    rc_frac=0.5,
    len_sd=0.0,
):
    """Generate reads and return (reads, ground_truth).

    Ground-truth rows are (name, ref_name, pos, strand, ref_span).
    """
    import math

    rng = random.Random(seed)
    sequences = [s for s in sequences if len(s[1]) >= read_len]
    if not sequences:
        raise ValueError("no reference sequence is long enough for the read length")

    reads = []
    ground_truth = []
    for i in range(n):
        span = read_len
        if len_sd > 0:
            span = max(
                50,
                round(read_len * math.exp(rng.gauss(-len_sd * len_sd / 2, len_sd))),
            )
        read = None
        seq_name = seq = pos = None
        for _ in range(200):
            seq_name, seq = sequences[lemire_randbelow(rng, len(sequences))]
            if len(seq) < span:
                continue
            pos = lemire_randbelow(rng, len(seq) - span + 1)
            candidate = seq[pos : pos + span]
            if "N" not in candidate and "n" not in candidate:
                read = candidate
                break
        if read is None:
            continue
        strand = "-" if rng.random() < rc_frac else "+"
        read = mutate(read, rng, error_rate, ins_rate, del_rate)
        if strand == "-":
            read = revcomp(read)
        reads.append((f">{i}", read))
        ground_truth.append((f"{i}", seq_name, pos, strand, span))
    return reads, ground_truth


def main():
    """Generate reads per the options and write them along with their truth."""
    parser = argparse.ArgumentParser(
        description="Generate random reads from a reference"
    )
    parser.add_argument("-r", required=True, help="Reference FASTA")
    parser.add_argument("-n", type=int, default=10000, help="Number of reads")
    parser.add_argument("-l", type=int, default=150, help="Read length")
    parser.add_argument("-e", type=float, default=0.0, help="Substitution rate")
    parser.add_argument("--ins", type=float, default=0.0, help="Insertion rate")
    parser.add_argument(
        "--del", dest="dele", type=float, default=0.0, help="Deletion rate"
    )
    parser.add_argument("--rc", type=float, default=0.5, help="Reverse-strand fraction")
    parser.add_argument(
        "--len-sd", type=float, default=0.0, help="Log-normal length sigma"
    )
    parser.add_argument("-o", default=None, help="Output FASTA")
    parser.add_argument("-s", type=int, default=42, help="Random seed")
    parser.add_argument("-c", type=float, default=0.0, help="Target coverage")
    parser.add_argument("--gt", default=None, help="Ground truth TSV output")
    args = parser.parse_args()

    sequences = read_fasta(args.r)
    total_len = sum(len(s) for _, s in sequences)

    # With a coverage given, the read count follows from coverage = n * read_len / total_len, so n = coverage * total_len / read_len.
    if args.c > 0:
        args.n = max(1, int(args.c * total_len / args.l))

    reads, ground_truth = generate_reads(
        sequences,
        args.n,
        args.l,
        args.e,
        args.s,
        args.ins,
        args.dele,
        args.rc,
        args.len_sd,
    )

    with contextlib.ExitStack() as stack:
        out = stack.enter_context(open(args.o, "w")) if args.o else sys.stdout
        for name, seq in reads:
            out.write(f"{name}\n{seq}\n")

    if args.gt:
        gt_path = args.gt
        with open(gt_path, "w") as f:
            f.write("read_name\tref_name\tpos\tstrand\tread_length\n")
            for name, ref_name, pos, strand, length in ground_truth:
                f.write(f"{name}\t{ref_name}\t{pos}\t{strand}\t{length}\n")


if __name__ == "__main__":
    main()
