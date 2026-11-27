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

from dna import read_fasta, revcomp

# Lemire draws are 64-bit, and a substitution picks one of the other NUM_BASES - 1 bases
RNG_WORD_BITS = 64
NUM_BASES = 4
# fraction of reads taken from the reverse strand by default
DEFAULT_RC_FRAC = 0.5
# log-normal read lengths are at least this long
MIN_READ_SPAN = 50
# N-free placements tried per read before giving up
MAX_PLACEMENT_TRIES = 200
# command-line defaults
DEFAULT_NUM_READS = 10000
DEFAULT_READ_LEN = 150
DEFAULT_SEED = 42


# A log-normal factor with mu = -sigma^2 / 2 has mean 1, so read lengths average the requested length.
LOGNORMAL_MEAN_DEN = 2


def lemire_randbelow(rng, bound):
    """Unbiased random integer in [0, n) by Lemire's method (Lemire 2019).

    ``x % n`` carries modulo bias whenever n does not divide the PRNG range.
    Lemire's scheme maps the 64-bit draw into [0, n) uniformly and redraws only on a rare threshold collision.
    A power-of-two n reduces to a masked draw.
    """
    if bound <= 0:
        raise ValueError("n must be >= 1")
    if bound & (bound - 1) == 0:
        return rng.getrandbits(bound.bit_length() - 1)
    threshold = (-bound) % bound
    mask64 = (1 << RNG_WORD_BITS) - 1
    while True:
        draw = rng.getrandbits(RNG_WORD_BITS)
        product = draw * bound
        if (product & mask64) >= threshold:
            return product >> RNG_WORD_BITS


def mutate(read, rng, sub, ins, dele):
    """Apply independent per-base substitutions, insertions and deletions."""
    if sub <= 0 and ins <= 0 and dele <= 0:
        return read
    out = []
    for base in read:
        # One roll picks deletion, substitution or keep. An insertion is drawn separately after the base.
        roll = rng.random()
        if roll < dele:
            continue
        if roll < dele + sub:
            alt = "ACGT".replace(base.upper(), "")
            out.append(
                alt[lemire_randbelow(rng, len(alt))]
                # An ambiguous base has no defined alternative, so any base may replace it.
                if len(alt) == NUM_BASES - 1
                else "ACGT"[lemire_randbelow(rng, NUM_BASES)]
            )
        else:
            out.append(base)
        if rng.random() < ins:
            out.append("ACGT"[lemire_randbelow(rng, NUM_BASES)])
    return "".join(out)


def generate_reads(
    sequences,
    num_reads,
    read_len,
    error_rate,
    seed,
    ins_rate=0.0,
    del_rate=0.0,
    rc_frac=DEFAULT_RC_FRAC,
    len_sd=0.0,
):
    """Generate reads and return (reads, ground_truth).

    Ground-truth rows are (name, ref_name, pos, strand, ref_span).
    """
    import math

    rng = random.Random(seed)
    sequences = [entry for entry in sequences if len(entry[1]) >= read_len]
    if not sequences:
        raise ValueError("no reference sequence is long enough for the read length")

    reads = []
    ground_truth = []
    for i in range(num_reads):
        span = read_len
        if len_sd > 0:
            span = max(
                MIN_READ_SPAN,
                round(read_len * math.exp(rng.gauss(-len_sd * len_sd / LOGNORMAL_MEAN_DEN, len_sd))),
            )
        read = None
        seq_name = seq = pos = None
        for _ in range(MAX_PLACEMENT_TRIES):
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
    parser.add_argument("-n", type=int, default=DEFAULT_NUM_READS, help="Number of reads")
    parser.add_argument("-l", type=int, default=DEFAULT_READ_LEN, help="Read length")
    parser.add_argument("-e", type=float, default=0.0, help="Substitution rate")
    parser.add_argument("--ins", type=float, default=0.0, help="Insertion rate")
    parser.add_argument(
        "--del", dest="dele", type=float, default=0.0, help="Deletion rate"
    )
    parser.add_argument("--rc", type=float, default=DEFAULT_RC_FRAC, help="Reverse-strand fraction")
    parser.add_argument(
        "--len-sd", type=float, default=0.0, help="Log-normal length sigma"
    )
    parser.add_argument("-o", default=None, help="Output FASTA")
    parser.add_argument("-s", type=int, default=DEFAULT_SEED, help="Random seed")
    parser.add_argument("-c", type=float, default=0.0, help="Target coverage")
    parser.add_argument("--gt", default=None, help="Ground truth TSV output")
    args = parser.parse_args()

    sequences = list(read_fasta(args.r))
    total_len = sum(len(seq) for _, seq in sequences)

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
        with open(gt_path, "w") as handle:
            handle.write("read_name\tref_name\tpos\tstrand\tread_length\n")
            for name, ref_name, pos, strand, length in ground_truth:
                handle.write(f"{name}\t{ref_name}\t{pos}\t{strand}\t{length}\n")


if __name__ == "__main__":
    main()
