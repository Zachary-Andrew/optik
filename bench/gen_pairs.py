#!/usr/bin/env python3
"""gen_pairs.py - simulate paired-end (FR) reads with a truth table.

A fragment whose length follows Normal(mean, sd), and is at least the read length, is drawn uniformly from the N-free reference on either strand.
R1 is its first `len` bases and R2 the reverse complement of its last `len` bases, as in an Illumina FR library.
Each read gets independent substitutions and indels (gen_reads.mutate).
The outputs are R1 and R2 FASTQ (names p<i>/1 and p<i>/2, qualities Q37 falling to Q20 over the last 20 bases) and a TSV with the columns

  name  ref  pos1  strand1  pos2  strand2  fragment

where pos is the 0-based leftmost reference position of the mate's origin and strand is the strand the mate aligns to.

Usage: python3 bench/gen_pairs.py -r ref.fa -n 100000 -o out_prefix
           [-l 150] [--mean 400] [--sd 50] [-e 0.005] [--ins 0.0005] [--del 0.0005] [-s 7]
"""

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_reads import (
    mutate,
    read_fasta,
    revcomp,
)

# Qualities are Phred+33 Q37 ("F"), except that the last QUAL_TAIL_LEN bases fall to Q20 (QUAL_TAIL_MIN_Q) in QUAL_TAIL_DROP steps.
PHRED_OFFSET = 33
HIGH_QUAL_CHAR = "F"
QUAL_TAIL_LEN = 20
QUAL_TAIL_MIN_Q = 20
QUAL_TAIL_DROP = 17
# command-line defaults: pairs, read length, fragment mean and sd, substitution, insertion and deletion rates, seed
DEFAULT_PAIRS = 100000
DEFAULT_READ_LEN = 150
DEFAULT_FRAGMENT_MEAN = 400
DEFAULT_FRAGMENT_SD = 50
DEFAULT_SUBST = 0.005
DEFAULT_INDEL = 0.0005
DEFAULT_SEED = 7
# references must hold a fragment this many sd above the mean
FRAGMENT_SD_SPAN = 4
# probability that a fragment is taken from the forward strand
FORWARD_STRAND_PROB = 0.5


def qual(length):
    """Phred+33 quality string of length n, Q37 ('F') falling to Q20 over the last 20 bases."""
    return "".join(
        HIGH_QUAL_CHAR if j < length - QUAL_TAIL_LEN else chr(PHRED_OFFSET + QUAL_TAIL_MIN_Q + (length - j) * QUAL_TAIL_DROP // QUAL_TAIL_LEN) for j in range(length)
    )


def main():
    """Write the R1, R2 and truth files for the given options."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-r", required=True)
    parser.add_argument("-n", type=int, default=DEFAULT_PAIRS)
    parser.add_argument("-l", type=int, default=DEFAULT_READ_LEN)
    parser.add_argument("--mean", type=int, default=DEFAULT_FRAGMENT_MEAN)
    parser.add_argument("--sd", type=int, default=DEFAULT_FRAGMENT_SD)
    parser.add_argument("-e", type=float, default=DEFAULT_SUBST)
    parser.add_argument("--ins", type=float, default=DEFAULT_INDEL)
    parser.add_argument("--del", dest="dele", type=float, default=DEFAULT_INDEL)
    parser.add_argument("-s", type=int, default=DEFAULT_SEED)
    parser.add_argument("-o", required=True, help="prefix: <o>_1.fq, <o>_2.fq, <o>.tsv")
    args = parser.parse_args()
    rng = random.Random(args.s)
    recs = [(name, seq.upper()) for name, seq in read_fasta(args.r) if len(seq) >= args.mean + FRAGMENT_SD_SPAN * args.sd]
    if not recs:
        sys.exit("no reference record is long enough for the fragment length")
    weights = [len(seq) for _, seq in recs]
    read_len = args.l
    with (
        open(args.o + "_1.fq", "w") as r1_file,
        open(args.o + "_2.fq", "w") as r2_file,
        open(args.o + ".tsv", "w") as truth_file,
    ):
        truth_file.write("name\tref\tpos1\tstrand1\tpos2\tstrand2\tfragment\n")
        i = 0
        while i < args.n:
            name, ref_seq = rng.choices(recs, weights)[0]
            frag_len = max(read_len, round(rng.gauss(args.mean, args.sd)))
            if frag_len > len(ref_seq):
                continue
            frag_start = rng.randrange(0, len(ref_seq) - frag_len + 1)
            frag = ref_seq[frag_start : frag_start + frag_len]
            if "N" in frag:
                continue
            plus = rng.random() < FORWARD_STRAND_PROB
            if not plus:
                frag = revcomp(frag)
            read1 = mutate(frag[:read_len], rng, args.e, args.ins, args.dele)
            read2 = mutate(revcomp(frag)[:read_len], rng, args.e, args.ins, args.dele)
            # R1 lies on the fragment's strand at its start, and R2 on the opposite strand at its end.
            pos1, st1 = (frag_start, "+") if plus else (frag_start + frag_len - read_len, "-")
            pos2, st2 = (frag_start + frag_len - read_len, "-") if plus else (frag_start, "+")
            r1_file.write(f"@p{i}/1\n{read1}\n+\n{qual(len(read1))}\n")
            r2_file.write(f"@p{i}/2\n{read2}\n+\n{qual(len(read2))}\n")
            truth_file.write(f"p{i}\t{name}\t{pos1}\t{st1}\t{pos2}\t{st2}\t{frag_len}\n")
            i += 1


if __name__ == "__main__":
    main()
