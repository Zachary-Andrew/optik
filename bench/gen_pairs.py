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
from gen_reads import mutate, read_fasta, revcomp


def qual(n):
    """Phred+33 quality string of length n, Q37 ('F') falling to Q20 over the last 20 bases."""
    return "".join(
        "F" if j < n - 20 else chr(33 + 20 + (n - j) * 17 // 20) for j in range(n)
    )


def main():
    """Write the R1, R2 and truth files for the given options."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-r", required=True)
    ap.add_argument("-n", type=int, default=100000)
    ap.add_argument("-l", type=int, default=150)
    ap.add_argument("--mean", type=int, default=400)
    ap.add_argument("--sd", type=int, default=50)
    ap.add_argument("-e", type=float, default=0.005)
    ap.add_argument("--ins", type=float, default=0.0005)
    ap.add_argument("--del", dest="dele", type=float, default=0.0005)
    ap.add_argument("-s", type=int, default=7)
    ap.add_argument("-o", required=True, help="prefix: <o>_1.fq, <o>_2.fq, <o>.tsv")
    a = ap.parse_args()
    rng = random.Random(a.s)
    recs = [(n, s.upper()) for n, s in read_fasta(a.r) if len(s) >= a.mean + 4 * a.sd]
    if not recs:
        sys.exit("no reference record is long enough for the fragment length")
    weights = [len(s) for _, s in recs]
    read_len = a.l
    with (
        open(a.o + "_1.fq", "w") as f1,
        open(a.o + "_2.fq", "w") as f2,
        open(a.o + ".tsv", "w") as ft,
    ):
        ft.write("name\tref\tpos1\tstrand1\tpos2\tstrand2\tfragment\n")
        i = 0
        while i < a.n:
            name, g = rng.choices(recs, weights)[0]
            frag_len = max(read_len, round(rng.gauss(a.mean, a.sd)))
            if frag_len > len(g):
                continue
            s = rng.randrange(0, len(g) - frag_len + 1)
            frag = g[s : s + frag_len]
            if "N" in frag:
                continue
            plus = rng.random() < 0.5
            if not plus:
                frag = revcomp(frag)
            r1 = mutate(frag[:read_len], rng, a.e, a.ins, a.dele)
            r2 = mutate(revcomp(frag)[:read_len], rng, a.e, a.ins, a.dele)
            # R1 lies on the fragment's strand at its start, and R2 on the opposite strand at its end.
            p1, st1 = (s, "+") if plus else (s + frag_len - read_len, "-")
            p2, st2 = (s + frag_len - read_len, "-") if plus else (s, "+")
            f1.write(f"@p{i}/1\n{r1}\n+\n{qual(len(r1))}\n")
            f2.write(f"@p{i}/2\n{r2}\n+\n{qual(len(r2))}\n")
            ft.write(f"p{i}\t{name}\t{p1}\t{st1}\t{p2}\t{st2}\t{frag_len}\n")
            i += 1


if __name__ == "__main__":
    main()
