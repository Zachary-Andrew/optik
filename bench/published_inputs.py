"""Input formatting of compare_tools_published_sets.py.

The fetched data are streamed into the formats every compared tool reads.
These are one upper-case sequence per line, ACGT query segments, joined mates, the shared SPSS of eulertigs or of GGCAT or Cuttlefish unitigs, and the sampled query k-mers.
"""

import gzip
import lzma
import math
import os
import random
import re

from dna import revcomp, write_fasta

# nucleotides, the base of log4, and how far below k the FMSI paper's SSHash m stays
NUM_BASES = 4
FMSI_M_BELOW_K = 2
# default seeds of the sampled positive and the random negative k-mers
DEFAULT_SAMPLE_SEED = 7
DEFAULT_RANDOM_SEED = 11


# Without a published m, SSHash runs at m = k / 2.
K_PER_MINIMIZER_LEN = 2
# With rc_half, every other sampled k-mer is reverse-complemented.
RC_EVERY = 2


def _open_text(path):
    """Open `path` as text, decompressing .xz and .gz by suffix.
    The caller closes the returned file, typically with a `with` block."""
    if path.endswith(".xz"):
        return lzma.open(path, "rt")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _records(path, limit_bases=None):
    """Stream (name, sequence) records of a plain, .gz or .xz FASTA or FASTQ file.

    Reading stops after limit_bases bases for smoke slices, and in that case a truncated compressed prefix is tolerated."""
    total = 0
    name, seq = "", []
    try:
        with _open_text(path) as handle:
            first = handle.read(1)
            if not first:
                return
            if first == "@":
                # A FASTQ record is a header, the sequence, '+' and the quality line.
                head = first + handle.readline()
                while head:
                    seq = handle.readline().strip()
                    handle.readline()
                    handle.readline()
                    yield head[1:].split()[0] if head[1:].split() else "", seq
                    total += len(seq)
                    if limit_bases and total >= limit_bases:
                        return
                    head = handle.readline()
                return
            tok = (first + handle.readline())[1:].split()
            name, seq, cur = (tok[0] if tok else ""), [], 0
            for line in handle:
                if line.startswith(">"):
                    yield name, "".join(seq)
                    total += cur
                    tok = line[1:].split()
                    name, seq, cur = (tok[0] if tok else ""), [], 0
                else:
                    seq.append(line.strip())
                    cur += len(seq[-1])
                    if limit_bases and total + cur >= limit_bases:
                        break  # a slice may end inside a long record such as a chromosome
            yield name, "".join(seq)
    except (EOFError, lzma.LZMAError, OSError):
        if not limit_bases:
            raise
        if seq:  # a truncated prefix yields what was read before the cut
            yield name, "".join(seq)


def join_mates(files, dst, limit_bases=None):
    """Join the two mates of each spot into one read, mate 1 then mate 2 as sequenced.

    This is SRA's spot, which fastq-dump writes without --split.
    The files come as <run>_1/<run>_2 pairs read in step, and mismatched mate names are an error."""
    if os.path.exists(dst):
        return dst
    pairs = {}
    for mate_file in files:
        match = re.match(r"(?P<run>.*)_(?P<mate>[12])\.fastq\.gz$", os.path.basename(mate_file))
        if not match:
            raise SystemExit(f"{mate_file}: not a _1/_2 mate file")
        pairs.setdefault(match.group("run"), {})[match.group("mate")] = mate_file
    spot_idx = 0
    with open(dst + ".part", "w") as out_handle:
        for run in sorted(pairs):
            mate1_file, mate2_file = pairs[run].get("1"), pairs[run].get("2")
            if not (mate1_file and mate2_file):
                raise SystemExit(f"{run}: mate file missing")
            for (name1, seq1), (name2, seq2) in zip(
                _records(mate1_file, limit_bases), _records(mate2_file, limit_bases), strict=False
            ):
                if name1.rsplit("/", 1)[0] != name2.rsplit("/", 1)[0]:
                    raise SystemExit(f"{run}: mates out of step ({name1} / {name2})")
                out_handle.write(f">{spot_idx}\n{seq1}{seq2}\n")
                spot_idx += 1
    os.replace(dst + ".part", dst)
    return dst


def spss_input(srcs, dst, kmer_len=None, limit_bases=None, delete_nonacgt=False):
    """Write one upper-case sequence per line.

    Records are split at every non-ACGT base and pieces shorter than k are dropped, as in scaling_bench.prep and _base.unwrap_fasta, since a k-mer spanning a removed base exists in no sequence.
    With delete_nonacgt those bases are deleted inside each record instead and its pieces joined, as in the SBWT paper's E. coli set.
    Returns (records, bases, k-mer windows at k)."""
    if os.path.exists(dst):
        num_records = num_bases = num_windows = 0
        with open(dst) as handle:
            for line in handle:
                if not line.startswith(">"):
                    num_records += 1
                    num_bases += len(line) - 1
                    num_windows += max(0, len(line) - 1 - (kmer_len or 1) + 1)
        return num_records, num_bases, num_windows
    num_records = num_bases = num_windows = 0
    with open(dst + ".part", "w") as out_handle:
        for src in srcs:
            for _, seq in _records(src, limit_bases):
                pieces = re.findall(r"[ACGT]+", seq.upper())
                for piece in ["".join(pieces)] if delete_nonacgt else pieces:
                    if kmer_len is None or len(piece) >= kmer_len:
                        out_handle.write(f">{num_records}\n{piece}\n")
                        num_records += 1
                        num_bases += len(piece)
                        num_windows += len(piece) - (kmer_len or 1) + 1
    os.replace(dst + ".part", dst)
    return num_records, num_bases, num_windows


def query_segments(src, dst, kmer_len, limit_bases=None):
    """Write a query file as ACGT runs of length k or more, one per record.

    Every tool reads the same file (_base.write_acgt_segments), since FMSI answers per base for a read that contains N.
    Returns the number of k-mer windows."""
    if os.path.exists(dst):
        with open(dst) as handle:
            return sum(max(0, len(line) - 1 - kmer_len + 1) for line in handle if not line.startswith(">"))
    num_segments = num_windows = 0
    with open(dst + ".part", "w") as out_handle:
        for _, seq in _records(src, limit_bases):
            for piece in re.findall(r"[ACGT]+", seq.upper()):
                if len(piece) >= kmer_len:
                    out_handle.write(f">s{num_segments}\n{piece}\n")
                    num_segments += 1
                    num_windows += len(piece) - kmer_len + 1
    os.replace(dst + ".part", dst)
    return num_windows


def reference_input(src, dst, limit_bases=None):
    """Write a reference with every header cut to its first token (validate_mapping.sanitize_ref), in upper case with one line per record."""
    if not os.path.exists(dst):
        with open(dst + ".part", "w") as out_handle:
            out_handle.writelines(
                f">{name or f'seq_{i}'}\n{seq.upper()}\n"
                for i, (name, seq) in enumerate(_records(src, limit_bases))
            )
        os.replace(dst + ".part", dst)
    return dst


def fmsi_sshash_m(kmer_len, kmers):
    """SSHash m used by the FMSI paper, min(ceil(log4 N) + 1, k - 2) for N k-mers in the set.

    The formula is in fmsi-supplement @ bd31808, experiments/*/getMinimizerLen4SSHash.py ("a formula from Giulio"), and its Snakefiles pass the result as -m."""
    return min(math.ceil(math.log(kmers, NUM_BASES)) + 1, kmer_len - FMSI_M_BELOW_K)


def sshash_m(dataset, kmer_len, kmers):
    """SSHash's -m for a dataset at k, since SSHash has no default.

    This is the m the dataset's paper gives, either as a number or as the paper's rule, and k/2 when the paper gives none."""
    if dataset.get("sshash_m_rule") == "fmsi":
        return fmsi_sshash_m(kmer_len, kmers)
    return dataset.get("sshash_m") or kmer_len // K_PER_MINIMIZER_LEN


def sample_kmers(spss, kmer_len, num_kmers, out, seed=DEFAULT_SAMPLE_SEED, rc_half=True):
    """Draw n k-mers uniformly from the k-mer windows of the SPSS, half of them reverse-complemented.

    SSHash 2026 sampled "half of which were reverse complemented".
    Two line-based passes keep any k-mer set out of memory."""
    if os.path.exists(out):
        return out
    with open(spss) as handle:
        lens = [
            max(0, len(line) - 1 - kmer_len + 1) for line in handle if not line.startswith(">")
        ]
    total = sum(lens)
    rng = random.Random(seed)
    picks = sorted(rng.randrange(total) for _ in range(min(num_kmers, total)))
    with open(spss) as handle, open(out + ".part", "w") as out_handle:
        base, pick_idx, rec = 0, 0, 0
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip()
            windows = lens[rec]
            while pick_idx < len(picks) and picks[pick_idx] < base + windows:
                offset = picks[pick_idx] - base
                kmer = seq[offset : offset + kmer_len]
                out_handle.write(f">{pick_idx}\n{revcomp(kmer) if rc_half and pick_idx % RC_EVERY else kmer}\n")
                pick_idx += 1
            base += windows
            rec += 1
    os.replace(out + ".part", out)
    return out


def windows_in(path, kmer_len):
    """Number of k-mer windows in a one-record-per-line query file, which every tool must report having queried."""
    with open(path) as handle:
        return sum(max(0, len(line) - 1 - kmer_len + 1) for line in handle if not line.startswith(">"))


def random_kmers(kmer_len, num_kmers, out, seed=DEFAULT_RANDOM_SEED):
    """Write n uniformly random k-mers to out once (seeded), and return out."""
    if not os.path.exists(out):
        rng = random.Random(seed)
        write_fasta(
            out, ((str(i), "".join(rng.choices("ACGT", k=kmer_len))) for i in range(num_kmers))
        )
    return out
