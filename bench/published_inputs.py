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

from small_common import rc, write_fasta


def _open_text(p):
    """Open p as text, decompressing .xz and .gz by suffix."""
    if p.endswith(".xz"):
        return lzma.open(p, "rt")
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    with open(p) as fh:
        return fh


def _records(path, limit_bases=None):
    """Stream (name, sequence) records of a plain, .gz or .xz FASTA or FASTQ file.

    Reading stops after limit_bases bases for smoke slices, and in that case a truncated compressed prefix is tolerated."""
    total = 0
    name, seq = "", []
    try:
        with _open_text(path) as f:
            first = f.read(1)
            if not first:
                return
            if first == "@":
                # A FASTQ record is a header, the sequence, '+' and the quality line.
                head = first + f.readline()
                while head:
                    sq = f.readline().strip()
                    f.readline()
                    f.readline()
                    yield head[1:].split()[0] if head[1:].split() else "", sq
                    total += len(sq)
                    if limit_bases and total >= limit_bases:
                        return
                    head = f.readline()
                return
            tok = (first + f.readline())[1:].split()
            name, seq, cur = (tok[0] if tok else ""), [], 0
            for line in f:
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
    for f in files:
        m = re.match(r"(.*)_([12])\.fastq\.gz$", os.path.basename(f))
        if not m:
            raise SystemExit(f"{f}: not a _1/_2 mate file")
        pairs.setdefault(m.group(1), {})[m.group(2)] = f
    n = 0
    with open(dst + ".part", "w") as o:
        for run in sorted(pairs):
            a, b = pairs[run].get("1"), pairs[run].get("2")
            if not (a and b):
                raise SystemExit(f"{run}: mate file missing")
            for (na, sa), (nb, sb) in zip(
                _records(a, limit_bases), _records(b, limit_bases), strict=False
            ):
                if na.rsplit("/", 1)[0] != nb.rsplit("/", 1)[0]:
                    raise SystemExit(f"{run}: mates out of step ({na} / {nb})")
                o.write(f">{n}\n{sa}{sb}\n")
                n += 1
    os.replace(dst + ".part", dst)
    return dst


def spss_input(srcs, dst, k=None, limit_bases=None, delete_nonacgt=False):
    """Write one upper-case sequence per line.

    Records are split at every non-ACGT base and pieces shorter than k are dropped, as in scaling_bench.prep and _base.unwrap_fasta, since a k-mer spanning a removed base exists in no sequence.
    With delete_nonacgt those bases are deleted inside each record instead and its pieces joined, as in the SBWT paper's E. coli set.
    Returns (records, bases, k-mer windows at k)."""
    if os.path.exists(dst):
        n = b = w = 0
        with open(dst) as fh:
            for line in fh:
                if not line.startswith(">"):
                    n += 1
                    b += len(line) - 1
                    w += max(0, len(line) - 1 - (k or 1) + 1)
        return n, b, w
    n = b = w = 0
    with open(dst + ".part", "w") as o:
        for src in srcs:
            for _, sq in _records(src, limit_bases):
                pieces = re.findall(r"[ACGT]+", sq.upper())
                for piece in ["".join(pieces)] if delete_nonacgt else pieces:
                    if k is None or len(piece) >= k:
                        o.write(f">{n}\n{piece}\n")
                        n += 1
                        b += len(piece)
                        w += len(piece) - (k or 1) + 1
    os.replace(dst + ".part", dst)
    return n, b, w


def query_segments(src, dst, k, limit_bases=None):
    """Write a query file as ACGT runs of length k or more, one per record.

    Every tool reads the same file (_base.write_acgt_segments), since FMSI answers per base for a read that contains N.
    Returns the number of k-mer windows."""
    if os.path.exists(dst):
        with open(dst) as fh:
            return sum(max(0, len(x) - 1 - k + 1) for x in fh if not x.startswith(">"))
    n = w = 0
    with open(dst + ".part", "w") as o:
        for _, sq in _records(src, limit_bases):
            for piece in re.findall(r"[ACGT]+", sq.upper()):
                if len(piece) >= k:
                    o.write(f">s{n}\n{piece}\n")
                    n += 1
                    w += len(piece) - k + 1
    os.replace(dst + ".part", dst)
    return w


def reference_input(src, dst, limit_bases=None):
    """Write a reference with every header cut to its first token (validate_mapping.sanitize_ref), in upper case with one line per record."""
    if not os.path.exists(dst):
        with open(dst + ".part", "w") as o:
            o.writelines(
                f">{name or f'seq_{i}'}\n{sq.upper()}\n"
                for i, (name, sq) in enumerate(_records(src, limit_bases))
            )
        os.replace(dst + ".part", dst)
    return dst


def fmsi_sshash_m(k, kmers):
    """SSHash m used by the FMSI paper, min(ceil(log4 N) + 1, k - 2) for N k-mers in the set.

    The formula is in fmsi-supplement @ bd31808, experiments/*/getMinimizerLen4SSHash.py ("a formula from Giulio"), and its Snakefiles pass the result as -m."""
    return min(math.ceil(math.log(kmers, 4)) + 1, k - 2)


def sshash_m(dd, k, kmers):
    """SSHash's -m for a dataset at k, since SSHash has no default.

    This is the m the dataset's paper gives, either as a number or as the paper's rule, and k/2 when the paper gives none."""
    if dd.get("sshash_m_rule") == "fmsi":
        return fmsi_sshash_m(k, kmers)
    return dd.get("sshash_m") or k // 2


def sample_kmers(spss, k, n, out, seed=7, rc_half=True):
    """Draw n k-mers uniformly from the k-mer windows of the SPSS, half of them reverse-complemented.

    SSHash 2026 sampled "half of which were reverse complemented".
    Two line-based passes keep any k-mer set out of memory."""
    if os.path.exists(out):
        return out
    with open(spss) as fh:
        lens = [
            max(0, len(line) - 1 - k + 1) for line in fh if not line.startswith(">")
        ]
    total = sum(lens)
    rng = random.Random(seed)
    picks = sorted(rng.randrange(total) for _ in range(min(n, total)))
    with open(spss) as f, open(out + ".part", "w") as o:
        base, j, rec = 0, 0, 0
        for line in f:
            if line.startswith(">"):
                continue
            s = line.strip()
            w = lens[rec]
            while j < len(picks) and picks[j] < base + w:
                p = picks[j] - base
                x = s[p : p + k]
                o.write(f">{j}\n{rc(x) if rc_half and j % 2 else x}\n")
                j += 1
            base += w
            rec += 1
    os.replace(out + ".part", out)
    return out


def windows_in(path, k):
    """Number of k-mer windows in a one-record-per-line query file, which every tool must report having queried."""
    with open(path) as fh:
        return sum(max(0, len(x) - 1 - k + 1) for x in fh if not x.startswith(">"))


def random_kmers(k, n, out, seed=11):
    """Write n uniformly random k-mers to out once (seeded), and return out."""
    if not os.path.exists(out):
        rng = random.Random(seed)
        write_fasta(
            out, ((str(i), "".join(rng.choices("ACGT", k=k))) for i in range(n))
        )
    return out
