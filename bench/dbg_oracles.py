"""Exact dBG output oracles shared by tool_bench.py and validate_dbg.py."""

from fastq import fastq_sequences


def filter_oracles(path, kmer_len, min_count):
    """Exact expected output k-mer sets computed from the raw input, one per filter model.

    kmer_filter holds the canonical k-mers seen at least c times (optik -c, Bifrost -r/-s).
    edge_filter holds the endpoints of canonical (k+1)-mers seen at least c times (Cuttlefish --read -c, which applies abundance to edges and drops a k-mer without a qualifying edge, such as one only ever seen as an isolated run of length k between Ns)."""
    import gzip
    import re
    from collections import Counter

    tbl = str.maketrans("ACGT", "TGCA")

    def can(kmer):
        """Canonical form of k-mer string x."""
        rev = kmer.translate(tbl)[::-1]
        return min(rev, kmer)

    kmer_counts, edge_counts = Counter(), Counter()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        first = handle.readline()
        fastq = first.startswith("@")
        lines = [first] + list(handle)
    seqs = (
        fastq_sequences(lines)
        if fastq
        else [text_line for text_line in lines if not text_line.startswith(">")]
    )
    if not fastq:  # rejoin wrapped FASTA records
        seqs, cur = [], []
        for text_line in lines:
            if text_line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(text_line.strip())
        if cur:
            seqs.append("".join(cur))
    for line in seqs:
        # Non-ACGT characters split a sequence, so no k-mer or edge spans them.
        for run in re.findall(r"[ACGT]+", line.upper()):
            for j in range(len(run) - kmer_len + 1):
                kmer_counts[can(run[j : j + kmer_len])] += 1
            for j in range(len(run) - kmer_len):
                edge_counts[can(run[j : j + kmer_len + 1])] += 1
    kmer_filter = {kmer for kmer, count in kmer_counts.items() if count >= min_count}
    edge_filter = set()
    # Cuttlefish's model keeps the k-mers at both ends of every solid edge.
    for edge, count in edge_counts.items():
        if count >= min_count:
            edge_filter.add(can(edge[:kmer_len]))
            edge_filter.add(can(edge[1:]))
    return kmer_filter, edge_filter
