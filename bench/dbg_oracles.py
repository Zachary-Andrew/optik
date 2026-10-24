"""Exact dBG output oracles shared by tool_bench.py and validate_dbg.py."""


def filter_oracles(path, k, c):
    """Exact expected output k-mer sets computed from the raw input, one per filter model.

    kmer_filter holds the canonical k-mers seen at least c times (optik -c, Bifrost -r/-s).
    edge_filter holds the endpoints of canonical (k+1)-mers seen at least c times (Cuttlefish --read -c, which applies abundance to edges and drops a k-mer without a qualifying edge, such as one only ever seen as an isolated run of length k between Ns)."""
    import gzip
    import re
    from collections import Counter

    tbl = str.maketrans("ACGT", "TGCA")

    def can(x):
        """Canonical form of k-mer string x."""
        r = x.translate(tbl)[::-1]
        return min(r, x)

    kc, ec = Counter(), Counter()
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        first = fh.readline()
        fastq = first.startswith("@")
        lines = [first] + list(fh)
    seqs = (
        lines[1::4]
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
        for run in re.findall(r"[ACGT]+", line.upper()):
            for j in range(len(run) - k + 1):
                kc[can(run[j : j + k])] += 1
            for j in range(len(run) - k):
                ec[can(run[j : j + k + 1])] += 1
    kf = {x for x, n in kc.items() if n >= c}
    ef = set()
    for e, n in ec.items():
        if n >= c:
            ef.add(can(e[:k]))
            ef.add(can(e[1:]))
    return kf, ef
