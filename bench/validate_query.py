#!/usr/bin/env python3
"""Validate optik query against a set-based k-mer oracle and SSHash.

Checks
  1. Oracle match.
     Every k-mer of the reference is found, and random k-mers absent from the reference are rejected.
  2. Positive count agreement.
     optik's query summary equals the count the oracle expects.
  3. Near-miss rejection.
     k-mers one substitution away from the reference that are not in the index are rejected.
  4. Streaming records.
     Reference windows from both strands with a few substitutions mix present and absent k-mers in each record, and their counts must match the oracle k-mer by k-mer.
  5. Cross-tool agreement.
     When SSHash is available, its reported counts on every set are compared as well.

Both optik index formats (build-spss and build-fm) are checked on every set.

Usage
    python validate_query.py -r ref.fa -k 31 [--sshash /path/to/sshash]

A JSON report goes to stdout and a human-readable summary to stderr.
"""

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

# k-mer helpers

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    """Reverse complement of seq."""
    return seq.translate(COMPLEMENT)[::-1]


def canonicalize(kmer: str) -> str:
    """Canonical form of kmer, the smaller of it and its reverse complement."""
    r = rc(kmer)
    return min(kmer, r)


def extract_kmers(seq: str, k: int):
    """Canonical k-mers of seq, case-insensitive, where any non-ACGT base breaks the run as it does in optik and SSHash."""
    kmers = []
    run = []
    for c in seq.upper() + "N":
        if c in "ACGT":
            run.append(c)
            continue
        r = "".join(run)
        kmers.extend(canonicalize(r[j : j + k]) for j in range(len(r) - k + 1))
        run = []
    return kmers


def read_fasta(path: str):
    """Read a FASTA file into a list of (name, sequence) records."""
    sequences = []
    name = None
    parts = []
    seq_idx = 0
    import gzip

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(parts)))
                tokens = line[1:].split()
                name = tokens[0] if tokens else f"seq_{seq_idx}"
                seq_idx += 1
                parts = []
            else:
                parts.append(line)
    if name is not None:
        sequences.append((name, "".join(parts)))
    return sequences


def build_oracle(ref_path: str, k: int):
    """Set of all canonical k-mers of the reference."""
    oracle = set()
    for _, seq in read_fasta(ref_path):
        oracle.update(extract_kmers(seq, k))
    return oracle


def random_dna(length: int, rng: random.Random) -> str:
    """Random ACGT string of the given length."""
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_single_kmer_fasta(kmers, path):
    """Write a FASTA file with one k-mer per record."""
    with open(path, "w") as f:
        f.writelines(f">q{i}\n{kmer}\n" for i, kmer in enumerate(kmers))


def run_optik_query(index_path: str, queries_fa: str, k: int, optik_bin: str):
    """Run optik query and parse its summary output."""
    cmd = [optik_bin, "query", "-i", index_path, "-q", queries_fa, "-k", str(k)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300, check=False)
    if proc.returncode != 0:
        return None, f"optik query failed: {proc.stderr[:500]}"

    result = {}
    for line in proc.stdout.split("\n"):
        line = line.strip()
        if "=" in line:
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.strip()
            # Drop the percentage annotation.
            if "(" in val:
                val = val[: val.index("(")].strip()
            try:
                result[key] = int(val)
            except ValueError:
                try:
                    result[key] = float(val)
                except ValueError:
                    result[key] = val
    return result, None


def parse_counts(text: str):
    """Parse 'num_kmers = N' and 'num_positive_kmers = N (x%)' lines into a dict."""
    out = {}
    for line in text.split("\n"):
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        val = val.strip().split("(")[0].strip()
        try:
            out[key.strip()] = int(val)
        except ValueError:
            pass
    return out


def validate(args):
    """Run every query check against the oracle and SSHash and return the report dict."""
    checks = {}
    k = args.k
    optik = args.optik
    tmpdir = tempfile.mkdtemp(prefix="vq_")

    try:
        # Oracle k-mer set
        oracle = build_oracle(args.ref, k)
        checks["oracle_size"] = len(oracle)

        # Controlled query sets
        rng = random.Random(args.seed)

        # Positive queries, k-mers taken from the reference
        pos_kmers = list(oracle)
        if len(pos_kmers) > args.max_query_kmers:
            pos_kmers = rng.sample(pos_kmers, args.max_query_kmers)

        # Negative queries, random k-mers absent from the reference
        neg_kmers = set()
        while len(neg_kmers) < min(len(pos_kmers), 5000):
            kmer = random_dna(k, rng)
            can = canonicalize(kmer)
            if can not in oracle:
                neg_kmers.add(can)
        neg_kmers = list(neg_kmers)

        # Near-miss negative queries, single substitutions of reference k-mers
        off1_kmers = set()
        bases = "ACGT"
        sample_pos = rng.sample(pos_kmers, min(2000, len(pos_kmers)))
        for kmer in sample_pos:
            for pos in range(k):
                for b in bases:
                    if b != kmer[pos]:
                        mutated = kmer[:pos] + b + kmer[pos + 1 :]
                        can = canonicalize(mutated)
                        if can not in oracle:
                            off1_kmers.add(can)
                            if len(off1_kmers) >= 2000:
                                break
                if len(off1_kmers) >= 2000:
                    break
            if len(off1_kmers) >= 2000:
                break
        off1_kmers = list(off1_kmers)

        # Query FASTA files
        pos_fa = os.path.join(tmpdir, "pos_kmers.fa")
        neg_fa = os.path.join(tmpdir, "neg_kmers.fa")
        off1_fa = os.path.join(tmpdir, "off1_kmers.fa")
        mixed_fa = os.path.join(tmpdir, "mixed_kmers.fa")

        # Probes are written in random orientation, half of them as the reverse complement of the canonical form, so a lookup that finds k-mers only in canonical orientation cannot pass.
        def orient(kms):
            """kms with each k-mer reverse-complemented with probability one half."""
            return [rc(x) if rng.random() < 0.5 else x for x in kms]

        make_single_kmer_fasta(orient(pos_kmers), pos_fa)
        make_single_kmer_fasta(orient(neg_kmers), neg_fa)
        make_single_kmer_fasta(orient(off1_kmers), off1_fa)
        make_single_kmer_fasta(orient(pos_kmers + neg_kmers), mixed_fa)

        # Streaming records are reference windows from both strands with a few substitutions, so each record mixes present and absent k-mers and exercises the multi-k-mer extension path rather than single-k-mer lookups.
        # The expected counts come from the oracle k-mer by k-mer.
        refs = [sq for _, sq in read_fasta(args.ref) if len(sq) >= 200]
        stream_fa = os.path.join(tmpdir, "stream_reads.fa")
        s_total = s_pos = 0
        with open(stream_fa, "w") as fh:
            for i in range(2000):
                sq = rng.choice(refs)
                p0 = rng.randrange(0, len(sq) - 150 + 1)
                r = list(sq[p0 : p0 + 150])
                for _ in range(rng.randrange(0, 4)):
                    j = rng.randrange(len(r))
                    r[j] = rng.choice("ACGT".replace(r[j].upper(), "") or "ACGT")
                r = "".join(r)
                if rng.random() < 0.5:
                    r = rc(r)
                fh.write(f">s{i}\n{r}\n")
                ks = extract_kmers(r, k)
                s_total += len(ks)
                s_pos += sum(1 for x in ks if x in oracle)

        sets = {
            "positive": (pos_fa, len(pos_kmers), len(pos_kmers)),
            "negative": (neg_fa, len(neg_kmers), 0),
            "off1_negative": (off1_fa, len(off1_kmers), 0),
            "mixed": (mixed_fa, len(pos_kmers) + len(neg_kmers), len(pos_kmers)),
            "stream": (stream_fa, s_total, s_pos),
        }

        # The official SSHash tool on every set, with its actual reported counts
        ss_counts = {}
        if args.sshash and os.path.isfile(args.sshash):
            ss_idx = os.path.join(tmpdir, "sshash.idx")
            m = min(k - 1, 15)
            b = subprocess.run(
                [
                    args.sshash,
                    "build",
                    "-i",
                    args.ref,
                    "-o",
                    ss_idx,
                    "-k",
                    str(k),
                    "-m",
                    str(m),
                    "-d",
                    tmpdir,
                ],
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            checks["sshash_build"] = {
                "pass": b.returncode == 0,
                "cmd": f"sshash build -i REF -k {k} -m {m}",
                "error": b.stderr[-300:] if b.returncode else "",
            }
            if b.returncode == 0:
                for name, (fa, _, _) in sets.items():
                    q = subprocess.run(
                        [args.sshash, "query", "-i", ss_idx, "-q", fa],
                        capture_output=True,
                        text=True,
                        timeout=600,
                        check=False,
                    )
                    ss_counts[name] = parse_counts(q.stdout)

        # optik with both index formats on every set, against the oracle and SSHash
        for fmt, build_cmd in (("spss", "build-spss"), ("fm", "build-fm")):
            idx_path = os.path.join(tmpdir, f"optik_{fmt}.idx")
            b = subprocess.run(
                [optik, build_cmd, "-i", args.ref, "-k", str(k), "-o", idx_path],
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            if b.returncode != 0:
                checks[f"{fmt}_build"] = {"pass": False, "error": b.stderr[-500:]}
                continue
            for name, (fa, n_kmers, expected) in sets.items():
                result, err = run_optik_query(idx_path, fa, k, optik)
                key = f"{fmt}_{name}_query"
                if err:
                    checks[key] = {"pass": False, "error": err}
                    continue
                got_pos = result.get("num_positive_kmers")
                got_n = result.get("num_kmers")
                ss = ss_counts.get(name)
                ok = got_pos == expected and got_n == n_kmers
                c = {
                    "expected_positive": expected,
                    "optik_positive": got_pos,
                    "optik_kmers": got_n,
                    "expected_kmers": n_kmers,
                }
                if ss is not None:
                    c["sshash_positive"] = ss.get("num_positive_kmers")
                    c["sshash_kmers"] = ss.get("num_kmers")
                    ok = (
                        ok
                        and ss.get("num_positive_kmers") == got_pos
                        and ss.get("num_kmers") == got_n
                    )
                c["pass"] = ok
                c["detail"] = (
                    "optik == oracle and == SSHash's reported counts"
                    if ss is not None
                    else "optik == oracle"
                )
                checks[key] = c

    finally:
        import shutil

        shutil.rmtree(tmpdir, ignore_errors=True)

    return checks


def main():
    """Run the query checks and print their report."""
    parser = argparse.ArgumentParser(description="Validate optik query correctness")
    parser.add_argument("-r", dest="ref", required=True, help="Reference FASTA")
    parser.add_argument("-k", type=int, required=True, help="K-mer size")
    parser.add_argument(
        "--optik", default="target/release/optik", help="Path to optik binary"
    )
    parser.add_argument(
        "--sshash", default=None, help="Path to SSHash binary (optional)"
    )
    parser.add_argument(
        "--max-query-kmers",
        type=int,
        default=50000,
        help="Max k-mers in positive query set",
    )
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    if not os.path.isfile(args.optik):
        print(f"ERROR: optik binary not found at {args.optik}", file=sys.stderr)
        sys.exit(1)

    checks = validate(args)

    # Overall verdict
    all_pass = all(
        c.get("pass", False)
        for c in checks.values()
        if isinstance(c, dict) and "pass" in c
    )
    overall = "PASS" if all_pass else "FAIL"

    # Human-readable summary
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  Query Validation: {overall}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    for name, check in checks.items():
        if not isinstance(check, dict) or "pass" not in check:
            continue
        status = "PASS" if check.get("pass", False) else "FAIL"
        line = f"  [{status}] {name}"
        if "expected" in check and "reported" in check:
            line += f"  (expected={check['expected']}, reported={check['reported']})"
        if "detail" in check:
            line += f"  -- {check['detail']}"
        if check.get("error"):
            line += f"  ERROR: {check['error']}"
        print(line, file=sys.stderr)
    print(file=sys.stderr)

    # JSON report
    report = {
        "overall": overall,
        "checks": checks,
    }
    print(json.dumps(report, indent=2))

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
