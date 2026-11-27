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

from diagnostics import DETAIL_CHARS
from dna import read_fasta, revcomp

# timeouts of the external commands, in seconds
QUERY_TIMEOUT_S = 300
BUILD_TIMEOUT_S = 1800
SSHASH_QUERY_TIMEOUT_S = 600
# most random absent k-mers in the negative set
MAX_NEGATIVE_KMERS = 5000
# reference k-mers mutated for the near-miss set, and the most near-miss k-mers kept
NEAR_MISS_SOURCE_KMERS = 2000
MAX_NEAR_MISS_KMERS = 2000
# probability that a probe or streaming record is written reverse-complemented
REVCOMP_PROB = 0.5
# streaming records: count, length, the shortest reference they are drawn from, and substitutions per record (exclusive upper bound)
STREAM_RECORDS = 2000
STREAM_READ_LEN = 150
STREAM_MIN_REF_LEN = 200
STREAM_MAX_SUBS = 4
# SSHash minimizer length, at most k - 1
SSHASH_MINIMIZER_LEN = 15
# argparse defaults
DEFAULT_MAX_QUERY_KMERS = 50000
DEFAULT_SEED = 42

# k-mer helpers


def canonicalize(kmer: str) -> str:
    """Canonical form of kmer, the smaller of it and its reverse complement."""
    rev = revcomp(kmer)
    return min(kmer, rev)


def extract_kmers(seq: str, kmer_len: int):
    """Canonical k-mers of seq, case-insensitive, where any non-ACGT base breaks the run as it does in optik and SSHash."""
    kmers = []
    run = []
    for base in seq.upper() + "N":
        if base in "ACGT":
            run.append(base)
            continue
        run_seq = "".join(run)
        kmers.extend(canonicalize(run_seq[j : j + kmer_len]) for j in range(len(run_seq) - kmer_len + 1))
        run = []
    return kmers


def build_oracle(ref_path: str, kmer_len: int):
    """Set of all canonical k-mers of the reference."""
    oracle = set()
    for _, seq in read_fasta(ref_path):
        oracle.update(extract_kmers(seq, kmer_len))
    return oracle


def random_dna(length: int, rng: random.Random) -> str:
    """Random ACGT string of the given length."""
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_single_kmer_fasta(kmers, path):
    """Write a FASTA file with one k-mer per record."""
    with open(path, "w") as handle:
        handle.writelines(f">q{i}\n{kmer}\n" for i, kmer in enumerate(kmers))


def run_optik_query(index_path: str, queries_fa: str, kmer_len: int, optik_bin: str):
    """Run optik query and parse its summary output."""
    cmd = [optik_bin, "query", "-i", index_path, "-q", queries_fa, "-k", str(kmer_len)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=QUERY_TIMEOUT_S, check=False)
    if proc.returncode != 0:
        return None, f"optik query failed: {proc.stderr[:DETAIL_CHARS]}"

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
            pass  # not a count line
    return out


def validate(args):
    """Run every query check against the oracle and SSHash and return the report dict."""
    checks = {}
    kmer_len = args.k
    optik = args.optik
    tmpdir = tempfile.mkdtemp(prefix="vq_")

    try:
        # Oracle k-mer set
        oracle = build_oracle(args.ref, kmer_len)
        checks["oracle_size"] = len(oracle)

        # Controlled query sets
        rng = random.Random(args.seed)

        # Positive queries, k-mers taken from the reference
        pos_kmers = list(oracle)
        if len(pos_kmers) > args.max_query_kmers:
            pos_kmers = rng.sample(pos_kmers, args.max_query_kmers)

        # Negative queries, random k-mers absent from the reference
        neg_kmers = set()
        while len(neg_kmers) < min(len(pos_kmers), MAX_NEGATIVE_KMERS):
            kmer = random_dna(kmer_len, rng)
            can = canonicalize(kmer)
            if can not in oracle:
                neg_kmers.add(can)
        neg_kmers = list(neg_kmers)

        # Near-miss negative queries, single substitutions of reference k-mers
        off1_kmers = set()
        bases = "ACGT"
        sample_pos = rng.sample(pos_kmers, min(NEAR_MISS_SOURCE_KMERS, len(pos_kmers)))
        for kmer in sample_pos:
            for pos in range(kmer_len):
                for base in bases:
                    if base != kmer[pos]:
                        mutated = kmer[:pos] + base + kmer[pos + 1 :]
                        can = canonicalize(mutated)
                        if can not in oracle:
                            off1_kmers.add(can)
                            if len(off1_kmers) >= MAX_NEAR_MISS_KMERS:
                                break
                if len(off1_kmers) >= MAX_NEAR_MISS_KMERS:
                    break
            if len(off1_kmers) >= MAX_NEAR_MISS_KMERS:
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
            return [revcomp(kmer) if rng.random() < REVCOMP_PROB else kmer for kmer in kms]

        make_single_kmer_fasta(orient(pos_kmers), pos_fa)
        make_single_kmer_fasta(orient(neg_kmers), neg_fa)
        make_single_kmer_fasta(orient(off1_kmers), off1_fa)
        make_single_kmer_fasta(orient(pos_kmers + neg_kmers), mixed_fa)

        # Streaming records are reference windows from both strands with a few substitutions, so each record mixes present and absent k-mers and exercises the multi-k-mer extension path rather than single-k-mer lookups.
        # The expected counts come from the oracle k-mer by k-mer.
        refs = [seq for _, seq in read_fasta(args.ref) if len(seq) >= STREAM_MIN_REF_LEN]
        stream_fa = os.path.join(tmpdir, "stream_reads.fa")
        s_total = s_pos = 0
        with open(stream_fa, "w") as handle:
            for i in range(STREAM_RECORDS):
                seq = rng.choice(refs)
                start = rng.randrange(0, len(seq) - STREAM_READ_LEN + 1)
                read = list(seq[start : start + STREAM_READ_LEN])
                for _ in range(rng.randrange(0, STREAM_MAX_SUBS)):
                    mut_pos = rng.randrange(len(read))
                    read[mut_pos] = rng.choice("ACGT".replace(read[mut_pos].upper(), "") or "ACGT")
                read = "".join(read)
                if rng.random() < REVCOMP_PROB:
                    read = revcomp(read)
                handle.write(f">s{i}\n{read}\n")
                kmers = extract_kmers(read, kmer_len)
                s_total += len(kmers)
                s_pos += sum(1 for kmer in kmers if kmer in oracle)

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
            minimizer_len = min(kmer_len - 1, SSHASH_MINIMIZER_LEN)
            base = subprocess.run(
                [
                    args.sshash,
                    "build",
                    "-i",
                    args.ref,
                    "-o",
                    ss_idx,
                    "-k",
                    str(kmer_len),
                    "-m",
                    str(minimizer_len),
                    "-d",
                    tmpdir,
                ],
                capture_output=True,
                text=True,
                timeout=BUILD_TIMEOUT_S,
                check=False,
            )
            checks["sshash_build"] = {
                "pass": base.returncode == 0,
                "cmd": f"sshash build -i REF -k {kmer_len} -m {minimizer_len}",
                "error": base.stderr[-DETAIL_CHARS:] if base.returncode else "",
            }
            if base.returncode == 0:
                for name, (query_fa, _, _) in sets.items():
                    proc = subprocess.run(
                        [args.sshash, "query", "-i", ss_idx, "-q", query_fa],
                        capture_output=True,
                        text=True,
                        timeout=SSHASH_QUERY_TIMEOUT_S,
                        check=False,
                    )
                    ss_counts[name] = parse_counts(proc.stdout)

        # optik with both index formats on every set, against the oracle and SSHash
        for fmt, build_cmd in (("spss", "build-spss"), ("fm", "build-fm")):
            idx_path = os.path.join(tmpdir, f"optik_{fmt}.idx")
            base = subprocess.run(
                [optik, build_cmd, "-i", args.ref, "-k", str(kmer_len), "-o", idx_path],
                capture_output=True,
                text=True,
                timeout=BUILD_TIMEOUT_S,
                check=False,
            )
            if base.returncode != 0:
                checks[f"{fmt}_build"] = {"pass": False, "error": base.stderr[-DETAIL_CHARS:]}
                continue
            for name, (query_fa, n_kmers, expected) in sets.items():
                result, err = run_optik_query(idx_path, query_fa, kmer_len, optik)
                key = f"{fmt}_{name}_query"
                if err:
                    checks[key] = {"pass": False, "error": err}
                    continue
                got_pos = result.get("num_positive_kmers")
                got_n = result.get("num_kmers")
                ss_count = ss_counts.get(name)
                passed = got_pos == expected and got_n == n_kmers
                check = {
                    "expected_positive": expected,
                    "optik_positive": got_pos,
                    "optik_kmers": got_n,
                    "expected_kmers": n_kmers,
                }
                if ss_count is not None:
                    check["sshash_positive"] = ss_count.get("num_positive_kmers")
                    check["sshash_kmers"] = ss_count.get("num_kmers")
                    passed = (
                        passed
                        and ss_count.get("num_positive_kmers") == got_pos
                        and ss_count.get("num_kmers") == got_n
                    )
                check["pass"] = passed
                check["detail"] = (
                    "optik == oracle and == SSHash's reported counts"
                    if ss_count is not None
                    else "optik == oracle"
                )
                checks[key] = check

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
        default=DEFAULT_MAX_QUERY_KMERS,
        help="Max k-mers in positive query set",
    )
    parser.add_argument("-s", "--seed", type=int, default=DEFAULT_SEED, help="Random seed")
    args = parser.parse_args()

    if not os.path.isfile(args.optik):
        print(f"ERROR: optik binary not found at {args.optik}", file=sys.stderr)
        sys.exit(1)

    checks = validate(args)

    # Overall verdict
    all_pass = all(
        check.get("pass", False)
        for check in checks.values()
        if isinstance(check, dict) and "pass" in check
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
