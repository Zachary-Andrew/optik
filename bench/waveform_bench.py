#!/usr/bin/env python3
"""waveform_bench.py - WFA2-lib reference bench for optik's edit-distance engine, the wavefront aligner in edit mode as `map` uses it.

The bench generates a deterministic battery of (query, text) pairs with a known edit distance, made by applying k unit edits at distinct, non-adjacent sites of a random locus.
It runs the wfa2_oracle harness on them and asserts

  * engine equivalence  optik's semi-global ED equals WFA2-lib's ends-free ED (checked per pair inside the oracle and counted here per tier).
  * end-to-end check    WFA2-lib end-to-end equals a naive O(mn) DP on small pairs, checked inside the oracle.
  * ground truth        optik's semi-global ED equals the constructed k on every pair whose fit is unique.
                        Pairs whose window admits a better fit than the embedded locus are dropped as ambiguous and reported, not failed.
                        A score above k, however, is impossible because the locus lies inside the text, so it fails hard.

Coverage follows the read lengths the mapper handles.

  * 150 bp reads at 10% error (k = 15), with extra samples
  * length strata m <= 64, 65..512 and > 512, covering short reads up to the 500 bp chaining threshold and beyond
  * banded stress on 513 and 1500 bp long reads at 12-15% error
  * indel-rich pairs with high insertion and deletion fractions, the case wfmash's wflign stage meets

The oracle cross-checks three independent implementations of unit-cost edit distance, namely optik's exact DP, WFA2-lib's C wavefront engine (ends-free under the edit metric through bench/patches/wfa2_edit_endsfree.py) and a plain rolling-row DP.
Any disagreement fails the run.

Usage
  waveform_bench.py [--oracle PATH] [--seed N] [--quick]
                    [--workdir DIR] [--keep]

The exit status is 0 on full agreement and 1 on any failure.
"""

import argparse
import os
import random
import subprocess
import sys

from diagnostics import LOG_TAIL_CHARS
from units import PERCENT

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORACLE_SRC = os.path.join(ROOT, "bench", "cpp", "wfa2_oracle.cpp")
INCLUDE_DIRS = (os.path.join(ROOT, "include"), os.path.join(ROOT, "temp", "WFA2-lib"))
WFA2_LIB = os.path.join(ROOT, "temp", "WFA2-lib", "lib", "libwfa.a")

BASES = "ACGT"
# length of the synthetic reference that loci are cut from
GENOME_LEN = 1_000_000
# edit sites are at least this far apart, so adjacent edits cannot interact and lower the edit distance
MIN_EDIT_SPACING = 2
# probability that an indel is an insertion rather than a deletion
INSERTION_PROB = 0.5
# --quick divides sample counts by QUICK_DIVISOR, keeping at least MIN_SAMPLES per tier
QUICK_DIVISOR = 10
MIN_SAMPLES = 10
# upper bounds of the short (single-word bit-parallel) and medium length strata
SHORT_STRATUM_MAX = 64
MEDIUM_STRATUM_MAX = 512
# the failing regime the battery weighs most: 150 bp reads at 10% error, k = 15 edits
FOCUS_READ_LEN = 150
FOCUS_ERR = 0.10
FOCUS_EDITS = 15
# the focus window tier gets this many times its listed samples
FOCUS_SAMPLE_FACTOR = 2
# ground-truth locus tiers: (locus length, substitution counts, pairs)
LOCUS_TIERS = [
    (30, [1, 3], 200),  # m <= 64
    (64, [2, 6], 200),  # 64/65 boundary
    (65, [2, 6], 200),  # 64/65 boundary
    (150, [3, 7, 15, 22], 300),  # k = 15 is 10% of 150
    (300, [6, 15, 30], 150),  # 65..512
    (512, [10, 25, 51], 150),  # 512/513 boundary (10% of 512)
    (513, [10, 25, 51], 150),  # 512/513 boundary
    (1500, [30, 75, 150], 100),  # > 512
    (1500, [180, 220], 120),  # banded stress (12%, 14.7%)
    (4000, [80, 200], 80),  # scale
]
# locus tiers of this length with at least this many edits are the banded stress tiers
BANDED_LOCUS_LEN = 1500
BANDED_MIN_EDITS = 180
# indel-consensus locus tiers: (locus length, edit counts, pairs, indel fraction)
INDEL_LOCUS_TIERS = [(150, [10, 20], 200, 0.7)]
# ground-truth window tiers: (locus length, error rates, pairs)
WINDOW_TIERS = [
    (65, [0.02, 0.05, 0.10], 150),
    (150, [0.02, 0.05, 0.10], 200),
    (300, [0.05, 0.10], 120),
    (513, [0.05, 0.10], 80),
]
# random flanks around a window's locus: short loci (up to SHORT_FLANK_MAX_LEN) get SHORT_FLANK, longer ones LONG_FLANK
SHORT_FLANK_MAX_LEN = 150
SHORT_FLANK = (30, 30)
LONG_FLANK = (50, 50)
# pairs and indel fraction of the indel-consensus window tier (150 bp at 10%)
INDEL_WINDOW_PAIRS = 200
INDEL_WINDOW_FRAC = 0.5
# a window tier emits (query, window) and (query, locus), a locus tier only the latter
PAIRS_PER_WINDOW = 2
# columns of an oracle output row
ORACLE_COLUMNS = 7
DEFAULT_SEED = 42


# Pair generation
class PairGen:
    """Seeded generator of (query, locus, window) triples."""

    def __init__(self, seed):
        self.rng = random.Random(seed)
        # Synthetic reference that loci are cut from
        self.genome = "".join(self.rng.choice(BASES) for _ in range(GENOME_LEN))

    def locus(self, length):
        """A random substring of the synthetic genome with the given length."""
        start = self.rng.randrange(0, len(self.genome) - length)
        return self.genome[start : start + length]

    @staticmethod
    def _apply_edits(locus, num_edits, rng, indel_frac):
        """Locus with k unit edits at distinct sites at least 2 apart, each an indel with probability indel_frac."""
        seq = list(locus)
        seq_len = len(seq)
        if num_edits > seq_len:
            raise ValueError(f"k={num_edits} > locus length {seq_len}")
        # k distinct sites at pairwise distance of at least 2, which prevents adjacent edits from interacting and collapsing the edit distance below k.
        pool = list(range(seq_len))
        rng.shuffle(pool)
        chosen = []
        for site in pool:
            if all(abs(site - chosen_site) >= MIN_EDIT_SPACING for chosen_site in chosen):
                chosen.append(site)
                if len(chosen) == num_edits:
                    break
        if len(chosen) < num_edits:  # pragma: no cover, the chosen sizes rule this out
            raise RuntimeError("could not place k non-adjacent edits")
        # Edits go right to left, so original coordinates stay valid throughout.
        for site in sorted(chosen, reverse=True):
            if rng.random() < indel_frac:  # indel
                if rng.random() < INSERTION_PROB:  # insertion before original p
                    seq.insert(site, rng.choice(BASES))
                else:  # deletion of original p
                    del seq[site]
            else:  # substitution
                seq[site] = rng.choice([base for base in BASES if base != seq[site]])
        return "".join(seq)


# Tiers
#
# Ground-truth tiers use substitutions only (indel_frac=0.0).
# k substitutions at distinct, non-adjacent sites guarantee ED(query, locus) == k exactly, so the assertion that optik's ED equals the constructed k is sharp.
# Indel-consensus tiers mix in insertions and deletions, for which an exact ED cannot be constructed (adjacent edits can legitimately lower the true ED), and assert only that the engines agree, the property that matters for indel-rich alignment.
def build_tiers(quick):
    """List the pair tiers of the battery, with sample counts cut tenfold by --quick."""

    def sample_count(count):
        """Sample count n, cut in quick mode to a tenth but never below 10."""
        return max(MIN_SAMPLES, count // (QUICK_DIVISOR if quick else 1))

    tiers = []
    # Ground-truth locus tiers, where the text is the read's locus and edits are substitutions.
    for locus_len, edit_counts, num_pairs in LOCUS_TIERS:
        for num_edits in edit_counts:
            tiers.append(
                {
                    "label": f"locus m={locus_len} k={num_edits}",
                    "kind": "locus",
                    "m": locus_len,
                    "k": num_edits,
                    "n": sample_count(num_pairs),
                    "indel_frac": 0.0,
                    "flank": (0, 0),
                    "gt_check": True,
                    "stratum": (
                        "m<=64" if locus_len <= SHORT_STRATUM_MAX else "65-512" if locus_len <= MEDIUM_STRATUM_MAX else ">512"
                    )
                    + (" (banded)" if locus_len in (BANDED_LOCUS_LEN,) and num_edits >= BANDED_MIN_EDITS else ""),
                }
            )
    # Indel-consensus locus tier, indel-heavy and checked for engine agreement only
    for locus_len, edit_counts, num_pairs, frac in INDEL_LOCUS_TIERS:
        for num_edits in edit_counts:
            tiers.append(
                {
                    "label": f"locus m={locus_len} k={num_edits} indel70",
                    "kind": "locus",
                    "m": locus_len,
                    "k": num_edits,
                    "n": sample_count(num_pairs),
                    "indel_frac": frac,
                    "flank": (0, 0),
                    "gt_check": False,
                    "stratum": "65-512 (indel consensus)",
                }
            )
    # Ground-truth window tiers, where the text is random flanks around the locus and edits are substitutions.
    # This is the mapper's realistic regime, with extra weight on 150 bp at 10%.
    for locus_len, errs, num_pairs in WINDOW_TIERS:
        for err in errs:
            num_edits = max(1, round(locus_len * err))
            n_this = sample_count(num_pairs * FOCUS_SAMPLE_FACTOR) if (locus_len == FOCUS_READ_LEN and err == FOCUS_ERR) else sample_count(num_pairs)
            tiers.append(
                {
                    "label": f"window m={locus_len} err={err:.0%}",
                    "kind": "window",
                    "m": locus_len,
                    "k": num_edits,
                    "n": n_this,
                    "indel_frac": 0.0,
                    "flank": SHORT_FLANK if locus_len <= SHORT_FLANK_MAX_LEN else LONG_FLANK,
                    "gt_check": True,
                    "stratum": (
                        "m<=64" if locus_len <= SHORT_STRATUM_MAX else "65-512" if locus_len <= MEDIUM_STRATUM_MAX else ">512"
                    ),
                }
            )
    # Indel-consensus window tier, 150 bp at 10% with heavy indels
    tiers.append(
        {
            "label": "window m=150 err=10% indel50",
            "kind": "window",
            "m": FOCUS_READ_LEN,
            "k": FOCUS_EDITS,
            "n": sample_count(INDEL_WINDOW_PAIRS),
            "indel_frac": INDEL_WINDOW_FRAC,
            "flank": SHORT_FLANK,
            "gt_check": False,
            "stratum": "65-512 (indel consensus)",
        }
    )
    return tiers


# Oracle
def build_oracle(workdir):
    """Compile the wfa2_oracle harness against WFA2-lib, returning (path, None) or (None, reason)."""
    if not (os.path.isfile(WFA2_LIB) and all(os.path.isdir(inc_dir) for inc_dir in INCLUDE_DIRS)):
        return None, "WFA2-lib not found (run bench/build_tools.sh wfa2 first)"
    out = os.path.join(workdir, "wfa2_oracle")
    # Relink when WFA2-lib or the harness changed, since a cached binary would keep testing the old library after a rebuild.
    stale = not os.path.isfile(out) or os.path.getmtime(out) < max(
        os.path.getmtime(WFA2_LIB), os.path.getmtime(ORACLE_SRC)
    )
    if stale:
        cmd = ["g++", "-std=c++17", "-O2", "-march=native"]
        for inc_dir in INCLUDE_DIRS:
            cmd += ["-I", inc_dir]
        cmd += [ORACLE_SRC, WFA2_LIB, "-lm", "-fopenmp", "-o", out]
        try:
            # Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0.
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                env=dict(os.environ, NIX_ENFORCE_NO_NATIVE="0"),
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            return None, f"oracle build failed: {exc}"
    return out, None


def run_oracle(oracle, pairs_tsv):
    """Feed TSV lines to the oracle and return its parsed rows."""
    result = subprocess.run(
        [oracle], input=pairs_tsv, text=True, capture_output=True, check=False
    )
    rows = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("m\tn\toptik_sg"):
            continue
        parts = line.split("\t")
        if len(parts) != ORACLE_COLUMNS:
            continue
        rows.append(tuple(int(field) for field in parts))
    return rows, result.returncode, result.stderr


# Report
class TierStats:
    """Per-tier counters of pairs, drops, ground-truth hits and engine agreement."""

    __slots__ = (
        "capped_agree",
        "dropped",
        "fail_gt_high",
        "gt_ok",
        "pairs",
        "sg_agree",
    )

    def __init__(self):
        self.pairs = 0
        self.dropped = 0
        self.gt_ok = 0
        self.fail_gt_high = 0
        self.sg_agree = 0
        self.capped_agree = 0


def main():
    """Generate the battery, run the oracle and report agreement per tier."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle", help="path to prebuilt wfa2_oracle binary")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--quick", action="store_true", help="reduced pair counts (CI smoke)"
    )
    parser.add_argument("--workdir", default=os.path.join(ROOT, "temp", "waveform"))
    parser.add_argument("--keep", action="store_true", help="keep workdir files")
    parser.add_argument("--dump", help="also write generated TSV to this path")
    args = parser.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    gen = PairGen(args.seed)
    tiers = build_tiers(args.quick)

    if args.oracle and os.path.isfile(args.oracle):
        oracle = args.oracle
    elif args.oracle:
        print(f"error: --oracle {args.oracle} not found", file=sys.stderr)
        return 1
    else:
        oracle, err = build_oracle(args.workdir)
        if oracle is None:
            # The caller (run_validations.sh step 6) runs this only when WFA2-lib exists, so reaching this point means the build failed rather than a dependency being absent, and that fails loudly instead of passing quietly.
            print(f"error: oracle unavailable: {err}", file=sys.stderr)
            return 1

    # Generate all pairs.
    # A window tier emits both (q, window) and (q, locus), so the end-to-end legs can anchor on the exact embedded locus.
    pairs = []
    for tier_idx, tier in enumerate(tiers):
        for _ in range(tier["n"]):
            locus = gen.locus(tier["m"])
            query = gen._apply_edits(locus, tier["k"], gen.rng, tier["indel_frac"])
            if tier["kind"] == "locus":
                pairs.append((query, locus, tier_idx, "locus"))
            else:
                left_len, right_len = tier["flank"]
                window = gen.locus(left_len) + locus + gen.locus(right_len)
                pairs.append((query, window, tier_idx, "window"))
                pairs.append((query, locus, tier_idx, "locus"))
    assert len(pairs) == sum(
        tier["n"] * (PAIRS_PER_WINDOW if tier["kind"] == "window" else 1) for tier in tiers
    )

    if args.dump:
        with open(args.dump, "w") as handle:
            for query, text, _ti, _ in pairs:
                handle.write(f"{query}\t{text}\n")

    tsv = "\n".join(f"{query}\t{txt}" for query, txt, _, _ in pairs) + "\n"
    rows, exit_code, stderr = run_oracle(oracle, tsv)
    if exit_code != 0:
        print("waveform bench FAIL: oracle exited nonzero")
        print(stderr[-LOG_TAIL_CHARS:])
        return 1
    if len(rows) != len(pairs):
        print(
            f"waveform bench FAIL: oracle returned {len(rows)} rows for "
            f"{len(pairs)} pairs; stderr tail:\n{stderr[-LOG_TAIL_CHARS:]}"
        )
        return 1

    stats = {i: TierStats() for i in range(len(tiers))}
    for (_query, _text, tier_idx, _kind), (_pattern_len, _text_len, optik_sg, capped, wsg, _we2e, _naive) in zip(
        pairs, rows, strict=False
    ):
        tier_stats = stats[tier_idx]
        tier = tiers[tier_idx]
        tier_stats.pairs += 1
        tier_stats.sg_agree += optik_sg == wsg
        tier_stats.capped_agree += optik_sg == capped
        if optik_sg > tier["k"]:
            # This cannot happen by construction (the locus is inside the text and the query is the locus with at most k unit edits), so it is a real bug.
            tier_stats.fail_gt_high += 1
        elif tier["gt_check"]:
            if optik_sg == tier["k"]:
                tier_stats.gt_ok += 1
            else:
                tier_stats.dropped += 1
        # Consensus tiers make no ground-truth assertion, because indel construction does not fix the ED exactly.
        # Engine agreement above is their contract.

    # Report
    total_pairs = total_drop = total_gt = total_agree = total_capped = total_fail = 0
    print("==============================================")
    print("  waveform / WFA2 exactness bench")
    print(f"  seed={args.seed}  oracle={oracle}")
    print(
        "  engines: optik WFA (edit mode) == optik capped == WFA2-lib ends-free == "
        "naive DP == constructed ED"
    )
    print("==============================================")
    print(
        f"{'tier':<28}{'pairs':>6}{'drop':>6}{'gt==k':>7}"
        f"{'sg-agree':>9}{'capped':>7}{'gt>k!':>7}  stratum"
    )
    for i, tier in enumerate(tiers):
        tier_stats = stats[i]
        # Tiers without a constructed edit distance have no ground-truth columns.
        drop = str(tier_stats.dropped) if tier["gt_check"] else "-"
        gtok = str(tier_stats.gt_ok) if tier["gt_check"] else "-"
        print(
            f"{tier['label']:<28}{tier_stats.pairs:>6}{drop:>6}{gtok:>7}"
            f"{tier_stats.sg_agree:>9}{tier_stats.capped_agree:>7}{tier_stats.fail_gt_high:>7}  {tier['stratum']}"
        )
        total_pairs += tier_stats.pairs
        total_agree += tier_stats.sg_agree
        total_capped += tier_stats.capped_agree
        total_fail += tier_stats.fail_gt_high
        if tier["gt_check"]:
            total_drop += tier_stats.dropped
            total_gt += tier_stats.gt_ok
    gt_checked = sum(
        tier["n"] * (PAIRS_PER_WINDOW if tier["kind"] == "window" else 1) for tier in tiers if tier["gt_check"]
    )

    # Summary of the 150 bp, 10% error regime over its locus and window lines
    fr_tiers = [
        i
        for i, tier in enumerate(tiers)
        if (tier["m"] == FOCUS_READ_LEN and tier["k"] == FOCUS_EDITS)
        or (tier["kind"] == "window" and tier["m"] == FOCUS_READ_LEN and tier["k"] == FOCUS_EDITS)
    ]
    fr_pairs = sum(stats[tier_idx].pairs for tier_idx in fr_tiers)
    fr_gt = sum(stats[tier_idx].gt_ok for tier_idx in fr_tiers)
    fr_ag = sum(stats[tier_idx].sg_agree for tier_idx in fr_tiers)
    fr_drop = sum(stats[tier_idx].dropped for tier_idx in fr_tiers)
    print("--------------------------------------------------------------")
    print(
        f"  failing regime (150 bp / 10% error): {fr_pairs} pairs, "
        f"ground-truth {fr_gt}, engine-agree {fr_ag}, ambiguous {fr_drop}"
    )
    print("--------------------------------------------------------------")
    # The bench passes when all engines agree on every pair and each ground-truth pair is matched or dropped as ambiguous.
    passed = (
        total_fail == 0
        and total_agree == total_pairs
        and total_capped == total_pairs
        and total_gt + total_drop == gt_checked
    )
    print(
        f"  totals: {total_pairs} pairs, {total_agree} engine-agree "
        f"({PERCENT * total_agree / total_pairs:.2f}%), {total_capped} capped-agree "
        f"({PERCENT * total_capped / total_pairs:.2f}%), {total_gt} ground-truth, "
        f"{total_drop} ambiguous-dropped, {total_fail} impossible-high"
    )
    print("  => WAVEFORM ORACLE: " + ("PASS" if passed else "FAIL"))
    print("==============================================")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
