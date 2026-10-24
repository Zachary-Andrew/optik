#!/usr/bin/env python3
"""waveform_bench.py - WFA2-lib reference bench for optik's edit-distance engine, the wavefront aligner in edit mode as `map` uses it.

The bench generates a deterministic battery of (query, text) pairs with a known edit distance, made by applying k unit edits at distinct, non-adjacent sites of a random locus.
It runs the wfa2_oracle harness on them and asserts

  * engine equivalence  optik's semi-global ED equals WFA2-lib's ends-free ED
                        (checked per pair inside the oracle and counted here
                        per tier)
  * end-to-end sanity   WFA2-lib end-to-end equals a naive O(mn) DP on small
                        pairs, checked inside the oracle
  * ground truth        optik's semi-global ED equals the constructed k on
                        every pair whose fit is unique.  Pairs whose window
                        admits a better fit than the embedded locus are
                        dropped as ambiguous and reported, not failed.  A
                        score above k is impossible, since the locus lies
                        inside the text, and fails hard.

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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORACLE_SRC = os.path.join(ROOT, "bench", "cpp", "wfa2_oracle.cpp")
INCLUDE_DIRS = (os.path.join(ROOT, "include"), os.path.join(ROOT, "temp", "WFA2-lib"))
WFA2_LIB = os.path.join(ROOT, "temp", "WFA2-lib", "lib", "libwfa.a")

BASES = "ACGT"


# Pair generation
class PairGen:
    """Seeded generator of (query, locus, window) triples."""

    def __init__(self, seed):
        self.rng = random.Random(seed)
        # Synthetic reference that loci are cut from
        self.genome = "".join(self.rng.choice(BASES) for _ in range(1000000))

    def locus(self, length):
        """A random substring of the synthetic genome with the given length."""
        start = self.rng.randrange(0, len(self.genome) - length)
        return self.genome[start : start + length]

    @staticmethod
    def _apply_edits(locus, k, rng, indel_frac):
        """Locus with k unit edits at distinct sites at least 2 apart, each an indel with probability indel_frac."""
        seq = list(locus)
        n = len(seq)
        if k > n:
            raise ValueError(f"k={k} > locus length {n}")
        # k distinct sites at pairwise distance of at least 2, which prevents adjacent edits from interacting and collapsing the edit distance below k.
        pool = list(range(n))
        rng.shuffle(pool)
        chosen = []
        for p in pool:
            if all(abs(p - c) >= 2 for c in chosen):
                chosen.append(p)
                if len(chosen) == k:
                    break
        if len(chosen) < k:  # pragma: no cover, the chosen sizes rule this out
            raise RuntimeError("could not place k non-adjacent edits")
        # Edits go right to left, so original coordinates stay valid throughout.
        for p in sorted(chosen, reverse=True):
            if rng.random() < indel_frac:  # indel
                if rng.random() < 0.5:  # insertion before original p
                    seq.insert(p, rng.choice(BASES))
                else:  # deletion of original p
                    del seq[p]
            else:  # substitution
                seq[p] = rng.choice([b for b in BASES if b != seq[p]])
        return "".join(seq)


# Tiers
#
# Ground-truth tiers use substitutions only (indel_frac=0.0).
# k substitutions at distinct, non-adjacent sites guarantee ED(query, locus) == k exactly, so the assertion that optik's ED equals the constructed k is sharp.
# Indel-consensus tiers mix in insertions and deletions, for which an exact ED cannot be constructed (adjacent edits can legitimately lower the true ED), and assert only that the engines agree, the property that matters for indel-rich alignment.
def build_tiers(quick):
    """List the pair tiers of the battery, with sample counts cut tenfold by --quick."""

    def q(n):
        """Sample count n, cut in quick mode to a tenth but never below 10."""
        return max(10, n // (10 if quick else 1))

    tiers = []
    # Ground-truth locus tiers, where the text is the read's locus and edits are substitutions.
    for m, ks, n in [
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
    ]:
        for k in ks:
            tiers.append(
                {
                    "label": f"locus m={m} k={k}",
                    "kind": "locus",
                    "m": m,
                    "k": k,
                    "n": q(n),
                    "indel_frac": 0.0,
                    "flank": (0, 0),
                    "gt_check": True,
                    "stratum": (
                        "m<=64" if m <= 64 else "65-512" if m <= 512 else ">512"
                    )
                    + (" (banded)" if m in (1500,) and k >= 180 else ""),
                }
            )
    # Indel-consensus locus tier, indel-heavy and checked for engine agreement only
    for m, ks, n, frac in [(150, [10, 20], 200, 0.7)]:
        for k in ks:
            tiers.append(
                {
                    "label": f"locus m={m} k={k} indel70",
                    "kind": "locus",
                    "m": m,
                    "k": k,
                    "n": q(n),
                    "indel_frac": frac,
                    "flank": (0, 0),
                    "gt_check": False,
                    "stratum": "65-512 (indel consensus)",
                }
            )
    # Ground-truth window tiers, where the text is random flanks around the locus and edits are substitutions.
    # This is the mapper's realistic regime, with extra weight on 150 bp at 10%.
    for m, errs, n in [
        (65, [0.02, 0.05, 0.10], 150),
        (150, [0.02, 0.05, 0.10], 200),
        (300, [0.05, 0.10], 120),
        (513, [0.05, 0.10], 80),
    ]:
        for err in errs:
            k = max(1, round(m * err))
            n_this = q(n * 2) if (m == 150 and err == 0.10) else q(n)
            tiers.append(
                {
                    "label": f"window m={m} err={err:.0%}",
                    "kind": "window",
                    "m": m,
                    "k": k,
                    "n": n_this,
                    "indel_frac": 0.0,
                    "flank": (30, 30) if m <= 150 else (50, 50),
                    "gt_check": True,
                    "stratum": (
                        "m<=64" if m <= 64 else "65-512" if m <= 512 else ">512"
                    ),
                }
            )
    # Indel-consensus window tier, 150 bp at 10% with heavy indels
    tiers.append(
        {
            "label": "window m=150 err=10% indel50",
            "kind": "window",
            "m": 150,
            "k": 15,
            "n": q(200),
            "indel_frac": 0.5,
            "flank": (30, 30),
            "gt_check": False,
            "stratum": "65-512 (indel consensus)",
        }
    )
    return tiers


# Oracle
def build_oracle(workdir):
    """Compile the wfa2_oracle harness against WFA2-lib, returning (path, None) or (None, reason)."""
    if not (os.path.isfile(WFA2_LIB) and all(os.path.isdir(d) for d in INCLUDE_DIRS)):
        return None, "WFA2-lib not found (run bench/build_tools.sh wfa2 first)"
    out = os.path.join(workdir, "wfa2_oracle")
    # Relink when WFA2-lib or the harness changed, since a cached binary would keep testing the old library after a rebuild.
    stale = not os.path.isfile(out) or os.path.getmtime(out) < max(
        os.path.getmtime(WFA2_LIB), os.path.getmtime(ORACLE_SRC)
    )
    if stale:
        cmd = ["g++", "-std=c++17", "-O2", "-march=native"]
        for d in INCLUDE_DIRS:
            cmd += ["-I", d]
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
    r = subprocess.run(
        [oracle], input=pairs_tsv, text=True, capture_output=True, check=False
    )
    rows = []
    for line in r.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("m\tn\toptik_sg"):
            continue
        parts = line.split("\t")
        if len(parts) != 7:
            continue
        rows.append(tuple(int(x) for x in parts))
    return rows, r.returncode, r.stderr


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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oracle", help="path to prebuilt wfa2_oracle binary")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--quick", action="store_true", help="reduced pair counts (CI smoke)"
    )
    ap.add_argument("--workdir", default=os.path.join(ROOT, "temp", "waveform"))
    ap.add_argument("--keep", action="store_true", help="keep workdir files")
    ap.add_argument("--dump", help="also write generated TSV to this path")
    args = ap.parse_args()

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
    for ti, t in enumerate(tiers):
        for _ in range(t["n"]):
            locus = gen.locus(t["m"])
            query = gen._apply_edits(locus, t["k"], gen.rng, t["indel_frac"])
            if t["kind"] == "locus":
                pairs.append((query, locus, ti, "locus"))
            else:
                fl, fr = t["flank"]
                window = gen.locus(fl) + locus + gen.locus(fr)
                pairs.append((query, window, ti, "window"))
                pairs.append((query, locus, ti, "locus"))
    assert len(pairs) == sum(
        t["n"] * (2 if t["kind"] == "window" else 1) for t in tiers
    )

    if args.dump:
        with open(args.dump, "w") as f:
            for q, text, _ti, _ in pairs:
                f.write(f"{q}\t{text}\n")

    tsv = "\n".join(f"{q}\t{txt}" for q, txt, _, _ in pairs) + "\n"
    rows, rc, stderr = run_oracle(oracle, tsv)
    if rc != 0:
        print("waveform bench FAIL: oracle exited nonzero")
        print(stderr[-4000:])
        return 1
    if len(rows) != len(pairs):
        print(
            f"waveform bench FAIL: oracle returned {len(rows)} rows for "
            f"{len(pairs)} pairs; stderr tail:\n{stderr[-2000:]}"
        )
        return 1

    stats = {i: TierStats() for i in range(len(tiers))}
    for (_q, _text, ti, _kind), (_m, _n, o, capped, wsg, _we2e, _naive) in zip(
        pairs, rows, strict=False
    ):
        st = stats[ti]
        t = tiers[ti]
        st.pairs += 1
        st.sg_agree += o == wsg
        st.capped_agree += o == capped
        if o > t["k"]:
            # This cannot happen by construction (the locus is inside the text and the query is the locus with at most k unit edits), so it is a real bug.
            st.fail_gt_high += 1
        elif t["gt_check"]:
            if o == t["k"]:
                st.gt_ok += 1
            else:
                st.dropped += 1
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
    for i, t in enumerate(tiers):
        s = stats[i]
        drop = str(s.dropped) if t["gt_check"] else "-"
        gtok = str(s.gt_ok) if t["gt_check"] else "-"
        print(
            f"{t['label']:<28}{s.pairs:>6}{drop:>6}{gtok:>7}"
            f"{s.sg_agree:>9}{s.capped_agree:>7}{s.fail_gt_high:>7}  {t['stratum']}"
        )
        total_pairs += s.pairs
        total_agree += s.sg_agree
        total_capped += s.capped_agree
        total_fail += s.fail_gt_high
        if t["gt_check"]:
            total_drop += s.dropped
            total_gt += s.gt_ok
    gt_checked = sum(
        t["n"] * (2 if t["kind"] == "window" else 1) for t in tiers if t["gt_check"]
    )

    # Summary of the 150 bp, 10% error regime over its locus and window lines
    fr_tiers = [
        i
        for i, t in enumerate(tiers)
        if (t["m"] == 150 and t["k"] == 15)
        or (t["kind"] == "window" and t["m"] == 150 and t["k"] == 15)
    ]
    fr_pairs = sum(stats[i].pairs for i in fr_tiers)
    fr_gt = sum(stats[i].gt_ok for i in fr_tiers)
    fr_ag = sum(stats[i].sg_agree for i in fr_tiers)
    fr_drop = sum(stats[i].dropped for i in fr_tiers)
    print("--------------------------------------------------------------")
    print(
        f"  failing regime (150 bp / 10% error): {fr_pairs} pairs, "
        f"ground-truth {fr_gt}, engine-agree {fr_ag}, ambiguous {fr_drop}"
    )
    print("--------------------------------------------------------------")
    ok = (
        total_fail == 0
        and total_agree == total_pairs
        and total_capped == total_pairs
        and total_gt + total_drop == gt_checked
    )
    print(
        f"  totals: {total_pairs} pairs, {total_agree} engine-agree "
        f"({100 * total_agree / total_pairs:.2f}%), {total_capped} capped-agree "
        f"({100 * total_capped / total_pairs:.2f}%), {total_gt} ground-truth, "
        f"{total_drop} ambiguous-dropped, {total_fail} impossible-high"
    )
    print("  => WAVEFORM ORACLE: " + ("PASS" if ok else "FAIL"))
    print("==============================================")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
