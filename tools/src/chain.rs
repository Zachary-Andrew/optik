//! Anchor chaining and long-read MAPQ.
//!
//! The chaining DP and gap cost are those of minimap2 (Li 2018), implemented independently.
//!
//! Anchors are exact k-mer matches (reference pos r, read pos q in the strand's frame) of one (reference, strand) group.
//! With anchors sorted by (r, q), the chain score is
//!   f(i) = max(k, max_j f(j) + min(qd, rd, k) - gap_cost(|rd - qd|))
//! over the previous `LOOKBACK` anchors j (up to `LOOKBACK_MAX` while no near-diagonal predecessor has been found) with 0 < qd, rd <= MAX_DIST and
//! |rd - qd| <= BANDWIDTH, where gap_cost(g) = 0.01*k*g + 0.5*log2(g).
//! The bounded lookback keeps the DP O(n * LOOKBACK) and its working set in L1 (anchors are 8-byte records scanned backwards contiguously).
//!
//! Chains are extracted best-first.
//! A backtrack stops at an anchor already claimed by a better chain and is scored relative to it, so every anchor belongs to one chain.

/// Predecessor anchors examined per anchor in the usual case.
pub(crate) const LOOKBACK: usize = 50;
/// The lookback extends to LOOKBACK_MAX while the best predecessor found so far needs a gap of more than GAP_SETTLED.
/// Inside a tandem array every read k-mer hits every copy, so the 50 nearest anchors can all lie on other diagonals and the same-diagonal predecessor is further back.
/// Stopping at 50 there made the chain jump one unit sideways, and the alignment then carried unit-sized false indels (0.7% of hg38 chr7 HiFi reads had a spurious indel of 50 bp or more, and 42% on a synthetic repeat genome).
pub(crate) const LOOKBACK_MAX: usize = 2000;
const GAP_SETTLED: i64 = 8;
/// Largest reference or read distance between chained anchors.
pub(crate) const MAX_DIST: i64 = 5000;
/// Largest diagonal shift |rd - qd| between chained anchors.
pub(crate) const BANDWIDTH: i64 = 500;
pub(crate) const MIN_CHAIN_SCORE: f64 = 40.0;
pub(crate) const MIN_CHAIN_ANCHORS: usize = 3;

/// minimap2's MAPQ coefficient and cap.
const MAPQ_COEF: f64 = 40.0;
const MAPQ_MAX: f64 = 60.0;
/// Anchor count at which the anchor-count penalty vanishes, 10 as in minimap2.
const FULL_WEIGHT_ANCHORS: f64 = 10.0;
/// Chain score at which the weak-chain penalty vanishes.
const FULL_WEIGHT_SCORE: f64 = 100.0;
/// Match score of the base-level alignment scores in [`MapqEvidence::dp`].
const MATCH_SCORE: f64 = 2.0;
/// BWA-MEM's phred-per-score-unit factor in its MAPQ cap (6.02 = 20 log10 2).
const BWA_PHRED_PER_SCORE: f64 = 6.02;

/// An exact k-mer match at reference position `r` and read position `q` in the strand's frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Anchor {
    pub r: u32,
    pub q: u32,
}

/// A chain of colinear anchors on one reference strand.
#[derive(Clone, Debug)]
pub(crate) struct Chain {
    pub rid: u32,
    pub rev: bool,
    pub score: f64,
    /// ascending in both r and q
    pub anchors: Vec<Anchor>,
}

impl Chain {
    /// Read positions of the first and last anchor.
    pub(crate) fn q_span(&self) -> (u32, u32) {
        (self.anchors[0].q, self.anchors[self.anchors.len() - 1].q)
    }
}

/// minimap2's gap cost of a diagonal shift `g`, 0.01 k g + 0.5 log2 g.
fn gap_cost(g: i64, k: f64) -> f64 {
    if g == 0 {
        0.0
    } else {
        0.01 * k * g as f64 + 0.5 * (g as f64).log2()
    }
}

/// Chains one (reference, strand) group of anchors into `out`.
///
/// `a` must be sorted by (r, q) and deduplicated.
/// Chains scoring below `min_score` or with fewer than `min_anchors` anchors are dropped.
pub(crate) fn chain_group(
    rid: u32,
    rev: bool,
    a: &[Anchor],
    k: usize,
    min_score: f64,
    min_anchors: usize,
    out: &mut Vec<Chain>,
) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let kf = k as f64;
    let mut f = vec![0f64; n];
    let mut p = vec![usize::MAX; n];
    for i in 0..n {
        let (ri, qi) = (a[i].r as i64, a[i].q as i64);
        let mut best = kf;
        let mut bp = usize::MAX;
        let mut best_g = i64::MAX;
        let lo = i.saturating_sub(LOOKBACK_MAX);
        for (n, j) in (lo..i).rev().enumerate() {
            if n >= LOOKBACK && best_g <= GAP_SETTLED {
                break;
            }
            let rd = ri - a[j].r as i64;
            if rd > MAX_DIST {
                break;
            } // sorted by r, so older ones are further away
            let qd = qi - a[j].q as i64;
            if rd <= 0 || qd <= 0 || qd > MAX_DIST {
                continue;
            }
            let g = (rd - qd).abs();
            if g > BANDWIDTH {
                continue;
            }
            let sc = f[j] + (qd.min(rd).min(k as i64)) as f64 - gap_cost(g, kf);
            if sc > best {
                best = sc;
                bp = j;
                best_g = g;
            }
        }
        f[i] = best;
        p[i] = bp;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&x, &y| f[y].total_cmp(&f[x]).then(x.cmp(&y)));
    let mut used = vec![false; n];
    for &end in &order {
        if used[end] {
            continue;
        }
        let mut idx = Vec::new();
        let mut i = end;
        let mut base = 0.0;
        loop {
            if used[i] {
                base = f[i];
                break;
            }
            used[i] = true;
            idx.push(i);
            if p[i] == usize::MAX {
                break;
            }
            i = p[i];
        }
        let score = f[end] - base;
        if score < min_score || idx.len() < min_anchors {
            continue;
        }
        idx.reverse();
        out.push(Chain {
            rid,
            rev,
            score,
            anchors: idx.iter().map(|&t| a[t]).collect(),
        });
    }
}

/// minimap2's MAPQ from the chain score `s1`, the best competing (query-overlapping) chain score `s2` and the anchor count.
pub(crate) fn mapq(s1: f64, s2: f64, n_anchors: usize) -> u8 {
    if s1 <= 0.0 {
        return 0;
    }
    let q =
        MAPQ_COEF * (1.0 - s2 / s1) * (n_anchors as f64 / FULL_WEIGHT_ANCHORS).min(1.0) * s1.ln();
    q.clamp(0.0, MAPQ_MAX) as u8
}

/// Evidence for one long-read alignment's MAPQ.
pub(crate) struct MapqEvidence {
    /// its chain score and the best competing (query-overlapping) chain's
    pub s1: f64,
    pub s2: f64,
    pub n_anchors: usize,
    /// sum of the read's kept (primary + supplementary) chain scores, and the read bases that only repetitive (dropped) minimizers cover
    pub sum_kept: f64,
    pub rep_len: u32,
    /// aligned identity, (columns - edits) / columns
    pub identity: f64,
    /// base-level alignment scores (match +2, mismatch -4, gap -(4 + 2 len)) of this alignment and of the best competing copy over the same read span, when competing copies were aligned
    pub dp: Option<(f64, f64)>,
}

/// Long-read MAPQ.
///
/// Concept from minimap2's MAPQ (Li 2018, 2021), formulated over optik's chains and alignments, with BWA-MEM's alignment-score cap (Li 2013).
///   pen  = min(min(1, s1/100) * uniq, min(1, n_anchors/10)),
///   uniq = sum_kept / (sum_kept + rep_len),
/// s2 floored at MIN_CHAIN_SCORE (a competitor below the chaining threshold was never extracted, not absent).
/// Without aligned competitors
///   MAPQ = 40 * identity * pen * (1 - s2/s1) * ln(s1);
/// with them, the chain and alignment ratios are combined,
///   x = (dp2/dp1) * (s2/s1),  MAPQ = 40 * identity * pen * (1 - x^2) * ln(dp1/2),
/// capped BWA-style by the alignment-score gap, 6.02 * identity^2 * (dp1 - dp2) / 2 with match score 2.
/// A weak chain cannot claim certainty, and a read lying mostly in dropped repeat seeds, with uniq << 1, is discounted, while a copy the alignment separates from its near-tied chain twin is not.
pub(crate) fn mapq_long(ev: &MapqEvidence) -> u8 {
    if ev.s1 <= 0.0 {
        return 0;
    }
    let s2 = ev.s2.max(MIN_CHAIN_SCORE);
    let uniq = if ev.sum_kept > 0.0 {
        ev.sum_kept / (ev.sum_kept + ev.rep_len as f64)
    } else {
        0.0
    };
    let pen_s1 = (ev.s1 / FULL_WEIGHT_SCORE).min(1.0) * uniq;
    let pen_cm = (ev.n_anchors as f64 / FULL_WEIGHT_ANCHORS).min(1.0);
    let pen = pen_s1.min(pen_cm);
    let xs = (s2 / ev.s1).min(1.0);
    let q = match ev.dp {
        Some((d1, d2)) if d1 > MATCH_SCORE => {
            let x = (d2.max(0.0) / d1).min(1.0) * xs;
            let q = MAPQ_COEF * ev.identity * pen * (1.0 - x * x) * (d1 / MATCH_SCORE).ln();
            q.min(
                BWA_PHRED_PER_SCORE * ev.identity * ev.identity * (d1 - d2).max(0.0) / MATCH_SCORE,
            )
        }
        _ => MAPQ_COEF * ev.identity * pen * (1.0 - xs) * ev.s1.ln(),
    };
    q.clamp(0.0, MAPQ_MAX) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colinear_anchors_form_one_chain() {
        let k = 15;
        let a: Vec<Anchor> = (0..40)
            .map(|i| Anchor {
                r: 1000 + i * 20,
                q: i * 20,
            })
            .collect();
        let mut out = Vec::new();
        chain_group(
            0,
            false,
            &a,
            k,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 40);
        assert!((out[0].score - (15.0 * 40.0)).abs() < 1e-9);
    }

    #[test]
    fn indel_costs_and_separate_loci_split() {
        let k = 15;
        let mut a: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                r: 100 + i * 20,
                q: i * 20,
            })
            .collect();
        // a 30 bp deletion in the read keeps one chain (gap 30 < bandwidth)
        a.extend((10..20).map(|i| Anchor {
            r: 100 + i * 20 + 30,
            q: i * 20,
        }));
        // a far locus (> MAX_DIST away) is a separate chain
        a.extend((0..10).map(|i| Anchor {
            r: 900_000 + i * 20,
            q: 400 + i * 20,
        }));
        a.sort();
        let mut out = Vec::new();
        chain_group(
            0,
            false,
            &a,
            k,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].anchors.len(), 20);
        assert!(out[0].score < 15.0 * 20.0 && out[0].score > 15.0 * 18.0);
        assert_eq!(out[1].anchors.len(), 10);
    }

    #[test]
    fn tandem_array_chain_stays_on_its_diagonal() {
        // A 3 kb read inside a 200-copy tandem array of a 25 bp unit, with one minimizer every 5 bp.
        // Each read k-mer hits all copies, so about 120 anchors share each reference position and the 50 nearest anchors of an anchor all lie on other diagonals.
        // The chain must follow one diagonal, without a unit-sized jump.
        let k = 15;
        let (unit, copies) = (25u32, 200u32);
        let (lo, hi) = (10_000u32, 10_000 + unit * copies);
        let mut a = Vec::new();
        for q in (0..3000u32).step_by(5) {
            for d in 0..copies {
                // read position q lies on copy diagonal d at reference lo + q + d * unit - 3000 while it is inside the array
                let r = (lo + q + d * unit).wrapping_sub(3000);
                if r >= lo && r < hi {
                    a.push(Anchor { r, q });
                }
            }
        }
        a.sort_unstable();
        a.dedup();
        let mut out = Vec::new();
        chain_group(
            0,
            false,
            &a,
            k,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        let best = out
            .iter()
            .max_by(|x, y| x.score.total_cmp(&y.score))
            .unwrap();
        let d0 = best.anchors[0].r as i64 - best.anchors[0].q as i64;
        assert!(
            best.anchors.iter().all(|x| x.r as i64 - x.q as i64 == d0),
            "chain jumped between tandem copies"
        );
        assert!(best.anchors.len() >= 300, "{} anchors", best.anchors.len());
    }

    #[test]
    fn mapq_long_behaviour() {
        let ev = |s1: f64, s2: f64, n: usize, rep: u32, dp| MapqEvidence {
            s1,
            s2,
            n_anchors: n,
            sum_kept: s1,
            rep_len: rep,
            identity: 0.95,
            dp,
        };
        // a strong unique chain is certain
        assert_eq!(mapq_long(&ev(3000.0, 0.0, 300, 0, None)), 60);
        // a tie is 0
        assert_eq!(mapq_long(&ev(3000.0, 3000.0, 300, 0, None)), 0);
        // A weak chain just above the threshold, for which the old formula said 60.
        // With s2 floored at the threshold and the weak-chain factor it is about 20.
        assert_eq!(mapq(60.0, 0.0, 4), 60);
        assert!(mapq_long(&ev(60.0, 0.0, 4, 0, None)) < 30);
        // a read mostly covered by dropped repeat seeds is discounted
        assert!(mapq_long(&ev(400.0, 0.0, 30, 14_000, None)) < 10);
        // aligned competitors one mismatch apart, a score gap of 6, cap MAPQ
        let q = mapq_long(&ev(3000.0, 0.0, 300, 0, Some((5000.0, 4994.0))));
        assert!(q > 10 && q < 20, "{q}");
        // Near-tied chains (s2/s1 = 0.97) whose alignments differ by 190 score units.
        // The chain ratio alone gives about 10, and combined about 25.
        let alone = mapq_long(&ev(14000.0, 13580.0, 1700, 0, None));
        let both = mapq_long(&ev(14000.0, 13580.0, 1700, 0, Some((29000.0, 28810.0))));
        assert!(alone < 15 && both > 20, "{alone} {both}");
        // identical copies with equal alignment scores stay a tie
        assert_eq!(
            mapq_long(&ev(14000.0, 14000.0, 1700, 0, Some((29000.0, 29000.0)))),
            0
        );
    }

    #[test]
    fn mapq_behaviour() {
        assert_eq!(mapq(300.0, 300.0, 30), 0);
        assert!(mapq(300.0, 0.0, 30) == 60);
        // A close competitor gives 40 * (1 - 280/300) * ln(300), about 15.
        assert!(mapq(300.0, 280.0, 30) > 5 && mapq(300.0, 280.0, 30) < 30);
    }
}
