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
/// Match score of the base-level alignment scores in [`MapqEvidence::align_scores`].
const MATCH_SCORE: f64 = 2.0;
/// BWA-MEM's phred-per-score-unit factor in its MAPQ cap (6.02 = 20 log10 2).
const BWA_PHRED_PER_SCORE: f64 = 6.02;
/// minimap2's chaining gap cost, GAP_COST_PER_K_BASE k g + GAP_COST_LOG_WEIGHT log2 g for a diagonal shift g.
const GAP_COST_PER_K_BASE: f64 = 0.01;
const GAP_COST_LOG_WEIGHT: f64 = 0.5;

/// An exact k-mer match at reference position `ref_pos` and read position `read_pos` in the strand's frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Anchor {
    pub ref_pos: u32,
    pub read_pos: u32,
}

/// A chain of colinear anchors on one reference strand.
#[derive(Clone, Debug)]
pub(crate) struct Chain {
    pub rid: u32,
    pub rev: bool,
    pub score: f64,
    /// Ascending in both r and q.
    pub anchors: Vec<Anchor>,
}

impl Chain {
    /// Read positions of the first and last anchor.
    pub(crate) fn q_span(&self) -> (u32, u32) {
        (
            self.anchors[0].read_pos,
            self.anchors[self.anchors.len() - 1].read_pos,
        )
    }
}

/// minimap2's gap cost of a diagonal shift `diag_shift`, 0.01 k g + 0.5 log2 g.
fn gap_cost(diag_shift: i64, kmer_len: f64) -> f64 {
    if diag_shift == 0 {
        0.0
    } else {
        GAP_COST_PER_K_BASE * kmer_len * diag_shift as f64
            + GAP_COST_LOG_WEIGHT * (diag_shift as f64).log2()
    }
}

/// Chains one (reference, strand) group of anchors into `out`.
///
/// `anchors` must be sorted by (r, q) and deduplicated.
/// Chains scoring below `min_score` or with fewer than `min_anchors` anchors are dropped.
pub(crate) fn chain_group(
    rid: u32,
    rev: bool,
    anchors: &[Anchor],
    kmer_len: usize,
    min_score: f64,
    min_anchors: usize,
    out: &mut Vec<Chain>,
) {
    let anchor_count = anchors.len();

    if anchor_count == 0 {
        return;
    }

    let kmer_len_f = kmer_len as f64;
    let mut scores = vec![0f64; anchor_count];
    let mut predecessors = vec![usize::MAX; anchor_count];

    for i in 0..anchor_count {
        let (ref_i, read_i) = (anchors[i].ref_pos as i64, anchors[i].read_pos as i64);
        // An anchor on its own scores k. Linking to a predecessor has to beat that fresh start.
        let mut best = kmer_len_f;
        let mut best_pred = usize::MAX;
        let mut best_g = i64::MAX;
        let window_start = i.saturating_sub(LOOKBACK_MAX);

        for (n, j) in (window_start..i).rev().enumerate() {
            // After LOOKBACK candidates a near-diagonal predecessor will not be beaten, so the scan stops early.
            if n >= LOOKBACK && best_g <= GAP_SETTLED {
                break;
            }

            let ref_dist = ref_i - anchors[j].ref_pos as i64;

            if ref_dist > MAX_DIST {
                break;
            } // sorted by r, so older ones are further away

            let read_dist = read_i - anchors[j].read_pos as i64;

            // A predecessor must lie strictly before on both sequences.
            if ref_dist <= 0 || read_dist <= 0 || read_dist > MAX_DIST {
                continue;
            }

            let diag_shift = (ref_dist - read_dist).abs();

            if diag_shift > BANDWIDTH {
                continue;
            }

            // As in minimap2, the predecessor adds at most k new matching bases and pays for the diagonal shift between the two anchors.
            let score = scores[j] + (read_dist.min(ref_dist).min(kmer_len as i64)) as f64
                - gap_cost(diag_shift, kmer_len_f);

            if score > best {
                best = score;
                best_pred = j;
                best_g = diag_shift;
            }
        }

        scores[i] = best;
        predecessors[i] = best_pred;
    }

    // Chains are traced back from the best-scoring ends first, and ties fall to the lower index so the output is deterministic.
    let mut order: Vec<usize> = (0..anchor_count).collect();

    order.sort_unstable_by(|&lhs, &rhs| scores[rhs].total_cmp(&scores[lhs]).then(lhs.cmp(&rhs)));
    let mut used = vec![false; anchor_count];

    for &end in &order {
        if used[end] {
            continue;
        }

        let mut idx = Vec::new();
        let mut node = end;
        let mut base = 0.0;

        // A traceback that reaches an anchor of an earlier chain stops there, and that anchor's score is subtracted, so the anchors are not shared and each chain keeps only its own part of the score.
        loop {
            if used[node] {
                base = scores[node];
                break;
            }

            used[node] = true;
            idx.push(node);

            if predecessors[node] == usize::MAX {
                break;
            }

            node = predecessors[node];
        }

        let score = scores[end] - base;

        if score < min_score || idx.len() < min_anchors {
            continue;
        }

        idx.reverse();
        out.push(Chain {
            rid,
            rev,
            score,
            anchors: idx.iter().map(|&anchor_idx| anchors[anchor_idx]).collect(),
        });
    }
}

/// minimap2's MAPQ from the chain score `chain_score`, the best competing (query-overlapping) chain score `competitor_score` and the anchor count.
pub(crate) fn mapq(chain_score: f64, competitor_score: f64, n_anchors: usize) -> u8 {
    if chain_score <= 0.0 {
        return 0;
    }

    let mapq = MAPQ_COEF
        * (1.0 - competitor_score / chain_score)
        * (n_anchors as f64 / FULL_WEIGHT_ANCHORS).min(1.0)
        * chain_score.ln();

    mapq.clamp(0.0, MAPQ_MAX) as u8
}

/// Evidence for one long-read alignment's MAPQ.
pub(crate) struct MapqEvidence {
    /// Its chain score and the best competing, query-overlapping chain's score.
    pub chain_score: f64,
    pub competitor_score: f64,
    pub n_anchors: usize,
    /// Sum of the read's kept chain scores, primary and supplementary, and the read bases that only repetitive, dropped minimizers cover.
    pub sum_kept: f64,
    pub rep_len: u32,
    /// Aligned identity, (columns - edits) / columns.
    pub identity: f64,
    /// Base-level alignment scores (match +2, mismatch -4, gap -(4 + 2 len)) of this alignment and of the best competing copy over the same read span, when competing copies were aligned.
    pub align_scores: Option<(f64, f64)>,
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
pub(crate) fn mapq_long(evidence: &MapqEvidence) -> u8 {
    if evidence.chain_score <= 0.0 {
        return 0;
    }

    let competitor_score = evidence.competitor_score.max(MIN_CHAIN_SCORE);
    let uniq = if evidence.sum_kept > 0.0 {
        evidence.sum_kept / (evidence.sum_kept + evidence.rep_len as f64)
    } else {
        0.0
    };
    let pen_s1 = (evidence.chain_score / FULL_WEIGHT_SCORE).min(1.0) * uniq;
    let pen_cm = (evidence.n_anchors as f64 / FULL_WEIGHT_ANCHORS).min(1.0);
    let pen = pen_s1.min(pen_cm);
    let chain_ratio = (competitor_score / evidence.chain_score).min(1.0);
    let mapq = match evidence.align_scores {
        Some((align_score, competitor_align_score)) if align_score > MATCH_SCORE => {
            let combined_ratio =
                (competitor_align_score.max(0.0) / align_score).min(1.0) * chain_ratio;
            let mapq_estimate = MAPQ_COEF
                * evidence.identity
                * pen
                * (1.0 - combined_ratio * combined_ratio)
                * (align_score / MATCH_SCORE).ln();

            mapq_estimate.min(
                BWA_PHRED_PER_SCORE
                    * evidence.identity
                    * evidence.identity
                    * (align_score - competitor_align_score).max(0.0)
                    / MATCH_SCORE,
            )
        }
        _ => MAPQ_COEF * evidence.identity * pen * (1.0 - chain_ratio) * evidence.chain_score.ln(),
    };

    mapq.clamp(0.0, MAPQ_MAX) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Anchors on one diagonal, 20 bases apart, chain into a single chain that scores k per anchor.
    #[test]
    fn colinear_anchors_form_one_chain() {
        let kmer_len = 15;
        let anchors: Vec<Anchor> = (0..40)
            .map(|i| Anchor {
                ref_pos: 1000 + i * 20,
                read_pos: i * 20,
            })
            .collect();
        let mut out = Vec::new();

        chain_group(
            0,
            false,
            &anchors,
            kmer_len,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 40);
        assert!((out[0].score - (15.0 * 40.0)).abs() < 1e-9);
    }

    /// An indel inside the bandwidth costs score but keeps one chain, while anchors at a distant locus form a chain of their own.
    #[test]
    fn indel_costs_and_separate_loci_split() {
        let kmer_len = 15;
        let mut anchors: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                ref_pos: 100 + i * 20,
                read_pos: i * 20,
            })
            .collect();

        // a 30 bp deletion in the read keeps one chain (gap 30 < bandwidth)
        anchors.extend((10..20).map(|i| Anchor {
            ref_pos: 100 + i * 20 + 30,
            read_pos: i * 20,
        }));
        // A far locus (> MAX_DIST away) is a separate chain.
        anchors.extend((0..10).map(|i| Anchor {
            ref_pos: 900_000 + i * 20,
            read_pos: 400 + i * 20,
        }));
        anchors.sort();
        let mut out = Vec::new();

        chain_group(
            0,
            false,
            &anchors,
            kmer_len,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].anchors.len(), 20);
        assert!(out[0].score < 15.0 * 20.0 && out[0].score > 15.0 * 18.0);
        assert_eq!(out[1].anchors.len(), 10);
    }

    /// Inside a tandem array the best chain follows one copy, even when every nearby anchor lies on another copy's diagonal.
    #[test]
    fn tandem_array_chain_stays_on_its_diagonal() {
        // A 3 kb read inside a 200-copy tandem array of a 25 bp unit, with one minimizer every 5 bp.
        // Each read k-mer hits all copies, so about 120 anchors share each reference position and the 50 nearest anchors of an anchor all lie on other diagonals.
        // The chain must follow one diagonal, without a unit-sized jump.
        let kmer_len = 15;
        let (unit, copies) = (25u32, 200u32);
        let (array_start, array_end) = (10_000u32, 10_000 + unit * copies);
        let mut anchors = Vec::new();

        for q in (0..3000u32).step_by(5) {
            for d in 0..copies {
                // Read position q lies on copy diagonal d at reference lo + q + d * unit - 3000 while it is inside the array.
                let ref_pos = (array_start + q + d * unit).wrapping_sub(3000);

                if ref_pos >= array_start && ref_pos < array_end {
                    anchors.push(Anchor {
                        ref_pos,
                        read_pos: q,
                    });
                }
            }
        }

        anchors.sort_unstable();
        anchors.dedup();
        let mut out = Vec::new();

        chain_group(
            0,
            false,
            &anchors,
            kmer_len,
            MIN_CHAIN_SCORE,
            MIN_CHAIN_ANCHORS,
            &mut out,
        );
        let best = out
            .iter()
            .max_by(|lhs, rhs| lhs.score.total_cmp(&rhs.score))
            .unwrap();
        let first_diag = best.anchors[0].ref_pos as i64 - best.anchors[0].read_pos as i64;

        assert!(
            best.anchors
                .iter()
                .all(|anchor| anchor.ref_pos as i64 - anchor.read_pos as i64 == first_diag),
            "chain jumped between tandem copies"
        );
        assert!(best.anchors.len() >= 300, "{} anchors", best.anchors.len());
    }

    /// The long-read MAPQ rewards strong unique chains and discounts ties, weak chains, repeat-covered reads and close aligned competitors.
    #[test]
    fn mapq_long_behaviour() {
        let evidence =
            |chain_score: f64, competitor_score: f64, n_anchors: usize, rep: u32, align_scores| {
                MapqEvidence {
                    chain_score,
                    competitor_score,
                    n_anchors,
                    sum_kept: chain_score,
                    rep_len: rep,
                    identity: 0.95,
                    align_scores,
                }
            };

        // A strong unique chain is certain.
        assert_eq!(mapq_long(&evidence(3000.0, 0.0, 300, 0, None)), 60);
        // A tie is 0.
        assert_eq!(mapq_long(&evidence(3000.0, 3000.0, 300, 0, None)), 0);
        // A weak chain just above the threshold, for which the old formula said 60.
        // With s2 floored at the threshold and the weak-chain factor it is about 20.
        assert_eq!(mapq(60.0, 0.0, 4), 60);
        assert!(mapq_long(&evidence(60.0, 0.0, 4, 0, None)) < 30);
        // A read mostly covered by dropped repeat seeds is discounted.
        assert!(mapq_long(&evidence(400.0, 0.0, 30, 14_000, None)) < 10);
        // Aligned competitors one mismatch apart, a score gap of 6, cap the MAPQ.
        let mapq = mapq_long(&evidence(3000.0, 0.0, 300, 0, Some((5000.0, 4994.0))));

        assert!(mapq > 10 && mapq < 20, "{mapq}");
        // Near-tied chains (s2/s1 = 0.97) whose alignments differ by 190 score units.
        // The chain ratio alone gives about 10, and combined about 25.
        let alone = mapq_long(&evidence(14000.0, 13580.0, 1700, 0, None));
        let both = mapq_long(&evidence(
            14000.0,
            13580.0,
            1700,
            0,
            Some((29000.0, 28810.0)),
        ));

        assert!(alone < 15 && both > 20, "{alone} {both}");
        // Identical copies with equal alignment scores stay a tie.
        assert_eq!(
            mapq_long(&evidence(
                14000.0,
                14000.0,
                1700,
                0,
                Some((29000.0, 29000.0))
            )),
            0
        );
    }

    /// minimap2's chain MAPQ gives 0 for a tie, 60 for a unique chain and an intermediate value for a close competitor.
    #[test]
    fn mapq_behaviour() {
        assert_eq!(mapq(300.0, 300.0, 30), 0);
        assert!(mapq(300.0, 0.0, 30) == 60);
        // A close competitor gives 40 * (1 - 280/300) * ln(300), about 15.
        assert!(mapq(300.0, 280.0, 30) > 5 && mapq(300.0, 280.0, 30) < 30);
    }
}
