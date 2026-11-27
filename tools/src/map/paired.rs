//! Paired-end mapping fits an insert-size model, rescues mates inside the insert window and selects both placements jointly.

use super::*;

// Paired-end reads, map -q R1 -q R2
//
// Each mate is first placed on its own exactly as a single-end read, with its best placement and its seeded competitors.
// The library's proper template range comes from the first batch's confident FR pairs, as quartiles +- 3 IQR the way BWA estimates it.
// A pair is then resolved as BWA-MEM resolves pairs, in likelihood form.
//   - The candidates are both mates' best placements and their competitors within 3 edits.
//     If no combination is a proper pair, mate rescue also aligns either mate inside the window its insert size allows next to the partner's best candidates.
//   - Every combination is weighted r1^(ed1 - best1) * r2^(ed2 - best2), with r the edit_ratio, times discordant_prior, worth 3.4 edits, when it is not a proper pair.
//     The most likely combination is chosen, and ties are broken by a hash of the pair's sequence.
//   - A mate's MAPQ is -10 log10 of the weight of the combinations that place it elsewhere, capped like the unpaired MAPQ by hidden-locus odds and seed evidence, so a mate in a repeat whose partner is unique gets the confidence the pair gives it.
//   - The SAM pair fields are 0x1, 0x2 for a proper pair, 0x8 and 0x20 for an unmapped or reverse mate, 0x40 / 0x80, RNEXT / PNEXT, a signed TLEN over the outer coordinates, positive on the leftmost mate, and the mate's CIGAR and MAPQ in MC:Z and MQ:i.
//     An unmapped mate is placed at its partner, as SAM spec 1.4 describes.
// Long mates of CHAIN_MIN_LEN or more keep their own placement, and their pair fields are filled the same way.

/// Prior weight of a placement the mate does not support, a discordant pair, in edits at the read's own error rate, r^3.4 for edit ratio r.
///
/// BWA-MEM's unpaired penalty is 17 score points, i.e. 3.4 mismatches at 5 points each.
/// A fixed 1e-3 was about 2 edits at 1% error.
/// Real HG002 mates whose partner was placed uniquely then went, as a discordant pair, to a paralog 2 mismatches closer through the sample's own SNPs, which the paralog shares, while GIAB confirms the proper placement.
/// Simulated reads have no variants, so they never show this.
pub(super) const UNPAIRED_EDITS: f64 = 3.4;

/// Pairs whose mates both reach this unpaired MAPQ (and have no repetitive seeds) train the insert model.
pub(super) const CONFIDENT_MAPQ: i32 = 30;

/// Template lengths at or beyond this are not a proper FR pair and stay out of the insert model.
pub(super) const MAX_TEMPLATE_LEN: i64 = 100_000;

/// Longest proper template while no insert model has been estimated yet.
pub(super) const FALLBACK_MAX_TLEN: i64 = 1000;

/// Proper template lengths used until the first batch yields an insert model.
pub(super) const FALLBACK_INSERT_MODEL: InsertModel = InsertModel {
    min_tlen: MIN_TEMPLATE_LEN,
    max_tlen: FALLBACK_MAX_TLEN,
};

/// The insert model needs at least this many confident pairs.
pub(super) const MIN_INSERT_PAIRS: usize = 30;

/// The proper range is [Q1 - IQR_FENCE IQR, Q3 + IQR_FENCE IQR] (Tukey's far-out fences), at least MIN_TEMPLATE_LEN.
pub(super) const LOWER_QUARTILE: f64 = 0.25;

pub(super) const UPPER_QUARTILE: f64 = 0.75;

pub(super) const IQR_FENCE: i64 = 3;

pub(super) const MIN_TEMPLATE_LEN: i64 = 1;

pub(super) fn discordant_prior(ratio: f64) -> f64 {
    ratio.powf(UNPAIRED_EDITS)
}

/// Proper template lengths [lo, hi] for FR orientation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InsertModel {
    pub min_tlen: i64,
    pub max_tlen: i64,
}

/// The insert model, [max(1, Q1 - 3 IQR), Q3 + 3 IQR] from template lengths of confidently placed pairs, or None below 30 pairs.
pub(crate) fn insert_model(mut tlens: Vec<i64>) -> Option<InsertModel> {
    if tlens.len() < MIN_INSERT_PAIRS {
        return None;
    }

    tlens.sort_unstable();
    let quantile = |frac: f64| tlens[((tlens.len() - 1) as f64 * frac).round() as usize];
    let (first_quartile, third_quartile) = (quantile(LOWER_QUARTILE), quantile(UPPER_QUARTILE));
    let iqr = (third_quartile - first_quartile).max(1);

    Some(InsertModel {
        min_tlen: (first_quartile - IQR_FENCE * iqr).max(MIN_TEMPLATE_LEN),
        max_tlen: third_quartile + IQR_FENCE * iqr,
    })
}

/// One mate after the independent pass.
pub(super) enum MateAln {
    Short(Option<ShortHit>),
    Lines(Vec<String>),
}

/// Places one mate on its own.
/// Long mates and short mates rescued on the chaining path keep their records, while short mates keep their hit for pairing.
pub(super) fn map_mate(
    read_name: &str,
    seq: &[u8],
    max_ed: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> MateAln {
    if seq.len() >= CHAIN_MIN_LEN {
        return MateAln::Lines(map_one_read_chain(
            read_name, seq, ctx, ref_names, LONG_CHAIN,
        ));
    }

    let hit = align_read_impl(seq, read_max_ed(max_ed, seq.len()), ctx);

    if hit.is_none() && seq.len() >= RESCUE_MIN_KMERS * ctx.kmer_len {
        let lines = map_one_read_chain(read_name, seq, ctx, ref_names, RESCUE_CHAIN);

        if !lines.is_empty() && !is_unmapped(&lines[0]) {
            return MateAln::Lines(lines);
        }
    }

    MateAln::Short(hit)
}

/// A placement in the form pairing uses, as reference, strand and [start, end).
#[derive(Clone, Copy, Debug)]
pub(super) struct Place {
    pub(super) rid: usize,
    pub(super) rev: bool,
    pub(super) start: i64,
    pub(super) end: i64,
}

/// Template length of an FR pair, from the forward mate's start to the reverse mate's end, or None if the two are not FR on one reference.
pub(super) fn fr_tlen(mate_a: Place, mate_b: Place) -> Option<i64> {
    if mate_a.rid != mate_b.rid || mate_a.rev == mate_b.rev {
        return None;
    }

    let (fwd_mate, rev_mate) = if mate_a.rev {
        (mate_b, mate_a)
    } else {
        (mate_a, mate_b)
    };

    if fwd_mate.start > rev_mate.end {
        return None;
    }

    Some(rev_mate.end - fwd_mate.start)
}

/// Whether two placements form a proper FR pair under the insert model.
pub(super) fn proper(mate_a: Place, mate_b: Place, model: &InsertModel) -> bool {
    fr_tlen(mate_a, mate_b).is_some_and(|tlen| tlen >= model.min_tlen && tlen <= model.max_tlen)
}

/// The placement of a short hit of a `len`-base read.
pub(super) fn hit_place(hit: &ShortHit, len: usize) -> Place {
    Place {
        rid: hit.rid,
        rev: hit.rev,
        start: hit.end - len as i64,
        end: hit.end,
    }
}

/// Primary-record fields pairing reads back from a SAM line.
pub(super) struct Prim {
    pub(super) mapped: bool,
    pub(super) rname: String,
    pub(super) pos: i64,
    pub(super) rev: bool,
    pub(super) mapq: i32,
    pub(super) cigar: String,
    pub(super) end: i64,
}

/// Reads the pairing fields back from a primary SAM record.
pub(super) fn parse_prim(line: &str) -> Prim {
    // Splitting stops at RNEXT, so the long SEQ and QUAL columns are never cut apart. Missing fields read as unmapped.
    let fields: Vec<&str> = line.splitn(SAM_COL_RNEXT + 1, '\t').collect();
    let flag: u16 = fields
        .get(SAM_COL_FLAG)
        .and_then(|field| field.parse().ok())
        .unwrap_or(SAM_FLAG_UNMAPPED);
    let pos: i64 = fields
        .get(SAM_COL_POS)
        .and_then(|field| field.parse().ok())
        .unwrap_or(0);
    let cigar = fields
        .get(SAM_COL_CIGAR)
        .copied()
        .unwrap_or("*")
        .to_string();

    Prim {
        mapped: flag & SAM_FLAG_UNMAPPED == 0,
        rname: fields
            .get(SAM_COL_RNAME)
            .copied()
            .unwrap_or("*")
            .to_string(),
        pos,
        rev: flag & SAM_FLAG_REVERSE != 0,
        mapq: fields
            .get(SAM_COL_MAPQ)
            .and_then(|field| field.parse().ok())
            .unwrap_or(0),
        end: pos - 1 + ref_span(&cigar),
        cigar,
    }
}

/// Reference bases a CIGAR covers (M, D, N, =, X).
pub(super) fn ref_span(cig: &str) -> i64 {
    let (mut run_len, mut span) = (0i64, 0i64);

    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        if matches!(byte, b'M' | b'D' | b'N' | b'=' | b'X') {
            span += run_len;
        }

        run_len = 0;
    }

    span
}

/// Aligns `seq`, on strand `rev`, inside the reference window where the partner's insert size puts it.
/// Returns the placement as a ShortHit without competitors yet, or None.
pub(super) fn rescue_in_window(
    seq: &[u8],
    rev_seq: &[u8],
    partner: Place,
    model: &InsertModel,
    max_ed: i32,
    ctx: &AlignCtx<'_>,
) -> Option<ShortHit> {
    let read_len = seq.len() as i64;
    let rlen = ctx.ref_seqs[partner.rid].len() as i64;
    let rev = !partner.rev;
    // With the partner forward the mate is reverse and ends in [p.start + lo, p.start + hi].
    // With the partner reverse the mate is forward and starts in [p.end - hi, p.end - lo].
    let (win_start, win_end) = if !partner.rev {
        (
            partner.start + model.min_tlen - read_len,
            partner.start + model.max_tlen,
        )
    } else {
        (
            partner.end - model.max_tlen,
            partner.end - model.min_tlen + read_len,
        )
    };
    let (win_start, win_end) = (win_start.max(0), win_end.min(rlen));

    if win_end - win_start < read_len {
        return None;
    }

    let oriented: &[u8] = if rev { rev_seq } else { seq };
    let sep = placement_sep(seq.len()) as i32;
    let text = &ctx.ref_seqs[partner.rid][win_start as usize..win_end as usize];
    let placement = crate::placements_x4(oriented, &[text], max_ed, sep)[0];

    if placement.best > max_ed {
        return None;
    }

    // The hit's window is re-centred on the placement found, as short_align expects.
    let end = win_start + placement.end as i64;
    let pad = max_ed.max(0) as i64 + WINDOW_PAD as i64;
    let (hit_start, hit_end) = ((end - read_len - pad).max(0), (end + pad).min(rlen));
    let (mut comp, mut comp_loc) = (Vec::new(), Vec::new());

    // The window's runner-up becomes the rescued hit's verified competitor.
    if placement.alt <= max_ed {
        comp.push(placement.alt);
        comp_loc.push((
            partner.rid as u32,
            rev,
            win_start + placement.alt_end as i64,
        ));
    }

    // A rescued mate has no seeds of its own, so its seed evidence is zero and seed_cap leaves it uncapped.
    Some(ShortHit {
        edit_dist: placement.best,
        comp,
        comp_loc,
        rep_frac: 0.0,
        seed_ev: [0; SEED_EV_FIELDS],
        gap: 0,
        rid: partner.rid,
        rev,
        win_start: hit_start as usize,
        win_end: hit_end as usize,
        end,
    })
}

/// A candidate placement of a mate, which is its best placement, a competitor its own seeds found, or one rescued next to the partner.
#[derive(Clone, Copy, Debug)]
pub(super) struct Cand {
    pub(super) edit_dist: i32,
    pub(super) rid: usize,
    pub(super) rev: bool,
    pub(super) end: i64,
    pub(super) rescued: bool,
}

/// Competitors further than this many edits from a mate's best placement are not pairing candidates (they weigh < 1e-7).
pub(super) const MATE_CAND_EDITS: i32 = 3;

/// Most pairing candidates per mate.
pub(super) const MAX_MATE_CANDS: usize = 64;

/// Pair combinations scoring within this relative margin of the best are ties.
pub(super) const PAIR_TIE_TOLERANCE: f64 = 1e-9;

/// The tie-breaking pair hash rotates mate 2's read hash by PAIR_HASH_ROTATE bits and folds itself by PAIR_HASH_FOLD bits.
pub(super) const PAIR_HASH_ROTATE: u32 = 17;

pub(super) const PAIR_HASH_FOLD: u32 = 29;

/// Smallest probability that a mate's chosen placement is wrong.
pub(super) const MATE_WRONG_PROB_FLOOR: f64 = 1e-12;

/// Mates are rescued unless the pair is proper and both are wrong with at most this probability.
pub(super) const RESCUE_WRONG_PROB: f64 = 1e-3;

/// Mate rescue aligns next to at most this many equally good partner candidates.
pub(super) const MAX_RESCUE_ANCHORS: usize = 4;

/// OPTIK_MAPQ_DEBUG writes this for a mate without a window MAPQ.
pub(super) const UNSET_DEBUG_MAPQ: i32 = -1;

/// Candidates of a short mate, the best placement and every competitor within MATE_CAND_EDITS of it, at most MAX_MATE_CANDS.
pub(super) fn mate_cands(hit: &ShortHit) -> Vec<Cand> {
    let mut cands = vec![Cand {
        edit_dist: hit.edit_dist,
        rid: hit.rid,
        rev: hit.rev,
        end: hit.end,
        rescued: false,
    }];

    for (&comp_ed, &(rid, rev, end)) in hit.comp.iter().zip(&hit.comp_loc) {
        if comp_ed <= hit.edit_dist + MATE_CAND_EDITS && cands.len() < MAX_MATE_CANDS {
            cands.push(Cand {
                edit_dist: comp_ed,
                rid: rid as usize,
                rev,
                end,
                rescued: false,
            });
        }
    }

    cands
}

/// Per-base floor on the edit rate for the sample's own differences from the reference.
/// A real genome's copy of a duplicated region can hold variants its reference copy lacks and a paralog shares, so a read may fit the paralog better by a few edits.
/// Simulated reads, drawn from the reference, never show it.
/// The floor was fitted on half of the real HG002 reads that their pair anchors, as the smallest value at which MAPQ >= 30 is calibrated with 23 wrong against 22.9 implied, and checked on the other half at 23 against 23.1.
/// Without it, reads with 0 edits at a paralog and 3 at their true copy were wrong at MAPQ 40-49 four times as often as implied.
pub(super) const SAMPLE_DIVERGENCE: f64 = 0.016;

/// edit_ratio's constants were fitted the same way as posterior_mapq_w's.
/// In order they are the prior edits added to the read's own, the substitution odds factor, the floor of 1 - e and the bounds of the ratio.
pub(super) const ER_PRIOR_EDITS: f64 = 1.2608;

pub(super) const ER_SUBST_ODDS: f64 = 0.8228;

pub(super) const ER_MIN_MATCH_RATE: f64 = 1e-3;

pub(super) const ER_MIN_RATIO: f64 = 1e-6;

pub(super) const ER_MAX_RATIO: f64 = 0.9;

/// The likelihood ratio of one more edit for a read with `edits` edits, i.e. the substitution odds 0.82 e / (1 - e) of the read's estimated error rate e = (ed + 1.26) / len + SAMPLE_DIVERGENCE, with constants fitted as in posterior_mapq_w.
pub(super) fn edit_ratio(edits: i32, len: usize) -> f64 {
    let err_rate = (edits as f64 + ER_PRIOR_EDITS) / len.max(1) as f64 + SAMPLE_DIVERGENCE;

    (ER_SUBST_ODDS * err_rate / (1.0 - err_rate).max(ER_MIN_MATCH_RATE))
        .clamp(ER_MIN_RATIO, ER_MAX_RATIO)
}

/// Whether two candidates are the same placement, with ends within `sep`.
pub(super) fn same_place(cand_a: &Cand, cand_b: &Cand, sep: i64) -> bool {
    cand_a.rid == cand_b.rid && cand_a.rev == cand_b.rev && (cand_a.end - cand_b.end).abs() <= sep
}

/// The ShortHit of a candidate.
/// The mate's own best is its hit, and any other candidate is re-verified in a window around its end (the aligner reports the same placement), so the record is built as for any short read.
pub(super) fn cand_hit(
    cand: &Cand,
    own: Option<&ShortHit>,
    seq: &[u8],
    rev_seq: &[u8],
    max_ed: i32,
    ctx: &AlignCtx<'_>,
) -> Option<ShortHit> {
    if let Some(own_hit) = own
        && own_hit.rid == cand.rid
        && own_hit.rev == cand.rev
        && own_hit.end == cand.end
    {
        return Some(clone_hit(own_hit));
    }

    let read_len = seq.len() as i64;
    let rlen = ctx.ref_seqs[cand.rid].len() as i64;
    // The window around the candidate's end leaves room for max_ed indels plus WINDOW_PAD bases on each side.
    let pad = max_ed.max(0) as i64 + WINDOW_PAD as i64;
    let (win_start, win_end) = (
        (cand.end - read_len - pad).max(0),
        (cand.end + pad).min(rlen),
    );

    if win_end - win_start < read_len {
        return None;
    }

    let oriented: &[u8] = if cand.rev { rev_seq } else { seq };
    let sep = placement_sep(seq.len()) as i32;
    let placement = crate::placements_x4(
        oriented,
        &[&ctx.ref_seqs[cand.rid][win_start as usize..win_end as usize]],
        max_ed,
        sep,
    )[0];

    if placement.best > max_ed {
        return None;
    }

    // The seed evidence stays with the mate's own hit, so a re-verified candidate keeps the same hidden-locus odds.
    Some(ShortHit {
        edit_dist: placement.best,
        comp: Vec::new(),
        comp_loc: Vec::new(),
        rep_frac: own.map_or(0.0, |own_hit| own_hit.rep_frac),
        seed_ev: own.map_or([0; SEED_EV_FIELDS], |own_hit| own_hit.seed_ev),
        gap: 0,
        rid: cand.rid,
        rev: cand.rev,
        win_start: win_start as usize,
        win_end: win_end as usize,
        end: win_start + placement.end as i64,
    })
}

/// A copy of a ShortHit, which holds vectors and so is not Copy.
pub(super) fn clone_hit(hit: &ShortHit) -> ShortHit {
    ShortHit {
        edit_dist: hit.edit_dist,
        comp: hit.comp.clone(),
        comp_loc: hit.comp_loc.clone(),
        rep_frac: hit.rep_frac,
        seed_ev: hit.seed_ev,
        gap: hit.gap,
        rid: hit.rid,
        rev: hit.rev,
        win_start: hit.win_start,
        win_end: hit.win_end,
        end: hit.end,
    }
}

/// Each mate's chosen hit with its pair MAPQ and whether rescue placed it, or None for an unplaced mate.
type ChosenMates = [Option<(ShortHit, i32, bool)>; MATES];

/// Pairs two independently placed mates (see the section comment), with a joint choice over both mates' candidates, mate rescue, pair MAPQ and records.
#[allow(clippy::too_many_arguments)]
/// Picks the most likely pair of placements from both mates' candidates, rescuing a mate inside the insert window when that helps.
/// Each chosen hit comes with its MAPQ before the cap and whether it was rescued, and the flag says whether the pick is a proper pair.
fn choose_pair(
    hits: &[Option<ShortHit>; MATES],
    seqs: &[(&[u8], &[u8]); MATES],
    len: &[usize; MATES],
    maxe: &[i32; MATES],
    model: &InsertModel,
    ctx: &AlignCtx<'_>,
) -> (ChosenMates, bool) {
    let (seq1, seq2) = (seqs[0].0, seqs[1].0);
    let as_place = |cand: &Cand, mate: usize| Place {
        rid: cand.rid,
        rev: cand.rev,
        start: cand.end - len[mate] as i64,
        end: cand.end,
    };
    let mut chosen: ChosenMates = [None, None];

    let mut cands: [Vec<Cand>; MATES] = [
        hits[0].as_ref().map_or(Vec::new(), mate_cands),
        hits[1].as_ref().map_or(Vec::new(), mate_cands),
    ];
    // Joint likelihood of every candidate combination, where an unmapped mate contributes a single "absent" option, giving the chosen combination and each mate's chance of being placed wrong.
    let resolve =
        |cands: &[Vec<Cand>; MATES]| -> ((Option<usize>, Option<usize>), [f64; MATES], bool) {
            let best_ed = [
                cands[0].iter().map(|cand| cand.edit_dist).min(),
                cands[1].iter().map(|cand| cand.edit_dist).min(),
            ];
            let ratio = [
                best_ed[0].map(|best| edit_ratio(best, len[0])),
                best_ed[1].map(|best| edit_ratio(best, len[1])),
            ];
            let opts = |mate: usize| -> Vec<Option<usize>> {
                if cands[mate].is_empty() {
                    vec![None]
                } else {
                    (0..cands[mate].len()).map(Some).collect()
                }
            };
            // Each candidate's likelihood relative to its mate's best is r^(extra edits), with r the per-edit likelihood ratio.
            let mut combos: Vec<(Option<usize>, Option<usize>, f64, bool)> = Vec::new();

            for cand1 in opts(0) {
                for cand2 in opts(1) {
                    let mut likelihood = 1.0;

                    if let Some(idx1) = cand1 {
                        likelihood *= ratio[0]
                            .expect("a chosen mate-1 candidate implies its ratio")
                            .powi(cands[0][idx1].edit_dist - best_ed[0].expect("a chosen mate-1 candidate implies a best distance"));
                    }

                    if let Some(idx2) = cand2 {
                        likelihood *= ratio[1]
                            .expect("a chosen mate-2 candidate implies its ratio")
                            .powi(cands[1][idx2].edit_dist - best_ed[1].expect("a chosen mate-2 candidate implies a best distance"));
                    }

                    let mut is_proper = false;

                    // A pair outside the insert model is still allowed, but its likelihood is scaled down by the discordant prior.
                    if let (Some(idx1), Some(idx2)) = (cand1, cand2) {
                        is_proper = proper(
                            as_place(&cands[0][idx1], 0),
                            as_place(&cands[1][idx2], 1),
                            model,
                        );
                        if !is_proper {
                            let worst_ratio =
                                ratio[0].unwrap_or(0.0).max(ratio[1].unwrap_or(0.0));

                            likelihood *= discordant_prior(worst_ratio);
                        }
                    }

                    combos.push((cand1, cand2, likelihood, is_proper));
                }
            }

            // Combinations within PAIR_TIE_TOLERANCE of the best are ties.
            // A hash of both sequences picks among them, so a rerun or another thread count gives the same choice.
            let top = combos.iter().map(|combo| combo.2).fold(0.0f64, f64::max);
            let tied: Vec<usize> = (0..combos.len())
                .filter(|&combo_idx| combos[combo_idx].2 >= top * (1.0 - PAIR_TIE_TOLERANCE))
                .collect();
            let mut pair_hash = read_hash(seq1) ^ read_hash(seq2).rotate_left(PAIR_HASH_ROTATE);

            pair_hash ^= pair_hash >> PAIR_HASH_FOLD;
            let pick = combos[tied[(pair_hash % tied.len() as u64) as usize]];
            // A mate's chance of being wrong is the likelihood share of combinations that put it somewhere else.
            let total: f64 = combos.iter().map(|combo| combo.2).sum();
            let mut wrong_prob = [1.0f64; MATES];

            for (mate, cand_idx) in [(0usize, pick.0), (1usize, pick.1)] {
                let Some(cand_idx) = cand_idx else { continue };
                let cand = cands[mate][cand_idx];
                let sep = placement_sep(len[mate]) as i64;
                // The combinations that place this mate elsewhere.
                let other: f64 = combos
                    .iter()
                    .filter(|combo| {
                        let mate_pick = if mate == 0 { combo.0 } else { combo.1 };

                        mate_pick.is_none_or(|pick_idx| {
                            !same_place(&cands[mate][pick_idx], &cand, sep)
                        })
                    })
                    .map(|combo| combo.2)
                    .sum();

                wrong_prob[mate] = (other / total).max(MATE_WRONG_PROB_FLOOR);
            }

            ((pick.0, pick.1), wrong_prob, pick.3)
        };
    let (mut pick, mut wrong_prob, mut is_proper) = resolve(&cands);

    // Mate rescue.
    // Unless the chosen pair is proper and both mates are confident, each mate is aligned inside the window its insert size allows next to the partner's best candidates, all equally good ones up to 4, since its true locus may have lost every seed.
    if !is_proper || wrong_prob.iter().any(|&prob| prob > RESCUE_WRONG_PROB) {
        let mut added = false;

        for x in 0..MATES {
            let anchor = 1 - x;

            if cands[anchor].is_empty() {
                continue;
            }

            let best = cands[anchor]
                .iter()
                .map(|cand| cand.edit_dist)
                .min()
                .expect("the anchor mate has at least one candidate");
            let anchors: Vec<Cand> = cands[anchor]
                .iter()
                .filter(|cand| cand.edit_dist == best)
                .take(MAX_RESCUE_ANCHORS)
                .copied()
                .collect();
            let (mate_seq, mate_rc) = seqs[x];
            let sep = placement_sep(mate_seq.len()) as i64;

            // The rescued placement joins the mate's candidates, and so does its runner-up if it is within MATE_CAND_EDITS.
            // Then the joint choice is made again over the larger set.
            for anchor_cand in anchors {
                if let Some(rescued) = rescue_in_window(
                    mate_seq,
                    mate_rc,
                    as_place(&anchor_cand, anchor),
                    model,
                    maxe[x],
                    ctx,
                ) {
                    let mut found = vec![Cand {
                        edit_dist: rescued.edit_dist,
                        rid: rescued.rid,
                        rev: rescued.rev,
                        end: rescued.end,
                        rescued: true,
                    }];

                    if let (Some(&comp_ed), Some(&(rid, rev, end))) =
                        (rescued.comp.first(), rescued.comp_loc.first())
                        && comp_ed <= rescued.edit_dist + MATE_CAND_EDITS
                    {
                        found.push(Cand {
                            edit_dist: comp_ed,
                            rid: rid as usize,
                            rev,
                            end,
                            rescued: true,
                        });
                    }

                    // A placement the mate already lists is not added twice.
                    for cand in found {
                        if !cands[x]
                            .iter()
                            .any(|existing| same_place(existing, &cand, sep))
                        {
                            cands[x].push(cand);
                            added = true;
                        }
                    }
                }
            }
        }

        if added {
            (pick, wrong_prob, is_proper) = resolve(&cands);
        }
    }

    // The pair MAPQ is the Phred-scaled wrong probability from the joint choice, and cand_hit realigns each chosen candidate to a full hit.
    let proper_pick = is_proper;

    for (mate, cand_idx) in [(0usize, pick.0), (1usize, pick.1)] {
        let Some(cand_idx) = cand_idx else { continue };
        let cand = cands[mate][cand_idx];
        let pair_mapq =
            ((-PHRED_SCALE * wrong_prob[mate].log10()).round() as i32).clamp(0, SAM_MAX_MAPQ);
        let (mate_seq, mate_rc) = seqs[mate];

        if let Some(hit) = cand_hit(
            &cand,
            hits[mate].as_ref(),
            mate_seq,
            mate_rc,
            maxe[mate],
            ctx,
        ) {
            chosen[mate] = Some((hit, pair_mapq, cand.rescued));
        }
    }

    (chosen, proper_pick)
}

/// Builds both mates' SAM records.
/// Placements are chosen jointly when both mates took the short-read path, and the pair fields are filled from the primaries written.
#[allow(clippy::too_many_arguments)]
pub(super) fn pair_records(
    name: &str,
    seq1: &[u8],
    seq2: &[u8],
    aln1: MateAln,
    aln2: MateAln,
    model: &InsertModel,
    max_ed: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> (Vec<String>, Vec<String>) {
    let (rc1, rc2) = (revcomp(seq1), revcomp(seq2));
    let seqs: [(&[u8], &[u8]); MATES] = [(seq1, &rc1), (seq2, &rc2)];
    let len = [seq1.len(), seq2.len()];
    let maxe = [
        read_max_ed(max_ed, seq1.len()),
        read_max_ed(max_ed, seq2.len()),
    ];
    let mut lines: [Option<Vec<String>>; MATES] = [None, None];
    let mut hits: [Option<ShortHit>; MATES] = [None, None];

    // A mate arrives either as finished long-read lines or as a short-read hit still to pair.
    for (i, aln) in [aln1, aln2].into_iter().enumerate() {
        match aln {
            MateAln::Lines(mate_lines) => lines[i] = Some(mate_lines),
            MateAln::Short(hit) => hits[i] = hit,
        }
    }

    // Mates already mapped by the long-read path come back as finished lines, and only short hits go through the joint choice.
    let (mut chosen, proper_pick) = if lines[0].is_none() && lines[1].is_none() {
        choose_pair(&hits, &seqs, &len, &maxe, model, ctx)
    } else {
        ([None, None], false)
    };

    // The hidden-locus odds and the seed-evidence cap of the unpaired MAPQ apply to the paired posterior too.
    // In a proper pair whose window next to the partner was aligned, a mate is wrong only if its true copy went unseen and either the partner's did too, or the true pair is discordant.
    // The first means the whole fragment came from another copy, with odds H_x H_p since the mates' seeds are drawn independently, and the second has odds H_x discordant_prior.
    // Where the mate's own hidden-locus odds bind, the window is aligned, and if the mate's placement is the best there, its hidden odds become H_x (H_p + discordant_prior) and the window's runner-up is a verified competitor.
    // Real HG002 pairs inside diverged L1 or segmental-duplication copies kept MAPQ 1-15 where BWA-MEM gives them 40 or more and GIAB confirms the placements, which affected 0.7% of chr20 SNPs.
    let mut alns: [Option<(ShortAln, i32)>; MATES] = [None, None];
    let mut hids = [f64::INFINITY; MATES];
    let mut window_mq: [Option<i32>; MATES] = [None, None];
    let mut debug_tags: [String; MATES] = [String::new(), String::new()]; // OPTIK_MAPQ_DEBUG pair inputs

    // First, each mate's MAPQ is the smallest of the pair posterior, its unpaired hidden-locus bound and its seed cap.
    for x in 0..MATES {
        if let Some((hit, pair_mapq, _)) = chosen[x].as_ref() {
            let aln = short_align(seqs[x].0, seqs[x].1, hit, ctx);

            hids[x] = hidden_odds(hit, aln.events, len[x], ctx.rho);
            let mapq_val = (*pair_mapq)
                .min(posterior_mapq_w(hit.edit_dist, &[], len[x], hids[x]))
                .min(seed_cap(hit));

            if mapq_debug() {
                debug_tags[x] = format!(
                    "pair={},hidden={},seed={}",
                    pair_mapq,
                    posterior_mapq_w(hit.edit_dist, &[], len[x], hids[x]),
                    seed_cap(hit)
                );
            }

            alns[x] = Some((aln, mapq_val));
        }
    }

    // Next, a mate held below its pair MAPQ in a proper pair gets its window next to the partner aligned.
    // When that alignment reproduces the same placement, the combined hidden odds described above replace its own.
    for x in 0..MATES {
        let (Some((hit, pair_mapq, _)), Some((partner_hit, _, _))) =
            (chosen[x].as_ref(), chosen[1 - x].as_ref())
        else {
            continue;
        };
        let mapq_val = alns[x].as_ref().map_or(0, |aln| aln.1);

        // Only a mate held below its pair bound can gain from rescue.
        if !proper_pick || mapq_val >= (*pair_mapq).min(seed_cap(hit)) {
            continue;
        }

        let partner = Place {
            rid: partner_hit.rid,
            rev: partner_hit.rev,
            start: partner_hit.end - len[1 - x] as i64,
            end: partner_hit.end,
        };
        let sep = placement_sep(len[x]) as i64;

        if let Some(rescued) = rescue_in_window(seqs[x].0, seqs[x].1, partner, model, maxe[x], ctx)
            && rescued.rid == hit.rid
            && rescued.rev == hit.rev
            && (rescued.end - hit.end).abs() <= sep
            // The rescue must find the very same placement. Its window competitors then stand in for the genome-wide ones.
            && rescued.edit_dist == hit.edit_dist
        {
            let comp: Vec<(i32, f64)> =
                rescued.comp.iter().map(|&comp_ed| (comp_ed, 1.0)).collect();
            let hid = hids[x] * (hids[1 - x] + discordant_prior(edit_ratio(hit.edit_dist, len[x])));

            window_mq[x] = Some(
                (*pair_mapq)
                    .min(posterior_mapq_w(hit.edit_dist, &comp, len[x], hid))
                    .min(seed_cap(hit)),
            );
        }
    }

    // The window result can only raise a mate's MAPQ.
    for x in 0..MATES {
        if let (Some(window_mapq), Some(aln)) = (window_mq[x], alns[x].as_mut()) {
            aln.1 = aln.1.max(window_mapq);
        }
    }

    // A mate placed by rescue has no evidence of its own.
    // It sits where its partner's placement puts it, so it is no more certain than the partner (a pair inside a segmental duplication whose partner is ambiguous was placed in the wrong copy at MAPQ 36-42).
    let mqs = [
        alns[0].as_ref().map(|aln| aln.1),
        alns[1].as_ref().map(|aln| aln.1),
    ];

    for x in 0..MATES {
        if let (Some((_, _, true)), Some(aln), Some(partner_mapq)) =
            (chosen[x].as_ref(), alns[x].as_mut(), mqs[1 - x])
        {
            aln.1 = aln.1.min(partner_mapq);
        }
    }

    // Pairing combines the mates' evidence but cannot create any, so a mate's MAPQ is at most the sum of both mates' unpaired MAPQs, as independent evidence adds in log odds.
    // Two mates that each tie alone lie in a repeat where every copy pair is proper, and the joint candidate lists are then incomplete.
    // Pairs inside chr7's centromeric satellite went to the wrong copy at MAPQ 10-14, 2.3x the implied error.
    let unpaired = |mate: usize| {
        hits[mate]
            .as_ref()
            .map_or(0, |hit| final_mapq(hit, hit.edit_dist, len[mate], ctx.rho))
    };
    let evidence = unpaired(0) + unpaired(1);

    for aln in alns.iter_mut().flatten() {
        aln.1 = aln.1.min(evidence);
    }

    if mapq_debug() {
        for x in 0..MATES {
            if !debug_tags[x].is_empty() {
                debug_tags[x] = format!(
                    "\tXP:Z:{},window={},proper={},sum={}",
                    debug_tags[x],
                    window_mq[x].unwrap_or(UNSET_DEBUG_MAPQ),
                    proper_pick as u8,
                    evidence
                );
            }
        }
    }

    // Each mate's record comes from the best source available, in this order.
    // Lines from the long-read path are kept as they are, then a paired short hit, then its unpaired hit, and finally an unmapped record.
    let mut out: [Vec<String>; MATES] = [Vec::new(), Vec::new()];

    for x in 0..MATES {
        let (mate_seq, mate_rc) = seqs[x];

        out[x] = if let Some(mate_lines) = lines[x].take() {
            mate_lines
        } else if let (Some((hit, _, _)), Some((aln, mapq_val))) =
            (chosen[x].take(), alns[x].take())
        {
            vec![
                short_format(
                    name, mate_seq, mate_rc, &hit, &aln, mapq_val, ctx, ref_names,
                ) + &debug_tags[x],
            ]
        } else if let Some(hit) = hits[x].as_ref() {
            vec![short_record(name, mate_seq, hit, ctx, ref_names, None)]
        } else {
            let seq_s = if mate_seq.is_empty() {
                "*"
            } else {
                std::str::from_utf8(mate_seq).unwrap_or("*")
            };

            vec![format!(
                "{name}\t{SAM_FLAG_UNMAPPED}\t*\t0\t0\t*\t*\t0\t0\t{seq_s}\t*"
            )]
        };
    }

    // The proper-pair flag is decided from the primaries as written, so it agrees with the records even where a mate's placement came from an unpaired path.
    let [mut lines1, mut lines2] = out;
    let (prim1, prim2) = (parse_prim(&lines1[0]), parse_prim(&lines2[0]));
    let proper_pair = prim1.mapped
        && prim2.mapped
        && {
            let place_of = |prim: &Prim| {
                ref_names
                    .iter()
                    .position(|ref_name| *ref_name == prim.rname)
                    .map(|rid| Place {
                        rid,
                        rev: prim.rev,
                        start: prim.pos - 1,
                        end: prim.end,
                    })
            };

            matches!((place_of(&prim1), place_of(&prim2)), (Some(place1), Some(place2)) if proper(place1, place2, model))
        };

    set_pair_fields(&mut lines1, &mut lines2, proper_pair);
    (lines1, lines2)
}

/// Fills the SAM pair fields of both mates' records, primary first.
pub(super) fn set_pair_fields(lines1: &mut [String], lines2: &mut [String], proper_pair: bool) {
    let (prim1, prim2) = (parse_prim(&lines1[0]), parse_prim(&lines2[0]));
    // Template length from the primaries' outer coordinates.
    let tlen = if prim1.mapped && prim2.mapped && prim1.rname == prim2.rname {
        let (start, end) = (prim1.pos.min(prim2.pos), prim1.end.max(prim2.end));

        end - start + 1
    } else {
        0
    };
    // The leftmost mate gets the positive TLEN. A tie goes to the forward mate, or to mate 1 when both share a strand.
    let first_is_left =
        prim1.pos < prim2.pos || (prim1.pos == prim2.pos && (prim1.rev == prim2.rev || !prim1.rev));

    for (mate_i, (lines, this_mate, other)) in [
        (&mut *lines1, &prim1, &prim2),
        (&mut *lines2, &prim2, &prim1),
    ]
    .into_iter()
    .enumerate()
    {
        // Supplementary records get mate fields too, but only the primary gets a TLEN.
        for (j, line) in lines.iter_mut().enumerate() {
            let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();

            if fields.len() < SAM_MANDATORY_COLUMNS {
                continue;
            }

            let mut flag: u16 = fields[SAM_COL_FLAG].parse().unwrap_or(0);

            flag |= SAM_FLAG_PAIRED
                | if mate_i == 0 {
                    SAM_FLAG_READ1
                } else {
                    SAM_FLAG_READ2
                };
            if proper_pair && this_mate.mapped && other.mapped {
                flag |= SAM_FLAG_PROPER_PAIR;
            }

            if !other.mapped {
                flag |= SAM_FLAG_MATE_UNMAPPED;
            } else if other.rev {
                flag |= SAM_FLAG_MATE_REVERSE;
            }

            let this_mapped = flag & SAM_FLAG_UNMAPPED == 0;

            if !this_mapped && other.mapped {
                // An unmapped mate takes its partner's position, per SAM 1.4.
                fields[SAM_COL_RNAME] = other.rname.clone();
                fields[SAM_COL_POS] = other.pos.to_string();
            }

            fields[SAM_COL_FLAG] = flag.to_string();

            if other.mapped {
                fields[SAM_COL_RNEXT] = if fields[SAM_COL_RNAME] == other.rname {
                    "=".into()
                } else {
                    other.rname.clone()
                };
                fields[SAM_COL_PNEXT] = other.pos.to_string();
            // Its mate is unplaced and borrowed this record's position, so RNEXT points back here.
            } else if this_mapped {
                fields[SAM_COL_RNEXT] = "=".into();
                fields[SAM_COL_PNEXT] = fields[SAM_COL_POS].clone();
            }

            let left = if mate_i == 0 {
                first_is_left
            } else {
                !first_is_left
            };

            fields[SAM_COL_TLEN] = if j == 0 && tlen != 0 {
                (if left { tlen } else { -tlen }).to_string()
            } else {
                "0".into()
            };
            if other.mapped {
                fields.push(format!("MC:Z:{}", other.cigar));
                fields.push(format!("MQ:i:{}", other.mapq));
            }

            *line = fields.join("\t");
        }
    }
}
