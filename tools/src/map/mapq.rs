//! Mapping quality combines the posterior over verified competitors with the odds that the true locus was never seeded.

use super::*;

/// Setting OPTIK_MAPQ_DEBUG adds an XD:Z tag of the MAPQ inputs to each short-read record.
/// The tag lists the edit distance, events, loci, minimizers, repetitive minimizers, votes for the best and the other loci, and every competitor with its edit distance.
pub(super) fn mapq_debug() -> bool {
    static D: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *D.get_or_init(|| std::env::var_os("OPTIK_MAPQ_DEBUG").is_some())
}

/// MAPQ of a short-read hit, the posterior over its verified competitors and the odds that its true locus was never a candidate (hidden_odds), capped by its seed evidence (seed_cap).
pub(super) fn final_mapq(hit: &ShortHit, events: i32, len: usize, rho: f64) -> i32 {
    let weighted: Vec<(i32, f64)> = hit.comp.iter().map(|&comp_ed| (comp_ed, 1.0)).collect();

    posterior_mapq_w(
        hit.edit_dist,
        &weighted,
        len,
        hidden_odds(hit, events, len, rho),
    )
    .min(seed_cap(hit))
}

/// MAPQ = -10 log10 P(best locus wrong), with P = (S + H) / (1 + S + H).
/// S sums r^(edits - best edits) over the verified distinct competitors, where r = edit_ratio is the likelihood ratio of one more edit.
/// H is hidden_odds, the odds that the true locus was never verified.
///
/// A competitor as good as the best is a tie with MAPQ 0.
///
/// The constants were fitted by maximum likelihood on simulated 100, 150 and 250 bp reads at 0.5-8% error.
/// Those reads came from 14 genomes, which were E. coli, Salmonella, four synthetic repeat genomes, hg38 chr1, chr2, chr7, chr17, chr19, chr21, chrX and the whole assembly.
/// They are judged on fresh simulated reads with bench/mapq_calibration.py --multi, on simulated pairs with paired_bench.py, and on real reads with bench/real_reads.py, using pair-anchored HG002 and E. coli mates and GIAB HG002 chr20 variant calling.
/// On none of them does a MAPQ bin above 0 have significantly more wrong reads than its MAPQs imply.
/// The real reads also exposed three cases that simulation missed.
/// Diverged Alu and L1 copies are handled by the competitor gap in hidden_odds, and satellites with few seeded minimizers by its seeded-share weight.
/// Then SAMPLE_DIVERGENCE covers variants of the sample that a paralog shares.
///
/// posterior_mapq_w takes a prior weight per competitor, which pairing uses to weigh a placement its mate does not support by discordant_prior.
/// A competitor at full weight as good as the best is a tie at 0, while a discounted one may even beat the best, which is then the mate-supported placement.
pub(super) fn posterior_mapq_w(best_ed: i32, competitors: &[(i32, f64)], len: usize, hidden: f64) -> i32 {
    if competitors
        .iter()
        .any(|&(comp_ed, weight)| comp_ed <= best_ed && weight >= 1.0)
    {
        return 0;
    }

    let ratio = edit_ratio(best_ed, len);
    let odds_sum: f64 = competitors
        .iter()
        .map(|&(comp_ed, weight)| weight * ratio.powi(comp_ed - best_ed))
        .sum::<f64>()
        + hidden;
    let p_wrong = (odds_sum / (1.0 + odds_sum)).max(P_WRONG_FLOOR);

    ((-PHRED_SCALE * p_wrong.log10()).round() as i32).clamp(0, SAM_MAX_MAPQ)
}

/// Smallest probability that the reported placement is wrong, which caps the posterior MAPQ before SAM_MAX_MAPQ.
pub(super) const P_WRONG_FLOOR: f64 = 1e-9;

// The HO_ constants are the coefficients of hidden_odds, fitted as its doc comment describes.

/// Seed votes of the best locus count up to this many.
pub(super) const HO_BEST_VOTES_CAP: u32 = 25;

/// The competitor term is weighted by 1 - HO_RHO_WEIGHT rho, times (1 - over-occurring share)^HO_REP_SHARE_EXP, times 1 - e^(-HO_VOTE_SATURATION votes).
pub(super) const HO_RHO_WEIGHT: f64 = 1.5407;

pub(super) const HO_REP_SHARE_EXP: f64 = 1.0366;

pub(super) const HO_VOTE_SATURATION: f64 = 8.1842;

/// Edits to the nearest verified competitor count up to this many.
pub(super) const HO_GAP_CAP: f64 = 4.2475;

/// Base log-odds and its per-vote slope.
pub(super) const HO_INTERCEPT: f64 = -2.8093;

pub(super) const HO_PER_VOTE: f64 = 0.5383;

/// Added when the candidate list was cut at MAX_SHORT_CANDIDATES.
pub(super) const HO_CUT_LIST: f64 = 0.485;

/// Slope in the alignment's event rate (mismatches, gap openings and clipped ends per base).
pub(super) const HO_EVENT_RATE: f64 = 71.0931;

/// Slope in the read's minimizer count, which counts up to HO_MINIMIZER_CAP.
pub(super) const HO_PER_MINIMIZER: f64 = 0.0527;

pub(super) const HO_MINIMIZER_CAP: u32 = 60;

/// Slope in the genome's repeat fraction rho.
pub(super) const HO_RHO: f64 = 9.0619;

/// Reads with over-occurring minimizers add their own intercept, a slope in the over-occurring share, and separate vote slopes below and above HO_REP_VOTE_SPLIT votes.
pub(super) const HO_REP_INTERCEPT: f64 = 2.8365;

pub(super) const HO_REP_SHARE: f64 = 5.2724;

pub(super) const HO_REP_VOTE_LOW: f64 = 0.1715;

pub(super) const HO_REP_VOTE_SPLIT: f64 = 10.0;

pub(super) const HO_REP_VOTE_HIGH: f64 = 0.0362;

/// The competitor term for reads with over-occurring minimizers, per edit to the nearest competitor and for when none lies within MAX_GAP.
pub(super) const HO_REP_PER_GAP: f64 = 2.5081;

pub(super) const HO_REP_NO_COMPETITOR: f64 = 9.7409;

/// Competitor term without over-occurring minimizers.
pub(super) const HO_PER_GAP: f64 = 1.6371;

pub(super) const HO_NO_COMPETITOR: f64 = 8.5511;

/// Log-odds are capped here before exponentiation.
pub(super) const HO_MAX_LOG_ODDS: f64 = 30.0;

/// Odds that the read's true locus was never a verified candidate, because its seeds were destroyed by read errors, filtered as over-occurring, or cut with the candidates beyond MAX_SHORT_CANDIDATES while another copy kept some.
/// Each seed vote of the best locus divides the odds by about 1.7, and a candidate list cut at MAX_SHORT_CANDIDATES multiplies them by about 1.6, since copies beyond the cut were never aligned.
/// Over-occurring minimizers multiply them by e^(2.8 + 5.3 x their fraction), and the best alignment's event rate f = events / len, counting mismatches, gap openings and clipped ends, by e^(71 f), since at high divergence the true copy is the one most likely to lose its seeds.
/// Each edit between the best placement and its nearest verified competitor divides them by about 5, or about 12 in a repeat family, up to about 4.2 edits.
/// Having no verified competitor within 12 edits divides them by e^8.6, or e^9.7 in a repeat family.
/// Both terms are weighted by the share of seeded minimizers and by the genome's repeat fraction rho, as described below.
/// Without the competitor term, reads of diverged Alu and L1 copies with 0-1 edits got MAPQ 7-14 on real HG002 reads, where BWA-MEM gives 60 and GIAB confirms the placement.
/// Like the posterior, the model was fitted by maximum likelihood on simulated reads of 14 genomes and checked on real reads with bench/real_reads.py.
pub(super) fn hidden_odds(hit: &ShortHit, events: i32, len: usize, rho: f64) -> f64 {
    let [loci, n_min, n_rep, n_best, _] = hit.seed_ev;
    let event_rate = events.max(0) as f64 / len.max(1) as f64;
    let best_votes = n_best.min(HO_BEST_VOTES_CAP) as f64;
    // The nearest verified competitor matters because a copy family that is diverged here, with every verified copy several edits away or none within 12, is unlikely to hold an unseen copy that fits the read as well.
    // That is evidence only as far as closer copies would have been seeded, so it is weighted by the share of the read's minimizers that seeded, (1 - over-occurring fraction)^1.04, and by 1 - 1.54 rho, since in a repeat-dense genome more copies are never seeded.
    // A rescued mate has no seed vote of its own, so the weight is zero for it.
    // Without the weights, reads of chr7's centromeric satellite with a few seeded minimizers and verified copies 10 or more edits away got MAPQ 10-17 and were wrong.
    let rep_share = n_rep as f64 / n_min.max(1) as f64;
    let seeded_weight = (1.0 - HO_RHO_WEIGHT * rho).max(0.0)
        * (1.0 - rep_share).powf(HO_REP_SHARE_EXP)
        * (1.0 - (-HO_VOTE_SATURATION * best_votes).exp());
    let none = hit.gap >= MAX_GAP;
    let gap = if none {
        0.0
    } else {
        (hit.gap as f64).min(HO_GAP_CAP)
    };
    let mut log_odds = HO_INTERCEPT - HO_PER_VOTE * best_votes
        + if loci > MAX_SHORT_CANDIDATES as u32 {
            HO_CUT_LIST
        } else {
            0.0
        }
        + HO_EVENT_RATE * event_rate
        - HO_PER_MINIMIZER * n_min.min(HO_MINIMIZER_CAP) as f64
        - HO_RHO * rho;

    if n_rep > 0 {
        // Over-occurring minimizers mean the read lies (partly) in a family of 1000 or more copies, of which only the seeded ones were verified.
        log_odds += HO_REP_INTERCEPT
            + HO_REP_SHARE * rep_share
            + HO_REP_VOTE_LOW * best_votes.min(HO_REP_VOTE_SPLIT)
            + HO_REP_VOTE_HIGH * (best_votes - HO_REP_VOTE_SPLIT).max(0.0)
            + seeded_weight
                * (-HO_REP_PER_GAP * gap + if none { -HO_REP_NO_COMPETITOR } else { 0.0 });
    } else {
        log_odds +=
            seeded_weight * (-HO_PER_GAP * gap + if none { -HO_NO_COMPETITOR } else { 0.0 });
    }

    log_odds.min(HO_MAX_LOG_ODDS).exp()
}

/// Cap from seed evidence alone, SEED_CAP_BASE + SEED_CAP_PER_VOTE per seed vote of the best locus.
/// A best locus with few seed votes may be a copy that only a read error made seedable, so one vote allows at most MAPQ 25.
/// Reads in repeat families with over-occurring minimizers get no separate cap, since `hidden_odds` models them, fitted on 14 simulated genomes and checked on real HG002 reads.
/// A hand-set cap held their MAPQ near 13 and cost variant recall.
/// A rescued mate without seeds of its own is not capped here, because its partner placed it.
pub(super) fn seed_cap(hit: &ShortHit) -> i32 {
    let [_, n_min, _, n_best, _] = hit.seed_ev;

    if n_min == 0 {
        return SAM_MAX_MAPQ;
    }

    const SEED_CAP_BASE: f64 = 20.0;
    const SEED_CAP_PER_VOTE: f64 = 5.0;
    let cap = SEED_CAP_BASE + SEED_CAP_PER_VOTE * n_best as f64;

    (cap.floor() as i32).clamp(0, SAM_MAX_MAPQ)
}
