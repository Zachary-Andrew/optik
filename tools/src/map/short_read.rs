//! Maps short reads.
//! Seed hits vote for diagonals, each candidate locus is verified with bit-parallel edit distance, and the best locus gets a gap-affine alignment.

use super::*;

/// Short-read mapping below CHAIN_MIN_LEN, with a single best record per read.
///
/// Seeds are hashed-order minimizers with strand.
/// Each reference hit votes for an alignment start on its diagonal, ref_pos - read_pos on the hit's strand.
/// Nearby starts merge into candidate loci, and each locus is verified once with the bit-parallel edit distance over a window of read_len + 2*pad around the diagonal, where pad = max_ed + 8, instead of re-aligning a window of 5x the read length for every seed hit.
/// MAPQ compares the best and second-best distinct loci and is 0 on a tie.
pub(super) fn map_one_read_short(
    read_name: &str,
    seq: &[u8],
    max_ed: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> String {
    let seq_s = if seq.is_empty() {
        "*"
    } else {
        std::str::from_utf8(seq).unwrap_or("*")
    };
    let unmapped = || format!("{read_name}\t{SAM_FLAG_UNMAPPED}\t*\t0\t0\t*\t*\t0\t0\t{seq_s}\t*");
    let Some(hit) = align_read_impl(seq, max_ed, ctx) else {
        return unmapped();
    };

    short_record(read_name, seq, &hit, ctx, ref_names, None)
}

/// A short read's alignment at a verified placement (`hit`), with an affine CIGAR in the placement's window, left-aligned indels and clipped ends.
pub(super) struct ShortAln {
    pub(super) cig: String,
    pub(super) offset: usize,
    pub(super) edit_dist: i32,
    pub(super) events: i32,
}

/// Aligns a short read at its verified placement, with a gap-affine CIGAR in the placement's window, left-aligned indels and low-scoring ends clipped.
pub(super) fn short_align(seq: &[u8], rev_seq: &[u8], hit: &ShortHit, ctx: &AlignCtx<'_>) -> ShortAln {
    let align_seq: &[u8] = if hit.rev { rev_seq } else { seq };
    let ref_window = &ctx.ref_seqs[hit.rid][hit.win_start..hit.win_end];
    // Gap-affine CIGAR from optik's own wavefront aligner with minimap2-like short-read penalties x=4, o=6, e=2.
    // The unit-cost scalar DP stays as the fallback and as the OPTIK_FORCE_SCALAR_CIGAR reference.
    let (cig, offset) = if force_scalar_cigar() {
        cigar(align_seq, ref_window)
    } else {
        crate::wfa_affine_cigar(align_seq, ref_window)
            .unwrap_or_else(|| cigar(align_seq, ref_window))
    };
    // Indels are left-aligned, the SAM/VCF convention that callers and pileups expect.
    let cig = left_align(&cig, align_seq, ref_window, offset);
    // Read ends that score below zero under extension scoring are soft-clipped, as the long-read flanks are and as minimap2 and BWA do.
    // An end-to-end CIGAR writes the last bases of a read that runs into an indel as mismatches, or starts the read inside deleted sequence, and variant callers count those bases as SNPs.
    let (cig, offset, edit_dist) = match clip_ends(&cig, align_seq, ref_window, offset) {
        Some((full, toff, edit_dist)) => (full, toff, edit_dist),
        None => {
            // NM is derived from the emitted CIGAR so the record is self-consistent.
            let edit_dist =
                cigar_ed_implied(&cig, align_seq, ref_window, offset).unwrap_or(hit.edit_dist);

            (cig, offset, edit_dist)
        }
    };
    let events = cigar_events(&cig, align_seq, ref_window, offset);

    ShortAln {
        cig,
        offset,
        edit_dist,
        events,
    }
}

/// Verified placement of one read.
pub(super) struct ShortHit {
    pub(super) edit_dist: i32,
    /// Distinct competitors' edit distances and the fraction of the read's minimizers dropped as over-occurring, both MAPQ inputs of final_mapq.
    pub(super) comp: Vec<i32>,
    /// Where each competitor lies, as (rid, rev, exclusive reference end).
    /// Pairing discounts the ones its mate does not support.
    pub(super) comp_loc: Vec<(u32, bool, i64)>,
    pub(super) rep_frac: f64,
    /// The seed evidence for hidden_odds is stored in this order.
    /// It counts the candidate loci before the MAX_SHORT_CANDIDATES cut, the read's minimizers and the over-occurring ones, then the seed votes of the best locus and of the strongest other locus.
    /// All are zero for a mate placed by rescue, which has no seeds of its own.
    pub(super) seed_ev: [u32; SEED_EV_FIELDS],
    /// Edits from the best placement to the nearest verified competitor, or MAX_GAP when none is that close, used by hidden_odds.
    /// It is 0 for a placement that is not the read's own best, such as a pair-chosen candidate or a rescued mate, whose own best elsewhere competes with it.
    pub(super) gap: i32,
    pub(super) rid: usize,
    pub(super) rev: bool,
    pub(super) win_start: usize,
    pub(super) win_end: usize,
    /// Exclusive reference end of the best placement.
    pub(super) end: i64,
}

/// Cap of ShortHit::gap, beyond which competitors over 12 edits away are not counted.
pub(super) const MAX_GAP: i32 = 13;

/// Cap on candidate loci verified per read, taken in order of seed support.
pub(super) const MAX_SHORT_CANDIDATES: usize = 128;

/// Bases of padding around a candidate window beyond the edit budget.
pub(super) const WINDOW_PAD: usize = 8;

/// ShortHit::seed_ev holds this many counts, in the order its field doc gives.
pub(super) const SEED_EV_FIELDS: usize = 5;

/// Placements whose ends differ by at most 1/PLACEMENT_SEP_FRACTION of the read, and at least PLACEMENT_SEP_MIN bases, are the same placement.
pub(super) const PLACEMENT_SEP_FRACTION: usize = 10;

pub(super) const PLACEMENT_SEP_MIN: usize = 10;

/// Placements whose reference ends differ by at most this, a tenth of the read and at least 10 bases, are the same placement.
/// The truth scoring uses the same read-origin tolerance.
pub(super) fn placement_sep(len: usize) -> usize {
    (len / PLACEMENT_SEP_FRACTION).max(PLACEMENT_SEP_MIN)
}

/// Finds and verifies a short read's best placement and its competitors, or returns None if no candidate locus aligns within `max_ed`.
pub(super) fn align_read_impl(seq: &[u8], max_ed: i32, ctx: &AlignCtx<'_>) -> Option<ShortHit> {
    let read_len = seq.len();

    if read_len < ctx.kmer_len {
        return None;
    }

    let pad = (max_ed.max(0) as usize) + WINDOW_PAD;
    // rid, rev, diagonal start.
    let (seed_hits, n_min, n_rep) = ctx.seeds.collect_with_rep(seq);
    let mut votes: Vec<(u32, bool, i64)> = seed_hits
        .iter()
        .map(|hit| (hit.rid, hit.rev != 0, hit.rpos as i64 - hit.qpos as i64))
        .collect();

    if votes.is_empty() {
        return None;
    }

    votes.sort_unstable();
    // Merge diagonals within `pad` of the cluster's first start.
    struct Locus {
        rid: u32,
        rev: bool,
        diag_lo: i64,
        diag_hi: i64,
        n_votes: usize,
    }
    let mut loci: Vec<Locus> = Vec::new();

    for &(rid, rev, diag) in &votes {
        match loci.last_mut() {
            Some(locus)
                if locus.rid == rid && locus.rev == rev && diag - locus.diag_lo <= pad as i64 =>
            {
                locus.diag_hi = diag;
                locus.n_votes += 1;
            }

            _ => loci.push(Locus {
                rid,
                rev,
                diag_lo: diag,
                diag_hi: diag,
                n_votes: 1,
            }),
        }
    }

    // Strongest loci first. The other keys only make ties deterministic.
    loci.sort_by(|lhs, rhs| {
        rhs.n_votes
            .cmp(&lhs.n_votes)
            .then(lhs.rid.cmp(&rhs.rid))
            .then(lhs.rev.cmp(&rhs.rev))
            .then(lhs.diag_lo.cmp(&rhs.diag_lo))
    });
    // The count before truncation feeds hidden_odds, since a cut list hides copies.
    let n_loci_all = loci.len() as u32;

    loci.truncate(MAX_SHORT_CANDIDATES);

    let rev_seq = revcomp(seq);
    // Verify every candidate locus up to MAX_SHORT_CANDIDATES, in order of seed support.
    // Seed support is not alignment quality.
    // At a few percent error the true copy of a repeat can keep fewer error-free minimizers than a diverged copy, so skipping weakly supported loci once a strong one aligns would pick the wrong copy at full confidence.
    // The loci share the read, so they are verified four at a time, bit-parallel with one SIMD lane per locus.
    let wins: Vec<(usize, usize)> = loci
        .iter()
        .map(|locus| {
            let rlen = ctx.ref_seqs[locus.rid as usize].len();
            // The window spans the cluster's diagonals plus pad on each side for indels.
            let win_start = (locus.diag_lo - pad as i64).max(0) as usize;
            let win_end =
                ((locus.diag_hi + read_len as i64 + pad as i64).max(0) as usize).min(rlen);

            (win_start, win_end)
        })
        .collect();
    // Two placements are distinct when their reference ends differ by more than the read-origin tolerance the truth scoring also uses.
    let sep = placement_sep(read_len) as i32;
    let mut pls = vec![
        crate::Placements {
            best: max_ed + 1,
            end: 0,
            alt: max_ed + 1,
            alt_end: 0
        };
        loci.len()
    ];

    // Loci are grouped by strand because one call aligns a single orientation of the read.
    for rev in [false, true] {
        let oriented: &[u8] = if rev { &rev_seq } else { seq };
        let idx: Vec<usize> = (0..loci.len())
            .filter(|&locus_idx| {
                loci[locus_idx].rev == rev && wins[locus_idx].0 < wins[locus_idx].1
            })
            .collect();

        for grp in idx.chunks(crate::PLACEMENT_LANES) {
            let texts: Vec<&[u8]> = grp
                .iter()
                .map(|&locus_idx| {
                    &ctx.ref_seqs[loci[locus_idx].rid as usize]
                        [wins[locus_idx].0..wins[locus_idx].1]
                })
                .collect();
            let out = crate::placements_x4(oriented, &texts, max_ed, sep);

            for (lane, &locus_idx) in grp.iter().enumerate() {
                pls[locus_idx] = out[lane];
            }
        }
    }

    // Every verified locus, in locus order.
    let found: Vec<usize> = (0..loci.len())
        .filter(|&locus_idx| pls[locus_idx].best <= max_ed)
        .collect();
    // The best placement has the lowest edit distance.
    // Equally good distinct placements, i.e. repeat copies, are chosen between pseudo-randomly, seeded by the read's sequence as BWA and minimap2 do, so a repeat family's reads spread over its copies instead of all landing on the first verified one.
    // Taking the first left 0.6% of the Salmonella LT2 genome under 10x at 30x.
    // The choice is deterministic at any thread count, and the MAPQ of a tie is 0 either way.
    let best_ed = found.iter().map(|&locus_idx| pls[locus_idx].best).min()?;
    let mut tied: Vec<usize> = Vec::new();

    for &locus_idx in &found {
        if pls[locus_idx].best != best_ed {
            continue;
        }

        let end = wins[locus_idx].0 as i64 + pls[locus_idx].end as i64;
        // Overlapping windows can report one placement twice. Such duplicates are not separate ties.
        let dup = tied.iter().any(|&tied_idx| {
            loci[tied_idx].rid == loci[locus_idx].rid
                && loci[tied_idx].rev == loci[locus_idx].rev
                && (wins[tied_idx].0 as i64 + pls[tied_idx].end as i64 - end).abs() <= sep as i64
        });

        if !dup {
            tied.push(locus_idx);
        }
    }

    // Sorting by position makes the hash pick independent of candidate order.
    tied.sort_by_key(|&locus_idx| {
        (
            loci[locus_idx].rid,
            wins[locus_idx].0 as i64 + pls[locus_idx].end as i64,
            loci[locus_idx].rev,
        )
    });
    let best_idx = if tied.len() == 1 {
        tied[0]
    } else {
        tied[(read_hash(seq) % tied.len() as u64) as usize]
    };
    let (edit_dist, rid, rev, (win_start, win_end)) = (
        pls[best_idx].best,
        loci[best_idx].rid,
        loci[best_idx].rev,
        wins[best_idx],
    );
    let best_end = win_start as i64 + pls[best_idx].end as i64;
    // Competitors are distinct placements.
    // Another window's best alignment counts unless it ends where the best one does, which is the same alignment seen from an overlapping window.
    // Every window's best alternative end counts too, because a read inside a tandem repeat aligns equally well one period over, often inside that same window.
    // Placements reported by several overlapping windows are merged by reference end.
    let mut alts: Vec<(u32, bool, i64, i32)> = Vec::new(); // (rid, rev, end, ed)

    for &locus_idx in &found {
        let placement = pls[locus_idx];
        let base = wins[locus_idx].0 as i64;
        let distinct = |end: i64| {
            loci[locus_idx].rid != rid
                || loci[locus_idx].rev != rev
                || (end - best_end).abs() > sep as i64
        };

        if locus_idx != best_idx && distinct(base + placement.end as i64) {
            alts.push((
                loci[locus_idx].rid,
                loci[locus_idx].rev,
                base + placement.end as i64,
                placement.best,
            ));
        }

        if placement.alt <= max_ed && distinct(base + placement.alt_end as i64) {
            alts.push((
                loci[locus_idx].rid,
                loci[locus_idx].rev,
                base + placement.alt_end as i64,
                placement.alt,
            ));
        }
    }

    // Sorted by place, near-identical ends sit together and merge into one competitor with the lower distance.
    alts.sort_unstable();
    let mut comp: Vec<i32> = Vec::new();
    let mut comp_loc: Vec<(u32, bool, i64)> = Vec::new();
    let mut last: Option<(u32, bool, i64)> = None;

    for &(alt_rid, alt_rev, alt_end, alt_ed) in &alts {
        match last {
            Some((last_rid, last_rev, last_end))
                if last_rid == alt_rid
                    && last_rev == alt_rev
                    && alt_end - last_end <= sep as i64 =>
            {
                let min_ed = comp.last_mut().expect("the match guard just read the last competitor");

                *min_ed = (*min_ed).min(alt_ed);
            }

            _ => {
                comp.push(alt_ed);
                comp_loc.push((alt_rid, alt_rev, alt_end));
                last = Some((alt_rid, alt_rev, alt_end));
            }
        }
    }

    // Seeds dropped as over-occurring, above the seed index's occurrence cap, are repeat copies that never became candidates.
    // A read whose minimizers are mostly repetitive, as in satellites or young repeat families, can look unique only because its competitors were filtered out, so the posterior is scaled by the fraction of the read's minimizers that still seeded.
    let rep_frac = if n_min > 0 {
        n_rep as f64 / n_min as f64
    } else {
        0.0
    };
    let n_best = loci[best_idx].n_votes as u32;
    let n_other = loci
        .iter()
        .enumerate()
        .filter(|&(locus_idx, _)| locus_idx != best_idx)
        .map(|(_, locus)| locus.n_votes as u32)
        .max()
        .unwrap_or(0);
    // gap is the edit margin to the nearest competitor, capped at MAX_GAP.
    let gap = comp
        .iter()
        .map(|&comp_ed| comp_ed - edit_dist)
        .min()
        .unwrap_or(MAX_GAP)
        .clamp(0, MAX_GAP);

    Some(ShortHit {
        edit_dist,
        comp,
        comp_loc,
        rep_frac,
        seed_ev: [n_loci_all, n_min, n_rep, n_best, n_other],
        gap,
        rid: rid as usize,
        rev,
        win_start,
        win_end,
        end: best_end,
    })
}
