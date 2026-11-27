//! Maps long reads and genome-scale queries by chaining anchors.
//! Gaps between anchors are filled with WFA or BiWFA, flanks are extended outwards, and clustered gaps are polished.
//! Chains that cover other parts of the query become supplementary records.

use super::*;

/// Reads at or above this length take the chaining path, with anchor chaining, gap-affine WFA between anchors, no -e cap, and primary plus supplementary records for query-disjoint chains, genome-scale queries included.
pub(super) const CHAIN_MIN_LEN: usize = 500;

/// Queries at least this long, such as assemblies and contigs, are mapped as genome-scale sequences, and record selection keeps every query-disjoint block.
pub(super) const GENOME_SCALE_LEN: usize = 100_000;

/// Gap-affine penalties x, o and e for chained long-read alignment.
pub(super) const CHAIN_PEN: (i32, i32, i32) = (4, 6, 2);

/// A flank whose affine score exceeds this many units per base is soft-clipped instead of aligned (random sequence scores ~3/base).
pub(super) const FLANK_MAX_PER_BASE: f64 = 1.5;

/// Flanks are extended outward from the end anchor in chunks of this many read bases, each aligned under its own FLANK_MAX_PER_BASE cap, and the first chunk over its cap soft-clips the rest in a Z-drop-style stop.
/// WFA memory grows with the square of the score, so one cap over a whole 40 kb flank meant gigabytes per thread, while a chunk cap keeps it to megabytes.
/// Chunks gallop, starting at FLANK_START_CHUNK, doubling after an aligned chunk up to FLANK_CHUNK and halving after a failed one down to FLANK_MIN_CHUNK.
/// The clip then lands near where the read leaves the reference, and a foreign tail costs only a few small capped alignments.
pub(super) const FLANK_CHUNK: usize = 1000;

pub(super) const FLANK_MIN_CHUNK: usize = 32;

pub(super) const FLANK_START_CHUNK: usize = 64;

/// A flank chunk's reference window is the chunk plus 1/FLANK_WINDOW_SLACK_DIV of it and FLANK_WINDOW_SLACK_BASES more, room for indels.
pub(super) const FLANK_WINDOW_SLACK_DIV: usize = 4;

pub(super) const FLANK_WINDOW_SLACK_BASES: usize = 32;

/// The reference window of an end flank spans this many times the read flank plus FLANK_REF_PAD bases, so the flank's indels fit.
pub(super) const FLANK_REF_SPAN_FACTOR: usize = 2;

pub(super) const FLANK_REF_PAD: usize = 64;

/// Extension scoring of one aligned column, used to cut flanks and read ends at their best-scoring prefix.
/// A match scores +2, a mismatch -4 and a gap of length len -(4 + 2 len), minimap2's defaults, with the open charged on the gap's first column.
pub(super) const EXT_MATCH: i64 = 2;

pub(super) const EXT_MISMATCH: i64 = -4;

pub(super) const EXT_GAP_OPEN: i64 = -4;

pub(super) const EXT_GAP_EXT: i64 = -2;

/// MAPQ evidence of the short-read path is put on minimap2's score scale: +MM2_MATCH_SCORE per read base and MM2_SCORE_PER_COST score units per unit of edit cost (a mismatch costs 4 here and 6 there).
pub(super) const MM2_MATCH_SCORE: f64 = 2.0;

pub(super) const MM2_SCORE_PER_COST: f64 = 1.5;

/// Extends a flank rightwards from an anchor, with `read_seq` and `ref_seq` starting right after it.
/// The flank is aligned in galloping chunks (FLANK_START_CHUNK, doubled with each aligned chunk up to FLANK_CHUNK and halved on a failed one down to FLANK_MIN_CHUNK) with the text ends free, each chunk under its own FLANK_MAX_PER_BASE cap.
/// The resulting path is then cut at its best-scoring prefix under extension scoring (match +2, mismatch -4, gap -(4 + 2*len)), the Z-drop and X-drop rule of minimap2 and BLAST, so bases past the point where the read leaves the reference are soft-clipped rather than aligned as noise.
/// Per-chunk caps bound the wavefront memory, since a single cap over a 40 kb flank meant GBs per thread.
/// Returns (ops, read bases used, reference bases used, chain-penalty cost including FLANK_MAX_PER_BASE per clipped base).
pub(super) fn extend_flank(read_seq: &[u8], ref_seq: &[u8]) -> (Vec<(u32, u8)>, usize, usize, i64) {
    let mut ops: Vec<(u32, u8)> = Vec::new();
    let (mut read_used, mut ref_used) = (0usize, 0usize);
    // The chunk size gallops, starting small, doubling after each aligned chunk and halving after a failed one.
    // A foreign tail then costs a few tiny capped alignments, while a long flank that truly matches reaches FLANK_CHUNK in four steps.
    const CHUNK_GROWTH: usize = 2;
    let mut chunk_len = FLANK_START_CHUNK;

    while read_used < read_seq.len() && ref_used < ref_seq.len() {
        // Near the contig end, a chunk longer than the remaining reference would have to absorb the overhang as an insertion and fail its cap, clipping end bases that do match.
        chunk_len = chunk_len
            .min(read_seq.len() - read_used)
            .min(ref_seq.len() - ref_used);
        // The text window covers the chunk plus 25% and 32 bases for indels.
        let span = (chunk_len + chunk_len / FLANK_WINDOW_SLACK_DIV + FLANK_WINDOW_SLACK_BASES)
            .min(ref_seq.len() - ref_used);
        let text = &ref_seq[ref_used..ref_used + span];

        match crate::wfa_affine_align(
            &read_seq[read_used..read_used + chunk_len],
            text,
            text.len(),
            CHAIN_PEN,
            (chunk_len as f64 * FLANK_MAX_PER_BASE) as i32,
        ) {
            Some(aln) => {
                push_cigar(&mut ops, &aln.cigar);
                ref_used += aln.tend;
                read_used += chunk_len;
                chunk_len = (chunk_len * CHUNK_GROWTH).min(FLANK_CHUNK);
            }
            None if chunk_len > FLANK_MIN_CHUNK => chunk_len /= CHUNK_GROWTH,
            None => break,
        }
    }

    // Exact matches past the last chunk, since a chunk that runs into an N gap or a foreign tail fails as a whole and drops its matching first bases.
    let mut run = 0u32;

    while read_used < read_seq.len()
        && ref_used < ref_seq.len()
        && ref_seq[ref_used] != b'N'
        && read_seq[read_used].eq_ignore_ascii_case(&ref_seq[ref_used])
    {
        read_used += 1;
        ref_used += 1;
        run += 1;
    }

    if run > 0 {
        push_op(&mut ops, run, b'M');
    }

    // Best-scoring prefix of the path.
    let same =
        |read_base: u8, ref_base: u8| read_base.eq_ignore_ascii_case(&ref_base) && ref_base != b'N';
    let (mut score, mut best) = (0i64, 0i64);
    let mut cut = (0usize, 0u32, 0usize, 0usize); // op index, bases of it, q, r
    let (mut read_i, mut ref_i) = (0usize, 0usize);

    for (i, &(run_len, cigar_op)) in ops.iter().enumerate() {
        match cigar_op {
            b'M' => {
                for t in 0..run_len {
                    score += if same(read_seq[read_i], ref_seq[ref_i]) {
                        EXT_MATCH
                    } else {
                        EXT_MISMATCH
                    };
                    read_i += 1;
                    ref_i += 1;

                    // The best prefix always ends on a match, since only a matching base raises score.
                    if score > best {
                        best = score;
                        cut = (i, t + 1, read_i, ref_i);
                    }
                }
            }
            b'I' => {
                score += EXT_GAP_OPEN + EXT_GAP_EXT * run_len as i64;
                read_i += run_len as usize;
            }
            _ => {
                score += EXT_GAP_OPEN + EXT_GAP_EXT * run_len as i64;
                ref_i += run_len as usize;
            }
        }
    }

    let (cut_op, cut_len, read_cut, ref_cut) = cut;
    let mut kept: Vec<(u32, u8)> = Vec::new();

    if best > 0 {
        kept.extend_from_slice(&ops[..cut_op]);
        kept.push((cut_len, b'M'));
    }

    // Chain-penalty cost of the kept path plus the clip.
    let (mut cost, mut read_i, mut ref_i) = (0i64, 0usize, 0usize);

    for &(run_len, cigar_op) in &kept {
        match cigar_op {
            b'M' => {
                for _ in 0..run_len {
                    if !same(read_seq[read_i], ref_seq[ref_i]) {
                        cost += CHAIN_PEN.0 as i64;
                    }

                    read_i += 1;
                    ref_i += 1;
                }
            }
            b'I' => {
                cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * run_len as i32) as i64;
                read_i += run_len as usize;
            }
            _ => {
                cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * run_len as i32) as i64;
                ref_i += run_len as usize;
            }
        }
    }

    // Each clipped read base is charged FLANK_MAX_PER_BASE, the most an aligned base may cost. Clipping is therefore never cheaper than a bad alignment.
    cost += ((read_seq.len() - read_cut) as f64 * FLANK_MAX_PER_BASE) as i64;
    (kept, read_cut, ref_cut, cost)
}

/// Windows longer than this many columns are not re-aligned by [`polish_gap_clusters`].
pub(super) const POLISH_MAX_WINDOW: usize = 4000;

/// Gap runs shorter than this many columns do not start a cluster.
pub(super) const POLISH_MIN_RUN: usize = 3;

/// Columns of exact context on each side of a re-aligned window.
pub(super) const POLISH_MARGIN: usize = 16;

/// Gap runs up to POLISH_SHORT_GAP columns join a cluster when the next one starts within max(POLISH_SHORT_GAP, POLISH_SHORT_GAP_REACH x gap).
/// Longer gap runs join within POLISH_LONG_GAP_REACH x gap instead.
pub(super) const POLISH_SHORT_GAP: usize = 64;

pub(super) const POLISH_SHORT_GAP_REACH: usize = 8;

pub(super) const POLISH_LONG_GAP_REACH: usize = 2;

/// Re-aligns clusters of nearby gaps freely.
/// A chain's alignment is forced through its exact anchors, and inside a tandem repeat the anchors can sit on a neighbouring diagonal, so the path contains indels that cancel.
/// Examples are 9D 8D 27I where 1D 9I 1I 1I fits, or a 497 bp D and I one repeat unit apart where the copy aligns with mismatches.
/// Gap runs closer than max(64, 8 x the larger gap), or 2 x for gaps over 64, form a cluster, since a 20D and a 20I 73 bp apart cost 92 where 1I 1D and 4 mismatches cost 36.
/// Windows over POLISH_MAX_WINDOW columns are left alone, because they lie in satellite arrays, where aligning them doubled the HiFi mapping time on hg38 for placements that stay ambiguous.
/// Each cluster's window, with POLISH_MARGIN columns of margin, is aligned end to end with gap-affine WFA and replaced when that scores lower.
/// Only gap runs of at least POLISH_MIN_RUN columns form clusters, and gaps inside a window are re-aligned with it.
/// The 1-2 bp indels of a noisy read lie a few dozen bases apart, and joining them re-aligned whole ONT reads in 4 kb windows at 2.3x their mapping time.
/// Returns the new ops and the cost saved.
pub(super) fn polish_gap_clusters(
    ops: &[(u32, u8)],
    read_seq: &[u8],
    ref_seq: &[u8],
    read_start: usize,
    ref_start: usize,
) -> (Vec<(u32, u8)>, i64) {
    // The ops are expanded to one byte per column, so windows can be cut and spliced by index.
    let mut cols: Vec<u8> = Vec::new();

    for &(run_len, cigar_op) in ops {
        cols.extend(std::iter::repeat_n(cigar_op, run_len as usize));
    }

    let num_cols = cols.len();
    // read_off[c] and ref_off[c] are the read and reference offsets at which column c starts.
    let (mut read_off, mut ref_off) = (vec![0usize; num_cols + 1], vec![0usize; num_cols + 1]);

    for i in 0..num_cols {
        read_off[i + 1] = read_off[i] + (cols[i] != b'D') as usize;
        ref_off[i + 1] = ref_off[i] + (cols[i] != b'I') as usize;
    }

    if read_start + read_off[num_cols] > read_seq.len()
        || ref_start + ref_off[num_cols] > ref_seq.len()
    {
        return (ops.to_vec(), 0);
    }

    // Only runs of one gap type that reach POLISH_MIN_RUN columns can seed a cluster.
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut col = 0;

    while col < num_cols {
        if cols[col] == b'M' {
            col += 1;
            continue;
        }

        let mut run_end = col;

        while run_end < num_cols && cols[run_end] == cols[col] {
            run_end += 1;
        }

        if run_end - col >= POLISH_MIN_RUN {
            runs.push((col, run_end));
        }

        col = run_end;
    }

    // The window's present cost uses the same penalties as the WFA call below, so the two scores compare directly.
    let base_ok = |base: u8| matches!(base, b'A' | b'C' | b'G' | b'T');
    let (mismatch_pen, gap_open_pen, gap_ext_pen) =
        (CHAIN_PEN.0 as i64, CHAIN_PEN.1 as i64, CHAIN_PEN.2 as i64);
    let mut out: Vec<u8> = Vec::with_capacity(num_cols);
    let (mut done, mut saved, mut run_idx) = (0usize, 0i64, 0usize);

    while run_idx < runs.len() {
        let (cluster_start, mut cluster_end) = runs[run_idx];
        let (mut last_run, mut maxg) = (run_idx, cluster_end - cluster_start);

        // The cluster grows while the next run starts within reach of the largest gap seen so far.
        while last_run + 1 < runs.len() {
            let (next_start, next_end) = runs[last_run + 1];
            let gap_len = maxg.max(next_end - next_start);

            if next_start - cluster_end
                > if gap_len <= POLISH_SHORT_GAP {
                    POLISH_SHORT_GAP.max(POLISH_SHORT_GAP_REACH * gap_len)
                } else {
                    POLISH_LONG_GAP_REACH * gap_len
                }
            {
                break;
            }
            last_run += 1;
            cluster_end = next_end;
            maxg = gap_len;
        }

        // Starting no earlier than `done` keeps this window clear of one that was already replaced.
        let (win_start, win_end) = (
            cluster_start.saturating_sub(POLISH_MARGIN).max(done),
            (cluster_end + POLISH_MARGIN).min(num_cols),
        );

        // A lone gap run has nothing to cancel against, so only clusters of two or more runs are re-aligned.
        // cost is the window's current price under CHAIN_PEN, which the re-alignment has to beat.
        if last_run > run_idx && win_end - win_start <= POLISH_MAX_WINDOW {
            let mut cost = 0i64;

            for c in win_start..win_end {
                match cols[c] {
                    b'M' => {
                        let (read_base, ref_base) = (
                            read_seq[read_start + read_off[c]],
                            ref_seq[ref_start + ref_off[c]],
                        );

                        if read_base != ref_base || !base_ok(read_base) || !base_ok(ref_base) {
                            cost += mismatch_pen;
                        }
                    }
                    gap_op => {
                        cost += gap_ext_pen;

                        if c == win_start || cols[c - 1] != gap_op {
                            cost += gap_open_pen;
                        }
                    }
                }
            }

            let (qsub, rsub) = (
                &read_seq[read_start + read_off[win_start]..read_start + read_off[win_end]],
                &ref_seq[ref_start + ref_off[win_start]..ref_start + ref_off[win_end]],
            );

            if !qsub.is_empty() && !rsub.is_empty() {
                // Capping WFA one below the present cost lets it give up as soon as it cannot improve on it.
                // Long windows go through BiWFA, whose memory grows linearly with the window instead of with its square.
                let cap = (cost - 1).min(i32::MAX as i64) as i32;
                let aln = if qsub.len().max(rsub.len()) > BIWFA_MIN_GAP {
                    crate::biwfa_align(qsub, rsub, CHAIN_PEN, cap)
                } else {
                    crate::wfa_affine_align(qsub, rsub, 0, CHAIN_PEN, cap)
                };

                if let Some(better) = aln.filter(|cand| (cand.score as i64) < cost) {
                    out.extend_from_slice(&cols[done..win_start]);
                    let mut run_len = 0usize;

                    // The WFA CIGAR distinguishes = and X, and both fold back into M columns here.
                    for byte in better.cigar.bytes() {
                        if byte.is_ascii_digit() {
                            run_len = push_digit(run_len, byte);
                            continue;
                        }

                        let cigar_op = if matches!(byte, b'=' | b'X') {
                            b'M'
                        } else {
                            byte
                        };

                        out.extend(std::iter::repeat_n(cigar_op, run_len));
                        run_len = 0;
                    }

                    // Moving done to win_end makes later windows start after the replaced columns.
                    saved += cost - better.score as i64;
                    done = win_end;
                }
            }
        }

        run_idx = last_run + 1;
    }

    if saved == 0 {
        return (ops.to_vec(), 0);
    }

    out.extend_from_slice(&cols[done..]);
    let mut res: Vec<(u32, u8)> = Vec::new();

    for &col_op in &out {
        push_op(&mut res, 1, col_op);
    }

    (res, saved)
}

/// minimap2's alignment scoring from v2.31 options.c, with `match_score` (minimap2's a), `mismatch` (b), two-piece affine gaps min(gap_open + gap_ext len, gap_open2 + gap_ext2 len), an N scoring -sc_ambi = -1, and the Z-drop threshold.
/// ksw2 counts a drop past zdrop + gap_ext2 per column of diagonal shift since the best score.
#[derive(Clone, Copy)]
pub(super) struct MmScore {
    pub(super) match_score: i64,
    pub(super) mismatch: i64,
    pub(super) gap_open: i64,
    pub(super) gap_ext: i64,
    pub(super) gap_open2: i64,
    pub(super) gap_ext2: i64,
    pub(super) zdrop: i64,
}

/// mm_mapopt_init for reads, where map-ont is the default.
pub(super) const MM_READS: MmScore = MmScore {
    match_score: 2,
    mismatch: 4,
    gap_open: 4,
    gap_ext: 2,
    gap_open2: 24,
    gap_ext2: 1,
    zdrop: 400,
};

/// -x asm20, used for genome-scale queries.
pub(super) const MM_ASM20: MmScore = MmScore {
    match_score: 1,
    mismatch: 4,
    gap_open: 6,
    gap_ext: 2,
    gap_open2: 26,
    gap_ext2: 1,
    zdrop: 200,
};

impl MmScore {
    /// Two-piece affine cost of a gap of `len`.
    pub(super) fn gap(&self, len: i64) -> i64 {
        (self.gap_open + self.gap_ext * len).min(self.gap_open2 + self.gap_ext2 * len)
    }
    /// Score of one aligned column, where an N scores -1.
    pub(super) fn base(&self, read_base: u8, ref_base: u8) -> i64 {
        /// Score of a column with an N on either side (minimap2 -sc_ambi = 1).
        const AMBIGUOUS_BASE_SCORE: i64 = -1;

        if read_base == b'N' || ref_base == b'N' {
            AMBIGUOUS_BASE_SCORE
        } else if read_base.eq_ignore_ascii_case(&ref_base) {
            self.match_score
        } else {
            -self.mismatch
        }
    }
}

/// Local clipping of an aligned path, as minimap2 treats a chain.
/// Returns the column ranges [start, end) worth keeping as records with their scores under `scoring`, minimap2's own scoring, where a gap column scores its increment of the two-piece gap cost, counted from the gap's first column in scan order.
/// The path is split wherever its score falls more than scoring.zdrop plus scoring.gap_ext2 per column of diagonal shift below the best so far.
/// Under this ksw2 Z-drop, one long indel never splits.
/// A divergent stretch inside the path, meaning anchors across inserted or rearranged sequence whose gaps are short enough to align uncapped, becomes a break between two records.
/// Each piece's ends are then cut at their best-scoring prefix and suffix, but never inward of its outermost anchors.
/// Between anchors the path stays global, as minimap2 fills a chain's gaps, and only what lies beyond them is clipped, as its end extension does.
/// A piece needs an anchor of its own, found through `anchor_q`, the ascending read offsets of the segment's anchors, and a score of at least one anchor, k a.
/// `read_off` and `ref_off` are the read and reference offsets before each column.
#[allow(clippy::too_many_arguments)]
pub(super) fn local_pieces(
    cols: &[u8],
    read_seq: &[u8],
    ref_seq: &[u8],
    read_off: &[usize],
    ref_off: &[usize],
    anchor_q: &[usize],
    kmer_len: usize,
    scoring: MmScore,
) -> Vec<(usize, usize, i64)> {
    let num_cols = cols.len();
    // Position of each gap column in its run, counted from either side with 1 as the first.
    let (mut gap_pos_fwd, mut gap_pos_back) = (vec![0i64; num_cols], vec![0i64; num_cols]);

    for c in 0..num_cols {
        if cols[c] != b'M' {
            gap_pos_fwd[c] = if c > 0 && cols[c - 1] == cols[c] {
                gap_pos_fwd[c - 1] + 1
            } else {
                1
            };
        }
    }

    for c in (0..num_cols).rev() {
        if cols[c] != b'M' {
            gap_pos_back[c] = if c + 1 < num_cols && cols[c + 1] == cols[c] {
                gap_pos_back[c + 1] + 1
            } else {
                1
            };
        }
    }

    let score = |col: usize, fwd: bool| -> i64 {
        if cols[col] == b'M' {
            return scoring.base(read_seq[read_off[col]], ref_seq[ref_off[col]]);
        }

        let gap_pos = if fwd {
            gap_pos_fwd[col]
        } else {
            gap_pos_back[col]
        };

        // A gap column scores gap(n) - gap(n - 1). Under minimap2's two-piece gaps the marginal cost depends on the position in the run.
        scoring.gap(gap_pos - 1) - scoring.gap(gap_pos)
    };
    // 1. Z-drop splits.
    let mut raw: Vec<(usize, usize)> = Vec::new();
    let (mut start, mut acc, mut best, mut best_end) = (0usize, 0i64, 0i64, 0usize);

    for j in 0..num_cols {
        acc += score(j, true);

        if acc > best {
            best = acc;
            best_end = j + 1;
        } else {
            // As in ksw2, the Z-drop allowance grows by gap_ext2 per diagonal of shift since the peak.
            let shift = (read_off[j + 1] - read_off[best_end])
                .abs_diff(ref_off[j + 1] - ref_off[best_end]) as i64;

            if best - acc > scoring.zdrop + scoring.gap_ext2 * shift {
                if best_end > start {
                    raw.push((start, best_end));
                }

                start = j + 1;
                acc = 0;
                best = 0;
                best_end = j + 1;
            }
        }
    }

    if num_cols > start {
        raw.push((start, num_cols));
    }

    // 2. Ends are cut at their best score, outside the outermost anchors.
    let col_of = |qoff: usize| read_off.partition_point(|&off| off < qoff);
    let mut out = Vec::new();

    for (piece_start, piece_end) in raw {
        let inside: Vec<usize> = anchor_q
            .iter()
            .copied()
            .filter(|&anchor_off| {
                anchor_off >= read_off[piece_start] && anchor_off + kmer_len <= read_off[piece_end]
            })
            .collect();
        let (Some(&first_anchor), Some(&last_anchor)) = (inside.first(), inside.last()) else {
            continue;
        };
        let (left_bound, right_bound) = (
            col_of(first_anchor).max(piece_start),
            col_of(last_anchor + kmer_len).min(piece_end),
        );
        // Right end, the best prefix score over [hi, e].
        let mut acc: i64 = (piece_start..right_bound).map(|col| score(col, true)).sum();
        let (mut end_col, mut end_score) = (right_bound, acc);

        for c in right_bound..piece_end {
            acc += score(c, true);

            if acc > end_score {
                end_score = acc;
                end_col = c + 1;
            }
        }

        // Left start, the best suffix score ending at ee over [s, lo].
        let mut acc: i64 = (left_bound..end_col)
            .rev()
            .map(|col| score(col, false))
            .sum();
        let (mut start_col, mut start_score) = (left_bound, acc);

        for c in (piece_start..left_bound).rev() {
            acc += score(c, false);

            if acc > start_score {
                start_score = acc;
                start_col = c;
            }
        }

        if start_score >= kmer_len as i64 * scoring.match_score {
            out.push((start_col, end_col, start_score));
        }
    }

    out
}

/// Genome-scale records scoring below this under MM_ASM20 are dropped, as for minimap2 -x asm20's minimal peak DP score (v2.31 options.c min_dp_max = 200, where hit.c drops a record whose dp_max is lower).
pub(super) const ASM_MIN_DP: i64 = 200;

/// A record's alignment score under minimap2 -x asm20's scoring, MM_ASM20.
/// local_pieces cuts records at their best-scoring ends with the same scoring, so the whole record's score stands for minimap2's peak, dp_max.
/// At 1 and -4 a column scores above zero only below 20% divergence.
/// Diverged homologous sequence, such as different genes of two strains, therefore scores under ASM_MIN_DP, and minimap2 -x asm20 reports no record for it.
pub(super) fn asm_score(cigar: &str, read_seq: &[u8], ref_seq: &[u8], pos: usize) -> i64 {
    let scoring = MM_ASM20;
    let (mut read_i, mut ref_i, mut run_len, mut tot) = (0usize, pos, 0usize, 0i64);

    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        match byte {
            // Clips score nothing. read_seq is the full read, which is why hard clips still advance read_i.
            b'S' | b'H' => read_i += run_len,
            b'M' | b'=' | b'X' => {
                for t in 0..run_len {
                    tot += scoring.base(read_seq[read_i + t], ref_seq[ref_i + t]);
                }

                read_i += run_len;
                ref_i += run_len;
            }
            b'I' => {
                tot -= scoring.gap(run_len as i64);
                read_i += run_len;
            }
            b'D' => {
                tot -= scoring.gap(run_len as i64);
                ref_i += run_len;
            }
            _ => {}
        }

        run_len = 0;
    }

    tot
}

/// One aligned piece of a chain, holding the SAM fields of a record.
pub(super) struct ChainAln {
    pub(super) rid: usize,
    pub(super) rev: bool,
    pub(super) pos: usize,
    pub(super) cigar: String,
    pub(super) edit_dist: i32,
    pub(super) read_start: usize,
    pub(super) read_end: usize,
    pub(super) n_anchors: usize,
    /// Gap-affine cost of the base-level alignment, clipped bases included.
    pub(super) cost: i64,
}

/// Gaps between anchors up to this long on both sides use the traceback WFA, whose score cap of 1.5 per base keeps its O(s^2) memory near 10 MB.
/// Longer gaps use BiWFA, with O(s) memory and the same optimum.
pub(super) const BIWFA_MIN_GAP: usize = 1000;

/// Gaps up to this long are aligned without a score cap and never split a chain (worst case s ~ 4 * 500, a few MB of wavefronts).
pub(super) const SPLIT_MIN_GAP: usize = 500;

/// A gap whose shorter side has at most this many bases, against a much longer side, is aligned as one long indel.
pub(super) const TINY_GAP_SIDE: usize = 32;

/// The tiny side is aligned against this many bases more than its own length at the near end of the long side.
pub(super) const TINY_GAP_WINDOW_PAD: usize = 64;

/// Keeps strictly non-overlapping anchors that lie inside both sequences, each contributing k exact matches.
fn plausible_anchors(chain: &crate::chain::Chain, read_seq: &[u8], rseq: &[u8], kmer_len: usize) -> Vec<(usize, usize)> {
    let mut anchors: Vec<(usize, usize)> = Vec::with_capacity(chain.anchors.len());

    for anchor in &chain.anchors {
        let (ref_pos, read_pos) = (anchor.ref_pos as usize, anchor.read_pos as usize);

        if let Some(&(prev_ref, prev_read)) = anchors.last()
            && (ref_pos < prev_ref + kmer_len || read_pos < prev_read + kmer_len)
        {
            continue;
        }

        if ref_pos + kmer_len > rseq.len() || read_pos + kmer_len > read_seq.len() {
            continue;
        }

        anchors.push((ref_pos, read_pos));
    }

    anchors
}

/// A run of anchors aligned end to end, with its path and gap-affine cost.
/// align_chain splits a chain into several of these wherever a gap between anchors is implausible.
struct AnchorSegment {
    first: usize,
    last: usize,
    ops: Vec<(u32, u8)>,
    cost: i64,
}

/// Aligns the gaps between consecutive anchors and starts a new segment wherever a gap fails its cost cap.
fn align_anchor_runs(anchors: &[(usize, usize)], read_seq: &[u8], rseq: &[u8], kmer_len: usize) -> Vec<AnchorSegment> {
    let mut segs: Vec<AnchorSegment> = Vec::new();
    let mut cur = AnchorSegment {
        first: 0,
        last: 0,
        ops: Vec::new(),
        cost: 0,
    };

    for w in 0..anchors.len() {
        let (ref_pos, read_pos) = anchors[w];

        // Each anchor is an exact k-mer match, so it contributes kmer_len M columns.
        push_op(&mut cur.ops, kmer_len as u32, b'M');
        cur.last = w;

        if w + 1 == anchors.len() {
            break;
        }

        let (next_ref, next_read) = anchors[w + 1];
        let (gap_read, gap_ref) = (
            &read_seq[read_pos + kmer_len..next_read],
            &rseq[ref_pos + kmer_len..next_ref],
        );

        // When only one side of the gap has bases, it is a pure insertion or deletion and needs no WFA call.
        match (gap_read.is_empty(), gap_ref.is_empty()) {
            (true, true) => {}
            (false, true) => {
                push_op(&mut cur.ops, gap_read.len() as u32, b'I');
                cur.cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * gap_read.len() as i32) as i64;
            }
            (true, false) => {
                push_op(&mut cur.ops, gap_ref.len() as u32, b'D');
                cur.cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * gap_ref.len() as i32) as i64;
            }
            _ if gap_read.len().min(gap_ref.len()) <= TINY_GAP_SIDE
                && gap_read.len().abs_diff(gap_ref.len()) > SPLIT_MIN_GAP =>
            {
                // One tiny side against a long one is a long deletion or insertion.
                // The short side is aligned ends-free against the near end of the long side, and the rest becomes one gap that left_align places.
                let (short, long, del) = if gap_read.len() < gap_ref.len() {
                    (gap_read, gap_ref, true)
                } else {
                    (gap_ref, gap_read, false)
                };
                let window_len = (short.len() + TINY_GAP_WINDOW_PAD).min(long.len());
                let aln = crate::wfa_affine_align(
                    short,
                    &long[..window_len],
                    window_len,
                    CHAIN_PEN,
                    NO_SCORE_CAP,
                );
                let (cig, tend, gap_score) = match aln {
                    Some(aln_res) => (aln_res.cigar, aln_res.tend, aln_res.score),
                    None => (String::new(), 0, 0),
                };
                let mut tmp: Vec<(u32, u8)> = Vec::new();

                if tend == 0 && !short.is_empty() {
                    // Fallback, with the short side as its own gap.
                    push_op(&mut tmp, short.len() as u32, if del { b'I' } else { b'D' });
                } else {
                    push_cigar(&mut tmp, &cig);

                    if !del {
                        // Aligned with the read's long side as text, so I and D swap.
                        for op_run in tmp.iter_mut() {
                            op_run.1 = match op_run.1 {
                                b'I' => b'D',
                                b'D' => b'I',
                                other => other,
                            };
                        }
                    }
                }

                for (run_len, cigar_op) in tmp {
                    push_op(&mut cur.ops, run_len, cigar_op);
                }

                let rest = long.len() - if tend == 0 { 0 } else { tend };

                push_op(&mut cur.ops, rest as u32, if del { b'D' } else { b'I' });
                cur.cost += gap_score as i64 + (CHAIN_PEN.1 + CHAIN_PEN.2 * rest as i32) as i64;
            }
            _ => {
                // Short gaps are always aligned, exactly and in tiny memory.
                // Their per-base cost is noisy, since a few indels in a 20 bp gap of a 10%-error read exceed any per-base cap, and splitting there would cut correct alignments.
                // Longer gaps are capped and split when the cap fails.
                // The cap allows the length difference as one indel, open plus extend per base, plus FLANK_MAX_PER_BASE over the part both sides share.
                // A 600 bp insertion with 80 reference bases around it then stays an insertion, written as one I in one record as minimap2 and SV callers expect, while two long sides that do not align still fail and split.
                // A flat 1.5 per base over the longer side split every insertion or deletion over about 500 bp, and on HG002 HiFi it moved 77 reads' primaries to the far side of a real insertion.
                let glen = gap_read.len().max(gap_ref.len());
                let (shared, diff) = (
                    gap_read.len().min(gap_ref.len()),
                    gap_read.len().abs_diff(gap_ref.len()),
                );
                let cap = if glen <= SPLIT_MIN_GAP {
                    NO_SCORE_CAP
                } else {
                    CHAIN_PEN.1
                        + CHAIN_PEN.2 * diff as i32
                        + (shared as f64 * FLANK_MAX_PER_BASE) as i32
                        + CHAIN_PEN.1
                };
                let aln = if glen > BIWFA_MIN_GAP {
                    crate::biwfa_align(gap_read, gap_ref, CHAIN_PEN, cap)
                } else {
                    crate::wfa_affine_align(gap_read, gap_ref, 0, CHAIN_PEN, cap)
                };

                // A gap that fails its cap closes the current segment, and the next one starts at the following anchor.
                match aln {
                    Some(aln) => {
                        push_cigar(&mut cur.ops, &aln.cigar);
                        cur.cost += aln.score as i64;
                    }
                    None => {
                        segs.push(std::mem::replace(
                            &mut cur,
                            AnchorSegment {
                                first: w + 1,
                                last: w + 1,
                                ops: Vec::new(),
                                cost: 0,
                            },
                        ));
                    }
                }
            }
        }
    }

    segs.push(cur);

    segs
}

/// Extends segment `seg_idx` into its flanks, up to the neighbouring segments, and clips it into its best-scoring pieces.
/// Each piece becomes one record, and a segment whose path no longer fits the sequences yields none.
#[allow(clippy::too_many_arguments)]
fn finish_segment(
    chain: &crate::chain::Chain,
    seg_idx: usize,
    segs: &[AnchorSegment],
    anchors: &[(usize, usize)],
    read_seq: &[u8],
    rseq: &[u8],
    kmer_len: usize,
    scoring: MmScore,
) -> Vec<ChainAln> {
    let seg = &segs[seg_idx];
    let i = seg_idx;
    let mut out = Vec::new();
    let (first_ref, first_read) = anchors[seg.first];
    let (last_ref, last_read) = anchors[seg.last];
    // Flank regions end at the neighbouring segment's anchor.
    let (q_lo, r_lo) = if i == 0 {
        (0, 0)
    } else {
        let (prev_ref, prev_read) = anchors[segs[i - 1].last];

        (prev_read + kmer_len, prev_ref + kmer_len)
    };
    let (q_hi, r_hi) = if i + 1 == segs.len() {
        (read_seq.len(), rseq.len())
    } else {
        let (next_ref, next_read) = anchors[segs[i + 1].first];

        (next_read, next_ref)
    };
    let mut ops: Vec<(u32, u8)> = Vec::new();
    let mut cost = seg.cost;
    // The left flank extends leftwards, i.e. rightwards over the reversed sequences.
    let left_len = first_read - q_lo;
    let ql_rev: Vec<u8> = read_seq[q_lo..first_read].iter().rev().copied().collect();
    let rl_rev: Vec<u8> = rseq[first_ref
        .saturating_sub(FLANK_REF_SPAN_FACTOR * left_len + FLANK_REF_PAD)
        .max(r_lo)..first_ref]
        .iter()
        .rev()
        .copied()
        .collect();
    let (fops, read_used, ref_used, fcost) = extend_flank(&ql_rev, &rl_rev);

    for &(run_len, cigar_op) in fops.iter().rev() {
        push_op(&mut ops, run_len, cigar_op);
    }

    cost += fcost;
    let (ref_start, clip_l) = (first_ref - ref_used, first_read - read_used);

    for &(run_len, cigar_op) in &seg.ops {
        push_op(&mut ops, run_len, cigar_op);
    }

    // Right flank.
    let (right_read, right_ref) = (last_read + kmer_len, last_ref + kmer_len);
    let right_len = q_hi - right_read;
    let (fops, read_used, _, fcost) = extend_flank(
        &read_seq[right_read..q_hi],
        &rseq[right_ref
            ..(right_ref + FLANK_REF_SPAN_FACTOR * right_len + FLANK_REF_PAD).min(r_hi)],
    );

    for (run_len, cigar_op) in fops {
        push_op(&mut ops, run_len, cigar_op);
    }

    cost += fcost;
    let clip_r = read_seq.len() - (right_read + read_used);
    let (ops, saved) = polish_gap_clusters(&ops, read_seq, rseq, clip_l, ref_start);

    cost -= saved;
    // Local clipping keeps the path's best-scoring pieces (local_pieces).
    let mut cols: Vec<u8> = Vec::new();

    for &(run_len, cigar_op) in &ops {
        cols.extend(std::iter::repeat_n(cigar_op, run_len as usize));
    }

    let num_cols = cols.len();
    let (mut read_off, mut ref_off) =
        (vec![clip_l; num_cols + 1], vec![ref_start; num_cols + 1]);

    for i in 0..num_cols {
        read_off[i + 1] = read_off[i] + (cols[i] != b'D') as usize;
        ref_off[i + 1] = ref_off[i] + (cols[i] != b'I') as usize;
    }

    // A path whose read span disagrees with the clips, or that runs off the reference, is dropped rather than reported.
    if read_off[num_cols] != read_seq.len() - clip_r || ref_off[num_cols] > rseq.len() {
        return out;
    }

    let anchors = &anchors[seg.first..=seg.last];
    let anchor_q: Vec<usize> = anchors
        .iter()
        .map(|&(_, anchor_read)| anchor_read)
        .collect();
    let pieces: Vec<(usize, usize)> = local_pieces(
        &cols, read_seq, rseq, &read_off, &ref_off, &anchor_q, kmer_len, scoring,
    )
    .into_iter()
    .map(|(piece_start, piece_end, _)| (piece_start, piece_end))
    .collect();
    // With no piece worth keeping, the whole path is kept as one piece, so the segment still gets a record.
    let pieces = if pieces.is_empty() {
        vec![(0, num_cols)]
    } else {
        pieces
    };

    // Columns outside every piece have their path cost replaced by the clip charge per read base, as extend_flank charges a clipped flank.
    if pieces != [(0, num_cols)] {
        let mut keep = vec![false; num_cols];

        for &(piece_start, piece_end) in &pieces {
            keep[piece_start..piece_end].fill(true);
        }

        let mut dropped = 0usize;

        for i in (0..num_cols).filter(|&i| !keep[i]) {
            cost -= match cols[i] {
                b'M' => {
                    if read_seq[read_off[i]].eq_ignore_ascii_case(&rseq[ref_off[i]])
                        && rseq[ref_off[i]] != b'N'
                    {
                        0
                    } else {
                        CHAIN_PEN.0 as i64
                    }
                }
                gap_op => {
                    CHAIN_PEN.2 as i64
                        + if i == 0 || cols[i - 1] != gap_op {
                            CHAIN_PEN.1 as i64
                        } else {
                            0
                        }
                }
            };
            // D columns consume no read base and add nothing to the clip charge.
            dropped += (cols[i] != b'D') as usize;
        }

        cost += (dropped as f64 * FLANK_MAX_PER_BASE) as i64;
    }

    // Each piece becomes one record, left-aligned and soft-clipped to the full read.
    for (piece_idx, &(piece_start, piece_end)) in pieces.iter().enumerate() {
        let (qs0, qe0, rs0) = (
            read_off[piece_start],
            read_off[piece_end],
            ref_off[piece_start],
        );
        let mut pops: Vec<(u32, u8)> = Vec::new();

        for &col_op in &cols[piece_start..piece_end] {
            push_op(&mut pops, 1, col_op);
        }

        let core: String = pops
            .iter()
            .map(|(run_len, cigar_op)| format!("{run_len}{}", *cigar_op as char))
            .collect();
        let core = left_align(&core, &read_seq[qs0..qe0], rseq, rs0);
        // cigar_ed_implied() fails only on inconsistent bookkeeping, in which case the records already built are returned.
        let Some(edit_dist) = cigar_ed_implied(&core, &read_seq[qs0..qe0], rseq, rs0) else {
            return out;
        };
        let mut cigar = String::new();

        if qs0 > 0 {
            cigar.push_str(&format!("{qs0}S"));
        }

        cigar.push_str(&core);

        if qe0 < read_seq.len() {
            cigar.push_str(&format!("{}S", read_seq.len() - qe0));
        }

        // Only anchors lying wholly inside the piece count toward its support.
        let n_anchors = anchors
            .iter()
            .filter(|&&(_, anchor_read)| anchor_read >= qs0 && anchor_read + kmer_len <= qe0)
            .count();

        out.push(ChainAln {
            rid: chain.rid as usize,
            rev: chain.rev,
            pos: rs0,
            cigar,
            edit_dist,
            read_start: qs0,
            read_end: qe0,
            n_anchors,
            // The segment's cost is charged once, to its first piece, so summing over records does not count it twice.
            cost: if piece_idx == 0 { cost } else { 0 },
        });
    }

    out
}

/// Aligns a chain with exact anchors, gap-affine WFA between consecutive anchors end to end under a FLANK_MAX_PER_BASE per-base cap, and flanks extended outward by `extend_flank`.
/// A gap that fails its cap is not a plausible alignment, because the sequence is unrelated or rearranged, so the chain is split there.
/// Each side becomes its own record that extends into the gap and soft-clips the rest, as minimap2 splits at a Z-drop.
/// `local_pieces` then clips each segment's path, cutting its ends at their best score and splitting where a divergent stretch drops the score by more than ZDROP.
/// Returns the records in read order, or nothing if every anchor is filtered out.
pub(super) fn align_chain(
    chain: &crate::chain::Chain,
    read_seq: &[u8],
    rseq: &[u8],
    kmer_len: usize,
    scoring: MmScore,
) -> Vec<ChainAln> {
    let anchors = plausible_anchors(chain, read_seq, rseq, kmer_len);

    if anchors.is_empty() {
        return Vec::new();
    }

    // Core alignments of anchor runs, split at implausible gaps.
    let segs = align_anchor_runs(&anchors, read_seq, rseq, kmer_len);

    (0..segs.len())
        .flat_map(|seg_idx| finish_segment(chain, seg_idx, &segs, &anchors, read_seq, rseq, kmer_len, scoring))
        .collect()
}

/// Ranking cost of a chain's alignment, summed over its segments.
/// Each segment's flank cost already charges the bases it clips.
pub(super) fn chain_cost(segs: &[ChainAln]) -> Option<i64> {
    if segs.is_empty() {
        None
    } else {
        Some(segs.iter().map(|aln| aln.cost).sum())
    }
}

/// Chaining thresholds.
/// Long reads use the defaults.
/// The short-read rescue of reads over -e accepts sparser chains, because few error-free minimizers survive in a high-error 150 bp read, which is also why minimap2's sr preset uses -n2 -m20.
/// It then accepts only alignments whose records together cover at least RESCUE_MIN_ALIGNED of the read with at most RESCUE_MAX_NM_RATE edits per aligned base.
#[derive(Clone, Copy)]
pub(super) struct ChainParams {
    pub(super) min_score: f64,
    pub(super) min_anchors: usize,
    pub(super) rescue: bool,
}

pub(super) const LONG_CHAIN: ChainParams = ChainParams {
    min_score: crate::chain::MIN_CHAIN_SCORE,
    min_anchors: crate::chain::MIN_CHAIN_ANCHORS,
    rescue: false,
};

pub(super) const RESCUE_CHAIN: ChainParams = ChainParams {
    min_score: 25.0,
    min_anchors: 2,
    rescue: true,
};

/// Rescue acceptance thresholds, the aligned read fraction and the edit density.
pub(super) const RESCUE_MIN_ALIGNED: f64 = 0.8;

pub(super) const RESCUE_MAX_NM_RATE: f64 = AUTO_ED_RATE;

/// Split rescue limits, each part's minimum length, the allowed query overlap between parts, and the most parts emitted.
pub(super) const RESCUE_MIN_PART: usize = 40;

pub(super) const RESCUE_PART_OVERLAP: usize = 10;

pub(super) const RESCUE_MAX_PARTS: usize = 3;

/// MAPQ ceiling of a rescued short read, set to 3.
/// The chain MAPQ is fitted to long reads with hundreds of anchors.
/// A rescued read is one the short-read aligner could not place within -e, its few anchors leave most repeat copies unchained, and such reads were wrong 54% of the time at chain MAPQ 11-35 on simulated reads at 0.5-8% error from hg38 and repeat-rich genomes.
pub(super) const RESCUE_MAPQ_MAX: u8 = 3;

/// Chains within this share of the best chain's score over the same read span are near-tied, and up to MAX_TIED_CHAINS of them are aligned and ranked by alignment cost.
pub(super) const NEAR_TIE_SHARE: f64 = 0.9;

pub(super) const MAX_TIED_CHAINS: usize = 5;

/// Records kept per read, primary plus supplementaries, and the share of the best chain's score a supplementary chain of a read needs.
pub(super) const MAX_READ_RECORDS: usize = 6;

pub(super) const SUPPLEMENTARY_MIN_SHARE: f64 = 0.3;

/// A genome-scale query keeps up to one record per this many query bases.
pub(super) const GENOME_BASES_PER_RECORD: usize = 10_000;

/// OPTIK_MAPQ_DEBUG prints at most this many chains per read.
pub(super) const DEBUG_MAX_CHAINS: usize = 8;

/// Spare capacity reserved beyond the read length when a CIGAR is expanded into columns or ops.
pub(super) const CIGAR_COLS_SLACK: usize = 16;

/// left_align stops shifting a gap once the shift adds more than this many mismatches.
pub(super) const SHIFT_MAX_EXTRA_MISMATCHES: i64 = 4;

/// left_align stops shifting a gap without a gain after this many columns.
pub(super) const SHIFT_MAX_NEUTRAL_COLS: usize = 256;

/// A chain's query span in forward read coordinates, which overlap tests need whatever the chain's strand.
fn forward_span(chain: &crate::chain::Chain, read_len: u32, kmer_len: u32) -> (u32, u32) {
    let (span_start, span_end) = chain.q_span();
    let (start, end) = (span_start, span_end + kmer_len);

    if chain.rev {
        (read_len - end, read_len - start)
    } else {
        (start, end)
    }
}

/// Two query spans overlap when they share more than half of the shorter one.
fn spans_overlap(span_a: (u32, u32), span_b: (u32, u32)) -> bool {
    // The spans overlap when the shared part exceeds half of the shorter one.
    const SHORTER_SPAN_SHARE_DEN: u32 = 2;
    let overlap_len = span_a
        .1
        .min(span_b.1)
        .saturating_sub(span_a.0.max(span_b.0));

    overlap_len * SHORTER_SPAN_SHARE_DEN > (span_a.1 - span_a.0).min(span_b.1 - span_b.0)
}

/// Aligns up to MAX_TIED_CHAINS near-tied chains over the primary's read span and moves the cheapest to the front.
/// Equal costs are broken pseudo-randomly by the read, and the result is the (chosen, runner-up) alignment scores used as MAPQ evidence.
fn resolve_near_ties(
    chains: &mut [crate::chain::Chain],
    read_name: &str,
    seq: &[u8],
    rev_seq: &[u8],
    rep_len: u32,
    ctx: &AlignCtx<'_>,
    mm_sc: MmScore,
) -> Option<(f64, f64)> {
    let (len, kmer_len) = (seq.len() as u32, ctx.kmer_len as u32);
    let forward_span_of = |chain: &crate::chain::Chain| forward_span(chain, len, kmer_len);
    let overlap = spans_overlap;
    let mut tie_dp: Option<(f64, f64)> = None;
    let top = chains[0].score;
    let span0 = forward_span_of(&chains[0]);
    // Near ties are chains within NEAR_TIE_SHARE of the top score that cover the same part of the read, so they compete for one placement.
    let cand: Vec<usize> = (0..chains.len())
        .filter(|&chain_idx| {
            chains[chain_idx].score >= NEAR_TIE_SHARE * top
                && overlap(forward_span_of(&chains[chain_idx]), span0)
        })
        .take(MAX_TIED_CHAINS)
        .collect();

    // OPTIK_DEBUG_CHAIN realigns the top chains just to print their cost.
    if env_flag("OPTIK_DEBUG_CHAIN", &DEBUG_CHAIN) {
        eprintln!(
            "[chain] {read_name} len={} rep_len={rep_len} chains={} near_tied={}",
            seq.len(),
            chains.len(),
            cand.len()
        );
        for chain in chains.iter().take(DEBUG_MAX_CHAINS) {
            let oriented: &[u8] = if chain.rev { rev_seq } else { seq };
            let cost = chain_cost(&align_chain(
                chain,
                oriented,
                &ctx.ref_seqs[chain.rid as usize],
                ctx.kmer_len,
                mm_sc,
            ));

            eprintln!(
                "[chain] {read_name} rid={} rev={} score={:.1} anchors={} r={}..{} q={:?} cost={:?}",
                chain.rid,
                chain.rev,
                chain.score,
                chain.anchors.len(),
                chain.anchors[0].ref_pos,
                chain.anchors[chain.anchors.len() - 1].ref_pos,
                forward_span_of(chain),
                cost
            );
        }
    }

    // Chain scores only count anchors, so the near ties are aligned and ranked by their real alignment cost.
    if cand.len() > 1 {
        let costs: Vec<(i64, usize)> = cand
            .iter()
            .filter_map(|&chain_idx| {
                let chain = &chains[chain_idx];
                let oriented: &[u8] = if chain.rev { rev_seq } else { seq };

                chain_cost(&align_chain(
                    chain,
                    oriented,
                    &ctx.ref_seqs[chain.rid as usize],
                    ctx.kmer_len,
                    mm_sc,
                ))
                .map(|cost| (cost, chain_idx))
            })
            .collect();

        if let Some(&(min_cost, _)) = costs.iter().min() {
            // Equal-cost copies (a read inside an identical repeat such as an rRNA operon) are chosen between pseudo-randomly, seeded by the read, as on the short-read path.
            // Always taking the first left the other copies without long-read coverage.
            let mut tied: Vec<usize> = costs
                .iter()
                .filter(|&&(cost, _)| cost == min_cost)
                .map(|&(_, chain_idx)| chain_idx)
                .collect();

            tied.sort_by_key(|&chain_idx| {
                (
                    chains[chain_idx].rid,
                    chains[chain_idx].anchors[0].ref_pos,
                    chains[chain_idx].rev,
                )
            });
            let chosen = tied[(read_hash(seq) % tied.len() as u64) as usize];
            // MAPQ evidence is the alignment score of the chosen copy and of the best other one, on minimap2's scale (match +2), from the same cost that ranked them, as 2 x read length - 1.5 x cost.
            // A mismatch costs 4 here and 6 score units there, and clipped bases are charged in the cost.
            let dp_score = |cost: i64| {
                MM2_MATCH_SCORE * seq.len() as f64 - MM2_SCORE_PER_COST * cost as f64
            };
            let runner_up = costs
                .iter()
                .filter(|entry| entry.1 != chosen)
                .map(|entry| entry.0)
                .min();

            if let Some(runner_cost) = runner_up {
                tie_dp = Some((dp_score(min_cost), dp_score(runner_cost)));
            }

            // The near tie with the lowest alignment cost is swapped to the front and becomes the primary.
            chains.swap(0, chosen);
        }
    }

    tie_dp
}

/// Picks the chains that become records, always starting with the primary.
fn select_kept_chains(chains: &[crate::chain::Chain], read_len: usize, genome_scale: bool, kmer_len: u32) -> Vec<usize> {
    let forward_span_of = |chain: &crate::chain::Chain| forward_span(chain, read_len as u32, kmer_len);
    // A genome-scale query can span many loci. Its record cap grows with read_len / GENOME_BASES_PER_RECORD.
    let max_kept = if genome_scale {
        (read_len / GENOME_BASES_PER_RECORD).max(MAX_READ_RECORDS)
    } else {
        MAX_READ_RECORDS
    };
    let mut kept: Vec<usize> = vec![0];

    for (chain_idx, chain) in chains.iter().enumerate().skip(1) {
        if kept.len() >= max_kept {
            break;
        }

        if !genome_scale && chain.score < SUPPLEMENTARY_MIN_SHARE * chains[0].score {
            break;
        }

        if kept
            .iter()
            .all(|&kept_idx| !spans_overlap(forward_span_of(chain), forward_span_of(&chains[kept_idx])))
        {
            kept.push(chain_idx);
        }
    }

    kept
}

/// Chaining-based mapping, for reads of at least CHAIN_MIN_LEN and for short reads over their edit budget.
/// Seeds become anchor chains per (reference, strand), then a primary and query-disjoint supplementary records.
pub(super) fn map_one_read_chain(
    read_name: &str,
    seq: &[u8],
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
    prm: ChainParams,
) -> Vec<String> {
    use crate::chain::{Anchor, Chain, MapqEvidence, chain_group, mapq, mapq_long};
    let seq_s = if seq.is_empty() {
        "*"
    } else {
        std::str::from_utf8(seq).unwrap_or("*")
    };
    let unmapped = || {
        vec![format!(
            "{read_name}\t{SAM_FLAG_UNMAPPED}\t*\t0\t0\t*\t*\t0\t0\t{seq_s}\t*"
        )]
    };
    // Long reads seed with repeat-run rescue through collect_chain and keep the repeat-covered length for MAPQ.
    // The short-read rescue keeps the seeds its fitted MAPQ cap was measured with.
    let (hits, rep_len) = if prm.rescue {
        (ctx.seeds.collect(seq), 0)
    } else {
        ctx.seeds.collect_chain(seq)
    };
    let build_chains = |mut hits: Vec<crate::SeedHit>| -> Vec<Chain> {
        let mut chains: Vec<Chain> = Vec::new();

        if hits.len() < prm.min_anchors {
            return chains;
        }

        // The sort groups hits by (rid, rev), and chain_group() runs once per group.
        hits.sort_unstable_by_key(|hit| (hit.rid, hit.rev, hit.rpos, hit.qpos));
        hits.dedup_by_key(|hit| (hit.rid, hit.rev, hit.rpos, hit.qpos));
        let mut i = 0;
        let mut grp: Vec<Anchor> = Vec::new();

        while i < hits.len() {
            let (rid, rev) = (hits[i].rid, hits[i].rev);

            grp.clear();

            while i < hits.len() && hits[i].rid == rid && hits[i].rev == rev {
                grp.push(Anchor {
                    ref_pos: hits[i].rpos,
                    read_pos: hits[i].qpos,
                });
                i += 1;
            }

            chain_group(
                rid,
                rev != 0,
                &grp,
                ctx.kmer_len,
                prm.min_score,
                prm.min_anchors,
                &mut chains,
            );
        }

        chains
    };
    let mut chains = build_chains(hits);

    // A read whose repeat runs left too few seeds for any chain is seeded again at the index's full occurrence cap, as minimap2 re-chains with a higher threshold.
    // rep_len keeps the repeat evidence of the first pass.
    if chains.is_empty() && !prm.rescue && rep_len > 0 {
        chains = build_chains(ctx.seeds.collect(seq));
    }

    if chains.is_empty() {
        return unmapped();
    }

    chains.sort_by(|lhs, rhs| {
        rhs.score
            .total_cmp(&lhs.score)
            .then(lhs.rid.cmp(&rhs.rid))
            .then(lhs.anchors[0].ref_pos.cmp(&rhs.anchors[0].ref_pos))
    });
    // Query spans in FORWARD read coordinates, for overlap tests.
    let len = seq.len() as u32;
    let kmer_len = ctx.kmer_len as u32;
    let fwd_span = |chain: &Chain| forward_span(chain, len, kmer_len);
    let overlap = spans_overlap;
    // Near-tied chains over the same read span are repeat copies that chain scores cannot separate, but base-level alignment can, so up to 5 are aligned and the lowest gap-affine cost becomes the primary.
    // Genome-scale queries, i.e. assemblies, are scored as minimap2 -x asm20 scores them, and reads with its defaults, the local clip and asm20's minimal score ASM_MIN_DP.
    let genome_scale = seq.len() >= GENOME_SCALE_LEN;
    let mm_sc = if genome_scale { MM_ASM20 } else { MM_READS };
    let rev_seq = revcomp(seq);
    // Best and second-best base-level cost over the primary's read span, kept as MAPQ evidence when near-tied copies were aligned.
    let tie_dp = resolve_near_ties(&mut chains, read_name, seq, &rev_seq, rep_len, ctx, mm_sc);
    // The primary and query-disjoint supplementaries are emitted.
    // Overlapping chains are secondary, not emitted, and feed MAPQ.
    // A read has one locus and a few chimeric parts, so query-disjoint chains within 30% of the best score are kept, at most 6.
    // A genome-scale query, such as an assembly in BWA-MEM's MG1655-vs-536 test, is a series of collinear blocks whose sizes bear no relation to the largest one.
    // Every query-disjoint chain that passed chaining's own thresholds is therefore kept, up to one per 10 kb of query, the way minimap2 keeps all query-disjoint primaries.
    let kept = select_kept_chains(&chains, seq.len(), genome_scale, ctx.kmer_len as u32);

    let sum_kept: f64 = kept.iter().map(|&kept_idx| chains[kept_idx].score).sum();
    let mut alns: Vec<(ChainAln, u8)> = Vec::new();

    for &chain_idx in &kept {
        let chain = &chains[chain_idx];
        let competitor_score = chains
            .iter()
            .enumerate()
            .filter(|(other_idx, other)| {
                *other_idx != chain_idx && overlap(fwd_span(other), fwd_span(chain))
            })
            .map(|(_, other)| other.score)
            .fold(0.0f64, f64::max);
        // competitor_score, the best overlapping rival, decides how unique this chain looks to mapq().
        let oriented: &[u8] = if chain.rev { &rev_seq } else { seq };
        let mut segs = align_chain(
            chain,
            oriented,
            &ctx.ref_seqs[chain.rid as usize],
            ctx.kmer_len,
            mm_sc,
        );

        if genome_scale {
            // Diverged homologous blocks go unreported, per minimap2 -x asm20's min_dp_max.
            segs.retain(|seg| {
                asm_score(&seg.cigar, oriented, &ctx.ref_seqs[seg.rid], seg.pos) >= ASM_MIN_DP
            });
            if segs.is_empty() {
                continue;
            }
        }

        let n_anchors: usize = segs.iter().map(|seg| seg.n_anchors).sum();
        let mapq_val = if prm.rescue {
            // Rescued reads use the chain MAPQ capped at RESCUE_MAPQ_MAX, the range its calibration covered.
            mapq(chain.score, competitor_score, n_anchors).min(RESCUE_MAPQ_MAX)
        } else {
            let (cols, edits) = segs
                .iter()
                .fold((0usize, 0i64), |(col_sum, edit_sum), seg| {
                    (
                        col_sum + cigar_cols(&seg.cigar),
                        edit_sum + seg.edit_dist as i64,
                    )
                });
            let identity = if cols > 0 {
                ((cols as i64 - edits).max(0) as f64) / cols as f64
            } else {
                0.0
            };

            mapq_long(&MapqEvidence {
                chain_score: chain.score,
                competitor_score,
                n_anchors,
                sum_kept,
                rep_len,
                identity,
                align_scores: if chain_idx == 0 { tie_dp } else { None },
            })
        };

        // The primary is the longest segment of the best chain, and a split chain's other segments become supplementaries.
        if alns.is_empty()
            && let Some(longest) = (0..segs.len()).max_by_key(|&seg_idx| {
                (
                    segs[seg_idx].read_end - segs[seg_idx].read_start,
                    std::cmp::Reverse(seg_idx),
                )
            })
        {
            segs.swap(0, longest);
        }

        for seg in segs {
            alns.push((seg, mapq_val));
        }
    }

    if alns.is_empty() {
        return unmapped();
    }

    if prm.rescue {
        // Either one near-end-to-end, low-divergence record or, for a read across a contig end or a breakpoint, a primary plus supplementaries whose query-disjoint parts are each long and low-divergence and together cover the read.
        let good = |aln: &ChainAln| {
            let aligned = (aln.read_end - aln.read_start) as f64;

            aln.edit_dist as f64 <= RESCUE_MAX_NM_RATE * aligned
        };
        let need = RESCUE_MIN_ALIGNED * seq.len() as f64;
        let primary = &alns[0].0;

        if (primary.read_end - primary.read_start) as f64 >= need && good(primary) {
            alns.truncate(1);
        } else {
            let read_len = seq.len();
            let fwd = |aln: &ChainAln| {
                if aln.rev {
                    (read_len - aln.read_end, read_len - aln.read_start)
                } else {
                    (aln.read_start, aln.read_end)
                }
            };
            let mut parts: Vec<(ChainAln, u8)> = Vec::new();

            for (aln, mapq_val) in alns.drain(..) {
                let (start, end) = fwd(&aln);

                if end - start < RESCUE_MIN_PART || !good(&aln) {
                    continue;
                }

                // Strictly query-disjoint, allowing a few bases of shared overlap.
                if parts.iter().any(|(part, _)| {
                    let (part_start, part_end) = fwd(part);

                    end.min(part_end).saturating_sub(start.max(part_start)) > RESCUE_PART_OVERLAP
                }) {
                    continue;
                }

                parts.push((aln, mapq_val));

                if parts.len() == RESCUE_MAX_PARTS {
                    break;
                }
            }

            // The parts' read intervals are merged in order, so a base covered twice counts once.
            let mut intervals: Vec<(usize, usize)> =
                parts.iter().map(|(aln, _)| fwd(aln)).collect();

            intervals.sort_unstable();
            let (mut covered, mut end) = (0usize, 0usize);

            for (iv_start, iv_end) in intervals {
                let clipped_start = iv_start.max(end);

                if iv_end > clipped_start {
                    covered += iv_end - clipped_start;
                    end = iv_end;
                }
            }

            // A split read needs at least two parts that together cover RESCUE_MIN_ALIGNED of it, and otherwise it stays unmapped.
            const MIN_SPLIT_PARTS: usize = 2;

            if parts.len() < MIN_SPLIT_PARTS || (covered as f64) < need {
                return unmapped();
            }

            alns = parts;
        }
    }

    // Every record lists the other parts in its SA tag, giving each one's rname, pos, strand, CIGAR, mapQ and NM as the SAM spec defines.
    let sa_entry = |aln: &ChainAln, mapq_val: u8| {
        format!(
            "{},{},{},{},{},{};",
            ref_names[aln.rid],
            aln.pos + 1,
            if aln.rev { '-' } else { '+' },
            aln.cigar,
            mapq_val,
            aln.edit_dist
        )
    };

    alns.iter()
        .enumerate()
        .map(|(i, (aln, mapq_val))| {
            // Records after the first get SAM_FLAG_SUPPLEMENTARY.
            let flag = (if aln.rev { SAM_FLAG_REVERSE } else { 0 })
                | (if i > 0 { SAM_FLAG_SUPPLEMENTARY } else { 0 });
            let sseq = if aln.rev {
                std::str::from_utf8(&rev_seq).unwrap_or("*")
            } else {
                seq_s
            };
            let sbytes: &[u8] = if aln.rev { &rev_seq } else { seq };
            let tags = diff_tag_fields(
                &aln.cigar,
                sbytes,
                &ctx.ref_seqs[aln.rid],
                aln.pos,
                ctx.cs_mode,
            );
            // Supplementary records hard-clip and keep only their aligned bases (the default of minimap2 and bwa), since a soft-clipped copy of the whole query per record makes a genome-scale query's SAM grow with records x query length.
            let (cig, sseq) = if i > 0 {
                hard_clip(&aln.cigar, sseq)
            } else {
                (aln.cigar.clone(), sseq)
            };
            let mut line = format!(
                "{read_name}\t{flag}\t{}\t{}\t{mapq_val}\t{cig}\t*\t0\t0\t{sseq}\t*\tNM:i:{}{tags}",
                ref_names[aln.rid],
                aln.pos + 1,
                aln.edit_dist
            );
            let sa_tag: String = alns
                .iter()
                .enumerate()
                .filter(|(other_idx, _)| *other_idx != i)
                .map(|(_, (other, omq))| sa_entry(other, *omq))
                .collect();

            if !sa_tag.is_empty() {
                line.push_str(&format!("\tSA:Z:{sa_tag}"));
            }

            line
        })
        .collect()
}
