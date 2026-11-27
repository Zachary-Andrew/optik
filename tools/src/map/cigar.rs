//! Builds and edits CIGARs.
//! The unit-cost DP produces them, and the helpers here count columns, left-align indels, clip low-scoring ends and convert soft clips to hard ones.

use super::*;

/// Traceback direction codes for the CIGAR DP matrix.
/// `back[i][j]` stores which DP recurrence produced the minimum score at cell (i, j).
pub(super) const DIR_MATCH: u8 = 0; // diagonal step, query[i] against text[j] as an M op, match or mismatch

pub(super) const DIR_DEL: u8 = 1; // vertical step, text[j] without a query counterpart, a D op

pub(super) const DIR_INS: u8 = 2; // horizontal step, query[i] without a text counterpart, an I op

/// Appends the ASCII digit `byte` to a run length being parsed from a CIGAR string.
pub(super) fn push_digit<T>(value: T, byte: u8) -> T
where
    T: From<u8> + std::ops::Mul<Output = T> + std::ops::Add<Output = T>,
{
    const DECIMAL_RADIX: u8 = 10;

    value * T::from(DECIMAL_RADIX) + T::from(byte - b'0')
}

/// Semi-global unit-cost alignment CIGAR by dynamic programming, with the query end to end and the text's ends free.
/// Returns (CIGAR, 0-based start offset in the text).
/// This is the fallback when the wavefront aligner gives up, and the reference path under OPTIK_FORCE_SCALAR_CIGAR.
pub(super) fn cigar(query: &[u8], text: &[u8]) -> (String, usize) {
    let query_len = query.len();
    let text_len = text.len();

    if query_len == 0 {
        return ("*".into(), 0);
    }

    // This is semi-global edit distance, with the query aligned end to end and free text on both sides.
    // Row 0 starts every query prefix as insertions, while column 0 costs nothing, so the query can begin anywhere in the text.
    let mut back = vec![vec![0u8; query_len + 1]; text_len + 1];

    back[0].fill(DIR_INS);

    let mut prev = (0..=query_len as i32).collect::<Vec<i32>>();

    let mut best_i = 0usize;
    let mut best_score = query_len as i32;

    for i in 1..=text_len {
        let mut cur = vec![0i32; query_len + 1];

        cur[0] = 0;
        back[i][0] = DIR_DEL;

        for j in 1..=query_len {
            let cost = if text[i - 1] == query[j - 1] { 0i32 } else { 1 };
            let diag = prev[j - 1] + cost;
            let from_up = prev[j] + 1;
            let left = cur[j - 1] + 1;

            // Ties go to the gaps, and between them to deletion, so the traceback is the same on every run.
            if diag < from_up && diag < left {
                cur[j] = diag;
                back[i][j] = DIR_MATCH;
            } else if from_up <= left {
                cur[j] = from_up;
                back[i][j] = DIR_DEL;
            } else {
                cur[j] = left;
                back[i][j] = DIR_INS;
            }
        }

        // The last row with the lowest cost is kept, which ends the alignment as far right as the score allows.
        if cur[query_len] <= best_score {
            best_score = cur[query_len];
            best_i = i;
        }

        prev = cur;
    }

    let mut ops: Vec<u8> = Vec::new();
    // The traceback walks back from the chosen end. Ops come out reversed.
    let (mut text_pos, mut query_pos) = (best_i as i32, query_len as i32);

    while query_pos > 0 {
        match back[text_pos as usize][query_pos as usize] {
            DIR_MATCH => {
                ops.push(b'M');
                text_pos -= 1;
                query_pos -= 1;
            }
            DIR_DEL => {
                ops.push(b'D');
                text_pos -= 1;
            }
            DIR_INS => {
                ops.push(b'I');
                query_pos -= 1;
            }
            _ => unreachable!(),
        }
    }

    let text_start = text_pos as usize;

    ops.reverse();

    // Deletions at either end only move the start in the text, so they are trimmed rather than written into the CIGAR.
    let lead = ops.iter().take_while(|&&cigar_op| cigar_op == b'D').count();
    let tail = ops
        .iter()
        .rev()
        .take_while(|&&cigar_op| cigar_op == b'D')
        .count();
    let core = &ops[lead..ops.len() - tail];

    let mut cigar = String::new();

    for run in core.chunk_by(|lhs, rhs| lhs == rhs) {
        cigar.push_str(&run.len().to_string());
        cigar.push(run[0] as char);
    }

    // An alignment made only of trimmed deletions still needs a valid CIGAR for its query.
    if cigar.is_empty() {
        cigar = format!("{}M", query_len);
    }

    (cigar, text_start)
}

/// Recomputes the edit distance a CIGAR implies against (query, text), where a mismatch inside M costs 1 and so does every I and D base.
/// Returns None if the CIGAR overshoots either string, which means inconsistent bookkeeping.
pub(super) fn cigar_ed_implied(cigar: &str, query: &[u8], text: &[u8], text_start: usize) -> Option<i32> {
    let (mut query_pos, mut text_pos) = (0usize, text_start);
    let mut edits = 0i32;
    let mut num: u64 = 0;

    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            num = push_digit(num, byte);
            continue;
        }

        if num == 0 {
            return None;
        }

        let run_len = num as usize;

        num = 0;

        match byte {
            b'M' => {
                for _ in 0..run_len {
                    if query_pos >= query.len() || text_pos >= text.len() {
                        return None;
                    }

                    // A non-ACGT base (N) never matches, as in SAM's NM.
                    if query[query_pos] != text[text_pos]
                        || !matches!(text[text_pos], b'A' | b'C' | b'G' | b'T')
                    {
                        edits += 1;
                    }

                    query_pos += 1;
                    text_pos += 1;
                }
            }
            b'I' => {
                for _ in 0..run_len {
                    if query_pos >= query.len() {
                        return None;
                    }

                    edits += 1;
                    query_pos += 1;
                }
            }
            b'D' => {
                for _ in 0..run_len {
                    if text_pos >= text.len() {
                        return None;
                    }

                    edits += 1;
                    text_pos += 1;
                }
            }
            _ => return None,
        }
    }

    Some(edits)
}

/// Alignment columns of a CIGAR (M, =, X, I, D).
pub(super) fn cigar_cols(cigar: &str) -> usize {
    let (mut run_len, mut tot) = (0usize, 0usize);

    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        if matches!(byte, b'M' | b'=' | b'X' | b'I' | b'D') {
            tot += run_len;
        }

        run_len = 0;
    }

    tot
}

/// Whether a SAM record has the unmapped flag set or no readable flag.
pub(super) fn is_unmapped(rec: &str) -> bool {
    rec.split('\t')
        .nth(1)
        .and_then(|flag_field| flag_field.parse::<u16>().ok())
        .is_none_or(|flag| flag & SAM_FLAG_UNMAPPED != 0)
}

/// Appends `count` copies of op to a run-length op list, merging with the tail.
pub(super) fn push_op(ops: &mut Vec<(u32, u8)>, count: u32, cigar_op: u8) {
    if count == 0 {
        return;
    }

    if let Some(last) = ops.last_mut()
        && last.1 == cigar_op
    {
        last.0 += count;
        return;
    }

    ops.push((count, cigar_op));
}

/// Appends every operation of a CIGAR string to a run-length op list.
pub(super) fn push_cigar(ops: &mut Vec<(u32, u8)>, cig: &str) {
    let mut run_len = 0u32;

    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
        } else {
            push_op(ops, run_len, byte);
            run_len = 0;
        }
    }
}

/// The CIGAR's leading and trailing hard-clip lengths, 0 where absent.
pub(super) fn hard_clips(cig: &str) -> (usize, usize) {
    let bytes = cig.as_bytes();
    let lead_digits = bytes
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    let lead = if lead_digits > 0 && bytes.get(lead_digits) == Some(&b'H') {
        cig[..lead_digits].parse().unwrap_or(0)
    } else {
        0
    };
    let mut trail = 0;

    // A lone clip such as 5H is the leading one. It must not count twice as trailing.
    if bytes.len() > 1 && bytes[bytes.len() - 1] == b'H' {
        let mut digit_start = bytes.len() - 1;

        while digit_start > 0 && bytes[digit_start - 1].is_ascii_digit() {
            digit_start -= 1;
        }

        if digit_start < bytes.len() - 1 && !(digit_start == 0 && lead > 0) {
            trail = cig[digit_start..bytes.len() - 1].parse().unwrap_or(0);
        }
    }

    (lead, trail)
}

/// `cig` with its end soft clips turned into hard clips, and `seq` cut to the bases the CIGAR still consumes.
pub(super) fn hard_clip<'a>(cig: &str, seq: &'a str) -> (String, &'a str) {
    let bytes = cig.as_bytes();
    let digits = |from: usize| {
        bytes[from..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count()
    };
    let mut lead = 0usize;
    let mut body = 0..bytes.len();
    let lead_digits = digits(0);

    if lead_digits > 0 && bytes.get(lead_digits) == Some(&b'S') {
        lead = cig[..lead_digits].parse().unwrap_or(0);
        body.start = lead_digits + 1;
    }

    let mut trail = 0usize;

    if bytes.last() == Some(&b'S') && body.end > body.start + 1 {
        let mut digit_start = body.end - 1;

        while digit_start > body.start && bytes[digit_start - 1].is_ascii_digit() {
            digit_start -= 1;
        }

        if digit_start < body.end - 1 {
            trail = cig[digit_start..body.end - 1].parse().unwrap_or(0);
            body.end = digit_start;
        }
    }

    // Nothing to clip, or clips longer than the read. Either way the input goes back unchanged.
    if lead == 0 && trail == 0 || lead + trail > seq.len() {
        return (cig.to_string(), seq);
    }

    let mut out = String::with_capacity(cig.len());

    if lead > 0 {
        out.push_str(&format!("{lead}H"));
    }

    out.push_str(&cig[body]);

    if trail > 0 {
        out.push_str(&format!("{trail}H"));
    }

    (out, &seq[lead..seq.len() - trail])
}

/// Left-aligns every indel of an M/I/D CIGAR of `query` against `text` from `toff`.
/// Each gap moves to the leftmost position within the match run before it at which the columns it passes have no more mismatches than before.
/// The gap's cost does not depend on its position, so the score stays the same or gets better.
/// This covers the normalization of GATK LeftAlignIndels and bcftools norm, shifting through a repeat, and also equal-score placements that move a mismatch.
/// Equivalent placements of an indel thus collapse onto the leftmost, where pileup callers, GIAB and BWA-MEM put it.
/// A mismatch next to the gap moves with it, so a read with a 1-base deletion and an A>T SNP in an A run is written as deletion then SNP, as GIAB and BWA-MEM write it.
/// Writing a mismatch beside the deletion at the SNP's place hid the SNP and cost 148 calls on HG002 chr20.
/// An unparseable CIGAR is returned unchanged.
pub(super) fn left_align(cig: &str, query: &[u8], text: &[u8], toff: usize) -> String {
    let mut ops: Vec<u8> = Vec::with_capacity(query.len() + CIGAR_COLS_SLACK);
    let mut run_len = 0usize;

    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        if !matches!(byte, b'M' | b'I' | b'D') {
            return cig.to_string();
        }

        // One op per column makes rotating a gap a slice operation.
        ops.extend(std::iter::repeat_n(byte, run_len));
        run_len = 0;
    }

    if !ops.iter().any(|&col_op| col_op != b'M') {
        return cig.to_string();
    }

    // Query and text position before each column, in one pass.
    // Shifting a gap only rearranges columns before its end, so the positions of every later column stay valid and the gaps are normalized left to right.
    let mut col_pos: Vec<(usize, usize)> = Vec::with_capacity(ops.len());

    {
        let (mut query_i, mut text_i) = (0usize, toff);

        for &col_op in &ops {
            col_pos.push((query_i, text_i));

            if col_op != b'D' {
                query_i += 1;
            }

            if col_op != b'I' {
                text_i += 1;
            }
        }

        if query_i > query.len() || text_i > text.len() {
            return cig.to_string();
        }
    }
    let base_ok = |base: u8| matches!(base, b'A' | b'C' | b'G' | b'T');
    let mut col = 0usize;

    while col < ops.len() {
        if ops[col] == b'M' {
            col += 1;
            continue;
        }

        let gap_op = ops[col];
        let mut gap_end = col;

        while gap_end < ops.len() && ops[gap_end] == gap_op {
            gap_end += 1;
        }

        let len = gap_end - col;
        // The run [lo, lo + len) starts at query qs and text ts.
        // Shifting it left by j re-pairs the j columns before it, so a deletion pairs q[qs - i] with t[ts + len - i] instead of t[ts - i], and an insertion pairs t[ts - i] with q[qs + len - i] instead of q[qs - i].
        // The gap cost is unchanged, so the leftmost shift whose re-paired columns have no more mismatches is an alignment of equal score.
        let (gap_col, (gap_query, gap_text)) = (col, col_pos[col]);
        let mismatch = |query_base: u8, text_base: u8| {
            query_base != text_base || !base_ok(query_base) || !base_ok(text_base)
        };
        let mut run = 0usize;

        // The gap may only move through the match run directly before it. Another gap stops it.
        while run < gap_col && ops[gap_col - 1 - run] == b'M' {
            run += 1;
        }

        let run = run.min(gap_query).min(gap_text - toff);
        let (mut delta, mut best) = (0i64, 0usize);

        for i in 1..=run {
            let (old, new) = if gap_op == b'D' {
                (
                    mismatch(query[gap_query - i], text[gap_text - i]),
                    mismatch(query[gap_query - i], text[gap_text + len - i]),
                )
            } else {
                (
                    mismatch(query[gap_query - i], text[gap_text - i]),
                    mismatch(query[gap_query + len - i], text[gap_text - i]),
                )
            };

            // delta tracks the mismatch change of the shift so far. Ties move left too.
            delta += new as i64 - old as i64;

            if delta <= 0 {
                best = i;
            }

            if delta > SHIFT_MAX_EXTRA_MISMATCHES || (i > SHIFT_MAX_NEUTRAL_COLS && delta > 0) {
                break;
            }
        }

        if best > 0 {
            // Rotate the `best` columns from before the gap to after it.
            ops[gap_col - best..gap_col + len].rotate_left(best);
        }

        col = gap_end;
    }

    let mut out = String::new();
    let mut run_start = 0;

    while run_start < ops.len() {
        let mut run_end = run_start;

        while run_end < ops.len() && ops[run_end] == ops[run_start] {
            run_end += 1;
        }

        out.push_str(&format!(
            "{}{}",
            run_end - run_start,
            ops[run_start] as char
        ));
        run_start = run_end;
    }

    out
}

/// Soft-clips the ends of an end-to-end alignment, an M/I/D CIGAR of `query` against `text` from `toff`, under the EXT_* scoring.
/// Each end is cut at its lowest-scoring prefix or suffix if that score is below zero, and only just before a matching column, so the kept core starts and ends with a match.
/// Returns the CIGAR with S ops, the new reference offset and the core's NM, or None when nothing is clipped or the core would keep under half the read.
pub(super) fn clip_ends(cig: &str, query: &[u8], text: &[u8], toff: usize) -> Option<(String, usize, i32)> {
    // Expand into columns of op and is_match.
    let mut cols: Vec<(u8, bool)> = Vec::with_capacity(query.len() + CIGAR_COLS_SLACK);
    let (mut query_i, mut text_i, mut run_len) = (0usize, toff, 0usize);

    // A CIGAR that runs past either sequence cannot be clipped, so it is rejected as it is expanded.
    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        for _ in 0..run_len {
            match byte {
                b'M' => {
                    if query_i >= query.len() || text_i >= text.len() {
                        return None;
                    }

                    let is_match = query[query_i] == text[text_i]
                        && matches!(text[text_i], b'A' | b'C' | b'G' | b'T');

                    cols.push((b'M', is_match));
                    query_i += 1;
                    text_i += 1;
                }
                b'I' => {
                    if query_i >= query.len() {
                        return None;
                    }

                    cols.push((b'I', false));
                    query_i += 1;
                }
                b'D' => {
                    if text_i >= text.len() {
                        return None;
                    }

                    cols.push((b'D', false));
                    text_i += 1;
                }
                _ => return None,
            }
        }

        run_len = 0;
    }

    // A CIGAR that stops short of the read's end is malformed.
    if query_i != query.len() {
        return None;
    }

    // A gap column pays the open cost only when it starts a new gap.
    let score = |col: (u8, bool), prev: Option<u8>| -> i64 {
        match col {
            (b'M', true) => EXT_MATCH,
            (b'M', false) => EXT_MISMATCH,
            (gap_op, _) => {
                if prev == Some(gap_op) {
                    EXT_GAP_EXT
                } else {
                    EXT_GAP_OPEN + EXT_GAP_EXT
                }
            }
        }
    };
    let num_cols = cols.len();
    // From the left, the cut goes after the lowest running score that a match follows, so the kept core starts on a match.
    let (mut cum, mut best, mut lcut) = (0i64, 0i64, 0usize);

    for c in 0..num_cols {
        cum += score(cols[c], if c > 0 { Some(cols[c - 1].0) } else { None });

        if cum < best && c + 1 < num_cols && cols[c + 1] == (b'M', true) {
            best = cum;
            lcut = c + 1;
        }
    }

    // The right end is scanned the same way, from the last column backwards.
    let (mut cum, mut best, mut rcut) = (0i64, 0i64, 0usize);

    for c in (0..num_cols).rev() {
        cum += score(
            cols[c],
            if c + 1 < num_cols {
                Some(cols[c + 1].0)
            } else {
                None
            },
        );
        if cum < best && c >= 1 && cols[c - 1] == (b'M', true) {
            best = cum;
            rcut = num_cols - c;
        }
    }

    if lcut == 0 && rcut == 0 {
        return None;
    }

    if lcut + rcut >= num_cols {
        return None;
    }

    let qlen = |span_cols: &[(u8, bool)]| span_cols.iter().filter(|col| col.0 != b'D').count();
    let tlen = |span_cols: &[(u8, bool)]| span_cols.iter().filter(|col| col.0 != b'I').count();
    let (left, core, right) = (
        &cols[..lcut],
        &cols[lcut..num_cols - rcut],
        &cols[num_cols - rcut..],
    );
    let (left_len, right_len, core_len) = (qlen(left), qlen(right), qlen(core));

    // A core shorter than half the read is rejected, because the clips would then outweigh the alignment.
    const MIN_CORE_SHARE_DEN: usize = 2;

    if core_len * MIN_CORE_SHARE_DEN < query.len() {
        return None;
    }

    let mut out = String::new();

    if left_len > 0 {
        out.push_str(&format!("{left_len}S"));
    }

    // The core is run-length encoded again, and its edit distance counts every mismatch and every gap base.
    let mut run_start = 0;
    let mut edit_dist = 0i32;

    while run_start < core.len() {
        let run_op = core[run_start].0;
        let mut run_end = run_start;

        while run_end < core.len() && core[run_end].0 == run_op {
            if run_op != b'M' || !core[run_end].1 {
                edit_dist += 1;
            }

            run_end += 1;
        }

        out.push_str(&format!("{}{}", run_end - run_start, run_op as char));
        run_start = run_end;
    }

    if right_len > 0 {
        out.push_str(&format!("{right_len}S"));
    }

    Some((out, toff + tlen(left), edit_dist))
}

/// Clipped bases that count as one more alignment event in [`cigar_events`].
pub(super) const CLIP_BASES_PER_EVENT: usize = 12;

/// Alignment events of a CIGAR (mismatches, gap openings, and clipped ends counted as one plus one per CLIP_BASES_PER_EVENT clipped bases), the error measure of the MAPQ cap.
/// A 10 bp indel is one event, not ten, because counting its length made reads carrying a true indel look like high-error reads and capped them near MAPQ 30, where pipelines filter.
pub(super) fn cigar_events(cig: &str, query: &[u8], text: &[u8], toff: usize) -> i32 {
    let (mut query_i, mut text_i, mut run_len, mut events) = (0usize, toff, 0usize, 0i32);

    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        match byte {
            b'M' => {
                for _ in 0..run_len {
                    if query_i < query.len()
                        && text_i < text.len()
                        && (query[query_i] != text[text_i]
                            || !matches!(text[text_i], b'A' | b'C' | b'G' | b'T'))
                    {
                        events += 1;
                    }

                    query_i += 1;
                    text_i += 1;
                }
            }
            // A gap of any length is one event.
            b'I' => {
                events += 1;
                query_i += run_len;
            }
            b'D' => {
                events += 1;
                text_i += run_len;
            }
            // A clipped end is one event plus one per 12 clipped bases, since a long clip is read sequence the placement does not explain.
            // An aligned part fitting a diverged repeat copy, while the whole read's locus was never seeded, was capped at MAPQ 41 or so.
            b'S' => {
                events += 1 + (run_len / CLIP_BASES_PER_EVENT) as i32;
                query_i += run_len;
            }
            _ => {}
        }

        run_len = 0;
    }

    events
}
