//! `optik map` maps reads with minimizer seeds.
//!
//! Short reads, under CHAIN_MIN_LEN bases, vote for diagonals with their minimizer hits.
//! Nearby votes merge into candidate loci, and each locus is verified by bit-parallel edit distance in include/bpm_edit.hpp, four loci per AVX2 register.
//! The best placement is aligned with gap-affine WFA, indels are left-aligned, and low-scoring ends are soft-clipped.
//! MAPQ is the posterior over the verified competitors plus the odds that the true locus was never seeded, computed by posterior_mapq_w and hidden_odds with constants fitted on simulated reads.
//!
//! Long reads and genome-scale queries have their anchors chained as minimap2 chains them in tools/src/chain.rs, the gaps between anchors are aligned with WFA or BiWFA, flanks are extended under a Z-drop-style rule, and query-disjoint chains become supplementary records (Li 2018, 2021).
//!
//! Paired-end reads come as two -q files.
//! Each mate is placed on its own, then the pair is chosen jointly as BWA-MEM pairs them (Li 2013), with mate rescue inside the insert-size window, as the "Paired-end reads" section explains.
//!
//! Records are SAM with NM, MD and optionally cs, written in input order at any thread count.

use crate::{SeedIdx, iterate_fasta, iterate_reads};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::time::Instant;

/// Traceback direction codes for the CIGAR DP matrix.
/// back[i][j] stores which DP recurrence produced the minimum score at cell (i,j).
const DIR_MATCH: u8 = 0; // diagonal step, query[i] against text[j] as an M op, match or mismatch
const DIR_DEL: u8 = 1; // vertical step, text[j] without a query counterpart, a D op
const DIR_INS: u8 = 2; // horizontal step, query[i] without a text counterpart, an I op

/// SAM flag bits (SAM specification, section 1.4).
const SAM_FLAG_PAIRED: u16 = 0x1;
const SAM_FLAG_PROPER_PAIR: u16 = 0x2;
const SAM_FLAG_UNMAPPED: u16 = 0x4;
const SAM_FLAG_MATE_UNMAPPED: u16 = 0x8;
const SAM_FLAG_REVERSE: u16 = 0x10;
const SAM_FLAG_MATE_REVERSE: u16 = 0x20;
const SAM_FLAG_READ1: u16 = 0x40;
const SAM_FLAG_READ2: u16 = 0x80;
const SAM_FLAG_SUPPLEMENTARY: u16 = 0x800;
const SAM_MAX_MAPQ: i32 = 60;

/// True when the WFA CIGAR must be bypassed in favour of the scalar unit-cost DP, which `OPTIK_FORCE_SCALAR_CIGAR` set to anything but "0" requests.
/// It is a safety switch and an A/B reference for the short-read path.
fn force_scalar_cigar() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| std::env::var_os("OPTIK_FORCE_SCALAR_CIGAR").is_some_and(|v| v != "0"))
}
/// Whether environment variable `name` is set, read once into `cell`, since a getenv per read showed in profiles.
fn env_flag(name: &'static str, cell: &'static std::sync::OnceLock<bool>) -> bool {
    *cell.get_or_init(|| std::env::var_os(name).is_some())
}
static DEBUG_CHAIN: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// Reverse complement of a DNA sequence, with any non-ACGT byte turned into N.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    // table lookup, no per-base branch (non-ACGT -> N)
    static RC: [u8; 256] = {
        let mut t = [b'N'; 256];
        t[b'A' as usize] = b'T';
        t[b'a' as usize] = b'T';
        t[b'C' as usize] = b'G';
        t[b'c' as usize] = b'G';
        t[b'G' as usize] = b'C';
        t[b'g' as usize] = b'C';
        t[b'T' as usize] = b'A';
        t[b't' as usize] = b'A';
        t
    };
    seq.iter().rev().map(|&b| RC[b as usize]).collect()
}

/// Semi-global unit-cost alignment CIGAR by dynamic programming, with the query end to end and the text's ends free.
/// Returns (CIGAR, 0-based start offset in the text).
/// This is the fallback when the wavefront aligner gives up, and the reference path under OPTIK_FORCE_SCALAR_CIGAR.
fn cigar(query: &[u8], text: &[u8]) -> (String, usize) {
    let q = query.len();
    let t = text.len();
    if q == 0 {
        return ("*".into(), 0);
    }

    let mut back = vec![vec![0u8; q + 1]; t + 1];
    back[0].fill(DIR_INS);

    let mut prev = (0..=q as i32).collect::<Vec<i32>>();

    let mut best_i = 0usize;
    let mut best_score = q as i32;

    for i in 1..=t {
        let mut cur = vec![0i32; q + 1];
        cur[0] = 0;
        back[i][0] = DIR_DEL;
        for j in 1..=q {
            let cost = if text[i - 1] == query[j - 1] { 0i32 } else { 1 };
            let diag = prev[j - 1] + cost;
            let up = prev[j] + 1;
            let left = cur[j - 1] + 1;
            if diag < up && diag < left {
                cur[j] = diag;
                back[i][j] = DIR_MATCH;
            } else if up <= left {
                cur[j] = up;
                back[i][j] = DIR_DEL;
            } else {
                cur[j] = left;
                back[i][j] = DIR_INS;
            }
        }
        if cur[q] <= best_score {
            best_score = cur[q];
            best_i = i;
        }
        prev = cur;
    }

    let mut ops: Vec<u8> = Vec::new();
    let (mut i, mut j) = (best_i as i32, q as i32);
    while j > 0 {
        match back[i as usize][j as usize] {
            DIR_MATCH => {
                ops.push(b'M');
                i -= 1;
                j -= 1;
            }
            DIR_DEL => {
                ops.push(b'D');
                i -= 1;
            }
            DIR_INS => {
                ops.push(b'I');
                j -= 1;
            }
            _ => unreachable!(),
        }
    }
    let text_start = i as usize;
    ops.reverse();

    let lead = ops.iter().take_while(|&&op| op == b'D').count();
    let tail = ops.iter().rev().take_while(|&&op| op == b'D').count();
    let core = &ops[lead..ops.len() - tail];

    let mut cigar = String::new();
    let mut n = 1u32;
    for w in core.windows(2) {
        if w[0] == w[1] {
            n += 1;
        } else {
            cigar.push_str(&n.to_string());
            cigar.push(w[0] as char);
            n = 1;
        }
    }
    if let Some(&last) = core.last() {
        cigar.push_str(&n.to_string());
        cigar.push(last as char);
    }

    if cigar.is_empty() {
        cigar = format!("{}M", q);
    }

    (cigar, text_start)
}

/// k chosen for noisy long reads, the same 15 as minimap2's map-ont.
const NOISY_K: usize = 15;
/// Long reads count as noisy above this per-base error, taken from base qualities, or below this share of their k = 21 minimizers found in the reference.
/// Both thresholds mark about 2% error, since (1 - 0.02)^21 = 0.65.
const NOISY_ERROR: f64 = 0.02;
const NOISY_PRESENCE: f64 = 0.65;
/// Reads sampled for the choice.
const NOISE_SAMPLE: usize = 1000;

/// What the first batch says about the reads' error rate.
/// They are noisy or accurate long reads, with the mean error from base qualities, mostly short reads, or long reads without qualities.
enum ReadNoise {
    Noisy(f64),
    Accurate(f64),
    Short,
    Unknown,
}

/// A read of a mapping batch, as name, seq, qual and, when paired, R2's seq and qual.
type MapRead = (String, Vec<u8>, Vec<u8>, Option<(Vec<u8>, Vec<u8>)>);

/// Whether the first batch is long reads and, if they carry base qualities, how noisy they are (mean per-base error 10^(-Q/10) over sampled long reads).
fn noisy_long_reads(batch: &[MapRead]) -> ReadNoise {
    let n_long = batch.iter().filter(|r| r.1.len() >= CHAIN_MIN_LEN).count();
    if n_long * 2 < batch.len() {
        return ReadNoise::Short; // mostly short reads keep the fitted short-read k
    }
    let long: Vec<_> = batch
        .iter()
        .filter(|r| r.1.len() >= CHAIN_MIN_LEN)
        .take(NOISE_SAMPLE)
        .collect();
    let (mut e, mut n) = (0.0f64, 0u64);
    for r in &long {
        if r.2.len() != r.1.len() {
            return ReadNoise::Unknown;
        } // FASTA carries no qualities
        for &q in &r.2 {
            e += 10f64.powf(-((q.saturating_sub(33)) as f64) / 10.0);
            n += 1;
        }
    }
    if n == 0 {
        return ReadNoise::Unknown;
    }
    let err = e / n as f64;
    if err > NOISY_ERROR {
        ReadNoise::Noisy(err)
    } else {
        ReadNoise::Accurate(err)
    }
}

/// Share of the sampled long reads' minimizers found in the index (or dropped there as too repetitive), about the share of error-free k-mers.
fn seed_presence(batch: &[MapRead], seeds: &SeedIdx, k: usize) -> f64 {
    let (mut found, mut total) = (0u64, 0u64);
    for r in batch
        .iter()
        .filter(|r| r.1.len() >= CHAIN_MIN_LEN)
        .take(NOISE_SAMPLE)
    {
        let (hits, n_min, n_rep) = seeds.collect_with_rep(&r.1);
        let len = r.1.len();
        let mut pos: Vec<usize> = hits
            .iter()
            .map(|h| {
                if h.rev != 0 {
                    len - k - h.qpos as usize
                } else {
                    h.qpos as usize
                }
            })
            .collect();
        pos.sort_unstable();
        pos.dedup();
        found += pos.len() as u64 + n_rep as u64;
        total += n_min as u64;
    }
    if total == 0 {
        1.0
    } else {
        found as f64 / total as f64
    }
}

/// The minimizer window for `-w 0`, read_len / 4 bounded below by 10 and above by max(k / 2, 10).
fn auto_tune_w(read_len: usize, k: usize) -> usize {
    let w = (k / 2).max(10);
    (read_len / 4).min(w).max(10)
}

/// Reads at or above this length take the chaining path, with anchor chaining, gap-affine WFA between anchors, no -e cap, and primary plus supplementary records for query-disjoint chains, genome-scale queries included.
const CHAIN_MIN_LEN: usize = 500;

/// Queries at least this long, such as assemblies and contigs, are mapped as genome-scale sequences, and record selection keeps every query-disjoint block.
const GENOME_SCALE_LEN: usize = 100_000;

/// Gap-affine penalties x, o and e for chained long-read alignment.
const CHAIN_PEN: (i32, i32, i32) = (4, 6, 2);

/// A flank whose affine score exceeds this many units per base is soft-clipped instead of aligned (random sequence scores ~3/base).
const FLANK_MAX_PER_BASE: f64 = 1.5;

/// Flanks are extended outward from the end anchor in chunks of this many read bases, each aligned under its own FLANK_MAX_PER_BASE cap, and the first chunk over its cap soft-clips the rest in a Z-drop-style stop.
/// WFA memory grows with the square of the score, so one cap over a whole 40 kb flank meant gigabytes per thread, while a chunk cap keeps it to megabytes.
/// Chunks gallop, starting at FLANK_START_CHUNK, doubling after an aligned chunk up to FLANK_CHUNK and halving after a failed one down to FLANK_MIN_CHUNK.
/// The clip then lands near where the read leaves the reference, and a foreign tail costs only a few small capped alignments.
const FLANK_CHUNK: usize = 1000;
const FLANK_MIN_CHUNK: usize = 32;
const FLANK_START_CHUNK: usize = 64;

/// Extension scoring of one aligned column, used to cut flanks and read ends at their best-scoring prefix.
/// A match scores +2, a mismatch -4 and a gap of length len -(4 + 2 len), minimap2's defaults, with the open charged on the gap's first column.
const EXT_MATCH: i64 = 2;
const EXT_MISMATCH: i64 = -4;
const EXT_GAP_OPEN: i64 = -4;
const EXT_GAP_EXT: i64 = -2;

/// Recompute the edit distance a CIGAR implies against (query, text), where a mismatch inside M costs 1 and so does every I and D base.
/// Returns None if the CIGAR overshoots either string, which means inconsistent bookkeeping.
fn cigar_ed_implied(cigar: &str, query: &[u8], text: &[u8], text_start: usize) -> Option<i32> {
    let (mut i, mut t) = (0usize, text_start);
    let mut ed = 0i32;
    let mut num: u64 = 0;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as u64;
            continue;
        }
        if num == 0 {
            return None;
        }
        let n = num as usize;
        num = 0;
        match ch {
            b'M' => {
                for _ in 0..n {
                    if i >= query.len() || t >= text.len() {
                        return None;
                    }
                    // a non-ACGT base (N) never matches, as in SAM's NM
                    if query[i] != text[t] || !matches!(text[t], b'A' | b'C' | b'G' | b'T') {
                        ed += 1;
                    }
                    i += 1;
                    t += 1;
                }
            }
            b'I' => {
                for _ in 0..n {
                    if i >= query.len() {
                        return None;
                    }
                    ed += 1;
                    i += 1;
                }
            }
            b'D' => {
                for _ in 0..n {
                    if t >= text.len() {
                        return None;
                    }
                    ed += 1;
                    t += 1;
                }
            }
            _ => return None,
        }
    }
    Some(ed)
}

/// Alignment columns of a CIGAR (M, =, X, I, D).
fn cigar_cols(cigar: &str) -> usize {
    let (mut n, mut tot) = (0usize, 0usize);
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        if matches!(ch, b'M' | b'=' | b'X' | b'I' | b'D') {
            tot += n;
        }
        n = 0;
    }
    tot
}

/// Whether a SAM record carries the unmapped flag or has no readable flag.
fn is_unmapped(rec: &str) -> bool {
    rec.split('\t')
        .nth(1)
        .and_then(|f| f.parse::<u16>().ok())
        .is_none_or(|f| f & SAM_FLAG_UNMAPPED != 0)
}

/// Append `n` copies of op to a run-length op list, merging with the tail.
fn push_op(ops: &mut Vec<(u32, u8)>, n: u32, op: u8) {
    if n == 0 {
        return;
    }
    if let Some(last) = ops.last_mut()
        && last.1 == op
    {
        last.0 += n;
        return;
    }
    ops.push((n, op));
}

/// Append every operation of a CIGAR string to a run-length op list.
fn push_cigar(ops: &mut Vec<(u32, u8)>, cig: &str) {
    let mut n = 0u32;
    for c in cig.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            push_op(ops, n, c);
            n = 0;
        }
    }
}

/// Extend a flank rightwards from an anchor, with `q` and `r` starting right after it.
/// The flank is aligned in galloping chunks (FLANK_START_CHUNK, doubled with each aligned chunk up to FLANK_CHUNK and halved on a failed one down to FLANK_MIN_CHUNK) with the text ends free, each chunk under its own FLANK_MAX_PER_BASE cap.
/// The resulting path is then cut at its best-scoring prefix under extension scoring (match +2, mismatch -4, gap -(4 + 2*len)), the Z-drop and X-drop rule of minimap2 and BLAST, so bases past the point where the read leaves the reference are soft-clipped rather than aligned as noise.
/// Per-chunk caps bound the wavefront memory, since a single cap over a 40 kb flank meant GBs per thread.
/// Returns (ops, read bases used, reference bases used, chain-penalty cost including FLANK_MAX_PER_BASE per clipped base).
fn extend_flank(q: &[u8], r: &[u8]) -> (Vec<(u32, u8)>, usize, usize, i64) {
    let mut ops: Vec<(u32, u8)> = Vec::new();
    let (mut qs, mut re) = (0usize, 0usize);
    // The chunk size gallops, starting small, doubling after each aligned chunk and halving after a failed one.
    // A foreign tail then costs a few tiny capped alignments, while a genuine long flank reaches FLANK_CHUNK in four steps.
    let mut cl = FLANK_START_CHUNK;
    while qs < q.len() && re < r.len() {
        // Near the contig end, a chunk longer than the remaining reference would have to absorb the overhang as an insertion and fail its cap, clipping genuine end bases.
        cl = cl.min(q.len() - qs).min(r.len() - re);
        // The text window covers the chunk plus 25% and 32 bases for indels.
        let span = (cl + cl / 4 + 32).min(r.len() - re);
        let text = &r[re..re + span];
        match crate::wfa_affine_align(
            &q[qs..qs + cl],
            text,
            text.len(),
            CHAIN_PEN,
            (cl as f64 * FLANK_MAX_PER_BASE) as i32,
        ) {
            Some(al) => {
                push_cigar(&mut ops, &al.cigar);
                re += al.tend;
                qs += cl;
                cl = (cl * 2).min(FLANK_CHUNK);
            }
            None if cl > FLANK_MIN_CHUNK => cl /= 2,
            None => break,
        }
    }
    // exact matches past the last chunk (a chunk that runs into an N gap or a foreign tail fails as a whole, dropping its genuine first bases)
    let mut run = 0u32;
    while qs < q.len() && re < r.len() && r[re] != b'N' && q[qs].eq_ignore_ascii_case(&r[re]) {
        qs += 1;
        re += 1;
        run += 1;
    }
    if run > 0 {
        push_op(&mut ops, run, b'M');
    }
    // best-scoring prefix of the path
    let same = |a: u8, b: u8| a.eq_ignore_ascii_case(&b) && b != b'N';
    let (mut sc, mut best) = (0i64, 0i64);
    let mut cut = (0usize, 0u32, 0usize, 0usize); // op index, bases of it, q, r
    let (mut qi, mut ri) = (0usize, 0usize);
    for (i, &(n, o)) in ops.iter().enumerate() {
        match o {
            b'M' => {
                for t in 0..n {
                    sc += if same(q[qi], r[ri]) {
                        EXT_MATCH
                    } else {
                        EXT_MISMATCH
                    };
                    qi += 1;
                    ri += 1;
                    if sc > best {
                        best = sc;
                        cut = (i, t + 1, qi, ri);
                    }
                }
            }
            b'I' => {
                sc += EXT_GAP_OPEN + EXT_GAP_EXT * n as i64;
                qi += n as usize;
            }
            _ => {
                sc += EXT_GAP_OPEN + EXT_GAP_EXT * n as i64;
                ri += n as usize;
            }
        }
    }
    let (ci, cn, qc, rc) = cut;
    let mut kept: Vec<(u32, u8)> = Vec::new();
    if best > 0 {
        kept.extend_from_slice(&ops[..ci]);
        kept.push((cn, b'M'));
    }
    // chain-penalty cost of the kept path plus the clip
    let (mut cost, mut qi, mut ri) = (0i64, 0usize, 0usize);
    for &(n, o) in &kept {
        match o {
            b'M' => {
                for _ in 0..n {
                    if !same(q[qi], r[ri]) {
                        cost += CHAIN_PEN.0 as i64;
                    }
                    qi += 1;
                    ri += 1;
                }
            }
            b'I' => {
                cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * n as i32) as i64;
                qi += n as usize;
            }
            _ => {
                cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * n as i32) as i64;
                ri += n as usize;
            }
        }
    }
    cost += ((q.len() - qc) as f64 * FLANK_MAX_PER_BASE) as i64;
    (kept, qc, rc, cost)
}

/// Windows longer than this many columns are not re-aligned by [`polish_gap_clusters`].
const POLISH_MAX_WINDOW: usize = 4000;
/// Gap runs shorter than this many columns do not start a cluster.
const POLISH_MIN_RUN: usize = 3;
/// Columns of exact context on each side of a re-aligned window.
const POLISH_MARGIN: usize = 16;

/// Re-align clusters of nearby gaps freely.
/// A chain's alignment is forced through its exact anchors, and inside a tandem repeat the anchors can sit on a neighbouring diagonal, so the path carries indels that cancel.
/// Examples are 9D 8D 27I where 1D 9I 1I 1I fits, or a 497 bp D and I one repeat unit apart where the copy aligns with mismatches.
/// Gap runs closer than max(64, 8 x the larger gap), or 2 x for gaps over 64, form a cluster, since a 20D and a 20I 73 bp apart cost 92 where 1I 1D and 4 mismatches cost 36.
/// Windows over POLISH_MAX_WINDOW columns are left alone, because they lie in satellite arrays, where aligning them doubled the HiFi mapping time on hg38 for placements that stay ambiguous.
/// Each cluster's window, with POLISH_MARGIN columns of margin, is aligned end to end with gap-affine WFA and replaced when that scores lower.
/// Only gap runs of at least POLISH_MIN_RUN columns form clusters, and gaps inside a window are re-aligned with it.
/// The 1-2 bp indels of a noisy read lie a few dozen bases apart, and joining them re-aligned whole ONT reads in 4 kb windows at 2.3x their mapping time.
/// Returns the new ops and the cost saved.
fn polish_gap_clusters(
    ops: &[(u32, u8)],
    q: &[u8],
    r: &[u8],
    qs: usize,
    rs: usize,
) -> (Vec<(u32, u8)>, i64) {
    let mut cols: Vec<u8> = Vec::new();
    for &(n, o) in ops {
        cols.extend(std::iter::repeat_n(o, n as usize));
    }
    let n = cols.len();
    let (mut qp, mut rp) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    for i in 0..n {
        qp[i + 1] = qp[i] + (cols[i] != b'D') as usize;
        rp[i + 1] = rp[i] + (cols[i] != b'I') as usize;
    }
    if qs + qp[n] > q.len() || rs + rp[n] > r.len() {
        return (ops.to_vec(), 0);
    }
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut i = 0;
    while i < n {
        if cols[i] == b'M' {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < n && cols[j] == cols[i] {
            j += 1;
        }
        if j - i >= POLISH_MIN_RUN {
            runs.push((i, j));
        }
        i = j;
    }
    let base_ok = |c: u8| matches!(c, b'A' | b'C' | b'G' | b'T');
    let (x, o, e) = (CHAIN_PEN.0 as i64, CHAIN_PEN.1 as i64, CHAIN_PEN.2 as i64);
    let mut out: Vec<u8> = Vec::with_capacity(n);
    let (mut done, mut saved, mut k) = (0usize, 0i64, 0usize);
    while k < runs.len() {
        let (s0, mut e0) = runs[k];
        let (mut j, mut maxg) = (k, e0 - s0);
        while j + 1 < runs.len() {
            let (s1, e1) = runs[j + 1];
            let g = maxg.max(e1 - s1);
            if s1 - e0 > if g <= 64 { 64usize.max(8 * g) } else { 2 * g } {
                break;
            }
            j += 1;
            e0 = e1;
            maxg = g;
        }
        let (ws, we) = (
            s0.saturating_sub(POLISH_MARGIN).max(done),
            (e0 + POLISH_MARGIN).min(n),
        );
        if j > k && we - ws <= POLISH_MAX_WINDOW {
            let mut cost = 0i64;
            for c in ws..we {
                match cols[c] {
                    b'M' => {
                        let (a, b) = (q[qs + qp[c]], r[rs + rp[c]]);
                        if a != b || !base_ok(a) || !base_ok(b) {
                            cost += x;
                        }
                    }
                    g => {
                        cost += e;
                        if c == ws || cols[c - 1] != g {
                            cost += o;
                        }
                    }
                }
            }
            let (qsub, rsub) = (&q[qs + qp[ws]..qs + qp[we]], &r[rs + rp[ws]..rs + rp[we]]);
            if !qsub.is_empty() && !rsub.is_empty() {
                let cap = (cost - 1).min(i32::MAX as i64) as i32;
                let al = if qsub.len().max(rsub.len()) > BIWFA_MIN_GAP {
                    crate::biwfa_align(qsub, rsub, CHAIN_PEN, cap)
                } else {
                    crate::wfa_affine_align(qsub, rsub, 0, CHAIN_PEN, cap)
                };
                if let Some(al) = al.filter(|al| (al.score as i64) < cost) {
                    out.extend_from_slice(&cols[done..ws]);
                    let mut m = 0usize;
                    for ch in al.cigar.bytes() {
                        if ch.is_ascii_digit() {
                            m = m * 10 + (ch - b'0') as usize;
                            continue;
                        }
                        let op = if matches!(ch, b'=' | b'X') { b'M' } else { ch };
                        out.extend(std::iter::repeat_n(op, m));
                        m = 0;
                    }
                    saved += cost - al.score as i64;
                    done = we;
                }
            }
        }
        k = j + 1;
    }
    if saved == 0 {
        return (ops.to_vec(), 0);
    }
    out.extend_from_slice(&cols[done..]);
    let mut res: Vec<(u32, u8)> = Vec::new();
    for &c in &out {
        push_op(&mut res, 1, c);
    }
    (res, saved)
}

/// minimap2's alignment scoring from v2.31 options.c, with match a, mismatch b, two-piece affine gaps min(q + e len, q2 + e2 len), an N scoring -sc_ambi = -1, and the Z-drop threshold.
/// ksw2 counts a drop past zdrop + e2 per column of diagonal shift since the best score.
#[derive(Clone, Copy)]
struct MmScore {
    a: i64,
    b: i64,
    q: i64,
    e: i64,
    q2: i64,
    e2: i64,
    zdrop: i64,
}
/// mm_mapopt_init for reads, where map-ont is the default
const MM_READS: MmScore = MmScore {
    a: 2,
    b: 4,
    q: 4,
    e: 2,
    q2: 24,
    e2: 1,
    zdrop: 400,
};
/// -x asm20, used for genome-scale queries
const MM_ASM20: MmScore = MmScore {
    a: 1,
    b: 4,
    q: 6,
    e: 2,
    q2: 26,
    e2: 1,
    zdrop: 200,
};

impl MmScore {
    /// Two-piece affine cost of a gap of `len`.
    fn gap(&self, len: i64) -> i64 {
        (self.q + self.e * len).min(self.q2 + self.e2 * len)
    }
    /// Score of one aligned column, where an N scores -1.
    fn base(&self, a: u8, b: u8) -> i64 {
        if a == b'N' || b == b'N' {
            -1
        } else if a.eq_ignore_ascii_case(&b) {
            self.a
        } else {
            -self.b
        }
    }
}

/// Local clipping of an aligned path, as minimap2 treats a chain.
/// Returns the column ranges [s, e) worth keeping as records with their scores under `sc`, minimap2's own scoring, where a gap column scores its increment of the two-piece gap cost, counted from the gap's first column in scan order.
/// The path is split wherever its score falls more than sc.zdrop plus sc.e2 per column of diagonal shift below the best so far.
/// Under this ksw2 Z-drop, one long indel never splits.
/// A divergent stretch inside the path, meaning anchors across inserted or rearranged sequence whose gaps are short enough to align uncapped, becomes a break between two records.
/// Each piece's ends are then cut at their best-scoring prefix and suffix, but never inward of its outermost anchors.
/// Between anchors the path stays global, as minimap2 fills a chain's gaps, and only what lies beyond them is clipped, as its end extension does.
/// A piece needs an anchor of its own, found through `anchor_q`, the ascending read offsets of the segment's anchors, and a score of at least one anchor, k a.
/// `qp` and `rp` are the read and reference offsets before each column.
#[allow(clippy::too_many_arguments)]
fn local_pieces(
    cols: &[u8],
    q: &[u8],
    r: &[u8],
    qp: &[usize],
    rp: &[usize],
    anchor_q: &[usize],
    k: usize,
    sc: MmScore,
) -> Vec<(usize, usize, i64)> {
    let n = cols.len();
    // position of each gap column in its run, from either side (1 = first)
    let (mut lf, mut lb) = (vec![0i64; n], vec![0i64; n]);
    for c in 0..n {
        if cols[c] != b'M' {
            lf[c] = if c > 0 && cols[c - 1] == cols[c] {
                lf[c - 1] + 1
            } else {
                1
            };
        }
    }
    for c in (0..n).rev() {
        if cols[c] != b'M' {
            lb[c] = if c + 1 < n && cols[c + 1] == cols[c] {
                lb[c + 1] + 1
            } else {
                1
            };
        }
    }
    let score = |c: usize, fwd: bool| -> i64 {
        if cols[c] == b'M' {
            return sc.base(q[qp[c]], r[rp[c]]);
        }
        let l = if fwd { lf[c] } else { lb[c] };
        sc.gap(l - 1) - sc.gap(l)
    };
    // 1. Z-drop splits
    let mut raw: Vec<(usize, usize)> = Vec::new();
    let (mut start, mut acc, mut best, mut bi) = (0usize, 0i64, 0i64, 0usize);
    for j in 0..n {
        acc += score(j, true);
        if acc > best {
            best = acc;
            bi = j + 1;
        } else {
            let shift = (qp[j + 1] - qp[bi]).abs_diff(rp[j + 1] - rp[bi]) as i64;
            if best - acc > sc.zdrop + sc.e2 * shift {
                if bi > start {
                    raw.push((start, bi));
                }
                start = j + 1;
                acc = 0;
                best = 0;
                bi = j + 1;
            }
        }
    }
    if n > start {
        raw.push((start, n));
    }
    // 2. ends cut at their best score, outside the outermost anchors
    let col_of = |qoff: usize| qp.partition_point(|&x| x < qoff);
    let mut out = Vec::new();
    for (s, e) in raw {
        let inside: Vec<usize> = anchor_q
            .iter()
            .copied()
            .filter(|&a| a >= qp[s] && a + k <= qp[e])
            .collect();
        let (Some(&fa), Some(&la)) = (inside.first(), inside.last()) else {
            continue;
        };
        let (lo, hi) = (col_of(fa).max(s), col_of(la + k).min(e));
        // Right end, the best prefix score over [hi, e]
        let mut acc: i64 = (s..hi).map(|c| score(c, true)).sum();
        let (mut ee, mut be) = (hi, acc);
        for c in hi..e {
            acc += score(c, true);
            if acc > be {
                be = acc;
                ee = c + 1;
            }
        }
        // Left start, the best suffix score ending at ee over [s, lo]
        let mut acc: i64 = (lo..ee).rev().map(|c| score(c, false)).sum();
        let (mut ss, mut bs) = (lo, acc);
        for c in (s..lo).rev() {
            acc += score(c, false);
            if acc > bs {
                bs = acc;
                ss = c;
            }
        }
        if bs >= k as i64 * sc.a {
            out.push((ss, ee, bs));
        }
    }
    out
}

/// Genome-scale records scoring below this under MM_ASM20 are dropped, as for minimap2 -x asm20's minimal peak DP score (v2.31 options.c min_dp_max = 200, where hit.c drops a record whose dp_max is lower).
const ASM_MIN_DP: i64 = 200;

/// A record's alignment score under minimap2 -x asm20's scoring, MM_ASM20.
/// local_pieces cuts records at their best-scoring ends with the same scoring, so the whole record's score stands for minimap2's peak, dp_max.
/// At 1 and -4 a column scores above zero only below 20% divergence.
/// Diverged homologous sequence, such as different genes of two strains, therefore scores under ASM_MIN_DP, and minimap2 -x asm20 reports no record for it.
fn asm_score(cigar: &str, q: &[u8], r: &[u8], pos: usize) -> i64 {
    let sc = MM_ASM20;
    let (mut qi, mut ri, mut n, mut tot) = (0usize, pos, 0usize, 0i64);
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        match ch {
            b'S' | b'H' => qi += n,
            b'M' | b'=' | b'X' => {
                for t in 0..n {
                    tot += sc.base(q[qi + t], r[ri + t]);
                }
                qi += n;
                ri += n;
            }
            b'I' => {
                tot -= sc.gap(n as i64);
                qi += n;
            }
            b'D' => {
                tot -= sc.gap(n as i64);
                ri += n;
            }
            _ => {}
        }
        n = 0;
    }
    tot
}

/// One aligned piece of a chain, holding the SAM fields of a record.
struct ChainAln {
    rid: usize,
    rev: bool,
    pos: usize,
    cigar: String,
    nm: i32,
    qs: usize,
    qe: usize,
    n_anchors: usize,
    /// gap-affine cost of the base-level alignment (+ clipped bases)
    cost: i64,
}

/// Gaps between anchors up to this long on both sides use the traceback WFA, whose score cap of 1.5 per base keeps its O(s^2) memory near 10 MB.
/// Longer gaps use BiWFA, with O(s) memory and the same optimum.
const BIWFA_MIN_GAP: usize = 1000;

/// Gaps up to this long are aligned without a score cap and never split a chain (worst case s ~ 4 * 500, a few MB of wavefronts).
const SPLIT_MIN_GAP: usize = 500;

/// Align a chain with exact anchors, gap-affine WFA between consecutive anchors end to end under a FLANK_MAX_PER_BASE per-base cap, and flanks extended outward by `extend_flank`.
/// A gap that fails its cap is not a plausible alignment, because the sequence is unrelated or rearranged, so the chain is split there.
/// Each side becomes its own record that extends into the gap and soft-clips the rest, as minimap2 splits at a Z-drop.
/// `local_pieces` then clips each segment's path, cutting its ends at their best score and splitting where a divergent stretch drops the score by more than ZDROP.
/// Returns the records in read order, or nothing if no anchor survives.
fn align_chain(
    c: &crate::chain::Chain,
    q: &[u8],
    rseq: &[u8],
    k: usize,
    sc: MmScore,
) -> Vec<ChainAln> {
    // Keep strictly non-overlapping anchors, each contributing k exact matches.
    let mut an: Vec<(usize, usize)> = Vec::with_capacity(c.anchors.len());
    for a in &c.anchors {
        let (r, qq) = (a.r as usize, a.q as usize);
        if let Some(&(pr, pq)) = an.last()
            && (r < pr + k || qq < pq + k)
        {
            continue;
        }
        if r + k > rseq.len() || qq + k > q.len() {
            continue;
        }
        an.push((r, qq));
    }
    if an.is_empty() {
        return Vec::new();
    }
    // core alignments of anchor runs, split at implausible gaps
    struct Seg {
        first: usize,
        last: usize,
        ops: Vec<(u32, u8)>,
        cost: i64,
    }
    let mut segs: Vec<Seg> = Vec::new();
    let mut cur = Seg {
        first: 0,
        last: 0,
        ops: Vec::new(),
        cost: 0,
    };
    for w in 0..an.len() {
        let (r, qq) = an[w];
        push_op(&mut cur.ops, k as u32, b'M');
        cur.last = w;
        if w + 1 == an.len() {
            break;
        }
        let (nr, nq) = an[w + 1];
        let (gq, gr) = (&q[qq + k..nq], &rseq[r + k..nr]);
        match (gq.is_empty(), gr.is_empty()) {
            (true, true) => {}
            (false, true) => {
                push_op(&mut cur.ops, gq.len() as u32, b'I');
                cur.cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * gq.len() as i32) as i64;
            }
            (true, false) => {
                push_op(&mut cur.ops, gr.len() as u32, b'D');
                cur.cost += (CHAIN_PEN.1 + CHAIN_PEN.2 * gr.len() as i32) as i64;
            }
            _ if gq.len().min(gr.len()) <= 32 && gq.len().abs_diff(gr.len()) > SPLIT_MIN_GAP => {
                // One tiny side against a long one is a long deletion or insertion.
                // The short side is aligned ends-free against the near end of the long side, and the rest becomes one gap that left_align places.
                let (short, long, del) = if gq.len() < gr.len() {
                    (gq, gr, true)
                } else {
                    (gr, gq, false)
                };
                let w = (short.len() + 64).min(long.len());
                let al = crate::wfa_affine_align(short, &long[..w], w, CHAIN_PEN, -1);
                let (cig, tend, sc) = match al {
                    Some(a) => (a.cigar, a.tend, a.score),
                    None => (String::new(), 0, 0),
                };
                let mut tmp: Vec<(u32, u8)> = Vec::new();
                if tend == 0 && !short.is_empty() {
                    // Fallback, with the short side as its own gap
                    push_op(&mut tmp, short.len() as u32, if del { b'I' } else { b'D' });
                } else {
                    push_cigar(&mut tmp, &cig);
                    if !del {
                        // aligned with the read's long side as text, so I and D swap
                        for o in tmp.iter_mut() {
                            o.1 = match o.1 {
                                b'I' => b'D',
                                b'D' => b'I',
                                x => x,
                            };
                        }
                    }
                }
                for (n, o) in tmp {
                    push_op(&mut cur.ops, n, o);
                }
                let rest = long.len() - if tend == 0 { 0 } else { tend };
                push_op(&mut cur.ops, rest as u32, if del { b'D' } else { b'I' });
                cur.cost += sc as i64 + (CHAIN_PEN.1 + CHAIN_PEN.2 * rest as i32) as i64;
            }
            _ => {
                // Short gaps are always aligned, exactly and in tiny memory.
                // Their per-base cost is noisy, since a few indels in a 20 bp gap of a 10%-error read exceed any per-base cap, and splitting there would cut genuine alignments.
                // Longer gaps are capped and split when the cap fails.
                // The cap allows the length difference as one indel, open plus extend per base, plus FLANK_MAX_PER_BASE over the part both sides share.
                // A 600 bp insertion with 80 reference bases around it then stays an insertion, written as one I in one record as minimap2 and SV callers expect, while two long sides that do not align still fail and split.
                // A flat 1.5 per base over the longer side split every insertion or deletion over about 500 bp, and on HG002 HiFi it moved 77 reads' primaries to the far side of a real insertion.
                let glen = gq.len().max(gr.len());
                let (shared, diff) = (gq.len().min(gr.len()), gq.len().abs_diff(gr.len()));
                let cap = if glen <= SPLIT_MIN_GAP {
                    -1
                } else {
                    CHAIN_PEN.1
                        + CHAIN_PEN.2 * diff as i32
                        + (shared as f64 * FLANK_MAX_PER_BASE) as i32
                        + CHAIN_PEN.1
                };
                let al = if glen > BIWFA_MIN_GAP {
                    crate::biwfa_align(gq, gr, CHAIN_PEN, cap)
                } else {
                    crate::wfa_affine_align(gq, gr, 0, CHAIN_PEN, cap)
                };
                match al {
                    Some(al) => {
                        push_cigar(&mut cur.ops, &al.cigar);
                        cur.cost += al.score as i64;
                    }
                    None => {
                        segs.push(std::mem::replace(
                            &mut cur,
                            Seg {
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
    let mut out = Vec::with_capacity(segs.len());
    for (i, sg) in segs.iter().enumerate() {
        let (r0, q0) = an[sg.first];
        let (rl, ql) = an[sg.last];
        // flank regions end at the neighbouring segment's anchor
        let (q_lo, r_lo) = if i == 0 {
            (0, 0)
        } else {
            let (pr, pq) = an[segs[i - 1].last];
            (pq + k, pr + k)
        };
        let (q_hi, r_hi) = if i + 1 == segs.len() {
            (q.len(), rseq.len())
        } else {
            let (nr, nq) = an[segs[i + 1].first];
            (nq, nr)
        };
        let mut ops: Vec<(u32, u8)> = Vec::new();
        let mut cost = sg.cost;
        // The left flank extends leftwards, i.e. rightwards over the reversed sequences.
        let lq = q0 - q_lo;
        let ql_rev: Vec<u8> = q[q_lo..q0].iter().rev().copied().collect();
        let rl_rev: Vec<u8> = rseq[r0.saturating_sub(2 * lq + 64).max(r_lo)..r0]
            .iter()
            .rev()
            .copied()
            .collect();
        let (fops, qc, rc, fcost) = extend_flank(&ql_rev, &rl_rev);
        for &(n, o) in fops.iter().rev() {
            push_op(&mut ops, n, o);
        }
        cost += fcost;
        let (ref_start, clip_l) = (r0 - rc, q0 - qc);
        for &(n, o) in &sg.ops {
            push_op(&mut ops, n, o);
        }
        // right flank
        let (qs, re) = (ql + k, rl + k);
        let rq = q_hi - qs;
        let (fops, qc, _, fcost) =
            extend_flank(&q[qs..q_hi], &rseq[re..(re + 2 * rq + 64).min(r_hi)]);
        for (n, o) in fops {
            push_op(&mut ops, n, o);
        }
        cost += fcost;
        let clip_r = q.len() - (qs + qc);
        let (ops, saved) = polish_gap_clusters(&ops, q, rseq, clip_l, ref_start);
        cost -= saved;
        // Local clipping keeps the path's best-scoring pieces (local_pieces).
        let mut cols: Vec<u8> = Vec::new();
        for &(n, o) in &ops {
            cols.extend(std::iter::repeat_n(o, n as usize));
        }
        let nc = cols.len();
        let (mut qp, mut rp) = (vec![clip_l; nc + 1], vec![ref_start; nc + 1]);
        for i in 0..nc {
            qp[i + 1] = qp[i] + (cols[i] != b'D') as usize;
            rp[i + 1] = rp[i] + (cols[i] != b'I') as usize;
        }
        if qp[nc] != q.len() - clip_r || rp[nc] > rseq.len() {
            continue;
        }
        let anchors = &an[sg.first..=sg.last];
        let anchor_q: Vec<usize> = anchors.iter().map(|&(_, aq)| aq).collect();
        let pieces: Vec<(usize, usize)> = local_pieces(&cols, q, rseq, &qp, &rp, &anchor_q, k, sc)
            .into_iter()
            .map(|(s, e, _)| (s, e))
            .collect();
        let pieces = if pieces.is_empty() {
            vec![(0, nc)]
        } else {
            pieces
        };
        // Columns outside every piece have their path cost replaced by the clip charge per read base, as extend_flank charges a clipped flank.
        if pieces != [(0, nc)] {
            let mut keep = vec![false; nc];
            for &(s, e) in &pieces {
                keep[s..e].fill(true);
            }
            let mut dropped = 0usize;
            for i in (0..nc).filter(|&i| !keep[i]) {
                cost -= match cols[i] {
                    b'M' => {
                        if q[qp[i]].eq_ignore_ascii_case(&rseq[rp[i]]) && rseq[rp[i]] != b'N' {
                            0
                        } else {
                            CHAIN_PEN.0 as i64
                        }
                    }
                    g => {
                        CHAIN_PEN.2 as i64
                            + if i == 0 || cols[i - 1] != g {
                                CHAIN_PEN.1 as i64
                            } else {
                                0
                            }
                    }
                };
                dropped += (cols[i] != b'D') as usize;
            }
            cost += (dropped as f64 * FLANK_MAX_PER_BASE) as i64;
        }
        for (pi, &(s, e)) in pieces.iter().enumerate() {
            let (qs0, qe0, rs0) = (qp[s], qp[e], rp[s]);
            let mut pops: Vec<(u32, u8)> = Vec::new();
            for &o in &cols[s..e] {
                push_op(&mut pops, 1, o);
            }
            let core: String = pops
                .iter()
                .map(|(n, o)| format!("{n}{}", *o as char))
                .collect();
            let core = left_align(&core, &q[qs0..qe0], rseq, rs0);
            let Some(nm) = cigar_ed_implied(&core, &q[qs0..qe0], rseq, rs0) else {
                continue;
            };
            let mut cigar = String::new();
            if qs0 > 0 {
                cigar.push_str(&format!("{qs0}S"));
            }
            cigar.push_str(&core);
            if qe0 < q.len() {
                cigar.push_str(&format!("{}S", q.len() - qe0));
            }
            let n_anchors = anchors
                .iter()
                .filter(|&&(_, aq)| aq >= qs0 && aq + k <= qe0)
                .count();
            out.push(ChainAln {
                rid: c.rid as usize,
                rev: c.rev,
                pos: rs0,
                cigar,
                nm,
                qs: qs0,
                qe: qe0,
                n_anchors,
                cost: if pi == 0 { cost } else { 0 },
            });
        }
    }
    out
}

/// Ranking cost of a chain's alignment, summed over its segments.
/// Each segment's flank cost already charges the bases it clips.
fn chain_cost(segs: &[ChainAln]) -> Option<i64> {
    if segs.is_empty() {
        None
    } else {
        Some(segs.iter().map(|a| a.cost).sum())
    }
}

/// Default per-read edit budget, used without `-e`, as a share of the read length.
/// A fixed cap rejects every read beyond it regardless of length, while divergence scales with length.
/// 15% of the read keeps alignments far from random-sequence identity, since an unrelated 150 bp read reaches a best semi-global identity of about 55-60%.
/// Candidates come only from exact seed hits, and 0 of 50k random reads map.
const AUTO_ED_RATE: f64 = 0.15;
/// Smallest automatic edit budget.
const MIN_AUTO_ED: i32 = 5;

/// The edit budget of a read of `len` bases, `max_ed` if given (>= 0) and otherwise AUTO_ED_RATE of the length, at least MIN_AUTO_ED.
fn read_max_ed(max_ed: i32, len: usize) -> i32 {
    if max_ed >= 0 {
        max_ed
    } else {
        ((len as f64 * AUTO_ED_RATE).ceil() as i32).max(MIN_AUTO_ED)
    }
}

/// Chaining thresholds.
/// Long reads use the defaults.
/// The short-read rescue of reads over -e accepts sparser chains, because few error-free minimizers survive in a high-error 150 bp read, which is also why minimap2's sr preset uses -n2 -m20.
/// It then requires the alignment itself to look genuine.
#[derive(Clone, Copy)]
struct ChainParams {
    min_score: f64,
    min_anchors: usize,
    rescue: bool,
}
const LONG_CHAIN: ChainParams = ChainParams {
    min_score: crate::chain::MIN_CHAIN_SCORE,
    min_anchors: crate::chain::MIN_CHAIN_ANCHORS,
    rescue: false,
};
const RESCUE_CHAIN: ChainParams = ChainParams {
    min_score: 25.0,
    min_anchors: 2,
    rescue: true,
};
/// Rescue acceptance thresholds, the aligned read fraction and the edit density
const RESCUE_MIN_ALIGNED: f64 = 0.8;
const RESCUE_MAX_NM_RATE: f64 = AUTO_ED_RATE;
/// Split rescue limits, each part's minimum length, the allowed query overlap between parts, and the most parts emitted
const RESCUE_MIN_PART: usize = 40;
const RESCUE_PART_OVERLAP: usize = 10;
const RESCUE_MAX_PARTS: usize = 3;
/// MAPQ ceiling of a rescued short read, set to 3.
/// The chain MAPQ is fitted to long reads with hundreds of anchors.
/// A rescued read is one the short-read aligner could not place within -e, its few anchors leave most repeat copies unchained, and such reads were wrong 54% of the time at chain MAPQ 11-35 on simulated reads at 0.5-8% error from hg38 and repeat-rich genomes.
const RESCUE_MAPQ_MAX: u8 = 3;

/// Chains within this share of the best chain's score over the same read span are near-tied, and up to MAX_TIED_CHAINS of them are aligned and ranked by alignment cost.
const NEAR_TIE_SHARE: f64 = 0.9;
const MAX_TIED_CHAINS: usize = 5;
/// Records kept per read, primary plus supplementaries, and the share of the best chain's score a supplementary chain of a read needs.
const MAX_READ_RECORDS: usize = 6;
const SUPPLEMENTARY_MIN_SHARE: f64 = 0.3;
/// A genome-scale query keeps up to one record per this many query bases.
const GENOME_BASES_PER_RECORD: usize = 10_000;

/// Chaining-based mapping, for reads of at least CHAIN_MIN_LEN and for short reads over their edit budget.
/// Seeds become anchor chains per (reference, strand), then a primary and query-disjoint supplementary records.
fn map_one_read_chain(
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
        hits.sort_unstable_by_key(|h| (h.rid, h.rev, h.rpos, h.qpos));
        hits.dedup_by_key(|h| (h.rid, h.rev, h.rpos, h.qpos));
        let mut i = 0;
        let mut grp: Vec<Anchor> = Vec::new();
        while i < hits.len() {
            let (rid, rev) = (hits[i].rid, hits[i].rev);
            grp.clear();
            while i < hits.len() && hits[i].rid == rid && hits[i].rev == rev {
                grp.push(Anchor {
                    r: hits[i].rpos,
                    q: hits[i].qpos,
                });
                i += 1;
            }
            chain_group(
                rid,
                rev != 0,
                &grp,
                ctx.k,
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
    chains.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.rid.cmp(&b.rid))
            .then(a.anchors[0].r.cmp(&b.anchors[0].r))
    });
    // query spans in FORWARD read coordinates, for overlap tests
    let len = seq.len() as u32;
    let kk = ctx.k as u32;
    let fwd_span = |c: &Chain| {
        let (a, b) = c.q_span();
        let (s, e) = (a, b + kk);
        if c.rev { (len - e, len - s) } else { (s, e) }
    };
    let overlap = |x: (u32, u32), y: (u32, u32)| {
        let o = x.1.min(y.1).saturating_sub(x.0.max(y.0));
        o * 2 > (x.1 - x.0).min(y.1 - y.0)
    };
    // Near-tied chains over the same read span are repeat copies that chain scores cannot separate, but base-level alignment can, so up to 5 are aligned and the lowest gap-affine cost becomes the primary.
    // Genome-scale queries, i.e. assemblies, are scored as minimap2 -x asm20 scores them, and reads with its defaults, the local clip and asm20's minimal score ASM_MIN_DP.
    let genome_scale = seq.len() >= GENOME_SCALE_LEN;
    let mm_sc = if genome_scale { MM_ASM20 } else { MM_READS };
    let rc = revcomp(seq);
    // Best and second-best base-level cost over the primary's read span, kept as MAPQ evidence when near-tied copies were aligned
    let mut tie_dp: Option<(f64, f64)> = None;
    {
        let top = chains[0].score;
        let span0 = fwd_span(&chains[0]);
        let cand: Vec<usize> = (0..chains.len())
            .filter(|&j| {
                chains[j].score >= NEAR_TIE_SHARE * top && overlap(fwd_span(&chains[j]), span0)
            })
            .take(MAX_TIED_CHAINS)
            .collect();
        if env_flag("OPTIK_DEBUG_CHAIN", &DEBUG_CHAIN) {
            eprintln!(
                "[chain] {read_name} len={} rep_len={rep_len} chains={} near_tied={}",
                seq.len(),
                chains.len(),
                cand.len()
            );
            for c in chains.iter().take(8) {
                let q: &[u8] = if c.rev { &rc } else { seq };
                let cost = chain_cost(&align_chain(
                    c,
                    q,
                    &ctx.ref_seqs[c.rid as usize],
                    ctx.k,
                    mm_sc,
                ));
                eprintln!(
                    "[chain] {read_name} rid={} rev={} score={:.1} anchors={} r={}..{} q={:?} cost={:?}",
                    c.rid,
                    c.rev,
                    c.score,
                    c.anchors.len(),
                    c.anchors[0].r,
                    c.anchors[c.anchors.len() - 1].r,
                    fwd_span(c),
                    cost
                );
            }
        }
        if cand.len() > 1 {
            let costs: Vec<(i64, usize)> = cand
                .iter()
                .filter_map(|&j| {
                    let c = &chains[j];
                    let q: &[u8] = if c.rev { &rc } else { seq };
                    chain_cost(&align_chain(
                        c,
                        q,
                        &ctx.ref_seqs[c.rid as usize],
                        ctx.k,
                        mm_sc,
                    ))
                    .map(|cost| (cost, j))
                })
                .collect();
            if let Some(&(min_cost, _)) = costs.iter().min() {
                // Equal-cost copies (a read inside an identical repeat such as an rRNA operon) are chosen between pseudo-randomly, seeded by the read, as on the short-read path.
                // Always taking the first left the other copies without long-read coverage.
                let mut tied: Vec<usize> = costs
                    .iter()
                    .filter(|&&(c, _)| c == min_cost)
                    .map(|&(_, j)| j)
                    .collect();
                tied.sort_by_key(|&j| (chains[j].rid, chains[j].anchors[0].r, chains[j].rev));
                let j = tied[(read_hash(seq) % tied.len() as u64) as usize];
                // MAPQ evidence is the alignment score of the chosen copy and of the best other one, on minimap2's scale (match +2), from the same cost that ranked them, as 2 x read length - 1.5 x cost.
                // A mismatch costs 4 here and 6 score units there, and clipped bases are charged in the cost.
                let dp = |cost: i64| 2.0 * seq.len() as f64 - 1.5 * cost as f64;
                let d2 = costs.iter().filter(|x| x.1 != j).map(|x| x.0).min();
                if let Some(c2) = d2 {
                    tie_dp = Some((dp(min_cost), dp(c2)));
                }
                chains.swap(0, j);
            }
        }
    }
    // The primary and query-disjoint supplementaries are emitted.
    // Overlapping chains are secondary, not emitted, and feed MAPQ.
    // A read has one locus and a few chimeric parts, so query-disjoint chains within 30% of the best score are kept, at most 6.
    // A genome-scale query, such as an assembly in BWA-MEM's MG1655-vs-536 test, is a series of collinear blocks whose sizes bear no relation to the largest one.
    // Every query-disjoint chain that passed chaining's own thresholds is therefore kept, up to one per 10 kb of query, the way minimap2 keeps all query-disjoint primaries.
    let max_kept = if genome_scale {
        (seq.len() / GENOME_BASES_PER_RECORD).max(MAX_READ_RECORDS)
    } else {
        MAX_READ_RECORDS
    };
    let mut kept: Vec<usize> = vec![0];
    for (ci, c) in chains.iter().enumerate().skip(1) {
        if kept.len() >= max_kept {
            break;
        }
        if !genome_scale && c.score < SUPPLEMENTARY_MIN_SHARE * chains[0].score {
            break;
        }
        if kept
            .iter()
            .all(|&j| !overlap(fwd_span(c), fwd_span(&chains[j])))
        {
            kept.push(ci);
        }
    }
    let sum_kept: f64 = kept.iter().map(|&j| chains[j].score).sum();
    let mut alns: Vec<(ChainAln, u8)> = Vec::new();
    for &ci in &kept {
        let c = &chains[ci];
        let s2 = chains
            .iter()
            .enumerate()
            .filter(|(j, o)| *j != ci && overlap(fwd_span(o), fwd_span(c)))
            .map(|(_, o)| o.score)
            .fold(0.0f64, f64::max);
        let q: &[u8] = if c.rev { &rc } else { seq };
        let mut segs = align_chain(c, q, &ctx.ref_seqs[c.rid as usize], ctx.k, mm_sc);
        if genome_scale {
            // diverged homologous blocks go unreported, per minimap2 -x asm20's min_dp_max
            segs.retain(|a| asm_score(&a.cigar, q, &ctx.ref_seqs[a.rid], a.pos) >= ASM_MIN_DP);
            if segs.is_empty() {
                continue;
            }
        }
        let n_anchors: usize = segs.iter().map(|a| a.n_anchors).sum();
        let mq = if prm.rescue {
            mapq(c.score, s2, n_anchors).min(RESCUE_MAPQ_MAX)
        } else {
            let (cols, edits) = segs.iter().fold((0usize, 0i64), |(c0, e0), a| {
                (c0 + cigar_cols(&a.cigar), e0 + a.nm as i64)
            });
            let identity = if cols > 0 {
                ((cols as i64 - edits).max(0) as f64) / cols as f64
            } else {
                0.0
            };
            mapq_long(&MapqEvidence {
                s1: c.score,
                s2,
                n_anchors,
                sum_kept,
                rep_len,
                identity,
                dp: if ci == 0 { tie_dp } else { None },
            })
        };
        // The primary is the longest segment of the best chain, and a split chain's other segments become supplementaries.
        if alns.is_empty()
            && let Some(p) =
                (0..segs.len()).max_by_key(|&i| (segs[i].qe - segs[i].qs, std::cmp::Reverse(i)))
        {
            segs.swap(0, p);
        }
        for a in segs {
            alns.push((a, mq));
        }
    }
    if alns.is_empty() {
        return unmapped();
    }
    if prm.rescue {
        // Either one near-end-to-end, low-divergence record or, for a read across a contig end or a breakpoint, a primary plus supplementaries whose query-disjoint parts are each long and low-divergence and together cover the read.
        let good = |a: &ChainAln| {
            let aligned = (a.qe - a.qs) as f64;
            a.nm as f64 <= RESCUE_MAX_NM_RATE * aligned
        };
        let need = RESCUE_MIN_ALIGNED * seq.len() as f64;
        let a0 = &alns[0].0;
        if (a0.qe - a0.qs) as f64 >= need && good(a0) {
            alns.truncate(1);
        } else {
            let n = seq.len();
            let fwd = |a: &ChainAln| {
                if a.rev {
                    (n - a.qe, n - a.qs)
                } else {
                    (a.qs, a.qe)
                }
            };
            let mut parts: Vec<(ChainAln, u8)> = Vec::new();
            for (a, mq) in alns.drain(..) {
                let (s, e) = fwd(&a);
                if e - s < RESCUE_MIN_PART || !good(&a) {
                    continue;
                }
                // strictly query-disjoint, allowing a few bases of shared overlap
                if parts.iter().any(|(p, _)| {
                    let (ps, pe) = fwd(p);
                    e.min(pe).saturating_sub(s.max(ps)) > RESCUE_PART_OVERLAP
                }) {
                    continue;
                }
                parts.push((a, mq));
                if parts.len() == RESCUE_MAX_PARTS {
                    break;
                }
            }
            let mut iv: Vec<(usize, usize)> = parts.iter().map(|(a, _)| fwd(a)).collect();
            iv.sort_unstable();
            let (mut covered, mut end) = (0usize, 0usize);
            for (s, e) in iv {
                let s = s.max(end);
                if e > s {
                    covered += e - s;
                    end = e;
                }
            }
            if parts.len() < 2 || (covered as f64) < need {
                return unmapped();
            }
            alns = parts;
        }
    }
    let sa_entry = |a: &ChainAln, mq: u8| {
        format!(
            "{},{},{},{},{},{};",
            ref_names[a.rid],
            a.pos + 1,
            if a.rev { '-' } else { '+' },
            a.cigar,
            mq,
            a.nm
        )
    };
    alns.iter()
        .enumerate()
        .map(|(i, (a, mq))| {
            let flag = (if a.rev { SAM_FLAG_REVERSE } else { 0 })
                | (if i > 0 { SAM_FLAG_SUPPLEMENTARY } else { 0 });
            let sseq = if a.rev {
                std::str::from_utf8(&rc).unwrap_or("*")
            } else {
                seq_s
            };
            let sbytes: &[u8] = if a.rev { &rc } else { seq };
            let tags = diff_tag_fields(&a.cigar, sbytes, &ctx.ref_seqs[a.rid], a.pos, ctx.cs);
            // Supplementary records hard-clip and carry only their aligned bases (the default of minimap2 and bwa), since a soft-clipped copy of the whole query per record makes a genome-scale query's SAM grow with records x query length.
            let (cig, sseq) = if i > 0 {
                hard_clip(&a.cigar, sseq)
            } else {
                (a.cigar.clone(), sseq)
            };
            let mut line = format!(
                "{read_name}\t{flag}\t{}\t{}\t{mq}\t{cig}\t*\t0\t0\t{sseq}\t*\tNM:i:{}{tags}",
                ref_names[a.rid],
                a.pos + 1,
                a.nm
            );
            let sa: String = alns
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, (o, omq))| sa_entry(o, *omq))
                .collect();
            if !sa.is_empty() {
                line.push_str(&format!("\tSA:Z:{sa}"));
            }
            line
        })
        .collect()
}

/// The CIGAR's leading and trailing hard-clip lengths, 0 where absent.
fn hard_clips(cig: &str) -> (usize, usize) {
    let b = cig.as_bytes();
    let n = b.iter().take_while(|c| c.is_ascii_digit()).count();
    let lead = if n > 0 && b.get(n) == Some(&b'H') {
        cig[..n].parse().unwrap_or(0)
    } else {
        0
    };
    let mut trail = 0;
    if b.len() > 1 && b[b.len() - 1] == b'H' {
        let mut j = b.len() - 1;
        while j > 0 && b[j - 1].is_ascii_digit() {
            j -= 1;
        }
        if j < b.len() - 1 && !(j == 0 && lead > 0) {
            trail = cig[j..b.len() - 1].parse().unwrap_or(0);
        }
    }
    (lead, trail)
}

/// `cig` with its end soft clips turned into hard clips, and `seq` cut to the bases the CIGAR still consumes.
fn hard_clip<'a>(cig: &str, seq: &'a str) -> (String, &'a str) {
    let b = cig.as_bytes();
    let digits = |from: usize| b[from..].iter().take_while(|c| c.is_ascii_digit()).count();
    let mut lead = 0usize;
    let mut body = 0..b.len();
    let n = digits(0);
    if n > 0 && b.get(n) == Some(&b'S') {
        lead = cig[..n].parse().unwrap_or(0);
        body.start = n + 1;
    }
    let mut trail = 0usize;
    if b.last() == Some(&b'S') && body.end > body.start + 1 {
        let mut j = body.end - 1;
        while j > body.start && b[j - 1].is_ascii_digit() {
            j -= 1;
        }
        if j < body.end - 1 {
            trail = cig[j..body.end - 1].parse().unwrap_or(0);
            body.end = j;
        }
    }
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

/// Short-read mapping below CHAIN_MIN_LEN, with a single best record per read.
///
/// Seeds are hashed-order minimizers with strand.
/// Each reference hit votes for an alignment start on its diagonal, ref_pos - read_pos on the hit's strand.
/// Nearby starts merge into candidate loci, and each locus is verified once with the bit-parallel edit distance over a window of read_len + 2*pad around the diagonal, where pad = max_ed + 8, instead of re-aligning a window of 5x the read length for every seed hit.
/// MAPQ compares the best and second-best distinct loci and is 0 on a tie.
fn map_one_read_short(
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
struct ShortAln {
    cig: String,
    offset: usize,
    nm: i32,
    events: i32,
}

/// Aligns a short read at its verified placement, with a gap-affine CIGAR in the placement's window, left-aligned indels and low-scoring ends clipped.
fn short_align(seq: &[u8], rc: &[u8], hit: &ShortHit, ctx: &AlignCtx<'_>) -> ShortAln {
    let align_seq: &[u8] = if hit.rev { rc } else { seq };
    let rw = &ctx.ref_seqs[hit.rid][hit.win_start..hit.win_end];
    // Gap-affine CIGAR from optik's own wavefront aligner with minimap2-like short-read penalties x=4, o=6, e=2.
    // The unit-cost scalar DP stays as the fallback and as the OPTIK_FORCE_SCALAR_CIGAR reference.
    let (cig, offset) = if force_scalar_cigar() {
        cigar(align_seq, rw)
    } else {
        crate::wfa_affine_cigar(align_seq, rw).unwrap_or_else(|| cigar(align_seq, rw))
    };
    // Indels left-aligned (the SAM/VCF convention callers and pileups expect)
    let cig = left_align(&cig, align_seq, rw, offset);
    // Read ends that score below zero under extension scoring are soft-clipped, as the long-read flanks are and as minimap2 and BWA do.
    // An end-to-end CIGAR writes the last bases of a read that runs into an indel as mismatches, or starts the read inside deleted sequence, and variant callers count those bases as SNPs.
    let (cig, offset, nm) = match clip_ends(&cig, align_seq, rw, offset) {
        Some((full, toff, nm)) => (full, toff, nm),
        None => {
            // NM is derived from the emitted CIGAR so the record is self-consistent.
            let nm = cigar_ed_implied(&cig, align_seq, rw, offset).unwrap_or(hit.ed);
            (cig, offset, nm)
        }
    };
    let events = cigar_events(&cig, align_seq, rw, offset);
    ShortAln {
        cig,
        offset,
        nm,
        events,
    }
}

/// The SAM record of a short read at a verified placement.
#[allow(clippy::too_many_arguments)]
fn short_format(
    read_name: &str,
    seq: &[u8],
    rc: &[u8],
    hit: &ShortHit,
    a: &ShortAln,
    mapq: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> String {
    let align_seq: &[u8] = if hit.rev { rc } else { seq };
    let flag = if hit.rev { SAM_FLAG_REVERSE } else { 0 };
    // SAM SEQ is the read as aligned to the forward reference strand.
    let out_seq = std::str::from_utf8(align_seq).unwrap_or("*");
    let tags = diff_tag_fields(
        &a.cig,
        align_seq,
        &ctx.ref_seqs[hit.rid],
        hit.win_start + a.offset,
        ctx.cs,
    );
    format!(
        "{read_name}\t{flag}\t{}\t{}\t{mapq}\t{}\t*\t0\t0\t{out_seq}\t*\tNM:i:{}{tags}",
        ref_names[hit.rid],
        hit.win_start + a.offset + 1,
        a.cig,
        a.nm
    )
}

/// The SAM record of a short read at its verified placement `hit`, with a MAPQ from its own evidence that `mapq_cap` caps.
fn short_record(
    read_name: &str,
    seq: &[u8],
    hit: &ShortHit,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
    mapq_cap: Option<i32>,
) -> String {
    let rc = revcomp(seq);
    let a = short_align(seq, &rc, hit, ctx);
    let mut mapq = final_mapq(hit, a.events, seq.len(), ctx.rho);
    if let Some(c) = mapq_cap {
        mapq = mapq.min(c);
    }
    let mut rec = short_format(read_name, seq, &rc, hit, &a, mapq, ctx, ref_names);
    if mapq_debug() {
        let c: Vec<String> = hit.comp.iter().map(|x| x.to_string()).collect();
        rec.push_str(&format!(
            "\tXD:Z:{},{},{},{},{},{},{},{}:{}",
            hit.ed,
            a.events,
            hit.seed_ev[0],
            hit.seed_ev[1],
            hit.seed_ev[2],
            hit.seed_ev[3],
            hit.seed_ev[4],
            c.len(),
            c.join(";")
        ));
    }
    rec
}

/// Left-align every indel of an M/I/D CIGAR of `q` against `t` from `toff`.
/// Each gap moves to the leftmost position within the match run before it at which the columns it passes have no more mismatches than before.
/// The gap's cost does not depend on its position, so the score stays the same or gets better.
/// This covers the normalization of GATK LeftAlignIndels and bcftools norm, shifting through a repeat, and also equal-score placements that move a mismatch.
/// Equivalent placements of an indel thus collapse onto the leftmost, where pileup callers, GIAB and BWA-MEM put it.
/// A mismatch next to the gap moves with it, so a read with a 1-base deletion and an A>T SNP in an A run is written as deletion then SNP, as GIAB and BWA-MEM write it.
/// Writing a mismatch beside the deletion at the SNP's place hid the SNP and cost 148 calls on HG002 chr20.
/// An unparseable CIGAR is returned unchanged.
fn left_align(cig: &str, q: &[u8], t: &[u8], toff: usize) -> String {
    let mut ops: Vec<u8> = Vec::with_capacity(q.len() + 16);
    let mut n = 0usize;
    for ch in cig.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        if !matches!(ch, b'M' | b'I' | b'D') {
            return cig.to_string();
        }
        ops.extend(std::iter::repeat_n(ch, n));
        n = 0;
    }
    if !ops.iter().any(|&o| o != b'M') {
        return cig.to_string();
    }
    // Query and text position before each column, in one pass.
    // Shifting a gap only rearranges columns before its end, so the positions of every later column stay valid and the gaps are normalized left to right.
    let mut p: Vec<(usize, usize)> = Vec::with_capacity(ops.len());
    {
        let (mut qi, mut ti) = (0usize, toff);
        for &o in &ops {
            p.push((qi, ti));
            if o != b'D' {
                qi += 1;
            }
            if o != b'I' {
                ti += 1;
            }
        }
        if qi > q.len() || ti > t.len() {
            return cig.to_string();
        }
    }
    let base_ok = |c: u8| matches!(c, b'A' | b'C' | b'G' | b'T');
    let mut a = 0usize;
    while a < ops.len() {
        if ops[a] == b'M' {
            a += 1;
            continue;
        }
        let op = ops[a];
        let mut b = a;
        while b < ops.len() && ops[b] == op {
            b += 1;
        }
        let len = b - a;
        // The run [lo, lo + len) starts at query qs and text ts.
        // Shifting it left by j re-pairs the j columns before it, so a deletion pairs q[qs - i] with t[ts + len - i] instead of t[ts - i], and an insertion pairs t[ts - i] with q[qs + len - i] instead of q[qs - i].
        // The gap cost is unchanged, so the leftmost shift whose re-paired columns have no more mismatches is an alignment of equal score.
        let (lo, (qs, ts)) = (a, p[a]);
        let mm = |x: u8, y: u8| x != y || !base_ok(x) || !base_ok(y);
        let mut run = 0usize;
        while run < lo && ops[lo - 1 - run] == b'M' {
            run += 1;
        }
        let run = run.min(qs).min(ts - toff);
        let (mut delta, mut best) = (0i64, 0usize);
        for i in 1..=run {
            let (old, new) = if op == b'D' {
                (mm(q[qs - i], t[ts - i]), mm(q[qs - i], t[ts + len - i]))
            } else {
                (mm(q[qs - i], t[ts - i]), mm(q[qs + len - i], t[ts - i]))
            };
            delta += new as i64 - old as i64;
            if delta <= 0 {
                best = i;
            }
            if delta > 4 || (i > 256 && delta > 0) {
                break;
            }
        }
        if best > 0 {
            // Rotate the `best` columns from before the gap to after it.
            ops[lo - best..lo + len].rotate_left(best);
        }
        a = b;
    }
    let mut out = String::new();
    let mut i = 0;
    while i < ops.len() {
        let mut j = i;
        while j < ops.len() && ops[j] == ops[i] {
            j += 1;
        }
        out.push_str(&format!("{}{}", j - i, ops[i] as char));
        i = j;
    }
    out
}

/// Write the read's base qualities into the QUAL column of its records.
/// They are reversed on reverse-strand records, whose SEQ is reverse-complemented, and FASTA reads keep '*'.
/// Supplementary records are hard-clipped, holding only the aligned bases in SEQ as minimap2 and BWA-MEM write them, so their QUAL is cut by the CIGAR's leading and trailing H the same way.
fn fill_qual(lines: &mut [String], qual: &[u8]) {
    if qual.is_empty() {
        return;
    }
    for line in lines.iter_mut() {
        let mut f: Vec<&str> = line.splitn(12, '\t').collect();
        if f.len() < 11 || f[9] == "*" {
            continue;
        }
        let flag: u16 = f[1].parse().unwrap_or(0);
        let oriented: Vec<u8> = if flag & SAM_FLAG_REVERSE != 0 {
            qual.iter().rev().copied().collect()
        } else {
            qual.to_vec()
        };
        let (lead, trail) = hard_clips(f[5]);
        if lead + trail > oriented.len() {
            continue;
        }
        let q: String = oriented[lead..oriented.len() - trail]
            .iter()
            .map(|&b| b as char)
            .collect();
        if q.len() != f[9].len() {
            continue;
        }
        f[10] = &q;
        let joined = f.join("\t");
        *line = joined;
    }
}

/// Deterministic per-read hash, 64-bit FNV-1a followed by MurmurHash3's fmix64 steps, that seeds the choice among equally good placements.
/// The same read thus gets the same copy at any thread count.
fn read_hash(seq: &[u8]) -> u64 {
    const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut h: u64 = FNV_OFFSET_BASIS;
    for &b in seq {
        h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h
}

/// Soft-clip the ends of an end-to-end alignment, an M/I/D CIGAR of `q` against `t` from `toff`, under the EXT_* scoring.
/// Each end is cut at its lowest-scoring prefix or suffix if that score is below zero, and only just before a matching column, so the kept core starts and ends with a match.
/// Returns the CIGAR with S ops, the new reference offset and the core's NM, or None when nothing is clipped or the core would keep under half the read.
fn clip_ends(cig: &str, q: &[u8], t: &[u8], toff: usize) -> Option<(String, usize, i32)> {
    // Expand into columns of op and is_match.
    let mut cols: Vec<(u8, bool)> = Vec::with_capacity(q.len() + 16);
    let (mut qi, mut ti, mut n) = (0usize, toff, 0usize);
    for ch in cig.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        for _ in 0..n {
            match ch {
                b'M' => {
                    if qi >= q.len() || ti >= t.len() {
                        return None;
                    }
                    let m = q[qi] == t[ti] && matches!(t[ti], b'A' | b'C' | b'G' | b'T');
                    cols.push((b'M', m));
                    qi += 1;
                    ti += 1;
                }
                b'I' => {
                    if qi >= q.len() {
                        return None;
                    }
                    cols.push((b'I', false));
                    qi += 1;
                }
                b'D' => {
                    if ti >= t.len() {
                        return None;
                    }
                    cols.push((b'D', false));
                    ti += 1;
                }
                _ => return None,
            }
        }
        n = 0;
    }
    if qi != q.len() {
        return None;
    }
    let score = |c: (u8, bool), prev: Option<u8>| -> i64 {
        match c {
            (b'M', true) => EXT_MATCH,
            (b'M', false) => EXT_MISMATCH,
            (op, _) => {
                if prev == Some(op) {
                    EXT_GAP_EXT
                } else {
                    EXT_GAP_OPEN + EXT_GAP_EXT
                }
            }
        }
    };
    let nc = cols.len();
    // left counts the columns to clip.
    let (mut cum, mut best, mut lcut) = (0i64, 0i64, 0usize);
    for c in 0..nc {
        cum += score(cols[c], if c > 0 { Some(cols[c - 1].0) } else { None });
        if cum < best && c + 1 < nc && cols[c + 1] == (b'M', true) {
            best = cum;
            lcut = c + 1;
        }
    }
    // right
    let (mut cum, mut best, mut rcut) = (0i64, 0i64, 0usize);
    for c in (0..nc).rev() {
        cum += score(
            cols[c],
            if c + 1 < nc {
                Some(cols[c + 1].0)
            } else {
                None
            },
        );
        if cum < best && c >= 1 && cols[c - 1] == (b'M', true) {
            best = cum;
            rcut = nc - c;
        }
    }
    if lcut == 0 && rcut == 0 {
        return None;
    }
    if lcut + rcut >= nc {
        return None;
    }
    let qlen = |cs: &[(u8, bool)]| cs.iter().filter(|c| c.0 != b'D').count();
    let tlen = |cs: &[(u8, bool)]| cs.iter().filter(|c| c.0 != b'I').count();
    let (left, core, right) = (&cols[..lcut], &cols[lcut..nc - rcut], &cols[nc - rcut..]);
    let (ql, qr, qc) = (qlen(left), qlen(right), qlen(core));
    if qc * 2 < q.len() {
        return None;
    }
    let mut out = String::new();
    if ql > 0 {
        out.push_str(&format!("{ql}S"));
    }
    let mut i = 0;
    let mut nm = 0i32;
    while i < core.len() {
        let op = core[i].0;
        let mut j = i;
        while j < core.len() && core[j].0 == op {
            if op != b'M' || !core[j].1 {
                nm += 1;
            }
            j += 1;
        }
        out.push_str(&format!("{}{}", j - i, op as char));
        i = j;
    }
    if qr > 0 {
        out.push_str(&format!("{qr}S"));
    }
    Some((out, toff + tlen(left), nm))
}

/// Clipped bases that count as one more alignment event in [`cigar_events`].
const CLIP_BASES_PER_EVENT: usize = 12;

/// Alignment events of a CIGAR (mismatches, gap openings, and clipped ends counted as one plus one per CLIP_BASES_PER_EVENT clipped bases), the error measure of the MAPQ cap.
/// A 10 bp indel is one event, not ten, because counting its length made reads carrying a true indel look like high-error reads and capped them near MAPQ 30, where pipelines filter.
fn cigar_events(cig: &str, q: &[u8], t: &[u8], toff: usize) -> i32 {
    let (mut qi, mut ti, mut n, mut ev) = (0usize, toff, 0usize, 0i32);
    for ch in cig.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        match ch {
            b'M' => {
                for _ in 0..n {
                    if qi < q.len()
                        && ti < t.len()
                        && (q[qi] != t[ti] || !matches!(t[ti], b'A' | b'C' | b'G' | b'T'))
                    {
                        ev += 1;
                    }
                    qi += 1;
                    ti += 1;
                }
            }
            b'I' => {
                ev += 1;
                qi += n;
            }
            b'D' => {
                ev += 1;
                ti += n;
            }
            // A clipped end is one event plus one per 12 clipped bases, since a long clip is read sequence the placement does not explain.
            // An aligned part fitting a diverged repeat copy, while the whole read's locus was never seeded, was capped at MAPQ 41 or so.
            b'S' => {
                ev += 1 + (n / CLIP_BASES_PER_EVENT) as i32;
                qi += n;
            }
            _ => {}
        }
        n = 0;
    }
    ev
}

/// With OPTIK_MAPQ_DEBUG set, short-read records carry an XD:Z tag of the MAPQ inputs (ed, events, loci, minimizers, repetitive minimizers, best-locus votes, other-locus votes, competitors with their eds).
fn mapq_debug() -> bool {
    static D: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *D.get_or_init(|| std::env::var_os("OPTIK_MAPQ_DEBUG").is_some())
}

/// MAPQ of a short-read hit, the posterior over its verified competitors and the odds that its true locus was never a candidate (hidden_odds), capped by its seed evidence (seed_cap).
fn final_mapq(hit: &ShortHit, events: i32, len: usize, rho: f64) -> i32 {
    let w: Vec<(i32, f64)> = hit.comp.iter().map(|&e| (e, 1.0)).collect();
    posterior_mapq_w(hit.ed, &w, len, hidden_odds(hit, events, len, rho)).min(seed_cap(hit))
}

/// Immutable lookup context shared by every read's alignment attempt.
/// `seeds` is the minimizer seed index of include/seed_index.hpp, whose hits carry the strand on which the read matches the reference.
struct AlignCtx<'a> {
    k: usize,
    seeds: &'a SeedIdx,
    ref_seqs: &'a [Vec<u8>],
    /// cs tag form, with MD always written
    cs: CsMode,
    /// the reference's minimizer repeat fraction from SeedIdx::repeat_fraction, a MAPQ input of hidden_odds
    rho: f64,
}

/// minimap2's cs tag, off by default, short with `:n` runs, or long with `=SEQ`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum CsMode {
    Off,
    Short,
    Long,
}

/// MD, per the SAM spec, and optionally minimap2's cs of one record.
/// The CIGAR's clips are skipped, and aligned columns are compared as NM counts them, by equal bases, where a non-ACGT reference base never matches.
/// Deletions are spelled from the reference and, in cs only, insertions from the read.
/// MD keeps a number between any two elements, as in `10A0^TT5`, the way samtools calmd writes it.
fn diff_tags(cig: &str, seq: &[u8], r: &[u8], pos: usize, cs: CsMode) -> (String, Option<String>) {
    let lc = |b: u8| (b as char).to_ascii_lowercase();
    let mut md = String::new();
    let mut run = 0usize; // MD matching run
    let mut csv = String::new();
    let mut ident = 0usize; // cs short identical run
    let mut ident_long = String::new(); // cs long identical bases
    let flush_cs = |csv: &mut String, ident: &mut usize, il: &mut String| {
        if *ident > 0 {
            csv.push_str(&format!(":{ident}"));
            *ident = 0;
        }
        if !il.is_empty() {
            csv.push('=');
            csv.push_str(il);
            il.clear();
        }
    };
    let (mut qi, mut ti, mut n) = (0usize, pos, 0usize);
    for ch in cig.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as usize;
            continue;
        }
        match ch {
            b'S' => qi += n,
            b'H' => {}
            b'M' | b'=' | b'X' => {
                for _ in 0..n {
                    let (qb, rb) = (
                        seq.get(qi).copied().unwrap_or(b'N'),
                        r.get(ti).copied().unwrap_or(b'N'),
                    );
                    if qb == rb && matches!(rb, b'A' | b'C' | b'G' | b'T') {
                        run += 1;
                        match cs {
                            CsMode::Short => ident += 1,
                            CsMode::Long => ident_long.push(qb as char),
                            CsMode::Off => {}
                        }
                    } else {
                        md.push_str(&run.to_string());
                        md.push(rb as char);
                        run = 0;
                        if cs != CsMode::Off {
                            flush_cs(&mut csv, &mut ident, &mut ident_long);
                            csv.push('*');
                            csv.push(lc(rb));
                            csv.push(lc(qb));
                        }
                    }
                    qi += 1;
                    ti += 1;
                }
            }
            b'I' => {
                if cs != CsMode::Off {
                    flush_cs(&mut csv, &mut ident, &mut ident_long);
                    csv.push('+');
                    let end = (qi + n).min(seq.len());
                    if qi < end {
                        for &b in &seq[qi..end] {
                            csv.push(lc(b));
                        }
                    }
                }
                qi += n;
            }
            b'D' | b'N' => {
                md.push_str(&run.to_string());
                md.push('^');
                run = 0;
                if cs != CsMode::Off {
                    flush_cs(&mut csv, &mut ident, &mut ident_long);
                    csv.push('-');
                }
                for j in ti..ti + n {
                    let rb = r.get(j).copied().unwrap_or(b'N');
                    md.push(rb as char);
                    if cs != CsMode::Off {
                        csv.push(lc(rb));
                    }
                }
                ti += n;
            }
            _ => {}
        }
        n = 0;
    }
    md.push_str(&run.to_string());
    if cs == CsMode::Off {
        return (md, None);
    }
    flush_cs(&mut csv, &mut ident, &mut ident_long);
    (md, Some(csv))
}

/// "\tMD:Z:..[\tcs:Z:..]" for a record.
fn diff_tag_fields(cig: &str, seq: &[u8], r: &[u8], pos: usize, cs: CsMode) -> String {
    let (md, c) = diff_tags(cig, seq, r, pos, cs);
    match c {
        Some(c) => format!("\tMD:Z:{md}\tcs:Z:{c}"),
        None => format!("\tMD:Z:{md}"),
    }
}

/// Verified placement of one read.
struct ShortHit {
    ed: i32,
    /// distinct competitors' edit distances and the fraction of the read's minimizers dropped as over-occurring, both MAPQ inputs of final_mapq
    comp: Vec<i32>,
    /// Where each competitor lies, as (rid, rev, exclusive reference end).
    /// Pairing discounts the ones its mate does not support.
    comp_loc: Vec<(u32, bool, i64)>,
    rep_frac: f64,
    /// Seed evidence for hidden_odds, namely candidate loci before the MAX_SHORT_CANDIDATES cut, the read's minimizers, the over-occurring ones, and the seed votes of the best locus and of the strongest other locus.
    /// All are zero for a mate placed by rescue, which has no seeds of its own.
    seed_ev: [u32; 5],
    /// Edits from the best placement to the nearest verified competitor, or MAX_GAP when none is that close, used by hidden_odds.
    /// It is 0 for a placement that is not the read's own best, such as a pair-chosen candidate or a rescued mate, whose own best elsewhere competes with it.
    gap: i32,
    rid: usize,
    rev: bool,
    win_start: usize,
    win_end: usize,
    /// exclusive reference end of the best placement
    end: i64,
}

/// Cap of ShortHit::gap, beyond which competitors over 12 edits away are not counted.
const MAX_GAP: i32 = 13;
/// Cap on candidate loci verified per read, taken in order of seed support.
const MAX_SHORT_CANDIDATES: usize = 128;
/// Bases of padding around a candidate window beyond the edit budget.
const WINDOW_PAD: usize = 8;

/// Placements whose reference ends differ by at most this, a tenth of the read and at least 10 bases, are the same placement.
/// The truth scoring uses the same read-origin tolerance.
fn placement_sep(len: usize) -> usize {
    (len / 10).max(10)
}

/// Finds and verifies a short read's best placement and its competitors, or returns None if no candidate locus aligns within `max_ed`.
fn align_read_impl(seq: &[u8], max_ed: i32, ctx: &AlignCtx<'_>) -> Option<ShortHit> {
    let l = seq.len();
    if l < ctx.k {
        return None;
    }
    let pad = (max_ed.max(0) as usize) + WINDOW_PAD;
    // rid, rev, diagonal start
    let (seed_hits, n_min, n_rep) = ctx.seeds.collect_with_rep(seq);
    let mut votes: Vec<(u32, bool, i64)> = seed_hits
        .iter()
        .map(|h| (h.rid, h.rev != 0, h.rpos as i64 - h.qpos as i64))
        .collect();
    if votes.is_empty() {
        return None;
    }
    votes.sort_unstable();
    // Merge diagonals within `pad` of the cluster's first start.
    struct Locus {
        rid: u32,
        rev: bool,
        lo: i64,
        hi: i64,
        n: usize,
    }
    let mut loci: Vec<Locus> = Vec::new();
    for &(rid, rev, st) in &votes {
        match loci.last_mut() {
            Some(c) if c.rid == rid && c.rev == rev && st - c.lo <= pad as i64 => {
                c.hi = st;
                c.n += 1;
            }
            _ => loci.push(Locus {
                rid,
                rev,
                lo: st,
                hi: st,
                n: 1,
            }),
        }
    }
    loci.sort_by(|a, b| {
        b.n.cmp(&a.n)
            .then(a.rid.cmp(&b.rid))
            .then(a.rev.cmp(&b.rev))
            .then(a.lo.cmp(&b.lo))
    });
    let n_loci_all = loci.len() as u32;
    loci.truncate(MAX_SHORT_CANDIDATES);

    let rc = revcomp(seq);
    // Verify every candidate locus up to MAX_SHORT_CANDIDATES, in order of seed support.
    // Seed support is not alignment quality.
    // At a few percent error the true copy of a repeat can keep fewer error-free minimizers than a diverged copy, so skipping weakly supported loci once a strong one aligns would pick the wrong copy at full confidence.
    // The loci share the read, so they are verified four at a time, bit-parallel with one SIMD lane per locus.
    let wins: Vec<(usize, usize)> = loci
        .iter()
        .map(|c| {
            let rlen = ctx.ref_seqs[c.rid as usize].len();
            let ws = (c.lo - pad as i64).max(0) as usize;
            let we = ((c.hi + l as i64 + pad as i64).max(0) as usize).min(rlen);
            (ws, we)
        })
        .collect();
    // Two placements are distinct when their reference ends differ by more than the read-origin tolerance the truth scoring also uses.
    let sep = placement_sep(l) as i32;
    let mut pls = vec![
        crate::Placements {
            best: max_ed + 1,
            end: 0,
            alt: max_ed + 1,
            alt_end: 0
        };
        loci.len()
    ];
    for rev in [false, true] {
        let q: &[u8] = if rev { &rc } else { seq };
        let idx: Vec<usize> = (0..loci.len())
            .filter(|&i| loci[i].rev == rev && wins[i].0 < wins[i].1)
            .collect();
        for grp in idx.chunks(crate::PLACEMENT_LANES) {
            let texts: Vec<&[u8]> = grp
                .iter()
                .map(|&i| &ctx.ref_seqs[loci[i].rid as usize][wins[i].0..wins[i].1])
                .collect();
            let out = crate::placements_x4(q, &texts, max_ed, sep);
            for (g, &i) in grp.iter().enumerate() {
                pls[i] = out[g];
            }
        }
    }
    // every verified locus, in locus order
    let found: Vec<usize> = (0..loci.len()).filter(|&i| pls[i].best <= max_ed).collect();
    // The best placement has the lowest edit distance.
    // Equally good distinct placements, i.e. repeat copies, are chosen between pseudo-randomly, seeded by the read's sequence as BWA and minimap2 do, so a repeat family's reads spread over its copies instead of all landing on the first verified one.
    // Taking the first left 0.6% of the Salmonella LT2 genome under 10x at 30x.
    // The choice is deterministic at any thread count, and the MAPQ of a tie is 0 either way.
    let best_ed = found.iter().map(|&i| pls[i].best).min()?;
    let mut tied: Vec<usize> = Vec::new();
    for &i in &found {
        if pls[i].best != best_ed {
            continue;
        }
        let end = wins[i].0 as i64 + pls[i].end as i64;
        let dup = tied.iter().any(|&j| {
            loci[j].rid == loci[i].rid
                && loci[j].rev == loci[i].rev
                && (wins[j].0 as i64 + pls[j].end as i64 - end).abs() <= sep as i64
        });
        if !dup {
            tied.push(i);
        }
    }
    tied.sort_by_key(|&i| {
        (
            loci[i].rid,
            wins[i].0 as i64 + pls[i].end as i64,
            loci[i].rev,
        )
    });
    let bi = if tied.len() == 1 {
        tied[0]
    } else {
        tied[(read_hash(seq) % tied.len() as u64) as usize]
    };
    let (ed, rid, rev, (ws, we)) = (pls[bi].best, loci[bi].rid, loci[bi].rev, wins[bi]);
    let best_end = ws as i64 + pls[bi].end as i64;
    // Competitors are distinct placements.
    // Another window's best alignment counts unless it ends where the best one does, which is the same alignment seen from an overlapping window.
    // Every window's best alternative end counts too, because a read inside a tandem repeat aligns equally well one period over, often inside that same window.
    // Placements reported by several overlapping windows are merged by reference end.
    let mut alts: Vec<(u32, bool, i64, i32)> = Vec::new(); // (rid, rev, end, ed)
    for &i in &found {
        let p = pls[i];
        let base = wins[i].0 as i64;
        let distinct =
            |e: i64| loci[i].rid != rid || loci[i].rev != rev || (e - best_end).abs() > sep as i64;
        if i != bi && distinct(base + p.end as i64) {
            alts.push((loci[i].rid, loci[i].rev, base + p.end as i64, p.best));
        }
        if p.alt <= max_ed && distinct(base + p.alt_end as i64) {
            alts.push((loci[i].rid, loci[i].rev, base + p.alt_end as i64, p.alt));
        }
    }
    alts.sort_unstable();
    let mut comp: Vec<i32> = Vec::new();
    let mut comp_loc: Vec<(u32, bool, i64)> = Vec::new();
    let mut last: Option<(u32, bool, i64)> = None;
    for &(r, v, e, d) in &alts {
        match last {
            Some((lr, lv, le)) if lr == r && lv == v && e - le <= sep as i64 => {
                let x = comp.last_mut().unwrap();
                *x = (*x).min(d);
            }
            _ => {
                comp.push(d);
                comp_loc.push((r, v, e));
                last = Some((r, v, e));
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
    let n_best = loci[bi].n as u32;
    let n_other = loci
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != bi)
        .map(|(_, c)| c.n as u32)
        .max()
        .unwrap_or(0);
    let gap = comp
        .iter()
        .map(|&d| d - ed)
        .min()
        .unwrap_or(MAX_GAP)
        .clamp(0, MAX_GAP);
    Some(ShortHit {
        ed,
        comp,
        comp_loc,
        rep_frac,
        seed_ev: [n_loci_all, n_min, n_rep, n_best, n_other],
        gap,
        rid: rid as usize,
        rev,
        win_start: ws,
        win_end: we,
        end: best_end,
    })
}

/// MAPQ = -10 log10 P(best locus wrong), with P = (S + H) / (1 + S + H), where
///   - S is the sum over verified distinct competitors of r^(their edits - best), with r = edit_ratio, the likelihood ratio of one more edit, and
///   - H is hidden_odds, the odds that the true locus was never verified.
///
/// A competitor as good as the best is a tie with MAPQ 0.
///
/// The constants were fitted by maximum likelihood on simulated reads of 14 genomes, namely E. coli, Salmonella, four synthetic repeat genomes, hg38 chr1, chr2, chr7, chr17, chr19, chr21, chrX and the whole assembly, at 100 / 150 / 250 bp and 0.5-8% error.
/// They are judged on fresh simulated reads with bench/mapq_calibration.py --multi, on simulated pairs with paired_bench.py, and on real reads with bench/real_reads.py, using pair-anchored HG002 and E. coli mates and GIAB HG002 chr20 variant calling.
/// On none of them does a MAPQ bin above 0 have significantly more wrong reads than its MAPQs imply.
/// The real reads showed what simulation could not, namely diverged Alu and L1 copies, which hidden_odds' competitor gap handles, satellites with few seeded minimizers, handled by its seeded-share weight, and the sample's own variants shared with a paralog, which SAMPLE_DIVERGENCE covers.
///
/// posterior_mapq_w takes a prior weight per competitor, which pairing uses to weigh a placement its mate does not support by discordant_prior.
/// A competitor at full weight as good as the best is a tie at 0, while a discounted one may even beat the best, which is then the mate-supported placement.
fn posterior_mapq_w(best_ed: i32, competitors: &[(i32, f64)], len: usize, hidden: f64) -> i32 {
    if competitors.iter().any(|&(e, w)| e <= best_ed && w >= 1.0) {
        return 0;
    }
    let r = edit_ratio(best_ed, len);
    let s: f64 = competitors
        .iter()
        .map(|&(e, w)| w * r.powi(e - best_ed))
        .sum::<f64>()
        + hidden;
    let p_wrong = (s / (1.0 + s)).max(1e-9);
    ((-10.0 * p_wrong.log10()).round() as i32).clamp(0, SAM_MAX_MAPQ)
}

/// Odds that the read's true locus was never a verified candidate, because its seeds were destroyed by read errors, filtered as over-occurring, or cut with the candidates beyond MAX_SHORT_CANDIDATES while another copy kept some.
/// Each seed vote of the best locus divides the odds by about 1.7, and a candidate list cut at MAX_SHORT_CANDIDATES multiplies them by about 1.6, since copies beyond the cut were never aligned.
/// Over-occurring minimizers multiply them by e^(2.8 + 5.3 x their fraction), and the best alignment's event rate f = events / len, counting mismatches, gap openings and clipped ends, by e^(71 f), since at high divergence the true copy is the one most likely to lose its seeds.
/// Each edit between the best placement and its nearest verified competitor divides them by about 5, or about 12 in a repeat family, up to about 4.2 edits.
/// Having no verified competitor within 12 edits divides them by e^8.6, or e^9.7 in a repeat family.
/// Both terms are weighted by the share of seeded minimizers and by the genome's repeat fraction rho, as described below.
/// Without the competitor term, reads of diverged Alu and L1 copies with 0-1 edits got MAPQ 7-14 on real HG002 reads, where BWA-MEM gives 60 and GIAB confirms the placement.
/// Fitted by maximum likelihood on simulated reads of 14 genomes and checked on real reads with bench/real_reads.py.
fn hidden_odds(h: &ShortHit, events: i32, len: usize, rho: f64) -> f64 {
    let [loci, n_min, n_rep, n_best, _] = h.seed_ev;
    let f = events.max(0) as f64 / len.max(1) as f64;
    let nb = n_best.min(25) as f64;
    // The nearest verified competitor matters because a copy family that is diverged here, with every verified copy several edits away or none within 12, is unlikely to hold an unseen copy that fits the read as well.
    // That is evidence only as far as closer copies would have been seeded, so it is weighted by the share of the read's minimizers that seeded, (1 - over-occurring fraction)^1.04, and by 1 - 1.54 rho, since in a repeat-dense genome more copies are never seeded.
    // A rescued mate, with no seed vote of its own, carries none.
    // Without the weights, reads of chr7's centromeric satellite with a few seeded minimizers and verified copies 10 or more edits away got MAPQ 10-17 and were wrong.
    let rf = n_rep as f64 / n_min.max(1) as f64;
    let w = (1.0 - 1.5407 * rho).max(0.0) * (1.0 - rf).powf(1.0366) * (1.0 - (-8.1842 * nb).exp());
    let none = h.gap >= MAX_GAP;
    let gap = if none {
        0.0
    } else {
        (h.gap as f64).min(4.2475)
    };
    let mut l = -2.8093 - 0.5383 * nb
        + if loci > MAX_SHORT_CANDIDATES as u32 {
            0.485
        } else {
            0.0
        }
        + 71.0931 * f
        - 0.0527 * n_min.min(60) as f64
        - 9.0619 * rho;
    if n_rep > 0 {
        // Over-occurring minimizers mean the read lies (partly) in a family of 1000 or more copies, of which only the seeded ones were verified.
        l += 2.8365
            + 5.2724 * rf
            + 0.1715 * nb.min(10.0)
            + 0.0362 * (nb - 10.0).max(0.0)
            + w * (-2.5081 * gap + if none { -9.7409 } else { 0.0 });
    } else {
        l += w * (-1.6371 * gap + if none { -8.5511 } else { 0.0 });
    }
    l.min(30.0).exp()
}

/// Cap from seed evidence alone, SEED_CAP_BASE + SEED_CAP_PER_VOTE per seed vote of the best locus.
/// A best locus with few seed votes may be a copy that only a read error made seedable, so one vote allows at most MAPQ 25.
/// Reads in repeat families with over-occurring minimizers get no separate cap, since `hidden_odds` models them, fitted on 14 simulated genomes and checked on real HG002 reads.
/// A hand-set cap held their MAPQ near 13 and cost variant recall.
/// A rescued mate without seeds of its own is not capped here, because its partner placed it.
fn seed_cap(h: &ShortHit) -> i32 {
    let [_, n_min, _, n_best, _] = h.seed_ev;
    if n_min == 0 {
        return SAM_MAX_MAPQ;
    }
    const SEED_CAP_BASE: f64 = 20.0;
    const SEED_CAP_PER_VOTE: f64 = 5.0;
    let c = SEED_CAP_BASE + SEED_CAP_PER_VOTE * n_best as f64;
    (c.floor() as i32).clamp(0, SAM_MAX_MAPQ)
}

/// SAM QNAME length limit.
const MAX_QNAME: usize = 254;

/// Parses a `-R` read-group header line (the bwa and minimap2 form, with `\t` written as the two characters or as a tab), returning the normalized line and its ID.
pub(crate) fn parse_read_group(rg: &str) -> Result<(String, String), String> {
    let line = rg.replace("\\t", "\t");
    if !line.starts_with("@RG\t") {
        return Err(format!(
            "-R must be a read-group line starting with '@RG\\t', got '{rg}'"
        ));
    }
    if line.contains('\n') || line.contains('\r') {
        return Err("-R may not contain a line break".into());
    }
    let id = line
        .split('\t')
        .skip(1)
        .find_map(|f| f.strip_prefix("ID:"))
        .filter(|id| !id.is_empty())
        .ok_or_else(|| format!("-R read group has no ID field: '{rg}'"))?;
    for f in line.split('\t').skip(1) {
        if f.len() < 4 || f.as_bytes()[2] != b':' {
            return Err(format!("-R field '{f}' is not TAG:VALUE"));
        }
    }
    Ok((line.clone(), id.to_string()))
}

/// Seed occurrence cap of the mapper's index, SeedIndex max_occ.
/// Minimizers occurring more often are dropped from short-read seeding.
const SEED_MAX_OCC: u32 = 1000;
/// Short reads over their edit budget are retried on the chaining path if they are at least this many k-mers long.
const RESCUE_MIN_KMERS: usize = 4;

/// Runs `optik map`, mapping single-end reads, or pairs when `reads2` is given, to the reference at `ref_path` and writing SAM to `output`, or stdout if None.
///
/// `k_auto` chooses k from the first batch, 15 for noisy long reads and `k` otherwise.
/// `w` 0 tunes the window from the read length.
/// `max_ed` < 0 uses the automatic edit budget.
/// `cmdline` goes into the @PG header line.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run(
    ref_path: &str,
    reads_path: &str,
    reads2: Option<&str>,
    k: usize,
    k_auto: bool,
    w: usize,
    max_ed: i32,
    output: Option<&str>,
    threads: usize,
    cs: CsMode,
    read_group: Option<&str>,
    cmdline: &str,
) -> Result<(), String> {
    // -R adds a header line plus an RG:Z tag on every record, as bwa and minimap2 -R do, since GATK and Picard require read groups.
    let rg = read_group.map(parse_read_group).transpose()?;
    if k_auto {
        eprintln!(
            "[map] loading reference {ref_path} (k chosen from the reads, default {k}; w={w})..."
        );
    } else {
        eprintln!("[map] loading reference {ref_path} (k={k} w={w})...");
    }
    let t0 = Instant::now();

    // Load the reference genome into memory as a list of named sequences.
    let mut ref_names: Vec<String> = Vec::new();
    let mut ref_seqs: Vec<Vec<u8>> = Vec::new();
    iterate_fasta(ref_path, |name, seq| {
        ref_names.push(name.to_string());
        // Normalized reference.
        // Soft-masked lowercase bases count as the same bases, since case-sensitive alignment scored every one as a mismatch, RNA U becomes T, and any other non-ACGT code becomes N.
        ref_seqs.push(
            seq.iter()
                .map(|&b| match b.to_ascii_uppercase() {
                    c @ (b'A' | b'C' | b'G' | b'T') => c,
                    b'U' => b'T',
                    _ => b'N',
                })
                .collect(),
        );
    })?;
    // SAM requires unique @SQ SN names, since a duplicate would make every record on it ambiguous, so the run refuses instead of writing an invalid header.
    {
        let mut seen = std::collections::HashSet::with_capacity(ref_names.len());
        for name in &ref_names {
            if !seen.insert(name.as_str()) {
                return Err(format!(
                    "duplicate reference sequence name '{name}' in {ref_path} \
                                    (SAM @SQ names must be unique)"
                ));
            }
        }
    }

    // Open output and write SAM header.
    let mut out: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(
            std::fs::File::create(path).map_err(|e| e.to_string())?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted").map_err(|e| e.to_string())?;
    for (i, name) in ref_names.iter().enumerate() {
        writeln!(out, "@SQ\tSN:{name}\tLN:{}", ref_seqs[i].len()).map_err(|e| e.to_string())?;
    }
    if let Some((line, _)) = &rg {
        writeln!(out, "{line}").map_err(|e| e.to_string())?;
    }
    // CL holds the command line, so the file records how it was made (tabs and line breaks would break the header line).
    let cl: String = cmdline
        .chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect();
    writeln!(
        out,
        "@PG\tID:optik\tPN:optik\tVN:{}\tCL:{cl}",
        env!("CARGO_PKG_VERSION")
    )
    .map_err(|e| e.to_string())?;
    let rg_tag: String = rg
        .as_ref()
        .map(|(_, id)| format!("\tRG:Z:{id}"))
        .unwrap_or_default();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| format!("failed to build {threads}-thread pool: {e}"))?;

    // Reads are streamed in bounded batches.
    // Each batch is mapped in parallel and written in input order, which is deterministic at any thread count and bounds memory by one batch of reads and records.
    // The index is built when the first batch arrives, so `-w 0` can be tuned from real read lengths.
    const BATCH: usize = 1 << 16;
    // name, seq, qual, and R2's seq and qual in paired mode
    type Rd = MapRead;
    let mut batch: Vec<Rd> = Vec::with_capacity(BATCH);
    let mut model: Option<InsertModel> = None;
    let mut state: Option<(usize, usize, SeedIdx)> = None; // (w, k, index)
    let (mut n_reads, mut n_mapped) = (0u64, 0u64);
    let mut io_err: Option<String> = None;
    let tmap = Instant::now();

    let mut flush = |batch: &mut Vec<Rd>,
                     state: &mut Option<(usize, usize, SeedIdx)>|
     -> Result<(), String> {
        if batch.is_empty() {
            return Ok(());
        }
        if state.is_none() {
            let w = if w == 0 {
                let avg = batch.iter().map(|(_, s, _, _)| s.len()).sum::<usize>() / batch.len();
                let tuned = auto_tune_w(avg, k);
                eprintln!("[map] auto-tuned w = {tuned}");
                tuned
            } else {
                w
            };
            // Minimizer MPHF seed index that uses the reference as its key store
            for (rid, seq) in ref_seqs.iter().enumerate() {
                if seq.len() >= (1usize << 31) {
                    return Err(format!("reference {} exceeds 2^31 bases", ref_names[rid]));
                }
            }
            let mut k = k;
            let mut why = String::new();
            if k_auto {
                // Noisy long reads seed better with shorter minimizers.
                // On pbsim2 R94 ONT reads of CHM13, k = 15 placed 149,895 reads at MAPQ >= 10 with 0 wrong and k = 21 placed 147,815 with 3 wrong.
                // HiFi-like reads lose nothing at k = 21, which runs 2.3x faster.
                match noisy_long_reads(batch) {
                    ReadNoise::Noisy(e) => {
                        k = NOISY_K;
                        why = format!("base qualities: {:.1}% error", 100.0 * e);
                    }
                    ReadNoise::Accurate(e) => {
                        why = format!("base qualities: {:.2}% error", 100.0 * e)
                    }
                    ReadNoise::Short => why = "short reads".into(),
                    ReadNoise::Unknown => {}
                }
            }
            let mut seeds = SeedIdx::build(&ref_seqs, k, w, SEED_MAX_OCC, threads.max(1))?;
            if k_auto && why.is_empty() {
                // Long reads without qualities are judged by the share of their minimizers the reference holds, about (1 - error)^k, which at k = 21 is 96% at 0.2% error and 34% at 5%.
                let f = seed_presence(batch, &seeds, k);
                if f < NOISY_PRESENCE {
                    k = NOISY_K;
                    why = format!("{:.0}% of long-read minimizers found at k = 21", 100.0 * f);
                    drop(seeds); // free the k = 21 index before building the k = 15 one
                    seeds = SeedIdx::build(&ref_seqs, k, w, SEED_MAX_OCC, threads.max(1))?;
                } else {
                    why = format!("{:.0}% of long-read minimizers found", 100.0 * f);
                }
            }
            if k_auto {
                eprintln!("[map] k = {k} ({why})");
            }
            eprintln!(
                "[map] {} minimizers indexed ({:.1} MB, repeat fraction {:.4}, chaining occurrence \
                       threshold {}) in {:.2}s",
                seeds.distinct(),
                seeds.bytes() as f64 / 1048576.0,
                seeds.repeat_fraction(),
                seeds.mid_occ(),
                t0.elapsed().as_secs_f64()
            );
            *state = Some((w, k, seeds));
        }
        let (_, k, seeds) = state.as_ref().unwrap();
        let k = *k;
        let ctx = AlignCtx {
            k,
            seeds,
            ref_seqs: &ref_seqs,
            cs,
            rho: seeds.repeat_fraction(),
        };
        let map_one = |read_name: &str, seq: &[u8]| -> Vec<String> {
            if seq.len() >= CHAIN_MIN_LEN {
                map_one_read_chain(read_name, seq, &ctx, &ref_names, LONG_CHAIN)
            } else {
                let rec = map_one_read_short(
                    read_name,
                    seq,
                    read_max_ed(max_ed, seq.len()),
                    &ctx,
                    &ref_names,
                );
                // Short reads beyond -e, i.e. high-error ones, get a second chance on the chaining path, which has no edit-distance cap.
                if is_unmapped(&rec) && seq.len() >= RESCUE_MIN_KMERS * ctx.k {
                    map_one_read_chain(read_name, seq, &ctx, &ref_names, RESCUE_CHAIN)
                } else {
                    vec![rec]
                }
            }
        };
        let recs: Vec<Vec<String>> = if reads2.is_none() {
            pool.install(|| {
                batch
                    .par_iter()
                    .map(|(name, seq, qual, _)| {
                        let mut lines = map_one(name, seq);
                        fill_qual(&mut lines, qual);
                        lines
                    })
                    .collect()
            })
        } else {
            // In paired mode both mates are placed on their own, then paired.
            let firsts: Vec<(MateAln, MateAln)> = pool.install(|| {
                batch
                    .par_iter()
                    .map(|(name, s1, _, m)| {
                        let s2 = &m.as_ref().unwrap().0;
                        (
                            map_mate(name, s1, max_ed, &ctx, &ref_names),
                            map_mate(name, s2, max_ed, &ctx, &ref_names),
                        )
                    })
                    .collect()
            });
            if model.is_none() {
                // template lengths of the first batch's confidently placed FR pairs (unpaired posterior >= 30 for both, no repetitive seeds)
                let sure = |h: &ShortHit, len: usize| {
                    h.rep_frac == 0.0 && final_mapq(h, h.ed, len, ctx.rho) >= 30
                };
                let tl: Vec<i64> = firsts
                    .iter()
                    .zip(batch.iter())
                    .filter_map(|((a, b), (_, s1, _, m))| {
                        let (MateAln::Short(Some(h1)), MateAln::Short(Some(h2))) = (a, b) else {
                            return None;
                        };
                        if !sure(h1, s1.len()) || !sure(h2, m.as_ref().unwrap().0.len()) {
                            return None;
                        }
                        let t = fr_tlen(
                            hit_place(h1, s1.len()),
                            hit_place(h2, m.as_ref().unwrap().0.len()),
                        )?;
                        (t > 0 && t < 100_000).then_some(t)
                    })
                    .collect();
                let n = tl.len();
                let est = insert_model(tl);
                let m = est.unwrap_or(InsertModel { lo: 1, hi: 1000 });
                if est.is_some() {
                    eprintln!(
                        "[map] paired: proper template length {}..{} (from {n} confident pairs)",
                        m.lo, m.hi
                    );
                } else {
                    eprintln!(
                        "[map] paired: only {n} confident FR pairs in the first batch; proper template \
                               length assumed {}..{}",
                        m.lo, m.hi
                    );
                }
                model = Some(m);
            }
            let m = model.unwrap();
            let pairs: Vec<(Vec<String>, Vec<String>)> = pool.install(|| {
                firsts
                    .into_par_iter()
                    .zip(batch.par_iter())
                    .map(|((a1, a2), (name, s1, q1, mt))| {
                        let (s2, q2) = mt.as_ref().unwrap();
                        let (mut l1, mut l2) =
                            pair_records(name, s1, s2, a1, a2, &m, max_ed, &ctx, &ref_names);
                        fill_qual(&mut l1, q1);
                        fill_qual(&mut l2, q2);
                        (l1, l2)
                    })
                    .collect()
            });
            pairs.into_iter().flat_map(|(a, b)| [a, b]).collect()
        };
        for lines in recs {
            n_reads += 1;
            // Count reads whose primary record, the first line, is mapped.
            let flag: u16 = lines
                .first()
                .and_then(|l| l.split('\t').nth(1))
                .and_then(|f| f.parse().ok())
                .unwrap_or(SAM_FLAG_UNMAPPED);
            if flag & SAM_FLAG_UNMAPPED == 0 {
                n_mapped += 1;
            }
            for line in lines {
                writeln!(out, "{line}{rg_tag}").map_err(|e| e.to_string())?;
            }
        }
        batch.clear();
        Ok(())
    };

    let mut qname_warned = false;
    // Every input record gets one SAM record, empty ones included.
    // Reads are aligned case-insensitively, with RNA U read as T.
    let mut qual_warned = false;
    let mut n_bad_qual = 0u64;
    let mut prep = |read_name: &str, seq: &[u8], qual: &[u8]| -> (String, Vec<u8>, Vec<u8>) {
        // SAM QUAL must be '*' or exactly as long as SEQ
        let qual: Vec<u8> = if qual.is_empty() {
            Vec::new()
        } else if qual.len() == seq.len() && qual.iter().all(|&q| (b'!'..=b'~').contains(&q)) {
            qual.to_vec()
        } else {
            // A wrong length, or a byte outside Phred+33 '!'..'~' that samtools rejects, gives QUAL '*' instead of a broken SAM.
            if !qual_warned {
                eprintln!(
                    "[map] warning: FASTQ record '{read_name}' has {} quality values for {} \
                               bases or a value outside '!'..'~'; its QUAL is written as '*'",
                    qual.len(),
                    seq.len()
                );
                qual_warned = true;
            }
            n_bad_qual += 1;
            Vec::new()
        };
        let seq: Vec<u8> = seq
            .iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'U' => b'T',
                c => c,
            })
            .collect();
        // SAM QNAME holds at most 254 characters, and samtools rejects longer ones.
        let name = if read_name.len() > MAX_QNAME {
            if !qname_warned {
                eprintln!(
                    "[map] warning: read names longer than {MAX_QNAME} characters \
                           are truncated (SAM QNAME limit)"
                );
                qname_warned = true;
            }
            let mut cut = MAX_QNAME;
            while !read_name.is_char_boundary(cut) {
                cut -= 1;
            }
            read_name[..cut].to_string()
        } else {
            read_name.to_string()
        };
        (name, seq, qual)
    };
    // In paired mode R2 is parsed on its own thread, but only while this thread waits for it.
    // This thread asks for the next R2_BLOCK records and blocks until they arrive, so at most one of the two reading threads runs at a time and the mapping pool keeps all -t threads (bench/thread_check.py).
    const R2_BLOCK: usize = 4096;
    type Rec = (String, Vec<u8>, Vec<u8>);
    type Msg = Result<Vec<Rec>, String>;
    let (want_tx, want_rx) = std::sync::mpsc::sync_channel::<()>(0);
    let (tx, rx) = std::sync::mpsc::sync_channel::<Msg>(0);
    let r2_thread = reads2.map(|p2| {
        let p2 = p2.to_string();
        std::thread::spawn(move || {
            let mut buf: Vec<Rec> = Vec::with_capacity(R2_BLOCK);
            let mut asked = want_rx.recv().is_ok();
            let mut stop = !asked;
            let r = iterate_reads(&p2, |n, s, q| {
                if stop {
                    return;
                }
                buf.push((n.to_string(), s.to_vec(), q.to_vec()));
                if buf.len() == R2_BLOCK {
                    let _ = tx.send(Ok(std::mem::replace(
                        &mut buf,
                        Vec::with_capacity(R2_BLOCK),
                    )));
                    asked = want_rx.recv().is_ok(); // wait for the next request or the end
                    stop = !asked;
                }
            });
            if !stop {
                // the last, partial block, then the error, or the end once the sender drops
                if !buf.is_empty() {
                    let _ = tx.send(Ok(buf));
                    let _ = want_rx.recv();
                }
                if let Err(e) = r {
                    let _ = tx.send(Err(e));
                }
            }
        })
    });
    // The R2 reader gets one block asked for at a time.
    // With -t >= 2 the next block is asked for as soon as one arrives, so R2 is parsed beside R1 by two reading threads while the mapping pool is idle.
    // An outstanding block is collected before every mapping batch, so the R2 reader is idle while the pool runs and the run never has more than -t threads busy.
    struct R2Feed {
        want: std::sync::mpsc::SyncSender<()>,
        rx: std::sync::mpsc::Receiver<Msg>,
        buf: std::collections::VecDeque<Rec>,
        pending: bool,
        end: bool,
        ahead: bool,
    }
    impl R2Feed {
        fn ask(&mut self) {
            if !self.pending && !self.end {
                if self.want.send(()).is_err() {
                    self.end = true;
                } else {
                    self.pending = true;
                }
            }
        }
        fn collect(&mut self) -> Option<String> {
            if !self.pending {
                return None;
            }
            self.pending = false;
            match self.rx.recv() {
                Ok(Ok(v)) => {
                    self.buf.extend(v);
                    None
                }
                Ok(Err(e)) => Some(e),
                Err(_) => {
                    self.end = true;
                    None
                }
            }
        }
        fn next(&mut self) -> Option<Result<Rec, String>> {
            if self.buf.is_empty() {
                self.ask();
                if let Some(e) = self.collect() {
                    return Some(Err(e));
                }
                if self.ahead {
                    self.ask();
                }
            }
            self.buf.pop_front().map(Ok)
        }
    }
    let mut r2 = R2Feed {
        want: want_tx,
        rx,
        buf: std::collections::VecDeque::new(),
        pending: false,
        end: false,
        ahead: threads > 1,
    };
    let mut n_rec = 0u64;
    iterate_reads(reads_path, |read_name, seq, qual| {
        if io_err.is_some() {
            return;
        }
        n_rec += 1;
        if reads2.is_some() {
            let (n2, s2, q2) = match r2.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => {
                    io_err = Some(e);
                    return;
                }
                None => {
                    io_err = Some(format!(
                        "the second reads file has fewer records than the first \
                                           (it ends before record {n_rec})"
                    ));
                    return;
                }
            };
            if mate_name(read_name.split_ascii_whitespace().next().unwrap_or(""))
                != mate_name(n2.split_ascii_whitespace().next().unwrap_or(""))
            {
                io_err = Some(format!(
                    "mate names differ at record {n_rec}: '{read_name}' vs '{n2}' \
                                       (R1 and R2 must list the pairs in the same order)"
                ));
                return;
            }
            let (name, seq, qual) = prep(mate_name(read_name), seq, qual);
            let (_, s2, q2) = prep(&n2, &s2, &q2);
            batch.push((name, seq, qual, Some((s2, q2))));
        } else {
            let (name, seq, qual) = prep(read_name, seq, qual);
            batch.push((name, seq, qual, None));
        }
        if batch.len() >= BATCH {
            // The R2 block being parsed is collected first, so no reader runs beside the pool.
            if reads2.is_some()
                && let Some(e) = r2.collect()
            {
                io_err = Some(e);
                return;
            }
            if let Err(e) = flush(&mut batch, &mut state) {
                io_err = Some(e);
            }
        }
    })?;
    if reads2.is_some() && io_err.is_none() {
        match r2.next() {
            Some(Ok((n2, _, _))) => {
                io_err = Some(format!(
                    "the second reads file has more records than the first (first extra: '{n2}')"
                ))
            }
            Some(Err(e)) => io_err = Some(e),
            None => {}
        }
    }
    let _ = r2.collect(); // a block still being parsed
    drop(r2); // the R2 reader stops at its next wait
    if let Some(t) = r2_thread {
        let _ = t.join();
    }
    if let Some(e) = io_err {
        return Err(e);
    }
    flush(&mut batch, &mut state)?;
    // ends the closure's borrows of the writer and counters, since it owns nothing to free
    #[allow(clippy::drop_non_drop)]
    drop(flush);
    out.flush().map_err(|e| e.to_string())?;

    #[allow(clippy::drop_non_drop)]
    drop(prep); // ends its borrows, as above
    let pct = 100.0 * n_mapped as f64 / n_reads.max(1) as f64;
    eprintln!(
        "[map] {n_reads} reads in {:.2}s -- {n_mapped} mapped ({pct:.1}%)",
        tmap.elapsed().as_secs_f64()
    );
    if n_bad_qual > 0 {
        eprintln!("[map] warning: {n_bad_qual} FASTQ records had unusable qualities (QUAL '*')");
    }
    Ok(())
}

// Paired-end reads, map -q R1 -q R2
//
// Each mate is first placed on its own exactly as a single-end read, with its best placement and its seeded competitors.
// The library's proper template range comes from the first batch's confident FR pairs, as quartiles +- 3 IQR the way BWA estimates it.
// A pair is then resolved as BWA-MEM resolves pairs, in likelihood form.
//   - The candidates are both mates' best placements and their competitors within 3 edits.
//     If no combination is a proper pair, mate rescue also aligns either mate inside the window its insert size allows next to the partner's best candidates.
//   - Every combination is weighted r1^(ed1 - best1) * r2^(ed2 - best2), with r the edit_ratio, times discordant_prior, worth 3.4 edits, when it is not a proper pair.
//     The most likely combination wins, with ties seeded by the pair's sequence.
//   - A mate's MAPQ is -10 log10 of the weight of the combinations that place it elsewhere, capped like the unpaired MAPQ by hidden-locus odds and seed evidence, so a mate in a repeat whose partner is unique gets the confidence the pair gives it.
//   - The SAM pair fields are 0x1, 0x2 for a proper pair, 0x8 and 0x20 for an unmapped or reverse mate, 0x40 / 0x80, RNEXT / PNEXT, a signed TLEN over the outer coordinates, positive on the leftmost mate, and the mate's CIGAR and MAPQ in MC:Z and MQ:i.
//     An unmapped mate is placed at its partner, as SAM spec 1.4 describes.
// Long mates of CHAIN_MIN_LEN or more keep their own placement, and their pair fields are filled the same way.

/// Prior weight of a placement the mate does not support, a discordant pair, in edits at the read's own error rate, r^3.4 for edit ratio r.
/// BWA-MEM's unpaired penalty is 17 score points, i.e. 3.4 mismatches at 5 points each.
/// A fixed 1e-3 was about 2 edits at 1% error.
/// Real HG002 mates whose partner was placed uniquely then went, as a discordant pair, to a paralog 2 mismatches closer through the sample's own SNPs, which the paralog shares, while GIAB confirms the proper placement.
/// Simulated reads carry no variants, so they never show this.
const UNPAIRED_EDITS: f64 = 3.4;
fn discordant_prior(r: f64) -> f64 {
    r.powf(UNPAIRED_EDITS)
}

/// Proper template lengths [lo, hi] for FR orientation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InsertModel {
    pub lo: i64,
    pub hi: i64,
}

/// The insert model, [max(1, Q1 - 3 IQR), Q3 + 3 IQR] from template lengths of confidently placed pairs, or None below 30 pairs.
pub(crate) fn insert_model(mut tl: Vec<i64>) -> Option<InsertModel> {
    if tl.len() < 30 {
        return None;
    }
    tl.sort_unstable();
    let q = |f: f64| tl[((tl.len() - 1) as f64 * f).round() as usize];
    let (q1, q3) = (q(0.25), q(0.75));
    let iqr = (q3 - q1).max(1);
    Some(InsertModel {
        lo: (q1 - 3 * iqr).max(1),
        hi: q3 + 3 * iqr,
    })
}

/// One mate after the independent pass.
enum MateAln {
    Short(Option<ShortHit>),
    Lines(Vec<String>),
}

/// Places one mate on its own.
/// Long mates and short mates rescued on the chaining path keep their records, while short mates keep their hit for pairing.
fn map_mate(
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
    if hit.is_none() && seq.len() >= RESCUE_MIN_KMERS * ctx.k {
        let l = map_one_read_chain(read_name, seq, ctx, ref_names, RESCUE_CHAIN);
        if !l.is_empty() && !is_unmapped(&l[0]) {
            return MateAln::Lines(l);
        }
    }
    MateAln::Short(hit)
}

/// A placement as pairing sees it, as reference, strand and [start, end).
#[derive(Clone, Copy, Debug)]
struct Place {
    rid: usize,
    rev: bool,
    start: i64,
    end: i64,
}

/// Template length of an FR pair, from the forward mate's start to the reverse mate's end, or None if the two are not FR on one reference.
fn fr_tlen(a: Place, b: Place) -> Option<i64> {
    if a.rid != b.rid || a.rev == b.rev {
        return None;
    }
    let (f, r) = if a.rev { (b, a) } else { (a, b) };
    if f.start > r.end {
        return None;
    }
    Some(r.end - f.start)
}

/// Whether two placements form a proper FR pair under the insert model.
fn proper(a: Place, b: Place, m: &InsertModel) -> bool {
    fr_tlen(a, b).is_some_and(|t| t >= m.lo && t <= m.hi)
}

/// The placement of a short hit of a `len`-base read.
fn hit_place(h: &ShortHit, len: usize) -> Place {
    Place {
        rid: h.rid,
        rev: h.rev,
        start: h.end - len as i64,
        end: h.end,
    }
}

/// Primary-record fields pairing reads back from a SAM line.
struct Prim {
    mapped: bool,
    rname: String,
    pos: i64,
    rev: bool,
    mapq: i32,
    cigar: String,
    end: i64,
}

/// Reads the pairing fields back from a primary SAM record.
fn parse_prim(line: &str) -> Prim {
    let f: Vec<&str> = line.splitn(7, '\t').collect();
    let flag: u16 = f
        .get(1)
        .and_then(|x| x.parse().ok())
        .unwrap_or(SAM_FLAG_UNMAPPED);
    let pos: i64 = f.get(3).and_then(|x| x.parse().ok()).unwrap_or(0);
    let cigar = f.get(5).copied().unwrap_or("*").to_string();
    Prim {
        mapped: flag & SAM_FLAG_UNMAPPED == 0,
        rname: f.get(2).copied().unwrap_or("*").to_string(),
        pos,
        rev: flag & SAM_FLAG_REVERSE != 0,
        mapq: f.get(4).and_then(|x| x.parse().ok()).unwrap_or(0),
        end: pos - 1 + ref_span(&cigar),
        cigar,
    }
}

/// Reference bases a CIGAR covers (M, D, N, =, X).
fn ref_span(cig: &str) -> i64 {
    let (mut n, mut span) = (0i64, 0i64);
    for ch in cig.bytes() {
        if ch.is_ascii_digit() {
            n = n * 10 + (ch - b'0') as i64;
            continue;
        }
        if matches!(ch, b'M' | b'D' | b'N' | b'=' | b'X') {
            span += n;
        }
        n = 0;
    }
    span
}

/// Align `seq`, on strand `rev`, inside the reference window where the partner's insert size puts it.
/// Returns the placement as a ShortHit without competitors yet, or None.
fn rescue_in_window(
    seq: &[u8],
    rc: &[u8],
    partner: Place,
    m: &InsertModel,
    max_ed: i32,
    ctx: &AlignCtx<'_>,
) -> Option<ShortHit> {
    let l = seq.len() as i64;
    let rlen = ctx.ref_seqs[partner.rid].len() as i64;
    let rev = !partner.rev;
    // With the partner forward the mate is reverse and ends in [p.start + lo, p.start + hi].
    // With the partner reverse the mate is forward and starts in [p.end - hi, p.end - lo].
    let (ws, we) = if !partner.rev {
        (partner.start + m.lo - l, partner.start + m.hi)
    } else {
        (partner.end - m.hi, partner.end - m.lo + l)
    };
    let (ws, we) = (ws.max(0), we.min(rlen));
    if we - ws < l {
        return None;
    }
    let q: &[u8] = if rev { rc } else { seq };
    let sep = placement_sep(seq.len()) as i32;
    let text = &ctx.ref_seqs[partner.rid][ws as usize..we as usize];
    let p = crate::placements_x4(q, &[text], max_ed, sep)[0];
    if p.best > max_ed {
        return None;
    }
    let end = ws + p.end as i64;
    let pad = max_ed.max(0) as i64 + WINDOW_PAD as i64;
    let (hs, he) = ((end - l - pad).max(0), (end + pad).min(rlen));
    let (mut comp, mut comp_loc) = (Vec::new(), Vec::new());
    if p.alt <= max_ed {
        comp.push(p.alt);
        comp_loc.push((partner.rid as u32, rev, ws + p.alt_end as i64));
    }
    Some(ShortHit {
        ed: p.best,
        comp,
        comp_loc,
        rep_frac: 0.0,
        seed_ev: [0; 5],
        gap: 0,
        rid: partner.rid,
        rev,
        win_start: hs as usize,
        win_end: he as usize,
        end,
    })
}

/// A candidate placement of a mate, which is its best placement, a competitor its own seeds found, or one rescued next to the partner.
#[derive(Clone, Copy, Debug)]
struct Cand {
    ed: i32,
    rid: usize,
    rev: bool,
    end: i64,
    rescued: bool,
}

/// Competitors further than this many edits from a mate's best placement are not pairing candidates (they weigh < 1e-7).
const MATE_CAND_EDITS: i32 = 3;
/// Most pairing candidates per mate.
const MAX_MATE_CANDS: usize = 64;

/// Candidates of a short mate, the best placement and every competitor within MATE_CAND_EDITS of it, at most MAX_MATE_CANDS.
fn mate_cands(h: &ShortHit) -> Vec<Cand> {
    let mut v = vec![Cand {
        ed: h.ed,
        rid: h.rid,
        rev: h.rev,
        end: h.end,
        rescued: false,
    }];
    for (&e, &(rid, rev, end)) in h.comp.iter().zip(&h.comp_loc) {
        if e <= h.ed + MATE_CAND_EDITS && v.len() < MAX_MATE_CANDS {
            v.push(Cand {
                ed: e,
                rid: rid as usize,
                rev,
                end,
                rescued: false,
            });
        }
    }
    v
}

/// Per-base floor on the edit rate for the sample's own differences from the reference.
/// A real genome's copy of a duplicated region can carry variants its reference copy lacks and a paralog shares, so a read may fit the paralog better by a few edits.
/// Simulated reads, drawn from the reference, never show it.
/// The floor was fitted on half of the real HG002 reads that their pair anchors, as the smallest value at which MAPQ >= 30 is calibrated with 23 wrong against 22.9 implied, and checked on the other half at 23 against 23.1.
/// Without it, reads with 0 edits at a paralog and 3 at their true copy were wrong at MAPQ 40-49 four times as often as implied.
const SAMPLE_DIVERGENCE: f64 = 0.016;

/// The likelihood ratio of one more edit for a read with `ed` edits, i.e. the substitution odds 0.82 e / (1 - e) of the read's estimated error rate e = (ed + 1.26) / len + SAMPLE_DIVERGENCE, with constants fitted as in posterior_mapq_w.
fn edit_ratio(ed: i32, len: usize) -> f64 {
    let e = (ed as f64 + 1.2608) / len.max(1) as f64 + SAMPLE_DIVERGENCE;
    (0.8228 * e / (1.0 - e).max(1e-3)).clamp(1e-6, 0.9)
}

/// Whether two candidates are the same placement, with ends within `sep`.
fn same_place(a: &Cand, b: &Cand, sep: i64) -> bool {
    a.rid == b.rid && a.rev == b.rev && (a.end - b.end).abs() <= sep
}

/// The ShortHit of a candidate.
/// The mate's own best is its hit, and any other candidate is re-verified in a window around its end (the aligner reports the same placement), so the record is built as for any short read.
fn cand_hit(
    c: &Cand,
    own: Option<&ShortHit>,
    seq: &[u8],
    rc: &[u8],
    max_ed: i32,
    ctx: &AlignCtx<'_>,
) -> Option<ShortHit> {
    if let Some(h) = own
        && h.rid == c.rid
        && h.rev == c.rev
        && h.end == c.end
    {
        return Some(clone_hit(h));
    }
    let l = seq.len() as i64;
    let rlen = ctx.ref_seqs[c.rid].len() as i64;
    let pad = max_ed.max(0) as i64 + WINDOW_PAD as i64;
    let (ws, we) = ((c.end - l - pad).max(0), (c.end + pad).min(rlen));
    if we - ws < l {
        return None;
    }
    let q: &[u8] = if c.rev { rc } else { seq };
    let sep = placement_sep(seq.len()) as i32;
    let p = crate::placements_x4(
        q,
        &[&ctx.ref_seqs[c.rid][ws as usize..we as usize]],
        max_ed,
        sep,
    )[0];
    if p.best > max_ed {
        return None;
    }
    Some(ShortHit {
        ed: p.best,
        comp: Vec::new(),
        comp_loc: Vec::new(),
        rep_frac: own.map_or(0.0, |h| h.rep_frac),
        seed_ev: own.map_or([0; 5], |h| h.seed_ev),
        gap: 0,
        rid: c.rid,
        rev: c.rev,
        win_start: ws as usize,
        win_end: we as usize,
        end: ws + p.end as i64,
    })
}

/// A copy of a ShortHit, which holds vectors and so is not Copy.
fn clone_hit(h: &ShortHit) -> ShortHit {
    ShortHit {
        ed: h.ed,
        comp: h.comp.clone(),
        comp_loc: h.comp_loc.clone(),
        rep_frac: h.rep_frac,
        seed_ev: h.seed_ev,
        gap: h.gap,
        rid: h.rid,
        rev: h.rev,
        win_start: h.win_start,
        win_end: h.win_end,
        end: h.end,
    }
}

/// Pair two independently placed mates (see the section comment), with a joint choice over both mates' candidates, mate rescue, pair MAPQ and records.
#[allow(clippy::too_many_arguments)]
fn pair_records(
    name: &str,
    s1: &[u8],
    s2: &[u8],
    a1: MateAln,
    a2: MateAln,
    m: &InsertModel,
    max_ed: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> (Vec<String>, Vec<String>) {
    let (rc1, rc2) = (revcomp(s1), revcomp(s2));
    let seqs: [(&[u8], &[u8]); 2] = [(s1, &rc1), (s2, &rc2)];
    let len = [s1.len(), s2.len()];
    let maxe = [read_max_ed(max_ed, s1.len()), read_max_ed(max_ed, s2.len())];
    let mut lines: [Option<Vec<String>>; 2] = [None, None];
    let mut hits: [Option<ShortHit>; 2] = [None, None];
    for (i, a) in [a1, a2].into_iter().enumerate() {
        match a {
            MateAln::Lines(l) => lines[i] = Some(l),
            MateAln::Short(h) => hits[i] = h,
        }
    }
    let as_place = |c: &Cand, i: usize| Place {
        rid: c.rid,
        rev: c.rev,
        start: c.end - len[i] as i64,
        end: c.end,
    };
    let mut chosen: [Option<(ShortHit, i32, bool)>; 2] = [None, None]; // hit, MAPQ before the cap, rescued
    let mut proper_pick = false;
    if lines[0].is_none() && lines[1].is_none() {
        let mut cands: [Vec<Cand>; 2] = [
            hits[0].as_ref().map_or(Vec::new(), mate_cands),
            hits[1].as_ref().map_or(Vec::new(), mate_cands),
        ];
        // Joint likelihood of every candidate combination, where an unmapped mate contributes a single "absent" option, giving the chosen combination and each mate's chance of being placed wrong
        let resolve = |cands: &[Vec<Cand>; 2]| -> ((Option<usize>, Option<usize>), [f64; 2], bool) {
            let best_ed = [
                cands[0].iter().map(|c| c.ed).min(),
                cands[1].iter().map(|c| c.ed).min(),
            ];
            let ratio = [
                best_ed[0].map(|e| edit_ratio(e, len[0])),
                best_ed[1].map(|e| edit_ratio(e, len[1])),
            ];
            let opts = |i: usize| -> Vec<Option<usize>> {
                if cands[i].is_empty() {
                    vec![None]
                } else {
                    (0..cands[i].len()).map(Some).collect()
                }
            };
            let mut combos: Vec<(Option<usize>, Option<usize>, f64, bool)> = Vec::new();
            for i in opts(0) {
                for j in opts(1) {
                    let mut lk = 1.0;
                    if let Some(i) = i {
                        lk *= ratio[0].unwrap().powi(cands[0][i].ed - best_ed[0].unwrap());
                    }
                    if let Some(j) = j {
                        lk *= ratio[1].unwrap().powi(cands[1][j].ed - best_ed[1].unwrap());
                    }
                    let mut pr = false;
                    if let (Some(i), Some(j)) = (i, j) {
                        pr = proper(as_place(&cands[0][i], 0), as_place(&cands[1][j], 1), m);
                        if !pr {
                            let rr = ratio[0].unwrap_or(0.0).max(ratio[1].unwrap_or(0.0));
                            lk *= discordant_prior(rr);
                        }
                    }
                    combos.push((i, j, lk, pr));
                }
            }
            let top = combos.iter().map(|c| c.2).fold(0.0f64, f64::max);
            let tied: Vec<usize> = (0..combos.len())
                .filter(|&k| combos[k].2 >= top * (1.0 - 1e-9))
                .collect();
            let mut hs = read_hash(s1) ^ read_hash(s2).rotate_left(17);
            hs ^= hs >> 29;
            let pick = combos[tied[(hs % tied.len() as u64) as usize]];
            let total: f64 = combos.iter().map(|c| c.2).sum();
            let mut pw = [1.0f64; 2];
            for (x, ci) in [(0usize, pick.0), (1usize, pick.1)] {
                let Some(ci) = ci else { continue };
                let c = cands[x][ci];
                let sep = placement_sep(len[x]) as i64;
                // the combinations that place this mate elsewhere
                let other: f64 = combos
                    .iter()
                    .filter(|cb| {
                        let k = if x == 0 { cb.0 } else { cb.1 };
                        k.is_none_or(|k| !same_place(&cands[x][k], &c, sep))
                    })
                    .map(|cb| cb.2)
                    .sum();
                pw[x] = (other / total).max(1e-12);
            }
            ((pick.0, pick.1), pw, pick.3)
        };
        let (mut pick, mut pw, mut pr) = resolve(&cands);
        // Mate rescue.
        // Unless the chosen pair is proper and both mates are confident, each mate is aligned inside the window its insert size allows next to the partner's best candidates, all equally good ones up to 4, since its true locus may have lost every seed.
        if !pr || pw.iter().any(|&p| p > 1e-3) {
            let mut added = false;
            for x in 0..2 {
                let anchor = 1 - x;
                if cands[anchor].is_empty() {
                    continue;
                }
                let best = cands[anchor].iter().map(|c| c.ed).min().unwrap();
                let anchors: Vec<Cand> = cands[anchor]
                    .iter()
                    .filter(|c| c.ed == best)
                    .take(4)
                    .copied()
                    .collect();
                let (s, rc) = seqs[x];
                let sep = placement_sep(s.len()) as i64;
                for ac in anchors {
                    if let Some(r) = rescue_in_window(s, rc, as_place(&ac, anchor), m, maxe[x], ctx)
                    {
                        let mut found = vec![Cand {
                            ed: r.ed,
                            rid: r.rid,
                            rev: r.rev,
                            end: r.end,
                            rescued: true,
                        }];
                        if let (Some(&e), Some(&(rid, rev, end))) =
                            (r.comp.first(), r.comp_loc.first())
                            && e <= r.ed + MATE_CAND_EDITS
                        {
                            found.push(Cand {
                                ed: e,
                                rid: rid as usize,
                                rev,
                                end,
                                rescued: true,
                            });
                        }
                        for c in found {
                            if !cands[x].iter().any(|y| same_place(y, &c, sep)) {
                                cands[x].push(c);
                                added = true;
                            }
                        }
                    }
                }
            }
            if added {
                (pick, pw, pr) = resolve(&cands);
            }
        }
        proper_pick = pr;
        for (x, ci) in [(0usize, pick.0), (1usize, pick.1)] {
            let Some(ci) = ci else { continue };
            let c = cands[x][ci];
            let q = ((-10.0 * pw[x].log10()).round() as i32).clamp(0, SAM_MAX_MAPQ);
            let (s, rc) = seqs[x];
            if let Some(h) = cand_hit(&c, hits[x].as_ref(), s, rc, maxe[x], ctx) {
                chosen[x] = Some((h, q, c.rescued));
            }
        }
    }
    // The hidden-locus odds and the seed-evidence cap of the unpaired MAPQ apply to the paired posterior too.
    // In a proper pair whose window next to the partner was aligned, a mate is wrong only if its true copy went unseen and either the partner's did too, or the true pair is discordant.
    // The first means the whole fragment came from another copy, with odds H_x H_p since the mates' seeds are drawn independently, and the second has odds H_x discordant_prior.
    // Where the mate's own hidden-locus odds bind, the window is aligned, and if the mate's placement is the best there, its hidden odds become H_x (H_p + discordant_prior) and the window's runner-up is a verified competitor.
    // Real HG002 pairs inside diverged L1 or segmental-duplication copies kept MAPQ 1-15 where BWA-MEM gives them 40 or more and GIAB confirms the placements, which affected 0.7% of chr20 SNPs.
    let mut alns: [Option<(ShortAln, i32)>; 2] = [None, None];
    let mut hids = [f64::INFINITY; 2];
    let mut window_mq: [Option<i32>; 2] = [None, None];
    let mut xp: [String; 2] = [String::new(), String::new()]; // OPTIK_MAPQ_DEBUG pair inputs
    for x in 0..2 {
        if let Some((h, q, _)) = chosen[x].as_ref() {
            let a = short_align(seqs[x].0, seqs[x].1, h, ctx);
            hids[x] = hidden_odds(h, a.events, len[x], ctx.rho);
            let mq = (*q)
                .min(posterior_mapq_w(h.ed, &[], len[x], hids[x]))
                .min(seed_cap(h));
            if mapq_debug() {
                xp[x] = format!(
                    "pair={},hidden={},seed={}",
                    q,
                    posterior_mapq_w(h.ed, &[], len[x], hids[x]),
                    seed_cap(h)
                );
            }
            alns[x] = Some((a, mq));
        }
    }
    for x in 0..2 {
        let (Some((h, q, _)), Some((p, _, _))) = (chosen[x].as_ref(), chosen[1 - x].as_ref())
        else {
            continue;
        };
        let mq = alns[x].as_ref().map_or(0, |a| a.1);
        if !proper_pick || mq >= (*q).min(seed_cap(h)) {
            continue;
        }
        let partner = Place {
            rid: p.rid,
            rev: p.rev,
            start: p.end - len[1 - x] as i64,
            end: p.end,
        };
        let sep = placement_sep(len[x]) as i64;
        if let Some(r) = rescue_in_window(seqs[x].0, seqs[x].1, partner, m, maxe[x], ctx)
            && r.rid == h.rid
            && r.rev == h.rev
            && (r.end - h.end).abs() <= sep
            && r.ed == h.ed
        {
            let comp: Vec<(i32, f64)> = r.comp.iter().map(|&e| (e, 1.0)).collect();
            let hid = hids[x] * (hids[1 - x] + discordant_prior(edit_ratio(h.ed, len[x])));
            window_mq[x] = Some(
                (*q).min(posterior_mapq_w(h.ed, &comp, len[x], hid))
                    .min(seed_cap(h)),
            );
        }
    }
    for x in 0..2 {
        if let (Some(w), Some(a)) = (window_mq[x], alns[x].as_mut()) {
            a.1 = a.1.max(w);
        }
    }
    // A mate placed by rescue has no evidence of its own.
    // It sits where its partner's placement puts it, so it is no more certain than the partner (a pair inside a segmental duplication whose partner is ambiguous was placed in the wrong copy at MAPQ 36-42).
    let mqs = [alns[0].as_ref().map(|a| a.1), alns[1].as_ref().map(|a| a.1)];
    for x in 0..2 {
        if let (Some((_, _, true)), Some(a), Some(p)) =
            (chosen[x].as_ref(), alns[x].as_mut(), mqs[1 - x])
        {
            a.1 = a.1.min(p);
        }
    }
    // Pairing combines the mates' evidence but cannot create any, so a mate's MAPQ is at most the sum of both mates' unpaired MAPQs, as independent evidence adds in log odds.
    // Two mates that each tie alone lie in a repeat where every copy pair is proper, and the joint candidate lists are then incomplete.
    // Pairs inside chr7's centromeric satellite went to the wrong copy at MAPQ 10-14, 2.3x the implied error.
    let unpaired = |x: usize| {
        hits[x]
            .as_ref()
            .map_or(0, |h| final_mapq(h, h.ed, len[x], ctx.rho))
    };
    let evidence = unpaired(0) + unpaired(1);
    for a in alns.iter_mut().flatten() {
        a.1 = a.1.min(evidence);
    }
    if mapq_debug() {
        for x in 0..2 {
            if !xp[x].is_empty() {
                xp[x] = format!(
                    "\tXP:Z:{},window={},proper={},sum={}",
                    xp[x],
                    window_mq[x].unwrap_or(-1),
                    proper_pick as u8,
                    evidence
                );
            }
        }
    }
    let mut out: [Vec<String>; 2] = [Vec::new(), Vec::new()];
    for x in 0..2 {
        let (s, rc) = seqs[x];
        out[x] = if let Some(l) = lines[x].take() {
            l
        } else if let (Some((h, _, _)), Some((a, mq))) = (chosen[x].take(), alns[x].take()) {
            vec![short_format(name, s, rc, &h, &a, mq, ctx, ref_names) + &xp[x]]
        } else if let Some(h) = hits[x].as_ref() {
            vec![short_record(name, s, h, ctx, ref_names, None)]
        } else {
            let seq_s = if s.is_empty() {
                "*"
            } else {
                std::str::from_utf8(s).unwrap_or("*")
            };
            vec![format!(
                "{name}\t{SAM_FLAG_UNMAPPED}\t*\t0\t0\t*\t*\t0\t0\t{seq_s}\t*"
            )]
        };
    }
    let [mut l1, mut l2] = out;
    let (p1, p2) = (parse_prim(&l1[0]), parse_prim(&l2[0]));
    let proper_pair = p1.mapped && p2.mapped && {
        let pl = |p: &Prim| {
            ref_names
                .iter()
                .position(|n| *n == p.rname)
                .map(|rid| Place {
                    rid,
                    rev: p.rev,
                    start: p.pos - 1,
                    end: p.end,
                })
        };
        matches!((pl(&p1), pl(&p2)), (Some(a), Some(b)) if proper(a, b, m))
    };
    set_pair_fields(&mut l1, &mut l2, proper_pair);
    (l1, l2)
}

/// Fill the SAM pair fields of both mates' records, primary first.
fn set_pair_fields(l1: &mut [String], l2: &mut [String], proper_pair: bool) {
    let (p1, p2) = (parse_prim(&l1[0]), parse_prim(&l2[0]));
    // template length from the primaries' outer coordinates
    let tlen = if p1.mapped && p2.mapped && p1.rname == p2.rname {
        let (s, e) = (p1.pos.min(p2.pos), p1.end.max(p2.end));
        e - s + 1
    } else {
        0
    };
    let first_is_left = p1.pos < p2.pos || (p1.pos == p2.pos && (p1.rev == p2.rev || !p1.rev));
    for (mate_i, (lines, me, other)) in [(&mut *l1, &p1, &p2), (&mut *l2, &p2, &p1)]
        .into_iter()
        .enumerate()
    {
        for (j, line) in lines.iter_mut().enumerate() {
            let mut f: Vec<String> = line.split('\t').map(str::to_string).collect();
            if f.len() < 11 {
                continue;
            }
            let mut flag: u16 = f[1].parse().unwrap_or(0);
            flag |= SAM_FLAG_PAIRED
                | if mate_i == 0 {
                    SAM_FLAG_READ1
                } else {
                    SAM_FLAG_READ2
                };
            if proper_pair && me.mapped && other.mapped {
                flag |= SAM_FLAG_PROPER_PAIR;
            }
            if !other.mapped {
                flag |= SAM_FLAG_MATE_UNMAPPED;
            } else if other.rev {
                flag |= SAM_FLAG_MATE_REVERSE;
            }
            let this_mapped = flag & SAM_FLAG_UNMAPPED == 0;
            if !this_mapped && other.mapped {
                // an unmapped mate takes its partner's position, per SAM 1.4
                f[2] = other.rname.clone();
                f[3] = other.pos.to_string();
            }
            f[1] = flag.to_string();
            if other.mapped {
                f[6] = if f[2] == other.rname {
                    "=".into()
                } else {
                    other.rname.clone()
                };
                f[7] = other.pos.to_string();
            } else if this_mapped {
                f[6] = "=".into();
                f[7] = f[3].clone();
            }
            let left = if mate_i == 0 {
                first_is_left
            } else {
                !first_is_left
            };
            f[8] = if j == 0 && tlen != 0 {
                (if left { tlen } else { -tlen }).to_string()
            } else {
                "0".into()
            };
            if other.mapped {
                f.push(format!("MC:Z:{}", other.cigar));
                f.push(format!("MQ:i:{}", other.mapq));
            }
            *line = f.join("\t");
        }
    }
}

/// Mate names, where a trailing /1 or /2 is not part of the read's name
fn mate_name(n: &str) -> &str {
    n.strip_suffix("/1")
        .or_else(|| n.strip_suffix("/2"))
        .unwrap_or(n)
}

#[cfg(test)]
mod tests {
    #[test]
    fn local_pieces_clips_and_splits() {
        use super::{MM_ASM20, MM_READS, asm_score, local_pieces};
        let offs = |cols: &[u8]| {
            let (mut qp, mut rp) = (vec![0usize], vec![0usize]);
            for &c in cols {
                qp.push(qp.last().unwrap() + (c != b'D') as usize);
                rp.push(rp.last().unwrap() + (c != b'I') as usize);
            }
            (qp, rp)
        };
        let k = 15;
        // a clean path is one piece over every column
        let q = rand_dna(400, 1);
        let cols = vec![b'M'; 400];
        let (qp, rp) = offs(&cols);
        let an = [0, 100, 200, 300];
        assert_eq!(
            local_pieces(&cols, &q, &q, &qp, &rp, &an, k, MM_READS),
            vec![(0, 400, 800)]
        );
        // a junk end past the last anchor is clipped back to the last match before it
        let mut r = q.clone();
        for i in 386..400 {
            if i % 7 != 0 {
                r[i] = if q[i] == b'A' { b'C' } else { b'A' };
            }
        }
        let p = local_pieces(&cols, &q, &r, &qp, &rp, &an, k, MM_READS);
        assert!(p.len() == 1 && p[0].0 == 0 && p[0].1 <= 386, "{p:?}");
        // One long deletion between anchors is neither a Z-drop nor clipped, so it stays one piece with all columns.
        let mut cols2 = vec![b'M'; 200];
        cols2.extend(std::iter::repeat_n(b'D', 300));
        cols2.extend(std::iter::repeat_n(b'M', 200));
        let r2: Vec<u8> = q[..200]
            .iter()
            .chain(rand_dna(300, 2).iter())
            .chain(q[200..].iter())
            .copied()
            .collect();
        let (qp2, rp2) = offs(&cols2);
        let p2 = local_pieces(
            &cols2,
            &q,
            &r2,
            &qp2,
            &rp2,
            &[0, 100, 250, 350],
            k,
            MM_READS,
        );
        assert!(p2.len() == 1 && p2[0].0 == 0 && p2[0].1 == 700, "{p2:?}");
        // 300 columns of unrelated sequence between two anchored blocks split into two pieces.
        let q3: Vec<u8> = q[..150]
            .iter()
            .chain(rand_dna(300, 3).iter())
            .chain(q[150..300].iter())
            .copied()
            .collect();
        let r3: Vec<u8> = q[..150]
            .iter()
            .chain(rand_dna(300, 4).iter())
            .chain(q[150..300].iter())
            .copied()
            .collect();
        let cols3 = vec![b'M'; 600];
        let (qp3, rp3) = offs(&cols3);
        let p3 = local_pieces(
            &cols3,
            &q3,
            &r3,
            &qp3,
            &rp3,
            &[0, 60, 120, 460, 520, 580],
            k,
            MM_READS,
        );
        assert!(
            p3.len() == 2 && p3[0].0 == 0 && p3[0].1 <= 154 && p3[1].0 >= 446 && p3[1].1 == 600,
            "{p3:?}"
        );
        // Under asm20 scoring the long deletion is still one piece, thanks to the two-piece gap's e2 allowance.
        let p4 = local_pieces(
            &cols2,
            &q,
            &r2,
            &qp2,
            &rp2,
            &[0, 100, 250, 350],
            k,
            MM_ASM20,
        );
        assert!(p4.len() == 1 && p4[0].0 == 0 && p4[0].1 == 700, "{p4:?}");
        // A clean 400 bp block scores 400, and the same block at 19% substitutions scores under 200.
        assert_eq!(asm_score("400M", &q, &q, 0), 400);
        let mut d = q.clone();
        for i in (0..400).filter(|i| i % 16 < 3) {
            d[i] = if q[i] == b'A' { b'C' } else { b'A' };
        }
        assert!(asm_score("400M", &q, &d, 0) < super::ASM_MIN_DP);
        assert_eq!(
            asm_score(
                "10S100M30D100M",
                q[..210].to_vec().as_slice(),
                &[&q[10..110], &d[..30], &q[110..210]].concat(),
                0
            ),
            200 - 26 - 30
        );
    }

    #[test]
    fn hard_clip_supplementary() {
        use super::hard_clip;
        assert_eq!(
            hard_clip("3S4M2S", "AAACCCCGG"),
            ("3H4M2H".to_string(), "CCCC")
        );
        assert_eq!(
            hard_clip("5M2I3M", "ACGTACGTAC"),
            ("5M2I3M".to_string(), "ACGTACGTAC")
        );
        assert_eq!(
            hard_clip("10S5M", "AAAAAAAAAACCCCC"),
            ("10H5M".to_string(), "CCCCC")
        );
        assert_eq!(
            hard_clip("5M1D2M12S", "CCCCCGGAAAAAAAAAAAA"),
            ("5M1D2M12H".to_string(), "CCCCCGG")
        );
    }

    /// QUAL follows SEQ on hard-clipped supplementary records, cut by the CIGAR's H in the record's orientation, with qualities reversed on reverse records.
    #[test]
    fn qual_of_hard_clipped_records() {
        use super::{fill_qual, hard_clips};
        assert_eq!(hard_clips("3H4M2H"), (3, 2));
        assert_eq!(hard_clips("10H5M"), (10, 0));
        assert_eq!(hard_clips("5M12H"), (0, 12));
        assert_eq!(hard_clips("3S4M2S"), (0, 0));
        // Read qualities are "ABCDEFGHI", and forward record 3H4M2H holds bases 3..7.
        let mut l = vec!["r\t2048\tc\t1\t60\t3H4M2H\t*\t0\t0\tCCCC\t*".to_string()];
        fill_qual(&mut l, b"ABCDEFGHI");
        assert!(l[0].ends_with("\tCCCC\tDEFG"), "{}", l[0]);
        // A reverse record has its qualities reversed ("IHGFEDCBA"), then cut 2H..1H.
        let mut l = vec!["r\t2064\tc\t1\t60\t2H6M1H\t*\t0\t0\tAAAAAA\t*".to_string()];
        fill_qual(&mut l, b"ABCDEFGHI");
        assert!(l[0].ends_with("\tAAAAAA\tGFEDCB"), "{}", l[0]);
        // An unclipped primary keeps the whole quality string.
        let mut l = vec!["r\t0\tc\t1\t60\t9M\t*\t0\t0\tAAAAAAAAA\t*".to_string()];
        fill_qual(&mut l, b"ABCDEFGHI");
        assert!(l[0].ends_with("\tABCDEFGHI"), "{}", l[0]);
    }

    use super::*;

    /// Insert model from quartiles +- 3 IQR, with none below 30 pairs.
    #[test]
    fn insert_model_from_quartiles() {
        assert_eq!(insert_model(vec![400; 29]), None);
        let tl: Vec<i64> = (0..100).map(|i| 300 + 2 * i).collect(); // 300..498
        let m = insert_model(tl).unwrap();
        // Q1 = 350 and Q3 = 448 by rounded ranks, so IQR = 98.
        assert_eq!((m.lo, m.hi), (350 - 294, 448 + 294));
    }

    /// FR template length, from forward start to reverse end on the same reference.
    #[test]
    fn fr_template_length() {
        let f = Place {
            rid: 0,
            rev: false,
            start: 100,
            end: 250,
        };
        let r = Place {
            rid: 0,
            rev: true,
            start: 350,
            end: 500,
        };
        assert_eq!(fr_tlen(f, r), Some(400));
        assert_eq!(fr_tlen(r, f), Some(400));
        assert_eq!(fr_tlen(f, Place { rev: false, ..r }), None); // FF
        assert_eq!(fr_tlen(f, Place { rid: 1, ..r }), None); // other contig
        assert_eq!(
            fr_tlen(
                Place {
                    start: 600,
                    end: 750,
                    ..f
                },
                r
            ),
            None
        ); // RF (outward)
        let m = InsertModel { lo: 200, hi: 600 };
        assert!(proper(f, r, &m));
        assert!(!proper(
            f,
            Place {
                start: 950,
                end: 1100,
                ..r
            },
            &m
        ));
    }

    #[test]
    fn mate_names() {
        assert_eq!(mate_name("r7/1"), "r7");
        assert_eq!(mate_name("r7/2"), "r7");
        assert_eq!(mate_name("r7"), "r7");
        assert_eq!(mate_name("r7/3"), "r7/3");
    }

    /// SAM pair fields, i.e. flags, RNEXT/PNEXT, signed TLEN, MC/MQ, and an unmapped mate placed at its partner.
    #[test]
    fn pair_fields() {
        let mut a = vec!["p\t0\tc\t101\t60\t150M\t*\t0\t0\tACGT\t*\tNM:i:0".to_string()];
        let mut b = vec!["p\t16\tc\t351\t40\t140M10S\t*\t0\t0\tACGT\t*\tNM:i:1".to_string()];
        set_pair_fields(&mut a, &mut b, true);
        let fa: Vec<&str> = a[0].split('\t').collect();
        let fb: Vec<&str> = b[0].split('\t').collect();
        assert_eq!(fa[1], (1 | 2 | 0x20 | 0x40).to_string());
        assert_eq!(fb[1], (1 | 2 | 16 | 0x80).to_string());
        assert_eq!((fa[6], fa[7], fa[8]), ("=", "351", "390")); // 101..490
        assert_eq!((fb[6], fb[7], fb[8]), ("=", "101", "-390"));
        assert!(
            a[0].ends_with("\tMC:Z:140M10S\tMQ:i:40") && b[0].ends_with("\tMC:Z:150M\tMQ:i:60")
        );
        // With one mate unmapped it takes the partner's position, with flags 0x8 / 0x4.
        let mut a = vec!["p\t0\tc\t101\t60\t150M\t*\t0\t0\tACGT\t*\tNM:i:0".to_string()];
        let mut u = vec!["p\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*".to_string()];
        set_pair_fields(&mut a, &mut u, false);
        let fa: Vec<&str> = a[0].split('\t').collect();
        let fu: Vec<&str> = u[0].split('\t').collect();
        assert_eq!(fa[1], (1 | 8 | 0x40).to_string());
        assert_eq!((fa[6], fa[7], fa[8]), ("=", "101", "0"));
        assert_eq!(fu[1], (1 | 4 | 0x80).to_string());
        assert_eq!(
            (fu[2], fu[3], fu[6], fu[7], fu[8]),
            ("c", "101", "=", "101", "0")
        );
    }

    /// Deterministic xorshift64 PRNG, avoiding an external dev-dependency.
    fn rng_state(seed: u64) -> u64 {
        let mut s = seed;
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    }
    fn next_rng(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }

    const BASES: &[u8] = b"ACGT";

    /// Random DNA window of `len` bases from a deterministic PRNG.
    fn rand_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut s = rng_state(seed);
        (0..len)
            .map(|_| BASES[(next_rng(&mut s) % 4) as usize])
            .collect()
    }

    /// Reconstruct the semi-global edit distance a CIGAR implies by walking it from (query=0, text=text_start).
    /// Returns None if the CIGAR overshoots.
    fn cigar_edit_distance(
        cigar: &str,
        query: &[u8],
        text: &[u8],
        text_start: usize,
    ) -> Option<usize> {
        let mut i = 0usize; // query offset
        let mut t = text_start; // text offset
        let mut ed = 0usize;
        let mut num = 0usize;
        for ch in cigar.bytes() {
            match ch {
                b'0'..=b'9' => num = num * 10 + (ch - b'0') as usize,
                _ => {
                    if num == 0 {
                        return None;
                    }
                    match ch {
                        b'M' => {
                            for _ in 0..num {
                                if i >= query.len() || t >= text.len() {
                                    return None;
                                }
                                if query[i] != text[t] {
                                    ed += 1;
                                }
                                i += 1;
                                t += 1;
                            }
                        }
                        b'I' => {
                            for _ in 0..num {
                                if i >= query.len() {
                                    return None;
                                }
                                ed += 1;
                                i += 1;
                            }
                        }
                        b'D' => {
                            for _ in 0..num {
                                if t >= text.len() {
                                    return None;
                                }
                                ed += 1;
                                t += 1;
                            }
                        }
                        _ => return None,
                    }
                    num = 0;
                }
            }
        }
        if i != query.len() {
            return None;
        }
        Some(ed)
    }

    /// Independent O(mn) semi-global edit distance (query end-to-end, text ends free).
    fn naive_semi_global_ed(q: &[u8], t: &[u8]) -> usize {
        let mut prev = vec![0usize; t.len() + 1]; // free text start
        for i in 1..=q.len() {
            let mut cur = vec![i; t.len() + 1];
            for j in 1..=t.len() {
                let sub = prev[j - 1] + (q[i - 1] != t[j - 1]) as usize;
                cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
            }
            prev = cur;
        }
        *prev.iter().min().unwrap() // free text end
    }

    #[test]
    fn left_align_normalizes_and_preserves() {
        // Random reads with planted indels in repeat contexts.
        // left_align keeps the same columns, i.e. the ops multiset, the mismatches and both strings fully consumed, and leaves no gap that could still move left.
        let expand = |c: &str| -> Vec<u8> {
            let mut v = Vec::new();
            let mut n = 0usize;
            for ch in c.bytes() {
                if ch.is_ascii_digit() {
                    n = n * 10 + (ch - b'0') as usize;
                } else {
                    v.extend(std::iter::repeat_n(ch, n));
                    n = 0;
                }
            }
            v
        };
        let mut st = rng_state(77);
        for case in 0..400 {
            // text with homopolymers / short tandem repeats
            let mut t = Vec::new();
            while t.len() < 200 {
                let unit: Vec<u8> = (0..(1 + next_rng(&mut st) % 3))
                    .map(|_| b"ACGT"[(next_rng(&mut st) % 4) as usize])
                    .collect();
                for _ in 0..(1 + next_rng(&mut st) % 4) {
                    t.extend_from_slice(&unit);
                }
            }
            let a = 20 + (next_rng(&mut st) % 40) as usize;
            let mut q: Vec<u8> = t[a..a + 120].to_vec();
            let at = 30 + (next_rng(&mut st) % 60) as usize;
            let l = 1 + (next_rng(&mut st) % 4) as usize;
            let del = case % 2 == 0;
            if del {
                q.drain(at..at + l);
            } else {
                let ins: Vec<u8> = q[at - l..at].to_vec();
                for (j, b) in ins.into_iter().enumerate() {
                    q.insert(at + j, b);
                }
            }
            let Some((cig, off)) = crate::wfa_affine_cigar(&q, &t) else {
                continue;
            };
            let out = left_align(&cig, &q, &t, off);
            let (e0, e1) = (expand(&cig), expand(&out));
            let count = |v: &[u8], o: u8| v.iter().filter(|&&x| x == o).count();
            for o in *b"MID" {
                assert_eq!(count(&e0, o), count(&e1, o), "case {case}: {cig} -> {out}");
            }
            assert_eq!(
                cigar_ed_implied(&cig, &q, &t, off),
                cigar_ed_implied(&out, &q, &t, off),
                "case {case}"
            );
            // No gap can still shift, i.e. the column before each gap run is not a match equal to the gap's last base.
            let (mut qi, mut ti) = (0usize, off);
            let mut i = 0;
            while i < e1.len() {
                let o = e1[i];
                if o == b'M' {
                    qi += 1;
                    ti += 1;
                    i += 1;
                    continue;
                }
                let mut j = i;
                while j < e1.len() && e1[j] == o {
                    j += 1;
                }
                let len = j - i;
                if i > 0 && e1[i - 1] == b'M' {
                    let (was, last) = if o == b'D' {
                        (t[ti - 1], t[ti + len - 1])
                    } else {
                        (q[qi - 1], q[qi + len - 1])
                    };
                    assert_ne!(last, was, "case {case}: gap not left-aligned in {out}");
                }
                if o == b'D' {
                    ti += len;
                } else {
                    qi += len;
                }
                i = j;
            }
        }
    }

    #[test]
    fn left_align_moves_gap_past_mismatch_in_run() {
        // HG002 chr20:59112821-59112830 holds a 1-base deletion and an A>T SNP in an A run.
        // The reference reads GATT A12 TGT and the read GATT A7 T A3 TGT.
        // The aligner may write the T as a mismatch beside the deletion, and normalization puts the deletion at the start of the run and the T at the SNP, as GIAB and BWA-MEM write it.
        let t = b"GATTAAAAAAAAAAAATGT";
        let q = b"GATTAAAAAAATAAATGT";
        assert_eq!(left_align("12M1D6M", q, t, 0), "4M1D14M");
        assert_eq!(
            cigar_ed_implied("12M1D6M", q, t, 0),
            cigar_ed_implied("4M1D14M", q, t, 0)
        );
    }

    #[test]
    fn polish_gap_clusters_replaces_cancelling_indels() {
        // The read is the reference with 10 bases inserted at 100.
        // A path forced through 5D ... 15I is replaced by a cheaper free alignment of the window that consumes the same bases with fewer edits.
        let mut st = rng_state(5);
        let r: Vec<u8> = (0..200)
            .map(|_| b"ACGT"[(next_rng(&mut st) % 4) as usize])
            .collect();
        let mut q = r[..100].to_vec();
        q.extend((0..10).map(|_| b"ACGT"[(next_rng(&mut st) % 4) as usize]));
        q.extend_from_slice(&r[100..]);
        let ops = vec![(95u32, b'M'), (5, b'D'), (10, b'M'), (15, b'I'), (90, b'M')];
        let (out, saved) = polish_gap_clusters(&ops, &q, &r, 0, 0);
        assert!(saved > 0, "{out:?}");
        let cig = |v: &[(u32, u8)]| {
            v.iter()
                .map(|(n, o)| format!("{n}{}", *o as char))
                .collect::<String>()
        };
        let qn: u32 = out.iter().filter(|o| o.1 != b'D').map(|o| o.0).sum();
        let rn: u32 = out.iter().filter(|o| o.1 != b'I').map(|o| o.0).sum();
        assert_eq!((qn as usize, rn as usize), (q.len(), r.len()));
        let (before, after) = (
            cigar_ed_implied(&cig(&ops), &q, &r, 0).unwrap(),
            cigar_ed_implied(&cig(&out), &q, &r, 0).unwrap(),
        );
        assert!(after < before, "{before} -> {after}: {}", cig(&out));
        // an alignment with a single gap is left alone
        let one = vec![(100u32, b'M'), (10, b'I'), (100, b'M')];
        assert_eq!(polish_gap_clusters(&one, &q, &r, 0, 0).1, 0);
        // 1-2 bp gaps from sequencing noise do not form clusters by themselves
        let noise = vec![
            (95u32, b'M'),
            (1, b'D'),
            (20, b'M'),
            (2, b'I'),
            (20, b'M'),
            (1, b'D'),
            (63, b'M'),
        ];
        assert_eq!(polish_gap_clusters(&noise, &q, &r, 0, 0).1, 0);
    }

    /// The scalar CIGAR DP (the fallback when the wavefront aligner gives up) on hand-checked cases.
    #[test]
    fn scalar_cigar_trivials() {
        // Exact window placement -> single M run, zero edits, start = 5.
        let window = b"AAAAACGGTACTTCCG".to_vec();
        let (sc, ss) = cigar(b"CGGTAC", &window);
        assert_eq!((sc.as_str(), ss), ("6M", 5));

        // For a query with no exact placement the CIGAR must consume the query and imply the optimal (minimum) semi-global edit distance.
        let read2 = b"ACGTACGT".to_vec();
        let (sc2, ss2) = cigar(&read2, &window);
        assert_eq!(
            cigar_edit_distance(&sc2, &read2, &window, ss2),
            Some(naive_semi_global_ed(&read2, &window)),
            "cigar={sc2} start={ss2}"
        );

        // Empty query is a special-cased "*".
        let (sc3, ss3) = cigar(b"", &window);
        assert_eq!((sc3.as_str(), ss3), ("*", 0));
    }
}
