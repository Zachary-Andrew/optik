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

mod sam;
mod cigar;
mod tune;
mod long_read;
mod short_read;
mod mapq;
mod paired;

pub(crate) use sam::*;
use cigar::*;
use tune::*;
pub(crate) use tune::{NOISY_ERROR, NOISY_K};
use long_read::*;
use short_read::*;
use mapq::*;
pub(crate) use paired::*;

/// Number of distinct byte values, the size of a byte lookup table.
const BYTE_VALUES: usize = 256;
/// Multiplier from a fraction to a percentage in log messages.
const PERCENT: f64 = 100.0;
/// Bytes per MB (MiB) in log messages.
const BYTES_PER_MB: f64 = 1048576.0;
/// The two mates of a read pair.
const MATES: usize = 2;
/// Score `-1` that callers of crate::wfa_affine_align pass for "no score cap".
const NO_SCORE_CAP: i32 = -1;

/// True when the WFA CIGAR must be bypassed in favour of the scalar unit-cost DP, which `OPTIK_FORCE_SCALAR_CIGAR` set to anything but "0" requests.
/// It is a safety switch and an A/B reference for the short-read path.
fn force_scalar_cigar() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| std::env::var_os("OPTIK_FORCE_SCALAR_CIGAR").is_some_and(|val| val != "0"))
}

/// Whether environment variable `name` is set, read once into `cell`, since a getenv per read showed in profiles.
fn env_flag(name: &'static str, cell: &'static std::sync::OnceLock<bool>) -> bool {
    *cell.get_or_init(|| std::env::var_os(name).is_some())
}
static DEBUG_CHAIN: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// Reverse complement of a DNA sequence, with any non-ACGT byte turned into N.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    // Table lookup with no per-base branch, turning non-ACGT into N.
    static RC: [u8; 256] = {
        let mut table = [b'N'; BYTE_VALUES];

        table[b'A' as usize] = b'T';
        table[b'a' as usize] = b'T';
        table[b'C' as usize] = b'G';
        table[b'c' as usize] = b'G';
        table[b'G' as usize] = b'C';
        table[b'g' as usize] = b'C';
        table[b'T' as usize] = b'A';
        table[b't' as usize] = b'A';
        table
    };
    seq.iter().rev().map(|&base| RC[base as usize]).collect()
}

/// Default per-read edit budget, used without `-e`, as a share of the read length.
/// A fixed cap rejects every read beyond it regardless of length, while divergence scales with length.
/// 15% of the read keeps alignments far from random-sequence identity, since an unrelated 150 bp read reaches a best semi-global identity of about 55-60%.
/// Candidates come only from exact seed hits, and 0 of 50k random reads map.
pub(crate) const AUTO_ED_RATE: f64 = 0.15;
/// Smallest automatic edit budget.
pub(crate) const MIN_AUTO_ED: i32 = 5;

/// The edit budget of a read of `len` bases, `max_ed` if given (>= 0) and otherwise AUTO_ED_RATE of the length, at least MIN_AUTO_ED.
fn read_max_ed(max_ed: i32, len: usize) -> i32 {
    if max_ed >= 0 {
        max_ed
    } else {
        ((len as f64 * AUTO_ED_RATE).ceil() as i32).max(MIN_AUTO_ED)
    }
}

/// Deterministic per-read hash, 64-bit FNV-1a followed by MurmurHash3's fmix64 steps, that seeds the choice among equally good placements.
/// The same read thus gets the same copy at any thread count.
fn read_hash(seq: &[u8]) -> u64 {
    const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    /// MurmurHash3 fmix64 finalizer shift and multiplier.
    const FMIX64_SHIFT: u32 = 33;
    const FMIX64_MULT: u64 = 0xff51_afd7_ed55_8ccd;
    let mut hash: u64 = FNV_OFFSET_BASIS;

    for &byte in seq {
        hash = (hash ^ byte as u64).wrapping_mul(FNV_PRIME);
    }

    hash ^= hash >> FMIX64_SHIFT;
    hash = hash.wrapping_mul(FMIX64_MULT);
    hash ^= hash >> FMIX64_SHIFT;
    hash
}

/// Immutable lookup context shared by every read's alignment attempt.
/// `seeds` is the minimizer seed index of include/seed_index.hpp, whose hits record the strand on which the read matches the reference.
struct AlignCtx<'a> {
    kmer_len: usize,
    seeds: &'a SeedIdx<'a>,
    ref_seqs: &'a [Vec<u8>],
    /// The cs tag form, with MD always written.
    cs_mode: CsMode,
    /// The reference's minimizer repeat fraction from SeedIdx::repeat_fraction, a MAPQ input of hidden_odds.
    rho: f64,
}

/// Seed occurrence cap of the mapper's index, SeedIndex max_occ.
/// Minimizers occurring more often are dropped from short-read seeding.
const SEED_MAX_OCC: u32 = 1000;
/// Short reads over their edit budget are retried on the chaining path if they are at least this many k-mers long.
const RESCUE_MIN_KMERS: usize = 4;

/// Runs `optik map`, mapping single-end reads, or pairs when `reads2` is given, to the reference at `ref_path` and writing SAM to `output`, or stdout if None.
///
/// `k_auto` chooses k from the first batch, 15 for noisy long reads and `kmer_len` otherwise.
/// `window` 0 tunes the window from the read length.
/// `max_ed` < 0 uses the automatic edit budget.
/// `cmdline` goes into the @PG header line.
#[allow(clippy::too_many_arguments)]
/// Maps one unpaired read.
/// Long reads go straight to the chaining path, while short ones are verified against candidate loci first.
fn map_single_read(read_name: &str, seq: &[u8], max_ed: i32, ctx: &AlignCtx<'_>, ref_names: &[String]) -> Vec<String> {
    if seq.len() >= CHAIN_MIN_LEN {
        map_one_read_chain(read_name, seq, ctx, ref_names, LONG_CHAIN)
    } else {
        let rec = map_one_read_short(
            read_name,
            seq,
            read_max_ed(max_ed, seq.len()),
            ctx,
            ref_names,
        );

        // Short reads beyond -e, i.e. high-error ones, get a second chance on the chaining path, which has no edit-distance cap.
        if is_unmapped(&rec) && seq.len() >= RESCUE_MIN_KMERS * ctx.kmer_len {
            map_one_read_chain(read_name, seq, ctx, ref_names, RESCUE_CHAIN)
        } else {
            vec![rec]
        }
    }
}

/// Estimates the insert-size model from the first batch, or falls back to FALLBACK_INSERT_MODEL when too few pairs qualify.
fn estimate_insert_model(firsts: &[(MateAln, MateAln)], batch: &[MapRead], ctx: &AlignCtx<'_>) -> InsertModel {
    // Template lengths of the first batch's confidently placed FR pairs, with an unpaired posterior >= 30 for both and no repetitive seeds.
    let sure = |hit: &ShortHit, len: usize| {
        hit.rep_frac == 0.0
            && final_mapq(hit, hit.edit_dist, len, ctx.rho) >= CONFIDENT_MAPQ
    };
    // Only pairs with both mates placed uniquely and confidently teach the insert model. Repeats would skew it.
    let tlens: Vec<i64> = firsts
        .iter()
        .zip(batch.iter())
        .filter_map(|((aln1, aln2), (_, seq1, _, mate))| {
            let (MateAln::Short(Some(hit1)), MateAln::Short(Some(hit2))) = (aln1, aln2)
            else {
                return None;
            };

            if !sure(hit1, seq1.len()) || !sure(hit2, mate.as_ref().expect("paired mode always reads a mate").0.len()) {
                return None;
            }

            let tlen = fr_tlen(
                hit_place(hit1, seq1.len()),
                hit_place(hit2, mate.as_ref().expect("paired mode always reads a mate").0.len()),
            )?;

            (tlen > 0 && tlen < MAX_TEMPLATE_LEN).then_some(tlen)
        })
        .collect();
    let n_pairs = tlens.len();
    let est = insert_model(tlens);
    let est_model = est.unwrap_or(FALLBACK_INSERT_MODEL);

    if est.is_some() {
        eprintln!(
            "[map] paired: proper template length {}..{} (from {n_pairs} confident pairs)",
            est_model.min_tlen, est_model.max_tlen
        );
    } else {
        eprintln!(
            "[map] paired: only {n_pairs} confident FR pairs in the first batch; proper template \
                   length assumed {}..{}",
            est_model.min_tlen, est_model.max_tlen
        );
    }

    est_model
}

/// R2 records are handed over in blocks of this many.
const R2_BLOCK: usize = 4096;
/// One R2 record as (name, sequence, qualities).
type R2Record = (String, Vec<u8>, Vec<u8>);
/// A block of R2 records, or the reader's error.
type R2Block = Result<Vec<R2Record>, String>;

/// Parses R2 on its own thread, one block per request.
/// The thread waits for a request on `want_rx` before parsing each block, so it never runs while the mapping pool does unless the caller asks ahead.
fn spawn_r2_reader(
    path2: &str,
    want_rx: std::sync::mpsc::Receiver<()>,
    block_tx: std::sync::mpsc::SyncSender<R2Block>,
) -> std::thread::JoinHandle<()> {
    let path2 = path2.to_string();

    std::thread::spawn(move || {
        let mut buf: Vec<R2Record> = Vec::with_capacity(R2_BLOCK);
        let mut asked = want_rx.recv().is_ok();
        // A dropped request channel before the first ask means the mapper never needs R2.
        let mut stop = !asked;
        let read_result = iterate_reads(&path2, |name, seq, qual| {
            if stop {
                return;
            }

            buf.push((name.to_string(), seq.to_vec(), qual.to_vec()));

            // A failed send means the mapper hung up early. Nothing is left to report to.
            if buf.len() == R2_BLOCK {
                let _ = block_tx.send(Ok(std::mem::replace(
                    &mut buf,
                    Vec::with_capacity(R2_BLOCK),
                )));

                asked = want_rx.recv().is_ok(); // wait for the next request or the end
                stop = !asked;
            }
        });

        if !stop {
            // The last, partial block, then the error, or the end once the sender drops.
            if !buf.is_empty() {
                let _ = block_tx.send(Ok(buf));
                let _ = want_rx.recv();
            }

            if let Err(err) = read_result {
                let _ = block_tx.send(Err(err));
            }
        }
    })
}

// The R2 reader gets one block asked for at a time.
// With -t >= 2 the next block is asked for as soon as one arrives, so R2 is parsed beside R1 by two reading threads while the mapping pool is idle.
// An outstanding block is collected before every mapping batch, so the R2 reader is idle while the pool runs and the run never has more than -t threads busy.
struct R2Feed {
    want: std::sync::mpsc::SyncSender<()>,
    block_rx: std::sync::mpsc::Receiver<R2Block>,
    buf: std::collections::VecDeque<R2Record>,
    pending: bool,
    end: bool,
    ahead: bool,
}

impl R2Feed {
    /// Asks the reader for the next block, unless one is already on its way or the file has ended.
    fn ask(&mut self) {
        if !self.pending && !self.end {
            if self.want.send(()).is_err() {
                self.end = true;
            } else {
                self.pending = true;
            }
        }
    }

    /// Waits for the block that was asked for, if any, and returns the reader's error when it hit one.
    fn collect(&mut self) -> Option<String> {
        if !self.pending {
            return None;
        }

        self.pending = false;

        match self.block_rx.recv() {
            Ok(Ok(block)) => {
                self.buf.extend(block);
                None
            }
            Ok(Err(err)) => Some(err),
            Err(_) => {
                self.end = true;
                None
            }
        }
    }

    /// Returns the next R2 record, fetching a new block when the buffer runs dry.
    fn next(&mut self) -> Option<Result<R2Record, String>> {
        if self.buf.is_empty() {
            self.ask();

            if let Some(err) = self.collect() {
                return Some(Err(err));
            }

            if self.ahead {
                self.ask();
            }
        }

        self.buf.pop_front().map(Ok)
    }
}

/// Cleans each input record before mapping, so every one of them yields exactly one valid SAM record.
/// It warns once per kind of problem and counts the records whose qualities had to be dropped.
#[derive(Default)]
struct ReadPrep {
    qname_warned: bool,
    qual_warned: bool,
    bad_qual_records: u64,
}

impl ReadPrep {
    /// Uppercases the bases (reading RNA U as T), truncates over-long names and drops unusable qualities.
    fn prepare(&mut self, read_name: &str, seq: &[u8], qual: &[u8]) -> (String, Vec<u8>, Vec<u8>) {
        // SAM QUAL must be '*' or exactly as long as SEQ.
        let qual: Vec<u8> = if qual.is_empty() {
            Vec::new()
        } else if qual.len() == seq.len()
            && qual
                .iter()
                .all(|&qual_byte| (b'!'..=b'~').contains(&qual_byte))
        {
            qual.to_vec()
        } else {
            // A wrong length, or a byte outside Phred+33 '!'..'~' that samtools rejects, gives QUAL '*' instead of a broken SAM.
            if !self.qual_warned {
                eprintln!(
                    "[map] warning: FASTQ record '{read_name}' has {} quality values for {} \
                               bases or a value outside '!'..'~'; its QUAL is written as '*'",
                    qual.len(),
                    seq.len()
                );
                self.qual_warned = true;
            }

            self.bad_qual_records += 1;
            Vec::new()
        };
        let seq: Vec<u8> = seq
            .iter()
            .map(|&base| match base.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect();
        // SAM QNAME holds at most 254 characters, and samtools rejects longer ones.
        let name = if read_name.len() > MAX_QNAME {
            if !self.qname_warned {
                eprintln!(
                    "[map] warning: read names longer than {MAX_QNAME} characters \
                           are truncated (SAM QNAME limit)"
                );
                self.qname_warned = true;
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
    }
}

/// Reads the reference into memory as named sequences, normalised for alignment.
/// Soft-masked lowercase bases count as the same bases, RNA U becomes T, and any other non-ACGT code becomes N.
/// Duplicate names are refused, because SAM @SQ names must be unique.
fn load_reference(ref_path: &str) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let mut ref_names: Vec<String> = Vec::new();
    let mut ref_seqs: Vec<Vec<u8>> = Vec::new();

    iterate_fasta(ref_path, |name, seq| {
        ref_names.push(name.to_string());
        // Case-sensitive alignment scored every soft-masked base as a mismatch, hence the uppercasing.
        ref_seqs.push(
            seq.iter()
                .map(|&base| match base.to_ascii_uppercase() {
                    base_upper @ (b'A' | b'C' | b'G' | b'T') => base_upper,
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

    Ok((ref_names, ref_seqs))
}

/// Opens the SAM destination, which is stdout when no -o path is given.
fn open_output(output: Option<&str>) -> Result<Box<dyn Write>, String> {
    Ok(match output {
        Some(path) => Box::new(BufWriter::new(
            std::fs::File::create(path).map_err(|err| err.to_string())?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    })
}

/// Writes the @HD, @SQ, optional @RG and @PG header lines.
fn write_sam_header(
    out: &mut dyn Write,
    ref_names: &[String],
    ref_seqs: &[Vec<u8>],
    rg_line: Option<&str>,
    cmdline: &str,
) -> Result<(), String> {
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted").map_err(|err| err.to_string())?;

    for (i, name) in ref_names.iter().enumerate() {
        writeln!(out, "@SQ\tSN:{name}\tLN:{}", ref_seqs[i].len()).map_err(|err| err.to_string())?;
    }

    if let Some(line) = rg_line {
        writeln!(out, "{line}").map_err(|err| err.to_string())?;
    }

    // CL holds the command line, so the file records how it was made (tabs and line breaks would break the header line).
    let cmd_clean: String = cmdline
        .chars()
        .map(|ch| {
            if ch == '\t' || ch == '\n' || ch == '\r' {
                ' '
            } else {
                ch
            }
        })
        .collect();

    writeln!(
        out,
        "@PG\tID:optik\tPN:optik\tVN:{}\tCL:{cmd_clean}",
        env!("CARGO_PKG_VERSION")
    )
    .map_err(|err| err.to_string())?;

    Ok(())
}

/// Builds the seed index when the first batch arrives, so `-w 0` and automatic k can look at real reads.
/// Returns the window, the k that was chosen and the index.
#[allow(clippy::too_many_arguments)]
fn build_seed_index<'r>(
    batch: &[MapRead],
    ref_names: &[String],
    ref_seqs: &'r [Vec<u8>],
    kmer_len: usize,
    k_auto: bool,
    window: usize,
    threads: usize,
    load_start: Instant,
) -> Result<(usize, usize, SeedIdx<'r>), String> {
    let window_size = if window == 0 {
        let avg = batch.iter().map(|(_, seq, _, _)| seq.len()).sum::<usize>() / batch.len();
        let tuned = auto_tune_w(avg, kmer_len);

        eprintln!("[map] auto-tuned w = {tuned}");
        tuned
    } else {
        window
    };

    // Minimizer MPHF seed index that uses the reference as its key store.
    for (rid, seq) in ref_seqs.iter().enumerate() {
        if seq.len() >= crate::MAX_REF_SEQ_LEN {
            return Err(format!("reference {} exceeds 2^31 bases", ref_names[rid]));
        }
    }

    let mut chosen_k = kmer_len;
    let mut why = String::new();

    if k_auto {
        // Noisy long reads seed better with shorter minimizers.
        // On pbsim2 R94 ONT reads of CHM13, k = 15 placed 149,895 reads at MAPQ >= 10 with 0 wrong and k = 21 placed 147,815 with 3 wrong.
        // HiFi-like reads lose nothing at k = 21, which runs 2.3x faster.
        match noisy_long_reads(batch) {
            ReadNoise::Noisy(err_rate) => {
                chosen_k = NOISY_K;
                why = format!("base qualities: {:.1}% error", PERCENT * err_rate);
            }
            ReadNoise::Accurate(err_rate) => {
                why = format!("base qualities: {:.2}% error", PERCENT * err_rate)
            }
            ReadNoise::Short => why = "short reads".into(),
            ReadNoise::Unknown => {}
        }
    }

    let mut seeds = SeedIdx::build(
        ref_seqs,
        chosen_k,
        window_size,
        SEED_MAX_OCC,
        threads.max(1),
    )?;

    if k_auto && why.is_empty() {
        // Long reads without qualities are judged by the share of their minimizers the reference holds, about (1 - error)^k, which at k = 21 is 96% at 0.2% error and 34% at 5%.
        let presence = seed_presence(batch, &seeds, chosen_k);

        if presence < NOISY_PRESENCE {
            let first_k = chosen_k;

            chosen_k = NOISY_K;
            why = format!(
                "{:.0}% of long-read minimizers found at k = {first_k}",
                PERCENT * presence
            );
            drop(seeds); // free the first index before building the one for noisy reads
            seeds = SeedIdx::build(
                ref_seqs,
                chosen_k,
                window_size,
                SEED_MAX_OCC,
                threads.max(1),
            )?;
        } else {
            why = format!("{:.0}% of long-read minimizers found", PERCENT * presence);
        }
    }

    if k_auto {
        eprintln!("[map] k = {chosen_k} ({why})");
    }

    // The summary line also reports the repeat fraction that hidden_odds uses and the occurrence cutoff of the chaining path.
    eprintln!(
        "[map] {} minimizers indexed ({:.1} MB, repeat fraction {:.4}, chaining occurrence \
               threshold {}) in {:.2}s",
        seeds.distinct(),
        seeds.bytes() as f64 / BYTES_PER_MB,
        seeds.repeat_fraction(),
        seeds.mid_occ(),
        load_start.elapsed().as_secs_f64()
    );
    Ok((window_size, chosen_k, seeds))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run(
    ref_path: &str,
    reads_path: &str,
    reads2: Option<&str>,
    kmer_len: usize,
    k_auto: bool,
    window: usize,
    max_ed: i32,
    output: Option<&str>,
    threads: usize,
    cs_mode: CsMode,
    read_group: Option<&str>,
    cmdline: &str,
) -> Result<(), String> {
    // -R adds a header line plus an RG:Z tag on every record, as bwa and minimap2 -R do, since GATK and Picard require read groups.
    let rg_line = read_group.map(parse_read_group).transpose()?;

    if k_auto {
        eprintln!(
            "[map] loading reference {ref_path} (k chosen from the reads, default {kmer_len}; w={window})..."
        );
    } else {
        eprintln!("[map] loading reference {ref_path} (k={kmer_len} w={window})...");
    }

    let load_start = Instant::now();

    let (ref_names, ref_seqs) = load_reference(ref_path)?;

    let mut out = open_output(output)?;

    write_sam_header(&mut out, &ref_names, &ref_seqs, rg_line.as_ref().map(|(line, _)| line.as_str()), cmdline)?;
    let rg_tag: String = rg_line
        .as_ref()
        .map(|(_, rg_id)| format!("\tRG:Z:{rg_id}"))
        .unwrap_or_default();

    // A private pool keeps -t exact even when rayon's global pool is sized differently.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|err| format!("failed to build {threads}-thread pool: {err}"))?;

    // Reads are streamed in bounded batches.
    // Each batch is mapped in parallel and written in input order, which is deterministic at any thread count and bounds memory by one batch of reads and records.
    // The index is built when the first batch arrives, so `-w 0` can be tuned from real read lengths.
    const BATCH: usize = 1 << 16;
    // Name, seq, qual, and R2's seq and qual in paired mode.
    type Rd = MapRead;
    let mut batch: Vec<Rd> = Vec::with_capacity(BATCH);
    let mut model: Option<InsertModel> = None;
    let mut state: Option<(usize, usize, SeedIdx<'_>)> = None; // (w, k, index)
    let (mut n_reads, mut n_mapped) = (0u64, 0u64);
    let mut io_err: Option<String> = None;
    let tmap = Instant::now();

    let mut flush = |batch: &mut Vec<Rd>| -> Result<(), String> {
        if batch.is_empty() {
            return Ok(());
        }

        // The first batch also picks k when it is automatic, so the index waits for it.
        if state.is_none() {
            state = Some(build_seed_index(batch, &ref_names, &ref_seqs, kmer_len, k_auto, window, threads, load_start)?);
        }

        let (_, state_k, seeds) = state.as_ref().expect("seed state is built before the first batch");
        let kmer_len = *state_k;
        let ctx = AlignCtx {
            kmer_len,
            seeds,
            ref_seqs: &ref_seqs,
            cs_mode,
            // rho scales the hidden-copy prior by how repetitive this reference is.
            rho: seeds.repeat_fraction(),
        };
        let map_one = |read_name: &str, seq: &[u8]| map_single_read(read_name, seq, max_ed, &ctx, &ref_names);
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
                    .map(|(name, seq1, _, mate)| {
                        let seq2 = &mate.as_ref().expect("paired mode always reads a mate").0;

                        (
                            map_mate(name, seq1, max_ed, &ctx, &ref_names),
                            map_mate(name, seq2, max_ed, &ctx, &ref_names),
                        )
                    })
                    .collect()
            });

            // The insert model comes from the first batch's confident pairs and is reused for the rest of the run.
            if model.is_none() {
                model = Some(estimate_insert_model(&firsts, batch, &ctx));
            }

            let pair_model = model.expect("the insert model is estimated before pairing starts");
            // The second pass pairs the mates in parallel, then the pairs are flattened back into R1, R2 order.
            let pairs: Vec<(Vec<String>, Vec<String>)> = pool.install(|| {
                firsts
                    .into_par_iter()
                    .zip(batch.par_iter())
                    .map(|((aln1, aln2), (name, seq1, qual1, mate))| {
                        let (seq2, qual2) = mate.as_ref().expect("paired mode always reads a mate");
                        let (mut lines1, mut lines2) = pair_records(
                            name,
                            seq1,
                            seq2,
                            aln1,
                            aln2,
                            &pair_model,
                            max_ed,
                            &ctx,
                            &ref_names,
                        );

                        fill_qual(&mut lines1, qual1);
                        fill_qual(&mut lines2, qual2);
                        (lines1, lines2)
                    })
                    .collect()
            });

            pairs
                .into_iter()
                .flat_map(|(lines1, lines2)| [lines1, lines2])
                .collect()
        };

        for lines in recs {
            n_reads += 1;
            // Count reads whose primary record, the first line, is mapped.
            let flag: u16 = lines
                .first()
                .and_then(|line| line.split('\t').nth(1))
                .and_then(|flag_field| flag_field.parse().ok())
                .unwrap_or(SAM_FLAG_UNMAPPED);

            if flag & SAM_FLAG_UNMAPPED == 0 {
                n_mapped += 1;
            }

            for line in lines {
                writeln!(out, "{line}{rg_tag}").map_err(|err| err.to_string())?;
            }
        }

        batch.clear();
        Ok(())
    };

    let mut read_prep = ReadPrep::default();
    let mut prep = |read_name: &str, seq: &[u8], qual: &[u8]| read_prep.prepare(read_name, seq, qual);

    // In paired mode R2 is parsed on its own thread, but only while this thread waits for it.
    // This thread asks for the next R2_BLOCK records and blocks until they arrive, so at most one of the two reading threads runs at a time and the mapping pool keeps all -t threads (bench/thread_check.py).
    // Zero-capacity channels make the R2 reader a rendezvous. It reads a block only when asked and hands it over directly.
    let (want_tx, want_rx) = std::sync::mpsc::sync_channel::<()>(0);
    let (block_tx, block_rx) = std::sync::mpsc::sync_channel::<R2Block>(0);
    let r2_thread = reads2.map(|path2| spawn_r2_reader(path2, want_rx, block_tx));
    let mut r2_feed = R2Feed {
        want: want_tx,
        block_rx,
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

        // Mates pair by record order, so the two files must have the same count.
        n_rec += 1;

        if reads2.is_some() {
            let (name2, seq2, qual2) = match r2_feed.next() {
                Some(Ok(rec)) => rec,
                Some(Err(err)) => {
                    io_err = Some(err);
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

            // Mates are matched by name with any /1 or /2 suffix removed, which catches R1 and R2 files that have drifted out of step.
            if mate_name(read_name.split_ascii_whitespace().next().unwrap_or(""))
                != mate_name(name2.split_ascii_whitespace().next().unwrap_or(""))
            {
                io_err = Some(format!(
                    "mate names differ at record {n_rec}: '{read_name}' vs '{name2}' \
                                       (R1 and R2 must list the pairs in the same order)"
                ));
                return;
            }

            let (name, seq, qual) = prep(mate_name(read_name), seq, qual);
            let (_, seq2_clean, qual2_clean) = prep(&name2, &seq2, &qual2);

            batch.push((name, seq, qual, Some((seq2_clean, qual2_clean))));
        } else {
            let (name, seq, qual) = prep(read_name, seq, qual);

            batch.push((name, seq, qual, None));
        }

        if batch.len() >= BATCH {
            // The R2 block being parsed is collected first, so no reader runs beside the pool.
            if reads2.is_some()
                && let Some(err) = r2_feed.collect()
            {
                io_err = Some(err);
                return;
            }

            if let Err(err) = flush(&mut batch) {
                io_err = Some(err);
            }
        }
    })?;
    if reads2.is_some() && io_err.is_none() {
        match r2_feed.next() {
            Some(Ok((extra_name, _, _))) => {
                io_err = Some(format!(
                    "the second reads file has more records than the first (first extra: '{extra_name}')"
                ))
            }
            Some(Err(err)) => io_err = Some(err),
            None => {}
        }
    }

    let _ = r2_feed.collect(); // a block still being parsed

    drop(r2_feed); // the R2 reader stops at its next wait

    if let Some(handle) = r2_thread {
        let _ = handle.join();
    }

    if let Some(err) = io_err {
        return Err(err);
    }

    flush(&mut batch)?;
    // Ends the closure's borrows of the writer and counters, since it owns nothing to free.
    #[allow(clippy::drop_non_drop)]
    drop(flush);
    out.flush().map_err(|err| err.to_string())?;

    #[allow(clippy::drop_non_drop)]
    drop(prep); // ends its borrows, as above
    let pct = PERCENT * n_mapped as f64 / n_reads.max(1) as f64;

    eprintln!(
        "[map] {n_reads} reads in {:.2}s -- {n_mapped} mapped ({pct:.1}%)",
        tmap.elapsed().as_secs_f64()
    );
    let n_bad_qual = read_prep.bad_qual_records;

    if n_bad_qual > 0 {
        eprintln!("[map] warning: {n_bad_qual} FASTQ records had unusable qualities (QUAL '*')");
    }

    Ok(())
}

#[cfg(test)]
mod tests;
