//! Looks at a sample of the reads before mapping starts.
//! It decides whether long reads are noisy enough to need shorter seeds, and it picks the minimizer window.

use super::*;

/// k chosen for noisy long reads, the same 15 as minimap2's map-ont.
pub(crate) const NOISY_K: usize = 15;

/// Long reads count as noisy above this per-base error, taken from base qualities, or below this share of their k = 21 minimizers found in the reference.
/// Both thresholds mark about 2% error, since (1 - 0.02)^21 = 0.65.
pub(crate) const NOISY_ERROR: f64 = 0.02;

pub(super) const NOISY_PRESENCE: f64 = 0.65;

/// Reads sampled for the choice.
pub(super) const NOISE_SAMPLE: usize = 1000;

/// What the first batch says about the reads' error rate.
/// They are noisy or accurate long reads, with the mean error from base qualities, mostly short reads, or long reads without qualities.
pub(super) enum ReadNoise {
    Noisy(f64),
    Accurate(f64),
    Short,
    Unknown,
}

/// A read of a mapping batch, as name, seq, qual and, when paired, R2's seq and qual.
pub(super) type MapRead = (String, Vec<u8>, Vec<u8>, Option<(Vec<u8>, Vec<u8>)>);

/// Whether the first batch is long reads and, if they have base qualities, how noisy they are (mean per-base error 10^(-Q/10) over sampled long reads).
pub(super) fn noisy_long_reads(batch: &[MapRead]) -> ReadNoise {
    // Long reads must make up at least half the batch.
    const LONG_SHARE_DEN: usize = 2;
    let n_long = batch
        .iter()
        .filter(|read| read.1.len() >= CHAIN_MIN_LEN)
        .count();

    if n_long * LONG_SHARE_DEN < batch.len() {
        return ReadNoise::Short; // mostly short reads keep the fitted short-read k
    }

    let long: Vec<_> = batch
        .iter()
        .filter(|read| read.1.len() >= CHAIN_MIN_LEN)
        .take(NOISE_SAMPLE)
        .collect();
    let (mut err_sum, mut base_count) = (0.0f64, 0u64);

    for read in &long {
        if read.2.len() != read.1.len() {
            return ReadNoise::Unknown;
        } // FASTA has no qualities

        for &qual in &read.2 {
            err_sum +=
                PHRED_SCALE.powf(-((qual.saturating_sub(PHRED_ASCII_OFFSET)) as f64) / PHRED_SCALE);
            base_count += 1;
        }
    }

    if base_count == 0 {
        return ReadNoise::Unknown;
    }

    let err = err_sum / base_count as f64;

    if err > NOISY_ERROR {
        ReadNoise::Noisy(err)
    } else {
        ReadNoise::Accurate(err)
    }
}

/// Share of the sampled long reads' minimizers found in the index (or dropped there as too repetitive), about the share of error-free k-mers.
pub(super) fn seed_presence(batch: &[MapRead], seeds: &SeedIdx<'_>, kmer_len: usize) -> f64 {
    let (mut found, mut total) = (0u64, 0u64);

    for read in batch
        .iter()
        .filter(|read| read.1.len() >= CHAIN_MIN_LEN)
        .take(NOISE_SAMPLE)
    {
        let (hits, n_min, n_rep) = seeds.collect_with_rep(&read.1);
        let len = read.1.len();
        let mut pos: Vec<usize> = hits
            .iter()
            .map(|hit| {
                // Reverse hits report positions on the reverse complement. Mapping them back lets one read position count once across strands.
                if hit.rev != 0 {
                    len - kmer_len - hit.qpos as usize
                } else {
                    hit.qpos as usize
                }
            })
            .collect();

        pos.sort_unstable();
        pos.dedup();
        // Over-occurring minimizers were present but filtered, so they count as found.
        found += pos.len() as u64 + n_rep as u64;
        total += n_min as u64;
    }

    // No sampled minimizers gives no evidence of noise.
    if total == 0 {
        1.0
    } else {
        found as f64 / total as f64
    }
}

/// The minimizer window for `-w 0`, read_len / 4 bounded below by 10 and above by max(k / 2, 10).
pub(super) fn auto_tune_w(read_len: usize, kmer_len: usize) -> usize {
    /// Smallest window the tuning picks.
    const MIN_AUTO_W: usize = 10;
    /// The window is this fraction (1/4) of the read length.
    const READ_LEN_PER_W: usize = 4;
    /// The window is at most half the minimizer length.
    const K_PER_MAX_W: usize = 2;
    let max_w = (kmer_len / K_PER_MAX_W).max(MIN_AUTO_W);

    (read_len / READ_LEN_PER_W).min(max_w).max(MIN_AUTO_W)
}
