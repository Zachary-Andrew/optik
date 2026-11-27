//! Index builders for SPSS, with `build-spss`, its alias `build` and `spss-update`, and for FM with `build-fm`.
//! Each prints tab-separated metrics on stdout, which the benchmark scripts parse.

use crate::{DynSpssIdx, FmIdx, RssWatermark, SpssIdx};
use std::time::Instant;

const BYTES_PER_MIB: f64 = (1 << 20) as f64;
/// Bits per byte, for index bits per k-mer.
const BITS_PER_BYTE: f64 = 8.0;

/// Builds and saves an FM membership index, the k-order BWT over the unitigs of include/fm_index.hpp.
///
/// The index is exact and its size is the BWT's rank support only.
/// Prints build_s, rss_mb, the file size idx_mb, bits_per, distinct_kmers, unitigs, super_len and bwt_len.
///
/// # Errors
///
/// Returns an error if the build or the write fails.
pub(crate) fn build_fm(
    input: &str,
    kmer_len: usize,
    min_count: u64,
    nthreads: u32,
    output: &str,
) -> Result<(), String> {
    let mut hwm = RssWatermark::new();
    let build_start = Instant::now();
    let idx = FmIdx::build_from_fasta(input, kmer_len, min_count, nthreads)
        .map_err(|err| format!("build_from_fasta failed: {err}"))?;
    let build_s = build_start.elapsed().as_secs_f64();

    // RSS is sampled before save() to leave the index file's write buffers out of the figure.
    hwm.sample();
    let rss_mb = hwm.net_mb();

    idx.save(output)?;

    // idx_mb reports the saved file's size and falls back to idx.bytes() only if the metadata call fails.
    let idx_mb = std::fs::metadata(output)
        .map(|meta| meta.len() as f64 / BYTES_PER_MIB)
        .unwrap_or(idx.bytes() as f64 / BYTES_PER_MIB);

    println!("build_s\t{build_s:.6}");
    println!("rss_mb\t{rss_mb:.3}");
    println!("idx_mb\t{idx_mb:.6}");
    println!("bits_per\t{:.2}", idx.bits_per_entry());
    println!("distinct_kmers\t{}", idx.size());
    println!("unitigs\t{}", idx.num_unitigs());
    println!("super_len\t{}", idx.super_len());
    println!("bwt_len\t{}", idx.bwt_len());

    Ok(())
}

/// Builds and saves a static SPSS index (include/spss_core.hpp).
///
/// `m_arg` 0 lets SpssCore::choose_m pick the minimizer length from the input.
/// `compact` counts the k-mers seen at least `min_count` times and indexes their maximal unitigs instead of the input records.
/// Prints build_s, rss_mb, the file size idx_mb, bits_per distinct k-mer, minimizer_len, distinct_kmers, k, m, superkmers and kmer_positions.
///
/// # Errors
///
/// Returns an error when the input has no k-mer of length `kmer_len`, or if the build or the write fails.
pub(crate) fn build_spss(
    input: &str,
    kmer_len: usize,
    m_arg: usize,
    compact: bool,
    min_count: u64,
    threads: usize,
    output: &str,
) -> Result<(), String> {
    let minimizer_len = m_arg;
    let mut hwm = RssWatermark::new();
    let build_start = Instant::now();
    let idx = SpssIdx::build_from_fasta(
        input,
        kmer_len,
        minimizer_len,
        compact,
        min_count.max(1),
        threads,
    )?;

    if idx.kmers() == 0 {
        return Err(format!(
            "no k-mers of length {kmer_len} in {} (empty input, or every record shorter than k)",
            input.replace('\n', ", ")
        ));
    }

    let build_s = build_start.elapsed().as_secs_f64();

    hwm.sample();
    idx.save(output)?;
    let file_bytes = std::fs::metadata(output)
        .map(|meta| meta.len())
        .unwrap_or(idx.bytes());
    let positions = idx.kmers();
    // Without --compact the records are indexed as given, so a k-mer that occurs several times in the input is stored several times.
    // A compact index holds each k-mer once, so its positions are already distinct.
    let distinct_kmers = if compact { positions } else { idx.distinct() };

    if distinct_kmers < positions {
        eprintln!(
            "[build-spss] warning: {} of {positions} indexed k-mer positions repeat a k-mer \
                   ({distinct_kmers} distinct); query --ids gives each distinct k-mer one id, but the ids \
                   are positions in [0, {positions}), not a bijection onto [0, {distinct_kmers}).  Use --compact for a duplicate-free index.",
            positions - distinct_kmers
        );
    }

    println!("build_s\t{build_s:.6}");
    println!("rss_mb\t{:.3}", hwm.net_mb());
    println!("idx_mb\t{:.6}", file_bytes as f64 / BYTES_PER_MIB);
    println!(
        "bits_per\t{:.2}",
        if distinct_kmers > 0 {
            file_bytes as f64 * BITS_PER_BYTE / distinct_kmers as f64
        } else {
            0.0
        }
    );
    println!("minimizer_len\t{}", idx.minimizer_len());
    println!("distinct_kmers\t{distinct_kmers}");
    println!("k\t{}", idx.kmer_len());
    println!("m\t{}", idx.minimizer_len());
    println!("superkmers\t{}", idx.superkmers());
    println!("kmer_positions\t{positions}");
    Ok(())
}

/// Applies the k-mer insertions in `adds` first and the removals in `removes` to an SPSS index and writes the compacted result to `output`.
///
/// Edits go to concurrent open-addressing sets over the static core, as include/dyn_spss.hpp describes, and compaction then rebuilds the maximal unitigs of the edited set.
/// Prints the k-mer counts before and after, the edit counts and the timings.
///
/// # Errors
///
/// Returns an error if the index or an edit file cannot be read, or the write fails.
pub(crate) fn update_spss(
    input: &str,
    adds: &[&str],
    removes: &[&str],
    threads: usize,
    output: &str,
) -> Result<(), String> {
    let start = Instant::now();
    let mut idx = DynSpssIdx::open(input)?;
    let before = idx.kmers();
    let mut inserted = 0u64;
    let mut erased = 0u64;

    // Adds apply before removes, so a k-mer named in both ends up removed.
    for add_path in adds {
        inserted += idx.apply_fasta(add_path, true)?;
    }

    for remove_path in removes {
        erased += idx.apply_fasta(remove_path, false)?;
    }

    let (added, deleted) = (idx.added(), idx.deleted());
    let t_edit = start.elapsed().as_secs_f64();

    idx.compact_save(output, threads)?;
    println!("kmers_before\t{before}");
    println!("inserted\t{inserted}");
    println!("erased\t{erased}");
    println!("delta_added\t{added}");
    println!("delta_deleted\t{deleted}");
    println!("kmers_after\t{}", idx.kmers());
    println!("edit_s\t{t_edit:.4}");
    println!("total_s\t{:.4}", start.elapsed().as_secs_f64());
    Ok(())
}
