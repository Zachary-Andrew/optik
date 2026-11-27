//! `optik query` looks up every k-mer of the query records in a saved SPSS or FM index and prints SSHash's query report.
//!
//! A k-mer window containing a non-ACGT character is invalid.
//! The C++ scanners restart at such characters, so records are passed whole.
//!
//! num_searches and num_extensions follow SSHash's streaming query, where a search is a full index lookup and an extension answers a k-mer from the previous one without a lookup.
//! Only the SPSS index streams, and the FM index reports num_searches equal to the valid k-mers.
//!
//! Each record, or each batch of records for SPSS, is scanned by one C++ call, i.e. one FFI crossing and no per-k-mer callback.
//! With -t > 1 the records are spread over a rayon thread pool.
//! Long records, such as a genome or an assembly, are first cut into window-exact slices.
//! A slice holds the k-mer windows [s, e) as bytes [s, e + k - 1), so every window belongs to exactly one slice and num_kmers, positive, negative and invalid are exact at any thread count.
//! num_searches and num_extensions can differ slightly with -t > 1, because each slice restarts the streaming cursor.

use std::time::Instant;

use rayon::prelude::*;

use crate::seqio::{batch_bytes_for, for_each_batch};
use crate::{FM_MAGIC, FmIdx, SPSS_MAGIC, SpssIdx, index_magic};

/// A rayon task gets at least this many query records when there are many records per thread.
const MIN_RECORDS_PER_TASK: usize = 256;
/// Multiplier from a fraction to a percentage in the report.
const PERCENT: f64 = 100.0;
/// Time units of the report's elapsed line.
const MS_PER_SEC: u64 = 1000;
const SECS_PER_MIN: u64 = 60;
const NS_PER_MS: f64 = 1e6;
/// Characters below U+0020 are JSON control characters and must be written as \u escapes.
const JSON_FIRST_PRINTABLE: u32 = 0x20;

/// A loaded SPSS or FM index.
enum LoadedIdx {
    Fm(FmIdx),
    Spss(SpssIdx),
}

impl LoadedIdx {
    /// Scans one record, looking up every valid k-mer, and returns the counters (valid, positive, searches, extensions) it contributed.
    fn scan_chunk(&self, chunk: &[u8], kmer_len: usize) -> (u64, u64, u64, u64) {
        match self {
            LoadedIdx::Fm(fm) => fm.scan_chunk(chunk, kmer_len),
            LoadedIdx::Spss(spss) => spss.scan_chunk(chunk),
        }
    }
}

/// Query counters of one record or batch, additive across records.
#[derive(Default, Clone, Copy)]
struct Rec {
    total_positions: u64,
    num_valid: u64,
    num_positive: u64,
    num_searches: u64,
    num_extensions: u64,
}

/// Counts one query record, with record_len - k + 1 total positions and one scan call (the scanner excludes non-ACGT windows).
fn process_record(seq: &[u8], kmer_len: usize, idx: &LoadedIdx) -> Rec {
    let total_positions = if seq.len() >= kmer_len {
        (seq.len() - kmer_len + 1) as u64
    } else {
        0
    };
    let (valid, positive, searches, extensions) = if total_positions > 0 {
        idx.scan_chunk(seq, kmer_len)
    } else {
        (0, 0, 0, 0)
    };

    Rec {
        total_positions,
        num_valid: valid,
        num_positive: positive,
        num_searches: searches,
        num_extensions: extensions,
    }
}

/// Runs `optik query`, looking up the k-mers of `queries` in the index at `index_path` and printing SSHash's query report on stdout, with a JSON summary on stderr.
///
/// The index format comes from the file's magic, either an SPSS index (OKSPSS) or a k-order BWT FM index (OKFM).
///
/// * `k_arg` is the `-k` given on the command line, if any.
///   Every index records its k, which applies when `-k` is absent.
/// * `threads` of 1 runs single-threaded, while more spread the records over a pool.
/// * `ids_out`, with an SPSS index, names a file that also receives each record's k-mer ids, one line per record (`--ids`).
/// * `verify` turns on the SPSS index checksum check at load.
///
/// # Errors
///
/// Returns an error for an unreadable or unsupported index, a `-k` that differs from the index's k (it would report every k-mer negative), `--ids` with an FM index, or an I/O failure.
pub(crate) fn query(
    index_path: &str,
    queries: &str,
    k_arg: Option<usize>,
    threads: usize,
    ids_out: Option<&str>,
    verify: bool,
) -> Result<(), String> {
    let k_check = |built: usize| -> Result<(), String> {
        match k_arg {
            Some(given_k) if given_k != built => Err(format!(
                "index was built with k={built} but -k {given_k} was given"
            )),
            _ => Ok(()),
        }
    };
    // FM_MAGIC or SPSS_MAGIC at the start of the file selects the loader, which is why query has no index-type flag.
    let idx: LoadedIdx = match index_magic(index_path)? {
        Some(FM_MAGIC) => {
            let fm = FmIdx::load(index_path)?;

            k_check(fm.kmer_len())?;
            LoadedIdx::Fm(fm)
        }
        Some(SPSS_MAGIC) => {
            let spss = SpssIdx::load(index_path, verify)?;

            k_check(spss.kmer_len())?;
            LoadedIdx::Spss(spss)
        }
        _ => {
            return Err(format!(
                "{index_path} is not an optik SPSS or FM index (bad magic); rebuild it with \
                 build-spss or build-fm"
            ));
        }
    };

    let kmer_len = match &idx {
        LoadedIdx::Fm(fm) => fm.kmer_len(),
        LoadedIdx::Spss(spss) => spss.kmer_len(),
    };

    // One pool serves the whole run, reused by every batch.
    let pool = if threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|err| err.to_string())?,
        )
    } else {
        None
    };

    let t_start = Instant::now();

    let add = |lhs: Rec, rhs: Rec| Rec {
        total_positions: lhs.total_positions + rhs.total_positions,
        num_valid: lhs.num_valid + rhs.num_valid,
        num_positive: lhs.num_positive + rhs.num_positive,
        num_searches: lhs.num_searches + rhs.num_searches,
        num_extensions: lhs.num_extensions + rhs.num_extensions,
    };
    // SPSS records go through the pipelined batch scan (misses of consecutive lookups overlap, across record boundaries too).
    let batch_rec = |recs: &[&[u8]], spss: &SpssIdx| -> Rec {
        let total_positions = recs
            .iter()
            .map(|record| {
                if record.len() >= kmer_len {
                    (record.len() - kmer_len + 1) as u64
                } else {
                    0
                }
            })
            .sum();
        let (valid, positive, searches, extensions) = spss.scan_many(recs);

        Rec {
            total_positions,
            num_valid: valid,
            num_positive: positive,
            num_searches: searches,
            num_extensions: extensions,
        }
    };

    const BATCH: usize = 4096; // records per SPSS scan_many() call

    // --ids writes one line per record, batch by batch, so memory stays bounded by one batch.
    // The ids and the counts come from the same pass over the queries, since a second pass would read a stdin input ('-q -') as empty and report 0 k-mers.
    use std::io::Write;
    let mut ids_w: Option<(std::io::BufWriter<std::fs::File>, &SpssIdx)> = match ids_out {
        None => None,
        Some(path) => {
            let LoadedIdx::Spss(spss) = &idx else {
                return Err("--ids needs an SPSS index (build-spss)".into());
            };

            Some((
                std::io::BufWriter::new(
                    std::fs::File::create(path).map_err(|err| format!("{path}: {err}"))?,
                ),
                spss,
            ))
        }
    };
    let id_line = |spss: &SpssIdx, record: &&[u8]| -> String {
        spss.ids(record)
            .iter()
            .map(|id| id.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    };

    // Records stream through seqio::for_each_batch in batches sized by the thread count, and a batch is split into work units processed in parallel.
    // A record, or a window-exact slice of one, goes to exactly one worker, so the aggregated Rec is exact.
    let mut out = Rec::default();

    for_each_batch(queries, batch_bytes_for(threads), |records| {
        if let Some((ids_writer, spss_ref)) = ids_w.as_mut() {
            let spss: &SpssIdx = spss_ref;
            let lines: Vec<String> = match &pool {
                Some(thread_pool) => thread_pool.install(|| {
                    records
                        .par_iter()
                        .map(|record| id_line(spss, record))
                        .collect()
                }),
                None => records.iter().map(|record| id_line(spss, record)).collect(),
            };

            for line in lines {
                writeln!(ids_writer, "{line}").map_err(|err| err.to_string())?;
            }
        }

        // Work units are whole records, except that with -t > 1 a record with more windows than one slice is cut into window-exact slices, since a single multi-Mbp record would otherwise run on one worker.
        let total_windows: usize = records
            .iter()
            .map(|record| record.len().saturating_sub(kmer_len - 1))
            .sum();

        // Each thread gets about 8 slices of at least MIN_SLICE_WINDOWS windows.
        const MIN_SLICE_WINDOWS: usize = 1 << 16;
        const SLICES_PER_THREAD: usize = 8;
        let slice_windows =
            (total_windows / (threads.max(1) * SLICES_PER_THREAD)).max(MIN_SLICE_WINDOWS);
        let mut work: Vec<&[u8]> = Vec::with_capacity(records.len());

        for &record in records {
            let nwin = record.len().saturating_sub(kmer_len - 1);

            if pool.is_none() || nwin <= slice_windows {
                work.push(record);
                continue;
            }

            let mut slice_start = 0;

            while slice_start < nwin {
                let slice_end = (slice_start + slice_windows).min(nwin);

                work.push(&record[slice_start..slice_end + kmer_len - 1]);
                slice_start = slice_end;
            }
        }

        let rec = match (&idx, &pool) {
            (LoadedIdx::Spss(spss), None) => work
                .chunks(BATCH)
                .map(|chunk| batch_rec(chunk, spss))
                .fold(Rec::default(), add),
            (LoadedIdx::Spss(spss), Some(thread_pool)) => {
                // Sliced records go out in small chunks, so their slices spread over the threads.
                let batch = if work.len() < threads * BATCH {
                    1
                } else {
                    BATCH
                };

                thread_pool.install(|| {
                    work.par_chunks(batch)
                        .map(|chunk| batch_rec(chunk, spss))
                        .reduce(Rec::default, add)
                })
            }
            // The single-threaded headline path has no pool and no task overhead.
            (_, None) => work
                .iter()
                .map(|record| process_record(record, kmer_len, &idx))
                .fold(Rec::default(), add),
            (_, Some(thread_pool)) => thread_pool.install(|| {
                work.par_iter()
                    // Batches smaller than threads * MIN_RECORDS_PER_TASK split per record so every thread gets work, and larger ones use coarser tasks.
                    .with_min_len(if work.len() < threads * MIN_RECORDS_PER_TASK {
                        1
                    } else {
                        MIN_RECORDS_PER_TASK
                    })
                    .map(|record| process_record(record, kmer_len, &idx))
                    .reduce(Rec::default, add)
            }),
        };

        out = add(out, rec);
        Ok(())
    })?;
    if let Some((mut ids_writer, _)) = ids_w {
        ids_writer.flush().map_err(|err| err.to_string())?;
    }

    let total_positions = out.total_positions;
    let num_valid = out.num_valid;
    let num_positive = out.num_positive;
    let num_searches = out.num_searches;
    let num_extensions = out.num_extensions;
    let num_invalid = total_positions - num_valid;
    let num_negative = num_valid - num_positive;

    let elapsed_ms = t_start.elapsed().as_millis() as u64;

    let pct = |part: u64, whole: u64| -> f64 {
        if whole > 0 {
            (part as f64 * PERCENT) / whole as f64
        } else {
            0.0
        }
    };

    // Machine-readable JSON summary on stderr.
    // Paths are JSON-escaped, since a quote or backslash in a file name would otherwise make the summary unparseable for the benchmark harnesses.
    eprintln!(
        "{{\"index_filename\":{},\"query_filename\":{},\"num_kmers\":{},\"num_positive_kmers\":{},\"num_negative_kmers\":{},\"num_invalid_kmers\":{},\"num_searches\":{},\"elapsed_millisec\":{}}}",
        json_str(index_path),
        json_str(queries),
        total_positions,
        num_positive,
        num_negative,
        num_invalid,
        num_searches,
        elapsed_ms
    );

    // Human-readable report on stdout, in SSHash's format.
    println!("==== query report:");
    println!("num_kmers = {total_positions}");
    println!(
        "num_positive_kmers = {num_positive} ({:.4}%)",
        pct(num_positive, total_positions)
    );
    println!(
        "num_negative_kmers = {num_negative} ({:.4}%)",
        pct(num_negative, total_positions)
    );
    println!(
        "num_invalid_kmers = {num_invalid} ({:.4}%)",
        pct(num_invalid, total_positions)
    );
    println!(
        "num_searches = {num_searches}/{num_valid} ({:.4}%)",
        pct(num_searches, num_valid)
    );
    println!(
        "num_extensions = {num_extensions}/{num_valid} ({:.4}%)",
        pct(num_extensions, num_valid)
    );

    let sec = elapsed_ms / MS_PER_SEC;
    let min = sec / SECS_PER_MIN;
    let ns_per_kmer = if total_positions > 0 {
        (elapsed_ms as f64 * NS_PER_MS) / total_positions as f64
    } else {
        0.0
    };

    println!("elapsed = {sec} sec / {min} min / {ns_per_kmer:.0} ns/kmer");

    Ok(())
}

/// A JSON string literal with quotes, backslashes and control bytes escaped.
fn json_str(raw: &str) -> String {
    const QUOTES: usize = 2;
    let mut escaped = String::with_capacity(raw.len() + QUOTES);

    escaped.push('"');

    for ch in raw.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            // JSON forbids raw control characters inside strings.
            ch if (ch as u32) < JSON_FIRST_PRINTABLE => {
                escaped.push_str(&format!("\\u{:04x}", ch as u32))
            }
            ch => escaped.push(ch),
        }
    }

    escaped.push('"');
    escaped
}

#[cfg(test)]
mod tests {
    use super::json_str;

    #[test]
    fn json_strings_escape_quotes_backslashes_and_control_bytes() {
        assert_eq!(json_str("plain"), "\"plain\"");
        assert_eq!(json_str("a\"b\\c"), "\"a\\\"b\\\\c\"");
        assert_eq!(json_str("tab\there"), "\"tab\\u0009here\"");
    }
}
