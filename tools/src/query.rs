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

/// A loaded SPSS or FM index.
enum LoadedIdx {
    Fm(FmIdx),
    Spss(SpssIdx),
}

impl LoadedIdx {
    /// Scans one record, looking up every valid k-mer, and returns the counters (valid, positive, searches, extensions) it contributed.
    fn scan_chunk(&self, chunk: &[u8], k: usize) -> (u64, u64, u64, u64) {
        match self {
            LoadedIdx::Fm(i) => i.scan_chunk(chunk, k),
            LoadedIdx::Spss(i) => i.scan_chunk(chunk),
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
fn process_record(seq: &[u8], k: usize, idx: &LoadedIdx) -> Rec {
    let total_positions = if seq.len() >= k {
        (seq.len() - k + 1) as u64
    } else {
        0
    };
    let (v, p, s, e) = if total_positions > 0 {
        idx.scan_chunk(seq, k)
    } else {
        (0, 0, 0, 0)
    };
    Rec {
        total_positions,
        num_valid: v,
        num_positive: p,
        num_searches: s,
        num_extensions: e,
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
            Some(k) if k != built => Err(format!(
                "index was built with k={built} but -k {k} was given"
            )),
            _ => Ok(()),
        }
    };
    let idx: LoadedIdx = match index_magic(index_path)? {
        Some(FM_MAGIC) => {
            let i = FmIdx::load(index_path)?;
            k_check(i.k())?;
            LoadedIdx::Fm(i)
        }
        Some(SPSS_MAGIC) => {
            let i = SpssIdx::load(index_path, verify)?;
            k_check(i.k())?;
            LoadedIdx::Spss(i)
        }
        _ => {
            return Err(format!(
                "{index_path} is not an optik SPSS or FM index (bad magic); rebuild it with \
                 build-spss or build-fm"
            ));
        }
    };

    let k = match &idx {
        LoadedIdx::Fm(i) => i.k(),
        LoadedIdx::Spss(i) => i.k(),
    };

    // One pool serves the whole run, reused by every batch.
    let pool = if threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|e| e.to_string())?,
        )
    } else {
        None
    };

    let t_start = Instant::now();

    let add = |a: Rec, b: Rec| Rec {
        total_positions: a.total_positions + b.total_positions,
        num_valid: a.num_valid + b.num_valid,
        num_positive: a.num_positive + b.num_positive,
        num_searches: a.num_searches + b.num_searches,
        num_extensions: a.num_extensions + b.num_extensions,
    };
    // SPSS records go through the pipelined batch scan (misses of consecutive lookups overlap, across record boundaries too).
    let batch_rec = |recs: &[&[u8]], si: &SpssIdx| -> Rec {
        let total_positions = recs
            .iter()
            .map(|r| {
                if r.len() >= k {
                    (r.len() - k + 1) as u64
                } else {
                    0
                }
            })
            .sum();
        let (v, p, s, e) = si.scan_many(recs);
        Rec {
            total_positions,
            num_valid: v,
            num_positive: p,
            num_searches: s,
            num_extensions: e,
        }
    };
    const BATCH: usize = 4096; // records per SPSS scan_many() call

    // --ids writes one line per record, batch by batch, so memory stays bounded by one batch.
    // The ids and the counts come from the same pass over the queries, since a second pass would read a stdin input ('-q -') as empty and report 0 k-mers.
    use std::io::Write;
    let mut ids_w: Option<(std::io::BufWriter<std::fs::File>, &SpssIdx)> = match ids_out {
        None => None,
        Some(path) => {
            let LoadedIdx::Spss(si) = &idx else {
                return Err("--ids needs an SPSS index (build-spss)".into());
            };
            Some((
                std::io::BufWriter::new(
                    std::fs::File::create(path).map_err(|e| format!("{path}: {e}"))?,
                ),
                si,
            ))
        }
    };
    let id_line = |si: &SpssIdx, r: &&[u8]| -> String {
        si.ids(r)
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    };

    // Records stream through seqio::for_each_batch in batches sized by the thread count, and a batch is split into work units processed in parallel.
    // A record, or a window-exact slice of one, goes to exactly one worker, so the aggregated Rec is exact.
    let mut out = Rec::default();
    for_each_batch(queries, batch_bytes_for(threads), |records| {
        if let Some((w, si)) = ids_w.as_mut() {
            let si: &SpssIdx = si;
            let lines: Vec<String> = match &pool {
                Some(pl) => pl.install(|| records.par_iter().map(|r| id_line(si, r)).collect()),
                None => records.iter().map(|r| id_line(si, r)).collect(),
            };
            for l in lines {
                writeln!(w, "{l}").map_err(|e| e.to_string())?;
            }
        }
        // Work units are whole records, except that with -t > 1 a record with more windows than one slice is cut into window-exact slices, since a single multi-Mbp record would otherwise run on one worker.
        let total_windows: usize = records.iter().map(|r| r.len().saturating_sub(k - 1)).sum();
        // about 8 slices per thread, each at least MIN_SLICE_WINDOWS
        const MIN_SLICE_WINDOWS: usize = 1 << 16;
        let slice_windows = (total_windows / (threads.max(1) * 8)).max(MIN_SLICE_WINDOWS);
        let mut work: Vec<&[u8]> = Vec::with_capacity(records.len());
        for &r in records {
            let nwin = r.len().saturating_sub(k - 1);
            if pool.is_none() || nwin <= slice_windows {
                work.push(r);
                continue;
            }
            let mut s = 0;
            while s < nwin {
                let e = (s + slice_windows).min(nwin);
                work.push(&r[s..e + k - 1]);
                s = e;
            }
        }
        let rec = match (&idx, &pool) {
            (LoadedIdx::Spss(si), None) => work
                .chunks(BATCH)
                .map(|c| batch_rec(c, si))
                .fold(Rec::default(), add),
            (LoadedIdx::Spss(si), Some(pl)) => {
                // small chunks when records were sliced, so slices spread out
                let batch = if work.len() < threads * BATCH {
                    1
                } else {
                    BATCH
                };
                pl.install(|| {
                    work.par_chunks(batch)
                        .map(|c| batch_rec(c, si))
                        .reduce(Rec::default, add)
                })
            }
            // The single-threaded headline path has no pool and no task overhead.
            (_, None) => work
                .iter()
                .map(|r| process_record(r, k, &idx))
                .fold(Rec::default(), add),
            (_, Some(pl)) => pl.install(|| {
                work.par_iter()
                    .with_min_len(if work.len() < threads * 256 { 1 } else { 256 })
                    .map(|r| process_record(r, k, &idx))
                    .reduce(Rec::default, add)
            }),
        };
        out = add(out, rec);
        Ok(())
    })?;
    if let Some((mut w, _)) = ids_w {
        w.flush().map_err(|e| e.to_string())?;
    }

    let total_positions = out.total_positions;
    let num_valid = out.num_valid;
    let num_positive = out.num_positive;
    let num_searches = out.num_searches;
    let num_extensions = out.num_extensions;
    let num_invalid = total_positions - num_valid;
    let num_negative = num_valid - num_positive;

    let elapsed_ms = t_start.elapsed().as_millis() as u64;

    let pct = |a: u64, b: u64| -> f64 {
        if b > 0 {
            (a as f64 * 100.0) / b as f64
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

    // Human-readable report on stdout matching sshash format.
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

    let sec = elapsed_ms / 1000;
    let min = sec / 60;
    let ns_per_kmer = if total_positions > 0 {
        (elapsed_ms as f64 * 1e6) / total_positions as f64
    } else {
        0.0
    };
    println!("elapsed = {sec} sec / {min} min / {ns_per_kmer:.0} ns/kmer");

    Ok(())
}

/// A JSON string literal with quotes, backslashes and control bytes escaped.
fn json_str(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}
