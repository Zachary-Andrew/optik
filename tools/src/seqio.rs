//! Streaming FASTA/FASTQ record reader for `optik query`.
//!
//! Records are delivered in bounded batches, so memory does not grow with the query file.
//! A plain file is memory-mapped.
//! Records whose sequence sits on one line are slices into the mapping without a copy, while multi-line FASTA sequences and FASTQ sequences split over several lines are concatenated into an owned buffer.
//! After each batch the consumed part of the mapping is released (MADV_DONTNEED), so resident memory stays near one batch instead of the whole file plus a slice per record.
//! Gzipped input and `-` (stdin) go through the C++ reader (include/fasta_reader.hpp), batched the same way.
use std::borrow::Cow;

use crate::iterate_fasta_all;

/// Sequence bytes per batch, large enough to keep every worker busy (the query slices long records itself) and small enough to bound resident memory.
pub(crate) const BATCH_BYTES: usize = 8 << 20;

/// Batch bytes per `query` worker, so that [`batch_bytes_for`] keeps 1 MiB of input resident on one thread and [`BATCH_BYTES`] on eight or more.
pub(crate) const BATCH_BYTES_PER_THREAD: usize = 1 << 20;

/// The batch size for `threads` workers, [`BATCH_BYTES_PER_THREAD`] each and at most [`BATCH_BYTES`].
pub(crate) fn batch_bytes_for(threads: usize) -> usize {
    (BATCH_BYTES_PER_THREAD * threads.max(1)).min(BATCH_BYTES)
}

/// Records per batch, which bounds the per-record bookkeeping when records are tiny (a million single k-mers would otherwise put about 250k records, and their slice vectors, in one batch).
pub(crate) const BATCH_RECORDS: usize = 1 << 15;

/// Page granularity for releasing consumed parts of the mapping.
const PAGE: usize = 4096;

/// The two bytes that start every gzip member (RFC 1952).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Calls `on_batch` with consecutive batches of record sequences, in file order.
///
/// Every record of the file is passed exactly once, empty records included.
/// A batch closes at `batch_bytes` of sequence or [`BATCH_RECORDS`] records.
///
/// # Errors
///
/// Returns an error if the file cannot be read, and stops at the first error `on_batch` returns.
pub(crate) fn for_each_batch<F>(
    path: &str,
    batch_bytes: usize,
    mut on_batch: F,
) -> Result<(), String>
where
    F: FnMut(&[&[u8]]) -> Result<(), String>,
{
    if path != "-" {
        let file = std::fs::File::open(path).map_err(|err| format!("{path}: {err}"))?;
        let len = file.metadata().map_err(|err| err.to_string())?.len();

        if len == 0 {
            return Ok(());
        }

        // SAFETY: the file is opened read-only and not modified while mapped.
        // The standard mmap caveat applies, since external truncation would fault.
        let mapping = unsafe { memmap2::Mmap::map(&file) }.map_err(|err| err.to_string())?;

        if !mapping.starts_with(&GZIP_MAGIC) {
            return mapped_batches(&mapping, batch_bytes, &mut on_batch);
        }
    }

    // For gzip or stdin the C++ reader decodes, and records are copied and batched.
    // Empty records are kept, as on the mapped path, because `query --ids` writes one line per record and dropping one would shift every later line onto the wrong read.
    let mut batch: Vec<Vec<u8>> = Vec::new();
    let mut bytes = 0usize;
    let mut err: Option<String> = None;

    // iterate_fasta_all has no way to stop the C++ reader early, so after a failed batch the remaining records are ignored.
    iterate_fasta_all(path, |_name, seq| {
        if err.is_some() {
            return;
        }

        bytes += seq.len();
        batch.push(seq.to_vec());

        // BATCH_RECORDS caps the record count as well as the bytes, so millions of tiny records still flush in bounded batches.
        if bytes >= batch_bytes || batch.len() >= BATCH_RECORDS {
            let seqs: Vec<&[u8]> = batch.iter().map(|record| record.as_slice()).collect();

            if let Err(batch_err) = on_batch(&seqs) {
                err = Some(batch_err);
            }

            batch.clear();
            bytes = 0;
        }
    })?;
    if let Some(first_err) = err {
        return Err(first_err);
    }

    if !batch.is_empty() {
        let seqs: Vec<&[u8]> = batch.iter().map(|record| record.as_slice()).collect();

        on_batch(&seqs)?;
    }

    Ok(())
}

/// [`for_each_batch`] over a mapped plain-text file.
fn mapped_batches<F>(
    mapping: &memmap2::Mmap,
    batch_bytes: usize,
    on_batch: &mut F,
) -> Result<(), String>
where
    F: FnMut(&[&[u8]]) -> Result<(), String>,
{
    // Single-line records are borrowed from the mmap, and only wrapped records are joined into an owned copy.
    let data: &[u8] = mapping;
    let mut records = Records::new(data);
    let mut released = 0usize;

    loop {
        let mut batch: Vec<Cow<'_, [u8]>> = Vec::new();
        let mut bytes = 0usize;

        while bytes < batch_bytes && batch.len() < BATCH_RECORDS {
            match records.next() {
                Some(record) => {
                    bytes += record.len();
                    batch.push(record);
                }
                None => break,
            }
        }

        if batch.is_empty() {
            return Ok(());
        }

        let seqs: Vec<&[u8]> = batch.iter().map(|record| record.as_ref()).collect();

        on_batch(&seqs)?;
        drop(seqs);
        drop(batch);
        // Every byte before the next unparsed record is consumed, so those pages are dropped from this process while staying in the page cache.
        let upto = records.pos / PAGE * PAGE;

        if upto > released {
            // SAFETY: read-only file mapping with no live slice into [released, upto), since the batch was dropped, and a later access would just fault the page back in from the file.
            unsafe {
                let _ = mapping.unchecked_advise_range(
                    memmap2::UncheckedAdvice::DontNeed,
                    released,
                    upto - released,
                );
            }
            released = upto;
        }
    }
}

/// Returns the line of `data` starting at `*pos`, without its newline and trailing CR, spaces or tabs, and advances `*pos` past the newline.
fn next_line<'a>(data: &'a [u8], pos: &mut usize) -> &'a [u8] {
    let start = *pos;
    let end = match memchr(b'\n', &data[start..]) {
        Some(offset) => start + offset,
        None => data.len(),
    };

    *pos = (end + 1).min(data.len());
    let mut trimmed_end = end;

    // Trailing CR, spaces and tabs are not sequence, as in the C++ reader.
    while trimmed_end > start && matches!(data[trimmed_end - 1], b'\r' | b' ' | b'\t') {
        trimmed_end -= 1;
    }

    &data[start..trimmed_end]
}

/// Position of the first `needle` in `haystack`.
fn memchr(needle: u8, haystack: &[u8]) -> Option<usize> {
    haystack.iter().position(|&byte| byte == needle)
}

/// Incremental record parser over a mapped FASTA or FASTQ file, yielding each record's sequence.
struct Records<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Records<'a> {
    fn new(data: &'a [u8]) -> Self {
        let mut pos = 0usize;

        // Anything before the first header is skipped.
        while pos < data.len() && data[pos] != b'>' && data[pos] != b'@' {
            next_line(data, &mut pos);
        }

        Records { data, pos }
    }
}

impl<'a> Iterator for Records<'a> {
    type Item = Cow<'a, [u8]>;

    fn next(&mut self) -> Option<Cow<'a, [u8]>> {
        let data = self.data;
        let pos = &mut self.pos;

        if *pos >= data.len() {
            return None;
        }

        let hdr = data[*pos];

        next_line(data, pos); // header

        if hdr == b'>' {
            let first = next_line_seq(data, pos, b'>');
            let mut owned: Option<Vec<u8>> = None;

            while *pos < data.len() && data[*pos] != b'>' {
                let line = next_line(data, pos);

                if line.is_empty() {
                    continue;
                }

                owned
                    .get_or_insert_with(|| first.to_vec())
                    .extend_from_slice(line);
            }

            Some(match owned {
                Some(seq) => Cow::Owned(seq),
                None => Cow::Borrowed(first),
            })
        } else {
            // FASTQ sequence lines run until '+', and as many quality characters follow it.
            let first = next_line(data, pos);
            let mut owned: Option<Vec<u8>> = None;

            while *pos < data.len() && data[*pos] != b'+' {
                let line = next_line(data, pos);

                owned
                    .get_or_insert_with(|| first.to_vec())
                    .extend_from_slice(line);
            }

            let slen = owned.as_ref().map_or(first.len(), |seq| seq.len());

            next_line(data, pos); // '+' line
            let mut qual_len = 0usize;

            while qual_len < slen && *pos < data.len() {
                qual_len += next_line(data, pos).len();
            }

            // Blank lines between records are tolerated.
            while *pos < data.len() && data[*pos] != b'@' && data[*pos] != b'>' {
                next_line(data, pos);
            }

            Some(match owned {
                Some(seq) => Cow::Owned(seq),
                None => Cow::Borrowed(first),
            })
        }
    }
}

/// First non-empty sequence line of a FASTA record, or an empty slice if the record has none before the next `hdr` line.
fn next_line_seq<'a>(data: &'a [u8], pos: &mut usize, hdr: u8) -> &'a [u8] {
    while *pos < data.len() && data[*pos] != hdr {
        let line = next_line(data, pos);

        if !line.is_empty() {
            return line;
        }
    }

    &data[0..0]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch_size_grows_with_threads_up_to_the_cap() {
        assert_eq!(batch_bytes_for(0), BATCH_BYTES_PER_THREAD);
        assert_eq!(batch_bytes_for(1), BATCH_BYTES_PER_THREAD);
        assert_eq!(batch_bytes_for(4), 4 * BATCH_BYTES_PER_THREAD);
        assert_eq!(batch_bytes_for(64), BATCH_BYTES);
    }

    #[test]
    fn next_line_trims_line_endings_and_advances() {
        let data = b"ACGT \r\nTT\t\nlast";
        let mut pos = 0;

        assert_eq!(next_line(data, &mut pos), b"ACGT");
        assert_eq!(next_line(data, &mut pos), b"TT");
        assert_eq!(next_line(data, &mut pos), b"last");
        assert_eq!(pos, data.len());
    }

    #[test]
    fn next_line_seq_skips_blank_lines_and_stops_at_the_next_header() {
        let data = b"\n\nACGT\n>next\nGG\n";
        let mut pos = 0;

        assert_eq!(next_line_seq(data, &mut pos, b'>'), b"ACGT");
        assert_eq!(next_line_seq(data, &mut pos, b'>'), b"");
        assert_eq!(data[pos], b'>');
    }
}
