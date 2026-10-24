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

/// Batch bytes per worker of `query`, whose batches are `batch_bytes_for()` of its thread count, so one thread keeps 1 MiB of input resident and eight or more keep BATCH_BYTES.
pub(crate) const BATCH_BYTES_PER_THREAD: usize = 1 << 20;

/// The batch size for `threads` workers, BATCH_BYTES_PER_THREAD each, at most BATCH_BYTES.
pub(crate) fn batch_bytes_for(threads: usize) -> usize {
    (BATCH_BYTES_PER_THREAD * threads.max(1)).min(BATCH_BYTES)
}

/// Records per batch, which bounds the per-record bookkeeping when records are tiny (a million single k-mers would otherwise put about 250k records, and their slice vectors, in one batch).
pub(crate) const BATCH_RECORDS: usize = 1 << 15;

/// Page granularity for releasing consumed parts of the mapping.
const PAGE: usize = 4096;

/// The two bytes that start every gzip member (RFC 1952).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Calls `f` with consecutive batches of record sequences, in file order.
///
/// Every record of the file is passed exactly once, empty records included.
/// A batch closes at `batch_bytes` of sequence or [`BATCH_RECORDS`] records.
///
/// # Errors
///
/// Returns an error if the file cannot be read, and stops at the first error `f` returns.
pub(crate) fn for_each_batch<F>(path: &str, batch_bytes: usize, mut f: F) -> Result<(), String>
where
    F: FnMut(&[&[u8]]) -> Result<(), String>,
{
    if path != "-" {
        let file = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
        let len = file.metadata().map_err(|e| e.to_string())?.len();
        if len == 0 {
            return Ok(());
        }
        // SAFETY: the file is opened read-only and not modified while mapped.
        // The standard mmap caveat applies, since external truncation would fault.
        let m = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| e.to_string())?;
        if !m.starts_with(&GZIP_MAGIC) {
            return mapped_batches(&m, batch_bytes, &mut f);
        }
    }
    // For gzip or stdin the C++ reader decodes, and records are copied and batched.
    // Empty records are kept, as on the mapped path, because `query --ids` writes one line per record and dropping one would shift every later line onto the wrong read.
    let mut batch: Vec<Vec<u8>> = Vec::new();
    let mut bytes = 0usize;
    let mut err: Option<String> = None;
    iterate_fasta_all(path, |_name, seq| {
        if err.is_some() {
            return;
        }
        bytes += seq.len();
        batch.push(seq.to_vec());
        if bytes >= batch_bytes || batch.len() >= BATCH_RECORDS {
            let s: Vec<&[u8]> = batch.iter().map(|v| v.as_slice()).collect();
            if let Err(e) = f(&s) {
                err = Some(e);
            }
            batch.clear();
            bytes = 0;
        }
    })?;
    if let Some(e) = err {
        return Err(e);
    }
    if !batch.is_empty() {
        let s: Vec<&[u8]> = batch.iter().map(|v| v.as_slice()).collect();
        f(&s)?;
    }
    Ok(())
}

/// [`for_each_batch`] over a mapped plain-text file.
fn mapped_batches<F>(m: &memmap2::Mmap, batch_bytes: usize, f: &mut F) -> Result<(), String>
where
    F: FnMut(&[&[u8]]) -> Result<(), String>,
{
    let d: &[u8] = m;
    let mut it = Records::new(d);
    let mut released = 0usize;
    loop {
        let mut batch: Vec<Cow<'_, [u8]>> = Vec::new();
        let mut bytes = 0usize;
        while bytes < batch_bytes && batch.len() < BATCH_RECORDS {
            match it.next() {
                Some(r) => {
                    bytes += r.len();
                    batch.push(r);
                }
                None => break,
            }
        }
        if batch.is_empty() {
            return Ok(());
        }
        let s: Vec<&[u8]> = batch.iter().map(|c| c.as_ref()).collect();
        f(&s)?;
        drop(s);
        drop(batch);
        // Every byte before the next unparsed record is consumed, so those pages are dropped from this process while staying in the page cache.
        let upto = it.pos / PAGE * PAGE;
        if upto > released {
            // SAFETY: read-only file mapping with no live slice into [released, upto), since the batch was dropped, and a later access would just fault the page back in from the file.
            unsafe {
                let _ = m.unchecked_advise_range(
                    memmap2::UncheckedAdvice::DontNeed,
                    released,
                    upto - released,
                );
            }
            released = upto;
        }
    }
}

/// Returns the line of `d` starting at `*pos`, without its newline and trailing CR, spaces or tabs, and advances `*pos` past the newline.
fn next_line<'a>(d: &'a [u8], pos: &mut usize) -> &'a [u8] {
    let start = *pos;
    let end = match memchr(b'\n', &d[start..]) {
        Some(i) => start + i,
        None => d.len(),
    };
    *pos = (end + 1).min(d.len());
    let mut e = end;
    // trailing CR / spaces / tabs are not sequence (matches the C++ reader)
    while e > start && matches!(d[e - 1], b'\r' | b' ' | b'\t') {
        e -= 1;
    }
    &d[start..e]
}

/// Position of the first `c` in `s`.
fn memchr(c: u8, s: &[u8]) -> Option<usize> {
    s.iter().position(|&b| b == c)
}

/// Incremental record parser over a mapped FASTA or FASTQ file, yielding each record's sequence.
struct Records<'a> {
    d: &'a [u8],
    pos: usize,
}

impl<'a> Records<'a> {
    fn new(d: &'a [u8]) -> Self {
        let mut pos = 0usize;
        // skip anything before the first header
        while pos < d.len() && d[pos] != b'>' && d[pos] != b'@' {
            next_line(d, &mut pos);
        }
        Records { d, pos }
    }
}

impl<'a> Iterator for Records<'a> {
    type Item = Cow<'a, [u8]>;

    fn next(&mut self) -> Option<Cow<'a, [u8]>> {
        let d = self.d;
        let pos = &mut self.pos;
        if *pos >= d.len() {
            return None;
        }
        let hdr = d[*pos];
        next_line(d, pos); // header
        if hdr == b'>' {
            let first = next_line_seq(d, pos, b'>');
            let mut owned: Option<Vec<u8>> = None;
            while *pos < d.len() && d[*pos] != b'>' {
                let l = next_line(d, pos);
                if l.is_empty() {
                    continue;
                }
                owned
                    .get_or_insert_with(|| first.to_vec())
                    .extend_from_slice(l);
            }
            Some(match owned {
                Some(v) => Cow::Owned(v),
                None => Cow::Borrowed(first),
            })
        } else {
            // FASTQ has sequence lines until '+', then as many quality chars
            let first = next_line(d, pos);
            let mut owned: Option<Vec<u8>> = None;
            while *pos < d.len() && d[*pos] != b'+' {
                let l = next_line(d, pos);
                owned
                    .get_or_insert_with(|| first.to_vec())
                    .extend_from_slice(l);
            }
            let slen = owned.as_ref().map_or(first.len(), |v| v.len());
            next_line(d, pos); // '+' line
            let mut q = 0usize;
            while q < slen && *pos < d.len() {
                q += next_line(d, pos).len();
            }
            // tolerate blank lines between records
            while *pos < d.len() && d[*pos] != b'@' && d[*pos] != b'>' {
                next_line(d, pos);
            }
            Some(match owned {
                Some(v) => Cow::Owned(v),
                None => Cow::Borrowed(first),
            })
        }
    }
}

/// First non-empty sequence line of a FASTA record, or an empty slice if the record has none before the next `hdr` line.
fn next_line_seq<'a>(d: &'a [u8], pos: &mut usize, hdr: u8) -> &'a [u8] {
    while *pos < d.len() && d[*pos] != hdr {
        let l = next_line(d, pos);
        if !l.is_empty() {
            return l;
        }
    }
    &d[0..0]
}
