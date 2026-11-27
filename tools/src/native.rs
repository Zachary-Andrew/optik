//! Wraps the C++ core in safe Rust types.
//! Each index (FM, SPSS, dynamic SPSS, de Bruijn and seed) owns its handle and frees it on drop.
//! The readers and aligners are wrapped in the same way.

use super::*;

/// FM k-mer index over the separator-delimited unitigs, after the FMSI concept of Sladky et al. and implemented in include/fm_index.hpp.
/// A separator after every unitig means no length-k window spans a boundary, so membership is an FM-index backward search that finds the k-mer iff the SA interval is non-empty.
/// Lookups are exact.
pub(crate) struct FmIdx {
    pub(super) ptr: *mut OkFmIndex,
}

// SAFETY: FmIdx owns its index outright, and the C++ object has no thread affinity, so moving it to another thread is sound.
unsafe impl Send for FmIdx {}

// SAFETY: every &self method passes a const handle, and FM lookups never write to the index, so concurrent readers cannot race.
unsafe impl Sync for FmIdx {}

impl FmIdx {
    /// Builds the index from the k-mers of a FASTA/FASTQ file, keeping those seen at least `min_count` times, so 1 keeps every k-mer.
    pub(crate) fn build_from_fasta(
        fasta_path: &str,
        kmer_len: usize,
        min_count: u64,
        nthreads: u32,
    ) -> Result<Self, String> {
        let cpath = CString::new(fasta_path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okfm_build_from_fasta(
                cpath.as_ptr(),
                kmer_len as u64,
                min_count,
                nthreads,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "unknown build error"));
        }

        Ok(FmIdx { ptr })
    }

    /// Deserializes an index from a file.
    pub(crate) fn load(path: &str) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okfm_load(
                cpath.as_ptr(),
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "unknown load error"));
        }

        Ok(FmIdx { ptr })
    }

    /// Serializes the index to a file.
    pub(crate) fn save(&self, path: &str) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: the handle is live, cpath is NUL-terminated, and errbuf.len() is the buffer's real size.
        let status = unsafe {
            okfm_save(
                self.ptr,
                cpath.as_ptr(),
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if status != 0 {
            return Err(ffi_error(&errbuf, "unknown save error"));
        }

        Ok(())
    }

    /// Looks up every k-mer of `chunk` in one C++ call, returning (valid, positive, searches, extensions), where the FM index has no extensions.
    pub(crate) fn scan_chunk(&self, chunk: &[u8], kmer_len: usize) -> (u64, u64, u64, u64) {
        let (mut valid, mut positive, mut searches) = (0u64, 0u64, 0u64);

        // SAFETY: the handle is live, chunk supplies its own pointer and length, and the three counters are locals that outlive the call.
        unsafe {
            okfm_scan_seq(
                self.ptr,
                chunk.as_ptr() as *const c_char,
                chunk.len() as u64,
                kmer_len as u64,
                &mut valid,
                &mut positive,
                &mut searches,
            );
        }
        (valid, positive, searches, 0)
    }

    /// K-mer length used when building.
    pub(crate) fn kmer_len(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_k only reads it.
        unsafe { okfm_k(self.ptr) as usize }
    }
    /// Number of distinct canonical k-mers indexed.
    pub(crate) fn size(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_size only reads it.
        unsafe { okfm_size(self.ptr) as usize }
    }
    /// Number of maximal unitigs indexed.
    pub(crate) fn num_unitigs(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_num_unitigs only reads it.
        unsafe { okfm_num_unitigs(self.ptr) as usize }
    }
    /// Total bases of the unitigs, separators excluded.
    pub(crate) fn super_len(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_super_len only reads it.
        unsafe { okfm_super_len(self.ptr) as usize }
    }
    /// Length of the BWT (superstring + separators + sentinel).
    pub(crate) fn bwt_len(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_bwt_len only reads it.
        unsafe { okfm_bwt_len(self.ptr) as usize }
    }
    /// Size of the saved index file, in bytes.
    pub(crate) fn bytes(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okfm_bytes only reads it.
        unsafe { okfm_bytes(self.ptr) as usize }
    }
    /// Average bits per stored k-mer.
    pub(crate) fn bits_per_entry(&self) -> f64 {
        // SAFETY: the handle is live for as long as self is, and okfm_bits_per_entry only reads it.
        unsafe { okfm_bits_per_entry(self.ptr) }
    }
}

impl Drop for FmIdx {
    fn drop(&mut self) {
        // SAFETY: ptr came from a successful build or load, and drop runs once, so the index is freed exactly once.
        unsafe { okfm_free(self.ptr) }
    }
}

/// Compacted de Bruijn graph built by the k-mer counter (src/ffi_dbg.cpp).
pub(crate) struct DbgIdx {
    pub(super) ptr: *mut OkDbg,
}

// SAFETY: DbgIdx owns its graph outright, and the C++ object has no thread affinity.
unsafe impl Send for DbgIdx {}

// SAFETY: only save and core_size take &self, and both pass a const handle.
// write_gfa marks visited nodes in the graph, which is why it takes &mut self.
unsafe impl Sync for DbgIdx {}

impl DbgIdx {
    /// Counts the canonical k-mers of `fasta_path` seen at least `min_count` times into a graph on `nthreads` threads, as src/kmer_count.cpp does.
    /// `keep_nodes` keeps the node table that the `-s` archive needs.
    /// Otherwise the unitigs may come from partitioned compaction.
    pub(crate) fn build_from_fasta(
        fasta_path: &str,
        kmer_len: usize,
        min_count: u64,
        nthreads: usize,
        keep_nodes: bool,
    ) -> Result<Self, String> {
        let cpath = CString::new(fasta_path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okdbg_build_from_fasta(
                cpath.as_ptr(),
                kmer_len as u64,
                min_count,
                nthreads.max(1) as u32,
                keep_nodes as c_int,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "unknown build error"));
        }

        Ok(DbgIdx { ptr })
    }

    /// Serializes the graph core to a kdbg archive.
    pub(crate) fn save(&self, path: &str) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;

        // SAFETY: the handle is live and cpath is a NUL-terminated string that outlives the call.
        if unsafe { okdbg_save(self.ptr, cpath.as_ptr()) } != 0 {
            return Err("dbg save failed".into());
        }

        Ok(())
    }

    /// Distinct canonical k-mers in the graph.
    pub(crate) fn size(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okdbg_core_size only reads it.
        unsafe { okdbg_core_size(self.ptr) as usize }
    }

    /// Writes the maximal unitigs as a GFA1 graph of S and L lines, or as unitig FASTA when `fasta` is set, and returns the summary statistics.
    pub(crate) fn write_gfa(
        &mut self,
        path: &str,
        nthreads: usize,
        fasta: bool,
    ) -> Result<GfaStats, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut stats = GfaStats::default();
        // SAFETY: &mut self gives this call the graph exclusively, cpath outlives it, and stats is a local the C++ side fills.
        let status = unsafe {
            okdbg_write_gfa(
                self.ptr,
                cpath.as_ptr(),
                nthreads.max(1) as u32,
                &mut stats,
                fasta as c_int,
            )
        };

        if status != 0 {
            return Err(format!(
                "failed to write {} to {path}",
                if fasta { "FASTA" } else { "GFA" }
            ));
        }

        Ok(stats)
    }
}

impl Drop for DbgIdx {
    fn drop(&mut self) {
        // SAFETY: ptr came from a successful build, and drop runs once, so the graph is freed exactly once.
        unsafe { okdbg_free(self.ptr) }
    }
}

/// Static SPSS core, exact k-mer membership over super-k-mer strings with a minimizer MPHF (spss_core.hpp).
pub(crate) struct SpssIdx {
    pub(super) ptr: *mut OkSpss,
}

// SAFETY: SpssIdx owns its index outright, and the C++ object has no thread affinity.
unsafe impl Send for SpssIdx {}

// SAFETY: every &self method passes a const handle, and SPSS scans keep their lane state on the stack, so concurrent readers share nothing mutable.
unsafe impl Sync for SpssIdx {}

impl SpssIdx {
    /// Builds from FASTA/FASTQ.
    /// `compact` counts the k-mers seen at least `min_count` times and indexes their maximal unitigs.
    /// Otherwise the input sequences are indexed as is.
    pub(crate) fn build_from_fasta(
        path: &str,
        kmer_len: usize,
        minimizer_len: usize,
        compact: bool,
        min_count: u64,
        nthreads: usize,
    ) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okspss_build_from_fasta(
                cpath.as_ptr(),
                kmer_len as u64,
                minimizer_len as u64,
                compact as c_int,
                min_count,
                nthreads.max(1) as u32,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "SPSS build failed"));
        }

        Ok(SpssIdx { ptr })
    }
    /// Maps an index, where `verify` first checks its checksum by reading the file once.
    pub(crate) fn load(path: &str, verify: bool) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okspss_load(
                cpath.as_ptr(),
                verify as c_int,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "SPSS load failed"));
        }

        Ok(SpssIdx { ptr })
    }
    /// Writes the index to `path` (atomically, src/ffi_spss.cpp).
    pub(crate) fn save(&self, path: &str) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;

        // SAFETY: the handle is live and cpath is a NUL-terminated string that outlives the call.
        if unsafe { okspss_save(self.ptr, cpath.as_ptr()) } != 0 {
            return Err(format!("failed to write {path}"));
        }

        Ok(())
    }
    /// Indexed k-mer positions (a repeated k-mer counts once per position).
    pub(crate) fn kmers(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okspss_kmers only reads it.
        unsafe { okspss_kmers(self.ptr) }
    }
    /// In-memory size of the index, in bytes.
    pub(crate) fn bytes(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okspss_bytes only reads it.
        unsafe { okspss_bytes(self.ptr) }
    }
    pub(crate) fn kmer_len(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okspss_k only reads it.
        unsafe { okspss_k(self.ptr) as usize }
    }
    pub(crate) fn minimizer_len(&self) -> usize {
        // SAFETY: the handle is live for as long as self is, and okspss_m only reads it.
        unsafe { okspss_m(self.ptr) as usize }
    }
    pub(crate) fn superkmers(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okspss_superkmers only reads it.
        unsafe { okspss_superkmers(self.ptr) }
    }
    /// Distinct k-mers (kmers() counts positions, which repeat duplicates).
    pub(crate) fn distinct(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okspss_distinct only reads it.
        unsafe { okspss_distinct(self.ptr) }
    }
    /// Dense k-mer identifier of every window of `seq`, -1 when absent or N.
    pub(crate) fn ids(&self, seq: &[u8]) -> Vec<i64> {
        if seq.len() < self.kmer_len() {
            return Vec::new();
        }

        let mut out = vec![0i64; seq.len() - self.kmer_len() + 1];

        // SAFETY: the handle is live, seq supplies its own pointer and length, and out has one slot for each of the seq.len() - k + 1 windows the C++ side writes.
        unsafe {
            okspss_ids(
                self.ptr,
                seq.as_ptr() as *const c_char,
                seq.len() as u64,
                out.as_mut_ptr(),
            );
        }
        out
    }
    /// Streaming lookup of every k-mer of `chunk`, returning (valid, positive, full lookups, extensions).
    pub(crate) fn scan_chunk(&self, chunk: &[u8]) -> (u64, u64, u64, u64) {
        let (mut valid, mut positive, mut searches, mut extensions) = (0u64, 0u64, 0u64, 0u64);

        // SAFETY: the handle is live, chunk supplies its own pointer and length, and the four counters are locals that outlive the call.
        unsafe {
            okspss_scan_seq(
                self.ptr,
                chunk.as_ptr() as *const c_char,
                chunk.len() as u64,
                &mut valid,
                &mut positive,
                &mut searches,
                &mut extensions,
            );
        }
        (valid, positive, searches, extensions)
    }

    /// Streaming lookup over a batch of records with cache misses overlapping across records, returning (valid, positive, searches, extensions).
    pub(crate) fn scan_many(&self, recs: &[&[u8]]) -> (u64, u64, u64, u64) {
        let ptrs: Vec<*const c_char> = recs
            .iter()
            .map(|record| record.as_ptr() as *const c_char)
            .collect();
        let lens: Vec<u64> = recs.iter().map(|record| record.len() as u64).collect();
        let (mut valid, mut positive, mut searches, mut extensions) = (0u64, 0u64, 0u64, 0u64);

        // SAFETY: ptrs and lens both hold recs.len() entries taken from live slices, and the counters are locals that outlive the call.
        unsafe {
            okspss_scan_many(
                self.ptr,
                ptrs.as_ptr(),
                lens.as_ptr(),
                recs.len() as u64,
                &mut valid,
                &mut positive,
                &mut searches,
                &mut extensions,
            );
        }
        (valid, positive, searches, extensions)
    }
}

impl Drop for SpssIdx {
    fn drop(&mut self) {
        // SAFETY: ptr came from a successful build or load, and drop runs once, so the index is freed exactly once.
        unsafe { okspss_free(self.ptr) }
    }
}

/// SPSS index opened for editing, i.e. the static core plus concurrent open-addressing edit sets (include/dyn_spss.hpp).
pub(crate) struct DynSpssIdx {
    pub(super) ptr: *mut OkDynSpss,
}

impl DynSpssIdx {
    /// Loads the index at `path`, checksum verified.
    pub(crate) fn open(path: &str) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: cpath and errbuf outlive the call, and errbuf.len() is the buffer's real size.
        let ptr = unsafe {
            okdyn_open(
                cpath.as_ptr(),
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if ptr.is_null() {
            return Err(ffi_error(&errbuf, "SPSS load failed"));
        }

        Ok(DynSpssIdx { ptr })
    }
    /// Inserts every k-mer of a FASTA/FASTQ file when `add` is set, or erases them, and returns how many k-mers changed membership.
    pub(crate) fn apply_fasta(&mut self, path: &str, add: bool) -> Result<u64, String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        // SAFETY: &mut self gives this call the edit sets exclusively, while cpath and errbuf outlive it.
        let changed = unsafe {
            okdyn_apply_fasta(
                self.ptr,
                cpath.as_ptr(),
                add as c_int,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };

        if changed < 0 {
            return Err(ffi_error(&errbuf, "failed to read k-mers"));
        }

        Ok(changed as u64)
    }
    /// Size of the insertion set.
    pub(crate) fn added(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okdyn_added only reads it.
        unsafe { okdyn_added(self.ptr) }
    }
    /// Size of the erasure set.
    pub(crate) fn deleted(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okdyn_deleted only reads it.
        unsafe { okdyn_deleted(self.ptr) }
    }
    /// k-mer positions of the static core.
    pub(crate) fn kmers(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okdyn_kmers only reads it.
        unsafe { okdyn_kmers(self.ptr) }
    }
    /// Compacts the edits into the static core and saves the result to `path`.
    pub(crate) fn compact_save(&mut self, path: &str, nthreads: usize) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|err| err.to_string())?;

        // SAFETY: &mut self gives this call the index exclusively, and cpath is a NUL-terminated string that outlives it.
        if unsafe { okdyn_compact_save(self.ptr, cpath.as_ptr(), nthreads.max(1) as u32) } != 0 {
            return Err(format!("failed to compact/write {path}"));
        }

        Ok(())
    }
}

impl Drop for DynSpssIdx {
    fn drop(&mut self) {
        // SAFETY: ptr came from a successful open, and drop runs once, so the index is freed exactly once.
        unsafe { okdyn_free(self.ptr) }
    }
}

// The first four bytes of each index file, read as a little-endian u32 (include/fm_index.hpp and src/ffi_spss.cpp write them).

pub(crate) const FM_MAGIC: u32 = u32::from_le_bytes(*b"OKFM"); // k-order BWT FM index

pub(crate) const SPSS_MAGIC: u32 = u32::from_le_bytes(*b"OKSP"); // SPSS index, "OKSPSS" padded to 8 bytes

/// Every index file starts with a 4-byte magic.
pub(crate) const INDEX_MAGIC_LEN: usize = 4;

/// Returns the first 4 bytes of an index file, its format magic, or None if the file is shorter than that.
/// An unreadable file is an error naming the path.
pub(crate) fn index_magic(path: &str) -> Result<Option<u32>, String> {
    use std::fs::File;
    use std::io::Read;
    let mut file = File::open(path).map_err(|err| format!("cannot open index {path}: {err}"))?;
    let mut magic = [0u8; INDEX_MAGIC_LEN];

    Ok(file
        .read_exact(&mut magic)
        .ok()
        .map(|_| u32::from_le_bytes(magic)))
}

/// Calls `f(name, seq)` for every record of a FASTA/FASTQ file with a non-empty sequence, read by include/fasta_reader.hpp with line breaks inside a sequence removed.
pub(crate) fn iterate_fasta<F: FnMut(&str, &[u8])>(path: &str, on_record: F) -> Result<(), String> {
    iterate_fasta_impl(path, on_record, false)
}

/// As [`iterate_fasta`], but records with an empty sequence are reported too.
pub(crate) fn iterate_fasta_all<F: FnMut(&str, &[u8])>(
    path: &str,
    on_record: F,
) -> Result<(), String> {
    iterate_fasta_impl(path, on_record, true)
}

/// What a record callback returns to keep the C++ reader going.
const RECORD_CONTINUE: c_int = 0;

/// What a record callback returns to stop the C++ reader, which then unwinds and reports failure.
const RECORD_STOP: c_int = 1;

/// The state a record callback shares with its reader call.
struct RecordVisit<'p, F> {
    on_record: F,
    path: &'p str,
    /// Records seen so far, which numbers a rejected record when `exact_count` holds.
    seen: u64,
    /// Whether `seen` counts every record of the file, which is true only when empty records are delivered too.
    exact_count: bool,
    /// Why the read was stopped, if a callback stopped it.
    error: Option<String>,
}

impl<F> RecordVisit<'_, F> {
    fn new(on_record: F, path: &str, exact_count: bool) -> RecordVisit<'_, F> {
        RecordVisit { on_record, path, seen: 0, exact_count, error: None }
    }

    /// The record name as text, or None after recording why it is rejected.
    /// Names must be valid UTF-8, and a name that is not is refused rather than altered, because SAM output allows only printable ASCII anyway.
    fn name<'n>(&mut self, name: &'n [u8]) -> Option<&'n str> {
        self.seen += 1;

        match std::str::from_utf8(name) {
            Ok(text) => Some(text),
            Err(_) => {
                let which = if self.exact_count {
                    format!("record {}", self.seen)
                } else {
                    "a record".to_string()
                };

                self.error = Some(format!(
                    "{}: {which} has a name that is not valid UTF-8 (\"{}\"). \
                     Rename it using ASCII characters, which SAM requires.",
                    self.path,
                    name.escape_ascii()
                ));
                None
            }
        }
    }

    /// The reader's outcome, where an error recorded by a callback takes precedence over the reader's own message.
    fn finish(self, status: c_int, errbuf: &[c_char]) -> Result<(), String> {
        if let Some(err) = self.error {
            return Err(err);
        }

        if status != 0 {
            // SAFETY: errbuf starts zeroed and the C++ side copies at most len - 1 bytes into it, so it stays NUL-terminated.
            let msg = unsafe { std::ffi::CStr::from_ptr(errbuf.as_ptr()) };

            return Err(msg.to_string_lossy().into_owned());
        }

        Ok(())
    }
}

/// A byte slice over `len` bytes at `ptr`, or an empty slice when `len` is 0, whatever `ptr` is.
///
/// # Safety
/// When `len` is nonzero, `ptr` must point at `len` readable bytes that stay valid and unchanged for `'a`.
unsafe fn raw_bytes<'a>(ptr: *const c_char, len: u64) -> &'a [u8] {
    if len == 0 {
        return &[];
    }

    // SAFETY: the caller guarantees that ptr points at len readable bytes valid for 'a, and len is nonzero here, so ptr is not null.
    unsafe { std::slice::from_raw_parts(ptr as *const u8, len as usize) }
}

/// Calls `f(name, seq, qual)` for every record, empty ones included.
/// `qual` is the FASTQ quality string, empty for FASTA.
/// A record whose name is not valid UTF-8 stops the read with an error that names it.
pub(crate) fn iterate_reads<F: FnMut(&str, &[u8], &[u8])>(path: &str, on_record: F) -> Result<(), String> {
    extern "C" fn read_record_callback<F: FnMut(&str, &[u8], &[u8])>(
        cookie: *mut c_void,
        name: *const c_char,
        nlen: u64,
        seq: *const c_char,
        slen: u64,
        qual: *const c_char,
        qlen: u64,
    ) -> c_int {
        // SAFETY: cookie is the `&mut visit` that iterate_reads passes to the reader, which outlives the call, and the reader calls back on this thread only.
        let visit = unsafe { &mut *(cookie as *mut RecordVisit<'_, F>) };
        // SAFETY: the reader passes each field as a pointer and length that stay valid until the callback returns.
        let (name_b, seq_b, qual_b) = unsafe { (raw_bytes(name, nlen), raw_bytes(seq, slen), raw_bytes(qual, qlen)) };
        let Some(name_s) = visit.name(name_b) else {
            return RECORD_STOP;
        };

        (visit.on_record)(name_s, seq_b, qual_b);
        RECORD_CONTINUE
    }
    let cpath = CString::new(path).map_err(|err| err.to_string())?;
    let mut errbuf = vec![0 as c_char; READER_ERRBUF_LEN];
    let mut visit = RecordVisit::new(on_record, path, true);
    // SAFETY: cpath and errbuf outlive the call, errbuf.len() is its real capacity, and the cookie points at `visit` with the type the callback casts it back to.
    let status = unsafe {
        ok_iterate_reads(
            cpath.as_ptr(),
            read_record_callback::<F>,
            &mut visit as *mut RecordVisit<'_, F> as *mut c_void,
            errbuf.as_mut_ptr(),
            errbuf.len() as u64,
        )
    };

    visit.finish(status, &errbuf)
}

/// Shared body of [`iterate_fasta`] and [`iterate_fasta_all`].
/// A record whose name is not valid UTF-8 stops the read with an error that names it.
pub(super) fn iterate_fasta_impl<F: FnMut(&str, &[u8])>(path: &str, on_record: F, keep_empty: bool) -> Result<(), String> {
    extern "C" fn fasta_record_callback<F: FnMut(&str, &[u8])>(
        cookie: *mut c_void,
        name: *const c_char,
        nlen: u64,
        seq: *const c_char,
        slen: u64,
    ) -> c_int {
        // SAFETY: cookie is the `&mut visit` that iterate_fasta_impl passes to the reader, which outlives the call, and the reader calls back on this thread only.
        let visit = unsafe { &mut *(cookie as *mut RecordVisit<'_, F>) };
        // SAFETY: the reader passes each field as a pointer and length that stay valid until the callback returns.
        let (name_b, seq_b) = unsafe { (raw_bytes(name, nlen), raw_bytes(seq, slen)) };
        let Some(name_s) = visit.name(name_b) else {
            return RECORD_STOP;
        };

        (visit.on_record)(name_s, seq_b);
        RECORD_CONTINUE
    }
    let cpath = CString::new(path).map_err(|err| err.to_string())?;
    let mut errbuf = vec![0 as c_char; READER_ERRBUF_LEN];
    // Without empty records the C++ side skips some, so only keep_empty gives an exact record number.
    let mut visit = RecordVisit::new(on_record, path, keep_empty);
    let iter = if keep_empty {
        ok_iterate_fasta_all
    } else {
        ok_iterate_fasta
    };
    // SAFETY: cpath and errbuf outlive the call, errbuf.len() is its real capacity, and the cookie points at `visit` with the type the callback casts it back to.
    let status = unsafe {
        iter(
            cpath.as_ptr(),
            fasta_record_callback::<F>,
            &mut visit as *mut RecordVisit<'_, F> as *mut c_void,
            errbuf.as_mut_ptr(),
            errbuf.len() as u64,
        )
    };

    visit.finish(status, &errbuf)
}

/// Mapper seed index, a minimizer MPHF with CSR postings and the reference as key store (seed_index.hpp).
/// It reads the reference sequences in place, so the `'r` borrow keeps them alive and unchanged for as long as the index exists.
pub(crate) struct SeedIdx<'r> {
    pub(super) ptr: *mut OkSeedIdx,
    pub(super) refs: std::marker::PhantomData<&'r [Vec<u8>]>,
}

// SAFETY: SeedIdx owns its index outright, and the 'r borrow keeps the references it reads alive wherever it moves.
unsafe impl Send for SeedIdx<'_> {}

// SAFETY: every &self method passes a const handle, and the hit buffers that collect returns are thread_local on the C++ side, so threads never share one.
unsafe impl Sync for SeedIdx<'_> {}

impl<'r> SeedIdx<'r> {
    /// Builds the seed index over `refs` with minimizer length `kmer_len`, window `window` and repeat cutoff `max_occ` (include/seed_index.hpp).
    pub(crate) fn build(
        refs: &'r [Vec<u8>],
        kmer_len: usize,
        window: usize,
        max_occ: u32,
        threads: usize,
    ) -> Result<Self, String> {
        // Seed postings pack a position in 31 bits, so an over-long reference is reported by name instead of failing opaquely or wrapping its positions.
        if let Some((seq_idx, long_ref)) = refs
            .iter()
            .enumerate()
            .find(|(_, ref_seq)| ref_seq.len() >= MAX_REF_SEQ_LEN)
        {
            return Err(format!(
                "reference sequence #{} is {} bases; map supports sequences \
                                shorter than 2^31 (2,147,483,648) bases -- split it",
                seq_idx + 1,
                long_ref.len()
            ));
        }

        let ptrs: Vec<*const c_char> = refs
            .iter()
            .map(|ref_seq| ref_seq.as_ptr() as *const c_char)
            .collect();
        let lens: Vec<u64> = refs.iter().map(|ref_seq| ref_seq.len() as u64).collect();
        // SAFETY: ptrs and lens hold refs.len() entries from live slices.
        // The index keeps reading those sequences after the call, which the 'r borrow on refs allows.
        let ptr = unsafe {
            okseed_build(
                ptrs.as_ptr(),
                lens.as_ptr(),
                refs.len() as u64,
                kmer_len as u32,
                window as u32,
                max_occ,
                threads as u32,
            )
        };

        if ptr.is_null() {
            return Err(
                "seed index build failed (more than 2^32 seed postings -- a reference of \
                        ~25 Gbp or more -- or out of memory)"
                    .into(),
            );
        }

        Ok(SeedIdx { ptr, refs: std::marker::PhantomData })
    }
    /// Distinct minimizers indexed.
    pub(crate) fn distinct(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okseed_distinct only reads it.
        unsafe { okseed_distinct(self.ptr) }
    }
    /// In-memory size of the index, in bytes.
    pub(crate) fn bytes(&self) -> u64 {
        // SAFETY: the handle is live for as long as self is, and okseed_bytes only reads it.
        unsafe { okseed_bytes(self.ptr) }
    }
    /// Fraction of the reference's minimizer occurrences that repeat an earlier one (1 - distinct / occurrences), 0 for a genome without repeats.
    /// A genome-level MAPQ input (map.rs, hidden_odds).
    pub(crate) fn repeat_fraction(&self) -> f64 {
        // SAFETY: the handle is live for as long as self is, and okseed_occurrences only reads it.
        let occ = unsafe { okseed_occurrences(self.ptr) };

        if occ == 0 {
            0.0
        } else {
            1.0 - self.distinct() as f64 / occ as f64
        }
    }
    /// Seed hits of a read on both strands, copied out of the C++ buffer.
    pub(crate) fn collect(&self, seq: &[u8]) -> Vec<SeedHit> {
        self.collect_with_rep(seq).0
    }
    /// Seed hits, the read's minimizer count, and how many of its minimizers were dropped as over-occurring (repetitive).
    pub(crate) fn collect_with_rep(&self, seq: &[u8]) -> (Vec<SeedHit>, u32, u32) {
        let mut hits_ptr: *const SeedHit = std::ptr::null();
        let (mut n_min, mut n_rep) = (0u32, 0u32);
        // SAFETY: the handle is live, seq supplies its own pointer and length, and the three out-parameters are locals.
        let hit_count = unsafe {
            okseed_collect(
                self.ptr,
                seq.as_ptr() as *const c_char,
                seq.len() as u64,
                &mut hits_ptr,
                &mut n_min,
                &mut n_rep,
            )
        };
        let hits = if hit_count == 0 || hits_ptr.is_null() {
            Vec::new()
        } else {
            // SAFETY: a non-null hits_ptr points at hit_count hits in this thread's buffer, which stays untouched until this thread collects again, and to_vec() copies them first.
            unsafe { std::slice::from_raw_parts(hits_ptr, hit_count as usize) }.to_vec()
        };

        (hits, n_min, n_rep)
    }
    /// The chaining path's occurrence threshold (seed_index.hpp mid_occ).
    pub(crate) fn mid_occ(&self) -> u32 {
        // SAFETY: the handle is live for as long as self is, and okseed_mid_occ only reads it.
        unsafe { okseed_mid_occ(self.ptr) }
    }
    /// Seed hits for chaining (repeat runs keep their least frequent minimizers) and how many read bases lie only under dropped minimizers.
    pub(crate) fn collect_chain(&self, seq: &[u8]) -> (Vec<SeedHit>, u32) {
        let mut hits_ptr: *const SeedHit = std::ptr::null();
        let mut rep_len = 0u32;
        // SAFETY: the handle is live, seq supplies its own pointer and length, and the two out-parameters are locals.
        let hit_count = unsafe {
            okseed_collect_chain(
                self.ptr,
                seq.as_ptr() as *const c_char,
                seq.len() as u64,
                &mut hits_ptr,
                &mut rep_len,
            )
        };
        let hits = if hit_count == 0 || hits_ptr.is_null() {
            Vec::new()
        } else {
            // SAFETY: as in collect_with_rep, the thread's own buffer holds hit_count hits until its next collect, and to_vec() copies them first.
            unsafe { std::slice::from_raw_parts(hits_ptr, hit_count as usize) }.to_vec()
        };

        (hits, rep_len)
    }
}

impl Drop for SeedIdx<'_> {
    fn drop(&mut self) {
        // SAFETY: ptr came from a successful build, and drop runs once, so the index is freed exactly once.
        unsafe { okseed_free(self.ptr) }
    }
}

/// Placements of one query in a text, i.e. the best edit distance capped at max_ed + 1 with its exclusive end offset, and the best distance among ends more than `sep` away, a distinct placement, with its end.
#[derive(Clone, Copy)]
pub(crate) struct Placements {
    pub best: i32,
    pub end: i32,
    pub alt: i32,
    pub alt_end: i32,
}

/// Texts [`placements_x4`] aligns at once, one per AVX2 lane.
pub(crate) const PLACEMENT_LANES: usize = 4;

/// Semi-global placements of one query in up to [`PLACEMENT_LANES`] texts, with the query end to end, the text ends free and unit costs, computed by the bit-parallel aligner of include/bpm_edit.hpp.
pub(crate) fn placements_x4(
    query: &[u8],
    texts: &[&[u8]],
    max_ed: i32,
    sep: i32,
) -> [Placements; PLACEMENT_LANES] {
    debug_assert!(texts.len() <= PLACEMENT_LANES);
    let ptrs: Vec<*const c_char> = texts
        .iter()
        .map(|text| text.as_ptr() as *const c_char)
        .collect();
    let lens: Vec<i32> = texts.iter().map(|text| text.len() as i32).collect();
    // Unused lanes keep max_ed + 1, the no-placement value callers already test for.
    let (mut best, mut best_end, mut alt, mut alt_end) = (
        [max_ed + 1; PLACEMENT_LANES],
        [0i32; PLACEMENT_LANES],
        [max_ed + 1; PLACEMENT_LANES],
        [0i32; PLACEMENT_LANES],
    );

    // SAFETY: ptrs and lens hold texts.len() entries, at most PLACEMENT_LANES of which are passed.
    // The four output arrays each have PLACEMENT_LANES slots, one per text lane.
    unsafe {
        ok_placements_x4(
            query.as_ptr() as *const c_char,
            query.len() as i32,
            ptrs.as_ptr(),
            lens.as_ptr(),
            texts.len().min(PLACEMENT_LANES) as i32,
            max_ed,
            sep,
            best.as_mut_ptr(),
            best_end.as_mut_ptr(),
            alt.as_mut_ptr(),
            alt_end.as_mut_ptr(),
        );
    }
    let mut out = [Placements {
        best: max_ed + 1,
        end: 0,
        alt: max_ed + 1,
        alt_end: 0,
    }; PLACEMENT_LANES];

    for i in 0..PLACEMENT_LANES {
        out[i] = Placements {
            best: best[i],
            end: best_end[i],
            alt: alt[i],
            alt_end: alt_end[i],
        };
    }

    out
}

/// Result of a gap-affine alignment, as CIGAR, score and text end.
pub(crate) struct AffineAln {
    pub cigar: String,
    pub score: i32,
    pub tend: usize,
}

/// A CIGAR operation takes at most this many bytes in the C++ writer's text (count digits plus the op).
pub(super) const CIGAR_BYTES_PER_OP: usize = 8;

/// Spare bytes beyond the per-op bound in the CIGAR buffers of the affine and BiWFA aligners and of the semi-global CIGAR.
pub(super) const AFFINE_CIGAR_SLACK: usize = 32;

pub(super) const SEMIGLOBAL_CIGAR_SLACK: usize = 16;

/// A buffer large enough for the CIGAR of aligning `qlen` against `tlen` bases, i.e. at most one operation per base, each at most 8 bytes, plus `slack`.
pub(super) fn cigar_buffer(qlen: usize, tlen: usize, slack: usize) -> Vec<u8> {
    vec![0u8; slack + CIGAR_BYTES_PER_OP * (qlen + tlen)]
}

/// Gap-affine alignment with the query end to end from text position 0, while the text may end anywhere in [tlen - tef, tlen].
/// A negative `max_score` leaves the score unbounded.
pub(crate) fn wfa_affine_align(
    query: &[u8],
    text: &[u8],
    tef: usize,
    pen: (i32, i32, i32),
    max_score: i32,
) -> Option<AffineAln> {
    let mut buf = cigar_buffer(query.len(), text.len(), AFFINE_CIGAR_SLACK);
    let (mut score, mut text_end) = (0i32, 0i32);
    // SAFETY: query and text supply their own lengths, buf.len() is the CIGAR buffer's real size, and the score and end are locals.
    let status = unsafe {
        ok_wfa_affine_align(
            query.as_ptr() as *const c_char,
            query.len() as i32,
            text.as_ptr() as *const c_char,
            text.len() as i32,
            tef as i32,
            pen.0,
            pen.1,
            pen.2,
            max_score,
            buf.as_mut_ptr() as *mut c_char,
            buf.len() as i32,
            &mut score,
            &mut text_end,
        )
    };

    if status != 0 {
        return None;
    }

    let cigar = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();

    Some(AffineAln {
        cigar,
        score,
        tend: text_end as usize,
    })
}

/// End-to-end gap-affine alignment in O(s) memory by BiWFA.
/// The result is the same as `wfa_affine_align(query, text, 0, ..)` with tend = text length.
pub(crate) fn biwfa_align(
    query: &[u8],
    text: &[u8],
    pen: (i32, i32, i32),
    max_score: i32,
) -> Option<AffineAln> {
    let mut buf = cigar_buffer(query.len(), text.len(), AFFINE_CIGAR_SLACK);
    let mut score = 0i32;
    // SAFETY: query and text supply their own lengths, buf.len() is the CIGAR buffer's real size, and score is a local.
    let status = unsafe {
        ok_biwfa_align(
            query.as_ptr() as *const c_char,
            query.len() as i32,
            text.as_ptr() as *const c_char,
            text.len() as i32,
            pen.0,
            pen.1,
            pen.2,
            max_score,
            buf.as_mut_ptr() as *mut c_char,
            buf.len() as i32,
            &mut score,
        )
    };

    if status != 0 {
        return None;
    }

    let cigar = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();

    Some(AffineAln {
        cigar,
        score,
        tend: text.len(),
    })
}

/// Penalties of [`wfa_affine_cigar`], as mismatch, gap open and gap extend.
pub(super) const CIGAR_PENALTIES: (i32, i32, i32) = (4, 6, 2);

/// Gap-affine semi-global CIGAR (include/wfa_affine.hpp) with [`CIGAR_PENALTIES`], returned as (cigar, 0-based text start), or None on failure.
pub(crate) fn wfa_affine_cigar(query: &[u8], text: &[u8]) -> Option<(String, usize)> {
    let mut buf = cigar_buffer(query.len(), text.len(), SEMIGLOBAL_CIGAR_SLACK);
    let mut score = 0i32;
    // SAFETY: query and text supply their own lengths, buf.len() is the CIGAR buffer's real size, and score is a local.
    let status = unsafe {
        ok_wfa_affine_cigar(
            query.as_ptr() as *const c_char,
            query.len() as i32,
            text.as_ptr() as *const c_char,
            text.len() as i32,
            CIGAR_PENALTIES.0,
            CIGAR_PENALTIES.1,
            CIGAR_PENALTIES.2,
            buf.as_mut_ptr() as *mut c_char,
            buf.len() as i32,
            &mut score,
        )
    };

    // A non-negative status is the text start rather than an error code.
    if status < 0 {
        return None;
    }

    let cig = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();

    Some((cig, status as usize))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A header byte that is valid Latin-1 (e with an acute accent) but never valid UTF-8 on its own.
    const LATIN1_E_ACUTE: u8 = 0xE9;

    /// Writes two records whose second name holds a Latin-1 byte, and returns the file's path.
    fn latin1_fasta(tag: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("optik_utf8_{tag}_{}.fa", std::process::id()));
        let mut body = b">ok\nACGT\n>caf".to_vec();

        body.push(LATIN1_E_ACUTE);
        body.extend_from_slice(b"\nACGT\n");
        std::fs::write(&path, body).expect("write the test FASTA");
        path
    }

    #[test]
    fn non_utf8_names_are_rejected_with_the_exact_bytes() {
        let path = latin1_fasta("reads");
        let path_s = path.to_str().expect("temp paths are UTF-8");
        let mut names = Vec::new();
        let err = iterate_reads(path_s, |name, _, _| names.push(name.to_string())).unwrap_err();

        assert_eq!(names, ["ok"], "the record before the bad one is still delivered");
        assert!(err.contains("record 2"), "{err}");
        assert!(err.contains("caf\\xe9"), "{err}");

        let mut fasta_names = Vec::new();
        let fasta_err = iterate_fasta(path_s, |name, _| fasta_names.push(name.to_string())).unwrap_err();

        assert_eq!(fasta_names, ["ok"]);
        assert!(fasta_err.contains("not valid UTF-8") && fasta_err.contains("caf\\xe9"), "{fasta_err}");
        std::fs::remove_file(&path).ok();
    }
}
