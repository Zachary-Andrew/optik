//! optik, a genomic sequence analysis toolkit.
//!
//! This file holds the command line, which checks arguments and dispatches, the C ABI declarations of the C++ core, and the safe Rust wrappers every subcommand uses.
//! The C++ core in include/ and src/ does the k-mer counting, minimizer selection, the bit-parallel and wavefront aligners, the de Bruijn graph and the SPSS and FM indexes.
//!
//! Subcommands
//!   tui        Interactive terminal UI
//!   map        Read mapping via minimizer seeds
//!   dbruijn    De Bruijn graph assembly
//!   build-spss Build the recommended SPSS membership index, also called `build`
//!   spss-update Add or remove k-mers of an SPSS index
//!   build-fm   Build the FM membership index over a k-order BWT
//!   query      Query an SPSS or FM index

mod chain;
mod map;
mod seqio;

mod build;
mod dbruijn;
mod query;
mod tui;

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::{c_char, c_int, c_void};

/// C++ OkFmIndex, the FM k-mer index of src/ffi_fm.cpp, seen as an opaque type.
pub(crate) struct OkFmIndex {
    _private: [u8; 0],
}
/// The compacted de Bruijn graph (OkDbg, src/ffi_dbg.cpp).
pub(crate) struct OkDbg {
    _private: [u8; 0],
}
/// The mapper's minimizer seed index (OkSeedIdx).
pub(crate) struct OkSeedIdx {
    _private: [u8; 0],
}
/// One seed hit, mirroring SeedIndex::Hit, where qpos is in the strand's read frame.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct SeedHit {
    pub rid: u32,
    pub rev: u32,
    pub rpos: u32,
    pub qpos: u32,
}
/// The static SPSS index (OkSpss, src/ffi_spss.cpp).
pub(crate) struct OkSpss {
    _private: [u8; 0],
}
/// An SPSS index with edit sets (OkDynSpss).
pub(crate) struct OkDynSpss {
    _private: [u8; 0],
}

// The C ABI of the C++ core in src/ffi*.cpp.
// improper_ctypes is allowed because the opaque handle structs above are not #[repr(C)], and they are only ever used behind pointers.

/// Size of the error buffers passed to the index functions.
const ERRBUF_LEN: usize = 4096;
/// Size of the error buffers passed to the FASTA/FASTQ readers.
const READER_ERRBUF_LEN: usize = 512;

#[allow(improper_ctypes)]
unsafe extern "C" {
    // Sequence readers and global settings (src/ffi.cpp)
    fn ok_iterate_fasta(
        path: *const c_char,
        cb: extern "C" fn(*mut c_void, *const c_char, u64, *const c_char, u64),
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    fn ok_iterate_fasta_all(
        path: *const c_char,
        cb: extern "C" fn(*mut c_void, *const c_char, u64, *const c_char, u64),
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    fn ok_iterate_reads(
        path: *const c_char,
        cb: extern "C" fn(*mut c_void, *const c_char, u64, *const c_char, u64, *const c_char, u64),
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    fn ok_set_num_threads(n: u32);
    fn ok_set_count_limits(max_mem: u64, tmp_dir: *const c_char, graph_mode: c_int);

    // FM k-mer index (src/ffi_fm.cpp, include/fm_index.hpp)
    fn okfm_build_from_fasta(
        fasta_path: *const c_char,
        k: u64,
        min_count: u64,
        nthreads: u32,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkFmIndex;
    fn okfm_load(path: *const c_char, errbuf: *mut c_char, errbuf_len: u64) -> *mut OkFmIndex;
    fn okfm_save(
        h: *const OkFmIndex,
        path: *const c_char,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    fn okfm_free(h: *mut OkFmIndex);
    fn okfm_k(h: *const OkFmIndex) -> u32;
    fn okfm_size(h: *const OkFmIndex) -> u64;
    fn okfm_num_unitigs(h: *const OkFmIndex) -> u64;
    fn okfm_super_len(h: *const OkFmIndex) -> u64;
    fn okfm_bwt_len(h: *const OkFmIndex) -> u64;
    fn okfm_bytes(h: *const OkFmIndex) -> u64;
    fn okfm_bits_per_entry(h: *const OkFmIndex) -> f64;
    fn okfm_scan_seq(
        h: *const OkFmIndex,
        seq: *const c_char,
        slen: u64,
        k: u64,
        num_valid: *mut u64,
        num_positive: *mut u64,
        num_searches: *mut u64,
    );

    // Compacted de Bruijn graph (src/ffi_dbg.cpp)
    fn okdbg_build_from_fasta(
        fasta_path: *const c_char,
        k: u64,
        min_count: u64,
        nthreads: u32,
        keep_nodes: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkDbg;
    fn okdbg_save(h: *const OkDbg, path: *const c_char) -> c_int;
    fn okdbg_free(h: *mut OkDbg);
    fn okdbg_core_size(h: *const OkDbg) -> u64;
    fn okdbg_write_gfa(
        h: *mut OkDbg,
        path: *const c_char,
        nthreads: u32,
        st: *mut GfaStats,
        fasta: c_int,
    ) -> c_int;

    // SPSS index, static and with edit sets (src/ffi_spss.cpp)
    fn okspss_build_from_fasta(
        path: *const c_char,
        k: u64,
        m: u64,
        compact: c_int,
        min_count: u64,
        nthreads: u32,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkSpss;
    fn okspss_free(h: *mut OkSpss);
    fn okspss_save(h: *const OkSpss, path: *const c_char) -> c_int;
    fn okspss_load(
        path: *const c_char,
        verify: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkSpss;
    fn okspss_kmers(h: *const OkSpss) -> u64;
    fn okspss_bytes(h: *const OkSpss) -> u64;
    fn okspss_k(h: *const OkSpss) -> u32;
    fn okspss_m(h: *const OkSpss) -> u32;
    fn okspss_superkmers(h: *const OkSpss) -> u64;
    fn okspss_distinct(h: *const OkSpss) -> u64;
    fn okspss_ids(h: *const OkSpss, seq: *const c_char, len: u64, out: *mut i64) -> u64;
    fn okspss_scan_seq(
        h: *const OkSpss,
        seq: *const c_char,
        len: u64,
        num_valid: *mut u64,
        num_positive: *mut u64,
        num_searches: *mut u64,
        num_extensions: *mut u64,
    );
    fn okdyn_open(path: *const c_char, errbuf: *mut c_char, errbuf_len: u64) -> *mut OkDynSpss;
    fn okdyn_free(h: *mut OkDynSpss);
    fn okdyn_apply_fasta(
        h: *mut OkDynSpss,
        path: *const c_char,
        add: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> i64;
    fn okdyn_added(h: *const OkDynSpss) -> u64;
    fn okdyn_deleted(h: *const OkDynSpss) -> u64;
    fn okdyn_compact_save(h: *mut OkDynSpss, path: *const c_char, nthreads: u32) -> c_int;
    fn okdyn_kmers(h: *const OkDynSpss) -> u64;
    fn okspss_scan_many(
        h: *const OkSpss,
        seqs: *const *const c_char,
        lens: *const u64,
        n: u64,
        num_valid: *mut u64,
        num_positive: *mut u64,
        num_searches: *mut u64,
        num_extensions: *mut u64,
    );
}

/// Summary returned by `okdbg_write_gfa`, mirroring OkDbgGfaStats.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct GfaStats {
    pub n_unitigs: u64,
    pub total_bp: u64,
    pub n_kmers: u64,
    pub max_len: u64,
    pub min_len: u64,
    pub n50: u64,
    pub n_links: u64,
    pub n_dangling: u64,
}

// Read mapper kernels (src/ffi_map.cpp)
#[allow(improper_ctypes)]
unsafe extern "C" {
    fn okseed_build(
        seqs: *const *const c_char,
        lens: *const u64,
        n: u64,
        k: u32,
        w: u32,
        max_occ: u32,
        nthreads: u32,
    ) -> *mut OkSeedIdx;
    fn okseed_free(h: *mut OkSeedIdx);
    fn okseed_distinct(h: *const OkSeedIdx) -> u64;
    fn okseed_bytes(h: *const OkSeedIdx) -> u64;
    fn okseed_occurrences(h: *const OkSeedIdx) -> u64;
    fn okseed_collect(
        h: *const OkSeedIdx,
        seq: *const c_char,
        len: u64,
        out: *mut *const SeedHit,
        n_min: *mut u32,
        n_rep: *mut u32,
    ) -> u64;
    fn okseed_collect_chain(
        h: *const OkSeedIdx,
        seq: *const c_char,
        len: u64,
        out: *mut *const SeedHit,
        rep_len: *mut u32,
    ) -> u64;
    fn okseed_mid_occ(h: *const OkSeedIdx) -> u32;
    fn ok_wfa_affine_align(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        tef: i32,
        x: i32,
        o: i32,
        e: i32,
        max_score: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
        out_tend: *mut i32,
    ) -> i32;
    fn ok_biwfa_align(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        x: i32,
        o: i32,
        e: i32,
        max_score: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
    ) -> i32;
    fn ok_placements_x4(
        query: *const c_char,
        qlen: i32,
        texts: *const *const c_char,
        tlens: *const i32,
        cnt: i32,
        max_ed: i32,
        sep: i32,
        out_best: *mut i32,
        out_end: *mut i32,
        out_alt: *mut i32,
        out_alt_end: *mut i32,
    );
    fn ok_wfa_affine_cigar(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        x: i32,
        o: i32,
        e: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
    ) -> i32;
}

/// Maximum supported k-mer size, since k-mers are packed 2 bits per base into a 64-bit word.
pub(crate) const MAX_K: usize = 32;

/// Copies the NUL-terminated message the C++ side wrote into `buf`, or `default` if there is none.
fn ffi_error(buf: &[u8], default: &str) -> String {
    CStr::from_bytes_until_nul(buf)
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|_| default.into())
}

/// Makes the OpenMP runtime follow the thread budget `threads`.
///
/// libgomp reads OMP_NUM_THREADS and OMP_WAIT_POLICY once, in its load-time constructor before main, so setting them here would change nothing.
/// When they do not already say what -t needs, the process re-executes itself once with them set, keeping the same arguments, standard streams and the rest of the environment.
/// OMP_NUM_THREADS = t becomes the default of a region entered from any thread.
/// OMP_WAIT_POLICY = PASSIVE (unless the user set one) makes idle OpenMP threads sleep immediately instead of spinning.
/// By default they spin after a region, and a phase starting its own threads then ran beside them, past the budget, as bench/thread_check.py found.
/// OPTIK_OMP_ENV marks the re-executed process, so it never loops.
fn omp_env_for_budget(threads: usize) {
    use std::os::unix::process::CommandExt;
    let want_n = threads.to_string();
    let n_ok = std::env::var("OMP_NUM_THREADS").ok().as_deref() == Some(want_n.as_str());
    let wait_ok = std::env::var_os("OMP_WAIT_POLICY").is_some();
    if (n_ok && wait_ok) || std::env::var_os("OPTIK_OMP_ENV").is_some() {
        return;
    }
    let Ok(exe) = std::env::current_exe() else {
        return;
    };
    let mut cmd = std::process::Command::new(exe);
    cmd.args(std::env::args_os().skip(1))
        .env("OMP_NUM_THREADS", &want_n)
        .env("OPTIK_OMP_ENV", "1");
    if !wait_ok {
        cmd.env("OMP_WAIT_POLICY", "PASSIVE");
    }
    // exec replaces this process under the same pid and returns only on failure.
    // The run then goes on as it is, with omp_set_num_threads still covering this thread.
    let err = cmd.exec();
    eprintln!(
        "[optik] warning: could not re-execute with the OpenMP environment ({err}); continuing"
    );
}

/// Default thread budget, half the logical CPUs and at least 1.
/// -t is a budget for everything a command keeps running at once, i.e. worker pools, input readers (a .gz is decoded inside its reader's thread, with no helper process) and the C++ core's OpenMP regions, which bench/thread_check.py measures.
/// The C++ core resolves 0 the same way (thread_budget.hpp).
fn default_threads() -> usize {
    let cpus = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1);
    (cpus / 2).max(1)
}

/// FM k-mer index over the separator-delimited unitigs, after the FMSI concept of Sladky et al. and implemented in include/fm_index.hpp.
/// A separator after every unitig means no length-k window spans a boundary, so membership is an FM-index backward search that finds the k-mer iff the SA interval is non-empty.
/// Exact.
pub(crate) struct FmIdx {
    ptr: *mut OkFmIndex,
}
unsafe impl Send for FmIdx {}
unsafe impl Sync for FmIdx {}

impl FmIdx {
    /// Builds the index from the k-mers of a FASTA/FASTQ file, keeping those seen at least `min_count` times, so 1 keeps every k-mer.
    pub(crate) fn build_from_fasta(
        fasta_path: &str,
        k: usize,
        min_count: u64,
        nthreads: u32,
    ) -> Result<Self, String> {
        let cpath = CString::new(fasta_path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        let ptr = unsafe {
            okfm_build_from_fasta(
                cpath.as_ptr(),
                k as u64,
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

    /// Deserialize an index from a file.
    pub(crate) fn load(path: &str) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
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

    /// Serialize the index to a file.
    pub(crate) fn save(&self, path: &str) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        let rc = unsafe {
            okfm_save(
                self.ptr,
                cpath.as_ptr(),
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };
        if rc != 0 {
            return Err(ffi_error(&errbuf, "unknown save error"));
        }
        Ok(())
    }

    /// Looks up every k-mer of `chunk` in one C++ call, returning (valid, positive, searches, extensions), where the FM index has no extensions.
    pub(crate) fn scan_chunk(&self, chunk: &[u8], k: usize) -> (u64, u64, u64, u64) {
        let (mut nv, mut np, mut ns) = (0u64, 0u64, 0u64);
        unsafe {
            okfm_scan_seq(
                self.ptr,
                chunk.as_ptr() as *const c_char,
                chunk.len() as u64,
                k as u64,
                &mut nv,
                &mut np,
                &mut ns,
            );
        }
        (nv, np, ns, 0)
    }

    /// K-mer length used when building.
    pub(crate) fn k(&self) -> usize {
        unsafe { okfm_k(self.ptr) as usize }
    }
    /// Number of distinct canonical k-mers indexed.
    pub(crate) fn size(&self) -> usize {
        unsafe { okfm_size(self.ptr) as usize }
    }
    /// Number of maximal unitigs indexed.
    pub(crate) fn num_unitigs(&self) -> usize {
        unsafe { okfm_num_unitigs(self.ptr) as usize }
    }
    /// Total bases of the unitigs, separators excluded.
    pub(crate) fn super_len(&self) -> usize {
        unsafe { okfm_super_len(self.ptr) as usize }
    }
    /// Length of the BWT (superstring + separators + sentinel).
    pub(crate) fn bwt_len(&self) -> usize {
        unsafe { okfm_bwt_len(self.ptr) as usize }
    }
    /// Size of the saved index file, in bytes.
    pub(crate) fn bytes(&self) -> usize {
        unsafe { okfm_bytes(self.ptr) as usize }
    }
    /// Average bits per stored k-mer.
    pub(crate) fn bits_per_entry(&self) -> f64 {
        unsafe { okfm_bits_per_entry(self.ptr) }
    }
}
impl Drop for FmIdx {
    fn drop(&mut self) {
        unsafe { okfm_free(self.ptr) }
    }
}

/// Compacted de Bruijn graph built by the k-mer counter (src/ffi_dbg.cpp).
pub(crate) struct DbgIdx {
    ptr: *mut OkDbg,
}
unsafe impl Send for DbgIdx {}
unsafe impl Sync for DbgIdx {}

impl DbgIdx {
    /// Counts the canonical k-mers of `fasta_path` seen at least `min_count` times into a graph on `nthreads` threads, as src/kmer_count.cpp does.
    /// `keep_nodes` keeps the node table that the `-s` archive needs.
    /// Otherwise the unitigs may come from partitioned compaction.
    pub(crate) fn build_from_fasta(
        fasta_path: &str,
        k: usize,
        min_count: u64,
        nthreads: usize,
        keep_nodes: bool,
    ) -> Result<Self, String> {
        let cpath = CString::new(fasta_path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        let ptr = unsafe {
            okdbg_build_from_fasta(
                cpath.as_ptr(),
                k as u64,
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

    /// Serialize the graph core to a kdbg archive.
    pub(crate) fn save(&self, path: &str) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        if unsafe { okdbg_save(self.ptr, cpath.as_ptr()) } != 0 {
            return Err("dbg save failed".into());
        }
        Ok(())
    }

    /// Distinct canonical k-mers in the graph.
    pub(crate) fn size(&self) -> usize {
        unsafe { okdbg_core_size(self.ptr) as usize }
    }

    /// Writes the maximal unitigs as a GFA1 graph of S and L lines, or as unitig FASTA when `fasta` is set, and returns the summary statistics.
    pub(crate) fn write_gfa(
        &self,
        path: &str,
        nthreads: usize,
        fasta: bool,
    ) -> Result<GfaStats, String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut st = GfaStats::default();
        let rc = unsafe {
            okdbg_write_gfa(
                self.ptr,
                cpath.as_ptr(),
                nthreads.max(1) as u32,
                &mut st,
                fasta as c_int,
            )
        };
        if rc != 0 {
            return Err(format!(
                "failed to write {} to {path}",
                if fasta { "FASTA" } else { "GFA" }
            ));
        }
        Ok(st)
    }
}
impl Drop for DbgIdx {
    fn drop(&mut self) {
        unsafe { okdbg_free(self.ptr) }
    }
}

/// Static SPSS core, exact k-mer membership over super-k-mer strings with a minimizer MPHF (spss_core.hpp).
pub(crate) struct SpssIdx {
    ptr: *mut OkSpss,
}
unsafe impl Send for SpssIdx {}
unsafe impl Sync for SpssIdx {}

impl SpssIdx {
    /// Builds from FASTA/FASTQ.
    /// `compact` counts the k-mers seen at least `min_count` times and indexes their maximal unitigs.
    /// Otherwise the input sequences are indexed as is.
    pub(crate) fn build_from_fasta(
        path: &str,
        k: usize,
        m: usize,
        compact: bool,
        min_count: u64,
        nthreads: usize,
    ) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        let ptr = unsafe {
            okspss_build_from_fasta(
                cpath.as_ptr(),
                k as u64,
                m as u64,
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
    /// Map an index, where `verify` first checks its checksum by reading the file once.
    pub(crate) fn load(path: &str, verify: bool) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
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
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        if unsafe { okspss_save(self.ptr, cpath.as_ptr()) } != 0 {
            return Err(format!("failed to write {path}"));
        }
        Ok(())
    }
    /// Indexed k-mer positions (a repeated k-mer counts once per position).
    pub(crate) fn kmers(&self) -> u64 {
        unsafe { okspss_kmers(self.ptr) }
    }
    /// In-memory size of the index, in bytes.
    pub(crate) fn bytes(&self) -> u64 {
        unsafe { okspss_bytes(self.ptr) }
    }
    pub(crate) fn k(&self) -> usize {
        unsafe { okspss_k(self.ptr) as usize }
    }
    pub(crate) fn m(&self) -> usize {
        unsafe { okspss_m(self.ptr) as usize }
    }
    pub(crate) fn superkmers(&self) -> u64 {
        unsafe { okspss_superkmers(self.ptr) }
    }
    /// Distinct k-mers (kmers() counts positions, which repeat duplicates).
    pub(crate) fn distinct(&self) -> u64 {
        unsafe { okspss_distinct(self.ptr) }
    }
    /// Dense k-mer identifier of every window of `seq`, -1 when absent or N.
    pub(crate) fn ids(&self, seq: &[u8]) -> Vec<i64> {
        if seq.len() < self.k() {
            return Vec::new();
        }
        let mut out = vec![0i64; seq.len() - self.k() + 1];
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
        let (mut v, mut p, mut s, mut e) = (0u64, 0u64, 0u64, 0u64);
        unsafe {
            okspss_scan_seq(
                self.ptr,
                chunk.as_ptr() as *const c_char,
                chunk.len() as u64,
                &mut v,
                &mut p,
                &mut s,
                &mut e,
            );
        }
        (v, p, s, e)
    }

    /// Streaming lookup over a batch of records with cache misses overlapping across records, returning (valid, positive, searches, extensions).
    pub(crate) fn scan_many(&self, recs: &[&[u8]]) -> (u64, u64, u64, u64) {
        let ptrs: Vec<*const c_char> = recs.iter().map(|r| r.as_ptr() as *const c_char).collect();
        let lens: Vec<u64> = recs.iter().map(|r| r.len() as u64).collect();
        let (mut v, mut p, mut s, mut e) = (0u64, 0u64, 0u64, 0u64);
        unsafe {
            okspss_scan_many(
                self.ptr,
                ptrs.as_ptr(),
                lens.as_ptr(),
                recs.len() as u64,
                &mut v,
                &mut p,
                &mut s,
                &mut e,
            );
        }
        (v, p, s, e)
    }
}
impl Drop for SpssIdx {
    fn drop(&mut self) {
        unsafe { okspss_free(self.ptr) }
    }
}

/// SPSS index opened for editing, i.e. the static core plus concurrent open-addressing edit sets (include/dyn_spss.hpp).
pub(crate) struct DynSpssIdx {
    ptr: *mut OkDynSpss,
}

impl DynSpssIdx {
    /// Loads the index at `path`, checksum verified.
    pub(crate) fn open(path: &str) -> Result<Self, String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
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
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        let mut errbuf = vec![0u8; ERRBUF_LEN];
        let n = unsafe {
            okdyn_apply_fasta(
                self.ptr,
                cpath.as_ptr(),
                add as c_int,
                errbuf.as_mut_ptr() as *mut c_char,
                errbuf.len() as u64,
            )
        };
        if n < 0 {
            return Err(ffi_error(&errbuf, "failed to read k-mers"));
        }
        Ok(n as u64)
    }
    /// Size of the insertion set.
    pub(crate) fn added(&self) -> u64 {
        unsafe { okdyn_added(self.ptr) }
    }
    /// Size of the erasure set.
    pub(crate) fn deleted(&self) -> u64 {
        unsafe { okdyn_deleted(self.ptr) }
    }
    /// k-mer positions of the static core.
    pub(crate) fn kmers(&self) -> u64 {
        unsafe { okdyn_kmers(self.ptr) }
    }
    /// Compacts the edits into the static core and saves the result to `path`.
    pub(crate) fn compact_save(&mut self, path: &str, nthreads: usize) -> Result<(), String> {
        let cpath = CString::new(path).map_err(|e| e.to_string())?;
        if unsafe { okdyn_compact_save(self.ptr, cpath.as_ptr(), nthreads.max(1) as u32) } != 0 {
            return Err(format!("failed to compact/write {path}"));
        }
        Ok(())
    }
}
impl Drop for DynSpssIdx {
    fn drop(&mut self) {
        unsafe { okdyn_free(self.ptr) }
    }
}

// The first four bytes of each index file, read as a little-endian u32 (include/fm_index.hpp and src/ffi_spss.cpp write them).

pub(crate) const FM_MAGIC: u32 = u32::from_le_bytes(*b"OKFM"); // k-order BWT FM index
pub(crate) const SPSS_MAGIC: u32 = u32::from_le_bytes(*b"OKSP"); // SPSS index, "OKSPSS" padded to 8 bytes

/// Returns the first 4 bytes of an index file, its format magic, or None if the file is shorter than that.
/// An unreadable file is an error naming the path.
pub(crate) fn index_magic(path: &str) -> Result<Option<u32>, String> {
    use std::fs::File;
    use std::io::Read;
    let mut f = File::open(path).map_err(|e| format!("cannot open index {path}: {e}"))?;
    let mut b = [0u8; 4];
    Ok(f.read_exact(&mut b).ok().map(|_| u32::from_le_bytes(b)))
}

/// Calls `f(name, seq)` for every record of a FASTA/FASTQ file with a non-empty sequence, read by include/fasta_reader.hpp with line breaks inside a sequence removed.
pub(crate) fn iterate_fasta<F: FnMut(&str, &[u8])>(path: &str, f: F) -> Result<(), String> {
    iterate_fasta_impl(path, f, false)
}

/// As [`iterate_fasta`], but records with an empty sequence are reported too.
pub(crate) fn iterate_fasta_all<F: FnMut(&str, &[u8])>(path: &str, f: F) -> Result<(), String> {
    iterate_fasta_impl(path, f, true)
}

/// Calls `f(name, seq, qual)` for every record, empty ones included.
/// `qual` is the FASTQ quality string, empty for FASTA.
pub(crate) fn iterate_reads<F: FnMut(&str, &[u8], &[u8])>(
    path: &str,
    mut f: F,
) -> Result<(), String> {
    extern "C" fn cb<F: FnMut(&str, &[u8], &[u8])>(
        cookie: *mut c_void,
        name: *const c_char,
        nlen: u64,
        seq: *const c_char,
        slen: u64,
        qual: *const c_char,
        qlen: u64,
    ) {
        let f = unsafe { &mut *(cookie as *mut F) };
        let name_s = unsafe {
            std::str::from_utf8_unchecked(std::slice::from_raw_parts(
                name as *const u8,
                nlen as usize,
            ))
        };
        let seq_b = unsafe { std::slice::from_raw_parts(seq as *const u8, slen as usize) };
        let qual_b: &[u8] = if qlen == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(qual as *const u8, qlen as usize) }
        };
        f(name_s, seq_b, qual_b);
    }
    let cpath = CString::new(path).map_err(|e| e.to_string())?;
    let mut errbuf = vec![0i8; READER_ERRBUF_LEN];
    let rc = unsafe {
        ok_iterate_reads(
            cpath.as_ptr(),
            cb::<F>,
            &mut f as *mut F as *mut c_void,
            errbuf.as_mut_ptr(),
            errbuf.len() as u64,
        )
    };
    if rc != 0 {
        let msg = unsafe {
            std::ffi::CStr::from_ptr(errbuf.as_ptr())
                .to_string_lossy()
                .into_owned()
        };
        return Err(msg);
    }
    Ok(())
}

/// Shared body of [`iterate_fasta`] and [`iterate_fasta_all`].
fn iterate_fasta_impl<F: FnMut(&str, &[u8])>(
    path: &str,
    mut f: F,
    keep_empty: bool,
) -> Result<(), String> {
    extern "C" fn cb<F: FnMut(&str, &[u8])>(
        cookie: *mut c_void,
        name: *const c_char,
        nlen: u64,
        seq: *const c_char,
        slen: u64,
    ) {
        let f = unsafe { &mut *(cookie as *mut F) };
        let name_s = unsafe {
            std::str::from_utf8_unchecked(std::slice::from_raw_parts(
                name as *const u8,
                nlen as usize,
            ))
        };
        let seq_b = unsafe { std::slice::from_raw_parts(seq as *const u8, slen as usize) };
        f(name_s, seq_b);
    }
    let cpath = CString::new(path).map_err(|e| e.to_string())?;
    let mut errbuf = vec![0i8; READER_ERRBUF_LEN];
    let rc = unsafe {
        let iter = if keep_empty {
            ok_iterate_fasta_all
        } else {
            ok_iterate_fasta
        };
        iter(
            cpath.as_ptr(),
            cb::<F>,
            &mut f as *mut F as *mut c_void,
            errbuf.as_mut_ptr(),
            errbuf.len() as u64,
        )
    };
    if rc != 0 {
        let msg = unsafe {
            std::ffi::CStr::from_ptr(errbuf.as_ptr())
                .to_string_lossy()
                .into_owned()
        };
        return Err(msg);
    }
    Ok(())
}

/// Mapper seed index, a minimizer MPHF with CSR postings and the reference as key store (seed_index.hpp).
/// It borrows the reference sequences, which must outlive it.
pub(crate) struct SeedIdx {
    ptr: *mut OkSeedIdx,
}
unsafe impl Send for SeedIdx {}
unsafe impl Sync for SeedIdx {}
impl SeedIdx {
    /// Builds the seed index over `refs` with minimizer length `k`, window `w` and repeat cutoff `max_occ` (include/seed_index.hpp).
    pub(crate) fn build(
        refs: &[Vec<u8>],
        k: usize,
        w: usize,
        max_occ: u32,
        threads: usize,
    ) -> Result<Self, String> {
        // Seed postings pack a position in 31 bits, so the limit is named rather than letting the run fail opaquely, or wrap positions as it did before this guard.
        if let Some((i, r)) = refs
            .iter()
            .enumerate()
            .find(|(_, r)| r.len() >= (1usize << 31))
        {
            return Err(format!(
                "reference sequence #{} is {} bases; map supports sequences \
                                shorter than 2^31 (2,147,483,648) bases -- split it",
                i + 1,
                r.len()
            ));
        }
        let ptrs: Vec<*const c_char> = refs.iter().map(|r| r.as_ptr() as *const c_char).collect();
        let lens: Vec<u64> = refs.iter().map(|r| r.len() as u64).collect();
        let ptr = unsafe {
            okseed_build(
                ptrs.as_ptr(),
                lens.as_ptr(),
                refs.len() as u64,
                k as u32,
                w as u32,
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
        Ok(SeedIdx { ptr })
    }
    /// Distinct minimizers indexed.
    pub(crate) fn distinct(&self) -> u64 {
        unsafe { okseed_distinct(self.ptr) }
    }
    /// In-memory size of the index, in bytes.
    pub(crate) fn bytes(&self) -> u64 {
        unsafe { okseed_bytes(self.ptr) }
    }
    /// Fraction of the reference's minimizer occurrences that repeat an earlier one (1 - distinct / occurrences), 0 for a genome without repeats.
    /// A genome-level MAPQ input (map.rs, hidden_odds).
    pub(crate) fn repeat_fraction(&self) -> f64 {
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
        let mut p: *const SeedHit = std::ptr::null();
        let (mut n_min, mut n_rep) = (0u32, 0u32);
        let n = unsafe {
            okseed_collect(
                self.ptr,
                seq.as_ptr() as *const c_char,
                seq.len() as u64,
                &mut p,
                &mut n_min,
                &mut n_rep,
            )
        };
        let hits = if n == 0 || p.is_null() {
            Vec::new()
        } else {
            unsafe { std::slice::from_raw_parts(p, n as usize) }.to_vec()
        };
        (hits, n_min, n_rep)
    }
    /// The chaining path's occurrence threshold (seed_index.hpp mid_occ).
    pub(crate) fn mid_occ(&self) -> u32 {
        unsafe { okseed_mid_occ(self.ptr) }
    }
    /// Seed hits for chaining (repeat runs keep their least frequent minimizers) and how many read bases lie only under dropped minimizers.
    pub(crate) fn collect_chain(&self, seq: &[u8]) -> (Vec<SeedHit>, u32) {
        let mut p: *const SeedHit = std::ptr::null();
        let mut rep_len = 0u32;
        let n = unsafe {
            okseed_collect_chain(
                self.ptr,
                seq.as_ptr() as *const c_char,
                seq.len() as u64,
                &mut p,
                &mut rep_len,
            )
        };
        let hits = if n == 0 || p.is_null() {
            Vec::new()
        } else {
            unsafe { std::slice::from_raw_parts(p, n as usize) }.to_vec()
        };
        (hits, rep_len)
    }
}
impl Drop for SeedIdx {
    fn drop(&mut self) {
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
) -> [Placements; 4] {
    debug_assert!(texts.len() <= PLACEMENT_LANES);
    let ptrs: Vec<*const c_char> = texts.iter().map(|t| t.as_ptr() as *const c_char).collect();
    let lens: Vec<i32> = texts.iter().map(|t| t.len() as i32).collect();
    let (mut b, mut e, mut a, mut ae) = ([max_ed + 1; 4], [0i32; 4], [max_ed + 1; 4], [0i32; 4]);
    unsafe {
        ok_placements_x4(
            query.as_ptr() as *const c_char,
            query.len() as i32,
            ptrs.as_ptr(),
            lens.as_ptr(),
            texts.len().min(PLACEMENT_LANES) as i32,
            max_ed,
            sep,
            b.as_mut_ptr(),
            e.as_mut_ptr(),
            a.as_mut_ptr(),
            ae.as_mut_ptr(),
        );
    }
    let mut out = [Placements {
        best: max_ed + 1,
        end: 0,
        alt: max_ed + 1,
        alt_end: 0,
    }; 4];
    for i in 0..PLACEMENT_LANES {
        out[i] = Placements {
            best: b[i],
            end: e[i],
            alt: a[i],
            alt_end: ae[i],
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

/// A buffer large enough for the CIGAR of aligning `qlen` against `tlen` bases, i.e. at most one operation per base, each at most 8 bytes, plus `slack`.
fn cigar_buffer(qlen: usize, tlen: usize, slack: usize) -> Vec<u8> {
    vec![0u8; slack + 8 * (qlen + tlen)]
}

/// Gap-affine alignment with the query end to end from text position 0, while the text may end anywhere in [tlen - tef, tlen].  max_score < 0 is unbounded.
pub(crate) fn wfa_affine_align(
    query: &[u8],
    text: &[u8],
    tef: usize,
    pen: (i32, i32, i32),
    max_score: i32,
) -> Option<AffineAln> {
    let mut buf = cigar_buffer(query.len(), text.len(), 32);
    let (mut sc, mut te) = (0i32, 0i32);
    let st = unsafe {
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
            &mut sc,
            &mut te,
        )
    };
    if st != 0 {
        return None;
    }
    let cigar = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();
    Some(AffineAln {
        cigar,
        score: sc,
        tend: te as usize,
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
    let mut buf = cigar_buffer(query.len(), text.len(), 32);
    let mut sc = 0i32;
    let st = unsafe {
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
            &mut sc,
        )
    };
    if st != 0 {
        return None;
    }
    let cigar = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();
    Some(AffineAln {
        cigar,
        score: sc,
        tend: text.len(),
    })
}

/// Penalties of [`wfa_affine_cigar`], as mismatch, gap open and gap extend.
const CIGAR_PENALTIES: (i32, i32, i32) = (4, 6, 2);

/// Gap-affine semi-global CIGAR (include/wfa_affine.hpp) with [`CIGAR_PENALTIES`], returned as (cigar, 0-based text start), or None on failure.
pub(crate) fn wfa_affine_cigar(query: &[u8], text: &[u8]) -> Option<(String, usize)> {
    let mut buf = cigar_buffer(query.len(), text.len(), 16);
    let mut score = 0i32;
    let st = unsafe {
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
    if st < 0 {
        return None;
    }
    let cig = CStr::from_bytes_until_nul(&buf)
        .ok()?
        .to_string_lossy()
        .into_owned();
    Some((cig, st as usize))
}

/// Peak resident-set growth since creation.
///
/// Reads VmRSS from /proc/self/status at creation and at each `sample()`, so the difference is the memory an operation added, without the process's baseline (libraries, runtime).
pub(crate) struct RssWatermark {
    before_kb: u64,
    peak_kb: u64,
}

impl RssWatermark {
    /// Starts tracking from the current resident set.
    pub(crate) fn new() -> Self {
        let r = Self::read_kb();
        RssWatermark {
            before_kb: r,
            peak_kb: r,
        }
    }

    /// Record the current RSS, updating the peak if it increased.
    pub(crate) fn sample(&mut self) {
        let r = Self::read_kb();
        if r > self.peak_kb {
            self.peak_kb = r;
        }
    }

    /// Net memory increase from start to peak, in megabytes.
    pub(crate) fn net_mb(&self) -> f64 {
        self.peak_kb.saturating_sub(self.before_kb) as f64 / 1024.0
    }

    /// Read the current RSS from the kernel's proc filesystem.
    fn read_kb() -> u64 {
        std::fs::read_to_string("/proc/self/status")
            .unwrap_or_default()
            .lines()
            .find(|l| l.starts_with("VmRSS:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    }
}

/// Parses a memory size, a number with an optional K/M/G/T suffix (powers of 1024, where a bare number is bytes).
/// None if malformed or not positive.
fn parse_size(v: &str) -> Option<u64> {
    // 8G, 8GB, 8GiB and 8g are the same size
    let t = v.trim();
    let t = t
        .strip_suffix("iB")
        .or_else(|| t.strip_suffix("ib"))
        .unwrap_or(t);
    let t = if t.len() > 1 && t[..t.len() - 1].ends_with(|c: char| "kmgtKMGT".contains(c)) {
        t.strip_suffix('B')
            .or_else(|| t.strip_suffix('b'))
            .unwrap_or(t)
    } else {
        t
    };
    let (num, mul) = match t.chars().last()?.to_ascii_uppercase() {
        'K' => (&t[..t.len() - 1], 1u64 << 10),
        'M' => (&t[..t.len() - 1], 1u64 << 20),
        'G' => (&t[..t.len() - 1], 1u64 << 30),
        'T' => (&t[..t.len() - 1], 1u64 << 40),
        _ => (t, 1u64),
    };
    let x: f64 = num.parse().ok()?;
    if !x.is_finite() || x <= 0.0 {
        return None;
    }
    let b = x * mul as f64;
    if b >= u64::MAX as f64 {
        return None;
    }
    Some(b as u64)
}

/// Whether two paths name the same existing file, by device and inode.
fn same_file(a: &str, b: &str) -> bool {
    use std::os::unix::fs::MetadataExt;
    match (std::fs::metadata(a), std::fs::metadata(b)) {
        (Ok(x), Ok(y)) => x.dev() == y.dev() && x.ino() == y.ino(),
        _ => false,
    }
}

/// The value flags and the switch flags subcommand `sub` accepts, or None for an unknown subcommand.
/// The TUI's forms are checked against this table.
pub(crate) fn subcommand_flags(
    sub: &str,
) -> Option<(&'static [&'static str], &'static [&'static str])> {
    let flags: (&'static [&'static str], &'static [&'static str]) = match sub {
        "map" => (
            &["-r", "-q", "-k", "-w", "-e", "-o", "-t", "-R"],
            &["--cs", "--cs=short", "--cs=long"],
        ),
        "dbruijn" => (
            &[
                "-i",
                "-k",
                "-c",
                "-t",
                "-o",
                "-s",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &[],
        ),
        "spss-update" => (&["-i", "--add", "--remove", "-t", "-o"], &[]),
        "build-spss" => (
            &[
                "-i",
                "-k",
                "-m",
                "-c",
                "-t",
                "-o",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &["--compact"],
        ),
        "build-fm" => (
            &[
                "-i",
                "-k",
                "-c",
                "-t",
                "-o",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &[],
        ),
        "query" => (&["-i", "-q", "-k", "-t", "--ids"], &["--no-verify"]),
        _ => return None,
    };
    Some(flags)
}

/// One subcommand's help: a one-line summary, the synopsis, its FLAGS block and examples.
struct CommandHelp {
    name: &'static str,
    summary: &'static str,
    synopsis: &'static str,
    flags: String,
    examples: &'static str,
}

/// The help of every subcommand, in the order the overview lists them.
fn command_helps() -> Vec<CommandHelp> {
    let mk = crate::MAX_K;
    let dk = crate::MAX_K - 1; // dbruijn's (k+1)-mer edge set must fit 64 bits
    vec![
        CommandHelp {
            name: "map",
            summary: "Map short, paired-end or long reads to a reference, writing SAM",
            synopsis: "optik map -r REF.fa -q READS.fq [-q MATES.fq] [options] [-o OUT.sam]",
            flags: "FLAGS (map)
  -r FILE    Reference genome in FASTA format
  -q FILE    Sequencing reads in FASTA or FASTQ format; give -q twice for
               paired-end reads (-q R1.fq -q R2.fq, same order, FR library):
               mate rescue, proper-pair flags, mate fields, TLEN, MC:Z
  -k INT     Seed (minimizer) length for map  [default: chosen from the
               reads: 21, or 15 for noisy long reads (over ~2% error, from
               base qualities, or from the share of read minimizers found
               in the reference for FASTA reads); -k fixes it]
  -w INT     Sliding window for minimizer selection (1..256, 0 = auto-tune
               from the read length)  [default: 10]
               Smaller window = more minimizer seeds, higher sensitivity
               Reads shorter than k+w-1 bases produce NO minimizers (all unmapped)
  -e INT     Maximum edit distance for short-read alignment (>= 0)
               [default: auto = 15% of the read length, at least 5]
               Higher = more tolerant of mismatches and indels
  -t INT     Number of threads for parallel mapping (0 = auto)  [default: half of logical CPUs]
  -o FILE    Write SAM to FILE  [default: standard output]
  -R STR     Read group line, e.g. '@RG\\tID:s1\\tSM:sample1' (as bwa and
               minimap2): written to the header, and RG:Z:<ID> added to
               every record (GATK / Picard require read groups)
  --cs       Also write minimap2's cs tag (short form, `:n` runs); --cs=long
               spells identical bases (`=ACGT`).  Every record carries
               NM and MD (SAM spec) either way.

NOTES (map)
  Reads must be at least k + w - 1 bases long for minimizer extraction.
  E.g., with the default k=21 w=10 each read needs >=30 bases; shorter reads
  always appear unmapped. Pre-extracted k-mers (single k-mer per
  sequence) should use the `query` subcommand instead.
  RNA U is read as T.  Reference names must be unique (SAM @SQ); read
  names over 254 characters are truncated (SAM QNAME limit)."
                .to_string(),
            examples: "  optik map -r genome.fa -q reads.fq -o reads.sam                 # single-end or long reads\n  optik map -r genome.fa -q R1.fq.gz -q R2.fq.gz -o pairs.sam     # paired-end\n  optik map -r genome.fa -q R1.fq -q R2.fq -R '@RG\\tID:s1\\tSM:s1' --cs -o pairs.sam",
        },
        CommandHelp {
            name: "dbruijn",
            summary: "Count k-mers and write the compacted de Bruijn graph (GFA with links, or unitig FASTA)",
            synopsis: "optik dbruijn -i SEQS.fa|.fq[.gz] [-i ...] [-k 31] [-c 2] [-o graph.gfa]",
            flags: format!(
                "FLAGS (dbruijn)
  -i FILE    Input sequencing reads in FASTA or FASTQ format (.gz, - for
               stdin); repeat for several files, e.g. -i R1.fq -i R2.fq
  -k INT     k-mer size for De Bruijn graph (1..{dk})  [default: 31]
               Larger k = more specific graph spurs
  -c INT     Minimum k-mer abundance to include (0 acts as 1)  [default: 2]
               Higher = filters more sequencing errors
  -t INT     Number of threads for counting + unitig walk (0 = auto)
               [default: half of logical CPUs]
  -o FILE    Output  [default: graph.gfa]: GFA1 (segments + links), or
               unitig FASTA if FILE ends in .fa/.fasta/.fna (BCALM-style
               headers `>id LN:i KC:i km:f`, what SSHash/SBWT/FMSI index)
  -s FILE    Also serialize the graph core to a kdbg archive (not with
               --graph part, which keeps no node table)
  --max-mem SIZE  Memory limit for k-mer counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   Spill directory  [default: the output file's directory]
  --graph MODE    How the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= 31) or table.  The
               output is identical either way"
            ),
            examples: "  optik dbruijn -i genome.fa -k 31 -c 1 -o graph.gfa          # a genome: keep every k-mer\n  optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa         # unitig FASTA instead of GFA\n  optik dbruijn -i R1.fq.gz -i R2.fq.gz -k 31 -o graph.gfa    # reads: -c 2 drops error k-mers",
        },
        CommandHelp {
            name: "build-spss",
            summary: "Build an exact SPSS k-mer index (membership and dense ids)",
            synopsis: "optik build-spss -i SEQS.fa [-i ...] [-k 31] [--compact [-c N]] -o INDEX.spss",
            flags: "FLAGS (build-spss)
  -i FILE    FASTA/FASTQ to index (an SPSS such as unitigs is size-optimal;
               any input gives exact membership); repeatable
  -k INT     k-mer length (1..32)  [default: 31]
  -m INT     Minimizer length (< k)  [default: 0 = auto: chosen from the
               input by an estimate of index size and lookup work per m;
               the index is exact at every m]
  --compact  Count k-mers (>= -c) and index their maximal unitigs instead of
               the raw input sequences
  -c INT     Minimum abundance with --compact  [default: 1]
  -t INT     Threads for the build, and for counting and compaction with
               --compact (0 = auto)  [default: half of logical CPUs]
  -o FILE    Output index file (required)
  --max-mem SIZE  --compact: memory limit for counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   --compact: spill directory  [default: the output file's directory]
  --graph MODE    --compact: how the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= 31) or table.  The
               output is identical either way"
                .to_string(),
            examples: "  optik build-spss -i unitigs.fa -k 31 -o genome.spss\n  optik build-spss --compact -i R1.fq.gz -i R2.fq.gz -c 2 -o reads.spss   # count, compact, index",
        },
        CommandHelp {
            name: "spss-update",
            summary: "Add or remove the k-mers of files in an SPSS index and write it re-compacted",
            synopsis: "optik spss-update -i INDEX.spss [--add FILE ...] [--remove FILE ...] -o OUT.spss",
            flags: "FLAGS (spss-update)
  -i FILE        SPSS index to edit
  --add FILE     Insert every k-mer of FILE (repeatable)
  --remove FILE  Erase every k-mer of FILE (repeatable; applied after adds)
  -t INT         Threads for compaction (0 = auto)
                   [default: half of logical CPUs]
  -o FILE        Write the compacted index (may equal -i)"
                .to_string(),
            examples: "  optik spss-update -i genome.spss --add new.fa --remove old.fa -o genome.spss",
        },
        CommandHelp {
            name: "build-fm",
            summary: "Build an exact FM k-mer index of the unitigs (the smallest index)",
            synopsis: "optik build-fm -i SEQS.fa [-i ...] [-k 31] [-c 1] -o INDEX.fm",
            flags: format!(
                "FLAGS (build-fm)
  -i FILE    Reference FASTA/FASTQ to build the FM index from; repeatable
  -k INT     k-mer length (1..{mk})  [default: 31]
  -c INT     Minimum k-mer abundance to include  [default: 1]
               Default keeps every distinct k-mer (exact membership)
  -t INT     Threads for counting, compaction and the index build (0 = auto)
               [default: half of logical CPUs]
  -o FILE    Write serialized index to file (required)
  --max-mem SIZE  Memory limit for k-mer counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   Spill directory  [default: the output file's directory]
  --graph MODE    How the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= 31) or table.  The
               output is identical either way"
            ),
            examples: "  optik build-fm -i genome.fa -k 31 -o genome.fm",
        },
        CommandHelp {
            name: "query",
            summary: "Look up every k-mer of a FASTA/FASTQ in an SPSS or FM index",
            synopsis: "optik query -i INDEX.spss|INDEX.fm -q QUERIES.fa|.fq[.gz] [-t N] [--ids IDS.txt]",
            flags: format!(
                "FLAGS (query)
  -i FILE    Saved SPSS or FM index (build-spss / build-fm)
  -q FILE    FASTA/FASTQ whose k-mers are looked up (k-mers or reads; .gz ok)
  -k INT     k-mer length (1..{mk}); must equal the index's  [default: the index's k]
  -t INT     Threads for the query phase (1 = single-threaded,
               >1 = parallel over records, long records split into
               window-exact slices; counts are identical at any -t;
               0 = auto)  [default: half of logical CPUs]
  --ids FILE SPSS index only: also write, per query record, the dense
               identifier in [0, n) of every k-mer window (-1 = absent),
               space-separated, one line per record (SSHash Lookup)
  --no-verify  SPSS index only: skip the checksum check at load (by default
               the whole index is checksummed once, so a truncated or
               corrupted file is an error instead of wrong answers)"
            ),
            examples: "  optik query -i genome.spss -q reads.fq                  # counts of present and absent k-mers\n  optik query -i genome.spss -q reads.fq --ids ids.txt    # also one id per k-mer window",
        },
    ]
}

/// The overview `optik` prints when no command or an unknown one is given:
/// what each command does, one example each and where to read more.
fn overview_text() -> String {
    let mut t = String::from(
        "optik -- genomic sequence analysis toolkit

USAGE
  optik <command> [options]
  optik <command> --help      that command's flags and examples
  optik --help                every command and flag
  optik tui                   interactive terminal UI (also plain `optik` in a terminal)
  optik --version

COMMANDS
",
    );
    for c in command_helps() {
        t.push_str(&format!("  {:<12}{}\n", c.name, c.summary));
    }
    t.push_str("  build       Alias of build-spss\n  tui         Interactive terminal UI over every command\n\nQUICK START\n");
    t.push_str(
        "  optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa   # graph of a genome, as unitigs
  optik build-spss -i unitigs.fa -o genome.spss          # exact k-mer index of the unitigs
  optik query -i genome.spss -q reads.fq                 # look up every k-mer of the reads
  optik map -r genome.fa -q reads.fq -o reads.sam        # map the reads to SAM
",
    );
    t.push_str("\nRun `optik <command> --help` for a command's flags and more examples.\n");
    t
}

/// The help of subcommand @sub, or None for an unknown one.
fn command_help_text(sub: &str) -> Option<String> {
    let c = command_helps().into_iter().find(|c| c.name == sub)?;
    Some(format!(
        "optik {} -- {}\n\nUSAGE\n  {}\n\n{}\n\nEXAMPLES\n{}\n",
        c.name, c.summary, c.synopsis, c.flags, c.examples
    ))
}

/// Prints the help text on stderr.
fn usage() {
    eprintln!("{}", usage_text());
}

/// The full help, the overview and then every command's flags.
fn usage_text() -> String {
    let mut t = overview_text();
    t.push('\n');
    for c in command_helps() {
        t.push('\n');
        t.push_str(&c.flags);
        t.push('\n');
    }
    t
}

/// The nearest known command to @word by edit distance, for "did you mean".
fn closest_command(word: &str) -> Option<&'static str> {
    fn dist(a: &str, b: &str) -> usize {
        let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
        let mut row: Vec<usize> = (0..=b.len()).collect();
        for i in 1..=a.len() {
            let mut prev = row[0];
            row[0] = i;
            for j in 1..=b.len() {
                let cur = row[j];
                row[j] = (row[j] + 1)
                    .min(row[j - 1] + 1)
                    .min(prev + usize::from(a[i - 1] != b[j - 1]));
                prev = cur;
            }
        }
        row[b.len()]
    }
    const NAMES: [&str; 8] = [
        "map",
        "dbruijn",
        "build-spss",
        "spss-update",
        "build-fm",
        "query",
        "build",
        "tui",
    ];
    NAMES
        .iter()
        .map(|n| (dist(word, n), *n))
        .filter(|(d, n)| *d <= 2.max(n.len() / 3))
        .min()
        .map(|(_, n)| n)
}

/// Report a command-line mistake in @sub: the message, the command's synopsis and where its full help is, then exit 1.
fn arg_error(sub: &str, msg: &str) -> ! {
    eprintln!("error: {msg}");
    if let Some(c) = command_helps().into_iter().find(|c| c.name == sub) {
        eprintln!("usage: {}", c.synopsis);
    }
    eprintln!("see: optik {sub} --help");
    std::process::exit(1)
}

fn main() {
    // Output written to a closed pipe, as in `optik query ... | head`, ends the process quietly, as for C tools.
    // Rust ignores SIGPIPE by default and println! then panics, which with panic = "abort" dumped core.
    // Any other failure to write stdout, such as a full disk, is an error rather than a crash.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let msg = info
            .payload()
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| info.payload().downcast_ref::<&str>().copied())
            .unwrap_or("");
        if msg.starts_with("failed printing to stdout") {
            eprintln!(
                "error: writing to stdout failed: {}",
                msg.trim_start_matches("failed printing to stdout: ")
            );
            std::process::exit(1);
        }
        default_hook(info);
    }));
    // glibc raises its mmap threshold, up to 32 MB, each time a large mmap'd block is freed.
    // Later large buffers, such as per-thread counting and walk arrays and query batches, then land on arena heaps and stay resident after they are freed, adding about 100 MB to the peak on a single-record genome.
    // A fixed 1 MiB threshold keeps large buffers mmap'd and returned to the OS on free, at the same measured speed.
    // Only glibc does this.
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        const MMAP_THRESHOLD_BYTES: c_int = 1 << 20;
        libc::mallopt(libc::M_MMAP_THRESHOLD, MMAP_THRESHOLD_BYTES);
    }
    let args: Vec<String> = std::env::args().collect();

    // Without arguments and without a terminal (scripts, pipes) the usage is printed instead of failing to start the TUI.
    if args.len() < 2 {
        use std::io::IsTerminal;
        if !std::io::stdout().is_terminal() || !std::io::stdin().is_terminal() {
            eprintln!("{}", overview_text());
            std::process::exit(1);
        }
    }
    // If no arguments or the first argument is "tui", start the terminal UI.
    // The TUI provides a graphical interface for all subcommands.
    if args.len() < 2 || args[1] == "tui" {
        if let Err(e) = tui::run() {
            eprintln!("TUI error: {e}");
            std::process::exit(1);
        }
        return;
    }

    // Help flags redirect to the usage text.
    if args[1] == "-V" || args[1] == "--version" {
        println!("optik {}", env!("CARGO_PKG_VERSION"));
        return;
    }
    if args[1] == "-h" || args[1] == "--help" || args[1] == "help" {
        // `optik help map` is `optik map --help`.
        match args
            .get(2)
            .and_then(|c| command_help_text(if c == "build" { "build-spss" } else { c }))
        {
            Some(t) => eprintln!("{t}"),
            None => usage(),
        }
        return;
    }

    // `build` is an alias of build-spss.
    let sub: &String = if args[1] == "build" {
        static SPSS: std::sync::OnceLock<String> = std::sync::OnceLock::new();
        SPSS.get_or_init(|| "build-spss".to_string())
    } else {
        &args[1]
    };
    let rest = &args[2..];
    // `optik <command> --help` prints that command's flags and examples.
    if rest.iter().any(|a| a == "-h" || a == "--help") {
        match command_help_text(sub) {
            Some(t) => eprintln!("{t}"),
            None => usage(),
        }
        return;
    }
    // A command given no arguments at all prints its help, since every command needs an input.
    if rest.is_empty()
        && let Some(t) = command_help_text(sub)
    {
        eprintln!("{t}");
        std::process::exit(1);
    }

    // Every argument must be a flag this subcommand knows (a typo such as "-K 25" or "--thread 4" would otherwise be ignored silently and the run would use a default), and every value flag needs its value.
    {
        let Some((vals, bools)) = subcommand_flags(sub) else {
            eprintln!("error: unknown command '{sub}'");
            if let Some(c) = closest_command(sub) {
                eprintln!("did you mean: optik {c}");
            }
            eprintln!("\n{}", overview_text());
            std::process::exit(1);
        };
        // A value flag given twice is an error, since the second would be ignored silently.
        // The exceptions take several values, namely spss-update --add/--remove and -i of the k-mer counting commands (several input files, such as paired-end R1 and R2).
        let multi = |f: &str| {
            f == "--add"
                || f == "--remove"
                || (f == "-i" && matches!(sub.as_str(), "dbruijn" | "build-spss" | "build-fm"))
                || (f == "-q" && sub == "map")
        };
        let mut seen: Vec<&str> = Vec::new();
        let mut i = 0;
        while i < rest.len() {
            let a = rest[i].as_str();
            if vals.contains(&a) {
                if i + 1 >= rest.len() {
                    arg_error(sub, &format!("{a} needs a value"));
                }
                if seen.contains(&a) && !multi(a) {
                    arg_error(sub, &format!("{a} given more than once"));
                }
                seen.push(a);
                i += 2;
            } else if bools.contains(&a) {
                i += 1;
            } else {
                arg_error(sub, &format!("unknown argument '{a}' for {sub}"));
            }
        }
    }

    // Simple flag parser that looks for -flag value pairs in the argument list
    let get = |flag: &str| -> Option<&str> {
        rest.windows(2)
            .find(|w| w[0] == flag)
            .map(|w| w[1].as_str())
    };
    fn required<'a>(
        sub: &str,
        get: &dyn Fn(&str) -> Option<&'a str>,
        flag: &str,
        msg: &str,
    ) -> &'a str {
        get(flag).unwrap_or_else(|| arg_error(sub, msg))
    }
    // Strict numeric flags.
    // A flag that is present but malformed or out of its type's range is an error, never a silent fall-back to the default.
    let num = |flag: &str, default: f64| -> f64 {
        match get(flag) {
            None => default,
            Some(v) => match v.parse::<f64>() {
                Ok(x) if x.is_finite() && x == x.trunc() => x,
                _ => {
                    eprintln!("error: {flag} expects an integer, got '{v}'");
                    std::process::exit(1);
                }
            },
        }
    };
    let in_range = |flag: &str, x: f64, lo: f64, hi: f64| -> f64 {
        if x < lo || x > hi {
            eprintln!("error: {flag} {x} out of range [{lo}, {hi}]");
            std::process::exit(1);
        }
        x
    };
    // map seeds default to k=21, since at 1% error most 31-mers of a 150 bp read contain an error (P(error-free 31-mer) = 0.73 vs 0.81), and reads with 5-6 spread edits keep no 31-mer seed at all
    let k_raw = num("-k", if sub == "map" { 21.0 } else { 31.0 });
    if k_raw < 1.0 || k_raw > MAX_K as f64 {
        eprintln!(
            "error: k={k_raw} is out of range [1, {MAX_K}] (k-mers are 2-bit packed in 64 bits)"
        );
        std::process::exit(1);
    }
    let k = k_raw as usize;
    // -w 0 auto-tunes, and the minimizer window ring holds at most 256 k-mers
    let w = in_range("-w", num("-w", 10.0), 0.0, 256.0) as usize;
    // map -e is an explicit fixed cap (>= 0) or, by default, auto at 15% of each read's length.
    let e: i32 = match get("-e") {
        None | Some("auto") => -1,
        Some(_) => in_range("-e", num("-e", 0.0), 0.0, 1_000_000.0) as i32,
    };
    let c = in_range("-c", num("-c", 2.0), 0.0, u32::MAX as f64) as u32;
    // build-fm and build-spss --compact take the same -c, whose default of 1 keeps every distinct k-mer for exact membership.
    let fm_min_count = in_range("-c", num("-c", 1.0), 0.0, u32::MAX as f64) as u64;
    let threads: usize = match get("-t") {
        None => default_threads(),
        Some(_) => match in_range("-t", num("-t", 0.0), 0.0, 4096.0) as usize {
            0 => default_threads(),
            t => t,
        },
    };
    // OpenMP regions in the C++ core follow -t too, though they default to all cores.
    // omp_set_num_threads covers regions entered from this thread, and the environment set by omp_env_for_budget covers every region's default and the idle policy.
    omp_env_for_budget(threads.max(1));
    unsafe { ok_set_num_threads(threads.max(1) as u32) };

    let m_val = in_range("-m", num("-m", 0.0), 0.0, 31.0) as usize;

    // k-mer counting memory limit and spill directory (dbruijn, build-fm, build-spss --compact).
    // Past the limit the counter's partition buffers go to disk, by default in the output's directory rather than /tmp, which is often RAM-backed (tmpfs).
    if matches!(sub.as_str(), "dbruijn" | "build-fm" | "build-spss") {
        let max_mem = match get("--max-mem") {
            None => 0u64,
            Some(v) => parse_size(v).unwrap_or_else(|| {
                eprintln!("error: --max-mem expects a size such as 8G, 512M or 4.5G, got '{v}'");
                std::process::exit(1);
            }),
        };
        const MIN_MAX_MEM: u64 = 64 << 20;
        if get("--max-mem").is_some() && max_mem < MIN_MAX_MEM {
            eprintln!("error: --max-mem must be at least 64M");
            std::process::exit(1);
        }
        let dir = match get("--tmp-dir") {
            Some(d) => d.to_string(),
            None => {
                let out = get("-o").unwrap_or(if sub == "dbruijn" { "graph.gfa" } else { "." });
                match std::path::Path::new(out).parent() {
                    Some(p) if !p.as_os_str().is_empty() => p.to_string_lossy().into_owned(),
                    _ => ".".to_string(),
                }
            }
        };
        if !std::path::Path::new(&dir).is_dir() {
            eprintln!("error: spill directory '{dir}' does not exist (--tmp-dir)");
            std::process::exit(1);
        }
        let cdir = CString::new(dir).unwrap_or_else(|_| {
            eprintln!("error: --tmp-dir contains a NUL byte");
            std::process::exit(1);
        });
        // --graph chooses how the unitigs are built (by default partitioned from k = 19 up and
        // the table below that), with identical output either way.
        let graph_mode: c_int = match get("--graph") {
            None | Some("auto") => 0,
            Some("part") | Some("partitioned") => 1,
            Some("table") => 2,
            Some(v) => {
                eprintln!("error: --graph expects auto, part or table, got '{v}'");
                std::process::exit(1);
            }
        };
        if graph_mode == 1 && get("-s").is_some() {
            eprintln!("error: --graph part keeps no node table, so it cannot write the -s archive");
            std::process::exit(1);
        }
        if graph_mode == 1 && k > 31 {
            eprintln!("error: --graph part needs k <= 31");
            std::process::exit(1);
        }
        unsafe { ok_set_count_limits(max_mem, cdir.as_ptr(), graph_mode) };
    }

    // An output path must not name an input, since creating the output truncates it before it is read.
    // spss-update -i X -o X is the supported in-place update, because the index is fully loaded before it is rewritten.
    {
        let outs: Vec<&str> = ["-o", "--ids", "-s"]
            .iter()
            .filter_map(|f| get(f))
            .collect();
        let mut ins: Vec<(&str, &str)> = Vec::new();
        for flag in ["-i", "-q", "-r", "--add", "--remove"] {
            for w in rest.windows(2).filter(|w| w[0] == flag) {
                if sub == "spss-update" && flag == "-i" {
                    continue;
                }
                ins.push((flag, w[1].as_str()));
            }
        }
        // Two outputs on one file would let the second overwrite the first.
        for (a, oa) in outs.iter().enumerate() {
            for ob in &outs[a + 1..] {
                if *oa != "/dev/null" && (oa == ob || same_file(oa, ob)) {
                    eprintln!("error: two outputs name the same file {oa}");
                    std::process::exit(1);
                }
            }
        }
        for o in &outs {
            for (flag, i) in &ins {
                if *i != "-" && same_file(o, i) {
                    eprintln!(
                        "error: output {o} is the {flag} input {i}; refusing to overwrite it"
                    );
                    std::process::exit(1);
                }
            }
        }
    }

    // Every -i of a counting command, joined by '\n' for the C++ reader to split, so a path containing a newline cannot be passed this way
    let inputs: String = {
        let v: Vec<&str> = rest
            .windows(2)
            .filter(|w| w[0] == "-i")
            .map(|w| w[1].as_str())
            .collect();
        if v.iter().any(|p| p.contains('\n')) {
            eprintln!("error: input paths may not contain a newline");
            std::process::exit(1);
        }
        if v.iter().filter(|p| **p == "-").count() > 1 {
            eprintln!("error: stdin (-) can be given as -i only once");
            std::process::exit(1);
        }
        v.join("\n")
    };
    let need_i = |msg: &str| -> &str {
        if inputs.is_empty() {
            arg_error(sub, msg)
        }
        inputs.as_str()
    };

    // Dispatch to the appropriate subcommand module.
    let result = match sub.as_str() {
"map"     => {
            // -q R1 -q R2 means paired-end, as minimap2 and BWA take two read files
            let qs: Vec<&str> = rest.windows(2).filter(|x| x[0] == "-q").map(|x| x[1].as_str()).collect();
            if qs.len() > 2 {
                eprintln!("error: map takes one -q (single-end) or two (paired-end R1 R2)");
                std::process::exit(1);
            }
            if qs.len() == 2 && qs.iter().all(|q| *q == "-") {
                eprintln!("error: stdin (-) can be only one of the two read files");
                std::process::exit(1);
            }
            map::run(required(sub, &get, "-r", "-r required"),
                               required(sub, &get, "-q", "-q required"), qs.get(1).copied(), k,
                               get("-k").is_none(), w, e,
                               get("-o"), threads,
                               if rest.iter().any(|x| x == "--cs=long") { map::CsMode::Long }
                               else if rest.iter().any(|x| x == "--cs" || x == "--cs=short") { map::CsMode::Short }
                               else { map::CsMode::Off },
                               get("-R"), &args.join(" "))
        }
        "dbruijn" => dbruijn::run(need_i("-i required"), k, c,
                                  get("-o").unwrap_or("graph.gfa"),
                                  get("-s"), threads),
        "spss-update" => {
            let all = |flag: &str| -> Vec<&str> {
                rest.windows(2).filter(|w| w[0] == flag).map(|w| w[1].as_str()).collect()
            };
            build::update_spss(required(sub, &get, "-i", "-i <index> required"),
                               &all("--add"), &all("--remove"), threads,
                               required(sub, &get, "-o", "-o <index> required"))
        }
        // -c filters counted k-mers, which only --compact does.
        // Without it the records are indexed as given and -c would be ignored silently.
        "build-spss" if get("-c").is_some() && !rest.iter().any(|a| a == "--compact") => Err(
            "-c (minimum abundance) applies only with --compact; without it every k-mer of the \
             input is indexed as given".into()),
        // The counter's memory limit likewise applies only to --compact, which counts k-mers.
        "build-spss" if (get("--max-mem").is_some() || get("--tmp-dir").is_some() || get("--graph").is_some())
                        && !rest.iter().any(|a| a == "--compact") => Err(
            "--max-mem / --tmp-dir / --graph apply only with --compact (the k-mer counter); without it the \
             records are indexed as given".into()),
        "build-spss" => build::build_spss(need_i("-i <ref.fa> required"),
                                          k, m_val, rest.iter().any(|a| a == "--compact"),
                                          fm_min_count, threads,
                                          required(sub, &get, "-o", "-o <index> required")),
        "build-fm"  => build::build_fm(need_i("-i <ref.fa> required"),
                                       k, fm_min_count, threads as u32,
                                       required(sub, &get, "-o", "-o <index.bin> required")),
        "query"     => query::query(required(sub, &get, "-i", "-i <index.bin> required"),
                                     required(sub, &get, "-q", "-q required"),
                                     get("-k").map(|_| k), threads,
                                     get("--ids"), !rest.iter().any(|a| a == "--no-verify")),
        _ => { usage(); std::process::exit(1); }
    };

    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod help_tests {
    const COMMANDS: [&str; 6] = [
        "map",
        "dbruijn",
        "build-spss",
        "build-fm",
        "query",
        "spss-update",
    ];

    #[test]
    fn each_command_help_lists_its_flags() {
        for sub in COMMANDS {
            let help = super::command_help_text(sub).unwrap_or_else(|| panic!("no help for {sub}"));
            let (vals, bools) = super::subcommand_flags(sub).unwrap();
            for f in vals.iter().chain(bools) {
                let name = f.split('=').next().unwrap();
                assert!(
                    help.lines().any(|l| l.trim_start().starts_with(name)),
                    "{sub} --help does not document {f}"
                );
            }
            assert!(
                help.contains("EXAMPLES") && help.contains(&format!("optik {sub} ")),
                "{sub}"
            );
        }
    }

    #[test]
    fn overview_names_every_command_and_full_help_has_every_block() {
        let overview = super::overview_text();
        let full = super::usage_text();
        for sub in COMMANDS {
            assert!(
                overview.contains(&format!("  {sub} ")),
                "overview lacks {sub}"
            );
            assert!(
                full.contains(&format!("FLAGS ({sub})")),
                "--help lacks {sub}"
            );
        }
    }

    #[test]
    fn typos_suggest_the_command() {
        assert_eq!(super::closest_command("qeury"), Some("query"));
        assert_eq!(super::closest_command("dbrujin"), Some("dbruijn"));
        assert_eq!(super::closest_command("buildfm"), Some("build-fm"));
        assert_eq!(super::closest_command("zzzzzzzz"), None);
    }
}

#[cfg(test)]
mod correctness_tests {
    use std::collections::HashSet;

    const K: usize = 31;

    /// Reverse complement with the C++ encoding (A=0, C=1, G=2, T=3, so the complement is XOR 3).
    fn rc(pk: u64, k: usize) -> u64 {
        let mut r = 0u64;
        for i in 0..k {
            let b = (pk >> (2 * i)) & 3;
            r |= (b ^ 3) << (2 * (k - 1 - i));
        }
        r
    }

    /// Packs a DNA string with the C++ encoding, without canonicalising it.
    fn pack(seq: &[u8]) -> u64 {
        let mut pk = 0u64;
        for &b in seq {
            pk = (pk << 2)
                | match b {
                    b'A' => 0u64,
                    b'C' => 1,
                    b'G' => 2,
                    _ => 3, // T=3 to match C++ base_encode
                };
        }
        pk
    }

    /// Spells a packed k-mer as bases.
    fn decode(pk: u64, k: usize) -> Vec<u8> {
        (0..k)
            .rev()
            .map(|i| b"ACGT"[((pk >> (2 * i)) & 3) as usize])
            .collect()
    }
    /// Membership of one k-mer through the scan path `optik query` uses.
    fn spss_has(idx: &crate::SpssIdx, pk: u64) -> bool {
        idx.scan_chunk(&decode(pk, K)).1 == 1
    }
    fn fm_has(idx: &crate::FmIdx, pk: u64) -> bool {
        idx.scan_chunk(&decode(pk, K), K).1 == 1
    }

    /// Writes a synthetic tandem-repeat FASTA (220 records, each a random 400 bp unit repeated 5 times) and returns its path and canonical k-mer set.
    fn synthetic_tandem_fasta(path: &std::path::Path, seed: u64) -> (String, HashSet<u64>) {
        let mut fa = String::new();
        let mut staged: HashSet<u64> = HashSet::new();
        let mut rng = seed;
        for rec in 0..220usize {
            let mut seed2 = Vec::with_capacity(400);
            for _ in 0..400 {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                seed2.push(match rng & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                });
            }
            let mut seq = Vec::with_capacity(2000);
            for _ in 0..5 {
                seq.extend_from_slice(&seed2);
            }
            for w in 0..=seq.len() - K {
                let pk = pack(&seq[w..w + K]);
                staged.insert(pk.min(rc(pk, K)));
            }
            fa.push_str(&format!(">r{rec}\n{}\n", String::from_utf8(seq).unwrap()));
        }
        std::fs::write(path, &fa).unwrap();
        (path.to_str().unwrap().to_string(), staged)
    }

    /// --max-mem sizes with K/M/G/T suffixes in powers of 1024, fractions, and the
    /// rejection of malformed or non-positive values.
    #[test]
    fn parse_size_units() {
        use crate::parse_size;
        assert_eq!(parse_size("1024"), Some(1024));
        assert_eq!(parse_size("64M"), Some(64 << 20));
        assert_eq!(parse_size("8g"), Some(8 << 30));
        assert_eq!(parse_size("1.5G"), Some(3 << 29));
        assert_eq!(parse_size("2T"), Some(2 << 40));
        assert_eq!(parse_size("8GB"), Some(8 << 30));
        assert_eq!(parse_size("8GiB"), Some(8 << 30));
        assert_eq!(parse_size("500mb"), Some(500 << 20));
        for bad in [
            "", "G", "GB", "-1G", "0", "abc", "1.2.3M", "NaNG", "infG", "8B",
        ] {
            assert_eq!(parse_size(bad), None, "{bad}");
        }
    }

    /// SPSS index (build-spss --compact), with exact membership on every positional k-mer of a tandem-repeat reference and no false positive on random k-mers, before and after a save/load round trip, and a distinct count equal to the k-mer set.
    #[test]
    fn spss_index_is_exact() {
        use crate::SpssIdx;
        let dir = std::env::temp_dir().join(format!("optik_spss_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let (fa_path, staged) = synthetic_tandem_fasta(&dir.join("sp.fa"), 0x9E3779B97F4A7C15u64);
        assert!(
            staged.len() >= 65_536,
            "test setup: too few distinct k-mers"
        );

        let check = |idx: &SpssIdx, staged: &HashSet<u64>| {
            assert_eq!(idx.k(), K);
            assert_eq!(idx.distinct() as usize, staged.len());
            for &ks in staged {
                assert!(spss_has(idx, ks), "covered k-mer {:x} not found", ks);
            }
            let mut neg: Vec<u64> = Vec::with_capacity(2000);
            let mut r = 0x123456789ABCDEF0u64;
            while neg.len() < 2000 {
                r ^= r << 13;
                r ^= r >> 7;
                r ^= r << 17;
                let pk = (r & ((1u64 << (2 * K)) - 1)) | 1;
                let c = pk.min(rc(pk, K));
                if !staged.contains(&c) {
                    neg.push(c);
                }
            }
            for &ks in &neg {
                assert!(!spss_has(idx, ks), "false positive {:x}", ks);
            }
        };

        let idx = SpssIdx::build_from_fasta(&fa_path, K, K / 2, true, 1, 0).expect("SPSS build");
        check(&idx, &staged);
        let bin = format!("{}/sp.spss", dir.to_str().unwrap());
        idx.save(&bin).expect("save");
        let idx2 = SpssIdx::load(&bin, true).expect("load");
        check(&idx2, &staged);

        std::fs::remove_dir_all(&dir).ok();
    }

    /// FM index with exact membership.
    /// Every k-mer of a synthetic reference must be found and random non-reference k-mers rejected, before and after a save/load round trip, and the saved file must be exactly bytes() long.
    #[test]
    fn fm_index_is_exact() {
        use crate::FmIdx;
        let dir = std::env::temp_dir().join(format!("optik_fm_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let (fa_path, staged) = synthetic_tandem_fasta(&dir.join("fm.fa"), 0xABCDEF0123456789u64);
        assert!(
            staged.len() >= 65_536,
            "test setup: too few distinct k-mers"
        );

        let check = |idx: &FmIdx, staged: &HashSet<u64>| {
            assert_eq!(idx.k(), K);
            assert_eq!(
                idx.size(),
                staged.len(),
                "indexed {} distinct k-mers, expected {}",
                idx.size(),
                staged.len()
            );
            for &ks in staged {
                assert!(fm_has(idx, ks), "covered k-mer {:x} not found", ks);
            }
            let mut neg: Vec<u64> = Vec::with_capacity(2000);
            let mut r = 0xFEDCBA9876543210u64;
            while neg.len() < 2000 {
                r ^= r << 13;
                r ^= r >> 7;
                r ^= r << 17;
                let pk = (r & ((1u64 << (2 * K)) - 1)) | 1;
                let c = pk.min(rc(pk, K));
                if !staged.contains(&c) {
                    neg.push(c);
                }
            }
            for &ks in &neg {
                assert!(!fm_has(idx, ks), "false positive {:x}", ks);
            }
        };

        let idx = FmIdx::build_from_fasta(&fa_path, K, 1, 0).expect("FM build");
        let est = idx.bytes() as f64;
        assert!(
            est < idx.size() as f64 * 8.0,
            "FM index {est:.0} B should stay below 8 b/k ({} kmers)",
            idx.size()
        );
        check(&idx, &staged);

        // The round trip's load path re-reads the masked-BWT rank arrays.
        let bin = format!("{}/fm.bin", dir.to_str().unwrap());
        idx.save(&bin).expect("save");
        assert_eq!(
            std::fs::metadata(&bin).map(|m| m.len() as usize).unwrap(),
            idx.bytes(),
            "serialized file size must equal bytes()"
        );
        let idx2 = FmIdx::load(&bin).expect("load");
        check(&idx2, &staged);

        std::fs::remove_dir_all(&dir).ok();
    }
}
