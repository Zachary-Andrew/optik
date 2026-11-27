//! Declares the raw C ABI of the C++ core (src/ffi*.cpp) and the opaque handle types it passes around.
//! Nothing outside native.rs should call these directly.

use super::*;

/// C++ OkFmIndex, the FM k-mer index of src/ffi_fm.cpp, seen as an opaque type.
pub(crate) struct OkFmIndex {
    pub(super) _private: [u8; 0],
}

/// The compacted de Bruijn graph (OkDbg, src/ffi_dbg.cpp).
pub(crate) struct OkDbg {
    pub(super) _private: [u8; 0],
}

/// The mapper's minimizer seed index (OkSeedIdx).
pub(crate) struct OkSeedIdx {
    pub(super) _private: [u8; 0],
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
    pub(super) _private: [u8; 0],
}

/// An SPSS index with edit sets (OkDynSpss).
pub(crate) struct OkDynSpss {
    pub(super) _private: [u8; 0],
}

// The C ABI of the C++ core in src/ffi*.cpp.
// improper_ctypes is allowed because the opaque handle structs above are not #[repr(C)], and they are only ever used behind pointers.

/// Size of the error buffers passed to the index functions.
pub(super) const ERRBUF_LEN: usize = 4096;

/// Size of the error buffers passed to the FASTA/FASTQ readers.
pub(super) const READER_ERRBUF_LEN: usize = 512;

#[allow(improper_ctypes)]
unsafe extern "C" {
    // Sequence readers and global settings (src/ffi.cpp)
    pub(super) fn ok_iterate_fasta(
        path: *const c_char,
        callback: extern "C" fn(*mut c_void, *const c_char, u64, *const c_char, u64) -> c_int,
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    pub(super) fn ok_iterate_fasta_all(
        path: *const c_char,
        callback: extern "C" fn(*mut c_void, *const c_char, u64, *const c_char, u64) -> c_int,
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    pub(super) fn ok_iterate_reads(
        path: *const c_char,
        callback: extern "C" fn(
            *mut c_void,
            *const c_char,
            u64,
            *const c_char,
            u64,
            *const c_char,
            u64,
        ) -> c_int,
        cookie: *mut c_void,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    pub(super) fn ok_set_num_threads(num_threads: u32);
    pub(super) fn ok_set_count_limits(max_mem: u64, tmp_dir: *const c_char, graph_mode: c_int);

    // FM k-mer index (src/ffi_fm.cpp, include/fm_index.hpp)
    pub(super) fn okfm_build_from_fasta(
        fasta_path: *const c_char,
        kmer_len: u64,
        min_count: u64,
        nthreads: u32,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkFmIndex;
    pub(super) fn okfm_load(path: *const c_char, errbuf: *mut c_char, errbuf_len: u64) -> *mut OkFmIndex;
    pub(super) fn okfm_save(
        handle: *const OkFmIndex,
        path: *const c_char,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> c_int;
    pub(super) fn okfm_free(handle: *mut OkFmIndex);
    pub(super) fn okfm_k(handle: *const OkFmIndex) -> u32;
    pub(super) fn okfm_size(handle: *const OkFmIndex) -> u64;
    pub(super) fn okfm_num_unitigs(handle: *const OkFmIndex) -> u64;
    pub(super) fn okfm_super_len(handle: *const OkFmIndex) -> u64;
    pub(super) fn okfm_bwt_len(handle: *const OkFmIndex) -> u64;
    pub(super) fn okfm_bytes(handle: *const OkFmIndex) -> u64;
    pub(super) fn okfm_bits_per_entry(handle: *const OkFmIndex) -> f64;
    pub(super) fn okfm_scan_seq(
        handle: *const OkFmIndex,
        seq: *const c_char,
        slen: u64,
        kmer_len: u64,
        num_valid: *mut u64,
        num_positive: *mut u64,
        num_searches: *mut u64,
    );

    // Compacted de Bruijn graph (src/ffi_dbg.cpp)
    pub(super) fn okdbg_build_from_fasta(
        fasta_path: *const c_char,
        kmer_len: u64,
        min_count: u64,
        nthreads: u32,
        keep_nodes: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkDbg;
    pub(super) fn okdbg_save(handle: *const OkDbg, path: *const c_char) -> c_int;
    pub(super) fn okdbg_free(handle: *mut OkDbg);
    pub(super) fn okdbg_core_size(handle: *const OkDbg) -> u64;
    pub(super) fn okdbg_write_gfa(
        handle: *mut OkDbg,
        path: *const c_char,
        nthreads: u32,
        stats: *mut GfaStats,
        fasta: c_int,
    ) -> c_int;

    // SPSS index, static and with edit sets (src/ffi_spss.cpp)
    pub(super) fn okspss_build_from_fasta(
        path: *const c_char,
        kmer_len: u64,
        minimizer_len: u64,
        compact: c_int,
        min_count: u64,
        nthreads: u32,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkSpss;
    pub(super) fn okspss_free(handle: *mut OkSpss);
    pub(super) fn okspss_save(handle: *const OkSpss, path: *const c_char) -> c_int;
    pub(super) fn okspss_load(
        path: *const c_char,
        verify: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> *mut OkSpss;
    pub(super) fn okspss_kmers(handle: *const OkSpss) -> u64;
    pub(super) fn okspss_bytes(handle: *const OkSpss) -> u64;
    pub(super) fn okspss_k(handle: *const OkSpss) -> u32;
    pub(super) fn okspss_m(handle: *const OkSpss) -> u32;
    pub(super) fn okspss_superkmers(handle: *const OkSpss) -> u64;
    pub(super) fn okspss_distinct(handle: *const OkSpss) -> u64;
    pub(super) fn okspss_ids(handle: *const OkSpss, seq: *const c_char, len: u64, out: *mut i64) -> u64;
    pub(super) fn okspss_scan_seq(
        handle: *const OkSpss,
        seq: *const c_char,
        len: u64,
        num_valid: *mut u64,
        num_positive: *mut u64,
        num_searches: *mut u64,
        num_extensions: *mut u64,
    );
    pub(super) fn okdyn_open(path: *const c_char, errbuf: *mut c_char, errbuf_len: u64) -> *mut OkDynSpss;
    pub(super) fn okdyn_free(handle: *mut OkDynSpss);
    pub(super) fn okdyn_apply_fasta(
        handle: *mut OkDynSpss,
        path: *const c_char,
        add: c_int,
        errbuf: *mut c_char,
        errbuf_len: u64,
    ) -> i64;
    pub(super) fn okdyn_added(handle: *const OkDynSpss) -> u64;
    pub(super) fn okdyn_deleted(handle: *const OkDynSpss) -> u64;
    pub(super) fn okdyn_compact_save(handle: *mut OkDynSpss, path: *const c_char, nthreads: u32) -> c_int;
    pub(super) fn okdyn_kmers(handle: *const OkDynSpss) -> u64;
    pub(super) fn okspss_scan_many(
        handle: *const OkSpss,
        seqs: *const *const c_char,
        lens: *const u64,
        num_seqs: u64,
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
    pub(super) fn okseed_build(
        seqs: *const *const c_char,
        lens: *const u64,
        num_seqs: u64,
        kmer_len: u32,
        window: u32,
        max_occ: u32,
        nthreads: u32,
    ) -> *mut OkSeedIdx;
    pub(super) fn okseed_free(handle: *mut OkSeedIdx);
    pub(super) fn okseed_distinct(handle: *const OkSeedIdx) -> u64;
    pub(super) fn okseed_bytes(handle: *const OkSeedIdx) -> u64;
    pub(super) fn okseed_occurrences(handle: *const OkSeedIdx) -> u64;
    pub(super) fn okseed_collect(
        handle: *const OkSeedIdx,
        seq: *const c_char,
        len: u64,
        out: *mut *const SeedHit,
        n_min: *mut u32,
        n_rep: *mut u32,
    ) -> u64;
    pub(super) fn okseed_collect_chain(
        handle: *const OkSeedIdx,
        seq: *const c_char,
        len: u64,
        out: *mut *const SeedHit,
        rep_len: *mut u32,
    ) -> u64;
    pub(super) fn okseed_mid_occ(handle: *const OkSeedIdx) -> u32;
    pub(super) fn ok_wfa_affine_align(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        tef: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
        max_score: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
        out_tend: *mut i32,
    ) -> i32;
    pub(super) fn ok_biwfa_align(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
        max_score: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
    ) -> i32;
    pub(super) fn ok_placements_x4(
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
    pub(super) fn ok_wfa_affine_cigar(
        query: *const c_char,
        qlen: i32,
        text: *const c_char,
        tlen: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
        out_cig: *mut c_char,
        cig_cap: i32,
        out_score: *mut i32,
    ) -> i32;
}

/// Copies the NUL-terminated message the C++ side wrote into `buf`, or `default` if there is none.
pub(super) fn ffi_error(buf: &[u8], default: &str) -> String {
    CStr::from_bytes_until_nul(buf)
        .map(|msg| msg.to_string_lossy().into_owned())
        .unwrap_or_else(|_| default.into())
}
