/*
 * ffi.cpp - C ABI for the helpers shared by every subcommand, namely FASTA/FASTQ iteration, the k-mer counter's memory limits and the OpenMP thread cap.
 *
 * Each index has its own bridge in ffi_spss.cpp, ffi_fm.cpp, ffi_dbg.cpp and ffi_map.cpp.
 * The Rust declarations of these functions are in tools/src/main.rs.
 *
 * By the error convention, a function that can fail returns 0 on success and -1 on failure, after copying at most errbuf_len - 1 bytes of the exception message into the zeroed buffer the caller passes.
 */
#include <cstddef>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstring>
#include <exception>
#include <string>

#include "../include/fasta_reader.hpp"
#include "../include/spill.hpp"

namespace {

/// run_reporting_errors() - run @body, returning 0, or -1 after an exception whose message is copied to @errbuf.
template <typename Body>
int run_reporting_errors(char* errbuf, uint64_t errbuf_len, Body&& body) {
    try {
        body();
        return 0;
    } catch (const std::exception& e) {
        if (errbuf && errbuf_len)
            std::strncpy(errbuf, e.what(), errbuf_len - 1);
        return -1;
    }
}

}

extern "C" {

/// Record callback taking the cookie, name, name length, sequence and sequence length.
typedef void (*SeqCb)(void*, const char*, uint64_t, const char*, uint64_t);

/// Record callback with the FASTQ quality string appended, empty for FASTA.
typedef void (*ReadCb)(void*, const char*, uint64_t, const char*, uint64_t, const char*, uint64_t);

/// ok_iterate_fasta() - call @cb for every record of @path that has a non-empty sequence.
int ok_iterate_fasta(const char* path, SeqCb cb, void* cookie,
                     char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& r) {
            cb(cookie,
               r.name.data(), (uint64_t)r.name.size(),
               r.seq.data(),  (uint64_t)r.seq.size());
        });
    });
}

/// ok_iterate_fasta_all() - as ok_iterate_fasta(), but records with an empty sequence are reported too.
int ok_iterate_fasta_all(const char* path, SeqCb cb, void* cookie,
                         char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& r) {
            cb(cookie,
               r.name.data(), (uint64_t)r.name.size(),
               r.seq.data(),  (uint64_t)r.seq.size());
        }, /*keep_empty=*/true);
    });
}

/// ok_iterate_reads() - every record, empty ones included, with its FASTQ quality string.
/// The mapper writes the qualities to SAM QUAL.
int ok_iterate_reads(const char* path, ReadCb cb, void* cookie,
                     char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& r) {
            cb(cookie,
               r.name.data(), (uint64_t)r.name.size(),
               r.seq.data(),  (uint64_t)r.seq.size(),
               r.qual.data(), (uint64_t)r.qual.size());
        }, /*keep_empty=*/true);
    });
}

/// ok_set_count_limits() - the counter's memory limit (@max_mem bytes, 0 = automatic), spill directory and graph mode (spill.hpp), used by every build that counts k-mers (dbruijn, build-fm, build-spss --compact).
void ok_set_count_limits(uint64_t max_mem, const char* tmp_dir, int graph_mode) {
    tpoptoa::count_limits().graph_mode = graph_mode;
    tpoptoa::count_limits().max_mem = max_mem;
    tpoptoa::count_limits().tmp_dir = tmp_dir ? tmp_dir : "";
}

/// ok_set_num_threads() - cap OpenMP regions started from the calling thread at @n, where 0 leaves the setting alone.
///
/// omp_set_num_threads() affects only the calling thread's regions.
/// The CLI also re-executes itself with OMP_NUM_THREADS set to -t, so regions entered from any thread respect it (omp_env_for_budget() in tools/src/main.rs).
void ok_set_num_threads(unsigned n) {
#ifdef _OPENMP
    if (n > 0) omp_set_num_threads(static_cast<int>(n));
#else
    (void)n;
#endif
}

}  // extern "C"
