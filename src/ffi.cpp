/**
 * @file
 * @brief C ABI shared by every subcommand, covering FASTA/FASTQ iteration, the counter's memory limits and the OpenMP thread cap.
 * The index-specific bridges are in ffi_spss.cpp, ffi_fm.cpp, ffi_dbg.cpp and ffi_map.cpp.
 * On the Rust side, tools/src/main.rs declares all of them.
 *
 * Failable calls return kStatusOk or kStatusFailed.
 * After a failure, the exception text (at most errbuf_len - 1 bytes) sits in the zeroed buffer the caller supplied.
 */
#include <cstddef>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>

#include "../include/fasta_reader.hpp"
#include "../include/spill.hpp"
#include "ffi_error.hpp"

namespace {

/// Thrown when a record callback returns nonzero, which unwinds the reader like any I/O error.
/// The caller already knows why it stopped, so the message only marks the cause.
struct ReadStopped : std::runtime_error {
    ReadStopped() : std::runtime_error("record callback stopped the read") {}
};

template <typename Body>
int run_reporting_errors(char* errbuf, uint64_t errbuf_len, Body&& body) {
    try {
        body();
        return tpoptoa::ffi::kStatusOk;
    } catch (const std::exception& err) {
        if (errbuf && errbuf_len)
            std::strncpy(errbuf, err.what(), errbuf_len - 1);

        return tpoptoa::ffi::kStatusFailed;
    }
}

}

extern "C" {

/// (cookie, name, name_len, seq, seq_len), returning 0 to go on or nonzero to stop the read.
typedef int (*SeqCb)(void*, const char*, uint64_t, const char*, uint64_t);

/// SeqCb plus (qual, qual_len), where the quality string is empty for FASTA.
typedef int (*ReadCb)(void*, const char*, uint64_t, const char*, uint64_t, const char*, uint64_t);

/// Skips records with an empty sequence.
int ok_iterate_fasta(const char* path, SeqCb callback, void* cookie,
                     char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& rec) {
            if (callback(cookie,
                         rec.name.data(), (uint64_t)rec.name.size(),
                         rec.seq.data(),  (uint64_t)rec.seq.size()))
                throw ReadStopped();
        });
    });
}

/// Like ok_iterate_fasta() but keeps empty records.
int ok_iterate_fasta_all(const char* path, SeqCb callback, void* cookie,
                         char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& rec) {
            if (callback(cookie,
                         rec.name.data(), (uint64_t)rec.name.size(),
                         rec.seq.data(),  (uint64_t)rec.seq.size()))
                throw ReadStopped();
        }, /*keep_empty=*/true);
    });
}

/// Keeps empty records and passes qualities along, which the mapper copies into SAM QUAL.
int ok_iterate_reads(const char* path, ReadCb callback, void* cookie,
                     char* errbuf, uint64_t errbuf_len) {
    return run_reporting_errors(errbuf, errbuf_len, [&] {
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& rec) {
            if (callback(cookie,
                         rec.name.data(), (uint64_t)rec.name.size(),
                         rec.seq.data(),  (uint64_t)rec.seq.size(),
                         rec.qual.data(), (uint64_t)rec.qual.size()))
                throw ReadStopped();
        }, /*keep_empty=*/true);
    });
}

/// Settings for every k-mer-counting build (spill.hpp). @p max_mem = 0 picks the limit automatically.
void ok_set_count_limits(uint64_t max_mem, const char* tmp_dir, int graph_mode) {
    tpoptoa::count_limits().graph_mode = graph_mode;
    tpoptoa::count_limits().max_mem = max_mem;
    tpoptoa::count_limits().tmp_dir = tmp_dir ? tmp_dir : "";
}

/// Caps OpenMP regions started from this thread, while 0 leaves the setting alone.
/// Because that only covers the calling thread, the CLI also re-execs itself with OMP_NUM_THREADS set.
/// See omp_env_for_budget() in tools/src/main.rs, which makes regions on other threads obey -t too.
void ok_set_num_threads(unsigned num_threads) {
#ifdef _OPENMP
    if (num_threads > 0) omp_set_num_threads(static_cast<int>(num_threads));
#else
    (void)num_threads;
#endif
}

}
