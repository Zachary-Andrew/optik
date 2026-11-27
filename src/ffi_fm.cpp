/**
 * @file
 * @brief The `okfm_*` C ABI behind build-fm and query on FM indexes.
 *
 * Building counts k-mers into a de Bruijn graph (kmer_count.cpp) and walks its maximal unitigs into one separator-joined text.
 * That text is then indexed with the k-order BWT in fm_index.hpp, and lookups are exact.
 */
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "../include/dbg_graph.hpp"
#include "../include/dbg_graph_ef.hpp"
#include "../include/dbg_part.hpp"
#include "../include/fm_index.hpp"
#include "kmer_count.hpp"
#include "../include/kmer_extractor.hpp"
#include "../include/phase_timer.hpp"
#include "ffi_error.hpp"

using tpoptoa::ffi::set_error;
using tpoptoa::ffi::kStatusFailed;
using tpoptoa::ffi::kStatusOk;

/// Text capacity is n_kmers * (1 + 1/kTextSlackDiv), with the extra reserved for separators.
constexpr std::size_t kTextSlackDiv = 4;

extern "C" {

struct OkFmIndex { tpoptoa::FmIndex* idx; bool batch; };

/// @p min_count = 1 keeps every k-mer.
/// Errors come back as nullptr plus a message in @p errbuf.
OkFmIndex* okfm_build_from_fasta(const char* fasta_path, uint64_t kmer_len_u64,
                                 uint64_t min_count, uint32_t nthreads,
                                 char* errbuf, uint64_t errbuf_len) {
    try {
        if (kmer_len_u64 < 1 || kmer_len_u64 > tpoptoa::kMaxK) {
            if (errbuf && errbuf_len)
                std::snprintf(errbuf, errbuf_len,
                              "k must be in [1, %zu] (got %llu)",
                              tpoptoa::kMaxK, (unsigned long long)kmer_len_u64);

            return nullptr;
        }

        tpoptoa::PhaseTimer timer;
        const unsigned kmer_len = static_cast<unsigned>(kmer_len_u64);
        // The FM index only needs the text, so the graph is scoped to the walk and freed before the BWT sort.
        // Otherwise both would be resident at the memory peak.
        std::vector<uint8_t> text;
        std::size_t n_kmers = 0, n_units = 0;

        {
            tpoptoa::DbgGraph graph(kmer_len);
            std::unique_ptr<tpoptoa::DbgGraphEF> ef_graph;       // set if the plain table won't fit
            std::unique_ptr<tpoptoa::PartUnitigs> part_unitigs;

            // The counter fills exactly one of the three graph forms, chosen by k and by whether the plain table fits the memory limit.
            tpoptoa::detail::count_into_graph(fasta_path, static_cast<std::size_t>(kmer_len),
                                              static_cast<uint32_t>(min_count), nthreads, graph, &ef_graph,
                                              tpoptoa::part_mode_wanted(kmer_len) ? &part_unitigs : nullptr);
            n_kmers = part_unitigs ? part_unitigs->core_size() : ef_graph ? ef_graph->core_size() : graph.core_size();

            if (n_kmers == 0) {
                set_error(errbuf, errbuf_len, "no k-mers extracted from input");
                return nullptr;
            }

            timer.lap("fm: count k-mers + graph");
            text.reserve(n_kmers + n_kmers / kTextSlackDiv);
            auto add = [&](const uint8_t* seq, uint32_t len, uint64_t) {
                tpoptoa::FmIndex::append_unitig(text, reinterpret_cast<const char*>(seq), len);
                ++n_units;
            };

            // Whichever graph was filled is walked into unitigs, and each one goes into the text followed by its separator.
            if (part_unitigs) part_unitigs->unitigs(add, nthreads);
            else if (ef_graph) ef_graph->unitigs(add, nthreads);
            else graph.unitigs(add, nthreads);
        }
        timer.lap("fm: unitig walk");
        auto* handle = new OkFmIndex();

        handle->idx = new tpoptoa::FmIndex(kmer_len, std::move(text), n_units, n_kmers, nthreads);
        // Batched lookups only help once the rank blocks outgrow the last-level cache, which batch_pays() decides from the index size.
        handle->batch = handle->idx->batch_pays();
        timer.lap("fm: FM index (k-order BWT, rank)");
        return handle;
    } catch (const std::exception& err) {
        set_error(errbuf, errbuf_len, err.what());
        return nullptr;
    } catch (...) {
        set_error(errbuf, errbuf_len, "build failed (allocation error)");
        return nullptr;
    }
}

/// Loads a saved FM index, or returns null with the reason in @p errbuf when the file is not a valid one.
OkFmIndex* okfm_load(const char* path, char* errbuf, uint64_t errbuf_len) {
    try {
        auto* idx = tpoptoa::FmIndex::load(path);

        if (!idx) {
            set_error(errbuf, errbuf_len,
                      "not a valid optik FM index (corrupt, truncated or written by an older optik); "
                      "rebuild it with build-fm");
            return nullptr;
        }

        return new OkFmIndex{idx, idx->batch_pays()};
    } catch (...) {
        set_error(errbuf, errbuf_len, "load failed (allocation error)");
        return nullptr;
    }
}

/// Writes the index to @p path and returns kStatusOk, or kStatusFailed with the reason in @p errbuf.
int okfm_save(const OkFmIndex* handle, const char* path,
              char* errbuf, uint64_t errbuf_len) {
    try {
        if (!handle->idx->save(path)) {
            set_error(errbuf, errbuf_len, "save failed (I/O error)");
            return kStatusFailed;
        }

        return kStatusOk;
    } catch (...) {
        set_error(errbuf, errbuf_len, "save failed (I/O error)");
        return kStatusFailed;
    }
}

/// Frees the handle together with its index, and accepts null.
void okfm_free(OkFmIndex* handle) { if (handle) { delete handle->idx; delete handle; } }

// Read-only figures for the summaries that build-fm and query print.
unsigned    okfm_k(const OkFmIndex* handle)           { return handle->idx->kmer_len(); }
uint64_t    okfm_size(const OkFmIndex* handle)        { return handle->idx->size(); }
uint64_t    okfm_num_unitigs(const OkFmIndex* handle) { return handle->idx->num_unitigs(); }
uint64_t    okfm_super_len(const OkFmIndex* handle)   { return handle->idx->super_len(); }
uint64_t    okfm_bwt_len(const OkFmIndex* handle)     { return handle->idx->bwt_len(); }
uint64_t    okfm_bytes(const OkFmIndex* handle)       { return handle->idx->bytes(); }
double      okfm_bits_per_entry(const OkFmIndex* handle) { return handle->idx->bits_per_entry(); }

/// Handles one record per FFI call.
/// @p num_searches follows SSHash's definition, counting lookups that aren't streaming extensions.
/// Because this scan never extends, it always equals @p num_valid.
void okfm_scan_seq(const OkFmIndex* handle, const char* seq, uint64_t slen, uint64_t kmer_len,
                   uint64_t* num_valid, uint64_t* num_positive, uint64_t* num_searches) {
    // Batching lets find_batch() overlap cache misses once the index outgrows the LLC.
    // The buffers are thread_local, so consecutive records reuse them.
    constexpr std::size_t kScanChunk = 4096;
    thread_local std::vector<tpoptoa::KmerWord> kmers;
    thread_local std::vector<uint8_t> hits;

    kmers.clear();
    uint64_t valid_count = 0, positive_count = 0;
    auto flush = [&] {
        hits.resize(kmers.size());

        // batch_pays() is false while the rank blocks fit in the LLC, where plain find() is faster than the batched lanes.
        if (handle->batch) {
            handle->idx->find_batch(kmers.data(), kmers.size(), hits.data());
        } else {
            for (std::size_t i = 0; i < kmers.size(); ++i) hits[i] = handle->idx->find(kmers[i]);
        }

        for (uint8_t hit : hits) positive_count += hit;

        kmers.clear();
    };

    tpoptoa::extract_kmers(seq, static_cast<std::size_t>(slen),
                           static_cast<std::size_t>(kmer_len),
                           [&](tpoptoa::KmerWord kmer) {
        ++valid_count;
        kmers.push_back(kmer);

        if (kmers.size() == kScanChunk) flush();
    });
    if (!kmers.empty()) flush();

    *num_valid = valid_count;
    *num_positive = positive_count;
    // Every valid k-mer counts as a search because FmIndex::find() never reuses the previous k-mer's interval.
    *num_searches = valid_count;
}

}  // extern "C"