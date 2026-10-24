/*
 * ffi_fm.cpp - C ABI of the FM k-mer index (build-fm, query), the okfm_* functions.
 *
 * The build path counts the canonical k-mers (minimizer-partitioned counter, kmer_count.cpp) into a de Bruijn graph with the observed (k+1)-mer edges, walk its maximal unitigs, and index the separator-delimited unitigs with an FM-index over the k-order BWT (fm_index.hpp).
 * Membership is exact.
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

extern "C" {

/// OkFmIndex - an FM index and whether its lookups go through find_batch().
struct OkFmIndex { tpoptoa::FmIndex* idx; bool batch; };

/// okfm_build_from_fasta() - build an FM index from the k-mers of a FASTA/FASTQ file.
/// @min_count filters the k-mers by abundance (1 = all).
/// Returns nullptr with a message in @errbuf on failure.
OkFmIndex* okfm_build_from_fasta(const char* fasta_path, uint64_t k,
                                 uint64_t min_count, uint32_t nthreads,
                                 char* errbuf, uint64_t errbuf_len) {
    try {
        if (k < 1 || k > tpoptoa::kMaxK) {
            if (errbuf && errbuf_len)
                std::snprintf(errbuf, errbuf_len,
                              "k must be in [1, %zu] (got %llu)",
                              tpoptoa::kMaxK, (unsigned long long)k);
            return nullptr;
        }
        tpoptoa::PhaseTimer pt;
        const unsigned K = static_cast<unsigned>(k);
        // Walk the maximal unitigs straight into the FM symbol text.
        // The The graph lives only for the walk, since the FM index needs nothing but the text, so the graph is freed before the BWT is sorted.
        std::vector<uint8_t> text;
        std::size_t n_kmers = 0, n_units = 0;
        {
            tpoptoa::DbgGraph g(K);
            std::unique_ptr<tpoptoa::DbgGraphEF> ef;       // when the node table would not fit the limit
            std::unique_ptr<tpoptoa::PartUnitigs> pu;      // partitioned compaction
            tpoptoa::detail::count_into_graph(fasta_path, static_cast<std::size_t>(K),
                                              static_cast<uint32_t>(min_count), nthreads, g, &ef,
                                              tpoptoa::part_mode_wanted(K) ? &pu : nullptr);
            n_kmers = pu ? pu->core_size() : ef ? ef->core_size() : g.core_size();
            if (n_kmers == 0) {
                set_error(errbuf, errbuf_len, "no k-mers extracted from input");
                return nullptr;
            }
            pt.lap("fm: count k-mers + graph");
            text.reserve(n_kmers + n_kmers / 4);
            auto add = [&](const uint8_t* seq, uint32_t len, uint64_t) {
                tpoptoa::FmIndex::append_unitig(text, reinterpret_cast<const char*>(seq), len);
                ++n_units;
            };
            if (pu) pu->unitigs(add, nthreads);
            else if (ef) ef->unitigs(add, nthreads);
            else g.unitigs(add, nthreads);
        }
        pt.lap("fm: unitig walk");
        auto* h = new OkFmIndex();
        h->idx = new tpoptoa::FmIndex(K, std::move(text), n_units, n_kmers, nthreads);
        h->batch = h->idx->batch_pays();
        pt.lap("fm: FM index (k-order BWT, rank)");
        return h;
    } catch (const std::exception& e) {
        set_error(errbuf, errbuf_len, e.what());
        return nullptr;
    } catch (...) {
        set_error(errbuf, errbuf_len, "build failed (allocation error)");
        return nullptr;
    }
}

/// okfm_load() - read and validate a saved index (FmIndex::load()).
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

/// okfm_save() - write the index to @path, returning 0 on success and -1 on an I/O error.
int okfm_save(const OkFmIndex* h, const char* path,
              char* errbuf, uint64_t errbuf_len) {
    try {
        if (!h->idx->save(path)) {
            set_error(errbuf, errbuf_len, "save failed (I/O error)");
            return -1;
        }
        return 0;
    } catch (...) {
        set_error(errbuf, errbuf_len, "save failed (I/O error)");
        return -1;
    }
}

void okfm_free(OkFmIndex* h) { if (h) { delete h->idx; delete h; } }

unsigned    okfm_k(const OkFmIndex* h)           { return h->idx->k(); }
uint64_t    okfm_size(const OkFmIndex* h)        { return h->idx->size(); }
uint64_t    okfm_num_unitigs(const OkFmIndex* h) { return h->idx->num_unitigs(); }
uint64_t    okfm_super_len(const OkFmIndex* h)   { return h->idx->super_len(); }
uint64_t    okfm_bwt_len(const OkFmIndex* h)     { return h->idx->bwt_len(); }
uint64_t    okfm_bytes(const OkFmIndex* h)       { return h->idx->bytes(); }
double      okfm_bits_per_entry(const OkFmIndex* h) { return h->idx->bits_per_entry(); }

/// okfm_scan_seq() - look up every k-mer of one record (non-ACGT characters restart the k-mer), one FFI call per record.
///
/// num_searches counts full index lookups, as SSHash counts them (a lookup that is not a streaming extension).
/// This scan has no streaming extension, so every valid k-mer is one full lookup.
void okfm_scan_seq(const OkFmIndex* h, const char* seq, uint64_t slen, uint64_t k,
                   uint64_t* num_valid, uint64_t* num_positive, uint64_t* num_searches) {
    // The k-mers are collected and looked up kScanChunk at a time, through find_batch() when the index is larger than the last-level cache, which overlaps their memory latency.
    // The buffers are per thread and reused across records.
    constexpr std::size_t kScanChunk = 4096;
    thread_local std::vector<tpoptoa::KmerWord> kms;
    thread_local std::vector<uint8_t> hits;
    kms.clear();
    uint64_t nvalid = 0, npos = 0;
    auto flush = [&] {
        hits.resize(kms.size());
        if (h->batch) {
            h->idx->find_batch(kms.data(), kms.size(), hits.data());
        } else {
            for (std::size_t i = 0; i < kms.size(); ++i) hits[i] = h->idx->find(kms[i]);
        }
        for (uint8_t x : hits) npos += x;
        kms.clear();
    };
    tpoptoa::extract_kmers(seq, static_cast<std::size_t>(slen),
                           static_cast<std::size_t>(k),
                           [&](tpoptoa::KmerWord km) {
        ++nvalid;
        kms.push_back(km);
        if (kms.size() == kScanChunk) flush();
    });
    if (!kms.empty()) flush();
    *num_valid = nvalid;
    *num_positive = npos;
    *num_searches = nvalid;   // every valid k-mer is one full lookup
}

}  // extern "C"