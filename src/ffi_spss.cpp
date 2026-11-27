/**
 * @file
 * @brief C ABI of the SPSS index behind build-spss, spss-update and query.
 * okspss_* functions cover the static index, and okdyn_* its editable form.
 *
 * An index file is the 8-byte magic "OKSPSS", the SpssCore::save() blob and an 8-byte checksum trailer.
 * Padding the magic to 8 bytes keeps every array in the blob 8-byte aligned, so a mapped file can be used in place.
 */
#include <cstdio>
#include <cstring>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>

#include "../include/fasta_reader.hpp"
#include "../include/hugepage.hpp"
#include "../include/dyn_spss.hpp"
#include "kmer_count.hpp"
#include "../include/spss_core.hpp"
#include "../include/dbg_graph_ef.hpp"
#include "../include/dbg_part.hpp"
#include "../include/atomic_file.hpp"
#include "../include/dbg_graph.hpp"
#include "../include/wyhash.hpp"
#include "ffi_error.hpp"
#include <algorithm>

using tpoptoa::ffi::set_error;

namespace {
constexpr std::size_t kMagicBytes = 8;

using tpoptoa::ffi::kStatusFailed;
using tpoptoa::ffi::kStatusOk;
/// okspss_ids result for a window that is absent or contains a non-ACGT base.
constexpr int64_t kAbsentId = -1;
constexpr std::size_t kChecksumBytes = 8;
constexpr char kSpssMagic[kMagicBytes] = {'O', 'K', 'S', 'P', 'S', 'S', 0, 0};

/// Detects bit rot and truncation, with no cryptographic claim.
///
/// Blocks are hashed in parallel, folding 16 bytes per wymix() step, and the block hashes are then folded in order.
/// Since it's part of the file format, changing it would invalidate every saved index.
uint64_t blob_checksum(const uint8_t* data, std::size_t num_bytes) {
    constexpr std::size_t kStepBytes = 2 * sizeof(uint64_t);   // two words per wymix()
    constexpr std::size_t kBlock = std::size_t(1) << 20;       // 1 MiB, the unit of parallel work
    const std::size_t num_blocks = (num_bytes + kBlock - 1) / kBlock;
    std::vector<uint64_t> block_hash(num_blocks);
    const uint64_t* secret = tpoptoa::kWyhashDefaultSecret;

    #pragma omp parallel for schedule(static)
    for (std::size_t b = 0; b < num_blocks; ++b) {
        const uint8_t* block_data = data + b * kBlock;
        const std::size_t len = std::min(kBlock, num_bytes - b * kBlock);
        // The block index seeds each digest, tying it to its position.
        uint64_t hash = secret[0] ^ b, lo_word, hi_word;
        std::size_t pos = 0;

        for (; pos + kStepBytes <= len; pos += kStepBytes) {
            std::memcpy(&lo_word, block_data + pos, sizeof(uint64_t)); std::memcpy(&hi_word, block_data + pos + sizeof(uint64_t), sizeof(uint64_t));
            hash = tpoptoa::wymix(lo_word ^ secret[1] ^ hash, hi_word ^ secret[2]);
        }

        // The tail is zero-padded to a full word, and folding in its length keeps a short block from colliding with one ending in zero bytes.
        uint8_t tail[kStepBytes] = {0};

        std::memcpy(tail, block_data + pos, len - pos);
        std::memcpy(&lo_word, tail, sizeof(uint64_t)); std::memcpy(&hi_word, tail + sizeof(uint64_t), sizeof(uint64_t));
        block_hash[b] = tpoptoa::wymix(lo_word ^ secret[1] ^ hash, hi_word ^ secret[3] ^ len);
    }

    uint64_t hash = secret[0] ^ num_bytes;

    // Block order makes the checksum independent of thread count.
    for (uint64_t block_digest : block_hash) hash = tpoptoa::wymix(hash ^ secret[1], block_digest ^ secret[2]);

    return hash;
}

/// Appends the checksum and writes the file through atomic_file.hpp.
/// That way an interrupted write never leaves a truncated index, and -i X -o X never truncates the mapped input.
bool write_index(const char* path, std::vector<uint8_t>& blob) {
    const uint64_t checksum = blob_checksum(blob.data() + kMagicBytes, blob.size() - kMagicBytes);
    uint8_t bytes[sizeof(uint64_t)]; std::memcpy(bytes, &checksum, sizeof(uint64_t));

    blob.insert(blob.end(), bytes, bytes + sizeof(uint64_t));
    return tpoptoa::atomic_write(path, [&](std::ofstream& out) {
        out.write(reinterpret_cast<const char*>(blob.data()), static_cast<std::streamsize>(blob.size()));
    });
}

}

extern "C" {

struct OkSpss {
    tpoptoa::SpssCore core;
    void* map = nullptr;          // when set, core's arrays point into this mapping
    std::size_t map_len = 0;

    ~OkSpss() { if (map) ::munmap(map, map_len); }
};

/// With @p compact == 0 the input sequences are indexed as given.
/// That is exact for any input and already size-optimal when the input is an SPSS such as unitigs.
/// Otherwise the k-mers are counted under min_count and compacted into maximal unitigs first.
OkSpss* okspss_build_from_fasta(const char* path, uint64_t kmer_len, uint64_t minimizer_len,
                                int compact, uint64_t min_count, unsigned nthreads,
                                char* errbuf, uint64_t errbuf_len) {
    auto fail = [&](const char* msg) -> OkSpss* {
        set_error(errbuf, errbuf_len, msg);
        return nullptr;
    };

    try {
        tpoptoa::PhaseTimer timer;
        std::vector<std::string> seqs;

        if (compact) {
            tpoptoa::DbgGraph graph(static_cast<unsigned>(kmer_len));
            std::unique_ptr<tpoptoa::DbgGraphEF> ef_graph;       // set if the plain table won't fit
            std::unique_ptr<tpoptoa::PartUnitigs> part_unitigs;

            tpoptoa::detail::count_into_graph(path, kmer_len, static_cast<uint32_t>(min_count),
                                              nthreads, graph, &ef_graph,
                                              tpoptoa::part_mode_wanted(static_cast<unsigned>(kmer_len)) ? &part_unitigs : nullptr);
            if ((part_unitigs ? part_unitigs->core_size() : ef_graph ? ef_graph->core_size() : graph.core_size()) == 0)
                return fail("no k-mers extracted from input");

            auto add = [&](const uint8_t* seq, uint32_t len, uint64_t) {
                seqs.emplace_back(reinterpret_cast<const char*>(seq), len);
            };

            if (part_unitigs) part_unitigs->unitigs(add, nthreads);
            else if (ef_graph) ef_graph->unitigs(add, nthreads);
            else graph.unitigs(add, nthreads);
        } else {
            tpoptoa::iterate_inputs(path, [&](const tpoptoa::SeqRecord& rec) {
                seqs.push_back(rec.seq);
            });
        }

        timer.lap("spss: read input");

        // m = 0 lets SpssCore::choose_m pick a value from the input, and any m gives an exact index.
        // m = 0 asks for the estimate-driven choice.
        if (minimizer_len == 0) minimizer_len = tpoptoa::SpssCore::choose_m(seqs, static_cast<unsigned>(kmer_len));

        timer.lap("spss: choose m");
        auto* handle = new OkSpss();

        if (!handle->core.build(seqs, static_cast<unsigned>(kmer_len), static_cast<unsigned>(minimizer_len), nthreads, &seqs)) {
            delete handle;
            return fail("SPSS build failed (check 1 <= m < k <= 32)");
        }

        return handle;
    } catch (const std::exception& err) {
        return fail(err.what());
    }
}

void okspss_free(OkSpss* handle) { delete handle; }

int okspss_save(const OkSpss* handle, const char* path) {
    try {
        std::vector<uint8_t> blob(kSpssMagic, kSpssMagic + kMagicBytes);

        handle->core.save(blob);
        return write_index(path, blob) ? kStatusOk : kStatusFailed;
    } catch (...) {
        return kStatusFailed;
    }
}

/// Maps the file and uses its big arrays in place, so opening costs only page-table setup plus whatever pages queries touch.
/// Passing @p verify != 0 checks the checksum first, which reads the whole file once.
OkSpss* okspss_load(const char* path, int verify, char* errbuf, uint64_t errbuf_len) {
    auto fail = [&](const char* msg) -> OkSpss* {
        set_error(errbuf, errbuf_len, msg);
        return nullptr;
    };
    const int file_desc = ::open(path, O_RDONLY);

    if (file_desc < 0) return fail("cannot open SPSS index");

    struct stat info;

    if (::fstat(file_desc, &info) != 0 || info.st_size < static_cast<off_t>(kMagicBytes)) { ::close(file_desc); return fail("not an SPSS index"); }

    const std::size_t len = static_cast<std::size_t>(info.st_size);
    void* mapping = ::mmap(nullptr, len, PROT_READ, MAP_PRIVATE, file_desc, 0);

    ::close(file_desc);

    if (mapping == MAP_FAILED) return fail("mmap failed");

    tpoptoa::advise_huge(mapping, len);   // file-backed THP, where the kernel supports it
    const auto* bytes = static_cast<const uint8_t*>(mapping);

    if (std::memcmp(bytes, kSpssMagic, kMagicBytes) != 0) {
        ::munmap(mapping, len);
        return fail("not an optik SPSS index (bad magic); rebuild it with build-spss");
    }

    // The 8-byte checksum follows the body. A file shorter than magic plus checksum cannot hold one and is rejected as truncated.
    if (len < kMagicBytes + kChecksumBytes) { ::munmap(mapping, len); return fail("corrupt SPSS index (truncated)"); }

    const std::size_t body = len - kMagicBytes - kChecksumBytes;

    if (verify) {
        uint64_t want; std::memcpy(&want, bytes + len - kChecksumBytes, kChecksumBytes);

        if (blob_checksum(bytes + kMagicBytes, body) != want) {
            ::munmap(mapping, len);
            return fail("corrupt SPSS index (checksum mismatch: truncated or modified file)");
        }
    }

    auto* handle = new OkSpss();

    // SpssCore borrows its arrays from the mapping, which okspss_free() unmaps together with the handle.
    handle->map = mapping; handle->map_len = len;

    if (handle->core.load(bytes + kMagicBytes, body, /*borrow=*/true) != body) {
        delete handle;
        return fail("corrupt SPSS index");
    }

    return handle;
}

uint64_t okspss_kmers(const OkSpss* handle) { return handle->core.stats().kmers; }
uint64_t okspss_bytes(const OkSpss* handle) { return handle->core.bytes(); }
uint32_t okspss_k(const OkSpss* handle) { return handle->core.kmer_len(); }
uint32_t okspss_m(const OkSpss* handle) { return handle->core.minimizer_len(); }
uint64_t okspss_superkmers(const OkSpss* handle) { return handle->core.stats().superkmers; }

/// Counts distinct canonical k-mers.
/// okspss_kmers instead counts positions, which repeat whenever the input does, since build-spss without --compact indexes records as given.
/// A position is a k-mer's first occurrence exactly when looking up that k-mer returns it.
uint64_t okspss_distinct(const OkSpss* handle) { return handle->core.count_distinct(); }

/// Writes each window's dense id in [0, num_kmers) to out[i], or kAbsentId when the window is missing or invalid.
/// @p out needs len - k + 1 entries, and that count is returned.
uint64_t okspss_ids(const OkSpss* handle, const char* seq, uint64_t len, int64_t* out) {
    const unsigned kmer_len = handle->core.kmer_len();

    if (len < kmer_len) return 0;

    const uint64_t mask = tpoptoa::kmer_mask(kmer_len);
    uint64_t fwd = 0; uint64_t valid = 0;

    for (uint64_t i = 0; i < len; ++i) {
        const uint8_t base = tpoptoa::base_encode(seq[i]);

        if (base == tpoptoa::kInvalidBase) { valid = 0; fwd = 0; }
        else { fwd = ((fwd << tpoptoa::kBitsPerBase) | base) & mask; ++valid; }

        // Every window position gets an output slot, and a window containing N is left at the absent id.
        if (i + 1 >= kmer_len) {
            int64_t kmer_ident = kAbsentId;

            if (valid >= kmer_len) {
                const uint64_t found = handle->core.kmer_id(fwd);

                if (found != tpoptoa::SpssCore::kNotFound) kmer_ident = static_cast<int64_t>(found);
            }

            out[i + 1 - kmer_len] = kmer_ident;
        }
    }

    return len - kmer_len + 1;
}

void okspss_scan_seq(const OkSpss* handle, const char* seq, uint64_t len,
                     uint64_t* num_valid, uint64_t* num_positive,
                     uint64_t* num_searches, uint64_t* num_extensions) {
    tpoptoa::SpssCore::ScanStats stats;

    handle->core.scan(seq, static_cast<std::size_t>(len), stats);
    *num_valid = stats.valid;
    *num_positive = stats.positive;
    *num_searches = stats.searches;
    *num_extensions = stats.extensions;
}

/// Handles many records per call, so cache misses overlap across records (SpssCore::scan_many()).
void okspss_scan_many(const OkSpss* handle, const char* const* seqs, const uint64_t* lens,
                      uint64_t count, uint64_t* num_valid, uint64_t* num_positive,
                      uint64_t* num_searches, uint64_t* num_extensions) {
    std::vector<std::size_t> lens_sz(lens, lens + count);
    tpoptoa::SpssCore::ScanStats stats;

    handle->core.scan_many(seqs, lens_sz.data(), static_cast<std::size_t>(count), stats);
    *num_valid = stats.valid;
    *num_positive = stats.positive;
    *num_searches = stats.searches;
    *num_extensions = stats.extensions;
}

// spss-update edits a loaded index through DynSpss (dyn_spss.hpp) and compacts it into a new static core on save.

struct OkDynSpss {
    void* map = nullptr;             // the loaded core still points into this
    std::size_t map_len = 0;
    std::unique_ptr<tpoptoa::DynSpss> dyn;

    ~OkDynSpss() { dyn.reset(); if (map) ::munmap(map, map_len); }
};

/// Always verifies the checksum, because an edit would otherwise bake corruption into the new file.
OkDynSpss* okdyn_open(const char* path, char* errbuf, uint64_t errbuf_len) {
    OkSpss* loaded = okspss_load(path, /*verify=*/1, errbuf, errbuf_len);

    if (!loaded) return nullptr;

    auto* handle = new OkDynSpss();

    handle->map = loaded->map; handle->map_len = loaded->map_len;
    loaded->map = nullptr;                                  // the new handle owns it now
    handle->dyn = std::make_unique<tpoptoa::DynSpss>(
        std::make_unique<tpoptoa::SpssCore>(std::move(loaded->core)));
    delete loaded;
    return handle;
}

void okdyn_free(OkDynSpss* handle) { delete handle; }

/// Inserts (@p add != 0) or erases every k-mer of the file.
/// It returns how many k-mers changed membership, or kStatusFailed with a message in @p errbuf.
int64_t okdyn_apply_fasta(OkDynSpss* handle, const char* path, int add,
                          char* errbuf, uint64_t errbuf_len) {
    try {
        int64_t changed = 0;
        const unsigned kmer_len = handle->dyn->kmer_len();

        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& rec) {
            uint64_t fwd = 0; std::size_t run_len = 0;
            const uint64_t mask = tpoptoa::kmer_mask(kmer_len);

            for (char ch : rec.seq) {
                const uint8_t base = tpoptoa::base_encode(ch);

                if (base == tpoptoa::kInvalidBase) { run_len = 0; fwd = 0; continue; }

                fwd = ((fwd << tpoptoa::kBitsPerBase) | base) & mask;

                // Only k-mers whose membership flips are counted.
                if (++run_len >= kmer_len) changed += add ? handle->dyn->insert(fwd) : handle->dyn->erase(fwd);
            }
        });
        return changed;
    } catch (const std::exception& err) {
        set_error(errbuf, errbuf_len, err.what());
        return kStatusFailed;
    }
}

uint64_t okdyn_added(const OkDynSpss* handle) { return handle->dyn->added(); }
uint64_t okdyn_deleted(const OkDynSpss* handle) { return handle->dyn->deleted(); }

/// Folds the edits into a fresh static core before writing.
int okdyn_compact_save(OkDynSpss* handle, const char* path, unsigned nthreads) {
    try {
        if (!handle->dyn->compact(nthreads)) return kStatusFailed;

        std::vector<uint8_t> blob(kSpssMagic, kSpssMagic + kMagicBytes);

        handle->dyn->core().save(blob);
        return write_index(path, blob) ? kStatusOk : kStatusFailed;
    } catch (...) {
        return kStatusFailed;
    }
}

uint64_t okdyn_kmers(const OkDynSpss* handle) { return handle->dyn->core().stats().kmers; }

}  // extern "C"
