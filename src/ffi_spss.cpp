/*
 * ffi_spss.cpp - C ABI of the SPSS index (build-spss, spss-update, query), with okspss_* for the static index, okdyn_* for its dynamic form.
 *
 * The index file holds the 8-byte magic "OKSPSS", the SpssCore::save() blob and an 8-byte checksum trailer.
 * The 8-byte magic keeps every array of the blob 8-byte aligned when the file is mapped, so the arrays are used in place.
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
// The 8-byte magic, padded with zeros so the blob after it stays 8-byte aligned.
constexpr std::size_t kMagicBytes = 8;
constexpr std::size_t kChecksumBytes = 8;
constexpr char kSpssMagic[kMagicBytes] = {'O', 'K', 'S', 'P', 'S', 'S', 0, 0};

/// blob_checksum() - corruption checksum of @n bytes at @p.
///
/// 1 MiB blocks are hashed in parallel (wymix() 128-bit multiply-fold per 16 bytes, wyhash.hpp) and the block hashes folded in order.
/// Not cryptographic, and it detects bit rot and truncation.
/// It is part of the file format, so it must not change.
uint64_t blob_checksum(const uint8_t* p, std::size_t n) {
    constexpr std::size_t kBlock = std::size_t(1) << 20;
    const std::size_t nb = (n + kBlock - 1) / kBlock;
    std::vector<uint64_t> block_hash(nb);
    const uint64_t* S = tpoptoa::kWyhashDefaultSecret;
    #pragma omp parallel for schedule(static)
    for (std::size_t b = 0; b < nb; ++b) {
        const uint8_t* q = p + b * kBlock;
        const std::size_t len = std::min(kBlock, n - b * kBlock);
        uint64_t h = S[0] ^ b, w0, w1;
        std::size_t i = 0;
        for (; i + 16 <= len; i += 16) {
            std::memcpy(&w0, q + i, 8); std::memcpy(&w1, q + i + 8, 8);
            h = tpoptoa::wymix(w0 ^ S[1] ^ h, w1 ^ S[2]);
        }
        uint8_t tail[16] = {0};
        std::memcpy(tail, q + i, len - i);
        std::memcpy(&w0, tail, 8); std::memcpy(&w1, tail + 8, 8);
        block_hash[b] = tpoptoa::wymix(w0 ^ S[1] ^ h, w1 ^ S[3] ^ len);
    }
    uint64_t h = S[0] ^ n;
    for (uint64_t x : block_hash) h = tpoptoa::wymix(h ^ S[1], x ^ S[2]);
    return h;
}

/// write_index() - append the checksum to @blob, the magic followed by the SpssCore blob, and write it to @path atomically through atomic_file.hpp.
/// An interrupted write then never leaves a truncated index, and an in-place update with -i X -o X never truncates the mapped input.
bool write_index(const char* path, std::vector<uint8_t>& blob) {
    const uint64_t cs = blob_checksum(blob.data() + kMagicBytes, blob.size() - kMagicBytes);
    uint8_t b[8]; std::memcpy(b, &cs, 8);
    blob.insert(blob.end(), b, b + 8);
    return tpoptoa::atomic_write(path, [&](std::ofstream& os) {
        os.write(reinterpret_cast<const char*>(blob.data()), static_cast<std::streamsize>(blob.size()));
    });
}

}

extern "C" {

/// OkSpss - a static index, built in memory or mapped from a file.
struct OkSpss {
    tpoptoa::SpssCore core;
    void* map = nullptr;          // zero-copy mapping backing core's arrays
    std::size_t map_len = 0;
    ~OkSpss() { if (map) ::munmap(map, map_len); }
};

/// okspss_build_from_fasta() - build from a FASTA/FASTQ.
/// With @compact == 0 the input sequences are indexed as given, which is exact for any input and size-optimal when the input is already an SPSS such as unitigs.
/// With compact != 0 the k-mers are counted under the min_count filter and compacted into maximal unitigs first.
OkSpss* okspss_build_from_fasta(const char* path, uint64_t k, uint64_t m,
                                int compact, uint64_t min_count, unsigned nthreads,
                                char* errbuf, uint64_t errbuf_len) {
    auto fail = [&](const char* msg) -> OkSpss* {
        set_error(errbuf, errbuf_len, msg);
        return nullptr;
    };
    try {
        tpoptoa::PhaseTimer pt;                         // OPTIK_TIMING=1
        std::vector<std::string> seqs;
        if (compact) {
            tpoptoa::DbgGraph g(static_cast<unsigned>(k));
            std::unique_ptr<tpoptoa::DbgGraphEF> ef;       // when the node table exceeds the limit
            std::unique_ptr<tpoptoa::PartUnitigs> pu;      // partitioned compaction
            tpoptoa::detail::count_into_graph(path, k, static_cast<uint32_t>(min_count),
                                              nthreads, g, &ef,
                                              tpoptoa::part_mode_wanted(static_cast<unsigned>(k)) ? &pu : nullptr);
            if ((pu ? pu->core_size() : ef ? ef->core_size() : g.core_size()) == 0)
                return fail("no k-mers extracted from input");
            auto add = [&](const uint8_t* s, uint32_t len, uint64_t) {
                seqs.emplace_back(reinterpret_cast<const char*>(s), len);
            };
            if (pu) pu->unitigs(add, nthreads);
            else if (ef) ef->unitigs(add, nthreads);
            else g.unitigs(add, nthreads);
        } else {
            tpoptoa::iterate_inputs(path, [&](const tpoptoa::SeqRecord& rec) {
                seqs.push_back(rec.seq);
            });
        }
        pt.lap("spss: read input");
        // m = 0 lets SpssCore::choose_m pick it from the input, exact at any m
        if (m == 0) m = tpoptoa::SpssCore::choose_m(seqs, static_cast<unsigned>(k));
        pt.lap("spss: choose m");
        auto* h = new OkSpss();
        if (!h->core.build(seqs, static_cast<unsigned>(k), static_cast<unsigned>(m), nthreads, &seqs)) {
            delete h;
            return fail("SPSS build failed (check 1 <= m < k <= 32)");
        }
        return h;
    } catch (const std::exception& e) {
        return fail(e.what());
    }
}

void okspss_free(OkSpss* h) { delete h; }

/// okspss_save() - write the index in the current format, returning 0 on success.
int okspss_save(const OkSpss* h, const char* path) {
    try {
        std::vector<uint8_t> blob(kSpssMagic, kSpssMagic + kMagicBytes);
        h->core.save(blob);
        return write_index(path, blob) ? 0 : -1;
    } catch (...) {
        return -1;
    }
}

/// okspss_load() - map an index file and use its big arrays in place, without a read loop or copy, so opening an index costs page-table setup plus the pages queries touch.
/// @verify != 0 first checks the checksum, which reads the whole file once.
OkSpss* okspss_load(const char* path, int verify, char* errbuf, uint64_t errbuf_len) {
    auto fail = [&](const char* msg) -> OkSpss* {
        set_error(errbuf, errbuf_len, msg);
        return nullptr;
    };
    const int fd = ::open(path, O_RDONLY);
    if (fd < 0) return fail("cannot open SPSS index");
    struct stat st;
    if (::fstat(fd, &st) != 0 || st.st_size < static_cast<off_t>(kMagicBytes)) { ::close(fd); return fail("not an SPSS index"); }
    const std::size_t len = static_cast<std::size_t>(st.st_size);
    void* m = ::mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (m == MAP_FAILED) return fail("mmap failed");
    tpoptoa::advise_huge(m, len);   // file THP where the kernel supports it
    const auto* p = static_cast<const uint8_t*>(m);
    if (std::memcmp(p, kSpssMagic, kMagicBytes) != 0) {
        ::munmap(m, len);
        return fail("not an optik SPSS index (bad magic); rebuild it with build-spss");
    }
    if (len < kMagicBytes + kChecksumBytes) { ::munmap(m, len); return fail("corrupt SPSS index (truncated)"); }
    const std::size_t body = len - kMagicBytes - kChecksumBytes;   // blob bytes after the magic
    if (verify) {
        uint64_t want; std::memcpy(&want, p + len - kChecksumBytes, kChecksumBytes);
        if (blob_checksum(p + kMagicBytes, body) != want) {
            ::munmap(m, len);
            return fail("corrupt SPSS index (checksum mismatch: truncated or modified file)");
        }
    }
    auto* h = new OkSpss();
    h->map = m; h->map_len = len;
    if (h->core.load(p + kMagicBytes, body, /*borrow=*/true) != body) {
        delete h;
        return fail("corrupt SPSS index");
    }
    return h;
}

uint64_t okspss_kmers(const OkSpss* h) { return h->core.stats().kmers; }
uint64_t okspss_bytes(const OkSpss* h) { return h->core.bytes(); }
uint32_t okspss_k(const OkSpss* h) { return h->core.k(); }
uint32_t okspss_m(const OkSpss* h) { return h->core.m(); }
uint64_t okspss_superkmers(const OkSpss* h) { return h->core.stats().superkmers; }

/// okspss_distinct() - distinct canonical k-mers.
/// okspss_kmers counts k-mer positions in the
/// spectrum-preserving string set, which repeats a k-mer when the input does (build-spss without --compact indexes records as given), where a position is the k-mer's first occurrence iff the lookup of its k-mer returns it.
uint64_t okspss_distinct(const OkSpss* h) { return h->core.count_distinct(); }

/// okspss_ids() - dictionary lookup of every k-mer window of @seq:
/// out[i] = dense k-mer
/// identifier in [0, num_kmers) of window i, or -1 when it is absent or holds a non-ACGT base.
/// out holds len - k + 1 entries, and that count is returned.
uint64_t okspss_ids(const OkSpss* h, const char* seq, uint64_t len, int64_t* out) {
    const unsigned k = h->core.k();
    if (len < k) return 0;
    const uint64_t mask = tpoptoa::kmer_mask(k);
    uint64_t f = 0; uint64_t valid = 0;
    for (uint64_t i = 0; i < len; ++i) {
        const uint8_t b = tpoptoa::base_encode(seq[i]);
        if (b == tpoptoa::kInvalidBase) { valid = 0; f = 0; }
        else { f = ((f << 2) | b) & mask; ++valid; }
        if (i + 1 >= k) {
            int64_t id = -1;
            if (valid >= k) {
                const uint64_t r = h->core.kmer_id(f);
                if (r != tpoptoa::SpssCore::kNotFound) id = static_cast<int64_t>(r);
            }
            out[i + 1 - k] = id;
        }
    }
    return len - k + 1;
}

/// okspss_scan_seq() - streaming lookup of every k-mer of one sequence (non-ACGT characters restart the k-mer).
void okspss_scan_seq(const OkSpss* h, const char* seq, uint64_t len,
                     uint64_t* num_valid, uint64_t* num_positive,
                     uint64_t* num_searches, uint64_t* num_extensions) {
    tpoptoa::SpssCore::ScanStats st;
    h->core.scan(seq, static_cast<std::size_t>(len), st);
    *num_valid = st.valid;
    *num_positive = st.positive;
    *num_searches = st.searches;
    *num_extensions = st.extensions;
}

/// okspss_scan_many() - streaming lookup over many records at once, so cache misses overlap across records (SpssCore::scan_many()).
void okspss_scan_many(const OkSpss* h, const char* const* seqs, const uint64_t* lens,
                      uint64_t n, uint64_t* num_valid, uint64_t* num_positive,
                      uint64_t* num_searches, uint64_t* num_extensions) {
    std::vector<std::size_t> l(lens, lens + n);
    tpoptoa::SpssCore::ScanStats st;
    h->core.scan_many(seqs, l.data(), static_cast<std::size_t>(n), st);
    *num_valid = st.valid;
    *num_positive = st.positive;
    *num_searches = st.searches;
    *num_extensions = st.extensions;
}

// Dynamic updates use the static core plus open-addressing edit sets (dyn_spss.hpp), compacted into a new static core on save.

/// OkDynSpss - an index opened for spss-update.
struct OkDynSpss {
    void* map = nullptr;             // backing mapping of the loaded core
    std::size_t map_len = 0;
    std::unique_ptr<tpoptoa::DynSpss> dyn;
    ~OkDynSpss() { dyn.reset(); if (map) ::munmap(map, map_len); }
};

/// okdyn_open() - load @path for editing after verifying its checksum.
OkDynSpss* okdyn_open(const char* path, char* errbuf, uint64_t errbuf_len) {
    OkSpss* s = okspss_load(path, /*verify=*/1, errbuf, errbuf_len);
    if (!s) return nullptr;
    auto* h = new OkDynSpss();
    h->map = s->map; h->map_len = s->map_len;
    s->map = nullptr;                                  // ownership moves to h
    h->dyn = std::make_unique<tpoptoa::DynSpss>(
        std::make_unique<tpoptoa::SpssCore>(std::move(s->core)));
    delete s;
    return h;
}

void okdyn_free(OkDynSpss* h) { delete h; }

/// okdyn_apply_fasta() - insert every k-mer of a FASTA/FASTQ when @add != 0, or erase them, returning the number of k-mers whose membership changed, or -1 with a message in @errbuf.
int64_t okdyn_apply_fasta(OkDynSpss* h, const char* path, int add,
                          char* errbuf, uint64_t errbuf_len) {
    try {
        int64_t changed = 0;
        const unsigned K = h->dyn->k();
        tpoptoa::iterate_sequences(path, [&](const tpoptoa::SeqRecord& rec) {
            uint64_t f = 0; std::size_t v = 0;
            const uint64_t mask = tpoptoa::kmer_mask(K);
            for (char c : rec.seq) {
                const uint8_t b = tpoptoa::base_encode(c);
                if (b == tpoptoa::kInvalidBase) { v = 0; f = 0; continue; }
                f = ((f << 2) | b) & mask;
                if (++v >= K) changed += add ? h->dyn->insert(f) : h->dyn->erase(f);
            }
        });
        return changed;
    } catch (const std::exception& e) {
        set_error(errbuf, errbuf_len, e.what());
        return -1;
    }
}

uint64_t okdyn_added(const OkDynSpss* h) { return h->dyn->added(); }
uint64_t okdyn_deleted(const OkDynSpss* h) { return h->dyn->deleted(); }

/// okdyn_compact_save() - fold the edits into a fresh static core and write it to @path.
/// Returns 0 on success.
int okdyn_compact_save(OkDynSpss* h, const char* path, unsigned nthreads) {
    try {
        if (!h->dyn->compact(nthreads)) return -1;
        std::vector<uint8_t> blob(kSpssMagic, kSpssMagic + kMagicBytes);
        h->dyn->core().save(blob);
        return write_index(path, blob) ? 0 : -1;
    } catch (...) {
        return -1;
    }
}

uint64_t okdyn_kmers(const OkDynSpss* h) { return h->dyn->core().stats().kmers; }

}  // extern "C"
