#pragma once

#include "atomic_file.hpp"

/*
 * dbg_graph.hpp - compacted de Bruijn graph held in one open-addressing table.
 *
 * Slots.
 * 12 bytes per node, five per 64-byte line.
 * The line holds the five key words (canonical k-mer + 1, 0 = empty) and then the five node words (abundance, the walker's two flags, and the in << 4 | out edge mask).
 * Everything a walk step needs is in the node's own line, so one probe (linear probing from a fast-range home line, load at most 0.8) answers "is this neighbour present, and what are its edges".
 * There is no sorted key array, rank table or per-node walker array, so the whole build takes about 16 bytes per k-mer.
 *
 * Counts.
 * The node word keeps 22 bits of abundance, and the rare abundances of kCountMax and above are kept exactly in a sorted side list.
 *
 * Build.
 * The counter inserts each solid k-mer straight into its slot (src/kmer_count.cpp).
 * Every k-mer is counted in exactly one minimizer partition, so one thread owns it, with a compare-and-swap on the key word followed by a plain store of the node word.
 * Observed edges go into the mask byte, and resolve_masks() then restricts them in place to present neighbours.
 *
 * Walk.
 * dbg_unitigs.hpp extracts the maximal unitigs in linear time without locks.
 * Slot positions depend on insertion timing, so every tie-break uses k-mer values, never slot indices.
 * Sorting the output by head k-mer makes it independent of the thread count.
 *
 * Lookups compare full keys, so membership is exact.
 *
 * Sources are the bidirected node-and-mask model of Conway and Bromage 2011 and Cuttlefish 2 (Khan et al. 2022), linear probing as in Knuth, TAOCP vol. 3, 6.4, and fast range reduction from Lemire 2019.
 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include <sys/mman.h>

#include "hugepage.hpp"
#include "kmer.hpp"
#include "thread_budget.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

struct DbgUnitigExtractor;                      // defined in dbg_unitigs.hpp

class DbgGraph {
public:
    static constexpr std::size_t kNpos = static_cast<std::size_t>(-1);
    static constexpr unsigned kSlotsPerLine = 5;
    /// Node word, with bits 0..21 abundance (kCountMax meaning look in the side list), bits 22..23 walker flags and bits 24..31 the edge mask in << 4 | out, in the canonical k-mer's frame.
    static constexpr uint32_t kCountMax = (uint32_t(1) << 22) - 1;
    static constexpr unsigned kFlagShift = 22, kMaskShift = 24;
    static constexpr unsigned kWalkFlags = 3u;            // the two walker flag bits, unshifted
    /// Edge mask with all 8 edges set, where every (k-1)-overlap is a candidate.
    static constexpr uint8_t kAllEdges = 0xFF;
    struct Line {
        uint64_t keys[kSlotsPerLine];    // canonical k-mer + 1, or 0 for an empty slot
        uint32_t nodes[kSlotsPerLine];   // node words
        uint32_t pad;
    };
    static_assert(sizeof(Line) == 64, "one line = five slots = one cache line");

    explicit DbgGraph(unsigned kk) : k_(kk) {}

    unsigned    k() const noexcept { return k_; }
    std::size_t core_size() const noexcept { return n_; }
    std::size_t capacity() const noexcept { return cap_; }

    /// lines_for() and table_bytes() - lines and bytes of the table init_table(@n) allocates, n * 1.25 + 16 slots for a load of at most 0.8.
    static constexpr std::size_t lines_for(std::size_t n) noexcept { return (n + n / 4 + 16 + kSlotsPerLine - 1) / kSlotsPerLine; }
    static constexpr uint64_t table_bytes(std::size_t n) noexcept { return uint64_t(lines_for(n)) * sizeof(Line); }

    /// init_table() - allocate an empty table for @n nodes.
    ///
    /// It uses anonymous mmap, so zero pages come from the OS on first touch and the allocation is O(1), and the range is marked MADV_HUGEPAGE so random probes of a multi-GB table hit the TLB instead of paying a page walk (hugepage.hpp).
    void init_table(std::size_t n) {
        nlines_ = lines_for(n);
        cap_ = nlines_ * kSlotsPerLine;
        lines_.reset(alloc_zeroed(nlines_), map_len(nlines_));
        n_ = 0;
        saturated_.clear();
    }

    /// rank() - slot of canonical k-mer @c, or kNpos.
    std::size_t rank(uint64_t c) const noexcept {
        if (!cap_) return kNpos;
        const uint64_t want = c + 1;
        for (std::size_t L = home(c);;) {
            const Line& ln = lines_[L];
            for (unsigned s = 0; s < kSlotsPerLine; ++s) {
                const uint64_t v = ln.keys[s];
                if (v == want) return L * kSlotsPerLine + s;
                if (v == 0) return kNpos;
            }
            if (++L == nlines_) L = 0;
        }
    }
    /// prefetch() and prefetch_w() - prefetch the home line of @c for reading or writing.
    void prefetch(uint64_t c) const noexcept { if (cap_) __builtin_prefetch(&lines_[home(c)]); }
    void prefetch_w(uint64_t c) const noexcept { if (cap_) __builtin_prefetch(&lines_[home(c)], 1); }

    /// claim() - slot of canonical k-mer @c, inserting it if absent.
    ///
    /// Thread-safe.
    /// The node count is not maintained here, since a shared counter would bounce one cache line between all inserting threads, so callers insert distinct keys and report the total with finish_insert().
    std::size_t claim(uint64_t c) noexcept {
        const uint64_t want = c + 1;
        for (std::size_t L = home(c);;) {
            Line& ln = lines_[L];
            for (unsigned s = 0; s < kSlotsPerLine; ++s) {
                uint64_t expect = 0;
                if (__atomic_compare_exchange_n(&ln.keys[s], &expect, want, false,
                                                __ATOMIC_RELAXED, __ATOMIC_RELAXED))
                    return L * kSlotsPerLine + s;
                if (expect == want) return L * kSlotsPerLine + s;
            }
            if (++L == nlines_) L = 0;
        }
    }
    /// put() - store the abundance @cnt and observed edges @raw of claimed slot @i.
    /// Thread-safe for distinct slots.
    void put(std::size_t i, uint32_t cnt, uint8_t raw) {
        if (cnt >= kCountMax) {
            std::lock_guard<std::mutex> l(saturated_mutex_);
            saturated_.push_back({i, cnt});
        }
        node_word(i) = std::min(cnt, kCountMax) | (uint32_t(raw) << kMaskShift);
    }
    /// finish_insert() - record that @n distinct nodes were claimed, after all inserts.
    void finish_insert(std::size_t n) {
        n_ = n;
        std::sort(saturated_.begin(), saturated_.end());
    }

    // Per-slot accessors.
    // Node words are read with relaxed atomic loads, because walkers OR their flags into the same word concurrently (these compile to plain loads on x86).
    bool     used(std::size_t i) const noexcept { return key_word(i) != 0; }
    uint64_t key(std::size_t i)  const noexcept { return key_word(i) - 1; }
    unsigned edge_mask(std::size_t i)  const noexcept { return node(i) >> kMaskShift; }
    uint64_t count(std::size_t i) const noexcept {
        const uint32_t c = node(i) & kCountMax;
        if (c < kCountMax) return c;
        const auto it = std::lower_bound(saturated_.begin(), saturated_.end(),
                                         std::pair<std::size_t, uint32_t>{i, 0});
        return it != saturated_.end() && it->first == i ? it->second : c;
    }

    /// flags_or(), flags(), clear_flags() - the walker flags of slot @i (2 bits, zero between walks), as an atomic OR returning the previous flags, a relaxed read, and a reset by the slot's one writer.
    unsigned flags_or(std::size_t i, unsigned f) const noexcept {
        return (__atomic_fetch_or(&node_word(i), uint32_t(f) << kFlagShift, __ATOMIC_ACQ_REL) >> kFlagShift) & kWalkFlags;
    }
    unsigned flags(std::size_t i) const noexcept { return (node(i) >> kFlagShift) & kWalkFlags; }
    void clear_flags(std::size_t i) const noexcept { node_word(i) &= ~(uint32_t(kWalkFlags) << kFlagShift); }

    /// aux_clear(), aux_set(), aux() - reuse of the node words after the graph's last walk.
    ///
    /// The GFA writer stores segment ids there for link lookup instead of building a separate index.
    /// aux_clear() zeroes every node word and marks the graph consumed, since abundances, masks and walker flags are gone, and save(), unitigs() and export_sorted() refuse it afterwards.
    void aux_clear(unsigned nthreads) {
        consumed_ = true;
        const std::size_t nlines = nlines_;
        const unsigned threads = nthreads ? nthreads : 1;
        tpoptoa::run_on_threads(threads, [&](unsigned t) {
            for (std::size_t l = nlines * t / threads; l < nlines * (t + 1) / threads; ++l)
                for (unsigned s = 0; s < kSlotsPerLine; ++s) lines_[l].nodes[s] = 0;
        });
    }
    void aux_set(std::size_t i, uint32_t v) noexcept { __atomic_store_n(&node_word(i), v, __ATOMIC_RELAXED); }
    uint32_t aux(std::size_t i) const noexcept { return node(i); }
    bool consumed() const noexcept { return consumed_; }

    /// out_mask() and in_mask() - edge masks of canonical k-mer @c in its canonical frame, or 0 if it is absent.
    unsigned out_mask(uint64_t c) const noexcept {
        const std::size_t r = rank(c);
        return r == kNpos ? 0u : edge_mask(r) & 0xFu;
    }
    unsigned in_mask(uint64_t c) const noexcept {
        const std::size_t r = rank(c);
        return r == kNpos ? 0u : edge_mask(r) >> 4;
    }

    /// masks_from_raw() - no-op for the interface shared with DbgGraphEF and PartUnitigs.
    /// It is called when every observed edge's neighbour is known to be present, so the stored edges already are the masks.
    void masks_from_raw() noexcept {}

    /// resolve_masks() - restrict every node's observed edges to present neighbours, in place.
    ///
    /// Software pipeline (Slotin, Algorithms for Modern Hardware).
    /// The neighbour keys of block b + 1 are computed and their home lines prefetched while block b is probed, so each probe's cache miss was issued one block earlier.
    /// Blocks are small enough that the prefetches in flight fit a core's line-fill buffers.
    /// A block reads its own observed edges before any of its masks is rewritten, and probes read only key words.
    void resolve_masks() {
        const unsigned K = k_;
        const uint64_t kmask = kmer_mask(K);
        const unsigned shK = 2u * (K - 1);
        auto canon = [&](uint64_t v) { const uint64_t r = reverse_complement(v, K); return v < r ? v : r; };
        constexpr std::size_t kBlock = 4, kSpan = 4096;     // slots per prefetch block / per OpenMP task
        constexpr unsigned kEdges = 8;                     // 4 out-edges, then 4 in-edges
        const std::size_t nspan = (cap_ + kSpan - 1) / kSpan;
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 4)
#endif
        for (std::size_t sp = 0; sp < nspan; ++sp) {
            const std::size_t lo = sp * kSpan, hi = std::min(cap_, lo + kSpan);
            // neighbour keys of the current and next block, with ~0 for no edge
            uint64_t neighbours[2][kBlock][kEdges];
            auto prep = [&](std::size_t b0, uint64_t (*t)[kEdges]) {
                for (std::size_t i = b0; i < std::min(hi, b0 + kBlock); ++i) {
                    uint64_t* r = t[i - b0];
                    for (unsigned b = 0; b < kEdges; ++b) r[b] = ~uint64_t(0);
                    if (!used(i)) continue;
                    const uint64_t c = key(i);
                    const unsigned raw = edge_mask(i);
                    for (unsigned b = 0; b < 4; ++b) {
                        if (raw & (1u << b)) {
                            r[b] = canon(((c << 2) | b) & kmask);
                            prefetch(r[b]);
                        }
                        if (raw & (1u << (b + 4))) {
                            r[b + 4] = canon(((uint64_t(b) << shK) | (c >> 2)) & kmask);
                            prefetch(r[b + 4]);
                        }
                    }
                }
            };
            prep(lo, neighbours[0]);
            for (std::size_t b0 = lo, cur = 0; b0 < hi; b0 += kBlock, cur ^= 1) {
                if (b0 + kBlock < hi) prep(b0 + kBlock, neighbours[cur ^ 1]);
                for (std::size_t i = b0; i < std::min(hi, b0 + kBlock); ++i) {
                    if (!used(i)) continue;
                    unsigned m = 0;
                    for (unsigned b = 0; b < kEdges; ++b)
                        if (neighbours[cur][i - b0][b] != ~uint64_t(0) && rank(neighbours[cur][i - b0][b]) != kNpos) m |= 1u << b;
                    node_word(i) = (node_word(i) & ((uint32_t(1) << kMaskShift) - 1)) | (uint32_t(m) << kMaskShift);
                }
            }
        }
    }

    /// unitigs() - cb(seq, len, coverage) for every maximal unitig in head k-mer order, the same for any thread count.
    /// dbg_unitigs.hpp defines it.
    template <typename CB>
    void unitigs(CB&& cb, unsigned nthreads = 0) const;
    /// unitigs_sorted() - batch form, calling sink(n, at) with at(i) -> {seq, len, cov} of the i-th unitig in output order (DbgUnitigExtractor::run_sorted).
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned nthreads = 0) const;

    /// export_sorted() - (canonical k-mer, abundance) pairs in ascending order, and their edge masks.
    void export_sorted(std::vector<std::pair<uint64_t, uint32_t>>& kv,
                       std::vector<uint8_t>& masks) const {
        if (consumed_) throw std::logic_error("DbgGraph: node words were reused (aux_clear)");
        std::vector<uint64_t> keys;
        keys.reserve(n_);
        for (std::size_t i = 0; i < cap_; ++i)
            if (used(i)) keys.push_back(key(i));
        std::sort(keys.begin(), keys.end());
        kv.resize(keys.size());
        masks.resize(keys.size());
        for (std::size_t i = 0; i < keys.size(); ++i) {
            const std::size_t r = rank(keys[i]);
            kv[i] = {keys[i], static_cast<uint32_t>(count(r))};
            masks[i] = static_cast<uint8_t>(edge_mask(r));
        }
    }

    /// save() - write the little-endian kdbg archive (dbruijn -s)
    ///   "OKKDBG" | u8 k | u64 n | n x u64 sorted canonical k-mers |
    ///   n x u32 abundances | n x u8 masks (in << 4 | out)
    /// Returns false on an I/O error or a consumed graph.
    bool save(const char* path) const;

    /// build_core() - build from given nodes and masks, with @pairs sorted by canonical k-mer and duplicate-free, masks[i] = in << 4 | out of pairs[i], already restricted to present neighbours.
    /// @pairs is released.
    void build_core(std::vector<std::pair<uint64_t, uint32_t>>&& pairs,
                    const std::vector<uint8_t>& masks) {
        if (masks.size() != pairs.size())
            throw std::invalid_argument("dbg: masks/pairs size mismatch");
        for (std::size_t i = 1; i < pairs.size(); ++i)
            if (pairs[i - 1].first >= pairs[i].first)
                throw std::invalid_argument("dbg: k-mers must be sorted and unique");
        init_table(pairs.size());
        for (std::size_t i = 0; i < pairs.size(); ++i)
            put(claim(pairs[i].first), pairs[i].second, masks[i]);
        finish_insert(pairs.size());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(pairs);
    }

    /// build_core() - node-centric build, where every (k-1)-overlap between two present k-mers is an edge, for k-mer sets without read windows (DynSpss::compact).
    /// @pairs must be sorted by canonical k-mer and duplicate-free.
    void build_core(std::vector<std::pair<uint64_t, uint32_t>>&& pairs) {
        init_table(pairs.size());
        for (const auto& kv : pairs)                // every overlap is a candidate edge
            put(claim(kv.first), kv.second, kAllEdges);
        finish_insert(pairs.size());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(pairs);
        resolve_masks();
    }

private:
    friend struct DbgUnitigExtractor;

    /// Seed of the home-line hash (2^64 / golden ratio).
    static constexpr uint64_t kHomeSeed = 0x9E3779B97F4A7C15ULL;

    /// home() - home line of @c: its hash mapped onto [0, nlines_) by fast range reduction (Lemire 2019).
    std::size_t home(uint64_t c) const noexcept {
        return static_cast<std::size_t>(
            (static_cast<__uint128_t>(mix64(c, kHomeSeed)) * nlines_) >> 64);
    }
    uint64_t& key_word(std::size_t i) const noexcept { return lines_[i / kSlotsPerLine].keys[i % kSlotsPerLine]; }
    uint32_t& node_word(std::size_t i) const noexcept { return lines_[i / kSlotsPerLine].nodes[i % kSlotsPerLine]; }
    uint32_t node(std::size_t i) const noexcept { return __atomic_load_n(&node_word(i), __ATOMIC_RELAXED); }

    static std::size_t map_len(std::size_t nlines) noexcept {
        return (nlines * sizeof(Line) + kHugePageBytes - 1) & ~(kHugePageBytes - 1);
    }
    static Line* alloc_zeroed(std::size_t nlines) {
        const std::size_t len = map_len(nlines);
        void* p = ::mmap(nullptr, len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (p == MAP_FAILED) throw std::bad_alloc();
        advise_huge(p, len);
        return static_cast<Line*>(p);
    }
    /// LinePtr - owner of the table's mmap.
    struct LinePtr {
        Line* p = nullptr;
        std::size_t len = 0;
        LinePtr() = default;
        LinePtr(const LinePtr&) = delete;
        LinePtr& operator=(const LinePtr&) = delete;
        ~LinePtr() { if (p) ::munmap(p, len); }
        void reset(Line* q, std::size_t l) { if (p) ::munmap(p, len); p = q; len = l; }
        Line& operator[](std::size_t i) const noexcept { return p[i]; }
    };

    unsigned k_;
    std::size_t nlines_ = 0;
    std::size_t cap_ = 0;             // nlines_ * kSlotsPerLine slots
    bool consumed_ = false;           // node words reused by aux_clear
    std::size_t n_ = 0;               // nodes, i.e. claimed slots
    LinePtr lines_;
    std::vector<std::pair<std::size_t, uint32_t>> saturated_;   // (slot, abundance >= kCountMax), sorted
    std::mutex saturated_mutex_;
};

inline bool DbgGraph::save(const char* path) const {
    if (consumed_) return false;
    std::vector<std::pair<uint64_t, uint32_t>> kv;
    std::vector<uint8_t> masks;
    export_sorted(kv, masks);
    return atomic_write(path, [&](std::ofstream& os) {
        os.write("OKKDBG", 6);
        os.put(static_cast<char>(k_));
        const uint64_t n = kv.size();
        os.write(reinterpret_cast<const char*>(&n), sizeof(n));
        for (const auto& p : kv) os.write(reinterpret_cast<const char*>(&p.first), 8);
        for (const auto& p : kv) os.write(reinterpret_cast<const char*>(&p.second), 4);
        os.write(reinterpret_cast<const char*>(masks.data()), static_cast<std::streamsize>(n));
    });
}

}

#include "dbg_unitigs.hpp"
