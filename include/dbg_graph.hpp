#pragma once

#include "atomic_file.hpp"

/**
 * @file
 * @brief Compacted de Bruijn graph held in one open-addressing table.
 *
 * dbg_graph.hpp - compacted de Bruijn graph held in one open-addressing table.
 *
 * @par Slots
 * Each node takes 12 bytes, and one 64-byte line holds five of them.
 * A line stores five key words (canonical k-mer + 1, with 0 meaning empty) followed by five node words.
 * Node words pack the abundance, the walker's two flags and the in << 4 | out edge mask.
 * Since everything a walk step needs sits in the node's own line, a single probe tells whether a neighbour exists and what its edges are.
 * Probing is linear from a fast-range home line, with the load kept at or below 0.8.
 * No sorted key array, rank table or per-node walker array exists, so a whole build costs about 16 bytes per k-mer.
 *
 * @par Counts
 * Node words keep 22 bits of abundance.
 * The rare abundances of kCountMax and above are stored exactly in a sorted side list.
 *
 * @par Build
 * The counter inserts each solid k-mer straight into its slot (src/kmer_count.cpp).
 * Each k-mer belongs to exactly one minimizer partition, so a single thread owns it.
 * That thread claims the key word with a compare-and-swap and then writes the node word with a plain store.
 * Observed edges go into the mask byte, and resolve_masks() later trims them to neighbours that exist.
 *
 * @par Walk
 * dbg_unitigs.hpp extracts the maximal unitigs in linear time without locks.
 * Slot positions depend on insertion timing, so every tie-break compares k-mer values rather than slot indices.
 * Sorting the output by head k-mer then makes it independent of the thread count.
 *
 * Lookups compare full keys keeping membership exact.
 *
 * The bidirected node-and-mask model follows Conway and Bromage 2011 and Cuttlefish 2 (Khan et al. 2022).
 * Linear probing is as it is in Knuth, TAOCP vol. 3, 6.4, and the fast range reduction comes from Lemire 2019.
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

#include "bits.hpp"
#include "hugepage.hpp"
#include "kmer.hpp"
#include "thread_budget.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

struct DbgUnitigExtractor;                      // defined in dbg_unitigs.hpp

class DbgGraph {
public:
    static constexpr std::size_t kNpos = static_cast<std::size_t>(-1);
    static constexpr std::size_t kCacheLineBytes = 64;   // x86-64, the unit one probe reads
    static constexpr unsigned kSlotsPerLine = 5;
    /// Node word layout, in the canonical k-mer's frame:
    ///   bits 0..21   abundance (kCountMax means "look in the side list")
    ///   bits 22..23  walker flags
    ///   bits 24..31  edge mask, in << 4 | out
    static constexpr unsigned kFlagShift = 22, kMaskShift = 24;
    static constexpr uint32_t kCountMax = (uint32_t(1) << kFlagShift) - 1;
    static constexpr unsigned kWalkFlags = 3u;            // the two walker flag bits, unshifted
    /// Every (k-1)-overlap is a candidate edge when all 8 bits are set.
    static constexpr uint8_t kAllEdges = 0xFF;
    static constexpr char kKdbgMagic[] = "OKKDBG";       // first bytes of a save() archive
    static constexpr std::size_t kKdbgMagicLen = sizeof(kKdbgMagic) - 1;
    static constexpr int kNoFd = -1;                      // anonymous mappings take fd -1

    struct Line {
        uint64_t keys[kSlotsPerLine];    // canonical k-mer + 1, or 0 for an empty slot
        uint32_t nodes[kSlotsPerLine];
        uint32_t pad;
    };
    static_assert(sizeof(Line) == kCacheLineBytes, "one line = five slots = one cache line");

    explicit DbgGraph(unsigned kmer_len) : kmer_len_(kmer_len) {}

    unsigned    kmer_len() const noexcept { return kmer_len_; }
    std::size_t core_size() const noexcept { return num_nodes_; }
    std::size_t capacity() const noexcept { return cap_; }

    /// The table holds the nodes plus a quarter (kSpareSlotDiv) and kSpareSlots more, which keeps the load at or below 0.8.
    static constexpr std::size_t kSpareSlotDiv = 4;
    static constexpr std::size_t kSpareSlots = 16;
    static constexpr std::size_t lines_for(std::size_t nodes) noexcept { return (nodes + nodes / kSpareSlotDiv + kSpareSlots + kSlotsPerLine - 1) / kSlotsPerLine; }
    static constexpr uint64_t table_bytes(std::size_t nodes) noexcept { return uint64_t(lines_for(nodes)) * sizeof(Line); }

    /// Anonymous mmap makes the allocation O(1), since zero pages only arrive on first touch.
    /// The range is also advised MADV_HUGEPAGE, so random probes into a multi-GB table hit the TLB instead of walking pages (hugepage.hpp).
    void init_table(std::size_t nodes) {
        nlines_ = lines_for(nodes);
        cap_ = nlines_ * kSlotsPerLine;
        lines_.reset(alloc_zeroed(nlines_), map_len(nlines_));
        num_nodes_ = 0;
        saturated_.clear();
    }

    /// Returns kNpos when the k-mer is absent.
    std::size_t rank(uint64_t canon_kmer) const noexcept {
        if (!cap_) return kNpos;

        const uint64_t want = canon_kmer + 1;

        for (std::size_t line_idx = home(canon_kmer);;) {
            const Line& line = lines_[line_idx];

            for (unsigned s = 0; s < kSlotsPerLine; ++s) {
                const uint64_t stored = line.keys[s];

                if (stored == want) return line_idx * kSlotsPerLine + s;

                if (stored == 0) return kNpos;
            }

            if (++line_idx == nlines_) line_idx = 0;
        }
    }
    void prefetch(uint64_t canon_kmer) const noexcept { if (cap_) __builtin_prefetch(&lines_[home(canon_kmer)]); }
    void prefetch_w(uint64_t canon_kmer) const noexcept { if (cap_) __builtin_prefetch(&lines_[home(canon_kmer)], 1); }

    /// Thread-safe, and inserts the k-mer if it's missing.
    /// The node count isn't updated here, because a shared counter would bounce one cache line between every inserting thread.
    /// Instead, callers insert distinct keys and report the total through finish_insert().
    std::size_t claim(uint64_t canon_kmer) noexcept {
        const uint64_t want = canon_kmer + 1;

        for (std::size_t line_idx = home(canon_kmer);;) {
            Line& line = lines_[line_idx];

            for (unsigned s = 0; s < kSlotsPerLine; ++s) {
                uint64_t expect = 0;

                if (__atomic_compare_exchange_n(&line.keys[s], &expect, want, false,
                                                __ATOMIC_RELAXED, __ATOMIC_RELAXED))
                    return line_idx * kSlotsPerLine + s;

                if (expect == want) return line_idx * kSlotsPerLine + s;
            }

            if (++line_idx == nlines_) line_idx = 0;
        }
    }
    /// Safe to call concurrently on distinct slots.
    void put(std::size_t slot, uint32_t cnt, uint8_t raw) {
        if (cnt >= kCountMax) {
            std::lock_guard<std::mutex> guard(saturated_mutex_);

            saturated_.push_back({slot, cnt});
        }

        node_word(slot) = std::min(cnt, kCountMax) | (uint32_t(raw) << kMaskShift);
    }
    void finish_insert(std::size_t nodes) {
        num_nodes_ = nodes;
        std::sort(saturated_.begin(), saturated_.end());
    }

    // Per-slot accessors.
    // Walkers OR their flags into node words concurrently, so reads use relaxed atomic loads, which are plain loads on x86 anyway.
    bool     used(std::size_t slot) const noexcept { return key_word(slot) != 0; }
    uint64_t key(std::size_t slot)  const noexcept { return key_word(slot) - 1; }
    unsigned edge_mask(std::size_t slot)  const noexcept { return node(slot) >> kMaskShift; }

    uint64_t count(std::size_t slot) const noexcept {
        const uint32_t stored = node(slot) & kCountMax;

        if (stored < kCountMax) return stored;

        const auto found = std::lower_bound(saturated_.begin(), saturated_.end(),
                                         std::pair<std::size_t, uint32_t>{slot, 0});

        return found != saturated_.end() && found->first == slot ? found->second : stored;
    }

    /// Each slot has 2 walker-flag bits, which are zero between walks.
    /// They support an atomic OR that returns the old value, a relaxed read, and a reset by the slot's single writer.
    unsigned flags_or(std::size_t slot, unsigned flag_bits) const noexcept {
        return (__atomic_fetch_or(&node_word(slot), uint32_t(flag_bits) << kFlagShift, __ATOMIC_ACQ_REL) >> kFlagShift) & kWalkFlags;
    }
    unsigned flags(std::size_t slot) const noexcept { return (node(slot) >> kFlagShift) & kWalkFlags; }
    void clear_flags(std::size_t slot) const noexcept { node_word(slot) &= ~(uint32_t(kWalkFlags) << kFlagShift); }

    /// Lets the node words be reused after the graph's last walk.
    ///
    /// The GFA writer stores segment ids there for link lookup, which saves building a separate index.
    /// aux_clear() zeroes every node word and marks the graph consumed, because abundances, masks and walker flags are gone.
    /// After that, save(), unitigs() and export_sorted() all refuse to run.
    void aux_clear(unsigned nthreads) {
        consumed_ = true;
        const std::size_t nlines = nlines_;
        const unsigned threads = nthreads ? nthreads : 1;

        tpoptoa::run_on_threads(threads, [&](unsigned thread_idx) {
            for (std::size_t l = nlines * thread_idx / threads; l < nlines * (thread_idx + 1) / threads; ++l)
                for (unsigned s = 0; s < kSlotsPerLine; ++s) lines_[l].nodes[s] = 0;
        });
    }
    void aux_set(std::size_t slot, uint32_t value) noexcept { __atomic_store_n(&node_word(slot), value, __ATOMIC_RELAXED); }
    uint32_t aux(std::size_t slot) const noexcept { return node(slot); }
    bool consumed() const noexcept { return consumed_; }

    /// Returns 0 for an absent k-mer.
    unsigned out_mask(uint64_t canon_kmer) const noexcept {
        const std::size_t slot = rank(canon_kmer);

        return slot == kNpos ? 0u : edge_mask(slot) & kOutEdgeMask;
    }
    unsigned in_mask(uint64_t canon_kmer) const noexcept {
        const std::size_t slot = rank(canon_kmer);

        return slot == kNpos ? 0u : edge_mask(slot) >> kOutEdgeBits;
    }

    /// Does nothing here, but DbgGraphEF and PartUnitigs share the interface.
    /// Callers use it when every observed neighbour is known to exist, so the stored edges are already final.
    void masks_from_raw() noexcept {}

    /// Trims every node's observed edges to neighbours that exist, in place.
    ///
    /// This is a software pipeline (Slotin, Algorithms for Modern Hardware).
    /// While block b is probed, the neighbour keys of block b + 1 are computed and their home lines prefetched, so each miss was issued a block earlier.
    /// Blocks stay small enough that the prefetches in flight fit a core's line-fill buffers.
    /// Each block reads its own observed edges before rewriting any mask, and probes only ever read key words.
    void resolve_masks() {
        const unsigned kmer_len = kmer_len_;
        const uint64_t kmask = kmer_mask(kmer_len);
        const unsigned first_base_shift = kBitsPerBase * (kmer_len - 1);
        auto canon = [&](uint64_t kmer) { const uint64_t rev_comp = reverse_complement(kmer, kmer_len); return kmer < rev_comp ? kmer : rev_comp; };
        constexpr std::size_t kBlock = 4, kSpan = 4096;     // slots per prefetch block / per OpenMP task
        constexpr unsigned kEdges = 8;                     // 4 out-edges, then 4 in-edges
        constexpr int kSpanSchedChunk = 4;                 // spans per OpenMP chunk
        constexpr std::size_t kPipelineStages = 2;          // the block being probed and the one being prefetched
        const std::size_t nspan = (cap_ + kSpan - 1) / kSpan;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, kSpanSchedChunk)
#endif
        for (std::size_t sp = 0; sp < nspan; ++sp) {
            const std::size_t span_lo = sp * kSpan, span_hi = std::min(cap_, span_lo + kSpan);
            // Neighbour keys of the current and next block, where ~0 marks a missing edge.
            uint64_t neighbours[kPipelineStages][kBlock][kEdges];
            auto prep = [&](std::size_t block_start, uint64_t (*block_neighbours)[kEdges]) {
                for (std::size_t slot = block_start; slot < std::min(span_hi, block_start + kBlock); ++slot) {
                    uint64_t* slot_neighbours = block_neighbours[slot - block_start];

                    // All ones can mark a missing edge, because no canonical k-mer is all ones.
                    for (unsigned b = 0; b < kEdges; ++b) slot_neighbours[b] = ~uint64_t(0);

                    if (!used(slot)) continue;

                    const uint64_t kmer = key(slot);
                    const unsigned raw = edge_mask(slot);

                    for (unsigned b = 0; b < kNumBases; ++b) {
                        if (raw & (1u << b)) {
                            slot_neighbours[b] = canon(((kmer << kBitsPerBase) | b) & kmask);
                            prefetch(slot_neighbours[b]);
                        }

                        if (raw & (1u << (b + kOutEdgeBits))) {
                            slot_neighbours[b + kOutEdgeBits] = canon(((uint64_t(b) << first_base_shift) | (kmer >> kBitsPerBase)) & kmask);
                            prefetch(slot_neighbours[b + kOutEdgeBits]);
                        }
                    }
                }
            };

            // With two buffers, the next block's neighbour keys are computed and prefetched while the current block's lookups resolve.
            prep(span_lo, neighbours[0]);

            for (std::size_t block_start = span_lo, cur = 0; block_start < span_hi; block_start += kBlock, cur ^= 1) {
                if (block_start + kBlock < span_hi) prep(block_start + kBlock, neighbours[cur ^ 1]);

                for (std::size_t slot = block_start; slot < std::min(span_hi, block_start + kBlock); ++slot) {
                    if (!used(slot)) continue;

                    unsigned mask = 0;

                    for (unsigned b = 0; b < kEdges; ++b)
                        if (neighbours[cur][slot - block_start][b] != ~uint64_t(0) && rank(neighbours[cur][slot - block_start][b]) != kNpos) mask |= 1u << b;

                    // An edge bit is cleared when its neighbour k-mer is missing from the table.
                    node_word(slot) = (node_word(slot) & ((uint32_t(1) << kMaskShift) - 1)) | (uint32_t(mask) << kMaskShift);
                }
            }
        }
    }

    /// Calls cb(seq, len, coverage) for every maximal unitig in head k-mer order, regardless of thread count.
    /// It is defined in dbg_unitigs.hpp.
    template <typename CB>
    void unitigs(CB&& callback, unsigned nthreads = 0) const;

    /// Batch form that calls sink(count, unitig_at), where unitig_at(i) yields {seq, len, cov} for the i-th unitig in output order.
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned nthreads = 0) const;

    /// Returns (canonical k-mer, abundance) pairs in ascending order, along with their edge masks.
    void export_sorted(std::vector<std::pair<uint64_t, uint32_t>>& pairs_out,
                       std::vector<uint8_t>& masks) const {
        if (consumed_) throw std::logic_error("DbgGraph: node words were reused (aux_clear)");

        std::vector<uint64_t> keys;

        keys.reserve(num_nodes_);

        for (std::size_t slot = 0; slot < cap_; ++slot)
            if (used(slot)) keys.push_back(key(slot));

        std::sort(keys.begin(), keys.end());
        pairs_out.resize(keys.size());
        masks.resize(keys.size());

        for (std::size_t i = 0; i < keys.size(); ++i) {
            const std::size_t slot = rank(keys[i]);

            pairs_out[i] = {keys[i], static_cast<uint32_t>(count(slot))};
            masks[i] = static_cast<uint8_t>(edge_mask(slot));
        }
    }

    /// Writes the little-endian kdbg archive of dbruijn -s:
    ///   "OKKDBG" | u8 k | u64 n | n x u64 sorted canonical k-mers | n x u32 abundances | n x u8 masks (in << 4 | out)
    /// It fails on an I/O error or on a consumed graph.
    bool save(const char* path) const;

    /// Builds from nodes whose masks are already trimmed to existing neighbours.
    /// @p pairs must be sorted and duplicate-free, and it is released afterwards.
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

    /// Node-centric build, where any (k-1)-overlap between two present k-mers counts as an edge.
    /// It suits k-mer sets that come without read windows, like DynSpss::compact().
    /// @p pairs must be sorted and duplicate-free.
    void build_core(std::vector<std::pair<uint64_t, uint32_t>>&& pairs) {
        init_table(pairs.size());

        for (const auto& pair : pairs)
            put(claim(pair.first), pair.second, kAllEdges);

        finish_insert(pairs.size());
        std::vector<std::pair<uint64_t, uint32_t>>().swap(pairs);
        resolve_masks();
    }

private:
    friend struct DbgUnitigExtractor;

    static constexpr uint64_t kHomeSeed = kGoldenRatio64;

    /// Maps the hash onto [0, nlines_) by fast range reduction (Lemire 2019).
    std::size_t home(uint64_t canon_kmer) const noexcept {
        return static_cast<std::size_t>(
            (static_cast<__uint128_t>(mix64(canon_kmer, kHomeSeed)) * nlines_) >> kWordBits);
    }
    uint64_t& key_word(std::size_t slot) const noexcept { return lines_[slot / kSlotsPerLine].keys[slot % kSlotsPerLine]; }
    uint32_t& node_word(std::size_t slot) const noexcept { return lines_[slot / kSlotsPerLine].nodes[slot % kSlotsPerLine]; }
    uint32_t node(std::size_t slot) const noexcept { return __atomic_load_n(&node_word(slot), __ATOMIC_RELAXED); }

    static std::size_t map_len(std::size_t nlines) noexcept {
        // Rounding up to whole huge pages lets advise_huge() back the entire table.
        return (nlines * sizeof(Line) + kHugePageBytes - 1) & ~(kHugePageBytes - 1);
    }
    static Line* alloc_zeroed(std::size_t nlines) {
        const std::size_t len = map_len(nlines);
        // Anonymous pages arrive zeroed and need no clearing pass.
        void* mem = ::mmap(nullptr, len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, kNoFd, 0);

        if (mem == MAP_FAILED) throw std::bad_alloc();

        advise_huge(mem, len);
        return static_cast<Line*>(mem);
    }
    struct LinePtr {
        Line* ptr = nullptr;
        std::size_t len = 0;

        LinePtr() = default;
        LinePtr(const LinePtr&) = delete;
        LinePtr& operator=(const LinePtr&) = delete;
        ~LinePtr() { if (ptr) ::munmap(ptr, len); }
        void reset(Line* new_ptr, std::size_t new_len) { if (ptr) ::munmap(ptr, len); ptr = new_ptr; len = new_len; }

        Line& operator[](std::size_t line_idx) const noexcept { return ptr[line_idx]; }
    };

    unsigned kmer_len_;
    std::size_t nlines_ = 0;
    std::size_t cap_ = 0;             // nlines_ * kSlotsPerLine slots
    bool consumed_ = false;           // set by aux_clear()
    std::size_t num_nodes_ = 0;       // claimed slots
    LinePtr lines_;
    std::vector<std::pair<std::size_t, uint32_t>> saturated_;   // (slot, abundance >= kCountMax), sorted
    std::mutex saturated_mutex_;
};

inline bool DbgGraph::save(const char* path) const {
    if (consumed_) return false;

    std::vector<std::pair<uint64_t, uint32_t>> pairs;
    std::vector<uint8_t> masks;

    export_sorted(pairs, masks);
    return atomic_write(path, [&](std::ofstream& out) {
        out.write(kKdbgMagic, kKdbgMagicLen);
        out.put(static_cast<char>(kmer_len_));
        const uint64_t count = pairs.size();

        out.write(reinterpret_cast<const char*>(&count), sizeof(count));

        for (const auto& pair : pairs) out.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));

        for (const auto& pair : pairs) out.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));

        out.write(reinterpret_cast<const char*>(masks.data()), static_cast<std::streamsize>(count));
    });
}

}

#include "dbg_unitigs.hpp"
