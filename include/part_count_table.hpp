#pragma once
/*
 * part_count_table.hpp - one thread's k-mer counting table for one minimizer partition at a time (the counting passes of src/kmer_count.cpp).
 *
 * Lines.
 * 64 bytes = 5 keys | 5 node words | u8 fill | u8 overflow | pad.
 * A node word is count << 8 | edge mask (in << 4 | out), the DbgGraph line layout (dbg_graph.hpp), so a probe reads a line of keys and their counts at once.
 *
 * Placement.
 * Within a level, a key goes to the first line with room in a window of kWindowLines lines starting at its home line, which is chosen by Fibonacci hashing (Knuth, TAOCP vol. 3, 6.4).
 * Lines fill in order and nothing is deleted, so a line with room ends a search, and the common insert or lookup reads one line.
 *
 * Levels.
 * When a key's whole window is full it sets its home line's overflow flag and moves one level down, and a lookup descends only past a full window whose home line is flagged.
 * This is the level structure of DeltaOA (delta_oa.hpp), after Farach-Colton, Krapivin and Kuszmaul (FOCS 2024).
 * Keys never move.
 *
 * Growth is the change made for counting.
 * DeltaOA's levels are fixed and shrink geometrically, whereas here a full table appends a level as large as all levels before it, so capacity doubles without rehashing or moving a key.
 * A partition dominated by one heavy minimizer (adapter read-through, poly-A) grows without a doubling table's rehash passes and copy.
 *
 * Sizing.
 * The caller sizes level 0 per partition from an estimate of its distinct k-mers (src/kmer_count.cpp uses a HyperLogLog pre-pass on large partitions, else the median distinct share of recent partitions), so nearly every partition stays in level 0.
 *
 * Reuse.
 * Storage is kept across partitions, where only used lines are cleared, and released past a size cap.
 * Cuttlefish 3 (Khan et al. 2025) also reuses worker-local maps across subgraphs and grows them from observed insertions.
 * The line layout, bounded line-linear placement, levels and append-only growth are this table's own.
 *
 * Counts.
 * The node word holds 24 bits, and a key whose count reaches kCountSat keeps the excess in an exact side map, and totals saturate at 2^32 - 1 as the counter's always have.
 *
 * The table is single-threaded, since one thread owns a partition, so DeltaOA's seqlocks are not needed.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

namespace tpoptoa {

class PartCountTable {
public:
    static constexpr unsigned kSlots = 5;              // keys per 64-byte line
    static constexpr std::size_t kWindowLines = 8;          // lines probed per level
    static constexpr uint32_t kNotFound = ~uint32_t(0);
    static constexpr unsigned kCountShift = 8;        // node word = count << 8 | edge mask
    static constexpr uint32_t kEdgeMask = 0xFFu;
    static constexpr uint32_t kCountSat = (uint32_t(1) << 24) - 1;
    static constexpr uint64_t kFibonacciMultiplier = 0x9E3779B97F4A7C15ULL;   // 2^64 / golden ratio
    static constexpr unsigned kLevelRotation = 13;    // hash bits rotated per level

    /// begin() - start a partition expected to hold about @est_distinct keys, with level 0 sized for a load factor of about 3/4.
    void begin(std::size_t est_distinct) {
        const std::size_t want = (est_distinct * 4 + 2 * kSlots) / (3 * kSlots) + 1;
        std::size_t lines = kWindowLines;
        while (lines < want) lines <<= 1;
        nlevels_ = 0;
        add_level(lines);
        size_ = 0;
        deep_keys_ = 0;
        excess_counts_.clear();
    }

    /// hash() - Fibonacci hash of @key, whose top bits a level takes.
    ///
    /// The caller hashes keys ahead and prefetches their home lines (group prefetching, Slotin, Algorithms for Modern Hardware), so the cache misses of a large partition overlap.
    static constexpr uint64_t hash(uint64_t key) noexcept { return key * kFibonacciMultiplier; }
    /// prefetch() - prefetch the level-0 home line of hash @h for writing.
    void prefetch(uint64_t h) const noexcept {
        const Level& L = levels_[0];
        __builtin_prefetch(&L.lines[home(h, 0, L.log2_lines)], 1);
    }

    /// add() - count @n occurrences of @key, OR-ing in edge bits @bits.
    void add(uint64_t key, uint8_t bits, uint32_t n = 1) { add_hashed(key, hash(key), bits, n); }
    /// add_hashed() - add() with the key's hash() already computed.
    void add_hashed(uint64_t key, uint64_t h, uint8_t bits, uint32_t n = 1) {
        for (unsigned i = 0;; ++i) {
            if (i == nlevels_) add_level(total_lines_);          // append, which doubles capacity
            Level& L = levels_[i];
            const std::size_t h0 = home(h, i, L.log2_lines), line_mask = L.n - 1;
            for (std::size_t w = 0; w < kWindowLines; ++w) {
                Line& l = L.lines[(h0 + w) & line_mask];
                const int s = slot_of(l, key);
                if (s >= 0) { bump(l, s, key, bits, n); return; }
                if (l.fill < kSlots) {                         // room here means the key is nowhere further
                    l.key[l.fill] = key;
                    l.node[l.fill] = bits;
                    bump(l, l.fill, key, 0, n);
                    ++l.fill;
                    ++size_;
                    if (i) ++deep_keys_;
                    return;
                }
            }
            L.lines[h0].overflow = 1;                              // window full, so go one level down
        }
    }

    /// find() - the count of @key including the side map, or kNotFound.
    uint32_t find(uint64_t key) const {
        const uint64_t h = hash(key);
        for (unsigned i = 0; i < nlevels_; ++i) {
            const Level& L = levels_[i];
            const std::size_t h0 = home(h, i, L.log2_lines), line_mask = L.n - 1;
            for (std::size_t w = 0; w < kWindowLines; ++w) {
                const Line& l = L.lines[(h0 + w) & line_mask];
                const int s = slot_of(l, key);
                if (s >= 0) return count(key, l.node[s]);
                if (l.fill < kSlots) return kNotFound;
            }
            if (!L.lines[h0].overflow) return kNotFound;
        }
        return kNotFound;
    }

    /// for_each() - f(key, count, edge mask) for every key, in storage order.
    template <typename F>
    void for_each(F&& f) const {
        for (unsigned i = 0; i < nlevels_; ++i) {
            const Level& L = levels_[i];
            for (std::size_t j = 0; j < L.n; ++j) {
                const Line& l = L.lines[j];
                for (unsigned s = 0; s < l.fill; ++s)
                    f(l.key[s], count(l.key[s], l.node[s]), static_cast<uint8_t>(l.node[s] & kEdgeMask));
            }
        }
    }

    std::size_t size() const noexcept { return size_; }
    /// deep() and levels() - for profiling, the keys placed below level 0 and the levels in use.
    std::size_t deep() const noexcept { return deep_keys_; }
    unsigned levels() const noexcept { return nlevels_; }
    /// bytes() - line storage held by all levels, used or kept for reuse.
    std::size_t bytes() const noexcept {
        std::size_t b = 0;
        for (const auto& L : levels_) b += L.cap * sizeof(Line);
        return b;
    }

    /// finish() - clear the partition's lines and keep the storage, or free it all if it exceeds @cap_bytes.
    void finish(std::size_t cap_bytes) {
        if (bytes() > cap_bytes) {
            levels_.clear();
        } else {
            for (unsigned i = 0; i < nlevels_; ++i)
                std::memset(static_cast<void*>(levels_[i].lines), 0, levels_[i].n * sizeof(Line));
        }
        nlevels_ = 0;
        total_lines_ = 0;
        size_ = 0;
        deep_keys_ = 0;
        excess_counts_.clear();
    }

private:
    struct alignas(64) Line {
        uint64_t key[kSlots];
        uint32_t node[kSlots];
        uint8_t fill, overflow;         // filled slots and the overflow flag of a home line
        uint8_t pad[2];
    };
    static_assert(sizeof(Line) == 64, "one cache line");

    struct Level {
        std::unique_ptr<Line[]> store;
        Line* lines = nullptr;
        std::size_t n = 0, cap = 0;       // lines in use / allocated
        unsigned log2_lines = 0;                  // log2(n)
    };

    /// add_level() - append level nlevels_ with @lines lines (a power of two >= kWindowLines), reusing that level's storage from an earlier partition.
    void add_level(std::size_t lines) {
        if (levels_.size() <= nlevels_) levels_.emplace_back();
        Level& L = levels_[nlevels_];
        if (L.cap < lines) {
            L.store.reset(new Line[lines]);
            L.cap = lines;
            std::memset(static_cast<void*>(L.store.get()), 0, lines * sizeof(Line));
        }
        L.lines = L.store.get();
        L.n = lines;
        L.log2_lines = 0;
        while ((std::size_t(1) << L.log2_lines) < lines) ++L.log2_lines;
        total_lines_ += lines;
        ++nlevels_;
    }

    /// home() - home line of hash @h on level @i: the top @log2_lines bits of the hash rotated left by kLevelRotation * i bits, so a key's home lines on different levels are independent.
    static constexpr std::size_t home(uint64_t h, unsigned i, unsigned log2_lines) noexcept {
        if (i) h = (h << (i * kLevelRotation & 63)) | (h >> ((64 - (i * kLevelRotation & 63)) & 63));
        return static_cast<std::size_t>(h >> (64 - log2_lines));
    }
    /// slot_of() - index of @key among the line's filled slots, or -1.
    /// All kSlots keys are compared without branches and the result is masked by the fill.
    static int slot_of(const Line& l, uint64_t key) noexcept {
        unsigned m = 0;
        for (unsigned s = 0; s < kSlots; ++s) m |= unsigned(l.key[s] == key) << s;
        m &= (1u << l.fill) - 1u;
        return m ? __builtin_ctz(m) : -1;
    }
    /// bump() - @n more occurrences of the key in slot @s, OR-ing in edge bits @bits.
    /// Past kCountSat the excess goes to the side map, which saturates at a total of 2^32 - 1.
    void bump(Line& l, int s, uint64_t key, uint8_t bits, uint32_t n) {
        uint32_t& node = l.node[s];
        node |= bits;
        const uint32_t room = kCountSat - (node >> kCountShift);
        if (n <= room) { node += n << kCountShift; return; }
        node += room << kCountShift;
        uint64_t& excess = excess_counts_[key];
        excess = std::min<uint64_t>(excess + (n - room), uint64_t(UINT32_MAX) - kCountSat);
    }
    /// count() - total count of @key whose node word is @node.
    uint32_t count(uint64_t key, uint32_t node) const {
        const uint32_t c = node >> kCountShift;
        if (c < kCountSat) return c;
        auto it = excess_counts_.find(key);
        return it == excess_counts_.end() ? c : static_cast<uint32_t>(std::min<uint64_t>(UINT32_MAX, kCountSat + it->second));
    }

    std::vector<Level> levels_;
    unsigned nlevels_ = 0;
    std::size_t total_lines_ = 0;
    std::size_t size_ = 0, deep_keys_ = 0;
    std::unordered_map<uint64_t, uint64_t> excess_counts_;
};

}
