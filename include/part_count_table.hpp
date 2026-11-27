#pragma once
/**
 * @file
 * @brief The k-mer counting table one thread uses for one minimizer partition at a time.
 * The counting passes in src/kmer_count.cpp are its only users.
 *
 * @par Lines
 * 64 bytes = 5 keys | 5 node words | u8 fill | u8 overflow | pad.
 * A node word is count << 8 | edge mask (in << 4 | out), the same layout DbgGraph uses (dbg_graph.hpp).
 * As a result, one probe reads a line of keys together with their counts.
 *
 * @par Placement
 * A key lands in the first line with room inside a window of kWindowLines lines that starts at its home line.
 * Fibonacci hashing picks the home line (Knuth, TAOCP vol. 3, 6.4).
 * Because lines fill in order and nothing is ever deleted, a line with room ends the search, so most operations read a single line.
 *
 * @par Levels
 * When a key's whole window is full, it flags its home line as overflowed and moves one level down.
 * Lookups therefore descend only past a full window whose home line has that flag set.
 * The scheme comes from DeltaOA (delta_oa.hpp), which in turn follows Farach-Colton, Krapivin and Kuszmaul (FOCS 2024), and keys never move.
 *
 * @par Growth
 * DeltaOA's levels are fixed and shrink geometrically.
 * Here, instead, a full table appends a level as large as everything before it, so capacity doubles without rehashing.
 * That helps partitions dominated by one heavy minimizer (adapter read-through, poly-A), which would otherwise pay for repeated rehash-and-copy passes.
 *
 * @par Sizing
 * The caller sizes level 0 from an estimate of the partition's distinct k-mers, so nearly every partition stays in level 0.
 * For large partitions src/kmer_count.cpp estimates with a HyperLogLog pre-pass, and otherwise uses the median distinct share of recent partitions.
 *
 * @par Reuse
 * Storage is reused across partitions, with only the used lines cleared, until it passes a size cap.
 * Cuttlefish 3 (Khan et al. 2025) similarly reuses worker-local maps and grows them from observed insertions.
 * The line layout, bounded placement, levels and append-only growth are specific to this table, though.
 *
 * @par Counts
 * A node word holds 24 bits of count, so anything past kCountSat goes to an exact side map.
 * Totals still saturate at 2^32 - 1, matching the counter's long-standing behaviour.
 *
 * One thread owns a partition, which is why none of DeltaOA's seqlocks are needed.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

#include "bits.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

class PartCountTable {
public:
    static constexpr std::size_t kCacheLineBytes = 64;   // x86-64, the unit one probe reads
    static constexpr unsigned kSlots = 5;                // keys per line
    /// Explicit padding after the two flag bytes, so a line's fields end on a 4-byte boundary.
    static constexpr std::size_t kLinePadBytes = 2;
    static constexpr std::size_t kWindowLines = 8;       // lines probed per level
    static constexpr uint32_t kNotFound = ~uint32_t(0);
    static constexpr unsigned kCountShift = 8;        // node word = count << 8 | edge mask
    static constexpr uint32_t kEdgeMask = 0xFFu;
    static constexpr uint32_t kCountSat = (uint32_t(1) << 24) - 1;
    static constexpr unsigned kLevelRotation = 13;    // hash bits rotated per level
    /// Level 0 is sized for a load factor of about kLoadNum / kLoadDen.
    static constexpr std::size_t kLoadNum = 3, kLoadDen = 4;
    static constexpr int kNoSlot = -1;                // slot_of() for a key the line doesn't hold

    /// Starts a partition expected to hold about @p est_distinct keys.
    void begin(std::size_t est_distinct) {
        constexpr std::size_t kRoundUpSlots = 2 * kSlots;   // pads the division so a partly filled line still counts
        const std::size_t want = (est_distinct * kLoadDen + kRoundUpSlots) / (kLoadNum * kSlots) + 1;
        std::size_t lines = kWindowLines;

        while (lines < want) lines <<= 1;

        nlevels_ = 0;
        add_level(lines);
        size_ = 0;
        deep_keys_ = 0;
        excess_counts_.clear();
    }

    /// Fibonacci hash, of which each level takes the top bits.
    /// Callers hash ahead and prefetch the home lines (group prefetching, Slotin), so a large partition's cache misses overlap.
    static constexpr uint64_t hash(uint64_t key) noexcept { return key * kGoldenRatio64; }

    void prefetch(uint64_t hash_val) const noexcept {
        const Level& level = levels_[0];

        prefetch_write(&level.lines[home(hash_val, 0, level.log2_lines)]);
    }

    /// Adds @p occurrences to @p key's count and ORs in edge bits @p bits.
    void add(uint64_t key, uint8_t bits, uint32_t occurrences = 1) { add_hashed(key, hash(key), bits, occurrences); }

    /// Same as add(), for callers that already computed hash().
    void add_hashed(uint64_t key, uint64_t hash_val, uint8_t bits, uint32_t occurrences = 1) {
        for (unsigned i = 0;; ++i) {
            if (i == nlevels_) add_level(total_lines_);          // appending doubles capacity

            Level& level = levels_[i];
            const std::size_t home_line = home(hash_val, i, level.log2_lines), line_mask = level.used_lines - 1;

            for (std::size_t w = 0; w < kWindowLines; ++w) {
                Line& line = level.lines[(home_line + w) & line_mask];
                const int slot = slot_of(line, key);

                if (slot >= 0) { bump(line, slot, key, bits, occurrences); return; }

                if (line.fill < kSlots) {                         // lines fill in order, so the key can't be further on
                    line.key[line.fill] = key;
                    line.node[line.fill] = bits;
                    bump(line, line.fill, key, 0, occurrences);
                    ++line.fill;
                    ++size_;

                    if (i) ++deep_keys_;

                    return;
                }
            }

            level.lines[home_line].overflow = 1;                              // window full, so the key drops a level
        }
    }

    /// Includes any excess held in the side map.
    uint32_t find(uint64_t key) const {
        const uint64_t hash_val = hash(key);

        for (unsigned i = 0; i < nlevels_; ++i) {
            const Level& level = levels_[i];
            const std::size_t home_line = home(hash_val, i, level.log2_lines), line_mask = level.used_lines - 1;

            for (std::size_t w = 0; w < kWindowLines; ++w) {
                const Line& line = level.lines[(home_line + w) & line_mask];
                const int slot = slot_of(line, key);

                if (slot >= 0) return count(key, line.node[slot]);

                if (line.fill < kSlots) return kNotFound;
            }

            if (!level.lines[home_line].overflow) return kNotFound;
        }

        return kNotFound;
    }

    /// Calls f(key, count, edge mask) in storage order.
    template <typename F>
    void for_each(F&& func) const {
        for (unsigned i = 0; i < nlevels_; ++i) {
            const Level& level = levels_[i];

            for (std::size_t j = 0; j < level.used_lines; ++j) {
                const Line& line = level.lines[j];

                for (unsigned s = 0; s < line.fill; ++s)
                    func(line.key[s], count(line.key[s], line.node[s]), static_cast<uint8_t>(line.node[s] & kEdgeMask));
            }
        }
    }

    std::size_t size() const noexcept { return size_; }
    /// Reports how many keys sit below level 0 and how many levels exist, for profiling.
    std::size_t deep() const noexcept { return deep_keys_; }
    unsigned levels() const noexcept { return nlevels_; }

    /// Counts storage kept for reuse as well as lines in use.
    std::size_t bytes() const noexcept {
        std::size_t total = 0;

        for (const auto& level : levels_) total += level.cap * sizeof(Line);

        return total;
    }

    /// Clears the partition's lines but keeps the storage, unless it exceeds @p cap_bytes.
    void finish(std::size_t cap_bytes) {
        if (bytes() > cap_bytes) {
            levels_.clear();
        } else {
            for (unsigned i = 0; i < nlevels_; ++i)
                std::memset(static_cast<void*>(levels_[i].lines), 0, levels_[i].used_lines * sizeof(Line));
        }

        nlevels_ = 0;
        total_lines_ = 0;
        size_ = 0;
        deep_keys_ = 0;
        excess_counts_.clear();
    }

private:
    struct alignas(kCacheLineBytes) Line {
        uint64_t key[kSlots];
        uint32_t node[kSlots];
        uint8_t fill, overflow;         // overflow is only meaningful on a home line
        uint8_t pad[kLinePadBytes];
    };
    static_assert(sizeof(Line) == kCacheLineBytes, "one cache line");

    struct Level {
        std::unique_ptr<Line[]> store;
        Line* lines = nullptr;
        std::size_t used_lines = 0, cap = 0;       // lines in use / allocated
        unsigned log2_lines = 0;
    };

    /// Appends a level of @p lines lines, a power of two no smaller than kWindowLines.
    /// Storage left over from an earlier partition at that depth gets reused.
    void add_level(std::size_t lines) {
        if (levels_.size() <= nlevels_) levels_.emplace_back();

        Level& level = levels_[nlevels_];

        if (level.cap < lines) {
            level.store.reset(new Line[lines]);
            level.cap = lines;
            std::memset(static_cast<void*>(level.store.get()), 0, lines * sizeof(Line));
        }

        level.lines = level.store.get();
        level.used_lines = lines;
        level.log2_lines = 0;

        while ((std::size_t(1) << level.log2_lines) < lines) ++level.log2_lines;

        total_lines_ += lines;
        ++nlevels_;
    }

    /// Takes the top @p log2_lines bits after rotating the hash by kLevelRotation per level.
    /// The rotation keeps a key's home lines on different levels independent of each other.
    static constexpr std::size_t home(uint64_t hash_val, unsigned level_idx, unsigned log2_lines) noexcept {
        hash_val = rotate_left(hash_val, level_idx * kLevelRotation);

        return static_cast<std::size_t>(hash_val >> (kWordBits - log2_lines));
    }
    /// Compares all kSlots keys without branches, then masks the result by the fill count.
    static int slot_of(const Line& line, uint64_t key) noexcept {
        unsigned match_mask = 0;

        for (unsigned s = 0; s < kSlots; ++s) match_mask |= unsigned(line.key[s] == key) << s;

        match_mask &= (1u << line.fill) - 1u;
        return match_mask ? __builtin_ctz(match_mask) : kNoSlot;
    }
    /// Past kCountSat the excess moves to the side map, which itself saturates at a total of 2^32 - 1.
    void bump(Line& line, int slot, uint64_t key, uint8_t bits, uint32_t occurrences) {
        uint32_t& node = line.node[slot];

        node |= bits;
        const uint32_t room = kCountSat - (node >> kCountShift);

        // Most counts fit in the 24-bit field, and excess_counts_ is touched only by k-mers past kCountSat.
        if (occurrences <= room) { node += occurrences << kCountShift; return; }

        node += room << kCountShift;
        uint64_t& excess = excess_counts_[key];

        excess = std::min<uint64_t>(excess + (occurrences - room), uint64_t(UINT32_MAX) - kCountSat);
    }
    uint32_t count(uint64_t key, uint32_t node) const {
        const uint32_t stored = node >> kCountShift;

        if (stored < kCountSat) return stored;

        // A saturated field defers to the side map for the rest of the count.
        auto found = excess_counts_.find(key);

        return found == excess_counts_.end() ? stored : static_cast<uint32_t>(std::min<uint64_t>(UINT32_MAX, kCountSat + found->second));
    }

    std::vector<Level> levels_;
    unsigned nlevels_ = 0;
    std::size_t total_lines_ = 0;
    std::size_t size_ = 0, deep_keys_ = 0;
    std::unordered_map<uint64_t, uint64_t> excess_counts_;
};

}
