#pragma once

/**
 * @file
 * @brief Window minimizers of a nucleotide sequence.
 *
 * Out of every w consecutive canonical k-mers, we keep the one with the smallest hash (Roberts et al. 2004).
 * Roughly 2/(w+1) of positions survive, yet any two sequences sharing w + k - 1 bases still share a minimizer.
 *
 * Ordering follows minimap2 (Li 2018), using an invertible integer hash rather than lexicographic order.
 * That choice avoids poly-A bias and the density loss of lexicographic minimizers.
 * Like minimap2, we keep the window in a ring buffer and rescan when the minimum drops out.
 * The rescan, however, uses conditional moves (Slotin's branch-free sliding minimum), so random hashes cause no mispredictions.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <functional>
#include <limits>

#include "kmer.hpp"
#include "kmer_extractor.hpp"

namespace tpoptoa {

/// Doubles as the ring-buffer size, so keep it a power of two.
inline constexpr std::size_t kMaxMinimizerWindow = 256;

struct MinimizerSeed {
    KmerWord kmer;   ///< canonical k-mer
    uint32_t pos;    ///< 0-based start offset of the k-mer in the sequence
    bool     rev;    ///< true when the sequence contains rc(kmer) at pos
};

/// Thomas Wang's 64-bit hash shifts, in application order (minimap2's hash64()).
inline constexpr unsigned kWangShift1 = 21, kWangShift2 = 24, kWangShift3a = 3, kWangShift3b = 8, kWangShift4 = 14,
                          kWangShift5a = 2, kWangShift5b = 4, kWangShift6 = 28, kWangShift7 = 31;

/// Being invertible, distinct k-mers can never tie.
constexpr inline uint64_t minimizer_order_hash(uint64_t key, uint64_t mask) noexcept {
    key = (~key + (key << kWangShift1)) & mask;
    key = key ^ (key >> kWangShift2);
    key = ((key + (key << kWangShift3a)) + (key << kWangShift3b)) & mask;
    key = key ^ (key >> kWangShift4);
    key = ((key + (key << kWangShift5a)) + (key << kWangShift5b)) & mask;
    key = key ^ (key >> kWangShift6);
    key = (key + (key << kWangShift7)) & mask;
    return key;
}

/// Reports each selected occurrence once, with offset and strand.
/// Non-ACGT characters reset the window, so segments shorter than @p window k-mers produce nothing.
/// Ties go to the leftmost k-mer, since only a strictly smaller hash replaces the minimum.
template <typename Callback>
void extract_minimizers(const char* seq, std::size_t len,
                        std::size_t kmer_len, std::size_t window,
                        Callback&& callback)
{
    if (kmer_len == 0 || window == 0 || window > kMaxMinimizerWindow || len < kmer_len + window - 1) return;

    const KmerWord mask = kmer_mask(kmer_len);
    const uint64_t hmask = mask;
    const unsigned rc_shift = static_cast<unsigned>(kBitsPerBase * (kmer_len - 1));
    constexpr std::size_t kRingMask = kMaxMinimizerWindow - 1;

    struct WindowEntry { uint64_t hash; KmerWord kmer; uint32_t pos; bool rev; };
    WindowEntry ring[kMaxMinimizerWindow];    // last w k-mers, slot = k-mer index mod ring size
    KmerWord fwd = 0, rev = 0;
    std::size_t valid = 0;                     // bases since the last reset
    std::size_t kmer_index = 0;
    uint32_t last_pos = UINT32_MAX;            // avoids reporting the same minimizer twice
    std::size_t min_index = 0;
    uint64_t min_hash = ~uint64_t(0);

    for (std::size_t i = 0; i < len; ++i) {
        const uint8_t base = base_encode(seq[i]);

        if (base == kInvalidBase) { valid = 0; kmer_index = 0; continue; }

        fwd = ((fwd << kBitsPerBase) | base) & mask;
        rev = ((rev >> kBitsPerBase) | (static_cast<KmerWord>(base ^ kComplementMask) << rc_shift)) & mask;

        if (++valid < kmer_len) continue;

        const bool is_rev = rev < fwd;
        const KmerWord can = is_rev ? rev : fwd;
        const uint64_t hash = minimizer_order_hash(can, hmask);
        const std::size_t cur = kmer_index++;

        ring[cur & kRingMask] = WindowEntry{hash, can, static_cast<uint32_t>(i + 1 - kmer_len), is_rev};

        if (cur == 0) { min_hash = hash; min_index = 0; }
        else if (hash < min_hash) { min_hash = hash; min_index = cur; }

        if (cur + 1 < window) continue;

        const std::size_t first = cur + 1 - window;

        if (min_index < first) {               // minimum slid out, which happens for ~2/(w+1) of k-mers
            min_hash = ~uint64_t(0); min_index = first;

            for (std::size_t t = first; t <= cur; ++t) {
                const uint64_t entry_hash = ring[t & kRingMask].hash;
                const bool is_less = entry_hash < min_hash;

                min_hash = is_less ? entry_hash : min_hash;
                min_index = is_less ? t : min_index;
            }
        }

        const WindowEntry& min_entry = ring[min_index & kRingMask];

        // Comparing pos against last_pos reports a minimizer once even when many consecutive windows share it.
        if (min_entry.pos != last_pos) {
            callback(MinimizerSeed{min_entry.kmer, min_entry.pos, min_entry.rev});
            last_pos = min_entry.pos;
        }
    }
}

template <typename Callback>
void extract_minimizers(const std::string& seq,
                        std::size_t kmer_len, std::size_t window,
                        Callback&& callback)
{
    extract_minimizers(seq.data(), seq.size(), kmer_len, window,
                       std::forward<Callback>(callback));
}

}
