#pragma once

/*
 * minimizer.hpp - window minimizers of a nucleotide sequence.
 *
 * A window is w consecutive canonical k-mers, and its minimizer is the k-mer that is smallest under a random order.
 * Keeping only minimizers samples about 2/(w+1) of the k-mer positions while guaranteeing that two sequences sharing w + k - 1 consecutive bases share a minimizer.
 *
 * Sources:
 *   Roberts et al. 2004 define window minimizers.
 *   Li 2018 (minimap2) orders k-mers by an invertible integer hash rather
 *   than lexicographically, which avoids the density loss and poly-A bias
 *   of the lexicographic order.  minimizer_order_hash() is that hash.
 *
 * minimap2 keeps the window minimum in a ring buffer and rescans it when the minimum leaves the window.
 * extract_minimizers() modifies this by doing the same with a branch-free rescan (the conditional moves follow the branch-free sliding minimum of Slotin, Algorithms for Modern Hardware), so random hashes do not cause data-dependent branch mispredictions.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <functional>
#include <limits>

#include "kmer.hpp"
#include "kmer_extractor.hpp"

namespace tpoptoa {

/// Largest window extract_minimizers() accepts.
/// It is also the ring-buffer size, so it must stay a power of two.
inline constexpr std::size_t kMaxMinimizerWindow = 256;

struct MinimizerSeed {
    KmerWord kmer;   ///< canonical k-mer
    uint32_t pos;    ///< 0-based start offset of the k-mer in the sequence
    bool     rev;    ///< true when the sequence carries rc(kmer) at pos
};

/// minimizer_order_hash() - invertible 64-bit mix restricted to @mask, used to order k-mers for minimizer selection.
///
/// Thomas Wang's integer hash, as minimap2 uses it (hash64() in sketch.c).
/// Invertible, so distinct k-mers never tie.
constexpr inline uint64_t minimizer_order_hash(uint64_t key, uint64_t mask) noexcept {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

/// extract_minimizers() - call @cb with the window minimizers of @seq.
///
/// Windows are @w consecutive canonical @k-mers of one segment free of non-ACGT characters.
/// Such a character ends the segment and resets the window, so no window spans an ambiguous base.
/// Within a window the k-mer with the smallest minimizer_order_hash() wins, the leftmost one on ties.
/// Each selected occurrence is reported once, with its sequence offset and strand.
/// Segments with fewer than @w k-mers yield no minimizer, and
/// @w must be at most kMaxMinimizerWindow.
///
/// The current minimum is replaced only when a strictly smaller hash enters, which keeps the leftmost tie.
/// The window is rescanned only when the minimum leaves it, which happens for about 2/(w+1) of the k-mers.
template <typename Callback>
void extract_minimizers(const char* seq, std::size_t len,
                        std::size_t k, std::size_t w,
                        Callback&& cb)
{
    if (k == 0 || w == 0 || w > kMaxMinimizerWindow || len < k + w - 1) return;
    const KmerWord mask = kmer_mask(k);
    const uint64_t hmask = mask;
    const unsigned rc_shift = static_cast<unsigned>(2 * (k - 1));
    constexpr std::size_t kRingMask = kMaxMinimizerWindow - 1;

    struct WindowEntry { uint64_t h; KmerWord kmer; uint32_t pos; bool rev; };
    WindowEntry ring[kMaxMinimizerWindow];    // the last w k-mers, by k-mer index mod ring size
    KmerWord fwd = 0, rev = 0;
    std::size_t valid = 0;                     // valid bases in current segment
    std::size_t kmer_index = 0;                // k-mers seen in current segment
    uint32_t last_pos = UINT32_MAX;            // position of the last reported minimizer
    std::size_t min_index = 0;                 // k-mer index of the window minimum
    uint64_t min_hash = ~uint64_t(0);

    for (std::size_t i = 0; i < len; ++i) {
        const uint8_t b = base_encode(seq[i]);
        if (b == kInvalidBase) { valid = 0; kmer_index = 0; continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = ((rev >> 2) | (static_cast<KmerWord>(b ^ 3u) << rc_shift)) & mask;
        if (++valid < k) continue;
        const bool is_rev = rev < fwd;
        const KmerWord can = is_rev ? rev : fwd;
        const uint64_t h = minimizer_order_hash(can, hmask);
        const std::size_t cur = kmer_index++;
        ring[cur & kRingMask] = WindowEntry{h, can, static_cast<uint32_t>(i + 1 - k), is_rev};
        if (cur == 0) { min_hash = h; min_index = 0; }
        else if (h < min_hash) { min_hash = h; min_index = cur; }
        if (cur + 1 < w) continue;
        const std::size_t first = cur + 1 - w;
        if (min_index < first) {               // minimum left the window, so rescan
            min_hash = ~uint64_t(0); min_index = first;
            for (std::size_t t = first; t <= cur; ++t) {
                const uint64_t ht = ring[t & kRingMask].h;
                const bool lt = ht < min_hash;
                min_hash = lt ? ht : min_hash;
                min_index = lt ? t : min_index;
            }
        }
        const WindowEntry& m = ring[min_index & kRingMask];
        if (m.pos != last_pos) {
            cb(MinimizerSeed{m.kmer, m.pos, m.rev});
            last_pos = m.pos;
        }
    }
}

template <typename Callback>
void extract_minimizers(const std::string& seq,
                        std::size_t k, std::size_t w,
                        Callback&& cb)
{
    extract_minimizers(seq.data(), seq.size(), k, w,
                       std::forward<Callback>(cb));
}

}
