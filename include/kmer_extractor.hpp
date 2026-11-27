#pragma once

/**
 * @file
 * @brief Rolling extraction of canonical k-mers.
 *
 * Both strands are tracked in a single pass at O(1) per base.
 * The forward word shifts left and takes the new base at the bottom.
 * Meanwhile, the reverse-complement word shifts right and takes the complemented base at bit 2(k-1).
 * Once k bases are in, the smaller of the two words is the canonical k-mer.
 */

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include "kmer.hpp"

namespace tpoptoa {

/// Calls @p callback with each canonical k-mer of @p seq in order.
/// N, IUPAC codes and other non-bases reset both words, so no k-mer crosses them.
template <typename Callback>
void extract_kmers(const char* seq, std::size_t len,
                   std::size_t kmer_len, Callback&& callback)
{
    if (len < kmer_len) return;

    KmerWord mask = kmer_mask(kmer_len);
    KmerWord fwd = 0;
    KmerWord rev = 0;
    std::size_t valid = 0;  // bases since the last restart
    unsigned rc_shift = static_cast<unsigned>(kBitsPerBase * (kmer_len - 1));

    for (std::size_t i = 0; i < len; ++i) {
        uint8_t fwd_base = base_encode(seq[i]);

        // An N restarts both rolling k-mers.
        if (fwd_base == kInvalidBase) {
            fwd = 0;
            rev = 0;
            valid = 0;
            continue;
        }

        fwd = ((fwd << kBitsPerBase) | fwd_base) & mask;

        uint8_t rev_base = base_complement(fwd_base);

        // The reverse complement grows from the high end as the forward k-mer grows from the low end.
        rev = (rev >> kBitsPerBase) | (static_cast<KmerWord>(rev_base) << rc_shift);
        rev &= mask;

        ++valid;

        if (valid >= kmer_len) {
            KmerWord can = (fwd < rev) ? fwd : rev;

            callback(can);
        }
    }
}

template <typename Callback>
void extract_kmers(const std::string& seq, std::size_t kmer_len, Callback&& callback) {
    extract_kmers(seq.data(), seq.size(), kmer_len, std::forward<Callback>(callback));
}

}
