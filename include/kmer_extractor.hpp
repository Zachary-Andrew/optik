#pragma once

/*
 * kmer_extractor.hpp - rolling extraction of canonical k-mers.
 *
 * One pass over the sequence keeps the forward k-mer and its reverse complement up to date in O(1) per base.
 * The forward word shifts left and takes the new base at the bottom, and the reverse-complement word shifts right and takes the complement of the new base at the top (bit 2(k-1)).
 * After k bases the two words are exact reverse complements, and the canonical k-mer is the smaller of the two.
 */

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include "kmer.hpp"

namespace tpoptoa {

/// extract_kmers() - call @cb with every canonical k-mer of @seq, in order.
///
/// A character that is not a base, such as N or an IUPAC code, restarts both rolling words, so no reported k-mer spans it.
template <typename Callback>
void extract_kmers(const char* seq, std::size_t len,
                   std::size_t k, Callback&& cb)
{
    if (len < k) return;

    KmerWord mask = kmer_mask(k);
    KmerWord fwd = 0;
    KmerWord rev = 0;
    std::size_t valid = 0;  // bases since the last restart
    unsigned rc_shift = static_cast<unsigned>(2 * (k - 1));

    for (std::size_t i = 0; i < len; ++i) {
        uint8_t fb = base_encode(seq[i]);

        if (fb == kInvalidBase) {
            fwd = 0;
            rev = 0;
            valid = 0;
            continue;
        }

        fwd = ((fwd << 2) | fb) & mask;

        uint8_t rb = base_complement(fb);
        rev = (rev >> 2) | (static_cast<KmerWord>(rb) << rc_shift);
        rev &= mask;

        ++valid;
        if (valid >= k) {
            KmerWord can = (fwd < rev) ? fwd : rev;
            cb(can);
        }
    }
}

template <typename Callback>
void extract_kmers(const std::string& seq, std::size_t k, Callback&& cb) {
    extract_kmers(seq.data(), seq.size(), k, std::forward<Callback>(cb));
}

}
