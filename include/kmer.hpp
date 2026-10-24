#pragma once

/*
 * kmer.hpp - 2-bit packed k-mer encoding.
 *
 * Each base is stored in 2 bits (A=0, C=1, G=2, T=3), so a k-mer of up to kMaxK = 32 bases fits one 64-bit word.
 * The first base of a k-mer occupies the highest 2 bits of the word, since extraction shifts left by 2 per base, so after k bases the first one has reached bits [2k-2, 2k).
 *
 * With this code the complement of a base is its XOR with 3, and the reverse complement of a packed k-mer is a word-level bit reversal (byte swap, then nibble and pair swaps) applied to 2-bit fields.
 */

#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace tpoptoa {

static constexpr std::size_t kMaxK = 32;

/// Code returned by base_encode() for any character other than A, C, G, T and U, such as N, IUPAC ambiguity codes and gaps.
/// Callers treat it as a break in the sequence, which no k-mer spans.
static constexpr uint8_t kInvalidBase = 0xFF;

using KmerWord = uint64_t;

/// base_encode() - 2-bit code of a nucleotide, case-insensitive.
///
/// RNA U is read as T.
/// Returns kInvalidBase for every other character.
constexpr inline uint8_t base_encode(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case 'U': case 'u': return 3;
        default: return kInvalidBase;
    }
}

/// base_complement() - Watson-Crick complement of a 2-bit base.
///
/// A=00 and T=11 are bitwise complements, as are C=01 and G=10, so the complement is XOR with 3.
constexpr inline uint8_t base_complement(uint8_t b) noexcept {
    return b ^ 3;
}

/// encode_kmer() - pack the first @k characters of @seq into @out.
///
/// The first base lands in the highest 2 bits, so "ACGT" packs to 0b00'01'10'11 = 0x1B.
/// Returns false, with @out partially written, if any of the k characters is not a valid base.
constexpr inline bool encode_kmer(const char* seq, std::size_t k, KmerWord& out) noexcept {
    out = 0;
    for (std::size_t i = 0; i < k; ++i) {
        uint8_t b = base_encode(seq[i]);
        if (b == kInvalidBase) return false;
        out = (out << 2) | b;
    }
    return true;
}

/// reverse_complement() - reverse complement of a packed k-mer, O(1).
///
/// Complementing every base is a bitwise NOT.
/// Reversing the 32 two-bit fields of the word is a byte swap followed by swapping the nibbles and then the 2-bit pairs inside each byte.
/// The k bases then sit in the top 2k bits and are shifted down.
/// Returns 0 for k = 0.
constexpr inline KmerWord reverse_complement(KmerWord kmer, std::size_t k) noexcept {
    uint64_t x = ~static_cast<uint64_t>(kmer);
    x = __builtin_bswap64(x);
    x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    return k == 0 ? 0 : static_cast<KmerWord>(x >> (64 - 2 * k));
}

/// kmer_mask() - mask of the low 2k bits of a KmerWord.
///
/// For k = kMaxK the mask is all ones, since a shift by 64 would be undefined.
/// Throws std::invalid_argument for k > kMaxK.
constexpr inline KmerWord kmer_mask(std::size_t k) {
    if (k > kMaxK)
        throw std::invalid_argument(
            "k=" + std::to_string(k) + " exceeds MAX_K=" +
            std::to_string(kMaxK) + "; 2k bits exceed 64-bit KmerWord");
    return (k < kMaxK) ? ((KmerWord(1) << (2 * k)) - 1) : ~KmerWord(0);
}

// The packing conventions above, checked at compile time.
static_assert(base_complement(0) == 3 && base_complement(1) == 2, "A pairs with T, C with G");
static_assert(kmer_mask(31) == (KmerWord(1) << 62) - 1 && kmer_mask(kMaxK) == ~KmerWord(0),
              "kmer_mask covers 2k bits");
static_assert(reverse_complement(0x01, 3) == 0x2F, "AAC reverse-complements to GTT");

}
