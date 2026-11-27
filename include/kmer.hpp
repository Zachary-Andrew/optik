#pragma once

/**
 * @file
 * @brief 2-bit packed k-mer encoding.
 *
 * Bases take two bits each (A=0, C=1, G=2, T=3), so up to kMaxK = 32 of them fit in a 64-bit word.
 * Extraction shifts left as bases arrive, which leaves the first base in the top occupied bits [2k-2, 2k).
 *
 * With this encoding, complementing a base is an XOR with 3.
 * Likewise, reverse-complementing a whole word is a bit reversal over 2-bit fields.
 */

#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace tpoptoa {

static constexpr unsigned kBitsPerBase = 2;
/// XOR with this complements a 2-bit base (A<->T, C<->G).
static constexpr uint8_t kComplementMask = 3;
/// Selects one 2-bit base code.
static constexpr unsigned kBaseMask = (1u << kBitsPerBase) - 1;
/// Code-to-character table, so kBaseChars[code] is the base letter.
inline constexpr char kBaseChars[] = "ACGT";
static constexpr unsigned kKmerWordBits = 64;
static constexpr std::size_t kMaxK = kKmerWordBits / kBitsPerBase;
/// The graph builders store (k+1)-mer edges, which need one base of headroom.
static constexpr std::size_t kDbgMaxK = kMaxK - 1;
static constexpr unsigned kNumBases = 4;
/// The 2-bit code of each base, in the order of kBaseChars.
static constexpr uint8_t kCodeA = 0, kCodeC = 1, kCodeG = 2, kCodeT = 3;
/// A k-mer can be read on either strand, so a non-palindromic one has two orientations.
static constexpr unsigned kStrands = 2;
static constexpr unsigned kBasesPerByte = 8 / kBitsPerBase;
/// Entries in a table indexed by a raw byte.
static constexpr int kByteValues = 256;
/// A graph edge mask keeps out-edges in its low nibble (one bit per next base) and in-edges in its high nibble.
static constexpr unsigned kOutEdgeBits = 4;
static constexpr unsigned kOutEdgeMask = (1u << kOutEdgeBits) - 1;
static constexpr unsigned kEdgeMaskBits = 2 * kOutEdgeBits;

/// Marks anything that isn't A/C/G/T/U, such as N, IUPAC codes or gaps.
/// Callers treat it as a sequence break.
static constexpr uint8_t kInvalidBase = 0xFF;

using KmerWord = uint64_t;

/// Case-insensitive, and RNA's U counts as T.
constexpr inline uint8_t base_encode(char ch) noexcept {
    switch (ch) {
        case 'A': case 'a': return kCodeA;
        case 'C': case 'c': return kCodeC;
        case 'G': case 'g': return kCodeG;
        case 'T': case 't': return kCodeT;
        case 'U': case 'u': return kCodeT;
        default: return kInvalidBase;
    }
}

/// A=00/T=11 and C=01/G=10 are bitwise opposites, hence the XOR.
constexpr inline uint8_t base_complement(uint8_t base) noexcept {
    return base ^ kComplementMask;
}

/// Packs the first @p kmer_len characters with the first base highest, so "ACGT" becomes 0b00'01'10'11 = 0x1B.
/// On an invalid base it returns false and leaves @p out half-written.
constexpr inline bool encode_kmer(const char* seq, std::size_t kmer_len, KmerWord& out) noexcept {
    out = 0;

    for (std::size_t i = 0; i < kmer_len; ++i) {
        uint8_t base = base_encode(seq[i]);

        if (base == kInvalidBase) return false;

        out = (out << kBitsPerBase) | base;
    }

    return true;
}

/// NOT complements every base at once.
/// Reversing the 32 two-bit fields then takes a byte swap followed by nibble and pair swaps within each byte.
/// Finally, the k bases sit at the top and one shift brings them down.
constexpr inline KmerWord reverse_complement(KmerWord kmer, std::size_t kmer_len) noexcept {
    uint64_t bits = ~static_cast<uint64_t>(kmer);

    bits = __builtin_bswap64(bits);
    constexpr unsigned kNibbleBits = 4;
    constexpr uint64_t kLowNibbles = 0x0F0F0F0F0F0F0F0FULL, kLowPairs = 0x3333333333333333ULL;
    bits = ((bits >> kNibbleBits) & kLowNibbles) | ((bits & kLowNibbles) << kNibbleBits);
    bits = ((bits >> kBitsPerBase) & kLowPairs) | ((bits & kLowPairs) << kBitsPerBase);
    return kmer_len == 0 ? 0 : static_cast<KmerWord>(bits >> (kKmerWordBits - kBitsPerBase * kmer_len));
}

/// Addressing for arrays of 2-bit base codes packed kBasesPerByte to a byte, first base in the low bits.
constexpr inline std::size_t packed_byte(uint64_t pos) noexcept { return static_cast<std::size_t>(pos / kBasesPerByte); }
constexpr inline unsigned packed_shift(uint64_t pos) noexcept { return static_cast<unsigned>(pos % kBasesPerByte) * kBitsPerBase; }
constexpr inline std::size_t packed_bytes_for(uint64_t bases) noexcept {
    return static_cast<std::size_t>((bases + kBasesPerByte - 1) / kBasesPerByte);
}

/// Low 2k bits set. k = kMaxK is special-cased because shifting by 64 is undefined behaviour.
constexpr inline KmerWord kmer_mask(std::size_t kmer_len) {
    if (kmer_len > kMaxK)
        throw std::invalid_argument(
            "k=" + std::to_string(kmer_len) + " exceeds MAX_K=" +
            std::to_string(kMaxK) + "; 2k bits exceed 64-bit KmerWord");

    return (kmer_len < kMaxK) ? ((KmerWord(1) << (kBitsPerBase * kmer_len)) - 1) : ~KmerWord(0);
}

static_assert(base_complement(0) == 3 && base_complement(1) == 2, "A pairs with T, C with G");
static_assert(kmer_mask(31) == (KmerWord(1) << 62) - 1 && kmer_mask(kMaxK) == ~KmerWord(0),
              "kmer_mask covers 2k bits");
static_assert(reverse_complement(0x01, 3) == 0x2F, "AAC reverse-complements to GTT");

}
