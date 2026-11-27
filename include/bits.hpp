#pragma once

/**
 * @file
 * @brief Addressing helpers for bit vectors packed into 64-bit words.
 *
 * Several indexes (Elias-Fano keys, MPHF pilots, rank blocks, packed strings) store bits LSB-first in uint64_t arrays.
 * These helpers give the word/offset arithmetic the indexes share a readable format.
 * For unsigned operands the compiler turns the divisions into shifts and masks, so nothing is lost over hand-written ones.
 */

#include <cstddef>
#include <cstdint>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

namespace tpoptoa {

inline constexpr unsigned kWordBits = 64;
inline constexpr unsigned kBitsPerByte = 8;
/// Radix of decimal digit strings, such as CIGAR run lengths and numeric environment values.
inline constexpr int kDecimalBase = 10;
inline constexpr uint64_t kTopBit = uint64_t(1) << (kWordBits - 1);

/// __builtin_prefetch's access hints, read or write.
inline constexpr int kPrefetchRead = 0, kPrefetchWrite = 1;

/// __builtin_prefetch's locality hints, from data used once (none) to data used again soon (high).
inline constexpr int kLocalityNone = 0, kLocalityLow = 1, kLocalityHigh = 3;

/// Prefetches @p ptr for a read that comes soon, keeping the line in every cache level.
inline void prefetch_read(const void* ptr) noexcept { __builtin_prefetch(ptr, kPrefetchRead, kLocalityHigh); }

/// Prefetches @p ptr for a read that is less likely to recur, so the line may leave the inner caches first.
inline void prefetch_read_low(const void* ptr) noexcept { __builtin_prefetch(ptr, kPrefetchRead, kLocalityLow); }

/// Prefetches @p ptr for a single read, without keeping the line once it is used.
inline void prefetch_read_once(const void* ptr) noexcept { __builtin_prefetch(ptr, kPrefetchRead, kLocalityNone); }

/// Prefetches @p ptr for a write that comes soon.
inline void prefetch_write(const void* ptr) noexcept { __builtin_prefetch(ptr, kPrefetchWrite, kLocalityHigh); }

/// Prefetches @p ptr for a single write, without keeping the line once it is used.
inline void prefetch_write_once(const void* ptr) noexcept { __builtin_prefetch(ptr, kPrefetchWrite, kLocalityNone); }

/// Index of the word that holds bit @p bit.
constexpr inline std::size_t word_index(uint64_t bit) noexcept { return static_cast<std::size_t>(bit / kWordBits); }

/// Position of bit @p bit within its word, counted from the least significant end.
constexpr inline unsigned bit_offset(uint64_t bit) noexcept { return static_cast<unsigned>(bit % kWordBits); }

/// First bit of word @p word.
constexpr inline uint64_t first_bit_of_word(uint64_t word) noexcept { return word * kWordBits; }

/// Words needed to hold @p bits bits.
constexpr inline std::size_t words_for_bits(uint64_t bits) noexcept {
    return static_cast<std::size_t>((bits + kWordBits - 1) / kWordBits);
}

/// Rotates @p value left by @p shift bits, wrapping any shift of a whole word or more.
constexpr inline uint64_t rotate_left(uint64_t value, unsigned shift) noexcept {
    shift %= kWordBits;
    return shift ? (value << shift) | (value >> (kWordBits - shift)) : value;
}

/// floor(log2(@p value)) for @p value > 0.
constexpr inline unsigned floor_log2(uint64_t value) noexcept {
    return kWordBits - 1 - static_cast<unsigned>(__builtin_clzll(value));
}

/// Position of the @p rank-th set bit of @p word, counting from 0.
/// With BMI2 a single PDEP does it, otherwise the lowest set bits are cleared one at a time.
inline unsigned select_in_word(uint64_t word, unsigned rank) noexcept {
#if defined(__BMI2__)
    return static_cast<unsigned>(__builtin_ctzll(_pdep_u64(uint64_t(1) << rank, word)));
#else
    for (unsigned t = 0; t < rank; ++t) word &= word - 1;

    return static_cast<unsigned>(__builtin_ctzll(word));
#endif
}

/// Bytes needed to hold @p bits bits, rounded down, as the size estimates use it.
constexpr inline uint64_t bits_to_bytes(uint64_t bits) noexcept { return bits / kBitsPerByte; }

static_assert(word_index(130) == 2 && bit_offset(130) == 2, "bit 130 is bit 2 of word 2");
static_assert(words_for_bits(64) == 1 && words_for_bits(65) == 2, "partial words round up");

}
