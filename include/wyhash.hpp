#pragma once

/**
 * @file
 * @brief 64-bit integer hashes.
 *
 * wymix() and its secret come from Wang Yi's public-domain wyhash (github.com/wangyi-fudan/wyhash).
 * Checksumming SPSS index files (ffi_spss.cpp) needs only the multiply-fold step, so the streaming and PRNG parts were left out.
 * Consequently, our checksum won't match reference wyhash output.
 *
 * mix64() is MurmurHash3's finalizer and serves as the seeded hash for every k-mer table.
 */

#include <cstdint>

#include "bits.hpp"

namespace tpoptoa {

/// XORs the high half of the 128-bit product into the low half.
/// Since each half depends on every input bit, one multiplication is enough for full avalanche.
constexpr inline uint64_t wymix(uint64_t lhs, uint64_t rhs) noexcept {
    unsigned __int128 product = lhs;

    product *= rhs;
    lhs = static_cast<uint64_t>(product);
    rhs = static_cast<uint64_t>(product >> kWordBits);
    return lhs ^ rhs;
}

/// The default secret of the reference wyhash.
inline constexpr uint64_t kWyhashDefaultSecret[4] = {
    UINT64_C(0xa0761d6478bd642f), UINT64_C(0xe7037ed1a0b428db),
    UINT64_C(0x8ebc6af09c88c6e3), UINT64_C(0x589965cc75374cc3),
};

/// 2^64 / phi, rounded to odd.
/// Multiplying by it scatters consecutive keys across the 64-bit range (Fibonacci hashing), and SplitMix64 uses it as its increment.
/// Index files store hashes seeded with this value, so changing it would break every saved index.
inline constexpr uint64_t kGoldenRatio64 = 0x9E3779B97F4A7C15ULL;

/// Computes fmix64(@p value ^ @p seed), so pass @p seed = 0 for plain fmix64.
/// The low bits make a safe table index even for keys that differ only in their high bits.
constexpr inline uint64_t mix64(uint64_t value, uint64_t seed) noexcept {
    // Appleby picked these by searching for the best avalanche behaviour.
    constexpr unsigned shift = 33;
    constexpr uint64_t first_multiplier = 0xff51afd7ed558ccdULL;
    constexpr uint64_t second_multiplier = 0xc4ceb9fe1a85ec53ULL;

    value ^= seed;
    value ^= value >> shift;
    value *= first_multiplier;
    value ^= value >> shift;
    value *= second_multiplier;
    value ^= value >> shift;
    return value;
}

// Saved indexes depend on this exact output.
static_assert(mix64(1, 0) == 0xb456bcfc34c2cb2cULL, "MurmurHash3 fmix64(1)");

}
