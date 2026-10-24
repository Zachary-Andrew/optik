#pragma once

/*
 * wyhash.hpp - 64-bit integer hashes.
 *
 * wymix() and the four secret constants are taken from wyhash by Wang Yi (public domain, https://github.com/wangyi-fudan/wyhash).
 * Only the multiply-fold mixer is kept.
 * The SPSS index checksum in ffi_spss.cpp folds file blocks with it, so the streaming and PRNG parts of wyhash are not included and the checksum is not output-compatible with the reference wyhash.
 *
 * mix64() is MurmurHash3's 64-bit finalizer, the seeded hash optik's tables use for k-mers and keys.
 */

#include <cstdint>

namespace tpoptoa {

/// wymix() - multiply-fold, the 128-bit product of @a and @b with the high half XORed into the low half.
///
/// Both halves of the product depend on every input bit, so folding them gives full avalanche from one multiplication.
constexpr inline uint64_t wymix(uint64_t a, uint64_t b) noexcept {
    unsigned __int128 r = a;
    r *= b;
    a = static_cast<uint64_t>(r);
    b = static_cast<uint64_t>(r >> 64);
    return a ^ b;
}

/// The default secret of the reference wyhash.
inline constexpr uint64_t kWyhashDefaultSecret[4] = {
    UINT64_C(0xa0761d6478bd642f), UINT64_C(0xe7037ed1a0b428db),
    UINT64_C(0x8ebc6af09c88c6e3), UINT64_C(0x589965cc75374cc3),
};

/// mix64() - MurmurHash3's 64-bit finalizer (fmix64) applied to @x ^ @seed.
constexpr inline uint64_t mix64(uint64_t x, uint64_t seed) noexcept {
    x ^= seed;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

// Hash tables built from mix64 are stored in index files, so its output is fixed.
static_assert(mix64(1, 0) == 0xb456bcfc34c2cb2cULL, "MurmurHash3 fmix64(1)");

}
