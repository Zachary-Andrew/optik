#pragma once

/**
 * @file
 * @brief CPU cache sizes read from sysfs (readable without root).
 *
 * FmIndex's batched lookups and the micro benchmark's cache/DRAM cases size themselves from these.
 * That way neither assumes one particular machine.
 */

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>

#include "bits.hpp"

namespace tpoptoa {

// Typical x86-64 sizes, used when sysfs has nothing (containers, odd kernels).
inline constexpr std::size_t kFallbackL1dBytes = std::size_t(32) << 10;
inline constexpr std::size_t kFallbackLlcBytes = std::size_t(32) << 20;
// Real CPUs list 4 or 5 entries (index0, index1, ...), so this is a generous upper bound.
inline constexpr int kMaxSysfsCacheEntries = 16;

/// Largest data/unified cache of @p cpu at @p level, where level 0 means the highest one present.
/// Returns 0 when sysfs doesn't list a match.
inline std::size_t cache_bytes(int cpu, int level) {
    std::size_t best = 0;
    int best_level = 0;

    for (int i = 0; i < kMaxSysfsCacheEntries; ++i) {
        const std::string dir = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/cache/index" +
                                std::to_string(i) + "/";
        std::ifstream level_file(dir + "level"), type_file(dir + "type"), size_file(dir + "size");
        int cache_level = 0;
        std::string type, size;

        if (!(level_file >> cache_level) || !(type_file >> type) || !(size_file >> size)) continue;

        if (type == "Instruction" || (level != 0 && cache_level != level)) continue;

        std::size_t bytes = std::strtoull(size.c_str(), nullptr, kDecimalBase);

        // sysfs prints sizes like "32K" or "32M", in binary units.
        constexpr std::size_t kKiB = 1024, kMiB = kKiB * 1024, kGiB = kMiB * 1024;
        const char unit = size.empty() ? ' ' : size.back();
        if (unit == 'K') bytes *= kKiB;
        else if (unit == 'M') bytes *= kMiB;
        else if (unit == 'G') bytes *= kGiB;

        if (cache_level > best_level || (cache_level == best_level && bytes > best)) { best = bytes; best_level = cache_level; }
    }

    return best;
}

/// Last-level cache of CPU 0, read once.
inline std::size_t last_level_cache_bytes() {
    static const std::size_t llc = [] {
        const std::size_t detected = cache_bytes(0, 0);

        return detected ? detected : kFallbackLlcBytes;
    }();

    return llc;
}

}
