#pragma once

/*
 * cache_info.hpp - CPU cache sizes from sysfs, which any user may read.
 *
 * Code that changes strategy with the working set, such as FmIndex's batched lookups and the micro benchmark's cache- and DRAM-sized cases, reads the sizes here instead of assuming one machine's caches.
 */

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>

namespace tpoptoa {

/// Cache sizes assumed when sysfs lists none, typical of current x86-64 CPUs.
inline constexpr std::size_t kFallbackL1dBytes = std::size_t(32) << 10;
inline constexpr std::size_t kFallbackLlcBytes = std::size_t(32) << 20;
/// Sysfs cache entries examined per CPU (index0, index1, ...).
inline constexpr int kMaxCacheIndex = 16;

/// cache_bytes() - size of the largest data or unified cache of @cpu at
/// @level, or of the highest level present when @level is 0, from /sys/devices/system/cpu/cpuN/cache.
/// Returns 0 when sysfs lists none.
inline std::size_t cache_bytes(int cpu, int level) {
    std::size_t best = 0;
    int best_level = 0;
    for (int i = 0; i < kMaxCacheIndex; ++i) {
        const std::string dir = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/cache/index" +
                                std::to_string(i) + "/";
        std::ifstream lf(dir + "level"), tf(dir + "type"), sf(dir + "size");
        int lv = 0;
        std::string type, size;
        if (!(lf >> lv) || !(tf >> type) || !(sf >> size)) continue;
        if (type == "Instruction" || (level != 0 && lv != level)) continue;
        std::size_t bytes = std::strtoull(size.c_str(), nullptr, 10);
        const char unit = size.empty() ? ' ' : size.back();
        if (unit == 'K') bytes <<= 10;
        else if (unit == 'M') bytes <<= 20;
        else if (unit == 'G') bytes <<= 30;
        if (lv > best_level || (lv == best_level && bytes > best)) { best = bytes; best_level = lv; }
    }
    return best;
}

/// last_level_cache_bytes() - the last-level cache of CPU 0, read once, or kFallbackLlcBytes when sysfs lists none.
inline std::size_t last_level_cache_bytes() {
    static const std::size_t llc = [] {
        const std::size_t b = cache_bytes(0, 0);
        return b ? b : kFallbackLlcBytes;
    }();
    return llc;
}

}
