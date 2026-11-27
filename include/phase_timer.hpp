#pragma once

/**
 * @file
 * @brief Per-phase wall time and memory on stderr.
 *
 * Setting OPTIK_TIMING turns it on.
 * Each lap() then prints the time since the last lap, the current RSS and the process peak.
 * That is how we tell which build phase is responsible for peak memory.
 */

#include <chrono>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>

namespace tpoptoa {

struct PhaseTimer {
    static constexpr double kBytesPerMiB = 1024.0 * 1024.0;
    static constexpr double kKiBPerMiB = 1024.0;   // ru_maxrss is in KiB on Linux

    bool enabled = std::getenv("OPTIK_TIMING") != nullptr;
    std::chrono::steady_clock::time_point last = std::chrono::steady_clock::now();

    void lap(const char* what) {
        if (!enabled) return;

        auto now = std::chrono::steady_clock::now();
        long pages = 0, resident = 0;  // statm fields, in pages

        if (std::FILE* statm = std::fopen("/proc/self/statm", "r")) {
            constexpr int kStatmFields = 2;   // total and resident pages

        if (std::fscanf(statm, "%ld %ld", &pages, &resident) != kStatmFields) resident = 0;

            std::fclose(statm);
        }

        struct rusage usage {};
        getrusage(RUSAGE_SELF, &usage);
        std::fprintf(stderr, "[timing] %-28s %.3f s  rss %.0f MB  peak %.0f MB\n", what,
                     std::chrono::duration<double>(now - last).count(),
                     resident * (sysconf(_SC_PAGESIZE) / kBytesPerMiB), usage.ru_maxrss / kKiBPerMiB);
        last = now;
    }
};

}
