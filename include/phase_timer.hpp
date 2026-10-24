#pragma once

/*
 * phase_timer.hpp - per-phase wall time and memory on stderr.
 *
 * Enabled by setting OPTIK_TIMING in the environment.
 * Each lap() prints the time since the previous lap, the current resident set (/proc/self/statm) and the process peak (getrusage), so a build's peak memory can be attributed to the phase that reached it.
 */

#include <chrono>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>

namespace tpoptoa {

struct PhaseTimer {
    static constexpr double kBytesPerMb = 1048576.0;
    static constexpr double kKbPerMb = 1024.0;

    bool on = std::getenv("OPTIK_TIMING") != nullptr;
    std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();

    /// lap() - report the phase that just ended, named @what.
    void lap(const char* what) {
        if (!on) return;
        auto n = std::chrono::steady_clock::now();
        long pages = 0, resident = 0;  // statm fields, in pages
        if (std::FILE* f = std::fopen("/proc/self/statm", "r")) {
            if (std::fscanf(f, "%ld %ld", &pages, &resident) != 2) resident = 0;
            std::fclose(f);
        }
        struct rusage ru {};
        getrusage(RUSAGE_SELF, &ru);  // ru_maxrss is in KB on Linux
        std::fprintf(stderr, "[timing] %-28s %.3f s  rss %.0f MB  peak %.0f MB\n", what,
                     std::chrono::duration<double>(n - t).count(),
                     resident * (sysconf(_SC_PAGESIZE) / kBytesPerMb), ru.ru_maxrss / kKbPerMb);
        t = n;
    }
};

}
