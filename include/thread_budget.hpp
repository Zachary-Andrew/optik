#pragma once

/**
 * @file
 * @brief The default worker-thread count.
 *
 * A thread budget counts everything a call keeps runnable at once: workers, input readers and OpenMP regions.
 * Readers count too, because each one inflates its own .gz input.
 * The CLI resolves -t the same way (default_threads() in tools/src/main.rs).
 * So a library caller passing 0 gets the CLI default rather than every core.
 */

#include <algorithm>
#include <exception>
#include <thread>
#include <vector>

namespace tpoptoa {

/// @p requested if positive, otherwise half the logical CPUs and at least 1.
inline unsigned resolve_threads(unsigned requested) {
    if (requested > 0) return requested;

    constexpr unsigned kLogicalCpusPerThread = 2;   // one thread per physical core on an SMT machine

    return std::max(1u, std::thread::hardware_concurrency() / kLogicalCpusPerThread);
}

/// Runs @p func(t) for t in [0, @p count) on separate threads, with the caller doing t = 0 itself.
/// Because the caller works instead of sitting in join(), a phase never has more than @p count runnable threads.
/// An exception from the caller's share is rethrown only after the others have joined.
template <typename F>
void run_on_threads(unsigned count, F&& func) {
    if (count == 0) return;

    std::vector<std::thread> threads;

    threads.reserve(count - 1);

    for (unsigned t = 1; t < count; ++t) threads.emplace_back([&func, t] { func(t); });

    std::exception_ptr err;

    try {
        func(0u);
    } catch (...) {
        err = std::current_exception();
    }
    for (auto& thread : threads) thread.join();

    if (err) std::rethrow_exception(err);
}

}
