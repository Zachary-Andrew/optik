#pragma once

/*
 * thread_budget.hpp - the default worker-thread count.
 *
 * A thread budget counts every thread a call keeps runnable at once, i.e. workers, input readers (each inflates its own .gz in-thread) and OpenMP regions.
 * The CLI resolves -t the same way (default_threads() in tools/src/main.rs), so a library caller passing 0 gets the CLI default rather than every core.
 * run_on_threads() keeps a parallel phase inside its budget by having the calling thread take one share of the work, so a phase of n shares starts n - 1 threads and no more than n are runnable, even while they start.
 */

#include <algorithm>
#include <exception>
#include <thread>
#include <vector>

namespace tpoptoa {

/// resolve_threads() - @n if positive, else half the logical CPUs (at least 1).
inline unsigned resolve_threads(unsigned n) {
    if (n > 0) return n;
    return std::max(1u, std::thread::hardware_concurrency() / 2);
}

/// run_on_threads() - call @fn(t) for every t in [0, @n), each on its own thread, the calling thread taking t = 0, and return once all calls have.
///
/// The caller works instead of waiting in join(), so the phase never has more than @n runnable threads.
/// An exception from the caller's share is rethrown after the other threads are joined.
/// @n == 0 calls nothing.
template <typename F>
void run_on_threads(unsigned n, F&& fn) {
    if (n == 0) return;
    std::vector<std::thread> th;
    th.reserve(n - 1);
    for (unsigned t = 1; t < n; ++t) th.emplace_back([&fn, t] { fn(t); });
    std::exception_ptr err;
    try {
        fn(0u);
    } catch (...) {
        err = std::current_exception();
    }
    for (auto& x : th) x.join();
    if (err) std::rethrow_exception(err);
}

}
