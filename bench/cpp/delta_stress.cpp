// delta_stress.cpp - DeltaOA single-threaded oracle and concurrent stress.
//
//   g++ -std=c++17 -O2 -march=native -pthread -I include bench/cpp/delta_stress.cpp && ./a.out
//   (add -fsanitize=thread -O1 -g for the race check)
//
// 1. Oracle, a random insert/erase/contains sequence against std::unordered_set, with enough keys to force several resizes.
// 2. Concurrency, where writer threads insert, erase and re-insert disjoint key ranges while reader threads continuously check that a STABLE key set is always present and an ABSENT set never is.
#include <atomic>
#include <cstdio>
#include <random>
#include <thread>
#include <unordered_set>
#include <vector>
#include "delta_oa.hpp"
using namespace tpoptoa;

int main(int argc, char** argv) {
    int fail = 0;
    const std::size_t N = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200000;
    {   // 1. oracle
        DeltaOA d(64);
        std::unordered_set<uint64_t> o;
        std::mt19937_64 rng(1);
        std::size_t bad = 0;
        for (std::size_t i = 0; i < 4 * N; ++i) {
            const uint64_t k = rng() % (N * 2);
            switch (rng() % 3) {
            case 0: bad += d.insert(k) != o.insert(k).second; break;
            case 1: bad += d.erase(k) != (o.erase(k) == 1); break;
            default: bad += d.contains(k) != (o.count(k) == 1);
            }
        }
        bad += d.size() != o.size();
        std::size_t seen = 0; d.for_each([&](uint64_t k) { seen += o.count(k); });
        bad += seen != o.size();
        std::printf("oracle: ops=%zu size=%zu bad=%zu bytes/key=%.1f\n", 4 * N, d.size(), bad,
                    d.bytes() / double(std::max<std::size_t>(1, d.size())));
        fail += bad != 0;
    }
    {   // 2. concurrency
        DeltaOA d(1024);
        const unsigned W = 4, R = 4;
        const uint64_t STABLE = 1ull << 40, ABSENT = 1ull << 41;
        for (uint64_t i = 0; i < N / 4; ++i) d.insert(STABLE + i);
        std::atomic<bool> stop{false};
        std::atomic<uint64_t> bad{0}, reads{0};
        std::vector<std::thread> th;
        for (unsigned r = 0; r < R; ++r)
            th.emplace_back([&, r] {
                std::mt19937_64 rng(100 + r); uint64_t n = 0;
                while (!stop.load(std::memory_order_relaxed)) {
                    if (!d.contains(STABLE + rng() % (N / 4))) bad.fetch_add(1);
                    if (d.contains(ABSENT + rng() % N)) bad.fetch_add(1);
                    n += 2;
                }
                reads.fetch_add(n);
            });
        std::vector<std::thread> wt;
        for (unsigned w = 0; w < W; ++w)
            wt.emplace_back([&, w] {
                const uint64_t base = (uint64_t(w) + 1) << 32;
                for (int round = 0; round < 3; ++round) {
                    for (uint64_t i = 0; i < N; ++i) if (!d.insert(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 0; i < N; i += 2) if (!d.erase(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 1; i < N; i += 2) if (!d.contains(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 0; i < N; i += 2) if (d.contains(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 1; i < N; i += 2) d.erase(base + i);
                }
            });
        for (auto& t : wt) t.join();
        stop.store(true);
        for (auto& t : th) t.join();
        const bool size_ok = d.size() == N / 4;
        std::printf("concurrency: reads=%llu bad=%llu final_size_ok=%d\n",
                    (unsigned long long)reads.load(), (unsigned long long)bad.load(), size_ok);
        fail += bad.load() != 0 || !size_ok;
    }
    std::printf(fail ? "delta_stress: %d FAILURE(S)\n" : "delta_stress: ALL PASS\n", fail);
    return fail != 0;
}
