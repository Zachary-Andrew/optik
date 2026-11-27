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
    const std::size_t num_ops = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200000;
    {   // 1. oracle
        DeltaOA delta(64);
        std::unordered_set<uint64_t> oracle;
        std::mt19937_64 rng(1);
        std::size_t bad = 0;
        for (std::size_t i = 0; i < 4 * num_ops; ++i) {
            const uint64_t key = rng() % (num_ops * 2);
            switch (rng() % 3) {
            case 0: bad += delta.insert(key) != oracle.insert(key).second; break;
            case 1: bad += delta.erase(key) != (oracle.erase(key) == 1); break;
            default: bad += delta.contains(key) != (oracle.count(key) == 1);
            }
        }
        bad += delta.size() != oracle.size();
        std::size_t seen = 0; delta.for_each([&](uint64_t key) { seen += oracle.count(key); });
        bad += seen != oracle.size();
        std::printf("oracle: ops=%zu size=%zu bad=%zu bytes/key=%.1f\n", 4 * num_ops, delta.size(), bad,
                    delta.bytes() / double(std::max<std::size_t>(1, delta.size())));
        fail += bad != 0;
    }
    {   // 2. concurrency
        DeltaOA delta(1024);
        const unsigned num_writers = 4, num_readers = 4;
        const uint64_t STABLE = 1ull << 40, ABSENT = 1ull << 41;
        for (uint64_t i = 0; i < num_ops / 4; ++i) delta.insert(STABLE + i);
        std::atomic<bool> stop{false};
        std::atomic<uint64_t> bad{0}, reads{0};
        std::vector<std::thread> reader_threads;
        for (unsigned r = 0; r < num_readers; ++r)
            reader_threads.emplace_back([&, r] {
                std::mt19937_64 rng(100 + r); uint64_t reads_done = 0;
                while (!stop.load(std::memory_order_relaxed)) {
                    if (!delta.contains(STABLE + rng() % (num_ops / 4))) bad.fetch_add(1);
                    if (delta.contains(ABSENT + rng() % num_ops)) bad.fetch_add(1);
                    reads_done += 2;
                }
                reads.fetch_add(reads_done);
            });
        std::vector<std::thread> writer_threads;
        for (unsigned w = 0; w < num_writers; ++w)
            writer_threads.emplace_back([&, w] {
                const uint64_t base = (uint64_t(w) + 1) << 32;
                for (int round = 0; round < 3; ++round) {
                    for (uint64_t i = 0; i < num_ops; ++i) if (!delta.insert(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 0; i < num_ops; i += 2) if (!delta.erase(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 1; i < num_ops; i += 2) if (!delta.contains(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 0; i < num_ops; i += 2) if (delta.contains(base + i)) bad.fetch_add(1);
                    for (uint64_t i = 1; i < num_ops; i += 2) delta.erase(base + i);
                }
            });
        for (auto& thread : writer_threads) thread.join();
        stop.store(true);
        for (auto& thread : reader_threads) thread.join();
        const bool size_ok = delta.size() == num_ops / 4;
        std::printf("concurrency: reads=%llu bad=%llu final_size_ok=%d\n",
                    (unsigned long long)reads.load(), (unsigned long long)bad.load(), size_ok);
        fail += bad.load() != 0 || !size_ok;
    }
    std::printf(fail ? "delta_stress: %d FAILURE(S)\n" : "delta_stress: ALL PASS\n", fail);
    return fail != 0;
}
