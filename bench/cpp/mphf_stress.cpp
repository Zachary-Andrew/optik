// mphf_stress.cpp - PilotMphf as a bijection onto [0, n), save/load and timing.
//   g++ -std=c++17 -O2 -march=native -I include bench/cpp/mphf_stress.cpp && ./a.out
#include <chrono>
#include <cstdio>
#include <random>
#include <unordered_set>
#include <vector>
#include "phobic_mphf.hpp"
using namespace tpoptoa;
int main() {
    int fail = 0;
    std::mt19937_64 rng(9);
    for (std::size_t num_keys : {1ul, 2ul, 7ul, 100ul, 1000ul, 100000ul, 2000000ul}) {
        std::unordered_set<uint64_t> key_set;
        while (key_set.size() < num_keys) key_set.insert(rng());
        std::vector<uint64_t> keys(key_set.begin(), key_set.end());
        PilotMphf mphf;
        auto start = std::chrono::steady_clock::now();
        if (!mphf.build(keys.data(), num_keys)) { std::printf("n=%zu build failed\n", num_keys); ++fail; continue; }
        double build_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::vector<uint8_t> seen(num_keys, 0);
        std::size_t bad = 0;
        for (uint64_t key : keys) { uint64_t rank = mphf(key); if (rank >= num_keys || seen[rank]++) ++bad; }
        std::vector<uint8_t> blob; mphf.save(blob);
        PilotMphf reloaded; std::size_t used = reloaded.load(blob.data(), blob.size());
        std::size_t diff = (used != blob.size());
        for (uint64_t key : keys) diff += mphf(key) != reloaded(key);
        // A borrowing load reads the pilots and the remap in place, and a copy of an owned or a borrowing function must keep answering through its own views.
        PilotMphf borrowed; diff += borrowed.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
        const PilotMphf owned_copy = mphf, borrowed_copy = borrowed;
        for (uint64_t key : keys) diff += (borrowed(key) != mphf(key)) + (owned_copy(key) != mphf(key)) + (borrowed_copy(key) != mphf(key));
        // lookup timing in random order
        std::shuffle(keys.begin(), keys.end(), rng);
        start = std::chrono::steady_clock::now();
        uint64_t acc = 0;
        for (uint64_t key : keys) acc += mphf(key);
        double lookup_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::printf("n=%-8zu bad=%zu reload_diff=%zu bits/key=%.2f build=%.3fs lookup=%.1fns (%llu)\n",
                    num_keys, bad, diff, mphf.bytes() * 8.0 / num_keys, build_secs, lookup_secs * 1e9 / num_keys, (unsigned long long)(acc & 1));
        fail += (bad != 0) + (diff != 0);
    }
    // Partitioned MPHF bijection, reload and parallel build speed
    for (std::size_t num_keys : {10ul, 300000ul, 2000000ul}) {
        std::unordered_set<uint64_t> key_set;
        while (key_set.size() < num_keys) key_set.insert(rng());
        std::vector<uint64_t> keys(key_set.begin(), key_set.end());
        Mphf mphf;
        auto start = std::chrono::steady_clock::now();
        if (!mphf.build(keys.data(), num_keys)) { std::printf("partitioned n=%zu build failed\n", num_keys); ++fail; continue; }
        double build_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::vector<uint8_t> seen(num_keys, 0);
        std::size_t bad = 0;
        for (uint64_t key : keys) { uint64_t rank = mphf(key); if (rank >= num_keys || seen[rank]++) ++bad; }
        std::vector<uint8_t> blob; mphf.save(blob);
        Mphf reloaded; std::size_t diff = reloaded.load(blob.data(), blob.size()) != blob.size();
        for (uint64_t key : keys) diff += mphf(key) != reloaded(key);
        Mphf borrowed; diff += borrowed.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
        for (uint64_t key : keys) diff += borrowed(key) != mphf(key);
        std::printf("partitioned n=%-8zu bad=%zu reload_diff=%zu bits/key=%.2f build=%.3fs\n",
                    num_keys, bad, diff, mphf.bytes() * 8.0 / num_keys, build_secs);
        fail += (bad != 0) + (diff != 0);
    }
    // duplicate keys must be rejected
    std::vector<uint64_t> dup_keys = {5, 7, 5};
    PilotMphf mphf; if (mphf.build(dup_keys.data(), dup_keys.size())) { std::printf("duplicates accepted\n"); ++fail; }
    std::printf(fail ? "mphf_stress: %d FAILURE(S)\n" : "mphf_stress: ALL PASS\n", fail);
    return fail != 0;
}
