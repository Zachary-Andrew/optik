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
    for (std::size_t n : {1ul, 2ul, 7ul, 100ul, 1000ul, 100000ul, 2000000ul}) {
        std::unordered_set<uint64_t> s;
        while (s.size() < n) s.insert(rng());
        std::vector<uint64_t> keys(s.begin(), s.end());
        PilotMphf f;
        auto t0 = std::chrono::steady_clock::now();
        if (!f.build(keys.data(), n)) { std::printf("n=%zu build failed\n", n); ++fail; continue; }
        double bt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::vector<uint8_t> seen(n, 0);
        std::size_t bad = 0;
        for (uint64_t k : keys) { uint64_t r = f(k); if (r >= n || seen[r]++) ++bad; }
        std::vector<uint8_t> blob; f.save(blob);
        PilotMphf g; std::size_t used = g.load(blob.data(), blob.size());
        std::size_t diff = (used != blob.size());
        for (uint64_t k : keys) diff += f(k) != g(k);
        // A borrowing load reads the pilots and the remap in place, and a copy of an owned or a borrowing function must keep answering through its own views.
        PilotMphf gb; diff += gb.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
        const PilotMphf fc = f, gc = gb;
        for (uint64_t k : keys) diff += (gb(k) != f(k)) + (fc(k) != f(k)) + (gc(k) != f(k));
        // lookup timing in random order
        std::shuffle(keys.begin(), keys.end(), rng);
        t0 = std::chrono::steady_clock::now();
        uint64_t acc = 0;
        for (uint64_t k : keys) acc += f(k);
        double lt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("n=%-8zu bad=%zu reload_diff=%zu bits/key=%.2f build=%.3fs lookup=%.1fns (%llu)\n",
                    n, bad, diff, f.bytes() * 8.0 / n, bt, lt * 1e9 / n, (unsigned long long)(acc & 1));
        fail += (bad != 0) + (diff != 0);
    }
    // Partitioned MPHF bijection, reload and parallel build speed
    for (std::size_t n : {10ul, 300000ul, 2000000ul}) {
        std::unordered_set<uint64_t> s;
        while (s.size() < n) s.insert(rng());
        std::vector<uint64_t> keys(s.begin(), s.end());
        Mphf f;
        auto t0 = std::chrono::steady_clock::now();
        if (!f.build(keys.data(), n)) { std::printf("partitioned n=%zu build failed\n", n); ++fail; continue; }
        double bt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::vector<uint8_t> seen(n, 0);
        std::size_t bad = 0;
        for (uint64_t k : keys) { uint64_t r = f(k); if (r >= n || seen[r]++) ++bad; }
        std::vector<uint8_t> blob; f.save(blob);
        Mphf g; std::size_t diff = g.load(blob.data(), blob.size()) != blob.size();
        for (uint64_t k : keys) diff += f(k) != g(k);
        Mphf gb; diff += gb.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
        for (uint64_t k : keys) diff += gb(k) != f(k);
        std::printf("partitioned n=%-8zu bad=%zu reload_diff=%zu bits/key=%.2f build=%.3fs\n",
                    n, bad, diff, f.bytes() * 8.0 / n, bt);
        fail += (bad != 0) + (diff != 0);
    }
    // duplicate keys must be rejected
    std::vector<uint64_t> d = {5, 7, 5};
    PilotMphf f; if (f.build(d.data(), d.size())) { std::printf("duplicates accepted\n"); ++fail; }
    std::printf(fail ? "mphf_stress: %d FAILURE(S)\n" : "mphf_stress: ALL PASS\n", fail);
    return fail != 0;
}
