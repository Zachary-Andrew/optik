// dyn_spss_stress.cpp - DynSpss (static SPSS core plus OA deltas, the `spss-update` path) against a hash-set oracle, covering insert/erase return values, delta sizes, and after compaction the rebuilt core's membership and scans.
//   g++ -std=c++17 -O2 -march=native -fopenmp -pthread -I include -I src 
//       bench/cpp/dyn_spss_stress.cpp src/kmer_count.cpp && ./a.out
#include <cstdio>
#include <random>
#include <unordered_set>
#include "dyn_spss.hpp"
using namespace tpoptoa;

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("FAIL %s:%d ", __FILE__, __LINE__); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

int main() {
    std::mt19937_64 rng(5);
    const unsigned k = 21, m = 11;
    std::vector<std::string> seqs;
    for (int i = 0; i < 40; ++i) {
        std::string s; for (int j = 0; j < 400; ++j) s.push_back("ACGT"[rng() % 4]);
        seqs.push_back(s);
    }
    auto core = std::make_unique<SpssCore>();
    if (!core->build(seqs, k, m)) { std::printf("build failed\n"); return 1; }
    DynSpss d(std::move(core));
    std::unordered_set<uint64_t> o;
    auto canon = [&](uint64_t x) { uint64_t r = SpssCore::rc(x, k); return x < r ? x : r; };
    std::vector<uint64_t> pool;
    for (auto& s : seqs) for (std::size_t i = 0; i + k <= s.size(); ++i) {
        uint64_t x; encode_kmer(s.data() + i, k, x); o.insert(canon(x)); pool.push_back(x);
    }
    for (int i = 0; i < 20000; ++i) pool.push_back(rng() & kmer_mask(k));
    std::unordered_set<uint64_t> base = o;            // static set of the current core
    auto check_core = [&](const char* tag) {
        std::size_t b = 0;
        for (uint64_t x : pool) b += d.core().contains(x) != (o.count(canon(x)) != 0);
        // streaming scan over the sequences
        for (int i = 0; i < 40; ++i) {
            SpssCore::ScanStats st; d.core().scan(seqs[i].data(), seqs[i].size(), st);
            uint64_t want = 0;
            for (std::size_t p = 0; p + k <= seqs[i].size(); ++p) { uint64_t x; encode_kmer(seqs[i].data() + p, k, x); want += o.count(canon(x)); }
            b += st.positive != want;
        }
        CHECK(b == 0, "%s: %zu mismatches", tag, b);
    };
    check_core("initial");
    for (int round = 0; round < 3; ++round) {
        for (int i = 0; i < 20000; ++i) {
            uint64_t x = pool[rng() % pool.size()];
            if (rng() & 1) { bool was = !o.count(canon(x)); CHECK(d.insert(x) == was, "insert"); o.insert(canon(x)); }
            else { bool was = o.erase(canon(x)) == 1; CHECK(d.erase(x) == was, "erase"); }
        }
        std::size_t want_added = 0, want_deleted = 0;
        for (uint64_t c : o) want_added += !base.count(c);
        for (uint64_t c : base) want_deleted += !o.count(c);
        CHECK(d.added() == want_added && d.deleted() == want_deleted,
              "delta sizes added=%zu/%zu deleted=%zu/%zu", d.added(), want_added,
              d.deleted(), want_deleted);
        CHECK(d.compact(4), "compact");
        CHECK(d.added() == 0 && d.deleted() == 0, "deltas not cleared");
        CHECK(d.core().stats().kmers == o.size(), "compacted core holds %llu k-mers, want %zu",
              (unsigned long long)d.core().stats().kmers, o.size());
        check_core("after compact");
        base = o;
    }
    std::printf(g_fail ? "dyn_spss_stress: %d FAILURE(S)\n" : "dyn_spss_stress: ALL PASS\n", g_fail);
    return g_fail != 0;
}
