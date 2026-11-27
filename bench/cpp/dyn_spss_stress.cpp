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
    const unsigned kmer_len = 21, minimizer_len = 11;
    std::vector<std::string> seqs;
    for (int i = 0; i < 40; ++i) {
        std::string seq; for (int j = 0; j < 400; ++j) seq.push_back("ACGT"[rng() % 4]);
        seqs.push_back(seq);
    }
    auto core = std::make_unique<SpssCore>();
    if (!core->build(seqs, kmer_len, minimizer_len)) { std::printf("build failed\n"); return 1; }
    DynSpss dyn(std::move(core));
    std::unordered_set<uint64_t> oracle;
    auto canon = [&](uint64_t kmer) { uint64_t rev = SpssCore::revcomp_packed(kmer, kmer_len); return kmer < rev ? kmer : rev; };
    std::vector<uint64_t> pool;
    for (auto& seq : seqs) for (std::size_t i = 0; i + kmer_len <= seq.size(); ++i) {
        uint64_t kmer; encode_kmer(seq.data() + i, kmer_len, kmer); oracle.insert(canon(kmer)); pool.push_back(kmer);
    }
    for (int i = 0; i < 20000; ++i) pool.push_back(rng() & kmer_mask(kmer_len));
    std::unordered_set<uint64_t> base = oracle;            // static set of the current core
    auto check_core = [&](const char* tag) {
        std::size_t mismatches = 0;
        for (uint64_t kmer : pool) mismatches += dyn.core().contains(kmer) != (oracle.count(canon(kmer)) != 0);
        // streaming scan over the sequences
        for (int i = 0; i < 40; ++i) {
            SpssCore::ScanStats stats; dyn.core().scan(seqs[i].data(), seqs[i].size(), stats);
            uint64_t want = 0;
            for (std::size_t p = 0; p + kmer_len <= seqs[i].size(); ++p) { uint64_t kmer; encode_kmer(seqs[i].data() + p, kmer_len, kmer); want += oracle.count(canon(kmer)); }
            mismatches += stats.positive != want;
        }
        CHECK(mismatches == 0, "%s: %zu mismatches", tag, mismatches);
    };
    check_core("initial");
    for (int round = 0; round < 3; ++round) {
        for (int i = 0; i < 20000; ++i) {
            uint64_t kmer = pool[rng() % pool.size()];
            if (rng() & 1) { bool was = !oracle.count(canon(kmer)); CHECK(dyn.insert(kmer) == was, "insert"); oracle.insert(canon(kmer)); }
            else { bool was = oracle.erase(canon(kmer)) == 1; CHECK(dyn.erase(kmer) == was, "erase"); }
        }
        std::size_t want_added = 0, want_deleted = 0;
        for (uint64_t kmer : oracle) want_added += !base.count(kmer);
        for (uint64_t kmer : base) want_deleted += !oracle.count(kmer);
        CHECK(dyn.added() == want_added && dyn.deleted() == want_deleted,
              "delta sizes added=%zu/%zu deleted=%zu/%zu", dyn.added(), want_added,
              dyn.deleted(), want_deleted);
        CHECK(dyn.compact(4), "compact");
        CHECK(dyn.added() == 0 && dyn.deleted() == 0, "deltas not cleared");
        CHECK(dyn.core().stats().kmers == oracle.size(), "compacted core holds %llu k-mers, want %zu",
              (unsigned long long)dyn.core().stats().kmers, oracle.size());
        check_core("after compact");
        base = oracle;
    }
    std::printf(g_fail ? "dyn_spss_stress: %d FAILURE(S)\n" : "dyn_spss_stress: ALL PASS\n", g_fail);
    return g_fail != 0;
}
