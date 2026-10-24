// minimizer_stress.cpp - brute-force oracle for extract_minimizers, checking true sequence offsets, N-segment resets, hashed order, strand flag.
//   g++ -std=c++17 -O2 -I include bench/cpp/minimizer_stress.cpp && ./a.out
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include <tuple>
#include "minimizer.hpp"
using namespace tpoptoa;
int main() {
    std::mt19937_64 rng(3); int bad = 0; long tot = 0;
    for (int it = 0; it < 3000; ++it) {
        std::size_t k = 3 + rng() % 29, w = 1 + rng() % 20, n = rng() % 400;
        std::string s(n, 'A');
        for (auto& c : s) { unsigned r = rng() % 100; c = r < 2 ? 'N' : "ACGTacgt"[rng() % 8]; if (r >= 2 && r < 30) c = 'A'; }
        std::vector<std::tuple<uint64_t,uint32_t,bool>> got, want;
        extract_minimizers(s.data(), s.size(), k, w, [&](MinimizerSeed m){ got.emplace_back(m.kmer, m.pos, m.rev); });
        // brute force
        KmerWord mask = kmer_mask(k);
        std::vector<long> seg; // per position, whether the k-mer is valid
        struct KI { uint64_t can, h; bool rev; bool ok; };
        std::vector<KI> km(n >= k ? n - k + 1 : 0);
        for (std::size_t p = 0; p + k <= n; ++p) {
            KmerWord f; if (!encode_kmer(s.data() + p, k, f)) { km[p].ok = false; continue; }
            KmerWord r = reverse_complement(f, k); bool rv = r < f; KmerWord c = rv ? r : f;
            km[p] = {c, minimizer_order_hash(c, mask), rv, true};
        }
        long last = -1;
        for (std::size_t p = 0; p + w <= km.size(); ++p) {
            bool all = true; for (std::size_t j = p; j < p + w; ++j) all &= km[j].ok;
            if (!all) continue;
            std::size_t best = p; for (std::size_t j = p; j < p + w; ++j) if (km[j].h < km[best].h) best = j;
            if ((long)best != last) { want.emplace_back(km[best].can, best, km[best].rev); last = best; }
        }
        tot += want.size();
        if (got != want) { if (++bad < 4) printf("mismatch k=%zu w=%zu n=%zu got=%zu want=%zu\n", k, w, n, got.size(), want.size()); }
    }
    printf("minimizer brute-force: %d bad cases, %ld minimizers checked\n", bad, tot);
    return bad != 0;
}
