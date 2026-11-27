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
        std::size_t kmer_len = 3 + rng() % 29, window = 1 + rng() % 20, seq_len = rng() % 400;
        std::string seq(seq_len, 'A');
        for (auto& ch : seq) { unsigned roll = rng() % 100; ch = roll < 2 ? 'N' : "ACGTacgt"[rng() % 8]; if (roll >= 2 && roll < 30) ch = 'A'; }
        std::vector<std::tuple<uint64_t,uint32_t,bool>> got, want;
        extract_minimizers(seq.data(), seq.size(), kmer_len, window, [&](MinimizerSeed seed){ got.emplace_back(seed.kmer, seed.pos, seed.rev); });
        // brute force
        KmerWord mask = kmer_mask(kmer_len);
        std::vector<long> seg; // per position, whether the k-mer is valid
        struct KmerInfo { uint64_t can, hash; bool rev; bool valid; };
        std::vector<KmerInfo> kmers(seq_len >= kmer_len ? seq_len - kmer_len + 1 : 0);
        for (std::size_t p = 0; p + kmer_len <= seq_len; ++p) {
            KmerWord fwd; if (!encode_kmer(seq.data() + p, kmer_len, fwd)) { kmers[p].valid = false; continue; }
            KmerWord rev = reverse_complement(fwd, kmer_len); bool is_rev = rev < fwd; KmerWord canon = is_rev ? rev : fwd;
            kmers[p] = {canon, minimizer_order_hash(canon, mask), is_rev, true};
        }
        long last = -1;
        for (std::size_t p = 0; p + window <= kmers.size(); ++p) {
            bool all = true; for (std::size_t j = p; j < p + window; ++j) all &= kmers[j].valid;
            if (!all) continue;
            std::size_t best = p; for (std::size_t j = p; j < p + window; ++j) if (kmers[j].hash < kmers[best].hash) best = j;
            if ((long)best != last) { want.emplace_back(kmers[best].can, best, kmers[best].rev); last = best; }
        }
        tot += want.size();
        if (got != want) { if (++bad < 4) printf("mismatch k=%zu w=%zu n=%zu got=%zu want=%zu\n", kmer_len, window, seq_len, got.size(), want.size()); }
    }
    printf("minimizer brute-force: %d bad cases, %ld minimizers checked\n", bad, tot);
    return bad != 0;
}
