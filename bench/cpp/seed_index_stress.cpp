// seed_index_stress - the parallel chunked SeedIndex build must equal a single serial pass over each reference (the oracle runs extract_minimizers on the whole sequence, sorted and deduplicated), at several thread counts, with N runs, repeats and chunk boundaries (references > 1 Mbp).
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "seed_index.hpp"

using namespace tpoptoa;

int main() {
    std::mt19937_64 rng(7);
    const char* B = "ACGT";
    std::vector<std::string> seqs;
    for (std::size_t L : {2500000ul, 1048600ul, 5000ul, 20ul}) {
        std::string s(L, 'A');
        for (auto& c : s) c = B[rng() & 3];
        for (int r = 0; r < 40 && L > 10000; ++r) {          // repeats and N runs
            std::size_t a = rng() % (L - 3000), b = rng() % (L - 3000);
            s.replace(b, 2000, s.substr(a, 2000));
            std::size_t n = rng() % (L - 50);
            for (int t = 0; t < int(rng() % 40); ++t) s[n + t] = 'N';
        }
        seqs.push_back(s);
    }
    std::vector<std::pair<const char*, std::size_t>> refs;
    for (auto& s : seqs) refs.push_back({s.data(), s.size()});
    int fail = 0;
    for (unsigned k : {15u, 21u}) for (unsigned w : {5u, 10u, 50u}) {
        std::map<uint64_t, std::vector<uint64_t>> oracle;
        for (std::size_t rid = 0; rid < refs.size(); ++rid)
            extract_minimizers(refs[rid].first, refs[rid].second, k, w, [&](const MinimizerSeed& m) {
                oracle[m.kmer].push_back((uint64_t(rid) << 33) | (uint64_t(m.pos) << 1) | m.rev);
            });
        std::size_t occ = 0;
        for (auto& [key, v] : oracle) { std::sort(v.begin(), v.end()); occ += v.size(); }
        for (unsigned th : {1u, 3u, 8u}) {
            SeedIndex idx;
            if (!idx.build(refs, k, w, 1u << 30, th)) { std::printf("build failed\n"); ++fail; continue; }
            std::size_t bad = idx.distinct() != oracle.size() || idx.occurrences() != occ;
            for (auto& [key, v] : oracle) {
                auto [p, n] = idx.find(key);
                bad += !p || n != v.size() || !std::equal(v.begin(), v.end(), p);
            }
            std::size_t fp = 0;
            for (int i = 0; i < 20000; ++i) {
                uint64_t x = rng() & kmer_mask(k), xr = reverse_complement(x, k);
                x = std::min(x, xr);
                fp += !oracle.count(x) && idx.find(x).first;
            }
            // pipelined collect() == naive per-minimizer find() in read order
            std::size_t cbad = 0;
            for (int i = 0; i < 300; ++i) {
                const auto& src = seqs[rng() % 2];
                std::size_t L = 50 + rng() % 400, p0 = rng() % (src.size() - L);
                std::string rd = src.substr(p0, L);
                for (auto& c : rd) if (rng() % 50 == 0) c = B[rng() & 3];
                std::vector<SeedIndex::Hit> got, want;
                idx.collect(rd.data(), rd.size(), got);
                extract_minimizers(rd.data(), rd.size(), k, w, [&](const MinimizerSeed& s) {
                    auto [pp, nn] = idx.find(s.kmer);
                    if (!pp) return;
                    for (uint32_t t = 0; t < nn; ++t) {
                        const bool rev = ((pp[t] & 1) != 0) != s.rev;
                        want.push_back({uint32_t(pp[t] >> 33), rev ? 1u : 0u,
                                        uint32_t((pp[t] >> 1) & 0xFFFFFFFFu),
                                        rev ? uint32_t(rd.size() - k - s.pos) : s.pos});
                    }
                });
                cbad += got.size() != want.size() || !std::equal(got.begin(), got.end(), want.begin(),
                    [](const auto& a, const auto& b) {
                        return a.rid == b.rid && a.rev == b.rev && a.rpos == b.rpos && a.qpos == b.qpos; });
            }
            std::printf("k=%u w=%-2u t=%u distinct=%zu occ=%zu bad=%zu fp=%zu collect_bad=%zu\n",
                        k, w, th, idx.distinct(), idx.occurrences(), bad, fp, cbad);
            fail += bad != 0 || fp != 0 || cbad != 0;
        }
    }
    std::printf(fail ? "seed_index_stress: %d FAILURE(S)\n" : "seed_index_stress: ALL PASS\n", fail);
    return fail != 0;
}
