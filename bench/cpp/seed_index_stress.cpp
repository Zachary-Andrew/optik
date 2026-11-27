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
    const char* bases = "ACGT";
    std::vector<std::string> seqs;
    for (std::size_t seq_len : {2500000ul, 1048600ul, 5000ul, 20ul}) {
        std::string seq(seq_len, 'A');
        for (auto& ch : seq) ch = bases[rng() & 3];
        for (int r = 0; r < 40 && seq_len > 10000; ++r) {          // repeats and N runs
            std::size_t src_pos = rng() % (seq_len - 3000), dst_pos = rng() % (seq_len - 3000);
            seq.replace(dst_pos, 2000, seq.substr(src_pos, 2000));
            std::size_t n_start = rng() % (seq_len - 50);
            for (int t = 0; t < int(rng() % 40); ++t) seq[n_start + t] = 'N';
        }
        seqs.push_back(seq);
    }
    std::vector<std::pair<const char*, std::size_t>> refs;
    for (auto& seq : seqs) refs.push_back({seq.data(), seq.size()});
    int fail = 0;
    for (unsigned kmer_len : {15u, 21u}) for (unsigned window : {5u, 10u, 50u}) {
        std::map<uint64_t, std::vector<uint64_t>> oracle;
        for (std::size_t rid = 0; rid < refs.size(); ++rid)
            extract_minimizers(refs[rid].first, refs[rid].second, kmer_len, window, [&](const MinimizerSeed& seed) {
                oracle[seed.kmer].push_back((uint64_t(rid) << 33) | (uint64_t(seed.pos) << 1) | seed.rev);
            });
        std::size_t occ = 0;
        for (auto& [key, posts] : oracle) { std::sort(posts.begin(), posts.end()); occ += posts.size(); }
        for (unsigned num_threads : {1u, 3u, 8u}) {
            SeedIndex idx;
            if (!idx.build(refs, kmer_len, window, 1u << 30, num_threads)) { std::printf("build failed\n"); ++fail; continue; }
            std::size_t bad = idx.distinct() != oracle.size() || idx.occurrences() != occ;
            for (auto& [key, posts] : oracle) {
                auto [posts_ptr, count] = idx.find(key);
                bad += !posts_ptr || count != posts.size() || !std::equal(posts.begin(), posts.end(), posts_ptr);
            }
            std::size_t false_pos = 0;
            for (int i = 0; i < 20000; ++i) {
                uint64_t kmer = rng() & kmer_mask(kmer_len), rev_comp = reverse_complement(kmer, kmer_len);
                kmer = std::min(kmer, rev_comp);
                false_pos += !oracle.count(kmer) && idx.find(kmer).first;
            }
            // pipelined collect() == naive per-minimizer find() in read order
            std::size_t cbad = 0;
            for (int i = 0; i < 300; ++i) {
                const auto& src = seqs[rng() % 2];
                std::size_t read_len = 50 + rng() % 400, read_start = rng() % (src.size() - read_len);
                std::string read = src.substr(read_start, read_len);
                for (auto& ch : read) if (rng() % 50 == 0) ch = bases[rng() & 3];
                std::vector<SeedIndex::Hit> got, want;
                idx.collect(read.data(), read.size(), got);
                extract_minimizers(read.data(), read.size(), kmer_len, window, [&](const MinimizerSeed& seed) {
                    auto [posts_ptr, count] = idx.find(seed.kmer);
                    if (!posts_ptr) return;
                    for (uint32_t t = 0; t < count; ++t) {
                        const bool rev = ((posts_ptr[t] & 1) != 0) != seed.rev;
                        want.push_back({uint32_t(posts_ptr[t] >> 33), rev ? 1u : 0u,
                                        uint32_t((posts_ptr[t] >> 1) & 0xFFFFFFFFu),
                                        rev ? uint32_t(read.size() - kmer_len - seed.pos) : seed.pos});
                    }
                });
                cbad += got.size() != want.size() || !std::equal(got.begin(), got.end(), want.begin(),
                    [](const auto& lhs, const auto& rhs) {
                        return lhs.rid == rhs.rid && lhs.rev == rhs.rev && lhs.rpos == rhs.rpos && lhs.qpos == rhs.qpos; });
            }
            std::printf("k=%u w=%-2u t=%u distinct=%zu occ=%zu bad=%zu fp=%zu collect_bad=%zu\n",
                        kmer_len, window, num_threads, idx.distinct(), idx.occurrences(), bad, false_pos, cbad);
            fail += bad != 0 || false_pos != 0 || cbad != 0;
        }
    }
    std::printf(fail ? "seed_index_stress: %d FAILURE(S)\n" : "seed_index_stress: ALL PASS\n", fail);
    return fail != 0;
}
