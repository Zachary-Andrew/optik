// kernels_stress.cpp - exactness of the scaling kernels against brute force.
//
//   g++ -std=c++17 -O2 -march=native -fopenmp -I include bench/cpp/kernels_stress.cpp && ./a.out
//
//   rc       constant-time reverse_complement vs a base-by-base loop, all k in 1..32, masked and unmasked words.
//   bpm      bit-parallel semi-global edit distance, scalar and four-lane AVX2, vs a naive DP over query lengths 1..300, i.e. 1..5 words, N, lowercase and junk bytes, embedded near-matches, 1..4 lanes of unequal length, empty texts, every cap.
//   fm       k-order-BWT FM index vs the brute-force window set, for k in {1, 2, 5, 16, 31, 32}, repeats, tiny and long unitigs, unitigs shorter than k.
//            It covers all windows of the text on both strands and random k-mers, in memory and after save/load, and BWT lengths at each multiple of 512, the rank superblock, up to 4 kb.
//   spss     the per-bucket distinct count vs brute force on inputs with duplicate k-mers and reverse-complement palindromes, before and after save/load, where the load path recounts by lookups.
//            kmer_id must be strand-independent and one id per distinct k-mer.
#include <algorithm>
#include <cstdio>
#include <map>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include <unistd.h>

#include "bpm_edit.hpp"
#include "fm_index.hpp"
#include "kmer.hpp"
#include "spss_core.hpp"

using namespace tpoptoa;

static long g_fail = 0;
static void expect(bool passed, const char* what, long& bad) { if (!passed) { ++bad; ++g_fail; (void)what; } }

static std::string rcs(const std::string& seq) {
    std::string rev(seq.rbegin(), seq.rend());
    for (char& ch : rev) ch = ch == 'A' ? 'T' : ch == 'C' ? 'G' : ch == 'G' ? 'C' : ch == 'T' ? 'A' : ch;
    return rev;
}
static std::string rnd(std::mt19937_64& rng, std::size_t len) {
    std::string seq(len, 'A');
    for (char& ch : seq) ch = "ACGT"[rng() & 3];
    return seq;
}

// rc
static void test_rc() {
    std::mt19937_64 rng(1);
    long bad = 0, tot = 0;
    for (std::size_t kmer_len = 1; kmer_len <= 32; ++kmer_len)
        for (int t = 0; t < 100000; ++t) {
            uint64_t kmer = rng();
            if (t & 1) kmer &= kmer_mask(kmer_len);
            uint64_t ref = 0, rest = kmer;
            for (std::size_t i = 0; i < kmer_len; ++i) { ref = (ref << 2) | ((rest & 3) ^ 3); rest >>= 2; }
            ++tot;
            expect(reverse_complement(kmer, kmer_len) == ref, "rc", bad);
        }
    std::printf("rc      %ld cases, %ld wrong\n", tot, bad);
}

// bpm
static bool isb(char ch) { return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't'; }
static int naive_ed(const std::string& query, const std::string& text) {
    const int query_len = static_cast<int>(query.size()), text_len = static_cast<int>(text.size());
    std::vector<int> prev_col(query_len + 1), cur_col(query_len + 1);
    for (int i = 0; i <= query_len; ++i) prev_col[i] = i;
    int best = prev_col[query_len];
    for (int j = 1; j <= text_len; ++j) {
        cur_col[0] = 0;
        for (int i = 1; i <= query_len; ++i) {
            const bool same = isb(query[i - 1]) && isb(text[j - 1]) && ((query[i - 1] | 32) == (text[j - 1] | 32));
            cur_col[i] = std::min({prev_col[i - 1] + (same ? 0 : 1), prev_col[i] + 1, cur_col[i - 1] + 1});
        }
        best = std::min(best, cur_col[query_len]);
        std::swap(prev_col, cur_col);
    }
    return best;
}
static void test_bpm() {
    std::mt19937_64 rng(7);
    const char* alphabet = "ACGTACGTACGTacgtNnX";
    long bad = 0, tot = 0;
    for (int it = 0; it < 60000; ++it) {
        const int query_len = 1 + static_cast<int>(rng() % (it % 3 == 0 ? 300 : 160));
        std::string query(query_len, 'A');
        for (char& ch : query) ch = alphabet[rng() % (it % 5 == 0 ? 19 : 4)];
        const int cnt = 1 + static_cast<int>(rng() % 4);
        std::string texts[4];
        const char* text_ptrs[4];
        int text_lens[4];
        for (int i = 0; i < cnt; ++i) {
            const int text_len = static_cast<int>(rng() % 320);
            texts[i].assign(text_len, 'A');
            if (text_len > query_len && (rng() & 1)) {                    // embedded near-match
                for (char& ch : texts[i]) ch = "ACGT"[rng() & 3];
                const int off = static_cast<int>(rng() % (text_len - query_len + 1));
                for (int x = 0; x < query_len; ++x) texts[i][off + x] = query[x];
                for (int e = static_cast<int>(rng() % 8); e--;) texts[i][rng() % text_len] = "ACGTN"[rng() % 5];
            } else {
                for (char& ch : texts[i]) ch = alphabet[rng() % (it % 7 == 0 ? 19 : 4)];
            }
            text_ptrs[i] = texts[i].data();
            text_lens[i] = text_len;
        }
        const int cap = static_cast<int>(rng() % 40);
        int out[4];
        bpm_semiglobal_ed_x4(query.data(), query_len, text_ptrs, text_lens, cnt, cap, out);
        for (int i = 0; i < cnt; ++i) {
            const int want = std::min(naive_ed(query, texts[i]), cap + 1);
            tot += 2;
            expect(bpm_semiglobal_ed(query.data(), query_len, text_ptrs[i], text_lens[i], cap) == want, "bpm", bad);
            expect(out[i] == want, "bpm_x4", bad);
        }
    }
    std::printf("bpm     %ld cases, %ld wrong\n", tot, bad);
}

// fm
static uint64_t pack(const std::string& seq) {
    uint64_t packed = 0;
    for (char ch : seq) packed = (packed << 2) | (ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2 : 3);
    return packed;
}
static void test_fm() {
    std::mt19937_64 rng(11);
    char fm_tmp[] = "/tmp/optik_kernels_fm_XXXXXX";
    { const int file_desc = mkstemp(fm_tmp); if (file_desc >= 0) close(file_desc); }
    long bad = 0, tot = 0;
    for (unsigned kmer_len : {1u, 2u, 5u, 16u, 31u, 32u}) {
        for (int rep = 0; rep < 6; ++rep) {
            std::vector<std::string> units;
            const std::string motif = rnd(rng, 3 + rng() % 40);
            const int num_units = 1 + static_cast<int>(rng() % 300);
            for (int u = 0; u < num_units; ++u) {
                const int kind = static_cast<int>(rng() % 6);
                std::string unit;
                if (kind == 0) unit = rnd(rng, rng() % (kmer_len + 1));                  // may be < k
                else if (kind == 1) { while (unit.size() < kmer_len + 50) unit += motif; }   // repeats
                else if (kind == 2) unit = std::string(kmer_len + rng() % 60, "ACGT"[rng() & 3]);
                else unit = rnd(rng, kmer_len + rng() % (rep == 5 ? 5000 : 200));
                units.push_back(unit);
            }
            std::unordered_set<uint64_t> truth;
            for (const auto& unit : units)
                for (std::size_t i = 0; i + kmer_len <= unit.size(); ++i) {
                    const std::string window = unit.substr(i, kmer_len);
                    truth.insert(pack(window));
                    truth.insert(pack(rcs(window)));
                }
            FmIndex fm_idx(kmer_len, units, truth.size(), 4);
            fm_idx.save(fm_tmp);
            FmIndex* loaded = FmIndex::load(fm_tmp);
            ++tot;
            expect(loaded != nullptr && loaded->bwt_len() == fm_idx.bwt_len(), "fm load", bad);
            for (const auto& kmer : truth) {
                tot += 2;
                expect(fm_idx.find(kmer), "fm+", bad);
                expect(loaded && loaded->find(kmer), "fm+ loaded", bad);
            }
            for (int t = 0; t < 20000; ++t) {
                const uint64_t kmer = rng() & kmer_mask(kmer_len);
                tot += 2;
                expect(fm_idx.find(kmer) == (truth.count(kmer) > 0), "fm?", bad);
                expect(loaded && loaded->find(kmer) == (truth.count(kmer) > 0), "fm? loaded", bad);
            }
            // find_batch() over present and random k-mers mixed, in a count that is not a multiple of the batch width, must give find()'s answers.
            std::vector<uint64_t> mix(truth.begin(), truth.end());
            for (int t = 0; t < 5003; ++t) mix.push_back(rng() & kmer_mask(kmer_len));
            std::shuffle(mix.begin(), mix.end(), rng);
            std::vector<uint8_t> hit(mix.size(), 2), hit_ld(mix.size(), 2);
            fm_idx.find_batch(mix.data(), mix.size(), hit.data());
            if (loaded) loaded->find_batch(mix.data(), mix.size(), hit_ld.data());
            for (std::size_t i = 0; i < mix.size(); ++i) {
                tot += 2;
                const uint8_t want = truth.count(mix[i]) > 0;
                expect(hit[i] == want, "fm batch", bad);
                expect(hit_ld[i] == want, "fm batch loaded", bad);
            }
            delete loaded;
        }
    }
    // BWT lengths at and around multiples of the rank block size, 448 bits, and of 512, the earlier superblock size.
    // Backward search starts with rank(n), which once read one superblock entry past the end when n was a multiple of 512, giving wrong answers or a crash.
    for (unsigned kmer_len : {16u, 31u})
        for (std::size_t n = 440; n <= 4100; ++n) {
            const std::size_t block_bits = FmRankBits::kBitsPerBlock;
            if ((n % 512 > 2 && n % 512 < 510) && (n % block_bits > 2 && n % block_bits < block_bits - 2)) continue;
            const std::vector<std::string> one{rnd(rng, n - 2)};       // n = len + X + $
            FmIndex fm_idx(kmer_len, one, 0, 1);
            ++tot;
            expect(fm_idx.bwt_len() == n, "fm len", bad);
            std::vector<uint64_t> win;
            for (std::size_t i = 0; i + kmer_len <= one[0].size(); i += 7) {
                ++tot;
                win.push_back(pack(one[0].substr(i, kmer_len)));
                expect(fm_idx.find(win.back()), "fm boundary", bad);
            }
            std::vector<uint8_t> hit(win.size(), 0);
            fm_idx.find_batch(win.data(), win.size(), hit.data());
            for (uint8_t hit_val : hit) { ++tot; expect(hit_val == 1, "fm batch boundary", bad); }
        }
    std::remove(fm_tmp);
    std::printf("fm      %ld cases, %ld wrong\n", tot, bad);
}

// spss
static void test_spss() {
    std::mt19937_64 rng(23);
    long bad = 0, tot = 0;
    for (int rep = 0; rep < 24; ++rep) {
        const unsigned kmer_len = rep % 3 == 0 ? 31 : rep % 3 == 1 ? 21 : 15;
        const unsigned minimizer_len = kmer_len / 2;
        std::vector<std::string> seqs;
        const std::string shared = rnd(rng, 400);
        for (int s = 0; s < 40; ++s) {
            std::string seq;
            for (int part = 0; part < 6; ++part) {
                switch (rng() % 5) {
                    case 0: seq += rnd(rng, 50 + rng() % 300); break;
                    case 1: { const std::string half = rnd(rng, kmer_len / 2 + 3 + rng() % 10); seq += half + rcs(half); break; }   // palindrome
                    case 2: seq += shared.substr(rng() % 200, 100 + rng() % 100); break;                    // duplicates
                    case 3: seq += rcs(shared.substr(rng() % 200, 100 + rng() % 100)); break;               // rc duplicates
                    default: seq += std::string(40 + rng() % 200, "ACGT"[rng() & 3]); break;                // homopolymer
                }
            }
            seqs.push_back(seq);
        }
        // brute distinct canonical k-mers
        std::map<std::string, int> canon;
        for (const auto& seq : seqs)
            for (std::size_t i = 0; i + kmer_len <= seq.size(); ++i) {
                const std::string window = seq.substr(i, kmer_len), rev = rcs(window);
                canon[std::min(window, rev)] = 0;
            }
        SpssCore idx;
        if (!idx.build(seqs, kmer_len, minimizer_len, 4)) { ++bad; ++g_fail; continue; }
        ++tot;
        expect(idx.count_distinct() == canon.size(), "distinct(build)", bad);
        std::vector<uint8_t> blob;
        idx.save(blob);
        SpssCore loaded;
        ++tot;
        expect(loaded.load(blob.data(), blob.size()) == blob.size() && loaded.count_distinct() == canon.size(),
               "distinct(load)", bad);
        // one strand-independent id per distinct k-mer
        std::map<std::string, std::set<uint64_t>> ids;
        for (const auto& seq : seqs)
            for (std::size_t i = 0; i + kmer_len <= seq.size(); ++i) {
                const std::string window = seq.substr(i, kmer_len), rev = rcs(window);
                ids[std::min(window, rev)].insert(idx.kmer_id(pack(window)));
                ids[std::min(window, rev)].insert(idx.kmer_id(pack(rev)));
            }
        std::set<uint64_t> all;
        for (const auto& entry : ids) {
            ++tot;
            expect(entry.second.size() == 1 && *entry.second.begin() != SpssCore::kNotFound, "one id", bad);
            all.insert(*entry.second.begin());
        }
        ++tot;
        expect(all.size() == ids.size(), "ids distinct", bad);
    }
    std::printf("spss    %ld cases, %ld wrong\n", tot, bad);
}

int main() {
    test_rc();
    test_bpm();
    test_fm();
    test_spss();
    std::printf("%s\n", g_fail ? "KERNELS STRESS FAILED" : "kernels_stress: ALL PASS");
    return g_fail ? 1 : 0;
}
