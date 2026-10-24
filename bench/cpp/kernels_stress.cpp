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
static void expect(bool ok, const char* what, long& bad) { if (!ok) { ++bad; ++g_fail; (void)what; } }

static std::string rcs(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
    return r;
}
static std::string rnd(std::mt19937_64& g, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[g() & 3];
    return s;
}

// rc
static void test_rc() {
    std::mt19937_64 g(1);
    long bad = 0, tot = 0;
    for (std::size_t k = 1; k <= 32; ++k)
        for (int t = 0; t < 100000; ++t) {
            uint64_t x = g();
            if (t & 1) x &= kmer_mask(k);
            uint64_t ref = 0, y = x;
            for (std::size_t i = 0; i < k; ++i) { ref = (ref << 2) | ((y & 3) ^ 3); y >>= 2; }
            ++tot;
            expect(reverse_complement(x, k) == ref, "rc", bad);
        }
    std::printf("rc      %ld cases, %ld wrong\n", tot, bad);
}

// bpm
static bool isb(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't'; }
static int naive_ed(const std::string& q, const std::string& t) {
    const int m = static_cast<int>(q.size()), n = static_cast<int>(t.size());
    std::vector<int> p(m + 1), c(m + 1);
    for (int i = 0; i <= m; ++i) p[i] = i;
    int best = p[m];
    for (int j = 1; j <= n; ++j) {
        c[0] = 0;
        for (int i = 1; i <= m; ++i) {
            const bool eq = isb(q[i - 1]) && isb(t[j - 1]) && ((q[i - 1] | 32) == (t[j - 1] | 32));
            c[i] = std::min({p[i - 1] + (eq ? 0 : 1), p[i] + 1, c[i - 1] + 1});
        }
        best = std::min(best, c[m]);
        std::swap(p, c);
    }
    return best;
}
static void test_bpm() {
    std::mt19937_64 g(7);
    const char* al = "ACGTACGTACGTacgtNnX";
    long bad = 0, tot = 0;
    for (int it = 0; it < 60000; ++it) {
        const int m = 1 + static_cast<int>(g() % (it % 3 == 0 ? 300 : 160));
        std::string q(m, 'A');
        for (char& c : q) c = al[g() % (it % 5 == 0 ? 19 : 4)];
        const int cnt = 1 + static_cast<int>(g() % 4);
        std::string ts[4];
        const char* tp[4];
        int ns[4];
        for (int i = 0; i < cnt; ++i) {
            const int n = static_cast<int>(g() % 320);
            ts[i].assign(n, 'A');
            if (n > m && (g() & 1)) {                    // embedded near-match
                for (char& c : ts[i]) c = "ACGT"[g() & 3];
                const int off = static_cast<int>(g() % (n - m + 1));
                for (int x = 0; x < m; ++x) ts[i][off + x] = q[x];
                for (int e = static_cast<int>(g() % 8); e--;) ts[i][g() % n] = "ACGTN"[g() % 5];
            } else {
                for (char& c : ts[i]) c = al[g() % (it % 7 == 0 ? 19 : 4)];
            }
            tp[i] = ts[i].data();
            ns[i] = n;
        }
        const int cap = static_cast<int>(g() % 40);
        int out[4];
        bpm_semiglobal_ed_x4(q.data(), m, tp, ns, cnt, cap, out);
        for (int i = 0; i < cnt; ++i) {
            const int want = std::min(naive_ed(q, ts[i]), cap + 1);
            tot += 2;
            expect(bpm_semiglobal_ed(q.data(), m, tp[i], ns[i], cap) == want, "bpm", bad);
            expect(out[i] == want, "bpm_x4", bad);
        }
    }
    std::printf("bpm     %ld cases, %ld wrong\n", tot, bad);
}

// fm
static uint64_t pack(const std::string& s) {
    uint64_t v = 0;
    for (char c : s) v = (v << 2) | (c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3);
    return v;
}
static void test_fm() {
    std::mt19937_64 g(11);
    char fm_tmp[] = "/tmp/optik_kernels_fm_XXXXXX";
    { const int fd = mkstemp(fm_tmp); if (fd >= 0) close(fd); }
    long bad = 0, tot = 0;
    for (unsigned k : {1u, 2u, 5u, 16u, 31u, 32u}) {
        for (int rep = 0; rep < 6; ++rep) {
            std::vector<std::string> units;
            const std::string motif = rnd(g, 3 + g() % 40);
            const int nu = 1 + static_cast<int>(g() % 300);
            for (int u = 0; u < nu; ++u) {
                const int kind = static_cast<int>(g() % 6);
                std::string s;
                if (kind == 0) s = rnd(g, g() % (k + 1));                  // may be < k
                else if (kind == 1) { while (s.size() < k + 50) s += motif; }   // repeats
                else if (kind == 2) s = std::string(k + g() % 60, "ACGT"[g() & 3]);
                else s = rnd(g, k + g() % (rep == 5 ? 5000 : 200));
                units.push_back(s);
            }
            std::unordered_set<uint64_t> truth;
            for (const auto& u : units)
                for (std::size_t i = 0; i + k <= u.size(); ++i) {
                    const std::string w = u.substr(i, k);
                    truth.insert(pack(w));
                    truth.insert(pack(rcs(w)));
                }
            FmIndex fm(k, units, truth.size(), 4);
            fm.save(fm_tmp);
            FmIndex* ld = FmIndex::load(fm_tmp);
            ++tot;
            expect(ld != nullptr && ld->bwt_len() == fm.bwt_len(), "fm load", bad);
            for (const auto& x : truth) {
                tot += 2;
                expect(fm.find(x), "fm+", bad);
                expect(ld && ld->find(x), "fm+ loaded", bad);
            }
            for (int t = 0; t < 20000; ++t) {
                const uint64_t x = g() & kmer_mask(k);
                tot += 2;
                expect(fm.find(x) == (truth.count(x) > 0), "fm?", bad);
                expect(ld && ld->find(x) == (truth.count(x) > 0), "fm? loaded", bad);
            }
            // find_batch() over present and random k-mers mixed, in a count that is not a multiple of the batch width, must give find()'s answers.
            std::vector<uint64_t> mix(truth.begin(), truth.end());
            for (int t = 0; t < 5003; ++t) mix.push_back(g() & kmer_mask(k));
            std::shuffle(mix.begin(), mix.end(), g);
            std::vector<uint8_t> hit(mix.size(), 2), hit_ld(mix.size(), 2);
            fm.find_batch(mix.data(), mix.size(), hit.data());
            if (ld) ld->find_batch(mix.data(), mix.size(), hit_ld.data());
            for (std::size_t i = 0; i < mix.size(); ++i) {
                tot += 2;
                const uint8_t want = truth.count(mix[i]) > 0;
                expect(hit[i] == want, "fm batch", bad);
                expect(hit_ld[i] == want, "fm batch loaded", bad);
            }
            delete ld;
        }
    }
    // BWT lengths at and around multiples of the rank block size, 448 bits, and of 512, the earlier superblock size.
    // Backward search starts with rank(n), which once read one superblock entry past the end when n was a multiple of 512, giving wrong answers or a crash.
    for (unsigned k : {16u, 31u})
        for (std::size_t n = 440; n <= 4100; ++n) {
            const std::size_t b = FmRankBits::kBitsPerBlock;
            if ((n % 512 > 2 && n % 512 < 510) && (n % b > 2 && n % b < b - 2)) continue;
            const std::vector<std::string> one{rnd(g, n - 2)};       // n = len + X + $
            FmIndex fm(k, one, 0, 1);
            ++tot;
            expect(fm.bwt_len() == n, "fm len", bad);
            std::vector<uint64_t> win;
            for (std::size_t i = 0; i + k <= one[0].size(); i += 7) {
                ++tot;
                win.push_back(pack(one[0].substr(i, k)));
                expect(fm.find(win.back()), "fm boundary", bad);
            }
            std::vector<uint8_t> hit(win.size(), 0);
            fm.find_batch(win.data(), win.size(), hit.data());
            for (uint8_t h : hit) { ++tot; expect(h == 1, "fm batch boundary", bad); }
        }
    std::remove(fm_tmp);
    std::printf("fm      %ld cases, %ld wrong\n", tot, bad);
}

// spss
static void test_spss() {
    std::mt19937_64 g(23);
    long bad = 0, tot = 0;
    for (int rep = 0; rep < 24; ++rep) {
        const unsigned k = rep % 3 == 0 ? 31 : rep % 3 == 1 ? 21 : 15;
        const unsigned m = k / 2;
        std::vector<std::string> seqs;
        const std::string shared = rnd(g, 400);
        for (int s = 0; s < 40; ++s) {
            std::string x;
            for (int part = 0; part < 6; ++part) {
                switch (g() % 5) {
                    case 0: x += rnd(g, 50 + g() % 300); break;
                    case 1: { const std::string y = rnd(g, k / 2 + 3 + g() % 10); x += y + rcs(y); break; }   // palindrome
                    case 2: x += shared.substr(g() % 200, 100 + g() % 100); break;                    // duplicates
                    case 3: x += rcs(shared.substr(g() % 200, 100 + g() % 100)); break;               // rc duplicates
                    default: x += std::string(40 + g() % 200, "ACGT"[g() & 3]); break;                // homopolymer
                }
            }
            seqs.push_back(x);
        }
        // brute distinct canonical k-mers
        std::map<std::string, int> canon;
        for (const auto& s : seqs)
            for (std::size_t i = 0; i + k <= s.size(); ++i) {
                const std::string w = s.substr(i, k), r = rcs(w);
                canon[std::min(w, r)] = 0;
            }
        SpssCore idx;
        if (!idx.build(seqs, k, m, 4)) { ++bad; ++g_fail; continue; }
        ++tot;
        expect(idx.count_distinct() == canon.size(), "distinct(build)", bad);
        std::vector<uint8_t> blob;
        idx.save(blob);
        SpssCore ld;
        ++tot;
        expect(ld.load(blob.data(), blob.size()) == blob.size() && ld.count_distinct() == canon.size(),
               "distinct(load)", bad);
        // one strand-independent id per distinct k-mer
        std::map<std::string, std::set<uint64_t>> ids;
        for (const auto& s : seqs)
            for (std::size_t i = 0; i + k <= s.size(); ++i) {
                const std::string w = s.substr(i, k), r = rcs(w);
                ids[std::min(w, r)].insert(idx.kmer_id(pack(w)));
                ids[std::min(w, r)].insert(idx.kmer_id(pack(r)));
            }
        std::set<uint64_t> all;
        for (const auto& kv : ids) {
            ++tot;
            expect(kv.second.size() == 1 && *kv.second.begin() != SpssCore::kNotFound, "one id", bad);
            all.insert(*kv.second.begin());
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
