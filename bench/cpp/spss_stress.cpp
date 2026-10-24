// spss_stress.cpp - SpssCore exactness against a hash-set oracle.
//
//   g++ -std=c++17 -O2 -march=native -I include bench/cpp/spss_stress.cpp && ./a.out [unitigs.fa]
//
// The random battery runs many small SPSS-like inputs (random sequences, repeats, low-complexity runs, N breaks) across k and m.
// Every input k-mer must be found in both orientations, random absent k-mers rejected, the streaming scan must agree with per-k-mer lookups, save/load must round-trip, and the dictionary (kmer_id / access) must be a dense bijection onto [0, n).
// With a FASTA argument the same checks run on that file, plus timings.
#include <chrono>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
#include "spss_core.hpp"
using namespace tpoptoa;

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("FAIL %s:%d ", __FILE__, __LINE__); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static std::vector<std::string> read_fa(const char* p) {
    std::vector<std::string> v; std::ifstream f(p); std::string l;
    while (std::getline(f, l)) { if (l.empty()) continue; if (l[0] == '>') v.emplace_back(); else if (!v.empty()) v.back() += l; }
    return v;
}

static void check(const std::vector<std::string>& seqs, unsigned k, unsigned m,
                  std::mt19937_64& rng, bool verbose) {
    SpssCore idx;
    auto t0 = std::chrono::steady_clock::now();
    if (!idx.build(seqs, k, m)) { CHECK(false, "build failed k=%u m=%u", k, m); return; }
    double bt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::unordered_set<uint64_t> truth;           // canonical k-mers
    std::vector<uint64_t> fw;                     // forward k-mers in order
    for (const auto& s : seqs) {
        uint64_t f = 0; std::size_t v = 0;
        for (char ch : s) {
            uint8_t b = base_encode(ch);
            if (b == 0xFF) { v = 0; f = 0; continue; }
            f = ((f << 2) | b) & kmer_mask(k);
            if (++v >= k) { fw.push_back(f); truth.insert(std::min(f, SpssCore::rc(f, k))); }
        }
    }
    std::size_t fn = 0, fp = 0;
    t0 = std::chrono::steady_clock::now();
    for (uint64_t f : fw) { fn += !idx.contains(f); fn += !idx.contains(SpssCore::rc(f, k)); }
    double lt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::size_t neg = 0, tries = 0;
    // A tiny k can have (almost) every k-mer present, so the attempts are bounded.
    while (neg < std::max<std::size_t>(fw.size(), 1000) && ++tries < 50 * (fw.size() + 1000)) {
        uint64_t y = rng() & kmer_mask(k);
        if (truth.count(std::min(y, SpssCore::rc(y, k)))) continue;
        fp += idx.contains(y); ++neg;
    }
    // streaming scan vs per-k-mer lookups on the sequences and their rc
    std::size_t scan_bad = 0;
    for (std::size_t i = 0; i < seqs.size() && i < 2000; ++i) {
        for (int o = 0; o < 2; ++o) {
            std::string s = seqs[i];
            if (o) { std::reverse(s.begin(), s.end()); for (auto& c : s) c = c=='A'?'T':c=='C'?'G':c=='G'?'C':c=='T'?'A':c; }
            SpssCore::ScanStats st; idx.scan(s.data(), s.size(), st);
            if (st.positive != st.valid) ++scan_bad;
        }
    }
    // Mixed scans, where mutated copies must count exactly the oracle positives, one record at a time and all at once through the pipelined scan_many
    std::vector<std::string> muts;
    SpssCore::ScanStats sum_one;
    for (std::size_t i = 0; i < seqs.size() && i < 500; ++i) {
        std::string s = seqs[i];
        for (auto& c : s) if (rng() % 20 == 0) c = "ACGTN"[rng() % 5];
        if (rng() & 1) { std::reverse(s.begin(), s.end()); for (auto& c : s) c = c=='A'?'T':c=='C'?'G':c=='G'?'C':c=='T'?'A':c; }
        uint64_t f = 0, want = 0, nv = 0; std::size_t v = 0;
        for (char ch : s) {
            uint8_t b = base_encode(ch);
            if (b == 0xFF) { v = 0; f = 0; continue; }
            f = ((f << 2) | b) & kmer_mask(k);
            if (++v >= k) { ++nv; want += truth.count(std::min(f, SpssCore::rc(f, k))) != 0; }
        }
        SpssCore::ScanStats st; idx.scan(s.data(), s.size(), st);
        if (st.positive != want || st.valid != nv) ++scan_bad;
        SpssCore::ScanStats ref; idx.scan_reference(s.data(), s.size(), ref);
        // records of <= kGroupMaxKmers windows are batch-looked-up with a search per k-mer, and longer ones stream exactly like the reference
        const bool tiny = s.size() >= k && s.size() - k + 1 <= SpssCore::kGroupMaxKmers;
        if (ref.positive != want || ref.valid != nv
            || (!tiny && (ref.searches != st.searches || ref.extensions != st.extensions))) ++scan_bad;
        sum_one.valid += st.valid; sum_one.positive += st.positive;
        muts.push_back(s);
        if (rng() % 3 == 0) {                               // single k-mer records
            for (std::size_t o = 0; o + k <= s.size() && o < 5 * k; o += k) muts.push_back(s.substr(o, k));
            for (std::size_t o = 0; o + k <= s.size() && o < 5 * k; o += k) {
                SpssCore::ScanStats t; idx.scan(muts[muts.size() - 1].data(), 0, t);
            }
        }
    }
    {
        std::vector<const char*> ps; std::vector<std::size_t> ls;
        SpssCore::ScanStats one;
        SpssCore::ScanStats ref;                           // single-cursor reference
        for (const auto& m2 : muts) {
            ps.push_back(m2.data()); ls.push_back(m2.size());
            idx.scan(m2.data(), m2.size(), one);
            if (m2.size() >= k && m2.size() - k + 1 <= SpssCore::kGroupMaxKmers) {
                SpssCore::ScanStats t; idx.scan_reference(m2.data(), m2.size(), t);
                ref.valid += t.valid; ref.positive += t.positive; ref.searches += t.valid;
            } else {
                idx.scan_reference(m2.data(), m2.size(), ref);
            }
        }
        // the interleaved lanes (scan_many) and the reference must agree on Every counter gives the same searches, extensions and answers per record.
        SpssCore::ScanStats all; idx.scan_many(ps.data(), ls.data(), ps.size(), all);
        if (all.valid != one.valid || all.positive != one.positive) ++scan_bad;
        if (all.valid != ref.valid || all.positive != ref.positive
            || all.searches != ref.searches || all.extensions != ref.extensions) ++scan_bad;
    }
    // Dictionary ids lie in [0, n), and access inverts lookup (in the reported orientation), absent k-mers have none, and with distinct input k-mers the ids are exactly {0, ..., n - 1}
    std::size_t dict_bad = 0;
    {
        const uint64_t n = idx.num_kmers();
        std::unordered_set<uint64_t> ids;
        for (uint64_t f : fw) {
            bool rev = false;
            const uint64_t id = idx.kmer_id(f, &rev);
            if (id == SpssCore::kNotFound || id >= n) { ++dict_bad; continue; }
            const uint64_t got = idx.access(id);
            if (got != (rev ? SpssCore::rc(f, k) : f)) ++dict_bad;
            ids.insert(id);
        }
        if (truth.size() == fw.size()) {               // every k-mer distinct
            if (n != fw.size() || ids.size() != n) ++dict_bad;
            for (uint64_t i = 0; i < n; ++i) {
                const uint64_t x = idx.access(i);
                if (idx.kmer_id(x) != i) { ++dict_bad; break; }
            }
        }
        for (int t = 0; t < 2000; ++t) {
            const uint64_t y = rng() & kmer_mask(k);
            if (!truth.count(std::min(y, SpssCore::rc(y, k))) && idx.kmer_id(y) != SpssCore::kNotFound) ++dict_bad;
        }
    }
    std::vector<uint8_t> blob; idx.save(blob);
    SpssCore ld; std::size_t used = ld.load(blob.data(), blob.size());
    std::size_t rt = used != blob.size();
    for (std::size_t i = 0; i < fw.size(); i += 7) rt += !ld.contains(fw[i]);
    // A borrowing load, the one query uses on a mapped file, answers the same with every array, the MPHF's and the Elias-Fano starts' included, read in place from the blob.
    SpssCore lb; rt += lb.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
    for (std::size_t i = 0; i < fw.size(); i += 7) rt += !lb.contains(fw[i]) || lb.kmer_id(fw[i]) != idx.kmer_id(fw[i]);
    for (int t = 0; t < 500; ++t) {
        const uint64_t y = rng() & kmer_mask(k);
        rt += lb.contains(y) != idx.contains(y);
    }
    CHECK(fn == 0 && fp == 0 && scan_bad == 0 && rt == 0 && dict_bad == 0,
          "k=%u m=%u FN=%zu FP=%zu scan_bad=%zu roundtrip=%zu dict_bad=%zu", k, m, fn, fp, scan_bad, rt, dict_bad);
    if (verbose) {
        const auto& S = idx.stats();
        std::printf("  k=%u m=%u kmers=%lu superkmers=%lu buckets=%lu ovf=%lu bits/kmer=%.2f build=%.2fs lookup=%.1fns FN=%zu FP=%zu\n",
                    k, m, (unsigned long)S.kmers, (unsigned long)S.superkmers, (unsigned long)S.buckets,
                    (unsigned long)S.ovf_buckets, idx.bytes() * 8.0 / S.kmers, bt, lt * 1e9 / (2 * fw.size()), fn, fp);
    }
}

int main(int argc, char** argv) {
    std::mt19937_64 rng(77);
    const char* A = "ACGT";
    for (int it = 0; it < 300; ++it) {
        unsigned k = 5 + rng() % 28, m = 2 + rng() % (k - 2);
        if (m > 31) m = 31;
        std::vector<std::string> seqs;
        int ns = 1 + rng() % 30;
        for (int i = 0; i < ns; ++i) {
            std::string s; std::size_t L = rng() % 300;
            for (std::size_t j = 0; j < L; ++j) {
                unsigned r = rng() % 100;
                s.push_back(r < 2 ? 'N' : r < 20 ? 'A' : A[rng() % 4]);
            }
            if (!seqs.empty() && rng() % 4 == 0) s += seqs[rng() % seqs.size()];   // repeats
            seqs.push_back(s);
        }
        check(seqs, k, m, rng, false);
    }
    std::printf("random battery done (%d failures)\n", g_fail);
    // Repeat families (human-like), with one unit copied hundreds to thousands of times with point mutations puts many super-k-mers behind a single minimizer, i.e. buckets beyond kSkewMin that the skew tables answer.
    // Exact duplicates of whole copies are included (repeated k-mers).
    uint64_t skew_seen = 0;
    for (int it = 0; it < 12; ++it) {
        const unsigned k = 31, m = (it % 3 == 0) ? 11 : (it % 3 == 1) ? 15 : 21;
        std::string unit;
        for (int j = 0; j < 200 + int(rng() % 200); ++j) unit.push_back(A[rng() % 4]);
        std::vector<std::string> seqs;
        const int copies = 300 + int(rng() % 2000);
        for (int c = 0; c < copies; ++c) {
            std::string s;
            for (int j = 0; j < int(rng() % 40); ++j) s.push_back(A[rng() % 4]);   // flank
            const unsigned div = 1 + rng() % 5;                                      // 1-5% diverged
            for (char ch : unit) s.push_back(rng() % 100 < div ? A[rng() % 4] : ch);
            for (int j = 0; j < int(rng() % 40); ++j) s.push_back(A[rng() % 4]);
            seqs.push_back(s);
            if (rng() % 10 == 0) seqs.push_back(s);                                   // exact duplicate
        }
        SpssCore probe_idx;
        if (probe_idx.build(seqs, k, m)) skew_seen += probe_idx.skew_entries();
        check(seqs, k, m, rng, false);
    }
    CHECK(skew_seen > 0, "repeat battery never reached the skew tables");
    std::printf("repeat-family battery done (%d failures, %llu k-mers via skew tables)\n",
                g_fail, (unsigned long long)skew_seen);
    if (argc > 1) {
        auto seqs = read_fa(argv[1]);
        for (unsigned m : {15u, 17u, 19u, 21u}) check(seqs, 31, m, rng, true);
    }
    std::printf(g_fail ? "spss_stress: %d FAILURE(S)\n" : "spss_stress: ALL PASS\n", g_fail);
    return g_fail != 0;
}
