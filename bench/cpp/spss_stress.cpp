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

static std::vector<std::string> read_fa(const char* path) {
    std::vector<std::string> seqs; std::ifstream file(path); std::string line;
    while (std::getline(file, line)) { if (line.empty()) continue; if (line[0] == '>') seqs.emplace_back(); else if (!seqs.empty()) seqs.back() += line; }
    return seqs;
}

static void check(const std::vector<std::string>& seqs, unsigned kmer_len, unsigned minimizer_len,
                  std::mt19937_64& rng, bool verbose) {
    SpssCore idx;
    auto start = std::chrono::steady_clock::now();
    if (!idx.build(seqs, kmer_len, minimizer_len)) { CHECK(false, "build failed k=%u m=%u", kmer_len, minimizer_len); return; }
    double build_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::unordered_set<uint64_t> truth;           // canonical k-mers
    std::vector<uint64_t> fwd_kmers;                     // forward k-mers in order
    for (const auto& seq : seqs) {
        uint64_t fwd = 0; std::size_t run_len = 0;
        for (char base_char : seq) {
            uint8_t base = base_encode(base_char);
            if (base == 0xFF) { run_len = 0; fwd = 0; continue; }
            fwd = ((fwd << 2) | base) & kmer_mask(kmer_len);
            if (++run_len >= kmer_len) { fwd_kmers.push_back(fwd); truth.insert(std::min(fwd, SpssCore::revcomp_packed(fwd, kmer_len))); }
        }
    }
    std::size_t false_neg = 0, false_pos = 0;
    start = std::chrono::steady_clock::now();
    for (uint64_t kmer : fwd_kmers) { false_neg += !idx.contains(kmer); false_neg += !idx.contains(SpssCore::revcomp_packed(kmer, kmer_len)); }
    double lookup_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::size_t neg = 0, tries = 0;
    // A tiny k can have (almost) every k-mer present, so the attempts are bounded.
    while (neg < std::max<std::size_t>(fwd_kmers.size(), 1000) && ++tries < 50 * (fwd_kmers.size() + 1000)) {
        uint64_t kmer = rng() & kmer_mask(kmer_len);
        if (truth.count(std::min(kmer, SpssCore::revcomp_packed(kmer, kmer_len)))) continue;
        false_pos += idx.contains(kmer); ++neg;
    }
    // streaming scan vs per-k-mer lookups on the sequences and their rc
    std::size_t scan_bad = 0;
    for (std::size_t i = 0; i < seqs.size() && i < 2000; ++i) {
        for (int o = 0; o < 2; ++o) {
            std::string seq = seqs[i];
            if (o) { std::reverse(seq.begin(), seq.end()); for (auto& ch : seq) ch = ch=='A'?'T':ch=='C'?'G':ch=='G'?'C':ch=='T'?'A':ch; }
            SpssCore::ScanStats stats; idx.scan(seq.data(), seq.size(), stats);
            if (stats.positive != stats.valid) ++scan_bad;
        }
    }
    // Mixed scans, where mutated copies must count exactly the oracle positives, one record at a time and all at once through the pipelined scan_many
    std::vector<std::string> muts;
    SpssCore::ScanStats sum_one;
    for (std::size_t i = 0; i < seqs.size() && i < 500; ++i) {
        std::string seq = seqs[i];
        for (auto& ch : seq) if (rng() % 20 == 0) ch = "ACGTN"[rng() % 5];
        if (rng() & 1) { std::reverse(seq.begin(), seq.end()); for (auto& ch : seq) ch = ch=='A'?'T':ch=='C'?'G':ch=='G'?'C':ch=='T'?'A':ch; }
        uint64_t fwd = 0, want = 0, num_valid = 0; std::size_t run_len = 0;
        for (char base_char : seq) {
            uint8_t base = base_encode(base_char);
            if (base == 0xFF) { run_len = 0; fwd = 0; continue; }
            fwd = ((fwd << 2) | base) & kmer_mask(kmer_len);
            if (++run_len >= kmer_len) { ++num_valid; want += truth.count(std::min(fwd, SpssCore::revcomp_packed(fwd, kmer_len))) != 0; }
        }
        SpssCore::ScanStats stats; idx.scan(seq.data(), seq.size(), stats);
        if (stats.positive != want || stats.valid != num_valid) ++scan_bad;
        SpssCore::ScanStats ref; idx.scan_reference(seq.data(), seq.size(), ref);
        // records of <= kGroupMaxKmers windows are batch-looked-up with a search per k-mer, and longer ones stream exactly like the reference
        const bool tiny = seq.size() >= kmer_len && seq.size() - kmer_len + 1 <= SpssCore::kGroupMaxKmers;
        if (ref.positive != want || ref.valid != num_valid
            || (!tiny && (ref.searches != stats.searches || ref.extensions != stats.extensions))) ++scan_bad;
        sum_one.valid += stats.valid; sum_one.positive += stats.positive;
        muts.push_back(seq);
        if (rng() % 3 == 0) {                               // single k-mer records
            for (std::size_t o = 0; o + kmer_len <= seq.size() && o < 5 * kmer_len; o += kmer_len) muts.push_back(seq.substr(o, kmer_len));
            for (std::size_t o = 0; o + kmer_len <= seq.size() && o < 5 * kmer_len; o += kmer_len) {
                SpssCore::ScanStats stats_tmp; idx.scan(muts[muts.size() - 1].data(), 0, stats_tmp);
            }
        }
    }
    {
        std::vector<const char*> ptrs; std::vector<std::size_t> lens;
        SpssCore::ScanStats one;
        SpssCore::ScanStats ref;                           // single-cursor reference
        for (const auto& mutated : muts) {
            ptrs.push_back(mutated.data()); lens.push_back(mutated.size());
            idx.scan(mutated.data(), mutated.size(), one);
            if (mutated.size() >= kmer_len && mutated.size() - kmer_len + 1 <= SpssCore::kGroupMaxKmers) {
                SpssCore::ScanStats ref_one; idx.scan_reference(mutated.data(), mutated.size(), ref_one);
                ref.valid += ref_one.valid; ref.positive += ref_one.positive; ref.searches += ref_one.valid;
            } else {
                idx.scan_reference(mutated.data(), mutated.size(), ref);
            }
        }
        // the interleaved lanes (scan_many) and the reference must agree on Every counter gives the same searches, extensions and answers per record.
        SpssCore::ScanStats all; idx.scan_many(ptrs.data(), lens.data(), ptrs.size(), all);
        if (all.valid != one.valid || all.positive != one.positive) ++scan_bad;
        if (all.valid != ref.valid || all.positive != ref.positive
            || all.searches != ref.searches || all.extensions != ref.extensions) ++scan_bad;
    }
    // Dictionary ids lie in [0, n), and access inverts lookup (in the reported orientation), absent k-mers have none, and with distinct input k-mers the ids are exactly {0, ..., n - 1}
    std::size_t dict_bad = 0;
    {
        const uint64_t num_kmers = idx.num_kmers();
        std::unordered_set<uint64_t> ids;
        for (uint64_t kmer : fwd_kmers) {
            bool rev = false;
            const uint64_t kmer_ident = idx.kmer_id(kmer, &rev);
            if (kmer_ident == SpssCore::kNotFound || kmer_ident >= num_kmers) { ++dict_bad; continue; }
            const uint64_t got = idx.access(kmer_ident);
            if (got != (rev ? SpssCore::revcomp_packed(kmer, kmer_len) : kmer)) ++dict_bad;
            ids.insert(kmer_ident);
        }
        if (truth.size() == fwd_kmers.size()) {               // every k-mer distinct
            if (num_kmers != fwd_kmers.size() || ids.size() != num_kmers) ++dict_bad;
            for (uint64_t i = 0; i < num_kmers; ++i) {
                const uint64_t kmer = idx.access(i);
                if (idx.kmer_id(kmer) != i) { ++dict_bad; break; }
            }
        }
        for (int t = 0; t < 2000; ++t) {
            const uint64_t kmer = rng() & kmer_mask(kmer_len);
            if (!truth.count(std::min(kmer, SpssCore::revcomp_packed(kmer, kmer_len))) && idx.kmer_id(kmer) != SpssCore::kNotFound) ++dict_bad;
        }
    }
    std::vector<uint8_t> blob; idx.save(blob);
    SpssCore loaded; std::size_t used = loaded.load(blob.data(), blob.size());
    std::size_t roundtrip_bad = used != blob.size();
    for (std::size_t i = 0; i < fwd_kmers.size(); i += 7) roundtrip_bad += !loaded.contains(fwd_kmers[i]);
    // A borrowing load, the one query uses on a mapped file, answers the same with every array, the MPHF's and the Elias-Fano starts' included, read in place from the blob.
    SpssCore borrowed; roundtrip_bad += borrowed.load(blob.data(), blob.size(), /*borrow=*/true) != blob.size();
    for (std::size_t i = 0; i < fwd_kmers.size(); i += 7) roundtrip_bad += !borrowed.contains(fwd_kmers[i]) || borrowed.kmer_id(fwd_kmers[i]) != idx.kmer_id(fwd_kmers[i]);
    for (int t = 0; t < 500; ++t) {
        const uint64_t kmer = rng() & kmer_mask(kmer_len);
        roundtrip_bad += borrowed.contains(kmer) != idx.contains(kmer);
    }
    CHECK(false_neg == 0 && false_pos == 0 && scan_bad == 0 && roundtrip_bad == 0 && dict_bad == 0,
          "k=%u m=%u FN=%zu FP=%zu scan_bad=%zu roundtrip=%zu dict_bad=%zu", kmer_len, minimizer_len, false_neg, false_pos, scan_bad, roundtrip_bad, dict_bad);
    if (verbose) {
        const auto& stats_ref = idx.stats();
        std::printf("  k=%u m=%u kmers=%lu superkmers=%lu buckets=%lu ovf=%lu bits/kmer=%.2f build=%.2fs lookup=%.1fns FN=%zu FP=%zu\n",
                    kmer_len, minimizer_len, (unsigned long)stats_ref.kmers, (unsigned long)stats_ref.superkmers, (unsigned long)stats_ref.buckets,
                    (unsigned long)stats_ref.ovf_buckets, idx.bytes() * 8.0 / stats_ref.kmers, build_secs, lookup_secs * 1e9 / (2 * fwd_kmers.size()), false_neg, false_pos);
    }
}

int main(int argc, char** argv) {
    std::mt19937_64 rng(77);
    const char* bases = "ACGT";
    for (int it = 0; it < 300; ++it) {
        unsigned kmer_len = 5 + rng() % 28, minimizer_len = 2 + rng() % (kmer_len - 2);
        if (minimizer_len > 31) minimizer_len = 31;
        std::vector<std::string> seqs;
        int num_seqs = 1 + rng() % 30;
        for (int i = 0; i < num_seqs; ++i) {
            std::string seq; std::size_t len = rng() % 300;
            for (std::size_t j = 0; j < len; ++j) {
                unsigned roll = rng() % 100;
                seq.push_back(roll < 2 ? 'N' : roll < 20 ? 'A' : bases[rng() % 4]);
            }
            if (!seqs.empty() && rng() % 4 == 0) seq += seqs[rng() % seqs.size()];   // repeats
            seqs.push_back(seq);
        }
        check(seqs, kmer_len, minimizer_len, rng, false);
    }
    std::printf("random battery done (%d failures)\n", g_fail);
    // Repeat families (human-like), with one unit copied hundreds to thousands of times with point mutations puts many super-k-mers behind a single minimizer, i.e. buckets beyond kSkewMin that the skew tables answer.
    // Exact duplicates of whole copies are included (repeated k-mers).
    uint64_t skew_seen = 0;
    for (int it = 0; it < 12; ++it) {
        const unsigned kmer_len = 31, minimizer_len = (it % 3 == 0) ? 11 : (it % 3 == 1) ? 15 : 21;
        std::string unit;
        for (int j = 0; j < 200 + int(rng() % 200); ++j) unit.push_back(bases[rng() % 4]);
        std::vector<std::string> seqs;
        const int copies = 300 + int(rng() % 2000);
        for (int c = 0; c < copies; ++c) {
            std::string seq;
            for (int j = 0; j < int(rng() % 40); ++j) seq.push_back(bases[rng() % 4]);   // flank
            const unsigned div = 1 + rng() % 5;                                      // 1-5% diverged
            for (char base_char : unit) seq.push_back(rng() % 100 < div ? bases[rng() % 4] : base_char);
            for (int j = 0; j < int(rng() % 40); ++j) seq.push_back(bases[rng() % 4]);
            seqs.push_back(seq);
            if (rng() % 10 == 0) seqs.push_back(seq);                                   // exact duplicate
        }
        SpssCore probe_idx;
        if (probe_idx.build(seqs, kmer_len, minimizer_len)) skew_seen += probe_idx.skew_entries();
        check(seqs, kmer_len, minimizer_len, rng, false);
    }
    CHECK(skew_seen > 0, "repeat battery never reached the skew tables");
    std::printf("repeat-family battery done (%d failures, %llu k-mers via skew tables)\n",
                g_fail, (unsigned long long)skew_seen);
    if (argc > 1) {
        auto seqs = read_fa(argv[1]);
        for (unsigned minimizer_len : {15u, 17u, 19u, 21u}) check(seqs, 31, minimizer_len, rng, true);
    }
    std::printf(g_fail ? "spss_stress: %d FAILURE(S)\n" : "spss_stress: ALL PASS\n", g_fail);
    return g_fail != 0;
}
