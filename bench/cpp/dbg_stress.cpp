// dbg_stress - exactness battery for the de Bruijn graph build path that `optik dbruijn` runs, from FASTA through count_kmers_with_edges and DbgGraph::build_core to unitigs, out_mask, in_mask and save.
//
// Ground truth is brute force over the input sequences, sharing no code with the counter or the walker.
//   * solid nodes = canonical k-mers with >= min_count occurrences, broken at N
//   * edge E -> X = both endpoints solid AND the (K+1)-window E+x observed in
//                   the input in either orientation, an edge-centric graph
//   * unitigs     = maximal compaction of that graph by dbg_oracle::unitigs
// Checks per case
//   1. counter output, where the k-mer set, abundances and canonical masks == truth, identical at 1/4/8 threads
//   2. mask queries (out_mask/in_mask) == truth, 0 for absent k-mers
//   3. unitig multiset == oracle, every solid k-mer in exactly one unitig, unitig coverage == summed abundances
//   4. emission order byte-identical at 1/4/8 threads
//   5. save() writes the documented header
// The node-centric build (DynSpss::compact, a k-mer set without read windows in which every (K-1)-overlap is an edge) checks masks and unitigs against the same oracle with an always-true edge predicate.
// Regimes are dense small K (4..7, with palindromes at even K, self-loops, hairpins, N runs, sub-K records) at min_count 1..3, and contiguous genomes with low-complexity runs (K = 15/21/31) plus error-bearing reads at min_count 2.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "dbg_graph.hpp"
#include "dbg_oracle.hpp"
#include "kmer_count.hpp"

using namespace tpoptoa;

static int g_failures = 0;
#define CHECK(cond, ...)                                                    \
    do { if (!(cond)) { ++g_failures;                                       \
        std::fprintf(stderr, "FAIL %s:%d ", __FILE__, __LINE__);            \
        std::fprintf(stderr, __VA_ARGS__); std::fprintf(stderr, "\n"); } }  \
    while (0)

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static uint64_t xrng() {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static int code(char base_char) {
    switch (base_char) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}

struct Truth {
    std::map<uint64_t, uint32_t> counts;          // solid canonical k-mer -> abundance
    std::unordered_set<uint64_t> windows;         // observed oriented (K+1)-mers
    std::map<uint64_t, uint8_t> masks;            // solid canonical -> in<<4|out
};

static Truth brute(const std::vector<std::string>& seqs, unsigned kmer_len, uint32_t minc) {
    using dbg_oracle::canon;
    using dbg_oracle::oracle_revcomp;
    Truth truth;
    std::map<uint64_t, uint32_t> all;
    const uint64_t kmask = (uint64_t(1) << (2 * kmer_len)) - 1;
    const uint64_t win_mask = kmer_len + 1 >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * (kmer_len + 1))) - 1;
    for (const std::string& seq : seqs) {
        uint64_t fwd = 0;
        std::size_t valid = 0;
        for (char base_char : seq) {
            const int code_val = code(base_char);
            if (code_val < 0) { fwd = 0; valid = 0; continue; }
            fwd = (fwd << 2) | unsigned(code_val);
            ++valid;
            if (valid >= kmer_len) ++all[canon(fwd & kmask, kmer_len)];
            if (valid >= kmer_len + 1) truth.windows.insert(fwd & win_mask);
        }
    }
    for (auto& entry : all) if (entry.second >= minc) truth.counts.insert(entry);
    auto solid = [&](uint64_t kmer) { return truth.counts.count(canon(kmer, kmer_len)) != 0; };
    auto win = [&](uint64_t window) { return truth.windows.count(window) || truth.windows.count(oracle_revcomp(window, kmer_len + 1)); };
    for (auto& entry : truth.counts) {
        const uint64_t canon_kmer = entry.first;
        unsigned in_mask = 0, out = 0;
        for (unsigned b = 0; b < 4; ++b) {
            if (solid(((canon_kmer << 2) | b) & kmask) && win((canon_kmer << 2) | b)) out |= 1u << b;
            if (solid((uint64_t(b) << (2 * (kmer_len - 1))) | (canon_kmer >> 2)) &&
                win((uint64_t(b) << (2 * kmer_len)) | canon_kmer)) in_mask |= 1u << b;
        }
        truth.masks[canon_kmer] = uint8_t((in_mask << 4) | out);
    }
    return truth;
}

static std::string write_fasta(const std::vector<std::string>& seqs) {
    const char* dir = std::getenv("TMPDIR");
    std::string path = std::string(dir && *dir ? dir : "/tmp") + "/dbg_stress_XXXXXX";
    std::vector<char> buf(path.begin(), path.end());
    buf.push_back('\0');
    const int file_desc = mkstemp(buf.data());
    if (file_desc < 0) { std::perror("mkstemp"); std::exit(2); }
    std::string body;
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        body += ">s" + std::to_string(i) + "\n";
        for (std::size_t p = 0; p < seqs[i].size(); p += 70)        // multi-line records
            body += seqs[i].substr(p, 70) + "\n";
    }
    if (write(file_desc, body.data(), body.size()) != static_cast<ssize_t>(body.size())) {
        std::perror("write"); std::exit(2);
    }
    close(file_desc);
    return std::string(buf.data());
}

static void run_case(const std::string& name, const std::vector<std::string>& seqs,
                     unsigned kmer_len, uint32_t minc) {
    const Truth truth = brute(seqs, kmer_len, minc);
    const std::string fasta_path = write_fasta(seqs);

    // 1. counter output, thread-count invariant
    std::vector<std::pair<uint64_t, uint32_t>> pairs;
    std::vector<uint8_t> masks;
    for (unsigned num_threads : {1u, 4u, 8u}) {
        std::vector<std::pair<uint64_t, uint32_t>> pairs_got;
        std::vector<uint8_t> masks_got;
        detail::count_kmers_with_edges(fasta_path.c_str(), kmer_len, minc, num_threads, pairs_got, masks_got);
        if (num_threads == 1) { pairs = pairs_got; masks = masks_got; continue; }
        CHECK(pairs_got == pairs && masks_got == masks, "[%s] counter output differs at %u threads",
              name.c_str(), num_threads);
    }
    std::remove(fasta_path.c_str());
    CHECK(pairs.size() == truth.counts.size(), "[%s] K=%u solid k-mers got=%zu want=%zu",
          name.c_str(), kmer_len, pairs.size(), truth.counts.size());
    CHECK(masks.size() == pairs.size(), "[%s] masks not index-aligned", name.c_str());
    {
        std::size_t rec_idx = 0, bad = 0;
        for (auto& entry : truth.counts) {
            if (rec_idx >= pairs.size()) break;
            if (pairs[rec_idx].first != entry.first || pairs[rec_idx].second != entry.second ||
                masks[rec_idx] != truth.masks.at(entry.first)) {
                if (++bad <= 3)
                    std::fprintf(stderr, "    [%s] i=%zu got (%llu,%u,%02x) want (%llu,%u,%02x)\n",
                                 name.c_str(), rec_idx, (unsigned long long)pairs[rec_idx].first,
                                 pairs[rec_idx].second, masks[rec_idx], (unsigned long long)entry.first,
                                 entry.second, truth.masks.at(entry.first));
            }
            ++rec_idx;
        }
        CHECK(bad == 0, "[%s] K=%u %zu counter records differ from truth", name.c_str(), kmer_len, bad);
    }
    if (truth.counts.empty()) {
        std::fprintf(stderr, "[%s] K=%u min_count=%u empty graph OK\n", name.c_str(), kmer_len, minc);
        return;
    }

    DbgGraph graph(kmer_len);
    graph.build_core(std::move(pairs), masks);
    CHECK(graph.core_size() == truth.counts.size(), "[%s] core_size", name.c_str());

    // 2. mask queries, including absent k-mers
    for (auto& entry : truth.masks) {
        const unsigned want = entry.second;
        CHECK(graph.out_mask(entry.first) == (want & 0xFu) && graph.in_mask(entry.first) == (want >> 4),
              "[%s] mask query mismatch key=%llu", name.c_str(), (unsigned long long)entry.first);
    }
    const uint64_t kmask = (uint64_t(1) << (2 * kmer_len)) - 1;
    for (int i = 0; i < 200; ++i) {
        const uint64_t canon_kmer = dbg_oracle::canon(xrng() & kmask, kmer_len);
        if (truth.counts.count(canon_kmer)) continue;
        CHECK(graph.out_mask(canon_kmer) == 0 && graph.in_mask(canon_kmer) == 0, "[%s] absent key has a mask",
              name.c_str());
    }

    // 3. unitigs vs the oracle, then 4. emission order across thread counts
    std::vector<std::pair<std::string, uint64_t>> emitted[3];
    const unsigned nts[3] = {1u, 4u, 8u};
    for (int i = 0; i < 3; ++i)
        graph.unitigs([&](const uint8_t* seq, uint32_t len, uint64_t cov) {
            emitted[i].push_back({std::string(reinterpret_cast<const char*>(seq), len), cov});
        }, nts[i]);
    CHECK(emitted[0] == emitted[1] && emitted[0] == emitted[2],
          "[%s] unitig emission differs across thread counts", name.c_str());

    std::set<uint64_t> nodes;
    for (auto& entry : truth.counts) nodes.insert(entry.first);
    auto win = [&](uint64_t from_kmer, uint64_t to_kmer) {
        const uint64_t window = (from_kmer << 2) | (to_kmer & 3u);
        return truth.windows.count(window) || truth.windows.count(dbg_oracle::oracle_revcomp(window, kmer_len + 1));
    };
    const std::multiset<std::string> want = dbg_oracle::unitigs(nodes, kmer_len, win);
    std::multiset<std::string> got;
    std::unordered_map<uint64_t, unsigned> seen;
    std::size_t cov_bad = 0;
    for (auto& [unitig, cov] : emitted[0]) {
        got.insert(dbg_oracle::norm(unitig));
        uint64_t fwd = 0, sum = 0;
        for (std::size_t i = 0; i < unitig.size(); ++i) {
            fwd = ((fwd << 2) | unsigned(code(unitig[i]))) & kmask;
            if (i + 1 < kmer_len) continue;
            const uint64_t canon_kmer = dbg_oracle::canon(fwd, kmer_len);
            ++seen[canon_kmer];
            auto found = truth.counts.find(canon_kmer);
            sum += found == truth.counts.end() ? 0 : found->second;
        }
        if (sum != cov) ++cov_bad;
    }
    CHECK(got == want, "[%s] K=%u unitig multiset mismatch got=%zu want=%zu",
          name.c_str(), kmer_len, got.size(), want.size());
    if (got != want) {
        std::multiset<std::string> miss = want, extra = got;
        for (auto& unitig : got)  { auto found = miss.find(unitig);  if (found != miss.end())  miss.erase(found); }
        for (auto& unitig : want) { auto found = extra.find(unitig); if (found != extra.end()) extra.erase(found); }
        int shown = 0;
        for (auto& unitig : miss)  { std::fprintf(stderr, "    MISSING %s\n", unitig.c_str()); if (++shown == 4) break; }
        shown = 0;
        for (auto& unitig : extra) { std::fprintf(stderr, "    EXTRA   %s\n", unitig.c_str()); if (++shown == 4) break; }
    }
    std::size_t dup = 0;
    for (auto& entry : seen) if (entry.second != 1 || !truth.counts.count(entry.first)) ++dup;
    CHECK(seen.size() == truth.counts.size() && dup == 0,
          "[%s] unitigs are not a partition of the solid k-mers (%zu covered, %zu bad)",
          name.c_str(), seen.size(), dup);
    CHECK(cov_bad == 0, "[%s] %zu unitigs with coverage != summed abundance",
          name.c_str(), cov_bad);

    // 5. save
    {
        const std::string out = write_fasta({}) + ".kdbg";
        CHECK(graph.save(out.c_str()), "[%s] save failed", name.c_str());
        std::FILE* file = std::fopen(out.c_str(), "rb");
        char hdr[7] = {0};
        const bool valid_header = file && std::fread(hdr, 1, 7, file) == 7 &&
                        std::memcmp(hdr, "OKKDBG", 6) == 0 && unsigned(uint8_t(hdr[6])) == kmer_len;
        if (file) std::fclose(file);
        CHECK(valid_header, "[%s] save header", name.c_str());
        std::remove(out.c_str());
        std::remove(out.substr(0, out.size() - 5).c_str());
    }
    std::fprintf(stderr, "[%s] K=%u min_count=%u nodes=%zu unitigs=%zu %s\n",
                 name.c_str(), kmer_len, minc, truth.counts.size(), got.size(),
                 got == want ? "OK" : "MISMATCH");
}

static void run_node_centric(const std::string& name, const std::set<uint64_t>& nodes,
                             unsigned kmer_len) {
    std::vector<std::pair<uint64_t, uint32_t>> pairs;
    for (uint64_t canon_kmer : nodes) pairs.push_back({canon_kmer, 1});
    DbgGraph graph(kmer_len);
    graph.build_core(std::move(pairs));
    const uint64_t kmask = (uint64_t(1) << (2 * kmer_len)) - 1;
    auto present = [&](uint64_t kmer) { return nodes.count(dbg_oracle::canon(kmer, kmer_len)) != 0; };
    std::size_t bad = 0;
    for (uint64_t canon_kmer : nodes) {
        unsigned in_mask = 0, out = 0;
        for (unsigned b = 0; b < 4; ++b) {
            if (present(((canon_kmer << 2) | b) & kmask)) out |= 1u << b;
            if (present((uint64_t(b) << (2 * (kmer_len - 1))) | (canon_kmer >> 2))) in_mask |= 1u << b;
        }
        if (graph.out_mask(canon_kmer) != out || graph.in_mask(canon_kmer) != in_mask) ++bad;
    }
    CHECK(bad == 0, "[%s] %zu node-centric masks differ", name.c_str(), bad);
    std::vector<std::string> emitted1, emitted8;
    graph.unitigs([&](const uint8_t* seq, uint32_t len, uint64_t) { emitted1.emplace_back((const char*)seq, len); }, 1);
    graph.unitigs([&](const uint8_t* seq, uint32_t len, uint64_t) { emitted8.emplace_back((const char*)seq, len); }, 8);
    CHECK(emitted1 == emitted8, "[%s] node-centric emission differs across thread counts", name.c_str());
    std::multiset<std::string> got;
    for (auto& unitig : emitted1) got.insert(dbg_oracle::norm(unitig));
    const auto want = dbg_oracle::unitigs(nodes, kmer_len, [](uint64_t, uint64_t) { return true; });
    CHECK(got == want, "[%s] K=%u node-centric unitigs got=%zu want=%zu",
          name.c_str(), kmer_len, got.size(), want.size());
    std::fprintf(stderr, "[%s] K=%u node-centric nodes=%zu unitigs=%zu %s\n", name.c_str(),
                 kmer_len, nodes.size(), got.size(), got == want ? "OK" : "MISMATCH");
}

static std::string rand_seq(std::size_t len) {
    std::string seq(len, 'A');
    for (auto& base_char : seq) base_char = "ACGT"[xrng() & 3];
    return seq;
}

static std::string genome(std::size_t len) {
    std::string genome_seq;
    genome_seq.reserve(len);
    while (genome_seq.size() < len) {
        if (xrng() % 23 == 0) genome_seq.append(3 + xrng() % 20, "ACGT"[xrng() & 3]);   // homopolymer
        else if (xrng() % 41 == 0) {                                          // tandem repeat
            const std::string unit = rand_seq(2 + xrng() % 5);
            for (int r = 0, repeats = 2 + int(xrng() % 6); r < repeats; ++r) genome_seq += unit;
        } else genome_seq.push_back("ACGT"[xrng() & 3]);
    }
    genome_seq.resize(len);
    return genome_seq;
}

int main() {
    std::fprintf(stderr, "== dbg_stress ==\n");
    int case_id = 0;

    // Dense small-K graphs hold most of the k-mer space, so palindromes at even K, self-loops, hairpins and pure cycles all occur.
    for (unsigned kmer_len = 4; kmer_len <= 7; ++kmer_len)
        for (uint32_t minc = 1; minc <= 3; ++minc)
            for (int rep = 0; rep < 3; ++rep, ++case_id) {
                std::vector<std::string> seqs;
                const int num_seqs = 5 + int(xrng() % 40);
                for (int i = 0; i < num_seqs; ++i) {
                    std::string seq = rand_seq(xrng() % (6 * kmer_len));          // includes sub-K records
                    if (seq.size() > kmer_len && xrng() % 4 == 0)
                        seq[xrng() % seq.size()] = 'N';                      // break the k-mer chain
                    if (xrng() % 7 == 0) seq = dbg_oracle::rc_str(seq);
                    seqs.push_back(seq);
                }
                run_case("dense" + std::to_string(case_id), seqs, kmer_len, minc);
            }

    // Degenerate inputs.
    run_case("empty", {}, 5, 1);
    run_case("subk", {"ACG", "TTA", ""}, 5, 1);
    run_case("allN", {"NNNNNNNNNNNN"}, 5, 1);
    run_case("homopolymer", {std::string(50, 'A')}, 5, 1);
    run_case("palindrome", {"ACGT", "ACGTACGT"}, 4, 1);
    run_case("nothing-solid", {"ACGTACGTTGCA"}, 5, 9);   // nothing solid

    // Abundances beyond 24 bits (formerly packed with the masks and saturated at 2^24 - 1), where node abundance and unitig coverage are exact up to the full 32 bits.
    {
        constexpr unsigned kmer_len = 5;
        const uint64_t kmer_a = 0x1B, kmer_b = 0x2C;                 // two canonical, unlinked 5-mers
        const uint32_t big[2] = {(1u << 24) + 5, 0xFFFFFFFFu};
        std::vector<std::pair<uint64_t, uint32_t>> pairs = {{kmer_a, big[0]}, {kmer_b, big[1]}};
        std::sort(pairs.begin(), pairs.end());
        std::vector<uint8_t> masks(2, 0);
        DbgGraph graph(kmer_len);
        graph.build_core(std::move(pairs), masks);
        std::vector<uint64_t> cov;
        graph.unitigs([&](const uint8_t*, uint32_t, uint64_t cov_val) { cov.push_back(cov_val); }, 1);
        std::sort(cov.begin(), cov.end());
        CHECK(cov.size() == 2 && cov[0] == big[0] && cov[1] == big[1],
              "[abundance-32bit] unitig coverage %llu/%llu, want %u/%u",
              cov.size() > 0 ? (unsigned long long)cov[0] : 0ull,
              cov.size() > 1 ? (unsigned long long)cov[1] : 0ull, big[0], big[1]);
        std::fprintf(stderr, "[abundance-32bit] coverage exact above 2^24\n");
    }

    // Contiguous genomes at realistic K, and error-bearing reads at min_count 2.
    for (unsigned kmer_len : {15u, 21u, 31u})
        for (int trial = 0; trial < 3; ++trial) {
            const std::string gnm = genome(20000 + xrng() % 20000);
            run_case("genome-k" + std::to_string(kmer_len) + "-" + std::to_string(trial), {gnm}, kmer_len, 1);
            std::vector<std::string> reads;
            for (int r = 0; r < 600; ++r) {
                const std::size_t len = 150, pos = xrng() % (gnm.size() - len);
                std::string read = gnm.substr(pos, len);
                for (auto& base_char : read) if (xrng() % 100 == 0) base_char = "ACGT"[xrng() & 3];
                if (xrng() & 1) read = dbg_oracle::rc_str(read);
                reads.push_back(read);
            }
            run_case("reads-k" + std::to_string(kmer_len) + "-" + std::to_string(trial), reads, kmer_len, 2);
        }

    // Node-centric builds on random dense k-mer sets and on a genome's k-mer set
    for (unsigned kmer_len = 4; kmer_len <= 9; ++kmer_len)
        for (int rep = 0; rep < 3; ++rep) {
            std::set<uint64_t> nodes;
            const uint64_t kmask = (uint64_t(1) << (2 * kmer_len)) - 1;
            const std::size_t num_nodes = 1 + xrng() % std::min<uint64_t>(kmask / 2, 3000);
            while (nodes.size() < num_nodes) nodes.insert(dbg_oracle::canon(xrng() & kmask, kmer_len));
            run_node_centric("nc-dense" + std::to_string(kmer_len) + "-" + std::to_string(rep), nodes, kmer_len);
        }
    for (unsigned kmer_len : {15u, 31u}) {
        const std::string gnm = genome(30000);
        std::set<uint64_t> nodes;
        const uint64_t kmask = kmer_len >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * kmer_len)) - 1;
        uint64_t fwd = 0;
        for (std::size_t i = 0; i < gnm.size(); ++i) {
            fwd = ((fwd << 2) | unsigned(code(gnm[i]))) & kmask;
            if (i + 1 >= kmer_len) nodes.insert(dbg_oracle::canon(fwd, kmer_len));
        }
        run_node_centric("nc-genome-k" + std::to_string(kmer_len), nodes, kmer_len);
    }

    std::fprintf(stderr, "%s (%d failures)\n", g_failures ? "FAILED" : "ALL OK", g_failures);
    return g_failures ? 1 : 0;
}
