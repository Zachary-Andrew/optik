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

static int code(char ch) {
    switch (ch) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}

struct Truth {
    std::map<uint64_t, uint32_t> counts;          // solid canonical k-mer -> abundance
    std::unordered_set<uint64_t> windows;         // observed oriented (K+1)-mers
    std::map<uint64_t, uint8_t> masks;            // solid canonical -> in<<4|out
};

static Truth brute(const std::vector<std::string>& seqs, unsigned K, uint32_t minc) {
    using dbg_oracle::canon;
    using dbg_oracle::rc;
    Truth t;
    std::map<uint64_t, uint32_t> all;
    const uint64_t km = (uint64_t(1) << (2 * K)) - 1;
    const uint64_t wm = K + 1 >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * (K + 1))) - 1;
    for (const std::string& s : seqs) {
        uint64_t f = 0;
        std::size_t valid = 0;
        for (char ch : s) {
            const int b = code(ch);
            if (b < 0) { f = 0; valid = 0; continue; }
            f = (f << 2) | unsigned(b);
            ++valid;
            if (valid >= K) ++all[canon(f & km, K)];
            if (valid >= K + 1) t.windows.insert(f & wm);
        }
    }
    for (auto& kv : all) if (kv.second >= minc) t.counts.insert(kv);
    auto solid = [&](uint64_t x) { return t.counts.count(canon(x, K)) != 0; };
    auto win = [&](uint64_t m) { return t.windows.count(m) || t.windows.count(rc(m, K + 1)); };
    for (auto& kv : t.counts) {
        const uint64_t c = kv.first;
        unsigned in = 0, out = 0;
        for (unsigned b = 0; b < 4; ++b) {
            if (solid(((c << 2) | b) & km) && win((c << 2) | b)) out |= 1u << b;
            if (solid((uint64_t(b) << (2 * (K - 1))) | (c >> 2)) &&
                win((uint64_t(b) << (2 * K)) | c)) in |= 1u << b;
        }
        t.masks[c] = uint8_t((in << 4) | out);
    }
    return t;
}

static std::string write_fasta(const std::vector<std::string>& seqs) {
    const char* dir = std::getenv("TMPDIR");
    std::string path = std::string(dir && *dir ? dir : "/tmp") + "/dbg_stress_XXXXXX";
    std::vector<char> buf(path.begin(), path.end());
    buf.push_back('\0');
    const int fd = mkstemp(buf.data());
    if (fd < 0) { std::perror("mkstemp"); std::exit(2); }
    std::string body;
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        body += ">s" + std::to_string(i) + "\n";
        for (std::size_t p = 0; p < seqs[i].size(); p += 70)        // multi-line records
            body += seqs[i].substr(p, 70) + "\n";
    }
    if (write(fd, body.data(), body.size()) != static_cast<ssize_t>(body.size())) {
        std::perror("write"); std::exit(2);
    }
    close(fd);
    return std::string(buf.data());
}

static void run_case(const std::string& name, const std::vector<std::string>& seqs,
                     unsigned K, uint32_t minc) {
    const Truth t = brute(seqs, K, minc);
    const std::string fa = write_fasta(seqs);

    // 1. counter output, thread-count invariant
    std::vector<std::pair<uint64_t, uint32_t>> pairs;
    std::vector<uint8_t> masks;
    for (unsigned nt : {1u, 4u, 8u}) {
        std::vector<std::pair<uint64_t, uint32_t>> p;
        std::vector<uint8_t> m;
        detail::count_kmers_with_edges(fa.c_str(), K, minc, nt, p, m);
        if (nt == 1) { pairs = p; masks = m; continue; }
        CHECK(p == pairs && m == masks, "[%s] counter output differs at %u threads",
              name.c_str(), nt);
    }
    std::remove(fa.c_str());
    CHECK(pairs.size() == t.counts.size(), "[%s] K=%u solid k-mers got=%zu want=%zu",
          name.c_str(), K, pairs.size(), t.counts.size());
    CHECK(masks.size() == pairs.size(), "[%s] masks not index-aligned", name.c_str());
    {
        std::size_t i = 0, bad = 0;
        for (auto& kv : t.counts) {
            if (i >= pairs.size()) break;
            if (pairs[i].first != kv.first || pairs[i].second != kv.second ||
                masks[i] != t.masks.at(kv.first)) {
                if (++bad <= 3)
                    std::fprintf(stderr, "    [%s] i=%zu got (%llu,%u,%02x) want (%llu,%u,%02x)\n",
                                 name.c_str(), i, (unsigned long long)pairs[i].first,
                                 pairs[i].second, masks[i], (unsigned long long)kv.first,
                                 kv.second, t.masks.at(kv.first));
            }
            ++i;
        }
        CHECK(bad == 0, "[%s] K=%u %zu counter records differ from truth", name.c_str(), K, bad);
    }
    if (t.counts.empty()) {
        std::fprintf(stderr, "[%s] K=%u min_count=%u empty graph OK\n", name.c_str(), K, minc);
        return;
    }

    DbgGraph g(K);
    g.build_core(std::move(pairs), masks);
    CHECK(g.core_size() == t.counts.size(), "[%s] core_size", name.c_str());

    // 2. mask queries, including absent k-mers
    for (auto& kv : t.masks) {
        const unsigned want = kv.second;
        CHECK(g.out_mask(kv.first) == (want & 0xFu) && g.in_mask(kv.first) == (want >> 4),
              "[%s] mask query mismatch key=%llu", name.c_str(), (unsigned long long)kv.first);
    }
    const uint64_t km = (uint64_t(1) << (2 * K)) - 1;
    for (int i = 0; i < 200; ++i) {
        const uint64_t c = dbg_oracle::canon(xrng() & km, K);
        if (t.counts.count(c)) continue;
        CHECK(g.out_mask(c) == 0 && g.in_mask(c) == 0, "[%s] absent key has a mask",
              name.c_str());
    }

    // 3. unitigs vs the oracle, then 4. emission order across thread counts
    std::vector<std::pair<std::string, uint64_t>> emitted[3];
    const unsigned nts[3] = {1u, 4u, 8u};
    for (int i = 0; i < 3; ++i)
        g.unitigs([&](const uint8_t* s, uint32_t len, uint64_t cov) {
            emitted[i].push_back({std::string(reinterpret_cast<const char*>(s), len), cov});
        }, nts[i]);
    CHECK(emitted[0] == emitted[1] && emitted[0] == emitted[2],
          "[%s] unitig emission differs across thread counts", name.c_str());

    std::set<uint64_t> nodes;
    for (auto& kv : t.counts) nodes.insert(kv.first);
    auto win = [&](uint64_t E, uint64_t X) {
        const uint64_t m = (E << 2) | (X & 3u);
        return t.windows.count(m) || t.windows.count(dbg_oracle::rc(m, K + 1));
    };
    const std::multiset<std::string> want = dbg_oracle::unitigs(nodes, K, win);
    std::multiset<std::string> got;
    std::unordered_map<uint64_t, unsigned> seen;
    std::size_t cov_bad = 0;
    for (auto& [u, cov] : emitted[0]) {
        got.insert(dbg_oracle::norm(u));
        uint64_t f = 0, sum = 0;
        for (std::size_t i = 0; i < u.size(); ++i) {
            f = ((f << 2) | unsigned(code(u[i]))) & km;
            if (i + 1 < K) continue;
            const uint64_t c = dbg_oracle::canon(f, K);
            ++seen[c];
            auto it = t.counts.find(c);
            sum += it == t.counts.end() ? 0 : it->second;
        }
        if (sum != cov) ++cov_bad;
    }
    CHECK(got == want, "[%s] K=%u unitig multiset mismatch got=%zu want=%zu",
          name.c_str(), K, got.size(), want.size());
    if (got != want) {
        std::multiset<std::string> miss = want, extra = got;
        for (auto& u : got)  { auto it = miss.find(u);  if (it != miss.end())  miss.erase(it); }
        for (auto& u : want) { auto it = extra.find(u); if (it != extra.end()) extra.erase(it); }
        int n = 0;
        for (auto& u : miss)  { std::fprintf(stderr, "    MISSING %s\n", u.c_str()); if (++n == 4) break; }
        n = 0;
        for (auto& u : extra) { std::fprintf(stderr, "    EXTRA   %s\n", u.c_str()); if (++n == 4) break; }
    }
    std::size_t dup = 0;
    for (auto& kv : seen) if (kv.second != 1 || !t.counts.count(kv.first)) ++dup;
    CHECK(seen.size() == t.counts.size() && dup == 0,
          "[%s] unitigs are not a partition of the solid k-mers (%zu covered, %zu bad)",
          name.c_str(), seen.size(), dup);
    CHECK(cov_bad == 0, "[%s] %zu unitigs with coverage != summed abundance",
          name.c_str(), cov_bad);

    // 5. save
    {
        const std::string out = write_fasta({}) + ".kdbg";
        CHECK(g.save(out.c_str()), "[%s] save failed", name.c_str());
        std::FILE* f = std::fopen(out.c_str(), "rb");
        char hdr[7] = {0};
        const bool ok = f && std::fread(hdr, 1, 7, f) == 7 &&
                        std::memcmp(hdr, "OKKDBG", 6) == 0 && unsigned(uint8_t(hdr[6])) == K;
        if (f) std::fclose(f);
        CHECK(ok, "[%s] save header", name.c_str());
        std::remove(out.c_str());
        std::remove(out.substr(0, out.size() - 5).c_str());
    }
    std::fprintf(stderr, "[%s] K=%u min_count=%u nodes=%zu unitigs=%zu %s\n",
                 name.c_str(), K, minc, t.counts.size(), got.size(),
                 got == want ? "OK" : "MISMATCH");
}

static void run_node_centric(const std::string& name, const std::set<uint64_t>& nodes,
                             unsigned K) {
    std::vector<std::pair<uint64_t, uint32_t>> pairs;
    for (uint64_t c : nodes) pairs.push_back({c, 1});
    DbgGraph g(K);
    g.build_core(std::move(pairs));
    const uint64_t km = (uint64_t(1) << (2 * K)) - 1;
    auto present = [&](uint64_t x) { return nodes.count(dbg_oracle::canon(x, K)) != 0; };
    std::size_t bad = 0;
    for (uint64_t c : nodes) {
        unsigned in = 0, out = 0;
        for (unsigned b = 0; b < 4; ++b) {
            if (present(((c << 2) | b) & km)) out |= 1u << b;
            if (present((uint64_t(b) << (2 * (K - 1))) | (c >> 2))) in |= 1u << b;
        }
        if (g.out_mask(c) != out || g.in_mask(c) != in) ++bad;
    }
    CHECK(bad == 0, "[%s] %zu node-centric masks differ", name.c_str(), bad);
    std::vector<std::string> e1, e8;
    g.unitigs([&](const uint8_t* s, uint32_t len, uint64_t) { e1.emplace_back((const char*)s, len); }, 1);
    g.unitigs([&](const uint8_t* s, uint32_t len, uint64_t) { e8.emplace_back((const char*)s, len); }, 8);
    CHECK(e1 == e8, "[%s] node-centric emission differs across thread counts", name.c_str());
    std::multiset<std::string> got;
    for (auto& u : e1) got.insert(dbg_oracle::norm(u));
    const auto want = dbg_oracle::unitigs(nodes, K, [](uint64_t, uint64_t) { return true; });
    CHECK(got == want, "[%s] K=%u node-centric unitigs got=%zu want=%zu",
          name.c_str(), K, got.size(), want.size());
    std::fprintf(stderr, "[%s] K=%u node-centric nodes=%zu unitigs=%zu %s\n", name.c_str(),
                 K, nodes.size(), got.size(), got == want ? "OK" : "MISMATCH");
}

static std::string rand_seq(std::size_t len) {
    std::string s(len, 'A');
    for (auto& ch : s) ch = "ACGT"[xrng() & 3];
    return s;
}

static std::string genome(std::size_t L) {
    std::string g;
    g.reserve(L);
    while (g.size() < L) {
        if (xrng() % 23 == 0) g.append(3 + xrng() % 20, "ACGT"[xrng() & 3]);   // homopolymer
        else if (xrng() % 41 == 0) {                                          // tandem repeat
            const std::string unit = rand_seq(2 + xrng() % 5);
            for (int r = 0, n = 2 + int(xrng() % 6); r < n; ++r) g += unit;
        } else g.push_back("ACGT"[xrng() & 3]);
    }
    g.resize(L);
    return g;
}

int main() {
    std::fprintf(stderr, "== dbg_stress ==\n");
    int id = 0;

    // Dense small-K graphs hold most of the k-mer space, so palindromes at even K, self-loops, hairpins and pure cycles all occur.
    for (unsigned K = 4; K <= 7; ++K)
        for (uint32_t minc = 1; minc <= 3; ++minc)
            for (int rep = 0; rep < 3; ++rep, ++id) {
                std::vector<std::string> seqs;
                const int n = 5 + int(xrng() % 40);
                for (int i = 0; i < n; ++i) {
                    std::string s = rand_seq(xrng() % (6 * K));          // includes sub-K records
                    if (s.size() > K && xrng() % 4 == 0)
                        s[xrng() % s.size()] = 'N';                      // break the k-mer chain
                    if (xrng() % 7 == 0) s = dbg_oracle::rc_str(s);
                    seqs.push_back(s);
                }
                run_case("dense" + std::to_string(id), seqs, K, minc);
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
        constexpr unsigned K = 5;
        const uint64_t a = 0x1B, b = 0x2C;                 // two canonical, unlinked 5-mers
        const uint32_t big[2] = {(1u << 24) + 5, 0xFFFFFFFFu};
        std::vector<std::pair<uint64_t, uint32_t>> pairs = {{a, big[0]}, {b, big[1]}};
        std::sort(pairs.begin(), pairs.end());
        std::vector<uint8_t> masks(2, 0);
        DbgGraph g(K);
        g.build_core(std::move(pairs), masks);
        std::vector<uint64_t> cov;
        g.unitigs([&](const uint8_t*, uint32_t, uint64_t c) { cov.push_back(c); }, 1);
        std::sort(cov.begin(), cov.end());
        CHECK(cov.size() == 2 && cov[0] == big[0] && cov[1] == big[1],
              "[abundance-32bit] unitig coverage %llu/%llu, want %u/%u",
              cov.size() > 0 ? (unsigned long long)cov[0] : 0ull,
              cov.size() > 1 ? (unsigned long long)cov[1] : 0ull, big[0], big[1]);
        std::fprintf(stderr, "[abundance-32bit] coverage exact above 2^24\n");
    }

    // Contiguous genomes at realistic K, and error-bearing reads at min_count 2.
    for (unsigned K : {15u, 21u, 31u})
        for (int trial = 0; trial < 3; ++trial) {
            const std::string gnm = genome(20000 + xrng() % 20000);
            run_case("genome-k" + std::to_string(K) + "-" + std::to_string(trial), {gnm}, K, 1);
            std::vector<std::string> reads;
            for (int r = 0; r < 600; ++r) {
                const std::size_t len = 150, pos = xrng() % (gnm.size() - len);
                std::string rd = gnm.substr(pos, len);
                for (auto& ch : rd) if (xrng() % 100 == 0) ch = "ACGT"[xrng() & 3];
                if (xrng() & 1) rd = dbg_oracle::rc_str(rd);
                reads.push_back(rd);
            }
            run_case("reads-k" + std::to_string(K) + "-" + std::to_string(trial), reads, K, 2);
        }

    // Node-centric builds on random dense k-mer sets and on a genome's k-mer set
    for (unsigned K = 4; K <= 9; ++K)
        for (int rep = 0; rep < 3; ++rep) {
            std::set<uint64_t> nodes;
            const uint64_t km = (uint64_t(1) << (2 * K)) - 1;
            const std::size_t n = 1 + xrng() % std::min<uint64_t>(km / 2, 3000);
            while (nodes.size() < n) nodes.insert(dbg_oracle::canon(xrng() & km, K));
            run_node_centric("nc-dense" + std::to_string(K) + "-" + std::to_string(rep), nodes, K);
        }
    for (unsigned K : {15u, 31u}) {
        const std::string gnm = genome(30000);
        std::set<uint64_t> nodes;
        const uint64_t km = K >= 32 ? ~uint64_t(0) : (uint64_t(1) << (2 * K)) - 1;
        uint64_t f = 0;
        for (std::size_t i = 0; i < gnm.size(); ++i) {
            f = ((f << 2) | unsigned(code(gnm[i]))) & km;
            if (i + 1 >= K) nodes.insert(dbg_oracle::canon(f, K));
        }
        run_node_centric("nc-genome-k" + std::to_string(K), nodes, K);
    }

    std::fprintf(stderr, "%s (%d failures)\n", g_failures ? "FAILED" : "ALL OK", g_failures);
    return g_failures ? 1 : 0;
}
