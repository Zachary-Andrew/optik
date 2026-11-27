// dbg_oracle.hpp - independent brute-force maximal-unitig oracle for the bidirected de Bruijn graph.
// It is test-only and shares no code with the extractor.
//
// The graph is a set of canonical k-mers plus an edge predicate over oriented k-mers, edge(E, X) with X = E[1..]+b.
// Node-centric graphs pass a predicate that only checks membership, while edge-centric graphs also require the (K+1)-mer window in either orientation.
//
// Maximal unitigs follow Conway & Bromage 2011 and Cuttlefish 2.
// Oriented E -> X is compactable iff E has exactly one successor X, X has exactly one predecessor, and canon(X) != canon(E).
// Each maximal compactable path is one unitig, and every k-mer lies in exactly one unitig.
// Palindromic k-mers, possible at even K, are singleton unitigs.
// Pure cycles are rotated to start at their smallest canonical k-mer in canonical orientation.
// Output strings are rc-normalized (min(s, rc(s))) for orientation-free comparison.
#pragma once
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

namespace dbg_oracle {
/// Largest k whose 2-bit k-mers fill a 64-bit word, and the number of bases.
constexpr unsigned kOracleMaxK = 32;
constexpr unsigned kOracleBases = 4;


inline uint64_t oracle_revcomp(uint64_t kmer, unsigned kmer_len) {
    uint64_t rev = 0;
    for (unsigned i = 0; i < kmer_len; ++i) { rev = (rev << 2) | ((kmer & 3u) ^ 3u); kmer >>= 2; }
    return rev;
}
inline uint64_t canon(uint64_t kmer, unsigned kmer_len) { uint64_t rev = oracle_revcomp(kmer, kmer_len); return kmer < rev ? kmer : rev; }
inline std::string rc_str(const std::string& seq) {
    std::string rev;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
        rev.push_back(*it == 'A' ? 'T' : *it == 'C' ? 'G' : *it == 'G' ? 'C' : 'A');
    return rev;
}
inline std::string norm(const std::string& seq) { std::string rev = rc_str(seq); return seq < rev ? seq : rev; }

using EdgeFn = std::function<bool(uint64_t from_kmer, uint64_t to_kmer)>;

inline std::multiset<std::string> unitigs(const std::set<uint64_t>& nodes,
                                          unsigned kmer_len, const EdgeFn& edge) {
    const uint64_t mask = kmer_len >= kOracleMaxK ? ~uint64_t(0) : ((uint64_t(1) << (2 * kmer_len)) - 1);
    auto succ = [&](uint64_t oriented) {
        std::vector<uint64_t> out;
        for (unsigned b = 0; b < kOracleBases; ++b) {
            uint64_t next_kmer = ((oriented << 2) | b) & mask;
            if (nodes.count(canon(next_kmer, kmer_len)) && edge(oriented, next_kmer)) out.push_back(next_kmer);
        }
        return out;
    };
    auto compact = [&](uint64_t oriented, uint64_t& next_kmer) {
        auto succs = succ(oriented);
        if (succs.size() != 1) return false;
        next_kmer = succs[0];
        if (canon(next_kmer, kmer_len) == canon(oriented, kmer_len)) return false;
        if (oriented == oracle_revcomp(oriented, kmer_len) || next_kmer == oracle_revcomp(next_kmer, kmer_len)) return false;   // palindromes stay single
        return succ(oracle_revcomp(next_kmer, kmer_len)).size() == 1;
    };
    auto dec = [&](uint64_t kmer) {
        std::string seq;
        for (int i = int(kmer_len) - 1; i >= 0; --i) seq.push_back("ACGT"[(kmer >> (2 * i)) & 3u]);
        return seq;
    };
    std::multiset<std::string> out;
    std::unordered_set<uint64_t> used;
    for (uint64_t canon_node : nodes) {
        if (used.count(canon_node)) continue;
        // Walk left to the start, detecting a pure cycle.
        uint64_t start = canon_node;
        bool cycle = false;
        for (;;) {
            uint64_t pred_rc;
            if (!compact(oracle_revcomp(start, kmer_len), pred_rc)) break;
            uint64_t pred = oracle_revcomp(pred_rc, kmer_len);
            if (canon(pred, kmer_len) == canon_node) { cycle = true; break; }
            start = pred;
        }
        if (cycle) {
            // rotate to the smallest canonical node, canonical orientation
            uint64_t smallest = canon_node, cur = canon_node;
            for (;;) {
                uint64_t next_kmer;
                if (!compact(cur, next_kmer)) break;          // defensive, since this is not a cycle
                if (canon(next_kmer, kmer_len) == canon_node) break;
                smallest = std::min(smallest, canon(next_kmer, kmer_len));
                cur = next_kmer;
            }
            start = smallest;
        }
        std::string seq = dec(start);
        used.insert(canon(start, kmer_len));
        uint64_t cur = start;
        for (;;) {
            uint64_t next_kmer;
            if (!compact(cur, next_kmer)) break;
            if (canon(next_kmer, kmer_len) == canon(start, kmer_len)) break;
            seq.push_back("ACGT"[next_kmer & 3u]);
            used.insert(canon(next_kmer, kmer_len));
            cur = next_kmer;
        }
        out.insert(norm(seq));
    }
    return out;
}

}
