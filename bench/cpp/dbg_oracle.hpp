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

inline uint64_t rc(uint64_t x, unsigned K) {
    uint64_t r = 0;
    for (unsigned i = 0; i < K; ++i) { r = (r << 2) | ((x & 3u) ^ 3u); x >>= 2; }
    return r;
}
inline uint64_t canon(uint64_t x, unsigned K) { uint64_t r = rc(x, K); return x < r ? x : r; }
inline std::string rc_str(const std::string& s) {
    std::string r;
    for (auto it = s.rbegin(); it != s.rend(); ++it)
        r.push_back(*it == 'A' ? 'T' : *it == 'C' ? 'G' : *it == 'G' ? 'C' : 'A');
    return r;
}
inline std::string norm(const std::string& s) { std::string r = rc_str(s); return s < r ? s : r; }

using EdgeFn = std::function<bool(uint64_t E, uint64_t X)>;

inline std::multiset<std::string> unitigs(const std::set<uint64_t>& nodes,
                                          unsigned K, const EdgeFn& edge) {
    const uint64_t mask = K >= 32 ? ~uint64_t(0) : ((uint64_t(1) << (2 * K)) - 1);
    auto succ = [&](uint64_t E) {
        std::vector<uint64_t> out;
        for (unsigned b = 0; b < 4; ++b) {
            uint64_t X = ((E << 2) | b) & mask;
            if (nodes.count(canon(X, K)) && edge(E, X)) out.push_back(X);
        }
        return out;
    };
    auto compact = [&](uint64_t E, uint64_t& X) {
        auto s = succ(E);
        if (s.size() != 1) return false;
        X = s[0];
        if (canon(X, K) == canon(E, K)) return false;
        if (E == rc(E, K) || X == rc(X, K)) return false;   // palindromes stay single
        return succ(rc(X, K)).size() == 1;
    };
    auto dec = [&](uint64_t x) {
        std::string s;
        for (int i = int(K) - 1; i >= 0; --i) s.push_back("ACGT"[(x >> (2 * i)) & 3u]);
        return s;
    };
    std::multiset<std::string> out;
    std::unordered_set<uint64_t> used;
    for (uint64_t c : nodes) {
        if (used.count(c)) continue;
        // Walk left to the start, detecting a pure cycle.
        uint64_t E = c;
        bool cycle = false;
        for (;;) {
            uint64_t Y;
            if (!compact(rc(E, K), Y)) break;
            uint64_t P = rc(Y, K);
            if (canon(P, K) == c) { cycle = true; break; }
            E = P;
        }
        if (cycle) {
            // rotate to the smallest canonical node, canonical orientation
            uint64_t mn = c, cur = c;
            for (;;) {
                uint64_t X;
                if (!compact(cur, X)) break;          // defensive, since this is not a cycle
                if (canon(X, K) == c) break;
                mn = std::min(mn, canon(X, K));
                cur = X;
            }
            E = mn;
        }
        std::string s = dec(E);
        used.insert(canon(E, K));
        uint64_t cur = E;
        for (;;) {
            uint64_t X;
            if (!compact(cur, X)) break;
            if (canon(X, K) == canon(E, K)) break;
            s.push_back("ACGT"[X & 3u]);
            used.insert(canon(X, K));
            cur = X;
        }
        out.insert(norm(s));
    }
    return out;
}

}
