#pragma once

#include "thread_budget.hpp"
#include "phase_timer.hpp"

/*
 * dbg_unitigs.hpp - maximal unitig extraction for DbgGraph and DbgGraphEF, included at the end of dbg_graph.hpp.
 *
 * Bidirected compaction after Conway and Bromage 2011 and Cuttlefish 2 (Khan et al. 2022).
 * Nodes are canonical k-mers, each with canonical in/out edge masks, which are symmetric, so an edge is visible from both of its endpoints.
 * An oriented k-mer E extends to X iff E has exactly one successor X, X has exactly one predecessor and canon(X) != canon(E).
 * Maximal unitigs are the maximal paths of such steps, and run_sorted() states the exact rules.
 *
 * Parallelism, lock-free.
 * Every oriented start is walked independently in O(N) total, and each unitig is emitted once, from its smaller oriented start, into a per-thread buffer.
 * Nodes are marked with relaxed atomic ORs, and the unmarked remainder, the pure cycles, is handled serially.
 * The output is sorted by head k-mer, so it is identical at every thread count.
 *
 * bench/validate_unitigs.py checks that the partition has no duplicate k-mers, the path, interior and maximality rules, exact L links, and that the unitig partition equals Cuttlefish 2's --read output at min_count 1.
 */

#include <algorithm>
#include <atomic>
#include <memory>
#include <string>

namespace tpoptoa {

struct DbgUnitigExtractor {
    static constexpr std::size_t kNpos = DbgGraph::kNpos;

    struct Result {
        uint64_t             head = 0;
        std::vector<uint8_t> seq;
        uint64_t             cov  = 0;
    };
    /// A thread's emitted unitigs, with bases appended to one arena plus one fixed record each, so there is no per-unitig allocation (at k = 15 a chromosome has tens of millions of unitigs a few dozen bases long).
    struct Rec { uint64_t head, cov, off; uint32_t len, tid; };
    /// Output order is by head k-mer, the primary key.
    /// One American-flag pass radix-buckets the records in place by the top kRadixBits bits of the head, with no second array, and the buckets are sorted independently in parallel.
    /// That costs O(n) moves plus small sorts instead of one O(n log n) sort of all unitigs.
    static constexpr unsigned kRadixBits = 12, kRadixBuckets = 1u << kRadixBits;
    /// Bases are kept 2-bit packed (A,C,G,T = 0..3, four per byte, with Rec::off counting bases), because the arena holds every unitig base until output, 3.4 G bases for a human genome at k = 31, and packing takes a quarter of the bytes.
    struct Out {
        std::vector<uint8_t> bases;
        std::size_t          nbases = 0;
        std::vector<Rec>     recs;
        void add(const Result& r, uint32_t tid) {
            recs.push_back(Rec{r.head, r.cov, nbases, static_cast<uint32_t>(r.seq.size()), tid});
            bases.resize((nbases + r.seq.size() + 3) / 4, 0);
            for (const uint8_t ch : r.seq) {
                const uint8_t c = ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2 : 3;
                bases[nbases >> 2] |= static_cast<uint8_t>(c << ((nbases & 3) * 2));
                ++nbases;
            }
        }
        uint8_t base(std::size_t i) const noexcept { return (bases[i >> 2] >> ((i & 3) * 2)) & 3u; }
    };

    // Linear-time bidirected unitig extraction over the canonical masks.
    //
    // The graph is the node set plus the canonical in/out masks, edge-centric in that a mask bit requires the (k+1)-mer to have been seen in either orientation.
    // The masks are symmetric by construction, so a walk sees each bidirected edge from both ends.
    //
    // Oriented k-mer E -> X is compactable iff E has exactly one successor X, X has exactly one predecessor, and canon(X) != canon(E).
    // Self-loops and hairpins end a unitig, and palindromic k-mers, possible only for even k, are singleton unitigs.
    // An oriented E is a start iff it cannot be compacted to the left.
    // Every non-cycle unitig has exactly two starts (its first k-mer and the reverse complement of its last), is walked from both and is emitted only from the smaller oriented start, so the output is a partition of the node set (each k-mer in exactly one unitig) at O(N) total work.
    // Nodes left unvisited lie on pure cycles and are emitted by a deterministic serial pass.

    /// UnitigView - one emitted unitig, a view into the walk's buffers.
    struct UnitigView { const uint8_t* seq; uint32_t len; uint64_t cov; };

    /// run_masks() - cb(seq, len, coverage) per unitig, in output order.
    template <typename Graph, typename CB>
    static void run_masks(const Graph& g, CB&& cb, unsigned nthreads) {
        run_sorted(g, [&](std::size_t n, const auto& at) {
            for (std::size_t i = 0; i < n; ++i) {
                const UnitigView u = at(i);
                cb(u.seq, u.len, u.cov);
            }
        }, nthreads);
    }

    /// run_sorted() - batch form, calling sink(n, at) once with at(i) the i-th unitig in output order and views valid during the call, so a writer can format them on several threads and still write them in order.
    ///
    /// Graph is DbgGraph, the open-addressing table, or DbgGraphEF with Elias-Fano keys and dense node words, which share one node interface.
    template <typename Graph, typename Sink>
    static void run_sorted(const Graph& g, Sink&& sink, unsigned nthreads) {
        if (g.consumed()) throw std::logic_error("DbgGraph: node words were reused (aux_clear)");
        const unsigned K = g.k();
        const uint64_t kmask = K >= 32 ? ~uint64_t(0)
                                       : ((uint64_t(1) << (2u * K)) - 1);
        PhaseTimer pt;
        // Nodes are the occupied slots of the graph's table.
        // In DbgGraph a node's key, masks, abundance and walker flags share one 64-byte line, so a step costs one probe.
        const std::size_t table_cap = g.cap_;
        const std::size_t N = g.n_;
        if (N == 0) return;
        nthreads = resolve_threads(nthreads);
        if (nthreads == 0) nthreads = 1;

        auto rc = [&](uint64_t v) { return reverse_complement(v, K); };
        auto canon = [&](uint64_t v) { const uint64_t r = rc(v); return v < r ? v : r; };
        auto key_of = [&](std::size_t i) { return g.key(i); };
        constexpr uint8_t kVisited = 1, kClaimed = 2;      // walker flag bits
        auto vis_or = [&](std::size_t i, uint8_t f) { return g.flags_or(i, f); };
        auto vis_get = [&](std::size_t i) { return g.flags(i); };
        const bool even_k = (K % 2) == 0;

        // The compaction test of oriented E at node i comes in two halves.
        // Its local half needs only node i's own slot and checks that E has one successor X, which is neither node i itself, as in a self-loop or hairpin, nor a palindrome, and that node i is not a palindrome.
        // For the remote half, X's slot shows whether X has exactly one predecessor, from the popcount of its in-mask, or of its out-mask when X is the reverse strand.
        // A walk loads X's slot anyway to take the step, so stepping costs no extra probe, and a neighbour is probed only by start tests on unvisited nodes.
        auto cand = [&](uint64_t E, std::size_t i, uint64_t& X) -> bool {
            const uint64_t c = key_of(i);
            if (even_k && c == rc(c)) return false;              // palindromic node
            const bool fwd = E == c;
            const unsigned mk = g.edge_mask(i);
            const unsigned m = fwd ? (mk & 0xFu) : (mk >> 4);
            if (__builtin_popcount(m) != 1) return false;
            const unsigned bb = fwd ? unsigned(__builtin_ctz(m)) : (unsigned(__builtin_ctz(m)) ^ 3u);
            X = ((E << 2) | bb) & kmask;
            if (canon(X) == c) return false;                     // self-loop / hairpin
            if (even_k && X == rc(X)) return false;              // palindrome
            return true;
        };
        auto accept = [&](uint64_t X, std::size_t xi) -> bool {
            const unsigned mk = g.edge_mask(xi);
            return __builtin_popcount(X == canon(X) ? (mk >> 4) : (mk & 0xFu)) == 1;
        };
        auto steps_ok = [&](uint64_t E, std::size_t i, uint64_t& X, std::size_t& xi) -> bool {
            if (!cand(E, i, X)) return false;
            xi = g.rank(canon(X));
            return xi != kNpos && accept(X, xi);
        };
        // Two-sided start test, which probes the predecessor, and its prefetch.
        auto is_start = [&](uint64_t E, std::size_t i) -> bool {
            uint64_t X; std::size_t xi;
            return !steps_ok(rc(E), i, X, xi);
        };
        const unsigned shK = 2u * (K - 1);
        auto prefetch_nbrs = [&](std::size_t i) {
            if (!g.used(i) || (vis_get(i) & kVisited)) return;
            const uint64_t c = key_of(i);
            const unsigned out = g.edge_mask(i) & 0xFu, in = g.edge_mask(i) >> 4;
            if (__builtin_popcount(out) == 1)
                g.prefetch(canon(((c << 2) | unsigned(__builtin_ctz(out))) & kmask));
            if (__builtin_popcount(in) == 1)
                g.prefetch(canon((uint64_t(__builtin_ctz(in)) << shK) | (c >> 2)));
        };
        // Choose the start rule that needs fewer probes on this graph.
        // Both find exactly the same starts, and the choice is global, so every start is found by one rule.
        // The two-sided test probes about 2 neighbours per node, while branch enumeration probes only the successors of branching orientations, which is far fewer on a genome at k = 31 and more on a dense k = 15 graph.
        // The choice compares the two probe counts on a sample of slots.
        bool branch_mode;
        {
            const std::size_t kStartRuleSample = std::min<std::size_t>(table_cap, std::size_t(1) << 18);
            uint64_t two = 0, br = 0;
            for (std::size_t i = 0; i < kStartRuleSample; ++i) {
                if (!g.used(i)) continue;
                const unsigned mk = g.edge_mask(i);
                two += 2;
                const unsigned o = __builtin_popcount(mk & 0xFu), n = __builtin_popcount(mk >> 4);
                br += (o >= 2 ? o : 0) + (n >= 2 ? n : 0);
            }
            branch_mode = br * 4 < two;      // branch probes under a quarter of two-sided ones
        }

        // Prefetch the home slots of the successors of node i's branching orientations, where the start scan probes next.
        auto prefetch_branch = [&](std::size_t i) {
            if (!g.used(i)) return;
            const uint64_t c = key_of(i);
            const unsigned mk = g.edge_mask(i);
            const unsigned out = mk & 0xFu, in = mk >> 4;
            if (__builtin_popcount(out) >= 2)
                for (unsigned t = 0; t < 4; ++t)
                    if (out & (1u << t)) g.prefetch(canon(((c << 2) | t) & kmask));
            if (__builtin_popcount(in) >= 2) {
                const uint64_t r = rc(c);
                for (unsigned t = 0; t < 4; ++t)
                    if (in & (1u << t)) g.prefetch(canon(((r << 2) | (t ^ 3u)) & kmask));
            }
        };

        // Walk right from oriented start E at slot i, appending bases to res, with one probe per step.
        // This is serial and used only for pure cycles.
        // Returns the last oriented k-mer.
        auto walk = [&](uint64_t E, std::size_t i, Result& res,
                        std::size_t* end_slot = nullptr) -> uint64_t {
            const uint64_t start_c = canon(E);
            uint64_t cur = E; std::size_t ci = i;
            res.head = start_c;
            res.seq.clear();
            for (int b = static_cast<int>(K) - 1; b >= 0; --b)
                res.seq.push_back(static_cast<uint8_t>("ACGT"[(E >> (2u * b)) & 3u]));
            res.cov = g.count(i);
            vis_or(i, kVisited);
            for (std::size_t steps = 0; steps < N; ++steps) {
                uint64_t X; std::size_t xi;
                if (!cand(cur, ci, X)) break;
                const uint64_t cx = canon(X);
                if (cx == start_c) break;                   // closed a cycle
                xi = g.rank(cx);
                if (xi == kNpos || !accept(X, xi)) break;
                res.seq.push_back(static_cast<uint8_t>("ACGT"[X & 3u]));
                res.cov += g.count(xi);
                vis_or(xi, kVisited);
                cur = X; ci = xi;
            }
            if (end_slot) *end_slot = ci;
            return cur;
        };

        pt.lap("walk: setup");
        std::vector<Out> per_thread(nthreads);
        // One walk per unitig emits while walking.
        // Both ends may race, and exactly one wins the claim on the end node with the smaller k-mer (a value, not a slot index, since slot order depends on insertion timing).
        // The kept sequence is normalized to start at the smaller of the unitig's two oriented starts, so the output does not depend on which end won.
        auto finish = [&](uint64_t E, std::size_t i, uint64_t T, std::size_t ti,
                          Result& res, Out& out, uint32_t tid) {
            const std::size_t claim = key_of(ti) < key_of(i) ? ti : i;
            if (vis_or(claim, kClaimed) & kClaimed) return;
            if (E > rc(T)) {                           // emitted from the other end
                std::reverse(res.seq.begin(), res.seq.end());
                for (auto& ch : res.seq)
                    ch = ch == 'A' ? 'T' : ch == 'C' ? 'G' : ch == 'G' ? 'C' : 'A';
                res.head = canon(rc(T));
            }
            out.add(res, tid);
        };

        // Interleaved walks.
        // A walk is a chain of dependent probes, since the next k-mer is known only once the current slot is loaded, so one walk per thread waits out one DRAM miss per step.
        // Each thread instead advances kInterleavedWalks walks round-robin.
        // In every round a lane probes the k-mer whose home slot it prefetched in the previous round and prefetches its next one, so that many misses are in flight per thread (group prefetching, Slotin, Algorithms for Modern Hardware).
        {
            constexpr unsigned kInterleavedWalks = 16;
            constexpr std::size_t kScanChunk = 1024;       // slots a thread takes per claim
            constexpr std::size_t kPrefetchAhead = 8;      // slots the start scan prefetches ahead
            std::atomic<std::size_t> next{0};
            struct Lane {
                uint64_t E, cur, X, start_c;
                std::size_t i, ci, steps;
                Result res;
                bool live = false;
            };
            tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                Out& out = per_thread[t];
                Lane lanes[kInterleavedWalks];
                // Start-scan state, the next slot s_i of the claimed range [.., s_hi) and the next orientation s_o of that slot
                std::size_t s_i = 0, s_hi = 0;
                unsigned s_o = 0;
                bool scan_done = false;
                // Next oriented start not already ended from its other end, or false when the table is exhausted.
                //
                // E is a start iff rc(E) cannot step or its one predecessor P has out-degree >= 2.
                // rc(E) cannot step when E's in-degree is not 1 or the step is a self-loop, hairpin or palindrome, all read from E's own slot.
                // The first condition is tested at E's slot and the second at P's, since every branching orientation enumerates its successors and those that pass the local test are starts.
                // The scan therefore probes a neighbour only at branch nodes, about 2 per unitig, rather than twice per node, and finds exactly the starts the plain two-sided test finds.
                std::vector<std::pair<uint64_t, std::size_t>> pend;
                auto next_start = [&](uint64_t& E, std::size_t& i) -> bool {
                    for (;;) {
                        while (!pend.empty()) {
                            const auto p = pend.back();
                            pend.pop_back();
                            if (vis_get(p.second) & kVisited) continue;
                            E = p.first; i = p.second;
                            return true;
                        }
                        if (s_i >= s_hi) {
                            if (scan_done) return false;
                            s_i = next.fetch_add(kScanChunk, std::memory_order_relaxed);
                            if (s_i >= table_cap) { scan_done = true; return false; }
                            s_hi = std::min(table_cap, s_i + kScanChunk);
                            s_o = 0;
                        }
                        if (!branch_mode) {                 // two-sided test
                            if (s_o == 0 && s_i + kPrefetchAhead < table_cap) prefetch_nbrs(s_i + kPrefetchAhead);
                            if (!g.used(s_i) || (vis_get(s_i) & kVisited)) { ++s_i; s_o = 0; continue; }
                            const uint64_t c = key_of(s_i);
                            const uint64_t rcc = rc(c);
                            const unsigned n_or = c == rcc ? 1u : 2u;
                            if (s_o >= n_or) { ++s_i; s_o = 0; continue; }
                            const uint64_t e = s_o == 0 ? c : rcc;
                            ++s_o;
                            if (!is_start(e, s_i)) continue;
                            E = e; i = s_i;
                            return true;
                        }
                        if (s_o == 0 && s_i + kPrefetchAhead < table_cap) prefetch_branch(s_i + kPrefetchAhead);
                        if (!g.used(s_i)) { ++s_i; s_o = 0; continue; }
                        const uint64_t c = key_of(s_i);
                        const uint64_t rcc = rc(c);
                        const unsigned n_or = c == rcc ? 1u : 2u;
                        if (s_o >= n_or) { ++s_i; s_o = 0; continue; }
                        const uint64_t e = s_o == 0 ? c : rcc;
                        ++s_o;
                        // Successors of a branching orientation, visited or not, since their starts are found here
                        {
                            const bool fwd = e == c;
                            const unsigned mk = g.edge_mask(s_i);
                            const unsigned m = fwd ? (mk & 0xFu) : (mk >> 4);
                            if (__builtin_popcount(m) >= 2) {
                                for (unsigned t = 0; t < 4; ++t) {
                                    if (!(m & (1u << t))) continue;
                                    const unsigned bb = fwd ? t : (t ^ 3u);
                                    const uint64_t X = ((e << 2) | bb) & kmask;
                                    const std::size_t xi = g.rank(canon(X));
                                    if (xi == kNpos || (vis_get(xi) & kVisited)) continue;
                                    uint64_t Y;
                                    if (cand(rc(X), xi, Y)) pend.push_back({X, xi});
                                }
                            }
                        }
                        // Visited nodes are never starts, since an interior node compacts on both sides and a walk's own start or end node is handled by that walk.
                        if (vis_get(s_i) & kVisited) continue;
                        uint64_t Y;
                        if (cand(rc(e), s_i, Y)) continue;          // it steps left, so it is not a start
                        E = e; i = s_i;
                        return true;
                    }
                };
                // Lane's next step, computing the successor and prefetching its home slot, or false if the walk ends at cur.
                auto advance = [&](Lane& L) -> bool {
                    if (L.steps >= N || !cand(L.cur, L.ci, L.X)) return false;
                    const uint64_t cx = canon(L.X);
                    if (cx == L.start_c) return false;              // closed a cycle
                    g.prefetch(cx);
                    return true;
                };
                auto start = [&](Lane& L) -> bool {
                    uint64_t E; std::size_t i;
                    while (next_start(E, i)) {
                        L.E = L.cur = E; L.i = L.ci = i; L.steps = 0;
                        L.start_c = canon(E);
                        L.res.head = L.start_c;
                        L.res.seq.clear();
                        for (int b = static_cast<int>(K) - 1; b >= 0; --b)
                            L.res.seq.push_back(static_cast<uint8_t>("ACGT"[(E >> (2u * b)) & 3u]));
                        L.res.cov = g.count(i);
                        vis_or(i, kVisited);
                        if (advance(L)) return L.live = true;
                        finish(L.E, L.i, L.cur, L.ci, L.res, out, t);  // single-node walk
                    }
                    return L.live = false;
                };
                unsigned live = 0;
                for (unsigned l = 0; l < kInterleavedWalks; ++l) live += start(lanes[l]);
                while (live) {
                    for (unsigned l = 0; l < kInterleavedWalks; ++l) {
                        Lane& L = lanes[l];
                        if (!L.live) continue;
                        const std::size_t xi = g.rank(canon(L.X));
                        if (xi == kNpos || !accept(L.X, xi)) {     // X has other predecessors
                            finish(L.E, L.i, L.cur, L.ci, L.res, out, t);
                            if (!start(L)) --live;
                            continue;
                        }
                        L.res.seq.push_back(static_cast<uint8_t>("ACGT"[L.X & 3u]));
                        L.res.cov += g.count(xi);
                        vis_or(xi, kVisited);
                        L.cur = L.X; L.ci = xi; ++L.steps;
                        if (advance(L)) continue;
                        finish(L.E, L.i, L.cur, L.ci, L.res, out, t);
                        if (!start(L)) --live;
                    }
                }
            });
        }
        // Pure cycles, where every node has in = out = 1 and no walk reached them.
        // One parallel sweep collects the unvisited nodes and clears every other node's walker flags, since the table is reusable and unitigs() may be called again.
        // Each cycle is then emitted from its smallest canonical k-mer in canonical orientation.
        // One pass over the cycle finds that k-mer while marking the nodes, and a second one emits.
        std::vector<std::vector<std::size_t>> left(nthreads);
        {
            constexpr std::size_t kSweepChunk = 1u << 16;       // slots per claim
            std::atomic<std::size_t> next{0};
            tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                for (;;) {
                    const std::size_t lo = next.fetch_add(kSweepChunk, std::memory_order_relaxed);
                    if (lo >= table_cap) break;
                    for (std::size_t i = lo; i < std::min(table_cap, lo + kSweepChunk); ++i) {
                        if (!g.used(i)) continue;
                        if (vis_get(i) & kVisited) g.clear_flags(i);
                        else left[t].push_back(i);
                    }
                }
            });
        }
        for (auto& v : left) {
            for (const std::size_t i : v) {
                if (vis_get(i) & kVisited) continue;
                Result tmp;
                uint64_t best = key_of(i);
                std::size_t best_i = i;
                {
                    uint64_t cur = best; std::size_t ci = i;
                    uint64_t X; std::size_t xi;
                    vis_or(i, kVisited);
                    for (std::size_t steps = 0; steps < N && steps_ok(cur, ci, X, xi); ++steps) {
                        if (xi == i) break;
                        vis_or(xi, kVisited);
                        if (canon(X) < best) { best = canon(X); best_i = xi; }
                        cur = X; ci = xi;
                    }
                }
                walk(best, best_i, tmp);
                per_thread[0].add(tmp, 0);
            }
        }
        for (auto& v : left)
            for (const std::size_t i : v) g.clear_flags(i);

        pt.lap("walk: walks");
        // the code order A < C < G < T is the ASCII order of the bases
        auto rec_less = [&](const Rec& a, const Rec& b) {
            if (a.head != b.head) return a.head < b.head;
            if (a.len != b.len) return a.len < b.len;
            const Out& oa = per_thread[a.tid];
            const Out& ob = per_thread[b.tid];
            for (uint32_t j = 0; j < a.len; ++j) {
                const uint8_t x = oa.base(a.off + j), y = ob.base(b.off + j);
                if (x != y) return x < y;
            }
            return false;
        };
        // Heads are canonical k-mers below 4^K, bucketed by their top kRadixBits bits.
        const unsigned rshift = 2 * K > kRadixBits ? 2 * K - kRadixBits : 0;
        const unsigned T = nthreads;
        // (1) Every thread buckets its own records in place with an American-flag pass, all
        //     threads at once, where cnt[t][b] is the size of thread t's bucket b.
        std::vector<std::vector<std::size_t>> cnt(T, std::vector<std::size_t>(kRadixBuckets + 1, 0));
        {
            tpoptoa::run_on_threads(T, [&](unsigned t) {
                auto& v = per_thread[t].recs;
                auto& c = cnt[t];
                for (const Rec& r : v) ++c[(r.head >> rshift) + 1];
                for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt) c[bkt + 1] += c[bkt];   // starts
                std::vector<std::size_t> nxt(c.begin(), c.end() - 1);
                for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt)
                    while (nxt[bkt] < c[bkt + 1]) {
                        const unsigned d = static_cast<unsigned>(v[nxt[bkt]].head >> rshift);
                        if (d == bkt) { ++nxt[bkt]; continue; }
                        std::swap(v[nxt[bkt]], v[nxt[d]++]);
                    }
            });
        }
        // (2) Bucket b's slots for thread t start at boff[b] plus the records of b in earlier threads.
        std::vector<std::size_t> boff(kRadixBuckets + 1, 0);
        for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt) {
            std::size_t n = 0;
            for (unsigned t = 0; t < T; ++t) n += cnt[t][bkt + 1] - cnt[t][bkt];
            boff[bkt + 1] = boff[bkt] + n;
        }
        const std::size_t nrec = boff[kRadixBuckets];
        // (3) copied thread by thread into an uninitialised small-page array (pages become
        //     resident as they fill), with each thread's vector released right after, so at
        //     most one thread's records exist twice
        SmallPageArray<Rec> recs_own(nrec);
        Rec* const recs = recs_own.data();
        {
            std::vector<std::size_t> fill(boff.begin(), boff.end() - 1);
            for (unsigned t = 0; t < T; ++t) {
                auto& v = per_thread[t].recs;
                std::atomic<unsigned> nb{0};
                tpoptoa::run_on_threads(T, [&](unsigned) {
                    for (unsigned bkt; (bkt = nb.fetch_add(1)) < kRadixBuckets;) {
                        const std::size_t lo = cnt[t][bkt], hi = cnt[t][bkt + 1];
                        std::copy(v.begin() + lo, v.begin() + hi, recs + fill[bkt]);
                        fill[bkt] += hi - lo;
                    }
                });
                std::vector<Rec>().swap(v);
            }
        }
        // Buckets are sorted independently.
        // Heads are distinct, as each k-mer lies in one unitig, so the order is total and equals that of one global sort.
        {
            std::atomic<unsigned> nb{0};
            tpoptoa::run_on_threads(nthreads, [&](unsigned) {
                for (unsigned bkt; (bkt = nb.fetch_add(1)) < kRadixBuckets;)
                    std::sort(recs + boff[bkt], recs + boff[bkt + 1], rec_less);
            });
        }
        pt.lap("walk: gather + sort");
        // A view unpacks into a buffer of the calling thread that stays valid until that thread's next at(), and every writer uses a view at once.
        sink(nrec, [&](std::size_t i) {
            const Rec& r = recs[i];
            thread_local std::vector<uint8_t> buf;
            buf.resize(r.len);
            const Out& o = per_thread[r.tid];
            for (uint32_t j = 0; j < r.len; ++j) buf[j] = static_cast<uint8_t>("ACGT"[o.base(r.off + j)]);
            return UnitigView{buf.data(), r.len, r.cov};
        });
        pt.lap("walk: emit (GFA S lines)");
    }
};

// DbgGraph member templates, defined here because they need the extractor.

template <typename CB>
void DbgGraph::unitigs(CB&& cb, unsigned nthreads) const {
    DbgUnitigExtractor::run_masks(*this, std::forward<CB>(cb), nthreads);
}

template <typename Sink>
void DbgGraph::unitigs_sorted(Sink&& sink, unsigned nthreads) const {
    DbgUnitigExtractor::run_sorted(*this, std::forward<Sink>(sink), nthreads);
}

}
