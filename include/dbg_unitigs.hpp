#pragma once

#include "thread_budget.hpp"
#include "phase_timer.hpp"

/**
 * @file
 * @brief Maximal unitig extraction for DbgGraph and DbgGraphEF.
 * dbg_graph.hpp includes it at the very end.
 *
 * The compaction is bidirected, after Conway and Bromage 2011 and Cuttlefish 2 (Khan et al. 2022).
 * Nodes are canonical k-mers with symmetric in/out edge masks, so each edge is visible from both endpoints.
 * An oriented k-mer E extends to X iff E has exactly one successor X, X has exactly one predecessor, and canon(X) != canon(E).
 * Maximal unitigs are the maximal paths of such steps, and run_sorted() spells out the exact rules.
 *
 * @par Parallelism, lock-free
 * Every oriented start is walked independently in O(N) total work.
 * Each unitig is emitted exactly once, from its smaller oriented start, into a per-thread buffer.
 * Nodes are marked with relaxed atomic ORs, and whatever stays unmarked (the pure cycles) is handled serially.
 * Sorting the output by head k-mer makes it identical at every thread count.
 *
 * bench/validate_unitigs.py checks the result against the path, interior and maximality rules, exact L links and duplicate k-mers.
 * It also confirms that the unitig partition equals Cuttlefish 2's --read output at min_count 1.
 */

#include <algorithm>
#include <atomic>
#include <memory>
#include <string>

namespace tpoptoa {

/// An orientation with this many successors or more is a branch, so no unitig continues through it.
inline constexpr unsigned kBranchDegree = 2;

struct DbgUnitigExtractor {
    static constexpr std::size_t kNpos = DbgGraph::kNpos;

    struct Result {
        uint64_t             head = 0;
        std::vector<uint8_t> seq;
        uint64_t             cov  = 0;
    };
    /// Bases go into one arena with a fixed record per unitig, so nothing is allocated per unitig.
    /// That matters because at k = 15 a chromosome has tens of millions of unitigs a few dozen bases long.
    struct Rec { uint64_t head, cov, off; uint32_t len, tid; };
    /// Output is ordered by head k-mer.
    /// One in-place American-flag pass buckets the records by the head's top kRadixBits bits, and the buckets are then sorted in parallel.
    /// That costs O(n) moves plus many small sorts, instead of one O(n log n) sort over every unitig.
    static constexpr unsigned kRadixBits = 12, kRadixBuckets = 1u << kRadixBits;

    /// Bases stay 2-bit packed (four per byte, with Rec::off counting bases) because the arena holds every unitig base until output.
    /// For a human genome at k = 31 that is 3.4 G bases, so packing saves three quarters of the memory.
    struct Out {
        std::vector<uint8_t> bases;
        std::size_t          nbases = 0;
        std::vector<Rec>     recs;

        void add(const Result& unitig, uint32_t tid) {
            recs.push_back(Rec{unitig.head, unitig.cov, nbases, static_cast<uint32_t>(unitig.seq.size()), tid});
            bases.resize((nbases + unitig.seq.size() + kBasesPerByte - 1) / kBasesPerByte, 0);

            for (const uint8_t base_char : unitig.seq) {
                const uint8_t code = base_encode(static_cast<char>(base_char));   // walks only ever emit A/C/G/T

                bases[nbases / kBasesPerByte] |= static_cast<uint8_t>(code << ((nbases % kBasesPerByte) * kBitsPerBase));
                ++nbases;
            }
        }
        uint8_t base(std::size_t base_idx) const noexcept {
            return (bases[base_idx / kBasesPerByte] >> ((base_idx % kBasesPerByte) * kBitsPerBase)) & kBaseMask;
        }
    };

    // Linear-time bidirected unitig extraction over the canonical masks.
    //
    // The graph is the node set plus the canonical in/out masks.
    // It is edge-centric, since a mask bit requires the (k+1)-mer to have been seen in either orientation.
    // Because the masks are symmetric by construction, a walk reaches each bidirected edge from both ends.
    //
    // Oriented k-mer E -> X is compactable iff E has exactly one successor X, X has exactly one predecessor, and canon(X) != canon(E).
    // Self-loops and hairpins end a unitig, while palindromic k-mers (possible only for even k) become singleton unitigs.
    // An oriented E is a start iff it can't be compacted to the left.
    // Every non-cycle unitig has exactly two starts: its first k-mer and the reverse complement of its last.
    // Both are walked, but only the smaller start emits, so each k-mer lands in exactly one unitig at O(N) total work.
    // Whatever stays unvisited lies on pure cycles, which a deterministic serial pass emits afterwards.

    /// Points into the walk's buffers rather than owning a copy.
    struct UnitigView { const uint8_t* seq; uint32_t len; uint64_t cov; };

    template <typename Graph, typename CB>
    static void run_masks(const Graph& graph, CB&& callback, unsigned nthreads) {
        run_sorted(graph, [&](std::size_t count, const auto& unitig_at) {
            for (std::size_t i = 0; i < count; ++i) {
                const UnitigView view = unitig_at(i);

                callback(view.seq, view.len, view.cov);
            }
        }, nthreads);
    }

    /// Calls sink(count, unitig_at) once, where unitig_at(i) returns the i-th unitig in output order.
    /// The views stay valid for the whole call, so a writer can format them on several threads and still write them in order.
    ///
    /// Graph may be DbgGraph (the open-addressing table) or DbgGraphEF (Elias-Fano keys), since both share one node interface.
    template <typename Graph, typename Sink>
    static void run_sorted(const Graph& graph, Sink&& sink, unsigned nthreads) {
        if (graph.consumed()) throw std::logic_error("DbgGraph: node words were reused (aux_clear)");

        const unsigned kmer_len = graph.kmer_len();
        const uint64_t kmask = kmer_mask(kmer_len);
        PhaseTimer timer;
        // Nodes are the occupied slots of the graph's table.
        // In DbgGraph a node's key, masks, abundance and walker flags share one 64-byte line, so each step costs a single probe.
        const std::size_t table_cap = graph.cap_;
        const std::size_t num_nodes = graph.num_nodes_;

        if (num_nodes == 0) return;

        nthreads = resolve_threads(nthreads);

        if (nthreads == 0) nthreads = 1;

        auto rev_comp = [&](uint64_t kmer) { return reverse_complement(kmer, kmer_len); };
        auto canon = [&](uint64_t kmer) { const uint64_t rev = rev_comp(kmer); return kmer < rev ? kmer : rev; };
        auto key_of = [&](std::size_t slot) { return graph.key(slot); };
        constexpr uint8_t kVisited = 1, kClaimed = 2;      // walker flag bits
        auto vis_or = [&](std::size_t slot, uint8_t flag_bits) { return graph.flags_or(slot, flag_bits); };
        auto vis_get = [&](std::size_t slot) { return graph.flags(slot); };
        const bool even_k = (kmer_len & 1u) == 0;

        // The compaction test of oriented E at node i comes in two halves.
        // The local half reads only node i's slot.
        // It checks that E has exactly one successor X, that X isn't node i itself (a self-loop or hairpin), and that neither is a palindrome.
        // The remote half reads X's slot to see whether X has exactly one predecessor, using its in-mask, or its out-mask when X is on the reverse strand.
        // A walk loads X's slot anyway to take the step, so stepping costs no extra probe.
        // Neighbours are therefore probed only by start tests on unvisited nodes.
        auto cand = [&](uint64_t oriented, std::size_t slot, uint64_t& succ) -> bool {
            const uint64_t canon_kmer = key_of(slot);

            if (even_k && canon_kmer == rev_comp(canon_kmer)) return false;              // palindromic node

            const bool fwd = oriented == canon_kmer;
            const unsigned edge_bits = graph.edge_mask(slot);
            const unsigned dir_mask = fwd ? (edge_bits & kOutEdgeMask) : (edge_bits >> kOutEdgeBits);

            if (__builtin_popcount(dir_mask) != 1) return false;

            const unsigned succ_base = fwd ? unsigned(__builtin_ctz(dir_mask)) : (unsigned(__builtin_ctz(dir_mask)) ^ kComplementMask);

            succ = ((oriented << kBitsPerBase) | succ_base) & kmask;

            if (canon(succ) == canon_kmer) return false;                     // self-loop / hairpin

            if (even_k && succ == rev_comp(succ)) return false;              // palindrome

            return true;
        };
        auto accept = [&](uint64_t succ, std::size_t succ_slot) -> bool {
            const unsigned edge_bits = graph.edge_mask(succ_slot);

            return __builtin_popcount(succ == canon(succ) ? (edge_bits >> kOutEdgeBits) : (edge_bits & kOutEdgeMask)) == 1;
        };
        auto steps_ok = [&](uint64_t oriented, std::size_t slot, uint64_t& succ, std::size_t& succ_slot) -> bool {
            if (!cand(oriented, slot, succ)) return false;

            succ_slot = graph.rank(canon(succ));
            return succ_slot != kNpos && accept(succ, succ_slot);
        };
        // The two-sided start test probes the predecessor, and this is its prefetch.
        auto is_start = [&](uint64_t oriented, std::size_t slot) -> bool {
            uint64_t succ; std::size_t succ_slot;

            return !steps_ok(rev_comp(oriented), slot, succ, succ_slot);
        };
        const unsigned first_base_shift = kBitsPerBase * (kmer_len - 1);
        auto prefetch_nbrs = [&](std::size_t slot) {
            if (!graph.used(slot) || (vis_get(slot) & kVisited)) return;

            const uint64_t canon_kmer = key_of(slot);
            const unsigned out = graph.edge_mask(slot) & kOutEdgeMask, in_edges = graph.edge_mask(slot) >> kOutEdgeBits;

            if (__builtin_popcount(out) == 1)
                graph.prefetch(canon(((canon_kmer << kBitsPerBase) | unsigned(__builtin_ctz(out))) & kmask));

            if (__builtin_popcount(in_edges) == 1)
                graph.prefetch(canon((uint64_t(__builtin_ctz(in_edges)) << first_base_shift) | (canon_kmer >> kBitsPerBase)));
        };
        // Pick whichever start rule needs fewer probes on this graph.
        // Both find exactly the same starts, and the choice is global, so each start is found by one rule.
        // The two-sided test probes about 2 neighbours per node.
        // Branch enumeration, by contrast, probes only the successors of branching orientations.
        // That is far fewer on a genome at k = 31, yet more on a dense k = 15 graph, so the decision compares probe counts on a sample of slots.
        bool branch_mode;

        {
            constexpr std::size_t kStartRuleSampleCap = std::size_t(1) << 18;
            const std::size_t kStartRuleSample = std::min<std::size_t>(table_cap, kStartRuleSampleCap);
            uint64_t two = 0, branch_probes = 0;

            for (std::size_t i = 0; i < kStartRuleSample; ++i) {
                if (!graph.used(i)) continue;

                const unsigned edge_bits = graph.edge_mask(i);

                two += kStrands;
                const unsigned out_deg = __builtin_popcount(edge_bits & kOutEdgeMask), in_deg = __builtin_popcount(edge_bits >> kOutEdgeBits);

                branch_probes += (out_deg >= kBranchDegree ? out_deg : 0) + (in_deg >= kBranchDegree ? in_deg : 0);
            }

            constexpr int kBranchProbeShareDiv = 4;      // branch probes must be under 1/4 of the two-sided ones

            branch_mode = branch_probes * kBranchProbeShareDiv < two;
        }

        // The start scan probes the successors of node i's branching orientations next, so their home slots are prefetched.
        auto prefetch_branch = [&](std::size_t slot) {
            if (!graph.used(slot)) return;

            const uint64_t canon_kmer = key_of(slot);
            const unsigned edge_bits = graph.edge_mask(slot);
            const unsigned out = edge_bits & kOutEdgeMask, in_edges = edge_bits >> kOutEdgeBits;

            if (__builtin_popcount(out) >= kBranchDegree)
                for (unsigned t = 0; t < kNumBases; ++t)
                    if (out & (1u << t)) graph.prefetch(canon(((canon_kmer << kBitsPerBase) | t) & kmask));

            if (__builtin_popcount(in_edges) >= kBranchDegree) {
                const uint64_t rev = rev_comp(canon_kmer);

                for (unsigned t = 0; t < kNumBases; ++t)
                    if (in_edges & (1u << t)) graph.prefetch(canon(((rev << kBitsPerBase) | (t ^ kComplementMask)) & kmask));
            }
        };

        // Walks right from oriented start E at slot i, appending bases to res, and returns the last oriented k-mer.
        // Only the serial pure-cycle pass uses it.
        auto walk = [&](uint64_t start_kmer, std::size_t slot, Result& res,
                        std::size_t* end_slot = nullptr) -> uint64_t {
            const uint64_t start_c = canon(start_kmer);
            uint64_t cur = start_kmer; std::size_t cur_slot = slot;

            res.head = start_c;
            res.seq.clear();

            for (int b = static_cast<int>(kmer_len) - 1; b >= 0; --b)
                res.seq.push_back(static_cast<uint8_t>(kBaseChars[(start_kmer >> (kBitsPerBase * b)) & kBaseMask]));

            res.cov = graph.count(slot);
            vis_or(slot, kVisited);

            for (std::size_t steps = 0; steps < num_nodes; ++steps) {
                uint64_t succ; std::size_t succ_slot;

                if (!cand(cur, cur_slot, succ)) break;

                const uint64_t canon_succ = canon(succ);

                if (canon_succ == start_c) break;                   // closed a cycle

                succ_slot = graph.rank(canon_succ);

                if (succ_slot == kNpos || !accept(succ, succ_slot)) break;

                res.seq.push_back(static_cast<uint8_t>(kBaseChars[succ & kBaseMask]));
                res.cov += graph.count(succ_slot);
                vis_or(succ_slot, kVisited);
                cur = succ; cur_slot = succ_slot;
            }

            if (end_slot) *end_slot = cur_slot;

            return cur;
        };

        timer.lap("walk: setup");
        std::vector<Out> per_thread(nthreads);
        // Each walk emits as it goes.
        // Both ends may race, and only the walk that claims the end node with the smaller k-mer emits the unitig.
        // Comparing k-mer values rather than slot indices matters, because slot order depends on insertion timing.
        // The kept sequence is then normalized to begin at the smaller oriented start, so the output never depends on which end claimed it first.
        auto finish = [&](uint64_t start_kmer, std::size_t start_slot, uint64_t tail_kmer, std::size_t tail_slot,
                          Result& res, Out& out, uint32_t tid) {
            const std::size_t claim = key_of(tail_slot) < key_of(start_slot) ? tail_slot : start_slot;

            if (vis_or(claim, kClaimed) & kClaimed) return;

            if (start_kmer > rev_comp(tail_kmer)) {                           // emitted from the other end
                std::reverse(res.seq.begin(), res.seq.end());

                for (auto& base_char : res.seq)
                    base_char = base_char == 'A' ? 'T' : base_char == 'C' ? 'G' : base_char == 'G' ? 'C' : 'A';

                res.head = canon(rev_comp(tail_kmer));
            }

            out.add(res, tid);
        };

        // Interleaved walks.
        // A walk is a chain of dependent probes, because the next k-mer is only known once the current slot loads.
        // A single walk per thread would therefore wait out one DRAM miss per step.
        // Instead each thread advances kInterleavedWalks walks round-robin.
        // In every round a lane probes the k-mer it prefetched last round and prefetches its next one, keeping many misses in flight (group prefetching, Slotin).
        {
            constexpr unsigned kInterleavedWalks = 16;
            constexpr std::size_t kScanChunk = 1024;       // slots a thread takes per claim
            constexpr std::size_t kPrefetchAhead = 8;      // slots the start scan prefetches ahead
            std::atomic<std::size_t> next{0};

            struct Lane {
                uint64_t start_kmer, cur, succ, start_c;
                std::size_t start_slot, cur_slot, steps;
                Result res;
                bool live = false;
            };
            tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                Out& out = per_thread[thread_idx];
                Lane lanes[kInterleavedWalks];
                // The start scan remembers the next slot s_i in its claimed range [.., s_hi) and the next orientation s_o of that slot.
                std::size_t s_i = 0, s_hi = 0;
                unsigned s_o = 0;
                bool scan_done = false;
                // Finds the next oriented start that wasn't already finished from its other end, or returns false once the table is exhausted.
                //
                // E is a start iff rc(E) can't step or its single predecessor P has out-degree >= 2.
                // rc(E) can't step when E's in-degree isn't 1 or the step is a self-loop, hairpin or palindrome, all visible in E's own slot.
                // The first condition is tested at E's slot and the second at P's, because each branching orientation enumerates its successors.
                // Successors that pass the local test are starts.
                // As a result, the scan probes neighbours only at branch nodes (about 2 per unitig) instead of twice per node, yet finds the same starts.
                std::vector<std::pair<uint64_t, std::size_t>> pend;
                auto next_start = [&](uint64_t& start_kmer, std::size_t& slot) -> bool {
                    for (;;) {
                        // Starts queued at an earlier branch node come first, and any that a walk has visited since are skipped.
                        while (!pend.empty()) {
                            const auto pending = pend.back();

                            pend.pop_back();

                            if (vis_get(pending.second) & kVisited) continue;

                            start_kmer = pending.first; slot = pending.second;
                            return true;
                        }

                        // Once its chunk is exhausted, the thread claims the next kScanChunk slots from the shared cursor.
                        if (s_i >= s_hi) {
                            if (scan_done) return false;

                            s_i = next.fetch_add(kScanChunk, std::memory_order_relaxed);

                            if (s_i >= table_cap) { scan_done = true; return false; }

                            s_hi = std::min(table_cap, s_i + kScanChunk);
                            s_o = 0;
                        }

                        // Without branch mode, each orientation of every unvisited slot is tested directly with is_start().
                        if (!branch_mode) {
                            if (s_o == 0 && s_i + kPrefetchAhead < table_cap) prefetch_nbrs(s_i + kPrefetchAhead);

                            if (!graph.used(s_i) || (vis_get(s_i) & kVisited)) { ++s_i; s_o = 0; continue; }

                            const uint64_t canon_kmer = key_of(s_i);
                            const uint64_t rcc = rev_comp(canon_kmer);
                            const unsigned n_or = canon_kmer == rcc ? 1u : kStrands;

                            if (s_o >= n_or) { ++s_i; s_o = 0; continue; }

                            const uint64_t oriented = s_o == 0 ? canon_kmer : rcc;

                            ++s_o;

                            if (!is_start(oriented, s_i)) continue;

                            start_kmer = oriented; slot = s_i;
                            return true;
                        }

                        // In branch mode, only branching orientations are expanded, and their successors are queued as starts.
                        if (s_o == 0 && s_i + kPrefetchAhead < table_cap) prefetch_branch(s_i + kPrefetchAhead);

                        if (!graph.used(s_i)) { ++s_i; s_o = 0; continue; }

                        const uint64_t canon_kmer = key_of(s_i);
                        const uint64_t rcc = rev_comp(canon_kmer);
                        const unsigned n_or = canon_kmer == rcc ? 1u : kStrands;

                        if (s_o >= n_or) { ++s_i; s_o = 0; continue; }

                        const uint64_t oriented = s_o == 0 ? canon_kmer : rcc;

                        ++s_o;
                        // Enumerate a branching orientation's successors whether visited or not, because this is where their starts get found.
                        {
                            const bool fwd = oriented == canon_kmer;
                            const unsigned edge_bits = graph.edge_mask(s_i);
                            const unsigned dir_mask = fwd ? (edge_bits & kOutEdgeMask) : (edge_bits >> kOutEdgeBits);

                            if (__builtin_popcount(dir_mask) >= kBranchDegree) {
                                for (unsigned t = 0; t < kNumBases; ++t) {
                                    if (!(dir_mask & (1u << t))) continue;

                                    const unsigned succ_base = fwd ? t : (t ^ kComplementMask);
                                    const uint64_t succ = ((oriented << kBitsPerBase) | succ_base) & kmask;
                                    const std::size_t succ_slot = graph.rank(canon(succ));

                                    if (succ_slot == kNpos || (vis_get(succ_slot) & kVisited)) continue;

                                    uint64_t pred;

                                    if (cand(rev_comp(succ), succ_slot, pred)) pend.push_back({succ, succ_slot});
                                }
                            }
                        }
                        // Visited nodes are never starts: interior nodes compact on both sides, and a walk handles its own start and end nodes.
                        if (vis_get(s_i) & kVisited) continue;

                        uint64_t pred;

                        if (cand(rev_comp(oriented), s_i, pred)) continue;          // it steps left, so it is not a start

                        start_kmer = oriented; slot = s_i;
                        return true;
                    }
                };
                // Computes the lane's successor and prefetches its home slot, or returns false when the walk ends at cur.
                auto advance = [&](Lane& lane) -> bool {
                    if (lane.steps >= num_nodes || !cand(lane.cur, lane.cur_slot, lane.succ)) return false;

                    const uint64_t canon_succ = canon(lane.succ);

                    if (canon_succ == lane.start_c) return false;              // closed a cycle

                    graph.prefetch(canon_succ);
                    return true;
                };
                auto start = [&](Lane& lane) -> bool {
                    uint64_t start_kmer; std::size_t slot;

                    while (next_start(start_kmer, slot)) {
                        lane.start_kmer = lane.cur = start_kmer; lane.start_slot = lane.cur_slot = slot; lane.steps = 0;
                        lane.start_c = canon(start_kmer);
                        lane.res.head = lane.start_c;
                        lane.res.seq.clear();

                        for (int b = static_cast<int>(kmer_len) - 1; b >= 0; --b)
                            lane.res.seq.push_back(static_cast<uint8_t>(kBaseChars[(start_kmer >> (kBitsPerBase * b)) & kBaseMask]));

                        lane.res.cov = graph.count(slot);
                        vis_or(slot, kVisited);

                        if (advance(lane)) return lane.live = true;

                        finish(lane.start_kmer, lane.start_slot, lane.cur, lane.cur_slot, lane.res, out, thread_idx);  // single-node walk
                    }

                    return lane.live = false;
                };
                unsigned live = 0;

                for (unsigned l = 0; l < kInterleavedWalks; ++l) live += start(lanes[l]);

                while (live) {
                    for (unsigned l = 0; l < kInterleavedWalks; ++l) {
                        Lane& lane = lanes[l];

                        if (!lane.live) continue;

                        const std::size_t succ_slot = graph.rank(canon(lane.succ));

                        if (succ_slot == kNpos || !accept(lane.succ, succ_slot)) {     // X has other predecessors
                            finish(lane.start_kmer, lane.start_slot, lane.cur, lane.cur_slot, lane.res, out, thread_idx);

                            if (!start(lane)) --live;

                            continue;
                        }

                        lane.res.seq.push_back(static_cast<uint8_t>(kBaseChars[lane.succ & kBaseMask]));
                        lane.res.cov += graph.count(succ_slot);
                        vis_or(succ_slot, kVisited);
                        lane.cur = lane.succ; lane.cur_slot = succ_slot; ++lane.steps;

                        if (advance(lane)) continue;

                        finish(lane.start_kmer, lane.start_slot, lane.cur, lane.cur_slot, lane.res, out, thread_idx);

                        if (!start(lane)) --live;
                    }
                }
            });
        }
        // Pure cycles, where every node has in = out = 1 and no walk ever reached them.
        // One parallel sweep collects the unvisited nodes and clears every other node's walker flags, since unitigs() may run again on the same table.
        // Each cycle is then emitted from its smallest canonical k-mer, in canonical orientation.
        // The first pass around a cycle finds that k-mer while marking nodes, and a second pass emits.
        std::vector<std::vector<std::size_t>> left(nthreads);

        {
            constexpr std::size_t kSweepChunk = 1u << 16;       // slots per claim
            std::atomic<std::size_t> next{0};

            tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                for (;;) {
                    const std::size_t chunk_start = next.fetch_add(kSweepChunk, std::memory_order_relaxed);

                    if (chunk_start >= table_cap) break;

                    for (std::size_t i = chunk_start; i < std::min(table_cap, chunk_start + kSweepChunk); ++i) {
                        if (!graph.used(i)) continue;

                        if (vis_get(i) & kVisited) graph.clear_flags(i);
                        else left[thread_idx].push_back(i);
                    }
                }
            });
        }
        // Cycles are rare and short next to the unitig walks, so they are handled on one thread, which also keeps their output order fixed.
        for (auto& thread_left : left) {
            for (const std::size_t slot : thread_left) {
                if (vis_get(slot) & kVisited) continue;

                Result tmp;
                uint64_t best = key_of(slot);
                std::size_t best_i = slot;

                {
                    uint64_t cur = best; std::size_t cur_slot = slot;
                    uint64_t succ; std::size_t succ_slot;

                    vis_or(slot, kVisited);

                    for (std::size_t steps = 0; steps < num_nodes && steps_ok(cur, cur_slot, succ, succ_slot); ++steps) {
                        if (succ_slot == slot) break;

                        vis_or(succ_slot, kVisited);

                        if (canon(succ) < best) { best = canon(succ); best_i = succ_slot; }

                        cur = succ; cur_slot = succ_slot;
                    }
                }
                // Starting at the smallest canonical k-mer fixes the cycle's rotation.
                walk(best, best_i, tmp);
                per_thread[0].add(tmp, 0);
            }
        }

        // The cycle nodes were skipped by the sweep's flag reset, so they are cleared here for the next call.
        for (auto& thread_left : left)
            for (const std::size_t slot : thread_left) graph.clear_flags(slot);

        timer.lap("walk: walks");
        // The code order A < C < G < T is also the ASCII order of the bases.
        auto rec_less = [&](const Rec& lhs, const Rec& rhs) {
            if (lhs.head != rhs.head) return lhs.head < rhs.head;

            if (lhs.len != rhs.len) return lhs.len < rhs.len;

            const Out& lhs_out = per_thread[lhs.tid];
            const Out& rhs_out = per_thread[rhs.tid];

            for (uint32_t j = 0; j < lhs.len; ++j) {
                const uint8_t lhs_base = lhs_out.base(lhs.off + j), rhs_base = rhs_out.base(rhs.off + j);

                if (lhs_base != rhs_base) return lhs_base < rhs_base;
            }

            return false;
        };
        // Heads are canonical k-mers below 4^K, bucketed by their top kRadixBits bits.
        const unsigned head_bits = kBitsPerBase * kmer_len;
        const unsigned rshift = head_bits > kRadixBits ? head_bits - kRadixBits : 0;
        const unsigned num_threads = nthreads;
        // (1) Every thread buckets its own records in place with an American-flag pass, and cnt[t][b] becomes the size of thread t's bucket b.
        std::vector<std::vector<std::size_t>> cnt(num_threads, std::vector<std::size_t>(kRadixBuckets + 1, 0));

        {
            tpoptoa::run_on_threads(num_threads, [&](unsigned thread_idx) {
                auto& thread_recs = per_thread[thread_idx].recs;
                auto& counts = cnt[thread_idx];

                for (const Rec& rec : thread_recs) ++counts[(rec.head >> rshift) + 1];

                for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt) counts[bkt + 1] += counts[bkt];   // starts

                std::vector<std::size_t> nxt(counts.begin(), counts.end() - 1);

                for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt)
                    while (nxt[bkt] < counts[bkt + 1]) {
                        const unsigned dest_bkt = static_cast<unsigned>(thread_recs[nxt[bkt]].head >> rshift);

                        if (dest_bkt == bkt) { ++nxt[bkt]; continue; }

                        std::swap(thread_recs[nxt[bkt]], thread_recs[nxt[dest_bkt]++]);
                    }
            });
        }
        // (2) Thread t's share of bucket b starts at boff[b] plus the records earlier threads put in b.
        std::vector<std::size_t> boff(kRadixBuckets + 1, 0);

        for (unsigned bkt = 0; bkt < kRadixBuckets; ++bkt) {
            std::size_t bucket_size = 0;

            for (unsigned t = 0; t < num_threads; ++t) bucket_size += cnt[t][bkt + 1] - cnt[t][bkt];

            boff[bkt + 1] = boff[bkt] + bucket_size;
        }

        const std::size_t nrec = boff[kRadixBuckets];
        // (3) Records are copied thread by thread into an uninitialised small-page array, whose pages become resident only as they fill.
        //     Each thread's vector is released right after its copy, so at most one thread's records ever exist twice.
        SmallPageArray<Rec> recs_own(nrec);
        Rec* const recs = recs_own.data();

        {
            std::vector<std::size_t> fill(boff.begin(), boff.end() - 1);

            for (unsigned t = 0; t < num_threads; ++t) {
                auto& thread_recs = per_thread[t].recs;
                std::atomic<unsigned> next_bkt{0};

                tpoptoa::run_on_threads(num_threads, [&](unsigned) {
                    for (unsigned bkt; (bkt = next_bkt.fetch_add(1)) < kRadixBuckets;) {
                        const std::size_t range_lo = cnt[t][bkt], range_hi = cnt[t][bkt + 1];

                        std::copy(thread_recs.begin() + range_lo, thread_recs.begin() + range_hi, recs + fill[bkt]);
                        fill[bkt] += range_hi - range_lo;
                    }
                });
                std::vector<Rec>().swap(thread_recs);
            }
        }
        // Buckets are sorted independently.
        // Since each k-mer lies in one unitig, heads are distinct, so the order is total and matches one global sort.
        {
            std::atomic<unsigned> next_bkt{0};

            tpoptoa::run_on_threads(nthreads, [&](unsigned) {
                for (unsigned bkt; (bkt = next_bkt.fetch_add(1)) < kRadixBuckets;)
                    std::sort(recs + boff[bkt], recs + boff[bkt + 1], rec_less);
            });
        }
        timer.lap("walk: gather + sort");
        // A view unpacks into the calling thread's buffer, which stays valid until that thread calls unitig_at() again.
        // Writers consume each view immediately, so that lifetime is enough.
        sink(nrec, [&](std::size_t rec_idx) {
            const Rec& rec = recs[rec_idx];
            thread_local std::vector<uint8_t> buf;

            buf.resize(rec.len);
            const Out& owner = per_thread[rec.tid];

            for (uint32_t j = 0; j < rec.len; ++j) buf[j] = static_cast<uint8_t>(kBaseChars[owner.base(rec.off + j)]);

            return UnitigView{buf.data(), rec.len, rec.cov};
        });
        timer.lap("walk: emit (GFA S lines)");
    }
};

// DbgGraph's member templates are defined here because they need the extractor.

template <typename CB>
void DbgGraph::unitigs(CB&& cb, unsigned nthreads) const {
    DbgUnitigExtractor::run_masks(*this, std::forward<CB>(cb), nthreads);
}

template <typename Sink>
void DbgGraph::unitigs_sorted(Sink&& sink, unsigned nthreads) const {
    DbgUnitigExtractor::run_sorted(*this, std::forward<Sink>(sink), nthreads);
}

}
