#pragma once

/*
 * dbg_part.hpp - partitioned compaction, which builds maximal unitigs without a global k-mer table.
 *
 * BCALM (Chikhi et al. 2014) compacts inside minimizer buckets and glues the pieces across buckets.
 * optik follows that idea in its own form, where the pieces are walked in each counted partition's cache-resident table and glued by matching the canonical (k+1)-mer keys of their open ends.
 *
 * The counter gives every canonical k-mer to exactly one minimizer partition, and consecutive k-mers that share a minimizer share a partition, so most edges of a unitig never leave their partition.
 * Once a partition is counted, its open-addressing table, sized for the L2 cache, holds every k-mer of it with its abundance and edges.
 * PartBuilder::partition() walks the partition's pieces of unitigs, called fragments, there on exact masks with one probe of a cache-resident table per step.
 * It applies the walker's compaction rule from dbg_unitigs.hpp, where E -> X compacts iff E has one successor X, X has one predecessor, and the step is no self-loop, hairpin or palindrome.
 *
 * A fragment end whose one neighbour lies in another partition is open.
 * The end's own half of the compaction test passed, out-degree 1 on the right or in-degree 1 on the left, and the other half belongs to the neighbour.
 * An open end is keyed by the canonical (k+1)-mer of its edge, which both sides compute alike, so two fragments are glued exactly when both sides of the edge are open, which is the full compaction test.
 * End keys are bucketed, spilled past the budget, and matched.
 * Chains of glued fragments are followed from their free ends, which pure cycles lack, and every unitig is assembled 2-bit packed from the fragment stream.
 *
 * The output equals the walker's byte for byte.
 * Each unitig starts at the smaller of its two oriented end k-mers, whose canonical form is the head, pure cycles start at their smallest canonical k-mer in canonical orientation, unitigs are sorted by head, KC is the sum of abundances, and the GFA link masks are the end k-mers' exact masks.
 *
 * With -c > 1 an edge between partitions is real only if both k-mers are solid, which neither partition can see alone.
 * In the counter's first pass every solid k-mer claims its cross-partition edges under the same keys.
 * An edge claimed once has a non-solid end, and its bit is cleared when the second pass reaches the claimer's partition.
 *
 * The fragment stream and end records spill like the counter's buffers.
 * In memory there is one partition table per thread, then two 32-bit words per fragment while chains are followed, holding the partner ends and later the unitig id and offset.
 * The assembly reads those words in fragment order, so past kLinkSpill fragments they are streamed back from the spill directory instead of being held next to the packed unitigs.
 * The unitigs are then sorted through an index, with heads read from the bases, and never copied.
 * Requires k <= 31, so that the edge key, a (k+1)-mer, fits 64 bits.
 */

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <malloc.h>
#include <vector>

#include "hugepage.hpp"
#include "kmer.hpp"
#include "phase_timer.hpp"
#include "spill.hpp"
#include "wyhash.hpp"
#include "thread_budget.hpp"

namespace tpoptoa {

/// PartNode - a solid k-mer of one partition, with its canonical key, abundance and in << 4 | out edge mask (exact inside the partition, while PartBuilder checks the cross-partition bits).
struct PartNode { uint64_t key; uint32_t cnt; uint32_t edge_mask; };

/// PartUnitigs - the unitigs of a partitioned compaction, sorted by head k-mer.
class PartUnitigs {
public:
    struct View { const uint8_t* seq; uint32_t len; uint64_t cov; };

    unsigned    k() const noexcept { return k_; }
    std::size_t core_size() const noexcept { return nkmers_; }
    std::size_t n_unitigs() const noexcept { return len_.size(); }
    bool        consumed() const noexcept { return false; }
    /// exit_mask() - unitig @i's exit masks, with the successors of its last k-mer in the low nibble and those of the reverse complement of its first in the high nibble.
    uint8_t     exit_mask(std::size_t i) const noexcept { return exit_masks_[order_[i]]; }

    /// unitigs_sorted() - calls sink(n, at), where at(i) is the i-th unitig in output order, a view into a buffer of the calling thread that stays valid until that thread's next at().
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned = 0) const {
        sink(len_.size(), [&](std::size_t i) {
            const std::size_t u = order_[i];
            // Output order is not storage order, so the fields of a unitig a few places ahead are fetched now, and its bases a few places later, once its offset is in cache.
            constexpr std::size_t kFieldsAhead = 16, kBasesAhead = 8;
            if (i + kFieldsAhead < order_.size()) {
                const uint32_t f = order_[i + kFieldsAhead];
                __builtin_prefetch(&off_[f]); __builtin_prefetch(&len_[f]); __builtin_prefetch(&cov_[f]);
            }
            if (i + kBasesAhead < order_.size()) __builtin_prefetch(bases_.data() + (off_[order_[i + kBasesAhead]] >> 2));
            const uint32_t len = len_[u];
            const uint64_t off = off_[u];
            thread_local std::vector<uint8_t> buf;
            buf.resize(len);
            // Four bases per byte, with unitigs starting byte-aligned, and lut.c[byte] spells the byte's four bases.
            // Built at compile time, so no first-call initialization check per unitig.
            struct Lut {
                char c[256][4] = {};
                constexpr Lut() {
                    constexpr char kBases[] = "ACGT";
                    for (int b = 0; b < 256; ++b)
                        for (int q = 0; q < 4; ++q) c[b][q] = kBases[(b >> (2 * q)) & 3];
                }
            };
            static constexpr Lut lut{};
            const uint8_t* src = bases_.data() + (off >> 2);
            uint32_t j = 0;
            for (; j + 4 <= len; j += 4) std::memcpy(buf.data() + j, lut.c[src[j >> 2]], 4);
            for (; j < len; ++j) buf[j] = static_cast<uint8_t>(lut.c[src[j >> 2]][j & 3]);
            return View{buf.data(), len, cov_[u]};
        });
    }
    // There is no node table.
    // The GFA writer takes the exit masks from exit_mask() and finds link targets through its own index of segment ends, and these members keep its graph-mode code compiling.
    static constexpr std::size_t kNpos = ~std::size_t(0);
    std::size_t rank(uint64_t) const noexcept { return kNpos; }
    void prefetch(uint64_t) const noexcept {}
    unsigned out_mask(uint64_t) const noexcept { return 0; }
    unsigned in_mask(uint64_t) const noexcept { return 0; }
    void aux_clear(unsigned) noexcept {}
    void aux_set(std::size_t, uint32_t) noexcept {}
    uint32_t aux(std::size_t) const noexcept { return 0; }

    template <typename CB>
    void unitigs(CB&& cb, unsigned nthreads = 0) const {
        unitigs_sorted([&](std::size_t n, const auto& at) {
            for (std::size_t i = 0; i < n; ++i) { const View u = at(i); cb(u.seq, u.len, u.cov); }
        }, nthreads);
    }

private:
    friend class PartBuilder;
    unsigned k_ = 0;
    uint64_t nkmers_ = 0;
    // Unitig u's bases sit at off_[u] (byte-aligned) and its fields are in assembly order, while order_ is the output order (by head k-mer).
    // The bases and fields are not copied into output order, since that copy was the compaction's memory peak.
    HugeVec<uint8_t> bases_;
    std::vector<uint32_t> order_;             // output position -> unitig
    std::vector<uint64_t> off_, cov_;
    std::vector<uint32_t> len_;
    std::vector<uint8_t> exit_masks_;
};

/// Smallest k at which --graph auto chooses partitioned compaction.
inline constexpr unsigned kPartModeMinK = 19;

/// part_mode_wanted() - whether a build at @k uses partitioned compaction.
///
/// The order of precedence is OPTIK_DBG_MODE=part|table, then the CLI's --graph part|table read through count_limits().graph_mode, then auto from kPartModeMinK up.
/// Never for k > 31.
///
/// The threshold comes from hg38 chr1 at 8 threads, measured at commit 9b6e0c6 with identical GFA either way.
///   k = 31  5.5 s / 1.00 GB  vs the table's 6.8 s / 3.84 GB
///   k = 21  7.3 s / 1.22 GB  vs 7.1 s / 3.74 GB
///   k = 15  21 s / 5.2 GB    vs 19 s / 5.3 GB (60M short unitigs)
inline bool part_mode_wanted(unsigned k) {
    if (k > 31) return false;
    if (const char* e = std::getenv("OPTIK_DBG_MODE")) return std::strcmp(e, "part") == 0;
    const int m = count_limits().graph_mode;
    if (m != 0) return m == 1;
    return k >= kPartModeMinK;
}

/// PartBuilder - builds PartUnitigs from the counter's partitions.
///
/// With @filter, claim() runs during the first counting pass and match_claims() after it.
/// partition() runs once per counted partition in the second pass, then finish().
class PartBuilder {
public:
    static constexpr uint32_t kNone = ~uint32_t(0);

    /// @K <= 31.
    /// @filter: min_count > 1, so cross-partition edges are claimed and checked.
    PartBuilder(unsigned K, unsigned nthreads, std::size_t nparts, bool filter,
                uint64_t per_thread_bytes, std::string dir)
        : k_(K), nthreads_(nthreads), nparts_(nparts), filter_(filter), budget_(per_thread_bytes), dir_(std::move(dir)),
          kmask_(kmer_mask(K)), frag_writers_(nthreads), scratch_(nthreads) {
        if (K > 31) throw std::invalid_argument("partitioned compaction needs k <= 31");
        ends_ = std::make_unique<spill::BucketSpill<EndRec>>(nthreads_, kEndBuckets, budget_, dir_);
        if (filter_) claims_ = std::make_unique<spill::BucketSpill<ClaimRec>>(nthreads_, kEndBuckets, budget_, dir_);
    }

    /// mask_edge() - the canonical (k+1)-mer of the edge of canonical k-mer
    /// @c through mask bit @bit (bits 0-3 are out bases and 4-7 in bases, in the canonical frame).
    uint64_t mask_edge(uint64_t c, unsigned bit) const noexcept {
        const uint64_t w = bit < 4 ? ((c << 2) | bit) : ((uint64_t(bit - 4) << (2 * k_)) | c);
        return canon1(w);
    }

    /// claim() - first counting pass, with filter only.
    /// A solid k-mer of partition @part has the cross-partition edge @edge_key.
    void claim(unsigned tid, uint32_t part, uint64_t edge_key) {
        claims_->add(tid, bucket(edge_key), ClaimRec{edge_key, part});
    }
    /// match_claims() - after the first pass, edges claimed only once have a non-solid end and go to their claimer's per-partition invalid list.
    void match_claims() {
        invalid_edges_ = std::make_unique<spill::BucketSpill<uint64_t>>(nthreads_, nparts_, budget_, dir_);
        std::atomic<std::size_t> nb{0};
        std::exception_ptr err;
        std::mutex mu;
        tpoptoa::run_on_threads(nthreads_, [&](unsigned t) {
            try {
                std::vector<ClaimRec> v;
                for (std::size_t b; (b = nb.fetch_add(1)) < kEndBuckets;) {
                    claims_->load(b, v);
                    std::sort(v.begin(), v.end(), [](const ClaimRec& x, const ClaimRec& y) {
                        const uint64_t a = x.key, c = y.key; return a < c; });
                    for (std::size_t i = 0; i < v.size();) {
                        std::size_t j = i + 1;
                        while (j < v.size() && v[j].key == v[i].key) ++j;
                        if (j - i == 1) invalid_edges_->add(t, v[i].part, v[i].key);
                        i = j;
                    }
                    std::vector<ClaimRec>().swap(v);
                }
            } catch (...) { std::lock_guard<std::mutex> l(mu); if (!err) err = std::current_exception(); }
        });
        if (err) std::rethrow_exception(err);
        claims_.reset();
    }

    /// partition() - second pass, walking the fragments of counted partition
    /// @part, whose solid k-mers are @nodes.
    void partition(unsigned tid, uint32_t part, std::vector<PartNode>& nodes);

    /// finish() - glue, assemble, orient and sort the fragments into the unitigs.
    std::unique_ptr<PartUnitigs> finish(PhaseTimer& pt);

private:
    static constexpr std::size_t kEndBuckets = 4096;              // end / claim buckets
    static constexpr uint64_t kBucketSeed = 0x51ED270B0A3C4E1DULL;
    static constexpr uint64_t kFibonacciMultiplier = 0x9E3779B97F4A7C15ULL;
    /// Fragment record header, a u32 k-mer count and a u64 coverage.
    static constexpr std::size_t kFragHeaderBytes = 12;
    /// FragWriter::nk value of a fragment whose k-mer count is in `big`.
    static constexpr uint16_t kNkInBig = 0xFFFFu;
    /// Link word 1 after chain walking, where bit 31 says the fragment is reversed in its unitig and bits 0..30 hold its base offset there.
    static constexpr uint32_t kReversedBit = 0x80000000u, kOffsetMask = 0x7FFFFFFFu;
    /// Chain walks advanced round-robin per thread (see finish()).
    static constexpr unsigned kInterleavedChains = 16;
    /// Fragments past which link words are streamed from disk in finish().
    static constexpr uint64_t kLinkSpill = uint64_t(1) << 24;
    /// Smallest per-partition index size, in slots.
    static constexpr std::size_t kMinIndexSlots = 16;
    struct __attribute__((packed)) EndRec { uint64_t key; uint32_t local; uint16_t tid; uint8_t side; };
    struct __attribute__((packed)) ClaimRec { uint64_t key; uint32_t part; };
    /// FragWriter - one thread's fragments, records of [u32 nk][u64 cov][packed bases] kept in memory up to the budget and then appended to its spill file in chunks.
    struct FragWriter {
        std::vector<uint8_t> mem;
        spill::File file;
        std::vector<std::pair<uint64_t, uint64_t>> chunks;
        std::vector<uint16_t> nk;              // k-mers per fragment, where kNkInBig points to `big`
        std::vector<std::pair<uint32_t, uint32_t>> big;   // local index and k-mers of fragments of at least kNkInBig
        std::vector<uint8_t> masks;            // exit masks, low nibble right and high nibble left
    };
    struct Scratch {
        std::vector<uint32_t> idx;             // open addressing, holding node + 1
        std::vector<uint8_t> vis;
        std::vector<uint8_t> seq;              // the fragment's base codes
        std::vector<uint64_t> inv;
    };

    std::size_t bucket(uint64_t key) const noexcept { return mix64(key, kBucketSeed) & (kEndBuckets - 1); }
    uint64_t rc(uint64_t v) const noexcept { return reverse_complement(v, k_); }
    uint64_t canon(uint64_t v) const noexcept { const uint64_t r = rc(v); return v < r ? v : r; }
    /// canon1() - canonical form of the (k+1)-mer @w.
    uint64_t canon1(uint64_t w) const noexcept {
        const uint64_t r = reverse_complement(w, k_ + 1); return w < r ? w : r;
    }

    void emit(unsigned tid, const std::vector<uint8_t>& seq, uint64_t cov, int open_l, uint64_t key_l,
              int open_r, uint64_t key_r, uint8_t masks);

    unsigned k_, nthreads_;
    std::size_t nparts_;
    bool filter_;
    uint64_t budget_;
    std::string dir_;
    uint64_t kmask_;
    std::vector<FragWriter> frag_writers_;
    std::vector<Scratch> scratch_;
    std::unique_ptr<spill::BucketSpill<EndRec>> ends_;
    std::unique_ptr<spill::BucketSpill<ClaimRec>> claims_;
    std::unique_ptr<spill::BucketSpill<uint64_t>> invalid_edges_;
};

/// emit() - store one fragment, with 2-bit codes @seq, coverage @cov and exit
/// @masks, together with its open ends' keys, spilling the thread's buffer past the budget.
inline void PartBuilder::emit(unsigned tid, const std::vector<uint8_t>& seq, uint64_t cov, int open_l,
                              uint64_t key_l, int open_r, uint64_t key_r, uint8_t masks) {
    FragWriter& w = frag_writers_[tid];
    const uint32_t local = static_cast<uint32_t>(w.nk.size());
    const uint32_t nk = static_cast<uint32_t>(seq.size() - k_ + 1);
    // Fragments are runs inside one partition and almost always short.
    w.nk.push_back(static_cast<uint16_t>(std::min<uint32_t>(nk, kNkInBig)));
    if (nk >= kNkInBig) w.big.push_back({local, nk});
    w.masks.push_back(masks);
    const std::size_t o = w.mem.size();
    const std::size_t nbytes = (seq.size() + 3) / 4;
    w.mem.resize(o + kFragHeaderBytes + nbytes, 0);
    std::memcpy(w.mem.data() + o, &nk, 4);
    std::memcpy(w.mem.data() + o + 4, &cov, 8);
    uint8_t* d = w.mem.data() + o + kFragHeaderBytes;
    for (std::size_t i = 0; i < seq.size(); ++i) d[i >> 2] |= static_cast<uint8_t>(seq[i] << ((i & 3) * 2));
    if (open_l) ends_->add(tid, bucket(key_l), EndRec{key_l, local, static_cast<uint16_t>(tid), 0});
    if (open_r) ends_->add(tid, bucket(key_r), EndRec{key_r, local, static_cast<uint16_t>(tid), 1});
    if (w.mem.size() > budget_) {
        w.file.open(dir_);
        std::vector<struct iovec> iov{{w.mem.data(), w.mem.size()}};
        w.chunks.push_back({w.file.size(), w.mem.size()});
        w.file.append(iov);
        w.mem.clear();
    }
}

inline void PartBuilder::partition(unsigned tid, uint32_t part, std::vector<PartNode>& nodes) {
    const std::size_t n = nodes.size();
    if (n == 0) return;
    Scratch& scr = scratch_[tid];
    const unsigned K = k_;
    const uint64_t kmask = kmask_;
    const bool even_k = (K % 2) == 0;
    // local index of the partition's solid nodes (cache-resident)
    std::size_t cap = kMinIndexSlots;
    while (cap < 2 * n) cap <<= 1;
    unsigned log2_cap = 0;
    while ((std::size_t(1) << log2_cap) < cap) ++log2_cap;
    scr.idx.assign(cap, 0);
    auto slot = [&](uint64_t key) {
        return static_cast<std::size_t>((key * kFibonacciMultiplier) >> (64 - log2_cap));
    };
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t h = slot(nodes[i].key);
        while (scr.idx[h]) h = (h + 1) & (cap - 1);
        scr.idx[h] = static_cast<uint32_t>(i + 1);
    }
    auto find = [&](uint64_t c) -> std::size_t {
        for (std::size_t h = slot(c);; h = (h + 1) & (cap - 1)) {
            const uint32_t v = scr.idx[h];
            if (!v) return SIZE_MAX;
            if (nodes[v - 1].key == c) return v - 1;
        }
    };
    const unsigned shK = 2u * (K - 1);
    auto nbr = [&](uint64_t c, unsigned bit) -> uint64_t {          // neighbour's canonical k-mer
        return canon(bit < 4 ? (((c << 2) | bit) & kmask) : ((uint64_t(bit - 4) << shK) | (c >> 2)));
    };
    // With -c > 1, cross-partition edges claimed once have a non-solid end.
    if (filter_) {
        invalid_edges_->load(part, scr.inv);
        if (!scr.inv.empty()) {
            std::sort(scr.inv.begin(), scr.inv.end());
            for (auto& nd : nodes)
                for (unsigned b = 0; b < 8; ++b)
                    if ((nd.edge_mask & (1u << b)) && find(nbr(nd.key, b)) == SIZE_MAX
                        && std::binary_search(scr.inv.begin(), scr.inv.end(), mask_edge(nd.key, b)))
                        nd.edge_mask &= ~(1u << b);
        }
    }
    // the walker's compaction test (dbg_unitigs.hpp), on node i
    auto cand = [&](uint64_t E, std::size_t i, uint64_t& X, unsigned& bb) -> bool {
        const uint64_t c = nodes[i].key;
        if (even_k && c == rc(c)) return false;
        const bool fwd = E == c;
        const unsigned mk = nodes[i].edge_mask;
        const unsigned m = fwd ? (mk & 0xFu) : ((mk >> 4) & 0xFu);
        if (__builtin_popcount(m) != 1) return false;
        bb = fwd ? unsigned(__builtin_ctz(m)) : (unsigned(__builtin_ctz(m)) ^ 3u);
        X = ((E << 2) | bb) & kmask;
        if (canon(X) == c) return false;
        if (even_k && X == rc(X)) return false;
        return true;
    };
    auto accept = [&](uint64_t X, std::size_t xi) -> bool {
        const unsigned mk = nodes[xi].edge_mask;
        return __builtin_popcount(X == canon(X) ? ((mk >> 4) & 0xFu) : (mk & 0xFu)) == 1;
    };
    // Successors of oriented k-mer E at node i, i.e. its exit mask
    auto succ = [&](uint64_t E, std::size_t i) -> unsigned {
        const uint64_t c = nodes[i].key;
        const unsigned mk = nodes[i].edge_mask;
        unsigned m = 0;
        if (E == c) m |= mk & 0xFu;
        if (rc(E) == c) {
            const unsigned in = (mk >> 4) & 0xFu;
            for (unsigned b = 0; b < 4; ++b) if (in & (1u << (b ^ 3u))) m |= 1u << b;
        }
        return m;
    };
    // internal step E -> X with both here and compactable
    auto step_in = [&](uint64_t E, std::size_t i, uint64_t& X, std::size_t& xi) -> bool {
        unsigned bb;
        if (!cand(E, i, X, bb)) return false;
        xi = find(canon(X));
        return xi != SIZE_MAX && accept(X, xi);
    };
    scr.vis.assign(n, 0);
    auto walk = [&](uint64_t E, std::size_t i, bool cycle) {
        // Fragment from oriented start E at node i, extended right as far as it compacts inside the partition
        scr.seq.clear();
        for (int b = static_cast<int>(K) - 1; b >= 0; --b) scr.seq.push_back(static_cast<uint8_t>((E >> (2u * b)) & 3u));
        uint64_t cov = nodes[i].cnt;
        scr.vis[i] = 1;
        const uint64_t start_c = canon(E);
        uint64_t cur = E; std::size_t ci = i;
        int open_r = 0; uint64_t key_r = 0;
        for (;;) {
            uint64_t X; unsigned bb;
            if (!cand(cur, ci, X, bb)) break;
            const uint64_t cx = canon(X);
            if (cx == start_c) break;                         // closed a (local) cycle
            const std::size_t xi = find(cx);
            if (xi == SIZE_MAX) {                             // the neighbour is elsewhere
                if (!cycle) { open_r = 1; key_r = canon1((cur << 2) | bb); }
                break;
            }
            if (!accept(X, xi)) break;
            scr.seq.push_back(static_cast<uint8_t>(X & 3u));
            cov += nodes[xi].cnt;
            scr.vis[xi] = 1;
            cur = X; ci = xi;
        }
        // The left end steps as rc(E).
        int open_l = 0; uint64_t key_l = 0;
        if (!cycle) {
            uint64_t Y; unsigned bb;
            if (cand(rc(E), i, Y, bb) && find(canon(Y)) == SIZE_MAX) { open_l = 1; key_l = canon1((rc(E) << 2) | bb); }
        }
        const uint8_t masks = static_cast<uint8_t>(succ(cur, ci) | (succ(rc(E), i) << 4));
        emit(tid, scr.seq, cov, open_l, key_l, open_r, key_r, masks);
    };
    // Fragment starts are oriented k-mers with no internal compactable predecessor.
    for (std::size_t i = 0; i < n; ++i) {
        const uint64_t c = nodes[i].key;
        const uint64_t rcc = rc(c);
        for (unsigned o = 0; o < (c == rcc ? 1u : 2u); ++o) {
            if (scr.vis[i]) break;
            const uint64_t E = o == 0 ? c : rcc;
            uint64_t Y; std::size_t yi;
            if (step_in(rc(E), i, Y, yi)) continue;           // compacts to the left here
            walk(E, i, false);
        }
    }
    // The rest lie on cycles inside the partition, the graph's pure cycles, and are emitted from their smallest canonical k-mer in canonical orientation, as the walker does.
    for (std::size_t i = 0; i < n; ++i) {
        if (scr.vis[i]) continue;
        uint64_t best = nodes[i].key; std::size_t best_i = i;
        uint64_t cur = nodes[i].key; std::size_t ci = i;
        for (std::size_t steps = 0; steps < n; ++steps) {
            uint64_t X; std::size_t xi;
            if (!step_in(cur, ci, X, xi) || xi == i) break;
            if (canon(X) < best) { best = canon(X); best_i = xi; }
            cur = X; ci = xi;
        }
        walk(best, best_i, true);                             // marks the whole cycle
    }
}

inline std::unique_ptr<PartUnitigs> PartBuilder::finish(PhaseTimer& pt) {
    const unsigned K = k_;
    const unsigned T = nthreads_;
    // Fragment ids of thread t follow those of threads < t.
    std::vector<uint64_t> base(T + 1, 0);
    for (unsigned t = 0; t < T; ++t) base[t + 1] = base[t] + frag_writers_[t].nk.size();
    const uint64_t nfrags = base[T];
    if (nfrags >= (uint64_t(1) << 31)) throw std::length_error("partitioned compaction: 2^31 fragments or more");
    HugeVec<uint16_t> frag_kmers16(nfrags);
    std::vector<std::pair<uint64_t, uint32_t>> big;            // (fragment, k-mers), sorted
    HugeVec<uint8_t> frag_masks(nfrags);
    for (unsigned t = 0; t < T; ++t) {
        std::copy(frag_writers_[t].nk.begin(), frag_writers_[t].nk.end(), frag_kmers16.begin() + static_cast<std::ptrdiff_t>(base[t]));
        for (const auto& b : frag_writers_[t].big) big.push_back({base[t] + b.first, b.second});
        std::vector<std::pair<uint32_t, uint32_t>>().swap(frag_writers_[t].big);
        std::copy(frag_writers_[t].masks.begin(), frag_writers_[t].masks.end(), frag_masks.begin() + static_cast<std::ptrdiff_t>(base[t]));
        std::vector<uint16_t>().swap(frag_writers_[t].nk);
        std::vector<uint8_t>().swap(frag_writers_[t].masks);
    }
    auto nk = [&](uint64_t f) -> uint64_t {
        const uint16_t v = frag_kmers16[f];
        if (v != kNkInBig) return v;
        return std::lower_bound(big.begin(), big.end(), std::pair<uint64_t, uint32_t>{f, 0})->second;
    };
    // Partner end of each fragment end (2f left, 2f+1 right)
    HugeVec<uint32_t> link(2 * nfrags, kNone);
    auto par = [&](auto&& fn) {
        std::exception_ptr err;
        std::mutex mu;
        tpoptoa::run_on_threads(T, [&](unsigned t) {
            try { fn(t); } catch (...) { std::lock_guard<std::mutex> l(mu); if (!err) err = std::current_exception(); }
        });
        if (err) std::rethrow_exception(err);
    };
    {
        std::atomic<std::size_t> nb{0};
        par([&](unsigned) {
            std::vector<EndRec> v;
            for (std::size_t b; (b = nb.fetch_add(1)) < kEndBuckets;) {
                if (b + T < kEndBuckets) ends_->will_need(b + T);
                ends_->load(b, v);
                std::sort(v.begin(), v.end(), [](const EndRec& x, const EndRec& y) {
                    const uint64_t a = x.key, c = y.key; return a < c; });
                for (std::size_t i = 0; i < v.size();) {
                    std::size_t j = i + 1;
                    while (j < v.size() && v[j].key == v[i].key) ++j;
                    if (j - i > 2) throw std::logic_error("partitioned compaction: an edge with three open ends");
                    if (j - i == 2) {
                        const uint32_t a = static_cast<uint32_t>((base[v[i].tid] + v[i].local) * 2 + v[i].side);
                        const uint32_t c = static_cast<uint32_t>((base[v[i + 1].tid] + v[i + 1].local) * 2 + v[i + 1].side);
                        link[a] = c; link[c] = a;
                    }
                    i = j;
                }
                std::vector<EndRec>().swap(v);
            }
        });
        ends_.reset();
    }
    malloc_trim(0);                  // the end-key buckets' freed pages
    if (pt.on) std::fprintf(stderr, "[timing] part: %llu fragments\n", static_cast<unsigned long long>(nfrags));
    pt.lap("part: glue (end keys matched)");

    // Chains are followed from free ends, where an end is entered and the fragment's other end exits, following partners until an end has none.
    // A chain is found from both of its free ends and kept from the smaller.
    // Its links are dependent loads, so a thread advances kInterleavedChains walks round-robin, and every walk prefetches its next link word one round ahead as the graph walker does, putting many misses in flight instead of one.
    std::vector<std::vector<uint32_t>> starts(T);
    {
        const uint64_t NE = 2 * nfrags;
        par([&](unsigned t) {
            const uint64_t lo = NE * t / T, hi = NE * (t + 1) / T;
            uint64_t next = lo;
            struct Lane { uint32_t start, cur; uint64_t steps; bool live; } L[kInterleavedChains];
            auto start = [&](Lane& l) -> bool {
                while (next < hi) {
                    const uint64_t e = next++;
                    if (link[e] != kNone) continue;
                    l.start = l.cur = static_cast<uint32_t>(e); l.steps = 0;
                    __builtin_prefetch(&link[e ^ 1u]);
                    return l.live = true;
                }
                return l.live = false;
            };
            unsigned live = 0;
            for (auto& l : L) live += start(l);
            while (live) {
                for (auto& l : L) {
                    if (!l.live) continue;
                    const uint32_t exit = l.cur ^ 1u;
                    const uint32_t nx = link[exit];
                    if (nx == kNone) {
                        if (l.start <= exit) starts[t].push_back(l.start);
                        if (!start(l)) --live;
                        continue;
                    }
                    if (++l.steps > nfrags) throw std::logic_error("partitioned compaction: a chain from a free end does not end");
                    l.cur = nx;
                    __builtin_prefetch(&link[nx ^ 1u]);
                }
            }
        });
    }
    // Unitig ids.
    // Each chain is walked again, and each of its fragments' link words is overwritten with the unitig id and the base offset ORed with the reversed bit.
    // A fragment belongs to one chain, so no other walk reads them.
    std::vector<uint64_t> sbase(T + 1, 0);
    for (unsigned t = 0; t < T; ++t) sbase[t + 1] = sbase[t] + starts[t].size();
    const uint64_t nchains = sbase[T];
    std::vector<uint32_t> ulen(nchains);                 // bases (< 2^31, checked)
    std::vector<uint8_t> unitig_exit_masks(nchains);
    std::vector<uint8_t> done((nfrags + 7) / 8, 0);
    par([&](unsigned t) {
        struct Lane { uint64_t uid, off; uint32_t cur; uint8_t mstart; bool live; } L[kInterleavedChains];
        std::size_t next = 0;
        auto start = [&](Lane& l) -> bool {
            if (next >= starts[t].size()) return l.live = false;
            l.uid = sbase[t] + next;
            l.cur = starts[t][next++];
            l.off = 0;
            const uint32_t f0 = l.cur >> 1;
            // Exiting through the entered end at the start gives the unitig's left exit.
            l.mstart = (l.cur & 1u) ? (frag_masks[f0] & 0xFu) : (frag_masks[f0] >> 4);
            __builtin_prefetch(&link[l.cur & ~1u]);
            return l.live = true;
        };
        unsigned live = 0;
        for (auto& l : L) live += start(l);
        while (live) {
            for (auto& l : L) {
                if (!l.live) continue;
                const uint32_t f = l.cur >> 1;
                const bool rev = l.cur & 1u;                   // entered at its right end
                const uint32_t exit = l.cur ^ 1u;
                const uint32_t nx = link[exit];
                link[2 * f] = static_cast<uint32_t>(l.uid);
                link[2 * f + 1] = static_cast<uint32_t>(l.off) | (rev ? kReversedBit : 0u);
                __atomic_fetch_or(&done[f >> 3], static_cast<uint8_t>(1u << (f & 7)), __ATOMIC_RELAXED);
                l.off += nk(f);
                if (nx == kNone) {
                    const uint8_t mend = (exit & 1u) ? (frag_masks[f] & 0xFu) : (frag_masks[f] >> 4);
                    unitig_exit_masks[l.uid] = static_cast<uint8_t>(mend | (l.mstart << 4));
                    if (l.off + K - 1 >= kReversedBit)
                        throw std::length_error("partitioned compaction: a unitig of 2^31 bases");
                    ulen[l.uid] = static_cast<uint32_t>(l.off + K - 1);
                    if (!start(l)) --live;
                    continue;
                }
                l.cur = nx;
                __builtin_prefetch(&link[nx & ~1u]);
                __builtin_prefetch(&frag_kmers16[nx >> 1]);
            }
        }
    });
    // Pure cycles, where every end is glued and no walk reached them, are rare and run serially.
    std::vector<uint64_t> cycles;
    for (uint64_t f = 0; f < nfrags; ++f) {
        if (done[f >> 3] & (1u << (f & 7))) continue;
        const uint64_t uid = ulen.size();
        uint64_t off = 0;
        uint32_t cur = static_cast<uint32_t>(2 * f);             // enter at the left end
        for (uint64_t steps = 0; steps <= nfrags; ++steps) {
            const uint32_t g = cur >> 1;
            const bool rev = cur & 1u;
            const uint32_t nx = link[cur ^ 1u];
            if (done[g >> 3] & (1u << (g & 7))) break;           // closed the cycle
            link[2 * g] = static_cast<uint32_t>(uid);
            link[2 * g + 1] = static_cast<uint32_t>(off) | (rev ? kReversedBit : 0u);
            done[g >> 3] |= static_cast<uint8_t>(1u << (g & 7));
            off += nk(g);
            if (nx == kNone) throw std::logic_error("partitioned compaction: an open chain left after the walks");
            cur = nx;
        }
        if (off + K - 1 >= kReversedBit) throw std::length_error("partitioned compaction: a unitig of 2^31 bases");
        ulen.push_back(static_cast<uint32_t>(off + K - 1));
        unitig_exit_masks.push_back(0);
        cycles.push_back(uid);
    }
    const uint64_t U = ulen.size();
    if (U >= (uint64_t(1) << 32)) throw std::length_error("partitioned compaction: 2^32 unitigs or more");
    HugeVec<uint16_t>().swap(frag_kmers16);
    HugeVec<uint8_t>().swap(frag_masks);
    std::vector<uint8_t>().swap(done);
    for (auto& v : starts) std::vector<uint32_t>().swap(v);
    malloc_trim(0);
    pt.lap("part: chains");

    // The link words now hold each fragment's (unitig, offset), and the assembly reads them in fragment order, the order of the fragment stream.
    // Past kLinkSpill fragments they are written to the spill directory and read back in slices beside the stream, so the 8 bytes per fragment are not held in memory with the bases.
    constexpr uint64_t kSlice = uint64_t(1) << 18;               // fragments per read
    spill::File link_file;
    const bool links_spilled = nfrags >= kLinkSpill || std::getenv("OPTIK_SPILL_TEST");   // tests spill at any size
    if (links_spilled) {
        link_file.open(dir_);
        std::vector<struct iovec> iov{{link.data(), 2 * nfrags * sizeof(uint32_t)}};
        link_file.append(iov);
        HugeVec<uint32_t>().swap(link);
        malloc_trim(0);
    }

    // Assembly writes unitig u's bases at uoff[u] (byte-aligned, so no two unitigs share a byte) and each fragment's bases at its offset, oriented.
    // Glued fragments overlap by k - 1 bases of equal value, which are written with atomic OR.
    auto out = std::make_unique<PartUnitigs>();
    out->k_ = K;
    std::vector<uint64_t> uoff(U + 1, 0);
    for (uint64_t u = 0; u < U; ++u) uoff[u + 1] = uoff[u] + ((ulen[u] + 3) & ~uint64_t(3));
    out->bases_.assign((uoff[U] / 4 + 7) & ~uint64_t(7), 0);       // whole 64-bit words
    uint8_t* arena = out->bases_.data();
    uint64_t* arena_words = reinterpret_cast<uint64_t*>(arena);
    std::vector<uint64_t> cov(U, 0);
    par([&](unsigned t) {
        FragWriter& w = frag_writers_[t];
        uint64_t j = 0;
        std::vector<uint8_t> buf;
        // fragment f's two link words, from memory or from the current slice
        std::vector<uint32_t> lw;
        uint64_t lw_lo = 0, lw_hi = 0;
        auto words = [&](uint64_t f) -> const uint32_t* {
            if (!links_spilled) return &link[2 * f];
            if (f >= lw_hi) {
                lw_lo = f;
                lw_hi = std::min(f + kSlice, base[t + 1]);
                lw.resize(2 * (lw_hi - lw_lo));
                link_file.read(2 * lw_lo * sizeof(uint32_t), lw.size() * sizeof(uint32_t),
                           reinterpret_cast<uint8_t*>(lw.data()));
            }
            return &lw[2 * (f - lw_lo)];
        };
        auto parse = [&](const uint8_t* p, std::size_t len) {
            std::size_t o = 0;
            while (o < len) {
                uint32_t n; uint64_t c;
                std::memcpy(&n, p + o, 4);
                std::memcpy(&c, p + o + 4, 8);
                const uint32_t L = n + K - 1;
                const uint8_t* d = p + o + kFragHeaderBytes;
                const uint64_t f = base[t] + j;
                const uint32_t* lk = words(f);
                const uint64_t uid = lk[0];
                const uint32_t w1 = lk[1];
                const bool rev = w1 & kReversedBit;
                const uint64_t at = uoff[uid] + (w1 & kOffsetMask);
                __atomic_fetch_add(&cov[uid], c, __ATOMIC_RELAXED);
                // 32 bases per 64-bit word take one atomic OR per word (glued fragments share their k - 1 overlapping bases, with equal values).
                uint64_t acc = 0, pos = at;
                for (uint32_t q = 0; q < L; ++q, ++pos) {
                    const uint32_t src = rev ? L - 1 - q : q;
                    uint64_t code = (d[src >> 2] >> ((src & 3) * 2)) & 3u;
                    if (rev) code ^= 3u;
                    acc |= code << ((pos & 31) * 2);
                    if ((pos & 31) == 31) {
                        if (acc) __atomic_fetch_or(&arena_words[pos >> 5], acc, __ATOMIC_RELAXED);
                        acc = 0;
                    }
                }
                if (acc) __atomic_fetch_or(&arena_words[(pos - 1) >> 5], acc, __ATOMIC_RELAXED);
                o += kFragHeaderBytes + (L + 3) / 4;
                ++j;
            }
        };
        for (const auto& ch : w.chunks) {
            buf.resize(ch.second);
            w.file.read(ch.first, ch.second, buf.data());
            parse(buf.data(), buf.size());
        }
        parse(w.mem.data(), w.mem.size());
        std::vector<uint8_t>().swap(w.mem);
    });
    frag_writers_.clear();
    HugeVec<uint32_t>().swap(link);
    if (links_spilled) link_file.punch(0, link_file.size());
    pt.lap("part: assemble");

    // orientation (the smaller oriented end k-mer first), head, cycles
    auto get = [&](uint64_t pos) -> uint8_t { return (arena[pos >> 2] >> ((pos & 3) * 2)) & 3u; };
    auto put = [&](uint64_t pos, uint8_t c) {
        uint8_t& b = arena[pos >> 2];
        b = static_cast<uint8_t>((b & ~(3u << ((pos & 3) * 2))) | (c << ((pos & 3) * 2)));
    };
    auto kmer_at = [&](uint64_t pos) {
        uint64_t v = 0;
        for (unsigned q = 0; q < K; ++q) v = (v << 2) | get(pos + q);
        return v;
    };
    par([&](unsigned t) {
        // Cycle ids are appended in increasing order, so the list is sorted.
        auto cy = std::lower_bound(cycles.begin(), cycles.end(), U * t / T);
        for (uint64_t u = U * t / T; u < U * (t + 1) / T; ++u) {
            const uint64_t a = uoff[u], L = ulen[u];
            if (cy != cycles.end() && *cy == u) {
                ++cy;
                // A circular sequence of n = L - K + 1 bases.
                // The walker starts at the smallest canonical k-mer in canonical orientation.
                const uint64_t n = L - K + 1;
                std::vector<uint8_t> c(n);
                for (uint64_t i = 0; i < n; ++i) c[i] = get(a + i);
                uint64_t fw = 0;
                for (unsigned q = 0; q < K; ++q) fw = (fw << 2) | c[q % n];
                uint64_t best = ~uint64_t(0), bi = 0; bool bfwd = true;
                for (uint64_t i = 0; i < n; ++i) {
                    if (i) fw = ((fw << 2) | c[(i + K - 1) % n]) & kmask_;
                    const uint64_t r = rc(fw);
                    const uint64_t cn = fw < r ? fw : r;
                    if (cn < best) { best = cn; bi = i; bfwd = fw <= r; }
                }
                // Forward reads c[bi], c[bi+1], ..., and reversed reads the complement walking back from the k-mer's last base.
                for (uint64_t q = 0; q < L; ++q) {
                    const uint8_t code = bfwd ? c[(bi + q) % n]
                                              : static_cast<uint8_t>(c[(bi + K - 1 + n * (L / n + 2) - q) % n] ^ 3u);
                    put(a + q, code);
                }
                const uint8_t first_succ = get(a + K - 1);           // the cycle continues into its start
                const uint8_t pre = get(a + L - K);                   // the last k-mer's first base
                unitig_exit_masks[u] = static_cast<uint8_t>((1u << first_succ) | ((1u << (pre ^ 3u)) << 4));
                continue;
            }
            const uint64_t first_kmer = kmer_at(a), last_kmer = kmer_at(a + L - K);
            if (first_kmer > rc(last_kmer)) {                                     // emitted from the other end
                for (uint64_t i = 0, j = L - 1; i < j; ++i, --j) {
                    const uint8_t x = get(a + i), y = get(a + j);
                    put(a + i, y ^ 3u); put(a + j, x ^ 3u);
                }
                if (L & 1) put(a + L / 2, get(a + L / 2) ^ 3u);
                unitig_exit_masks[u] = static_cast<uint8_t>((unitig_exit_masks[u] >> 4) | (unitig_exit_masks[u] << 4));
            }
        }
    });
    // Sorted by head.
    // Oriented, a unitig's head is the canonical form of its first k-mer (a path starts at its smaller oriented end k-mer, a cycle at its smallest canonical k-mer), so a head is read from the bases where needed instead of being kept per unitig, through an index bucketed by the heads' top bits, whose buckets are sorted in parallel with their heads gathered just for them.
    auto head_of = [&](uint64_t u) { return canon(kmer_at(uoff[u])); };
    uint64_t nkm = 0;
    for (uint64_t u = 0; u < U; ++u) nkm += ulen[u] - K + 1;
    out->nkmers_ = nkm;
    {
        const unsigned kRadixBits = 12, shift = 2 * K > kRadixBits ? 2 * K - kRadixBits : 0;
        const std::size_t NBK = std::size_t(1) << std::min(kRadixBits, 2 * K);
        std::vector<uint16_t> top(U);
        std::vector<std::vector<uint64_t>> tc(T, std::vector<uint64_t>(NBK, 0));
        par([&](unsigned t) {
            for (uint64_t u = U * t / T; u < U * (t + 1) / T; ++u) {
                top[u] = static_cast<uint16_t>(head_of(u) >> shift);
                ++tc[t][top[u]];
            }
        });
        std::vector<uint64_t> cnt(NBK + 1, 0);
        for (std::size_t b = 0; b < NBK; ++b) {
            cnt[b + 1] = cnt[b];
            for (unsigned t = 0; t < T; ++t) cnt[b + 1] += tc[t][b];
        }
        std::vector<std::vector<uint64_t>>().swap(tc);
        out->order_.resize(U);
        std::vector<uint64_t> fill(cnt.begin(), cnt.end() - 1);
        for (uint64_t u = 0; u < U; ++u) out->order_[fill[top[u]]++] = static_cast<uint32_t>(u);
        std::vector<uint16_t>().swap(top);
        std::atomic<std::size_t> nb{0};
        par([&](unsigned) {
            std::vector<std::pair<uint64_t, uint32_t>> v;
            for (std::size_t b; (b = nb.fetch_add(1)) < NBK;) {
                const auto lo = out->order_.begin() + static_cast<std::ptrdiff_t>(cnt[b]);
                const auto hi = out->order_.begin() + static_cast<std::ptrdiff_t>(cnt[b + 1]);
                v.clear();
                for (auto it = lo; it != hi; ++it) v.push_back({head_of(*it), *it});
                std::sort(v.begin(), v.end());
                for (std::size_t q = 0; q < v.size(); ++q) lo[static_cast<std::ptrdiff_t>(q)] = v[q].second;
            }
        });
    }
    std::vector<uint64_t>().swap(cycles);
    // The fields stay in assembly order, with order_ mapping output positions to them, and off_ takes the unitigs' base offsets.
    uoff.pop_back();
    out->off_ = std::move(uoff);
    out->cov_ = std::move(cov);
    out->len_ = std::move(ulen);
    out->exit_masks_ = std::move(unitig_exit_masks);
    pt.lap("part: orient + sort");
    return out;
}

}
