#pragma once

/**
 * @file
 * @brief Partitioned compaction, which builds maximal unitigs without a global k-mer table.
 *
 * BCALM (Chikhi et al. 2014) compacts inside minimizer buckets and then glues the pieces across buckets.
 * We follow the same idea in our own form.
 * Pieces are walked in each counted partition's cache-resident table and glued by matching the canonical (k+1)-mer keys of their open ends.
 *
 * The counter assigns every canonical k-mer to exactly one minimizer partition, and consecutive k-mers sharing a minimizer share a partition.
 * Most edges of a unitig therefore never leave their partition.
 * Once a partition is counted, its L2-sized table holds every one of its k-mers with abundance and edges.
 * PartBuilder::partition() walks the partition's unitig pieces, called fragments, on exact masks with one cache-resident probe per step.
 * It applies the walker's compaction rule from dbg_unitigs.hpp: E -> X compacts iff E has one successor X, X has one predecessor, and the step is no self-loop, hairpin or palindrome.
 *
 * A fragment end whose single neighbour is in another partition is called open.
 * Its own half of the compaction test passed (out-degree 1 on the right, or in-degree 1 on the left), while the other half belongs to the neighbour.
 * Both sides key an open end by the canonical (k+1)-mer of its edge.
 * Two fragments are therefore glued exactly when both sides of the edge are open, which completes the compaction test.
 * End keys are bucketed, spilled past the budget, and matched.
 * Chains of glued fragments are then followed from their free ends (which pure cycles lack), and every unitig is assembled 2-bit packed from the fragment stream.
 *
 * The output equals the walker's byte for byte.
 * Each unitig starts at the smaller of its two oriented end k-mers, whose canonical form is the head.
 * Pure cycles start at their smallest canonical k-mer in canonical orientation, and all unitigs are sorted by head.
 * KC is the sum of abundances, and the GFA link masks are the end k-mers' exact masks.
 *
 * With -c > 1, an edge between partitions is real only if both k-mers are solid, which neither partition can check alone.
 * So in the counter's first pass every solid k-mer claims its cross-partition edges under the same keys.
 * An edge claimed only once has a non-solid end, and its bit is cleared when the second pass reaches the claimer's partition.
 *
 * The fragment stream and end records spill just like the counter's buffers.
 * Memory holds one partition table per thread, then two 32-bit words per fragment while chains are followed.
 * Those words first hold the partner ends and later the unitig id and offset.
 * Assembly reads them in fragment order, so past kLinkSpill fragments they are streamed back from disk instead of sitting next to the packed unitigs.
 * Finally the unitigs are sorted through an index, reading heads from the bases, so nothing is copied.
 * k must be at most 31, because the edge key is a (k+1)-mer that has to fit 64 bits.
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

#include "bits.hpp"
#include "hugepage.hpp"
#include "kmer.hpp"
#include "phase_timer.hpp"
#include "spill.hpp"
#include "wyhash.hpp"
#include "thread_budget.hpp"

namespace tpoptoa {

/// A solid k-mer of one partition, with its canonical key, abundance and in << 4 | out edge mask (exact inside the partition, while PartBuilder checks the cross-partition bits).
struct PartNode { uint64_t key; uint32_t cnt; uint32_t edge_mask; };

/// The unitigs of a partitioned compaction, sorted by head k-mer.
class PartUnitigs {
public:
    struct View { const uint8_t* seq; uint32_t len; uint64_t cov; };

    unsigned    kmer_len() const noexcept { return kmer_len_; }
    std::size_t core_size() const noexcept { return nkmers_; }
    std::size_t n_unitigs() const noexcept { return len_.size(); }
    bool        consumed() const noexcept { return false; }
    /// Unitig @p unitig_idx's exit masks, with the successors of its last k-mer in the low nibble and those of the reverse complement of its first in the high nibble.
    uint8_t     exit_mask(std::size_t unitig_idx) const noexcept { return exit_masks_[order_[unitig_idx]]; }

    /// Calls sink(count, unitig_at), where unitig_at(i) is the i-th unitig in output order, a view into a buffer of the calling thread that stays valid until that thread's next unitig_at().
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned = 0) const {
        sink(len_.size(), [&](std::size_t out_pos) {
            const std::size_t unitig = order_[out_pos];
            // Output order is not storage order, so the fields of a unitig a few places ahead are fetched now, and its bases a few places later, once its offset is in cache.
            constexpr std::size_t kFieldsAhead = 16, kBasesAhead = 8;

            if (out_pos + kFieldsAhead < order_.size()) {
                const uint32_t ahead = order_[out_pos + kFieldsAhead];

                __builtin_prefetch(&off_[ahead]); __builtin_prefetch(&len_[ahead]); __builtin_prefetch(&cov_[ahead]);
            }

            if (out_pos + kBasesAhead < order_.size()) __builtin_prefetch(bases_.data() + packed_byte(off_[order_[out_pos + kBasesAhead]]));

            const uint32_t len = len_[unitig];
            const uint64_t off = off_[unitig];
            thread_local std::vector<uint8_t> buf;

            buf.resize(len);
            // Four bases per byte, with unitigs starting byte-aligned, and lut.c[byte] spells the byte's four bases.
            // The table is built at compile time, so no unitig pays a first-call initialization check.
            struct Lut {
                char spell[kByteValues][kBasesPerByte] = {};

                constexpr Lut() {

                    for (int b = 0; b < kByteValues; ++b)
                        for (int q = 0; q < int(kBasesPerByte); ++q) spell[b][q] = kBaseChars[(b >> (kBitsPerBase * q)) & kBaseMask];
                }
            };
            static constexpr Lut lut{};
            const uint8_t* src = bases_.data() + packed_byte(off);
            uint32_t base_idx = 0;

            for (; base_idx + kBasesPerByte <= len; base_idx += kBasesPerByte) std::memcpy(buf.data() + base_idx, lut.spell[src[packed_byte(base_idx)]], kBasesPerByte);

            for (; base_idx < len; ++base_idx) buf[base_idx] = static_cast<uint8_t>(lut.spell[src[packed_byte(base_idx)]][base_idx % kBasesPerByte]);

            return View{buf.data(), len, cov_[unitig]};
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
    void unitigs(CB&& callback, unsigned nthreads = 0) const {
        unitigs_sorted([&](std::size_t count, const auto& unitig_at) {
            for (std::size_t i = 0; i < count; ++i) { const View view = unitig_at(i); callback(view.seq, view.len, view.cov); }
        }, nthreads);
    }

private:
    friend class PartBuilder;
    unsigned kmer_len_ = 0;
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

/// Whether a build at @p kmer_len uses partitioned compaction.
///
/// The order of precedence is OPTIK_DBG_MODE=part|table, then the CLI's --graph part|table read through count_limits().graph_mode, then auto from kPartModeMinK up.
/// Never for k > 31.
///
/// The threshold comes from hg38 chr1 at 8 threads, measured at commit 9b6e0c6 with identical GFA either way.
///   k = 31  5.5 s / 1.00 GB  vs the table's 6.8 s / 3.84 GB
///   k = 21  7.3 s / 1.22 GB  vs 7.1 s / 3.74 GB
///   k = 15  21 s / 5.2 GB    vs 19 s / 5.3 GB (60M short unitigs)
inline bool part_mode_wanted(unsigned kmer_len) {
    if (kmer_len > kDbgMaxK) return false;

    if (const char* env = std::getenv("OPTIK_DBG_MODE")) return std::strcmp(env, "part") == 0;

    const int mode = count_limits().graph_mode;

    if (mode != 0) return mode == 1;

    return kmer_len >= kPartModeMinK;
}

/// Builds PartUnitigs from the counter's partitions.
///
/// With @p filter, claim() runs during the first counting pass and match_claims() after it.
/// partition() runs once per counted partition in the second pass, then finish().
class PartBuilder {
public:
    static constexpr uint32_t kNone = ~uint32_t(0);

    /// @p kmer_len <= 31.
    /// @p filter is set when min_count > 1, in which case cross-partition edges are claimed and checked.
    PartBuilder(unsigned kmer_len, unsigned nthreads, std::size_t nparts, bool filter,
                uint64_t per_thread_bytes, std::string dir)
        : kmer_len_(kmer_len), nthreads_(nthreads), nparts_(nparts), filter_(filter), budget_(per_thread_bytes), dir_(std::move(dir)),
          kmask_(kmer_mask(kmer_len)), frag_writers_(nthreads), scratch_(nthreads) {
        if (kmer_len > kDbgMaxK) throw std::invalid_argument("partitioned compaction needs k <= 31");

        ends_ = std::make_unique<spill::BucketSpill<EndRec>>(nthreads_, kEndBuckets, budget_, dir_);

        if (filter_) claims_ = std::make_unique<spill::BucketSpill<ClaimRec>>(nthreads_, kEndBuckets, budget_, dir_);
    }

    /// The canonical (k+1)-mer of the edge of canonical k-mer @p canon_kmer through mask bit @p bit (bits 0-3 are out bases and 4-7 in bases, in the canonical frame).
    uint64_t mask_edge(uint64_t canon_kmer, unsigned bit) const noexcept {
        const uint64_t edge_kmer = bit < kOutEdgeBits ? ((canon_kmer << kBitsPerBase) | bit) : ((uint64_t(bit - kOutEdgeBits) << (kBitsPerBase * kmer_len_)) | canon_kmer);

        return canon1(edge_kmer);
    }

    /// First counting pass, with filter only.
    /// A solid k-mer of partition @p part has the cross-partition edge @p edge_key.
    void claim(unsigned tid, uint32_t part, uint64_t edge_key) {
        claims_->add(tid, bucket(edge_key), ClaimRec{edge_key, part});
    }
    /// After the first pass, edges claimed only once have a non-solid end and go to their claimer's per-partition invalid list.
    void match_claims() {
        invalid_edges_ = std::make_unique<spill::BucketSpill<uint64_t>>(nthreads_, nparts_, budget_, dir_);
        std::atomic<std::size_t> next_bucket{0};
        std::exception_ptr err;
        std::mutex err_mutex;

        tpoptoa::run_on_threads(nthreads_, [&](unsigned thread_idx) {
            try {
                std::vector<ClaimRec> claims;

                for (std::size_t bucket; (bucket = next_bucket.fetch_add(1)) < kEndBuckets;) {
                    claims_->load(bucket, claims);
                    std::sort(claims.begin(), claims.end(), [](const ClaimRec& lhs, const ClaimRec& rhs) {
                        const uint64_t lhs_key = lhs.key, rhs_key = rhs.key; return lhs_key < rhs_key; });
                    for (std::size_t run_start = 0; run_start < claims.size();) {
                        std::size_t run_end = run_start + 1;

                        while (run_end < claims.size() && claims[run_end].key == claims[run_start].key) ++run_end;

                        // A solid edge is claimed from both of its k-mers. A lone claim therefore means the far k-mer fell below the count cutoff.
                        if (run_end - run_start == 1) invalid_edges_->add(thread_idx, claims[run_start].part, claims[run_start].key);

                        run_start = run_end;
                    }

                    std::vector<ClaimRec>().swap(claims);
                }
            } catch (...) { std::lock_guard<std::mutex> guard(err_mutex); if (!err) err = std::current_exception(); }
        });
        if (err) std::rethrow_exception(err);

        claims_.reset();
    }

    /// Second pass, walking the fragments of counted partition @p part, whose solid k-mers are @p nodes.
    void partition(unsigned tid, uint32_t part, std::vector<PartNode>& nodes);

    /// Glue, assemble, orient and sort the fragments into the unitigs.
    std::unique_ptr<PartUnitigs> finish(PhaseTimer& timer);

private:
    static constexpr std::size_t kEndBuckets = 4096;              // end / claim buckets
    static constexpr uint64_t kBucketSeed = 0x51ED270B0A3C4E1DULL;
    static constexpr uint64_t kFibonacciMultiplier = kGoldenRatio64;
    /// Fragment f has ends 2f (left) and 2f + 1 (right).
    static constexpr uint64_t kEndsPerFragment = 2;
    /// After the walks, a fragment's two end slots are reused for its unitig id and its base offset.
    static constexpr uint64_t kLinkWordsPerFragment = kEndsPerFragment;
    /// An open end pairs with exactly one other, so an end key shared by more ends means a broken partition.
    static constexpr std::size_t kEndsPerEdge = 2;
    /// The local node index keeps its load at or below one half.
    static constexpr std::size_t kIndexLoadDen = 2;
    /// Fragment record header, a u32 k-mer count and a u64 coverage.
    static constexpr std::size_t kFragHeaderBytes = 12;
    /// Byte offset of the coverage in a fragment header, after the u32 k-mer count.
    static constexpr std::size_t kFragCovOffset = sizeof(uint32_t);
    /// Fragment and unitig ids are 31- and 32-bit, so these counts are out of range.
    static constexpr uint64_t kMaxFragments = uint64_t(1) << 31;
    static constexpr uint64_t kMaxUnitigs = uint64_t(1) << 32;
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
    /// One thread's fragments, records of [u32 nk][u64 cov][packed bases] kept in memory up to the budget and then appended to its spill file in chunks.
    struct FragWriter {
        std::vector<uint8_t> mem;
        spill::File file;
        std::vector<std::pair<uint64_t, uint64_t>> chunks;
        std::vector<uint16_t> nkmers;              // k-mers per fragment, where kNkInBig points to `big`
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
    uint64_t revcomp(uint64_t kmer) const noexcept { return reverse_complement(kmer, kmer_len_); }
    uint64_t canon(uint64_t kmer) const noexcept { const uint64_t rev = revcomp(kmer); return kmer < rev ? kmer : rev; }

    /// Canonical form of the (k+1)-mer @p kmer1.
    uint64_t canon1(uint64_t kmer1) const noexcept {
        const uint64_t rev = reverse_complement(kmer1, kmer_len_ + 1); return kmer1 < rev ? kmer1 : rev;
    }

    void emit(unsigned tid, const std::vector<uint8_t>& seq, uint64_t cov, int open_l, uint64_t key_l,
              int open_r, uint64_t key_r, uint8_t masks);

    unsigned kmer_len_, nthreads_;
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

/// Store one fragment, with 2-bit codes @p seq, coverage @p cov and exit @p masks, together with its open ends' keys, spilling the thread's buffer past the budget.
inline void PartBuilder::emit(unsigned tid, const std::vector<uint8_t>& seq, uint64_t cov, int open_l,
                              uint64_t key_l, int open_r, uint64_t key_r, uint8_t masks) {
    FragWriter& writer = frag_writers_[tid];
    const uint32_t local = static_cast<uint32_t>(writer.nkmers.size());
    const uint32_t nkmers = static_cast<uint32_t>(seq.size() - kmer_len_ + 1);

    // Fragments are runs inside one partition and almost always short.
    writer.nkmers.push_back(static_cast<uint16_t>(std::min<uint32_t>(nkmers, kNkInBig)));

    // A fragment longer than 65535 k-mers overflows the 16-bit length field and records its true length in a side list.
    if (nkmers >= kNkInBig) writer.big.push_back({local, nkmers});

    writer.masks.push_back(masks);
    const std::size_t start = writer.mem.size();
    const std::size_t nbytes = packed_bytes_for(seq.size());

    writer.mem.resize(start + kFragHeaderBytes + nbytes, 0);
    std::memcpy(writer.mem.data() + start, &nkmers, sizeof(uint32_t));
    std::memcpy(writer.mem.data() + start + kFragCovOffset, &cov, sizeof(uint64_t));
    uint8_t* dst = writer.mem.data() + start + kFragHeaderBytes;

    for (std::size_t i = 0; i < seq.size(); ++i) dst[packed_byte(i)] |= static_cast<uint8_t>(seq[i] << packed_shift(i));

    // Filing open ends by key lets the glue step pair them across partitions.
    if (open_l) ends_->add(tid, bucket(key_l), EndRec{key_l, local, static_cast<uint16_t>(tid), 0});

    if (open_r) ends_->add(tid, bucket(key_r), EndRec{key_r, local, static_cast<uint16_t>(tid), 1});

    // Over budget, fragments move to this thread's spill file.
    if (writer.mem.size() > budget_) {
        writer.file.open(dir_);
        std::vector<struct iovec> iov{{writer.mem.data(), writer.mem.size()}};

        writer.chunks.push_back({writer.file.size(), writer.mem.size()});
        writer.file.append(iov);
        writer.mem.clear();
    }
}

inline void PartBuilder::partition(unsigned tid, uint32_t part, std::vector<PartNode>& nodes) {
    const std::size_t count = nodes.size();

    if (count == 0) return;

    Scratch& scr = scratch_[tid];
    const unsigned kmer_len = kmer_len_;
    const uint64_t kmask = kmask_;
    const bool even_k = (kmer_len & 1u) == 0;
    // Local index of the partition's solid nodes, which stays cache-resident.
    std::size_t cap = kMinIndexSlots;

    while (cap < kIndexLoadDen * count) cap <<= 1;

    unsigned log2_cap = 0;

    while ((std::size_t(1) << log2_cap) < cap) ++log2_cap;

    scr.idx.assign(cap, 0);
    // Multiplying by the golden-ratio constant spreads clustered keys, and the top log2_cap bits of the product pick the slot.
    auto slot = [&](uint64_t key) {
        return static_cast<std::size_t>((key * kFibonacciMultiplier) >> (kWordBits - log2_cap));
    };

    for (std::size_t i = 0; i < count; ++i) {
        std::size_t slot_idx = slot(nodes[i].key);

        while (scr.idx[slot_idx]) slot_idx = (slot_idx + 1) & (cap - 1);

        // Storing index + 1 leaves 0 for an empty slot.
        scr.idx[slot_idx] = static_cast<uint32_t>(i + 1);
    }

    auto find = [&](uint64_t key) -> std::size_t {
        for (std::size_t slot_idx = slot(key);; slot_idx = (slot_idx + 1) & (cap - 1)) {
            const uint32_t entry = scr.idx[slot_idx];

            if (!entry) return SIZE_MAX;

            if (nodes[entry - 1].key == key) return entry - 1;
        }
    };
    const unsigned first_base_shift = kBitsPerBase * (kmer_len - 1);
    auto nbr = [&](uint64_t key, unsigned bit) -> uint64_t {          // neighbour's canonical k-mer
        return canon(bit < kOutEdgeBits ? (((key << kBitsPerBase) | bit) & kmask) : ((uint64_t(bit - kOutEdgeBits) << first_base_shift) | (key >> kBitsPerBase)));
    };

    // With -c > 1, cross-partition edges claimed once have a non-solid end.
    if (filter_) {
        invalid_edges_->load(part, scr.inv);

        if (!scr.inv.empty()) {
            std::sort(scr.inv.begin(), scr.inv.end());

            for (auto& node : nodes)
                for (unsigned b = 0; b < kEdgeMaskBits; ++b)
                    if ((node.edge_mask & (1u << b)) && find(nbr(node.key, b)) == SIZE_MAX
                        && std::binary_search(scr.inv.begin(), scr.inv.end(), mask_edge(node.key, b)))
                        node.edge_mask &= ~(1u << b);
        }
    }

    // The walker's compaction test (dbg_unitigs.hpp), applied to node i.
    auto cand = [&](uint64_t oriented, std::size_t node, uint64_t& succ, unsigned& succ_base) -> bool {
        const uint64_t canon_kmer = nodes[node].key;

        if (even_k && canon_kmer == revcomp(canon_kmer)) return false;

        const bool fwd = oriented == canon_kmer;
        const unsigned edge_bits = nodes[node].edge_mask;
        const unsigned dir_mask = fwd ? (edge_bits & kOutEdgeMask) : ((edge_bits >> kOutEdgeBits) & kOutEdgeMask);

        if (__builtin_popcount(dir_mask) != 1) return false;

        succ_base = fwd ? unsigned(__builtin_ctz(dir_mask)) : (unsigned(__builtin_ctz(dir_mask)) ^ kComplementMask);
        succ = ((oriented << kBitsPerBase) | succ_base) & kmask;

        if (canon(succ) == canon_kmer) return false;

        if (even_k && succ == revcomp(succ)) return false;

        return true;
    };
    auto accept = [&](uint64_t succ, std::size_t succ_node) -> bool {
        const unsigned edge_bits = nodes[succ_node].edge_mask;

        return __builtin_popcount(succ == canon(succ) ? ((edge_bits >> kOutEdgeBits) & kOutEdgeMask) : (edge_bits & kOutEdgeMask)) == 1;
    };
    // Successors of oriented k-mer E at node i, i.e. its exit mask.
    auto succ = [&](uint64_t oriented, std::size_t node) -> unsigned {
        const uint64_t canon_kmer = nodes[node].key;
        const unsigned edge_bits = nodes[node].edge_mask;
        unsigned exits = 0;

        if (oriented == canon_kmer) exits |= edge_bits & kOutEdgeMask;

        if (revcomp(oriented) == canon_kmer) {
            const unsigned in_edges = (edge_bits >> kOutEdgeBits) & kOutEdgeMask;

            for (unsigned b = 0; b < kNumBases; ++b) if (in_edges & (1u << (b ^ kComplementMask))) exits |= 1u << b;
        }

        return exits;
    };
    // Internal step E -> X, with both k-mers in this partition and the step compactable.
    auto step_in = [&](uint64_t oriented, std::size_t node, uint64_t& succ, std::size_t& succ_node) -> bool {
        unsigned succ_base;

        if (!cand(oriented, node, succ, succ_base)) return false;

        succ_node = find(canon(succ));
        return succ_node != SIZE_MAX && accept(succ, succ_node);
    };

    scr.vis.assign(count, 0);
    auto walk = [&](uint64_t start_kmer, std::size_t node, bool cycle) {
        // Fragment from oriented start E at node i, extended right as far as it compacts inside the partition.
        scr.seq.clear();

        for (int b = static_cast<int>(kmer_len) - 1; b >= 0; --b) scr.seq.push_back(static_cast<uint8_t>((start_kmer >> (kBitsPerBase * b)) & kBaseMask));

        uint64_t cov = nodes[node].cnt;

        scr.vis[node] = 1;
        const uint64_t start_c = canon(start_kmer);
        uint64_t cur = start_kmer; std::size_t cur_node = node;
        int open_r = 0; uint64_t key_r = 0;

        for (;;) {
            uint64_t succ; unsigned succ_base;

            if (!cand(cur, cur_node, succ, succ_base)) break;

            const uint64_t canon_succ = canon(succ);

            if (canon_succ == start_c) break;                         // closed a (local) cycle

            const std::size_t succ_node = find(canon_succ);

            if (succ_node == SIZE_MAX) {                             // the neighbour is elsewhere
                if (!cycle) { open_r = 1; key_r = canon1((cur << kBitsPerBase) | succ_base); }

                break;
            }

            if (!accept(succ, succ_node)) break;

            scr.seq.push_back(static_cast<uint8_t>(succ & kBaseMask));
            cov += nodes[succ_node].cnt;
            scr.vis[succ_node] = 1;
            cur = succ; cur_node = succ_node;
        }

        // The left end steps as rc(E).
        int open_l = 0; uint64_t key_l = 0;

        if (!cycle) {
            uint64_t pred; unsigned pred_base;

            if (cand(revcomp(start_kmer), node, pred, pred_base) && find(canon(pred)) == SIZE_MAX) { open_l = 1; key_l = canon1((revcomp(start_kmer) << kBitsPerBase) | pred_base); }
        }

        const uint8_t masks = static_cast<uint8_t>(succ(cur, cur_node) | (succ(revcomp(start_kmer), node) << kOutEdgeBits));

        emit(tid, scr.seq, cov, open_l, key_l, open_r, key_r, masks);
    };

    // Fragment starts are oriented k-mers with no internal compactable predecessor.
    for (std::size_t i = 0; i < count; ++i) {
        const uint64_t canon_kmer = nodes[i].key;
        const uint64_t rcc = revcomp(canon_kmer);

        for (unsigned o = 0; o < (canon_kmer == rcc ? 1u : kStrands); ++o) {
            if (scr.vis[i]) break;

            const uint64_t oriented = o == 0 ? canon_kmer : rcc;
            uint64_t pred; std::size_t pred_node;

            if (step_in(revcomp(oriented), i, pred, pred_node)) continue;           // compacts to the left here

            walk(oriented, i, false);
        }
    }

    // The rest lie on cycles inside the partition, the graph's pure cycles, and are emitted from their smallest canonical k-mer in canonical orientation, as the walker does.
    for (std::size_t i = 0; i < count; ++i) {
        if (scr.vis[i]) continue;

        uint64_t best = nodes[i].key; std::size_t best_i = i;
        uint64_t cur = nodes[i].key; std::size_t cur_node = i;

        for (std::size_t steps = 0; steps < count; ++steps) {
            uint64_t succ; std::size_t succ_node;

            if (!step_in(cur, cur_node, succ, succ_node) || succ_node == i) break;

            if (canon(succ) < best) { best = canon(succ); best_i = succ_node; }

            cur = succ; cur_node = succ_node;
        }

        walk(best, best_i, true);                             // marks the whole cycle
    }
}

inline std::unique_ptr<PartUnitigs> PartBuilder::finish(PhaseTimer& timer) {
    const unsigned kmer_len = kmer_len_;
    const unsigned num_threads = nthreads_;
    // Fragment ids of thread t follow those of threads < t.
    std::vector<uint64_t> base(num_threads + 1, 0);

    for (unsigned t = 0; t < num_threads; ++t) base[t + 1] = base[t] + frag_writers_[t].nkmers.size();

    const uint64_t nfrags = base[num_threads];

    if (nfrags >= kMaxFragments) throw std::length_error("partitioned compaction: 2^31 fragments or more");

    HugeVec<uint16_t> frag_kmers16(nfrags);
    std::vector<std::pair<uint64_t, uint32_t>> big;            // (fragment, k-mers), sorted
    HugeVec<uint8_t> frag_masks(nfrags);

    for (unsigned t = 0; t < num_threads; ++t) {
        std::copy(frag_writers_[t].nkmers.begin(), frag_writers_[t].nkmers.end(), frag_kmers16.begin() + static_cast<std::ptrdiff_t>(base[t]));

        for (const auto& big_entry : frag_writers_[t].big) big.push_back({base[t] + big_entry.first, big_entry.second});

        std::vector<std::pair<uint32_t, uint32_t>>().swap(frag_writers_[t].big);
        std::copy(frag_writers_[t].masks.begin(), frag_writers_[t].masks.end(), frag_masks.begin() + static_cast<std::ptrdiff_t>(base[t]));
        std::vector<uint16_t>().swap(frag_writers_[t].nkmers);
        std::vector<uint8_t>().swap(frag_writers_[t].masks);
    }

    auto frag_nkmers = [&](uint64_t frag) -> uint64_t {
        const uint16_t stored = frag_kmers16[frag];

        if (stored != kNkInBig) return stored;

        return std::lower_bound(big.begin(), big.end(), std::pair<uint64_t, uint32_t>{frag, 0})->second;
    };
    // Partner end of each fragment end (2f left, 2f+1 right).
    HugeVec<uint32_t> link(kEndsPerFragment * nfrags, kNone);
    auto par = [&](auto&& func) {
        std::exception_ptr err;
        std::mutex err_mutex;

        tpoptoa::run_on_threads(num_threads, [&](unsigned thread_idx) {
            try { func(thread_idx); } catch (...) { std::lock_guard<std::mutex> guard(err_mutex); if (!err) err = std::current_exception(); }
        });
        if (err) std::rethrow_exception(err);
    };

    {
        std::atomic<std::size_t> next_bucket{0};

        par([&](unsigned) {
            std::vector<EndRec> ends;

            for (std::size_t bucket; (bucket = next_bucket.fetch_add(1)) < kEndBuckets;) {
                // Read ahead the bucket this thread takes next.
                if (bucket + num_threads < kEndBuckets) ends_->will_need(bucket + num_threads);

                ends_->load(bucket, ends);
                std::sort(ends.begin(), ends.end(), [](const EndRec& lhs, const EndRec& rhs) {
                    const uint64_t lhs_key = lhs.key, rhs_key = rhs.key; return lhs_key < rhs_key; });
                for (std::size_t run_start = 0; run_start < ends.size();) {
                    std::size_t run_end = run_start + 1;

                    while (run_end < ends.size() && ends[run_end].key == ends[run_start].key) ++run_end;

                    if (run_end - run_start > kEndsPerEdge) throw std::logic_error("partitioned compaction: an edge with three open ends");

                    // Two ends sharing a key face each other across a partition boundary.
                    // An end key that occurs once has no partner across partitions and remains a free end of a chain.
                    if (run_end - run_start == kEndsPerEdge) {
                        const uint32_t end_a = static_cast<uint32_t>((base[ends[run_start].tid] + ends[run_start].local) * kEndsPerFragment + ends[run_start].side);
                        const uint32_t end_b = static_cast<uint32_t>((base[ends[run_start + 1].tid] + ends[run_start + 1].local) * kEndsPerFragment + ends[run_start + 1].side);

                        link[end_a] = end_b; link[end_b] = end_a;
                    }

                    run_start = run_end;
                }

                std::vector<EndRec>().swap(ends);
            }
        });
        ends_.reset();
    }
    malloc_trim(0);                  // the end-key buckets' freed pages

    if (timer.enabled) std::fprintf(stderr, "[timing] part: %llu fragments\n", static_cast<unsigned long long>(nfrags));

    timer.lap("part: glue (end keys matched)");

    // Chains are followed from free ends, where an end is entered and the fragment's other end exits, following partners until an end has none.
    // A chain is found from both of its free ends and kept from the smaller.
    // Its links are dependent loads, so a thread advances kInterleavedChains walks round-robin, and every walk prefetches its next link word one round ahead as the graph walker does, putting many misses in flight instead of one.
    std::vector<std::vector<uint32_t>> starts(num_threads);

    {
        const uint64_t num_ends = kEndsPerFragment * nfrags;

        par([&](unsigned thread_idx) {
            const uint64_t range_lo = num_ends * thread_idx / num_threads, range_hi = num_ends * (thread_idx + 1) / num_threads;
            uint64_t next = range_lo;

            struct Lane { uint32_t start, cur; uint64_t steps; bool live; } lanes[kInterleavedChains];
            auto start = [&](Lane& lane) -> bool {
                while (next < range_hi) {
                    const uint64_t end_idx = next++;

                    // Only a free end, one without a partner, can start a chain.
                    if (link[end_idx] != kNone) continue;

                    lane.start = lane.cur = static_cast<uint32_t>(end_idx); lane.steps = 0;
                    __builtin_prefetch(&link[end_idx ^ 1u]);
                    return lane.live = true;
                }

                return lane.live = false;
            };
            unsigned live = 0;

            for (auto& lane : lanes) live += start(lane);

            while (live) {
                for (auto& lane : lanes) {
                    if (!lane.live) continue;

                    const uint32_t exit = lane.cur ^ 1u;
                    const uint32_t next_end = link[exit];

                    if (next_end == kNone) {
                        // Both free ends reach this point for the same chain, and only the walk from the smaller one records it.
                        if (lane.start <= exit) starts[thread_idx].push_back(lane.start);

                        if (!start(lane)) --live;

                        continue;
                    }

                    // A walk from a free end visits each fragment once at most, so more steps than fragments means the links are corrupt.
                    if (++lane.steps > nfrags) throw std::logic_error("partitioned compaction: a chain from a free end does not end");

                    lane.cur = next_end;
                    __builtin_prefetch(&link[next_end ^ 1u]);
                }
            }
        });
    }
    // Unitig ids.
    // Each chain is walked again, and each of its fragments' link words is overwritten with the unitig id and the base offset ORed with the reversed bit.
    // A fragment belongs to one chain, so no other walk reads them.
    std::vector<uint64_t> sbase(num_threads + 1, 0);

    for (unsigned t = 0; t < num_threads; ++t) sbase[t + 1] = sbase[t] + starts[t].size();

    const uint64_t nchains = sbase[num_threads];
    std::vector<uint32_t> ulen(nchains);                 // bases (< 2^31, checked)
    std::vector<uint8_t> unitig_exit_masks(nchains);
    std::vector<uint8_t> done((nfrags + kBitsPerByte - 1) / kBitsPerByte, 0);

    par([&](unsigned thread_idx) {
        struct Lane { uint64_t uid, off; uint32_t cur; uint8_t mstart; bool live; } lanes[kInterleavedChains];
        std::size_t next = 0;
        auto start = [&](Lane& lane) -> bool {
            if (next >= starts[thread_idx].size()) return lane.live = false;

            lane.uid = sbase[thread_idx] + next;
            lane.cur = starts[thread_idx][next++];
            lane.off = 0;
            const uint32_t first_frag = lane.cur >> 1;

            // Exiting through the entered end at the start gives the unitig's left exit.
            lane.mstart = (lane.cur & 1u) ? (frag_masks[first_frag] & kOutEdgeMask) : (frag_masks[first_frag] >> kOutEdgeBits);
            __builtin_prefetch(&link[lane.cur & ~1u]);
            return lane.live = true;
        };
        unsigned live = 0;

        for (auto& lane : lanes) live += start(lane);

        while (live) {
            for (auto& lane : lanes) {
                if (!lane.live) continue;

                const uint32_t frag = lane.cur >> 1;
                const bool rev = lane.cur & 1u;                   // entered at its right end
                const uint32_t exit = lane.cur ^ 1u;
                const uint32_t next_end = link[exit];

                // The end slots now hold the fragment's unitig and offset.
                link[kLinkWordsPerFragment * frag] = static_cast<uint32_t>(lane.uid);
                link[kLinkWordsPerFragment * frag + 1] = static_cast<uint32_t>(lane.off) | (rev ? kReversedBit : 0u);
                __atomic_fetch_or(&done[frag / kBitsPerByte], static_cast<uint8_t>(1u << (frag % kBitsPerByte)), __ATOMIC_RELAXED);
                lane.off += frag_nkmers(frag);

                // A free exit closes the chain. Its exit edges come from the last fragment in walk orientation.
                if (next_end == kNone) {
                    const uint8_t mend = (exit & 1u) ? (frag_masks[frag] & kOutEdgeMask) : (frag_masks[frag] >> kOutEdgeBits);

                    unitig_exit_masks[lane.uid] = static_cast<uint8_t>(mend | (lane.mstart << kOutEdgeBits));

                    if (lane.off + kmer_len - 1 >= kReversedBit)
                        throw std::length_error("partitioned compaction: a unitig of 2^31 bases");

                    ulen[lane.uid] = static_cast<uint32_t>(lane.off + kmer_len - 1);

                    if (!start(lane)) --live;

                    continue;
                }

                lane.cur = next_end;
                __builtin_prefetch(&link[next_end & ~1u]);
                __builtin_prefetch(&frag_kmers16[next_end >> 1]);
            }
        }
    });
    // Pure cycles, where every end is glued and no walk reached them, are rare and run serially.
    std::vector<uint64_t> cycles;

    for (uint64_t f = 0; f < nfrags; ++f) {
        if (done[f / kBitsPerByte] & (1u << (f % kBitsPerByte))) continue;

        const uint64_t uid = ulen.size();
        uint64_t off = 0;
        uint32_t cur = static_cast<uint32_t>(f * kEndsPerFragment);             // enter at the left end

        for (uint64_t steps = 0; steps <= nfrags; ++steps) {
            const uint32_t frag = cur >> 1;
            const bool rev = cur & 1u;
            const uint32_t next_end = link[cur ^ 1u];

            if (done[frag / kBitsPerByte] & (1u << (frag % kBitsPerByte))) break;           // closed the cycle

            link[kLinkWordsPerFragment * frag] = static_cast<uint32_t>(uid);
            link[kLinkWordsPerFragment * frag + 1] = static_cast<uint32_t>(off) | (rev ? kReversedBit : 0u);
            done[frag / kBitsPerByte] |= static_cast<uint8_t>(1u << (frag % kBitsPerByte));
            off += frag_nkmers(frag);

            // Fragments the free-end walks never reached can only lie on closed cycles, where every exit has a partner.
            if (next_end == kNone) throw std::logic_error("partitioned compaction: an open chain left after the walks");

            cur = next_end;
        }

        if (off + kmer_len - 1 >= kReversedBit) throw std::length_error("partitioned compaction: a unitig of 2^31 bases");

        ulen.push_back(static_cast<uint32_t>(off + kmer_len - 1));
        unitig_exit_masks.push_back(0);
        cycles.push_back(uid);
    }

    const uint64_t num_unitigs = ulen.size();

    if (num_unitigs >= kMaxUnitigs) throw std::length_error("partitioned compaction: 2^32 unitigs or more");

    HugeVec<uint16_t>().swap(frag_kmers16);
    HugeVec<uint8_t>().swap(frag_masks);
    std::vector<uint8_t>().swap(done);

    for (auto& thread_starts : starts) std::vector<uint32_t>().swap(thread_starts);

    malloc_trim(0);
    timer.lap("part: chains");

    // The link words now hold each fragment's (unitig, offset), and the assembly reads them in fragment order, the order of the fragment stream.
    // Past kLinkSpill fragments they are written to the spill directory and read back in slices beside the stream, so the 8 bytes per fragment are not held in memory with the bases.
    constexpr uint64_t kSlice = uint64_t(1) << 18;               // fragments per read
    spill::File link_file;
    const bool links_spilled = nfrags >= kLinkSpill || std::getenv("OPTIK_SPILL_TEST");   // tests spill at any size

    if (links_spilled) {
        link_file.open(dir_);
        std::vector<struct iovec> iov{{link.data(), kLinkWordsPerFragment * nfrags * sizeof(uint32_t)}};

        link_file.append(iov);
        HugeVec<uint32_t>().swap(link);
        malloc_trim(0);
    }

    // Assembly writes unitig u's bases at uoff[u] (byte-aligned, so no two unitigs share a byte) and each fragment's bases at its offset, oriented.
    // Glued fragments overlap by k - 1 bases of equal value, which are written with atomic OR.
    auto out = std::make_unique<PartUnitigs>();

    out->kmer_len_ = kmer_len;
    std::vector<uint64_t> uoff(num_unitigs + 1, 0);

    for (uint64_t u = 0; u < num_unitigs; ++u) uoff[u + 1] = uoff[u] + ((ulen[u] + kBasesPerByte - 1) & ~uint64_t(kBasesPerByte - 1));

    out->bases_.assign((packed_byte(uoff[num_unitigs]) + sizeof(uint64_t) - 1) & ~uint64_t(sizeof(uint64_t) - 1), 0);       // whole 64-bit words
    uint8_t* arena = out->bases_.data();
    uint64_t* arena_words = reinterpret_cast<uint64_t*>(arena);
    std::vector<uint64_t> cov(num_unitigs, 0);

    par([&](unsigned thread_idx) {
        FragWriter& writer = frag_writers_[thread_idx];
        uint64_t frag_idx = 0;
        std::vector<uint8_t> buf;
        // Fragment f's two link words, taken from memory or from the current slice.
        std::vector<uint32_t> slice_words;
        uint64_t lw_lo = 0, lw_hi = 0;
        auto words = [&](uint64_t frag) -> const uint32_t* {
            if (!links_spilled) return &link[kLinkWordsPerFragment * frag];

            // Spilled link words come back a slice at a time beside the fragment stream.
            if (frag >= lw_hi) {
                lw_lo = frag;
                lw_hi = std::min(frag + kSlice, base[thread_idx + 1]);
                slice_words.resize(kLinkWordsPerFragment * (lw_hi - lw_lo));
                link_file.read(kLinkWordsPerFragment * lw_lo * sizeof(uint32_t), slice_words.size() * sizeof(uint32_t),
                           reinterpret_cast<uint8_t*>(slice_words.data()));
            }

            return &slice_words[kLinkWordsPerFragment * (frag - lw_lo)];
        };
        auto parse = [&](const uint8_t* data, std::size_t len) {
            std::size_t offset = 0;

            while (offset < len) {
                uint32_t nkmers; uint64_t frag_cov;

                std::memcpy(&nkmers, data + offset, sizeof(uint32_t));
                std::memcpy(&frag_cov, data + offset + kFragCovOffset, sizeof(uint64_t));
                const uint32_t frag_len = nkmers + kmer_len - 1;
                const uint8_t* packed = data + offset + kFragHeaderBytes;
                const uint64_t frag = base[thread_idx] + frag_idx;
                const uint32_t* link_words = words(frag);
                const uint64_t uid = link_words[0];
                const uint32_t offset_word = link_words[1];
                const bool rev = offset_word & kReversedBit;
                // The low bits of the offset word give the fragment's position in its unitig, and the top bit asks for it to be written reverse-complemented.
                const uint64_t dest = uoff[uid] + (offset_word & kOffsetMask);

                __atomic_fetch_add(&cov[uid], frag_cov, __ATOMIC_RELAXED);
                // 32 bases per 64-bit word take one atomic OR per word (glued fragments share their k - 1 overlapping bases, with equal values).
                constexpr uint64_t kBasesPerWord = kWordBits / kBitsPerBase;
                uint64_t acc = 0, pos = dest;

                for (uint32_t q = 0; q < frag_len; ++q, ++pos) {
                    const uint32_t src = rev ? frag_len - 1 - q : q;
                    uint64_t code = (packed[packed_byte(src)] >> packed_shift(src)) & kBaseMask;

                    if (rev) code ^= kComplementMask;

                    acc |= code << ((pos % kBasesPerWord) * kBitsPerBase);

                    if (pos % kBasesPerWord == kBasesPerWord - 1) {
                        if (acc) __atomic_fetch_or(&arena_words[pos / kBasesPerWord], acc, __ATOMIC_RELAXED);

                        acc = 0;
                    }
                }

                if (acc) __atomic_fetch_or(&arena_words[(pos - 1) / kBasesPerWord], acc, __ATOMIC_RELAXED);

                offset += kFragHeaderBytes + packed_bytes_for(frag_len);
                ++frag_idx;
            }
        };

        // Spilled chunks precede the in-memory tail, matching fragment numbering.
        for (const auto& chunk : writer.chunks) {
            buf.resize(chunk.second);
            writer.file.read(chunk.first, chunk.second, buf.data());
            parse(buf.data(), buf.size());
        }

        parse(writer.mem.data(), writer.mem.size());
        std::vector<uint8_t>().swap(writer.mem);
    });
    frag_writers_.clear();
    HugeVec<uint32_t>().swap(link);

    if (links_spilled) link_file.punch(0, link_file.size());

    timer.lap("part: assemble");

    // Orientation, with the smaller oriented end k-mer first, then heads and cycles.
    auto get = [&](uint64_t pos) -> uint8_t { return (arena[packed_byte(pos)] >> packed_shift(pos)) & kBaseMask; };
    auto put = [&](uint64_t pos, uint8_t code) {
        uint8_t& byte = arena[packed_byte(pos)];

        byte = static_cast<uint8_t>((byte & ~(kBaseMask << packed_shift(pos))) | (code << packed_shift(pos)));
    };
    auto kmer_at = [&](uint64_t pos) {
        uint64_t kmer = 0;

        for (unsigned q = 0; q < kmer_len; ++q) kmer = (kmer << kBitsPerBase) | get(pos + q);

        return kmer;
    };

    par([&](unsigned thread_idx) {
        // Cycle ids are appended in increasing order, so the list is sorted.
        auto next_cycle = std::lower_bound(cycles.begin(), cycles.end(), num_unitigs * thread_idx / num_threads);

        for (uint64_t u = num_unitigs * thread_idx / num_threads; u < num_unitigs * (thread_idx + 1) / num_threads; ++u) {
            const uint64_t unit_start = uoff[u], unit_len = ulen[u];

            if (next_cycle != cycles.end() && *next_cycle == u) {
                ++next_cycle;
                // A circular sequence of n = L - K + 1 bases.
                // The walker starts at the smallest canonical k-mer in canonical orientation.
                const uint64_t cycle_len = unit_len - kmer_len + 1;
                std::vector<uint8_t> codes(cycle_len);

                for (uint64_t i = 0; i < cycle_len; ++i) codes[i] = get(unit_start + i);

                uint64_t fwd_kmer = 0;

                for (unsigned q = 0; q < kmer_len; ++q) fwd_kmer = (fwd_kmer << kBitsPerBase) | codes[q % cycle_len];

                uint64_t best = ~uint64_t(0), best_idx = 0; bool bfwd = true;

                for (uint64_t i = 0; i < cycle_len; ++i) {
                    if (i) fwd_kmer = ((fwd_kmer << kBitsPerBase) | codes[(i + kmer_len - 1) % cycle_len]) & kmask_;

                    const uint64_t rev = revcomp(fwd_kmer);
                    const uint64_t canon_kmer = fwd_kmer < rev ? fwd_kmer : rev;

                    if (canon_kmer < best) { best = canon_kmer; best_idx = i; bfwd = fwd_kmer <= rev; }
                }

                // Forward reads c[bi], c[bi+1], ..., and reversed reads the complement walking back from the k-mer's last base.
                // The reverse index adds whole cycles first, so it stays positive before the modulo.
                constexpr uint64_t kWrapCycles = 2;

                for (uint64_t q = 0; q < unit_len; ++q) {
                    const uint8_t code = bfwd ? codes[(best_idx + q) % cycle_len]
                                              : static_cast<uint8_t>(codes[(best_idx + kmer_len - 1 + cycle_len * (unit_len / cycle_len + kWrapCycles) - q) % cycle_len] ^ kComplementMask);

                    put(unit_start + q, code);
                }

                const uint8_t first_succ = get(unit_start + kmer_len - 1);           // the cycle continues into its start
                const uint8_t pre = get(unit_start + unit_len - kmer_len);                   // the last k-mer's first base

                unitig_exit_masks[u] = static_cast<uint8_t>((1u << first_succ) | ((1u << (pre ^ kComplementMask)) << kOutEdgeBits));
                continue;
            }

            const uint64_t first_kmer = kmer_at(unit_start), last_kmer = kmer_at(unit_start + unit_len - kmer_len);

            if (first_kmer > revcomp(last_kmer)) {                                     // emitted from the other end
                uint64_t i = 0, j = unit_len - 1;

                for (; i < j; ++i, --j) {
                    const uint8_t left_base = get(unit_start + i), right_base = get(unit_start + j);

                    put(unit_start + i, right_base ^ kComplementMask); put(unit_start + j, left_base ^ kComplementMask);
                }

                // An odd length leaves the middle base, where the indices meet, still to complement.
                if (i == j) put(unit_start + i, get(unit_start + i) ^ kComplementMask);

                unitig_exit_masks[u] = static_cast<uint8_t>((unitig_exit_masks[u] >> kOutEdgeBits) | (unitig_exit_masks[u] << kOutEdgeBits));
            }
        }
    });
    // Sorted by head.
    // Oriented, a unitig's head is the canonical form of its first k-mer (a path starts at its smaller oriented end k-mer, a cycle at its smallest canonical k-mer), so a head is read from the bases where needed instead of being kept per unitig, through an index bucketed by the heads' top bits, whose buckets are sorted in parallel with their heads gathered just for them.
    auto head_of = [&](uint64_t unitig) { return canon(kmer_at(uoff[unitig])); };
    uint64_t nkm = 0;

    for (uint64_t u = 0; u < num_unitigs; ++u) nkm += ulen[u] - kmer_len + 1;

    out->nkmers_ = nkm;
    {
        const unsigned kRadixBits = 12, key_bits = kBitsPerBase * kmer_len, shift = key_bits > kRadixBits ? key_bits - kRadixBits : 0;
        const std::size_t NBK = std::size_t(1) << std::min(kRadixBits, key_bits);
        std::vector<uint16_t> top(num_unitigs);
        std::vector<std::vector<uint64_t>> thread_counts(num_threads, std::vector<uint64_t>(NBK, 0));

        // The first pass buckets every unitig by the top kRadixBits of its head, counting per thread so no counter is shared.
        par([&](unsigned thread_idx) {
            for (uint64_t u = num_unitigs * thread_idx / num_threads; u < num_unitigs * (thread_idx + 1) / num_threads; ++u) {
                top[u] = static_cast<uint16_t>(head_of(u) >> shift);
                ++thread_counts[thread_idx][top[u]];
            }
        });
        std::vector<uint64_t> cnt(NBK + 1, 0);

        for (std::size_t b = 0; b < NBK; ++b) {
            cnt[b + 1] = cnt[b];

            for (unsigned t = 0; t < num_threads; ++t) cnt[b + 1] += thread_counts[t][b];
        }

        std::vector<std::vector<uint64_t>>().swap(thread_counts);
        out->order_.resize(num_unitigs);
        std::vector<uint64_t> fill(cnt.begin(), cnt.end() - 1);

        // The scatter is serial, which keeps unitigs of one bucket in assembly order before the bucket sort.
        for (uint64_t u = 0; u < num_unitigs; ++u) out->order_[fill[top[u]]++] = static_cast<uint32_t>(u);

        std::vector<uint16_t>().swap(top);
        std::atomic<std::size_t> next_bucket{0};

        par([&](unsigned) {
            std::vector<std::pair<uint64_t, uint32_t>> heads;

            for (std::size_t bucket; (bucket = next_bucket.fetch_add(1)) < NBK;) {
                const auto range_begin = out->order_.begin() + static_cast<std::ptrdiff_t>(cnt[bucket]);
                const auto range_end = out->order_.begin() + static_cast<std::ptrdiff_t>(cnt[bucket + 1]);

                heads.clear();

                // Heads are rebuilt from the bases only for the bucket being sorted, so all of them are never held at once.
                for (auto entry_it = range_begin; entry_it != range_end; ++entry_it) heads.push_back({head_of(*entry_it), *entry_it});

                std::sort(heads.begin(), heads.end());

                for (std::size_t q = 0; q < heads.size(); ++q) range_begin[static_cast<std::ptrdiff_t>(q)] = heads[q].second;
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
    timer.lap("part: orient + sort");
    return out;
}

}
