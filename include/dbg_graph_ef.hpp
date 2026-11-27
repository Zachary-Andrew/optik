#pragma once

/**
 * @file
 * @brief Compacted de Bruijn graph over Elias-Fano coded keys.
 * It exists for k-mer sets whose open-addressing table (DbgGraph, about 16 bytes per k-mer) won't fit in memory.
 * A whole human genome, for instance, has about 2.5e9 distinct 31-mers, which would need 40 GB of table.
 *
 * @par Keys
 * The sorted canonical k-mers form an Elias-Fano sequence (Elias 1974, Fano 1971), with select sampling as in Vigna 2013.
 * Each key keeps l = floor(log2(4^k / n)) packed low bits, while the high parts form a unary bitvector of about 2n bits.
 * That comes to about 2 + log2(4^k / n) bits per k-mer, or 32 for the human genome at k = 31.
 * rank(c) finds c's high-part bucket with a sampled select0 and then compares the bucket's few low parts, so membership is exact.
 * key(i), in turn, is a sampled select1.
 * Because node i is the i-th smallest k-mer, iteration comes out sorted.
 *
 * @par Nodes
 * Each k-mer gets one 32-bit word laid out as in DbgGraph: 22 bits of abundance, the walker's two flags and the in << 4 | out edge mask.
 * Abundances that don't fit go to a sorted side list.
 *
 * Overall the graph takes about 8 bytes per k-mer against the table's 16.
 * The price is 2-3 dependent cache misses per lookup instead of one.
 * For that reason the counter only picks it when the table won't fit the memory limit (src/kmer_count.cpp).
 *
 * The node interface matches DbgGraph's, so the unitig walker (dbg_unitigs.hpp) and the GFA writer (src/ffi_dbg.cpp) run on either.
 * Their output is identical, which bench/spill_check.py verifies.
 */

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include "bits.hpp"
#include "dbg_graph.hpp"
#include "hugepage.hpp"
#include "thread_budget.hpp"

#if defined(__BMI2__)
#include <immintrin.h>
#endif

namespace tpoptoa {

class DbgGraphEF {
public:
    static constexpr std::size_t kNpos = DbgGraph::kNpos;
    static constexpr uint32_t kCountMax = DbgGraph::kCountMax;
    static constexpr unsigned kFlagShift = DbgGraph::kFlagShift, kMaskShift = DbgGraph::kMaskShift;
    static constexpr unsigned kWalkFlags = DbgGraph::kWalkFlags;

    explicit DbgGraphEF(unsigned kmer_len) : kmer_len_(kmer_len) {}

    unsigned    kmer_len() const noexcept { return kmer_len_; }
    std::size_t core_size() const noexcept { return num_nodes_; }
    std::size_t capacity() const noexcept { return cap_; }

    /// Covers keys, node words and select samples.
    static constexpr uint64_t bytes_for(std::size_t nodes, unsigned kmer_len) noexcept {
        if (nodes == 0) return 0;

        const unsigned low_width = low_bits(nodes, kmer_len);
        const uint64_t high_len = nodes + ((universe(kmer_len) - 1) >> low_width) + kHighBitsSlack;
        const uint64_t low_bytes = bits_to_bytes(nodes * uint64_t(low_width)), high_bytes = bits_to_bytes(high_len);
        const uint64_t sample_bytes = (high_bytes / kSelect0Sample + nodes / kSelect1Sample) * sizeof(uint64_t);

        return low_bytes + high_bytes + nodes * sizeof(uint32_t) + sample_bytes + kScalarBytes;
    }

    // Building takes three steps: begin(n), then put() for every node (any order, any thread), then finish_put().

    void begin(std::size_t nodes) {
        num_nodes_ = cap_ = 0;
        reserve_ = nodes;
        nlow_ = nodes ? low_bits(nodes, kmer_len_) : 0;
        const uint64_t hbits = nodes + ((universe(kmer_len_) - 1) >> nlow_) + kHighBitsSlack;

        low_parts_.assign(words_for_bits(nodes * nlow_) + kPadWords, 0);
        high_bits_.assign(hbits / kWordBits + kPadWords, 0);
        node_words_.assign(nodes, 0);
        select1_samples_.clear();
        saturated_.clear();
        last_ = 0;
    }
    /// Stores @p key, abundance @p cnt and observed edges @p raw for one node.
    /// Concurrent calls are safe for distinct @p node_idx, provided keys increase with @p node_idx.
    void put(std::size_t node_idx, uint64_t key, uint32_t cnt, uint8_t raw) {
        if (nlow_) {
            const uint64_t lowv = key & ((uint64_t(1) << nlow_) - 1), bit = node_idx * nlow_;

            __atomic_fetch_or(&low_parts_[word_index(bit)], lowv << bit_offset(bit), __ATOMIC_RELAXED);

            if (bit_offset(bit) + nlow_ > kWordBits) __atomic_fetch_or(&low_parts_[word_index(bit) + 1], lowv >> (kWordBits - bit_offset(bit)), __ATOMIC_RELAXED);
        }

        const uint64_t high_pos = (key >> nlow_) + node_idx;

        __atomic_fetch_or(&high_bits_[word_index(high_pos)], uint64_t(1) << bit_offset(high_pos), __ATOMIC_RELAXED);

        if (cnt >= kCountMax) { std::lock_guard<std::mutex> guard(saturated_mutex_); saturated_.push_back({node_idx, cnt}); }

        node_words_[node_idx] = std::min(cnt, kCountMax) | (uint32_t(raw) << kMaskShift);
    }
    /// Takes the select1 samples once every node is in, and checks the key order on the way.
    /// Keys that weren't distinct and increasing throw std::logic_error.
    void finish_put() {
        num_nodes_ = reserve_;
        select1_samples_.clear();
        uint64_t ones = 0;

        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {
            const uint64_t word = high_bits_[wi];
            const uint64_t ones_in_word = static_cast<uint64_t>(__builtin_popcountll(word));

            for (uint64_t sample_rank = (ones + kSelect1Sample - 1) / kSelect1Sample * kSelect1Sample; sample_rank < ones + ones_in_word; sample_rank += kSelect1Sample)
                select1_samples_.push_back(first_bit_of_word(uint64_t(wi)) + select_in_word(word, static_cast<unsigned>(sample_rank - ones)));

            ones += ones_in_word;
        }

        if (ones != num_nodes_) throw std::logic_error("DbgGraphEF: keys were not distinct and increasing");

        last_ = num_nodes_ ? key(num_nodes_ - 1) : 0;
        finish();
    }
    void finish() {
        if (num_nodes_ != reserve_) throw std::logic_error("DbgGraphEF: fewer keys appended than reserved");

        cap_ = num_nodes_;
        // Record the position of every kSelect0Sample-th zero in the high bits.
        select0_samples_.clear();
        uint64_t zeros = 0;

        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {
            const uint64_t zero_bits = ~high_bits_[wi];
            const uint64_t zeros_in_word = static_cast<uint64_t>(__builtin_popcountll(zero_bits));

            for (uint64_t sample_rank = (zeros + kSelect0Sample - 1) / kSelect0Sample * kSelect0Sample; sample_rank < zeros + zeros_in_word; sample_rank += kSelect0Sample)
                select0_samples_.push_back(first_bit_of_word(uint64_t(wi)) + select_in_word(zero_bits, static_cast<unsigned>(sample_rank - zeros)));

            zeros += zeros_in_word;
        }

        nzeros_ = zeros;
        std::sort(saturated_.begin(), saturated_.end());
    }

    /// Returns kNpos when the k-mer is absent.
    std::size_t rank(uint64_t canon_kmer) const noexcept {
        if (num_nodes_ == 0 || canon_kmer > last_) return kNpos;

        const uint64_t high = canon_kmer >> nlow_;
        uint64_t pos;

        if (high == 0) pos = 0;
        else {
            // In Elias-Fano, bucket h starts right after the h-th zero of the high bits, which select0 locates from a sampled position.
            if (high - 1 >= nzeros_) return kNpos;

            pos = select0(high - 1) + 1;
        }

        uint64_t idx = pos - high;
        const uint64_t low_part = canon_kmer & low_mask();

        // Each set bit in the high array is one key of the bucket, in low-bit order. The scan can stop at the first larger key.
        while (idx < num_nodes_ && ((high_bits_[word_index(pos)] >> bit_offset(pos)) & 1)) {
            const uint64_t stored_low = low(idx);

            if (stored_low == low_part) return idx;

            if (stored_low > low_part) return kNpos;

            ++idx; ++pos;
        }

        return kNpos;
    }
    uint64_t key(std::size_t node_idx) const noexcept { return ((select1(node_idx) - node_idx) << nlow_) | low(node_idx); }
    bool used(std::size_t node_idx) const noexcept { return node_idx < num_nodes_; }

    /// Prefetches only the select0 sample rank() starts from.
    /// The high and low words depend on that sample's value, and guessing their position from the k-mer alone misses by about sqrt(n) bits.
    void prefetch(uint64_t canon_kmer) const noexcept {
        const uint64_t high = canon_kmer >> nlow_;

        if (high && (high - 1) / kSelect0Sample < select0_samples_.size()) __builtin_prefetch(&select0_samples_[(high - 1) / kSelect0Sample]);
    }
    void prefetch_w(uint64_t canon_kmer) const noexcept { prefetch(canon_kmer); }

    unsigned edge_mask(std::size_t node_idx) const noexcept { return node(node_idx) >> kMaskShift; }

    uint64_t count(std::size_t node_idx) const noexcept {
        const uint32_t stored = node(node_idx) & kCountMax;

        if (stored < kCountMax) return stored;

        const auto found = std::lower_bound(saturated_.begin(), saturated_.end(), std::pair<std::size_t, uint32_t>{node_idx, 0});

        return found != saturated_.end() && found->first == node_idx ? found->second : stored;
    }
    unsigned flags_or(std::size_t node_idx, unsigned flag_bits) const noexcept {
        return (__atomic_fetch_or(&node_word(node_idx), uint32_t(flag_bits) << kFlagShift, __ATOMIC_ACQ_REL) >> kFlagShift) & kWalkFlags;
    }
    unsigned flags(std::size_t node_idx) const noexcept { return (node(node_idx) >> kFlagShift) & kWalkFlags; }
    void clear_flags(std::size_t node_idx) const noexcept { node_word(node_idx) &= ~(uint32_t(kWalkFlags) << kFlagShift); }

    /// Lets the node words be reused after the last walk, as DbgGraph::aux_clear() describes.
    void aux_clear(unsigned nthreads) {
        consumed_ = true;
        const unsigned threads = nthreads ? nthreads : 1;

        tpoptoa::run_on_threads(threads, [&](unsigned thread_idx) {
            std::fill(node_words_.begin() + num_nodes_ * thread_idx / threads, node_words_.begin() + num_nodes_ * (thread_idx + 1) / threads, 0u);
        });
    }
    void aux_set(std::size_t node_idx, uint32_t value) noexcept { __atomic_store_n(&node_word(node_idx), value, __ATOMIC_RELAXED); }
    uint32_t aux(std::size_t node_idx) const noexcept { return node(node_idx); }
    bool consumed() const noexcept { return consumed_; }

    unsigned out_mask(uint64_t canon_kmer) const noexcept {
        const std::size_t node_idx = rank(canon_kmer);

        return node_idx == kNpos ? 0u : edge_mask(node_idx) & kOutEdgeMask;
    }
    unsigned in_mask(uint64_t canon_kmer) const noexcept {
        const std::size_t node_idx = rank(canon_kmer);

        return node_idx == kNpos ? 0u : edge_mask(node_idx) >> kOutEdgeBits;
    }
    void masks_from_raw() noexcept {}

    /// Drops observed edges whose neighbour isn't in the graph.
    /// Each node touches only its own mask and its neighbours' keys, so nodes run in parallel safely.
    void resolve_masks() {
        const unsigned kmer_len = kmer_len_;
        const uint64_t kmask = kmer_mask(kmer_len);
        const unsigned first_base_shift = kBitsPerBase * (kmer_len - 1);
        auto canon = [&](uint64_t kmer) { const uint64_t rev_comp = reverse_complement(kmer, kmer_len); return kmer < rev_comp ? kmer : rev_comp; };
        const int64_t num_nodes = static_cast<int64_t>(num_nodes_);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, kResolveChunk)
#endif
        for (int64_t ii = 0; ii < num_nodes; ++ii) {
            const std::size_t node_idx = static_cast<std::size_t>(ii);
            const uint64_t kmer = key(node_idx);
            const unsigned raw = edge_mask(node_idx);
            unsigned mask = 0;

            for (unsigned b = 0; b < kNumBases; ++b) {
                // resolve_masks() clears an edge bit when the neighbour it points to fell below the count cutoff.
                if ((raw & (1u << b)) && rank(canon(((kmer << kBitsPerBase) | b) & kmask)) != kNpos) mask |= 1u << b;

                if ((raw & (1u << (b + kOutEdgeBits))) && rank(canon(((uint64_t(b) << first_base_shift) | (kmer >> kBitsPerBase)) & kmask)) != kNpos)
                    mask |= 1u << (b + kOutEdgeBits);
            }

            node_word(node_idx) = (node_word(node_idx) & ((uint32_t(1) << kMaskShift) - 1)) | (uint32_t(mask) << kMaskShift);
        }
    }

    template <typename CB>
    void unitigs(CB&& callback, unsigned nthreads = 0) const {
        DbgUnitigExtractor::run_masks(*this, std::forward<CB>(callback), nthreads);
    }
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned nthreads = 0) const {
        DbgUnitigExtractor::run_sorted(*this, std::forward<Sink>(sink), nthreads);
    }

    /// Node order already is key order, so no sort is needed.
    void export_sorted(std::vector<std::pair<uint64_t, uint32_t>>& pairs_out, std::vector<uint8_t>& masks) const {
        if (consumed_) throw std::logic_error("DbgGraphEF: node words were reused (aux_clear)");

        pairs_out.resize(num_nodes_);
        masks.resize(num_nodes_);

        for (std::size_t i = 0; i < num_nodes_; ++i) {
            pairs_out[i] = {key(i), static_cast<uint32_t>(count(i))};
            masks[i] = static_cast<uint8_t>(edge_mask(i));
        }
    }
    /// Writes the same OKKDBG archive as DbgGraph::save(), streamed straight from the sorted keys.
    bool save(const char* path) const {
        if (consumed_) return false;

        return atomic_write(path, [&](std::ofstream& out) {
            out.write(DbgGraph::kKdbgMagic, DbgGraph::kKdbgMagicLen);
            out.put(static_cast<char>(kmer_len_));
            const uint64_t node_count = num_nodes_;

            out.write(reinterpret_cast<const char*>(&node_count), sizeof(node_count));
            // The layout matches DbgGraph::save() so that either graph loads either file.
            std::vector<char> buf;
            auto flush = [&] { out.write(buf.data(), static_cast<std::streamsize>(buf.size())); buf.clear(); };

            for (std::size_t i = 0; i < num_nodes_; ++i) {
                const uint64_t kmer = key(i);

                buf.insert(buf.end(), reinterpret_cast<const char*>(&kmer), reinterpret_cast<const char*>(&kmer) + sizeof(kmer));

                if (buf.size() >= kSaveBufferBytes) flush();
            }

            for (std::size_t i = 0; i < num_nodes_; ++i) {
                const uint32_t abundance = static_cast<uint32_t>(count(i));

                buf.insert(buf.end(), reinterpret_cast<const char*>(&abundance), reinterpret_cast<const char*>(&abundance) + sizeof(abundance));

                if (buf.size() >= kSaveBufferBytes) flush();
            }

            for (std::size_t i = 0; i < num_nodes_; ++i) {
                buf.push_back(static_cast<char>(edge_mask(i)));

                if (buf.size() >= kSaveBufferBytes) flush();
            }

            flush();
        });
    }

private:
    friend struct DbgUnitigExtractor;
    static constexpr uint64_t kSelect0Sample = 64, kSelect1Sample = 256;   // one sample per this many zeros / ones
    static constexpr int kResolveChunk = 1 << 14;                          // nodes per OpenMP task in resolve_masks()
    static constexpr std::size_t kSaveBufferBytes = std::size_t(1) << 20;
    static constexpr uint64_t kScalarBytes = 64;      // rough size of the members, for bytes()
    /// low() reads two words at a time, so each bit array is padded by this many words.
    static constexpr std::size_t kPadWords = 2;
    /// The high bitvector gets a couple of spare bits past n + (universe >> l), so select never runs off its end.
    static constexpr uint64_t kHighBitsSlack = 2;

    /// 4^k distinct k-mers, saturating at k = kMaxK.
    static constexpr uint64_t universe(unsigned kmer_len) noexcept { return kmer_len >= kMaxK ? ~uint64_t(0) : uint64_t(1) << (kBitsPerBase * kmer_len); }

    /// floor(log2(universe / n)), capped one short of a full word.
    static constexpr unsigned low_bits(std::size_t nodes, unsigned kmer_len) noexcept {
        unsigned width = 0;
        const uint64_t ratio = universe(kmer_len) / nodes;

        while (width < kWordBits - 1 && (uint64_t(1) << (width + 1)) <= ratio) ++width;

        return width;
    }
    uint64_t low_mask() const noexcept { return nlow_ ? (uint64_t(1) << nlow_) - 1 : 0; }

    uint64_t low(uint64_t node_idx) const noexcept {
        if (!nlow_) return 0;

        const uint64_t bit = node_idx * nlow_;
        const uint64_t lo_word = low_parts_[word_index(bit)], hi_word = low_parts_[word_index(bit) + 1];
        const unsigned shift = bit_offset(bit);
        const uint64_t bits = shift ? (lo_word >> shift) | (hi_word << (kWordBits - shift)) : lo_word;

        return bits & low_mask();
    }
    uint64_t select0(uint64_t rank) const noexcept {          // finds the rank-th zero
        const uint64_t sample_pos = select0_samples_[rank / kSelect0Sample];
        uint64_t rem = rank % kSelect0Sample;
        uint64_t word_idx = word_index(sample_pos);
        uint64_t word = ~high_bits_[word_idx] & (~uint64_t(0) << bit_offset(sample_pos));

        for (;;) {
            const uint64_t set_bits = static_cast<uint64_t>(__builtin_popcountll(word));

            if (rem < set_bits) return first_bit_of_word(word_idx) + select_in_word(word, static_cast<unsigned>(rem));

            rem -= set_bits;
            word = ~high_bits_[++word_idx];
        }
    }
    uint64_t select1(uint64_t rank) const noexcept {          // finds the rank-th one
        const uint64_t sample_pos = select1_samples_[rank / kSelect1Sample];
        uint64_t remaining = rank % kSelect1Sample;
        uint64_t word_idx = word_index(sample_pos);
        uint64_t word = high_bits_[word_idx] & (~uint64_t(0) << bit_offset(sample_pos));

        for (;;) {
            const uint64_t set_bits = static_cast<uint64_t>(__builtin_popcountll(word));

            if (remaining < set_bits) return first_bit_of_word(word_idx) + select_in_word(word, static_cast<unsigned>(remaining));

            remaining -= set_bits;
            word = high_bits_[++word_idx];
        }
    }
    uint32_t& node_word(std::size_t node_idx) const noexcept { return const_cast<uint32_t&>(node_words_[node_idx]); }
    uint32_t node(std::size_t node_idx) const noexcept { return __atomic_load_n(&node_word(node_idx), __ATOMIC_RELAXED); }

    unsigned kmer_len_;
    std::size_t num_nodes_ = 0, cap_ = 0, reserve_ = 0;
    uint64_t nlow_ = 0, last_ = 0, nzeros_ = 0;
    bool consumed_ = false;
    HugeVec<uint64_t> low_parts_, high_bits_;
    HugeVec<uint32_t> node_words_;
    HugeVec<uint64_t> select0_samples_, select1_samples_;       // huge pages here too, for the TLB
    std::vector<std::pair<std::size_t, uint32_t>> saturated_;
    std::mutex saturated_mutex_;
};

}
