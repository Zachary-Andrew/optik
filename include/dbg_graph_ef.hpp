#pragma once

/*
 * dbg_graph_ef.hpp - compacted de Bruijn graph over Elias-Fano coded keys, for k-mer sets whose open-addressing table (DbgGraph, about 16 bytes per k-mer) does not fit in memory.
 * A whole human genome has about 2.5e9 distinct 31-mers, 40 GB of table.
 *
 * Keys.
 * The sorted canonical k-mers as an Elias-Fano sequence (Elias 1974, Fano 1971; select sampling as in Vigna 2013), with l = floor(log2(4^k / n)) low bits per key, packed, and the high parts as a unary bitvector of about 2n bits, so about 2 + log2(4^k / n) bits per k-mer (32 for the human genome at k = 31).
 * rank(c), the node index of k-mer c or kNpos, finds c's high-part bucket with a sampled select0 and compares the bucket's few low parts, giving exact membership with no false positives.
 * key(i) is a sampled select1.
 * Node i is the i-th smallest k-mer, so iteration is in sorted order.
 *
 * Nodes.
 * One 32-bit word per k-mer, laid out as in DbgGraph, holds 22 bits of abundance, with larger ones in a sorted side list, the walker's two flags and the in << 4 | out edge mask.
 *
 * About 8 bytes per k-mer in all, against the table's 16, at 2-3 dependent cache misses per lookup instead of one.
 * The counter chooses it only when the table would not fit the memory limit (src/kmer_count.cpp), so every graph that fits keeps the faster table.
 *
 * The interface is DbgGraph's node interface, so the unitig walker (dbg_unitigs.hpp) and the GFA writer (src/ffi_dbg.cpp) run on either, and their output is identical (bench/spill_check.py compares them).
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

    explicit DbgGraphEF(unsigned kk) : k_(kk) {}

    unsigned    k() const noexcept { return k_; }
    std::size_t core_size() const noexcept { return n_; }
    std::size_t capacity() const noexcept { return cap_; }

    /// bytes_for() - bytes of a graph of @n k-mers, covering keys, node words and select samples.
    static constexpr uint64_t bytes_for(std::size_t n, unsigned k) noexcept {
        if (n == 0) return 0;
        const unsigned l = low_bits(n, k);
        const uint64_t hi = n + ((universe(k) - 1) >> l) + 2;
        return n * uint64_t(l) / 8 + hi / 8 + n * 4 + hi / 8 / kSelect0Sample * 8 + n / kSelect1Sample * 8 + 64;
    }

    // Build with begin(n), then put() for every node in any order from any thread, then finish_put().

    /// begin() - start a build of @n nodes.
    void begin(std::size_t n) {
        n_ = cap_ = 0;
        reserve_ = n;
        nlow_ = n ? low_bits(n, k_) : 0;
        const uint64_t hbits = n + ((universe(k_) - 1) >> nlow_) + 2;
        low_parts_.assign((n * nlow_ + 63) / 64 + 2, 0);
        high_bits_.assign(hbits / 64 + 2, 0);
        node_words_.assign(n, 0);
        select1_samples_.clear();
        saturated_.clear();
        last_ = 0;
    }
    /// put() - node @i gets @key, abundance @cnt and observed edges @raw.
    /// Thread-safe for distinct @i, with the caller guaranteeing that keys increase with i.
    void put(std::size_t i, uint64_t key, uint32_t cnt, uint8_t raw) {
        if (nlow_) {
            const uint64_t lowv = key & ((uint64_t(1) << nlow_) - 1), bit = i * nlow_;
            __atomic_fetch_or(&low_parts_[bit >> 6], lowv << (bit & 63), __ATOMIC_RELAXED);
            if ((bit & 63) + nlow_ > 64) __atomic_fetch_or(&low_parts_[(bit >> 6) + 1], lowv >> (64 - (bit & 63)), __ATOMIC_RELAXED);
        }
        const uint64_t hp = (key >> nlow_) + i;
        __atomic_fetch_or(&high_bits_[hp >> 6], uint64_t(1) << (hp & 63), __ATOMIC_RELAXED);
        if (cnt >= kCountMax) { std::lock_guard<std::mutex> lk(saturated_mutex_); saturated_.push_back({i, cnt}); }
        node_words_[i] = std::min(cnt, kCountMax) | (uint32_t(raw) << kMaskShift);
    }
    /// finish_put() - after put() of every node 0..n-1, take the select1 samples from the high bits and check the key order.
    /// Throws std::logic_error if the keys were not distinct and increasing.
    void finish_put() {
        n_ = reserve_;
        select1_samples_.clear();
        uint64_t ones = 0;
        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {
            const uint64_t w = high_bits_[wi];
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(w));
            for (uint64_t r = (ones + kSelect1Sample - 1) / kSelect1Sample * kSelect1Sample; r < ones + c; r += kSelect1Sample)
                select1_samples_.push_back((uint64_t(wi) << 6) + select_in_word(w, static_cast<unsigned>(r - ones)));
            ones += c;
        }
        if (ones != n_) throw std::logic_error("DbgGraphEF: keys were not distinct and increasing");
        last_ = n_ ? key(n_ - 1) : 0;
        finish();
    }
    /// finish() - take the select0 samples and sort the side list.
    void finish() {
        if (n_ != reserve_) throw std::logic_error("DbgGraphEF: fewer keys appended than reserved");
        cap_ = n_;
        // select0 samples hold the position of every kSelect0Sample-th zero of the high bits
        select0_samples_.clear();
        uint64_t zeros = 0;
        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {
            const uint64_t z = ~high_bits_[wi];
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(z));
            for (uint64_t r = (zeros + kSelect0Sample - 1) / kSelect0Sample * kSelect0Sample; r < zeros + c; r += kSelect0Sample)
                select0_samples_.push_back((uint64_t(wi) << 6) + select_in_word(z, static_cast<unsigned>(r - zeros)));
            zeros += c;
        }
        nzeros_ = zeros;
        std::sort(saturated_.begin(), saturated_.end());
    }

    /// rank() - node index of canonical k-mer @c, or kNpos.
    std::size_t rank(uint64_t c) const noexcept {
        if (n_ == 0 || c > last_) return kNpos;
        const uint64_t h = c >> nlow_;
        uint64_t pos;
        if (h == 0) pos = 0;
        else {
            if (h - 1 >= nzeros_) return kNpos;
            pos = select0(h - 1) + 1;
        }
        uint64_t idx = pos - h;
        const uint64_t lc = c & low_mask();
        while (idx < n_ && ((high_bits_[pos >> 6] >> (pos & 63)) & 1)) {
            const uint64_t lw = low(idx);
            if (lw == lc) return idx;
            if (lw > lc) return kNpos;
            ++idx; ++pos;
        }
        return kNpos;
    }
    uint64_t key(std::size_t i) const noexcept { return ((select1(i) - i) << nlow_) | low(i); }
    bool used(std::size_t i) const noexcept { return i < n_; }
    /// prefetch() - start the first load rank(@c) needs, its select0 sample.
    /// The high and low words depend on that sample's value, and estimating their position from c alone is off by about sqrt(n) bits, so they are not prefetched.
    void prefetch(uint64_t c) const noexcept {
        const uint64_t h = c >> nlow_;
        if (h && (h - 1) / kSelect0Sample < select0_samples_.size()) __builtin_prefetch(&select0_samples_[(h - 1) / kSelect0Sample]);
    }
    void prefetch_w(uint64_t c) const noexcept { prefetch(c); }

    unsigned edge_mask(std::size_t i) const noexcept { return node(i) >> kMaskShift; }
    uint64_t count(std::size_t i) const noexcept {
        const uint32_t c = node(i) & kCountMax;
        if (c < kCountMax) return c;
        const auto it = std::lower_bound(saturated_.begin(), saturated_.end(), std::pair<std::size_t, uint32_t>{i, 0});
        return it != saturated_.end() && it->first == i ? it->second : c;
    }
    unsigned flags_or(std::size_t i, unsigned f) const noexcept {
        return (__atomic_fetch_or(&node_word(i), uint32_t(f) << kFlagShift, __ATOMIC_ACQ_REL) >> kFlagShift) & kWalkFlags;
    }
    unsigned flags(std::size_t i) const noexcept { return (node(i) >> kFlagShift) & kWalkFlags; }
    void clear_flags(std::size_t i) const noexcept { node_word(i) &= ~(uint32_t(kWalkFlags) << kFlagShift); }

    /// aux_clear(), aux_set(), aux() - reuse of the node words after the last walk, as DbgGraph::aux_clear() describes.
    void aux_clear(unsigned nthreads) {
        consumed_ = true;
        const unsigned threads = nthreads ? nthreads : 1;
        tpoptoa::run_on_threads(threads, [&](unsigned t) {
            std::fill(node_words_.begin() + n_ * t / threads, node_words_.begin() + n_ * (t + 1) / threads, 0u);
        });
    }
    void aux_set(std::size_t i, uint32_t v) noexcept { __atomic_store_n(&node_word(i), v, __ATOMIC_RELAXED); }
    uint32_t aux(std::size_t i) const noexcept { return node(i); }
    bool consumed() const noexcept { return consumed_; }

    unsigned out_mask(uint64_t c) const noexcept {
        const std::size_t r = rank(c);
        return r == kNpos ? 0u : edge_mask(r) & 0xFu;
    }
    unsigned in_mask(uint64_t c) const noexcept {
        const std::size_t r = rank(c);
        return r == kNpos ? 0u : edge_mask(r) >> 4;
    }
    void masks_from_raw() noexcept {}

    /// resolve_masks() - restrict every node's observed edges to present neighbours, in place.
    /// Each node reads only its own mask and its neighbours' keys.
    void resolve_masks() {
        const unsigned K = k_;
        const uint64_t kmask = kmer_mask(K);
        const unsigned shK = 2u * (K - 1);
        auto canon = [&](uint64_t v) { const uint64_t r = reverse_complement(v, K); return v < r ? v : r; };
        const int64_t N = static_cast<int64_t>(n_);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, kResolveChunk)
#endif
        for (int64_t ii = 0; ii < N; ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            const uint64_t c = key(i);
            const unsigned raw = edge_mask(i);
            unsigned m = 0;
            for (unsigned b = 0; b < 4; ++b) {
                if ((raw & (1u << b)) && rank(canon(((c << 2) | b) & kmask)) != kNpos) m |= 1u << b;
                if ((raw & (1u << (b + 4))) && rank(canon(((uint64_t(b) << shK) | (c >> 2)) & kmask)) != kNpos)
                    m |= 1u << (b + 4);
            }
            node_word(i) = (node_word(i) & ((uint32_t(1) << kMaskShift) - 1)) | (uint32_t(m) << kMaskShift);
        }
    }

    template <typename CB>
    void unitigs(CB&& cb, unsigned nthreads = 0) const {
        DbgUnitigExtractor::run_masks(*this, std::forward<CB>(cb), nthreads);
    }
    template <typename Sink>
    void unitigs_sorted(Sink&& sink, unsigned nthreads = 0) const {
        DbgUnitigExtractor::run_sorted(*this, std::forward<Sink>(sink), nthreads);
    }

    /// export_sorted() - (k-mer, abundance) pairs and edge masks in key order, which node order already is.
    void export_sorted(std::vector<std::pair<uint64_t, uint32_t>>& kv, std::vector<uint8_t>& masks) const {
        if (consumed_) throw std::logic_error("DbgGraphEF: node words were reused (aux_clear)");
        kv.resize(n_);
        masks.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            kv[i] = {key(i), static_cast<uint32_t>(count(i))};
            masks[i] = static_cast<uint8_t>(edge_mask(i));
        }
    }
    /// save() - the OKKDBG archive of DbgGraph::save(), streamed from the sorted keys.
    bool save(const char* path) const {
        if (consumed_) return false;
        return atomic_write(path, [&](std::ofstream& os) {
            os.write("OKKDBG", 6);
            os.put(static_cast<char>(k_));
            const uint64_t n = n_;
            os.write(reinterpret_cast<const char*>(&n), sizeof(n));
            std::vector<char> buf;
            auto flush = [&] { os.write(buf.data(), static_cast<std::streamsize>(buf.size())); buf.clear(); };
            for (std::size_t i = 0; i < n_; ++i) {
                const uint64_t x = key(i);
                buf.insert(buf.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 8);
                if (buf.size() >= kSaveBufferBytes) flush();
            }
            for (std::size_t i = 0; i < n_; ++i) {
                const uint32_t x = static_cast<uint32_t>(count(i));
                buf.insert(buf.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 4);
                if (buf.size() >= kSaveBufferBytes) flush();
            }
            for (std::size_t i = 0; i < n_; ++i) {
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

    /// universe() - 4^k, the number of distinct k-mers (saturating at k = 32).

    static constexpr uint64_t universe(unsigned k) noexcept { return k >= 32 ? ~uint64_t(0) : uint64_t(1) << (2 * k); }
    /// low_bits() - Elias-Fano low-part width floor(log2(universe / n)), at most 63.
    static constexpr unsigned low_bits(std::size_t n, unsigned k) noexcept {
        unsigned l = 0;
        const uint64_t q = universe(k) / n;
        while (l < 63 && (uint64_t(2) << l) <= q) ++l;
        return l;
    }
    uint64_t low_mask() const noexcept { return nlow_ ? (uint64_t(1) << nlow_) - 1 : 0; }
    /// low() - low part of key @i.
    uint64_t low(uint64_t i) const noexcept {
        if (!nlow_) return 0;
        const uint64_t bit = i * nlow_;
        const uint64_t w0 = low_parts_[bit >> 6], w1 = low_parts_[(bit >> 6) + 1];
        const unsigned sh = bit & 63;
        const uint64_t v = sh ? (w0 >> sh) | (w1 << (64 - sh)) : w0;
        return v & low_mask();
    }
    /// select_in_word() - position of the @r-th set bit of @w (0-based), by one PDEP with BMI2, else clearing r low set bits.
    static unsigned select_in_word(uint64_t w, unsigned r) noexcept {
#if defined(__BMI2__)
        return static_cast<unsigned>(__builtin_ctzll(_pdep_u64(uint64_t(1) << r, w)));
#else
        for (unsigned t = 0; t < r; ++t) w &= w - 1;
        return static_cast<unsigned>(__builtin_ctzll(w));
#endif
    }
    uint64_t select0(uint64_t r) const noexcept {          // position of the r-th zero
        const uint64_t sp = select0_samples_[r / kSelect0Sample];
        uint64_t rem = r % kSelect0Sample;
        uint64_t wi = sp >> 6;
        uint64_t w = ~high_bits_[wi] & (~uint64_t(0) << (sp & 63));
        for (;;) {
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(w));
            if (rem < c) return (wi << 6) + select_in_word(w, static_cast<unsigned>(rem));
            rem -= c;
            w = ~high_bits_[++wi];
        }
    }
    uint64_t select1(uint64_t i) const noexcept {          // position of the i-th one
        const uint64_t sp = select1_samples_[i / kSelect1Sample];
        uint64_t r = i % kSelect1Sample;
        uint64_t wi = sp >> 6;
        uint64_t w = high_bits_[wi] & (~uint64_t(0) << (sp & 63));
        for (;;) {
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(w));
            if (r < c) return (wi << 6) + select_in_word(w, static_cast<unsigned>(r));
            r -= c;
            w = high_bits_[++wi];
        }
    }
    uint32_t& node_word(std::size_t i) const noexcept { return const_cast<uint32_t&>(node_words_[i]); }
    uint32_t node(std::size_t i) const noexcept { return __atomic_load_n(&node_word(i), __ATOMIC_RELAXED); }

    unsigned k_;
    std::size_t n_ = 0, cap_ = 0, reserve_ = 0;
    uint64_t nlow_ = 0, last_ = 0, nzeros_ = 0;
    bool consumed_ = false;
    HugeVec<uint64_t> low_parts_, high_bits_;
    HugeVec<uint32_t> node_words_;
    HugeVec<uint64_t> select0_samples_, select1_samples_;       // select samples, on huge pages too (TLB)
    std::vector<std::pair<std::size_t, uint32_t>> saturated_;
    std::mutex saturated_mutex_;
};

}
