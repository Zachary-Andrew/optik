#pragma once

/*
 * phobic_mphf.hpp - minimal perfect hash functions over 64-bit keys, and the Elias-Fano sequence they use.
 *
 * The sources are PTHash (Pibiri and Trani 2021) and PHOBIC (Hermann et al. 2024).
 * The keys are split into many small buckets, and every bucket gets a "pilot" p such that all of its keys land on free slots of a table of size T = n / alpha via
 *
 *     pos(key) = fastrange(mix(h(key) ^ mix(p)), T)
 *
 * A lookup is one hash, one pilot read and one position computation.
 * The few slots at or past n that alpha < 1 leaves are remapped into the holes below n.
 * fastrange replaces modulo everywhere (Lemire 2019).
 *
 * Adaptations in optik
 *   - Skewed bucket assignment as in PTHash sends 60% of the keys to the first 30% of the buckets, the dense ones.
 *   - Dense and sparse buckets have far apart pilot ranges, so each part stores its pilots at its own fixed bit width, read with one unaligned 64-bit load.
 *   - The remap of the alpha slack is non-decreasing, so EfSeq codes it as Elias-Fano, which makes a low alpha cheap, with fewer pilot trials and smaller pilots.
 *   - Keys are expected to be well mixed already, as hashed minimizers are.
 *     They are re-hashed with a per-build seed so a failed build can reseed.
 *   - Mphf partitions the keys and builds one PilotMphf per partition on its own thread, as PHOBIC's partitioned construction does.
 *
 * The pilot of a bucket only locates keys together with their hash, which relates it to the tiny pointers of Bender et al. 2025.
 *
 * A key outside the build set gets an arbitrary value in [0, n), so callers must verify, as the SPSS core does by comparing against its string.
 */

#include "thread_budget.hpp"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <thread>
#include <vector>

#if defined(__BMI2__)
#include <immintrin.h>
#endif

namespace tpoptoa {

/// WordView - a read-only run of 64-bit words held elsewhere, by a vector of the owning object after a build, or by a mapped index file after a borrowing load, which then costs no copy and no second resident set.
struct WordView {
    const uint64_t* p = nullptr;
    std::size_t n = 0;
    const uint64_t& operator[](std::size_t i) const noexcept { return p[i]; }
    std::size_t size() const noexcept { return n; }
    bool empty() const noexcept { return n == 0; }
    const uint64_t* data() const noexcept { return p; }
    const uint64_t* begin() const noexcept { return p; }
    const uint64_t* end() const noexcept { return p + n; }
    static WordView of(const std::vector<uint64_t>& v) noexcept { return {v.data(), v.size()}; }
};

/// load_words() - read a u64 count and that many words at @p + @o into @view, pointing into @p when @borrow (so @p must outlive the view and be 8-byte aligned) and otherwise copying into @own.
/// Returns false past @sz.
inline bool load_words(const uint8_t* p, std::size_t sz, std::size_t& o, bool borrow,
                       std::vector<uint64_t>& own, WordView& view) {
    uint64_t m;
    if (o + 8 > sz) return false;
    std::memcpy(&m, p + o, 8);
    o += 8;
    if (m > (sz - o) / 8) return false;
    if (borrow) {
        own.clear();
        own.shrink_to_fit();
        view = {reinterpret_cast<const uint64_t*>(p + o), static_cast<std::size_t>(m)};
    } else {
        own.resize(m);
        if (m) std::memcpy(own.data(), p + o, m * 8);
        view = WordView::of(own);
    }
    o += m * 8;
    return true;
}

/// Elias-Fano sequence of non-decreasing integers with O(1) access through sampled select (Elias 1974, Fano 1971; select sampling as in Vigna, Quasi-succinct indices, WSDM 2013).
/// Every S-th set bit of the high bitvector has its position stored.
/// access(i) jumps to the sample, popcount-scans whole words and finishes with an in-word select, using BMI2 pdep and tzcnt when available, as algorithmica.org/hpc "bit manipulation" describes.
class EfSeq {
public:
    static constexpr uint64_t kSampleRate = 256;

    /// build() - encode the non-decreasing @v, every value below @universe.
    void build(const std::vector<uint64_t>& v, uint64_t universe) {
        n_ = v.size();
        nlow_ = 0;
        if (n_ && universe / n_ > 1) while ((uint64_t(2) << nlow_) <= universe / n_) ++nlow_;
        low_own_.assign((n_ * nlow_ + 63) / 64 + 2, 0);
        const uint64_t hbits = n_ + (universe >> nlow_) + 2;
        high_own_.assign(hbits / 64 + 2, 0);
        one_own_.clear();
        for (uint64_t i = 0; i < n_; ++i) {
            if (nlow_) {
                const uint64_t lowv = v[i] & ((uint64_t(1) << nlow_) - 1), bit = i * nlow_;
                low_own_[bit >> 6] |= lowv << (bit & 63);
                if ((bit & 63) + nlow_ > 64) low_own_[(bit >> 6) + 1] |= lowv >> (64 - (bit & 63));
            }
            const uint64_t hp = (v[i] >> nlow_) + i;
            high_own_[hp >> 6] |= uint64_t(1) << (hp & 63);
            if (i % kSampleRate == 0) one_own_.push_back(hp);
        }
        borrowed_ = false;
        adopt_owned();
        build_zero_samples();
    }

    /// count_le() - number of values <= @x.
    ///
    /// The bucket of x's high part starts after the (x >> l)-th zero of the high bits, found by sampled select0, and its few ones are compared on their low bits in O(1) expected time instead of a binary search over operator[].
    uint64_t count_le(uint64_t x) const noexcept {
        if (n_ == 0) return 0;
        const uint64_t h = x >> nlow_;
        uint64_t pos;                                   // first position of bucket h
        if (h == 0) pos = 0;
        else {
            if (h - 1 >= nzeros_) return n_;            // past every bucket
            pos = select0(h - 1) + 1;
        }
        uint64_t idx = pos - h;                         // values with high part < h
        const uint64_t lowx = x & ((uint64_t(1) << nlow_) - 1);
        while (idx < n_ && ((high_bits_[pos >> 6] >> (pos & 63)) & 1)) {
            if (low(idx) > lowx) break;
            ++idx; ++pos;
        }
        return idx;
    }

    /// operator[]() - the @i-th value.
    uint64_t operator[](uint64_t i) const noexcept {
        return ((select1(i) - i) << nlow_) | low(i);
    }

    std::size_t bytes() const noexcept {
        return low_parts_.size() * 8 + high_bits_.size() * 8 + one_samples_.size() * 8 + 16;
    }
    uint64_t size() const noexcept { return n_; }
    std::size_t mem_bytes() const noexcept { return bytes() + zero_samples_.size() * 8; }
    /// Structural check after load that the arrays can hold n_ values of l_ low bits and every sample points into the high bits (a truncated or corrupted blob is rejected instead of read out of bounds).
    bool sane() const noexcept {
        if (nlow_ > 63 || low_parts_.size() < (n_ * nlow_ + 63) / 64 + 1 || one_samples_.size() < (n_ + kSampleRate - 1) / kSampleRate) return false;
        if (n_ && high_bits_.empty()) return false;
        for (uint64_t sp : one_samples_) if ((sp >> 6) >= high_bits_.size()) return false;
        uint64_t ones = 0;
        for (uint64_t w : high_bits_) ones += static_cast<uint64_t>(__builtin_popcountll(w));
        return ones == n_;
    }
    /// save() and load() - little-endian serialization.
    /// load() returns the bytes consumed, or 0 for a malformed blob.
    /// The select0 samples are rebuilt rather than stored.
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t v) { uint8_t b[8]; std::memcpy(b, &v, 8); out.insert(out.end(), b, b + 8); };
        put(n_); put(nlow_);
        for (const auto* vec : {&low_parts_, &high_bits_, &one_samples_}) { put(vec->size()); for (uint64_t x : *vec) put(x); }
    }
    /// With @borrow the three arrays are used in place, so @p must stay mapped, and 8-byte aligned, for the lifetime of this object.
    std::size_t load(const uint8_t* p, std::size_t sz, bool borrow = false) {
        std::size_t o = 0;
        auto get = [&](uint64_t& v) { if (o + 8 > sz) return false; std::memcpy(&v, p + o, 8); o += 8; return true; };
        if (!get(n_) || !get(nlow_)) return 0;
        if (!load_words(p, sz, o, borrow, low_own_, low_parts_) || !load_words(p, sz, o, borrow, high_own_, high_bits_)
            || !load_words(p, sz, o, borrow, one_own_, one_samples_)) return 0;
        borrowed_ = borrow;
        build_zero_samples();                            // derived, not stored
        return o;
    }

    EfSeq() = default;
    EfSeq(EfSeq&&) noexcept = default;               // moved vectors keep their buffers
    EfSeq& operator=(EfSeq&&) noexcept = default;
    /// A copy owns copies of the arrays and points its views at them, unless the source borrows a mapping, whose views stay valid as they are.
    EfSeq(const EfSeq& o) { *this = o; }
    EfSeq& operator=(const EfSeq& o) {
        if (this == &o) return *this;
        n_ = o.n_; nlow_ = o.nlow_; nzeros_ = o.nzeros_; borrowed_ = o.borrowed_;
        low_own_ = o.low_own_; high_own_ = o.high_own_; one_own_ = o.one_own_;
        zero_samples_ = o.zero_samples_;
        if (borrowed_) { low_parts_ = o.low_parts_; high_bits_ = o.high_bits_; one_samples_ = o.one_samples_; }
        else adopt_owned();
        return *this;
    }

private:
    static unsigned select_in_word(uint64_t w, unsigned r) noexcept {   // r-th set bit
#if defined(__BMI2__)
        return static_cast<unsigned>(__builtin_ctzll(_pdep_u64(uint64_t(1) << r, w)));
#else
        for (unsigned t = 0; t < r; ++t) w &= w - 1;
        return static_cast<unsigned>(__builtin_ctzll(w));
#endif
    }
    uint64_t low(uint64_t i) const noexcept {
        if (!nlow_) return 0;
        const uint64_t bit = i * nlow_;
        uint64_t w;
        std::memcpy(&w, reinterpret_cast<const uint8_t*>(low_parts_.data()) + (bit >> 3), 8);
        return (w >> (bit & 7)) & ((uint64_t(1) << nlow_) - 1);
    }
    /// zero_samples_[t] = position of the (t * S)-th zero of high_bits_ (over the positions that can hold buckets, [0, n + (max value >> l) + 1)).
    void build_zero_samples() {
        zero_samples_.clear();
        nzeros_ = 0;
        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {   // a word at a time
            const uint64_t z = ~high_bits_[wi];
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(z));
            // next sampled rank, a multiple of S, that falls in this word
            for (uint64_t r = (nzeros_ + kSampleRate - 1) / kSampleRate * kSampleRate; r < nzeros_ + c; r += kSampleRate)
                zero_samples_.push_back((uint64_t(wi) << 6)
                                    + select_in_word(z, static_cast<unsigned>(r - nzeros_)));
            nzeros_ += c;
        }
    }
    uint64_t select0(uint64_t r) const noexcept {        // position of the r-th zero
        const uint64_t sp = zero_samples_[r / kSampleRate];
        uint64_t rem = r % kSampleRate;
        uint64_t wi = sp >> 6;
        uint64_t w = ~high_bits_[wi] & (~uint64_t(0) << (sp & 63));
        for (;;) {
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(w));
            if (rem < c) return (wi << 6) + select_in_word(w, static_cast<unsigned>(rem));
            rem -= c;
            w = ~high_bits_[++wi];
        }
    }
    uint64_t select1(uint64_t i) const noexcept {
        const uint64_t sp = one_samples_[i / kSampleRate];
        uint64_t r = i % kSampleRate;                          // ones still to skip after sp
        uint64_t wi = sp >> 6;
        uint64_t w = high_bits_[wi] & (~uint64_t(0) << (sp & 63));
        for (;;) {
            const uint64_t c = static_cast<uint64_t>(__builtin_popcountll(w));
            if (r < c) return (wi << 6) + select_in_word(w, static_cast<unsigned>(r));
            r -= c;
            w = high_bits_[++wi];
        }
    }
    uint64_t n_ = 0, nlow_ = 0;
    // The arrays every lookup reads, as views of the owned vectors after a build or of a mapped file after a borrowing load.
    std::vector<uint64_t> low_own_, high_own_, one_own_;
    WordView low_parts_, high_bits_, one_samples_;
    bool borrowed_ = false;
    std::vector<uint64_t> zero_samples_;                    // select0 samples (derived)
    void adopt_owned() noexcept {
        low_parts_ = WordView::of(low_own_);
        high_bits_ = WordView::of(high_own_);
        one_samples_ = WordView::of(one_own_);
    }
    uint64_t nzeros_ = 0;
};

/// PilotMphf - one minimal perfect hash function mapping n keys to [0, n).
class PilotMphf {
public:
    PilotMphf() = default;

    /// build() - build over @n distinct keys, returning false on duplicate keys or when kMaxBuildAttempts seeds all fail.
    ///
    /// @c sets the bucket count, c * n / log2(n), and @alpha is the load of the table before remapping.
    /// The defaults c = 6 and alpha = 0.85 gave the smallest total size with the split pilot widths and Elias-Fano remap when they were chosen, on 0.54M random keys, at 4.7 bits per key against 6.9 for the earlier c = 8, alpha = 0.97 with a 32-bit remap, at the same build time.
    /// Mphf passes its own alpha of 0.94.
    bool build(const uint64_t* keys, std::size_t n, double c = 6.0,
               double alpha = 0.85) {
        n_ = n;
        dense_own_.clear(); sparse_own_.clear();
        borrowed_ = false;
        adopt_owned();
        if (n == 0) { table_size_ = 0; nbuckets_ = 0; dense_width_ = sparse_width_ = 0; return true; }
        for (uint64_t attempt = 0; attempt < kMaxBuildAttempts; ++attempt) {
            seed_ = kSeedStep * (attempt + 1) ^ kSeedBase;
            if (try_build(keys, n, c, alpha)) return true;
        }
        return false;
    }

    /// operator()() - rank in [0, n) of a member @key.
    uint64_t operator()(uint64_t key) const noexcept {
        const uint64_t h = mix(key ^ seed_);
        const uint64_t p = pilot(bucket_of(h));
        const uint64_t pos = fastrange(mix(h ^ mix(p + kPilotSalt)), table_size_);
        return pos < n_ ? pos : remap_[pos - n_];         // Elias-Fano access
    }

    /// prefetch() - prefetch the pilot a later operator()(@key) will read.
    void prefetch(uint64_t key) const noexcept {
        const uint64_t h = mix(key ^ seed_);
        const uint64_t b = bucket_of(h);
        if (b < ndense_) __builtin_prefetch(dense_pilots_.data() + ((b * dense_width_) >> 6), 0, 1);
        else __builtin_prefetch(sparse_pilots_.data() + (((b - ndense_) * sparse_width_) >> 6), 0, 1);
    }

    std::size_t size() const noexcept { return n_; }
    std::size_t bytes() const noexcept {
        return dense_pilots_.size() * 8 + sparse_pilots_.size() * 8 + remap_.bytes() + 64;
    }

    /// save() - little-endian, self-delimited serialization.
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t v) { uint8_t b[8]; std::memcpy(b, &v, 8); out.insert(out.end(), b, b + 8); };
        put(n_); put(table_size_); put(nbuckets_); put(ndense_); put(seed_); put(dense_width_); put(sparse_width_);
        put(dense_pilots_.size()); for (uint64_t w : dense_pilots_) put(w);
        put(sparse_pilots_.size()); for (uint64_t w : sparse_pilots_) put(w);
        remap_.save(out);
    }
    /// load() - returns the bytes consumed, or 0 for a malformed blob.
    std::size_t load(const uint8_t* p, std::size_t sz, bool borrow = false) {
        std::size_t o = 0;
        auto get = [&](uint64_t& v) { if (o + 8 > sz) return false; std::memcpy(&v, p + o, 8); o += 8; return true; };
        if (!get(n_) || !get(table_size_) || !get(nbuckets_) || !get(ndense_) || !get(seed_) || !get(dense_width_) || !get(sparse_width_)) return 0;
        if (!load_words(p, sz, o, borrow, dense_own_, dense_pilots_)
            || !load_words(p, sz, o, borrow, sparse_own_, sparse_pilots_)) return 0;
        borrowed_ = borrow;
        const std::size_t used = remap_.load(p + o, sz - o, borrow);
        if (!used) return 0;
        return o + used;
    }

private:
    static constexpr uint64_t kMaxBuildAttempts = 16;   // seeds tried before build() gives up
    static constexpr uint64_t kSeedStep = 0x9E3779B97F4A7C15ULL;
    static constexpr uint64_t kSeedBase = 0xD6E8FEB86659FD93ULL;
    static constexpr uint64_t kPilotSalt = 0x5851F42D4C957F2DULL;   // added to a pilot before mixing
    static constexpr uint64_t kMaxPilot = uint64_t(1) << 24;         // a bucket needing more reseeds
    static constexpr double kDenseKeyShare = 0.6;      // share of keys in the dense buckets
    static constexpr double kDenseBucketShare = 0.3;   // share of buckets that are dense

    static constexpr uint64_t mix(uint64_t x) noexcept {           // MurmurHash3 fmix64
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;
        return x;
    }
    static constexpr uint64_t fastrange(uint64_t h, uint64_t n) noexcept {
        return static_cast<uint64_t>((static_cast<__uint128_t>(h) * n) >> 64);
    }
    /// bucket_of() - skewed split, where kDenseKeyShare of the hash space goes to the first ndense_ (kDenseBucketShare) buckets.
    uint64_t bucket_of(uint64_t h) const noexcept {
        static constexpr uint64_t kDenseHashCutoff = static_cast<uint64_t>(kDenseKeyShare * 18446744073709551615.0);
        return h < kDenseHashCutoff ? fastrange(h, ndense_) : ndense_ + fastrange(h, nbuckets_ - ndense_);
    }
    /// read_packed() and write_packed() - element @idx of a packed array of
    /// @w-bit values (one unaligned 64-bit load to read).
    static uint64_t read_packed(const WordView& v, uint64_t idx, uint64_t w) noexcept {
        if (w == 0) return 0;
        const uint64_t bit = idx * w;
        uint64_t x;
        std::memcpy(&x, reinterpret_cast<const uint8_t*>(v.data()) + (bit >> 3), 8);
        return (x >> (bit & 7)) & ((uint64_t(1) << w) - 1);
    }
    static void write_packed(std::vector<uint64_t>& v, uint64_t idx, uint64_t w, uint64_t val) noexcept {
        if (w == 0) return;
        const uint64_t bit = idx * w;
        v[bit >> 6] |= val << (bit & 63);
        if ((bit & 63) + w > 64) v[(bit >> 6) + 1] |= val >> (64 - (bit & 63));
    }
    uint64_t pilot(uint64_t b) const noexcept {
        return b < ndense_ ? read_packed(dense_pilots_, b, dense_width_) : read_packed(sparse_pilots_, b - ndense_, sparse_width_);
    }

    /// try_build() - one build attempt with the current seed_, which fails when a duplicate turns up or a bucket needs a pilot above kMaxPilot.
    bool try_build(const uint64_t* keys, std::size_t n, double c, double alpha) {
        table_size_ = std::max<uint64_t>(n, static_cast<uint64_t>(std::ceil(n / alpha)));
        const double lg = std::max(1.0, std::log2(static_cast<double>(n)));
        nbuckets_ = std::max<uint64_t>(1, static_cast<uint64_t>(std::ceil(c * n / lg)));
        ndense_ = std::max<uint64_t>(1, static_cast<uint64_t>(kDenseBucketShare * nbuckets_));
        if (ndense_ >= nbuckets_) ndense_ = nbuckets_ - 1 ? nbuckets_ - 1 : 1;
        if (nbuckets_ == 1) ndense_ = 1;

        // hashes, grouped by bucket with a counting sort
        std::vector<uint64_t> h(n);
        for (std::size_t i = 0; i < n; ++i) h[i] = mix(keys[i] ^ seed_);
        std::vector<uint32_t> cnt(nbuckets_ + 1, 0);
        std::vector<uint64_t> bucket_ids(n);
        for (std::size_t i = 0; i < n; ++i) { bucket_ids[i] = bucket_of(h[i]); ++cnt[bucket_ids[i] + 1]; }
        for (uint64_t b = 0; b < nbuckets_; ++b) cnt[b + 1] += cnt[b];
        std::vector<uint64_t> grouped(n);
        {
            std::vector<uint32_t> w(cnt.begin(), cnt.end() - 1);
            for (std::size_t i = 0; i < n; ++i) grouped[w[bucket_ids[i]]++] = h[i];
        }
        // duplicate keys or hashes can only meet in one bucket, so check there (buckets hold ~log2(n)/c keys) instead of sorting all n hashes
        uint32_t max_sz = 0;
        for (uint64_t b = 0; b < nbuckets_; ++b) {
            const uint32_t lo = cnt[b], hi = cnt[b + 1];
            max_sz = std::max(max_sz, hi - lo);
            for (uint32_t x = lo + 1; x < hi; ++x)
                for (uint32_t y = lo; y < x; ++y)
                    if (grouped[x] == grouped[y]) return false;
        }
        // Largest buckets first, by a counting sort on size that is stable, so ties keep bucket order as the former stable_sort did
        std::vector<uint32_t> order(nbuckets_);
        {
            std::vector<uint32_t> start(max_sz + 2, 0);
            for (uint64_t b = 0; b < nbuckets_; ++b) ++start[max_sz - (cnt[b + 1] - cnt[b]) + 1];
            for (uint32_t z = 0; z <= max_sz; ++z) start[z + 1] += start[z];
            for (uint64_t b = 0; b < nbuckets_; ++b)
                order[start[max_sz - (cnt[b + 1] - cnt[b])]++] = static_cast<uint32_t>(b);
        }

        std::vector<uint64_t> taken((table_size_ + 63) / 64, 0);
        std::vector<uint64_t> pilots(nbuckets_, 0);
        std::vector<uint64_t> pos;
        uint64_t max_pilot = 0;
        for (uint32_t b : order) {
            const uint32_t lo = cnt[b], hi = cnt[b + 1];
            if (lo == hi) continue;
            pos.resize(hi - lo);
            uint64_t p = 0;
            for (;; ++p) {
                if (p > kMaxPilot) return false;                   // reseed
                const uint64_t mp = mix(p + kPilotSalt);
                bool ok = true;
                for (uint32_t j = lo; j < hi && ok; ++j) {
                    const uint64_t q = fastrange(mix(grouped[j] ^ mp), table_size_);
                    if ((taken[q >> 6] >> (q & 63)) & 1) { ok = false; break; }
                    for (uint32_t t = lo; t < j; ++t) if (pos[t - lo] == q) { ok = false; break; }
                    pos[j - lo] = q;
                }
                if (ok) break;
            }
            for (uint32_t j = lo; j < hi; ++j) taken[pos[j - lo] >> 6] |= uint64_t(1) << (pos[j - lo] & 63);
            pilots[b] = p;
            max_pilot = std::max(max_pilot, p);
        }
        (void)max_pilot;
        // pilots in two fixed-width parts, plus padding for the unaligned load
        auto width_of = [](uint64_t m) { uint64_t w = 0; while ((uint64_t(1) << w) <= m) ++w; return w; };
        uint64_t m1 = 0, m2 = 0;
        for (uint64_t b = 0; b < nbuckets_; ++b) (b < ndense_ ? m1 : m2) = std::max(b < ndense_ ? m1 : m2, pilots[b]);
        dense_width_ = width_of(m1); sparse_width_ = width_of(m2);
        dense_own_.assign((ndense_ * dense_width_ + 63) / 64 + 2, 0);
        sparse_own_.assign(((nbuckets_ - ndense_) * sparse_width_ + 63) / 64 + 2, 0);
        for (uint64_t b = 0; b < nbuckets_; ++b) {
            if (b < ndense_) write_packed(dense_own_, b, dense_width_, pilots[b]);
            else write_packed(sparse_own_, b - ndense_, sparse_width_, pilots[b]);
        }
        borrowed_ = false;
        adopt_owned();
        // Remap the slack slots >= n onto the holes < n.
        // Holes are assigned in increasing order, and untaken slots repeat the previous value since they are never queried, so the table is non-decreasing and suits Elias-Fano.
        std::vector<uint64_t> remap_values(table_size_ - n_, 0);
        uint64_t hole = 0, last = 0;
        auto is_taken = [&](uint64_t q) { return (taken[q >> 6] >> (q & 63)) & 1; };
        for (uint64_t q = n_; q < table_size_; ++q) {
            if (is_taken(q)) {
                while (is_taken(hole)) ++hole;
                last = hole++;
            }
            remap_values[q - n_] = last;
        }
        remap_.build(remap_values, n_);
        return true;
    }

    uint64_t n_ = 0, table_size_ = 0, nbuckets_ = 0, ndense_ = 1, seed_ = 0, dense_width_ = 0, sparse_width_ = 0;
    // pilots of the dense buckets [0, ndense_) and of the rest, each packed at its own width, read through views of the owned vectors after a build or of a mapped file after a borrowing load
    std::vector<uint64_t> dense_own_, sparse_own_;
    WordView dense_pilots_, sparse_pilots_;
    bool borrowed_ = false;
    EfSeq remap_;
    void adopt_owned() noexcept {
        dense_pilots_ = WordView::of(dense_own_);
        sparse_pilots_ = WordView::of(sparse_own_);
    }

public:
    PilotMphf(PilotMphf&&) noexcept = default;       // moved vectors keep their buffers
    PilotMphf& operator=(PilotMphf&&) noexcept = default;
    /// A copy owns copies of the pilots and points its views at them, unless the source borrows a mapping, whose views stay valid as they are.
    PilotMphf(const PilotMphf& o) { *this = o; }
    PilotMphf& operator=(const PilotMphf& o) {
        if (this == &o) return *this;
        n_ = o.n_; table_size_ = o.table_size_; nbuckets_ = o.nbuckets_; ndense_ = o.ndense_;
        seed_ = o.seed_; dense_width_ = o.dense_width_; sparse_width_ = o.sparse_width_;
        dense_own_ = o.dense_own_; sparse_own_ = o.sparse_own_; borrowed_ = o.borrowed_;
        remap_ = o.remap_;
        if (borrowed_) { dense_pilots_ = o.dense_pilots_; sparse_pilots_ = o.sparse_pilots_; }
        else adopt_owned();
        return *this;
    }
};

/// Mphf - partitioned minimal perfect hash (the partitioned construction of PHOBIC and partitioned PTHash).
///
/// Keys are split by an independent hash into partitions of about kPartKeys keys, and the partitions build their PilotMphfs in parallel, since the pilot search is sequential within one.
/// rank = offset[part] + local rank.
/// The partition table is small enough to stay in cache, so a lookup is still one pilot access, and each build's working set (taken bitmap, bucket arrays) stays about L2-sized.
class Mphf {
public:
    static constexpr std::size_t kPartKeys = 250000;

    /// build() - build over @n distinct keys on @nthreads threads, where 0 means resolve_threads(0).
    /// Returns false on duplicate keys.
    bool build(const uint64_t* keys, std::size_t n, double c = 6.0, double alpha = 0.94,
               unsigned nthreads = 0) {
        n_ = n;
        parts_.clear(); part_offset_.clear();
        nparts_ = std::max<std::size_t>(1, (n + kPartKeys - 1) / kPartKeys);
        if (n == 0) { part_offset_.assign(2, 0); parts_.resize(1); return parts_[0].build(keys, 0, c, alpha); }
        std::vector<std::vector<uint64_t>> part_keys(nparts_);
        for (std::size_t i = 0; i < n; ++i) part_keys[part_of(keys[i])].push_back(keys[i]);
        return build_parts(part_keys, c, alpha, nthreads);
    }

    /// build() - the same function built from keys held in several vectors, without a concatenated copy, for key sets too large to hold twice.
    /// The mapper's seed index over a human genome has 415M keys in 3.3 GB.
    /// Each vector is freed once its keys are distributed.
    bool build(std::vector<std::vector<uint64_t>>& spans, double c = 6.0, double alpha = 0.94,
               unsigned nthreads = 0) {
        std::size_t n = 0;
        for (const auto& v : spans) n += v.size();
        n_ = n;
        parts_.clear(); part_offset_.clear();
        nparts_ = std::max<std::size_t>(1, (n + kPartKeys - 1) / kPartKeys);
        if (n == 0) { part_offset_.assign(2, 0); parts_.resize(1); return parts_[0].build(nullptr, 0, c, alpha); }
        std::vector<std::size_t> cnt(nparts_, 0);
        for (const auto& v : spans) for (uint64_t k : v) ++cnt[part_of(k)];
        std::vector<std::vector<uint64_t>> part_keys(nparts_);
        for (std::size_t p = 0; p < nparts_; ++p) part_keys[p].reserve(cnt[p]);
        for (auto& v : spans) {
            for (uint64_t k : v) part_keys[part_of(k)].push_back(k);
            std::vector<uint64_t>().swap(v);
        }
        return build_parts(part_keys, c, alpha, nthreads);
    }

private:
    bool build_parts(std::vector<std::vector<uint64_t>>& part_keys, double c, double alpha, unsigned nthreads) {
        parts_.resize(nparts_);
        part_offset_.assign(nparts_ + 1, 0);
        for (std::size_t p = 0; p < nparts_; ++p) part_offset_[p + 1] = part_offset_[p] + part_keys[p].size();
        nthreads = resolve_threads(nthreads);
        nthreads = std::max(1u, std::min<unsigned>(nthreads, static_cast<unsigned>(nparts_)));
        std::atomic<std::size_t> next{0};
        std::atomic<bool> ok{true};
        tpoptoa::run_on_threads(nthreads, [&](unsigned) {
            for (;;) {
                const std::size_t p = next.fetch_add(1);
                if (p >= nparts_) break;
                if (!parts_[p].build(part_keys[p].data(), part_keys[p].size(), c, alpha)) ok = false;
            }
        });
        return ok.load();
    }

public:
    uint64_t operator()(uint64_t key) const noexcept {
        const std::size_t p = part_of(key);
        return part_offset_[p] + parts_[p](key);
    }
    void prefetch(uint64_t key) const noexcept { parts_[part_of(key)].prefetch(key); }
    std::size_t size() const noexcept { return n_; }
    std::size_t bytes() const noexcept {
        std::size_t b = part_offset_.size() * 8 + 32;
        for (const auto& x : parts_) b += x.bytes();
        return b;
    }
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t v) { uint8_t b[8]; std::memcpy(b, &v, 8); out.insert(out.end(), b, b + 8); };
        put(n_); put(parts_.size());   // 0 for a never-built (empty) function
        for (const auto& x : parts_) x.save(out);
    }
    /// load() - read a save() blob, with @borrow as PilotMphf::load() takes it.
    std::size_t load(const uint8_t* p, std::size_t sz, bool borrow = false) {
        std::size_t o = 0;
        auto get = [&](uint64_t& v) { if (o + 8 > sz) return false; std::memcpy(&v, p + o, 8); o += 8; return true; };
        uint64_t n, P;
        if (!get(n) || !get(P) || P > sz) return 0;
        n_ = n; nparts_ = std::max<uint64_t>(P, 1);
        parts_.assign(P, PilotMphf());
        part_offset_.assign(P + 1, 0);
        for (std::size_t i = 0; i < P; ++i) {
            const std::size_t used = parts_[i].load(p + o, sz - o, borrow);
            if (!used) return 0;
            o += used;
            part_offset_[i + 1] = part_offset_[i] + parts_[i].size();
        }
        return part_offset_[P] == n_ ? o : 0;
    }

private:
    static constexpr uint64_t kPartMult1 = 0xD6E8FEB86659FD93ULL, kPartMult2 = 0x9E3779B97F4A7C15ULL;

    /// part_of() - partition of @key, by a hash independent of PilotMphf's.
    std::size_t part_of(uint64_t key) const noexcept {
        uint64_t h = key * kPartMult1;
        h ^= h >> 32;
        return static_cast<std::size_t>((static_cast<__uint128_t>(h * kPartMult2) * nparts_) >> 64);
    }
    std::size_t n_ = 0, nparts_ = 1;
    std::vector<PilotMphf> parts_;
    std::vector<uint64_t> part_offset_;
};

}
