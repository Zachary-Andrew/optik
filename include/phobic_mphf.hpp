#pragma once

/**
 * @file
 * @brief Minimal perfect hash functions over 64-bit keys, and the Elias-Fano sequence they use.
 *
 * The design follows PTHash (Pibiri and Trani 2021) and PHOBIC (Hermann et al. 2024).
 * Keys are split into many small buckets, and each bucket gets a "pilot" p chosen so that all of its keys land on free slots of a table of size T = n / alpha:
 *
 *     pos(key) = fastrange(mix(h(key) ^ mix(p)), T)
 *
 * A lookup is thus one hash, one pilot read and one position computation.
 * Because alpha < 1, a few slots at or past n get used, and those are remapped into the holes below n.
 * fastrange replaces modulo throughout (Lemire 2019).
 *
 * What we changed:
 *   - Skewed bucket assignment, as in PTHash, sends kDenseKeyShare of the keys to the first kDenseBucketShare of the buckets.
 *   - Dense and sparse buckets have very different pilot ranges, so each part stores its pilots at its own fixed width, read with one unaligned load.
 *   - The remap table is non-decreasing, so EfSeq stores it as Elias-Fano.
 *     That makes a low alpha cheap, which in turn means fewer pilot trials and smaller pilots.
 *   - Keys are assumed to be well mixed already, as hashed minimizers are.
 *     They are still re-hashed with a per-build seed, so that a failed build can reseed.
 *   - Mphf partitions the keys and builds one PilotMphf per partition on its own thread, like PHOBIC's partitioned construction.
 *
 * A pilot only locates keys together with their hash, as the data access method we use was originally inspired by the tiny pointers of Bender et al. 2025.
 *
 * A key outside the build set maps to an arbitrary value in [0, n).
 * Callers must therefore verify membership themselves, as the SPSS core does by comparing against its string.
 */

#include "bits.hpp"
#include "thread_budget.hpp"
#include "wyhash.hpp"
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

/// Read-only view of 64-bit words owned elsewhere.
/// After a build they belong to a vector of the owning object, and after a borrowing load to a mapped index file, which avoids a copy and a second resident set.
struct WordView {
    const uint64_t* ptr = nullptr;
    std::size_t count = 0;

    const uint64_t& operator[](std::size_t idx) const noexcept { return ptr[idx]; }
    std::size_t size() const noexcept { return count; }
    bool empty() const noexcept { return count == 0; }
    const uint64_t* data() const noexcept { return ptr; }
    const uint64_t* begin() const noexcept { return ptr; }
    const uint64_t* end() const noexcept { return ptr + count; }
    static WordView from_vector(const std::vector<uint64_t>& words) noexcept { return {words.data(), words.size()}; }
};

/// Reads a u64 count and then that many words.
/// With @p borrow the view points into @p buf, which must then outlive it and be 8-byte aligned, otherwise the words are copied into @p own.
/// Running past @p buf_len returns false.
inline bool load_words(const uint8_t* buf, std::size_t buf_len, std::size_t& offset, bool borrow,
                       std::vector<uint64_t>& own, WordView& view) {
    uint64_t count;

    if (offset + sizeof(uint64_t) > buf_len) return false;

    std::memcpy(&count, buf + offset, sizeof(uint64_t));
    offset += sizeof(uint64_t);

    if (count > (buf_len - offset) / sizeof(uint64_t)) return false;

    if (borrow) {
        own.clear();
        own.shrink_to_fit();
        view = {reinterpret_cast<const uint64_t*>(buf + offset), static_cast<std::size_t>(count)};
    } else {
        own.resize(count);

        if (count) std::memcpy(own.data(), buf + offset, count * sizeof(uint64_t));

        view = WordView::from_vector(own);
    }

    offset += count * sizeof(uint64_t);
    return true;
}

/// Elias-Fano sequence of non-decreasing integers with O(1) access (Elias 1974, Fano 1971).
/// Select sampling follows Vigna, "Quasi-succinct indices" (WSDM 2013): every kSampleRate-th set bit of the high bitvector has its position stored.
/// access(i) jumps to the nearest sample, popcount-scans whole words and finishes with select_in_word().
class EfSeq {
public:
    static constexpr uint64_t kSampleRate = 256;
    static constexpr std::size_t kPadWords = 2;      // low() reads two words at a time
    static constexpr uint64_t kHighBitsSlack = 2;    // spare high bits so select never runs off the end
    static constexpr std::size_t kScalarBytes = 16;  // rough size of the members, for bytes()

    void build(const std::vector<uint64_t>& values, uint64_t universe) {
        num_values_ = values.size();
        nlow_ = 0;

        // Elias-Fano keeps floor(log2(universe / n)) low bits per value, which minimises total space.
        if (num_values_ && universe / num_values_ > 1) while ((uint64_t(1) << (nlow_ + 1)) <= universe / num_values_) ++nlow_;

        low_own_.assign(words_for_bits(num_values_ * nlow_) + kPadWords, 0);
        const uint64_t hbits = num_values_ + (universe >> nlow_) + kHighBitsSlack;

        high_own_.assign(hbits / kWordBits + kPadWords, 0);
        one_own_.clear();

        for (uint64_t i = 0; i < num_values_; ++i) {
            if (nlow_) {
                const uint64_t lowv = values[i] & ((uint64_t(1) << nlow_) - 1), bit = i * nlow_;

                low_own_[word_index(bit)] |= lowv << bit_offset(bit);

                if (bit_offset(bit) + nlow_ > kWordBits) low_own_[word_index(bit) + 1] |= lowv >> (kWordBits - bit_offset(bit));
            }

            // Value i's high part, offset by i, gives a distinct set bit in the unary high array.
            const uint64_t high_pos = (values[i] >> nlow_) + i;

            high_own_[word_index(high_pos)] |= uint64_t(1) << bit_offset(high_pos);

            // Every kSampleRate-th one is sampled to speed up select1.
            if (i % kSampleRate == 0) one_own_.push_back(high_pos);
        }

        borrowed_ = false;
        adopt_owned();
        build_zero_samples();
    }

    /// Counts the values <= @p value.
    ///
    /// The bucket of @p value's high part starts after the (@p value >> l)-th zero of the high bits, which sampled select0 finds.
    /// Its few entries are then compared on their low bits, in O(1) expected time rather than a binary search over operator[].
    uint64_t count_le(uint64_t value) const noexcept {
        if (num_values_ == 0) return 0;

        const uint64_t high = value >> nlow_;
        uint64_t pos;                                   // where bucket h starts in the high bits

        if (high == 0) pos = 0;
        else {
            if (high - 1 >= nzeros_) return num_values_;            // past every bucket

            pos = select0(high - 1) + 1;
        }

        uint64_t idx = pos - high;                         // how many values have a high part below h
        const uint64_t lowx = value & ((uint64_t(1) << nlow_) - 1);

        while (idx < num_values_ && ((high_bits_[word_index(pos)] >> bit_offset(pos)) & 1)) {
            if (low(idx) > lowx) break;

            ++idx; ++pos;
        }

        return idx;
    }

    uint64_t operator[](uint64_t idx) const noexcept {
        return ((select1(idx) - idx) << nlow_) | low(idx);
    }

    std::size_t bytes() const noexcept {
        return low_parts_.size() * sizeof(uint64_t) + high_bits_.size() * sizeof(uint64_t) + one_samples_.size() * sizeof(uint64_t) + kScalarBytes;
    }
    uint64_t size() const noexcept { return num_values_; }
    std::size_t mem_bytes() const noexcept { return bytes() + zero_samples_.size() * sizeof(uint64_t); }

    /// Checks after a load that the arrays hold num_values_ values of nlow_ bits and that every sample points into the high bits.
    /// A truncated or corrupted blob is rejected here instead of being read out of bounds later.
    bool sane() const noexcept {
        if (nlow_ > kWordBits - 1 || low_parts_.size() < words_for_bits(num_values_ * nlow_) + 1 || one_samples_.size() < (num_values_ + kSampleRate - 1) / kSampleRate) return false;

        if (num_values_ && high_bits_.empty()) return false;

        for (uint64_t sample_pos : one_samples_) if (word_index(sample_pos) >= high_bits_.size()) return false;

        uint64_t ones = 0;

        for (uint64_t word : high_bits_) ones += static_cast<uint64_t>(__builtin_popcountll(word));

        return ones == num_values_;
    }
    /// Little-endian serialization, in which the select0 samples are rebuilt on load rather than stored.
    /// load() returns the bytes consumed, or 0 for a malformed blob.
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t value) { uint8_t bytes[sizeof(uint64_t)]; std::memcpy(bytes, &value, sizeof(uint64_t)); out.insert(out.end(), bytes, bytes + sizeof(uint64_t)); };

        put(num_values_); put(nlow_);

        for (const auto* vec : {&low_parts_, &high_bits_, &one_samples_}) { put(vec->size()); for (uint64_t word : *vec) put(word); }
    }
    /// With @p borrow the three arrays are used in place, so @p buf must stay mapped and 8-byte aligned for this object's lifetime.
    std::size_t load(const uint8_t* buf, std::size_t buf_len, bool borrow = false) {
        std::size_t offset = 0;
        auto get = [&](uint64_t& value) { if (offset + sizeof(uint64_t) > buf_len) return false; std::memcpy(&value, buf + offset, sizeof(uint64_t)); offset += sizeof(uint64_t); return true; };

        if (!get(num_values_) || !get(nlow_)) return 0;

        if (!load_words(buf, buf_len, offset, borrow, low_own_, low_parts_) || !load_words(buf, buf_len, offset, borrow, high_own_, high_bits_)
            || !load_words(buf, buf_len, offset, borrow, one_own_, one_samples_)) return 0;

        borrowed_ = borrow;
        build_zero_samples();
        return offset;
    }

    EfSeq() = default;
    EfSeq(EfSeq&&) noexcept = default;               // moved vectors keep their buffers
    EfSeq& operator=(EfSeq&&) noexcept = default;
    /// Copies take their own arrays and repoint the views, unless the source borrows a mapping, in which case the views stay valid as they are.
    EfSeq(const EfSeq& other) { *this = other; }
    EfSeq& operator=(const EfSeq& other) {
        if (this == &other) return *this;

        num_values_ = other.num_values_; nlow_ = other.nlow_; nzeros_ = other.nzeros_; borrowed_ = other.borrowed_;
        low_own_ = other.low_own_; high_own_ = other.high_own_; one_own_ = other.one_own_;
        zero_samples_ = other.zero_samples_;

        if (borrowed_) { low_parts_ = other.low_parts_; high_bits_ = other.high_bits_; one_samples_ = other.one_samples_; }
        else adopt_owned();

        return *this;
    }

private:
    uint64_t low(uint64_t idx) const noexcept {
        if (!nlow_) return 0;

        const uint64_t bit = idx * nlow_;
        uint64_t word;

        std::memcpy(&word, reinterpret_cast<const uint8_t*>(low_parts_.data()) + bits_to_bytes(bit), sizeof(uint64_t));
        return (word >> (bit % kBitsPerByte)) & ((uint64_t(1) << nlow_) - 1);
    }
    /// zero_samples_[t] holds the position of the (t * kSampleRate)-th zero of high_bits_.
    /// Only positions in [0, n + (max value >> l) + 1) can hold buckets, so only they count.
    void build_zero_samples() {
        zero_samples_.clear();
        nzeros_ = 0;

        for (std::size_t wi = 0; wi < high_bits_.size(); ++wi) {
            const uint64_t zero_bits = ~high_bits_[wi];
            const uint64_t zeros_in_word = static_cast<uint64_t>(__builtin_popcountll(zero_bits));

            // Find the next sampled rank (a multiple of kSampleRate) that falls in this word.
            for (uint64_t sample_rank = (nzeros_ + kSampleRate - 1) / kSampleRate * kSampleRate; sample_rank < nzeros_ + zeros_in_word; sample_rank += kSampleRate)
                zero_samples_.push_back(first_bit_of_word(uint64_t(wi))
                                    + select_in_word(zero_bits, static_cast<unsigned>(sample_rank - nzeros_)));

            nzeros_ += zeros_in_word;
        }
    }
    uint64_t select0(uint64_t rank) const noexcept {        // finds the rank-th zero
        const uint64_t sample_pos = zero_samples_[rank / kSampleRate];
        uint64_t rem = rank % kSampleRate;
        uint64_t word_idx = word_index(sample_pos);
        uint64_t word = ~high_bits_[word_idx] & (~uint64_t(0) << bit_offset(sample_pos));

        for (;;) {
            const uint64_t set_bits = static_cast<uint64_t>(__builtin_popcountll(word));

            if (rem < set_bits) return first_bit_of_word(word_idx) + select_in_word(word, static_cast<unsigned>(rem));

            rem -= set_bits;
            word = ~high_bits_[++word_idx];
        }
    }
    uint64_t select1(uint64_t idx) const noexcept {
        const uint64_t sample_pos = one_samples_[idx / kSampleRate];
        uint64_t remaining = idx % kSampleRate;                          // ones left to skip past the sample
        uint64_t word_idx = word_index(sample_pos);
        uint64_t word = high_bits_[word_idx] & (~uint64_t(0) << bit_offset(sample_pos));

        for (;;) {
            const uint64_t set_bits = static_cast<uint64_t>(__builtin_popcountll(word));

            if (remaining < set_bits) return first_bit_of_word(word_idx) + select_in_word(word, static_cast<unsigned>(remaining));

            remaining -= set_bits;
            word = high_bits_[++word_idx];
        }
    }
    uint64_t num_values_ = 0, nlow_ = 0;
    // Every lookup reads these views, which point at the owned vectors after a build or at a mapped file after a borrowing load.
    std::vector<uint64_t> low_own_, high_own_, one_own_;
    WordView low_parts_, high_bits_, one_samples_;
    bool borrowed_ = false;
    std::vector<uint64_t> zero_samples_;                    // always derived, never stored

    void adopt_owned() noexcept {
        low_parts_ = WordView::from_vector(low_own_);
        high_bits_ = WordView::from_vector(high_own_);
        one_samples_ = WordView::from_vector(one_own_);
    }
    uint64_t nzeros_ = 0;
};

/// Minimal perfect hash function mapping n keys onto [0, n).
class PilotMphf {
public:
    /// Bucket factor c (giving c n / log2 n buckets) and table load alpha before remapping, as build() explains.
    static constexpr double kDefaultBucketFactor = 6.0;
    static constexpr double kDefaultAlpha = 0.85;

    PilotMphf() = default;

    /// Fails on duplicate keys, or when all kMaxBuildAttempts seeds fail.
    ///
    /// @p bucket_factor (c) sets the bucket count to c * n / log2(n), and @p alpha is the table load before remapping.
    /// When chosen, c = 6 and alpha = 0.85 gave the smallest total size with split pilot widths and the Elias-Fano remap.
    /// On 0.54M random keys that was 4.7 bits per key, against 6.9 for the earlier c = 8, alpha = 0.97 with a 32-bit remap, at the same build time.
    /// Mphf, however, passes its own kDefaultAlpha.
    bool build(const uint64_t* keys, std::size_t count, double bucket_factor = kDefaultBucketFactor,
               double alpha = kDefaultAlpha) {
        num_keys_ = count;

        dense_own_.clear(); sparse_own_.clear();
        borrowed_ = false;
        adopt_owned();

        if (count == 0) { table_size_ = 0; nbuckets_ = 0; dense_width_ = sparse_width_ = 0; return true; }

        for (uint64_t attempt = 0; attempt < kMaxBuildAttempts; ++attempt) {
            seed_ = kSeedStep * (attempt + 1) ^ kSeedBase;

            if (try_build(keys, count, bucket_factor, alpha)) return true;
        }

        return false;
    }

    /// Only meaningful for keys in the build set.
    uint64_t operator()(uint64_t key) const noexcept {
        const uint64_t hash = mix(key ^ seed_);
        const uint64_t pilot_val = pilot(bucket_of(hash));
        const uint64_t pos = fastrange(mix(hash ^ mix(pilot_val + kPilotSalt)), table_size_);

        return pos < num_keys_ ? pos : remap_[pos - num_keys_];
    }

    /// Prefetches the pilot that a later operator()(@p key) will read.
    void prefetch(uint64_t key) const noexcept {
        const uint64_t hash = mix(key ^ seed_);
        const uint64_t bucket = bucket_of(hash);

        if (bucket < ndense_) prefetch_read_low(dense_pilots_.data() + word_index(bucket * dense_width_));
        else prefetch_read_low(sparse_pilots_.data() + word_index((bucket - ndense_) * sparse_width_));
    }

    std::size_t size() const noexcept { return num_keys_; }

    std::size_t bytes() const noexcept {
        return dense_pilots_.size() * sizeof(uint64_t) + sparse_pilots_.size() * sizeof(uint64_t) + remap_.bytes() + kScalarBytes;
    }

    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t value) { uint8_t bytes[sizeof(uint64_t)]; std::memcpy(bytes, &value, sizeof(uint64_t)); out.insert(out.end(), bytes, bytes + sizeof(uint64_t)); };

        put(num_keys_); put(table_size_); put(nbuckets_); put(ndense_); put(seed_); put(dense_width_); put(sparse_width_);
        put(dense_pilots_.size()); for (uint64_t word : dense_pilots_) put(word);
        put(sparse_pilots_.size()); for (uint64_t word : sparse_pilots_) put(word);
        remap_.save(out);
    }
    /// Returns 0 for a malformed blob.
    std::size_t load(const uint8_t* buf, std::size_t buf_len, bool borrow = false) {
        std::size_t offset = 0;
        auto get = [&](uint64_t& value) { if (offset + sizeof(uint64_t) > buf_len) return false; std::memcpy(&value, buf + offset, sizeof(uint64_t)); offset += sizeof(uint64_t); return true; };

        if (!get(num_keys_) || !get(table_size_) || !get(nbuckets_) || !get(ndense_) || !get(seed_) || !get(dense_width_) || !get(sparse_width_)) return 0;

        if (!load_words(buf, buf_len, offset, borrow, dense_own_, dense_pilots_)
            || !load_words(buf, buf_len, offset, borrow, sparse_own_, sparse_pilots_)) return 0;

        borrowed_ = borrow;
        const std::size_t used = remap_.load(buf + offset, buf_len - offset, borrow);

        if (!used) return 0;

        return offset + used;
    }

private:
    static constexpr uint64_t kMaxBuildAttempts = 16;   // seeds tried before build() gives up
    static constexpr uint64_t kSeedStep = kGoldenRatio64;
    static constexpr uint64_t kSeedBase = 0xD6E8FEB86659FD93ULL;
    static constexpr uint64_t kPilotSalt = 0x5851F42D4C957F2DULL;   // added to a pilot before mixing
    static constexpr uint64_t kMaxPilot = uint64_t(1) << 24;         // a bucket needing a larger pilot forces a reseed
    static constexpr double kDenseKeyShare = 0.6;      // share of keys in the dense buckets
    static constexpr double kDenseBucketShare = 0.3;   // share of buckets that are dense
    static constexpr std::size_t kPadWords = 2;      // pilot reads may span two words
    static constexpr std::size_t kScalarBytes = 64;  // rough size of the members, for bytes()

    static constexpr uint64_t mix(uint64_t value) noexcept {
        return mix64(value, 0);
    }
    static constexpr uint64_t fastrange(uint64_t hash, uint64_t range) noexcept {
        return static_cast<uint64_t>((static_cast<__uint128_t>(hash) * range) >> kWordBits);
    }
    /// Skewed split that sends kDenseKeyShare of the hash space to the first ndense_ buckets.
    uint64_t bucket_of(uint64_t hash) const noexcept {
        static constexpr uint64_t kDenseHashCutoff = static_cast<uint64_t>(kDenseKeyShare * static_cast<double>(UINT64_MAX));

        return hash < kDenseHashCutoff ? fastrange(hash, ndense_) : ndense_ + fastrange(hash, nbuckets_ - ndense_);
    }
    /// Reads one @p width-bit value with a single unaligned 64-bit load.
    static uint64_t read_packed(const WordView& words, uint64_t idx, uint64_t width) noexcept {
        if (width == 0) return 0;

        const uint64_t bit = idx * width;
        uint64_t bits;

        std::memcpy(&bits, reinterpret_cast<const uint8_t*>(words.data()) + bits_to_bytes(bit), sizeof(uint64_t));
        return (bits >> (bit % kBitsPerByte)) & ((uint64_t(1) << width) - 1);
    }
    static void write_packed(std::vector<uint64_t>& words, uint64_t idx, uint64_t width, uint64_t val) noexcept {
        if (width == 0) return;

        const uint64_t bit = idx * width;

        words[word_index(bit)] |= val << bit_offset(bit);

        if (bit_offset(bit) + width > kWordBits) words[word_index(bit) + 1] |= val >> (kWordBits - bit_offset(bit));
    }
    uint64_t pilot(uint64_t bucket) const noexcept {
        return bucket < ndense_ ? read_packed(dense_pilots_, bucket, dense_width_) : read_packed(sparse_pilots_, bucket - ndense_, sparse_width_);
    }

    /// One attempt with the current seed_.
    /// It fails when a duplicate turns up or a bucket needs a pilot above kMaxPilot.
    bool try_build(const uint64_t* keys, std::size_t count, double bucket_factor, double alpha) {
        table_size_ = std::max<uint64_t>(count, static_cast<uint64_t>(std::ceil(count / alpha)));
        const double log_n = std::max(1.0, std::log2(static_cast<double>(count)));

        nbuckets_ = std::max<uint64_t>(1, static_cast<uint64_t>(std::ceil(bucket_factor * count / log_n)));
        ndense_ = std::max<uint64_t>(1, static_cast<uint64_t>(kDenseBucketShare * nbuckets_));

        if (ndense_ >= nbuckets_) ndense_ = nbuckets_ - 1 ? nbuckets_ - 1 : 1;

        if (nbuckets_ == 1) ndense_ = 1;

        // Group the hashes by bucket with a counting sort.
        std::vector<uint64_t> hashes(count);

        for (std::size_t i = 0; i < count; ++i) hashes[i] = mix(keys[i] ^ seed_);

        std::vector<uint32_t> cnt(nbuckets_ + 1, 0);
        std::vector<uint64_t> bucket_ids(count);

        for (std::size_t i = 0; i < count; ++i) { bucket_ids[i] = bucket_of(hashes[i]); ++cnt[bucket_ids[i] + 1]; }

        for (uint64_t b = 0; b < nbuckets_; ++b) cnt[b + 1] += cnt[b];

        std::vector<uint64_t> grouped(count);

        {
            std::vector<uint32_t> cursor(cnt.begin(), cnt.end() - 1);

            for (std::size_t i = 0; i < count; ++i) grouped[cursor[bucket_ids[i]]++] = hashes[i];
        }
        // Duplicate keys or hashes can only meet within one bucket of about log2(n)/c keys.
        // Checking there is far cheaper than sorting all n hashes.
        uint32_t max_sz = 0;

        for (uint64_t b = 0; b < nbuckets_; ++b) {
            const uint32_t bucket_lo = cnt[b], bucket_hi = cnt[b + 1];

            max_sz = std::max(max_sz, bucket_hi - bucket_lo);

            for (uint32_t x = bucket_lo + 1; x < bucket_hi; ++x)
                for (uint32_t y = bucket_lo; y < x; ++y)
                    if (grouped[x] == grouped[y]) return false;
        }

        // Process the largest buckets first.
        // A stable counting sort on size keeps ties in bucket order, exactly as the former stable_sort did.
        std::vector<uint32_t> order(nbuckets_);

        {
            // Sizes run 0..max_sz, and one more slot holds the running total for the prefix sum.
            constexpr uint32_t kPrefixSumExtraSlots = 2;
            std::vector<uint32_t> start(max_sz + kPrefixSumExtraSlots, 0);

            // A counting sort orders buckets largest first. Big buckets are hardest to place and go while the table is emptiest.
            for (uint64_t b = 0; b < nbuckets_; ++b) ++start[max_sz - (cnt[b + 1] - cnt[b]) + 1];

            for (uint32_t z = 0; z <= max_sz; ++z) start[z + 1] += start[z];

            for (uint64_t b = 0; b < nbuckets_; ++b)
                order[start[max_sz - (cnt[b + 1] - cnt[b])]++] = static_cast<uint32_t>(b);
        }

        std::vector<uint64_t> taken(words_for_bits(table_size_), 0);
        std::vector<uint64_t> pilots(nbuckets_, 0);
        std::vector<uint64_t> pos;
        uint64_t max_pilot = 0;

        for (uint32_t bucket : order) {
            const uint32_t bucket_lo = cnt[bucket], bucket_hi = cnt[bucket + 1];

            if (bucket_lo == bucket_hi) continue;

            pos.resize(bucket_hi - bucket_lo);
            uint64_t pilot_val = 0;

            for (;; ++pilot_val) {
                if (pilot_val > kMaxPilot) return false;

                // Try pilots in turn until every key of the bucket lands on a free, distinct slot.
                const uint64_t pilot_hash = mix(pilot_val + kPilotSalt);
                bool fits = true;

                for (uint32_t j = bucket_lo; j < bucket_hi && fits; ++j) {
                    const uint64_t slot = fastrange(mix(grouped[j] ^ pilot_hash), table_size_);

                    if ((taken[word_index(slot)] >> bit_offset(slot)) & 1) { fits = false; break; }

                    // Two keys of one bucket may not share a slot either.
                    for (uint32_t t = bucket_lo; t < j; ++t) if (pos[t - bucket_lo] == slot) { fits = false; break; }

                    pos[j - bucket_lo] = slot;
                }

                if (fits) break;
            }

            for (uint32_t j = bucket_lo; j < bucket_hi; ++j) taken[word_index(pos[j - bucket_lo])] |= uint64_t(1) << bit_offset(pos[j - bucket_lo]);

            pilots[bucket] = pilot_val;
            max_pilot = std::max(max_pilot, pilot_val);
        }

        (void)max_pilot;
        // Pilots are stored in two fixed-width parts, plus padding for the unaligned load.
        auto width_of = [](uint64_t max_val) { uint64_t width = 0; while ((uint64_t(1) << width) <= max_val) ++width; return width; };
        uint64_t dense_max = 0, sparse_max = 0;

        for (uint64_t b = 0; b < nbuckets_; ++b) (b < ndense_ ? dense_max : sparse_max) = std::max(b < ndense_ ? dense_max : sparse_max, pilots[b]);

        dense_width_ = width_of(dense_max); sparse_width_ = width_of(sparse_max);
        dense_own_.assign(words_for_bits(ndense_ * dense_width_) + kPadWords, 0);
        sparse_own_.assign(words_for_bits((nbuckets_ - ndense_) * sparse_width_) + kPadWords, 0);

        for (uint64_t b = 0; b < nbuckets_; ++b) {
            if (b < ndense_) write_packed(dense_own_, b, dense_width_, pilots[b]);
            else write_packed(sparse_own_, b - ndense_, sparse_width_, pilots[b]);
        }

        borrowed_ = false;
        adopt_owned();
        // Remap the slack slots >= n onto the holes below n.
        // Holes are assigned in increasing order, and untaken slots repeat the previous value because no lookup reads them.
        // The table therefore stays non-decreasing, which is what Elias-Fano needs.
        std::vector<uint64_t> remap_values(table_size_ - num_keys_, 0);
        uint64_t hole = 0, last = 0;
        auto is_taken = [&](uint64_t slot) { return (taken[word_index(slot)] >> bit_offset(slot)) & 1; };

        for (uint64_t q = num_keys_; q < table_size_; ++q) {
            if (is_taken(q)) {
                while (is_taken(hole)) ++hole;

                last = hole++;
            }

            remap_values[q - num_keys_] = last;
        }

        remap_.build(remap_values, num_keys_);
        return true;
    }

    uint64_t num_keys_ = 0, table_size_ = 0, nbuckets_ = 0, ndense_ = 1, seed_ = 0, dense_width_ = 0, sparse_width_ = 0;
    // Pilots of the dense buckets [0, ndense_) and of the rest, each packed at its own width.
    // Like EfSeq's arrays, they are views of the owned vectors after a build or of a mapped file after a borrowing load.
    std::vector<uint64_t> dense_own_, sparse_own_;
    WordView dense_pilots_, sparse_pilots_;
    bool borrowed_ = false;
    EfSeq remap_;

    void adopt_owned() noexcept {
        dense_pilots_ = WordView::from_vector(dense_own_);
        sparse_pilots_ = WordView::from_vector(sparse_own_);
    }

public:
    PilotMphf(PilotMphf&&) noexcept = default;       // moved vectors keep their buffers
    PilotMphf& operator=(PilotMphf&&) noexcept = default;
    /// Copies take their own pilots and repoint the views, unless the source borrows a mapping.
    PilotMphf(const PilotMphf& other) { *this = other; }
    PilotMphf& operator=(const PilotMphf& other) {
        if (this == &other) return *this;

        num_keys_ = other.num_keys_; table_size_ = other.table_size_; nbuckets_ = other.nbuckets_; ndense_ = other.ndense_;
        seed_ = other.seed_; dense_width_ = other.dense_width_; sparse_width_ = other.sparse_width_;
        dense_own_ = other.dense_own_; sparse_own_ = other.sparse_own_; borrowed_ = other.borrowed_;
        remap_ = other.remap_;

        if (borrowed_) { dense_pilots_ = other.dense_pilots_; sparse_pilots_ = other.sparse_pilots_; }
        else adopt_owned();

        return *this;
    }
};

/// Partitioned minimal perfect hash, after the partitioned constructions of PHOBIC and PTHash.
///
/// An independent hash splits the keys into partitions of about kPartKeys keys.
/// Because the pilot search is sequential within a partition, the partitions build their PilotMphfs in parallel.
/// A key's rank is offset[part] + its local rank.
/// The partition table stays in cache, so a lookup still costs one pilot access, and each build's working set stays about L2-sized.
class Mphf {
public:
    static constexpr std::size_t kPartKeys = 250000;
    /// Per-partition settings, with a higher table load than PilotMphf's default.
    static constexpr double kDefaultBucketFactor = 6.0;
    static constexpr double kDefaultAlpha = 0.94;

    /// Fails on duplicate keys.
    bool build(const uint64_t* keys, std::size_t count, double bucket_factor = kDefaultBucketFactor, double alpha = kDefaultAlpha,
               unsigned nthreads = 0) {
        num_keys_ = count;

        parts_.clear(); part_offset_.clear();
        nparts_ = std::max<std::size_t>(1, (count + kPartKeys - 1) / kPartKeys);

        if (count == 0) { part_offset_.assign(kSinglePartOffsets, 0); parts_.resize(1); return parts_[0].build(keys, 0, bucket_factor, alpha); }

        std::vector<std::vector<uint64_t>> part_keys(nparts_);

        for (std::size_t i = 0; i < count; ++i) part_keys[part_of(keys[i])].push_back(keys[i]);

        return build_parts(part_keys, bucket_factor, alpha, nthreads);
    }

    /// Builds the same function from keys spread over several vectors, skipping the concatenated copy.
    /// That matters for key sets too large to hold twice, such as the mapper's 415M human-genome keys (3.3 GB).
    /// Each vector is freed as soon as its keys are distributed.
    bool build(std::vector<std::vector<uint64_t>>& spans, double bucket_factor = kDefaultBucketFactor, double alpha = kDefaultAlpha,
               unsigned nthreads = 0) {
        std::size_t count = 0;

        for (const auto& span : spans) count += span.size();

        num_keys_ = count;
        parts_.clear(); part_offset_.clear();
        nparts_ = std::max<std::size_t>(1, (count + kPartKeys - 1) / kPartKeys);

        if (count == 0) { part_offset_.assign(kSinglePartOffsets, 0); parts_.resize(1); return parts_[0].build(nullptr, 0, bucket_factor, alpha); }

        std::vector<std::size_t> cnt(nparts_, 0);

        // Counting first lets every partition reserve its exact size, so the distribution below never reallocates.
        for (const auto& span : spans) for (uint64_t key : span) ++cnt[part_of(key)];

        std::vector<std::vector<uint64_t>> part_keys(nparts_);

        for (std::size_t p = 0; p < nparts_; ++p) part_keys[p].reserve(cnt[p]);

        for (auto& span : spans) {
            for (uint64_t key : span) part_keys[part_of(key)].push_back(key);

            std::vector<uint64_t>().swap(span);
        }

        return build_parts(part_keys, bucket_factor, alpha, nthreads);
    }

private:
    bool build_parts(std::vector<std::vector<uint64_t>>& part_keys, double bucket_factor, double alpha, unsigned nthreads) {
        parts_.resize(nparts_);
        part_offset_.assign(nparts_ + 1, 0);

        // A key's global value is its partition's offset plus its value inside the partition, so the offsets are prefix sums of the partition sizes.
        for (std::size_t p = 0; p < nparts_; ++p) part_offset_[p + 1] = part_offset_[p] + part_keys[p].size();

        nthreads = resolve_threads(nthreads);
        nthreads = std::max(1u, std::min<unsigned>(nthreads, static_cast<unsigned>(nparts_)));
        std::atomic<std::size_t> next{0};
        std::atomic<bool> all_built{true};

        // Partitions are independent, so threads claim them one at a time from a shared counter.
        tpoptoa::run_on_threads(nthreads, [&](unsigned) {
            for (;;) {
                const std::size_t part = next.fetch_add(1);

                if (part >= nparts_) break;

                if (!parts_[part].build(part_keys[part].data(), part_keys[part].size(), bucket_factor, alpha)) all_built = false;
            }
        });
        return all_built.load();
    }

public:
    uint64_t operator()(uint64_t key) const noexcept {
        const std::size_t part = part_of(key);

        // A partition's offset turns its local rank into a global one.
        return part_offset_[part] + parts_[part](key);
    }
    void prefetch(uint64_t key) const noexcept { parts_[part_of(key)].prefetch(key); }
    std::size_t size() const noexcept { return num_keys_; }

    std::size_t bytes() const noexcept {
        std::size_t total = part_offset_.size() * sizeof(uint64_t) + kScalarBytes;

        for (const auto& part_fn : parts_) total += part_fn.bytes();

        return total;
    }
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t value) { uint8_t bytes[sizeof(uint64_t)]; std::memcpy(bytes, &value, sizeof(uint64_t)); out.insert(out.end(), bytes, bytes + sizeof(uint64_t)); };

        put(num_keys_); put(parts_.size());   // an empty, never-built function stores 0 here

        for (const auto& part_fn : parts_) part_fn.save(out);
    }
    /// @p borrow works as in PilotMphf::load().
    std::size_t load(const uint8_t* buf, std::size_t buf_len, bool borrow = false) {
        std::size_t offset = 0;
        auto get = [&](uint64_t& value) { if (offset + sizeof(uint64_t) > buf_len) return false; std::memcpy(&value, buf + offset, sizeof(uint64_t)); offset += sizeof(uint64_t); return true; };
        uint64_t count, num_parts;

        if (!get(count) || !get(num_parts) || num_parts > buf_len) return 0;

        num_keys_ = count; nparts_ = std::max<uint64_t>(num_parts, 1);
        parts_.assign(num_parts, PilotMphf());
        part_offset_.assign(num_parts + 1, 0);

        for (std::size_t i = 0; i < num_parts; ++i) {
            const std::size_t used = parts_[i].load(buf + offset, buf_len - offset, borrow);

            if (!used) return 0;

            offset += used;
            part_offset_[i + 1] = part_offset_[i] + parts_[i].size();
        }

        return part_offset_[num_parts] == num_keys_ ? offset : 0;
    }

private:
    static constexpr uint64_t kPartMult1 = 0xD6E8FEB86659FD93ULL, kPartMult2 = kGoldenRatio64;
    static constexpr unsigned kPartMixShift = 32;    // folds the hash's high half into its low half
    static constexpr std::size_t kSinglePartOffsets = 2;   // start and end of the only partition
    static constexpr std::size_t kScalarBytes = 32;  // rough size of the members, for bytes()

    /// Uses a hash independent of PilotMphf's, so partitioning doesn't bias the buckets.
    std::size_t part_of(uint64_t key) const noexcept {
        uint64_t hash = key * kPartMult1;

        hash ^= hash >> kPartMixShift;
        return static_cast<std::size_t>((static_cast<__uint128_t>(hash * kPartMult2) * nparts_) >> kWordBits);
    }
    std::size_t num_keys_ = 0, nparts_ = 1;
    std::vector<PilotMphf> parts_;
    std::vector<uint64_t> part_offset_;
};

}
