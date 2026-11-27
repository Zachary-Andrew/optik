#pragma once

/**
 * @file
 * @brief SpssCore serialisation.
 *
 * save() writes the little-endian blob that ffi_spss.cpp wraps with its magic and checksum.
 * load() reads it back, either copying the arrays or borrowing them from a mapped file.
 * spss_core.hpp includes this file at the end, so callers never include it directly.
 */

#include "spss_core.hpp"

namespace tpoptoa {

inline void SpssCore::save(std::vector<uint8_t>& out) const {
    auto put = [&](uint64_t value) { uint8_t bytes[sizeof(uint64_t)]; std::memcpy(bytes, &value, sizeof(uint64_t)); out.insert(out.end(), bytes, bytes + sizeof(uint64_t)); };
    auto putv = [&](const uint64_t* words, std::size_t count) {
        put(count);
        const std::size_t start = out.size();

        out.resize(start + count * sizeof(uint64_t));

        if (count) std::memcpy(out.data() + start, words, count * sizeof(uint64_t));
    };

    put(kmer_len_); put(minimizer_len_); put(lo_bits_); put(rec_bits_); put(n_starts_);
    put(stats_.kmers); put(stats_.bases); put(stats_.seqs); put(stats_.superkmers);
    put(stats_.buckets); put(stats_.ovf_buckets);
    putv(str_view_, str_words_); putv(rec_view_, rec_words_); putv(ovf_view_, ovf_words_);
    // The skew tables' layout parameters, then the tables, 8-byte aligned like every array above.
    put(kSkewMin); put(kFpBits); put(kMaxClass);

    for (unsigned c = 0; c <= kMaxClass; ++c) {
        put(skew_cap_[c]);
        putv(skew_tables_[c], skew_words_[c]);
    }

    mphf_.save(out);
    ef_starts_.save(out);
}

inline std::size_t SpssCore::load(const uint8_t* buf, std::size_t buf_len, bool borrow) {
    has_distinct_ = false;
    std::size_t offset = 0;
    auto get = [&](uint64_t& value) { if (offset + sizeof(uint64_t) > buf_len) return false; std::memcpy(&value, buf + offset, sizeof(uint64_t)); offset += sizeof(uint64_t); return true; };
    // With `borrow`, an array points straight into the buffer (a memory-mapped file) instead of being copied, so the caller must keep the buffer alive.
    auto getv = [&](std::vector<uint64_t>& owned, const uint64_t*& view, std::size_t& count) {
        // A count larger than the bytes left marks a truncated or corrupt file.
        uint64_t word_count; if (!get(word_count) || word_count > (buf_len - offset) / sizeof(uint64_t)) return false;

        if (borrow) {
            owned.clear(); owned.shrink_to_fit();
            view = reinterpret_cast<const uint64_t*>(buf + offset);
        } else {
            owned.resize(word_count); std::memcpy(owned.data(), buf + offset, word_count * sizeof(uint64_t)); view = owned.data();
        }

        count = word_count; offset += word_count * sizeof(uint64_t); return true; };
    uint64_t kmer_len, minimizer_len, low_bits, record_bits, num_starts;

    // load() reads fields in the order save() wrote them, and any short read makes the whole load fail.
    if (!get(kmer_len) || !get(minimizer_len) || !get(low_bits) || !get(record_bits) || !get(num_starts)) return 0;

    if (!get(stats_.kmers) || !get(stats_.bases) || !get(stats_.seqs) || !get(stats_.superkmers)
        || !get(stats_.buckets) || !get(stats_.ovf_buckets)) return 0;

    if (!getv(str_, str_view_, str_words_) || !getv(rec_, rec_view_, rec_words_)
        || !getv(ovf_, ovf_view_, ovf_words_)) return 0;

    skew_own_.assign(kMaxClass + 1, {});
    skew_cap_.assign(kMaxClass + 1, 0);
    skew_tables_.assign(kMaxClass + 1, nullptr);
    skew_words_.assign(kMaxClass + 1, 0);
    uint64_t smin, fpb, maxc;

    if (!get(smin) || !get(fpb) || !get(maxc)) return 0;

    // A file written with other skew-table parameters cannot be read by this layout, so it is rejected as invalid.
    if (smin != kSkewMin || fpb != kFpBits || maxc != kMaxClass) return 0;

    for (unsigned c = 0; c <= kMaxClass; ++c) {
        if (!get(skew_cap_[c])) return 0;

        if (!getv(skew_own_[c], skew_tables_[c], skew_words_[c])) return 0;

        // A table must hold its capacity's slots, or the file is corrupt.
        if (skew_cap_[c] && skew_words_[c] * kWordBits < skew_cap_[c] * (c + kFpBits)) return 0;
    }

    const std::size_t used = mphf_.load(buf + offset, buf_len - offset, borrow);

    if (!used) return 0;

    offset += used;
    const std::size_t ef_used = ef_starts_.load(buf + offset, buf_len - offset, borrow);

    if (!ef_used || ef_starts_.size() != num_starts || !ef_starts_.sane()) return 0;

    offset += ef_used;
    kmer_len_ = unsigned(kmer_len); minimizer_len_ = unsigned(minimizer_len); window_ = kmer_len_ - minimizer_len_ + 1;
    lo_bits_ = unsigned(low_bits); rec_bits_ = unsigned(record_bits);
    n_starts_ = num_starts;
    kmask_ = kmer_mask(kmer_len_); mmask_ = kmer_mask(minimizer_len_);
    return offset;
}

}
