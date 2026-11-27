#pragma once

/**
 * @file
 * @brief Bit-parallel semi-global edit distance.
 *
 * The bit-vector recurrence is Myers 1999, and the multi-word block form is Hyyro 2003.
 * Variable names map onto Myers' notation: match_bits is Peq, pos_vert/neg_vert are Pv/Mv, pos_horiz/neg_horiz are Ph/Mh, and vert_x/horiz_x are Xv/Xh.
 *
 * The query is aligned end to end while the text's ends stay free.
 * In other words, the result is the minimum of ed(query, t) over all text substrings t, at unit costs.
 * Each DP column is held as vertical delta bit-vectors with one bit per query row, so a text base advances a whole word of rows in a few operations.
 * The cost is O(n * ceil(m / 64)) regardless of distance.
 * A wavefront aligner with free text starts, by comparison, has to carry all n start diagonals through every score.
 *
 * @par Free text start
 * Row 0 is 0 in every column, which amounts to a horizontal input of 0 at the top.
 * The score at row m starts at m in column 0 and then moves by the horizontal delta out of row m, and the answer is its minimum.
 * Padding rows past m in the last word are harmless, because carries only travel toward higher rows.
 *
 * Bases A/C/G/T in either case match only themselves.
 * Any other byte matches nothing, not even N against N, which mirrors the wavefront aligner this code replaced for candidate verification.
 *
 * Our addition is bpm_semiglobal_ed_x4(), which runs four texts against one query in the four 64-bit lanes of an AVX2 register.
 */

#include <cstdint>
#include <vector>
#include "bits.hpp"
#include "kmer.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace tpoptoa {

/// Texts bpm_semiglobal_ed_x4() aligns at once, one per 64-bit AVX2 lane.
inline constexpr int kBpmLanes = 4;

namespace bpm_detail {

/// Code for every byte other than A, C, G or T.
/// Its match vector stays all zero, so it never matches.
inline constexpr uint8_t kOther = 4;
inline constexpr int kPeqRows = 5;                 // A, C, G, T, kOther
inline constexpr int kDeltaRows = 2;               // the vertical deltas Pv and Mv, one row each
inline constexpr std::size_t kAvx2Align = 32;      // register width in bytes

/// Maps A/C/G/T in either case to 0..3 and every other byte to kOther.
/// Because it's built at compile time, the table is placed in read-only data and needs no first-call check.
struct BaseCode {
    uint8_t codes[kByteValues] = {};

    constexpr BaseCode() {
        for (int i = 0; i < kByteValues; ++i) codes[i] = kOther;

        codes['A'] = codes['a'] = kCodeA; codes['C'] = codes['c'] = kCodeC;
        codes['G'] = codes['g'] = kCodeG; codes['T'] = codes['t'] = kCodeT;
    }
};
inline constexpr BaseCode kBaseCode{};

inline const BaseCode& base_code() {
    return kBaseCode;
}

}

/// Returns min(semi-global edit distance, @p max_ed + 1).
///
/// When @p prof is given, it receives the DP's last row as @p text_len ints.
/// Entry j is then the query's distance to the best text substring ending at position j (inclusive).
inline int bpm_semiglobal_ed(const char* query, int query_len, const char* text, int text_len, int max_ed,
                             int* prof = nullptr) {
    using namespace bpm_detail;

    if (query_len <= 0) return 0;

    const BaseCode& code = base_code();
    const int num_words = static_cast<int>(words_for_bits(static_cast<uint64_t>(query_len)));
    thread_local std::vector<uint64_t> buf;

    buf.assign(static_cast<std::size_t>(num_words) * (kPeqRows + kDeltaRows), 0);
    uint64_t* peq = buf.data();                 // kPeqRows x W, with the kOther row left at 0
    uint64_t* pos_vert = peq + kPeqRows * num_words;
    uint64_t* neg_vert = pos_vert + num_words;

    for (int i = 0; i < query_len; ++i) {
        const uint8_t base = code.codes[static_cast<unsigned char>(query[i])];

        if (base < kOther) peq[base * num_words + word_index(i)] |= uint64_t(1) << bit_offset(i);
    }

    for (int w = 0; w < num_words; ++w) pos_vert[w] = ~uint64_t(0);

    const uint64_t last = uint64_t(1) << bit_offset(query_len - 1);
    int score = query_len, best = query_len;

    for (int j = 0; j < text_len; ++j) {
        const uint64_t* match_row = peq + code.codes[static_cast<unsigned char>(text[j])] * num_words;
        // Myers' recurrence is applied 64 rows at a time, and hin passes the horizontal delta from one word to the next.
        int hin = 0;                            // free text start, so row 0 is 0

        for (int w = 0; w < num_words; ++w) {
            uint64_t match_bits = match_row[w];
            const uint64_t pos_vert_word = pos_vert[w], neg_vert_word = neg_vert[w];
            const uint64_t hneg = hin < 0 ? 1u : 0u;
            const uint64_t vert_x = match_bits | neg_vert_word;

            match_bits |= hneg;
            // The addition propagates match runs through the vertical deltas in one operation.
            const uint64_t horiz_x = (((match_bits & pos_vert_word) + pos_vert_word) ^ pos_vert_word) | match_bits;
            uint64_t pos_horiz = neg_vert_word | ~(horiz_x | pos_vert_word);
            uint64_t neg_horiz = pos_vert_word & horiz_x;
            // The last word's delta leaves at the query's final row, not at bit 63.
            const uint64_t out_bit = w == num_words - 1 ? last : kTopBit;
            const int hout = (pos_horiz & out_bit) ? 1 : (neg_horiz & out_bit) ? -1 : 0;

            pos_horiz = (pos_horiz << 1) | (hin > 0 ? 1u : 0u);
            neg_horiz = (neg_horiz << 1) | hneg;
            pos_vert[w] = neg_horiz | ~(vert_x | pos_horiz);
            neg_vert[w] = pos_horiz & vert_x;
            hin = hout;
        }

        // The bottom-row score moves by the delta that left the last word.
        score += hin;

        if (score < best) best = score;

        if (prof) prof[j] = score;
    }

    return best > max_ed ? max_ed + 1 : best;
}

/// The best placement of a query in one text, plus the best distinct alternative.
///
/// bpm_placements() reads both from the last-row profile @p prof.
/// First comes the best edit distance with its leftmost end.
/// The alternative is the best distance among ends more than @p sep away, such as one tandem-repeat period over.
/// Distances are capped at @p max_ed + 1, and ends are exclusive text offsets (end = j + 1).
struct BpmPlacements { int best, end, alt, alt_end; };
inline BpmPlacements bpm_placements(const int* prof, int text_len, int sep, int max_ed) {
    BpmPlacements placements{max_ed + 1, 0, max_ed + 1, 0};

    for (int j = 0; j < text_len; ++j)
        if (prof[j] < placements.best) { placements.best = prof[j]; placements.end = j + 1; }

    if (placements.best > max_ed) { placements.best = max_ed + 1; return placements; }

    for (int j = 0; j < text_len; ++j) {
        const int offset = j + 1 - placements.end;

        if ((offset > sep || offset < -sep) && prof[j] < placements.alt) { placements.alt = prof[j]; placements.alt_end = j + 1; }
    }

    if (placements.alt > max_ed) placements.alt = max_ed + 1;

    return placements;
}

/// Gives out[i] = bpm_semiglobal_ed(q, m, t[i], n[i], max_ed) for up to kBpmLanes texts at once.
///
/// With AVX2, each text runs in its own 64-bit lane of the same instruction stream while sharing the query's match vectors.
/// Lanes stay independent because words within a lane are linked only through the +-1 horizontal delta, never an add carry.
/// Once a lane passes the end of its text, it stops updating its minimum.
/// Optionally, prof[i] receives text i's last-row profile of n[i] ints.
/// Without AVX2, the texts run one after another.
inline void bpm_semiglobal_ed_x4(const char* query, int query_len, const char* const* texts, const int* text_lens,
                                 int cnt, int max_ed, int* out, int* const* prof = nullptr) {
#ifdef __AVX2__
    using namespace bpm_detail;
    constexpr int kLanes = kBpmLanes;

    if (query_len <= 0 || cnt <= 1) {
        for (int i = 0; i < cnt; ++i)
            out[i] = bpm_semiglobal_ed(query, query_len, texts[i], text_lens[i], max_ed, prof ? prof[i] : nullptr);

        return;
    }

    const BaseCode& code = base_code();
    const int num_words = static_cast<int>(words_for_bits(static_cast<uint64_t>(query_len)));
    thread_local std::vector<uint64_t> peq_buf;
    thread_local std::vector<uint64_t> delta_buf;   // Pv then Mv, W vectors of kLanes lanes each

    peq_buf.assign(static_cast<std::size_t>(num_words) * kPeqRows, 0);
    delta_buf.resize(static_cast<std::size_t>(num_words) * kDeltaRows * kLanes);
    uint64_t* peq = peq_buf.data();             // kPeqRows x W, with the kOther row left at 0
    auto vec = [&](int word_idx) { return reinterpret_cast<__m256i*>(delta_buf.data() + kLanes * word_idx); };

    for (int i = 0; i < query_len; ++i) {
        const uint8_t base = code.codes[static_cast<unsigned char>(query[i])];

        if (base < kOther) peq[base * num_words + word_index(i)] |= uint64_t(1) << bit_offset(i);
    }

    const __m256i ones = _mm256_set1_epi64x(-1), one = _mm256_set1_epi64x(1);
    const __m256i hbit = _mm256_set1_epi64x(static_cast<long long>(kTopBit));
    const __m256i lastbit = _mm256_set1_epi64x(static_cast<long long>(uint64_t(1) << bit_offset(query_len - 1)));

    for (int w = 0; w < num_words; ++w) {
        _mm256_storeu_si256(vec(w), ones);
        _mm256_storeu_si256(vec(num_words + w), _mm256_setzero_si256());
    }

    const char* lane_text[kLanes];
    int len[kLanes], nmax = 0;

    for (int i = 0; i < kLanes; ++i) {
        lane_text[i] = i < cnt ? texts[i] : texts[0];
        len[i] = i < cnt ? text_lens[i] : 0;
        nmax = len[i] > nmax ? len[i] : nmax;
    }

    // _mm256_set_epi64x takes its four lanes from the highest down, so here and in the gather below the lane indices are spelled out to keep the values in registers.
    static_assert(kLanes == 4, "the lane sets in this function name each of the four lanes");
    const __m256i lenv = _mm256_set_epi64x(len[3], len[2], len[1], len[0]);
    __m256i score = _mm256_set1_epi64x(query_len), best = score;

    for (int j = 0; j < nmax; ++j) {
        // Past a lane's end the code is kOther, so the padding matches nothing.
        unsigned lane_code[kLanes];

        for (int i = 0; i < kLanes; ++i)
            lane_code[i] = j < len[i] ? code.codes[static_cast<unsigned char>(lane_text[i][j])] : unsigned(kOther);

        __m256i hpos = _mm256_setzero_si256(), hneg = _mm256_setzero_si256();

        for (int w = 0; w < num_words; ++w) {
            __m256i match_bits = _mm256_set_epi64x(
                static_cast<long long>(peq[lane_code[3] * num_words + w]), static_cast<long long>(peq[lane_code[2] * num_words + w]),
                static_cast<long long>(peq[lane_code[1] * num_words + w]), static_cast<long long>(peq[lane_code[0] * num_words + w]));
            const __m256i pos_vert_word = _mm256_loadu_si256(vec(w)), neg_vert_word = _mm256_loadu_si256(vec(num_words + w));
            const __m256i hneg1 = _mm256_and_si256(hneg, one);
            const __m256i vert_x = _mm256_or_si256(match_bits, neg_vert_word);

            match_bits = _mm256_or_si256(match_bits, hneg1);
            const __m256i horiz_x = _mm256_or_si256(
                _mm256_xor_si256(_mm256_add_epi64(_mm256_and_si256(match_bits, pos_vert_word), pos_vert_word), pos_vert_word), match_bits);
            __m256i pos_horiz = _mm256_or_si256(neg_vert_word, _mm256_xor_si256(_mm256_or_si256(horiz_x, pos_vert_word), ones));
            __m256i neg_horiz = _mm256_and_si256(pos_vert_word, horiz_x);
            const __m256i out_bit = w == num_words - 1 ? lastbit : hbit;
            const __m256i out_pos = _mm256_cmpeq_epi64(_mm256_and_si256(pos_horiz, out_bit), out_bit);
            const __m256i out_neg = _mm256_cmpeq_epi64(_mm256_and_si256(neg_horiz, out_bit), out_bit);

            pos_horiz = _mm256_or_si256(_mm256_slli_epi64(pos_horiz, 1), _mm256_and_si256(hpos, one));
            neg_horiz = _mm256_or_si256(_mm256_slli_epi64(neg_horiz, 1), hneg1);
            _mm256_storeu_si256(vec(w), _mm256_or_si256(neg_horiz, _mm256_xor_si256(_mm256_or_si256(vert_x, pos_horiz), ones)));
            _mm256_storeu_si256(vec(num_words + w), _mm256_and_si256(pos_horiz, vert_x));
            hpos = out_pos; hneg = out_neg;
        }

        // The score moves by hout, which is +1 where hpos = -1 and -1 where hneg = -1 (the masks are exclusive).
        score = _mm256_add_epi64(_mm256_sub_epi64(score, hpos), hneg);
        // Lanes past their text's end no longer update best.
        const __m256i active = _mm256_cmpgt_epi64(lenv, _mm256_set1_epi64x(j));
        const __m256i lower = _mm256_and_si256(_mm256_cmpgt_epi64(best, score), active);

        best = _mm256_blendv_epi8(best, score, lower);

        if (prof) {
            alignas(kAvx2Align) long long lane_score[kLanes];
            _mm256_store_si256(reinterpret_cast<__m256i*>(lane_score), score);

            for (int i = 0; i < cnt; ++i)
                if (j < len[i]) prof[i][j] = static_cast<int>(lane_score[i]);
        }
    }

    alignas(kAvx2Align) long long lane_best[kLanes];
    _mm256_store_si256(reinterpret_cast<__m256i*>(lane_best), best);

    for (int i = 0; i < cnt; ++i) out[i] = lane_best[i] > max_ed ? max_ed + 1 : static_cast<int>(lane_best[i]);

#else
    for (int i = 0; i < cnt; ++i)
        out[i] = bpm_semiglobal_ed(query, query_len, texts[i], text_lens[i], max_ed, prof ? prof[i] : nullptr);
#endif
}

}
