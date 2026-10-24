#pragma once

/*
 * bpm_edit.hpp - bit-parallel semi-global edit distance.
 *
 * Sources are Myers 1999 for the bit-vector recurrence and Hyyro 2003 for its multi-word block form.
 * Variable names follow Myers' notation, with Peq for the match vectors, Pv/Mv for the vertical +1/-1 deltas, Ph/Mh for the horizontal deltas, and Xv/Xh.
 *
 * The query is aligned end to end and the text's ends are free, so the result is the minimum over text substrings t of ed(query, t) at unit costs.
 * Column j of the DP matrix is kept as vertical delta bit-vectors, one bit per query row and 64 rows per word, and a text base advances all words with a handful of word operations.
 * The cost is O(n * ceil(m / 64)) whatever the distance, where a wavefront aligner with free text starts must carry all n start diagonals through each score.
 *
 * Free text start.
 * Row 0 is 0 in every column, i.e. a horizontal input of 0 at the top.
 * The score at row m starts at m in column 0 and moves by the horizontal delta out of row m in each column, and the answer is its minimum.
 * Rows past m in the last word never influence rows <= m (carries in the addition only move toward higher rows), so the padding is inert.
 *
 * Bases A/C/G/T, in either case, match only themselves, and any other byte in the query or the text matches nothing, not even N against N, as in the wavefront aligner this replaced for candidate verification.
 *
 * As a modification, bpm_semiglobal_ed_x4() runs four texts against one query in the four 64-bit lanes of an AVX2 register.
 */

#include <cstdint>
#include <vector>

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace tpoptoa {

namespace bpm_detail {

/// Code of every byte that is not A, C, G or T.
/// Its match vector is all zero, so it matches nothing.
inline constexpr uint8_t kOther = 4;
/// Match-vector rows for A, C, G, T and kOther.
inline constexpr int kPeqRows = 5;

/// BaseCode - byte -> 0..3 for A/C/G/T in either case, and kOther for the rest.
/// Built at compile time, so the table sits in read-only data with no first-call initialization check.
struct BaseCode {
    uint8_t c[256] = {};
    constexpr BaseCode() {
        for (int i = 0; i < 256; ++i) c[i] = kOther;
        c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
        c['G'] = c['g'] = 2; c['T'] = c['t'] = 3;
    }
};
inline constexpr BaseCode kBaseCode{};

inline const BaseCode& base_code() {
    return kBaseCode;
}

}

/// bpm_semiglobal_ed() - semi-global edit distance of the @m-base query @q in the @n-base text @t, capped, returning min(distance, @max_ed + 1).
///
/// @prof, if given, receives @n ints, where prof[j] is the distance of the query against the best text substring ending at position j (inclusive), i.e. the DP's last row.
inline int bpm_semiglobal_ed(const char* q, int m, const char* t, int n, int max_ed,
                             int* prof = nullptr) {
    using namespace bpm_detail;
    if (m <= 0) return 0;
    const BaseCode& code = base_code();
    const int W = (m + 63) / 64;                // query words
    thread_local std::vector<uint64_t> buf;
    buf.assign(static_cast<std::size_t>(W) * (kPeqRows + 2), 0);
    uint64_t* peq = buf.data();                 // kPeqRows x W, with the kOther row left at 0
    uint64_t* Pv = peq + kPeqRows * W;
    uint64_t* Mv = Pv + W;
    for (int i = 0; i < m; ++i) {
        const uint8_t c = code.c[static_cast<unsigned char>(q[i])];
        if (c < kOther) peq[c * W + (i >> 6)] |= uint64_t(1) << (i & 63);
    }
    for (int w = 0; w < W; ++w) Pv[w] = ~uint64_t(0);
    const uint64_t last = uint64_t(1) << ((m - 1) & 63);
    int score = m, best = m;
    for (int j = 0; j < n; ++j) {
        const uint64_t* eq = peq + code.c[static_cast<unsigned char>(t[j])] * W;
        int hin = 0;                            // free text start, so row 0 is 0
        for (int w = 0; w < W; ++w) {
            uint64_t Eq = eq[w];
            const uint64_t pv = Pv[w], mv = Mv[w];
            const uint64_t hneg = hin < 0 ? 1u : 0u;
            const uint64_t Xv = Eq | mv;
            Eq |= hneg;
            const uint64_t Xh = (((Eq & pv) + pv) ^ pv) | Eq;
            uint64_t Ph = mv | ~(Xh | pv);
            uint64_t Mh = pv & Xh;
            const uint64_t hb = w == W - 1 ? last : (uint64_t(1) << 63);
            const int hout = (Ph & hb) ? 1 : (Mh & hb) ? -1 : 0;
            Ph = (Ph << 1) | (hin > 0 ? 1u : 0u);
            Mh = (Mh << 1) | hneg;
            Pv[w] = Mh | ~(Xv | Ph);
            Mv[w] = Ph & Xv;
            hin = hout;
        }
        score += hin;
        if (score < best) best = score;
        if (prof) prof[j] = score;
    }
    return best > max_ed ? max_ed + 1 : best;
}

/// BpmPlacements - the best placement of a query in one text and the best distinct alternative.
///
/// bpm_placements() fills it from the last-row profile @prof: the best edit distance and its leftmost end, and the best edit distance among ends more than @sep positions away from that end, a distinct placement such as one tandem-repeat period over.
/// Distances are capped at @max_ed + 1, and ends are exclusive text offsets, end = j + 1.
struct BpmPlacements { int best, end, alt, alt_end; };
inline BpmPlacements bpm_placements(const int* prof, int n, int sep, int max_ed) {
    BpmPlacements r{max_ed + 1, 0, max_ed + 1, 0};
    for (int j = 0; j < n; ++j)
        if (prof[j] < r.best) { r.best = prof[j]; r.end = j + 1; }
    if (r.best > max_ed) { r.best = max_ed + 1; return r; }
    for (int j = 0; j < n; ++j) {
        const int d = j + 1 - r.end;
        if ((d > sep || d < -sep) && prof[j] < r.alt) { r.alt = prof[j]; r.alt_end = j + 1; }
    }
    if (r.alt > max_ed) r.alt = max_ed + 1;
    return r;
}

/// bpm_semiglobal_ed_x4() - four texts against one query at once, with out[i] = bpm_semiglobal_ed(q, m, t[i], n[i], max_ed) for i < @cnt, where
/// @cnt <= 4.
///
/// With AVX2 each text runs in its own 64-bit lane of the same instruction stream.
/// The query's match vectors are shared, and the words of one lane are linked only through the +-1 horizontal delta, never through an add carry, so lanes stay independent.
/// A lane past the end of its text stops updating its minimum.
/// @prof, if given, receives in prof[i] text i's last-row profile of n[i] ints.
/// Without AVX2 the texts run one by one.
inline void bpm_semiglobal_ed_x4(const char* q, int m, const char* const* t, const int* n,
                                 int cnt, int max_ed, int* out, int* const* prof = nullptr) {
#ifdef __AVX2__
    using namespace bpm_detail;
    constexpr int kLanes = 4;
    if (m <= 0 || cnt <= 1) {
        for (int i = 0; i < cnt; ++i)
            out[i] = bpm_semiglobal_ed(q, m, t[i], n[i], max_ed, prof ? prof[i] : nullptr);
        return;
    }
    const BaseCode& code = base_code();
    const int W = (m + 63) / 64;                // query words
    thread_local std::vector<uint64_t> peq_buf;
    thread_local std::vector<uint64_t> delta_buf;   // Pv then Mv, W vectors of kLanes lanes each
    peq_buf.assign(static_cast<std::size_t>(W) * kPeqRows, 0);
    delta_buf.resize(static_cast<std::size_t>(W) * 2 * kLanes);
    uint64_t* peq = peq_buf.data();             // kPeqRows x W, with the kOther row left at 0
    auto vec = [&](int i) { return reinterpret_cast<__m256i*>(delta_buf.data() + kLanes * i); };
    for (int i = 0; i < m; ++i) {
        const uint8_t c = code.c[static_cast<unsigned char>(q[i])];
        if (c < kOther) peq[c * W + (i >> 6)] |= uint64_t(1) << (i & 63);
    }
    const __m256i ones = _mm256_set1_epi64x(-1), one = _mm256_set1_epi64x(1);
    const __m256i hbit = _mm256_set1_epi64x(static_cast<long long>(uint64_t(1) << 63));
    const __m256i lastbit = _mm256_set1_epi64x(static_cast<long long>(uint64_t(1) << ((m - 1) & 63)));
    for (int w = 0; w < W; ++w) {
        _mm256_storeu_si256(vec(w), ones);
        _mm256_storeu_si256(vec(W + w), _mm256_setzero_si256());
    }
    const char* tp[kLanes];
    int len[kLanes], nmax = 0;
    for (int i = 0; i < kLanes; ++i) {
        tp[i] = i < cnt ? t[i] : t[0];
        len[i] = i < cnt ? n[i] : 0;
        nmax = len[i] > nmax ? len[i] : nmax;
    }
    const __m256i lenv = _mm256_set_epi64x(len[3], len[2], len[1], len[0]);
    __m256i score = _mm256_set1_epi64x(m), best = score;
    for (int j = 0; j < nmax; ++j) {
        // per-lane text code, kOther past a lane's end so it matches nothing
        unsigned cj[kLanes];
        for (int i = 0; i < kLanes; ++i)
            cj[i] = j < len[i] ? code.c[static_cast<unsigned char>(tp[i][j])] : unsigned(kOther);
        __m256i hpos = _mm256_setzero_si256(), hneg = _mm256_setzero_si256();
        for (int w = 0; w < W; ++w) {
            __m256i Eq = _mm256_set_epi64x(
                static_cast<long long>(peq[cj[3] * W + w]), static_cast<long long>(peq[cj[2] * W + w]),
                static_cast<long long>(peq[cj[1] * W + w]), static_cast<long long>(peq[cj[0] * W + w]));
            const __m256i pv = _mm256_loadu_si256(vec(w)), mv = _mm256_loadu_si256(vec(W + w));
            const __m256i hneg1 = _mm256_and_si256(hneg, one);
            const __m256i Xv = _mm256_or_si256(Eq, mv);
            Eq = _mm256_or_si256(Eq, hneg1);
            const __m256i Xh = _mm256_or_si256(
                _mm256_xor_si256(_mm256_add_epi64(_mm256_and_si256(Eq, pv), pv), pv), Eq);
            __m256i Ph = _mm256_or_si256(mv, _mm256_xor_si256(_mm256_or_si256(Xh, pv), ones));
            __m256i Mh = _mm256_and_si256(pv, Xh);
            const __m256i hb = w == W - 1 ? lastbit : hbit;
            const __m256i op = _mm256_cmpeq_epi64(_mm256_and_si256(Ph, hb), hb);
            const __m256i on = _mm256_cmpeq_epi64(_mm256_and_si256(Mh, hb), hb);
            Ph = _mm256_or_si256(_mm256_slli_epi64(Ph, 1), _mm256_and_si256(hpos, one));
            Mh = _mm256_or_si256(_mm256_slli_epi64(Mh, 1), hneg1);
            _mm256_storeu_si256(vec(w), _mm256_or_si256(Mh, _mm256_xor_si256(_mm256_or_si256(Xv, Ph), ones)));
            _mm256_storeu_si256(vec(W + w), _mm256_and_si256(Ph, Xv));
            hpos = op; hneg = on;
        }
        // score += hout, +1 where hpos = -1 and -1 where hneg = -1, with exclusive masks
        score = _mm256_add_epi64(_mm256_sub_epi64(score, hpos), hneg);
        const __m256i active = _mm256_cmpgt_epi64(lenv, _mm256_set1_epi64x(j));
        const __m256i lower = _mm256_and_si256(_mm256_cmpgt_epi64(best, score), active);
        best = _mm256_blendv_epi8(best, score, lower);
        if (prof) {
            alignas(32) long long sc[kLanes];
            _mm256_store_si256(reinterpret_cast<__m256i*>(sc), score);
            for (int i = 0; i < cnt; ++i)
                if (j < len[i]) prof[i][j] = static_cast<int>(sc[i]);
        }
    }
    alignas(32) long long b[kLanes];
    _mm256_store_si256(reinterpret_cast<__m256i*>(b), best);
    for (int i = 0; i < cnt; ++i) out[i] = b[i] > max_ed ? max_ed + 1 : static_cast<int>(b[i]);
#else
    for (int i = 0; i < cnt; ++i)
        out[i] = bpm_semiglobal_ed(q, m, t[i], n[i], max_ed, prof ? prof[i] : nullptr);
#endif
}

}
