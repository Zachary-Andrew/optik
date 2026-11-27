#pragma once

/**
 * @file
 * @brief Exact gap-affine wavefront alignment, WFA, and its O(s) memory form, BiWFA.
 *
 * WFA comes from Marco-Sola, Moure, Moreto and Espinosa 2021, and BiWFA from Marco-Sola, Eizenga, Guarracino, Paten, Garrison and Moreto 2023.
 * Both are written from scratch here.
 * bench/cpp/wfa_affine_stress.cpp and bench/cpp/wfa2_oracle.cpp validate them against an independent Gotoh DP and against WFA2-lib.
 *
 * The score is minimised: a match costs 0, a mismatch x, and a gap of length L costs o + e * L.
 * The pattern (a read of length m) is aligned end to end, whereas the text (a reference window of length n) may have free leading and trailing parts.
 * So an alignment may start at text offset h in [0, text_begin_free] and end at h in [n - text_end_free, n], matching WFA2-lib's ends-free semantics.
 *
 * @par Wavefronts
 * For score s and diagonal k = h - v, with h the text position and v the pattern position, M/I/D hold the furthest-reaching text offset h.
 *   I[s][k] = max(M[s-o-e][k-1], I[s-e][k-1]) + 1     gap in pattern
 *   D[s][k] = max(M[s-o-e][k+1], D[s-e][k+1])         gap in text
 *   M[s][k] = max(M[s-x][k] + 1, I[s][k], D[s][k])
 * and then M is extended along exact matches.
 * Each score's wavefront is a contiguous int32 array over [lo, hi], so diagonals are dense, SIMD-friendly and need no hashing.
 * In the code, M/I/D are Wave::match_offs/ins_offs/del_offs, [lo, hi] is [diag_lo, diag_hi], and x, o, e are AffinePenalties::mismatch, gap_open, gap_extend.
 *
 * Extension, where we depart from WFA2-lib.
 * Pattern and text are copied into byte buffers padded with mutually distinct sentinels, and N maps to a different byte in each.
 * A match run is then found with one AVX2 32-byte compare, movemask and ctz per 32 bases, with no bounds or N checks.
 * WFA2-lib, for comparison, compares 8 bytes per step.
 * A 2-bit packed variant (32 bases per 64-bit XOR) measured slower, because extracting 32 bases at an arbitrary offset costs two loads and a 128-bit shift per step.
 *
 * WfaAffine keeps every wavefront for the traceback.
 * Its O(s^2) memory is fine for read-length problems, since the mapper aligns between chain anchors where s stays small.
 * Gaps between anchors longer than about 1 kb go to BiWfa instead.
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include "bits.hpp"
#include "kmer.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace tpoptoa {

struct AffinePenalties {
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 2;
};

/// Gap state at an alignment end, used by BiWFA's split subproblems so a gap crossing a breakpoint pays its open cost once.
/// kText is a gap consuming only text (a CIGAR D), while kPat consumes only pattern (a CIGAR I).
enum class GapEnd : uint8_t { kNone = 0, kText = 1, kPat = 2 };

/// kNoScore marks a result that exceeded max_score, and kNoScoreCap is the max_score meaning "unbounded".
constexpr int kNoScore = -1;
constexpr int kNoScoreCap = -1;
/// With AVX2, extend() compares kAvx2Bytes per step, and the bytes are all equal when the movemask is kAllBytesEqual.
/// Otherwise it compares kScalarCompareBytes.
constexpr int kAvx2Bytes = 32;
constexpr uint32_t kAllBytesEqual = 0xFFFFFFFFu;
constexpr int kScalarCompareBytes = 8;
/// ctz of a byte-mismatch mask, shifted right by this, gives the index of the first differing byte.
constexpr int kBitsPerByteLog2 = 3;
/// extend_x8() works on kAvx2IntLanes int lanes and gathers kGatherBases bases per lane.
/// It finds the lowest set bit through the float layout, hence the mantissa width, exponent bias and exponent mask.
constexpr int kAvx2IntLanes = 8;
constexpr int kGatherBases = 4;
constexpr int kFloatMantissaBits = 23;
constexpr int kFloatExponentBias = 127;
constexpr int kFloatExponentMask = 0xFF;
/// BiWfa reads five wavefront sources, since a mismatch comes from M, an insertion from M and I, and a deletion from M and D.
/// Together they cut a wave into at most kMaxCuts segments (two per source plus the wave's ends).
constexpr int kSrcMismatch = 0;
constexpr int kSrcInsOpen = 1;
constexpr int kSrcInsExtend = 2;
constexpr int kSrcDelOpen = 3;
constexpr int kSrcDelExtend = 4;
constexpr int kNumSources = 5;
constexpr int kMaxCuts = 12;
constexpr int kDefaultBaseScore = 512;
/// load() keeps separate tables for the pattern (non-ACGT becomes 'N') and the text ('n').
constexpr int kNumLoadTables = 2;

struct WfaAffineResult {
    int score = kNoScore;          // stays kNoScore when max_score was exceeded
    int text_begin = 0;      // first text base of the alignment
    int text_end = 0;        // one past the last text base
    std::string cigar;       // M for matches, X for mismatches, plus I/D
};

class BiWfa;

class WfaAffine {
    friend class BiWfa;
public:
    explicit WfaAffine(AffinePenalties config = {}) : pen_(config) {}

    /// The score comes back as kNoScore when it exceeds @p max_score, and a negative @p max_score means unbounded.
    WfaAffineResult align(const char* pattern, int pattern_len, const char* text, int text_len,
                          int text_begin_free, int text_end_free,
                          int max_score = kNoScoreCap, bool want_cigar = true) {
        return run(pattern, pattern_len, text, text_len, text_begin_free, text_end_free, max_score,
                   want_cigar, GapEnd::kNone, GapEnd::kNone);
    }

    /// End-to-end alignment with gap-state constraints at either end.
    /// `start_free` lets the alignment begin inside a gap of that kind whose open cost was already paid, though it may also begin normally.
    /// `end_in`, by contrast, forces it to end in a gap of that kind and pays the open here.
    WfaAffineResult align_global(const char* pattern, int pattern_len, const char* text, int text_len,
                                 int max_score, GapEnd start_free, GapEnd end_in) {
        return run(pattern, pattern_len, text, text_len, 0, 0, max_score, true, start_free, end_in);
    }

    AffinePenalties penalties() const noexcept { return pen_; }

private:
    WfaAffineResult run(const char* pattern, int pattern_len, const char* text, int text_len,
                        int text_begin_free, int text_end_free, int max_score,
                        bool want_cigar, GapEnd start_free, GapEnd end_in) {
        WfaAffineResult result;

        load(pattern, pattern_len, pattern_buf_, 'N', '!');      // N and padding never match text
        load(text, text_len, text_buf_, 'n', '?');
        pattern_len_ = pattern_len; text_len_ = text_len;
        text_begin_free_ = std::min(std::max(text_begin_free, 0), text_len);
        text_end_free_ = std::min(std::max(text_end_free, 0), text_len);

        // A forced end gap needs a sequence for it to consume.
        if ((end_in == GapEnd::kPat && pattern_len == 0) || (end_in == GapEnd::kText && text_len == 0)) return result;

        waves_.clear();
        reset_arena();
        const int mismatch = pen_.mismatch, open_extend = pen_.gap_open + pen_.gap_extend, gap_extend = pen_.gap_extend;
        // Termination diagonals, where v == m and h in [n - tef, n], so k = h - m.
        const int kend_lo = text_len_ - text_end_free_ - pattern_len_, kend_hi = text_len_ - pattern_len_;

        for (int score = 0;; ++score) {
            Wave wave;

            if (score == 0) {
                // Every allowed start diagonal k = h in [0, tbf], at offset h.
                wave.diag_lo = 0; wave.diag_hi = text_begin_free_;
                wave.match_offs = alloc(wave.width());
                int* match_offs = wave.match_offs;

                for (int diag = wave.diag_lo; diag <= wave.diag_hi; ++diag) match_offs[diag - wave.diag_lo] = diag;

                // Only split subproblems begin inside an already-open gap.
                if (start_free == GapEnd::kText) { wave.ins_offs = alloc(wave.width()); wave.ins_offs[0] = 0; }

                if (start_free == GapEnd::kPat)  { wave.del_offs = alloc(wave.width()); wave.del_offs[0] = 0; }
            } else {
                compute(score, mismatch, open_extend, gap_extend, wave);
            }

            // WFA's extend step slides each new wave along exact matches.
            waves_.push_back(wave);

            if (wave.diag_lo <= wave.diag_hi) {
                int* match_offs = wave.match_offs;

                for (int diag = wave.diag_lo; diag <= wave.diag_hi; ++diag) {
                    int& offset = match_offs[diag - wave.diag_lo];

                    if (offset >= 0) offset = extend(offset, diag);
                }

                // Done once an allowed end diagonal reaches the pattern's end in the required component.
                const int end_lo = std::max(wave.diag_lo, kend_lo), end_hi = std::min(wave.diag_hi, kend_hi);
                const int* endc = end_in == GapEnd::kNone ? wave.match_offs
                                : end_in == GapEnd::kText ? wave.ins_offs : wave.del_offs;

                for (int diag = end_lo; endc && diag <= end_hi; ++diag) {
                    const int offset = endc[diag - wave.diag_lo];

                    if (offset >= 0 && offset - diag == pattern_len_) {
                        result.score = score; result.text_end = offset;

                        if (want_cigar) backtrace(result, score, diag, end_in);

                        return result;
                    }
                }
            }

            // max_score stops the wavefront early on alignments the caller would discard anyway.
            if (max_score >= 0 && score >= max_score) break;
        }

        return result;
    }

    /// Far enough below any real offset that adding penalties can never wrap.
    static constexpr int kUnreachable = INT_MIN / 4;
    static constexpr std::size_t kRetainCells = std::size_t(8) << 20;   // 32 MiB of offsets
    static constexpr std::size_t kBlockCells = std::size_t(1) << 20;    // 4 MiB blocks

    /// M/I/D point into the arena, or are null when that component is absent.
    struct Wave {
        int diag_lo = 1, diag_hi = 0;
        int* match_offs = nullptr; int* ins_offs = nullptr; int* del_offs = nullptr;
        int width() const noexcept { return diag_hi - diag_lo + 1; }
    };

    /// Block arena that keeps pointers stable, never copies on growth and reuses blocks across calls.
    /// Long-lived (thread_local) aligners keep at most kRetainCells between calls, so one rare high-score alignment can't pin its peak footprint for the rest of the run.
    struct Block { std::unique_ptr<int[]> cells; std::size_t cap; };
    std::vector<Block> blocks_;
    std::size_t block_index_ = 0, used_ = 0;

    void reset_arena() {
        std::size_t keep = 0, cells = 0;

        // reset() keeps the first blocks up to kRetainCells for the next alignment and frees anything beyond.
        while (keep < blocks_.size() && cells + blocks_[keep].cap <= kRetainCells) cells += blocks_[keep++].cap;

        blocks_.resize(keep);
        block_index_ = 0; used_ = 0;
    }
    int* alloc(int width) {
        const std::size_t count = static_cast<std::size_t>(width);

        // A block too small for the wave slice is skipped rather than split, since slices must be contiguous.
        while (block_index_ < blocks_.size() && used_ + count > blocks_[block_index_].cap) { ++block_index_; used_ = 0; }

        if (block_index_ == blocks_.size()) {
            const std::size_t cap = std::max(count, kBlockCells);

            blocks_.push_back(Block{std::unique_ptr<int[]>(new int[cap]), cap});
            used_ = 0;
        }

        int* slice = blocks_[block_index_].cells.get() + used_;

        used_ += count;
        // Fresh cells start unreachable until a recurrence feeds them.
        std::fill(slice, slice + count, kUnreachable);
        return slice;
    }
    const Wave* get(int score) const noexcept {
        if (score < 0 || score >= static_cast<int>(waves_.size())) return nullptr;

        const Wave& wave = waves_[score];

        return wave.diag_lo <= wave.diag_hi ? &wave : nullptr;
    }
    static int offset_at(const Wave* wave, const int* cls, int diag) noexcept {
        if (!wave || !cls || diag < wave->diag_lo || diag > wave->diag_hi) return kUnreachable;

        return cls[diag - wave->diag_lo];
    }

    /// Branch-free dst[k] = max(dst[k], src[k + shift] + add) over the overlap.
    void fold(int* dst, int diag_lo, int diag_hi, const Wave* src_wave, const int* cls, int shift, int add) {
        if (!src_wave || !cls) return;

        // Only diagonals present in both waves are folded, and shift maps source diagonal k onto target k + shift.
        const int from = std::max(diag_lo, src_wave->diag_lo - shift), to = std::min(diag_hi, src_wave->diag_hi - shift);

        if (from > to) return;

        const int* __restrict src = cls + (from + shift - src_wave->diag_lo);
        int* __restrict out = dst + (from - diag_lo);     // distinct waves never overlap
        const int len = to - from + 1;

        for (int i = 0; i < len; ++i) {
            const int cand = src[i] + add;

            out[i] = cand > out[i] ? cand : out[i];
        }
    }

    /// Applies the recurrences from the file comment, but doesn't extend yet.
    void compute(int score, int mismatch, int open_extend, int gap_extend, Wave& wave) {
        const Wave* from_mismatch = get(score - mismatch);
        const Wave* moe = get(score - open_extend);
        const Wave* from_extend = get(score - gap_extend);
        const Wave* from_del_extend = from_extend;
        const bool haveI = moe || (from_extend && from_extend->ins_offs);
        const bool haveD = moe || (from_del_extend && from_del_extend->del_offs);
        int diag_lo = INT_MAX, diag_hi = INT_MIN;

        // A mismatch stays on diagonal k while a gap moves to k - 1 or k + 1, which widens the new wave by one on each side.
        if (from_mismatch)  { diag_lo = std::min(diag_lo, from_mismatch->diag_lo);      diag_hi = std::max(diag_hi, from_mismatch->diag_hi); }

        if (moe) { diag_lo = std::min(diag_lo, moe->diag_lo - 1); diag_hi = std::max(diag_hi, moe->diag_hi + 1); }

        if (from_extend && from_extend->ins_offs) { diag_lo = std::min(diag_lo, from_extend->diag_lo + 1); diag_hi = std::max(diag_hi, from_extend->diag_hi + 1); }

        if (from_del_extend && from_del_extend->del_offs) { diag_lo = std::min(diag_lo, from_del_extend->diag_lo - 1); diag_hi = std::max(diag_hi, from_del_extend->diag_hi - 1); }

        if (diag_lo > diag_hi) return;

        // Diagonals outside the matrix are clipped.
        diag_lo = std::max(diag_lo, -pattern_len_); diag_hi = std::min(diag_hi, text_len_);

        if (diag_lo > diag_hi) return;

        wave.diag_lo = diag_lo; wave.diag_hi = diag_hi;
        const int width = diag_hi - diag_lo + 1;

        wave.match_offs = alloc(width);

        if (haveI) wave.ins_offs = alloc(width);

        if (haveD) wave.del_offs = alloc(width);

        int* match_offs = wave.match_offs;
        int* ins_offs = haveI ? wave.ins_offs : nullptr;
        int* del_offs = haveD ? wave.del_offs : nullptr;

        // I[k] = max(M[s-oe][k-1], I[s-e][k-1]) + 1
        if (haveI) {
            fold(ins_offs, diag_lo, diag_hi, moe, moe ? moe->match_offs : nullptr, -1, 1);
            fold(ins_offs, diag_lo, diag_hi, from_extend, from_extend ? from_extend->ins_offs : nullptr, -1, 1);
        }

        // D[k] = max(M[s-oe][k+1], D[s-e][k+1])
        if (haveD) {
            fold(del_offs, diag_lo, diag_hi, moe, moe ? moe->match_offs : nullptr, +1, 0);
            fold(del_offs, diag_lo, diag_hi, from_del_extend, from_del_extend ? from_del_extend->del_offs : nullptr, +1, 0);
        }

        // M[k] = M[s-x][k] + 1
        fold(match_offs, diag_lo, diag_hi, from_mismatch, from_mismatch ? from_mismatch->match_offs : nullptr, 0, 1);
        // A cell is valid when 0 <= v = h-k <= m and h <= n, and after that check M = max(M, I, D).
        // The pass is branch-free, with one specialised loop per combination of present sources so each one vectorises.
        const int text_len = text_len_, pattern_len = pattern_len_;

        int* __restrict match_row = match_offs;

        if (haveI && haveD) {
            int* __restrict ins_row = ins_offs; int* __restrict del_row = del_offs;

            for (int i = 0; i < width; ++i) {
                const int diag = diag_lo + i;
                int match_val = match_row[i], ins_val = ins_row[i], del_val = del_row[i];

                match_val = (match_val > text_len || match_val - diag > pattern_len) ? kUnreachable : match_val;
                ins_val = (ins_val > text_len || ins_val - diag > pattern_len) ? kUnreachable : ins_val;
                del_val = (del_val - diag > pattern_len || del_val - diag < 0) ? kUnreachable : del_val;
                ins_row[i] = ins_val; del_row[i] = del_val;
                match_val = ins_val > match_val ? ins_val : match_val;
                match_val = del_val > match_val ? del_val : match_val;
                match_row[i] = match_val < 0 ? kUnreachable : match_val;
            }
        } else {
            // With only one gap source present, haveI and haveD are the same for every cell, so the compiler can hoist both tests out of the loop.
            for (int i = 0; i < width; ++i) {
                const int diag = diag_lo + i;
                int match_val = match_row[i];

                match_val = (match_val > text_len || match_val - diag > pattern_len) ? kUnreachable : match_val;

                // Present gap components are clipped and folded into M.
                if (haveI) {
                    int ins_val = ins_offs[i];

                    ins_val = (ins_val > text_len || ins_val - diag > pattern_len) ? kUnreachable : ins_val;
                    ins_offs[i] = ins_val;
                    match_val = ins_val > match_val ? ins_val : match_val;
                }

                if (haveD) {
                    int del_val = del_offs[i];

                    del_val = (del_val - diag > pattern_len || del_val - diag < 0) ? kUnreachable : del_val;
                    del_offs[i] = del_val;
                    match_val = del_val > match_val ? del_val : match_val;
                }

                match_row[i] = match_val < 0 ? kUnreachable : match_val;
            }
        }
    }

    static constexpr int kPad = 64;           // sentinel bytes after a buffer, at least one AVX2 load
    static constexpr int kAsciiCaseOffset = 'a' - 'A';

    /// Copies @p seq with ACGT upper-cased, any other base turned into @p nch, and kPad bytes of `pad` appended.
    /// Pattern and text use different N and pad bytes, so a compare can never match an N or run off either end.
    /// That is what lets extend() skip bounds and N checks entirely.
    static void load(const char* seq, int len, std::vector<char>& out, char nch, char pad) {
        // A byte-indexed table keeps this branch-free, because a switch mispredicts on random bases.
        // Building it at compile time also spares every load a first-call initialisation check.
        struct Lut {
            char table[kNumLoadTables][kByteValues] = {};

            constexpr Lut() {
                constexpr char kBases[] = "ACGT";

                for (int w = 0; w < kNumLoadTables; ++w) {
                    for (int i = 0; i < kByteValues; ++i) table[w][i] = w ? 'n' : 'N';

                    for (int b = 0; b < int(kNumBases); ++b) {
                        table[w][static_cast<unsigned char>(kBases[b])] = kBases[b];
                        table[w][static_cast<unsigned char>(kBases[b] + kAsciiCaseOffset)] = kBases[b];
                    }
                }
            }
        };
        static constexpr Lut lut{};
        // Ambiguous bases map to N in the pattern but n in the text, and the differing bytes mean two ambiguous bases never match.
        const char* tab = lut.table[nch == 'n'];

        out.resize(static_cast<std::size_t>(len) + kPad);
        char* dst = out.data();

        for (int i = 0; i < len; ++i) dst[i] = tab[static_cast<unsigned char>(seq[i])];

        std::memset(out.data() + len, pad, kPad);
    }
    /// Follows exact matches along @p diag, kAvx2Bytes bases per step with AVX2 or kScalarCompareBytes otherwise.
    int extend(int offset, int diag) const noexcept {
        const char* pattern_ptr = pattern_buf_.data() + (offset - diag);
        const char* text_ptr = text_buf_.data() + offset;

#ifdef __AVX2__
        for (;;) {
            const __m256i pattern_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pattern_ptr));
            const __m256i text_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(text_ptr));
            const uint32_t eq_mask = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(pattern_vec, text_vec)));

            // Distinct pad bytes guarantee a mismatch before either buffer ends.
            if (eq_mask != kAllBytesEqual) return offset + __builtin_ctz(~eq_mask);

            pattern_ptr += kAvx2Bytes; text_ptr += kAvx2Bytes; offset += kAvx2Bytes;
        }
#else
        for (;;) {
            uint64_t pattern_word, text_word;

            std::memcpy(&pattern_word, pattern_ptr, kScalarCompareBytes); std::memcpy(&text_word, text_ptr, kScalarCompareBytes);
            const uint64_t diff = pattern_word ^ text_word;

            // The lowest set bit of the XOR marks the first unequal byte.
            if (diff) return offset + (__builtin_ctzll(diff) >> kBitsPerByteLog2);

            pattern_ptr += kScalarCompareBytes; text_ptr += kScalarCompareBytes; offset += kScalarCompareBytes;
        }
#endif
    }

    void backtrace(WfaAffineResult& result, int score, int diag, GapEnd end_in) const {
        std::string ops;                          // built backwards, reversed at the end

        enum { kInM, kInI, kInD } state = end_in == GapEnd::kText ? kInI : end_in == GapEnd::kPat ? kInD : kInM;
        int offset = offset_at(get(score), state == kInM ? waves_[score].match_offs : state == kInI ? waves_[score].ins_offs : waves_[score].del_offs, diag);
        const int mismatch = pen_.mismatch, open_extend = pen_.gap_open + pen_.gap_extend, gap_extend = pen_.gap_extend;

        enum { kViaMismatch = 0, kViaIns = 1, kViaDel = 2, kViaStart = 3 };   // how M was reached

        while (true) {
            if (state == kInM) {
                int from, how;                    // how is one of the kVia* values

                if (score == 0) { from = diag; how = kViaStart; }
                else {
                    const Wave* cur = get(score);
                    const Wave* from_mismatch = get(score - mismatch);
                    const int mxv = offset_at(from_mismatch, from_mismatch ? from_mismatch->match_offs : nullptr, diag);
                    const int via_mismatch = mxv >= 0 ? mxv + 1 : kUnreachable;
                    const int via_ins = offset_at(cur, cur->ins_offs, diag), via_del = offset_at(cur, cur->del_offs, diag);
                    // The M cell is the max of its three sources, and whichever source equals it is the predecessor on the path.
                    const int base = std::max(via_mismatch, std::max(via_ins, via_del));

                    from = base;
                    how = (base == via_mismatch) ? kViaMismatch : (base == via_ins ? kViaIns : kViaDel);
                }

                // Extend covered the offsets in between with exact matches.
                for (int t = offset; t > from; --t) ops.push_back('M');

                offset = from;

                if (how == kViaStart) break;

                if (how == kViaMismatch) { ops.push_back('X'); offset -= 1; score -= mismatch; continue; }

                state = (how == kViaIns) ? kInI : kInD;
                continue;
            }

            if (score == 0) break;                    // inside a free start gap, at (0, 0)

            if (state == kInI) {                       // gap in pattern, consuming text
                const Wave* moe = get(score - open_extend);

                ops.push_back('D');
                const int prev = offset - 1;

                // Opened here when the neighbouring M cell holds the predecessor offset. Extended otherwise.
                if (offset_at(moe, moe ? moe->match_offs : nullptr, diag - 1) == prev) { score -= open_extend; state = kInM; }
                else { score -= gap_extend; state = kInI; }

                diag -= 1; offset = prev;
                continue;
            }

            const Wave* moe = get(score - open_extend);        // kInD is a gap in the text, consuming pattern

            ops.push_back('I');

            // A text gap runs the same test on the other neighbour.
            if (offset_at(moe, moe ? moe->match_offs : nullptr, diag + 1) == offset) { score -= open_extend; state = kInM; }
            else { score -= gap_extend; state = kInD; }

            diag += 1;
        }

        result.text_begin = offset;
        std::reverse(ops.begin(), ops.end());
        result.cigar.clear();

        for (std::size_t run_start = 0; run_start < ops.size();) {
            std::size_t run_end = run_start;

            while (run_end < ops.size() && ops[run_end] == ops[run_start]) ++run_end;

            result.cigar += std::to_string(run_end - run_start);
            result.cigar += ops[run_start];
            run_start = run_end;
        }
    }

    AffinePenalties pen_;
    std::vector<char> pattern_buf_, text_buf_;
    int pattern_len_ = 0, text_len_ = 0, text_begin_free_ = 0, text_end_free_ = 0;
    std::vector<Wave> waves_;
};

/// Optimal end-to-end gap-affine alignment in O(s) memory (BiWFA, Marco-Sola et al. 2023).
///
/// A forward WFA from (0, 0) and a reverse WFA from (m, n) on the reversed sequences take turns, each keeping only its last max(x, o+e) + 1 wavefronts.
/// Forward diagonal k meets reverse diagonal (n - m) - k, and the paths join there once forward offset + reverse offset >= n.
/// An M-M join scores s_f + s_r, while a join inside a gap of one kind scores s_f + s_r - o, since only one open is paid.
/// The search runs until no later pair can beat the best join, and then recurses on the two halves.
/// After a gap join, the left half must end in that gap and the right half may start inside it with the open already paid.
/// Halves scoring at most base_score go to WfaAffine::align_global, whose O(s^2) memory is then bounded by base_score^2.
class BiWfa {
public:
    /// Lower values recurse deeper, and tests go all the way down to 0.
    explicit BiWfa(AffinePenalties config = {}, int base_score = kDefaultBaseScore) : pen_(config), base_(config),
        ring_size_(std::max(config.mismatch, config.gap_open + config.gap_extend) + 1), base_score_(base_score) {}

    /// The score comes back as kNoScore if it exceeds a non-negative @p max_score.
    WfaAffineResult align(const char* pattern, int pattern_len, const char* text, int text_len,
                          int max_score = kNoScoreCap) {
        WfaAffineResult result;

        pattern_.assign(pattern, pattern + pattern_len);
        text_.assign(text, text + text_len);
        ops_.clear();
        // solve() appends each half's ops in order, building the full alignment in ops_.
        const int score = solve(0, pattern_len, 0, text_len, GapEnd::kNone, GapEnd::kNone, max_score);

        if (score < 0) return result;

        result.score = score; result.text_begin = 0; result.text_end = text_len;

        for (std::size_t run_start = 0; run_start < ops_.size();) {
            std::size_t run_end = run_start;

            while (run_end < ops_.size() && ops_[run_end] == ops_[run_start]) ++run_end;

            result.cigar += std::to_string(run_end - run_start);
            result.cigar += ops_[run_start];
            run_start = run_end;
        }

        return result;
    }

    AffinePenalties penalties() const noexcept { return pen_; }

private:
    static constexpr int kUnreachable = INT_MIN / 4;

    /// One wavefront of a windowed engine.
    struct Wave {
        int score = kNoScore, diag_lo = 1, diag_hi = 0;
        int amax = INT_MIN;                // largest antidiagonal h + v over M/I/D
        // Storage is kept across the steps that reuse this ring slot, so a resize only fills the growth.
        // hasI and hasD record whether the I and D components exist at this score.
        std::vector<int> match_offs, ins_offs, del_offs;
        bool hasI = false, hasD = false;
    };

    /// Keeps only the last W waves.
    struct Engine {
        const char* pattern = nullptr; const char* text = nullptr;   // padded copies
        int pattern_len = 0, text_len = 0, mismatch = 0, open_extend = 0, gap_extend = 0, ring_size = 0;
        GapEnd start = GapEnd::kNone;     // free start gap, whose open was paid elsewhere
        GapEnd forced = GapEnd::kNone;    // the alignment must start with a gap of this kind
        int last_score = kNoScore;
        std::vector<Wave> ring;

        const Wave* get(int score) const noexcept {
            if (score < 0 || score > last_score) return nullptr;

            const Wave& wave = ring[static_cast<std::size_t>(score % ring_size)];

            return (wave.score == score && wave.diag_lo <= wave.diag_hi) ? &wave : nullptr;
        }
        int extend(int offset, int diag) const noexcept {
            const char* pattern_ptr = pattern + (offset - diag);
            const char* text_ptr = text + offset;

#ifdef __AVX2__
            for (;;) {                      // the distinct pads guarantee this stops
                const __m256i pattern_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pattern_ptr));
                const __m256i text_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(text_ptr));
                const uint32_t eq_mask = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(pattern_vec, text_vec)));

                // Distinct pad bytes guarantee a mismatch before either buffer ends.
                if (eq_mask != kAllBytesEqual) return offset + __builtin_ctz(~eq_mask);

                pattern_ptr += kAvx2Bytes; text_ptr += kAvx2Bytes; offset += kAvx2Bytes;
            }
#else
            for (;;) {
                uint64_t pattern_word, text_word;

                std::memcpy(&pattern_word, pattern_ptr, kScalarCompareBytes); std::memcpy(&text_word, text_ptr, kScalarCompareBytes);
                const uint64_t diff = pattern_word ^ text_word;

                // The lowest set bit of the XOR marks the first unequal byte.
                if (diff) return offset + (__builtin_ctzll(diff) >> kBitsPerByteLog2);

                pattern_ptr += kScalarCompareBytes; text_ptr += kScalarCompareBytes; offset += kScalarCompareBytes;
            }
#endif
        }
        /// Fills one segment [a, b] of a wave from sources indexed by k.
        /// They are q[0] M(s-x), q[1] M(s-o-e) for I, q[2] I(s-e), q[3] M(s-o-e) for D and q[4] D(s-e).
        /// Being branch-free, the kernel vectorises.
        template <bool HI, bool HD>
        static void kernel(const int* const* sources, int* __restrict match_out, int* __restrict ins_out, int* __restrict del_out,
                           int seg_lo, int seg_hi, int pattern_len, int text_len) noexcept {
            const int* __restrict src_mismatch = sources[kSrcMismatch];

            const int* __restrict pio = sources[kSrcInsOpen];
            const int* __restrict pie = sources[kSrcInsExtend];
            const int* __restrict pdo = sources[kSrcDelOpen];
            const int* __restrict pde = sources[kSrcDelExtend];

            for (int diag = seg_lo; diag <= seg_hi; ++diag) {
                int match_val = src_mismatch[diag] + 1;

                match_val = (match_val > text_len || match_val - diag > pattern_len) ? kUnreachable : match_val;

                // Compiling gap components in only when present gives each combination a branch-free loop.
                if constexpr (HI) {
                    int ins_val = (pio[diag] > pie[diag] ? pio[diag] : pie[diag]) + 1;

                    ins_val = (ins_val > text_len || ins_val - diag > pattern_len) ? kUnreachable : ins_val;
                    ins_out[diag] = ins_val;
                    match_val = match_val > ins_val ? match_val : ins_val;
                }

                if constexpr (HD) {
                    int del_val = pdo[diag] > pde[diag] ? pdo[diag] : pde[diag];

                    del_val = (del_val - diag > pattern_len || del_val - diag < 0) ? kUnreachable : del_val;
                    del_out[diag] = del_val;
                    match_val = match_val > del_val ? match_val : del_val;
                }

                match_out[diag] = match_val < 0 ? kUnreachable : match_val;
            }
        }
#ifdef __AVX2__
        /// Extends kAvx2IntLanes diagonals at once.
        /// This is the profile's hot loop: at 10% error a match run averages about 9 bases, so a per-diagonal 32-byte compare was mostly overhead.
        /// The first kGatherBases bases of each diagonal are gathered from P and T and compared in one step.
        /// The number of leading equal bytes is the lowest set bit of the mismatch mask, read from the float exponent of that isolated bit, divided by kBitsPerByte.
        /// Only diagonals whose gathered bases all match continue in extend(), while invalid (negative) offsets gather from 0 and stay unchanged.
        /// The return value is the first unprocessed diagonal, leaving a tail shorter than kAvx2IntLanes.
        int extend_x8(int* match_offs, int diag_lo, int width, int& amax) const noexcept {
            const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            const __m256i neg1 = _mm256_set1_epi32(-1), zero = _mm256_setzero_si256();
            const __m256i bias = _mm256_set1_epi32(kFloatExponentBias), mexp = _mm256_set1_epi32(kFloatExponentMask);
            __m256i vmax = _mm256_set1_epi32(INT_MIN);
            int group_start = 0;

            for (; group_start + kAvx2IntLanes <= width; group_start += kAvx2IntLanes) {
                const __m256i off = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(match_offs + group_start));
                const __m256i diag_vec = _mm256_add_epi32(_mm256_set1_epi32(diag_lo + group_start), lane);
                const __m256i valid = _mm256_cmpgt_epi32(off, neg1);
                const __m256i text_off = _mm256_and_si256(valid, off);
                const __m256i pattern_off = _mm256_and_si256(valid, _mm256_sub_epi32(off, diag_vec));
                const __m256i pattern_bytes = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pattern), pattern_off, 1);
                const __m256i text_bytes = _mm256_i32gather_epi32(reinterpret_cast<const int*>(text), text_off, 1);
                const __m256i mismatch_bytes = _mm256_xor_si256(_mm256_cmpeq_epi8(pattern_bytes, text_bytes), neg1);    // 0xFF marks each mismatching byte
                const __m256i low = _mm256_and_si256(mismatch_bytes, _mm256_sub_epi32(zero, mismatch_bytes));  // lowest set bit
                __m256i pos = _mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(low)), kFloatMantissaBits);

                pos = _mm256_sub_epi32(_mm256_and_si256(pos, mexp), bias);              // bit index, or -bias when no bit is set
                const __m256i full = _mm256_cmpeq_epi32(mismatch_bytes, zero);                     // every gathered base matched
                __m256i eqc = _mm256_srli_epi32(_mm256_max_epi32(pos, zero), kBitsPerByteLog2);

                eqc = _mm256_blendv_epi8(eqc, _mm256_set1_epi32(kGatherBases), full);
                const __m256i new_off = _mm256_blendv_epi8(off, _mm256_add_epi32(off, eqc), valid);

                _mm256_storeu_si256(reinterpret_cast<__m256i*>(match_offs + group_start), new_off);
                int more = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_and_si256(full, valid)));

                // Lanes whose whole gathered block matched may extend further, so they finish with the scalar extend().
                while (more) {
                    const int lane_idx = __builtin_ctz(static_cast<unsigned>(more));

                    more &= more - 1;
                    match_offs[group_start + lane_idx] = extend(match_offs[group_start + lane_idx], diag_lo + group_start + lane_idx);
                }

                const __m256i fin = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(match_offs + group_start));
                // The anti-diagonal h + v = 2h - k measures how far a wave has come, and BiWFA uses its maximum to skip wave pairs that cannot meet yet.
                const __m256i anti = _mm256_sub_epi32(_mm256_add_epi32(fin, fin), diag_vec);

                vmax = _mm256_max_epi32(vmax, _mm256_blendv_epi8(_mm256_set1_epi32(INT_MIN), anti, valid));
            }

            alignas(kAvx2Bytes) int tmp[kAvx2IntLanes];
            _mm256_store_si256(reinterpret_cast<__m256i*>(tmp), vmax);

            for (int t = 0; t < kAvx2IntLanes; ++t) amax = std::max(amax, tmp[t]);

            return group_start;
        }
#endif
        void init() {
            ring.assign(static_cast<std::size_t>(ring_size), Wave{});
            Wave& wave = ring[0];

            wave.score = 0; wave.diag_lo = 0; wave.diag_hi = 0;
            wave.match_offs.assign(1, 0);

            if (start == GapEnd::kText) { wave.ins_offs.assign(1, 0); wave.hasI = true; }

            if (start == GapEnd::kPat)  { wave.del_offs.assign(1, 0); wave.hasD = true; }

            if (forced == GapEnd::kNone) wave.match_offs[0] = extend(0, 0);

            const int start_offset = wave.match_offs[0];

            wave.amax = start_offset + (start_offset - wave.diag_lo);   // h + v, with v = h - k
            last_score = 0;
        }
        const Wave& step() {
            const int score = ++last_score;
            // The affine recurrences reach back at most ring_size scores. Wave storage is therefore indexed modulo ring_size.
            Wave& wave = ring[static_cast<std::size_t>(score % ring_size)];

            wave.score = score; wave.hasI = wave.hasD = false; wave.diag_lo = 1; wave.diag_hi = 0;
            wave.amax = INT_MIN;
            // With a forced start gap, the score-0 "M" is a virtual origin from which only the forced kind of gap may open.
            const bool virt = forced != GapEnd::kNone;
            const Wave* from_mismatch = (virt && score - mismatch == 0) ? nullptr : get(score - mismatch);
            const Wave* moe = get(score - open_extend);
            const Wave* moeI = (virt && score - open_extend == 0 && forced != GapEnd::kText) ? nullptr : moe;
            const Wave* moeD = (virt && score - open_extend == 0 && forced != GapEnd::kPat) ? nullptr : moe;
            const Wave* from_extend = get(score - gap_extend);
            const bool ext_ins = from_extend && from_extend->hasI, ext_del = from_extend && from_extend->hasD;
            int diag_lo = INT_MAX, diag_hi = INT_MIN;

            // As in the forward aligner, a mismatch source keeps its diagonals and each gap source widens the range by one.
            if (from_mismatch)   { diag_lo = std::min(diag_lo, from_mismatch->diag_lo);       diag_hi = std::max(diag_hi, from_mismatch->diag_hi); }

            if (moeI) { diag_lo = std::min(diag_lo, moeI->diag_lo + 1); diag_hi = std::max(diag_hi, moeI->diag_hi + 1); }

            if (moeD) { diag_lo = std::min(diag_lo, moeD->diag_lo - 1); diag_hi = std::max(diag_hi, moeD->diag_hi - 1); }

            if (ext_ins)   { diag_lo = std::min(diag_lo, from_extend->diag_lo + 1);   diag_hi = std::max(diag_hi, from_extend->diag_hi + 1); }

            if (ext_del)   { diag_lo = std::min(diag_lo, from_extend->diag_lo - 1);   diag_hi = std::max(diag_hi, from_extend->diag_hi - 1); }

            diag_lo = std::max(diag_lo, -pattern_len); diag_hi = std::min(diag_hi, text_len);

            if (diag_lo > diag_hi) return wave;

            wave.diag_lo = diag_lo; wave.diag_hi = diag_hi;
            const int width = diag_hi - diag_lo + 1;
            const bool haveI = moeI || ext_ins, haveD = moeD || ext_del;

            wave.match_offs.resize(static_cast<std::size_t>(width));

            if (haveI) wave.ins_offs.resize(static_cast<std::size_t>(width));

            if (haveD) wave.del_offs.resize(static_cast<std::size_t>(width));

            wave.hasI = haveI; wave.hasD = haveD;
            int* match_offs = wave.match_offs.data();
            int* ins_offs = haveI ? wave.ins_offs.data() : nullptr;
            int* del_offs = haveD ? wave.del_offs.data() : nullptr;

            // One fused pass per wave, because separate fills, five shifted-max folds and a validity pass cost about 7 sweeps in the profile.
            // I comes from M(s-o-e)[k-1] + 1 and I(s-e)[k-1] + 1, D from M(s-o-e)[k+1] and D(s-e)[k+1], and M from M(s-x)[k] + 1.
            // Missing sources read as kUnreachable, offsets past either sequence are invalidated, and finally M = max(M, I, D).
            struct Src { const int* row = nullptr; int first_diag = 1, last_diag = 0; };   // row[k] is valid for k in [first_diag, last_diag]
            auto src = [](const Wave* src_wave, const std::vector<int>& cells, int shift) {
                Src result;

                if (src_wave && src_wave->diag_lo <= src_wave->diag_hi && !cells.empty()) {
                    result.row = cells.data() - src_wave->diag_lo + shift; result.first_diag = src_wave->diag_lo - shift; result.last_diag = src_wave->diag_hi - shift;
                }

                return result;
            };
            const Src srcs[kNumSources] = {src(from_mismatch, from_mismatch ? from_mismatch->match_offs : wave.match_offs, 0),
                                 src(moeI, moeI ? moeI->match_offs : wave.match_offs, -1), src(ext_ins ? from_extend : nullptr, ext_ins ? from_extend->ins_offs : wave.match_offs, -1),
                                 src(moeD, moeD ? moeD->match_offs : wave.match_offs, +1), src(ext_del ? from_extend : nullptr, ext_del ? from_extend->del_offs : wave.match_offs, +1)};
            // Missing sources, and sources outside their range, read from a row of kUnreachable.
            thread_local std::vector<int> negrow;

            if (negrow.size() < static_cast<std::size_t>(width)) negrow.assign(static_cast<std::size_t>(width), kUnreachable);

            const int* unreachable_row = negrow.data() - diag_lo;                    // valid for k in [diag_lo, diag_hi]
            // Cut the wave at every source range end, so inside each segment every source is either wholly present or wholly absent.
            // The kernel then needs no branches.
            int cut[kMaxCuts], num_cuts = 0;

            cut[num_cuts++] = diag_lo; cut[num_cuts++] = diag_hi + 1;

            for (const Src& source : srcs)
                if (source.row) {
                    if (source.first_diag > diag_lo && source.first_diag <= diag_hi) cut[num_cuts++] = source.first_diag;

                    if (source.last_diag + 1 > diag_lo && source.last_diag + 1 <= diag_hi) cut[num_cuts++] = source.last_diag + 1;
                }

            for (int x = 1; x < num_cuts; ++x)                           // at most kMaxCuts, so insertion sort is fine
                for (int y = x; y > 0 && cut[y - 1] > cut[y]; --y) std::swap(cut[y - 1], cut[y]);

            num_cuts = static_cast<int>(std::unique(cut, cut + num_cuts) - cut);

            for (int c = 0; c + 1 < num_cuts; ++c) {
                const int seg_lo = cut[c], seg_hi = cut[c + 1] - 1;
                const int* sources[kNumSources];

                for (int t = 0; t < kNumSources; ++t)
                    sources[t] = (srcs[t].row && srcs[t].first_diag <= seg_lo && seg_hi <= srcs[t].last_diag) ? srcs[t].row : unreachable_row;

                if (haveI && haveD)  kernel<true, true>(sources, match_offs - diag_lo, ins_offs - diag_lo, del_offs - diag_lo, seg_lo, seg_hi, pattern_len, text_len);
                else if (haveI)      kernel<true, false>(sources, match_offs - diag_lo, ins_offs - diag_lo, nullptr, seg_lo, seg_hi, pattern_len, text_len);
                else if (haveD)      kernel<false, true>(sources, match_offs - diag_lo, nullptr, del_offs - diag_lo, seg_lo, seg_hi, pattern_len, text_len);
                else                 kernel<false, false>(sources, match_offs - diag_lo, nullptr, nullptr, seg_lo, seg_hi, pattern_len, text_len);
            }

            // Extend along matches.
            // Since M >= I and M >= D, M alone bounds the antidiagonal.
            int first_scalar = 0;
#ifdef __AVX2__
            first_scalar = extend_x8(match_offs, diag_lo, width, wave.amax);

#endif
            for (int i = first_scalar; i < width; ++i) {
                if (match_offs[i] < 0) continue;

                const int diag = diag_lo + i;
                const int offset = extend(match_offs[i], diag);

                match_offs[i] = offset;
                wave.amax = std::max(wave.amax, offset + (offset - diag));   // h + v, with v = h - k
            }

            return wave;
        }
    };

    /// Same as WfaAffine::load() but for a slice, optionally reversed.
    static void pad(const char* seq, int len, bool rev, std::vector<char>& out, char nch, char padc) {
        out.resize(static_cast<std::size_t>(len) + WfaAffine::kPad);

        for (int i = 0; i < len; ++i) {
            char ch = seq[rev ? len - 1 - i : i];

            if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - WfaAffine::kAsciiCaseOffset);

            out[static_cast<std::size_t>(i)] = (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') ? ch : nch;
        }

        std::memset(out.data() + len, padc, WfaAffine::kPad);
    }

    /// Returns the optimal score (kNoScore above max_score) and reports the join as forward join_pattern, join_text and bkind.
    int search(int pat_begin, int pat_end, int text_begin, int text_end, GapEnd sfree, GapEnd ein, int max_score,
               int& join_pattern, int& join_text, GapEnd& bkind) {
        const int pattern_len = pat_end - pat_begin, text_len = text_end - text_begin;

        pad(pattern_.data() + pat_begin, pattern_len, false, fwd_pattern_, 'N', '!');
        pad(text_.data() + text_begin, text_len, false, fwd_text_, 'n', '?');
        pad(pattern_.data() + pat_begin, pattern_len, true, rev_pattern_, 'N', '!');
        pad(text_.data() + text_begin, text_len, true, rev_text_, 'n', '?');
        Engine fwd, rev;

        for (Engine* engine : {&fwd, &rev}) {
            engine->pattern_len = pattern_len; engine->text_len = text_len; engine->mismatch = pen_.mismatch; engine->open_extend = pen_.gap_open + pen_.gap_extend; engine->gap_extend = pen_.gap_extend; engine->ring_size = ring_size_;
        }

        fwd.pattern = fwd_pattern_.data(); fwd.text = fwd_text_.data(); fwd.start = sfree;
        rev.pattern = rev_pattern_.data(); rev.text = rev_text_.data(); rev.forced = ein;
        fwd.init(); rev.init();
        int best = INT_MAX;
        const int diag_sum = text_len - pattern_len;
        // Pairs one new wave `wave_a` (direction `a_fwd`) against every live wave of `other`.
        // A join needs h_f + h_r >= n on one diagonal, which means antidiagonals (2h - k) summing to at least m + n.
        // Waves whose maxima fall short can't meet, so the scan is skipped until the fronts get close.
        auto check = [&](const Wave& wave_a, bool a_fwd, const Engine& other) {
            if (wave_a.diag_lo > wave_a.diag_hi) return;

            for (int score_b = std::max(0, other.last_score - ring_size_ + 1); score_b <= other.last_score; ++score_b) {
                const Wave* wave_b = other.get(score_b);

                if (!wave_b) continue;

                if (static_cast<long>(wave_a.amax) + wave_b->amax < static_cast<long>(pattern_len) + text_len) continue;

                const bool skipM = (a_fwd ? (rev.forced != GapEnd::kNone && score_b == 0)
                                          : (rev.forced != GapEnd::kNone && wave_a.score == 0));

                for (int diag = wave_a.diag_lo; diag <= wave_a.diag_hi; ++diag) {
                    // Forward diagonal k meets reverse diagonal (n - m) - k.
                    const int diag_b = diag_sum - diag;

                    if (diag_b < wave_b->diag_lo || diag_b > wave_b->diag_hi) continue;

                    const std::size_t idx_a = static_cast<std::size_t>(diag - wave_a.diag_lo);
                    const std::size_t idx_b = static_cast<std::size_t>(diag_b - wave_b->diag_lo);
                    const int fwd_diag = a_fwd ? diag : diag_b;
                    auto join = [&](int off_a, int off_b, int cost, GapEnd kind) {
                        // Offsets summing to at least n cover the whole text on this diagonal.
                        if (off_a < 0 || off_b < 0 || off_a + off_b < text_len || cost >= best) return;

                        const int fwd_off = a_fwd ? off_a : off_b;

                        best = cost; join_text = fwd_off; join_pattern = fwd_off - fwd_diag; bkind = kind;
                    };

                    if (!skipM) join(wave_a.match_offs[idx_a], wave_b->match_offs[idx_b], wave_a.score + score_b, GapEnd::kNone);

                    if (wave_a.hasI && wave_b->hasI)
                        // Forward and reverse halves each paid the gap open, and a join inside a gap subtracts one of those penalties.
                        join(wave_a.ins_offs[idx_a], wave_b->ins_offs[idx_b], wave_a.score + score_b - pen_.gap_open, GapEnd::kText);

                    if (wave_a.hasD && wave_b->hasD)
                        join(wave_a.del_offs[idx_a], wave_b->del_offs[idx_b], wave_a.score + score_b - pen_.gap_open, GapEnd::kPat);
                }
            }
        };

        check(fwd.ring[0], true, rev);

        for (;;) {
            // No join found from here on can score below this.
            const int bound = std::min(fwd.last_score, rev.last_score) + 1 + std::max(fwd.last_score, rev.last_score) - (ring_size_ - 1) - pen_.gap_open;

            if (best != INT_MAX && bound >= best) break;

            if (max_score >= 0 && bound > max_score) return kNoScore;

            // The lagging direction steps next.
            if (fwd.last_score <= rev.last_score) check(fwd.step(), true, rev);
            else            check(rev.step(), false, fwd);
        }

        if (max_score >= 0 && best > max_score) return kNoScore;

        return best;
    }

    /// Appends ops for the sub-rectangle under the given gap constraints.
    int solve(int pat_begin, int pat_end, int text_begin, int text_end, GapEnd sfree, GapEnd ein, int max_score) {
        // Every alignment costs at most a mismatch per diagonal step plus two gaps.
        // That finite bound lets an infeasible constraint end the search instead of looping.
        const int full_pattern_len = pat_end - pat_begin, full_text_len = text_end - text_begin;
        constexpr int kBoundGaps = 2;
        const int bound = pen_.mismatch * std::min(full_pattern_len, full_text_len) + kBoundGaps * (pen_.gap_open + pen_.gap_extend * std::max(full_pattern_len, full_text_len)) + pen_.gap_open;

        if (max_score < 0 || max_score > bound) max_score = bound;

        int join_pattern = 0, join_text = 0; GapEnd join_kind = GapEnd::kNone;
        const int score = search(pat_begin, pat_end, text_begin, text_end, sfree, ein, max_score, join_pattern, join_text, join_kind);

        if (score < 0) return kNoScore;

        const int pattern_len = pat_end - pat_begin, text_len = text_end - text_begin;
        // A join at a corner would make one half the whole problem again, and the recursion would never terminate.
        const bool corner = (join_pattern == 0 && join_text == 0) || (join_pattern == pattern_len && join_text == text_len);

        if (score <= base_score_ || corner) {
            WfaAffineResult result = base_.align_global(pattern_.data() + pat_begin, pattern_len, text_.data() + text_begin, text_len,
                                                   score, sfree, ein);

            if (result.score < 0) return kNoScore;

            append(result.cigar);
            return result.score;
        }

        // Both halves are solved recursively, with the left half ending in the join's gap and the right half starting inside it.
        if (solve(pat_begin, pat_begin + join_pattern, text_begin, text_begin + join_text, sfree, join_kind, kNoScoreCap) < 0) return kNoScore;

        if (solve(pat_begin + join_pattern, pat_end, text_begin + join_text, text_end, join_kind, ein, kNoScoreCap) < 0) return kNoScore;

        return score;
    }

    void append(const std::string& cig) {
        int run_len = 0;

        for (char ch : cig) {
            if (ch >= '0' && ch <= '9') run_len = run_len * kDecimalBase + (ch - '0');
            else { ops_.append(static_cast<std::size_t>(run_len), ch); run_len = 0; }
        }
    }

    AffinePenalties pen_;
    WfaAffine base_;
    int ring_size_, base_score_;
    std::vector<char> pattern_, text_, fwd_pattern_, fwd_text_, rev_pattern_, rev_text_;
    std::string ops_;
};

}
