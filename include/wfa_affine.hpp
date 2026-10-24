#pragma once

/*
 * wfa_affine.hpp - exact gap-affine wavefront alignment, WFA, and its O(s) memory form, BiWFA.
 *
 * The sources are WFA by Marco-Sola, Moure, Moreto and Espinosa 2021 and BiWFA by Marco-Sola, Eizenga, Guarracino, Paten, Garrison and Moreto 2023.
 * Both are implemented here from scratch and validated against an independent Gotoh DP and WFA2-lib by bench/cpp/wfa_affine_stress.cpp and bench/cpp/wfa2_oracle.cpp.
 *
 * Score model, minimised, with mismatch x, a gap of length L costing o + e * L, and match 0.
 * The pattern, a read of length m, is aligned end to end, while the text, a reference window of length n, may have free leading and trailing parts.
 * The alignment may thus start at text offset h in [0, text_begin_free] and end at h in [n - text_end_free, n], WFA2-lib's ends-free semantics.
 *
 * Wavefronts.
 * For score s and diagonal k = h - v, with h the text position and v the pattern position, M/I/D hold the furthest-reaching text offset h.
 *   I[s][k] = max(M[s-o-e][k-1], I[s-e][k-1]) + 1     gap in pattern
 *   D[s][k] = max(M[s-o-e][k+1], D[s-e][k+1])         gap in text
 *   M[s][k] = max(M[s-x][k] + 1, I[s][k], D[s][k])
 * then M is extended along exact matches.
 * Each score's wavefront is a contiguous int32 array over [lo, hi], with dense diagonals, SIMD-friendly and without hashing.
 *
 * Extension, the modification from WFA2-lib.
 * Pattern and text are copied into byte buffers padded with mutually distinct sentinels, with N mapped to a distinct byte in each, so a match run is found with one AVX2 32-byte compare + movemask + ctz per 32 bases and no bounds or N checks.
 * WFA2-lib compares 8 bytes per step.
 * A 2-bit packed variant with 32 bases per 64-bit XOR measured slower, since extracting 32 bases at an arbitrary offset costs two loads and a 128-bit shift per sequence per step.
 *
 * WfaAffine keeps every wavefront for the traceback.
 * Its O(s^2) memory suits read-length problems, since the mapper aligns between chain anchors and s stays small.
 * BiWfa takes the gaps between anchors over 1 kb.
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace tpoptoa {

struct AffinePenalties {
    int x = 4;   // mismatch
    int o = 6;   // gap open
    int e = 2;   // gap extend
};

/// GapEnd - gap state at an alignment end for BiWFA's split subproblems, so a gap that crosses a breakpoint is charged its open cost once.
/// kText is a gap consuming only text, a CIGAR D, and kPat one consuming only pattern, a CIGAR I.
enum class GapEnd : uint8_t { kNone = 0, kText = 1, kPat = 2 };

struct WfaAffineResult {
    int score = -1;          // -1 when max_score was exceeded
    int text_begin = 0;      // first text base of the alignment
    int text_end = 0;        // one past the last text base
    std::string cigar;       // M/X/I/D ops ('=' matches as 'M', mismatches 'X')
};

class BiWfa;

class WfaAffine {
    friend class BiWfa;
public:
    explicit WfaAffine(AffinePenalties p = {}) : pen_(p) {}

    /// align() - align pattern[0..m) against text[0..n) with the given free text ends.
    /// The result score is -1 if it exceeds @max_score, and @max_score < 0 means unbounded.
    WfaAffineResult align(const char* pattern, int m, const char* text, int n,
                          int text_begin_free, int text_end_free,
                          int max_score = -1, bool want_cigar = true) {
        return run(pattern, m, text, n, text_begin_free, text_end_free, max_score,
                   want_cigar, GapEnd::kNone, GapEnd::kNone);
    }

    /// align_global() - end-to-end alignment with gap-state constraints at the ends.
    /// With `start_free` the alignment may begin inside a gap of that kind whose open cost is already paid, or begin normally.
    /// With `end_in` it must end in a gap of that kind, paying the open here.
    WfaAffineResult align_global(const char* pattern, int m, const char* text, int n,
                                 int max_score, GapEnd start_free, GapEnd end_in) {
        return run(pattern, m, text, n, 0, 0, max_score, true, start_free, end_in);
    }

    AffinePenalties penalties() const noexcept { return pen_; }

private:
    WfaAffineResult run(const char* pattern, int m, const char* text, int n,
                        int text_begin_free, int text_end_free, int max_score,
                        bool want_cigar, GapEnd start_free, GapEnd end_in) {
        WfaAffineResult r;
        load(pattern, m, pattern_buf_, 'N', '!');      // N and padding never match text
        load(text, n, text_buf_, 'n', '?');
        m_ = m; n_ = n;
        text_begin_free_ = std::min(std::max(text_begin_free, 0), n);
        text_end_free_ = std::min(std::max(text_end_free, 0), n);
        // a forced end gap needs a sequence for it to consume
        if ((end_in == GapEnd::kPat && m == 0) || (end_in == GapEnd::kText && n == 0)) return r;
        wf_.clear();
        reset_arena();
        const int x = pen_.x, oe = pen_.o + pen_.e, e = pen_.e;
        // Termination diagonals, where v == m and h in [n - tef, n], so k = h - m
        const int kend_lo = n_ - text_end_free_ - m_, kend_hi = n_ - m_;
        for (int s = 0;; ++s) {
            Wave w;
            if (s == 0) {
                // every allowed start diagonal k = h in [0, tbf], offset h
                w.lo = 0; w.hi = text_begin_free_;
                w.M = alloc(w.width());
                int* M = w.M;
                for (int k = w.lo; k <= w.hi; ++k) M[k - w.lo] = k;
                // begin inside an already-open gap, which only split subproblems do
                if (start_free == GapEnd::kText) { w.I = alloc(w.width()); w.I[0] = 0; }
                if (start_free == GapEnd::kPat)  { w.D = alloc(w.width()); w.D[0] = 0; }
            } else {
                compute(s, x, oe, e, w);
            }
            wf_.push_back(w);
            if (w.lo <= w.hi) {
                int* M = w.M;
                for (int k = w.lo; k <= w.hi; ++k) {
                    int& h = M[k - w.lo];
                    if (h >= 0) h = extend(h, k);
                }
                const int a2 = std::max(w.lo, kend_lo), b2 = std::min(w.hi, kend_hi);
                const int* endc = end_in == GapEnd::kNone ? w.M
                                : end_in == GapEnd::kText ? w.I : w.D;
                for (int k = a2; endc && k <= b2; ++k) {
                    const int h = endc[k - w.lo];
                    if (h >= 0 && h - k == m_) {
                        r.score = s; r.text_end = h;
                        if (want_cigar) backtrace(r, s, k, end_in);
                        return r;
                    }
                }
            }
            if (max_score >= 0 && s >= max_score) break;
        }
        return r;
    }

    /// Offset of an unreachable cell, far enough below any real offset that adding penalties never wraps.
    static constexpr int kUnreachable = INT_MIN / 4;
    static constexpr std::size_t kRetainCells = std::size_t(8) << 20;   // 32 MiB of offsets
    static constexpr std::size_t kBlockCells = std::size_t(1) << 20;    // 4 MiB blocks
    /// A wavefront over diagonals [lo, hi], with M/I/D pointing into the arena or null when absent.
    struct Wave {
        int lo = 1, hi = 0;
        int* M = nullptr; int* I = nullptr; int* D = nullptr;
        int width() const noexcept { return hi - lo + 1; }
    };

    /// Block arena for the wavefronts, which keeps pointers stable, never copies on growth and reuses blocks across calls.
    /// Long-lived (thread_local) aligners keep at most kRetainCells between calls, so a rare high-score alignment does not pin its peak footprint for the rest of the run.
    struct Block { std::unique_ptr<int[]> p; std::size_t cap; };
    std::vector<Block> blocks_;
    std::size_t block_index_ = 0, used_ = 0;
    void reset_arena() {
        std::size_t keep = 0, cells = 0;
        while (keep < blocks_.size() && cells + blocks_[keep].cap <= kRetainCells) cells += blocks_[keep++].cap;
        blocks_.resize(keep);
        block_index_ = 0; used_ = 0;
    }
    int* alloc(int width) {
        const std::size_t n = static_cast<std::size_t>(width);
        while (block_index_ < blocks_.size() && used_ + n > blocks_[block_index_].cap) { ++block_index_; used_ = 0; }
        if (block_index_ == blocks_.size()) {
            const std::size_t cap = std::max(n, kBlockCells);
            blocks_.push_back(Block{std::unique_ptr<int[]>(new int[cap]), cap});
            used_ = 0;
        }
        int* p = blocks_[block_index_].p.get() + used_;
        used_ += n;
        std::fill(p, p + n, kUnreachable);
        return p;
    }
    const Wave* get(int s) const noexcept {
        if (s < 0 || s >= static_cast<int>(wf_.size())) return nullptr;
        const Wave& w = wf_[s];
        return w.lo <= w.hi ? &w : nullptr;
    }
    static int at(const Wave* w, const int* cls, int k) noexcept {
        if (!w || !cls || k < w->lo || k > w->hi) return kUnreachable;
        return cls[k - w->lo];
    }

    /// fold() - dst[k] = max(dst[k], src[k + shift] + add) over the overlap, branch-free.
    void fold(int* dst, int lo, int hi, const Wave* sw, const int* cls, int shift, int add) {
        if (!sw || !cls) return;
        const int a = std::max(lo, sw->lo - shift), b = std::min(hi, sw->hi - shift);
        if (a > b) return;
        const int* __restrict src = cls + (a + shift - sw->lo);
        int* __restrict d = dst + (a - lo);     // distinct waves never overlap
        const int len = b - a + 1;
        for (int i = 0; i < len; ++i) {
            const int v = src[i] + add;
            d[i] = v > d[i] ? v : d[i];
        }
    }

    /// compute() - the wavefront of score @s from the earlier ones (the recurrences in the file comment), before extension.
    void compute(int s, int x, int oe, int e, Wave& w) {
        const Wave* mx = get(s - x);
        const Wave* moe = get(s - oe);
        const Wave* ie = get(s - e);
        const Wave* de = ie;
        const bool haveI = moe || (ie && ie->I);
        const bool haveD = moe || (de && de->D);
        int lo = INT_MAX, hi = INT_MIN;
        if (mx)  { lo = std::min(lo, mx->lo);      hi = std::max(hi, mx->hi); }
        if (moe) { lo = std::min(lo, moe->lo - 1); hi = std::max(hi, moe->hi + 1); }
        if (ie && ie->I) { lo = std::min(lo, ie->lo + 1); hi = std::max(hi, ie->hi + 1); }
        if (de && de->D) { lo = std::min(lo, de->lo - 1); hi = std::max(hi, de->hi - 1); }
        if (lo > hi) return;
        lo = std::max(lo, -m_); hi = std::min(hi, n_);
        if (lo > hi) return;
        w.lo = lo; w.hi = hi;
        const int W = hi - lo + 1;
        w.M = alloc(W);
        if (haveI) w.I = alloc(W);
        if (haveD) w.D = alloc(W);
        int* M = w.M;
        int* I = haveI ? w.I : nullptr;
        int* D = haveD ? w.D : nullptr;
        // I[k] = max(M[s-oe][k-1], I[s-e][k-1]) + 1
        if (haveI) {
            fold(I, lo, hi, moe, moe ? moe->M : nullptr, -1, 1);
            fold(I, lo, hi, ie, ie ? ie->I : nullptr, -1, 1);
        }
        // D[k] = max(M[s-oe][k+1], D[s-e][k+1])
        if (haveD) {
            fold(D, lo, hi, moe, moe ? moe->M : nullptr, +1, 0);
            fold(D, lo, hi, de, de ? de->D : nullptr, +1, 0);
        }
        // M[k] = M[s-x][k] + 1
        fold(M, lo, hi, mx, mx ? mx->M : nullptr, 0, 1);
        // validity, 0 <= v = h-k <= m and h <= n, then M = max(M, I, D).
        // It is branch-free, with one specialised loop per class set so each vectorises
        const int n = n_, m = m_;
        int* __restrict Mr = M;
        if (haveI && haveD) {
            int* __restrict Ir = I; int* __restrict Dr = D;
            for (int i = 0; i < W; ++i) {
                const int k = lo + i;
                int mv = Mr[i], iv = Ir[i], dv = Dr[i];
                mv = (mv > n || mv - k > m) ? kUnreachable : mv;
                iv = (iv > n || iv - k > m) ? kUnreachable : iv;
                dv = (dv - k > m || dv - k < 0) ? kUnreachable : dv;
                Ir[i] = iv; Dr[i] = dv;
                mv = iv > mv ? iv : mv;
                mv = dv > mv ? dv : mv;
                Mr[i] = mv < 0 ? kUnreachable : mv;
            }
        } else {
            for (int i = 0; i < W; ++i) {
                const int k = lo + i;
                int mv = Mr[i];
                mv = (mv > n || mv - k > m) ? kUnreachable : mv;
                if (haveI) {
                    int iv = I[i];
                    iv = (iv > n || iv - k > m) ? kUnreachable : iv;
                    I[i] = iv;
                    mv = iv > mv ? iv : mv;
                }
                if (haveD) {
                    int dv = D[i];
                    dv = (dv - k > m || dv - k < 0) ? kUnreachable : dv;
                    D[i] = dv;
                    mv = dv > mv ? dv : mv;
                }
                Mr[i] = mv < 0 ? kUnreachable : mv;
            }
        }
    }

    static constexpr int kPad = 64;           // sentinel bytes after a buffer, >= one AVX2 load
    static constexpr int kAsciiCaseOffset = 'a' - 'A';
    /// load() - byte copy of @s with ACGT upper-cased and any other base turned into @nch, followed by kPad bytes of `pad`.
    /// Pattern and text use different N and pad bytes, so a compare can never match an N or run past either end, and extend needs no bounds or N checks.
    static void load(const char* s, int len, std::vector<char>& out, char nch, char pad) {
        // Branch-free through a 256-entry table, since a switch mispredicts on random bases Built at compile time, so no first-call initialization check per load.
        struct Lut {
            char t[2][256] = {};
            constexpr Lut() {
                constexpr char kBases[] = "ACGT";
                for (int w = 0; w < 2; ++w) {
                    for (int i = 0; i < 256; ++i) t[w][i] = w ? 'n' : 'N';
                    for (int b = 0; b < 4; ++b) {
                        t[w][static_cast<unsigned char>(kBases[b])] = kBases[b];
                        t[w][static_cast<unsigned char>(kBases[b] + kAsciiCaseOffset)] = kBases[b];
                    }
                }
            }
        };
        static constexpr Lut lut{};
        const char* tab = lut.t[nch == 'n'];
        out.resize(static_cast<std::size_t>(len) + kPad);
        char* o = out.data();
        for (int i = 0; i < len; ++i) o[i] = tab[static_cast<unsigned char>(s[i])];
        std::memset(out.data() + len, pad, kPad);
    }
    /// Extend offset h on diagonal k along exact matches, 32 bases per step with an AVX2 byte compare or 8 with a scalar 64-bit XOR.
    int extend(int h, int k) const noexcept {
        const char* P = pattern_buf_.data() + (h - k);
        const char* T = text_buf_.data() + h;
#ifdef __AVX2__
        for (;;) {
            const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(P));
            const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(T));
            const uint32_t eq = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)));
            if (eq != 0xFFFFFFFFu) return h + __builtin_ctz(~eq);
            P += 32; T += 32; h += 32;
        }
#else
        for (;;) {
            uint64_t a, b;
            std::memcpy(&a, P, 8); std::memcpy(&b, T, 8);
            const uint64_t x = a ^ b;
            if (x) return h + (__builtin_ctzll(x) >> 3);
            P += 8; T += 8; h += 8;
        }
#endif
    }

    void backtrace(WfaAffineResult& r, int s, int k, GapEnd end_in) const {
        std::string ops;                          // reversed
        enum { kInM, kInI, kInD } st = end_in == GapEnd::kText ? kInI : end_in == GapEnd::kPat ? kInD : kInM;
        int h = at(get(s), st == kInM ? wf_[s].M : st == kInI ? wf_[s].I : wf_[s].D, k);
        const int x = pen_.x, oe = pen_.o + pen_.e, e = pen_.e;
        while (true) {
            if (st == kInM) {
                int from, how;                    // 0 X, 1 I, 2 D, 3 start
                if (s == 0) { from = k; how = 3; }
                else {
                    const Wave* cur = get(s);
                    const Wave* mx = get(s - x);
                    const int mxv = at(mx, mx ? mx->M : nullptr, k);
                    const int Xv = mxv >= 0 ? mxv + 1 : kUnreachable;
                    const int Iv = at(cur, cur->I, k), Dv = at(cur, cur->D, k);
                    const int base = std::max(Xv, std::max(Iv, Dv));
                    from = base;
                    how = (base == Xv) ? 0 : (base == Iv ? 1 : 2);
                }
                for (int t = h; t > from; --t) ops.push_back('M');
                h = from;
                if (how == 3) break;
                if (how == 0) { ops.push_back('X'); h -= 1; s -= x; continue; }
                st = (how == 1) ? kInI : kInD;
                continue;
            }
            if (s == 0) break;                    // inside a free start gap, at (0, 0)
            if (st == kInI) {                       // gap in pattern, consuming text
                const Wave* moe = get(s - oe);
                ops.push_back('D');
                const int prev = h - 1;
                if (at(moe, moe ? moe->M : nullptr, k - 1) == prev) { s -= oe; st = kInM; }
                else { s -= e; st = kInI; }
                k -= 1; h = prev;
                continue;
            }
            const Wave* moe = get(s - oe);        // kInD is a gap in the text, consuming pattern
            ops.push_back('I');
            if (at(moe, moe ? moe->M : nullptr, k + 1) == h) { s -= oe; st = kInM; }
            else { s -= e; st = kInD; }
            k += 1;
        }
        r.text_begin = h;
        std::reverse(ops.begin(), ops.end());
        r.cigar.clear();
        for (std::size_t i = 0; i < ops.size();) {
            std::size_t j = i;
            while (j < ops.size() && ops[j] == ops[i]) ++j;
            r.cigar += std::to_string(j - i);
            r.cigar += ops[i];
            i = j;
        }
    }

    AffinePenalties pen_;
    std::vector<char> pattern_buf_, text_buf_;
    int m_ = 0, n_ = 0, text_begin_free_ = 0, text_end_free_ = 0;
    std::vector<Wave> wf_;
};

/// BiWfa - optimal gap-affine end-to-end alignment in O(s) memory, the BiWFA of Marco-Sola et al. 2023.
///
/// A forward WFA from (0, 0) and a reverse WFA from (m, n), run on the reversed sequences, advance alternately, each keeping only the last max(x, o+e) + 1 wavefronts.
/// Forward diagonal k meets reverse diagonal (n - m) - k, and when forward offset + reverse offset >= n the two paths join there, M-M at s_f + s_r, or inside a gap of the same kind at s_f + s_r - o, paying one open instead of two.
/// The search continues until no later pair can beat the best join, then recurses on the two halves.
/// A gap join makes the left half end in that gap and lets the right half start inside it with the open already paid.
/// Halves whose score is at most base_score are solved by WfaAffine::align_global, the traceback WFA, whose O(s^2) memory is then bounded by base_score^2.
class BiWfa {
public:
    /// Halves scoring at most `base_score` are solved by traceback WFA in memory of about base_score^2.
    /// Lower values recurse deeper, down to 0 in tests.
    explicit BiWfa(AffinePenalties p = {}, int base_score = 512) : pen_(p), base_(p),
        W_(std::max(p.x, p.o + p.e) + 1), base_score_(base_score) {}

    /// align() - optimal end-to-end alignment, whose score is -1 if it exceeds the non-negative @max_score.
    WfaAffineResult align(const char* pattern, int m, const char* text, int n,
                          int max_score = -1) {
        WfaAffineResult r;
        P_.assign(pattern, pattern + m);
        T_.assign(text, text + n);
        ops_.clear();
        const int sc = solve(0, m, 0, n, GapEnd::kNone, GapEnd::kNone, max_score);
        if (sc < 0) return r;
        r.score = sc; r.text_begin = 0; r.text_end = n;
        // run-length encode the op string
        for (std::size_t i = 0; i < ops_.size();) {
            std::size_t j = i;
            while (j < ops_.size() && ops_[j] == ops_[i]) ++j;
            r.cigar += std::to_string(j - i);
            r.cigar += ops_[i];
            i = j;
        }
        return r;
    }

    AffinePenalties penalties() const noexcept { return pen_; }

private:
    static constexpr int kUnreachable = INT_MIN / 4;

    /// One wavefront of a windowed engine.
    struct Wave {
        int s = -1, lo = 1, hi = 0;
        int amax = INT_MIN;                // max antidiagonal h + v over M/I/D
        // Storage is kept across the steps that reuse this ring slot, so a resize fills only the growth.
        // hasI and hasD say whether the I and D components exist at this score.
        std::vector<int> M, I, D;
        bool hasI = false, hasD = false;
    };

    /// Windowed WFA over (P, T) from (0, 0), keeping the last W waves.
    struct Engine {
        const char* P = nullptr; const char* T = nullptr;   // padded buffers
        int m = 0, n = 0, x = 0, oe = 0, e = 0, W = 0;
        GapEnd start = GapEnd::kNone;     // free start gap, paid elsewhere
        GapEnd forced = GapEnd::kNone;    // must start with a gap of this kind
        int s = -1;                      // last computed score
        std::vector<Wave> ring;

        const Wave* get(int sc) const noexcept {
            if (sc < 0 || sc > s) return nullptr;
            const Wave& w = ring[static_cast<std::size_t>(sc % W)];
            return (w.s == sc && w.lo <= w.hi) ? &w : nullptr;
        }
        int extend(int h, int k) const noexcept {
            const char* p = P + (h - k);
            const char* t = T + h;
#ifdef __AVX2__
            for (;;) {                      // distinct pads stop the run
                const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
                const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(t));
                const uint32_t eq = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)));
                if (eq != 0xFFFFFFFFu) return h + __builtin_ctz(~eq);
                p += 32; t += 32; h += 32;
            }
#else
            for (;;) {
                uint64_t a, b;
                std::memcpy(&a, p, 8); std::memcpy(&b, t, 8);
                const uint64_t d = a ^ b;
                if (d) return h + (__builtin_ctzll(d) >> 3);
                p += 8; t += 8; h += 8;
            }
#endif
        }
        /// One segment [a, b] of a wave, all sources indexed by k (q[0] M(s-x), q[1] M(s-o-e) for I, q[2] I(s-e), q[3] M(s-o-e) for D, q[4] D(s-e)):
        /// branch-free, so the compiler vectorizes it.
        template <bool HI, bool HD>
        static void kernel(const int* const* q, int* __restrict Mo, int* __restrict Io, int* __restrict Do,
                           int a, int b, int m, int n) noexcept {
            const int* __restrict pm = q[0];
            const int* __restrict pio = q[1];
            const int* __restrict pie = q[2];
            const int* __restrict pdo = q[3];
            const int* __restrict pde = q[4];
            for (int k = a; k <= b; ++k) {
                int mv = pm[k] + 1;
                mv = (mv > n || mv - k > m) ? kUnreachable : mv;
                if constexpr (HI) {
                    int iv = (pio[k] > pie[k] ? pio[k] : pie[k]) + 1;
                    iv = (iv > n || iv - k > m) ? kUnreachable : iv;
                    Io[k] = iv;
                    mv = mv > iv ? mv : iv;
                }
                if constexpr (HD) {
                    int dv = pdo[k] > pde[k] ? pdo[k] : pde[k];
                    dv = (dv - k > m || dv - k < 0) ? kUnreachable : dv;
                    Do[k] = dv;
                    mv = mv > dv ? mv : dv;
                }
                Mo[k] = mv < 0 ? kUnreachable : mv;
            }
        }
#ifdef __AVX2__
        /// Extend 8 diagonals at once.
        /// This is the profile's hot loop, since at 10% error a match run averages about 9 bases and a per-diagonal 32-byte compare was mostly overhead.
        /// The first 4 bases of each diagonal are gathered from P and T and compared in one step, and the count of leading equal bytes is the position of the lowest set bit of the mismatch mask (read from the float exponent of that isolated bit), divided by 8.
        /// Only diagonals whose 4 bases all match continue in extend().
        /// Invalid offsets (< 0) gather from 0 and are left unchanged.
        /// Returns the first diagonal not processed (a tail of
        /// < 8), with amax updated for the processed ones.
        int extend_x8(int* M, int lo, int Wd, int& amax) const noexcept {
            const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            const __m256i neg1 = _mm256_set1_epi32(-1), zero = _mm256_setzero_si256();
            const __m256i bias = _mm256_set1_epi32(127), mexp = _mm256_set1_epi32(0xFF);
            __m256i vmax = _mm256_set1_epi32(INT_MIN);
            int i = 0;
            for (; i + 8 <= Wd; i += 8) {
                const __m256i off = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(M + i));
                const __m256i kv = _mm256_add_epi32(_mm256_set1_epi32(lo + i), lane);
                const __m256i valid = _mm256_cmpgt_epi32(off, neg1);
                const __m256i h = _mm256_and_si256(valid, off);
                const __m256i v = _mm256_and_si256(valid, _mm256_sub_epi32(off, kv));
                const __m256i pv = _mm256_i32gather_epi32(reinterpret_cast<const int*>(P), v, 1);
                const __m256i tv = _mm256_i32gather_epi32(reinterpret_cast<const int*>(T), h, 1);
                const __m256i ne = _mm256_xor_si256(_mm256_cmpeq_epi8(pv, tv), neg1);    // 0xFF per mismatch byte
                const __m256i low = _mm256_and_si256(ne, _mm256_sub_epi32(zero, ne));  // lowest set bit
                __m256i pos = _mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(low)), 23);
                pos = _mm256_sub_epi32(_mm256_and_si256(pos, mexp), bias);              // bit index, -127 if none
                const __m256i full = _mm256_cmpeq_epi32(ne, zero);                     // all 4 bases equal
                __m256i eqc = _mm256_srli_epi32(_mm256_max_epi32(pos, zero), 3);
                eqc = _mm256_blendv_epi8(eqc, _mm256_set1_epi32(4), full);
                const __m256i nh = _mm256_blendv_epi8(off, _mm256_add_epi32(off, eqc), valid);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(M + i), nh);
                int more = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_and_si256(full, valid)));
                while (more) {
                    const int t = __builtin_ctz(static_cast<unsigned>(more));
                    more &= more - 1;
                    M[i + t] = extend(M[i + t], lo + i + t);
                }
                const __m256i fin = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(M + i));
                const __m256i anti = _mm256_sub_epi32(_mm256_add_epi32(fin, fin), kv);
                vmax = _mm256_max_epi32(vmax, _mm256_blendv_epi8(_mm256_set1_epi32(INT_MIN), anti, valid));
            }
            alignas(32) int tmp[8];
            _mm256_store_si256(reinterpret_cast<__m256i*>(tmp), vmax);
            for (int t = 0; t < 8; ++t) amax = std::max(amax, tmp[t]);
            return i;
        }
#endif
        void init() {
            ring.assign(static_cast<std::size_t>(W), Wave{});
            Wave& w = ring[0];
            w.s = 0; w.lo = 0; w.hi = 0;
            w.M.assign(1, 0);
            if (start == GapEnd::kText) { w.I.assign(1, 0); w.hasI = true; }
            if (start == GapEnd::kPat)  { w.D.assign(1, 0); w.hasD = true; }
            if (forced == GapEnd::kNone) w.M[0] = extend(0, 0);
            w.amax = 2 * w.M[0];
            s = 0;
        }
        /// Compute and extend the wave of score s + 1.
        const Wave& step() {
            const int sc = ++s;
            Wave& w = ring[static_cast<std::size_t>(sc % W)];
            w.s = sc; w.hasI = w.hasD = false; w.lo = 1; w.hi = 0;
            w.amax = INT_MIN;
            // With a forced start gap the score-0 "M" is the virtual origin, from which only a gap of the forced kind may open.
            const bool virt = forced != GapEnd::kNone;
            const Wave* mx = (virt && sc - x == 0) ? nullptr : get(sc - x);
            const Wave* moe = get(sc - oe);
            const Wave* moeI = (virt && sc - oe == 0 && forced != GapEnd::kText) ? nullptr : moe;
            const Wave* moeD = (virt && sc - oe == 0 && forced != GapEnd::kPat) ? nullptr : moe;
            const Wave* ge = get(sc - e);
            const bool gI = ge && ge->hasI, gD = ge && ge->hasD;
            int lo = INT_MAX, hi = INT_MIN;
            if (mx)   { lo = std::min(lo, mx->lo);       hi = std::max(hi, mx->hi); }
            if (moeI) { lo = std::min(lo, moeI->lo + 1); hi = std::max(hi, moeI->hi + 1); }
            if (moeD) { lo = std::min(lo, moeD->lo - 1); hi = std::max(hi, moeD->hi - 1); }
            if (gI)   { lo = std::min(lo, ge->lo + 1);   hi = std::max(hi, ge->hi + 1); }
            if (gD)   { lo = std::min(lo, ge->lo - 1);   hi = std::max(hi, ge->hi - 1); }
            lo = std::max(lo, -m); hi = std::min(hi, n);
            if (lo > hi) return w;
            w.lo = lo; w.hi = hi;
            const int Wd = hi - lo + 1;
            const bool haveI = moeI || gI, haveD = moeD || gD;
            w.M.resize(static_cast<std::size_t>(Wd));
            if (haveI) w.I.resize(static_cast<std::size_t>(Wd));
            if (haveD) w.D.resize(static_cast<std::size_t>(Wd));
            w.hasI = haveI; w.hasD = haveD;
            int* M = w.M.data();
            int* I = haveI ? w.I.data() : nullptr;
            int* D = haveD ? w.D.data() : nullptr;
            // One fused pass per wave, where separate fills, five shifted-max folds and a validity pass took about 7 sweeps of the wave in the profile.
            // I comes from M(s-o-e)[k-1] + 1 and I(s-e)[k-1] + 1, D from M(s-o-e)[k+1] and D(s-e)[k+1], and M from M(s-x)[k] + 1.
            // Each source is read where it exists and reads as kUnreachable elsewhere, offsets past either sequence are invalidated, and M = max(M, I, D).
            struct Src { const int* p = nullptr; int a = 1, b = 0; };   // p[k] valid for k in [a, b]
            auto src = [](const Wave* wv, const std::vector<int>& v, int shift) {
                Src r;
                if (wv && wv->lo <= wv->hi && !v.empty()) {
                    r.p = v.data() - wv->lo + shift; r.a = wv->lo - shift; r.b = wv->hi - shift;
                }
                return r;
            };
            const Src srcs[5] = {src(mx, mx ? mx->M : w.M, 0),
                                 src(moeI, moeI ? moeI->M : w.M, -1), src(gI ? ge : nullptr, gI ? ge->I : w.M, -1),
                                 src(moeD, moeD ? moeD->M : w.M, +1), src(gD ? ge : nullptr, gD ? ge->D : w.M, +1)};
            // missing sources, and sources outside their range, read a kUnreachable row
            thread_local std::vector<int> negrow;
            if (negrow.size() < static_cast<std::size_t>(Wd)) negrow.assign(static_cast<std::size_t>(Wd), kUnreachable);
            const int* nr = negrow.data() - lo;                    // nr[k] valid for k in [lo, hi]
            // Segments between source range ends, inside each of which every source is wholly present or wholly absent, so the kernel has no branches
            int cut[12], nc = 0;
            cut[nc++] = lo; cut[nc++] = hi + 1;
            for (const Src& x : srcs)
                if (x.p) {
                    if (x.a > lo && x.a <= hi) cut[nc++] = x.a;
                    if (x.b + 1 > lo && x.b + 1 <= hi) cut[nc++] = x.b + 1;
                }
            for (int x = 1; x < nc; ++x)                           // at most 12 cuts, so insertion sort
                for (int y = x; y > 0 && cut[y - 1] > cut[y]; --y) std::swap(cut[y - 1], cut[y]);
            nc = static_cast<int>(std::unique(cut, cut + nc) - cut);
            for (int c = 0; c + 1 < nc; ++c) {
                const int a = cut[c], b = cut[c + 1] - 1;
                const int* q[5];
                for (int t = 0; t < 5; ++t)
                    q[t] = (srcs[t].p && srcs[t].a <= a && b <= srcs[t].b) ? srcs[t].p : nr;
                if (haveI && haveD)  kernel<true, true>(q, M - lo, I - lo, D - lo, a, b, m, n);
                else if (haveI)      kernel<true, false>(q, M - lo, I - lo, nullptr, a, b, m, n);
                else if (haveD)      kernel<false, true>(q, M - lo, nullptr, D - lo, a, b, m, n);
                else                 kernel<false, false>(q, M - lo, nullptr, nullptr, a, b, m, n);
            }
            // Extend along matches.
            // M >= I and M >= D, so M bounds the antidiagonal.
            int i0 = 0;
#ifdef __AVX2__
            i0 = extend_x8(M, lo, Wd, w.amax);
#endif
            for (int i = i0; i < Wd; ++i) {
                if (M[i] < 0) continue;
                const int k = lo + i;
                const int h = extend(M[i], k);
                M[i] = h;
                w.amax = std::max(w.amax, 2 * h - k);
            }
            return w;
        }
    };

    /// pad() - WfaAffine::load() for a slice, optionally reversed, giving @len bases of @s upper-cased with others turned into @nch, then kPad bytes of @padc.
    static void pad(const char* s, int len, bool rev, std::vector<char>& out, char nch, char padc) {
        out.resize(static_cast<std::size_t>(len) + WfaAffine::kPad);
        for (int i = 0; i < len; ++i) {
            char c = s[rev ? len - 1 - i : i];
            if (c >= 'a' && c <= 'z') c = static_cast<char>(c - WfaAffine::kAsciiCaseOffset);
            out[static_cast<std::size_t>(i)] = (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? c : nch;
        }
        std::memset(out.data() + len, padc, WfaAffine::kPad);
    }

    /// search() - breakpoint search on P_[p0,p1) x T_[t0,t1).
    /// Returns the optimal score, -1 above max_score, and the join as forward v, h and kind.
    int search(int p0, int p1, int t0, int t1, GapEnd sfree, GapEnd ein, int max_score,
               int& bv, int& bh, GapEnd& bkind) {
        const int m = p1 - p0, n = t1 - t0;
        pad(P_.data() + p0, m, false, fp_, 'N', '!');
        pad(T_.data() + t0, n, false, ft_, 'n', '?');
        pad(P_.data() + p0, m, true, rp_, 'N', '!');
        pad(T_.data() + t0, n, true, rt_, 'n', '?');
        Engine F, R;
        for (Engine* E : {&F, &R}) {
            E->m = m; E->n = n; E->x = pen_.x; E->oe = pen_.o + pen_.e; E->e = pen_.e; E->W = W_;
        }
        F.P = fp_.data(); F.T = ft_.data(); F.start = sfree;
        R.P = rp_.data(); R.T = rt_.data(); R.forced = ein;
        F.init(); R.init();
        int best = INT_MAX;
        const int dk = n - m;
        // Pair one new wave of `A` (A_is_fwd) against every live wave of `B`.
        // A join needs h_f + h_r >= n on one diagonal, i.e. antidiagonals (2h - k) summing to at least m + n, so waves whose maxima fall short cannot meet, which skips the scan until the fronts are close.
        auto check = [&](const Wave& a, bool a_fwd, const Engine& B) {
            if (a.lo > a.hi) return;
            for (int sb = std::max(0, B.s - W_ + 1); sb <= B.s; ++sb) {
                const Wave* b = B.get(sb);
                if (!b) continue;
                if (static_cast<long>(a.amax) + b->amax < static_cast<long>(m) + n) continue;
                const bool skipM = (a_fwd ? (R.forced != GapEnd::kNone && sb == 0)
                                          : (R.forced != GapEnd::kNone && a.s == 0));
                for (int k = a.lo; k <= a.hi; ++k) {
                    const int kb = dk - k;
                    if (kb < b->lo || kb > b->hi) continue;
                    const std::size_t ia = static_cast<std::size_t>(k - a.lo);
                    const std::size_t ib = static_cast<std::size_t>(kb - b->lo);
                    const int kf = a_fwd ? k : kb;              // forward diagonal
                    auto join = [&](int ha, int hb, int cost, GapEnd kind) {
                        if (ha < 0 || hb < 0 || ha + hb < n || cost >= best) return;
                        const int h = a_fwd ? ha : hb;          // the forward offset
                        best = cost; bh = h; bv = h - kf; bkind = kind;
                    };
                    if (!skipM) join(a.M[ia], b->M[ib], a.s + sb, GapEnd::kNone);
                    if (a.hasI && b->hasI)
                        join(a.I[ia], b->I[ib], a.s + sb - pen_.o, GapEnd::kText);
                    if (a.hasD && b->hasD)
                        join(a.D[ia], b->D[ib], a.s + sb - pen_.o, GapEnd::kPat);
                }
            }
        };
        check(F.ring[0], true, R);
        for (;;) {
            // lower bound of any join found from here on
            const int lb = std::min(F.s, R.s) + 1 + std::max(F.s, R.s) - (W_ - 1) - pen_.o;
            if (best != INT_MAX && lb >= best) break;
            if (max_score >= 0 && lb > max_score) return -1;
            if (F.s <= R.s) check(F.step(), true, R);
            else            check(R.step(), false, F);
        }
        if (max_score >= 0 && best > max_score) return -1;
        return best;
    }

    /// Align P_[p0,p1) x T_[t0,t1) under the gap constraints, appending ops.
    int solve(int p0, int p1, int t0, int t1, GapEnd sfree, GapEnd ein, int max_score) {
        // A finite bound, since every alignment costs at most a mismatch per diagonal step plus two gaps, so an infeasible constraint ends the search.
        const int m0 = p1 - p0, n0 = t1 - t0;
        const int bound = pen_.x * std::min(m0, n0) + 2 * (pen_.o + pen_.e * std::max(m0, n0)) + pen_.o;
        if (max_score < 0 || max_score > bound) max_score = bound;
        int bv = 0, bh = 0; GapEnd bk = GapEnd::kNone;
        const int sc = search(p0, p1, t0, t1, sfree, ein, max_score, bv, bh, bk);
        if (sc < 0) return -1;
        const int m = p1 - p0, n = t1 - t0;
        const bool corner = (bv == 0 && bh == 0) || (bv == m && bh == n);
        if (sc <= base_score_ || corner) {
            WfaAffineResult r = base_.align_global(P_.data() + p0, m, T_.data() + t0, n,
                                                   sc, sfree, ein);
            if (r.score < 0) return -1;
            append(r.cigar);
            return r.score;
        }
        if (solve(p0, p0 + bv, t0, t0 + bh, sfree, bk, -1) < 0) return -1;
        if (solve(p0 + bv, p1, t0 + bh, t1, bk, ein, -1) < 0) return -1;
        return sc;
    }

    void append(const std::string& cig) {
        int n = 0;
        for (char c : cig) {
            if (c >= '0' && c <= '9') n = n * 10 + (c - '0');
            else { ops_.append(static_cast<std::size_t>(n), c); n = 0; }
        }
    }

    AffinePenalties pen_;
    WfaAffine base_;
    int W_, base_score_;
    std::vector<char> P_, T_, fp_, ft_, rp_, rt_;
    std::string ops_;
};

}
