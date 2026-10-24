// wfa_affine_stress.cpp - optik's own gap-affine WFA against an independent Gotoh DP.
//
//   g++ -std=c++17 -O2 -march=native -I include bench/cpp/wfa_affine_stress.cpp && ./a.out
//
// Random pattern/text pairs with substitutions, indels, N bases and varied ends-free ranges and penalties are checked for these properties.
//   * score == Gotoh DP optimum under identical ends-free semantics
//   * the CIGAR consumes the whole pattern, spans [text_begin, text_end) of the text within the free ranges, and re-scores to exactly `score`
//   * the max_score cut-off returns -1 exactly when the optimum exceeds it
//   * end-to-end alignment starting inside an open gap == a constrained Gotoh DP
//   * BiWFA with O(s) memory == Gotoh optimum with a valid CIGAR, at recursion thresholds 0 / 8 / 512, and == the traceback WFA on 2-20 kb pairs.
#include <algorithm>
#include <climits>
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include "wfa_affine.hpp"
using namespace tpoptoa;

static int gotoh(const std::string& p, const std::string& t, AffinePenalties pen, int tbf, int tef) {
    const int m = p.size(), n = t.size(), INF = INT_MAX / 4;
    auto eq = [](char a, char b) { return a == b && a != 'N'; };
    std::vector<int> M(n + 1), X(n + 1), Y(n + 1), pM(n + 1), pX(n + 1), pY(n + 1);
    // X is a gap in the pattern, consuming text, and Y a gap in the text, consuming pattern
    for (int j = 0; j <= n; ++j) { pM[j] = j <= tbf ? 0 : INF; pX[j] = j <= tbf ? INF : INF; pY[j] = INF; }
    // allow leading text gap beyond tbf as a real gap from column tbf
    for (int j = 1; j <= n; ++j) if (j > tbf) { pX[j] = std::min(pM[j - 1] + pen.o + pen.e, pX[j - 1] + pen.e); pM[j] = pX[j]; }
    for (int i = 1; i <= m; ++i) {
        M[0] = INF; X[0] = INF; Y[0] = std::min(pM[0] + pen.o + pen.e, pY[0] + pen.e);
        M[0] = Y[0];
        for (int j = 1; j <= n; ++j) {
            X[j] = std::min(M[j - 1] + pen.o + pen.e, X[j - 1] + pen.e);
            Y[j] = std::min(pM[j] + pen.o + pen.e, pY[j] + pen.e);
            const int d = pM[j - 1] + (eq(p[i - 1], t[j - 1]) ? 0 : pen.x);
            M[j] = std::min(d, std::min(X[j], Y[j]));
        }
        std::swap(M, pM); std::swap(X, pX); std::swap(Y, pY);
    }
    int best = INF;
    for (int j = std::max(0, n - tef); j <= n; ++j) best = std::min(best, pM[j]);
    return best;
}

static int rescore(const std::string& cig, const std::string& p, const std::string& t, int tb,
                   AffinePenalties pen, int& pconsumed, int& tend) {
    int i = 0, j = tb, sc = 0; std::size_t q = 0;
    while (q < cig.size()) {
        int len = 0; while (isdigit(cig[q])) len = len * 10 + (cig[q++] - '0');
        const char op = cig[q++];
        if (op == 'M') { for (int r = 0; r < len; ++r, ++i, ++j) if (p[i] != t[j] || p[i] == 'N') return -9; }
        else if (op == 'X') { for (int r = 0; r < len; ++r, ++i, ++j) sc += pen.x; }
        else if (op == 'I') { sc += pen.o + pen.e * len; i += len; }
        else if (op == 'D') { sc += pen.o + pen.e * len; j += len; }
    }
    pconsumed = i; tend = j;
    return sc;
}

/// Independent O(mn) semi-global unit-cost edit distance with the pattern end to end and the text ends free.
/// A pattern 'N' matches nothing.
static int naive_semi_global_ed(const std::string& p, const std::string& t) {
    std::vector<int> prev(t.size() + 1, 0), cur(t.size() + 1);
    for (std::size_t i = 1; i <= p.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (std::size_t j = 1; j <= t.size(); ++j) {
            const int sub = prev[j - 1] + ((p[i - 1] != t[j - 1] || p[i - 1] == 'N') ? 1 : 0);
            cur[j] = std::min({sub, prev[j] + 1, cur[j - 1] + 1});
        }
        std::swap(prev, cur);
    }
    return *std::min_element(prev.begin(), prev.end());
}

/// Global Gotoh with gap-state constraints, where start_free means it may begin inside an open gap of that kind at no open cost, 1 consuming text and 2 consuming pattern.
/// end_in requires ending with a gap of that kind.
/// Full matrices.
static int gotoh_c(const std::string& p, const std::string& t, AffinePenalties pen,
                   int start_free, int end_in) {
    const int m = p.size(), n = t.size(), INF = INT_MAX / 4;
    auto eq = [](char a, char b) { return a == b && a != 'N'; };
    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, INF)), X = M, Y = M;
    M[0][0] = 0;
    if (start_free == 1) X[0][0] = 0;
    if (start_free == 2) Y[0][0] = 0;
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= n; ++j) {
            if (j > 0) X[i][j] = std::min(X[i][j], std::min(M[i][j - 1] + pen.o + pen.e, X[i][j - 1] + pen.e));
            if (i > 0) Y[i][j] = std::min(Y[i][j], std::min(M[i - 1][j] + pen.o + pen.e, Y[i - 1][j] + pen.e));
            if (i > 0 && j > 0) M[i][j] = std::min(M[i][j], M[i - 1][j - 1] + (eq(p[i - 1], t[j - 1]) ? 0 : pen.x));
            if (i || j) M[i][j] = std::min(M[i][j], std::min(X[i][j], Y[i][j]));
        }
    const int v = end_in == 1 ? X[m][n] : end_in == 2 ? Y[m][n] : M[m][n];
    return v >= INF ? -1 : v;
}

static std::string mutate(const std::string& p, double er, std::mt19937_64& rng) {
    const char* A = "ACGT";
    std::string core;
    for (char c : p) {
        double u = std::uniform_real_distribution<double>(0, 1)(rng);
        if (u < er / 3) continue;
        if (u < 2 * er / 3) { core.push_back(A[rng() % 4]); core.push_back(c); continue; }
        if (u < er) { core.push_back(A[rng() % 4]); continue; }
        core.push_back(c);
    }
    return core;
}

int main() {
    std::mt19937_64 rng(11);
    int fail = 0; long cases = 0;
    const char* A = "ACGT";
    for (int it = 0; it < 20000; ++it) {
        AffinePenalties pen;
        pen.x = 1 + rng() % 6; pen.o = rng() % 8; pen.e = 1 + rng() % 4;
        const int m = 1 + rng() % 120;
        std::string p; for (int i = 0; i < m; ++i) p.push_back(rng() % 50 == 0 ? 'N' : A[rng() % 4]);
        // text = flank + mutated pattern + flank
        std::string core;
        const double er = (rng() % 4) * 0.05;
        for (char c : p) {
            double u = std::uniform_real_distribution<double>(0, 1)(rng);
            if (u < er / 3) continue;
            if (u < 2 * er / 3) { core.push_back(A[rng() % 4]); core.push_back(c); continue; }
            if (u < er) { core.push_back(A[rng() % 4]); continue; }
            core.push_back(c);
        }
        std::string l, r;
        const int fl = rng() % 30, fr = rng() % 30;
        for (int i = 0; i < fl; ++i) l.push_back(A[rng() % 4]);
        for (int i = 0; i < fr; ++i) r.push_back(A[rng() % 4]);
        const std::string t = (rng() % 8 == 0) ? std::string() + A[rng() % 4] : l + core + r;
        const int n = t.size();
        const int tbf = (rng() % 3 == 0) ? 0 : rng() % (n + 1);
        const int tef = (rng() % 3 == 0) ? 0 : rng() % (n + 1);
        const int want = gotoh(p, t, pen, tbf, tef);
        WfaAffine w(pen);
        auto res = w.align(p.data(), m, t.data(), n, tbf, tef);
        ++cases;
        int pc = 0, te = 0;
        const int rs = rescore(res.cigar, p, t, res.text_begin, pen, pc, te);
        bool ok = res.score == want && rs == want && pc == m && te == res.text_end
                  && res.text_begin <= tbf && res.text_end >= n - tef && res.text_end <= n;
        // max_score cut-off
        if (want > 0) {
            auto cut = w.align(p.data(), m, t.data(), n, tbf, tef, want - 1, false);
            ok = ok && cut.score == -1;
            auto at = w.align(p.data(), m, t.data(), n, tbf, tef, want, false);
            ok = ok && at.score == want;
        }
        if (!ok && ++fail <= 5)
            std::printf("FAIL m=%d n=%d x=%d o=%d e=%d tbf=%d tef=%d want=%d got=%d rescore=%d pc=%d te=%d/%d cig=%s\n",
                        m, n, pen.x, pen.o, pen.e, tbf, tef, want, res.score, rs, pc, te, res.text_end, res.cigar.c_str());
    }
    // edit mode (x=1, o=0, e=1, text ends fully free, capped) == the scalar semi-global edit-distance DP used by the mapper before
    {
        WfaAffine ed(AffinePenalties{1, 0, 1});
        int bad = 0;
        for (int it = 0; it < 20000; ++it) {
            const int m = 1 + rng() % 150;
            std::string p; for (int i = 0; i < m; ++i) p.push_back(rng() % 60 == 0 ? 'N' : A[rng() % 4]);
            std::string t; const int n = rng() % 260;
            for (int i = 0; i < n; ++i) t.push_back(A[rng() % 4]);
            if (n > m && rng() % 2) t.replace(rng() % (n - m + 1), m, p);   // plant the read
            const int cap = rng() % 20;
            const int want = std::min(naive_semi_global_ed(p, t), cap + 1);
            auto r = ed.align(p.data(), m, t.data(), n, n, n, cap, false);
            const int got = r.score < 0 ? cap + 1 : r.score;
            if (got != want && ++bad <= 3) std::printf("EDIT FAIL m=%d n=%d cap=%d want=%d got=%d\n", m, n, cap, want, got);
        }
        fail += bad; cases += 20000;
    }
    // A free start gap with an unconstrained end equals constrained Gotoh.
    // Forced end gaps are exact only where a breakpoint search proved the state reachable, since furthest-reaching dominance fails at the sequence boundary for an arbitrary forced end.
    // The BiWFA cases below cover them, as base 0 routes every half through them.
    {
        int bad = 0;
        for (int it = 0; it < 20000; ++it) {
            AffinePenalties pen;
            pen.x = 1 + rng() % 6; pen.o = rng() % 8; pen.e = 1 + rng() % 4;
            const int m = rng() % 60;
            std::string p; for (int i = 0; i < m; ++i) p.push_back(rng() % 40 == 0 ? 'N' : A[rng() % 4]);
            std::string t = mutate(p, (rng() % 5) * 0.08, rng);
            if (rng() % 6 == 0) t.clear();
            const int sf = rng() % 3, ei = 0;
            const int want = gotoh_c(p, t, pen, sf, ei);
            WfaAffine w(pen);
            const int cap = 4 * (m + (int)t.size()) * (pen.x + pen.o + pen.e) + 8;
            auto r = w.align_global(p.data(), m, t.data(), t.size(), cap,
                                    static_cast<GapEnd>(sf), static_cast<GapEnd>(ei));
            bool ok = r.score == want;
            if (ok && want >= 0) {
                int pc = 0, te = 0;
                int rs = rescore(r.cigar, p, t, 0, pen, pc, te);
                // A free start gap is not charged its open, so it is credited back.
                if (sf && !r.cigar.empty()) {
                    std::size_t q = 0; while (isdigit(r.cigar[q])) ++q;
                    const char op0 = r.cigar[q];
                    if ((sf == 1 && op0 == 'D') || (sf == 2 && op0 == 'I')) rs -= pen.o;
                }
                const char last = r.cigar.empty() ? 0 : r.cigar.back();
                ok = rs == want && pc == m && te == (int)t.size()
                     && (ei == 0 || (ei == 1 && last == 'D') || (ei == 2 && last == 'I'));
            }
            if (!ok && ++bad <= 3)
                std::printf("CONSTRAINED FAIL m=%d n=%zu x=%d o=%d e=%d sf=%d ei=%d want=%d got=%d cig=%s p=%s t=%s\n",
                            m, t.size(), pen.x, pen.o, pen.e, sf, ei, want, r.score, r.cigar.c_str(), p.c_str(), t.c_str());
        }
        fail += bad; cases += 20000;
    }
    // BiWFA == Gotoh optimum, with small thresholds forcing deep recursion
    {
        int bad = 0;
        for (int it = 0; it < 30000; ++it) {
            AffinePenalties pen;
            pen.x = 1 + rng() % 6; pen.o = rng() % 8; pen.e = 1 + rng() % 4;
            const int m = rng() % 160;
            std::string p; for (int i = 0; i < m; ++i) p.push_back(rng() % 50 == 0 ? 'N' : A[rng() % 4]);
            std::string t = mutate(p, (rng() % 6) * 0.06, rng);
            if (rng() % 10 == 0) { t.clear(); for (int i = 0, L = rng() % 80; i < L; ++i) t.push_back(A[rng() % 4]); }
            const int want = gotoh(p, t, pen, 0, 0);
            const int base = (int[]){0, 8, 512}[it % 3];
            BiWfa b(pen, base);
            auto r = b.align(p.data(), m, t.data(), t.size());
            int pc = 0, te = 0;
            const int rs = r.score >= 0 ? rescore(r.cigar, p, t, 0, pen, pc, te) : -1;
            const bool ok = r.score == want && rs == want && pc == m && te == (int)t.size();
            if (!ok && ++bad <= 3)
                std::printf("BIWFA FAIL base=%d m=%d n=%zu x=%d o=%d e=%d want=%d got=%d rescore=%d pc=%d te=%d cig=%s\n",
                            base, m, t.size(), pen.x, pen.o, pen.e, want, r.score, rs, pc, te, r.cigar.c_str());
            // score cap
            if (ok && want > 0) {
                BiWfa c(pen, base);
                if (c.align(p.data(), m, t.data(), t.size(), want - 1).score != -1 && ++bad <= 3)
                    std::printf("BIWFA CAP FAIL m=%d want=%d\n", m, want);
            }
        }
        fail += bad; cases += 30000;
    }
    // On large pairs BiWFA must equal traceback WFA with a valid CIGAR.
    {
        int bad = 0;
        for (int it = 0; it < 40; ++it) {
            AffinePenalties pen;                       // mapper penalties
            const int m = 2000 + rng() % 18000;
            std::string p; for (int i = 0; i < m; ++i) p.push_back(A[rng() % 4]);
            std::string t = mutate(p, 0.02 + (rng() % 15) * 0.02, rng);
            WfaAffine w(pen);
            const int want = w.align(p.data(), m, t.data(), t.size(), 0, 0, -1, false).score;
            BiWfa b(pen);
            auto r = b.align(p.data(), m, t.data(), t.size());
            int pc = 0, te = 0;
            const int rs = rescore(r.cigar, p, t, 0, pen, pc, te);
            if ((r.score != want || rs != want || pc != m || te != (int)t.size()) && ++bad <= 3)
                std::printf("BIWFA LARGE FAIL m=%d n=%zu want=%d got=%d rescore=%d\n", m, t.size(), want, r.score, rs);
        }
        fail += bad; cases += 40;
    }
    if (fail) std::printf("wfa_affine_stress: %d/%ld FAILURE(S)\n", fail, cases);
    else std::printf("wfa_affine_stress: ALL PASS (%ld cases)\n", cases);
    return fail != 0;
}
