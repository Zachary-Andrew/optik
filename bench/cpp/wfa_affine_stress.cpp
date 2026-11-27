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

static int gotoh(const std::string& pattern, const std::string& text, AffinePenalties pen, int tbf, int tef) {
    const int pattern_len = pattern.size(), text_len = text.size(), INF = INT_MAX / 4;
    auto same = [](char lhs, char rhs) { return lhs == rhs && lhs != 'N'; };
    std::vector<int> match_row(text_len + 1), gap_text_row(text_len + 1), gap_pattern_row(text_len + 1), prev_match(text_len + 1), prev_gap_text(text_len + 1), prev_gap_pattern(text_len + 1);
    // X is a gap in the pattern, consuming text, and Y a gap in the text, consuming pattern
    for (int j = 0; j <= text_len; ++j) { prev_match[j] = j <= tbf ? 0 : INF; prev_gap_text[j] = j <= tbf ? INF : INF; prev_gap_pattern[j] = INF; }
    // allow leading text gap beyond tbf as a real gap from column tbf
    for (int j = 1; j <= text_len; ++j) if (j > tbf) { prev_gap_text[j] = std::min(prev_match[j - 1] + pen.gap_open + pen.gap_extend, prev_gap_text[j - 1] + pen.gap_extend); prev_match[j] = prev_gap_text[j]; }
    for (int i = 1; i <= pattern_len; ++i) {
        match_row[0] = INF; gap_text_row[0] = INF; gap_pattern_row[0] = std::min(prev_match[0] + pen.gap_open + pen.gap_extend, prev_gap_pattern[0] + pen.gap_extend);
        match_row[0] = gap_pattern_row[0];
        for (int j = 1; j <= text_len; ++j) {
            gap_text_row[j] = std::min(match_row[j - 1] + pen.gap_open + pen.gap_extend, gap_text_row[j - 1] + pen.gap_extend);
            gap_pattern_row[j] = std::min(prev_match[j] + pen.gap_open + pen.gap_extend, prev_gap_pattern[j] + pen.gap_extend);
            const int diag = prev_match[j - 1] + (same(pattern[i - 1], text[j - 1]) ? 0 : pen.mismatch);
            match_row[j] = std::min(diag, std::min(gap_text_row[j], gap_pattern_row[j]));
        }
        std::swap(match_row, prev_match); std::swap(gap_text_row, prev_gap_text); std::swap(gap_pattern_row, prev_gap_pattern);
    }
    int best = INF;
    for (int j = std::max(0, text_len - tef); j <= text_len; ++j) best = std::min(best, prev_match[j]);
    return best;
}

static int rescore(const std::string& cig, const std::string& pattern, const std::string& text, int text_begin,
                   AffinePenalties pen, int& pconsumed, int& tend) {
    int pat_pos = 0, text_pos = text_begin, score = 0; std::size_t cig_pos = 0;
    while (cig_pos < cig.size()) {
        int len = 0; while (isdigit(cig[cig_pos])) len = len * 10 + (cig[cig_pos++] - '0');
        const char cig_op = cig[cig_pos++];
        if (cig_op == 'M') { for (int r = 0; r < len; ++r, ++pat_pos, ++text_pos) if (pattern[pat_pos] != text[text_pos] || pattern[pat_pos] == 'N') return -9; }
        else if (cig_op == 'X') { for (int r = 0; r < len; ++r, ++pat_pos, ++text_pos) score += pen.mismatch; }
        else if (cig_op == 'I') { score += pen.gap_open + pen.gap_extend * len; pat_pos += len; }
        else if (cig_op == 'D') { score += pen.gap_open + pen.gap_extend * len; text_pos += len; }
    }
    pconsumed = pat_pos; tend = text_pos;
    return score;
}

/// Independent O(mn) semi-global unit-cost edit distance with the pattern end to end and the text ends free.
/// A pattern 'N' matches nothing.
static int naive_semi_global_ed(const std::string& pattern, const std::string& text) {
    std::vector<int> prev(text.size() + 1, 0), cur(text.size() + 1);
    for (std::size_t i = 1; i <= pattern.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (std::size_t j = 1; j <= text.size(); ++j) {
            const int sub = prev[j - 1] + ((pattern[i - 1] != text[j - 1] || pattern[i - 1] == 'N') ? 1 : 0);
            cur[j] = std::min({sub, prev[j] + 1, cur[j - 1] + 1});
        }
        std::swap(prev, cur);
    }
    return *std::min_element(prev.begin(), prev.end());
}

/// Global Gotoh with gap-state constraints, where start_free means it may begin inside an open gap of that kind at no open cost, 1 consuming text and 2 consuming pattern.
/// end_in requires ending with a gap of that kind.
/// Full matrices.
static int gotoh_c(const std::string& pattern, const std::string& text, AffinePenalties pen,
                   int start_free, int end_in) {
    const int pattern_len = pattern.size(), text_len = text.size(), INF = INT_MAX / 4;
    auto same = [](char lhs, char rhs) { return lhs == rhs && lhs != 'N'; };
    std::vector<std::vector<int>> match_dp(pattern_len + 1, std::vector<int>(text_len + 1, INF)), gap_text_dp = match_dp, gap_pattern_dp = match_dp;
    match_dp[0][0] = 0;
    if (start_free == 1) gap_text_dp[0][0] = 0;
    if (start_free == 2) gap_pattern_dp[0][0] = 0;
    for (int i = 0; i <= pattern_len; ++i)
        for (int j = 0; j <= text_len; ++j) {
            if (j > 0) gap_text_dp[i][j] = std::min(gap_text_dp[i][j], std::min(match_dp[i][j - 1] + pen.gap_open + pen.gap_extend, gap_text_dp[i][j - 1] + pen.gap_extend));
            if (i > 0) gap_pattern_dp[i][j] = std::min(gap_pattern_dp[i][j], std::min(match_dp[i - 1][j] + pen.gap_open + pen.gap_extend, gap_pattern_dp[i - 1][j] + pen.gap_extend));
            if (i > 0 && j > 0) match_dp[i][j] = std::min(match_dp[i][j], match_dp[i - 1][j - 1] + (same(pattern[i - 1], text[j - 1]) ? 0 : pen.mismatch));
            if (i || j) match_dp[i][j] = std::min(match_dp[i][j], std::min(gap_text_dp[i][j], gap_pattern_dp[i][j]));
        }
    const int result = end_in == 1 ? gap_text_dp[pattern_len][text_len] : end_in == 2 ? gap_pattern_dp[pattern_len][text_len] : match_dp[pattern_len][text_len];
    return result >= INF ? -1 : result;
}

static std::string mutate(const std::string& pattern, double err_rate, std::mt19937_64& rng) {
    const char* bases = "ACGT";
    std::string core;
    for (char ch : pattern) {
        double roll = std::uniform_real_distribution<double>(0, 1)(rng);
        if (roll < err_rate / 3) continue;
        if (roll < 2 * err_rate / 3) { core.push_back(bases[rng() % 4]); core.push_back(ch); continue; }
        if (roll < err_rate) { core.push_back(bases[rng() % 4]); continue; }
        core.push_back(ch);
    }
    return core;
}

int main() {
    std::mt19937_64 rng(11);
    int fail = 0; long cases = 0;
    const char* bases = "ACGT";
    for (int it = 0; it < 20000; ++it) {
        AffinePenalties pen;
        pen.mismatch = 1 + rng() % 6; pen.gap_open = rng() % 8; pen.gap_extend = 1 + rng() % 4;
        const int pattern_len = 1 + rng() % 120;
        std::string pattern; for (int i = 0; i < pattern_len; ++i) pattern.push_back(rng() % 50 == 0 ? 'N' : bases[rng() % 4]);
        // text = flank + mutated pattern + flank
        std::string core;
        const double err_rate = (rng() % 4) * 0.05;
        for (char ch : pattern) {
            double roll = std::uniform_real_distribution<double>(0, 1)(rng);
            if (roll < err_rate / 3) continue;
            if (roll < 2 * err_rate / 3) { core.push_back(bases[rng() % 4]); core.push_back(ch); continue; }
            if (roll < err_rate) { core.push_back(bases[rng() % 4]); continue; }
            core.push_back(ch);
        }
        std::string left_flank, right_flank;
        const int left_len = rng() % 30, right_len = rng() % 30;
        for (int i = 0; i < left_len; ++i) left_flank.push_back(bases[rng() % 4]);
        for (int i = 0; i < right_len; ++i) right_flank.push_back(bases[rng() % 4]);
        const std::string text = (rng() % 8 == 0) ? std::string() + bases[rng() % 4] : left_flank + core + right_flank;
        const int text_len = text.size();
        const int tbf = (rng() % 3 == 0) ? 0 : rng() % (text_len + 1);
        const int tef = (rng() % 3 == 0) ? 0 : rng() % (text_len + 1);
        const int want = gotoh(pattern, text, pen, tbf, tef);
        WfaAffine aligner(pen);
        auto res = aligner.align(pattern.data(), pattern_len, text.data(), text_len, tbf, tef);
        ++cases;
        int pat_consumed = 0, text_end = 0;
        const int rescored = rescore(res.cigar, pattern, text, res.text_begin, pen, pat_consumed, text_end);
        bool passed = res.score == want && rescored == want && pat_consumed == pattern_len && text_end == res.text_end
                  && res.text_begin <= tbf && res.text_end >= text_len - tef && res.text_end <= text_len;
        // max_score cut-off
        if (want > 0) {
            auto cut = aligner.align(pattern.data(), pattern_len, text.data(), text_len, tbf, tef, want - 1, false);
            passed = passed && cut.score == -1;
            auto at_cap = aligner.align(pattern.data(), pattern_len, text.data(), text_len, tbf, tef, want, false);
            passed = passed && at_cap.score == want;
        }
        if (!passed && ++fail <= 5)
            std::printf("FAIL m=%d n=%d x=%d o=%d e=%d tbf=%d tef=%d want=%d got=%d rescore=%d pc=%d te=%d/%d cig=%s\n",
                        pattern_len, text_len, pen.mismatch, pen.gap_open, pen.gap_extend, tbf, tef, want, res.score, rescored, pat_consumed, text_end, res.text_end, res.cigar.c_str());
    }
    // edit mode (x=1, o=0, e=1, text ends fully free, capped) == the scalar semi-global edit-distance DP used by the mapper before
    {
        WfaAffine edit_aligner(AffinePenalties{1, 0, 1});
        int bad = 0;
        for (int it = 0; it < 20000; ++it) {
            const int pattern_len = 1 + rng() % 150;
            std::string pattern; for (int i = 0; i < pattern_len; ++i) pattern.push_back(rng() % 60 == 0 ? 'N' : bases[rng() % 4]);
            std::string text; const int text_len = rng() % 260;
            for (int i = 0; i < text_len; ++i) text.push_back(bases[rng() % 4]);
            if (text_len > pattern_len && rng() % 2) text.replace(rng() % (text_len - pattern_len + 1), pattern_len, pattern);   // plant the read
            const int cap = rng() % 20;
            const int want = std::min(naive_semi_global_ed(pattern, text), cap + 1);
            auto result = edit_aligner.align(pattern.data(), pattern_len, text.data(), text_len, text_len, text_len, cap, false);
            const int got = result.score < 0 ? cap + 1 : result.score;
            if (got != want && ++bad <= 3) std::printf("EDIT FAIL m=%d n=%d cap=%d want=%d got=%d\n", pattern_len, text_len, cap, want, got);
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
            pen.mismatch = 1 + rng() % 6; pen.gap_open = rng() % 8; pen.gap_extend = 1 + rng() % 4;
            const int pattern_len = rng() % 60;
            std::string pattern; for (int i = 0; i < pattern_len; ++i) pattern.push_back(rng() % 40 == 0 ? 'N' : bases[rng() % 4]);
            std::string text = mutate(pattern, (rng() % 5) * 0.08, rng);
            if (rng() % 6 == 0) text.clear();
            const int start_free = rng() % 3, end_in = 0;
            const int want = gotoh_c(pattern, text, pen, start_free, end_in);
            WfaAffine aligner(pen);
            const int cap = 4 * (pattern_len + (int)text.size()) * (pen.mismatch + pen.gap_open + pen.gap_extend) + 8;
            auto result = aligner.align_global(pattern.data(), pattern_len, text.data(), text.size(), cap,
                                    static_cast<GapEnd>(start_free), static_cast<GapEnd>(end_in));
            bool passed = result.score == want;
            if (passed && want >= 0) {
                int pat_consumed = 0, text_end = 0;
                int rescored = rescore(result.cigar, pattern, text, 0, pen, pat_consumed, text_end);
                // A free start gap is not charged its open, so it is credited back.
                if (start_free && !result.cigar.empty()) {
                    std::size_t cig_pos = 0; while (isdigit(result.cigar[cig_pos])) ++cig_pos;
                    const char op0 = result.cigar[cig_pos];
                    if ((start_free == 1 && op0 == 'D') || (start_free == 2 && op0 == 'I')) rescored -= pen.gap_open;
                }
                const char last = result.cigar.empty() ? 0 : result.cigar.back();
                passed = rescored == want && pat_consumed == pattern_len && text_end == (int)text.size()
                     && (end_in == 0 || (end_in == 1 && last == 'D') || (end_in == 2 && last == 'I'));
            }
            if (!passed && ++bad <= 3)
                std::printf("CONSTRAINED FAIL m=%d n=%zu x=%d o=%d e=%d sf=%d ei=%d want=%d got=%d cig=%s p=%s t=%s\n",
                            pattern_len, text.size(), pen.mismatch, pen.gap_open, pen.gap_extend, start_free, end_in, want, result.score, result.cigar.c_str(), pattern.c_str(), text.c_str());
        }
        fail += bad; cases += 20000;
    }
    // BiWFA == Gotoh optimum, with small thresholds forcing deep recursion
    {
        int bad = 0;
        for (int it = 0; it < 30000; ++it) {
            AffinePenalties pen;
            pen.mismatch = 1 + rng() % 6; pen.gap_open = rng() % 8; pen.gap_extend = 1 + rng() % 4;
            const int pattern_len = rng() % 160;
            std::string pattern; for (int i = 0; i < pattern_len; ++i) pattern.push_back(rng() % 50 == 0 ? 'N' : bases[rng() % 4]);
            std::string text = mutate(pattern, (rng() % 6) * 0.06, rng);
            if (rng() % 10 == 0) { text.clear(); for (int i = 0, rand_len = rng() % 80; i < rand_len; ++i) text.push_back(bases[rng() % 4]); }
            const int want = gotoh(pattern, text, pen, 0, 0);
            const int base = (int[]){0, 8, 512}[it % 3];
            BiWfa biwfa(pen, base);
            auto result = biwfa.align(pattern.data(), pattern_len, text.data(), text.size());
            int pat_consumed = 0, text_end = 0;
            const int rescored = result.score >= 0 ? rescore(result.cigar, pattern, text, 0, pen, pat_consumed, text_end) : -1;
            const bool passed = result.score == want && rescored == want && pat_consumed == pattern_len && text_end == (int)text.size();
            if (!passed && ++bad <= 3)
                std::printf("BIWFA FAIL base=%d m=%d n=%zu x=%d o=%d e=%d want=%d got=%d rescore=%d pc=%d te=%d cig=%s\n",
                            base, pattern_len, text.size(), pen.mismatch, pen.gap_open, pen.gap_extend, want, result.score, rescored, pat_consumed, text_end, result.cigar.c_str());
            // score cap
            if (passed && want > 0) {
                BiWfa capped(pen, base);
                if (capped.align(pattern.data(), pattern_len, text.data(), text.size(), want - 1).score != -1 && ++bad <= 3)
                    std::printf("BIWFA CAP FAIL m=%d want=%d\n", pattern_len, want);
            }
        }
        fail += bad; cases += 30000;
    }
    // On large pairs BiWFA must equal traceback WFA with a valid CIGAR.
    {
        int bad = 0;
        for (int it = 0; it < 40; ++it) {
            AffinePenalties pen;                       // mapper penalties
            const int pattern_len = 2000 + rng() % 18000;
            std::string pattern; for (int i = 0; i < pattern_len; ++i) pattern.push_back(bases[rng() % 4]);
            std::string text = mutate(pattern, 0.02 + (rng() % 15) * 0.02, rng);
            WfaAffine aligner(pen);
            const int want = aligner.align(pattern.data(), pattern_len, text.data(), text.size(), 0, 0, -1, false).score;
            BiWfa biwfa(pen);
            auto result = biwfa.align(pattern.data(), pattern_len, text.data(), text.size());
            int pat_consumed = 0, text_end = 0;
            const int rescored = rescore(result.cigar, pattern, text, 0, pen, pat_consumed, text_end);
            if ((result.score != want || rescored != want || pat_consumed != pattern_len || text_end != (int)text.size()) && ++bad <= 3)
                std::printf("BIWFA LARGE FAIL m=%d n=%zu want=%d got=%d rescore=%d\n", pattern_len, text.size(), want, result.score, rescored);
        }
        fail += bad; cases += 40;
    }
    if (fail) std::printf("wfa_affine_stress: %d/%ld FAILURE(S)\n", fail, cases);
    else std::printf("wfa_affine_stress: ALL PASS (%ld cases)\n", cases);
    return fail != 0;
}
