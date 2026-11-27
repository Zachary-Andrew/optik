// wfa_affine_bench.cpp - exact gap-affine pairwise alignment with optik's own wavefront aligners against WFA2-lib, used only as an external reference with heuristics off.
//
// Methods, all exact, end-to-end and with identical penalties
//   optik-wfa      WfaAffine, traceback WFA, O(s^2) memory (Marco-Sola 2021)
//   optik-biwfa    BiWfa, bidirectional WFA, O(s) memory (Marco-Sola 2023)
//   wfa2-high      WFA2-lib MemoryHigh     (classic WFA)
//   wfa2-ultralow  WFA2-lib MemoryUltralow (BiWFA)
// Every alignment's CIGAR is re-scored against both sequences, and a CIGAR that does not consume them exactly or re-scores differently counts as invalid.
//
// Usage:
//   wfa_affine_bench                          self-check on built-in pair sets, all methods, scores must agree
//   wfa_affine_bench METHOD PAIRS.tsv OUT     one method in its own process, so
//                                             peak RSS is that method's.  PAIRS
//                                             has "pattern<TAB>text" lines, OUT
//                                             gets one score per line and stdout
//                                             "align_seconds S invalid N pairs P"
// Build (WFA2-lib from temp/):
//   g++ -std=c++17 -O3 -march=native -I include -I temp/WFA2-lib
//     bench/cpp/wfa_affine_bench.cpp temp/WFA2-lib/bindings/cpp/WFAligner.cpp
//     temp/WFA2-lib/lib/libwfa.a -fopenmp -o wfab
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include "wfa_affine.hpp"

/// rescore() result for a CIGAR that does not consume exactly the pattern and the text.
constexpr int kInvalidCigar = -1;
/// Each simulated pair is seeded with len * kSeedLenMul + err * kSeedErrScale, and its edits are deletions, insertions or substitutions with equal odds.
constexpr int kSeedLenMul = 7;
constexpr double kSeedErrScale = 100;
constexpr double kEditKinds = 3;
/// Argument count of the file mode (method, pair file, score output) and the index of the output path.
constexpr int kFileModeArgc = 4;
constexpr int kOutPathArg = 3;
/// Microseconds per second.
constexpr double kUsPerSec = 1e6;
#include "bindings/cpp/WFAligner.hpp"
using namespace tpoptoa;

static constexpr AffinePenalties PEN{4, 6, 2};

/// Re-score a CIGAR whose ops are M match, X mismatch, I pattern-only and D text-only, with WFA2-lib's '=' read as M.
/// Returns -1 if the CIGAR does not consume both sequences exactly or an M covers a mismatch.
static int rescore(const std::string& cig, const std::string& pattern, const std::string& text) {
    std::size_t pat_pos = 0, text_pos = 0, cig_pos = 0; int score = 0;
    while (cig_pos < cig.size()) {
        long run_len = 0;
        if (!isdigit(static_cast<unsigned char>(cig[cig_pos]))) run_len = 1;     // WFA2 run-less op
        while (cig_pos < cig.size() && isdigit(static_cast<unsigned char>(cig[cig_pos]))) run_len = run_len * 10 + (cig[cig_pos++] - '0');
        const char cig_op = cig[cig_pos++];
        if (cig_op == 'M' || cig_op == '=') {
            for (long r = 0; r < run_len; ++r, ++pat_pos, ++text_pos)
                if (pat_pos >= pattern.size() || text_pos >= text.size() || pattern[pat_pos] != text[text_pos]) return kInvalidCigar;
        } else if (cig_op == 'X') {
            for (long r = 0; r < run_len; ++r, ++pat_pos, ++text_pos) {
                if (pat_pos >= pattern.size() || text_pos >= text.size()) return kInvalidCigar;
                score += PEN.mismatch;
            }
        } else if (cig_op == 'I') { score += PEN.gap_open + PEN.gap_extend * run_len; pat_pos += run_len; }
        else if (cig_op == 'D') { score += PEN.gap_open + PEN.gap_extend * run_len; text_pos += run_len; }
        else return kInvalidCigar;
    }
    return (pat_pos == pattern.size() && text_pos == text.size()) ? score : kInvalidCigar;
}

struct Aligner {
    virtual ~Aligner() = default;
    /// score and CIGAR (M/X/I/D)
    virtual int align(const std::string& pattern, const std::string& text, std::string& cig) = 0;
};
struct OptikWfa : Aligner {
    WfaAffine aligner{PEN};
    int align(const std::string& pattern, const std::string& text, std::string& cig) override {
        auto result = aligner.align(pattern.data(), pattern.size(), text.data(), text.size(), 0, 0);
        cig = result.cigar; return result.score;
    }
};
struct OptikBiWfa : Aligner {
    BiWfa aligner{PEN};
    int align(const std::string& pattern, const std::string& text, std::string& cig) override {
        auto result = aligner.align(pattern.data(), pattern.size(), text.data(), text.size());
        cig = result.cigar; return result.score;
    }
};
struct Wfa2 : Aligner {
    wfa::WFAlignerGapAffine aligner;
    explicit Wfa2(wfa::WFAligner::MemoryModel mem_model)
        : aligner(PEN.mismatch, PEN.gap_open, PEN.gap_extend, wfa::WFAligner::Alignment, mem_model) { aligner.setHeuristicNone(); }
    int align(const std::string& pattern, const std::string& text, std::string& cig) override {
        aligner.alignEnd2End(pattern.data(), pattern.size(), text.data(), text.size());
        cig = aligner.getCIGAR(true);       // '=' / X / I (pattern-only) / D (text-only)
        return -aligner.getAlignmentScore();
    }
};
static std::unique_ptr<Aligner> make(const std::string& method) {
    if (method == "optik-wfa") return std::make_unique<OptikWfa>();
    if (method == "optik-biwfa") return std::make_unique<OptikBiWfa>();
    if (method == "wfa2-high") return std::make_unique<Wfa2>(wfa::WFAligner::MemoryHigh);
    if (method == "wfa2-ultralow") return std::make_unique<Wfa2>(wfa::WFAligner::MemoryUltralow);
    return nullptr;
}

static void mutate_pairs(int len, double err, int count, std::vector<std::string>& patterns, std::vector<std::string>& texts) {
    std::mt19937_64 rng(len * kSeedLenMul + int(err * kSeedErrScale));
    const char* bases = "ACGT";
    for (int i = 0; i < count; ++i) {
        std::string pattern; for (int j = 0; j < len; ++j) pattern.push_back(bases[rng() % kNumBases]);
        std::string text;
        for (char base_char : pattern) {
            double roll = std::uniform_real_distribution<double>(0, 1)(rng);
            if (roll < err / kEditKinds) continue;
            if (roll < 2 * err / kEditKinds) { text.push_back(bases[rng() % kNumBases]); text.push_back(base_char); continue; }
            if (roll < err) { text.push_back(bases[rng() % kNumBases]); continue; }
            text.push_back(base_char);
        }
        patterns.push_back(pattern); texts.push_back(text);
    }
}

int main(int argc, char** argv) {
    if (argc == kFileModeArgc) {
        auto aligner = make(argv[1]);
        if (!aligner) { std::fprintf(stderr, "unknown method %s\n", argv[1]); return 2; }
        std::vector<std::string> patterns, texts;
        std::ifstream in_file(argv[2]);
        std::string line;
        while (std::getline(in_file, line)) {
            const auto tab = line.find('\t');
            if (tab == std::string::npos) continue;
            patterns.push_back(line.substr(0, tab)); texts.push_back(line.substr(tab + 1));
        }
        std::vector<int> scores(patterns.size());
        std::vector<std::string> cigars(patterns.size());
        auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < patterns.size(); ++i) scores[i] = aligner->align(patterns[i], texts[i], cigars[i]);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        int invalid = 0;
        std::ofstream out(argv[kOutPathArg]);
        for (std::size_t i = 0; i < patterns.size(); ++i) {
            invalid += rescore(cigars[i], patterns[i], texts[i]) != scores[i];
            out << scores[i] << '\n';
        }
        std::printf("align_seconds %.6f invalid %d pairs %zu\n", secs, invalid, patterns.size());
        return 0;
    }
    // self-check
    const char* methods[] = {"optik-wfa", "optik-biwfa", "wfa2-high", "wfa2-ultralow"};
    int bad = 0;
    struct Cfg { int len; double err; int num_pairs; };
    for (Cfg cfg : {Cfg{150, 0.01, 20000}, Cfg{150, 0.05, 20000}, Cfg{1000, 0.05, 2000},
                  Cfg{10000, 0.05, 50}, Cfg{10000, 0.10, 20}}) {
        std::vector<std::string> patterns, texts;
        mutate_pairs(cfg.len, cfg.err, cfg.num_pairs, patterns, texts);
        std::vector<int> ref;
        std::printf("len=%-6d err=%.2f pairs=%-6d", cfg.len, cfg.err, cfg.num_pairs);
        for (const char* method : methods) {
            auto aligner = make(method);
            std::vector<int> scores(cfg.num_pairs); std::vector<std::string> cigars(cfg.num_pairs);
            auto start = std::chrono::steady_clock::now();
            for (int i = 0; i < cfg.num_pairs; ++i) scores[i] = aligner->align(patterns[i], texts[i], cigars[i]);
            const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            int diff = 0, inv = 0;
            if (ref.empty()) ref = scores;
            for (int i = 0; i < cfg.num_pairs; ++i) { diff += scores[i] != ref[i]; inv += rescore(cigars[i], patterns[i], texts[i]) != scores[i]; }
            bad += diff + inv;
            std::printf("  %s %.2f us%s", method, secs * kUsPerSec / cfg.num_pairs,
                        (diff || inv) ? (" [DIFF " + std::to_string(diff) + " INVALID " + std::to_string(inv) + "]").c_str() : "");
        }
        std::printf("\n");
    }
    std::printf(bad ? "wfa_affine_bench: %d SCORE/CIGAR DIFFERENCES\n"
                    : "wfa_affine_bench: scores identical, CIGARs valid (optik WFA/BiWFA, WFA2 WFA/BiWFA)\n", bad);
    return bad != 0;
}
