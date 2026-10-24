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
#include "bindings/cpp/WFAligner.hpp"
using namespace tpoptoa;

static constexpr AffinePenalties PEN{4, 6, 2};

/// Re-score a CIGAR whose ops are M match, X mismatch, I pattern-only and D text-only, with WFA2-lib's '=' read as M.
/// Returns -1 if the CIGAR does not consume both sequences exactly or an M covers a mismatch.
static int rescore(const std::string& cig, const std::string& p, const std::string& t) {
    std::size_t i = 0, j = 0, q = 0; int sc = 0;
    while (q < cig.size()) {
        long n = 0;
        if (!isdigit(static_cast<unsigned char>(cig[q]))) n = 1;     // WFA2 run-less op
        while (q < cig.size() && isdigit(static_cast<unsigned char>(cig[q]))) n = n * 10 + (cig[q++] - '0');
        const char op = cig[q++];
        if (op == 'M' || op == '=') {
            for (long r = 0; r < n; ++r, ++i, ++j)
                if (i >= p.size() || j >= t.size() || p[i] != t[j]) return -1;
        } else if (op == 'X') {
            for (long r = 0; r < n; ++r, ++i, ++j) {
                if (i >= p.size() || j >= t.size()) return -1;
                sc += PEN.x;
            }
        } else if (op == 'I') { sc += PEN.o + PEN.e * n; i += n; }
        else if (op == 'D') { sc += PEN.o + PEN.e * n; j += n; }
        else return -1;
    }
    return (i == p.size() && j == t.size()) ? sc : -1;
}

struct Aligner {
    virtual ~Aligner() = default;
    /// score and CIGAR (M/X/I/D)
    virtual int align(const std::string& p, const std::string& t, std::string& cig) = 0;
};
struct OptikWfa : Aligner {
    WfaAffine w{PEN};
    int align(const std::string& p, const std::string& t, std::string& cig) override {
        auto r = w.align(p.data(), p.size(), t.data(), t.size(), 0, 0);
        cig = r.cigar; return r.score;
    }
};
struct OptikBiWfa : Aligner {
    BiWfa b{PEN};
    int align(const std::string& p, const std::string& t, std::string& cig) override {
        auto r = b.align(p.data(), p.size(), t.data(), t.size());
        cig = r.cigar; return r.score;
    }
};
struct Wfa2 : Aligner {
    wfa::WFAlignerGapAffine a;
    explicit Wfa2(wfa::WFAligner::MemoryModel mm)
        : a(PEN.x, PEN.o, PEN.e, wfa::WFAligner::Alignment, mm) { a.setHeuristicNone(); }
    int align(const std::string& p, const std::string& t, std::string& cig) override {
        a.alignEnd2End(p.data(), p.size(), t.data(), t.size());
        cig = a.getCIGAR(true);       // '=' / X / I (pattern-only) / D (text-only)
        return -a.getAlignmentScore();
    }
};
static std::unique_ptr<Aligner> make(const std::string& m) {
    if (m == "optik-wfa") return std::make_unique<OptikWfa>();
    if (m == "optik-biwfa") return std::make_unique<OptikBiWfa>();
    if (m == "wfa2-high") return std::make_unique<Wfa2>(wfa::WFAligner::MemoryHigh);
    if (m == "wfa2-ultralow") return std::make_unique<Wfa2>(wfa::WFAligner::MemoryUltralow);
    return nullptr;
}

static void mutate_pairs(int len, double err, int n, std::vector<std::string>& P, std::vector<std::string>& T) {
    std::mt19937_64 rng(len * 7 + int(err * 100));
    const char* A = "ACGT";
    for (int i = 0; i < n; ++i) {
        std::string p; for (int j = 0; j < len; ++j) p.push_back(A[rng() % 4]);
        std::string t;
        for (char ch : p) {
            double u = std::uniform_real_distribution<double>(0, 1)(rng);
            if (u < err / 3) continue;
            if (u < 2 * err / 3) { t.push_back(A[rng() % 4]); t.push_back(ch); continue; }
            if (u < err) { t.push_back(A[rng() % 4]); continue; }
            t.push_back(ch);
        }
        P.push_back(p); T.push_back(t);
    }
}

int main(int argc, char** argv) {
    if (argc == 4) {
        auto al = make(argv[1]);
        if (!al) { std::fprintf(stderr, "unknown method %s\n", argv[1]); return 2; }
        std::vector<std::string> P, T;
        std::ifstream in(argv[2]);
        std::string line;
        while (std::getline(in, line)) {
            const auto tab = line.find('\t');
            if (tab == std::string::npos) continue;
            P.push_back(line.substr(0, tab)); T.push_back(line.substr(tab + 1));
        }
        std::vector<int> sc(P.size());
        std::vector<std::string> cg(P.size());
        auto t0 = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < P.size(); ++i) sc[i] = al->align(P[i], T[i], cg[i]);
        const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        int invalid = 0;
        std::ofstream out(argv[3]);
        for (std::size_t i = 0; i < P.size(); ++i) {
            invalid += rescore(cg[i], P[i], T[i]) != sc[i];
            out << sc[i] << '\n';
        }
        std::printf("align_seconds %.6f invalid %d pairs %zu\n", s, invalid, P.size());
        return 0;
    }
    // self-check
    const char* methods[] = {"optik-wfa", "optik-biwfa", "wfa2-high", "wfa2-ultralow"};
    int bad = 0;
    struct Cfg { int len; double err; int n; };
    for (Cfg c : {Cfg{150, 0.01, 20000}, Cfg{150, 0.05, 20000}, Cfg{1000, 0.05, 2000},
                  Cfg{10000, 0.05, 50}, Cfg{10000, 0.10, 20}}) {
        std::vector<std::string> P, T;
        mutate_pairs(c.len, c.err, c.n, P, T);
        std::vector<int> ref;
        std::printf("len=%-6d err=%.2f pairs=%-6d", c.len, c.err, c.n);
        for (const char* m : methods) {
            auto al = make(m);
            std::vector<int> sc(c.n); std::vector<std::string> cg(c.n);
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < c.n; ++i) sc[i] = al->align(P[i], T[i], cg[i]);
            const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            int diff = 0, inv = 0;
            if (ref.empty()) ref = sc;
            for (int i = 0; i < c.n; ++i) { diff += sc[i] != ref[i]; inv += rescore(cg[i], P[i], T[i]) != sc[i]; }
            bad += diff + inv;
            std::printf("  %s %.2f us%s", m, s * 1e6 / c.n,
                        (diff || inv) ? (" [DIFF " + std::to_string(diff) + " INVALID " + std::to_string(inv) + "]").c_str() : "");
        }
        std::printf("\n");
    }
    std::printf(bad ? "wfa_affine_bench: %d SCORE/CIGAR DIFFERENCES\n"
                    : "wfa_affine_bench: scores identical, CIGARs valid (optik WFA/BiWFA, WFA2 WFA/BiWFA)\n", bad);
    return bad != 0;
}
