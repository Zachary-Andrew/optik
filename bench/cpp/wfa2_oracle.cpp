// wfa2_oracle.cpp - external exact-alignment oracle for optik's edit-distance engine.
//
// For (query, text) pairs read from stdin as TSV lines (`query<TAB>text`, upper-case ACGT), it computes the unit-cost edit distance by independent legs.
//
//   optik_sg     optik's own wavefront aligner (WfaAffine, wfa_affine.hpp) with edit penalties x=1 o=0 e=1, the query end-to-end and the text ends free, uncapped.
//                The mapper aligns with this engine under affine penalties, and its candidate verification uses the bit-parallel kernel of bpm_edit.hpp, which kernels_stress.cpp checks against a naive DP.
//
//   optik_capped that engine with the mapper's score cap set to the true
//                distance, where the capped path must return the same score.
//
//   wfa2_sg      WFA2-lib ends-free alignment under the edit metric.
//                It needs the vendored ends-free boundary-key fix (bench/patches/wfa2_edit_endsfree.py), because upstream WFA2 silently disables ends-free initialization for the edit and indel metrics, where match==0.
//                The score is read from cigar->score, which classic WFA sets to wf_score for edit, and not from the CIGAR string, since upstream writes the free text ends as I operations that a naive re-scoring would over-count.
//
//   wfa2_e2e     WFA2-lib end-to-end alignment under the edit metric with wf_heuristic_none, exact and without heuristic pruning.
//                It applies only when the text is the query's locus, so the driver feeds locus-length texts for this leg.
//
//   naive_e2e    a 10-line O(m*n) rolling-row unit-cost Levenshtein DP, written independently, as the final arbiter for small pairs (m*n <= 4e6, else -1).
//
// The harness fails hard with exit 1 on any disagreement between optik_sg and wfa2_sg, between optik_sg and optik_capped, and between wfa2_e2e and naive_e2e when computed.
// The driver (waveform_bench.py) also checks every pair against the known ground-truth edit distance and ties the legs together per tier.
//
// This "waveform reference" bench confirms independently, with an entirely different code base (WFA2-lib's C wavefront engine plus a naive DP), the exact scores the mapper depends on, including the 150 bp / 10% error regime where wavefront heuristics once produced wrong scores.
//
// Usage:  wfa2_oracle < pairs.tsv
// Each output line is m<TAB>n<TAB>optik_sg<TAB>optik_capped<TAB>wfa2_sg<TAB>wfa2_e2e<TAB>naive_e2e.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "wfa_affine.hpp"                  // tpoptoa::WfaAffine, optik's aligner
// WFA2-lib is a C library whose headers lack extern "C" guards.
extern "C" {
#include "wavefront/wfa.h"                 // wavefront_aligner_t
#include "wavefront/wavefront_align.h"     // wavefront_align
#include "wavefront/wavefront_aligner.h"   // wavefront_aligner_new/delete/reap

/// optik's WFA in edit mode: unit mismatch, no gap open, unit gap extend.
static constexpr tpoptoa::AffinePenalties kEditPenalties{1, 0, 1};
/// The naive O(m n) arbiter runs only up to this many DP cells.
/// Pairs it skips are marked kNotComputed.
static constexpr long kNaiveMaxCells = 4000000L;
static constexpr int kNotComputed = -1;
}

/// Normalize a sequence to upper case with non-ACGT mapped to 'A'.
/// The generators emit only ACGT, and an unknown character mismatches everything in both engines.
static std::string norm(const char* seq, size_t len) {
    std::string out;
    out.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        unsigned char ch = (unsigned char)seq[i];
        if (ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
            ch = (unsigned char)toupper(ch);
        if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') ch = 'A';
        out.push_back((char)ch);
    }
    return out;
}

/// Independent O(m*n) unit-cost end-to-end Levenshtein over rolling rows.
static int naive_end2end_ed(const char* query, int query_len, const char* text, int text_len) {
    if (query_len == 0) return text_len;
    if (text_len == 0) return query_len;
    std::vector<int> prev(text_len + 1), cur(text_len + 1);
    for (int j = 0; j <= text_len; ++j) prev[j] = j;
    for (int i = 1; i <= query_len; ++i) {
        cur[0] = i;
        for (int j = 1; j <= text_len; ++j) {
            int best = prev[j] + 1;                                    // del
            if (cur[j - 1] + 1 < best) best = cur[j - 1] + 1;             // ins
            if (prev[j - 1] + (query[i - 1] == text[j - 1] ? 0 : 1) < best)   // sub
                best = prev[j - 1] + (query[i - 1] == text[j - 1] ? 0 : 1);
            cur[j] = best;
        }
        std::swap(prev, cur);
    }
    return prev[text_len];
}

int main(void) {
    // One aligner is reused across pairs, configured exactly for both legs.
    //   heuristic none  no wf-adaptive, xdrop, zdrop or banded shortcuts
    //   memory med      piggyback wavefronts, fully supported for ends-free (upstream warns that BiWFA and ultralow ends-free are untested, so this oracle stays conservative)
    //   edit metric     unit costs, match 0, mismatch 1, indel 1
    //   scope           compute_alignment, storing the score in cigar->score
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric          = edit;
    attributes.alignment_scope          = compute_alignment;
    attributes.heuristic.strategy       = wf_heuristic_none;
    attributes.memory_mode              = wavefront_memory_med;
    attributes.system.verbose           = 0;
    attributes.system.max_num_threads   = 1;

    wavefront_aligner_t* const aligner = wavefront_aligner_new(&attributes);
    if (aligner == NULL) {
        fprintf(stderr, "[wfa2_oracle] wavefront_aligner_new failed\n");
        return 2;
    }

    char* line = NULL;
    size_t cap = 0;
    long n_pairs = 0, n_sg_mismatch = 0, n_woah_mismatch = 0, n_e2e_mismatch = 0;
    printf("m\tn\toptik_sg\toptik_capped\twfa2_sg\twfa2_e2e\tnaive_e2e\n");

    while (getline(&line, &cap, stdin) >= 0) {
        size_t llen = strlen(line);
        while (llen > 0 && (line[llen - 1] == '\n' || line[llen - 1] == '\r'))
            line[--llen] = '\0';
        if (llen == 0) continue;

        char* tab = strchr(line, '\t');
        if (tab == NULL) continue;  // malformed line, skipped since the driver controls the input
        *tab = '\0';
        const size_t mlen = (size_t)(tab - line);
        const size_t tlen = llen - mlen - 1;
        if (mlen == 0 || tlen == 0) continue;

        std::string query = norm(line, mlen);
        std::string text = norm(tab + 1, tlen);
        const int query_len = (int)query.size();
        const int text_len = (int)text.size();

        // optik's own WFA in edit mode, with the query end-to-end and the text ends free
        static tpoptoa::WfaAffine optik_ed(kEditPenalties);
        const int optik_sg = optik_ed.align(query.data(), query_len, text.data(), text_len, text_len, text_len, tpoptoa::kNoScoreCap, false).score;
        // The same engine under the mapper's score cap, set to the true distance
        const int woah_sg = optik_ed.align(query.data(), query_len, text.data(), text_len, text_len, text_len, optik_sg, false).score;

        // WFA2 ends-free (semi-global), with the query fully consumed and the text free on both sides, the fit semantics of optik's D[0][j]=0 and min_j D[m][j].
        // The aligner's stored score is read, not the CIGAR.
        wavefront_aligner_set_alignment_free_ends(aligner, 0, 0, text_len, text_len);
        wavefront_align(aligner, query.data(), query_len, text.data(), text_len);
        const int wfa2_sg = aligner->cigar->score;

        // WFA2 end-to-end, which for locus-length texts is the true edit distance
        wavefront_aligner_set_alignment_end_to_end(aligner);
        wavefront_align(aligner, query.data(), query_len, text.data(), text_len);
        const int wfa2_e2e = aligner->cigar->score;

        // Naive DP arbiter for small pairs
        const int naive_e2e =
            ((long)query_len * (long)text_len <= kNaiveMaxCells) ? naive_end2end_ed(query.data(), query_len, text.data(), text_len) : kNotComputed;

        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", query_len, text_len, optik_sg, woah_sg,
               wfa2_sg, wfa2_e2e, naive_e2e);

        if (optik_sg != wfa2_sg) {
            ++n_sg_mismatch;
            fprintf(stderr, "[wfa2_oracle] SG MISMATCH m=%d n=%d optik=%d wfa2=%d q=%.24s t=%.24s\n",
                    query_len, text_len, optik_sg, wfa2_sg, query.c_str(), text.c_str());
        }
        if (optik_sg != woah_sg) {
            ++n_woah_mismatch;
            fprintf(stderr, "[wfa2_oracle] CAPPED MISMATCH m=%d n=%d optik=%d capped=%d q=%.24s t=%.24s\n",
                    query_len, text_len, optik_sg, woah_sg, query.c_str(), text.c_str());
        }
        if (naive_e2e >= 0 && wfa2_e2e != naive_e2e) {
            ++n_e2e_mismatch;
            fprintf(stderr, "[wfa2_oracle] E2E MISMATCH m=%d n=%d wfa2=%d naive=%d\n",
                    query_len, text_len, wfa2_e2e, naive_e2e);
        }
        ++n_pairs;

        wavefront_aligner_reap(aligner);
    }

    free(line);
    wavefront_aligner_delete(aligner);

    fprintf(stderr, "[wfa2_oracle] %ld pairs, sg_mismatch=%ld woah_mismatch=%ld e2e_mismatch=%ld\n",
            n_pairs, n_sg_mismatch, n_woah_mismatch, n_e2e_mismatch);
    return (n_sg_mismatch == 0 && n_woah_mismatch == 0 &&
            n_e2e_mismatch == 0) ? 0 : 1;
}