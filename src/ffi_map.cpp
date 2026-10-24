/*
 * ffi_map.cpp - C ABI of the read mapper's native parts, i.e. the minimizer seed index (seed_index.hpp), bit-parallel candidate verification (bpm_edit.hpp) and the gap-affine wavefront aligners (wfa_affine.hpp).
 * The Rust side is tools/src/map.rs.
 */
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "../include/bpm_edit.hpp"
#include "../include/seed_index.hpp"
#include "../include/wfa_affine.hpp"

namespace {

/// Texts bpm_semiglobal_ed_x4() aligns at once, one per AVX2 lane.
constexpr int kBpmLanes = 4;

/// copy_cigar() - copy a CIGAR into a bounded C buffer, always NUL-terminating.
/// Returns the bytes copied, excluding the NUL, or -1 if the buffer is unusable.
int copy_cigar(const std::string& cigar, char* out_cig, int cig_cap) {
    if (!out_cig || cig_cap <= 0) return -1;
    const std::size_t space = static_cast<std::size_t>(cig_cap) - 1;
    const std::size_t n = std::min(cigar.size(), space);
    std::memcpy(out_cig, cigar.data(), n);
    out_cig[n] = '\0';
    return static_cast<int>(n);
}

/// fold_mismatches() - @cigar with its X (mismatch) runs merged into the neighbouring M runs, as SAM's M covers matches and mismatches alike.
std::string fold_mismatches(const std::string& cigar) {
    std::string cig;
    long run = 0; char op = 0;
    std::size_t i = 0;
    auto flush = [&]() { if (run) { cig += std::to_string(run); cig += op; } };
    while (i < cigar.size()) {
        long n = 0;
        while (i < cigar.size() && cigar[i] >= '0' && cigar[i] <= '9') n = n * 10 + (cigar[i++] - '0');
        char c = cigar[i++];
        if (c == 'X') c = 'M';
        if (c == op) run += n; else { flush(); op = c; run = n; }
    }
    flush();
    return cig;
}

}

extern "C" {

/// OkSeedIdx - minimizer seed index over reference sequences.
/// The sequences are not copied, so they must outlive the handle.
struct OkSeedIdx { tpoptoa::SeedIndex idx; };

/// okseed_build() - SeedIndex::build() over @n sequences, or nullptr on failure.
OkSeedIdx* okseed_build(const char* const* seqs, const uint64_t* lens, uint64_t n,
                        uint32_t k, uint32_t w, uint32_t max_occ, uint32_t nthreads) {
    try {
        std::vector<std::pair<const char*, std::size_t>> refs;
        for (uint64_t i = 0; i < n; ++i) refs.push_back({seqs[i], static_cast<std::size_t>(lens[i])});
        auto* h = new OkSeedIdx();
        if (!h->idx.build(refs, k, w, max_occ, nthreads)) { delete h; return nullptr; }
        return h;
    } catch (...) {
        return nullptr;
    }
}
void okseed_free(OkSeedIdx* h) { delete h; }
uint64_t okseed_distinct(const OkSeedIdx* h) { return h->idx.distinct(); }
uint64_t okseed_bytes(const OkSeedIdx* h) { return h->idx.bytes(); }
uint64_t okseed_occurrences(const OkSeedIdx* h) { return h->idx.occurrences(); }

/// okseed_collect() - seed hits of one read, with *out pointing at a thread-local buffer valid until the next call on this thread.
/// Returns the number of hits.
/// When not null, *n_min and *n_rep get the read's minimizer count and how many were dropped as over-occurring.
uint64_t okseed_collect(const OkSeedIdx* h, const char* seq, uint64_t len,
                        const tpoptoa::SeedIndex::Hit** out, uint32_t* n_min, uint32_t* n_rep) {
    thread_local std::vector<tpoptoa::SeedIndex::Hit> buf;
    try {
        h->idx.collect(seq, static_cast<std::size_t>(len), buf, n_min, n_rep);
    } catch (...) {
        buf.clear();
    }
    *out = buf.data();
    return buf.size();
}

/// okseed_collect_chain() - seed hits of one read for chaining (SeedIndex::collect_chain(), where repeat runs keep their least frequent minimizers), and *rep_len gets the read bases covered only by dropped minimizers.
/// Same buffer contract as okseed_collect.
uint64_t okseed_collect_chain(const OkSeedIdx* h, const char* seq, uint64_t len,
                              const tpoptoa::SeedIndex::Hit** out, uint32_t* rep_len) {
    thread_local std::vector<tpoptoa::SeedIndex::Hit> buf;
    try {
        h->idx.collect_chain(seq, static_cast<std::size_t>(len), buf, rep_len);
    } catch (...) {
        buf.clear();
    }
    *out = buf.data();
    return buf.size();
}
uint32_t okseed_mid_occ(const OkSeedIdx* h) { return h->idx.mid_occ(); }

/// ok_placements_x4() - placements of one query, end to end, in @cnt <= 4 texts with free ends, at unit costs capped at max_ed + 1.
/// Per text it gives the best edit distance with its leftmost exclusive end, and the best distance among ends more than `sep` away, i.e. a distinct placement such as one tandem-repeat period over.
/// It runs bit-parallel with one SIMD lane per text through bpm_edit.hpp, in O(t_len * ceil(q_len / 64)) whatever the distance.
void ok_placements_x4(const char* query, int32_t q_len,
                      const char* const* texts, const int32_t* t_lens, int32_t cnt,
                      int32_t max_ed, int32_t sep, int32_t* out_best, int32_t* out_end,
                      int32_t* out_alt, int32_t* out_alt_end) {
    const int c = std::min<int32_t>(std::max<int32_t>(cnt, 0), kBpmLanes);
    try {
        thread_local std::vector<int> prof;
        int n[kBpmLanes], o[kBpmLanes], off[kBpmLanes];
        std::size_t tot = 0;
        for (int i = 0; i < c; ++i) { n[i] = t_lens[i]; off[i] = static_cast<int>(tot); tot += static_cast<std::size_t>(n[i]); }
        prof.resize(tot + 1);
        int* pp[kBpmLanes];
        for (int i = 0; i < c; ++i) pp[i] = prof.data() + off[i];
        tpoptoa::bpm_semiglobal_ed_x4(query, q_len, texts, n, c, max_ed, o, pp);
        for (int i = 0; i < c; ++i) {
            const tpoptoa::BpmPlacements pl = tpoptoa::bpm_placements(pp[i], n[i], sep, max_ed);
            out_best[i] = o[i];
            out_end[i] = pl.end;
            out_alt[i] = pl.alt;
            out_alt_end[i] = pl.alt_end;
        }
    } catch (...) {
        for (int i = 0; i < c; ++i) { out_best[i] = out_alt[i] = max_ed + 1; out_end[i] = out_alt_end[i] = 0; }
    }
}

/// ok_wfa_affine_align() - gap-affine alignment, query end to end from text position 0, text end free within its last `text_end_free` bases.
/// CIGAR keeps M for matches and mismatches.
/// Returns 0 on success with the score and text end written out, or -1 on failure or when the score cap is exceeded.
int32_t ok_wfa_affine_align(const char* query, int32_t q_len,
                            const char* text, int32_t t_len, int32_t text_end_free,
                            int32_t x, int32_t o, int32_t e, int32_t max_score,
                            char* out_cig, int32_t cig_cap, int32_t* out_score,
                            int32_t* out_tend) {
    try {
        thread_local tpoptoa::WfaAffine aligner;
        tpoptoa::AffinePenalties pen; pen.x = x; pen.o = o; pen.e = e;
        if (aligner.penalties().x != x || aligner.penalties().o != o || aligner.penalties().e != e)
            aligner = tpoptoa::WfaAffine(pen);
        auto r = aligner.align(query, q_len, text, t_len, 0, text_end_free, max_score);
        if (r.score < 0) return -1;
        const std::string cig = fold_mismatches(r.cigar);
        if (copy_cigar(cig, out_cig, cig_cap) < 0) return -1;
        *out_score = r.score; *out_tend = r.text_end;
        return 0;
    } catch (...) {
        return -1;
    }
}

/// ok_biwfa_align() - gap-affine end-to-end alignment in O(s) memory (BiWfa).
/// CIGAR keeps M for matches and mismatches.
/// Returns 0 on success, -1 on failure / score cap exceeded.
int32_t ok_biwfa_align(const char* query, int32_t q_len, const char* text, int32_t t_len,
                       int32_t x, int32_t o, int32_t e, int32_t max_score,
                       char* out_cig, int32_t cig_cap, int32_t* out_score) {
    try {
        thread_local tpoptoa::BiWfa aligner;
        if (aligner.penalties().x != x || aligner.penalties().o != o || aligner.penalties().e != e)
            aligner = tpoptoa::BiWfa(tpoptoa::AffinePenalties{x, o, e});
        auto r = aligner.align(query, q_len, text, t_len, max_score);
        if (r.score < 0) return -1;
        const std::string cig = fold_mismatches(r.cigar);
        if (copy_cigar(cig, out_cig, cig_cap) < 0) return -1;
        *out_score = r.score;
        return 0;
    } catch (...) {
        return -1;
    }
}

/// ok_wfa_affine_cigar() - gap-affine semi-global alignment, with the query end to end and both text ends free.
/// Mismatch runs are folded into 'M' for SAM-compatible M/I/D.
/// Returns the 0-based text start, or -1 on error, and stores the affine score in *out_score.
int32_t ok_wfa_affine_cigar(const char* query, int32_t q_len,
                            const char* text, int32_t t_len,
                            int32_t x, int32_t o, int32_t e,
                            char* out_cig, int32_t cig_cap, int32_t* out_score) {
    try {
        thread_local tpoptoa::WfaAffine aligner;
        tpoptoa::AffinePenalties pen; pen.x = x; pen.o = o; pen.e = e;
        if (aligner.penalties().x != x || aligner.penalties().o != o || aligner.penalties().e != e)
            aligner = tpoptoa::WfaAffine(pen);
        auto r = aligner.align(query, q_len, text, t_len, t_len, t_len);
        if (r.score < 0) return -1;
        const std::string cig = fold_mismatches(r.cigar);
        if (copy_cigar(cig, out_cig, cig_cap) < 0) return -1;
        if (out_score) *out_score = r.score;
        return r.text_begin;
    } catch (...) {
        return -1;
    }
}

}  // extern "C"
