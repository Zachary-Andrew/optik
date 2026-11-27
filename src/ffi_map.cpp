/**
 * @file
 * @brief C ABI for the read mapper's native parts.
 * These are the minimizer seed index (seed_index.hpp), bit-parallel candidate verification (bpm_edit.hpp) and the gap-affine wavefront aligners (wfa_affine.hpp).
 * The matching Rust code is tools/src/map.rs.
 */
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "../include/bpm_edit.hpp"
#include "../include/seed_index.hpp"
#include "../include/wfa_affine.hpp"
#include "ffi_error.hpp"

using tpoptoa::ffi::kStatusFailed;
using tpoptoa::ffi::kStatusOk;
using tpoptoa::kBpmLanes;

namespace {

/// Copies a CIGAR into a bounded C buffer and always NUL-terminates it.
/// The result is the byte count without the NUL, or kStatusFailed if the buffer is unusable.
int copy_cigar(const std::string& cigar, char* out_cig, int cig_cap) {
    if (!out_cig || cig_cap <= 0) return kStatusFailed;

    const std::size_t space = static_cast<std::size_t>(cig_cap) - 1;
    const std::size_t copy_len = std::min(cigar.size(), space);

    std::memcpy(out_cig, cigar.data(), copy_len);
    out_cig[copy_len] = '\0';
    return static_cast<int>(copy_len);
}

/// Merges X (mismatch) runs into the neighbouring M runs, since SAM's M covers both matches and mismatches.
std::string fold_mismatches(const std::string& cigar) {
    std::string cig;
    long run = 0; char cur_op = 0;
    std::size_t pos = 0;
    auto flush = [&]() { if (run) { cig += std::to_string(run); cig += cur_op; } };

    while (pos < cigar.size()) {
        long run_len = 0;

        while (pos < cigar.size() && cigar[pos] >= '0' && cigar[pos] <= '9') run_len = run_len * tpoptoa::kDecimalBase + (cigar[pos++] - '0');

        char op_char = cigar[pos++];

        if (op_char == 'X') op_char = 'M';

        if (op_char == cur_op) run += run_len; else { flush(); cur_op = op_char; run = run_len; }
    }

    flush();
    return cig;
}

}

extern "C" {

/// Minimizer seed index over reference sequences.
/// Nothing is copied, so the sequences must outlive the handle.
struct OkSeedIdx { tpoptoa::SeedIndex idx; };

/// Returns nullptr if SeedIndex::build() fails.
OkSeedIdx* okseed_build(const char* const* seqs, const uint64_t* lens, uint64_t num_seqs,
                        uint32_t kmer_len, uint32_t window, uint32_t max_occ, uint32_t nthreads) {
    try {
        std::vector<std::pair<const char*, std::size_t>> refs;

        for (uint64_t i = 0; i < num_seqs; ++i) refs.push_back({seqs[i], static_cast<std::size_t>(lens[i])});

        auto* handle = new OkSeedIdx();

        if (!handle->idx.build(refs, kmer_len, window, max_occ, nthreads)) { delete handle; return nullptr; }

        return handle;
    } catch (...) {
        return nullptr;
    }
}
void okseed_free(OkSeedIdx* handle) { delete handle; }
uint64_t okseed_distinct(const OkSeedIdx* handle) { return handle->idx.distinct(); }
uint64_t okseed_bytes(const OkSeedIdx* handle) { return handle->idx.bytes(); }
uint64_t okseed_occurrences(const OkSeedIdx* handle) { return handle->idx.occurrences(); }

/// Returns the hit count and points *out at a thread-local buffer that stays valid until this thread calls again.
/// If provided, *n_min receives the read's minimizer count and *n_rep how many of them were dropped as over-occurring.
uint64_t okseed_collect(const OkSeedIdx* handle, const char* seq, uint64_t len,
                        const tpoptoa::SeedIndex::Hit** out, uint32_t* n_min, uint32_t* n_rep) {
    thread_local std::vector<tpoptoa::SeedIndex::Hit> buf;

    try {
        handle->idx.collect(seq, static_cast<std::size_t>(len), buf, n_min, n_rep);
    } catch (...) {
        buf.clear();
    }
    *out = buf.data();
    return buf.size();
}

/// Collects hits for chaining, where repeat runs keep their least frequent minimizers (SeedIndex::collect_chain()).
/// Meanwhile *rep_len receives the read bases covered only by dropped minimizers.
/// The buffer contract matches okseed_collect().
uint64_t okseed_collect_chain(const OkSeedIdx* handle, const char* seq, uint64_t len,
                              const tpoptoa::SeedIndex::Hit** out, uint32_t* rep_len) {
    thread_local std::vector<tpoptoa::SeedIndex::Hit> buf;

    try {
        handle->idx.collect_chain(seq, static_cast<std::size_t>(len), buf, rep_len);
    } catch (...) {
        buf.clear();
    }
    *out = buf.data();
    return buf.size();
}
uint32_t okseed_mid_occ(const OkSeedIdx* handle) { return handle->idx.mid_occ(); }

/// Places one query end to end in up to kBpmLanes free-ended texts, at unit costs capped at max_ed + 1.
/// For each text it reports the best distance with its leftmost exclusive end, plus the best distance among ends more than `sep` away.
/// That second figure catches distinct placements, such as one tandem-repeat period over.
/// Since it runs bit-parallel with one SIMD lane per text (bpm_edit.hpp), the cost is O(t_len * ceil(q_len / 64)) at any distance.
void ok_placements_x4(const char* query, int32_t q_len,
                      const char* const* texts, const int32_t* t_lens, int32_t cnt,
                      int32_t max_ed, int32_t sep, int32_t* out_best, int32_t* out_end,
                      int32_t* out_alt, int32_t* out_alt_end) {
    const int lanes_used = std::min<int32_t>(std::max<int32_t>(cnt, 0), kBpmLanes);

    try {
        thread_local std::vector<int> prof;
        int text_lens[kBpmLanes], best_eds[kBpmLanes], off[kBpmLanes];
        std::size_t tot = 0;

        // All lanes share one profile buffer, each at its own offset.
        for (int i = 0; i < lanes_used; ++i) { text_lens[i] = t_lens[i]; off[i] = static_cast<int>(tot); tot += static_cast<std::size_t>(text_lens[i]); }

        prof.resize(tot + 1);
        int* profiles[kBpmLanes];

        for (int i = 0; i < lanes_used; ++i) profiles[i] = prof.data() + off[i];

        tpoptoa::bpm_semiglobal_ed_x4(query, q_len, texts, text_lens, lanes_used, max_ed, best_eds, profiles);

        for (int i = 0; i < lanes_used; ++i) {
            // The profile holds the edit distance at every text end. Its best end and the runner-up further than sep away become the placements.
            const tpoptoa::BpmPlacements placements = tpoptoa::bpm_placements(profiles[i], text_lens[i], sep, max_ed);

            out_best[i] = best_eds[i];
            out_end[i] = placements.end;
            out_alt[i] = placements.alt;
            out_alt_end[i] = placements.alt_end;
        }
    } catch (...) {
        // On failure every lane reports no placement within max_ed.
        for (int i = 0; i < lanes_used; ++i) { out_best[i] = out_alt[i] = max_ed + 1; out_end[i] = out_alt_end[i] = 0; }
    }
}

/// Gap-affine alignment of the whole query from text position 0, with the text end free within its last `text_end_free` bases.
/// On success the score and text end are written out.
/// Failures and exceeded score caps both return kStatusFailed.
int32_t ok_wfa_affine_align(const char* query, int32_t q_len,
                            const char* text, int32_t t_len, int32_t text_end_free,
                            int32_t mismatch, int32_t gap_open, int32_t gap_extend, int32_t max_score,
                            char* out_cig, int32_t cig_cap, int32_t* out_score,
                            int32_t* out_tend) {
    try {
        thread_local tpoptoa::WfaAffine aligner;
        tpoptoa::AffinePenalties pen; pen.mismatch = mismatch; pen.gap_open = gap_open; pen.gap_extend = gap_extend;

        if (aligner.penalties().mismatch != mismatch || aligner.penalties().gap_open != gap_open || aligner.penalties().gap_extend != gap_extend)
            aligner = tpoptoa::WfaAffine(pen);

        auto result = aligner.align(query, q_len, text, t_len, 0, text_end_free, max_score);

        if (result.score < 0) return kStatusFailed;

        const std::string cig = fold_mismatches(result.cigar);

        if (copy_cigar(cig, out_cig, cig_cap) < 0) return kStatusFailed;

        *out_score = result.score; *out_tend = result.text_end;
        return kStatusOk;
    } catch (...) {
        return kStatusFailed;
    }
}

/// End-to-end gap-affine alignment in O(s) memory (BiWfa).
/// Failures and exceeded score caps both return kStatusFailed.
int32_t ok_biwfa_align(const char* query, int32_t q_len, const char* text, int32_t t_len,
                       int32_t mismatch, int32_t gap_open, int32_t gap_extend, int32_t max_score,
                       char* out_cig, int32_t cig_cap, int32_t* out_score) {
    try {
        thread_local tpoptoa::BiWfa aligner;

        if (aligner.penalties().mismatch != mismatch || aligner.penalties().gap_open != gap_open || aligner.penalties().gap_extend != gap_extend)
            aligner = tpoptoa::BiWfa(tpoptoa::AffinePenalties{mismatch, gap_open, gap_extend});

        auto result = aligner.align(query, q_len, text, t_len, max_score);

        if (result.score < 0) return kStatusFailed;

        const std::string cig = fold_mismatches(result.cigar);

        if (copy_cigar(cig, out_cig, cig_cap) < 0) return kStatusFailed;

        *out_score = result.score;
        return kStatusOk;
    } catch (...) {
        return kStatusFailed;
    }
}

/// Semi-global gap-affine alignment with the whole query and both text ends free.
/// It returns the 0-based text start, or kStatusFailed on error, and stores the affine score in *out_score.
int32_t ok_wfa_affine_cigar(const char* query, int32_t q_len,
                            const char* text, int32_t t_len,
                            int32_t mismatch, int32_t gap_open, int32_t gap_extend,
                            char* out_cig, int32_t cig_cap, int32_t* out_score) {
    try {
        thread_local tpoptoa::WfaAffine aligner;
        tpoptoa::AffinePenalties pen; pen.mismatch = mismatch; pen.gap_open = gap_open; pen.gap_extend = gap_extend;

        if (aligner.penalties().mismatch != mismatch || aligner.penalties().gap_open != gap_open || aligner.penalties().gap_extend != gap_extend)
            aligner = tpoptoa::WfaAffine(pen);

        auto result = aligner.align(query, q_len, text, t_len, t_len, t_len);

        if (result.score < 0) return kStatusFailed;

        const std::string cig = fold_mismatches(result.cigar);

        if (copy_cigar(cig, out_cig, cig_cap) < 0) return kStatusFailed;

        if (out_score) *out_score = result.score;

        return result.text_begin;
    } catch (...) {
        return kStatusFailed;
    }
}

}  // extern "C"
