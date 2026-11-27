#pragma once

/**
 * @file
 * @brief Minimizer-length selection for SpssCore.
 *
 * estimate_m() predicts an index's size and lookup cost at one minimizer length by sampling a quarter of the minimizer space.
 * choose_m() walks candidate lengths with those estimates and keeps the cheapest.
 * spss_core.hpp includes this file at the end, so callers never include it directly.
 */

#include "spss_core.hpp"

namespace tpoptoa {

inline SpssCore::MEstimate SpssCore::estimate_m(const std::vector<std::string>& seqs, unsigned kmer_len, unsigned minimizer_len) {
    constexpr unsigned kSampleShift = 2;                 // 1/4 of the minimizers
    constexpr double kMeanTriesDen = 2.0;                // a lookup tries (n + 1) / 2 of a bucket's n records on average
    constexpr uint64_t kPartSeed = 0x7A3C1D, kSampleSeed = 0x5A17E5ULL;
    const unsigned window = kmer_len - minimizer_len + 1;
    const uint64_t mmask = kmer_mask(minimizer_len);
    uint64_t total_bases = 0;

    for (const auto& seq : seqs) total_bases += seq.size();

    struct KeyLen { uint64_t key; uint32_t len; };
    // Per-thread buffers are split into kParts minimizer-hash partitions, so the per-bucket aggregation below runs one partition per thread.
    constexpr unsigned kParts = 64;
    constexpr unsigned kPartBits = 6;                   // log2(kParts)
    auto part_of = [](uint64_t key) { return static_cast<unsigned>(mix64(key, kPartSeed) >> (kWordBits - kPartBits)); };   // top 6 bits, so kParts = 64
    unsigned nthr = 1;
#ifdef _OPENMP
    nthr = static_cast<unsigned>(omp_get_max_threads());
#endif
    std::vector<std::vector<std::vector<KeyLen>>> per(nthr, std::vector<std::vector<KeyLen>>(kParts));
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
    std::vector<uint64_t> hashes, canon_mers;                       // per-thread, reused
    std::vector<std::size_t> window_deque;

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, kScanSchedChunk)
#endif
    for (int64_t si = 0; si < static_cast<int64_t>(seqs.size()); ++si) {
#ifdef _OPENMP
        auto& out = per[static_cast<std::size_t>(omp_get_thread_num())];
#else
        auto& out = per[0];
#endif
        const std::string& seq = seqs[static_cast<std::size_t>(si)];
        std::size_t run_start = 0;

        while (run_start < seq.size()) {                          // each ACGT run of length >= k
            while (run_start < seq.size() && base_encode(seq[run_start]) == kInvalidBase) ++run_start;

            std::size_t run_end = run_start;

            while (run_end < seq.size() && base_encode(seq[run_end]) != kInvalidBase) ++run_end;

            if (run_end - run_start >= kmer_len) {
                const std::size_t num_mers = (run_end - run_start) - minimizer_len + 1, num_kmers = (run_end - run_start) - kmer_len + 1;

                hashes.resize(num_mers); canon_mers.resize(num_mers); window_deque.resize(num_mers);
                uint64_t fwd = 0, rev = 0;

                // First, every m-mer of the run gets its canonical form and hash.
                for (std::size_t t = 0; t < run_end - run_start; ++t) {
                    const uint64_t base = base_encode(seq[run_start + t]);

                    fwd = ((fwd << kBitsPerBase) | base) & mmask;
                    rev = (rev >> kBitsPerBase) | ((base ^ kComplementMask) << (kBitsPerBase * (minimizer_len - 1)));

                    if (t + 1 >= minimizer_len) {
                        const uint64_t canon = fwd < rev ? fwd : rev;

                        canon_mers[t + 1 - minimizer_len] = canon; hashes[t + 1 - minimizer_len] = mhash(canon);
                    }
                }

                // Then a monotone deque yields each k-mer's window minimizer, and consecutive k-mers that share one form a super-k-mer.
                std::size_t deque_head = 0, deque_tail = 0, next_in = 0, cur_q = ~std::size_t(0);
                uint32_t cur_len = 0;
                // A super-k-mer is sampled when its minimizer's seeded hash has the top kSampleShift bits clear, so the same quarter of minimizers is kept everywhere.
                auto emit = [&](std::size_t mer_idx, uint32_t len) {
                    if ((mix64(canon_mers[mer_idx], kSampleSeed) >> (kWordBits - kSampleShift)) == 0)
                        out[part_of(canon_mers[mer_idx])].push_back({canon_mers[mer_idx], len});
                };

                for (std::size_t st = 0; st < num_kmers; ++st) {
                    while (next_in <= st + (window - 1)) {
                        // A larger hash behind a new m-mer can never be a minimum again.
                        while (deque_tail > deque_head && hashes[window_deque[deque_tail - 1]] > hashes[next_in]) --deque_tail;

                        window_deque[deque_tail++] = next_in++;
                    }

                    // Dropping m-mers that slid out leaves the minimum at the front.
                    while (window_deque[deque_head] < st) ++deque_head;

                    // A new minimizer closes the current super-k-mer.
                    if (window_deque[deque_head] != cur_q) {
                        if (cur_len) emit(cur_q, cur_len);

                        cur_q = window_deque[deque_head]; cur_len = 0;
                    }

                    ++cur_len;
                }

                if (cur_len) emit(cur_q, cur_len);
            }

            run_start = run_end;
        }
    }
    }
    // Counts from the sample are scaled by 2^kSampleShift to estimate the full index.
    constexpr double scale = double(1u << kSampleShift);
    const unsigned lo_bits = bits_for(total_bases + 1), sk_bits = lo_bits + 1;

    struct Agg { double keys = 0, ovf_recs = 0, ovf_cnt_bits = 0, skew_bits = 0,
                 kmers = 0, scan = 0;
                 uint64_t max_cnt = 1; };
    std::vector<Agg> agg(kParts);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int p = 0; p < int(kParts); ++p) {
        std::vector<KeyLen> all;

        for (auto& thread_buf : per) { all.insert(all.end(), thread_buf[p].begin(), thread_buf[p].end()); std::vector<KeyLen>().swap(thread_buf[p]); }

        // Sorting by key gathers each minimizer's super-k-mers into one run, which is the bucket the index would build.
        std::sort(all.begin(), all.end(), [](const KeyLen& lhs, const KeyLen& rhs) { return lhs.key < rhs.key; });
        Agg& part_agg = agg[p];

        for (std::size_t run_start = 0; run_start < all.size();) {
            std::size_t run_end = run_start;
            uint64_t run_kmers = 0;

            while (run_end < all.size() && all[run_end].key == all[run_start].key) run_kmers += all[run_end++].len;

            const uint64_t run_count = run_end - run_start;

            part_agg.keys += 1; part_agg.kmers += double(run_kmers); part_agg.max_cnt = std::max(part_agg.max_cnt, run_count);

            // Every multi-record bucket keeps its record list in the overflow stream, and heavy ones add a skew table on top.
            if (run_count > 1) { part_agg.ovf_recs += double(run_count); part_agg.ovf_cnt_bits += gamma_len(run_count - 1); }

            if (run_count > kSkewMin) {
                part_agg.skew_bits += double(run_kmers) * kSkewLoadNum / kSkewLoadDen * double(cls_of(run_count) + kFpBits);   // load 0.9
                part_agg.scan += double(run_kmers);                       // one probe
            } else {
                part_agg.scan += double(run_kmers) * double(run_count + 1) / kMeanTriesDen; // mean records tried, (n + 1) / 2
            }

            run_start = run_end;
        }
    }

    double keys = 0, ovf_recs = 0, ovf_cnt_bits = 0, skew_bits = 0, kmers = 0, scan = 0;
    uint64_t max_cnt = 1;

    for (const Agg& part_agg : agg) {
        keys += part_agg.keys; ovf_recs += part_agg.ovf_recs; ovf_cnt_bits += part_agg.ovf_cnt_bits;
        skew_bits += part_agg.skew_bits; kmers += part_agg.kmers; scan += part_agg.scan;
        max_cnt = std::max(max_cnt, part_agg.max_cnt);
    }

    (void)max_cnt;
    // Record words must hold both a string position and an overflow offset. The record width follows whichever needs more bits.
    const double ovf_total = scale * (ovf_cnt_bits + ovf_recs * sk_bits);
    const unsigned rec_bits = std::max(sk_bits, bits_for(static_cast<uint64_t>(ovf_total) + 1)) + 1;
    const double est_keys = keys * scale;
    const double bits = est_keys * rec_bits + ovf_total + skew_bits * scale + est_keys * kMphfBitsPerKey;

    return {minimizer_len, bits, kmers > 0 ? scan / kmers : 0.0,
            est_keys, est_keys * rec_bits, ovf_total, skew_bits * scale, est_keys * kMphfBitsPerKey, kmers * scale};
}

inline unsigned SpssCore::choose_m(const std::vector<std::string>& seqs, unsigned kmer_len) {
    constexpr unsigned kFallbackMDiv = 2;   // below kMinChooseK, m is half of k
    const unsigned half = std::max(1u, kmer_len / kFallbackMDiv);

    if (kmer_len < kMinChooseK) return half;

    uint64_t bases = 0;

    for (const auto& seq : seqs) bases += seq.size();

    unsigned log4_bases = 0;

    while (log4_bases < kMaxLog4Bases && (uint64_t(1) << (kBitsPerBase * log4_bases)) < bases) ++log4_bases;   // ceil(log4(bases))

    const unsigned m_lo = std::max(kMinCandidateM, half > kCandidatesBelowHalf ? half - kCandidatesBelowHalf : 1u), m_hi = kmer_len - kMinMBelowK;
    const unsigned start_m = std::min(std::max(m_lo, log4_bases + 1), m_hi);
    constexpr double kLookupWeight = 1.0;
    std::vector<std::pair<unsigned, double>> evaluated;     // (m, objective)
    auto eval = [&](unsigned cand_m) -> double {
        for (const auto& entry : evaluated) if (entry.first == cand_m) return entry.second;

        const MEstimate est = estimate_m(seqs, kmer_len, cand_m);
        const double obj = (est.kmers > 0 ? est.bits / est.kmers : 0.0) + kLookupWeight * est.cost;

        evaluated.push_back({cand_m, obj});
        return obj;
    };
    unsigned best_m = start_m;
    double best = eval(start_m);

    // Walk each way while it improves, through one non-improving step, since the objective is not unimodal at short k (on chr1 at k = 15, m = 12 is worse than both 11 and 13, and 13 is the best).
    for (int dir : {-1, +1}) {
        double run_best = eval(start_m);
        int misses = 0;

        for (int m = int(start_m) + dir; m >= int(m_lo) && m <= int(m_hi); m += dir) {
            const double objective = eval(unsigned(m));

            if (objective < run_best) { run_best = objective; misses = 0; }
            else if (++misses == kMaxNonImprovingSteps) break;

            if (objective < best) { best = objective; best_m = unsigned(m); }
        }
    }

    return best_m;
}

}
