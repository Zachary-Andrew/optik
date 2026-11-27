#pragma once

/**
 * @file
 * @brief SpssCore, static exact k-mer membership and dictionary over a spectrum-preserving string set (SPSS), such as maximal unitigs.
 *
 * The design follows SSHash (Pibiri 2022).
 * SSHash buckets super-k-mer strings by minimizer, hashes the minimizers with a minimal perfect hash and confirms every lookup against the strings.
 * We assemble it from our own parts (phobic_mphf.hpp and Elias-Fano sequence starts) and depart from SSHash where noted below.
 *
 * @par Layout
 *   str_      Every sequence at 2 bits per base, most significant first inside each 64-bit word (kBasesPerStrWord bases per word).
 *             A k-mer at any offset then costs two word loads and one 128-bit shift.
 *   super-k-mers
 *             Maximal runs of consecutive k-mers in one sequence that share the same minimizer occurrence.
 *             A k-mer's minimizer is the m-mer whose canonical form hashes smallest, with ties going leftmost.
 *             As a consequence, a k-mer and its reverse complement select the same m-mer, mirrored.
 *   mphf_     Partitioned Mphf over the distinct canonical minimizers, which gives each one its bucket.
 *   rec_      One fixed-width record per bucket, read in a single cache-line access, holding
 *               q     the global offset of the minimizer occurrence
 *               near  1 if a k-mer containing that occurrence could reach past its sequence's ends
 *             A bucket with several super-k-mers (a repeated minimizer) sets kOverflowFlag instead and points into ovf_.
 *             That stream holds an Elias-gamma count followed by records in the same layout.
 *
 * @par Lookup
 * An exact lookup takes one MPHF access, one record access and one string access.
 * It starts from x's canonical minimizer and that minimizer's offset p in x.
 * For each record in the bucket, the only candidate start of x is s = q - p when x is forward, or s = q - (k - m - p) when it is reverse complemented.
 * Comparing the string confirms the candidate and also rejects keys the MPHF was never built on.
 * Near-boundary records additionally check that [s, s + k) stays inside one sequence.
 * If the minimum hash occurs at several offsets of x (the same canonical m-mer twice, say), every tied offset is tried.
 * SSHash scans the whole super-k-mer at this point, whereas storing the minimizer's offset lets us jump straight to the one candidate window.
 *
 * Heavy buckets and skew tables.
 * A minimizer shared by many super-k-mers (up to 164k in human chr1 repeats) would force a positive lookup to scan its whole list.
 * Buckets with more than kSkewMin records therefore put their k-mers into an open-addressing table, one per size class c covering buckets of up to 2^c records.
 * Each bit-packed slot of c + kFpBits bits holds a fingerprint and the index of the k-mer's record within the bucket.
 * Slots are probed linearly from a fast-range home slot at a load of about 0.9.
 * A probe chain thus stays within one or two cache lines, and the fingerprint skips nearly every wrong slot without touching the string.
 * Even so, the record found is verified against the string, so membership stays exact.
 * SSHash's skew index does the same job with one MPHF per size class.
 *
 * @par Sequence bounds
 * Sequence starts are stored as an Elias-Fano sequence, so finding the sequence that holds an offset takes one rank.
 *
 * Streaming with scan and scan_many.
 * Consecutive query k-mers usually continue the same string match.
 * So the next k-mer is first compared against the string at s + 1 (forward) or s - 1 (reverse), within the current sequence and without any hashing.
 * SSHash calls this an "extension".
 */

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bits.hpp"
#include "kmer.hpp"
#include "par_sort.hpp"
#include "phase_timer.hpp"
#include "phobic_mphf.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

class SpssCore {
public:
    static constexpr uint64_t kNotFound = ~uint64_t(0);

    struct Stats { uint64_t kmers = 0, bases = 0, seqs = 0, superkmers = 0,
                   buckets = 0, ovf_buckets = 0; };

    unsigned kmer_len() const noexcept { return kmer_len_; }
    unsigned minimizer_len() const noexcept { return minimizer_len_; }
    const Stats& stats() const noexcept { return stats_; }

    std::size_t bytes() const noexcept {
        std::size_t skew = 0;

        for (std::size_t words : skew_words_) skew += words * sizeof(uint64_t);

        return (str_words_ + rec_words_ + ovf_words_) * sizeof(uint64_t) + skew + mphf_.bytes() + ef_starts_.mem_bytes() + kScalarBytes;
    }

    /// Estimated size and lookup cost of an index at minimizer length m (estimate_m()).
    struct MEstimate {
        unsigned minimizer_len; double bits; double cost;
        double keys, records, overflow, skew, mphf;     // components, keys as a count and the others in bits
        double kmers;                                   // estimated k-mer positions
    };
    /// Estimated size and lookup cost of the index of @p seqs at minimizer length @p minimizer_len, without building it.
    ///
    /// All super-k-mers are found (the same sliding-window minimum as build()), but only those whose minimizer falls in a 1/4 hash sample of the minimizer space are kept.
    /// Each occurrence of a sampled minimizer is seen, so its bucket (records, k-mers) is exact, and the totals are scaled by 4.
    /// The bits cover the parts that depend on m, i.e. bucket records, overflow lists, skew tables and the MPHF (about kMphfBitsPerKey per minimizer), while the strings and sequence starts do not change with m.
    /// The overflow stream holds the records of all multi-record buckets, skew buckets included.
    /// The cost field is the mean number of records scanned per positive lookup in buckets below the skew threshold, and a skew bucket costs one table probe.
    static MEstimate estimate_m(const std::vector<std::string>& seqs, unsigned kmer_len, unsigned minimizer_len);

    /// Default minimizer length for @p seqs at @p kmer_len, the m that minimizes the estimated bits per k-mer + kLookupWeight x records scanned per lookup, as estimate_m() predicts them.
    ///
    /// The numbers below were measured when the rule was chosen.
    /// A smaller m shrinks the index but fills buckets, so lookups scan more records.
    /// On E. coli at k = 21, m = 10 gives the smallest index but queries 2.5x slower than m = 12.
    /// A weight of 1 bit per record scanned favours speed, since one record is about 24 ns of a search's DRAM chain.
    /// On chr1 at k = 21, m = 15 queries reads 37% faster than the size optimum m = 13 for 15% more space, and the rule picks m = 17 at k = 31 and m = 13 at k = 15.
    /// The index is exact and lossless at every m, where m only trades size against lookup work, so the sampled estimate can only make the default less than optimal, never wrong.
    /// Candidates are walked down and up from log4(bases) + 1, near the optimum on chr21 and chr1 at k = 15..31, through one non-improving step, within [max(5, k/2 - 4), k - 2], which keeps a window of at least three m-mers.
    /// k < 12 keeps k/2.
    static unsigned choose_m(const std::vector<std::string>& seqs, unsigned kmer_len);

    /// Build from ACGT sequences, where non-ACGT characters split a sequence.
    /// Returns false for invalid k/m or a failed MPHF build.
    ///
    /// Exact membership needs no uniqueness, but dense k-mer ids (kmer_id()) assume each canonical k-mer occurs once, as in an SPSS.
    /// @p nthreads sets the build workers (0 means one for the chunked passes, while the MPHF resolves 0 to half the CPUs).
    /// @p release (optional) is the caller's input vector, freed as soon as it is packed, since only the packing pass reads it.
    /// At short k the input is tens of millions of small strings (2.4 GB for chr1 at k = 15).
    bool build(const std::vector<std::string>& seqs, unsigned kmer_len, unsigned minimizer_len,
               unsigned nthreads = 0, std::vector<std::string>* release = nullptr) {
        if (kmer_len == 0 || kmer_len > kMaxK || minimizer_len == 0 || minimizer_len >= kmer_len || minimizer_len > kMaxMinimizerLen) return false;

        kmer_len_ = kmer_len; minimizer_len_ = minimizer_len; window_ = kmer_len - minimizer_len + 1;
        kmask_ = kmer_mask(kmer_len); mmask_ = kmer_mask(minimizer_len);
        stats_ = Stats{};
        has_distinct_ = false;
        str_.clear();
        PhaseTimer timer;                                  // OPTIK_TIMING=1

        // 1. Pack all ACGT runs of length >= k.
        uint64_t num_bases = 0;
        std::vector<uint64_t> starts;
        auto push_base = [&](uint8_t base) {
            if (num_bases % kBasesPerStrWord == 0) str_.push_back(0);

            // Bases are packed from each word's high bits down. A k-mer read across a word boundary then has its bases in order.
            str_.back() |= uint64_t(base) << (kWordBits - kBitsPerBase * (num_bases % kBasesPerStrWord + 1));
            ++num_bases;
        };

        for (const std::string& seq : seqs) {
            std::size_t run_start = 0, seq_len = seq.size();

            while (run_start < seq_len) {
                while (run_start < seq_len && base_encode(seq[run_start]) == kInvalidBase) ++run_start;

                std::size_t run_end = run_start;

                while (run_end < seq_len && base_encode(seq[run_end]) != kInvalidBase) ++run_end;

                // N splits a run because no stored k-mer may contain it, and runs shorter than k are dropped as k-mer free.
                if (run_end - run_start >= kmer_len) {
                    starts.push_back(num_bases);

                    for (std::size_t t = run_start; t < run_end; ++t) push_base(base_encode(seq[t]));

                    ++stats_.seqs;
                    stats_.kmers += (run_end - run_start) - kmer_len + 1;
                }

                run_start = run_end;
            }
        }

        // A final sentinel start makes every run's end the next run's start.
        starts.push_back(num_bases);
        stats_.bases = num_bases;

        if (release) std::vector<std::string>().swap(*release);   // seqs is unused from here

        str_.push_back(0); str_.push_back(0);          // padding for 2-word loads
        lo_bits_ = bits_for(num_bases + 1);
        // Sequence starts are coded as Elias-Fano, at about 2 + log2(bases / sequences) bits each instead of lo_bits_.
        // At short k the unitigs are short and numerous, and on chr21 at k = 15 fixed-width starts would take 5.7 of 31 bits per k-mer.
        ef_starts_.build(starts, num_bases + 1);
        n_starts_ = starts.size();
        adopt_owned();                                  // views for get_mer below
        timer.lap("spss: pack");

        // 2. Super-k-mers per sequence, from a sliding-window minimum over m-mers.
        // q is the absolute offset of the minimizer occurrence, and near says a k-mer anchored at q could reach past its sequence's ends.
        // s0 and len, the super-k-mer's first k-mer and k-mer count, are needed only during the build, since the skew tables enumerate exactly the k-mers of each record.
        // An SK is 16 bytes, the minimizer followed by one word packing q in 40 bits, q - s0 <= k - m in 6 bits, len <= k - m + 1 in 6 bits and near in 1 bit.
        constexpr unsigned kQBits = 40, kDeltaShift = 40, LEN_SHIFT = 46, NEAR_SHIFT = 52;
        constexpr uint64_t kField6Mask = 63u;

        struct SuperKmer {
            uint64_t key, packed;
            uint64_t min_pos() const noexcept { return packed & ((uint64_t(1) << kQBits) - 1); }
            uint64_t first_kmer_pos() const noexcept { return min_pos() - ((packed >> kDeltaShift) & kField6Mask); }
            uint32_t len() const noexcept { return static_cast<uint32_t>((packed >> LEN_SHIFT) & kField6Mask); }
            uint32_t near() const noexcept { return static_cast<uint32_t>((packed >> NEAR_SHIFT) & 1u); }
        };
        if (num_bases >= (uint64_t(1) << kQBits))
            throw std::length_error("SPSS build: input of 2^40 bases or more");

        std::vector<SuperKmer> sks;
        // Any k-mer containing the m-mer at q starts in [q-(k-m), q] and ends by q+k, so the occurrence is flagged if that span can leave [a, b).
        auto make_sk = [&](uint64_t key, uint64_t mer_pos, uint64_t seq_start, uint64_t seq_end, uint64_t first_kmer, uint32_t len) {
            const uint32_t near = (mer_pos < seq_start + (kmer_len_ - minimizer_len_)) || (mer_pos + kmer_len_ > seq_end);

            return SuperKmer{key, mer_pos | ((mer_pos - first_kmer) << kDeltaShift) | (uint64_t(len) << LEN_SHIFT) | (uint64_t(near) << NEAR_SHIFT)};
        };
        // One sequence's super-k-mers, appended to out with reused buffers.
        auto seq_sks = [&](std::size_t seq_idx, std::vector<SuperKmer>& out, std::vector<uint64_t>& hashes,
                           std::vector<uint64_t>& canon_mers, std::vector<uint64_t>& window_deque) {
            const uint64_t seq_start = starts[seq_idx], seq_end = starts[seq_idx + 1];
            const uint64_t num_mers = (seq_end - seq_start) - minimizer_len_ + 1;

            hashes.resize(num_mers); canon_mers.resize(num_mers); window_deque.resize(num_mers);

            for (uint64_t t = 0; t < num_mers; ++t) {
                const uint64_t mer = get_mer(seq_start + t, minimizer_len_);
                const uint64_t canon = std::min(mer, revcomp_packed(mer, minimizer_len_));

                canon_mers[t] = canon; hashes[t] = mhash(canon);
            }

            // Monotone deque of m-mer indices, with the leftmost minimum at the front.
            std::size_t deque_head = 0, deque_tail = 0;
            uint64_t cur_q = ~uint64_t(0); uint32_t cur_len = 0;
            const uint64_t num_kmers = (seq_end - seq_start) - kmer_len_ + 1;
            uint64_t next_in = 0;

            for (uint64_t s = 0; s < num_kmers; ++s) {
                while (next_in <= s + (kmer_len_ - minimizer_len_)) {             // push m-mers of window
                    while (deque_tail > deque_head && hashes[window_deque[deque_tail - 1]] > hashes[next_in]) --deque_tail;

                    window_deque[deque_tail++] = next_in++;
                }

                while (window_deque[deque_head] < s) ++deque_head;                        // expire

                const uint64_t min_idx = window_deque[deque_head];

                if (min_idx != cur_q) {
                    if (cur_len) out.push_back(make_sk(canon_mers[cur_q], seq_start + cur_q, seq_start, seq_end, seq_start + s - cur_len, cur_len));

                    cur_q = min_idx; cur_len = 0;
                }

                ++cur_len;
            }

            if (cur_len) out.push_back(make_sk(canon_mers[cur_q], seq_start + cur_q, seq_start, seq_end, seq_start + num_kmers - cur_len, cur_len));
        };

        {
            // Contiguous chunks of sequences, concatenated in chunk order, give the same vector a serial pass builds.
            const std::size_t nseq = starts.size() - 1;
            const unsigned num_threads = nthreads ? nthreads : 1;
            const std::size_t nch = std::min<std::size_t>(nseq, std::size_t(num_threads) * kChunksPerThread);
            std::vector<std::vector<SuperKmer>> part(nch);
#ifdef _OPENMP
            #pragma omp parallel num_threads(num_threads)
#endif
            {
                std::vector<uint64_t> hashes, canon_mers, window_deque;

#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (int64_t c = 0; c < int64_t(nch); ++c) {
                    const std::size_t seq_lo = nseq * std::size_t(c) / nch, seq_hi = nseq * std::size_t(c + 1) / nch;

                    part[c].reserve((starts[seq_hi] - starts[seq_lo]) / std::max(1u, window_ / kSuperKmerLenPerWindowDen) + kSuperKmerSlack);

                    // Each thread reuses its scratch vectors across sequences, and the inner loop never allocates.
                    for (std::size_t si = seq_lo; si < seq_hi; ++si) seq_sks(si, part[c], hashes, canon_mers, window_deque);
                }
            }
            std::size_t tot = 0;

            for (const auto& chunk : part) tot += chunk.size();

            sks.reserve(tot);

            // Freeing each chunk right after copying it holds peak memory near one copy of the super-k-mers instead of two.
            for (auto& chunk : part) { sks.insert(sks.end(), chunk.begin(), chunk.end()); std::vector<SuperKmer>().swap(chunk); }
        }
        stats_.superkmers = sks.size();
        timer.lap("spss: super-k-mers");

        // 3. Build the MPHF over the distinct minimizers and bucket the super-k-mers.
        // Since a minimizer occurrence ends one super-k-mer, (key, q) is unique, so the order is total and the parallel sort reproduces the serial one.
        // The sort runs in place, without par_sort's merge buffer of 1.3 GB on chr1 at k=15: one American-flag pass buckets by the top bits of the minimizer, and the buckets are then sorted in parallel.
        radix_sort_inplace(sks, [&](const SuperKmer& sk) { return sk.key; }, kBitsPerBase * minimizer_len_,
                           [](const SuperKmer& lhs, const SuperKmer& rhs) {
                               return lhs.key != rhs.key ? lhs.key < rhs.key : lhs.min_pos() < rhs.min_pos(); }, nthreads);
        timer.lap("spss: sort");
        // Chunks of whole buckets for the parallel passes below, which avoid per-bucket offset arrays next to sks (8 bytes per bucket each, and on HPRC 3 such arrays pushed the build over a 24 GiB cap).
        const std::vector<std::size_t> cuts = bucket_cuts(sks, std::size_t(nthreads ? nthreads : 1) * kCutsPerThread);
        const int64_t nch = static_cast<int64_t>(cuts.size()) - 1;

        distinct_ = count_distinct_buckets(sks, cuts);
        has_distinct_ = true;
        timer.lap("spss: distinct");
        std::size_t nkeys = 0;

        for (std::size_t i = 0; i < sks.size(); ++i) nkeys += i == 0 || sks[i].key != sks[i - 1].key;

        stats_.buckets = nkeys;
        {
            std::vector<uint64_t> keys;

            keys.reserve(nkeys);

            for (std::size_t i = 0; i < sks.size(); ++i)
                if (i == 0 || sks[i].key != sks[i - 1].key) keys.push_back(sks[i].key);

            // At alpha 0.94 about 6% of lookups hit the Elias-Fano remap, against about 15% at the size-optimal 0.85.
            // The higher load costs only about 0.02 bits per k-mer more.
            if (!keys.empty() && !mphf_.build(keys.data(), keys.size(), kMphfC, kMphfAlpha, nthreads)) return false;
        }
        timer.lap("spss: mphf");

        // The record layout holds a super-k-mer record or a flagged overflow offset.
        {
            const unsigned sk_bits = lo_bits_ + 1;
            uint64_t ovf_total = 0;

            for (std::size_t run_start = 0; run_start < sks.size();) {
                std::size_t run_end = run_start;

                while (run_end < sks.size() && sks[run_end].key == sks[run_start].key) ++run_end;

                if (run_end - run_start > 1) ovf_total += gamma_len(run_end - run_start - 1) + (run_end - run_start) * sk_bits;

                run_start = run_end;
            }

            rec_bits_ = std::max(sk_bits, bits_for(ovf_total + 1)) + 1;
        }
        rec_.assign(words_for_bits(nkeys * rec_bits_) + kPadWords, 0);
        ovf_.clear();
        // Each overflow list's bit offset is a prefix sum over the chunks, then over the buckets inside each chunk as it is written.
        // Neighbouring fields can share a word, so bits are OR-ed in atomically.
        auto ovf_len = [&](std::size_t count) -> uint64_t {
            return count > 1 ? gamma_len(count - 1) + count * (lo_bits_ + 1) : 0;
        };
        std::vector<uint64_t> ovf_at(nch + 1, 0);

#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads ? nthreads : 1) schedule(dynamic, 1)
#endif
        for (int64_t c = 0; c < nch; ++c) {
            uint64_t bits = 0;

            for (std::size_t run_start = cuts[c], run_end; run_start < cuts[c + 1]; run_start = run_end) {
                for (run_end = run_start + 1; run_end < cuts[c + 1] && sks[run_end].key == sks[run_start].key; ++run_end) {}

                bits += ovf_len(run_end - run_start);
            }

            ovf_at[c + 1] = bits;
        }

        // The first pass only sized each chunk's overflow lists, so this prefix sum turns the sizes into start offsets.
        for (int64_t c = 0; c < nch; ++c) ovf_at[c + 1] += ovf_at[c];

        const uint64_t ovf_nbits = ovf_at[nch];
        std::vector<uint64_t> ovf_bits_vec(words_for_bits(ovf_nbits) + kPadWords, 0);
        auto put_bits = [](uint64_t* words, uint64_t pos, uint64_t val, unsigned nbits) {
            if (!nbits) return;

            __atomic_fetch_or(&words[word_index(pos)], val << bit_offset(pos), __ATOMIC_RELAXED);

            if (bit_offset(pos) + nbits > kWordBits)
                __atomic_fetch_or(&words[word_index(pos) + 1], val >> (kWordBits - bit_offset(pos)), __ATOMIC_RELAXED);
        };
        uint64_t n_ovf = 0;

#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads ? nthreads : 1) schedule(dynamic, 1) reduction(+ : n_ovf)
#endif
        for (int64_t c = 0; c < nch; ++c) {
            uint64_t bit_pos = ovf_at[c];

            for (std::size_t run_start = cuts[c], run_end; run_start < cuts[c + 1]; run_start = run_end) {
                for (run_end = run_start + 1; run_end < cuts[c + 1] && sks[run_end].key == sks[run_start].key; ++run_end) {}

                const uint64_t bucket = mphf_(sks[run_start].key);

                // A bucket with a single super-k-mer stores its position in the record itself, with the low flag bit clear.
                if (run_end - run_start == 1) {
                    put_bits(rec_.data(), bucket * rec_bits_, pack(sks[run_start].min_pos(), sks[run_start].near()) << kEntryPayloadShift, rec_bits_);
                    continue;
                }

                // Otherwise the record holds the list's bit offset with the flag set, and the list is a gamma-coded count followed by one entry per super-k-mer.
                ++n_ovf;
                put_bits(rec_.data(), bucket * rec_bits_, (bit_pos << 1) | 1, rec_bits_);
                put_bits(ovf_bits_vec.data(), bit_pos, gamma_word(run_end - run_start - 1), gamma_len(run_end - run_start - 1));
                bit_pos += gamma_len(run_end - run_start - 1);

                for (std::size_t t = run_start; t < run_end; ++t) {
                    put_bits(ovf_bits_vec.data(), bit_pos, pack(sks[t].min_pos(), sks[t].near()), lo_bits_ + 1);
                    bit_pos += lo_bits_ + 1;
                }
            }
        }

        stats_.ovf_buckets += n_ovf;
        ovf_ = std::move(ovf_bits_vec);
        // The two zero words are the kPadWords padding, since get_bits() loads two words at a time.
        ovf_.push_back(0); ovf_.push_back(0);
        adopt_owned();
        timer.lap("spss: records");
        build_skew(sks, nthreads);
        timer.lap("spss: skew");
        return true;
    }

private:
    /// Distinct canonical k-mers of the sorted super-k-mers, one bucket at a time.
    /// Every occurrence of a canonical k-mer on either strand has the same minimizer, since the m-mer hash is a bijection and the minimum hash names one m-mer value, and so lies in that minimizer's bucket.
    /// A small sort per bucket replaces a full lookup per stored k-mer.
    template <typename SKV>
    uint64_t count_distinct_buckets(const SKV& sks, const std::vector<std::size_t>& cuts) const {
        const int64_t nch = static_cast<int64_t>(cuts.size()) - 1;
        uint64_t distinct_count = 0;
#ifdef _OPENMP
        #pragma omp parallel reduction(+ : distinct_count)
#endif
        {
            std::vector<uint64_t> buf;

#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1)
#endif
            for (int64_t c = 0; c < nch; ++c)
                for (std::size_t group_start = cuts[c], group_end; group_start < cuts[c + 1]; group_start = group_end) {
                    buf.clear();

                    // A group is the super-k-mers of one minimizer bucket, and all of their canonical k-mers are collected for sorting.
                    for (group_end = group_start; group_end < cuts[c + 1] && sks[group_end].key == sks[group_start].key; ++group_end)
                        for (uint32_t u = 0; u < sks[group_end].len(); ++u) {
                            const uint64_t kmer = get_mer(sks[group_end].first_kmer_pos() + u, kmer_len_), rev_comp = revcomp_packed(kmer, kmer_len_);

                            buf.push_back(kmer < rev_comp ? kmer : rev_comp);
                        }

                    if (buf.size() > 1) std::sort(buf.begin(), buf.end());

                    // After the sort equal k-mers are adjacent, which makes the number of run starts the distinct count.
                    for (std::size_t i = 0; i < buf.size(); ++i) distinct_count += i == 0 || buf[i] != buf[i - 1];
                }
        }
        return distinct_count;
    }

    /// Up to @p nch + 1 cut points into the key-sorted @p sks, each at the start of a bucket, which split it into chunks of whole buckets for the parallel passes.
    template <typename SKV>
    static std::vector<std::size_t> bucket_cuts(const SKV& sks, std::size_t nch) {
        const std::size_t count = sks.size();
        std::vector<std::size_t> cut{0};

        for (std::size_t c = 1; c < nch; ++c) {
            std::size_t cut_pos = std::max(cut.back(), count / nch * c + count % nch * c / nch);

            while (cut_pos > 0 && cut_pos < count && sks[cut_pos].key == sks[cut_pos - 1].key) ++cut_pos;

            if (cut_pos > cut.back() && cut_pos < count) cut.push_back(cut_pos);
        }

        cut.push_back(count);
        return cut;
    }

    /// Skew tables over the k-mers of every bucket with more than kSkewMin records.
    /// There is one table per size class, filled sequentially per class with the classes in parallel.
    /// Bucket b of class c inserts each distinct canonical k-mer of its records once, as (fingerprint, t) with t its earliest record.
    /// A lookup returns the first matching entry along the probe chain, which is that same record, so dropping the repeats changes no answer, and it keeps a k-mer repeated n times (low-complexity runs) from building an n-long probe chain with quadratic insert time.
    template <typename SKV>
    void build_skew(const SKV& sks, unsigned nthreads) {
        skew_own_.assign(kMaxClass + 1, {});
        skew_cap_.assign(kMaxClass + 1, 0);
        skew_tables_.assign(kMaxClass + 1, nullptr);
        struct Range { std::size_t begin, end; uint64_t distinct; };     // n = distinct k-mers
        std::vector<std::vector<Range>> by_class(kMaxClass + 1);
        std::vector<uint64_t> entries(kMaxClass + 1, 0);

        // A bucket's k-mers as (canonical, record, order) triples, enumerated by record and then offset, with each canonical k-mer's repeats dropped.
        struct KmerOcc { uint64_t kmer; uint32_t rec, order; };
        auto distinct_of = [&](std::size_t first, std::size_t last, std::vector<KmerOcc>& out) {
            out.clear();

            for (std::size_t t = first; t < last; ++t)
                for (uint32_t u = 0; u < sks[t].len(); ++u) {
                    const uint64_t kmer = get_mer(sks[t].first_kmer_pos() + u, kmer_len_), rev_comp = revcomp_packed(kmer, kmer_len_);

                    out.push_back({kmer < rev_comp ? kmer : rev_comp, static_cast<uint32_t>(t - first), static_cast<uint32_t>(out.size())});
                }

            // Sorting by k-mer and then by order puts each k-mer's earliest occurrence first, so unique() keeps that one.
            // A final sort by order restores the enumeration order.
            std::sort(out.begin(), out.end(), [](const KmerOcc& lhs, const KmerOcc& rhs) { return lhs.kmer != rhs.kmer ? lhs.kmer < rhs.kmer : lhs.order < rhs.order; });
            out.erase(std::unique(out.begin(), out.end(), [](const KmerOcc& lhs, const KmerOcc& rhs) { return lhs.kmer == rhs.kmer; }), out.end());
            std::sort(out.begin(), out.end(), [](const KmerOcc& lhs, const KmerOcc& rhs) { return lhs.order < rhs.order; });
        };
        // Only buckets with more than kSkewMin records get skew entries, so the light ones skip this work.
        std::vector<Range> heavy;

        for (std::size_t run_start = 0; run_start < sks.size();) {
            std::size_t run_end = run_start;

            while (run_end < sks.size() && sks[run_end].key == sks[run_start].key) ++run_end;

            if (run_end - run_start > kSkewMin) heavy.push_back({run_start, run_end, 0});

            run_start = run_end;
        }
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            std::vector<KmerOcc> kmers;

#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1)
#endif
            for (std::size_t h = 0; h < heavy.size(); ++h) {
                distinct_of(heavy[h].begin, heavy[h].end, kmers);
                heavy[h].distinct = kmers.size();
            }
        }
        // Heavy buckets are grouped by size class, and each class shares one table sized for all its distinct k-mers.
        for (const Range& range : heavy) {
            const unsigned size_class = cls_of(range.end - range.begin);

            by_class[size_class].push_back(range);
            entries[size_class] += range.distinct;
        }

        std::vector<unsigned> todo;

        skew_entries_ = 0;

        for (unsigned c = 0; c <= kMaxClass; ++c) {
            skew_entries_ += entries[c];

            if (!entries[c]) continue;

            // A load of about 0.9 with fast-range home slots, so no power-of-two rounding.
            skew_cap_[c] = entries[c] + entries[c] / kSkewSlackDiv + 1;
            skew_own_[c].assign(words_for_bits(skew_cap_[c] * (c + kFpBits)) + kPadWords, 0);
            todo.push_back(c);
        }

        (void)nthreads;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1)
#endif
        for (std::size_t ti = 0; ti < todo.size(); ++ti) {
            const unsigned size_class = todo[ti];
            const unsigned slot_bits = size_class + kFpBits;
            std::vector<uint64_t>& table_words = skew_own_[size_class];
            const uint64_t cap = skew_cap_[size_class];
            std::vector<KmerOcc> kmers;

            for (const Range& range : by_class[size_class]) {
                distinct_of(range.begin, range.end, kmers);

                for (const KmerOcc& entry : kmers) {
                    const uint64_t hash = skew_hash(entry.kmer);
                    // The fingerprint occupies the high bits above the record index, which leaves the all-zero word free to mean an empty slot.
                    const uint64_t val = (skew_fp(hash) << size_class) | entry.rec;

                    // Collisions probe linearly and wrap at the end.
                    for (uint64_t slot = home(hash, cap);; slot = (slot + 1 == cap) ? 0 : slot + 1) {
                        if (get_bits(table_words.data(), slot * slot_bits, slot_bits) == 0) { set_bits(table_words.data(), slot * slot_bits, val, slot_bits); break; }
                    }
                }
            }
        }

        skew_words_.assign(kMaxClass + 1, 0);

        for (unsigned c = 0; c <= kMaxClass; ++c) {
            // Lookups go through these views, which point at the owned vectors after build() and into the mapping after load().
            skew_tables_[c] = skew_own_[c].data();
            skew_words_[c] = skew_own_[c].size();
        }
    }

public:

    /// Exact membership of k-mer @p kmer, given in forward packing, either strand.
    /// Returns the string offset of the matching occurrence, or kNotFound.
    uint64_t lookup(uint64_t kmer) const noexcept {
        if (stats_.buckets == 0) return kNotFound;

        MinInfo min_info = minimizer(kmer);

        return probe(kmer, min_info, nullptr);
    }

    bool contains(uint64_t kmer) const noexcept { return lookup(kmer) != kNotFound; }

    // Dictionary operations, SSHash's Lookup and Access.
    // String i occupies [start_i, start_i+1) and holds len_i - k + 1 k-mers, so the k-mer at offset s of string i has the dense identifier s - (k - 1) * i.
    // Identifiers cover [0, n) and are a bijection onto the indexed k-mers when the input holds each k-mer once (unitigs / an SPSS, as SSHash also requires).

    /// Dense identifier in [0, num_kmers) of k-mer @p kmer (either strand), or kNotFound.
    /// *@p rev is set when @p kmer is the reverse complement of the stored k-mer.
    uint64_t kmer_id(uint64_t kmer, bool* rev = nullptr) const noexcept {
        if (stats_.buckets == 0) return kNotFound;

        const uint64_t offset = probe(kmer, minimizer(kmer), rev);

        if (offset == kNotFound) return kNotFound;

        return offset - static_cast<uint64_t>(kmer_len_ - 1) * seq_index(offset);
    }

    /// The stored k-mer (forward packing, as it appears in the string) with identifier @p kmer_ident < num_kmers().
    uint64_t access(uint64_t kmer_ident) const noexcept {
        // The largest i with start_i - (k-1) i <= id.
        std::size_t lo_idx = 0, hi_idx = n_starts_ - 1;

        while (hi_idx - lo_idx > 1) {
            const std::size_t mid = lo_idx + (hi_idx - lo_idx) / kBisect;
            const uint64_t before = start_at(mid) - static_cast<uint64_t>(kmer_len_ - 1) * mid;

            if (before <= kmer_ident) lo_idx = mid; else hi_idx = mid;
        }

        return get_mer(kmer_ident + static_cast<uint64_t>(kmer_len_ - 1) * lo_idx, kmer_len_);
    }

    uint64_t num_kmers() const noexcept { return stats_.kmers; }
    /// k-mers indexed through skew tables by the last build(), or 0 after load.
    uint64_t skew_entries() const noexcept { return skew_entries_; }

    /// Distinct canonical k-mers, where a stored k-mer counts once, at the first occurrence a lookup returns.
    /// The count is known from build().
    /// After a load it walks the sequences directly, so no id has to be inverted with a binary search over the starts.
    uint64_t count_distinct() const noexcept {
        if (has_distinct_) return distinct_;

        if (stats_.buckets == 0) return 0;

        const int64_t num_seqs = static_cast<int64_t>(n_starts_) - 1;
        uint64_t distinct_count = 0;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, kDistinctSchedChunk) reduction(+ : distinct_count)
#endif
        for (int64_t si = 0; si < num_seqs; ++si) {
            const uint64_t seq_start = start_at(si);
            const uint64_t seq_end = start_at(si + 1);

            for (uint64_t s = seq_start; s + kmer_len_ <= seq_end; ++s)
                distinct_count += probe(get_mer(s, kmer_len_), minimizer(get_mer(s, kmer_len_)), nullptr) == s;
        }

        return distinct_count;
    }

    /// @p func(canonical k-mer) for every stored k-mer, in sequence order.
    template <typename F> void for_each_kmer(F&& func) const {
        for (std::size_t si = 0; si + 1 < n_starts_; ++si) {
            const uint64_t seq_start = start_at(si);
            const uint64_t seq_end = start_at(si + 1);

            for (uint64_t s = seq_start; s + kmer_len_ <= seq_end; ++s) {
                const uint64_t kmer = get_mer(s, kmer_len_), rev = revcomp_packed(kmer, kmer_len_);

                func(kmer < rev ? kmer : rev);
            }
        }
    }

    /// Streaming lookup of every k-mer of the query records, counted as SSHash's query counts them.
    ///
    /// Long records use one cursor.
    /// An extension (the next string k-mer continues the previous match) is one string compare and hashes nothing.
    /// A full search needs the minimizer of the current k-mer.
    /// The canonical m-mers roll per base (cheap), and their hashes are computed lazily, only for window positions not hashed yet, so a run of extensions costs no hashing and consecutive searches hash one m-mer each.
    /// The argmin and tie set come from a branch-free scan over the window's W hashes, since a monotone queue mispredicts on random hashes.
    ///
    /// Tiny records (<= kGroupMaxKmers k-mers, such as single k-mer queries) have no stream to extend along, so their k-mers are looked up in groups of kGroup with staged prefetching (all pilot lines, then all bucket records, then all string lines), and the three dependent misses of independent lookups overlap instead of serialising.
    struct ScanStats { uint64_t valid = 0, positive = 0, searches = 0, extensions = 0; };
    static constexpr unsigned kGroup = 16, kGroupMaxKmers = 4;

    void scan(const char* seq, std::size_t len, ScanStats& out) const noexcept;

    void scan_many(const char* const* seqs, const std::size_t* lens, std::size_t nseq,
                   ScanStats& out) const noexcept {
        GroupBuf group;
        // Longer records wait here for a lane (local, since scan_many runs concurrently on one shared index).
        std::vector<std::pair<const char*, std::size_t>> long_recs;

        long_recs.reserve(nseq);

        for (std::size_t r = 0; r < nseq; ++r) {
            if (lens[r] >= kmer_len_ && lens[r] - kmer_len_ + 1 <= kGroupMaxKmers) {
                // Collect every valid k-mer of the tiny record.
                uint64_t fwd = 0; std::size_t run_len = 0;

                for (std::size_t i = 0; i < lens[r]; ++i) {
                    const uint8_t base = base_encode(seqs[r][i]);

                    if (base == kInvalidBase) { run_len = 0; fwd = 0; continue; }

                    fwd = ((fwd << kBitsPerBase) | base) & kmask_;

                    if (++run_len >= kmer_len_) {
                        group.kmers[group.count++] = fwd;

                        if (group.count == kGroup) flush_group(group, out);
                    }
                }
            } else {
                long_recs.push_back({seqs[r], lens[r]});
            }
        }

        flush_group(group, out);

        if (long_recs.empty()) return;

        if (stats_.buckets == 0) {                    // empty index, which counts but never hits
            for (const auto& rec : long_recs) scan_one(rec.first, rec.second, out);
        } else {
            scan_lanes(long_recs.data(), long_recs.size(), out);
        }
    }

    /// Reference single-cursor scan, one record and one dependent miss chain at a time.
    /// scan_many runs records through scan_lanes instead, and the stress tests check that both give identical counts.
    void scan_reference(const char* seq, std::size_t len, ScanStats& out) const noexcept;

private:
    struct GroupBuf { uint64_t kmers[kGroup]; unsigned count = 0; };

    /// Look up the buffered independent k-mers with staged prefetching.
    void flush_group(GroupBuf& group, ScanStats& out) const noexcept;

    // Interleaved streaming scan.
    // A full search is a chain of dependent misses, from MPHF pilot to bucket record to overflow list to string, so one cursor per thread waits out one DRAM miss at a time.
    // scan_lanes advances kLanes records round-robin, and each search is split into stages that prefetch the line the next stage reads, so up to kLanes misses are in flight.
    // Per record the work and counts match scan_one's, with the same extension rule, bucket cache and probe.
    // Only the order across records changes, and the counts are sums.
    static constexpr unsigned kLanes = 16;
    static constexpr uint64_t kProbeChunk = 32;         // overflow records decoded per chunk
    /// Largest number of m-mers in a k-mer (k - m + 1 with k <= 32, m >= 1), also the size of the m-mer rings below, so it must be a power of two.
    static constexpr unsigned kMaxWindowMmers = 32;
    /// A minimum hash that ties is kept at both its first and its last position in the window.
    static constexpr uint32_t kTiedMinPositions = 2;
    static constexpr uint64_t kRingMask = kMaxWindowMmers - 1;
    /// A lane moves through these stages, each named after what it has just prefetched, while kStageNone has nothing pending.
    static constexpr unsigned kStageNone = 0;
    static constexpr unsigned kStagePilot = 1;
    static constexpr unsigned kStageRecord = 2;
    static constexpr unsigned kStageOverflow = 3;
    static constexpr unsigned kStageString = 4;

    /// A k-mer's canonical minimizer and the offsets in the k-mer where its minimum hash occurs, npos in number.
    struct MinInfo { uint64_t key; uint32_t npos; uint8_t pos[kMaxWindowMmers]; };
    struct Lane {
        const char* seq = nullptr; std::size_t len = 0, pos = 0;
        uint64_t fwd = 0, fwd_mer = 0, rev_mer = 0; std::size_t valid = 0;
        uint64_t midx = 0, hashed = 0, best = 0, best_lo = 0, best_hi = 0; bool have_min = false;
        uint64_t cur = kNotFound; bool cur_rev = false; uint64_t seq_lo = 0, seq_hi = 0;
        uint64_t cached_key = ~uint64_t(0), cached_entry = 0; bool cached_ok = false;
        // Pending search stage, where 0 is none, 1 pilot prefetched, 2 record prefetched, 3 overflow head prefetched and 4 string prefetched.
        unsigned stage = kStageNone; bool fresh = false;
        MinInfo min_info; uint64_t kmer = 0, entry = 0;
        uint64_t rcn[kMaxWindowMmers], ring_hashes[kMaxWindowMmers];   // rings of canonical m-mers and hashes
    };
    void prefetch_window(uint64_t mer_pos) const noexcept;
    void lane_reset(Lane& lane, const char* seq, std::size_t len) const noexcept;
    /// Run lane @p lane's record forward until it needs a full search whose data is not cached (returning true, with a stage pending and its first line prefetched) or the record ends (false).
    /// Mirrors scan_one's loop.
    bool lane_run(Lane& lane, ScanStats& out) const noexcept;
    /// One stage of lane @p lane's pending search, true while still pending.
    bool lane_stage(Lane& lane, ScanStats& out) const noexcept;
    void scan_lanes(const std::pair<const char*, std::size_t>* recs, std::size_t count,
                    ScanStats& out) const noexcept {
        Lane lanes[kLanes];
        std::size_t next = 0;
        unsigned live = 0;
        bool lane_on[kLanes] = {};
        // Start a lane on the next record that needs a search, or return false when no records are left.
        auto feed = [&](Lane& lane) -> bool {
            while (next < count) {
                lane_reset(lane, recs[next].first, recs[next].second);
                ++next;

                if (lane_run(lane, out)) return true;
            }

            return false;
        };

        for (unsigned l = 0; l < kLanes; ++l) { lane_on[l] = feed(lanes[l]); live += lane_on[l]; }

        while (live) {
            for (unsigned l = 0; l < kLanes; ++l) {
                if (!lane_on[l]) continue;

                Lane& lane = lanes[l];

                if (lane_stage(lane, out)) continue;           // still pending

                if (lane_run(lane, out)) continue;             // next search staged

                if (!feed(lane)) { lane_on[l] = false; --live; }
            }
        }
    }

    /// Single-cursor streaming scan of one record (see scan_many()).
    void scan_one(const char* seq, std::size_t len, ScanStats& out) const noexcept;

public:
    /// Serialize the index, as load() reads it back.
    void save(std::vector<uint8_t>& out) const;
    /// Read a save() blob, returning the bytes consumed, or 0 for a malformed one.
    /// With borrow = true every array is used in place, the MPHF's and the Elias-Fano starts' included (zero-copy, so `buf` must stay mapped and 8-byte aligned for the lifetime of this object), and otherwise they are copied.
    std::size_t load(const uint8_t* buf, std::size_t buf_len, bool borrow = false);

    /// Kept as a member for existing callers, while kmer.hpp's reverse_complement() does the work.
    static constexpr uint64_t revcomp_packed(uint64_t kmer, unsigned width) noexcept {
        return reverse_complement(kmer, width);
    }

private:

    static constexpr unsigned bits_for(uint64_t limit) noexcept {   // bits to represent [0, v)
        unsigned bits = 0;

        while ((uint64_t(1) << bits) < limit) ++bits;

        return bits ? bits : 1;
    }
    /// MurmurHash3 fmix64, the hash that orders canonical m-mers.
    /// Index files depend on it.
    static constexpr uint64_t mhash(uint64_t value) noexcept {
        return mix64(value, 0);
    }
    /// A record stores q above one near bit.
    static constexpr unsigned kNearBits = 1;
    static constexpr uint64_t pack(uint64_t mer_pos, uint32_t near) noexcept {
        return (mer_pos << kNearBits) | near;
    }
    /// A bucket entry is either a record payload or an overflow-stream offset, both shifted above kOverflowFlag.
    static constexpr uint64_t kOverflowFlag = 1;
    static constexpr unsigned kEntryPayloadShift = 1;
    static constexpr bool entry_overflows(uint64_t entry) noexcept { return (entry & kOverflowFlag) != 0; }
    static constexpr uint64_t entry_ovf_offset(uint64_t entry) noexcept { return entry >> kEntryPayloadShift; }
    static constexpr uint64_t entry_mer_pos(uint64_t entry) noexcept { return entry >> (kEntryPayloadShift + kNearBits); }
    /// Bases per word of the packed string, most significant first.
    static constexpr unsigned kBasesPerStrWord = kWordBits / kBitsPerBase;
    /// Elias-gamma code (Elias 1975) of n >= 1 in the LSB-first bit streams.
    /// That is L = floor(log2 n) zero bits, a one, and then the L low bits of n, for 2L + 1 bits in all.
    /// It codes overflow-bucket record counts, mostly 2-4 records in 1-3 bits, where a fixed field took 15-17 bits.
    static constexpr unsigned gamma_len(uint64_t value) noexcept {
        return floor_log2(value) + (floor_log2(value) + 1u);   // the zero prefix, then the value's own bits
    }
    static constexpr uint64_t gamma_word(uint64_t value) noexcept {
        const unsigned log_len = floor_log2(value);

        return (uint64_t(1) << log_len) | ((value & ((uint64_t(1) << log_len) - 1)) << (log_len + 1));
    }
    /// Record count of the overflow entry at bit offset @p off, and the width of its count field.
    uint64_t ovf_count(uint64_t off, unsigned& width) const noexcept {
        const uint64_t word = get_bits(ovf_view_, off, kWordBits);
        const unsigned zeros = static_cast<unsigned>(__builtin_ctzll(word));

        width = zeros + (zeros + 1u);   // Elias-gamma length, the zero prefix and then the value's bits
        return ((uint64_t(1) << zeros) | ((word >> (zeros + 1)) & ((uint64_t(1) << zeros) - 1))) + 1;   // n = cnt - 1
    }
    static uint64_t get_bits(const uint64_t* words, uint64_t pos, unsigned nbits) noexcept {
        const uint64_t lo_word = words[word_index(pos)], hi_word = words[word_index(pos) + 1];
        const unsigned shift = bit_offset(pos);
        uint64_t bits = shift ? ((lo_word >> shift) | (hi_word << (kWordBits - shift))) : lo_word;

        return nbits == kWordBits ? bits : (bits & ((uint64_t(1) << nbits) - 1));
    }
    /// Reads the @p width-mer (@p width <= 32) starting at global base offset @p offset.
    uint64_t get_mer(uint64_t offset, unsigned width) const noexcept {
        constexpr unsigned kPairBits = 2 * kWordBits;
        const uint64_t word_idx = offset / kBasesPerStrWord;
        const __uint128_t pair = ((__uint128_t)str_view_[word_idx] << kWordBits) | str_view_[word_idx + 1];

        return static_cast<uint64_t>((pair << (kBitsPerBase * (offset % kBasesPerStrWord))) >> (kPairBits - kBitsPerBase * width));
    }

    /// The minimizer of k-mer @p kmer and its tied offsets.
    MinInfo minimizer(uint64_t kmer) const noexcept {
        const uint64_t rev_comp = revcomp_packed(kmer, kmer_len_);
        const unsigned window = kmer_len_ - minimizer_len_ + 1;
        uint64_t hashes[kMaxWindowMmers] = {}, canon_mers[kMaxWindowMmers] = {};

        for (unsigned p = 0; p < window; ++p) {
            const uint64_t mer = (kmer >> (kBitsPerBase * (kmer_len_ - minimizer_len_ - p))) & mmask_;
            const uint64_t mer_rc = (rev_comp >> (kBitsPerBase * p)) & mmask_;     // rc of the m-mer at p

            canon_mers[p] = mer < mer_rc ? mer : mer_rc;
            hashes[p] = mhash(canon_mers[p]);
        }

        uint64_t best = ~uint64_t(0); unsigned best_pos = 0;

        for (unsigned p = 0; p < window; ++p) {
            const bool is_less = hashes[p] < best;

            best = is_less ? hashes[p] : best;
            best_pos = is_less ? p : best_pos;
        }

        MinInfo min_info;

        min_info.key = canon_mers[best_pos];
        min_info.npos = 0;

        for (unsigned p = best_pos; p < window; ++p) {
            min_info.pos[min_info.npos] = static_cast<uint8_t>(p);
            min_info.npos += hashes[p] == best;
        }

        return min_info;
    }

    /// Try every candidate window of @p kmer's bucket and return the string offset, or kNotFound.
    uint64_t probe(uint64_t kmer, const MinInfo& min_info, bool* rev_out) const noexcept {
        return probe_entry(kmer, min_info, bucket_entry(min_info.key), rev_out);
    }
    /// True iff bucket record @p entry belongs to minimizer @p key.
    /// Every occurrence stored in a bucket spells its canonical minimizer, so one m-mer compare rejects keys the MPHF maps to a foreign bucket.
    /// The string doubles as the key store and needs no fingerprint bits.
    bool entry_has_key(uint64_t entry, uint64_t key) const noexcept {
        uint64_t mer_pos;

        if (!entry_overflows(entry)) mer_pos = entry_mer_pos(entry);
        else {
            unsigned count_width;

            ovf_count(entry >> 1, count_width);
            mer_pos = get_bits(ovf_view_, (entry >> 1) + count_width, lo_bits_ + 1) >> 1;
        }

        const uint64_t mer = get_mer(mer_pos, minimizer_len_);

        return (mer < revcomp_packed(mer, minimizer_len_) ? mer : revcomp_packed(mer, minimizer_len_)) == key;
    }
    /// The record of minimizer @p key's bucket.
    uint64_t bucket_entry(uint64_t key) const noexcept {
        return get_bits(rec_view_, mphf_(key) * rec_bits_, rec_bits_);
    }
    /// probe() with the bucket record @p entry already read.
    uint64_t probe_entry(uint64_t kmer, const MinInfo& min_info, uint64_t entry,
                         bool* rev_out) const noexcept {
        const uint64_t rev_comp = revcomp_packed(kmer, kmer_len_);
        auto try_rec = [&](uint64_t rec) -> uint64_t {
            const uint64_t mer_pos = rec >> 1;
            const bool near = rec & 1;
            // An exact string match at s is a hit unless the window crosses into a neighbouring sequence, which only near-boundary occurrences can do.
            // The bounds of q's sequence, computed once per record and checked against every candidate offset of the record.
            uint64_t seq_lo = 0, seq_hi = 0;
            bool have = false;
            auto inside = [&](uint64_t start) {
                if (!near) return true;

                if (!have) { seq_bounds(mer_pos, seq_lo, seq_hi); have = true; }

                return start >= seq_lo && start + kmer_len_ <= seq_hi;
            };
            // The candidate offsets are the same set for x and rc(x), since the minimizer positions mirror, so the smallest matching offset, rather than the first one tried, gives the answer.
            // That makes it, and with it --ids, strand-independent when a record holds a k-mer twice (a hairpin with x and rc(x) a few bases apart).
            uint64_t best = kNotFound;
            bool best_rev = false;

            for (unsigned t = 0; t < min_info.npos; ++t) {
                const uint64_t mer_off = min_info.pos[t];

                if (mer_pos >= mer_off) {                                    // x forward
                    const uint64_t start = mer_pos - mer_off;

                    if (start < best && get_mer(start, kmer_len_) == kmer && inside(start)) { best = start; best_rev = false; }
                }

                const uint64_t rev_off = kmer_len_ - minimizer_len_ - mer_off;                 // x reverse

                if (mer_pos >= rev_off) {
                    const uint64_t start = mer_pos - rev_off;

                    if (start < best && get_mer(start, kmer_len_) == rev_comp && inside(start)) { best = start; best_rev = true; }
                }
            }

            if (best != kNotFound && rev_out) *rev_out = best_rev;

            return best;
        };

        if (!(entry & 1)) return try_rec(entry >> 1);

        uint64_t off = entry >> 1;
        unsigned count_width;
        const uint64_t cnt = ovf_count(off, count_width);

        off += count_width;
        const unsigned skb = lo_bits_ + 1;

        if (cnt > kSkewMin) {
            // In a heavy bucket the k-mer's record index comes from its size class's open-addressing table (fingerprint-filtered linear probing).
            const unsigned size_class = cls_of(cnt);
            const uint64_t cap = skew_cap_[size_class];
            const uint64_t* table = skew_tables_[size_class];
            const unsigned slot_bits = size_class + kFpBits;
            const uint64_t hash = skew_hash(kmer < rev_comp ? kmer : rev_comp);
            const uint64_t fingerprint = skew_fp(hash), imask = (uint64_t(1) << size_class) - 1;

            for (uint64_t slot = home(hash, cap);; slot = (slot + 1 == cap) ? 0 : slot + 1) {
                const uint64_t slot_val = get_bits(table, slot * slot_bits, slot_bits);

                if (slot_val == 0) return kNotFound;                  // end of the chain

                if ((slot_val >> size_class) == fingerprint && (slot_val & imask) < cnt) {
                    const uint64_t offset = try_rec(get_bits(ovf_view_, off + (slot_val & imask) * skb, skb));

                    if (offset != kNotFound) return offset;
                }
            }
        }

        // A bucket's records point into the string at unrelated places, so checking them one after another serialises one cache miss each.
        // Records are decoded a chunk at a time and, when there is more than one, the string lines of the whole chunk are prefetched before the first is compared, so their misses overlap (in the same order, with the same answer, the first record that holds the k-mer).
        constexpr uint64_t chunk = kProbeChunk;
        uint64_t recs[chunk];

        for (uint64_t base = 0; base < cnt; base += chunk) {
            const uint64_t count = std::min(chunk, cnt - base);

            for (uint64_t i = 0; i < count; ++i) recs[i] = get_bits(ovf_view_, off + (base + i) * skb, skb);

            if (count > 1)
                for (uint64_t i = 0; i < count; ++i) {
                    const uint64_t mer_pos = recs[i] >> 1;

                    __builtin_prefetch(&str_view_[(mer_pos > kmer_len_ ? mer_pos - kmer_len_ : 0) / kBasesPerStrWord]);
                    __builtin_prefetch(&str_view_[(mer_pos + kmer_len_) / kBasesPerStrWord]);
                }

            for (uint64_t i = 0; i < count; ++i) {
                const uint64_t offset = try_rec(recs[i]);

                if (offset != kNotFound) return offset;
            }
        }

        return kNotFound;
    }

    unsigned kmer_len_ = 0, minimizer_len_ = 0, window_ = 0;
    uint64_t kmask_ = 0, mmask_ = 0;
    unsigned lo_bits_ = 1, rec_bits_ = kInitialRecBits;
    std::size_t n_starts_ = 0;     // number of packed sequence starts (seqs + 1)
    // Sequence starts as an Elias-Fano sequence.
    EfSeq ef_starts_;
    uint64_t start_at(std::size_t idx) const noexcept { return ef_starts_[idx]; }

    /// Index of the sequence containing string offset @p pos.
    std::size_t seq_index(uint64_t pos) const noexcept {
        return static_cast<std::size_t>(ef_starts_.count_le(pos)) - 1;
    }
    /// [@p seq_lo, @p seq_hi) of the sequence holding base offset @p pos, one Elias-Fano rank and two accesses.
    void seq_bounds(uint64_t pos, uint64_t& seq_lo, uint64_t& seq_hi) const noexcept {
        const std::size_t seq_idx = seq_index(pos);

        seq_lo = start_at(seq_idx);
        seq_hi = start_at(seq_idx + 1);
    }
    // Skew tables for heavy buckets.
    // These three parameters are stored in the index file and checked at load.
    static constexpr uint64_t kSkewMin = 64;        // records above which a bucket gets a skew table
    static constexpr unsigned kFpBits = 5;          // fingerprint bits per skew slot
    static constexpr unsigned kMaxClass = 40;       // largest size class, buckets of up to 2^40 records
    // MPHF build parameters, plus the size per key that estimate_m() assumes.
    static constexpr double kMphfC = 6.0, kMphfAlpha = 0.94;
    static constexpr double kMphfBitsPerKey = 5.1;
    static constexpr std::size_t kScalarBytes = 64;   // rough size of the members, for bytes()
    /// Skew tables get kSkewLoadNum / kSkewLoadDen slots per entry (a load of 0.9), and their capacity adds 1/kSkewSlackDiv of the entries.
    static constexpr double kSkewLoadNum = 10.0;
    static constexpr double kSkewLoadDen = 9.0;
    static constexpr uint64_t kSkewSlackDiv = 9;
    /// choose_m() parameters.
    /// Below kMinChooseK it keeps m = k/2.
    /// Otherwise candidates start from ceil(log4 bases), capped at kMaxLog4Bases, and stay within [max(kMinCandidateM, k/2 - kCandidatesBelowHalf), k - kMinMBelowK].
    /// The walk stops after kMaxNonImprovingSteps steps without improvement.
    static constexpr unsigned kMinChooseK = 12;
    static constexpr unsigned kMaxLog4Bases = 31;
    static constexpr unsigned kMinCandidateM = 5;
    static constexpr unsigned kCandidatesBelowHalf = 4;
    static constexpr unsigned kMinMBelowK = 2;
    static constexpr int kMaxNonImprovingSteps = 2;
    /// m must stay below k, and k is at most 32.
    static constexpr unsigned kMaxMinimizerLen = 31;
    /// Packing-pass chunks per thread, and bucket cuts per thread for the parallel passes.
    static constexpr std::size_t kChunksPerThread = 16;
    static constexpr std::size_t kCutsPerThread = 64;
    /// Spare super-k-mer slots reserved per packing chunk.
    static constexpr std::size_t kSuperKmerSlack = 16;
    /// Super-k-mers average about half a window in length, which sizes the reservation per sequence chunk.
    static constexpr unsigned kSuperKmerLenPerWindowDen = 2;
    /// Bisection halves the search range each step.
    static constexpr std::size_t kBisect = 2;
    static constexpr std::size_t kPadWords = 2;   // get_bits() loads two words at a time
    /// OpenMP dynamic-schedule chunks for the sequence scan and the distinct-count pass.
    static constexpr int kScanSchedChunk = 64;
    static constexpr int kDistinctSchedChunk = 1024;
    static constexpr unsigned kInitialRecBits = 3;
    /// Fingerprints are drawn from the hash bits under this mask.
    static constexpr uint64_t kFpHashMask = 0xFFFF;
    /// SplitMix64 finalizer constants for skew_hash(), which index files depend on.
    /// mhash(), by contrast, is plain fmix64 (mix64() in wyhash.hpp).
    static constexpr unsigned kSplitmixShift1 = 31;
    static constexpr uint64_t kSplitmixMul1 = 0xbf58476d1ce4e5b9ULL;
    static constexpr unsigned kSplitmixShift2 = 27;
    static constexpr uint64_t kSplitmixMul2 = 0x94d049bb133111ebULL;
    static constexpr unsigned kSplitmixShift3 = 31;
    uint64_t skew_entries_ = 0;
    uint64_t distinct_ = 0;           // distinct canonical k-mers, as counted by build()
    bool has_distinct_ = false;
    std::vector<std::vector<uint64_t>> skew_own_;
    std::vector<uint64_t> skew_cap_;
    std::vector<const uint64_t*> skew_tables_;
    std::vector<std::size_t> skew_words_;

    /// The bits needed for an index below @p cnt.
    static constexpr unsigned cls_of(uint64_t cnt) noexcept {
        return cnt <= 1 ? 1u : floor_log2(cnt - 1) + 1;
    }
    /// SplitMix64's finalizer applied to @p value ^ kGoldenRatio64.
    static constexpr uint64_t skew_hash(uint64_t value) noexcept {
        value ^= kGoldenRatio64;
        value ^= value >> kSplitmixShift1; value *= kSplitmixMul1;
        value ^= value >> kSplitmixShift2; value *= kSplitmixMul2;
        value ^= value >> kSplitmixShift3;
        return value;
    }
    /// Returns a value in [1, 2^kFpBits), because 0 marks an empty slot.
    static constexpr uint64_t skew_fp(uint64_t hash) noexcept {
        return 1 + (hash & kFpHashMask) % ((uint64_t(1) << kFpBits) - 1);
    }
    static constexpr uint64_t home(uint64_t hash, uint64_t cap) noexcept {   // fast range reduction
        return static_cast<uint64_t>((static_cast<__uint128_t>(hash) * cap) >> kWordBits);
    }
    static void set_bits(uint64_t* words, uint64_t pos, uint64_t val, unsigned nbits) noexcept {
        words[word_index(pos)] |= val << bit_offset(pos);

        if (bit_offset(pos) + nbits > kWordBits) words[word_index(pos) + 1] |= val >> (kWordBits - bit_offset(pos));
    }

    Stats stats_;
    // Owned storage from a build or copying load, and the views every query uses.
    // A borrowing load points the views into the mapped file instead.
    std::vector<uint64_t> str_, rec_, ovf_;
    const uint64_t *str_view_ = nullptr, *rec_view_ = nullptr, *ovf_view_ = nullptr;
    std::size_t str_words_ = 0, rec_words_ = 0, ovf_words_ = 0;

    void adopt_owned() noexcept {
        str_view_ = str_.data(); str_words_ = str_.size();
        rec_view_ = rec_.data(); rec_words_ = rec_.size();
        ovf_view_ = ovf_.data(); ovf_words_ = ovf_.size();
    }
    Mphf mphf_;
};

}

// These headers hold the definitions of the larger members.
#include "spss_choose_m.hpp"
#include "spss_io.hpp"
#include "spss_scan.hpp"
