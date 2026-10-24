#pragma once

/*
 * spss_core.hpp - SpssCore, static exact k-mer membership and dictionary over a spectrum-preserving string set, or SPSS, such as maximal unitigs.
 *
 * The design follows SSHash (Pibiri 2022), which buckets super-k-mer strings by minimizer, hashes the minimizers with a minimal perfect hash, and confirms lookups against the strings.
 * optik builds it from its own parts, phobic_mphf.hpp and Elias-Fano sequence starts, and departs from SSHash where noted below.
 *
 * Layout.
 *   str_      every sequence at 2 bits per base, most significant first inside each 64-bit word, with base i at bits 62 - 2 * (i % 32) of word i / 32.
 *             A k-mer at any offset then costs two word loads and one 128-bit shift.
 *   super-k-mers
 *             maximal runs of consecutive k-mers of one sequence that share
 *             the same minimizer occurrence.  A k-mer's minimizer is its
 *             m-mer whose canonical form has the smallest hash, the leftmost
 *             on ties, so a k-mer and its reverse complement select the same
 *             m-mer, mirrored.
 *   mphf_     partitioned Mphf over the distinct canonical minimizers, giving each its bucket.
 *   rec_      one fixed-width record per bucket, read in one cache-line access, holding
 *               q     the global offset of the minimizer occurrence
 *               near  1 if a k-mer containing that occurrence could extend past its sequence's ends
 *             A bucket holding several super-k-mers, i.e. a repeated
 *             minimizer, sets a flag bit and points into ovf_, a stream of an
 *             Elias-gamma count followed by records in the same layout.
 *
 * Lookup.
 * An exact lookup takes one MPHF access, one record access and one
 *   string access.  Compute x's canonical minimizer and its offset p in x.
 *   For each record of the bucket the only candidate start of x is s = q - p
 *   when x is forward, or s = q - (k - m - p) when it is reverse
 *   complemented.  Comparing the string confirms the candidate and also
 *   rejects keys the MPHF was not built on, and near-boundary records also
 *   check that [s, s + k) stays in one sequence.  When the minimum hash occurs
 *   at several offsets of x, as with the same canonical m-mer twice, every
 *   tied offset is tried.  SSHash scans the whole super-k-mer, whereas storing
 *   the minimizer's offset lets optik jump straight to the one candidate window.
 *
 * Heavy buckets and skew tables.
 * A minimizer shared by many super-k-mers, as
 *   in repeats with up to 164k on human chr1, would make a positive lookup scan
 *   its whole list.  Buckets with more than kSkewMin records get their k-mers
 *   in an open-addressing table, one per size class c holding buckets of up to
 *   2^c records.  Its bit-packed slots of c + kFpBits bits hold a fingerprint
 *   and the index of the k-mer's record in the bucket, probed linearly at a
 *   load of about 0.9 from a fast-range home slot.  A probe chain stays within one
 *   or two cache lines, the fingerprint skips nearly every wrong slot
 *   without touching the string, and the record found is verified against
 *   the string, so membership stays exact.  SSHash's skew index serves the
 *   same role with one MPHF per size class.
 *
 * Sequence bounds.
 * Sequence starts are stored as an Elias-Fano sequence,
 *   so the sequence holding an offset is one rank.
 *
 * Streaming with scan and scan_many.
 * Consecutive query k-mers usually
 *   continue the same string match, so the next k-mer is first compared
 *   against the string at s + 1 going forward or s - 1 going in reverse inside
 *   the current sequence, with no hashing at all.  SSHash calls this an
 *   "extension".
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

    unsigned k() const noexcept { return k_; }
    unsigned m() const noexcept { return m_; }
    const Stats& stats() const noexcept { return stats_; }

    std::size_t bytes() const noexcept {
        std::size_t skew = 0;
        for (std::size_t n : skew_words_) skew += n * 8;
        return (str_words_ + rec_words_ + ovf_words_) * 8 + skew + mphf_.bytes() + ef_starts_.mem_bytes() + 64;
    }

    /// MEstimate - estimated size and lookup cost of an index at minimizer length m (estimate_m()).
    struct MEstimate {
        unsigned m; double bits; double cost;
        double keys, records, overflow, skew, mphf;     // components, keys as a count and the others in bits
        double kmers;                                   // estimated k-mer positions
    };
    /// estimate_m() - estimated size and lookup cost of the index of @seqs at minimizer length @m, without building it.
    ///
    /// All super-k-mers are found (the same sliding-window minimum as build()), but only those whose minimizer falls in a 1/4 hash sample of the minimizer space are kept.
    /// Each occurrence of a sampled minimizer is seen, so its bucket (records, k-mers) is exact, and the totals are scaled by 4.
    /// The bits cover the parts that depend on m, i.e. bucket records, overflow lists, skew tables and the MPHF (about kMphfBitsPerKey per minimizer), while the strings and sequence starts do not change with m.
    /// The overflow stream holds the records of all multi-record buckets, skew buckets included.
    /// cost is the mean number of records scanned per positive lookup in buckets below the skew threshold, and a skew bucket costs one table probe.
    static MEstimate estimate_m(const std::vector<std::string>& seqs, unsigned k, unsigned m) {
        constexpr unsigned kSampleShift = 2;                 // 1/4 of the minimizers
        constexpr uint64_t kPartSeed = 0x7A3C1D, kSampleSeed = 0x5A17E5ULL;
        const unsigned w = k - m + 1;
        const uint64_t mmask = kmer_mask(m);
        uint64_t total_bases = 0;
        for (const auto& s : seqs) total_bases += s.size();
        struct KL { uint64_t key; uint32_t len; };
        // per-thread buffers split into kParts minimizer-hash partitions, so the per-bucket aggregation below runs one partition per thread
        constexpr unsigned kParts = 64;
        auto part_of = [](uint64_t key) { return static_cast<unsigned>(mix64(key, kPartSeed) >> 58); };   // top 6 bits, so kParts = 64
        unsigned nthr = 1;
#ifdef _OPENMP
        nthr = static_cast<unsigned>(omp_get_max_threads());
#endif
        std::vector<std::vector<std::vector<KL>>> per(nthr, std::vector<std::vector<KL>>(kParts));
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
        std::vector<uint64_t> hq, cq;                       // per-thread, reused
        std::vector<std::size_t> dq;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 64)
#endif
        for (int64_t si = 0; si < static_cast<int64_t>(seqs.size()); ++si) {
#ifdef _OPENMP
            auto& out = per[static_cast<std::size_t>(omp_get_thread_num())];
#else
            auto& out = per[0];
#endif
            const std::string& s = seqs[static_cast<std::size_t>(si)];
            std::size_t i = 0;
            while (i < s.size()) {                          // each ACGT run of length >= k
                while (i < s.size() && base_encode(s[i]) == kInvalidBase) ++i;
                std::size_t j = i;
                while (j < s.size() && base_encode(s[j]) != kInvalidBase) ++j;
                if (j - i >= k) {
                    const std::size_t nm = (j - i) - m + 1, nk = (j - i) - k + 1;
                    hq.resize(nm); cq.resize(nm); dq.resize(nm);
                    uint64_t f = 0, r = 0;
                    for (std::size_t t = 0; t < j - i; ++t) {
                        const uint64_t b = base_encode(s[i + t]);
                        f = ((f << 2) | b) & mmask;
                        r = (r >> 2) | ((b ^ 3u) << (2 * (m - 1)));
                        if (t + 1 >= m) {
                            const uint64_t c = f < r ? f : r;
                            cq[t + 1 - m] = c; hq[t + 1 - m] = mhash(c);
                        }
                    }
                    std::size_t qh = 0, qt = 0, next_in = 0, cur_q = ~std::size_t(0);
                    uint32_t cur_len = 0;
                    auto emit = [&](std::size_t q, uint32_t len) {
                        if ((mix64(cq[q], kSampleSeed) >> (64 - kSampleShift)) == 0)
                            out[part_of(cq[q])].push_back({cq[q], len});
                    };
                    for (std::size_t st = 0; st < nk; ++st) {
                        while (next_in <= st + (w - 1)) {
                            while (qt > qh && hq[dq[qt - 1]] > hq[next_in]) --qt;
                            dq[qt++] = next_in++;
                        }
                        while (dq[qh] < st) ++qh;
                        if (dq[qh] != cur_q) {
                            if (cur_len) emit(cur_q, cur_len);
                            cur_q = dq[qh]; cur_len = 0;
                        }
                        ++cur_len;
                    }
                    if (cur_len) emit(cur_q, cur_len);
                }
                i = j;
            }
        }
        }
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
            std::vector<KL> all;
            for (auto& t : per) { all.insert(all.end(), t[p].begin(), t[p].end()); std::vector<KL>().swap(t[p]); }
            std::sort(all.begin(), all.end(), [](const KL& x, const KL& y) { return x.key < y.key; });
            Agg& g = agg[p];
            for (std::size_t a = 0; a < all.size();) {
                std::size_t b = a;
                uint64_t km = 0;
                while (b < all.size() && all[b].key == all[a].key) km += all[b++].len;
                const uint64_t c = b - a;
                g.keys += 1; g.kmers += double(km); g.max_cnt = std::max(g.max_cnt, c);
                // Every multi-record bucket keeps its record list in the overflow stream, and heavy ones add a skew table on top.
                if (c > 1) { g.ovf_recs += double(c); g.ovf_cnt_bits += gamma_len(c - 1); }
                if (c > kSkewMin) {
                    g.skew_bits += double(km) * 10.0 / 9.0 * double(cls_of(c) + kFpBits);   // load 0.9
                    g.scan += double(km);                       // one probe
                } else {
                    g.scan += double(km) * double(c + 1) / 2.0; // mean records tried
                }
                a = b;
            }
        }
        double keys = 0, ovf_recs = 0, ovf_cnt_bits = 0, skew_bits = 0, kmers = 0, scan = 0;
        uint64_t max_cnt = 1;
        for (const Agg& g : agg) {
            keys += g.keys; ovf_recs += g.ovf_recs; ovf_cnt_bits += g.ovf_cnt_bits;
            skew_bits += g.skew_bits; kmers += g.kmers; scan += g.scan;
            max_cnt = std::max(max_cnt, g.max_cnt);
        }
        (void)max_cnt;
        const double ovf_total = scale * (ovf_cnt_bits + ovf_recs * sk_bits);
        const unsigned rec_bits = std::max(sk_bits, bits_for(static_cast<uint64_t>(ovf_total) + 1)) + 1;
        const double D = keys * scale;
        const double bits = D * rec_bits + ovf_total + skew_bits * scale + D * kMphfBitsPerKey;
        return {m, bits, kmers > 0 ? scan / kmers : 0.0,
                D, D * rec_bits, ovf_total, skew_bits * scale, D * kMphfBitsPerKey, kmers * scale};
    }

    /// choose_m() - default minimizer length for @seqs at @k, the m that minimizes the estimated bits per k-mer + kLookupWeight x records scanned per lookup, as estimate_m() predicts them.
    ///
    /// The numbers below were measured when the rule was chosen.
    /// A smaller m shrinks the index but fills buckets, so lookups scan more records.
    /// On E. coli at k = 21, m = 10 gives the smallest index but queries 2.5x slower than m = 12.
    /// A weight of 1 bit per record scanned favours speed, since one record is about 24 ns of a search's DRAM chain.
    /// On chr1 at k = 21, m = 15 queries reads 37% faster than the size optimum m = 13 for 15% more space, and the rule picks m = 17 at k = 31 and m = 13 at k = 15.
    /// The index is exact and lossless at every m, where m only trades size against lookup work, so the sampled estimate can only make the default less than optimal, never wrong.
    /// Candidates are walked down and up from log4(bases) + 1, near the optimum on chr21 and chr1 at k = 15..31, through one non-improving step, within [max(5, k/2 - 4), k - 2], which keeps a window of at least three m-mers.
    /// k < 12 keeps k/2.
    static unsigned choose_m(const std::vector<std::string>& seqs, unsigned k) {
        const unsigned half = std::max(1u, k / 2);
        if (k < 12) return half;
        uint64_t bases = 0;
        for (const auto& s : seqs) bases += s.size();
        unsigned lg = 0;
        while (lg < 31 && (uint64_t(1) << (2 * lg)) < bases) ++lg;   // ceil(log4(bases))
        const unsigned lo = std::max(5u, half > 4 ? half - 4 : 1u), hi = k - 2;
        const unsigned m0 = std::min(std::max(lo, lg + 1), hi);
        constexpr double kLookupWeight = 1.0;
        std::vector<std::pair<unsigned, double>> ev;     // (m, objective)
        auto eval = [&](unsigned m) -> double {
            for (const auto& e : ev) if (e.first == m) return e.second;
            const MEstimate e = estimate_m(seqs, k, m);
            const double obj = (e.kmers > 0 ? e.bits / e.kmers : 0.0) + kLookupWeight * e.cost;
            ev.push_back({m, obj});
            return obj;
        };
        unsigned best_m = m0;
        double best = eval(m0);
        // Walk each way while it improves, through one non-improving step, since the objective is not unimodal at short k (on chr1 at k = 15, m = 12 is worse than both 11 and 13, and 13 is the best).
        for (int dir : {-1, +1}) {
            double run_best = eval(m0);
            int misses = 0;
            for (int m = int(m0) + dir; m >= int(lo) && m <= int(hi); m += dir) {
                const double o = eval(unsigned(m));
                if (o < run_best) { run_best = o; misses = 0; }
                else if (++misses == 2) break;
                if (o < best) { best = o; best_m = unsigned(m); }
            }
        }
        return best_m;
    }

    /// build() - build from ACGT sequences, where non-ACGT characters split a sequence.
    /// Returns false for invalid k/m or a failed MPHF build.
    ///
    /// Exact membership needs no uniqueness, but dense k-mer ids (kmer_id()) assume each canonical k-mer occurs once, as in an SPSS.
    /// @nthreads sets the build workers (0 means one for the chunked passes, while the MPHF resolves 0 to half the CPUs).
    /// @release (optional) is the caller's input vector, freed as soon as it is packed, since only the packing pass reads it.
    /// At short k the input is tens of millions of small strings (2.4 GB for chr1 at k = 15).
    bool build(const std::vector<std::string>& seqs, unsigned k, unsigned m,
               unsigned nthreads = 0, std::vector<std::string>* release = nullptr) {
        if (k == 0 || k > 32 || m == 0 || m >= k || m > 31) return false;
        k_ = k; m_ = m; w_ = k - m + 1;
        kmask_ = kmer_mask(k); mmask_ = kmer_mask(m);
        stats_ = Stats{};
        has_distinct_ = false;
        str_.clear();
        PhaseTimer pt;                                  // OPTIK_TIMING=1

        // 1. pack all ACGT runs of length >= k
        uint64_t nb = 0;
        std::vector<uint64_t> starts;
        auto push_base = [&](uint8_t b) {
            if ((nb & 31) == 0) str_.push_back(0);
            str_.back() |= uint64_t(b) << (62 - 2 * (nb & 31));
            ++nb;
        };
        for (const std::string& s : seqs) {
            std::size_t i = 0, n = s.size();
            while (i < n) {
                while (i < n && base_encode(s[i]) == kInvalidBase) ++i;
                std::size_t j = i;
                while (j < n && base_encode(s[j]) != kInvalidBase) ++j;
                if (j - i >= k) {
                    starts.push_back(nb);
                    for (std::size_t t = i; t < j; ++t) push_base(base_encode(s[t]));
                    ++stats_.seqs;
                    stats_.kmers += (j - i) - k + 1;
                }
                i = j;
            }
        }
        starts.push_back(nb);
        stats_.bases = nb;
        if (release) std::vector<std::string>().swap(*release);   // seqs is unused from here
        str_.push_back(0); str_.push_back(0);          // padding for 2-word loads
        lo_bits_ = bits_for(nb + 1);
        // Sequence starts are coded as Elias-Fano, at about 2 + log2(bases / sequences) bits each instead of lo_bits_.
        // At short k the unitigs are short and numerous, and on chr21 at k = 15 fixed-width starts would take 5.7 of 31 bits per k-mer.
        ef_starts_.build(starts, nb + 1);
        n_starts_ = starts.size();
        adopt_owned();                                  // views for get_mer below
        pt.lap("spss: pack");

        // 2. Super-k-mers per sequence, from a sliding-window minimum over m-mers.
        // q is the absolute offset of the minimizer occurrence, and near says a k-mer anchored at q could reach past its sequence's ends.
        // s0 and len, the super-k-mer's first k-mer and k-mer count, are needed only during the build, since the skew tables enumerate exactly the k-mers of each record.
        // An SK is 16 bytes, the minimizer followed by one word packing q in 40 bits, q - s0 <= k - m in 6 bits, len <= k - m + 1 in 6 bits and near in 1 bit.
        constexpr unsigned kQBits = 40, kDeltaShift = 40, LEN_SHIFT = 46, NEAR_SHIFT = 52;
        constexpr uint64_t kField6Mask = 63u;
        struct SK {
            uint64_t key, w;
            uint64_t q() const noexcept { return w & ((uint64_t(1) << kQBits) - 1); }
            uint64_t s0() const noexcept { return q() - ((w >> kDeltaShift) & kField6Mask); }
            uint32_t len() const noexcept { return static_cast<uint32_t>((w >> LEN_SHIFT) & kField6Mask); }
            uint32_t near() const noexcept { return static_cast<uint32_t>((w >> NEAR_SHIFT) & 1u); }
        };
        if (nb >= (uint64_t(1) << kQBits))
            throw std::length_error("SPSS build: input of 2^40 bases or more");
        std::vector<SK> sks;
        // Any k-mer containing the m-mer at q starts in [q-(k-m), q] and ends by q+k, so the occurrence is flagged if that span can leave [a, b).
        auto mk = [&](uint64_t key, uint64_t q, uint64_t a, uint64_t b, uint64_t s0, uint32_t len) {
            const uint32_t near = (q < a + (k_ - m_)) || (q + k_ > b);
            return SK{key, q | ((q - s0) << kDeltaShift) | (uint64_t(len) << LEN_SHIFT) | (uint64_t(near) << NEAR_SHIFT)};
        };
        // one sequence's super-k-mers, appended to out with reused buffers
        auto seq_sks = [&](std::size_t si, std::vector<SK>& out, std::vector<uint64_t>& hq,
                           std::vector<uint64_t>& cq, std::vector<uint64_t>& dq) {
            const uint64_t a = starts[si], b = starts[si + 1];
            const uint64_t nm = (b - a) - m_ + 1;
            hq.resize(nm); cq.resize(nm); dq.resize(nm);
            for (uint64_t t = 0; t < nm; ++t) {
                const uint64_t y = get_mer(a + t, m_);
                const uint64_t c = std::min(y, rc(y, m_));
                cq[t] = c; hq[t] = mhash(c);
            }
            // monotone deque of m-mer indices with the leftmost minimum at the front
            std::size_t qh = 0, qt = 0;
            uint64_t cur_q = ~uint64_t(0); uint32_t cur_len = 0;
            const uint64_t nk = (b - a) - k_ + 1;
            uint64_t next_in = 0;
            for (uint64_t s = 0; s < nk; ++s) {
                while (next_in <= s + (k_ - m_)) {             // push m-mers of window
                    while (qt > qh && hq[dq[qt - 1]] > hq[next_in]) --qt;
                    dq[qt++] = next_in++;
                }
                while (dq[qh] < s) ++qh;                        // expire
                const uint64_t q = dq[qh];
                if (q != cur_q) {
                    if (cur_len) out.push_back(mk(cq[cur_q], a + cur_q, a, b, a + s - cur_len, cur_len));
                    cur_q = q; cur_len = 0;
                }
                ++cur_len;
            }
            if (cur_len) out.push_back(mk(cq[cur_q], a + cur_q, a, b, a + nk - cur_len, cur_len));
        };
        {
            // Contiguous chunks of sequences, concatenated in chunk order, give the same vector a serial pass builds.
            const std::size_t nseq = starts.size() - 1;
            const unsigned T = nthreads ? nthreads : 1;
            const std::size_t nch = std::min<std::size_t>(nseq, std::size_t(T) * 16);
            std::vector<std::vector<SK>> part(nch);
#ifdef _OPENMP
            #pragma omp parallel num_threads(T)
#endif
            {
                std::vector<uint64_t> hq, cq, dq;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (int64_t c = 0; c < int64_t(nch); ++c) {
                    const std::size_t lo = nseq * std::size_t(c) / nch, hi = nseq * std::size_t(c + 1) / nch;
                    part[c].reserve((starts[hi] - starts[lo]) / std::max(1u, w_ / 2) + 16);
                    for (std::size_t si = lo; si < hi; ++si) seq_sks(si, part[c], hq, cq, dq);
                }
            }
            std::size_t tot = 0;
            for (const auto& v : part) tot += v.size();
            sks.reserve(tot);
            for (auto& v : part) { sks.insert(sks.end(), v.begin(), v.end()); std::vector<SK>().swap(v); }
        }
        stats_.superkmers = sks.size();
        pt.lap("spss: super-k-mers");

        // 3. MPHF over distinct minimizers, bucket the super-k-mers
        // (key, q) is unique (a minimizer occurrence ends one super-k-mer), so the order is total and the parallel sort's result is the serial one In place, without par_sort's merge buffer (1.3 GB on chr1 at k=15), one American-flag pass buckets by the top bits of the minimizer, then the buckets are sorted in parallel.
        radix_sort_inplace(sks, [&](const SK& x) { return x.key; }, 2 * m_,
                           [](const SK& x, const SK& y) {
                               return x.key != y.key ? x.key < y.key : x.q() < y.q(); }, nthreads);
        pt.lap("spss: sort");
        // Chunks of whole buckets for the parallel passes below, which avoid per-bucket offset arrays next to sks (8 bytes per bucket each, and on HPRC 3 such arrays pushed the build over a 24 GiB cap).
        const std::vector<std::size_t> cuts = bucket_cuts(sks, std::size_t(nthreads ? nthreads : 1) * 64);
        const int64_t nch = static_cast<int64_t>(cuts.size()) - 1;
        distinct_ = count_distinct_buckets(sks, cuts);
        has_distinct_ = true;
        pt.lap("spss: distinct");
        std::size_t nkeys = 0;
        for (std::size_t i = 0; i < sks.size(); ++i) nkeys += i == 0 || sks[i].key != sks[i - 1].key;
        stats_.buckets = nkeys;
        {
            std::vector<uint64_t> keys;
            keys.reserve(nkeys);
            for (std::size_t i = 0; i < sks.size(); ++i)
                if (i == 0 || sks[i].key != sks[i - 1].key) keys.push_back(sks[i].key);
            // At alpha 0.94 about 6% of lookups hit the Elias-Fano remap.
            // (against about 15% at the size-optimal 0.85) for about 0.02 bits per k-mer more
            if (!keys.empty() && !mphf_.build(keys.data(), keys.size(), kMphfC, kMphfAlpha, nthreads)) return false;
        }
        pt.lap("spss: mphf");

        // The record layout holds a super-k-mer record or a flagged overflow offset.
        {
            const unsigned sk_bits = lo_bits_ + 1;
            uint64_t ovf_total = 0;
            for (std::size_t i = 0; i < sks.size();) {
                std::size_t j = i;
                while (j < sks.size() && sks[j].key == sks[i].key) ++j;
                if (j - i > 1) ovf_total += gamma_len(j - i - 1) + (j - i) * sk_bits;
                i = j;
            }
            rec_bits_ = std::max(sk_bits, bits_for(ovf_total + 1)) + 1;
        }
        rec_.assign((nkeys * rec_bits_ + 63) / 64 + 2, 0);
        ovf_.clear();
        // Each overflow list's bit offset is a prefix sum over the chunks, then over the buckets inside each chunk as it is written.
        // Neighbouring fields can share a word, so bits are OR-ed in atomically.
        auto ovf_len = [&](std::size_t n) -> uint64_t {
            return n > 1 ? gamma_len(n - 1) + n * (lo_bits_ + 1) : 0;
        };
        std::vector<uint64_t> ovf_at(nch + 1, 0);
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads ? nthreads : 1) schedule(dynamic, 1)
#endif
        for (int64_t c = 0; c < nch; ++c) {
            uint64_t bits = 0;
            for (std::size_t i = cuts[c], j; i < cuts[c + 1]; i = j) {
                for (j = i + 1; j < cuts[c + 1] && sks[j].key == sks[i].key; ++j) {}
                bits += ovf_len(j - i);
            }
            ovf_at[c + 1] = bits;
        }
        for (int64_t c = 0; c < nch; ++c) ovf_at[c + 1] += ovf_at[c];
        const uint64_t ovf_nbits = ovf_at[nch];
        std::vector<uint64_t> ovf_bits_vec((ovf_nbits + 63) / 64 + 2, 0);
        auto put_bits = [](uint64_t* v, uint64_t pos, uint64_t val, unsigned nbits) {
            if (!nbits) return;
            __atomic_fetch_or(&v[pos >> 6], val << (pos & 63), __ATOMIC_RELAXED);
            if ((pos & 63) + nbits > 64)
                __atomic_fetch_or(&v[(pos >> 6) + 1], val >> (64 - (pos & 63)), __ATOMIC_RELAXED);
        };
        uint64_t n_ovf = 0;
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads ? nthreads : 1) schedule(dynamic, 1) reduction(+ : n_ovf)
#endif
        for (int64_t c = 0; c < nch; ++c) {
            uint64_t at = ovf_at[c];
            for (std::size_t i = cuts[c], j; i < cuts[c + 1]; i = j) {
                for (j = i + 1; j < cuts[c + 1] && sks[j].key == sks[i].key; ++j) {}
                const uint64_t bucket = mphf_(sks[i].key);
                if (j - i == 1) {
                    put_bits(rec_.data(), bucket * rec_bits_, pack(sks[i].q(), sks[i].near()) << 1, rec_bits_);
                    continue;
                }
                ++n_ovf;
                put_bits(rec_.data(), bucket * rec_bits_, (at << 1) | 1, rec_bits_);
                put_bits(ovf_bits_vec.data(), at, gamma_word(j - i - 1), gamma_len(j - i - 1));
                at += gamma_len(j - i - 1);
                for (std::size_t t = i; t < j; ++t) {
                    put_bits(ovf_bits_vec.data(), at, pack(sks[t].q(), sks[t].near()), lo_bits_ + 1);
                    at += lo_bits_ + 1;
                }
            }
        }
        stats_.ovf_buckets += n_ovf;
        ovf_ = std::move(ovf_bits_vec);
        ovf_.push_back(0); ovf_.push_back(0);
        adopt_owned();
        pt.lap("spss: records");
        build_skew(sks, nthreads);
        pt.lap("spss: skew");
        return true;
    }

private:
    /// count_distinct_buckets() - distinct canonical k-mers of the sorted super-k-mers, one bucket at a time.
    /// Every occurrence of a canonical k-mer on either strand has the same minimizer, since the m-mer hash is a bijection and the minimum hash names one m-mer value, and so lies in that minimizer's bucket.
    /// A small sort per bucket replaces a full lookup per stored k-mer.
    template <typename SKV>
    uint64_t count_distinct_buckets(const SKV& sks, const std::vector<std::size_t>& cuts) const {
        const int64_t nch = static_cast<int64_t>(cuts.size()) - 1;
        uint64_t d = 0;
#ifdef _OPENMP
        #pragma omp parallel reduction(+ : d)
#endif
        {
            std::vector<uint64_t> buf;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1)
#endif
            for (int64_t c = 0; c < nch; ++c)
                for (std::size_t g = cuts[c], h; g < cuts[c + 1]; g = h) {
                    buf.clear();
                    for (h = g; h < cuts[c + 1] && sks[h].key == sks[g].key; ++h)
                        for (uint32_t u = 0; u < sks[h].len(); ++u) {
                            const uint64_t x = get_mer(sks[h].s0() + u, k_), xr = rc(x, k_);
                            buf.push_back(x < xr ? x : xr);
                        }
                    if (buf.size() > 1) std::sort(buf.begin(), buf.end());
                    for (std::size_t i = 0; i < buf.size(); ++i) d += i == 0 || buf[i] != buf[i - 1];
                }
        }
        return d;
    }

    /// bucket_cuts() - up to @nch + 1 cut points into the key-sorted @sks, each at the start of a bucket, which split it into chunks of whole buckets for the parallel passes.
    template <typename SKV>
    static std::vector<std::size_t> bucket_cuts(const SKV& sks, std::size_t nch) {
        const std::size_t n = sks.size();
        std::vector<std::size_t> cut{0};
        for (std::size_t c = 1; c < nch; ++c) {
            std::size_t i = std::max(cut.back(), n / nch * c + n % nch * c / nch);
            while (i > 0 && i < n && sks[i].key == sks[i - 1].key) ++i;
            if (i > cut.back() && i < n) cut.push_back(i);
        }
        cut.push_back(n);
        return cut;
    }

    /// build_skew() - skew tables over the k-mers of every bucket with more than kSkewMin records.
    /// There is one table per size class, filled sequentially per class with the classes in parallel.
    /// Bucket b of class c inserts each distinct canonical k-mer of its records once, as (fingerprint, t) with t its earliest record.
    /// A lookup returns the first matching entry along the probe chain, which is that same record, so dropping the repeats changes no answer, and it keeps a k-mer repeated n times (low-complexity runs) from building an n-long probe chain with quadratic insert time.
    template <typename SKV>
    void build_skew(const SKV& sks, unsigned nthreads) {
        skew_own_.assign(kMaxClass + 1, {});
        skew_cap_.assign(kMaxClass + 1, 0);
        skew_tables_.assign(kMaxClass + 1, nullptr);
        struct Range { std::size_t i, j; uint64_t n; };     // n = distinct k-mers
        std::vector<std::vector<Range>> by_class(kMaxClass + 1);
        std::vector<uint64_t> entries(kMaxClass + 1, 0);
        // A bucket's k-mers as (canonical, record, order) triples, enumerated by record and then offset, with each canonical k-mer's repeats dropped
        struct E { uint64_t x; uint32_t t, p; };
        auto distinct_of = [&](std::size_t i, std::size_t j, std::vector<E>& v) {
            v.clear();
            for (std::size_t t = i; t < j; ++t)
                for (uint32_t u = 0; u < sks[t].len(); ++u) {
                    const uint64_t x = get_mer(sks[t].s0() + u, k_), xr = rc(x, k_);
                    v.push_back({x < xr ? x : xr, static_cast<uint32_t>(t - i), static_cast<uint32_t>(v.size())});
                }
            std::sort(v.begin(), v.end(), [](const E& a, const E& b) { return a.x != b.x ? a.x < b.x : a.p < b.p; });
            v.erase(std::unique(v.begin(), v.end(), [](const E& a, const E& b) { return a.x == b.x; }), v.end());
            std::sort(v.begin(), v.end(), [](const E& a, const E& b) { return a.p < b.p; });
        };
        std::vector<Range> heavy;
        for (std::size_t i = 0; i < sks.size();) {
            std::size_t j = i;
            while (j < sks.size() && sks[j].key == sks[i].key) ++j;
            if (j - i > kSkewMin) heavy.push_back({i, j, 0});
            i = j;
        }
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            std::vector<E> v;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 1)
#endif
            for (std::size_t h = 0; h < heavy.size(); ++h) {
                distinct_of(heavy[h].i, heavy[h].j, v);
                heavy[h].n = v.size();
            }
        }
        for (const Range& rg : heavy) {
            const unsigned c = cls_of(rg.j - rg.i);
            by_class[c].push_back(rg);
            entries[c] += rg.n;
        }
        std::vector<unsigned> todo;
        skew_entries_ = 0;
        for (unsigned c = 0; c <= kMaxClass; ++c) {
            skew_entries_ += entries[c];
            if (!entries[c]) continue;
            // A load of about 0.9 with fast-range home slots, so no power-of-two rounding
            skew_cap_[c] = entries[c] + entries[c] / 9 + 1;
            skew_own_[c].assign((skew_cap_[c] * (c + kFpBits) + 63) / 64 + 2, 0);
            todo.push_back(c);
        }
        (void)nthreads;
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1)
#endif
        for (std::size_t ti = 0; ti < todo.size(); ++ti) {
            const unsigned c = todo[ti];
            const unsigned w = c + kFpBits;
            std::vector<uint64_t>& W = skew_own_[c];
            const uint64_t cap = skew_cap_[c];
            std::vector<E> v;
            for (const Range& rg : by_class[c]) {
                distinct_of(rg.i, rg.j, v);
                for (const E& e : v) {
                    const uint64_t h = skew_hash(e.x);
                    const uint64_t val = (skew_fp(h) << c) | e.t;
                    for (uint64_t sl = home(h, cap);; sl = (sl + 1 == cap) ? 0 : sl + 1) {
                        if (get_bits(W.data(), sl * w, w) == 0) { set_bits(W.data(), sl * w, val, w); break; }
                    }
                }
            }
        }
        skew_words_.assign(kMaxClass + 1, 0);
        for (unsigned c = 0; c <= kMaxClass; ++c) {
            skew_tables_[c] = skew_own_[c].data();
            skew_words_[c] = skew_own_[c].size();
        }
    }

public:

    /// lookup() - exact membership of k-mer @x, given in forward packing, either strand.
    /// Returns the string offset of the matching occurrence, or kNotFound.
    uint64_t lookup(uint64_t x) const noexcept {
        if (stats_.buckets == 0) return kNotFound;
        MinInfo mi = minimizer(x);
        return probe(x, mi, nullptr);
    }

    bool contains(uint64_t x) const noexcept { return lookup(x) != kNotFound; }

    // Dictionary operations, SSHash's Lookup and Access.
    // String i occupies [start_i, start_i+1) and holds len_i - k + 1 k-mers, so the k-mer at offset s of string i has the dense identifier s - (k - 1) * i.
    // Identifiers cover [0, n) and are a bijection onto the indexed k-mers when the input holds each k-mer once (unitigs / an SPSS, as SSHash also requires).

    /// kmer_id() - dense identifier in [0, num_kmers) of k-mer @x (either strand), or kNotFound.
    /// *@rev: x is the reverse complement of the stored k-mer.
    uint64_t kmer_id(uint64_t x, bool* rev = nullptr) const noexcept {
        if (stats_.buckets == 0) return kNotFound;
        const uint64_t s = probe(x, minimizer(x), rev);
        if (s == kNotFound) return kNotFound;
        return s - static_cast<uint64_t>(k_ - 1) * seq_index(s);
    }

    /// access() - the stored k-mer (forward packing, as it appears in the string) with identifier @id < num_kmers().
    uint64_t access(uint64_t id) const noexcept {
        // largest i with start_i - (k-1) i <= id
        std::size_t a = 0, b = n_starts_ - 1;
        while (b - a > 1) {
            const std::size_t mid = (a + b) / 2;
            const uint64_t before = start_at(mid) - static_cast<uint64_t>(k_ - 1) * mid;
            if (before <= id) a = mid; else b = mid;
        }
        return get_mer(id + static_cast<uint64_t>(k_ - 1) * a, k_);
    }

    uint64_t num_kmers() const noexcept { return stats_.kmers; }
    /// k-mers indexed through skew tables by the last build(), or 0 after load.
    uint64_t skew_entries() const noexcept { return skew_entries_; }

    /// count_distinct() - distinct canonical k-mers, where a stored k-mer counts once, at the first occurrence a lookup returns.
    /// The count is known from build().
    /// After a load it walks the sequences directly, so no id has to be inverted with a binary search over the starts.
    uint64_t count_distinct() const noexcept {
        if (has_distinct_) return distinct_;
        if (stats_.buckets == 0) return 0;
        const int64_t ns = static_cast<int64_t>(n_starts_) - 1;
        uint64_t d = 0;
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1024) reduction(+ : d)
#endif
        for (int64_t si = 0; si < ns; ++si) {
            const uint64_t a = start_at(si);
            const uint64_t b = start_at(si + 1);
            for (uint64_t s = a; s + k_ <= b; ++s)
                d += probe(get_mer(s, k_), minimizer(get_mer(s, k_)), nullptr) == s;
        }
        return d;
    }

    /// for_each_kmer() - @f(canonical k-mer) for every stored k-mer, in sequence order.
    template <typename F> void for_each_kmer(F&& f) const {
        for (std::size_t si = 0; si + 1 < n_starts_; ++si) {
            const uint64_t a = start_at(si);
            const uint64_t b = start_at(si + 1);
            for (uint64_t s = a; s + k_ <= b; ++s) {
                const uint64_t x = get_mer(s, k_), r = rc(x, k_);
                f(x < r ? x : r);
            }
        }
    }

    /// scan() and scan_many() - streaming lookup of every k-mer of the query records, counted as SSHash's query counts them.
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

    void scan(const char* seq, std::size_t len, ScanStats& out) const noexcept {
        scan_many(&seq, &len, 1, out);
    }

    void scan_many(const char* const* seqs, const std::size_t* lens, std::size_t nseq,
                   ScanStats& out) const noexcept {
        GroupBuf g;
        // Longer records wait here for a lane (local, since scan_many runs concurrently on one shared index).
        std::vector<std::pair<const char*, std::size_t>> lq;
        lq.reserve(nseq);
        for (std::size_t r = 0; r < nseq; ++r) {
            if (lens[r] >= k_ && lens[r] - k_ + 1 <= kGroupMaxKmers) {
                // collect every valid k-mer of the tiny record
                uint64_t f = 0; std::size_t v = 0;
                for (std::size_t i = 0; i < lens[r]; ++i) {
                    const uint8_t b = base_encode(seqs[r][i]);
                    if (b == kInvalidBase) { v = 0; f = 0; continue; }
                    f = ((f << 2) | b) & kmask_;
                    if (++v >= k_) {
                        g.x[g.n++] = f;
                        if (g.n == kGroup) flush_group(g, out);
                    }
                }
            } else {
                lq.push_back({seqs[r], lens[r]});
            }
        }
        flush_group(g, out);
        if (lq.empty()) return;
        if (stats_.buckets == 0) {                    // empty index, which counts but never hits
            for (const auto& x : lq) scan_one(x.first, x.second, out);
        } else {
            scan_lanes(lq.data(), lq.size(), out);
        }
    }

    /// Reference single-cursor scan, one record and one dependent miss chain at a time.
    /// scan_many runs records through scan_lanes instead, and the stress tests check that both give identical counts.
    void scan_reference(const char* seq, std::size_t len, ScanStats& out) const noexcept {
        scan_one(seq, len, out);
    }

private:
    struct GroupBuf { uint64_t x[kGroup]; unsigned n = 0; };

    /// flush_group() - look up the buffered independent k-mers with staged prefetching.
    void flush_group(GroupBuf& g, ScanStats& out) const noexcept {
        if (g.n == 0) return;
        const unsigned n = g.n;
        g.n = 0;
        out.valid += n;
        out.searches += n;
        if (stats_.buckets == 0) return;
        MinInfo mi[kGroup];
        uint64_t bi[kGroup], en[kGroup];
        for (unsigned i = 0; i < n; ++i) {
            mi[i] = minimizer(g.x[i]);
            mphf_.prefetch(mi[i].key);
        }
        for (unsigned i = 0; i < n; ++i) {
            bi[i] = mphf_(mi[i].key);
            __builtin_prefetch(reinterpret_cast<const uint8_t*>(rec_view_) + ((bi[i] * rec_bits_) >> 3), 0, 3);
        }
        for (unsigned i = 0; i < n; ++i) {
            en[i] = get_bits(rec_view_, bi[i] * rec_bits_, rec_bits_);
            if (!(en[i] & 1)) {
                const uint64_t q = en[i] >> 2;
                __builtin_prefetch(str_view_ + ((q >= k_ ? q - k_ : 0) >> 5), 0, 3);
                __builtin_prefetch(str_view_ + (q >> 5) + 1, 0, 3);
            } else {
                __builtin_prefetch(reinterpret_cast<const uint8_t*>(ovf_view_) + ((en[i] >> 1) >> 3), 0, 3);
            }
        }
        for (unsigned i = 0; i < n; ++i) {
            const bool hit = entry_has_key(en[i], mi[i].key)
                             && probe_entry(g.x[i], mi[i], en[i], nullptr) != kNotFound;
            out.positive += hit;
        }
    }

    // Interleaved streaming scan.
    // A full search is a chain of dependent misses, from MPHF pilot to bucket record to overflow list to string, so one cursor per thread waits out one DRAM miss at a time.
    // scan_lanes advances kLanes records round-robin, and each search is split into stages that prefetch the line the next stage reads, so up to kLanes misses are in flight.
    // Per record the work and counts match scan_one's, with the same extension rule, bucket cache and probe.
    // Only the order across records changes, and the counts are sums.
    static constexpr unsigned kLanes = 16;
    static constexpr uint64_t kProbeChunk = 32;         // overflow records decoded per chunk
    /// Largest number of m-mers in a k-mer (k - m + 1 with k <= 32, m >= 1), also the size of the m-mer rings below, so it must be a power of two.
    static constexpr unsigned kMaxWindowMmers = 32;
    static constexpr uint64_t kRingMask = kMaxWindowMmers - 1;
    /// MinInfo - a k-mer's canonical minimizer and the offsets in the k-mer where its minimum hash occurs, npos in number.
    struct MinInfo { uint64_t key; uint32_t npos; uint8_t pos[kMaxWindowMmers]; };
    struct Lane {
        const char* seq = nullptr; std::size_t len = 0, i = 0;
        uint64_t fwd = 0, mf = 0, mr = 0; std::size_t valid = 0;
        uint64_t midx = 0, hashed = 0, best = 0, bl = 0, br = 0; bool have_min = false;
        uint64_t cur = kNotFound; bool cur_rev = false; uint64_t seq_lo = 0, seq_hi = 0;
        uint64_t cached_key = ~uint64_t(0), cached_entry = 0; bool cached_ok = false;
        // Pending search stage, where 0 is none, 1 pilot prefetched, 2 record prefetched, 3 overflow head prefetched and 4 string prefetched
        unsigned stage = 0; bool fresh = false;
        MinInfo mi; uint64_t x = 0, entry = 0;
        uint64_t rcn[kMaxWindowMmers], rh[kMaxWindowMmers];   // rings of canonical m-mers and hashes
    };
    void prefetch_window(uint64_t q) const noexcept {
        __builtin_prefetch(str_view_ + ((q >= k_ ? q - k_ : 0) >> 5), 0, 3);
        __builtin_prefetch(str_view_ + ((q + k_) >> 5), 0, 3);
    }
    void lane_reset(Lane& L, const char* s, std::size_t n) const noexcept {
        L.seq = s; L.len = n; L.i = 0;
        L.fwd = L.mf = L.mr = 0; L.valid = 0; L.midx = L.hashed = 0; L.have_min = false;
        L.cur = kNotFound; L.cached_key = ~uint64_t(0); L.cached_ok = false; L.stage = 0;
    }
    /// lane_run() - run lane @L's record forward until it needs a full search whose data is not cached (returning true, with a stage pending and its first line prefetched) or the record ends (false).
    /// Mirrors scan_one's loop.
    bool lane_run(Lane& L, ScanStats& out) const noexcept {
        const unsigned W = k_ - m_ + 1;
        const unsigned mshift = 2 * (m_ - 1);
        while (L.i < L.len) {
            const uint8_t b = base_encode(L.seq[L.i++]);
            if (b == kInvalidBase) { L.valid = 0; L.fwd = 0; L.midx = L.hashed = 0; L.have_min = false; L.cur = kNotFound; continue; }
            L.fwd = ((L.fwd << 2) | b) & kmask_;
            L.mf = ((L.mf << 2) | b) & mmask_;
            L.mr = (L.mr >> 2) | (uint64_t(b ^ 3u) << mshift);
            ++L.valid;
            if (L.valid >= m_) { L.rcn[L.midx & kRingMask] = L.mf < L.mr ? L.mf : L.mr; ++L.midx; }
            if (L.valid < k_) continue;
            ++out.valid;
            if (L.cur != kNotFound) {
                const bool room = L.cur_rev ? (L.cur > L.seq_lo) : (L.cur + 1 + k_ <= L.seq_hi);
                const uint64_t ns = L.cur_rev ? L.cur - 1 : L.cur + 1;
                if (room && get_mer(ns, k_) == (L.cur_rev ? rc(L.fwd, k_) : L.fwd)) {
                    L.cur = ns; ++out.positive; ++out.extensions;
                    continue;
                }
            }
            const uint64_t first = L.midx - W;
            if (L.hashed < first) { L.hashed = first; L.have_min = false; }
            if (L.have_min && L.bl < first) L.have_min = false;
            if (!L.have_min) {
                for (; L.hashed < L.midx; ++L.hashed) L.rh[L.hashed & kRingMask] = mhash(L.rcn[L.hashed & kRingMask]);
                uint64_t best = ~uint64_t(0), bl = first, br = first;
                for (uint64_t t = first; t < L.midx; ++t) {
                    const uint64_t h = L.rh[t & kRingMask];
                    const bool lt = h < best, le = h <= best;
                    best = lt ? h : best;
                    bl = lt ? t : bl;
                    br = le ? t : br;
                }
                L.best = best; L.bl = bl; L.br = br;
                L.have_min = true;
            } else {
                for (; L.hashed < L.midx; ++L.hashed) {
                    const uint64_t h = mhash(L.rcn[L.hashed & kRingMask]);
                    L.rh[L.hashed & kRingMask] = h;
                    if (h < L.best) { L.best = h; L.bl = L.br = L.hashed; }
                    else if (h == L.best) L.br = L.hashed;
                }
            }
            L.mi.key = L.rcn[L.bl & kRingMask];
            L.mi.pos[0] = static_cast<uint8_t>(L.bl - first);
            L.mi.pos[1] = static_cast<uint8_t>(L.br - first);
            L.mi.npos = (L.br != L.bl) ? 2 : 1;
            L.x = L.fwd;
            ++out.searches;
            if (L.mi.key != L.cached_key) {                // a new bucket, so stage it
                L.fresh = true;
                mphf_.prefetch(L.mi.key);
                L.stage = 1;
                return true;
            }
            if (!L.cached_ok) { L.cur = kNotFound; continue; }   // minimizer absent
            L.fresh = false;                                // bucket known, so only the string remains
            L.entry = L.cached_entry;
            if (!(L.entry & 1)) { prefetch_window(L.entry >> 2); L.stage = 4; }
            else {
                __builtin_prefetch(reinterpret_cast<const uint8_t*>(ovf_view_) + ((L.entry >> 1) >> 3), 0, 3);
                L.stage = 3;
            }
            return true;
        }
        return false;
    }
    /// lane_stage() - one stage of lane @L's pending search, true while still pending.
    bool lane_stage(Lane& L, ScanStats& out) const noexcept {
        switch (L.stage) {
        case 1: {
            const uint64_t bkt = mphf_(L.mi.key);
            L.entry = bkt;                                  // bucket id until stage 2
            __builtin_prefetch(reinterpret_cast<const uint8_t*>(rec_view_) + ((bkt * rec_bits_) >> 3), 0, 3);
            L.stage = 2;
            return true;
        }
        case 2: {
            L.entry = get_bits(rec_view_, L.entry * rec_bits_, rec_bits_);
            if (!(L.entry & 1)) { prefetch_window(L.entry >> 2); L.stage = 4; }
            else {
                __builtin_prefetch(reinterpret_cast<const uint8_t*>(ovf_view_) + ((L.entry >> 1) >> 3), 0, 3);
                L.stage = 3;
            }
            return true;
        }
        case 3: {                                           // first overflow record
            // Prefetching more records' windows measured no faster, since a probe usually stops early and extra prefetches crowd the fill buffers.
            unsigned cw;
            ovf_count(L.entry >> 1, cw);
            prefetch_window(get_bits(ovf_view_, (L.entry >> 1) + cw, lo_bits_ + 1) >> 1);
            L.stage = 4;
            return true;
        }
        default: {                                          // stage 4 compares and probes
            L.stage = 0;
            if (L.fresh) {
                L.cached_key = L.mi.key;
                L.cached_entry = L.entry;
                L.cached_ok = entry_has_key(L.entry, L.mi.key);
                if (!L.cached_ok) { L.cur = kNotFound; return false; }
            }
            bool rev = false;
            const uint64_t s = probe_entry(L.x, L.mi, L.entry, &rev);
            if (s == kNotFound) { L.cur = kNotFound; return false; }
            ++out.positive;
            L.cur = s; L.cur_rev = rev;
            seq_bounds(s, L.seq_lo, L.seq_hi);
            return false;
        }
        }
    }
    void scan_lanes(const std::pair<const char*, std::size_t>* recs, std::size_t n,
                    ScanStats& out) const noexcept {
        Lane lanes[kLanes];
        std::size_t next = 0;
        unsigned live = 0;
        bool on[kLanes] = {};
        // Start a lane on the next record that needs a search, or return false when no records are left.
        auto feed = [&](Lane& L) -> bool {
            while (next < n) {
                lane_reset(L, recs[next].first, recs[next].second);
                ++next;
                if (lane_run(L, out)) return true;
            }
            return false;
        };
        for (unsigned l = 0; l < kLanes; ++l) { on[l] = feed(lanes[l]); live += on[l]; }
        while (live) {
            for (unsigned l = 0; l < kLanes; ++l) {
                if (!on[l]) continue;
                Lane& L = lanes[l];
                if (lane_stage(L, out)) continue;           // still pending
                if (lane_run(L, out)) continue;             // next search staged
                if (!feed(L)) { on[l] = false; --live; }
            }
        }
    }

    /// scan_one() - single-cursor streaming scan of one record (see scan_many()).
    void scan_one(const char* seq, std::size_t len, ScanStats& out) const noexcept {
        if (stats_.buckets == 0) {                    // empty index, which counts but never hits
            for (std::size_t i = 0, v = 0; i < len; ++i) {
                if (base_encode(seq[i]) == kInvalidBase) { v = 0; continue; }
                if (++v >= k_) { ++out.valid; ++out.searches; }
            }
            return;
        }
        const unsigned W = k_ - m_ + 1;           // m-mers per k-mer (<= 32)
        const unsigned mshift = 2 * (m_ - 1);
        uint64_t rcn[kMaxWindowMmers], rh[kMaxWindowMmers];   // rings of canonical m-mers and hashes
        uint64_t fwd = 0, mf = 0, mr = 0;
        std::size_t valid = 0;
        uint64_t midx = 0, hashed = 0;             // m-mers seen and hashed in the segment
        // Incremental window minimum over hashed m-mers at absolute indices, with the leftmost and rightmost positions of the minimum hash
        uint64_t best = ~uint64_t(0), bl = 0, br = 0;
        bool have_min = false;
        uint64_t cur = kNotFound; bool cur_rev = false; uint64_t seq_lo = 0, seq_hi = 0;
        uint64_t cached_key = ~uint64_t(0), cached_entry = 0;
        bool cached_ok = false;
        for (std::size_t i = 0; i < len; ++i) {
            const uint8_t b = base_encode(seq[i]);
            if (b == kInvalidBase) { valid = 0; fwd = 0; midx = hashed = 0; have_min = false; cur = kNotFound; continue; }
            fwd = ((fwd << 2) | b) & kmask_;
            mf = ((mf << 2) | b) & mmask_;
            mr = (mr >> 2) | (uint64_t(b ^ 3u) << mshift);
            ++valid;
            if (valid >= m_) { rcn[midx & kRingMask] = mf < mr ? mf : mr; ++midx; }
            if (valid < k_) continue;
            ++out.valid;
            if (cur != kNotFound) {
                const bool room = cur_rev ? (cur > seq_lo) : (cur + 1 + k_ <= seq_hi);
                const uint64_t ns = cur_rev ? cur - 1 : cur + 1;
                if (room && get_mer(ns, k_) == (cur_rev ? rc(fwd, k_) : fwd)) {
                    cur = ns; ++out.positive; ++out.extensions;
                    continue;
                }
            }
            // A full search hashes the window's not-yet-hashed m-mers and updates the window minimum incrementally, rescanning (branch-free) only when the minimum left the window or hashing restarted after a run.
            const uint64_t first = midx - W;
            if (hashed < first) { hashed = first; have_min = false; }
            if (have_min && bl < first) have_min = false;
            if (!have_min) {
                for (; hashed < midx; ++hashed) rh[hashed & kRingMask] = mhash(rcn[hashed & kRingMask]);
                best = ~uint64_t(0); bl = br = first;
                for (uint64_t t = first; t < midx; ++t) {
                    const uint64_t h = rh[t & kRingMask];
                    const bool lt = h < best, le = h <= best;
                    best = lt ? h : best;
                    bl = lt ? t : bl;
                    br = le ? t : br;
                }
                have_min = true;
            } else {
                for (; hashed < midx; ++hashed) {
                    const uint64_t h = mhash(rcn[hashed & kRingMask]);
                    rh[hashed & kRingMask] = h;
                    if (h < best) { best = h; bl = br = hashed; }
                    else if (h == best) br = hashed;
                }
            }
            MinInfo mi;
            mi.key = rcn[bl & kRingMask];
            mi.pos[0] = static_cast<uint8_t>(bl - first);
            mi.pos[1] = static_cast<uint8_t>(br - first);
            mi.npos = (br != bl) ? 2 : 1;
            ++out.searches;
            if (mi.key != cached_key) {
                cached_key = mi.key;
                cached_entry = bucket_entry(mi.key);
                cached_ok = entry_has_key(cached_entry, mi.key);
            }
            if (!cached_ok) {                                  // minimizer absent
                cur = kNotFound; continue;
            }
            bool rev = false;
            const uint64_t s = probe_entry(fwd, mi, cached_entry, &rev);
            if (s == kNotFound) { cur = kNotFound; continue; }
            ++out.positive;
            cur = s; cur_rev = rev;
            seq_bounds(s, seq_lo, seq_hi);
        }
    }

public:
    /// save() - serialize the index, as load() reads it back.
    void save(std::vector<uint8_t>& out) const {
        auto put = [&](uint64_t v) { uint8_t b[8]; std::memcpy(b, &v, 8); out.insert(out.end(), b, b + 8); };
        auto putv = [&](const uint64_t* v, std::size_t n) {
            put(n);
            const std::size_t o = out.size();
            out.resize(o + n * 8);
            if (n) std::memcpy(out.data() + o, v, n * 8);
        };
        put(k_); put(m_); put(lo_bits_); put(rec_bits_); put(n_starts_);
        put(stats_.kmers); put(stats_.bases); put(stats_.seqs); put(stats_.superkmers);
        put(stats_.buckets); put(stats_.ovf_buckets);
        putv(str_view_, str_words_); putv(rec_view_, rec_words_); putv(ovf_view_, ovf_words_);
        // The skew tables' layout parameters, then the tables, 8-byte aligned like every array above.
        put(kSkewMin); put(kFpBits); put(kMaxClass);
        for (unsigned c = 0; c <= kMaxClass; ++c) {
            put(skew_cap_[c]);
            putv(skew_tables_[c], skew_words_[c]);
        }
        mphf_.save(out);
        ef_starts_.save(out);
    }
    /// load() - read a save() blob, returning the bytes consumed, or 0 for a malformed one.
    /// With borrow = true every array is used in place, the MPHF's and the Elias-Fano starts' included (zero-copy, so `p` must stay mapped and 8-byte aligned for the lifetime of this object), and otherwise they are copied.
    std::size_t load(const uint8_t* p, std::size_t sz, bool borrow = false) {
        has_distinct_ = false;
        std::size_t o = 0;
        auto get = [&](uint64_t& v) { if (o + 8 > sz) return false; std::memcpy(&v, p + o, 8); o += 8; return true; };
        auto getv = [&](std::vector<uint64_t>& v, const uint64_t*& view, std::size_t& vn) {
            uint64_t n; if (!get(n) || n > (sz - o) / 8) return false;
            if (borrow) {
                v.clear(); v.shrink_to_fit();
                view = reinterpret_cast<const uint64_t*>(p + o);
            } else {
                v.resize(n); std::memcpy(v.data(), p + o, n * 8); view = v.data();
            }
            vn = n; o += n * 8; return true; };
        uint64_t k, m, lb, rb, ns;
        if (!get(k) || !get(m) || !get(lb) || !get(rb) || !get(ns)) return 0;
        if (!get(stats_.kmers) || !get(stats_.bases) || !get(stats_.seqs) || !get(stats_.superkmers)
            || !get(stats_.buckets) || !get(stats_.ovf_buckets)) return 0;
        if (!getv(str_, str_view_, str_words_) || !getv(rec_, rec_view_, rec_words_)
            || !getv(ovf_, ovf_view_, ovf_words_)) return 0;
        skew_own_.assign(kMaxClass + 1, {});
        skew_cap_.assign(kMaxClass + 1, 0);
        skew_tables_.assign(kMaxClass + 1, nullptr);
        skew_words_.assign(kMaxClass + 1, 0);
        uint64_t smin, fpb, maxc;
        if (!get(smin) || !get(fpb) || !get(maxc)) return 0;
        if (smin != kSkewMin || fpb != kFpBits || maxc != kMaxClass) return 0;
        for (unsigned c = 0; c <= kMaxClass; ++c) {
            if (!get(skew_cap_[c])) return 0;
            if (!getv(skew_own_[c], skew_tables_[c], skew_words_[c])) return 0;
            // a table must hold its capacity's slots, or the file is corrupt
            if (skew_cap_[c] && skew_words_[c] * 64 < skew_cap_[c] * (c + kFpBits)) return 0;
        }
        const std::size_t used = mphf_.load(p + o, sz - o, borrow);
        if (!used) return 0;
        o += used;
        const std::size_t eu = ef_starts_.load(p + o, sz - o, borrow);
        if (!eu || ef_starts_.size() != ns || !ef_starts_.sane()) return 0;
        o += eu;
        k_ = unsigned(k); m_ = unsigned(m); w_ = k_ - m_ + 1;
        lo_bits_ = unsigned(lb); rec_bits_ = unsigned(rb);
        n_starts_ = ns;
        kmask_ = kmer_mask(k_); mmask_ = kmer_mask(m_);
        return o;
    }

    /// rc() - reverse complement of a 2-bit packed @w-mer (w <= 32), branch-free.
    static constexpr uint64_t rc(uint64_t x, unsigned w) noexcept {
        x = ~x;
        x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
        x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
        x = __builtin_bswap64(x);
        return x >> (64 - 2 * w);
    }

private:

    static constexpr unsigned bits_for(uint64_t v) noexcept {   // bits to represent [0, v)
        unsigned b = 0;
        while ((uint64_t(1) << b) < v) ++b;
        return b ? b : 1;
    }
    /// mhash() - MurmurHash3 fmix64, the hash that orders canonical m-mers.
    /// Index files depend on it.
    static constexpr uint64_t mhash(uint64_t c) noexcept {
        c ^= c >> 33; c *= 0xff51afd7ed558ccdULL;
        c ^= c >> 33; c *= 0xc4ceb9fe1a85ec53ULL;
        c ^= c >> 33;
        return c;
    }
    /// pack() - a record's payload, q << 1 | near.
    static constexpr uint64_t pack(uint64_t q, uint32_t near) noexcept {
        return (q << 1) | near;
    }
    /// gamma_len() and gamma_word() - Elias-gamma code (Elias 1975) of n >= 1 in the LSB-first bit streams, i.e. L = floor(log2 n) zero bits, a one, then the L low bits of n (2L + 1 bits).
    /// It codes the overflow-bucket record counts, mostly 2-4 records in 1-3 bits, where a fixed field took bits(max count) = 15-17 bits.
    static constexpr unsigned gamma_len(uint64_t n) noexcept {
        return 2u * (63u - static_cast<unsigned>(__builtin_clzll(n))) + 1u;
    }
    static constexpr uint64_t gamma_word(uint64_t n) noexcept {
        const unsigned L = 63u - static_cast<unsigned>(__builtin_clzll(n));
        return (uint64_t(1) << L) | ((n & ((uint64_t(1) << L) - 1)) << (L + 1));
    }
    /// ovf_count() - record count of the overflow entry at bit offset @off, and the width of its count field.
    uint64_t ovf_count(uint64_t off, unsigned& width) const noexcept {
        const uint64_t w = get_bits(ovf_view_, off, 64);
        const unsigned L = static_cast<unsigned>(__builtin_ctzll(w));
        width = 2u * L + 1u;
        return ((uint64_t(1) << L) | ((w >> (L + 1)) & ((uint64_t(1) << L) - 1))) + 1;   // n = cnt - 1
    }
    static uint64_t get_bits(const uint64_t* v, uint64_t pos, unsigned nbits) noexcept {
        const uint64_t w0 = v[pos >> 6], w1 = v[(pos >> 6) + 1];
        const unsigned sh = pos & 63;
        uint64_t x = sh ? ((w0 >> sh) | (w1 << (64 - sh))) : w0;
        return nbits == 64 ? x : (x & ((uint64_t(1) << nbits) - 1));
    }
    /// get_mer() - the @w-mer (w <= 32) starting at global base offset @s.
    uint64_t get_mer(uint64_t s, unsigned w) const noexcept {
        const uint64_t i = s >> 5;
        const __uint128_t v = ((__uint128_t)str_view_[i] << 64) | str_view_[i + 1];
        return static_cast<uint64_t>((v << (2 * (s & 31))) >> (128 - 2 * w));
    }

    /// minimizer() - the minimizer of k-mer @x and its tied offsets.
    MinInfo minimizer(uint64_t x) const noexcept {
        const uint64_t xr = rc(x, k_);
        const unsigned W = k_ - m_ + 1;
        uint64_t h[kMaxWindowMmers] = {}, c[kMaxWindowMmers] = {};
        for (unsigned p = 0; p < W; ++p) {
            const uint64_t y = (x >> (2 * (k_ - m_ - p))) & mmask_;
            const uint64_t yr = (xr >> (2 * p)) & mmask_;     // rc of the m-mer at p
            c[p] = y < yr ? y : yr;
            h[p] = mhash(c[p]);
        }
        uint64_t best = ~uint64_t(0); unsigned bp = 0;
        for (unsigned p = 0; p < W; ++p) {
            const bool lt = h[p] < best;
            best = lt ? h[p] : best;
            bp = lt ? p : bp;
        }
        MinInfo mi;
        mi.key = c[bp];
        mi.npos = 0;
        for (unsigned p = bp; p < W; ++p) {
            mi.pos[mi.npos] = static_cast<uint8_t>(p);
            mi.npos += h[p] == best;
        }
        return mi;
    }

    /// probe() - try every candidate window of @x's bucket and return the string offset, or kNotFound.
    uint64_t probe(uint64_t x, const MinInfo& mi, bool* rev_out) const noexcept {
        return probe_entry(x, mi, bucket_entry(mi.key), rev_out);
    }
    /// entry_has_key() - true iff bucket record @entry belongs to minimizer
    /// @key.
    /// Every occurrence stored in a bucket spells its canonical minimizer, so one m-mer compare rejects keys the MPHF maps to a foreign bucket.
    /// The string doubles as the key store and needs no fingerprint bits.
    bool entry_has_key(uint64_t entry, uint64_t key) const noexcept {
        uint64_t q;
        if (!(entry & 1)) q = entry >> 2;
        else {
            unsigned cw;
            ovf_count(entry >> 1, cw);
            q = get_bits(ovf_view_, (entry >> 1) + cw, lo_bits_ + 1) >> 1;
        }
        const uint64_t y = get_mer(q, m_);
        return (y < rc(y, m_) ? y : rc(y, m_)) == key;
    }
    /// bucket_entry() - the record of minimizer @key's bucket.
    uint64_t bucket_entry(uint64_t key) const noexcept {
        return get_bits(rec_view_, mphf_(key) * rec_bits_, rec_bits_);
    }
    /// probe_entry() - probe() with the bucket record @e already read.
    uint64_t probe_entry(uint64_t x, const MinInfo& mi, uint64_t e,
                         bool* rev_out) const noexcept {
        const uint64_t xr = rc(x, k_);
        auto try_rec = [&](uint64_t r) -> uint64_t {
            const uint64_t q = r >> 1;
            const bool near = r & 1;
            // An exact string match at s is a hit unless the window crosses into a neighbouring sequence, which only near-boundary occurrences can do.
            // the bounds of q's sequence, computed once per record (every candidate offset of the record is checked against them)
            uint64_t lo = 0, hi = 0;
            bool have = false;
            auto inside = [&](uint64_t s) {
                if (!near) return true;
                if (!have) { seq_bounds(q, lo, hi); have = true; }
                return s >= lo && s + k_ <= hi;
            };
            // The candidate offsets are the same set for x and rc(x), since the minimizer positions mirror, so the smallest matching offset, rather than the first one tried, gives the answer.
            // That makes it, and with it --ids, strand-independent when a record holds a k-mer twice (a hairpin with x and rc(x) a few bases apart).
            uint64_t best = kNotFound;
            bool best_rev = false;
            for (unsigned t = 0; t < mi.npos; ++t) {
                const uint64_t p = mi.pos[t];
                if (q >= p) {                                    // x forward
                    const uint64_t s = q - p;
                    if (s < best && get_mer(s, k_) == x && inside(s)) { best = s; best_rev = false; }
                }
                const uint64_t pr = k_ - m_ - p;                 // x reverse
                if (q >= pr) {
                    const uint64_t s = q - pr;
                    if (s < best && get_mer(s, k_) == xr && inside(s)) { best = s; best_rev = true; }
                }
            }
            if (best != kNotFound && rev_out) *rev_out = best_rev;
            return best;
        };
        if (!(e & 1)) return try_rec(e >> 1);
        uint64_t off = e >> 1;
        unsigned cw;
        const uint64_t cnt = ovf_count(off, cw);
        off += cw;
        const unsigned skb = lo_bits_ + 1;
        if (cnt > kSkewMin) {
            // In a heavy bucket the k-mer's record index comes from its size class's open-addressing table (fingerprint-filtered linear probing).
            const unsigned c = cls_of(cnt);
            const uint64_t cap = skew_cap_[c];
            const uint64_t* W = skew_tables_[c];
            const unsigned w = c + kFpBits;
            const uint64_t h = skew_hash(x < xr ? x : xr);
            const uint64_t fp = skew_fp(h), imask = (uint64_t(1) << c) - 1;
            for (uint64_t sl = home(h, cap);; sl = (sl + 1 == cap) ? 0 : sl + 1) {
                const uint64_t v = get_bits(W, sl * w, w);
                if (v == 0) return kNotFound;                  // end of the chain
                if ((v >> c) == fp && (v & imask) < cnt) {
                    const uint64_t s = try_rec(get_bits(ovf_view_, off + (v & imask) * skb, skb));
                    if (s != kNotFound) return s;
                }
            }
        }
        // A bucket's records point into the string at unrelated places, so checking them one after another serialises one cache miss each.
        // Records are decoded a chunk at a time and, when there is more than one, the string lines of the whole chunk are prefetched before the first is compared, so their misses overlap (in the same order, with the same answer, the first record that holds the k-mer).
        constexpr uint64_t chunk = kProbeChunk;
        uint64_t recs[chunk];
        for (uint64_t base = 0; base < cnt; base += chunk) {
            const uint64_t n = std::min(chunk, cnt - base);
            for (uint64_t i = 0; i < n; ++i) recs[i] = get_bits(ovf_view_, off + (base + i) * skb, skb);
            if (n > 1)
                for (uint64_t i = 0; i < n; ++i) {
                    const uint64_t q = recs[i] >> 1;
                    __builtin_prefetch(&str_view_[(q > k_ ? q - k_ : 0) >> 5]);
                    __builtin_prefetch(&str_view_[(q + k_) >> 5]);
                }
            for (uint64_t i = 0; i < n; ++i) {
                const uint64_t s = try_rec(recs[i]);
                if (s != kNotFound) return s;
            }
        }
        return kNotFound;
    }

    unsigned k_ = 0, m_ = 0, w_ = 0;
    uint64_t kmask_ = 0, mmask_ = 0;
    unsigned lo_bits_ = 1, rec_bits_ = 3;
    std::size_t n_starts_ = 0;     // number of packed sequence starts (seqs + 1)
    // sequence starts as an Elias-Fano sequence
    EfSeq ef_starts_;
    uint64_t start_at(std::size_t i) const noexcept { return ef_starts_[i]; }

    /// seq_index() - index of the sequence containing string offset @pos.
    std::size_t seq_index(uint64_t pos) const noexcept {
        return static_cast<std::size_t>(ef_starts_.count_le(pos)) - 1;
    }
    /// seq_bounds() - [@lo, @hi) of the sequence holding base offset @pos, one Elias-Fano rank and two accesses.
    void seq_bounds(uint64_t pos, uint64_t& lo, uint64_t& hi) const noexcept {
        const std::size_t a = seq_index(pos);
        lo = start_at(a);
        hi = start_at(a + 1);
    }
    // Skew tables for heavy buckets, whose three parameters are stored in the index file and checked at load.
    static constexpr uint64_t kSkewMin = 64;        // records above which a bucket gets a skew table
    static constexpr unsigned kFpBits = 5;          // fingerprint bits per skew slot
    static constexpr unsigned kMaxClass = 40;       // largest size class, buckets of up to 2^40 records
    // MPHF parameters (Mphf::build) and its size per key for estimate_m()
    static constexpr double kMphfC = 6.0, kMphfAlpha = 0.94;
    static constexpr double kMphfBitsPerKey = 5.1;
    uint64_t skew_entries_ = 0;
    uint64_t distinct_ = 0;           // distinct canonical k-mers counted by build
    bool has_distinct_ = false;
    std::vector<std::vector<uint64_t>> skew_own_;
    std::vector<uint64_t> skew_cap_;
    std::vector<const uint64_t*> skew_tables_;
    std::vector<std::size_t> skew_words_;
    /// cls_of() - size class of a bucket of @cnt records, i.e. the bits for an index < cnt.
    static constexpr unsigned cls_of(uint64_t cnt) noexcept {
        return cnt <= 1 ? 1u : 64u - static_cast<unsigned>(__builtin_clzll(cnt - 1));
    }
    /// skew_hash() - SplitMix64's finalizer applied to @c ^ 2^64/phi.
    static constexpr uint64_t skew_hash(uint64_t c) noexcept {
        c ^= 0x9E3779B97F4A7C15ULL;
        c ^= c >> 31; c *= 0xbf58476d1ce4e5b9ULL;
        c ^= c >> 27; c *= 0x94d049bb133111ebULL;
        c ^= c >> 31;
        return c;
    }
    /// skew_fp() - fingerprint of hash @h, in [1, 2^kFpBits); 0 marks an empty slot.
    static constexpr uint64_t skew_fp(uint64_t h) noexcept {
        return 1 + (h & 0xFFFF) % ((uint64_t(1) << kFpBits) - 1);
    }
    static constexpr uint64_t home(uint64_t h, uint64_t cap) noexcept {   // fast range
        return static_cast<uint64_t>((static_cast<__uint128_t>(h) * cap) >> 64);
    }
    static void set_bits(uint64_t* v, uint64_t pos, uint64_t val, unsigned nbits) noexcept {
        v[pos >> 6] |= val << (pos & 63);
        if ((pos & 63) + nbits > 64) v[(pos >> 6) + 1] |= val >> (64 - (pos & 63));
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
