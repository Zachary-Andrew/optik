#pragma once

/**
 * @file
 * @brief SeedIndex, which maps minimizers to reference positions for the mapper.
 *
 *   keys      distinct canonical minimizers, in minimizer.hpp's hashed order
 *   mphf_     partitioned Mphf (phobic_mphf.hpp) that maps each key to a rank r
 *   off_[r]   uint32 CSR offsets into post_
 *   post_     postings (see make_posting()), grouped by rank and kept in reference order within a rank
 *
 * The reference itself serves as the key store.
 * A query minimizer's rank is accepted only if the canonical reference k-mer at its first posting equals it.
 * So no keys are stored, and a foreign minimizer costs one MPHF access plus one reference read.
 *
 * collect() drops minimizers occurring more than max_occ times, like minimap2's -f filter (Li 2018).
 * collect_chain(), on the other hand, rescues some of them (Li 2021).
 * Both seed a whole read in one call and emit strand-aware hits.
 */

#include "bits.hpp"
#include "thread_budget.hpp"
#include <algorithm>
#include <atomic>
#include <malloc.h>
#include <cstdint>
#include <string>
#include <thread>
#include <vector>

#include "hugepage.hpp"
#include "phase_timer.hpp"
#include "kmer.hpp"
#include "minimizer.hpp"
#include "phobic_mphf.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

class SeedIndex {
public:
    /// Strand 1 means the read is reverse complemented, and qpos is measured in that strand's frame.
    /// The alignment-start diagonal is therefore rpos - qpos.
    struct Hit { uint32_t rid, rev, rpos, qpos; };

    /// @p refs must stay alive and unmodified for as long as the index exists.
    /// The build fails for a reference of kMaxRefLen bases or more, for kMaxPostings or more postings, or when the MPHF build fails.
    ///
    /// Although the build is parallel, it produces exactly what a serial build would:
    ///  1. Each reference is cut into tasks of kExtractChunk window starts.
    ///     A task emits only the minimizer positions it owns, sorted into kNpart hash partitions.
    ///  2. A single thread collects each partition's distinct keys, which works because a key never spans partitions.
    ///  3. One MPHF then covers all keys, built partition by partition in parallel.
    ///  4. Finally, a counting sort on rank fills the CSR, again per partition.
    bool build(const std::vector<std::pair<const char*, std::size_t>>& refs,
               unsigned kmer_len, unsigned window, uint32_t max_occ, unsigned nthreads = 0) {
        kmer_len_ = kmer_len; window_ = window; max_occ_ = max_occ; refs_ = refs;
        PhaseTimer timer;

        for (const auto& ref : refs) if (ref.second >= kMaxRefLen) return false;

        nthreads = resolve_threads(nthreads);

        struct Task { uint32_t rid; std::size_t win_start, win_end; };     // window starts [start, end)
        constexpr std::size_t kExtractChunk = std::size_t(1) << 20;
        std::vector<Task> tasks;

        for (std::size_t rid = 0; rid < refs.size(); ++rid) {
            const std::size_t len = refs[rid].second;

            if (len < kmer_len + window - 1) continue;

            const std::size_t nwin = len - (kmer_len + window - 1) + 1;

            for (std::size_t chunk_start = 0; chunk_start < nwin; chunk_start += kExtractChunk)
                tasks.push_back({static_cast<uint32_t>(rid), chunk_start, std::min(nwin, chunk_start + kExtractChunk)});
        }

        nthreads = std::max(1u, std::min<unsigned>(nthreads, static_cast<unsigned>(std::max<std::size_t>(tasks.size(), 1))));
        // 1. Every task owns the minimizer positions [start, end), and the last task of a reference owns through its end.
        //    It scans windows starting in [start - (w-1), end), since any of those can select an owned position.
        //    Emitting only owned positions means no seed appears twice and none goes missing.
        //    Per-task buffers keep postings in reference order.
        //    Only the 8-byte posting is kept, because key_of() can re-read the key from the reference whenever needed.
        //    Carrying the key as well doubled the largest buffers, and hg38's 535M postings then peaked at 25 GB for a 6 GB index.
        //    Buckets are sized from the expected minimizer density, which keeps the doubling slack small.
        std::vector<std::vector<std::vector<uint64_t>>> task_buckets(tasks.size(), std::vector<std::vector<uint64_t>>(kNpart));

        parallel(nthreads, tasks.size(), [&](unsigned, std::size_t task_idx) {
            const Task& task = tasks[task_idx];
            const char* seq = refs[task.rid].first;
            const std::size_t len = refs[task.rid].second;
            const std::size_t nwin = len - (kmer_len + window - 1) + 1;
            const std::size_t scan_start = task.win_start >= window - 1 ? task.win_start - (window - 1) : 0;
            const std::size_t end = task.win_end + (window - 1) + (kmer_len - 1);   // exclusive base end
            const std::size_t own_end = task.win_end == nwin ? len : task.win_end;
            auto& out = task_buckets[task_idx];
            const std::size_t expect = kMinimizerDensityNumerator * (task.win_end - task.win_start) / (window + 1) / kNpart;

            for (auto& bucket : out) bucket.reserve(expect + expect / kBucketSlackDiv + kBucketSlack);

            extract_minimizers(seq + scan_start, end - scan_start, kmer_len, window, [&](const MinimizerSeed& seed) {
                const uint64_t pos = scan_start + seed.pos;

                if (pos < task.win_start || pos >= own_end) return;

                out[part_of(seed.kmer)].push_back(make_posting(task.rid, pos, seed.rev));
            });
        });
        timer.lap("seed: extract minimizers");

        // 2. Each partition's distinct keys are read straight from the per-task buckets, in reference order, through a cache-sized hash set.
        //    Sorting is unnecessary, because the MPHF ranks do the grouping below.
        //    There is no gathered copy either, because freed small buckets stay resident in the allocator and a copy cost its full size again (4.3 GB on hg38).
        auto each = [&](std::size_t part, auto&& func) {
            for (const auto& task : task_buckets) for (const uint64_t post : task[part]) func(post);
        };
        // each() serves the passes that re-read keys.
        // A bucket's postings sit about w/2 x kNpart bases apart, so every key read misses cache.
        // To hide that, the reference bytes are prefetched kPrefetchDistance postings ahead.
        auto each_key = [&](std::size_t part, auto&& func) {
            for (const auto& task : task_buckets) {
                const auto& bucket = task[part];
                const std::size_t count = bucket.size();

                for (std::size_t i = 0; i < std::min(count, kPrefetchDistance); ++i) prefetch_ref(bucket[i]);

                for (std::size_t i = 0; i < count; ++i) {
                    if (i + kPrefetchDistance < count) prefetch_ref(bucket[i + kPrefetchDistance]);

                    func(bucket[i]);
                }
            }
        };
        std::vector<std::size_t> psize(kNpart, 0);
        std::vector<std::vector<uint64_t>> pkeys(kNpart);

        parallel(nthreads, kNpart, [&](unsigned, std::size_t part) {
            std::size_t tot = 0;

            for (auto& task : task_buckets) tot += task[part].size();

            psize[part] = tot;
            std::size_t cap = kMinKeySetSlots;

            while (cap < tot * kKeySetSlotsPerKey) cap <<= 1;

            std::vector<uint64_t> set(cap, ~uint64_t(0));     // ~0 is free, since no canonical k-mer with k <= 32 reaches it
            auto& distinct_keys = pkeys[part];

            each_key(part, [&](uint64_t post) {
                const uint64_t key = key_of(post);
                // part_of() used the top bits, which every key in this partition shares, so an independent hash is needed here.
                std::size_t slot = static_cast<std::size_t>(mix64(key, 0)) & (cap - 1);

                for (;;) {
                    if (set[slot] == key) break;

                    if (set[slot] == ~uint64_t(0)) { set[slot] = key; distinct_keys.push_back(key); break; }

                    slot = (slot + 1) & (cap - 1);
                }
            });
            distinct_keys.shrink_to_fit();       // doubling slack across all partitions reached 3.3 GB on hg38
        });
        num_keys_ = 0;

        for (const auto& keys : pkeys) num_keys_ += keys.size();

        timer.lap("seed: distinct keys");

        // Build speed takes priority over MPHF space here, because the 8-byte postings dominate the index anyway.
        // When this was chosen, c = 16 built 17.5M keys in 1.1 s per thread at 8 bits per key, versus 4.2 s and 4.8 bits per key at the size-tuned c = 6.
        // The per-partition key lists go in directly, without a concatenated copy.
        if (num_keys_ && !mphf_.build(pkeys, kMphfC, kMphfAlpha, nthreads)) return false;

        timer.lap("seed: mphf");

        // 3. CSR by counting sort on rank.
        //    Each posting's rank is computed once and kept for the scatter pass.
        //    Computing it re-reads the key from the reference, with the MPHF line prefetched kPrefetchDistance postings ahead.
        //    The count and scatter passes prefetch their random off_ and post_ lines in the same way (Slotin, "Prefetching").
        //    The scatter uses off_ itself as the write cursor and shifts it back afterwards, since a separate cursor took 1.7 GB on hg38.
        //    All postings of a key are in one partition and are scattered in reference order, so each posting list ends up sorted.
        // Because off_ is 32-bit, kMaxPostings or more postings (a reference near 25 Gbp) would wrap silently, so the build refuses instead.
        {
            std::size_t total = 0;

            for (std::size_t part_size : psize) total += part_size;

            if (total >= kMaxPostings) return false;
        }
        off_.assign(num_keys_ + 1, 0);
        std::vector<std::vector<uint32_t>> prank(kNpart);

        parallel(nthreads, kNpart, [&](unsigned, std::size_t part) {
            auto& ranks = prank[part];

            ranks.resize(psize[part]);
            std::size_t cursor = 0;

            for (const auto& task : task_buckets) {
                const auto& bucket = task[part];
                const std::size_t count = bucket.size();
                // Each key is re-read exactly once.
                // Reference bytes are prefetched 2 * kPrefetchDistance postings ahead, and the key itself is read kPrefetchDistance ahead.
                // Its MPHF line is prefetched at that point, and the key then waits in a ring until it's used.
                uint64_t ring[kPrefetchDistance];
                constexpr std::size_t kRefPrefetchDistance = 2 * kPrefetchDistance;

                for (std::size_t i = 0; i < std::min(count, kRefPrefetchDistance); ++i) prefetch_ref(bucket[i]);

                for (std::size_t i = 0; i < std::min(count, kPrefetchDistance); ++i) {
                    ring[i] = key_of(bucket[i]);
                    mphf_.prefetch(ring[i]);
                }

                for (std::size_t i = 0; i < count; ++i) {
                    const uint64_t key = ring[i % kPrefetchDistance];

                    if (i + kRefPrefetchDistance < count) prefetch_ref(bucket[i + kRefPrefetchDistance]);

                    if (i + kPrefetchDistance < count) {
                        ring[i % kPrefetchDistance] = key_of(bucket[i + kPrefetchDistance]);
                        mphf_.prefetch(ring[i % kPrefetchDistance]);
                    }

                    // The MPHF rank of each posting's key decides which list it joins.
                    ranks[cursor++] = static_cast<uint32_t>(mphf_(key));
                }
            }

            for (std::size_t i = 0; i < ranks.size(); ++i) {
                if (i + kPrefetchDistance < ranks.size()) prefetch_write_once(&off_[ranks[i + kPrefetchDistance] + 1]);

                // Counting into off_[r + 1] makes the prefix sum below yield each list's start in off_[r].
                ++off_[ranks[i] + 1];
            }
        });
        for (std::size_t r = 0; r < num_keys_; ++r) off_[r + 1] += off_[r];

        post_.resize(off_[num_keys_]);
        parallel(nthreads, kNpart, [&](unsigned, std::size_t part) {
            const auto& ranks = prank[part];
            std::size_t cursor = 0;

            each(part, [&](uint64_t post) {
                if (cursor + kPrefetchDistance < ranks.size()) prefetch_write_once(&off_[ranks[cursor + kPrefetchDistance]]);

                // Scattering partitions in reference order leaves every posting list sorted by position with no extra sort.
                post_[off_[ranks[cursor]]++] = post;
                ++cursor;
            });
            std::vector<uint32_t>().swap(prank[part]);
        });
        // off_[r] now holds the end of rank r's list, so everything shifts back by one.
        for (std::size_t r = num_keys_; r > 0; --r) off_[r] = off_[r - 1];

        off_[0] = 0;
        task_buckets.clear();
        // The allocator keeps freed small buckets, so they are handed back explicitly.
        // Otherwise the mapper would run with 4.3 GB of them resident on hg38.
        malloc_trim(0);
        timer.lap("seed: csr fill");
        set_mid_occ();
        return true;
    }

    /// collect_chain() skips minimizers above this occurrence count, which a kMidOccFrac fraction of distinct minimizers exceed.
    /// The count is clamped to [kMinMidOcc, max_occ].
    /// Deriving it from the index's own distribution lets it track the genome's repeat content (minimap2's mid_occ, Li 2021).
    static constexpr double kMidOccFrac = 2e-4;
    static constexpr uint32_t kMinMidOcc = 10;
    /// Within a run of over-threshold minimizers, collect_chain() keeps the least frequent one per kOccDist read bases, as long as it occurs at most kMaxMaxOcc times.
    static constexpr uint32_t kOccDist = 500;
    static constexpr uint32_t kMaxMaxOcc = 4095;
    uint32_t mid_occ() const noexcept { return mid_occ_; }

    /// Minimizers occurring at most mid_occ() times seed exactly as in collect().
    /// A run of over-threshold minimizers, meaning the read lies in a repeat, isn't dropped whole, though.
    /// Between the flanking kept minimizers (or read ends) at ps and pe, the round((pe - ps) / kOccDist) least frequent ones still seed if they occur at most kMaxMaxOcc times.
    /// A read inside a satellite array thus keeps anchors at its true copy and at its competitors, and MAPQ can account for both.
    /// *rep_len receives the read bases covered only by dropped minimizers, as the union of their [pos, pos + k) intervals.
    void collect_chain(const char* seq, std::size_t len, std::vector<Hit>& out,
                       uint32_t* rep_len = nullptr) const {
        out.clear();

        if (rep_len) *rep_len = 0;

        if (num_keys_ == 0) return;

        thread_local std::vector<MinimizerSeed> seeds;
        thread_local std::vector<uint64_t> ranks;
        thread_local std::vector<uint32_t> occ;
        thread_local std::vector<uint8_t> keep;
        thread_local std::vector<uint32_t> cand;

        seeds.clear();
        extract_minimizers(seq, len, kmer_len_, window_, [&](const MinimizerSeed& seed) { seeds.push_back(seed); });
        const std::size_t count = seeds.size();

        resolve(seeds, ranks);
        occ.assign(count, 0);
        keep.assign(count, 0);

        for (std::size_t i = 0; i < count; ++i) {
            const uint32_t post_lo = range_lo(ranks[i]), post_hi = range_hi(ranks[i]);

            // The MPHF also maps keys it was never built on. Comparing the first posting's reference k-mer confirms the minimizer is in the index.
            if (post_lo == post_hi || !verify(post_[post_lo], seeds[i].kmer)) continue;

            occ[i] = post_hi - post_lo;
            // Minimizers above mid_occ_ are dropped at first, and the loop below restores the least frequent one in each uncovered run.
            keep[i] = occ[i] <= mid_occ_;
        }

        // Find runs of over-threshold minimizers, bounded by kept ones or by the read ends.
        std::size_t run_start = 0;

        while (run_start < count) {
            if (occ[run_start] <= mid_occ_) { ++run_start; continue; }

            std::size_t run_end = run_start;

            while (run_end < count && occ[run_end] > mid_occ_) ++run_end;

            const std::size_t prev_pos = run_start > 0 ? seeds[run_start - 1].pos : 0;
            const std::size_t next_pos = run_end < count ? seeds[run_end].pos : len;
            const std::size_t want = static_cast<std::size_t>(double(next_pos - prev_pos) / kOccDist + kWantRoundingBias);

            if (want > 0) {
                cand.clear();

                for (std::size_t t = run_start; t < run_end; ++t)
                    if (occ[t] <= kMaxMaxOcc) cand.push_back(static_cast<uint32_t>(t));

                if (cand.size() > want) {
                    // Keep the `want` least frequent ones, breaking ties toward the earlier ones in the read.
                    std::nth_element(cand.begin(), cand.begin() + want, cand.end(), [&](uint32_t lhs, uint32_t rhs) {
                        return occ[lhs] != occ[rhs] ? occ[lhs] < occ[rhs] : lhs < rhs;
                    });
                    cand.resize(want);
                }

                for (uint32_t seed_idx : cand) keep[seed_idx] = 1;
            }

            run_start = run_end;
        }

        uint32_t rep_bases = 0;
        std::size_t cov_end = 0;                              // union of dropped [pos, pos + k)

        for (std::size_t t = 0; t < count; ++t) {
            if (occ[t] == 0) continue;

            if (!keep[t]) {
                const std::size_t cov_start = std::max<std::size_t>(seeds[t].pos, cov_end), cov_stop = seeds[t].pos + kmer_len_;

                if (cov_stop > cov_start) { rep_bases += static_cast<uint32_t>(cov_stop - cov_start); cov_end = cov_stop; }

                continue;
            }

            const uint32_t post_lo = range_lo(ranks[t]), post_hi = range_hi(ranks[t]);

            emit(seeds[t], post_lo, post_hi, len, out);
        }

        if (rep_len) *rep_len = rep_bases;
    }

    std::size_t distinct() const noexcept { return num_keys_; }
    std::size_t occurrences() const noexcept { return post_.size(); }

    std::size_t bytes() const noexcept {
        return mphf_.bytes() + off_.size() * sizeof(uint32_t) + post_.size() * sizeof(uint64_t);
    }

    /// Returns {nullptr, 0} for an absent key.
    std::pair<const uint64_t*, uint32_t> find(uint64_t key) const noexcept {
        if (num_keys_ == 0) return {nullptr, 0};

        const uint64_t rank = mphf_(key);
        const uint32_t post_lo = off_[rank], post_hi = off_[rank + 1];

        if (post_lo == post_hi) return {nullptr, 0};

        if (!verify(post_[post_lo], key)) return {nullptr, 0};

        return {post_.data() + post_lo, post_hi - post_lo};
    }

    /// Gives the same hits as calling find() per minimizer in read order, on both strands.
    /// Each find() is a chain of four dependent misses: pilot, off_[r], post_[lo] and the reference k-mer.
    /// So instead, all minimizers advance stage by stage, with each stage's lines prefetched for all of them first (group prefetching, Slotin).
    /// That exposes about 4 miss latencies per read rather than about 4 per minimizer.
    /// Optionally, @p n_min receives the read's minimizer count and @p n_rep how many were dropped for occurring more than max_occ times.
    /// A read whose seeds are mostly dropped lies in a repeat whose other copies were never candidates, which is why callers want that count.
    void collect(const char* seq, std::size_t len, std::vector<Hit>& out,
                 uint32_t* n_min = nullptr, uint32_t* n_rep = nullptr) const {
        out.clear();

        if (n_min) *n_min = 0;

        if (n_rep) *n_rep = 0;

        if (num_keys_ == 0) return;

        thread_local std::vector<MinimizerSeed> seeds;
        thread_local std::vector<uint64_t> ranks;

        seeds.clear();
        extract_minimizers(seq, len, kmer_len_, window_, [&](const MinimizerSeed& seed) { seeds.push_back(seed); });
        const std::size_t count = seeds.size();

        resolve(seeds, ranks);
        uint32_t rep = 0;

        for (std::size_t i = 0; i < count; ++i) {
            const MinimizerSeed& seed = seeds[i];
            const uint32_t post_lo = range_lo(ranks[i]), post_hi = range_hi(ranks[i]);

            if (post_lo == post_hi || !verify(post_[post_lo], seed.kmer)) continue;

            // Skipped but counted, because hidden_odds weighs their share.
            if (post_hi - post_lo > max_occ_) { ++rep; continue; }

            emit(seed, post_lo, post_hi, len, out);
        }

        if (n_min) *n_min = static_cast<uint32_t>(count);

        if (n_rep) *n_rep = rep;
    }

private:
    /// Fills ranks[i] with minimizer i's packed posting range (see pack_range()), staged and prefetched as collect() describes.
    void resolve(const std::vector<MinimizerSeed>& seeds, std::vector<uint64_t>& ranks) const {
        const std::size_t count = seeds.size();

        ranks.resize(count);

        // Each of the four passes prefetches what the next pass reads, letting a read's cache misses overlap instead of happening one after another.
        for (std::size_t i = 0; i < count; ++i) mphf_.prefetch(seeds[i].kmer);

        for (std::size_t i = 0; i < count; ++i) {
            ranks[i] = mphf_(seeds[i].kmer);
            __builtin_prefetch(&off_[ranks[i]]);
        }

        for (std::size_t i = 0; i < count; ++i) {
            const uint32_t post_lo = off_[ranks[i]];

            ranks[i] = pack_range(post_lo, off_[ranks[i] + 1]);
            __builtin_prefetch(&post_[post_lo]);
        }

        for (std::size_t i = 0; i < count; ++i) {
            const uint32_t post_lo = range_lo(ranks[i]), post_hi = range_hi(ranks[i]);

            if (post_lo == post_hi) continue;

            const uint64_t post = post_[post_lo];
            const char* ref_ptr = refs_[posting_ref(post)].first + posting_pos(post);

            __builtin_prefetch(ref_ptr);
            __builtin_prefetch(ref_ptr + kmer_len_ - 1);
        }
    }

    void emit(const MinimizerSeed& seed, uint32_t post_lo, uint32_t post_hi, std::size_t len, std::vector<Hit>& out) const {
        for (uint32_t t = post_lo; t < post_hi; ++t) {
            const uint64_t post = post_[t];
            const uint32_t rid = static_cast<uint32_t>(posting_ref(post));
            const uint32_t rpos = static_cast<uint32_t>(posting_pos(post));
            // A hit is reverse when the read and the reference saw the minimizer on opposite strands.
            const bool rev = posting_reverse(post) != seed.rev;
            // Chaining works on the reverse-complemented read for reverse hits.
            const uint32_t read_pos = rev ? static_cast<uint32_t>(len - kmer_len_ - seed.pos) : seed.pos;

            out.push_back({rid, rev ? 1u : 0u, rpos, read_pos});
        }
    }

    /// Works from a histogram of occurrence counts, because a per-key array would be huge for hg38's hundreds of millions of keys.
    void set_mid_occ() {
        const uint32_t cap = std::max<uint32_t>(max_occ_, 1);

        mid_occ_ = std::min(kMinMidOcc, cap);

        if (num_keys_ == 0) return;

        // Bins 0..cap hold exact counts, and bin cap + 1 collects everything above.
        std::vector<uint64_t> hist(static_cast<std::size_t>(cap) + kHistogramOverflowBins, 0);

        for (std::size_t r = 0; r < num_keys_; ++r) ++hist[std::min<uint32_t>(off_[r + 1] - off_[r], cap + 1)];

        const uint64_t above = static_cast<uint64_t>(kMidOccFrac * double(num_keys_));
        uint64_t cum = 0;
        uint32_t threshold = cap + 1;

        // Find the smallest t with at most `above` keys occurring more than t times.
        while (threshold > 0 && cum + hist[threshold] <= above) { cum += hist[threshold]; --threshold; }

        mid_occ_ = std::max(std::min(kMinMidOcc, cap), std::min(threshold, cap));
    }

    /// Confirms that the reference spells @p key at posting @p post.
    bool verify(uint64_t post, uint64_t key) const noexcept {
        const auto& ref = refs_[posting_ref(post)];
        const std::size_t pos = posting_pos(post);
        uint64_t kmer;

        if (!encode_kmer(ref.first + pos, kmer_len_, kmer)) return false;

        const uint64_t rev_comp = reverse_complement(kmer, kmer_len_);

        return (kmer < rev_comp ? kmer : rev_comp) == key;
    }

    void prefetch_ref(uint64_t post) const noexcept {
        const char* ref_ptr = refs_[posting_ref(post)].first + posting_pos(post);

        prefetch_read_once(ref_ptr);
        prefetch_read_once(ref_ptr + kmer_len_ - 1);    // a k-mer may cross a line
    }

    /// Reads the minimizer key back out of the reference.
    ///
    /// A posting always points at k valid bases, because minimizers never span an N.
    /// That allows decoding without a table or validity check, since ((c >> 1) ^ (c >> 2)) & 3 maps A, C, G, T to 0..3 in either case.
    uint64_t key_of(uint64_t post) const noexcept {
        const auto& ref = refs_[posting_ref(post)];
        const auto* bases = reinterpret_cast<const unsigned char*>(ref.first + posting_pos(post));
        uint64_t kmer = 0;

        for (unsigned i = 0; i < kmer_len_; ++i)
            kmer = (kmer << kBitsPerBase) | (((bases[i] >> kAsciiDecodeShift1) ^ (bases[i] >> kAsciiDecodeShift2)) & kComplementMask);

        const uint64_t rev_comp = reverse_complement(kmer, kmer_len_);

        return kmer < rev_comp ? kmer : rev_comp;
    }

    /// A posting stores the reference id above bit kRidShift, the position in bits 1..32 and the strand in bit 0.
    static constexpr unsigned kPosShift = 1;
    static constexpr unsigned kRidShift = 33;
    static constexpr uint64_t kPosMask = 0xFFFFFFFFu;
    static constexpr uint64_t kStrandBit = 1;
    /// A rank's posting range is packed as (lo << kRangeShift) | hi.
    static constexpr unsigned kRangeShift = 32;
    /// Shifts that turn an ASCII base into its 2-bit code (see key_of()).
    static constexpr unsigned kAsciiDecodeShift1 = 1, kAsciiDecodeShift2 = 2;
    /// A window of w k-mers yields about 2/(w + 1) minimizers per position.
    static constexpr std::size_t kMinimizerDensityNumerator = 2;
    /// The per-partition key set keeps at most half its slots full.
    static constexpr std::size_t kKeySetSlotsPerKey = 2;
    static constexpr std::size_t kHistogramOverflowBins = 2;   // bin 0 plus the "above cap" bin
    static constexpr std::size_t kPrefetchDistance = 16;   // postings ahead
    static constexpr std::size_t kMinKeySetSlots = 64;
    static constexpr double kMphfC = 16.0, kMphfAlpha = 0.94;

    static constexpr unsigned kNpartBits = 8;
    /// Reference positions get 31 bits and posting offsets 32, so anything at or past these is out of range.
    static constexpr std::size_t kMaxRefLen = std::size_t(1) << 31;
    static constexpr std::size_t kMaxPostings = std::size_t(1) << 32;
    /// Partition buckets reserve an eighth (kBucketSlackDiv) more than the expected minimizers, plus kBucketSlack.
    static constexpr std::size_t kBucketSlackDiv = 8;
    static constexpr std::size_t kBucketSlack = 8;
    /// Added before truncating a repeat run's wanted seed count, so that exact halves round down.
    static constexpr double kWantRoundingBias = 0.499;
    static constexpr std::size_t kNpart = std::size_t(1) << kNpartBits;

    static constexpr uint64_t make_posting(uint64_t rid, uint64_t pos, bool rev) noexcept {
        return (rid << kRidShift) | (pos << kPosShift) | (rev ? kStrandBit : 0);
    }
    static constexpr uint64_t posting_ref(uint64_t post) noexcept { return post >> kRidShift; }
    static constexpr std::size_t posting_pos(uint64_t post) noexcept { return (post >> kPosShift) & kPosMask; }
    static constexpr bool posting_reverse(uint64_t post) noexcept { return (post & kStrandBit) != 0; }
    static constexpr uint64_t pack_range(uint32_t lo, uint32_t hi) noexcept { return (uint64_t(lo) << kRangeShift) | hi; }
    static constexpr uint32_t range_lo(uint64_t range) noexcept { return static_cast<uint32_t>(range >> kRangeShift); }
    static constexpr uint32_t range_hi(uint64_t range) noexcept { return static_cast<uint32_t>(range); }

    /// Takes the top kNpartBits of the key's Fibonacci hash.
    static constexpr std::size_t part_of(uint64_t key) noexcept {
        return static_cast<std::size_t>((key * kGoldenRatio64) >> (kWordBits - kNpartBits));
    }
    /// Runs f(thread, item) for every item in [0, @p count), handing items out through an atomic counter.
    /// The calling thread takes part as well.
    template <typename F>
    static void parallel(unsigned nthreads, std::size_t count, F&& func) {
        std::atomic<std::size_t> next{0};
        auto work = [&](unsigned thread_idx) {
            for (std::size_t item; (item = next.fetch_add(1, std::memory_order_relaxed)) < count;) func(thread_idx, item);
        };
        std::vector<std::thread> threads;

        for (unsigned t = 1; t < nthreads; ++t) threads.emplace_back(work, t);

        work(0);

        for (auto& thread : threads) thread.join();
    }

    unsigned kmer_len_ = 0, window_ = 0;
    uint32_t max_occ_ = 0;
    uint32_t mid_occ_ = 0;
    std::size_t num_keys_ = 0;
    std::vector<std::pair<const char*, std::size_t>> refs_;
    Mphf mphf_;
    HugeVec<uint32_t> off_;       // read at random per seed, hence huge pages
    HugeVec<uint64_t> post_;
};

}
