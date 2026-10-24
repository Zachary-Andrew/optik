#pragma once

/*
 * seed_index.hpp - SeedIndex, minimizer -> reference positions for the mapper.
 *
 *   keys      distinct canonical minimizers in minimizer.hpp's hashed order
 *   mphf_     partitioned Mphf of phobic_mphf.hpp, mapping each key to a rank r
 *   off_[r]   uint32 CSR offsets into post_
 *   post_     postings (rid << kRidShift) | (pos << 1) | rev, grouped by rank and in reference order within a rank
 *
 * The reference is the key store.
 * A query minimizer's rank is accepted only if the canonical reference k-mer at its first posting equals it, so no keys are stored and a foreign minimizer costs one MPHF access plus one reference read.
 * collect() drops minimizers occurring more than max_occ times, as minimap2's -f filter drops repetitive minimizers (Li 2018), while collect_chain() keeps some of them (Li 2021).
 *
 * collect() and collect_chain() seed a whole read in one call and emit strand-aware hits.
 */

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

namespace tpoptoa {

class SeedIndex {
public:
    /// Hit - one seed hit, with reference id, strand (1 = read reverse complemented), reference position, and position in the read in that strand's frame, so the alignment-start diagonal is rpos - qpos.
    struct Hit { uint32_t rid, rev, rpos, qpos; };

    /// build() - index the minimizers of @refs with k-mer length @k and window @w.
    /// @refs must stay alive and unmodified for the lifetime of the index.
    /// Returns false for a reference of 2^31 bases or more, more than 2^32 - 1 postings, or a failed MPHF build.
    ///
    /// The build is parallel, with @nthreads = 0 resolved by resolve_threads(), and gives the same result as a serial build.
    ///  1. Each reference is cut into tasks of kExtractChunk window starts.
    ///     A task owns the minimizer positions in its range and emits only those, into kNpart hash partitions.
    ///  2. One thread collects the distinct keys of each partition, since distinct keys never span partitions.
    ///  3. One MPHF covers all keys, built in partitions in parallel.
    ///  4. A counting sort on rank fills the CSR, per partition.
    bool build(const std::vector<std::pair<const char*, std::size_t>>& refs,
               unsigned k, unsigned w, uint32_t max_occ, unsigned nthreads = 0) {
        k_ = k; w_ = w; max_occ_ = max_occ; refs_ = refs;
        PhaseTimer pt;
        for (const auto& r : refs) if (r.second >= (std::size_t(1) << 31)) return false;
        nthreads = resolve_threads(nthreads);

        struct Task { uint32_t rid; std::size_t s, e; };     // window starts [s, e)
        constexpr std::size_t kExtractChunk = std::size_t(1) << 20;
        std::vector<Task> tasks;
        for (std::size_t rid = 0; rid < refs.size(); ++rid) {
            const std::size_t len = refs[rid].second;
            if (len < k + w - 1) continue;
            const std::size_t nwin = len - (k + w - 1) + 1;
            for (std::size_t s = 0; s < nwin; s += kExtractChunk)
                tasks.push_back({static_cast<uint32_t>(rid), s, std::min(nwin, s + kExtractChunk)});
        }
        nthreads = std::max(1u, std::min<unsigned>(nthreads, static_cast<unsigned>(std::max<std::size_t>(tasks.size(), 1))));
        // 1. Every task owns the minimizer positions [s, e), and the last task of a reference owns through its end.
        //    It scans the windows starting in [s - (w-1), e), which can select an owned position, and emits only owned positions, so no seed is emitted twice and none is lost.
        //    Per-task buffers keep postings in reference order.
        //    Only the 8-byte posting is kept, since key_of re-reads its key, the reference k-mer it points at, when needed.
        //    Carrying the key doubled the build's largest buffers, and hg38 with 535M postings peaked at 25 GB for a 6 GB index.
        //    Buckets are sized from the expected density 2 / (w + 1), so doubling slack stays small.
        std::vector<std::vector<std::vector<uint64_t>>> task_buckets(tasks.size(), std::vector<std::vector<uint64_t>>(kNpart));
        parallel(nthreads, tasks.size(), [&](unsigned, std::size_t i) {
            const Task& tk = tasks[i];
            const char* seq = refs[tk.rid].first;
            const std::size_t len = refs[tk.rid].second;
            const std::size_t nwin = len - (k + w - 1) + 1;
            const std::size_t s0 = tk.s >= w - 1 ? tk.s - (w - 1) : 0;
            const std::size_t end = tk.e + (w - 1) + (k - 1);   // exclusive base end
            const std::size_t own_end = tk.e == nwin ? len : tk.e;
            auto& out = task_buckets[i];
            const std::size_t expect = 2 * (tk.e - tk.s) / (w + 1) / kNpart;
            for (auto& b : out) b.reserve(expect + expect / 8 + 8);
            extract_minimizers(seq + s0, end - s0, k, w, [&](const MinimizerSeed& m) {
                const uint64_t pos = s0 + m.pos;
                if (pos < tk.s || pos >= own_end) return;
                out[part_of(m.kmer)].push_back((uint64_t(tk.rid) << kRidShift) | (pos << 1) | uint64_t(m.rev));
            });
        });
        pt.lap("seed: extract minimizers");

        // 2. Distinct keys per partition are read straight from the per-task buckets in task order, which is reference order, through a cache-sized hash set.
        //    No sort is needed, since grouping is done by MPHF rank below.
        //    There is no gathered copy, because freed small buckets stay resident in the allocator and a copy cost its full size again, 4.3 GB on hg38.
        auto each = [&](std::size_t p, auto&& f) {
            for (const auto& v : task_buckets) for (const uint64_t x : v[p]) f(x);
        };
        // each() serves the passes that re-read keys.
        // A bucket's postings are about w/2 x kNpart bases apart, so every key read misses cache, and the reference bytes are prefetched kPrefetchDistance postings ahead.
        auto each_key = [&](std::size_t p, auto&& f) {
            for (const auto& v : task_buckets) {
                const auto& a = v[p];
                const std::size_t m = a.size();
                for (std::size_t i = 0; i < std::min(m, kPrefetchDistance); ++i) prefetch_ref(a[i]);
                for (std::size_t i = 0; i < m; ++i) {
                    if (i + kPrefetchDistance < m) prefetch_ref(a[i + kPrefetchDistance]);
                    f(a[i]);
                }
            }
        };
        std::vector<std::size_t> psize(kNpart, 0);
        std::vector<std::vector<uint64_t>> pkeys(kNpart);
        parallel(nthreads, kNpart, [&](unsigned, std::size_t p) {
            std::size_t tot = 0;
            for (auto& v : task_buckets) tot += v[p].size();
            psize[p] = tot;
            std::size_t cap = kMinKeySetSlots;
            while (cap < tot * 2) cap <<= 1;
            std::vector<uint64_t> set(cap, ~uint64_t(0));     // keys < 2^64 - 1, which canonical k-mers with k <= 32 are
            auto& distinct_keys = pkeys[p];
            each_key(p, [&](uint64_t post) {
                const uint64_t key = key_of(post);
                // MurmurHash3 fmix64, a hash independent of part_of(), whose top bits all keys here share
                uint64_t z = key ^ (key >> 33);
                z *= 0xff51afd7ed558ccdULL; z ^= z >> 33;
                z *= 0xc4ceb9fe1a85ec53ULL; z ^= z >> 33;
                std::size_t h = static_cast<std::size_t>(z) & (cap - 1);
                for (;;) {
                    if (set[h] == key) break;
                    if (set[h] == ~uint64_t(0)) { set[h] = key; distinct_keys.push_back(key); break; }
                    h = (h + 1) & (cap - 1);
                }
            });
            distinct_keys.shrink_to_fit();       // doubling slack x 256 partitions reached 3.3 GB on hg38
        });
        n_ = 0;
        for (const auto& v : pkeys) n_ += v.size();
        pt.lap("seed: distinct keys");
        // Build speed is preferred over MPHF space here, since the 8-byte postings dominate the index.
        // When this was chosen, c = 16 built 17.5M keys in 1.1 s per thread at 8 bits per key, against 4.2 s and 4.8 bits per key at the size-tuned c = 6.
        // The per-partition key lists go in directly, with no concatenated copy.
        if (n_ && !mphf_.build(pkeys, kMphfC, kMphfAlpha, nthreads)) return false;
        pt.lap("seed: mphf");

        // 3. CSR by counting sort on rank.
        //    Each posting's rank is computed once, re-reading the key from the reference with its MPHF line prefetched kPrefetchDistance postings ahead, and kept for the scatter pass.
        //    The count and scatter passes prefetch their random off_ and post_ lines likewise (Slotin, Algorithms for Modern Hardware, "Prefetching").
        //    The scatter uses off_ itself as the write cursor and shifts it back afterwards, where a cursor copy took 1.7 GB on hg38.
        //    All postings of a key live in one partition and are scattered in reference order, so each posting list is sorted exactly as before.
        // off_ is 32-bit, and more than 2^32 - 1 postings, from a reference of about 25 Gbp, would wrap the offsets silently, so the build refuses instead.
        {
            std::size_t total = 0;
            for (std::size_t t : psize) total += t;
            if (total >= (std::size_t(1) << 32)) return false;
        }
        off_.assign(n_ + 1, 0);
        std::vector<std::vector<uint32_t>> prank(kNpart);
        parallel(nthreads, kNpart, [&](unsigned, std::size_t p) {
            auto& ranks = prank[p];
            ranks.resize(psize[p]);
            std::size_t x = 0;
            for (const auto& v : task_buckets) {
                const auto& a = v[p];
                const std::size_t m = a.size();
                // Each key is re-read once.
                // Reference bytes are prefetched 2 * kPrefetchDistance postings ahead and the key read kPrefetchDistance ahead (its MPHF line then prefetched), after which it waits in a ring until its use.
                uint64_t ring[kPrefetchDistance];
                for (std::size_t i = 0; i < std::min(m, 2 * kPrefetchDistance); ++i) prefetch_ref(a[i]);
                for (std::size_t i = 0; i < std::min(m, kPrefetchDistance); ++i) {
                    ring[i] = key_of(a[i]);
                    mphf_.prefetch(ring[i]);
                }
                for (std::size_t i = 0; i < m; ++i) {
                    const uint64_t key = ring[i % kPrefetchDistance];
                    if (i + 2 * kPrefetchDistance < m) prefetch_ref(a[i + 2 * kPrefetchDistance]);
                    if (i + kPrefetchDistance < m) {
                        ring[i % kPrefetchDistance] = key_of(a[i + kPrefetchDistance]);
                        mphf_.prefetch(ring[i % kPrefetchDistance]);
                    }
                    ranks[x++] = static_cast<uint32_t>(mphf_(key));
                }
            }
            for (std::size_t i = 0; i < ranks.size(); ++i) {
                if (i + kPrefetchDistance < ranks.size()) __builtin_prefetch(&off_[ranks[i + kPrefetchDistance] + 1], 1, 0);
                ++off_[ranks[i] + 1];
            }
        });
        for (std::size_t r = 0; r < n_; ++r) off_[r + 1] += off_[r];
        post_.resize(off_[n_]);
        parallel(nthreads, kNpart, [&](unsigned, std::size_t p) {
            const auto& ranks = prank[p];
            std::size_t x = 0;
            each(p, [&](uint64_t post) {
                if (x + kPrefetchDistance < ranks.size()) __builtin_prefetch(&off_[ranks[x + kPrefetchDistance]], 1, 0);
                post_[off_[ranks[x]]++] = post;
                ++x;
            });
            std::vector<uint32_t>().swap(prank[p]);
        });
        // off_[r] now holds the end of rank r's list, so shift back by one.
        for (std::size_t r = n_; r > 0; --r) off_[r] = off_[r - 1];
        off_[0] = 0;
        task_buckets.clear();
        // The buckets were small allocations the allocator keeps after free, so they are handed back, or the mapper would run with 4.3 GB of them resident on hg38.
        malloc_trim(0);
        pt.lap("seed: csr fill");
        set_mid_occ();
        return true;
    }

    /// mid_occ() - occurrence threshold of the chaining path (collect_chain), the occurrence count exceeded by a kMidOccFrac fraction of the distinct minimizers, clamped to [kMinMidOcc, max_occ].
    /// It comes from the index's own distribution, so it follows the genome's repeat content instead of a fixed constant (the concept is minimap2's mid_occ, Li 2021 btab705).
    static constexpr double kMidOccFrac = 2e-4;
    static constexpr uint32_t kMinMidOcc = 10;
    /// collect_chain keeps, in each run of over-threshold minimizers, the least frequent one per kOccDist read bases, if it occurs at most kMaxMaxOcc times, as btab705 rescues repeat seeds.
    static constexpr uint32_t kOccDist = 500;
    static constexpr uint32_t kMaxMaxOcc = 4095;
    uint32_t mid_occ() const noexcept { return mid_occ_; }

    /// collect_chain() - seed hits for chaining.
    /// Minimizers occurring at most mid_occ() times seed as in collect().
    /// A run of consecutive over-threshold minimizers, i.e. a repeat the read lies in, is not dropped whole.
    /// Between the flanking kept minimizers or read ends at positions ps and pe, the round((pe - ps) / kOccDist) least frequent ones that occur at most kMaxMaxOcc times still seed.
    /// A read inside a satellite array thus keeps anchors at its true copy and at its competitors, which MAPQ then sees.
    /// *rep_len gets the read bases that only dropped minimizers cover, the union of their [pos, pos + k) intervals.
    void collect_chain(const char* seq, std::size_t len, std::vector<Hit>& out,
                       uint32_t* rep_len = nullptr) const {
        out.clear();
        if (rep_len) *rep_len = 0;
        if (n_ == 0) return;
        thread_local std::vector<MinimizerSeed> ms;
        thread_local std::vector<uint64_t> ranks;
        thread_local std::vector<uint32_t> occ;
        thread_local std::vector<uint8_t> keep;
        thread_local std::vector<uint32_t> cand;
        ms.clear();
        extract_minimizers(seq, len, k_, w_, [&](const MinimizerSeed& s) { ms.push_back(s); });
        const std::size_t n = ms.size();
        resolve(ms, ranks);
        occ.assign(n, 0);
        keep.assign(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            const uint32_t lo = uint32_t(ranks[i] >> 32), hi = uint32_t(ranks[i]);
            if (lo == hi || !verify(post_[lo], ms[i].kmer)) continue;
            occ[i] = hi - lo;
            keep[i] = occ[i] <= mid_occ_;
        }
        // runs of over-threshold minimizers, bounded by kept ones or read ends
        std::size_t i = 0;
        while (i < n) {
            if (occ[i] <= mid_occ_) { ++i; continue; }
            std::size_t j = i;
            while (j < n && occ[j] > mid_occ_) ++j;          // [i, j) is one run
            const std::size_t ps = i > 0 ? ms[i - 1].pos : 0;
            const std::size_t pe = j < n ? ms[j].pos : len;
            const std::size_t want = static_cast<std::size_t>(double(pe - ps) / kOccDist + 0.499);
            if (want > 0) {
                cand.clear();
                for (std::size_t t = i; t < j; ++t)
                    if (occ[t] <= kMaxMaxOcc) cand.push_back(static_cast<uint32_t>(t));
                if (cand.size() > want) {
                    // the `want` least frequent, with ties going to the earlier ones in the read
                    std::nth_element(cand.begin(), cand.begin() + want, cand.end(), [&](uint32_t a, uint32_t b) {
                        return occ[a] != occ[b] ? occ[a] < occ[b] : a < b;
                    });
                    cand.resize(want);
                }
                for (uint32_t t : cand) keep[t] = 1;
            }
            i = j;
        }
        uint32_t rl = 0;
        std::size_t cov_end = 0;                              // union of dropped [pos, pos + k)
        for (std::size_t t = 0; t < n; ++t) {
            if (occ[t] == 0) continue;
            if (!keep[t]) {
                const std::size_t a = std::max<std::size_t>(ms[t].pos, cov_end), b = ms[t].pos + k_;
                if (b > a) { rl += static_cast<uint32_t>(b - a); cov_end = b; }
                continue;
            }
            const uint32_t lo = uint32_t(ranks[t] >> 32), hi = uint32_t(ranks[t]);
            emit(ms[t], lo, hi, len, out);
        }
        if (rep_len) *rep_len = rl;
    }

    std::size_t distinct() const noexcept { return n_; }
    std::size_t occurrences() const noexcept { return post_.size(); }
    std::size_t bytes() const noexcept {
        return mphf_.bytes() + off_.size() * 4 + post_.size() * 8;
    }

    /// find() - postings of canonical minimizer @key, or {nullptr, 0} if absent.
    std::pair<const uint64_t*, uint32_t> find(uint64_t key) const noexcept {
        if (n_ == 0) return {nullptr, 0};
        const uint64_t r = mphf_(key);
        const uint32_t lo = off_[r], hi = off_[r + 1];
        if (lo == hi) return {nullptr, 0};
        if (!verify(post_[lo], key)) return {nullptr, 0};
        return {post_.data() + lo, hi - lo};
    }

    /// collect() - seed hits of a read on both strands, identical to calling find() per minimizer in read order.
    /// Every find() is a chain of four dependent misses (pilot, off_[r], post_[lo], reference k-mer), so the read's minimizers are resolved stage by stage with each stage's lines prefetched for all of them first (group prefetching, Slotin), giving about 4 exposed miss latencies per read instead of about 4 per minimizer.
    /// @n_min and @n_rep (optional) receive the read's minimizers and how many of them occur in the reference but more than max_occ times (dropped as repetitive, since a read whose seeds are mostly dropped lies in a repeat whose other copies were never candidates).
    void collect(const char* seq, std::size_t len, std::vector<Hit>& out,
                 uint32_t* n_min = nullptr, uint32_t* n_rep = nullptr) const {
        out.clear();
        if (n_min) *n_min = 0;
        if (n_rep) *n_rep = 0;
        if (n_ == 0) return;
        thread_local std::vector<MinimizerSeed> ms;
        thread_local std::vector<uint64_t> ranks;
        ms.clear();
        extract_minimizers(seq, len, k_, w_, [&](const MinimizerSeed& s) { ms.push_back(s); });
        const std::size_t n = ms.size();
        resolve(ms, ranks);
        uint32_t rep = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const MinimizerSeed& s = ms[i];
            const uint32_t lo = uint32_t(ranks[i] >> 32), hi = uint32_t(ranks[i]);
            if (lo == hi || !verify(post_[lo], s.kmer)) continue;
            if (hi - lo > max_occ_) { ++rep; continue; }
            emit(s, lo, hi, len, out);
        }
        if (n_min) *n_min = static_cast<uint32_t>(n);
        if (n_rep) *n_rep = rep;
    }

private:
    /// resolve() - ranks[i] = (lo << 32) | hi, the posting range of minimizer i, resolved stage by stage with each stage's lines prefetched for all minimizers first (see collect()).
    void resolve(const std::vector<MinimizerSeed>& ms, std::vector<uint64_t>& ranks) const {
        const std::size_t n = ms.size();
        ranks.resize(n);
        for (std::size_t i = 0; i < n; ++i) mphf_.prefetch(ms[i].kmer);
        for (std::size_t i = 0; i < n; ++i) {
            ranks[i] = mphf_(ms[i].kmer);
            __builtin_prefetch(&off_[ranks[i]]);
        }
        for (std::size_t i = 0; i < n; ++i) {
            const uint32_t lo = off_[ranks[i]];
            ranks[i] = (uint64_t(lo) << 32) | off_[ranks[i] + 1];        // (lo, hi)
            __builtin_prefetch(&post_[lo]);
        }
        for (std::size_t i = 0; i < n; ++i) {
            const uint32_t lo = uint32_t(ranks[i] >> 32), hi = uint32_t(ranks[i]);
            if (lo == hi) continue;
            const uint64_t p = post_[lo];
            const char* r = refs_[p >> kRidShift].first + ((p >> 1) & kPosMask);
            __builtin_prefetch(r);
            __builtin_prefetch(r + k_ - 1);
        }
    }

    /// emit() - strand-aware hits of minimizer @s over postings [@lo, @hi).
    void emit(const MinimizerSeed& s, uint32_t lo, uint32_t hi, std::size_t len, std::vector<Hit>& out) const {
        for (uint32_t t = lo; t < hi; ++t) {
            const uint64_t e = post_[t];
            const uint32_t rid = static_cast<uint32_t>(e >> kRidShift);
            const uint32_t rpos = static_cast<uint32_t>((e >> 1) & kPosMask);
            const bool rev = ((e & 1) != 0) != s.rev;
            const uint32_t qp = rev ? static_cast<uint32_t>(len - k_ - s.pos) : s.pos;
            out.push_back({rid, rev ? 1u : 0u, rpos, qp});
        }
    }

    /// set_mid_occ() - mid_occ_ from a histogram of the keys' occurrence counts, with no per-key array, since hg38 has hundreds of millions of keys.
    void set_mid_occ() {
        const uint32_t cap = std::max<uint32_t>(max_occ_, 1);
        mid_occ_ = std::min(kMinMidOcc, cap);
        if (n_ == 0) return;
        std::vector<uint64_t> hist(static_cast<std::size_t>(cap) + 2, 0);
        for (std::size_t r = 0; r < n_; ++r) ++hist[std::min<uint32_t>(off_[r + 1] - off_[r], cap + 1)];
        const uint64_t above = static_cast<uint64_t>(kMidOccFrac * double(n_));
        uint64_t cum = 0;
        uint32_t t = cap + 1;
        // the smallest t with at most `above` keys occurring more than t times
        while (t > 0 && cum + hist[t] <= above) { cum += hist[t]; --t; }
        mid_occ_ = std::max(std::min(kMinMidOcc, cap), std::min(t, cap));
    }

    /// verify() - the reference as key store, where posting @p must spell @key.
    bool verify(uint64_t p, uint64_t key) const noexcept {
        const auto& ref = refs_[p >> kRidShift];
        const std::size_t pos = (p >> 1) & kPosMask;
        uint64_t x;
        if (!encode_kmer(ref.first + pos, k_, x)) return false;
        const uint64_t xr = reverse_complement(x, k_);
        return (x < xr ? x : xr) == key;
    }

    /// prefetch_ref() - prefetch the reference bytes posting @p points at.
    void prefetch_ref(uint64_t p) const noexcept {
        const char* s = refs_[p >> kRidShift].first + ((p >> 1) & kPosMask);
        __builtin_prefetch(s, 0, 0);
        __builtin_prefetch(s + k_ - 1, 0, 0);    // a k-mer may cross a line
    }

    /// key_of() - the canonical k-mer posting @p points at (its minimizer key).
    ///
    /// A posting always points at k bases of A/C/G/T (minimizers never span an N), so each base is coded without a table or validity check, since ((c >> 1) ^ (c >> 2)) & 3 gives A 0, C 1, G 2, T 3 for upper- and lower-case ASCII, encode_kmer's code.
    uint64_t key_of(uint64_t p) const noexcept {
        const auto& ref = refs_[p >> kRidShift];
        const auto* s = reinterpret_cast<const unsigned char*>(ref.first + ((p >> 1) & kPosMask));
        uint64_t x = 0;
        for (unsigned i = 0; i < k_; ++i) x = (x << 2) | (((s[i] >> 1) ^ (s[i] >> 2)) & 3u);
        const uint64_t xr = reverse_complement(x, k_);
        return x < xr ? x : xr;
    }

    /// Posting layout, with the reference id above bit kRidShift, the position in bits 1..32 and the strand in bit 0.
    static constexpr unsigned kRidShift = 33;
    static constexpr uint64_t kPosMask = 0xFFFFFFFFu;
    static constexpr std::size_t kPrefetchDistance = 16;   // postings ahead
    static constexpr std::size_t kMinKeySetSlots = 64;
    static constexpr double kMphfC = 16.0, kMphfAlpha = 0.94;

    static constexpr std::size_t kNpart = 256;
    /// part_of() - hash partition of @key: the top 8 bits (kNpart = 256) of its Fibonacci hash.
    static constexpr std::size_t part_of(uint64_t key) noexcept {
        return static_cast<std::size_t>((key * 0x9E3779B97F4A7C15ULL) >> 56);
    }
    /// parallel() - run f(thread, item) over [0, @n) with dynamic (atomic counter) scheduling, the calling thread included.
    template <typename F>
    static void parallel(unsigned nthreads, std::size_t n, F&& f) {
        std::atomic<std::size_t> next{0};
        auto work = [&](unsigned t) {
            for (std::size_t i; (i = next.fetch_add(1, std::memory_order_relaxed)) < n;) f(t, i);
        };
        std::vector<std::thread> th;
        for (unsigned t = 1; t < nthreads; ++t) th.emplace_back(work, t);
        work(0);
        for (auto& x : th) x.join();
    }

    unsigned k_ = 0, w_ = 0;
    uint32_t max_occ_ = 0;
    uint32_t mid_occ_ = 0;
    std::size_t n_ = 0;
    std::vector<std::pair<const char*, std::size_t>> refs_;
    Mphf mphf_;
    HugeVec<uint32_t> off_;       // accessed at random per seed, so on huge pages
    HugeVec<uint64_t> post_;
};

}
