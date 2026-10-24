#pragma once

/*
 * par_sort.hpp - deterministic parallel sorts for the index builders.
 *
 * Both sorts give exactly std::sort's result for a strict total order (no two elements compare equal), so a build is byte-identical at every thread count.
 *
 * radix_sort_inplace() distributes with one American-flag pass (in-place MSD radix, where each element is swapped straight into its bucket with no second array), as in McIlroy, Bostic and McIlroy, "Engineering radix sort" (Computing Systems 6(1), 1993).
 * Only the first pass is radix, and buckets are then finished with the comparison sort, one bucket per OpenMP task.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

namespace tpoptoa {

/// Inputs smaller than this are sorted with one std::sort, since below it the thread start-up costs more than the sort.
inline constexpr std::size_t kParSortMinElements = std::size_t(1) << 16;

/// par_sort() - sort @v by @cmp on @nthreads threads.
///
/// The vector is cut into @nthreads equal slices, one sorted per thread, then adjacent runs are merged pairwise in log2(nthreads) parallel rounds.
template <typename V, typename C>
void par_sort(V& v, C cmp, unsigned nthreads) {
    const std::size_t n = v.size();
    const unsigned threads = nthreads ? nthreads : 1;
    if (threads < 2 || n < kParSortMinElements) { std::sort(v.begin(), v.end(), cmp); return; }
    std::vector<std::size_t> cut(threads + 1);
    for (unsigned t = 0; t <= threads; ++t) cut[t] = n * t / threads;
#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (int t = 0; t < int(threads); ++t) std::sort(v.begin() + cut[t], v.begin() + cut[t + 1], cmp);
    for (unsigned width = 1; width < threads; width *= 2) {
        const int pairs = int((threads + 2 * width - 1) / (2 * width));
#ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
        for (int p = 0; p < pairs; ++p) {
            const unsigned first = unsigned(p) * 2 * width;
            const unsigned mid = std::min(first + width, threads);
            const unsigned end = std::min(first + 2 * width, threads);
            if (mid < end)
                std::inplace_merge(v.begin() + cut[first], v.begin() + cut[mid], v.begin() + cut[end], cmp);
        }
    }
}

/// radix_sort_inplace() - sort @v by @cmp, distributing first on the top kRadixBits of an unsigned @key_bits-bit primary key @key(x).
///
/// @cmp must order by key(x) first, and the result then equals std::sort(cmp).
template <typename V, typename Key, typename C>
void radix_sort_inplace(V& v, Key key, unsigned key_bits, C cmp, unsigned nthreads) {
    constexpr unsigned kRadixBits = 12, kRadixBuckets = 1u << kRadixBits;
    constexpr int kBucketsPerTask = 16;  // OpenMP dynamic chunk for the bucket sorts
    const std::size_t n = v.size();
    const unsigned threads = nthreads ? nthreads : 1;
    if (n < kParSortMinElements) { std::sort(v.begin(), v.end(), cmp); return; }
    const unsigned shift = key_bits > kRadixBits ? key_bits - kRadixBits : 0;
    std::vector<std::size_t> bucket_start(kRadixBuckets + 1, 0);
    {
        // per-thread histograms over equal slices, then one prefix sum
        std::vector<std::vector<std::size_t>> histogram(threads, std::vector<std::size_t>(kRadixBuckets, 0));
#ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
        for (int t = 0; t < int(threads); ++t)
            for (std::size_t i = n * std::size_t(t) / threads; i < n * std::size_t(t + 1) / threads; ++i)
                ++histogram[t][static_cast<std::size_t>(key(v[i]) >> shift)];
        for (unsigned b = 0; b < kRadixBuckets; ++b) {
            std::size_t count = 0;
            for (unsigned t = 0; t < threads; ++t) count += histogram[t][b];
            bucket_start[b + 1] = bucket_start[b] + count;
        }
    }
    {
        // American-flag permutation, swapping each element into its bucket's next free slot until every bucket holds only its own elements
        std::vector<std::size_t> next_free(bucket_start.begin(), bucket_start.end() - 1);
        for (unsigned b = 0; b < kRadixBuckets; ++b)
            while (next_free[b] < bucket_start[b + 1]) {
                const unsigned dest = static_cast<unsigned>(key(v[next_free[b]]) >> shift);
                if (dest == b) { ++next_free[b]; continue; }
                std::swap(v[next_free[b]], v[next_free[dest]++]);
            }
    }
#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(dynamic, kBucketsPerTask)
#endif
    for (int b = 0; b < int(kRadixBuckets); ++b)
        std::sort(v.begin() + bucket_start[b], v.begin() + bucket_start[b + 1], cmp);
}

}
