#pragma once

/**
 * @file
 * @brief Deterministic parallel sorts for the index builders.
 *
 * Given a strict total order, both produce exactly what std::sort would.
 * Builds therefore stay byte-identical regardless of thread count.
 *
 * radix_sort_inplace() starts with one American-flag pass (McIlroy, Bostic and McIlroy, "Engineering radix sort", 1993).
 * That is an in-place MSD radix pass, where elements are swapped straight into their buckets without a second array.
 * Afterwards, each bucket is finished by a comparison sort running as its own OpenMP task.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

namespace tpoptoa {

/// Below this size, starting threads costs more than the sort itself.
inline constexpr std::size_t kParSortMinElements = std::size_t(1) << 16;

/// Sorts equal slices in parallel, then merges neighbouring runs pairwise over log2(nthreads) rounds.
template <typename V, typename C>
void par_sort(V& vec, C cmp, unsigned nthreads) {
    const std::size_t size = vec.size();
    const unsigned threads = nthreads ? nthreads : 1;

    constexpr unsigned kMergeWays = 2;   // each merge round joins neighbouring pairs of sorted runs

    if (threads < kMergeWays || size < kParSortMinElements) { std::sort(vec.begin(), vec.end(), cmp); return; }

    std::vector<std::size_t> cut(threads + 1);

    // Threads first sort equal slices, then neighbouring sorted runs are merged pairwise with the run width doubling each round.
    for (unsigned t = 0; t <= threads; ++t) cut[t] = size * t / threads;

#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (int t = 0; t < int(threads); ++t) std::sort(vec.begin() + cut[t], vec.begin() + cut[t + 1], cmp);

    for (unsigned width = 1; width < threads; width *= kMergeWays) {
        const int pairs = int((threads + kMergeWays * width - 1) / (kMergeWays * width));

#ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
        for (int p = 0; p < pairs; ++p) {
            const unsigned first = unsigned(p) * kMergeWays * width;
            const unsigned mid = std::min(first + width, threads);
            const unsigned end = std::min(first + kMergeWays * width, threads);

            if (mid < end)
                std::inplace_merge(vec.begin() + cut[first], vec.begin() + cut[mid], vec.begin() + cut[end], cmp);
        }
    }
}

/// Buckets on the top kRadixBits of the @p key_bits-bit primary key @p key(x) before comparison-sorting.
/// For that reason @p cmp has to order by key(x) first, otherwise the buckets would come out in the wrong order.
template <typename V, typename Key, typename C>
void radix_sort_inplace(V& vec, Key key, unsigned key_bits, C cmp, unsigned nthreads) {
    constexpr unsigned kRadixBits = 12, kRadixBuckets = 1u << kRadixBits;
    constexpr int kBucketsPerTask = 16;  // OpenMP dynamic chunk for the bucket sorts
    const std::size_t size = vec.size();
    const unsigned threads = nthreads ? nthreads : 1;

    if (size < kParSortMinElements) { std::sort(vec.begin(), vec.end(), cmp); return; }

    const unsigned shift = key_bits > kRadixBits ? key_bits - kRadixBits : 0;
    std::vector<std::size_t> bucket_start(kRadixBuckets + 1, 0);

    {
        // Each thread counts its own slice, and a prefix sum then turns the counts into bucket starts.
        std::vector<std::vector<std::size_t>> histogram(threads, std::vector<std::size_t>(kRadixBuckets, 0));

#ifdef _OPENMP
        #pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
        for (int t = 0; t < int(threads); ++t)
            for (std::size_t i = size * std::size_t(t) / threads; i < size * std::size_t(t + 1) / threads; ++i)
                ++histogram[t][static_cast<std::size_t>(key(vec[i]) >> shift)];

        for (unsigned b = 0; b < kRadixBuckets; ++b) {
            std::size_t count = 0;

            for (unsigned t = 0; t < threads; ++t) count += histogram[t][b];

            bucket_start[b + 1] = bucket_start[b] + count;
        }
    }
    {
        // Swap every element into the next free slot of its bucket until each bucket is pure.
        std::vector<std::size_t> next_free(bucket_start.begin(), bucket_start.end() - 1);

        for (unsigned b = 0; b < kRadixBuckets; ++b)
            while (next_free[b] < bucket_start[b + 1]) {
                const unsigned dest = static_cast<unsigned>(key(vec[next_free[b]]) >> shift);

                if (dest == b) { ++next_free[b]; continue; }

                std::swap(vec[next_free[b]], vec[next_free[dest]++]);
            }
    }
#ifdef _OPENMP
    #pragma omp parallel for num_threads(threads) schedule(dynamic, kBucketsPerTask)
#endif
    for (int b = 0; b < int(kRadixBuckets); ++b)
        std::sort(vec.begin() + bucket_start[b], vec.begin() + bucket_start[b + 1], cmp);
}

}
