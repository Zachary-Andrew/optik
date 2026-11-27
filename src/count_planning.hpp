#pragma once

/**
 * @file
 * @brief Decisions the k-mer counter makes before pass 1 starts.
 *
 * choose_log_partitions() sizes the partition count from the input size.
 * HeavyMinimizers samples the input for m-mers that would swamp one partition and ranks them last.
 * warn_if_spill_disk_short() checks whether the spill directory could hold a worst-case spill.
 * The m-mer hash and base table are defined here too, because pass 1 in kmer_count.cpp uses the same ones.
 * Only kmer_count.cpp includes this file.
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "../include/bits.hpp"
#include "../include/fasta_reader.hpp"
#include "../include/kmer.hpp"
#include "../include/part_count_table.hpp"
#include "../include/spill.hpp"
#include "../include/wyhash.hpp"

namespace tpoptoa {
namespace detail {
namespace count_planning {

constexpr uint64_t kCtEmpty = ~uint64_t(0);    // never a canonical k-mer, since k <= 32
constexpr uint64_t kPartitionSeed = kGoldenRatio64;   // with ct_hash(), maps a minimizer to its partition
constexpr unsigned kMaxEnvLogPartitions = 24;         // ceiling for OPTIK_LOGP

// Minimizers and heavy-minimizer sampling.
constexpr std::size_t kSampleTableSlack = 64;
/// An m-mer counts as heavy once it minimizes more k-mers than this many average partitions.
constexpr uint64_t kHeavyPartitionFactor = 4;

// Input size estimates, in bases per file byte, for gzip FASTQ, gzip FASTA, plain FASTQ and plain FASTA.
constexpr double kGzFastqBasesPerByte = 1.6;
constexpr double kGzFastaBasesPerByte = 3.5;
constexpr double kFastqBasesPerByte = 0.45;
constexpr double kFastaBasesPerByte = 1.0;
/// Spilling writes about this many bytes per input base (measured on the CF2 gut reads).
constexpr double kSpillBytesPerBase = 1.6;

/// Free disk space is reported in decimal GB, as df does.
constexpr double kBytesPerGB = 1e9;

/// Orders m-mers and assigns partitions.
constexpr inline uint64_t ct_hash(uint64_t value) noexcept {
    return mix64(value, 0);
}

/// base_encode() as a 256-entry table.
/// The pass-1 scan reads every input base once through it, where a switch per base measured about 1% of the run on reads.
struct BaseLut {
    uint8_t codes[kByteValues];

    constexpr BaseLut() : codes{} {
        for (int c = 0; c < kByteValues; ++c) codes[c] = base_encode(static_cast<char>(c));
    }
};
inline constexpr BaseLut kBaseLut{};

/// log2 of the partition count for the inputs @p fasta_path, given as '\n'-joined paths.
///
/// About 64k k-mer occurrences per partition keep each private counting table near L2 size.
/// Every input file counts, gzip ones at about 4x their compressed size.
/// There are at most 2^16 partitions, since a 3.1 Gbp human genome at 2^14 had about 190k occurrences per partition, tables past L2 and a 45 s sizing pass.
/// `OPTIK_LOGP=<n>` sets the count to 2^n directly as a profiling aid.
inline unsigned choose_log_partitions(const char* fasta_path) {
    constexpr unsigned kMinLogPartitions = 6, kMaxLogPartitions = 16;
    constexpr std::size_t kInputBytesPerPartition = 65536;
    constexpr std::size_t kGzipExpansion = 4;
    unsigned logP = kMinLogPartitions;

    {
        std::size_t bytes = 0;

        for (const std::string& path : split_inputs(fasta_path)) {
            struct stat info;

            if (::stat(path.c_str(), &info) != 0 || info.st_size <= 0) continue;

            const bool is_gzip = path.size() > kGzSuffixLen && (path.compare(path.size() - kGzSuffixLen, kGzSuffixLen, kGzSuffix) == 0 ||
                                             path.compare(path.size() - kGzSuffixLen, kGzSuffixLen, kGzSuffixUpper) == 0);

            // Gzipped input counts at its rough uncompressed size.
            bytes += static_cast<std::size_t>(info.st_size) * (is_gzip ? kGzipExpansion : 1);
        }

        const std::size_t want = bytes / kInputBytesPerPartition;

        // About kInputBytesPerPartition of input per partition keeps each partition's table cache-sized.
        while (logP < kMaxLogPartitions && (std::size_t(1) << logP) < want) ++logP;
    }
    if (const char* env = std::getenv("OPTIK_LOGP")) logP = std::clamp(unsigned(std::atoi(env)), kMinLogPartitions, kMaxEnvLogPartitions);

    return logP;
}

/// The m-mer order of the counter, with heavy minimizers ranked last.
///
/// A few m-mers can be the minimizer of a large share of all k-mers, through adapter read-through such as Nextera or TruSeq or through poly-A.
/// The k-mers they minimize all go to one partition, which one thread then counts alone with a table and bin buffer of GBs.
/// sample() finds them in the start of the input, and order() ranks them after every other m-mer while keeping their own order, so a k-mer with any other m-mer in its window takes that one.
/// The order stays a fixed function of the m-mer, so all occurrences of a k-mer still land in one partition, and the result does not depend on the partitioning at all.
class HeavyMinimizers {
public:
    /// Rank of canonical m-mer @p mer, given by ct_hash(), a full finalizer, with heavy m-mers last.
    /// A single-multiply order measured slower, since a less random order gives shorter super-k-mers and more window rescans.
    uint64_t order(uint64_t mer) const {
        const uint64_t hash = ct_hash(mer) & ~kHeavyBit;

        return n_heavy_ && contains(mer) ? hash | kHeavyBit : hash;
    }

    /// Find the heavy minimizers of k-mer length @p kmer_len and minimizer length @p minimizer_len in up to 8 Mbp (kSample) from the start of the first input, against @p num_parts partitions.
    /// Rounds repeat because the k-mers moved off a heavy m-mer can make their next one heavy.
    void sample(const char* fasta_path, unsigned kmer_len, unsigned minimizer_len, std::size_t num_parts);

private:
    static constexpr uint64_t kHeavyBit = kTopBit;   // set in a heavy m-mer's order, ranking it last
    static constexpr unsigned kHeavyRounds = 3;                // sample passes (see sample())
    static constexpr std::size_t kMaxHeavy = 256;              // heavy m-mers kept at most
    static constexpr uint64_t kMinHeavyCount = 64;             // a heavy minimizer must minimize at least this many k-mers
    static constexpr std::size_t kHeavySetSlotsPerKey = 2;

    std::size_t slot(uint64_t mer) const {                      // Fibonacci hashing on the top heavy_lg_ bits
        return static_cast<std::size_t>((mer * kPartitionSeed) >> (kWordBits - heavy_lg_));
    }
    bool contains(uint64_t mer) const {
        for (std::size_t slot_idx = slot(mer);; slot_idx = (slot_idx + 1) & (heavy_.size() - 1)) {
            if (heavy_[slot_idx] == mer) return true;

            if (heavy_[slot_idx] == kCtEmpty) return false;
        }
    }

    static constexpr unsigned kHeavyInitialLog2 = 6;                      // 64 slots to start with
    std::vector<uint64_t> heavy_ = std::vector<uint64_t>(std::size_t(1) << kHeavyInitialLog2, kCtEmpty);   // open addressing, canonical m-mers
    unsigned heavy_lg_ = kHeavyInitialLog2;                               // log2(heavy_.size())
    std::size_t n_heavy_ = 0;
};

void HeavyMinimizers::sample(const char* fasta_path, unsigned kmer_len, unsigned minimizer_len, std::size_t num_parts) {
    const uint64_t mmask = kmer_mask(minimizer_len);
    const unsigned mrc_shift = kBitsPerBase * (minimizer_len - 1);
    const unsigned window = kmer_len - minimizer_len + 1;
    const std::vector<std::string> ins = split_inputs(fasta_path);
    std::vector<std::string> sample;
    std::size_t got = 0;
    // 8 Mbp from the start of the first input, where a long record such as a whole chromosome contributes its first bases only.
    constexpr std::size_t kSample = std::size_t(8) << 20;    // bases

    struct Enough {};

    if (ins[0] != "-") {
        try {
            iterate_sequences(ins[0].c_str(), [&](const SeqRecord& rec) {
                if (rec.seq.size() >= kmer_len) {
                    sample.push_back(rec.seq.substr(0, kSample - got));
                    got += sample.back().size();
                }

                if (got >= kSample) throw Enough{};
            }, false, kSample);
        } catch (const Enough&) {}
    }

    // The window's m-mer orders and m-mers sit in a ring, m-mer t at t & 63, and a window holds W <= 32 of them.
    // The minimizers' k-mer counts go to the partition counting table of part_count_table.hpp, sized for one minimizer per (W + 1) / 2 k-mers.
    constexpr std::size_t kRing = 64;

    static_assert(kRing >= kMaxK, "a window of a k <= 32 k-mer");
    std::vector<uint64_t> orders(kRing), cmv(kRing);
    PartCountTable cnt;

    for (unsigned round = 0; round < kHeavyRounds && !sample.empty(); ++round) {
        // Random minimizers have density 2 / (w + 1), so that many of the sampled k-mers start a new super-k-mer.
        constexpr uint64_t kMinimizerDensityNum = 2;

        cnt.begin(got * kMinimizerDensityNum / (window + 1) + kSampleTableSlack);
        uint64_t total = 0;
        // Runs of k-mers with one minimizer get one table update each.
        uint64_t run_cm = 0;
        uint32_t run_n = 0;
        auto flush = [&] { if (run_n) cnt.add(run_cm, 0, run_n); run_n = 0; };

        for (const std::string& seq : sample) {
            uint64_t fwd_mer = 0, rev_mer = 0;
            std::size_t mvalid = 0, mer_count = 0;
            // The window minimum is kept as the pass-1 scan keeps it, rescanned only when it leaves the window.
            uint64_t best = ~uint64_t(0);
            std::size_t best_i = 0;

            for (std::size_t i = 0; i < seq.size(); ++i) {
                const uint8_t base = kBaseLut.codes[uint8_t(seq[i])];

                if (base == kInvalidBase) { mvalid = 0; mer_count = 0; continue; }

                fwd_mer = ((fwd_mer << kBitsPerBase) | base) & mmask;
                rev_mer = ((rev_mer >> kBitsPerBase) | (static_cast<uint64_t>(base ^ kComplementMask) << mrc_shift)) & mmask;

                if (++mvalid < minimizer_len) continue;

                const uint64_t mer = fwd_mer < rev_mer ? fwd_mer : rev_mer;

                orders[mer_count & (kRing - 1)] = order(mer);
                cmv[mer_count & (kRing - 1)] = mer;
                ++mer_count;

                if (mvalid < kmer_len) continue;

                const std::size_t last = mer_count - 1, window_start = last + 1 - window;

                // Rescan only when the window is new or its minimum slid out. Otherwise one compare with the newest m-mer suffices.
                if (mer_count == window || best_i < window_start) {
                    best = ~uint64_t(0);

                    for (std::size_t t = window_start; t <= last; ++t)
                        if (orders[t & (kRing - 1)] < best) { best = orders[t & (kRing - 1)]; best_i = t; }
                } else if (orders[last & (kRing - 1)] < best) {
                    best = orders[last & (kRing - 1)]; best_i = last;
                }

                const uint64_t minimizer = cmv[best_i & (kRing - 1)];

                // A change of minimizer ends a run of k-mers sharing one.
                if (minimizer != run_cm) flush();

                run_cm = minimizer;
                ++run_n;
                ++total;
            }
        }

        flush();
        // A heavy m-mer is the minimizer of more k-mers than four average partitions.
        const uint64_t thr = std::max<uint64_t>(kMinHeavyCount, total * kHeavyPartitionFactor / num_parts);
        std::vector<std::pair<uint64_t, uint64_t>> add;

        cnt.for_each([&](uint64_t mer, uint32_t count, uint8_t) {
            if (count > thr && !contains(mer)) add.push_back({count, mer});
        });
        cnt.finish(0);

        if (add.empty()) break;

        std::sort(add.rbegin(), add.rend());

        for (const auto& [count, mer] : add) {
            if (n_heavy_ >= kMaxHeavy) break;

            if ((n_heavy_ + 1) * kHeavySetSlotsPerKey > heavy_.size()) {        // keep the set at most half full
                constexpr std::size_t kHeavyGrowth = 2;
                std::vector<uint64_t> old(heavy_.size() * kHeavyGrowth, kCtEmpty);

                old.swap(heavy_);
                ++heavy_lg_;

                for (uint64_t mer : old) if (mer != kCtEmpty) {
                    std::size_t slot_idx = slot(mer);

                    while (heavy_[slot_idx] != kCtEmpty) slot_idx = (slot_idx + 1) & (heavy_.size() - 1);

                    heavy_[slot_idx] = mer;
                }
            }

            std::size_t slot_idx = slot(mer);

            while (heavy_[slot_idx] != kCtEmpty) slot_idx = (slot_idx + 1) & (heavy_.size() - 1);

            heavy_[slot_idx] = mer;
            ++n_heavy_;
        }
    }

    if (n_heavy_ && std::getenv("OPTIK_TIMING"))
        std::fprintf(stderr, "[timing] %zu heavy minimizers ranked last (sample of %zu bases)\n", n_heavy_, got);
}

/// Warn when the spill directory @p sdir may not hold what the counter could spill past @p buf_budget bytes of buffers.
///
/// A spilling count writes about 1.6 bytes per input base, counting 2-bit super-k-mers with their k - 1 overlaps and record headers, as measured on the CF2 gut reads.
/// Bases are estimated per input kind.
/// A shortfall gives a warning with the numbers rather than an error, since the estimate is rough and ENOSPC still stops the run with an error.
inline void warn_if_spill_disk_short(const char* fasta_path, uint64_t buf_budget, const std::string& sdir) {
    double bases = 0;

    for (const std::string& path : split_inputs(fasta_path)) {
        struct stat info;

        if (::stat(path.c_str(), &info) != 0 || info.st_size <= 0) continue;

        std::string lower = path;

        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char ch) { return std::tolower(ch); });
        const bool is_gzip = lower.size() > kGzSuffixLen && lower.compare(lower.size() - kGzSuffixLen, kGzSuffixLen, kGzSuffix) == 0;
        const bool is_fastq = lower.find(".fq") != std::string::npos || lower.find(".fastq") != std::string::npos;

        // File size converts to bases by format. FASTQ spends half its bytes on qualities, and gzip shrinks both.
        bases += double(info.st_size) * (is_gzip ? (is_fastq ? kGzFastqBasesPerByte : kGzFastaBasesPerByte) : (is_fastq ? kFastqBasesPerByte : kFastaBasesPerByte));
    }

    const double need = kSpillBytesPerBase * bases - double(buf_budget);
    const uint64_t avail = spill::free_bytes(sdir);

    // This is only a warning, because the bases-per-byte estimate is rough and a full disk still makes the spill write fail with ENOSPC.
    if (need > 0 && avail && double(avail) < need)
        std::fprintf(stderr, "[count] warning: the k-mer buffers may spill about %.0f GB to %s, which has "
                     "%.0f GB free; set --tmp-dir to a larger disk or raise --max-mem\n",
                     need / kBytesPerGB, sdir.c_str(), double(avail) / kBytesPerGB);
}

}
}
}
