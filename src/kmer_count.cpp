/**
 * @file
 * @brief Canonical k-mer counting with inline edge recording, as declared in kmer_count.hpp.
 *
 * Sources
 *   - Minimizer-partitioned super-k-mer counting from KMC 2 (Deorowicz et al. 2015).
 *     Unlike KMC we keep partitions in memory, and spill.hpp writes them to disk only past the memory budget.
 *   - The edge-centric graph of Cuttlefish 2 (Khan et al. 2022), where a k-mer's edges are the (k+1)-mers observed around it.
 *   - HyperLogLog (Flajolet et al. 2007) for sizing large partitions' tables.
 *   - Fibonacci hashing (Knuth, TAOCP vol. 3) and group prefetching (Slotin, Algorithms for Modern Hardware).
 *
 * Phases of count_into_graph()
 *   1. Heavy minimizers.
 *      Sampling the input finds m-mers that would minimize a large share of k-mers, as adapter read-through or poly-A do.
 *      Those are then ranked after every other m-mer.
 *   2. Pass 1.
 *      Readers cut the input into chunks of kChunk k-mer starts, and workers split each chunk into super-k-mers (maximal runs sharing a minimizer).
 *      Each super-k-mer is appended 2-bit packed, flanking bases included, to per-thread buffers of its partition.
 *      Because the flanks are the observed (k+1)-mers of the first and last k-mer, no separate (k+1)-mer pass is needed.
 *      Past the memory budget, buffers spill to disk.
 *   3. Pass 2.
 *      One thread counts each partition in a private, reused, cache-sized table (part_count_table.hpp).
 *      A sizing run first counts the solid k-mers, and an insert run then places them in the graph.
 *      That graph is the DbgGraph table, the compact DbgGraphEF, or dbg_part.hpp's partitioned compaction.
 *   4. With min_count > 1, edge masks are finally trimmed to solid neighbours.
 */
#include "../include/bits.hpp"
#include "../include/thread_budget.hpp"
#include "kmer_count.hpp"

#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <thread>
#include <unordered_map>

#include <malloc.h>
#include <sys/stat.h>

#include "../include/fasta_reader.hpp"
#include "../include/kmer.hpp"
#include "../include/dbg_graph.hpp"
#include "../include/dbg_graph_ef.hpp"
#include "../include/dbg_part.hpp"
#include "../include/part_count_table.hpp"
#include "../include/spill.hpp"
#include "../include/wyhash.hpp"
#include "count_planning.hpp"

namespace tpoptoa {
namespace detail {

using namespace count_planning;

namespace {

// Records and partitions.
constexpr std::size_t kChunk = 1u << 18;       // k-mer starts per work item
/// The record header's length field has 10 bits, so longer super-k-mers are cut.
constexpr std::size_t kMaxRecordBases = 1000;
/// Each record starts with a u16 header packed as nbases | flank flags << kHdrFlagShift | flank codes << kHdrCodeShift.
constexpr unsigned kHdrLenMask = 0x3FFu, kHdrFlagShift = 10, kHdrCodeShift = 12;
constexpr unsigned kLeftFlank = 1u, kRightFlank = 2u, kHdrFlagMask = kLeftFlank | kRightFlank;
constexpr std::size_t kHdrBytes = 2;
/// Spilled solid lists use LEB128 varints, seven payload bits per byte with the high bit set while more bytes follow.
constexpr unsigned kVarintPayloadBits = 7;
constexpr uint8_t kVarintMore = 0x80, kVarintPayloadMask = 0x7F;
constexpr std::size_t kSpillHdrBytes = 16;            // on-disk SpillHdr
constexpr std::size_t kKRecBytes = 13;                // on-disk packed KRec
/// Bins group runs of partitions, at most 2^kMaxLogBins of them.
constexpr unsigned kMaxLogBins = 10;

/// The minimizer length is k minus this, capped at the largest supported length.
constexpr unsigned kMinimizerShortfall = 3;

// Threads and work queues.
constexpr std::size_t kQueuedChunksPerThread = 4;
constexpr unsigned kThreadsPerReader = 4;
constexpr std::size_t kBatchSlack = 4096;    // spare bytes in a read batch beyond one work item

// Memory plan.
/// With auto memory, buffers get 1/kAutoBufferShareDen of the limit.
/// With --max-mem they get kExplicitBufferShareNum/kExplicitBufferShareDen of it instead.
constexpr uint64_t kAutoBufferShareDen = 2;
constexpr uint64_t kExplicitBufferShareNum = 8;
constexpr uint64_t kExplicitBufferShareDen = 10;
/// Partitioned compaction keeps only half of --max-mem for buffers, and at most kPartAutoBufferCap by default.
constexpr uint64_t kPartBufferShareDen = 2;
constexpr uint64_t kPartAutoBufferCap = uint64_t(1) << 30;   // 1 GiB
constexpr uint64_t kMinThreadBufferBytes = uint64_t(4) << 20;   // 4 MiB, below which spilling starts
/// Headroom of 5% for the per-thread counting tables and the reader.
constexpr uint64_t kHeadroomDiv = 20;
/// A graph needs 20% working room beyond its own bytes.
constexpr uint64_t kWorkingRoomDiv = 5;
constexpr int kMaxSpillRounds = 8;
/// Fragment and end-key streams get 1/kStreamShareDiv of the free memory per thread, clamped to [kMinStreamBuffer, kMaxStreamBuffer].
constexpr uint64_t kStreamShareDiv = 16;
constexpr uint64_t kMinStreamBuffer = uint64_t(4) << 20;    // 4 MiB
constexpr uint64_t kMaxStreamBuffer = uint64_t(16) << 20;   // 16 MiB
/// Solid-record key buckets get half of the free memory, with at least kMinBucketBuffer per thread.
constexpr uint64_t kBucketRoomShareDen = 2;
constexpr uint64_t kMinBucketBuffer = uint64_t(16) << 20;   // 16 MiB
/// Buffers larger than these are freed after each partition rather than reused.
constexpr std::size_t kMaxKeptDecodeBytes = std::size_t(64) << 20;
constexpr std::size_t kMaxKeptNodes = std::size_t(1) << 20;

// Table sizing.
constexpr unsigned kMaxEfBucketBits = 12;    // leading key bits the compact graph buckets on
constexpr unsigned kLevelSlots = 8;          // OPTIK_TIMING level histogram, last slot collects deeper levels
/// Before any partition is measured, assume this distinct share, and afterwards use the median of the last kRatioHistory.
constexpr double kInitialDistinctRatio = 0.125;
constexpr std::size_t kRatioHistory = 32;
/// HyperLogLog (Flajolet et al. 2007) uses alpha_m = kHllAlphaNum / (1 + kHllAlphaDen / m).
/// Below kHllLinearCountingCut * m it falls back to linear counting.
constexpr double kHllAlphaNum = 0.7213;
constexpr double kHllAlphaDen = 1.079;
constexpr double kHllLinearCountingCut = 2.5;
constexpr double kHllMargin = 1.05;
constexpr double kRatioMargin = 1.5;          // safety factor on a size estimated from recent shares
constexpr std::size_t kTableSlack = 16;
/// Small partitions say little about the distinct share, so they don't update the history.
constexpr std::size_t kMinRatioOccurrences = 1024;

// Messages.
constexpr std::size_t kMsgBufLen = 256;
/// Memory is reported in binary MB and GB.
constexpr double kBytesPerMiB = 1024.0 * 1024.0;
constexpr double kBytesPerGiB = 1024.0 * 1024.0 * 1024.0;

struct Chunk {
    std::shared_ptr<const std::string> seq;
    std::size_t begin_pos, end_pos;                        // k-mer start positions [lo, hi)
};

/// Bounded multi-producer multi-consumer work queue under a mutex, with one push and pop per kChunk k-mers, far off the per-k-mer hot path.
class ChunkQueue {
public:
    explicit ChunkQueue(std::size_t cap) : cap_(cap) {}

    void push(Chunk chunk) {
        std::unique_lock<std::mutex> lock(mutex_);

        // Blocking on a full queue caps the chunks held in memory.
        not_full_.wait(lock, [&] { return queue_.size() < cap_; });
        queue_.push_back(std::move(chunk));
        not_empty_.notify_one();
    }
    bool pop(Chunk& chunk) {
        std::unique_lock<std::mutex> lock(mutex_);

        not_empty_.wait(lock, [&] { return !queue_.empty() || done_; });

        // An empty queue after close() tells the worker to stop.
        if (queue_.empty()) return false;

        chunk = std::move(queue_.front());
        queue_.pop_front();
        not_full_.notify_one();
        return true;
    }
    void close() {
        std::lock_guard<std::mutex> lock(mutex_);

        done_ = true;
        not_empty_.notify_all();
    }
private:
    std::mutex mutex_;
    std::condition_variable not_empty_, not_full_;
    std::deque<Chunk> queue_;
    std::size_t cap_;
    bool done_ = false;
};

}

// All occurrences of a k-mer share its minimizer, so they meet in exactly one partition.
// Pass 1 stores about (k + #kmers - 1) / 4 + 2 bytes per super-k-mer instead of a slot per distinct k-mer in one shared table, and pass 2's inserts hit L2 instead of DRAM.
namespace {

/// Append-only byte buffer of packed super-k-mer records.
///
/// Each record opens with a little-endian u16 header, followed by ceil(nbases / 4) bytes of 2-bit bases.
/// The header packs nbases in 10 bits, since records are capped at kMaxRecordBases, then flank flags << 10 (bit 0 left, bit 1 right) and flank codes << 12 (left | right << 2).
/// It is kept at two bytes rather than four because on read sets the spill consists mostly of these records.
struct SkBuf {
    std::vector<uint8_t> bytes;
    std::size_t kmers = 0;
};

/// A solid k-mer of a counted partition.
struct LocalSlot { uint64_t key; uint32_t cnt; uint32_t edge_mask; };

/// A run of packed super-k-mer records (in memory or read back from disk).
struct Span { const uint8_t* data; std::size_t nbytes; };
/// Header of one partition buffer in a spill file.
struct SpillHdr { uint32_t part, pad; uint64_t nbytes; };
static_assert(sizeof(SpillHdr) == kSpillHdrBytes, "spill header layout");
/// A spill-file extent holding records of one bin (one thread, one flush).
struct Extent { uint64_t off, len; };
/// A solid k-mer on its way into the compact graph (key-ordered buckets), packed into 13 bytes in memory and in the spill files.
struct __attribute__((packed)) KRec { uint64_t key; uint32_t cnt; uint8_t edge_mask; };
static_assert(sizeof(KRec) == kKRecBytes, "packed record");

}

/// `OPTIK_SPILL_TEST=<bytes>` replaces a per-thread buffer share, so tests can force many small spills on small inputs.
inline uint64_t spill_share_or_test(uint64_t share) {
    const char* env = std::getenv("OPTIK_SPILL_TEST");

    return env ? std::strtoull(env, nullptr, kDecimalBase) : share;
}

/// Warns when -c > 1 meets single-copy input such as a genome, an assembly or unitigs, where nearly every k-mer occurs once and the cutoff silently keeps only repeats.
/// Reads have coverage well above 1, so fewer than 1.5 occurrences per distinct k-mer with under half of them solid says the input is not reads.
inline void warn_if_single_copy_input(uint32_t min_count, uint64_t n_occ, uint64_t n_distinct, std::size_t n_solid) {
    constexpr uint64_t kReadCoverageExtraDen = 2;   // reads exceed 1 + 1/2 occurrences per distinct k-mer
    constexpr uint64_t kMinSolidShareDen = 2;       // and keep at least half of their k-mers solid
    constexpr double kPercent = 100.0;

    if (min_count <= 1 || n_distinct == 0) return;

    if (n_occ >= n_distinct + n_distinct / kReadCoverageExtraDen || n_solid >= n_distinct / kMinSolidShareDen) return;

    std::fprintf(stderr, "[count] warning: mean k-mer coverage is %.2fx, so -c %u kept %zu of %llu "
                 "distinct k-mers (%.1f%%); for a genome or assembly use -c 1\n",
                 double(n_occ) / double(n_distinct), min_count, n_solid,
                 static_cast<unsigned long long>(n_distinct), kPercent * double(n_solid) / double(n_distinct));
}

/// How much memory the counter may use, and how much of it the partition buffers get.
struct CountMemoryPlan {
    bool explicit_limit;          // true when --max-mem set the budget
    uint64_t limit;               // buffers plus graph must stay under this
    uint64_t buffer_budget;       // all threads' partition buffers together
    uint64_t per_thread_buffer;   // one worker's share before it spills
};

/// Works out the memory plan for @p nthreads workers.
///
/// The partition buffers grow with the input (coverage times genome), while the graph grows only with its distinct solid k-mers.
/// Below the budget nothing is spilled, so the counter behaves exactly as if spilling did not exist.
inline CountMemoryPlan plan_count_memory(unsigned nthreads, bool partitioned) {
    CountMemoryPlan plan{};
    const uint64_t budget = spill::memory_budget(plan.explicit_limit);

    // With auto memory, buffers get half of what is available and the whole of it bounds buffers plus graph.
    // With --max-mem, both numbers come from the limit.
    plan.limit = plan.explicit_limit ? budget : budget * kAutoBufferShareDen;
    plan.buffer_budget = plan.explicit_limit ? budget / kExplicitBufferShareDen * kExplicitBufferShareNum : budget;

    // Partitioned compaction holds little after counting, and spilling on fast storage costs little (a human genome took 86 s spilled against 90 s in memory).
    // Its buffers therefore get a smaller share, half of --max-mem and by default at most 1 GB.
    if (partitioned) plan.buffer_budget = plan.explicit_limit ? budget / kPartBufferShareDen : std::min<uint64_t>(budget, kPartAutoBufferCap);

    plan.per_thread_buffer = std::max<uint64_t>(plan.buffer_budget / nthreads, kMinThreadBufferBytes);

    plan.per_thread_buffer = spill_share_or_test(plan.per_thread_buffer);

    return plan;
}

void count_into_graph(const char* fasta_path, std::size_t kmer_len_in,
                      uint32_t min_count, unsigned nthreads, DbgGraph& graph,
                      std::unique_ptr<DbgGraphEF>* compact, std::unique_ptr<PartUnitigs>* part) {
    const unsigned kmer_len = static_cast<unsigned>(kmer_len_in);

    nthreads = resolve_threads(nthreads);

    if (nthreads == 0) nthreads = 1;

    const uint64_t kmask = kmer_mask(kmer_len);
    const unsigned rc_shift = kBitsPerBase * (kmer_len - 1);

    if (min_count == 0) min_count = 1;

    PhaseTimer timer;

    // The minimizer length is long enough to spread partitions and short enough for long super-k-mers (window w = K - M + 1).
    constexpr unsigned kMaxMinimizerLen = 13;
    const unsigned minimizer_len = kmer_len <= kMinimizerShortfall ? kmer_len : std::min(kMaxMinimizerLen, kmer_len - kMinimizerShortfall);
    const uint64_t mmask = kmer_mask(minimizer_len);
    const unsigned mrc_shift = kBitsPerBase * (minimizer_len - 1);
    const unsigned window = kmer_len - minimizer_len + 1;

    const unsigned logP = choose_log_partitions(fasta_path);
    const std::size_t num_parts = std::size_t(1) << logP;

    std::vector<std::vector<SkBuf>> bufs(nthreads, std::vector<SkBuf>(num_parts));
    // Heavy minimizers are ranked after every other m-mer (HeavyMinimizers).
    HeavyMinimizers heavy;

    heavy.sample(fasta_path, kmer_len, minimizer_len, num_parts);
    auto order_of = [&](uint64_t mer) { return heavy.order(mer); };

    timer.lap("count: heavy minimizers");

    ChunkQueue queue(kQueuedChunksPerThread * nthreads);

    // Each worker keeps its partition buffers in memory up to its share of the budget and spills past it (spill.hpp).
    const CountMemoryPlan plan = plan_count_memory(nthreads, part != nullptr);
    const bool explicit_limit = plan.explicit_limit;
    const uint64_t limit = plan.limit, buf_budget = plan.buffer_budget, per_thread = plan.per_thread_buffer;

    // Bins are runs of P / B consecutive partitions, the unit read back at once.
    const unsigned logB = std::min(logP, kMaxLogBins);
    const std::size_t num_bins = std::size_t(1) << logB;
    const unsigned bin_shift = logP - logB;
    const std::string sdir = spill::spill_dir();

    warn_if_spill_disk_short(fasta_path, buf_budget, sdir);
    std::vector<spill::File> files(nthreads);
    std::vector<std::vector<std::vector<Extent>>> ext(nthreads, std::vector<std::vector<Extent>>(num_bins));
    std::vector<uint64_t> held(nthreads, 0);          // buffer capacity per thread
    std::atomic<bool> failed{false};
    std::mutex err_mu;
    std::string err_msg;
    std::once_flag spill_note;
    auto fail_with = [&](const std::exception& err) {
        // The message is written before the flag is raised, since readers check the flag and then read the message.
        std::lock_guard<std::mutex> lock(err_mu);

        if (!failed.load()) { err_msg = err.what(); failed.store(true, std::memory_order_release); }
    };
    // Write thread t's buffers of partitions [p0, p1) to its spill file and free them.
    // Records keep their partition, and extents are indexed by bin.
    auto spill_range = [&](unsigned thread_idx, std::size_t part_begin, std::size_t part_end) {
        std::call_once(spill_note, [&] {
            std::fprintf(stderr, "[count] k-mer buffers exceed the memory budget (%.2f GB%s): "
                         "spilling them to %s\n", buf_budget / kBytesPerGiB,
                         explicit_limit ? ", --max-mem" : ", auto", sdir.c_str());
            if (spill::ram_backed(sdir))
                std::fprintf(stderr, "[count] warning: %s is RAM-backed (tmpfs); spilling there "
                             "saves no memory -- set --tmp-dir to a disk\n", sdir.c_str());
        });
        files[thread_idx].open(sdir);
        auto& mine = bufs[thread_idx];
        std::size_t nonempty = 0;

        for (std::size_t p = part_begin; p < part_end; ++p) nonempty += !mine[p].bytes.empty();

        if (nonempty == 0) return;

        std::vector<SpillHdr> hdr; hdr.reserve(nonempty);
        constexpr std::size_t kIovecsPerPartition = 2;   // its SpillHdr and its bytes
        std::vector<struct iovec> iov; iov.reserve(kIovecsPerPartition * nonempty);
        uint64_t off = files[thread_idx].size();
        std::size_t cur_bin = SIZE_MAX;

        for (std::size_t p = part_begin; p < part_end; ++p) {
            auto& bytes_buf = mine[p].bytes;

            if (bytes_buf.empty()) continue;

            const std::size_t bin = p >> bin_shift;

            if (bin != cur_bin) { ext[thread_idx][bin].push_back({off, 0}); cur_bin = bin; }

            hdr.push_back({static_cast<uint32_t>(p), 0, bytes_buf.size()});
            iov.push_back({&hdr.back(), sizeof(SpillHdr)});
            iov.push_back({bytes_buf.data(), bytes_buf.size()});
            ext[thread_idx][bin].back().len += sizeof(SpillHdr) + bytes_buf.size();
            off += sizeof(SpillHdr) + bytes_buf.size();
        }

        files[thread_idx].append(iov);

        for (std::size_t p = part_begin; p < part_end; ++p) {
            held[thread_idx] -= mine[p].bytes.capacity();
            std::vector<uint8_t>().swap(mine[p].bytes);   // kmers stays, since it counts spilled ones too
        }
    };

    // Pass 1 puts super-k-mers into partition buffers.
    // `single` processes that one chunk on the calling thread and returns (with -t 1 the reader does the work itself, so one thread runs, not two).
    auto worker = [&](unsigned tid, const Chunk* single = nullptr) {
        auto& mine = bufs[tid];
        uint64_t& hold = held[tid];
        Chunk chunk;
        std::vector<uint64_t> mer_hashes;         // canonical m-mer hashes of the run
        std::vector<uint8_t> base_codes;          // 2-bit codes of the chunk's scanned bases
        uint8_t packed[packed_bytes_for(kMaxRecordBases) + kHdrBytes];   // one packed record
        auto body = [&](const Chunk& chunk) {
            if (failed.load(std::memory_order_relaxed)) return;   // drain the queue

            const std::string& seq = *chunk.seq;
            const std::size_t seq_len = seq.size();
            const std::size_t scan0 = chunk.begin_pos > 0 ? chunk.begin_pos - 1 : 0;
            const std::size_t scan_end = std::min(seq_len, chunk.end_pos + kmer_len);   // exclusive

            base_codes.resize(scan_end - scan0);
            // State of the current super-k-mer.
            std::size_t sk_first = 0, sk_last = 0;   // start positions of first/last k-mer
            uint64_t sk_min = 0;
            bool open = false;
            auto emit = [&]() {
                if (!open) return;

                const std::size_t rec_start = sk_first, rec_end = sk_last + kmer_len;         // bases [a, e)
                const std::size_t num_bases = rec_end - rec_start;
                SkBuf& part_buf = mine[ct_hash(sk_min ^ kPartitionSeed) & (num_parts - 1)];
                uint8_t flags = 0, codes = 0;

                if (rec_start > 0) { const uint8_t base = kBaseLut.codes[uint8_t(seq[rec_start - 1])]; if (base != kInvalidBase) { flags |= kLeftFlank; codes |= base; } }

                if (rec_end < seq_len) { const uint8_t base = kBaseLut.codes[uint8_t(seq[rec_end])]; if (base != kInvalidBase) { flags |= kRightFlank; codes |= base << kBitsPerBase; } }

                // The record is packed on the stack from the codes the scan already computed, then appended in one copy without a zero fill.
                const std::size_t nby = packed_bytes_for(num_bases);
                const uint16_t header = static_cast<uint16_t>(num_bases | (unsigned(flags) << kHdrFlagShift) | (unsigned(codes) << kHdrCodeShift));

                packed[0] = static_cast<uint8_t>(header);                  // low byte first
                packed[1] = static_cast<uint8_t>(header >> kBitsPerByte);
                const uint8_t* codes_ptr = base_codes.data() + (rec_start - scan0);
                std::size_t pos = 0, out_pos = kHdrBytes;

                for (; pos + kBasesPerByte <= num_bases; pos += kBasesPerByte) {
                    uint8_t byte = 0;

                    for (unsigned r = 0; r < kBasesPerByte; ++r) byte |= static_cast<uint8_t>(codes_ptr[pos + r] << (r * kBitsPerBase));

                    packed[out_pos++] = byte;
                }

                if (pos < num_bases) {
                    uint8_t last_byte = 0;

                    for (unsigned r = 0; pos < num_bases; ++pos, ++r) last_byte |= static_cast<uint8_t>(codes_ptr[pos] << (r * kBitsPerBase));

                    packed[out_pos++] = last_byte;
                }

                const std::size_t cap0 = part_buf.bytes.capacity();

                part_buf.bytes.insert(part_buf.bytes.end(), packed, packed + kHdrBytes + nby);

                if (part_buf.bytes.capacity() != cap0) hold += part_buf.bytes.capacity() - cap0;

                part_buf.kmers += sk_last - sk_first + 1;
                open = false;
            };
            uint64_t fwd_mer = 0, rev_mer = 0;            // rolling m-mer on both strands
            std::size_t mvalid = 0;             // valid bases in the current run
            std::size_t run0 = scan0;           // start of the current ACGT run

            mer_hashes.clear();
            // Window minimum over the W m-mers of the current k-mer.
            uint64_t best = ~uint64_t(0);
            std::size_t best_i = 0;             // index in mh of the minimum

            for (std::size_t i = scan0; i < scan_end; ++i) {
                const uint8_t base = kBaseLut.codes[uint8_t(seq[i])];

                base_codes[i - scan0] = base;

                if (base == kInvalidBase) {
                    emit();
                    mvalid = 0; mer_hashes.clear(); run0 = i + 1;
                    continue;
                }

                fwd_mer = ((fwd_mer << kBitsPerBase) | base) & mmask;
                rev_mer = ((rev_mer >> kBitsPerBase) | (static_cast<uint64_t>(base ^ kComplementMask) << mrc_shift)) & mmask;
                ++mvalid;

                if (mvalid < minimizer_len) continue;

                mer_hashes.push_back(order_of(fwd_mer < rev_mer ? fwd_mer : rev_mer));   // m-mer ending at i

                if (mvalid < kmer_len) continue;

                const std::size_t start = i + 1 - kmer_len;

                if (start >= chunk.end_pos) break;

                // The m-mers of this k-mer are mh[q .. q + W - 1], with q = start - run0.
                const std::size_t window_start = start - run0, last = window_start + window - 1;
                const uint64_t new_hash = mer_hashes[last];

                if (start == run0 || best_i < window_start) {         // first k-mer, or the minimum left, so rescan
                    best = ~uint64_t(0);

                    for (std::size_t t = window_start; t <= last; ++t) {
                        const bool is_less = mer_hashes[t] < best;     // leftmost minimum, so ties keep the first

                        best = is_less ? mer_hashes[t] : best;
                        best_i = is_less ? t : best_i;
                    }
                } else if (new_hash < best) {
                    best = new_hash; best_i = last;
                }

                if (start < chunk.begin_pos) continue;               // context k-mer only

                if (open && best == sk_min) {
                    sk_last = start;
                } else {
                    emit();
                    open = true; sk_first = sk_last = start; sk_min = best;
                }

                if (sk_last - sk_first + kmer_len >= kMaxRecordBases) emit();   // cap record length
            }

            emit();

            if (hold > per_thread) {
                try { spill_range(tid, 0, num_parts); } catch (const std::exception& err) { fail_with(err); }
            }
        };

        if (single) { body(*single); return; }

        while (queue.pop(chunk)) body(chunk);
    };

    std::vector<std::thread> workers;

    workers.reserve(nthreads);
    // Several input files are read at once, each by its own reader with its own batch, since one gzip stream decodes only so fast and reading a multi-file read set file after file left the workers waiting on one decoder.
    // Record order cannot change the counts, because workers already append to partitions in any order, so the output is the same.
    // Readers come out of -t, as Cuttlefish 3's decoders come out of its --threads.
    // Of the t threads, r = t/4 are readers, at most one per file, each inflating its own .gz in-thread through fasta_reader.hpp, and t - r are pass-1 super-k-mer workers.
    const std::vector<std::string> inputs = split_inputs(fasta_path);
    const unsigned nreaders = static_cast<unsigned>(
        std::min<std::size_t>(inputs.size(), std::max(1u, nthreads / kThreadsPerReader)));
    // With -t 1 there are no worker threads, and the one reader processes its chunks itself.
    const bool inline_work = nthreads == 1;
    const unsigned nworkers = inline_work ? 0 : nthreads > nreaders ? nthreads - nreaders : 1;

    for (unsigned t = 0; t < nworkers; ++t) workers.emplace_back(worker, t, nullptr);

    // Short records (reads) are batched into one work item joined by 'N', which resets the scanner and blanks flanks, so the result is identical to separate records.
    auto push_seq = [&](std::shared_ptr<const std::string> seq_ptr) {
        const std::size_t nstarts = seq_ptr->size() - kmer_len + 1;

        // Long records split into kChunk-start chunks that share one sequence buffer.
        for (std::size_t chunk_start = 0; chunk_start < nstarts; chunk_start += kChunk) {
            const Chunk chunk{seq_ptr, chunk_start, std::min(nstarts, chunk_start + kChunk)};

            // With -t 1 the reader counts each chunk itself.
            if (inline_work) worker(0, &chunk);
            else queue.push(chunk);
        }
    };
    std::atomic<std::size_t> next_input{0};
    std::exception_ptr read_err;
    std::mutex read_mu;
    std::atomic<bool> read_failed{false};
    auto reader = [&]() {
        std::string batch;

        batch.reserve(kChunk + kBatchSlack);
        auto flush_batch = [&]() {
            if (batch.size() >= kmer_len) push_seq(std::make_shared<const std::string>(std::move(batch)));

            batch = std::string();
            batch.reserve(kChunk + kBatchSlack);
        };

        try {
            // Readers claim whole files because a gzip stream decodes on one thread.
            for (std::size_t input_idx; (input_idx = next_input.fetch_add(1)) < inputs.size();) {
                iterate_sequences(inputs[input_idx].c_str(), [&](const SeqRecord& rec) {
                    if (failed.load(std::memory_order_acquire))      // a spill write failed, so stop reading
                        throw std::runtime_error(err_msg);

                    if (read_failed.load(std::memory_order_relaxed))   // another reader failed
                        throw std::runtime_error("input read stopped");

                    if (rec.seq.size() < kmer_len) return;

                    if (rec.seq.size() >= kChunk) {
                        push_seq(std::make_shared<const std::string>(rec.seq));
                        return;
                    }

                    if (!batch.empty()) batch.push_back('N');

                    batch += rec.seq;

                    if (batch.size() >= kChunk) flush_batch();
                });
            }

            flush_batch();
        } catch (...) {
            std::lock_guard<std::mutex> lock(read_mu);

            if (!read_err) read_err = std::current_exception();

            read_failed.store(true);
        }
    };

    // A reader error (missing file, I/O error) must not unwind past running worker threads (std::terminate), so the queue is closed and the workers joined before rethrowing.
    {
        std::vector<std::thread> readers;

        for (unsigned r = 1; r < nreaders; ++r) readers.emplace_back(reader);

        reader();

        for (auto& reader_thread : readers) reader_thread.join();
    }
    if (read_err) {
        queue.close();

        for (auto& worker_thread : workers) worker_thread.join();

        std::rethrow_exception(read_err);
    }

    queue.close();

    for (auto& worker_thread : workers) worker_thread.join();

    if (failed) throw std::runtime_error(err_msg);

    timer.lap("count: parse + super-k-mers");

    // Pass 2 counts each partition in a private cache-sized table.
    // It runs twice, first only to learn how many solid k-mers each partition has (which sizes the graph table exactly), then to insert the solid k-mers straight into their graph slots.
    // Counting a partition is cheap, since it is L2-resident.
    // The alternative, keeping all partitions' solid lists, gathering, sorting and building a separate solid-key set, held several copies of each k-mer.
    // The second pass frees partition buffers as it goes.
    //
    // A partition's records are spans of its in-memory buffers plus, after a spill, the parts read back from disk.
    // Without a spill the passes hand out single partitions.
    // With one they hand out bins, each read back with a few large preads (the next bins' extents are prefetched), so a worker holds one bin of disk records at a time.
    bool spilled = false;

    for (auto& file : files) spilled = spilled || file.is_open();

    // Spilled buffers were freed into malloc's free lists (many small vectors, below the mmap threshold), so those pages are handed back to the OS and the resident set is what is still held.
    if (spilled) malloc_trim(0);

    auto bin_bytes = [&](std::size_t bin) {
        uint64_t total = 0;

        for (unsigned t = 0; t < nthreads; ++t) for (const Extent& extent : ext[t][bin]) total += extent.len;

        return total;
    };
    uint64_t n_distinct = 0, n_occ = 0;        // from the sizing pass, declared early for count_pass
    // With the compact DbgGraphEF graph, solid k-mers go to key-ordered buckets.
    bool ef_mode = false;
    const unsigned ef_bits = std::min(kMaxEfBucketBits, kBitsPerBase * kmer_len);
    const unsigned ef_shift = kBitsPerBase * kmer_len - ef_bits;
    std::unique_ptr<spill::BucketSpill<KRec>> ksp;
    // Partitioned compaction (dbg_part.hpp) uses no global table at all.
    const bool part_mode = part != nullptr;
    std::unique_ptr<PartBuilder> part_builder;
    // Work is handed out largest first by k-mer occurrences (bins after a spill, otherwise partitions).
    // A few partitions are far heavier than the rest (adapter read-through, poly-A), and one taken last would leave a single thread working while the others wait.
    std::vector<std::size_t> order(spilled ? num_bins : num_parts);

    {
        std::vector<uint64_t> weights(order.size(), 0);

        for (unsigned w = 0; w < nthreads; ++w)
            for (std::size_t p = 0; p < num_parts; ++p) weights[spilled ? p >> bin_shift : p] += bufs[w][p].kmers;

        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;

        std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) { return weights[lhs] > weights[rhs]; });
    }
    // With -c > 1 and partitioned compaction, the first pass already has each partition's solid k-mers with their in-partition edge masks, which is all the compaction needs (it applies the cross-partition claims itself).
    // The first pass therefore saves them, in memory or after a spill appended to a per-thread file, and the compaction reads them back instead of reading every super-k-mer again and recounting every occurrence.
    // On reads the solid k-mers are a small fraction of the occurrences.
    const bool save_solid = part_mode && min_count > 1;

    // Spilled lists are sorted by key and written compactly as varint key deltas, varint counts and one mask byte (about 5 bytes a k-mer instead of 16).
    struct SolidRef { uint64_t off = 0, count = 0, bytes = 0; uint32_t thread_idx = 0; };
    std::vector<SolidRef> sref(save_solid ? num_parts : 0);
    std::vector<std::vector<PartNode>> smem(save_solid && !spilled ? num_parts : 0);
    std::vector<spill::File> sfile(save_solid && spilled ? nthreads : 0);
    std::vector<std::vector<uint8_t>> swbuf(save_solid && spilled ? nthreads : 0);
    constexpr std::size_t kSwbuf = std::size_t(8) << 20;
    auto flush_solid = [&](unsigned thread_idx) {
        if (swbuf.empty() || swbuf[thread_idx].empty()) return;

        sfile[thread_idx].open(sdir);
        std::vector<struct iovec> iov{{swbuf[thread_idx].data(), swbuf[thread_idx].size()}};

        sfile[thread_idx].append(iov);
        swbuf[thread_idx].clear();
    };
    // This is LEB128, seven bits per byte with the high bit set while more bytes follow.
    auto put_varint = [](std::vector<uint8_t>& out, uint64_t value) {
        while (value >= kVarintMore) { out.push_back(static_cast<uint8_t>(value | kVarintMore)); value >>= kVarintPayloadBits; }

        out.push_back(static_cast<uint8_t>(value));
    };
    auto save_nodes = [&](unsigned thread_idx, std::size_t part, std::vector<PartNode>& nodes) {
        if (nodes.empty()) return;

        // Without a spill, the list is kept in memory as it is.
        if (!spilled) { smem[part] = nodes; return; }

        std::sort(nodes.begin(), nodes.end(), [](const PartNode& lhs, const PartNode& rhs) { return lhs.key < rhs.key; });
        std::vector<uint8_t>& out = swbuf[thread_idx];
        const uint64_t start = out.size();
        uint64_t prev = 0;

        // Sorted keys make the deltas small, so most of them fit in one or two varint bytes.
        for (const PartNode& node : nodes) {
            put_varint(out, node.key - prev);
            prev = node.key;
            put_varint(out, node.cnt);
            out.push_back(static_cast<uint8_t>(node.edge_mask));
        }

        sref[part] = SolidRef{sfile[thread_idx].size() + start, nodes.size(), out.size() - start, thread_idx};

        if (out.size() >= kSwbuf) flush_solid(thread_idx);
    };
    auto count_pass = [&](bool insert) {
        std::atomic<std::size_t> next{0};
        std::atomic<std::size_t> total{0};
        std::atomic<uint64_t> all_distinct{0}, all_occ{0};
        std::atomic<uint64_t> max_tab{0}, max_bin{0}, max_occ{0};   // OPTIK_TIMING report
        std::atomic<uint64_t> all_deep{0};
        std::atomic<uint64_t> all_lv[kLevelSlots] = {};

        tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
          try {
            // One table per thread is reused across partitions.
            // Only the slots a partition touched are visited and reset afterwards, so a partition costs time in its occurrences rather than in the table capacity.
            PartCountTable table;                 // reused across partitions
            double ratio = kInitialDistinctRatio;             // distinct / occurrences of the last partition
            std::array<double, kRatioHistory> rat{};     // recent ratios of ordinary partitions
            std::size_t nrat = 0;
            constexpr unsigned kHllP = 12;    // 2^12 HyperLogLog registers
            std::vector<uint8_t> hll;
            std::vector<LocalSlot> solid;
            std::size_t mine = 0;
            uint64_t my_distinct = 0, my_occ = 0, my_max_bin = 0, my_max_occ = 0, my_max_tab = 0, my_deep = 0;
            uint64_t my_lv[kLevelSlots] = {0};
            std::vector<Span> spans;
            std::vector<PartNode> pnodes;
            bool release_solid = false;
            constexpr std::size_t kBigTabBytes = std::size_t(128) << 20;
            constexpr std::size_t kBigBin = std::size_t(64) << 20;  // bytes
            auto count_partition = [&](std::size_t part, const Span* spans_ptr, std::size_t num_spans) {
                std::size_t occ = 0;

                for (unsigned w = 0; w < nthreads; ++w) occ += bufs[w][part].kmers;

                if (occ == 0) return;

                my_occ += occ;
                my_max_occ = std::max<uint64_t>(my_max_occ, occ);
                // f(canonical, edge bits) for every k-mer occurrence of the partition's records.
                auto each_kmer = [&](auto&& func) {
                    for (std::size_t si = 0; si < num_spans; ++si) {
                        const uint8_t* const span_data = spans_ptr[si].data;
                        const std::size_t by_n = spans_ptr[si].nbytes;
                        std::size_t offset = 0;

                        while (offset < by_n) {
                            const unsigned header = span_data[offset] | (unsigned(span_data[offset + 1]) << kBitsPerByte);
                            const std::size_t num_bases = header & kHdrLenMask;
                            const uint8_t flags = static_cast<uint8_t>((header >> kHdrFlagShift) & kHdrFlagMask);
                            const uint8_t codes = static_cast<uint8_t>(header >> kHdrCodeShift);
                            const uint8_t* packed = span_data + offset + kHdrBytes;
                            auto base = [&](std::size_t pos) -> uint8_t { return (packed[packed_byte(pos)] >> packed_shift(pos)) & kBaseMask; };
                            uint64_t fwd = 0, rev = 0;

                            for (std::size_t t = 0; t < num_bases; ++t) {
                                const uint8_t cur_base = base(t);

                                fwd = ((fwd << kBitsPerBase) | cur_base) & kmask;
                                rev = ((rev >> kBitsPerBase) | (static_cast<uint64_t>(cur_base ^ kComplementMask) << rc_shift)) & kmask;

                                if (t + 1 < kmer_len) continue;

                                const std::size_t kmer_start = t + 1 - kmer_len;
                                // At either end of a record, the neighbour base comes from the flank code in its header, so edges across super-k-mer boundaries are kept.
                                const uint8_t prev_base = kmer_start > 0 ? base(kmer_start - 1) : ((flags & kLeftFlank) ? (codes & kBaseMask) : kInvalidBase);
                                const uint8_t next_base = t + 1 < num_bases ? base(t + 1) : ((flags & kRightFlank) ? ((codes >> kBitsPerBase) & kBaseMask) : kInvalidBase);
                                const uint64_t cano = fwd < rev ? fwd : rev;
                                uint8_t in_b = 0, out_b = 0;

                                // Edges are stored relative to the canonical orientation, so a k-mer read on the reverse strand swaps and complements its neighbours.
                                // A palindrome takes both branches.
                                if (fwd == cano) {
                                    if (next_base != kInvalidBase) out_b |= 1u << next_base;

                                    if (prev_base != kInvalidBase) in_b  |= 1u << prev_base;
                                }

                                if (rev == cano) {
                                    if (next_base != kInvalidBase) in_b  |= 1u << (next_base ^ kComplementMask);

                                    if (prev_base != kInvalidBase) out_b |= 1u << (prev_base ^ kComplementMask);
                                }

                                func(cano, static_cast<uint8_t>((in_b << kOutEdgeBits) | out_b));
                            }

                            offset += kHdrBytes + packed_bytes_for(num_bases);
                        }
                    }
                };
                // Level 0 is sized from the expected distinct k-mers, so the table stays one level (part_count_table.hpp).
                // A large partition, among the first handed out and often one heavy minimizer, gets a HyperLogLog estimate from a pre-pass over its k-mers (Flajolet et al. 2007, with 4 KB of registers in L1 and a standard error of 1.04 / sqrt(2^12), about 1.6%).
                // The rest use the distinct share this thread last measured, with a margin.
                constexpr std::size_t kHllMinOccurrences = std::size_t(1) << 21;
                constexpr uint64_t kHllSeed = 0x5851F42D4C957F2DULL;
                double est;

                if (occ > kHllMinOccurrences) {
                    hll.assign(std::size_t(1) << kHllP, 0);
                    each_kmer([&](uint64_t cano, uint8_t) {
                        const uint64_t hash = ct_hash(cano ^ kHllSeed);
                        const std::size_t reg = hash >> (kWordBits - kHllP);
                        const uint8_t leading = static_cast<uint8_t>(__builtin_clzll((hash << kHllP) | (uint64_t(1) << (kHllP - 1))) + 1);

                        if (leading > hll[reg]) hll[reg] = leading;
                    });
                    double sum = 0; std::size_t zeros = 0;

                    for (uint8_t reg_val : hll) { sum += std::ldexp(1.0, -int(reg_val)); zeros += reg_val == 0; }

                    constexpr double mreg = double(std::size_t(1) << kHllP);

                    // The raw estimate is alpha_m * m^2 / sum 2^-M[j] with alpha_m = 0.7213 / (1 + 1.079 / m), switching to linear counting below 2.5 m (Flajolet et al. 2007).
                    est = kHllAlphaNum / (1 + kHllAlphaDen / mreg) * mreg * mreg / sum;

                    if (est < kHllLinearCountingCut * mreg && zeros) est = mreg * std::log(mreg / double(zeros));   // linear counting

                    est *= kHllMargin;                                // 5% margin
                } else {
                    // The median of this thread's recent ratios is used, so a single highly repetitive partition (ratio about 1e-4) cannot shrink the next ordinary one's table.
                    double med = kInitialDistinctRatio;

                    if (nrat) {
                        std::array<double, kRatioHistory> tmp;
                        const std::size_t count = std::min<std::size_t>(nrat, tmp.size());

                        std::copy_n(rat.begin(), count, tmp.begin());
                        constexpr std::size_t kMedianDiv = 2;   // the median sits at the middle index
                        const std::size_t middle = count / kMedianDiv;

                        std::nth_element(tmp.begin(), tmp.begin() + middle, tmp.begin() + count);
                        med = tmp[middle];
                    }

                    est = static_cast<double>(occ) * med * kRatioMargin;
                }

                table.begin(static_cast<std::size_t>(est) + kTableSlack);
                {
                    // Group prefetching keeps a ring of kRing k-mers whose lines are in flight.
                    // Each new k-mer's lines are prefetched and the oldest is inserted, and the order does not change the counts.
                    constexpr std::size_t kRing = 16;            // a power of two
                    uint64_t ring_keys[kRing], ring_hashes[kRing];
                    uint8_t ring_bits[kRing];
                    std::size_t written = 0, consumed = 0;                 // FIFO counters

                    each_kmer([&](uint64_t cano, uint8_t bits) {
                        const uint64_t hash = PartCountTable::hash(cano);

                        table.prefetch(hash);

                        // Inserting kRing k-mers behind the prefetch hides each line's miss.
                        if (written - consumed == kRing) {
                            const std::size_t slot = consumed++ & (kRing - 1);

                            table.add_hashed(ring_keys[slot], ring_hashes[slot], ring_bits[slot]);
                        }

                        const std::size_t slot = written++ & (kRing - 1);

                        ring_keys[slot] = cano; ring_hashes[slot] = hash; ring_bits[slot] = bits;
                    });
                    // The ring's last entries are flushed once the partition ends.
                    while (consumed < written) {
                        const std::size_t slot = consumed++ & (kRing - 1);

                        table.add_hashed(ring_keys[slot], ring_hashes[slot], ring_bits[slot]);
                    }
                }
                const std::size_t num_distinct = table.size();

                my_max_tab = std::max<uint64_t>(my_max_tab, table.bytes());
                my_deep += table.deep();
                my_lv[std::min(table.levels(), kLevelSlots - 1)]++;
                ratio = std::min(1.0, static_cast<double>(num_distinct) / static_cast<double>(occ));

                // Mid-sized partitions add to the distinct-share history used to size later tables. Large ones are sized by their own HyperLogLog pass.
                if (occ >= kMinRatioOccurrences && occ <= kHllMinOccurrences) rat[nrat++ % rat.size()] = ratio;

                my_distinct += num_distinct;
                solid.clear();

                if (part_mode) {
                    // The partition's table gives exact in-partition masks, where a neighbour counted here below min_count is no edge.
                    // Cross-partition edges are claimed in the first pass and handed over in the second.
                    const unsigned first_base_shift = kBitsPerBase * (kmer_len - 1);
                    auto nbr = [&](uint64_t kmer, unsigned bit) {
                        const uint64_t neighbour = bit < kOutEdgeBits ? (((kmer << kBitsPerBase) | bit) & kmask) : ((uint64_t(bit - kOutEdgeBits) << first_base_shift) | (kmer >> kBitsPerBase));
                        const uint64_t rev = reverse_complement(neighbour, kmer_len);

                        return neighbour < rev ? neighbour : rev;
                    };
                    const bool keep = insert || save_solid;

                    if (keep) pnodes.clear();

                    table.for_each([&](uint64_t key, uint32_t cnt, uint8_t edge_mask) {
                        if (cnt < min_count) return;

                        ++mine;
                        uint32_t mask = edge_mask;

                        if (min_count > 1)
                            for (unsigned b = 0; b < kEdgeMaskBits; ++b) {
                                if (!(mask & (1u << b))) continue;

                                const uint32_t nbr_count = table.find(nbr(key, b));

                                // A neighbour inside this partition is looked up directly, and one outside it is claimed for the cross-partition match.
                                if (nbr_count != PartCountTable::kNotFound) { if (nbr_count < min_count) mask &= ~(1u << b); }
                                else if (!insert) part_builder->claim(thread_idx, static_cast<uint32_t>(part), part_builder->mask_edge(key, b));
                            }

                        if (keep) pnodes.push_back(PartNode{key, cnt, mask});
                    });
                    // The first pass only records the partition's solid list, and the second pass compacts it into fragments.
                    if (insert) part_builder->partition(thread_idx, static_cast<uint32_t>(part), pnodes);
                    else if (save_solid) save_nodes(thread_idx, part, pnodes);
                } else {
                    table.for_each([&](uint64_t key, uint32_t cnt, uint8_t edge_mask) {
                        if (cnt < min_count) return;

                        ++mine;

                        if (insert) solid.push_back(LocalSlot{key, cnt, edge_mask});
                    });
                }

                // A heavy partition (one minimizer shared by many k-mers, as with adapter read-through or poly-A) grows this thread's table far past the usual size.
                // If it were kept for reuse every thread would end the pass holding its largest, so past the cap it is released, after its last use below, since solid is still needed there.
                if (table.bytes() > kBigTabBytes) release_solid = true;

                table.finish(kBigTabBytes);
                // One thread owns each k-mer through its partition, so claiming takes one CAS on the key word and then a plain store of the node word.
                // A locked CAS keeps later loads from running ahead, so each claim would wait out its own DRAM miss.
                // The home line kPrefetchAhead claims ahead is prefetched for writing instead, and the misses overlap.
                constexpr std::size_t kPrefetchAhead = 16;

                // The compact graph is built later from key-ordered buckets, picked by the top key bits.
                if (ef_mode) {
                    for (const LocalSlot& entry : solid)
                        ksp->add(thread_idx, entry.key >> ef_shift, KRec{entry.key, entry.cnt, static_cast<uint8_t>(entry.edge_mask)});
                } else {
                    for (std::size_t j = 0; j < solid.size(); ++j) {
                        if (j + kPrefetchAhead < solid.size()) graph.prefetch_w(solid[j + kPrefetchAhead].key);

                        graph.put(graph.claim(solid[j].key), solid[j].cnt,
                                  static_cast<uint8_t>(solid[j].edge_mask));
                    }
                }

                // An oversized partition's scratch memory is released.
                if (release_solid) {
                    std::vector<LocalSlot>().swap(solid);
                    std::vector<PartNode>().swap(pnodes);
                    release_solid = false;
                }
            };
            auto mem_spans = [&](std::size_t part) {
                for (unsigned w = 0; w < nthreads; ++w)
                    if (!bufs[w][part].bytes.empty())
                        spans.push_back({bufs[w][part].bytes.data(), bufs[w][part].bytes.size()});
            };
            auto free_mem = [&](std::size_t part) {
                if (!insert && !save_solid) return;   // a later pass reads them again

                for (unsigned w = 0; w < nthreads; ++w)           // free as we go
                    std::vector<uint8_t>().swap(bufs[w][part].bytes);
            };

            if (!spilled) {
                for (std::size_t idx, part; (idx = next.fetch_add(1)) < num_parts && (part = order[idx], true);) {
                    spans.clear();
                    mem_spans(part);
                    count_partition(part, spans.data(), spans.size());
                    free_mem(part);
                }
            } else {
                const std::size_t parts_per_bin = num_parts >> logB;   // partitions per bin
                std::vector<uint8_t> binbuf;
                std::vector<std::vector<Span>> disk(parts_per_bin);

                for (std::size_t idx, bin; (idx = next.fetch_add(1)) < num_bins && (bin = order[idx], true);) {
                    if (failed.load(std::memory_order_relaxed)) break;

                    // Start the reads of the bins the other workers take next.
                    if (const std::size_t next_idx = idx + nthreads; next_idx < num_bins)
                        for (unsigned w = 0; w < nthreads; ++w)
                            for (const Extent& extent : ext[w][order[next_idx]]) files[w].will_need(extent.off, extent.len);

                    binbuf.resize(bin_bytes(bin));
                    my_max_bin = std::max<uint64_t>(my_max_bin, binbuf.size());
                    uint64_t filled = 0;

                    // Large sequential reads gather every worker's extents for this bin.
                    for (unsigned w = 0; w < nthreads; ++w)
                        for (const Extent& extent : ext[w][bin]) {
                            files[w].read(extent.off, extent.len, binbuf.data() + filled);
                            filled += extent.len;
                        }

                    for (auto& bin_spans : disk) bin_spans.clear();

                    for (uint64_t offset = 0; offset < filled;) {
                        SpillHdr header;

                        std::memcpy(&header, binbuf.data() + offset, sizeof header);

                        // A record outside this bin or past the data means a damaged spill file.
                        if ((header.part >> bin_shift) != bin || offset + sizeof header + header.nbytes > filled)
                            throw std::runtime_error("spill file corrupt (bad record header)");

                        disk[header.part - (bin << bin_shift)].push_back({binbuf.data() + offset + sizeof header, header.nbytes});
                        offset += sizeof header + header.nbytes;
                    }

                    for (std::size_t j = 0; j < parts_per_bin; ++j) {
                        const std::size_t part = (bin << bin_shift) + j;

                        // A partition is its spilled spans plus what stayed in memory.
                        spans.assign(disk[j].begin(), disk[j].end());
                        mem_spans(part);
                        count_partition(part, spans.data(), spans.size());
                        free_mem(part);
                    }

                    // This is the last pass over the bin, so its disk blocks go back.
                    if (insert || save_solid)
                        for (unsigned w = 0; w < nthreads; ++w)
                            for (const Extent& extent : ext[w][bin]) files[w].punch(extent.off, extent.len);

                    if (binbuf.capacity() > kBigBin) std::vector<uint8_t>().swap(binbuf);
                }
            }

            if (!insert && save_solid) flush_solid(thread_idx);

            all_deep.fetch_add(my_deep);

            for (int i = 0; i < int(kLevelSlots); ++i) all_lv[i].fetch_add(my_lv[i]);

            {
                // Each thread folds its maxima into the shared ones with a compare-and-swap loop, because std::atomic has no fetch_max.
                uint64_t value = my_max_tab, current = max_tab.load();

                while (value > current && !max_tab.compare_exchange_weak(current, value)) {}

                value = my_max_bin; current = max_bin.load();

                while (value > current && !max_bin.compare_exchange_weak(current, value)) {}

                value = my_max_occ; current = max_occ.load();

                while (value > current && !max_occ.compare_exchange_weak(current, value)) {}
            }
            total.fetch_add(mine);
            all_distinct.fetch_add(my_distinct);
            all_occ.fetch_add(my_occ);
          } catch (const std::exception& err) { fail_with(err); }
        });
        // Worker errors are rethrown here on the calling thread after every worker has joined.
        if (failed) throw std::runtime_error(err_msg);

        // Sizing-pass totals size the graph table and drive the low-coverage warning.
        if (!insert) { n_distinct = all_distinct.load(); n_occ = all_occ.load(); }

        if (std::getenv("OPTIK_TIMING"))
            std::fprintf(stderr, "[timing] count pass: P = 2^%u, largest partition %llu occurrences, "
                         "per-thread table up to %.0f MB, bin buffer up to %.0f MB\n", logP,
                         static_cast<unsigned long long>(max_occ.load()), max_tab.load() / kBytesPerMiB,
                         max_bin.load() / kBytesPerMiB);

        if (std::getenv("OPTIK_TIMING")) {
            std::fprintf(stderr, "[timing] count pass: %llu of %llu distinct keys below level 0; partitions by levels:",
                         static_cast<unsigned long long>(all_deep.load()), static_cast<unsigned long long>(all_distinct.load()));
            for (int i = 1; i < int(kLevelSlots); ++i) std::fprintf(stderr, " %d:%llu", i, static_cast<unsigned long long>(all_lv[i].load()));

            std::fprintf(stderr, "\n");
        }

        return total.load();
    };

    if (part_mode) {
        // At -c 1 every counted k-mer is solid and no sizing pass is needed.
        // At -c > 1 the first pass claims the cross-partition edges.
        const uint64_t res = spill::resident_bytes();
        const uint64_t room = limit > res ? limit - res : 0;
        // Fragment and end-key streams keep 4-16 MB per thread in memory and the rest on disk, read back once and sequentially by bucket.
        uint64_t per = std::clamp<uint64_t>(room / kStreamShareDiv / nthreads, kMinStreamBuffer, kMaxStreamBuffer);

        per = spill_share_or_test(per);
        part_builder = std::make_unique<PartBuilder>(kmer_len, nthreads, num_parts, min_count > 1, per, sdir);

        if (min_count > 1) {
            const std::size_t num_solid = count_pass(false);

            warn_if_single_copy_input(min_count, n_occ, n_distinct, num_solid);

            timer.lap("count: partitions (first pass: edge claims)");
            part_builder->match_claims();
            timer.lap("count: claims matched");
        }

        if (!save_solid) {
            count_pass(true);
        } else {
            // The saved solid lists, largest first, each compacted by one thread.
            bufs.clear();
            malloc_trim(0);
            std::vector<std::size_t> porder(num_parts);

            for (std::size_t p = 0; p < num_parts; ++p) porder[p] = p;

            auto nsolid = [&](std::size_t part) { return spilled ? sref[part].count : smem[part].size(); };

            std::stable_sort(porder.begin(), porder.end(),
                             [&](std::size_t lhs, std::size_t rhs) { return nsolid(lhs) > nsolid(rhs); });
            std::atomic<std::size_t> next{0};

            tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                try {
                    std::vector<PartNode> nodes;
                    std::vector<uint8_t> enc;

                    for (std::size_t idx; (idx = next.fetch_add(1)) < num_parts;) {
                        if (failed.load(std::memory_order_relaxed)) break;

                        const std::size_t part = porder[idx];

                        if (nsolid(part) == 0) break;            // sorted, so the rest are empty

                        if (spilled) {
                            const SolidRef& ref = sref[part];

                            enc.resize(ref.bytes);
                            sfile[ref.thread_idx].read(ref.off, ref.bytes, enc.data());
                            nodes.resize(ref.count);
                            const uint8_t* cursor = enc.data();
                            const uint8_t* const end_ptr = cursor + ref.bytes;
                            auto get = [&]() {
                                uint64_t value = 0;

                                for (unsigned shift = 0; cursor < end_ptr; shift += kVarintPayloadBits) {
                                    const uint8_t byte = *cursor++;

                                    value |= uint64_t(byte & kVarintPayloadMask) << shift;

                                    if (!(byte & kVarintMore)) return value;
                                }

                                throw std::runtime_error("spill file corrupt (solid list)");
                            };
                            uint64_t key = 0;

                            // A running sum undoes the delta encoding.
                            for (PartNode& node : nodes) {
                                key += get();
                                node.key = key;
                                node.cnt = static_cast<uint32_t>(get());

                                if (cursor >= end_ptr) throw std::runtime_error("spill file corrupt (solid list)");

                                node.edge_mask = *cursor++;
                            }

                            sfile[ref.thread_idx].punch(ref.off, ref.bytes);          // read once, so give it back

                            if (enc.capacity() > kMaxKeptDecodeBytes) std::vector<uint8_t>().swap(enc);
                        } else {
                            // Unspilled lists change hands by swap.
                            nodes.swap(smem[part]);
                            std::vector<PartNode>().swap(smem[part]);
                        }

                        part_builder->partition(thread_idx, static_cast<uint32_t>(part), nodes);

                        if (nodes.capacity() > kMaxKeptNodes) std::vector<PartNode>().swap(nodes);
                    }
                } catch (const std::exception& err) { fail_with(err); }
            });
            if (failed) throw std::runtime_error(err_msg);
        }

        bufs.clear();
        // The partition buffers were many small vectors, and their freed pages sit in malloc's free lists until handed back.
        malloc_trim(0);
        timer.lap("count: partitions (compaction)");
        *part = part_builder->finish(timer);
        return;
    }

    const std::size_t n_solid = count_pass(false);

    warn_if_single_copy_input(min_count, n_occ, n_distinct, n_solid);

    timer.lap("count: partitions (sizing)");

    // The graph is allocated before the buffers are released, since freed buffer memory is not reused for the graph's pages.
    // If the resident set plus the graph plus the bins read back would pass the limit, buffers are spilled, from the last partitions back, until it fits.
    // This is measured rather than estimated, because the resident set includes malloc overhead and what fragmentation keeps.
    // A graph larger than the limit by itself cannot be helped by spilling.
    {
        const uint64_t gbytes = DbgGraph::table_bytes(n_solid);
        auto max_bin = [&] {
            uint64_t largest = 0;

            for (std::size_t b = 0; b < num_bins; ++b) largest = std::max(largest, bin_bytes(b));

            return largest;
        };
        // 5% headroom for the per-thread counting tables and the reader.
        auto need = [&] {
            const uint64_t need_bytes = gbytes + spill::resident_bytes() + (spilled ? max_bin() * nthreads : 0);

            return need_bytes + need_bytes / kHeadroomDiv;
        };
        uint64_t needed = need();

        if (timer.enabled)
            std::fprintf(stderr, "[timing] count: memory plan: graph %.0f MB + resident %.0f MB "
                         "(+ bins read back) = %.0f MB vs limit %.0f MB (spilled %s)\n",
                         gbytes / kBytesPerMiB, spill::resident_bytes() / kBytesPerMiB, needed / kBytesPerMiB,
                         limit / kBytesPerMiB, spilled ? "yes" : "no");

        // The table does not fit, so the compact graph is used (Elias-Fano keys, about 8 bytes per k-mer, and OPTIK_DBG_COMPACT=1 forces it for tests).
        const uint64_t efbytes = DbgGraphEF::bytes_for(n_solid, kmer_len);

        if (compact && (std::getenv("OPTIK_DBG_COMPACT") || gbytes + gbytes / kWorkingRoomDiv > limit)) {
            ef_mode = true;

            if (efbytes + efbytes / kWorkingRoomDiv > limit) {
                char msg[kMsgBufLen];

                std::snprintf(msg, sizeof msg, "the compact graph of %zu solid k-mers needs %.2f GB "
                              "(%.2f GB + 20%% working room), more than the memory limit (%.2f GB%s)", n_solid,
                              (efbytes + efbytes / kWorkingRoomDiv) / kBytesPerGiB, efbytes / kBytesPerGiB,
                              limit / kBytesPerGiB, explicit_limit ? ", --max-mem" : ", available");
                // A user-given --max-mem is enforced as an error. The automatic limit from MemAvailable only prints a warning.
                if (explicit_limit) throw std::runtime_error(std::string(msg) + "; raise --max-mem or -c");

                std::fprintf(stderr, "[count] warning: %s\n", msg);
            }

            std::fprintf(stderr, "[count] %zu solid k-mers: the table would need %.2f GB, so the graph is "
                         "built compact (%.2f GB, Elias-Fano keys)\n", n_solid, gbytes / kBytesPerGiB,
                         efbytes / kBytesPerGiB);
            // The solid records wait in key buckets, in memory up to what the compact graph leaves of the limit and on disk past it.
            const uint64_t res = spill::resident_bytes();
            const uint64_t room = limit > res + efbytes ? limit - res - efbytes : 0;
            uint64_t per = std::max<uint64_t>(room / kBucketRoomShareDen / nthreads, kMinBucketBuffer);

            per = spill_share_or_test(per);
            ksp = std::make_unique<spill::BucketSpill<KRec>>(nthreads, std::size_t(1) << ef_bits, per, sdir);
        } else if (gbytes > limit) {
            char msg[kMsgBufLen];

            std::snprintf(msg, sizeof msg, "the graph of %zu solid k-mers needs %.2f GB, more than the "
                          "memory limit (%.2f GB%s)", n_solid, gbytes / kBytesPerGiB,
                          limit / kBytesPerGiB, explicit_limit ? ", --max-mem" : ", available");
            // A user-given --max-mem is enforced as an error. The automatic limit from MemAvailable only prints a warning.
            if (explicit_limit) throw std::runtime_error(std::string(msg) + "; raise --max-mem or -c");

            std::fprintf(stderr, "[count] warning: %s\n", msg);
        } else {
            std::size_t cut = num_parts;

            for (int round = 0; needed > limit && cut > 0 && round < kMaxSpillRounds; ++round) {
                // Free at least the excess capacity, starting from the end.
                const std::size_t range_end = cut;
                uint64_t freed = 0;

                while (cut > 0 && freed < needed - limit) {
                    --cut;

                    for (unsigned t = 0; t < nthreads; ++t) freed += bufs[t][cut].bytes.capacity();
                }

                tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                    try { spill_range(thread_idx, cut, range_end); } catch (const std::exception& err) { fail_with(err); }
                });
                if (failed) throw std::runtime_error(err_msg);

                // RSS is sampled again because glibc may keep freed buffers mapped.
                spilled = true;
                malloc_trim(0);
                needed = need();
            }

            if (cut < num_parts) timer.lap("count: spill before insert");
        }
    }
    if (ef_mode) {
        count_pass(true);
        bufs.clear();

        if (ksp->spilled()) malloc_trim(0);

        timer.lap("count: partitions (insert, key buckets)");
        // Buckets come in key order, and bucket b's nodes start after every smaller bucket's, so the buckets are loaded, sorted and written into the Elias-Fano arrays in parallel, the keys being distinct.
        auto ef_graph = std::make_unique<DbgGraphEF>(kmer_len);

        ef_graph->begin(n_solid);
        const std::size_t num_buckets = std::size_t(1) << ef_bits;
        std::vector<std::size_t> base(num_buckets + 1, 0);

        // A prefix sum over bucket sizes gives each bucket its first node index, and threads then fill buckets independently.
        for (std::size_t b = 0; b < num_buckets; ++b) base[b + 1] = base[b] + ksp->size(b);

        if (base[num_buckets] != n_solid) throw std::logic_error("compact graph: bucket sizes do not add up");

        {
            std::atomic<std::size_t> next_bucket{0};

            tpoptoa::run_on_threads(nthreads, [&](unsigned) {
                try {
                    std::vector<KRec> recs;

                    for (std::size_t bucket; (bucket = next_bucket.fetch_add(1)) < num_buckets;) {
                        if (bucket + nthreads < num_buckets) ksp->will_need(bucket + nthreads);

                        ksp->load(bucket, recs);
                        std::sort(recs.begin(), recs.end(), [](const KRec& lhs, const KRec& rhs) {
                            const uint64_t lhs_key = lhs.key, rhs_key = rhs.key; return lhs_key < rhs_key; });
                        for (std::size_t j = 0; j < recs.size(); ++j)
                            // base[bucket] + j is final because keys are distinct and the bucket is sorted, so no second pass renumbers nodes.
                            ef_graph->put(base[bucket] + j, recs[j].key, recs[j].cnt, recs[j].edge_mask);

                        std::vector<KRec>().swap(recs);
                    }
                } catch (const std::exception& err) { fail_with(err); }
            });
            if (failed) throw std::runtime_error(err_msg);
        }
        ef_graph->finish_put();
        ksp.reset();
        timer.lap("count: compact graph (Elias-Fano)");

        // Edges to neighbours below the cutoff can only go once every node is known.
        if (min_count > 1) ef_graph->resolve_masks();

        *compact = std::move(ef_graph);
        timer.lap("count: edge masks");
        return;
    }

    graph.init_table(n_solid);
    count_pass(true);
    graph.finish_insert(n_solid);     // solid k-mers are distinct
    bufs.clear();
    timer.lap("count: partitions (insert)");

    // Edges whose neighbour is solid.
    // At min_count 1 every observed edge already is one, since both k-mers of an observed (K+1)-window were counted and so are both present.
    if (min_count > 1) graph.resolve_masks();
    else graph.masks_from_raw();

    timer.lap("count: edge masks");
}

void count_kmers_with_edges(const char* fasta_path, std::size_t kmer_len,
                            uint32_t min_count, unsigned nthreads,
                            std::vector<std::pair<uint64_t, uint32_t>>& out_kmers,
                            std::vector<uint8_t>& out_masks) {
    DbgGraph graph(static_cast<unsigned>(kmer_len));

    count_into_graph(fasta_path, kmer_len, min_count, nthreads, graph);
    graph.export_sorted(out_kmers, out_masks);
}

}
}
