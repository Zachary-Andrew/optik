/*
 * kmer_count.cpp - canonical k-mer counting with inline edge recording, as declared in kmer_count.hpp.
 *
 * Sources
 *   - minimizer-partitioned counting of super-k-mers from KMC 2 (Deorowicz et al. 2015), where optik keeps the partitions in memory and spill.hpp writes them to disk only past the memory budget
 *   - the edge-centric graph of Cuttlefish 2 (Khan et al. 2022), where each k-mer's edges are the (k+1)-mers observed around its occurrences
 *   - HyperLogLog (Flajolet et al. 2007) to size the tables of large partitions, Fibonacci hashing (Knuth, TAOCP vol. 3) and group prefetching (Slotin, Algorithms for Modern Hardware)
 *
 * Phases of count_into_graph()
 *   1. Heavy minimizers.
 *      A sample of the input finds the m-mers that would be the minimizer of a large share of the k-mers, as adapter read-through or poly-A would, and ranks them after every other m-mer.
 *   2. Pass 1.
 *      Readers cut the input into chunks of kChunk k-mer start positions, and workers split each chunk into super-k-mers, the maximal runs of k-mers sharing a minimizer.
 *      Each is appended 2-bit packed, with its flanking bases, to per-thread buffers of the minimizer's partition.
 *      The flanks are the observed (k+1)-mers of the first and last k-mer, so no separate (k+1)-mer pass is needed.
 *      Past the memory budget, buffers spill to disk.
 *   3. Pass 2.
 *      One thread counts each partition in a private, reused, cache-sized table from part_count_table.hpp.
 *      A sizing run counts the solid k-mers, then an insert run places them in the graph, which is the DbgGraph table, the compact DbgGraphEF, or dbg_part.hpp's partitioned compaction.
 *   4. With min_count > 1, edge masks are restricted to solid neighbours.
 */
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

namespace tpoptoa {
namespace detail {

namespace {

constexpr uint64_t kCtEmpty = ~uint64_t(0);   // never a canonical k-mer (k <= 32)
constexpr std::size_t kChunk = 1u << 18;       // k-mer starts per work item

/// Super-k-mer records are cut at this many bases (the header's length field has 10 bits).
constexpr std::size_t kMaxRecordBases = 1000;
/// Record header (u16), nbases | flank flags << kFlagShift | flank codes
/// << CODE_SHIFT.
constexpr unsigned kHdrLenMask = 0x3FFu, kHdrFlagShift = 10, kHdrCodeShift = 12;
constexpr std::size_t kHdrBytes = 2;
/// Seed that maps a minimizer to its partition (with ct_hash()).
constexpr uint64_t kPartitionSeed = 0x9E3779B97F4A7C15ULL;

/// ct_hash() - MurmurHash3's fmix64, which orders m-mers and assigns partitions.
constexpr inline uint64_t ct_hash(uint64_t x) noexcept {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

/// BaseLut - base_encode() as a 256-entry table.
/// The pass-1 scan reads every input base once through it, where a switch per base measured about 1% of the run on reads.
struct BaseLut {
    uint8_t v[256];
    constexpr BaseLut() : v{} {
        for (int c = 0; c < 256; ++c) v[c] = base_encode(static_cast<char>(c));
    }
};
constexpr BaseLut kBaseLut{};

struct Chunk {
    std::shared_ptr<const std::string> seq;
    std::size_t lo, hi;                        // k-mer start positions [lo, hi)
};

/// ChunkQueue - bounded multi-producer multi-consumer work queue under a mutex, with one push and pop per kChunk k-mers, far off the per-k-mer hot path.
class ChunkQueue {
public:
    explicit ChunkQueue(std::size_t cap) : cap_(cap) {}
    void push(Chunk c) {
        std::unique_lock<std::mutex> l(mu_);
        not_full_.wait(l, [&] { return q_.size() < cap_; });
        q_.push_back(std::move(c));
        not_empty_.notify_one();
    }
    bool pop(Chunk& c) {
        std::unique_lock<std::mutex> l(mu_);
        not_empty_.wait(l, [&] { return !q_.empty() || done_; });
        if (q_.empty()) return false;
        c = std::move(q_.front());
        q_.pop_front();
        not_full_.notify_one();
        return true;
    }
    void close() {
        std::lock_guard<std::mutex> l(mu_);
        done_ = true;
        not_empty_.notify_all();
    }
private:
    std::mutex mu_;
    std::condition_variable not_empty_, not_full_;
    std::deque<Chunk> q_;
    std::size_t cap_;
    bool done_ = false;
};

}

// All occurrences of a k-mer share its minimizer, so they meet in exactly one partition.
// Pass 1 stores about (k + #kmers - 1) / 4 + 2 bytes per super-k-mer instead of a slot per distinct k-mer in one shared table, and pass 2's inserts hit L2 instead of DRAM.
namespace {

/// SkBuf - append-only byte buffer of packed super-k-mer records.
///   Each record has a little-endian u16 header, nbases (10 bits, since
///   records are capped at kMaxRecordBases) | flank flags << 10 (bit 0 left,
///   bit 1 right) | flank codes << 12 (left | right << 2), followed by
///   ceil(nbases / 4) bytes of 2-bit bases.  The header has two bytes, not
///   four, because on read sets the spill is mostly these records.
struct SkBuf {
    std::vector<uint8_t> bytes;
    std::size_t kmers = 0;
};

/// LocalSlot - a solid k-mer of a counted partition.
struct LocalSlot { uint64_t key; uint32_t cnt; uint32_t edge_mask; };

/// Span - a run of packed super-k-mer records (in memory or read back from disk).
struct Span { const uint8_t* d; std::size_t n; };
/// SpillHdr - header of one partition buffer in a spill file.
struct SpillHdr { uint32_t p, pad; uint64_t nbytes; };
static_assert(sizeof(SpillHdr) == 16, "spill header layout");
/// Extent - a spill-file extent holding records of one bin (one thread, one flush).
struct Extent { uint64_t off, len; };
/// KRec - a solid k-mer on its way into the compact graph (key-ordered buckets), packed into 13 bytes in memory and in the spill files.
struct __attribute__((packed)) KRec { uint64_t key; uint32_t cnt; uint8_t edge_mask; };
static_assert(sizeof(KRec) == 13, "packed record");

/// choose_log_partitions() - log2 of the partition count for the inputs
/// @fasta_path, given as '\n'-joined paths.
///
/// About 64k k-mer occurrences per partition keep each private counting table near L2 size.
/// Every input file counts, gzip ones at about 4x their compressed size.
/// There are at most 2^16 partitions, since a 3.1 Gbp human genome at 2^14 had about 190k occurrences per partition, tables past L2 and a 45 s sizing pass.
/// OPTIK_LOGP=<n> sets the count to 2^n directly as a profiling aid.
unsigned choose_log_partitions(const char* fasta_path) {
    constexpr unsigned kMinLogPartitions = 6, kMaxLogPartitions = 16;
    constexpr std::size_t kInputBytesPerPartition = 65536;
    constexpr std::size_t kGzipExpansion = 4;
    unsigned logP = kMinLogPartitions;
    {
        std::size_t bytes = 0;
        for (const std::string& f : split_inputs(fasta_path)) {
            struct stat st;
            if (::stat(f.c_str(), &st) != 0 || st.st_size <= 0) continue;
            const bool gz = f.size() > 3 && (f.compare(f.size() - 3, 3, ".gz") == 0 ||
                                             f.compare(f.size() - 3, 3, ".GZ") == 0);
            bytes += static_cast<std::size_t>(st.st_size) * (gz ? kGzipExpansion : 1);
        }
        const std::size_t want = bytes / kInputBytesPerPartition;
        while (logP < kMaxLogPartitions && (std::size_t(1) << logP) < want) ++logP;
    }
    if (const char* e = std::getenv("OPTIK_LOGP")) logP = std::clamp(unsigned(std::atoi(e)), 6u, 24u);
    return logP;
}

/// HeavyMinimizers - the m-mer order of the counter, with heavy minimizers ranked last.
///
/// A few m-mers can be the minimizer of a large share of all k-mers, through adapter read-through such as Nextera or TruSeq or through poly-A.
/// The k-mers they win all go to one partition, which one thread then counts alone with a table and bin buffer of GBs.
/// sample() finds them in the start of the input, and order() ranks them after every other m-mer while keeping their own order, so a k-mer with any other m-mer in its window takes that one.
/// The order stays a fixed function of the m-mer, so all occurrences of a k-mer still land in one partition, and the result does not depend on the partitioning at all.
class HeavyMinimizers {
public:
    /// order() - rank of canonical m-mer @cm, given by ct_hash(), a full finalizer, with heavy m-mers last.
    /// A single-multiply order measured slower, since a less random order gives shorter super-k-mers and more window rescans.
    uint64_t order(uint64_t cm) const {
        const uint64_t h = ct_hash(cm) & ~kHeavyBit;
        return n_heavy_ && contains(cm) ? h | kHeavyBit : h;
    }

    /// sample() - find the heavy minimizers of k-mer length @K and minimizer length @M in up to 8 Mbp (kSample) from the start of the first input, against @P partitions.
    /// Rounds repeat because the k-mers moved off a heavy m-mer can make their next one heavy.
    void sample(const char* fasta_path, unsigned K, unsigned M, std::size_t P);

private:
    static constexpr uint64_t kHeavyBit = uint64_t(1) << 63;   // set in a heavy m-mer's order, ranking it last
    static constexpr unsigned kHeavyRounds = 3;                // sample passes (see sample())
    static constexpr std::size_t kMaxHeavy = 256;              // heavy m-mers kept at most
    static constexpr uint64_t kMinHeavyCount = 64;             // k-mers a heavy minimizer must win at least

    std::size_t slot(uint64_t cm) const {                      // Fibonacci hashing on the top heavy_lg_ bits
        return static_cast<std::size_t>((cm * kPartitionSeed) >> (64 - heavy_lg_));
    }
    bool contains(uint64_t cm) const {
        for (std::size_t h = slot(cm);; h = (h + 1) & (heavy_.size() - 1)) {
            if (heavy_[h] == cm) return true;
            if (heavy_[h] == kCtEmpty) return false;
        }
    }

    std::vector<uint64_t> heavy_ = std::vector<uint64_t>(64, kCtEmpty);   // open addressing, canonical m-mers
    unsigned heavy_lg_ = 6;                                               // log2(heavy_.size())
    std::size_t n_heavy_ = 0;
};

void HeavyMinimizers::sample(const char* fasta_path, unsigned K, unsigned M, std::size_t P) {
    const uint64_t mmask = kmer_mask(M);
    const unsigned mrc_shift = 2u * (M - 1);
    const unsigned W = K - M + 1;
    const std::vector<std::string> ins = split_inputs(fasta_path);
    std::vector<std::string> sample;
    std::size_t got = 0;
    // 8 Mbp from the start of the first input (a long record, e.g. a whole chromosome, contributes its first bases only)
    constexpr std::size_t kSample = std::size_t(8) << 20;    // bases
    struct Enough {};
    if (ins[0] != "-") {
        try {
            iterate_sequences(ins[0].c_str(), [&](const SeqRecord& rec) {
                if (rec.seq.size() >= K) {
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
    static_assert(kRing >= 32, "a window of a k <= 32 k-mer");
    std::vector<uint64_t> oh(kRing), cmv(kRing);
    PartCountTable cnt;
    for (unsigned round = 0; round < kHeavyRounds && !sample.empty(); ++round) {
        cnt.begin(got * 2 / (W + 1) + 64);
        uint64_t total = 0;
        // Runs of k-mers with one minimizer get one table update each.
        uint64_t run_cm = 0;
        uint32_t run_n = 0;
        auto flush = [&] { if (run_n) cnt.add(run_cm, 0, run_n); run_n = 0; };
        for (const std::string& r : sample) {
            uint64_t mf = 0, mr = 0;
            std::size_t mvalid = 0, n = 0;
            // window minimum kept as the pass-1 scan keeps it, rescanned only when it leaves
            uint64_t best = ~uint64_t(0);
            std::size_t best_i = 0;
            for (std::size_t i = 0; i < r.size(); ++i) {
                const uint8_t b = kBaseLut.v[uint8_t(r[i])];
                if (b == kInvalidBase) { mvalid = 0; n = 0; continue; }
                mf = ((mf << 2) | b) & mmask;
                mr = ((mr >> 2) | (static_cast<uint64_t>(b ^ 3u) << mrc_shift)) & mmask;
                if (++mvalid < M) continue;
                const uint64_t cm = mf < mr ? mf : mr;
                oh[n & (kRing - 1)] = order(cm);
                cmv[n & (kRing - 1)] = cm;
                ++n;
                if (mvalid < K) continue;
                const std::size_t last = n - 1, q = last + 1 - W;
                if (n == W || best_i < q) {
                    best = ~uint64_t(0);
                    for (std::size_t t = q; t <= last; ++t)
                        if (oh[t & (kRing - 1)] < best) { best = oh[t & (kRing - 1)]; best_i = t; }
                } else if (oh[last & (kRing - 1)] < best) {
                    best = oh[last & (kRing - 1)]; best_i = last;
                }
                const uint64_t m = cmv[best_i & (kRing - 1)];
                if (m != run_cm) flush();
                run_cm = m;
                ++run_n;
                ++total;
            }
        }
        flush();
        // A heavy m-mer is the minimizer of more k-mers than four average partitions.
        const uint64_t thr = std::max<uint64_t>(kMinHeavyCount, total * 4 / P);
        std::vector<std::pair<uint64_t, uint64_t>> add;
        cnt.for_each([&](uint64_t cm, uint32_t c, uint8_t) {
            if (c > thr && !contains(cm)) add.push_back({c, cm});
        });
        cnt.finish(0);
        if (add.empty()) break;
        std::sort(add.rbegin(), add.rend());
        for (const auto& [c, cm] : add) {
            if (n_heavy_ >= kMaxHeavy) break;
            if ((n_heavy_ + 1) * 2 > heavy_.size()) {        // grow the set (load <= 1/2)
                std::vector<uint64_t> old(heavy_.size() * 2, kCtEmpty);
                old.swap(heavy_);
                ++heavy_lg_;
                for (uint64_t x : old) if (x != kCtEmpty) {
                    std::size_t h = slot(x);
                    while (heavy_[h] != kCtEmpty) h = (h + 1) & (heavy_.size() - 1);
                    heavy_[h] = x;
                }
            }
            std::size_t h = slot(cm);
            while (heavy_[h] != kCtEmpty) h = (h + 1) & (heavy_.size() - 1);
            heavy_[h] = cm;
            ++n_heavy_;
        }
    }
    if (n_heavy_ && std::getenv("OPTIK_TIMING"))
        std::fprintf(stderr, "[timing] %zu heavy minimizers ranked last (sample of %zu bases)\n", n_heavy_, got);
}

/// warn_if_spill_disk_short() - warn when the spill directory @sdir may not hold what the counter could spill past @buf_budget bytes of buffers.
///
/// A spilling count writes about 1.6 bytes per input base, counting 2-bit super-k-mers with their k - 1 overlaps and record headers, as measured on the CF2 gut reads.
/// Bases are estimated per input kind.
/// A shortfall gives a warning with the numbers rather than an error, since the estimate is rough and ENOSPC still stops cleanly.
void warn_if_spill_disk_short(const char* fasta_path, uint64_t buf_budget, const std::string& sdir) {
    double bases = 0;
    for (const std::string& f : split_inputs(fasta_path)) {
        struct stat st;
        if (::stat(f.c_str(), &st) != 0 || st.st_size <= 0) continue;
        std::string n = f;
        std::transform(n.begin(), n.end(), n.begin(), [](unsigned char c) { return std::tolower(c); });
        const bool gz = n.size() > 3 && n.compare(n.size() - 3, 3, ".gz") == 0;
        const bool fq = n.find(".fq") != std::string::npos || n.find(".fastq") != std::string::npos;
        bases += double(st.st_size) * (gz ? (fq ? 1.6 : 3.5) : (fq ? 0.45 : 1.0));
    }
    const double need = 1.6 * bases - double(buf_budget);
    const uint64_t avail = spill::free_bytes(sdir);
    if (need > 0 && avail && double(avail) < need)
        std::fprintf(stderr, "[count] warning: the k-mer buffers may spill about %.0f GB to %s, which has "
                     "%.0f GB free; set --tmp-dir to a larger disk or raise --max-mem\n",
                     need / 1e9, sdir.c_str(), double(avail) / 1e9);
}

}

void count_into_graph(const char* fasta_path, std::size_t kk,
                      uint32_t min_count, unsigned nthreads, DbgGraph& graph,
                      std::unique_ptr<DbgGraphEF>* compact, std::unique_ptr<PartUnitigs>* part) {
    const unsigned K = static_cast<unsigned>(kk);
    nthreads = resolve_threads(nthreads);
    if (nthreads == 0) nthreads = 1;
    const uint64_t kmask = kmer_mask(K);
    const unsigned rc_shift = 2u * (K - 1);
    if (min_count == 0) min_count = 1;
    PhaseTimer pt;

    // The minimizer length is long enough to spread partitions and short enough for long super-k-mers (window w = K - M + 1).
    constexpr unsigned kMaxMinimizerLen = 13;
    const unsigned M = K <= 3 ? K : std::min(kMaxMinimizerLen, K - 3);
    const uint64_t mmask = kmer_mask(M);
    const unsigned mrc_shift = 2u * (M - 1);
    const unsigned W = K - M + 1;

    const unsigned logP = choose_log_partitions(fasta_path);
    const std::size_t P = std::size_t(1) << logP;

    std::vector<std::vector<SkBuf>> bufs(nthreads, std::vector<SkBuf>(P));
    // Heavy minimizers are ranked after every other m-mer (HeavyMinimizers).
    HeavyMinimizers heavy;
    heavy.sample(fasta_path, K, M, P);
    auto order_of = [&](uint64_t cm) { return heavy.order(cm); };
    pt.lap("count: heavy minimizers");

    ChunkQueue queue(4 * nthreads);

    // Memory budget and disk spill, as spill.hpp implements them.
    // The partition buffers grow with the input, coverage times genome, and the graph only with its distinct solid k-mers.
    // Each worker keeps its buffers in memory up to its share of the budget.
    // Past it the worker appends them to its own spill file and starts over, so the counting tables below see the same records either way.
    // Below the budget nothing is written and the counter runs exactly as before.
    bool explicit_limit = false;
    const uint64_t budget = spill::memory_budget(explicit_limit);
    // With auto, buffers get half of the available memory and all of it is the limit that buffers plus graph must stay under.
    // With --max-mem both come from the limit.
    const uint64_t limit = explicit_limit ? budget : budget * 2;
    // Partitioned compaction holds little after counting, and on fast storage a spill costs next to nothing, with a human genome taking 86 s spilled against 90 s in memory.
    // Its partition buffers therefore keep a smaller share, half of --max-mem and by default at most 1 GB.
    uint64_t buf_budget = explicit_limit ? budget / 10 * 8 : budget;
    if (part) buf_budget = explicit_limit ? budget / 2 : std::min<uint64_t>(budget, uint64_t(1) << 30);
    // OPTIK_SPILL_TEST=<bytes> sets the per-thread share directly (tests force many small spills on small inputs)
    uint64_t per_thread = std::max<uint64_t>(buf_budget / nthreads, uint64_t(4) << 20);
    if (const char* e = std::getenv("OPTIK_SPILL_TEST")) per_thread = std::strtoull(e, nullptr, 10);
    // Bins are runs of P / B consecutive partitions, the unit read back at once.
    const unsigned logB = std::min(logP, 10u);
    const std::size_t B = std::size_t(1) << logB;
    const unsigned bin_shift = logP - logB;
    const std::string sdir = spill::spill_dir();
    warn_if_spill_disk_short(fasta_path, buf_budget, sdir);
    std::vector<spill::File> files(nthreads);
    std::vector<std::vector<std::vector<Extent>>> ext(nthreads, std::vector<std::vector<Extent>>(B));
    std::vector<uint64_t> held(nthreads, 0);          // buffer capacity per thread
    std::atomic<bool> failed{false};
    std::mutex err_mu;
    std::string err_msg;
    std::once_flag spill_note;
    auto fail_with = [&](const std::exception& e) {
        // the message is written before the flag is raised (readers check the flag, then read the message)
        std::lock_guard<std::mutex> l(err_mu);
        if (!failed.load()) { err_msg = e.what(); failed.store(true, std::memory_order_release); }
    };
    // Write thread t's buffers of partitions [p0, p1) to its spill file and free them.
    // Records keep their partition, and extents are indexed by bin.
    auto spill_range = [&](unsigned t, std::size_t p0, std::size_t p1) {
        std::call_once(spill_note, [&] {
            std::fprintf(stderr, "[count] k-mer buffers exceed the memory budget (%.2f GB%s): "
                         "spilling them to %s\n", buf_budget / 1073741824.0,
                         explicit_limit ? ", --max-mem" : ", auto", sdir.c_str());
            if (spill::ram_backed(sdir))
                std::fprintf(stderr, "[count] warning: %s is RAM-backed (tmpfs); spilling there "
                             "saves no memory -- set --tmp-dir to a disk\n", sdir.c_str());
        });
        files[t].open(sdir);
        auto& mine = bufs[t];
        std::size_t n = 0;
        for (std::size_t p = p0; p < p1; ++p) n += !mine[p].bytes.empty();
        if (n == 0) return;
        std::vector<SpillHdr> hdr; hdr.reserve(n);
        std::vector<struct iovec> iov; iov.reserve(2 * n);
        uint64_t off = files[t].size();
        std::size_t cur_bin = SIZE_MAX;
        for (std::size_t p = p0; p < p1; ++p) {
            auto& by = mine[p].bytes;
            if (by.empty()) continue;
            const std::size_t b = p >> bin_shift;
            if (b != cur_bin) { ext[t][b].push_back({off, 0}); cur_bin = b; }
            hdr.push_back({static_cast<uint32_t>(p), 0, by.size()});
            iov.push_back({&hdr.back(), sizeof(SpillHdr)});
            iov.push_back({by.data(), by.size()});
            ext[t][b].back().len += sizeof(SpillHdr) + by.size();
            off += sizeof(SpillHdr) + by.size();
        }
        files[t].append(iov);
        for (std::size_t p = p0; p < p1; ++p) {
            held[t] -= mine[p].bytes.capacity();
            std::vector<uint8_t>().swap(mine[p].bytes);   // kmers stays, since it counts spilled ones too
        }
    };

    // Pass 1 puts super-k-mers into partition buffers.
    // `single` processes that one chunk on the calling thread and returns (with -t 1 the reader does the work itself, so one thread runs, not two).
    auto worker = [&](unsigned tid, const Chunk* single = nullptr) {
        auto& mine = bufs[tid];
        uint64_t& hold = held[tid];
        Chunk c;
        std::vector<uint64_t> mh;         // canonical m-mer hashes of the run
        std::vector<uint8_t> cd;          // 2-bit codes of the chunk's scanned bases
        uint8_t pk[(kMaxRecordBases + 3) / 4 + kHdrBytes];   // one packed record
        auto body = [&](const Chunk& c) {
            if (failed.load(std::memory_order_relaxed)) return;   // drain the queue
            const std::string& s = *c.seq;
            const std::size_t n = s.size();
            const std::size_t scan0 = c.lo > 0 ? c.lo - 1 : 0;
            const std::size_t scan_end = std::min(n, c.hi + K);   // exclusive
            cd.resize(scan_end - scan0);
            // state of the current super-k-mer
            std::size_t sk_first = 0, sk_last = 0;   // start positions of first/last k-mer
            uint64_t sk_min = 0;
            bool open = false;
            auto emit = [&]() {
                if (!open) return;
                const std::size_t a = sk_first, e = sk_last + K;         // bases [a, e)
                const std::size_t nb = e - a;
                SkBuf& sb = mine[ct_hash(sk_min ^ kPartitionSeed) & (P - 1)];
                uint8_t flags = 0, codes = 0;
                if (a > 0) { const uint8_t b = kBaseLut.v[uint8_t(s[a - 1])]; if (b != kInvalidBase) { flags |= 1; codes |= b; } }
                if (e < n) { const uint8_t b = kBaseLut.v[uint8_t(s[e])]; if (b != kInvalidBase) { flags |= 2; codes |= b << 2; } }
                // The record is packed on the stack from the codes the scan already computed, then appended in one copy without a zero fill.
                const std::size_t nby = (nb + 3) / 4;
                const uint16_t hd = static_cast<uint16_t>(nb | (unsigned(flags) << kHdrFlagShift) | (unsigned(codes) << kHdrCodeShift));
                pk[0] = static_cast<uint8_t>(hd & 0xFF);
                pk[1] = static_cast<uint8_t>(hd >> 8);
                const uint8_t* q = cd.data() + (a - scan0);
                std::size_t t = 0, j = kHdrBytes;
                for (; t + 4 <= nb; t += 4)
                    pk[j++] = static_cast<uint8_t>(q[t] | (q[t + 1] << 2) | (q[t + 2] << 4) | (q[t + 3] << 6));
                if (t < nb) {
                    uint8_t x = 0;
                    for (unsigned r = 0; t < nb; ++t, ++r) x |= static_cast<uint8_t>(q[t] << (2 * r));
                    pk[j++] = x;
                }
                const std::size_t cap0 = sb.bytes.capacity();
                sb.bytes.insert(sb.bytes.end(), pk, pk + kHdrBytes + nby);
                if (sb.bytes.capacity() != cap0) hold += sb.bytes.capacity() - cap0;
                sb.kmers += sk_last - sk_first + 1;
                open = false;
            };
            uint64_t mf = 0, mr = 0;            // rolling m-mer on both strands
            std::size_t mvalid = 0;             // valid bases in the current run
            std::size_t run0 = scan0;           // start of the current ACGT run
            mh.clear();
            // window minimum over the W m-mers of the current k-mer
            uint64_t best = ~uint64_t(0);
            std::size_t best_i = 0;             // index in mh of the minimum
            for (std::size_t i = scan0; i < scan_end; ++i) {
                const uint8_t b = kBaseLut.v[uint8_t(s[i])];
                cd[i - scan0] = b;
                if (b == kInvalidBase) {
                    emit();
                    mvalid = 0; mh.clear(); run0 = i + 1;
                    continue;
                }
                mf = ((mf << 2) | b) & mmask;
                mr = ((mr >> 2) | (static_cast<uint64_t>(b ^ 3u) << mrc_shift)) & mmask;
                ++mvalid;
                if (mvalid < M) continue;
                mh.push_back(order_of(mf < mr ? mf : mr));   // m-mer ending at i
                if (mvalid < K) continue;
                const std::size_t start = i + 1 - K;
                if (start >= c.hi) break;
                // m-mers of this k-mer are mh[q .. q + W - 1], with q = start - run0
                const std::size_t q = start - run0, last = q + W - 1;
                const uint64_t hn = mh[last];
                if (start == run0 || best_i < q) {         // first k-mer, or the minimum left, so rescan
                    best = ~uint64_t(0);
                    for (std::size_t t = q; t <= last; ++t) {
                        const bool lt = mh[t] < best;     // leftmost minimum, so ties keep the first
                        best = lt ? mh[t] : best;
                        best_i = lt ? t : best_i;
                    }
                } else if (hn < best) {
                    best = hn; best_i = last;
                }
                if (start < c.lo) continue;               // context k-mer only
                if (open && best == sk_min) {
                    sk_last = start;
                } else {
                    emit();
                    open = true; sk_first = sk_last = start; sk_min = best;
                }
                if (sk_last - sk_first + K >= kMaxRecordBases) emit();   // cap record length
            }
            emit();
            if (hold > per_thread) {
                try { spill_range(tid, 0, P); } catch (const std::exception& e) { fail_with(e); }
            }
        };
        if (single) { body(*single); return; }
        while (queue.pop(c)) body(c);
    };

    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    // Several input files are read at once, each by its own reader with its own batch, since one gzip stream decodes only so fast and reading a multi-file read set file after file left the workers waiting on one decoder.
    // Record order cannot change the counts, because workers already append to partitions in any order, so the output is the same.
    // Readers come out of -t, as Cuttlefish 3's decoders come out of its --threads.
    // Of the t threads, r = t/4 are readers, at most one per file, each inflating its own .gz in-thread through fasta_reader.hpp, and t - r are pass-1 super-k-mer workers.
    const std::vector<std::string> inputs = split_inputs(fasta_path);
    const unsigned nreaders = static_cast<unsigned>(
        std::min<std::size_t>(inputs.size(), std::max(1u, nthreads / 4)));
    // With -t 1 there are no worker threads, and the one reader processes its chunks itself.
    const bool inline_work = nthreads == 1;
    const unsigned nworkers = inline_work ? 0 : nthreads > nreaders ? nthreads - nreaders : 1;
    for (unsigned t = 0; t < nworkers; ++t) workers.emplace_back(worker, t, nullptr);

    // Short records (reads) are batched into one work item joined by 'N', which resets the scanner and blanks flanks, so the result is identical to separate records.
    auto push_seq = [&](std::shared_ptr<const std::string> sp) {
        const std::size_t nstarts = sp->size() - K + 1;
        for (std::size_t lo = 0; lo < nstarts; lo += kChunk) {
            const Chunk ch{sp, lo, std::min(nstarts, lo + kChunk)};
            if (inline_work) worker(0, &ch);
            else queue.push(ch);
        }
    };
    std::atomic<std::size_t> next_input{0};
    std::exception_ptr read_err;
    std::mutex read_mu;
    std::atomic<bool> read_failed{false};
    auto reader = [&]() {
        std::string batch;
        batch.reserve(kChunk + 4096);
        auto flush_batch = [&]() {
            if (batch.size() >= K) push_seq(std::make_shared<const std::string>(std::move(batch)));
            batch = std::string();
            batch.reserve(kChunk + 4096);
        };
        try {
            for (std::size_t f; (f = next_input.fetch_add(1)) < inputs.size();) {
                iterate_sequences(inputs[f].c_str(), [&](const SeqRecord& rec) {
                    if (failed.load(std::memory_order_acquire))      // a spill write failed, so stop reading
                        throw std::runtime_error(err_msg);
                    if (read_failed.load(std::memory_order_relaxed))   // another reader failed
                        throw std::runtime_error("input read stopped");
                    if (rec.seq.size() < K) return;
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
            std::lock_guard<std::mutex> l(read_mu);
            if (!read_err) read_err = std::current_exception();
            read_failed.store(true);
        }
    };
    // A reader error (missing file, I/O error) must not unwind past running worker threads (std::terminate), so the queue is closed and the workers joined before rethrowing.
    {
        std::vector<std::thread> readers;
        for (unsigned r = 1; r < nreaders; ++r) readers.emplace_back(reader);
        reader();
        for (auto& r : readers) r.join();
    }
    if (read_err) {
        queue.close();
        for (auto& w : workers) w.join();
        std::rethrow_exception(read_err);
    }
    queue.close();
    for (auto& w : workers) w.join();
    if (failed) throw std::runtime_error(err_msg);
    pt.lap("count: parse + super-k-mers");

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
    for (auto& f : files) spilled = spilled || f.is_open();
    // Spilled buffers were freed into malloc's free lists (many small vectors, below the mmap threshold), so those pages are handed back to the OS and the resident set is what is still held.
    if (spilled) malloc_trim(0);
    auto bin_bytes = [&](std::size_t b) {
        uint64_t n = 0;
        for (unsigned t = 0; t < nthreads; ++t) for (const Extent& e : ext[t][b]) n += e.len;
        return n;
    };
    uint64_t n_distinct = 0, n_occ = 0;        // from the sizing pass, declared early for count_pass
    // With the compact DbgGraphEF graph, solid k-mers go to key-ordered buckets.
    bool ef_mode = false;
    const unsigned ef_bits = std::min(12u, 2 * K);
    const unsigned ef_shift = 2 * K - ef_bits;
    std::unique_ptr<spill::BucketSpill<KRec>> ksp;
    // Partitioned compaction (dbg_part.hpp) uses no global table at all.
    const bool part_mode = part != nullptr;
    std::unique_ptr<PartBuilder> pb;
    // Work is handed out largest first by k-mer occurrences (bins after a spill, otherwise partitions).
    // A few partitions are far heavier than the rest (adapter read-through, poly-A), and one taken last would leave a single thread working while the others wait.
    std::vector<std::size_t> order(spilled ? B : P);
    {
        std::vector<uint64_t> wt(order.size(), 0);
        for (unsigned w = 0; w < nthreads; ++w)
            for (std::size_t p = 0; p < P; ++p) wt[spilled ? p >> bin_shift : p] += bufs[w][p].kmers;
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return wt[a] > wt[b]; });
    }
    // With -c > 1 and partitioned compaction, the first pass already has each partition's solid k-mers with their in-partition edge masks, which is all the compaction needs (it applies the cross-partition claims itself).
    // The first pass therefore saves them, in memory or after a spill appended to a per-thread file, and the compaction reads them back instead of reading every super-k-mer again and recounting every occurrence.
    // On reads the solid k-mers are a small fraction of the occurrences.
    const bool save_solid = part_mode && min_count > 1;
    // Spilled lists are sorted by key and written compactly as varint key deltas, varint counts and one mask byte (about 5 bytes a k-mer instead of 16).
    struct SolidRef { uint64_t off = 0, n = 0, bytes = 0; uint32_t t = 0; };
    std::vector<SolidRef> sref(save_solid ? P : 0);
    std::vector<std::vector<PartNode>> smem(save_solid && !spilled ? P : 0);
    std::vector<spill::File> sfile(save_solid && spilled ? nthreads : 0);
    std::vector<std::vector<uint8_t>> swbuf(save_solid && spilled ? nthreads : 0);
    constexpr std::size_t kSwbuf = std::size_t(8) << 20;
    auto flush_solid = [&](unsigned t) {
        if (swbuf.empty() || swbuf[t].empty()) return;
        sfile[t].open(sdir);
        std::vector<struct iovec> iov{{swbuf[t].data(), swbuf[t].size()}};
        sfile[t].append(iov);
        swbuf[t].clear();
    };
    auto put_varint = [](std::vector<uint8_t>& o, uint64_t x) {
        while (x >= 0x80) { o.push_back(static_cast<uint8_t>(x | 0x80)); x >>= 7; }
        o.push_back(static_cast<uint8_t>(x));
    };
    auto save_nodes = [&](unsigned t, std::size_t p, std::vector<PartNode>& v) {
        if (v.empty()) return;
        if (!spilled) { smem[p] = v; return; }
        std::sort(v.begin(), v.end(), [](const PartNode& a, const PartNode& b) { return a.key < b.key; });
        std::vector<uint8_t>& o = swbuf[t];
        const uint64_t at = o.size();
        uint64_t prev = 0;
        for (const PartNode& nd : v) {
            put_varint(o, nd.key - prev);
            prev = nd.key;
            put_varint(o, nd.cnt);
            o.push_back(static_cast<uint8_t>(nd.edge_mask));
        }
        sref[p] = SolidRef{sfile[t].size() + at, v.size(), o.size() - at, t};
        if (o.size() >= kSwbuf) flush_solid(t);
    };
    auto count_pass = [&](bool insert) {
        std::atomic<std::size_t> next{0};
        std::atomic<std::size_t> total{0};
        std::atomic<uint64_t> all_distinct{0}, all_occ{0};
        std::atomic<uint64_t> max_tab{0}, max_bin{0}, max_occ{0};   // OPTIK_TIMING report
        std::atomic<uint64_t> all_deep{0};
        std::atomic<uint64_t> all_lv[8] = {};
        tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
          try {
            // One table per thread is reused across partitions.
            // Only the slots a partition touched are visited and reset afterwards, so a partition costs time in its occurrences rather than in the table capacity.
            PartCountTable T;                 // reused across partitions
            double ratio = 0.125;             // distinct / occurrences of the last partition
            std::array<double, 32> rat{};     // recent ratios of ordinary partitions
            std::size_t nrat = 0;
            constexpr unsigned kHllP = 12;    // 2^12 HyperLogLog registers
            std::vector<uint8_t> hll;
            std::vector<LocalSlot> solid;
            std::size_t mine = 0;
            uint64_t my_distinct = 0, my_occ = 0, my_max_bin = 0, my_max_occ = 0, my_max_tab = 0, my_deep = 0;
            uint64_t my_lv[8] = {0};
            std::vector<Span> spans;
            std::vector<PartNode> pnodes;
            bool release_solid = false;
            constexpr std::size_t kBigTabBytes = std::size_t(128) << 20;
            constexpr std::size_t kBigBin = std::size_t(64) << 20;  // bytes
            auto count_partition = [&](std::size_t p, const Span* sp, std::size_t ns) {
                std::size_t occ = 0;
                for (unsigned w = 0; w < nthreads; ++w) occ += bufs[w][p].kmers;
                if (occ == 0) return;
                my_occ += occ;
                my_max_occ = std::max<uint64_t>(my_max_occ, occ);
                // f(canonical, edge bits) for every k-mer occurrence of the partition's records
                auto each_kmer = [&](auto&& f) {
                    for (std::size_t si = 0; si < ns; ++si) {
                        const uint8_t* const by = sp[si].d;
                        const std::size_t by_n = sp[si].n;
                        std::size_t o = 0;
                        while (o < by_n) {
                            const unsigned hd = by[o] | (unsigned(by[o + 1]) << 8);
                            const std::size_t nb = hd & kHdrLenMask;
                            const uint8_t flags = static_cast<uint8_t>((hd >> kHdrFlagShift) & 3u);
                            const uint8_t codes = static_cast<uint8_t>(hd >> kHdrCodeShift);
                            const uint8_t* d = by + o + kHdrBytes;
                            auto base = [&](std::size_t t) -> uint8_t { return (d[t >> 2] >> ((t & 3) * 2)) & 3u; };
                            uint64_t fwd = 0, rev = 0;
                            for (std::size_t t = 0; t < nb; ++t) {
                                const uint8_t b = base(t);
                                fwd = ((fwd << 2) | b) & kmask;
                                rev = ((rev >> 2) | (static_cast<uint64_t>(b ^ 3u) << rc_shift)) & kmask;
                                if (t + 1 < K) continue;
                                const std::size_t st = t + 1 - K;
                                const uint8_t pb = st > 0 ? base(st - 1) : ((flags & 1) ? (codes & 3u) : kInvalidBase);
                                const uint8_t nx = t + 1 < nb ? base(t + 1) : ((flags & 2) ? ((codes >> 2) & 3u) : kInvalidBase);
                                const uint64_t cano = fwd < rev ? fwd : rev;
                                uint8_t in_b = 0, out_b = 0;
                                if (fwd == cano) {
                                    if (nx != kInvalidBase) out_b |= 1u << nx;
                                    if (pb != kInvalidBase) in_b  |= 1u << pb;
                                }
                                if (rev == cano) {
                                    if (nx != kInvalidBase) in_b  |= 1u << (nx ^ 3u);
                                    if (pb != kInvalidBase) out_b |= 1u << (pb ^ 3u);
                                }
                                f(cano, static_cast<uint8_t>((in_b << 4) | out_b));
                            }
                            o += kHdrBytes + (nb + 3) / 4;
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
                        const uint64_t h = ct_hash(cano ^ kHllSeed);
                        const std::size_t r = h >> (64 - kHllP);
                        const uint8_t z = static_cast<uint8_t>(__builtin_clzll((h << kHllP) | (uint64_t(1) << (kHllP - 1))) + 1);
                        if (z > hll[r]) hll[r] = z;
                    });
                    double sum = 0; std::size_t zeros = 0;
                    for (uint8_t v : hll) { sum += std::ldexp(1.0, -int(v)); zeros += v == 0; }
                    constexpr double mreg = double(std::size_t(1) << kHllP);
                    // the raw estimate alpha_m * m^2 / sum 2^-M[j] with alpha_m = 0.7213 / (1 + 1.079 / m), and linear counting below 2.5 m (Flajolet et al. 2007)
                    est = 0.7213 / (1 + 1.079 / mreg) * mreg * mreg / sum;
                    if (est < 2.5 * mreg && zeros) est = mreg * std::log(mreg / double(zeros));   // linear counting
                    est *= 1.05;                                // 5% margin
                } else {
                    // The median of this thread's recent ratios is used, so a single highly repetitive partition (ratio about 1e-4) cannot shrink the next ordinary one's table.
                    double med = 0.125;
                    if (nrat) {
                        std::array<double, 32> tmp;
                        const std::size_t n = std::min<std::size_t>(nrat, tmp.size());
                        std::copy_n(rat.begin(), n, tmp.begin());
                        std::nth_element(tmp.begin(), tmp.begin() + n / 2, tmp.begin() + n);
                        med = tmp[n / 2];
                    }
                    est = static_cast<double>(occ) * med * 1.5;
                }
                T.begin(static_cast<std::size_t>(est) + 16);
                {
                    // Group prefetching keeps a ring of kRing k-mers whose lines are in flight.
                    // Each new k-mer's lines are prefetched and the oldest is inserted, and the order does not change the counts.
                    constexpr std::size_t kRing = 16;            // a power of two
                    uint64_t rk[kRing], rh[kRing];
                    uint8_t rb[kRing];
                    std::size_t wr = 0, rd = 0;                 // FIFO counters
                    each_kmer([&](uint64_t cano, uint8_t bits) {
                        const uint64_t h = PartCountTable::hash(cano);
                        T.prefetch(h);
                        if (wr - rd == kRing) {
                            const std::size_t j = rd++ & (kRing - 1);
                            T.add_hashed(rk[j], rh[j], rb[j]);
                        }
                        const std::size_t j = wr++ & (kRing - 1);
                        rk[j] = cano; rh[j] = h; rb[j] = bits;
                    });
                    while (rd < wr) {
                        const std::size_t j = rd++ & (kRing - 1);
                        T.add_hashed(rk[j], rh[j], rb[j]);
                    }
                }
                const std::size_t nd = T.size();
                my_max_tab = std::max<uint64_t>(my_max_tab, T.bytes());
                my_deep += T.deep();
                my_lv[std::min(T.levels(), 7u)]++;
                ratio = std::min(1.0, static_cast<double>(nd) / static_cast<double>(occ));
                if (occ >= 1024 && occ <= kHllMinOccurrences) rat[nrat++ % rat.size()] = ratio;
                my_distinct += nd;
                solid.clear();
                if (part_mode) {
                    // The partition's table gives exact in-partition masks, where a neighbour counted here below min_count is no edge.
                    // Cross-partition edges are claimed in the first pass and handed over in the second.
                    const unsigned shK = 2u * (K - 1);
                    auto nbr = [&](uint64_t c, unsigned b) {
                        const uint64_t x = b < 4 ? (((c << 2) | b) & kmask) : ((uint64_t(b - 4) << shK) | (c >> 2));
                        const uint64_t r = reverse_complement(x, K);
                        return x < r ? x : r;
                    };
                    const bool keep = insert || save_solid;
                    if (keep) pnodes.clear();
                    T.for_each([&](uint64_t key, uint32_t cnt, uint8_t edge_mask) {
                        if (cnt < min_count) return;
                        ++mine;
                        uint32_t m = edge_mask;
                        if (min_count > 1)
                            for (unsigned b = 0; b < 8; ++b) {
                                if (!(m & (1u << b))) continue;
                                const uint32_t y = T.find(nbr(key, b));
                                if (y != PartCountTable::kNotFound) { if (y < min_count) m &= ~(1u << b); }
                                else if (!insert) pb->claim(t, static_cast<uint32_t>(p), pb->mask_edge(key, b));
                            }
                        if (keep) pnodes.push_back(PartNode{key, cnt, m});
                    });
                    if (insert) pb->partition(t, static_cast<uint32_t>(p), pnodes);
                    else if (save_solid) save_nodes(t, p, pnodes);
                } else {
                    T.for_each([&](uint64_t key, uint32_t cnt, uint8_t edge_mask) {
                        if (cnt < min_count) return;
                        ++mine;
                        if (insert) solid.push_back(LocalSlot{key, cnt, edge_mask});
                    });
                }
                // A heavy partition (one minimizer shared by many k-mers, as with adapter read-through or poly-A) grows this thread's table far past the usual size.
                // If it were kept for reuse every thread would end the pass holding its largest, so past the cap it is released, after its last use below, since solid is still needed there.
                if (T.bytes() > kBigTabBytes) release_solid = true;
                T.finish(kBigTabBytes);
                // One thread owns each k-mer through its partition, so claiming takes one CAS on the key word and then a plain store of the node word.
                // A locked CAS keeps later loads from running ahead, so each claim would wait out its own DRAM miss.
                // The home line kPrefetchAhead claims ahead is prefetched for writing instead, and the misses overlap.
                constexpr std::size_t kPrefetchAhead = 16;
                if (ef_mode) {
                    for (const LocalSlot& sl : solid)
                        ksp->add(t, sl.key >> ef_shift, KRec{sl.key, sl.cnt, static_cast<uint8_t>(sl.edge_mask)});
                } else {
                    for (std::size_t j = 0; j < solid.size(); ++j) {
                        if (j + kPrefetchAhead < solid.size()) graph.prefetch_w(solid[j + kPrefetchAhead].key);
                        graph.put(graph.claim(solid[j].key), solid[j].cnt,
                                  static_cast<uint8_t>(solid[j].edge_mask));
                    }
                }
                if (release_solid) {
                    std::vector<LocalSlot>().swap(solid);
                    std::vector<PartNode>().swap(pnodes);
                    release_solid = false;
                }
            };
            auto mem_spans = [&](std::size_t p) {
                for (unsigned w = 0; w < nthreads; ++w)
                    if (!bufs[w][p].bytes.empty())
                        spans.push_back({bufs[w][p].bytes.data(), bufs[w][p].bytes.size()});
            };
            auto free_mem = [&](std::size_t p) {
                if (!insert && !save_solid) return;   // a later pass reads them again
                for (unsigned w = 0; w < nthreads; ++w)           // free as we go
                    std::vector<uint8_t>().swap(bufs[w][p].bytes);
            };
            if (!spilled) {
                for (std::size_t i, p; (i = next.fetch_add(1)) < P && (p = order[i], true);) {
                    spans.clear();
                    mem_spans(p);
                    count_partition(p, spans.data(), spans.size());
                    free_mem(p);
                }
            } else {
                const std::size_t PB = P >> logB;   // partitions per bin
                std::vector<uint8_t> binbuf;
                std::vector<std::vector<Span>> disk(PB);
                for (std::size_t i, b; (i = next.fetch_add(1)) < B && (b = order[i], true);) {
                    if (failed.load(std::memory_order_relaxed)) break;
                    // Start the reads of the bins the other workers take next.
                    if (const std::size_t ni = i + nthreads; ni < B)
                        for (unsigned w = 0; w < nthreads; ++w)
                            for (const Extent& e : ext[w][order[ni]]) files[w].will_need(e.off, e.len);
                    binbuf.resize(bin_bytes(b));
                    my_max_bin = std::max<uint64_t>(my_max_bin, binbuf.size());
                    uint64_t at = 0;
                    for (unsigned w = 0; w < nthreads; ++w)
                        for (const Extent& e : ext[w][b]) {
                            files[w].read(e.off, e.len, binbuf.data() + at);
                            at += e.len;
                        }
                    for (auto& v : disk) v.clear();
                    for (uint64_t o = 0; o < at;) {
                        SpillHdr h;
                        std::memcpy(&h, binbuf.data() + o, sizeof h);
                        if ((h.p >> bin_shift) != b || o + sizeof h + h.nbytes > at)
                            throw std::runtime_error("spill file corrupt (bad record header)");
                        disk[h.p - (b << bin_shift)].push_back({binbuf.data() + o + sizeof h, h.nbytes});
                        o += sizeof h + h.nbytes;
                    }
                    for (std::size_t j = 0; j < PB; ++j) {
                        const std::size_t p = (b << bin_shift) + j;
                        spans.assign(disk[j].begin(), disk[j].end());
                        mem_spans(p);
                        count_partition(p, spans.data(), spans.size());
                        free_mem(p);
                    }
                    // This is the last pass over the bin, so its disk blocks go back.
                    if (insert || save_solid)
                        for (unsigned w = 0; w < nthreads; ++w)
                            for (const Extent& e : ext[w][b]) files[w].punch(e.off, e.len);
                    if (binbuf.capacity() > kBigBin) std::vector<uint8_t>().swap(binbuf);
                }
            }
            if (!insert && save_solid) flush_solid(t);
            all_deep.fetch_add(my_deep);
            for (int i = 0; i < 8; ++i) all_lv[i].fetch_add(my_lv[i]);
            {
                uint64_t v = my_max_tab, c = max_tab.load();
                while (v > c && !max_tab.compare_exchange_weak(c, v)) {}
                v = my_max_bin; c = max_bin.load();
                while (v > c && !max_bin.compare_exchange_weak(c, v)) {}
                v = my_max_occ; c = max_occ.load();
                while (v > c && !max_occ.compare_exchange_weak(c, v)) {}
            }
            total.fetch_add(mine);
            all_distinct.fetch_add(my_distinct);
            all_occ.fetch_add(my_occ);
          } catch (const std::exception& e) { fail_with(e); }
        });
        if (failed) throw std::runtime_error(err_msg);
        if (!insert) { n_distinct = all_distinct.load(); n_occ = all_occ.load(); }
        if (std::getenv("OPTIK_TIMING"))
            std::fprintf(stderr, "[timing] count pass: P = 2^%u, largest partition %llu occurrences, "
                         "per-thread table up to %.0f MB, bin buffer up to %.0f MB\n", logP,
                         static_cast<unsigned long long>(max_occ.load()), max_tab.load() / 1048576.0,
                         max_bin.load() / 1048576.0);
        if (std::getenv("OPTIK_TIMING")) {
            std::fprintf(stderr, "[timing] count pass: %llu of %llu distinct keys below level 0; partitions by levels:",
                         static_cast<unsigned long long>(all_deep.load()), static_cast<unsigned long long>(all_distinct.load()));
            for (int i = 1; i < 8; ++i) std::fprintf(stderr, " %d:%llu", i, static_cast<unsigned long long>(all_lv[i].load()));
            std::fprintf(stderr, "\n");
        }
        return total.load();
    };
    if (part_mode) {
        // At -c 1 every counted k-mer is solid and no sizing pass is needed.
        // At -c > 1 the first pass claims the cross-partition edges.
        const uint64_t res = spill::resident_bytes();
        const uint64_t room = limit > res ? limit - res : 0;
        // Fragment and end-key streams keep 4-16 MB per thread in memory and the rest on disk, read back once and sequentially by bucket
        uint64_t per = std::clamp<uint64_t>(room / 16 / nthreads, uint64_t(4) << 20, uint64_t(16) << 20);
        if (const char* e = std::getenv("OPTIK_SPILL_TEST")) per = std::strtoull(e, nullptr, 10);
        pb = std::make_unique<PartBuilder>(K, nthreads, P, min_count > 1, per, sdir);
        if (min_count > 1) {
            const std::size_t ns = count_pass(false);
            if (n_distinct > 0 && n_occ < n_distinct + n_distinct / 2 && ns < n_distinct / 2)
                std::fprintf(stderr, "[count] warning: mean k-mer coverage is %.2fx, so -c %u kept %zu of %llu "
                             "distinct k-mers (%.1f%%); for a genome or assembly use -c 1\n",
                             double(n_occ) / double(n_distinct), min_count, ns,
                             static_cast<unsigned long long>(n_distinct), 100.0 * double(ns) / double(n_distinct));
            pt.lap("count: partitions (first pass: edge claims)");
            pb->match_claims();
            pt.lap("count: claims matched");
        }
        if (!save_solid) {
            count_pass(true);
        } else {
            // the saved solid lists, largest first, each compacted by one thread
            bufs.clear();
            malloc_trim(0);
            std::vector<std::size_t> porder(P);
            for (std::size_t p = 0; p < P; ++p) porder[p] = p;
            auto nsolid = [&](std::size_t p) { return spilled ? sref[p].n : smem[p].size(); };
            std::stable_sort(porder.begin(), porder.end(),
                             [&](std::size_t a, std::size_t b) { return nsolid(a) > nsolid(b); });
            std::atomic<std::size_t> next{0};
            tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                try {
                    std::vector<PartNode> nodes;
                    std::vector<uint8_t> enc;
                    for (std::size_t i; (i = next.fetch_add(1)) < P;) {
                        if (failed.load(std::memory_order_relaxed)) break;
                        const std::size_t p = porder[i];
                        if (nsolid(p) == 0) break;            // sorted, so the rest are empty
                        if (spilled) {
                            const SolidRef& r = sref[p];
                            enc.resize(r.bytes);
                            sfile[r.t].read(r.off, r.bytes, enc.data());
                            nodes.resize(r.n);
                            const uint8_t* q = enc.data();
                            const uint8_t* const qe = q + r.bytes;
                            auto get = [&]() {
                                uint64_t x = 0;
                                for (unsigned sh = 0; q < qe; sh += 7) {
                                    const uint8_t c = *q++;
                                    x |= uint64_t(c & 0x7F) << sh;
                                    if (!(c & 0x80)) return x;
                                }
                                throw std::runtime_error("spill file corrupt (solid list)");
                            };
                            uint64_t key = 0;
                            for (PartNode& nd : nodes) {
                                key += get();
                                nd.key = key;
                                nd.cnt = static_cast<uint32_t>(get());
                                if (q >= qe) throw std::runtime_error("spill file corrupt (solid list)");
                                nd.edge_mask = *q++;
                            }
                            sfile[r.t].punch(r.off, r.bytes);          // read once, so give it back
                            if (enc.capacity() > (std::size_t(64) << 20)) std::vector<uint8_t>().swap(enc);
                        } else {
                            nodes.swap(smem[p]);
                            std::vector<PartNode>().swap(smem[p]);
                        }
                        pb->partition(t, static_cast<uint32_t>(p), nodes);
                        if (nodes.capacity() > (std::size_t(1) << 20)) std::vector<PartNode>().swap(nodes);
                    }
                } catch (const std::exception& e) { fail_with(e); }
            });
            if (failed) throw std::runtime_error(err_msg);
        }
        bufs.clear();
        // The partition buffers were many small vectors, and their freed pages sit in malloc's free lists until handed back.
        malloc_trim(0);
        pt.lap("count: partitions (compaction)");
        *part = pb->finish(pt);
        return;
    }
    const std::size_t n_solid = count_pass(false);
    // -c > 1 on single-copy input such as a genome, an assembly or unitigs silently keeps only repeats, since nearly every k-mer occurs once.
    // Reads have coverage well above 1, so a warning is printed when the mean occurrence per distinct k-mer says the input is not reads.
    if (min_count > 1 && n_distinct > 0 && n_occ < n_distinct + n_distinct / 2 &&
        n_solid < n_distinct / 2)
        std::fprintf(stderr, "[count] warning: mean k-mer coverage is %.2fx, so -c %u kept %zu of %llu "
                     "distinct k-mers (%.1f%%); for a genome or assembly use -c 1\n",
                     double(n_occ) / double(n_distinct), min_count, n_solid,
                     static_cast<unsigned long long>(n_distinct), 100.0 * double(n_solid) / double(n_distinct));
    pt.lap("count: partitions (sizing)");

    // The graph is allocated before the buffers are released, since freed buffer memory is not reused for the graph's pages.
    // If the resident set plus the graph plus the bins read back would pass the limit, buffers are spilled, from the last partitions back, until it fits.
    // This is measured rather than estimated, because the resident set includes malloc overhead and what fragmentation keeps.
    // A graph larger than the limit by itself cannot be helped by spilling.
    {
        const uint64_t gbytes = DbgGraph::table_bytes(n_solid);
        auto max_bin = [&] {
            uint64_t m = 0;
            for (std::size_t b = 0; b < B; ++b) m = std::max(m, bin_bytes(b));
            return m;
        };
        // 5% headroom for the per-thread counting tables and the reader
        auto need = [&] {
            const uint64_t n = gbytes + spill::resident_bytes() + (spilled ? max_bin() * nthreads : 0);
            return n + n / 20;
        };
        uint64_t nd = need();
        if (pt.on)
            std::fprintf(stderr, "[timing] count: memory plan: graph %.0f MB + resident %.0f MB "
                         "(+ bins read back) = %.0f MB vs limit %.0f MB (spilled %s)\n",
                         gbytes / 1048576.0, spill::resident_bytes() / 1048576.0, nd / 1048576.0,
                         limit / 1048576.0, spilled ? "yes" : "no");
        // The table does not fit, so the compact graph is used (Elias-Fano keys, about 8 bytes per k-mer, and OPTIK_DBG_COMPACT=1 forces it for tests).
        const uint64_t efbytes = DbgGraphEF::bytes_for(n_solid, K);
        if (compact && (std::getenv("OPTIK_DBG_COMPACT") || gbytes + gbytes / 5 > limit)) {
            ef_mode = true;
            if (efbytes + efbytes / 5 > limit) {
                char msg[256];
                std::snprintf(msg, sizeof msg, "the compact graph of %zu solid k-mers needs %.2f GB "
                              "(%.2f GB + 20%% working room), more than the memory limit (%.2f GB%s)", n_solid,
                              (efbytes + efbytes / 5) / 1073741824.0, efbytes / 1073741824.0,
                              limit / 1073741824.0, explicit_limit ? ", --max-mem" : ", available");
                if (explicit_limit) throw std::runtime_error(std::string(msg) + "; raise --max-mem or -c");
                std::fprintf(stderr, "[count] warning: %s\n", msg);
            }
            std::fprintf(stderr, "[count] %zu solid k-mers: the table would need %.2f GB, so the graph is "
                         "built compact (%.2f GB, Elias-Fano keys)\n", n_solid, gbytes / 1073741824.0,
                         efbytes / 1073741824.0);
            // the solid records wait in key buckets, in memory up to what the compact graph leaves of the limit, on disk past it
            const uint64_t res = spill::resident_bytes();
            const uint64_t room = limit > res + efbytes ? limit - res - efbytes : 0;
            uint64_t per = std::max<uint64_t>(room / 2 / nthreads, uint64_t(16) << 20);
            if (const char* e = std::getenv("OPTIK_SPILL_TEST")) per = std::strtoull(e, nullptr, 10);
            ksp = std::make_unique<spill::BucketSpill<KRec>>(nthreads, std::size_t(1) << ef_bits, per, sdir);
        } else if (gbytes > limit) {
            char msg[256];
            std::snprintf(msg, sizeof msg, "the graph of %zu solid k-mers needs %.2f GB, more than the "
                          "memory limit (%.2f GB%s)", n_solid, gbytes / 1073741824.0,
                          limit / 1073741824.0, explicit_limit ? ", --max-mem" : ", available");
            if (explicit_limit) throw std::runtime_error(std::string(msg) + "; raise --max-mem or -c");
            std::fprintf(stderr, "[count] warning: %s\n", msg);
        } else {
            std::size_t cut = P;
            for (int round = 0; nd > limit && cut > 0 && round < 8; ++round) {
                // free at least the excess capacity, starting from the end
                const std::size_t hi = cut;
                uint64_t freed = 0;
                while (cut > 0 && freed < nd - limit) {
                    --cut;
                    for (unsigned t = 0; t < nthreads; ++t) freed += bufs[t][cut].bytes.capacity();
                }
                tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                    try { spill_range(t, cut, hi); } catch (const std::exception& e) { fail_with(e); }
                });
                if (failed) throw std::runtime_error(err_msg);
                spilled = true;
                malloc_trim(0);
                nd = need();
            }
            if (cut < P) pt.lap("count: spill before insert");
        }
    }
    if (ef_mode) {
        count_pass(true);
        bufs.clear();
        if (ksp->spilled()) malloc_trim(0);
        pt.lap("count: partitions (insert, key buckets)");
        // Buckets come in key order, and bucket b's nodes start after every smaller bucket's, so the buckets are loaded, sorted and written into the Elias-Fano arrays in parallel, the keys being distinct.
        auto ef = std::make_unique<DbgGraphEF>(K);
        ef->begin(n_solid);
        const std::size_t NB = std::size_t(1) << ef_bits;
        std::vector<std::size_t> base(NB + 1, 0);
        for (std::size_t b = 0; b < NB; ++b) base[b + 1] = base[b] + ksp->size(b);
        if (base[NB] != n_solid) throw std::logic_error("compact graph: bucket sizes do not add up");
        {
            std::atomic<std::size_t> nb{0};
            tpoptoa::run_on_threads(nthreads, [&](unsigned) {
                try {
                    std::vector<KRec> v;
                    for (std::size_t b; (b = nb.fetch_add(1)) < NB;) {
                        if (b + nthreads < NB) ksp->will_need(b + nthreads);
                        ksp->load(b, v);
                        std::sort(v.begin(), v.end(), [](const KRec& x, const KRec& y) {
                            const uint64_t a = x.key, c = y.key; return a < c; });
                        for (std::size_t j = 0; j < v.size(); ++j)
                            ef->put(base[b] + j, v[j].key, v[j].cnt, v[j].edge_mask);
                        std::vector<KRec>().swap(v);
                    }
                } catch (const std::exception& e) { fail_with(e); }
            });
            if (failed) throw std::runtime_error(err_msg);
        }
        ef->finish_put();
        ksp.reset();
        pt.lap("count: compact graph (Elias-Fano)");
        if (min_count > 1) ef->resolve_masks();
        *compact = std::move(ef);
        pt.lap("count: edge masks");
        return;
    }
    graph.init_table(n_solid);
    count_pass(true);
    graph.finish_insert(n_solid);     // solid k-mers are distinct
    bufs.clear();
    pt.lap("count: partitions (insert)");

    // Edges whose neighbour is solid.
    // At min_count 1 every observed edge already is one, since both k-mers of an observed (K+1)-window were counted and so are both present.
    if (min_count > 1) graph.resolve_masks();
    else graph.masks_from_raw();
    pt.lap("count: edge masks");
}

void count_kmers_with_edges(const char* fasta_path, std::size_t kk,
                            uint32_t min_count, unsigned nthreads,
                            std::vector<std::pair<uint64_t, uint32_t>>& out_kmers,
                            std::vector<uint8_t>& out_masks) {
    DbgGraph g(static_cast<unsigned>(kk));
    count_into_graph(fasta_path, kk, min_count, nthreads, g);
    g.export_sorted(out_kmers, out_masks);
}

}
}
