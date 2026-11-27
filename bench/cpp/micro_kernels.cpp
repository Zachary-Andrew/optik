// micro_kernels.cpp - per-operation hardware counters of optik's hot kernels, measured with the unprivileged nanoBench-style harness (nanobench.hpp).
// Each kernel prints one JSON line on stdout with cycles, instructions and branch, L1D, LLC and dTLB misses per operation (median and minimum of 31 repetitions of the difference between n and 2n operations).
//
// The calibration kernels come first and check the method rather than optik.
//   nop              one NOP per op, 1 instruction and a fraction of a cycle
//   add_dependent    one dependent ADD per op, 1 instruction and about 1 cycle
//   add_independent  4 independent ADDs per op, 4 instructions and about 1 cycle
//   (the calibration kernels are straight-line blocks of 64 operations)
//   chase_l1         dependent loads in a ring half the size of the L1D, i.e. the L1 load-to-use latency
//   chase_dram       dependent loads in a random ring of at least 1 GB and 32 times the last-level cache, i.e. DRAM latency, with about 1 L1D and 1 LLC miss per op
// optik's kernels follow (k = 31), each at a cache-resident and a DRAM-resident size (_cache and _dram, where --quick runs only _cache).
// Every size comes from the cache sizes of the CPU the run is pinned to, read from sysfs, which any user may read, so the cases mean the same on any machine.
//   spss_lookup_pos / _neg   SpssCore::contains on random present or absent k-mers
//   spss_scan_read           streaming scan of 150 bp reads, per k-mer
//   mphf_lookup              minimizer MPHF rank
//   ef_select / ef_rank      Elias-Fano access and count_le
//   fm_find                  FM backward search in both orientations
//   fm_find_batch            the same k-mers through FmIndex::find_batch()
//   bpm_150                  bit-parallel semi-global edit distance, 150 bp read against a 200 bp window
//   wfa_150                  gap-affine WFA traceback on the same pair (1% errors)
//
// Usage: micro_kernels [--quick]
#include "nanobench.hpp"

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "bpm_edit.hpp"
#include "cache_info.hpp"
#include "fm_index.hpp"
#include "kmer.hpp"
#include "phobic_mphf.hpp"
#include "spss_core.hpp"
#include "wfa_affine.hpp"

using namespace tpoptoa;

static std::string rand_dna(std::mt19937_64& rng, std::size_t len) {
    static constexpr char bases[] = "ACGT";
    std::string seq(len, 'A');
    for (auto& ch : seq) ch = bases[rng() & 3];
    return seq;
}

/// The dram cases size from the cache but at most from physical memory divided by this, keeping the largest, the 32x cache ring, within an eighth of RAM.
static constexpr std::size_t kDramMemoryShare = 256;
/// Index builds are untimed and run after the pin, which their threads inherit, so one thread does the same work without queueing threads on the one CPU.
static constexpr unsigned kBuildThreads = 1;
/// Operations per straight-line asm block.
/// The asm strings repeat each instruction this many times (.rept 64), so the two numbers must match.
static constexpr uint64_t kUnroll = 64;
/// Operations per measurement call of the calibration kernels, the L1 chase and the DRAM chase, and the DRAM chase's repetitions and warm-ups.
static constexpr uint64_t kCalibrationOps = 1 << 20;
static constexpr uint64_t kChaseL1Ops = 1 << 18;
static constexpr uint64_t kChaseDramOps = 1 << 16;
static constexpr int kChaseDramReps = 11;
static constexpr int kChaseDramWarmup = 2;
/// Pointer-chase rings use one pointer per cache line.
/// The L1 ring takes 1/kL1RingShare of L1D, while the DRAM ring takes at least kMinDramRingBytes and kDramRingLlcMultiple times the cache.
static constexpr std::size_t kCacheLineBytes = 64;   // x86-64
static constexpr std::size_t kWordsPerLine = kCacheLineBytes / sizeof(uint64_t);
static constexpr std::size_t kL1RingShare = 2;
static constexpr std::size_t kMinDramRingBytes = std::size_t(1) << 30;
static constexpr std::size_t kDramRingLlcMultiple = 32;
/// "cache" sizes: a genome of 1/kCacheGenomeLlcShare of the LLC and keys of max(Mbp / kCacheKeysMbpDiv, LLC / kCacheKeysLlcDiv).
static constexpr std::size_t kCacheGenomeLlcShare = 8;
static constexpr std::size_t kCacheKeysMbpDiv = 16;
static constexpr std::size_t kCacheKeysLlcDiv = 32;
/// "dram" sizes: genome max(kDramGenomeMbp Mbp, kDramGenomeLlcMultiple LLC), keys max(kDramKeysMbp Mbp, kDramKeysLlcMultiple LLC), FM pieces max(kDramFmPiecesMin, kDramFmLlcMultiple LLC / Mbp).
static constexpr std::size_t kDramGenomeMbp = 256;
static constexpr std::size_t kDramGenomeLlcMultiple = 8;
static constexpr std::size_t kDramKeysMbp = 64;
static constexpr std::size_t kDramKeysLlcMultiple = 2;
static constexpr std::size_t kDramFmPiecesMin = 96;
static constexpr std::size_t kDramFmLlcMultiple = 3;
/// The genome is generated in pieces of this many bases.
static constexpr std::size_t kPieceBases = 1000000;
/// SPSS minimizer length for genomes above kLargeGenomeBases, and below.
static constexpr std::size_t kLargeGenomeBases = std::size_t(16) << 20;
static constexpr unsigned kLargeGenomeM = 19;
static constexpr unsigned kSmallGenomeM = 17;
/// Decimal megabytes in the size reports.
static constexpr double kBytesPerMB = 1e6;
/// Operations per measurement call of the lookups, MPHF and Elias-Fano queries, and FM searches.
static constexpr uint64_t kLookupOps = 1 << 14;
static constexpr uint64_t kMphfOps = 1 << 15;
static constexpr uint64_t kEfOps = 1 << 15;
static constexpr uint64_t kFmOps = 1 << 11;
/// The read scans use kNumReads reads of kReadLen bases, each with one substitution per kSubstOneIn bases.
/// scan_many takes kScanBatch reads per call, and each measurement call covers kScanReadsPerCall reads.
static constexpr int kNumReads = 1 << 16;
static constexpr std::size_t kReadLen = 150;
static constexpr uint64_t kSubstOneIn = 100;
static constexpr uint64_t kScanBatch = 4096;
static constexpr uint64_t kScanReadsPerCall = 512;
/// Elias-Fano test values grow by 1 to kMaxEfGap per entry.
static constexpr uint64_t kMaxEfGap = 64;
/// The alignment kernels align a read taken at kAlignReadOffset of a kAlignTextLen text.
/// That read has an edit every kAlignEditEvery bases, starting at kAlignEditOffset, and the BPM edit cap and operations per call follow.
static constexpr int kAlignTextLen = 200;
static constexpr int kAlignReadOffset = 25;
static constexpr int kAlignReadLen = 150;
static constexpr int kAlignEditEvery = 100;
static constexpr int kAlignEditOffset = 3;
static constexpr int kBpmMaxEd = 30;
static constexpr uint64_t kBpmOps = 1 << 12;
static constexpr uint64_t kWfaOps = 1 << 11;

static void emit(const nanobench::Result& result) { std::printf("%s\n", nanobench::to_json(result).c_str()); std::fflush(stdout); }

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "--quick";
    const int cpu = nanobench::pin_here();
    std::size_t l1d = cache_bytes(cpu, 1), llc = cache_bytes(cpu, 0);
    if (!l1d) l1d = kFallbackL1dBytes;
    if (!llc) llc = kFallbackLlcBytes;
    std::fprintf(stderr, "[micro] pinned to cpu %d, L1D %zu KB, last-level cache %zu MB\n", cpu,
                 l1d / 1024, llc / (1024 * 1024));
    // The dram cases grow with the cache, but the largest stays within an eighth of physical memory, so a machine with a large cache and little RAM still runs them.
    const std::size_t phys = std::size_t(sysconf(_SC_PHYS_PAGES)) * std::size_t(sysconf(_SC_PAGESIZE));
    const std::size_t dram_llc = std::min(llc, std::max(kFallbackLlcBytes, phys / kDramMemoryShare));
    std::mt19937_64 rng(1);

    // Calibration.
    // Straight-line blocks of 64 operations, with one dec and jnz per block as nanoBench unrolls, put the loop overhead at 2/64 instructions per operation, so the expected counts are exact.
    // nop is 1 instruction, a dependent add 1 instruction and 1 cycle, and four independent adds 4 instructions and about 1 cycle on four ALUs.
    // n must be a multiple of 64.
    emit(nanobench::measure("nop", [](uint64_t count) {
        uint64_t blocks = count / kUnroll;
        asm volatile("1:\n\t.rept 64\n\tnop\n\t.endr\n\tdec %0\n\tjnz 1b" : "+r"(blocks) :: "cc");
    }, kCalibrationOps));
    emit(nanobench::measure("add_dependent", [](uint64_t count) {
        uint64_t blocks = count / kUnroll, acc = 1;
        asm volatile("1:\n\t.rept 64\n\tadd %1, %1\n\t.endr\n\tdec %0\n\tjnz 1b" : "+r"(blocks), "+r"(acc) :: "cc");
        nanobench::keep(acc);
    }, kCalibrationOps));
    emit(nanobench::measure("add_independent", [](uint64_t count) {
        uint64_t blocks = count / kUnroll, acc_a = 1, acc_b = 2, acc_c = 3, acc_d = 4;
        asm volatile("1:\n\t.rept 64\n\tadd %1, %1\n\tadd %2, %2\n\tadd %3, %3\n\tadd %4, %4\n\t.endr\n\t"
                     "dec %0\n\tjnz 1b" : "+r"(blocks), "+r"(acc_a), "+r"(acc_b), "+r"(acc_c), "+r"(acc_d) :: "cc");
        nanobench::keep(acc_a + acc_b + acc_c + acc_d);
    }, kCalibrationOps));
    auto ring = [&](std::size_t bytes) {
        const std::size_t num_lines = bytes / kCacheLineBytes;
        std::vector<uint64_t> perm(num_lines);
        for (std::size_t i = 0; i < num_lines; ++i) perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), rng);
        std::vector<uint64_t> next(num_lines * kWordsPerLine, 0);              // one pointer per 64-byte line
        for (std::size_t i = 0; i < num_lines; ++i) next[perm[i] * kWordsPerLine] = perm[(i + 1) % num_lines] * kWordsPerLine;
        return next;
    };
    {
        auto l1_ring = ring(l1d / kL1RingShare);
        const uint64_t* base1 = l1_ring.data();
        emit(nanobench::measure("chase_l1", [&](uint64_t count) {
            uint64_t blocks = count / kUnroll, cursor = 0;
            asm volatile("1:\n\t.rept 64\n\tmov (%2,%1,8), %1\n\t.endr\n\tdec %0\n\tjnz 1b"
                         : "+r"(blocks), "+r"(cursor) : "r"(base1) : "cc", "memory");
            nanobench::keep(cursor);
        }, kChaseL1Ops));
    }
    if (!quick) {
        // The chase continues where the previous call stopped, so every load touches a line not seen in this measurement.
        // Restarting at the ring's head made all repetitions walk the same 4 MB, which a 32 MB L3 then held, giving 62 cycles instead of DRAM latency.
        auto dram_ring = ring(std::max(kMinDramRingBytes, kDramRingLlcMultiple * dram_llc));
        const uint64_t* base2 = dram_ring.data();
        uint64_t dram_cursor = 0;
        emit(nanobench::measure("chase_dram", [&](uint64_t count) {
            uint64_t blocks = count / kUnroll, cursor = dram_cursor;
            asm volatile("1:\n\t.rept 64\n\tmov (%2,%1,8), %1\n\t.endr\n\tdec %0\n\tjnz 1b"
                         : "+r"(blocks), "+r"(cursor) : "r"(base2) : "cc", "memory");
            dram_cursor = cursor;
        }, kChaseDramOps, kChaseDramReps, kChaseDramWarmup));
    }

    // optik kernels.
    // Each structure is measured at two sizes, "cache", an eighth of the last-level cache so it stays resident, and "dram", well past it as real indexes are, like the chr1 SPSS at about 200 MB.
    // On a 32 MB L3 these are a 4 Mbp genome with 1M keys and a 256 Mbp genome with 64M keys.
    // Query positions continue across calls, so repetitions do not touch the same lines again.
    constexpr unsigned kmer_len = 31;
    const uint64_t mask = kmer_mask(kmer_len);
    // fm_pieces is how many 1 Mbp pieces the FM index takes, all of the cache genome and three times the last-level cache for dram.
    struct Size { const char* tag; std::size_t genome, keys, fm_pieces; };
    // The genome is generated in 1 Mbp pieces, so the cache case keeps at least one.
    constexpr std::size_t mbp = std::size_t(1) << 20;
    const std::size_t cache_genome = std::max(mbp, llc / kCacheGenomeLlcShare);
    std::vector<Size> sizes = {{"cache", cache_genome, std::max(mbp / kCacheKeysMbpDiv, llc / kCacheKeysLlcDiv), cache_genome / mbp}};
    if (!quick)
        sizes.push_back({"dram", std::max(kDramGenomeMbp * mbp, kDramGenomeLlcMultiple * dram_llc), std::max(kDramKeysMbp * mbp, kDramKeysLlcMultiple * dram_llc),
                         std::max<std::size_t>(kDramFmPiecesMin, kDramFmLlcMultiple * dram_llc / mbp)});
    for (const Size& size_cfg : sizes)
        std::fprintf(stderr, "[micro] %s: %zu Mbp genome, %zu keys, FM over %zu Mbp\n", size_cfg.tag,
                     size_cfg.genome / mbp, size_cfg.keys, size_cfg.fm_pieces);
    for (const Size& size_cfg : sizes) {
        const std::string suffix = std::string("_") + size_cfg.tag;
        std::vector<std::string> seqs;
        for (std::size_t done = 0; done < size_cfg.genome; done += kPieceBases) seqs.push_back(rand_dna(rng, kPieceBases));
        auto pack = [&](const std::string& seq, std::size_t start) {
            uint64_t packed = 0;
            for (unsigned i = 0; i < kmer_len; ++i) packed = (packed << 2) | base_encode(seq[start + i]);
            return packed & mask;
        };
        constexpr std::size_t num_queries = 1u << 22;                 // more queries than any run uses
        std::vector<uint64_t> pos(num_queries), neg(num_queries);
        for (auto& query : pos) { const auto& seq = seqs[rng() % seqs.size()]; query = pack(seq, rng() % (seq.size() - kmer_len)); }
        for (auto& query : neg) query = rng() & mask;
        {
            SpssCore idx;
            idx.build(seqs, kmer_len, size_cfg.genome > kLargeGenomeBases ? kLargeGenomeM : kSmallGenomeM, kBuildThreads);
            std::fprintf(stderr, "[micro] spss%s: %.1f MB\n", suffix.c_str(), idx.bytes() / kBytesPerMB);
            uint64_t query_pos = 0;
            emit(nanobench::measure("spss_lookup_pos" + suffix, [&](uint64_t count) {
                uint64_t hits = 0;
                for (uint64_t i = 0; i < count; ++i) hits += idx.contains(pos[query_pos++ % num_queries]);
                nanobench::keep(hits);
            }, kLookupOps));
            emit(nanobench::measure("spss_lookup_neg" + suffix, [&](uint64_t count) {
                uint64_t hits = 0;
                for (uint64_t i = 0; i < count; ++i) hits += idx.contains(neg[query_pos++ % num_queries]);
                nanobench::keep(hits);
            }, kLookupOps));
            // 150 bp reads with 1% substitutions through scan_many, as query batches them
            std::vector<std::string> reads;
            for (int i = 0; i < kNumReads; ++i) {
                const auto& seq = seqs[rng() % seqs.size()];
                std::string read = seq.substr(rng() % (seq.size() - kReadLen), kReadLen);
                for (auto& ch : read) if (rng() % kSubstOneIn == 0) ch = "ACGT"[rng() & 3];
                reads.push_back(read);
            }
            std::vector<const char*> ptr; std::vector<std::size_t> len;
            for (auto& read : reads) { ptr.push_back(read.data()); len.push_back(read.size()); }
            const uint64_t per_read = kReadLen - kmer_len + 1;
            std::size_t read_pos = 0;
            emit(nanobench::measure("spss_scan_read" + suffix, [&](uint64_t count) {
                SpssCore::ScanStats stats;
                uint64_t left = count / per_read;
                while (left) {
                    const std::size_t batch = std::min<uint64_t>({left, kScanBatch, uint64_t(reads.size() - read_pos)});
                    idx.scan_many(ptr.data() + read_pos, len.data() + read_pos, batch, stats);
                    left -= batch; read_pos = (read_pos + batch) % reads.size();
                }
                nanobench::keep(stats.positive);
            }, per_read * kScanReadsPerCall));
        }
        {
            std::vector<uint64_t> keys(size_cfg.keys);
            for (auto& key : keys) key = rng();
            std::sort(keys.begin(), keys.end());
            keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
            Mphf mphf;
            mphf.build(keys.data(), keys.size(), Mphf::kDefaultBucketFactor, Mphf::kDefaultAlpha, kBuildThreads);
            std::vector<uint64_t> probe(num_queries);
            for (auto& key : probe) key = keys[rng() % keys.size()];
            uint64_t key_pos = 0;
            emit(nanobench::measure("mphf_lookup" + suffix, [&](uint64_t count) {
                uint64_t hits = 0;
                for (uint64_t i = 0; i < count; ++i) hits += mphf(probe[key_pos++ % num_queries]);
                nanobench::keep(hits);
            }, kMphfOps));
            std::vector<uint64_t> values(size_cfg.keys);
            uint64_t acc = 0;
            for (auto& value : values) { acc += 1 + rng() % kMaxEfGap; value = acc; }
            EfSeq ef_seq;
            ef_seq.build(values, acc + 1);
            std::vector<uint64_t> positions(num_queries), val(num_queries);
            for (auto& pos_val : positions) pos_val = rng() % values.size();
            for (auto& query_val : val) query_val = rng() % acc;
            uint64_t ef_pos = 0;
            emit(nanobench::measure("ef_select" + suffix, [&](uint64_t count) {
                uint64_t total = 0;
                for (uint64_t i = 0; i < count; ++i) total += ef_seq[positions[ef_pos++ % num_queries]];
                nanobench::keep(total);
            }, kEfOps));
            emit(nanobench::measure("ef_rank" + suffix, [&](uint64_t count) {
                uint64_t total = 0;
                for (uint64_t i = 0; i < count; ++i) total += ef_seq.count_le(val[ef_pos++ % num_queries]);
                nanobench::keep(total);
            }, kEfOps));
        }
        {
            std::vector<std::string> unitigs(seqs.begin(), seqs.begin() + std::min(seqs.size(), size_cfg.fm_pieces));
            FmIndex fm_idx(kmer_len, unitigs, 0, kBuildThreads);
            std::fprintf(stderr, "[micro] fm%s: %.1f MB\n", suffix.c_str(), fm_idx.bytes() / kBytesPerMB);
            std::vector<uint64_t> queries(num_queries);
            for (auto& query : queries) { const auto& unitig = unitigs[rng() % unitigs.size()]; query = pack(unitig, rng() % (unitig.size() - kmer_len)); }
            uint64_t fm_pos = 0;
            emit(nanobench::measure("fm_find" + suffix, [&](uint64_t count) {
                uint64_t hits = 0;
                for (uint64_t i = 0; i < count; ++i) hits += fm_idx.find(queries[fm_pos++ % num_queries]);
                nanobench::keep(hits);
            }, kFmOps));
            // The same lookups through find_batch(), whose searches overlap their misses.
            std::vector<uint8_t> hits(num_queries);
            uint64_t batch_pos = 0;
            emit(nanobench::measure("fm_find_batch" + suffix, [&](uint64_t count) {
                if (batch_pos + count > num_queries) batch_pos = 0;
                fm_idx.find_batch(queries.data() + batch_pos, count, hits.data() + batch_pos);
                batch_pos += count;
                nanobench::keep(hits[batch_pos - 1]);
            }, kFmOps));
        }
    }
    {
        // The alignment kernels are compute-bound, since the pair stays in L1.
        const std::string text = rand_dna(rng, kAlignTextLen);
        std::string read = text.substr(kAlignReadOffset, kAlignReadLen);
        for (int i = 0; i < kAlignReadLen; i += kAlignEditEvery) read[i + kAlignEditOffset] = read[i + kAlignEditOffset] == 'A' ? 'C' : 'A';
        emit(nanobench::measure("bpm_150", [&](uint64_t count) {
            int total = 0;
            for (uint64_t i = 0; i < count; ++i) total += bpm_semiglobal_ed(read.data(), kAlignReadLen, text.data(), kAlignTextLen, kBpmMaxEd);
            nanobench::keep(total);
        }, kBpmOps));
        AffinePenalties pen;   // the default 4 mismatch, 6 gap open, 2 gap extend
        WfaAffine aligner(pen);
        emit(nanobench::measure("wfa_150", [&](uint64_t count) {
            int total = 0;
            for (uint64_t i = 0; i < count; ++i) total += aligner.align(read.data(), kAlignReadLen, text.data(), kAlignTextLen, kAlignTextLen, kAlignTextLen).score;
            nanobench::keep(total);
        }, kWfaOps));
    }
    return 0;
}
