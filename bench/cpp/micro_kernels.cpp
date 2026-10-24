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

static std::string rand_dna(std::mt19937_64& r, std::size_t n) {
    static constexpr char B[] = "ACGT";
    std::string s(n, 'A');
    for (auto& c : s) c = B[r() & 3];
    return s;
}

/// The dram cases size from the cache but at most from physical memory divided by this, keeping the largest, the 32x cache ring, within an eighth of RAM.
static constexpr std::size_t kDramMemoryShare = 256;
/// Index builds are untimed and run after the pin, which their threads inherit, so one thread does the same work without queueing threads on the one CPU.
static constexpr unsigned kBuildThreads = 1;

static void emit(const nb::Result& r) { std::printf("%s\n", nb::to_json(r).c_str()); std::fflush(stdout); }

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "--quick";
    const int cpu = nb::pin_here();
    std::size_t l1d = cache_bytes(cpu, 1), llc = cache_bytes(cpu, 0);
    if (!l1d) l1d = kFallbackL1dBytes;
    if (!llc) llc = kFallbackLlcBytes;
    std::fprintf(stderr, "[micro] pinned to cpu %d, L1D %zu KB, last-level cache %zu MB\n", cpu,
                 l1d >> 10, llc >> 20);
    // The dram cases grow with the cache, but the largest stays within an eighth of physical memory, so a machine with a large cache and little RAM still runs them.
    const std::size_t phys = std::size_t(sysconf(_SC_PHYS_PAGES)) * std::size_t(sysconf(_SC_PAGESIZE));
    const std::size_t dram_llc = std::min(llc, std::max(kFallbackLlcBytes, phys / kDramMemoryShare));
    std::mt19937_64 rng(1);

    // Calibration.
    // Straight-line blocks of 64 operations, with one dec and jnz per block as nanoBench unrolls, put the loop overhead at 2/64 instructions per operation, so the expected counts are exact.
    // nop is 1 instruction, a dependent add 1 instruction and 1 cycle, and four independent adds 4 instructions and about 1 cycle on four ALUs.
    // n must be a multiple of 64.
    emit(nb::measure("nop", [](uint64_t n) {
        uint64_t blocks = n / 64;
        asm volatile("1:\n\t.rept 64\n\tnop\n\t.endr\n\tdec %0\n\tjnz 1b" : "+r"(blocks) :: "cc");
    }, 1 << 20));
    emit(nb::measure("add_dependent", [](uint64_t n) {
        uint64_t blocks = n / 64, a = 1;
        asm volatile("1:\n\t.rept 64\n\tadd %1, %1\n\t.endr\n\tdec %0\n\tjnz 1b" : "+r"(blocks), "+r"(a) :: "cc");
        nb::keep(a);
    }, 1 << 20));
    emit(nb::measure("add_independent", [](uint64_t n) {
        uint64_t blocks = n / 64, a = 1, b = 2, c = 3, d = 4;
        asm volatile("1:\n\t.rept 64\n\tadd %1, %1\n\tadd %2, %2\n\tadd %3, %3\n\tadd %4, %4\n\t.endr\n\t"
                     "dec %0\n\tjnz 1b" : "+r"(blocks), "+r"(a), "+r"(b), "+r"(c), "+r"(d) :: "cc");
        nb::keep(a + b + c + d);
    }, 1 << 20));
    auto ring = [&](std::size_t bytes) {
        const std::size_t n = bytes / 64;
        std::vector<uint64_t> perm(n);
        for (std::size_t i = 0; i < n; ++i) perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), rng);
        std::vector<uint64_t> next(n * 8, 0);              // one pointer per 64-byte line
        for (std::size_t i = 0; i < n; ++i) next[perm[i] * 8] = perm[(i + 1) % n] * 8;
        return next;
    };
    {
        auto r1 = ring(l1d / 2);
        const uint64_t* base1 = r1.data();
        emit(nb::measure("chase_l1", [&](uint64_t n) {
            uint64_t blocks = n / 64, p = 0;
            asm volatile("1:\n\t.rept 64\n\tmov (%2,%1,8), %1\n\t.endr\n\tdec %0\n\tjnz 1b"
                         : "+r"(blocks), "+r"(p) : "r"(base1) : "cc", "memory");
            nb::keep(p);
        }, 1 << 18));
    }
    if (!quick) {
        // The chase continues where the previous call stopped, so every load touches a line not seen in this measurement.
        // Restarting at the ring's head made all repetitions walk the same 4 MB, which a 32 MB L3 then held, giving 62 cycles instead of DRAM latency.
        auto r2 = ring(std::max(std::size_t(1) << 30, 32 * dram_llc));
        const uint64_t* base2 = r2.data();
        uint64_t p2 = 0;
        emit(nb::measure("chase_dram", [&](uint64_t n) {
            uint64_t blocks = n / 64, p = p2;
            asm volatile("1:\n\t.rept 64\n\tmov (%2,%1,8), %1\n\t.endr\n\tdec %0\n\tjnz 1b"
                         : "+r"(blocks), "+r"(p) : "r"(base2) : "cc", "memory");
            p2 = p;
        }, 1 << 16, 11, 2));
    }

    // optik kernels.
    // Each structure is measured at two sizes, "cache", an eighth of the last-level cache so it stays resident, and "dram", well past it as real indexes are, like the chr1 SPSS at about 200 MB.
    // On a 32 MB L3 these are a 4 Mbp genome with 1M keys and a 256 Mbp genome with 64M keys.
    // Query positions continue across calls, so repetitions do not touch the same lines again.
    constexpr unsigned K = 31;
    const uint64_t mask = kmer_mask(K);
    // fm_pieces is how many 1 Mbp pieces the FM index takes, all of the cache genome and three times the last-level cache for dram.
    struct Size { const char* tag; std::size_t genome, keys, fm_pieces; };
    // The genome is generated in 1 Mbp pieces, so the cache case keeps at least one.
    constexpr std::size_t mbp = std::size_t(1) << 20;
    const std::size_t cache_genome = std::max(mbp, llc / 8);
    std::vector<Size> sizes = {{"cache", cache_genome, std::max(mbp / 16, llc / 32), cache_genome / mbp}};
    if (!quick)
        sizes.push_back({"dram", std::max(256 * mbp, 8 * dram_llc), std::max(64 * mbp, 2 * dram_llc),
                         std::max<std::size_t>(96, 3 * dram_llc / mbp)});
    for (const Size& sz : sizes)
        std::fprintf(stderr, "[micro] %s: %zu Mbp genome, %zu keys, FM over %zu Mbp\n", sz.tag,
                     sz.genome >> 20, sz.keys, sz.fm_pieces);
    for (const Size& sz : sizes) {
        const std::string T = std::string("_") + sz.tag;
        std::vector<std::string> seqs;
        for (std::size_t done = 0; done < sz.genome; done += 1000000) seqs.push_back(rand_dna(rng, 1000000));
        auto pack = [&](const std::string& s, std::size_t at) {
            uint64_t v = 0;
            for (unsigned i = 0; i < K; ++i) v = (v << 2) | base_encode(s[at + i]);
            return v & mask;
        };
        constexpr std::size_t NQ = 1u << 22;                 // more queries than any run uses
        std::vector<uint64_t> pos(NQ), neg(NQ);
        for (auto& x : pos) { const auto& s = seqs[rng() % seqs.size()]; x = pack(s, rng() % (s.size() - K)); }
        for (auto& x : neg) x = rng() & mask;
        {
            SpssCore idx;
            idx.build(seqs, K, sz.genome > (16u << 20) ? 19 : 17, kBuildThreads);
            std::fprintf(stderr, "[micro] spss%s: %.1f MB\n", T.c_str(), idx.bytes() / 1e6);
            uint64_t qp = 0;
            emit(nb::measure("spss_lookup_pos" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += idx.contains(pos[qp++ % NQ]);
                nb::keep(h);
            }, 1 << 14));
            emit(nb::measure("spss_lookup_neg" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += idx.contains(neg[qp++ % NQ]);
                nb::keep(h);
            }, 1 << 14));
            // 150 bp reads with 1% substitutions through scan_many, as query batches them
            std::vector<std::string> reads;
            for (int i = 0; i < (1 << 16); ++i) {
                const auto& s = seqs[rng() % seqs.size()];
                std::string q = s.substr(rng() % (s.size() - 150), 150);
                for (auto& c : q) if (rng() % 100 == 0) c = "ACGT"[rng() & 3];
                reads.push_back(q);
            }
            std::vector<const char*> ptr; std::vector<std::size_t> len;
            for (auto& r : reads) { ptr.push_back(r.data()); len.push_back(r.size()); }
            const uint64_t per_read = 150 - K + 1;
            std::size_t rp = 0;
            emit(nb::measure("spss_scan_read" + T, [&](uint64_t n) {
                SpssCore::ScanStats st;
                uint64_t left = n / per_read;
                while (left) {
                    const std::size_t b = std::min<uint64_t>({left, uint64_t(4096), uint64_t(reads.size() - rp)});
                    idx.scan_many(ptr.data() + rp, len.data() + rp, b, st);
                    left -= b; rp = (rp + b) % reads.size();
                }
                nb::keep(st.positive);
            }, per_read * 512));
        }
        {
            std::vector<uint64_t> keys(sz.keys);
            for (auto& x : keys) x = rng();
            std::sort(keys.begin(), keys.end());
            keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
            Mphf f;
            f.build(keys.data(), keys.size(), 6.0, 0.94, kBuildThreads);
            std::vector<uint64_t> probe(NQ);
            for (auto& x : probe) x = keys[rng() % keys.size()];
            uint64_t kp = 0;
            emit(nb::measure("mphf_lookup" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += f(probe[kp++ % NQ]);
                nb::keep(h);
            }, 1 << 15));
            std::vector<uint64_t> v(sz.keys);
            uint64_t acc = 0;
            for (auto& x : v) { acc += 1 + rng() % 64; x = acc; }
            EfSeq ef;
            ef.build(v, acc + 1);
            std::vector<uint64_t> at(NQ), val(NQ);
            for (auto& x : at) x = rng() % v.size();
            for (auto& x : val) x = rng() % acc;
            uint64_t ep = 0;
            emit(nb::measure("ef_select" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += ef[at[ep++ % NQ]];
                nb::keep(h);
            }, 1 << 15));
            emit(nb::measure("ef_rank" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += ef.count_le(val[ep++ % NQ]);
                nb::keep(h);
            }, 1 << 15));
        }
        {
            std::vector<std::string> unitigs(seqs.begin(), seqs.begin() + std::min(seqs.size(), sz.fm_pieces));
            FmIndex fm(K, unitigs, 0, kBuildThreads);
            std::fprintf(stderr, "[micro] fm%s: %.1f MB\n", T.c_str(), fm.bytes() / 1e6);
            std::vector<uint64_t> q(NQ);
            for (auto& x : q) { const auto& s = unitigs[rng() % unitigs.size()]; x = pack(s, rng() % (s.size() - K)); }
            uint64_t fp = 0;
            emit(nb::measure("fm_find" + T, [&](uint64_t n) {
                uint64_t h = 0;
                for (uint64_t i = 0; i < n; ++i) h += fm.find(q[fp++ % NQ]);
                nb::keep(h);
            }, 1 << 11));
            // The same lookups through find_batch(), whose searches overlap their misses.
            std::vector<uint8_t> hits(NQ);
            uint64_t bp = 0;
            emit(nb::measure("fm_find_batch" + T, [&](uint64_t n) {
                if (bp + n > NQ) bp = 0;
                fm.find_batch(q.data() + bp, n, hits.data() + bp);
                bp += n;
                nb::keep(hits[bp - 1]);
            }, 1 << 11));
        }
    }
    {
        // The alignment kernels are compute-bound, since the pair stays in L1.
        const std::string t = rand_dna(rng, 200);
        std::string qr = t.substr(25, 150);
        for (int i = 0; i < 150; i += 100) qr[i + 3] = qr[i + 3] == 'A' ? 'C' : 'A';
        emit(nb::measure("bpm_150", [&](uint64_t n) {
            int h = 0;
            for (uint64_t i = 0; i < n; ++i) h += bpm_semiglobal_ed(qr.data(), 150, t.data(), 200, 30);
            nb::keep(h);
        }, 1 << 12));
        AffinePenalties pen; pen.x = 4; pen.o = 6; pen.e = 2;
        WfaAffine w(pen);
        emit(nb::measure("wfa_150", [&](uint64_t n) {
            int h = 0;
            for (uint64_t i = 0; i < n; ++i) h += w.align(qr.data(), 150, t.data(), 200, 200, 200).score;
            nb::keep(h);
        }, 1 << 11));
    }
    return 0;
}
