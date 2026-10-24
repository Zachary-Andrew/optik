#pragma once

#include "thread_budget.hpp"
#include "atomic_file.hpp"

/*
 * fm_index.hpp - FmIndex, exact k-mer membership by FM-index backward search over the separator-delimited unitigs.
 *
 * The source is FMSI (Sladky, Vesely and Brinda 2025), which indexes a masked superstring of a k-mer set, i.e. overlapping k-mers merged plus a mask bitvector, under an FM-index (Burrows and Wheeler 1994; Ferragina and Manzini 2000).
 * optik modifies it by indexing the maximal unitigs concatenated with a separator X instead, which is simpler but longer than an optimized masked superstring.
 *
 * The unitigs are maximal paths of the de Bruijn graph built from the k-mer set, so every length-k window inside a unitig is a distinct present k-mer, and the canonical forms of all unitig windows are exactly the canonical k-mer set.
 * X between unitigs means no length-k window spans a unitig boundary, so backward search of an ACGT-only query never matches a junction, and a k-mer is present iff its SA interval is non-empty.
 * There is no position array, mask or minimizer index, so the size is only the BWT's rank support, a few bits per k-mer.
 *
 * Alphabet as byte values
 *     $ = 0  unique terminator, smallest, at the very end
 *     X = 1  separator after every unitig
 *     A = 2, C = 3, G = 4, T = 5
 *
 * Only k-mers are searched, so the BWT is the k-order BWT (Schindler 1997), whose rows are text positions sorted by their first k symbols only.
 * Backward search of a pattern of length <= k is exact on it (Culpepper, Petri and Puglisi 2012), so no full suffix array is needed.
 * A window is cut at the first separator, since no ACGT pattern reaches past X.
 * A row's sort key is then one 64-bit word of packed, zero-padded bases plus the window length, and shorter windows sort ahead among equal padded words because X < A.
 * Rows are sorted in passes over ranges of the key's top kBucketBits bits to bound memory, each pass a parallel counting scatter plus small per-bucket sorts, with ties broken by text position so the BWT is deterministic.
 * The BWT is stored as four rank bitvectors, one per base, on huge pages for the random rank probes, with a cache-friendly layout after Movi (Zakeri et al. 2024).
 *
 * A query searches both orientations, because the set is canonical but unitigs are stored oriented, so a k-mer's strand may appear either way.
 * Lookups are exact, with no false positives or negatives.
 */

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <memory>
#include <vector>

#include <atomic>
#include <thread>

#include "cache_info.hpp"
#include "hugepage.hpp"
#include "kmer.hpp"
#include "wyhash.hpp"

namespace tpoptoa {

/// FmRankBits - bitvector with rank in one cache line per query.
///
/// The bits live in 64-byte blocks, each holding the ones before it and the next kBitsPerBlock = 448 bits as seven words, so rank() reads only the line that holds its position (Slotin, Algorithms for Modern Hardware, "Cache lines" and "Memory-level parallelism").
/// The division by 448 is by a constant, which the compiler turns into a multiply.
/// Since 448 = 7 x 64, word i of the flat LSB-first layout, the one the file stores, is blocks[i / 7].w[i % 7], so building, loading and saving need no copy.
struct FmRankBits {
    static constexpr unsigned kWordsPerBlock = 7;
    static constexpr uint64_t kBitsPerBlock = 64 * kWordsPerBlock;
    struct alignas(64) Block {
        uint64_t before;                   // ones in all earlier blocks
        uint64_t w[kWordsPerBlock];        // the block's bits, LSB-first
    };
    static_assert(sizeof(Block) == 64, "a block is one cache line");

    uint64_t n = 0;
    // The blocks are a slice of FmIndex's one rank store, which keeps all four bases on 2 MiB pages (hugepage.hpp) in a single allocation, so huge-page rounding is paid once.
    // A backward search does about 2k random rank probes, and with 4 KiB pages most of them also missed the TLB (bench/cpp/micro_kernels.cpp counted 32 dTLB misses per find on a 48 MB index).
    // One block past the last full one serves rank(n) when n is a multiple of kBitsPerBlock.
    Block* blocks = nullptr;
    uint64_t nblocks = 0;

    /// blocks_for() - blocks of a bitvector of @nbits bits.
    static constexpr uint64_t blocks_for(uint64_t nbits) noexcept { return nbits / kBitsPerBlock + 1; }

    void clear() noexcept { n = 0; blocks = nullptr; nblocks = 0; }

    /// init() - a bitvector of @nbits bits over @store, blocks_for(@nbits) zeroed blocks that outlive it.
    void init(uint64_t nbits, Block* store) noexcept {
        n = nbits;
        blocks = store;
        nblocks = blocks_for(nbits);
    }

    /// nwords() - words of the flat layout, ceil(n / 64).
    uint64_t nwords() const noexcept { return (n + 63) >> 6; }
    /// word() and word_ref() - word @i of the flat layout.
    uint64_t word(uint64_t i) const noexcept { return blocks[i / kWordsPerBlock].w[i % kWordsPerBlock]; }
    uint64_t& word_ref(uint64_t i) noexcept { return blocks[i / kWordsPerBlock].w[i % kWordsPerBlock]; }

    /// set() - set the bit at @pos during the build, calling build_rank() once after all sets.
    void set(uint64_t pos) noexcept { word_ref(pos >> 6) |= uint64_t(1) << (pos & 63); }

    /// build_rank() - fill every block's count of earlier ones.
    /// Called after all set() calls and by load().
    void build_rank() noexcept {
        uint64_t acc = 0;
        for (uint64_t i = 0; i < nblocks; ++i) {
            Block& b = blocks[i];
            b.before = acc;
            for (unsigned j = 0; j < kWordsPerBlock; ++j) acc += __builtin_popcountll(b.w[j]);
        }
    }

    /// block_of() - the block rank(@pos) reads, for prefetching.
    const Block* block_of(uint64_t pos) const noexcept { return &blocks[pos / kBitsPerBlock]; }

    /// rank() - number of set bits in positions [0, @pos), @pos <= n.
    ///
    /// Branch-free over the block's seven words, where words before @pos's count whole, its own word counts below @pos, and later ones not at all.
    uint64_t rank(uint64_t pos) const noexcept {
        const Block& b = blocks[pos / kBitsPerBlock];
        const unsigned off = static_cast<unsigned>(pos % kBitsPerBlock);
        const unsigned wi = off >> 6;
        const uint64_t part = (uint64_t(1) << (off & 63)) - 1;
        uint64_t res = b.before;
        for (unsigned j = 0; j < kWordsPerBlock; ++j) {
            const uint64_t m = j < wi ? ~uint64_t(0) : j == wi ? part : 0;
            res += __builtin_popcountll(b.w[j] & m);
        }
        return res;
    }
};

class FmIndex {
public:
    static constexpr uint32_t kMagic = 0x4D464B4Fu;  // file bytes "OKFM" (LE)

    static constexpr uint8_t kSymSentinel = 0;  // '$' (unique, last)
    static constexpr uint8_t kSymSep      = 1;  // 'X' separator
    static constexpr uint8_t kSymBase0    = 2;  // 'A'
    static constexpr int kNumBases = 4;
    static constexpr int kAlphabetSize = kSymBase0 + kNumBases;   // $, X, A, C, G, T
    /// Entries of the C array, one per symbol plus C[kAlphabetSize] = text length, so the rows of symbol c are [C[c], C[c + 1]).
    static constexpr int kCEntries = kAlphabetSize + 1;

    FmIndex() = default;
    FmIndex(FmIndex&&) = default;
    FmIndex& operator=(FmIndex&&) = default;
    FmIndex(const FmIndex&) = delete;
    FmIndex& operator=(const FmIndex&) = delete;

    /// FmIndex() - build from the unitigs of a k-mer set.
    ///
    /// @unitigs: every length-k window inside each unitig is a distinct present k-mer and, jointly, the windows reproduce the canonical k-mer set exactly, the de Bruijn unitig property.
    /// @n_kmers is the number of distinct canonical k-mers, kept for size() to report.
    FmIndex(unsigned kk, const std::vector<std::string>& unitigs,
            std::size_t n_kmers, unsigned nthreads = 0) {
        std::vector<uint8_t> text;
        for (const auto& u : unitigs) append_unitig(text, u.data(), u.size());
        build(kk, std::move(text), unitigs.size(), n_kmers, nthreads);
    }

    /// FmIndex() - build from a symbol text made with append_unitig() (no sentinel yet).
    FmIndex(unsigned kk, std::vector<uint8_t>&& text, std::size_t n_unitigs,
            std::size_t n_kmers, unsigned nthreads = 0) {
        build(kk, std::move(text), n_unitigs, n_kmers, nthreads);
    }

    /// append_unitig() - append one unitig and its separator X to a symbol text.
    /// Any byte other than A/C/G in either case is stored as T.
    static void append_unitig(std::vector<uint8_t>& text, const char* seq, std::size_t len) {
        const std::size_t o = text.size();
        text.resize(o + len + 1);
        for (std::size_t i = 0; i < len; ++i) {
            uint8_t base;
            switch (seq[i]) {
                case 'A': case 'a': base = 0; break;
                case 'C': case 'c': base = 1; break;
                case 'G': case 'g': base = 2; break;
                default:            base = 3; break;
            }
            text[o + i] = static_cast<uint8_t>(kSymBase0 + base);
        }
        text[o + len] = kSymSep;
    }

    unsigned     k()        const noexcept { return k_; }
    std::size_t  size()     const noexcept { return n_kmers_; }
    std::size_t  num_unitigs() const noexcept { return n_unitigs_; }
    std::size_t  super_len()  const noexcept { return super_len_; }
    std::size_t  bwt_len()    const noexcept { return bwt_len_; }

    /// find() - exact membership of k-mer @km in either orientation.
    bool find(KmerWord km) const noexcept {
        if (k_ == 0 || bwt_len_ == 0) return false;
        if (backward_search(km)) return true;
        const KmerWord r = reverse_complement(km, k_);
        return r != km && backward_search(r);
    }

    /// Searches find_batch() keeps in flight at once.
    /// Each waits on a cache miss per step, so this many overlap their misses.
    static constexpr unsigned kFindBatch = 16;

    /// batch_pays() - whether find_batch() is faster than find() one k-mer at a time, i.e. whether the rank blocks exceed the last-level cache.
    /// Inside the cache a rank costs a few cycles and the lane bookkeeping outweighs the overlap, so on hg38 chr1 (132 MB index, 32 MB L3) the batch made queries 1.1x faster than single lookups with the same layout, and on a 12 MB index 1.4x slower.
    bool batch_pays() const noexcept {
        std::size_t rank_bytes = 0;
        for (const FmRankBits& r : ranks_) rank_bytes += r.nblocks * sizeof(FmRankBits::Block);
        return rank_bytes > last_level_cache_bytes();
    }

    /// find_batch() - find() for each of @n k-mers at @km, writing 1 or 0 to
    /// @hit[i].
    ///
    /// The searches run in lockstep, kFindBatch at a time, so the rank misses of different k-mers overlap instead of each search waiting out its own (Slotin, Algorithms for Modern Hardware, "Memory-level parallelism" and "Prefetching").
    /// Every step first prefetches the two blocks each active search will read, then advances them all.
    /// A search that ends frees its lane for the next k-mer, and a k-mer absent in its given orientation is searched again as its reverse complement, so the answers equal find().
    void find_batch(const KmerWord* km, std::size_t n, uint8_t* hit) const noexcept {
        if (k_ == 0 || bwt_len_ == 0) { std::fill(hit, hit + n, uint8_t(0)); return; }
        struct Lane { KmerWord pat; uint64_t lo, hi; std::size_t idx; unsigned step; bool rc; };
        Lane lane[kFindBatch];
        unsigned active = 0;
        std::size_t next = 0;
        const unsigned K = k_;
        auto start = [&](Lane& l, std::size_t i, bool rc) {
            l.idx = i; l.rc = rc; l.pat = rc ? reverse_complement(km[i], K) : km[i];
            l.lo = 1; l.hi = bwt_len_; l.step = 0;
        };
        while (active < kFindBatch && next < n) start(lane[active++], next++, false);
        while (active) {
            for (unsigned j = 0; j < active; ++j) {
                const Lane& l = lane[j];
                const FmRankBits& r = ranks_[(l.pat >> (2u * l.step)) & 3u];
                __builtin_prefetch(r.block_of(l.lo - 1));
                __builtin_prefetch(r.block_of(l.hi));
            }
            for (unsigned j = 0; j < active;) {
                Lane& l = lane[j];
                const uint8_t sym = kSymBase0 + static_cast<uint8_t>((l.pat >> (2u * l.step)) & 3u);
                const uint64_t rl = rank_sym(sym, l.lo - 1);
                const uint64_t rh = rank_sym(sym, l.hi);
                l.lo = C_[sym] + rl + 1;
                l.hi = C_[sym] + rh;
                ++l.step;
                const bool miss = l.lo > l.hi;
                if (!miss && l.step < K) { ++j; continue; }
                if (miss && !l.rc && reverse_complement(km[l.idx], K) != km[l.idx]) {
                    start(l, l.idx, true);
                    ++j;
                    continue;
                }
                hit[l.idx] = miss ? 0 : 1;
                if (next < n) { start(l, next++, false); ++j; }
                else lane[j] = lane[--active];
            }
        }
    }

    /// backward_search() - true iff the oriented packed k-mer @km occurs in the text.
    bool backward_search(KmerWord km) const noexcept {
        const unsigned K = k_;
        uint64_t lo = 1, hi = bwt_len_;   // 1-indexed SA interval [lo, hi]
        // Packing puts the first base in the high bits and the last (3') base in the LSB, so backward search consumes bits from offset 0 upward, position K-1 first and position 0 last.
        for (unsigned i = 0; i < K; ++i) {
            const unsigned b = static_cast<unsigned>((km >> (2u * i)) & 3u);
            const uint8_t sym = kSymBase0 + static_cast<uint8_t>(b);
            const uint64_t rl = rank_sym(sym, lo - 1);
            const uint64_t rh = rank_sym(sym, hi);
            lo = C_[sym] + rl + 1;
            hi = C_[sym] + rh;
            if (lo > hi) return false;
        }
        return lo <= hi;
    }

    /// bytes() - size of the saved index file.
    std::size_t bytes() const noexcept {
        std::size_t sz = 4;                    // magic
        sz += 5 * 8;                           // k, n_kmers_, n_unitigs_, super_len_, bwt_len_
        sz += kCEntries * 8;                   // C_
        for (int b = 0; b < kNumBases; ++b) {
            sz += 8 + 8;                       // stored nbits, nwords per rank
            sz += ranks_[b].nwords() * 8;      // packed bits words
        }
        sz += 8;                               // checksum
        return sz;
    }

    double bits_per_entry() const noexcept {
        return n_kmers_ > 0 ? (bytes() * 8.0) / n_kmers_ : 0.0;
    }

    // File format, little-endian
    //
    //   [0..4)   kMagic "OKFM"
    //   [..)     u64 k, n_kmers, U_, super_len_, N_
    //   [..)     kCEntries = 7 u64 values of C_
    //   [..)     for A, C, G and T, u64 nbits, u64 nwords, nwords x u64 bits
    //   [..)     u64 checksum of every u64 above
    //
    // Load validates the structure, checking the k range, that every rank vector holds N_ bits with nothing past N_, that C_ is monotone up to N_, that each base's ones equal its C_ count, that the checksum matches and that the file is consumed exactly, so a truncated, corrupted or bit-flipped index is an error instead of out-of-bounds reads or silently wrong answers.
    /// Sum - streaming checksum over u64 words, which detects corruption only and is not cryptographic.
    /// Its constants are part of the file format.
    struct Sum {
        uint64_t h = kMagic;
        void add(uint64_t v) noexcept { h = mix64(h ^ v, 0x9E3779B97F4A7C15ULL) + 0x632BE59BD9B4E019ULL; }
    };
    /// save() - write the index to @path through atomic_write(), so a failed write, a full disk included, leaves no file and returns false.
    bool save(const char* path) const noexcept {
      try {
        return atomic_write(path, [&](std::ofstream& os) {
        Sum sum;
        auto w64 = [&](uint64_t v) { sum.add(v); os.write(reinterpret_cast<const char*>(&v), 8); };
        os.write(reinterpret_cast<const char*>(&kMagic), 4);
        w64(k_); w64(n_kmers_); w64(n_unitigs_); w64(super_len_); w64(bwt_len_);
        for (int i = 0; i < kCEntries; ++i) w64(C_[i]);
        for (int b = 0; b < kNumBases; ++b) {
            const FmRankBits& r = ranks_[b];
            w64(r.n);
            w64(r.nwords());
            for (uint64_t i = 0; i < r.nwords(); ++i) {
                const uint64_t x = r.word(i);
                sum.add(x);
                os.write(reinterpret_cast<const char*>(&x), 8);
            }
        }
        const uint64_t cs = sum.h;
        os.write(reinterpret_cast<const char*>(&cs), 8);
        });
      } catch (...) { return false; }
    }

    /// load() - read and validate an index file, returning nullptr if it is missing, malformed, truncated or fails its checksum.
    /// The caller owns the result.
    static FmIndex* load(const char* path) {
        std::ifstream is(path, std::ios::binary);
        if (!is) return nullptr;
        // Total serialized size, for exact-size validation below.
        is.seekg(0, std::ios::end);
        const std::streamoff fsize = is.tellg();
        if (fsize < 0) return nullptr;
        is.seekg(0, std::ios::beg);
        Sum sum;
        auto r64 = [&]() -> uint64_t {
            uint64_t v = 0; is.read(reinterpret_cast<char*>(&v), 8); sum.add(v); return v;
        };
        uint32_t magic;
        is.read(reinterpret_cast<char*>(&magic), 4);
        if (!is || magic != kMagic) return nullptr;
        std::unique_ptr<FmIndex> idx(new FmIndex());
        const uint64_t k = r64();
        idx->n_kmers_  = static_cast<std::size_t>(r64());
        idx->n_unitigs_        = static_cast<std::size_t>(r64());
        idx->super_len_= static_cast<std::size_t>(r64());
        idx->bwt_len_        = static_cast<std::size_t>(r64());
        for (int i = 0; i < kCEntries; ++i) idx->C_[i] = r64();
        if (!is || k < 1 || k > kMaxK) return nullptr;
        idx->k_ = static_cast<unsigned>(k);
        const uint64_t N = idx->bwt_len_;
        if (idx->C_[0] != 0 || idx->C_[kAlphabetSize] != N) return nullptr;
        for (int i = 0; i < kAlphabetSize; ++i) if (idx->C_[i] > idx->C_[i + 1]) return nullptr;
        // Each rank is stored as u64 nbits, u64 nwords and nwords u64 words.
        // The real remaining file size rejects truncated files before the payload is allocated.
        // seekg/tellg gives it accurately for files, unlike std::filebuf::in_avail, which only reflects the read-ahead buffer.
        // The four rank bitvectors must fit in what is left of the file before their one store is allocated, so a corrupt length cannot request an outsized allocation.
        if (uint64_t(kNumBases) * ((N + 63) / 64) > static_cast<uint64_t>((fsize - is.tellg()) / 8)) return nullptr;
        idx->alloc_ranks(N);
        for (int b = 0; b < kNumBases; ++b) {
            FmRankBits& r = idx->ranks_[b];
            const uint64_t rn = r64();
            const uint64_t nw = r64();
            if (!is || rn != N || nw != (N + 63) / 64) return nullptr;
            const std::streamoff remain = fsize - is.tellg();
            if (nw > static_cast<uint64_t>(remain / 8)) return nullptr;
            // The flat words go straight into the blocks, seven per block.
            for (uint64_t w = 0; w < nw; w += FmRankBits::kWordsPerBlock) {
                const uint64_t cnt = std::min<uint64_t>(FmRankBits::kWordsPerBlock, nw - w);
                uint64_t* dst = r.blocks[w / FmRankBits::kWordsPerBlock].w;
                is.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(cnt * 8));
                if (!is) return nullptr;
                for (uint64_t j = 0; j < cnt; ++j) sum.add(dst[j]);
            }
            if ((N & 63) && nw && (r.word(nw - 1) >> (N & 63))) return nullptr;   // bits past N
            r.build_rank();
            const uint64_t ones = r.rank(N);
            const uint8_t sym = static_cast<uint8_t>(kSymBase0 + b);
            if (ones != idx->C_[sym + 1] - idx->C_[sym]) return nullptr;
        }
        uint64_t cs = 0;
        is.read(reinterpret_cast<char*>(&cs), 8);
        if (!is || cs != sum.h) return nullptr;
        if (is.tellg() != fsize) return nullptr;                 // trailing bytes
        return idx.release();
    }

private:
    /// build() - build from text = u0 X u1 X ... u{U-1} X.
    /// The sentinel $ is appended here.
    void build(unsigned kk, std::vector<uint8_t>&& text, std::size_t n_unitigs,
               std::size_t n_kmers, unsigned nthreads) {
        k_ = kk;
        n_kmers_ = n_kmers;
        n_unitigs_ = n_unitigs;
        super_len_ = text.size() - n_unitigs;
        text.push_back(kSymSentinel);
        bwt_len_ = text.size();
        nthreads = resolve_threads(nthreads);
        if (nthreads == 0) nthreads = 1;
        // C_[c] = number of symbols strictly smaller than c (the BWT is a permutation of the text).
        uint64_t counts[kCEntries] = {0, 0, 0, 0, 0, 0, 0};
        for (uint8_t c : text) ++counts[c];
        uint64_t front = 0;
        for (int c = 0; c < kCEntries; ++c) { C_[c] = front; front += counts[c]; }
        alloc_ranks(bwt_len_);
        build_kbwt(text, nthreads);
        for (int b = 0; b < kNumBases; ++b) ranks_[b].build_rank();
    }

    /// Row - one row of the k-order BWT, whose key is the window's bases packed 2 bits each, cut at the first non-base and left-aligned in 2k bits with zero padding.
    /// aux = window length << kAuxLenShift | position << kAuxPosShift | BWT symbol.
    struct Row { uint64_t key, aux; };
    static constexpr unsigned kAuxLenShift = 40, kAuxPosShift = 3;
    static constexpr uint64_t kAuxSymMask = 7u;

    /// build_kbwt() - sort the rows, i.e. text positions, by key, length and position, and emit their BWT symbols straight into the rank bitvectors.
    /// Row 0 belongs to the sentinel, the smallest window.
    void build_kbwt(const std::vector<uint8_t>& text, unsigned T) {
        const std::size_t n = text.size();
        const unsigned K = k_;
        const uint64_t mask = K >= 32 ? ~uint64_t(0) : ((uint64_t(1) << (2u * K)) - 1);
        constexpr unsigned kBucketBits = 16;                      // rows bucketed by the key's top bits
        const unsigned shift = 2u * K > kBucketBits ? 2u * K - kBucketBits : 0;
        constexpr std::size_t kNbuckets = std::size_t(1) << kBucketBits;
        auto is_base = [&](std::size_t i) { return text[i] >= kSymBase0; };
        // Windows of a position range in order, as f(pos, key, len).
        auto windows = [&](std::size_t lo, std::size_t hi, auto&& f) {
            uint64_t key = 0;
            unsigned m = 0;
            bool fresh = true;
            for (std::size_t i = lo; i < hi; ++i) {
                if (!is_base(i)) { key = 0; m = 0; fresh = true; f(i, key, m); continue; }
                if (fresh) {                                   // read the window
                    key = 0; m = 0;
                    while (m < K && i + m < n && is_base(i + m)) {
                        key = (key << 2) | uint64_t(text[i + m] - kSymBase0);
                        ++m;
                    }
                    key <<= 2u * (K - m);
                    fresh = false;
                } else if (m == K && is_base(i + K - 1)) {     // roll in one base
                    key = ((key << 2) | uint64_t(text[i + K - 1] - kSymBase0)) & mask;
                } else {                                       // window shrinks
                    key = (key << 2) & mask;
                    --m;
                }
                f(i, key, m);
            }
        };
        // Per-thread histograms of the buckets over fixed position ranges (the sentinel row n - 1 is emitted first, outside the sort).
        const std::size_t body = n - 1;
        std::vector<std::size_t> cut(T + 1);
        for (unsigned t = 0; t <= T; ++t) cut[t] = body * t / T;
        std::vector<std::vector<uint64_t>> hist(T, std::vector<uint64_t>(kNbuckets, 0));
        {
            tpoptoa::run_on_threads(T, [&](unsigned t) {
                auto& h = hist[t];
                windows(cut[t], cut[t + 1], [&](std::size_t, uint64_t key, unsigned) { ++h[key >> shift]; });
            });
        }
        std::vector<uint64_t> bsize(kNbuckets, 0);
        for (unsigned t = 0; t < T; ++t)
            for (std::size_t b = 0; b < kNbuckets; ++b) bsize[b] += hist[t][b];

        uint64_t row = 0;
        auto emit = [&](uint8_t sym) {
            if (sym >= kSymBase0) ranks_[sym - kSymBase0].set(row);
            ++row;
        };
        emit(n >= 2 ? text[n - 2] : kSymSentinel);   // the sentinel's row, since $ sorts first

        // Passes over bucket ranges holding about PASS rows each, at least 2^24 rows and at most 16 passes.
        const std::size_t PASS = std::max<std::size_t>(std::size_t(1) << 24, body / 16);
        std::vector<Row> rows;
        std::vector<std::size_t> off(kNbuckets + 1), toff(T);
        for (std::size_t b0 = 0; b0 < kNbuckets;) {
            std::size_t b1 = b0, cnt = 0;
            while (b1 < kNbuckets && (b1 == b0 || cnt + bsize[b1] <= PASS)) cnt += bsize[b1++];
            if (cnt == 0) { b0 = b1; continue; }
            rows.resize(cnt);
            off[b0] = 0;
            for (std::size_t b = b0; b < b1; ++b) off[b + 1] = off[b] + bsize[b];
            // In the scatter thread t writes bucket b at off[b] + (rows of b in threads < t).
            {
                tpoptoa::run_on_threads(T, [&](unsigned t) {
                    std::vector<std::size_t> pos(b1 - b0);
                    for (std::size_t b = b0; b < b1; ++b) {
                        std::size_t o = off[b];
                        for (unsigned u = 0; u < t; ++u) o += hist[u][b];
                        pos[b - b0] = o;
                    }
                    windows(cut[t], cut[t + 1], [&](std::size_t i, uint64_t key, unsigned m) {
                        const std::size_t b = key >> shift;
                        if (b < b0 || b >= b1) return;
                        const uint8_t bw = i == 0 ? kSymSentinel : text[i - 1];
                        rows[pos[b - b0]++] = Row{key, (uint64_t(m) << kAuxLenShift) | (uint64_t(i) << kAuxPosShift) | bw};
                    });
                });
            }
            // sort every bucket, whose rows share their top bits
            {
                std::atomic<std::size_t> next{b0};
                tpoptoa::run_on_threads(T, [&](unsigned) {
                    for (std::size_t b; (b = next.fetch_add(1)) < b1;) {
                        if (bsize[b] < 2) continue;
                        std::sort(rows.begin() + off[b], rows.begin() + off[b + 1],
                                  [](const Row& x, const Row& y) {
                                      return x.key != y.key ? x.key < y.key : x.aux < y.aux;
                                  });
                    }
                });
            }
            for (const Row& r : rows) emit(static_cast<uint8_t>(r.aux & kAuxSymMask));
            b0 = b1;
        }
    }

    uint64_t rank_sym(uint8_t sym, uint64_t pos) const noexcept {
        return ranks_[sym - kSymBase0].rank(pos);
    }

    unsigned     k_ = 0;
    std::size_t  n_kmers_ = 0;
    std::size_t  n_unitigs_ = 0;
    std::size_t  super_len_ = 0;
    std::size_t  bwt_len_ = 0;
    uint64_t     C_[kCEntries] = {0, 0, 0, 0, 0, 0, 0};   // FM-index C array, counting symbols smaller than c
    FmRankBits   ranks_[kNumBases];
    HugeVec<FmRankBits::Block> rank_store_;   // the blocks of all four ranks_, one allocation

    /// alloc_ranks() - zeroed blocks for four @n-bit rank bitvectors in one allocation, handed out to ranks_.
    void alloc_ranks(uint64_t n) {
        const uint64_t nb = FmRankBits::blocks_for(n);
        rank_store_.assign(kNumBases * nb, FmRankBits::Block{});
        for (int b = 0; b < kNumBases; ++b) ranks_[b].init(n, rank_store_.data() + b * nb);
    }
};

}