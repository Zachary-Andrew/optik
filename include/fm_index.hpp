#pragma once

#include "bits.hpp"
#include "thread_budget.hpp"
#include "atomic_file.hpp"

/**
 * @file
 * @brief FmIndex, exact k-mer membership by FM-index backward search over separator-delimited unitigs.
 *
 * The idea comes from FMSI (Sladky, Vesely and Brinda 2025), which puts a masked superstring of a k-mer set under an FM-index.
 * A masked superstring merges overlapping k-mers and adds a mask bitvector (FM-index: Burrows and Wheeler 1994, Ferragina and Manzini 2000).
 * We index the maximal unitigs joined by a separator X instead, which is simpler, though longer than an optimised masked superstring.
 *
 * Unitigs are maximal de Bruijn paths, so every length-k window inside one is a distinct k-mer from the set.
 * Together, the canonical forms of all unitig windows are exactly the canonical k-mer set.
 * The X between unitigs keeps any length-k window from spanning a boundary, so an ACGT-only query never matches a junction.
 * A k-mer is therefore present iff its SA interval is non-empty.
 * With no position array, mask or minimizer index, the size is just the BWT's rank support, a few bits per k-mer.
 *
 * Alphabet as byte values
 *     $ = 0  unique terminator, smallest, at the very end
 *     X = 1  separator after every unitig
 *     A = 2, C = 3, G = 4, T = 5
 *
 * Since only k-mers are ever searched, the BWT is the k-order BWT (Schindler 1997), whose rows are sorted by their first k symbols only.
 * Backward search for patterns of length <= k is still exact on it (Culpepper, Petri and Puglisi 2012), so no full suffix array is needed.
 * Each window is cut at the first separator, because no ACGT pattern reaches past X.
 * A row's sort key then fits in one 64-bit word of packed, zero-padded bases plus the window length.
 * Among equal padded words, shorter windows sort first, which matches X < A.
 * To bound memory, rows are sorted in passes over ranges of the key's top kBucketBits bits.
 * Each pass is a parallel counting scatter followed by small per-bucket sorts, with ties broken by text position so the BWT is deterministic.
 * The BWT is stored as four rank bitvectors, one per base, on huge pages for the random rank probes and with a cache-friendly layout after Movi (Zakeri et al. 2024).
 *
 * Queries search both orientations, because the set is canonical while unitigs are stored oriented.
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

/// Bitvector whose rank() touches one cache line.
///
/// Each 64-byte block holds the count of earlier ones plus the next kBitsPerBlock bits as seven words.
/// rank() therefore reads only the line containing its position (Slotin, "Cache lines" and "Memory-level parallelism").
/// Dividing by kBitsPerBlock is a division by a constant, which the compiler turns into a multiply.
/// Because kBitsPerBlock is a whole number of words, word i of the flat file layout is blocks[i / 7].words[i % 7].
/// Building, loading and saving thus need no copy.
struct FmRankBits {
    static constexpr std::size_t kCacheLineBytes = 64;   // x86-64, and one Block fills exactly one
    static constexpr unsigned kWordsPerBlock = 7;
    static constexpr uint64_t kBitsPerBlock = kWordBits * kWordsPerBlock;

    struct alignas(kCacheLineBytes) Block {
        uint64_t before;                   // ones in all earlier blocks
        uint64_t words[kWordsPerBlock];    // LSB-first
    };
    static_assert(sizeof(Block) == kCacheLineBytes, "a block is one cache line");

    uint64_t bit_count = 0;
    // The blocks are a slice of FmIndex's single rank store, which keeps all four bases on 2 MiB pages, so huge-page rounding is paid once.
    // That matters because a backward search makes about 2k random rank probes.
    // With 4 KiB pages most of them also missed the TLB, and bench/cpp/micro_kernels.cpp counted 32 dTLB misses per find on a 48 MB index.
    // One extra block past the last full one serves rank(n) when n is a multiple of kBitsPerBlock.
    Block* blocks = nullptr;
    uint64_t nblocks = 0;

    static constexpr uint64_t blocks_for(uint64_t nbits) noexcept { return nbits / kBitsPerBlock + 1; }

    void clear() noexcept { bit_count = 0; blocks = nullptr; nblocks = 0; }

    /// @p store must hold blocks_for(@p nbits) zeroed blocks and outlive this bitvector.
    void init(uint64_t nbits, Block* store) noexcept {
        bit_count = nbits;
        blocks = store;
        nblocks = blocks_for(nbits);
    }

    uint64_t nwords() const noexcept { return words_for_bits(bit_count); }
    uint64_t word(uint64_t word_idx) const noexcept { return blocks[word_idx / kWordsPerBlock].words[word_idx % kWordsPerBlock]; }
    uint64_t& word_ref(uint64_t word_idx) noexcept { return blocks[word_idx / kWordsPerBlock].words[word_idx % kWordsPerBlock]; }

    /// Build-time only, and build_rank() must follow once all bits are set.
    void set(uint64_t pos) noexcept { word_ref(word_index(pos)) |= uint64_t(1) << bit_offset(pos); }

    /// Fills in every block's count of earlier ones, after the last set() and from load().
    void build_rank() noexcept {
        uint64_t acc = 0;

        for (uint64_t i = 0; i < nblocks; ++i) {
            Block& block = blocks[i];

            block.before = acc;

            for (unsigned j = 0; j < kWordsPerBlock; ++j) acc += __builtin_popcountll(block.words[j]);
        }
    }

    /// Exposed so callers can prefetch the block rank(@p pos) will read.
    const Block* block_of(uint64_t pos) const noexcept { return &blocks[pos / kBitsPerBlock]; }

    /// Counts set bits in [0, @p pos), for @p pos <= n.
    ///
    /// It is branch-free over the block's seven words.
    /// Words before @p pos count in full, @p pos's own word counts only the bits below it, and later words count not at all.
    uint64_t rank(uint64_t pos) const noexcept {
        const Block& block = blocks[pos / kBitsPerBlock];
        const unsigned off = static_cast<unsigned>(pos % kBitsPerBlock);
        const unsigned word_in_block = word_index(off);
        const uint64_t part = (uint64_t(1) << bit_offset(off)) - 1;
        uint64_t res = block.before;

        for (unsigned j = 0; j < kWordsPerBlock; ++j) {
            const uint64_t mask = j < word_in_block ? ~uint64_t(0) : j == word_in_block ? part : 0;

            res += __builtin_popcountll(block.words[j] & mask);
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
    static constexpr int kNumBases = static_cast<int>(tpoptoa::kNumBases);
    static constexpr int kAlphabetSize = kSymBase0 + kNumBases;   // $, X, A, C, G, T
    /// One C-array entry per symbol plus C[kAlphabetSize] = text length, so symbol c owns rows [C[c], C[c + 1]).
    static constexpr int kCEntries = kAlphabetSize + 1;
    /// Field counts of the file layout described above save().
    static constexpr int kMagicBytes = static_cast<int>(sizeof(kMagic));
    static constexpr int kSerialWordBytes = static_cast<int>(sizeof(uint64_t));
    static constexpr int kHeaderWords = 5;
    static constexpr int kRankHeaderWords = 2;
    /// Checksum update is state = mix64(state ^ word, kGoldenRatio64) + kChecksumAdd, and the constant is part of the file format.
    static constexpr uint64_t kChecksumAdd = 0x632BE59BD9B4E019ULL;
    /// Construction passes hold about max(kMinPassRows, rows / kMaxPasses) rows each.
    static constexpr std::size_t kMinPassRows = std::size_t(1) << 24;
    static constexpr std::size_t kMaxPasses = 16;

    FmIndex() = default;
    FmIndex(FmIndex&&) = default;
    FmIndex& operator=(FmIndex&&) = default;
    FmIndex(const FmIndex&) = delete;
    FmIndex& operator=(const FmIndex&) = delete;

    /// @p unitigs must have the de Bruijn unitig property: every window is a distinct k-mer, and together they reproduce the canonical set exactly.
    /// @p n_kmers only feeds size().
    FmIndex(unsigned kmer_len, const std::vector<std::string>& unitigs,
            std::size_t n_kmers, unsigned nthreads = 0) {
        std::vector<uint8_t> text;

        for (const auto& unitig : unitigs) append_unitig(text, unitig.data(), unitig.size());

        build(kmer_len, std::move(text), unitigs.size(), n_kmers, nthreads);
    }

    /// Takes a symbol text built with append_unitig(), without the sentinel.
    FmIndex(unsigned kmer_len, std::vector<uint8_t>&& text, std::size_t n_unitigs,
            std::size_t n_kmers, unsigned nthreads = 0) {
        build(kmer_len, std::move(text), n_unitigs, n_kmers, nthreads);
    }

    /// Appends a unitig followed by its separator X.
    /// Bytes other than A/C/G (in either case) are stored as T.
    static void append_unitig(std::vector<uint8_t>& text, const char* seq, std::size_t len) {
        const std::size_t start = text.size();

        text.resize(start + len + 1);

        for (std::size_t i = 0; i < len; ++i) {
            uint8_t base;

            switch (seq[i]) {
                case 'A': case 'a': base = kCodeA; break;
                case 'C': case 'c': base = kCodeC; break;
                case 'G': case 'g': base = kCodeG; break;
                default:            base = kCodeT; break;
            }

            // Bases sit above the separator and sentinel in symbol order.
            text[start + i] = static_cast<uint8_t>(kSymBase0 + base);
        }

        // The separator keeps any k-mer window from spanning two unitigs.
        text[start + len] = kSymSep;
    }

    unsigned     kmer_len()        const noexcept { return kmer_len_; }
    std::size_t  size()     const noexcept { return n_kmers_; }
    std::size_t  num_unitigs() const noexcept { return n_unitigs_; }
    std::size_t  super_len()  const noexcept { return super_len_; }
    std::size_t  bwt_len()    const noexcept { return bwt_len_; }

    /// Checks both orientations.
    bool find(KmerWord kmer) const noexcept {
        if (kmer_len_ == 0 || bwt_len_ == 0) return false;

        if (backward_search(kmer)) return true;

        const KmerWord rev_comp = reverse_complement(kmer, kmer_len_);

        return rev_comp != kmer && backward_search(rev_comp);
    }

    /// Searches find_batch() keeps in flight, each of which waits on a cache miss per step.
    static constexpr unsigned kFindBatch = 16;

    /// True when the rank blocks outgrow the last-level cache, which is when find_batch() beats find().
    /// Inside the cache a rank costs a few cycles, so the lane bookkeeping outweighs the overlap.
    /// With the same layout, hg38 chr1 (132 MB index, 32 MB L3) ran 1.1x faster batched, while a 12 MB index ran 1.4x slower.
    bool batch_pays() const noexcept {
        std::size_t rank_bytes = 0;

        for (const FmRankBits& rank_bits : ranks_) rank_bytes += rank_bits.nblocks * sizeof(FmRankBits::Block);

        return rank_bytes > last_level_cache_bytes();
    }

    /// Same answers as calling find() on each k-mer, written to @p hit[i] as 1 or 0.
    ///
    /// The searches run in lockstep, kFindBatch at a time, so rank misses from different k-mers overlap (Slotin, "Memory-level parallelism").
    /// Every step first prefetches the two blocks each active search will read, and only then advances them all.
    /// A finished search frees its lane for the next k-mer.
    /// A k-mer absent in its given orientation is searched again as its reverse complement, which keeps the answers identical to find().
    void find_batch(const KmerWord* kmers, std::size_t count, uint8_t* hit) const noexcept {
        if (kmer_len_ == 0 || bwt_len_ == 0) { std::fill(hit, hit + count, uint8_t(0)); return; }

        struct Lane { KmerWord pat; uint64_t sa_lo, sa_hi; std::size_t idx; unsigned step; bool is_rc; };
        Lane lanes[kFindBatch];
        unsigned active = 0;
        std::size_t next = 0;
        const unsigned kmer_len = kmer_len_;
        auto start = [&](Lane& lane, std::size_t kmer_idx, bool is_rc) {
            lane.idx = kmer_idx; lane.is_rc = is_rc; lane.pat = is_rc ? reverse_complement(kmers[kmer_idx], kmer_len) : kmers[kmer_idx];
            lane.sa_lo = 1; lane.sa_hi = bwt_len_; lane.step = 0;
        };

        while (active < kFindBatch && next < count) start(lanes[active++], next++, false);

        while (active) {
            for (unsigned j = 0; j < active; ++j) {
                const Lane& lane = lanes[j];
                const FmRankBits& rank_bits = ranks_[(lane.pat >> (kBitsPerBase * lane.step)) & kBaseMask];

                __builtin_prefetch(rank_bits.block_of(lane.sa_lo - 1));
                __builtin_prefetch(rank_bits.block_of(lane.sa_hi));
            }

            for (unsigned lane_idx = 0; lane_idx < active;) {
                Lane& lane = lanes[lane_idx];
                const uint8_t sym = kSymBase0 + static_cast<uint8_t>((lane.pat >> (kBitsPerBase * lane.step)) & kBaseMask);
                const uint64_t rank_lo = rank_sym(sym, lane.sa_lo - 1);
                const uint64_t rank_hi = rank_sym(sym, lane.sa_hi);

                // This is one backward-search step, mapping both ends of the interval through LF.
                lane.sa_lo = smaller_counts_[sym] + rank_lo + 1;
                lane.sa_hi = smaller_counts_[sym] + rank_hi;
                ++lane.step;
                const bool miss = lane.sa_lo > lane.sa_hi;

                if (!miss && lane.step < kmer_len) { ++lane_idx; continue; }

                // A palindromic k-mer is its own reverse complement, so it is not searched a second time.
                if (miss && !lane.is_rc && reverse_complement(kmers[lane.idx], kmer_len) != kmers[lane.idx]) {
                    start(lane, lane.idx, true);
                    ++lane_idx;
                    continue;
                }

                hit[lane.idx] = miss ? 0 : 1;

                // With no k-mers left, the last active lane moves into the finished one's slot, which keeps the active lanes packed.
                if (next < count) { start(lane, next++, false); ++lane_idx; }
                else lanes[lane_idx] = lanes[--active];
            }
        }
    }

    /// Searches one orientation only.
    bool backward_search(KmerWord kmer) const noexcept {
        const unsigned kmer_len = kmer_len_;
        uint64_t sa_lo = 1, sa_hi = bwt_len_;   // 1-indexed SA interval [lo, hi]

        // Packing puts the first base in the high bits and the last (3') base in the LSB.
        // Backward search therefore consumes bits from offset 0 upward, position K-1 first and position 0 last.
        for (unsigned i = 0; i < kmer_len; ++i) {
            const unsigned base = static_cast<unsigned>((kmer >> (kBitsPerBase * i)) & kBaseMask);
            const uint8_t sym = kSymBase0 + static_cast<uint8_t>(base);
            const uint64_t rank_lo = rank_sym(sym, sa_lo - 1);
            const uint64_t rank_hi = rank_sym(sym, sa_hi);

            sa_lo = smaller_counts_[sym] + rank_lo + 1;
            sa_hi = smaller_counts_[sym] + rank_hi;

            if (sa_lo > sa_hi) return false;
        }

        return sa_lo <= sa_hi;
    }

    std::size_t bytes() const noexcept {
        std::size_t total = kMagicBytes;

        total += kHeaderWords * kSerialWordBytes;
        total += kCEntries * kSerialWordBytes;

        for (int b = 0; b < kNumBases; ++b) {
            total += kRankHeaderWords * kSerialWordBytes;
            total += ranks_[b].nwords() * kSerialWordBytes;
        }

        total += kSerialWordBytes;                               // trailing checksum
        return total;
    }

    double bits_per_entry() const noexcept {
        return n_kmers_ > 0 ? (bytes() * double(kBitsPerByte)) / n_kmers_ : 0.0;
    }

    // File format, little-endian
    //
    //   [0..4)   kMagic "OKFM"
    //   [..)     u64 k, n_kmers, n_unitigs_, super_len_, bwt_len_
    //   [..)     kCEntries = 7 u64 values of smaller_counts_
    //   [..)     for A, C, G and T, u64 nbits, u64 nwords, nwords x u64 bits
    //   [..)     u64 checksum of every u64 above
    //
    // load() validates everything it can: the k range, that each rank vector holds exactly bwt_len_ bits, and that smaller_counts_ is monotone up to bwt_len_.
    // It also checks that each base's ones match its smaller_counts_ entry, that the checksum agrees and that the file is consumed exactly.
    // So a truncated, corrupted or bit-flipped index fails with an error instead of causing out-of-bounds reads or silently wrong answers.
    /// Streaming, non-cryptographic checksum over u64 words, whose constants belong to the file format.
    struct Sum {
        uint64_t state = kMagic;
        void add(uint64_t word) noexcept { state = mix64(state ^ word, kGoldenRatio64) + kChecksumAdd; }
    };
    /// Goes through atomic_write(), so a failed write (a full disk included) leaves no file behind.
    bool save(const char* path) const noexcept {
      try {
        return atomic_write(path, [&](std::ofstream& out) {
        Sum sum;
        // Every word written also feeds the running checksum stored at the end.
        auto w64 = [&](uint64_t word) { sum.add(word); out.write(reinterpret_cast<const char*>(&word), kSerialWordBytes); };

        out.write(reinterpret_cast<const char*>(&kMagic), kMagicBytes);
        w64(kmer_len_); w64(n_kmers_); w64(n_unitigs_); w64(super_len_); w64(bwt_len_);

        for (int i = 0; i < kCEntries; ++i) w64(smaller_counts_[i]);

        // Only the raw bits of the four rank arrays are written. load() rebuilds the per-block counts with build_rank().
        for (int b = 0; b < kNumBases; ++b) {
            const FmRankBits& rank_bits = ranks_[b];

            w64(rank_bits.bit_count);
            w64(rank_bits.nwords());

            for (uint64_t i = 0; i < rank_bits.nwords(); ++i) {
                const uint64_t word = rank_bits.word(i);

                sum.add(word);
                out.write(reinterpret_cast<const char*>(&word), kSerialWordBytes);
            }
        }

        const uint64_t checksum = sum.state;

        out.write(reinterpret_cast<const char*>(&checksum), kSerialWordBytes);
        });
      } catch (...) { return false; }
    }

    /// Returns nullptr for a missing, malformed or truncated file, or one that fails its checksum.
    /// The caller owns the result.
    static FmIndex* load(const char* path) {
        std::ifstream in(path, std::ios::binary);

        if (!in) return nullptr;

        // Needed below to insist on an exact file size.
        in.seekg(0, std::ios::end);
        const std::streamoff fsize = in.tellg();

        if (fsize < 0) return nullptr;

        in.seekg(0, std::ios::beg);
        Sum sum;
        auto r64 = [&]() -> uint64_t {
            uint64_t word = 0; in.read(reinterpret_cast<char*>(&word), kSerialWordBytes); sum.add(word); return word;
        };
        uint32_t magic;

        in.read(reinterpret_cast<char*>(&magic), kMagicBytes);

        if (!in || magic != kMagic) return nullptr;

        std::unique_ptr<FmIndex> idx(new FmIndex());
        const uint64_t kmer_len = r64();

        idx->n_kmers_  = static_cast<std::size_t>(r64());
        idx->n_unitigs_        = static_cast<std::size_t>(r64());
        idx->super_len_= static_cast<std::size_t>(r64());
        idx->bwt_len_        = static_cast<std::size_t>(r64());

        for (int i = 0; i < kCEntries; ++i) idx->smaller_counts_[i] = r64();

        if (!in || kmer_len < 1 || kmer_len > kMaxK) return nullptr;

        idx->kmer_len_ = static_cast<unsigned>(kmer_len);
        const uint64_t bwt_len = idx->bwt_len_;

        // The C array must start at 0, end at the BWT length and never decrease. Anything else is corrupt.
        if (idx->smaller_counts_[0] != 0 || idx->smaller_counts_[kAlphabetSize] != bwt_len) return nullptr;

        for (int i = 0; i < kAlphabetSize; ++i) if (idx->smaller_counts_[i] > idx->smaller_counts_[i + 1]) return nullptr;

        // Each rank is stored as u64 nbits, u64 nwords and then nwords u64 words.
        // Checking against the real remaining file size rejects truncated files before anything is allocated.
        // seekg/tellg reports that size accurately, whereas std::filebuf::in_avail only reflects the read-ahead buffer.
        // All four rank bitvectors must fit in what's left of the file, so a corrupt length can't request an outsized allocation.
        if (uint64_t(kNumBases) * words_for_bits(bwt_len) > static_cast<uint64_t>((fsize - in.tellg()) / kSerialWordBytes)) return nullptr;

        idx->alloc_ranks(bwt_len);

        for (int b = 0; b < kNumBases; ++b) {
            FmRankBits& rank_bits = idx->ranks_[b];
            const uint64_t stored_bits = r64();
            const uint64_t stored_words = r64();

            if (!in || stored_bits != bwt_len || stored_words != words_for_bits(bwt_len)) return nullptr;

            const std::streamoff remain = fsize - in.tellg();

            if (stored_words > static_cast<uint64_t>(remain / kSerialWordBytes)) return nullptr;

            // Flat words go straight into the blocks, kWordsPerBlock at a time.
            for (uint64_t word_start = 0; word_start < stored_words; word_start += FmRankBits::kWordsPerBlock) {
                const uint64_t cnt = std::min<uint64_t>(FmRankBits::kWordsPerBlock, stored_words - word_start);
                uint64_t* dst = rank_bits.blocks[word_start / FmRankBits::kWordsPerBlock].words;

                in.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(cnt * kSerialWordBytes));

                if (!in) return nullptr;

                for (uint64_t j = 0; j < cnt; ++j) sum.add(dst[j]);
            }

            if (bit_offset(bwt_len) && stored_words && (rank_bits.word(stored_words - 1) >> bit_offset(bwt_len))) return nullptr;   // stray bits past N

            rank_bits.build_rank();
            const uint64_t ones = rank_bits.rank(bwt_len);
            const uint8_t sym = static_cast<uint8_t>(kSymBase0 + b);

            if (ones != idx->smaller_counts_[sym + 1] - idx->smaller_counts_[sym]) return nullptr;
        }

        uint64_t checksum = 0;

        in.read(reinterpret_cast<char*>(&checksum), kSerialWordBytes);

        if (!in || checksum != sum.state) return nullptr;

        if (in.tellg() != fsize) return nullptr;                 // trailing garbage

        return idx.release();
    }

private:
    /// Expects text = u0 X u1 X ... u{U-1} X and appends the sentinel $ itself.
    void build(unsigned kmer_len, std::vector<uint8_t>&& text, std::size_t n_unitigs,
               std::size_t n_kmers, unsigned nthreads) {
        kmer_len_ = kmer_len;

        n_kmers_ = n_kmers;
        n_unitigs_ = n_unitigs;
        super_len_ = text.size() - n_unitigs;
        text.push_back(kSymSentinel);
        bwt_len_ = text.size();
        nthreads = resolve_threads(nthreads);

        if (nthreads == 0) nthreads = 1;

        // smaller_counts_[c] counts symbols strictly smaller than c, which works because the BWT is a permutation of the text.
        uint64_t counts[kCEntries] = {0, 0, 0, 0, 0, 0, 0};

        for (uint8_t sym : text) ++counts[sym];

        uint64_t front = 0;

        for (int c = 0; c < kCEntries; ++c) { smaller_counts_[c] = front; front += counts[c]; }

        alloc_ranks(bwt_len_);
        build_kbwt(text, nthreads);

        for (int b = 0; b < kNumBases; ++b) ranks_[b].build_rank();
    }

    /// One row of the k-order BWT.
    /// Its key packs the window's bases (cut at the first non-base) left-aligned in 2k bits with zero padding.
    /// aux packs window length << kAuxLenShift | position << kAuxPosShift | BWT symbol.
    struct Row { uint64_t key, aux; };
    static constexpr unsigned kAuxLenShift = 40, kAuxPosShift = 3;
    static constexpr uint64_t kAuxSymMask = 7u;

    /// Sorts text positions by key, length and position, emitting their BWT symbols straight into the rank bitvectors.
    /// Row 0 belongs to the sentinel, which is the smallest window.
    void build_kbwt(const std::vector<uint8_t>& text, unsigned num_threads) {
        const std::size_t text_len = text.size();
        const unsigned kmer_len = kmer_len_;
        const uint64_t mask = kmer_mask(kmer_len);
        constexpr unsigned kBucketBits = 16;
        const unsigned shift = kBitsPerBase * kmer_len > kBucketBits ? kBitsPerBase * kmer_len - kBucketBits : 0;
        constexpr std::size_t kNbuckets = std::size_t(1) << kBucketBits;
        auto is_base = [&](std::size_t pos) { return text[pos] >= kSymBase0; };
        // Visits the windows of a position range in order, calling f(pos, key, len).
        auto windows = [&](std::size_t range_lo, std::size_t range_hi, auto&& emit_window) {
            uint64_t key = 0;
            unsigned win_len = 0;
            bool fresh = true;

            for (std::size_t i = range_lo; i < range_hi; ++i) {
                if (!is_base(i)) { key = 0; win_len = 0; fresh = true; emit_window(i, key, win_len); continue; }

                if (fresh) {
                    key = 0; win_len = 0;

                    while (win_len < kmer_len && i + win_len < text_len && is_base(i + win_len)) {
                        key = (key << kBitsPerBase) | uint64_t(text[i + win_len] - kSymBase0);
                        ++win_len;
                    }

                    key <<= kBitsPerBase * (kmer_len - win_len);
                    fresh = false;
                } else if (win_len == kmer_len && is_base(i + kmer_len - 1)) {     // roll in one base
                    key = ((key << kBitsPerBase) | uint64_t(text[i + kmer_len - 1] - kSymBase0)) & mask;
                } else {                                       // window shrinks
                    key = (key << kBitsPerBase) & mask;
                    --win_len;
                }

                emit_window(i, key, win_len);
            }
        };
        // Per-thread bucket histograms over fixed position ranges.
        // The sentinel row (n - 1) is emitted first, outside the sort.
        const std::size_t body = text_len - 1;
        std::vector<std::size_t> cut(num_threads + 1);

        for (unsigned t = 0; t <= num_threads; ++t) cut[t] = body * t / num_threads;

        std::vector<std::vector<uint64_t>> hist(num_threads, std::vector<uint64_t>(kNbuckets, 0));

        {
            tpoptoa::run_on_threads(num_threads, [&](unsigned thread_idx) {
                auto& thread_hist = hist[thread_idx];

                windows(cut[thread_idx], cut[thread_idx + 1], [&](std::size_t, uint64_t key, unsigned) { ++thread_hist[key >> shift]; });
            });
        }
        std::vector<uint64_t> bsize(kNbuckets, 0);

        for (unsigned t = 0; t < num_threads; ++t)
            for (std::size_t b = 0; b < kNbuckets; ++b) bsize[b] += hist[t][b];

        uint64_t row = 0;
        auto emit = [&](uint8_t sym) {
            if (sym >= kSymBase0) ranks_[sym - kSymBase0].set(row);

            ++row;
        };

        // $ sorts first, so its row comes out first, carrying the symbol just before it.
        const std::size_t sentinel_pos = text_len - 1;
        emit(sentinel_pos >= 1 ? text[sentinel_pos - 1] : kSymSentinel);

        // Each pass covers a bucket range of about rows_per_pass rows, bounded by kMinPassRows and kMaxPasses.
        const std::size_t rows_per_pass = std::max<std::size_t>(kMinPassRows, body / kMaxPasses);
        std::vector<Row> rows;
        std::vector<std::size_t> off(kNbuckets + 1), toff(num_threads);

        for (std::size_t pass_lo = 0; pass_lo < kNbuckets;) {
            std::size_t pass_hi = pass_lo, cnt = 0;

            while (pass_hi < kNbuckets && (pass_hi == pass_lo || cnt + bsize[pass_hi] <= rows_per_pass)) cnt += bsize[pass_hi++];

            if (cnt == 0) { pass_lo = pass_hi; continue; }

            rows.resize(cnt);
            off[pass_lo] = 0;

            for (std::size_t b = pass_lo; b < pass_hi; ++b) off[b + 1] = off[b] + bsize[b];

            // During the scatter, thread t writes bucket b starting at off[b] plus the rows earlier threads put in b.
            {
                tpoptoa::run_on_threads(num_threads, [&](unsigned thread_idx) {
                    std::vector<std::size_t> pos(pass_hi - pass_lo);

                    for (std::size_t b = pass_lo; b < pass_hi; ++b) {
                        std::size_t slot = off[b];

                        for (unsigned u = 0; u < thread_idx; ++u) slot += hist[u][b];

                        pos[b - pass_lo] = slot;
                    }

                    windows(cut[thread_idx], cut[thread_idx + 1], [&](std::size_t text_pos, uint64_t key, unsigned win_len) {
                        const std::size_t bucket = key >> shift;

                        if (bucket < pass_lo || bucket >= pass_hi) return;

                        const uint8_t bwt_sym = text_pos == 0 ? kSymSentinel : text[text_pos - 1];

                        rows[pos[bucket - pass_lo]++] = Row{key, (uint64_t(win_len) << kAuxLenShift) | (uint64_t(text_pos) << kAuxPosShift) | bwt_sym};
                    });
                });
            }
            // Rows within a bucket share their top bits, so each bucket sorts independently.
            {
                constexpr std::size_t kMinSortedEntries = 2;   // a single entry is already in order
                std::atomic<std::size_t> next{pass_lo};

                tpoptoa::run_on_threads(num_threads, [&](unsigned) {
                    for (std::size_t bucket; (bucket = next.fetch_add(1)) < pass_hi;) {
                        if (bsize[bucket] < kMinSortedEntries) continue;

                        std::sort(rows.begin() + off[bucket], rows.begin() + off[bucket + 1],
                                  [](const Row& lhs, const Row& rhs) {
                                      return lhs.key != rhs.key ? lhs.key < rhs.key : lhs.aux < rhs.aux;
                                  });
                    }
                });
            }
            for (const Row& row_entry : rows) emit(static_cast<uint8_t>(row_entry.aux & kAuxSymMask));

            pass_lo = pass_hi;
        }
    }

    uint64_t rank_sym(uint8_t sym, uint64_t pos) const noexcept {
        return ranks_[sym - kSymBase0].rank(pos);
    }

    unsigned     kmer_len_ = 0;
    std::size_t  n_kmers_ = 0;
    std::size_t  n_unitigs_ = 0;
    std::size_t  super_len_ = 0;
    std::size_t  bwt_len_ = 0;
    uint64_t     smaller_counts_[kCEntries] = {};   // the FM-index C array
    FmRankBits   ranks_[kNumBases];
    HugeVec<FmRankBits::Block> rank_store_;   // one allocation backing all four ranks_

    void alloc_ranks(uint64_t nbits) {
        const uint64_t blocks_per_rank = FmRankBits::blocks_for(nbits);

        rank_store_.assign(kNumBases * blocks_per_rank, FmRankBits::Block{});

        for (int b = 0; b < kNumBases; ++b) ranks_[b].init(nbits, rank_store_.data() + b * blocks_per_rank);
    }
};

}
