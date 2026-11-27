/**
 * @file
 * @brief The `okdbg_*` C ABI behind dbruijn.
 *
 * okdbg_build_from_fasta() counts k-mers into one of three graph forms: the DbgGraph table, DbgGraphEF's compact keys, or PartUnitigs from partitioned compaction.
 * okdbg_write_gfa() then writes the maximal unitigs either as GFA1 (S and L lines) or as unitig FASTA with BCALM-style headers.
 * Output is identical at every thread count and across all three graph forms.
 */
#include "../include/bits.hpp"
#include "../include/thread_budget.hpp"
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "../include/dbg_graph.hpp"
#include "../include/dbg_graph_ef.hpp"
#include "../include/dbg_part.hpp"
#include "../include/spill.hpp"
#include "kmer_count.hpp"
#include "ffi_error.hpp"
#include "../include/kmer.hpp"

using tpoptoa::ffi::kStatusFailed;
using tpoptoa::ffi::kStatusOk;

constexpr std::size_t kGfaBufSlack = 4096;   // headroom GfaOut keeps past its flush threshold
/// Fmt formats a decimal u64, a fixed-point double and the "<k-1>M" overlap CIGAR in buffers of these sizes.
constexpr std::size_t kU64DigitsBuf = 24;
constexpr std::size_t kDoubleDigitsBuf = 64;
constexpr std::size_t kOverlapCigarBuf = 16;
constexpr int kKmDecimals = 2;               // decimals of km:f
/// Byte lengths of two fixed GFA pieces, an orientation with its tabs ("\t+\t") and an S-line tag prefix (" LN:i:").
constexpr std::size_t kOrientFieldLen = 3;
constexpr std::size_t kTagPrefixLen = 6;
/// Lines are formatted in nthreads * kSlicesPerThread slices, of which nthreads * kWindowPerThread are in flight.
constexpr std::size_t kSlicesPerThread = 16;
constexpr std::size_t kWindowPerThread = 4;
/// Every segment has two oriented ends, so link endpoints stay below 2n.
constexpr uint64_t kOrientations = 2;
/// OpenMP chunk sizes of the per-end loops and of the node-word sweep.
constexpr int kEndSchedChunk = 256;
constexpr int kNodeSweepChunk = 65536;
/// The in-graph entry index is faster once the ends are at least an eighth (kInGraphMinShare) of the nodes.
constexpr uint64_t kInGraphMinShare = 8;
/// The separate entry index gets kEntrySlotsPerEnd slots per end, plus kEntryExtraNum/kEntryExtraDen more and kEntrySlackSlots.
constexpr uint64_t kEntrySlotsPerEnd = 2;
constexpr uint64_t kEntryExtraNum = 6;
constexpr uint64_t kEntryExtraDen = 7;
constexpr uint64_t kEntrySlackSlots = 16;

using tpoptoa::ffi::set_error;

namespace {

/// Buffered GFA writer.
/// Small pieces collect in a large buffer that is flushed with fwrite.
/// block() writes whole chunks produced by Fmt, in order, which avoids one fprintf parse per line.
class GfaOut {
public:
    explicit GfaOut(std::FILE* file) : file_(file) { buf_.reserve(kFlushBytes + kGfaBufSlack); }

    ~GfaOut() { flush(); }
    void put(const char* text, std::size_t len) {
        buf_.append(text, len);

        if (buf_.size() >= kFlushBytes) flush();
    }
    void put(const char* text) { put(text, std::strlen(text)); }

    bool flush() {
        // GfaOut records the first write failure instead of throwing, and okdbg_write_gfa checks it once after the last block.
        if (!buf_.empty() && std::fwrite(buf_.data(), 1, buf_.size(), file_) != buf_.size()) write_ok_ = false;

        buf_.clear();
        return write_ok_;
    }
    void block(const std::string& text) {
        // Buffered small writes are flushed first so the large block lands after them in file order.
        flush();

        if (!text.empty() && std::fwrite(text.data(), 1, text.size(), file_) != text.size()) write_ok_ = false;
    }
private:
    static constexpr std::size_t kFlushBytes = std::size_t(8) << 20;
    std::FILE* file_;
    std::string buf_;
    bool write_ok_ = true;
};

/// Append-only formatter, one per thread chunk, whose text GfaOut later emits in order.
/// Numbers go through std::to_chars, which is exact, locale-free and matches %.2f for doubles.
struct Fmt {
    std::string text;
    void put(char ch) { text.push_back(ch); }
    void put(const char* chars, std::size_t len) { text.append(chars, len); }

    void num(uint64_t value) {
        char digits[kU64DigitsBuf];
        const auto result = std::to_chars(digits, digits + sizeof(digits), value);

        text.append(digits, result.ptr);
    }
    void fixed2(double value) {
        char digits[kDoubleDigitsBuf];
        const auto result = std::to_chars(digits, digits + sizeof(digits), value, std::chars_format::fixed, kKmDecimals);

        text.append(digits, result.ptr);
    }
};

}

extern "C" {

/// A built graph in one of three forms.
/// The table is always allocated.
/// The compact graph (dbg_graph_ef.hpp) replaces it when the table won't fit the memory limit, and partitioned unitigs (dbg_part.hpp) skip it altogether.
struct OkDbg { tpoptoa::DbgGraph* graph; std::unique_ptr<tpoptoa::DbgGraphEF> ef_graph;
               std::unique_ptr<tpoptoa::PartUnitigs> part_unitigs; };
uint64_t okdbg_core_size(const OkDbg* handle);

/// @p keep_nodes keeps the node table for the -s archive, which rules out partitioned compaction.
/// Failures return nullptr with a message in @p errbuf.
OkDbg* okdbg_build_from_fasta(const char* fasta_path, uint64_t kmer_len,
                              uint64_t min_count, unsigned nthreads, int keep_nodes,
                              char* errbuf, uint64_t errbuf_len) {
    try {
        auto* handle = new OkDbg();

        handle->graph = new tpoptoa::DbgGraph(static_cast<unsigned>(kmer_len));
        // Lock-free counting in a single parse, straight into the graph's table.
        // Edge masks come from the (K+1)-windows observed around each occurrence and are trimmed to solid neighbours.
        // The result is the edge-centric graph over the min_count-solid k-mers.
        tpoptoa::detail::count_into_graph(fasta_path,
            static_cast<std::size_t>(kmer_len), static_cast<uint32_t>(min_count),
            nthreads, *handle->graph, &handle->ef_graph,
            !keep_nodes && tpoptoa::part_mode_wanted(static_cast<unsigned>(kmer_len)) ? &handle->part_unitigs : nullptr);
        if (okdbg_core_size(handle) == 0) {
            delete handle->graph; delete handle;
            set_error(errbuf, errbuf_len, "no k-mers extracted from input");
            return nullptr;
        }

        return handle;
    } catch (const std::exception& err) {
        set_error(errbuf, errbuf_len, err.what());
        return nullptr;
    }
}

/// Fails for a graph without a node table, and on I/O errors.
int okdbg_save(const OkDbg* handle, const char* path) {
    try {
        if (handle->part_unitigs) return kStatusFailed;                    // built without -s, so there's no node table

        return (handle->ef_graph ? handle->ef_graph->save(path) : handle->graph->save(path)) ? kStatusOk : kStatusFailed;
    } catch (...) {
        return kStatusFailed;
    }
}

void okdbg_free(OkDbg* handle) { delete handle->graph; delete handle; }

uint64_t okdbg_core_size(const OkDbg* handle)  {
    return handle->part_unitigs ? handle->part_unitigs->core_size() : handle->ef_graph ? handle->ef_graph->core_size() : handle->graph->core_size();
}

typedef void (*DbgUnitigCb)(void* cookie, const char* seq, uint32_t len,
                            uint64_t coverage);

struct OkDbgGfaStats {
    uint64_t n_unitigs;
    uint64_t total_bp;      // sum of segment lengths, with overlaps counted twice
    uint64_t n_kmers;       // sum of len - k + 1 over segments
    uint64_t max_len, min_len, n50;
    uint64_t n_links;       // distinct L lines written
    uint64_t n_dangling;    // edges whose target isn't a segment end, where nonzero means the unitigs aren't a maximal partition
};

}  // extern "C"

namespace {

/// Segment lengths for the GFA summary.
/// Lengths below kHistLimit go into a histogram and only longer ones are listed, so memory doesn't grow by a word per segment.
struct SegmentLengths {
    static constexpr std::size_t kHistLimit = std::size_t(1) << 16;
    std::vector<uint64_t> hist = std::vector<uint64_t>(kHistLimit, 0);
    std::vector<uint32_t> big;

    void add(uint32_t len) {
        if (len < hist.size()) ++hist[len];
        else big.push_back(len);
    }
};

/// Fills count, total, min, max and N50 of the segment lengths, sorting @p lens.big on the way.
/// A null @p stats makes it a no-op.
void fill_gfa_stats(OkDbgGfaStats* stats, SegmentLengths& lens, std::size_t n_seg, uint64_t num_kmers,
                    uint64_t n_links, uint64_t dangling) {
    if (!stats) return;

    std::vector<uint64_t>& len_hist = lens.hist;
    std::vector<uint32_t>& len_big = lens.big;

    std::memset(stats, 0, sizeof(*stats));
    stats->n_unitigs = n_seg;
    stats->n_kmers = num_kmers;
    stats->n_links = n_links;
    stats->n_dangling = dangling;

    if (n_seg) {
        // Lengths are walked in decreasing order, so the sorted long ones come first and the histogram follows.
        std::sort(len_big.begin(), len_big.end(), std::greater<uint32_t>());
        uint64_t tot = 0;

        for (const uint32_t seg_len : len_big) tot += seg_len;

        for (std::size_t seg_len = 0; seg_len < len_hist.size(); ++seg_len) tot += len_hist[seg_len] * seg_len;

        stats->total_bp = tot;
        std::size_t hist_start = 0;

        while (hist_start < len_hist.size() && !len_hist[hist_start]) ++hist_start;

        std::size_t hist_end = len_hist.size();

        while (hist_end > 0 && !len_hist[hist_end - 1]) --hist_end;

        stats->max_len = len_big.empty() ? hist_end - 1 : len_big.front();
        stats->min_len = hist_start < len_hist.size() ? hist_start : len_big.back();
        uint64_t cum = 0;
        bool found = false;

        for (const uint32_t seg_len : len_big) {
            cum += seg_len;

            // N50 is reached once the lengths seen so far cover at least half the total.
            if (cum >= tot - cum) { stats->n50 = seg_len; found = true; break; }
        }

        for (std::size_t seg_len = hist_end; !found && seg_len-- > hist_start;) {
            cum += len_hist[seg_len] * seg_len;

            if (len_hist[seg_len] && cum >= tot - cum) { stats->n50 = seg_len; found = true; }
        }
    }
}

/// A link (from, tgt), where from = segment << 1 | leaves reversed and tgt = segment << 1 | entered reversed.
using Link = std::pair<uint64_t, uint64_t>;

/// Writes the sorted, duplicate-free union of the per-thread link lists and returns the line count.
///
/// The lists are cut into slices by source end, and each slice is merged and formatted by one thread before being written in order.
/// Duplicates share their source, so they always meet in the same slice, and the text equals what one global merge would produce.
uint64_t write_link_lines(GfaOut& out, const std::vector<std::vector<Link>>& part, std::size_t n_seg,
                          unsigned nthreads, unsigned kmer_len) {
    uint64_t n_links = 0;
    char overlap_cigar[kOverlapCigarBuf];
    const auto ovl = std::to_chars(overlap_cigar, overlap_cigar + sizeof(overlap_cigar), kmer_len - 1);
    const uint64_t span = uint64_t(n_seg) * kOrientations;
    const std::size_t num_slices = std::max<std::size_t>(1, std::size_t(nthreads) * kSlicesPerThread);
    const std::size_t window_slices = std::size_t(nthreads) * kWindowPerThread;
    std::vector<Fmt> buf(window_slices);
    const std::size_t num_parts = part.size();

    for (std::size_t win_start = 0; win_start < num_slices; win_start += window_slices) {
        const std::size_t win_end = std::min(num_slices, win_start + window_slices);
        uint64_t win_links = 0;

#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1) reduction(+ : win_links)
#endif
        for (int64_t sl = int64_t(win_start); sl < int64_t(win_end); ++sl) {
            Fmt& fmt = buf[std::size_t(sl) - win_start];

            fmt.text.clear();
            // A slice is a range of source ends, located in each sorted thread list by binary search.
            const uint64_t flo = span * uint64_t(sl) / num_slices, fhi = span * uint64_t(sl + 1) / num_slices;
            std::vector<std::size_t> cursor(num_parts), end(num_parts);

            for (std::size_t t = 0; t < num_parts; ++t) {
                cursor[t] = std::lower_bound(part[t].begin(), part[t].end(), Link{flo, 0}) - part[t].begin();
                end[t] = std::lower_bound(part[t].begin(), part[t].end(), Link{fhi, 0}) - part[t].begin();
            }

            bool have_prev = false;
            Link prev{0, 0};

            // A k-way merge over the thread lists yields the slice's links in order, and a repeat of the previous link is skipped.
            for (;;) {
                std::size_t best = num_parts;

                for (std::size_t t = 0; t < num_parts; ++t)
                    if (cursor[t] < end[t] && (best == num_parts || part[t][cursor[t]] < part[best][cursor[best]]))
                        best = t;

                if (best == num_parts) break;

                const Link link = part[best][cursor[best]++];

                if (have_prev && link == prev) continue;

                // An end id holds the segment in its high bits and the orientation in bit 0, and every link overlaps by k - 1 bases.
                prev = link; have_prev = true;
                fmt.put("L\t", sizeof("L\t") - 1); fmt.num(link.first >> 1);
                fmt.put((link.first & 1) ? "\t-\t" : "\t+\t", kOrientFieldLen);
                fmt.num(link.second >> 1);
                fmt.put((link.second & 1) ? "\t-\t" : "\t+\t", kOrientFieldLen);
                fmt.put(overlap_cigar, static_cast<std::size_t>(ovl.ptr - overlap_cigar));
                fmt.put("M\n", sizeof("M\n") - 1);
                ++win_links;
            }
        }

        n_links += win_links;

        for (std::size_t c = win_start; c < win_end; ++c) out.block(buf[c - win_start].text);
    }

    return n_links;
}

}

/// Shared by every graph form, with the output documented at okdbg_write_gfa().
template <typename G>
static int write_gfa_impl(G& graph, const char* path, unsigned nthreads, OkDbgGfaStats* stats, int fasta)
{
    try {
        const unsigned kmer_len = graph.kmer_len();
        const uint64_t mask = tpoptoa::kmer_mask(kmer_len);
        std::FILE* file = std::fopen(path, "w");

        if (!file) return kStatusFailed;

        auto enc = [](uint8_t ch) -> uint64_t { return tpoptoa::base_encode(static_cast<char>(ch)); };   // unitigs are pure ACGT
        auto pack = [&](const uint8_t* seq) {
            uint64_t packed = 0;

            for (unsigned i = 0; i < kmer_len; ++i) packed = (packed << tpoptoa::kBitsPerBase) | enc(seq[i]);

            return packed;
        };
        auto rev_comp = [&](uint64_t kmer) { return tpoptoa::reverse_complement(kmer, kmer_len); };

        // Maps oriented end k-mers to (segment id << 1 | entered_reverse):
        //   the first k-mer of S enters the segment in '+'
        //   rc(last k-mer of S) enters it in '-'
        tpoptoa::PhaseTimer timer;

        struct Ends { uint64_t first, last; };
        std::vector<Ends> ends;                 // GFA only, because links are matched on segment ends
        SegmentLengths seg_lengths;
        std::size_t n_seg = 0;
        uint64_t num_kmers = 0;
        GfaOut out(file);

        if (!fasta) out.put("H\tVN:Z:1.0\n");

        nthreads = tpoptoa::resolve_threads(nthreads);
        // S lines come from chunks of the sorted unitigs, formatted on all threads and written in order, so the file never depends on the thread count.
        // Working through a window of chunks at a time also bounds the text held in memory.
        graph.unitigs_sorted([&](std::size_t count, const auto& unitig_at) {
            n_seg = count;

            if (!fasta) ends.resize(count);

            constexpr std::size_t kSLineChunk = std::size_t(1) << 15;
            const std::size_t nch = (count + kSLineChunk - 1) / kSLineChunk;
            const std::size_t window_chunks = std::size_t(nthreads) * kWindowPerThread;
            std::vector<Fmt> buf(window_chunks);
            std::vector<std::vector<uint32_t>> clen(window_chunks);    // each chunk's segment lengths

            // Chunks are formatted in parallel a window at a time, then written in id order, so the output is the same at every thread count while memory stays bounded by the window.
            for (std::size_t win_start = 0; win_start < nch; win_start += window_chunks) {
                const std::size_t win_end = std::min(nch, win_start + window_chunks);
                uint64_t win_kmers = 0;

#ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1) reduction(+ : win_kmers)
#endif
                for (int64_t c = int64_t(win_start); c < int64_t(win_end); ++c) {
                    Fmt& fmt = buf[std::size_t(c) - win_start];

                    fmt.text.clear();
                    std::vector<uint32_t>& chunk_lens = clen[std::size_t(c) - win_start];

                    chunk_lens.clear();
                    const std::size_t chunk_lo = std::size_t(c) * kSLineChunk, chunk_hi = std::min(count, chunk_lo + kSLineChunk);

                    for (std::size_t id = chunk_lo; id < chunk_hi; ++id) {
                        const auto unitig = unitig_at(id);
                        const uint8_t* seq = unitig.seq;
                        const uint32_t len = unitig.len;
                        const uint64_t kmers = len >= kmer_len ? len - kmer_len + 1 : 0;

                        // SSHash, SBWT and FMSI expect BCALM-style tags in the FASTA header. GFA puts them after the sequence.
                        if (fasta) {
                            fmt.put('>'); fmt.num(id);
                            fmt.put(" LN:i:", kTagPrefixLen); fmt.num(len);
                            fmt.put(" KC:i:", kTagPrefixLen); fmt.num(unitig.cov);
                            fmt.put(" km:f:", kTagPrefixLen); fmt.fixed2(kmers ? static_cast<double>(unitig.cov) / kmers : 0.0);
                            fmt.put('\n');
                            fmt.put(reinterpret_cast<const char*>(seq), len);
                            fmt.put('\n');
                        } else {
                            fmt.put("S\t", sizeof("S\t") - 1); fmt.num(id); fmt.put('\t');
                            fmt.put(reinterpret_cast<const char*>(seq), len);
                            fmt.put("\tLN:i:", kTagPrefixLen); fmt.num(len);
                            fmt.put("\tKC:i:", kTagPrefixLen); fmt.num(unitig.cov);
                            fmt.put("\tkm:f:", kTagPrefixLen); fmt.fixed2(kmers ? static_cast<double>(unitig.cov) / kmers : 0.0);
                            fmt.put('\n');
                        }

                        // GFA output also needs each segment's first and last k-mer, because the links are built from them afterwards.
                        if (!fasta) ends[id] = Ends{pack(seq), pack(seq + len - kmer_len)};

                        chunk_lens.push_back(len);
                        win_kmers += kmers;
                    }
                }

                num_kmers += win_kmers;

                for (std::size_t c = win_start; c < win_end; ++c) {
                    out.block(buf[c - win_start].text);

                    for (const uint32_t seg_len : clen[c - win_start]) seg_lengths.add(seg_len);
                }
            }
        }, nthreads);
        timer.lap("gfa: unitigs + S lines");

        if (fasta) {
            const bool wrote = out.flush();
            const bool closed = std::fclose(file) == 0 && wrote;

            fill_gfa_stats(stats, seg_lengths, n_seg, num_kmers, 0, 0);
            return closed ? kStatusOk : kStatusFailed;
        }

        // Successor bases of an oriented k-mer E, read from the canonical masks.
        auto succ = [&](uint64_t oriented) -> unsigned {
            const uint64_t rev = rev_comp(oriented);
            const uint64_t canon_kmer = oriented < rev ? oriented : rev;
            unsigned succ_mask = 0;

            if (oriented == canon_kmer) succ_mask |= graph.out_mask(canon_kmer);

            if (rev == canon_kmer) {                       // E is on the reverse strand
                const unsigned in_edges = graph.in_mask(canon_kmer);

                for (unsigned b = 0; b < tpoptoa::kNumBases; ++b)
                    if (in_edges & (1u << (b ^ tpoptoa::kComplementMask))) succ_mask |= 1u << b;
            }

            return succ_mask;
        };
        auto canon = [&](uint64_t kmer) { const uint64_t rev = rev_comp(kmer); return kmer < rev ? kmer : rev; };

        nthreads = tpoptoa::resolve_threads(nthreads);
        const uint64_t n_ends = ends.size();
        const int64_t num_ends_signed = static_cast<int64_t>(n_ends);

        // 1. Successor masks of every segment end, which are the only masks the links need.
        //    The low nibble is i's exit going forward through its last k-mer, and the high nibble its exit going reversed through rc of its first.
        std::vector<uint8_t> exit_mask(n_ends);
        constexpr bool kIsPartitioned = std::is_same_v<G, tpoptoa::PartUnitigs>;

        if constexpr (kIsPartitioned) {
            // Partitioned compaction already kept every segment's exit masks.
            for (int64_t i = 0; i < num_ends_signed; ++i) exit_mask[i] = graph.exit_mask(static_cast<std::size_t>(i));
        } else {
            // The loop steps one block of kBlk ends at a time, in chunks of kEndSchedChunk blocks.
            // A coarser chunk (65536 iterations, 4M ends) left chr1's 2.4M ends on one thread and a human genome's 71M on 9 of 16.
            constexpr int64_t kBlk = 64;

#ifdef _OPENMP
            #pragma omp parallel for num_threads(nthreads) schedule(dynamic, kEndSchedChunk)
#endif
            for (int64_t block_start = 0; block_start < num_ends_signed; block_start += kBlk) {
                const int64_t block_end = std::min(num_ends_signed, block_start + kBlk);

                for (int64_t i = block_start; i < block_end; ++i) {
                    graph.prefetch(canon(ends[i].last));
                    graph.prefetch(canon(ends[i].first));
                }

                for (int64_t i = block_start; i < block_end; ++i)
                    exit_mask[i] = static_cast<uint8_t>(succ(ends[i].last) | (succ(rev_comp(ends[i].first)) << tpoptoa::kOutEdgeBits));
            }
        }

        timer.lap("gfa: links: end masks");

        std::vector<std::vector<Link>> part(nthreads);
        std::vector<uint64_t> dang(nthreads, 0);

        if constexpr (kIsPartitioned) {
            // Partitioned unitigs have no node table, but both ends of every link are segment ends.
            // So links come from matching the ends' exit windows instead.
            // End A's exit (k+1)-mer W_A and end B's W_B traverse the same edge iff W_B = rc(W_A), and a palindromic W is its own reverse (a self-mirrored link).
            // The link is then A -> flip(B), recorded from whichever end holds the canonical representative, as below.
            // Windows are grouped by canonical key in passes over ranges of a key hash, which bounds memory without indexing the end k-mers.
            // Each exit window stores its canonical key plus from << 1 | (window == key).
            // W_B = rc(W_A) iff those low bits differ, or the key is its own reverse complement, in which case every exit pairs with every other.
            struct XRec { uint64_t key, from_bits; };
            auto canon1 = [&](uint64_t window) { const uint64_t rev = tpoptoa::reverse_complement(window, kmer_len + 1); return window < rev ? window : rev; };
            constexpr std::size_t kExitBuckets = 4096;
            constexpr uint64_t kExitBucketSeed = 0x6A09E667F3BCC908ULL;
            // Each thread keeps 64 MB in memory and spills the rest.
            constexpr uint64_t kExitMemoryPerThread = uint64_t(64) << 20;
            tpoptoa::spill::BucketSpill<XRec> exit_spill(nthreads, kExitBuckets, kExitMemoryPerThread, tpoptoa::spill::spill_dir());

            {
                std::vector<std::exception_ptr> errs(nthreads);

                tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                    try {
                        // Each thread takes a contiguous slice of the ends, and every successor of each oriented end becomes one exit record.
                        for (uint64_t i = n_ends * thread_idx / nthreads; i < n_ends * (thread_idx + 1) / nthreads; ++i)
                            for (unsigned o = 0; o < kOrientations; ++o) {
                                const uint64_t oriented = o == 0 ? ends[i].last : rev_comp(ends[i].first);
                                const unsigned succ_mask = (exit_mask[i] >> (tpoptoa::kOutEdgeBits * o)) & tpoptoa::kOutEdgeMask;

                                for (unsigned b = 0; b < tpoptoa::kNumBases; ++b) {
                                    if (!(succ_mask & (1u << b))) continue;

                                    const uint64_t window = (oriented << tpoptoa::kBitsPerBase) | b;
                                    const uint64_t key = canon1(window);

                                    exit_spill.add(thread_idx, tpoptoa::mix64(key, kExitBucketSeed) & (kExitBuckets - 1),
                                           // The record packs end i with its orientation and a bit for an already canonical window.
                                           XRec{key, (((i << 1) | o) << 1) | (window == key ? 1u : 0u)});
                                }
                            }
                    } catch (...) { errs[thread_idx] = std::current_exception(); }
                });
                for (auto& err : errs) if (err) std::rethrow_exception(err);
            }
            {
                std::atomic<std::size_t> next_bucket{0};
                std::vector<std::exception_ptr> errs(nthreads);

                tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                    try {
                        std::vector<XRec> all;

                        // Threads then claim whole buckets, and inside a bucket the records that share a key are the candidate pairs.
                        for (std::size_t bucket; (bucket = next_bucket.fetch_add(1)) < kExitBuckets;) {
                            exit_spill.load(bucket, all);
                            std::sort(all.begin(), all.end(), [](const XRec& lhs, const XRec& rhs) {
                                return lhs.key != rhs.key ? lhs.key < rhs.key : lhs.from_bits < rhs.from_bits; });
                            for (std::size_t run_start = 0; run_start < all.size();) {
                                std::size_t run_end = run_start + 1;

                                while (run_end < all.size() && all[run_end].key == all[run_start].key) ++run_end;

                                const bool pal = tpoptoa::reverse_complement(all[run_start].key, kmer_len + 1) == all[run_start].key;

                                for (std::size_t x = run_start; x < run_end; ++x) {
                                    bool any = false;

                                    for (std::size_t y = run_start; y < run_end; ++y) {
                                        if (!pal && (all[y].from_bits & 1) == (all[x].from_bits & 1)) continue;

                                        any = true;
                                        const uint64_t from = all[x].from_bits >> 1, tgt = (all[y].from_bits >> 1) ^ 1;

                                        // A link and its reverse describe the same edge, so only the smaller of the two is kept.
                                        if (std::make_pair(from, tgt) <= std::make_pair(tgt ^ 1, from ^ 1))
                                            part[thread_idx].push_back({from, tgt});
                                    }

                                    // An exit that pairs with no end leaves the graph, which the statistics count as a dangling edge.
                                    if (!any) ++dang[thread_idx];
                                }

                                run_start = run_end;
                            }
                        }
                    } catch (...) { errs[thread_idx] = std::current_exception(); }
                });
                for (auto& err : errs) if (err) std::rethrow_exception(err);
            }
            // The L-line merge below expects each thread's list to be sorted.
            {
                tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) { std::sort(part[thread_idx].begin(), part[thread_idx].end()); });
            }
        } else {
            // 2. Where a link lands, meaning the segment an oriented k-mer X enters.
            //    The in-graph form handles fewer than 2^kIdBits segments that make up at least an eighth of the nodes.
            //    Once the walk is done the node words are free, so each boundary node gets its segment id plus four flags there.
            //    The flags say is-first, first-is-forward, is-last and last-is-forward.
            //    A lookup is then one graph probe with no extra index, which saved 2.7 GB on chr1 at k=15.
            //    Otherwise a separate open-addressing entry index is built.
            constexpr uint64_t kIdBits = 28;
            // Clearing the node words sweeps the whole table, which only pays when segments are a sizable share of the nodes.
            // On chr1, k=15 has 60M segments over 112M nodes, while k=31 has 2.4M over 204M, where the separate index needs only 0.1 GB and builds faster.
            bool in_graph = !kIsPartitioned && n_ends < (uint64_t(1) << kIdBits) && n_ends * kInGraphMinShare >= graph.core_size();

            // OPTIK_DBG_LINKS=graph|index forces a path, so tests can compare both.
            if (const char* env = std::getenv("OPTIK_DBG_LINKS")) {
                if (std::strcmp(env, "graph") == 0) in_graph = !kIsPartitioned && n_ends < (uint64_t(1) << kIdBits);
                else if (std::strcmp(env, "index") == 0) in_graph = false;
            }

            // Flags stored above the segment id in a boundary node's word.
            constexpr uint32_t kIsFirst = 1u << kIdBits, kFirstForward = kIsFirst << 1, kIsLast = kIsFirst << 2, kLastForward = kIsFirst << 3;

            struct ESlot { uint64_t key1, val; };
            std::vector<ESlot> etab;
            std::size_t ecap = 0;
            auto ehome = [&](uint64_t kmer) {
                return static_cast<std::size_t>(
                    (static_cast<__uint128_t>(tpoptoa::mix64(kmer, tpoptoa::kGoldenRatio64)) * ecap) >> tpoptoa::kWordBits);
            };

            // When the ends are at least an eighth of the nodes, their entries are stored in the node words, and a separate hash table serves otherwise.
            if (in_graph) {
                graph.aux_clear(nthreads);
                constexpr int64_t kBlk = 64;

    #ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(dynamic, kEndSchedChunk)
    #endif
                for (int64_t block_start = 0; block_start < num_ends_signed; block_start += kBlk) {
                    const int64_t block_end = std::min(num_ends_signed, block_start + kBlk);

                    // Prefetching the whole block before any write overlaps its misses.
                    for (int64_t i = block_start; i < block_end; ++i) {
                        graph.prefetch(canon(ends[i].first));
                        graph.prefetch(canon(ends[i].last));
                    }

                    for (int64_t i = block_start; i < block_end; ++i) {
                        // Each node belongs to exactly one segment, so only that segment writes it (from both ends when it's a single k-mer).
                        const uint64_t first_kmer = ends[i].first, last_kmer = ends[i].last;
                        const uint64_t canon_first = canon(first_kmer), canon_last = canon(last_kmer);
                        uint32_t first_word = uint32_t(i) | kIsFirst | (first_kmer == canon_first ? kFirstForward : 0u);
                        const uint32_t last_word = uint32_t(i) | kIsLast | (last_kmer == canon_last ? kLastForward : 0u);

                        if (canon_first == canon_last) first_word |= last_word;

                        graph.aux_set(graph.rank(canon_first), first_word);

                        if (canon_first != canon_last) graph.aux_set(graph.rank(canon_last), last_word);
                    }
                }
            } else {
                ecap = n_ends * kEntrySlotsPerEnd + n_ends * kEntryExtraNum / kEntryExtraDen + kEntrySlackSlots;
                etab.assign(ecap, ESlot{0, 0});
                // CAS-claimed slots may hold equal keys in any order along a chain.
                // Lookups visit every match and each end's links get ordered below, so the output never depends on that order.
                auto eput = [&](uint64_t key, uint64_t value) {
                    for (std::size_t slot = ehome(key);; slot = (slot + 1 == ecap) ? 0 : slot + 1) {
                        uint64_t expect = 0;

                        if (__atomic_compare_exchange_n(&etab[slot].key1, &expect, key + 1, false,
                                                        __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                            etab[slot].val = value;
                            return;
                        }
                    }
                };
    #ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(static, kNodeSweepChunk)
    #endif
                for (int64_t i = 0; i < num_ends_signed; ++i) {
                    eput(ends[i].first, uint64_t(i) << 1);
                    eput(rev_comp(ends[i].last), (uint64_t(i) << 1) | 1);
                }
            }

            timer.lap(in_graph ? "gfa: links: entry lookup (graph ids)" : "gfa: links: entry lookup (index)");
            // Collects the tgt values (segment << 1 | entered reversed) of every segment entered through X, up to kMaxEntries.
            // An empty result marks a dangling edge.
            constexpr unsigned kMaxEntries = 8;
            auto entries = [&](uint64_t succ, uint64_t* out, unsigned& count) -> void {
                count = 0;

                if (in_graph) {
                    const uint64_t canon_kmer = canon(succ);
                    const std::size_t node = graph.rank(canon_kmer);

                    if (node == G::kNpos) return;

                    const uint32_t word = graph.aux(node);
                    const uint64_t sid = word & ((uint32_t(1) << kIdBits) - 1);

                    // A first k-mer equal to X enters the segment forward.
                    if ((word & kIsFirst) && ((word & kFirstForward) ? succ == canon_kmer : succ == rev_comp(canon_kmer))) out[count++] = sid << 1;

                    // rc(last k-mer) == X, in other words last == rc(X), enters it reversed.
                    if ((word & kIsLast) && ((word & kLastForward) ? succ == rev_comp(canon_kmer) : succ == canon_kmer)) out[count++] = (sid << 1) | 1;

                    return;
                }

                for (std::size_t slot = ehome(succ); etab[slot].key1; slot = (slot + 1 == ecap) ? 0 : slot + 1)
                    if (etab[slot].key1 == succ + 1 && count < kMaxEntries) out[count++] = etab[slot].val;
            };

            // 3. Each segment end's out-edges, split over threads.
            //    Every link is recorded from its source end, and ends are visited in increasing order.
            //    So once an end's own few links are ordered, each thread's list is sorted by construction.
            {
                tpoptoa::run_on_threads(nthreads, [&](unsigned thread_idx) {
                    auto& links = part[thread_idx];
                    const uint64_t range_lo = n_ends * thread_idx / nthreads, range_hi = n_ends * (thread_idx + 1) / nthreads;

                    links.reserve((range_hi - range_lo) * kOrientations);
                    constexpr uint64_t kBlk = 64;

                    for (uint64_t block_start = range_lo; block_start < range_hi; block_start += kBlk) {
                        const uint64_t block_end = std::min(range_hi, block_start + kBlk);

                        for (uint64_t i = block_start; i < block_end; ++i)          // prefetch all targets first
                            for (unsigned o = 0; o < kOrientations; ++o) {
                                const uint64_t oriented = o == 0 ? ends[i].last : rev_comp(ends[i].first);
                                const unsigned succ_mask = (exit_mask[i] >> (tpoptoa::kOutEdgeBits * o)) & tpoptoa::kOutEdgeMask;

                                for (unsigned b = 0; b < tpoptoa::kNumBases; ++b) {
                                    if (!(succ_mask & (1u << b))) continue;

                                    const uint64_t succ = ((oriented << tpoptoa::kBitsPerBase) | b) & mask;

                                    if (in_graph) graph.prefetch(canon(succ));
                                    else __builtin_prefetch(&etab[ehome(succ)]);
                                }
                            }

                        for (uint64_t i = block_start; i < block_end; ++i)
                            for (unsigned o = 0; o < kOrientations; ++o) {
                                const uint64_t oriented = o == 0 ? ends[i].last : rev_comp(ends[i].first);
                                const unsigned succ_mask = (exit_mask[i] >> (tpoptoa::kOutEdgeBits * o)) & tpoptoa::kOutEdgeMask;
                                const uint64_t from = (i << 1) | o;
                                const std::size_t sorted_from = links.size();

                                for (unsigned b = 0; b < tpoptoa::kNumBases; ++b) {
                                    if (!(succ_mask & (1u << b))) continue;

                                    const uint64_t succ = ((oriented << tpoptoa::kBitsPerBase) | b) & mask;
                                    uint64_t targets[kMaxEntries]; unsigned num_targets;

                                    entries(succ, targets, num_targets);

                                    // This successor is not a segment end in the graph. The link would dangle, and n_dangling counts it instead.
                                    if (!num_targets) { ++dang[thread_idx]; continue; }

                                    for (unsigned q = 0; q < num_targets; ++q) {
                                        // (from -> tgt) and (flip(tgt) -> flip(from)) are the same bidirected edge.
                                        // Each edge is found from both ends, so only the end holding the smaller form records it, and a self-mirrored edge has just one form anyway.
                                        const uint64_t tgt = targets[q], rfrom = tgt ^ 1, rto = from ^ 1;

                                        if (std::make_pair(from, tgt) <= std::make_pair(rfrom, rto))
                                            links.push_back({from, tgt});
                                    }
                                }

                                // Each oriented end adds only a few links, so an insertion sort keeps the list ordered as it grows.
                                for (std::size_t a2 = sorted_from + 1; a2 < links.size(); ++a2)
                                    for (std::size_t b2 = a2; b2 > sorted_from && links[b2] < links[b2 - 1]; --b2)
                                        std::swap(links[b2], links[b2 - 1]);
                            }
                    }

                    // Construction already sorted the list in O(n), and one linear pass confirms it, with a full sort as a fallback.
                    if (!std::is_sorted(links.begin(), links.end()))
                        std::sort(links.begin(), links.end());
                });
            }
            std::vector<ESlot>().swap(etab);
        }

        timer.lap("gfa: links: find + sort");
        std::vector<uint8_t>().swap(exit_mask);
        std::vector<Ends>().swap(ends);
        uint64_t dangling = 0;

        for (uint64_t count : dang) dangling += count;

        const uint64_t n_links = write_link_lines(out, part, n_seg, nthreads, kmer_len);
        const bool wrote = out.flush();
        const bool closed = std::fclose(file) == 0 && wrote;

        timer.lap("gfa: links + L lines");

        fill_gfa_stats(stats, seg_lengths, n_seg, num_kmers, n_links, dangling);
        return closed ? kStatusOk : kStatusFailed;
    } catch (...) {
        return kStatusFailed;
    }
}

extern "C" {

/// Writes the maximal unitigs as GFA1, an S line per unitig and an L line per bidirected edge.
///
///     S <id> <seq> LN:i KC:i km:f
///     L <a> <+/-> <b> <+/-> <k-1>M
///
/// KC sums the abundances of the unitig's k-mers, while km is their mean.
/// Links come from the graph's own canonical masks, so they are exactly the edges the unitigs were compacted over.
/// Segments appear in the walker's deterministic order.
///
/// With @p fasta != 0 it writes unitig FASTA instead, the input SSHash, SBWT and FMSI pipelines expect.
/// Those records have BCALM-style headers (">id LN:i:.. KC:i:.. km:f:..") and no links.
/// I/O and allocation failures return kStatusFailed.
int okdbg_write_gfa(OkDbg* handle, const char* path, unsigned nthreads,
                    OkDbgGfaStats* stats, int fasta) {
    return handle->part_unitigs ? write_gfa_impl(*handle->part_unitigs, path, nthreads, stats, fasta)
         : handle->ef_graph ? write_gfa_impl(*handle->ef_graph, path, nthreads, stats, fasta)
                 : write_gfa_impl(*handle->graph, path, nthreads, stats, fasta);
}

}  // extern "C"
