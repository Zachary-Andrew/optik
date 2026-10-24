/*
 * ffi_dbg.cpp - C ABI of the de Bruijn graph builder (dbruijn), the okdbg_* functions.
 *
 * okdbg_build_from_fasta() counts k-mers into one of three graph forms (DbgGraph table, DbgGraphEF compact keys, or PartUnitigs from partitioned compaction).
 * okdbg_write_gfa() writes the maximal unitigs as GFA1 (S and L lines, GFA spec) or as unitig FASTA with BCALM-style headers.
 * The output is identical at every thread count and in every graph form.
 */
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

using tpoptoa::ffi::set_error;

namespace {

/// GfaOut - buffered GFA text writer.
/// Small pieces collect in a large buffer flushed with fwrite, and block() writes the text that Fmt chunks produce, in order, instead of one fprintf parse per line.
class GfaOut {
public:
    explicit GfaOut(std::FILE* f) : f_(f) { buf_.reserve(kFlushBytes + 4096); }
    ~GfaOut() { flush(); }
    void put(const char* s, std::size_t n) {
        buf_.append(s, n);
        if (buf_.size() >= kFlushBytes) flush();
    }
    void put(const char* s) { put(s, std::strlen(s)); }
    bool flush() {
        if (!buf_.empty() && std::fwrite(buf_.data(), 1, buf_.size(), f_) != buf_.size()) ok_ = false;
        buf_.clear();
        return ok_;
    }
    /// block() - write a pre-formatted block, keeping buffer order.
    void block(const std::string& t) {
        flush();
        if (!t.empty() && std::fwrite(t.data(), 1, t.size(), f_) != t.size()) ok_ = false;
    }
private:
    static constexpr std::size_t kFlushBytes = std::size_t(8) << 20;   // flush threshold, bytes
    std::FILE* f_;
    std::string buf_;
    bool ok_ = true;
};

/// Fmt - append-only text formatter, one per thread chunk, whose string GfaOut emits in order.
/// Numbers go through std::to_chars, which is exact, locale-free and identical to %.2f for doubles.
struct Fmt {
    std::string b;
    void put(char c) { b.push_back(c); }
    void put(const char* s, std::size_t n) { b.append(s, n); }
    void num(uint64_t v) {
        char t[24];
        const auto r = std::to_chars(t, t + sizeof(t), v);
        b.append(t, r.ptr);
    }
    void fixed2(double v) {
        char t[64];
        const auto r = std::to_chars(t, t + sizeof(t), v, std::chars_format::fixed, 2);
        b.append(t, r.ptr);
    }
};

}

extern "C" {

/// OkDbg - a built graph in one of three forms.
/// g is the open-addressing table, ef the compact graph the counter builds instead when the table would not fit the memory limit (dbg_graph_ef.hpp), and pu the unitigs of partitioned compaction (dbg_part.hpp).
/// g is always allocated, while ef and pu are set only when used.
struct OkDbg { tpoptoa::DbgGraph* g; std::unique_ptr<tpoptoa::DbgGraphEF> ef;
               std::unique_ptr<tpoptoa::PartUnitigs> pu; };
uint64_t okdbg_core_size(const OkDbg* h);

/// okdbg_build_from_fasta() - count the k-mers of @fasta_path (abundance at least @min_count) straight into a graph.
/// @keep_nodes keeps the node table (needed for the -s archive), which rules out partitioned compaction.
/// Returns nullptr with a message in @errbuf on failure.
OkDbg* okdbg_build_from_fasta(const char* fasta_path, uint64_t k,
                              uint64_t min_count, unsigned nthreads, int keep_nodes,
                              char* errbuf, uint64_t errbuf_len) {
    try {
        auto* h = new OkDbg();
        h->g = new tpoptoa::DbgGraph(static_cast<unsigned>(k));
        // Single-parse lock-free counting straight into the graph's table.
        // Edge masks come from the observed (K+1)-windows around each occurrence, restricted to solid neighbours, which gives the edge-centric graph over the min_count-solid k-mers.
        tpoptoa::detail::count_into_graph(fasta_path,
            static_cast<std::size_t>(k), static_cast<uint32_t>(min_count),
            nthreads, *h->g, &h->ef,
            // the node table is kept when its archive is asked for (-s)
            !keep_nodes && tpoptoa::part_mode_wanted(static_cast<unsigned>(k)) ? &h->pu : nullptr);
        if (okdbg_core_size(h) == 0) {
            delete h->g; delete h;
            set_error(errbuf, errbuf_len, "no k-mers extracted from input");
            return nullptr;
        }
        return h;
    } catch (const std::exception& e) {
        set_error(errbuf, errbuf_len, e.what());
        return nullptr;
    }
}

/// okdbg_save() - write the kdbg archive (DbgGraph::save()); -1 for a graph without a node table or on an I/O error.
int okdbg_save(const OkDbg* h, const char* path) {
    try {
        if (h->pu) return -1;                    // no node table (built without -s)
        return (h->ef ? h->ef->save(path) : h->g->save(path)) ? 0 : -1;
    } catch (...) {
        return -1;
    }
}

void okdbg_free(OkDbg* h) { delete h->g; delete h; }

uint64_t okdbg_core_size(const OkDbg* h)  {
    return h->pu ? h->pu->core_size() : h->ef ? h->ef->core_size() : h->g->core_size();
}

typedef void (*DbgUnitigCb)(void* cookie, const char* seq, uint32_t len,
                            uint64_t coverage);

/// Summary of a written GFA, which okdbg_write_gfa fills.
struct OkDbgGfaStats {
    uint64_t n_unitigs;
    uint64_t total_bp;      // sum of segment lengths, with overlaps counted twice
    uint64_t n_kmers;       // sum of len - k + 1 over segments
    uint64_t max_len, min_len, n50;
    uint64_t n_links;       // distinct L lines written
    uint64_t n_dangling;    // graph edges whose target is not a segment end, a validation
                            // signal, since nonzero means the unitig decomposition is not maximal or not a partition
};

}  // extern "C"

namespace {

/// SegmentLengths - segment lengths for the GFA summary, as a histogram where lengths below kHistLimit are counted, longer ones listed, so memory does not grow by one word per segment.
struct SegmentLengths {
    static constexpr std::size_t kHistLimit = std::size_t(1) << 16;
    std::vector<uint64_t> hist = std::vector<uint64_t>(kHistLimit, 0);
    std::vector<uint32_t> big;
    void add(uint32_t len) {
        if (len < hist.size()) ++hist[len];
        else big.push_back(len);
    }
};

/// fill_gfa_stats() - the summary of a written graph, i.e. the count, total, minimum, maximum and N50 of the segment lengths.
/// Sorts @lens.big.
/// Does nothing if @st is null.
void fill_gfa_stats(OkDbgGfaStats* st, SegmentLengths& lens, std::size_t n_seg, uint64_t nk,
                    uint64_t n_links, uint64_t dangling) {
    if (!st) return;
    std::vector<uint64_t>& len_hist = lens.hist;
    std::vector<uint32_t>& len_big = lens.big;
    std::memset(st, 0, sizeof(*st));
    st->n_unitigs = n_seg;
    st->n_kmers = nk;
    st->n_links = n_links;
    st->n_dangling = dangling;
    if (n_seg) {
        // Lengths in decreasing order, the long ones sorted and then the histogram
        std::sort(len_big.begin(), len_big.end(), std::greater<uint32_t>());
        uint64_t tot = 0;
        for (const uint32_t l : len_big) tot += l;
        for (std::size_t l = 0; l < len_hist.size(); ++l) tot += len_hist[l] * l;
        st->total_bp = tot;
        std::size_t lo = 0;
        while (lo < len_hist.size() && !len_hist[lo]) ++lo;
        std::size_t hi = len_hist.size();
        while (hi > 0 && !len_hist[hi - 1]) --hi;
        st->max_len = len_big.empty() ? hi - 1 : len_big.front();
        st->min_len = lo < len_hist.size() ? lo : len_big.back();
        uint64_t cum = 0;
        bool found = false;
        for (const uint32_t l : len_big) {
            cum += l;
            if (cum * 2 >= tot) { st->n50 = l; found = true; break; }
        }
        for (std::size_t l = hi; !found && l-- > lo;) {
            cum += len_hist[l] * l;
            if (len_hist[l] && cum * 2 >= tot) { st->n50 = l; found = true; }
        }
    }
}

/// A link (from, tgt), where from = segment << 1 | leaves reversed and tgt = segment << 1 | entered reversed.
using Link = std::pair<uint64_t, uint64_t>;

/// write_link_lines() - write the L lines, the sorted, duplicate-free union of the per-thread sorted link lists @part.
/// Returns the lines written.
///
/// The lists are cut into slices by source end.
/// A slice merges its part of all lists and is formatted by one thread, and the slices are written in order.
/// Repeats share their source, so they meet in one slice, and the text equals one global merge's.
uint64_t write_link_lines(GfaOut& out, const std::vector<std::vector<Link>>& part, std::size_t n_seg,
                          unsigned nthreads, unsigned K) {
    uint64_t n_links = 0;
    char ov[16];
    const auto ovl = std::to_chars(ov, ov + sizeof(ov), K - 1);
    const uint64_t span = uint64_t(n_seg) * 2;                  // from < 2n
    const std::size_t NS = std::max<std::size_t>(1, std::size_t(nthreads) * 16);
    const std::size_t WIN = std::size_t(nthreads) * 4;
    std::vector<Fmt> buf(WIN);
    const std::size_t P = part.size();
    for (std::size_t w0 = 0; w0 < NS; w0 += WIN) {
        const std::size_t w1 = std::min(NS, w0 + WIN);
        uint64_t wl = 0;
#ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1) reduction(+ : wl)
#endif
        for (int64_t sl = int64_t(w0); sl < int64_t(w1); ++sl) {
            Fmt& o = buf[std::size_t(sl) - w0];
            o.b.clear();
            const uint64_t flo = span * uint64_t(sl) / NS, fhi = span * uint64_t(sl + 1) / NS;
            std::vector<std::size_t> at(P), end(P);
            for (std::size_t t = 0; t < P; ++t) {
                at[t] = std::lower_bound(part[t].begin(), part[t].end(), Link{flo, 0}) - part[t].begin();
                end[t] = std::lower_bound(part[t].begin(), part[t].end(), Link{fhi, 0}) - part[t].begin();
            }
            bool have_prev = false;
            Link prev{0, 0};
            for (;;) {
                std::size_t best = P;
                for (std::size_t t = 0; t < P; ++t)
                    if (at[t] < end[t] && (best == P || part[t][at[t]] < part[best][at[best]]))
                        best = t;
                if (best == P) break;
                const Link l = part[best][at[best]++];
                if (have_prev && l == prev) continue;
                prev = l; have_prev = true;
                o.put("L\t", 2); o.num(l.first >> 1);
                o.put((l.first & 1) ? "\t-\t" : "\t+\t", 3);
                o.num(l.second >> 1);
                o.put((l.second & 1) ? "\t-\t" : "\t+\t", 3);
                o.put(ov, static_cast<std::size_t>(ovl.ptr - ov));
                o.put("M\n", 2);
                ++wl;
            }
        }
        n_links += wl;
        for (std::size_t c = w0; c < w1; ++c) out.block(buf[c - w0].b);
    }
    return n_links;
}

}

/// write_gfa_impl() - the GFA and FASTA writer over any graph form, described at okdbg_write_gfa().
template <typename G>
static int write_gfa_impl(G& g, const char* path, unsigned nthreads, OkDbgGfaStats* st, int fasta)
{
    try {
        const unsigned K = g.k();
        const uint64_t mask = K >= 32 ? ~uint64_t(0)
                                      : ((uint64_t(1) << (2u * K)) - 1);
        std::FILE* f = std::fopen(path, "w");
        if (!f) return -1;

        auto enc = [](uint8_t c) -> uint64_t {
            switch (c) { case 'A': return 0; case 'C': return 1;
                         case 'G': return 2; default: return 3; }
        };
        auto pack = [&](const uint8_t* s) {
            uint64_t v = 0;
            for (unsigned i = 0; i < K; ++i) v = (v << 2) | enc(s[i]);
            return v;
        };
        auto rc = [&](uint64_t v) { return tpoptoa::reverse_complement(v, K); };

        // Oriented end k-mers -> (segment id << 1 | entered_reverse)
        //   first k-mer of S   enters the segment in '+'
        //   rc(last k-mer of S) enters the segment in '-'
        tpoptoa::PhaseTimer pt;
        struct Ends { uint64_t first, last; };
        std::vector<Ends> ends;                 // GFA only, since the links match segment ends
        SegmentLengths seg_lengths;             // for the statistics
        std::size_t n_seg = 0;
        uint64_t nk = 0;
        GfaOut out(f);
        if (!fasta) out.put("H\tVN:Z:1.0\n");
        nthreads = tpoptoa::resolve_threads(nthreads);
        // S lines come from chunks of the sorted unitigs, formatted on all threads and written in order, so the file is the same at any thread count.
        // Working through a window of chunks at a time bounds the text held in memory.
        g.unitigs_sorted([&](std::size_t n, const auto& at) {
            n_seg = n;
            if (!fasta) ends.resize(n);
            constexpr std::size_t kSLineChunk = std::size_t(1) << 15;
            const std::size_t nch = (n + kSLineChunk - 1) / kSLineChunk;
            const std::size_t WIN = std::size_t(nthreads) * 4;
            std::vector<Fmt> buf(WIN);
            std::vector<std::vector<uint32_t>> clen(WIN);    // each chunk's segment lengths
            for (std::size_t w0 = 0; w0 < nch; w0 += WIN) {
                const std::size_t w1 = std::min(nch, w0 + WIN);
                uint64_t wk = 0;
#ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1) reduction(+ : wk)
#endif
                for (int64_t c = int64_t(w0); c < int64_t(w1); ++c) {
                    Fmt& o = buf[std::size_t(c) - w0];
                    o.b.clear();
                    std::vector<uint32_t>& cl = clen[std::size_t(c) - w0];
                    cl.clear();
                    const std::size_t lo = std::size_t(c) * kSLineChunk, hi = std::min(n, lo + kSLineChunk);
                    for (std::size_t id = lo; id < hi; ++id) {
                        const auto u = at(id);
                        const uint8_t* seq = u.seq;
                        const uint32_t len = u.len;
                        const uint64_t kmers = len >= K ? len - K + 1 : 0;
                        if (fasta) {
                            o.put('>'); o.num(id);
                            o.put(" LN:i:", 6); o.num(len);
                            o.put(" KC:i:", 6); o.num(u.cov);
                            o.put(" km:f:", 6); o.fixed2(kmers ? static_cast<double>(u.cov) / kmers : 0.0);
                            o.put('\n');
                            o.put(reinterpret_cast<const char*>(seq), len);
                            o.put('\n');
                        } else {
                            o.put("S\t", 2); o.num(id); o.put('\t');
                            o.put(reinterpret_cast<const char*>(seq), len);
                            o.put("\tLN:i:", 6); o.num(len);
                            o.put("\tKC:i:", 6); o.num(u.cov);
                            o.put("\tkm:f:", 6); o.fixed2(kmers ? static_cast<double>(u.cov) / kmers : 0.0);
                            o.put('\n');
                        }
                        if (!fasta) ends[id] = Ends{pack(seq), pack(seq + len - K)};
                        cl.push_back(len);
                        wk += kmers;
                    }
                }
                nk += wk;
                for (std::size_t c = w0; c < w1; ++c) {
                    out.block(buf[c - w0].b);
                    for (const uint32_t l : clen[c - w0]) seg_lengths.add(l);
                }
            }
        }, nthreads);
        pt.lap("gfa: unitigs + S lines");
        if (fasta) {
            const bool wrote = out.flush();
            const bool ok = std::fclose(f) == 0 && wrote;
            fill_gfa_stats(st, seg_lengths, n_seg, nk, 0, 0);
            return ok ? 0 : -1;
        }

        // Successor bases of an ORIENTED k-mer E per the canonical masks.
        auto succ = [&](uint64_t E) -> unsigned {
            const uint64_t r = rc(E);
            const uint64_t c = E < r ? E : r;
            unsigned m = 0;
            if (E == c) m |= g.out_mask(c);
            if (r == c) {                       // E is the reverse strand
                const unsigned in = g.in_mask(c);
                for (unsigned b = 0; b < 4; ++b)
                    if (in & (1u << (b ^ 3u))) m |= 1u << b;
            }
            return m;
        };
        auto canon = [&](uint64_t x) { const uint64_t r = rc(x); return x < r ? x : r; };
        nthreads = tpoptoa::resolve_threads(nthreads);
        const uint64_t n_ends = ends.size();
        const int64_t ne = static_cast<int64_t>(n_ends);

        // 1. Successor masks of every segment end, the only masks the links need.
        //    The low nibble is the exit of i forward through its last k-mer, and the high nibble the exit of i reversed through rc of its first.
        std::vector<uint8_t> exit_mask(n_ends);
        constexpr bool kIsPartitioned = std::is_same_v<G, tpoptoa::PartUnitigs>;
        if constexpr (kIsPartitioned) {
            // partitioned compaction kept every segment's exit masks
            for (int64_t i = 0; i < ne; ++i) exit_mask[i] = g.exit_mask(static_cast<std::size_t>(i));
        } else {
            // The loop steps one block of kBlk ends, in chunks of 256 blocks.
            // A chunk of 65536 iterations was 4M ends, so chr1's 2.4M ends ran on one thread and a human genome's 71M on 9 of 16.
            constexpr int64_t kBlk = 64;
#ifdef _OPENMP
            #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 256)
#endif
            for (int64_t b0 = 0; b0 < ne; b0 += kBlk) {
                const int64_t b1 = std::min(ne, b0 + kBlk);
                for (int64_t i = b0; i < b1; ++i) {
                    g.prefetch(canon(ends[i].last));
                    g.prefetch(canon(ends[i].first));
                }
                for (int64_t i = b0; i < b1; ++i)
                    exit_mask[i] = static_cast<uint8_t>(succ(ends[i].last) | (succ(rc(ends[i].first)) << 4));
            }
        }
        pt.lap("gfa: links: end masks");

        std::vector<std::vector<Link>> part(nthreads);
        std::vector<uint64_t> dang(nthreads, 0);
        if constexpr (kIsPartitioned) {
            // Partitioned unitigs have no node table, and both ends of every link are segment ends, so links come from matching the ends' exit windows.
            // End A's exit (k+1)-mer W_A and end B's W_B are the two traversals of one edge iff W_B = rc(W_A), with a palindromic W its own reverse (the self-mirrored link).
            // A -> flip(B) is then the link, kept from the end whose (from, tgt) is the canonical representative, as below.
            // Windows are grouped by canonical key in passes over ranges of a key hash, so memory stays bounded without an index of the end k-mers.
            // An exit window holds its canonical key and from << 1 | (window == key).
            // W_B = rc(W_A) iff the bits differ, or the key is its own reverse complement, in which case every exit of it pairs with every other.
            struct XRec { uint64_t key, fb; };
            auto canon1 = [&](uint64_t w) { const uint64_t r = tpoptoa::reverse_complement(w, K + 1); return w < r ? w : r; };
            constexpr std::size_t kExitBuckets = 4096;
            constexpr uint64_t kExitBucketSeed = 0x6A09E667F3BCC908ULL;
            // 64 MB per thread in memory, the rest in the spill directory
            constexpr uint64_t kExitMemoryPerThread = uint64_t(64) << 20;
            tpoptoa::spill::BucketSpill<XRec> xs(nthreads, kExitBuckets, kExitMemoryPerThread, tpoptoa::spill::spill_dir());
            {
                std::vector<std::exception_ptr> errs(nthreads);
                tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                    try {
                        for (uint64_t i = n_ends * t / nthreads; i < n_ends * (t + 1) / nthreads; ++i)
                            for (unsigned o = 0; o < 2; ++o) {
                                const uint64_t E = o == 0 ? ends[i].last : rc(ends[i].first);
                                const unsigned sm = (exit_mask[i] >> (4 * o)) & 0xFu;
                                for (unsigned b = 0; b < 4; ++b) {
                                    if (!(sm & (1u << b))) continue;
                                    const uint64_t w = (E << 2) | b;
                                    const uint64_t key = canon1(w);
                                    xs.add(t, tpoptoa::mix64(key, kExitBucketSeed) & (kExitBuckets - 1),
                                           XRec{key, (((i << 1) | o) << 1) | (w == key ? 1u : 0u)});
                                }
                            }
                    } catch (...) { errs[t] = std::current_exception(); }
                });
                for (auto& e : errs) if (e) std::rethrow_exception(e);
            }
            {
                std::atomic<std::size_t> nb{0};
                std::vector<std::exception_ptr> errs(nthreads);
                tpoptoa::run_on_threads(nthreads, [&](unsigned q) {
                    try {
                        std::vector<XRec> all;
                        for (std::size_t bk; (bk = nb.fetch_add(1)) < kExitBuckets;) {
                            xs.load(bk, all);
                            std::sort(all.begin(), all.end(), [](const XRec& x, const XRec& y) {
                                return x.key != y.key ? x.key < y.key : x.fb < y.fb; });
                            for (std::size_t i = 0; i < all.size();) {
                                std::size_t j = i + 1;
                                while (j < all.size() && all[j].key == all[i].key) ++j;
                                const bool pal = tpoptoa::reverse_complement(all[i].key, K + 1) == all[i].key;
                                for (std::size_t x = i; x < j; ++x) {
                                    bool any = false;
                                    for (std::size_t y = i; y < j; ++y) {
                                        if (!pal && (all[y].fb & 1) == (all[x].fb & 1)) continue;
                                        any = true;
                                        const uint64_t from = all[x].fb >> 1, tgt = (all[y].fb >> 1) ^ 1;
                                        if (std::make_pair(from, tgt) <= std::make_pair(tgt ^ 1, from ^ 1))
                                            part[q].push_back({from, tgt});
                                    }
                                    if (!any) ++dang[q];
                                }
                                i = j;
                            }
                        }
                    } catch (...) { errs[q] = std::current_exception(); }
                });
                for (auto& e : errs) if (e) std::rethrow_exception(e);
            }
            // each thread's list sorted (the L-line merge below expects it)
            {
                tpoptoa::run_on_threads(nthreads, [&](unsigned t) { std::sort(part[t].begin(), part[t].end()); });
            }
        } else {
            // 2. Where a link lands, i.e. the segment that an oriented k-mer X enters.
            //    The in-graph form serves fewer than 2^28 segments that make up at least 1/8 of the nodes.
            //    The node words, free once the walk is done, are cleared and each boundary node gets its 28-bit segment id and the flags is-first, first-is-forward, is-last and last-is-forward.
            //    A lookup is then one graph probe and no index is allocated, saving 2.7 GB on chr1 at k=15.
            //    Otherwise a separate open-addressing entry index is built.
            constexpr uint64_t kIdBits = 28;
            // Clearing the node words sweeps the whole table, which only pays when segments are a sizable share of the nodes.
            // On chr1, k=15 has 60M segments over 112M nodes, and k=31 has 2.4M over 204M, where the separate index takes 0.1 GB and builds faster.
            bool in_graph = !kIsPartitioned && n_ends < (uint64_t(1) << kIdBits) && n_ends * 8 >= g.core_size();
            // OPTIK_DBG_LINKS=graph|index forces a path, so tests can compare both.
            if (const char* e = std::getenv("OPTIK_DBG_LINKS")) {
                if (std::strcmp(e, "graph") == 0) in_graph = !kIsPartitioned && n_ends < (uint64_t(1) << kIdBits);
                else if (std::strcmp(e, "index") == 0) in_graph = false;
            }
            // flags above the segment id in a boundary node's word
            constexpr uint32_t kIsFirst = 1u << 28, kFirstForward = 1u << 29, kIsLast = 1u << 30, kLastForward = 1u << 31;
            struct ESlot { uint64_t key1, val; };
            std::vector<ESlot> etab;
            std::size_t ecap = 0;
            auto ehome = [&](uint64_t x) {
                return static_cast<std::size_t>(
                    (static_cast<__uint128_t>(tpoptoa::mix64(x, 0x9E3779B97F4A7C15ULL)) * ecap) >> 64);
            };
            if (in_graph) {
                g.aux_clear(nthreads);
                constexpr int64_t kBlk = 64;
    #ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 256)
    #endif
                for (int64_t b0 = 0; b0 < ne; b0 += kBlk) {
                    const int64_t b1 = std::min(ne, b0 + kBlk);
                    for (int64_t i = b0; i < b1; ++i) {
                        g.prefetch(canon(ends[i].first));
                        g.prefetch(canon(ends[i].last));
                    }
                    for (int64_t i = b0; i < b1; ++i) {
                        // One node per k-mer and one segment per node, so a node is written only by its own segment (both ends when the segment is a single k-mer).
                        const uint64_t F = ends[i].first, L = ends[i].last;
                        const uint64_t cf = canon(F), cl = canon(L);
                        uint32_t wf = uint32_t(i) | kIsFirst | (F == cf ? kFirstForward : 0u);
                        const uint32_t wl = uint32_t(i) | kIsLast | (L == cl ? kLastForward : 0u);
                        if (cf == cl) wf |= wl;
                        g.aux_set(g.rank(cf), wf);
                        if (cf != cl) g.aux_set(g.rank(cl), wl);
                    }
                }
            } else {
                ecap = n_ends * 2 + n_ends * 6 / 7 + 16;
                etab.assign(ecap, ESlot{0, 0});
                // In CAS-claimed slots, equal keys may sit in any order along a chain.
                // Lookups visit every match and each end's links are ordered below, so the output does not depend on it.
                auto eput = [&](uint64_t x, uint64_t v) {
                    for (std::size_t i = ehome(x);; i = (i + 1 == ecap) ? 0 : i + 1) {
                        uint64_t expect = 0;
                        if (__atomic_compare_exchange_n(&etab[i].key1, &expect, x + 1, false,
                                                        __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                            etab[i].val = v;
                            return;
                        }
                    }
                };
    #ifdef _OPENMP
                #pragma omp parallel for num_threads(nthreads) schedule(static, 65536)
    #endif
                for (int64_t i = 0; i < ne; ++i) {
                    eput(ends[i].first, uint64_t(i) << 1);
                    eput(rc(ends[i].last), (uint64_t(i) << 1) | 1);
                }
            }
            pt.lap(in_graph ? "gfa: links: entry lookup (graph ids)" : "gfa: links: entry lookup (index)");
            // tgt values, segment << 1 | entered reversed, of every segment entered through oriented k-mer X, at most kMaxEntries.
            // n = 0 marks a dangling edge.
            constexpr unsigned kMaxEntries = 8;
            auto entries = [&](uint64_t X, uint64_t* out, unsigned& n) -> void {
                n = 0;
                if (in_graph) {
                    const uint64_t c = canon(X);
                    const std::size_t xi = g.rank(c);
                    if (xi == G::kNpos) return;
                    const uint32_t w = g.aux(xi);
                    const uint64_t sid = w & ((uint32_t(1) << kIdBits) - 1);
                    // first k-mer == X means it is entered forward
                    if ((w & kIsFirst) && ((w & kFirstForward) ? X == c : X == rc(c))) out[n++] = sid << 1;
                    // rc(last k-mer) == X, i.e. last == rc(X), means it is entered reversed
                    if ((w & kIsLast) && ((w & kLastForward) ? X == rc(c) : X == c)) out[n++] = (sid << 1) | 1;
                    return;
                }
                for (std::size_t j = ehome(X); etab[j].key1; j = (j + 1 == ecap) ? 0 : j + 1)
                    if (etab[j].key1 == X + 1 && n < kMaxEntries) out[n++] = etab[j].val;
            };

            // 3. Each segment end's out-edges, split over threads.
            //    Every link is recorded from its source end, and ends are visited in increasing order, so each thread's list is sorted by construction once an end's own few links are ordered.
            {
                tpoptoa::run_on_threads(nthreads, [&](unsigned t) {
                    auto& links = part[t];
                    const uint64_t lo = n_ends * t / nthreads, hi = n_ends * (t + 1) / nthreads;
                    links.reserve((hi - lo) * 2);
                    constexpr uint64_t kBlk = 64;
                    for (uint64_t b0 = lo; b0 < hi; b0 += kBlk) {
                        const uint64_t b1 = std::min(hi, b0 + kBlk);
                        for (uint64_t i = b0; i < b1; ++i)          // prefetch every target
                            for (unsigned o = 0; o < 2; ++o) {
                                const uint64_t E = o == 0 ? ends[i].last : rc(ends[i].first);
                                const unsigned sm = (exit_mask[i] >> (4 * o)) & 0xFu;
                                for (unsigned b = 0; b < 4; ++b) {
                                    if (!(sm & (1u << b))) continue;
                                    const uint64_t X = ((E << 2) | b) & mask;
                                    if (in_graph) g.prefetch(canon(X));
                                    else __builtin_prefetch(&etab[ehome(X)]);
                                }
                            }
                        for (uint64_t i = b0; i < b1; ++i)
                            for (unsigned o = 0; o < 2; ++o) {
                                // exit k-mer of segment i traversed in orientation o
                                const uint64_t E = o == 0 ? ends[i].last : rc(ends[i].first);
                                const unsigned sm = (exit_mask[i] >> (4 * o)) & 0xFu;
                                const uint64_t from = (i << 1) | o;
                                const std::size_t g0 = links.size();
                                for (unsigned b = 0; b < 4; ++b) {
                                    if (!(sm & (1u << b))) continue;
                                    const uint64_t X = ((E << 2) | b) & mask;
                                    uint64_t tg[kMaxEntries]; unsigned nt;
                                    entries(X, tg, nt);
                                    if (!nt) { ++dang[t]; continue; }
                                    for (unsigned q = 0; q < nt; ++q) {
                                        // Canonical representative of the bidirected
                                        // edge, (from -> tgt) == (flip(tgt) -> flip(from)).
                                        // Each edge is found from both of its ends, and only the end whose (from, tgt) is the representative records it (a self-mirrored edge has one form and is found once).
                                        const uint64_t tgt = tg[q], rfrom = tgt ^ 1, rto = from ^ 1;
                                        if (std::make_pair(from, tgt) <= std::make_pair(rfrom, rto))
                                            links.push_back({from, tgt});
                                    }
                                }
                                for (std::size_t a2 = g0 + 1; a2 < links.size(); ++a2)
                                    for (std::size_t b2 = a2; b2 > g0 && links[b2] < links[b2 - 1]; --b2)
                                        std::swap(links[b2], links[b2 - 1]);
                            }
                    }
                    // sorted by construction in O(n) without a sort, and checked in one linear pass, with a full sort as the fallback
                    if (!std::is_sorted(links.begin(), links.end()))
                        std::sort(links.begin(), links.end());
                });
            }
            std::vector<ESlot>().swap(etab);
        }
        pt.lap("gfa: links: find + sort");
        // The entry index is not needed past this point.
        std::vector<uint8_t>().swap(exit_mask);
        std::vector<Ends>().swap(ends);
        uint64_t dangling = 0;
        for (uint64_t d : dang) dangling += d;
        const uint64_t n_links = write_link_lines(out, part, n_seg, nthreads, K);
        const bool wrote = out.flush();
        const bool ok = std::fclose(f) == 0 && wrote;
        pt.lap("gfa: links + L lines");

        fill_gfa_stats(st, seg_lengths, n_seg, nk, n_links, dangling);
        return ok ? 0 : -1;
    } catch (...) {
        return -1;
    }
}

extern "C" {

/// okdbg_write_gfa() - write the maximal unitigs of @h to @path as GFA1
///   S <id> <seq> LN:i KC:i km:f   KC is the summed k-mer abundance and km the
///                                 mean abundance per k-mer
///   L <a> <+/-> <b> <+/-> <k-1>M  one line per bidirected edge
/// Links come from the graph's own canonical in/out masks, so they are exactly the edge set the unitigs were compacted over.
/// Segments are written in the walker's deterministic order.
///
/// With @fasta != 0 it writes unitig FASTA instead, which SSHash, SBWT and FMSI pipelines index, with BCALM-style headers ">id LN:i:.. KC:i:.. km:f:.." and no links.
/// Fills @st if given.
/// Returns 0 on success, -1 on I/O or allocation failure.
int okdbg_write_gfa(OkDbg* h, const char* path, unsigned nthreads,
                    OkDbgGfaStats* st, int fasta) {
    return h->pu ? write_gfa_impl(*h->pu, path, nthreads, st, fasta)
         : h->ef ? write_gfa_impl(*h->ef, path, nthreads, st, fasta)
                 : write_gfa_impl(*h->g, path, nthreads, st, fasta);
}

}  // extern "C"
