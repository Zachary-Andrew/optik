#pragma once

/*
 * kmer_count.hpp - canonical k-mer counting with inline edge recording.
 *
 * It is shared by every build that counts k-mers, namely dbruijn, build-fm and build-spss --compact.
 * The counter records, next to each k-mer's abundance, which of its 8 possible (k+1)-mer extensions were seen in the input.
 * The de Bruijn graph is then the edge-centric graph of Cuttlefish 2 (Khan et al. 2022), with edges kept only between solid k-mers.
 */

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "../include/phase_timer.hpp"

namespace tpoptoa {

class DbgGraph;
class DbgGraphEF;
class PartUnitigs;

namespace detail {

/// count_into_graph() - count the canonical k-mers of a FASTA/FASTQ file straight into @graph.
///
/// The table of @graph is sized here.
/// Each k-mer seen at least @min_count times, i.e. a solid k-mer, gets its abundance and observed edges.
/// The edge masks are then restricted to solid neighbours.
/// @nthreads == 0 means resolve_threads(0), half the logical CPUs.
///
/// With @compact given, a set whose table exceeds the memory limit (or any set when OPTIK_DBG_COMPACT=1) is built into *@compact instead (DbgGraphEF, dbg_graph_ef.hpp) and @graph stays empty.
/// With @part given, no graph is built and the unitigs come from partitioned compaction (dbg_part.hpp).
void count_into_graph(const char* fasta_path, std::size_t kk,
                      uint32_t min_count, unsigned nthreads, DbgGraph& graph,
                      std::unique_ptr<DbgGraphEF>* compact = nullptr,
                      std::unique_ptr<PartUnitigs>* part = nullptr);

/// count_kmers_with_edges() - count_into_graph() exported as sorted arrays, for validation and tests.
///
///   out_kmers[i] = (canonical k-mer, abundance >= min_count)
///   out_masks[i] = in_mask << 4 | out_mask, in the canonical k-mer's frame,
///                  keeping only edges whose (k+1)-mer was observed in either
///                  orientation and whose neighbour k-mer is itself solid
///
/// @nthreads == 0 means half the logical CPUs.
void count_kmers_with_edges(const char* fasta_path, std::size_t kk,
                            uint32_t min_count, unsigned nthreads,
                            std::vector<std::pair<uint64_t, uint32_t>>& out_kmers,
                            std::vector<uint8_t>& out_masks);

}
}
