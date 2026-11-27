#pragma once

/**
 * @file
 * @brief Canonical k-mer counting with inline edge recording.
 *
 * Every build that counts k-mers uses this: dbruijn, build-fm and build-spss --compact.
 * Alongside each abundance, the counter records which of the 8 possible (k+1)-mer extensions appeared.
 * The graph therefore comes out edge-centric, as in Cuttlefish 2 (Khan et al. 2022), with edges only between solid k-mers.
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

/// Counts the canonical k-mers of a FASTA/FASTQ file straight into @p graph, sizing its table here.
/// K-mers seen at least @p min_count times keep their abundance and their edges to other solid k-mers.
///
/// When @p compact is given and the table won't fit the memory limit, the graph goes into *@p compact instead (dbg_graph_ef.hpp).
/// OPTIK_DBG_COMPACT=1 forces that path, and @p graph then stays empty.
/// With @p part, no graph is built at all, because the unitigs come from partitioned compaction (dbg_part.hpp).
void count_into_graph(const char* fasta_path, std::size_t kmer_len,
                      uint32_t min_count, unsigned nthreads, DbgGraph& graph,
                      std::unique_ptr<DbgGraphEF>* compact = nullptr,
                      std::unique_ptr<PartUnitigs>* part = nullptr);

/// Does the same counting as count_into_graph() but returns sorted arrays, for validation and tests.
///   out_kmers[i] = (canonical k-mer, abundance)
///   out_masks[i] = in_mask << 4 | out_mask in the canonical frame, solid neighbours only
void count_kmers_with_edges(const char* fasta_path, std::size_t kmer_len,
                            uint32_t min_count, unsigned nthreads,
                            std::vector<std::pair<uint64_t, uint32_t>>& out_kmers,
                            std::vector<uint8_t>& out_masks);

}
}
