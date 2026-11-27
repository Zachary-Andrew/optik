//! `optik dbruijn` builds a compacted de Bruijn graph from reads or genomes.
//!
//! It is a thin driver over the C++ graph builder, the `okdbg_*` functions of src/ffi_dbg.cpp.
//!   1. `okdbg_build_from_fasta` counts canonical k-mers in src/kmer_count.cpp, recording the observed (k+1)-mers around each occurrence as in/out edge masks, i.e. the edge-centric graph of Cuttlefish 2 over the abundance-filtered k-mers, into a node table or straight into partitioned compaction.
//!   2. The maximal unitigs are walked in linear time, each k-mer in exactly one unitig, sorted so the output is identical at any thread count (include/dbg_unitigs.hpp, include/dbg_part.hpp).
//!   3. `okdbg_write_gfa` writes one `S` line per unitig with `LN:i`, the summed k-mer abundance `KC:i` and the mean abundance `km:f`, and one `L` line per bidirected edge with a (k-1)M overlap, derived from the same edge masks the unitigs were compacted over.
//!      It can write unitig FASTA instead.
//!
//! `bench/validate_unitigs.py` validates the output, requiring a k-mer set equal to an oracle, no duplicated k-mers, maximality, exact L links, and unitig-partition equality with Cuttlefish (`--read`) at min_count 1.
//! For min_count > 1 optik filters k-mers (nodes) while Cuttlefish filters (k+1)-mers (edges), so the partitions may differ by design.
use std::time::Instant;

use crate::DbgIdx;

/// True if `path` names unitig FASTA output, ending in .fa, .fasta or .fna, rather than GFA.
fn fasta_out(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();

    [".fa", ".fasta", ".fna"]
        .iter()
        .any(|ext| lower.ends_with(ext))
}

/// Builds the compacted de Bruijn graph of `input` and writes it to `output`.
///
/// * `input` holds the FASTA/FASTQ paths, several joined by `'\n'` for a repeated `-i`.
/// * `kmer_len` is the k-mer length, 1..=31, since the (k+1)-mer edges must fit 64 bits.
/// * `min_count` is the lowest k-mer abundance kept.
/// * `output` receives GFA1, or unitig FASTA when the name ends in .fa/.fasta/.fna.
/// * `save` is the optional `-s` path of a kdbg archive of the graph core.
/// * `threads` bounds the threads for counting, walking and writing.
///
/// Prints a summary of unitigs, links and length statistics on stdout.
/// A graph where no k-mer reaches `min_count` is written empty with a warning.
///
/// # Errors
///
/// Returns an error for k outside 1..=31, an input holding no k-mer of length `kmer_len`, or any build or I/O failure.
pub(crate) fn run(
    input: &str,
    kmer_len: usize,
    min_count: u32,
    output: &str,
    save: Option<&str>,
    threads: usize,
) -> Result<(), String> {
    if kmer_len == 0 || kmer_len > crate::DBG_MAX_K {
        return Err(format!(
            "dbruijn requires 1 <= k <= 31 (k+1-mer edge set needs 2*(k+1) <= 64 bits), got k={kmer_len}"
        ));
    }

    eprintln!(
        "[dbruijn] counting k-mers from {} (k={kmer_len}, min_count={min_count})...",
        input.replace('\n', ", ")
    );
    let start = Instant::now();

    let mut graph = match DbgIdx::build_from_fasta(
        input,
        kmer_len,
        min_count as u64,
        threads,
        // Saving needs the node table, which rules out partitioned compaction.
        save.is_some(),
    ) {
        Ok(built) => built,
        // No k-mer reaches the cutoff, so the graph is empty, which is a valid result (header-only GFA) and not a failure.
        // An input without any k-mer at all is a user error, as for build-spss.
        Err(err) if err.contains("no k-mers extracted") => {
            // Several -i inputs arrive joined by '\n'.
            // stdin cannot be read twice, so a stdin input is taken to have had k-mers.
            let mut any = false;

            for input_path in input.split('\n') {
                if input_path == "-" {
                    any = true;
                    continue;
                }

                // A record holds a k-mer as soon as it has a run of k unambiguous bases.
                crate::iterate_fasta(input_path, |_, seq| {
                    let mut run = 0usize;

                    for &base in seq {
                        run = if matches!(
                            base.to_ascii_uppercase(),
                            b'A' | b'C' | b'G' | b'T' | b'U'
                        ) {
                            run + 1
                        } else {
                            0
                        };
                        if run >= kmer_len {
                            any = true;
                            break;
                        }
                    }
                })?;
            }

            if !any {
                return Err(format!(
                    "no k-mers of length {kmer_len} in {} (empty input, or every record shorter than k)",
                    input.replace('\n', ", ")
                ));
            }

            // The empty graph is an empty FASTA, or a GFA with only its header line.
            std::fs::write(
                output,
                if fasta_out(output) {
                    ""
                } else {
                    "H\tVN:Z:1.0\n"
                },
            )
            .map_err(|err| err.to_string())?;
            eprintln!(
                "[dbruijn] warning: no k-mer reaches min_count={min_count}; wrote an empty graph"
            );
            // Metric lines still print so scripts parsing them see zeros.
            println!("input_kmers\t0");
            println!("unitigs\t0");
            return Ok(());
        }
        Err(err) => return Err(err),
    };
    let kmer_count = graph.size();

    eprintln!("[dbruijn] {kmer_count} k-mers after filtering");

    if let Some(path) = save {
        graph.save(path)?;
        eprintln!("[dbruijn] graph core saved -> {path}");
    }

    let stats = graph.write_gfa(output, threads, fasta_out(output))?;

    eprintln!(
        "[dbruijn] {} unitigs, {} links -> {output} with {threads} thread(s) in {:.2}s",
        stats.n_unitigs,
        stats.n_links,
        start.elapsed().as_secs_f64()
    );
    if stats.n_dangling != 0 {
        // A graph edge whose target is not a segment end means the unitig set is not a maximal partition, which is reported rather than hidden.
        eprintln!(
            "[dbruijn] WARNING: {} edge(s) point into a segment interior (non-maximal unitigs)",
            stats.n_dangling
        );
    }

    // Deterministic report on stdout.
    // `total_length` counts every segment base (k-1 overlaps between adjacent unitigs are counted twice), so it is not a genome length estimate.
    // Tab-separated metrics on stdout, one per line, for the TUI and scripts.
    if stats.n_unitigs > 0 {
        println!("unitigs\t{}", stats.n_unitigs);
        println!("links\t{}", stats.n_links);
        println!("kmers\t{}", stats.n_kmers);
        println!("total_length\t{}", stats.total_bp);
        println!("max_unitig\t{}", stats.max_len);
        println!("min_unitig\t{}", stats.min_len);
        println!(
            "avg_unitig\t{:.1}",
            stats.total_bp as f64 / stats.n_unitigs as f64
        );
        println!("n50\t{}", stats.n50);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::fasta_out;

    #[test]
    fn fasta_extensions_pick_unitig_fasta_and_everything_else_gfa() {
        assert!(fasta_out("unitigs.fa"));
        assert!(fasta_out("UNITIGS.FASTA"));
        assert!(fasta_out("out.fna"));
        assert!(!fasta_out("graph.gfa"));
        assert!(!fasta_out("graph.fa.gfa"));
    }
}
