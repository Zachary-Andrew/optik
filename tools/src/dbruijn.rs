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
    let l = path.to_ascii_lowercase();
    [".fa", ".fasta", ".fna"].iter().any(|e| l.ends_with(e))
}

/// Builds the compacted de Bruijn graph of `input` and writes it to `output`.
///
/// * `input` holds the FASTA/FASTQ paths, several joined by `'\n'` for a repeated `-i`.
/// * `k` is the k-mer length, 1..=31, since the (k+1)-mer edges must fit 64 bits.
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
/// Returns an error for k outside 1..=31, an input holding no k-mer of length `k`, or any build or I/O failure.
pub(crate) fn run(
    input: &str,
    k: usize,
    min_count: u32,
    output: &str,
    save: Option<&str>,
    threads: usize,
) -> Result<(), String> {
    if k == 0 || k >= 32 {
        return Err(format!(
            "dbruijn requires 1 <= k <= 31 (k+1-mer edge set needs 2*(k+1) <= 64 bits), got k={k}"
        ));
    }
    eprintln!(
        "[dbruijn] counting k-mers from {} (k={k}, min_count={min_count})...",
        input.replace('\n', ", ")
    );
    let t0 = Instant::now();

    let graph = match DbgIdx::build_from_fasta(input, k, min_count as u64, threads, save.is_some())
    {
        Ok(g) => g,
        // No k-mer reaches the cutoff, so the graph is empty, which is a valid result (header-only GFA) and not a failure.
        // An input without any k-mer at all is a user error, as for build-spss.
        Err(e) if e.contains("no k-mers extracted") => {
            // Several -i inputs arrive joined by '\n'.
            // stdin cannot be read twice, so a stdin input is taken to have had k-mers.
            let mut any = false;
            for f in input.split('\n') {
                if f == "-" {
                    any = true;
                    continue;
                }
                crate::iterate_fasta(f, |_, s| {
                    let mut run = 0usize;
                    for &b in s {
                        run = if matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U')
                        {
                            run + 1
                        } else {
                            0
                        };
                        if run >= k {
                            any = true;
                            break;
                        }
                    }
                })?;
            }
            if !any {
                return Err(format!(
                    "no k-mers of length {k} in {} (empty input, or every record shorter than k)",
                    input.replace('\n', ", ")
                ));
            }
            std::fs::write(
                output,
                if fasta_out(output) {
                    ""
                } else {
                    "H\tVN:Z:1.0\n"
                },
            )
            .map_err(|e| e.to_string())?;
            eprintln!(
                "[dbruijn] warning: no k-mer reaches min_count={min_count}; wrote an empty graph"
            );
            println!("input_kmers\t0");
            println!("unitigs\t0");
            return Ok(());
        }
        Err(e) => return Err(e),
    };
    let n = graph.size();
    eprintln!("[dbruijn] {n} k-mers after filtering");

    if let Some(path) = save {
        graph.save(path)?;
        eprintln!("[dbruijn] graph core saved -> {path}");
    }

    let st = graph.write_gfa(output, threads, fasta_out(output))?;
    eprintln!(
        "[dbruijn] {} unitigs, {} links -> {output} with {threads} thread(s) in {:.2}s",
        st.n_unitigs,
        st.n_links,
        t0.elapsed().as_secs_f64()
    );
    if st.n_dangling != 0 {
        // A graph edge whose target is not a segment end means the unitig set is not a maximal partition, which is reported rather than hidden.
        eprintln!(
            "[dbruijn] WARNING: {} edge(s) point into a segment interior (non-maximal unitigs)",
            st.n_dangling
        );
    }

    // Deterministic report on stdout.
    // `total_length` counts every segment base (k-1 overlaps between adjacent unitigs are counted twice), so it is not a genome length estimate.
    if st.n_unitigs > 0 {
        println!("unitigs\t{}", st.n_unitigs);
        println!("links\t{}", st.n_links);
        println!("kmers\t{}", st.n_kmers);
        println!("total_length\t{}", st.total_bp);
        println!("max_unitig\t{}", st.max_len);
        println!("min_unitig\t{}", st.min_len);
        println!(
            "avg_unitig\t{:.1}",
            st.total_bp as f64 / st.n_unitigs as f64
        );
        println!("n50\t{}", st.n50);
    }
    Ok(())
}
