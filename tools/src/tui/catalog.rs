//! Describes every tool the TUI can launch, grouped by category, along with each parameter its form offers.

use super::*;

/// The groups the tool-selection screen shows.
#[derive(Clone, Copy, PartialEq)]
pub(super) enum ToolCategory {
    Alignment,
    Assembly,
    Index,
}

impl ToolCategory {
    pub(super) fn name(self) -> &'static str {
        match self {
            ToolCategory::Alignment => "Alignment",
            ToolCategory::Assembly => "Assembly",
            ToolCategory::Index => "Indexing",
        }
    }
    pub(super) fn color(self) -> Color {
        match self {
            ToolCategory::Alignment => Color::Cyan,
            ToolCategory::Assembly => Color::Green,
            ToolCategory::Index => Color::Magenta,
        }
    }
}

/// One tool of the selection screen and the parameters its form asks for.
#[derive(Clone, Copy)]
pub(super) struct ToolInfo {
    pub(super) name: &'static str,
    /// The optik subcommand the tool runs.
    pub(super) cmd: &'static str,
    pub(super) icon: &'static str,
    pub(super) desc: &'static str,
    pub(super) category: ToolCategory,
    pub(super) params: &'static [ParamDef],
}

/// One form field, with its label, the command-line flag it fills, its kind, a default value and help text.
/// An empty value is not passed to the command, and a flag ending in '=' is joined to its value as one argument (--cs=long).
#[derive(Clone, Copy)]
pub(super) struct ParamDef {
    pub(super) name: &'static str,
    pub(super) flag: &'static str,
    pub(super) kind: ParamKind,
    pub(super) default: &'static str,
    pub(super) desc: &'static str,
    /// The command fails without this field, so the form marks it with `*` and refuses to run until it is filled.
    pub(super) required: bool,
}

/// What a field holds, i.e. a path to read through the file browser, a path to write, a number, free text, or an on/off switch that passes its flag alone.
#[derive(Clone, Copy, PartialEq)]
pub(super) enum ParamKind {
    FileRead,
    FileWrite,
    Int,
    Text,
    Flag,
}

/// The value a Flag field holds when it is on.
/// Off is the empty string.
pub(super) const FLAG_ON: &str = "on";

pub(super) const THREADS: ParamDef = ParamDef {
    name: "Threads",
    flag: "-t",
    kind: ParamKind::Int,
    default: "",
    desc: "Worker threads; empty = half of the logical CPUs",
    required: false,
};

/// The k-mer counting limits shared by dbruijn, build-fm and build-spss --compact.
pub(super) const MAX_MEM: ParamDef = ParamDef {
    name: "Memory limit",
    flag: "--max-mem",
    kind: ParamKind::Text,
    default: "",
    desc: "Counting memory, e.g. 8G or 500M; past it buffers spill to disk; empty = half the available memory",
    required: false,
};

pub(super) const TMP_DIR: ParamDef = ParamDef {
    name: "Spill directory",
    flag: "--tmp-dir",
    kind: ParamKind::Text,
    default: "",
    desc: "Directory for spilled buffers; empty = the output file's directory",
    required: false,
};

pub(super) const GRAPH_MODE: ParamDef = ParamDef {
    name: "Graph mode",
    flag: "--graph",
    kind: ParamKind::Text,
    default: "",
    desc: "auto, part or table; identical output, part uses about 1/4 of the memory; empty = auto",
    required: false,
};

/// Every tool the TUI offers, in display order and grouped by category.
pub(super) const TOOLS: &[ToolInfo] = &[
    ToolInfo {
        name: "Map",
        cmd: "map",
        icon: " M ",
        desc: "Read mapping against reference genome using minimizer seeds; for pre-extracted k-mers use Query tool instead",
        category: ToolCategory::Alignment,
        params: &[
            ParamDef {
                name: "Reference FASTA",
                flag: "-r",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Reference genome in FASTA format",
                required: true,
            },
            ParamDef {
                name: "Reads FASTA/FASTQ",
                flag: "-q",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Sequencing reads in FASTA or FASTQ format (R1 of a pair)",
                required: true,
            },
            ParamDef {
                name: "Mate reads (R2)",
                flag: "-q",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Optional: R2 of paired-end reads (FR library, same order as R1)",
                required: false,
            },
            ParamDef {
                name: "k-mer length",
                flag: "-k",
                kind: ParamKind::Int,
                default: "",
                desc: "Seed length (1-32); empty = chosen from the reads, 21 or 15 for noisy long reads",
                required: false,
            },
            ParamDef {
                name: "Window size",
                flag: "-w",
                kind: ParamKind::Int,
                default: "10",
                desc: "Sliding window for minimizer selection, smaller = more seeds (1-256, 0 = auto-tune from the read length); reads shorter than k+w-1 bases yield NO minimizers (all unmapped)",
                required: false,
            },
            ParamDef {
                name: "Max edit distance",
                flag: "-e",
                kind: ParamKind::Int,
                default: "auto",
                desc: "Short-read edit budget; auto = 15% of the read length, at least 5",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "Read group",
                flag: "-R",
                kind: ParamKind::Text,
                default: "",
                desc: "@RG line such as @RG\\tID:s1\\tSM:sample1; adds RG:Z to every record",
                required: false,
            },
            ParamDef {
                name: "cs tag",
                flag: "--cs=",
                kind: ParamKind::Text,
                default: "",
                desc: "short or long to add minimap2's cs tag; empty = NM and MD only",
                required: false,
            },
            ParamDef {
                name: "Output SAM",
                flag: "-o",
                kind: ParamKind::FileWrite,
                default: "",
                desc: "Write SAM output to file (stdout if omitted)",
                required: false,
            },
        ],
    },
    ToolInfo {
        name: "Dbruijn",
        cmd: "dbruijn",
        icon: " D ",
        desc: "De Bruijn graph assembly: k-mer counting, filtering, unitig output",
        category: ToolCategory::Assembly,
        params: &[
            ParamDef {
                name: "Input FASTA/FASTQ",
                flag: "-i",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Input sequencing reads in FASTA or FASTQ format",
                required: true,
            },
            ParamDef {
                name: "k-mer length",
                flag: "-k",
                kind: ParamKind::Int,
                default: "31",
                desc: "K-mer size for De Bruijn graph, larger = more specific (1-31)",
                required: false,
            },
            ParamDef {
                name: "Min k-mer count",
                flag: "-c",
                kind: ParamKind::Int,
                default: "2",
                desc: "Minimum k-mer abundance to include, higher = less noise; use 1 for a genome",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "Output GFA",
                flag: "-o",
                kind: ParamKind::FileWrite,
                default: "graph.gfa",
                desc: "GFA graph, or unitig FASTA if the name ends in .fa/.fasta/.fna",
                required: false,
            },
            ParamDef {
                name: "Graph archive",
                flag: "-s",
                kind: ParamKind::FileWrite,
                default: "",
                desc: "Optional kdbg archive of the graph core (not with graph mode part)",
                required: false,
            },
            MAX_MEM,
            TMP_DIR,
            GRAPH_MODE,
        ],
    },
    ToolInfo {
        name: "Build SPSS",
        cmd: "build-spss",
        icon: " B ",
        desc: "Build and save an exact SPSS k-mer membership index",
        category: ToolCategory::Index,
        params: &[
            ParamDef {
                name: "Input FASTA/FASTQ",
                flag: "-i",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Sequences to index (unitigs give the smallest index; any input is exact)",
                required: true,
            },
            ParamDef {
                name: "k-mer length",
                flag: "-k",
                kind: ParamKind::Int,
                default: "31",
                desc: "K-mer length for indexing (1-32)",
                required: false,
            },
            ParamDef {
                name: "Minimizer length",
                flag: "-m",
                kind: ParamKind::Int,
                default: "0",
                desc: "0 = auto (estimated from the input; exact at any m); must be < k",
                required: false,
            },
            ParamDef {
                name: "Compact",
                flag: "--compact",
                kind: ParamKind::Flag,
                default: "",
                desc: "Count k-mers and index their maximal unitigs instead of the input sequences",
                required: false,
            },
            ParamDef {
                name: "Min k-mer count",
                flag: "-c",
                kind: ParamKind::Int,
                default: "",
                desc: "With Compact only: minimum abundance; empty = 1, every k-mer",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "Output index",
                flag: "-o",
                kind: ParamKind::FileWrite,
                default: "index.spss",
                desc: "Path to write the serialized index",
                required: true,
            },
            MAX_MEM,
            TMP_DIR,
            GRAPH_MODE,
        ],
    },
    ToolInfo {
        name: "Build FM",
        cmd: "build-fm",
        icon: " F ",
        desc: "Build and save an FM k-mer membership index (k-order BWT over the unitigs; exact)",
        category: ToolCategory::Index,
        params: &[
            ParamDef {
                name: "Input FASTA/FASTQ",
                flag: "-i",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Sequences whose k-mers are indexed",
                required: true,
            },
            ParamDef {
                name: "k-mer length",
                flag: "-k",
                kind: ParamKind::Int,
                default: "31",
                desc: "K-mer length for indexing (1-32)",
                required: false,
            },
            ParamDef {
                name: "Min k-mer count",
                flag: "-c",
                kind: ParamKind::Int,
                default: "",
                desc: "Minimum abundance; empty = 1, every k-mer (exact membership)",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "Output index",
                flag: "-o",
                kind: ParamKind::FileWrite,
                default: "index.fm",
                desc: "Path to write the serialized index",
                required: true,
            },
            MAX_MEM,
            TMP_DIR,
            GRAPH_MODE,
        ],
    },
    ToolInfo {
        name: "Query",
        cmd: "query",
        icon: " Q ",
        desc: "Query a saved SPSS or FM index with the k-mers of a FASTA/FASTQ file (exact membership counts)",
        category: ToolCategory::Index,
        params: &[
            ParamDef {
                name: "Index file",
                flag: "-i",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Saved index from Build SPSS or Build FM",
                required: true,
            },
            ParamDef {
                name: "Query FASTA/FASTQ",
                flag: "-q",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Sequences (k-mers or reads) whose k-mers are looked up in the index",
                required: true,
            },
            ParamDef {
                name: "k-mer length",
                flag: "-k",
                kind: ParamKind::Int,
                default: "",
                desc: "Empty = the index's own k (recorded in the file); a different k is an error",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "k-mer ids output",
                flag: "--ids",
                kind: ParamKind::FileWrite,
                default: "",
                desc: "SPSS index only: write each record's k-mer ids, -1 for absent",
                required: false,
            },
            ParamDef {
                name: "Skip checksum",
                flag: "--no-verify",
                kind: ParamKind::Flag,
                default: "",
                desc: "SPSS index only: skip the whole-file checksum at load",
                required: false,
            },
        ],
    },
    ToolInfo {
        name: "Update",
        cmd: "spss-update",
        icon: " U ",
        desc: "Add or remove the k-mers of a file in an SPSS index and write the compacted index",
        category: ToolCategory::Index,
        params: &[
            ParamDef {
                name: "Index file",
                flag: "-i",
                kind: ParamKind::FileRead,
                default: "",
                desc: "SPSS index to edit",
                required: true,
            },
            ParamDef {
                name: "Add FASTA/FASTQ",
                flag: "--add",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Insert every k-mer of this file",
                required: false,
            },
            ParamDef {
                name: "Remove FASTA/FASTQ",
                flag: "--remove",
                kind: ParamKind::FileRead,
                default: "",
                desc: "Erase every k-mer of this file, after the additions",
                required: false,
            },
            THREADS,
            ParamDef {
                name: "Output index",
                flag: "-o",
                kind: ParamKind::FileWrite,
                default: "",
                desc: "Path of the compacted index, which may be the edited index itself",
                required: true,
            },
        ],
    },
];
