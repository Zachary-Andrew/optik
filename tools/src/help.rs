//! Produces the command-line help.
//! Each command lists its flags and examples, and a mistyped command gets a did-you-mean suggestion.

/// The value flags and the switch flags subcommand `sub` accepts, or None for an unknown subcommand.
/// The TUI's forms are checked against this table.
pub(crate) fn subcommand_flags(
    sub: &str,
) -> Option<(&'static [&'static str], &'static [&'static str])> {
    let flags: (&'static [&'static str], &'static [&'static str]) = match sub {
        "map" => (
            &["-r", "-q", "-k", "-w", "-e", "-o", "-t", "-R"],
            &["--cs", "--cs=short", "--cs=long"],
        ),
        "dbruijn" => (
            &[
                "-i",
                "-k",
                "-c",
                "-t",
                "-o",
                "-s",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &[],
        ),
        "spss-update" => (&["-i", "--add", "--remove", "-t", "-o"], &[]),
        "build-spss" => (
            &[
                "-i",
                "-k",
                "-m",
                "-c",
                "-t",
                "-o",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &["--compact"],
        ),
        "build-fm" => (
            &[
                "-i",
                "-k",
                "-c",
                "-t",
                "-o",
                "--max-mem",
                "--tmp-dir",
                "--graph",
            ],
            &[],
        ),
        "query" => (&["-i", "-q", "-k", "-t", "--ids"], &["--no-verify"]),
        _ => return None,
    };

    Some(flags)
}

/// Help for one subcommand, made of a one-line summary, the synopsis, its FLAGS block and some examples.
pub(super) struct CommandHelp {
    pub(super) name: &'static str,
    pub(super) summary: &'static str,
    pub(super) synopsis: &'static str,
    pub(super) flags: String,
    pub(super) examples: &'static str,
}

/// The help of every subcommand, in the order the overview lists them.
pub(super) fn command_helps() -> Vec<CommandHelp> {
    // Defaults and limits come from the constants the commands use, so the help cannot drift from them.
    const PERCENT: f64 = 100.0;
    let max_k = crate::MAX_K;
    let dbg_max_k = crate::DBG_MAX_K;
    let default_k = crate::DEFAULT_K as usize;
    let map_k = crate::DEFAULT_MAP_K as usize;
    let noisy_k = crate::map::NOISY_K;
    let noisy_pct = crate::map::NOISY_ERROR * PERCENT;
    let window = crate::DEFAULT_WINDOW as usize;
    let max_window = crate::MAX_WINDOW as usize;
    let auto_ed_pct = crate::map::AUTO_ED_RATE * PERCENT;
    let min_auto_ed = crate::map::MIN_AUTO_ED;
    let min_read_len = map_k + window - 1;
    let max_qname = crate::map::MAX_QNAME;
    let dbg_min_count = crate::DEFAULT_DBG_MIN_COUNT as u64;
    let fm_min_count = crate::DEFAULT_FM_MIN_COUNT as u64;

    vec![
        CommandHelp {
            name: "map",
            summary: "Map short, paired-end or long reads to a reference, writing SAM",
            synopsis: "optik map -r REF.fa -q READS.fq [-q MATES.fq] [options] [-o OUT.sam]",
            flags: format!(
                "FLAGS (map)
  -r FILE    Reference genome in FASTA format
  -q FILE    Sequencing reads in FASTA or FASTQ format; give -q twice for
               paired-end reads (-q R1.fq -q R2.fq, same order, FR library):
               mate rescue, proper-pair flags, mate fields, TLEN, MC:Z
  -k INT     Seed (minimizer) length for map  [default: chosen from the
               reads: {map_k}, or {noisy_k} for noisy long reads (over ~{noisy_pct:.0}% error, from
               base qualities, or from the share of read minimizers found
               in the reference for FASTA reads); -k fixes it]
  -w INT     Sliding window for minimizer selection (1..{max_window}, 0 = auto-tune
               from the read length)  [default: {window}]
               Smaller window = more minimizer seeds, higher sensitivity
               Reads shorter than k+w-1 bases produce NO minimizers (all unmapped)
  -e INT     Maximum edit distance for short-read alignment (>= 0)
               [default: auto = {auto_ed_pct:.0}% of the read length, at least {min_auto_ed}]
               Higher = more tolerant of mismatches and indels
  -t INT     Number of threads for parallel mapping (0 = auto)  [default: half of logical CPUs]
  -o FILE    Write SAM to FILE  [default: standard output]
  -R STR     Read group line, e.g. '@RG\\tID:s1\\tSM:sample1' (as bwa and
               minimap2): written to the header, and RG:Z:<ID> added to
               every record (GATK / Picard require read groups)
  --cs       Also write minimap2's cs tag (short form, `:n` runs); --cs=long
               spells identical bases (`=ACGT`).  Every record includes
               NM and MD (SAM spec) either way.

NOTES (map)
  Reads must be at least k + w - 1 bases long for minimizer extraction.
  E.g., with the default k={map_k} w={window} each read needs >={min_read_len} bases; shorter reads
  always appear unmapped. Pre-extracted k-mers (single k-mer per
  sequence) should use the `query` subcommand instead.
  RNA U is read as T.  Reference names must be unique (SAM @SQ); read
  names over {max_qname} characters are truncated (SAM QNAME limit)."
            ),
            examples: "  optik map -r genome.fa -q reads.fq -o reads.sam                 # single-end or long reads\n  optik map -r genome.fa -q R1.fq.gz -q R2.fq.gz -o pairs.sam     # paired-end\n  optik map -r genome.fa -q R1.fq -q R2.fq -R '@RG\\tID:s1\\tSM:s1' --cs -o pairs.sam",
        },
        CommandHelp {
            name: "dbruijn",
            summary: "Count k-mers and write the compacted de Bruijn graph (GFA with links, or unitig FASTA)",
            synopsis: "optik dbruijn -i SEQS.fa|.fq[.gz] [-i ...] [-k 31] [-c 2] [-o graph.gfa]",
            flags: format!(
                "FLAGS (dbruijn)
  -i FILE    Input sequencing reads in FASTA or FASTQ format (.gz, - for
               stdin); repeat for several files, e.g. -i R1.fq -i R2.fq
  -k INT     k-mer size for De Bruijn graph (1..{dbg_max_k})  [default: {default_k}]
               Larger k = more specific graph spurs
  -c INT     Minimum k-mer abundance to include (0 acts as 1)  [default: {dbg_min_count}]
               Higher = filters more sequencing errors
  -t INT     Number of threads for counting + unitig walk (0 = auto)
               [default: half of logical CPUs]
  -o FILE    Output  [default: graph.gfa]: GFA1 (segments + links), or
               unitig FASTA if FILE ends in .fa/.fasta/.fna (BCALM-style
               headers `>id LN:i KC:i km:f`, what SSHash/SBWT/FMSI index)
  -s FILE    Also serialize the graph core to a kdbg archive (not with
               --graph part, which keeps no node table)
  --max-mem SIZE  Memory limit for k-mer counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   Spill directory  [default: the output file's directory]
  --graph MODE    How the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= {dbg_max_k}) or table.  The
               output is identical either way"
            ),
            examples: "  optik dbruijn -i genome.fa -k 31 -c 1 -o graph.gfa          # a genome: keep every k-mer\n  optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa         # unitig FASTA instead of GFA\n  optik dbruijn -i R1.fq.gz -i R2.fq.gz -k 31 -o graph.gfa    # reads: -c 2 drops error k-mers",
        },
        CommandHelp {
            name: "build-spss",
            summary: "Build an exact SPSS k-mer index (membership and dense ids)",
            synopsis: "optik build-spss -i SEQS.fa [-i ...] [-k 31] [--compact [-c N]] -o INDEX.spss",
            flags: format!(
                "FLAGS (build-spss)
  -i FILE    FASTA/FASTQ to index (an SPSS such as unitigs is size-optimal;
               any input gives exact membership); repeatable
  -k INT     k-mer length (1..{max_k})  [default: {default_k}]
  -m INT     Minimizer length (< k)  [default: 0 = auto: chosen from the
               input by an estimate of index size and lookup work per m;
               the index is exact at every m]
  --compact  Count k-mers (>= -c) and index their maximal unitigs instead of
               the raw input sequences
  -c INT     Minimum abundance with --compact  [default: {fm_min_count}]
  -t INT     Threads for the build, and for counting and compaction with
               --compact (0 = auto)  [default: half of logical CPUs]
  -o FILE    Output index file (required)
  --max-mem SIZE  --compact: memory limit for counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   --compact: spill directory  [default: the output file's directory]
  --graph MODE    --compact: how the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= {dbg_max_k}) or table.  The
               output is identical either way"
            ),
            examples: "  optik build-spss -i unitigs.fa -k 31 -o genome.spss\n  optik build-spss --compact -i R1.fq.gz -i R2.fq.gz -c 2 -o reads.spss   # count, compact, index",
        },
        CommandHelp {
            name: "spss-update",
            summary: "Add or remove the k-mers of files in an SPSS index and write it re-compacted",
            synopsis: "optik spss-update -i INDEX.spss [--add FILE ...] [--remove FILE ...] -o OUT.spss",
            flags: "FLAGS (spss-update)
  -i FILE        SPSS index to edit
  --add FILE     Insert every k-mer of FILE (repeatable)
  --remove FILE  Erase every k-mer of FILE (repeatable; applied after adds)
  -t INT         Threads for compaction (0 = auto)
                   [default: half of logical CPUs]
  -o FILE        Write the compacted index (may equal -i)"
                .to_string(),
            examples: "  optik spss-update -i genome.spss --add new.fa --remove old.fa -o genome.spss",
        },
        CommandHelp {
            name: "build-fm",
            summary: "Build an exact FM k-mer index of the unitigs (the smallest index)",
            synopsis: "optik build-fm -i SEQS.fa [-i ...] [-k 31] [-c 1] -o INDEX.fm",
            flags: format!(
                "FLAGS (build-fm)
  -i FILE    Reference FASTA/FASTQ to build the FM index from; repeatable
  -k INT     k-mer length (1..{max_k})  [default: {default_k}]
  -c INT     Minimum k-mer abundance to include  [default: {fm_min_count}]
               Default keeps every distinct k-mer (exact membership)
  -t INT     Threads for counting, compaction and the index build (0 = auto)
               [default: half of logical CPUs]
  -o FILE    Write serialized index to file (required)
  --max-mem SIZE  Memory limit for k-mer counting, e.g. 8G or 500M  [default:
               half of the available memory, cgroup limits included].  Past
               it the counter's partition buffers spill to disk (same
               result); the graph itself must fit
  --tmp-dir DIR   Spill directory  [default: the output file's directory]
  --graph MODE    How the unitigs are built: auto (default: partitioned from
               k = 19 up, the k-mer table below), part (partitioned
               compaction: fragments compacted per counted partition and
               glued, ~1/4 of the memory; any k <= {dbg_max_k}) or table.  The
               output is identical either way"
            ),
            examples: "  optik build-fm -i genome.fa -k 31 -o genome.fm",
        },
        CommandHelp {
            name: "query",
            summary: "Look up every k-mer of a FASTA/FASTQ in an SPSS or FM index",
            synopsis: "optik query -i INDEX.spss|INDEX.fm -q QUERIES.fa|.fq[.gz] [-t N] [--ids IDS.txt]",
            flags: format!(
                "FLAGS (query)
  -i FILE    Saved SPSS or FM index (build-spss / build-fm)
  -q FILE    FASTA/FASTQ whose k-mers are looked up (k-mers or reads; .gz ok)
  -k INT     k-mer length (1..{max_k}); must equal the index's  [default: the index's k]
  -t INT     Threads for the query phase (1 = single-threaded,
               >1 = parallel over records, long records split into
               window-exact slices; counts are identical at any -t;
               0 = auto)  [default: half of logical CPUs]
  --ids FILE SPSS index only: also write, per query record, the dense
               identifier in [0, n) of every k-mer window (-1 = absent),
               space-separated, one line per record (SSHash Lookup)
  --no-verify  SPSS index only: skip the checksum check at load (by default
               the whole index is checksummed once, so a truncated or
               corrupted file is an error instead of wrong answers)"
            ),
            examples: "  optik query -i genome.spss -q reads.fq                  # counts of present and absent k-mers\n  optik query -i genome.spss -q reads.fq --ids ids.txt    # also one id per k-mer window",
        },
    ]
}

/// The overview `optik` prints when no command or an unknown one is given:
/// What each command does, with one example each and where to read more.
pub(super) fn overview_text() -> String {
    let mut text = String::from(
        "optik -- genomic sequence analysis toolkit

USAGE
  optik <command> [options]
  optik <command> --help      that command's flags and examples
  optik --help                every command and flag
  optik tui                   interactive terminal UI (also plain `optik` in a terminal)
  optik --version

COMMANDS
",
    );

    // Command names are padded to one column, so the summaries line up.
    const COMMAND_NAME_W: usize = 12;

    for help in command_helps() {
        text.push_str(&format!("  {:<COMMAND_NAME_W$}{}\n", help.name, help.summary));
    }

    text.push_str("  build       Alias of build-spss\n  tui         Interactive terminal UI over every command\n\nQUICK START\n");
    text.push_str(
        "  optik dbruijn -i genome.fa -k 31 -c 1 -o unitigs.fa   # graph of a genome, as unitigs
  optik build-spss -i unitigs.fa -o genome.spss          # exact k-mer index of the unitigs
  optik query -i genome.spss -q reads.fq                 # look up every k-mer of the reads
  optik map -r genome.fa -q reads.fq -o reads.sam        # map the reads to SAM
",
    );
    text.push_str("\nRun `optik <command> --help` for a command's flags and more examples.\n");
    text
}

/// The help of subcommand `sub`, or `None` for an unknown one.
pub(super) fn command_help_text(sub: &str) -> Option<String> {
    let help = command_helps().into_iter().find(|help| help.name == sub)?;

    Some(format!(
        "optik {} -- {}\n\nUSAGE\n  {}\n\n{}\n\nEXAMPLES\n{}\n",
        help.name, help.summary, help.synopsis, help.flags, help.examples
    ))
}

/// Prints the help text on stderr.
pub(super) fn usage() {
    eprintln!("{}", usage_text());
}

/// The full help, the overview and then every command's flags.
pub(super) fn usage_text() -> String {
    let mut text = overview_text();

    text.push('\n');

    for help in command_helps() {
        text.push('\n');
        text.push_str(&help.flags);
        text.push('\n');
    }

    text
}

/// The nearest known command to `word` by edit distance, for "did you mean".
pub(super) fn closest_command(word: &str) -> Option<&'static str> {
    fn dist(left: &str, right: &str) -> usize {
        let (left_chars, right_chars): (Vec<char>, Vec<char>) =
            (left.chars().collect(), right.chars().collect());
        let mut row: Vec<usize> = (0..=right_chars.len()).collect();

        for i in 1..=left_chars.len() {
            let mut prev = row[0];

            row[0] = i;

            for j in 1..=right_chars.len() {
                let cur = row[j];

                row[j] = (row[j] + 1)
                    .min(row[j - 1] + 1)
                    .min(prev + usize::from(left_chars[i - 1] != right_chars[j - 1]));
                prev = cur;
            }
        }

        row[right_chars.len()]
    }
    const NAMES: [&str; 8] = [
        "map",
        "dbruijn",
        "build-spss",
        "spss-update",
        "build-fm",
        "query",
        "build",
        "tui",
    ];
    /// A typo matches a subcommand within this many edits, or a third of the name's length if that is more.
    const MAX_TYPO_EDITS: usize = 2;
    const NAME_CHARS_PER_TYPO: usize = 3;
    NAMES
        .iter()
        .map(|name| (dist(word, name), *name))
        .filter(|(distance, name)| {
            *distance <= MAX_TYPO_EDITS.max(name.len() / NAME_CHARS_PER_TYPO)
        })
        .min()
        .map(|(_, name)| name)
}

/// Reports a command-line mistake in `sub` with the command's synopsis and a pointer to its full help, then exits with status 1.
pub(super) fn arg_error(sub: &str, msg: &str) -> ! {
    eprintln!("error: {msg}");

    if let Some(help) = command_helps().into_iter().find(|help| help.name == sub) {
        eprintln!("usage: {}", help.synopsis);
    }

    eprintln!("see: optik {sub} --help");
    std::process::exit(1)
}

#[cfg(test)]
mod help_tests {
    const COMMANDS: [&str; 6] = [
        "map",
        "dbruijn",
        "build-spss",
        "build-fm",
        "query",
        "spss-update",
    ];

    /// Every flag a command accepts starts a line of its help, and the help ends with examples.
    #[test]
    fn each_command_help_lists_its_flags() {
        for sub in COMMANDS {
            let help = super::command_help_text(sub).unwrap_or_else(|| panic!("no help for {sub}"));
            let (vals, bools) = super::subcommand_flags(sub).unwrap();

            for flag in vals.iter().chain(bools) {
                let name = flag.split('=').next().unwrap();

                assert!(
                    help.lines().any(|line| line.trim_start().starts_with(name)),
                    "{sub} --help does not document {flag}"
                );
            }

            assert!(
                help.contains("EXAMPLES") && help.contains(&format!("optik {sub} ")),
                "{sub}"
            );
        }
    }

    /// The overview names every command, and the full help has each command's flag block.
    #[test]
    fn overview_names_every_command_and_full_help_has_every_block() {
        let overview = super::overview_text();
        let full = super::usage_text();

        for sub in COMMANDS {
            assert!(
                overview.contains(&format!("  {sub} ")),
                "overview lacks {sub}"
            );
            assert!(
                full.contains(&format!("FLAGS ({sub})")),
                "--help lacks {sub}"
            );
        }
    }

    /// Near misses get the right suggestion, while a word far from every command gets none.
    #[test]
    fn typos_suggest_the_command() {
        assert_eq!(super::closest_command("qeury"), Some("query"));
        assert_eq!(super::closest_command("dbrujin"), Some("dbruijn"));
        assert_eq!(super::closest_command("buildfm"), Some("build-fm"));
        assert_eq!(super::closest_command("zzzzzzzz"), None);
    }
}
