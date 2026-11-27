//! optik, a genomic sequence analysis toolkit.
//!
//! This file holds the command line, which checks arguments and dispatches, the C ABI declarations of the C++ core, and the safe Rust wrappers every subcommand uses.
//! The C++ core in include/ and src/ does the k-mer counting, minimizer selection, the bit-parallel and wavefront aligners, the de Bruijn graph and the SPSS and FM indexes.
//!
//! Subcommands
//!   tui        Interactive terminal UI
//!   map        Read mapping via minimizer seeds
//!   dbruijn    De Bruijn graph assembly
//!   build-spss Build the recommended SPSS membership index, also called `build`
//!   spss-update Add or remove k-mers of an SPSS index
//!   build-fm   Build the FM membership index over a k-order BWT
//!   query      Query an SPSS or FM index

mod chain;
mod map;
mod seqio;

mod build;
mod dbruijn;
mod query;
mod tui;

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::{c_char, c_int, c_void};

mod ffi;
mod native;
mod util;
mod help;

use ffi::*;
use native::*;
use util::*;
use help::*;

/// Maximum supported k-mer size, since k-mers are packed 2 bits per base into a 64-bit word.
pub(crate) const MAX_K: usize = 32;
/// dbruijn's largest k, since its (k+1)-mer edge set must fit 2 bits per base in 64 bits.
pub(crate) const DBG_MAX_K: usize = MAX_K - 1;
/// The seed index stores 31-bit reference positions, so no reference sequence may reach this length.
pub(crate) const MAX_REF_SEQ_LEN: usize = 1 << 31;

/// Command-line defaults and limits.
/// Default -k of map, whose seeds survive 1% read error better at 21 than at 31.
const DEFAULT_MAP_K: f64 = 21.0;
/// Default -k of every other subcommand.
const DEFAULT_K: f64 = 31.0;
/// Default -w, the minimizer window, and its largest value, the size of the window ring.
const DEFAULT_WINDOW: f64 = 10.0;
const MAX_WINDOW: f64 = 256.0;
/// map's -e when absent or "auto": the edit cap is chosen per read.
const AUTO_EDIT_CAP: i32 = -1;
/// Largest explicit -e.
const MAX_EDIT_CAP: f64 = 1_000_000.0;
/// Default -c of dbruijn (k-mers seen once are dropped as read errors) and of build-fm / build-spss --compact (every distinct k-mer is kept).
const DEFAULT_DBG_MIN_COUNT: f64 = 2.0;
const DEFAULT_FM_MIN_COUNT: f64 = 1.0;
/// The largest -t accepted, while -t 0 (or no -t) uses the default thread budget.
const AUTO_THREADS: usize = 0;
const MAX_THREADS: f64 = 4096.0;
/// The largest -m accepted, while -m 0 (or no -m) lets build-spss choose the minimizer length.
const AUTO_MINIMIZER_LEN: f64 = 0.0;
const MAX_MINIMIZER_LEN: f64 = 31.0;
/// Without --max-mem, the counter picks its memory limit itself.
const AUTO_MAX_MEM: u64 = 0;
/// --graph modes passed to ok_set_count_limits (include/spill.hpp CountLimits::graph_mode).
const GRAPH_MODE_AUTO: c_int = 0;
const GRAPH_MODE_PARTITIONED: c_int = 1;
const GRAPH_MODE_TABLE: c_int = 2;
/// Size suffixes K, M, G and T of --max-mem, in bytes (binary units).
const KIB: u64 = 1 << 10;
const MIB: u64 = 1 << 20;
const GIB: u64 = 1 << 30;
const TIB: u64 = 1 << 40;
/// KB per MB in memory reports.
const KB_PER_MB: f64 = 1024.0;

/// Arguments a value flag spans, the flag itself and the value after it.
const FLAG_AND_VALUE: usize = 2;

/// Every value given to `flag` in an argument list of -flag value pairs, in order of appearance.
fn flag_values<'a>(rest: &'a [String], flag: &str) -> impl Iterator<Item = &'a str> {
    rest.windows(FLAG_AND_VALUE).filter_map(move |pair| match pair {
        [name, value] if name == flag => Some(value.as_str()),
        _ => None,
    })
}

/// Runs the terminal UI, exiting with status 1 if it fails.
fn start_tui() {
    if let Err(err) = tui::run() {
        eprintln!("TUI error: {err}");
        std::process::exit(1);
    }
}

fn main() {
    // Output written to a closed pipe, as in `optik query ... | head`, ends the process quietly, as for C tools.
    // Rust ignores SIGPIPE by default and println! then panics, which with panic = "abort" dumped core.
    // Any other failure to write stdout, such as a full disk, is an error rather than a crash.
    // SAFETY: this runs first in main, before any other thread exists, and SIG_DFL installs no Rust handler.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    let default_hook = std::panic::take_hook();

    std::panic::set_hook(Box::new(move |info| {
        let msg = info
            .payload()
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| info.payload().downcast_ref::<&str>().copied())
            .unwrap_or("");

        if msg.starts_with("failed printing to stdout") {
            eprintln!(
                "error: writing to stdout failed: {}",
                msg.trim_start_matches("failed printing to stdout: ")
            );
            std::process::exit(1);
        }

        default_hook(info);
    }));
    // glibc raises its mmap threshold, up to 32 MB, each time a large mmap'd block is freed.
    // Later large buffers, such as per-thread counting and walk arrays and query batches, then land on arena heaps and stay resident after they are freed, adding about 100 MB to the peak on a single-record genome.
    // A fixed 1 MiB threshold keeps large buffers mmap'd and returned to the OS on free, at the same measured speed.
    // Only glibc does this.
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: mallopt only changes a glibc tuning value, and it runs before any thread allocates.
    unsafe {
        const MMAP_THRESHOLD_BYTES: c_int = 1 << 20;
        libc::mallopt(libc::M_MMAP_THRESHOLD, MMAP_THRESHOLD_BYTES);
    }
    let args: Vec<String> = std::env::args().collect();
    // The program name comes first, then the command and its own arguments.
    let (command, rest) = match args.as_slice() {
        [_, command, rest @ ..] => (command.as_str(), rest),
        _ => {
            use std::io::IsTerminal;

            // Without arguments and without a terminal (scripts, pipes) the usage is printed instead of failing to start the TUI.
            if !std::io::stdout().is_terminal() || !std::io::stdin().is_terminal() {
                eprintln!("{}", overview_text());
                std::process::exit(1);
            }

            start_tui();
            return;
        }
    };

    // "tui" as the command starts the terminal UI too.
    if command == "tui" {
        start_tui();
        return;
    }

    // Help flags redirect to the usage text.
    if command == "-V" || command == "--version" {
        println!("optik {}", env!("CARGO_PKG_VERSION"));
        return;
    }

    if command == "-h" || command == "--help" || command == "help" {
        // `optik help map` is `optik map --help`.
        match rest
            .first()
            .and_then(|cmd| command_help_text(if cmd == "build" { "build-spss" } else { cmd }))
        {
            Some(help_text) => eprintln!("{help_text}"),
            None => usage(),
        }

        return;
    }

    // `build` is an alias of build-spss.
    let sub: &str = if command == "build" { "build-spss" } else { command };

    // `optik <command> --help` prints that command's flags and examples.
    if rest.iter().any(|arg| arg == "-h" || arg == "--help") {
        match command_help_text(sub) {
            Some(help_text) => eprintln!("{help_text}"),
            None => usage(),
        }

        return;
    }

    // A command given no arguments at all prints its help, since every command needs an input.
    if rest.is_empty()
        && let Some(help_text) = command_help_text(sub)
    {
        eprintln!("{help_text}");
        std::process::exit(1);
    }

    // Every argument must be a flag this subcommand knows (a typo such as "-K 25" or "--thread 4" would otherwise be ignored silently and the run would use a default), and every value flag needs its value.
    {
        let Some((vals, bools)) = subcommand_flags(sub) else {
            eprintln!("error: unknown command '{sub}'");

            if let Some(suggestion) = closest_command(sub) {
                eprintln!("did you mean: optik {suggestion}");
            }

            eprintln!("\n{}", overview_text());
            std::process::exit(1);
        };
        // A value flag given twice is an error, since the second would be ignored silently.
        // Only spss-update --add/--remove and -i of the k-mer counting commands take several values.
        // The latter accepts several input files, such as paired-end R1 and R2.
        // Repeating any other flag is an error, because the last value would otherwise win silently.
        let multi = |flag: &str| {
            flag == "--add"
                || flag == "--remove"
                || (flag == "-i" && matches!(sub, "dbruijn" | "build-spss" | "build-fm"))
                || (flag == "-q" && sub == "map")
        };
        let mut seen: Vec<&str> = Vec::new();
        let mut i = 0;

        while i < rest.len() {
            let arg = rest[i].as_str();

            if vals.contains(&arg) {
                if i + 1 >= rest.len() {
                    arg_error(sub, &format!("{arg} needs a value"));
                }

                if seen.contains(&arg) && !multi(arg) {
                    arg_error(sub, &format!("{arg} given more than once"));
                }

                seen.push(arg);
                i += FLAG_AND_VALUE;
            } else if bools.contains(&arg) {
                i += 1;
            } else {
                arg_error(sub, &format!("unknown argument '{arg}' for {sub}"));
            }
        }
    }

    // Flag lookup over the -flag value pairs of the argument list.
    let get = |flag: &str| -> Option<&str> { flag_values(rest, flag).next() };

    fn required<'a>(
        sub: &str,
        get: &dyn Fn(&str) -> Option<&'a str>,
        flag: &str,
        msg: &str,
    ) -> &'a str {
        get(flag).unwrap_or_else(|| arg_error(sub, msg))
    }
    // Strict numeric flags.
    // A flag that is present but malformed or out of its type's range is an error, never a silent fall-back to the default.
    let num = |flag: &str, default: f64| -> f64 {
        match get(flag) {
            None => default,
            Some(raw) => match raw.parse::<f64>() {
                Ok(value) if value.is_finite() && value == value.trunc() => value,
                _ => {
                    eprintln!("error: {flag} expects an integer, got '{raw}'");
                    std::process::exit(1);
                }
            },
        }
    };
    let in_range = |flag: &str, value: f64, min: f64, max: f64| -> f64 {
        if value < min || value > max {
            eprintln!("error: {flag} {value} out of range [{min}, {max}]");
            std::process::exit(1);
        }

        value
    };
    // map seeds default to k=21, since at 1% error most 31-mers of a 150 bp read contain an error (P(error-free 31-mer) = 0.73 vs 0.81), and reads with 5-6 spread edits keep no 31-mer seed at all.
    let k_raw = num(
        "-k",
        if sub == "map" {
            DEFAULT_MAP_K
        } else {
            DEFAULT_K
        },
    );

    if k_raw < 1.0 || k_raw > MAX_K as f64 {
        eprintln!(
            "error: k={k_raw} is out of range [1, {MAX_K}] (k-mers are 2-bit packed in 64 bits)"
        );
        std::process::exit(1);
    }

    let kmer_len = k_raw as usize;
    // -w 0 auto-tunes, and the minimizer window ring holds at most 256 k-mers.
    let window = in_range("-w", num("-w", DEFAULT_WINDOW), 0.0, MAX_WINDOW) as usize;
    // map -e is an explicit fixed cap (>= 0) or, by default, auto at 15% of each read's length.
    let max_edits: i32 = match get("-e") {
        None | Some("auto") => AUTO_EDIT_CAP,
        Some(_) => in_range("-e", num("-e", 0.0), 0.0, MAX_EDIT_CAP) as i32,
    };
    let min_count = in_range("-c", num("-c", DEFAULT_DBG_MIN_COUNT), 0.0, u32::MAX as f64) as u32;
    // build-fm and build-spss --compact take the same -c, whose default of 1 keeps every distinct k-mer for exact membership.
    let fm_min_count = in_range("-c", num("-c", DEFAULT_FM_MIN_COUNT), 0.0, u32::MAX as f64) as u64;
    let threads: usize = match get("-t") {
        None => default_threads(),
        Some(_) => {
            match in_range("-t", num("-t", AUTO_THREADS as f64), 0.0, MAX_THREADS) as usize {
                AUTO_THREADS => default_threads(),
                requested => requested,
            }
        }
    };

    // OpenMP regions in the C++ core follow -t too, though they default to all cores.
    // omp_set_num_threads covers regions entered from this thread, and the environment set by omp_env_for_budget covers every region's default and the idle policy.
    omp_env_for_budget(threads.max(1));
    // SAFETY: the call takes a plain integer and runs before any OpenMP region starts.
    unsafe { ok_set_num_threads(threads.max(1) as u32) };

    let m_val = in_range("-m", num("-m", AUTO_MINIMIZER_LEN), 0.0, MAX_MINIMIZER_LEN) as usize;

    // k-mer counting memory limit and spill directory (dbruijn, build-fm, build-spss --compact).
    // Past the limit the counter's partition buffers go to disk, by default in the output's directory rather than /tmp, which is often RAM-backed (tmpfs).
    if matches!(sub, "dbruijn" | "build-fm" | "build-spss") {
        let max_mem = match get("--max-mem") {
            None => AUTO_MAX_MEM,
            Some(raw) => parse_size(raw).unwrap_or_else(|| {
                eprintln!("error: --max-mem expects a size such as 8G, 512M or 4.5G, got '{raw}'");
                std::process::exit(1);
            }),
        };

        const MIN_MAX_MEM: u64 = 64 << 20;

        if get("--max-mem").is_some() && max_mem < MIN_MAX_MEM {
            eprintln!("error: --max-mem must be at least 64M");
            std::process::exit(1);
        }

        let dir = match get("--tmp-dir") {
            Some(tmp_dir) => tmp_dir.to_string(),
            None => {
                // Spill files default to the output's directory, which is likely on a disk with room for results.
                let out = get("-o").unwrap_or(if sub == "dbruijn" { "graph.gfa" } else { "." });

                match std::path::Path::new(out).parent() {
                    Some(parent) if !parent.as_os_str().is_empty() => {
                        parent.to_string_lossy().into_owned()
                    }
                    _ => ".".to_string(),
                }
            }
        };

        if !std::path::Path::new(&dir).is_dir() {
            eprintln!("error: spill directory '{dir}' does not exist (--tmp-dir)");
            std::process::exit(1);
        }

        let cdir = CString::new(dir).unwrap_or_else(|_| {
            eprintln!("error: --tmp-dir contains a NUL byte");
            std::process::exit(1);
        });
        // --graph chooses how the unitigs are built (by default partitioned from k = 19 up and the table below that), with identical output either way.
        let graph_mode: c_int = match get("--graph") {
            None | Some("auto") => GRAPH_MODE_AUTO,
            Some("part") | Some("partitioned") => GRAPH_MODE_PARTITIONED,
            Some("table") => GRAPH_MODE_TABLE,
            Some(other) => {
                eprintln!("error: --graph expects auto, part or table, got '{other}'");
                std::process::exit(1);
            }
        };

        if graph_mode == GRAPH_MODE_PARTITIONED && get("-s").is_some() {
            eprintln!("error: --graph part keeps no node table, so it cannot write the -s archive");
            std::process::exit(1);
        }

        if graph_mode == GRAPH_MODE_PARTITIONED && kmer_len > DBG_MAX_K {
            eprintln!("error: --graph part needs k <= 31");
            std::process::exit(1);
        }

        // SAFETY: cdir is a NUL-terminated CString, and the C++ side copies it into its own std::string before returning.
        unsafe { ok_set_count_limits(max_mem, cdir.as_ptr(), graph_mode) };
    }

    // An output path must not name an input, since creating the output truncates it before it is read.
    // spss-update -i X -o X is the supported in-place update, because the index is fully loaded before it is rewritten.
    {
        let outs: Vec<&str> = ["-o", "--ids", "-s"]
            .iter()
            .filter_map(|flag| get(flag))
            .collect();
        let mut ins: Vec<(&str, &str)> = Vec::new();

        for flag in ["-i", "-q", "-r", "--add", "--remove"] {
            for value in flag_values(rest, flag) {
                if sub == "spss-update" && flag == "-i" {
                    continue;
                }

                ins.push((flag, value));
            }
        }

        // Two outputs on one file would let the second overwrite the first.
        for (out_idx, out_a) in outs.iter().enumerate() {
            for out_b in &outs[out_idx + 1..] {
                if *out_a != "/dev/null" && (out_a == out_b || same_file(out_a, out_b)) {
                    eprintln!("error: two outputs name the same file {out_a}");
                    std::process::exit(1);
                }
            }
        }

        for out_path in &outs {
            for (flag, in_path) in &ins {
                if *in_path != "-" && same_file(out_path, in_path) {
                    eprintln!(
                        "error: output {out_path} is the {flag} input {in_path}; refusing to overwrite it"
                    );
                    std::process::exit(1);
                }
            }
        }
    }

    // Every -i of a counting command, joined by '\n' for the C++ reader to split, so a path containing a newline cannot be passed this way.
    let inputs: String = {
        let input_paths: Vec<&str> = flag_values(rest, "-i").collect();

        if input_paths.iter().any(|path| path.contains('\n')) {
            eprintln!("error: input paths may not contain a newline");
            std::process::exit(1);
        }

        if input_paths.iter().filter(|path| **path == "-").count() > 1 {
            eprintln!("error: stdin (-) can be given as -i only once");
            std::process::exit(1);
        }

        input_paths.join("\n")
    };
    let need_i = |msg: &str| -> &str {
        if inputs.is_empty() {
            arg_error(sub, msg)
        }

        inputs.as_str()
    };

    let result = match sub {
"map"     => {
            // -q R1 -q R2 means paired-end, as minimap2 and BWA take two read files.
            // One read file is single-end and two are paired-end R1 and R2.
            const PAIRED_READ_FILES: usize = 2;
            let read_files: Vec<&str> = flag_values(rest, "-q").collect();

            if read_files.len() > PAIRED_READ_FILES {
                eprintln!("error: map takes one -q (single-end) or two (paired-end R1 R2)");
                std::process::exit(1);
            }

            // R1 and R2 cannot both be read from a single stdin stream.
            if read_files.len() == PAIRED_READ_FILES && read_files.iter().all(|path| *path == "-") {
                eprintln!("error: stdin (-) can be only one of the two read files");
                std::process::exit(1);
            }

            map::run(required(sub, &get, "-r", "-r required"),
                               required(sub, &get, "-q", "-q required"), read_files.get(1).copied(), kmer_len,
                               get("-k").is_none(), window, max_edits,
                               get("-o"), threads,
                               if rest.iter().any(|arg| arg == "--cs=long") { map::CsMode::Long }
                               else if rest.iter().any(|arg| arg == "--cs" || arg == "--cs=short") { map::CsMode::Short }
                               else { map::CsMode::Off },

                               get("-R"), &args.join(" "))
        }
        "dbruijn" => dbruijn::run(need_i("-i required"), kmer_len, min_count,
                                  get("-o").unwrap_or("graph.gfa"),
                                  get("-s"), threads),
        "spss-update" => {
            let all = |flag: &str| -> Vec<&str> { flag_values(rest, flag).collect() };

            build::update_spss(required(sub, &get, "-i", "-i <index> required"),
                               &all("--add"), &all("--remove"), threads,
                               required(sub, &get, "-o", "-o <index> required"))
        }
        // -c filters counted k-mers, which only --compact does.
        // Without it the records are indexed as given and -c would be ignored silently.
        "build-spss" if get("-c").is_some() && !rest.iter().any(|arg| arg == "--compact") => Err(
            "-c (minimum abundance) applies only with --compact; without it every k-mer of the \
             input is indexed as given".into()),
        // The counter's memory limit likewise applies only to --compact, which counts k-mers.
        "build-spss" if (get("--max-mem").is_some() || get("--tmp-dir").is_some() || get("--graph").is_some())
                        && !rest.iter().any(|arg| arg == "--compact") => Err(
            "--max-mem / --tmp-dir / --graph apply only with --compact (the k-mer counter); without it the \
             records are indexed as given".into()),
        "build-spss" => build::build_spss(need_i("-i <ref.fa> required"),
                                          kmer_len, m_val, rest.iter().any(|arg| arg == "--compact"),
                                          fm_min_count, threads,
                                          required(sub, &get, "-o", "-o <index> required")),
        "build-fm"  => build::build_fm(need_i("-i <ref.fa> required"),
                                       kmer_len, fm_min_count, threads as u32,
                                       required(sub, &get, "-o", "-o <index.bin> required")),
        "query"     => query::query(required(sub, &get, "-i", "-i <index.bin> required"),
                                     required(sub, &get, "-q", "-q required"),
                                     get("-k").map(|_| kmer_len), threads,
                                     get("--ids"), !rest.iter().any(|arg| arg == "--no-verify")),
        _ => { usage(); std::process::exit(1); }
    };

    if let Err(err) = result {
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod correctness_tests {
    use std::collections::HashSet;

    const TEST_KMER_LEN: usize = 31;

    /// Reverse complement with the C++ encoding (A=0, C=1, G=2, T=3, so the complement is XOR 3).
    fn revcomp(packed: u64, kmer_len: usize) -> u64 {
        let mut rev = 0u64;

        for i in 0..kmer_len {
            let base = (packed >> (2 * i)) & 3;

            rev |= (base ^ 3) << (2 * (kmer_len - 1 - i));
        }

        rev
    }

    /// Packs a DNA string with the C++ encoding, without canonicalising it.
    fn pack(seq: &[u8]) -> u64 {
        let mut packed = 0u64;

        for &base in seq {
            packed = (packed << 2)
                | match base {
                    b'A' => 0u64,
                    b'C' => 1,
                    b'G' => 2,
                    _ => 3, // T=3 to match C++ base_encode
                };
        }

        packed
    }

    /// Spells a packed k-mer as bases.
    fn decode(packed: u64, kmer_len: usize) -> Vec<u8> {
        (0..kmer_len)
            .rev()
            .map(|i| b"ACGT"[((packed >> (2 * i)) & 3) as usize])
            .collect()
    }
    /// Membership of one k-mer through the scan path `optik query` uses.
    fn spss_has(idx: &crate::SpssIdx, packed: u64) -> bool {
        idx.scan_chunk(&decode(packed, TEST_KMER_LEN)).1 == 1
    }
    fn fm_has(idx: &crate::FmIdx, packed: u64) -> bool {
        idx.scan_chunk(&decode(packed, TEST_KMER_LEN), TEST_KMER_LEN)
            .1
            == 1
    }

    /// Writes a synthetic tandem-repeat FASTA (220 records, each a random 400 bp unit repeated 5 times) and returns its path and canonical k-mer set.
    fn synthetic_tandem_fasta(path: &std::path::Path, seed: u64) -> (String, HashSet<u64>) {
        let mut fasta = String::new();
        let mut staged: HashSet<u64> = HashSet::new();
        let mut rng = seed;

        for rec in 0..220usize {
            // Repeating one random 400-mer five times per record makes every k-mer recur and exercises duplicate handling in both index builders.
            let mut seed2 = Vec::with_capacity(400);

            for _ in 0..400 {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                seed2.push(match rng & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                });
            }

            let mut seq = Vec::with_capacity(2000);

            for _ in 0..5 {
                seq.extend_from_slice(&seed2);
            }

            for w in 0..=seq.len() - TEST_KMER_LEN {
                let packed = pack(&seq[w..w + TEST_KMER_LEN]);

                staged.insert(packed.min(revcomp(packed, TEST_KMER_LEN)));
            }

            fasta.push_str(&format!(">r{rec}\n{}\n", String::from_utf8(seq).unwrap()));
        }

        std::fs::write(path, &fasta).unwrap();
        (path.to_str().unwrap().to_string(), staged)
    }

    /// --max-mem sizes with K/M/G/T suffixes in powers of 1024, fractions, and the rejection of malformed or non-positive values.
    #[test]
    fn parse_size_units() {
        use crate::parse_size;
        assert_eq!(parse_size("1024"), Some(1024));
        assert_eq!(parse_size("64M"), Some(64 << 20));
        assert_eq!(parse_size("8g"), Some(8 << 30));
        assert_eq!(parse_size("1.5G"), Some(3 << 29));
        assert_eq!(parse_size("2T"), Some(2 << 40));
        assert_eq!(parse_size("8GB"), Some(8 << 30));
        assert_eq!(parse_size("8GiB"), Some(8 << 30));
        assert_eq!(parse_size("500mb"), Some(500 << 20));

        for bad in [
            "", "G", "GB", "-1G", "0", "abc", "1.2.3M", "NaNG", "infG", "8B",
        ] {
            assert_eq!(parse_size(bad), None, "{bad}");
        }
    }

    /// SPSS index (build-spss --compact), with exact membership on every positional k-mer of a tandem-repeat reference and no false positive on random k-mers, before and after a save/load round trip, and a distinct count equal to the k-mer set.
    #[test]
    fn spss_index_is_exact() {
        use crate::SpssIdx;
        let dir = std::env::temp_dir().join(format!("optik_spss_{}", std::process::id()));

        std::fs::create_dir_all(&dir).unwrap();

        let (fa_path, staged) = synthetic_tandem_fasta(&dir.join("sp.fa"), 0x9E3779B97F4A7C15u64);

        assert!(
            staged.len() >= 65_536,
            "test setup: too few distinct k-mers"
        );

        let check = |idx: &SpssIdx, staged: &HashSet<u64>| {
            assert_eq!(idx.kmer_len(), TEST_KMER_LEN);
            assert_eq!(idx.distinct() as usize, staged.len());

            for &kmer in staged {
                assert!(spss_has(idx, kmer), "covered k-mer {:x} not found", kmer);
            }

            let mut neg: Vec<u64> = Vec::with_capacity(2000);
            // Negative k-mers come from a xorshift64 generator, and any that happen to be in the set are skipped.
            let mut rng = 0x123456789ABCDEF0u64;

            while neg.len() < 2000 {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                let packed = (rng & ((1u64 << (2 * TEST_KMER_LEN)) - 1)) | 1;
                let canonical = packed.min(revcomp(packed, TEST_KMER_LEN));

                if !staged.contains(&canonical) {
                    neg.push(canonical);
                }
            }

            for &kmer in &neg {
                assert!(!spss_has(idx, kmer), "false positive {:x}", kmer);
            }
        };

        let idx = SpssIdx::build_from_fasta(&fa_path, TEST_KMER_LEN, TEST_KMER_LEN / 2, true, 1, 0)
            .expect("SPSS build");

        check(&idx, &staged);
        let bin = format!("{}/sp.spss", dir.to_str().unwrap());

        idx.save(&bin).expect("save");
        let idx2 = SpssIdx::load(&bin, true).expect("load");

        check(&idx2, &staged);

        std::fs::remove_dir_all(&dir).ok();
    }

    /// FM index with exact membership.
    /// Every k-mer of a synthetic reference must be found and random non-reference k-mers rejected, before and after a save/load round trip, and the saved file must be exactly bytes() long.
    #[test]
    fn fm_index_is_exact() {
        use crate::FmIdx;
        let dir = std::env::temp_dir().join(format!("optik_fm_{}", std::process::id()));

        std::fs::create_dir_all(&dir).unwrap();

        let (fa_path, staged) = synthetic_tandem_fasta(&dir.join("fm.fa"), 0xABCDEF0123456789u64);

        assert!(
            staged.len() >= 65_536,
            "test setup: too few distinct k-mers"
        );

        let check = |idx: &FmIdx, staged: &HashSet<u64>| {
            assert_eq!(idx.kmer_len(), TEST_KMER_LEN);
            assert_eq!(
                idx.size(),
                staged.len(),
                "indexed {} distinct k-mers, expected {}",
                idx.size(),
                staged.len()
            );
            for &kmer in staged {
                assert!(fm_has(idx, kmer), "covered k-mer {:x} not found", kmer);
            }

            let mut neg: Vec<u64> = Vec::with_capacity(2000);
            // Negative k-mers come from a xorshift64 generator, and any that happen to be in the set are skipped.
            let mut rng = 0xFEDCBA9876543210u64;

            while neg.len() < 2000 {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                let packed = (rng & ((1u64 << (2 * TEST_KMER_LEN)) - 1)) | 1;
                let canonical = packed.min(revcomp(packed, TEST_KMER_LEN));

                if !staged.contains(&canonical) {
                    neg.push(canonical);
                }
            }

            for &kmer in &neg {
                assert!(!fm_has(idx, kmer), "false positive {:x}", kmer);
            }
        };

        let idx = FmIdx::build_from_fasta(&fa_path, TEST_KMER_LEN, 1, 0).expect("FM build");
        let est = idx.bytes() as f64;

        assert!(
            est < idx.size() as f64 * 8.0,
            "FM index {est:.0} B should stay below 8 b/k ({} kmers)",
            idx.size()
        );
        check(&idx, &staged);

        // The round trip's load path re-reads the k-order BWT's rank arrays.
        let bin = format!("{}/fm.bin", dir.to_str().unwrap());

        idx.save(&bin).expect("save");
        assert_eq!(
            std::fs::metadata(&bin)
                .map(|meta| meta.len() as usize)
                .unwrap(),
            idx.bytes(),
            "serialized file size must equal bytes()"
        );
        let idx2 = FmIdx::load(&bin).expect("load");

        check(&idx2, &staged);

        std::fs::remove_dir_all(&dir).ok();
    }
}
