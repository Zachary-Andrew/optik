//! `optik tui` is an interactive terminal front end.
//!
//! The user picks a tool, fills in its parameters with a file browser for paths, and runs it.
//! The TUI starts the same optik binary with the matching subcommand and flags as a child process, streams its output, and shows the result, formatting SAM records and the tab-separated metrics the subcommands print.
//! It is built on ratatui and crossterm.

use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::process::Child;
use std::sync::mpsc::Receiver;
use std::time::Instant;

use crossterm::{
    event::{self, Event, KeyCode, KeyEvent, KeyEventKind},
    execute,
    terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode},
};
use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span, Text},
    widgets::{
        Block, BorderType, Borders, Cell, Clear, List, ListItem, ListState, Paragraph, Row, Table,
        TableState,
    },
};

/// The groups the tool-selection screen shows.
#[derive(Clone, Copy, PartialEq)]
enum ToolCategory {
    Alignment,
    Assembly,
    Index,
}

impl ToolCategory {
    fn name(self) -> &'static str {
        match self {
            ToolCategory::Alignment => "Alignment",
            ToolCategory::Assembly => "Assembly",
            ToolCategory::Index => "Indexing",
        }
    }
    fn color(self) -> Color {
        match self {
            ToolCategory::Alignment => Color::Cyan,
            ToolCategory::Assembly => Color::Green,
            ToolCategory::Index => Color::Magenta,
        }
    }
}

/// One tool of the selection screen and the parameters its form asks for.
#[derive(Clone, Copy)]
struct ToolInfo {
    name: &'static str,
    /// The optik subcommand the tool runs.
    cmd: &'static str,
    icon: &'static str,
    desc: &'static str,
    category: ToolCategory,
    params: &'static [ParamDef],
}

/// One form field, with its label, the command-line flag it fills, its kind, a default value and help text.
/// An empty value is not passed to the command, and a flag ending in '=' is joined to its value as one argument (--cs=long).
#[derive(Clone, Copy)]
struct ParamDef {
    name: &'static str,
    flag: &'static str,
    kind: ParamKind,
    default: &'static str,
    desc: &'static str,
    /// The command fails without this field, so the form marks it with `*` and refuses to run until it is filled.
    required: bool,
}

/// What a field holds, i.e. a path to read through the file browser, a path to write, a number, free text, or an on/off switch that passes its flag alone.
#[derive(Clone, Copy, PartialEq)]
enum ParamKind {
    FileRead,
    FileWrite,
    Int,
    Text,
    Flag,
}

/// The value a Flag field holds when it is on.
/// Off is the empty string.
const FLAG_ON: &str = "on";

const THREADS: ParamDef = ParamDef {
    name: "Threads",
    flag: "-t",
    kind: ParamKind::Int,
    default: "",
    desc: "Worker threads; empty = half of the logical CPUs",
    required: false,
};

/// The k-mer counting limits shared by dbruijn, build-fm and build-spss --compact.
const MAX_MEM: ParamDef = ParamDef {
    name: "Memory limit",
    flag: "--max-mem",
    kind: ParamKind::Text,
    default: "",
    desc: "Counting memory, e.g. 8G or 500M; past it buffers spill to disk; empty = half the available memory",
    required: false,
};
const TMP_DIR: ParamDef = ParamDef {
    name: "Spill directory",
    flag: "--tmp-dir",
    kind: ParamKind::Text,
    default: "",
    desc: "Directory for spilled buffers; empty = the output file's directory",
    required: false,
};
const GRAPH_MODE: ParamDef = ParamDef {
    name: "Graph mode",
    flag: "--graph",
    kind: ParamKind::Text,
    default: "",
    desc: "auto, part or table; identical output, part uses about 1/4 of the memory; empty = auto",
    required: false,
};

/// SAM flag bits the results screen reads.
const SAM_FLAG_UNMAPPED: i32 = 0x4;
const SAM_FLAG_REVERSE: i32 = 0x10;

/// Every tool the TUI offers, in display order and grouped by category.
const TOOLS: &[ToolInfo] = &[
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

/// The screen shown, with its state.
enum Screen {
    ToolSelect,
    FileBrowser {
        tool_idx: usize,
        param_idx: usize,
        values: Vec<String>,
        focus: usize,
        cwd: PathBuf,
        entries: Vec<(String, bool)>,
        cursor: usize,
    },
    ToolForm {
        tool_idx: usize,
        values: Vec<String>,
        focus: usize,
    },
    Running {
        tool_idx: usize,
        output: Vec<String>,
        stderr: Vec<String>,
        start: Instant,
        rx_stdout: Option<Receiver<String>>,
        rx_stderr: Option<Receiver<String>>,
        child: Option<Child>,
    },
    Results {
        tool_idx: usize,
        output: Vec<String>,
        exit_code: i32,
        scroll: usize,
        formatted: Vec<Line<'static>>,
        /// The width `formatted` was wrapped to.
        formatted_w: usize,
    },
}

/// TUI state, i.e. the current screen, the selection on the tool screen, the terminal size of the last frame and the help, notice and browser settings every screen shares.
pub(crate) struct App {
    screen: Screen,
    selected_cat: usize,
    selected_in_cat: usize,
    term_w: u16,
    term_h: u16,
    /// The `?` overlay is open; any key closes it.
    show_help: bool,
    /// First row of the help overlay shown, moved by the arrow and page keys.
    help_scroll: u16,
    /// A message for the form, such as the required fields still empty.
    notice: Option<String>,
    /// The equivalent command line of the run in progress or just finished.
    cmd_line: String,
    /// File-name suffixes the browser lists, empty for every file.
    browse_exts: &'static [&'static str],
    /// The browser lists every file regardless of browse_exts (key `a`).
    browse_all: bool,
}

/// Terminal size assumed before the first frame is drawn.
const DEFAULT_TERM_W: u16 = 80;
const DEFAULT_TERM_H: u16 = 24;
/// How often the running screen polls the child's output and redraws.
const RUN_POLL_INTERVAL_MS: u64 = 50;

impl App {
    fn new() -> Self {
        App {
            screen: Screen::ToolSelect,
            selected_cat: 0,
            selected_in_cat: 0,
            term_w: DEFAULT_TERM_W,
            term_h: DEFAULT_TERM_H,
            show_help: false,
            help_scroll: 0,
            notice: None,
            cmd_line: String::new(),
            browse_exts: &[],
            browse_all: false,
        }
    }
}

/// Sequence files a FASTA/FASTQ field accepts, plain or gzipped.
const SEQ_EXTS: &[&str] = &[
    ".fa",
    ".fasta",
    ".fna",
    ".fq",
    ".fastq",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fq.gz",
    ".fastq.gz",
];
/// Index files the query and update forms read.
const INDEX_EXTS: &[&str] = &[".spss", ".fm"];

/// The suffixes the browser should list for field @param of @tool, empty for any file (output paths, and fields that take no file).
fn wanted_exts(tool: &ToolInfo, param: &ParamDef) -> &'static [&'static str] {
    match param.kind {
        ParamKind::FileRead
            if param.flag == "-i" && matches!(tool.cmd, "query" | "spss-update") =>
        {
            INDEX_EXTS
        }
        ParamKind::FileRead => SEQ_EXTS,
        _ => &[],
    }
}

/// shell_quote() - @s as one shell word, quoted only when it needs to be.
fn shell_quote(s: &str) -> String {
    if !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || "-_./=:@+,%".contains(c))
    {
        s.to_string()
    } else {
        format!("'{}'", s.replace('\'', "'\\''"))
    }
}

/// The command line a form runs, `optik <cmd> <flags>`, quoted for a shell.
fn command_line(tool: &ToolInfo, values: &[String]) -> String {
    std::iter::once("optik".to_string())
        .chain(std::iter::once(tool.cmd.to_string()))
        .chain(build_args(tool, values).iter().map(|a| shell_quote(a)))
        .collect::<Vec<_>>()
        .join(" ")
}

/// The required fields of a form that are still empty, by name.
fn missing_required(tool: &ToolInfo, values: &[String]) -> Vec<&'static str> {
    tool.params
        .iter()
        .zip(values)
        .filter(|(p, v)| p.required && v.trim().is_empty())
        .map(|(p, _)| p.name)
        .collect()
}

/// Runs the TUI until the user quits, restoring the terminal afterwards, even after a panic.
pub(crate) fn run() -> Result<(), String> {
    let prev_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        use std::io::Write;
        let _ = disable_raw_mode();
        let mut stdout = io::stdout();
        let _ = execute!(stdout, LeaveAlternateScreen);
        let _ = stdout.flush();
        prev_hook(info);
    }));

    enable_raw_mode().map_err(|e| e.to_string())?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen).map_err(|e| e.to_string())?;
    let mut terminal = ratatui::Terminal::new(ratatui::backend::CrosstermBackend::new(stdout))
        .map_err(|e| e.to_string())?;
    terminal.clear().map_err(|e| e.to_string())?;
    let mut app = App::new();
    let res = run_app(&mut terminal, &mut app);
    disable_raw_mode().map_err(|e| e.to_string())?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen).map_err(|e| e.to_string())?;
    terminal.show_cursor().map_err(|e| e.to_string())?;
    res
}

/// The draw / handle-event loop.
fn run_app<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    app: &mut App,
) -> Result<(), String> {
    loop {
        let size = terminal.size().map_err(|e| e.to_string())?;
        app.term_w = size.width;
        app.term_h = size.height;
        terminal
            .draw(|f| render(f, &mut *app))
            .map_err(|e| e.to_string())?;
        if !handle_event(app)? {
            return Ok(());
        }
    }
}

/// The tools grouped by category, as (category, indices into TOOLS).
fn grouped_tools() -> Vec<(ToolCategory, Vec<usize>)> {
    let mut groups: Vec<(ToolCategory, Vec<usize>)> = Vec::new();
    for (i, tool) in TOOLS.iter().enumerate() {
        match groups.last_mut() {
            Some((cat, indices)) if *cat == tool.category => indices.push(i),
            _ => groups.push((tool.category, vec![i])),
        }
    }
    groups
}

/// Number of categories.
fn cat_count() -> usize {
    grouped_tools().len()
}

/// Number of tools in category `cat_idx`.
fn tools_in_cat(cat_idx: usize) -> usize {
    grouped_tools()
        .get(cat_idx)
        .map(|(_, v)| v.len())
        .unwrap_or(0)
}

/// Index into TOOLS of tool `in_cat` of category `cat_idx` (0 if out of range).
fn global_tool_idx(cat_idx: usize, in_cat: usize) -> usize {
    grouped_tools()
        .get(cat_idx)
        .and_then(|(_, v)| v.get(in_cat).copied())
        .unwrap_or(0)
}

/// The file browser's entries for `cwd`, i.e. "..", then directories, then files, each sorted with hidden names skipped, as (name, is_dir) pairs.
/// With @exts non-empty only files ending in one of them are listed.
fn refresh_entries(cwd: &Path, exts: &[&str]) -> Vec<(String, bool)> {
    let mut entries: Vec<(String, bool)> = Vec::new();
    if cwd.parent().is_some() {
        entries.push(("..".to_string(), true));
    }
    if let Ok(rd) = std::fs::read_dir(cwd) {
        let mut dirs: Vec<String> = Vec::new();
        let mut files: Vec<String> = Vec::new();
        for entry in rd.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if name.starts_with('.') {
                continue;
            }
            if let Ok(ft) = entry.file_type() {
                if ft.is_dir() {
                    dirs.push(name);
                } else if exts.is_empty()
                    || exts.iter().any(|e| name.to_ascii_lowercase().ends_with(e))
                {
                    files.push(name);
                }
            }
        }
        dirs.sort();
        files.sort();
        for d in dirs {
            entries.push((d, true));
        }
        for f in files {
            entries.push((f, false));
        }
    }
    entries
}

/// Lines of a tool's output the TUI keeps.
/// Everything is still read (the child must not block on a full pipe), but a map run without an output file streams the whole SAM here, and holding it all would take as much memory as the SAM itself.
const MAX_KEPT_LINES: usize = 200_000;

/// Keeps `line` if fewer than MAX_KEPT_LINES are kept, and adds one truncation notice when the cap is reached.
fn push_capped(v: &mut Vec<String>, line: String) {
    if v.len() < MAX_KEPT_LINES {
        v.push(line);
    } else if v.len() == MAX_KEPT_LINES {
        v.push(format!(
            "... output truncated: the TUI keeps the first {MAX_KEPT_LINES} lines \
                        (set an output file to keep everything)"
        ));
    }
}

/// Collects the running child's output and, once it has exited, moves to the results screen.
fn process_running(app: &mut App) {
    let prev = std::mem::replace(&mut app.screen, Screen::ToolSelect);
    app.screen = match prev {
        Screen::Running {
            tool_idx,
            mut output,
            mut stderr,
            start,
            rx_stdout,
            rx_stderr,
            mut child,
        } => {
            if let Some(rx) = &rx_stdout {
                while let Ok(line) = rx.try_recv() {
                    push_capped(&mut output, line);
                }
            }
            if let Some(rx) = &rx_stderr {
                while let Ok(line) = rx.try_recv() {
                    push_capped(&mut stderr, line);
                }
            }
            if let Some(c) = &mut child {
                match c.try_wait() {
                    Ok(Some(status)) => {
                        if let Some(rx) = &rx_stdout {
                            while let Ok(line) = rx.try_recv() {
                                push_capped(&mut output, line);
                            }
                        }
                        if let Some(rx) = &rx_stderr {
                            while let Ok(line) = rx.try_recv() {
                                push_capped(&mut stderr, line);
                            }
                        }
                        let exit_code = status.code().unwrap_or(-1);
                        let mut combined = output;
                        combined.extend(stderr);
                        Screen::Results {
                            tool_idx,
                            output: combined,
                            exit_code,
                            scroll: 0,
                            formatted: Vec::new(),
                            formatted_w: 0,
                        }
                    }
                    Ok(None) => Screen::Running {
                        tool_idx,
                        output,
                        stderr,
                        start,
                        rx_stdout,
                        rx_stderr,
                        child,
                    },
                    Err(e) => {
                        output.push(format!("Process error: {e}"));
                        let mut combined = output;
                        combined.extend(stderr);
                        Screen::Results {
                            tool_idx,
                            output: combined,
                            exit_code: -1,
                            scroll: 0,
                            formatted: Vec::new(),
                            formatted_w: 0,
                        }
                    }
                }
            } else {
                Screen::Running {
                    tool_idx,
                    output,
                    stderr,
                    start,
                    rx_stdout,
                    rx_stderr,
                    child,
                }
            }
        }
        other => other,
    };
}

/// Handles one key event, or a poll timeout while a tool runs, and returns false to quit.
fn handle_event(app: &mut App) -> Result<bool, String> {
    if matches!(app.screen, Screen::Running { .. }) {
        process_running(app);
        if matches!(app.screen, Screen::Running { .. }) {
            let timeout = std::time::Duration::from_millis(RUN_POLL_INTERVAL_MS);
            if event::poll(timeout).map_err(|e| e.to_string())?
                && let Event::Key(key) = event::read().map_err(|e| e.to_string())?
                && key.kind == KeyEventKind::Press
                && matches!(key.code, KeyCode::Char('q') | KeyCode::Esc)
            {
                // Dropping the handle would leave the child to die of SIGPIPE at its next write, so it is stopped explicitly.
                if let Screen::Running { child: Some(c), .. } = &mut app.screen {
                    let _ = c.kill();
                    let _ = c.wait();
                }
                app.screen = Screen::ToolSelect;
                app.selected_cat = 0;
                app.selected_in_cat = 0;
            }
        }
        return Ok(true);
    }

    if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
        if key.kind != KeyEventKind::Press {
            return Ok(true);
        }
        // The arrow and page keys scroll the help overlay, and any other key closes it.
        // F1 opens it anywhere, and `?` everywhere except while a form field takes typed text.
        if app.show_help {
            let page = app.term_h.saturating_sub(4).max(1);
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    app.help_scroll = app.help_scroll.saturating_sub(1)
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    app.help_scroll = app.help_scroll.saturating_add(1)
                }
                KeyCode::PageUp => app.help_scroll = app.help_scroll.saturating_sub(page),
                KeyCode::PageDown => app.help_scroll = app.help_scroll.saturating_add(page),
                KeyCode::Home => app.help_scroll = 0,
                KeyCode::End => app.help_scroll = u16::MAX,
                _ => {
                    app.show_help = false;
                    app.help_scroll = 0;
                }
            }
            return Ok(true);
        }
        let typing = matches!(&app.screen, Screen::ToolForm { tool_idx, focus, .. }
            if *focus < TOOLS[*tool_idx].params.len()
                && TOOLS[*tool_idx].params[*focus].kind != ParamKind::Flag);
        if key.code == KeyCode::F(1) || (key.code == KeyCode::Char('?') && !typing) {
            app.show_help = true;
            return Ok(true);
        }
        app.notice = None;
        let app_browse_all = app.browse_all;
        let app_browse_exts = app.browse_exts;
        match &mut app.screen {
            Screen::ToolSelect => return Ok(handle_tool_select(app, key)),
            Screen::FileBrowser {
                tool_idx,
                param_idx,
                values,
                focus,
                cwd,
                entries,
                cursor,
            } => match key.code {
                KeyCode::Esc => {
                    app.screen = Screen::ToolForm {
                        tool_idx: *tool_idx,
                        values: values.clone(),
                        focus: *focus,
                    };
                }
                KeyCode::Up | KeyCode::Char('k') => {
                    *cursor = cursor.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    *cursor = (*cursor + 1).min(entries.len().saturating_sub(1));
                }
                KeyCode::Enter => {
                    if let Some((ref name, is_dir)) = entries.get(*cursor).cloned() {
                        if is_dir {
                            if name == ".." {
                                if let Some(parent) = cwd.parent() {
                                    *cwd = parent.to_path_buf();
                                    *entries = refresh_entries(
                                        cwd,
                                        if app_browse_all { &[] } else { app_browse_exts },
                                    );
                                    *cursor = 0;
                                }
                            } else {
                                let new_path = cwd.join(name);
                                if new_path.is_dir() {
                                    *cwd = new_path;
                                    *entries = refresh_entries(
                                        cwd,
                                        if app_browse_all { &[] } else { app_browse_exts },
                                    );
                                    *cursor = 0;
                                }
                            }
                        } else {
                            let full_path = cwd.join(name);
                            values[*param_idx] = full_path.to_string_lossy().to_string();
                            app.screen = Screen::ToolForm {
                                tool_idx: *tool_idx,
                                values: values.clone(),
                                focus: *focus,
                            };
                        }
                    }
                }
                KeyCode::Char('q') => {
                    app.screen = Screen::ToolForm {
                        tool_idx: *tool_idx,
                        values: values.clone(),
                        focus: *focus,
                    };
                }
                KeyCode::Char('a') => {
                    app.browse_all = !app_browse_all;
                    *entries =
                        refresh_entries(cwd, if app.browse_all { &[] } else { app_browse_exts });
                    *cursor = 0;
                }
                _ => {}
            },
            Screen::Results { scroll, .. } => match key.code {
                KeyCode::Char('q') | KeyCode::Esc => {
                    app.screen = Screen::ToolSelect;
                    app.selected_cat = 0;
                    app.selected_in_cat = 0;
                }
                KeyCode::Up | KeyCode::Char('k') => {
                    *scroll = scroll.saturating_sub(3);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    *scroll = scroll.saturating_add(3);
                }
                KeyCode::PageUp => {
                    *scroll = scroll.saturating_sub(app.term_h as usize);
                }
                KeyCode::PageDown => {
                    *scroll = scroll.saturating_add(app.term_h as usize);
                }
                KeyCode::Home => {
                    *scroll = 0;
                }
                KeyCode::End => {
                    *scroll = usize::MAX;
                }
                _ => {}
            },
            Screen::ToolForm {
                tool_idx,
                values,
                focus,
            } => {
                let tool = &TOOLS[*tool_idx];
                let nparams = tool.params.len();
                match key.code {
                    KeyCode::Esc => {
                        app.screen = Screen::ToolSelect;
                        app.selected_cat = 0;
                        app.selected_in_cat = 0;
                    }
                    // no j/k here, since the fields take typed paths and values
                    KeyCode::Tab | KeyCode::Down => {
                        *focus = (*focus + 1) % (nparams + 1);
                    }
                    KeyCode::BackTab | KeyCode::Up => {
                        *focus = focus.saturating_sub(1);
                    }
                    KeyCode::Enter if *focus == nparams => {
                        // Run only with every required field filled, else name the missing ones and move to the first.
                        let missing = missing_required(tool, values);
                        if let Some(first) = tool
                            .params
                            .iter()
                            .zip(values.iter())
                            .position(|(p, v)| p.required && v.trim().is_empty())
                        {
                            app.notice = Some(format!("Fill in first: {}", missing.join(", ")));
                            *focus = first;
                        } else {
                            app.cmd_line = command_line(tool, values);
                            let args = build_args(tool, values);
                            let full_args: Vec<String> =
                                std::iter::once(tool.cmd.to_string()).chain(args).collect();
                            let exe = std::env::current_exe().unwrap_or_else(|_| "optik".into());
                            let exe_s = exe.to_string_lossy().to_string();
                            match std::process::Command::new(&exe_s)
                                .args(&full_args)
                                .stdout(std::process::Stdio::piped())
                                .stderr(std::process::Stdio::piped())
                                .spawn()
                            {
                                Ok(mut child) => {
                                    let (tx_stdout, rx_stdout) = std::sync::mpsc::channel();
                                    let (tx_stderr, rx_stderr) = std::sync::mpsc::channel();
                                    if let Some(out) = child.stdout.take() {
                                        std::thread::spawn(move || {
                                            for line in io::BufReader::new(out).lines() {
                                                if let Ok(l) = line
                                                    && tx_stdout.send(l).is_err()
                                                {
                                                    break;
                                                }
                                            }
                                        });
                                    }
                                    if let Some(err) = child.stderr.take() {
                                        std::thread::spawn(move || {
                                            for line in io::BufReader::new(err).lines() {
                                                if let Ok(l) = line
                                                    && tx_stderr.send(l).is_err()
                                                {
                                                    break;
                                                }
                                            }
                                        });
                                    }
                                    app.screen = Screen::Running {
                                        tool_idx: *tool_idx,
                                        output: Vec::new(),
                                        stderr: Vec::new(),
                                        start: Instant::now(),
                                        rx_stdout: Some(rx_stdout),
                                        rx_stderr: Some(rx_stderr),
                                        child: Some(child),
                                    };
                                }
                                Err(e) => {
                                    app.screen = Screen::Results {
                                        tool_idx: *tool_idx,
                                        output: vec![format!("Error: failed to start tool: {e}")],
                                        exit_code: -1,
                                        scroll: 0,
                                        formatted: Vec::new(),
                                        formatted_w: 0,
                                    };
                                }
                            }
                        }
                    }
                    KeyCode::Enter
                        if *focus < nparams && tool.params[*focus].kind != ParamKind::Flag =>
                    {
                        if matches!(
                            tool.params[*focus].kind,
                            ParamKind::FileRead | ParamKind::FileWrite
                        ) {
                            let cwd = if values[*focus].is_empty() {
                                std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
                            } else {
                                let p = Path::new(&values[*focus]);
                                if p.is_dir() {
                                    p.to_path_buf()
                                } else {
                                    p.parent().map(|p| p.to_path_buf()).unwrap_or_else(|| {
                                        std::env::current_dir()
                                            .unwrap_or_else(|_| PathBuf::from("."))
                                    })
                                }
                            };
                            // The browser lists the files this field expects, `a` shows all.
                            let exts = wanted_exts(tool, &tool.params[*focus]);
                            app.browse_exts = exts;
                            app.browse_all = false;
                            let entries = refresh_entries(&cwd, exts);
                            app.screen = Screen::FileBrowser {
                                tool_idx: *tool_idx,
                                param_idx: *focus,
                                values: values.clone(),
                                focus: *focus,
                                cwd,
                                entries,
                                cursor: 0,
                            };
                        }
                    }
                    KeyCode::Enter | KeyCode::Char(' ')
                        if *focus < nparams && tool.params[*focus].kind == ParamKind::Flag =>
                    {
                        let v = &mut values[*focus];
                        *v = if v.is_empty() {
                            FLAG_ON.to_string()
                        } else {
                            String::new()
                        };
                    }
                    KeyCode::Char(c)
                        if *focus < nparams && tool.params[*focus].kind != ParamKind::Flag =>
                    {
                        values[*focus].push(c);
                    }
                    KeyCode::Backspace if *focus < nparams => {
                        values[*focus].pop();
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }
    Ok(true)
}

/// Opens the parameter form of tool `tool_idx` with its default values.
fn select_tool(app: &mut App, tool_idx: usize) {
    if tool_idx < TOOLS.len() {
        let tool = &TOOLS[tool_idx];
        let values: Vec<String> = tool.params.iter().map(|p| p.default.to_string()).collect();
        let groups = grouped_tools();
        for (gi, (_, indices)) in groups.iter().enumerate() {
            if let Some(ci) = indices.iter().position(|&ti| ti == tool_idx) {
                app.selected_cat = gi;
                app.selected_in_cat = ci;
                break;
            }
        }
        app.screen = Screen::ToolForm {
            tool_idx,
            values,
            focus: 0,
        };
    }
}

/// The tool whose one-letter icon is `c` (case-insensitive), the keyboard shortcut on the selection screen.
fn tool_by_icon_char(c: char) -> Option<usize> {
    let c_upper = c.to_ascii_uppercase();
    TOOLS.iter().position(|t| {
        let icon = t.icon.trim();
        icon.len() == 1
            && icon
                .chars()
                .next()
                .map(|ic| ic.to_ascii_uppercase() == c_upper)
                .unwrap_or(false)
    })
}

/// Returns false to quit, and the caller restores the terminal.
fn handle_tool_select(app: &mut App, key: KeyEvent) -> bool {
    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => return false,
        KeyCode::Up | KeyCode::Char('k') => {
            let nc = cat_count();
            if nc > 0 {
                app.selected_cat = app.selected_cat.saturating_sub(1);
                let nic = tools_in_cat(app.selected_cat);
                if app.selected_in_cat >= nic {
                    app.selected_in_cat = nic.saturating_sub(1);
                }
            }
        }
        KeyCode::Down | KeyCode::Char('j') => {
            let nc = cat_count();
            if nc > 0 {
                app.selected_cat = (app.selected_cat + 1).min(nc - 1);
                let nic = tools_in_cat(app.selected_cat);
                if app.selected_in_cat >= nic {
                    app.selected_in_cat = nic.saturating_sub(1);
                }
            }
        }
        KeyCode::Left | KeyCode::Char('h') => {
            let nic = tools_in_cat(app.selected_cat);
            if nic > 0 {
                app.selected_in_cat = app.selected_in_cat.saturating_sub(1);
            }
        }
        KeyCode::Right | KeyCode::Char('l') => {
            let nic = tools_in_cat(app.selected_cat);
            if nic > 0 {
                app.selected_in_cat = (app.selected_in_cat + 1).min(nic - 1);
            }
        }
        KeyCode::Enter => {
            let tool_idx = global_tool_idx(app.selected_cat, app.selected_in_cat);
            select_tool(app, tool_idx);
        }
        KeyCode::Char(c) => {
            if let Some(tool_idx) = tool_by_icon_char(c) {
                select_tool(app, tool_idx);
            }
        }
        _ => {}
    }
    true
}

/// The command-line arguments of a filled-in form, one per non-empty field: a flag and its value, a joined `flag=value`, or a switched-on flag alone.
fn build_args(tool: &ToolInfo, values: &[String]) -> Vec<String> {
    let mut args = Vec::new();
    for (i, param) in tool.params.iter().enumerate() {
        let v = values[i].trim();
        if v.is_empty() {
            continue;
        }
        if param.kind == ParamKind::Flag {
            args.push(param.flag.to_string());
        } else if param.flag.ends_with('=') {
            args.push(format!("{}{v}", param.flag));
        } else {
            args.push(param.flag.to_string());
            args.push(v.to_string());
        }
    }
    args
}

/// Draws the current screen.
fn render(f: &mut Frame, app: &mut App) {
    match &mut app.screen {
        Screen::ToolSelect => render_tool_select(f, app),
        Screen::FileBrowser {
            tool_idx,
            param_idx,
            cwd,
            entries,
            cursor,
            ..
        } => render_file_browser(f, *tool_idx, *param_idx, cwd, entries, *cursor),
        Screen::ToolForm {
            tool_idx,
            values,
            focus,
        } => render_tool_form(f, *tool_idx, values, *focus, app.notice.as_deref()),
        Screen::Running {
            tool_idx,
            output,
            stderr,
            start,
            ..
        } => render_running(f, *tool_idx, output, stderr, *start, &app.cmd_line),
        Screen::Results {
            tool_idx,
            output,
            exit_code,
            scroll,
            formatted,
            formatted_w,
        } => render_results(
            f,
            *tool_idx,
            output,
            *exit_code,
            scroll,
            formatted,
            formatted_w,
            &app.cmd_line,
        ),
    }
    if app.show_help {
        let mut scroll = app.help_scroll;
        render_help(f, app, &mut scroll);
        app.help_scroll = scroll;
    }
}

/// The title bar shown on every screen.
fn header() -> Block<'static> {
    Block::default()
        .title(Line::from(vec![
            Span::styled(
                " Optik ",
                Style::default()
                    .fg(Color::Black)
                    .bg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" "),
            Span::styled(
                "Genomic Sequence Analysis Toolkit",
                Style::default().fg(Color::Cyan),
            ),
            Span::raw(concat!("  v", env!("CARGO_PKG_VERSION"))),
        ]))
        .title_alignment(ratatui::layout::Alignment::Left)
        .borders(Borders::BOTTOM)
        .border_style(Style::default().fg(Color::Cyan))
}

/// header_rows() - Rows the header takes at @width, three, or more when the command line @cmd needs them.
fn header_rows(cmd: &str, width: u16) -> u16 {
    2 + cmd_rows(cmd, width).len().max(1) as u16
}

/// The command line of a run as wrapped `$ optik ...` rows for the header, none when @cmd is empty.
fn cmd_rows(cmd: &str, width: u16) -> Vec<Line<'static>> {
    if cmd.is_empty() {
        return Vec::new();
    }
    let line = Line::from(vec![
        Span::styled("$ ", Style::default().fg(Color::DarkGray)),
        Span::styled(cmd.to_string(), Style::default().fg(Color::Cyan)),
    ]);
    wrap_line(line, width.saturating_sub(2) as usize, 2)
}

/// Draws the title bar in @area, with the run's command line @cmd under the title when it is not empty.
fn render_header(f: &mut Frame, area: Rect, cmd: &str) {
    f.render_widget(header(), area);
    let rows = cmd_rows(cmd, area.width);
    if !rows.is_empty() {
        let body = Rect {
            x: area.x + 1,
            y: area.y + 1,
            width: area.width.saturating_sub(2),
            height: area.height.saturating_sub(2),
        };
        f.render_widget(Paragraph::new(rows), body);
    }
}

/// footer_rows() - Rows the key-help footer @text takes at @width.
fn footer_rows(text: &str, width: u16) -> u16 {
    if text.chars().count() + 2 <= width as usize {
        1
    } else {
        wrap_text(text, width.saturating_sub(2) as usize, 0, Style::default()).len() as u16
    }
}

/// Draws the key-help footer @text in @area, as one rule-line title when it fits and wrapped centered rows when it does not.
fn render_footer(f: &mut Frame, area: Rect, text: &str) {
    if area.height <= 1 {
        f.render_widget(footer(text), area);
        return;
    }
    let style = Style::default()
        .fg(Color::Black)
        .bg(Color::DarkGray)
        .add_modifier(Modifier::BOLD);
    let rows = wrap_text(text, area.width.saturating_sub(2) as usize, 0, style);
    f.render_widget(
        Paragraph::new(rows).alignment(ratatui::layout::Alignment::Center),
        area,
    );
}

/// The vertical split of a screen with @cmd in the header and @foot as the footer: header, body, footer.
fn screen_layout(area: Rect, cmd: &str, foot: &str) -> std::rc::Rc<[Rect]> {
    Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(header_rows(cmd, area.width)),
            Constraint::Min(1),
            Constraint::Length(footer_rows(foot, area.width)),
        ])
        .split(area)
}

/// The key-help line at the bottom of a screen.
fn footer(text: &str) -> Block<'static> {
    Block::default()
        .title(Line::from(Span::styled(
            format!(" {} ", text),
            Style::default()
                .fg(Color::Black)
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        )))
        .title_alignment(ratatui::layout::Alignment::Center)
        .borders(Borders::TOP)
        .border_style(Style::default().fg(Color::DarkGray))
}

/// Narrowest and widest tool card on the selection screen, in columns.
const MIN_CARD_W: u16 = 22;
const MAX_CARD_W: u16 = 42;

/// Draws the tool-selection screen.
fn render_tool_select(f: &mut Frame, app: &App) {
    let area = f.area();
    let foot = "^/v Group | <-/-> Tool | Enter or M/D/B/F/Q/U Select | ? Help | q Quit";
    let layout = screen_layout(area, "", foot);

    render_header(f, layout[0], "");
    f.render_widget(Clear, layout[1]);
    let inner = layout[1];
    let groups = grouped_tools();
    let gap = 2u16;
    let cat_gap = 1u16;

    // Each group's card geometry: (cards per row, card width, card height, rows the group takes).
    let geometry: Vec<(u16, u16, u16, u16)> = groups
        .iter()
        .map(|(_, indices)| {
            // Cards wrap onto further rows when a group has more than fit side by side.
            let n_tools = indices.len() as u16;
            let per_row = ((inner.width + gap) / (MIN_CARD_W + gap)).clamp(1, n_tools.max(1));
            let total_gap = gap * per_row.saturating_sub(1);
            let avail_w = inner.width.saturating_sub(total_gap);
            let card_w = (avail_w / per_row).clamp(MIN_CARD_W, MAX_CARD_W);
            // Cards are as tall as the longest wrapped description of the group needs.
            let inner_w = card_w.saturating_sub(2) as usize;
            let desc_rows = indices
                .iter()
                .map(|&ti| wrap_text(TOOLS[ti].desc, inner_w, 0, Style::default()).len())
                .max()
                .unwrap_or(1);
            let card_h = (4 + desc_rows as u16).max(6);
            let rows = 1 + n_tools.div_ceil(per_row) * card_h + cat_gap;
            (per_row, card_w, card_h, rows)
        })
        .collect();
    // When the groups are taller than the screen, the view starts at the first group that keeps the selected one fully on screen, so every tool can be reached and seen.
    let sel = app.selected_cat.min(groups.len().saturating_sub(1));
    let mut first = 0;
    while first < sel && geometry[first..=sel].iter().map(|g| g.3).sum::<u16>() > inner.height {
        first += 1;
    }
    let mut y = inner.y;
    if first > 0 {
        f.render_widget(
            Paragraph::new(Span::styled(
                format!("  ^ {first} more group(s) above"),
                Style::default().fg(Color::DarkGray),
            )),
            Rect { height: 1, ..inner },
        );
        y += 1;
    }

    for (gi, (cat, indices)) in groups.iter().enumerate().skip(first) {
        if y >= inner.y + inner.height {
            break;
        }
        let (per_row, card_w, card_h, rows) = geometry[gi];

        let is_sel_cat = gi == app.selected_cat;
        let hdr_style = if is_sel_cat {
            Style::default()
                .fg(cat.color())
                .add_modifier(Modifier::BOLD | Modifier::REVERSED)
        } else {
            Style::default()
                .fg(cat.color())
                .add_modifier(Modifier::BOLD)
        };
        let cat_hdr = Paragraph::new(Line::from(vec![
            Span::raw("  "),
            Span::styled(format!(" {} ", cat.name()), hdr_style),
        ]));
        let ha = Rect {
            x: inner.x,
            y,
            width: inner.width,
            height: 1,
        };
        f.render_widget(cat_hdr, ha);

        for (ci, &tool_idx) in indices.iter().enumerate() {
            let (row, col) = (ci as u16 / per_row, ci as u16 % per_row);
            let card_area = Rect {
                x: inner.x + col * (card_w + gap),
                y: y + 1 + row * card_h,
                width: card_w,
                height: card_h,
            };
            if card_area.right() > inner.x + inner.width
                || card_area.bottom() > inner.y + inner.height
            {
                continue;
            }
            let selected = is_sel_cat && ci == app.selected_in_cat;
            render_tool_card(f, card_area, tool_idx, selected);
        }

        y += rows;
    }

    render_footer(f, layout[2], foot);
}

/// Draws one tool's card on the selection screen.
fn render_tool_card(f: &mut Frame, area: Rect, tool_idx: usize, selected: bool) {
    let tool = &TOOLS[tool_idx];
    let bc = if selected {
        tool.category.color()
    } else {
        Color::DarkGray
    };
    let border_style = if selected {
        Style::default().fg(bc).add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(bc)
    };
    let block = Block::default()
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(border_style);
    let inner = block.inner(area);
    f.render_widget(block, area);

    let inner_w = inner.width as usize;
    let name = format!("{} {}", tool.icon, tool.name);
    let mut lines = wrap_text(
        &name,
        inner_w,
        0,
        Style::default()
            .fg(tool.category.color())
            .add_modifier(Modifier::BOLD),
    );
    lines.push(Line::from(Span::raw("")));
    lines.extend(wrap_text(
        tool.desc,
        inner_w,
        0,
        Style::default().fg(Color::Gray),
    ));

    f.render_widget(Paragraph::new(lines), inner);
}

#[allow(clippy::too_many_arguments)]
/// Draws the file browser over the form.
fn render_file_browser(
    f: &mut Frame,
    tool_idx: usize,
    param_idx: usize,
    cwd: &Path,
    entries: &[(String, bool)],
    cursor: usize,
) {
    let tool = &TOOLS[tool_idx];
    let area = f.area();
    let foot = "^/v Move  |  Enter Select/Open dir  |  a All files  |  ? Help  |  Esc Back";
    let layout = screen_layout(area, "", foot);
    render_header(f, layout[0], "");
    f.render_widget(Clear, layout[1]);

    let param = &tool.params[param_idx];
    let kind_tag = match param.kind {
        ParamKind::FileRead => "Select File",
        ParamKind::FileWrite => "Select Output Path",
        _ => "Select File",
    };

    let title_fmt = format!(
        " {} {} -- {} {} ",
        tool.icon, tool.name, kind_tag, param.name
    );

    let fb_area = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(10),
            Constraint::Percentage(80),
            Constraint::Percentage(10),
        ])
        .split(layout[1])[1];

    let path_str = cwd.to_string_lossy().to_string();
    let block = Block::default()
        .title(title_fmt)
        .title_alignment(ratatui::layout::Alignment::Left)
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(Style::default().fg(Color::Cyan));

    let inner = block.inner(fb_area);
    f.render_widget(block, fb_area);

    let path_style = Style::default()
        .fg(Color::Cyan)
        .add_modifier(Modifier::BOLD);
    let path_lines = wrap_text(
        &format!(" {}", path_str),
        inner.width as usize,
        1,
        path_style,
    );
    let path_h = (path_lines.len() as u16)
        .min(inner.height.saturating_sub(1))
        .max(1);
    let path_area = Rect {
        x: inner.x,
        y: inner.y,
        width: inner.width,
        height: path_h,
    };
    f.render_widget(Paragraph::new(path_lines), path_area);

    let list_y = inner.y + path_h;
    let list_h = inner.height.saturating_sub(path_h);
    let list_area = Rect {
        x: inner.x,
        y: list_y,
        width: inner.width,
        height: list_h,
    };

    let items: Vec<ListItem> = entries
        .iter()
        .enumerate()
        .map(|(i, (name, is_dir))| {
            let prefix = if *is_dir { " > " } else { "   " };
            let suffix = if *is_dir { "/" } else { "" };
            let style = if i == cursor {
                Style::default().add_modifier(Modifier::REVERSED)
            } else if *is_dir {
                Style::default().fg(Color::Yellow)
            } else {
                Style::default().fg(Color::White)
            };
            // A long name continues on further rows under its first character.
            ListItem::new(wrap_text(
                &format!("{}{}{}", prefix, name, suffix),
                list_area.width as usize,
                prefix.len(),
                style,
            ))
        })
        .collect();

    if !items.is_empty() {
        let list =
            List::new(items).highlight_style(Style::default().add_modifier(Modifier::REVERSED));
        let mut state = ListState::default().with_selected(Some(cursor));
        f.render_stateful_widget(list, list_area, &mut state);
    }

    render_footer(f, layout[2], foot);
}

/// Draws a tool's parameter form.
fn render_tool_form(
    f: &mut Frame,
    tool_idx: usize,
    values: &[String],
    focus: usize,
    notice: Option<&str>,
) {
    let tool = &TOOLS[tool_idx];
    let area = f.area();
    let foot = "^/v Move  |  type a value  |  Enter Browse/Toggle/Run  |  F1 Help  |  Esc Back";
    let layout = screen_layout(area, "", foot);
    render_header(f, layout[0], "");
    f.render_widget(Clear, layout[1]);

    let content = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(10),
            Constraint::Percentage(80),
            Constraint::Percentage(10),
        ])
        .split(layout[1])[1];
    let width = content.width as usize;
    let nparams = tool.params.len();
    // Each block below starts where the one above ended, so wrapped text never overlaps.
    let mut y = content.y;

    place_rows(
        f,
        content,
        layout[1],
        wrap_text(
            &format!(" {} {}", tool.icon, tool.desc),
            width,
            4,
            Style::default().fg(Color::Gray),
        ),
        &mut y,
    );

    // Column widths inside the table border: name, value and a type tag that always fits.
    let inner_w = content.width.saturating_sub(2);
    let tag_w: u16 = 10;
    let name_w = inner_w * 2 / 5;
    let val_w = inner_w.saturating_sub(name_w + tag_w + 2);
    let mut rows = Vec::new();
    let mut rows_h: u16 = 0;
    for (i, p) in tool.params.iter().enumerate() {
        let is_focus = focus == i;
        let disp = if p.kind == ParamKind::Flag {
            if values[i].is_empty() { "off" } else { "on" }.to_string()
        } else if values[i].is_empty() {
            if is_focus {
                " |".into()
            } else if p.required {
                "(required)".into()
            } else if p.default.is_empty() && p.kind == ParamKind::FileRead {
                "(none)".into()
            } else if p.default.is_empty() {
                "(auto)".into()
            } else {
                format!("(default: {})", p.default)
            }
        } else {
            values[i].clone()
        };
        let val_style = if is_focus {
            Style::default()
                .fg(Color::Black)
                .bg(Color::Cyan)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::default()
        };
        let name_style = if is_focus {
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::default().fg(Color::White)
        };
        let tag_style = if is_focus {
            Style::default().fg(Color::Cyan)
        } else {
            Style::default().fg(Color::DarkGray)
        };
        let tag = match p.kind {
            ParamKind::FileRead => " [file]",
            ParamKind::FileWrite => " [out]",
            ParamKind::Int => " [int]",
            ParamKind::Text => " [text]",
            ParamKind::Flag => " [on/off]",
        };
        let name_lines = wrap_text(
            &format!(" {}{}", p.name, if p.required { " *" } else { "" }),
            name_w as usize,
            1,
            name_style,
        );
        let val_lines = wrap_text(&format!(" {}", disp), val_w as usize, 1, val_style);
        let h = name_lines.len().max(val_lines.len()) as u16;
        rows_h += h;
        rows.push(
            Row::new(vec![
                Cell::from(Text::from(name_lines)),
                Cell::from(Text::from(val_lines)),
                Cell::from(Line::from(Span::styled(tag, tag_style))),
            ])
            .height(h),
        );
    }

    let title_fmt = format!(" {} {} ", tool.icon, tool.name);
    let table = Table::new(
        rows,
        vec![
            Constraint::Length(name_w),
            Constraint::Length(val_w),
            Constraint::Length(tag_w),
        ],
    )
    .block(
        Block::default()
            .title(title_fmt)
            .title_alignment(ratatui::layout::Alignment::Left)
            .borders(Borders::ALL)
            .border_type(BorderType::Rounded)
            .border_style(Style::default().fg(tool.category.color())),
    );
    let sel = if nparams > 0 {
        Some(focus.min(nparams.saturating_sub(1)))
    } else {
        None
    };
    let mut state = TableState::default().with_selected(sel);

    // The field help, the command preview and the notice are wrapped first, so the table gets the rows they and the Run button leave, and scrolls to the focused field when taller.
    let info_lines = if focus < nparams {
        let p = &tool.params[focus];
        let text = if p.default.is_empty() {
            format!("  {}", p.desc)
        } else {
            format!("  {} (default: {})", p.desc, p.default)
        };
        wrap_text(&text, width, 2, Style::default().fg(Color::DarkGray))
    } else {
        Vec::new()
    };
    let cmd = Line::from(vec![
        Span::styled("  $ ", Style::default().fg(Color::DarkGray)),
        Span::styled(command_line(tool, values), Style::default().fg(Color::Cyan)),
    ]);
    let cmd_lines = wrap_line(cmd, width, 4);
    let notice_lines = notice
        .map(|msg| {
            wrap_text(
                &format!("  {msg}"),
                width,
                2,
                Style::default().fg(Color::Red).add_modifier(Modifier::BOLD),
            )
        })
        .unwrap_or_default();
    let below_need =
        (1 + info_lines.len() + cmd_lines.len() + notice_lines.len()) as u16 + RUN_BUTTON_H;
    let below = content.bottom().saturating_sub(y);
    let table_h = (rows_h + 2)
        .min(below.saturating_sub(below_need))
        .max(MIN_TABLE_H);
    let table_area = Rect {
        x: content.x,
        y,
        width: content.width,
        height: table_h,
    }
    .intersection(layout[1]);
    f.render_stateful_widget(&table, table_area, &mut state);
    y += table_h + 1;

    // The field help, then the command the form would run, then a notice such as missing required fields
    place_rows(f, content, layout[1], info_lines, &mut y);
    place_rows(f, content, layout[1], cmd_lines, &mut y);
    place_rows(f, content, layout[1], notice_lines, &mut y);
    let btn_area = Rect {
        x: content.x + (content.width / 2).saturating_sub(10),
        y,
        width: 20.min(content.width),
        height: RUN_BUTTON_H,
    }
    .intersection(layout[1]);
    let is_run = focus >= nparams;
    let (btn_fg, btn_bg) = if is_run {
        (Color::Black, Color::Green)
    } else {
        (Color::Green, Color::Reset)
    };
    let btn = Block::default()
        .title(" Run ")
        .title_alignment(ratatui::layout::Alignment::Center)
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(
            Style::default()
                .fg(Color::Green)
                .add_modifier(Modifier::BOLD),
        )
        .style(Style::default().bg(btn_bg).fg(btn_fg));
    f.render_widget(btn, btn_area);

    render_footer(f, layout[2], foot);
}

/// Rows of the form's Run button, border included
const RUN_BUTTON_H: u16 = 3;
/// Fewest rows the form's table keeps, its border and one field
const MIN_TABLE_H: u16 = 3;

/// place_rows() - Draw @lines in the columns of @col from row @y down, clipped to @bound, and move @y below them.
fn place_rows(f: &mut Frame, col: Rect, bound: Rect, lines: Vec<Line<'static>>, y: &mut u16) {
    let h = lines.len() as u16;
    let r = Rect {
        x: col.x,
        y: *y,
        width: col.width,
        height: h,
    }
    .intersection(bound);
    f.render_widget(Paragraph::new(lines), r);
    *y += h;
}

/// Draws the running screen with the elapsed time and the tail of the output.
fn render_running(
    f: &mut Frame,
    tool_idx: usize,
    output: &[String],
    stderr: &[String],
    start: Instant,
    cmd: &str,
) {
    let tool = &TOOLS[tool_idx];
    let area = f.area();
    let foot = "Running...  q/Esc stop the run and go back  |  ? Help";
    let layout = screen_layout(area, cmd, foot);
    render_header(f, layout[0], cmd);
    f.render_widget(Clear, layout[1]);

    let elapsed = start.elapsed().as_secs_f64();
    const SPINNER: &str = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏";
    let sp_i = (start.elapsed().as_millis() / 100 % 10) as usize;
    let sp_char = SPINNER.chars().nth(sp_i).unwrap_or(' ');

    let title = format!(
        " {} {}  |  {:.1}s  |  out: {}  |  err: {}  |  {} ",
        tool.icon,
        tool.name,
        elapsed,
        output.len(),
        stderr.len(),
        sp_char,
    );

    let block = Block::default()
        .title(Line::from(Span::styled(
            &title,
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        )))
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(Style::default().fg(Color::Yellow));

    let inner = block.inner(layout[1]);
    f.render_widget(block, layout[1]);

    let max_w = inner.width as usize;
    let has_err = !stderr.is_empty();

    let err_h = if has_err {
        (inner.height / 4).clamp(2, 6)
    } else {
        0
    };
    let out_h = inner.height.saturating_sub(err_h);

    // Only the tail is shown, since formatting all kept lines, up to MAX_KEPT_LINES of SAM, on every 50 ms redraw would stall the screen.
    let tail = output.len().saturating_sub(4 * out_h as usize + 64)..output.len();
    let formatted = format_tool_output(output, max_w, tail);
    let out_start = formatted.len().saturating_sub(out_h as usize);
    let out_lines: Vec<Line> = formatted.into_iter().skip(out_start).collect();

    let out_area = Rect {
        x: inner.x,
        y: inner.y,
        width: inner.width,
        height: out_h,
    };
    f.render_widget(Paragraph::new(Text::from(out_lines)), out_area);

    if has_err && err_h > 0 {
        // The last lines that fill the pane once wrapped, so a long line is never cut off.
        let err_start = stderr.len().saturating_sub(err_h as usize);
        let sep_w = inner.width.saturating_sub(12) as usize;
        let sep = format!(" {} stderr ", "─".repeat(sep_w));
        let mut err_items: Vec<Line> = vec![Line::from(vec![
            Span::raw("  "),
            Span::styled(
                sep,
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            ),
        ])];
        let mut err_rows: Vec<Line<'static>> = Vec::new();
        for l in stderr[err_start..].iter() {
            let c = if l.starts_with("error") || l.starts_with("Error") {
                Color::Red
            } else if l.contains("phase") {
                Color::Cyan
            } else {
                Color::Yellow
            };
            err_rows.extend(wrap_text(l, max_w, 2, Style::default().fg(c)));
        }
        let keep = (err_h as usize).saturating_sub(1);
        err_items.extend(err_rows.split_off(err_rows.len().saturating_sub(keep)));

        let err_area = Rect {
            x: inner.x,
            y: inner.y + out_h,
            width: inner.width,
            height: err_h,
        };
        f.render_widget(Paragraph::new(Text::from(err_items)), err_area);
    }

    render_footer(f, layout[2], foot);
}

/// The results section a metric key belongs to, for the keys the subcommands print (build-spss, build-fm, dbruijn and spss-update).
fn classify_section(key: &str) -> &'static str {
    match key {
        "build_s" | "rss_mb" => "Build",
        "idx_mb" | "bits_per" | "distinct_kmers" | "k" | "m" | "minimizer_len" | "superkmers"
        | "kmer_positions" | "super_len" | "bwt_len" => "Index",
        "unitigs" | "links" | "kmers" | "total_length" | "max_unitig" | "min_unitig"
        | "avg_unitig" | "n50" => "Graph",
        "kmers_before" | "inserted" | "erased" | "delta_added" | "delta_deleted"
        | "kmers_after" | "edit_s" | "total_s" => "Update",
        _ => "Results",
    }
}

/// `n` with thousands separators ("1,234,567").
fn fmt_thousands(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// The display label of a metric key.
fn pretty_key(key: &str) -> String {
    match key {
        "build_s" => "Build time".into(),
        "rss_mb" => "Peak RSS".into(),
        "idx_mb" => "Index size".into(),
        "bits_per" => "Bits per k-mer".into(),
        "edit_s" => "Edit time".into(),
        "total_s" => "Total time".into(),
        _ => key.to_string(),
    }
}

/// A " ── name ───" section rule, `max_w` columns wide.
fn section_header(name: &str, max_w: usize) -> Line<'static> {
    // " ── " + name + " " + dashes fits max_w without wrapping
    let dash_len = max_w.saturating_sub(name.chars().count() + 6);
    let dashes = "─".repeat(dash_len);
    let title = name.to_string();
    Line::from(vec![
        Span::styled(
            " ── ",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(
            title,
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(format!(" {}", dashes), Style::default().fg(Color::DarkGray)),
    ])
}

/// One "label: value" results line.
fn render_kv_line(key: &str, value: &str, val_color: Color) -> Line<'static> {
    let label = pretty_key(key);
    let pad = 18usize.saturating_sub(label.len());
    Line::from(vec![
        Span::styled(
            format!("  {}", label),
            Style::default()
                .fg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled(":".to_string(), Style::default().fg(Color::DarkGray)),
        Span::styled(" ".repeat(pad), Style::default()),
        Span::styled(value.to_string(), Style::default().fg(val_color)),
    ])
}

/// The value of metric `key` with its unit, and counts with thousands separators.
fn format_metric(key: &str, v: &str) -> String {
    match key {
        "build_s" | "edit_s" | "total_s" => format!("{v} s"),
        "rss_mb" | "idx_mb" => format!("{v} MB"),
        "bits_per" => format!("{v} bits/k-mer"),
        "distinct_kmers" | "kmer_positions" | "superkmers" | "unitigs" | "links" | "kmers"
        | "total_length" | "max_unitig" | "min_unitig" | "n50" | "super_len" | "bwt_len"
        | "kmers_before" | "inserted" | "erased" | "delta_added" | "delta_deleted"
        | "kmers_after" => v
            .parse::<u64>()
            .map(fmt_thousands)
            .unwrap_or_else(|_| v.to_string()),
        _ => v.to_string(),
    }
}

/// The tab-separated metrics a subcommand printed, grouped into sections.
fn format_metrics(tsv: &[(String, String)], max_w: usize) -> Vec<Line<'static>> {
    let mut lines: Vec<Line<'static>> = Vec::new();
    let mut sections: Vec<(&str, Vec<(String, String)>)> = Vec::new();
    for (k, v) in tsv {
        // one block per section, in first-seen order
        let sec = classify_section(k);
        if let Some((_, pairs)) = sections.iter_mut().find(|(s, _)| *s == sec) {
            pairs.push((k.clone(), v.clone()));
        } else {
            sections.push((sec, vec![(k.clone(), v.clone())]));
        }
    }
    for (si, (sec_name, pairs)) in sections.iter().enumerate() {
        if si > 0 {
            lines.push(Line::from(Span::raw("")));
        }
        lines.push(section_header(sec_name, max_w));
        let sec_color = match *sec_name {
            "Build" | "Graph" => Color::Green,
            "Index" => Color::Magenta,
            "Update" => Color::Yellow,
            _ => Color::White,
        };
        for (k, v) in pairs {
            // Wrapped values continue under the value column.
            lines.extend(wrap_line(
                render_kv_line(k, &format_metric(k, v.trim()), sec_color),
                max_w,
                21,
            ));
        }
    }
    lines
}

/// `s` cut to at most `max` characters, ending in ".." when cut (char-safe, since tool output and file names may hold multi-byte text).
fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() <= max {
        return s.to_string();
    }
    let keep = max.saturating_sub(2);
    let mut t: String = s.chars().take(keep).collect();
    t.push_str("..");
    t
}

/// Columns a tab takes when output is wrapped for the screen
const TAB_COLS: usize = 2;

/// wrap_line() - A styled line broken into rows of at most @width characters.
/// @line: the line, whose spans keep their styles across the breaks.
/// @width: the widest row, in characters.
/// @indent: spaces that start every row after the first.
///
/// A row breaks after its last space when it has one, and inside a word only when the word alone is wider than a row, so no text is ever cut off.
/// Tabs become TAB_COLS spaces.
fn wrap_line(line: Line<'static>, width: usize, indent: usize) -> Vec<Line<'static>> {
    let cells: Vec<(char, Style)> = line
        .spans
        .iter()
        .flat_map(|sp| {
            let st = sp.style;
            sp.content
                .chars()
                .flat_map(move |c| {
                    let n = if c == '\t' { TAB_COLS } else { 1 };
                    std::iter::repeat_n(if c == '\t' { ' ' } else { c }, n)
                })
                .map(move |c| (c, st))
                .collect::<Vec<_>>()
        })
        .collect();
    let width = width.max(1);
    let indent = indent.min(width / 2);
    let mut rows: Vec<Vec<(char, Style)>> = Vec::new();
    let mut i = 0;
    while i < cells.len() || rows.is_empty() {
        let room = if rows.is_empty() {
            width
        } else {
            width - indent
        };
        let rest = cells.len() - i;
        let take = if rest <= room {
            rest
        } else {
            // The last space inside the row, or a hard break when the row holds one word
            match cells[i..i + room + 1].iter().rposition(|(c, _)| *c == ' ') {
                Some(sp) if sp > 0 => sp,
                _ => room,
            }
        };
        let mut row: Vec<(char, Style)> = Vec::with_capacity(take + indent);
        if !rows.is_empty() {
            row.extend(std::iter::repeat_n((' ', Style::default()), indent));
        }
        row.extend_from_slice(&cells[i..i + take]);
        rows.push(row);
        i += take;
        // The space a row broke at is not carried to the start of the next one.
        while i < cells.len() && cells[i].0 == ' ' && take < rest {
            i += 1;
        }
        if rest == 0 {
            break;
        }
    }
    rows.into_iter()
        .map(|row| {
            let mut spans: Vec<Span<'static>> = Vec::new();
            for (c, st) in row {
                match spans.last_mut() {
                    Some(last) if last.style == st => last.content.to_mut().push(c),
                    _ => spans.push(Span::styled(c.to_string(), st)),
                }
            }
            Line::from(spans)
        })
        .collect()
}

/// Rows of plain text @s wrapped like wrap_line(), all in @style.
fn wrap_text(s: &str, width: usize, indent: usize, style: Style) -> Vec<Line<'static>> {
    wrap_line(
        Line::from(Span::styled(s.to_string(), style)),
        width,
        indent,
    )
}

/// Whether the first lines of `output` look like SAM (header lines, or a header followed by an 11-column record).
fn is_sam_output(output: &[String]) -> bool {
    let sam_headers = ["@HD\t", "@SQ\t", "@PG\t", "@RG\t", "@CO\t"];
    let header_count = output
        .iter()
        .take(20)
        .filter(|l| sam_headers.iter().any(|h| l.starts_with(h)))
        .count();
    if header_count >= 2 {
        return true;
    }
    if output.len() > 2 && header_count >= 1 {
        let non_header = output.iter().skip(2).find(|l| !l.starts_with('@'));
        if let Some(line) = non_header {
            return line.split('\t').count() >= 11;
        }
    }
    false
}

/// SAM lines in `range`, each record shortened to its main fields.
fn format_sam_lines(output: &[String], range: std::ops::Range<usize>) -> Vec<Line<'static>> {
    let mut out: Vec<Line<'static>> = Vec::new();
    for l in output.iter().skip(range.start).take(range.len()) {
        if l.starts_with('[') || l.starts_with("===") {
            continue;
        }
        if l.starts_with('@') {
            out.push(Line::from(Span::styled(
                l.clone(),
                Style::default().fg(Color::Cyan),
            )));
            continue;
        }
        let cols: Vec<&str> = l.split('\t').collect();
        if cols.len() >= 11 {
            let qname = cols[0];
            let flags: i32 = cols[1].parse().unwrap_or(0);
            let rname = if cols[2] == "*" { "*" } else { cols[2] };
            let pos = cols[3];
            let mapq = cols[4];
            let cigar = cols[5];
            let rnext = cols[6];
            let pnext = cols[7];
            let seq_preview = if cols[9].len() > 12 && cols[9].is_ascii() {
                format!("{}..{}", &cols[9][..6], &cols[9][cols[9].len() - 4..])
            } else {
                truncate(cols[9], 12)
            };

            let flag_str = if flags & SAM_FLAG_UNMAPPED != 0 {
                "UNMAPPED"
            } else if flags & SAM_FLAG_REVERSE != 0 {
                "REVERSE"
            } else {
                "FORWARD"
            };

            let info = format!(
                " {}  {}  {}:{}  {}  {}  {}  seq={}",
                qname,
                flag_str,
                rname,
                pos,
                mapq,
                cigar,
                if rnext != "*" {
                    format!("next={}:{}", rnext, pnext)
                } else {
                    String::new()
                },
                seq_preview
            );
            out.push(Line::from(Span::styled(
                info,
                Style::default().fg(Color::White),
            )));
        } else {
            out.push(Line::from(Span::styled(
                l.clone(),
                Style::default().fg(Color::DarkGray),
            )));
        }
    }
    out
}

/// Plain output lines, errors in red.
fn format_plain_lines(lines: &[String], max_w: usize) -> Vec<Line<'static>> {
    let mut out: Vec<Line<'static>> = Vec::new();
    for l in lines {
        let c = if l.starts_with("error") || l.starts_with("Error") {
            Color::Red
        } else {
            Color::White
        };
        out.extend(wrap_text(l, max_w, 2, Style::default().fg(c)));
    }
    out
}

/// A tool's output for the results screen, either SAM records or metric sections followed by the remaining plain lines.
fn format_tool_output(
    output: &[String],
    max_w: usize,
    range: std::ops::Range<usize>,
) -> Vec<Line<'static>> {
    let sample = output.iter().take(20).cloned().collect::<Vec<_>>();

    if is_sam_output(&sample) {
        return format_sam_lines(output, range)
            .into_iter()
            .flat_map(|l| wrap_line(l, max_w, 2))
            .collect();
    }

    let mut tsv_pairs: Vec<(String, String)> = Vec::new();
    let mut raw: Vec<String> = Vec::new();

    for line in output.iter().skip(range.start).take(range.len()) {
        if let Some(tab) = line.find('\t') {
            let key = line[..tab].to_string();
            let val = line[tab + 1..].to_string();
            if key.chars().all(|c| c.is_alphanumeric() || c == '_') {
                tsv_pairs.push((key, val));
                continue;
            }
        }
        raw.push(line.clone());
    }

    let mut formatted: Vec<Line<'static>> = Vec::new();

    if !tsv_pairs.is_empty() {
        formatted.extend(format_metrics(&tsv_pairs, max_w));
        if !raw.is_empty() {
            formatted.push(Line::from(Span::raw("")));
            formatted.push(section_header("Output", max_w));
            formatted.extend(format_plain_lines(&raw, max_w));
        }
    } else {
        formatted.extend(format_plain_lines(&raw, max_w));
    }

    formatted
}

/// Draws the results screen, with `formatted` caching the formatted output.
#[allow(clippy::too_many_arguments)]
fn render_results(
    f: &mut Frame,
    tool_idx: usize,
    output: &[String],
    exit_code: i32,
    scroll: &mut usize,
    formatted: &mut Vec<Line<'static>>,
    formatted_w: &mut usize,
    cmd: &str,
) {
    let tool = &TOOLS[tool_idx];
    let area = f.area();
    // The footer carries the scroll position, so its height is planned for the widest one.
    let n = output.len().max(1) * 4;
    let foot = format!("  Esc/q Back  ^/v Scroll  ? Help  100% (ln {n}-{n}/{n})  ");
    let layout = screen_layout(area, cmd, &foot);
    render_header(f, layout[0], cmd);
    f.render_widget(Clear, layout[1]);

    let ok = exit_code == 0;
    let (status_color, status_word) = if ok {
        (Color::Green, "Completed")
    } else {
        (Color::Red, "Failed")
    };
    let title_fmt = format!(
        " {} {} -- {} (exit {}) ",
        tool.icon, tool.name, status_word, exit_code
    );

    let block = Block::default()
        .title(Line::from(Span::styled(
            title_fmt,
            Style::default()
                .fg(status_color)
                .add_modifier(Modifier::BOLD),
        )))
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(Style::default().fg(status_color));
    let inner = block.inner(layout[1]);
    f.render_widget(block, layout[1]);

    let max_w = inner.width as usize;

    // The formatted rows are cached, and formatted again when the width changes, since they are wrapped to it.
    if formatted.is_empty() || *formatted_w != max_w {
        *formatted = format_tool_output(output, max_w, 0..output.len());
        *formatted_w = max_w;
    }

    let display_h = inner.height as usize;
    if *scroll > formatted.len().saturating_sub(display_h) {
        *scroll = formatted.len().saturating_sub(display_h);
    }
    let end = (*scroll + display_h).min(formatted.len());
    let visible: Vec<Line> = formatted[*scroll..end].to_vec();
    f.render_widget(Paragraph::new(Text::from(visible)), inner);

    let total = formatted.len();
    let scroll_pct = if total > display_h {
        let pct = ((*scroll + display_h) * 100) / total.max(1);
        format!(
            "  Esc/q Back  ^/v Scroll  ? Help  {}% (ln {}-{}/{})  ",
            pct.min(100),
            *scroll + 1,
            end,
            total
        )
    } else {
        "  Esc/q Back  ? Help  ".to_string()
    };
    render_footer(f, layout[2], &scroll_pct);
}

/// Narrowest description column of the help overlay, below which descriptions go under their labels.
const HELP_MIN_DESC_W: usize = 20;

/// help_entry() - One help overlay entry, a label in @label_w columns and its @desc wrapped beside it.
/// @label: the label's spans, padded to @label_w by the caller.
/// @width: the overlay's inner width.
///
/// When the pane is too narrow for a description column, the description wraps on its own rows under the label.
fn help_entry(
    label: Vec<Span<'static>>,
    label_w: usize,
    desc: &str,
    width: usize,
) -> Vec<Line<'static>> {
    if width < label_w + HELP_MIN_DESC_W {
        let mut out = wrap_line(Line::from(label), width, 4);
        out.extend(wrap_text(
            &format!("    {desc}"),
            width,
            4,
            Style::default(),
        ));
        return out;
    }
    let mut spans = label;
    spans.push(Span::raw(desc.to_string()));
    wrap_line(Line::from(spans), width, label_w)
}

/// Draws the help overlay over the current screen: the keys, and for a form or browser the tool's fields with their flags, plus where to read more.
/// Every entry wraps to the pane's width, and @scroll, clamped here, is the first row shown when the help is taller than the pane.
fn render_help(f: &mut Frame, app: &App, scroll: &mut u16) {
    let area = f.area();
    let w = area.width.saturating_sub(4).min(100);
    let pane = Rect {
        x: area.x + (area.width - w) / 2,
        y: area.y + 1,
        width: w,
        height: area.height.saturating_sub(2),
    };
    let inner = w.saturating_sub(2) as usize;
    let bold = Style::default().add_modifier(Modifier::BOLD);
    let key_w = 12;
    let key = |k: &str, what: &str| {
        help_entry(
            vec![Span::styled(
                format!(" {k:<11}"),
                Style::default().fg(Color::Yellow),
            )],
            key_w,
            what,
            inner,
        )
    };
    let mut lines = vec![Line::from(Span::styled(" Keys", bold))];
    for (k, what) in [
        ("^/v Tab", "move between groups, tools and fields"),
        (
            "Enter",
            "choose a tool, browse for a file, toggle a switch, or Run",
        ),
        ("type", "fill the focused field"),
        ("a", "show every file in the file browser"),
        ("q / Esc", "go back, or stop a run"),
        (
            "? / F1",
            "open this help, where ^/v and PgUp/PgDn scroll and any other key closes it",
        ),
    ] {
        lines.extend(key(k, what));
    }
    lines.push(Line::from(""));
    let tool_idx = match &app.screen {
        Screen::ToolForm { tool_idx, .. }
        | Screen::FileBrowser { tool_idx, .. }
        | Screen::Running { tool_idx, .. }
        | Screen::Results { tool_idx, .. } => Some(*tool_idx),
        Screen::ToolSelect => None,
    };
    match tool_idx {
        Some(ti) => {
            let t = &TOOLS[ti];
            lines.extend(wrap_text(
                &format!(" {} (optik {}): {}", t.name, t.cmd, t.desc),
                inner,
                2,
                bold,
            ));
            lines.push(Line::from(Span::styled(
                "  * required, red",
                Style::default().fg(Color::DarkGray),
            )));
            for p in t.params {
                let flag = p.flag.trim_end_matches('=');
                let name = format!("{}{}", p.name, if p.required { " *" } else { "" });
                let color = if p.required { Color::Red } else { Color::Cyan };
                lines.extend(help_entry(
                    vec![
                        Span::styled(format!("  {flag:<10}"), Style::default().fg(color)),
                        Span::styled(format!("{name:<22}"), bold),
                    ],
                    34,
                    p.desc,
                    inner,
                ));
            }
            lines.push(Line::from(""));
            lines.extend(wrap_text(
                &format!(" Every flag, with examples: optik {} --help", t.cmd),
                inner,
                1,
                Style::default(),
            ));
        }
        None => {
            lines.push(Line::from(Span::styled(" Tools", bold)));
            for t in TOOLS {
                lines.extend(help_entry(
                    vec![Span::styled(
                        format!("  {:<12}", t.name),
                        Style::default().fg(Color::Cyan),
                    )],
                    14,
                    t.desc,
                    inner,
                ));
            }
            lines.push(Line::from(""));
            lines.extend(wrap_text(
                " A run shows its `$ optik ...` command, which works outside the TUI too.",
                inner,
                1,
                Style::default(),
            ));
            lines.push(Line::from(" optik --help lists every command."));
        }
    }
    let visible = pane.height.saturating_sub(2) as usize;
    let max_scroll = lines.len().saturating_sub(visible) as u16;
    *scroll = (*scroll).min(max_scroll);
    let title = if max_scroll > 0 {
        format!(
            " Help  rows {}-{} of {}, ^/v PgUp/PgDn scroll, any other key closes ",
            *scroll + 1,
            (*scroll as usize + visible).min(lines.len()),
            lines.len()
        )
    } else {
        " Help  (any key closes) ".to_string()
    };
    f.render_widget(Clear, pane);
    f.render_widget(
        Paragraph::new(lines).scroll((*scroll, 0)).block(
            Block::default()
                .title(title)
                .borders(Borders::ALL)
                .border_type(BorderType::Rounded)
                .border_style(Style::default().fg(Color::Yellow)),
        ),
        pane,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kept_output_is_capped() {
        let mut v = Vec::new();
        for i in 0..MAX_KEPT_LINES + 10 {
            push_capped(&mut v, i.to_string());
        }
        assert_eq!(v.len(), MAX_KEPT_LINES + 1);
        assert!(v.last().unwrap().contains("truncated"));
        assert_eq!(v[MAX_KEPT_LINES - 1], (MAX_KEPT_LINES - 1).to_string());
    }

    #[test]
    fn truncate_is_char_safe() {
        assert_eq!(truncate("abc", 5), "abc");
        assert_eq!(truncate("MAPQ ≥ 30 → wrong", 8), "MAPQ ≥..");
        assert_eq!(truncate("─────", 4), "──..");
    }

    #[test]
    fn every_tool_runs_a_real_subcommand() {
        for t in TOOLS {
            assert!(
                matches!(
                    t.cmd,
                    "map" | "dbruijn" | "build-spss" | "build-fm" | "query" | "spss-update"
                ),
                "{}",
                t.cmd
            );
        }
    }

    #[test]
    fn paired_map_passes_both_mates() {
        let m = TOOLS.iter().find(|t| t.cmd == "map").unwrap();
        let mut vals: Vec<String> = m.params.iter().map(|p| p.default.to_string()).collect();
        vals[0] = "ref.fa".into();
        vals[1] = "r1.fq".into();
        vals[2] = "r2.fq".into();
        let args = build_args(m, &vals);
        assert_eq!(args.iter().filter(|a| *a == "-q").count(), 2, "{args:?}");
    }

    #[test]
    fn map_defaults_leave_k_to_the_reads() {
        // map chooses k from the reads only when -k is absent, so a default here would pin it to one value.
        let m = TOOLS.iter().find(|t| t.cmd == "map").unwrap();
        let vals: Vec<String> = m.params.iter().map(|p| p.default.to_string()).collect();
        assert!(!build_args(m, &vals).iter().any(|a| a == "-k"));
    }

    #[test]
    fn flags_and_joined_values_build_single_arguments() {
        let q = TOOLS.iter().find(|t| t.cmd == "query").unwrap();
        let mut vals: Vec<String> = q.params.iter().map(|p| p.default.to_string()).collect();
        let nv = q
            .params
            .iter()
            .position(|p| p.flag == "--no-verify")
            .unwrap();
        assert!(!build_args(q, &vals).iter().any(|a| a == "--no-verify"));
        vals[nv] = FLAG_ON.into();
        let args = build_args(q, &vals);
        let at = args.iter().position(|a| a == "--no-verify").unwrap();
        assert!(
            args.get(at + 1).is_none_or(|a| a.starts_with('-')),
            "{args:?}"
        );

        let m = TOOLS.iter().find(|t| t.cmd == "map").unwrap();
        let mut vals: Vec<String> = m.params.iter().map(|p| p.default.to_string()).collect();
        let cs = m.params.iter().position(|p| p.flag == "--cs=").unwrap();
        vals[cs] = "long".into();
        assert!(build_args(m, &vals).iter().any(|a| a == "--cs=long"));
    }

    #[test]
    fn every_form_flag_is_accepted_by_its_subcommand() {
        for t in TOOLS {
            let (vals, bools) = crate::subcommand_flags(t.cmd).expect(t.cmd);
            for p in t.params {
                let ok = match p.kind {
                    ParamKind::Flag => bools.contains(&p.flag),
                    _ if p.flag.ends_with('=') => ["short", "long"]
                        .iter()
                        .all(|v| bools.contains(&format!("{}{v}", p.flag).as_str())),
                    _ => vals.contains(&p.flag),
                };
                assert!(ok, "{} {}", t.cmd, p.flag);
            }
        }
    }

    /// The text of every cell of a frame drawn by `draw` at `w` x `h`.
    fn draw_text(w: u16, h: u16, draw: impl FnOnce(&mut Frame)) -> String {
        let mut term = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
        term.draw(draw).unwrap();
        term.backend()
            .buffer()
            .content()
            .iter()
            .map(|c| c.symbol())
            .collect()
    }

    #[test]
    fn selection_screen_shows_every_tool() {
        let app = App::new();
        let text = draw_text(160, 60, |f| render_tool_select(f, &app));
        for t in TOOLS {
            assert!(text.contains(t.name), "{}", t.name);
        }
    }

    #[test]
    fn every_form_draws_on_a_small_terminal() {
        for (ti, t) in TOOLS.iter().enumerate() {
            let vals: Vec<String> = t.params.iter().map(|p| p.default.to_string()).collect();
            for focus in 0..=t.params.len() {
                draw_text(80, 24, |f| render_tool_form(f, ti, &vals, focus, None));
            }
        }
    }

    #[test]
    fn required_fields_match_the_cli() {
        // The flags each subcommand refuses to run without, from main's dispatch.
        let cli: &[(&str, &[&str])] = &[
            ("map", &["-r", "-q"]),
            ("dbruijn", &["-i"]),
            ("build-spss", &["-i", "-o"]),
            ("build-fm", &["-i", "-o"]),
            ("query", &["-i", "-q"]),
            ("spss-update", &["-i", "-o"]),
        ];
        for t in TOOLS {
            let want = cli.iter().find(|(c, _)| *c == t.cmd).expect(t.cmd).1;
            let mut got: Vec<&str> = t
                .params
                .iter()
                .filter(|p| p.required)
                .map(|p| p.flag)
                .collect();
            got.sort_unstable();
            let mut want = want.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "{}", t.cmd);
        }
    }

    #[test]
    fn form_marks_required_fields_and_previews_the_command() {
        let ti = TOOLS.iter().position(|t| t.cmd == "query").unwrap();
        let t = &TOOLS[ti];
        let mut vals: Vec<String> = t.params.iter().map(|p| p.default.to_string()).collect();
        assert!(!missing_required(t, &vals).is_empty());
        for (p, v) in t.params.iter().zip(vals.iter_mut()) {
            if p.flag == "-i" {
                *v = "idx.bin".into();
            } else if p.flag == "-q" {
                *v = "reads.fq".into();
            }
        }
        assert!(missing_required(t, &vals).is_empty());
        let text = draw_text(80, 24, |f| {
            render_tool_form(f, ti, &vals, 0, Some("Fill in first: x"))
        });
        assert!(text.contains(" *"), "required marker");
        assert!(text.contains("$ optik query"), "preview");
        assert!(text.contains("Fill in first"), "notice");
        assert!(command_line(t, &vals).contains("-i idx.bin"));
    }

    #[test]
    fn help_overlay_draws_on_every_screen() {
        let mut app = App::new();
        app.show_help = true;
        let text = draw_text(80, 24, |f| render(f, &mut app));
        assert!(text.contains("Help"), "{text}");
        let vals: Vec<String> = TOOLS[0]
            .params
            .iter()
            .map(|p| p.default.to_string())
            .collect();
        app.screen = Screen::ToolForm {
            tool_idx: 0,
            values: vals,
            focus: 0,
        };
        let text = draw_text(80, 24, |f| render(f, &mut app));
        assert!(text.contains(TOOLS[0].cmd), "{text}");
    }

    #[test]
    fn results_screen_shows_the_command_that_ran() {
        let mut app = App::new();
        app.cmd_line = "optik query -i idx.bin -q reads.fq".into();
        app.screen = Screen::Results {
            tool_idx: 0,
            output: vec!["done".into()],
            exit_code: 0,
            scroll: 0,
            formatted: Vec::new(),
            formatted_w: 0,
        };
        let text = draw_text(80, 24, |f| render(f, &mut app));
        assert!(text.contains("$ optik query -i idx.bin"), "{text}");
    }

    /// The rows of a frame drawn by `draw` at `w` x `h`, with trailing spaces removed.
    fn draw_rows(w: u16, h: u16, draw: impl FnOnce(&mut Frame)) -> Vec<String> {
        let text = draw_text(w, h, draw);
        let cells: Vec<char> = text.chars().collect();
        cells
            .chunks(w as usize)
            .map(|r| r.iter().collect::<String>().trim_end().to_string())
            .collect()
    }

    /// The words of `s`, for checking that wrapped text lost nothing.
    fn words(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_string).collect()
    }

    #[test]
    fn wrap_line_keeps_every_word_within_the_width() {
        let text = "query --ids gives each distinct k-mer one id, but the ids are positions in [0, 12345678), averyveryverylongwordthatmustbreakinside";
        for width in [8usize, 20, 33, 80] {
            let rows = wrap_text(text, width, 2, Style::default());
            for r in &rows {
                assert!(r.width() <= width, "row wider than {width}: {r:?}");
            }
            let joined: String = rows
                .iter()
                .map(|l| {
                    l.spans
                        .iter()
                        .map(|s| s.content.to_string())
                        .collect::<String>()
                })
                .collect::<Vec<_>>()
                .join(" ");
            assert_eq!(
                words(&joined).concat(),
                words(text).concat(),
                "width {width}"
            );
        }
        assert_eq!(wrap_text("", 10, 2, Style::default()).len(), 1);
    }

    #[test]
    fn dbruijn_is_selected_with_d() {
        let d = tool_by_icon_char('d').map(|i| TOOLS[i].cmd);
        assert_eq!(d, Some("dbruijn"));
        assert_eq!(tool_by_icon_char('g'), None);
    }

    #[test]
    fn long_output_lines_wrap_on_the_results_screen() {
        // The build-spss duplicate warning is longer than an 80-column screen.
        let warning = "[build-spss] warning: 120 of 4554269 indexed k-mer positions repeat a k-mer (4554149 distinct); query --ids gives each distinct k-mer one id, but the ids are positions in [0, 4554269), not a bijection onto [0, 4554149).  Use --compact for a duplicate-free index.";
        let mut app = App::new();
        app.screen = Screen::Results {
            tool_idx: TOOLS.iter().position(|t| t.cmd == "build-spss").unwrap(),
            output: vec!["build_s\t0.25".into(), warning.into()],
            exit_code: 0,
            scroll: 0,
            formatted: Vec::new(),
            formatted_w: 0,
        };
        let rows = draw_rows(80, 24, |f| render(f, &mut app));
        let shown = rows.join(" ").replace('│', " ");
        for w in ["bijection", "duplicate-free", "index."] {
            assert!(shown.contains(w), "{w} missing:\n{}", rows.join("\n"));
        }
    }

    #[test]
    fn form_shows_a_long_field_description_in_full() {
        let ti = TOOLS.iter().position(|t| t.cmd == "map").unwrap();
        let t = &TOOLS[ti];
        let focus = t.params.iter().position(|p| p.flag == "-w").unwrap();
        let vals: Vec<String> = t.params.iter().map(|p| p.default.to_string()).collect();
        let rows = draw_rows(80, 30, |f| render_tool_form(f, ti, &vals, focus, None));
        let shown = words(&rows.join(" ").replace('│', " ")).join(" ");
        let want = words(t.params[focus].desc).join(" ");
        assert!(
            shown.contains(&want),
            "description cut:\n{}",
            rows.join("\n")
        );
    }

    #[test]
    fn help_overlay_scrolls_to_its_last_row() {
        let mut app = App::new();
        let ti = TOOLS.iter().position(|t| t.cmd == "map").unwrap();
        let vals: Vec<String> = TOOLS[ti]
            .params
            .iter()
            .map(|p| p.default.to_string())
            .collect();
        app.screen = Screen::ToolForm {
            tool_idx: ti,
            values: vals,
            focus: 0,
        };
        app.show_help = true;
        app.help_scroll = u16::MAX;
        let rows = draw_rows(80, 24, |f| render(f, &mut app));
        assert!(
            rows.join("\n").contains("optik map --help"),
            "{}",
            rows.join("\n")
        );
        assert!(app.help_scroll > 0 && app.help_scroll < u16::MAX);
    }

    #[test]
    fn a_long_command_wraps_in_the_header() {
        let mut app = App::new();
        app.cmd_line = format!(
            "optik map -r {} -q reads.fq",
            "/data/".repeat(20) + "ref.fa"
        );
        app.screen = Screen::Results {
            tool_idx: 0,
            output: vec!["done".into()],
            exit_code: 0,
            scroll: 0,
            formatted: Vec::new(),
            formatted_w: 0,
        };
        let rows = draw_rows(80, 24, |f| render(f, &mut app));
        let shown = rows.join("");
        assert!(shown.contains("ref.fa -q reads.fq"), "{}", rows.join("\n"));
    }

    #[test]
    fn selection_screen_scrolls_to_the_selected_group() {
        let mut app = App::new();
        let groups = grouped_tools();
        for (gi, (_, indices)) in groups.iter().enumerate() {
            app.selected_cat = gi;
            let text = draw_text(80, 24, |f| render_tool_select(f, &app));
            for &ti in indices {
                assert!(
                    text.contains(TOOLS[ti].name),
                    "group {gi}: {} not drawn",
                    TOOLS[ti].name
                );
            }
        }
    }

    #[test]
    fn metrics_of_every_subcommand_get_units() {
        let tsv: Vec<(String, String)> = [
            ("build_s", "1.5"),
            ("idx_mb", "2"),
            ("kmers_after", "1234567"),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
        let text: String = format_metrics(&tsv, 60)
            .iter()
            .map(|l| {
                l.spans
                    .iter()
                    .map(|s| s.content.as_ref())
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            text.contains("1.5 s") && text.contains("2 MB") && text.contains("1,234,567"),
            "{text}"
        );
        assert!(text.contains("Update"), "{text}");
    }

    #[test]
    fn query_defaults_leave_k_to_the_index() {
        // A fixed default -k made every TUI query of an index built at another k fail, and empty values are not passed.
        let q = TOOLS.iter().find(|t| t.name == "Query").unwrap();
        let vals: Vec<String> = q.params.iter().map(|p| p.default.to_string()).collect();
        let args = build_args(q, &vals);
        assert!(!args.iter().any(|a| a == "-k"), "{args:?}");
    }
}
