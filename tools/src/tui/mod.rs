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

use crate::map::{
    SAM_COL_CIGAR, SAM_COL_FLAG, SAM_COL_MAPQ, SAM_COL_PNEXT, SAM_COL_POS, SAM_COL_QNAME, SAM_COL_RNAME, SAM_COL_RNEXT,
    SAM_COL_SEQ, SAM_FLAG_REVERSE, SAM_FLAG_UNMAPPED, SAM_MANDATORY_COLUMNS,
};

mod catalog;
mod input;
mod render;
mod format;

use catalog::*;
use input::*;
use render::*;
use format::*;

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
    /// The `?` overlay is open, and any key closes it.
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

/// Runs the TUI until the user quits, restoring the terminal afterwards, even after a panic.
pub(crate) fn run() -> Result<(), String> {
    let prev_hook = std::panic::take_hook();

    // The hook leaves raw mode and the alternate screen before the default hook prints, otherwise the panic message would be invisible.
    std::panic::set_hook(Box::new(move |info| {
        use std::io::Write;
        let _ = disable_raw_mode();
        let mut stdout = io::stdout();
        let _ = execute!(stdout, LeaveAlternateScreen);
        let _ = stdout.flush();

        prev_hook(info);
    }));

    enable_raw_mode().map_err(|err| err.to_string())?;
    let mut stdout = io::stdout();

    execute!(stdout, EnterAlternateScreen).map_err(|err| err.to_string())?;
    let mut terminal = ratatui::Terminal::new(ratatui::backend::CrosstermBackend::new(stdout))
        .map_err(|err| err.to_string())?;

    terminal.clear().map_err(|err| err.to_string())?;
    let mut app = App::new();
    // res is returned only after the terminal is restored, which lets an error message print on the normal screen.
    let res = run_app(&mut terminal, &mut app);

    disable_raw_mode().map_err(|err| err.to_string())?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen).map_err(|err| err.to_string())?;
    terminal.show_cursor().map_err(|err| err.to_string())?;
    res
}

/// The draw / handle-event loop.
fn run_app<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    app: &mut App,
) -> Result<(), String> {
    loop {
        let size = terminal.size().map_err(|err| err.to_string())?;

        app.term_w = size.width;
        app.term_h = size.height;
        terminal
            .draw(|frame| render(frame, &mut *app))
            .map_err(|err| err.to_string())?;

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
        .map(|(_, tool_indices)| tool_indices.len())
        .unwrap_or(0)
}

/// Index into TOOLS of tool `in_cat` of category `cat_idx` (0 if out of range).
fn global_tool_idx(cat_idx: usize, in_cat: usize) -> usize {
    grouped_tools()
        .get(cat_idx)
        .and_then(|(_, tool_indices)| tool_indices.get(in_cat).copied())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests;
