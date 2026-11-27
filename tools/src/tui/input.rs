//! Handles key presses and the file browser's state.
//! Once a form is filled in, this module also turns it into the command line that gets run.

use super::*;

/// How often the running screen polls the child's output and redraws.
const RUN_POLL_INTERVAL_MS: u64 = 50;

/// Exit code shown when a tool was killed by a signal (no status code) or failed to start.
const NO_EXIT_CODE: i32 = -1;

/// Rows of the help pane that are not scrolling text, subtracted from a page step.
const HELP_CHROME_ROWS: u16 = 4;

/// Lines moved by one Up or Down key in the output view.
const SCROLL_STEP: usize = 3;

/// Sequence files a FASTA/FASTQ field accepts, plain or gzipped.
pub(super) const SEQ_EXTS: &[&str] = &[
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
pub(super) const INDEX_EXTS: &[&str] = &[".spss", ".fm"];

/// The suffixes the browser should list for field `param` of `tool`, empty for any file (output paths, and fields that take no file).
pub(super) fn wanted_exts(tool: &ToolInfo, param: &ParamDef) -> &'static [&'static str] {
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

/// `s` as one shell word, quoted only when it needs to be.
pub(super) fn shell_quote(word: &str) -> String {
    if !word.is_empty()
        && word
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || "-_./=:@+,%".contains(ch))
    {
        word.to_string()
    } else {
        format!("'{}'", word.replace('\'', "'\\''"))
    }
}

/// The command line a form runs, `optik <cmd> <flags>`, quoted for a shell.
pub(super) fn command_line(tool: &ToolInfo, values: &[String]) -> String {
    std::iter::once("optik".to_string())
        .chain(std::iter::once(tool.cmd.to_string()))
        .chain(build_args(tool, values).iter().map(|arg| shell_quote(arg)))
        .collect::<Vec<_>>()
        .join(" ")
}

/// The required fields of a form that are still empty, by name.
pub(super) fn missing_required(tool: &ToolInfo, values: &[String]) -> Vec<&'static str> {
    tool.params
        .iter()
        .zip(values)
        .filter(|(param, value)| param.required && value.trim().is_empty())
        .map(|(param, _)| param.name)
        .collect()
}

/// The file browser's entries for `cwd`, i.e. "..", then directories, then files, each sorted with hidden names skipped, as (name, is_dir) pairs.
/// With `exts` non-empty only files ending in one of them are listed.
pub(super) fn refresh_entries(cwd: &Path, exts: &[&str]) -> Vec<(String, bool)> {
    let mut entries: Vec<(String, bool)> = Vec::new();

    if cwd.parent().is_some() {
        entries.push(("..".to_string(), true));
    }

    if let Ok(dir_entries) = std::fs::read_dir(cwd) {
        let mut dirs: Vec<String> = Vec::new();
        let mut files: Vec<String> = Vec::new();

        for entry in dir_entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();

            // Dotfiles are skipped, as most file pickers do.
            if name.starts_with('.') {
                continue;
            }

            if let Ok(file_type) = entry.file_type() {
                if file_type.is_dir() {
                    dirs.push(name);
                } else if exts.is_empty()
                    // An empty extension list comes from the a toggle and lists every file.
                    || exts
                        .iter()
                        .any(|ext| name.to_ascii_lowercase().ends_with(ext))
                {
                    files.push(name);
                }
            }
        }

        dirs.sort();
        files.sort();

        for dir_name in dirs {
            entries.push((dir_name, true));
        }

        for file_name in files {
            entries.push((file_name, false));
        }
    }

    entries
}

/// Lines of a tool's output the TUI keeps.
/// Everything is still read (the child must not block on a full pipe), but a map run without an output file streams the whole SAM here, and holding it all would take as much memory as the SAM itself.
pub(super) const MAX_KEPT_LINES: usize = 200_000;

/// Keeps `line` if fewer than MAX_KEPT_LINES are kept, and adds one truncation notice when the cap is reached.
pub(super) fn push_capped(kept: &mut Vec<String>, line: String) {
    if kept.len() < MAX_KEPT_LINES {
        kept.push(line);
    } else if kept.len() == MAX_KEPT_LINES {
        kept.push(format!(
            "... output truncated: the TUI keeps the first {MAX_KEPT_LINES} lines \
                        (set an output file to keep everything)"
        ));
    }
}

/// Collects the running child's output and, once it has exited, moves to the results screen.
pub(super) fn process_running(app: &mut App) {
    // The screen is taken out by value, because moving to Results consumes the output buffers and the child handle.
    // mem::replace takes ownership of the Running state, so its output vectors move into Results without copying.
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
            if let Some(receiver) = &rx_stdout {
                while let Ok(line) = receiver.try_recv() {
                    push_capped(&mut output, line);
                }
            }

            if let Some(receiver) = &rx_stderr {
                while let Ok(line) = receiver.try_recv() {
                    push_capped(&mut stderr, line);
                }
            }

            if let Some(running_child) = &mut child {
                match running_child.try_wait() {
                    // The pipes are drained once more after the exit, since the reader threads may have queued the last lines after the previous poll.
                    Ok(Some(status)) => {
                        if let Some(receiver) = &rx_stdout {
                            while let Ok(line) = receiver.try_recv() {
                                push_capped(&mut output, line);
                            }
                        }

                        if let Some(receiver) = &rx_stderr {
                            while let Ok(line) = receiver.try_recv() {
                                push_capped(&mut stderr, line);
                            }
                        }

                        let exit_code = status.code().unwrap_or(NO_EXIT_CODE);
                        // stdout comes before stderr in Results, so the metric lines appear before the progress log.
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
                    // Still running, so the screen is put back as it was.
                    Ok(None) => Screen::Running {
                        tool_idx,
                        output,
                        stderr,
                        start,
                        rx_stdout,
                        rx_stderr,
                        child,
                    },
                    Err(err) => {
                        output.push(format!("Process error: {err}"));
                        let mut combined = output;

                        combined.extend(stderr);
                        Screen::Results {
                            tool_idx,
                            output: combined,
                            exit_code: NO_EXIT_CODE,
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
/// While a tool runs, picks up its output and lets q or Esc stop it.
fn poll_running(app: &mut App) -> Result<(), String> {
    process_running(app);

    if matches!(app.screen, Screen::Running { .. }) {
        let timeout = std::time::Duration::from_millis(RUN_POLL_INTERVAL_MS);

        if event::poll(timeout).map_err(|err| err.to_string())?
            && let Event::Key(key) = event::read().map_err(|err| err.to_string())?
            && key.kind == KeyEventKind::Press
            && matches!(key.code, KeyCode::Char('q') | KeyCode::Esc)
        {
            // Dropping the handle would leave the child to die of SIGPIPE at its next write, so it is stopped explicitly.
            if let Screen::Running {
                child: Some(running_child),
                ..
            } = &mut app.screen
            {
                let _ = running_child.kill();
                let _ = running_child.wait();
            }

            app.screen = Screen::ToolSelect;
            app.selected_cat = 0;
            app.selected_in_cat = 0;
        }
    }

    Ok(())
}

/// Scrolls the help overlay with the arrow and page keys, and closes it on any other key.
fn handle_help_key(app: &mut App, key: KeyEvent) {
    let page = app.term_h.saturating_sub(HELP_CHROME_ROWS).max(1);

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
}

/// Moves through the file browser, enters directories and returns the picked path to the form.
fn handle_browser_key(app: &mut App, key: KeyEvent) {
    // The filter settings are copied out first, because the screen is borrowed mutably below.
    let app_browse_all = app.browse_all;
    let app_browse_exts = app.browse_exts;

    let Screen::FileBrowser {
        tool_idx,
        param_idx,
        values,
        focus,
        cwd,
        entries,
        cursor,
    } = &mut app.screen else {
        return;
    };

    match key.code {
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
        // Enter opens a directory, or picks a file and fills the field with its full path.
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
        // The a key switches between the parameter's usual file extensions and every file.
        KeyCode::Char('a') => {
            app.browse_all = !app_browse_all;
            *entries =
                refresh_entries(cwd, if app.browse_all { &[] } else { app_browse_exts });
            *cursor = 0;
        }
        _ => {}
    }
}

/// Scrolls the results view, or goes back to the tool picker on q or Esc.
fn handle_results_key(app: &mut App, key: KeyEvent) {
    let Screen::Results { scroll, .. } = &mut app.screen else {
        return;
    };

    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => {
            app.screen = Screen::ToolSelect;
            app.selected_cat = 0;
            app.selected_in_cat = 0;
        }
        KeyCode::Up | KeyCode::Char('k') => {
            *scroll = scroll.saturating_sub(SCROLL_STEP);
        }
        KeyCode::Down | KeyCode::Char('j') => {
            *scroll = scroll.saturating_add(SCROLL_STEP);
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
    }
}

/// Edits the form's fields, opens the file browser for path fields and runs the tool from the Run button.
fn handle_form_key(app: &mut App, key: KeyEvent) {
    let Screen::ToolForm {
        tool_idx,
        values,
        focus,
    } = &mut app.screen else {
        return;
    };

    let tool = &TOOLS[*tool_idx];
    let nparams = tool.params.len();

    match key.code {
        KeyCode::Esc => {
            app.screen = Screen::ToolSelect;
            app.selected_cat = 0;
            app.selected_in_cat = 0;
        }
        // No j/k here, since the fields take typed paths and values.
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
                .position(|(param, value)| param.required && value.trim().is_empty())
            {
                app.notice = Some(format!("Fill in first: {}", missing.join(", ")));
                *focus = first;
            } else {
                app.cmd_line = command_line(tool, values);
                let args = build_args(tool, values);
                let full_args: Vec<String> =
                    std::iter::once(tool.cmd.to_string()).chain(args).collect();
                // The tool runs as a child of this same binary, so the TUI always calls the build it belongs to.
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

                        // Each pipe gets its own reader thread and channel, so a busy stderr cannot block stdout and the UI thread never waits on either.
                        if let Some(out) = child.stdout.take() {
                            std::thread::spawn(move || {
                                for line in io::BufReader::new(out).lines() {
                                    if let Ok(text) = line
                                        && tx_stdout.send(text).is_err()
                                    {
                                        break;
                                    }
                                }
                            });
                        }

                        if let Some(err) = child.stderr.take() {
                            std::thread::spawn(move || {
                                for line in io::BufReader::new(err).lines() {
                                    if let Ok(text) = line
                                        && tx_stderr.send(text).is_err()
                                    {
                                        break;
                                    }
                                }
                            });
                        }

                        // The form switches to the running screen at once, and process_running() collects what the reader threads send.
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
                    // A tool that cannot start shows its error on the results screen, like any failed run.
                    Err(err) => {
                        app.screen = Screen::Results {
                            tool_idx: *tool_idx,
                            output: vec![format!("Error: failed to start tool: {err}")],
                            exit_code: NO_EXIT_CODE,
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
                // The browser opens in the field's directory or next to its file, and otherwise in the working directory.
                let cwd = if values[*focus].is_empty() {
                    std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
                } else {
                    let current = Path::new(&values[*focus]);

                    if current.is_dir() {
                        current.to_path_buf()
                    } else {
                        current
                            .parent()
                            .map(|parent| parent.to_path_buf())
                            .unwrap_or_else(|| {
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
            let value = &mut values[*focus];

            *value = if value.is_empty() {
                FLAG_ON.to_string()
            } else {
                String::new()
            };
        }

        KeyCode::Char(ch)
            if *focus < nparams && tool.params[*focus].kind != ParamKind::Flag =>
        {
            values[*focus].push(ch);
        }

        KeyCode::Backspace if *focus < nparams => {
            values[*focus].pop();
        }
        _ => {}
    }
}

pub(super) fn handle_event(app: &mut App) -> Result<bool, String> {
    if matches!(app.screen, Screen::Running { .. }) {
        poll_running(app)?;
        return Ok(true);
    }

    if let Event::Key(key) = event::read().map_err(|err| err.to_string())? {
        if key.kind != KeyEventKind::Press {
            return Ok(true);
        }

        // The arrow and page keys scroll the help overlay, and any other key closes it.
        // F1 opens it anywhere, and `?` everywhere except while a form field takes typed text.
        if app.show_help {
            handle_help_key(app, key);
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

        match &app.screen {
            Screen::ToolSelect => return Ok(handle_tool_select(app, key)),
            Screen::FileBrowser { .. } => handle_browser_key(app, key),
            Screen::Results { .. } => handle_results_key(app, key),
            Screen::ToolForm { .. } => handle_form_key(app, key),
            _ => {}
        }
    }

    Ok(true)
}

/// Opens the parameter form of tool `tool_idx` with its default values.
pub(super) fn select_tool(app: &mut App, tool_idx: usize) {
    if tool_idx < TOOLS.len() {
        let tool = &TOOLS[tool_idx];
        let values: Vec<String> = tool
            .params
            .iter()
            .map(|param| param.default.to_string())
            .collect();
        let groups = grouped_tools();

        for (group_idx, (_, indices)) in groups.iter().enumerate() {
            if let Some(pos_in_group) = indices
                .iter()
                .position(|&candidate_idx| candidate_idx == tool_idx)
            {
                app.selected_cat = group_idx;
                app.selected_in_cat = pos_in_group;
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

/// The tool whose one-letter icon is `icon_char` (case-insensitive), the keyboard shortcut on the selection screen.
pub(super) fn tool_by_icon_char(icon_char: char) -> Option<usize> {
    let c_upper = icon_char.to_ascii_uppercase();

    TOOLS.iter().position(|tool| {
        let icon = tool.icon.trim();

        icon.len() == 1
            && icon
                .chars()
                .next()
                .map(|icon_ch| icon_ch.to_ascii_uppercase() == c_upper)
                .unwrap_or(false)
    })
}

/// Returns false to quit, and the caller restores the terminal.
pub(super) fn handle_tool_select(app: &mut App, key: KeyEvent) -> bool {
    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => return false,
        // Up and Down move between categories and keep the column inside the new category, while Left and Right move within one.
        KeyCode::Up | KeyCode::Char('k') => {
            let num_cats = cat_count();

            if num_cats > 0 {
                // selected_in_cat is clamped because the new group may have fewer tools.
                app.selected_cat = app.selected_cat.saturating_sub(1);
                let nic = tools_in_cat(app.selected_cat);

                if app.selected_in_cat >= nic {
                    app.selected_in_cat = nic.saturating_sub(1);
                }
            }
        }
        KeyCode::Down | KeyCode::Char('j') => {
            let num_cats = cat_count();

            if num_cats > 0 {
                app.selected_cat = (app.selected_cat + 1).min(num_cats - 1);
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
        // Each tool's icon letter opens its form directly.
        KeyCode::Char(ch) => {
            if let Some(tool_idx) = tool_by_icon_char(ch) {
                select_tool(app, tool_idx);
            }
        }
        _ => {}
    }

    true
}

/// The command-line arguments of a filled-in form, one per non-empty field: a flag and its value, a joined `flag=value`, or a switched-on flag alone.
pub(super) fn build_args(tool: &ToolInfo, values: &[String]) -> Vec<String> {
    let mut args = Vec::new();

    for (i, param) in tool.params.iter().enumerate() {
        let value = values[i].trim();

        if value.is_empty() {
            continue;
        }

        if param.kind == ParamKind::Flag {
            args.push(param.flag.to_string());
        } else if param.flag.ends_with('=') {
            args.push(format!("{}{value}", param.flag));
        } else {
            args.push(param.flag.to_string());
            args.push(value.to_string());
        }
    }

    args
}
