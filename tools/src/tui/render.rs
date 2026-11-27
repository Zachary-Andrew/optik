//! Draws every screen of the TUI, from the tool picker and file browser to the form, the running view, results and help.
//! The header and footer they share are drawn here too.

use super::*;

/// Cells taken by a bordered block's two borders (left and right, or top and bottom).
const BORDER_CELLS: u16 = 2;

/// Rows of a tool card besides its description, and the card's minimum height.
const CARD_CHROME_ROWS: u16 = 4;
const MIN_CARD_H: u16 = 6;

/// The form and run screens put the content between two side margins, sized in percent of the width.
const SIDE_MARGIN_PCT: u16 = 10;
const CONTENT_PCT: u16 = 80;

/// Hanging indent of wrapped descriptions and command lines.
const WRAP_INDENT: usize = 4;

/// Width of the "$ " prompt before a command line, which its wrapped lines are indented by.
const PROMPT_CELLS: usize = 2;

/// The form table gives its type-tag column a fixed width, the name column 2/5 of the rest, and leaves gaps between the three columns.
const TAG_COL_W: u16 = 10;
const NAME_COL_SHARE_NUM: u16 = 2;
const NAME_COL_SHARE_DEN: u16 = 5;
const TABLE_COL_GAPS: u16 = 2;

/// Width of the centred run button.
const RUN_BUTTON_W: u16 = 20;

/// The spinner advances one of its SPINNER_FRAMES frames every SPINNER_FRAME_MS milliseconds.
const SPINNER_FRAME_MS: u128 = 100;
const SPINNER_FRAMES: u128 = 10;

/// The stderr pane takes 1/ERR_PANE_DIV of the output area, between MIN_ERR_ROWS and MAX_ERR_ROWS.
const ERR_PANE_DIV: u16 = 4;
const MIN_ERR_ROWS: u16 = 2;
const MAX_ERR_ROWS: u16 = 6;

/// A redraw formats FORMAT_LINES_PER_ROW lines per visible row plus FORMAT_SLACK_LINES, which is enough to fill the pane after wrapping.
const FORMAT_LINES_PER_ROW: usize = 4;
const FORMAT_SLACK_LINES: usize = 64;

/// Cells of the stderr separator taken by its " stderr " label and frame.
const STDERR_LABEL_CELLS: u16 = 12;

/// The footer's width is planned for line numbers up to this many times the output length (wrapped rows per line).
const FOOTER_ROWS_PER_LINE: usize = 4;

/// Multiplier from a fraction to a percentage, and the largest percentage.
const PERCENT: usize = 100;

/// The help pane leaves some cells free around itself, stops at a largest width, and pads key, parameter and tool labels to fixed widths.
const HELP_PANE_MARGIN: u16 = 4;
const HELP_PANE_MAX_W: u16 = 100;
const HELP_KEY_LABEL_W: usize = 12;
const HELP_PARAM_LABEL_W: usize = 34;
const HELP_TOOL_LABEL_W: usize = 14;

/// A key label starts with this many spaces, and the keys fill the rest of HELP_KEY_LABEL_W.
const HELP_KEY_LEAD: usize = 1;
const HELP_KEY_FIELD_W: usize = HELP_KEY_LABEL_W - HELP_KEY_LEAD;

/// Parameter and tool labels start with the two spaces written before them.
const HELP_ITEM_INDENT: usize = 2;

/// A parameter label holds the flag in a fixed field and the name in the rest.
const HELP_FLAG_FIELD_W: usize = 10;
const HELP_PARAM_NAME_W: usize = HELP_PARAM_LABEL_W - HELP_ITEM_INDENT - HELP_FLAG_FIELD_W;
const HELP_TOOL_NAME_W: usize = HELP_TOOL_LABEL_W - HELP_ITEM_INDENT;

/// A centred element splits the leftover width between its two sides.
const SIDES: u16 = 2;

/// Draws the current screen.
pub(super) fn render(frame: &mut Frame, app: &mut App) {
    // Results needs mutable access to clamp its scroll and to cache formatted lines for the current width.
    match &mut app.screen {
        Screen::ToolSelect => render_tool_select(frame, app),
        Screen::FileBrowser {
            tool_idx,
            param_idx,
            cwd,
            entries,
            cursor,
            ..
        } => render_file_browser(frame, *tool_idx, *param_idx, cwd, entries, *cursor),
        Screen::ToolForm {
            tool_idx,
            values,
            focus,
        } => render_tool_form(frame, *tool_idx, values, *focus, app.notice.as_deref()),
        Screen::Running {
            tool_idx,
            output,
            stderr,
            start,
            ..
        } => render_running(frame, *tool_idx, output, stderr, *start, &app.cmd_line),
        Screen::Results {
            tool_idx,
            output,
            exit_code,
            scroll,
            formatted,
            formatted_w,
        } => render_results(
            frame,
            *tool_idx,
            output,
            *exit_code,
            scroll,
            formatted,
            formatted_w,
            &app.cmd_line,
        ),
    }

    // The help overlay is drawn last, so it covers whichever screen lies under it.
    if app.show_help {
        let mut scroll = app.help_scroll;

        render_help(frame, app, &mut scroll);
        app.help_scroll = scroll;
    }
}

/// The title bar shown on every screen.
pub(super) fn header() -> Block<'static> {
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

/// Rows the header takes at `width`, three, or more when the command line `cmd` needs them.
pub(super) fn header_rows(cmd: &str, width: u16) -> u16 {
    BORDER_CELLS + cmd_rows(cmd, width).len().max(1) as u16
}

/// The command line of a run as wrapped `$ optik ...` rows for the header, none when `cmd` is empty.
pub(super) fn cmd_rows(cmd: &str, width: u16) -> Vec<Line<'static>> {
    if cmd.is_empty() {
        return Vec::new();
    }

    let line = Line::from(vec![
        Span::styled("$ ", Style::default().fg(Color::DarkGray)),
        Span::styled(cmd.to_string(), Style::default().fg(Color::Cyan)),
    ]);

    wrap_line(
        line,
        width.saturating_sub(BORDER_CELLS) as usize,
        PROMPT_CELLS,
    )
}

/// Draws the title bar in `area`, with the run's command line `cmd` under the title when it is not empty.
pub(super) fn render_header(frame: &mut Frame, area: Rect, cmd: &str) {
    frame.render_widget(header(), area);
    let rows = cmd_rows(cmd, area.width);

    if !rows.is_empty() {
        let body = Rect {
            x: area.x + 1,
            y: area.y + 1,
            width: area.width.saturating_sub(BORDER_CELLS),
            height: area.height.saturating_sub(BORDER_CELLS),
        };

        frame.render_widget(Paragraph::new(rows), body);
    }
}

/// Rows the key-help footer `text` takes at `width`.
pub(super) fn footer_rows(text: &str, width: u16) -> u16 {
    // footer() pads the text with one space on each side.
    const FOOTER_PAD_CELLS: usize = 2;

    if text.chars().count() + FOOTER_PAD_CELLS <= width as usize {
        1
    } else {
        wrap_text(
            text,
            width.saturating_sub(BORDER_CELLS) as usize,
            0,
            Style::default(),
        )
        .len() as u16
    }
}

/// Draws the key-help footer `text` in `area`, as one rule-line title when it fits and wrapped centered rows when it does not.
pub(super) fn render_footer(frame: &mut Frame, area: Rect, text: &str) {
    if area.height <= 1 {
        frame.render_widget(footer(text), area);
        return;
    }

    let style = Style::default()
        .fg(Color::Black)
        .bg(Color::DarkGray)
        .add_modifier(Modifier::BOLD);
    let rows = wrap_text(
        text,
        area.width.saturating_sub(BORDER_CELLS) as usize,
        0,
        style,
    );

    frame.render_widget(
        Paragraph::new(rows).alignment(ratatui::layout::Alignment::Center),
        area,
    );
}

/// The vertical split of a screen with `cmd` in the header and `foot` as the footer: header, body, footer.
pub(super) fn screen_layout(area: Rect, cmd: &str, foot: &str) -> ScreenAreas {
    let [header, body, footer] = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(header_rows(cmd, area.width)),
            Constraint::Min(1),
            Constraint::Length(footer_rows(foot, area.width)),
        ])
        .areas(area);

    ScreenAreas { header, body, footer }
}

/// The three bands of every screen, from top to bottom.
pub(super) struct ScreenAreas {
    pub(super) header: Rect,
    pub(super) body: Rect,
    pub(super) footer: Rect,
}

/// The middle column of `body`, between the two side margins.
fn content_column(body: Rect) -> Rect {
    let [_, content, _] = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(SIDE_MARGIN_PCT),
            Constraint::Percentage(CONTENT_PCT),
            Constraint::Percentage(SIDE_MARGIN_PCT),
        ])
        .areas(body);

    content
}

/// The key-help line at the bottom of a screen.
pub(super) fn footer(text: &str) -> Block<'static> {
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
pub(super) const MIN_CARD_W: u16 = 22;

pub(super) const MAX_CARD_W: u16 = 42;

/// Draws the tool-selection screen.
pub(super) fn render_tool_select(frame: &mut Frame, app: &App) {
    let area = frame.area();
    let foot = "^/v Group | <-/-> Tool | Enter or M/D/B/F/Q/U Select | ? Help | q Quit";
    let layout = screen_layout(area, "", foot);

    render_header(frame, layout.header, "");
    frame.render_widget(Clear, layout.body);
    let inner = layout.body;
    let groups = grouped_tools();
    const CARD_GAP: u16 = 2;
    let cat_gap = 1u16;

    // For each group, work out how many cards fit per row, how big they are, and how many rows the group takes.
    let geometry: Vec<(u16, u16, u16, u16)> = groups
        .iter()
        .map(|(_, indices)| {
            // Cards wrap onto further rows when a group has more than fit side by side.
            let n_tools = indices.len() as u16;
            let per_row = ((inner.width + CARD_GAP) / (MIN_CARD_W + CARD_GAP)).clamp(1, n_tools.max(1));
            let total_gap = CARD_GAP * per_row.saturating_sub(1);
            let avail_w = inner.width.saturating_sub(total_gap);
            let card_w = (avail_w / per_row).clamp(MIN_CARD_W, MAX_CARD_W);
            // Cards are as tall as the longest wrapped description of the group needs.
            let inner_w = card_w.saturating_sub(BORDER_CELLS) as usize;
            let desc_rows = indices
                .iter()
                .map(|&tool_idx| {
                    wrap_text(TOOLS[tool_idx].desc, inner_w, 0, Style::default()).len()
                })
                .max()
                .unwrap_or(1);
            let card_h = (CARD_CHROME_ROWS + desc_rows as u16).max(MIN_CARD_H);
            let rows = 1 + n_tools.div_ceil(per_row) * card_h + cat_gap;

            (per_row, card_w, card_h, rows)
        })
        .collect();
    // When the groups are taller than the screen, the view starts at the first group that keeps the selected one fully on screen, so every tool can be reached and seen.
    let sel = app.selected_cat.min(groups.len().saturating_sub(1));
    // The list scrolls by whole groups, so the first group shown moves down until the selected one fits.
    let mut first = 0;

    // first advances until the groups from first to sel fit in the inner height.
    while first < sel
        && geometry[first..=sel]
            .iter()
            .map(|group_geometry| group_geometry.3)
            .sum::<u16>()
            > inner.height
    {
        first += 1;
    }

    let mut cursor_y = inner.y;

    if first > 0 {
        frame.render_widget(
            Paragraph::new(Span::styled(
                format!("  ^ {first} more group(s) above"),
                Style::default().fg(Color::DarkGray),
            )),
            Rect { height: 1, ..inner },
        );
        cursor_y += 1;
    }

    // Drawing stops at the bottom edge. The loop above already guarantees the selected group is visible.
    for (group_idx, (cat, indices)) in groups.iter().enumerate().skip(first) {
        if cursor_y >= inner.y + inner.height {
            break;
        }

        let (per_row, card_w, card_h, rows) = geometry[group_idx];

        let is_sel_cat = group_idx == app.selected_cat;
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
        let header_area = Rect {
            x: inner.x,
            y: cursor_y,
            width: inner.width,
            height: 1,
        };

        frame.render_widget(cat_hdr, header_area);

        // Cards fill each group row by row, and a card that would cross the pane edge is skipped rather than cut.
        for (card_idx, &tool_idx) in indices.iter().enumerate() {
            let (row, col) = (card_idx as u16 / per_row, card_idx as u16 % per_row);
            let card_area = Rect {
                x: inner.x + col * (card_w + CARD_GAP),
                y: cursor_y + 1 + row * card_h,
                width: card_w,
                height: card_h,
            };

            if card_area.right() > inner.x + inner.width
                || card_area.bottom() > inner.y + inner.height
            {
                continue;
            }

            let selected = is_sel_cat && card_idx == app.selected_in_cat;

            render_tool_card(frame, card_area, tool_idx, selected);
        }

        cursor_y += rows;
    }

    render_footer(frame, layout.footer, foot);
}

/// Draws one tool's card on the selection screen.
pub(super) fn render_tool_card(frame: &mut Frame, area: Rect, tool_idx: usize, selected: bool) {
    let tool = &TOOLS[tool_idx];
    let border_color = if selected {
        tool.category.color()
    } else {
        Color::DarkGray
    };
    let border_style = if selected {
        Style::default()
            .fg(border_color)
            .add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(border_color)
    };
    let block = Block::default()
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(border_style);
    let inner = block.inner(area);

    frame.render_widget(block, area);

    let inner_w = inner.width as usize;
    let name = format!("{} {}", tool.icon, tool.name);
    // The name wraps above the description, with a blank row between them.
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

    frame.render_widget(Paragraph::new(lines), inner);
}

#[allow(clippy::too_many_arguments)]
/// Draws the file browser over the form.
pub(super) fn render_file_browser(
    frame: &mut Frame,
    tool_idx: usize,
    param_idx: usize,
    cwd: &Path,
    entries: &[(String, bool)],
    cursor: usize,
) {
    let tool = &TOOLS[tool_idx];
    let area = frame.area();
    let foot = "^/v Move  |  Enter Select/Open dir  |  a All files  |  ? Help  |  Esc Back";
    let layout = screen_layout(area, "", foot);

    render_header(frame, layout.header, "");
    frame.render_widget(Clear, layout.body);

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

    let fb_area = content_column(layout.body);

    let path_str = cwd.to_string_lossy().to_string();
    let block = Block::default()
        .title(title_fmt)
        .title_alignment(ratatui::layout::Alignment::Left)
        .borders(Borders::ALL)
        .border_type(BorderType::Rounded)
        .border_style(Style::default().fg(Color::Cyan));

    let inner = block.inner(fb_area);

    frame.render_widget(block, fb_area);

    let path_style = Style::default()
        .fg(Color::Cyan)
        .add_modifier(Modifier::BOLD);
    let path_lines = wrap_text(
        &format!(" {}", path_str),
        inner.width as usize,
        1,
        path_style,
    );
    // A deep path wraps onto several rows, but it always leaves at least one row for the listing.
    let path_h = (path_lines.len() as u16)
        .min(inner.height.saturating_sub(1))
        .max(1);
    let path_area = Rect {
        x: inner.x,
        y: inner.y,
        width: inner.width,
        height: path_h,
    };

    frame.render_widget(Paragraph::new(path_lines), path_area);

    let list_y = inner.y + path_h;
    let list_h = inner.height.saturating_sub(path_h);
    let list_area = Rect {
        x: inner.x,
        y: list_y,
        width: inner.width,
        height: list_h,
    };

    // Directories get a > marker and a trailing slash as well as the yellow colour, so they are readable on terminals without colour.
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

        frame.render_stateful_widget(list, list_area, &mut state);
    }

    render_footer(frame, layout.footer, foot);
}

/// Draws a tool's parameter form.
pub(super) fn render_tool_form(
    frame: &mut Frame,
    tool_idx: usize,
    values: &[String],
    focus: usize,
    notice: Option<&str>,
) {
    let tool = &TOOLS[tool_idx];
    let area = frame.area();
    let foot = "^/v Move  |  type a value  |  Enter Browse/Toggle/Run  |  F1 Help  |  Esc Back";
    let layout = screen_layout(area, "", foot);

    render_header(frame, layout.header, "");
    frame.render_widget(Clear, layout.body);

    let content = content_column(layout.body);
    let width = content.width as usize;
    let nparams = tool.params.len();
    // Each block below starts where the one above ended, so wrapped text never overlaps.
    let mut cursor_y = content.y;

    place_rows(
        frame,
        content,
        layout.body,
        wrap_text(
            &format!(" {} {}", tool.icon, tool.desc),
            width,
            WRAP_INDENT,
            Style::default().fg(Color::Gray),
        ),
        &mut cursor_y,
    );

    // Inside the table border, the name and value columns share the width left after a type tag that always fits.
    let inner_w = content.width.saturating_sub(BORDER_CELLS);
    let tag_w: u16 = TAG_COL_W;
    let name_w = inner_w * NAME_COL_SHARE_NUM / NAME_COL_SHARE_DEN;
    let val_w = inner_w.saturating_sub(name_w + tag_w + TABLE_COL_GAPS);
    let mut rows = Vec::new();
    let mut rows_h: u16 = 0;

    for (i, param) in tool.params.iter().enumerate() {
        let is_focus = focus == i;
        // An empty field shows what the tool will use in its place, and the focused one shows a cursor instead.
        let disp = if param.kind == ParamKind::Flag {
            if values[i].is_empty() { "off" } else { "on" }.to_string()
        } else if values[i].is_empty() {
            if is_focus {
                " |".into()
            } else if param.required {
                "(required)".into()
            } else if param.default.is_empty() && param.kind == ParamKind::FileRead {
                "(none)".into()
            } else if param.default.is_empty() {
                "(auto)".into()
            } else {
                format!("(default: {})", param.default)
            }
        } else {
            values[i].clone()
        };
        // The focused row is highlighted in all three of its columns.
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
        // The tag shows whether Enter opens the file browser, a text field or a toggle.
        let tag = match param.kind {
            ParamKind::FileRead => " [file]",
            ParamKind::FileWrite => " [out]",
            ParamKind::Int => " [int]",
            ParamKind::Text => " [text]",
            ParamKind::Flag => " [on/off]",
        };
        let name_lines = wrap_text(
            &format!(" {}{}", param.name, if param.required { " *" } else { "" }),
            name_w as usize,
            1,
            name_style,
        );
        let val_lines = wrap_text(&format!(" {}", disp), val_w as usize, 1, val_style);
        // A row is as tall as its longer wrapped cell, and the total decides whether the table has to scroll.
        let row_height = name_lines.len().max(val_lines.len()) as u16;

        rows_h += row_height;
        rows.push(
            Row::new(vec![
                Cell::from(Text::from(name_lines)),
                Cell::from(Text::from(val_lines)),
                Cell::from(Line::from(Span::styled(tag, tag_style))),
            ])
            .height(row_height),
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
    // When focus is on the Run button the last parameter row remains selected, so the table does not jump back to the top.
    let sel = if nparams > 0 {
        Some(focus.min(nparams.saturating_sub(1)))
    } else {
        None
    };
    let mut state = TableState::default().with_selected(sel);

    // The field help, the command preview and the notice are wrapped first, so the table gets the rows they and the Run button leave, and scrolls to the focused field when taller.
    let info_lines = if focus < nparams {
        let param = &tool.params[focus];
        let text = if param.default.is_empty() {
            format!("  {}", param.desc)
        } else {
            format!("  {} (default: {})", param.desc, param.default)
        };

        wrap_text(&text, width, MESSAGE_WRAP_INDENT, Style::default().fg(Color::DarkGray))
    } else {
        Vec::new()
    };
    // The preview shows the exact optik command the form will run.
    let cmd = Line::from(vec![
        Span::styled("  $ ", Style::default().fg(Color::DarkGray)),
        Span::styled(command_line(tool, values), Style::default().fg(Color::Cyan)),
    ]);
    let cmd_lines = wrap_line(cmd, width, WRAP_INDENT);
    let notice_lines = notice
        .map(|msg| {
            wrap_text(
                &format!("  {msg}"),
                width,
                MESSAGE_WRAP_INDENT,
                Style::default().fg(Color::Red).add_modifier(Modifier::BOLD),
            )
        })
        .unwrap_or_default();
    let below_need =
        (1 + info_lines.len() + cmd_lines.len() + notice_lines.len()) as u16 + RUN_BUTTON_H;
    let below = content.bottom().saturating_sub(cursor_y);
    // The table gets whatever height the help, command and notice lines leave, with MIN_TABLE_H as the floor.
    let table_h = (rows_h + BORDER_CELLS)
        .min(below.saturating_sub(below_need))
        .max(MIN_TABLE_H);
    let table_area = Rect {
        x: content.x,
        y: cursor_y,
        width: content.width,
        height: table_h,
    }
    .intersection(layout.body);
    frame.render_stateful_widget(&table, table_area, &mut state);
    cursor_y += table_h + 1;

    // The field help, then the command the form would run, then a notice such as missing required fields.
    place_rows(frame, content, layout.body, info_lines, &mut cursor_y);
    place_rows(frame, content, layout.body, cmd_lines, &mut cursor_y);
    place_rows(frame, content, layout.body, notice_lines, &mut cursor_y);
    let btn_area = Rect {
        x: content.x + (content.width / SIDES).saturating_sub(RUN_BUTTON_W / SIDES),
        y: cursor_y,
        width: RUN_BUTTON_W.min(content.width),
        height: RUN_BUTTON_H,
    }
    .intersection(layout.body);
    // Any focus at or past nparams is the Run button.
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

    frame.render_widget(btn, btn_area);

    render_footer(frame, layout.footer, foot);
}

/// Rows of the form's Run button, border included.
pub(super) const RUN_BUTTON_H: u16 = 3;

/// Fewest rows the form's table keeps, its border and one field.
pub(super) const MIN_TABLE_H: u16 = 3;

/// Draws `lines` in the columns of `col` from row `y` down, clipped to `bound`, and moves `y` below them.
pub(super) fn place_rows(
    frame: &mut Frame,
    col: Rect,
    bound: Rect,
    lines: Vec<Line<'static>>,
    cursor_y: &mut u16,
) {
    let height = lines.len() as u16;
    let rect = Rect {
        x: col.x,
        y: *cursor_y,
        width: col.width,
        height,
    }
    .intersection(bound);
    frame.render_widget(Paragraph::new(lines), rect);
    *cursor_y += height;
}

/// Draws the running screen with the elapsed time and the tail of the output.
pub(super) fn render_running(
    frame: &mut Frame,
    tool_idx: usize,
    output: &[String],
    stderr: &[String],
    start: Instant,
    cmd: &str,
) {
    let tool = &TOOLS[tool_idx];
    let area = frame.area();
    let foot = "Running...  q/Esc stop the run and go back  |  ? Help";
    let layout = screen_layout(area, cmd, foot);

    render_header(frame, layout.header, cmd);
    frame.render_widget(Clear, layout.body);

    let elapsed = start.elapsed().as_secs_f64();

    // The spinner advances with wall time rather than with redraws, so it keeps an even pace.
    const SPINNER: &str = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏";
    let sp_i = (start.elapsed().as_millis() / SPINNER_FRAME_MS % SPINNER_FRAMES) as usize;
    let sp_char = SPINNER.chars().nth(sp_i).unwrap_or(' ');

    // The out and err counts in the title still grow after the panes have filled.
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

    let inner = block.inner(layout.body);

    frame.render_widget(block, layout.body);

    let max_w = inner.width as usize;
    let has_err = !stderr.is_empty();

    // Standard error gets its own pane only when there is any, sized between MIN_ERR_ROWS and MAX_ERR_ROWS.
    let err_h = if has_err {
        (inner.height / ERR_PANE_DIV).clamp(MIN_ERR_ROWS, MAX_ERR_ROWS)
    } else {
        0
    };
    let out_h = inner.height.saturating_sub(err_h);

    // Only the tail is shown, since formatting all kept lines, up to MAX_KEPT_LINES of SAM, on every 50 ms redraw would stall the screen.
    let tail = output
        .len()
        .saturating_sub(FORMAT_LINES_PER_ROW * out_h as usize + FORMAT_SLACK_LINES)
        ..output.len();
    let formatted = format_tool_output(output, max_w, tail);
    let out_start = formatted.len().saturating_sub(out_h as usize);
    let out_lines: Vec<Line> = formatted.into_iter().skip(out_start).collect();

    let out_area = Rect {
        x: inner.x,
        y: inner.y,
        width: inner.width,
        height: out_h,
    };

    frame.render_widget(Paragraph::new(Text::from(out_lines)), out_area);

    if has_err && err_h > 0 {
        // The last lines that fill the pane once wrapped, so a long line is never cut off.
        let err_start = stderr.len().saturating_sub(err_h as usize);
        let sep_w = inner.width.saturating_sub(STDERR_LABEL_CELLS) as usize;
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

        for line in stderr[err_start..].iter() {
            // Lines starting with error are red, phase timings from PhaseTimer are cyan and other progress is yellow.
            let color = if line.starts_with("error") || line.starts_with("Error") {
                Color::Red
            } else if line.contains("phase") {
                Color::Cyan
            } else {
                Color::Yellow
            };

            err_rows.extend(wrap_text(line, max_w, MESSAGE_WRAP_INDENT, Style::default().fg(color)));
        }

        // One row is taken by the stderr separator, and the newest wrapped rows fill the rest.
        let keep = (err_h as usize).saturating_sub(1);

        err_items.extend(err_rows.split_off(err_rows.len().saturating_sub(keep)));

        let err_area = Rect {
            x: inner.x,
            y: inner.y + out_h,
            width: inner.width,
            height: err_h,
        };

        frame.render_widget(Paragraph::new(Text::from(err_items)), err_area);
    }

    render_footer(frame, layout.footer, foot);
}

/// Draws the results screen, with `formatted` caching the formatted output.
#[allow(clippy::too_many_arguments)]
pub(super) fn render_results(
    frame: &mut Frame,
    tool_idx: usize,
    output: &[String],
    exit_code: i32,
    scroll: &mut usize,
    formatted: &mut Vec<Line<'static>>,
    formatted_w: &mut usize,
    cmd: &str,
) {
    let tool = &TOOLS[tool_idx];
    let area = frame.area();
    // The footer shows the scroll position, so its height is planned for the widest one.
    let max_line_no = output.len().max(1) * FOOTER_ROWS_PER_LINE;
    let foot = format!(
        "  Esc/q Back  ^/v Scroll  ? Help  100% (ln {max_line_no}-{max_line_no}/{max_line_no})  "
    );
    let layout = screen_layout(area, cmd, &foot);

    render_header(frame, layout.header, cmd);
    frame.render_widget(Clear, layout.body);

    // The title colour gives the exit status at a glance, green for 0 and red otherwise.
    let succeeded = exit_code == 0;
    let (status_color, status_word) = if succeeded {
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
    let inner = block.inner(layout.body);

    frame.render_widget(block, layout.body);

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

    frame.render_widget(Paragraph::new(Text::from(visible)), inner);

    let total = formatted.len();
    let scroll_pct = if total > display_h {
        let pct = ((*scroll + display_h) * PERCENT) / total.max(1);

        format!(
            "  Esc/q Back  ^/v Scroll  ? Help  {}% (ln {}-{}/{})  ",
            pct.min(PERCENT),
            *scroll + 1,
            end,
            total
        )
    } else {
        "  Esc/q Back  ? Help  ".to_string()
    };

    render_footer(frame, layout.footer, &scroll_pct);
}

/// Narrowest description column of the help overlay, below which descriptions go under their labels.
pub(super) const HELP_MIN_DESC_W: usize = 20;

/// One help overlay entry, a label in `label_w` columns and its `desc` wrapped beside it.
///
/// * `label`: the label's spans, padded to `label_w` by the caller.
/// * `width`: the overlay's inner width.
///
/// When the pane is too narrow for a description column, the description wraps on its own rows under the label.
pub(super) fn help_entry(
    label: Vec<Span<'static>>,
    label_w: usize,
    desc: &str,
    width: usize,
) -> Vec<Line<'static>> {
    if width < label_w + HELP_MIN_DESC_W {
        let mut out = wrap_line(Line::from(label), width, WRAP_INDENT);

        out.extend(wrap_text(
            &format!("    {desc}"),
            width,
            WRAP_INDENT,
            Style::default(),
        ));
        return out;
    }

    let mut spans = label;

    spans.push(Span::raw(desc.to_string()));
    wrap_line(Line::from(spans), width, label_w)
}

/// Draws the help overlay over the current screen.
/// It always lists the keys, and on a form or in the browser it also shows the tool's fields with their flags and where to read more.
/// Every entry wraps to the pane's width, and `scroll`, clamped here, is the first row shown when the help is taller than the pane.
pub(super) fn render_help(frame: &mut Frame, app: &App, scroll: &mut u16) {
    let area = frame.area();
    let pane_w = area
        .width
        .saturating_sub(HELP_PANE_MARGIN)
        .min(HELP_PANE_MAX_W);
    // The pane is centred horizontally and leaves one row free above it.
    let pane = Rect {
        x: area.x + (area.width - pane_w) / SIDES,
        y: area.y + 1,
        width: pane_w,
        height: area.height.saturating_sub(BORDER_CELLS),
    };
    let inner = pane_w.saturating_sub(BORDER_CELLS) as usize;
    let bold = Style::default().add_modifier(Modifier::BOLD);
    let key_w = HELP_KEY_LABEL_W;
    let key = |keys: &str, what: &str| {
        help_entry(
            vec![Span::styled(
                format!(" {keys:<HELP_KEY_FIELD_W$}"),
                Style::default().fg(Color::Yellow),
            )],
            key_w,
            what,
            inner,
        )
    };
    // The key list comes first, because it applies on every screen.
    let mut lines = vec![Line::from(Span::styled(" Keys", bold))];

    for (keys, what) in [
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
        lines.extend(key(keys, what));
    }

    lines.push(Line::from(""));
    let tool_idx = match &app.screen {
        Screen::ToolForm { tool_idx, .. }
        | Screen::FileBrowser { tool_idx, .. }
        | Screen::Running { tool_idx, .. }
        | Screen::Results { tool_idx, .. } => Some(*tool_idx),
        Screen::ToolSelect => None,
    };

    // On a tool's screens the help describes that tool's fields, while on the tool list it describes every tool.
    match tool_idx {
        Some(selected_idx) => {
            let tool = &TOOLS[selected_idx];

            lines.extend(wrap_text(
                &format!(" {} (optik {}): {}", tool.name, tool.cmd, tool.desc),
                inner,
                MESSAGE_WRAP_INDENT,
                bold,
            ));
            lines.push(Line::from(Span::styled(
                "  * required, red",
                Style::default().fg(Color::DarkGray),
            )));
            // Required fields are marked with * and drawn in red, as on the form.
            // Flags such as --cs= show without their trailing =.
            for param in tool.params {
                let flag = param.flag.trim_end_matches('=');
                let name = format!("{}{}", param.name, if param.required { " *" } else { "" });
                let color = if param.required {
                    Color::Red
                } else {
                    Color::Cyan
                };

                lines.extend(help_entry(
                    vec![
                        Span::styled(format!("  {flag:<HELP_FLAG_FIELD_W$}"), Style::default().fg(color)),
                        Span::styled(format!("{name:<HELP_PARAM_NAME_W$}"), bold),
                    ],
                    HELP_PARAM_LABEL_W,
                    param.desc,
                    inner,
                ));
            }

            lines.push(Line::from(""));
            lines.extend(wrap_text(
                &format!(" Every flag, with examples: optik {} --help", tool.cmd),
                inner,
                1,
                Style::default(),
            ));
        }
        None => {
            lines.push(Line::from(Span::styled(" Tools", bold)));

            for tool in TOOLS {
                lines.extend(help_entry(
                    vec![Span::styled(
                        format!("  {:<HELP_TOOL_NAME_W$}", tool.name),
                        Style::default().fg(Color::Cyan),
                    )],
                    HELP_TOOL_LABEL_W,
                    tool.desc,
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

    // The scroll offset is clamped here, because the help's height depends on the wrap width, which the input handler does not know.
    let visible = pane.height.saturating_sub(BORDER_CELLS) as usize;
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

    frame.render_widget(Clear, pane);
    frame.render_widget(
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
