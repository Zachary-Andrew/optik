//! Formats tool output for the results view.
//! Metrics become aligned tables, SAM records are pretty-printed, and long lines are wrapped or truncated to fit.

use super::*;

/// Wrapped messages, errors and field help continue under their text, after its two-space lead.
pub(super) const MESSAGE_WRAP_INDENT: usize = 2;

/// Digits per thousands group.
const THOUSANDS_GROUP: usize = 3;

/// A section header spends these cells on " ── " before its name and a space after it.
const SECTION_HEADER_CHROME: usize = 6;

/// Width key labels are padded to, and the hanging indent of wrapped key-value lines (under the value column).
const KV_LABEL_W: usize = 18;
const KV_WRAP_INDENT: usize = 21;

/// Characters of the ".." that truncate appends.
const ELLIPSIS_LEN: usize = 2;

/// Sniffing for SAM inspects this many leading lines, accepts on enough header lines alone, and skips a few lines before looking for the first record.
const SAM_SNIFF_LINES: usize = 20;
const SAM_SNIFF_MIN_HEADERS: usize = 2;
const SAM_SNIFF_SKIP_LINES: usize = 2;

/// SEQ longer than SEQ_PREVIEW_MAX is shown as its first SEQ_PREVIEW_HEAD and last SEQ_PREVIEW_TAIL bases.
const SEQ_PREVIEW_MAX: usize = 12;
const SEQ_PREVIEW_HEAD: usize = 6;
const SEQ_PREVIEW_TAIL: usize = 4;

/// The results section a metric key belongs to, for the keys the subcommands print (build-spss, build-fm, dbruijn and spss-update).
pub(super) fn classify_section(key: &str) -> &'static str {
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

/// `value` with thousands separators ("1,234,567").
pub(super) fn fmt_thousands(value: u64) -> String {
    let digits = value.to_string();
    let mut result = String::with_capacity(digits.len() + digits.len() / THOUSANDS_GROUP);

    for (i, digit) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(THOUSANDS_GROUP) {
            result.push(',');
        }

        result.push(digit);
    }

    result
}

/// The display label of a metric key.
pub(super) fn pretty_key(key: &str) -> String {
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
pub(super) fn section_header(name: &str, max_w: usize) -> Line<'static> {
    // The rule " ── " + name + " " + dashes fits max_w without wrapping.
    let dash_len = max_w.saturating_sub(name.chars().count() + SECTION_HEADER_CHROME);
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
pub(super) fn render_kv_line(key: &str, value: &str, val_color: Color) -> Line<'static> {
    let label = pretty_key(key);
    let pad = KV_LABEL_W.saturating_sub(label.len());

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
pub(super) fn format_metric(key: &str, value: &str) -> String {
    match key {
        "build_s" | "edit_s" | "total_s" => format!("{value} s"),
        "rss_mb" | "idx_mb" => format!("{value} MB"),
        "bits_per" => format!("{value} bits/k-mer"),
        "distinct_kmers" | "kmer_positions" | "superkmers" | "unitigs" | "links" | "kmers"
        | "total_length" | "max_unitig" | "min_unitig" | "n50" | "super_len" | "bwt_len"
        | "kmers_before" | "inserted" | "erased" | "delta_added" | "delta_deleted"
        | "kmers_after" => value
            .parse::<u64>()
            .map(fmt_thousands)
            .unwrap_or_else(|_| value.to_string()),
        _ => value.to_string(),
    }
}

/// The tab-separated metrics a subcommand printed, grouped into sections.
pub(super) fn format_metrics(tsv: &[(String, String)], max_w: usize) -> Vec<Line<'static>> {
    let mut lines: Vec<Line<'static>> = Vec::new();
    let mut sections: Vec<(&str, Vec<(String, String)>)> = Vec::new();

    for (key, value) in tsv {
        // One block per section, in first-seen order.
        let sec = classify_section(key);

        if let Some((_, pairs)) = sections.iter_mut().find(|(section, _)| *section == sec) {
            pairs.push((key.clone(), value.clone()));
        } else {
            sections.push((sec, vec![(key.clone(), value.clone())]));
        }
    }

    for (section_idx, (sec_name, pairs)) in sections.iter().enumerate() {
        if section_idx > 0 {
            lines.push(Line::from(Span::raw("")));
        }

        lines.push(section_header(sec_name, max_w));
        let sec_color = match *sec_name {
            "Build" | "Graph" => Color::Green,
            "Index" => Color::Magenta,
            "Update" => Color::Yellow,
            _ => Color::White,
        };

        for (key, value) in pairs {
            // Wrapped values continue under the value column.
            lines.extend(wrap_line(
                render_kv_line(key, &format_metric(key, value.trim()), sec_color),
                max_w,
                KV_WRAP_INDENT,
            ));
        }
    }

    lines
}

/// `text` cut to at most `max` characters, ending in ".." when cut (char-safe, since tool output and file names may hold multi-byte text).
pub(super) fn truncate(text: &str, max: usize) -> String {
    if text.chars().count() <= max {
        return text.to_string();
    }

    let keep = max.saturating_sub(ELLIPSIS_LEN);
    let mut truncated: String = text.chars().take(keep).collect();

    truncated.push_str("..");
    truncated
}

/// Columns a tab takes when output is wrapped for the screen.
pub(super) const TAB_COLS: usize = 2;

/// A styled line broken into rows of at most `width` characters.
///
/// * `line`: the line, whose spans keep their styles across the breaks.
/// * `width`: the widest row, in characters.
/// * `indent`: spaces that start every row after the first.
///
/// A row breaks after its last space when it has one, and inside a word only when the word alone is wider than a row, so no text is ever cut off.
/// Tabs become TAB_COLS spaces.
pub(super) fn wrap_line(line: Line<'static>, width: usize, indent: usize) -> Vec<Line<'static>> {
    let cells: Vec<(char, Style)> = line
        .spans
        .iter()
        .flat_map(|span| {
            let style = span.style;

            span.content
                .chars()
                .flat_map(move |ch| {
                    let repeat = if ch == '\t' { TAB_COLS } else { 1 };

                    std::iter::repeat_n(if ch == '\t' { ' ' } else { ch }, repeat)
                })
                .map(move |ch| (ch, style))
                .collect::<Vec<_>>()
        })
        .collect();
    let width = width.max(1);
    // An indent may take at most half the row, so every row keeps room for text.
    const MAX_INDENT_SHARE_DEN: usize = 2;
    let indent = indent.min(width / MAX_INDENT_SHARE_DEN);
    let mut rows: Vec<Vec<(char, Style)>> = Vec::new();
    let mut cell_pos = 0;

    while cell_pos < cells.len() || rows.is_empty() {
        let room = if rows.is_empty() {
            width
        } else {
            width - indent
        };
        let rest = cells.len() - cell_pos;
        let take = if rest <= room {
            rest
        } else {
            // The last space inside the row, or a hard break when the row holds one word.
            match cells[cell_pos..cell_pos + room + 1]
                .iter()
                .rposition(|(ch, _)| *ch == ' ')
            {
                Some(space_at) if space_at > 0 => space_at,
                _ => room,
            }
        };
        let mut row: Vec<(char, Style)> = Vec::with_capacity(take + indent);

        if !rows.is_empty() {
            row.extend(std::iter::repeat_n((' ', Style::default()), indent));
        }

        row.extend_from_slice(&cells[cell_pos..cell_pos + take]);
        rows.push(row);
        cell_pos += take;

        // The space a row broke at is not carried to the start of the next one.
        while cell_pos < cells.len() && cells[cell_pos].0 == ' ' && take < rest {
            cell_pos += 1;
        }

        if rest == 0 {
            break;
        }
    }

    rows.into_iter()
        .map(|row| {
            let mut spans: Vec<Span<'static>> = Vec::new();

            for (ch, style) in row {
                match spans.last_mut() {
                    // Runs of cells with the same style become one Span instead of one per character.
                    Some(last) if last.style == style => last.content.to_mut().push(ch),
                    _ => spans.push(Span::styled(ch.to_string(), style)),
                }
            }

            Line::from(spans)
        })
        .collect()
}

/// Rows of plain text `s` wrapped like wrap_line(), all in `style`.
pub(super) fn wrap_text(text: &str, width: usize, indent: usize, style: Style) -> Vec<Line<'static>> {
    wrap_line(
        Line::from(Span::styled(text.to_string(), style)),
        width,
        indent,
    )
}

/// Whether the first lines of `output` look like SAM (header lines, or a header followed by an 11-column record).
pub(super) fn is_sam_output(output: &[String]) -> bool {
    let sam_headers = ["@HD\t", "@SQ\t", "@PG\t", "@RG\t", "@CO\t"];
    let header_count = output
        .iter()
        .take(SAM_SNIFF_LINES)
        .filter(|line| sam_headers.iter().any(|header| line.starts_with(header)))
        .count();

    if header_count >= SAM_SNIFF_MIN_HEADERS {
        return true;
    }

    // One @-line could be anything. A record with all 11 mandatory columns after the preamble is required too.
    if output.len() > SAM_SNIFF_SKIP_LINES && header_count >= 1 {
        let non_header = output
            .iter()
            .skip(SAM_SNIFF_SKIP_LINES)
            .find(|line| !line.starts_with('@'));

        if let Some(line) = non_header {
            return line.split('\t').count() >= SAM_MANDATORY_COLUMNS;
        }
    }

    false
}

/// SAM lines in `range`, each record shortened to its main fields.
pub(super) fn format_sam_lines(output: &[String], range: std::ops::Range<usize>) -> Vec<Line<'static>> {
    let mut out: Vec<Line<'static>> = Vec::new();

    for line in output.iter().skip(range.start).take(range.len()) {
        // Progress lines tagged [name] and === banners belong to the log, so the SAM view leaves them out.
        if line.starts_with('[') || line.starts_with("===") {
            continue;
        }

        if line.starts_with('@') {
            out.push(Line::from(Span::styled(
                line.clone(),
                Style::default().fg(Color::Cyan),
            )));
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();

        if cols.len() >= SAM_MANDATORY_COLUMNS {
            let qname = cols[SAM_COL_QNAME];
            let flags: u16 = cols[SAM_COL_FLAG].parse().unwrap_or(0);
            let rname = if cols[SAM_COL_RNAME] == "*" {
                "*"
            } else {
                cols[SAM_COL_RNAME]
            };
            let pos = cols[SAM_COL_POS];
            let mapq = cols[SAM_COL_MAPQ];
            let cigar = cols[SAM_COL_CIGAR];
            let rnext = cols[SAM_COL_RNEXT];
            let pnext = cols[SAM_COL_PNEXT];
            // A long sequence keeps only its head and tail, so each record still fits on one row.
            let seq_preview =
                if cols[SAM_COL_SEQ].len() > SEQ_PREVIEW_MAX && cols[SAM_COL_SEQ].is_ascii() {
                    format!(
                        "{}..{}",
                        &cols[SAM_COL_SEQ][..SEQ_PREVIEW_HEAD],
                        &cols[SAM_COL_SEQ][cols[SAM_COL_SEQ].len() - SEQ_PREVIEW_TAIL..]
                    )
                } else {
                    truncate(cols[SAM_COL_SEQ], SEQ_PREVIEW_MAX)
                };

            // The orientation word stands in for the raw flag value, which is hard to read at a glance.
            // UNMAPPED is tested first because the reverse bit of an unmapped record has no meaning.
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
                line.clone(),
                Style::default().fg(Color::DarkGray),
            )));
        }
    }

    out
}

/// Plain output lines, errors in red.
pub(super) fn format_plain_lines(lines: &[String], max_w: usize) -> Vec<Line<'static>> {
    let mut out: Vec<Line<'static>> = Vec::new();

    for line in lines {
        let color = if line.starts_with("error") || line.starts_with("Error") {
            Color::Red
        } else {
            Color::White
        };

        out.extend(wrap_text(line, max_w, MESSAGE_WRAP_INDENT, Style::default().fg(color)));
    }

    out
}

/// A tool's output for the results screen, either SAM records or metric sections followed by the remaining plain lines.
pub(super) fn format_tool_output(
    output: &[String],
    max_w: usize,
    range: std::ops::Range<usize>,
) -> Vec<Line<'static>> {
    let sample = output
        .iter()
        .take(SAM_SNIFF_LINES)
        .cloned()
        .collect::<Vec<_>>();

    if is_sam_output(&sample) {
        return format_sam_lines(output, range)
            .into_iter()
            .flat_map(|line| wrap_line(line, max_w, MESSAGE_WRAP_INDENT))
            .collect();
    }

    // Output that is not SAM is split into key-tab-value metrics, shown as a table, and everything else, shown as plain lines after it.
    let mut tsv_pairs: Vec<(String, String)> = Vec::new();
    let mut raw: Vec<String> = Vec::new();

    for line in output.iter().skip(range.start).take(range.len()) {
        if let Some(tab) = line.find('\t') {
            let key = line[..tab].to_string();
            let val = line[tab + 1..].to_string();

            // A metric line is an identifier, a tab and a value, as the subcommands print them, while other tabbed lines are shown raw.
            if key.chars().all(|ch| ch.is_alphanumeric() || ch == '_') {
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
