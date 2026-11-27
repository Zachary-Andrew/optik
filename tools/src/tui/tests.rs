use super::*;

#[test]
fn kept_output_is_capped() {
    let mut kept = Vec::new();

    for i in 0..MAX_KEPT_LINES + 10 {
        push_capped(&mut kept, i.to_string());
    }

    assert_eq!(kept.len(), MAX_KEPT_LINES + 1);
    assert!(kept.last().unwrap().contains("truncated"));
    assert_eq!(kept[MAX_KEPT_LINES - 1], (MAX_KEPT_LINES - 1).to_string());
}

#[test]
fn truncate_is_char_safe() {
    assert_eq!(truncate("abc", 5), "abc");
    assert_eq!(truncate("MAPQ ≥ 30 → wrong", 8), "MAPQ ≥..");
    assert_eq!(truncate("─────", 4), "──..");
}

#[test]
fn every_tool_runs_a_real_subcommand() {
    for tool in TOOLS {
        assert!(
            matches!(
                tool.cmd,
                "map" | "dbruijn" | "build-spss" | "build-fm" | "query" | "spss-update"
            ),
            "{}",
            tool.cmd
        );
    }
}

#[test]
fn paired_map_passes_both_mates() {
    let map_tool = TOOLS.iter().find(|tool| tool.cmd == "map").unwrap();
    let mut vals: Vec<String> = map_tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();

    vals[0] = "ref.fa".into();
    vals[1] = "r1.fq".into();
    vals[2] = "r2.fq".into();
    let args = build_args(map_tool, &vals);

    assert_eq!(
        args.iter().filter(|arg| *arg == "-q").count(),
        2,
        "{args:?}"
    );
}

#[test]
fn map_defaults_leave_k_to_the_reads() {
    // map chooses k from the reads only when -k is absent, so a default here would pin it to one value.
    let map_tool = TOOLS.iter().find(|tool| tool.cmd == "map").unwrap();
    let vals: Vec<String> = map_tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();

    assert!(!build_args(map_tool, &vals).iter().any(|arg| arg == "-k"));
}

#[test]
fn flags_and_joined_values_build_single_arguments() {
    let query_tool = TOOLS.iter().find(|tool| tool.cmd == "query").unwrap();
    let mut vals: Vec<String> = query_tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();
    let no_verify_idx = query_tool
        .params
        .iter()
        .position(|param| param.flag == "--no-verify")
        .unwrap();

    assert!(
        !build_args(query_tool, &vals)
            .iter()
            .any(|arg| arg == "--no-verify")
    );
    vals[no_verify_idx] = FLAG_ON.into();
    let args = build_args(query_tool, &vals);
    let flag_pos = args.iter().position(|arg| arg == "--no-verify").unwrap();

    assert!(
        args.get(flag_pos + 1)
            .is_none_or(|arg| arg.starts_with('-')),
        "{args:?}"
    );

    let map_tool = TOOLS.iter().find(|tool| tool.cmd == "map").unwrap();
    let mut vals: Vec<String> = map_tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();
    let cs_idx = map_tool
        .params
        .iter()
        .position(|param| param.flag == "--cs=")
        .unwrap();

    vals[cs_idx] = "long".into();
    assert!(
        build_args(map_tool, &vals)
            .iter()
            .any(|arg| arg == "--cs=long")
    );
}

#[test]
fn every_form_flag_is_accepted_by_its_subcommand() {
    for tool in TOOLS {
        let (vals, bools) = crate::subcommand_flags(tool.cmd).expect(tool.cmd);

        for param in tool.params {
            let accepted = match param.kind {
                ParamKind::Flag => bools.contains(&param.flag),
                _ if param.flag.ends_with('=') => ["short", "long"]
                    .iter()
                    .all(|mode| bools.contains(&format!("{}{mode}", param.flag).as_str())),
                _ => vals.contains(&param.flag),
            };

            assert!(accepted, "{} {}", tool.cmd, param.flag);
        }
    }
}

/// The text of every cell of a frame drawn by `draw` at `width` x `height`.
fn draw_text(width: u16, height: u16, draw: impl FnOnce(&mut Frame)) -> String {
    let mut term =
        ratatui::Terminal::new(ratatui::backend::TestBackend::new(width, height)).unwrap();

    term.draw(draw).unwrap();
    term.backend()
        .buffer()
        .content()
        .iter()
        .map(|cell| cell.symbol())
        .collect()
}

#[test]
fn selection_screen_shows_every_tool() {
    let app = App::new();
    let text = draw_text(160, 60, |frame| render_tool_select(frame, &app));

    for tool in TOOLS {
        assert!(text.contains(tool.name), "{}", tool.name);
    }
}

#[test]
fn every_form_draws_on_a_small_terminal() {
    for (tool_idx, tool) in TOOLS.iter().enumerate() {
        let vals: Vec<String> = tool
            .params
            .iter()
            .map(|param| param.default.to_string())
            .collect();

        for focus in 0..=tool.params.len() {
            draw_text(80, 24, |frame| {
                render_tool_form(frame, tool_idx, &vals, focus, None)
            });
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

    for tool in TOOLS {
        let want = cli
            .iter()
            .find(|(cmd, _)| *cmd == tool.cmd)
            .expect(tool.cmd)
            .1;
        let mut got: Vec<&str> = tool
            .params
            .iter()
            .filter(|param| param.required)
            .map(|param| param.flag)
            .collect();

        got.sort_unstable();
        let mut want = want.to_vec();

        want.sort_unstable();
        assert_eq!(got, want, "{}", tool.cmd);
    }
}

#[test]
fn form_marks_required_fields_and_previews_the_command() {
    let tool_idx = TOOLS.iter().position(|tool| tool.cmd == "query").unwrap();
    let tool = &TOOLS[tool_idx];
    let mut vals: Vec<String> = tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();

    assert!(!missing_required(tool, &vals).is_empty());

    for (param, value) in tool.params.iter().zip(vals.iter_mut()) {
        if param.flag == "-i" {
            *value = "idx.bin".into();
        } else if param.flag == "-q" {
            *value = "reads.fq".into();
        }
    }

    assert!(missing_required(tool, &vals).is_empty());
    let text = draw_text(80, 24, |frame| {
        render_tool_form(frame, tool_idx, &vals, 0, Some("Fill in first: x"))
    });

    assert!(text.contains(" *"), "required marker");
    assert!(text.contains("$ optik query"), "preview");
    assert!(text.contains("Fill in first"), "notice");
    assert!(command_line(tool, &vals).contains("-i idx.bin"));
}

#[test]
fn help_overlay_draws_on_every_screen() {
    let mut app = App::new();

    app.show_help = true;
    let text = draw_text(80, 24, |frame| render(frame, &mut app));

    assert!(text.contains("Help"), "{text}");
    let vals: Vec<String> = TOOLS[0]
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();

    app.screen = Screen::ToolForm {
        tool_idx: 0,
        values: vals,
        focus: 0,
    };
    let text = draw_text(80, 24, |frame| render(frame, &mut app));

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
    let text = draw_text(80, 24, |frame| render(frame, &mut app));

    assert!(text.contains("$ optik query -i idx.bin"), "{text}");
}

/// The rows of a frame drawn by `draw` at `width` x `height`, with trailing spaces removed.
fn draw_rows(width: u16, height: u16, draw: impl FnOnce(&mut Frame)) -> Vec<String> {
    let text = draw_text(width, height, draw);
    let cells: Vec<char> = text.chars().collect();

    cells
        .chunks(width as usize)
        .map(|row| row.iter().collect::<String>().trim_end().to_string())
        .collect()
}

/// The words of `text`, for checking that wrapped text lost nothing.
fn words(text: &str) -> Vec<String> {
    text.split_whitespace().map(str::to_string).collect()
}

#[test]
fn wrap_line_keeps_every_word_within_the_width() {
    let text = "query --ids gives each distinct k-mer one id, but the ids are positions in [0, 12345678), averyveryverylongwordthatmustbreakinside";

    for width in [8usize, 20, 33, 80] {
        let rows = wrap_text(text, width, 2, Style::default());

        for row in &rows {
            assert!(row.width() <= width, "row wider than {width}: {row:?}");
        }

        let joined: String = rows
            .iter()
            .map(|line| {
                line.spans
                    .iter()
                    .map(|span| span.content.to_string())
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
    let selected = tool_by_icon_char('d').map(|tool_idx| TOOLS[tool_idx].cmd);

    assert_eq!(selected, Some("dbruijn"));
    assert_eq!(tool_by_icon_char('g'), None);
}

#[test]
fn long_output_lines_wrap_on_the_results_screen() {
    // The build-spss duplicate warning is longer than an 80-column screen.
    let warning = "[build-spss] warning: 120 of 4554269 indexed k-mer positions repeat a k-mer (4554149 distinct); query --ids gives each distinct k-mer one id, but the ids are positions in [0, 4554269), not a bijection onto [0, 4554149).  Use --compact for a duplicate-free index.";
    let mut app = App::new();

    app.screen = Screen::Results {
        tool_idx: TOOLS
            .iter()
            .position(|tool| tool.cmd == "build-spss")
            .unwrap(),
        output: vec!["build_s\t0.25".into(), warning.into()],
        exit_code: 0,
        scroll: 0,
        formatted: Vec::new(),
        formatted_w: 0,
    };
    let rows = draw_rows(80, 24, |frame| render(frame, &mut app));
    let shown = rows.join(" ").replace('│', " ");

    for word in ["bijection", "duplicate-free", "index."] {
        assert!(shown.contains(word), "{word} missing:\n{}", rows.join("\n"));
    }
}

#[test]
fn form_shows_a_long_field_description_in_full() {
    let tool_idx = TOOLS.iter().position(|tool| tool.cmd == "map").unwrap();
    let tool = &TOOLS[tool_idx];
    let focus = tool
        .params
        .iter()
        .position(|param| param.flag == "-w")
        .unwrap();
    let vals: Vec<String> = tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();
    let rows = draw_rows(80, 30, |frame| {
        render_tool_form(frame, tool_idx, &vals, focus, None)
    });
    let shown = words(&rows.join(" ").replace('│', " ")).join(" ");
    let want = words(tool.params[focus].desc).join(" ");

    assert!(
        shown.contains(&want),
        "description cut:\n{}",
        rows.join("\n")
    );
}

#[test]
fn help_overlay_scrolls_to_its_last_row() {
    let mut app = App::new();
    let tool_idx = TOOLS.iter().position(|tool| tool.cmd == "map").unwrap();
    let vals: Vec<String> = TOOLS[tool_idx]
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();

    app.screen = Screen::ToolForm {
        tool_idx,
        values: vals,
        focus: 0,
    };
    app.show_help = true;
    app.help_scroll = u16::MAX;
    let rows = draw_rows(80, 24, |frame| render(frame, &mut app));

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
    let rows = draw_rows(80, 24, |frame| render(frame, &mut app));
    let shown = rows.join("");

    assert!(shown.contains("ref.fa -q reads.fq"), "{}", rows.join("\n"));
}

#[test]
fn selection_screen_scrolls_to_the_selected_group() {
    let mut app = App::new();
    let groups = grouped_tools();

    for (group_idx, (_, indices)) in groups.iter().enumerate() {
        app.selected_cat = group_idx;
        let text = draw_text(80, 24, |frame| render_tool_select(frame, &app));

        for &tool_idx in indices {
            assert!(
                text.contains(TOOLS[tool_idx].name),
                "group {group_idx}: {} not drawn",
                TOOLS[tool_idx].name
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
    .map(|(key, value)| (key.to_string(), value.to_string()))
    .collect();
    let text: String = format_metrics(&tsv, 60)
        .iter()
        .map(|line| {
            line.spans
                .iter()
                .map(|span| span.content.as_ref())
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
    let query_tool = TOOLS.iter().find(|tool| tool.name == "Query").unwrap();
    let vals: Vec<String> = query_tool
        .params
        .iter()
        .map(|param| param.default.to_string())
        .collect();
    let args = build_args(query_tool, &vals);

    assert!(!args.iter().any(|arg| arg == "-k"), "{args:?}");
}
