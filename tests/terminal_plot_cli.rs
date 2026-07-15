//! Integration tests for `trajectory --plot` (MBA-1320): the inline terminal chart.
//!
//! Conventions this locks down (see CLI_USAGE.md "Terminal Chart (`--plot`)"):
//! - The flag is additive: omitted, every byte of `-o table` output before where the
//!   chart would go is unchanged from a pre-MBA-1320 run.
//! - Bare `--plot` resolves to the braille renderer (clap `default_missing_value`);
//!   `--plot ascii` is the explicit ASCII fallback. Both live only under `-o table` —
//!   JSON/CSV stay machine-readable and don't grow the chart text.
//! - The canvas-primitive/render_chart unit tests live in src/terminal_plot.rs itself;
//!   this file exercises the real, wired-up CLI end to end, including one byte-exact
//!   "full deterministic solve" chart golden per the ticket.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

/// A short, deterministic, flat-fire shot: short --max-range keeps the solve (and thus
/// this file's golden) small, per the ticket's "keep goldens small".
const BASE_ARGS: &[&str] = &[
    "trajectory", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5", "--max-range", "100",
];

const CHART_MARKER: &str = "Drop vs Range:";

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().expect("run ballistics")
}

fn run_ok(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "command failed ({:?}): {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).to_string()
}

/// The exact two-chart block (drop vs. range, then lateral drift vs. range) for
/// BASE_ARGS + `--plot`, captured once from a real run and pinned here as a golden.
/// Any change to render_chart's formatting, the trajectory solve, or the point data fed
/// into the charts must be a deliberate re-generation of this constant, not a silent
/// diff.
const EXPECTED_CHART: &str = r#"Drop vs Range:
┌ drop (yd) — y:[1.60, 1.67] ────────────────────────────────────────────┐
│⠉⠀⠁⠀⠁⠀⠀⠐⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⢀│
└ x:[0.00, 100.00] ──────────────────────────────────────────────────────┘

Lateral Drift vs Range:
┌ drift (yd) — y:[-0.50, 0.50] ──────────────────────────────────────────┐
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠉⠀⠁⠀⠁⠀⠀⠈⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠈⠀⠈│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
└ x:[0.00, 100.00] ──────────────────────────────────────────────────────┘
"#;

#[test]
fn full_deterministic_solve_chart_is_byte_exact() {
    let mut args = BASE_ARGS.to_vec();
    args.push("--plot");
    let stdout = run_ok(&args);
    let idx = stdout
        .find(CHART_MARKER)
        .expect("--plot must print a \"Drop vs Range:\" section");
    let chart_section = &stdout[idx..];
    assert_eq!(chart_section, EXPECTED_CHART);
}

#[test]
fn bare_plot_defaults_to_braille() {
    let mut args = BASE_ARGS.to_vec();
    args.push("--plot");
    let stdout = run_ok(&args);
    let idx = stdout.find(CHART_MARKER).expect("chart present");
    let chart_section = &stdout[idx..];
    // A braille canvas necessarily contains non-ASCII glyphs (U+2800..=U+28FF);
    // the ASCII fallback never does.
    assert!(
        chart_section.chars().any(|c| c as u32 >= 0x2800 && c as u32 <= 0x28FF),
        "bare --plot should render braille dots:\n{chart_section}"
    );
}

#[test]
fn plot_ascii_uses_only_ascii_inside_the_canvas() {
    // The FRAME (┌─┐│└─┘ and the axis text's em dash) is plain Unicode box-drawing
    // regardless of style — those glyphs have near-universal font coverage, unlike the
    // braille block, so there's no reason to ASCII-ify them too. Only the interior DOT
    // CANVAS rows switch between braille dots and '*'/' '; scope this assertion to just
    // those rows (between the top and bottom border lines) rather than the whole block.
    let mut args = BASE_ARGS.to_vec();
    args.extend(["--plot", "ascii"]);
    let stdout = run_ok(&args);
    let idx = stdout.find(CHART_MARKER).expect("chart present");
    let chart_section = &stdout[idx..];
    let canvas_lines: Vec<&str> = chart_section
        .lines()
        .filter(|line| line.starts_with('\u{2502}')) // '│'
        .collect();
    assert!(!canvas_lines.is_empty(), "expected at least one canvas row");
    for line in &canvas_lines {
        // Strip the leading/trailing '│' frame characters themselves (non-ASCII) before
        // checking the interior.
        let interior = line.trim_matches('\u{2502}');
        assert!(
            interior.is_ascii(),
            "--plot ascii canvas row must be pure ASCII: {line:?}"
        );
    }
    // And it should actually have plotted something (not an all-blank canvas).
    assert!(chart_section.contains('*'));
}

#[test]
fn omitting_plot_leaves_output_unchanged_up_to_where_the_chart_would_start() {
    let without_plot = run_ok(BASE_ARGS);
    assert!(
        !without_plot.contains(CHART_MARKER),
        "no --plot flag must never print a chart"
    );

    let mut with_plot_args = BASE_ARGS.to_vec();
    with_plot_args.push("--plot");
    let with_plot = run_ok(&with_plot_args);
    let idx = with_plot.find(CHART_MARKER).expect("chart present");
    // Everything BEFORE the chart marker (summary box, notes, impact/ground messages)
    // must be unchanged whether or not --plot was passed — additive, no reformatting of
    // the pre-existing output. `with_plot[..idx]` includes one extra blank line (the
    // `println!()` spacer main.rs emits right before "Drop vs Range:"), which
    // `without_plot` naturally doesn't have anything after to produce — trim_end() on
    // both sides normalizes that away so this compares the actual content, not
    // incidental trailing blank-line bookkeeping.
    assert_eq!(without_plot.trim_end(), with_plot[..idx].trim_end());
}

#[test]
fn plot_does_not_affect_json_output() {
    let mut args = BASE_ARGS.to_vec();
    args.extend(["-o", "json", "--plot"]);
    let stdout = run_ok(&args);
    assert!(
        !stdout.contains(CHART_MARKER),
        "--plot must not touch JSON output (only -o table renders a chart)"
    );
    let _: serde_json::Value = serde_json::from_str(&stdout).expect("still valid JSON");
}

#[test]
fn plot_does_not_affect_csv_output() {
    let mut args = BASE_ARGS.to_vec();
    args.extend(["-o", "csv", "--plot", "--full"]);
    let stdout = run_ok(&args);
    assert!(
        !stdout.contains(CHART_MARKER),
        "--plot must not touch CSV output (only -o table renders a chart)"
    );
}

#[test]
fn invalid_plot_value_is_rejected_at_parse_time() {
    let mut args = BASE_ARGS.to_vec();
    args.extend(["--plot", "bogus"]);
    let out = run(&args);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("invalid value"), "stderr was: {stderr}");
}
