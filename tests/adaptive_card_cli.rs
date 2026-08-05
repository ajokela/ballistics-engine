//! MBA-1351 (0.33.0 decision-support Plan B Task 12): the `adaptive-card` CLI surface over
//! Task 11's `ballistics_engine::card::adaptive_card` engine.
//!
//! Every fixture here is the same representative metric .308-class load `src/card.rs`'s own
//! `adaptive_card_tests` and `src/hold_curve.rs`'s own tests use (`-v 800 -b 0.223 -m 10.9
//! -d 7.82 --drag-model g7`), run under `--units metric` so `--start`/`--end`/`--anchor`
//! line up directly with the engine's native metre grid -- this file's numbers were
//! independently confirmed against the real binary before being pinned (see Task 12's
//! report for the exact commands run).
//!
//! Domains are kept small (200-500 m or a narrower slice of it), per the brief's own
//! fixture-size guidance, so the whole file runs quickly even though (unlike `dial-plan`)
//! every test here solves a real trajectory.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

/// The shared load fixture, as CLI args, with `$1` left as a placeholder for the
/// subcommand-specific tail (domain/budget/output flags).
const LOAD_ARGS: &[&str] = &[
    "-v", "800", "-b", "0.223", "-m", "10.9", "-d", "7.82", "--drag-model", "g7",
    "--sight-height", "45", "--zero-distance", "100", "--temperature", "15",
    "--pressure", "1013.25", "--humidity", "50", "--wind-speed", "3", "--wind-direction", "90",
];

fn run(tail: &[&str]) -> (String, String, bool) {
    let mut args: Vec<&str> = vec!["--units", "metric", "adaptive-card"];
    args.extend_from_slice(LOAD_ARGS);
    args.extend_from_slice(tail);
    let output = Command::new(bin()).args(&args).output().expect("run ballistics");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

/// A fresh, uniquely-named scratch `$HOME` per call -- safe under `cargo test`'s default
/// parallel execution (mirrors `tests/dial_plan_cli.rs`'s own helper).
fn tempfile_dir(tag: &str) -> PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "adaptive-card-cli-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn save_profile(home: &Path, name: &str, extra: &[&str]) {
    let mut args: Vec<&str> = vec!["--units", "metric", "profile", "save", name];
    args.extend_from_slice(LOAD_ARGS);
    args.extend_from_slice(extra);
    let output = Command::new(bin()).env("HOME", home).args(&args).output().expect("profile save");
    assert!(output.status.success(), "profile save failed: {}", String::from_utf8_lossy(&output.stderr));
}

fn run_with_home(home: &Path, tail: &[&str]) -> (String, String, bool) {
    let mut args: Vec<&str> = vec!["--units", "metric", "adaptive-card", "--profile"];
    let output = {
        args.push("fixture");
        args.extend_from_slice(tail);
        Command::new(bin()).env("HOME", home).args(&args).output().expect("run ballistics")
    };
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

// ---- (a) happy path table: footer states budget met + worst error; row count <= max ----

#[test]
fn happy_path_table_states_budget_met_and_worst_error_within_row_cap() {
    let (table, stderr, ok) = run(&["--start", "200", "--end", "500"]);
    assert!(ok, "stderr: {stderr}");
    assert!(table.contains("budget met: yes"), "{table}");
    assert!(table.contains("worst error: elevation"), "{table}");
    assert!(table.contains("verification grid: 0.9144 m step"), "{table}");
    assert!(table.contains("rows: "), "{table}");

    // Row count <= --max-rows (default 25), read back from the footer's own "rows: N of 25
    // max" line rather than re-deriving it, so this assertion fails if that line's shape
    // ever changes out from under it.
    let rows_line = table.lines().find(|l| l.starts_with("rows: ")).expect("rows line");
    let n: usize = rows_line
        .trim_start_matches("rows: ")
        .split_whitespace()
        .next()
        .expect("row count")
        .parse()
        .expect("row count must be a number");
    assert!(n <= 25, "{rows_line}");
    // The domain is small and smooth (per the brief's own fixture-size guidance): the fixed
    // seed is 2 rows (both ends), so a handful more from real refinement is expected, but
    // nowhere near the 25-row cap.
    assert!(n <= 12, "expected a small card, got {n} rows: {table}");
}

// ---- (b) -o json carries schema_version/method/5 assumptions + budget_met ----

#[test]
fn json_output_carries_schema_method_assumptions_and_budget_met() {
    let (stdout, stderr, ok) = run(&["--start", "200", "--end", "500", "-o", "json"]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["method"], "greedy_worst_point_insertion_on_holdcurve_grid_v1");
    assert_eq!(v["assumptions"].as_array().expect("assumptions array").len(), 5);
    assert_eq!(v["budget_met"], true);
    assert_eq!(v["rows_capped"], false);
    assert!(v["worst_elevation_error"].as_f64().unwrap() >= 0.0);
    assert!(v["worst_windage_error"].as_f64().unwrap() >= 0.0);
    assert_eq!(v["verification_grid_step_m"].as_f64().unwrap(), 0.9144);

    // JSON is the report VERBATIM: range stays in the engine's native metres even under
    // `--units metric` (== same numeric space here, but the point is it is NOT re-derived
    // from `--units` the way the table/csv columns are) -- both domain ends are rows.
    let rows = v["rows"].as_array().expect("rows array");
    assert!(rows.len() >= 2);
    assert_eq!(rows.first().unwrap()["range"].as_f64().unwrap(), 200.0);
    assert_eq!(rows.last().unwrap()["range"].as_f64().unwrap(), 500.0);
    // Fields this engine never populates still serialize as explicit JSON null, not as an
    // absent key.
    assert!(rows[0]["come_up"].is_null());
    assert!(rows[0]["lead_adj"].is_null());
}

// ---- (c) budget below the click floor: exit 0, footer says budget NOT met with the
// measured floor (an honest report is a successful run) ----

#[test]
fn budget_below_click_floor_reports_not_met_with_the_measured_floor() {
    let home = tempfile_dir("floor");
    save_profile(&home, "fixture", &["--elevation-click", "0.1mil"]);

    // A small domain and a generous --max-rows so the search stops at the irreducible
    // half-click floor, not the row cap -- isolating exactly the property under test,
    // the same separation Task 11's own `quantization_floor_is_honest_and_terminates`
    // test makes.
    let (stdout, stderr, ok) = run_with_home(
        &home,
        &[
            "--start", "300", "--end", "320", "--max-rows", "500",
            "--elevation-budget", "0.001mil", "--windage-budget", "0.001mil",
            "-o", "json",
        ],
    );
    assert!(ok, "an honest budget-not-met report is still exit 0; stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["budget_met"], false);
    assert_eq!(v["rows_capped"], false, "the row cap must not be what stopped this search");
    let worst = v["worst_elevation_error"].as_f64().unwrap();
    // Half of 0.1mil, measured -- close to but never over the floor Task 11's engine pins.
    assert!((0.04..=0.0500001).contains(&worst), "worst_elevation_error {worst} not at the half-click floor");

    let (table, stderr, ok) = run_with_home(
        &home,
        &[
            "--start", "300", "--end", "320", "--max-rows", "500",
            "--elevation-budget", "0.001mil", "--windage-budget", "0.001mil",
        ],
    );
    assert!(ok, "stderr: {stderr}");
    assert!(table.contains("budget met: no"), "{table}");
    assert!(table.contains("worst error: elevation 0.05"), "footer must show the measured floor: {table}");

    let _ = std::fs::remove_dir_all(&home);
}

// ---- (d) --anchor outside --start/--end -> usage error naming the flag ----

#[test]
fn anchor_outside_domain_is_a_usage_error_naming_the_flag() {
    let (_stdout, stderr, ok) =
        run(&["--start", "200", "--end", "500", "--anchor", "100"]);
    assert!(!ok, "an out-of-domain anchor must be a usage error, not a silent drop");
    assert!(stderr.contains("--anchor"), "{stderr}");
    assert!(stderr.contains("--start"), "{stderr}");
    assert!(stderr.contains("--end"), "{stderr}");
}

// ---- --elevation-budget/--windage-budget without a unit suffix -> a clap usage error
// naming the flag, the same as every other angular flag on dial-plan/profile save (final
// branch review N-3: these two used to fall back to a bare "needs a unit suffix" message
// with nothing saying which flag) ----

#[test]
fn elevation_budget_without_a_unit_suffix_is_a_usage_error_naming_the_flag() {
    let (_stdout, stderr, ok) =
        run(&["--start", "200", "--end", "500", "--elevation-budget", "0.1"]);
    assert!(!ok, "a budget with no unit suffix must be a usage error");
    assert!(stderr.contains("--elevation-budget"), "{stderr}");
}

#[test]
fn windage_budget_without_a_unit_suffix_is_a_usage_error_naming_the_flag() {
    let (_stdout, stderr, ok) =
        run(&["--start", "200", "--end", "500", "--windage-budget", "0.1"]);
    assert!(!ok, "a budget with no unit suffix must be a usage error");
    assert!(stderr.contains("--windage-budget"), "{stderr}");
}

// ---- (e) -o pdf writes a non-empty file ----

#[cfg(feature = "pdf")]
#[test]
fn pdf_output_writes_a_non_empty_file() {
    let out_path = std::env::temp_dir().join(format!(
        "adaptive_card_mba1351_{}_{}.pdf",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let (_stdout, stderr, ok) = run(&[
        "--start", "200", "--end", "500", "-o", "pdf",
        "--output-file", out_path.to_str().unwrap(),
    ]);
    assert!(ok, "stderr: {stderr}");
    let bytes = std::fs::read(&out_path).expect("PDF written");
    assert!(!bytes.is_empty(), "PDF output must not be empty");
    assert_eq!(&bytes[..5], b"%PDF-", "output must start with a PDF header");
    let _ = std::fs::remove_file(&out_path);
}

#[cfg(feature = "pdf")]
#[test]
fn pdf_output_requires_output_file() {
    let (_stdout, stderr, ok) = run(&["--start", "200", "--end", "500", "-o", "pdf"]);
    assert!(!ok);
    assert!(stderr.contains("--output-file"), "{stderr}");
}

// ---- come_up_row_line trap (Task 9 review M1): adaptive rows always carry `come_up:
// None`, which panics `come_up_row_line` on ANY row (its unwrap runs before the i==0
// em-dash branch). adaptive-card's table renderer deliberately reuses `range_table_row_line`
// instead (adaptive rows are shaped exactly like range-table's: come_up/lead_adj always
// None, everything else Some), which never reads `.come_up` at all -- this test pins that
// the FIRST row (and the whole table) renders without panicking, as a regression guard in
// case a future refactor ever routes adaptive rows through the come-ups renderer instead. ----

#[test]
fn first_row_renders_without_panicking() {
    let (table, stderr, ok) = run(&["--start", "200", "--end", "220"]);
    assert!(ok, "stderr: {stderr}");
    assert!(!stderr.contains("panicked"), "stderr: {stderr}");
    // The first data row (range 200) must actually be present, not merely "didn't crash".
    let first_data_line = table
        .lines()
        .find(|l| l.trim_start().starts_with('│') && l.contains("200"))
        .expect("first row must render");
    assert!(first_data_line.contains('│'), "{first_data_line}");
}

// ---- supporting coverage: CSV is well-formed and the footer goes to stderr, not stdout ----

#[test]
fn csv_output_is_clean_and_the_footer_goes_to_stderr() {
    let (stdout, stderr, ok) = run(&["--start", "200", "--end", "500", "-o", "csv"]);
    assert!(ok, "stderr: {stderr}");
    let mut lines = stdout.lines();
    assert_eq!(
        lines.next().unwrap(),
        "range_m,drop_mm,drop_mil,wind_mm,wind_mil,velocity_m/s,energy_J,tof_s"
    );
    for line in lines {
        assert_eq!(line.split(',').count(), 8, "{line}");
    }
    // The footer's prose must not leak into the CSV data channel...
    assert!(!stdout.contains("budget met"), "stdout: {stdout}");
    // ...but the user still sees it, on stderr.
    assert!(stderr.contains("budget met: yes"), "stderr: {stderr}");
}

// ---- supporting coverage: default budget == half the profile's click, exactly ----

#[test]
fn default_budget_matches_half_the_profiles_click_exactly() {
    let home = tempfile_dir("default-budget");
    save_profile(&home, "fixture", &["--elevation-click", "0.1mil", "--windage-click", "0.2mil"]);

    let (default_json, stderr, ok) =
        run_with_home(&home, &["--start", "200", "--end", "500", "-o", "json"]);
    assert!(ok, "stderr: {stderr}");
    let (explicit_json, stderr, ok) = run_with_home(
        &home,
        &[
            "--start", "200", "--end", "500",
            "--elevation-budget", "0.05mil", "--windage-budget", "0.1mil",
            "-o", "json",
        ],
    );
    assert!(ok, "stderr: {stderr}");

    let default_v: serde_json::Value = serde_json::from_str(&default_json).unwrap();
    let explicit_v: serde_json::Value = serde_json::from_str(&explicit_json).unwrap();
    assert_eq!(default_v["rows"].as_array().unwrap().len(), explicit_v["rows"].as_array().unwrap().len());
    assert_eq!(
        default_v["worst_elevation_error"].as_f64().unwrap(),
        explicit_v["worst_elevation_error"].as_f64().unwrap()
    );
    // The elevation assertion above cannot discriminate the WINDAGE default on its own -- bind
    // that axis too, so a swapped or mis-derived windage default (e.g. reusing elevation's
    // click instead of windage's) fails this test instead of shipping silently.
    assert_eq!(
        default_v["worst_windage_error"].as_f64().unwrap(),
        explicit_v["worst_windage_error"].as_f64().unwrap()
    );

    let _ = std::fs::remove_dir_all(&home);
}

// ---- supporting coverage: --max-rows caps and the footer admits it ----

#[test]
fn max_rows_caps_and_the_footer_admits_it() {
    let (table, stderr, ok) = run(&[
        "--start", "200", "--end", "500", "--max-rows", "2",
        "--elevation-budget", "0.001mil", "--windage-budget", "0.001mil",
    ]);
    assert!(ok, "an honest capped report is still exit 0; stderr: {stderr}");
    assert!(table.contains("budget met: no"), "{table}");
    assert!(table.contains("rows: 2 of 2 max (row cap reached)"), "{table}");
}

// ---- review fix I-1: a fractional row range must print the SAME value on every format,
// not get silently rounded to a whole display unit on table/csv while json/pdf keep the
// exact value. Reproduces the review's own case: `--anchor 412.5` on a 200-500 m metric card.

#[test]
fn fractional_anchor_range_agrees_across_json_table_and_csv() {
    let tail = ["--start", "200", "--end", "500", "--anchor", "412.5"];

    let (json_out, stderr, ok) = {
        let mut t = tail.to_vec();
        t.extend(["-o", "json"]);
        run(&t)
    };
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&json_out).expect("json");
    let ranges: Vec<f64> = v["rows"].as_array().unwrap().iter().map(|r| r["range"].as_f64().unwrap()).collect();
    assert!(ranges.contains(&412.5), "json must carry the exact anchor: {ranges:?}");

    let (table, stderr, ok) = run(&tail);
    assert!(ok, "stderr: {stderr}");
    assert!(
        table.contains("│ 412.5 │"),
        "table must print the exact anchor 412.5, not a rounded '412': {table}"
    );
    assert!(
        !table.contains("│   412 │"),
        "table must not silently round the 412.5 m anchor to a whole metre: {table}"
    );

    let (csv, stderr, ok) = {
        let mut t = tail.to_vec();
        t.extend(["-o", "csv"]);
        run(&t)
    };
    assert!(ok, "stderr: {stderr}");
    assert!(
        csv.lines().any(|l| l.starts_with("412.5,")),
        "csv must print the exact anchor 412.5, not a rounded '412': {csv}"
    );
    assert!(
        !csv.lines().any(|l| l.starts_with("412,")),
        "csv must not silently round the 412.5 m anchor to a whole metre: {csv}"
    );
}

/// The near-integer rule's OTHER branch: a range that lands within floating-point noise of a
/// whole display unit (both domain ends here) must still print as a bare integer, not grow a
/// spurious ".0" -- the I-1 fix must not regress the common, previously-correct case.
#[test]
fn whole_unit_ranges_still_print_as_bare_integers() {
    let (table, stderr, ok) = run(&["--start", "200", "--end", "500"]);
    assert!(ok, "stderr: {stderr}");
    // A range cell formatted as "200.0" would render "│ 200.0 │" (6-wide, one decimal), not
    // "│   200 │" (6-wide, bare integer) -- the positive match below is mutually exclusive
    // with the regression this test guards against, so no separate negative check is needed.
    assert!(table.contains("│   200 │"), "domain start must print as a bare integer: {table}");
    assert!(table.contains("│   500 │"), "domain end must print as a bare integer: {table}");

    let (csv, stderr, ok) = run(&["--start", "200", "--end", "500", "-o", "csv"]);
    assert!(ok, "stderr: {stderr}");
    assert!(csv.lines().any(|l| l.starts_with("200,")), "{csv}");
    assert!(csv.lines().any(|l| l.starts_with("500,")), "{csv}");
}
