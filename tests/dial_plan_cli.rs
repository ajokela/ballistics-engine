//! MBA-1348 Plan B Task 6: the `dial-plan` CLI surface over `plan_corrections` (`src/optic.rs`).
//!
//! `dial-plan` runs no solve at all -- the TRUE angular correction is a direct CLI input, not
//! something derived from a trajectory -- so every fixture here is pure, instant arithmetic
//! and the whole file runs in a fraction of a second.
//!
//! Two ways to supply the optic are exercised: inline flags (`--elevation-click` and the rest
//! of the turret/hold flags, mirroring `profile save`'s own set) for most tests, and
//! `--profile NAME` (via a real saved profile under an isolated `$HOME`, matching
//! `tests/profile_turret_carry_forward.rs`'s hermetic pattern) for the profile-specific tests.
//!
//! A recurring theme below: 0.1 mil (and 1/4 MOA) click graduations are not exactly
//! representable in binary floating point, so an "exact click alignment" case like `2.3mil` on
//! `0.1mil` clicks does not produce a bit-for-bit `0.0` residual, only one within a few
//! `f64::EPSILON` of it. Assertions on residual/hold magnitudes use a small tolerance for
//! exactly this reason instead of `==`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(bin()).args(args).output().expect("run ballistics");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

/// A fresh, uniquely-named scratch `$HOME` per call -- safe under `cargo test`'s default
/// parallel execution (mirrors `tests/profile_turret_carry_forward.rs`'s own helper).
fn tempfile_dir(tag: &str) -> PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "dial-plan-cli-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn save_profile(home: &Path, name: &str, extra: &[&str]) {
    let mut args: Vec<&str> =
        vec!["profile", "save", name, "-v", "2700", "-b", "0.475", "-m", "175", "-d", "0.308"];
    args.extend_from_slice(extra);
    let output = Command::new(bin()).env("HOME", home).args(&args).output().expect("profile save");
    assert!(
        output.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn profile_path(home: &Path, name: &str) -> PathBuf {
    home.join(".ballistics").join("profiles").join(format!("{name}.json"))
}

fn run_with_home(home: &Path, args: &[&str]) -> (String, String, bool) {
    let output = Command::new(bin()).env("HOME", home).args(args).output().expect("run ballistics");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

// ---- (a) exact-click alignment: `plans[0]` is `dial_all`, residual is (within float noise) 0,
// and the table names the dial instruction ----

/// A correction that is a whole number of 0.1 mil clicks (23 of them), with elevation travel
/// declared generously enough that `dial_all` can affirm its own feasibility promise. No hold
/// bounds are declared, so `hold_all`/`hybrid` (whose own promise needs a verified hold) rank
/// behind it regardless of their residual -- feasibility is the PRIMARY ranking key
/// (`plan_corrections`' own doc comment) -- which is exactly the property this test pins.
#[test]
fn exact_click_alignment_ranks_dial_all_first_with_near_zero_residual() {
    let base = [
        "--units", "metric", "dial-plan", "--elevation", "2.3mil", "--range", "600",
        "--elevation-click", "0.1mil", "--travel-up", "30mil", "--travel-down", "5mil",
    ];

    let mut json_args = base.to_vec();
    json_args.extend(["-o", "json"]);
    let (stdout, stderr, ok) = run(&json_args);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["plans"][0]["strategy"], "dial_all");
    assert_eq!(v["plans"][0]["feasible"], true);
    let residual = v["plans"][0]["residual_linear_at_range_m"].as_f64().unwrap();
    assert!(residual.abs() < 1e-6, "expected ~0 residual, got {residual}");
    let elevation_instr = &v["plans"][0]["instructions"][0];
    assert_eq!(elevation_instr["axis"], "elevation");
    assert_eq!(elevation_instr["direction"], "up");
    assert_eq!(elevation_instr["target_clicks_from_zero"], 23);

    let (table, stderr, ok) = run(&base);
    assert!(ok, "stderr: {stderr}");
    assert!(
        table.contains("dial UP 23 clicks"),
        "table must render the dial instruction: {table}"
    );
}

// ---- (b) a declared `clicks_per_revolution` shows "rev N + M" phrasing ----

#[test]
fn clicks_per_revolution_annotation_shows_rev_phrasing() {
    let (table, stderr, ok) = run(&[
        "--units", "metric", "dial-plan", "--elevation", "2.3mil", "--range", "600",
        "--elevation-click", "0.1mil", "--travel-up", "30mil", "--travel-down", "5mil",
        "--clicks-per-rev", "10",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(
        table.contains("rev 2 + 3"),
        "23 clicks at 10 clicks/rev must render as 2 revolutions + 3 clicks: {table}"
    );
}

// ---- (c) infeasible travel is a successful (exit 0) analysis that names the violation ----

/// Elevation needs 5mil but only 1mil of travel is declared in each direction -- `dial_all`
/// cannot affirm its whole-correction promise and is reported INFEASIBLE, but the command
/// itself still exits 0 (concluding "this doesn't fit" is not a usage error). Generous hold
/// bounds are declared so `hold_all`/`hybrid` are feasible, isolating the one violation this
/// test is about instead of drowning it in unrelated `NoHoldBoundData` noise.
#[test]
fn infeasible_travel_still_exits_zero_and_names_the_violation() {
    let base = [
        "--units", "metric", "dial-plan", "--elevation", "5mil", "--range", "100",
        "--elevation-click", "0.1mil", "--travel-up", "1mil", "--travel-down", "1mil",
        "--hold-up", "10mil", "--hold-down", "10mil", "--hold-left", "10mil", "--hold-right",
        "10mil",
    ];

    let (table, stderr, ok) = run(&base);
    assert!(ok, "an infeasibility analysis must still exit 0: stderr={stderr}");
    assert!(table.contains("INFEASIBLE"), "{table}");
    assert!(table.contains("TravelExceeded"), "the violation must be named: {table}");

    let mut json_args = base.to_vec();
    json_args.extend(["-o", "json"]);
    let (stdout, stderr, ok) = run(&json_args);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    let dial_all = v["plans"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["strategy"] == "dial_all")
        .expect("a dial_all plan is always present");
    assert_eq!(dial_all["feasible"], false);
    let kinds: Vec<&str> =
        dial_all["limits_hit"].as_array().unwrap().iter().map(|l| l["kind"].as_str().unwrap()).collect();
    assert!(kinds.contains(&"travel_exceeded"), "{kinds:?}");
}

// ---- (d) `-o json` is `serde_json::to_string_pretty(&DialPlanReportV1)`, verbatim ----

/// Independently computes the SAME report through the library `dial-plan` itself calls
/// (identical inputs -- same deterministic floating-point arithmetic, so this is a byte-exact
/// comparison, not an approximate one) and asserts the CLI's stdout equals its pretty-printed
/// serialization exactly, with nothing else printed around it.
#[test]
fn json_output_is_the_verbatim_pretty_printed_report() {
    use ballistics_engine::adjustment::{ClickBase, ClickValue};
    use ballistics_engine::optic::{
        plan_corrections, AngularCorrection, OpticProfile, Preferences, TravelLimits,
    };

    let optic = OpticProfile {
        elevation_click: ClickValue { size: 0.1, base: ClickBase::Mil },
        windage_click: ClickValue { size: 0.1, base: ClickBase::Mil },
        clicks_per_revolution: None,
        zero_stop: false,
        elevation_travel: Some(TravelLimits { up_mil: 30.0, down_mil: 5.0 }),
        windage_travel: None,
        turret_state: None,
        reticle_hold_bounds: None,
    };
    let correction = AngularCorrection { elevation_mil: 2.3, windage_mil: 0.0 };
    let expected = plan_corrections(correction, &optic, 600.0, 1.0, 1.0, &Preferences::default())
        .expect("valid inputs");
    let expected_json = serde_json::to_string_pretty(&expected).unwrap();

    let (stdout, stderr, ok) = run(&[
        "--units", "metric", "dial-plan", "--elevation", "2.3mil", "--range", "600",
        "--elevation-click", "0.1mil", "--travel-up", "30mil", "--travel-down", "5mil", "-o",
        "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert_eq!(
        stdout.trim_end(),
        expected_json,
        "the CLI must print serde_json::to_string_pretty(&report) verbatim, nothing wrapped"
    );

    // The specific fields the task calls out, asserted directly (redundant with the full
    // equality above, but pins the intent even if that stricter check is ever loosened).
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["method"], "dial_space_quantization_v1");
    assert_eq!(v["assumptions"].as_array().unwrap().len(), 5);
}

// ---- (e) missing both --profile and --elevation-click names both flags ----

#[test]
fn missing_both_profile_and_elevation_click_names_both_flags() {
    let (_, stderr, ok) = run(&["dial-plan", "--elevation", "1mil", "--range", "100"]);
    assert!(!ok, "an optic-less dial-plan must be rejected");
    assert!(stderr.contains("--profile"), "{stderr}");
    assert!(stderr.contains("--elevation-click"), "{stderr}");
}

// ---- (f) the three strategies render distinctly, and windage falls back to the elevation
// click graduation when --windage-click is omitted ----

/// A correction that is NOT aligned to the click graduation on either axis, with travel and
/// hold bounds generous enough that all three strategies are feasible -- isolating the
/// per-axis INSTRUCTION LINE text itself (direction, delta clicks, hold) as the thing that
/// must differ, rather than merely differing because one plan is infeasible. `dial_all` and
/// `hybrid` dial the identical click count on each axis (same underlying `AxisComputation`,
/// see `src/optic.rs`) -- it is specifically the hold-mil portion of the line that must set
/// them apart, since `dial_all` always holds 0 and `hybrid` holds the leftover here.
///
/// `--windage-click` is deliberately never given: windage (1.16mil / 0.1mil = 11.6, rounding
/// to 12 clicks) only comes out right if it fell back to the elevation graduation, matching
/// `profile save`'s own precedent (`ProfileData::optic_profile`).
#[test]
fn three_strategies_render_distinct_instruction_lines() {
    let (table, stderr, ok) = run(&[
        "--units", "metric", "dial-plan", "--elevation", "2.34mil", "--windage", "1.16mil",
        "--range", "100", "--elevation-click", "0.1mil", "--travel-up", "30mil", "--travel-down",
        "5mil", "--windage-travel-left", "10mil", "--windage-travel-right", "10mil", "--hold-up",
        "10mil", "--hold-down", "10mil", "--hold-left", "10mil", "--hold-right", "10mil",
    ]);
    assert!(ok, "stderr: {stderr}");

    let elevation_lines: Vec<&str> =
        table.lines().filter(|l| l.trim_start().starts_with("elevation: dial")).collect();
    assert_eq!(elevation_lines.len(), 3, "one elevation instruction line per plan: {table}");
    let mut unique = elevation_lines.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(
        unique.len(),
        3,
        "all three strategies must render distinct elevation instruction lines: {elevation_lines:?}"
    );
    // Windage fell back to the 0.1mil elevation click: 1.16 / 0.1 = 11.6, nearest click 12.
    assert!(
        table.contains("windage: dial RIGHT 12 clicks"),
        "windage must fall back to --elevation-click's graduation: {table}"
    );
}

// ---- Beyond the brief's list: properties the KEY FACTS section calls out by name ----

#[test]
fn at_least_one_correction_axis_is_required() {
    let (_, stderr, ok) =
        run(&["dial-plan", "--range", "100", "--elevation-click", "0.1mil"]);
    assert!(!ok, "neither --elevation nor --windage given must be rejected");
    assert!(stderr.contains("--elevation"), "{stderr}");
    assert!(stderr.contains("--windage"), "{stderr}");
}

#[test]
fn csv_and_pdf_are_rejected() {
    for format in ["csv", "pdf"] {
        let (_, stderr, ok) = run(&[
            "dial-plan", "--elevation", "1mil", "--range", "100", "--elevation-click", "0.1mil",
            "-o", format,
        ]);
        assert!(!ok, "dial-plan -o {format} was accepted");
        assert!(stderr.contains("no "), "{format}: {stderr}");
    }
}

#[test]
fn profile_and_inline_optic_flags_conflict() {
    let home = tempfile_dir("conflict");
    save_profile(&home, "rifle", &["--elevation-click", "0.1mil"]);

    let (_, stderr, ok) = run_with_home(
        &home,
        &[
            "dial-plan", "--profile", "rifle", "--elevation", "1mil", "--range", "100",
            "--elevation-click", "0.2mil",
        ],
    );
    assert!(!ok, "combining --profile with an inline optic flag must be rejected");
    assert!(stderr.contains("--profile"), "{stderr}");
}

#[test]
fn profile_without_optic_data_is_a_named_error() {
    let home = tempfile_dir("no-optic");
    save_profile(&home, "bare", &[]);

    let (_, stderr, ok) = run_with_home(
        &home,
        &["dial-plan", "--profile", "bare", "--elevation", "1mil", "--range", "100"],
    );
    assert!(!ok, "a profile with no saved click graduation must be rejected");
    assert!(stderr.contains("bare"), "{stderr}");
    assert!(stderr.contains("--elevation-click"), "{stderr}");
}

/// The module's own hand-derived worked example (`cf_dial_space_worked_example`,
/// `src/optic.rs`), reached through `--profile` instead of a direct library call: CF 0.98 on
/// 0.1 mil clicks needs a dial target of `5.0 / 0.98 = 5.10204...` mil, quantizing to 51
/// clicks, which execute `51 * 0.1 * 0.98 = 4.998` true mil -- proving the CLI actually reads
/// `elevation_cf` off the saved profile and feeds it to `plan_corrections`, not a silently
/// assumed 1.0. `profile save` has no `--elevation-cf` flag (it is derived by `tall-target`
/// and hand-edited in, per that field's own doc comment), so the CF is patched into the saved
/// JSON directly, exactly as a shooter following that workflow would.
#[test]
fn profile_supplies_optic_and_tracking_cf_matches_the_worked_example() {
    let home = tempfile_dir("cf-worked-example");
    save_profile(&home, "rifle", &["--elevation-click", "0.1mil"]);

    let path = profile_path(&home, "rifle");
    let mut stored: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    stored["elevation_cf"] = serde_json::json!(0.98);
    std::fs::write(&path, stored.to_string()).unwrap();

    let (stdout, stderr, ok) = run_with_home(
        &home,
        &[
            "--units", "metric", "dial-plan", "--profile", "rifle", "--elevation", "5mil",
            "--range", "100", "-o", "json",
        ],
    );
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let hybrid = v["plans"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["strategy"] == "hybrid")
        .expect("a hybrid plan is always present");
    let elevation = &hybrid["instructions"][0];
    assert_eq!(elevation["target_clicks_from_zero"], 51);
    let dial_mil_true = elevation["dial_mil_true"].as_f64().unwrap();
    assert!((dial_mil_true - 4.998).abs() < 1e-9, "{dial_mil_true}");
    let hold_mil = elevation["hold_mil"].as_f64().unwrap();
    assert!((hold_mil - 0.002).abs() < 1e-9, "{hold_mil}");
    let residual_mil = elevation["residual_mil"].as_f64().unwrap();
    assert!(residual_mil.abs() < 1e-9, "hybrid must reconstruct exactly: {residual_mil}");
}
