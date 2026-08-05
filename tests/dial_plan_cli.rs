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
    assert!(table.contains("travel_exceeded"), "the violation must be named: {table}");

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

// ---- Review fix round 1: I1/I2/I3 (task-6-review.md) ----

/// I1: the ten inline-flag -> `OpticProfile` field mappings (`--travel-up`/`--travel-down`,
/// `--windage-travel-left`/`--windage-travel-right`, `--hold-up`/`--hold-down`/`--hold-left`/
/// `--hold-right`, `--turret-elev`/`--turret-wind`), pinned with TEN PAIRWISE-DISTINCT sentinel
/// magnitudes (8.1..8.8 mil for the eight travel/hold flags, 1.9/2.0 mil -- kept small so
/// `OpticProfile::validate`'s "turret state must be within its own declared travel" check still
/// accepts them -- for the two turret flags) rather than the original suite's symmetric
/// `10mil`-everywhere fixtures, which cannot distinguish a transposed mapping from a correct one.
///
/// Read-back technique: a correction of +-50mil on either axis is far larger than any of the
/// 8.1..8.8 mil bounds, so it unconditionally clamps on `dial_all` (travel) and violates on
/// `hold_all` (hold) regardless of which sentinel ends up where -- so instead of constructing a
/// borderline correction that only clamps under the CORRECT mapping (fragile, and doesn't
/// distinguish "clamped" from "clamped at the wrong number"), each probe reads the reported
/// `available_mil` (the exact declared bound the check consulted) straight out of the JSON and
/// asserts it against the ONE sentinel value that must have produced it. A transposition
/// anywhere in the ten mappings makes some `available_mil` read back a DIFFERENT sentinel than
/// asserted, failing the test. `turret_elev`/`turret_wind` are pinned the same way, via the
/// identity `target_clicks_from_zero - delta_clicks == state_clicks` (true regardless of
/// travel-clamping, so it needs no separate unclamped fixture) recovering the exact state-click
/// count `TurretState` was built from.
///
/// Four probes (elevation up/down, windage right/left) are needed because each declared bound
/// only participates in the check for corrections of its own sign -- one probe alone cannot
/// exercise all four travel bounds or all four hold bounds.
#[test]
fn inline_optic_flags_map_to_their_named_fields_with_distinct_sentinels() {
    // Ten pairwise-distinct sentinels. Expected destinations (verified by hand against
    // src/optic.rs and reconfirmed by this test): travel_up -> elevation TravelLimits::up_mil,
    // travel_down -> down_mil; windage_travel_left -> windage TravelLimits::down_mil (the
    // crossed mapping), windage_travel_right -> up_mil; hold_up/-down/-left/-right -> the
    // identically-named HoldBounds fields (plan_corrections does its own, separately-tested,
    // correction-space/reticle-space crossing on top of these); turret_elev -> TurretState
    // ::elevation_mil, turret_wind -> ::windage_mil.
    const TRAVEL_UP: &str = "8.1mil";
    const TRAVEL_DOWN: &str = "8.2mil";
    const WINDAGE_LEFT: &str = "8.3mil";
    const WINDAGE_RIGHT: &str = "8.4mil";
    const HOLD_UP: &str = "8.5mil";
    const HOLD_DOWN: &str = "8.6mil";
    const HOLD_LEFT: &str = "8.7mil";
    const HOLD_RIGHT: &str = "8.8mil";
    const TURRET_ELEV: &str = "1.9mil"; // must stay inside +-min(travel_up, travel_down)
    const TURRET_WIND: &str = "2.0mil"; // must stay inside +-min(windage_left, windage_right)
    const TURRET_ELEV_CLICKS: i64 = 19; // 1.9mil / 0.1mil
    const TURRET_WIND_CLICKS: i64 = 20; // 2.0mil / 0.1mil

    let base: Vec<&str> = vec![
        "--units", "metric", "dial-plan", "--elevation-click", "0.1mil",
        "--travel-up", TRAVEL_UP, "--travel-down", TRAVEL_DOWN,
        "--windage-travel-left", WINDAGE_LEFT, "--windage-travel-right", WINDAGE_RIGHT,
        "--hold-up", HOLD_UP, "--hold-down", HOLD_DOWN, "--hold-left", HOLD_LEFT,
        "--hold-right", HOLD_RIGHT, "--turret-elev", TURRET_ELEV, "--turret-wind", TURRET_WIND,
        "--range", "100", "-o", "json",
    ];

    // Runs `base` plus one huge (+-50mil) correction on `axis_flag`, and returns
    // (dial_all's violation on `axis_name`, hold_all's violation on `axis_name`, both plans'
    // full instruction pairs) so a single probe can check the travel bound, the hold bound, AND
    // (via both instructions, always present regardless of which axis carries the correction)
    // both turret-state mappings at once.
    let probe = |axis_flag: &str, magnitude: &str, axis_name: &str| -> serde_json::Value {
        let mut args = base.clone();
        args.push(axis_flag);
        args.push(magnitude);
        let (stdout, stderr, ok) = run(&args);
        assert!(ok, "probe {axis_flag}={magnitude}: stderr={stderr}");
        let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");

        let turret_check = |plan: &serde_json::Value| {
            for (idx, expected_state) in
                [(0, TURRET_ELEV_CLICKS), (1, TURRET_WIND_CLICKS)]
            {
                let instr = &plan["instructions"][idx];
                let target = instr["target_clicks_from_zero"].as_i64().unwrap();
                let delta = instr["delta_clicks"].as_i64().unwrap();
                assert_eq!(
                    target - delta,
                    expected_state,
                    "probe {axis_flag}={magnitude}, instruction {idx}: target({target}) - \
                     delta({delta}) should recover the turret_elev/turret_wind sentinel's own \
                     click count ({expected_state}), proving TurretState's field got the right \
                     flag's value: {plan}"
                );
            }
        };

        let dial_all = v["plans"].as_array().unwrap().iter().find(|p| p["strategy"] == "dial_all")
            .expect("dial_all plan present");
        let hold_all = v["plans"].as_array().unwrap().iter().find(|p| p["strategy"] == "hold_all")
            .expect("hold_all plan present");
        turret_check(dial_all);
        turret_check(hold_all);

        let dial_available = dial_all["limits_hit"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["axis"] == axis_name)
            .unwrap_or_else(|| panic!("probe {axis_flag}={magnitude}: no dial_all violation on {axis_name}: {dial_all}"))
            ["available_mil"]
            .clone();
        let hold_available = hold_all["limits_hit"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["axis"] == axis_name)
            .unwrap_or_else(|| panic!("probe {axis_flag}={magnitude}: no hold_all violation on {axis_name}: {hold_all}"))
            ["available_mil"]
            .clone();
        serde_json::json!({"dial_available_mil": dial_available, "hold_available_mil": hold_available})
    };

    let sentinel = |suffixed: &str| -> f64 { suffixed.trim_end_matches("mil").parse().unwrap() };

    // Probe A: elevation UP -- pins --travel-up into TravelLimits::up_mil (not --travel-down's
    // value) via dial_all, and pins --hold-down into the bound a POSITIVE elevation hold
    // consumes (the crossed correction-space/reticle-space rule) via hold_all.
    let a = probe("--elevation", "50mil", "elevation");
    assert_eq!(a["dial_available_mil"].as_f64().unwrap(), sentinel(TRAVEL_UP), "{a}");
    assert_eq!(a["hold_available_mil"].as_f64().unwrap(), sentinel(HOLD_DOWN), "{a}");

    // Probe B: elevation DOWN -- pins --travel-down and --hold-up (the negative-correction
    // bound), the opposite pair from probe A.
    let b = probe("--elevation", "-50mil", "elevation");
    assert_eq!(b["dial_available_mil"].as_f64().unwrap(), sentinel(TRAVEL_DOWN), "{b}");
    assert_eq!(b["hold_available_mil"].as_f64().unwrap(), sentinel(HOLD_UP), "{b}");

    // Probe C: windage RIGHT -- pins --windage-travel-right into TravelLimits::up_mil (the
    // crossed windage travel mapping I1 specifically calls out) and --hold-left.
    let c = probe("--windage", "50mil", "windage");
    assert_eq!(c["dial_available_mil"].as_f64().unwrap(), sentinel(WINDAGE_RIGHT), "{c}");
    assert_eq!(c["hold_available_mil"].as_f64().unwrap(), sentinel(HOLD_LEFT), "{c}");

    // Probe D: windage LEFT -- pins --windage-travel-left into TravelLimits::down_mil and
    // --hold-right, the opposite pair from probe C.
    let d = probe("--windage", "-50mil", "windage");
    assert_eq!(d["dial_available_mil"].as_f64().unwrap(), sentinel(WINDAGE_LEFT), "{d}");
    assert_eq!(d["hold_available_mil"].as_f64().unwrap(), sentinel(HOLD_RIGHT), "{d}");
}

/// I2a: `--prefer-hold` observably reorders `plans[0]`. `hold_all` and `hybrid` are BOTH
/// feasible and tied at EXACTLY `0.0` residual here (no travel declared makes `dial_all`
/// infeasible via `NoTravelData`, taking it out of contention regardless of preference), so the
/// only thing left to break the tie is `Preferences::prefer_hold` -- the default (`false`)
/// ranks `hybrid` first, `--prefer-hold` reverses that to `hold_all`. Fault-injection verified:
/// see the task-6 fix-round report for the substitution that makes this fail.
#[test]
fn prefer_hold_flips_the_top_ranked_strategy_on_a_tie() {
    let base = [
        "--units", "metric", "dial-plan", "--elevation-click", "0.1mil", "--elevation", "2.3mil",
        "--range", "100", "--hold-up", "10mil", "--hold-down", "10mil", "--hold-left", "10mil",
        "--hold-right", "10mil", "-o", "json",
    ];

    let (stdout, stderr, ok) = run(&base);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(v["plans"][0]["strategy"], "hybrid", "default (no --prefer-hold): {v}");

    let mut with_flag = base.to_vec();
    with_flag.push("--prefer-hold");
    let (stdout, stderr, ok) = run(&with_flag);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(v["plans"][0]["strategy"], "hold_all", "--prefer-hold must flip the tie: {v}");
}

/// I2b: `--max-hold` observably tightens the effective hold bound and can flip `feasible`.
/// `--hold-up/-down/-left/-right 10mil` alone comfortably covers the 2.3mil elevation hold
/// `hold_all` needs, so it is feasible with no violations; adding `--max-hold 0.05mil` caps the
/// SAME hold at 0.05mil regardless of the declared 10mil bound, so `hold_all` becomes infeasible
/// with a `hold_bound_exceeded` violation naming the tighter `available_mil: 0.05`, exactly as
/// `plan_corrections`' `Preferences::max_hold_mil` doc promises. Fault-injection verified: see
/// the task-6 fix-round report.
#[test]
fn max_hold_can_force_hold_all_infeasible() {
    let base = [
        "--units", "metric", "dial-plan", "--elevation-click", "0.1mil", "--elevation", "2.3mil",
        "--range", "100", "--travel-up", "30mil", "--travel-down", "5mil", "--hold-up", "10mil",
        "--hold-down", "10mil", "--hold-left", "10mil", "--hold-right", "10mil", "-o", "json",
    ];

    let (stdout, stderr, ok) = run(&base);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let hold_all = v["plans"].as_array().unwrap().iter().find(|p| p["strategy"] == "hold_all")
        .expect("hold_all plan present");
    assert_eq!(hold_all["feasible"], true, "without --max-hold, 10mil covers a 2.3mil hold: {v}");
    assert!(hold_all["limits_hit"].as_array().unwrap().is_empty(), "{v}");

    let mut with_cap = base.to_vec();
    with_cap.extend(["--max-hold", "0.05mil"]);
    let (stdout, stderr, ok) = run(&with_cap);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let hold_all = v["plans"].as_array().unwrap().iter().find(|p| p["strategy"] == "hold_all")
        .expect("hold_all plan present");
    assert_eq!(hold_all["feasible"], false, "--max-hold 0.05mil must cap below the 2.3mil hold: {v}");
    let violation = hold_all["limits_hit"]
        .as_array()
        .unwrap()
        .iter()
        .find(|l| l["axis"] == "elevation")
        .expect("elevation hold_bound_exceeded violation");
    assert_eq!(violation["kind"], "hold_bound_exceeded");
    let available = violation["available_mil"].as_f64().unwrap();
    assert!((available - 0.05).abs() < 1e-9, "{violation}");
}

/// I3: `--range` under `--units imperial` (yards) must convert to metres before reaching
/// `plan_corrections` -- the earlier suite only ever exercised `--units metric`, where that
/// conversion is the identity and a dropped call would still pass. `600` yards is exactly
/// `548.64` metres (`600 * 0.9144`); this test asserts `range_m` in the report equals that
/// converted value (not the raw `600`), AND that `dial_all`'s reported
/// `residual_linear_at_range_m` is consistent with the SAME converted range (recomputed by hand
/// from the response's own angular residuals via the small-angle formula `plan_corrections`
/// itself uses) -- so a dropped, inverted, or wrong-direction conversion is caught either way.
/// Fault-injection verified: see the task-6 fix-round report.
#[test]
fn imperial_range_converts_to_metres_and_scales_the_linear_residual() {
    let (stdout, stderr, ok) = run(&[
        "--units", "imperial", "dial-plan", "--elevation-click", "0.1mil", "--elevation",
        "2.34mil", "--range", "600", "--travel-up", "30mil", "--travel-down", "5mil", "-o",
        "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();

    let range_m = v["range_m"].as_f64().unwrap();
    assert!((range_m - 548.64).abs() < 1e-9, "600 yards must convert to 548.64 m: {range_m}");

    let dial_all = v["plans"].as_array().unwrap().iter().find(|p| p["strategy"] == "dial_all")
        .expect("dial_all plan present");
    let e_res = dial_all["instructions"][0]["residual_mil"].as_f64().unwrap();
    let w_res = dial_all["instructions"][1]["residual_mil"].as_f64().unwrap();
    let expected_linear =
        ((e_res / 1000.0 * range_m).powi(2) + (w_res / 1000.0 * range_m).powi(2)).sqrt();
    let reported_linear = dial_all["residual_linear_at_range_m"].as_f64().unwrap();
    assert!(
        (reported_linear - expected_linear).abs() < 1e-9,
        "reported {reported_linear} vs recomputed-from-converted-range {expected_linear}"
    );
    // A dropped conversion (range_m left as the raw 600) would report ~0.024 here instead --
    // pin the actually-converted number directly too, not just internal self-consistency.
    assert!(
        (reported_linear - 0.0219456).abs() < 1e-6,
        "expected ~0.0219456 m (0.04mil / 1000 * 548.64m): {reported_linear}"
    );
}
