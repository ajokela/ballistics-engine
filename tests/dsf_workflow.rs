//! End-to-end DSF workflow story (MBA-1357 Task 3): the full two-stage Applied
//! Ballistics-style truing narrative, driven entirely through the real `ballistics`
//! binary — never a synthetic fixture. Six beats:
//!
//! 1. MV calibration precedes DSF: `true-velocity` trues a muzzle velocity from a
//!    supersonic-range observation; the still-supersonic-near-muzzle profile that
//!    carries it has its `dsf` observation rejected with the exact text steering the
//!    shooter back to `true-velocity`.
//! 2. DSF accumulation: several well-separated transonic/subsonic observations each
//!    grow the table.
//! 3. Supersede: an observation within 0.05 Mach of an existing point replaces it.
//! 4. The 6-point cap: a 7th observation is rejected outright, naming `--clear-dsf`.
//! 5. The drop-only invariant, proven by a REAL second solve (not a synthetic
//!    fixture) across both the full-JSON and sampled-CSV trajectory outputs, with
//!    JSON purity (no DSF text, numeric drop fields genuinely differ).
//! 6. `--clear-dsf` restores the untrued output byte-for-byte.
//!
//! This complements, rather than re-derives, `tests/dsf_profile_cli.rs` (individual
//! `dsf`/auto-apply/`profile show` behaviors) and `tests/truing_analysis_cli.rs`
//! (`true-velocity`/`plan-truing` internals) — every observed-drop value below is
//! derived dynamically from a real baseline solve, never a hand-picked constant, so the
//! suite stays correct if the underlying physics ever shifts.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation
//! points `HOME` at a private temp dir (same pattern as the other CLI test files).

use std::path::Path;
use std::process::{Command, Output};

use serde_json::Value;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

/// LOS height (world frame, yd) for every profile in this suite: the default bore
/// height (60 in — no `--bore-height` flag on any `profile save` below) plus the fixed
/// 2.0 in `--sight-height` every profile is given. Same convention documented in
/// `tests/dsf_profile_cli.rs`.
const LOS_YD: f64 = (60.0 + 2.0) * 0.0254 / 0.9144;

fn cli() -> Command {
    Command::new(BIN)
}

fn tempfile_dir() -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "dsf-workflow-test-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn profile_path(home: &Path, name: &str) -> std::path::PathBuf {
    home.join(".ballistics").join("profiles").join(format!("{name}.json"))
}

fn profile_json(home: &Path, name: &str) -> Value {
    serde_json::from_str(&std::fs::read_to_string(profile_path(home, name)).unwrap()).unwrap()
}

fn parse_json(out: &Output) -> Value {
    serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("invalid JSON ({e}): {}", String::from_utf8_lossy(&out.stdout)))
}

/// Save a synthetic profile. Every scenario below derives its own predicted drop
/// dynamically from THIS profile's own real solve — never a hardcoded number — so the
/// exact velocity/BC/zero-distance combination only needs to be internally consistent,
/// not physically realistic.
#[allow(
    clippy::too_many_arguments,
    reason = "flat load parameters mirror `profile save`'s own flat CLI surface"
)]
fn save_profile_inner(
    home: &Path,
    name: &str,
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: &str,
    zero_distance: f64,
    clear_dsf: bool,
) -> Output {
    let mut cmd = cli();
    cmd.env("HOME", home).args([
        "profile",
        "save",
        name,
        "-v",
        &velocity.to_string(),
        "-b",
        &bc.to_string(),
        "-m",
        &mass.to_string(),
        "-d",
        &diameter.to_string(),
        "--drag-model",
        drag_model,
        "--zero-distance",
        &zero_distance.to_string(),
        "--sight-height",
        "2.0",
    ]);
    if clear_dsf {
        cmd.arg("--clear-dsf");
    }
    cmd.output().expect("spawn profile save")
}

#[allow(clippy::too_many_arguments, reason = "see save_profile_inner")]
fn save_profile(
    home: &Path,
    name: &str,
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: &str,
    zero_distance: f64,
) -> Output {
    save_profile_inner(home, name, velocity, bc, mass, diameter, drag_model, zero_distance, false)
}

#[allow(clippy::too_many_arguments, reason = "see save_profile_inner")]
fn save_profile_clear_dsf(
    home: &Path,
    name: &str,
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: &str,
    zero_distance: f64,
) -> Output {
    save_profile_inner(home, name, velocity, bc, mass, diameter, drag_model, zero_distance, true)
}

/// `trajectory --saved-profile NAME --max-range R [--sample-trajectory] -o FMT --full`.
fn full_trajectory_out(home: &Path, name: &str, max_range: f64, sampled: bool, output_fmt: &str) -> Output {
    let mut cmd = cli();
    cmd.env("HOME", home)
        .args(["trajectory", "--saved-profile", name])
        .arg("--max-range")
        .arg(max_range.to_string());
    if sampled {
        cmd.arg("--sample-trajectory");
    }
    cmd.args(["-o", output_fmt, "--full"]);
    cmd.output().expect("spawn trajectory --full")
}

/// `trajectory --saved-profile NAME --max-range R` (default table output, no `--full`) —
/// used to check the "DSF table active (...)" note, which is table-output-only.
fn trajectory_table_out(home: &Path, name: &str, max_range: f64) -> Output {
    cli()
        .env("HOME", home)
        .args(["trajectory", "--saved-profile", name])
        .arg("--max-range")
        .arg(max_range.to_string())
        .output()
        .expect("spawn trajectory table")
}

fn full_trajectory_json(home: &Path, name: &str, max_range: f64) -> Value {
    let out = full_trajectory_out(home, name, max_range, false, "json");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    parse_json(&out)
}

/// `dsf --saved-profile NAME --range R --observed-drop=VALUEin`. The `=`-joined form is
/// required (not a separate `--observed-drop` `VALUEin` pair) so clap accepts a
/// negative `VALUE` as part of one token instead of mistaking it for another flag.
fn run_dsf(home: &Path, name: &str, range_yd: f64, observed_drop_in: f64) -> Output {
    cli()
        .env("HOME", home)
        .args(["dsf", "--saved-profile", name])
        .arg("--range")
        .arg(format!("{range_yd:.2}"))
        .arg(format!("--observed-drop={observed_drop_in:.4}in"))
        .output()
        .expect("spawn dsf")
}

/// Linear interpolation of downrange height `y` (yd) at distance `target_z_yd`,
/// mirroring the engine's own `interpolate_position_and_velocity` bracket-and-lerp
/// approach (src/main.rs) over the SAME `trajectory` JSON array `dsf`'s own internal
/// solve reads from — this is what lets every observed-drop value below be derived from
/// a real solve instead of a hand-picked constant.
fn interp_y_at(points: &[Value], target_z_yd: f64) -> f64 {
    for w in points.windows(2) {
        let z1 = w[0]["z"].as_f64().unwrap();
        let z2 = w[1]["z"].as_f64().unwrap();
        if z2 >= target_z_yd {
            let y1 = w[0]["y"].as_f64().unwrap();
            let y2 = w[1]["y"].as_f64().unwrap();
            if (z2 - z1).abs() < 1e-9 {
                return y2;
            }
            let t = (target_z_yd - z1) / (z2 - z1);
            return y1 + t * (y2 - y1);
        }
    }
    panic!("target range {target_z_yd} yd is beyond the baseline trajectory's solved reach");
}

/// Predicted drop (in), the SAME `LOS - y` convention `dsf`'s own derivation step uses
/// (src/main.rs, `solve_profile_for_dsf`/`interpolate_position_and_velocity`). May be
/// negative when the trajectory arcs above the line of sight (mid-flight, before the
/// profile's own zero range) — that's a real, valid observation, not an error.
fn predicted_drop_in_at(points: &[Value], target_z_yd: f64) -> f64 {
    (LOS_YD - interp_y_at(points, target_z_yd)) * 36.0
}

// ---------------------------------------------------------------------------------
// Scenario 1: MV calibration precedes DSF.
// ---------------------------------------------------------------------------------

/// `true-velocity` (offline, single-observation legacy schema — no `--observed`) trues
/// a muzzle velocity from one supersonic-range drop reading; that trued velocity is
/// then carried into a saved profile (same load), which — still comfortably supersonic
/// at a short range — has its `dsf` observation rejected outright with the exact text
/// steering the shooter back to `true-velocity`, and nothing gets saved.
#[test]
fn stage_one_mv_calibration_precedes_dsf_and_supersonic_gate_rejects() {
    let home = tempfile_dir();

    let tv_out = cli()
        .env("HOME", &home)
        .args([
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "-b",
            "0.27",
            "--drag-model",
            "g7",
            "-m",
            "140",
            "-d",
            "0.264",
            "--offline",
            "-o",
            "json",
        ])
        .output()
        .expect("spawn true-velocity");
    assert!(tv_out.status.success(), "{}", String::from_utf8_lossy(&tv_out.stderr));
    let tv_json = parse_json(&tv_out);
    let trued_velocity = tv_json["effective_velocity"]
        .as_f64()
        .expect("effective_velocity field");
    assert!(
        trued_velocity.is_finite() && (500.0..5000.0).contains(&trued_velocity),
        "trued velocity should be a sane rifle muzzle velocity: {trued_velocity}"
    );
    assert_eq!(tv_json["velocity_unit"], "fps");

    // Carry the trued velocity into a saved profile using the SAME load true-velocity
    // just solved for — still comfortably supersonic near the muzzle.
    let save = save_profile(&home, "mv-calibrated", trued_velocity, 0.27, 140.0, 0.264, "g7", 300.0);
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let dsf_out = run_dsf(&home, "mv-calibrated", 100.0, 0.5);
    assert!(!dsf_out.status.success(), "a short-range observation on this load must still be supersonic");
    let stderr = String::from_utf8_lossy(&dsf_out.stderr);
    assert!(
        stderr.contains("error: observation is supersonic (Mach")
            && stderr.contains(
                "calibrate muzzle velocity first (true-velocity), then collect DSF points at Mach <= 1.2"
            ),
        "{stderr}"
    );

    let saved = profile_json(&home, "mv-calibrated");
    assert!(saved.get("dsf_points").is_none(), "rejected observation must not be saved: {saved}");
}

// ---------------------------------------------------------------------------------
// Scenario 2: DSF accumulation.
// ---------------------------------------------------------------------------------

/// Three observations at Mach values comfortably more than 0.05 apart (50/400/950 yd on
/// a long-zero profile whose muzzle Mach already sits just under the 1.2 ceiling — see
/// module docs above) each Add a new point: the table grows 1 -> 2 -> 3, and every
/// observed drop is derived from the profile's own real baseline solve.
#[test]
fn dsf_accumulates_multiple_observations_growing_the_table() {
    let home = tempfile_dir();
    let save = save_profile(&home, "dsf-accum", 1300.0, 0.4, 168.0, 0.308, "g1", 1500.0);
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let baseline = full_trajectory_json(&home, "dsf-accum", 1600.0);
    let points = baseline["trajectory"].as_array().expect("trajectory array");

    let ranges = [50.0, 400.0, 950.0];
    for (i, &range_yd) in ranges.iter().enumerate() {
        let predicted = predicted_drop_in_at(points, range_yd);
        assert!(predicted.abs() > 1.0, "predicted drop should be a real, non-trivial value: {predicted}");
        let observed = predicted * 1.1;
        let out = run_dsf(&home, "dsf-accum", range_yd, observed);
        assert!(out.status.success(), "range {range_yd}: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(stdout.contains("Added DSF point"), "{stdout}");
        assert!(!stdout.contains("Superseded"), "{stdout}");
        let expected_count = i + 1;
        assert!(stdout.contains(&format!("({expected_count} points, Mach")), "{stdout}");
    }

    let saved = profile_json(&home, "dsf-accum");
    assert_eq!(saved["dsf_points"].as_array().expect("dsf_points").len(), 3);
}

// ---------------------------------------------------------------------------------
// Scenario 3: supersede within Mach tolerance.
// ---------------------------------------------------------------------------------

/// A second observation (90 yd) at Mach ~0.04 from an existing point's Mach (50 yd) —
/// inside the 0.05 supersede tolerance — replaces it outright: the table's point count
/// stays at 1, stdout announces "Superseded", and the OLD dsf ratio is gone from the
/// saved profile (a DIFFERENT ratio than the original makes the replacement provable,
/// not just "same value re-saved").
#[test]
fn dsf_supersedes_a_point_within_mach_tolerance() {
    let home = tempfile_dir();
    let save = save_profile(&home, "dsf-supersede", 1300.0, 0.4, 168.0, 0.308, "g1", 1500.0);
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let baseline = full_trajectory_json(&home, "dsf-supersede", 1600.0);
    let points = baseline["trajectory"].as_array().unwrap();

    let first_predicted = predicted_drop_in_at(points, 50.0);
    let first = run_dsf(&home, "dsf-supersede", 50.0, first_predicted * 1.1);
    assert!(first.status.success(), "{}", String::from_utf8_lossy(&first.stderr));
    assert!(String::from_utf8_lossy(&first.stdout).contains("Added DSF point"));

    let first_dsf_value = profile_json(&home, "dsf-supersede")["dsf_points"][0]["dsf"]
        .as_f64()
        .unwrap();

    let second_predicted = predicted_drop_in_at(points, 90.0);
    let second = run_dsf(&home, "dsf-supersede", 90.0, second_predicted * 1.3);
    assert!(second.status.success(), "{}", String::from_utf8_lossy(&second.stderr));
    let stdout = String::from_utf8_lossy(&second.stdout);
    assert!(stdout.contains("Superseded DSF point"), "{stdout}");
    assert!(stdout.contains("(1 points, Mach"), "table must still hold exactly 1 point: {stdout}");

    let saved_after = profile_json(&home, "dsf-supersede");
    let after_points = saved_after["dsf_points"].as_array().unwrap();
    assert_eq!(after_points.len(), 1, "supersede must not grow the table");
    let second_dsf_value = after_points[0]["dsf"].as_f64().unwrap();
    assert!(
        (second_dsf_value - 1.3).abs() < 0.02,
        "surviving point must carry the NEW ratio: {second_dsf_value}"
    );
    assert!(
        (first_dsf_value - second_dsf_value).abs() > 0.1,
        "old ratio must be gone: first={first_dsf_value} second={second_dsf_value}"
    );
}

// ---------------------------------------------------------------------------------
// Scenario 4: the 6-point cap.
// ---------------------------------------------------------------------------------

/// Six observations, each at a Mach comfortably more than 0.05 from all the others,
/// each Add cleanly (table grows 1..6). A 7th observation — landing in the widest
/// remaining Mach gap, itself still more than 0.05 from every existing point — is
/// rejected outright, mentioning `--clear-dsf`, with a nonzero exit code, and the saved
/// table is left at exactly 6 points, unmodified.
#[test]
fn dsf_table_caps_at_six_points_and_seventh_errors() {
    let home = tempfile_dir();
    let save = save_profile(&home, "dsf-cap", 1300.0, 0.4, 168.0, 0.308, "g1", 1500.0);
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let baseline = full_trajectory_json(&home, "dsf-cap", 1600.0);
    let points = baseline["trajectory"].as_array().unwrap();

    let ranges = [50.0, 200.0, 400.0, 650.0, 950.0, 1300.0];
    for (i, &range_yd) in ranges.iter().enumerate() {
        let predicted = predicted_drop_in_at(points, range_yd);
        let out = run_dsf(&home, "dsf-cap", range_yd, predicted * 1.1);
        assert!(out.status.success(), "range {range_yd}: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(stdout.contains("Added DSF point"), "range {range_yd}: {stdout}");
        let expected_count = i + 1;
        assert!(stdout.contains(&format!("({expected_count} points, Mach")), "{stdout}");
    }

    let saved = profile_json(&home, "dsf-cap");
    assert_eq!(saved["dsf_points"].as_array().unwrap().len(), 6);

    // Distinct range, Mach comfortably inside the widest existing gap (between the
    // first two additions above) — a clean cap rejection, not a supersede.
    let seventh_predicted = predicted_drop_in_at(points, 120.0);
    let seventh = run_dsf(&home, "dsf-cap", 120.0, seventh_predicted * 1.1);
    assert!(!seventh.status.success(), "7th observation must be rejected outright");
    let stderr = String::from_utf8_lossy(&seventh.stderr);
    assert!(
        stderr.contains("already holds the maximum 6 points") && stderr.contains("--clear-dsf"),
        "{stderr}"
    );

    let saved_after = profile_json(&home, "dsf-cap");
    assert_eq!(
        saved_after["dsf_points"].as_array().unwrap().len(),
        6,
        "rejected 7th observation must not modify the saved table"
    );
}

// ---------------------------------------------------------------------------------
// Scenarios 5 + 6: the drop-only invariant, end-to-end through a real second solve —
// and `--clear-dsf` restoring the untrued output exactly.
// ---------------------------------------------------------------------------------

/// Scenario 5: with an active one-point DSF table, `trajectory --saved-profile`'s
/// velocity/energy/time-of-flight/max_range stay byte-identical to the pre-DSF
/// baseline (a REAL second solve through the binary, not `apply_dsf` called on a
/// synthetic fixture — closing Task 1's deferred review Minor), in both the full-JSON
/// `.points` array and the `--sample-trajectory` CSV `.sampled_points` array; JSON
/// output carries no DSF text anywhere. The trajectory's own final (ground-impact)
/// point necessarily sits at a Mach at or below the derivation point's Mach — squarely
/// in `factor_at`'s flat-clamp region (`src/truing_dsf.rs`) — so ITS drop must scale by
/// EXACTLY the derived ratio, not just "some nonzero change".
///
/// Scenario 6: `profile save --clear-dsf` (same load, unchanged) then reproduces the
/// ORIGINAL pre-DSF baseline byte-for-byte.
#[test]
fn drop_only_invariant_end_to_end_then_clear_dsf_restores_exactly() {
    let home = tempfile_dir();
    let save = save_profile(&home, "dsf-invariant", 1300.0, 0.4, 168.0, 0.308, "g1", 300.0);
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let baseline_full = full_trajectory_json(&home, "dsf-invariant", 1500.0);
    let baseline_sampled_out = full_trajectory_out(&home, "dsf-invariant", 1500.0, true, "csv");
    assert!(baseline_sampled_out.status.success(), "{}", String::from_utf8_lossy(&baseline_sampled_out.stderr));
    let baseline_sampled_csv = String::from_utf8_lossy(&baseline_sampled_out.stdout).into_owned();

    // Derive a DSF point from a genuinely subsonic (< 1000 fps), below-LOS point (real
    // positive drop) — same convention as tests/dsf_profile_cli.rs's pick_subsonic_point.
    let (range_yd, predicted_drop_in) = {
        let points = baseline_full["trajectory"].as_array().unwrap();
        points
            .iter()
            .find_map(|p| {
                let v = p["velocity"].as_f64()?;
                if v < 1000.0 {
                    let z = p["z"].as_f64()?;
                    let y = p["y"].as_f64()?;
                    Some((z, (LOS_YD - y) * 36.0))
                } else {
                    None
                }
            })
            .expect("a subsonic point must exist on this profile's trajectory")
    };
    assert!(predicted_drop_in > 0.0, "expected a real positive drop: {predicted_drop_in}");
    let dsf_out = run_dsf(&home, "dsf-invariant", range_yd, predicted_drop_in * 1.3);
    assert!(dsf_out.status.success(), "{}", String::from_utf8_lossy(&dsf_out.stderr));
    assert!(String::from_utf8_lossy(&dsf_out.stdout).contains("Added DSF point"));

    // --- Full JSON: top-level scalars untouched ---
    let corrected_full = full_trajectory_json(&home, "dsf-invariant", 1500.0);
    for key in [
        "units",
        "max_range",
        "max_height",
        "time_of_flight",
        "impact_velocity",
        "impact_energy",
        "stability_coefficient",
        "spin_drift",
        "legend",
    ] {
        assert_eq!(baseline_full[key], corrected_full[key], "field {key} must stay identical");
    }

    let base_points = baseline_full["trajectory"].as_array().unwrap();
    let corr_points = corrected_full["trajectory"].as_array().unwrap();
    assert_eq!(base_points.len(), corr_points.len());
    let mut any_drop_changed = false;
    for (bp, cp) in base_points.iter().zip(corr_points.iter()) {
        assert_eq!(bp["time"], cp["time"], "time must stay identical");
        assert_eq!(bp["velocity"], cp["velocity"], "velocity must stay identical");
        assert_eq!(bp["energy"], cp["energy"], "energy must stay identical");
        assert_eq!(bp["z"], cp["z"], "downrange must stay identical");
        assert_eq!(bp["x"], cp["x"], "lateral position must stay identical");
        let by = bp["y"].as_f64().unwrap();
        let cy = cp["y"].as_f64().unwrap();
        if (by - cy).abs() > 1e-6 {
            any_drop_changed = true;
        }
    }
    assert!(any_drop_changed, "at least one point's drop must change after DSF auto-apply");

    // The trajectory's final point (ground impact) is at or beyond the derivation
    // range, hence at or below the derivation point's Mach — flat-clamped to EXACTLY
    // the derived 1.3 ratio (src/truing_dsf.rs's factor_at). This is the "expected
    // factor direction" check, made exact rather than a vague "it changed".
    let last_base_y = base_points.last().unwrap()["y"].as_f64().unwrap();
    let last_corr_y = corr_points.last().unwrap()["y"].as_f64().unwrap();
    let last_base_drop = LOS_YD - last_base_y;
    let last_corr_drop = LOS_YD - last_corr_y;
    assert!(last_base_drop > 0.0, "the final point must be a real, positive drop below LOS: {last_base_drop}");
    assert!(
        (last_corr_drop / last_base_drop - 1.3).abs() < 0.02,
        "the flat-clamped final point's drop must scale by exactly the derived 1.3 ratio: \
         base_drop={last_base_drop} corr_drop={last_corr_drop} ratio={}",
        last_corr_drop / last_base_drop
    );

    // JSON purity: no DSF text anywhere in the raw stdout.
    let corrected_full_raw =
        String::from_utf8_lossy(&full_trajectory_out(&home, "dsf-invariant", 1500.0, false, "json").stdout)
            .into_owned();
    assert!(!corrected_full_raw.to_lowercase().contains("dsf"), "{corrected_full_raw}");

    // --- Sampled CSV (--sample-trajectory): same invariant, separate array ---
    let corrected_sampled_out = full_trajectory_out(&home, "dsf-invariant", 1500.0, true, "csv");
    assert!(corrected_sampled_out.status.success());
    let corrected_sampled_csv = String::from_utf8_lossy(&corrected_sampled_out.stdout).into_owned();
    let mut base_lines = baseline_sampled_csv.lines();
    let mut corr_lines = corrected_sampled_csv.lines();
    let header = base_lines.next().expect("csv header");
    assert_eq!(corr_lines.next(), Some(header), "csv header must be unchanged");
    assert!(header.starts_with("distance_yd,drop_in,drift_in,velocity_fps"), "{header}");

    let mut any_sampled_drop_changed = false;
    let mut last_base_row: Option<Vec<f64>> = None;
    let mut last_corr_row: Option<Vec<f64>> = None;
    for (b, c) in base_lines.zip(corr_lines) {
        let bcols: Vec<f64> = b.split(',').map(|s| s.parse().unwrap()).collect();
        let ccols: Vec<f64> = c.split(',').map(|s| s.parse().unwrap()).collect();
        assert!((bcols[0] - ccols[0]).abs() < 1e-9, "distance must be byte-identical");
        assert!((bcols[2] - ccols[2]).abs() < 1e-9, "drift must be byte-identical");
        assert!((bcols[3] - ccols[3]).abs() < 1e-9, "velocity must be byte-identical");
        assert!((bcols[4] - ccols[4]).abs() < 1e-9, "energy must be byte-identical");
        assert!((bcols[5] - ccols[5]).abs() < 1e-9, "time must be byte-identical");
        if (bcols[1] - ccols[1]).abs() > 1e-6 {
            any_sampled_drop_changed = true;
        }
        last_base_row = Some(bcols);
        last_corr_row = Some(ccols);
    }
    assert!(any_sampled_drop_changed, "sampled CSV drop must also bend");
    let last_base_drop_in = last_base_row.expect("at least one sampled row")[1];
    let last_corr_drop_in = last_corr_row.expect("at least one sampled row")[1];
    assert!(last_base_drop_in > 0.0, "final sampled row must be a real positive drop: {last_base_drop_in}");
    assert!(
        (last_corr_drop_in / last_base_drop_in - 1.3).abs() < 0.02,
        "final sampled row must also scale by exactly the derived 1.3 ratio: \
         base={last_base_drop_in} corr={last_corr_drop_in}"
    );
    assert!(!corrected_sampled_csv.to_lowercase().contains("dsf"));

    // Table output: the note IS present (table-output-only).
    let table_out = trajectory_table_out(&home, "dsf-invariant", 1500.0);
    assert!(table_out.status.success());
    assert!(String::from_utf8_lossy(&table_out.stdout).contains("DSF table active (1 points"));

    // --- Scenario 6: --clear-dsf restores the untrued output EXACTLY ---
    let clear = save_profile_clear_dsf(&home, "dsf-invariant", 1300.0, 0.4, 168.0, 0.308, "g1", 300.0);
    assert!(clear.status.success(), "{}", String::from_utf8_lossy(&clear.stderr));
    let cleared_profile = profile_json(&home, "dsf-invariant");
    assert!(cleared_profile.get("dsf_points").is_none(), "{cleared_profile}");

    let restored_full = full_trajectory_json(&home, "dsf-invariant", 1500.0);
    assert_eq!(
        restored_full, baseline_full,
        "full JSON trajectory must be byte-for-byte identical to the pre-DSF baseline"
    );

    let restored_table_out = trajectory_table_out(&home, "dsf-invariant", 1500.0);
    assert!(restored_table_out.status.success());
    assert!(!String::from_utf8_lossy(&restored_table_out.stdout).contains("DSF table active"));
}
