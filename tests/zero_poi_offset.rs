//! MBA-1359: deliberate POI offset at the zero range (Kestrel ZH/ZO semantics).
//!
//! `--zero-poi-up` / `--zero-poi-right` (and the `zero_poi_up_m`/`zero_poi_right_m`
//! solve-json request fields / profile fields) record that the rifle is deliberately zeroed
//! off — e.g. 0.1 in high at 100 yd — and shift the whole solution by the equivalent
//! angular bias (offset / zero distance) applied POST-solve to the zero elevation and
//! azimuth. Pins:
//!   (a) absent/0.0 flags are byte-identical to a run without them;
//!   (b) the vertical bias is exactly offset/zero-distance and scales angularly downrange;
//!   (c) a right offset shifts windage sign-correctly (positive = impacts right);
//!   (d) the bias is a constant additive on the solved angle, identical flat vs inclined;
//!   (e) profile fields round-trip and the CLI flags override them.

use serde_json::Value;
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const BASE_ARGS: &[&str] = &[
    "trajectory",
    "--velocity",
    "2700",
    "--bc",
    "0.475",
    "--mass",
    "168",
    "--diameter",
    "0.308",
    "--max-range",
    "300",
    "--auto-zero",
    "100",
];

fn run_ok(args: &[&str]) -> std::process::Output {
    let output = Command::new(BIN).args(args).output().expect("run binary");
    assert!(
        output.status.success(),
        "command failed: {:?}\nstderr: {}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

/// Linearly interpolate a full-point JSON field at downrange distance `at` (yards).
/// Full-point axes: `z` = downrange, `y` = height, `x` = lateral (all yards).
fn interpolate_at(doc: &Value, field: &str, at: f64) -> f64 {
    let points = doc["trajectory"].as_array().expect("trajectory points");
    let mut prev: Option<(f64, f64)> = None;
    for p in points {
        let z = p["z"].as_f64().expect("z");
        let v = p[field].as_f64().expect("field");
        if let Some((pz, pv)) = prev {
            if pz <= at && z >= at && z > pz {
                return pv + (v - pv) * (at - pz) / (z - pz);
            }
        }
        prev = Some((z, v));
    }
    panic!("no bracketing points around {at} yd for {field}");
}

const OFFSET_IN: f64 = 0.1; // --zero-poi-up test offset, inches
const ZERO_DISTANCE_M: f64 = 91.44; // 100 yd
/// The angular bias the offset must produce: offset / zero distance, radians.
const EXPECTED_BIAS_RAD: f64 = OFFSET_IN * 0.0254 / ZERO_DISTANCE_M;

// (a) explicit 0.0 offsets are byte-identical to a run without the flags.
#[test]
fn zero_valued_flags_are_byte_identical() {
    let mut with_flags: Vec<&str> = BASE_ARGS.to_vec();
    with_flags.extend(["--zero-poi-up", "0", "--zero-poi-right", "0", "-o", "json"]);
    let mut without: Vec<&str> = BASE_ARGS.to_vec();
    without.extend(["-o", "json"]);
    assert_eq!(
        run_ok(&with_flags).stdout,
        run_ok(&without).stdout,
        "0.0-valued zero POI flags must not change a single output byte"
    );
}

// (b) engine level: the solved zero angle shifts by EXACTLY offset/zero-distance.
#[test]
fn vertical_bias_is_exactly_offset_over_zero_distance() {
    let baseline = ballistics_engine::BallisticInputs {
        bc_value: 0.475,
        ..Default::default()
    };
    let mut biased = baseline.clone();
    biased.zero_poi_vertical_m = OFFSET_IN * 0.0254;

    let sight_height = baseline.sight_height;
    let angle_baseline =
        ballistics_engine::calculate_zero_angle(baseline, ZERO_DISTANCE_M, sight_height)
            .expect("baseline zero");
    let angle_biased =
        ballistics_engine::calculate_zero_angle(biased, ZERO_DISTANCE_M, sight_height)
            .expect("biased zero");

    // Post-solve application means the difference is the analytic bias to within the
    // solver's own bisection tolerance reproducibility (the trials are bias-free, so both
    // searches walk the identical bracket sequence — the delta is exact).
    assert!(
        (angle_biased - angle_baseline - EXPECTED_BIAS_RAD).abs() < 1e-12,
        "expected bias {EXPECTED_BIAS_RAD}, got {}",
        angle_biased - angle_baseline
    );
}

// (b) CLI level: +0.1 in up at 100 yd impacts ~0.1 in high at 100 yd and ~2x that at 200 yd
// (angular scaling).
#[test]
fn vertical_offset_shifts_impact_high_and_scales_downrange() {
    let mut baseline_args: Vec<&str> = BASE_ARGS.to_vec();
    baseline_args.extend(["--full", "-o", "json"]);
    let baseline: Value =
        serde_json::from_slice(&run_ok(&baseline_args).stdout).expect("baseline json");

    let mut biased_args: Vec<&str> = BASE_ARGS.to_vec();
    biased_args.extend(["--zero-poi-up", "0.1", "--full", "-o", "json"]);
    let biased: Value = serde_json::from_slice(&run_ok(&biased_args).stdout).expect("json");

    let dy_100 =
        interpolate_at(&biased, "y", 100.0) - interpolate_at(&baseline, "y", 100.0);
    let dy_200 =
        interpolate_at(&biased, "y", 200.0) - interpolate_at(&baseline, "y", 200.0);
    let expected_100_yd = OFFSET_IN / 36.0; // 0.1 in in yards
    assert!(
        (dy_100 - expected_100_yd).abs() < expected_100_yd * 0.05,
        "impact at 100 yd should be ~0.1 in high, got {} yd",
        dy_100
    );
    assert!(
        (dy_200 - 2.0 * expected_100_yd).abs() < expected_100_yd * 0.15,
        "angular bias should double the shift by 200 yd, got {} yd",
        dy_200
    );
}

// (c) a right offset shifts lateral impact right (+x in the full-point frame) by the
// offset at the zero range.
#[test]
fn right_offset_shifts_windage_right_at_zero_range() {
    let mut baseline_args: Vec<&str> = BASE_ARGS.to_vec();
    baseline_args.extend(["--full", "-o", "json"]);
    let baseline: Value =
        serde_json::from_slice(&run_ok(&baseline_args).stdout).expect("baseline json");

    let mut biased_args: Vec<&str> = BASE_ARGS.to_vec();
    biased_args.extend(["--zero-poi-right", "0.5", "--full", "-o", "json"]);
    let biased: Value = serde_json::from_slice(&run_ok(&biased_args).stdout).expect("json");

    let dx_100 =
        interpolate_at(&biased, "x", 100.0) - interpolate_at(&baseline, "x", 100.0);
    let expected_yd = 0.5 / 36.0; // 0.5 in in yards, positive = right
    assert!(
        (dx_100 - expected_yd).abs() < expected_yd * 0.05,
        "impact at 100 yd should be ~0.5 in RIGHT, got {} yd",
        dx_100
    );
}

// (d) the vertical bias is a constant additive on the solved zero angle — identical for a
// flat and an inclined shot (the zero itself always solves level; MBA-1344 convention).
#[test]
fn bias_is_identical_flat_and_inclined() {
    let echo = |extra: &[&str]| -> f64 {
        let mut args: Vec<&str> = BASE_ARGS.to_vec();
        args.extend_from_slice(extra);
        args.extend(["-o", "json"]);
        let doc: Value = serde_json::from_slice(&run_ok(&args).stdout).expect("json");
        doc["zero_angle_degrees"].as_f64().expect("zero echo")
    };

    let flat_delta = echo(&["--zero-poi-up", "0.1"]) - echo(&[]);
    let inclined_delta = echo(&["--zero-poi-up", "0.1", "--shooting-angle", "30"])
        - echo(&["--shooting-angle", "30"]);
    let expected_deg = EXPECTED_BIAS_RAD.to_degrees();

    assert!(
        (flat_delta - expected_deg).abs() < 1e-9,
        "flat bias {flat_delta} deg != expected {expected_deg} deg"
    );
    assert!(
        (flat_delta - inclined_delta).abs() < 1e-9,
        "bias must be identical flat ({flat_delta}) vs inclined ({inclined_delta})"
    );
}

// Wire path (solve-json v1): the new request fields decode (NOT unknown_field), apply, and
// explicit 0.0 stays byte-identical to omitting them.
#[test]
fn solve_json_wire_fields_decode_and_apply() {
    let request_json = |poi: Option<(f64, f64)>| -> String {
        let shot = match poi {
            Some((up, right)) => format!(
                r#"{{"max_range_m": 300.0, "zero_distance_m": 91.44,
                    "zero_poi_up_m": {up}, "zero_poi_right_m": {right}}}"#
            ),
            None => r#"{"max_range_m": 300.0, "zero_distance_m": 91.44}"#.to_string(),
        };
        format!(
            r#"{{
                "schema_version": 1,
                "projectile": {{"mass_kg": 0.01089, "diameter_m": 0.007823,
                                "drag_model": "G1", "ballistic_coefficient": 0.475}},
                "rifle": {{"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05}},
                "shot": {shot},
                "atmosphere": {{}},
                "wind": {{}},
                "solver": {{}},
                "effects": {{}},
                "sampling": {{"interval_m": 25.0}}
            }}"#
        )
    };

    let solve = |json: &str| -> ballistics_engine::solve_json::SolveSuccessV1 {
        let request = ballistics_engine::solve_json::decode_solve_request_v1(json)
            .expect("request with zero POI fields must decode (not unknown_field)");
        ballistics_engine::solve_v1(request).expect("solve")
    };

    let baseline = solve(&request_json(None));
    let biased = solve(&request_json(Some((0.00254, 0.0127)))); // 0.1 in up, 0.5 in right

    let sample_near = |s: &ballistics_engine::solve_json::SolveSuccessV1, at: f64| {
        s.samples
            .iter()
            .min_by(|a, b| {
                (a.distance_m - at)
                    .abs()
                    .total_cmp(&(b.distance_m - at).abs())
            })
            .cloned()
            .expect("sample")
    };
    let b0 = sample_near(&baseline, 91.44);
    let b1 = sample_near(&biased, 91.44);
    // drop_m positive = below LOS: a deliberately-high zero LOWERS drop by ~0.00254 m.
    let drop_delta = b0.drop_m - b1.drop_m;
    assert!(
        (drop_delta - 0.00254).abs() < 0.0005,
        "expected ~2.54 mm high at the zero range, got {} m",
        drop_delta
    );
    // windage_m positive = right.
    let windage_delta = b1.windage_m - b0.windage_m;
    assert!(
        (windage_delta - 0.0127).abs() < 0.0015,
        "expected ~12.7 mm right at the zero range, got {} m",
        windage_delta
    );

    // Explicit 0.0 == omitted, byte-identical response.
    let explicit_zero = solve(&request_json(Some((0.0, 0.0))));
    assert_eq!(
        serde_json::to_string(&baseline).unwrap(),
        serde_json::to_string(&explicit_zero).unwrap(),
        "explicit 0.0 zero POI fields must be byte-identical to omitting them"
    );
}

// (e) profile round-trip: stored SI fields apply, and explicit CLI flags override them.
#[test]
fn profile_fields_round_trip_and_cli_overrides() {
    static N: AtomicU32 = AtomicU32::new(0);
    let home = std::env::temp_dir().join(format!(
        "bx-zero-poi-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&home).unwrap();

    // Save a plain profile, then hand-edit the stored JSON to add the SI offset fields —
    // exactly the documented source for them besides `.a7p --zero-click` import.
    let save = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "profile", "save", "poi-test", "--velocity", "2700", "--bc", "0.475",
            "--drag-model", "g1", "--mass", "168", "--diameter", "0.308",
        ])
        .output()
        .expect("profile save");
    assert!(save.status.success());
    let path = home.join(".ballistics").join("profiles").join("poi-test.json");
    let mut profile: Value =
        serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    profile["zero_poi_up_m"] = Value::from(0.00254); // 0.1 in high
    std::fs::write(&path, serde_json::to_string_pretty(&profile).unwrap()).unwrap();

    let echo = |extra: &[&str]| -> f64 {
        let output = Command::new(BIN)
            .env("HOME", &home)
            .args([
                "trajectory",
                "--saved-profile",
                "poi-test",
                "--max-range",
                "300",
                "--auto-zero",
                "100",
                "-o",
                "json",
            ])
            .args(extra)
            .output()
            .expect("trajectory");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let doc: Value = serde_json::from_slice(&output.stdout).expect("json");
        doc["zero_angle_degrees"].as_f64().expect("zero echo")
    };

    // The stored 0.1-in-high offset must bias the profile solve...
    let profile_delta = echo(&[]) - echo(&["--zero-poi-up", "0"]);
    let expected_deg = EXPECTED_BIAS_RAD.to_degrees();
    assert!(
        (profile_delta - expected_deg).abs() < 1e-9,
        "stored profile offset must bias the zero: got {profile_delta} deg, expected {expected_deg}"
    );
    // ...and an explicit CLI value must override the stored one (0.2 in, not 0.1+0.2).
    let override_delta = echo(&["--zero-poi-up", "0.2"]) - echo(&["--zero-poi-up", "0"]);
    assert!(
        (override_delta - 2.0 * expected_deg).abs() < 1e-9,
        "CLI flag must override (not add to) the stored profile offset: got {override_delta} deg"
    );

    std::fs::remove_dir_all(&home).ok();
}
