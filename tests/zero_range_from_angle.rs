//! MBA-1402: `zero --from-angle` solves the zero RANGE a stored bore angle produces — the
//! exact inverse of the default `zero --target-distance` mode.
//!
//! Two properties are load-bearing here:
//!
//! 1. **Round trip.** Solving a distance to an angle, then feeding that angle back through
//!    the inverse, must return (approximately) the original distance. This is the strongest
//!    available correctness check for an inverse solver — a wrong root, a unit-conversion
//!    slip, or a sign error would all show up as a large round-trip miss.
//! 2. **The FAR crossing.** A flat-fire trajectory crosses a level line of sight twice (once
//!    ascending near the muzzle, once descending past the apex); "zero range" always means
//!    the far one. A test that only checked "some crossing was found" would pass even if the
//!    solver returned the near crossing and silently mis-zeroed every caller.

use serde_json::Value;
use std::process::Command;

fn run_zero(args: &[&str]) -> Value {
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(args)
        .output()
        .expect("run zero command");
    assert!(
        output.status.success(),
        "zero command failed (args={args:?}): {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "zero stdout is not JSON (args={args:?}): {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

fn solve_angle_imperial(target_distance_yd: f64) -> f64 {
    let result = run_zero(&[
        "zero",
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--target-distance",
        &target_distance_yd.to_string(),
        "--sight-height",
        "2",
        "--output",
        "json",
    ]);
    result["zero_angle_degrees"].as_f64().unwrap()
}

fn solve_range_from_angle_imperial(angle_deg: f64) -> f64 {
    let result = run_zero(&[
        "zero",
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--from-angle",
        &angle_deg.to_string(),
        "--sight-height",
        "2",
        "--output",
        "json",
    ]);
    result["zero_range"].as_f64().unwrap()
}

fn solve_angle_metric(target_distance_m: f64) -> f64 {
    let result = run_zero(&[
        "--units",
        "metric",
        "zero",
        "--velocity",
        "822.96",
        "--bc",
        "0.475",
        "--mass",
        "10.88621688",
        "--diameter",
        "7.8232",
        "--target-distance",
        &target_distance_m.to_string(),
        "--sight-height",
        "50.8",
        "--output",
        "json",
    ]);
    result["zero_angle_degrees"].as_f64().unwrap()
}

fn solve_range_from_angle_metric(angle_deg: f64) -> f64 {
    let result = run_zero(&[
        "--units",
        "metric",
        "zero",
        "--velocity",
        "822.96",
        "--bc",
        "0.475",
        "--mass",
        "10.88621688",
        "--diameter",
        "7.8232",
        "--from-angle",
        &angle_deg.to_string(),
        "--sight-height",
        "50.8",
        "--output",
        "json",
    ]);
    result["zero_range"].as_f64().unwrap()
}

/// Round trip is only as precise as the forward solver's own convergence granularity
/// (linear interpolation between adaptive-step sample points feeds both the forward
/// angle search and this inverse's far-crossing bisection). Observed worst case across the
/// distances below is ~0.35% (at the shortest, near-tangent 100 yd zero); 2% leaves ample
/// margin while still catching a wrong-root/unit/sign bug outright (which misses by tens of
/// percent, not fractions of one).
const ROUND_TRIP_TOLERANCE_FRACTION: f64 = 0.02;

fn assert_round_trips(original: f64, recovered: f64, label: &str) {
    let tolerance = (original * ROUND_TRIP_TOLERANCE_FRACTION).max(0.5);
    assert!(
        (recovered - original).abs() < tolerance,
        "{label}: round trip did not recover the original distance: original={original} \
         recovered={recovered} tolerance={tolerance}"
    );
}

#[test]
fn round_trip_imperial_several_distances() {
    for &distance_yd in &[100.0, 150.0, 200.0, 300.0, 500.0] {
        let angle = solve_angle_imperial(distance_yd);
        let recovered = solve_range_from_angle_imperial(angle);
        assert_round_trips(distance_yd, recovered, &format!("imperial {distance_yd} yd"));
    }
}

#[test]
fn round_trip_metric_several_distances() {
    for &distance_m in &[91.44, 137.16, 182.88, 274.32, 457.2] {
        let angle = solve_angle_metric(distance_m);
        let recovered = solve_range_from_angle_metric(angle);
        assert_round_trips(distance_m, recovered, &format!("metric {distance_m} m"));
    }
}

/// A slower bullet (more arc) zeroed at 100 yd puts its near crossing at ~11 yd and its far
/// crossing at ~100 yd -- a large, unambiguous separation. If the inverse solver picked the
/// FIRST sign change instead of the LAST, this would return ~11 yd, not ~100 yd.
#[test]
fn inverse_returns_the_far_crossing_not_the_near_one() {
    let angle = {
        let result = run_zero(&[
            "zero",
            "--velocity",
            "1000",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--target-distance",
            "100",
            "--sight-height",
            "2",
            "--output",
            "json",
        ]);
        result["zero_angle_degrees"].as_f64().unwrap()
    };

    let result = run_zero(&[
        "zero",
        "--velocity",
        "1000",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--from-angle",
        &angle.to_string(),
        "--sight-height",
        "2",
        "--output",
        "json",
    ]);
    let recovered = result["zero_range"].as_f64().unwrap();

    // The far crossing is close to the original 100 yd target.
    assert!(
        (recovered - 100.0).abs() < 2.0,
        "expected the far crossing near 100 yd, got {recovered} yd"
    );
    // The near crossing sits around 11 yd for this load -- nowhere close to the far one.
    // A solver that returned the FIRST crossing instead of the LAST would land here.
    assert!(
        recovered > 50.0,
        "inverse appears to have returned the NEAR crossing ({recovered} yd), not the far one"
    );
}

/// `--target-distance` and `--from-angle` are mutually exclusive and jointly required: clap
/// must reject both omitted and both supplied, before any solve runs.
#[test]
fn target_distance_and_from_angle_are_mutually_exclusive_and_jointly_required() {
    let common = [
        "zero",
        "--velocity",
        "2700",
        "--bc",
        "0.475",
        "--mass",
        "168",
        "--diameter",
        "0.308",
        "--output",
        "json",
    ];

    // Neither given.
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(common)
        .output()
        .expect("run zero with neither flag");
    assert!(
        !output.status.success(),
        "expected failure with neither --target-distance nor --from-angle"
    );

    // Both given.
    let mut both: Vec<&str> = common.to_vec();
    both.extend(["--target-distance", "100", "--from-angle", "0.1"]);
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(&both)
        .output()
        .expect("run zero with both flags");
    assert!(
        !output.status.success(),
        "expected failure with both --target-distance and --from-angle"
    );
}
