//! MBA-1402: `zero --from-angle` solves the zero RANGE(S) a stored bore angle produces.
//!
//! Corrected by the Tier 2 whole-branch review (C2): this is **not** a single-valued inverse
//! of the default `zero --target-distance` mode. A rifle sighted above the bore generally
//! crosses a level line of sight TWICE on a rising shot (once ascending, close to the muzzle
//! — the "near zero"; once descending past the apex — the "far zero"). Both are equally real
//! zero ranges for the SAME bore angle: this is the classic 25/300-yard battle-zero
//! relationship. `zero --target-distance`'s forward solve returns whichever root matches the
//! distance it was asked to hit — the near root for a short/flat zero, the far root for a
//! conventional one — so `--from-angle` cannot know in advance which one is "the" answer and
//! reports both instead.
//!
//! The old version of this test suite asserted the ORIGINAL distance equals the FAR crossing
//! specifically. That happened to hold at every distance/velocity combination the old suite
//! tried (100-500 yd at 2700 fps), but is false in general: at 25 yd (2700 fps) the original
//! distance is the NEAR crossing and the far crossing is ~287 yd; at 50 yd it is ~152 yd; and
//! even an ordinary 100 yd zero at a faster 3400 fps load returns a far crossing ~18% high
//! (its near crossing is the one that matches). The honest invariant, proven below, is that
//! the original distance appears as ONE of the two returned crossings — not that it always
//! matches a specific one of the two.

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

fn solve_angle_imperial(velocity_fps: f64, target_distance_yd: f64) -> f64 {
    let result = run_zero(&[
        "zero",
        "--velocity",
        &velocity_fps.to_string(),
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

/// Both crossings the stored angle implies (display units, imperial). Either may be absent
/// (Tier 2 review C2): `near_zero_range`/`far_zero_range` are additive, skip-when-absent JSON
/// keys, mirroring MBA-1402's own `skip_serializing_if` pattern.
fn solve_crossings_imperial(velocity_fps: f64, angle_deg: f64) -> (Option<f64>, Option<f64>) {
    let result = run_zero(&[
        "zero",
        "--velocity",
        &velocity_fps.to_string(),
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
    (
        result["near_zero_range"].as_f64(),
        result["far_zero_range"].as_f64(),
    )
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

fn solve_crossings_metric(angle_deg: f64) -> (Option<f64>, Option<f64>) {
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
    (
        result["near_zero_range"].as_f64(),
        result["far_zero_range"].as_f64(),
    )
}

/// Round trip is only as precise as the forward solver's own convergence granularity
/// (linear interpolation between adaptive-step sample points feeds both the forward angle
/// search and this inverse's crossing bisection), amplified by near-tangency at short zeros
/// (the descending/ascending slope through the line of sight flattens as the zero distance
/// shrinks, so the same absolute height error maps to a larger downrange error). Observed
/// worst case across the distances/velocities below is a few percent at the shortest,
/// near-tangent zeros; 5% leaves margin while still catching a wrong-root/unit/sign bug
/// outright (which misses by tens to hundreds of percent, not single-digit fractions).
const ROUND_TRIP_TOLERANCE_FRACTION: f64 = 0.05;

/// Assert `original` appears as ONE of the (up to two) returned crossings, within tolerance.
/// This is the honest invariant for a bore angle that generally implies two zeros — NOT that
/// `original` equals a specific one of the two (see module doc comment).
fn assert_original_is_one_of_the_crossings(
    original: f64,
    near: Option<f64>,
    far: Option<f64>,
    label: &str,
) {
    let tolerance = (original * ROUND_TRIP_TOLERANCE_FRACTION).max(0.5);
    let matches_near = near.is_some_and(|n| (n - original).abs() < tolerance);
    let matches_far = far.is_some_and(|f| (f - original).abs() < tolerance);
    assert!(
        matches_near || matches_far,
        "{label}: neither returned crossing recovered the original distance: \
         original={original} near={near:?} far={far:?} tolerance={tolerance}"
    );
}

#[test]
fn round_trip_imperial_several_distances_at_2700_fps() {
    for &distance_yd in &[25.0, 50.0, 100.0, 200.0, 300.0, 500.0] {
        let angle = solve_angle_imperial(2700.0, distance_yd);
        let (near, far) = solve_crossings_imperial(2700.0, angle);
        assert_original_is_one_of_the_crossings(
            distance_yd,
            near,
            far,
            &format!("imperial {distance_yd} yd @ 2700 fps"),
        );
    }
}

/// A second, faster muzzle velocity (Tier 2 review C2): even an ordinary 100 yd zero
/// mis-recovers under the old "always the far crossing" assumption once the load is fast
/// enough that 100 yd is the load's NEAR crossing (its far crossing lands ~18% high).
#[test]
fn round_trip_imperial_several_distances_at_3400_fps() {
    for &distance_yd in &[25.0, 50.0, 100.0, 200.0, 300.0, 500.0] {
        let angle = solve_angle_imperial(3400.0, distance_yd);
        let (near, far) = solve_crossings_imperial(3400.0, angle);
        assert_original_is_one_of_the_crossings(
            distance_yd,
            near,
            far,
            &format!("imperial {distance_yd} yd @ 3400 fps"),
        );
    }
}

#[test]
fn round_trip_metric_several_distances() {
    for &distance_m in &[22.86, 45.72, 91.44, 182.88, 274.32, 457.2] {
        let angle = solve_angle_metric(distance_m);
        let (near, far) = solve_crossings_metric(angle);
        assert_original_is_one_of_the_crossings(distance_m, near, far, &format!("metric {distance_m} m"));
    }
}

/// A short/flat zero (25 yd) has its near crossing at ~25 yd and its far crossing far
/// downrange (~287 yd for this load) -- the classic 25/300-yard battle-zero pairing. Both
/// must be reported, and the near one must be the one close to 25 yd (not confused with the
/// far one, and not silently dropped).
#[test]
fn short_zero_reports_both_the_near_and_a_much_farther_crossing() {
    let angle = solve_angle_imperial(2700.0, 25.0);
    let (near, far) = solve_crossings_imperial(2700.0, angle);

    let near = near.expect("a 25 yd zero's near crossing must be within the solved envelope");
    let far = far.expect("a 25 yd zero's far crossing must be within the solved envelope");

    assert!(
        (near - 25.0).abs() < 2.0,
        "expected the near crossing close to 25 yd, got {near} yd"
    );
    // The far crossing is far downrange -- nowhere close to 25 yd (the classic 25/300
    // relationship; observed ~287 yd for this load).
    assert!(
        far > 150.0,
        "expected a far crossing well beyond 25 yd, got {far} yd"
    );
}

/// `--output table` prints both crossings, clearly labelled.
#[test]
fn table_output_labels_both_crossings() {
    let angle = solve_angle_imperial(2700.0, 25.0);
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args([
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
            &angle.to_string(),
            "--sight-height",
            "2",
            "--output",
            "table",
        ])
        .output()
        .expect("run zero --from-angle table output");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Near Zero"), "{stdout}");
    assert!(stdout.contains("Far Zero"), "{stdout}");
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
