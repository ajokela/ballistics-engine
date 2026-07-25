//! Integration tests for MBA-1397: `--pressure-type <absolute|qnh>` (and the paired
//! `--zero-pressure-type` auto-zero override), the sea-level-corrected altimeter-setting
//! (QNH) pressure mode.
//!
//! Hard requirement: with the toggle absent or explicitly `absolute`, every surface must be
//! byte-identical to pre-MBA-1397 output. The `*_byte_identical_to_omitted` tests below
//! compare raw stdout bytes (not just parsed field subsets) to prove that.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

fn get_cli_binary() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("ballistics");
    if !path.exists() {
        path.pop();
        path.pop();
        path.push("release");
        path.push("ballistics");
    }
    path
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("failed to execute ballistics binary")
}

fn run_json(args: &[&str]) -> Value {
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

const BASE_TRAJECTORY_ARGS: &[&str] = &[
    "trajectory",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--units",
    "metric",
    "--max-range",
    "300",
    "--ignore-ground-impact",
    "-o",
    "json",
    "--altitude",
    "1500",
    "--pressure",
    "1030.0",
];

const BASE_ZERO_ARGS: &[&str] = &[
    "zero",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--units",
    "metric",
    "--target-distance",
    "300",
    "-o",
    "json",
    "--altitude",
    "1500",
    "--pressure",
    "1030.0",
];

fn with_extra<'a>(base: &[&'a str], extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = base.to_vec();
    args.extend_from_slice(extra);
    args
}

// ---- Hard requirement: byte-identical when absent/absolute ---------------------------------

#[test]
fn trajectory_omitted_pressure_type_is_byte_identical_to_explicit_absolute() {
    let omitted = run(&with_extra(BASE_TRAJECTORY_ARGS, &["--full"]));
    let explicit = run(&with_extra(
        BASE_TRAJECTORY_ARGS,
        &["--full", "--pressure-type", "absolute"],
    ));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --pressure-type must be byte-identical to --pressure-type absolute"
    );
}

#[test]
fn zero_omitted_pressure_type_is_byte_identical_to_explicit_absolute() {
    let omitted = run(BASE_ZERO_ARGS);
    let explicit = run(&with_extra(BASE_ZERO_ARGS, &["--pressure-type", "absolute"]));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --pressure-type must be byte-identical to --pressure-type absolute on zero"
    );
}

#[test]
fn trajectory_auto_zero_omitted_zero_pressure_type_is_byte_identical_to_explicit_absolute() {
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "1500",
        "--pressure",
        "1013.25",
        "--auto-zero",
        "300",
        "--zero-pressure",
        "1030.0",
    ];
    let omitted = run(args);
    let explicit = run(&with_extra(args, &["--zero-pressure-type", "absolute"]));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --zero-pressure-type must be byte-identical to --zero-pressure-type absolute"
    );
}

// ---- QNH must actually change the physics (reduced, not raw, pressure reaches the solve) ---

/// Hand-computed pin (see also src/atmosphere.rs's
/// `qnh_reduction_matches_hand_computed_value_at_nonzero_altitude`): reducing a QNH of
/// 1030.0 hPa at 1500 m gives 859.575312392644... hPa station pressure -- strictly, and
/// substantially, lower than the raw 1030.0 hPa. Lower pressure means less dense air, less
/// drag, and a HIGHER retained impact velocity at a fixed range.
#[test]
fn trajectory_qnh_mode_yields_a_higher_impact_velocity_than_absolute_at_the_same_reading() {
    let absolute = run_json(BASE_TRAJECTORY_ARGS);
    let qnh = run_json(&with_extra(BASE_TRAJECTORY_ARGS, &["--pressure-type", "qnh"]));

    let v_absolute = absolute["impact_velocity"].as_f64().unwrap();
    let v_qnh = qnh["impact_velocity"].as_f64().unwrap();
    assert!(
        v_qnh > v_absolute + 10.0,
        "QNH-reduced (lower) pressure must retain velocity better than treating the same \
         reading as absolute station pressure: absolute={v_absolute} qnh={v_qnh}"
    );
}

#[test]
fn zero_qnh_mode_yields_a_smaller_zero_angle_than_absolute_at_the_same_reading() {
    let absolute = run_json(BASE_ZERO_ARGS);
    let qnh = run_json(&with_extra(BASE_ZERO_ARGS, &["--pressure-type", "qnh"]));

    let angle_absolute = absolute["zero_angle_degrees"].as_f64().unwrap();
    let angle_qnh = qnh["zero_angle_degrees"].as_f64().unwrap();
    // Less dense air (QNH correctly reduced) means less drop over the same zero distance,
    // so a SMALLER compensating zero angle is needed.
    assert!(
        angle_qnh < angle_absolute,
        "QNH-reduced pressure must need a smaller zero angle: absolute={angle_absolute} \
         qnh={angle_qnh}"
    );
    assert!((angle_absolute - angle_qnh).abs() > 1e-4);
}

/// `--zero-pressure-type` is independent of the shot-day `--pressure-type`: this test uses a
/// shot day at the plain ICAO default (so the flown trajectory is unaffected) and only
/// changes how the auto-zero SOLVE interprets `--zero-pressure`. The two runs must therefore
/// differ (proving the zero-day reduction actually ran) even though the difference is a
/// second-order effect on the muzzle angle rather than a first-order density change.
#[test]
fn zero_pressure_type_qnh_changes_the_auto_zero_result_vs_absolute() {
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "1500",
        "--pressure",
        "1013.25",
        "--auto-zero",
        "300",
        "--zero-pressure",
        "1030.0",
    ];
    let absolute = run_json(&with_extra(args, &["--zero-pressure-type", "absolute"]));
    let qnh = run_json(&with_extra(args, &["--zero-pressure-type", "qnh"]));
    let v_absolute = absolute["impact_velocity"].as_f64().unwrap();
    let v_qnh = qnh["impact_velocity"].as_f64().unwrap();
    assert_ne!(
        v_absolute, v_qnh,
        "zero-day QNH reduction must change the solved zero angle, and therefore the flown \
         trajectory, versus treating the same zero-day reading as absolute"
    );
}

#[test]
fn qnh_pressure_at_sea_level_matches_absolute_exactly() {
    // At sea level (altitude 0) the QNH reduction is the identity: QNH IS station pressure.
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "0",
        "--pressure",
        "1005.0",
    ];
    let absolute = run(args);
    let qnh = run(&with_extra(args, &["--pressure-type", "qnh"]));
    assert!(absolute.status.success());
    assert!(qnh.status.success());
    assert_eq!(
        absolute.stdout, qnh.stdout,
        "at sea level, QNH and absolute must be byte-identical (QNH reduces to itself)"
    );
}
