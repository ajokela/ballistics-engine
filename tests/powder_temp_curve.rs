//! Integration tests for `--powder-temp-curve`: a measured powder-temperature ->
//! muzzle-velocity table that supersedes the linear `--powder-temp-sensitivity` model.
//! The muzzle velocity is interpolated from the table at the ambient `--temperature`,
//! clamped at the endpoints (no extrapolation).
//!
//! The first trajectory point's `velocity` is the muzzle velocity, so we read that.

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

const CURVE: &str = "40:2620,70:2700,100:2760";

/// Muzzle velocity (fps) the engine uses at ambient temperature `temp_f`, with the
/// given extra args. Uses -v 2700 as a nominal (the curve overrides it when present).
fn muzzle_velocity(temp_f: &str, extra: &[&str]) -> f64 {
    let mut args: Vec<&str> = vec![
        "trajectory", "-v", "2700", "-b", "0.19", "-m", "77", "-d", "0.224", "--drag-model", "g7",
        "--max-range", "5", "--sight-height", "2.48", "--temperature", temp_f, "--pressure",
        "29.92", "--full", "-o", "json",
    ];
    args.extend_from_slice(extra);
    let out = Command::new(get_cli_binary())
        .args(&args)
        .output()
        .expect("run");
    assert!(
        out.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let json: Value = serde_json::from_str(&String::from_utf8_lossy(&out.stdout)).expect("json");
    json["trajectory"][0]["velocity"].as_f64().expect("v0")
}

fn metric_linear_muzzle_velocity(sensitivity: Option<&str>) -> f64 {
    let mut args = vec![
        "--units",
        "metric",
        "trajectory",
        "-v",
        "800",
        "-b",
        "0.5",
        "-m",
        "10",
        "-d",
        "7.62",
        "--max-range",
        "5",
        "--temperature",
        "40",
        "--pressure",
        "1013.25",
        "--powder-temp",
        "20",
        "--use-powder-sensitivity",
        "--full",
        "-o",
        "json",
    ];
    if let Some(sensitivity) = sensitivity {
        args.extend(["--powder-temp-sensitivity", sensitivity]);
    }
    let out = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run");
    assert!(
        out.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let json: Value = serde_json::from_str(&String::from_utf8_lossy(&out.stdout)).expect("json");
    json["trajectory"][0]["velocity"].as_f64().expect("v0")
}

/// The curve reproduces its measured points exactly.
#[test]
fn curve_hits_measured_points() {
    for (t, expected) in [("40", 2620.0), ("70", 2700.0), ("100", 2760.0)] {
        let v = muzzle_velocity(t, &["--powder-temp-curve", CURVE]);
        assert!(
            (v - expected).abs() < 0.5,
            "at {t} F expected {expected}, got {v}"
        );
    }
}

/// Between points it interpolates linearly (55 F is the midpoint of 40 and 70).
#[test]
fn curve_interpolates_between_points() {
    let v = muzzle_velocity("55", &["--powder-temp-curve", CURVE]);
    assert!((v - 2660.0).abs() < 0.5, "55 F should interpolate to ~2660, got {v}");
}

/// Outside the measured range the velocity is CLAMPED to the endpoint (no extrapolation).
#[test]
fn curve_clamps_outside_range() {
    let cold = muzzle_velocity("10", &["--powder-temp-curve", CURVE]);
    let hot = muzzle_velocity("130", &["--powder-temp-curve", CURVE]);
    assert!((cold - 2620.0).abs() < 0.5, "10 F should clamp to 2620, got {cold}");
    assert!((hot - 2760.0).abs() < 0.5, "130 F should clamp to 2760, got {hot}");
}

/// Without a curve, `-v` is used verbatim regardless of temperature (backward compatible).
#[test]
fn no_curve_uses_velocity_verbatim() {
    let v = muzzle_velocity("100", &[]);
    assert!((v - 2700.0).abs() < 0.5, "no curve should keep -v 2700, got {v}");
}

#[test]
fn explicit_metric_linear_sensitivity_is_distinct_from_omitted_default() {
    let omitted = metric_linear_muzzle_velocity(None);
    let explicit_one = metric_linear_muzzle_velocity(Some("1.0"));

    assert!(
        (omitted - 810.9728).abs() < 1e-6,
        "the omitted metric value must retain the 1 fps/degree-F default: {omitted}"
    );
    assert!(
        (explicit_one - 820.0).abs() < 1e-6,
        "an explicit 1.0 m/s/degree-C sensitivity over 20 C must add 20 m/s: {explicit_one}"
    );

    let imperial_omitted =
        muzzle_velocity("68", &["--use-powder-sensitivity", "--powder-temp", "32"]);
    let imperial_explicit = muzzle_velocity(
        "68",
        &[
            "--use-powder-sensitivity",
            "--powder-temp",
            "32",
            "--powder-temp-sensitivity",
            "1.0",
        ],
    );
    assert!((imperial_omitted - 2736.0).abs() < 1e-6);
    assert_eq!(imperial_explicit.to_bits(), imperial_omitted.to_bits());
}

/// A curve with fewer than two points is rejected.
#[test]
fn curve_requires_two_points() {
    let out = Command::new(get_cli_binary())
        .args([
            "trajectory", "-v", "2700", "-b", "0.19", "-m", "77", "-d", "0.224", "--drag-model",
            "g7", "--max-range", "5", "--temperature", "70", "--powder-temp-curve", "40:2620",
        ])
        .output()
        .expect("run");
    assert!(!out.status.success(), "single-point curve should be rejected");
}

/// With `--powder-temp` given, the curve is interpolated at the POWDER temperature,
/// decoupled from `--temperature` (which drives air density only). MBA follow-up.
#[test]
fn curve_uses_powder_temp_decoupled_from_air() {
    // Hot air (100 F), cold powder (40 F) -> velocity from the 40 F curve point.
    let v = muzzle_velocity("100", &["--powder-temp-curve", CURVE, "--powder-temp", "40"]);
    assert!((v - 2620.0).abs() < 0.5, "cold powder in hot air -> 2620, got {v}");
    // Cold air (40 F), hot powder (100 F) -> 2760.
    let v = muzzle_velocity("40", &["--powder-temp-curve", CURVE, "--powder-temp", "100"]);
    assert!((v - 2760.0).abs() < 0.5, "hot powder in cold air -> 2760, got {v}");
}

/// Without `--powder-temp`, the curve still uses the air temperature (backward compatible).
#[test]
fn curve_defaults_to_air_temp_without_powder_temp() {
    let v = muzzle_velocity("100", &["--powder-temp-curve", CURVE]);
    assert!((v - 2760.0).abs() < 0.5, "no --powder-temp -> curve at air temp 100 = 2760, got {v}");
}
