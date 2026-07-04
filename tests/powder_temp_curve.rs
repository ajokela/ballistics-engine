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
