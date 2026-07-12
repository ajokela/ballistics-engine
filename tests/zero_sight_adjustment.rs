//! Reported sight adjustment must use the angle between the bore and the target line of sight.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const TARGET_DISTANCE_YARDS: f64 = 100.0;
const YARDS_TO_METERS: f64 = 0.9144;

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

fn run_zero(target_height_yards: f64) -> Value {
    let output = Command::new(get_cli_binary())
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
            "--target-distance",
            "100",
            "--target-height",
            &target_height_yards.to_string(),
            "--sight-height",
            "2",
            "--temperature",
            "59",
            "--pressure",
            "29.92",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute zero command");
    assert!(
        output.status.success(),
        "zero command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("valid zero JSON")
}

fn run_zero_metric() -> Value {
    let output = Command::new(get_cli_binary())
        .args([
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
            "91.44",
            "--target-height",
            "0",
            "--sight-height",
            "50.8",
            "--temperature",
            "15",
            "--pressure",
            "1013.207888",
            "--output",
            "json",
        ])
        .output()
        .expect("run metric zero command");
    assert!(
        output.status.success(),
        "metric zero command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("valid metric zero JSON")
}

#[test]
fn horizontal_los_sight_adjustment_equals_solved_zero_angle() {
    let result = run_zero(0.0);
    let zero_angle_moa = result["zero_angle_moa"].as_f64().unwrap();
    let sight_adjustment_moa = result["sight_adjustment_moa"].as_f64().unwrap();

    assert!(
        (sight_adjustment_moa - zero_angle_moa).abs() < 1e-12,
        "horizontal LOS adjustment double-counted sight height: adjustment={sight_adjustment_moa} zero_angle={zero_angle_moa}"
    );
}

#[test]
fn elevated_target_sight_adjustment_removes_los_slope() {
    let target_height_yards = 0.25;
    let result = run_zero(target_height_yards);
    let zero_angle_moa = result["zero_angle_moa"].as_f64().unwrap();
    let actual = result["sight_adjustment_moa"].as_f64().unwrap();
    let expected = zero_angle_moa - target_height_yards / TARGET_DISTANCE_YARDS * 3437.75;

    assert!(
        (actual - expected).abs() < 1e-12,
        "sight adjustment used the wrong LOS frame: actual={actual} expected={expected}"
    );
}

#[test]
fn point_blank_range_ignores_the_initial_bore_to_sight_offset() {
    let result = run_zero(0.0);
    let point_blank_range = result["point_blank_range"].as_f64().unwrap();
    let solver_horizon_yards = TARGET_DISTANCE_YARDS * 3.0;

    assert!(
        point_blank_range > TARGET_DISTANCE_YARDS && point_blank_range < solver_horizon_yards,
        "point-blank range did not use the post-zero band exit: pbr={point_blank_range} yd zero={TARGET_DISTANCE_YARDS} yd horizon={solver_horizon_yards} yd"
    );
}

#[test]
fn zero_json_distances_follow_selected_units() {
    let imperial = run_zero(0.0);
    let metric = run_zero_metric();

    assert_eq!(imperial["units"], "imperial");
    assert_eq!(metric["units"], "metric");
    for field in [
        "zero_angle_degrees",
        "zero_angle_moa",
        "zero_angle_mrad",
        "sight_adjustment_moa",
    ] {
        let imperial_angle = imperial[field].as_f64().unwrap();
        let metric_angle = metric[field].as_f64().unwrap();
        assert!(
            (imperial_angle - metric_angle).abs() < 1e-12,
            "{field} changed across equivalent unit systems"
        );
    }
    for field in ["max_ordinate", "point_blank_range"] {
        let yards = imperial[field].as_f64().unwrap();
        let meters = metric[field].as_f64().unwrap();
        assert!(
            (yards * YARDS_TO_METERS - meters).abs() < 1e-9,
            "{field} did not convert consistently: {yards} yd versus {meters} m"
        );
    }
}
