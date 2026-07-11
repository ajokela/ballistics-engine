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
    let target_distance_m = TARGET_DISTANCE_YARDS * YARDS_TO_METERS;
    let solver_horizon_m = target_distance_m * 3.0;

    assert!(
        point_blank_range > target_distance_m && point_blank_range < solver_horizon_m,
        "point-blank range did not use the post-zero band exit: pbr={point_blank_range} m zero={target_distance_m} m horizon={solver_horizon_m} m"
    );
}
