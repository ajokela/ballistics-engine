//! Auto-zero must solve with the same velocity-dependent BC configuration as the flight.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const ZERO_RANGE_YARDS: f64 = 600.0;
const EXPECTED_ZERO_HEIGHT_YARDS: f64 = (5.0 * 12.0 + 2.0) / 36.0;

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

fn zero_miss_inches(velocity: &str, bc: &str, mass: &str, diameter: &str, bc_option: &str) -> f64 {
    let output = Command::new(get_cli_binary())
        .args([
            "trajectory",
            "--velocity",
            velocity,
            "--bc",
            bc,
            "--mass",
            mass,
            "--diameter",
            diameter,
            "--drag-model",
            "g7",
            "--auto-zero",
            "600",
            "--max-range",
            "650",
            "--bore-height",
            "5",
            "--sight-height",
            "2",
            "--ignore-ground-impact",
            "--time-step",
            "0.005",
            "--full",
            "--output",
            "json",
            bc_option,
        ])
        .output()
        .expect("Failed to execute auto-zero trajectory");
    assert!(
        output.status.success(),
        "auto-zero trajectory failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let result: Value = serde_json::from_slice(&output.stdout).expect("valid trajectory JSON");
    let points = result["trajectory"].as_array().unwrap();
    let height = points
        .windows(2)
        .find_map(|window| {
            let z1 = window[0]["z"].as_f64().unwrap();
            let z2 = window[1]["z"].as_f64().unwrap();
            if z2 >= ZERO_RANGE_YARDS && z2 > z1 {
                let y1 = window[0]["y"].as_f64().unwrap();
                let y2 = window[1]["y"].as_f64().unwrap();
                let fraction = (ZERO_RANGE_YARDS - z1) / (z2 - z1);
                Some(y1 + fraction * (y2 - y1))
            } else {
                None
            }
        })
        .expect("trajectory should cross the requested zero range");

    (height - EXPECTED_ZERO_HEIGHT_YARDS) * 36.0
}

#[test]
fn estimator_segments_are_shared_by_zero_and_flight() {
    let miss = zero_miss_inches("2650", "0.19", "77", "0.224", "--use-bc-segments");
    assert!(
        miss.abs() < 0.01,
        "estimator-segmented flight missed its requested zero by {miss:.4} inches"
    );
}

#[test]
fn cluster_bc_is_shared_by_zero_and_flight() {
    let miss = zero_miss_inches("3600", "0.075", "40", "0.20", "--use-cluster-bc");
    assert!(
        miss.abs() < 0.01,
        "cluster-BC flight missed its requested zero by {miss:.4} inches"
    );
}
