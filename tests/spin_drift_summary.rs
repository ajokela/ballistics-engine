//! The summary must report the same canonical Litz drift applied to trajectory points.

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

fn get_cli_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/debug/ballistics")
}

#[test]
fn spin_drift_summary_matches_final_lateral_trajectory() {
    let output = Command::new(get_cli_binary())
        .args([
            "--units",
            "metric",
            "trajectory",
            "-v",
            "800",
            "-b",
            "0.24",
            "-m",
            "11.34",
            "-d",
            "7.82",
            "--bullet-length",
            "33",
            "--twist-rate",
            "254",
            "--enable-spin-drift",
            "--max-range",
            "1500",
            "--ignore-ground-impact",
            "--temperature",
            "0",
            "--pressure",
            "700",
            "--humidity",
            "50",
            "--altitude",
            "3000",
            "--full",
            "-o",
            "json",
        ])
        .output()
        .expect("run spin-drift trajectory");
    assert!(
        output.status.success(),
        "trajectory failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let json: Value = serde_json::from_slice(&output.stdout).expect("valid trajectory JSON");
    let summary = json["spin_drift"].as_f64().expect("summary spin drift");
    let points = json["trajectory"].as_array().expect("full trajectory");
    let final_lateral = points
        .last()
        .and_then(|point| point["x"].as_f64())
        .expect("final lateral position");

    assert!(
        final_lateral > 1.0,
        "fixture must produce a material right-hand spin drift"
    );
    assert!(
        (summary - final_lateral).abs() < 1e-9,
        "summary drift diverged from canonical trajectory drift: summary={summary} m, final lateral={final_lateral} m"
    );
}
