//! `zero` must solve targets below the shooter, including ones deeper than the solver's
//! default 100 m ground stop (the handbook's "250 yards below the shooter" example).

use serde_json::Value;
use std::path::PathBuf;
use std::process::Command;

const TARGET_DISTANCE_YARDS: f64 = 500.0;

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
            "0.462",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--target-distance",
            &TARGET_DISTANCE_YARDS.to_string(),
            &format!("--target-height={target_height_yards}"),
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute zero command");
    assert!(
        output.status.success(),
        "zero --target-height {target_height_yards} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("valid zero JSON")
}

fn line_of_sight_degrees(target_height_yards: f64) -> f64 {
    (target_height_yards / TARGET_DISTANCE_YARDS)
        .atan()
        .to_degrees()
}

/// Bore angle above the line of sight to the target, in degrees.
fn hold_degrees(target_height_yards: f64) -> f64 {
    run_zero(target_height_yards)["zero_angle_degrees"]
        .as_f64()
        .unwrap()
        - line_of_sight_degrees(target_height_yards)
}

#[test]
fn targets_below_the_shooter_get_a_bore_below_the_line_of_sight_to_them() {
    let level_hold = hold_degrees(0.0);

    // -5 yd and -100 yd sit above the default ground stop; -250 yd (about 229 m) is below it.
    for height in [-5.0, -100.0, -250.0] {
        let angle = run_zero(height)["zero_angle_degrees"].as_f64().unwrap();
        let hold = angle - line_of_sight_degrees(height);

        assert!(
            angle < 0.0,
            "a target {height} yd down needs a downward bore: angle={angle} deg"
        );
        // The bore still points above the line of sight to cover drop, by about the level
        // hold (rifleman's rule: within a few percent here, more for the longer slant range).
        assert!(
            hold > 0.0 && (hold - level_hold).abs() / level_hold < 0.05,
            "hold over the line of sight out of range at {height} yd: hold={hold} deg \
             level={level_hold} deg"
        );
    }
}

#[test]
fn a_deep_target_and_its_mirror_above_need_similar_holds() {
    // Gravity acts across the line of sight with the same cosine uphill and downhill, so the
    // two holds differ only by second-order effects. This catches a solution that lands on
    // the ground stop instead of on the target.
    let up = hold_degrees(250.0);
    let down = hold_degrees(-250.0);
    assert!(
        (up - down).abs() / up < 0.01,
        "uphill and downhill holds at 26.6 deg diverge: up={up} deg down={down} deg"
    );
}
