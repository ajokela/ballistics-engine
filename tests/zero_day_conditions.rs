//! Integration tests for the `--auto-zero` zero-day condition overrides
//! (`--zero-velocity`, `--zero-temperature`, `--zero-pressure`, `--zero-humidity`,
//! `--zero-altitude`).
//!
//! These let the zero ANGLE be solved under the conditions/velocity the rifle was
//! actually zeroed in, while the trajectory runs under the current shot-day
//! conditions — the classic "zeroed on a cold day, shooting on a warm day" shift.
//!
//! JSON axis convention (imperial): `z` = downrange (yards), `y` = height (yards).
//! The zero angle tilts the whole arc, so a larger zero angle raises `y` at every
//! downrange point.

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

/// Run a 300 yd trajectory zeroed at 100 yd and return (downrange_yd, height_yd) points.
fn run_traj(extra: &[&str]) -> Vec<(f64, f64)> {
    let mut args: Vec<&str> = vec![
        "trajectory",
        "-v",
        "2650",
        "-b",
        "0.19",
        "-m",
        "77",
        "-d",
        "0.224",
        "--drag-model",
        "g7",
        "--max-range",
        "300",
        "--sight-height",
        "2.48",
        "--temperature",
        "59",
        "--pressure",
        "29.92",
        "--humidity",
        "50",
        "--auto-zero",
        "100",
        "--full",
        "-o",
        "json",
    ];
    args.extend_from_slice(extra);

    let output = Command::new(get_cli_binary())
        .args(&args)
        .output()
        .expect("Failed to execute command");
    assert!(
        output.status.success(),
        "trajectory command should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).expect("valid JSON");
    json["trajectory"]
        .as_array()
        .expect("trajectory array")
        .iter()
        .map(|p| (p["z"].as_f64().unwrap(), p["y"].as_f64().unwrap()))
        .collect()
}

/// Interpolate height (yd) at a given downrange (yd).
fn height_at(pts: &[(f64, f64)], zyd: f64) -> f64 {
    for w in pts.windows(2) {
        let (z1, y1) = w[0];
        let (z2, y2) = w[1];
        if z2 >= zyd && z2 > z1 {
            let t = (zyd - z1) / (z2 - z1);
            return y1 + t * (y2 - y1);
        }
    }
    pts.last().unwrap().1
}

/// Supplying zero-day conditions equal to the shot-day conditions must reproduce the
/// plain auto-zero result exactly — the override path is a no-op when nothing differs,
/// so existing callers are unaffected.
#[test]
fn zero_day_equal_to_shot_day_is_a_noop() {
    let baseline = run_traj(&[]);
    let explicit = run_traj(&[
        "--zero-velocity",
        "2650",
        "--zero-temperature",
        "59",
        "--zero-pressure",
        "29.92",
        "--zero-humidity",
        "50",
    ]);
    assert_eq!(
        baseline.len(),
        explicit.len(),
        "same number of trajectory points"
    );
    for zyd in [0.0, 50.0, 100.0, 200.0, 300.0] {
        let diff = (height_at(&baseline, zyd) - height_at(&explicit, zyd)).abs();
        // Well under a thousandth of an inch (1 yd = 36 in).
        assert!(
            diff < 1e-6,
            "explicit zero-day==shot-day should match baseline at {zyd} yd (diff {diff} yd)"
        );
    }
}

/// Zeroing under a COLDER, SLOWER load (more drop) bakes in more elevation. Firing the
/// current warmer/faster (flatter) load then prints HIGH — height must exceed the
/// baseline zero at every downrange point past the muzzle.
#[test]
fn zeroed_cold_and_slow_prints_high() {
    let baseline = run_traj(&[]);
    let cold = run_traj(&["--zero-velocity", "2500", "--zero-temperature", "0"]);
    for zyd in [50.0, 100.0, 200.0, 300.0] {
        let b = height_at(&baseline, zyd);
        let c = height_at(&cold, zyd);
        assert!(
            c > b,
            "cold/slow zero should fly higher than baseline at {zyd} yd (cold {c} <= base {b})"
        );
    }
}

/// Zeroing under a HOTTER, FASTER load (flatter) bakes in less elevation, so the current
/// load prints LOW relative to the baseline zero.
#[test]
fn zeroed_hot_and_fast_prints_low() {
    let baseline = run_traj(&[]);
    let hot = run_traj(&["--zero-velocity", "2800", "--zero-temperature", "100"]);
    for zyd in [50.0, 100.0, 200.0, 300.0] {
        let b = height_at(&baseline, zyd);
        let h = height_at(&hot, zyd);
        assert!(
            h < b,
            "hot/fast zero should fly lower than baseline at {zyd} yd (hot {h} >= base {b})"
        );
    }
}
