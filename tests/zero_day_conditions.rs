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

/// Run a deliberately temperature-sensitive load with a 500 yd zero. The stated 2650 fps
/// velocity is referenced to 70 F, so the 110 F shot-day velocity resolves to 2850 fps.
fn run_linear_powder_zero(extra: &[&str]) -> (Vec<(f64, f64)>, f64) {
    let mut args = vec![
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
        "--temperature",
        "110",
        "--pressure",
        "29.92",
        "--powder-temp",
        "70",
        "--use-powder-sensitivity",
        "--powder-temp-sensitivity",
        "5",
        "--auto-zero",
        "500",
        "--max-range",
        "600",
        "--bore-height",
        "5",
        "--sight-height",
        "2.48",
        "--ignore-ground-impact",
        "--full",
        "-o",
        "json",
    ];
    args.extend_from_slice(extra);

    let output = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run linear-powder auto-zero trajectory");
    assert!(
        output.status.success(),
        "trajectory command should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: Value = serde_json::from_slice(&output.stdout).expect("valid trajectory JSON");
    let trajectory = json["trajectory"].as_array().expect("trajectory array");
    let muzzle_velocity = trajectory[0]["velocity"].as_f64().expect("muzzle velocity");
    let points = trajectory
        .iter()
        .map(|point| {
            (
                point["z"].as_f64().expect("downrange"),
                point["y"].as_f64().expect("height"),
            )
        })
        .collect();
    (points, muzzle_velocity)
}

/// Run a curve-resolved load with 32 F air and ammunition held at 68 F. The curve gives
/// 2706 fps at the explicit powder temperature.
fn run_curve_powder_zero(extra: &[&str]) -> (Vec<(f64, f64)>, f64) {
    let mut args = vec![
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
        "--temperature",
        "32",
        "--powder-temp",
        "68",
        "--powder-temp-curve",
        "32:2650,77:2720",
        "--auto-zero",
        "100",
        "--max-range",
        "300",
        "--bore-height",
        "5",
        "--sight-height",
        "2.48",
        "--ignore-ground-impact",
        "--full",
        "-o",
        "json",
    ];
    args.extend_from_slice(extra);

    let output = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run powder-curve auto-zero trajectory");
    assert!(
        output.status.success(),
        "trajectory command should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: Value = serde_json::from_slice(&output.stdout).expect("valid trajectory JSON");
    let trajectory = json["trajectory"].as_array().expect("trajectory array");
    let muzzle_velocity = trajectory[0]["velocity"].as_f64().expect("muzzle velocity");
    let points = trajectory
        .iter()
        .map(|point| {
            (
                point["z"].as_f64().expect("downrange"),
                point["y"].as_f64().expect("height"),
            )
        })
        .collect();
    (points, muzzle_velocity)
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

/// Regression (Bero/David): extreme zero-day atmosphere at a SHORT (100 yd) zero must
/// still change the trajectory. A 100 yd zero is nearly density-independent so the effect
/// is small, but the zero solve's convergence must be fine enough that it doesn't quantize
/// two very different atmospheres to a byte-identical zero angle (it used to converge at
/// only 1 mm of height, which rounded the effect to zero at short zeros).
#[test]
fn zero_day_atmosphere_resolves_at_short_zero() {
    let base = run_traj(&[]); // run_traj already zeroes at 100 yd
    let extreme = run_traj(&[
        "--zero-temperature", "20", "--zero-pressure", "20", "--zero-humidity", "90",
        "--zero-altitude", "12000",
    ]);
    // Measure near the far end of run_traj's 300 yd range.
    let far = 290.0;
    let diff = (height_at(&base, far) - height_at(&extreme, far)).abs();
    assert!(
        diff > 0.001, // yards (~0.036"); the real effect here is ~0.1"
        "extreme zero-day atmosphere must change the come-up at {far} yd, got diff {diff} yd"
    );
}

#[test]
fn linear_powder_model_is_shared_by_auto_zero_and_flight() {
    let (implicit, implicit_velocity) = run_linear_powder_zero(&[]);
    let (explicit, explicit_velocity) = run_linear_powder_zero(&["--zero-velocity", "2850"]);

    assert!((implicit_velocity - 2850.0).abs() < 1e-6);
    assert!((explicit_velocity - 2850.0).abs() < 1e-6);
    let implicit_height = height_at(&implicit, 500.0);
    let explicit_height = height_at(&explicit, 500.0);
    let sight_line_height = (5.0 * 12.0 + 2.48) / 36.0;
    assert!(
        (implicit_height - explicit_height).abs() < 1e-6,
        "implicit linear zero must match explicit 2850 fps zero: implicit={implicit_height} yd, explicit={explicit_height} yd"
    );
    assert!(
        (implicit_height - sight_line_height).abs() < 1e-4,
        "temperature-resolved flight missed its 500 yd zero: height={implicit_height} yd, sight line={sight_line_height} yd"
    );
}

#[test]
fn linear_powder_zero_uses_zero_day_temperature() {
    // At 30 F the zero-day velocity is 2650 + 5 * (30 - 70) = 2450 fps, while the
    // shot-day trajectory remains at 2850 fps under its 110 F conditions.
    let (implicit, _) = run_linear_powder_zero(&["--zero-temperature", "30"]);
    let (explicit, _) =
        run_linear_powder_zero(&["--zero-temperature", "30", "--zero-velocity", "2450"]);

    let implicit_height = height_at(&implicit, 500.0);
    let explicit_height = height_at(&explicit, 500.0);
    assert!(
        (implicit_height - explicit_height).abs() < 1e-6,
        "linear zero must use the zero-day temperature: implicit={implicit_height} yd, explicit={explicit_height} yd"
    );
}

#[test]
fn auto_zero_inherits_explicit_shot_day_powder_temperature() {
    let (inherited, inherited_velocity) = run_curve_powder_zero(&[]);
    let (explicit, explicit_velocity) = run_curve_powder_zero(&["--zero-powder-temp", "68"]);

    assert!((inherited_velocity - 2706.0).abs() < 1e-6);
    assert!((explicit_velocity - 2706.0).abs() < 1e-6);
    for distance in [100.0, 200.0, 300.0] {
        let inherited_height = height_at(&inherited, distance);
        let explicit_height = height_at(&explicit, distance);
        assert!(
            (inherited_height - explicit_height).abs() < 1e-6,
            "omitting --zero-powder-temp must inherit 68 F at {distance} yd: inherited={inherited_height}, explicit={explicit_height}"
        );
    }
}

#[test]
fn explicit_zero_air_temperature_remains_the_curve_fallback() {
    let (implicit, _) = run_curve_powder_zero(&["--zero-temperature", "77"]);
    let (explicit, _) =
        run_curve_powder_zero(&["--zero-temperature", "77", "--zero-powder-temp", "77"]);

    for distance in [100.0, 200.0, 300.0] {
        let implicit_height = height_at(&implicit, distance);
        let explicit_height = height_at(&explicit, distance);
        assert!(
            (implicit_height - explicit_height).abs() < 1e-6,
            "explicit zero air must remain the powder-curve fallback at {distance} yd: implicit={implicit_height}, explicit={explicit_height}"
        );
    }
}
