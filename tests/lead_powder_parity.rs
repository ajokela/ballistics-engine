//! Integration tests for `lead` powder-temperature parity (MBA-1325).
//!
//! `lead` gained the same `--use-powder-sensitivity` / `--powder-temp-sensitivity` /
//! `--powder-temp` / `--powder-temp-curve` flags as `trajectory`, plumbed identically
//! into `BallisticInputs` (both commands' inputs flow through the SAME
//! `cli_api::TrajectorySolver::new`, which resolves the powder correction). A
//! powder-corrected `lead` run must therefore reproduce exactly the same holds as
//! manually resolving the muzzle velocity via `trajectory` first and feeding it back
//! in with `-v` (no powder flags) — the "extract resolved velocity via the trajectory
//! path" technique from tests/powder_temp_curve.rs.

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

fn run_json(args: &[&str]) -> Value {
    let out = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run");
    assert!(
        out.status.success(),
        "command failed ({:?}): {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_str(&String::from_utf8_lossy(&out.stdout)).expect("json")
}

/// Muzzle velocity (fps) `trajectory` resolves for a given powder configuration — the
/// first trajectory point's velocity, same technique as tests/powder_temp_curve.rs.
fn resolved_velocity(temp_f: &str, extra: &[&str]) -> f64 {
    let mut args: Vec<&str> = vec![
        "trajectory", "-v", "2700", "-b", "0.19", "-m", "77", "-d", "0.224", "--drag-model", "g7",
        "--max-range", "5", "--sight-height", "2.48", "--temperature", temp_f, "--full", "-o",
        "json",
    ];
    args.extend_from_slice(extra);
    let json = run_json(&args);
    json["trajectory"][0]["velocity"].as_f64().expect("v0")
}

const CURVE: &str = "40:2620,70:2700,100:2760";

fn base_lead_args<'a>(temp_f: &'a str, velocity: &'a str) -> Vec<&'a str> {
    vec![
        "lead", "-v", velocity, "-b", "0.19", "-m", "77", "-d", "0.224", "--drag-model", "g7",
        "--sight-height", "2.48", "--temperature", temp_f, "--target-speed", "5", "--start",
        "300", "--end", "300", "-o", "json",
    ]
}

/// `lead --powder-temp-curve` at 55F (the midpoint of 40/70 -> ~2660 fps, NOT a
/// measured curve point, so this is a genuine correction) reproduces the exact same
/// lead solution as `lead -v 2660...` (no curve) once the velocity is resolved.
#[test]
fn lead_with_powder_curve_matches_manually_resolved_velocity() {
    let resolved = resolved_velocity("55", &["--powder-temp-curve", CURVE]);
    assert!(
        (resolved - 2660.0).abs() < 0.5,
        "sanity: 55F should interpolate to ~2660 fps, got {resolved}"
    );

    let mut with_curve_args = base_lead_args("55", "2700");
    with_curve_args.extend(["--powder-temp-curve", CURVE]);
    let with_curve = run_json(&with_curve_args);

    let resolved_str = resolved.to_string();
    let manual = run_json(&base_lead_args("55", &resolved_str));

    for field in ["lead_mil", "lead_moa", "lead", "tof_s", "intercept_range", "iterations"] {
        assert_eq!(
            with_curve["rows"][0][field], manual["rows"][0][field],
            "field '{field}' differs between curve-corrected and manually-resolved lead runs"
        );
    }
}

/// Same parity check for the linear `--use-powder-sensitivity` model (not just the
/// data-driven curve): 68F air, 32F reference, 1.0 fps/F sensitivity -> 68-32=36 F
/// delta -> resolved velocity = 2700 + 36 = 2736 fps.
#[test]
fn lead_with_linear_powder_sensitivity_matches_manually_resolved_velocity() {
    let resolved = resolved_velocity(
        "68",
        &[
            "--use-powder-sensitivity",
            "--powder-temp",
            "32",
            "--powder-temp-sensitivity",
            "1.0",
        ],
    );
    assert!(
        (resolved - 2736.0).abs() < 1e-6,
        "sanity: 68F with 32F reference at 1.0 fps/F should resolve to 2736 fps, got {resolved}"
    );

    let mut with_sensitivity_args = base_lead_args("68", "2700");
    with_sensitivity_args.extend([
        "--use-powder-sensitivity",
        "--powder-temp",
        "32",
        "--powder-temp-sensitivity",
        "1.0",
    ]);
    let with_sensitivity = run_json(&with_sensitivity_args);

    let resolved_str = resolved.to_string();
    let manual = run_json(&base_lead_args("68", &resolved_str));

    for field in ["lead_mil", "lead_moa", "lead", "tof_s"] {
        assert_eq!(
            with_sensitivity["rows"][0][field], manual["rows"][0][field],
            "field '{field}' differs between sensitivity-corrected and manually-resolved lead runs"
        );
    }
}

/// Omitting every powder flag must use `-v` verbatim (regression risk: zero) — same
/// invariant `trajectory` already guarantees (tests/powder_temp_curve.rs
/// `no_curve_uses_velocity_verbatim`), now also true through `lead`.
#[test]
fn lead_with_no_powder_flags_uses_velocity_verbatim() {
    let v = resolved_velocity("100", &[]);
    assert!(
        (v - 2700.0).abs() < 1e-6,
        "no powder flags should keep -v 2700 verbatim regardless of --temperature, got {v}"
    );
}

/// `lead --help` advertises the new flags (documentation/plumbing sanity check).
#[test]
fn lead_help_lists_new_powder_flags() {
    let out = Command::new(get_cli_binary())
        .args(["lead", "--help"])
        .output()
        .expect("run");
    let help = String::from_utf8_lossy(&out.stdout);
    for flag in [
        "--use-powder-sensitivity",
        "--powder-temp-sensitivity",
        "--powder-temp",
        "--powder-temp-curve",
    ] {
        assert!(help.contains(flag), "`lead --help` is missing {flag}:\n{help}");
    }
}
