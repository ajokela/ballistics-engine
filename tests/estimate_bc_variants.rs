//! Tests for multi-variant BC estimation: G1 / G7 drag models, fit against either a
//! drop curve or a velocity-retention curve.
//!
//! The library-level tests are round-trips: generate synthetic data from a *known* BC with
//! the engine's own forward model, then confirm the estimator recovers that BC — for both
//! drag models and both fit bases. The CLI test exercises the wired-up `estimate-bc` command
//! end to end.

use ballistics_engine::{
    estimate_bc_fit, BallisticInputs, BcFitMode, DragModel, TrajectoryPoint, TrajectorySolver,
};

// 77 gr .224 at 2650 fps, expressed in SI (m/s, kg, m).
const V: f64 = 2650.0 * 0.3048;
const M: f64 = 77.0 * 0.00006479891;
const D: f64 = 0.224 * 0.0254;
// 300 / 500 / 700 yd in meters.
const RANGES_M: [f64; 3] = [300.0 * 0.9144, 500.0 * 0.9144, 700.0 * 0.9144];

/// Forward model: solve a flat-fire trajectory for a known BC and drag model.
fn solve_ref(bc: f64, model: DragModel) -> Vec<TrajectoryPoint> {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = V;
    inputs.bc_value = bc;
    inputs.bc_type = model;
    inputs.bullet_mass = M;
    inputs.bullet_diameter = D;
    let mut solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
    solver.set_max_range(RANGES_M[2] * 1.5);
    solver.solve().expect("reference solve").points
}

/// Interpolate drop (m, `-y`) or speed (m/s) at a downrange distance (`position.x`).
fn interp(pts: &[TrajectoryPoint], target: f64, drop: bool) -> f64 {
    let val = |p: &TrajectoryPoint| {
        if drop {
            -p.position.y
        } else {
            p.velocity_magnitude
        }
    };
    for i in 1..pts.len() {
        if pts[i].position.x >= target {
            let p1 = &pts[i - 1];
            let p2 = &pts[i];
            let t = (target - p1.position.x) / (p2.position.x - p1.position.x);
            return val(p1) + t * (val(p2) - val(p1));
        }
    }
    val(pts.last().unwrap())
}

fn drop_data(bc: f64, model: DragModel) -> Vec<(f64, f64)> {
    let pts = solve_ref(bc, model);
    RANGES_M.iter().map(|&r| (r, interp(&pts, r, true))).collect()
}

fn velocity_data(bc: f64, model: DragModel) -> Vec<(f64, f64)> {
    let pts = solve_ref(bc, model);
    RANGES_M
        .iter()
        .map(|&r| (r, interp(&pts, r, false)))
        .collect()
}

/// A G7 BC must be recovered from a drop curve it generated.
#[test]
fn g7_recovered_from_drop() {
    let data = drop_data(0.19, DragModel::G7);
    let est = estimate_bc_fit(V, M, D, &data, DragModel::G7, BcFitMode::Drop).unwrap();
    assert!(
        (est.bc - 0.19).abs() < 0.01,
        "G7 drop-fit should recover ~0.19, got {}",
        est.bc
    );
}

/// A G7 BC must also be recovered from a velocity-retention curve it generated.
#[test]
fn g7_recovered_from_velocity() {
    let data = velocity_data(0.19, DragModel::G7);
    let est = estimate_bc_fit(V, M, D, &data, DragModel::G7, BcFitMode::Velocity).unwrap();
    assert!(
        (est.bc - 0.19).abs() < 0.01,
        "G7 velocity-fit should recover ~0.19, got {}",
        est.bc
    );
}

/// A G1 BC must be recovered from its own drop curve.
#[test]
fn g1_recovered_from_drop() {
    let data = drop_data(0.372, DragModel::G1);
    let est = estimate_bc_fit(V, M, D, &data, DragModel::G1, BcFitMode::Drop).unwrap();
    assert!(
        (est.bc - 0.372).abs() < 0.01,
        "G1 drop-fit should recover ~0.372, got {}",
        est.bc
    );
}

/// Fitting G7-generated data with the *wrong* model (G1) yields a different, larger BC —
/// G1 BCs run roughly 2x the G7 BC for a boat-tail. This guards against the two drag models
/// being silently conflated.
#[test]
fn g1_and_g7_fits_differ() {
    let data = drop_data(0.19, DragModel::G7);
    let g7 = estimate_bc_fit(V, M, D, &data, DragModel::G7, BcFitMode::Drop)
        .unwrap()
        .bc;
    let g1 = estimate_bc_fit(V, M, D, &data, DragModel::G1, BcFitMode::Drop)
        .unwrap()
        .bc;
    assert!(
        g1 > g7 * 1.5,
        "G1 fit ({g1}) should be markedly larger than G7 fit ({g7})"
    );
}

/// Drop-fit and velocity-fit of the same bullet must agree on the BC.
#[test]
fn drop_and_velocity_fits_agree() {
    let drop = estimate_bc_fit(
        V,
        M,
        D,
        &drop_data(0.22, DragModel::G7),
        DragModel::G7,
        BcFitMode::Drop,
    )
    .unwrap()
    .bc;
    let vel = estimate_bc_fit(
        V,
        M,
        D,
        &velocity_data(0.22, DragModel::G7),
        DragModel::G7,
        BcFitMode::Velocity,
    )
    .unwrap()
    .bc;
    assert!(
        (drop - vel).abs() < 0.01,
        "drop-fit ({drop}) and velocity-fit ({vel}) should agree"
    );
}

// ---- CLI end-to-end ----

use std::path::PathBuf;
use std::process::Command;

fn cli() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("target");
    p.push("debug");
    p.push("ballistics");
    if !p.exists() {
        p.pop();
        p.pop();
        p.push("release");
        p.push("ballistics");
    }
    p
}

/// `estimate-bc --drag-model both` on drop data prints both a G1 and a G7 variant, with the
/// G7 BC smaller than the G1 BC (Bero's data: a 77 gr .224 fits ~0.19 G7 / ~0.37 G1).
#[test]
fn cli_estimate_bc_both_models_json() {
    let out = Command::new(cli())
        .args([
            "estimate-bc",
            "-v",
            "2650",
            "-m",
            "77",
            "-d",
            "0.224",
            "--data",
            "300,29.0;500,89.9;700,204.6",
            "--drag-model",
            "both",
            "-o",
            "json",
        ])
        .output()
        .expect("run estimate-bc");
    assert!(
        out.status.success(),
        "estimate-bc should succeed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&out.stdout).expect("valid JSON");
    let variants = json["variants"].as_array().expect("variants array");
    assert_eq!(variants.len(), 2, "G1 + G7 = two variants");

    let bc_of = |model: &str| -> f64 {
        variants
            .iter()
            .find(|v| v["drag_model"] == model)
            .and_then(|v| v["estimated_bc"].as_f64())
            .unwrap_or_else(|| panic!("no {model} variant"))
    };
    let g1 = bc_of("G1");
    let g7 = bc_of("G7");
    assert!(
        (0.34..0.40).contains(&g1),
        "G1 BC for this bullet should be ~0.37, got {g1}"
    );
    assert!(
        (0.17..0.21).contains(&g7),
        "G7 BC for this bullet should be ~0.19, got {g7}"
    );
    assert!(g1 > g7, "G1 BC ({g1}) should exceed G7 BC ({g7})");
}
