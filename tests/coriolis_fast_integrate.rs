// MBA-957 regression guard: fast_integrate (the Monte Carlo / Python-binding path via
// solve_trajectory_for_monte_carlo) previously applied NO Coriolis, while the cli_api
// TrajectorySolver applies the validated -2 Omega x v (MBA-938, checked vs py_ballisticcalc).
// This test pins the two integrators to AGREE on the Coriolis lateral deflection: same sign
// (Northern hemisphere drifts RIGHT / +Z, Southern flips LEFT) and comparable magnitude.
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

const TARGET_M: f64 = 731.5; // ~800 yd — long enough that Coriolis deflection is well-resolved

fn shot(latitude: f64, enable_coriolis: bool) -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.target_distance = TARGET_M;
    i.muzzle_velocity = 823.0; // m/s
    i.bullet_mass = 168.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // m
    i.bullet_length = 1.215 * 0.0254;
    i.caliber_inches = 0.308;
    i.weight_grains = 168.0;
    i.bc_value = 0.475;
    i.bc_type = DragModel::G1;
    i.muzzle_angle = 0.012; // radians — enough elevation to carry to ~731 m
    i.sight_height = 0.05;
    i.muzzle_height = 0.0;
    i.temperature = 15.0;
    i.pressure = 1013.25;
    i.enable_coriolis = enable_coriolis;
    i.latitude = Some(latitude);
    i.azimuth_angle = 0.0; // 0 = North shot
    i
}

/// Lateral (Z) drift at the target from the fast/MC integrator.
fn fast_drift(latitude: f64, coriolis: bool) -> f64 {
    solve_trajectory_for_monte_carlo(&shot(latitude, coriolis))
        .expect("fast_integrate solve reached target")
        .wind_drift
}

/// Lateral (Z) drift at the target from the validated cli_api solver.
fn cli_drift(latitude: f64, coriolis: bool) -> f64 {
    let i = shot(latitude, coriolis);
    let mut solver =
        TrajectorySolver::new(i.clone(), WindConditions::default(), AtmosphericConditions::default());
    solver.set_max_range(TARGET_M);
    let r = solver.solve().expect("cli_api solve");
    r.points.last().unwrap().position.z
}

#[test]
fn fast_integrate_coriolis_matches_cli_api() {
    // Isolate the Coriolis contribution: (coriolis ON) - (coriolis OFF) for each integrator,
    // which cancels any baseline lateral and any integrator-specific zero offset.
    let fast_defl = fast_drift(45.0, true) - fast_drift(45.0, false);
    let cli_defl = cli_drift(45.0, true) - cli_drift(45.0, false);

    // 1. The fix actually added Coriolis to the fast path (was ~0 before).
    assert!(
        fast_defl.abs() > 1e-4,
        "fast_integrate Coriolis deflection should be non-trivial, got {fast_defl:.6} m"
    );
    // 2. Northern hemisphere North shot drifts RIGHT (+Z) in BOTH solvers.
    assert!(cli_defl > 0.0, "cli Coriolis deflection should be +Z (right), got {cli_defl:.6} m");
    assert!(fast_defl > 0.0, "fast Coriolis deflection should be +Z (right), got {fast_defl:.6} m");
    // 3. The two independent integrators agree on magnitude (loose tolerance: RK4 fixed-step
    //    fast path vs RK45 adaptive cli path differ slightly in TOF, hence the deflection).
    let rel = (fast_defl - cli_defl).abs() / cli_defl;
    assert!(
        rel < 0.30,
        "fast vs cli Coriolis deflection disagree by {:.1}% (fast={:.4} m, cli={:.4} m)",
        rel * 100.0,
        fast_defl,
        cli_defl
    );
}

#[test]
fn fast_integrate_coriolis_flips_with_hemisphere() {
    // Southern hemisphere flips the lateral Coriolis drift to LEFT (-Z).
    let north = fast_drift(45.0, true) - fast_drift(45.0, false);
    let south = fast_drift(-45.0, true) - fast_drift(-45.0, false);
    assert!(north > 0.0, "N hemisphere should drift right (+Z), got {north:.6} m");
    assert!(south < 0.0, "S hemisphere should drift left (-Z), got {south:.6} m");
}
