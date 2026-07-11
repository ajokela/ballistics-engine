//! MBA-1134 parity test: the empirical Litz spin-drift model must be the single canonical
//! lateral-drift path across all three solver families. A fixed .308 175gr shot fed through
//!   (a) cli_api::TrajectorySolver::solve()             (validated RK45 solver)
//!   (b) monte_carlo::solve_trajectory_for_monte_carlo  (fast_integrate, fixed-step RK4)
//!   (c) fast_trajectory::fast_integrate_with_segments  (RK45 via integrate_trajectory)
//! must agree on the final lateral (McCoy Z, windage) within < 0.3 in (0.00762 m). This locks the
//! three families together so no path silently uses a different (or stacked) spin-drift model.

use ballistics_engine::atmosphere::{calculate_atmosphere, resolve_station_conditions};
use ballistics_engine::fast_trajectory::{
    fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
};
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo;
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

const TARGET_M: f64 = 914.4; // 1000 yd
const TOL_M: f64 = 0.00762; // 0.3 in

/// Fixed parity bullet: .308, 175 gr, 1:10 right-hand twist, ~2625 fps, G7 BC, sea-level
/// standard atmosphere, no wind. Spin drift + advanced effects on.
fn parity_inputs() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0; // m/s (~2625 fps)
    i.bullet_mass = 175.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // m
    i.bullet_length = 1.24 * 0.0254; // m (175 SMK-class)
    i.caliber_inches = 0.308;
    i.weight_grains = 175.0;
    i.bc_value = 0.243;
    i.bc_type = DragModel::G7;
    i.twist_rate = 10.0; // 1:10"
    i.is_twist_right = true;
    i.muzzle_angle = 0.0; // flat launch; identical for all three paths
    i.target_distance = TARGET_M;
    // Sea-level standard atmosphere (matches AtmosphericConditions::default()).
    i.altitude = 0.0;
    i.temperature = 15.0;
    i.pressure = 1013.25;
    i.humidity = 0.5;
    // Canonical spin-drift + advanced effects.
    i.use_enhanced_spin_drift = true;
    i.enable_advanced_effects = true;
    i
}

/// (a) The validated cli_api solver: final lateral is the last trajectory point's McCoy Z.
fn lateral_cli_api(i: &BallisticInputs) -> f64 {
    let mut solver = TrajectorySolver::new(
        i.clone(),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(TARGET_M);
    let r = solver.solve().expect("cli_api solve");
    r.points.last().unwrap().position.z
}

/// (b) Monte-Carlo kernel: `wind_drift` is the total lateral (final integrated lateral + the Litz
/// spin-drift post-process). With no wind/Coriolis/Magnus the integrated lateral is ~0, so this is
/// the spin drift.
fn lateral_monte_carlo(i: &BallisticInputs) -> f64 {
    let out = solve_trajectory_for_monte_carlo(i).expect("monte carlo solve");
    out.wind_drift
}

/// (c) fast_integrate_with_segments: build the same physical shot the Monte-Carlo path builds, then
/// read the endpoint lateral (which now carries the canonical Litz post-process).
fn lateral_fast_segments(i: &BallisticInputs) -> f64 {
    let (resolved_temp_c, resolved_pressure_hpa) =
        resolve_station_conditions(i.temperature, i.pressure, i.altitude);
    let (air_density, _) = calculate_atmosphere(
        i.altitude,
        Some(resolved_temp_c),
        Some(resolved_pressure_hpa),
        i.humidity_percent(),
    );
    let base_ratio = air_density / 1.225;

    let v = i.muzzle_velocity;
    let a = i.muzzle_angle;
    let sight_position = i.muzzle_height + i.sight_height;
    let initial_state = [0.0, i.muzzle_height, 0.0, v * a.cos(), v * a.sin(), 0.0];
    let params = FastIntegrationParams {
        initial_state,
        t_span: (0.0, 30.0),
        horiz: TARGET_M,
        vert: sight_position,
        atmo_params: (i.altitude, resolved_temp_c, resolved_pressure_hpa, base_ratio),
        atmo_sock: None,
    };
    let sol = fast_integrate_with_segments(i, vec![], params);
    assert!(sol.success, "fast_integrate_with_segments should succeed");
    let n = sol.t.len();
    sol.y[2][n - 1] // McCoy Z (lateral)
}

/// (d) plain fast_integrate — the SINGLE-SHOT fast path (solve_trajectory_rust / the API).
/// MBA-1134 regression guard: this path lost spin drift when the in-integration term was removed
/// and the Litz post-process was only added to the MC / segmented paths. It must now agree too.
fn lateral_fast_integrate(i: &BallisticInputs) -> f64 {
    let (resolved_temp_c, resolved_pressure_hpa) =
        resolve_station_conditions(i.temperature, i.pressure, i.altitude);
    let (air_density, _) = calculate_atmosphere(
        i.altitude,
        Some(resolved_temp_c),
        Some(resolved_pressure_hpa),
        i.humidity_percent(),
    );
    let base_ratio = air_density / 1.225;

    let v = i.muzzle_velocity;
    let a = i.muzzle_angle;
    let sight_position = i.muzzle_height + i.sight_height;
    let initial_state = [0.0, i.muzzle_height, 0.0, v * a.cos(), v * a.sin(), 0.0];
    let params = FastIntegrationParams {
        initial_state,
        t_span: (0.0, 30.0),
        horiz: TARGET_M,
        vert: sight_position,
        atmo_params: (i.altitude, resolved_temp_c, resolved_pressure_hpa, base_ratio),
        atmo_sock: None,
    };
    let sol = fast_integrate(i, &WindSock::new(vec![]), params);
    assert!(sol.success, "fast_integrate should succeed");
    let n = sol.t.len();
    sol.y[2][n - 1] // McCoy Z (lateral)
}

#[test]
fn litz_spin_drift_parity_across_three_solvers() {
    let i = parity_inputs();

    let a = lateral_cli_api(&i);
    let b = lateral_monte_carlo(&i);
    let c = lateral_fast_segments(&i);
    let d = lateral_fast_integrate(&i);

    let to_in = |m: f64| m / 0.0254;
    println!("MBA-1134 spin-drift parity (.308 175gr, 1:10 RH, 1000 yd, no wind):");
    println!("  (a) cli_api TrajectorySolver : {:.6} m ({:.4} in)", a, to_in(a));
    println!("  (b) monte_carlo fast_integrate: {:.6} m ({:.4} in)", b, to_in(b));
    println!("  (c) fast_integrate_w_segments : {:.6} m ({:.4} in)", c, to_in(c));
    println!("  (d) plain fast_integrate      : {:.6} m ({:.4} in)", d, to_in(d));

    // Right-hand twist => positive (rightward) drift on all four.
    assert!(a > 0.0, "cli_api drift should be positive (RH twist), got {a}");
    assert!(b > 0.0, "monte_carlo drift should be positive (RH twist), got {b}");
    assert!(c > 0.0, "fast_segments drift should be positive (RH twist), got {c}");
    assert!(d > 0.0, "fast_integrate drift should be positive (RH twist), got {d}");

    for (name, x, y) in [
        ("cli_api vs monte_carlo", a, b),
        ("cli_api vs fast_segments", a, c),
        ("cli_api vs fast_integrate", a, d),
        ("monte_carlo vs fast_segments", b, c),
        ("monte_carlo vs fast_integrate", b, d),
        ("fast_segments vs fast_integrate", c, d),
    ] {
        assert!(
            (x - y).abs() < TOL_M,
            "{name} lateral disagree: {:.4} in vs {:.4} in (Δ {:.4} in > 0.3 in)",
            to_in(x),
            to_in(y),
            to_in((x - y).abs())
        );
    }
}

#[test]
fn litz_spin_drift_flips_with_twist_hand_all_paths() {
    // A left-hand twist must mirror the drift on every path (canonical sign convention).
    let mut right = parity_inputs();
    right.is_twist_right = true;
    let mut left = parity_inputs();
    left.is_twist_right = false;

    for (name, f) in [
        ("cli_api", lateral_cli_api as fn(&BallisticInputs) -> f64),
        ("monte_carlo", lateral_monte_carlo),
        ("fast_segments", lateral_fast_segments),
        ("fast_integrate", lateral_fast_integrate),
    ] {
        let r = f(&right);
        let l = f(&left);
        assert!(
            r > 0.0 && l < 0.0,
            "{name}: RH twist should drift +Z and LH -Z, got R={r}, L={l}"
        );
        assert!(
            (r + l).abs() < TOL_M,
            "{name}: RH and LH drift should be equal and opposite, got R={r}, L={l}"
        );
    }
}
