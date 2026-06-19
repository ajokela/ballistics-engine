//! MBA-959: verify aerodynamic jump is wired into the solver as an opt-in,
//! default-off muzzle launch-angle perturbation.

use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions,
};
use std::f64::consts::PI;

/// Build a solver for a representative .308-class load. `crosswind_mps` is applied
/// as a pure lateral (McCoy Z) wind so it drives aerodynamic jump.
fn solve(enable_aj: bool, crosswind_mps: f64) -> ballistics_engine::TrajectoryResult {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 800.0; // m/s
    inputs.muzzle_angle = 0.01; // small positive elevation so it carries downrange
    inputs.target_distance = 500.0;
    inputs.twist_rate = 11.0; // 1:11"
    inputs.is_twist_right = true;
    inputs.enable_aerodynamic_jump = enable_aj;

    // direction = PI/2 -> wind is purely lateral (crosswind), no head/tail component.
    let wind = WindConditions {
        speed: crosswind_mps,
        direction: PI / 2.0,
    };
    let atmo = AtmosphericConditions::default();
    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(600.0);
    solver.solve().expect("solve should succeed")
}

fn vertical_at_500(r: &ballistics_engine::TrajectoryResult) -> f64 {
    r.position_at_range(500.0)
        .expect("trajectory should reach 500 m")
        .y
}

#[test]
fn disabled_by_default_reports_no_jump() {
    let r = solve(false, 10.0);
    assert!(
        r.aerodynamic_jump.is_none(),
        "AJ must be None when the flag is off"
    );
}

#[test]
fn enabled_reports_jump_components() {
    let r = solve(true, 10.0);
    let aj = r
        .aerodynamic_jump
        .expect("AJ components must be present when enabled");
    assert!(
        aj.vertical_jump_moa.is_finite() && aj.vertical_jump_moa.abs() > 0.0,
        "crosswind should produce a nonzero vertical jump, got {}",
        aj.vertical_jump_moa
    );
}

#[test]
fn jump_shifts_vertical_impact_with_crosswind() {
    let off = vertical_at_500(&solve(false, 10.0));
    let on = vertical_at_500(&solve(true, 10.0));
    assert!(
        (on - off).abs() > 1e-6,
        "enabling AJ with a crosswind should shift the vertical impact (off={off}, on={on})"
    );
}

const MS_TO_MPH: f64 = 2.236_936_292_054_4;

#[test]
fn litz_magnitude_matches_engine_sg_exactly() {
    // The engine's vertical jump must equal Litz's Y = (0.01*Sg - 0.0024*L + 0.032)
    // * crosswind_mph, computed from the engine's OWN Miller Sg.
    use ballistics_engine::aerodynamic_jump::litz_crosswind_jump_moa;
    use ballistics_engine::stability::compute_stability_coefficient;

    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 790.0;
    inputs.bullet_diameter = 0.00782; // ~.308"
    inputs.bullet_length = 0.0312; // ~4.0 cal
    inputs.bullet_mass = 0.01134; // ~175 gr
    inputs.twist_rate = 11.0;
    inputs.is_twist_right = true;
    inputs.enable_aerodynamic_jump = true;
    inputs.muzzle_angle = 0.01;

    let atmo = AtmosphericConditions::default();
    let cw_mps = 4.4704_f64; // 10 mph
    let wind = WindConditions {
        speed: cw_mps,
        direction: -PI / 2.0, // from the right
    };

    let sg = compute_stability_coefficient(
        &inputs,
        (atmo.altitude, atmo.temperature, atmo.pressure, 0.0),
    );
    let length_cal = inputs.bullet_length / inputs.bullet_diameter;
    let cw_from_right_mph = -(cw_mps * (-PI / 2.0_f64).sin()) * MS_TO_MPH; // = +10 mph
    let expected = litz_crosswind_jump_moa(sg, length_cal, cw_from_right_mph, true);

    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(600.0);
    let r = solver.solve().unwrap();
    let aj = r.aerodynamic_jump.expect("AJ present when enabled");

    assert!(
        (aj.vertical_jump_moa - expected).abs() < 1e-9,
        "engine vertical jump {} != Litz {}",
        aj.vertical_jump_moa,
        expected
    );
    // Plausible Litz magnitude for a .308 at 10 mph crosswind (a few tenths MOA).
    assert!(
        (0.2..0.8).contains(&aj.vertical_jump_moa.abs()),
        "magnitude {} outside plausible Litz range",
        aj.vertical_jump_moa
    );
    // Physical anchor: a wind from the right pushes the bullet LEFT (z < 0)...
    let z = r.position_at_range(500.0).unwrap().z;
    assert!(z < 0.0, "wind from the right should drift the bullet left, z={z}");
    // ...and for a right twist that crosswind jumps the impact UP.
    assert!(aj.vertical_jump_moa > 0.0, "right twist + wind from right -> up");
}

#[test]
fn aj_direction_flips_with_wind_side_and_twist() {
    let v = |dir: f64, right_twist: bool| -> f64 {
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 790.0;
        inputs.twist_rate = 11.0;
        inputs.is_twist_right = right_twist;
        inputs.enable_aerodynamic_jump = true;
        inputs.muzzle_angle = 0.01;
        let wind = WindConditions {
            speed: 4.4704,
            direction: dir,
        };
        let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        s.set_max_range(600.0);
        s.solve()
            .unwrap()
            .aerodynamic_jump
            .unwrap()
            .vertical_jump_moa
    };
    let from_right = v(-PI / 2.0, true);
    let from_left = v(PI / 2.0, true);
    assert!(
        from_right > 0.0 && from_left < 0.0,
        "right twist: wind from right -> up, from left -> down (R={from_right}, L={from_left})"
    );
    assert!(
        (from_right + from_left).abs() < 1e-9,
        "left/right crosswind jumps should be symmetric"
    );
    // Left-hand twist reverses the jump direction.
    assert!(
        v(-PI / 2.0, false) < 0.0,
        "left twist reverses the jump direction"
    );
}

#[test]
fn zeroing_ignores_aerodynamic_jump() {
    // The zero must be found on the bare bore so AJ is an additive POI shift, never
    // absorbed by the zero search — even when zeroing with a crosswind present.
    use ballistics_engine::calculate_zero_angle_with_conditions;
    let make = |aj: bool| {
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 790.0;
        inputs.twist_rate = 11.0;
        inputs.enable_aerodynamic_jump = aj;
        inputs
    };
    let wind = WindConditions {
        speed: 4.4704, // 10 mph crosswind present during zeroing
        direction: PI / 2.0,
    };
    let atmo = AtmosphericConditions::default();
    let z_off =
        calculate_zero_angle_with_conditions(make(false), 200.0, 0.0, wind.clone(), atmo.clone())
            .unwrap();
    let z_on =
        calculate_zero_angle_with_conditions(make(true), 200.0, 0.0, wind, atmo).unwrap();
    assert!(
        (z_on - z_off).abs() < 1e-12,
        "zero angle must not depend on AJ (on={z_on}, off={z_off})"
    );
}

#[test]
fn aj_affects_only_vertical_not_windage() {
    // AJ is a vertical effect; it must NOT change the lateral wind drift. This proves it
    // doesn't double-count with the integrator's crosswind response (which is purely
    // lateral for a point-mass solve): AJ adds the missing vertical, leaving windage intact.
    let solve_xyz = |aj: bool| {
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 790.0;
        inputs.twist_rate = 11.0;
        inputs.enable_aerodynamic_jump = aj;
        inputs.muzzle_angle = 0.01;
        let wind = WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
        };
        let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        s.set_max_range(600.0);
        let p = s.solve().unwrap().position_at_range(500.0).unwrap();
        (p.y, p.z)
    };
    let (y_off, z_off) = solve_xyz(false);
    let (y_on, z_on) = solve_xyz(true);
    let dz = (z_on - z_off).abs();
    let dy = (y_on - y_off).abs();
    // The vertical shift is the real AJ effect (~0.1 m at 500 m).
    assert!(dy > 1e-2, "AJ should shift the vertical impact (dy={dy})");
    // Windage is essentially untouched: the only lateral change is the negligible
    // second-order coupling from cos(elevation) (~1e-6 m), NOT a duplicated wind drift.
    // It must be orders of magnitude smaller than the vertical shift.
    assert!(
        dz < 1e-4 && dz < 1e-3 * dy,
        "AJ must not double-count wind drift: lateral change {dz} not negligible vs vertical {dy}"
    );
}

#[test]
fn nan_twist_is_guarded_and_does_not_poison_trajectory() {
    // A NaN twist must not slip past the guard (NaN <= 0.0 is false) and NaN-out
    // the launch angle. AJ must be suppressed and the trajectory stay finite.
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 800.0;
    inputs.muzzle_angle = 0.01;
    inputs.target_distance = 500.0;
    inputs.enable_aerodynamic_jump = true;
    inputs.twist_rate = f64::NAN;
    let wind = WindConditions {
        speed: 10.0,
        direction: PI / 2.0,
    };
    let mut solver = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
    solver.set_max_range(600.0);
    let r = solver.solve().expect("solve should succeed even with NaN twist");
    assert!(
        r.aerodynamic_jump.is_none(),
        "AJ must be suppressed for a non-finite twist"
    );
    assert!(
        vertical_at_500(&r).is_finite(),
        "trajectory must stay finite when twist is NaN"
    );
}

#[test]
fn no_wind_no_tipoff_jump_is_negligible() {
    // With zero crosswind and zero tip-off yaw, the jump should be ~0, so the
    // trajectory must stay numerically indistinguishable from the disabled case.
    let off = vertical_at_500(&solve(false, 0.0));
    let on = vertical_at_500(&solve(true, 0.0));
    assert!(
        (on - off).abs() < 1e-6,
        "AJ with no crosswind/tip-off must not perturb the trajectory (off={off}, on={on})"
    );
}
