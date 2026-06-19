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
