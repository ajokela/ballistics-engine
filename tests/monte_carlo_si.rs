// Regression guard for the monte_carlo SI fix: solve_trajectory_for_monte_carlo
// now takes SI-canonical BallisticInputs (meters, m/s, kg, radians) and must
// agree with the cli_api TrajectorySolver on the same physical shot. Guards the
// SI unit conversions, the atmo_params ordering, and the base_ratio derivation
// (any of which, if regressed, makes the trajectory diverge or collapse).
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

fn si_bullet() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.target_distance = 457.2; // meters (500 yd)
    i.muzzle_velocity = 823.0; // m/s
    i.bullet_mass = 168.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // meters
    i.bullet_length = 1.215 * 0.0254;
    i.caliber_inches = 0.308;
    i.weight_grains = 168.0;
    i.bc_value = 0.475;
    i.bc_type = DragModel::G1;
    i.muzzle_angle = 0.006; // radians
    i.sight_height = 0.05; // meters
    i.muzzle_height = 0.0;
    i.temperature = 15.0; // Celsius
    i.pressure = 1013.25; // hPa
    i
}

#[test]
fn monte_carlo_agrees_with_cli_api_on_si_inputs() {
    let i = si_bullet();
    let mc = solve_trajectory_for_monte_carlo(&i).expect("monte carlo solve");

    let mut solver = TrajectorySolver::new(
        i.clone(),
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(457.2);
    let r = solver.solve().expect("cli_api solve");
    let last = r.points.last().unwrap();

    // monte_carlo must reach the target (not collapse at ~33 m) with physical energy.
    assert!(
        mc.distance > 450.0,
        "monte_carlo should reach ~457 m, got {} m",
        mc.distance
    );
    assert!(
        mc.energy > 1000.0 && mc.energy < 4000.0,
        "monte_carlo energy should be physical (~1700 J), got {} J",
        mc.energy
    );

    // The two independent solvers must agree on impact velocity and energy.
    let dvel = (mc.velocity - last.velocity_magnitude).abs() / last.velocity_magnitude;
    let den = (mc.energy - last.kinetic_energy).abs() / last.kinetic_energy;
    assert!(
        dvel < 0.05,
        "velocity disagreement {:.1}% (mc={:.1}, cli={:.1})",
        dvel * 100.0,
        mc.velocity,
        last.velocity_magnitude
    );
    assert!(
        den < 0.05,
        "energy disagreement {:.1}% (mc={:.0}, cli={:.0})",
        den * 100.0,
        mc.energy,
        last.kinetic_energy
    );
}
