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

    // monte_carlo must terminate on the requested target plane (not collapse at ~33 m or retain
    // the first fixed-step sample beyond the target) with physical energy.
    assert!(
        (mc.distance - i.target_distance).abs() < 1e-10,
        "monte_carlo should stop at {} m, got {} m",
        i.target_distance,
        mc.distance,
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

    // Both paths use the canonical horizontal line of sight at
    // muzzle_height + sight_height. Positive drop is below that line.
    let cli_at_target = r
        .position_at_range(i.target_distance)
        .expect("CLI trajectory reaches target distance");
    let cli_drop = i.muzzle_height + i.sight_height - cli_at_target.y;
    assert!(
        (mc.drop - cli_drop).abs() < 0.005,
        "drop-frame disagreement: mc={:.4}m cli={:.4}m",
        mc.drop,
        cli_drop
    );
}

#[test]
fn sight_height_offsets_line_of_sight_not_launch_position() {
    let mut bore_sight = si_bullet();
    bore_sight.muzzle_height = 1.5;
    bore_sight.sight_height = 0.0;
    let bore_drop = solve_trajectory_for_monte_carlo(&bore_sight)
        .expect("zero-sight-height Monte Carlo solve")
        .drop;

    let mut raised_sight = bore_sight;
    raised_sight.sight_height = 0.05;
    let raised_drop = solve_trajectory_for_monte_carlo(&raised_sight)
        .expect("raised-sight Monte Carlo solve")
        .drop;

    assert!(
        (raised_drop - bore_drop - 0.05).abs() < 1e-6,
        "raising the LOS by 5 cm must add 5 cm of reported drop: bore={bore_drop:.6}, \
         raised={raised_drop:.6}"
    );
}
