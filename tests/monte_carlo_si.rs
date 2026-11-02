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
    BallisticInputs {
        target_distance: 457.2,             // meters (500 yd)
        muzzle_velocity: 823.0,             // m/s
        bullet_mass: 168.0 * 0.00006479891, // kg
        bullet_diameter: 0.308 * 0.0254,    // meters
        bullet_length: 1.215 * 0.0254,
        caliber_inches: 0.308,
        weight_grains: 168.0,
        bc_value: 0.475,
        bc_type: DragModel::G1,
        muzzle_angle: 0.006, // radians
        sight_height: 0.05,  // meters
        muzzle_height: 0.0,
        temperature: 15.0, // Celsius
        pressure: 1013.25, // hPa
        ..BallisticInputs::default()
    }
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

    // The two independent solvers fly the same air (the station at the default 50 %), so only
    // their integrators differ. Measured on this shot (0.50.0's R5): velocity 2.1e-6 and energy
    // 4.3e-6 relative, drop 1.15e-4 m. The bounds below are about 10x, 10x and 2.6x those. The
    // lean Monte Carlo of 0.49.0 flew the station's sound speed dry and missed by 2.3e-4,
    // 4.7e-4 and 4.5e-4 m, so each bound catches it; the 5 % (and 5 mm) this test allowed
    // before passed it by two orders of magnitude.
    let dvel = (mc.velocity - last.velocity_magnitude).abs() / last.velocity_magnitude;
    let den = (mc.energy - last.kinetic_energy).abs() / last.kinetic_energy;
    assert!(
        dvel < 2e-5,
        "velocity disagreement {dvel:e} relative (mc={}, cli={})",
        mc.velocity,
        last.velocity_magnitude
    );
    assert!(
        den < 4e-5,
        "energy disagreement {den:e} relative (mc={}, cli={})",
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
        (mc.drop - cli_drop).abs() < 3e-4,
        "drop-frame disagreement: mc={:.6}m cli={:.6}m",
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
