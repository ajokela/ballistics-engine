// MBA-958: the Magnus / yaw-of-repose Sg now carries the canonical linear Miller density
// correction (T/T0)*(P0/P) = rho0/rho, matching the spin-drift Sg (MBA-942) and stability.rs.
// Magnus is off by default, so this guards that enabling it produces a finite, non-zero vertical
// effect that responds to altitude (the correction is no-op at sea-level standard, factor 1.0).
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

fn magnus_vertical_delta(pressure_hpa: f64, temp_c: f64) -> f64 {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 800.0; // m/s
    inputs.bc_value = 0.5;
    inputs.bc_type = DragModel::G7;
    inputs.bullet_mass = 168.0 * 0.00006479891; // kg
    inputs.bullet_diameter = 0.308 * 0.0254; // m
    inputs.bullet_length = 1.24 * 0.0254;
    inputs.caliber_inches = 0.308;
    inputs.weight_grains = 168.0;
    inputs.twist_rate = 12.0;
    inputs.is_twist_right = true;
    let mut a = AtmosphericConditions::default();
    a.pressure = pressure_hpa;
    a.temperature = temp_c;
    let endpoint_y = |enable_magnus| {
        let mut run_inputs = inputs.clone();
        run_inputs.enable_magnus = enable_magnus;
        let mut solver = TrajectorySolver::new(run_inputs, WindConditions::default(), a.clone());
        solver.set_max_range(800.0);
        solver
            .solve()
            .expect("solve")
            .points
            .last()
            .unwrap()
            .position
            .y
    };

    endpoint_y(true) - endpoint_y(false)
}

#[test]
fn magnus_enabled_vertical_finite_and_altitude_responsive() {
    let sea_level = magnus_vertical_delta(1013.25, 15.0); // standard -> correction factor = 1.0
    let altitude = magnus_vertical_delta(696.9, -4.8); // ~10 kft standard

    assert!(
        sea_level.is_finite() && sea_level < -1e-5,
        "sea-level right-hand Magnus vertical delta should be finite and down, got {sea_level}"
    );
    assert!(
        altitude.is_finite(),
        "altitude Magnus vertical delta should be finite, got {altitude}"
    );
    // The density-corrected Sg largely offsets thinner air in the Magnus force, so check a
    // scale-independent response rather than an absolute displacement tied to the old yaw model.
    let relative_change = (sea_level - altitude).abs() / sea_level.abs();
    assert!(
        relative_change > 0.01,
        "Magnus should respond to altitude: sea_level={sea_level}, altitude={altitude}, relative_change={relative_change}"
    );
}
