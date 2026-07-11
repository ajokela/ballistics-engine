// MBA-958: the Magnus / yaw-of-repose Sg now carries the canonical linear Miller density
// correction (T/T0)*(P0/P) = rho0/rho, matching the spin-drift Sg (MBA-942) and stability.rs.
// Magnus is off by default, so this guards that enabling it produces a finite, non-zero lateral
// effect that responds to altitude (the correction is no-op at sea-level standard, factor 1.0).
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

fn magnus_lateral(pressure_hpa: f64, temp_c: f64) -> f64 {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0; // m/s
    i.bc_value = 0.5;
    i.bc_type = DragModel::G7;
    i.bullet_mass = 168.0 * 0.00006479891; // kg
    i.bullet_diameter = 0.308 * 0.0254; // m
    i.bullet_length = 1.24 * 0.0254;
    i.caliber_inches = 0.308;
    i.weight_grains = 168.0;
    i.twist_rate = 12.0;
    i.is_twist_right = true;
    i.enable_magnus = true;
    let mut a = AtmosphericConditions::default();
    a.pressure = pressure_hpa;
    a.temperature = temp_c;
    let mut s = TrajectorySolver::new(i, WindConditions::default(), a);
    s.set_max_range(800.0);
    let r = s.solve().expect("solve");
    r.points.last().unwrap().position.z // lateral (McCoy Z)
}

#[test]
fn magnus_enabled_lateral_finite_and_altitude_responsive() {
    let sea_level = magnus_lateral(1013.25, 15.0); // standard -> density correction factor = 1.0
    let altitude = magnus_lateral(696.9, -4.8); // ~10 kft standard

    assert!(
        sea_level.is_finite() && sea_level.abs() > 1e-5,
        "sea-level Magnus lateral should be finite and non-zero, got {sea_level}"
    );
    assert!(
        altitude.is_finite(),
        "altitude Magnus lateral should be finite, got {altitude}"
    );
    // The density-corrected Sg largely offsets thinner air in the Magnus force, so check a
    // scale-independent response rather than an absolute displacement tied to the old yaw model.
    let relative_change = (sea_level - altitude).abs() / sea_level.abs();
    assert!(
        relative_change > 0.01,
        "Magnus should respond to altitude: sea_level={sea_level}, altitude={altitude}, relative_change={relative_change}"
    );
}
