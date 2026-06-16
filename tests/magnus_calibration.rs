// Magnus force must be a small secondary effect: lateral drift far below the
// gyroscopic spin drift, and below ~2 inches at 1000 yd.
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions,
};

fn lateral_drift_m(enable_magnus: bool, enable_spin: bool) -> f64 {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 823.0; // m/s (~2700 fps)
    inputs.bullet_mass = 168.0 * 0.00006479891; // kg
    inputs.bullet_diameter = 0.308 * 0.0254; // m
    inputs.bullet_length = 1.215 * 0.0254; // m
    inputs.bc_value = 0.475;
    inputs.twist_rate = 12.0;
    inputs.is_twist_right = true;
    inputs.enable_magnus = enable_magnus;
    inputs.use_enhanced_spin_drift = enable_spin;
    // McCoy frame: Z is lateral
    let mut solver =
        TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    solver.set_max_range(914.4); // 1000 yd
    let r = solver.solve().unwrap();
    r.points.last().unwrap().position.z
}

#[test]
fn magnus_is_small_secondary_effect() {
    let magnus = lateral_drift_m(true, false).abs();
    let spin = lateral_drift_m(false, true).abs();
    assert!(magnus > 0.0, "Magnus should produce some drift, got {magnus}");
    assert!(
        magnus < 0.05,
        "Magnus drift should be < ~2in (0.05m) at 1000yd, got {magnus} m"
    );
    assert!(
        magnus < spin,
        "Magnus drift ({magnus} m) should be far below spin drift ({spin} m)"
    );
}
