// Magnus force from yaw of repose must be a small vertical secondary effect.
use ballistics_engine::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};

fn endpoint_m(enable_magnus: bool, enable_spin: bool) -> nalgebra::Vector3<f64> {
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
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(914.4); // 1000 yd
    let r = solver.solve().unwrap();
    r.points.last().unwrap().position
}

#[test]
fn magnus_is_small_secondary_effect() {
    let baseline = endpoint_m(false, false);
    let magnus = endpoint_m(true, false) - baseline;
    let spin = endpoint_m(false, true) - baseline;
    assert!(
        magnus.y < 0.0,
        "right-hand Magnus should add a small downward displacement, got {magnus:?}"
    );
    assert!(
        magnus.y.abs() < 0.05,
        "Magnus displacement should be < ~2in (0.05m) at 1000yd, got {magnus:?}"
    );
    assert!(
        magnus.z.abs() < 1e-12,
        "Magnus must not masquerade as lateral spin drift, got {magnus:?}"
    );
    assert!(
        spin.z.abs() > magnus.z.abs(),
        "the separate Litz model must own lateral spin drift, got {spin:?}"
    );
}
