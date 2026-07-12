// Magnus force from yaw of repose must be a small vertical secondary effect.
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};

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

fn fast_endpoint_m(
    enable_magnus: bool,
    use_litz: bool,
    is_twist_right: bool,
    with_segments: bool,
) -> nalgebra::Vector3<f64> {
    use ballistics_engine::fast_trajectory::{
        fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
    };
    use ballistics_engine::wind::WindSock;

    let inputs = BallisticInputs {
        muzzle_velocity: 823.0,
        bullet_mass: 168.0 * 0.00006479891,
        bullet_diameter: 0.308 * 0.0254,
        bullet_length: 1.215 * 0.0254,
        caliber_inches: 0.308,
        weight_grains: 168.0,
        bc_value: 0.475,
        bc_type: DragModel::G1,
        twist_rate: 12.0,
        is_twist_right,
        enable_magnus,
        use_enhanced_spin_drift: use_litz,
        enable_advanced_effects: false,
        ground_threshold: -100.0,
        ..BallisticInputs::default()
    };
    let elevation = 0.02_f64;
    let params = FastIntegrationParams {
        horiz: 1_000.0,
        vert: 0.0,
        initial_state: [
            0.0,
            0.0,
            0.0,
            inputs.muzzle_velocity * elevation.cos(),
            inputs.muzzle_velocity * elevation.sin(),
            0.0,
        ],
        t_span: (0.0, 5.0),
        atmo_params: (0.0, 15.0, 1013.25, 1.0),
        atmo_sock: None,
    };
    let solution = if with_segments {
        fast_integrate_with_segments(&inputs, vec![], params)
    } else {
        fast_integrate(&inputs, &WindSock::new(vec![]), params)
    };
    let last = solution.t.len() - 1;
    nalgebra::Vector3::new(
        solution.y[0][last],
        solution.y[1][last],
        solution.y[2][last],
    )
}

#[test]
fn plain_fast_magnus_matches_segmented_twist_response() {
    for is_twist_right in [true, false] {
        let plain_off = fast_endpoint_m(false, false, is_twist_right, false);
        let plain_on = fast_endpoint_m(true, false, is_twist_right, false);
        let segmented_off = fast_endpoint_m(false, false, is_twist_right, true);
        let segmented_on = fast_endpoint_m(true, false, is_twist_right, true);
        let plain_delta = plain_on - plain_off;
        let segmented_delta = segmented_on - segmented_off;

        if is_twist_right {
            assert!(
                plain_delta.y < -1e-5,
                "right twist must move down: {plain_delta:?}"
            );
        } else {
            assert!(
                plain_delta.y > 1e-5,
                "left twist must move up: {plain_delta:?}"
            );
        }
        assert!(
            plain_delta.z.abs() < 1e-12,
            "Magnus must remain vertical: {plain_delta:?}"
        );
        assert!(
            (plain_delta.y - segmented_delta.y).abs() < 5e-6,
            "plain/segmented Magnus mismatch: plain={plain_delta:?}, segmented={segmented_delta:?}"
        );
    }

    let litz_off = fast_endpoint_m(false, true, true, false);
    let litz_with_magnus = fast_endpoint_m(true, true, true, false);
    assert!(
        (litz_with_magnus - litz_off).norm() < 1e-12,
        "Litz mode must suppress explicit Magnus"
    );
}
