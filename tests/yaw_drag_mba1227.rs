//! MBA-1227: quadratic yaw drag must be ADDITIVE per McCoy (CD = CD0 + CD_delta2 * delta^2),
//! not the old multiplicative Cd*(1 + delta^2), which implied CD_delta2 == CD0 — an order of
//! magnitude below literature (~4-20 per rad^2 for spitzer rifle bullets; 7.62 M80 ~= 9.6).
//!
//! The term lives in the RK4/RK45 generic integrator kernel (`derivatives::compute_derivatives`),
//! the only solver family that models tip-off yaw, so these tests drive the kernel directly.

use ballistics_engine::derivatives::compute_derivatives;
use ballistics_engine::BallisticInputs;
use nalgebra::Vector3;

/// Drag-axis acceleration magnitude for a straight 700 m/s shot at sea level.
fn drag_accel(tipoff_yaw_rad: f64, decay_m: f64, cd_delta2: f64, pos_x: f64) -> f64 {
    let inputs = BallisticInputs {
        tipoff_yaw: tipoff_yaw_rad,
        tipoff_decay_distance: decay_m,
        cd_delta2,
        ..BallisticInputs::default()
    };
    let d = compute_derivatives(
        Vector3::new(pos_x, 0.0, 0.0),
        Vector3::new(700.0, 0.0, 0.0),
        &inputs,
        Vector3::zeros(),
        (1.225, 340.0, 0.0, 0.0),
        inputs.bc_value,
        None,
        0.0,
        None,
    );
    d[3].abs()
}

/// With zero tip-off yaw (the default), cd_delta2 must be completely inert.
#[test]
fn zero_yaw_leaves_cd_delta2_inert() {
    assert_eq!(drag_accel(0.0, 0.0, 7.5, 0.0), drag_accel(0.0, 0.0, 20.0, 0.0));
}

/// A 2-degree constant yaw with the literature default must raise drag by a few
/// percent (the old multiplicative form produced ~0.12%, an order of magnitude low).
#[test]
fn two_degree_yaw_raises_drag_by_literature_sane_percent() {
    let base = drag_accel(0.0, 0.0, 7.5, 0.0);
    let yawed = drag_accel(2.0_f64.to_radians(), 0.0, 7.5, 0.0);
    let ratio = yawed / base - 1.0;
    assert!(
        (0.005..0.10).contains(&ratio),
        "2-deg yaw drag increase {ratio:.5} outside the literature-sane band \
         (old buggy behavior would give ~0.0012)"
    );
}

/// The additive term is exactly linear in CD_delta2: doubling the coefficient
/// must exactly double the extra drag (to FP tolerance).
#[test]
fn extra_drag_is_exactly_linear_in_cd_delta2() {
    let base = drag_accel(0.0, 0.0, 7.5, 0.0);
    let extra1 = drag_accel(2.0_f64.to_radians(), 0.0, 7.5, 0.0) - base;
    let extra2 = drag_accel(2.0_f64.to_radians(), 0.0, 15.0, 0.0) - base;
    assert!(extra1 > 0.0);
    let ratio = extra2 / extra1;
    assert!(
        (ratio - 2.0).abs() < 1e-9,
        "extra drag must be exactly linear in cd_delta2 (got ratio {ratio})"
    );
}

/// The exponential tip-off decay applies to delta, so the extra drag (~delta^2)
/// at x = 2*decay_distance must scale by exactly e^-4.
#[test]
fn decay_wiring_scales_extra_drag_by_exp_minus_four() {
    let delta = 2.0_f64.to_radians();
    let base_muzzle = drag_accel(0.0, 50.0, 7.5, 0.0);
    let base_100 = drag_accel(0.0, 50.0, 7.5, 100.0);
    let extra_muzzle = drag_accel(delta, 50.0, 7.5, 0.0) - base_muzzle;
    let extra_100 = drag_accel(delta, 50.0, 7.5, 100.0) - base_100;
    let ratio = extra_100 / extra_muzzle;
    let expected = (-4.0_f64).exp();
    assert!(
        (ratio - expected).abs() / expected < 1e-6,
        "decay scaling {ratio:.8} != e^-4 ({expected:.8})"
    );
}
