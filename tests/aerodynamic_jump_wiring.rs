//! MBA-959: verify aerodynamic jump is wired into the solver as an opt-in,
//! default-off muzzle launch-angle perturbation.

use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions,
};
use std::f64::consts::PI;

/// Build a solver for a representative .308-class load. `crosswind_mps` is applied
/// as a pure lateral (McCoy Z) wind so it drives aerodynamic jump.
fn solve(enable_aj: bool, crosswind_mps: f64) -> ballistics_engine::TrajectoryResult {
    let inputs = BallisticInputs {
        muzzle_velocity: 800.0, // m/s
        muzzle_angle: 0.01, // small positive elevation so it carries downrange
        target_distance: 500.0,
        twist_rate: 11.0, // 1:11"
        is_twist_right: true,
        enable_aerodynamic_jump: enable_aj,
        ..BallisticInputs::default()
    };

    // direction = PI/2 -> wind is purely lateral (crosswind), no head/tail component.
    let wind = WindConditions {
        speed: crosswind_mps,
        direction: PI / 2.0,
        vertical_speed: 0.0,
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

const MS_TO_MPH: f64 = 2.236_936_292_054_4;

#[test]
fn litz_magnitude_matches_engine_sg_exactly() {
    // The engine's vertical jump must equal Litz's Y = (0.01*Sg - 0.0024*L + 0.032)
    // * crosswind_mph, computed from the engine's OWN Miller Sg.
    use ballistics_engine::aerodynamic_jump::litz_crosswind_jump_moa;
    use ballistics_engine::stability::compute_stability_coefficient;

    let inputs = BallisticInputs {
        muzzle_velocity: 790.0,
        bullet_diameter: 0.00782, // ~.308"
        bullet_length: 0.0312, // ~4.0 cal
        bullet_mass: 0.01134, // ~175 gr
        twist_rate: 11.0,
        is_twist_right: true,
        enable_aerodynamic_jump: true,
        muzzle_angle: 0.01,
        ..BallisticInputs::default()
    };

    let atmo = AtmosphericConditions::default();
    let cw_mps = 4.4704_f64; // 10 mph
    let wind = WindConditions {
        speed: cw_mps,
        direction: PI / 2.0, // from the right (0=headwind, +90deg=from right)
        vertical_speed: 0.0,
    };

    let sg = compute_stability_coefficient(
        &inputs,
        (atmo.altitude, atmo.temperature, atmo.pressure, 0.0),
    );
    let length_cal = inputs.bullet_length / inputs.bullet_diameter;
    let cw_from_right_mph = cw_mps * (PI / 2.0_f64).sin() * MS_TO_MPH; // = +10 mph
    let expected = litz_crosswind_jump_moa(sg, length_cal, cw_from_right_mph, true);

    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(600.0);
    let r = solver.solve().unwrap();
    let aj = r.aerodynamic_jump.expect("AJ present when enabled");

    assert!(
        (aj.vertical_jump_moa - expected).abs() < 1e-9,
        "engine vertical jump {} != Litz {}",
        aj.vertical_jump_moa,
        expected
    );
    // Plausible Litz magnitude for a .308 at 10 mph crosswind (a few tenths MOA).
    assert!(
        (0.2..0.8).contains(&aj.vertical_jump_moa.abs()),
        "magnitude {} outside plausible Litz range",
        aj.vertical_jump_moa
    );
    // Physical anchor: a wind from the right pushes the bullet LEFT (z < 0)...
    let z = r.position_at_range(500.0).unwrap().z;
    assert!(z < 0.0, "wind from the right should drift the bullet left, z={z}");
    // ...and for a right twist that crosswind jumps the impact UP.
    assert!(aj.vertical_jump_moa > 0.0, "right twist + wind from right -> up");
}

#[test]
fn segmented_muzzle_wind_matches_scalar_aerodynamic_jump() {
    let inputs = BallisticInputs {
        muzzle_velocity: 790.0,
        twist_rate: 11.0,
        is_twist_right: true,
        enable_aerodynamic_jump: true,
        muzzle_angle: 0.01,
        ..BallisticInputs::default()
    };
    let atmosphere = AtmosphericConditions::default();

    let mut scalar_solver = TrajectorySolver::new(
        inputs.clone(),
        WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
            vertical_speed: 0.0,
        },
        atmosphere.clone(),
    );
    scalar_solver.set_max_range(600.0);
    let scalar_jump = scalar_solver
        .solve()
        .unwrap()
        .aerodynamic_jump
        .expect("scalar crosswind must produce aerodynamic jump")
        .vertical_jump_moa;

    // A conflicting 20 mph scalar wind from the left must be ignored while segments are present.
    let mut segmented_solver = TrajectorySolver::new(
        inputs,
        WindConditions {
            speed: 8.9408,
            direction: -PI / 2.0,
            vertical_speed: 0.0,
        },
        atmosphere,
    );
    segmented_solver.set_max_range(600.0);
    segmented_solver.set_wind_segments(vec![
        ballistics_engine::wind::WindSegment::new(16.09344, 90.0, 100.0),
        ballistics_engine::wind::WindSegment::new(32.18688, 270.0, 5000.0),
    ]);
    let segmented_jump = segmented_solver
        .solve()
        .unwrap()
        .aerodynamic_jump
        .expect("segmented muzzle crosswind must produce aerodynamic jump")
        .vertical_jump_moa;

    assert!(
        segmented_jump > 0.0,
        "muzzle wind from the right must jump up"
    );
    assert!(
        (segmented_jump - scalar_jump).abs() < 1e-9,
        "equivalent scalar and segmented muzzle wind must produce the same jump: scalar={scalar_jump}, segmented={segmented_jump}"
    );
}

#[test]
fn aj_direction_flips_with_wind_side_and_twist() {
    let v = |dir: f64, right_twist: bool| -> f64 {
        let inputs = BallisticInputs {
            muzzle_velocity: 790.0,
            twist_rate: 11.0,
            is_twist_right: right_twist,
            enable_aerodynamic_jump: true,
            muzzle_angle: 0.01,
            ..BallisticInputs::default()
        };
        let wind = WindConditions {
            speed: 4.4704,
            direction: dir,
            vertical_speed: 0.0,
        };
        let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        s.set_max_range(600.0);
        s.solve()
            .unwrap()
            .aerodynamic_jump
            .unwrap()
            .vertical_jump_moa
    };
    // Corrected wind-FROM convention: +90deg (+PI/2) = from the right, -90deg = from the left.
    let from_right = v(PI / 2.0, true);
    let from_left = v(-PI / 2.0, true);
    assert!(
        from_right > 0.0 && from_left < 0.0,
        "right twist: wind from right -> up, from left -> down (R={from_right}, L={from_left})"
    );
    assert!(
        (from_right + from_left).abs() < 1e-9,
        "left/right crosswind jumps should be symmetric"
    );
    // Left-hand twist reverses the jump direction.
    assert!(
        v(PI / 2.0, false) < 0.0,
        "left twist reverses the jump direction"
    );
}

#[test]
fn zeroing_ignores_aerodynamic_jump() {
    // The zero must be found on the bare bore so AJ is an additive POI shift, never
    // absorbed by the zero search — even when zeroing with a crosswind present.
    use ballistics_engine::calculate_zero_angle_with_conditions;
    let make = |aj: bool| BallisticInputs {
        muzzle_velocity: 790.0,
        twist_rate: 11.0,
        enable_aerodynamic_jump: aj,
        ..BallisticInputs::default()
    };
    let wind = WindConditions {
        speed: 4.4704, // 10 mph crosswind present during zeroing
        direction: PI / 2.0,
        vertical_speed: 0.0,
    };
    let atmo = AtmosphericConditions::default();
    let z_off =
        calculate_zero_angle_with_conditions(make(false), 200.0, 0.0, wind.clone(), atmo.clone())
            .unwrap();
    let z_on =
        calculate_zero_angle_with_conditions(make(true), 200.0, 0.0, wind, atmo).unwrap();
    assert!(
        (z_on - z_off).abs() < 1e-12,
        "zero angle must not depend on AJ (on={z_on}, off={z_off})"
    );
}

#[test]
fn aj_affects_only_vertical_not_windage() {
    // AJ is a vertical effect; it must NOT change the lateral wind drift. This proves it
    // doesn't double-count with the integrator's crosswind response (which is purely
    // lateral for a point-mass solve): AJ adds the missing vertical, leaving windage intact.
    let solve_xyz = |aj: bool| {
        let inputs = BallisticInputs {
            muzzle_velocity: 790.0,
            twist_rate: 11.0,
            enable_aerodynamic_jump: aj,
            muzzle_angle: 0.01,
            ..BallisticInputs::default()
        };
        let wind = WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
            vertical_speed: 0.0,
        };
        let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
        s.set_max_range(600.0);
        let p = s.solve().unwrap().position_at_range(500.0).unwrap();
        (p.y, p.z)
    };
    let (y_off, z_off) = solve_xyz(false);
    let (y_on, z_on) = solve_xyz(true);
    let dz = (z_on - z_off).abs();
    let dy = (y_on - y_off).abs();
    // The vertical shift is the real AJ effect (~0.1 m at 500 m).
    assert!(dy > 1e-2, "AJ should shift the vertical impact (dy={dy})");
    // Windage is essentially untouched: the only lateral change is the negligible
    // second-order coupling from cos(elevation) (~1e-6 m), NOT a duplicated wind drift.
    // It must be orders of magnitude smaller than the vertical shift.
    assert!(
        dz < 1e-4 && dz < 1e-3 * dy,
        "AJ must not double-count wind drift: lateral change {dz} not negligible vs vertical {dy}"
    );
}

#[test]
fn fast_path_launch_offset_sign_disabled_and_flips() {
    // MBA-959: the engine fast-integrate path (used by the binding's fast_integrate + MC).
    use ballistics_engine::fast_trajectory::aerodynamic_jump_launch_offset_rad;
    let mut inputs = BallisticInputs {
        muzzle_velocity: 800.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.0312,
        bullet_mass: 0.01134,
        twist_rate: 11.0,
        is_twist_right: true,
        wind_speed: 4.4704, // m/s (~10 mph)
        wind_angle: PI / 2.0, // BallisticInputs convention: 90deg = from the right
        ..BallisticInputs::default()
    };
    let atmo = (0.0, 15.0, 1013.25, 1.0);

    inputs.enable_aerodynamic_jump = false;
    assert_eq!(aerodynamic_jump_launch_offset_rad(&inputs, atmo), 0.0);

    inputs.enable_aerodynamic_jump = true;
    let from_right = aerodynamic_jump_launch_offset_rad(&inputs, atmo);
    assert!(
        from_right > 0.0,
        "right twist + wind from the right -> up (positive), got {from_right}"
    );
    // Reverse the wind side -> opposite sign, same magnitude.
    inputs.wind_angle = -PI / 2.0;
    assert!((aerodynamic_jump_launch_offset_rad(&inputs, atmo) + from_right).abs() < 1e-12);
    // Flip the twist hand -> opposite sign.
    inputs.wind_angle = PI / 2.0;
    inputs.is_twist_right = false;
    assert!((aerodynamic_jump_launch_offset_rad(&inputs, atmo) + from_right).abs() < 1e-12);
}

#[test]
fn fast_path_direct_atmosphere_uses_density_for_aerodynamic_jump() {
    use ballistics_engine::fast_trajectory::aerodynamic_jump_launch_offset_rad;

    let inputs = BallisticInputs {
        muzzle_velocity: 800.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.0312,
        bullet_mass: 0.01134,
        twist_rate: 11.0,
        is_twist_right: true,
        wind_speed: 4.4704,
        wind_angle: PI / 2.0,
        enable_aerodynamic_jump: true,
        ..BallisticInputs::default()
    };

    let standard = aerodynamic_jump_launch_offset_rad(&inputs, (0.0, 15.0, 1013.25, 1.0));
    let direct_standard_density =
        aerodynamic_jump_launch_offset_rad(&inputs, (1.225, 340.0, 0.0, 0.0));
    assert!(standard > 0.0);
    assert!(
        (direct_standard_density - standard).abs() < 1e-12,
        "direct standard density must preserve AJ: standard={standard}, direct={direct_standard_density}"
    );

    let thin_density = 0.9;
    let equivalent_pressure = 1013.25 * thin_density / 1.225;
    let direct_thin = aerodynamic_jump_launch_offset_rad(&inputs, (thin_density, 340.0, 0.0, 0.0));
    let standard_thin =
        aerodynamic_jump_launch_offset_rad(&inputs, (0.0, 15.0, equivalent_pressure, 1.0));
    assert!(
        (direct_thin - standard_thin).abs() < 1e-12,
        "direct density correction must match equivalent standard atmosphere: direct={direct_thin}, standard={standard_thin}"
    );
}

#[test]
fn fast_paths_apply_aerodynamic_jump() {
    use ballistics_engine::fast_trajectory::{
        fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
    };
    use ballistics_engine::wind::WindSock;

    fn y_at_500(enable_aj: bool, with_segments: bool, atmo_params: (f64, f64, f64, f64)) -> f64 {
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            bullet_diameter: 0.00782,
            bullet_length: 0.0312,
            bullet_mass: 0.01134,
            twist_rate: 11.0,
            is_twist_right: true,
            wind_speed: 4.4704,
            wind_angle: PI / 2.0, // from the right
            enable_aerodynamic_jump: enable_aj,
            ..BallisticInputs::default()
        };
        let v = inputs.muzzle_velocity;
        let elev = 0.02_f64;
        let initial_state = [0.0, 0.0, 0.0, v * elev.cos(), v * elev.sin(), 0.0];
        let params = FastIntegrationParams {
            horiz: 600.0,
            vert: 0.0,
            initial_state,
            t_span: (0.0, 3.0),
            atmo_params,
            atmo_sock: None,
        };
        let sol = if with_segments {
            fast_integrate_with_segments(&inputs, vec![], params)
        } else {
            fast_integrate(&inputs, &WindSock::new(vec![]), params)
        };
        // FastSolution.y is column-major [6][n_points]: y[0]=x series, y[1]=vertical series.
        let xs = &sol.y[0];
        let ys = &sol.y[1];
        let target = 457.2; // ~500 yd downrange
        let mut out = *ys.last().unwrap();
        for i in 1..xs.len() {
            if xs[i] >= target {
                let f = (target - xs[i - 1]) / (xs[i] - xs[i - 1]);
                out = ys[i - 1] + f * (ys[i] - ys[i - 1]);
                break;
            }
        }
        out
    }

    for (atmo_name, atmo_params) in [
        ("standard", (0.0, 15.0, 1013.25, 1.0)),
        ("direct", (1.225, 340.0, 0.0, 0.0)),
    ] {
        for (path, with_segments) in [("plain", false), ("segmented", true)] {
            let off = y_at_500(false, with_segments, atmo_params);
            let on = y_at_500(true, with_segments, atmo_params);
            assert!(
                (on - off).abs() > 1e-3,
                "AJ should shift the {path} fast-path vertical position in {atmo_name} mode (off={off}, on={on})"
            );
            assert!(
                on > off,
                "right twist + wind from the right should raise the {path} impact in {atmo_name} mode (off={off}, on={on})"
            );
        }
    }
}

#[test]
fn segmented_fast_path_records_the_same_aj_launch_state_as_plain_fast() {
    use ballistics_engine::fast_trajectory::{
        fast_integrate, fast_integrate_with_segments, FastIntegrationParams,
    };
    use ballistics_engine::wind::WindSock;

    let mut inputs = BallisticInputs {
        muzzle_velocity: 800.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.0312,
        bullet_mass: 0.01134,
        twist_rate: 11.0,
        is_twist_right: true,
        wind_speed: 4.4704,
        wind_angle: PI / 2.0,
        enable_aerodynamic_jump: true,
        ..BallisticInputs::default()
    };
    let raw_state = [0.0, 0.0, 0.0, 800.0, 0.0, 0.0];
    let params = || FastIntegrationParams {
        horiz: 20.0,
        vert: 0.0,
        initial_state: raw_state,
        t_span: (0.0, 0.1),
        atmo_params: (0.0, 15.0, 1013.25, 1.0),
        atmo_sock: None,
    };

    let plain = fast_integrate(&inputs, &WindSock::new(vec![]), params());
    let segmented = fast_integrate_with_segments(&inputs, vec![], params());
    for component in 3..6 {
        assert!(
            (plain.y[component][0] - segmented.y[component][0]).abs() < 1e-12,
            "AJ launch component {component} differs: plain={} segmented={}",
            plain.y[component][0],
            segmented.y[component][0]
        );
    }
    assert!(
        (plain.y[4][0] - raw_state[4]).abs() > 1e-9,
        "enabled AJ must measurably rotate the launch elevation"
    );

    inputs.enable_aerodynamic_jump = false;
    let disabled = fast_integrate_with_segments(&inputs, vec![], params());
    for (component, expected) in disabled.y[3..6].iter().zip(&raw_state[3..6]) {
        assert_eq!(component[0].to_bits(), expected.to_bits());
    }
}

#[test]
fn nan_twist_is_guarded_and_does_not_poison_trajectory() {
    // A NaN twist must not slip past the guard (NaN <= 0.0 is false) and NaN-out
    // the launch angle. AJ must be suppressed and the trajectory stay finite.
    let inputs = BallisticInputs {
        muzzle_velocity: 800.0,
        muzzle_angle: 0.01,
        target_distance: 500.0,
        enable_aerodynamic_jump: true,
        twist_rate: f64::NAN,
        ..BallisticInputs::default()
    };
    let wind = WindConditions {
        speed: 10.0,
        direction: PI / 2.0,
        vertical_speed: 0.0,
    };
    let mut solver = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
    solver.set_max_range(600.0);
    let r = solver.solve().expect("solve should succeed even with NaN twist");
    assert!(
        r.aerodynamic_jump.is_none(),
        "AJ must be suppressed for a non-finite twist"
    );
    assert!(
        vertical_at_500(&r).is_finite(),
        "trajectory must stay finite when twist is NaN"
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
