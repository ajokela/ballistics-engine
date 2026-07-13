// MBA-728 analytic acceptance gate for vertical wind. Drag is direction-agnostic, so a pure
// updraft `w` must deflect the shot VERTICALLY by the same amount (within 5%) that an equal
// pure crosswind `w` deflects it LATERALLY, at the same ranges.
use ballistics_engine::fast_trajectory::{fast_integrate, FastIntegrationParams};
use ballistics_engine::wind::{WindSegment, WindSock};
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectoryResult, TrajectorySolver,
    WindConditions,
};
use std::f64::consts::FRAC_PI_2;

fn base() -> BallisticInputs {
    let mut i = BallisticInputs::default();
    i.muzzle_velocity = 800.0;
    i.bc_value = 0.5;
    i.bc_type = DragModel::G7;
    i.bullet_mass = 0.0109;
    i.bullet_diameter = 0.00782;
    i.bullet_length = 0.0309;
    i.sight_height = 0.05;
    i.use_rk4 = true;
    i
}

fn solve(inputs: BallisticInputs, wind: WindConditions, max_range: f64) -> TrajectoryResult {
    let mut s = TrajectorySolver::new(inputs, wind, AtmosphericConditions::default());
    s.set_max_range(max_range);
    s.set_time_step(0.001);
    s.solve().expect("solve")
}

fn solve_segments(
    inputs: BallisticInputs,
    segments: Vec<WindSegment>,
    max_range: f64,
) -> TrajectoryResult {
    let mut s = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    s.set_max_range(max_range);
    s.set_time_step(0.001);
    s.set_wind_segments(segments);
    s.solve().expect("solve")
}

/// (y, z) at range x, linearly interpolated between the bracketing points.
fn yz_at(r: &TrajectoryResult, x: f64) -> (f64, f64) {
    let pts = &r.points;
    for i in 1..pts.len() {
        if pts[i].position.x >= x {
            let (p1, p2) = (&pts[i - 1], &pts[i]);
            let dx = p2.position.x - p1.position.x;
            let t = if dx.abs() < 1e-12 { 0.0 } else { (x - p1.position.x) / dx };
            return (
                p1.position.y + t * (p2.position.y - p1.position.y),
                p1.position.z + t * (p2.position.z - p1.position.z),
            );
        }
    }
    panic!("never reached {x}");
}

#[test]
fn vertical_and_lateral_wind_deflect_symmetrically() {
    let calm = solve(base(), WindConditions::default(), 700.0);
    let crosswind = solve(
        base(),
        WindConditions {
            speed: 5.0,
            direction: FRAC_PI_2,
            ..Default::default()
        },
        700.0,
    );
    let updraft = solve(
        base(),
        WindConditions {
            vertical_speed: 5.0,
            ..Default::default()
        },
        700.0,
    );

    for range in [300.0, 600.0] {
        let (y_calm, z_calm) = yz_at(&calm, range);
        let (y_updraft, _) = yz_at(&updraft, range);
        let (_, z_cross) = yz_at(&crosswind, range);

        let dy = y_updraft - y_calm;
        let dz = z_cross - z_calm;

        // Locked contract: an updraft raises POI.
        assert!(dy > 0.0, "updraft must raise POI at {range} m, got dy={dy}");
        // Measured (not assumed) sign of the from-right crosswind's lateral deflection;
        // this documents the existing convention rather than asserting an assumption.
        assert!(
            dz < 0.0,
            "from-right crosswind should drift the shot left (z<0) at {range} m, got dz={dz}"
        );

        let err = (dy.abs() - dz.abs()).abs();
        let bound = 0.05 * dz.abs();
        println!(
            "range={range} m: dy={dy:.6} dz={dz:.6} err={err:.6} ({:.3}% of |dz|)",
            100.0 * err / dz.abs()
        );
        assert!(
            err <= bound,
            "drag-symmetry violated at {range} m: |dy|={:.6}, |dz|={:.6}, err={:.6} > 5% bound {:.6}",
            dy.abs(),
            dz.abs(),
            err,
            bound
        );
    }
}

#[test]
fn segment_boundary_delays_vertical_deflection() {
    let calm = solve(base(), WindConditions::default(), 700.0);
    // Vertical wind only turns on at 300 m: [0-300) vertical 0, [300, 1000) vertical 5.
    let segments = vec![
        WindSegment {
            speed_kmh: 0.0,
            angle_deg: 0.0,
            until_m: 300.0,
            vertical_mps: 0.0,
        },
        WindSegment {
            speed_kmh: 0.0,
            angle_deg: 0.0,
            until_m: 1000.0,
            vertical_mps: 5.0,
        },
    ];
    let segmented = solve_segments(base(), segments, 700.0);
    let updraft_full = solve(
        base(),
        WindConditions {
            vertical_speed: 5.0,
            ..Default::default()
        },
        700.0,
    );

    let (y_calm_200, _) = yz_at(&calm, 200.0);
    let (y_seg_200, _) = yz_at(&segmented, 200.0);
    let dy_200 = y_seg_200 - y_calm_200;
    assert!(
        dy_200.abs() < 0.01,
        "no vertical wind should have applied before the 300 m boundary: dy(200)={dy_200}"
    );

    let (y_calm_600, _) = yz_at(&calm, 600.0);
    let (y_seg_600, _) = yz_at(&segmented, 600.0);
    let (y_full_600, _) = yz_at(&updraft_full, 600.0);
    let dy_seg_600 = y_seg_600 - y_calm_600;
    let dy_full_600 = y_full_600 - y_calm_600;
    println!("segment dy(600)={dy_seg_600:.6}, full-range updraft dy(600)={dy_full_600:.6}");
    assert!(
        dy_seg_600 > 0.0,
        "vertical wind after the 300 m boundary should raise POI: dy(600)={dy_seg_600}"
    );
    assert!(
        dy_seg_600 < dy_full_600,
        "partial-coverage updraft deflection {dy_seg_600} should be less than the full-range updraft {dy_full_600}"
    );
}

#[test]
fn updraft_still_raises_poi_on_an_incline() {
    let mut incline = base();
    incline.shooting_angle = 15f64.to_radians();

    let calm_incline = solve(incline.clone(), WindConditions::default(), 700.0);
    let updraft_incline = solve(
        incline,
        WindConditions {
            vertical_speed: 5.0,
            ..Default::default()
        },
        700.0,
    );

    let (y_calm, _) = yz_at(&calm_incline, 400.0);
    let (y_wind, _) = yz_at(&updraft_incline, 400.0);
    assert!(
        y_wind > y_calm,
        "updraft should still raise POI on a 15deg incline at 400 m: calm={y_calm} wind={y_wind}"
    );
}

#[test]
fn wind_conditions_default_has_zero_vertical_speed() {
    assert_eq!(WindConditions::default().vertical_speed, 0.0);
}

// MBA-728 coverage gap: `wind_shear::apply_boundary_layer_shear` preserves the vertical (Y)
// wind component unscaled while scaling horizontal (X/Z) wind. `TrajectorySolver::solve()`
// (used by every other test in this file) never actually calls that function -- its
// `enable_wind_shear` path applies an independent, hand-rolled equivalent inline in
// cli_api.rs's `get_wind_at_altitude`. The only public entry point that exercises
// `apply_boundary_layer_shear` itself is the fast-integration kernel
// (`fast_trajectory::fast_integrate`, via `compute_derivatives`), so this integration test
// drives that path directly instead.
fn fast_base_inputs(enable_wind_shear: bool, wind_shear_model: &str) -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: 800.0,
        bc_value: 0.5,
        bc_type: DragModel::G7,
        bullet_mass: 0.0109,
        bullet_diameter: 0.00782,
        bullet_length: 0.0309,
        enable_wind_shear,
        wind_shear_model: wind_shear_model.to_string(),
        ground_threshold: -1000.0,
        ..BallisticInputs::default()
    }
}

/// (y, z) at range 400 m for a fast-kernel run. The launch is a high-arc shot (~6.9 deg
/// elevation) so the trajectory climbs well above the ~10 m boundary-layer reference height
/// where shear actually departs from ratio == 1.0 -- a flat-fire shot never leaves the
/// near-ground "full wind" band, so shear would be a silent no-op and could not discriminate
/// the pass-through contract either way.
fn fast_yz_at_400(
    enable_wind_shear: bool,
    wind_shear_model: &str,
    wind_segments: Vec<WindSegment>,
) -> (f64, f64) {
    let inputs = fast_base_inputs(enable_wind_shear, wind_shear_model);
    let elevation = 0.12_f64;
    let solution = fast_integrate(
        &inputs,
        &WindSock::new(wind_segments),
        FastIntegrationParams {
            horiz: 400.0,
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
        },
    );
    let last = solution.t.len() - 1;
    assert_eq!(
        solution.y[0][last].to_bits(),
        400.0_f64.to_bits(),
        "fast_integrate should land exactly on the 400 m target plane"
    );
    (solution.y[1][last], solution.y[2][last])
}

#[test]
fn updraft_deflection_is_unaffected_by_shear() {
    // Same base wind as the vertical/lateral symmetry test above: a 5 m/s updraft. A modest
    // 3 m/s crosswind from the right is added so shear -- which only scales horizontal wind --
    // has something to act on; the two runs otherwise differ ONLY in enable_wind_shear.
    let wind_segments = vec![WindSegment {
        speed_kmh: 3.0 * 3.6,
        angle_deg: 90.0,
        until_m: 10_000.0,
        vertical_mps: 5.0,
    }];

    let (y_calm, z_calm) = fast_yz_at_400(false, "none", vec![]);
    let (y_no_shear, z_no_shear) = fast_yz_at_400(false, "none", wind_segments.clone());
    let (y_shear, z_shear) = fast_yz_at_400(true, "power_law", wind_segments);

    let dy_no_shear = y_no_shear - y_calm;
    let dy_shear = y_shear - y_calm;
    let dz_no_shear = z_no_shear - z_calm;
    let dz_shear = z_shear - z_calm;
    println!(
        "dy_no_shear={dy_no_shear:.6} dy_shear={dy_shear:.6} dz_no_shear={dz_no_shear:.6} dz_shear={dz_shear:.6}"
    );

    // Shear must not touch vertical: the updraft-caused vertical deflection at 400 m must
    // agree within 2% whether shear is on or off.
    let dy_err = (dy_shear - dy_no_shear).abs();
    let dy_bound = 0.02 * dy_no_shear.abs();
    assert!(
        dy_err <= dy_bound,
        "shear must not affect vertical deflection: dy_no_shear={dy_no_shear:.6}, dy_shear={dy_shear:.6}, err={dy_err:.6} > 2% bound {dy_bound:.6}"
    );

    // But shear DID scale the horizontal crosswind aloft, so the lateral deflection must
    // differ well beyond float noise -- proving shear was actually active in this run.
    let dz_diff = (dz_shear - dz_no_shear).abs();
    assert!(
        dz_diff > 1e-3,
        "shear should have changed the lateral drift (was it actually active?): dz_no_shear={dz_no_shear:.6}, dz_shear={dz_shear:.6}, diff={dz_diff:.6}"
    );
}
