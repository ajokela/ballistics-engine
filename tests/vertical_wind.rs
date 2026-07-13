// MBA-728 analytic acceptance gate for vertical wind. Drag is direction-agnostic, so a pure
// updraft `w` must deflect the shot VERTICALLY by the same amount (within 5%) that an equal
// pure crosswind `w` deflects it LATERALLY, at the same ranges.
use ballistics_engine::wind::WindSegment;
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
