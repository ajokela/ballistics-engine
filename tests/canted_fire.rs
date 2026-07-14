// MBA-1286 analytic gate: zero at 100 m un-canted, fire at 10 degrees of cant; the POI
// error vs the un-canted shot must match the small-angle canted-rifle prediction
//   dz = (D - sh) * sin(g)      (right for clockwise cant)
//   dy = -(D - sh) * (1-cos(g)) (low)
// where D is the MEASURED height the zero elevation adds at range R and sh is the
// sight height (the bore itself swings -sh*sin(g), partially offsetting the leak).
use ballistics_engine::{
    calculate_zero_angle, AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver,
    WindConditions,
};

fn base() -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: 800.0,
        bc_value: 0.5,
        bc_type: DragModel::G7,
        bullet_mass: 0.0109,
        bullet_diameter: 0.00782,
        bullet_length: 0.0309,
        sight_height: 0.05,
        twist_rate: 10.0,
        use_rk4: true,
        ..Default::default()
    }
}

fn solve(inputs: BallisticInputs, max_range: f64) -> ballistics_engine::TrajectoryResult {
    let mut s = TrajectorySolver::new(inputs, WindConditions::default(), AtmosphericConditions::default());
    s.set_max_range(max_range);
    s.set_time_step(0.001);
    s.solve().expect("solve")
}

fn yz_at(r: &ballistics_engine::TrajectoryResult, x: f64) -> (f64, f64) {
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
fn canted_fire_matches_small_angle_prediction() {
    let cant = 10f64.to_radians();
    let sh = base().sight_height;

    // Zero at 100 m, un-canted (zero_height 0 = bore-line datum used by the zero solver).
    // Relies on the engine zeroing un-canted (cant=0 forced in calculate_zero_angle);
    // that invariant is independently verified by cli_api.rs cant_tests::zero_angle_is_independent_of_cant.
    let zero_angle = calculate_zero_angle(base(), 100.0, 0.0).expect("zero");

    // Flat reference (no elevation), zeroed reference, and canted shot.
    let mut flat = base();
    flat.muzzle_angle = 0.0;
    let mut zeroed = base();
    zeroed.muzzle_angle = zero_angle;
    let mut canted = zeroed.clone();
    canted.cant_angle = cant;

    let r_flat = solve(flat, 700.0);
    let r_zeroed = solve(zeroed, 700.0);
    let r_canted = solve(canted, 700.0);

    for range in [300.0, 600.0] {
        let (y_flat, _) = yz_at(&r_flat, range);
        let (y_zeroed, z_zeroed) = yz_at(&r_zeroed, range);
        let (y_canted, z_canted) = yz_at(&r_canted, range);

        let d = y_zeroed - y_flat; // measured elevation-correction height at this range
        assert!(d > 0.1, "sanity: zero elevation must add height at {range} m (d={d})");

        let dz = z_canted - z_zeroed;
        let dy = y_canted - y_zeroed;
        let expected_dz = (d - sh) * cant.sin();
        let expected_dy = -(d - sh) * (1.0 - cant.cos());

        assert!(
            (dz - expected_dz).abs() <= 0.05 * expected_dz.abs(),
            "horizontal cant error at {range} m: got {dz:.4}, expected {expected_dz:.4} (5%)"
        );
        assert!(
            (dy - expected_dy).abs() <= 0.10 * expected_dy.abs(),
            "vertical cant error at {range} m: got {dy:.4}, expected {expected_dy:.4} (10%)"
        );
    }
}
