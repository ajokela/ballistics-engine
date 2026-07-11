#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain::ranged;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(velocity) = ranged(&mut u, 300.0, 1000.0) else { return };

    // Flat fire, near-drag-free: compare the DROP at a fixed downrange x against the exact
    // vacuum parabola drop = 0.5*g*(x/v)^2. bc_value=100 makes retardation negligible, so
    // the engine's drop should match to a small band. This is a ground-truth analytic check
    // measured at a fixed x (interpolated) — independent of the solver's ground/max-range
    // termination, so unlike a time-of-flight comparison it is well-defined.
    let mut b = BallisticInputs::default();
    b.bc_value = 100.0;            // near-drag-free: huge retardation denominator
    b.bc_type = DragModel::G1;
    b.bullet_mass = 0.01;
    b.bullet_diameter = 0.0077;
    b.bullet_length = 0.03;
    b.muzzle_velocity = velocity;
    b.muzzle_angle = 0.0;         // flat fire: y=0 at muzzle, drop measured off the bore line
    b.target_distance = 600.0;
    b.twist_rate = 10.0;
    b.use_rk4 = true;
    b.enable_coriolis = false;
    b.enable_magnus = false;
    b.use_enhanced_spin_drift = false;
    b.latitude = None;

    let solver = TrajectorySolver::new(b, WindConditions::default(), AtmosphericConditions::default());
    let Ok(r) = solver.solve() else { return };

    // (a) No lateral motion without wind/spin/Coriolis — exact.
    if let Some(last) = r.points.last() {
        assert!(last.position.z.abs() < 1e-6, "lateral drift in a no-cross-effect shot: {}", last.position.z);
    }

    // (b) Vacuum-drop check: interpolate y at a downrange x well inside the flight.
    let x = 300.0_f64; // meters, < target_distance and < the default max_range
    let mut y_at = None;
    for i in 1..r.points.len() {
        if r.points[i].position.x >= x {
            let p1 = &r.points[i - 1];
            let p2 = &r.points[i];
            let dx = p2.position.x - p1.position.x;
            let t = if dx.abs() < 1e-9 { 0.0 } else { (x - p1.position.x) / dx };
            y_at = Some(p1.position.y + t * (p2.position.y - p1.position.y));
            break;
        }
    }
    if let Some(y) = y_at {
        let g = 9.80665;
        let t_flight = x / velocity;                 // flat fire, cos(0)=1
        let drop_vac = 0.5 * g * t_flight * t_flight; // exact vacuum drop off the bore line
        let drop_eng = -y;                            // y is negative (below bore)
        // Near-drag-free: engine drop should sit within 15% + 2 cm of the vacuum drop.
        // Residual drag only slightly increases flight time -> slightly more drop; a value
        // far outside this band signals a gravity/integration error.
        let tol = 0.15 * drop_vac + 0.02;
        assert!((drop_eng - drop_vac).abs() <= tol,
            "near-vacuum drop {drop_eng:.4} m at {x} m deviates from vacuum {drop_vac:.4} m (v={velocity})");
    }
});
