#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain::valid_inputs;

/// Vertical position (y, McCoy: negative = drop) interpolated at downrange x=`range`.
/// None if the trajectory never reaches `range` (short shot) or the solve errors.
fn drop_at(inputs: &BallisticInputs, range: f64) -> Option<f64> {
    let solver = TrajectorySolver::new(inputs.clone(), WindConditions::default(), AtmosphericConditions::default());
    let r = solver.solve().ok()?;
    let pts = &r.points;
    for i in 1..pts.len() {
        if pts[i].position.x >= range {
            let p1 = &pts[i - 1];
            let p2 = &pts[i];
            let dx = p2.position.x - p1.position.x;
            if dx.abs() < 1e-9 { return Some(p1.position.y); }
            let t = (range - p1.position.x) / dx;
            return Some(p1.position.y + t * (p2.position.y - p1.position.y));
        }
    }
    None
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(mut base) = valid_inputs(&mut u) else { return };
    // Flat fire so the monotonicity properties are well-defined: a high launch angle
    // puts the apex far downrange, making "farther => more drop" false in the ascending
    // arc (a harness artifact, not an engine bug). Angle isn't the variable under test
    // here — BC and range are — so pin it near-horizontal and vary everything else.
    base.muzzle_angle = 0.0;

    // Compare drop at a range both shots plausibly reach.
    let range = (base.target_distance * 0.6).min(800.0);

    // Property A: higher BC => equal-or-less drop (|y| smaller) at fixed range.
    let mut hi_bc = base.clone();
    hi_bc.bc_value = base.bc_value * 1.5;
    if let (Some(y_lo), Some(y_hi)) = (drop_at(&base, range), drop_at(&hi_bc, range)) {
        // y is negative going down; higher BC should not drop MORE (more negative).
        let tol = 1e-3; // meters, float slack
        assert!(y_hi >= y_lo - tol,
            "higher BC dropped more at {range} m: base_y={y_lo}, hi_bc_y={y_hi}");
    }

    // Property B: farther downrange => equal-or-more drop (monotone descent after apex).
    // Compare two ranges past a modest distance to stay beyond the ascending arc.
    let near = 200.0_f64.min(base.target_distance * 0.4);
    let far = (near * 2.0).min(base.target_distance * 0.9);
    if far > near + 1.0 {
        if let (Some(y_near), Some(y_far)) = (drop_at(&base, near), drop_at(&base, far)) {
            let tol = 1e-2;
            assert!(y_far <= y_near + tol,
                "trajectory rose from {near}->{far} m past apex: y_near={y_near}, y_far={y_far}");
        }
    }
});
