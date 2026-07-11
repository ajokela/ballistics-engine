#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{
    calculate_zero_angle, AtmosphericConditions, TrajectorySolver, WindConditions,
};
use ballistics_engine_fuzz::domain::valid_inputs;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(inputs) = valid_inputs(&mut u) else { return };
    let target_distance = inputs.target_distance;
    let target_height = 0.0; // zero at the same height as the bore line

    // The solver must either cleanly report "un-zeroable" (Err) or return a FINITE,
    // physically-plausible launch angle. NaN/Inf/absurd angle = bug.
    match calculate_zero_angle(inputs.clone(), target_distance, target_height) {
        Err(_) => { /* clean rejection is acceptable */ }
        Ok(angle) => {
            assert!(angle.is_finite(), "zero angle not finite: {angle}");
            // The zero search brackets [0, 0.785] rad (it expands high_angle up to 45°),
            // so any returned angle is <= ~0.785. Flag only truly out-of-range values.
            assert!(angle.abs() < 0.8, "absurd zero angle {angle} rad for {target_distance} m");

            // Round-trip: firing at the solved angle should land ~on the target height at
            // target_distance (the solver's own convergence contract). Reproduce the zero
            // search's own integration settings (time_step 0.001, max_range = 2*distance) so
            // the round-trip trajectory matches what the search converged on.
            let mut z = inputs.clone();
            z.muzzle_angle = angle;
            z.enable_aerodynamic_jump = false;
            let mut solver = TrajectorySolver::new(z, WindConditions::default(), AtmosphericConditions::default());
            solver.set_max_range(target_distance * 2.0);
            solver.set_time_step(0.001);
            if let Ok(r) = solver.solve() {
                // interpolate y at target_distance
                let mut y_at = None;
                for i in 1..r.points.len() {
                    if r.points[i].position.x >= target_distance {
                        let p1 = &r.points[i - 1];
                        let p2 = &r.points[i];
                        let dx = p2.position.x - p1.position.x;
                        let t = if dx.abs() < 1e-9 { 0.0 } else { (target_distance - p1.position.x) / dx };
                        y_at = Some(p1.position.y + t * (p2.position.y - p1.position.y));
                        break;
                    }
                }
                if let Some(y) = y_at {
                    // Convergence tolerance: the zero search aims for the target height.
                    // 0.05 m is generous slack over the solver's internal precision.
                    assert!((y - target_height).abs() < 0.05,
                        "zero did not round-trip: y={y} at {target_distance} m, wanted {target_height} (angle={angle})");
                }
            }
        }
    }
});
