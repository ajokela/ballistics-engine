#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};
use ballistics_engine_fuzz::domain::{ranged, valid_inputs};

/// Lateral (z) position of the last trajectory point, i.e. total windage at impact.
fn lateral_at_impact(inputs: &BallisticInputs, wind: WindConditions) -> Option<f64> {
    let solver = TrajectorySolver::new(inputs.clone(), wind, AtmosphericConditions::default());
    let r = solver.solve().ok()?;
    r.points.last().map(|p| p.position.z)
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(mut base) = valid_inputs(&mut u) else { return };

    // Isolate wind: disable spin drift, Magnus, Coriolis so lateral motion is wind-only.
    base.enable_coriolis = false;
    base.enable_magnus = false;
    base.use_enhanced_spin_drift = false;
    base.latitude = None;

    // Property A: determinism — identical inputs produce bit-identical lateral impact.
    let z1 = lateral_at_impact(&base, WindConditions::default());
    let z2 = lateral_at_impact(&base, WindConditions::default());
    assert_eq!(z1.map(f64::to_bits), z2.map(f64::to_bits), "solve is non-deterministic");

    // Property B: zero wind (+ no spin/Magnus/Coriolis) => ~zero lateral.
    if let Some(z0) = lateral_at_impact(&base, WindConditions::default()) {
        assert!(z0.abs() < 1e-3, "nonzero drift with zero wind/spin/Coriolis: {z0}");
    }

    // Property C: mirrored crosswind => mirrored drift.
    let Ok(speed) = ranged(&mut u, 1.0, 15.0) else { return };
    let from_right = WindConditions { speed, direction: std::f64::consts::FRAC_PI_2, vertical_speed: 0.0 };      // PI/2
    let from_left = WindConditions { speed, direction: 3.0 * std::f64::consts::FRAC_PI_2, vertical_speed: 0.0 }; // 3PI/2
    if let (Some(zr), Some(zl)) = (lateral_at_impact(&base, from_right), lateral_at_impact(&base, from_left)) {
        // Opposite crosswinds => opposite-sign, equal-magnitude drift.
        let tol = 1e-2 + 0.02 * zr.abs();
        assert!((zr + zl).abs() <= tol, "crosswind not antisymmetric: from_right={zr}, from_left={zl}");
    }
});
