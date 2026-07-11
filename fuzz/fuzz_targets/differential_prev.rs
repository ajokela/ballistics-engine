#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{
    AtmosphericConditions as Atmo, BallisticInputs as Inp, TrajectorySolver as Solver,
    WindConditions as Wind,
};
use ballistics_engine_fuzz::domain::valid_inputs;

use engine_prev::{
    AtmosphericConditions as PAtmo, BallisticInputs as PInp, TrajectorySolver as PSolver,
    WindConditions as PWind,
};

/// Sign of drop response to a +BC bump, in the CURRENT engine.
fn response_sign_current(base: &Inp) -> Option<i8> {
    let range = (base.target_distance * 0.6).min(800.0);
    let y0 = drop_at_current(base, range)?;
    let mut hi = base.clone();
    hi.bc_value *= 1.5;
    let y1 = drop_at_current(&hi, range)?;
    Some(sign(y1 - y0))
}

fn response_sign_prev(base: &PInp) -> Option<i8> {
    let range = (base.target_distance * 0.6).min(800.0);
    let y0 = drop_at_prev(base, range)?;
    let mut hi = base.clone();
    hi.bc_value *= 1.5;
    let y1 = drop_at_prev(&hi, range)?;
    Some(sign(y1 - y0))
}

fn sign(x: f64) -> i8 { if x > 1e-4 { 1 } else if x < -1e-4 { -1 } else { 0 } }

fn drop_at_current(inputs: &Inp, range: f64) -> Option<f64> {
    let r = Solver::new(inputs.clone(), Wind::default(), Atmo::default()).solve().ok()?;
    interp_y(r.points.iter().map(|p| (p.position.x, p.position.y)), range)
}
fn drop_at_prev(inputs: &PInp, range: f64) -> Option<f64> {
    let r = PSolver::new(inputs.clone(), PWind::default(), PAtmo::default()).solve().ok()?;
    interp_y(r.points.iter().map(|p| (p.position.x, p.position.y)), range)
}

fn interp_y<I: Iterator<Item = (f64, f64)>>(pts: I, range: f64) -> Option<f64> {
    let v: Vec<(f64, f64)> = pts.collect();
    for i in 1..v.len() {
        if v[i].0 >= range {
            let (x1, y1) = v[i - 1];
            let (x2, y2) = v[i];
            let dx = x2 - x1;
            if dx.abs() < 1e-9 { return Some(y1); }
            let t = (range - x1) / dx;
            return Some(y1 + t * (y2 - y1));
        }
    }
    None
}

/// Mirror the physically-meaningful scalar fields into the 0.21.5 input struct.
/// Uses `..Default::default()` for the (identical) remaining fields.
fn to_prev(c: &Inp) -> PInp {
    PInp {
        bc_value: c.bc_value,
        bullet_mass: c.bullet_mass,
        bullet_diameter: c.bullet_diameter,
        bullet_length: c.bullet_length,
        muzzle_velocity: c.muzzle_velocity,
        muzzle_angle: c.muzzle_angle,
        target_distance: c.target_distance,
        twist_rate: c.twist_rate,
        temperature: c.temperature,
        pressure: c.pressure,
        altitude: c.altitude,
        humidity: c.humidity,
        use_rk4: true,
        ..Default::default()
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(cur) = valid_inputs(&mut u) else { return };
    let prev = to_prev(&cur);

    // RELATIONSHIP MODE: both versions must agree on the SIGN of the drop response
    // to a BC increase. Absolute values legitimately differ across 0.21->0.22.
    if let (Some(sc), Some(sp)) = (response_sign_current(&cur), response_sign_prev(&prev)) {
        // We assert only that the two engine versions AGREE on the sign of the
        // BC->drop response (relationship-mode), not the direction itself.
        // Disagreement on a NON-zero sign is a regression in a qualitative property.
        if sc != 0 && sp != 0 {
            assert_eq!(sc, sp,
                "BC->drop response sign diverged: current={sc}, prev(0.21.5)={sp}");
        }
    }
});
