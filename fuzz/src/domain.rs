//! Shared byte->input generators and invariant assertions for the fuzz harnesses.

use arbitrary::{Result, Unstructured};
use ballistics_engine::{BallisticInputs, DragModel, TrajectoryResult};

/// Map raw bytes to a finite f64 in `[lo, hi]` (deterministic).
pub fn ranged(u: &mut Unstructured, lo: f64, hi: f64) -> Result<f64> {
    let bits: u64 = u.arbitrary()?;
    let frac = (bits as f64) / (u64::MAX as f64);
    Ok(lo + frac * (hi - lo))
}

/// An adversarial f64: NaN, infinities, zero, negatives, huge magnitudes, or an
/// arbitrary bit pattern. Used for the PHYSICS fields in `hostile_inputs`.
pub fn wild(u: &mut Unstructured) -> Result<f64> {
    Ok(match u.int_in_range(0u8..=7)? {
        0 => f64::NAN,
        1 => f64::INFINITY,
        2 => f64::NEG_INFINITY,
        3 => 0.0,
        4 => -ranged(u, 0.0, 1.0e6)?,
        5 => ranged(u, 0.0, 1.0e9)?,
        6 => ranged(u, -1.0e300, 1.0e300)?,
        _ => f64::from_bits(u.arbitrary()?),
    })
}

/// A ballistically plausible input. Every field is finite and in a real-world
/// envelope so invariant/differential/solver harnesses spend budget on meaningful
/// solves rather than rejected garbage.
pub fn valid_inputs(u: &mut Unstructured) -> Result<BallisticInputs> {
    let mut b = BallisticInputs::default();
    b.bc_value = ranged(u, 0.1, 1.2)?;
    b.bc_type = if u.arbitrary()? { DragModel::G7 } else { DragModel::G1 };
    b.bullet_mass = ranged(u, 0.001, 0.05)?;      // kg (~15..770 gr)
    b.bullet_diameter = ranged(u, 0.004, 0.014)?; // m (.17...55 cal)
    b.bullet_length = ranged(u, 0.010, 0.050)?;   // m
    b.muzzle_velocity = ranged(u, 200.0, 1500.0)?;// m/s
    b.muzzle_angle = ranged(u, -0.05, 0.30)?;     // rad
    b.target_distance = ranged(u, 50.0, 2000.0)?; // m (loop-count bound)
    b.twist_rate = ranged(u, 5.0, 20.0)?;         // in/turn, > 0
    b.temperature = ranged(u, -40.0, 50.0)?;      // C
    b.pressure = ranged(u, 800.0, 1100.0)?;       // hPa
    b.altitude = ranged(u, 0.0, 4000.0)?;         // m
    b.humidity = ranged(u, 0.0, 1.0)?;            // FRACTION
    b.use_rk4 = true;
    Ok(b)
}

/// Physics fields are `wild` (NaN/Inf/negative/extreme) to test clean rejection;
/// loop-count fields (`target_distance`) stay bounded so a hostile input can't turn
/// into a multi-minute solve. Runtime termination is covered by `solver_zero` +
/// libFuzzer's own `-timeout`.
pub fn hostile_inputs(u: &mut Unstructured) -> Result<BallisticInputs> {
    let mut b = BallisticInputs::default();
    b.bc_value = wild(u)?;
    b.bullet_mass = wild(u)?;
    b.bullet_diameter = wild(u)?;
    b.bullet_length = wild(u)?;
    b.muzzle_velocity = wild(u)?;
    b.muzzle_angle = wild(u)?;
    b.twist_rate = wild(u)?;
    b.temperature = wild(u)?;
    b.pressure = wild(u)?;
    b.altitude = wild(u)?;
    b.humidity = wild(u)?;
    b.target_distance = ranged(u, 10.0, 3000.0)?; // bounded on purpose
    b.use_rk4 = true;
    Ok(b)
}

/// Every scalar output is finite and physically non-negative; every trajectory
/// point is finite. Panics (= a libFuzzer crash) on violation.
pub fn assert_finite_sane(r: &TrajectoryResult) {
    assert!(r.max_range.is_finite() && r.max_range >= 0.0, "bad max_range {}", r.max_range);
    assert!(r.max_height.is_finite(), "max_height not finite");
    assert!(r.time_of_flight.is_finite() && r.time_of_flight >= 0.0, "bad time_of_flight {}", r.time_of_flight);
    assert!(r.impact_velocity.is_finite() && r.impact_velocity >= 0.0, "bad impact_velocity {}", r.impact_velocity);
    assert!(r.impact_energy.is_finite() && r.impact_energy >= 0.0, "bad impact_energy {}", r.impact_energy);
    for (i, p) in r.points.iter().enumerate() {
        assert!(p.position.x.is_finite() && p.position.y.is_finite() && p.position.z.is_finite(), "point {i} position not finite");
        assert!(p.velocity_magnitude.is_finite() && p.velocity_magnitude >= 0.0, "point {i} bad velocity {}", p.velocity_magnitude);
        assert!(p.time.is_finite() && p.time >= 0.0, "point {i} bad time {}", p.time);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arbitrary::Unstructured;

    #[test]
    fn valid_inputs_are_finite_and_positive() {
        let bytes = vec![0xA5u8; 4096];
        let mut u = Unstructured::new(&bytes);
        let b = valid_inputs(&mut u).unwrap();
        assert!(b.bc_value.is_finite() && b.bc_value > 0.0);
        assert!(b.bullet_mass.is_finite() && b.bullet_mass > 0.0);
        assert!(b.muzzle_velocity.is_finite() && b.muzzle_velocity > 0.0);
        assert!(b.twist_rate.is_finite() && b.twist_rate > 0.0);
        assert!(b.target_distance.is_finite() && b.target_distance > 0.0);
    }

    #[test]
    fn valid_inputs_loop_bounds_are_capped() {
        let bytes = vec![0xFFu8; 4096];
        let mut u = Unstructured::new(&bytes);
        let b = valid_inputs(&mut u).unwrap();
        assert!(b.target_distance <= 2000.0);
    }
}
