//! Shared byte->input generators and invariant assertions for the fuzz harnesses.

use arbitrary::{Result, Unstructured};
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, DragModel, TrajectoryResult,
    HUMIDITY_SENSOR_TOLERANCE_PERCENT,
};

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

/// A launch angle that is either non-finite or finite within the forward-fire domain. Keeping
/// finite angles inside +/- 80 degrees makes `assert_finite_sane`'s non-negative downrange
/// invariant meaningful while still exercising the solve boundary's NaN/Inf rejection.
fn hostile_launch_angle(u: &mut Unstructured) -> Result<f64> {
    Ok(match u.int_in_range(0u8..=3)? {
        0 => f64::NAN,
        1 => f64::INFINITY,
        2 => f64::NEG_INFINITY,
        _ => ranged(u, -1.396_263_401_595_463_6, 1.396_263_401_595_463_6)?,
    })
}

/// A relative humidity, percent, that every surface accepts: anywhere in the sensor band
/// (`-HUMIDITY_SENSOR_TOLERANCE_PERCENT..=100 + HUMIDITY_SENSOR_TOLERANCE_PERCENT`, read as 0 or
/// 100 outside 0-100), with the band's edges, 0 and 100 drawn on purpose.
pub fn band_humidity_percent(u: &mut Unstructured) -> Result<f64> {
    let band = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    Ok(match u.int_in_range(0u8..=9)? {
        0 => -band,
        1 => 0.0,
        2 => 100.0,
        3 => 100.0 + band,
        4 => ranged(u, -band, 0.0)?,
        5 => ranged(u, 100.0, 100.0 + band)?,
        _ => ranged(u, 0.0, 100.0)?,
    })
}

/// A relative humidity, percent, at or just beyond the sensor band's edges: the last values
/// accepted and the first refused, and a sweep a little past each side.
pub fn band_edge_humidity_percent(u: &mut Unstructured) -> Result<f64> {
    let band = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    let below = -band;
    let above = 100.0 + band;
    Ok(match u.int_in_range(0u8..=5)? {
        0 => below,
        1 => above,
        // One ulp beyond each edge: refused.
        2 => f64::from_bits(below.to_bits() + 1),
        3 => f64::from_bits(above.to_bits() + 1),
        4 => ranged(u, below - 2.0, below + 2.0)?,
        _ => ranged(u, above - 2.0, above + 2.0)?,
    })
}

/// The station air a draw describes: its own temperature, pressure, humidity and altitude. A
/// harness must fly this, not `AtmosphericConditions::default()`: `TrajectorySolver` sets the
/// inputs' humidity from its atmosphere, so a default atmosphere flies 15 C, 1013.25 hPa and
/// 50 % whatever was drawn, and the drawn air never reaches the validator or the physics.
pub fn station_air(inputs: &BallisticInputs) -> AtmosphericConditions {
    AtmosphericConditions {
        temperature: inputs.temperature,
        pressure: inputs.pressure,
        humidity: inputs.humidity_percent,
        altitude: inputs.altitude,
    }
}

/// The air `robustness_inputs` flies a hostile draw in. Each mode reaches different checks:
/// - the default air, so the non-air fields' validation is reached even when the drawn air would
///   be refused first;
/// - the drawn air itself, wild temperature, pressure, altitude and humidity together;
/// - a plausible station with the drawn (wild) humidity, so a hostile humidity reaches the band
///   rule and the station-air construction instead of losing to a refused temperature;
/// - a plausible station at a band edge, just inside or just outside;
/// - a hot, saturated, low-pressure corner, where saturated air stops existing.
pub fn hostile_air(
    u: &mut Unstructured,
    inputs: &BallisticInputs,
) -> Result<AtmosphericConditions> {
    let plausible = |u: &mut Unstructured, humidity: f64| -> Result<AtmosphericConditions> {
        Ok(AtmosphericConditions {
            temperature: ranged(u, -40.0, 50.0)?,
            pressure: ranged(u, 700.0, 1100.0)?,
            humidity,
            altitude: ranged(u, 0.0, 3000.0)?,
        })
    };
    Ok(match u.int_in_range(0u8..=4)? {
        0 => AtmosphericConditions::default(),
        1 => station_air(inputs),
        2 => plausible(u, inputs.humidity_percent)?,
        3 => {
            let humidity = band_edge_humidity_percent(u)?;
            plausible(u, humidity)?
        }
        _ => AtmosphericConditions {
            temperature: ranged(u, 80.0, 110.0)?,
            pressure: ranged(u, 300.0, 1013.25)?,
            humidity: ranged(u, 90.0, 100.0 + HUMIDITY_SENSOR_TOLERANCE_PERCENT)?,
            altitude: ranged(u, 0.0, 3000.0)?,
        },
    })
}

/// A ballistically plausible input. Every field is finite and in a real-world
/// envelope so invariant/differential/solver harnesses spend budget on meaningful
/// solves rather than rejected garbage.
pub fn valid_inputs(u: &mut Unstructured) -> Result<BallisticInputs> {
    Ok(BallisticInputs {
        bc_value: ranged(u, 0.1, 1.2)?,
        bc_type: if u.arbitrary()? {
            DragModel::G7
        } else {
            DragModel::G1
        },
        bullet_mass: ranged(u, 0.001, 0.05)?, // kg (~15..770 gr)
        bullet_diameter: ranged(u, 0.004, 0.014)?, // m (.17...55 cal)
        bullet_length: ranged(u, 0.010, 0.050)?, // m
        muzzle_velocity: ranged(u, 200.0, 1500.0)?, // m/s
        muzzle_angle: ranged(u, -0.05, 0.30)?, // rad
        target_distance: ranged(u, 50.0, 2000.0)?, // m (loop-count bound)
        twist_rate: ranged(u, 5.0, 20.0)?,    // in/turn, > 0
        temperature: ranged(u, -40.0, 50.0)?, // C
        pressure: ranged(u, 800.0, 1100.0)?,  // hPa
        altitude: ranged(u, 0.0, 4000.0)?,    // m
        humidity_percent: band_humidity_percent(u)?, // PERCENT, the sensor band included
        use_rk4: true,
        ..Default::default()
    })
}

/// Physics fields are `wild` (NaN/Inf/negative/extreme) to test clean rejection;
/// loop-count fields (`target_distance`) stay bounded so a hostile input can't turn
/// into a multi-minute solve. Runtime termination is covered by `solver_zero` +
/// libFuzzer's own `-timeout`.
pub fn hostile_inputs(u: &mut Unstructured) -> Result<BallisticInputs> {
    Ok(BallisticInputs {
        bc_value: wild(u)?,
        bullet_mass: wild(u)?,
        bullet_diameter: wild(u)?,
        bullet_length: wild(u)?,
        muzzle_velocity: wild(u)?,
        muzzle_angle: hostile_launch_angle(u)?,
        twist_rate: wild(u)?,
        temperature: wild(u)?,
        pressure: wild(u)?,
        altitude: wild(u)?,
        humidity_percent: wild(u)?,
        target_distance: ranged(u, 10.0, 3000.0)?, // bounded on purpose
        use_rk4: true,
        ..Default::default()
    })
    .and_then(|mut inputs| {
        // Drawn last, so every field above decodes as it always did (the frozen reproducers
        // in tests/ decode through this function): in a third of draws the humidity becomes a
        // sensor-band edge or a value just beyond one, which `wild` almost never produces.
        if u.int_in_range(0u8..=2)? == 2 {
            inputs.humidity_percent = band_edge_humidity_percent(u)?;
        }
        Ok(inputs)
    })
}

/// Every scalar output is finite and physically non-negative; every trajectory
/// point is finite. Panics (= a libFuzzer crash) on violation.
pub fn assert_finite_sane(r: &TrajectoryResult) {
    assert!(
        r.max_range.is_finite() && r.max_range >= 0.0,
        "bad max_range {}",
        r.max_range
    );
    assert!(r.max_height.is_finite(), "max_height not finite");
    assert!(
        r.time_of_flight.is_finite() && r.time_of_flight >= 0.0,
        "bad time_of_flight {}",
        r.time_of_flight
    );
    assert!(
        r.impact_velocity.is_finite() && r.impact_velocity >= 0.0,
        "bad impact_velocity {}",
        r.impact_velocity
    );
    assert!(
        r.impact_energy.is_finite() && r.impact_energy >= 0.0,
        "bad impact_energy {}",
        r.impact_energy
    );
    for (i, p) in r.points.iter().enumerate() {
        assert!(
            p.position.x.is_finite() && p.position.y.is_finite() && p.position.z.is_finite(),
            "point {i} position not finite"
        );
        assert!(
            p.velocity_magnitude.is_finite() && p.velocity_magnitude >= 0.0,
            "point {i} bad velocity {}",
            p.velocity_magnitude
        );
        assert!(
            p.time.is_finite() && p.time >= 0.0,
            "point {i} bad time {}",
            p.time
        );
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
    fn valid_inputs_stay_inside_the_sensor_band_and_reach_its_edges() {
        let band = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
        let mut seen = std::collections::BTreeSet::new();
        for seed in 0u8..=255 {
            let bytes: Vec<u8> = (0..4096)
                .map(|i| seed.wrapping_mul(31).wrapping_add(i as u8))
                .collect();
            let mut u = Unstructured::new(&bytes);
            let humidity = valid_inputs(&mut u).unwrap().humidity_percent;
            assert!((-band..=100.0 + band).contains(&humidity), "{humidity}");
            seen.insert(humidity.to_bits());
        }
        for edge in [-band, 0.0, 100.0, 100.0 + band] {
            assert!(seen.contains(&edge.to_bits()), "edge {edge} never drawn");
        }
    }

    /// The point of `station_air`: a harness that flies it delivers the drawn humidity to the
    /// solver. A NaN humidity in an otherwise plausible draw is refused, naming the humidity.
    /// (`AtmosphericConditions::default()`, which every robustness and invariant harness flew
    /// until 0.50.0, flies 50 % whatever was drawn, and the same draw solves.)
    #[test]
    fn the_drawn_humidity_reaches_the_solver() {
        use ballistics_engine::{TrajectorySolver, WindConditions};
        let bytes = vec![0xA5u8; 4096];
        let mut u = Unstructured::new(&bytes);
        let mut inputs = valid_inputs(&mut u).unwrap();
        inputs.target_distance = 300.0;
        assert!(
            TrajectorySolver::new(
                inputs.clone(),
                WindConditions::default(),
                station_air(&inputs)
            )
            .solve()
            .is_ok(),
            "the plausible draw solves in its own air"
        );
        inputs.humidity_percent = f64::NAN;
        let refused = TrajectorySolver::new(
            inputs.clone(),
            WindConditions::default(),
            station_air(&inputs),
        )
        .solve()
        .expect_err("a NaN humidity reached the solver and was refused");
        assert!(
            refused.to_string().contains("humidity"),
            "the refusal names the humidity: {refused}"
        );
    }

    #[test]
    fn valid_inputs_loop_bounds_are_capped() {
        let bytes = vec![0xFFu8; 4096];
        let mut u = Unstructured::new(&bytes);
        let b = valid_inputs(&mut u).unwrap();
        assert!(b.target_distance <= 2000.0);
    }
}
