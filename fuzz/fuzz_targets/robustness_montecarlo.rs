#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::{
    run_monte_carlo_with_wind, MonteCarloParams, WindConditions, TARGET_NOT_REACHED_SENTINEL_M,
};
use ballistics_engine_fuzz::domain::{ranged, valid_inputs};

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(inputs) = valid_inputs(&mut u) else { return };
    // Realistic positive std-devs. (The engine already rejects negative std-devs via
    // Normal::new(...).map_err(..)? = Err, verified; and a very wide bc/velocity std-dev
    // would sample non-physical params (negative BC), which is Finding #1's validation
    // gap, not a new MC bug. So we fuzz within plausible spreads and assert MC output
    // soundness: every realistic MC run yields finite ranges/velocities/positions.)
    let Ok(vsd) = ranged(&mut u, 0.0, 15.0) else { return };
    let Ok(asd) = ranged(&mut u, 0.0, 0.005) else { return };
    let Ok(bsd) = ranged(&mut u, 0.0, 0.03) else { return };
    let Ok(wsd) = ranged(&mut u, 0.0, 5.0) else { return };
    let params = MonteCarloParams {
        num_simulations: 16,
        velocity_std_dev: vsd,
        angle_std_dev: asd,
        bc_std_dev: bsd,
        wind_speed_std_dev: wsd,
        target_distance: Some(inputs.target_distance),
        base_wind_speed: 0.0,
        base_wind_direction: 0.0,
        azimuth_std_dev: 0.0,
    };
    if let Ok(results) = run_monte_carlo_with_wind(inputs, WindConditions::default(), params) {
        for r in &results.ranges {
            assert!(r.is_finite(), "MC range not finite: {r}");
        }
        for v in &results.impact_velocities {
            assert!(v.is_finite() && *v >= 0.0, "MC impact velocity bad: {v}");
        }
        for p in &results.impact_positions {
            // A miss is the sentinel row; every other component must be finite.
            let is_sentinel = p.y == TARGET_NOT_REACHED_SENTINEL_M;
            assert!(is_sentinel || (p.x.is_finite() && p.y.is_finite() && p.z.is_finite()),
                "MC impact position not finite: {:?}", p);
        }
    }
});
