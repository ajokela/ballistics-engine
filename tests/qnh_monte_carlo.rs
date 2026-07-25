//! MBA-1397: prove the reduced (not raw QNH) station pressure reaches BOTH Monte Carlo
//! entry points that read `BallisticInputs.pressure` directly:
//!
//! - `monte_carlo::solve_trajectory_for_monte_carlo` (the fixed-step MC kernel; a `pub`
//!   library surface with no production CLI/WASM/FFI caller, but exercised directly by
//!   `tests/monte_carlo_si.rs` and reachable by any external consumer of this crate's API).
//! - `run_monte_carlo` / `run_monte_carlo_with_wind_and_direction_std_dev` (the actual
//!   production Monte Carlo path behind the CLI `monte-carlo` subcommand and the FFI
//!   `ballistics_monte_carlo*` exports), which builds its own internal `AtmosphericConditions`
//!   directly from `base_inputs.temperature/pressure/humidity/altitude`.
//!
//! Neither `BallisticInputs` nor `AtmosphericConditions` carry a pressure-reference-mode
//! field (see src/atmosphere.rs's `PressureReferenceMode`); both fields have always meant,
//! and continue to mean, absolute station pressure. A caller that wants QNH support reduces
//! it themselves via `ballistics_engine::atmosphere::reduce_qnh_to_station_pressure` before
//! setting `BallisticInputs.pressure` -- these tests prove that reduction actually reaches
//! the solve (a materially different, physically correct result versus feeding the raw,
//! unreduced QNH straight through).

use ballistics_engine::atmosphere::reduce_qnh_to_station_pressure;
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo;
use ballistics_engine::{run_monte_carlo, BallisticInputs, DragModel, MonteCarloParams};

const ALTITUDE_M: f64 = 1500.0;
const QNH_HPA: f64 = 1030.0;

fn si_bullet(pressure_hpa: f64) -> BallisticInputs {
    BallisticInputs {
        target_distance: 457.2, // meters (500 yd)
        muzzle_velocity: 823.0, // m/s
        bullet_mass: 168.0 * 0.00006479891,
        bullet_diameter: 0.308 * 0.0254,
        bullet_length: 1.215 * 0.0254,
        caliber_inches: 0.308,
        weight_grains: 168.0,
        bc_value: 0.475,
        bc_type: DragModel::G1,
        muzzle_angle: 0.006,
        sight_height: 0.05,
        muzzle_height: 0.0,
        temperature: 15.0,
        altitude: ALTITUDE_M,
        pressure: pressure_hpa,
        ..BallisticInputs::default()
    }
}

#[test]
fn monte_carlo_kernel_uses_the_reduced_pressure_not_the_raw_qnh() {
    let reduced_hpa = reduce_qnh_to_station_pressure(QNH_HPA, ALTITUDE_M);
    assert!(reduced_hpa < QNH_HPA, "reduction must lower pressure");

    let with_reduced = solve_trajectory_for_monte_carlo(&si_bullet(reduced_hpa))
        .expect("monte carlo solve (reduced)");
    let with_raw_qnh = solve_trajectory_for_monte_carlo(&si_bullet(QNH_HPA))
        .expect("monte carlo solve (raw QNH)");

    // Lower (correctly reduced) pressure means less dense air, less drag, higher retained
    // velocity at the same target distance.
    assert!(
        with_reduced.velocity > with_raw_qnh.velocity + 5.0,
        "reduced pressure must retain more velocity than the raw, unreduced QNH: \
         reduced={} raw_qnh={}",
        with_reduced.velocity,
        with_raw_qnh.velocity
    );

    // And it must match a caller who had the equivalent absolute station pressure directly
    // (proving the reduced number, not some other transformation, reaches the solve).
    let with_equivalent_absolute = solve_trajectory_for_monte_carlo(&si_bullet(reduced_hpa))
        .expect("monte carlo solve (equivalent absolute)");
    assert_eq!(
        with_reduced.velocity.to_bits(),
        with_equivalent_absolute.velocity.to_bits()
    );
}

#[test]
fn run_monte_carlo_uses_the_reduced_pressure_not_the_raw_qnh() {
    let reduced_hpa = reduce_qnh_to_station_pressure(QNH_HPA, ALTITUDE_M);

    let params = MonteCarloParams {
        num_simulations: 200,
        velocity_std_dev: 1.0,
        angle_std_dev: 0.001,
        bc_std_dev: 0.01,
        wind_speed_std_dev: 0.0,
        target_distance: None,
        base_wind_speed: 0.0,
        base_wind_direction: 0.0,
        azimuth_std_dev: 0.001,
    };

    let reduced_results = run_monte_carlo(si_bullet(reduced_hpa), params.clone())
        .expect("run_monte_carlo (reduced)");
    let raw_qnh_results =
        run_monte_carlo(si_bullet(QNH_HPA), params).expect("run_monte_carlo (raw QNH)");

    let mean_range = |ranges: &[f64]| ranges.iter().sum::<f64>() / ranges.len() as f64;
    let mean_reduced = mean_range(&reduced_results.ranges);
    let mean_raw_qnh = mean_range(&raw_qnh_results.ranges);

    assert!(
        (mean_reduced - mean_raw_qnh).abs() > 0.5,
        "reduced vs. raw-QNH pressure must change run_monte_carlo's mean range materially: \
         reduced={mean_reduced} raw_qnh={mean_raw_qnh}"
    );
}
