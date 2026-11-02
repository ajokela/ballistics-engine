//! What humid air costs a flight (R4, review HR-2), in pairs: `fast_integrate`, the lean Monte
//! Carlo trial and `TrajectorySolver` RK4 and RK45, each at 0 % and 50 % relative humidity.
//!
//! Every kernel looks up its air at every step through `air_at`. The governing air's
//! water-vapour mole fraction, and a downrange zone's density ratio, are fixed for the flight,
//! so the engine computes them once, when the station or the zones are built; a step evaluates
//! no saturation pressure. Each kernel is one criterion group with its dry and humid runs side
//! by side, so the pair shows directly that a real humidity adds no work, and a kernel that
//! regressed to per-step vapour work shows up against a saved baseline.
//!
//! Compare two commits with a criterion baseline (release profile, nothing else running):
//!
//! ```text
//! git switch <before> && cargo bench --bench humidity_air_bench -- --save-baseline before
//! git switch <after>  && cargo bench --bench humidity_air_bench -- --baseline before
//! ```
//!
//! Every case flies one shot: 850 m/s, a G1 BC of 0.5, 1000 m flat, from a sea-level station
//! at 30 °C and 1013.25 hPa (`BallisticInputs::default()` otherwise). Each iteration builds its
//! own station, so the per-flight work is counted too.

use ballistics_engine::fast_trajectory::{fast_integrate, FastIntegrationParams};
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo;
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, KernelAtmosphere, RelativeHumidity, StationAir,
    TrajectorySolver, WindConditions,
};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

/// The relative humidities of each pair, percent: dry and the default.
const HUMIDITIES: [f64; 2] = [0.0, 50.0];

const MUZZLE_VELOCITY_MPS: f64 = 850.0;
const RANGE_M: f64 = 1000.0;
const TEMPERATURE_C: f64 = 30.0;
const PRESSURE_HPA: f64 = 1013.25;

fn inputs(humidity_percent: f64, use_rk4: bool, use_adaptive_rk45: bool) -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bc_value: 0.5,
        target_distance: RANGE_M,
        temperature: TEMPERATURE_C,
        pressure: PRESSURE_HPA,
        altitude: 0.0,
        humidity_percent,
        use_rk4,
        use_adaptive_rk45,
        ..BallisticInputs::default()
    }
}

/// `fast_integrate` over its own station, built each flight.
fn fast_integrate_flight(inputs: &BallisticInputs, humidity_percent: f64) -> f64 {
    let station = StationAir::new(
        0.0,
        TEMPERATURE_C,
        PRESSURE_HPA,
        RelativeHumidity::try_percent(humidity_percent).expect("a valid humidity"),
    )
    .expect("air that exists");
    let solution = fast_integrate(
        inputs,
        &WindSock::new(Vec::new()),
        FastIntegrationParams {
            horiz: RANGE_M,
            vert: 0.0,
            initial_state: [0.0, 0.0, 0.0, MUZZLE_VELOCITY_MPS, 0.0, 0.0],
            t_span: (0.0, 10.0),
            atmosphere: KernelAtmosphere::Station(station),
            atmo_sock: None,
        },
    );
    assert!(solution.success, "fast_integrate flew");
    solution.t.last().copied().unwrap_or_default()
}

/// One lean Monte Carlo trial (the path an MC of N trials runs N times).
fn lean_monte_carlo_trial(inputs: &BallisticInputs) -> f64 {
    let trial = solve_trajectory_for_monte_carlo(inputs).expect("the trial flew");
    trial.drop + trial.time + trial.velocity
}

/// `TrajectorySolver::solve` with the integrator `inputs` selects.
fn solver_flight(inputs: BallisticInputs, humidity_percent: f64) -> f64 {
    let atmosphere = AtmosphericConditions {
        temperature: TEMPERATURE_C,
        pressure: PRESSURE_HPA,
        humidity: humidity_percent,
        altitude: 0.0,
    };
    let mut solver = TrajectorySolver::new(inputs, WindConditions::default(), atmosphere);
    solver.set_max_range(RANGE_M + 100.0);
    let result = solver.solve().expect("the solver flew");
    result.max_range + result.time_of_flight + result.impact_velocity
}

fn bench_fast_integrate(c: &mut Criterion) {
    let mut group = c.benchmark_group("humid_air/fast_integrate");
    for humidity in HUMIDITIES {
        let inputs = inputs(humidity, true, false);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("rh{humidity}")),
            &humidity,
            |b, &humidity| b.iter(|| fast_integrate_flight(black_box(&inputs), humidity)),
        );
    }
    group.finish();
}

fn bench_lean_monte_carlo(c: &mut Criterion) {
    let mut group = c.benchmark_group("humid_air/lean_mc_trial");
    for humidity in HUMIDITIES {
        let inputs = inputs(humidity, true, true);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("rh{humidity}")),
            &humidity,
            |b, _| b.iter(|| lean_monte_carlo_trial(black_box(&inputs))),
        );
    }
    group.finish();
}

fn bench_solver(c: &mut Criterion) {
    for (name, use_adaptive_rk45) in [("solver_rk4", false), ("solver_rk45", true)] {
        let mut group = c.benchmark_group(format!("humid_air/{name}"));
        for humidity in HUMIDITIES {
            let inputs = inputs(humidity, true, use_adaptive_rk45);
            group.bench_with_input(
                BenchmarkId::from_parameter(format!("rh{humidity}")),
                &humidity,
                |b, &humidity| b.iter(|| solver_flight(black_box(inputs.clone()), humidity)),
            );
        }
        group.finish();
    }
}

criterion_group!(
    benches,
    bench_fast_integrate,
    bench_lean_monte_carlo,
    bench_solver
);
criterion_main!(benches);
