//! Humidity consistency through the PUBLIC API (humidity-consistency train, plan C0b and C2):
//! every kernel family responds to relative humidity as `TrajectorySolver` does.
//!
//! The in-crate class guard (`src/cli_api/humidity_consistency_tests.rs`) proves every family
//! composes the reference's air at every step. This file checks what a library caller sees:
//!
//! - **Trajectory parity.** `try_integrate_trajectory` RK4 and `TrajectorySolver` RK4, at the
//!   same fixed step, every effect off and the same launch state, agree to round-off at RH 0, 50
//!   and 100, on a flat shot and on a 30-degree lob.
//! - **Endpoint response.** Each family's RH response at a fixed downrange,
//!   `R_F = endpoint_F(RH) - endpoint_F(0)` in height and in speed, must be at least half the
//!   reference's (the sensitivity floor) and within 2 % of it (integrator truncation cancels in
//!   the response). A companion shows the response is far above the families' integrator
//!   differences, so neither check is vacuous. `TrajectorySolver` Monte Carlo, whose trial
//!   outcome is a deviation from its baseline rather than a point at the target, is held to the
//!   solver itself instead: a zero-spread trial is the solver's flight bit for bit, dry and
//!   humid.
//! - **Zones.** One zone equal to the station, and an empty zone set, are bitwise the station;
//!   a humid zone over a dry station moves every zoned family as it moves the reference.

use ballistics_engine::atmosphere::{AtmoSegment, AtmoSock};
use ballistics_engine::derivatives::compute_derivatives;
use ballistics_engine::fast_trajectory::{
    fast_integrate, fast_integrate_with_segments, FastIntegrationParams, FastSolution,
};
use ballistics_engine::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence;
use ballistics_engine::trajectory_integration::{
    try_integrate_trajectory, try_solve_trajectory_rust, TrajectoryParams,
};
use ballistics_engine::wind::WindSock;
use ballistics_engine::{
    run_monte_carlo_with_station_presence, AtmosphericConditions, BallisticInputs, DragModel,
    KernelAtmosphere, MonteCarloParams, RelativeHumidity, StationAir, StationAtmospherePresence,
    TrajectoryResult, TrajectorySolver, WindConditions,
};
use nalgebra::Vector3;

// =================================================================================================
// The load and the air
// =================================================================================================

const MASS_KG: f64 = 0.01134;
const DIAMETER_M: f64 = 0.00782;
const LENGTH_M: f64 = 0.031;
const MUZZLE_VELOCITY_MPS: f64 = 823.0;
const BC_G7: f64 = 0.243;
const TWIST_IN: f64 = 10.0;
const FLAT_ANGLE_RAD: f64 = 0.002;
const LOB_ANGLE_RAD: f64 = std::f64::consts::FRAC_PI_6;
const GRAINS_TO_KG: f64 = 0.00006479891;

#[derive(Clone, Copy, Debug)]
struct Station {
    altitude_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
}

const SEA_LEVEL: Station = Station {
    altitude_m: 0.0,
    temperature_c: 15.0,
    pressure_hpa: 1013.25,
};
const HOT: Station = Station {
    altitude_m: 0.0,
    temperature_c: 35.0,
    pressure_hpa: 1013.25,
};
const HIGH: Station = Station {
    altitude_m: 1500.0,
    temperature_c: 25.0,
    pressure_hpa: 830.0,
};
/// The `tests/atmo_sock_zones.rs` station, at which a zone of the same T/P/RH reproduces the
/// station's density ratio to the bit.
const ZONE_STATION: Station = Station {
    altitude_m: 0.0,
    temperature_c: 20.0,
    pressure_hpa: 1000.0,
};

impl Station {
    fn conditions(self, humidity_percent: f64) -> AtmosphericConditions {
        AtmosphericConditions {
            temperature: self.temperature_c,
            pressure: self.pressure_hpa,
            humidity: humidity_percent,
            altitude: self.altitude_m,
        }
    }

    /// What the kernel families are given for this station at `humidity_percent`: the station
    /// itself, humidity included. (Written first against the 0.49.0 API, this was the legacy
    /// `(alt, T, P, ratio)` tuple carrying the station's density ratio at the RH, the best a
    /// caller could do: the tuple had no slot for the RH, and every kernel family failed.)
    fn kernel_atmosphere(self, humidity_percent: f64) -> KernelAtmosphere {
        KernelAtmosphere::Station(
            StationAir::new(
                self.altitude_m,
                self.temperature_c,
                self.pressure_hpa,
                RelativeHumidity::try_percent(humidity_percent).unwrap(),
            )
            .unwrap(),
        )
    }
}

fn load(station: Station, muzzle_angle_rad: f64, humidity_percent: f64) -> BallisticInputs {
    BallisticInputs {
        bc_value: BC_G7,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: TWIST_IN,
        caliber_inches: DIAMETER_M / 0.0254,
        weight_grains: MASS_KG / GRAINS_TO_KG,
        muzzle_angle: muzzle_angle_rad,
        sight_height: 0.0,
        muzzle_height: 0.0,
        altitude: station.altitude_m,
        temperature: station.temperature_c,
        pressure: station.pressure_hpa,
        humidity_percent,
        ground_threshold: -10_000.0,
        ..BallisticInputs::default()
    }
}

fn launch_state(muzzle_angle_rad: f64) -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * muzzle_angle_rad.cos(),
        MUZZLE_VELOCITY_MPS * muzzle_angle_rad.sin(),
        0.0,
    ]
}

fn trajectory_params(
    station: Station,
    humidity_percent: f64,
    zones: Option<Vec<AtmoSegment>>,
    target_distance_m: f64,
) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: MASS_KG,
        bc: BC_G7,
        drag_model: DragModel::G7,
        wind_segments: Vec::new(),
        atmosphere: station.kernel_atmosphere(humidity_percent),
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: station.altitude_m,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: TWIST_IN,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -10_000.0,
        atmo_sock: zones.map(AtmoSock::new),
    }
}

fn fast_params(
    station: Station,
    humidity_percent: f64,
    zones: Option<Vec<AtmoSegment>>,
    horiz_m: f64,
) -> FastIntegrationParams {
    FastIntegrationParams {
        horiz: horiz_m,
        vert: 0.0,
        initial_state: launch_state(FLAT_ANGLE_RAD),
        t_span: (0.0, 30.0),
        atmosphere: station.kernel_atmosphere(humidity_percent),
        atmo_sock: zones.map(AtmoSock::new),
    }
}

fn solver(
    station: Station,
    humidity_percent: f64,
    zones: Option<Vec<AtmoSegment>>,
    inputs: BallisticInputs,
    max_range_m: f64,
) -> TrajectorySolver {
    let mut solver = TrajectorySolver::new_with_station_presence(
        inputs,
        WindConditions::default(),
        station.conditions(humidity_percent),
        StationAtmospherePresence::BOTH,
    );
    solver.set_max_range(max_range_m);
    if let Some(zones) = zones {
        solver.set_atmo_segments(zones);
    }
    solver
}

// =================================================================================================
// Trajectory parity: kernel RK4 vs TrajectorySolver RK4 at a matched fixed step
// =================================================================================================

/// `(time, x, y, z, speed)` at every RK4 step of `TrajectorySolver`, before its terminal point.
fn solver_rk4_steps(station: Station, humidity_percent: f64, angle: f64) -> Vec<[f64; 5]> {
    let inputs = BallisticInputs {
        use_rk4: true,
        use_adaptive_rk45: false,
        ..load(station, angle, humidity_percent)
    };
    let result: TrajectoryResult = solver(station, humidity_percent, None, inputs, 1.0e6)
        .solve()
        .expect("reference solve");
    result
        .points
        .iter()
        .map(|p| {
            [
                p.time,
                p.position.x,
                p.position.y,
                p.position.z,
                p.velocity_magnitude,
            ]
        })
        .collect()
}

fn kernel_rk4_steps(
    station: Station,
    humidity_percent: f64,
    angle: f64,
    seconds: f64,
) -> Vec<[f64; 5]> {
    try_integrate_trajectory(
        launch_state(angle),
        (0.0, seconds),
        trajectory_params(station, humidity_percent, None, 1.0e6),
        "RK4",
        1e-6,
        0.001,
    )
    .expect("valid kernel input")
    .iter()
    .map(|(t, s)| [*t, s[0], s[1], s[2], Vector3::new(s[3], s[4], s[5]).norm()])
    .collect()
}

#[test]
fn kernel_rk4_matches_trajectory_solver_rk4_to_round_off_at_every_humidity() {
    let mut failures = Vec::new();
    for station in [SEA_LEVEL, HIGH] {
        for (angle, seconds) in [(FLAT_ANGLE_RAD, 2.0), (LOB_ANGLE_RAD, 20.0)] {
            for humidity in [0.0, 50.0, 100.0] {
                let reference = solver_rk4_steps(station, humidity, angle);
                let kernel = kernel_rk4_steps(station, humidity, angle, seconds);
                // Both step at 1 ms from the same state; compare every step both took, short
                // of the kernel's final (shortened) step.
                let steps = (kernel.len() - 1).min(reference.len() - 1);
                assert!(steps as f64 > 0.9 * seconds / 0.001, "{steps} steps");
                let mut worst = (0.0_f64, 0usize, 0usize);
                for step in 0..steps {
                    for component in 0..5 {
                        let (a, b) = (kernel[step][component], reference[step][component]);
                        let error = (a - b).abs() / b.abs().max(1.0);
                        if error > worst.0 {
                            worst = (error, step, component);
                        }
                    }
                }
                if worst.0 > 1e-12 {
                    let (_, step, component) = worst;
                    failures.push(format!(
                        "{station:?} angle {angle} RH {humidity}: step {step} component \
                         {component}: kernel {} vs TrajectorySolver {} (relative {:.3e}; the \
                         height gap there is {:.4} m)",
                        kernel[step][component],
                        reference[step][component],
                        worst.0,
                        kernel[step][2] - reference[step][2]
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// =================================================================================================
// Endpoint response: every family responds to RH as TrajectorySolver does
// =================================================================================================

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Family {
    /// The reference: `TrajectorySolver`, fixed-step RK4.
    TrajectorySolverRk4,
    TrajectorySolverRk45,
    TryIntegrateTrajectoryRk4,
    TryIntegrateTrajectoryRk45,
    TrySolveTrajectoryRust,
    FastIntegrate,
    FastIntegrateWithSegments,
    LeanMonteCarlo,
}

const FAMILIES: [Family; 8] = [
    Family::TrajectorySolverRk4,
    Family::TrajectorySolverRk45,
    Family::TryIntegrateTrajectoryRk4,
    Family::TryIntegrateTrajectoryRk45,
    Family::TrySolveTrajectoryRust,
    Family::FastIntegrate,
    Family::FastIntegrateWithSegments,
    Family::LeanMonteCarlo,
];

/// Height and speed where the flight crosses the target plane.
#[derive(Clone, Copy, Debug)]
struct Endpoint {
    height_m: f64,
    speed_mps: f64,
}

fn fast_endpoint(solution: &FastSolution) -> Endpoint {
    assert!(solution.success);
    let last = solution.t.len() - 1;
    Endpoint {
        height_m: solution.y[1][last],
        speed_mps: Vector3::new(
            solution.y[3][last],
            solution.y[4][last],
            solution.y[5][last],
        )
        .norm(),
    }
}

fn endpoint(family: Family, station: Station, humidity: f64, target_m: f64) -> Endpoint {
    let inputs = load(station, FLAT_ANGLE_RAD, humidity);
    let kernel = |method: &str| {
        let states = try_integrate_trajectory(
            launch_state(FLAT_ANGLE_RAD),
            (0.0, 30.0),
            trajectory_params(station, humidity, None, target_m),
            method,
            1e-6,
            if method == "RK4" { 0.001 } else { 0.01 },
        )
        .expect("valid kernel input");
        let (_, last) = states.last().unwrap();
        assert!(
            (last[0] - target_m).abs() < 1e-9,
            "{family:?} stopped short"
        );
        Endpoint {
            height_m: last[1],
            speed_mps: Vector3::new(last[3], last[4], last[5]).norm(),
        }
    };
    match family {
        Family::TrajectorySolverRk4 | Family::TrajectorySolverRk45 => {
            let inputs = BallisticInputs {
                use_rk4: true,
                use_adaptive_rk45: family == Family::TrajectorySolverRk45,
                ..inputs
            };
            let result = solver(station, humidity, None, inputs, target_m)
                .solve()
                .expect("reference solve");
            let last = result.points.last().unwrap();
            assert!((last.position.x - target_m).abs() < 1e-9);
            Endpoint {
                height_m: last.position.y,
                speed_mps: last.velocity_magnitude,
            }
        }
        Family::TryIntegrateTrajectoryRk4 => kernel("RK4"),
        Family::TryIntegrateTrajectoryRk45 => kernel("RK45"),
        Family::TrySolveTrajectoryRust => {
            let points = try_solve_trajectory_rust(
                launch_state(FLAT_ANGLE_RAD),
                (0.0, 30.0),
                MASS_KG,
                BC_G7,
                DragModel::G7,
                Vec::new(),
                station.kernel_atmosphere(humidity),
                None,
                false,
                false,
                false,
                "RK45".to_string(),
                1e-6,
                0.01,
                target_m,
            )
            .expect("valid kernel input");
            let last = points.last().unwrap();
            assert!((last["x"] - target_m).abs() < 1e-9);
            Endpoint {
                height_m: last["y"],
                speed_mps: Vector3::new(last["vx"], last["vy"], last["vz"]).norm(),
            }
        }
        Family::FastIntegrate => fast_endpoint(&fast_integrate(
            &inputs,
            &WindSock::new(Vec::new()),
            fast_params(station, humidity, None, target_m),
        )),
        Family::FastIntegrateWithSegments => {
            let solution = fast_integrate_with_segments(
                &inputs,
                Vec::new(),
                fast_params(station, humidity, None, target_m),
            );
            // The RK45 path saves ~50 range-spaced states; the target crossing is the last.
            fast_endpoint(&solution)
        }
        Family::LeanMonteCarlo => {
            let output = solve_trajectory_for_monte_carlo_with_station_presence(
                &BallisticInputs {
                    target_distance: target_m,
                    ..inputs
                },
                StationAtmospherePresence::BOTH,
            )
            .expect("lean Monte Carlo solve");
            Endpoint {
                // The line of sight is at the bore here, so the drop is minus the height.
                height_m: -output.drop,
                speed_mps: output.velocity,
            }
        }
    }
}

/// A flat shot to 1000 m, and the same load to 1200 m, where it is transonic.
const TARGETS_M: [f64; 2] = [1000.0, 1200.0];
/// RH 50 and 100 at 15 C, and 90 % at 35 C, where humid air is furthest from dry.
const HUMID_CASES: [(Station, f64); 3] = [(SEA_LEVEL, 50.0), (SEA_LEVEL, 100.0), (HOT, 90.0)];

#[test]
fn every_family_responds_to_humidity_as_trajectory_solver_does() {
    let mut failures = Vec::new();
    for target_m in TARGETS_M {
        for (station, humidity) in HUMID_CASES {
            let reference_dry = endpoint(Family::TrajectorySolverRk4, station, 0.0, target_m);
            let reference_humid =
                endpoint(Family::TrajectorySolverRk4, station, humidity, target_m);
            let reference = [
                reference_humid.height_m - reference_dry.height_m,
                reference_humid.speed_mps - reference_dry.speed_mps,
            ];
            for family in FAMILIES {
                let dry = endpoint(family, station, 0.0, target_m);
                let humid = endpoint(family, station, humidity, target_m);
                let response = [
                    humid.height_m - dry.height_m,
                    humid.speed_mps - dry.speed_mps,
                ];
                for (axis, (r, r_ref)) in ["height", "speed"]
                    .iter()
                    .zip(response.iter().zip(reference))
                {
                    let floor = r.abs() >= 0.5 * r_ref.abs();
                    let agrees = (r - r_ref).abs() <= 0.02 * r_ref.abs();
                    if !(floor && agrees) {
                        failures.push(format!(
                            "{family:?} {target_m} m {station:?} RH {humidity}: {axis} response \
                             {r:+.6} vs reference {r_ref:+.6} (ratio {:.4})",
                            r / r_ref
                        ));
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn the_humidity_response_dwarfs_the_families_integrator_differences() {
    // Non-vacuity for the response test: the reference's response must stand at least 10x
    // above how far the families' RH-0 endpoints already sit from the reference's (the
    // integrator differences that cancel in a response), and be a real distance and speed.
    for target_m in TARGETS_M {
        for (station, humidity) in HUMID_CASES {
            let dry = endpoint(Family::TrajectorySolverRk4, station, 0.0, target_m);
            let humid = endpoint(Family::TrajectorySolverRk4, station, humidity, target_m);
            let response = [
                (humid.height_m - dry.height_m).abs(),
                (humid.speed_mps - dry.speed_mps).abs(),
            ];
            assert!(
                response[0] > 0.005 && response[1] > 0.05,
                "{target_m} m {station:?} RH {humidity}: response {response:?} is too small to test"
            );
            for family in FAMILIES {
                let other = endpoint(family, station, 0.0, target_m);
                let gap = [
                    (other.height_m - dry.height_m).abs(),
                    (other.speed_mps - dry.speed_mps).abs(),
                ];
                for axis in 0..2 {
                    assert!(
                        response[axis] >= 10.0 * gap[axis],
                        "{family:?} {target_m} m {station:?} RH {humidity}: response \
                         {response:?} vs integrator gap {gap:?}"
                    );
                }
            }
        }
    }
}

/// `TrajectorySolver` Monte Carlo flies each trial through `TrajectorySolver` in the station's
/// air, humidity included. Its trial outcome is not a point at the target (that is a deviation
/// from the baseline, 0 with every spread zero), so instead of the endpoint above it is held to
/// the solver itself: one trial with every spread zero records exactly the range and impact
/// speed `TrajectorySolver` gives the same inputs, bit for bit, at RH 0 and humid, so its RH
/// response is the reference's. A Monte Carlo-only air (a per-batch cache, a trial built
/// without the humidity) fails here.
#[test]
fn trajectory_solver_monte_carlo_responds_to_humidity_as_its_solver_does() {
    for target_m in TARGETS_M {
        for (station, humidity) in HUMID_CASES {
            let mut responses = Vec::new();
            for rh in [0.0, humidity] {
                let inputs = BallisticInputs {
                    use_rk4: true,
                    ..load(station, FLAT_ANGLE_RAD, rh)
                };
                let params = MonteCarloParams {
                    num_simulations: 1,
                    velocity_std_dev: 0.0,
                    angle_std_dev: 0.0,
                    bc_std_dev: 0.0,
                    wind_speed_std_dev: 0.0,
                    target_distance: Some(target_m),
                    base_wind_speed: 0.0,
                    base_wind_direction: 0.0,
                    azimuth_std_dev: 0.0,
                };
                let monte_carlo = run_monte_carlo_with_station_presence(
                    inputs.clone(),
                    WindConditions::default(),
                    params,
                    0.0,
                    Some(0x5EED),
                    StationAtmospherePresence::BOTH,
                )
                .expect("Monte Carlo");
                // The sampler's own solve: the inputs' air, to twice the larger of the target
                // and 1000 m.
                let mut reference = TrajectorySolver::new_with_station_presence(
                    inputs,
                    WindConditions::default(),
                    station.conditions(rh),
                    StationAtmospherePresence::BOTH,
                );
                reference.set_max_range(target_m.max(1000.0) * 2.0);
                let reference = reference.solve().expect("reference solve");
                let what = format!("{target_m} m {station:?} RH {rh}");
                assert_eq!(monte_carlo.ranges.len(), 1, "{what}: the one trial flew");
                assert_eq!(
                    [
                        monte_carlo.ranges[0].to_bits(),
                        monte_carlo.impact_velocities[0].to_bits()
                    ],
                    [
                        reference.max_range.to_bits(),
                        reference.impact_velocity.to_bits()
                    ],
                    "{what}: Monte Carlo trial (range {}, speed {}) vs TrajectorySolver (range \
                     {}, speed {})",
                    monte_carlo.ranges[0],
                    monte_carlo.impact_velocities[0],
                    reference.max_range,
                    reference.impact_velocity
                );
                responses.push(monte_carlo.impact_velocities[0]);
            }
            assert!(
                (responses[1] - responses[0]).abs() > 0.05,
                "{target_m} m {station:?} RH {humidity}: the humidity moved the trial's impact \
                 speed by only {:+.6} m/s",
                responses[1] - responses[0]
            );
        }
    }
}

// =================================================================================================
// Zones
// =================================================================================================

/// The zoned families' endpoint (height, speed) at 1000 m, with `zones` over the station.
fn zoned_endpoints(
    station: Station,
    humidity: f64,
    zones: Option<Vec<AtmoSegment>>,
) -> Vec<(&'static str, [u64; 2])> {
    let target_m = 1000.0;
    let inputs = load(station, FLAT_ANGLE_RAD, humidity);
    let bits = |e: Endpoint| [e.height_m.to_bits(), e.speed_mps.to_bits()];
    let mut out = Vec::new();
    for (name, rk45) in [
        ("TrajectorySolver RK4", false),
        ("TrajectorySolver RK45", true),
    ] {
        let result = solver(
            station,
            humidity,
            zones.clone(),
            BallisticInputs {
                use_rk4: true,
                use_adaptive_rk45: rk45,
                ..inputs.clone()
            },
            target_m,
        )
        .solve()
        .expect("reference solve");
        let last = result.points.last().unwrap();
        out.push((
            name,
            bits(Endpoint {
                height_m: last.position.y,
                speed_mps: last.velocity_magnitude,
            }),
        ));
    }
    for method in ["RK4", "RK45"] {
        let states = try_integrate_trajectory(
            launch_state(FLAT_ANGLE_RAD),
            (0.0, 30.0),
            trajectory_params(station, humidity, zones.clone(), target_m),
            method,
            1e-6,
            if method == "RK4" { 0.001 } else { 0.01 },
        )
        .expect("valid kernel input");
        let (_, last) = states.last().unwrap();
        out.push((
            if method == "RK4" {
                "try_integrate_trajectory RK4"
            } else {
                "try_integrate_trajectory RK45"
            },
            bits(Endpoint {
                height_m: last[1],
                speed_mps: Vector3::new(last[3], last[4], last[5]).norm(),
            }),
        ));
    }
    out.push((
        "fast_integrate",
        bits(fast_endpoint(&fast_integrate(
            &inputs,
            &WindSock::new(Vec::new()),
            fast_params(station, humidity, zones.clone(), target_m),
        ))),
    ));
    out.push((
        "fast_integrate_with_segments",
        bits(fast_endpoint(&fast_integrate_with_segments(
            &inputs,
            Vec::new(),
            fast_params(station, humidity, zones.clone(), target_m),
        ))),
    ));
    // compute_derivatives at three states (its accelerations, not an endpoint).
    let sock = zones.map(AtmoSock::new);
    for (name, x, y) in [
        ("compute_derivatives at 0 m", 0.0, 0.0),
        ("compute_derivatives at 500 m", 500.0, 20.0),
        ("compute_derivatives at 900 m", 900.0, -40.0),
    ] {
        let d = compute_derivatives(
            Vector3::new(x, y, 0.0),
            Vector3::new(650.0, 1.0, 0.0),
            &inputs,
            Vector3::zeros(),
            &station.kernel_atmosphere(humidity),
            BC_G7,
            None,
            0.0,
            sock.as_ref(),
        );
        out.push((name, [d[3].to_bits(), d[4].to_bits()]));
    }
    out
}

#[test]
fn one_zone_equal_to_the_station_is_the_station_bit_for_bit() {
    for humidity in [0.0, 50.0, 100.0] {
        let station = ZONE_STATION;
        let zone = vec![(station.temperature_c, station.pressure_hpa, humidity, 1.0e6)];
        assert_eq!(
            zoned_endpoints(station, humidity, Some(zone)),
            zoned_endpoints(station, humidity, None),
            "RH {humidity}"
        );
    }
}

#[test]
fn an_empty_zone_set_is_the_station_bit_for_bit() {
    for (station, humidity) in [(ZONE_STATION, 50.0), (HIGH, 0.0), (HOT, 100.0)] {
        assert_eq!(
            zoned_endpoints(station, humidity, Some(Vec::new())),
            zoned_endpoints(station, humidity, None),
            "{station:?} RH {humidity}"
        );
    }
}

#[test]
fn a_humid_zone_moves_every_zoned_family_as_it_moves_the_reference() {
    // A zone at the station's T/P but RH 100, over a station at RH 0: every family's response
    // to the zone is the reference's within 2 %.
    let station = ZONE_STATION;
    let zone = vec![(station.temperature_c, station.pressure_hpa, 100.0, 1.0e6)];
    let target_m = 1000.0;
    let decode = |bits: [u64; 2]| [f64::from_bits(bits[0]), f64::from_bits(bits[1])];
    let dry = zoned_endpoints(station, 0.0, None);
    let zoned = zoned_endpoints(station, 0.0, Some(zone));
    let reference = {
        let (d, z) = (decode(dry[0].1), decode(zoned[0].1));
        [z[0] - d[0], z[1] - d[1]]
    };
    assert!(reference[0].abs() > 0.005 && reference[1].abs() > 0.05);
    let mut failures = Vec::new();
    // The compute_derivatives entries are accelerations, not endpoints; skip them here.
    for ((name, d), (_, z)) in dry
        .iter()
        .zip(&zoned)
        .filter(|((name, _), _)| !name.starts_with("compute_derivatives"))
    {
        let (d, z) = (decode(*d), decode(*z));
        for axis in 0..2 {
            let response = z[axis] - d[axis];
            if (response - reference[axis]).abs() > 0.02 * reference[axis].abs() {
                failures.push(format!(
                    "{name} at {target_m} m: zone response {response:+.6} vs reference {:+.6}",
                    reference[axis]
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// =================================================================================================
// The sensor-tolerance band (S5b, owner decision 2026-11-01)
// =================================================================================================

/// Readings within 5 points of the range, each with the end it reads as.
const SENSOR_READINGS: [(f64, f64); 6] = [
    (100.000_000_1, 100.0),
    (101.7, 100.0),
    (105.0, 100.0),
    (-0.0, 0.0),
    (-2.0, 0.0),
    (-5.0, 0.0),
];

#[test]
fn every_family_flies_a_sensor_reading_as_the_end_it_reads_as_bit_for_bit() {
    let bits = |e: Endpoint| [e.height_m.to_bits(), e.speed_mps.to_bits()];
    for (reading, end) in SENSOR_READINGS {
        for family in FAMILIES {
            assert_eq!(
                bits(endpoint(family, SEA_LEVEL, reading, 1000.0)),
                bits(endpoint(family, SEA_LEVEL, end, 1000.0)),
                "{family:?}: station RH {reading} must fly {end}"
            );
        }
        // The zoned families, the station and a zone each at the reading.
        let station = ZONE_STATION;
        let zone = |humidity: f64| {
            vec![
                (station.temperature_c, station.pressure_hpa, 50.0, 400.0),
                (30.0, 990.0, humidity, 1.0e6),
            ]
        };
        assert_eq!(
            zoned_endpoints(station, reading, Some(zone(reading))),
            zoned_endpoints(station, end, Some(zone(end))),
            "zoned families at RH {reading}"
        );
    }
}
