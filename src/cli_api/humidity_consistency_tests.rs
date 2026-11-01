//! The humidity class guard (humidity-consistency train, plan C0b and C2): every solver family
//! flies `TrajectorySolver`'s humid air.
//!
//! `TrajectorySolver` is the reference. Its per-step air is [`air_at_with_gravity`] over its
//! station ([`StationAir`], humidity included) and its downrange zones. The test-only trace
//! [`air_trace`] records every air composition a path makes: the base it composed from, the
//! relative humidity it gave the sound speed, and the air it got. Each test runs a path with the
//! trace on and checks what it recorded.
//!
//! - **One air.** Every family, at RH 50 and 100, with and without humid zones, composes exactly
//!   the reference's air at every step, to the bit ([`one_air_report`]).
//! - **The dry control.** A deliberately dry path must FAIL that check, by more than 1e-4 in the
//!   sound speed, so no tolerance can ever hide a path that goes dry.
//! - **Per step.** At fixed states, the four per-step sites (`TrajectorySolver`,
//!   `compute_derivatives`, the RK4/RK45 kernel, the `fast_integrate` substep) each compose the
//!   air [`air_at_with_gravity`] gives at that state, over a grid of stations, RH, zones,
//!   heights, downrange and gravity.
//! - **Response.** The drag and the Magnus accelerations respond to humidity at the same state
//!   exactly as `TrajectorySolver`'s do, to 1e-12 relative: the Mach every site uses is the
//!   reference's.
//! - **Lean Monte Carlo.** Its reported Mach divides by the sound speed it flew.
//! - **Aerodynamic jump.** The fast path's launch offset is `TrajectorySolver`'s, from the same
//!   station air, at every RH. It is humidity-free by design: the published Miller T/P density
//!   correction is dry.

use super::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};
use crate::atmosphere::air_trace::{traced, AirTraceEntry};
use crate::atmosphere::{
    air_at_with_gravity, calculate_air_density_cimp, get_local_atmosphere_humid_with_gravity,
    AtmoSegment, AtmoSock, KernelAtmosphere, LocalAir, RelativeHumidity, StationAir,
    StationAtmospherePresence,
};
use crate::gravity::Gravity;
use crate::trajectory_integration::{
    build_inputs, compute_derivatives_vec_with_gravity, integrate_trajectory,
    try_integrate_trajectory, try_solve_trajectory_rust, TrajectoryParams,
};
use crate::wind::{WindSegment, WindSock};
use crate::DragModel;
use nalgebra::{Vector3, Vector6};

// =================================================================================================
// The load and the air
// =================================================================================================

const MASS_KG: f64 = 0.01134;
const DIAMETER_M: f64 = 0.00782;
const LENGTH_M: f64 = 0.031;
const MUZZLE_VELOCITY_MPS: f64 = 823.0;
const BC_G7: f64 = 0.243;
const TWIST_IN: f64 = 10.0;
const MUZZLE_ANGLE_RAD: f64 = 0.002;
const MOA_PER_RAD: f64 = 3437.7467707849;

/// The air at the shooter.
#[derive(Clone, Copy, Debug)]
struct Station {
    altitude_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
    humidity_percent: f64,
}

impl Station {
    const fn new(altitude_m: f64, temperature_c: f64, pressure_hpa: f64) -> Self {
        Self {
            altitude_m,
            temperature_c,
            pressure_hpa,
            humidity_percent: 0.0,
        }
    }

    fn at_humidity(self, humidity_percent: f64) -> Self {
        Self {
            humidity_percent,
            ..self
        }
    }

    /// The validated station air, exactly as `TrajectorySolver` resolves it for these values
    /// (with `StationAtmospherePresence::BOTH`, so they fly as given).
    fn air(self) -> StationAir {
        StationAir::new(
            self.altitude_m,
            self.temperature_c,
            self.pressure_hpa,
            RelativeHumidity::try_percent(self.humidity_percent).unwrap(),
        )
        .unwrap()
    }

    fn conditions(self) -> AtmosphericConditions {
        AtmosphericConditions {
            temperature: self.temperature_c,
            pressure: self.pressure_hpa,
            humidity: self.humidity_percent,
            altitude: self.altitude_m,
        }
    }
}

const SEA_LEVEL: Station = Station::new(0.0, 15.0, 1013.25);
const HOT: Station = Station::new(0.0, 35.0, 1013.25);
const HIGH: Station = Station::new(1500.0, 25.0, 830.0);

/// Two downrange zones whose temperature, pressure and RH (20 %, 80 %) differ from each other
/// and from the station, so a recorded base names its zone unambiguously.
fn humid_zones(station: Station) -> Vec<AtmoSegment> {
    vec![
        (
            station.temperature_c + 2.0,
            station.pressure_hpa - 3.0,
            20.0,
            300.0,
        ),
        (
            station.temperature_c - 3.0,
            station.pressure_hpa + 2.0,
            80.0,
            100_000.0,
        ),
    ]
}

/// One atmosphere a family is flown in: a station, optional zones and the gravity.
#[derive(Clone, Debug)]
struct Scenario {
    station: Station,
    zones: Option<Vec<AtmoSegment>>,
    gravity: Gravity,
}

impl Scenario {
    fn new(station: Station, zones: Option<Vec<AtmoSegment>>) -> Self {
        Self {
            station,
            zones,
            gravity: Gravity::EARTH,
        }
    }

    fn sock(&self) -> Option<AtmoSock> {
        self.zones.clone().map(AtmoSock::new)
    }

    /// The reference: `TrajectorySolver`'s air at one state.
    fn reference_air(&self, downrange_m: f64, altitude_m: f64) -> LocalAir {
        air_at_with_gravity(
            &KernelAtmosphere::Station(self.station.air()),
            self.sock().as_ref(),
            downrange_m,
            altitude_m,
            self.gravity,
        )
    }

    /// The atmosphere the kernel families are given: the station itself, humidity included.
    /// (Written first against the 0.49.0 API, this was the legacy `(alt, T, P, ratio)` tuple
    /// with the station's density ratio at its RH, the best a caller could do: the tuple had no
    /// slot for the RH, and every kernel family failed these tests.)
    fn kernel_atmosphere(&self) -> KernelAtmosphere {
        KernelAtmosphere::Station(self.station.air())
    }

    fn solver(&self, inputs: BallisticInputs) -> TrajectorySolver {
        let mut solver = TrajectorySolver::new_with_station_presence_and_gravity(
            inputs,
            WindConditions::default(),
            self.station.conditions(),
            StationAtmospherePresence::BOTH,
            self.gravity,
        );
        if let Some(zones) = &self.zones {
            solver.set_atmo_segments(zones.clone());
        }
        solver
    }
}

/// The station, then the station with humid zones, at each RH.
fn scenarios(humidities: &[f64]) -> Vec<Scenario> {
    let mut out = Vec::new();
    for station in [SEA_LEVEL, HOT, HIGH] {
        for &humidity in humidities {
            let station = station.at_humidity(humidity);
            out.push(Scenario::new(station, None));
            out.push(Scenario::new(station, Some(humid_zones(station))));
        }
    }
    out
}

/// A .308 175 gr-class load at the station, every effect off, SI and imperial fields in step.
fn load(station: Station) -> BallisticInputs {
    BallisticInputs {
        bc_value: BC_G7,
        bc_type: DragModel::G7,
        bullet_mass: MASS_KG,
        muzzle_velocity: MUZZLE_VELOCITY_MPS,
        bullet_diameter: DIAMETER_M,
        bullet_length: LENGTH_M,
        twist_rate: TWIST_IN,
        caliber_inches: DIAMETER_M / 0.0254,
        weight_grains: MASS_KG / crate::constants::GRAINS_TO_KG,
        muzzle_angle: MUZZLE_ANGLE_RAD,
        sight_height: 0.0,
        muzzle_height: 0.0,
        altitude: station.altitude_m,
        temperature: station.temperature_c,
        pressure: station.pressure_hpa,
        // BallisticInputs.humidity is a FRACTION.
        humidity: station.humidity_percent / 100.0,
        target_distance: 600.0,
        ..BallisticInputs::default()
    }
}

fn launch_state() -> [f64; 6] {
    [
        0.0,
        0.0,
        0.0,
        MUZZLE_VELOCITY_MPS * MUZZLE_ANGLE_RAD.cos(),
        MUZZLE_VELOCITY_MPS * MUZZLE_ANGLE_RAD.sin(),
        0.0,
    ]
}

fn trajectory_params(
    scenario: &Scenario,
    target_distance_m: f64,
    enable_magnus: bool,
) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: MASS_KG,
        bc: BC_G7,
        drag_model: DragModel::G7,
        wind_segments: Vec::new(),
        atmosphere: scenario.kernel_atmosphere(),
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus,
        enable_coriolis: false,
        target_distance_m,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: scenario.station.altitude_m,
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
        ground_threshold: -1000.0,
        atmo_sock: scenario.sock(),
    }
}

fn fast_params(scenario: &Scenario, horiz_m: f64) -> crate::fast_trajectory::FastIntegrationParams {
    crate::fast_trajectory::FastIntegrationParams {
        horiz: horiz_m,
        vert: 0.0,
        initial_state: launch_state(),
        t_span: (0.0, 30.0),
        atmosphere: scenario.kernel_atmosphere(),
        atmo_sock: scenario.sock(),
    }
}

// =================================================================================================
// The families, and what each composed
// =================================================================================================

macro_rules! solver_families {
    ($($(#[$doc:meta])* $family:ident),+ $(,)?) => {
        /// Every solver family that flies air. The list is the enum: adding a variant here adds
        /// it to [`SolverFamily::ALL`], and [`run_family`]'s exhaustive match then refuses to
        /// compile until it has a runner.
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        enum SolverFamily {
            $($(#[$doc])* $family),+
        }

        impl SolverFamily {
            const ALL: &'static [SolverFamily] = &[$(SolverFamily::$family),+];
        }
    };
}

solver_families!(
    /// `TrajectorySolver::solve`, Euler.
    TrajectorySolverEuler,
    /// `TrajectorySolver::solve`, RK4.
    TrajectorySolverRk4,
    /// `TrajectorySolver::solve`, RK45.
    TrajectorySolverRk45,
    /// A zero through `TrajectorySolver` (`calculate_zero_angle_with_station_presence`).
    TrajectorySolverZero,
    /// `derivatives::compute_derivatives` at several states.
    ComputeDerivatives,
    /// `trajectory_integration::integrate_trajectory`, RK4.
    IntegrateTrajectoryRk4,
    /// `trajectory_integration::integrate_trajectory`, RK45.
    IntegrateTrajectoryRk45,
    /// `trajectory_integration::try_integrate_trajectory`, RK4.
    TryIntegrateTrajectoryRk4,
    /// `trajectory_integration::try_integrate_trajectory`, RK45.
    TryIntegrateTrajectoryRk45,
    /// `trajectory_integration::try_solve_trajectory_rust` (no zones).
    TrySolveTrajectoryRust,
    /// `trajectory_integration::solve_trajectory_rust` (no zones).
    SolveTrajectoryRust,
    /// `fast_trajectory::fast_integrate`.
    FastIntegrate,
    /// `fast_trajectory::fast_integrate_with_segments`.
    FastIntegrateWithSegments,
    /// `monte_carlo::solve_trajectory_for_monte_carlo` (no zones; presence-blind).
    LeanMonteCarlo,
    /// `monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence` (no zones).
    LeanMonteCarloWithStationPresence,
);

/// Run `family` in `scenario` with the trace on and return every air composition it made.
/// `None` when the family takes no downrange zones and the scenario has some.
fn run_family(family: SolverFamily, scenario: &Scenario) -> Option<Vec<AirTraceEntry>> {
    let station = scenario.station;
    let solve = |use_rk4: bool, use_adaptive_rk45: bool| {
        let mut solver = scenario.solver(BallisticInputs {
            use_rk4,
            use_adaptive_rk45,
            ..load(station)
        });
        solver.set_max_range(600.0);
        traced(|| solver.solve().expect("reference solve")).1
    };
    let entries = match family {
        SolverFamily::TrajectorySolverEuler => solve(false, false),
        SolverFamily::TrajectorySolverRk4 => solve(true, false),
        SolverFamily::TrajectorySolverRk45 => solve(true, true),
        SolverFamily::TrajectorySolverZero => {
            if scenario.zones.is_some() {
                return None;
            }
            traced(|| {
                crate::calculate_zero_angle_with_station_presence(
                    load(station),
                    300.0,
                    0.0,
                    WindConditions::default(),
                    station.conditions(),
                    StationAtmospherePresence::BOTH,
                )
                .expect("zero")
            })
            .1
        }
        SolverFamily::ComputeDerivatives => {
            let inputs = load(station);
            let sock = scenario.sock();
            traced(|| {
                for (x, y) in [(0.0, 0.0), (99.0, 4.0), (450.0, 120.0), (700.0, -30.0)] {
                    crate::derivatives::compute_derivatives(
                        Vector3::new(x, y, 0.0),
                        Vector3::new(700.0, 2.0, 0.5),
                        &inputs,
                        Vector3::zeros(),
                        &scenario.kernel_atmosphere(),
                        BC_G7,
                        None,
                        0.0,
                        sock.as_ref(),
                    );
                }
            })
            .1
        }
        SolverFamily::IntegrateTrajectoryRk4 | SolverFamily::IntegrateTrajectoryRk45 => {
            let method = if family == SolverFamily::IntegrateTrajectoryRk4 {
                "RK4"
            } else {
                "RK45"
            };
            traced(|| {
                integrate_trajectory(
                    launch_state(),
                    (0.0, 2.0),
                    trajectory_params(scenario, 600.0, false),
                    method,
                    1e-6,
                    0.001,
                )
            })
            .1
        }
        SolverFamily::TryIntegrateTrajectoryRk4 | SolverFamily::TryIntegrateTrajectoryRk45 => {
            let method = if family == SolverFamily::TryIntegrateTrajectoryRk4 {
                "RK4"
            } else {
                "RK45"
            };
            traced(|| {
                try_integrate_trajectory(
                    launch_state(),
                    (0.0, 2.0),
                    trajectory_params(scenario, 600.0, false),
                    method,
                    1e-6,
                    0.001,
                )
                .expect("valid kernel input")
            })
            .1
        }
        SolverFamily::TrySolveTrajectoryRust | SolverFamily::SolveTrajectoryRust => {
            if scenario.zones.is_some() {
                return None;
            }
            let run = || {
                let arguments = (
                    launch_state(),
                    (0.0, 2.0),
                    MASS_KG,
                    BC_G7,
                    DragModel::G7,
                    Vec::<WindSegment>::new(),
                    scenario.kernel_atmosphere(),
                );
                if family == SolverFamily::TrySolveTrajectoryRust {
                    try_solve_trajectory_rust(
                        arguments.0,
                        arguments.1,
                        arguments.2,
                        arguments.3,
                        arguments.4,
                        arguments.5,
                        arguments.6,
                        None,
                        false,
                        false,
                        false,
                        "RK45".to_string(),
                        1e-6,
                        0.01,
                        600.0,
                    )
                    .expect("valid kernel input")
                } else {
                    crate::trajectory_integration::solve_trajectory_rust(
                        arguments.0,
                        arguments.1,
                        arguments.2,
                        arguments.3,
                        arguments.4,
                        arguments.5,
                        arguments.6,
                        None,
                        false,
                        false,
                        false,
                        "RK45".to_string(),
                        1e-6,
                        0.01,
                        600.0,
                    )
                }
            };
            traced(run).1
        }
        SolverFamily::FastIntegrate => {
            let inputs = load(station);
            traced(|| {
                crate::fast_trajectory::fast_integrate(
                    &inputs,
                    &WindSock::new(Vec::new()),
                    fast_params(scenario, 600.0),
                )
            })
            .1
        }
        SolverFamily::FastIntegrateWithSegments => {
            let inputs = load(station);
            traced(|| {
                crate::fast_trajectory::fast_integrate_with_segments(
                    &inputs,
                    Vec::new(),
                    fast_params(scenario, 600.0),
                )
            })
            .1
        }
        SolverFamily::LeanMonteCarlo | SolverFamily::LeanMonteCarloWithStationPresence => {
            if scenario.zones.is_some() {
                return None;
            }
            let inputs = load(station);
            traced(|| {
                if family == SolverFamily::LeanMonteCarlo {
                    crate::monte_carlo::solve_trajectory_for_monte_carlo(&inputs)
                } else {
                    crate::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence(
                        &inputs,
                        StationAtmospherePresence::BOTH,
                    )
                }
                .expect("lean Monte Carlo solve")
            })
            .1
        }
    };
    Some(entries)
}

/// A deliberately DRY path: the per-step lookup the kernels made on 0.49.0, which composed the
/// station's (or the zone's) base with the right density ratio but gave the sound speed 0 % RH.
/// It must FAIL [`one_air_report`].
fn run_dry_control(scenario: &Scenario) -> Vec<AirTraceEntry> {
    let station = scenario.station.air();
    let sock = scenario.sock();
    traced(|| {
        for (downrange_m, height_m) in [(0.0, 0.0), (99.0, 4.0), (450.0, 120.0), (700.0, -30.0)] {
            let (base_temp_c, base_press_hpa, base_ratio) = match &sock {
                Some(sock) => {
                    let (t, p, h) = sock.atmo_for_range(downrange_m);
                    (t, p, calculate_air_density_cimp(t, p, h) / 1.225)
                }
                None => (
                    station.temperature_c(),
                    station.pressure_hpa(),
                    station.density_ratio(),
                ),
            };
            get_local_atmosphere_humid_with_gravity(
                station.altitude_m() + height_m,
                station.altitude_m(),
                base_temp_c,
                base_press_hpa,
                base_ratio,
                RelativeHumidity::DRY.percent(),
                scenario.gravity,
            );
        }
    })
    .1
}

/// How a path's recorded air compares with the reference's.
#[derive(Debug, Default)]
struct OneAirReport {
    entries: usize,
    mismatches: usize,
    max_relative_sound_speed_error: f64,
    max_relative_density_error: f64,
    first_mismatch: Option<String>,
}

/// A base the reference can compose from: the station, or one zone, with a downrange inside it.
struct ReferenceBase {
    temperature_c: f64,
    pressure_hpa: f64,
    density_ratio: f64,
    humidity_percent: f64,
    downrange_m: f64,
}

fn reference_bases(scenario: &Scenario) -> Vec<ReferenceBase> {
    let station = scenario.station.air();
    match &scenario.zones {
        None => vec![ReferenceBase {
            temperature_c: station.temperature_c(),
            pressure_hpa: station.pressure_hpa(),
            density_ratio: station.density_ratio(),
            humidity_percent: station.humidity().percent(),
            downrange_m: 0.0,
        }],
        Some(zones) => {
            let mut sorted = zones.clone();
            sorted.sort_by(|a, b| a.3.total_cmp(&b.3));
            let mut start_m = 0.0;
            sorted
                .iter()
                .map(
                    |&(temperature_c, pressure_hpa, humidity_percent, until_m)| {
                        let base = ReferenceBase {
                            temperature_c,
                            pressure_hpa,
                            density_ratio: calculate_air_density_cimp(
                                temperature_c,
                                pressure_hpa,
                                humidity_percent,
                            ) / 1.225,
                            humidity_percent,
                            downrange_m: start_m,
                        };
                        start_m = until_m;
                        base
                    },
                )
                .collect()
        }
    }
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    if actual == expected {
        0.0
    } else {
        ((actual - expected) / expected).abs()
    }
}

/// Check every recorded composition against the reference at the same state: the lapse starts
/// at the station, the base is the station's (or the active zone's) with ITS density ratio and
/// ITS relative humidity, and the air is `air_at`'s, to the bit.
fn one_air_report(scenario: &Scenario, entries: &[AirTraceEntry]) -> OneAirReport {
    let station = scenario.station.air();
    let bases = reference_bases(scenario);
    let mut report = OneAirReport {
        entries: entries.len(),
        ..OneAirReport::default()
    };
    for (index, entry) in entries.iter().enumerate() {
        let base = bases.iter().find(|base| {
            base.temperature_c.to_bits() == entry.base_temp_c.to_bits()
                && base.pressure_hpa.to_bits() == entry.base_press_hpa.to_bits()
        });
        let problem = match base {
            None => Some(format!(
                "entry {index}: base {} C / {} hPa is neither the station nor a zone",
                entry.base_temp_c, entry.base_press_hpa
            )),
            Some(base) => {
                let expected = scenario.reference_air(base.downrange_m, entry.altitude_m);
                report.max_relative_sound_speed_error = report.max_relative_sound_speed_error.max(
                    relative_error(entry.air.speed_of_sound_mps, expected.speed_of_sound_mps),
                );
                report.max_relative_density_error = report.max_relative_density_error.max(
                    relative_error(entry.air.density_kg_m3, expected.density_kg_m3),
                );
                let air_bits = |air: &LocalAir| {
                    [
                        air.density_kg_m3.to_bits(),
                        air.speed_of_sound_mps.to_bits(),
                        air.temperature_k.to_bits(),
                        air.pressure_pa.to_bits(),
                    ]
                };
                if entry.base_alt_m.to_bits() != station.altitude_m().to_bits() {
                    Some(format!(
                        "entry {index}: lapse from {} m, station at {} m",
                        entry.base_alt_m,
                        station.altitude_m()
                    ))
                } else if entry.gravity_mps2.to_bits() != scenario.gravity.magnitude().to_bits() {
                    Some(format!(
                        "entry {index}: gravity {}, scenario {}",
                        entry.gravity_mps2,
                        scenario.gravity.magnitude()
                    ))
                } else if entry.base_ratio.to_bits() != base.density_ratio.to_bits() {
                    Some(format!(
                        "entry {index}: density ratio {} where the reference has {}",
                        entry.base_ratio, base.density_ratio
                    ))
                } else if entry.humidity_percent.to_bits() != base.humidity_percent.to_bits() {
                    Some(format!(
                        "entry {index}: sound speed at {} % RH where the reference flies {} % \
                         ({} vs {} m/s)",
                        entry.humidity_percent,
                        base.humidity_percent,
                        entry.air.speed_of_sound_mps,
                        expected.speed_of_sound_mps
                    ))
                } else if air_bits(&entry.air) != air_bits(&expected) {
                    Some(format!(
                        "entry {index}: air {:?} where the reference has {expected:?}",
                        entry.air
                    ))
                } else {
                    None
                }
            }
        };
        if let Some(problem) = problem {
            report.mismatches += 1;
            report.first_mismatch.get_or_insert(problem);
        }
    }
    report
}

/// The humidities a family is held to the reference at, and every scenario it runs in.
const HUMID: [f64; 2] = [50.0, 100.0];

#[test]
fn every_family_flies_the_reference_air_at_every_step() {
    let mut failures = Vec::new();
    for &family in SolverFamily::ALL {
        let mut ran = 0;
        for scenario in scenarios(&HUMID) {
            let Some(entries) = run_family(family, &scenario) else {
                continue;
            };
            ran += 1;
            let report = one_air_report(&scenario, &entries);
            if report.entries == 0 || report.mismatches > 0 {
                failures.push(format!(
                    "{family:?} at {:?} {}: {} of {} steps differ (max relative c error {:.3e}, \
                     density {:.3e}); first: {}",
                    scenario.station,
                    if scenario.zones.is_some() {
                        "with zones RH 20/80"
                    } else {
                        "no zones"
                    },
                    report.mismatches,
                    report.entries,
                    report.max_relative_sound_speed_error,
                    report.max_relative_density_error,
                    report
                        .first_mismatch
                        .unwrap_or_else(|| "no air recorded".to_string())
                ));
            }
        }
        assert!(ran > 0, "{family:?} ran in no scenario");
    }
    assert!(
        failures.is_empty(),
        "{} family/scenario runs do not fly TrajectorySolver's air:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn every_family_gives_the_sound_speed_the_station_or_zone_humidity() {
    // The coarse form of the one-air test: the set of RH a family handed the sound speed must
    // be exactly the station's (no zones) or the zones' (with zones), never 0 % or anything else.
    let mut failures = Vec::new();
    for &family in SolverFamily::ALL {
        for scenario in scenarios(&HUMID) {
            let Some(entries) = run_family(family, &scenario) else {
                continue;
            };
            let mut used: Vec<f64> = entries.iter().map(|e| e.humidity_percent).collect();
            used.sort_by(f64::total_cmp);
            used.dedup();
            let allowed: Vec<f64> = match &scenario.zones {
                None => vec![scenario.station.humidity_percent],
                Some(zones) => zones.iter().map(|zone| zone.2).collect(),
            };
            if used.is_empty() || used.iter().any(|rh| !allowed.contains(rh)) {
                failures.push(format!(
                    "{family:?} at {:?} zones={}: sound speed given RH {used:?}, \
                     allowed {allowed:?}",
                    scenario.station,
                    scenario.zones.is_some()
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn the_kernel_atmosphere_wins_over_the_inputs_own_station_fields() {
    // A library caller's BallisticInputs may still carry another station (here the sea-level
    // defaults, dry). The kernel atmosphere is the air that flies, measured from ITS station's
    // altitude, exactly as TrajectorySolver flies its own: the inputs' altitude, temperature,
    // pressure and humidity must not leak into any step.
    for station in [HIGH.at_humidity(50.0), HOT.at_humidity(100.0)] {
        for zones in [None, Some(humid_zones(station))] {
            let scenario = Scenario::new(station, zones);
            let stray = BallisticInputs {
                altitude: 0.0,
                temperature: 15.0,
                pressure: 1013.25,
                humidity: 0.0,
                ..load(station)
            };
            let sock = scenario.sock();
            let runs = [
                (
                    "compute_derivatives",
                    traced(|| {
                        for (x, y) in [(0.0, 0.0), (99.0, 4.0), (450.0, 120.0)] {
                            crate::derivatives::compute_derivatives(
                                Vector3::new(x, y, 0.0),
                                Vector3::new(700.0, 2.0, 0.5),
                                &stray,
                                Vector3::zeros(),
                                &scenario.kernel_atmosphere(),
                                BC_G7,
                                None,
                                0.0,
                                sock.as_ref(),
                            );
                        }
                    })
                    .1,
                ),
                (
                    "fast_integrate",
                    traced(|| {
                        crate::fast_trajectory::fast_integrate(
                            &stray,
                            &WindSock::new(Vec::new()),
                            fast_params(&scenario, 600.0),
                        )
                    })
                    .1,
                ),
                (
                    "fast_integrate_with_segments",
                    traced(|| {
                        crate::fast_trajectory::fast_integrate_with_segments(
                            &stray,
                            Vec::new(),
                            fast_params(&scenario, 600.0),
                        )
                    })
                    .1,
                ),
            ];
            for (name, entries) in runs {
                let report = one_air_report(&scenario, &entries);
                assert!(
                    report.entries > 0 && report.mismatches == 0,
                    "{name} at {station:?} zones={}: {report:?}",
                    scenario.zones.is_some()
                );
                // The one-air check takes each step's altitude as given; the heights must be
                // measured from the station too. compute_derivatives ran at known heights; the
                // flights launch at the station and stay within a few metres of it.
                if name == "compute_derivatives" {
                    let heights: Vec<u64> =
                        entries.iter().map(|e| e.altitude_m.to_bits()).collect();
                    let expected: Vec<u64> = [0.0, 4.0, 120.0]
                        .iter()
                        .map(|y| (station.altitude_m + y).to_bits())
                        .collect();
                    assert_eq!(heights, expected, "{name} at {station:?}: step altitudes");
                } else {
                    assert_eq!(
                        entries[0].altitude_m.to_bits(),
                        station.altitude_m.to_bits(),
                        "{name} at {station:?}: launch altitude"
                    );
                    assert!(
                        entries
                            .iter()
                            .all(|e| (e.altitude_m - station.altitude_m).abs() < 50.0),
                        "{name} at {station:?}: a step left the station's altitude band"
                    );
                }
            }
        }
    }
}

#[test]
fn a_dry_control_path_fails_the_one_air_check_by_more_than_any_tolerance() {
    // The separation companion: a path that flies the right density but dry sound must be
    // caught, and by a margin (> 1e-4 relative in c) that no plausible tolerance would absorb.
    for scenario in scenarios(&HUMID) {
        let report = one_air_report(&scenario, &run_dry_control(&scenario));
        assert!(report.entries > 0);
        assert_eq!(
            report.mismatches, report.entries,
            "every dry step must be caught at {:?}",
            scenario.station
        );
        assert!(
            report.max_relative_sound_speed_error > 1e-4,
            "dry sound speed is only {:.3e} off at {:?}",
            report.max_relative_sound_speed_error,
            scenario.station
        );
    }
}

#[test]
fn every_family_is_dry_consistent_at_zero_humidity_too() {
    // At RH 0 every family was already TrajectorySolver's air (the S1 pins hold that to the
    // bit); the one-air check must agree, or it is checking the wrong thing.
    for &family in SolverFamily::ALL {
        let scenario = Scenario::new(SEA_LEVEL, None);
        let entries = run_family(family, &scenario).expect("no zones");
        let report = one_air_report(&scenario, &entries);
        assert!(
            report.entries > 0 && report.mismatches == 0,
            "{family:?} at RH 0: {report:?}"
        );
    }
}

// =================================================================================================
// Per step: the four per-step sites at fixed states
// =================================================================================================

/// The four places a per-step acceleration is computed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PerStepSite {
    /// `TrajectorySolver::calculate_acceleration`, the reference.
    TrajectorySolver,
    /// `derivatives::compute_derivatives_with_gravity`.
    ComputeDerivatives,
    /// `trajectory_integration::compute_derivatives_vec_with_gravity` behind `build_inputs`.
    TrajectoryIntegration,
    /// The `fast_integrate` substep.
    FastIntegrate,
}

const SITES: [PerStepSite; 4] = [
    PerStepSite::TrajectorySolver,
    PerStepSite::ComputeDerivatives,
    PerStepSite::TrajectoryIntegration,
    PerStepSite::FastIntegrate,
];

/// The acceleration `site` computes at one state, with the air compositions it made.
fn per_step(
    site: PerStepSite,
    scenario: &Scenario,
    inputs: &BallisticInputs,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
) -> (Vector3<f64>, Vec<AirTraceEntry>) {
    // Every site sees the inputs exactly as TrajectorySolver normalized them.
    let solver = scenario.solver(inputs.clone());
    let inputs = solver.inputs.clone();
    let sock = scenario.sock();
    match site {
        PerStepSite::TrajectorySolver => {
            let station = solver.station_air().unwrap();
            traced(|| {
                solver.calculate_acceleration(&position, &velocity, &Vector3::zeros(), &station)
            })
        }
        PerStepSite::ComputeDerivatives => {
            let (derivative, entries) = traced(|| {
                crate::derivatives::compute_derivatives_with_gravity(
                    position,
                    velocity,
                    &inputs,
                    Vector3::zeros(),
                    &scenario.kernel_atmosphere(),
                    inputs.bc_value,
                    None,
                    0.0,
                    sock.as_ref(),
                    scenario.gravity,
                )
            });
            (
                Vector3::new(derivative[3], derivative[4], derivative[5]),
                entries,
            )
        }
        PerStepSite::TrajectoryIntegration => {
            let params = trajectory_params(scenario, 1000.0, inputs.enable_magnus);
            let built = build_inputs(&params, inputs.muzzle_velocity);
            let state = Vector6::new(
                position.x, position.y, position.z, velocity.x, velocity.y, velocity.z,
            );
            let (derivative, entries) = traced(|| {
                compute_derivatives_vec_with_gravity(&state, 0.0, &params, &built, scenario.gravity)
            });
            (
                Vector3::new(derivative[3], derivative[4], derivative[5]),
                entries,
            )
        }
        PerStepSite::FastIntegrate => {
            // fast_integrate's own setup: the inputs carrying the station it was given.
            let station = scenario.station.air();
            let mut effective = inputs.clone();
            effective.altitude = station.altitude_m();
            effective.temperature = station.temperature_c();
            effective.pressure = station.pressure_hpa();
            let shape = crate::transonic_drag::resolve_projectile_shape(
                effective.bullet_model.as_deref(),
                effective.caliber_inches,
                effective.weight_grains,
                "G7",
            );
            let state = [
                position.x, position.y, position.z, velocity.x, velocity.y, velocity.z,
            ];
            let wind = WindSock::new(Vec::new());
            let (derivative, entries) = traced(|| {
                crate::fast_trajectory::compute_derivatives_with_gravity(
                    &state,
                    &effective,
                    &wind,
                    &scenario.kernel_atmosphere(),
                    &effective.bc_type,
                    shape,
                    effective.bc_value,
                    false,
                    false,
                    None,
                    None,
                    sock.as_ref(),
                    scenario.gravity,
                )
            });
            (
                Vector3::new(derivative[3], derivative[4], derivative[5]),
                entries,
            )
        }
    }
}

#[test]
fn every_per_step_site_composes_air_at_at_the_same_state() {
    let gravities = [Gravity::EARTH, Gravity::new(1.62).unwrap()];
    let mut failures = Vec::new();
    let mut checked = 0;
    for temperature_c in [-10.0, 15.0, 35.0] {
        for (altitude_m, pressure_hpa) in [(1500.0, 850.0), (0.0, 1013.25)] {
            for humidity in [0.0, 50.0, 100.0] {
                let station =
                    Station::new(altitude_m, temperature_c, pressure_hpa).at_humidity(humidity);
                for zones in [None, Some(humid_zones(station))] {
                    for gravity in gravities {
                        let scenario = Scenario {
                            station,
                            zones: zones.clone(),
                            gravity,
                        };
                        let inputs = load(station);
                        for height_m in [0.0, 50.0, 1000.0, -200.0] {
                            for downrange_m in [0.0, 99.0, 500.0] {
                                let position = Vector3::new(downrange_m, height_m, 0.3);
                                let velocity = Vector3::new(700.0, 3.0, 0.5);
                                let expected =
                                    scenario.reference_air(downrange_m, altitude_m + height_m);
                                for site in SITES {
                                    let (_, entries) =
                                        per_step(site, &scenario, &inputs, position, velocity);
                                    checked += 1;
                                    let matches = entries.len() == 1
                                        && entries[0].air.density_kg_m3.to_bits()
                                            == expected.density_kg_m3.to_bits()
                                        && entries[0].air.speed_of_sound_mps.to_bits()
                                            == expected.speed_of_sound_mps.to_bits()
                                        && entries[0].air.temperature_k.to_bits()
                                            == expected.temperature_k.to_bits();
                                    if !matches && failures.len() < 40 {
                                        failures.push(format!(
                                            "{site:?} {station:?} zones={} g={} h={height_m} \
                                             x={downrange_m}: {:?} vs air_at {expected:?}",
                                            zones.is_some(),
                                            gravity.magnitude(),
                                            entries.iter().map(|e| e.air).collect::<Vec<_>>()
                                        ));
                                    } else if !matches {
                                        failures.push(String::new());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(checked, 3 * 2 * 3 * 2 * 2 * 4 * 3 * 4);
    assert!(
        failures.is_empty(),
        "{} of {checked} per-step lookups differ from air_at; first ones:\n{}",
        failures.len(),
        failures
            .iter()
            .filter(|f| !f.is_empty())
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

// =================================================================================================
// Response: drag and Magnus respond to RH exactly as TrajectorySolver's do
// =================================================================================================

/// `site`'s aerodynamic acceleration (drag, plus Magnus when on) at one state: the total minus
/// the shot-frame gravity.
fn aerodynamic(
    site: PerStepSite,
    scenario: &Scenario,
    inputs: &BallisticInputs,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
) -> Vector3<f64> {
    per_step(site, scenario, inputs, position, velocity).0
        - scenario.gravity.in_frame(inputs.shooting_angle)
}

/// The humidity response `|a(RH)| / |a(0)|` of `site`'s drag (Magnus off) or its Magnus
/// acceleration alone (`a(Magnus on) - a(Magnus off)`).
fn humidity_response(
    site: PerStepSite,
    station: Station,
    humidity: f64,
    magnus: bool,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
) -> f64 {
    let acceleration = |humidity: f64| {
        let station = station.at_humidity(humidity);
        let scenario = Scenario::new(station, None);
        let off = load(station);
        let drag = aerodynamic(site, &scenario, &off, position, velocity);
        if !magnus {
            return drag;
        }
        let on = BallisticInputs {
            enable_magnus: true,
            use_enhanced_spin_drift: false,
            ..off
        };
        aerodynamic(site, &scenario, &on, position, velocity) - drag
    };
    acceleration(humidity).norm() / acceleration(0.0).norm()
}

#[test]
fn drag_and_magnus_respond_to_humidity_exactly_as_the_reference_does() {
    let mut failures = Vec::new();
    for station in [SEA_LEVEL, HOT] {
        for humidity in HUMID {
            // Supersonic (Magnus C_Np rises with Mach), transonic, and subsonic (C_Np flat).
            for speed in [700.0, 400.0, 250.0] {
                for height_m in [0.0, 500.0] {
                    let position = Vector3::new(100.0, height_m, 0.0);
                    let velocity = Vector3::new(speed, 1.0, 0.2);
                    for magnus in [false, true] {
                        let reference = humidity_response(
                            PerStepSite::TrajectorySolver,
                            station,
                            humidity,
                            magnus,
                            position,
                            velocity,
                        );
                        // Non-vacuity: the reference responds, except subsonic Magnus, whose C_Np
                        // is flat in Mach and whose force is density-free (the yaw of repose
                        // goes as 1/rho), so it barely moves with RH at all.
                        if !(magnus && speed < 300.0) {
                            assert!(
                                (reference - 1.0).abs() > 1e-6,
                                "the reference must respond: {reference} ({} v={speed})",
                                if magnus { "Magnus" } else { "drag" }
                            );
                        }
                        for site in SITES {
                            let response = humidity_response(
                                site, station, humidity, magnus, position, velocity,
                            );
                            if (response - reference).abs() > 1e-12 * reference {
                                failures.push(format!(
                                    "{site:?} {} {station:?} RH {humidity} v={speed} \
                                     h={height_m}: response {response} vs reference {reference}",
                                    if magnus { "Magnus" } else { "drag" }
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// =================================================================================================
// The lean Monte Carlo reports Mach against the sound speed it flew (H13)
// =================================================================================================

#[test]
fn lean_monte_carlo_mach_divides_by_the_sound_speed_it_flew() {
    for station in [SEA_LEVEL, HOT, HIGH] {
        for humidity in [0.0, 50.0, 100.0] {
            let station = station.at_humidity(humidity);
            let inputs = load(station);
            let (output, entries) = traced(|| {
                crate::monte_carlo::solve_trajectory_for_monte_carlo_with_station_presence(
                    &inputs,
                    StationAtmospherePresence::BOTH,
                )
                .expect("lean Monte Carlo solve")
            });
            // The first lookup is the launch state, at the station itself.
            let launch = entries.first().expect("the lean Monte Carlo flew some air");
            assert_eq!(launch.altitude_m.to_bits(), station.altitude_m.to_bits());
            let flown = launch.air.speed_of_sound_mps;
            assert_eq!(
                flown.to_bits(),
                station.air().speed_of_sound_mps().to_bits(),
                "{station:?}: flew c = {flown} at the muzzle, station c = {}",
                station.air().speed_of_sound_mps()
            );
            assert_eq!(
                output.mach.to_bits(),
                (output.velocity / flown).to_bits(),
                "{station:?}: reported Mach {} is not speed / flown c {}",
                output.mach,
                output.velocity / flown
            );
        }
    }
}

// =================================================================================================
// Aerodynamic jump: the fast path's launch offset is TrajectorySolver's, from the same station
// =================================================================================================

#[test]
fn fast_path_aerodynamic_jump_is_the_reference_jump_from_the_same_station() {
    for station in [SEA_LEVEL, HOT, HIGH] {
        let mut offsets = Vec::new();
        for humidity in [0.0, 50.0, 100.0] {
            let station = station.at_humidity(humidity);
            let scenario = Scenario::new(station, None);
            let inputs = BallisticInputs {
                enable_aerodynamic_jump: true,
                wind_speed: 4.0,
                wind_angle: std::f64::consts::FRAC_PI_2,
                ..load(station)
            };
            let mut solver = TrajectorySolver::new_with_station_presence(
                inputs.clone(),
                WindConditions {
                    speed: 4.0,
                    direction: std::f64::consts::FRAC_PI_2,
                    vertical_speed: 0.0,
                },
                station.conditions(),
                StationAtmospherePresence::BOTH,
            );
            solver.set_max_range(600.0);
            let reference = solver
                .aerodynamic_jump_components()
                .expect("jump components")
                .vertical_jump_moa
                / MOA_PER_RAD;
            let fast = crate::fast_trajectory::aerodynamic_jump_launch_offset_rad(
                &solver.inputs,
                &scenario.kernel_atmosphere(),
            );
            assert_eq!(
                fast.to_bits(),
                reference.to_bits(),
                "{station:?}: fast {fast} vs TrajectorySolver {reference}"
            );
            offsets.push(fast.to_bits());
        }
        // The published Miller density correction is dry (T/P only), by decision: the jump does
        // not move with RH on either path.
        offsets.dedup();
        assert_eq!(offsets.len(), 1, "{station:?}: jump moved with RH");
    }
}
