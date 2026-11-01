//! H5 builder pins (humidity-consistency train, plan C3): every builder of a `BallisticInputs`
//! leaves the relative humidity in `BallisticInputs.humidity_percent` as the PERCENT the rest of
//! the engine flies, exactly.
//!
//! 0.49.0's `BallisticInputs.humidity` was a 0-1 FRACTION that some builders filled with a
//! percent (read back through the clamping `humidity_percent()` method as 100 %) and others
//! round-tripped through `/ 100` and `* 100`. The field is now a percent, so a percent passes
//! through unconverted. These tests were written first against that API, where every builder
//! but the two fraction writers (solve v1, the WEZ echo) and the default failed.
//!
//! - **Probes.** 7 % and 55 % do not survive a fraction: 0.07 x 100 is 7.000000000000001.
//! - **The solver syncs.** Every `TrajectorySolver` constructor sets the inputs' humidity from
//!   its `AtmosphericConditions`, the air the solve flies.
//! - **Writers.** The hold-curve and DSF builders, the FFI Monte Carlo, solve v1, the kernel's
//!   derivative inputs and the fast path's station copy each write the percent.
//! - **Readers.** The Monte Carlo sampler, the lean Monte Carlo and the WEZ attribution each
//!   read it as a percent.
//!
//! The WASM and CLI writers are inline in their command handlers; each writes the same percent
//! variable into both structs, the solver sync makes the inputs' copy match regardless, and
//! `tests/humidity_guard.rs` refuses a hand conversion anywhere in `src/`.

use super::{
    AtmosphericConditions, BallisticInputs, MonteCarloParams, MonteCarloTrialSampler,
    TrajectorySolver, WindConditions,
};
use crate::atmosphere::air_trace::traced;
use crate::atmosphere::{
    KernelAtmosphere, RelativeHumidity, StationAir, StationAtmospherePresence, StationReading,
    DEFAULT_HUMIDITY_PERCENT,
};
use crate::trajectory_integration::TrajectoryParams;
use crate::DragModel;

/// 7 % and 55 % are the probes: 0.07 x 100 and 0.55 x 100 are not 7 and 55 in binary floating
/// point, so a builder that carries the humidity through a fraction anywhere cannot return them
/// exactly. 0, 30 and 100 are exact either way.
const HUMIDITIES: [f64; 5] = [0.0, 7.0, 30.0, 55.0, 100.0];

fn assert_percent(builder: &str, humidity: f64, got: f64) {
    assert_eq!(
        got.to_bits(),
        humidity.to_bits(),
        "{builder} at {humidity} %: BallisticInputs carries {got:?} %"
    );
}

/// A .308 load at a 1500 m, 25 C, 830 hPa station at `humidity` percent.
fn load(humidity: f64) -> BallisticInputs {
    BallisticInputs {
        bc_value: 0.243,
        bc_type: DragModel::G7,
        bullet_mass: 0.01134,
        muzzle_velocity: 823.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.031,
        muzzle_angle: 0.002,
        target_distance: 500.0,
        altitude: 1_500.0,
        temperature: 25.0,
        pressure: 830.0,
        humidity_percent: humidity,
        ..BallisticInputs::default()
    }
}

fn station(humidity: f64) -> StationAir {
    StationAir::new(
        1_500.0,
        25.0,
        830.0,
        RelativeHumidity::try_percent(humidity).unwrap(),
    )
    .unwrap()
}

#[test]
fn the_default_inputs_carry_the_default_percent() {
    assert_percent(
        "BallisticInputs::default",
        DEFAULT_HUMIDITY_PERCENT,
        BallisticInputs::default().humidity_percent,
    );
    assert_eq!(
        AtmosphericConditions::default().humidity.to_bits(),
        DEFAULT_HUMIDITY_PERCENT.to_bits()
    );
}

#[test]
fn every_trajectory_solver_constructor_syncs_the_inputs_humidity_from_its_atmosphere() {
    // The atmosphere is what the solver flies, so it wins; the inputs then say the same thing.
    // An out-of-range value is copied as it is: judging it is the solve's job.
    for humidity in HUMIDITIES.into_iter().chain([150.0, -5.0]) {
        let atmosphere = AtmosphericConditions {
            humidity,
            ..AtmosphericConditions::default()
        };
        let inputs = || load(80.0);
        let wind = WindConditions::default;
        let solvers = [
            (
                "new",
                TrajectorySolver::new(inputs(), wind(), atmosphere.clone()),
            ),
            (
                "new_with_station_presence(NEITHER)",
                TrajectorySolver::new_with_station_presence(
                    inputs(),
                    wind(),
                    atmosphere.clone(),
                    StationAtmospherePresence::NEITHER,
                ),
            ),
            (
                "new_with_station_presence(BOTH)",
                TrajectorySolver::new_with_station_presence(
                    inputs(),
                    wind(),
                    atmosphere.clone(),
                    StationAtmospherePresence::BOTH,
                ),
            ),
            (
                "new_with_resolved_station_atmosphere",
                TrajectorySolver::new_with_resolved_station_atmosphere(
                    inputs(),
                    wind(),
                    atmosphere.clone(),
                ),
            ),
        ];
        for (name, solver) in solvers {
            assert_percent(name, humidity, solver.inputs.humidity_percent);
        }
    }
}

#[test]
fn the_monte_carlo_sampler_flies_the_inputs_percent() {
    let params = MonteCarloParams {
        num_simulations: 1,
        target_distance: Some(300.0),
        ..MonteCarloParams::default()
    };
    for humidity in HUMIDITIES {
        let sampler = MonteCarloTrialSampler::new(
            load(humidity),
            &WindConditions::default(),
            &params,
            0.0,
            StationAtmospherePresence::BOTH,
        )
        .expect("baseline solve");
        assert_percent(
            "MonteCarloTrialSampler",
            humidity,
            sampler.atmosphere.humidity,
        );
    }
}

#[test]
fn the_lean_monte_carlo_flies_the_inputs_percent_exactly() {
    for humidity in HUMIDITIES {
        let inputs = load(humidity);
        let (result, trace) =
            traced(|| crate::monte_carlo::solve_trajectory_for_monte_carlo(&inputs));
        result.expect("lean Monte Carlo solve");
        assert!(!trace.is_empty(), "the lean Monte Carlo composes air");
        for entry in &trace {
            assert_percent(
                "solve_trajectory_for_monte_carlo",
                humidity,
                entry.humidity_percent,
            );
        }
    }
}

#[test]
fn the_hold_curve_builder_writes_the_percent_into_both_structs() {
    for humidity in HUMIDITIES {
        let (inputs, _wind, atmosphere) = crate::hold_curve::build_trajectory_components(
            823.0,
            0.243,
            0.01134,
            0.00782,
            DragModel::G7,
            0.05,
            15.0,
            1013.25,
            humidity,
            0.0,
            0.0,
            0.0,
            1_000.0,
            10.0,
            None,
            None,
            0.0,
            0.0,
            0.0,
        );
        assert_percent(
            "hold_curve::build_trajectory_components",
            humidity,
            inputs.humidity_percent,
        );
        assert_eq!(atmosphere.humidity.to_bits(), humidity.to_bits());
    }
}

#[test]
fn the_dsf_builder_writes_the_percent_into_both_structs() {
    for humidity in HUMIDITIES {
        let dsf = crate::truing_dsf::DsfSolveInputs {
            muzzle_velocity_fps: 2700.0,
            ballistic_coefficient: 0.243,
            drag_model: DragModel::G7,
            mass_gr: 175.0,
            diameter_in: 0.308,
            sight_height_in: 1.5,
            temperature_f: StationReading::supplied(59.0),
            pressure_inhg: StationReading::supplied(29.92),
            humidity_pct: humidity,
            altitude_ft: 0.0,
            wind_speed_mps: None,
            wind_direction_rad: None,
            shooting_angle_rad: None,
            zero_poi_vertical_m: None,
            zero_poi_horizontal_m: None,
            sight_offset_lateral_m: None,
            bc_reference_standard: None,
            use_bc_segments: false,
            bc_segments: None,
            custom_drag_table: None,
            zero_distance_yd: None,
            bore_height_m: 1.524,
        };
        let flight = crate::truing_dsf::dsf_flight(&dsf, 1_000.0);
        assert_percent(
            "truing_dsf::dsf_flight",
            humidity,
            flight.inputs.humidity_percent,
        );
        assert_eq!(flight.atmosphere.humidity.to_bits(), humidity.to_bits());
    }
}

#[cfg(feature = "ffi")]
#[test]
fn the_ffi_monte_carlo_builder_writes_the_callers_percent() {
    use crate::ffi::{FFIAtmosphericConditions, FFIBallisticInputs};
    let inputs = FFIBallisticInputs {
        muzzle_velocity: 800.0,
        muzzle_angle: 0.0,
        bc_value: 0.5,
        bullet_mass: 0.01,
        bullet_diameter: 0.00762,
        bc_type: 0,
        sight_height: 0.05,
        target_distance: 300.0,
        temperature: 15.0,
        twist_rate: 12.0,
        is_twist_right: 1,
        shooting_angle: 0.0,
        altitude: 0.0,
        latitude: f64::NAN,
        azimuth_angle: 0.0,
        use_rk4: 1,
        use_adaptive_rk45: 0,
        enable_wind_shear: 0,
        enable_trajectory_sampling: 0,
        sample_interval: 10.0,
        enable_pitch_damping: 0,
        enable_precession_nutation: 0,
        enable_spin_drift: 0,
        enable_magnus: 0,
        enable_coriolis: 0,
        shot_azimuth: 0.0,
        cant_angle: 0.0,
        zero_poi_vertical: 0.0,
        zero_poi_horizontal: 0.0,
        sight_offset_lateral: 0.0,
    };
    for humidity in HUMIDITIES {
        let atmosphere = FFIAtmosphericConditions {
            temperature: 20.0,
            pressure: 1000.0,
            humidity,
            altitude: 300.0,
        };
        let (built, presence) = crate::ffi::monte_carlo_base_inputs(
            &inputs,
            Some(&atmosphere),
            StationAtmospherePresence::BOTH,
        );
        assert_percent(
            "ffi::monte_carlo_base_inputs",
            humidity,
            built.humidity_percent,
        );
        assert_eq!(presence, StationAtmospherePresence::BOTH);
    }
    // A null atmosphere carries no reading: the inputs keep the default percent.
    let (built, presence) =
        crate::ffi::monte_carlo_base_inputs(&inputs, None, StationAtmospherePresence::BOTH);
    assert_percent(
        "ffi::monte_carlo_base_inputs(null)",
        DEFAULT_HUMIDITY_PERCENT,
        built.humidity_percent,
    );
    assert_eq!(presence, StationAtmospherePresence::NEITHER);
}

#[test]
fn the_solve_v1_builder_writes_the_wire_fraction_as_one_percent() {
    // The SI wire's relative_humidity is a fraction; the engine structs both carry the one
    // percent RelativeHumidity::try_fraction makes of it.
    for humidity in HUMIDITIES {
        let fraction = RelativeHumidity::try_percent(humidity).unwrap().fraction();
        let request = crate::solve_json::decode_solve_request_v1(&format!(
            r#"{{"schema_version": 1,
                "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782,
                                "drag_model": "G7", "ballistic_coefficient": 0.243}},
                "rifle": {{"muzzle_velocity_mps": 823.0}},
                "shot": {{"max_range_m": 500.0}},
                "atmosphere": {{"relative_humidity": {fraction:?}}},
                "wind": {{"speed_mps": 0.0, "direction_from_rad": 0.0}},
                "solver": {{"method": "rk45"}},
                "effects": {{}},
                "sampling": {{"interval_m": 100.0}}}}"#
        ))
        .expect("valid request");
        let prepared = crate::solve_v1::prepare_request(&request).expect("prepared");
        let percent = RelativeHumidity::try_fraction(fraction).unwrap().percent();
        assert_percent(
            "solve_v1::prepare_request",
            percent,
            prepared.inputs.humidity_percent,
        );
        assert_eq!(prepared.atmosphere.humidity.to_bits(), percent.to_bits());
    }
}

fn trajectory_params(atmosphere: KernelAtmosphere) -> TrajectoryParams {
    TrajectoryParams {
        mass_kg: 0.01134,
        bc: 0.243,
        drag_model: DragModel::G7,
        wind_segments: Vec::new(),
        atmosphere,
        omega_vector: None,
        enable_spin_drift: false,
        enable_magnus: false,
        enable_coriolis: false,
        target_distance_m: 500.0,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: 1_500.0,
        is_twist_right: true,
        shooting_angle: 0.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.031,
        twist_rate: 10.0,
        custom_drag_table: None,
        cd_scale: 1.0,
        custom_drag_reference: false,
        bc_segments: None,
        use_bc_segments: false,
        bc_segments_data: None,
        ground_threshold: -1000.0,
        atmo_sock: None,
    }
}

#[test]
fn the_kernel_derivative_inputs_carry_the_station_percent() {
    for humidity in HUMIDITIES {
        let params = trajectory_params(KernelAtmosphere::Station(station(humidity)));
        let inputs = crate::trajectory_integration::build_inputs(&params, 823.0);
        assert_percent(
            "trajectory_integration::build_inputs",
            humidity,
            inputs.humidity_percent,
        );
    }
    // Uniform air has no station: the default percent.
    let params = trajectory_params(KernelAtmosphere::uniform(1.1, 335.0).unwrap());
    let inputs = crate::trajectory_integration::build_inputs(&params, 823.0);
    assert_percent(
        "trajectory_integration::build_inputs(uniform)",
        DEFAULT_HUMIDITY_PERCENT,
        inputs.humidity_percent,
    );
}

#[test]
fn the_fast_path_copies_the_station_percent_into_its_inputs() {
    for humidity in HUMIDITIES {
        let mut inputs = load(80.0);
        crate::fast_trajectory::station_into_inputs(
            &KernelAtmosphere::Station(station(humidity)),
            &mut inputs,
        );
        assert_percent(
            "fast_trajectory::station_into_inputs",
            humidity,
            inputs.humidity_percent,
        );
    }
    // Uniform air has no station and leaves the inputs as the caller set them.
    let mut inputs = load(80.0);
    crate::fast_trajectory::station_into_inputs(
        &KernelAtmosphere::uniform(1.1, 335.0).unwrap(),
        &mut inputs,
    );
    assert_percent(
        "fast_trajectory::station_into_inputs(uniform)",
        80.0,
        inputs.humidity_percent,
    );
}

#[test]
fn the_wez_attribution_echoes_the_inputs_percent_as_the_wire_fraction() {
    for humidity in HUMIDITIES {
        let resolved =
            crate::wez::wez_resolved_request(&load(humidity), &WindConditions::default(), 2_000.0)
                .expect("a built-in drag model");
        assert_eq!(
            resolved.atmosphere.relative_humidity.to_bits(),
            RelativeHumidity::try_percent(humidity)
                .unwrap()
                .fraction()
                .to_bits(),
            "{humidity} %"
        );
    }
}
