//! Pins the deliberate station-atmosphere divergence between the two public solve surfaces.
//!
//! The C ABI (`src/ffi.rs`) resolves `FFIAtmosphericConditions` through the legacy default
//! sentinels: at a nonzero altitude, 15 °C / 1013.25 hPa mean "omitted, use the ICAO standard
//! here". solve-json (`src/solve_v1.rs`) builds through
//! `TrajectorySolver::new_with_resolved_station_atmosphere` and trusts a supplied station
//! reading as given. Neither is wrong — the C struct has no presence channel and the JSON DTO
//! does (MBA-1397) — but they return different numbers for the same nominal inputs, so the
//! behaviour is pinned here rather than left to be rediscovered.
//!
//! If one of these assertions starts failing, the fix is NOT to relax it: either a surface's
//! resolution mode changed (which silently moves every existing caller's numbers and needs to
//! be a deliberate, announced break), or the documentation in `src/ffi.rs` and
//! `AtmosphereV1` is now wrong and must be updated with it.

#![cfg(feature = "ffi")]

use ballistics_engine::atmosphere::resolve_station_conditions;
use ballistics_engine::ffi::{
    ballistics_calculate_zero_angle, FFIAtmosphericConditions, FFIBallisticInputs,
    FFIWindConditions,
};
use ballistics_engine::{
    calculate_zero_angle_with_resolved_conditions, AtmosphericConditions, BallisticInputs,
    DragModel, WindConditions,
};

/// 2000 m: high enough that the ICAO lapse is worth real drop, and well past the `> 1 m`
/// guard that keeps the sentinels from firing at sea level.
const ALTITUDE_M: f64 = 2000.0;
const ZERO_DISTANCE_M: f64 = 300.0;
const SIGHT_HEIGHT_M: f64 = 0.05;
/// Exactly the sentinel values -- this is the whole point of the fixture.
const SENTINEL_TEMP_C: f64 = 15.0;
const SENTINEL_PRESSURE_HPA: f64 = 1013.25;
const HUMIDITY_PERCENT: f64 = 50.0;

const MOA_PER_RAD: f64 = 3437.7467707849;

/// 168 gr .308 at 2700 fps -- the fixture every other zero test in this suite uses.
fn native_inputs() -> BallisticInputs {
    BallisticInputs {
        muzzle_velocity: 823.0,
        bc_value: 0.475,
        bc_type: DragModel::G1,
        bullet_mass: 0.01088,
        bullet_diameter: 0.00782,
        sight_height: SIGHT_HEIGHT_M,
        altitude: ALTITUDE_M,
        ..Default::default()
    }
}

/// The same rifle as `native_inputs`, expressed over the C ABI. `convert_inputs` starts from
/// `BallisticInputs::default()` and overwrites exactly these fields, so the two agree.
fn ffi_inputs() -> FFIBallisticInputs {
    FFIBallisticInputs {
        muzzle_velocity: 823.0,
        muzzle_angle: 0.0,
        bc_value: 0.475,
        bullet_mass: 0.01088,
        bullet_diameter: 0.00782,
        bc_type: 0, // G1
        sight_height: SIGHT_HEIGHT_M,
        target_distance: ZERO_DISTANCE_M,
        temperature: SENTINEL_TEMP_C,
        twist_rate: BallisticInputs::default().twist_rate,
        is_twist_right: 1,
        shooting_angle: 0.0,
        altitude: ALTITUDE_M,
        latitude: f64::NAN,
        azimuth_angle: 0.0,
        use_rk4: 1,
        use_adaptive_rk45: 0,
        enable_wind_shear: 0,
        enable_trajectory_sampling: 0,
        sample_interval: 0.0,
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
    }
}

fn ffi_atmosphere(temperature_c: f64, pressure_hpa: f64) -> FFIAtmosphericConditions {
    FFIAtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity: HUMIDITY_PERCENT,
        altitude: ALTITUDE_M,
    }
}

fn ffi_zero_angle(temperature_c: f64, pressure_hpa: f64) -> f64 {
    let inputs = ffi_inputs();
    let wind = FFIWindConditions {
        speed: 0.0,
        direction: 0.0,
        vertical_speed: 0.0,
    };
    let atmosphere = ffi_atmosphere(temperature_c, pressure_hpa);
    unsafe { ballistics_calculate_zero_angle(&inputs, &wind, &atmosphere, ZERO_DISTANCE_M) }
}

/// The solve-json / bridge behaviour: whatever station reading it is handed is authoritative.
fn authoritative_zero_angle(temperature_c: f64, pressure_hpa: f64) -> f64 {
    calculate_zero_angle_with_resolved_conditions(
        native_inputs(),
        ZERO_DISTANCE_M,
        SIGHT_HEIGHT_M,
        WindConditions::default(),
        AtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity: HUMIDITY_PERCENT,
            altitude: ALTITUDE_M,
        },
    )
    .expect("zero angle solves")
}

/// The C ABI reads the sentinels as "omitted": the answer it gives for 15 °C / 1013.25 hPa at
/// altitude is bit-for-bit the answer for the ICAO conditions there, not for the literal
/// numbers passed in.
#[test]
fn c_abi_reinterprets_sentinel_station_conditions_as_icao_at_altitude() {
    let (icao_temp_c, icao_pressure_hpa) =
        resolve_station_conditions(SENTINEL_TEMP_C, SENTINEL_PRESSURE_HPA, ALTITUDE_M);

    // Guard the fixture itself: if the sentinel bands ever stop firing here, every assertion
    // below would pass vacuously.
    assert!(
        (icao_temp_c - SENTINEL_TEMP_C).abs() > 1.0
            && (icao_pressure_hpa - SENTINEL_PRESSURE_HPA).abs() > 1.0,
        "fixture no longer trips the sentinels: resolved {icao_temp_c} °C / {icao_pressure_hpa} hPa"
    );

    let via_sentinels = ffi_zero_angle(SENTINEL_TEMP_C, SENTINEL_PRESSURE_HPA);
    let via_icao_spelled_out = ffi_zero_angle(icao_temp_c, icao_pressure_hpa);

    assert!(via_sentinels.is_finite(), "C ABI zero angle returned NaN");
    assert_eq!(
        via_sentinels, via_icao_spelled_out,
        "the C ABI is supposed to resolve sentinel station conditions to the ICAO atmosphere \
         at {ALTITUDE_M} m ({icao_temp_c} °C / {icao_pressure_hpa} hPa); it no longer does"
    );
}

/// The two surfaces disagree, on purpose, and by an amount a shooter would see. Pinned as a
/// floor rather than an exact value so ordinary physics work does not churn this test.
#[test]
fn c_abi_and_solve_json_disagree_on_sentinel_station_conditions() {
    let c_abi = ffi_zero_angle(SENTINEL_TEMP_C, SENTINEL_PRESSURE_HPA);
    let solve_json = authoritative_zero_angle(SENTINEL_TEMP_C, SENTINEL_PRESSURE_HPA);

    let delta_moa = (c_abi - solve_json).abs() * MOA_PER_RAD;
    assert!(
        delta_moa > 0.2,
        "the documented C-ABI-vs-solve-json divergence has collapsed to {delta_moa:.4} MOA \
         (C ABI {c_abi} rad, solve-json {solve_json} rad). Either a surface's station-atmosphere \
         resolution changed -- which moves every existing caller's numbers -- or the docs in \
         src/ffi.rs and on AtmosphereV1 need updating."
    );
}

/// Outside the sentinel bands the divergence vanishes: a C caller who states a real station
/// reading gets exactly what solve-json would give it. This is what makes the sentinel an
/// omission channel rather than a second physics model.
#[test]
fn the_two_surfaces_agree_once_the_station_reading_escapes_the_sentinel_bands() {
    // 0.2 °C / 1 hPa off the defaults -- just past the 0.1 °C and 0.5 hPa tolerances.
    let temperature_c = SENTINEL_TEMP_C + 0.2;
    let pressure_hpa = SENTINEL_PRESSURE_HPA - 1.0;

    let c_abi = ffi_zero_angle(temperature_c, pressure_hpa);
    let solve_json = authoritative_zero_angle(temperature_c, pressure_hpa);

    assert_eq!(
        c_abi, solve_json,
        "an explicitly-stated station reading must be authoritative on BOTH surfaces"
    );
}
