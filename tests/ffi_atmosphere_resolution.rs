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

// ---------------------------------------------------------------------------------------
// 0.49.0: the presence-aware exports
// ---------------------------------------------------------------------------------------
//
// `FFIAtmosphericConditions` still cannot say which of its doubles were supplied, so the
// struct exports above keep the legacy rule. The `*_with_station_presence` exports take the
// presence beside it as bits. Every legacy number below was captured at ab4eb0c; every
// new-export number was captured from ab4eb0c with the two sentinels compiled out, where a
// supplied value flies as given.

mod station_presence_exports {
    use ballistics_engine::ffi::{
        ballistics_calculate_trajectory, ballistics_calculate_trajectory_with_station_presence,
        ballistics_calculate_zero_angle, ballistics_calculate_zero_angle_with_station_presence,
        ballistics_density_altitude_altitude_m, ballistics_density_altitude_pressure_hpa,
        ballistics_density_altitude_temperature_c, ballistics_free_monte_carlo_results,
        ballistics_free_trajectory_result, ballistics_monte_carlo,
        ballistics_monte_carlo_with_station_presence, ballistics_reduce_qnh_pressure,
        FFIAtmosphericConditions, FFIBallisticInputs, FFIMonteCarloParams,
        BALLISTICS_STATION_PRESSURE_SUPPLIED, BALLISTICS_STATION_TEMPERATURE_SUPPLIED,
    };
    use std::ptr;

    const T: u32 = BALLISTICS_STATION_TEMPERATURE_SUPPLIED;
    const P: u32 = BALLISTICS_STATION_PRESSURE_SUPPLIED;

    /// The C ABI rifle the pinned rows were captured with, with `temperature` echoed into the
    /// inputs.
    fn inputs(altitude: f64, temperature: f64) -> FFIBallisticInputs {
        FFIBallisticInputs {
            muzzle_velocity: 823.0,
            muzzle_angle: 0.0,
            bc_value: 0.475,
            bullet_mass: 0.0109,
            bullet_diameter: 0.00782,
            bc_type: 0,
            sight_height: 0.0381,
            target_distance: 600.0,
            temperature,
            twist_rate: 10.0,
            is_twist_right: 1,
            shooting_angle: 0.0,
            altitude,
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
        }
    }

    fn air(t: f64, p: f64, alt: f64) -> FFIAtmosphericConditions {
        FFIAtmosphericConditions {
            temperature: t,
            pressure: p,
            humidity: 50.0,
            altitude: alt,
        }
    }

    fn legacy_zero(t: f64, p: f64, alt: f64) -> f64 {
        unsafe {
            ballistics_calculate_zero_angle(&inputs(alt, t), ptr::null(), &air(t, p, alt), 600.0)
        }
    }

    fn presence_zero(t: f64, p: f64, alt: f64, flags: u32) -> f64 {
        unsafe {
            ballistics_calculate_zero_angle_with_station_presence(
                &inputs(alt, t),
                ptr::null(),
                &air(t, p, alt),
                600.0,
                ptr::null(),
                ptr::null(),
                0,
                1.0,
                flags,
            )
        }
    }

    /// One pinned row: (label, (t, p, alt), flags, legacy export, new export).
    type PinnedRow = (&'static str, (f64, f64, f64), u32, f64, f64);

    /// The pinned zero angles (rad, 600 m), row for row.
    #[test]
    fn zero_export_matches_the_pinned_rows() {
        let qnh = ballistics_reduce_qnh_pressure(1031.0, 150.0);
        let da_t = ballistics_density_altitude_temperature_c(1000.0, 15.0);
        let da_p = ballistics_density_altitude_pressure_hpa(1000.0, 15.0);
        let da_a = ballistics_density_altitude_altitude_m(1000.0, 15.0);
        let rows: [PinnedRow; 7] = [
            (
                "abs 1013.0 at 150 m",
                (20.0, 1013.0, 150.0),
                P,
                0.0062198638916015625,
                0.006264114379882812,
            ),
            (
                "QNH helper 1031 at 150 m",
                (20.0, qnh, 150.0),
                P,
                0.0062198638916015625,
                0.00626373291015625,
            ),
            (
                "15 C at 1000 m",
                (15.0, 900.0, 1000.0),
                T,
                0.006062698364257813,
                0.0060211181640625,
            ),
            (
                "DA helper 1000 m, 15 C",
                (da_t, da_p, da_a),
                T,
                0.006103515625,
                0.006068801879882813,
            ),
            (
                "omitted at 150 m",
                (15.0, 1013.25, 150.0),
                0,
                0.00626373291015625,
                0.00626373291015625,
            ),
            (
                "omitted at 0.5 m",
                (15.0, 1013.25, 0.5),
                0,
                0.006302261352539062,
                0.006302261352539062,
            ),
            (
                "abs 1012.0 at 150 m",
                (20.0, 1012.0, 150.0),
                P,
                0.006261444091796875,
                0.006261444091796875,
            ),
        ];
        for (label, (t, p, alt), flags, legacy, new) in rows {
            assert_eq!(
                legacy_zero(t, p, alt).to_bits(),
                legacy.to_bits(),
                "{label}: legacy export"
            );
            assert_eq!(
                presence_zero(t, p, alt, flags).to_bits(),
                new.to_bits(),
                "{label}: new export, flags {flags}"
            );
            // Flags 0 IS the legacy export, whatever the inputs: compared with the pinned
            // legacy number, not just with the legacy export.
            assert_eq!(
                presence_zero(t, p, alt, 0).to_bits(),
                legacy.to_bits(),
                "{label}: new export, flags 0"
            );
        }
        // Supplying the in-band pressure moves the answer off the legacy one.
        assert_ne!(
            presence_zero(20.0, 1013.0, 150.0, P).to_bits(),
            legacy_zero(20.0, 1013.0, 150.0).to_bits()
        );
    }

    #[test]
    fn unknown_bits_and_a_scale_without_a_deck_are_refused() {
        for flags in [4, 8, 1 << 31, T | P | 4] {
            assert!(
                presence_zero(20.0, 1013.0, 150.0, flags).is_nan(),
                "flags {flags}"
            );
        }
        assert!(presence_zero(20.0, 1013.0, 150.0, T | P).is_finite());
        let scaled = unsafe {
            ballistics_calculate_zero_angle_with_station_presence(
                &inputs(150.0, 20.0),
                ptr::null(),
                &air(20.0, 1013.0, 150.0),
                600.0,
                ptr::null(),
                ptr::null(),
                0,
                1.05,
                P,
            )
        };
        assert!(scaled.is_nan(), "cd_scale without a deck must be refused");
        // A deck is validated as by the _scaled export.
        let mach = [0.5, 1.0, 2.0];
        let cd = [0.2, 0.4, 0.3];
        let with_deck = unsafe {
            ballistics_calculate_zero_angle_with_station_presence(
                &inputs(150.0, 20.0),
                ptr::null(),
                &air(20.0, 1013.0, 150.0),
                600.0,
                mach.as_ptr(),
                cd.as_ptr(),
                3,
                1.0,
                P,
            )
        };
        assert!(with_deck.is_finite());
        let bad_deck = unsafe {
            ballistics_calculate_zero_angle_with_station_presence(
                &inputs(150.0, 20.0),
                ptr::null(),
                &air(20.0, 1013.0, 150.0),
                600.0,
                mach.as_ptr(),
                cd.as_ptr(),
                1,
                1.0,
                P,
            )
        };
        assert!(bad_deck.is_nan());
    }

    /// The capture fixture: a 0.002 rad shot sampled every 50 m.
    fn trajectory_inputs() -> FFIBallisticInputs {
        FFIBallisticInputs {
            muzzle_angle: 0.002,
            temperature: 20.0,
            altitude: 150.0,
            enable_trajectory_sampling: 1,
            sample_interval: 50.0,
            ..inputs(150.0, 20.0)
        }
    }

    /// (impact velocity, max height, sample 6 velocity, sample 6 Mach)
    fn trajectory_headline(flags: Option<u32>) -> Option<[f64; 4]> {
        let atmosphere = air(20.0, 1013.0, 150.0);
        unsafe {
            let result = match flags {
                None => ballistics_calculate_trajectory(
                    &trajectory_inputs(),
                    ptr::null(),
                    &atmosphere,
                    600.0,
                    1.0,
                ),
                Some(flags) => ballistics_calculate_trajectory_with_station_presence(
                    &trajectory_inputs(),
                    ptr::null(),
                    &atmosphere,
                    600.0,
                    1.0,
                    ptr::null(),
                    ptr::null(),
                    0,
                    1.0,
                    flags,
                ),
            };
            if result.is_null() {
                return None;
            }
            let r = &*result;
            let samples =
                std::slice::from_raw_parts(r.sampled_points, r.sampled_point_count as usize);
            let headline = [
                r.impact_velocity,
                r.max_height,
                samples[6].velocity_mps,
                samples[6].mach,
            ];
            ballistics_free_trajectory_result(result);
            Some(headline)
        }
    }

    #[test]
    fn trajectory_export_is_the_legacy_one_with_flags_0_and_flies_a_supplied_pressure() {
        // captured at ab4eb0c (legacy export, in-band 1013.0 hPa at 150 m)
        let pinned = [
            0x407f0676b1cf0f7a,
            0x3fc08c94333f41cc,
            0x40843afdb0c2e019,
            0x3ffe1f1ff823bb7a,
        ];
        let bits = |h: [f64; 4]| h.map(f64::to_bits);
        assert_eq!(bits(trajectory_headline(None).unwrap()), pinned);
        assert_eq!(bits(trajectory_headline(Some(0)).unwrap()), pinned);
        let supplied = trajectory_headline(Some(P)).unwrap();
        assert_ne!(bits(supplied), pinned);
        // 1013 hPa is denser than the ICAO 995 hPa: slower at impact.
        assert!(supplied[0] < f64::from_bits(pinned[0]));
        assert!(trajectory_headline(Some(4)).is_none());
    }

    /// `ballistics_monte_carlo` is unseeded, so the parameters carry no dispersion at all:
    /// every trial is then the baseline, and counts and means are exact.
    fn zero_sigma_monte_carlo(flags: Option<u32>) -> Option<(i32, f64, f64)> {
        let params = FFIMonteCarloParams {
            num_simulations: 8,
            velocity_std_dev: 0.0,
            angle_std_dev: 0.0,
            bc_std_dev: 0.0,
            wind_speed_std_dev: 0.0,
            target_distance: 500.0,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0,
        };
        let atmosphere = air(20.0, 1013.0, 150.0);
        unsafe {
            let results = match flags {
                None => ballistics_monte_carlo(&trajectory_inputs(), &atmosphere, &params),
                Some(flags) => ballistics_monte_carlo_with_station_presence(
                    &trajectory_inputs(),
                    &atmosphere,
                    &params,
                    0.0,
                    flags,
                ),
            };
            if results.is_null() {
                return None;
            }
            let r = &*results;
            let out = (r.num_results, r.mean_range, r.mean_impact_velocity);
            ballistics_free_monte_carlo_results(results);
            Some(out)
        }
    }

    #[test]
    fn monte_carlo_export_is_the_legacy_one_with_flags_0_and_flies_a_supplied_pressure() {
        // captured at ab4eb0c
        let pinned = (8, 0x407fe680b099e436_u64, 0x4080d4ad5eca7c79_u64);
        let bits = |(n, range, velocity): (i32, f64, f64)| (n, range.to_bits(), velocity.to_bits());
        assert_eq!(bits(zero_sigma_monte_carlo(None).unwrap()), pinned);
        assert_eq!(bits(zero_sigma_monte_carlo(Some(0)).unwrap()), pinned);
        let supplied = zero_sigma_monte_carlo(Some(P)).unwrap();
        assert_eq!(supplied.0, 8);
        assert_ne!(supplied.2.to_bits(), pinned.2);
        assert!(zero_sigma_monte_carlo(Some(4)).is_none());
    }
}
