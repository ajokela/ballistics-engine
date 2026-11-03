//! A humidity that is not a number never becomes steam.
//!
//! The engine's raw humidity primitives take a bare `f64` percent. They clamp a number outside
//! 0-100 into it, a total-function guard for unvalidated callers, but before 0.50.0 they also
//! turned a NaN into pure water vapour: `f64::clamp` lets a NaN through, and the `.min(1.0)`
//! cap on the vapour mole fraction then returned 1.0 for it (`NaN.min(1.0) == 1.0`). Air at
//! 15 °C and 1013.25 hPa came out at 0.81 kg/m³ and 387 m/s. The same cap did it above water's
//! critical point, where the saturation pressure itself is NaN, even at 0 %.
//!
//! The one vapour-fraction cap now keeps a NaN a NaN, so every primitive answers NaN for a NaN
//! humidity and for air above the critical point, and every finite answer is bit for bit what
//! it was: the cap differs from `.min(1.0)` for a NaN only.

use super::*;

/// The saturated air these tests compare against: what a NaN used to fly as.
fn saturated_station() -> (f64, f64) {
    calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 100.0)
}

#[test]
fn a_nan_humidity_is_nan_air_through_every_raw_primitive() {
    let nan = f64::NAN;
    let (density, sound) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), nan);
    assert!(
        density.is_nan() && sound.is_nan(),
        "calculate_atmosphere at a NaN humidity is ({density}, {sound}), not NaN; saturated air \
         is {:?}",
        saturated_station()
    );
    // The ICAO-standard branch of the same function, with no overrides.
    let (density, sound) = calculate_atmosphere(1000.0, None, None, nan);
    assert!(density.is_nan() && sound.is_nan(), "({density}, {sound})");

    let cimp = calculate_air_density_cimp(15.0, 1013.25, nan);
    assert!(cimp.is_nan(), "calculate_air_density_cimp: {cimp}");
    let cipm = calculate_air_density_cipm(15.0, 1013.25, nan);
    assert!(cipm.is_nan(), "calculate_air_density_cipm: {cipm}");
    let sound = moist_speed_of_sound(288.15, 101_325.0, nan);
    assert!(sound.is_nan(), "moist_speed_of_sound: {sound}");

    let (density, sound) = get_local_atmosphere_humid(1000.0, 0.0, 15.0, 1013.25, 1.0, nan);
    assert!(
        density.is_nan() && sound.is_nan(),
        "get_local_atmosphere_humid: ({density}, {sound})"
    );
    let (density, sound) = get_local_atmosphere_humid_with_gravity(
        1000.0,
        0.0,
        15.0,
        1013.25,
        1.0,
        nan,
        Gravity::EARTH,
    );
    assert!(density.is_nan() && sound.is_nan(), "({density}, {sound})");

    let altitude = density_matched_altitude_m(15.0, 1013.25, nan);
    assert!(altitude.is_nan(), "density_matched_altitude_m: {altitude}");
}

#[test]
fn a_nan_zone_humidity_flies_nan_air_at_every_height_not_steam() {
    let station = KernelAtmosphere::Station(
        StationAir::new(0.0, 15.0, 1013.25, RelativeHumidity::DEFAULT).expect("ordinary air"),
    );
    // AtmoSock::new does not validate; AtmoSock::validate refuses this zone, but air_at, as a
    // per-step lookup, never calls it.
    let zones = AtmoSock::new(vec![(15.0, 1013.25, f64::NAN, 1_000.0)]);
    assert!(
        zones.validate().is_err(),
        "the zone is refused where it is validated"
    );
    for altitude_m in [0.0, 1_000.0, 5_000.0] {
        let air = air_at(&station, Some(&zones), 10.0, altitude_m);
        assert!(
            air.density_kg_m3.is_nan() && air.speed_of_sound_mps.is_nan(),
            "air_at({altitude_m} m) in a NaN zone is {air:?}, not NaN"
        );
    }
    // A finite out-of-range zone humidity is still clamped, as it always was.
    let wet = AtmoSock::new(vec![(15.0, 1013.25, 150.0, 1_000.0)]);
    let saturated = AtmoSock::new(vec![(15.0, 1013.25, 100.0, 1_000.0)]);
    assert_eq!(
        air_at(&station, Some(&wet), 10.0, 1_000.0),
        air_at(&station, Some(&saturated), 10.0, 1_000.0)
    );
}

#[test]
fn above_waters_critical_point_the_vapour_model_answers_nan_not_steam() {
    // 400 °C is above 647.096 K, where the Wagner-Pruss saturation pressure is undefined.
    // Before 0.50.0 the cap read its NaN as pure steam, 0.293 kg/m³ and 591 m/s, at any
    // humidity, 0 % included.
    for humidity in [0.0, 50.0] {
        let (density, sound) = calculate_atmosphere(0.0, Some(400.0), Some(1013.25), humidity);
        assert!(
            density.is_nan() && sound.is_nan(),
            "{humidity} % at 400 °C: ({density}, {sound})"
        );
    }
    // Every surface that builds a station refuses it, as before.
    assert_eq!(
        StationAir::new(0.0, 400.0, 1013.25, RelativeHumidity::DRY)
            .unwrap_err()
            .reason(),
        AtmosphereErrorReason::ImpossibleAir
    );
}

#[test]
fn the_vapour_fraction_cap_matches_min_one_for_every_number() {
    let samples = [
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        1e-300,
        0.25,
        0.999_999_999_999_999_9,
        1.0,
        1.000_000_000_000_000_2,
        2.0,
        1e300,
        f64::INFINITY,
        f64::NEG_INFINITY,
        -1.0,
    ];
    for x in samples {
        assert_eq!(
            capped_vapor_mole_fraction(x).to_bits(),
            x.min(1.0).to_bits(),
            "the cap must be bit-identical to min(1.0) for {x:e}"
        );
    }
    assert!(capped_vapor_mole_fraction(f64::NAN).is_nan());
}
