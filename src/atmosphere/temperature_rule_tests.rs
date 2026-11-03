//! One air-temperature rule, and where the air of a flight reaches absolute zero.
//!
//! A station temperature is refused on its own only when it is not finite or at or below
//! absolute zero ([`check_station_temperature_c`]); the impossible-air rule is the only other
//! reason. Dry air holds no vapour at any temperature, so it flies above water's critical point
//! until its CIPM-2007 density gives out. And a flight is refused where the lapse from its base
//! temperature takes the air to 0 K ([`absolute_zero_altitudes_m`], [`air_aloft_refusal`]),
//! which these tests hold to the very temperatures [`air_at`] flies.

use super::*;

fn dry_station(altitude_m: f64, temperature_k: f64, pressure_hpa: f64) -> StationAir {
    StationAir::new(
        altitude_m,
        temperature_k - KELVIN_OFFSET,
        pressure_hpa,
        RelativeHumidity::DRY,
    )
    .expect("a dry station above absolute zero flies")
}

/// The temperature, K, the lapse of [`air_at`] carries a base to at `altitude_m`.
fn lapsed_temperature_k(base_temp_k: f64, base_alt_m: f64, altitude_m: f64) -> f64 {
    local_temp_pressure_density(
        altitude_m,
        base_alt_m,
        base_temp_k - KELVIN_OFFSET,
        1013.25,
        1.0,
    )
    .0
}

#[test]
fn the_one_rule_refuses_only_a_temperature_that_is_not_finite_or_at_absolute_zero() {
    for temperature_c in [
        -273.149_999,
        -200.0,
        -150.0,
        -100.000_001,
        100.000_001,
        250.0,
        1e6,
    ] {
        assert!(
            check_station_temperature_c(temperature_c).is_ok(),
            "{temperature_c} °C is a temperature (extrapolated outside -100..100 °C)"
        );
    }
    for (temperature_c, reason) in [
        (-273.15, AtmosphereErrorReason::AtOrBelowAbsoluteZero),
        (-300.0, AtmosphereErrorReason::AtOrBelowAbsoluteZero),
        (f64::NAN, AtmosphereErrorReason::NotFinite),
        (f64::INFINITY, AtmosphereErrorReason::NotFinite),
        (f64::NEG_INFINITY, AtmosphereErrorReason::NotFinite),
    ] {
        let error = check_station_temperature_c(temperature_c).unwrap_err();
        assert_eq!(error.reason(), reason, "{temperature_c}");
        assert_eq!(error.field(), "temperature_c");
    }
    // StationAir::new applies exactly this rule to the temperature (then impossible air).
    for temperature_c in [-273.15, -200.0, 250.0, f64::NAN] {
        let station = StationAir::new(0.0, temperature_c, 1013.25, RelativeHumidity::DRY);
        assert_eq!(
            station.is_ok(),
            check_station_temperature_c(temperature_c).is_ok(),
            "{temperature_c} °C"
        );
    }
}

#[test]
fn a_flag_refuses_in_its_own_unit_exactly_what_the_rule_refuses() {
    for (value, fahrenheit) in [
        (-459.67, true),
        (-459.66, true),
        (-273.15, false),
        (-273.149, false),
        (482.0, true),
        (250.0, false),
        (f64::NAN, true),
        (f64::INFINITY, false),
    ] {
        let celsius = if fahrenheit {
            (value - 32.0) * 5.0 / 9.0
        } else {
            value
        };
        assert_eq!(
            check_temperature_flag("--temperature", value, fahrenheit).is_ok(),
            check_station_temperature_c(celsius).is_ok(),
            "{value} (fahrenheit: {fahrenheit})"
        );
    }
    assert_eq!(
        check_temperature_flag("--zero-temperature", -273.15, false).unwrap_err(),
        "--zero-temperature -273.15 C is at or below absolute zero (-273.15 C)"
    );
    assert_eq!(
        check_temperature_flag("--temperature", f64::NAN, true).unwrap_err(),
        "--temperature NaN F is not a finite number"
    );
}

/// Dry air above water's critical point flies until its CIPM-2007 compressibility gives no
/// positive density: about 90,000 K at 1013.25 hPa (Z reaches 0 at 90,181 K there).
#[test]
fn dry_air_flies_above_the_critical_point_until_its_density_model_gives_out() {
    for temperature_c in [373.95, 400.0, 1_000.0, 10_000.0, 89_000.0] {
        let station = StationAir::new(0.0, temperature_c, 1013.25, RelativeHumidity::DRY)
            .unwrap_or_else(|error| panic!("dry air at {temperature_c} °C: {error}"));
        assert!(station.density_kg_m3() > 0.0 && station.density_kg_m3().is_finite());
        assert_eq!(
            station.speed_of_sound_mps().to_bits(),
            (GAMMA * R_AIR * (temperature_c + KELVIN_OFFSET))
                .sqrt()
                .to_bits()
        );
    }
    for temperature_c in [91_000.0, 1e6, 1e300] {
        let error = StationAir::new(0.0, temperature_c, 1013.25, RelativeHumidity::DRY)
            .expect_err("no positive density");
        assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
        assert_eq!(error.field(), "temperature_c");
    }
    // Humid air above the critical point cannot exist at any humidity above 0.
    let error = StationAir::new(
        0.0,
        374.0,
        1013.25,
        RelativeHumidity::try_percent(0.001).unwrap(),
    )
    .expect_err("humid air above the critical point");
    assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
    assert_eq!(error.field(), "temperature_c");
    // Zones are held to the same rule.
    assert!(AtmoSock::new(vec![(400.0, 1013.25, 0.0, 100.0)])
        .validate()
        .is_ok());
    assert!(AtmoSock::new(vec![(400.0, 1013.25, 1.0, 100.0)])
        .validate()
        .is_err());
}

/// Each edge of the band is where the lapse [`air_at`] flies reaches 0 K: the air a metre inside
/// is above absolute zero and has a finite density, the air a metre outside has none.
#[test]
fn the_absolute_zero_altitudes_are_where_the_flown_lapse_reaches_zero() {
    let mut ceilings = 0;
    let mut floors = 0;
    for (base_temp_k, base_alt_m) in [
        (10.0, 0.0),
        (50.0, 0.0),
        (288.15, 0.0),
        (330.0, -400.0),
        (5.0, 1_500.0),
        (60.0, 10_500.0),
        (4.0, 25_000.0),
        (12.0, 40_000.0),
        (150.0, 60_000.0),
    ] {
        let (floor, ceiling) = absolute_zero_altitudes_m(base_temp_k, base_alt_m);
        assert!(
            floor < base_alt_m && base_alt_m < ceiling,
            "{base_temp_k} K at {base_alt_m} m"
        );
        for (edge, inside, outside) in [
            (ceiling, ceiling - 1.0, ceiling + 1.0),
            (floor, floor + 1.0, floor - 1.0),
        ] {
            if !edge.is_finite() {
                continue;
            }
            let t_inside = lapsed_temperature_k(base_temp_k, base_alt_m, inside);
            let t_outside = lapsed_temperature_k(base_temp_k, base_alt_m, outside);
            assert!(
                t_inside > 0.0,
                "{base_temp_k} K at {base_alt_m} m: {t_inside} K at {inside} m"
            );
            assert!(
                t_outside.is_nan() || t_outside <= 0.0,
                "{base_temp_k} K at {base_alt_m} m: {t_outside} K at {outside} m"
            );
            let station = KernelAtmosphere::Station(dry_station(base_alt_m, base_temp_k, 1013.25));
            let air = air_at(&station, None, 0.0, inside);
            assert!(
                air.density_kg_m3.is_finite() && air.speed_of_sound_mps.is_finite(),
                "{air:?}"
            );
            let air = air_at(&station, None, 0.0, outside);
            assert!(
                !air.speed_of_sound_mps.is_finite() || !air.density_kg_m3.is_finite(),
                "{air:?}"
            );
        }
        ceilings += usize::from(ceiling.is_finite());
        floors += usize::from(floor.is_finite());
    }
    assert!(
        ceilings >= 8 && floors >= 2,
        "{ceilings} ceilings, {floors} floors"
    );
    // The figures the docs quote: a 10 K station at sea level reaches 0 K 1.54 km up, and 15 °C
    // at sea level 183 km up (the last ICAO layer's lapse continues above 84.852 km).
    let (_, ten_kelvin) = absolute_zero_altitudes_m(10.0, 0.0);
    assert!((ten_kelvin - 1_538.83).abs() < 0.01, "{ten_kelvin}");
    let (standard_floor, standard) = absolute_zero_altitudes_m(288.15, 0.0);
    assert!((183_000.0..184_000.0).contains(&standard), "{standard}");
    assert_eq!(standard_floor, f64::NEG_INFINITY);
}

#[test]
fn a_failed_step_is_blamed_on_the_air_only_when_its_reach_crosses_absolute_zero() {
    let station = KernelAtmosphere::Station(dry_station(0.0, 10.0, 10.0));
    let (_, ceiling) = absolute_zero_altitudes_m(10.0, 0.0);
    let refusal = air_aloft_refusal(&station, None, 500.0, ceiling - 1.0, 2.0)
        .expect("the step's reach crosses the ceiling");
    assert_eq!(
        refusal.reason(),
        AtmosphereErrorReason::AirAloftBelowAbsoluteZero
    );
    assert_eq!(refusal.field(), "temperature_c");
    assert_eq!(refusal.value().to_bits(), (10.0 - KELVIN_OFFSET).to_bits());
    assert_eq!(refusal.altitude_m(), Some(ceiling));
    assert_eq!(
        refusal.to_string(),
        format!(
            "temperature_c = {}: the air aloft would fall below absolute zero: from this \
             temperature the standard lapse reaches 0 K at 1539 m above mean sea level, which \
             this flight reaches (check the temperature and its unit)",
            10.0 - KELVIN_OFFSET
        )
    );
    // A step that stays below the ceiling failed for another reason.
    assert!(air_aloft_refusal(&station, None, 500.0, ceiling - 10.0, 2.0).is_none());
    // Uniform air does not lapse, and a non-finite start says nothing.
    let uniform = KernelAtmosphere::uniform(1.0, 330.0).unwrap();
    assert!(air_aloft_refusal(&uniform, None, 0.0, 1e9, 1.0).is_none());
    assert!(air_aloft_refusal(&station, None, f64::NAN, ceiling, 1.0).is_none());
    // A zone's own temperature governs the air it covers: a 300 K station whose first 1 km is a
    // 10 K zone is refused there, and the air beyond, the station's, is not.
    let warm = KernelAtmosphere::Station(dry_station(0.0, 300.0, 1013.25));
    let zones = AtmoSock::new(vec![(10.0 - KELVIN_OFFSET, 10.0, 0.0, 1_000.0)]);
    zones.validate().unwrap();
    let in_zone = air_aloft_refusal(&warm, Some(&zones), 500.0, ceiling, 1.0).expect("the zone");
    assert_eq!(in_zone.value().to_bits(), (10.0 - KELVIN_OFFSET).to_bits());
    // Past the last zone the last zone still governs (a sock has no end), so test the station
    // alone beyond a sock that covers nothing.
    assert!(air_aloft_refusal(&warm, None, 5_000.0, ceiling, 1.0).is_none());
}
