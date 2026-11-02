//! Unit tests for the single humidity source of truth in `atmosphere.rs`: the validated
//! [`RelativeHumidity`], [`StationAir`], [`KernelAtmosphere`] and [`air_at`].
//!
//! Two kinds of test live here.
//! - Validation and conversion: every constructor's bounds, NaN and the infinities, exact 0 and
//!   100, the fraction round trips, and impossible air refused by exactly the predicate the
//!   `atmosphere.density_altitude` service applies.
//! - Bit-identity of the refactor: the merged `vapor_mole_fraction`, `StationAir::new` against
//!   `calculate_atmosphere`, and `air_at` against both the legacy `get_local_atmosphere_humid`
//!   wrapper and a verbatim copy of 0.49.0's `TrajectorySolver::calculate_acceleration` air
//!   block. Those comparisons are `f64::to_bits`, never a tolerance.

use super::*;

// =================================================================================================
// Grids
// =================================================================================================

/// Station temperatures (C) spanning the plan's -40..45 C grid plus the hot and cold edges.
const TEMPS_C: [f64; 12] = [
    -40.0, -12.3, -10.0, 0.0, 7.9, 15.0, 21.7, 30.0, 33.3, 35.0, 45.0, 60.0,
];
/// Station pressures (hPa) spanning 500..1100 hPa.
const PRESSURES_HPA: [f64; 6] = [500.0, 700.0, 850.0, 1000.0, 1013.25, 1100.0];
/// Relative humidities (percent), both exact ends included.
const HUMIDITIES: [f64; 7] = [0.0, 12.5, 37.5, 50.0, 62.5, 90.0, 100.0];
/// Query altitudes relative to the station (m): at, just above, well above, and below it.
const ALT_OFFSETS_M: [f64; 6] = [0.0, 0.5, 50.0, 1_000.0, 12_000.0, -200.0];

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn local_bits(air: LocalAir) -> [u64; 4] {
    [
        air.density_kg_m3.to_bits(),
        air.speed_of_sound_mps.to_bits(),
        air.temperature_k.to_bits(),
        air.pressure_pa.to_bits(),
    ]
}

fn rh(percent: f64) -> RelativeHumidity {
    RelativeHumidity::try_percent(percent).expect("valid test humidity")
}

/// The difference in units in the last place between two finite, same-sign doubles.
fn ulps_apart(a: f64, b: f64) -> u64 {
    a.to_bits().abs_diff(b.to_bits())
}

// =================================================================================================
// RelativeHumidity
// =================================================================================================

#[test]
fn relative_humidity_constants_are_the_documented_default_and_dry_air() {
    assert_eq!(DEFAULT_HUMIDITY_PERCENT.to_bits(), 50.0_f64.to_bits());
    assert_eq!(
        RelativeHumidity::DEFAULT.percent().to_bits(),
        DEFAULT_HUMIDITY_PERCENT.to_bits()
    );
    assert_eq!(
        RelativeHumidity::DEFAULT.fraction().to_bits(),
        0.5_f64.to_bits()
    );
    assert_eq!(RelativeHumidity::DRY.percent().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        RelativeHumidity::DRY.fraction().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(rh(50.0), RelativeHumidity::DEFAULT);
    assert_eq!(rh(0.0), RelativeHumidity::DRY);
}

#[test]
fn try_percent_accepts_the_closed_range_and_stores_the_value_verbatim() {
    for percent in [0.0, 1e-300, 0.0001, 12.5, 37.5, 50.0, 99.9999, 100.0] {
        let humidity = RelativeHumidity::try_percent(percent).expect("in range");
        assert_eq!(humidity.percent().to_bits(), percent.to_bits(), "{percent}");
        assert_eq!(
            humidity.fraction().to_bits(),
            (percent / 100.0).to_bits(),
            "{percent}"
        );
    }
    // Exact ends, by bits.
    assert_eq!(rh(0.0).percent().to_bits(), 0.0_f64.to_bits());
    assert_eq!(rh(100.0).percent().to_bits(), 100.0_f64.to_bits());
    assert_eq!(rh(100.0).fraction().to_bits(), 1.0_f64.to_bits());
}

#[test]
fn negative_zero_percent_is_dry_air_with_positive_zero_bits() {
    let humidity = RelativeHumidity::try_percent(-0.0).expect("-0.0 is zero");
    assert_eq!(humidity.percent().to_bits(), 0.0_f64.to_bits());
    assert_eq!(humidity, RelativeHumidity::DRY);
    let humidity = RelativeHumidity::try_fraction(-0.0).expect("-0.0 is zero");
    assert_eq!(humidity.percent().to_bits(), 0.0_f64.to_bits());
}

#[test]
fn try_percent_rejects_out_of_range_and_non_finite_values_without_clamping() {
    // S5b: a reading within 5 points of the range is saturated or dry air (the band tests
    // above); only a value beyond the band is refused, and it is refused unchanged.
    let out_of_range = [-5.000_000_1, 105.000_000_1, 150.0, 1e300];
    for percent in out_of_range {
        let error = RelativeHumidity::try_percent(percent).expect_err("out of range");
        assert_eq!(error.field(), "humidity_percent", "{percent}");
        assert_eq!(error.value().to_bits(), percent.to_bits(), "{percent}");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::HumidityPercentOutOfRange,
            "{percent}"
        );
    }
    for percent in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = RelativeHumidity::try_percent(percent).expect_err("non-finite");
        assert_eq!(error.field(), "humidity_percent");
        assert_eq!(error.reason(), AtmosphereErrorReason::NotFinite);
        assert!(
            error.value().is_nan() || error.value() == percent,
            "the error echoes the rejected value"
        );
    }
}

#[test]
fn try_fraction_accepts_the_closed_unit_range_and_scales_by_exactly_100() {
    // The scaling is `fraction * 100.0`, the very conversion the SI wire applies today
    // (solve_v1's relative_humidity * 100.0 and the density_altitude service's), so a later
    // migration of those sites through try_fraction is bit-identical.
    for k in 0..=1000 {
        let fraction = f64::from(k) / 1000.0;
        let humidity = RelativeHumidity::try_fraction(fraction).expect("in range");
        assert_eq!(
            humidity.percent().to_bits(),
            (fraction * 100.0).to_bits(),
            "{fraction}"
        );
    }
    assert_eq!(
        RelativeHumidity::try_fraction(0.5).unwrap(),
        RelativeHumidity::try_percent(50.0).unwrap()
    );
    assert_eq!(
        RelativeHumidity::try_fraction(1.0)
            .unwrap()
            .percent()
            .to_bits(),
        100.0_f64.to_bits()
    );
    assert_eq!(
        RelativeHumidity::try_fraction(0.0)
            .unwrap()
            .percent()
            .to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        RelativeHumidity::try_fraction(0.375).unwrap(),
        RelativeHumidity::try_percent(37.5).unwrap()
    );
}

#[test]
fn try_fraction_rejects_out_of_range_and_non_finite_values_without_clamping() {
    for fraction in [1.5, -0.1, -0.050_000_1, 1.050_000_1, 50.0, 100.0] {
        let error = RelativeHumidity::try_fraction(fraction).expect_err("out of range");
        assert_eq!(error.field(), "relative_humidity", "{fraction}");
        assert_eq!(error.value().to_bits(), fraction.to_bits(), "{fraction}");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::HumidityFractionOutOfRange,
            "{fraction}"
        );
    }
    for fraction in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = RelativeHumidity::try_fraction(fraction).expect_err("non-finite");
        assert_eq!(error.field(), "relative_humidity");
        assert_eq!(error.reason(), AtmosphereErrorReason::NotFinite);
    }
}

#[test]
fn fraction_round_trips_are_exact_at_the_ends_and_within_one_ulp_elsewhere() {
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let back = RelativeHumidity::try_fraction(fraction).unwrap().fraction();
        assert_eq!(back.to_bits(), fraction.to_bits(), "{fraction}");
    }
    for k in 0..=10_000 {
        let fraction = f64::from(k) / 10_000.0;
        let back = RelativeHumidity::try_fraction(fraction).unwrap().fraction();
        assert!(
            ulps_apart(back, fraction) <= 1,
            "fraction {fraction} came back as {back}"
        );
        let percent = f64::from(k) / 100.0;
        let back = RelativeHumidity::try_fraction(rh(percent).fraction())
            .unwrap()
            .percent();
        assert!(
            ulps_apart(back, percent) <= 1,
            "percent {percent} came back as {back}"
        );
    }
}

// =================================================================================================
// The sensor-tolerance band (owner decision 2026-11-01, step S5b)
// =================================================================================================

/// Percents a humidity SENSOR can report around saturated air: within 5 points above 100,
/// they are saturated air, stored as exactly 100.
const SATURATED_READINGS: [f64; 6] = [
    100.000_000_000_001,
    100.000_000_1,
    100.0001,
    101.7,
    104.999_999_9,
    105.0,
];
/// ...and around dry air: within 5 points below 0 (and a negative zero), dry air, stored as
/// exactly +0.0.
const DRY_READINGS: [f64; 7] = [-0.0, -1e-300, -0.0001, -2.0, -3.0, -4.999_999_9, -5.0];

#[test]
fn the_sensor_band_is_five_points_either_side_in_both_units() {
    assert_eq!(
        HUMIDITY_SENSOR_TOLERANCE_PERCENT.to_bits(),
        5.0_f64.to_bits()
    );
    // The band edges are the decimal literals themselves, so 105.0, -5.0, 1.05 and -0.05 are
    // inside and the next double out is not.
    assert_eq!(
        [HUMIDITY_PERCENT_BAND.0, HUMIDITY_PERCENT_BAND.1].map(f64::to_bits),
        [-5.0_f64, 105.0].map(f64::to_bits)
    );
    assert_eq!(
        [HUMIDITY_FRACTION_BAND.0, HUMIDITY_FRACTION_BAND.1].map(f64::to_bits),
        [-0.05_f64, 1.05].map(f64::to_bits)
    );
    for edge in [-5.0_f64, 105.0] {
        assert!(RelativeHumidity::try_percent(edge).is_ok(), "{edge}");
        let beyond = f64::from_bits(if edge > 0.0 {
            edge.to_bits() + 1
        } else {
            (-edge).to_bits().wrapping_add(1) | (1 << 63)
        });
        assert!(RelativeHumidity::try_percent(beyond).is_err(), "{beyond}");
    }
    for edge in [-0.05_f64, 1.05] {
        assert!(RelativeHumidity::try_fraction(edge).is_ok(), "{edge}");
        let beyond = f64::from_bits(if edge > 0.0 {
            edge.to_bits() + 1
        } else {
            (-edge).to_bits().wrapping_add(1) | (1 << 63)
        });
        assert!(RelativeHumidity::try_fraction(beyond).is_err(), "{beyond}");
    }
}

/// S5b: the two other places a humidity's range is needed take it from the one rule rather than
/// writing a second one. The MCP schema publishes `RelativeHumidity::FRACTION_BAND`, which is
/// the band `try_fraction` enforces; and the perturbation layer's physical domain is
/// `is_physical_fraction`, the fractions the band leaves where it found them, which is exactly
/// `0.0..=1.0` (the domain that layer had before the band, so its derivatives are unchanged).
#[test]
fn the_published_fraction_band_and_the_physical_domain_come_from_the_one_rule() {
    assert_eq!(
        [
            RelativeHumidity::FRACTION_BAND.0,
            RelativeHumidity::FRACTION_BAND.1
        ]
        .map(f64::to_bits),
        [HUMIDITY_FRACTION_BAND.0, HUMIDITY_FRACTION_BAND.1].map(f64::to_bits)
    );
    let below_one = f64::from_bits(1.0_f64.to_bits() - 1);
    let above_one = f64::from_bits(1.0_f64.to_bits() + 1);
    for fraction in [0.0, -0.0, 1e-300, 0.25, 0.5, below_one, 1.0] {
        assert!(
            RelativeHumidity::is_physical_fraction(fraction),
            "{fraction}"
        );
    }
    for fraction in [
        above_one,
        1.001,
        1.03,
        1.05,
        -1e-300,
        -0.001,
        -0.03,
        -0.05,
        1.06,
        1.5,
        50.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(
            !RelativeHumidity::is_physical_fraction(fraction),
            "{fraction}"
        );
    }
    // Exactly the closed unit interval, on a sweep across both band edges.
    for step in -1_200..=1_200 {
        let fraction = f64::from(step) * 1e-3;
        assert_eq!(
            RelativeHumidity::is_physical_fraction(fraction),
            (0.0..=1.0).contains(&fraction),
            "{fraction}"
        );
    }
}

#[test]
fn try_percent_reads_a_sensor_reading_within_five_points_as_saturated_or_dry_air() {
    for percent in SATURATED_READINGS {
        let humidity = RelativeHumidity::try_percent(percent)
            .unwrap_or_else(|e| panic!("{percent} is a saturated-air reading: {e}"));
        assert_eq!(
            humidity.percent().to_bits(),
            100.0_f64.to_bits(),
            "{percent}"
        );
        assert_eq!(
            humidity.fraction().to_bits(),
            1.0_f64.to_bits(),
            "{percent}"
        );
        assert_eq!(humidity, rh(100.0), "{percent}");
    }
    for percent in DRY_READINGS {
        let humidity = RelativeHumidity::try_percent(percent)
            .unwrap_or_else(|e| panic!("{percent} is a dry-air reading: {e}"));
        assert_eq!(humidity.percent().to_bits(), 0.0_f64.to_bits(), "{percent}");
        assert_eq!(humidity, RelativeHumidity::DRY, "{percent}");
    }
    // The exact ends and the inside of the range are untouched, bit for bit.
    for percent in [0.0, 1e-300, 0.0001, 7.0, 55.0, 99.9999, 100.0] {
        assert_eq!(
            rh(percent).percent().to_bits(),
            percent.to_bits(),
            "{percent}"
        );
    }
}

#[test]
fn try_percent_refuses_a_reading_beyond_the_sensor_band_unchanged() {
    for percent in [
        105.000_000_1,
        105.000_000_000_001,
        106.0,
        150.0,
        1e300,
        -5.000_000_1,
        -5.000_000_000_001,
        -20.0,
        -1e300,
    ] {
        let error = RelativeHumidity::try_percent(percent).expect_err("beyond the band");
        assert_eq!(error.field(), "humidity_percent", "{percent}");
        assert_eq!(error.value().to_bits(), percent.to_bits(), "{percent}");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::HumidityPercentOutOfRange,
            "{percent}"
        );
    }
    for percent in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = RelativeHumidity::try_percent(percent).expect_err("non-finite");
        assert_eq!(error.reason(), AtmosphereErrorReason::NotFinite);
    }
}

#[test]
fn try_fraction_reads_a_sensor_reading_within_five_points_as_saturated_or_dry_air() {
    // 1.0 is exactly 100 %, and a fraction up to 0.05 above it is saturated air.
    for fraction in [1.0, 1.000_000_000_000_001, 1.0001, 1.03, 1.049_999_9, 1.05] {
        let humidity = RelativeHumidity::try_fraction(fraction)
            .unwrap_or_else(|e| panic!("{fraction} is a saturated-air reading: {e}"));
        assert_eq!(
            humidity.percent().to_bits(),
            100.0_f64.to_bits(),
            "{fraction}"
        );
        assert_eq!(humidity, rh(100.0), "{fraction}");
    }
    for fraction in [-0.0, -1e-300, -0.0001, -0.03, -0.049_999_9, -0.05] {
        let humidity = RelativeHumidity::try_fraction(fraction)
            .unwrap_or_else(|e| panic!("{fraction} is a dry-air reading: {e}"));
        assert_eq!(
            humidity.percent().to_bits(),
            0.0_f64.to_bits(),
            "{fraction}"
        );
        assert_eq!(humidity, RelativeHumidity::DRY, "{fraction}");
    }
}

#[test]
fn try_fraction_refuses_a_reading_beyond_the_sensor_band_unchanged() {
    // A percent sent where a fraction belongs (50.0, 100.0) stays an error.
    for fraction in [
        1.050_000_1,
        1.06,
        1.5,
        50.0,
        100.0,
        -0.050_000_1,
        -0.1,
        -1.0,
    ] {
        let error = RelativeHumidity::try_fraction(fraction).expect_err("beyond the band");
        assert_eq!(error.field(), "relative_humidity", "{fraction}");
        assert_eq!(error.value().to_bits(), fraction.to_bits(), "{fraction}");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::HumidityFractionOutOfRange,
            "{fraction}"
        );
    }
    for fraction in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error = RelativeHumidity::try_fraction(fraction).expect_err("non-finite");
        assert_eq!(error.reason(), AtmosphereErrorReason::NotFinite);
    }
}

#[test]
fn a_reading_in_the_band_is_the_same_station_and_zone_air_as_the_end_it_reads_as() {
    // The snapped humidity is the end itself, so the station air (density, sound speed, ratio)
    // and a zone built at the reading are the end's to the bit.
    for (reading, end) in [(101.7, 100.0), (105.0, 100.0), (-2.0, 0.0), (-5.0, 0.0)] {
        for (altitude_m, temp_c, pressure_hpa) in [
            (0.0, 15.0, 1013.25),
            (1_500.0, 30.0, 830.0),
            (0.0, -10.0, 1030.0),
        ] {
            let at_reading = StationAir::new(altitude_m, temp_c, pressure_hpa, rh(reading))
                .expect("a reading in the band is a station");
            let at_end = StationAir::new(altitude_m, temp_c, pressure_hpa, rh(end)).unwrap();
            assert_eq!(at_reading, at_end, "{reading} -> {end} at {temp_c} C");
        }
        assert!(AtmoSock::new(vec![(15.0, 1013.25, reading, 300.0)])
            .validate()
            .is_ok());
    }
    // The band does not relax impossible air: saturated air at the boiling point read as 102 %
    // is refused exactly as 100 % is, at the station and in a zone.
    for reading in [100.0, 101.7, 105.0] {
        let error = StationAir::new(0.0, 100.0, 1013.25, rh(reading)).expect_err("steam");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::ImpossibleAir,
            "{reading}"
        );
        let error = AtmoSock::new(vec![(100.0, 1013.25, reading, 300.0)])
            .validate()
            .expect_err("steam");
        assert_eq!(
            error.reason(),
            AtmosphereErrorReason::ImpossibleAir,
            "{reading}"
        );
    }
}

// =================================================================================================
// AtmosphereError
// =================================================================================================

#[test]
fn atmosphere_error_names_its_field_and_value_and_can_be_relabelled_at_a_boundary() {
    let error = RelativeHumidity::try_percent(150.0).unwrap_err();
    let text = error.to_string();
    assert!(text.contains("humidity_percent"), "{text}");
    assert!(text.contains("150"), "{text}");
    assert!(text.contains("0") && text.contains("100"), "{text}");

    let relabelled = error.clone().at_field("$.atmosphere.relative_humidity");
    assert_eq!(relabelled.field(), "$.atmosphere.relative_humidity");
    assert_eq!(relabelled.value().to_bits(), error.value().to_bits());
    assert_eq!(relabelled.reason(), error.reason());

    // It is a std error, so `?` can lift it into Box<dyn Error> and friends.
    let boxed: Box<dyn std::error::Error> = Box::new(relabelled);
    assert!(boxed
        .to_string()
        .starts_with("$.atmosphere.relative_humidity"));
}

// =================================================================================================
// One vapour mole fraction (H24 merge)
// =================================================================================================

/// 0.49.0's inline vapour mole fraction from `moist_speed_of_sound`, verbatim, before its
/// `.min(1.0)` cap.
fn inline_mole_fraction_0_49(temp_k: f64, pressure_pa: f64, humidity_percent: f64) -> f64 {
    let humidity_clamped = humidity_percent.clamp(0.0, 100.0);
    let temp_c = temp_k - 273.15;
    let p_sv_hpa = enhanced_saturation_vapor_pressure(temp_k);
    let f = enhanced_enhancement_factor(pressure_pa, temp_c);
    let vapor_pressure_pa = humidity_clamped / 100.0 * f * p_sv_hpa * 100.0;
    vapor_pressure_pa / pressure_pa.max(f64::MIN_POSITIVE)
}

/// 0.49.0's `moist_speed_of_sound`, verbatim.
fn moist_speed_of_sound_0_49(temp_k: f64, pressure_pa: f64, humidity_percent: f64) -> f64 {
    let mole_fraction_vapor =
        inline_mole_fraction_0_49(temp_k, pressure_pa, humidity_percent).min(1.0);
    let gamma_moist = GAMMA * (1.0 - mole_fraction_vapor * 0.062);
    let r_moist = R_AIR * (1.0 + 0.378 * mole_fraction_vapor);
    (gamma_moist * r_moist * temp_k).sqrt()
}

#[test]
fn merged_vapor_mole_fraction_is_bit_identical_to_the_old_inline_copy() {
    let mut checked = 0usize;
    for t_step in 0..=85 {
        let temp_c = -40.0 + f64::from(t_step);
        let temp_k = temp_c + 273.15;
        for p_step in 0..=60 {
            let pressure_pa = (500.0 + 10.0 * f64::from(p_step)) * 100.0;
            for h_step in 0..=20 {
                let humidity = 5.0 * f64::from(h_step);
                let merged = moist_air_vapor_mole_fraction(temp_k, pressure_pa, humidity);
                let inline = inline_mole_fraction_0_49(temp_k, pressure_pa, humidity);
                assert_eq!(
                    merged.to_bits(),
                    inline.to_bits(),
                    "x_v at {temp_c} C, {pressure_pa} Pa, {humidity}%"
                );
                assert_eq!(
                    moist_speed_of_sound(temp_k, pressure_pa, humidity).to_bits(),
                    moist_speed_of_sound_0_49(temp_k, pressure_pa, humidity).to_bits(),
                    "c at {temp_c} C, {pressure_pa} Pa, {humidity}%"
                );
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 86 * 61 * 21);
}

#[test]
fn merged_moist_speed_of_sound_keeps_the_old_total_function_guards() {
    // Out-of-range RH, a zero or negative pressure and NaN inputs reach the legacy public
    // function unvalidated; the merge must keep its clamp, its pressure floor and its x_v cap
    // bit-for-bit.
    let cases = [
        (288.15, 101_325.0, -5.0),
        (288.15, 101_325.0, 150.0),
        (288.15, 0.0, 50.0),
        (288.15, -5.0, 50.0),
        (373.15, 101_325.0, 100.0),
        (500.0, 1_000.0, 100.0),
        (700.0, 101_325.0, 0.0),
        (288.15, 101_325.0, f64::NAN),
        (f64::NAN, 101_325.0, 50.0),
        (288.15, f64::NAN, 50.0),
    ];
    for (temp_k, pressure_pa, humidity) in cases {
        let merged = moist_speed_of_sound(temp_k, pressure_pa, humidity);
        let old = moist_speed_of_sound_0_49(temp_k, pressure_pa, humidity);
        assert_eq!(
            merged.to_bits(),
            old.to_bits(),
            "c at {temp_k} K, {pressure_pa} Pa, {humidity}%: merged {merged} old {old}"
        );
    }
}

#[test]
fn dry_sound_speed_is_the_engines_single_dry_constant() {
    for temp_k in [223.15, 250.0, 288.15, 303.15, 318.15] {
        for pressure_pa in [50_000.0, 101_325.0, 110_000.0] {
            assert_eq!(
                moist_speed_of_sound(temp_k, pressure_pa, 0.0).to_bits(),
                (GAMMA * R_AIR * temp_k).sqrt().to_bits(),
                "{temp_k} K {pressure_pa} Pa"
            );
        }
    }
}

// =================================================================================================
// StationAir
// =================================================================================================

#[test]
fn station_air_new_is_calculate_atmosphere_bit_for_bit() {
    for altitude_m in [-100.0, 0.0, 1_500.0, 4_000.0] {
        for temp_c in TEMPS_C {
            for pressure_hpa in PRESSURES_HPA {
                for humidity in HUMIDITIES {
                    let Ok(station) =
                        StationAir::new(altitude_m, temp_c, pressure_hpa, rh(humidity))
                    else {
                        // Only impossible air may be refused on this grid.
                        let x_v =
                            water_vapor_mole_fraction_uncapped(temp_c, pressure_hpa, humidity);
                        assert!(x_v >= 1.0, "{temp_c} C {pressure_hpa} hPa {humidity}%");
                        continue;
                    };
                    let (density, sound) = calculate_atmosphere(
                        altitude_m,
                        Some(temp_c),
                        Some(pressure_hpa),
                        humidity,
                    );
                    assert_eq!(
                        bits(&[
                            station.density_kg_m3(),
                            station.speed_of_sound_mps(),
                            station.density_ratio(),
                        ]),
                        bits(&[density, sound, density / 1.225]),
                        "{altitude_m} m {temp_c} C {pressure_hpa} hPa {humidity}%"
                    );
                    assert_eq!(
                        bits(&[
                            station.altitude_m(),
                            station.temperature_c(),
                            station.pressure_hpa(),
                            station.humidity().percent(),
                        ]),
                        bits(&[altitude_m, temp_c, pressure_hpa, humidity])
                    );
                }
            }
        }
    }
}

#[test]
fn station_air_accepts_exact_zero_and_one_hundred_percent() {
    for humidity in [0.0, 100.0] {
        let station = StationAir::new(0.0, 15.0, 1013.25, rh(humidity)).expect("possible air");
        assert_eq!(station.humidity().percent().to_bits(), humidity.to_bits());
    }
}

#[test]
fn station_air_refuses_non_finite_and_unphysical_state_naming_the_field() {
    let humidity = RelativeHumidity::DEFAULT;
    let cases: [(f64, f64, f64, &str, AtmosphereErrorReason); 11] = [
        (
            f64::NAN,
            15.0,
            1013.25,
            "altitude_m",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            f64::INFINITY,
            15.0,
            1013.25,
            "altitude_m",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            0.0,
            f64::NAN,
            1013.25,
            "temperature_c",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            0.0,
            f64::NEG_INFINITY,
            1013.25,
            "temperature_c",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            0.0,
            -273.15,
            1013.25,
            "temperature_c",
            AtmosphereErrorReason::AtOrBelowAbsoluteZero,
        ),
        (
            0.0,
            -300.0,
            1013.25,
            "temperature_c",
            AtmosphereErrorReason::AtOrBelowAbsoluteZero,
        ),
        (
            0.0,
            15.0,
            f64::NAN,
            "pressure_hpa",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            0.0,
            15.0,
            f64::INFINITY,
            "pressure_hpa",
            AtmosphereErrorReason::NotFinite,
        ),
        (
            0.0,
            15.0,
            0.0,
            "pressure_hpa",
            AtmosphereErrorReason::NonPositivePressure,
        ),
        (
            0.0,
            15.0,
            -0.0,
            "pressure_hpa",
            AtmosphereErrorReason::NonPositivePressure,
        ),
        (
            0.0,
            15.0,
            -1013.25,
            "pressure_hpa",
            AtmosphereErrorReason::NonPositivePressure,
        ),
    ];
    for (altitude_m, temp_c, pressure_hpa, field, reason) in cases {
        let error = StationAir::new(altitude_m, temp_c, pressure_hpa, humidity)
            .expect_err("unphysical state");
        assert_eq!(error.field(), field, "{altitude_m} {temp_c} {pressure_hpa}");
        assert_eq!(
            error.reason(),
            reason,
            "{altitude_m} {temp_c} {pressure_hpa}"
        );
        let error = StationAir::with_density_ratio(altitude_m, temp_c, pressure_hpa, 1.0, humidity)
            .expect_err("unphysical state");
        assert_eq!(error.field(), field, "{altitude_m} {temp_c} {pressure_hpa}");
        assert_eq!(
            error.reason(),
            reason,
            "{altitude_m} {temp_c} {pressure_hpa}"
        );
    }
    // Just above absolute zero is a state the model can evaluate.
    assert!(StationAir::new(0.0, -273.0, 1013.25, RelativeHumidity::DRY).is_ok());
}

#[test]
fn station_air_refuses_the_impossible_air_that_density_altitude_refuses() {
    // The two cases the plan names, refused at every RH-carrying constructor.
    for (temp_c, pressure_hpa) in [(100.0, 1013.25), (95.0, 500.0)] {
        let error = StationAir::new(0.0, temp_c, pressure_hpa, rh(100.0)).expect_err("steam");
        assert_eq!(error.field(), "humidity_percent");
        assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
        assert_eq!(error.value().to_bits(), 100.0_f64.to_bits());
        let error = StationAir::with_density_ratio(0.0, temp_c, pressure_hpa, 0.6, rh(100.0))
            .expect_err("steam");
        assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
    }

    // Parity with the service's predicate (atmosphere_service.rs: NaN or >= 1 is refused) over
    // a grid that crosses the boundary, up past water's critical point.
    let temps_c = [
        -40.0, 0.0, 20.0, 40.0, 60.0, 80.0, 90.0, 95.0, 99.0, 100.0, 101.0, 110.0, 150.0, 200.0,
        300.0, 373.0, 374.0, 380.0, 500.0,
    ];
    let pressures_hpa = [300.0, 500.0, 700.0, 850.0, 1013.25, 1100.0, 2000.0];
    let humidities = [0.0, 0.5, 1.0, 50.0, 99.0, 100.0];
    let (mut refused, mut accepted) = (0, 0);
    for temp_c in temps_c {
        for pressure_hpa in pressures_hpa {
            for humidity in humidities {
                let x_v = water_vapor_mole_fraction_uncapped(temp_c, pressure_hpa, humidity);
                let service_refuses = x_v.is_nan() || x_v >= 1.0;
                match StationAir::new(0.0, temp_c, pressure_hpa, rh(humidity)) {
                    Ok(_) => {
                        assert!(
                            !service_refuses,
                            "accepted impossible air: {temp_c} C {pressure_hpa} hPa {humidity}%"
                        );
                        accepted += 1;
                    }
                    Err(error) => {
                        assert!(
                            service_refuses,
                            "refused possible air: {temp_c} C {pressure_hpa} hPa {humidity}%: \
                             {error}"
                        );
                        assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
                        // Above water's critical point the vapour formula itself is undefined
                        // (NaN at any RH), which is the temperature's fault, not the humidity's.
                        let field = if x_v.is_nan() {
                            "temperature_c"
                        } else {
                            "humidity_percent"
                        };
                        assert_eq!(error.field(), field, "{temp_c} C {pressure_hpa} hPa");
                        refused += 1;
                    }
                }
            }
        }
    }
    assert!(
        refused > 20 && accepted > 20,
        "the grid crosses the boundary"
    );
}

#[test]
fn with_density_ratio_keeps_the_callers_ratio_verbatim() {
    for ratio in [1e-9, 0.5, 0.95, 1.0, 1.25, 1.999_999] {
        for humidity in [0.0, 37.5, 100.0] {
            let station = StationAir::with_density_ratio(250.0, 21.0, 980.0, ratio, rh(humidity))
                .expect("valid ratio");
            assert_eq!(station.density_ratio().to_bits(), ratio.to_bits());
            assert_eq!(station.density_kg_m3().to_bits(), (ratio * 1.225).to_bits());
            // The station sound speed is the RH-carrying one, the same as `new` would give.
            let reference = StationAir::new(250.0, 21.0, 980.0, rh(humidity)).unwrap();
            assert_eq!(
                station.speed_of_sound_mps().to_bits(),
                reference.speed_of_sound_mps().to_bits()
            );
            assert_eq!(station.humidity(), rh(humidity));
            assert_eq!(
                bits(&[
                    station.altitude_m(),
                    station.temperature_c(),
                    station.pressure_hpa()
                ]),
                bits(&[250.0, 21.0, 980.0])
            );
        }
    }
}

#[test]
fn with_density_ratio_refuses_ratios_outside_zero_to_two() {
    for ratio in [
        0.0,
        -0.0,
        -1.0,
        2.0,
        2.5,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        let error =
            StationAir::with_density_ratio(0.0, 15.0, 1013.25, ratio, RelativeHumidity::DRY)
                .expect_err("bad ratio");
        assert_eq!(error.field(), "density_ratio", "{ratio}");
        let reason = if ratio.is_finite() {
            AtmosphereErrorReason::DensityRatioOutOfRange
        } else {
            AtmosphereErrorReason::NotFinite
        };
        assert_eq!(error.reason(), reason, "{ratio}");
    }
}

#[test]
fn the_solver_station_is_station_air_new_and_refuses_what_0_49_clamped() {
    // TrajectorySolver resolves its station through StationAir::new (plan C4), the interim
    // crate-private `legacy_clamped` constructor is gone. A valid station is still exactly
    // calculate_atmosphere at the given humidity, to the bit; the inputs 0.49.0 flew clamped or
    // capped (an RH beyond the sensor band around [0, 100], impossible air, a pressure of 0) are
    // refused instead.
    for (altitude_m, temp_c, pressure_hpa, humidity) in [
        (0.0, 15.0, 1013.25, 50.0),
        (1_500.0, 25.0, 830.0, 0.0),
        (0.0, 30.0, 1013.25, 100.0),
    ] {
        let station = StationAir::new(altitude_m, temp_c, pressure_hpa, rh(humidity))
            .expect("a valid station");
        let (density, sound) =
            calculate_atmosphere(altitude_m, Some(temp_c), Some(pressure_hpa), humidity);
        assert_eq!(
            bits(&[
                station.density_kg_m3(),
                station.speed_of_sound_mps(),
                station.density_ratio()
            ]),
            bits(&[density, sound, density / 1.225]),
            "{temp_c} C {pressure_hpa} hPa {humidity}%"
        );
        assert_eq!(station.humidity().percent().to_bits(), humidity.to_bits());
    }
    for humidity in [150.0, -5.000_000_1, f64::NAN] {
        let error = RelativeHumidity::try_percent(humidity).expect_err("not a humidity");
        assert_eq!(error.field(), "humidity_percent", "{humidity}");
    }
    let error = StationAir::new(0.0, 100.0, 1013.25, rh(100.0)).expect_err("impossible air");
    assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
    let error = StationAir::new(0.0, 15.0, 0.0, rh(50.0)).expect_err("no pressure");
    assert_eq!(error.reason(), AtmosphereErrorReason::NonPositivePressure);
}

// =================================================================================================
// KernelAtmosphere
// =================================================================================================

#[test]
fn uniform_accepts_todays_direct_mode_range() {
    for (density, sound) in [(1e-9, 200.000_001), (1.225, 340.294), (1.999_999, 1_500.0)] {
        match KernelAtmosphere::uniform(density, sound).expect("valid direct air") {
            KernelAtmosphere::Uniform {
                density_kg_m3,
                speed_of_sound_mps,
            } => {
                assert_eq!(density_kg_m3.to_bits(), density.to_bits());
                assert_eq!(speed_of_sound_mps.to_bits(), sound.to_bits());
            }
            other => panic!("expected Uniform, got {other:?}"),
        }
    }
}

#[test]
fn uniform_rejects_what_the_direct_mode_kernels_reject() {
    let bad_density = [0.0, -0.0, -1.0, 2.0, 2.5, f64::NAN, f64::INFINITY];
    for density in bad_density {
        let error = KernelAtmosphere::uniform(density, 340.0).expect_err("bad density");
        assert_eq!(error.field(), "density_kg_m3", "{density}");
        let reason = if density.is_finite() {
            AtmosphereErrorReason::UniformDensityOutOfRange
        } else {
            AtmosphereErrorReason::NotFinite
        };
        assert_eq!(error.reason(), reason, "{density}");
    }
    for sound in [200.0, 199.0, 0.0, -340.0, f64::NAN, f64::INFINITY] {
        let error = KernelAtmosphere::uniform(1.2, sound).expect_err("bad sound speed");
        assert_eq!(error.field(), "speed_of_sound_mps", "{sound}");
        let reason = if sound.is_finite() {
            AtmosphereErrorReason::UniformSpeedOfSoundOutOfRange
        } else {
            AtmosphereErrorReason::NotFinite
        };
        assert_eq!(error.reason(), reason, "{sound}");
    }
}

#[test]
fn from_legacy_params_maps_the_direct_sentinel_to_uniform() {
    let atmosphere =
        KernelAtmosphere::from_legacy_params((1.1, 330.0, 0.0, 0.0), RelativeHumidity::DEFAULT)
            .expect("direct sentinel");
    assert_eq!(atmosphere, KernelAtmosphere::uniform(1.1, 330.0).unwrap());
    // -0.0 sentinels compare equal to 0.0, as the kernels' `== 0.0` tests do.
    let atmosphere =
        KernelAtmosphere::from_legacy_params((1.1, 330.0, -0.0, -0.0), RelativeHumidity::DRY)
            .expect("direct sentinel");
    assert_eq!(atmosphere, KernelAtmosphere::uniform(1.1, 330.0).unwrap());
    // A sentinel with out-of-range direct values is refused, not reinterpreted.
    for params in [
        (2.5, 340.0, 0.0, 0.0),
        (1.2, 150.0, 0.0, 0.0),
        (0.0, 340.0, 0.0, 0.0),
    ] {
        assert!(
            KernelAtmosphere::from_legacy_params(params, RelativeHumidity::DRY).is_err(),
            "{params:?}"
        );
    }
}

#[test]
fn from_legacy_params_maps_standard_mode_to_a_station_with_the_ratio_verbatim() {
    let humidity = rh(37.5);
    let atmosphere = KernelAtmosphere::from_legacy_params((1_500.0, 25.0, 830.0, 0.93), humidity)
        .expect("standard mode");
    assert_eq!(
        atmosphere,
        KernelAtmosphere::Station(
            StationAir::with_density_ratio(1_500.0, 25.0, 830.0, 0.93, humidity).unwrap()
        )
    );
    // A nonpositive ratio means "not supplied" and falls back to 1.0 (MBA-1157).
    for missing in [0.0, -0.0, -1.0, -1e300] {
        let KernelAtmosphere::Station(station) =
            KernelAtmosphere::from_legacy_params((0.0, 15.0, 1013.25, missing), humidity)
                .expect("missing ratio")
        else {
            panic!("standard mode must give a station");
        };
        assert_eq!(
            station.density_ratio().to_bits(),
            1.0_f64.to_bits(),
            "{missing}"
        );
    }
    // A supplied ratio at or above 2 is refused, as fast_integrate refuses it today.
    for ratio in [2.0, 3.0, f64::NAN, f64::INFINITY] {
        let error = KernelAtmosphere::from_legacy_params((0.0, 15.0, 1013.25, ratio), humidity)
            .expect_err("overlarge ratio");
        assert_eq!(error.field(), "density_ratio", "{ratio}");
    }
    // A zero station pressure that is not the direct sentinel is refused.
    let error = KernelAtmosphere::from_legacy_params((0.0, 15.0, 0.0, 1.0), humidity)
        .expect_err("zero pressure");
    assert_eq!(error.field(), "pressure_hpa");
    // Impossible air is refused through the station constructor.
    let error = KernelAtmosphere::from_legacy_params((0.0, 100.0, 1013.25, 0.6), rh(100.0))
        .expect_err("steam");
    assert_eq!(error.reason(), AtmosphereErrorReason::ImpossibleAir);
}

// =================================================================================================
// air_at
// =================================================================================================

#[test]
fn air_at_uniform_returns_the_callers_air_everywhere_and_ignores_zones() {
    let atmosphere = KernelAtmosphere::uniform(1.1, 330.0).unwrap();
    let zones = AtmoSock::new(vec![
        (35.0, 900.0, 90.0, 100.0),
        (-5.0, 1050.0, 10.0, 500.0),
    ]);
    for zone in [None, Some(&zones)] {
        for downrange_m in [0.0, 99.0, 500.0, 2_000.0] {
            for altitude_m in [-200.0, 0.0, 1_000.0, 12_000.0] {
                let air = air_at(&atmosphere, zone, downrange_m, altitude_m);
                assert_eq!(air.density_kg_m3.to_bits(), 1.1_f64.to_bits());
                assert_eq!(air.speed_of_sound_mps.to_bits(), 330.0_f64.to_bits());
                assert_eq!(
                    air.temperature_k.to_bits(),
                    (330.0 * 330.0 / (GAMMA * R_AIR)).to_bits()
                );
                assert_eq!(air.pressure_pa.to_bits(), 0.0_f64.to_bits());
            }
        }
    }
}

#[test]
fn air_at_station_is_the_legacy_humid_lookup_bit_for_bit() {
    for base_alt in [0.0, 1_500.0] {
        for temp_c in [-10.0, 15.0, 35.0] {
            for pressure_hpa in [850.0, 1013.25] {
                for humidity in [0.0, 50.0, 100.0] {
                    let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity))
                        .expect("possible air");
                    let atmosphere = KernelAtmosphere::Station(station);
                    for offset in ALT_OFFSETS_M {
                        let altitude_m = base_alt + offset;
                        let (density, sound) = get_local_atmosphere_humid(
                            altitude_m,
                            base_alt,
                            temp_c,
                            pressure_hpa,
                            station.density_ratio(),
                            humidity,
                        );
                        // The lapse at the station's own vapour mole fraction, held aloft
                        // (S8: the virtual-temperature lapse; dry at 0 %).
                        let (temp_k, pressure_pa, density_again) =
                            local_temp_pressure_density_with_gravity(
                                altitude_m,
                                base_alt,
                                temp_c,
                                pressure_hpa,
                                station.density_ratio(),
                                governing_vapor_mole_fraction(temp_c, pressure_hpa, humidity),
                                Gravity::EARTH,
                            );
                        assert_eq!(density.to_bits(), density_again.to_bits());
                        for downrange_m in [0.0, 99.0, 500.0] {
                            let air = air_at(&atmosphere, None, downrange_m, altitude_m);
                            assert_eq!(
                                local_bits(air),
                                [
                                    density.to_bits(),
                                    sound.to_bits(),
                                    temp_k.to_bits(),
                                    pressure_pa.to_bits()
                                ],
                                "{base_alt} m {temp_c} C {pressure_hpa} hPa {humidity}% \
                                 +{offset} m"
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn air_at_with_gravity_matches_the_legacy_gravity_aware_lookup() {
    let moon = Gravity::new(1.62).unwrap();
    let zero = Gravity::new(0.0).unwrap();
    let station = StationAir::new(300.0, 20.0, 950.0, rh(62.5)).unwrap();
    let atmosphere = KernelAtmosphere::Station(station);
    for gravity in [Gravity::EARTH, moon, zero] {
        for offset in ALT_OFFSETS_M {
            let altitude_m = 300.0 + offset;
            let (density, sound) = get_local_atmosphere_humid_with_gravity(
                altitude_m,
                300.0,
                20.0,
                950.0,
                station.density_ratio(),
                62.5,
                gravity,
            );
            let air = air_at_with_gravity(&atmosphere, None, 0.0, altitude_m, gravity);
            assert_eq!(
                [
                    air.density_kg_m3.to_bits(),
                    air.speed_of_sound_mps.to_bits()
                ],
                [density.to_bits(), sound.to_bits()],
                "g {} +{offset} m",
                gravity.magnitude()
            );
        }
    }
    // The public entry is the Earth-gravity case.
    for offset in ALT_OFFSETS_M {
        assert_eq!(
            local_bits(air_at(&atmosphere, None, 10.0, 300.0 + offset)),
            local_bits(air_at_with_gravity(
                &atmosphere,
                None,
                10.0,
                300.0 + offset,
                Gravity::EARTH
            ))
        );
    }
}

/// 0.49.0's `TrajectorySolver::calculate_acceleration` air block, verbatim but for `self`:
/// the zone (or station) base, then `get_local_atmosphere_humid_with_gravity`.
#[allow(clippy::too_many_arguments)]
fn trajectory_solver_air_0_49(
    atmo_sock: Option<&AtmoSock>,
    position_x: f64,
    local_alt: f64,
    station_altitude: f64,
    base_temp_c: f64,
    base_press_hpa: f64,
    station_ratio: f64,
    station_humidity: f64,
    gravity: Gravity,
) -> (f64, f64) {
    let (drag_base_temp_c, drag_base_press_hpa, drag_base_ratio, drag_humidity_percent) =
        if let Some(sock) = atmo_sock {
            let (zone_temp_c, zone_press_hpa, zone_humidity) = sock.atmo_for_range(position_x);
            let zone_base_ratio =
                calculate_air_density_cimp(zone_temp_c, zone_press_hpa, zone_humidity) / 1.225;
            (zone_temp_c, zone_press_hpa, zone_base_ratio, zone_humidity)
        } else {
            (base_temp_c, base_press_hpa, station_ratio, station_humidity)
        };
    get_local_atmosphere_humid_with_gravity(
        local_alt,
        station_altitude,
        drag_base_temp_c,
        drag_base_press_hpa,
        drag_base_ratio,
        drag_humidity_percent,
        gravity,
    )
}

#[test]
fn air_at_reproduces_the_0_49_trajectory_solver_air_with_and_without_zones() {
    // Zone RH differs from the station's, the zones cross the downrange grid, and one zone
    // carries an out-of-range RH that the internal guards still clamp until AtmoSock::validate
    // exists.
    let zones = AtmoSock::new(vec![
        (20.0, 1000.0, 20.0, 100.0),
        (30.0, 990.0, 80.0, 450.0),
        (5.0, 1020.0, 120.0, 800.0),
    ]);
    let moon = Gravity::new(1.62).unwrap();
    for base_alt in [0.0, 1_500.0] {
        for temp_c in [-10.0, 15.0, 35.0] {
            for pressure_hpa in [850.0, 1013.25] {
                for humidity in [0.0, 50.0, 100.0] {
                    let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity))
                        .expect("possible air");
                    let atmosphere = KernelAtmosphere::Station(station);
                    for sock in [None, Some(&zones)] {
                        for gravity in [Gravity::EARTH, moon] {
                            for downrange_m in [0.0, 99.0, 100.0, 500.0, 2_000.0] {
                                for offset in [0.0, 50.0, 1_000.0, -200.0] {
                                    let altitude_m = base_alt + offset;
                                    let (density, sound) = trajectory_solver_air_0_49(
                                        sock,
                                        downrange_m,
                                        altitude_m,
                                        base_alt,
                                        temp_c,
                                        pressure_hpa,
                                        station.density_ratio(),
                                        humidity,
                                        gravity,
                                    );
                                    let air = air_at_with_gravity(
                                        &atmosphere,
                                        sock,
                                        downrange_m,
                                        altitude_m,
                                        gravity,
                                    );
                                    assert_eq!(
                                        [
                                            air.density_kg_m3.to_bits(),
                                            air.speed_of_sound_mps.to_bits()
                                        ],
                                        [density.to_bits(), sound.to_bits()],
                                        "{base_alt} m {temp_c} C {pressure_hpa} hPa \
                                         {humidity}% zones={} x={downrange_m} +{offset} m",
                                        sock.is_some()
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn air_at_zone_base_is_the_cipm_ratio_at_the_zone_to_the_bit() {
    // A dense sweep of single-zone socks, so a last-bit change in the zone density ratio (a
    // reordered division, a different standard-density constant) cannot hide behind a lucky
    // handful of zones. Compared against the 0.49.0 solver block, which divides the CIPM
    // density by 1.225.
    let station = StationAir::new(0.0, 15.0, 1013.25, RelativeHumidity::DEFAULT).unwrap();
    let atmosphere = KernelAtmosphere::Station(station);
    let mut checked = 0usize;
    for t_step in 0..=10 {
        let zone_temp_c = -10.0 + 5.0 * f64::from(t_step);
        for p_step in 0..=8 {
            let zone_press_hpa = 850.0 + 25.0 * f64::from(p_step);
            for h_step in 0..=8 {
                let zone_humidity = 12.5 * f64::from(h_step);
                let zones = AtmoSock::new(vec![(zone_temp_c, zone_press_hpa, zone_humidity, 1e9)]);
                for offset in [0.0, 120.0, -60.0] {
                    let (density, sound) = trajectory_solver_air_0_49(
                        Some(&zones),
                        250.0,
                        offset,
                        0.0,
                        15.0,
                        1013.25,
                        station.density_ratio(),
                        50.0,
                        Gravity::EARTH,
                    );
                    let air = air_at(&atmosphere, Some(&zones), 250.0, offset);
                    assert_eq!(
                        [
                            air.density_kg_m3.to_bits(),
                            air.speed_of_sound_mps.to_bits()
                        ],
                        [density.to_bits(), sound.to_bits()],
                        "zone {zone_temp_c} C {zone_press_hpa} hPa {zone_humidity}% +{offset} m"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert_eq!(checked, 11 * 9 * 9 * 3);
}

#[test]
fn air_at_selects_the_zone_base_by_downrange_and_lapses_from_the_station_altitude() {
    let station = StationAir::new(500.0, 10.0, 950.0, rh(30.0)).unwrap();
    let atmosphere = KernelAtmosphere::Station(station);
    let zones = AtmoSock::new(vec![(25.0, 940.0, 70.0, 300.0), (0.0, 960.0, 5.0, 900.0)]);
    for (downrange_m, (zt, zp, zh)) in [
        (0.0, (25.0, 940.0, 70.0)),
        (299.0, (25.0, 940.0, 70.0)),
        (300.0, (0.0, 960.0, 5.0)),
        (5_000.0, (0.0, 960.0, 5.0)),
    ] {
        let ratio = calculate_air_density_cimp(zt, zp, zh) / 1.225;
        for altitude_m in [500.0, 800.0, 300.0] {
            let (density, sound) = get_local_atmosphere_humid(altitude_m, 500.0, zt, zp, ratio, zh);
            let air = air_at(&atmosphere, Some(&zones), downrange_m, altitude_m);
            assert_eq!(
                [
                    air.density_kg_m3.to_bits(),
                    air.speed_of_sound_mps.to_bits()
                ],
                [density.to_bits(), sound.to_bits()],
                "x={downrange_m} alt={altitude_m}"
            );
        }
    }
}

#[test]
fn an_empty_zone_set_means_the_station() {
    let empty = AtmoSock::new(Vec::new());
    for humidity in [0.0, 50.0, 100.0] {
        let station = StationAir::new(1_500.0, 25.0, 830.0, rh(humidity)).unwrap();
        let atmosphere = KernelAtmosphere::Station(station);
        for downrange_m in [0.0, 500.0] {
            for altitude_m in [1_500.0, 2_500.0, 1_300.0] {
                assert_eq!(
                    local_bits(air_at(&atmosphere, Some(&empty), downrange_m, altitude_m)),
                    local_bits(air_at(&atmosphere, None, downrange_m, altitude_m)),
                    "{humidity}% x={downrange_m} alt={altitude_m}"
                );
            }
        }
    }
}

#[test]
fn air_at_flies_the_station_humidity_into_the_sound_speed() {
    // Density is held by the station ratio; the RH reaches the sound speed through the one
    // vapour mole fraction, so a humid station is faster than a dry one at the same state.
    let dry =
        StationAir::with_density_ratio(0.0, 30.0, 1013.25, 1.0, RelativeHumidity::DRY).unwrap();
    let humid = StationAir::with_density_ratio(0.0, 30.0, 1013.25, 1.0, rh(100.0)).unwrap();
    let dry_air = air_at(&KernelAtmosphere::Station(dry), None, 0.0, 0.0);
    let humid_air = air_at(&KernelAtmosphere::Station(humid), None, 0.0, 0.0);
    assert_eq!(
        dry_air.density_kg_m3.to_bits(),
        humid_air.density_kg_m3.to_bits()
    );
    assert!(humid_air.speed_of_sound_mps > dry_air.speed_of_sound_mps * 1.004);
    assert_eq!(
        dry_air.speed_of_sound_mps.to_bits(),
        (GAMMA * R_AIR * dry_air.temperature_k).sqrt().to_bits()
    );
    // At the station altitude the local sound speed is the station's own.
    assert_eq!(
        humid_air.speed_of_sound_mps.to_bits(),
        humid.speed_of_sound_mps().to_bits()
    );
}

#[test]
fn atmo_sock_validate_accepts_valid_and_empty_zone_sets() {
    assert!(AtmoSock::new(Vec::new()).validate().is_ok());
    assert!(AtmoSock::new(vec![
        (15.0, 1013.25, 0.0, 300.0),
        (-5.0, 900.0, 100.0, 1_000.0),
        (35.0, 1013.25, 50.0, f64::INFINITY),
    ])
    .validate()
    .is_ok());
}

#[test]
fn atmo_sock_validate_refuses_each_bad_zone_field_by_name_without_clamping() {
    use AtmosphereErrorReason as Reason;
    let cases: [(AtmoSegment, &str, Reason); 12] = [
        (
            (15.0, 1013.25, 150.0, 300.0),
            "humidity_percent",
            Reason::HumidityPercentOutOfRange,
        ),
        (
            (15.0, 1013.25, 105.000_000_1, 300.0),
            "humidity_percent",
            Reason::HumidityPercentOutOfRange,
        ),
        (
            (15.0, 1013.25, -5.000_000_1, 300.0),
            "humidity_percent",
            Reason::HumidityPercentOutOfRange,
        ),
        (
            (15.0, 1013.25, f64::NAN, 300.0),
            "humidity_percent",
            Reason::NotFinite,
        ),
        (
            (-273.15, 1013.25, 50.0, 300.0),
            "temperature_c",
            Reason::AtOrBelowAbsoluteZero,
        ),
        (
            (f64::NAN, 1013.25, 50.0, 300.0),
            "temperature_c",
            Reason::NotFinite,
        ),
        (
            (15.0, 0.0, 50.0, 300.0),
            "pressure_hpa",
            Reason::NonPositivePressure,
        ),
        (
            (15.0, f64::INFINITY, 50.0, 300.0),
            "pressure_hpa",
            Reason::NotFinite,
        ),
        // The density-altitude service's impossible air: saturated at 100 C and 1013.25 hPa.
        (
            (100.0, 1013.25, 100.0, 300.0),
            "humidity_percent",
            Reason::ImpossibleAir,
        ),
        (
            (15.0, 1013.25, 50.0, 0.0),
            "until_distance_m",
            Reason::NonPositiveZoneThreshold,
        ),
        (
            (15.0, 1013.25, 50.0, -5.0),
            "until_distance_m",
            Reason::NonPositiveZoneThreshold,
        ),
        (
            (15.0, 1013.25, 50.0, f64::NAN),
            "until_distance_m",
            Reason::NonPositiveZoneThreshold,
        ),
    ];
    for (zone, field, reason) in cases {
        let error = AtmoSock::new(vec![zone]).validate().unwrap_err();
        assert_eq!(
            error.field(),
            format!("atmo_sock.zones[0].{field}"),
            "{zone:?}"
        );
        assert_eq!(error.reason(), reason, "{zone:?}");
        assert!(
            error.value().to_bits() == zone.0.to_bits()
                || error.value().to_bits() == zone.1.to_bits()
                || error.value().to_bits() == zone.2.to_bits()
                || error.value().to_bits() == zone.3.to_bits(),
            "the refused value is reported unchanged: {error}"
        );
    }
    // The index is the zone's place in threshold order: a valid near zone comes first.
    let error = AtmoSock::new(vec![
        (15.0, 1013.25, 150.0, 900.0),
        (15.0, 1013.25, 50.0, 300.0),
    ])
    .validate()
    .unwrap_err();
    assert_eq!(error.field(), "atmo_sock.zones[1].humidity_percent");
    assert_eq!(
        error.to_string(),
        "atmo_sock.zones[1].humidity_percent = 150: relative humidity must be a percent from 0 \
         to 100"
    );
}

#[test]
fn an_empty_atmo_sock_lookup_falls_back_at_the_one_default_humidity() {
    // air_at never asks an empty sock (it means the station); the raw lookup's own fallback is
    // the engine's default humidity, not a 0 % of its own.
    let (_, _, humidity) = AtmoSock::new(Vec::new()).atmo_for_range(250.0);
    assert_eq!(humidity.to_bits(), DEFAULT_HUMIDITY_PERCENT.to_bits());
}

#[test]
#[allow(deprecated)] // the dry legacy lookup is what is under test
fn the_deprecated_dry_lookup_is_the_shared_composition_at_zero_humidity() {
    // One dry constant: get_local_atmosphere is get_local_atmosphere_humid at 0 %, to the bit,
    // over heights below, at and far above the base, so a dry path can only differ from a humid
    // one by its humidity.
    for (base_alt, base_t, base_p, ratio) in [
        (0.0, 15.0, 1013.25, 1.0),
        (1_500.0, 25.0, 830.0, 0.79),
        (500.0, -10.0, 950.0, 1.05),
    ] {
        for offset in ALT_OFFSETS_M {
            let altitude = base_alt + offset;
            let dry = get_local_atmosphere(altitude, base_alt, base_t, base_p, ratio);
            let humid0 = get_local_atmosphere_humid(altitude, base_alt, base_t, base_p, ratio, 0.0);
            assert_eq!(
                [dry.0.to_bits(), dry.1.to_bits()],
                [humid0.0.to_bits(), humid0.1.to_bits()],
                "{base_alt} m {base_t} C {base_p} hPa, query {altitude} m"
            );
        }
    }
}

#[test]
fn a_terminal_humidity_flag_parses_through_the_one_validator() {
    // The WASM terminal's parser for every humidity flag (plan C4): a percent in [0, 100], a
    // sensor reading within 5 points of it read as the end (S5b), anything farther out refused
    // (never clamped), with the flag named in the error.
    for (raw, percent) in [("0", 0.0), ("100", 100.0), ("37.5", 37.5), (" 50 ", 50.0)] {
        assert_eq!(
            parse_humidity_percent_arg("--humidity", raw)
                .unwrap()
                .to_bits(),
            f64::to_bits(percent),
            "{raw}"
        );
    }
    // -0.0 is dry air with positive-zero bits, as everywhere else.
    assert_eq!(
        parse_humidity_percent_arg("--humidity", "-0")
            .unwrap()
            .to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(
        parse_humidity_percent_arg("--humidity", "150").unwrap_err(),
        "--humidity = 150: relative humidity must be a percent from 0 to 100"
    );
    assert_eq!(
        parse_humidity_percent_arg("--zero-humidity", "-5.0000001").unwrap_err(),
        "--zero-humidity = -5.0000001: relative humidity must be a percent from 0 to 100"
    );
    // A sensor reading within 5 points of the range is the end it reads as, to the bit.
    for (raw, end) in [
        ("100.0001", 100.0),
        ("101.7", 100.0),
        ("105", 100.0),
        ("-0.0001", 0.0),
        ("-2", 0.0),
        ("-5", 0.0),
    ] {
        assert_eq!(
            parse_humidity_percent_arg("--humidity", raw)
                .unwrap_or_else(|e| panic!("{raw}: {e}"))
                .to_bits(),
            f64::to_bits(end),
            "{raw}"
        );
    }
    for raw in ["105.0000001", "-5.0000001", "1e3"] {
        let error = parse_humidity_percent_arg("--humidity", raw).unwrap_err();
        assert!(
            error.starts_with("--humidity = ") && error.contains("percent from 0 to 100"),
            "{raw}: {error}"
        );
    }
    for raw in ["nan", "inf", "-inf"] {
        let error = parse_humidity_percent_arg("--humidity", raw).unwrap_err();
        assert!(error.ends_with("must be finite"), "{raw}: {error}");
    }
    let error = parse_humidity_percent_arg("--humidity", "fifty").unwrap_err();
    assert!(
        error.contains("--humidity = fifty: not a number"),
        "{error}"
    );
}
