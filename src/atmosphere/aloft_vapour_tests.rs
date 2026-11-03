//! The aloft vapour model.
//!
//! Until 0.50.0 the air aloft modelled its water vapour two ways inside one function. The density
//! scaled the base density by `(T_0 · P) / (P_0 · T)`, which keeps the base air's vapour mole
//! fraction `x_v`; the sound speed kept the base RELATIVE humidity and re-derived `x_v` at the
//! local temperature and pressure; and the hydrostatic lapse was dry. 0.50.0 holds the governing
//! air's `x_v` for all three (the well-mixed assumption), with a virtual-temperature lapse: the dry
//! exponent scaled by the moist-air molar-mass ratio `M_m / M_a = 1 - x_v (1 - M_v / M_a)`, the
//! same factor the CIPM-2007 density applies.
//!
//! The governing air is the station's, or inside a downrange [`AtmoSock`] zone the zone's own:
//! that is the base `air_at` lapses from, the one whose density ratio it already scales.
//!
//! - At 0 % `x_v` is 0 and `M_m = M_a`, so dry air is bit-identical to the model before
//!   ([`dry_air_aloft_is_bit_identical_to_the_model_before`]).
//! - At a humidity above 0 the sound speed, the density and the lapse all carry the base `x_v`,
//!   checked against the formulas and against an independent hand calculation.

use super::*;

fn rh(percent: f64) -> RelativeHumidity {
    RelativeHumidity::try_percent(percent).expect("valid test humidity")
}

/// The governing (base) air's vapour mole fraction, through the very conversions the station
/// sound speed uses (`calculate_atmosphere`: K = C + 273.15, Pa = hPa · 100), capped at 1.
fn base_x_v(temp_c: f64, press_hpa: f64, humidity_percent: f64) -> f64 {
    moist_air_vapor_mole_fraction(temp_c + 273.15, press_hpa * 100.0, humidity_percent).min(1.0)
}

/// The first-order moist-air sound speed at a GIVEN vapour mole fraction.
fn sound_speed_at(temp_k: f64, x_v: f64) -> f64 {
    (GAMMA * (1.0 - 0.062 * x_v) * R_AIR * (1.0 + 0.378 * x_v) * temp_k).sqrt()
}

/// The vapour mole fraction a first-order moist sound speed implies at `temp_k`: the root near
/// zero of `(1 - 0.062 x)(1 + 0.378 x) = c² / (1.4 · 287.0531 · T)`, that is of
/// `a x² - b x + q = 0` with `b = 0.316`, `a = 0.062 · 0.378` and `q = c² / (1.4 R T) - 1`, in
/// the cancellation-free form `2q / (b + sqrt(b² - 4aq))`.
fn x_v_from_sound_speed(sound_speed_mps: f64, temp_k: f64) -> f64 {
    let q = sound_speed_mps * sound_speed_mps / (GAMMA * R_AIR * temp_k) - 1.0;
    let (b, a) = (0.378 - 0.062, 0.062 * 0.378);
    2.0 * q / (b + (b * b - 4.0 * a * q).sqrt())
}

/// The moist-air molar-mass ratio `M_m / M_a` at `x_v`.
fn molar_mass_ratio(x_v: f64) -> f64 {
    1.0 - x_v * (1.0 - M_V / M_A)
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    ((actual - expected) / expected).abs()
}

/// Humid stations: `(altitude m, temperature C, pressure hPa)`.
const STATIONS: [(f64, f64, f64); 4] = [
    (0.0, 15.0, 1013.25),
    (0.0, 35.0, 1013.25),
    (1_500.0, 25.0, 830.0),
    (300.0, -10.0, 980.0),
];
const HUMID: [f64; 4] = [12.5, 50.0, 90.0, 100.0];
/// Heights above (and one below) the station, m. All far enough from it that `ln(P / P_0)`
/// carries at least 13 significant digits.
const OFFSETS_M: [f64; 6] = [50.0, 1_000.0, 3_000.0, 9_000.0, 14_000.0, -200.0];

fn gravities() -> [Gravity; 2] {
    [Gravity::EARTH, Gravity::new(1.62).expect("lunar gravity")]
}

// =================================================================================================
// At 0 %: bit-identical to the model before the aloft change
// =================================================================================================

/// The lapse before the aloft change, verbatim: the ICAO layers in geopotential height with the dry
/// `R_AIR`.
fn layers_before_0_50_0(
    base_geopotential_m: f64,
    target_geopotential_m: f64,
    mut temp_k: f64,
    mut pressure_hpa: f64,
    gravity: Gravity,
) -> (f64, f64) {
    fn segment(
        base_temp_k: f64,
        base_pressure_hpa: f64,
        height_diff: f64,
        lapse_rate: f64,
        gravity: Gravity,
    ) -> (f64, f64) {
        let temp_k = base_temp_k + lapse_rate * height_diff;
        let pressure_hpa = if lapse_rate.abs() < 1e-10 {
            base_pressure_hpa * (-gravity.magnitude() * height_diff / (R_AIR * base_temp_k)).exp()
        } else {
            let temp_ratio = temp_k / base_temp_k;
            base_pressure_hpa * temp_ratio.powf(-gravity.magnitude() / (lapse_rate * R_AIR))
        };
        (temp_k, pressure_hpa)
    }
    let mut current_alt = base_geopotential_m;
    if target_geopotential_m > current_alt {
        while current_alt < target_geopotential_m {
            let layer_index = ICAO_LAYERS
                .iter()
                .rposition(|layer| current_alt >= layer.base_altitude)
                .unwrap_or(0);
            let segment_end = ICAO_LAYERS
                .get(layer_index + 1)
                .map_or(target_geopotential_m, |next| {
                    target_geopotential_m.min(next.base_altitude)
                });
            (temp_k, pressure_hpa) = segment(
                temp_k,
                pressure_hpa,
                segment_end - current_alt,
                ICAO_LAYERS[layer_index].lapse_rate,
                gravity,
            );
            current_alt = segment_end;
        }
    } else {
        while current_alt > target_geopotential_m {
            let layer_index = ICAO_LAYERS
                .iter()
                .rposition(|layer| current_alt > layer.base_altitude)
                .unwrap_or(0);
            let segment_end = if layer_index == 0 {
                target_geopotential_m
            } else {
                target_geopotential_m.max(ICAO_LAYERS[layer_index].base_altitude)
            };
            (temp_k, pressure_hpa) = segment(
                temp_k,
                pressure_hpa,
                segment_end - current_alt,
                ICAO_LAYERS[layer_index].lapse_rate,
                gravity,
            );
            current_alt = segment_end;
        }
    }
    (temp_k, pressure_hpa)
}

/// The per-step air before the aloft change, verbatim: the dry lapse, the density from the base
/// ratio, and `moist_speed_of_sound` at the LOCAL temperature and pressure and the base relative
/// humidity.
#[allow(clippy::too_many_arguments)]
fn air_before_0_50_0(
    altitude_m: f64,
    base_alt: f64,
    base_temp_c: f64,
    base_press_hpa: f64,
    base_ratio: f64,
    humidity_percent: f64,
    gravity: Gravity,
) -> LocalAir {
    let base_temp_k = base_temp_c + 273.15;
    let (temp_k, pressure_hpa) = layers_before_0_50_0(
        geometric_to_geopotential_height_m(base_alt),
        geometric_to_geopotential_height_m(altitude_m),
        base_temp_k,
        base_press_hpa,
        gravity,
    );
    let density_ratio = base_ratio * (base_temp_k * pressure_hpa) / (base_press_hpa * temp_k);
    let pressure_pa = pressure_hpa * 100.0;
    LocalAir {
        density_kg_m3: density_ratio * 1.225,
        speed_of_sound_mps: moist_speed_of_sound(temp_k, pressure_pa, humidity_percent),
        temperature_k: temp_k,
        pressure_pa,
    }
}

fn local_bits(air: LocalAir) -> [u64; 4] {
    [
        air.density_kg_m3.to_bits(),
        air.speed_of_sound_mps.to_bits(),
        air.temperature_k.to_bits(),
        air.pressure_pa.to_bits(),
    ]
}

#[test]
fn dry_air_aloft_is_bit_identical_to_the_model_before() {
    // Dry air over station temperatures from -30 to 40 C and pressures from 700 to 1050 hPa,
    // heights from 10 m below the station to 12 km above it, and Earth's and the Moon's gravity,
    // at two station altitudes, through the station and through a dry zone; plus the legacy
    // wrapper that shares the composition.
    let mut checked = 0usize;
    for base_alt in [0.0, 1_500.0] {
        for temp_c in [-30.0, 0.0, 15.0, 40.0] {
            for pressure_hpa in [700.0, 850.0, 1013.25, 1050.0] {
                let station =
                    StationAir::new(base_alt, temp_c, pressure_hpa, RelativeHumidity::DRY)
                        .expect("dry air exists");
                let atmosphere = KernelAtmosphere::Station(station);
                let zone = (temp_c + 4.0, pressure_hpa - 6.0, 0.0, 1e9);
                let zones = AtmoSock::new(vec![zone]);
                let zone_ratio = calculate_air_density_cimp(zone.0, zone.1, 0.0) / 1.225;
                for offset in [0.0, 10.0, -10.0, 500.0, 3_000.0, 12_000.0] {
                    let altitude_m = base_alt + offset;
                    for gravity in gravities() {
                        let before = air_before_0_50_0(
                            altitude_m,
                            base_alt,
                            temp_c,
                            pressure_hpa,
                            station.density_ratio(),
                            0.0,
                            gravity,
                        );
                        let what = format!(
                            "{base_alt} m {temp_c} C {pressure_hpa} hPa +{offset} m g {}",
                            gravity.magnitude()
                        );
                        assert_eq!(
                            local_bits(air_at_with_gravity(
                                &atmosphere,
                                None,
                                0.0,
                                altitude_m,
                                gravity
                            )),
                            local_bits(before),
                            "station {what}"
                        );
                        let (density, sound) = get_local_atmosphere_humid_with_gravity(
                            altitude_m,
                            base_alt,
                            temp_c,
                            pressure_hpa,
                            station.density_ratio(),
                            0.0,
                            gravity,
                        );
                        assert_eq!(
                            [density.to_bits(), sound.to_bits()],
                            [
                                before.density_kg_m3.to_bits(),
                                before.speed_of_sound_mps.to_bits()
                            ],
                            "legacy wrapper {what}"
                        );
                        let zone_before = air_before_0_50_0(
                            altitude_m, base_alt, zone.0, zone.1, zone_ratio, 0.0, gravity,
                        );
                        assert_eq!(
                            local_bits(air_at_with_gravity(
                                &atmosphere,
                                Some(&zones),
                                10.0,
                                altitude_m,
                                gravity
                            )),
                            local_bits(zone_before),
                            "dry zone {what}"
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    assert_eq!(checked, 2 * 4 * 4 * 6 * 2);
}

// =================================================================================================
// Above 0 %: one vapour mole fraction for the sound speed, the density and the lapse
// =================================================================================================

#[test]
fn the_sound_speed_aloft_holds_the_governing_airs_vapour_mole_fraction() {
    let mut largest_old_model_error = 0.0_f64;
    for (base_alt, temp_c, pressure_hpa) in STATIONS {
        for humidity in HUMID {
            let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity))
                .expect("possible air");
            let atmosphere = KernelAtmosphere::Station(station);
            let x_v = base_x_v(temp_c, pressure_hpa, humidity);
            assert!(x_v > 0.0);
            for gravity in gravities() {
                for offset in OFFSETS_M {
                    let altitude_m = base_alt + offset;
                    let air = air_at_with_gravity(&atmosphere, None, 0.0, altitude_m, gravity);
                    let expected = sound_speed_at(air.temperature_k, x_v);
                    assert!(
                        relative_error(air.speed_of_sound_mps, expected) <= 1e-14,
                        "{base_alt} m {temp_c} C {pressure_hpa} hPa {humidity}% +{offset} m g {}: \
                         c {} where the station x_v {x_v} gives {expected}",
                        gravity.magnitude(),
                        air.speed_of_sound_mps
                    );
                    // The model before the aloft change re-derived x_v at the local state, so it
                    // fails this check by far more than any tolerance.
                    let old = air_before_0_50_0(
                        altitude_m,
                        base_alt,
                        temp_c,
                        pressure_hpa,
                        station.density_ratio(),
                        humidity,
                        gravity,
                    );
                    largest_old_model_error = largest_old_model_error.max(relative_error(
                        old.speed_of_sound_mps,
                        sound_speed_at(old.temperature_k, x_v),
                    ));
                }
            }
        }
    }
    assert!(
        largest_old_model_error > 1e-3,
        "the separation companion must fail by more than 1e-3: {largest_old_model_error:e}"
    );
}

#[test]
fn the_density_and_the_sound_speed_aloft_carry_the_same_vapour_mole_fraction() {
    for (base_alt, temp_c, pressure_hpa) in STATIONS {
        for humidity in HUMID {
            let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity))
                .expect("possible air");
            let atmosphere = KernelAtmosphere::Station(station);
            let x_v = base_x_v(temp_c, pressure_hpa, humidity);
            // The station's own composition: P / (rho T), the CIPM density at the station x_v.
            let station_air = air_at(&atmosphere, None, 0.0, base_alt);
            let composition =
                station_air.pressure_pa / (station_air.density_kg_m3 * station_air.temperature_k);
            for offset in OFFSETS_M {
                let air = air_at(&atmosphere, None, 0.0, base_alt + offset);
                let what =
                    format!("{base_alt} m {temp_c} C {pressure_hpa} hPa {humidity}% +{offset} m");
                // The density aloft keeps the station's composition, so its x_v...
                assert!(
                    relative_error(
                        air.pressure_pa / (air.density_kg_m3 * air.temperature_k),
                        composition
                    ) <= 1e-13,
                    "density composition {what}"
                );
                // ...and the sound speed aloft carries that same x_v.
                let x_sound = x_v_from_sound_speed(air.speed_of_sound_mps, air.temperature_k);
                assert!(
                    relative_error(x_sound, x_v) <= 1e-9,
                    "sound speed x_v {x_sound} where the density carries {x_v}: {what}"
                );
            }
        }
    }
}

/// `ln(P / P_0)` of the air at `altitude_m` over a station at `(base_alt, T, P)`.
fn ln_pressure_ratio(
    base_alt: f64,
    temp_c: f64,
    pressure_hpa: f64,
    humidity: f64,
    altitude_m: f64,
    gravity: Gravity,
) -> f64 {
    let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity)).unwrap();
    let air = air_at_with_gravity(
        &KernelAtmosphere::Station(station),
        None,
        0.0,
        altitude_m,
        gravity,
    );
    (air.pressure_pa / (pressure_hpa * 100.0)).ln()
}

#[test]
fn the_lapse_scales_ln_p_by_the_moist_molar_mass_ratio() {
    // ln(P/P0) = -(g M_m / (R* ...)) ...: at the same temperature profile the moist exponent is
    // the dry one times M_m / M_a = R_AIR / R_m, on a lapse layer, on the isothermal layer, and
    // across a layer boundary.
    let cases = [
        // Troposphere (lapse layer).
        (0.0, 15.0, 1013.25, 50.0, 1_000.0),
        (0.0, 15.0, 1013.25, 50.0, 3_000.0),
        (0.0, 35.0, 1013.25, 90.0, 3_000.0),
        (1_500.0, 25.0, 830.0, 100.0, 1_450.0),
        (1_500.0, 25.0, 830.0, 100.0, 2_500.0),
        // The isothermal layer, 11-20 km geopotential: a warm, humid station inside it.
        (11_500.0, 20.0, 250.0, 60.0, 15_000.0),
        (11_500.0, 20.0, 250.0, 60.0, 11_900.0),
        // Across the tropopause.
        (0.0, 30.0, 1013.25, 80.0, 14_000.0),
    ];
    for (base_alt, temp_c, pressure_hpa, humidity, altitude_m) in cases {
        let x_v = base_x_v(temp_c, pressure_hpa, humidity);
        for gravity in gravities() {
            let moist = ln_pressure_ratio(
                base_alt,
                temp_c,
                pressure_hpa,
                humidity,
                altitude_m,
                gravity,
            );
            let dry = ln_pressure_ratio(base_alt, temp_c, pressure_hpa, 0.0, altitude_m, gravity);
            assert!(
                relative_error(moist / dry, molar_mass_ratio(x_v)) <= 1e-12,
                "{base_alt} m {temp_c} C {pressure_hpa} hPa {humidity}% to {altitude_m} m g {}: \
                 ln ratio {} where M_m/M_a = {}",
                gravity.magnitude(),
                moist / dry,
                molar_mass_ratio(x_v)
            );
        }
    }
}

#[test]
fn one_layer_of_the_lapse_scales_ln_p_by_r_air_over_r_m() {
    // The layer integrator itself, at the moist gas constant R_m = R_AIR / (M_m / M_a): the
    // temperature is the dry one to the bit, and ln(P / P0) is the dry value times R_AIR / R_m,
    // on a lapse layer (up and down), on an isothermal layer and on a warming layer.
    assert_eq!(
        (R_AIR / moist_molar_mass_ratio(0.0)).to_bits(),
        R_AIR.to_bits()
    );
    for x_v in [0.0021, 0.00845, 0.0502, 0.2] {
        let gas_constant = R_AIR / moist_molar_mass_ratio(x_v);
        assert!(
            gas_constant > R_AIR,
            "moist air is lighter: a larger gas constant"
        );
        for gravity in gravities() {
            for (lapse_rate, height_diff) in [
                (-0.0065, 1_000.0),
                (-0.0065, -300.0),
                (0.0, 2_500.0),
                (0.0, -800.0),
                (0.001, 4_000.0),
            ] {
                let (t_moist, p_moist) = integrate_local_atmosphere_segment(
                    288.15,
                    1013.25,
                    height_diff,
                    lapse_rate,
                    gas_constant,
                    gravity,
                );
                let (t_dry, p_dry) = integrate_local_atmosphere_segment(
                    288.15,
                    1013.25,
                    height_diff,
                    lapse_rate,
                    R_AIR,
                    gravity,
                );
                assert_eq!(t_moist.to_bits(), t_dry.to_bits());
                let ratio = (p_moist / 1013.25).ln() / (p_dry / 1013.25).ln();
                assert!(
                    relative_error(ratio, R_AIR / gas_constant) <= 1e-12,
                    "x_v {x_v} lapse {lapse_rate} dh {height_diff} g {}: {ratio} vs {}",
                    gravity.magnitude(),
                    R_AIR / gas_constant
                );
            }
        }
    }
}

#[test]
fn a_zone_governs_the_vapour_aloft_with_its_own_mole_fraction() {
    // The zone's air is the base air_at lapses from (its density ratio already was), so its own
    // x_v is the one held aloft, for the sound speed and the lapse.
    let station = StationAir::new(200.0, 10.0, 990.0, rh(20.0)).unwrap();
    let atmosphere = KernelAtmosphere::Station(station);
    let zones = AtmoSock::new(vec![(30.0, 985.0, 85.0, 400.0), (5.0, 995.0, 0.0, 1e9)]);
    let dry_zone_twin = AtmoSock::new(vec![(30.0, 985.0, 0.0, 400.0), (5.0, 995.0, 0.0, 1e9)]);
    let x_zone = base_x_v(30.0, 985.0, 85.0);
    for offset in [150.0, 1_000.0, 2_500.0] {
        let altitude_m = 200.0 + offset;
        let air = air_at(&atmosphere, Some(&zones), 100.0, altitude_m);
        let expected = sound_speed_at(air.temperature_k, x_zone);
        assert!(
            relative_error(air.speed_of_sound_mps, expected) <= 1e-14,
            "zone +{offset} m: c {} where the zone x_v gives {expected}",
            air.speed_of_sound_mps
        );
        let dry = air_at(&atmosphere, Some(&dry_zone_twin), 100.0, altitude_m);
        let ln_moist = (air.pressure_pa / 98_500.0).ln();
        let ln_dry = (dry.pressure_pa / 98_500.0).ln();
        assert!(
            relative_error(ln_moist / ln_dry, molar_mass_ratio(x_zone)) <= 1e-12,
            "zone +{offset} m: ln ratio {}",
            ln_moist / ln_dry
        );
        // Past the first zone the dry zone governs: dry air, whatever the station's humidity.
        let beyond = air_at(&atmosphere, Some(&zones), 600.0, altitude_m);
        assert_eq!(
            beyond.speed_of_sound_mps.to_bits(),
            (GAMMA * R_AIR * beyond.temperature_k).sqrt().to_bits()
        );
    }
}

#[test]
fn at_the_station_height_the_air_is_the_station_air() {
    // The unified model changes nothing at the governing air's own height.
    for (base_alt, temp_c, pressure_hpa) in STATIONS {
        for humidity in [0.0, 12.5, 50.0, 90.0, 100.0] {
            let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity)).unwrap();
            for gravity in gravities() {
                let air = air_at_with_gravity(
                    &KernelAtmosphere::Station(station),
                    None,
                    0.0,
                    base_alt,
                    gravity,
                );
                assert_eq!(
                    air.speed_of_sound_mps.to_bits(),
                    station.speed_of_sound_mps().to_bits()
                );
                assert_eq!(air.temperature_k.to_bits(), (temp_c + 273.15).to_bits());
                assert_eq!(air.pressure_pa.to_bits(), (pressure_hpa * 100.0).to_bits());
            }
        }
    }
}

// =================================================================================================
// An independent hand calculation
// =================================================================================================

/// `air_at` against a hand calculation from the published formulas, written independently of
/// the engine's code and evaluated in Python doubles: the station x_v from CIPM-2007
/// (enhancement factor, Wagner-Pruss saturation pressure), T from the ICAO lapse rate in
/// geopotential height, P from the virtual-temperature lapse `P0 (T/T0)^(-g M_m / (L R_d M_a))`
/// (or `P0 exp(-g dh M_m / (R_d M_a T0))` on the isothermal layer), the density as the station's
/// CIPM density times `(T0 P) / (P0 T)`, and `c = sqrt(1.4 (1 - 0.062 x_v) 287.0531 (1 + 0.378
/// x_v) T)` at the station x_v.
#[test]
fn the_moist_lapse_matches_a_hand_calculation() {
    // (station alt m, T C, P hPa, RH %, query alt m) -> (x_v, T K, P Pa, rho kg/m3, c m/s)
    #[allow(clippy::type_complexity)]
    let cases: [((f64, f64, f64, f64, f64), (f64, f64, f64, f64, f64)); 5] = [
        (
            (0.0, 15.0, 1013.25, 50.0, 1000.0),
            (
                0.008449919693477451,
                281.6510223716947,
                89910.71626853914,
                1.1090273440008176,
                336.8833060538346,
            ),
        ),
        (
            (0.0, 15.0, 1013.25, 50.0, 3000.0),
            (
                0.008449919693477451,
                268.65919845164115,
                70203.66851104423,
                0.9078209920327416,
                329.0218059249786,
            ),
        ),
        (
            (0.0, 35.0, 1013.25, 90.0, 1000.0),
            (
                0.050223405365983464,
                301.6510223716947,
                90778.88172414222,
                1.028908684360444,
                350.9167495141029,
            ),
        ),
        (
            (0.0, 35.0, 1013.25, 90.0, 3000.0),
            (
                0.050223405365983464,
                288.65919845164115,
                72344.17798562118,
                0.8568701044639647,
                343.2767562867101,
            ),
        ),
        (
            (11500.0, 20.0, 250.0, 60.0, 15000.0),
            (
                0.05623210279727311,
                293.15,
                16799.04055999137,
                0.19542304087937543,
                346.2573166522401,
            ),
        ),
    ];
    for ((base_alt, temp_c, pressure_hpa, humidity, altitude_m), (x_v, t_k, p_pa, rho, c)) in cases
    {
        let what =
            format!("{base_alt} m {temp_c} C {pressure_hpa} hPa {humidity}% to {altitude_m} m");
        assert!(
            relative_error(base_x_v(temp_c, pressure_hpa, humidity), x_v) <= 1e-12,
            "x_v {what}"
        );
        let station = StationAir::new(base_alt, temp_c, pressure_hpa, rh(humidity)).unwrap();
        let air = air_at(&KernelAtmosphere::Station(station), None, 0.0, altitude_m);
        for (name, got, want) in [
            ("T", air.temperature_k, t_k),
            ("P", air.pressure_pa, p_pa),
            ("rho", air.density_kg_m3, rho),
            ("c", air.speed_of_sound_mps, c),
        ] {
            assert!(
                relative_error(got, want) <= 1e-12,
                "{name} {what}: engine {got:?}, hand calculation {want:?}"
            );
        }
    }
}
