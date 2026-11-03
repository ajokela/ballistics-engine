//! The governing air is computed once per base air, never per step.
//!
//! A flight's governing air is the station's, or inside a downrange [`AtmoSock`] zone the zone's
//! own. Its water-vapour mole fraction `x_v` (held at every height since the aloft change) and, for
//! a zone, its CIPM-2007 density ratio are fixed for the whole flight. They are computed when the
//! station or the zones are built and read at every step, so [`air_at`] evaluates no Wagner-Pruss
//! saturation pressure.
//!
//! - **No per-step work.** [`air_at`] over a built station, through either constructor, with no
//!   zones, an empty sock, valid zones or unvalidated ones, evaluates no saturation pressure at
//!   any point, at every humidity, dry included
//!   ([`air_at_evaluates_no_saturation_pressure_at_any_point`]).
//! - **The same bits.** The held values give exactly the air that a composition computing
//!   them afresh from the raw base gives, and the legacy wrapper's too, which still computes
//!   them on every call ([`the_held_air_is_the_air_a_fresh_composition_computes`]).
//! - **The same Debug.** The held values are a cache, not state: [`StationAir`]'s and
//!   [`AtmoSock`]'s `Debug` output names the fields it named before they held it
//!   ([`the_held_values_leave_debug_output_unchanged`]).

use super::*;

fn rh(percent: f64) -> RelativeHumidity {
    RelativeHumidity::try_percent(percent).expect("valid test humidity")
}

/// Stations through both constructors: [`StationAir::new`] over a grid of altitudes,
/// temperatures, pressures and humidities (dry, partial, the default and saturated), and
/// [`StationAir::with_density_ratio`], whose density ratio is the caller's.
fn stations() -> Vec<StationAir> {
    let mut out = Vec::new();
    for altitude_m in [0.0, 1_500.0] {
        for temperature_c in [-30.0, 15.0, 35.0] {
            for pressure_hpa in [700.0, 1013.25] {
                for humidity_percent in [0.0, 37.5, 50.0, 100.0] {
                    out.push(
                        StationAir::new(
                            altitude_m,
                            temperature_c,
                            pressure_hpa,
                            rh(humidity_percent),
                        )
                        .expect("air that exists"),
                    );
                    out.push(
                        StationAir::with_density_ratio(
                            altitude_m,
                            temperature_c,
                            pressure_hpa,
                            0.93,
                            rh(humidity_percent),
                        )
                        .expect("air that exists"),
                    );
                }
            }
        }
    }
    out
}

/// No zones, an empty sock, valid zones at several humidities, and zones that
/// [`AtmoSock::validate`] refuses or would snap (a NaN, 150 % and -3 % humidity, and a NaN
/// threshold), which the per-step lookup still composes without validating.
fn zone_sets() -> Vec<Option<AtmoSock>> {
    vec![
        None,
        Some(AtmoSock::new(Vec::new())),
        Some(AtmoSock::new(vec![
            (17.0, 1010.0, 20.0, 300.0),
            (12.0, 1015.0, 80.0, 900.0),
            (5.0, 990.0, 0.0, f64::INFINITY),
        ])),
        Some(AtmoSock::new(vec![(30.0, 985.0, 100.0, 400.0)])),
        Some(AtmoSock::new(vec![
            (15.0, 1013.25, f64::NAN, 200.0),
            (15.0, 1013.25, 150.0, 400.0),
            (15.0, 1013.25, -3.0, 600.0),
            (15.0, 1013.25, 50.0, f64::NAN),
        ])),
    ]
}

/// Downrange distances, zone thresholds and NaN included.
const DOWNRANGE_M: [f64; 8] = [-1.0, 0.0, 199.9, 200.0, 300.0, 650.0, 1e4, f64::NAN];

/// Heights above the station, m.
const HEIGHTS_M: [f64; 6] = [0.0, 10.0, -10.0, 500.0, 3_000.0, 12_000.0];

fn gravities() -> [Gravity; 2] {
    [Gravity::EARTH, Gravity::new(1.62).expect("lunar gravity")]
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
fn air_at_evaluates_no_saturation_pressure_at_any_point() {
    // Everything is built first: the station's held x_v, and each zone's held x_v and density
    // ratio, are the construction's work, counted outside the lookups.
    let stations = stations();
    let zone_sets = zone_sets();
    let before = air_trace::saturation_pressures();
    let mut lookups = 0u64;
    for station in &stations {
        let atmosphere = KernelAtmosphere::Station(*station);
        for zones in &zone_sets {
            for downrange_m in DOWNRANGE_M {
                for height_m in HEIGHTS_M {
                    for gravity in gravities() {
                        std::hint::black_box(air_at_with_gravity(
                            &atmosphere,
                            zones.as_ref(),
                            downrange_m,
                            station.altitude_m() + height_m,
                            gravity,
                        ));
                        lookups += 1;
                    }
                }
            }
        }
    }
    let evaluated = air_trace::saturation_pressures() - before;
    assert_eq!(
        evaluated, 0,
        "{lookups} per-step air lookups evaluated {evaluated} Wagner-Pruss saturation \
         pressures; the governing air's x_v and a zone's density ratio are fixed for the flight"
    );
}

#[test]
fn the_held_air_is_the_air_a_fresh_composition_computes() {
    // The reference recomputes every governing value from the raw base on every call: the zone
    // the lookup selects (atmo_for_range, the public lookup), the zone's CIPM density ratio, and
    // the governing x_v through the one formula. The held values must give the same air to the
    // bit, at every point the per-step lookup can be asked.
    let mut checked = 0usize;
    for station in stations() {
        let atmosphere = KernelAtmosphere::Station(station);
        for zones in zone_sets() {
            for downrange_m in DOWNRANGE_M {
                let (base_temp_c, base_press_hpa, base_ratio, base_humidity_percent) =
                    match zones.as_ref().filter(|sock| !sock.is_empty()) {
                        Some(sock) => {
                            let (t, p, h) = sock.atmo_for_range(downrange_m);
                            (t, p, calculate_air_density_cimp(t, p, h) / 1.225, h)
                        }
                        None => (
                            station.temperature_c(),
                            station.pressure_hpa(),
                            station.density_ratio(),
                            station.humidity().percent(),
                        ),
                    };
                let x_v = governing_vapor_mole_fraction(
                    base_temp_c,
                    base_press_hpa,
                    base_humidity_percent,
                );
                for height_m in HEIGHTS_M {
                    let altitude_m = station.altitude_m() + height_m;
                    for gravity in gravities() {
                        let (temperature_k, pressure_pa, density_kg_m3) =
                            local_temp_pressure_density_with_gravity(
                                altitude_m,
                                station.altitude_m(),
                                base_temp_c,
                                base_press_hpa,
                                base_ratio,
                                x_v,
                                gravity,
                            );
                        let fresh = LocalAir {
                            density_kg_m3,
                            speed_of_sound_mps: moist_sound_speed_at_mole_fraction(
                                temperature_k,
                                x_v,
                            ),
                            temperature_k,
                            pressure_pa,
                        };
                        let held = air_at_with_gravity(
                            &atmosphere,
                            zones.as_ref(),
                            downrange_m,
                            altitude_m,
                            gravity,
                        );
                        let what = format!(
                            "{station:?} zones {zones:?} at {downrange_m} m downrange, \
                             +{height_m} m, g {}",
                            gravity.magnitude()
                        );
                        assert_eq!(local_bits(held), local_bits(fresh), "{what}");
                        let (legacy_density, legacy_sound) =
                            get_local_atmosphere_humid_with_gravity(
                                altitude_m,
                                station.altitude_m(),
                                base_temp_c,
                                base_press_hpa,
                                base_ratio,
                                base_humidity_percent,
                                gravity,
                            );
                        assert_eq!(
                            [legacy_density.to_bits(), legacy_sound.to_bits()],
                            [
                                held.density_kg_m3.to_bits(),
                                held.speed_of_sound_mps.to_bits()
                            ],
                            "legacy wrapper: {what}"
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    assert_eq!(
        checked,
        stations().len() * zone_sets().len() * DOWNRANGE_M.len() * HEIGHTS_M.len() * 2
    );
}

/// The types' `Debug` before they held the vapour, derived from exactly the fields they had then.
mod before_the_held_vapour {
    #![allow(dead_code)] // read only through `Debug`

    use super::{AtmoSegment, RelativeHumidity};

    #[derive(Debug)]
    pub(super) struct StationAir {
        pub altitude_m: f64,
        pub temperature_c: f64,
        pub pressure_hpa: f64,
        pub humidity: RelativeHumidity,
        pub density_ratio: f64,
        pub density_kg_m3: f64,
        pub speed_of_sound_mps: f64,
    }

    #[derive(Debug)]
    pub(super) struct AtmoSock {
        pub segments: Vec<AtmoSegment>,
    }

    #[derive(Debug)]
    pub(super) enum KernelAtmosphere {
        Station(StationAir),
    }
}

#[test]
fn the_held_values_leave_debug_output_unchanged() {
    for station in stations() {
        let before = before_the_held_vapour::StationAir {
            altitude_m: station.altitude_m(),
            temperature_c: station.temperature_c(),
            pressure_hpa: station.pressure_hpa(),
            humidity: station.humidity(),
            density_ratio: station.density_ratio(),
            density_kg_m3: station.density_kg_m3(),
            speed_of_sound_mps: station.speed_of_sound_mps(),
        };
        assert_eq!(format!("{station:?}"), format!("{before:?}"));
        assert_eq!(format!("{station:#?}"), format!("{before:#?}"));
        let atmosphere = KernelAtmosphere::Station(station);
        let atmosphere_before = before_the_held_vapour::KernelAtmosphere::Station(before);
        assert_eq!(format!("{atmosphere:?}"), format!("{atmosphere_before:?}"));
        assert_eq!(
            format!("{atmosphere:#?}"),
            format!("{atmosphere_before:#?}")
        );
    }
    for zones in zone_sets().into_iter().flatten() {
        let before = before_the_held_vapour::AtmoSock {
            segments: zones.segments.clone(),
        };
        assert_eq!(format!("{zones:?}"), format!("{before:?}"));
        assert_eq!(format!("{zones:#?}"), format!("{before:#?}"));
    }
}
