//! Describe an atmosphere the way a solve resolves it: the service behind the
//! `atmosphere.density_altitude` bridge command.
//!
//! # Why this exists
//!
//! Density altitude was reachable only as one field inside a DOPE card, so an embedding app
//! wanting to show it had to reimplement the formula — and an external Android app did,
//! then hand-maintained the copy with its own regression tests (2026-09). A copy of a formula
//! is a copy that drifts. This service answers the question directly, from the same code the
//! solver and the card use.
//!
//! # One atmosphere schema, decoded and resolved one way
//!
//! The request carries the SAME `atmosphere` object a solve-json v1 request does
//! ([`AtmosphereV1`]). It is decoded by the solve's own shape validator and resolved by the
//! solve's own resolver (`solve_v1::resolve_atmosphere`), so every error the solve would give
//! for that object comes back identically — same code, same `$.atmosphere.*` path — as do the
//! ICAO-standard defaults for whatever is omitted (each announced in `assumptions`) and the
//! QNH-to-station reduction. An app can send the atmosphere it is about to solve with and get
//! the density altitude of precisely that air.
//!
//! # Then checks the solve does not make (yet — MBA-1594)
//!
//! The solve's resolver only requires temperature and pressure to be positive, which admits air
//! that cannot exist: 5 K, a pressure altitude of 30 km. (A humidity whose vapor pressure
//! exceeds the air's total pressure the solve refuses too: both ask the one rule,
//! `StationAir::new`, through `solve_v1::resolved_station_air`, which blames such air in the
//! order of the checks below, so whatever the solve refuses this service refuses with the same
//! code and path, and with the same message unless one of these checks refused it first.)
//! Those are also exactly what the two likeliest
//! unit mistakes produce — °C sent as `temperature_k`, hPa sent as `pressure_pa` — so a
//! figure computed from them is a wrong answer delivered confidently. This service refuses
//! them instead, with the same `invalid_value` envelope and the path of the input to blame:
//!
//! - `temperature_k` outside 173.15–373.15 K (−100 °C to 100 °C), where the engine's
//!   saturation-vapor and compressibility formulas are defined;
//! - station pressure whose pressure altitude is outside −5 km to 11 km — the standard
//!   atmosphere's troposphere, which is all the closed-form arithmetic here covers;
//! - a relative humidity that would need more vapor pressure than the total pressure
//!   (`StationAir::new`'s refusal, the solve's own);
//! - either density altitude outside −5 km to 11 km.
//!
//! None of these rejects a condition anyone has fired a rifle in.
//!
//! # Two density altitudes, by name
//!
//! "Density altitude" means two different things in the field, and one number under one name
//! is how the confusion this module answers arose. The response therefore never offers a bare
//! `density_altitude`; it offers both, and a caller has to choose:
//!
//! - **`faa_rule`** — [`crate::atmosphere::faa_rule_density_altitude_ft`]: NWS pressure
//!   altitude plus the FAA's "120 ft per °C" rule of thumb. Humidity-free. This is the formula
//!   the engine's DOPE card header uses, and the exact inverse of the density-altitude ENTRY
//!   mode.
//! - **`density_matched`** — [`crate::atmosphere::density_matched_altitude_m`]: the ISA
//!   altitude whose standard density equals this air's actual density, humidity included. The
//!   textbook definition, and what humidity-aware tools report — the National Weather
//!   Service's own density-altitude calculator included.
//!
//! They coincide at ICAO standard sea-level dry air; anywhere else a match is coincidence.
//! They differ by 33.1 m at 15 °C, 1013.25 hPa and 50% RH (humidity), and by about 22 m at
//! 30 °C in perfectly DRY air (the FAA rule's straight-line approximation). The difference
//! between them is therefore NOT a humidity correction, and must not be presented as one.
//!
//! `air_density_kg_m3` is reported as well: the CIPM-2007 density the solver itself uses at
//! the muzzle, which answers "how thin is the air" with no definition to argue about.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::atmosphere::RelativeHumidity;
use crate::solve_json::{
    AtmosphereV1, PressureReferenceV1, ResolvedAtmosphereV1, SolveErrorCodeV1,
    SolveErrorEnvelopeV1, SolveErrorV1, SolveNoticeV1,
};
use crate::solve_v1::{KELVIN_OFFSET_C, PASCALS_PER_HECTOPASCAL};

const METERS_PER_FOOT: f64 = 0.3048;
/// The wire's one humidity field: a FRACTION, as in solve-json v1.
const RELATIVE_HUMIDITY_PATH: &str = "$.atmosphere.relative_humidity";

/// −100 °C. The engine's saturation-vapor-pressure and CIPM compressibility formulas clamp
/// their temperature here, so below it the density would mix a clamped and an unclamped term.
/// The one bound every surface's impossible-air refusal blames the temperature beyond.
const MIN_TEMPERATURE_K: f64 = crate::atmosphere::HUMID_AIR_MODEL_MIN_K;
/// 100 °C. Above the boiling point the saturation-vapor-pressure formula stops describing a
/// stable atmosphere, and above water's critical point (647 K) it is not defined at all.
const MAX_TEMPERATURE_K: f64 = crate::atmosphere::HUMID_AIR_MODEL_MAX_K;
/// The lowest altitude solve-json v1 accepts, and the bottom of the ICAO table.
const LOWEST_ALTITUDE_M: f64 = crate::atmosphere::LOWEST_PRESSURE_ALTITUDE_M;
/// The top of the standard atmosphere's troposphere: every closed-form altitude here is the
/// troposphere formula, and it keeps answering, wrongly, above this. The same band the shared
/// impossible-air refusal blames a pressure outside of.
const TROPOSPHERE_TOP_M: f64 = crate::atmosphere::TROPOSPHERE_TOP_M;

/// `atmosphere.density_altitude` request.
///
/// `atmosphere` is exactly the solve-json v1 `atmosphere` object: SI units
/// (`altitude_m`, `temperature_k`, `pressure_pa`), `relative_humidity` as a FRACTION 0..1,
/// and an optional `pressure_reference` of `absolute` (the default) or `qnh`. Every field is
/// optional; an omitted one resolves exactly as a solve would resolve it. Decode with
/// [`decode_density_altitude_request_v1`] to get the solve's error envelope for a malformed
/// object.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityAltitudeRequestV1 {
    pub atmosphere: AtmosphereV1,
}

/// One altitude in both units a shooter meets it in. Density altitude is conventionally
/// quoted in feet even by shooters who range in meters, so both travel together rather than
/// asking every caller to convert.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AltitudeV1 {
    pub m: f64,
    pub ft: f64,
}

impl AltitudeV1 {
    fn from_ft(ft: f64) -> Self {
        Self {
            m: ft * METERS_PER_FOOT,
            ft,
        }
    }

    fn from_m(m: f64) -> Self {
        Self {
            m,
            ft: m / METERS_PER_FOOT,
        }
    }
}

/// The two density altitudes, each under its own name. See the module docs for why there is
/// no single `density_altitude` value.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityAltitudesV1 {
    /// NWS pressure altitude plus the FAA 120 ft/°C rule: humidity-free, the formula the DOPE
    /// card header uses, and the exact inverse of density-altitude entry.
    pub faa_rule: AltitudeV1,
    /// ISA altitude of equal actual air density, humidity included.
    pub density_matched: AltitudeV1,
}

/// `atmosphere.density_altitude` result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DensityAltitudeResponseV1 {
    /// The atmosphere these figures describe, AFTER resolution: `pressure_pa` is station
    /// pressure (already reduced if a QNH was supplied), and every omitted field carries the
    /// value it was defaulted to. `pressure_reference` echoes the mode that was SENT, exactly
    /// as the solve's `resolved_request` does — so this is a report, not a request: re-sent
    /// with `qnh`, the pressure would be reduced twice.
    pub atmosphere: ResolvedAtmosphereV1,
    /// NWS pressure altitude of the resolved station pressure.
    pub pressure_altitude: AltitudeV1,
    pub density_altitude: DensityAltitudesV1,
    /// CIPM-2007 humid-air density, kg/m³ — the density the solver uses at the muzzle.
    pub air_density_kg_m3: f64,
    /// Every default applied and any QNH-to-station reduction made while resolving the
    /// atmosphere, in the solve's own notice format.
    pub assumptions: Vec<SolveNoticeV1>,
}

/// Decode an `atmosphere.density_altitude` request from JSON the way `solve` decodes its
/// `atmosphere`: the solve's own shape validator runs first, so an unknown or unit-suffixed
/// field, an explicit `null`, a string where a number belongs or an unrecognized
/// `pressure_reference` is refused with the solve's code at the solve's `$.atmosphere.*` path.
pub fn decode_density_altitude_request_v1(
    value: &Value,
) -> Result<DensityAltitudeRequestV1, SolveErrorEnvelopeV1> {
    let root = crate::solve_json::require_object(value, "$")?;
    crate::solve_json::validate_members(root, "$", &["atmosphere"], &["atmosphere"])?;
    crate::solve_json::validate_atmosphere(&root["atmosphere"])?;
    serde_json::from_value(value.clone()).map_err(|error| invalid_value("$", error.to_string()))
}

/// Resolve `request.atmosphere` exactly as a solve would, check that the air it describes can
/// exist and lies where these formulas hold, then report its pressure altitude, both density
/// altitudes and its air density.
///
/// Errors are the solve's own validation envelope, located at `$.atmosphere.*`.
pub fn density_altitude_v1(
    request: &DensityAltitudeRequestV1,
) -> Result<DensityAltitudeResponseV1, SolveErrorEnvelopeV1> {
    let mut assumptions = Vec::new();
    let atmosphere = crate::solve_v1::resolve_atmosphere(&request.atmosphere, &mut assumptions)?;

    // Temperature first: above water's critical point the vapor formula is NaN, which the
    // vapor check below would otherwise read as "fine". A defaulted temperature is ICAO at an
    // altitude the solve accepts, always inside this range, so only a supplied one can fail.
    if crate::atmosphere::outside_humid_air_model_range(atmosphere.temperature_k) {
        return Err(invalid_value(
            "$.atmosphere.temperature_k",
            format!(
                "{} K is outside {MIN_TEMPERATURE_K}–{MAX_TEMPERATURE_K} K (−100 °C to 100 °C), \
                 the range the humid-air model behind these figures is defined for; \
                 temperature_k is kelvin (°C + 273.15)",
                echo(atmosphere.temperature_k)
            ),
        ));
    }

    // The same conversions the solve applies before handing the atmosphere to the engine
    // (solve_v1.rs): Kelvin to Celsius, Pa to hPa, and the wire's humidity fraction through the
    // one converter, `RelativeHumidity::try_fraction`. The resolver has already read a sensor
    // reading within the band as the end it reads as, and refused anything beyond it, so this
    // cannot fail; it is mapped rather than unwrapped all the same.
    let temperature_c = atmosphere.temperature_k - KELVIN_OFFSET_C;
    let temperature_f = temperature_c * 9.0 / 5.0 + 32.0;
    let station_pressure_hpa = atmosphere.pressure_pa / PASCALS_PER_HECTOPASCAL;
    let humidity = RelativeHumidity::try_fraction(atmosphere.relative_humidity)
        .map_err(|error| invalid_value(RELATIVE_HUMIDITY_PATH, error.to_string()))?;

    // Pressure before vapor: hPa sent as Pa is the likelier mistake, and at 1013 Pa a warm
    // humid day would otherwise be blamed on the humidity.
    let pressure_altitude_ft = crate::atmosphere::nws_pressure_altitude_ft(station_pressure_hpa);
    let pressure_altitude = AltitudeV1::from_ft(pressure_altitude_ft);
    if !within_troposphere(pressure_altitude.m) {
        return Err(pressure_outside_troposphere(
            &request.atmosphere,
            &atmosphere,
            pressure_altitude.m,
        ));
    }

    // Whether this air can exist is the one rule's to decide, `StationAir::new`'s, which every
    // solve applies to the same air through this same function: a humidity needing more
    // water-vapour pressure than the air's total pressure is refused by both, with one code,
    // path and message. Air the checks above refused first is refused by the solve too, at the
    // same path: that function blames impossible air in their order (the temperature outside
    // this range, then the reading that took the pressure out of the troposphere). The station
    // it builds is the solver's muzzle air, so its density is the one reported below.
    let station = crate::solve_v1::resolved_station_air(&request.atmosphere, &atmosphere)?;
    debug_assert_eq!(station.humidity(), humidity);

    let faa_rule = AltitudeV1::from_ft(crate::atmosphere::faa_rule_density_altitude_ft(
        station_pressure_hpa,
        temperature_f,
    ));
    let density_matched = AltitudeV1::from_m(crate::atmosphere::density_matched_altitude_m(
        temperature_c,
        station_pressure_hpa,
        humidity.percent(),
    ));
    let air_density_kg_m3 = station.density_kg_m3();

    for (name, altitude) in [("faa_rule", faa_rule), ("density_matched", density_matched)] {
        if !within_troposphere(altitude.m) {
            return Err(invalid_value(
                "$.atmosphere",
                format!(
                    "the {name} density altitude of these conditions is {:.0} m; density \
                     altitude is computed only between {LOWEST_ALTITUDE_M} m and \
                     {TROPOSPHERE_TOP_M} m, the standard atmosphere's troposphere. Check that \
                     temperature_k is kelvin and pressure_pa is pascals",
                    altitude.m
                ),
            ));
        }
    }
    // Every input above is bounded, so this cannot fire; it is here so that if a bound is ever
    // loosened, the failure is a refusal rather than `ok: true` with a `null` (serde_json
    // writes NaN as null) or a negative-zero density.
    let figures = [
        pressure_altitude.m,
        pressure_altitude.ft,
        faa_rule.m,
        faa_rule.ft,
        density_matched.m,
        density_matched.ft,
    ];
    if !(figures.iter().all(|x| x.is_finite())
        && air_density_kg_m3.is_finite()
        && air_density_kg_m3 > 0.0)
    {
        return Err(invalid_value(
            "$.atmosphere",
            "these conditions do not produce a finite, positive air density",
        ));
    }

    Ok(DensityAltitudeResponseV1 {
        atmosphere,
        pressure_altitude,
        density_altitude: DensityAltitudesV1 {
            faa_rule,
            density_matched,
        },
        air_density_kg_m3,
        assumptions,
    })
}

/// The refusal for a station pressure outside the troposphere, blaming the input that is
/// actually wrong and quoting what the caller SENT rather than a pressure derived from it.
fn pressure_outside_troposphere(
    sent: &AtmosphereV1,
    resolved: &ResolvedAtmosphereV1,
    pressure_altitude_m: f64,
) -> SolveErrorEnvelopeV1 {
    let scope = format!(
        "is a pressure altitude of {pressure_altitude_m:.0} m; density altitude is computed only \
         between {LOWEST_ALTITUDE_M} m and {TROPOSPHERE_TOP_M} m, the standard atmosphere's \
         troposphere"
    );
    let station = echo(resolved.pressure_pa);
    let altitude = echo(resolved.altitude_m);
    // The reading to blame is the one the solve's impossible-air refusal blames for the same
    // station, by the one rule: a QNH that would itself be a troposphere pressure at sea level
    // is a plausible reading, so the altitude it was reduced to is what took it out of range
    // (feet sent as meters, say); an omitted pressure is the altitude's too.
    let path = match crate::solve_v1::station_pressure_origin(sent) {
        crate::atmosphere::StationPressureOrigin::Reading => "$.atmosphere.pressure_pa",
        crate::atmosphere::StationPressureOrigin::Standard
        | crate::atmosphere::StationPressureOrigin::Qnh => "$.atmosphere.altitude_m",
    };
    match sent.pressure_pa {
        Some(qnh) if sent.pressure_reference == Some(PressureReferenceV1::Qnh) => invalid_value(
            path,
            format!(
                "a QNH of {} Pa, reduced to a station pressure of {station} Pa at \
                     {altitude} m, {scope}. pressure_pa is pascals (hPa × 100) and altitude_m \
                     is meters",
                echo(qnh)
            ),
        ),
        Some(_) => invalid_value(
            path,
            format!(
                "a station pressure of {station} Pa {scope}. pressure_pa is pascals (hPa × 100)"
            ),
        ),
        None => invalid_value(
            path,
            format!(
                "the ICAO standard pressure at {altitude} m, {station} Pa, {scope}. altitude_m \
                 is meters"
            ),
        ),
    }
}

/// A caller's number as it should appear in a message: plain for ordinary magnitudes, and in
/// exponent form for the absurd ones, which `{}` would otherwise print as 300 digits.
fn echo(value: f64) -> String {
    crate::solve_v1::echo_si(value)
}

fn within_troposphere(altitude_m: f64) -> bool {
    (LOWEST_ALTITUDE_M..=TROPOSPHERE_TOP_M).contains(&altitude_m)
}

fn invalid_value(path: &str, message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(
        SolveErrorV1::new(SolveErrorCodeV1::InvalidValue, message).at_path(path),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn describe(atmosphere: AtmosphereV1) -> DensityAltitudeResponseV1 {
        density_altitude_v1(&DensityAltitudeRequestV1 { atmosphere }).expect("valid atmosphere")
    }

    fn refusal(atmosphere: AtmosphereV1) -> SolveErrorV1 {
        density_altitude_v1(&DensityAltitudeRequestV1 { atmosphere })
            .expect_err("atmosphere should be refused")
            .error
    }

    fn at(temperature_c: f64, pressure_hpa: f64, relative_humidity: f64) -> AtmosphereV1 {
        AtmosphereV1 {
            altitude_m: Some(0.0),
            temperature_k: Some(temperature_c + KELVIN_OFFSET_C),
            pressure_pa: Some(pressure_hpa * PASCALS_PER_HECTOPASCAL),
            relative_humidity: Some(relative_humidity),
            ..AtmosphereV1::default()
        }
    }

    #[test]
    fn icao_standard_dry_air_is_zero_by_both_definitions() {
        let r = describe(at(15.0, 1013.25, 0.0));
        assert!(r.density_altitude.faa_rule.m.abs() < 1e-9, "{r:?}");
        assert!(r.density_altitude.density_matched.m.abs() < 1e-9, "{r:?}");
        assert!(r.pressure_altitude.m.abs() < 1e-9, "{r:?}");
    }

    /// The reported case: 0 m by the FAA rule where a humidity-aware tool read 36 m.
    /// Density-matched gives 33.1 m; the remaining few meters are that tool's own method.
    ///
    /// 33.10 m was recomputed independently of this crate, from the published CIPM-2007
    /// coefficients and the ISA troposphere, and agreed to 0.01 m (2026-09-27) — so the pin is
    /// tight on purpose: a drift here is a change to the solver's air density.
    #[test]
    fn humidity_moves_the_density_matched_figure_and_not_the_faa_one() {
        let r = describe(at(15.0, 1013.25, 0.5));
        assert!(r.density_altitude.faa_rule.m.abs() < 1e-9, "{r:?}");
        let matched = r.density_altitude.density_matched.m;
        assert!(
            (matched - 33.10).abs() < 0.2,
            "50% RH at ICAO standard should be 33.1 m density-matched, got {matched}"
        );
    }

    #[test]
    fn density_matched_rises_monotonically_with_humidity() {
        let mut last = f64::NEG_INFINITY;
        for rh in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let m = describe(at(15.0, 1013.25, rh))
                .density_altitude
                .density_matched
                .m;
            assert!(m > last, "RH {rh} gave {m}, not above {last}");
            last = m;
        }
    }

    /// The two definitions part company in DRY air away from standard, which is why their
    /// difference must never be labeled a humidity correction. Ideal-gas ISA says 30 °C dry at
    /// 1013.25 hPa is ~525.5 m, CIPM-2007 (whose compressibility factor also moves with
    /// temperature) 526.9 m, and the FAA rule 548.6 m.
    #[test]
    fn the_definitions_differ_in_dry_air_away_from_standard() {
        let r = describe(at(30.0, 1013.25, 0.0));
        let faa = r.density_altitude.faa_rule.m;
        let matched = r.density_altitude.density_matched.m;
        assert!((faa - 548.64).abs() < 0.05, "FAA rule at 30 C dry: {faa}");
        assert!(
            (matched - 526.92).abs() < 0.2,
            "density-matched at 30 C dry: {matched}"
        );
        assert!(
            faa - matched > 15.0,
            "expected the FAA rule to read ~22 m high, got {faa} vs {matched}"
        );
    }

    /// Pins the "bit-identical to what the DOPE card has always printed" claim against the
    /// formula as it was written out in `pdf_dope_card` before it moved to `atmosphere` —
    /// comparing against the function the service calls would test nothing.
    #[test]
    fn faa_rule_is_bit_identical_to_the_card_formula_it_replaced() {
        for (t_c, p_hpa) in [
            (15.0, 1013.25),
            (30.0, 1013.25),
            (15.0, 950.0),
            (0.0, 1013.25),
            (-10.0, 850.0),
            (41.3, 702.7),
        ] {
            let r = describe(at(t_c, p_hpa, 0.5));
            // The formula is what is pinned, so feed it the inputs AFTER the service's own
            // round trip through kelvin and pascals (41.3 °C does not survive that exactly).
            let t_c = (t_c + KELVIN_OFFSET_C) - KELVIN_OFFSET_C;
            let p_hpa = (p_hpa * PASCALS_PER_HECTOPASCAL) / PASCALS_PER_HECTOPASCAL;
            let t_f = t_c * 9.0 / 5.0 + 32.0;
            let pressure_alt = 145_366.45 * (1.0 - (p_hpa / 1013.25_f64).powf(0.190_284));
            let isa_temp_f = 59.0 - (pressure_alt / 1000.0) * 3.57;
            let expected = pressure_alt + (120.0 * 5.0 / 9.0) * (t_f - isa_temp_f);
            assert_eq!(
                r.density_altitude.faa_rule.ft.to_bits(),
                expected.to_bits(),
                "{t_c} C {p_hpa} hPa"
            );
        }
    }

    /// `air_density_kg_m3` is the solver's own muzzle density, not a second model.
    #[test]
    fn air_density_is_the_solvers_density() {
        let r = describe(at(15.0, 1013.25, 0.5));
        let (solver_density, _) =
            crate::atmosphere::calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 50.0);
        assert_eq!(r.air_density_kg_m3.to_bits(), solver_density.to_bits());
        assert!(
            (r.air_density_kg_m3 - 1.2216).abs() < 0.0002,
            "{}",
            r.air_density_kg_m3
        );
    }

    /// Humidity is a FRACTION on this wire, as in solve-json v1. 50 means 5000%, and is refused
    /// with the solve's own path rather than silently read as 50%.
    #[test]
    fn humidity_in_percent_is_refused_at_its_path() {
        let error = refusal(at(15.0, 1013.25, 50.0));
        assert_eq!(error.path(), Some("$.atmosphere.relative_humidity"));
    }

    /// Omitting everything is legal, resolves to ICAO standard at sea level with 50% RH, and
    /// SAYS so — the same four assumptions a solve would announce.
    #[test]
    fn omitted_fields_resolve_like_a_solve_and_are_announced() {
        let r = describe(AtmosphereV1::default());
        assert_eq!(r.atmosphere.altitude_m, 0.0);
        assert!((r.atmosphere.temperature_k - 288.15).abs() < 1e-9);
        assert!((r.atmosphere.pressure_pa - 101_325.0).abs() < 1e-6);
        assert_eq!(r.atmosphere.relative_humidity, 0.5);
        let paths: Vec<_> = r
            .assumptions
            .iter()
            .filter_map(|a| a.path.as_deref())
            .collect();
        for expected in [
            "$.atmosphere.altitude_m",
            "$.atmosphere.temperature_k",
            "$.atmosphere.pressure_pa",
            "$.atmosphere.relative_humidity",
        ] {
            assert!(
                paths.contains(&expected),
                "missing assumption at {expected}: {paths:?}"
            );
        }
    }

    /// A QNH is reduced to station pressure BEFORE either density altitude is computed —
    /// density altitude of the sea, rather than of the firing point, is the MBA-643 defect.
    #[test]
    fn a_qnh_is_reduced_to_station_pressure_first() {
        let altitude_m = 1500.0;
        let mut qnh = at(15.0, 1013.25, 0.0);
        qnh.altitude_m = Some(altitude_m);
        qnh.pressure_reference = Some(PressureReferenceV1::Qnh);
        let r = describe(qnh);
        assert!(
            r.atmosphere.pressure_pa < 90_000.0,
            "QNH 1013.25 at 1500 m should reduce to ~845 hPa station, got {} Pa",
            r.atmosphere.pressure_pa
        );
        assert!(
            (r.pressure_altitude.m - altitude_m).abs() < 30.0,
            "pressure altitude of a standard QNH should be ~the station altitude: {}",
            r.pressure_altitude.m
        );
        assert!(r
            .assumptions
            .iter()
            .any(|a| a.path.as_deref() == Some("$.atmosphere.pressure_pa")));
    }

    #[test]
    fn meters_and_feet_agree() {
        let r = describe(at(30.0, 950.0, 0.5));
        for altitude in [
            r.pressure_altitude,
            r.density_altitude.faa_rule,
            r.density_altitude.density_matched,
        ] {
            assert!(
                (altitude.m - altitude.ft * METERS_PER_FOOT).abs() < 1e-9,
                "{altitude:?}"
            );
        }
    }

    /// °C sent as `temperature_k` is the likeliest temperature mistake, and it used to come
    /// back `ok` with a density of 23 kg/m³. So did 400 K, and 1e-300 K came back with a
    /// `null` density.
    #[test]
    fn a_temperature_no_air_can_have_is_refused_at_temperature_k() {
        for kelvin in [15.0, 5.0, 1e-300, 173.0, 373.2, 400.0, 1e300] {
            let mut atmosphere = at(15.0, 1013.25, 0.5);
            atmosphere.temperature_k = Some(kelvin);
            let error = refusal(atmosphere);
            assert_eq!(
                error.path(),
                Some("$.atmosphere.temperature_k"),
                "{kelvin} K: {error:?}"
            );
        }
    }

    /// hPa sent as `pressure_pa` — a pressure altitude of about 30 km — is refused at the
    /// pressure, not blamed on the humidity even on a day warm enough that the vapor check
    /// would also trip.
    #[test]
    fn a_pressure_outside_the_troposphere_is_refused_at_pressure_pa() {
        for pascals in [1013.25, 101.325, 29.92, 1e-300, 250_000.0, 1e300] {
            let mut atmosphere = at(25.0, 1013.25, 0.9);
            atmosphere.pressure_pa = Some(pascals);
            let error = refusal(atmosphere);
            assert_eq!(
                error.path(),
                Some("$.atmosphere.pressure_pa"),
                "{pascals} Pa: {error:?}"
            );
        }
    }

    /// With the pressure defaulted, a station above the troposphere is the altitude's fault.
    #[test]
    fn a_defaulted_pressure_above_the_troposphere_is_refused_at_altitude_m() {
        let error = refusal(AtmosphereV1 {
            altitude_m: Some(20_000.0),
            ..AtmosphereV1::default()
        });
        assert_eq!(error.path(), Some("$.atmosphere.altitude_m"));
    }

    /// A plausible QNH reduced at an impossible altitude (feet sent as meters) is the
    /// altitude's fault; an implausible QNH is the pressure's. Either way the message quotes
    /// the QNH the caller sent, not the station pressure derived from it.
    #[test]
    fn a_qnh_refusal_blames_the_input_that_is_wrong_and_quotes_what_was_sent() {
        let mut feet_as_meters = at(15.0, 1013.25, 0.5);
        feet_as_meters.altitude_m = Some(36_000.0);
        feet_as_meters.pressure_reference = Some(PressureReferenceV1::Qnh);
        let error = refusal(feet_as_meters);
        assert_eq!(error.path(), Some("$.atmosphere.altitude_m"), "{error:?}");
        assert!(
            error.message.contains("a QNH of 101325 Pa"),
            "{}",
            error.message
        );

        let mut hpa_as_pa = at(15.0, 1013.25, 0.5);
        hpa_as_pa.altitude_m = Some(1500.0);
        hpa_as_pa.pressure_pa = Some(1013.25);
        hpa_as_pa.pressure_reference = Some(PressureReferenceV1::Qnh);
        let error = refusal(hpa_as_pa);
        assert_eq!(error.path(), Some("$.atmosphere.pressure_pa"), "{error:?}");
        assert!(
            error.message.contains("a QNH of 1013.25 Pa"),
            "{}",
            error.message
        );
    }

    /// An absurd input is echoed in exponent form, not as hundreds of digits.
    #[test]
    fn absurd_inputs_are_echoed_compactly() {
        let mut atmosphere = at(15.0, 1013.25, 0.5);
        atmosphere.temperature_k = Some(1e300);
        let error = refusal(atmosphere);
        assert!(error.message.starts_with("1e300 K"), "{}", error.message);
        assert!(error.message.len() < 300, "{}", error.message);
    }

    /// Saturated air at 100 °C needs more vapor pressure than a standard sea-level
    /// atmosphere's total pressure: that air does not exist, and computing it as pure steam is
    /// not an answer.
    #[test]
    fn humidity_beyond_the_total_pressure_is_refused_at_relative_humidity() {
        let error = refusal(at(100.0, 1013.25, 1.0));
        assert_eq!(
            error.path(),
            Some("$.atmosphere.relative_humidity"),
            "{error:?}"
        );
        // The same air at 90% is merely unusual.
        describe(at(100.0, 1013.25, 0.9));
        // And the hottest, wettest sea-level day on record is fine.
        describe(at(45.0, 1013.25, 1.0));
    }

    /// Anything accepted is finite, with a positive density: sweep the whole accepted domain.
    #[test]
    fn every_accepted_atmosphere_reports_finite_figures() {
        let mut accepted = 0;
        for t in (0..=20).map(|i| MIN_TEMPERATURE_K + 10.0 * f64::from(i)) {
            for p in (0..=16).map(|i| 20_000.0 + 10_000.0 * f64::from(i)) {
                for rh in [0.0, 0.3, 0.7, 1.0] {
                    let atmosphere = AtmosphereV1 {
                        altitude_m: Some(0.0),
                        temperature_k: Some(t),
                        pressure_pa: Some(p),
                        relative_humidity: Some(rh),
                        ..AtmosphereV1::default()
                    };
                    let Ok(r) = density_altitude_v1(&DensityAltitudeRequestV1 { atmosphere })
                    else {
                        continue;
                    };
                    accepted += 1;
                    assert!(r.air_density_kg_m3.is_finite() && r.air_density_kg_m3 > 0.0);
                    for altitude in [
                        r.pressure_altitude,
                        r.density_altitude.faa_rule,
                        r.density_altitude.density_matched,
                    ] {
                        assert!(altitude.m.is_finite() && altitude.ft.is_finite(), "{r:?}");
                        assert!(within_troposphere(altitude.m), "{r:?}");
                    }
                }
            }
        }
        assert!(accepted > 300, "only {accepted} accepted");
    }

    /// Malformed objects are refused by the solve's own shape validator, with its codes and
    /// paths — the same answer `solve` gives for the same `atmosphere`.
    #[test]
    fn malformed_requests_get_the_solves_codes_and_paths() {
        let cases = [
            (
                json!({"atmosphere": {"temperature_c": 30.0}}),
                SolveErrorCodeV1::UnknownField,
                "$.atmosphere.temperature_c",
            ),
            (
                json!({"atmosphere": {"temperature_k": null}}),
                SolveErrorCodeV1::InvalidValue,
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"atmosphere": {"pressure_pa": "101325"}}),
                SolveErrorCodeV1::InvalidValue,
                "$.atmosphere.pressure_pa",
            ),
            (
                json!({"atmosphere": {"pressure_reference": "altimeter"}}),
                SolveErrorCodeV1::InvalidValue,
                "$.atmosphere.pressure_reference",
            ),
            (json!({}), SolveErrorCodeV1::MissingField, "$.atmosphere"),
            (
                json!({"atmosphere": {}, "units": "si"}),
                SolveErrorCodeV1::UnknownField,
                "$.units",
            ),
            (json!([{}]), SolveErrorCodeV1::InvalidValue, "$"),
        ];
        for (request, code, path) in cases {
            let error = decode_density_altitude_request_v1(&request)
                .expect_err("malformed")
                .error;
            assert_eq!(error.code, code, "{request}: {error:?}");
            assert_eq!(error.path(), Some(path), "{request}: {error:?}");
        }
        assert!(decode_density_altitude_request_v1(&json!({"atmosphere": {}})).is_ok());
    }

    /// The air outside this command's own model range that a solve refuses: the same field is
    /// to blame in both. Saturated air at 380 K, 400 K at 90 %, -103 °C at a pressure
    /// of 0.0001 Pa, and a station above the troposphere because its pressure was omitted or
    /// was a plausible QNH reduced 20 km up.
    #[test]
    fn a_solve_blames_impossible_air_on_the_reading_this_command_blames() {
        let solve_path = |atmosphere: Value| {
            let mut request = json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0},
                "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
                "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
                "solver": {"method": "rk45", "time_step_s": 0.001},
                "effects": {"magnus": false, "coriolis": false, "enhanced_spin_drift": false},
                "sampling": {"interval_m": 50.0}
            });
            request["atmosphere"] = atmosphere;
            let request = crate::solve_json::decode_solve_request_v1(&request.to_string())
                .expect("a well-formed solve request");
            let error = crate::solve_v1::prepare_request(&request)
                .expect_err("the solve refuses this air")
                .error;
            (error.code, error.path().map(str::to_string))
        };
        for (atmosphere, path) in [
            (
                json!({"temperature_k": 380.0, "pressure_pa": 101325.0, "relative_humidity": 1.0}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"temperature_k": 400.0, "pressure_pa": 101325.0, "relative_humidity": 0.9}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"temperature_k": 170.0, "pressure_pa": 0.0001, "relative_humidity": 1.0}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"altitude_m": 20000.0, "temperature_k": 310.0, "relative_humidity": 1.0}),
                "$.atmosphere.altitude_m",
            ),
            (
                json!({"altitude_m": 20000.0, "temperature_k": 310.0, "pressure_pa": 101325.0,
                       "pressure_reference": "qnh", "relative_humidity": 1.0}),
                "$.atmosphere.altitude_m",
            ),
        ] {
            let request = decode_density_altitude_request_v1(&json!({ "atmosphere": atmosphere }))
                .expect("a well-formed atmosphere");
            let here = density_altitude_v1(&request).expect_err("refused").error;
            assert_eq!(here.path(), Some(path), "{atmosphere}: {here:?}");
            assert_eq!(
                solve_path(atmosphere.clone()),
                (here.code, here.path().map(str::to_string)),
                "{atmosphere}: the solve blames another field"
            );
        }
    }

    /// Whenever a solve refuses an `atmosphere` object, this command refuses the same object
    /// with the same code and path, across every way the object can be built: each
    /// field sent or omitted, an absolute pressure or a QNH, temperatures and pressures from
    /// absurd through ordinary to absurd again, and both edges of this command's own model
    /// range to the bit. The message is the same word for word too, except where this
    /// command's own range checks (the temperature range, the troposphere) refuse the air first.
    #[test]
    fn every_atmosphere_a_solve_refuses_is_refused_here_with_the_same_code_and_path() {
        let up = |x: f64| f64::from_bits(x.to_bits() + 1);
        let down = |x: f64| f64::from_bits(x.to_bits() - 1);
        let temperatures = [
            None,
            Some(1e-300),
            Some(5.0),
            Some(100.0),
            Some(down(MIN_TEMPERATURE_K)),
            Some(MIN_TEMPERATURE_K),
            Some(250.0),
            Some(300.0),
            Some(330.0),
            Some(360.0),
            Some(MAX_TEMPERATURE_K),
            Some(up(MAX_TEMPERATURE_K)),
            Some(380.0),
            Some(400.0),
            Some(500.0),
            Some(647.0),
            Some(700.0),
            Some(1e300),
        ];
        let pressures = [
            None,
            Some(1e-300),
            Some(1e-4),
            Some(0.01),
            Some(29.92),
            Some(1013.25),
            Some(5_000.0),
            Some(20_000.0),
            Some(50_000.0),
            Some(101_325.0),
            Some(250_000.0),
            Some(1e300),
        ];
        let altitudes = [
            None,
            Some(0.0),
            Some(1500.0),
            Some(12_000.0),
            Some(20_000.0),
            Some(40_000.0),
        ];
        let humidities = [None, Some(0.0), Some(0.5), Some(0.9), Some(1.0), Some(1.06)];
        let references = [None, Some("absolute"), Some("qnh")];

        let base = json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0},
            "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
            "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
            "solver": {"method": "rk45", "time_step_s": 0.001},
            "effects": {"magnus": false, "coriolis": false, "enhanced_spin_drift": false},
            "sampling": {"interval_m": 50.0}
        });
        let mut compared = 0usize;
        let mut blamed = std::collections::BTreeMap::<String, usize>::new();
        for temperature in temperatures {
            for pressure in pressures {
                for altitude in altitudes {
                    for humidity in humidities {
                        for reference in references {
                            let mut atmosphere = serde_json::Map::new();
                            for (name, value) in [
                                ("temperature_k", temperature),
                                ("pressure_pa", pressure),
                                ("altitude_m", altitude),
                                ("relative_humidity", humidity),
                            ] {
                                if let Some(value) = value {
                                    atmosphere.insert(name.to_string(), json!(value));
                                }
                            }
                            if let Some(reference) = reference {
                                atmosphere
                                    .insert("pressure_reference".to_string(), json!(reference));
                            }
                            let atmosphere = Value::Object(atmosphere);
                            let mut solve = base.clone();
                            solve["atmosphere"] = atmosphere.clone();
                            let solve =
                                crate::solve_json::decode_solve_request_v1(&solve.to_string())
                                    .and_then(|request| crate::solve_v1::prepare_request(&request));
                            let Err(solve) = solve else { continue };
                            let solve = solve.error;
                            let here = decode_density_altitude_request_v1(
                                &json!({ "atmosphere": atmosphere }),
                            )
                            .and_then(|request| density_altitude_v1(&request))
                            .expect_err("a solve refuses this atmosphere, so this command must")
                            .error;
                            assert_eq!(
                                (here.code, here.path()),
                                (solve.code, solve.path()),
                                "{atmosphere}\n  solve: {}\n  here:  {}",
                                solve.message,
                                here.message
                            );
                            if here.message != solve.message {
                                assert!(
                                    here.message.contains("the range the humid-air model")
                                        || here
                                            .message
                                            .contains("density altitude is computed only between"),
                                    "{atmosphere}: the messages differ, but not because this \
                                     command's own range check refused first\n  solve: {}\n  \
                                     here:  {}",
                                    solve.message,
                                    here.message
                                );
                            }
                            compared += 1;
                            if solve.message.contains("no such air exists") {
                                *blamed
                                    .entry(solve.path().unwrap_or_default().to_string())
                                    .or_default() += 1;
                            }
                        }
                    }
                }
            }
        }
        // The grid reaches every blame the impossible-air rule makes.
        for path in [
            "$.atmosphere.relative_humidity",
            "$.atmosphere.temperature_k",
            "$.atmosphere.pressure_pa",
            "$.atmosphere.altitude_m",
        ] {
            assert!(
                blamed.get(path).copied().unwrap_or(0) > 0,
                "no impossible air blamed at {path}: {blamed:?}"
            );
        }
        assert!(compared > 1000, "{compared} refusals compared");
    }
}
