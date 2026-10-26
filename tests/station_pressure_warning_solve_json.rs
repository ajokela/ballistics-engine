//! `station_pressure_implausible_for_altitude`: an explicit absolute `pressure_pa` that cannot
//! be a station reading at `altitude_m`.
//!
//! An explicit absolute pressure is authoritative, so the solve flies it as sent. A bridge app
//! that always sent 101325 Pa beside the real `altitude_m` of Cusco (3310.128 m) flew sea-level
//! air at 3310 m with nothing in the response to say so. The warning is that signal. It is a
//! notice only: the solve still uses the pressure as sent.
//!
//! The band is 0.85 to 1.08 times the ICAO standard pressure at the altitude
//! (`atmosphere::station_pressure_plausibility`). These cases pin where it fires and where it
//! must not:
//!
//!   - 101325 Pa at Cusco is 1.503 times standard: warned, high side, naming the placeholder
//!     and QNH causes;
//!   - a real Cusco reading (68830 Pa) is close to standard: no warning;
//!   - the same reading with `altitude_m` omitted or 0 is 0.679 times standard: warned, low
//!     side, naming the altitude and the unit but not QNH;
//!   - a QNH, sent as `"qnh"`, and an omitted pressure are never warned;
//!   - an explicit `"absolute"` reference is warned exactly like the default.

use ballistics_engine::solve_json::{decode_solve_request_v1, SolveNoticeV1, SolveSuccessV1};
use ballistics_engine::solve_v1;
use ballistics_engine::solve_v1::WARNING_STATION_PRESSURE_IMPLAUSIBLE_FOR_ALTITUDE as CODE;

const CUSCO_M: f64 = 3310.128;

/// A .308 175 gr to 600 m, zeroed at 100 m, at 16 °C and 52% RH. `atmosphere_extra` holds the
/// fields that vary: altitude, pressure and pressure reference.
fn request_json(atmosphere_extra: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.0315,
                            "drag_model": "G7", "ballistic_coefficient": 0.243}},
            "rifle": {{"muzzle_velocity_mps": 792.48, "sight_height_m": 0.0508,
                       "twist_rate_m_per_turn": 0.254}},
            "shot": {{"max_range_m": 600.0, "zero_distance_m": 100.0}},
            "atmosphere": {{"temperature_k": 289.15, "relative_humidity": 0.52{atmosphere_extra}}},
            "wind": {{}},
            "solver": {{"method": "rk45"}},
            "effects": {{}},
            "sampling": {{"interval_m": 100.0}}
        }}"#
    )
}

fn solve(atmosphere_extra: &str) -> SolveSuccessV1 {
    let request = decode_solve_request_v1(&request_json(atmosphere_extra)).expect("decodes");
    solve_v1(request).expect("solves")
}

fn pressure_warnings(response: &SolveSuccessV1) -> Vec<&SolveNoticeV1> {
    response
        .warnings
        .iter()
        .filter(|n| n.code == CODE)
        .collect()
}

fn only_pressure_warning(response: &SolveSuccessV1) -> &SolveNoticeV1 {
    let warnings = pressure_warnings(response);
    assert_eq!(warnings.len(), 1, "{:?}", response.warnings);
    assert_eq!(
        warnings[0].path.as_deref(),
        Some("$.atmosphere.pressure_pa")
    );
    warnings[0]
}

#[test]
fn a_sea_level_placeholder_at_cusco_is_warned_and_flown_as_sent() {
    let response = solve(&format!(
        r#", "altitude_m": {CUSCO_M}, "pressure_pa": 101325.0"#
    ));
    let message = &only_pressure_warning(&response).message;
    for expected in [
        "pressure_pa 101325 Pa",
        "1.503 times",
        "altitude_m 3310.128 m",
        "67402.29 Pa",
        "placeholder",
        "\"pressure_reference\": \"qnh\"",
        "not in meters",
        "used pressure_pa as sent",
    ] {
        assert!(
            message.contains(expected),
            "{expected:?} missing: {message}"
        );
    }
    assert!(!message.contains("defaulted"), "{message}");
    // Notice only: the pressure is flown as sent.
    assert_eq!(response.resolved_request.atmosphere.pressure_pa, 101_325.0);
}

#[test]
fn an_explicit_absolute_reference_is_warned_like_the_default() {
    let response = solve(&format!(
        r#", "altitude_m": {CUSCO_M}, "pressure_pa": 101325.0, "pressure_reference": "absolute""#
    ));
    only_pressure_warning(&response);
}

#[test]
fn a_real_cusco_station_reading_is_not_warned() {
    let response = solve(&format!(
        r#", "altitude_m": {CUSCO_M}, "pressure_pa": 68830.0"#
    ));
    assert!(
        pressure_warnings(&response).is_empty(),
        "{:?}",
        response.warnings
    );
}

#[test]
fn a_station_reading_at_an_omitted_or_zero_altitude_is_warned_on_the_low_side() {
    let omitted = solve(r#", "pressure_pa": 68830.0"#);
    let message = &only_pressure_warning(&omitted).message;
    for expected in [
        "0.679 times",
        "altitude_m 0 m",
        "101325.00 Pa",
        "altitude_m was omitted and defaulted to 0 m",
        "below the 0.85 to 1.08 band",
        "not the station elevation",
        "not in pascals",
    ] {
        assert!(
            message.contains(expected),
            "{expected:?} missing: {message}"
        );
    }
    assert!(
        !message.contains("qnh"),
        "a low reading is never a QNH: {message}"
    );

    let zero = solve(r#", "altitude_m": 0.0, "pressure_pa": 68830.0"#);
    let message = &only_pressure_warning(&zero).message;
    assert!(message.contains("0.679 times"), "{message}");
    assert!(
        !message.contains("defaulted"),
        "altitude_m was sent: {message}"
    );
}

#[test]
fn a_qnh_and_an_omitted_pressure_are_never_warned() {
    let qnh = solve(&format!(
        r#", "altitude_m": {CUSCO_M}, "pressure_pa": 103300.0, "pressure_reference": "qnh""#
    ));
    assert!(pressure_warnings(&qnh).is_empty(), "{:?}", qnh.warnings);

    let omitted = solve(&format!(r#", "altitude_m": {CUSCO_M}"#));
    assert!(
        pressure_warnings(&omitted).is_empty(),
        "{:?}",
        omitted.warnings
    );
}

#[cfg(feature = "bridge")]
#[test]
fn the_bridge_solve_carries_the_warning() {
    use ballistics_engine::bridge::{bridge_call, BRIDGE_API_VERSION};
    let request: serde_json::Value = serde_json::from_str(&request_json(&format!(
        r#", "altitude_m": {CUSCO_M}, "pressure_pa": 101325.0"#
    )))
    .unwrap();
    let envelope = serde_json::json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "solve",
        "request": request,
    });
    let out: serde_json::Value = serde_json::from_str(&bridge_call(&envelope.to_string())).unwrap();
    assert_eq!(out["ok"], true, "{out}");
    let codes: Vec<&str> = out["result"]["warnings"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|w| w["code"].as_str())
        .collect();
    assert_eq!(codes, vec![CODE], "{out}");
}
