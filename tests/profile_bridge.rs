//! Integration tests for the bridge's profile commands (C2): `profile.validate`,
//! `profile.normalize`, and — when the `profile-import` feature is compiled in —
//! `profile.import_a7p`. The import tests reuse the same spec-derived synthetic `.a7p`
//! encoders as `tests/a7p_import_cli.rs` (the repo has no binary `.a7p` fixture on
//! purpose; documents are synthesized from the wire spec, exactly as the parser's own
//! unit tests do) and assert the bridge surfaces the identical mapping the CLI performs.

#![cfg(feature = "bridge")]

use ballistics_engine::bridge::{bridge_call, BRIDGE_API_VERSION};
use serde_json::{json, Value};

fn call(envelope: Value) -> Value {
    let raw = bridge_call(&envelope.to_string());
    serde_json::from_str(&raw).expect("bridge output must be valid JSON")
}

/// A saved-profile document as it exists on disk (`~/.ballistics/profiles/*.json`),
/// carrying the major optional surfaces plus an unknown key a future engine might write.
fn fixture_profile() -> Value {
    json!({
        "name": "bridge-fixture",
        "velocity": 762.0,
        "bc": 0.243,
        "mass": 11.33980925,
        "diameter": 7.8232,
        "drag_model": "G7",
        "units": "metric",
        "temperature": 15.0,
        "pressure": 1013.25,
        "humidity": 50.0,
        "altitude": 0.0,
        "bc_segments": [
            {"bc": 0.243, "velocity_mps": 792.0},
            {"bc": 0.230, "velocity_mps": 400.0}
        ],
        "dsf_points": [{"mach": 0.9, "dsf": 1.04}],
        "elevation_cf": 0.97,
        "windage_cf": 1.02,
        "elevation_click": "0.1mil",
        "zero_sets": [
            {"name": "suppressed", "zero_distance": 200.0, "poi_up_mil": -0.3,
             "poi_right_mil": 0.1, "notes": "suppressed load"}
        ],
        "reticle": {
            "name": "mil-grid",
            "focal_plane": "ffp",
            "reference_magnification": 10.0,
            "marks": [{"down_mil": 0.0, "right_mil": 0.0, "kind": "center"}]
        },
        "future_field_this_engine_never_wrote": true
    })
}

#[test]
fn validate_accepts_the_fixture_and_normalizes_it() {
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.validate",
        "request": fixture_profile(),
    }));
    assert_eq!(out["ok"], true, "{out}");
    let result = &out["result"];
    assert_eq!(result["valid"], true, "{result}");
    assert_eq!(result["warnings"].as_array().unwrap().len(), 0, "{result}");

    // The normalized document is this engine's own serialization: stored fields survive,
    // defaults fill in, and the unknown key is dropped (the CLI's own load-then-save
    // behavior — documented, not silent: validate exists precisely to preview this).
    let normalized = &result["normalized"];
    assert_eq!(normalized["name"], "bridge-fixture");
    assert_eq!(normalized["bc_segments"][0]["velocity_mps"], 792.0);
    assert_eq!(normalized["zero_sets"][0]["name"], "suppressed");
    assert_eq!(normalized["reticle"]["marks"][0]["kind"], "center");
    assert_eq!(normalized["elevation_click"], "0.1mil");
    assert!(normalized.get("future_field_this_engine_never_wrote").is_none());
}

#[test]
fn validate_reports_the_cli_load_gates_as_warnings() {
    let mut profile = fixture_profile();
    profile["elevation_cf"] = json!(0.4); // outside the (0.5, 1.5) tracking-CF band
    profile["units"] = json!("nautical");
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.validate",
        "request": profile,
    }));
    assert_eq!(out["ok"], true, "{out}");
    assert_eq!(out["result"]["valid"], false, "{out}");
    let warnings: Vec<String> =
        serde_json::from_value(out["result"]["warnings"].clone()).unwrap();
    assert!(warnings.iter().any(|w| w.contains("elevation_cf")), "{warnings:?}");
    assert!(warnings.iter().any(|w| w.contains("unsupported units")), "{warnings:?}");
}

#[test]
fn normalize_round_trips_and_is_idempotent() {
    let first = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.normalize",
        "request": fixture_profile(),
    }));
    assert_eq!(first["ok"], true, "{first}");
    let normalized = first["result"]["profile"].clone();
    assert_eq!(normalized["velocity"], 762.0);
    assert_eq!(normalized["dsf_points"][0]["dsf"], 1.04);

    // Feeding the engine's own output back in must be the identity.
    let second = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.normalize",
        "request": normalized.clone(),
    }));
    assert_eq!(second["ok"], true, "{second}");
    assert_eq!(second["result"]["profile"], normalized);
}

#[test]
fn a_non_profile_document_is_an_invalid_request() {
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.validate",
        "request": {"name": "missing-everything"},
    }));
    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(out["error"]["code"], "invalid_request");
}

/// Feature-detection contract: `profile.import_a7p` appears in `meta.capabilities`
/// exactly when the feature is compiled in, and is `unknown_command` when it is not.
#[test]
fn import_command_presence_matches_compiled_features() {
    let caps = call(json!({"api_version": BRIDGE_API_VERSION, "command": "meta.capabilities"}));
    let commands: Vec<String> =
        serde_json::from_value(caps["result"]["commands"].clone()).unwrap();
    assert_eq!(
        commands.contains(&"profile.import_a7p".to_string()),
        cfg!(feature = "profile-import")
    );
    assert!(commands.contains(&"profile.validate".to_string()));
    assert!(commands.contains(&"profile.normalize".to_string()));
}

#[cfg(not(feature = "profile-import"))]
#[test]
fn import_without_the_feature_is_unknown_command() {
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "profile.import_a7p",
        "request": {"a7p_base64": "AAAA"},
    }));
    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(out["error"]["code"], "unknown_command");
}

#[cfg(feature = "profile-import")]
mod import {
    use super::*;

    // Spec-derived encoders (same shape as tests/a7p_import_cli.rs and the parser's own
    // unit tests — synthesized documents, no binary fixture in the repo).
    fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                return;
            }
            out.push(byte | 0x80);
        }
    }
    fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
        enc_varint(u64::from(number) << 3, out);
        enc_varint(value as u64, out);
    }
    fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(s.len() as u64, out);
        out.extend_from_slice(s.as_bytes());
    }
    fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        enc_varint((u64::from(number) << 3) | 2, out);
        enc_varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    /// The same synthetic 6.5 CM G7 export as `tests/a7p_import_cli.rs`, plus an unknown
    /// future field (number 63) the report must surface.
    fn synthetic_a7p() -> Vec<u8> {
        let mut p = Vec::new();
        enc_str(1, "BRIDGE 6.5CM", &mut p);
        enc_str(3, "140 ELD-M", &mut p);
        enc_i32(9, 50, &mut p); // 50 mm sight height
        enc_i32(10, 800, &mut p); // 8.00 in twist
        enc_i32(11, 8230, &mut p); // 823.0 m/s
        enc_i32(15, 15, &mut p);
        enc_i32(16, 10132, &mut p); // 1013.2 hPa
        enc_i32(17, 50, &mut p);
        enc_i32(20, 264, &mut p); // 0.264 in
        enc_i32(21, 1400, &mut p); // 140 gr
        enc_i32(22, 1360, &mut p); // 1.36 in
        enc_i32(24, 1, &mut p); // G7
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed); // one distance: 100.00 m (zero index defaults 0)
        enc_bytes(26, &packed, &mut p);
        let mut row = Vec::new();
        enc_i32(1, 3260, &mut row); // G7 BC 0.326
        enc_i32(2, 8230, &mut row);
        enc_bytes(27, &row, &mut p);
        enc_i32(63, 42, &mut p); // unknown future field, must be reported
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        ballistics_engine::profile_import::wrap_payload(&payload)
    }

    fn to_base64(bytes: &[u8]) -> String {
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for chunk in bytes.chunks(3) {
            let b = [chunk[0], *chunk.get(1).unwrap_or(&0), *chunk.get(2).unwrap_or(&0)];
            let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
            out.push(char::from(ALPHABET[(n >> 18 & 63) as usize]));
            out.push(char::from(ALPHABET[(n >> 12 & 63) as usize]));
            out.push(if chunk.len() > 1 { char::from(ALPHABET[(n >> 6 & 63) as usize]) } else { '=' });
            out.push(if chunk.len() > 2 { char::from(ALPHABET[(n & 63) as usize]) } else { '=' });
        }
        out
    }

    #[test]
    fn import_maps_the_document_end_to_end_and_matches_the_library_mapping() {
        let file = synthetic_a7p();
        let out = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(&file)},
        }));
        assert_eq!(out["ok"], true, "{out}");
        let result = &out["result"];
        let profile = &result["profile"];

        // The mapped fields, as the CLI import would store them (metric, SI-pinned).
        assert_eq!(profile["name"], "BRIDGE 6.5CM");
        assert_eq!(profile["units"], "metric");
        assert_eq!(profile["drag_model"], "G7");
        assert_eq!(profile["velocity"], 823.0);
        assert_eq!(profile["bc"], 0.326);
        assert_eq!(profile["bullet_name"], "140 ELD-M");
        assert!((profile["mass"].as_f64().unwrap() - 140.0 * 0.06479891).abs() < 1e-9);
        assert!((profile["diameter"].as_f64().unwrap() - 0.264 * 25.4).abs() < 1e-9);
        assert_eq!(profile["sight_height"], 50.0);
        assert_eq!(profile["zero_distance"], 100.0);
        assert_eq!(profile["twist_right"], true);

        // The report surfaces everything the CLI's report carries.
        assert!(result["mapped"].as_array().unwrap().iter().any(|row| row[0] == "b_weight"));
        assert!(!result["unmapped"].as_array().unwrap().is_empty());
        let unknown = result["unknown_fields"].as_array().unwrap();
        assert!(
            unknown.iter().any(|u| u["context"] == "Profile" && u["number"] == 63),
            "{unknown:?}"
        );

        // Bridge result must equal the library mapping verbatim — one mapping, two doors.
        let doc = ballistics_engine::profile_import::parse_a7p(&file).unwrap();
        let reference =
            ballistics_engine::profile_import::map_a7p_to_profile(&doc, None, None).unwrap();
        let mut reference_profile = serde_json::to_value(&reference.profile).unwrap();
        let mut bridge_profile = profile.clone();
        // `created` is a wall-clock timestamp; everything else must be identical.
        reference_profile.as_object_mut().unwrap().remove("created");
        bridge_profile.as_object_mut().unwrap().remove("created");
        assert_eq!(bridge_profile, reference_profile);

        // And the imported document passes profile.validate cleanly.
        let validated = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.validate",
            "request": profile,
        }));
        assert_eq!(validated["ok"], true, "{validated}");
        assert_eq!(validated["result"]["valid"], true, "{validated}");
    }

    #[test]
    fn zero_click_is_honored_like_the_cli_flag() {
        // Zeroing state: zero_x = -20_000 (20 clicks RIGHT after upstream's negation),
        // zero_y = 10_000 (10 clicks UP) — the sign conventions proven in the library's
        // own mapping tests.
        let mut p = Vec::new();
        enc_str(1, "ZEROED", &mut p);
        enc_i32(7, -20_000, &mut p);
        enc_i32(8, 10_000, &mut p);
        enc_i32(11, 7920, &mut p);
        enc_i32(20, 338, &mut p);
        enc_i32(21, 3000, &mut p);
        enc_i32(24, 0, &mut p);
        let mut row = Vec::new();
        enc_i32(1, 7160, &mut row);
        enc_i32(2, 7920, &mut row);
        enc_bytes(27, &row, &mut p);
        let mut packed = Vec::new();
        enc_varint(10_000, &mut packed); // 100.00 m zero
        enc_bytes(26, &packed, &mut p);
        let mut payload = Vec::new();
        enc_bytes(1, &p, &mut payload);
        let file = ballistics_engine::profile_import::wrap_payload(&payload);

        let out = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(&file), "zero_click": "0.1mil"},
        }));
        assert_eq!(out["ok"], true, "{out}");
        let profile = &out["result"]["profile"];
        // up: 10 clicks * 0.1/1000 rad * 100 m = 0.1 m; right: 20 clicks -> 0.2 m.
        assert!((profile["zero_poi_up_m"].as_f64().unwrap() - 0.1).abs() < 1e-12);
        assert!((profile["zero_poi_right_m"].as_f64().unwrap() - 0.2).abs() < 1e-12);
        assert_eq!(profile["elevation_click"], "0.1mil");
        assert_eq!(profile["zero_sets"][0]["name"], "a7p-zero");

        let bad = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(&file), "zero_click": "not-a-click"},
        }));
        assert_eq!(bad["error"]["code"], "invalid_request", "{bad}");
    }

    #[test]
    fn md5_mismatch_warns_by_default_and_rejects_under_strict() {
        let mut file = synthetic_a7p();
        // Flip one hex digit of the MD5 envelope prefix: payload parses, checksum fails.
        file[0] = if file[0] == b'0' { b'1' } else { b'0' };

        let lenient = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(&file)},
        }));
        assert_eq!(lenient["ok"], true, "{lenient}");
        let warnings: Vec<String> =
            serde_json::from_value(lenient["result"]["warnings"].clone()).unwrap();
        assert!(
            warnings.iter().any(|w| w.contains("checksum mismatch")),
            "{warnings:?}"
        );

        let strict = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(&file), "strict": true},
        }));
        assert_eq!(strict["ok"], false, "{strict}");
        assert_eq!(strict["error"]["code"], "command_failed");
        assert!(
            strict["error"]["message"]
                .as_str()
                .unwrap()
                .contains("checksum mismatch"),
            "{strict}"
        );
    }

    #[test]
    fn garbage_and_bad_base64_are_clean_errors() {
        let garbage = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": to_base64(b"not a profile")},
        }));
        assert_eq!(garbage["ok"], false, "{garbage}");
        assert_eq!(garbage["error"]["code"], "command_failed");
        assert!(
            garbage["error"]["message"].as_str().unwrap().contains("a7p"),
            "{garbage}"
        );

        let bad_text = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": "@@@not base64@@@"},
        }));
        assert_eq!(bad_text["error"]["code"], "invalid_request", "{bad_text}");

        let unknown_key = call(json!({
            "api_version": BRIDGE_API_VERSION,
            "command": "profile.import_a7p",
            "request": {"a7p_base64": "AAAA", "surprise": 1},
        }));
        assert_eq!(unknown_key["error"]["code"], "invalid_request", "{unknown_key}");
    }
}
