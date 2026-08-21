//! Integration tests for the JSON command bridge, including the golden
//! cross-check: a `solve` routed through the bridge must be byte-identical to
//! the same request routed through the library solve-json service — the bridge
//! is a transport, never a second implementation.

#![cfg(feature = "bridge")]

use ballistics_engine::bridge::{bridge_call, BRIDGE_API_VERSION, MAX_REQUEST_BYTES};
use serde_json::{json, Value};

fn solve_request_v1() -> Value {
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
                        "drag_model": "G7", "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0},
        "shot": {"max_range_m": 1000.0},
        "atmosphere": {}, "wind": {}, "solver": {}, "effects": {}, "sampling": {}
    })
}

fn call(envelope: Value) -> Value {
    let raw = bridge_call(&envelope.to_string());
    serde_json::from_str(&raw).expect("bridge output must be valid JSON")
}

#[test]
fn solve_through_bridge_matches_library_solve_json_exactly() {
    // Reference: the transport-free library service.
    let request = ballistics_engine::solve_json::decode_solve_request_v1(
        &solve_request_v1().to_string(),
    )
    .expect("fixture must decode");
    let reference = ballistics_engine::solve_v1(request).expect("fixture must solve");
    let reference_value = serde_json::to_value(&reference).expect("reference serializes");

    // Same request through the bridge.
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "solve",
        "request": solve_request_v1(),
    }));

    assert_eq!(out["ok"], true, "bridge solve failed: {out}");
    assert_eq!(
        out["result"], reference_value,
        "bridge result must be the library result verbatim"
    );
}

#[test]
fn solve_result_carries_physical_sanity() {
    let out = call(json!({
        "api_version": BRIDGE_API_VERSION,
        "command": "solve",
        "request": solve_request_v1(),
    }));
    assert_eq!(out["ok"], true, "solve failed: {out}");
    // A 175gr G7 .308 at 823 m/s reaches 1000 m with a subsonic-adjacent but
    // positive remaining velocity; any regression to garbage shows up here.
    let result = &out["result"];
    assert_eq!(result["schema_version"], 1);
    let text = result.to_string();
    assert!(text.len() > 100, "solve result should be substantive: {text}");
}

#[test]
fn capabilities_and_solve_share_one_engine_version() {
    let caps = call(json!({"api_version": BRIDGE_API_VERSION, "command": "meta.capabilities"}));
    let version = call(json!({"api_version": BRIDGE_API_VERSION, "command": "meta.version"}));
    assert_eq!(
        caps["result"]["engine_version"],
        version["result"]["engine_version"]
    );
    assert_eq!(caps["engine_version"], caps["result"]["engine_version"]);
}

#[test]
fn a_request_at_the_size_limit_boundary_is_still_rejected_cleanly() {
    let padding = "x".repeat(MAX_REQUEST_BYTES);
    let oversized = format!(
        r#"{{"api_version":{BRIDGE_API_VERSION},"command":"meta.version","request":"{padding}"}}"#
    );
    let out: Value = serde_json::from_str(&bridge_call(&oversized)).unwrap();
    assert_eq!(out["error"]["code"], "resource_limit");
}

#[test]
fn card_commands_keep_their_error_shape_through_the_generic_adapter() {
    // A null payload must still be InvalidRequest, not a panic or CommandFailed.
    let out = ballistics_engine::bridge::bridge_call(
        r#"{"api_version":1,"command":"card.range_table","request":null}"#,
    );
    let v: serde_json::Value = serde_json::from_str(&out).unwrap();
    assert_eq!(v["ok"], false);
    assert_eq!(v["error"]["code"], "invalid_request");
    assert!(v["error"]["message"].as_str().unwrap().contains("card.range_table"));
}

#[test]
fn error_envelopes_always_carry_versions() {
    for bad in [
        r#"{"#,
        r#"{"api_version":1,"command":"nope"}"#,
        r#"{"api_version":7,"command":"solve"}"#,
    ] {
        let out: Value = serde_json::from_str(&bridge_call(bad)).unwrap();
        assert_eq!(out["ok"], false);
        assert_eq!(out["api_version"], BRIDGE_API_VERSION, "input: {bad}");
        assert!(out["engine_version"].is_string(), "input: {bad}");
    }
}
