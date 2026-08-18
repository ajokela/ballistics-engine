//! Versioned JSON command bridge for embedded (mobile/FFI) consumers.
//!
//! One entry point, [`bridge_call`], accepts a JSON envelope and returns a JSON
//! envelope. Request semantics live in the transport-free library services
//! (starting with [`crate::solve_v1()`]); this module contains only the envelope
//! contract, command dispatch, and panic containment. The C ABI wrapper lives in
//! [`crate::bridge::ffi`] (feature `ffi`).
//!
//! ## Envelope contract (v1)
//!
//! Request:
//! ```json
//! { "api_version": 1, "command": "solve", "request": { ... } }
//! ```
//!
//! Success response:
//! ```json
//! { "ok": true, "api_version": 1, "engine_version": "0.33.1",
//!   "command": "solve", "result": { ... } }
//! ```
//!
//! Error response (always in-band; the bridge never signals failure any other way):
//! ```json
//! { "ok": false, "api_version": 1, "engine_version": "0.33.1",
//!   "error": { "code": "command_failed", "message": "...", "details": { ... } } }
//! ```
//!
//! Compatibility policy: the envelope itself rejects unknown fields (a caller that
//! misspells `command` should hear about it), while inner `request` payloads follow
//! each command's own schema discipline (e.g. `solve` uses the solve-json v1
//! decoder, which also rejects unknown fields with location info). New commands
//! and new OPTIONAL response fields may appear within api_version 1; anything that
//! would break an existing well-formed caller bumps `BRIDGE_API_VERSION`. Callers
//! feature-detect with `meta.capabilities` instead of sniffing versions.

#[cfg(feature = "ffi")]
pub mod ffi;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::panic::{catch_unwind, AssertUnwindSafe};

/// Bridge envelope version. Bumped only for breaking envelope changes.
pub const BRIDGE_API_VERSION: u32 = 1;

/// Hard cap on request size, matching the solve-json transport.
pub const MAX_REQUEST_BYTES: usize = 1024 * 1024;

const ENGINE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Commands available in this build, in dispatch order.
/// `meta.capabilities` reports exactly this list so apps can feature-detect.
fn command_names() -> Vec<&'static str> {
    vec![
        "meta.capabilities",
        "meta.version",
        "solve",
        "card.come_ups",
        "card.range_table",
        "card.wind",
    ]
}

fn compiled_features() -> Vec<&'static str> {
    [
        ("pdf", cfg!(feature = "pdf")),
        ("profile-import", cfg!(feature = "profile-import")),
        ("online", cfg!(feature = "online")),
    ]
    .iter()
    .filter(|(_, enabled)| *enabled)
    .map(|(name, _)| *name)
    .collect()
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct BridgeRequest {
    api_version: u32,
    command: String,
    #[serde(default)]
    request: Value,
}

/// Machine-readable bridge error codes. Distinct from any command's own error
/// vocabulary: a `command_failed` carries the command's typed error in `details`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BridgeErrorCode {
    InvalidJson,
    UnsupportedApiVersion,
    UnknownCommand,
    InvalidRequest,
    ResourceLimit,
    CommandFailed,
    InternalError,
}

fn success(command: &str, result: Value) -> String {
    serialize_envelope(&json!({
        "ok": true,
        "api_version": BRIDGE_API_VERSION,
        "engine_version": ENGINE_VERSION,
        "command": command,
        "result": result,
    }))
}

fn error(code: BridgeErrorCode, message: impl Into<String>, details: Option<Value>) -> String {
    let mut error = json!({
        "code": code,
        "message": message.into(),
    });
    if let Some(details) = details {
        error["details"] = details;
    }
    serialize_envelope(&json!({
        "ok": false,
        "api_version": BRIDGE_API_VERSION,
        "engine_version": ENGINE_VERSION,
        "error": error,
    }))
}

/// Serialization of the envelope itself must not be able to fail the bridge:
/// fall back to a hand-written internal_error document.
fn serialize_envelope(value: &Value) -> String {
    serde_json::to_string(value).unwrap_or_else(|_| {
        format!(
            r#"{{"ok":false,"api_version":{BRIDGE_API_VERSION},"engine_version":"{ENGINE_VERSION}","error":{{"code":"internal_error","message":"bridge response serialization failed"}}}}"#
        )
    })
}

/// Process one bridge exchange. Never panics; every failure mode is an in-band
/// error envelope. This is the function the C ABI wraps.
pub fn bridge_call(request_json: &str) -> String {
    let guarded = catch_unwind(AssertUnwindSafe(|| dispatch(request_json)));
    guarded.unwrap_or_else(|_| {
        error(
            BridgeErrorCode::InternalError,
            "bridge command failed unexpectedly",
            None,
        )
    })
}

fn dispatch(request_json: &str) -> String {
    if request_json.len() > MAX_REQUEST_BYTES {
        return error(
            BridgeErrorCode::ResourceLimit,
            format!("bridge request exceeds the {MAX_REQUEST_BYTES}-byte limit"),
            None,
        );
    }

    let request: BridgeRequest = match serde_json::from_str(request_json) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidJson,
                format!("bridge request is not a valid envelope: {err}"),
                None,
            )
        }
    };

    if request.api_version != BRIDGE_API_VERSION {
        return error(
            BridgeErrorCode::UnsupportedApiVersion,
            format!(
                "unsupported api_version {}; this build speaks {BRIDGE_API_VERSION}",
                request.api_version
            ),
            None,
        );
    }

    match request.command.as_str() {
        "meta.capabilities" => success(
            "meta.capabilities",
            json!({
                "engine_version": ENGINE_VERSION,
                "bridge_api_version": BRIDGE_API_VERSION,
                "commands": command_names(),
                "features": compiled_features(),
                "solve_schema_version": crate::solve_json::SOLVE_JSON_SCHEMA_VERSION_V1,
            }),
        ),
        "meta.version" => success(
            "meta.version",
            json!({ "engine_version": ENGINE_VERSION }),
        ),
        "solve" => run_solve(&request.request),
        "card.come_ups" => run_card(&request.request, "card.come_ups", crate::card_service::come_ups_v1),
        "card.range_table" => {
            run_card(&request.request, "card.range_table", crate::card_service::range_table_v1)
        }
        "card.wind" => run_card(&request.request, "card.wind", crate::card_service::wind_card_v1),
        other => error(
            BridgeErrorCode::UnknownCommand,
            format!(
                "unknown command '{other}'; this build supports: {}",
                command_names().join(", ")
            ),
            None,
        ),
    }
}

/// `solve` delegates verbatim to the solve-json v1 service. The inner request is
/// re-serialized and run through [`crate::solve_json::decode_solve_request_v1`] so
/// callers get the exact same schema validation (unknown-field rejection, explicit
/// SI units, typed error locations) as the CLI `solve-json` transport.
fn run_solve(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'solve' requires a request payload (solve-json v1 document)",
            None,
        );
    }
    let inner_text = match serde_json::to_string(inner) {
        Ok(text) => text,
        Err(err) => {
            return error(
                BridgeErrorCode::InternalError,
                format!("failed to re-serialize solve request: {err}"),
                None,
            )
        }
    };

    let request = match crate::solve_json::decode_solve_request_v1(&inner_text) {
        Ok(request) => request,
        Err(envelope) => return command_error("solve request rejected", &envelope),
    };

    match crate::solve_v1(request) {
        Ok(successful) => match serde_json::to_value(&successful) {
            Ok(result) => success("solve", result),
            Err(err) => error(
                BridgeErrorCode::InternalError,
                format!("failed to serialize solve result: {err}"),
                None,
            ),
        },
        Err(envelope) => command_error("solve failed", &envelope),
    }
}

/// Shared runner for the three card commands: decode the typed request (rejecting
/// unknown fields), run the transport-free service, serialize the typed response.
fn run_card(
    inner: &Value,
    command: &'static str,
    service: fn(
        &crate::card_service::CardRequestV1,
    ) -> Result<crate::card_service::CardResponseV1, crate::card_service::CardServiceError>,
) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            format!("'{command}' requires a request payload (card v1 document)"),
            None,
        );
    }
    let request: crate::card_service::CardRequestV1 =
        match serde_json::from_value(inner.clone()) {
            Ok(request) => request,
            Err(err) => {
                return error(
                    BridgeErrorCode::InvalidRequest,
                    format!("{command} request rejected: {err}"),
                    None,
                )
            }
        };
    match service(&request) {
        Ok(response) => match serde_json::to_value(&response) {
            Ok(result) => success(command, result),
            Err(err) => error(
                BridgeErrorCode::InternalError,
                format!("failed to serialize {command} result: {err}"),
                None,
            ),
        },
        Err(err) => error(
            BridgeErrorCode::CommandFailed,
            format!("{command} failed: {err}"),
            None,
        ),
    }
}

/// Wrap a command's own typed error envelope losslessly in `details`.
fn command_error<E: Serialize>(message: &str, typed: &E) -> String {
    let details = serde_json::to_value(typed).ok();
    error(BridgeErrorCode::CommandFailed, message, details)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(value: Value) -> Value {
        let raw = bridge_call(&value.to_string());
        serde_json::from_str(&raw).expect("bridge output must be valid JSON")
    }

    #[test]
    fn capabilities_reports_commands_and_versions() {
        let out = call(json!({"api_version": 1, "command": "meta.capabilities"}));
        assert_eq!(out["ok"], true);
        assert_eq!(out["api_version"], 1);
        assert_eq!(out["result"]["engine_version"], ENGINE_VERSION);
        let commands: Vec<String> =
            serde_json::from_value(out["result"]["commands"].clone()).unwrap();
        assert!(commands.contains(&"solve".to_string()));
        assert!(commands.contains(&"meta.capabilities".to_string()));
    }

    #[test]
    fn invalid_json_is_an_envelope_not_a_panic() {
        let out: Value = serde_json::from_str(&bridge_call("{not json")).unwrap();
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_json");
    }

    #[test]
    fn unknown_envelope_field_is_rejected() {
        let out = call(json!({"api_version": 1, "command": "meta.version", "extra": 1}));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_json");
    }

    #[test]
    fn unknown_command_lists_supported_ones() {
        let out = call(json!({"api_version": 1, "command": "card.pdf"}));
        assert_eq!(out["error"]["code"], "unknown_command");
        assert!(out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("meta.capabilities"));
    }

    #[test]
    fn wrong_api_version_is_rejected() {
        let out = call(json!({"api_version": 99, "command": "meta.version"}));
        assert_eq!(out["error"]["code"], "unsupported_api_version");
    }

    #[test]
    fn oversize_request_is_a_resource_limit() {
        let big = format!(
            r#"{{"api_version":1,"command":"meta.version","request":"{}"}}"#,
            "x".repeat(MAX_REQUEST_BYTES)
        );
        let out: Value = serde_json::from_str(&bridge_call(&big)).unwrap();
        assert_eq!(out["error"]["code"], "resource_limit");
    }

    #[test]
    fn solve_without_payload_is_invalid_request() {
        let out = call(json!({"api_version": 1, "command": "solve"}));
        assert_eq!(out["error"]["code"], "invalid_request");
    }

    #[test]
    fn solve_with_bad_schema_carries_typed_details() {
        let out = call(json!({
            "api_version": 1,
            "command": "solve",
            "request": {"schema_version": 1, "unknown_field": true}
        }));
        assert_eq!(out["error"]["code"], "command_failed");
        // The solve-json envelope rides along losslessly.
        assert_eq!(out["error"]["details"]["status"], "error");
    }
}
