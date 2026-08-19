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
    let mut names = vec![
        "meta.capabilities",
        "meta.version",
        "solve",
        "card.come_ups",
        "card.range_table",
        "card.wind",
    ];
    // Listed ONLY when compiled in (mirroring compiled_features) so apps feature-detect
    // the command list instead of probing for unknown_command. Each conditional push sits
    // at its dispatch position so this list stays in dispatch order, as documented above.
    #[cfg(feature = "pdf")]
    names.push("card.pdf");
    names.extend(["profile.validate", "profile.normalize"]);
    #[cfg(feature = "profile-import")]
    names.push("profile.import_a7p");
    // Filesystem-backed (BC5D tables are loaded from caller-supplied paths), so absent on
    // wasm32 — the same "list only what this build can run" rule as profile.import_a7p.
    #[cfg(not(target_arch = "wasm32"))]
    names.push("bc5d.info");
    names
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
        #[cfg(feature = "pdf")]
        "card.pdf" => run_card_pdf(&request.request),
        "profile.validate" => run_profile_validate(&request.request),
        "profile.normalize" => run_profile_normalize(&request.request),
        #[cfg(feature = "profile-import")]
        "profile.import_a7p" => run_profile_import_a7p(&request.request),
        #[cfg(not(target_arch = "wasm32"))]
        "bc5d.info" => run_bc5d_info(&request.request),
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

/// Hard cap on the PDF `card.pdf` will hand back, measured on the RAW document (the
/// base64 text in the response is ~4/3 of it, so this bounds a ~5.6 MiB response body).
///
/// Every dope card carries a ~815 KiB floor: the two Liberation Sans faces
/// `pdf_dope_card` embeds. Rows are cheap on top of that (~0.5 KiB each — a 300-row,
/// 4-page card is ~950 KiB), so this cap leaves room for roughly 7,000 rows / 70 pages and
/// only trips on a request that asked for something absurd: a hair-fine `step` over a long
/// domain, or a header/footer label of hundreds of KiB (the `pdf` block's strings are drawn
/// verbatim on every page). In those cases a typed refusal the caller can act on beats
/// pushing a multi-megabyte base64 string through an embedded FFI hop.
#[cfg(feature = "pdf")]
pub const MAX_PDF_BYTES: usize = 4 * 1024 * 1024;

/// The over-cap envelope for a generated PDF, or `None` when it fits. Split out so the
/// boundary itself is unit-testable at exactly [`MAX_PDF_BYTES`] and one byte past it.
#[cfg(feature = "pdf")]
fn pdf_over_cap_error(byte_length: usize) -> Option<String> {
    (byte_length > MAX_PDF_BYTES).then(|| {
        error(
            BridgeErrorCode::ResourceLimit,
            format!(
                "generated dope card is {byte_length} bytes; the limit is {MAX_PDF_BYTES} \
                 (coarsen step, shorten the range domain, or shorten the pdf labels)"
            ),
            None,
        )
    })
}

/// `card.pdf`: the printable dope card, as base64. The request is the SAME
/// [`crate::card_service::CardRequestV1`] the on-screen card commands take — an app stores
/// one request per saved card and replays it here — with the optional presentation-only
/// `pdf` block for the header/footer labels, font size, and the Lead column's target speed.
///
/// The printed Range/Drop/Wind figures are the rows `card.range_table` returns for that
/// same request, from the same computation (see `card_service::range_table_rows`), so a
/// reprint cannot disagree with the card the shooter already read.
///
/// Result: `{ "pdf_base64": ..., "byte_length": <raw PDF bytes>, "page_count": ... }`.
/// `byte_length` describes the DECODED document, not the base64 text. Documents over
/// [`MAX_PDF_BYTES`] are refused with `resource_limit` instead of returned.
///
/// Errors follow the sibling card commands exactly: a malformed payload is
/// `invalid_request`, anything the service rejects (including an out-of-band
/// `pdf.font_scale`) is `command_failed` with the service's own message.
#[cfg(feature = "pdf")]
fn run_card_pdf(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'card.pdf' requires a request payload (card v1 document)",
            None,
        );
    }
    let request: crate::card_service::CardRequestV1 = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("card.pdf request rejected: {err}"),
                None,
            )
        }
    };
    let card = match crate::card_service::pdf_card_v1(&request) {
        Ok(card) => card,
        Err(err) => {
            return error(
                BridgeErrorCode::CommandFailed,
                format!("card.pdf failed: {err}"),
                None,
            )
        }
    };
    let byte_length = card.pdf_bytes.len();
    if let Some(envelope) = pdf_over_cap_error(byte_length) {
        return envelope;
    }
    success(
        "card.pdf",
        json!({
            "pdf_base64": encode_base64(&card.pdf_bytes),
            "byte_length": byte_length,
            "page_count": card.page_count,
        }),
    )
}

/// RFC 4648 standard-alphabet base64 encoder with padding, for `card.pdf`.
///
/// Hand-rolled for the same reason as `decode_base64` below (plain text, not a doc link:
/// that function is gated on `profile-import`, which a pdf-only build need not enable): no
/// direct base64 dependency
/// exists in `Cargo.toml`, and adding one for twenty lines of arithmetic would ride along on
/// all thirteen release targets.
#[cfg(feature = "pdf")]
fn encode_base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let triple = (u32::from(chunk[0]) << 16)
            | (u32::from(chunk.get(1).copied().unwrap_or(0)) << 8)
            | u32::from(chunk.get(2).copied().unwrap_or(0));
        out.push(char::from(ALPHABET[(triple >> 18) as usize & 63]));
        out.push(char::from(ALPHABET[(triple >> 12) as usize & 63]));
        // A 1- or 2-byte tail pads rather than encoding the zero bits it never carried.
        out.push(if chunk.len() > 1 {
            char::from(ALPHABET[(triple >> 6) as usize & 63])
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            char::from(ALPHABET[triple as usize & 63])
        } else {
            '='
        });
    }
    out
}

/// Wrap a command's own typed error envelope losslessly in `details`.
fn command_error<E: Serialize>(message: &str, typed: &E) -> String {
    let details = serde_json::to_value(typed).ok();
    error(BridgeErrorCode::CommandFailed, message, details)
}

/// Shared decode for the two profile document commands: the inner request IS a
/// [`crate::profile::ProfileData`] JSON document (the exact schema of
/// `~/.ballistics/profiles/*.json` — same field names, same defaults, unknown keys
/// tolerated, exactly as the CLI loads it).
fn decode_profile_document(
    inner: &Value,
    command: &'static str,
) -> Result<crate::profile::ProfileData, String> {
    if inner.is_null() {
        return Err(error(
            BridgeErrorCode::InvalidRequest,
            format!("'{command}' requires a request payload (a ProfileData JSON document)"),
            None,
        ));
    }
    serde_json::from_value(inner.clone()).map_err(|err| {
        error(
            BridgeErrorCode::InvalidRequest,
            format!("{command} request is not a ProfileData document: {err}"),
            None,
        )
    })
}

/// `profile.validate`: parse a ProfileData document and run the cheap invariants the CLI
/// applies when loading/saving a profile (units string, MBA-1358 tracking-CF band,
/// MBA-1348 turret/optic assembly + validation including click-value parse) — see
/// [`crate::profile::ProfileData::validation_warnings`]. No new physics checks. Result:
/// `{ "valid": bool, "warnings": [..], "normalized": <the profile re-serialized by this
/// engine> }` — `valid` is simply `warnings.is_empty()`.
fn run_profile_validate(inner: &Value) -> String {
    let profile = match decode_profile_document(inner, "profile.validate") {
        Ok(profile) => profile,
        Err(envelope) => return envelope,
    };
    let warnings = profile.validation_warnings();
    match serde_json::to_value(&profile) {
        Ok(normalized) => success(
            "profile.validate",
            json!({
                "valid": warnings.is_empty(),
                "warnings": warnings,
                "normalized": normalized,
            }),
        ),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize normalized profile: {err}"),
            None,
        ),
    }
}

/// `profile.normalize`: parse a ProfileData document and hand it back re-serialized by
/// THIS engine — the supported way for an app to round-trip a stored blob through a newer
/// engine version (unknown keys are tolerated on input and dropped on output; defaults
/// fill in; `skip_serializing_if` keys stay absent — the same round-trip a CLI
/// load-then-save performs). Result: `{ "profile": <re-serialized ProfileData> }`.
fn run_profile_normalize(inner: &Value) -> String {
    let profile = match decode_profile_document(inner, "profile.normalize") {
        Ok(profile) => profile,
        Err(envelope) => return envelope,
    };
    match serde_json::to_value(&profile) {
        Ok(normalized) => success("profile.normalize", json!({ "profile": normalized })),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize normalized profile: {err}"),
            None,
        ),
    }
}

/// `bc5d.info` request payload: the filesystem path of a downloaded BC5D table.
#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Bc5dInfoRequest {
    path: String,
}

/// `bc5d.info`: open, parse, and CRC-verify a BC5D correction table at a caller-supplied
/// path — the exact same load-with-verification (`bc_table_5d::path_cache::load_verified`)
/// the `solve`/card `bc5d_table_path` fields use, so "info says valid" and "the solve will
/// accept it" cannot drift apart. Lets an app validate a table right after downloading it.
///
/// Result on success: `{ "valid": true, "crc_ok": true, ... }` with the identifying
/// metadata the header carries (format version, caliber, generator API version,
/// generation timestamp, per-axis bin counts, total cells, weight/velocity coverage).
/// A missing, unreadable, corrupt, or non-BC5D file is a `command_failed` envelope with
/// a human-readable message (`invalid_request` when the payload itself is malformed).
///
/// `caliber` is the raw header value (an `f32`, so a .308 table reports 0.30799998) and
/// `caliber_key` is that same value as the 3-digit BC5D key — EXACTLY the integer the
/// solve/card caliber guard compares (`Bc5dTable::ensure_caliber_matches`). An app can
/// therefore pre-check a downloaded table itself with
/// `round(bullet_diameter_inches * 1000) == caliber_key` and show its own friendly
/// message instead of provoking the `command_failed`. This command deliberately does NOT
/// take a caliber: it describes a file, and only the surfaces that have a shot in hand
/// enforce the match.
#[cfg(not(target_arch = "wasm32"))]
fn run_bc5d_info(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'bc5d.info' requires a request payload ({\"path\": ...})",
            None,
        );
    }
    let request: Bc5dInfoRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("bc5d.info request rejected: {err}"),
                None,
            )
        }
    };

    let table = match crate::bc_table_5d::path_cache::load_verified(std::path::Path::new(
        &request.path,
    )) {
        Ok(table) => table,
        Err(err) => {
            return error(
                BridgeErrorCode::CommandFailed,
                format!("bc5d.info: not a usable BC5D table: {err}"),
                None,
            )
        }
    };

    let (weight, bc, muzzle_vel, current_vel, drag_types) = table.bin_counts();
    let (weight_lo, weight_hi) = table.weight_range();
    let (vel_lo, vel_hi) = table.velocity_range();
    success(
        "bc5d.info",
        json!({
            // Reaching here means the magic, format version, dimensions, AND the stored
            // CRC32 all checked out — crc_ok is not a separate weaker probe.
            "valid": true,
            "crc_ok": true,
            "format_version": table.version(),
            "caliber": table.caliber(),
            // The integer the caliber guard actually compares (see the doc comment).
            "caliber_key": table.caliber_key(),
            "api_version": table.api_version(),
            "generated_timestamp": table.timestamp(),
            // Axis order matches the on-disk layout: [drag_type][weight][bc][mv][cv].
            "bins": {
                "weight": weight,
                "bc": bc,
                "muzzle_velocity": muzzle_vel,
                "current_velocity": current_vel,
                "drag_types": drag_types,
            },
            "total_cells": table.total_cells(),
            "weight_range_grains": [weight_lo, weight_hi],
            "velocity_range_fps": [vel_lo, vel_hi],
        }),
    )
}

/// Hard cap on the DECODED `.a7p` payload accepted by `profile.import_a7p`. Real files
/// are a few KiB; this exists purely as a resource bound (the request envelope's own
/// [`MAX_REQUEST_BYTES`] already caps the base64 text).
#[cfg(feature = "profile-import")]
pub const MAX_A7P_DECODED_BYTES: usize = 1024 * 1024;

/// `profile.import_a7p` request payload. `zero_click` mirrors the CLI's `--zero-click`
/// (the source device's click graduation, e.g. `"0.1mil"`), enabling the same optional
/// zero_x/zero_y click-count conversion; omitted keeps the CLI's default behavior (the
/// counts are reported as unmapped). `strict` mirrors `--strict`: reject the file on an
/// MD5 envelope mismatch instead of importing with a warning.
#[cfg(feature = "profile-import")]
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProfileImportA7pRequest {
    a7p_base64: String,
    #[serde(default)]
    zero_click: Option<String>,
    #[serde(default)]
    strict: bool,
}

/// `profile.import_a7p`: run the cleanroom `.a7p` parser + the SAME
/// [`crate::profile_import::map_a7p_to_profile`] mapping the CLI's `profile import` uses,
/// on a base64-supplied file. Result: `{ "profile": <ProfileData>, "warnings": [..],
/// "mapped": [[source, raw, converted, destination], ..], "unmapped": [[field, why], ..],
/// "unknown_fields": [{context, number}, ..] }` — the full import report, nothing
/// silently dropped (`unmapped` includes the unknown-field entries too, exactly as the
/// CLI prints them; `unknown_fields` additionally lists the parser-level unknowns in
/// structured form). The profile name is derived from the file (sanitized); renaming is
/// the caller's business — there is no name override here.
#[cfg(feature = "profile-import")]
fn run_profile_import_a7p(inner: &Value) -> String {
    use crate::profile_import::{map_a7p_to_profile, parse_a7p, EnvelopeStatus};

    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'profile.import_a7p' requires a request payload ({\"a7p_base64\": ...})",
            None,
        );
    }
    let request: ProfileImportA7pRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("profile.import_a7p request rejected: {err}"),
                None,
            )
        }
    };

    let zero_click = match request.zero_click.as_deref() {
        Some(raw) => match crate::adjustment::parse_click_value(raw) {
            Ok(click) => Some(click),
            Err(err) => {
                return error(
                    BridgeErrorCode::InvalidRequest,
                    format!("profile.import_a7p zero_click: {err}"),
                    None,
                )
            }
        },
        None => None,
    };

    let bytes = match decode_base64(&request.a7p_base64) {
        Ok(bytes) => bytes,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("profile.import_a7p a7p_base64: {err}"),
                None,
            )
        }
    };
    if bytes.len() > MAX_A7P_DECODED_BYTES {
        return error(
            BridgeErrorCode::ResourceLimit,
            format!(
                "decoded .a7p payload is {} bytes; the limit is {MAX_A7P_DECODED_BYTES}",
                bytes.len()
            ),
            None,
        );
    }

    let doc = match parse_a7p(&bytes) {
        Ok(doc) => doc,
        Err(err) => {
            return error(
                BridgeErrorCode::CommandFailed,
                format!("not a usable .a7p file: {err}"),
                None,
            )
        }
    };
    // Same refusal (and message) as the CLI's --strict; without it the mismatch becomes
    // a warning in the report, also exactly as the CLI behaves.
    if request.strict {
        if let EnvelopeStatus::Mismatch { expected, actual } = &doc.envelope {
            return error(
                BridgeErrorCode::CommandFailed,
                format!(
                    "checksum mismatch (file says {expected}, payload hashes to {actual}) — refusing under strict"
                ),
                None,
            );
        }
    }

    let outcome = match map_a7p_to_profile(&doc, None, zero_click) {
        Ok(outcome) => outcome,
        Err(err) => return error(BridgeErrorCode::CommandFailed, err, None),
    };
    let unknown_fields: Vec<Value> = doc
        .unknown_fields
        .iter()
        .map(|u| json!({ "context": u.context, "number": u.number }))
        .collect();
    match serde_json::to_value(&outcome.profile) {
        Ok(profile) => success(
            "profile.import_a7p",
            json!({
                "profile": profile,
                "warnings": outcome.report.warnings,
                "mapped": outcome.report.mapped,
                "unmapped": outcome.report.unmapped,
                "unknown_fields": unknown_fields,
            }),
        ),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize imported profile: {err}"),
            None,
        ),
    }
}

/// Minimal strict RFC 4648 standard-alphabet base64 decoder for `profile.import_a7p`.
///
/// Hand-rolled rather than a new dependency, deliberately: the crate already carries its
/// own cleanroom MD5 (`profile_import::md5`) and statistical constants for the same
/// thirteen-platform reasons, `Cargo.toml` has no direct base64 dependency today, and the
/// input here is a few KiB. Strict: rejects any character outside `A-Za-z0-9+/`, `=`
/// anywhere but as final padding, and lengths of form 4n+1.
#[cfg(feature = "profile-import")]
fn decode_base64(input: &str) -> Result<Vec<u8>, String> {
    fn sextet(c: u8) -> Result<u32, String> {
        match c {
            b'A'..=b'Z' => Ok(u32::from(c - b'A')),
            b'a'..=b'z' => Ok(u32::from(c - b'a') + 26),
            b'0'..=b'9' => Ok(u32::from(c - b'0') + 52),
            b'+' => Ok(62),
            b'/' => Ok(63),
            _ => Err(format!("invalid base64 character {:?}", char::from(c))),
        }
    }
    let bytes = input.as_bytes();
    let data = match bytes {
        [rest @ .., b'=', b'='] => rest,
        [rest @ .., b'='] => rest,
        _ => bytes,
    };
    if data.contains(&b'=') {
        return Err("'=' is only valid as trailing padding".to_string());
    }
    if data.len() % 4 == 1 {
        return Err("base64 text has an impossible length (4n+1 data characters)".to_string());
    }
    let mut out = Vec::with_capacity(data.len() / 4 * 3 + 2);
    let mut acc: u32 = 0;
    let mut bits: u32 = 0;
    for &c in data {
        acc = (acc << 6) | sextet(c)?;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((acc >> bits) as u8);
        }
    }
    Ok(out)
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
        // Was `card.pdf` until that became a real (pdf-gated) command; use a name no build
        // can ever dispatch so this test means the same thing in every feature set.
        let out = call(json!({"api_version": 1, "command": "card.semaphore"}));
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
    fn profile_commands_without_payload_are_invalid_requests() {
        for command in ["profile.validate", "profile.normalize"] {
            let out = call(json!({"api_version": 1, "command": command}));
            assert_eq!(out["error"]["code"], "invalid_request", "{command}: {out}");
            assert!(
                out["error"]["message"]
                    .as_str()
                    .unwrap()
                    .contains("ProfileData"),
                "{command}: {out}"
            );
        }
    }

    #[test]
    fn capabilities_lists_profile_commands_and_gates_import_on_the_feature() {
        let out = call(json!({"api_version": 1, "command": "meta.capabilities"}));
        let commands: Vec<String> =
            serde_json::from_value(out["result"]["commands"].clone()).unwrap();
        assert!(commands.contains(&"profile.validate".to_string()));
        assert!(commands.contains(&"profile.normalize".to_string()));
        assert_eq!(
            commands.contains(&"profile.import_a7p".to_string()),
            cfg!(feature = "profile-import"),
            "profile.import_a7p must be listed exactly when compiled in"
        );
        assert_eq!(
            commands.contains(&"bc5d.info".to_string()),
            cfg!(not(target_arch = "wasm32")),
            "bc5d.info must be listed exactly when the build has filesystem access"
        );
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn bc5d_info_without_payload_or_with_missing_file_fails_cleanly() {
        let out = call(json!({"api_version": 1, "command": "bc5d.info"}));
        assert_eq!(out["error"]["code"], "invalid_request", "{out}");

        let out = call(json!({
            "api_version": 1,
            "command": "bc5d.info",
            "request": {"path": "/nonexistent/bc5d_308.bin"}
        }));
        assert_eq!(out["error"]["code"], "command_failed", "{out}");
        assert!(
            out["error"]["message"]
                .as_str()
                .unwrap()
                .contains("not a usable BC5D table"),
            "{out}"
        );
    }

    #[cfg(feature = "profile-import")]
    #[test]
    fn base64_decoder_round_trips_and_rejects_garbage() {
        // RFC 4648 test vectors.
        for (text, bytes) in [
            ("", &b""[..]),
            ("Zg==", b"f"),
            ("Zm8=", b"fo"),
            ("Zm9v", b"foo"),
            ("Zm9vYg==", b"foob"),
            ("Zm9vYmE=", b"fooba"),
            ("Zm9vYmFy", b"foobar"),
        ] {
            assert_eq!(decode_base64(text).unwrap(), bytes, "{text}");
        }
        assert!(decode_base64("Zm9v\n").is_err(), "whitespace is rejected");
        assert!(decode_base64("Zg=X").is_err(), "inner padding is rejected");
        assert!(decode_base64("Z").is_err(), "4n+1 length is rejected");
        assert!(decode_base64("Zm9v!").is_err(), "non-alphabet byte is rejected");
    }

    /// `card.pdf` must be listed exactly when the `pdf` feature is compiled in, and be an
    /// honest `unknown_command` otherwise — the same rule `profile.import_a7p` follows. The
    /// pdf-absent half of this only runs under `--no-default-features --features bridge`.
    #[test]
    fn capabilities_gates_card_pdf_on_the_pdf_feature() {
        let out = call(json!({"api_version": 1, "command": "meta.capabilities"}));
        let commands: Vec<String> =
            serde_json::from_value(out["result"]["commands"].clone()).unwrap();
        assert_eq!(
            commands.contains(&"card.pdf".to_string()),
            cfg!(feature = "pdf"),
            "card.pdf must be listed exactly when compiled in: {out}"
        );
        let features: Vec<String> =
            serde_json::from_value(out["result"]["features"].clone()).unwrap();
        assert_eq!(
            features.contains(&"pdf".to_string()),
            cfg!(feature = "pdf"),
            "the command list and the feature list must agree: {out}"
        );
    }

    #[cfg(not(feature = "pdf"))]
    #[test]
    fn card_pdf_is_an_unknown_command_without_the_pdf_feature() {
        let out = call(json!({
            "api_version": 1,
            "command": "card.pdf",
            "request": {
                "muzzle_velocity": 2600.0, "ballistic_coefficient": 0.243,
                "mass": 175.0, "diameter": 0.308,
                "zero_distance": 100.0, "start": 100.0, "end": 300.0, "step": 100.0
            }
        }));
        assert_eq!(out["error"]["code"], "unknown_command", "{out}");
    }

    /// The `pdf` presentation block must survive a build that cannot render it: an app
    /// stores one request per card and replays it against `card.range_table` too, so a
    /// pdf-less engine has to ACCEPT the field rather than reject it as unknown.
    #[test]
    fn the_pdf_presentation_block_is_accepted_by_the_on_screen_card_in_every_build() {
        let out = call(json!({
            "api_version": 1,
            "command": "card.range_table",
            "request": {
                "muzzle_velocity": 2600.0, "ballistic_coefficient": 0.243,
                "mass": 175.0, "diameter": 0.308,
                "zero_distance": 100.0, "start": 100.0, "end": 300.0, "step": 100.0,
                "pdf": {"title": "Stored Card", "target_speed": 8.0, "font_preset": "large"}
            }
        }));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["kind"], "range_table", "{out}");
    }

    #[cfg(feature = "pdf")]
    #[test]
    fn card_pdf_without_payload_is_invalid_request() {
        let out = call(json!({"api_version": 1, "command": "card.pdf"}));
        assert_eq!(out["error"]["code"], "invalid_request", "{out}");
        assert!(
            out["error"]["message"].as_str().unwrap().contains("card v1 document"),
            "{out}"
        );
    }

    /// The output cap's boundary, both sides. Generating a genuinely over-cap dope card
    /// would take tens of thousands of rows, so the predicate is tested directly — see
    /// `pdf_over_cap_error`'s own comment.
    #[cfg(feature = "pdf")]
    #[test]
    fn pdf_output_cap_refuses_only_over_the_limit() {
        assert!(pdf_over_cap_error(0).is_none());
        assert!(
            pdf_over_cap_error(MAX_PDF_BYTES).is_none(),
            "a document exactly at the cap fits"
        );
        let envelope: Value =
            serde_json::from_str(&pdf_over_cap_error(MAX_PDF_BYTES + 1).expect("over cap"))
                .unwrap();
        assert_eq!(envelope["ok"], false);
        assert_eq!(envelope["error"]["code"], "resource_limit");
        assert!(
            envelope["error"]["message"].as_str().unwrap().contains("dope card"),
            "{envelope}"
        );
    }

    #[cfg(feature = "pdf")]
    #[test]
    fn base64_encoder_matches_the_rfc_4648_vectors() {
        for (bytes, text) in [
            (&b""[..], ""),
            (b"f", "Zg=="),
            (b"fo", "Zm8="),
            (b"foo", "Zm9v"),
            (b"foob", "Zm9vYg=="),
            (b"fooba", "Zm9vYmE="),
            (b"foobar", "Zm9vYmFy"),
        ] {
            assert_eq!(encode_base64(bytes), text, "{bytes:?}");
        }
        // Full-byte-range coverage: the >> 18 / >> 12 / >> 6 masking must not sign- or
        // width-mangle a high byte, which is most of a PDF's content.
        assert_eq!(encode_base64(&[0xff, 0xff, 0xff]), "////");
        assert_eq!(encode_base64(&[0x00, 0x00, 0x00]), "AAAA");
        assert_eq!(encode_base64(&[0xfb, 0xff, 0xbf]), "+/+/");
    }

    /// The encoder and the (profile-import) decoder must be inverses — the property that
    /// makes `pdf_base64` a lossless transport for a binary document.
    #[cfg(all(feature = "pdf", feature = "profile-import"))]
    #[test]
    fn base64_encode_decode_round_trips_arbitrary_bytes() {
        for len in 0..=32usize {
            let bytes: Vec<u8> = (0..len).map(|i| (i as u8).wrapping_mul(37).wrapping_add(11)).collect();
            let decoded = decode_base64(&encode_base64(&bytes)).expect("own output decodes");
            assert_eq!(decoded, bytes, "len {len}");
        }
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
