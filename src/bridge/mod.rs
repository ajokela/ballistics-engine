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
//!
//! ## `.a7p` interop: `profile.import_a7p` and `profile.export_a7p`
//!
//! Both directions of the ArcherBC2 format, and both are LOSSY in a way the caller is
//! required to surface. Import reports `unmapped` — what the file held that a
//! `ProfileData` cannot. Export (MBA-1556) reports `not_carried` — every `ProfileData`
//! field the file cannot hold, listed unconditionally with a `populated` flag, plus the
//! pre-filtered `dropped_fields` for the ones this profile really did lose.
//!
//! The asymmetry is deliberate. On import the user is gaining a profile and can see what
//! arrived; on export they are handing a rifle to somebody else and cannot see what
//! arrives at the far end. An app that renders "exported" without rendering the drop list
//! has told a shooter something untrue, so the list travels in the result rather than
//! being something a caller has to ask for.
//!
//! Export never extends the format to make room: `.a7p` belongs to somebody else, and
//! smuggling our fields into unused field numbers would produce files Archer's own tools
//! misread.
//!
//! ## The `reticle.*` family
//!
//! `reticle.describe`, `reticle.hold` and `reticle.import` put the engine's reticle stack
//! (MBA-1361, MBA-1440, MBA-1544) in reach of an app. Before MBA-1558 none of it was: the
//! bridge named no reticle command, and since the apps vendor only `ballistics_bridge_call`
//! and friends, a capability the bridge does not name does not exist for them.
//!
//! - `reticle.hold` is the one the family exists for — given a firing solution ALREADY
//!   reduced to angles and the optic's current magnification, it reports which mark to hold
//!   on. It runs no physics: the angles are inputs, and `crate::reticle` keeps its
//!   no-physics property precisely because nothing here re-derives them.
//! - `reticle.describe` resolves a reticle, supplied inline or built from a generator, and
//!   optionally reports every mark's TRUE angular position at a magnification. This is what
//!   a picker and a drawing are built on. It carries `focal_plane` and
//!   `magnification_dependent` so a UI knows whether a magnification control changes
//!   anything.
//! - `reticle.import` reads the two third-party formats the crate already parses (Ventum
//!   JSON, `.reticle` XML), each with its report. Both are TEXT and travel inline, so
//!   unlike `profile.import_a7p` there is no base64 step.
//! - `reticle.catalog` lists the NAMED reticles this build can produce
//!   (`crate::reticle_catalog`, MBA-1545), which is what an app picker is populated from —
//!   hardcoding ids in an app means a picker that offers what a pinned engine refuses. Each
//!   entry carries its `source`, and that provenance reaches the wire on purpose: a wrong
//!   subtension is a hold wrong by a mark and looks entirely normal, so what makes it
//!   traceable belongs where the person reading the reticle can see it.
//!
//! Two shapes are shared across the family and are worth stating once. A reticle is named
//! by EXACTLY ONE of `reticle` (a full description), `generator`, or `catalog` (an id from
//! `reticle.catalog`) — a request naming more than one is refused rather than resolved by
//! precedence, because supplying more than one is a caller bug and picking a winner hides
//! it. An unknown catalog id is refused for the same reason and never substituted: a
//! shooter handed a different reticle is being shown holds for glass they are not looking
//! through. And because every generator returns FFP,
//! `focal_plane` / `reference_magnification` may be supplied ALONGSIDE a generator to make
//! an SFP reticle; the same keys beside a full `reticle` are refused, since that
//! description already carries its own.
//!
//! Errors follow `true.dsf`'s convention: `error.code` stays `command_failed` and a stable
//! `reason` rides in `error.details` (one per `ReticleError` variant) with the offending
//! numbers beside it. Additive within api_version 1, and listed by `meta.capabilities`.
//!
//! ## The `true.*` truing family
//!
//! `true.fit`, `true.wind`, `true.tall_target`, `true.dsf`, `true.plan`, and `true.dial_plan`
//! expose the engine's truing methods. All six are unconditional (no filesystem access, so
//! all six are present on wasm32), but they are not otherwise uniform:
//!
//! - `true.fit` (joint MV+BC truing) is backed by the uncertainty solver, so its result
//!   always carries `approximation` — a required enum that is either `Available` with
//!   intervals for both muzzle velocity and BC, or `Unavailable` with a reason, never simply
//!   absent. There is deliberately no command that returns a bare truing point estimate.
//! - `true.wind` (effective crosswind from an observed miss) is the one exception to that
//!   guarantee: `solve_wind_truing` has no uncertainty model, so its result is a bare point
//!   value with no interval. Callers must not present it with `true.fit`'s confidence.
//!   `true.wind` is also the one command whose wire shape is SI throughout (`range_m`,
//!   `miss_right_m`, `sigma_m`) while every other command, including the rest of this
//!   family, is imperial; apps convert at the boundary for `true.wind` specifically.
//! - `true.tall_target` returns a scope's tracking correction factor from a tall-target
//!   test.
//! - `true.dsf` derives a single Mach-keyed drop-scale-factor point from an observed
//!   transonic drop; it never persists into a profile's DSF table, which is the caller's
//!   job. It established this family's structured-`error.details` convention: a
//!   machine-readable `reason` (`invalid_input`, `supersonic`, `out_of_range`,
//!   `degenerate_drop`, `forward_model`) alongside the message; `error.code` stays
//!   `command_failed` for all commands, so existing callers are unaffected.
//! - `true.plan` recommends which candidate ranges to shoot for a joint MV/BC truing
//!   experiment (`crate::truing_plan::plan_truing_experiment_v1`, wired directly — no new
//!   service function). Its error also carries structured `error.details`, under the same
//!   `reason` key as every other command in this family (`invalid_request`,
//!   `insufficient_reachable_candidates`, `no_feasible_design`) plus the `rejected_candidates`
//!   diagnostics the typed error itself carries.
//! - `true.dial_plan` turns a TRUE angular correction into ranked dial/hold/hybrid
//!   execution plans for an INLINE optic (`crate::truing_service::dial_plan_v1`, wrapping
//!   `crate::optic::plan_corrections`). Unlike the CLI's `dial-plan --profile` mode, there
//!   is no profile-loading path here — the optic is supplied inline in the request, since a
//!   saved-profile filesystem read must not enter this bridge. Its error also carries
//!   structured `error.details`: a stable `reason` per `OpticError` variant.
//!
//! None of this needed a `BRIDGE_API_VERSION` bump: the six commands are additive within
//! api_version 1, and `meta.capabilities` lists all six for feature detection.

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
    #[cfg(feature = "profile-export")]
    names.push("profile.export_a7p");
    // Unconditional: the reticle stack is pure geometry with no filesystem access, so all
    // three are present on wasm32 like the `true.*` family.
    names.extend([
        "reticle.catalog",
        "reticle.describe",
        "reticle.hold",
        "reticle.holds",
        "reticle.import",
    ]);
    names.extend([
        "true.fit",
        "true.wind",
        "true.tall_target",
        "true.dsf",
        "true.plan",
        "true.dial_plan",
    ]);
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
        ("profile-export", cfg!(feature = "profile-export")),
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
        "meta.version" => success("meta.version", json!({ "engine_version": ENGINE_VERSION })),
        "solve" => run_solve(&request.request),
        "card.come_ups" => run_service(
            &request.request,
            "card.come_ups",
            crate::card_service::come_ups_v1,
        ),
        "card.range_table" => run_service(
            &request.request,
            "card.range_table",
            crate::card_service::range_table_v1,
        ),
        "card.wind" => run_service(
            &request.request,
            "card.wind",
            crate::card_service::wind_card_v1,
        ),
        #[cfg(feature = "pdf")]
        "card.pdf" => run_card_pdf(&request.request),
        "profile.validate" => run_profile_validate(&request.request),
        "profile.normalize" => run_profile_normalize(&request.request),
        #[cfg(feature = "profile-import")]
        "profile.import_a7p" => run_profile_import_a7p(&request.request),
        #[cfg(feature = "profile-export")]
        "profile.export_a7p" => run_profile_export_a7p(&request.request),
        "reticle.catalog" => run_reticle_catalog(),
        "reticle.describe" => run_reticle_describe(&request.request),
        "reticle.hold" => run_reticle_hold(&request.request),
        "reticle.holds" => run_reticle_holds(&request.request),
        "reticle.import" => run_reticle_import(&request.request),
        "true.fit" => run_service(
            &request.request,
            "true.fit",
            crate::truing_uncertainty::run_uncertainty_truing_v1,
        ),
        "true.wind" => run_service(
            &request.request,
            "true.wind",
            crate::truing_wind::solve_wind_truing,
        ),
        "true.tall_target" => run_service(
            &request.request,
            "true.tall_target",
            crate::truing_service::tall_target_v1,
        ),
        "true.dsf" => run_service_detailed(
            &request.request,
            "true.dsf",
            crate::truing_service::derive_dsf_point_v1,
            crate::truing_service::DsfServiceErrorV1::failure_details,
        ),
        "true.plan" => run_service_detailed(
            &request.request,
            "true.plan",
            crate::truing_plan::plan_truing_experiment_v1,
            crate::truing_plan::TruingPlanErrorV1::failure_details,
        ),
        "true.dial_plan" => run_service_detailed(
            &request.request,
            "true.dial_plan",
            crate::truing_service::dial_plan_v1,
            crate::optic::OpticError::failure_details,
        ),
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

/// Shared adapter for every command backed by a transport-free service: null-check the
/// payload, deserialize the request, call the service, serialize the response. The error
/// mapping is fixed so all commands report failures identically.
fn run_service<Req, Resp, E, F>(inner: &Value, command: &'static str, service: F) -> String
where
    Req: serde::de::DeserializeOwned,
    Resp: serde::Serialize,
    E: std::fmt::Display,
    F: FnOnce(&Req) -> Result<Resp, E>,
{
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            format!("'{command}' requires a request payload"),
            None,
        );
    }
    let request: Req = match serde_json::from_value(inner.clone()) {
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

/// [`run_service`] for a service whose error carries a machine-readable reason.
///
/// `details` lands in `error.details` so a wizard can branch on "supersonic" or
/// "out_of_range" instead of pattern-matching prose. `code` stays `command_failed`, so
/// existing callers are unaffected.
///
/// Callers: `true.dsf`, `true.plan`, and `true.dial_plan`, each of whose service error
/// carries a machine-readable reason worth surfacing in `error.details`.
fn run_service_detailed<Req, Resp, E, F, D>(
    inner: &Value,
    command: &'static str,
    service: F,
    details: D,
) -> String
where
    Req: serde::de::DeserializeOwned,
    Resp: serde::Serialize,
    E: std::fmt::Display,
    F: FnOnce(&Req) -> Result<Resp, E>,
    D: FnOnce(&E) -> Option<Value>,
{
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            format!("'{command}' requires a request payload"),
            None,
        );
    }
    let request: Req = match serde_json::from_value(inner.clone()) {
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
        Err(err) => {
            let d = details(&err);
            error(
                BridgeErrorCode::CommandFailed,
                format!("{command} failed: {err}"),
                d,
            )
        }
    }
}

/// Hard cap on the PDF `card.pdf` will hand back, measured on the RAW document (the
/// base64 text in the response is ~4/3 of it, so this bounds a ~5.6 MiB response body).
///
/// Every dope card carries a ~815 KiB floor: the two Liberation Sans faces
/// `pdf_dope_card` embeds. Rows are cheap on top of that (~0.5 KiB each — a 300-row,
/// 4-page card is ~950 KiB).
///
/// This is the BACKSTOP, not the first line: the row set is refused on its own row and page
/// count before any document exists (`card_service::MAX_PDF_ROWS` / `MAX_PDF_PAGES`), because
/// measuring bytes means having already built and paginated them. What survives that check
/// and still lands here is a card made huge by its LABELS — the `pdf` block's strings are
/// drawn verbatim on every page — and for those a typed refusal the caller can act on beats
/// pushing a multi-megabyte base64 string through an embedded FFI hop.
#[cfg(feature = "pdf")]
pub const MAX_PDF_BYTES: usize = 4 * 1024 * 1024;

/// The over-cap envelope for a generated PDF, or `None` when it fits. Split out so the
/// boundary itself is unit-testable at exactly [`MAX_PDF_BYTES`] and one byte past it.
///
/// The message states what is true of the document — its size, and how many rows and pages
/// it holds. It deliberately does NOT advise "coarsen the step" or "shorten the range
/// domain": for a saved card those are immutable (there is no editor for a snapshot's
/// domain), so naming them told the one user who ever sees this message to do something
/// impossible.
#[cfg(feature = "pdf")]
fn pdf_over_cap_error(byte_length: usize, row_count: usize, page_count: usize) -> Option<String> {
    (byte_length > MAX_PDF_BYTES).then(|| {
        error(
            BridgeErrorCode::ResourceLimit,
            format!(
                "generated dope card is {byte_length} bytes; the limit is {MAX_PDF_BYTES} \
                 ({row_count} rows, {page_count} pages)"
            ),
            None,
        )
    })
}

/// The one `card.pdf`-only key on the request: the rows to print, instead of solving. Not a
/// field on [`crate::card_service::CardRequestV1`], because it is not part of a saved card —
/// it is the card's stored RESPONSE, attached at export time — and a stored request must stay
/// replayable against `card.range_table` unchanged. Removed from the payload before the card
/// request is decoded, so `deny_unknown_fields` still governs everything else (a
/// `stored_cards` typo is an honest `invalid_request`).
#[cfg(feature = "pdf")]
const STORED_CARD_KEY: &str = "stored_card";

/// `card.pdf`: the printable dope card, as base64. The request is the SAME
/// [`crate::card_service::CardRequestV1`] the on-screen card commands take — an app stores
/// one request per saved card and replays it here — with the optional presentation-only
/// `pdf` block for the header/footer labels, font size, and the Lead column's target speed,
/// plus one `card.pdf`-only key:
///
/// * `stored_card` (optional): `{ "card": <a stored card.range_table result, verbatim>,
///   "engine_version": "0.34.1", "bc5d_table_version": "2.5.0" }`. Supply it and this
///   command PRINTS THOSE ROWS: no zero solve, no trajectory, and `bc5d_table_path` is never
///   opened, so a saved card reprints identically after an engine bump, after the correction
///   table at that path is overwritten in place, and even after it is deleted. The footer's
///   `BC:` is the stored card's own `bc_for_solve`, and its `Engine:`/`Table:` are the two
///   provenance strings, so paper and screen can be reconciled afterwards.
/// * Omit it (or send `null`, which means the same thing) and the rows are solved here, from
///   the same `card_service::range_table_rows` call `card.range_table` makes — the
///   pre-existing behaviour, unchanged.
///
/// This surface prints a range-table card and says so. A request carrying a wind card's
/// `wind_speeds`/`wind_angles_deg`, or a `stored_card` of another kind, is REFUSED: an `ok`
/// response whose defining field was silently ignored is worse than no response.
///
/// Result: `{ "pdf_base64": ..., "byte_length": <raw PDF bytes>, "page_count": ...,
/// "row_count": ..., "kind": "range_table", "source": "solve" | "stored_rows",
/// "unprintable_title_chars": "" }`.
/// `byte_length` describes the DECODED document, not the base64 text; `source` lets a caller
/// verify it got a reprint rather than a re-solve. `unprintable_title_chars` is normally
/// empty and names the characters of `pdf.title` the card font could not draw when it is not
/// — the card still printed, with a visible stand-in for each of them, but a caller that
/// accepts any card name should warn rather than hand over an untitled card. A card too big
/// to print is refused with
/// `resource_limit` — on its row/page count first (`card_service::MAX_PDF_ROWS` /
/// `MAX_PDF_PAGES`), and on [`MAX_PDF_BYTES`] as the backstop.
///
/// Other errors follow the sibling card commands exactly: a malformed payload is
/// `invalid_request`, anything the service rejects (including an out-of-band
/// `pdf.font_scale`) is `command_failed` with the service's own message.
#[cfg(feature = "pdf")]
fn run_card_pdf(inner: &Value) -> String {
    use crate::card_service::CardServiceError;

    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'card.pdf' requires a request payload (card v1 document)",
            None,
        );
    }
    let mut payload = inner.clone();
    let stored_value = payload
        .as_object_mut()
        .and_then(|object| object.remove(STORED_CARD_KEY))
        .filter(|value| !value.is_null());
    let request: crate::card_service::CardRequestV1 = match serde_json::from_value(payload) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("card.pdf request rejected: {err}"),
                None,
            )
        }
    };
    let stored: Option<crate::card_service::StoredCardV1> = match stored_value {
        Some(value) => match serde_json::from_value(value) {
            Ok(stored) => Some(stored),
            Err(err) => {
                return error(
                    BridgeErrorCode::InvalidRequest,
                    format!("card.pdf {STORED_CARD_KEY} rejected: {err}"),
                    None,
                )
            }
        },
        None => None,
    };

    let card = match crate::card_service::pdf_card_v1(&request, stored.as_ref()) {
        Ok(card) => card,
        // A card too large to print is a resource refusal, not a command failure: same code
        // the byte cap below reports, so a caller has one condition to handle.
        Err(err @ CardServiceError::TooLarge(_)) => {
            return error(
                BridgeErrorCode::ResourceLimit,
                format!("card.pdf refused: {err}"),
                None,
            )
        }
        Err(err) => {
            return error(
                BridgeErrorCode::CommandFailed,
                format!("card.pdf failed: {err}"),
                None,
            )
        }
    };
    let byte_length = card.pdf_bytes.len();
    if let Some(envelope) = pdf_over_cap_error(byte_length, card.row_count, card.page_count) {
        return envelope;
    }
    let mut response = json!({
        "pdf_base64": encode_base64(&card.pdf_bytes),
        "byte_length": byte_length,
        "page_count": card.page_count,
        "row_count": card.row_count,
        "kind": crate::card_service::PDF_CARD_KIND,
        "source": card.source.as_str(),
        "unprintable_title_chars": card.unprintable_title_chars,
    });
    // MBA-1477: additive, and spelled exactly as `card.range_table`'s own block, so an app
    // that already reads a truncation off the on-screen card reads the printed one the same
    // way. Absent on a card that reached every row it asked for, which is every card the
    // previous shape could describe at all — a truncated document used to come back here
    // indistinguishable from a complete one, with only a smaller `row_count` to notice.
    if let (Some(truncation), Some(object)) = (card.truncation, response.as_object_mut()) {
        object.insert(
            "truncation".to_string(),
            json!({
                "requested_end": truncation.requested_end,
                "last_row": truncation.last_row,
                "reach": truncation.reach,
            }),
        );
    }
    success("card.pdf", response)
}

/// RFC 4648 standard-alphabet base64 encoder with padding, for `card.pdf` and
/// `profile.export_a7p` — hence the two-feature gate, which must list every feature that
/// has a binary to hand out or the function goes missing from exactly the build that
/// needs it.
///
/// Hand-rolled for the same reason as `decode_base64` below (plain text, not a doc link:
/// that function is gated on `profile-import`, which a pdf-only build need not enable): no
/// direct base64 dependency
/// exists in `Cargo.toml`, and adding one for twenty lines of arithmetic would ride along on
/// all thirteen release targets.
#[cfg(any(feature = "pdf", feature = "profile-export"))]
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

    let table =
        match crate::bc_table_5d::path_cache::load_verified(std::path::Path::new(&request.path)) {
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

/// `profile.export_a7p`: encode a ProfileData document as an ArcherBC2 `.a7p` file,
/// with the full account of what could not travel. Request payload is the ProfileData
/// document itself, exactly like `profile.validate`/`profile.normalize`.
///
/// Result: `{ "a7p_base64": "..", "byte_length": N, "carried_fields": [..],
/// "not_carried": [{field, populated, value, reason}, ..], "dropped_fields": [..],
/// "warnings": [..] }`.
///
/// `not_carried` is the reason this command has a shape of its own rather than
/// returning bare bytes, and it is NOT optional output. `.a7p` is a third-party format
/// built around one rifle, one load and one device; a saved profile holds a good deal it
/// has no slot for, and an app that shows a shooter "exported" without showing them what
/// their friend will not receive has told them something untrue. `not_carried` lists
/// every unsupported field unconditionally, with `populated` saying whether THIS profile
/// actually lost anything there, and `dropped_fields` is the pre-filtered subset that
/// did — the list to put in front of a person. `carried_fields` is the other half of the
/// same partition, so a caller can account for the whole document.
///
/// Errors carry a structured `reason` in `error.details` under the `true.*` convention
/// (`error.code` stays `command_failed`): `units`, `drag_model`, or `field` plus the
/// offending field name.
#[cfg(feature = "profile-export")]
fn run_profile_export_a7p(inner: &Value) -> String {
    use crate::profile_export::{export_a7p, A7pExportError, CARRIED_FIELDS};

    let profile = match decode_profile_document(inner, "profile.export_a7p") {
        Ok(profile) => profile,
        Err(envelope) => return envelope,
    };
    let export = match export_a7p(&profile) {
        Ok(export) => export,
        Err(err) => {
            let details = match &err {
                A7pExportError::Units(_) => json!({ "reason": "units" }),
                A7pExportError::DragModel(model) => {
                    json!({ "reason": "drag_model", "drag_model": model })
                }
                A7pExportError::Field { field, .. } => {
                    json!({ "reason": "field", "field": field })
                }
            };
            return error(
                BridgeErrorCode::CommandFailed,
                err.to_string(),
                Some(details),
            );
        }
    };
    let not_carried = match serde_json::to_value(&export.not_carried) {
        Ok(value) => value,
        Err(err) => {
            return error(
                BridgeErrorCode::InternalError,
                format!("failed to serialize the not-carried list: {err}"),
                None,
            )
        }
    };
    success(
        "profile.export_a7p",
        json!({
            "a7p_base64": encode_base64(&export.bytes),
            "byte_length": export.bytes.len(),
            "carried_fields": CARRIED_FIELDS,
            "not_carried": not_carried,
            "dropped_fields": export.dropped_fields(),
            "warnings": export.warnings,
        }),
    )
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

// ---------------------------------------------------------------------------
// reticle.* — MBA-1558
// ---------------------------------------------------------------------------

/// How a request names the reticle it is asking about: either the full
/// [`crate::reticle::ReticleDescription`] inline, or a generator to build one from.
///
/// Exactly one, and the requests below check that rather than silently preferring one —
/// a caller that sends both has a bug, and picking a winner for them hides it.
#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum ReticleGenerator {
    /// [`crate::reticle::ReticleDescription::mil_grid`].
    MilGrid { spacing_mil: f64, extent_mil: f64 },
    /// [`crate::reticle::ReticleDescription::tree`].
    Tree {
        rows: usize,
        row_spacing_mil: f64,
        spread_step_mil: f64,
    },
    /// [`crate::reticle::ReticleDescription::bdc_from_drops`]. Pairs are
    /// `[range_metres, drop_mil]`, and the generator labels each mark `"<range> m"` —
    /// the range unit is the wire's, not the shooter's, and an app showing those labels
    /// converts them itself.
    Bdc { drops: Vec<(f64, f64)> },
}

/// Optional focal-plane overrides for a GENERATED reticle.
///
/// Every generator returns FFP with `reference_magnification` 1.0, and the generator docs
/// say callers wanting SFP "set those two fields afterwards". On the bridge there is no
/// afterwards — the caller never holds the struct — so the overrides have to travel with
/// the request or the generator path could only ever produce FFP reticles, which is half
/// the optics on the market.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticlePlaneOverride {
    #[serde(default)]
    focal_plane: Option<crate::reticle::FocalPlane>,
    #[serde(default)]
    reference_magnification: Option<f64>,
}

/// Build the description a request is about, from whichever source it supplied.
///
/// `Err` is a finished error envelope, matching how the other commands' helpers report.
fn resolve_reticle(
    command: &'static str,
    reticle: Option<crate::reticle::ReticleDescription>,
    generator: Option<ReticleGenerator>,
    catalog: Option<String>,
    plane: &ReticlePlaneOverride,
) -> Result<crate::reticle::ReticleDescription, String> {
    use crate::reticle::ReticleDescription;

    // Exactly one source. Counted rather than matched pair-by-pair, so adding a fourth
    // source later cannot quietly reintroduce a precedence rule.
    let named = [reticle.is_some(), generator.is_some(), catalog.is_some()]
        .iter()
        .filter(|named| **named)
        .count();
    if named > 1 {
        return Err(error(
            BridgeErrorCode::InvalidRequest,
            format!(
                "{command} takes exactly one of 'reticle', 'generator' or 'catalog' — \
                 supplying more than one is a caller bug rather than a preference"
            ),
            None,
        ));
    }

    // Resolved before the match below so a catalog id reaches it as a full description.
    let reticle = match (reticle, catalog) {
        (reticle, None) => reticle,
        (_, Some(id)) => match crate::reticle_catalog::by_id(&id) {
            Some(Ok(built)) => Some(built),
            Some(Err(err)) => return Err(reticle_error_envelope(command, &err)),
            // Not a default: a caller asking for a specific reticle and silently getting a
            // different one would be shown holds for glass they are not looking through.
            None => {
                return Err(error(
                    BridgeErrorCode::InvalidRequest,
                    format!(
                        "{command}: no catalog reticle with id {id:?}; \
                         'reticle.catalog' lists what this build has"
                    ),
                    None,
                ))
            }
        },
    };

    let mut described = match (reticle, generator) {
        (Some(_), Some(_)) => unreachable!("the count above rejected more than one source"),
        (None, None) => {
            return Err(error(
                BridgeErrorCode::InvalidRequest,
                format!(
                    "{command} requires one of 'reticle' (a full description), \
                     'generator' or 'catalog'"
                ),
                None,
            ))
        }
        (Some(reticle), None) => {
            if plane.focal_plane.is_some() || plane.reference_magnification.is_some() {
                return Err(error(
                    BridgeErrorCode::InvalidRequest,
                    format!(
                        "{command}: 'focal_plane' and 'reference_magnification' override a \
                         GENERATED reticle; a 'reticle' or 'catalog' one already carries its own"
                    ),
                    None,
                ));
            }
            reticle
        }
        (None, Some(generator)) => {
            let built = match generator {
                ReticleGenerator::MilGrid {
                    spacing_mil,
                    extent_mil,
                } => ReticleDescription::mil_grid(spacing_mil, extent_mil),
                ReticleGenerator::Tree {
                    rows,
                    row_spacing_mil,
                    spread_step_mil,
                } => ReticleDescription::tree(rows, row_spacing_mil, spread_step_mil),
                ReticleGenerator::Bdc { drops } => ReticleDescription::bdc_from_drops(&drops),
            };
            match built {
                Ok(built) => built,
                Err(err) => return Err(reticle_error_envelope(command, &err)),
            }
        }
    };

    if let Some(focal_plane) = plane.focal_plane {
        described.focal_plane = focal_plane;
    }
    if let Some(reference_magnification) = plane.reference_magnification {
        described.reference_magnification = reference_magnification;
    }
    Ok(described)
}

/// A [`crate::reticle::ReticleError`] as this family's structured error envelope.
///
/// Follows the `true.*` convention established by `true.dsf`: `error.code` stays
/// `command_failed`, and a stable machine-readable `reason` rides in `error.details`
/// alongside the human message, so a front end can render its own wording without
/// parsing prose. The reasons are the error's own variants, which is what makes them
/// stable — a new variant is a new reason, never a re-spelling of an old one.
fn reticle_error_envelope(command: &'static str, err: &crate::reticle::ReticleError) -> String {
    use crate::reticle::ReticleError;

    let mut details = json!({ "reason": match err {
        ReticleError::NonPositiveMagnification { .. } => "non_positive_magnification",
        ReticleError::NonPositiveReferenceMagnification { .. } => {
            "non_positive_reference_magnification"
        }
        ReticleError::NoMarks => "no_marks",
        ReticleError::TooManyMarks { .. } => "too_many_marks",
        ReticleError::NonFiniteMark { .. } => "non_finite_mark",
        ReticleError::NonFiniteHold { .. } => "non_finite_hold",
        ReticleError::InvalidGeneratorParameter { .. } => "invalid_generator_parameter",
        ReticleError::InvalidSpec(_) => "invalid_spec",
    }});

    // The numbers the variant carries, so a caller can point at the offending input
    // without re-parsing the sentence.
    match err {
        ReticleError::NonPositiveMagnification { magnification } => {
            details["magnification"] = json!(magnification);
        }
        ReticleError::NonPositiveReferenceMagnification {
            reference_magnification,
        } => {
            details["reference_magnification"] = json!(reference_magnification);
        }
        ReticleError::TooManyMarks { count, max } => {
            details["count"] = json!(count);
            details["max"] = json!(max);
        }
        ReticleError::NonFiniteMark { index } => {
            details["index"] = json!(index);
        }
        ReticleError::InvalidGeneratorParameter {
            parameter,
            value,
            rule,
        } => {
            details["parameter"] = json!(parameter);
            details["value"] = json!(value);
            details["rule"] = json!(rule);
        }
        ReticleError::NoMarks
        | ReticleError::NonFiniteHold { .. }
        | ReticleError::InvalidSpec(_) => {}
    }
    details["command"] = json!(command);

    error(
        BridgeErrorCode::CommandFailed,
        err.to_string(),
        Some(details),
    )
}

/// The nearest mark a hold landed by, reported so a caller does not index back into the
/// description and re-apply the scale itself.
///
/// BOTH positions are given because they answer different questions and are not the same
/// number on an SFP optic: `nominal` is the mark as authored (what a reticle diagram
/// prints), `true_angular` is where it actually sits at this magnification (what the hold
/// is measured against).
fn nearest_mark_value(
    reticle: &crate::reticle::ReticleDescription,
    hold: &crate::reticle::ReticleHold,
) -> Value {
    match hold
        .nearest_mark
        .and_then(|index| reticle.marks.get(index).map(|mark| (index, mark)))
    {
        Some((index, mark)) => json!({
            "index": index,
            "kind": mark.kind.as_str(),
            "label": mark.label,
            "nominal": { "down_mil": mark.down_mil, "right_mil": mark.right_mil },
            "true_angular": {
                "down_mil": mark.down_mil * hold.mark_scale,
                "right_mil": mark.right_mil * hold.mark_scale,
            },
        }),
        None => Value::Null,
    }
}

/// `reticle.catalog`: the named reticles this build can produce.
///
/// Result: `{ "reticles": [{ "id", "display_name", "source", "notes", "mark_count" }, ..] }`.
/// An app populates a picker from this rather than hardcoding ids, so a reticle added to a
/// later engine appears without an app release — and one REMOVED does not leave a picker
/// offering something the engine will refuse.
///
/// `source` is carried onto the wire deliberately. A wrong subtension is a hold that is
/// wrong by a mark and looks entirely normal, so the provenance that makes it traceable
/// belongs where the person looking at the reticle can see it, not only in the crate.
fn run_reticle_catalog() -> String {
    let reticles: Vec<Value> = crate::reticle_catalog::catalog()
        .into_iter()
        .map(|entry| {
            // Listed entries are built by `by_id` under test, so this cannot be the place a
            // caller learns an entry is broken; report what is known rather than failing the
            // whole listing over one bad row.
            let mark_count = crate::reticle_catalog::by_id(entry.id)
                .and_then(|built| built.ok())
                .map(|reticle| reticle.marks.len());
            json!({
                "id": entry.id,
                "display_name": entry.display_name,
                "source": entry.source,
                "notes": entry.notes,
                "mark_count": mark_count,
            })
        })
        .collect();
    success("reticle.catalog", json!({ "reticles": reticles }))
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticleDescribeRequest {
    #[serde(default)]
    reticle: Option<crate::reticle::ReticleDescription>,
    #[serde(default)]
    generator: Option<ReticleGenerator>,
    /// A `reticle.catalog` id. The third and simplest way to name a reticle,
    /// and the one an app picker uses.
    #[serde(default)]
    catalog: Option<String>,
    #[serde(default)]
    focal_plane: Option<crate::reticle::FocalPlane>,
    #[serde(default)]
    reference_magnification: Option<f64>,
    /// When present, the result also carries every mark's TRUE angular position at this
    /// magnification. Absent means "just the description".
    #[serde(default)]
    magnification: Option<f64>,
}

/// `reticle.describe`: resolve a reticle — supplied or generated — and report it.
///
/// Result: `{ "reticle": <ReticleDescription>, "focal_plane": "FFP"|"SFP",
/// "magnification_dependent": bool, "mark_count": n, "scaled": {...}|null }`.
///
/// This is the command a picker and a drawing are built on. `scaled` is present only when
/// the request supplied a magnification, and carries `mark_scale` plus the marks in true
/// angular space — for an FFP reticle that is the nominal marks and a scale of exactly
/// 1.0, which is worth returning anyway so a caller has one code path.
fn run_reticle_describe(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'reticle.describe' requires a request payload ({\"reticle\": ...} or {\"generator\": ...})",
            None,
        );
    }
    let request: ReticleDescribeRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("reticle.describe request rejected: {err}"),
                None,
            )
        }
    };

    let plane = ReticlePlaneOverride {
        focal_plane: request.focal_plane,
        reference_magnification: request.reference_magnification,
    };
    let reticle = match resolve_reticle(
        "reticle.describe",
        request.reticle,
        request.generator,
        request.catalog,
        &plane,
    ) {
        Ok(reticle) => reticle,
        Err(envelope) => return envelope,
    };

    // Validate even when no magnification was asked for: a description that cannot be
    // held on is not a description worth handing back as if it were usable.
    if let Err(err) = reticle.validate() {
        return reticle_error_envelope("reticle.describe", &err);
    }

    let scaled = match request.magnification {
        Some(magnification) => match reticle.scaled_marks(magnification) {
            Ok(marks) => json!({
                "magnification": magnification,
                "mark_scale": reticle.mark_scale(magnification),
                "marks": marks
                    .iter()
                    .map(|mark| json!({
                        "down_mil": mark.down_mil,
                        "right_mil": mark.right_mil,
                    }))
                    .collect::<Vec<_>>(),
            }),
            Err(err) => return reticle_error_envelope("reticle.describe", &err),
        },
        None => Value::Null,
    };

    match serde_json::to_value(&reticle) {
        Ok(described) => success(
            "reticle.describe",
            json!({
                "reticle": described,
                "focal_plane": reticle.focal_plane.label(),
                "magnification_dependent": reticle.focal_plane.is_magnification_dependent(),
                "mark_count": reticle.marks.len(),
                "scaled": scaled,
            }),
        ),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize reticle description: {err}"),
            None,
        ),
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticleHoldRequest {
    #[serde(default)]
    reticle: Option<crate::reticle::ReticleDescription>,
    #[serde(default)]
    generator: Option<ReticleGenerator>,
    /// A `reticle.catalog` id. The third and simplest way to name a reticle,
    /// and the one an app picker uses.
    #[serde(default)]
    catalog: Option<String>,
    #[serde(default)]
    focal_plane: Option<crate::reticle::FocalPlane>,
    #[serde(default)]
    reference_magnification: Option<f64>,
    /// Angular drop below the line of sight, positive = below: the come-up the shooter
    /// would otherwise dial.
    drop_mil: f64,
    /// Angular wind deflection, positive = the bullet goes RIGHT.
    #[serde(default)]
    wind_mil: f64,
    /// The optic's CURRENT magnification. Required on every focal plane — an FFP hold does
    /// not depend on it, but zero magnification is not a physical optic and
    /// `hold_point_in_reticle` rejects it on both planes rather than masking a caller bug.
    magnification: f64,
}

/// `reticle.hold`: which mark to hold on.
///
/// The command this family exists for. Given a firing solution already reduced to angles
/// (`drop_mil` / `wind_mil`) and the optic's current magnification, place it in the
/// reticle and report the nearest mark.
///
/// ⚠️ THE ANGLES ARE INPUTS, NOT A SOLVE. This command runs no physics — it is the mark
/// search around [`crate::reticle::hold_point_in_reticle`], which is why the module can
/// keep its no-physics property. A caller gets `drop_mil` and `wind_mil` from `solve` (or
/// a card row) and passes them in; nothing here re-derives them, and a stale angle in
/// produces a confident hold out.
///
/// Result: `{ "hold": <ReticleHold>, "nearest_mark": {...}|null }`. `hold.off_reticle`
/// is the field a UI must not ignore: true means the solution has run off the marked part
/// of the reticle and there is nothing honest to hold on.
fn run_reticle_hold(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'reticle.hold' requires a request payload ({reticle|generator, drop_mil, magnification})",
            None,
        );
    }
    let request: ReticleHoldRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("reticle.hold request rejected: {err}"),
                None,
            )
        }
    };

    let plane = ReticlePlaneOverride {
        focal_plane: request.focal_plane,
        reference_magnification: request.reference_magnification,
    };
    let reticle = match resolve_reticle(
        "reticle.hold",
        request.reticle,
        request.generator,
        request.catalog,
        &plane,
    ) {
        Ok(reticle) => reticle,
        Err(envelope) => return envelope,
    };

    let hold = match crate::reticle::hold_point_in_reticle(
        request.drop_mil,
        request.wind_mil,
        request.magnification,
        &reticle,
    ) {
        Ok(hold) => hold,
        Err(err) => return reticle_error_envelope("reticle.hold", &err),
    };

    let nearest = nearest_mark_value(&reticle, &hold);
    match serde_json::to_value(&hold) {
        Ok(hold_value) => success(
            "reticle.hold",
            json!({
                "hold": hold_value,
                "nearest_mark": nearest,
            }),
        ),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize reticle hold: {err}"),
            None,
        ),
    }
}

/// Most holds `reticle.holds` will place in one call.
///
/// A sampled trajectory is the thing this exists for and runs to a hundred-odd rows; the
/// cap is far above that and exists so a malformed request cannot ask for unbounded work.
pub const MAX_RETICLE_HOLDS: usize = 4096;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticleHoldsEntry {
    drop_mil: f64,
    #[serde(default)]
    wind_mil: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticleHoldsRequest {
    #[serde(default)]
    reticle: Option<crate::reticle::ReticleDescription>,
    #[serde(default)]
    generator: Option<ReticleGenerator>,
    #[serde(default)]
    catalog: Option<String>,
    #[serde(default)]
    focal_plane: Option<crate::reticle::FocalPlane>,
    #[serde(default)]
    reference_magnification: Option<f64>,
    magnification: f64,
    /// The firing solutions to place, in the caller's own order, which the result
    /// preserves.
    holds: Vec<ReticleHoldsEntry>,
}

/// `reticle.holds`: [`run_reticle_hold`] for a whole table, in one call.
///
/// A holdover COLUMN needs a hold per row, and a sampled trajectory is a hundred-odd rows.
/// Done one at a time that is a hundred round trips through the FFI with a JSON encode and
/// decode each; done in the app it is the nearest-mark search re-implemented on two
/// platforms, which is the second spelling `AdjustmentConversion` exists to argue against.
/// So it is one call: the reticle is resolved and validated ONCE, and the search runs per
/// row inside the engine.
///
/// Result: `{ "mark_scale": s, "holds": [ <the same shape reticle.hold returns>, .. ] }`,
/// in request order. `mark_scale` is hoisted because it is a property of the reticle and
/// the magnification, identical for every row, and repeating it per row would invite a
/// caller to wonder when it might differ.
///
/// The whole request fails on the first unusable hold rather than returning a list with a
/// hole in it: a column that silently skips a row is worse than one that does not draw.
fn run_reticle_holds(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'reticle.holds' requires a request payload ({reticle|generator|catalog, magnification, holds})",
            None,
        );
    }
    let request: ReticleHoldsRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("reticle.holds request rejected: {err}"),
                None,
            )
        }
    };

    if request.holds.len() > MAX_RETICLE_HOLDS {
        return error(
            BridgeErrorCode::ResourceLimit,
            format!(
                "reticle.holds was given {} holds; the limit is {MAX_RETICLE_HOLDS}",
                request.holds.len()
            ),
            None,
        );
    }

    let plane = ReticlePlaneOverride {
        focal_plane: request.focal_plane,
        reference_magnification: request.reference_magnification,
    };
    let reticle = match resolve_reticle(
        "reticle.holds",
        request.reticle,
        request.generator,
        request.catalog,
        &plane,
    ) {
        Ok(reticle) => reticle,
        Err(envelope) => return envelope,
    };

    let mut placed = Vec::with_capacity(request.holds.len());
    let mut scale = 1.0;
    for entry in &request.holds {
        let hold = match crate::reticle::hold_point_in_reticle(
            entry.drop_mil,
            entry.wind_mil,
            request.magnification,
            &reticle,
        ) {
            Ok(hold) => hold,
            // Defensive, and unreachable from the wire: a bad reticle and a non-positive
            // magnification are both settled before this loop, and the only error left is
            // a non-finite hold, which JSON cannot carry (serde_json refuses an
            // out-of-range literal, and infinity encodes as null). Kept because this is
            // also reachable from the Rust API, and because a silently skipped row would
            // be worse than a refused request.
            Err(err) => return reticle_error_envelope("reticle.holds", &err),
        };
        scale = hold.mark_scale;
        let nearest = nearest_mark_value(&reticle, &hold);
        match serde_json::to_value(&hold) {
            Ok(mut value) => {
                // Hoisted to the top level; see the doc comment.
                if let Some(object) = value.as_object_mut() {
                    object.remove("mark_scale");
                }
                placed.push(json!({ "hold": value, "nearest_mark": nearest }));
            }
            Err(err) => {
                return error(
                    BridgeErrorCode::InternalError,
                    format!("failed to serialize reticle hold: {err}"),
                    None,
                )
            }
        }
    }

    success(
        "reticle.holds",
        json!({ "mark_scale": scale, "holds": placed }),
    )
}

/// Largest reticle document `reticle.import` will parse.
///
/// The envelope cap ([`MAX_REQUEST_BYTES`]) already bounds the whole request, but a
/// document-specific cap gives a caller a message naming the document rather than the
/// envelope, and bounds the parse before it starts. Both formats are TEXT and travel
/// inline — unlike `.a7p`, which is binary and needs base64 — so there is no inflation
/// factor to leave headroom for.
pub const MAX_RETICLE_DOCUMENT_BYTES: usize = 512 * 1024;

#[derive(Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
enum ReticleDocumentFormat {
    /// The Ventum JSON reticle spec (MBA-1440).
    Ventum,
    /// The third-party `.reticle` XML drawing format (MBA-1544).
    ReticleXml,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReticleImportRequest {
    format: ReticleDocumentFormat,
    /// The document itself, as text.
    document: String,
}

/// `reticle.import`: turn a third-party reticle document into a
/// [`crate::reticle::ReticleDescription`].
///
/// Result: `{ "reticle": <ReticleDescription>, "format": "ventum"|"reticle_xml",
/// "report": {...} }`. The report is format-specific and is NOT decoration: both importers
/// drop elements they cannot turn into hold points, and a caller that ignores the report
/// shows a shooter a reticle missing marks their drawing had. The same obligation
/// `profile.import_a7p` discharges with its `unmapped` list.
fn run_reticle_import(inner: &Value) -> String {
    if inner.is_null() {
        return error(
            BridgeErrorCode::InvalidRequest,
            "'reticle.import' requires a request payload ({\"format\": ..., \"document\": ...})",
            None,
        );
    }
    let request: ReticleImportRequest = match serde_json::from_value(inner.clone()) {
        Ok(request) => request,
        Err(err) => {
            return error(
                BridgeErrorCode::InvalidRequest,
                format!("reticle.import request rejected: {err}"),
                None,
            )
        }
    };

    if request.document.len() > MAX_RETICLE_DOCUMENT_BYTES {
        return error(
            BridgeErrorCode::ResourceLimit,
            format!(
                "reticle document is {} bytes; the limit is {MAX_RETICLE_DOCUMENT_BYTES}",
                request.document.len()
            ),
            None,
        );
    }

    let (reticle, report, format) = match request.format {
        ReticleDocumentFormat::Ventum => {
            match crate::reticle_import::import_ventum_reticle_with_report(&request.document) {
                Ok((reticle, report)) => (
                    reticle,
                    json!({
                        "dropped_elements": report.dropped_elements,
                        "dropped_element_types": report
                            .dropped_element_types
                            .iter()
                            .map(|(tag, count)| json!({ "tag": tag, "count": count }))
                            .collect::<Vec<_>>(),
                        // Arcs are REPORTED, never imported as marks: a horseshoe is a
                        // shape a shooter indexes on, and which of its points counts as an
                        // aiming point is the caller's decision, not this bridge's. Apex
                        // and the two tips are the points the importer resolves.
                        "arcs": report
                            .arcs
                            .iter()
                            .map(|arc| json!({
                                "center": { "down_mil": arc.center.down_mil, "right_mil": arc.center.right_mil },
                                "radius_mil": arc.radius_mil,
                                "start_degrees": arc.start_degrees,
                                "end_degrees": arc.end_degrees,
                                "sweep_degrees": arc.sweep_degrees,
                                "start_tip": { "down_mil": arc.start_tip.down_mil, "right_mil": arc.start_tip.right_mil },
                                "apex": { "down_mil": arc.apex.down_mil, "right_mil": arc.apex.right_mil },
                                "end_tip": { "down_mil": arc.end_tip.down_mil, "right_mil": arc.end_tip.right_mil },
                            }))
                            .collect::<Vec<_>>(),
                        // Both shortfalls the report declares about itself, carried rather
                        // than dropped: `arcs.len() + arcs_unresolved` is the arc tally, and
                        // `circle_repeats_unreadable` says the element count is low by an
                        // unknown amount.
                        "arcs_unresolved": report.arcs_unresolved,
                        "circle_repeats_unreadable": report.circle_repeats_unreadable,
                    }),
                    "ventum",
                ),
                Err(err) => return reticle_error_envelope("reticle.import", &err),
            }
        }
        ReticleDocumentFormat::ReticleXml => {
            match crate::reticle_document_import::import_reticle_document_with_report(
                &request.document,
            ) {
                Ok((reticle, report)) => (
                    reticle,
                    json!({
                        "drawing_elements_dropped": report.drawing_elements_dropped,
                        "dropped_element_tags": report
                            .dropped_element_tags
                            .iter()
                            .map(|(tag, count)| json!({ "tag": tag, "count": count }))
                            .collect::<Vec<_>>(),
                        "canvas_mil": report
                            .canvas_mil
                            .map(|(x, y)| json!({ "size_x_mil": x, "size_y_mil": y })),
                        "zero_in_canvas_mil": report
                            .zero_in_canvas_mil
                            .map(|(x, y)| json!({ "zero_x_mil": x, "zero_y_mil": y })),
                    }),
                    "reticle_xml",
                ),
                Err(err) => return reticle_error_envelope("reticle.import", &err),
            }
        }
    };

    match serde_json::to_value(&reticle) {
        Ok(imported) => success(
            "reticle.import",
            json!({
                "reticle": imported,
                "format": format,
                "mark_count": reticle.marks.len(),
                "report": report,
            }),
        ),
        Err(err) => error(
            BridgeErrorCode::InternalError,
            format!("failed to serialize imported reticle: {err}"),
            None,
        ),
    }
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
            commands.contains(&"profile.export_a7p".to_string()),
            cfg!(feature = "profile-export"),
            "profile.export_a7p must be listed exactly when compiled in"
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
        assert!(
            decode_base64("Zm9v!").is_err(),
            "non-alphabet byte is rejected"
        );
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
            out["error"]["message"]
                .as_str()
                .unwrap()
                .contains("card v1 document"),
            "{out}"
        );
    }

    /// The output cap's boundary, both sides. Generating a genuinely over-cap dope card
    /// would take tens of thousands of rows, so the predicate is tested directly — see
    /// `pdf_over_cap_error`'s own comment.
    #[cfg(feature = "pdf")]
    #[test]
    fn pdf_output_cap_refuses_only_over_the_limit() {
        assert!(pdf_over_cap_error(0, 0, 0).is_none());
        assert!(
            pdf_over_cap_error(MAX_PDF_BYTES, 6, 1).is_none(),
            "a document exactly at the cap fits"
        );
        let envelope: Value =
            serde_json::from_str(&pdf_over_cap_error(MAX_PDF_BYTES + 1, 6, 1).expect("over cap"))
                .unwrap();
        assert_eq!(envelope["ok"], false);
        assert_eq!(envelope["error"]["code"], "resource_limit");
        let message = envelope["error"]["message"].as_str().unwrap();
        assert!(message.contains("dope card"), "{envelope}");
        // What is true of the document, not advice about controls a saved card lacks.
        assert!(message.contains("6 rows"), "{envelope}");
        assert!(message.contains("1 pages"), "{envelope}");
        for absent in ["coarsen", "shorten"] {
            assert!(!message.contains(absent), "{envelope}");
        }
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
            let bytes: Vec<u8> = (0..len)
                .map(|i| (i as u8).wrapping_mul(37).wrapping_add(11))
                .collect();
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

    // -----------------------------------------------------------------------
    // reticle.* — MBA-1558
    // -----------------------------------------------------------------------

    /// A 0.5 mil grid out to 5 mil, as a full description rather than a generator, so a
    /// test can assert against marks it wrote itself.
    fn hash_reticle() -> Value {
        json!({
            "name": "test cross",
            "focal_plane": "ffp",
            "reference_magnification": 1.0,
            "marks": [
                {"down_mil": 0.0, "right_mil": 0.0, "kind": "center"},
                {"down_mil": 2.0, "right_mil": 0.0, "kind": "hash", "label": "2 mil"},
                {"down_mil": 4.0, "right_mil": 0.0, "kind": "hash"},
                {"down_mil": 0.0, "right_mil": 2.0, "kind": "hash"},
            ],
        })
    }

    #[test]
    fn capabilities_lists_the_reticle_family() {
        // The whole point of MBA-1558: an app feature-detects through this list, so a
        // command that dispatches but is not listed is a command no app will call.
        let out = call(json!({"api_version": 1, "command": "meta.capabilities"}));
        let commands: Vec<String> =
            serde_json::from_value(out["result"]["commands"].clone()).unwrap();
        for command in ["reticle.describe", "reticle.hold", "reticle.import"] {
            assert!(
                commands.contains(&command.to_string()),
                "{command} dispatches but meta.capabilities does not list it"
            );
        }
    }

    #[test]
    fn hold_reports_the_nearest_mark_and_its_label() {
        let out = call(json!({
            "api_version": 1,
            "command": "reticle.hold",
            "request": {
                "reticle": hash_reticle(),
                "drop_mil": 2.1,
                "wind_mil": 0.0,
                "magnification": 10.0,
            },
        }));
        assert_eq!(out["ok"], true, "{out}");
        // The hold coordinates are the inputs: this command places a solution, it does
        // not compute one.
        assert_eq!(out["result"]["hold"]["down_mil"], 2.1);
        assert_eq!(out["result"]["hold"]["right_mil"], 0.0);
        assert_eq!(out["result"]["hold"]["off_reticle"], false);
        assert_eq!(out["result"]["nearest_mark"]["index"], 1);
        assert_eq!(out["result"]["nearest_mark"]["label"], "2 mil");
        assert_eq!(out["result"]["nearest_mark"]["kind"], "hash");
    }

    #[test]
    fn an_sfp_reticle_moves_its_marks_and_the_hold_follows() {
        // The reason `magnification` is required on every plane. The same solution finds a
        // DIFFERENT mark on an SFP optic depending on the zoom ring, and an app that
        // ignored this would draw a confident hold on the wrong hash.
        let mut sfp = hash_reticle();
        sfp["focal_plane"] = json!("sfp");
        sfp["reference_magnification"] = json!(10.0);

        let at_reference = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"reticle": sfp, "drop_mil": 2.0, "magnification": 10.0},
        }));
        assert_eq!(at_reference["result"]["hold"]["mark_scale"], 1.0);
        assert_eq!(at_reference["result"]["nearest_mark"]["index"], 1);
        assert_eq!(
            at_reference["result"]["nearest_mark"]["true_angular"]["down_mil"],
            2.0
        );

        // Halve the magnification and every subtension doubles, so the 2 mil hash now
        // sits at 4 mil true and the 4 mil hash at 8.
        let zoomed_out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"reticle": sfp, "drop_mil": 2.0, "magnification": 5.0},
        }));
        assert_eq!(zoomed_out["result"]["hold"]["mark_scale"], 2.0);
        assert_eq!(
            zoomed_out["result"]["nearest_mark"]["true_angular"]["down_mil"], 0.0,
            "at 5x the nearest mark to a 2 mil hold is now CENTER, not the 2 mil hash"
        );
        assert_eq!(zoomed_out["result"]["nearest_mark"]["index"], 0);
        // And the nominal position is unchanged, which is why both are reported.
        assert_eq!(
            zoomed_out["result"]["nearest_mark"]["nominal"]["down_mil"],
            0.0
        );

        // A mark whose coordinates are NOT zero, because the assertions above cannot tell
        // a scaled position from an unscaled one: 0.0 * 2.0 is 0.0. Proved by mutation --
        // dropping the scale from `true_angular` left every case above green.
        //
        // At 5x on a 10x-reference SFP reticle the 2 mil hash sits at 4 mil true, so a
        // 4 mil hold lands exactly on it. The two positions differ by the scale, which is
        // the whole reason both are reported.
        let on_the_scaled_hash = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"reticle": sfp, "drop_mil": 4.0, "magnification": 5.0},
        }));
        assert_eq!(
            on_the_scaled_hash["result"]["nearest_mark"]["nominal"]["down_mil"], 2.0,
            "as authored, the mark is at 2 mil"
        );
        assert_eq!(
            on_the_scaled_hash["result"]["nearest_mark"]["true_angular"]["down_mil"], 4.0,
            "at 5x it subtends 4 mil, which is what the hold is measured against"
        );
        assert_eq!(
            on_the_scaled_hash["result"]["hold"]["nearest_mark_distance_mil"], 0.0,
            "the hold lands exactly on it"
        );
    }

    #[test]
    fn a_hold_off_the_marked_part_says_so() {
        // The field a UI must not ignore. 40 mil is far outside a reticle whose lowest
        // mark is at 4.
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"reticle": hash_reticle(), "drop_mil": 40.0, "magnification": 10.0},
        }));
        assert_eq!(out["ok"], true);
        assert_eq!(out["result"]["hold"]["off_reticle"], true);
    }

    #[test]
    fn describe_builds_from_a_generator_and_reports_the_plane() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {"generator": {"kind": "mil_grid", "spacing_mil": 1.0, "extent_mil": 5.0}},
        }));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["focal_plane"], "FFP");
        assert_eq!(out["result"]["magnification_dependent"], false);
        assert!(out["result"]["mark_count"].as_u64().unwrap() > 1);
        // No magnification asked for, so no scaled block — not an empty one.
        assert!(out["result"]["scaled"].is_null());
    }

    #[test]
    fn a_generated_reticle_can_be_made_sfp() {
        // Every generator returns FFP and the library's answer is "set those two fields
        // afterwards". On the bridge there is no afterwards, so the overrides travel with
        // the request — without them the generator path could only ever make FFP.
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {
                "generator": {"kind": "mil_grid", "spacing_mil": 1.0, "extent_mil": 5.0},
                "focal_plane": "sfp",
                "reference_magnification": 12.0,
                "magnification": 6.0,
            },
        }));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["focal_plane"], "SFP");
        assert_eq!(out["result"]["magnification_dependent"], true);
        assert_eq!(out["result"]["scaled"]["mark_scale"], 2.0);
    }

    #[test]
    fn a_bdc_ladder_carries_its_range_labels() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {"generator": {"kind": "bdc", "drops": [[300.0, 1.2], [400.0, 2.4]]}},
        }));
        assert_eq!(out["ok"], true, "{out}");
        let marks = out["result"]["reticle"]["marks"]
            .as_array()
            .unwrap()
            .clone();
        let labels: Vec<&str> = marks.iter().filter_map(|m| m["label"].as_str()).collect();
        // The generator labels in METRES, which is the wire's unit and not the shooter's.
        assert_eq!(labels, vec!["300 m", "400 m"]);
    }

    #[test]
    fn naming_a_reticle_twice_is_refused_rather_than_resolved() {
        // Supplying both is a caller bug. Picking a winner would hide it.
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {
                "reticle": hash_reticle(),
                "generator": {"kind": "mil_grid", "spacing_mil": 1.0, "extent_mil": 5.0},
                "drop_mil": 1.0,
                "magnification": 10.0,
            },
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_request");
        assert!(out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("exactly one"));
    }

    #[test]
    fn a_catalog_id_is_the_third_way_to_name_a_reticle() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"catalog": "mil-dot", "drop_mil": 3.0, "magnification": 10.0},
        }));
        assert_eq!(out["ok"], true, "{out}");
        // A mil-dot has a dot exactly 3 mil down, so the hold lands on it.
        assert_eq!(out["result"]["hold"]["nearest_mark_distance_mil"], 0.0);
        assert_eq!(out["result"]["nearest_mark"]["nominal"]["down_mil"], 3.0);
        assert_eq!(out["result"]["nearest_mark"]["label"], "3 mil down");
    }

    #[test]
    fn an_unknown_catalog_id_is_refused_rather_than_substituted() {
        // Handing back a different reticle would show a shooter holds for glass they are
        // not looking through.
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"catalog": "tremor-9000", "drop_mil": 3.0, "magnification": 10.0},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_request");
        assert!(out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("reticle.catalog"));
    }

    #[test]
    fn naming_a_reticle_three_ways_is_refused_too() {
        // The count-based check exists so a third source could not reintroduce precedence.
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {
                "catalog": "mil-dot",
                "generator": {"kind": "mil_grid", "spacing_mil": 1.0, "extent_mil": 5.0},
            },
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_request");
    }

    #[test]
    fn the_catalog_lists_what_the_other_commands_will_accept() {
        // The contract an app picker is built on: every listed id must resolve, or the
        // picker offers something the engine refuses.
        let out = call(json!({"api_version": 1, "command": "reticle.catalog"}));
        assert_eq!(out["ok"], true, "{out}");
        let reticles = out["result"]["reticles"].as_array().unwrap().clone();
        assert!(!reticles.is_empty());
        for entry in reticles {
            let id = entry["id"].as_str().unwrap();
            // Provenance travels to the wire on purpose — see reticle_catalog's header.
            assert!(
                !entry["source"].as_str().unwrap().is_empty(),
                "{id} has no source"
            );
            assert!(
                !entry["notes"].as_str().unwrap().is_empty(),
                "{id} has no notes"
            );
            assert!(
                entry["mark_count"].as_u64().unwrap() > 0,
                "{id} has no marks"
            );

            let held = call(json!({
                "api_version": 1, "command": "reticle.hold",
                "request": {"catalog": id, "drop_mil": 1.0, "magnification": 10.0},
            }));
            assert_eq!(held["ok"], true, "listed id {id} does not resolve: {held}");
        }
    }

    #[test]
    fn a_batch_of_holds_comes_back_in_request_order() {
        // The column this exists for: one call per table, not one per row.
        let out = call(json!({
            "api_version": 1, "command": "reticle.holds",
            "request": {
                "catalog": "mil-dot",
                "magnification": 10.0,
                "holds": [
                    {"drop_mil": 1.0},
                    {"drop_mil": 3.0, "wind_mil": 0.0},
                    {"drop_mil": 2.0},
                ],
            },
        }));
        assert_eq!(out["ok"], true, "{out}");
        let holds = out["holds"].clone();
        assert!(holds.is_null(), "holds live under result, not the envelope");
        let holds = out["result"]["holds"].as_array().unwrap().clone();
        assert_eq!(holds.len(), 3);
        // Order is the caller's, NOT sorted — a table row must line up with its hold.
        assert_eq!(holds[0]["nearest_mark"]["nominal"]["down_mil"], 1.0);
        assert_eq!(holds[1]["nearest_mark"]["nominal"]["down_mil"], 3.0);
        assert_eq!(holds[2]["nearest_mark"]["nominal"]["down_mil"], 2.0);
        // mark_scale is hoisted: one reticle at one magnification has exactly one.
        assert_eq!(out["result"]["mark_scale"], 1.0);
        assert!(holds[0]["hold"]["mark_scale"].is_null());
    }

    #[test]
    fn a_batch_agrees_with_the_single_hold_command() {
        // The batch must not be a second implementation. If these ever disagree, one of
        // them is doing its own mark search.
        let single = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"catalog": "mil-dot", "drop_mil": 2.4, "wind_mil": 0.7, "magnification": 10.0},
        }));
        let batch = call(json!({
            "api_version": 1, "command": "reticle.holds",
            "request": {
                "catalog": "mil-dot", "magnification": 10.0,
                "holds": [{"drop_mil": 2.4, "wind_mil": 0.7}],
            },
        }));
        assert_eq!(single["ok"], true);
        assert_eq!(batch["ok"], true);
        assert_eq!(
            single["result"]["nearest_mark"],
            batch["result"]["holds"][0]["nearest_mark"]
        );
        assert_eq!(
            single["result"]["hold"]["off_reticle"],
            batch["result"]["holds"][0]["hold"]["off_reticle"]
        );
        assert_eq!(
            single["result"]["hold"]["nearest_mark_distance_mil"],
            batch["result"]["holds"][0]["hold"]["nearest_mark_distance_mil"]
        );
    }

    #[test]
    fn a_malformed_hold_entry_refuses_the_whole_batch() {
        // A column that silently skips a row is worse than one that does not draw, so the
        // request fails rather than returning a list with a hole in it.
        //
        // Note WHICH failure this is. The per-row arm in run_reticle_holds is defensive
        // and unreachable from the wire: `hold_point_in_reticle` rejects a bad reticle and
        // a non-positive magnification, both of which are settled before the loop, and its
        // only remaining error is a non-finite hold — which JSON cannot deliver.
        // serde_json refuses an out-of-range literal outright ("number out of range" for
        // 1e999) and `json!(f64::INFINITY)` encodes as null, so a non-finite drop is
        // stopped by the transport. `NonFiniteHold` is reachable from the Rust API and not
        // from here; an earlier version of this test asserted that reason and failed,
        // which is how the unreachability was found.
        for bad in [
            json!(null),
            json!("3.0"),
            json!({"drop_mil": 1.0, "nonsense": 2}),
        ] {
            let out = call(json!({
                "api_version": 1, "command": "reticle.holds",
                "request": {
                    "catalog": "mil-dot", "magnification": 10.0,
                    "holds": [{"drop_mil": 1.0}, bad],
                },
            }));
            assert_eq!(out["ok"], false, "accepted a malformed entry: {out}");
            assert_eq!(out["error"]["code"], "invalid_request");
        }
    }

    #[test]
    fn a_bad_magnification_fails_the_batch_before_any_row() {
        // The reachable whole-batch failure: settled once, up front, not per row.
        let out = call(json!({
            "api_version": 1, "command": "reticle.holds",
            "request": {
                "catalog": "mil-dot", "magnification": 0.0,
                "holds": [{"drop_mil": 1.0}, {"drop_mil": 2.0}],
            },
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(
            out["error"]["details"]["reason"],
            "non_positive_magnification"
        );
    }

    #[test]
    fn an_empty_batch_is_an_empty_list_not_an_error() {
        // A trajectory with no samples is a legitimate thing to ask about.
        let out = call(json!({
            "api_version": 1, "command": "reticle.holds",
            "request": {"catalog": "mil-dot", "magnification": 10.0, "holds": []},
        }));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["holds"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn too_many_holds_is_a_resource_limit() {
        let holds: Vec<Value> = (0..MAX_RETICLE_HOLDS + 1)
            .map(|i| json!({"drop_mil": (i % 5) as f64}))
            .collect();
        let out = call(json!({
            "api_version": 1, "command": "reticle.holds",
            "request": {"catalog": "mil-dot", "magnification": 10.0, "holds": holds},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "resource_limit");
    }

    #[test]
    fn capabilities_lists_the_catalog_command() {
        let out = call(json!({"api_version": 1, "command": "meta.capabilities"}));
        let commands: Vec<String> =
            serde_json::from_value(out["result"]["commands"].clone()).unwrap();
        assert!(commands.contains(&"reticle.catalog".to_string()));
        assert!(commands.contains(&"reticle.holds".to_string()));
    }

    #[test]
    fn naming_no_reticle_at_all_is_refused() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"drop_mil": 1.0, "magnification": 10.0},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_request");
    }

    #[test]
    fn plane_overrides_beside_a_full_reticle_are_refused() {
        // A supplied description already carries its own plane; silently overriding it
        // would let a request contradict the document it sent.
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {"reticle": hash_reticle(), "focal_plane": "sfp"},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "invalid_request");
    }

    #[test]
    fn a_reticle_error_carries_a_stable_reason_and_its_numbers() {
        // The `true.dsf` convention: code stays command_failed, the machine-readable
        // reason is in details, and the offending value rides along so a caller can point
        // at the input without re-parsing the sentence.
        let out = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {"reticle": hash_reticle(), "drop_mil": 1.0, "magnification": 0.0},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "command_failed");
        assert_eq!(
            out["error"]["details"]["reason"],
            "non_positive_magnification"
        );
        assert_eq!(out["error"]["details"]["magnification"], 0.0);
        assert_eq!(out["error"]["details"]["command"], "reticle.hold");
    }

    #[test]
    fn a_reticle_with_no_marks_is_refused_by_describe_too() {
        // Validation is not deferred to the hold: handing back an unusable description as
        // though it were usable is the failure this guards.
        let mut empty = hash_reticle();
        empty["marks"] = json!([]);
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {"reticle": empty},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["details"]["reason"], "no_marks");
    }

    #[test]
    fn an_unreadable_generator_parameter_names_the_parameter() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.describe",
            "request": {"generator": {"kind": "mil_grid", "spacing_mil": 0.0, "extent_mil": 5.0}},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(
            out["error"]["details"]["reason"],
            "invalid_generator_parameter"
        );
        assert_eq!(out["error"]["details"]["parameter"], "spacing");
    }

    #[test]
    fn importing_a_ventum_document_returns_the_reticle_and_its_report() {
        // Bero's own MBR dot-tree spec, the reticle that produced MBA-1440 — three rows
        // each stamped by a mirrored repeat, expanding to 78 holdable marks. Using the
        // real thing rather than a toy is deliberate: an earlier version of this test sent
        // a document with an invented key, and because it accepted either outcome it
        // passed while reaching nothing.
        let document = r#"{"name":"MBR","plane":"ffp","unit":"mil","spec":[
            {"type":"dot","y":4,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":9,"mirror":true}},
            {"type":"dot","y":8,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":13,"mirror":true}},
            {"type":"dot","y":12,"x":1,"r":0.12,"repeat":{"axis":"x","step":1,"n":17,"mirror":true}}]}"#;

        let out = call(json!({
            "api_version": 1, "command": "reticle.import",
            "request": {"format": "ventum", "document": document},
        }));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["format"], "ventum");
        assert_eq!(out["result"]["mark_count"], 78);
        assert_eq!(out["result"]["reticle"]["name"], "MBR");
        assert_eq!(out["result"]["reticle"]["focal_plane"], "ffp");
        // Nothing was dropped, and the report says so rather than being absent.
        assert_eq!(out["result"]["report"]["dropped_elements"], 0);
        assert_eq!(out["result"]["report"]["arcs_unresolved"], 0);

        // And the imported reticle is immediately usable by the command it exists for.
        let held = call(json!({
            "api_version": 1, "command": "reticle.hold",
            "request": {
                "reticle": out["result"]["reticle"],
                "drop_mil": 8.0,
                "wind_mil": 2.0,
                "magnification": 10.0,
            },
        }));
        assert_eq!(held["ok"], true, "{held}");
        assert_eq!(held["result"]["hold"]["off_reticle"], false);
        assert_eq!(held["result"]["nearest_mark"]["nominal"]["down_mil"], 8.0);
        assert_eq!(held["result"]["nearest_mark"]["nominal"]["right_mil"], 2.0);
        assert_eq!(held["result"]["hold"]["nearest_mark_distance_mil"], 0.0);
    }

    #[test]
    fn a_document_that_is_not_the_named_format_fails_as_a_command_failure() {
        let out = call(json!({
            "api_version": 1, "command": "reticle.import",
            "request": {"format": "reticle_xml", "document": "{\"this\": \"is json\"}"},
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "command_failed");
    }

    #[test]
    fn an_oversize_reticle_document_is_a_resource_limit_not_a_parse_failure() {
        // Named as the document's limit rather than the envelope's, and refused before
        // the parse rather than during it.
        let out = call(json!({
            "api_version": 1, "command": "reticle.import",
            "request": {
                "format": "ventum",
                "document": "x".repeat(MAX_RETICLE_DOCUMENT_BYTES + 1),
            },
        }));
        assert_eq!(out["ok"], false);
        assert_eq!(out["error"]["code"], "resource_limit");
        assert!(out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("reticle document"));
    }

    #[test]
    fn every_reticle_command_rejects_an_unknown_request_field() {
        for command in ["reticle.describe", "reticle.hold", "reticle.import"] {
            let out = call(json!({
                "api_version": 1,
                "command": command,
                "request": {"nonsense": 1},
            }));
            assert_eq!(out["ok"], false, "{command} accepted an unknown field");
            assert_eq!(out["error"]["code"], "invalid_request", "{command}");
        }
    }

    #[test]
    fn every_reticle_command_rejects_a_missing_payload() {
        for command in ["reticle.describe", "reticle.hold", "reticle.import"] {
            let out = call(json!({"api_version": 1, "command": command}));
            assert_eq!(out["ok"], false, "{command} accepted a null payload");
            assert_eq!(out["error"]["code"], "invalid_request", "{command}");
        }
    }

    /// The bridge's own end-to-end proof: a ProfileData in, a real `.a7p` out, and the
    /// file it hands back is one `profile.import_a7p` accepts. Round-trip fidelity of the
    /// individual fields is the encoder's own test; what is proved here is that the two
    /// bridge commands are actually inverse over the wire, base64 and all.
    #[cfg(all(feature = "profile-export", feature = "profile-import"))]
    #[test]
    fn export_a7p_produces_a_file_import_a7p_reads_back() {
        let profile = json!({
            "name": "bridge-export",
            "velocity": 792.0,
            "bc": 0.381,
            "mass": 19.439673,
            "diameter": 8.5852,
            "drag_model": "G7",
            "twist_rate": 254.0,
            "sight_height": 90.0,
            "zero_distance": 100.0,
            "units": "metric",
            "temperature": 15.0,
            "pressure": 1000.0,
            "humidity": 50.0,
            "altitude": 0.0,
            "bullet_name": "300GR OTM",
            "twist_right": false,
            "bullet_length": 45.72,
            "elevation_cf": 0.97,
            "dsf_points": [{"mach": 0.9, "dsf": 1.04}]
        });

        let out = call(json!({
            "api_version": 1, "command": "profile.export_a7p", "request": profile
        }));
        assert_eq!(out["ok"], true, "{out}");
        let result = &out["result"];
        assert!(result["byte_length"].as_u64().unwrap() > 32);

        // The honesty contract: the two fields this profile set that .a7p cannot take
        // are named, and nothing that DID travel is named alongside them.
        let dropped: Vec<String> =
            serde_json::from_value(result["dropped_fields"].clone()).unwrap();
        assert!(dropped.contains(&"elevation_cf".to_string()), "{out}");
        assert!(dropped.contains(&"dsf_points".to_string()), "{out}");
        assert!(!dropped.contains(&"velocity".to_string()), "{out}");
        // ...and the full list names the unsupported fields this profile left empty too,
        // so a caller can tell "had none" from "lost it".
        let not_carried = result["not_carried"].as_array().unwrap();
        let reticle = not_carried
            .iter()
            .find(|e| e["field"] == "reticle")
            .expect("reticle is listed even though this profile has none");
        assert_eq!(reticle["populated"], false, "{out}");
        assert!(reticle["reason"].as_str().unwrap().len() > 10, "{out}");

        let back = call(json!({
            "api_version": 1,
            "command": "profile.import_a7p",
            "request": { "a7p_base64": result["a7p_base64"].clone(), "strict": true }
        }));
        // `strict` refuses on an envelope mismatch, so this passing is also the proof
        // that the MD5 prefix the encoder wrote is correct.
        assert_eq!(back["ok"], true, "{back}");
        let reimported = &back["result"]["profile"];
        assert_eq!(reimported["name"], "bridge-export", "{back}");
        assert_eq!(reimported["drag_model"], "G7", "{back}");
        assert_eq!(reimported["twist_right"], false, "{back}");
    }

    /// A drag model the format cannot name is refused with the structured `reason` the
    /// `true.*` family established, not written out as something else.
    #[cfg(feature = "profile-export")]
    #[test]
    fn export_a7p_refuses_an_unrepresentable_drag_model_with_a_reason() {
        let out = call(json!({
            "api_version": 1,
            "command": "profile.export_a7p",
            // Complete in every other respect, so the refusal can only be about the
            // drag model: the format also REQUIRES a bullet length and a zero
            // distance, and those refusals are checked in the encoder's own tests.
            "request": {
                "name": "g5-load", "velocity": 792.0, "bc": 0.3, "mass": 19.4,
                "diameter": 8.5852, "drag_model": "G5", "units": "metric",
                "bullet_length": 45.72, "zero_distance": 100.0,
                "temperature": 15.0, "pressure": 1000.0, "humidity": 50.0
            }
        }));
        assert_eq!(out["error"]["code"], "command_failed", "{out}");
        assert_eq!(out["error"]["details"]["reason"], "drag_model", "{out}");
        assert_eq!(out["error"]["details"]["drag_model"], "G5", "{out}");
    }

    #[cfg(feature = "profile-export")]
    #[test]
    fn export_a7p_without_a_payload_is_an_invalid_request() {
        let out = call(json!({"api_version": 1, "command": "profile.export_a7p"}));
        assert_eq!(out["error"]["code"], "invalid_request", "{out}");
        assert!(
            out["error"]["message"]
                .as_str()
                .unwrap()
                .contains("ProfileData"),
            "{out}"
        );
    }
}
