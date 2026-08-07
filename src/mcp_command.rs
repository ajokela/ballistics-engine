//! Process transport for the `ballistics mcp` command.
//!
//! This implements a small subset of the Model Context Protocol (MCP) over the stdio
//! transport: JSON-RPC 2.0 messages, one per line, on standard input and standard output.
//! Only what two tools ("solve" and "engine_info") require is implemented — this is
//! deliberately not a general-purpose MCP SDK.
//!
//! Framing follows the stdio transport section of the MCP specification: messages are
//! JSON-RPC requests, notifications, or responses, delimited by newlines, and must not
//! contain embedded newlines. The server must not write anything to stdout that is not a
//! valid MCP message; it may log freely to stderr (this implementation does not).
//!
//! Request semantics for the "solve" tool remain in the transport-free library service,
//! [`ballistics_engine::solve_v1()`], reusing the same decoding entry point as the
//! `solve-json` command (see [`crate::solve_json_command`]).
//!
//! ## Error-reporting convention (documented, not part of the MCP spec)
//!
//! Two different kinds of "the tool call failed" exist, and this server reports them two
//! different ways:
//!
//! - **Malformed tool arguments** — the `tools/call` envelope itself is malformed (missing
//!   `name`/`arguments`), an unknown tool is named, or `solve`'s arguments fail solve-json v1
//!   *structural* validation (unknown/missing fields, wrong JSON types, an unsupported
//!   `schema_version`) — are reported as a JSON-RPC protocol error, code `-32602` (Invalid
//!   params). The MCP tools specification's own error-handling example uses this code for an
//!   unknown tool name, and lists "Invalid arguments" under the same protocol-error category.
//! - **A well-formed request the engine could not solve** — solve-json v1 *semantic* failures
//!   (out-of-range values, conflicting fields caught during resolution, the sample-count
//!   resource limit) and genuine solve failures — are reported as a successful `tools/call`
//!   result with `isError: true`, whose text content is the solve-json v1 error envelope. This
//!   matches the MCP tools specification's "Tool Execution Errors" category (API failures,
//!   invalid input data, business logic errors).
//!
//! Concretely: anything [`ballistics_engine::solve_json::decode_solve_request_v1`] rejects is
//! `-32602`; anything [`ballistics_engine::solve_v1()`] rejects after that is `isError: true`.

use ballistics_engine::solve_json::{decode_solve_request_v1, DRAG_MODEL_WIRE_NAMES_V1};
use serde_json::{json, Value};
use std::io::{self, BufRead, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};

/// Maximum bytes accepted for one JSON-RPC message line, including its trailing newline.
///
/// There is no MCP-specified limit, so this reuses the spirit (and the exact value) of
/// `solve-json`'s `MAX_INPUT_BYTES` (MBA-1299's allocation-bounding convention): a hostile or
/// buggy client sending an unbounded line without a newline must never grow this process's
/// memory without bound. A line at or under the limit is accepted; the first byte over it makes
/// the whole line rejected as oversized, without ever buffering more than one internal read
/// chunk beyond the limit.
const MAX_MCP_LINE_BYTES: usize = 1024 * 1024;

const JSON_RPC_PARSE_ERROR: i64 = -32700;
const JSON_RPC_INVALID_REQUEST: i64 = -32600;
const JSON_RPC_METHOD_NOT_FOUND: i64 = -32601;
const JSON_RPC_INVALID_PARAMS: i64 = -32602;
const JSON_RPC_INTERNAL_ERROR: i64 = -32603;
/// MCP lifecycle: request received before a successful `initialize`. The MCP spec
/// mandates the rejection but reserves no code for it; -32002 follows the LSP
/// ServerNotInitialized convention MCP implementations commonly borrow.
const JSON_RPC_SERVER_NOT_INITIALIZED: i64 = -32002;

/// Fallback `protocolVersion` used only when a client's `initialize` request omits it or sends
/// a non-string value. This server otherwise always echoes the client's requested
/// `protocolVersion` back unchanged (see [`handle_initialize`]) rather than asserting a fixed
/// version of its own: the MCP spec version churns, this server has no version-gated behavior,
/// and disagreeing with a client that sent a well-formed request would only break well-behaved
/// clients. `2025-06-18` was the current MCP protocol revision at the time this was written.
const FALLBACK_PROTOCOL_VERSION: &str = "2025-06-18";

const EXIT_SUCCESS: i32 = 0;
const EXIT_IO_ERROR: i32 = 1;

/// Run the MCP stdio server loop on the process standard streams until EOF or a fatal I/O error.
pub(crate) fn run_stdio() -> i32 {
    let stdin = io::stdin();
    let stdout = io::stdout();
    run_session(stdin.lock(), stdout.lock())
}

/// One outcome of reading a single bounded, newline-delimited message from the transport.
#[derive(Debug, PartialEq, Eq)]
enum LineOutcome {
    /// The stream ended with no further message (EOF at a message boundary).
    Eof,
    /// A complete line was read and is valid UTF-8.
    Line(String),
    /// A line exceeded [`MAX_MCP_LINE_BYTES`] before a newline (or EOF) was found.
    TooLong,
    /// A complete line was read but was not valid UTF-8. JSON-RPC messages must be UTF-8.
    InvalidUtf8,
}

fn run_session<R: BufRead, W: Write>(mut reader: R, mut writer: W) -> i32 {
    // MCP lifecycle state: set once an `initialize` request succeeds (MBA-1337 m3).
    let mut initialized = false;
    loop {
        let outcome = match read_bounded_line(&mut reader) {
            Ok(outcome) => outcome,
            Err(error) => {
                // The transport itself failed (not merely malformed content). Report it and
                // stop; unlike a malformed message, a broken pipe will not recover on retry.
                let response = error_response(
                    Value::Null,
                    JSON_RPC_INTERNAL_ERROR,
                    &format!("failed to read from stdin: {error}"),
                    None,
                );
                let _ = write_message(&mut writer, &response);
                return EXIT_IO_ERROR;
            }
        };

        let response = match outcome {
            LineOutcome::Eof => return EXIT_SUCCESS,
            LineOutcome::TooLong => Some(error_response(
                Value::Null,
                JSON_RPC_PARSE_ERROR,
                &format!("message exceeds the {MAX_MCP_LINE_BYTES}-byte line limit"),
                None,
            )),
            LineOutcome::InvalidUtf8 => Some(error_response(
                Value::Null,
                JSON_RPC_PARSE_ERROR,
                "message is not valid UTF-8",
                None,
            )),
            LineOutcome::Line(line) => dispatch_line(&line, &mut initialized),
        };

        if let Some(response) = response {
            if !write_message(&mut writer, &response) {
                return EXIT_IO_ERROR;
            }
        }
    }
}

/// Read one bounded, newline-delimited message from `reader`.
///
/// This intentionally does not use [`BufRead::read_line`], which allocates without bound while
/// scanning for a newline. Instead it grows a buffer chunk by chunk via `fill_buf`/`consume`,
/// refusing to keep buffering once [`MAX_MCP_LINE_BYTES`] is exceeded — the transient overshoot
/// is bounded by one internal buffer chunk, not by the attacker-controlled line length.
fn read_bounded_line<R: BufRead>(reader: &mut R) -> io::Result<LineOutcome> {
    let mut buf: Vec<u8> = Vec::new();
    let mut too_long = false;

    loop {
        let mut newline_found = false;
        let consume_len;
        {
            let available = reader.fill_buf()?;
            if available.is_empty() {
                // Genuine EOF: the underlying reader has no more bytes to give us.
                break;
            }
            consume_len = match available.iter().position(|&byte| byte == b'\n') {
                Some(newline_at) => {
                    newline_found = true;
                    accumulate(&mut buf, &mut too_long, &available[..newline_at]);
                    newline_at + 1
                }
                None => {
                    accumulate(&mut buf, &mut too_long, available);
                    available.len()
                }
            };
        }
        reader.consume(consume_len);
        if newline_found {
            return Ok(finish_line(buf, too_long));
        }
    }

    if buf.is_empty() && !too_long {
        Ok(LineOutcome::Eof)
    } else {
        // A final message with no trailing newline before EOF is accepted rather than dropped.
        Ok(finish_line(buf, too_long))
    }
}

fn accumulate(buf: &mut Vec<u8>, too_long: &mut bool, chunk: &[u8]) {
    if *too_long {
        return;
    }
    if buf.len() + chunk.len() > MAX_MCP_LINE_BYTES {
        *too_long = true;
        buf.clear();
        buf.shrink_to_fit();
    } else {
        buf.extend_from_slice(chunk);
    }
}

fn finish_line(buf: Vec<u8>, too_long: bool) -> LineOutcome {
    if too_long {
        return LineOutcome::TooLong;
    }
    match String::from_utf8(buf) {
        Ok(text) => LineOutcome::Line(text),
        Err(_) => LineOutcome::InvalidUtf8,
    }
}

fn write_message<W: Write>(writer: &mut W, value: &Value) -> bool {
    let mut payload = serde_json::to_vec(value).unwrap_or_else(|_| {
        br#"{"jsonrpc":"2.0","id":null,"error":{"code":-32603,"message":"internal error"}}"#
            .to_vec()
    });
    payload.push(b'\n');
    writer.write_all(&payload).and_then(|()| writer.flush()).is_ok()
}

/// A JSON-RPC protocol-level error: an unknown method, malformed params, or a contained panic.
struct RpcError {
    code: i64,
    message: String,
    data: Option<Value>,
}

impl RpcError {
    fn new(code: i64, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            data: None,
        }
    }

    fn with_data(code: i64, message: impl Into<String>, data: Value) -> Self {
        Self {
            code,
            message: message.into(),
            data: Some(data),
        }
    }
}

fn success_response(id: Value, result: Value) -> Value {
    json!({"jsonrpc": "2.0", "id": id, "result": result})
}

fn error_response(id: Value, code: i64, message: &str, data: Option<Value>) -> Value {
    let mut error = json!({"code": code, "message": message});
    if let Some(data) = data {
        error["data"] = data;
    }
    json!({"jsonrpc": "2.0", "id": id, "error": error})
}

/// Parse and dispatch one line. Returns `None` exactly when the line was a well-formed JSON-RPC
/// notification (no `id` member): notifications never receive a response, per JSON-RPC 2.0,
/// even when the notified method is unrecognized. Any line that cannot be parsed as a
/// structurally valid JSON-RPC message (bad JSON, not an object, wrong `jsonrpc` version, a
/// missing or non-string `method`, or an `id` of a disallowed type) always gets a response, with
/// `id: null` when no valid id could be recovered — this is the sense in which "malformed lines
/// produce JSON-RPC errors, never process exit."
fn dispatch_line(line: &str, initialized: &mut bool) -> Option<Value> {
    let value: Value = match serde_json::from_str(line) {
        Ok(value) => value,
        Err(error) => {
            return Some(error_response(
                Value::Null,
                JSON_RPC_PARSE_ERROR,
                &format!("parse error: {error}"),
                None,
            ));
        }
    };

    let Some(object) = value.as_object() else {
        return Some(error_response(
            Value::Null,
            JSON_RPC_INVALID_REQUEST,
            "invalid request: expected a JSON object",
            None,
        ));
    };

    let id_present = object.contains_key("id");
    let id_value = object.get("id").cloned().unwrap_or(Value::Null);
    // JSON-RPC 2.0: a Number id SHOULD NOT contain fractional parts — reject 1.5 but
    // keep accepting integral-valued floats like 2.0, which the spec permits and which
    // this server always accepted (MBA-1337 m1).
    let id_type_ok = match &id_value {
        Value::String(_) | Value::Null => true,
        Value::Number(n) => {
            n.is_i64()
                || n.is_u64()
                || n.as_f64().is_some_and(|f| f.is_finite() && f.fract() == 0.0)
        }
        _ => false,
    };

    let jsonrpc_ok = object.get("jsonrpc").and_then(Value::as_str) == Some("2.0");
    let method = object.get("method").and_then(Value::as_str).map(str::to_string);

    if !jsonrpc_ok || method.is_none() || !id_type_ok {
        let response_id = if id_type_ok { id_value } else { Value::Null };
        return Some(error_response(
            response_id,
            JSON_RPC_INVALID_REQUEST,
            "invalid request: expected a JSON-RPC 2.0 request or notification",
            None,
        ));
    }

    let method = method.expect("checked above");
    let params = object.get("params").cloned();

    if !id_present {
        // A well-formed notification. `notifications/initialized` is accepted as a no-op, and so
        // is every other notification method: there is no response channel to report anything
        // on, so there is nothing further to do.
        return None;
    }

    // MCP lifecycle (MBA-1337 m3): until an `initialize` request has succeeded, only
    // `initialize` and `ping` are served. Notifications already returned above (no
    // response channel), so this gates requests only.
    if !*initialized && method != "initialize" && method != "ping" {
        return Some(error_response(
            id_value,
            JSON_RPC_SERVER_NOT_INITIALIZED,
            "server not initialized: send an initialize request first",
            None,
        ));
    }

    let response = dispatch_guarded(id_value, || handle_request(&method, params.as_ref()));
    if method == "initialize" && response.get("result").is_some() {
        *initialized = true;
    }
    Some(response)
}

/// Run `call` behind [`catch_unwind`] and turn its outcome into a response envelope.
///
/// A contained panic anywhere in request handling (including deep inside the engine) becomes an
/// internal-error response, never a dead server process — mirroring `solve_json_command`'s
/// precedent for the one-shot `solve-json` transport. `AssertUnwindSafe` is appropriate here for
/// the same reason it is there: `call` closes over process-local, per-message state that is
/// discarded after this one response is built.
fn dispatch_guarded<F>(id: Value, call: F) -> Value
where
    F: FnOnce() -> Result<Value, RpcError>,
{
    match catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(result)) => success_response(id, result),
        Ok(Err(rpc_error)) => error_response(id, rpc_error.code, &rpc_error.message, rpc_error.data),
        Err(_) => error_response(id, JSON_RPC_INTERNAL_ERROR, "internal error", None),
    }
}

fn handle_request(method: &str, params: Option<&Value>) -> Result<Value, RpcError> {
    match method {
        "initialize" => Ok(handle_initialize(params)),
        "ping" => Ok(json!({})),
        "tools/list" => Ok(handle_tools_list()),
        "tools/call" => handle_tools_call(params),
        _ => Err(RpcError::new(
            JSON_RPC_METHOD_NOT_FOUND,
            format!("method not found: {method}"),
        )),
    }
}

fn handle_initialize(params: Option<&Value>) -> Value {
    let protocol_version = params
        .and_then(|params| params.get("protocolVersion"))
        .and_then(Value::as_str)
        .unwrap_or(FALLBACK_PROTOCOL_VERSION);
    json!({
        "protocolVersion": protocol_version,
        "capabilities": {
            "tools": {}
        },
        "serverInfo": {
            "name": "ballistics-engine",
            "version": env!("CARGO_PKG_VERSION")
        }
    })
}

fn handle_tools_list() -> Value {
    json!({
        "tools": [solve_tool_definition(), engine_info_tool_definition()]
    })
}

fn solve_tool_definition() -> Value {
    json!({
        "name": "solve",
        "description": "Compute a full ballistics trajectory solve. Arguments ARE a solve-json \
            v1 request object (schema_version 1; SI units throughout — see \
            docs/SOLVE_JSON_V1.md in the ballistics-engine repository for the full field \
            reference). The result's text content is the solve-json v1 response JSON: a success \
            envelope with resolved_request/summary/samples when isError is false, or a \
            solve-json v1 error envelope when isError is true. Arguments that do not even form a \
            structurally valid solve-json v1 request (unknown or missing fields, wrong JSON \
            types, an unsupported schema_version) are instead rejected as a JSON-RPC \
            'Invalid params' protocol error.",
        "inputSchema": solve_input_schema(),
    })
}

fn engine_info_tool_definition() -> Value {
    json!({
        "name": "engine_info",
        "description": "Report static information about this ballistics-engine build: its crate \
            version, the drag models the solve-json v1 contract accepts, and the crate feature \
            flags this binary was compiled with. Takes no arguments.",
        "inputSchema": {
            "type": "object",
            "properties": {},
            "additionalProperties": false
        },
    })
}

/// A JSON Schema for [`ballistics_engine::solve_json::SolveRequestV1`], kept consistent with
/// that type's `deny_unknown_fields` shape and the field reference in `docs/SOLVE_JSON_V1.md`.
/// This schema is advisory for MCP clients; the "solve" tool itself always re-validates
/// arguments with the same `decode_solve_request_v1` entry point the `solve-json` command uses,
/// so this schema drifting would be a documentation problem, not a validation-bypass one.
fn solve_input_schema() -> Value {
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": [
            "schema_version", "projectile", "rifle", "shot", "atmosphere", "wind", "solver",
            "effects", "sampling"
        ],
        "properties": {
            "schema_version": {"type": "integer", "const": 1},
            "projectile": {
                "type": "object",
                "additionalProperties": false,
                "required": ["mass_kg", "diameter_m", "drag_model", "ballistic_coefficient"],
                "properties": {
                    "mass_kg": {"type": "number"},
                    "diameter_m": {"type": "number"},
                    "length_m": {"type": "number"},
                    "drag_model": {"type": "string", "enum": DRAG_MODEL_WIRE_NAMES_V1},
                    "ballistic_coefficient": {"type": "number"}
                }
            },
            "rifle": {
                "type": "object",
                "additionalProperties": false,
                "required": ["muzzle_velocity_mps"],
                "properties": {
                    "muzzle_velocity_mps": {"type": "number"},
                    "sight_height_m": {"type": "number"},
                    "muzzle_height_m": {"type": "number"},
                    "twist_rate_m_per_turn": {"type": "number"},
                    "twist_direction": {"type": "string", "enum": ["left", "right"]},
                    "sight_offset_lateral_m": {"type": "number"}
                }
            },
            "shot": {
                "type": "object",
                "additionalProperties": false,
                "required": ["max_range_m"],
                "properties": {
                    "max_range_m": {"type": "number"},
                    "zero_distance_m": {"type": "number"},
                    "muzzle_angle_rad": {"type": "number"},
                    "aim_azimuth_rad": {"type": "number"},
                    "shot_azimuth_rad": {"type": "number"},
                    "shooting_angle_rad": {"type": "number"},
                    "cant_angle_rad": {"type": "number"},
                    "target_height_m": {"type": "number"},
                    "ground_threshold_m": {"type": "number"},
                    "zero_poi_up_m": {"type": "number"},
                    "zero_poi_right_m": {"type": "number"},
                    "drops_reference": {"type": "string", "enum": ["los", "target"]}
                }
            },
            "atmosphere": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "altitude_m": {"type": "number"},
                    "temperature_k": {"type": "number"},
                    "pressure_pa": {"type": "number"},
                    "relative_humidity": {"type": "number"},
                    "latitude_rad": {"type": "number"}
                }
            },
            "wind": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "speed_mps": {"type": "number"},
                    "direction_from_rad": {"type": "number"},
                    "vertical_speed_mps": {"type": "number"},
                    "segments": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "additionalProperties": false,
                            "required": ["until_distance_m", "speed_mps", "direction_from_rad"],
                            "properties": {
                                "until_distance_m": {"type": "number"},
                                "speed_mps": {"type": "number"},
                                "direction_from_rad": {"type": "number"},
                                "vertical_speed_mps": {"type": "number"}
                            }
                        }
                    },
                    "wind_reference": {"type": "string", "enum": ["shooter", "compass"]}
                }
            },
            "solver": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "method": {"type": "string", "enum": ["euler", "rk4", "rk45"]},
                    "time_step_s": {"type": "number"}
                }
            },
            "effects": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "magnus": {"type": "boolean"},
                    "coriolis": {"type": "boolean"},
                    "enhanced_spin_drift": {"type": "boolean"}
                }
            },
            "sampling": {
                "type": "object",
                "additionalProperties": false,
                "properties": {
                    "interval_m": {"type": "number"}
                }
            },
            // Optional reticle hold request (MBA-1361), mirroring solve-json's
            // ReticleRequestV1. Without this the top-level additionalProperties:false above
            // rejects a `reticle` block that the solve-json CLI path accepts. The envelope
            // is strict; `description` stays permissive about extra keys, matching
            // ReticleRequestV1 (front ends carry render metadata inside the description).
            "reticle": {
                "type": "object",
                "additionalProperties": false,
                "required": ["range_m", "magnification", "description"],
                "properties": {
                    "range_m": {"type": "number"},
                    "magnification": {"type": "number"},
                    "description": {"type": "object"}
                }
            }
        }
    })
}

fn handle_tools_call(params: Option<&Value>) -> Result<Value, RpcError> {
    let params = params
        .ok_or_else(|| RpcError::new(JSON_RPC_INVALID_PARAMS, "tools/call requires params"))?;
    let name = params
        .get("name")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            RpcError::new(JSON_RPC_INVALID_PARAMS, "tools/call params.name must be a string")
        })?;
    let arguments = params
        .get("arguments")
        .cloned()
        .unwrap_or_else(|| Value::Object(serde_json::Map::new()));

    match name {
        "solve" => call_solve_tool(arguments),
        // engine_info advertises additionalProperties:false and "takes no arguments";
        // enforce it (MBA-1337 m2) instead of silently discarding whatever arrives.
        "engine_info" => match &arguments {
            Value::Object(map) if map.is_empty() => Ok(call_engine_info_tool()),
            Value::Object(map) => Err(RpcError::new(
                JSON_RPC_INVALID_PARAMS,
                format!(
                    "engine_info takes no arguments; unexpected argument(s): {}",
                    map.keys().cloned().collect::<Vec<_>>().join(", ")
                ),
            )),
            _ => Err(RpcError::new(
                JSON_RPC_INVALID_PARAMS,
                "engine_info arguments must be an object (or omitted)",
            )),
        },
        other => Err(RpcError::new(
            JSON_RPC_INVALID_PARAMS,
            format!("unknown tool: {other}"),
        )),
    }
}

fn call_solve_tool(arguments: Value) -> Result<Value, RpcError> {
    // Reuse solve-json v1's own request decoder rather than duplicating its validation: the
    // tools/call arguments object IS a solve-json v1 request object.
    let request_json = serde_json::to_string(&arguments).map_err(|error| {
        RpcError::new(
            JSON_RPC_INTERNAL_ERROR,
            format!("failed to re-encode tool arguments as JSON: {error}"),
        )
    })?;

    let request = decode_solve_request_v1(&request_json).map_err(|envelope| {
        let data = serde_json::to_value(&envelope.error).unwrap_or(Value::Null);
        RpcError::with_data(
            JSON_RPC_INVALID_PARAMS,
            "solve arguments failed solve-json v1 validation",
            data,
        )
    })?;

    let (is_error, payload) = match ballistics_engine::solve_v1(request) {
        Ok(success) => (false, serde_json::to_value(&success)),
        Err(envelope) => (true, serde_json::to_value(&envelope)),
    };
    let payload = payload.map_err(|error| {
        RpcError::new(
            JSON_RPC_INTERNAL_ERROR,
            format!("failed to encode solve result: {error}"),
        )
    })?;
    let text = serde_json::to_string(&payload).unwrap_or_else(|_| "null".to_string());

    Ok(json!({
        "content": [{"type": "text", "text": text}],
        "isError": is_error
    }))
}

fn call_engine_info_tool() -> Value {
    let mut features: Vec<&str> = Vec::new();
    if cfg!(feature = "online") {
        features.push("online");
    }
    if cfg!(feature = "pdf") {
        features.push("pdf");
    }
    if cfg!(feature = "profile-import") {
        features.push("profile-import");
    }
    if cfg!(feature = "jemalloc") {
        features.push("jemalloc");
    }
    if cfg!(feature = "mimalloc") {
        features.push("mimalloc");
    }
    if cfg!(feature = "validation") {
        features.push("validation");
    }

    let info = json!({
        "engine_version": env!("CARGO_PKG_VERSION"),
        "drag_models": DRAG_MODEL_WIRE_NAMES_V1,
        "features": features
    });
    let text = serde_json::to_string(&info).unwrap_or_else(|_| "null".to_string());

    json!({
        "content": [{"type": "text", "text": text}],
        "isError": false
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Dispatch one line in an already-initialized session — the common case for
    /// request tests. Lifecycle-specific tests drive dispatch_line directly.
    fn dispatch(line: &str) -> Value {
        let mut initialized = true;
        dispatch_line(line, &mut initialized).expect("expected a response for this line")
    }

    #[test]
    fn read_bounded_line_reads_a_single_terminated_line() {
        let mut input: &[u8] = b"{\"a\":1}\n{\"b\":2}\n";
        let first = read_bounded_line(&mut input).expect("read first line");
        assert_eq!(first, LineOutcome::Line("{\"a\":1}".to_string()));
        let second = read_bounded_line(&mut input).expect("read second line");
        assert_eq!(second, LineOutcome::Line("{\"b\":2}".to_string()));
        let eof = read_bounded_line(&mut input).expect("read eof");
        assert_eq!(eof, LineOutcome::Eof);
    }

    #[test]
    fn read_bounded_line_accepts_a_final_unterminated_line() {
        let mut input: &[u8] = b"{\"a\":1}";
        let line = read_bounded_line(&mut input).expect("read final line");
        assert_eq!(line, LineOutcome::Line("{\"a\":1}".to_string()));
        let eof = read_bounded_line(&mut input).expect("read eof");
        assert_eq!(eof, LineOutcome::Eof);
    }

    #[test]
    fn read_bounded_line_rejects_an_oversized_line_without_unbounded_buffering() {
        let mut oversized = vec![b'a'; MAX_MCP_LINE_BYTES + 1];
        oversized.push(b'\n');
        oversized.extend_from_slice(b"{\"next\":true}\n");
        let mut input: &[u8] = &oversized;
        let outcome = read_bounded_line(&mut input).expect("read oversized line");
        assert_eq!(outcome, LineOutcome::TooLong);
        // The stream must resynchronize at the next newline: the following message reads clean.
        let next = read_bounded_line(&mut input).expect("read next line");
        assert_eq!(next, LineOutcome::Line("{\"next\":true}".to_string()));
    }

    #[test]
    fn read_bounded_line_accepts_exactly_the_limit() {
        let mut exact = vec![b'a'; MAX_MCP_LINE_BYTES];
        exact.push(b'\n');
        let mut input: &[u8] = &exact;
        let outcome = read_bounded_line(&mut input).expect("read exact-limit line");
        match outcome {
            LineOutcome::Line(text) => assert_eq!(text.len(), MAX_MCP_LINE_BYTES),
            other => panic!("expected an exact-limit line to be accepted, got {other:?}"),
        }
    }

    #[test]
    fn read_bounded_line_reports_invalid_utf8() {
        let mut input: &[u8] = &[0xff, 0xfe, b'\n'];
        let outcome = read_bounded_line(&mut input).expect("read invalid utf-8 line");
        assert_eq!(outcome, LineOutcome::InvalidUtf8);
    }

    #[test]
    fn malformed_json_is_a_parse_error() {
        let response = dispatch("not json");
        assert_eq!(response["error"]["code"], JSON_RPC_PARSE_ERROR);
        assert_eq!(response["id"], Value::Null);
    }

    #[test]
    fn non_object_message_is_an_invalid_request() {
        let response = dispatch("[1,2,3]");
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_REQUEST);
    }

    #[test]
    fn wrong_jsonrpc_version_is_an_invalid_request() {
        let response = dispatch(r#"{"jsonrpc":"1.0","id":1,"method":"ping"}"#);
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_REQUEST);
        assert_eq!(response["id"], json!(1));
    }

    #[test]
    fn well_formed_notification_gets_no_response() {
        let mut initialized = true;
        assert_eq!(
            dispatch_line(
                r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#,
                &mut initialized
            ),
            None
        );
    }

    #[test]
    fn fractional_numeric_id_is_invalid_request() {
        // JSON-RPC 2.0: Number ids SHOULD NOT contain fractional parts (MBA-1337 m1).
        let response = dispatch(r#"{"jsonrpc":"2.0","id":1.5,"method":"ping"}"#);
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_REQUEST);
        assert_eq!(response["id"], Value::Null);
    }

    #[test]
    fn integral_valued_float_id_stays_accepted() {
        // 2.0 has no fractional part: the spec permits it and this server always
        // accepted it — the m1 tightening must not regress that.
        let response = dispatch(r#"{"jsonrpc":"2.0","id":2.0,"method":"ping"}"#);
        assert_eq!(response["result"], json!({}));
        assert_eq!(response["id"], json!(2.0));
    }

    #[test]
    fn engine_info_rejects_unexpected_arguments() {
        // The schema advertises additionalProperties:false — enforce it (MBA-1337 m2).
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"engine_info","arguments":{"foo":1}}}"#,
        );
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_PARAMS);
        assert!(response["error"]["message"].as_str().unwrap().contains("foo"));
    }

    #[test]
    fn engine_info_rejects_non_object_arguments() {
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"engine_info","arguments":42}}"#,
        );
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_PARAMS);
    }

    #[test]
    fn requests_before_initialize_are_rejected_except_ping() {
        // MCP lifecycle (MBA-1337 m3): a fresh session serves only initialize + ping.
        let mut initialized = false;
        let rejected = dispatch_line(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#,
            &mut initialized,
        )
        .unwrap();
        assert_eq!(rejected["error"]["code"], JSON_RPC_SERVER_NOT_INITIALIZED);
        assert!(!initialized);

        let ping = dispatch_line(
            r#"{"jsonrpc":"2.0","id":2,"method":"ping"}"#,
            &mut initialized,
        )
        .unwrap();
        assert_eq!(ping["result"], json!({}));
        assert!(!initialized, "ping must not count as initialization");

        let init = dispatch_line(
            r#"{"jsonrpc":"2.0","id":3,"method":"initialize","params":{}}"#,
            &mut initialized,
        )
        .unwrap();
        assert!(init["result"].is_object());
        assert!(initialized);

        let listed = dispatch_line(
            r#"{"jsonrpc":"2.0","id":4,"method":"tools/list"}"#,
            &mut initialized,
        )
        .unwrap();
        assert!(listed["result"]["tools"].is_array());
    }

    #[test]
    fn unknown_method_is_method_not_found() {
        let response = dispatch(r#"{"jsonrpc":"2.0","id":7,"method":"nope"}"#);
        assert_eq!(response["error"]["code"], JSON_RPC_METHOD_NOT_FOUND);
        assert_eq!(response["id"], json!(7));
    }

    #[test]
    fn ping_returns_an_empty_result() {
        let response = dispatch(r#"{"jsonrpc":"2.0","id":1,"method":"ping"}"#);
        assert_eq!(response["result"], json!({}));
    }

    #[test]
    fn initialize_echoes_the_requested_protocol_version() {
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2099-01-01"}}"#,
        );
        assert_eq!(response["result"]["protocolVersion"], "2099-01-01");
        assert_eq!(response["result"]["serverInfo"]["name"], "ballistics-engine");
        assert_eq!(
            response["result"]["serverInfo"]["version"],
            env!("CARGO_PKG_VERSION")
        );
        assert!(response["result"]["capabilities"]["tools"].is_object());
    }

    #[test]
    fn initialize_without_a_protocol_version_falls_back() {
        let response = dispatch(r#"{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}"#);
        assert_eq!(
            response["result"]["protocolVersion"],
            FALLBACK_PROTOCOL_VERSION
        );
    }

    #[test]
    fn tools_list_reports_solve_and_engine_info() {
        let response = dispatch(r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#);
        let tools = response["result"]["tools"].as_array().expect("tools array");
        let names: Vec<&str> = tools
            .iter()
            .map(|tool| tool["name"].as_str().expect("tool name"))
            .collect();
        assert_eq!(names, vec!["solve", "engine_info"]);
        for tool in tools {
            assert_eq!(tool["inputSchema"]["type"], "object");
        }
        assert_eq!(
            tools[0]["inputSchema"]["properties"]["projectile"]["properties"]["drag_model"]
                ["enum"],
            json!(DRAG_MODEL_WIRE_NAMES_V1)
        );
    }

    #[test]
    fn tools_call_unknown_tool_is_invalid_params() {
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"nope"}}"#,
        );
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_PARAMS);
    }

    #[test]
    fn tools_call_engine_info_reports_version_and_drag_models() {
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"engine_info"}}"#,
        );
        assert_eq!(response["result"]["isError"], false);
        let text = response["result"]["content"][0]["text"]
            .as_str()
            .expect("text content");
        let info: Value = serde_json::from_str(text).expect("engine_info text is JSON");
        assert_eq!(info["engine_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(info["drag_models"], json!(DRAG_MODEL_WIRE_NAMES_V1));
        assert!(info["features"].is_array());
    }

    #[test]
    fn tools_call_solve_with_malformed_arguments_is_invalid_params() {
        let response = dispatch(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"solve","arguments":{}}}"#,
        );
        assert_eq!(response["error"]["code"], JSON_RPC_INVALID_PARAMS);
        assert!(response["error"]["data"]["code"].is_string());
    }

    #[test]
    fn panics_in_request_handling_become_internal_errors_not_process_death() {
        let response = dispatch_guarded(json!(9), || panic!("injected panic payload"));
        assert_eq!(response["error"]["code"], JSON_RPC_INTERNAL_ERROR);
        assert_eq!(response["id"], json!(9));
        let serialized = response.to_string();
        assert!(!serialized.contains("injected panic payload"));
    }
}
