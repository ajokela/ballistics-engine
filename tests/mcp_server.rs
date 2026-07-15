//! Integration tests for the `ballistics mcp` command: a Model Context Protocol server over the
//! stdio JSON-RPC 2.0 transport (MBA-1319).
//!
//! These tests spawn the real `ballistics` binary (the `CARGO_BIN_EXE_ballistics` pattern used
//! by the sibling `solve-json` transport tests in `tests/solve_json_cli.rs`) and drive a full
//! session over its stdin/stdout, rather than calling the server loop as a library function.

use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

/// Must match `MAX_MCP_LINE_BYTES` in `src/mcp_command.rs`.
const MAX_MCP_LINE_BYTES: usize = 1024 * 1024;

fn solve_request_arguments() -> Value {
    json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 0.01134,
            "diameter_m": 0.00782,
            "length_m": 0.031,
            "drag_model": "G7",
            "ballistic_coefficient": 0.243
        },
        "rifle": {
            "muzzle_velocity_mps": 823.0,
            "sight_height_m": 0.05,
            "muzzle_height_m": 1.2,
            "twist_rate_m_per_turn": 0.2032,
            "twist_direction": "right"
        },
        "shot": {
            "max_range_m": 50.0,
            "muzzle_angle_rad": 0.001,
            "aim_azimuth_rad": 0.0,
            "shot_azimuth_rad": 0.0,
            "shooting_angle_rad": 0.0,
            "cant_angle_rad": 0.0,
            "target_height_m": 1.25,
            "ground_threshold_m": -100.0
        },
        "atmosphere": {
            "altitude_m": 0.0,
            "temperature_k": 293.15,
            "pressure_pa": 100000.0,
            "relative_humidity": 0.4
        },
        "wind": {
            "speed_mps": 0.0,
            "direction_from_rad": 0.0,
            "vertical_speed_mps": 0.0
        },
        "solver": {
            "method": "rk45",
            "time_step_s": 0.001
        },
        "effects": {
            "magnus": false,
            "coriolis": false,
            "enhanced_spin_drift": false
        },
        "sampling": {
            "interval_m": 10.0
        }
    })
}

/// A live `ballistics mcp` subprocess with its stdio pipes wired for line-at-a-time driving.
struct McpSession {
    child: Child,
    stdin: Option<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl McpSession {
    fn spawn() -> Self {
        let mut child = Command::new(env!("CARGO_BIN_EXE_ballistics"))
            .arg("mcp")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("spawn ballistics mcp");
        let stdin = child.stdin.take().expect("piped stdin");
        let stdout = BufReader::new(child.stdout.take().expect("piped stdout"));
        Self {
            child,
            stdin: Some(stdin),
            stdout,
        }
    }

    fn send(&mut self, message: &Value) {
        let mut line = serde_json::to_vec(message).expect("serialize JSON-RPC message");
        line.push(b'\n');
        self.send_raw(&line);
    }

    fn send_raw(&mut self, raw: &[u8]) {
        let stdin = self.stdin.as_mut().expect("stdin still open");
        stdin.write_all(raw).expect("write to ballistics mcp stdin");
        stdin.flush().expect("flush ballistics mcp stdin");
    }

    /// Read exactly one JSON-RPC response line. Panics on EOF: callers expect a response here.
    fn recv(&mut self) -> Value {
        let mut line = String::new();
        let bytes_read = self
            .stdout
            .read_line(&mut line)
            .expect("read a response line from ballistics mcp stdout");
        assert!(
            bytes_read > 0,
            "expected a JSON-RPC response line but the server closed stdout"
        );
        serde_json::from_str(line.trim_end_matches('\n'))
            .unwrap_or_else(|error| panic!("response line was not JSON ({error}): {line:?}"))
    }

    /// Close stdin (signaling EOF to the server) and wait for the process to exit.
    fn close_and_wait(mut self) -> std::process::ExitStatus {
        self.stdin = None; // drop the ChildStdin, closing the pipe
        let status = self.child.wait().expect("wait for ballistics mcp to exit");
        // Read stderr fully so an unexpectedly noisy server surfaces in the failure output.
        if let Some(mut stderr) = self.child.stderr.take() {
            use std::io::Read;
            let mut buffer = String::new();
            let _ = stderr.read_to_string(&mut buffer);
            assert!(
                buffer.is_empty(),
                "ballistics mcp wrote unexpected stderr output: {buffer}"
            );
        }
        status
    }
}

fn initialize(session: &mut McpSession, id: i64) -> Value {
    session.send(&json!({
        "jsonrpc": "2.0",
        "id": id,
        "method": "initialize",
        "params": {
            "protocolVersion": "2025-06-18",
            "capabilities": {},
            "clientInfo": {"name": "mcp-server-test", "version": "0.0.0"}
        }
    }));
    session.recv()
}

fn send_initialized_notification(session: &mut McpSession) {
    session.send(&json!({
        "jsonrpc": "2.0",
        "method": "notifications/initialized"
    }));
}

#[test]
fn full_session_initialize_list_call_and_clean_eof() {
    let mut session = McpSession::spawn();

    // 1. initialize: protocolVersion is echoed back, serverInfo identifies the engine.
    let initialize_response = initialize(&mut session, 1);
    assert_eq!(initialize_response["jsonrpc"], "2.0");
    assert_eq!(initialize_response["id"], 1);
    assert_eq!(
        initialize_response["result"]["protocolVersion"],
        "2025-06-18"
    );
    assert_eq!(
        initialize_response["result"]["serverInfo"]["name"],
        "ballistics-engine"
    );
    assert!(initialize_response["result"]["capabilities"]["tools"].is_object());

    // 2. notifications/initialized: accepted, no-op, and — being a notification — produces no
    // response line at all. We prove that by immediately sending a request and observing that
    // its response is the very next line on stdout.
    send_initialized_notification(&mut session);

    // 3. tools/list: both tools are present with an object-typed inputSchema.
    session.send(&json!({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}));
    let tools_response = session.recv();
    assert_eq!(tools_response["id"], 2);
    let tools = tools_response["result"]["tools"]
        .as_array()
        .expect("tools/list result.tools must be an array");
    let names: Vec<&str> = tools
        .iter()
        .map(|tool| tool["name"].as_str().expect("tool.name is a string"))
        .collect();
    assert_eq!(names, vec!["solve", "engine_info"]);
    for tool in tools {
        assert_eq!(
            tool["inputSchema"]["type"], "object",
            "tool {:?} must publish an object JSON Schema",
            tool["name"]
        );
        assert!(tool["description"].as_str().is_some_and(|d| !d.is_empty()));
    }
    let solve_schema = &tools[0]["inputSchema"];
    let required = solve_schema["required"]
        .as_array()
        .expect("solve inputSchema must list required top-level fields");
    for field in [
        "schema_version",
        "projectile",
        "rifle",
        "shot",
        "atmosphere",
        "wind",
        "solver",
        "effects",
        "sampling",
    ] {
        assert!(
            required.contains(&json!(field)),
            "solve inputSchema.required is missing {field:?}"
        );
    }

    // 4. tools/call "solve" with a valid solve-json v1 request: a real trajectory comes back.
    session.send(&json!({
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {"name": "solve", "arguments": solve_request_arguments()}
    }));
    let solve_response = session.recv();
    assert_eq!(solve_response["id"], 3);
    assert_eq!(solve_response["result"]["isError"], false);
    let content_text = solve_response["result"]["content"][0]["text"]
        .as_str()
        .expect("tool result content[0].text must be a string");
    let solve_payload: Value =
        serde_json::from_str(content_text).expect("solve tool text content must be JSON");
    assert_eq!(solve_payload["schema_version"], 1);
    assert_eq!(solve_payload["status"], "ok");
    assert!(solve_payload["summary"]["actual_range_m"]
        .as_f64()
        .is_some_and(|range| range > 0.0));
    let samples = solve_payload["samples"]
        .as_array()
        .expect("a real trajectory result must include samples");
    assert!(!samples.is_empty());

    // 5. tools/call "engine_info": no arguments, reports the crate version and drag models.
    session.send(&json!({
        "jsonrpc": "2.0",
        "id": 4,
        "method": "tools/call",
        "params": {"name": "engine_info", "arguments": {}}
    }));
    let engine_info_response = session.recv();
    assert_eq!(engine_info_response["result"]["isError"], false);
    let engine_info_text = engine_info_response["result"]["content"][0]["text"]
        .as_str()
        .expect("engine_info text content");
    let engine_info: Value =
        serde_json::from_str(engine_info_text).expect("engine_info text content must be JSON");
    assert_eq!(engine_info["engine_version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(engine_info["drag_models"], json!(["G1", "G6", "G7", "G8"]));

    // 6. tools/call "solve" with structurally invalid arguments: rejected as Invalid params
    // (-32602), per this server's documented split between malformed-arguments protocol errors
    // and well-formed-but-unsolvable tool-execution errors.
    session.send(&json!({
        "jsonrpc": "2.0",
        "id": 5,
        "method": "tools/call",
        "params": {"name": "solve", "arguments": {"schema_version": 1}}
    }));
    let invalid_params_response = session.recv();
    assert_eq!(invalid_params_response["id"], 5);
    assert!(invalid_params_response.get("result").is_none());
    assert_eq!(invalid_params_response["error"]["code"], -32602);
    assert_eq!(
        invalid_params_response["error"]["data"]["code"],
        "missing_field"
    );

    // 7. An unknown top-level method: Method not found (-32601).
    session.send(&json!({"jsonrpc": "2.0", "id": 6, "method": "not/a/real/method"}));
    let unknown_method_response = session.recv();
    assert_eq!(unknown_method_response["id"], 6);
    assert_eq!(unknown_method_response["error"]["code"], -32601);

    // 8. EOF: closing stdin ends the session cleanly with exit code 0, no stderr noise.
    let status = session.close_and_wait();
    assert_eq!(status.code(), Some(0));
}

#[test]
fn tools_call_solve_execution_failure_is_a_tool_error_result_not_a_protocol_error() {
    // A structurally valid solve-json v1 request that the engine cannot actually solve
    // (an absurd muzzle velocity) surfaces as isError: true in a normal JSON-RPC *result*,
    // not as a JSON-RPC protocol error — this is the other half of this server's documented
    // invalid-params-vs-tool-error split.
    let mut session = McpSession::spawn();
    let _ = initialize(&mut session, 1);
    send_initialized_notification(&mut session);

    let mut arguments = solve_request_arguments();
    arguments["rifle"]["muzzle_velocity_mps"] = json!(1.0e308);
    session.send(&json!({
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/call",
        "params": {"name": "solve", "arguments": arguments}
    }));
    let response = session.recv();
    assert!(response.get("error").is_none());
    assert_eq!(response["result"]["isError"], true);
    let text = response["result"]["content"][0]["text"]
        .as_str()
        .expect("tool error content is a text block");
    let error_payload: Value = serde_json::from_str(text).expect("error content is JSON");
    assert_eq!(error_payload["status"], "error");
    assert_eq!(error_payload["error"]["code"], "solve_failed");

    let status = session.close_and_wait();
    assert_eq!(status.code(), Some(0));
}

#[test]
fn malformed_json_line_gets_a_parse_error_and_the_session_keeps_going() {
    let mut session = McpSession::spawn();

    session.send_raw(b"this is not json at all\n");
    let parse_error_response = session.recv();
    assert!(parse_error_response.get("result").is_none());
    assert_eq!(parse_error_response["error"]["code"], -32700);
    assert_eq!(parse_error_response["id"], Value::Null);

    // The malformed line must not have killed the server: a subsequent well-formed request still
    // gets a normal response.
    session.send(&json!({"jsonrpc": "2.0", "id": 1, "method": "ping"}));
    let ping_response = session.recv();
    assert_eq!(ping_response["id"], 1);
    assert_eq!(ping_response["result"], json!({}));

    let status = session.close_and_wait();
    assert_eq!(status.code(), Some(0));
}

#[test]
fn an_oversized_line_is_rejected_without_killing_the_session() {
    let mut session = McpSession::spawn();

    let mut oversized = vec![b' '; MAX_MCP_LINE_BYTES + 1];
    oversized.push(b'\n');
    session.send_raw(&oversized);
    let response = session.recv();
    assert_eq!(response["error"]["code"], -32700);
    assert_eq!(response["id"], Value::Null);

    session.send(&json!({"jsonrpc": "2.0", "id": 1, "method": "ping"}));
    let ping_response = session.recv();
    assert_eq!(ping_response["result"], json!({}));

    let status = session.close_and_wait();
    assert_eq!(status.code(), Some(0));
}

#[test]
fn eof_with_no_input_at_all_exits_zero_cleanly() {
    let session = McpSession::spawn();
    let status = session.close_and_wait();
    assert_eq!(status.code(), Some(0));
}
