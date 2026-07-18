use serde_json::{json, Value};
use std::io::Write;
use std::process::{Command, Output, Stdio};

const MAX_INPUT_BYTES: usize = 1024 * 1024;

fn request_value() -> Value {
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

fn encode(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).expect("request fixture must serialize")
}

fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn ballistics");
    child
        .stdin
        .take()
        .expect("piped stdin")
        .write_all(input)
        .expect("write request");
    child.wait_with_output().expect("wait for ballistics")
}

fn envelope(output: &Output) -> Value {
    assert!(
        output.stdout.ends_with(b"\n"),
        "envelope needs one terminator"
    );
    assert_eq!(
        output.stdout.iter().filter(|&&byte| byte == b'\n').count(),
        1,
        "stdout must contain exactly one compact JSON envelope: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    serde_json::from_slice(&output.stdout).expect("stdout must be one JSON envelope")
}

fn assert_exit(output: &Output, expected: i32) {
    assert_eq!(
        output.status.code(),
        Some(expected),
        "unexpected exit status; stderr={} stdout={}",
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
}

#[test]
fn success_is_one_clean_envelope_and_units_are_ignored() {
    let mut request = request_value();
    // Euler historically emitted debug lines; the process protocol must stay clean for every
    // solver method.
    request["solver"]["method"] = json!("euler");
    let input = encode(&request);
    let metric = run(&["--units", "metric", "solve-json"], &input);
    let imperial = run(&["--units", "imperial", "solve-json"], &input);

    for output in [&metric, &imperial] {
        assert_exit(output, 0);
        assert!(
            output.stderr.is_empty(),
            "stderr must stay quiet on success"
        );
        let wire = envelope(output);
        assert_eq!(wire["schema_version"], 1);
        assert_eq!(wire["status"], "ok");
        assert!(wire["samples"]
            .as_array()
            .is_some_and(|rows| !rows.is_empty()));
    }
    assert_eq!(
        metric.stdout, imperial.stdout,
        "the explicit-SI protocol cannot depend on global --units"
    );
}

#[test]
fn named_conformance_fixture_requests_stay_clean_over_the_cli_transport() {
    for (name, request, success) in [
        (
            "calm-air",
            include_str!("fixtures/solve_json_v1/calm-air.request.json"),
            include_str!("fixtures/solve_json_v1/calm-air.success.json"),
        ),
        (
            "crosswind",
            include_str!("fixtures/solve_json_v1/crosswind.request.json"),
            include_str!("fixtures/solve_json_v1/crosswind.success.json"),
        ),
        (
            "segmented-wind",
            include_str!("fixtures/solve_json_v1/segmented-wind.request.json"),
            include_str!("fixtures/solve_json_v1/segmented-wind.success.json"),
        ),
        (
            "zeroed",
            include_str!("fixtures/solve_json_v1/zeroed.request.json"),
            include_str!("fixtures/solve_json_v1/zeroed.success.json"),
        ),
    ] {
        let output = run(&["solve-json"], request.as_bytes());
        assert_exit(&output, 0);
        assert!(
            output.stderr.is_empty(),
            "{name}: stderr must stay quiet on success"
        );
        let wire = envelope(&output);
        let expected: Value =
            serde_json::from_str(success).expect("checked success fixture must be valid JSON");
        assert_eq!(wire["schema_version"], 1, "{name}");
        assert_eq!(wire["status"], "ok", "{name}");
        assert_eq!(
            wire["summary"]["termination"], expected["summary"]["termination"],
            "{name}"
        );
        assert_eq!(
            wire["samples"].as_array().map(Vec::len),
            expected["samples"].as_array().map(Vec::len),
            "{name}: sample cardinality diverged from the checked fixture"
        );
        // The resolved request is a deterministic echo of the input plus documented defaults for
        // every pair except `zeroed`, whose effective muzzle angle is solver-derived and is pinned
        // with numeric tolerance by the contract tests instead.
        if name != "zeroed" {
            assert_eq!(
                wire["resolved_request"], expected["resolved_request"],
                "{name}: resolved request diverged from the checked fixture"
            );
        }
    }
}

#[test]
fn malformed_json_and_invalid_utf8_report_source_locations() {
    for input in [
        b"{\n  \"schema_version\": 1,\n  \"projectile\": ]\n}".as_slice(),
        b"{\n  \"schema_version\": \xff\n}".as_slice(),
    ] {
        let output = run(&["solve-json"], input);
        assert_exit(&output, 2);
        assert!(output.stderr.is_empty());
        let wire = envelope(&output);
        assert_eq!(wire["error"]["code"], "invalid_json");
        assert!(wire["error"]["line"].as_u64().is_some_and(|line| line > 0));
        assert!(wire["error"]["column"]
            .as_u64()
            .is_some_and(|column| column > 0));
        assert!(wire["error"]["path"].is_null());
    }
}

#[test]
fn schema_shape_and_semantic_failures_use_request_exit_status() {
    let mut cases = Vec::new();

    let mut unsupported = request_value();
    unsupported["schema_version"] = json!(2);
    cases.push((
        unsupported,
        "unsupported_schema_version",
        "$.schema_version",
    ));

    let mut unknown = request_value();
    unknown["unexpected"] = json!(true);
    cases.push((unknown, "unknown_field", "$.unexpected"));

    let mut missing = request_value();
    missing.as_object_mut().expect("object").remove("sampling");
    cases.push((missing, "missing_field", "$.sampling"));

    let mut invalid = request_value();
    invalid["projectile"]["mass_kg"] = json!(-1.0);
    cases.push((invalid, "invalid_value", "$.projectile.mass_kg"));

    let mut conflicting = request_value();
    conflicting["shot"]["zero_distance_m"] = json!(25.0);
    cases.push((conflicting, "conflicting_fields", "$.shot"));

    for (request, code, path) in cases {
        let output = run(&["solve-json"], &encode(&request));
        assert_exit(&output, 2);
        assert!(output.stderr.is_empty());
        let wire = envelope(&output);
        assert_eq!(wire["error"]["code"], code);
        assert_eq!(wire["error"]["path"], path);
    }
}

#[test]
fn concatenated_requests_are_rejected_as_one_malformed_exchange() {
    let request = encode(&request_value());
    let mut input = request.clone();
    input.extend_from_slice(&request);
    let output = run(&["solve-json"], &input);
    assert_exit(&output, 2);
    assert!(output.stderr.is_empty());
    assert_eq!(envelope(&output)["error"]["code"], "invalid_json");
}

#[test]
fn input_limit_accepts_exactly_one_mib_and_rejects_one_byte_more() {
    let mut exact = encode(&request_value());
    assert!(exact.len() < MAX_INPUT_BYTES);
    exact.resize(MAX_INPUT_BYTES, b' ');

    let accepted = run(&["solve-json"], &exact);
    assert_exit(&accepted, 0);
    assert_eq!(envelope(&accepted)["status"], "ok");

    exact.push(b' ');
    let rejected = run(&["solve-json"], &exact);
    assert_exit(&rejected, 3);
    assert!(rejected.stderr.is_empty());
    assert_eq!(envelope(&rejected)["error"]["code"], "resource_limit");
}

#[test]
fn service_sample_limit_uses_resource_exit_status() {
    let mut request = request_value();
    request["sampling"]["interval_m"] = json!(0.001);
    let output = run(&["solve-json"], &encode(&request));
    assert_exit(&output, 3);
    assert!(output.stderr.is_empty());
    let wire = envelope(&output);
    assert_eq!(wire["error"]["code"], "resource_limit");
    assert_eq!(wire["error"]["path"], "$.sampling.interval_m");
}

#[test]
fn engine_failure_uses_solve_exit_status() {
    let mut request = request_value();
    request["rifle"]["muzzle_velocity_mps"] = json!(1.0e308);
    let output = run(&["solve-json"], &encode(&request));
    assert_exit(&output, 3);
    assert!(output.stderr.is_empty());
    assert_eq!(envelope(&output)["error"]["code"], "solve_failed");
}

#[cfg(unix)]
#[test]
fn stdin_read_failure_is_one_io_error_envelope() {
    let (reader, writer) =
        std::os::unix::net::UnixStream::pair().expect("create stdin socket pair");
    reader
        .set_nonblocking(true)
        .expect("make stdin socket nonblocking");
    let reader: std::os::fd::OwnedFd = reader.into();
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .arg("solve-json")
        .stdin(Stdio::from(reader))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("run ballistics with unreadable stdin");
    drop(writer);
    assert_exit(&output, 1);
    assert!(output.stderr.is_empty());
    assert_eq!(envelope(&output)["error"]["code"], "io_error");
}

#[test]
fn command_help_does_not_read_the_protocol_stream() {
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(["solve-json", "--help"])
        .output()
        .expect("run solve-json help");
    assert_exit(&output, 0);
    assert!(output.stdout.starts_with(b"Solve one explicit-SI"));
}
