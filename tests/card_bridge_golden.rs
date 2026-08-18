//! Golden cross-check: the bridge's card commands must agree row-for-row with the
//! CLI card commands they mirror. The CLI is the reference implementation; any
//! divergence here means `card_service` drifted from `main.rs` and one of them
//! is lying to shooters.
//!
//! Requires both the `bridge` and `cli` features (the default set), because it
//! executes the actual `ballistics` binary via CARGO_BIN_EXE.

#![cfg(all(feature = "bridge", feature = "cli"))]

use serde_json::{json, Value};
use std::process::Command;

const TOL: f64 = 1e-9;

fn cli(args: &[&str]) -> Value {
    let out = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(args)
        .output()
        .expect("run ballistics CLI");
    assert!(
        out.status.success(),
        "CLI failed: {}\n{}",
        String::from_utf8_lossy(&out.stderr),
        String::from_utf8_lossy(&out.stdout)
    );
    serde_json::from_slice(&out.stdout).expect("CLI -o json output must parse")
}

fn bridge(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    let raw = ballistics_engine::bridge::bridge_call(&envelope.to_string());
    let out: Value = serde_json::from_str(&raw).expect("bridge output parses");
    assert_eq!(out["ok"], true, "bridge {command} failed: {out}");
    out["result"].clone()
}

fn assert_close(a: &Value, b: &Value, label: &str) {
    match (a.as_f64(), b.as_f64()) {
        (Some(x), Some(y)) => {
            let scale = x.abs().max(y.abs()).max(1.0);
            assert!(
                (x - y).abs() <= TOL * scale,
                "{label}: CLI {x} vs bridge {y}"
            );
        }
        _ => assert_eq!(a, b, "{label}: non-numeric mismatch"),
    }
}

/// The shared fixture load: a 175gr .308 at 2600 fps, G7 0.243, 100 yd zero,
/// full-value 10 mph wind, 100-600 yd in 100s, MIL adjustments.
fn fixture_request() -> Value {
    json!({
        "units": "imperial",
        "muzzle_velocity": 2600.0,
        "ballistic_coefficient": 0.243,
        "mass": 175.0,
        "diameter": 0.308,
        "drag_model": "g7",
        "sight_height": 1.5,
        "zero_distance": 100.0,
        "wind_speed": 10.0,
        "wind_direction_deg": 90.0,
        "start": 100.0,
        "end": 600.0,
        "step": 100.0,
        "adjustment_unit": "mil"
    })
}

fn fixture_cli_args<'a>(extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = vec![
        "-v", "2600", "-b", "0.243", "-m", "175", "-d", "0.308",
        "--drag-model", "g7",
        "--sight-height", "1.5",
        "--wind-speed", "10", "--wind-direction", "90",
        "--start", "100", "--end", "600", "--step", "100",
        "--adjustment-unit", "mil",
        "-o", "json",
    ];
    args.extend_from_slice(extra);
    args
}

#[test]
fn come_ups_matches_cli() {
    let cli_out = cli(
        &[&["come-ups", "--zero-distance", "100"][..], &fixture_cli_args(&[])[..]].concat(),
    );
    let bridge_out = bridge("card.come_ups", fixture_request());

    let cli_rows = cli_out["data"].as_array().expect("CLI data rows");
    let bridge_rows = bridge_out["rows"].as_array().expect("bridge rows");
    assert_eq!(cli_rows.len(), bridge_rows.len(), "row count");
    assert!(!cli_rows.is_empty(), "fixture must produce rows");

    for (c, b) in cli_rows.iter().zip(bridge_rows) {
        assert_close(&c["range"], &b["range"], "range");
        assert_close(&c["drop"], &b["drop_adj"], "drop_adj");
        assert_close(&c["come_up"], &b["come_up"], "come_up");
        assert_close(&c["velocity"], &b["velocity"], "velocity");
        assert_close(&c["energy"], &b["energy"], "energy");
        assert_close(&c["time"], &b["time"], "time");
    }
}

#[test]
fn range_table_matches_cli() {
    let cli_out = cli(
        &[&["range-table", "--zero-distance", "100"][..], &fixture_cli_args(&[])[..]].concat(),
    );
    let bridge_out = bridge("card.range_table", fixture_request());

    let cli_rows = cli_out["data"].as_array().expect("CLI data rows");
    let bridge_rows = bridge_out["rows"].as_array().expect("bridge rows");
    assert_eq!(cli_rows.len(), bridge_rows.len(), "row count");
    assert!(!cli_rows.is_empty(), "fixture must produce rows");

    for (c, b) in cli_rows.iter().zip(bridge_rows) {
        assert_close(&c["range"], &b["range"], "range");
        assert_close(&c["drop"], &b["drop_linear"], "drop_linear");
        assert_close(&c["drop_adj"], &b["drop_adj"], "drop_adj");
        assert_close(&c["wind_drift"], &b["wind_linear"], "wind_linear");
        assert_close(&c["wind_adj"], &b["wind_adj"], "wind_adj");
        assert_close(&c["velocity"], &b["velocity"], "velocity");
        assert_close(&c["energy"], &b["energy"], "energy");
        assert_close(&c["time"], &b["time"], "time");
    }
}

#[test]
fn wind_card_matches_cli() {
    let cli_out = cli(&[
        "wind-card",
        "--zero-distance", "100",
        "-v", "2600", "-b", "0.243", "-m", "175", "-d", "0.308",
        "--drag-model", "g7",
        "--sight-height", "1.5",
        "--start", "100", "--end", "600", "--step", "100",
        "--wind-speeds", "5,10,15",
        "--adjustment-unit", "mil",
        "-o", "json",
    ]);
    let mut request = fixture_request();
    request["wind_speeds"] = json!([5.0, 10.0, 15.0]);
    let bridge_out = bridge("card.wind", request);

    // CLI wind-card JSON: {"cards":[{"wind_direction":..,"data":[{"range":..,"drifts":[..]}]}]}
    // or a flat {"data": ...} shape for the single default angle — accept either.
    let cli_rows = cli_out["data"]
        .as_array()
        .or_else(|| cli_out["cards"][0]["data"].as_array())
        .expect("CLI wind-card rows");
    let bridge_rows = bridge_out["rows"].as_array().expect("bridge rows");
    assert_eq!(cli_rows.len(), bridge_rows.len(), "row count");
    assert!(!cli_rows.is_empty(), "fixture must produce rows");

    // The CLI emits one "wind_<speed>" key per swept speed; the bridge emits an
    // ordered wind_columns array aligned with the request's wind_speeds.
    let speeds = [5.0_f64, 10.0, 15.0];
    for (c, b) in cli_rows.iter().zip(bridge_rows) {
        assert_close(&c["range"], &b["range"], "range");
        let bridge_drifts = b["wind_columns"].as_array().expect("bridge wind_columns");
        assert_eq!(speeds.len(), bridge_drifts.len(), "column count");
        for (i, speed) in speeds.iter().enumerate() {
            let key = format!("wind_{speed}");
            let cli_cell = &c[&key];
            assert!(!cli_cell.is_null(), "CLI row missing {key}: {c}");
            assert_close(cli_cell, &bridge_drifts[i], &format!("drift col {key}"));
        }
    }
}
