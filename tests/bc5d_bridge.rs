//! BC5D offline correction over the mobile bridge (solve `corrections.bc5d_table_path`,
//! card `bc_segments`/`bc5d_table_path`, and `bc5d.info`).
//!
//! Fixture provenance: every table here is SYNTHESIZED in-test in the BC5D v2 binary
//! layout (same builder shape as `tests/bc_table_adjustment.rs`) — no real production
//! `.bin` is committed. The single-plane fixture keys its correction off the BC axis
//! (0.4 -> 0.8, 0.5 -> 1.2), so a 0.4-BC load gets a strong, easily asserted 0.8
//! correction across the whole velocity envelope.
//!
//! The golden test executes the actual CLI (`trajectory --bc-table-dir`) against the
//! bridge `solve` with `bc5d_table_path` on the same load and requires the sampled
//! drops to agree — the CLI remains the reference implementation for BC5D semantics.

#![cfg(all(feature = "bridge", feature = "cli"))]

use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[cfg(feature = "pdf")]
#[path = "support/pdf_text.rs"]
mod pdf_text;

// ---------------------------------------------------------------------------
// Synthesized BC5D v2 fixture (see module doc for provenance)
// ---------------------------------------------------------------------------

fn crc32_ieee(data: &[u8]) -> u32 {
    let mut crc = 0xffff_ffff_u32;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            let mask = 0_u32.wrapping_sub(crc & 1);
            crc = (crc >> 1) ^ (0xedb8_8320 & mask);
        }
    }
    !crc
}

fn push_f32s(output: &mut Vec<u8>, values: &[f32]) {
    for value in values {
        output.extend_from_slice(&value.to_le_bytes());
    }
}

/// One G1 plane, correction keyed off the BC axis: BC 0.4 -> 0.8, BC 0.5 -> 1.2.
fn bc5d_fixture() -> Vec<u8> {
    let weight_bins = [168.0_f32];
    let bc_bins = [0.4_f32, 0.5_f32];
    let muzzle_bins = [2_800.0_f32];
    let current_bins = [2_000.0_f32];
    let data = [0.8_f32, 1.2_f32];

    let mut checksum_data = Vec::new();
    push_f32s(&mut checksum_data, &weight_bins);
    push_f32s(&mut checksum_data, &bc_bins);
    push_f32s(&mut checksum_data, &muzzle_bins);
    push_f32s(&mut checksum_data, &current_bins);
    push_f32s(&mut checksum_data, &data);

    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"BC5D");
    bytes.extend_from_slice(&2_u32.to_le_bytes()); // format version
    bytes.extend_from_slice(&0.308_f32.to_le_bytes()); // caliber
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // flags
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // padding
    bytes.extend_from_slice(&(weight_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(bc_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(muzzle_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(current_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&1_u32.to_le_bytes()); // drag types
    bytes.extend_from_slice(&0_u64.to_le_bytes()); // timestamp
    bytes.extend_from_slice(&crc32_ieee(&checksum_data).to_le_bytes());
    let mut api_version = [0_u8; 16];
    api_version[..4].copy_from_slice(b"test");
    bytes.extend_from_slice(&api_version);
    bytes.extend_from_slice(&[0_u8; 12]);
    bytes.extend_from_slice(&checksum_data);
    assert_eq!(bytes.len(), 108);
    bytes
}

fn unique_temp_dir(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bc5d-bridge-{label}-{}-{nonce}",
        std::process::id()
    ))
}

/// Write the fixture as `bc5d_308.bin` in a fresh temp dir; returns (dir, file path).
fn write_fixture(label: &str, bytes: &[u8]) -> (PathBuf, PathBuf) {
    let dir = unique_temp_dir(label);
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("bc5d_308.bin");
    std::fs::write(&file, bytes).unwrap();
    (dir, file)
}

// ---------------------------------------------------------------------------
// Bridge helpers
// ---------------------------------------------------------------------------

fn bridge_raw(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    let raw = ballistics_engine::bridge::bridge_call(&envelope.to_string());
    serde_json::from_str(&raw).expect("bridge output parses")
}

fn bridge_ok(command: &str, request: Value) -> Value {
    let out = bridge_raw(command, request);
    assert_eq!(out["ok"], true, "bridge {command} failed: {out}");
    out["result"].clone()
}

/// The shared card load: 168 gr .308 at 2800 fps, G1 0.4 — the fixture's hot cell.
fn card_request(extra: Value) -> Value {
    let mut request = json!({
        "units": "imperial",
        "muzzle_velocity": 2800.0,
        "ballistic_coefficient": 0.4,
        "mass": 168.0,
        "diameter": 0.308,
        "drag_model": "g1",
        "sight_height": 1.5,
        "zero_distance": 100.0,
        "start": 100.0,
        "end": 600.0,
        "step": 100.0,
        "adjustment_unit": "mil"
    });
    for (key, value) in extra.as_object().expect("extra must be an object") {
        request[key] = value.clone();
    }
    request
}

/// The same load as `card_request`, as a solve-json v1 document (SI units).
fn solve_request(corrections: Option<Value>) -> Value {
    let mut request = json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 168.0 * 0.000_064_798_91,
            "diameter_m": 0.308 * 0.0254,
            "drag_model": "G1",
            "ballistic_coefficient": 0.4
        },
        "rifle": {
            "muzzle_velocity_mps": 2800.0 * 0.3048,
            "sight_height_m": 1.5 * 0.0254
        },
        "shot": {
            "max_range_m": 600.0 * 0.9144,
            "zero_distance_m": 100.0 * 0.9144,
            // CLI-parity zero convention: the CLI's auto-zero converges the bullet onto
            // the LINE OF SIGHT at the zero range, i.e. onto the world-vertical height
            // of the sight above the (0-height) muzzle. solve-json's target_height_m
            // defaults to 0 (return to bore height), so state the LOS height explicitly.
            "target_height_m": 1.5 * 0.0254
        },
        // The CLI's imperial defaults, stated explicitly in SI: 59 F and 29.92 inHg
        // (via the CLI's own 33.8639 hPa/inHg constant), 50 % RH, sea level.
        "atmosphere": {
            "altitude_m": 0.0,
            "temperature_k": 288.15,
            "pressure_pa": 29.92 * 33.8639 * 100.0,
            "relative_humidity": 0.5
        },
        "wind": {},
        "solver": {},
        "effects": {},
        "sampling": {"interval_m": 100.0}
    });
    if let Some(corrections) = corrections {
        request["corrections"] = corrections;
    }
    request
}

fn last_row_drop_adj(result: &Value) -> f64 {
    result["rows"]
        .as_array()
        .expect("rows")
        .last()
        .expect("at least one row")["drop_adj"]
        .as_f64()
        .expect("drop_adj")
}

// ---------------------------------------------------------------------------
// bc5d.info
// ---------------------------------------------------------------------------

#[test]
fn bc5d_info_reports_table_metadata() {
    let (dir, file) = write_fixture("info", &bc5d_fixture());
    let result = bridge_ok("bc5d.info", json!({"path": file.to_str().unwrap()}));
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(result["valid"], true);
    assert_eq!(result["crc_ok"], true);
    assert_eq!(result["format_version"], 2);
    assert!((result["caliber"].as_f64().unwrap() - 0.308).abs() < 1e-4);
    assert_eq!(result["api_version"], "test");
    assert_eq!(result["bins"]["weight"], 1);
    assert_eq!(result["bins"]["bc"], 2);
    assert_eq!(result["bins"]["muzzle_velocity"], 1);
    assert_eq!(result["bins"]["current_velocity"], 1);
    assert_eq!(result["bins"]["drag_types"], 1);
    assert_eq!(result["total_cells"], 2);
    assert_eq!(result["weight_range_grains"][0], 168.0);
    assert_eq!(result["velocity_range_fps"][0], 2000.0);
}

#[test]
fn bc5d_info_rejects_a_corrupt_crc_with_a_clean_envelope() {
    let mut bytes = bc5d_fixture();
    let last = bytes.len() - 1;
    bytes[last] ^= 0xFF; // corrupt the data section; header CRC no longer matches
    let (dir, file) = write_fixture("info-corrupt", &bytes);
    let out = bridge_raw("bc5d.info", json!({"path": file.to_str().unwrap()}));
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    let message = out["error"]["message"].as_str().unwrap();
    assert!(
        message.contains("Checksum mismatch"),
        "corruption must be named as a checksum failure: {message}"
    );
}

#[test]
fn bc5d_info_rejects_wrong_magic() {
    let mut bytes = bc5d_fixture();
    bytes[..4].copy_from_slice(b"NOPE");
    let (dir, file) = write_fixture("info-magic", &bytes);
    let out = bridge_raw("bc5d.info", json!({"path": file.to_str().unwrap()}));
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    assert!(
        out["error"]["message"].as_str().unwrap().contains("magic"),
        "{out}"
    );
}

// ---------------------------------------------------------------------------
// Cards: bc5d_table_path and bc_segments
// ---------------------------------------------------------------------------

#[test]
fn card_rows_change_when_a_bc5d_table_applies() {
    let (dir, file) = write_fixture("card-table", &bc5d_fixture());
    let baseline = bridge_ok("card.range_table", card_request(json!({})));
    let corrected = bridge_ok(
        "card.range_table",
        card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );

    let baseline_rows = baseline["rows"].as_array().unwrap();
    let corrected_rows = corrected["rows"].as_array().unwrap();
    assert_eq!(baseline_rows.len(), corrected_rows.len());
    assert!(!baseline_rows.is_empty());

    // The 0.8 correction degrades BC 0.4 -> 0.32: strictly more drop at distance.
    let baseline_drop = last_row_drop_adj(&baseline);
    let corrected_drop = last_row_drop_adj(&corrected);
    assert!(
        corrected_drop - baseline_drop > 0.3,
        "a 0.8 BC correction must add substantial dial-up at 600 yd: \
         baseline {baseline_drop} mil vs corrected {corrected_drop} mil"
    );

    // Same schedule on every card surface: come-ups and the wind card apply it too.
    let come_ups_corrected = bridge_ok(
        "card.come_ups",
        card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(
        (last_row_drop_adj(&come_ups_corrected) - corrected_drop).abs() < 1e-9,
        "come-ups and range-table must dial the same corrected elevation"
    );
}

#[test]
fn card_wind_card_applies_the_schedule() {
    let (dir, file) = write_fixture("card-wind", &bc5d_fixture());
    let extra = json!({"wind_speeds": [10.0]});
    let baseline = bridge_ok("card.wind", card_request(extra.clone()));
    let corrected = bridge_ok(
        "card.wind",
        card_request(json!({
            "wind_speeds": [10.0],
            "bc5d_table_path": file.to_str().unwrap()
        })),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    let baseline_drift = baseline["rows"].as_array().unwrap().last().unwrap()["wind_columns"][0]
        .as_f64()
        .unwrap();
    let corrected_drift = corrected["rows"].as_array().unwrap().last().unwrap()["wind_columns"][0]
        .as_f64()
        .unwrap();
    assert!(
        (corrected_drift - baseline_drift).abs() > 0.05,
        "a degraded BC must change 600 yd wind drift: {baseline_drift} vs {corrected_drift}"
    );
}

#[test]
fn explicit_bc_segments_win_over_the_table() {
    let (dir, file) = write_fixture("card-precedence", &bc5d_fixture());
    // A schedule deliberately DIFFERENT from anything the table would generate
    // (table: bc 0.32 everywhere; explicit: 0.45), covering the whole envelope.
    let segments = json!([{"velocity_min": 500.0, "velocity_max": 4000.0, "bc": 0.45}]);
    let with_both = bridge_ok(
        "card.range_table",
        card_request(json!({
            "bc_segments": segments.clone(),
            "bc5d_table_path": file.to_str().unwrap()
        })),
    );
    let segments_only = bridge_ok(
        "card.range_table",
        card_request(json!({"bc_segments": segments})),
    );
    let table_only = bridge_ok(
        "card.range_table",
        card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(
        serde_json::to_string(&with_both).unwrap(),
        serde_json::to_string(&segments_only).unwrap(),
        "explicit bc_segments must fully override the table"
    );
    assert!(
        (last_row_drop_adj(&with_both) - last_row_drop_adj(&table_only)).abs() > 0.3,
        "precedence test is vacuous if the two schedules produce the same card"
    );
}

#[test]
fn card_rejects_malformed_bc_segments_and_bad_tables() {
    // velocity_min >= velocity_max
    let out = bridge_raw(
        "card.range_table",
        card_request(json!({
            "bc_segments": [{"velocity_min": 2000.0, "velocity_max": 2000.0, "bc": 0.4}]
        })),
    );
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    assert!(
        out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("velocity_min must be < velocity_max"),
        "{out}"
    );

    // Non-positive BC
    let out = bridge_raw(
        "card.range_table",
        card_request(json!({
            "bc_segments": [{"velocity_min": 1000.0, "velocity_max": 2000.0, "bc": 0.0}]
        })),
    );
    assert!(
        out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("bc must be > 0"),
        "{out}"
    );

    // Corrupt table file: clean command_failed, not a panic or a silent no-op.
    let mut bytes = bc5d_fixture();
    let last = bytes.len() - 1;
    bytes[last] ^= 0xFF;
    let (dir, file) = write_fixture("card-corrupt", &bytes);
    let out = bridge_raw(
        "card.range_table",
        card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    assert!(
        out["error"]["message"]
            .as_str()
            .unwrap()
            .contains("Checksum mismatch"),
        "{out}"
    );
}

// ---------------------------------------------------------------------------
// Solve: corrections.bc5d_table_path
// ---------------------------------------------------------------------------

#[test]
fn solve_applies_the_table_to_zero_and_flight() {
    let (dir, file) = write_fixture("solve-table", &bc5d_fixture());
    let baseline = bridge_ok("solve", solve_request(None));
    let corrected = bridge_ok(
        "solve",
        solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()}))),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    // The corrections block is echoed on the resolved request (round-trip contract).
    assert_eq!(
        corrected["resolved_request"]["corrections"]["bc5d_table_path"],
        json!(file.to_str().unwrap())
    );
    assert!(baseline["resolved_request"].get("corrections").is_none());

    // Zero solved WITH the schedule: the degraded BC needs more launch angle.
    let baseline_angle = baseline["resolved_request"]["shot"]["muzzle_angle_rad"]
        .as_f64()
        .unwrap();
    let corrected_angle = corrected["resolved_request"]["shot"]["muzzle_angle_rad"]
        .as_f64()
        .unwrap();
    assert!(
        corrected_angle > baseline_angle,
        "zero must be solved with the corrected (worse) BC: {baseline_angle} vs {corrected_angle}"
    );

    // And the flight: strictly more drop at the terminal sample.
    let terminal_drop = |out: &Value| {
        out["samples"].as_array().unwrap().last().unwrap()["drop_m"]
            .as_f64()
            .unwrap()
    };
    let baseline_drop = terminal_drop(&baseline);
    let corrected_drop = terminal_drop(&corrected);
    assert!(
        corrected_drop - baseline_drop > 0.05,
        "a 0.8 BC correction must add drop at 600 yd: {baseline_drop} m vs {corrected_drop} m"
    );

    // G1 request against a G1 table: no coercion warning.
    assert!(
        !corrected["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w["code"] == "bc5d_drag_model_coerced"),
        "{corrected}"
    );
}

#[test]
fn solve_warns_when_an_aux_drag_model_is_coerced_to_g1() {
    let (dir, file) = write_fixture("solve-coerce", &bc5d_fixture());
    let mut request = solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()})));
    request["projectile"]["drag_model"] = json!("G5");
    let corrected = bridge_ok("solve", request);
    std::fs::remove_dir_all(&dir).unwrap();

    let warnings = corrected["warnings"].as_array().unwrap();
    assert!(
        warnings
            .iter()
            .any(|w| w["code"] == "bc5d_drag_model_coerced"
                && w["path"] == "$.corrections.bc5d_table_path"),
        "a G5 request against a G1/G7 table must carry the coercion warning: {warnings:?}"
    );
}

/// MBA-1597: a correction table adjusts a G1/G7 BC along the flight, and a caller-supplied
/// `projectile.drag_table` replaces the G-curve that BC belongs to, so the pair is refused —
/// even with a valid table on disk, which is never opened.
#[test]
fn solve_refuses_a_correction_table_alongside_a_drag_table() {
    let (dir, file) = write_fixture("solve-with-drag-table", &bc5d_fixture());
    let mut request = solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()})));
    request["projectile"]["drag_table"] = json!({
        "kind": "reference",
        "points": [{"mach": 0.0, "cd": 0.25}, {"mach": 3.0, "cd": 0.25}]
    });
    let out = bridge_raw("solve", request);
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["code"], "conflicting_fields",
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"], "$.corrections.bc5d_table_path",
        "{out}"
    );
}

/// MBA-1598: a correction table builds its own BC ladder from `ballistic_coefficient`, so a
/// caller's `projectile.bc_segments` beside it is refused. The refusal comes before the table
/// is opened: a corrupt table on disk would otherwise answer with its own CRC error at
/// `$.corrections.bc5d_table_path`, and a valid one gets the same refusal.
#[test]
fn solve_refuses_a_correction_table_alongside_a_bc_ladder() {
    let mut corrupt = bc5d_fixture();
    let last = corrupt.len() - 1;
    corrupt[last] ^= 0xFF;
    for (label, bytes) in [
        ("solve-with-ladder", bc5d_fixture()),
        ("solve-with-ladder-corrupt", corrupt),
    ] {
        let (dir, file) = write_fixture(label, &bytes);
        let mut request = solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()})));
        request["projectile"]["bc_segments"] = json!([
            {"velocity_min_mps": 600.0, "velocity_max_mps": 999.0, "ballistic_coefficient": 0.42},
            {"velocity_min_mps": 0.0, "velocity_max_mps": 600.0, "ballistic_coefficient": 0.38}
        ]);
        let out = bridge_raw("solve", request);
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(out["ok"], false, "{label}: {out}");
        assert_eq!(
            out["error"]["details"]["error"]["code"], "conflicting_fields",
            "{label}: {out}"
        );
        assert_eq!(
            out["error"]["details"]["error"]["path"], "$.projectile.bc_segments",
            "{label}: {out}"
        );
        assert_eq!(
            out["error"]["details"]["error"]["message"],
            "projectile.bc_segments cannot be combined with corrections.bc5d_table_path: the \
             correction table builds its own BC ladder from ballistic_coefficient",
            "{label}: {out}"
        );
    }
}

/// MBA-1598: with a ladder, a drag table and a correction table all present, the older drag
/// table check answers first, at the correction table's path, as docs/SOLVE_JSON_V1.md says.
#[test]
fn solve_refuses_a_ladder_drag_table_and_correction_table_at_the_older_path() {
    let (dir, file) = write_fixture("solve-with-all-three", &bc5d_fixture());
    let mut request = solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()})));
    request["projectile"]["drag_table"] = json!({
        "kind": "reference",
        "points": [{"mach": 0.0, "cd": 0.25}, {"mach": 3.0, "cd": 0.25}]
    });
    request["projectile"]["bc_segments"] = json!([
        {"velocity_min_mps": 0.0, "velocity_max_mps": 999.0, "ballistic_coefficient": 0.42}
    ]);
    let out = bridge_raw("solve", request);
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["code"], "conflicting_fields",
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"], "$.corrections.bc5d_table_path",
        "{out}"
    );
}

#[test]
fn solve_rejects_corrupt_and_missing_tables_with_typed_envelopes() {
    // Corrupt CRC -> invalid_value at the corrections path.
    let mut bytes = bc5d_fixture();
    let last = bytes.len() - 1;
    bytes[last] ^= 0xFF;
    let (dir, file) = write_fixture("solve-corrupt", &bytes);
    let out = bridge_raw(
        "solve",
        solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()}))),
    );
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(out["ok"], false, "{out}");
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["code"], "invalid_value",
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"], "$.corrections.bc5d_table_path",
        "{out}"
    );

    // Missing file -> io_error at the same path.
    let out = bridge_raw(
        "solve",
        solve_request(Some(
            json!({"bc5d_table_path": "/nonexistent/bc5d_308.bin"}),
        )),
    );
    assert_eq!(
        out["error"]["details"]["error"]["code"], "io_error",
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"], "$.corrections.bc5d_table_path",
        "{out}"
    );

    // Unknown members inside the corrections block are a schema rejection.
    let out = bridge_raw(
        "solve",
        solve_request(Some(json!({"bc5d_table_path": "x", "extra": 1}))),
    );
    assert_eq!(
        out["error"]["details"]["error"]["code"], "unknown_field",
        "{out}"
    );

    // A non-string path is a shape error, not a coerced value.
    let out = bridge_raw("solve", solve_request(Some(json!({"bc5d_table_path": 42}))));
    assert_eq!(
        out["error"]["details"]["error"]["code"], "invalid_value",
        "{out}"
    );
}

#[test]
fn solve_resolved_request_roundtrip_reapplies_not_compounds() {
    let (dir, file) = write_fixture("solve-roundtrip", &bc5d_fixture());
    let request = ballistics_engine::solve_json::decode_solve_request_v1(
        &solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()}))).to_string(),
    )
    .expect("decode");
    let first = ballistics_engine::solve_v1(request).expect("first solve");
    let rebuilt: ballistics_engine::solve_json::SolveRequestV1 = (&first.resolved_request).into();
    assert!(
        rebuilt.corrections.is_some(),
        "the corrections block must survive the resolved->request round-trip"
    );
    let second = ballistics_engine::solve_v1(rebuilt).expect("second solve");
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(
        serde_json::to_value(&first.resolved_request).unwrap(),
        serde_json::to_value(&second.resolved_request).unwrap(),
        "resolved request changed after a round-trip"
    );
    assert_eq!(
        first.samples, second.samples,
        "a round-tripped correction must re-apply identically, not compound"
    );
}

// ---------------------------------------------------------------------------
// Golden CLI parity: trajectory --bc-table-dir vs bridge solve bc5d_table_path
// ---------------------------------------------------------------------------

/// Run the reference CLI with sampled output and parse (distance_yd, drop_in) rows.
fn cli_sampled_drops(table_dir: Option<&Path>) -> Vec<(f64, f64)> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_ballistics"));
    command.args([
        "trajectory",
        "-v",
        "2800",
        "-b",
        "0.4",
        "-m",
        "168",
        "-d",
        "0.308",
        "--drag-model",
        "g1",
        "--sight-height",
        "1.5",
        "--auto-zero",
        "100",
        "--max-range",
        "600",
        "--sample-trajectory",
        "--sample-interval",
        "100", // meters, matching the solve request's interval_m
        "-o",
        "csv",
        "--full",
    ]);
    if let Some(dir) = table_dir {
        command.arg("--bc-table-dir").arg(dir);
    }
    let output = command.output().expect("run ballistics CLI");
    assert!(
        output.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("utf8");
    let mut rows = Vec::new();
    for line in stdout.lines().skip(1) {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 6 {
            continue;
        }
        let distance_yd: f64 = fields[0].parse().expect("distance field");
        let drop_in: f64 = fields[1].parse().expect("drop field");
        rows.push((distance_yd, drop_in));
    }
    assert!(
        rows.len() >= 4,
        "expected several sampled CLI rows:\n{stdout}"
    );
    rows
}

/// Bridge solve for the same load; returns (distance_m, drop_m) samples.
fn bridge_sampled_drops(table_path: Option<&Path>) -> Vec<(f64, f64)> {
    let corrections = table_path.map(|p| json!({"bc5d_table_path": p.to_str().unwrap()}));
    let result = bridge_ok("solve", solve_request(corrections));
    result["samples"]
        .as_array()
        .expect("samples")
        .iter()
        .map(|s| {
            (
                s["distance_m"].as_f64().unwrap(),
                s["drop_m"].as_f64().unwrap(),
            )
        })
        .collect()
}

#[test]
fn golden_cli_trajectory_and_bridge_solve_agree_on_bc5d_drops() {
    let (dir, _file) = write_fixture("golden", &bc5d_fixture());
    let table_file = dir.join("bc5d_308.bin");

    let cli_base = cli_sampled_drops(None);
    let cli_corrected = cli_sampled_drops(Some(&dir));
    let bridge_base = bridge_sampled_drops(None);
    let bridge_corrected = bridge_sampled_drops(Some(&table_file));
    std::fs::remove_dir_all(&dir).unwrap();

    // Match bridge samples to CLI rows by downrange distance (CLI prints yards; the
    // sample grid itself is metric on both sides). Returns the aligned pairs so the
    // delta comparison below reuses the exact same rows.
    let aligned =
        |cli: &[(f64, f64)], bridge: &[(f64, f64)], label: &str| -> Vec<(f64, f64, f64)> {
            let mut pairs = Vec::new();
            for &(distance_yd, cli_drop_in) in cli {
                let distance_m = distance_yd * 0.9144;
                let Some(&(_, bridge_drop_m)) =
                    bridge.iter().find(|(d, _)| (d - distance_m).abs() < 0.5)
                else {
                    continue;
                };
                pairs.push((distance_m, cli_drop_in, bridge_drop_m * 39.3701));
            }
            assert!(
                pairs.len() >= 4,
                "{label}: too few aligned samples to be meaningful ({})",
                pairs.len()
            );
            pairs
        };

    let base_pairs = aligned(&cli_base, &bridge_base, "no table");
    let corrected_pairs = aligned(&cli_corrected, &bridge_corrected, "with BC5D table");

    // Direct agreement, every aligned range: the CLI prints drops at 0.01 in
    // resolution, so 0.05 in of headroom is rounding plus nothing.
    for (label, pairs) in [
        ("no table", &base_pairs),
        ("with BC5D table", &corrected_pairs),
    ] {
        for &(distance_m, cli_in, bridge_in) in pairs.iter() {
            assert!(
                (bridge_in - cli_in).abs() < 0.05,
                "{label}: drop mismatch at {distance_m:.0} m: CLI {cli_in} in vs bridge {bridge_in} in"
            );
        }
    }

    // The table must have MOVED both implementations, and by the same amount, at the
    // farthest range both sampled in both runs — the actual BC5D-parity assertion
    // (any shared systematics cancel in the delta).
    let (far_m, cli_base_in, bridge_base_in) = *base_pairs.last().unwrap();
    let &(corr_far_m, cli_corr_in, bridge_corr_in) = corrected_pairs
        .iter()
        .find(|(d, _, _)| (d - far_m).abs() < 0.5)
        .expect("corrected run must sample the same far range");
    let cli_delta = cli_corr_in - cli_base_in;
    let bridge_delta = bridge_corr_in - bridge_base_in;
    assert!(
        cli_delta > 2.0,
        "the 0.8 correction must visibly add drop at {corr_far_m:.0} m (CLI delta {cli_delta} in)"
    );
    assert!(
        (cli_delta - bridge_delta).abs() < 0.05,
        "CLI and bridge must attribute the same drop change to the table at {far_m:.0} m: \
         CLI {cli_delta} in vs bridge {bridge_delta} in"
    );
}

// ---------------------------------------------------------------------------
// MBA-1603: `explain` carries the correction table with ProjectileDrag
// ---------------------------------------------------------------------------

fn resolved_request(request: Value) -> ballistics_engine::solve_json::ResolvedSolveRequestV1 {
    let decoded = ballistics_engine::solve_json::decode_solve_request_v1(&request.to_string())
        .expect("request decodes");
    ballistics_engine::solve_v1(decoded)
        .expect("request solves")
        .resolved_request
}

/// Before MBA-1603 no group swapped `corrections`, so a table's whole effect sat in the
/// interaction remainder with `skipped_axes` empty: an unexplained difference, unexplained.
#[test]
fn explain_attributes_a_correction_table_to_projectile_drag() {
    use ballistics_engine::explain::explain_difference;
    use ballistics_engine::perturbation::InputGroup;

    let (dir, file) = write_fixture("explain-attribution", &bc5d_fixture());
    let plain = resolved_request(solve_request(None));
    let corrected = resolved_request(solve_request(Some(
        json!({"bc5d_table_path": file.to_str().unwrap()}),
    )));
    let report = explain_difference(&plain, &corrected, &[300.0, 500.0]);
    std::fs::remove_dir_all(&dir).unwrap();
    let report = report.expect("explain runs");

    assert!(report.skipped_axes.is_empty(), "{:?}", report.skipped_axes);
    for row in &report.rows {
        let drag = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::ProjectileDrag)
            .unwrap()
            .delta;
        assert!(
            row.total.drop_m > 1e-3,
            "the 0.8 correction must add drop at {} m: {:?}",
            row.range_m,
            row.total
        );
        assert!(
            (drag.drop_m - row.total.drop_m).abs() < 1e-9,
            "range {}: ProjectileDrag {} vs total {} (remainder {})",
            row.range_m,
            drag.drop_m,
            row.total.drop_m,
            row.interaction_remainder.drop_m
        );
        assert!(row.interaction_remainder.drop_m.abs() < 1e-9);
    }
}

/// The swap used to write the OTHER request's diameter onto a request still carrying its own
/// .308 table, which the caliber check refuses, so the whole comparison aborted with
/// `invalid_value` at `bc5d_table_path`. The table now leaves before the diameter arrives.
#[test]
fn explain_compares_a_corrected_308_with_a_224_instead_of_aborting() {
    use ballistics_engine::explain::explain_difference;

    let (dir, file) = write_fixture("explain-caliber", &bc5d_fixture());
    let corrected_308 = resolved_request(solve_request(Some(
        json!({"bc5d_table_path": file.to_str().unwrap()}),
    )));
    let mut plain_224 = solve_request(None);
    plain_224["projectile"]["diameter_m"] = json!(0.224 * 0.0254);
    let plain_224 = resolved_request(plain_224);
    let report = explain_difference(&corrected_308, &plain_224, &[300.0]);
    std::fs::remove_dir_all(&dir).unwrap();

    let report = report.expect("a foreign diameter must never meet the table mid-swap");
    assert!(report.rows[0].interaction_remainder.drop_m.abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// MBA-1600: card.* refuses a BC ladder it could not fly as written
// ---------------------------------------------------------------------------

/// The refusal `command` gives for `extra` on the shared card load, checked to be a typed
/// `command_failed` envelope rather than a card.
fn card_refusal(command: &str, extra: Value) -> String {
    let out = bridge_raw(command, card_request(extra));
    assert_eq!(out["ok"], false, "{command} must refuse this ladder: {out}");
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    out["error"]["message"].as_str().unwrap().to_string()
}

fn band(velocity_min: f64, velocity_max: f64, bc: f64) -> Value {
    json!({"velocity_min": velocity_min, "velocity_max": velocity_max, "bc": bc})
}

/// The engine's lookup flies the FIRST listed band containing a speed and builds no blend
/// across an overlap, so an overlapping ladder used to fly differently by listing order alone
/// (0.9 mil apart at 1000 yd for this load), and a nested band was silently dead.
#[test]
fn card_refuses_overlapping_bc_segments_in_any_order() {
    const RULE: &str = "bands may touch or leave gaps but must not overlap";
    let cases = [
        (
            json!({"bc_segments": [band(2000.0, 3500.0, 0.5), band(1000.0, 2200.0, 0.3)]}),
            "bc_segments[0] (2000 to 3500 fps) overlaps bc_segments[1] (1000 to 2200 fps)",
        ),
        (
            json!({"bc_segments": [band(1000.0, 2200.0, 0.3), band(2000.0, 3500.0, 0.5)]}),
            "bc_segments[1] (2000 to 3500 fps) overlaps bc_segments[0] (1000 to 2200 fps)",
        ),
        (
            json!({"bc_segments": [band(1000.0, 3500.0, 0.5), band(1800.0, 2200.0, 0.3)]}),
            "bc_segments[1] (1800 to 2200 fps) overlaps bc_segments[0] (1000 to 3500 fps)",
        ),
        (
            json!({"bc_segments": [band(1000.0, 2000.0, 0.4), band(1000.0, 1500.0, 0.3)]}),
            "bc_segments[1] (1000 to 1500 fps) overlaps bc_segments[0] (1000 to 2000 fps)",
        ),
        (
            json!({
                "velocity_unit": "mps",
                "muzzle_velocity": 853.0,
                "bc_segments": [band(600.0, 1100.0, 0.5), band(300.0, 650.0, 0.3)]
            }),
            "bc_segments[0] (600 to 1100 m/s) overlaps bc_segments[1] (300 to 650 m/s)",
        ),
    ];
    for (extra, expected) in cases {
        assert_eq!(
            card_refusal("card.range_table", extra),
            format!("card.range_table failed: invalid request: {expected}: {RULE}")
        );
    }

    // Every surface that flies the ladder refuses it the same way.
    let overlap = [band(2000.0, 3500.0, 0.5), band(1000.0, 2200.0, 0.3)];
    let mut surfaces = vec![
        ("card.come_ups", json!({"bc_segments": overlap})),
        (
            "card.wind",
            json!({"wind_speeds": [10.0], "bc_segments": overlap}),
        ),
    ];
    if cfg!(feature = "pdf") {
        surfaces.push(("card.pdf", json!({"bc_segments": overlap})));
    }
    for (command, extra) in surfaces {
        let message = card_refusal(command, extra);
        assert!(
            message.starts_with(&format!(
                "{command} failed: invalid request: bc_segments[0] (2000 to 3500 fps) overlaps"
            )),
            "{message}"
        );
    }
}

/// The tightening stops at overlaps: card's own conventions (any order, touching edges, gaps,
/// `[]` and `null` as no ladder) are unchanged, and so are the cards they give.
#[test]
fn card_still_accepts_touching_gapped_and_reordered_bc_segments() {
    let touching = bridge_ok(
        "card.range_table",
        card_request(json!({
            "bc_segments": [band(2000.0, 3500.0, 0.5), band(1000.0, 2000.0, 0.3)]
        })),
    );
    let touching_reversed = bridge_ok(
        "card.range_table",
        card_request(json!({
            "bc_segments": [band(1000.0, 2000.0, 0.3), band(2000.0, 3500.0, 0.5)]
        })),
    );
    // The rows, not the whole response: `bc_ladder` echoes the bands in request order
    // (MBA-1601), so the two responses differ there and only there.
    assert_eq!(
        serde_json::to_string(&touching["rows"]).unwrap(),
        serde_json::to_string(&touching_reversed["rows"]).unwrap()
    );

    bridge_ok(
        "card.range_table",
        card_request(json!({
            "bc_segments": [
                band(2400.0, 3500.0, 0.5),
                band(1500.0, 2000.0, 0.3),
                band(0.0, 1200.0, 0.25)
            ]
        })),
    );

    let plain =
        serde_json::to_string(&bridge_ok("card.range_table", card_request(json!({})))).unwrap();
    for empty in [json!([]), Value::Null] {
        let out = bridge_ok(
            "card.range_table",
            card_request(json!({ "bc_segments": empty })),
        );
        assert_eq!(serde_json::to_string(&out).unwrap(), plain);
    }
}

/// A band BC is held to the bound solve-json puts on any BC. 1e6 used to fly a vacuum card,
/// and 1e-9 failed as an unexplained `zero solve failed`.
#[test]
fn card_bounds_each_band_bc_like_solve_json() {
    let cases = [
        (
            json!([band(0.0, 5000.0, 1e6)]),
            "bc_segments[0]: bc must be between 0.0001 and 100, got 1000000",
        ),
        (
            json!([band(0.0, 5000.0, 1e-9)]),
            "bc_segments[0]: bc must be between 0.0001 and 100, got 0.000000001",
        ),
        (
            json!([band(2000.0, 5000.0, 0.4), band(0.0, 500.0, 9.99e-5)]),
            "bc_segments[1]: bc must be between 0.0001 and 100, got 0.0000999",
        ),
    ];
    for (ladder, expected) in cases {
        assert_eq!(
            card_refusal("card.range_table", json!({ "bc_segments": ladder })),
            format!("card.range_table failed: invalid request: {expected}")
        );
    }

    // Both edges are inclusive, as in solve-json.
    for ladder in [
        json!([band(0.0, 5000.0, 100.0)]),
        json!([band(2000.0, 5000.0, 0.4), band(0.0, 500.0, 1e-4)]),
    ] {
        bridge_ok(
            "card.range_table",
            card_request(json!({ "bc_segments": ladder })),
        );
    }
}

/// The lookup scans the ladder on every derivative evaluation, so solve time grew with the
/// band count and nothing capped it (9000 bands took about 100 times as long as one).
#[test]
fn card_caps_bc_segments_at_the_solve_json_band_limit() {
    use ballistics_engine::solve_json::MAX_BC_SEGMENTS_V1;

    // Contiguous, fastest first, 0 to 4000 fps, every band at the scalar's 0.4.
    let ladder = |n: usize| -> Value {
        (0..n)
            .map(|i| {
                band(
                    4000.0 * (n - i - 1) as f64 / n as f64,
                    4000.0 * (n - i) as f64 / n as f64,
                    0.4,
                )
            })
            .collect()
    };

    bridge_ok(
        "card.range_table",
        card_request(json!({ "bc_segments": ladder(MAX_BC_SEGMENTS_V1) })),
    );
    assert_eq!(
        card_refusal(
            "card.range_table",
            json!({ "bc_segments": ladder(MAX_BC_SEGMENTS_V1 + 1) })
        ),
        format!(
            "card.range_table failed: invalid request: bc_segments must have at most {} bands, \
             got {}",
            MAX_BC_SEGMENTS_V1,
            MAX_BC_SEGMENTS_V1 + 1
        )
    );
}

// ---------------------------------------------------------------------------
// MBA-1601: a card states the BC ladder that flew, not only the scalar
// ---------------------------------------------------------------------------

/// The velocity-varying correction: G1 plane, BC 0.4 corrected 0.8 at 1000 fps and 1.0 at
/// 2800 fps (BC 0.5: 1.0 and 1.2). A 0.4 load at 2800 fps keeps its scalar at the muzzle,
/// while every band below it flies a smaller BC.
fn bc5d_velocity_fixture() -> Vec<u8> {
    let weight_bins = [168.0_f32];
    let bc_bins = [0.4_f32, 0.5_f32];
    let muzzle_bins = [2_800.0_f32];
    let current_bins = [1_000.0_f32, 2_800.0_f32];
    let data = [0.8_f32, 1.0_f32, 1.0_f32, 1.2_f32];

    let mut checksum_data = Vec::new();
    push_f32s(&mut checksum_data, &weight_bins);
    push_f32s(&mut checksum_data, &bc_bins);
    push_f32s(&mut checksum_data, &muzzle_bins);
    push_f32s(&mut checksum_data, &current_bins);
    push_f32s(&mut checksum_data, &data);

    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"BC5D");
    bytes.extend_from_slice(&2_u32.to_le_bytes()); // format version
    bytes.extend_from_slice(&0.308_f32.to_le_bytes()); // caliber
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // flags
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // padding
    bytes.extend_from_slice(&(weight_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(bc_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(muzzle_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&(current_bins.len() as u32).to_le_bytes());
    bytes.extend_from_slice(&1_u32.to_le_bytes()); // drag types
    bytes.extend_from_slice(&0_u64.to_le_bytes()); // timestamp
    bytes.extend_from_slice(&crc32_ieee(&checksum_data).to_le_bytes());
    let mut api_version = [0_u8; 16];
    api_version[..4].copy_from_slice(b"test");
    bytes.extend_from_slice(&api_version);
    bytes.extend_from_slice(&[0_u8; 12]);
    bytes.extend_from_slice(&checksum_data);
    // 80 header bytes plus 10 f32 cells: one current bin and two data cells more than
    // `bc5d_fixture`.
    assert_eq!(bytes.len(), 120);
    bytes
}

/// Under a gapless ladder the scalar never flies (the rows are identical at 0.4 and 0.2), yet
/// `bc_for_solve` was the only BC a card response stated.
#[test]
fn card_responses_carry_the_ladder_that_flew() {
    let ladder = json!([band(2000.0, 4000.0, 0.45), band(500.0, 2000.0, 0.42)]);
    for (command, extra) in [
        ("card.range_table", json!({})),
        ("card.come_ups", json!({})),
        ("card.wind", json!({"wind_speeds": [10.0]})),
    ] {
        let mut request = card_request(extra);
        request["bc_segments"] = ladder.clone();
        let out = bridge_ok(command, request.clone());
        assert_eq!(out["bc_for_solve"], json!(0.4), "{command}: {out}");
        assert_eq!(
            out["bc_ladder"]["source"], "bc_segments",
            "{command}: {out}"
        );
        assert_eq!(out["bc_ladder"]["segments"], ladder, "{command}: {out}");
        assert_eq!(out["bc_ladder"]["interior_gap"], false, "{command}: {out}");

        request["ballistic_coefficient"] = json!(0.2);
        let other_scalar = bridge_ok(command, request);
        assert_eq!(
            out["rows"], other_scalar["rows"],
            "{command}: the scalar must not fly under a gapless ladder"
        );
    }

    // No ladder, no block: those responses are byte-identical to before.
    for (command, extra) in [
        ("card.range_table", json!({})),
        ("card.come_ups", json!({})),
        ("card.wind", json!({"wind_speeds": [10.0]})),
    ] {
        let out = bridge_ok(command, card_request(extra));
        assert!(out.get("bc_ladder").is_none(), "{command}: {out}");
    }
}

#[test]
fn an_interior_gap_is_reported_because_the_scalar_flies_there() {
    let ladder = json!([band(2300.0, 4000.0, 0.45), band(500.0, 1800.0, 0.42)]);
    let at = |bc: f64| {
        bridge_ok(
            "card.range_table",
            card_request(json!({"bc_segments": ladder, "ballistic_coefficient": bc})),
        )
    };
    let (published, halved) = (at(0.4), at(0.2));
    assert_eq!(published["bc_ladder"]["interior_gap"], true, "{published}");
    assert_ne!(
        published["rows"], halved["rows"],
        "the scalar flies in the 1800-2300 fps gap, so it must move the rows"
    );
}

#[test]
fn a_correction_table_card_reports_its_generated_ladder() {
    let (dir, file) = write_fixture("card-ladder", &bc5d_velocity_fixture());
    let table = file.to_str().unwrap();
    let out = bridge_ok(
        "card.range_table",
        card_request(json!({ "bc5d_table_path": table })),
    );
    let metric = bridge_ok(
        "card.range_table",
        card_request(json!({
            "units": "metric",
            "muzzle_velocity": 853.44,
            "mass": 10.886,
            "diameter": 7.8232,
            "sight_height": 38.1,
            "bc5d_table_path": table
        })),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    let ladder = &out["bc_ladder"];
    assert_eq!(ladder["source"], "bc5d_table", "{out}");
    assert_eq!(ladder["interior_gap"], false, "{out}");
    let segments = ladder["segments"].as_array().unwrap();
    assert_eq!(segments.len(), 26, "{out}");
    assert_eq!(segments[0]["velocity_max"], json!(2800.0), "{out}");
    assert_eq!(segments[25]["velocity_min"], json!(500.0), "{out}");
    // The muzzle correction is 1.0, so the scalar is the published 0.4 and every band flies
    // less: the scalar was never the card's BC.
    let scalar = out["bc_for_solve"].as_f64().unwrap();
    assert_eq!(scalar, 0.4, "{out}");
    for segment in segments {
        assert!(segment["bc"].as_f64().unwrap() < scalar - 1e-3, "{segment}");
    }

    // The same ladder in the request's own velocity unit: the top band ends at the metric
    // muzzle velocity, round-tripped.
    let segments = metric["bc_ladder"]["segments"].as_array().unwrap();
    assert_eq!(segments.len(), 26, "{metric}");
    assert!(
        (segments[0]["velocity_max"].as_f64().unwrap() - 853.44).abs() < 1e-9,
        "{metric}"
    );
}

/// The printed footer of a correction-table card states the ladder the table generated, not the
/// scalar (MBA-1601). The apps send `bc5d_table_path`, so this is the footer on their printed
/// cards. Under the velocity fixture the ladder runs 0.320 to 0.398 and the 0.400 scalar never
/// flies, yet the footer used to read `BC:0.400 (g1)`.
#[cfg(feature = "pdf")]
#[test]
fn a_correction_table_cards_pdf_footer_states_the_generated_ladder() {
    let (dir, file) = write_fixture("card-ladder-pdf", &bc5d_velocity_fixture());
    let request: ballistics_engine::card_service::CardRequestV1 = serde_json::from_value(
        card_request(json!({ "bc5d_table_path": file.to_str().unwrap() })),
    )
    .unwrap();
    let card = ballistics_engine::card_service::pdf_card_v1(&request, None);
    std::fs::remove_dir_all(&dir).unwrap();

    let card = card.expect("a correction-table card must print");
    let printed = pdf_text::tokens(&card.pdf_bytes, "correction-table card");
    assert!(
        printed.iter().any(|t| t == "BC:0.320-0.398"),
        "the footer must state the generated ladder's range: {printed:?}"
    );
    assert!(printed.iter().any(|t| t == "ladder)"), "{printed:?}");
    assert!(
        !printed.iter().any(|t| t == "BC:0.400"),
        "the scalar never flew on this card: {printed:?}"
    );
}
