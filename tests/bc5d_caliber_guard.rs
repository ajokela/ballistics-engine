//! BC5D caliber identity: a table may only correct the caliber its own header declares.
//!
//! Nothing used to bind a table's CONTENT to the shot it was applied to. The file is
//! chosen by name (`bc5d_<key>.bin`) or handed over as a caller-supplied path, the CRC
//! only proves the bytes are intact, and every lookup clamps to the edge bins — so a
//! `.224` table applied to a `.308` shot returned a full, plausible-looking ladder that
//! silently biased every row (measured: segment BCs 0.4710..0.5114 where the correct
//! table gives 0.4989..0.5072 for a 175 gr G1 0.505 at 2600 fps; a `.308` table on a
//! `.243` shot came out -17.9 % on effective BC). A wrong-caliber table is worse than no
//! table, so it is now refused everywhere a shot is in hand.
//!
//! Fixture provenance: every table here is SYNTHESIZED in-test in the BC5D v2 binary
//! layout (the same builder shape as `tests/bc5d_bridge.rs`) — no real production `.bin`
//! is committed. Correction is keyed off the BC axis (0.4 -> 0.8, 0.5 -> 1.2), so the
//! 0.4-BC load below gets a strong 0.8 correction across the whole velocity envelope.

#![cfg(all(feature = "bridge", feature = "cli"))]

use serde_json::{json, Value};
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

// ---------------------------------------------------------------------------
// Synthesized BC5D v2 fixture, parameterized by header caliber
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
/// `caliber` is what the HEADER declares — the whole point of these tests.
fn bc5d_fixture(caliber: f32) -> Vec<u8> {
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
    bytes.extend_from_slice(&caliber.to_le_bytes());
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
    std::env::temp_dir().join(format!("bc5d-guard-{label}-{}-{nonce}", std::process::id()))
}

/// Write `bytes` as `file_name` in a fresh temp dir; returns (dir, file path).
fn write_fixture(label: &str, file_name: &str, bytes: &[u8]) -> (PathBuf, PathBuf) {
    let dir = unique_temp_dir(label);
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join(file_name);
    std::fs::write(&file, bytes).unwrap();
    (dir, file)
}

// ---------------------------------------------------------------------------
// Bridge helpers (same envelope shape as tests/bc5d_bridge.rs)
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

/// The shared load: 168 gr .308 at 2800 fps, G1 0.4 — the fixture's hot cell.
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

/// The same load in metric: 10.886 g, 7.8232 mm (= 0.308 in), 853.44 m/s.
fn metric_card_request(extra: Value) -> Value {
    let mut request = json!({
        "units": "metric",
        "muzzle_velocity": 2800.0 * 0.3048,
        "ballistic_coefficient": 0.4,
        "mass": 10.886,
        "diameter": 7.8232,
        "drag_model": "g1",
        "sight_height": 38.1,
        "zero_distance": 91.44,
        "start": 91.44,
        "end": 548.64,
        "step": 91.44,
        "adjustment_unit": "mil"
    });
    for (key, value) in extra.as_object().expect("extra must be an object") {
        request[key] = value.clone();
    }
    request
}

/// The imperial load as a solve-json v1 document (SI units), as in tests/bc5d_bridge.rs.
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
            "target_height_m": 1.5 * 0.0254
        },
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

fn drop_adjustments(result: &Value) -> Vec<f64> {
    result["rows"]
        .as_array()
        .expect("rows")
        .iter()
        .map(|row| row["drop_adj"].as_f64().expect("drop_adj"))
        .collect()
}

// ---------------------------------------------------------------------------
// A MATCHING table is untouched: the correction values are bit-identical
// ---------------------------------------------------------------------------

/// Golden ladder for the fixture at 168 gr / BC 0.4 / 2800 fps, captured from the engine
/// BEFORE the caliber guard existed (`origin/main`, commit d78a64b) by running the four
/// "matching table" tests in this file against that tree — they pass there unchanged,
/// which is the before/after equality claim. The guard is a gate, not a transform: for a
/// matching table every velocity band and every corrected BC must come back with the
/// identical BITS.
const GOLDEN_LADDER: [(f64, f64); 26] = [
    (2700.0, 2800.0),
    (2500.0, 2700.0),
    (2300.0, 2500.0),
    (2100.0, 2300.0),
    (2000.0, 2100.0),
    (1900.0, 2000.0),
    (1800.0, 1900.0),
    (1700.0, 1800.0),
    (1600.0, 1700.0),
    (1500.0, 1600.0),
    (1400.0, 1500.0),
    (1350.0, 1400.0),
    (1300.0, 1350.0),
    (1250.0, 1300.0),
    (1200.0, 1250.0),
    (1150.0, 1200.0),
    (1100.0, 1150.0),
    (1050.0, 1100.0),
    (1000.0, 1050.0),
    (950.0, 1000.0),
    (900.0, 950.0),
    (850.0, 900.0),
    (800.0, 850.0),
    (700.0, 800.0),
    (600.0, 700.0),
    (500.0, 600.0),
];

/// The corrected BC every band of that ladder carries: `0.4 * (0.8_f32 as f64)`, the
/// fixture's uniform correction, to the bit.
const GOLDEN_SEGMENT_BC: f64 = 0.320_000_004_768_371_6;

/// Card elevation (mil) at 100..600 yd with the matching table applied, captured from the
/// same pre-guard commit. Compared with `==` on purpose: these are the numbers the guard
/// must not move.
const GOLDEN_CARD_DROP_ADJ: [f64; 6] = [
    2.593_140_741_718_406_5e-5,
    0.558_401_068_207_466_7,
    1.384_719_268_983_751_5,
    2.407_150_226_123_941,
    3.645_404_844_306_932_5,
    5.143_799_279_278_284_5,
];

/// Terminal sampled drop (m) for the solve path with the matching table, same provenance.
const GOLDEN_SOLVE_TERMINAL_DROP_M: f64 = 2.822_094_036_583_238;

#[test]
fn matching_table_ladder_is_bit_identical() {
    let table = ballistics_engine::bc_table_5d::Bc5dTable::from_bytes(&bc5d_fixture(0.308))
        .expect("fixture parses");
    let segments = table
        .generate_segments(0.4, "G1", 168.0, Some(2800.0))
        .expect("the fixture carries a real correction");

    assert_eq!(
        segments.len(),
        GOLDEN_LADDER.len(),
        "the ladder gained or lost bands"
    );
    for (segment, &(velocity_min, velocity_max)) in segments.iter().zip(&GOLDEN_LADDER) {
        assert_eq!(segment.velocity_min.to_bits(), velocity_min.to_bits());
        assert_eq!(segment.velocity_max.to_bits(), velocity_max.to_bits());
        assert_eq!(
            segment.bc_value.to_bits(),
            GOLDEN_SEGMENT_BC.to_bits(),
            "corrected BC moved in band {}..{}: {} vs golden {}",
            segment.velocity_min,
            segment.velocity_max,
            segment.bc_value,
            GOLDEN_SEGMENT_BC
        );
    }
}

#[test]
fn matching_table_card_and_solve_are_bit_identical() {
    let (dir, file) = write_fixture("match", "bc5d_308.bin", &bc5d_fixture(0.308));
    let card = bridge_ok(
        "card.range_table",
        card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );
    let solve = bridge_ok(
        "solve",
        solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()}))),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    let adjustments = drop_adjustments(&card);
    assert_eq!(adjustments.len(), GOLDEN_CARD_DROP_ADJ.len());
    for (actual, expected) in adjustments.iter().zip(&GOLDEN_CARD_DROP_ADJ) {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "corrected card row moved: {actual} vs golden {expected}"
        );
    }

    let terminal = solve["samples"]
        .as_array()
        .expect("samples")
        .last()
        .expect("at least one sample")["drop_m"]
        .as_f64()
        .expect("drop_m");
    assert_eq!(
        terminal.to_bits(),
        GOLDEN_SOLVE_TERMINAL_DROP_M.to_bits(),
        "corrected solve moved: {terminal} vs golden {GOLDEN_SOLVE_TERMINAL_DROP_M}"
    );
}

/// A metric request states its diameter in MILLIMETRES; the guard must convert before
/// comparing, or every metric caller would be refused their own caliber's table.
#[test]
fn metric_diameter_matches_the_same_table() {
    let (dir, file) = write_fixture("metric", "bc5d_308.bin", &bc5d_fixture(0.308));
    let baseline = bridge_ok("card.range_table", metric_card_request(json!({})));
    let corrected = bridge_ok(
        "card.range_table",
        metric_card_request(json!({"bc5d_table_path": file.to_str().unwrap()})),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    let baseline_drop = *drop_adjustments(&baseline).last().unwrap();
    let corrected_drop = *drop_adjustments(&corrected).last().unwrap();
    // The 0.8 correction degrades BC 0.4 -> 0.32, so the far row needs more dial-up.
    assert!(
        corrected_drop > baseline_drop + 0.3,
        "a 7.8232 mm request must still be corrected by the .308 table: \
         baseline {baseline_drop} mil vs corrected {corrected_drop} mil"
    );
}

// ---------------------------------------------------------------------------
// A MISMATCHED table is refused, on every surface
// ---------------------------------------------------------------------------

/// The message every surface must carry: both calibers, named.
const EXPECTED_TEXT: &str = "table is for 0.224, shot is 0.308";

#[test]
fn card_commands_refuse_a_foreign_caliber_table() {
    let (dir, file) = write_fixture("card-mismatch", "bc5d_224.bin", &bc5d_fixture(0.224));
    let path = file.to_str().unwrap().to_owned();

    for command in ["card.range_table", "card.come_ups", "card.wind"] {
        let out = bridge_raw(
            command,
            card_request(json!({"bc5d_table_path": path, "wind_speeds": [10.0]})),
        );
        assert_eq!(out["ok"], false, "{command} accepted a .224 table: {out}");
        assert_eq!(out["error"]["code"], "command_failed", "{out}");
        let message = out["error"]["message"].as_str().unwrap();
        assert!(
            message.contains(EXPECTED_TEXT),
            "{command} must name both calibers: {message}"
        );
    }

    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn solve_refuses_a_foreign_caliber_table() {
    let (dir, file) = write_fixture("solve-mismatch", "bc5d_224.bin", &bc5d_fixture(0.224));
    let out = bridge_raw(
        "solve",
        solve_request(Some(json!({"bc5d_table_path": file.to_str().unwrap()}))),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(out["ok"], false, "solve accepted a .224 table: {out}");
    assert_eq!(out["error"]["code"], "command_failed", "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["code"], "invalid_value",
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"], "$.corrections.bc5d_table_path",
        "{out}"
    );
    let message = out["error"]["details"]["error"]["message"]
        .as_str()
        .unwrap();
    assert!(message.contains(EXPECTED_TEXT), "{message}");
}

/// The refusal must be a REFUSAL, not a downgrade: an uncorrected answer would be
/// indistinguishable from a corrected one to the caller, which is how the original defect
/// stayed invisible.
#[test]
fn a_refused_table_never_yields_an_uncorrected_answer() {
    let (dir, file) = write_fixture("no-downgrade", "bc5d_224.bin", &bc5d_fixture(0.224));
    let path = file.to_str().unwrap().to_owned();
    let uncorrected = bridge_ok("card.range_table", card_request(json!({})));
    let mismatched = bridge_raw(
        "card.range_table",
        card_request(json!({"bc5d_table_path": path})),
    );
    std::fs::remove_dir_all(&dir).unwrap();

    assert!(uncorrected["rows"].as_array().unwrap().len() > 1);
    assert!(
        mismatched.get("result").is_none(),
        "a mismatched table must not produce rows at all: {mismatched}"
    );
}

/// `bc5d.info` describes a file rather than judging a shot, so it still succeeds — but it
/// must report the identity the guard compares, so an app can pre-check and say something
/// friendlier than the refusal.
#[test]
fn bc5d_info_reports_the_key_the_guard_compares() {
    let (dir, file) = write_fixture("info", "bc5d_224.bin", &bc5d_fixture(0.224));
    let result = bridge_ok("bc5d.info", json!({"path": file.to_str().unwrap()}));
    std::fs::remove_dir_all(&dir).unwrap();

    assert_eq!(result["valid"], true);
    assert_eq!(
        result["caliber_key"], 224,
        "the pre-check key must be the guard's integer: {result}"
    );
    // The raw f32 header value stays as it was (0.224 is not exact in f32).
    assert!((result["caliber"].as_f64().unwrap() - 0.224).abs() < 1e-4, "{result}");

    // What an app's pre-check looks like: same rule, same verdict as the guard.
    let shot_key = (0.308_f64 * 1000.0).round() as i64;
    assert_ne!(shot_key, result["caliber_key"].as_i64().unwrap());
}

// ---------------------------------------------------------------------------
// CLI: --bc-table-dir picks by filename, so a mislabeled file must be fatal
// ---------------------------------------------------------------------------

#[test]
fn cli_refuses_a_mislabeled_table_file() {
    // The bytes say .224; the filename says .308, which is what --bc-table-dir looks up.
    let (dir, _file) = write_fixture("cli-mislabeled", "bc5d_308.bin", &bc5d_fixture(0.224));
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args([
            "trajectory", "-v", "2800", "-b", "0.4", "-m", "168", "-d", "0.308", "--drag-model",
            "g1", "--max-range", "300", "-o", "json",
        ])
        .arg("--bc-table-dir")
        .arg(&dir)
        .output()
        .expect("run ballistics CLI");
    std::fs::remove_dir_all(&dir).unwrap();

    assert!(
        !output.status.success(),
        "a mislabeled table must fail the run, not be skipped: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(EXPECTED_TEXT) && stderr.contains("--bc-table-dir"),
        "the CLI error must name the flag and both calibers: {stderr}"
    );
}

/// Same directory, correct filename: the CLI still applies the table it was given, so the
/// guard has not broken the ordinary `--bc-table-dir` path.
#[test]
fn cli_still_applies_a_correctly_labeled_table() {
    let (dir, _file) = write_fixture("cli-ok", "bc5d_308.bin", &bc5d_fixture(0.308));
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args([
            "trajectory", "-v", "2800", "-b", "0.4", "-m", "168", "-d", "0.308", "--drag-model",
            "g1", "--max-range", "300", "-o", "json",
        ])
        .arg("--bc-table-dir")
        .arg(&dir)
        .output()
        .expect("run ballistics CLI");
    std::fs::remove_dir_all(&dir).unwrap();

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("BC5D Table: Loaded caliber 0.308"),
        "the matching table must still load: {stderr}"
    );
    assert!(
        stderr.contains("BC5D Table: Generated"),
        "the matching table must still generate its ladder: {stderr}"
    );
}
