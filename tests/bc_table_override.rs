//! Manual BC segments must replace both the segment ladder and scalar fallback from a table.

use serde_json::Value;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn get_cli_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/debug/ballistics")
}

fn unique_temp_dir() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "ballistics-bc-table-override-{}-{nonce}",
        std::process::id()
    ))
}

fn push_f32(output: &mut Vec<u8>, value: f32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn flat_table_fixture() -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"BCCR");
    bytes.extend_from_slice(&1_u32.to_le_bytes());
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // flags
    for dimension in [1_u32; 5] {
        bytes.extend_from_slice(&dimension.to_le_bytes());
    }
    bytes.extend_from_slice(&0_u64.to_le_bytes()); // timestamp
    bytes.extend_from_slice(&0_u32.to_le_bytes()); // legacy no-checksum sentinel
    bytes.extend_from_slice(&[0_u8; 16]);
    push_f32(&mut bytes, 0.4); // BC bin
    push_f32(&mut bytes, 168.0); // mass bin (grains)
    push_f32(&mut bytes, 1.0); // length bin (inches)
    push_f32(&mut bytes, 2_800.0); // velocity bin (fps)
    push_f32(&mut bytes, 0.8); // correction
    assert_eq!(bytes.len(), 80);
    bytes
}

fn run(table: Option<&Path>, manual_segments: bool) -> Output {
    let mut command = Command::new(get_cli_binary());
    command.args([
        "trajectory",
        "-v",
        "2800",
        "-b",
        "0.4",
        "--bc-adjustment",
        "1.1",
        "-m",
        "168",
        "-d",
        "0.308",
        "--bullet-length",
        "1.0",
        "--drag-model",
        "g1",
        "--max-range",
        "1000",
        "--ignore-ground-impact",
        "-o",
        "json",
    ]);
    if manual_segments {
        // The gap between these bands exercises the scalar fallback. Velocities outside the
        // global coverage clamp to the nearest band and would not expose the leaked correction.
        command.args([
            "--bc-segment",
            "500:600:0.2",
            "--bc-segment",
            "3000:3100:0.2",
        ]);
    }
    if let Some(path) = table {
        command.arg("--bc-table").arg(path);
    }
    command.output().expect("run ballistics")
}

fn result(output: Output) -> (f64, f64, String) {
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: Value = serde_json::from_slice(&output.stdout).expect("valid trajectory JSON");
    (
        json["impact_velocity"].as_f64().expect("impact velocity"),
        json["time_of_flight"].as_f64().expect("time of flight"),
        String::from_utf8(output.stderr).expect("UTF-8 diagnostics"),
    )
}

#[test]
fn manual_bc_segments_fully_override_flat_bc_table() {
    let directory = unique_temp_dir();
    fs::create_dir_all(&directory).unwrap();
    let table_path = directory.join("table.bin");
    fs::write(&table_path, flat_table_fixture()).unwrap();

    let manual = result(run(None, true));
    let combined = result(run(Some(&table_path), true));
    let table_only = result(run(Some(&table_path), false));
    fs::remove_dir_all(&directory).unwrap();

    assert!(
        table_only.2.contains("Muzzle correction=0.8000"),
        "fixture correction was not loaded:\n{}",
        table_only.2
    );
    assert!(
        (table_only.0 - manual.0).abs() > 10.0,
        "fixture must materially change the scalar BC"
    );
    assert!(
        (combined.0 - manual.0).abs() < 1e-6,
        "flat table leaked into manual fallback: manual={} fps, combined={} fps",
        manual.0,
        combined.0
    );
    assert!(
        (combined.1 - manual.1).abs() < 1e-9,
        "flat table changed manual time of flight: manual={} s, combined={} s",
        manual.1,
        combined.1
    );
}
