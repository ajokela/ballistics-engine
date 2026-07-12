use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn get_cli_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/debug/ballistics")
}

fn unique_temp_dir(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("ballistics-{label}-{}-{nonce}", std::process::id()))
}

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

fn bc5d_fixture() -> Vec<u8> {
    let weight_bins = [168.0_f32];
    let bc_bins = [0.4_f32, 0.5_f32];
    let muzzle_bins = [2_800.0_f32];
    let current_bins = [2_000.0_f32];
    // Correction depends only on the BC axis: raw 0.4 -> 0.8, adjusted 0.5 -> 1.2.
    let data = [0.8_f32, 1.2_f32];

    let mut checksum_data = Vec::new();
    push_f32s(&mut checksum_data, &weight_bins);
    push_f32s(&mut checksum_data, &bc_bins);
    push_f32s(&mut checksum_data, &muzzle_bins);
    push_f32s(&mut checksum_data, &current_bins);
    push_f32s(&mut checksum_data, &data);

    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"BC5D");
    bytes.extend_from_slice(&2_u32.to_le_bytes());
    bytes.extend_from_slice(&0.308_f32.to_le_bytes());
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

fn run_with_table(path: &Path) -> Output {
    Command::new(get_cli_binary())
        .args([
            "trajectory",
            "-v",
            "2800",
            "-b",
            "0.4",
            "--bc-adjustment",
            "1.25",
            "-m",
            "168",
            "-d",
            "0.308",
            "--bullet-length",
            "1.2",
            "--drag-model",
            "g1",
            "--max-range",
            "5",
            "-o",
            "json",
        ])
        .arg("--bc-table-dir")
        .arg(path)
        .output()
        .expect("run ballistics")
}

fn ladder_command(metric: bool, velocity: &str) -> Command {
    let mut command = Command::new(get_cli_binary());
    if metric {
        command.args(["--units", "metric"]);
    }
    command.arg("trajectory");
    if metric {
        command.args([
            "-v",
            velocity,
            "-b",
            "0.4",
            "-m",
            "10.886",
            "-d",
            "7.8232",
            "--bullet-length",
            "30.48",
        ]);
    } else {
        command.args([
            "-v",
            velocity,
            "-b",
            "0.4",
            "-m",
            "168",
            "-d",
            "0.308",
            "--bullet-length",
            "1.2",
        ]);
    }
    command.args(["--drag-model", "g1", "--max-range", "5", "-o", "json"]);
    command
}

fn print_ladder(path: &Path, metric: bool, velocity: &str) -> Output {
    ladder_command(metric, velocity)
        .arg("--bc-table-dir")
        .arg(path)
        .arg("--print-bc-segments")
        .output()
        .expect("print BC5D ladder")
}

fn ladder_segments(output: &Output) -> Vec<String> {
    assert!(
        output.status.success(),
        "table command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stderr)
        .lines()
        .filter_map(|line| line.trim().strip_prefix("--bc-segment "))
        .map(str::to_owned)
        .collect()
}

fn assert_ladder_round_trips(directory: &Path, metric: bool, velocity: &str, expected_top: &str) {
    let printed = print_ladder(directory, metric, velocity);
    let segments = ladder_segments(&printed);
    assert!(!segments.is_empty(), "no ladder was printed");
    assert!(
        segments[0].starts_with(expected_top),
        "top band lost its narrow boundary: {}",
        segments[0]
    );
    for segment in &segments {
        let mut fields = segment.split(':');
        let min: f64 = fields.next().unwrap().parse().unwrap();
        let max: f64 = fields.next().unwrap().parse().unwrap();
        assert!(min < max, "printed degenerate segment: {segment}");
    }

    let mut pasted = ladder_command(metric, velocity);
    for segment in &segments {
        pasted.arg("--bc-segment").arg(segment);
    }
    let output = pasted.output().expect("paste BC5D ladder");
    assert!(
        output.status.success(),
        "printed ladder was rejected: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_adjusted_axis_used(output: Output) {
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("BC5D Table: BC range 0.60000 - 0.60000"),
        "segment ladder did not use adjusted BC 0.5:\n{stderr}"
    );
    assert!(
        stderr.contains("BC5D Table: Muzzle correction=1.2000 for BC=0.500"),
        "muzzle lookup did not use adjusted BC 0.5:\n{stderr}"
    );
    assert!(
        stderr.contains("table-corrected=0.6000, factor=1.2000"),
        "scalar fallback disagrees with the segment ladder:\n{stderr}"
    );
    assert!(!stderr.contains("Muzzle correction=0.8000"));
}

#[test]
fn bc_adjustment_and_bc5d_share_the_effective_bc_axis() {
    let directory = unique_temp_dir("bc5d-adjustment");
    fs::create_dir_all(&directory).unwrap();
    fs::write(directory.join("bc5d_308.bin"), bc5d_fixture()).unwrap();
    let output = run_with_table(&directory);
    fs::remove_dir_all(&directory).unwrap();

    assert_adjusted_axis_used(output);
}

#[test]
fn printed_imperial_ladder_preserves_narrow_muzzle_band() {
    let directory = unique_temp_dir("bc5d-print-imperial");
    fs::create_dir_all(&directory).unwrap();
    fs::write(directory.join("bc5d_308.bin"), bc5d_fixture()).unwrap();
    assert_ladder_round_trips(&directory, false, "2700.3", "2700:2700.3:");
    fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn printed_metric_ladder_preserves_narrow_muzzle_band() {
    let directory = unique_temp_dir("bc5d-print-metric");
    fs::create_dir_all(&directory).unwrap();
    fs::write(directory.join("bc5d_308.bin"), bc5d_fixture()).unwrap();
    assert_ladder_round_trips(&directory, true, "823", "822.96:823:");
    fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn printed_metric_ladder_survives_display_unit_float_collision() {
    let directory = unique_temp_dir("bc5d-print-metric-ulp");
    fs::create_dir_all(&directory).unwrap();
    fs::write(directory.join("bc5d_308.bin"), bc5d_fixture()).unwrap();
    assert_ladder_round_trips(
        &directory,
        true,
        "518.1600000000001",
        "518.1600000000001:518.1600000000002:",
    );
    fs::remove_dir_all(&directory).unwrap();
}
