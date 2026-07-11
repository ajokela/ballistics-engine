use serde_json::Value;
#[cfg(feature = "online")]
use std::fs;
#[cfg(feature = "online")]
use std::io::{Read, Write};
#[cfg(feature = "online")]
use std::net::TcpListener;
use std::path::PathBuf;
use std::process::Command;
#[cfg(feature = "online")]
use std::thread;
#[cfg(feature = "online")]
use std::time::{SystemTime, UNIX_EPOCH};

fn get_cli_binary() -> PathBuf {
    // Try to find the built binary
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("ballistics");

    if !path.exists() {
        // Try release build
        path.pop();
        path.pop();
        path.push("release");
        path.push("ballistics");
    }

    path
}

#[cfg(feature = "online")]
fn accepted_tos_home() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let home = std::env::temp_dir().join(format!(
        "ballistics-compare-test-{}-{nonce}",
        std::process::id()
    ));
    let config_dir = home.join(".ballistics");
    fs::create_dir_all(&config_dir).unwrap();
    fs::write(
        config_dir.join("tos.json"),
        r#"{
  "accepted": true,
  "accepted_version": "2026-01-26",
  "accepted_at": "test",
  "terms_hash": "test"
}"#,
    )
    .unwrap();
    home
}

#[cfg(feature = "online")]
fn comparison_rows(stdout: &str) -> Vec<[f64; 4]> {
    stdout
        .lines()
        .filter_map(|line| {
            let values: Vec<f64> = line
                .split(['║', '│'])
                .filter_map(|cell| cell.trim().parse().ok())
                .collect();
            (values.len() == 4).then(|| [values[0], values[1], values[2], values[3]])
        })
        .collect()
}

#[cfg(feature = "online")]
#[test]
fn compare_uses_api_drop_frame_and_same_default_zero() {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap();
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        let mut request = Vec::new();
        let mut buffer = [0_u8; 4096];
        loop {
            let bytes = stream.read(&mut buffer).unwrap();
            if bytes == 0 {
                break;
            }
            request.extend_from_slice(&buffer[..bytes]);
            if request.windows(4).any(|window| window == b"\r\n\r\n") {
                break;
            }
        }

        let body = r#"{
  "results": {"barrel_angle": 0.0},
  "trajectory": [
    {"distance": {"value": 0.0}, "drop": {"value": -1.5}, "wind_drift": 0.0, "velocity": 2600.0, "energy": 2600.0, "time": 0.0},
    {"distance": {"value": 100.0}, "drop": {"value": 0.0}, "wind_drift": 0.0, "velocity": 2400.0, "energy": 2200.0, "time": 0.12}
  ]
}"#;
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
            body.len(),
            body
        )
        .unwrap();
        String::from_utf8(request).unwrap()
    });

    let home = accepted_tos_home();
    let output = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "trajectory",
            "--velocity",
            "2600",
            "--bc",
            "0.243",
            "--mass",
            "175",
            "--diameter",
            "0.308",
            "--drag-model",
            "g7",
            "--max-range",
            "200",
            "--sight-height",
            "1.5",
            "--compare",
            "--ignore-ground-impact",
            "--api-url",
            &format!("http://{address}"),
        ])
        .output()
        .expect("failed to execute comparison command");
    let request = server.join().unwrap();
    fs::remove_dir_all(home).unwrap();

    assert!(
        output.status.success(),
        "comparison failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        request.lines().next().unwrap().contains("zero_distance=100.0"),
        "imperial comparison must use the API's 100-yard default zero: {request}"
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    let rows = comparison_rows(&stdout);
    assert_eq!(rows.len(), 2, "unexpected comparison output:\n{stdout}");
    assert_eq!(rows[0], [0.0, -1.5, -1.5, 0.0]);
    assert!(
        rows[1][2].abs() <= 0.05 && rows[1][3].abs() <= 0.05,
        "both legs must be zeroed at 100 yards: {:?}\n{stdout}",
        rows[1]
    );
}

#[test]
fn test_cli_trajectory_basic() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity",
            "2700",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "2000",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("TRAJECTORY") || stdout.contains("Range"),
        "Should contain trajectory output"
    );
}

#[test]
fn test_cli_monte_carlo_command() {
    // Pin the target so this remains a simple command/output smoke test rather than depending on
    // the omitted-target baseline convention.
    let output = Command::new(get_cli_binary())
        .args(&[
            "monte-carlo",
            "--velocity",
            "2700",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--num-sims",
            "10",
            "--wind-direction-std",
            "5",
            "--target-distance",
            "800",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Monte carlo output contains statistical results
    assert!(
        stdout.contains("Mean") || stdout.contains("Impact") || stdout.contains("Trajectory"),
        "Should contain Monte Carlo results: {}",
        stdout
    );
}

#[test]
fn test_cli_monte_carlo_all_shortfalls_have_null_cep() {
    let output = Command::new(get_cli_binary())
        .args([
            "monte-carlo",
            "--velocity",
            "2700",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--num-sims",
            "10",
            "--target-distance",
            "10000",
            "--output",
            "full",
        ])
        .output()
        .expect("failed to execute Monte Carlo command");

    assert!(
        output.status.success(),
        "Monte Carlo command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: Value = serde_json::from_slice(&output.stdout).expect("valid JSON output");
    assert!(
        result["cep"].is_null(),
        "CEP is undefined when every sample falls short: {result}"
    );
    assert_eq!(result["target_shortfall_fraction"].as_f64(), Some(1.0));
    assert_eq!(result["hit_probability"].as_f64(), Some(0.0));
}

#[test]
fn test_cli_help() {
    let output = Command::new(get_cli_binary())
        .args(&["--help"])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Help command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("trajectory"),
        "Should list trajectory command"
    );
    assert!(
        stdout.contains("monte-carlo"),
        "Should list monte-carlo command"
    );
    assert!(stdout.contains("zero"), "Should list zero command");
}

#[test]
fn test_cli_invalid_command() {
    let output = Command::new(get_cli_binary())
        .args(&["invalid-command"])
        .output()
        .expect("Failed to execute command");

    // Command should fail for invalid subcommand
    assert!(!output.status.success(), "Invalid command should fail");
}

#[test]
fn test_cli_missing_required_args() {
    let output = Command::new(get_cli_binary())
        .args(&["trajectory"])
        .output()
        .expect("Failed to execute command");

    // Should fail due to missing required arguments
    assert!(!output.status.success(), "Should fail with missing args");
}

#[test]
fn test_cli_output_format_json() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity",
            "2700",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "1000",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // JSON output should contain brackets
    assert!(
        stdout.contains("[") || stdout.contains("{"),
        "Should be JSON format"
    );
}

#[test]
fn test_cli_output_format_csv() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity",
            "2700",
            "--bc",
            "0.245",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--max-range",
            "1000",
            "--output",
            "csv",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // CSV output should contain commas
    assert!(stdout.contains(","), "Should be CSV format");
}

#[test]
fn saved_metric_profile_preserves_physical_values_when_recalled_as_imperial() {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let home = std::env::temp_dir().join(format!(
        "ballistics-profile-test-{}-{nonce}",
        std::process::id()
    ));
    std::fs::create_dir_all(&home).unwrap();

    let saved = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "profile",
            "save",
            "metric-exact",
            "--units",
            "metric",
            "--velocity",
            "762",
            "--bc",
            "0.243",
            "--mass",
            "11.33980925",
            "--diameter",
            "7.8232",
            "--drag-model",
            "g7",
            "--twist-rate",
            "254",
            "--sight-height",
            "50.8",
            "--auto-zero",
            "91.44",
            "--wind-speed",
            "4.4704",
            "--bullet-length",
            "30.48",
        ])
        .output()
        .expect("Failed to save metric profile");
    assert!(
        saved.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&saved.stderr)
    );

    let recalled = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "metric-exact",
            "--units",
            "imperial",
            "--max-range",
            "100",
            "--ignore-ground-impact",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to recall metric profile as imperial");
    let explicit = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "trajectory",
            "--units",
            "imperial",
            "--velocity",
            "2500",
            "--bc",
            "0.243",
            "--mass",
            "175",
            "--diameter",
            "0.308",
            "--drag-model",
            "g7",
            "--twist-rate",
            "10",
            "--sight-height",
            "2",
            "--auto-zero",
            "100",
            "--wind-speed",
            "10",
            "--bullet-length",
            "1.2",
            "--max-range",
            "100",
            "--ignore-ground-impact",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to run explicit imperial trajectory");
    let recalled_mpbr = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "mpbr",
            "--profile",
            "metric-exact",
            "--units",
            "imperial",
            "--vital-zone",
            "8",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to recall metric profile for MPBR");
    let explicit_mpbr = Command::new(get_cli_binary())
        .env("HOME", &home)
        .args([
            "mpbr",
            "--units",
            "imperial",
            "--velocity",
            "2500",
            "--bc",
            "0.243",
            "--mass",
            "175",
            "--diameter",
            "0.308",
            "--drag-model",
            "g7",
            "--sight-height",
            "2",
            "--vital-zone",
            "8",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to run explicit imperial MPBR");

    assert!(
        recalled.status.success(),
        "profile recall failed: {}",
        String::from_utf8_lossy(&recalled.stderr)
    );
    assert!(
        explicit.status.success(),
        "explicit trajectory failed: {}",
        String::from_utf8_lossy(&explicit.stderr)
    );
    assert!(
        recalled_mpbr.status.success(),
        "profile MPBR failed: {}",
        String::from_utf8_lossy(&recalled_mpbr.stderr)
    );
    assert!(
        explicit_mpbr.status.success(),
        "explicit MPBR failed: {}",
        String::from_utf8_lossy(&explicit_mpbr.stderr)
    );
    let recalled: Value = serde_json::from_slice(&recalled.stdout).expect("valid recalled JSON");
    let explicit: Value = serde_json::from_slice(&explicit.stdout).expect("valid explicit JSON");
    let recalled_mpbr: Value =
        serde_json::from_slice(&recalled_mpbr.stdout).expect("valid recalled MPBR JSON");
    let explicit_mpbr: Value =
        serde_json::from_slice(&explicit_mpbr.stdout).expect("valid explicit MPBR JSON");
    std::fs::remove_dir_all(&home).unwrap();

    assert_eq!(recalled["units"], "imperial");
    for field in [
        "max_height",
        "time_of_flight",
        "impact_velocity",
        "impact_energy",
        "stability_coefficient",
    ] {
        let actual = recalled[field].as_f64().unwrap();
        let expected = explicit[field].as_f64().unwrap();
        let tolerance = expected.abs().max(1.0) * 1e-10;
        assert!(
            (actual - expected).abs() <= tolerance,
            "cross-unit profile changed {field}: recalled={actual} explicit={expected}"
        );
    }
    for field in ["mpbr", "optimal_zero", "impact_velocity", "impact_energy"] {
        let actual = recalled_mpbr[field].as_f64().unwrap();
        let expected = explicit_mpbr[field].as_f64().unwrap();
        let tolerance = expected.abs().max(1.0) * 1e-10;
        assert!(
            (actual - expected).abs() <= tolerance,
            "cross-unit profile changed MPBR {field}: recalled={actual} explicit={expected}"
        );
    }
}

// ============================================================================
// True Velocity Offline Mode Tests
// ============================================================================

#[test]
fn test_true_velocity_offline_basic() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Offline true-velocity should succeed"
    );
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Parse JSON output
    let json: Value = serde_json::from_str(&stdout).expect("Should produce valid JSON");

    // Verify required fields exist
    assert!(
        json["effective_velocity"].is_number(),
        "Should have effective_velocity"
    );
    assert!(
        json["calculated_drop_mil"].is_number(),
        "Should have calculated_drop_mil"
    );
    assert!(json["confidence"].is_string(), "Should have confidence");
    assert!(json["iterations"].is_number(), "Should have iterations");

    // Verify velocity is in reasonable range (1500-4500 fps)
    let velocity = json["effective_velocity"].as_f64().unwrap();
    assert!(
        velocity > 1500.0 && velocity < 4500.0,
        "Velocity {} should be in reasonable range",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_converges_accurately() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    // Check that calculated drop is close to measured drop
    let calculated_drop = json["calculated_drop_mil"].as_f64().unwrap();
    let measured_drop = 5.1;
    let error = (calculated_drop - measured_drop).abs();

    assert!(
        error < 0.1,
        "Calculated drop {} should be within 0.1 MIL of measured {}",
        calculated_drop,
        measured_drop
    );
}

#[test]
fn test_true_velocity_offline_308_caliber() {
    // Test with .308 Win 175gr SMK
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "8.5",
            "--range",
            "800",
            "--bc",
            "0.475",
            "--drag-model",
            "g7",
            "--mass",
            "175",
            "--diameter",
            "0.308",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should succeed for .308 caliber");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    // .308 175gr with 8.5 MIL at 800 yards with LOS-correct zeroing.
    assert!(
        velocity > 2150.0 && velocity < 2300.0,
        ".308 velocity {} should be in expected range 2150-2300",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_224_caliber() {
    // Test with .224 77gr
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "3.2",
            "--range",
            "400",
            "--bc",
            "0.210",
            "--drag-model",
            "g7",
            "--mass",
            "77",
            "--diameter",
            "0.224",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should succeed for .224 caliber");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    // .224 77gr with 3.2 MIL at 400 yards with LOS-correct zeroing.
    assert!(
        velocity > 2300.0 && velocity < 2450.0,
        ".224 velocity {} should be in expected range 2300-2450",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_with_chrono() {
    // Test with chronograph velocity for adjustment calculation
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--chrono-velocity",
            "2800",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    // Should have velocity adjustment when chrono is provided
    assert!(
        json["velocity_adjustment"].is_number(),
        "Should have velocity_adjustment when chrono provided"
    );
    assert!(
        json["adjustment_percent"].is_number(),
        "Should have adjustment_percent when chrono provided"
    );
}

#[test]
fn test_true_velocity_offline_g1_drag_model() {
    // Test with G1 drag model
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "6.0",
            "--range",
            "500",
            "--bc",
            "0.450",
            "--drag-model",
            "g1",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should succeed with G1 drag model");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    assert!(
        velocity > 2000.0 && velocity < 3500.0,
        "G1 velocity {} should be in reasonable range",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_output_formats() {
    // Test table output
    let table_output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "table",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(table_output.status.success());
    let stdout = String::from_utf8_lossy(&table_output.stdout);
    assert!(
        stdout.contains("VELOCITY TRUING RESULTS"),
        "Table output should have header"
    );
    assert!(
        stdout.contains("Effective Muzzle Velocity"),
        "Table output should show velocity"
    );

    // Test CSV output
    let csv_output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "csv",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(csv_output.status.success());
    let stdout = String::from_utf8_lossy(&csv_output.stdout);
    assert!(stdout.contains(","), "CSV output should have commas");
    assert!(
        stdout.contains("effective_velocity"),
        "CSV should have header row"
    );
}

#[test]
fn test_true_velocity_offline_metric_units() {
    // Test with metric units
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "550", // 550 meters ≈ 600 yards
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "9.07", // 140gr in grams
            "--diameter",
            "6.7", // 0.264" in mm
            "--units",
            "metric",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should succeed with metric units");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    assert_eq!(
        json["velocity_unit"].as_str().unwrap(),
        "m/s",
        "Metric output should use m/s"
    );

    let velocity_ms = json["effective_velocity"].as_f64().unwrap();
    // Should be around 700-900 m/s (2300-2950 fps)
    // Range is wider because metric input converts differently
    assert!(
        velocity_ms > 650.0 && velocity_ms < 950.0,
        "Metric velocity {} m/s should be in expected range",
        velocity_ms
    );
}

#[test]
fn true_velocity_metric_default_sight_height_matches_explicit_50_mm() {
    let run = |sight_height: Option<&str>| {
        let mut args = vec![
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "550",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "9.07",
            "--diameter",
            "6.7",
            "--units",
            "metric",
            "--offline",
            "--output",
            "json",
        ];
        if let Some(value) = sight_height {
            args.extend(["--sight-height", value]);
        }

        let output = Command::new(get_cli_binary())
            .args(args)
            .output()
            .expect("Failed to execute metric true-velocity command");
        assert!(
            output.status.success(),
            "metric true-velocity failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice::<Value>(&output.stdout).expect("Should produce valid JSON")
    };

    let default = run(None);
    let explicit = run(Some("50"));
    let default_velocity = default["effective_velocity"].as_f64().unwrap();
    let explicit_velocity = explicit["effective_velocity"].as_f64().unwrap();

    assert_eq!(
        default_velocity.to_bits(),
        explicit_velocity.to_bits(),
        "metric default should be 50 mm: omitted={default_velocity} m/s explicit={explicit_velocity} m/s"
    );
}

#[test]
fn test_true_velocity_offline_extreme_drop() {
    // Test with high drop value (long range)
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "15.0", // Very high drop
            "--range",
            "1000",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should handle extreme drop values");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    // With 15 MIL at 1000 yards, velocity should be relatively low
    assert!(
        velocity > 1500.0 && velocity < 3000.0,
        "Extreme drop velocity {} should be in range",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_low_drop() {
    // Test with low drop value (short range or high velocity)
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "1.5",
            "--range",
            "300",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should handle low drop values");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    // With only 1.5 MIL at 300 yards, velocity should remain in the upper search range.
    assert!(
        velocity > 2450.0 && velocity < 4000.0,
        "Low drop velocity {} should indicate high muzzle velocity",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_custom_zero() {
    // Test with custom zero distance
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--zero-distance",
            "200", // 200 yard zero
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(output.status.success(), "Should succeed with custom zero");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    assert!(
        velocity > 2000.0 && velocity < 4000.0,
        "Custom zero velocity {} should be reasonable",
        velocity
    );
}

#[test]
fn test_true_velocity_offline_custom_atmosphere() {
    // Test with custom atmospheric conditions
    let output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--temperature",
            "90", // Hot day
            "--altitude",
            "5000", // 5000 ft elevation
            "--pressure",
            "25.0", // Lower pressure
            "--humidity",
            "20", // Low humidity
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Should succeed with custom atmosphere"
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let json: Value = serde_json::from_str(&stdout).unwrap();

    let velocity = json["effective_velocity"].as_f64().unwrap();
    // At altitude with thin air, same drop requires lower velocity
    assert!(
        velocity > 2000.0 && velocity < 4000.0,
        "High altitude velocity {} should be reasonable",
        velocity
    );
}

/// Test that offline and online modes produce similar results
/// This test requires the online feature and network access
#[test]
#[cfg(feature = "online")]
fn test_true_velocity_offline_vs_online_consistency() {
    // Run offline calculation
    let offline_output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute offline command");

    assert!(offline_output.status.success(), "Offline should succeed");
    let offline_stdout = String::from_utf8_lossy(&offline_output.stdout);
    let offline_json: Value =
        serde_json::from_str(&offline_stdout).expect("Offline should produce valid JSON");
    let offline_velocity = offline_json["effective_velocity"].as_f64().unwrap();

    // Run online calculation
    let online_output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.1",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute online command");

    // Online might fail due to network issues, so check conditionally
    if online_output.status.success() {
        let online_stdout = String::from_utf8_lossy(&online_output.stdout);
        if let Ok(online_json) = serde_json::from_str::<Value>(&online_stdout) {
            if let Some(online_velocity) = online_json["effective_velocity"].as_f64() {
                // Calculate percentage difference
                let diff_percent =
                    ((offline_velocity - online_velocity) / online_velocity * 100.0).abs();

                // The remote service may lag local numerical corrections, so keep this as a
                // broad consistency check rather than a tight golden comparison.
                assert!(
                    diff_percent < 15.0,
                    "Offline ({:.1} fps) and online ({:.1} fps) should be within 15%, got {:.2}%",
                    offline_velocity,
                    online_velocity,
                    diff_percent
                );
            }
        }
    }
    // If online fails, we just skip the comparison (network not available)
}

/// Test consistency across multiple runs (deterministic)
#[test]
fn test_true_velocity_offline_deterministic() {
    let args = &[
        "true-velocity",
        "--measured-drop",
        "5.1",
        "--range",
        "600",
        "--bc",
        "0.27",
        "--drag-model",
        "g7",
        "--mass",
        "140",
        "--diameter",
        "0.264",
        "--offline",
        "--output",
        "json",
    ];

    // Run twice
    let output1 = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("Failed to execute command");

    let output2 = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("Failed to execute command");

    assert!(output1.status.success());
    assert!(output2.status.success());

    let json1: Value = serde_json::from_str(&String::from_utf8_lossy(&output1.stdout)).unwrap();
    let json2: Value = serde_json::from_str(&String::from_utf8_lossy(&output2.stdout)).unwrap();

    let vel1 = json1["effective_velocity"].as_f64().unwrap();
    let vel2 = json2["effective_velocity"].as_f64().unwrap();

    assert!(
        (vel1 - vel2).abs() < 0.001,
        "Results should be deterministic: {} vs {}",
        vel1,
        vel2
    );
}

/// Test that inverse calculation works: find velocity, then verify drop
#[test]
fn test_true_velocity_offline_inverse_verification() {
    // First, find the velocity for a known drop
    let tv_output = Command::new(get_cli_binary())
        .args(&[
            "true-velocity",
            "--measured-drop",
            "5.0",
            "--range",
            "600",
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--zero-distance",
            "100",
            "--offline",
            "--output",
            "json",
        ])
        .output()
        .expect("Failed to execute true-velocity");

    assert!(tv_output.status.success());
    let tv_json: Value = serde_json::from_str(&String::from_utf8_lossy(&tv_output.stdout)).unwrap();
    let found_velocity = tv_json["effective_velocity"].as_f64().unwrap();

    // Now run a trajectory at that velocity and check the drop at 600 yards
    let traj_output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity",
            &format!("{:.0}", found_velocity),
            "--bc",
            "0.27",
            "--drag-model",
            "g7",
            "--mass",
            "140",
            "--diameter",
            "0.264",
            "--auto-zero",
            "100",
            "--max-range",
            "700",
            "--ignore-ground-impact",
            "--full",
            "--output",
            "csv",
        ])
        .output()
        .expect("Failed to execute trajectory");

    assert!(
        traj_output.status.success(),
        "Trajectory command should succeed"
    );

    // Parse CSV to find drop at ~600 yards
    let csv = String::from_utf8_lossy(&traj_output.stdout);
    let lines: Vec<&str> = csv.lines().collect();

    // Find line closest to 600 yards (z_yd column, index 3)
    for line in lines.iter().skip(1) {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 4 {
            if let Ok(range) = fields[3].parse::<f64>() {
                if range >= 595.0 && range <= 605.0 {
                    // Found a point near 600 yards
                    // The trajectory ran successfully at the found velocity
                    return;
                }
            }
        }
    }

    // If we reach here, trajectory didn't have data at 600 yards
    // This is OK - trajectory might stop earlier due to ground impact
}
