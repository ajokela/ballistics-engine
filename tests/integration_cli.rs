use std::process::Command;
use std::path::PathBuf;

fn get_cli_binary() -> PathBuf {
    // Try to find the built binary
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("ballistics-cli");
    
    if !path.exists() {
        // Try release build
        path.pop();
        path.pop();
        path.push("release");
        path.push("ballistics-cli");
    }
    
    path
}

#[test]
fn test_cli_trajectory_basic() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity", "2700",
            "--bc", "0.475",
            "--mass", "168",
            "--diameter", "0.308",
            "--max-range", "2000"
        ])
        .output()
        .expect("Failed to execute command");
    
    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("TRAJECTORY") || stdout.contains("Range"), 
            "Should contain trajectory output");
}



#[test]
fn test_cli_monte_carlo_command() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "monte-carlo",
            "--velocity", "2700",
            "--bc", "0.475",
            "--mass", "168",
            "--diameter", "0.308",
            "--num-sims", "10"
        ])
        .output()
        .expect("Failed to execute command");
    
    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Monte carlo output contains statistical results
    assert!(stdout.contains("Mean") || stdout.contains("Impact") || stdout.contains("Trajectory"), 
            "Should contain Monte Carlo results: {}", stdout);
}

#[test]
fn test_cli_help() {
    let output = Command::new(get_cli_binary())
        .args(&[ "--help"])
        .output()
        .expect("Failed to execute command");
    
    assert!(output.status.success(), "Help command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("trajectory"), "Should list trajectory command");
    assert!(stdout.contains("monte-carlo"), "Should list monte-carlo command");
    assert!(stdout.contains("info"), "Should list info command");
}

#[test]
fn test_cli_invalid_command() {
    let output = Command::new(get_cli_binary())
        .args(&[ "invalid-command"])
        .output()
        .expect("Failed to execute command");
    
    // Command should fail for invalid subcommand
    assert!(!output.status.success(), "Invalid command should fail");
}

#[test]
fn test_cli_missing_required_args() {
    let output = Command::new(get_cli_binary())
        .args(&[ "trajectory"])
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
            "--velocity", "2700",
            "--bc", "0.475",
            "--mass", "168",
            "--diameter", "0.308",
            "--max-range", "1000",
            "--output", "json"
        ])
        .output()
        .expect("Failed to execute command");
    
    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // JSON output should contain brackets
    assert!(stdout.contains("[") || stdout.contains("{"), "Should be JSON format");
}

#[test]
fn test_cli_output_format_csv() {
    let output = Command::new(get_cli_binary())
        .args(&[
            "trajectory",
            "--velocity", "2700",
            "--bc", "0.245",
            "--mass", "168",
            "--diameter", "0.308",
            "--max-range", "1000",
            "--output", "csv"
        ])
        .output()
        .expect("Failed to execute command");
    
    assert!(output.status.success(), "Command should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    // CSV output should contain commas
    assert!(stdout.contains(","), "Should be CSV format");
}