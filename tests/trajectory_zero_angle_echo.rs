//! MBA-1402: `trajectory` echoes the auto-zero-solved bore angle (degrees) in its summary,
//! JSON, and CSV output -- additive only. When auto-zero did not run (a bare `--angle`), the
//! new key/row/line must be absent and every other byte of the output must be unchanged.

use serde_json::Value;
use std::process::Command;

fn run_trajectory(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(args)
        .output()
        .expect("run trajectory command")
}

const BASE_ARGS: &[&str] = &[
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
    "300",
];

#[test]
fn json_echoes_zero_angle_only_when_auto_zero_ran() {
    // No auto-zero: the key must be absent entirely (skip_serializing_if), not null.
    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--angle", "1.0", "--output", "json"]);
    let output = run_trajectory(&args);
    assert!(output.status.success(), "{:?}", output);
    let doc: Value = serde_json::from_slice(&output.stdout).expect("valid JSON");
    assert!(
        doc.get("zero_angle_degrees").is_none(),
        "zero_angle_degrees must be absent without --auto-zero: {doc}"
    );

    // With auto-zero: the key is present and matches the angle `zero` itself would solve for
    // the same load/target-distance.
    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--auto-zero", "200", "--output", "json"]);
    let output = run_trajectory(&args);
    assert!(output.status.success(), "{:?}", output);
    let doc: Value = serde_json::from_slice(&output.stdout).expect("valid JSON");
    let echoed = doc["zero_angle_degrees"]
        .as_f64()
        .expect("zero_angle_degrees present and numeric when --auto-zero ran");

    let zero_output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args([
            "zero",
            "--velocity",
            "2700",
            "--bc",
            "0.475",
            "--mass",
            "168",
            "--diameter",
            "0.308",
            "--target-distance",
            "200",
            "--output",
            "json",
        ])
        .output()
        .expect("run zero command");
    let zero_doc: Value = serde_json::from_slice(&zero_output.stdout).expect("valid JSON");
    let expected = zero_doc["zero_angle_degrees"].as_f64().unwrap();

    assert!(
        (echoed - expected).abs() < 1e-9,
        "trajectory's echoed auto-zero angle ({echoed}) disagrees with zero's own solve \
         ({expected}) for the same load/target-distance"
    );
}

#[test]
fn csv_summary_echoes_zero_angle_only_when_auto_zero_ran() {
    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--angle", "1.0", "--output", "csv"]);
    let output = run_trajectory(&args);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        !stdout.contains("zero_angle_degrees"),
        "bare --angle CSV summary must not mention zero_angle_degrees:\n{stdout}"
    );

    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--auto-zero", "200", "--output", "csv"]);
    let output = run_trajectory(&args);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.lines().any(|line| line.starts_with("zero_angle_degrees,")),
        "auto-zero CSV summary must include a zero_angle_degrees row:\n{stdout}"
    );
}

#[test]
fn table_prints_zero_angle_only_when_auto_zero_ran() {
    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--angle", "1.0"]);
    let output = run_trajectory(&args);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        !stdout.contains("Zero Angle:"),
        "bare --angle table must not print a Zero Angle row:\n{stdout}"
    );

    let mut args: Vec<&str> = BASE_ARGS.to_vec();
    args.extend(["--auto-zero", "200"]);
    let output = run_trajectory(&args);
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("Zero Angle:"),
        "auto-zero table must print a Zero Angle row:\n{stdout}"
    );
}

/// Additive-only contract (MBA-1402): with neither the echo (--auto-zero absent) nor the
/// inverse mode requested, `trajectory` output must be BYTE-IDENTICAL to before this feature
/// existed. Pinned here for JSON and CSV; `tests/legacy_trajectory_json_contract.rs` already
/// pins the exact top-level JSON key set for a no-auto-zero run and continues to pass
/// unmodified, which is the strongest form of this guarantee for JSON.
#[test]
fn no_auto_zero_output_is_byte_identical_across_repeated_runs() {
    for format in ["json", "csv"] {
        let mut args: Vec<&str> = BASE_ARGS.to_vec();
        args.extend(["--angle", "1.0", "--output", format]);
        let first = run_trajectory(&args);
        let second = run_trajectory(&args);
        assert!(first.status.success() && second.status.success());
        assert_eq!(
            first.stdout, second.stdout,
            "{format} output was not deterministic/byte-identical across runs"
        );
    }
}
