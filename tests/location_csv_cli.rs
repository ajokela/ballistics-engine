//! MBA-1425: a supplied `--location`/`--site` (and `--profile`/`--profile-row`) must either
//! load or stop the run — never silently fall back to defaults.
//!
//! Three silences existed, in ascending order of visibility:
//!   1. `--location <FILE>` with no `--site`: the file was never read at all — no message,
//!      no error, default atmosphere. (This one cost real debugging time during MBA-1421,
//!      helped along by the near-miss `--location-name` flag, which is only a PDF header
//!      override.)
//!   2. an unreadable path: a stderr warning, then the run continued on defaults.
//!   3. a site name not present in the file: same warn-and-continue.
//!
//! All three produced a clean, plausible trajectory computed against air the user never
//! chose. 1 is now a clap pairing error; 2 and 3 are hard errors.

use std::io::Write;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn base_args() -> Vec<&'static str> {
    vec![
        "trajectory", "-b", "0.243", "--drag-model", "g7", "-v", "2700", "-m", "175", "-d",
        "0.308", "--max-range", "300", "-o", "json",
    ]
}

fn location_csv(name: &str) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!("bx_1425_{name}.csv"));
    let mut f = std::fs::File::create(&path).expect("create csv");
    writeln!(f, "NAME,LATITUDE,ALTITUDE,TEMPERATURE,PRESSURE,HUMIDITY").unwrap();
    writeln!(f, "TestSite,45.0,5000,59,24.90,50").unwrap();
    path
}

#[test]
fn location_without_site_is_a_usage_error_not_a_silent_no_op() {
    let csv = location_csv("pair_a");
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--location", csv.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success(), "must not run with half the pair");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--site"),
        "the error must name the missing half: {stderr}"
    );
}

/// The pairing is deliberately ASYMMETRIC: the file without a selector was the silent
/// trap, but a standalone selector has a legitimate pre-existing job — it labels the PDF
/// dope card (rifle name / location line). Requiring the file here would have broken that.
#[test]
fn a_standalone_site_stays_legal_for_pdf_labeling() {
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--site", "TestSite"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "standalone --site is a display label and must keep working: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn an_unreadable_location_file_stops_the_run() {
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--location", "/nonexistent/bx_1425.csv", "--site", "X"])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "a typo'd path must not produce a default-atmosphere trajectory"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("bx_1425.csv"),
        "the error must name the path the user typed: {stderr}"
    );
}

#[test]
fn a_missing_site_row_stops_the_run_and_names_the_row() {
    let csv = location_csv("missing_row");
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--location", csv.to_str().unwrap(), "--site", "NoSuchSite"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("NoSuchSite"), "must name the row: {stderr}");
}

#[test]
fn profile_csv_without_profile_row_is_a_usage_error() {
    let csv = location_csv("profile_pair");
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--profile", csv.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--profile-row"));
}

/// And the same asymmetry for the profile pair.
#[test]
fn a_standalone_profile_row_stays_legal_for_pdf_labeling() {
    let out = Command::new(BIN)
        .args(base_args())
        .args(["--profile-row", "R700"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "standalone --profile-row is a display label: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// The control that MBA-1421's bad experiment lacked: the valid pair must still load AND must
/// actually change the physics versus no location at all.
#[test]
fn the_valid_pair_still_loads_and_reaches_the_atmosphere() {
    let csv = location_csv("happy");
    let with = Command::new(BIN)
        .args(base_args())
        .args(["--location", csv.to_str().unwrap(), "--site", "TestSite"])
        .output()
        .unwrap();
    assert!(with.status.success(), "{}", String::from_utf8_lossy(&with.stderr));
    let without = Command::new(BIN).args(base_args()).output().unwrap();

    let iv = |o: &std::process::Output| -> f64 {
        serde_json::from_slice::<serde_json::Value>(&o.stdout).unwrap()["impact_velocity"]
            .as_f64()
            .unwrap()
    };
    assert!(
        (iv(&with) - iv(&without)).abs() > 1.0,
        "the 5000 ft site must move the physics"
    );
}
