//! R3 (whole-branch review C2, VAL-1): every numeric cell `trajectory` reads from a
//! `--location` or `--profile` CSV row is read one way.
//!
//! - A cell is trimmed, so `" 70"` is 70, as it is on `--humidity`.
//! - An absent column, an empty cell or a whitespace-only cell is absent: the value falls back
//!   to the next source (a saved profile, then the default).
//! - Any other cell that is not a finite number stops the run with an error naming the column
//!   and quoting the cell. Before 0.50.0 it was read as absent, so `"65%"` in a HUMIDITY column
//!   silently flew 50 %, and `"5g"` in a TEMPERATURE column flew 59 °F.
//! - A HUMIDITY cell goes through the one humidity validator, and the stderr `Effective values`
//!   line echoes the humidity the run flies.

#![cfg(feature = "cli")]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const LOAD: &[&str] = &[
    "trajectory",
    "-v",
    "2700",
    "-b",
    "0.243",
    "-m",
    "168",
    "-d",
    "0.308",
    "--drag-model",
    "g7",
    "--max-range",
    "500",
    "-o",
    "json",
];

/// A private scratch directory, removed with everything in it when the test drops it.
struct PrivateDir(PathBuf);

impl Drop for PrivateDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn private_dir() -> PrivateDir {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-csv-cells-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    PrivateDir(dir)
}

/// A one-row CSV in `dir`: the header, then the row, named `S1` in its first column.
fn csv(dir: &PrivateDir, file: &str, header: &str, cells: &str) -> PathBuf {
    let path = dir.0.join(file);
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "NAME,{header}").unwrap();
    writeln!(f, "S1,{cells}").unwrap();
    path
}

fn run(args: &[&str]) -> Output {
    Command::new(BIN)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn with_location(location: &Path, extra: &[&str]) -> Vec<String> {
    LOAD.iter()
        .map(|s| s.to_string())
        .chain(["--location".into(), location.display().to_string()])
        .chain(["--site".into(), "S1".into()])
        .chain(extra.iter().map(|s| s.to_string()))
        .collect()
}

fn with_profile(profile: &Path, extra: &[&str]) -> Vec<String> {
    LOAD.iter()
        .map(|s| s.to_string())
        .chain(["--profile".into(), profile.display().to_string()])
        .chain(["--profile-row".into(), "S1".into()])
        .chain(extra.iter().map(|s| s.to_string()))
        .collect()
}

fn run_owned(args: &[String]) -> Output {
    run(&args.iter().map(String::as_str).collect::<Vec<_>>())
}

fn stdout_of(output: &Output, label: &str) -> String {
    assert!(
        output.status.success(),
        "{label} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn refusal_of(output: &Output, label: &str) -> String {
    assert!(
        !output.status.success(),
        "{label} must be refused, but it ran: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.status.code(), Some(1), "{label}: a refusal exits 1");
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Cells that are present but are no finite number, each of which a sloppy export or a typo
/// really produces.
const BAD_CELLS: [&str; 5] = ["abc", "65%", "8O", "5g", "nan"];

/// Every numeric `--location` column `trajectory` reads, with an ordinary value for it.
const LOCATION_COLUMNS: [(&str, &str); 11] = [
    ("TARGET_TEMP", "70"),
    ("TEMPERATURE", "70"),
    ("TEMP", "70"),
    ("PRESSURE", "29.5"),
    ("ALTITUDE", "1000"),
    ("ALT", "1000"),
    ("HUMIDITY", "70"),
    ("DA", "3000"),
    ("DENSITY_ALTITUDE", "3000"),
    ("WIND_DIR", "90"),
    ("WIND_DIRECTION", "90"),
];

/// Every numeric `--profile` column `trajectory` reads, with an ordinary value for it.
const PROFILE_COLUMNS: [(&str, &str); 14] = [
    ("VELOCITY", "2650"),
    ("MV", "2650"),
    ("MUZZLE_VELOCITY", "2650"),
    ("VELOCITY_ADJ", "10"),
    ("VEL_ADJ", "10"),
    ("BC", "0.25"),
    ("BC_ADJ", "0.98"),
    ("BC_ADJUSTMENT", "0.98"),
    ("ZERO_TEMP", "40"),
    ("ZERO_ALT", "500"),
    ("ZERO_RANGE", "200"),
    ("ZERO_DISTANCE", "200"),
    ("V_OFFSET_MIL", "0.2"),
    ("H_OFFSET_MIL", "0.1"),
];

#[test]
fn an_unreadable_location_cell_is_refused_naming_its_column_and_quoting_it() {
    let dir = private_dir();
    for (column, _) in LOCATION_COLUMNS {
        for bad in BAD_CELLS {
            let path = csv(&dir, "loc.csv", column, bad);
            let stderr = refusal_of(
                &run_owned(&with_location(&path, &[])),
                &format!("--location {column} = {bad:?}"),
            );
            assert!(
                stderr.contains(&format!("--location CSV column {column} = '{bad}'")),
                "the refusal names the column and quotes the cell: {stderr}"
            );
        }
    }
}

#[test]
fn an_unreadable_profile_cell_is_refused_naming_its_column_and_quoting_it() {
    let dir = private_dir();
    for (column, _) in PROFILE_COLUMNS {
        for bad in BAD_CELLS {
            let path = csv(&dir, "gun.csv", column, bad);
            let stderr = refusal_of(
                &run_owned(&with_profile(&path, &[])),
                &format!("--profile {column} = {bad:?}"),
            );
            assert!(
                stderr.contains(&format!("--profile CSV column {column} = '{bad}'")),
                "the refusal names the column and quotes the cell: {stderr}"
            );
        }
    }
}

#[test]
fn a_humidity_cell_that_is_no_number_states_the_humidity_rule() {
    let dir = private_dir();
    let path = csv(&dir, "loc.csv", "HUMIDITY", "65%");
    let stderr = refusal_of(&run_owned(&with_location(&path, &[])), "HUMIDITY 65%");
    assert!(
        stderr.contains("--location CSV column HUMIDITY = '65%': not a number")
            && stderr.contains("relative humidity is a percent from 0 to 100"),
        "{stderr}"
    );
}

#[test]
fn a_padded_cell_is_the_number_it_holds() {
    let dir = private_dir();
    for (column, value) in LOCATION_COLUMNS {
        let padded = csv(&dir, "padded.csv", column, &format!("  {value} "));
        let plain = csv(&dir, "plain.csv", column, value);
        assert_eq!(
            stdout_of(&run_owned(&with_location(&padded, &[])), column),
            stdout_of(&run_owned(&with_location(&plain, &[])), column),
            "--location {column} \"  {value} \" flies as {value}"
        );
    }
    for (column, value) in PROFILE_COLUMNS {
        let padded = csv(&dir, "padded.csv", column, &format!(" {value}  "));
        let plain = csv(&dir, "plain.csv", column, value);
        assert_eq!(
            stdout_of(&run_owned(&with_profile(&padded, &[])), column),
            stdout_of(&run_owned(&with_profile(&plain, &[])), column),
            "--profile {column} \" {value}  \" flies as {value}"
        );
    }
}

#[test]
fn an_empty_or_blank_cell_is_absent() {
    let dir = private_dir();
    let nothing = csv(&dir, "none.csv", "COMMENT", "x");
    let baseline = stdout_of(&run_owned(&with_location(&nothing, &[])), "no columns");
    for (column, _) in LOCATION_COLUMNS {
        for blank in ["", "   "] {
            let path = csv(&dir, "blank.csv", column, blank);
            assert_eq!(
                stdout_of(&run_owned(&with_location(&path, &[])), column),
                baseline,
                "--location {column} {blank:?} is absent"
            );
        }
    }
}

#[test]
fn the_effective_values_line_echoes_the_humidity_the_run_flies() {
    let dir = private_dir();
    for (cell, flown) in [
        ("65", "65.0"),
        ("102", "100.0"),
        ("-2", "0.0"),
        ("", "50.0"),
    ] {
        let path = csv(
            &dir,
            "loc.csv",
            "TEMPERATURE,HUMIDITY",
            &format!("95,{cell}"),
        );
        let output = run_owned(&with_location(&path, &[]));
        assert!(output.status.success(), "HUMIDITY {cell:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let line = stderr
            .lines()
            .find(|line| line.contains("altitude=") && line.contains("temp="))
            .unwrap_or_else(|| panic!("no location line in {stderr}"));
        assert!(
            line.ends_with(&format!("humidity={flown}%")),
            "HUMIDITY {cell:?} flies {flown} %: {line}"
        );
    }
}
