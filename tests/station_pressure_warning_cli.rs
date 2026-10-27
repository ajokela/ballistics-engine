//! The native CLI's station-pressure plausibility warning, run through the real binary.
//!
//! `trajectory` (shot day and zero day) and `zero` print a `warning:` line on stderr when the
//! absolute station pressure they fly is outside 0.85 to 1.08 times the ICAO standard pressure
//! at the altitude (`atmosphere::station_pressure_plausibility`). The message text itself is
//! unit-tested in `src/main.rs`; these tests pin the wiring at each call site and its gates:
//!
//!   - a weather report's altimeter setting, 30.50 inHg, typed as station pressure at Cusco
//!     (10860 ft) warns on the `trajectory` shot day, on `zero`, and on the `trajectory` zero
//!     day when it arrives through `--zero-pressure`;
//!   - the same pair declared as a QNH (`--pressure-type qnh`, `--zero-pressure-type qnh`)
//!     does not warn;
//!   - a `--density-altitude` run does not warn, because the density altitude supplies the
//!     pressure, and neither does a zero day that inherits that pressure;
//!   - the default pressure at 10860 ft does not warn, because the CLI flies the ICAO pressure
//!     there instead.

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

/// Present in every form of the warning, shot day and zero day, either side of the band.
const WARNING_MARKER: &str = "times the ICAO standard pressure at";

/// A .308 168 gr to 300 yd, imperial (the default unit system).
const TRAJECTORY: &[&str] = &[
    "trajectory",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--max-range",
    "300",
    "-o",
    "json",
];

/// The same rifle, solving a 100 yd zero.
const ZERO: &[&str] = &[
    "zero",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--target-distance",
    "100",
    "-o",
    "json",
];

/// Cusco, with a weather report's altimeter setting typed as station pressure.
const CUSCO_WITH_QNH_AS_STATION: &[&str] = &["--altitude", "10860", "--pressure", "30.50"];

/// Cusco with a real station reading on the shot day (1.02 times standard), zeroed at 100 yd
/// on a day whose `--zero-pressure` is the altimeter setting typed as station pressure.
const CUSCO_ZERO_DAY_WITH_QNH_AS_STATION: &[&str] = &[
    "--altitude",
    "10860",
    "--pressure",
    "20.33",
    "--auto-zero",
    "100",
    "--zero-pressure",
    "30.50",
];

/// Runs `base` + each slice of `extra`, asserts success, and returns stderr.
fn stderr_of(base: &[&str], extra: &[&[&str]]) -> String {
    let mut args = base.to_vec();
    for more in extra {
        args.extend_from_slice(more);
    }
    let output: Output = Command::new(BIN)
        .args(&args)
        .output()
        .expect("spawn ballistics");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "ballistics {args:?} failed: {stderr}"
    );
    stderr
}

/// The stderr lines that carry the warning.
fn warnings(stderr: &str) -> Vec<&str> {
    stderr
        .lines()
        .filter(|line| line.contains(WARNING_MARKER))
        .collect()
}

/// The one warning line `stderr` must carry.
fn single_warning(stderr: &str) -> &str {
    let lines = warnings(stderr);
    assert_eq!(
        lines.len(),
        1,
        "exactly one warning expected; stderr: {stderr}"
    );
    lines[0]
}

fn assert_no_warning(stderr: &str, case: &str) {
    assert!(
        warnings(stderr).is_empty(),
        "{case} must not warn; stderr: {stderr}"
    );
}

// ---- Must warn ------------------------------------------------------------------------------

#[test]
fn trajectory_shot_day_warns_for_a_qnh_typed_as_station_pressure() {
    let stderr = stderr_of(TRAJECTORY, &[CUSCO_WITH_QNH_AS_STATION]);
    let line = single_warning(&stderr);
    assert!(
        line.starts_with("warning: station pressure 30.50 inHg is 1.532 times"),
        "{line}"
    );
    assert!(line.contains("at 10860 ft (19.90 inHg)"), "{line}");
    assert!(line.contains("add --pressure-type qnh"), "{line}");
}

#[test]
fn zero_warns_for_a_qnh_typed_as_station_pressure() {
    let stderr = stderr_of(ZERO, &[CUSCO_WITH_QNH_AS_STATION]);
    let line = single_warning(&stderr);
    assert!(
        line.starts_with("warning: station pressure 30.50 inHg is 1.532 times"),
        "{line}"
    );
    assert!(line.contains("add --pressure-type qnh"), "{line}");
}

#[test]
fn trajectory_zero_day_warns_for_a_qnh_typed_as_zero_pressure() {
    let stderr = stderr_of(TRAJECTORY, &[CUSCO_ZERO_DAY_WITH_QNH_AS_STATION]);
    // The shot-day reading is real, so the zero day's is the only warning.
    let line = single_warning(&stderr);
    assert!(
        line.starts_with("warning: zero-day station pressure 30.50 inHg is 1.532 times"),
        "{line}"
    );
    assert!(line.contains("add --zero-pressure-type qnh"), "{line}");
    // No --zero-altitude was given, so the zero day flies at --altitude and the text says so
    // rather than pointing at a flag the user never typed.
    assert!(
        line.contains("The zero day uses --altitude because --zero-altitude was not given."),
        "{line}"
    );
    assert!(
        line.contains("check that --altitude is the station elevation"),
        "{line}"
    );
}

// ---- Must not warn --------------------------------------------------------------------------

#[test]
fn a_declared_qnh_is_not_warned() {
    let qnh: &[&str] = &["--pressure-type", "qnh"];
    assert_no_warning(
        &stderr_of(TRAJECTORY, &[CUSCO_WITH_QNH_AS_STATION, qnh]),
        "trajectory with --pressure-type qnh",
    );
    assert_no_warning(
        &stderr_of(ZERO, &[CUSCO_WITH_QNH_AS_STATION, qnh]),
        "zero with --pressure-type qnh",
    );
    assert_no_warning(
        &stderr_of(
            TRAJECTORY,
            &[
                CUSCO_ZERO_DAY_WITH_QNH_AS_STATION,
                &["--zero-pressure-type", "qnh"],
            ],
        ),
        "trajectory zero day with --zero-pressure-type qnh",
    );
}

#[test]
fn a_density_altitude_run_is_not_warned() {
    let da: &[&str] = &["--density-altitude", "10000"];
    assert_no_warning(
        &stderr_of(TRAJECTORY, &[CUSCO_WITH_QNH_AS_STATION, da]),
        "trajectory with --density-altitude",
    );
    assert_no_warning(
        &stderr_of(ZERO, &[CUSCO_WITH_QNH_AS_STATION, da]),
        "zero with --density-altitude",
    );
    // A zero day that brings only its own altitude inherits the pressure the density altitude
    // supplied, which is not an explicit pressure either, so it is not checked against
    // --zero-altitude.
    assert_no_warning(
        &stderr_of(
            TRAJECTORY,
            &[da, &["--auto-zero", "100", "--zero-altitude", "0"]],
        ),
        "trajectory zero day inheriting a --density-altitude pressure",
    );
}

#[test]
fn the_default_pressure_at_altitude_is_not_warned() {
    let altitude: &[&str] = &["--altitude", "10860"];
    assert_no_warning(
        &stderr_of(TRAJECTORY, &[altitude]),
        "trajectory with the default pressure",
    );
    assert_no_warning(
        &stderr_of(ZERO, &[altitude]),
        "zero with the default pressure",
    );
    // A zero day that brings only its own altitude inherits the omitted pressure with its
    // presence: still not a reading, still flown as the ICAO pressure there, so not warned.
    // (Since 0.49.0 the gate is presence, not the old value-based sentinel: a TYPED 29.92 inHg
    // at altitude is a reading and is warned; tests/station_presence_cli.rs pins that.)
    assert_no_warning(
        &stderr_of(
            TRAJECTORY,
            &[altitude, &["--auto-zero", "100", "--zero-altitude", "500"]],
        ),
        "trajectory zero day with the default pressure",
    );
}
