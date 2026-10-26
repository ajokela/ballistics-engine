//! MBA-1593: the CLI's dope cards state, in their header, the station air their rows were
//! solved in — on both CLI PDF surfaces, `trajectory -o pdf` and `adaptive-card -o pdf`. The
//! bridge's `card.pdf` is pinned the same way in `tests/card_pdf_bridge.rs`.
//!
//! Until 0.49.0 the header printed the air the flags STATED: the sea-level defaults when a
//! pressure or temperature was left out at altitude (which the solve replaces with ICAO air at
//! the site), and the QNH itself when one was declared (which the solve reduces to station
//! pressure). The expected headers here are derived from the air the solve flies,
//! `atmosphere::resolve_station_conditions_with_presence` with the run's presence (0.49.0: a
//! typed value flies as typed, an omitted one keeps the legacy rule), and the main ones are also
//! pinned as literals.

#![cfg(feature = "pdf")]

#[path = "support/pdf_text.rs"]
mod pdf_text;

use std::path::PathBuf;
use std::process::Command;

use ballistics_engine::atmosphere::{
    resolve_station_conditions_with_presence, PressureReferenceMode, StationAtmospherePresence,
};
use pdf_text::{contains_phrase, header_air_phrase, tokens};

/// The MBA-1593 report's metric load, as `trajectory` flags.
const TRAJECTORY_LOAD: &[&str] = &[
    "trajectory",
    "-v",
    "823",
    "-b",
    "0.475",
    "-m",
    "10.9",
    "-d",
    "7.82",
    "--max-range",
    "300",
    "--auto-zero",
    "100",
    "--humidity",
    "0",
    "--sample-trajectory",
];

/// A .308 load for the imperial runs, as `trajectory` flags.
const TRAJECTORY_LOAD_IMPERIAL: &[&str] = &[
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
    "--auto-zero",
    "100",
    "--humidity",
    "0",
    "--sample-trajectory",
];

/// The metric load as `adaptive-card` flags.
const ADAPTIVE_LOAD: &[&str] = &[
    "adaptive-card",
    "-v",
    "823",
    "-b",
    "0.475",
    "-m",
    "10.9",
    "-d",
    "7.82",
    "--zero-distance",
    "100",
    "--humidity",
    "0",
    "--start",
    "100",
    "--end",
    "600",
];

fn scratch_pdf(label: &str) -> PathBuf {
    let tag: String = label
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    std::env::temp_dir().join(format!("bx_header_air_{}_{tag}.pdf", std::process::id()))
}

/// Run the CLI with `units`, `load` and `air`, and return the PDF it wrote.
fn card(units: &str, load: &[&str], air: &[&str], label: &str) -> Vec<u8> {
    let out = scratch_pdf(label);
    let output = Command::new(env!("CARGO_BIN_EXE_ballistics"))
        .args(["--units", units])
        .args(load)
        .args(air)
        .args(["-o", "pdf", "--output-file", out.to_str().unwrap()])
        .output()
        .expect("run the ballistics CLI");
    assert!(
        output.status.success(),
        "{label}: the CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = std::fs::read(&out).unwrap_or_else(|_| panic!("{label}: no PDF written"));
    let _ = std::fs::remove_file(&out);
    bytes
}

fn assert_header_says(bytes: &[u8], label: &str, phrase: &str) {
    assert!(
        contains_phrase(bytes, label, phrase),
        "{label}: the header must read {phrase:?}: {:?}",
        tokens(bytes, label)
    );
}

/// The header for the air a run flies, from what it hands the solve (°C, hPa, m) and which of
/// the pair it typed.
fn flown_header(
    temperature_c: f64,
    pressure_hpa: f64,
    altitude_m: f64,
    mode: PressureReferenceMode,
    presence: StationAtmospherePresence,
) -> String {
    let (temperature_c, pressure_hpa) = resolve_station_conditions_with_presence(
        temperature_c,
        pressure_hpa,
        altitude_m,
        mode,
        presence,
    );
    header_air_phrase(temperature_c, pressure_hpa)
}

const NEITHER: StationAtmospherePresence = StationAtmospherePresence::NEITHER;
const BOTH: StationAtmospherePresence = StationAtmospherePresence::BOTH;
const PRESSURE_ONLY: StationAtmospherePresence = StationAtmospherePresence::new(false, true);

#[test]
fn trajectory_header_states_icao_air_when_none_is_supplied_at_altitude() {
    let want = flown_header(
        15.0,
        1013.25,
        1500.0,
        PressureReferenceMode::Absolute,
        NEITHER,
    );
    assert_eq!(want, "DA:4919 ft Pressure:24.97/846 Temp:41");
    let bytes = card(
        "metric",
        TRAJECTORY_LOAD,
        &["--altitude", "1500"],
        "omitted",
    );
    assert_header_says(&bytes, "trajectory, omitted air at 1500 m", &want);
}

/// Imperial runs hand the solve inHg and feet converted to hPa and metres; the header must
/// still land on the ICAO air at 5000 ft, not the 29.92 inHg / 59 °F defaults.
#[test]
fn trajectory_header_states_icao_air_when_none_is_supplied_at_altitude_imperial() {
    let want = flown_header(
        15.0,
        29.92 * 33.8639,
        5000.0 * 0.3048,
        PressureReferenceMode::Absolute,
        NEITHER,
    );
    assert_eq!(want, "DA:4998 ft Pressure:24.90/843 Temp:41");
    let bytes = card(
        "imperial",
        TRAJECTORY_LOAD_IMPERIAL,
        &["--altitude", "5000"],
        "omitted imperial",
    );
    assert_header_says(&bytes, "trajectory, omitted air at 5000 ft", &want);
}

#[test]
fn trajectory_header_states_the_temperature_flown_when_only_pressure_is_supplied() {
    let want = flown_header(
        15.0,
        845.6,
        1500.0,
        PressureReferenceMode::Absolute,
        PRESSURE_ONLY,
    );
    assert_eq!(want, "DA:4919 ft Pressure:24.97/846 Temp:41");
    let air = ["--altitude", "1500", "--pressure", "845.6"];
    let bytes = card("metric", TRAJECTORY_LOAD, &air, "station pressure");
    assert_header_says(&bytes, "trajectory, 845.6 hPa and no temperature", &want);
}

/// The ticket's QNH repro. Its expected 6088.6 ft is the DA of 845.6 hPa at 15 °C. Before
/// 0.49.0 a typed `--temperature 15` at altitude was read as the CLI's default-temperature
/// sentinel and flew the ICAO 5.25 °C, so this header said about 4919 ft. A typed value now
/// flies as typed (the station-sentinel rule), so the run flies 15 °C and the header says the
/// ticket's 6089 ft.
#[test]
fn trajectory_header_states_station_pressure_for_a_declared_qnh() {
    let want = flown_header(15.0, 1013.25, 1500.0, PressureReferenceMode::Qnh, BOTH);
    assert_eq!(want, "DA:6089 ft Pressure:24.97/846 Temp:59");
    let air = [
        "--altitude",
        "1500",
        "--temperature",
        "15",
        "--pressure",
        "1013.25",
        "--pressure-type",
        "qnh",
    ];
    let bytes = card("metric", TRAJECTORY_LOAD, &air, "qnh standard");
    assert_header_says(&bytes, "trajectory, QNH 1013.25 hPa at 1500 m", &want);

    // A non-standard QNH and temperature: the header prints the reduced station pressure,
    // not the 1030 hPa typed.
    let want = flown_header(20.0, 1030.0, 1500.0, PressureReferenceMode::Qnh, BOTH);
    assert_eq!(want, "DA:6126 ft Pressure:25.40/860 Temp:68");
    let air = [
        "--altitude",
        "1500",
        "--temperature",
        "20",
        "--pressure",
        "1030",
        "--pressure-type",
        "qnh",
    ];
    let bytes = card("metric", TRAJECTORY_LOAD, &air, "qnh 1030");
    assert_header_says(&bytes, "trajectory, QNH 1030 hPa at 1500 m", &want);
}

/// Where the solve flies exactly what it was handed, the header is what it always was, digit
/// for digit. 987.5 hPa and 30.245 inHg are readings whose unit round trip lands on the other
/// side of the printed rounding, so they fail if the header re-derives a figure it should keep.
/// A density-altitude run hands the solve already-resolved air and is unchanged too.
#[test]
fn trajectory_header_prints_sea_level_and_fully_stated_air_as_before() {
    for (air, want) in [
        (&[][..], "DA:0 ft Pressure:29.92/1013 Temp:59"),
        (
            &[
                "--altitude",
                "1500",
                "--pressure",
                "845.6",
                "--temperature",
                "5.25",
            ][..],
            "DA:4918 ft Pressure:24.97/846 Temp:41",
        ),
        (
            &["--pressure", "987.5", "--temperature", "20.5"][..],
            "DA:1539 ft Pressure:29.16/987 Temp:69",
        ),
        (
            &["--density-altitude", "1500"][..],
            "DA:4921 ft Pressure:24.97/845 Temp:41",
        ),
    ] {
        let label = format!("trajectory {air:?}");
        let bytes = card("metric", TRAJECTORY_LOAD, air, &label);
        assert_header_says(&bytes, &label, want);
    }

    let air = ["--pressure", "30.245", "--temperature", "59.5"];
    let bytes = card(
        "imperial",
        TRAJECTORY_LOAD_IMPERIAL,
        &air,
        "imperial 30.245 inHg",
    );
    assert_header_says(
        &bytes,
        "trajectory, 30.245 inHg at sea level",
        "DA:-336 ft Pressure:30.25/1024 Temp:60",
    );
}

#[test]
fn adaptive_card_header_states_icao_air_when_none_is_supplied_at_altitude() {
    let want = flown_header(
        15.0,
        1013.25,
        1500.0,
        PressureReferenceMode::Absolute,
        NEITHER,
    );
    assert_eq!(want, "DA:4919 ft Pressure:24.97/846 Temp:41");
    let bytes = card(
        "metric",
        ADAPTIVE_LOAD,
        &["--altitude", "1500"],
        "adaptive omitted",
    );
    assert_header_says(&bytes, "adaptive-card, omitted air at 1500 m", &want);
}

/// `adaptive-card` reduces a declared QNH before it solves, so its header already showed the
/// station pressure; it is the temperature it printed wrong (the 15 °C default, where the
/// solve flew the ICAO 5.25 °C), giving `DA:6089 ft` until 0.49.0.
#[test]
fn adaptive_card_header_states_the_temperature_flown_for_a_declared_qnh() {
    let want = flown_header(
        15.0,
        1013.25,
        1500.0,
        PressureReferenceMode::Qnh,
        PRESSURE_ONLY,
    );
    assert_eq!(want, "DA:4919 ft Pressure:24.97/846 Temp:41");
    let air = [
        "--altitude",
        "1500",
        "--pressure",
        "1013.25",
        "--pressure-type",
        "qnh",
    ];
    let bytes = card("metric", ADAPTIVE_LOAD, &air, "adaptive qnh");
    assert_header_says(&bytes, "adaptive-card, QNH 1013.25 hPa at 1500 m", &want);
}

#[test]
fn adaptive_card_header_prints_sea_level_and_fully_stated_air_as_before() {
    for (air, want) in [
        (&[][..], "DA:0 ft Pressure:29.92/1013 Temp:59"),
        (
            &[
                "--altitude",
                "1500",
                "--pressure",
                "845.6",
                "--temperature",
                "5.25",
            ][..],
            "DA:4918 ft Pressure:24.97/846 Temp:41",
        ),
        (
            &[
                "--altitude",
                "1500",
                "--pressure",
                "1030",
                "--pressure-type",
                "qnh",
                "--temperature",
                "20",
            ][..],
            "DA:6126 ft Pressure:25.40/860 Temp:68",
        ),
    ] {
        let label = format!("adaptive-card {air:?}");
        let bytes = card("metric", ADAPTIVE_LOAD, air, &label);
        assert_header_says(&bytes, &label, want);
    }
}
