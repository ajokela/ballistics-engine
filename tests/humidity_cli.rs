//! Humidity-consistency train, plan C4, at the CLI: one validator, one default, every flag.
//!
//! - Every `--humidity` is a PERCENT that is range-checked by the parser; 0.49.0 left
//!   `wind-card` and `compare` unchecked, so 150 flew clamped to 100 %.
//! - A humidity read from a `--location` CSV or a saved profile is validated too, and the
//!   error names where it came from. A refused stored humidity says to edit the profile's file,
//!   which it names, and never costs the profile its data: `profile save` reads the stored file
//!   for its carry-forward without the load-time field checks, so a re-save keeps every field it
//!   cannot express (gate G17). Nor does `profile save` overwrite a stored profile it cannot
//!   read, or carry forward a tracking CF that every later load would refuse.
//! - `monte-carlo` gains `--humidity` (0.49.0 always flew 50 %), and it reaches every mode:
//!   the fixed-count run, `--adaptive` and `--wez`.
//! - An omitted humidity is the one default: byte-for-byte the output of `--humidity 50`.
//! - Air that cannot exist is refused by the solver's one rule on every command that takes a
//!   station temperature.
//!
//! - The sensor band (S5b, owner decision 2026-11-01): a reading within 5 points of the range
//!   (101.7 in fog, -2 in desert air), from a flag, a `--location` CSV or a saved profile, is
//!   the saturated or dry air it reads, byte for byte the output of the end itself. Only a value
//!   beyond the band is refused.
//!
//! The clap defaults themselves are pinned to `DEFAULT_HUMIDITY_PERCENT` inside the binary's
//! own tests (`src/main.rs`, `humidity_flag_tests`).

#![cfg(feature = "cli")]

use std::io::Write;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn run(args: &[&str]) -> Output {
    Command::new(BIN)
        .args(args)
        .output()
        .expect("spawn ballistics")
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
        "{label} must be refused, but it ran:\n{}",
        String::from_utf8_lossy(&output.stdout)
    );
    format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn with<'a>(base: &[&'a str], extra: &[&'a str]) -> Vec<&'a str> {
    base.iter().chain(extra).copied().collect()
}

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
    "500",
    "-o",
    "json",
];

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
    "200",
];

const WIND_CARD: &[&str] = &[
    "wind-card",
    "-v",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
    "-b",
    "0.5",
    "--zero-distance",
    "100",
    "--end",
    "600",
];

const COMPARE: &[&str] = &[
    "compare",
    "--load",
    "175 SMK:g7:0.243:175:2650",
    "--load",
    "168 ELD-M:g7:0.523:168:2700",
    "--zero-distance",
    "100",
    "--end",
    "500",
    "--step",
    "100",
];

/// A seeded, small fixed-count Monte Carlo: reproducible run to run.
const MONTE_CARLO: &[&str] = &[
    "--units",
    "metric",
    "monte-carlo",
    "-v",
    "800",
    "-b",
    "0.5",
    "-m",
    "10",
    "-d",
    "7.62",
    "-n",
    "40",
    "--target-distance",
    "300",
    "--seed",
    "7",
    "-o",
    "json",
];

/// The same load under `--adaptive`, kept small.
const MONTE_CARLO_ADAPTIVE: &[&str] = &[
    "--units",
    "metric",
    "monte-carlo",
    "-v",
    "800",
    "-b",
    "0.5",
    "-m",
    "10",
    "-d",
    "7.62",
    "--target-distance",
    "300",
    "--adaptive",
    "--seed",
    "7",
    "--min-samples",
    "40",
    "--max-samples",
    "80",
    "--mc-batch-size",
    "40",
    "-o",
    "json",
];

/// A two-step WEZ sweep (deterministic: it seeds itself per range step).
const MONTE_CARLO_WEZ: &[&str] = &[
    "monte-carlo",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "-n",
    "40",
    "--wez",
    "--target-size",
    "6x6",
    "--wez-start",
    "200",
    "--wez-end",
    "300",
    "--wez-step",
    "100",
    "-o",
    "full",
];

#[test]
fn wind_card_and_compare_range_check_their_humidity_like_every_other_command() {
    for (label, base) in [("wind-card", WIND_CARD), ("compare", COMPARE)] {
        // Beyond the 5-point sensor band (S5b: -5 and 100.0001 are readings, flown below).
        for bad in ["150", "-5.0000001", "105.0000001", "nan"] {
            // The `=` form, so a negative value reaches the humidity parser as a value.
            let flag = format!("--humidity={bad}");
            let text = refusal_of(
                &run(&with(base, &[flag.as_str()])),
                &format!("{label} {flag}"),
            );
            assert!(text.contains("humidity"), "{label} {bad}: {text}");
        }
        for good in ["0", "100"] {
            stdout_of(
                &run(&with(base, &["--humidity", good])),
                &format!("{label} --humidity {good}"),
            );
        }
    }
}

#[test]
fn monte_carlo_refuses_an_out_of_range_humidity() {
    for base in [MONTE_CARLO, MONTE_CARLO_ADAPTIVE, MONTE_CARLO_WEZ] {
        for bad in ["150", "-5.0000001"] {
            let flag = format!("--humidity={bad}");
            let text = refusal_of(
                &run(&with(base, &[flag.as_str()])),
                &format!("monte-carlo {flag}"),
            );
            // The parser's own range refusal, not an unknown-flag error.
            assert!(
                text.contains("humidity") && text.contains("percent from 0 to 100"),
                "{text}"
            );
        }
    }
}

#[test]
fn monte_carlo_humidity_reaches_every_mode() {
    for (label, base) in [
        ("fixed-count", MONTE_CARLO),
        ("--adaptive", MONTE_CARLO_ADAPTIVE),
        ("--wez", MONTE_CARLO_WEZ),
    ] {
        let dry = stdout_of(&run(&with(base, &["--humidity", "0"])), label);
        let saturated = stdout_of(&run(&with(base, &["--humidity", "100"])), label);
        assert_ne!(
            dry, saturated,
            "{label}: --humidity 0 and 100 fly the same air, so the flag never reached it"
        );
    }
}

#[test]
fn an_omitted_humidity_is_byte_for_byte_humidity_50() {
    for (label, base) in [
        ("trajectory", TRAJECTORY),
        ("zero", ZERO),
        ("wind-card", WIND_CARD),
        ("compare", COMPARE),
        ("monte-carlo", MONTE_CARLO),
        ("monte-carlo --adaptive", MONTE_CARLO_ADAPTIVE),
        ("monte-carlo --wez", MONTE_CARLO_WEZ),
    ] {
        let omitted = stdout_of(&run(base), label);
        let typed = stdout_of(&run(&with(base, &["--humidity", "50"])), label);
        assert_eq!(omitted, typed, "{label}");
    }
}

// --- impossible air -----------------------------------------------------------------------------

const COME_UPS: &[&str] = &[
    "come-ups",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--zero-distance",
    "100",
    "--start",
    "200",
    "--end",
    "300",
    "--step",
    "100",
];

const ESTIMATE_BC: &[&str] = &[
    "estimate-bc",
    "-v",
    "2700",
    "-m",
    "168",
    "-d",
    "0.308",
    "--data",
    "100,2.5;200,10.2;300,23.5",
];

const OPTIMAL_ZERO: &[&str] = &[
    "optimal-zero",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--target",
    "200:10",
    "--target",
    "400:12",
    "--target",
    "600:18",
];

const MPBR: &[&str] = &[
    "mpbr",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--vital-zone",
    "8",
];

const MARK_TO_RANGE: &[&str] = &[
    "mark-to-range",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--mark",
    "2.0",
    "--zero-distance",
    "100",
];

/// 212 F at 29.92 inHg and 100 % is saturated air at the boiling point, which cannot exist; every
/// command refuses it with the solver's one rule (`StationAir::new`), naming the field, and flies
/// it at 50 %. The searches (`estimate-bc`, `optimal-zero`) used to drop every failed trial and
/// blame the data or the fitted zero instead. `monte-carlo` has no `--temperature`, so it always
/// flies 15 C, where every humidity exists. (`lead` reports a failed solve per row, by design,
/// and names the same refusal in each row.)
#[test]
fn every_command_with_a_temperature_refuses_impossible_air() {
    for (label, base) in [
        ("trajectory", TRAJECTORY),
        ("zero", ZERO),
        ("wind-card", WIND_CARD),
        ("compare", COMPARE),
        ("come-ups", COME_UPS),
        ("estimate-bc", ESTIMATE_BC),
        ("optimal-zero", OPTIMAL_ZERO),
        ("mpbr", MPBR),
        ("mark-to-range", MARK_TO_RANGE),
    ] {
        let text = refusal_of(
            &run(&with(base, &["--temperature", "212", "--humidity", "100"])),
            &format!("{label} at 212 F / 100 %"),
        );
        assert!(
            text.contains("atmosphere.humidity = 100: no such air exists"),
            "{label}: refused with StationAir's own reason, naming the field: {text}"
        );
        stdout_of(
            &run(&with(base, &["--temperature", "212", "--humidity", "50"])),
            &format!("{label} at 212 F / 50 %"),
        );
    }
}

// --- humidity read from a file ------------------------------------------------------------------

/// A private scratch directory, removed with everything in it when the test drops it.
struct PrivateDir(std::path::PathBuf);

impl PrivateDir {
    fn path(&self) -> &std::path::Path {
        &self.0
    }
}

impl Drop for PrivateDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn private_dir(tag: &str) -> PrivateDir {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-humidity-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    PrivateDir(dir)
}

/// A `--location` CSV with one site at the given HUMIDITY, and the directory that holds it.
fn location_csv(humidity: &str) -> (PrivateDir, std::path::PathBuf) {
    let dir = private_dir("csv");
    let path = dir.path().join("sites.csv");
    let mut f = std::fs::File::create(&path).expect("create csv");
    writeln!(f, "NAME,ALTITUDE,TEMPERATURE,PRESSURE,HUMIDITY").unwrap();
    writeln!(f, "TestSite,0,59,29.92,{humidity}").unwrap();
    (dir, path)
}

#[test]
fn a_location_csv_humidity_outside_zero_to_one_hundred_is_refused_naming_its_source() {
    // Beyond the 5-point sensor band (S5b: -5 and 100.5 are readings, flown in the band test).
    for bad in ["150", "-5.0000001", "105.0000001"] {
        let (_dir, csv) = location_csv(bad);
        let text = refusal_of(
            &run(&with(
                TRAJECTORY,
                &["--location", csv.to_str().unwrap(), "--site", "TestSite"],
            )),
            &format!("--location HUMIDITY {bad}"),
        );
        assert!(
            text.contains("HUMIDITY") && text.contains("--location"),
            "the error names the CSV column it came from: {text}"
        );
    }
    // A valid CSV humidity flies, and is the same air as the flag.
    let (_dir, csv) = location_csv("80");
    let from_csv = stdout_of(
        &run(&with(
            TRAJECTORY,
            &["--location", csv.to_str().unwrap(), "--site", "TestSite"],
        )),
        "--location HUMIDITY 80",
    );
    let from_flag = stdout_of(
        &run(&with(
            TRAJECTORY,
            &[
                "--location",
                csv.to_str().unwrap(),
                "--site",
                "TestSite",
                "--humidity",
                "80",
            ],
        )),
        "--humidity 80",
    );
    assert_eq!(from_csv, from_flag);
}

/// `$HOME/.ballistics/profiles/humid.json` under a private home.
fn humid_profile_path(home: &std::path::Path) -> std::path::PathBuf {
    home.join(".ballistics").join("profiles").join("humid.json")
}

fn home_with_profile(humidity: f64) -> PrivateDir {
    let home = private_dir("home");
    let path = humid_profile_path(home.path());
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(
        &path,
        serde_json::json!({
            "name": "humid", "velocity": 2700.0, "bc": 0.243, "mass": 175.0,
            "diameter": 0.308, "drag_model": "G7", "humidity": humidity
        })
        .to_string(),
    )
    .unwrap();
    home
}

/// Runs the binary with `HOME` pointed at `home`.
fn run_at_home(home: &std::path::Path, args: &[&str]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

const SAVED_PROFILE_TRAJECTORY: &[&str] = &[
    "trajectory",
    "--saved-profile",
    "humid",
    "--max-range",
    "500",
    "-o",
    "json",
];

#[test]
fn a_saved_profile_humidity_outside_zero_to_one_hundred_is_refused_naming_its_source() {
    for bad in [150.0, -5.000_000_1] {
        let home = home_with_profile(bad);
        let output = run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY);
        let text = refusal_of(&output, &format!("saved profile humidity {bad}"));
        assert!(
            text.contains("profile") && text.contains("humidity"),
            "the error names the profile field: {text}"
        );
    }
    let home = home_with_profile(100.0);
    let output = run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY);
    stdout_of(&output, "a saved profile at 100 %");
}

// --- a refused stored value never costs a profile its data (gate G17) ---------------------------

/// `profile save humid` with the five mandatory fields; the caller adds `--humidity`.
const SAVE_HUMID: &[&str] = &[
    "profile",
    "save",
    "humid",
    "--velocity",
    "2700",
    "--bc",
    "0.243",
    "--mass",
    "175",
    "--diameter",
    "0.308",
];

/// Every field `profile save` cannot express from its flags alone, so carries forward from the
/// stored file: the DSF table, the zero POI and lateral sight offset, the tracking CFs, the zero
/// sets, the reticle, the click graduations and the twelve turret/hold fields. Distinct values
/// that pass every check (the shapes of `src/profile.rs`'s wire-shape fixture).
fn carried_fields() -> serde_json::Value {
    serde_json::json!({
        "dsf_points": [{"mach": 0.9, "dsf": 1.04}, {"mach": 1.1, "dsf": 0.98}],
        "zero_poi_up_m": 0.01,
        "zero_poi_right_m": -0.02,
        "sight_offset_lateral_m": 0.03,
        "elevation_cf": 1.05,
        "windage_cf": 0.97,
        "zero_sets": [
            {"name": "suppressed", "zero_distance": 200.0, "poi_up_mil": -0.3,
             "poi_right_mil": 0.1, "notes": "suppressed load"}
        ],
        "reticle": {
            "name": "mil-grid 0.5/10",
            "focal_plane": "ffp",
            "reference_magnification": 10.0,
            "marks": [
                {"down_mil": 0.0, "right_mil": 0.0, "kind": "center"},
                {"down_mil": 1.0, "right_mil": 0.0, "kind": "hash", "label": "1.0"}
            ]
        },
        "elevation_click": "0.1mil",
        "windage_click": "0.25moa",
        "clicks_per_revolution": 100,
        "zero_stop": true,
        "elevation_travel_up_mil": 26.0,
        "elevation_travel_down_mil": 4.0,
        "windage_travel_left_mil": 12.0,
        "windage_travel_right_mil": 12.0,
        "turret_elevation_dialed_mil": 5.4,
        "turret_windage_dialed_mil": -0.2,
        "hold_bound_up_mil": 5.0,
        "hold_bound_down_mil": 10.0,
        "hold_bound_left_mil": 6.0,
        "hold_bound_right_mil": 6.0
    })
}

fn stored_profile(home: &std::path::Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(humid_profile_path(home)).unwrap()).unwrap()
}

fn store_profile(home: &std::path::Path, profile: &serde_json::Value) {
    std::fs::write(
        humid_profile_path(home),
        serde_json::to_string_pretty(profile).unwrap(),
    )
    .unwrap();
}

/// A private home whose `humid` profile holds every carried field, as `profile save` itself
/// wrote it (a valid re-save at 60 % carried them forward), and that stored profile.
fn home_with_carried_fields() -> (PrivateDir, serde_json::Value) {
    let home = private_dir("carry");
    stdout_of(
        &run_at_home(home.path(), &with(SAVE_HUMID, &["--humidity", "50"])),
        "the first save",
    );
    let mut profile = stored_profile(home.path());
    for (key, value) in carried_fields().as_object().unwrap() {
        profile[key] = value.clone();
    }
    store_profile(home.path(), &profile);
    stdout_of(
        &run_at_home(home.path(), &with(SAVE_HUMID, &["--humidity", "60"])),
        "a valid re-save at 60 %",
    );
    let baseline = stored_profile(home.path());
    assert_eq!(baseline["humidity"], serde_json::json!(60.0));
    for key in carried_fields().as_object().unwrap().keys() {
        assert!(
            !baseline[key].is_null(),
            "a valid re-save carries {key} forward: {baseline}"
        );
    }
    stdout_of(
        &run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY),
        "the stored profile flies",
    );
    (home, baseline)
}

#[test]
fn re_saving_a_profile_whose_stored_humidity_is_refused_keeps_every_carried_field() {
    for bad in [150.0, -5.000_000_1] {
        let (home, baseline) = home_with_carried_fields();
        let mut profile = baseline.clone();
        profile["humidity"] = serde_json::json!(bad);
        store_profile(home.path(), &profile);
        refusal_of(
            &run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY),
            &format!("a stored humidity of {bad}"),
        );
        // The user re-saves at a valid humidity: the stored value the load refused is replaced
        // by the flag's, and nothing else the profile stored is lost.
        stdout_of(
            &run_at_home(home.path(), &with(SAVE_HUMID, &["--humidity", "50"])),
            &format!("the re-save over a stored {bad}"),
        );
        let resaved = stored_profile(home.path());
        assert_eq!(resaved["humidity"], serde_json::json!(50.0));
        for key in carried_fields().as_object().unwrap().keys() {
            assert_eq!(
                resaved[key], baseline[key],
                "RH {bad}: the re-save keeps the stored {key}"
            );
        }
        stdout_of(
            &run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY),
            "the re-saved profile flies",
        );
    }
}

#[test]
fn a_refused_stored_humidity_says_to_edit_the_profile_file_it_names() {
    let home = home_with_profile(150.0);
    let text = refusal_of(
        &run_at_home(home.path(), SAVED_PROFILE_TRAJECTORY),
        "a stored humidity of 150",
    );
    let path = humid_profile_path(home.path());
    let advice = format!(
        "profile 'humid' field humidity = 150: relative humidity must be a percent from 0 to 100; \
         edit the humidity stored in the profile to a percent from 0 to 100 (the profile is \
         stored at {})",
        path.display()
    );
    assert!(
        text.contains(&advice),
        "the refusal names the field, the rule and the file to edit:\n{text}"
    );
}

#[test]
fn profile_save_never_overwrites_a_stored_profile_it_cannot_read() {
    let (home, _) = home_with_carried_fields();
    let path = humid_profile_path(home.path());
    let mut damaged = std::fs::read_to_string(&path).unwrap();
    damaged.truncate(damaged.len() / 2);
    std::fs::write(&path, &damaged).unwrap();
    let text = refusal_of(
        &run_at_home(home.path(), &with(SAVE_HUMID, &["--humidity", "50"])),
        "a re-save over an unreadable profile",
    );
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        damaged,
        "the stored file is untouched"
    );
    assert!(
        text.contains(&path.display().to_string()) && text.contains("profile delete humid"),
        "the refusal names the file and how to start over: {text}"
    );
}

#[test]
fn profile_save_never_carries_forward_a_tracking_cf_every_load_would_refuse() {
    for (key, bad) in [("elevation_cf", 1.6), ("windage_cf", 0.4)] {
        let (home, mut profile) = home_with_carried_fields();
        profile[key] = serde_json::json!(bad);
        store_profile(home.path(), &profile);
        let before = std::fs::read_to_string(humid_profile_path(home.path())).unwrap();
        let text = refusal_of(
            &run_at_home(home.path(), &with(SAVE_HUMID, &["--humidity", "50"])),
            &format!("a re-save carrying {key} {bad}"),
        );
        assert_eq!(
            std::fs::read_to_string(humid_profile_path(home.path())).unwrap(),
            before,
            "{key} {bad}: the stored file is untouched"
        );
        assert!(
            text.contains(key)
                && text.contains(&humid_profile_path(home.path()).display().to_string()),
            "the refusal names the field and the file: {text}"
        );
    }
}

// --- the sensor band (S5b) ----------------------------------------------------------------------

/// Readings within 5 points of the range, each with the end it reads as.
const SENSOR_READINGS: [(&str, &str); 6] = [
    ("100.0000001", "100"),
    ("101.7", "100"),
    ("105", "100"),
    ("-0", "0"),
    ("-2", "0"),
    ("-5", "0"),
];

/// `trajectory --auto-zero` with a zero-day humidity.
const AUTO_ZERO: &[&str] = &[
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
    "500",
    "--auto-zero",
    "200",
    "-o",
    "json",
];

#[test]
fn every_command_flies_a_sensor_reading_as_the_end_it_reads_as() {
    for (label, base, flag) in [
        ("trajectory", TRAJECTORY, "--humidity"),
        ("trajectory --auto-zero", AUTO_ZERO, "--zero-humidity"),
        ("zero", ZERO, "--humidity"),
        ("wind-card", WIND_CARD, "--humidity"),
        ("compare", COMPARE, "--humidity"),
        ("come-ups", COME_UPS, "--humidity"),
        ("estimate-bc", ESTIMATE_BC, "--humidity"),
        ("mpbr", MPBR, "--humidity"),
        ("mark-to-range", MARK_TO_RANGE, "--humidity"),
        ("monte-carlo", MONTE_CARLO, "--humidity"),
        ("monte-carlo --adaptive", MONTE_CARLO_ADAPTIVE, "--humidity"),
        ("monte-carlo --wez", MONTE_CARLO_WEZ, "--humidity"),
    ] {
        for (reading, end) in SENSOR_READINGS {
            let at = |value: &str| {
                let typed = format!("{flag}={value}");
                stdout_of(
                    &run(&with(base, &[typed.as_str()])),
                    &format!("{label} {typed}"),
                )
            };
            assert_eq!(at(reading), at(end), "{label} {flag}={reading}");
        }
        // A negative reading needs no `=`: the flag takes it as its value.
        let spaced = stdout_of(
            &run(&with(base, &[flag, "-2"])),
            &format!("{label} {flag} -2"),
        );
        let end = format!("{flag}=0");
        assert_eq!(
            spaced,
            stdout_of(
                &run(&with(base, &[end.as_str()])),
                &format!("{label} {end}")
            ),
            "{label} {flag} -2"
        );
    }
}

#[test]
fn a_location_csv_or_saved_profile_reading_flies_as_the_end_it_reads_as() {
    for (reading, end) in SENSOR_READINGS {
        let (_dir, csv) = location_csv(reading);
        let from_csv = stdout_of(
            &run(&with(
                TRAJECTORY,
                &["--location", csv.to_str().unwrap(), "--site", "TestSite"],
            )),
            &format!("--location HUMIDITY {reading}"),
        );
        let (_end_dir, end_csv) = location_csv(end);
        let at_end = stdout_of(
            &run(&with(
                TRAJECTORY,
                &[
                    "--location",
                    end_csv.to_str().unwrap(),
                    "--site",
                    "TestSite",
                ],
            )),
            &format!("--location HUMIDITY {end}"),
        );
        assert_eq!(from_csv, at_end, "--location HUMIDITY {reading}");

        let reading_home = home_with_profile(reading.parse().unwrap());
        let end_home = home_with_profile(end.parse().unwrap());
        assert_eq!(
            stdout_of(
                &run_at_home(reading_home.path(), SAVED_PROFILE_TRAJECTORY),
                &format!("a saved profile at {reading}")
            ),
            stdout_of(
                &run_at_home(end_home.path(), SAVED_PROFILE_TRAJECTORY),
                &format!("a saved profile at {end}")
            ),
            "saved profile humidity {reading}"
        );
    }
}

#[test]
fn profile_save_stores_a_sensor_reading_as_the_end_it_reads_as() {
    for (reading, end) in SENSOR_READINGS {
        let home = private_dir("save");
        let flag = format!("--humidity={reading}");
        stdout_of(
            &run_at_home(home.path(), &with(SAVE_HUMID, &[flag.as_str()])),
            &format!("profile save {flag}"),
        );
        assert_eq!(
            stored_profile(home.path())["humidity"]
                .as_f64()
                .unwrap()
                .to_bits(),
            end.parse::<f64>().unwrap().to_bits(),
            "profile save {flag}"
        );
    }
}
