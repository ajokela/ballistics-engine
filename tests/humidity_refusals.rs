//! Humidity-consistency train, review step R2: one way to refuse a humidity or air that cannot
//! exist, and every refusal readable.
//!
//! - **Display, never Debug.** Every error that reaches a CLI user prints its message: no Rust
//!   struct dump (`Error: BallisticsError { message: ... }`), no quoted and escaped string
//!   (`Error: "..."`), no wire code through `{:?}` (`SolveFailed`).
//! - **The cause first.** Air that cannot exist is refused up front, before any per-row or
//!   per-candidate work, so the first line names it in the surface's own words (`--humidity`,
//!   `--zero-humidity`, the `--location` CSV column, the saved profile's field), and a command
//!   whose every row would fail exits non-zero instead of printing a table of dashes.
//! - **The band, stated truthfully.** A humidity refusal states the range it enforces: 0 to 100
//!   with a sensor reading up to [`HUMIDITY_SENSOR_TOLERANCE_PERCENT`] outside read as the end
//!   (and the same band in fraction units on the SI wire). The text is derived here from the
//!   constant, so a change to the band and a stale message cannot both pass.
//! - **One machine-readable envelope.** On every bridge command an out-of-band humidity and
//!   impossible air carry `error.details.error.code == "invalid_value"` at the path of the
//!   field to blame in that command's own request.
//! - **Profiles.** `profile.validate`, `.a7p` export and import, and the CLI's `profile save`
//!   and `profile import` apply the impossible-air rule a flight does; `profile show` displays
//!   a profile whose stored humidity it would refuse to fly, with the refusal as a warning.

#![cfg(feature = "cli")]

use ballistics_engine::HUMIDITY_SENSOR_TOLERANCE_PERCENT;
use serde_json::{json, Value};
use std::io::Write;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

// --- the words every humidity refusal states ----------------------------------------------------

/// The accepted range of a PERCENT, as every percent refusal states it, from the constant.
fn percent_band() -> String {
    let t = HUMIDITY_SENSOR_TOLERANCE_PERCENT;
    format!(
        "a percent from 0 to 100 (a sensor reading from {} to {} is read as 0 or 100)",
        -t,
        100.0 + t
    )
}

/// The accepted range of a FRACTION (the SI wire's `relative_humidity`), from the constant.
fn fraction_band() -> String {
    let t = HUMIDITY_SENSOR_TOLERANCE_PERCENT / 100.0;
    format!(
        "a fraction from 0 to 1 (0.5 = 50 %; a sensor reading from {} to {} is read as 0 or 1)",
        -t,
        1.0 + t
    )
}

/// Any Rust Debug rendering a user could be shown in place of a message.
fn debug_forms(text: &str) -> Vec<String> {
    let mut found = Vec::new();
    for needle in [
        "{ message:",
        "BallisticsError",
        "AtmosphereError {",
        "TruingPlanErrorV1 {",
        "RejectedTruingCandidateV1",
        "SolveFailed",
        "InvalidValue",
        "HumidityPercentOutOfRange",
        "Error: \"",
    ] {
        if text.contains(needle) {
            found.push(needle.to_string());
        }
    }
    for line in text.lines() {
        if line.starts_with("Error: \"") || line.starts_with("Error: Error(") {
            found.push(line.to_string());
        }
    }
    found
}

// --- CLI plumbing -------------------------------------------------------------------------------

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
        "bx-refusals-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    PrivateDir(dir)
}

/// Runs the binary with `HOME` at a fresh private directory, so no saved profile or ToS state of
/// the machine's user can reach the run.
fn run_in(home: &std::path::Path, args: &[&str]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn run(args: &[&str]) -> Output {
    let home = private_dir("home");
    run_in(home.path(), args)
}

fn with<'a>(base: &[&'a str], extra: &[&'a str]) -> Vec<&'a str> {
    base.iter().chain(extra).copied().collect()
}

/// The refusal a CLI run printed: it must have exited 1, printed no Debug form anywhere, and
/// put the cause on stderr's first line after any note of which inputs were read: `Error: `
/// and then `cause`.
fn assert_cli_refusal(label: &str, output: &Output, cause: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(
        output.status.code(),
        Some(1),
        "{label}: refused with exit status 1\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    // A note that says which inputs were read (`Loaded location 'X' from loc.csv`) may come
    // first; no progress line of the work itself may, and no other line.
    let mut lines = stderr.lines();
    let mut first = lines.next().unwrap_or_default();
    while first.starts_with("Loaded ") {
        assert!(
            !first.contains('"'),
            "{label}: a path is printed as itself, never quoted through Debug: {first}"
        );
        first = lines.next().unwrap_or_default();
    }
    assert!(
        first.starts_with(&format!("Error: {cause}")),
        "{label}: stderr's first line leads with the cause\n  expected prefix: Error: {cause}\n  \
         got: {first}\nstderr:\n{stderr}"
    );
    let debug = debug_forms(&format!("{stdout}{stderr}"));
    assert!(
        debug.is_empty(),
        "{label}: printed a Debug form {debug:?}\nstderr:\n{stderr}"
    );
}

const LOAD: &[&str] = &["-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82"];
const BULLET: &[&str] = &[
    "--bc",
    "0.475",
    "--drag-model",
    "g1",
    "--velocity",
    "823",
    "-m",
    "10.9",
    "-d",
    "7.82",
];

/// Every CLI command that takes a station `--temperature` and `--humidity`, metric, at air that
/// flies (15 C, 50 %) when the impossible pair is not appended.
fn cli_commands(dir: &std::path::Path) -> Vec<(&'static str, Vec<String>)> {
    let scenarios = scenario_file(dir);
    let reticle = reticle_file(dir);
    let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    let cat = |a: &[&str], b: &[&str]| s(&with(a, b));
    vec![
        (
            "trajectory",
            cat(&["trajectory"], &with(LOAD, &["--max-range", "500"])),
        ),
        (
            "zero",
            cat(&["zero"], &with(BULLET, &["--target-distance", "300"])),
        ),
        (
            "mpbr",
            cat(&["mpbr"], &with(BULLET, &["--vital-zone", "20"])),
        ),
        (
            "come-ups",
            cat(
                &["come-ups"],
                &with(BULLET, &["--zero-distance", "300", "--end", "600"]),
            ),
        ),
        (
            "lead",
            cat(&["lead"], &with(BULLET, &["--target-speed", "4"])),
        ),
        (
            "wind-card",
            cat(&["wind-card"], &with(BULLET, &["--zero-distance", "300"])),
        ),
        (
            "range-table",
            cat(&["range-table"], &with(BULLET, &["--zero-distance", "300"])),
        ),
        (
            "true-velocity",
            s(&[
                "true-velocity",
                "--measured-drop",
                "8",
                "--range",
                "500",
                "--bc",
                "0.475",
                "--mass",
                "10.9",
                "--diameter",
                "7.82",
                "--offline",
                "--zero-distance",
                "300",
            ]),
        ),
        (
            "true-wind",
            cat(
                &["true-wind", "--miss", "500:15"],
                &with(BULLET, &["--twist-rate", "254", "--zero-distance", "300"]),
            ),
        ),
        (
            "plan-truing",
            s(&[
                "plan-truing",
                "--measurement-resolution",
                "0.25",
                "--candidate-ranges",
                "300,500,700",
                "--velocity",
                "823",
                "--bc",
                "0.475",
                "--mass",
                "10.9",
                "--diameter",
                "7.82",
                "--zero-distance",
                "300",
            ]),
        ),
        (
            "compare",
            s(&[
                "compare",
                "--load",
                "A:g1:0.475:10.9:823:7.82",
                "--load",
                "B:g1:0.500:11.66:838:7.82",
                "--zero-distance",
                "300",
            ]),
        ),
        (
            "estimate-bc",
            s(&[
                "estimate-bc",
                "--velocity",
                "823",
                "--mass",
                "10.9",
                "--diameter",
                "7.82",
                "--data",
                "300,0.85;500,2.7",
                "--zero-range",
                "100",
            ]),
        ),
        (
            "mark-to-range",
            cat(
                &["mark-to-range", "--mark", "2.0", "--mark", "5.0"],
                &with(LOAD, &["--zero-distance", "300"]),
            ),
        ),
        (
            "bdc-match",
            cat(
                &[
                    "bdc-match",
                    "--mark-range",
                    "2:300",
                    "--mark-range",
                    "4:500",
                    "--mark-range",
                    "6:700",
                ],
                &with(LOAD, &["--reference-mag", "12", "--zero-distance", "300"]),
            ),
        ),
        (
            "optimal-zero",
            cat(
                &[
                    "optimal-zero",
                    "--target",
                    "200:10",
                    "--target",
                    "400:12",
                    "--target",
                    "600:18",
                ],
                LOAD,
            ),
        ),
        (
            "adaptive-card",
            cat(
                &["adaptive-card"],
                &with(
                    LOAD,
                    &["--zero-distance", "300", "--start", "100", "--end", "600"],
                ),
            ),
        ),
        (
            "hold-corridor",
            cat(
                &[
                    "hold-corridor",
                    "--scenarios",
                    scenarios.to_str().unwrap(),
                    "--ranges",
                    "200,400,600",
                    "--target",
                    "rect:12x18",
                ],
                &with(LOAD, &["--zero-distance", "300"]),
            ),
        ),
        (
            "reticle hold",
            cat(
                &[
                    "reticle",
                    "hold",
                    "--reticle-json",
                    reticle.to_str().unwrap(),
                    "--range",
                    "300",
                    "--drag-model",
                    "g1",
                ],
                &with(LOAD, &["--zero-distance", "100"]),
            ),
        ),
    ]
}

/// A plain mil cross, as `reticle generate` writes it, for `reticle hold --range`.
fn reticle_file(dir: &std::path::Path) -> std::path::PathBuf {
    let output = Command::new(BIN)
        .args(["reticle", "generate", "mil-grid", "-o", "json"])
        .output()
        .expect("spawn ballistics");
    assert!(output.status.success(), "reticle generate");
    let path = dir.join("reticle.json");
    std::fs::write(&path, &output.stdout).unwrap();
    path
}

fn scenario_file(dir: &std::path::Path) -> std::path::PathBuf {
    let path = dir.join("scenarios.json");
    std::fs::write(
        &path,
        r#"{"version":1,"nominal":"low","scenarios":[{"name":"low","segments":["4:90:1000"]},{"name":"high","segments":["14:90:1000"]}]}"#,
    )
    .unwrap();
    path
}

/// The impossible-air refusal every metric `--humidity` command gives at 100 C and 100 % under
/// the default 1013.25 hPa: the flag, the value, the air it was refused in, and the reason.
const CLI_IMPOSSIBLE: &str = "--humidity = 100 at 100 °C and 1013.25 hPa: no such air exists";

#[test]
fn every_cli_command_refuses_impossible_air_up_front_with_the_cause_first() {
    let dir = private_dir("scen");
    for (label, base) in cli_commands(dir.path()) {
        let base: Vec<&str> = base.iter().map(String::as_str).collect();
        let metric = with(&["-u", "metric"], &base);
        let fine = run(&with(&metric, &["--temperature", "15", "--humidity", "50"]));
        assert!(
            fine.status.success(),
            "{label}: the fixture flies at 15 C / 50 %: {}",
            String::from_utf8_lossy(&fine.stderr)
        );
        let refused = run(&with(
            &metric,
            &["--temperature", "100", "--humidity", "100"],
        ));
        assert_cli_refusal(
            &format!("{label} at 100 C / 100 %"),
            &refused,
            CLI_IMPOSSIBLE,
        );
        // Half the humidity in the same heat is air that exists, and flies.
        let half = run(&with(
            &metric,
            &["--temperature", "100", "--humidity", "50"],
        ));
        assert!(
            half.status.success(),
            "{label}: 100 C at 50 % exists: {}",
            String::from_utf8_lossy(&half.stderr)
        );
    }
}

/// A station above the troposphere because no pressure was given is the altitude's fault: every
/// command that takes a station `--altitude` says to check the altitude's unit, never the unit of
/// a `--pressure` nobody typed (gate R2), whether it resolves the station itself (`trajectory`)
/// or hands the presence to the engine (`zero`, the cards, the truing commands). 20 km of
/// standard atmosphere is about 55 hPa, which saturated air at 40 C exceeds.
#[test]
fn every_cli_command_blames_an_omitted_pressure_above_the_troposphere_on_the_altitude() {
    let dir = private_dir("alt");
    for (label, base) in cli_commands(dir.path()) {
        let base: Vec<&str> = base.iter().map(String::as_str).collect();
        let metric = with(&["-u", "metric"], &base);
        let output = run(&with(
            &metric,
            &[
                "--altitude",
                "20000",
                "--temperature",
                "40",
                "--humidity",
                "100",
            ],
        ));
        assert_cli_refusal(
            &format!("{label} at 20 km"),
            &output,
            "--humidity = 100 at 40 °C and ",
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("no --pressure was given")
                && stderr.contains("check the altitude's unit")
                && !stderr.contains("check the pressure's unit"),
            "{label}: blames the altitude\n{stderr}"
        );
    }
    // A QNH is a sea-level figure: one that is itself plausible, reduced to a station 20 km up,
    // is the altitude's fault too.
    let output = run(&with(
        &["-u", "metric", "trajectory"],
        &with(
            LOAD,
            &[
                "--max-range",
                "500",
                "--altitude",
                "20000",
                "--pressure",
                "1013.25",
                "--pressure-type",
                "qnh",
                "--temperature",
                "40",
                "--humidity",
                "100",
            ],
        ),
    ));
    assert_cli_refusal(
        "trajectory QNH at 20 km",
        &output,
        "--humidity = 100 at 40 °C and ",
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--pressure is a QNH") && stderr.contains("check the altitude's unit"),
        "{stderr}"
    );
    // The zero day resolves its own station: its own altitude is to blame. No `--zero-pressure`
    // was given, so it inherits the shot day's, which was not given either.
    let output = run(&with(
        &["-u", "metric", "trajectory"],
        &with(
            LOAD,
            &[
                "--max-range",
                "500",
                "--auto-zero",
                "100",
                "--zero-altitude",
                "20000",
                "--zero-temperature",
                "40",
                "--zero-humidity",
                "100",
            ],
        ),
    ));
    assert_cli_refusal(
        "trajectory zero day at 20 km",
        &output,
        "--zero-humidity = 100 at 40 °C and ",
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("no --pressure was given") && stderr.contains("check the altitude's unit"),
        "{stderr}"
    );
}

/// The zero day names each reading where it came from: a humidity the zero day inherits from
/// `--humidity` is `--humidity`, even when another zero-day flag is set (gate R2).
#[test]
fn the_zero_day_names_an_inherited_humidity_by_the_flag_it_came_from() {
    let output = run(&with(
        &["-u", "metric", "trajectory"],
        &with(
            LOAD,
            &[
                "--max-range",
                "500",
                "--auto-zero",
                "100",
                "--zero-temperature",
                "100",
                "--humidity",
                "100",
            ],
        ),
    ));
    assert_cli_refusal(
        "trajectory --humidity 100 at --zero-temperature 100",
        &output,
        "--humidity = 100 at 100 °C and 1013.25 hPa: no such air exists",
    );
}

#[test]
fn the_zero_day_air_is_refused_under_its_own_flag() {
    let output = run(&with(
        &["-u", "metric", "trajectory"],
        &with(
            LOAD,
            &[
                "--max-range",
                "500",
                "--auto-zero",
                "100",
                "--zero-temperature",
                "100",
                "--zero-humidity",
                "100",
            ],
        ),
    ));
    assert_cli_refusal(
        "trajectory --zero-humidity 100 at --zero-temperature 100",
        &output,
        "--zero-humidity = 100 at 100 °C and 1013.25 hPa: no such air exists",
    );
}

fn location_csv(dir: &std::path::Path, temperature: &str, humidity: &str) -> std::path::PathBuf {
    let path = dir.join("loc.csv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "SITE,TEMPERATURE,PRESSURE,HUMIDITY,ALTITUDE").unwrap();
    writeln!(f, "TestSite,{temperature},1013.25,{humidity},0").unwrap();
    path
}

#[test]
fn a_location_csv_air_is_refused_naming_its_column() {
    let dir = private_dir("loc");
    let csv = location_csv(dir.path(), "15", "150");
    let base = with(
        &["-u", "metric", "trajectory"],
        &with(LOAD, &["--max-range", "500"]),
    );
    let output = run(&with(
        &base,
        &["--location", csv.to_str().unwrap(), "--site", "TestSite"],
    ));
    assert_cli_refusal(
        "--location HUMIDITY 150",
        &output,
        &format!(
            "--location CSV column HUMIDITY = 150: relative humidity must be {}",
            percent_band()
        ),
    );
    let csv = location_csv(dir.path(), "100", "100");
    let output = run(&with(
        &base,
        &["--location", csv.to_str().unwrap(), "--site", "TestSite"],
    ));
    assert_cli_refusal(
        "--location 100 C / 100 %",
        &output,
        "--location CSV column HUMIDITY = 100 at 100 °C and 1013.25 hPa: no such air exists",
    );
}

fn profile_path(home: &std::path::Path, name: &str) -> std::path::PathBuf {
    home.join(".ballistics")
        .join("profiles")
        .join(format!("{name}.json"))
}

fn write_profile(home: &std::path::Path, name: &str, air: Value) {
    let mut profile = json!({
        "name": name, "velocity": 823.0, "bc": 0.475, "mass": 10.9, "diameter": 7.82,
        "drag_model": "G1", "units": "metric", "zero_distance": 100.0,
        "temperature_supplied": true, "pressure_supplied": true
    });
    for (key, value) in air.as_object().unwrap() {
        profile[key] = value.clone();
    }
    let path = profile_path(home, name);
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(path, profile.to_string()).unwrap();
}

#[test]
fn a_saved_profile_whose_air_cannot_exist_is_refused_naming_its_field() {
    let home = private_dir("prof");
    write_profile(
        home.path(),
        "boil",
        json!({"temperature": 100.0, "pressure": 1013.25, "humidity": 100.0}),
    );
    let cause = "profile 'boil' field humidity = 100 at 100 °C and 1013.25 hPa: no such air exists";
    let output = run_in(
        home.path(),
        &[
            "-u",
            "metric",
            "trajectory",
            "--saved-profile",
            "boil",
            "--max-range",
            "300",
        ],
    );
    assert_cli_refusal("trajectory --saved-profile boil", &output, cause);
    let output = run_in(
        home.path(),
        &[
            "-u",
            "metric",
            "dsf",
            "--saved-profile",
            "boil",
            "--range",
            "300",
            "--observed-drop",
            "2.0mil",
        ],
    );
    assert_cli_refusal("dsf --saved-profile boil", &output, cause);
    let output = run_in(
        home.path(),
        &[
            "-u",
            "metric",
            "plan-truing",
            "--profile",
            "boil",
            "--measurement-resolution",
            "0.25",
            "--candidate-ranges",
            "300,500,700",
        ],
    );
    assert_cli_refusal("plan-truing --profile boil", &output, cause);
}

#[test]
fn lead_exits_non_zero_when_its_air_cannot_exist() {
    let output = run(&with(
        &["-u", "metric", "lead"],
        &with(
            BULLET,
            &[
                "--target-speed",
                "4",
                "-o",
                "json",
                "--temperature",
                "100",
                "--humidity",
                "100",
            ],
        ),
    ));
    assert_cli_refusal("lead -o json at 100 C / 100 %", &output, CLI_IMPOSSIBLE);
}

#[test]
fn every_cli_humidity_refusal_states_the_band() {
    let base = with(
        &["-u", "metric", "trajectory"],
        &with(LOAD, &["--max-range", "300"]),
    );
    for (raw, reason) in [
        (
            "106",
            format!("relative humidity must be {}", percent_band()),
        ),
        (
            "-5.5",
            format!("relative humidity must be {}", percent_band()),
        ),
        (
            "nan",
            format!(
                "relative humidity must be a finite {}",
                &percent_band()[2..]
            ),
        ),
        (
            "abc",
            format!("not a number; relative humidity is {}", percent_band()),
        ),
    ] {
        let flag = format!("--humidity={raw}");
        let output = run(&with(&base, &[flag.as_str()]));
        assert_eq!(output.status.code(), Some(2), "clap refuses {raw}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(&reason),
            "--humidity {raw}: states the band\n  expected: {reason}\n  got: {stderr}"
        );
    }
}

// --- request files, solve-json, MCP -------------------------------------------------------------

fn solve_request(temperature_k: f64, pressure_pa: f64, relative_humidity: f64) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01134, "diameter_m": 0.00782, "drag_model": "G7",
                       "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0},
        "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
        "atmosphere": {"temperature_k": temperature_k, "pressure_pa": pressure_pa,
                       "relative_humidity": relative_humidity},
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": "rk45", "time_step_s": 0.001},
        "effects": {"magnus": false, "coriolis": false, "enhanced_spin_drift": false},
        "sampling": {"interval_m": 50.0}
    })
}

/// `relative_humidity = 1` at 373.15 K and 100000 Pa: the boiling point at 1000 hPa is below
/// 100 C, so saturated air there needs more vapour pressure than the air has.
fn impossible_solve() -> Value {
    solve_request(373.15, 100_000.0, 1.0)
}

const WIRE_IMPOSSIBLE: &str = "relative_humidity = 1 at 373.15 K and 100000 Pa: no such air exists";

#[test]
fn the_request_file_commands_print_the_wire_code_by_its_name() {
    let dir = private_dir("req");
    let path = dir.path().join("impossible.json");
    std::fs::write(&path, impossible_solve().to_string()).unwrap();
    let file = path.to_str().unwrap();
    for (label, args) in [
        (
            "explain",
            vec!["explain", "--a", file, "--b", file, "--ranges", "200"],
        ),
        (
            "tolerance",
            vec![
                "tolerance",
                "--request",
                file,
                "--range",
                "200",
                "--target",
                "rect:0.5x0.75",
                "--axis",
                "wind-speed",
                "--domain",
                "wind-speed=0:20",
            ],
        ),
        (
            "error-budget",
            vec![
                "error-budget",
                "--request",
                file,
                "--ranges",
                "200",
                "--sigma",
                "wind-speed=1.5",
            ],
        ),
    ] {
        let output = run(&args);
        assert_cli_refusal(
            label,
            &output,
            &format!("{file}: invalid_value at $.atmosphere.relative_humidity: {WIRE_IMPOSSIBLE}"),
        );
    }
}

/// A domain whose edge the solve refuses is refused naming the code by its wire name, never the
/// Rust variant (`solve failed (InvalidValue)`).
///
/// R3 (review C3): a humidity domain beyond 0 to 1 no longer reaches the solve at all. It is
/// refused up front for leaving the axis's physical range, naming the fraction it must be, the
/// way the humidity axis's unit is named everywhere now; before, `relative-humidity=0:1.2`
/// reached the solve and came back as `solve failed (invalid_value): ... outside the
/// perturbation domain [0, 1]`. A temperature domain reaching below 0 K still reaches the
/// solve, so it carries the wire-code check.
#[test]
fn tolerance_names_a_refused_solve_by_its_wire_code() {
    let dir = private_dir("tol");
    let path = dir.path().join("ok.json");
    std::fs::write(&path, solve_request(288.15, 101_325.0, 0.5).to_string()).unwrap();
    let tolerance = |axis: &str, domain: &str| {
        run(&[
            "tolerance",
            "--request",
            path.to_str().unwrap(),
            "--range",
            "200",
            "--target",
            "rect:0.5x0.75",
            "--axis",
            axis,
            "--domain",
            domain,
        ])
    };
    assert_cli_refusal(
        "tolerance --domain relative-humidity=0:1.2",
        &tolerance("relative-humidity", "relative-humidity=0:1.2"),
        "invalid tolerance-envelope domain for axis RelativeHumidity: the domain must lie \
         within the axis's physical range, a fraction from 0 to 1 (0.5 = 50 %)",
    );
    assert_cli_refusal(
        "tolerance --domain temperature=-10:400",
        &tolerance("temperature", "temperature=-10:400"),
        "solve failed (invalid_value): ",
    );
}

/// Every solve-json error code displays as its wire name, the name a reader branches on.
#[test]
fn a_solve_error_code_displays_as_its_wire_name() {
    use ballistics_engine::solve_json::SolveErrorCodeV1 as Code;
    for code in [
        Code::InvalidJson,
        Code::UnsupportedSchemaVersion,
        Code::UnknownField,
        Code::MissingField,
        Code::InvalidValue,
        Code::ConflictingFields,
        Code::ResourceLimit,
        Code::SolveFailed,
        Code::IoError,
        Code::InternalError,
    ] {
        assert_eq!(
            serde_json::to_value(code).unwrap(),
            code.to_string(),
            "{code:?}"
        );
    }
}

fn solve_json(request: &Value) -> (Option<i32>, Value) {
    let mut child = Command::new(BIN)
        .arg("solve-json")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(request.to_string().as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    (
        output.status.code(),
        serde_json::from_slice(&output.stdout).expect("solve-json writes JSON"),
    )
}

#[test]
fn solve_json_refuses_impossible_air_as_invalid_value_at_the_field_to_blame() {
    for (request, path, cause) in [
        (
            impossible_solve(),
            "$.atmosphere.relative_humidity",
            WIRE_IMPOSSIBLE.to_string(),
        ),
        // hPa sent as pascals: the pressure is what is wrong, so it is what is blamed.
        (
            solve_request(293.15, 1013.25, 0.5),
            "$.atmosphere.pressure_pa",
            "relative_humidity = 0.5 at 293.15 K and 1013.25 Pa: no such air exists".to_string(),
        ),
        (
            solve_request(293.15, 101_325.0, 1.06),
            "$.atmosphere.relative_humidity",
            format!(
                "relative_humidity = 1.06: relative humidity must be {}",
                fraction_band()
            ),
        ),
    ] {
        let (code, envelope) = solve_json(&request);
        assert_ne!(code, Some(0), "{envelope}");
        assert_eq!(envelope["error"]["code"], "invalid_value", "{envelope}");
        assert_eq!(envelope["error"]["path"], path, "{envelope}");
        let message = envelope["error"]["message"].as_str().unwrap_or_default();
        assert!(
            message.starts_with(&cause),
            "{message}\n  expected: {cause}"
        );
    }
}

#[test]
fn the_mcp_solve_tool_refuses_impossible_air_with_the_same_envelope() {
    use std::io::{BufRead, BufReader};
    let mut child = Command::new(BIN)
        .arg("mcp")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let mut stdout = BufReader::new(child.stdout.take().unwrap());
    let mut send = |message: Value| {
        writeln!(stdin, "{message}").unwrap();
        stdin.flush().unwrap();
    };
    let mut recv = || {
        let mut line = String::new();
        stdout.read_line(&mut line).unwrap();
        serde_json::from_str::<Value>(&line).unwrap()
    };
    send(
        json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
        "protocolVersion": "2025-06-18", "capabilities": {},
        "clientInfo": {"name": "r2", "version": "0"}}}),
    );
    recv();
    send(json!({"jsonrpc": "2.0", "method": "notifications/initialized"}));
    send(json!({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                "params": {"name": "solve", "arguments": impossible_solve()}}));
    let response = recv();
    drop(stdin); // EOF: the server exits cleanly
    let _ = child.wait();
    assert_eq!(response["result"]["isError"], true, "{response}");
    let payload: Value = serde_json::from_str(
        response["result"]["content"][0]["text"]
            .as_str()
            .unwrap_or_default(),
    )
    .unwrap();
    assert_eq!(payload["error"]["code"], "invalid_value", "{payload}");
    assert_eq!(
        payload["error"]["path"], "$.atmosphere.relative_humidity",
        "{payload}"
    );
    assert!(
        payload["error"]["message"]
            .as_str()
            .unwrap_or_default()
            .starts_with(WIRE_IMPOSSIBLE),
        "{payload}"
    );
}

// --- the bridge ---------------------------------------------------------------------------------

#[cfg(feature = "bridge")]
mod bridge {
    use super::*;

    fn call(command: &str, request: Value) -> Value {
        let envelope = json!({"api_version": 1, "command": command, "request": request});
        serde_json::from_str(&ballistics_engine::bridge::bridge_call(
            &envelope.to_string(),
        ))
        .unwrap()
    }

    /// Asserts the one envelope: `command_failed`, and the solve-json error in `details` with
    /// code `invalid_value` at `path`, its message leading with `cause`. Returns `details`.
    fn assert_envelope(label: &str, out: &Value, path: &str, cause: &str) -> Value {
        assert_eq!(out["ok"], false, "{label}: refused: {out}");
        assert_eq!(out["error"]["code"], "command_failed", "{label}: {out}");
        let details = &out["error"]["details"];
        assert_eq!(details["error"]["code"], "invalid_value", "{label}: {out}");
        assert_eq!(details["error"]["path"], path, "{label}: {out}");
        let inner = details["error"]["message"].as_str().unwrap_or_default();
        assert!(
            inner.starts_with(cause),
            "{label}: the inner message leads with the cause\n  expected: {cause}\n  got: {inner}"
        );
        // `solve` and `atmosphere.density_altitude` keep their documented outer messages
        // ("solve failed", "... request rejected"); every other command's leads with the cause.
        let outer = out["error"]["message"].as_str().unwrap_or_default();
        if !label.starts_with("solve") && !label.starts_with("atmosphere.") {
            assert!(
                outer.contains(cause),
                "{label}: the outer message carries the cause: {outer}"
            );
        }
        assert!(debug_forms(outer).is_empty(), "{label}: {outer}");
        details.clone()
    }

    fn card(command: &str, temperature_f: f64, humidity: f64) -> Value {
        let mut request = json!({
            "units": "imperial", "muzzle_velocity": 2600.0, "ballistic_coefficient": 0.243,
            "mass": 175.0, "diameter": 0.308, "drag_model": "g7", "zero_distance": 100.0,
            "start": 100.0, "end": 300.0, "step": 100.0,
            "temperature": temperature_f, "pressure": 29.92, "humidity": humidity
        });
        if command == "card.wind" {
            request["wind_speeds"] = json!([5.0, 10.0]);
        }
        request
    }

    fn model(temperature_f: f64, humidity: f64) -> Value {
        json!({
            "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
            "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
            "zero_distance_yd": 100.0, "sight_height_in": 2.0,
            "temperature_f": temperature_f, "pressure_inhg": 29.92,
            "humidity_pct": humidity, "altitude_ft": 0.0
        })
    }

    fn truing(command: &str, temperature_f: f64, humidity: f64) -> Value {
        match command {
            "true.fit" => json!({
                "model": model(temperature_f, humidity), "drop_unit": "mil",
                "observations": [
                    {"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1},
                    {"range_yd": 600.0, "drop": 4.8626, "sigma": 0.1},
                    {"range_yd": 800.0, "drop": 7.6598, "sigma": 0.1}
                ],
                "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                           "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
                "predictions": []
            }),
            "true.dsf" => {
                let mut model = model(temperature_f, humidity);
                model["muzzle_velocity_fps"] = json!(1600.0);
                json!({"model": model, "range_yd": 200.0, "observed_drop": 12.0,
                       "drop_unit": "mil"})
            }
            "true.plan" => json!({
                "model": model(temperature_f, humidity),
                "candidate_ranges_yd": [200.0, 400.0, 600.0, 800.0],
                "observation_count": 3, "minimum_separation_yd": 100.0,
                "measurement_sigma_1sd": 0.1, "drop_unit": "mil"
            }),
            "true.wind" => json!({
                "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
                "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
                "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
                "sight_height_in": 2.0, "temperature_f": temperature_f, "pressure_inhg": 29.92,
                "humidity_pct": humidity, "altitude_ft": 0.0,
                "twist": {"rate_in": 11.0, "right_hand": true},
                "earth": null, "called_crosswind_mph": null
            }),
            other => panic!("{other}"),
        }
    }

    #[test]
    fn solve_and_density_altitude_share_the_envelope() {
        for command in ["solve", "atmosphere.density_altitude"] {
            let request = |atmosphere: Value| match command {
                "solve" => {
                    let mut request = solve_request(0.0, 0.0, 0.0);
                    request["atmosphere"] = atmosphere;
                    request
                }
                _ => json!({"atmosphere": atmosphere}),
            };
            assert_envelope(
                &format!("{command} impossible"),
                &call(command, request(impossible_solve()["atmosphere"].clone())),
                "$.atmosphere.relative_humidity",
                WIRE_IMPOSSIBLE,
            );
            assert_envelope(
                &format!("{command} 1.06"),
                &call(
                    command,
                    request(json!({"temperature_k": 293.15, "pressure_pa": 101325.0,
                                   "relative_humidity": 1.06})),
                ),
                "$.atmosphere.relative_humidity",
                &format!(
                    "relative_humidity = 1.06: relative humidity must be {}",
                    fraction_band()
                ),
            );
        }
    }

    /// One rule decides which reading impossible air is blamed on (gate R2): the temperature
    /// when it is outside 173.15-373.15 K; then, at a station pressure outside the troposphere,
    /// the pressure when it was sent as a station pressure (or as a QNH that is itself outside
    /// it) and the altitude when it was omitted or was a plausible QNH; otherwise the humidity.
    /// That is the order `atmosphere.density_altitude`'s own range checks blame in, so the two
    /// commands give the same code and path for every one, as its documentation says.
    #[test]
    fn solve_and_density_altitude_blame_the_same_reading_for_the_same_air() {
        for (atmosphere, path) in [
            (
                json!({"temperature_k": 380.0, "pressure_pa": 101325.0,
                       "relative_humidity": 1.0}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"temperature_k": 400.0, "pressure_pa": 101325.0,
                       "relative_humidity": 0.9}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"temperature_k": 170.0, "pressure_pa": 0.0001,
                       "relative_humidity": 1.0}),
                "$.atmosphere.temperature_k",
            ),
            (
                json!({"temperature_k": 373.15, "pressure_pa": 100000.0,
                       "relative_humidity": 1.0}),
                "$.atmosphere.relative_humidity",
            ),
            (
                json!({"temperature_k": 293.15, "pressure_pa": 1013.25,
                       "relative_humidity": 0.5}),
                "$.atmosphere.pressure_pa",
            ),
            (
                json!({"altitude_m": 20000.0, "temperature_k": 310.0, "relative_humidity": 1.0}),
                "$.atmosphere.altitude_m",
            ),
            (
                json!({"altitude_m": 20000.0, "temperature_k": 310.0, "pressure_pa": 101325.0,
                       "pressure_reference": "qnh", "relative_humidity": 1.0}),
                "$.atmosphere.altitude_m",
            ),
            (
                json!({"altitude_m": 1500.0, "temperature_k": 293.15, "pressure_pa": 1013.25,
                       "pressure_reference": "qnh", "relative_humidity": 0.5}),
                "$.atmosphere.pressure_pa",
            ),
        ] {
            let mut solve = solve_request(0.0, 0.0, 0.0);
            solve["atmosphere"] = atmosphere.clone();
            for (command, out) in [
                ("solve", call("solve", solve)),
                (
                    "atmosphere.density_altitude",
                    call(
                        "atmosphere.density_altitude",
                        json!({ "atmosphere": atmosphere }),
                    ),
                ),
            ] {
                assert_eq!(out["ok"], false, "{command} {atmosphere}: {out}");
                assert_eq!(out["error"]["code"], "command_failed", "{command}: {out}");
                let error = &out["error"]["details"]["error"];
                assert_eq!(
                    error["code"], "invalid_value",
                    "{command} {atmosphere}: {out}"
                );
                assert_eq!(error["path"], path, "{command} {atmosphere}: {out}");
            }
        }
        // The solve says which unit to check: the altitude's, for a pressure nobody sent.
        let mut solve = solve_request(0.0, 0.0, 0.0);
        solve["atmosphere"] =
            json!({"altitude_m": 20000.0, "temperature_k": 310.0, "relative_humidity": 1.0});
        let out = call("solve", solve);
        let message = out["error"]["details"]["error"]["message"]
            .as_str()
            .unwrap_or_default();
        assert!(
            message.starts_with("relative_humidity = 1 at 310 K and ")
                && message.contains("no pressure_pa was given")
                && message.contains("check the altitude's unit"),
            "{message}"
        );
    }

    /// Every other wire command blames the same reading for the same air, at its own member: a
    /// temperature outside 173.15-373.15 K is the temperature's fault, and an omitted pressure
    /// that leaves the station above the troposphere is the altitude's (gate R2).
    #[test]
    fn every_wire_command_blames_the_reading_solve_blames() {
        let hot = "humidity = 100 at 230 °F and 29.92 inHg: no such air exists";
        assert_envelope(
            "card.come_ups 230 F / 100 %",
            &call("card.come_ups", card("card.come_ups", 230.0, 100.0)),
            "$.temperature",
            hot,
        );
        let mut high = card("card.come_ups", 104.0, 100.0);
        high.as_object_mut().unwrap().remove("pressure");
        high["altitude"] = json!(20000.0);
        let details = assert_envelope(
            "card.come_ups at 20 km, no pressure",
            &call("card.come_ups", high),
            "$.altitude",
            "humidity = 100 at 104 °F and ",
        );
        let message = details["error"]["message"].as_str().unwrap_or_default();
        assert!(
            message.contains("no pressure was given")
                && message.contains("check the altitude's unit"),
            "{message}"
        );
        for (command, at) in [("true.fit", "$.model."), ("true.wind", "$.")] {
            assert_envelope(
                &format!("{command} 230 F / 100 %"),
                &call(command, truing(command, 230.0, 100.0)),
                &format!("{at}temperature_f"),
                "humidity_pct = 100 at 230 °F and 29.92 inHg: no such air exists",
            );
        }
        let mut high = truing("true.fit", 104.0, 100.0);
        high["model"]
            .as_object_mut()
            .unwrap()
            .remove("pressure_inhg");
        high["model"]["altitude_ft"] = json!(65617.0);
        assert_envelope(
            "true.fit at 20 km, no pressure",
            &call("true.fit", high),
            "$.model.altitude_ft",
            "humidity_pct = 100 at 104 °F and ",
        );
    }

    /// `profile.export_a7p`'s message is the refusal itself, cause first, with the field named
    /// once (gate R2 note: it read `humidity: profile 'boil' field humidity = 100 ...`).
    #[cfg(feature = "profile-export")]
    #[test]
    fn profile_export_a7p_names_the_refused_field_once() {
        let out = call("profile.export_a7p", stored_profile(100.0, 100.0));
        let message = out["error"]["message"].as_str().unwrap_or_default();
        assert!(
            message.starts_with(
                "profile 'boil' field humidity = 100 at 100 °C and 1013.25 hPa: no such air exists"
            ),
            "{message}"
        );
    }

    /// A stored card is reprinted only for air its request could fly: `card.pdf` refuses the
    /// request's impossible air before it reads the stored card, like every card command.
    #[cfg(feature = "pdf")]
    #[test]
    fn card_pdf_refuses_a_reprint_whose_request_air_cannot_exist() {
        let table = call("card.range_table", card("card.range_table", 212.0, 50.0));
        assert_eq!(table["ok"], true, "{table}");
        let mut reprint = card("card.pdf", 212.0, 50.0);
        reprint["stored_card"] = json!({"card": table["result"]});
        let out = call("card.pdf", reprint.clone());
        assert_eq!(out["ok"], true, "212 F at 50 % reprints: {out}");
        reprint["humidity"] = json!(100.0);
        assert_envelope(
            "card.pdf reprint at 212 F / 100 %",
            &call("card.pdf", reprint),
            "$.humidity",
            "humidity = 100 at 212 °F and 29.92 inHg: no such air exists",
        );
    }

    #[test]
    fn every_card_command_refuses_with_the_envelope_at_its_humidity() {
        let mut commands = vec!["card.come_ups", "card.range_table", "card.wind"];
        if cfg!(feature = "pdf") {
            commands.push("card.pdf");
        }
        for command in commands {
            assert_envelope(
                &format!("{command} 212 F / 100 %"),
                &call(command, card(command, 212.0, 100.0)),
                "$.humidity",
                "humidity = 100 at 212 °F and 29.92 inHg: no such air exists",
            );
            assert_envelope(
                &format!("{command} 106 %"),
                &call(command, card(command, 59.0, 106.0)),
                "$.humidity",
                &format!(
                    "humidity = 106: relative humidity must be {}",
                    percent_band()
                ),
            );
        }
    }

    #[test]
    fn every_truing_command_refuses_with_the_envelope_at_its_humidity() {
        for (command, path, field) in [
            ("true.fit", "$.model.humidity_pct", "humidity_pct"),
            ("true.dsf", "$.model.humidity_pct", "humidity_pct"),
            ("true.plan", "$.model.humidity_pct", "humidity_pct"),
            ("true.wind", "$.humidity_pct", "humidity_pct"),
        ] {
            let details = assert_envelope(
                &format!("{command} 212 F / 100 %"),
                &call(command, truing(command, 212.0, 100.0)),
                path,
                &format!("{field} = 100 at 212 °F and 29.92 inHg: no such air exists"),
            );
            let flat = match command {
                "true.dsf" => Some("invalid_input"),
                "true.plan" => Some("invalid_request"),
                _ => None,
            };
            if let Some(reason) = flat {
                // The family's flat reason stays: impossible air is a bad input, never a
                // shortage of reachable ranges or a forward-model failure.
                assert_eq!(details["reason"], reason, "{command}: {details}");
            }
            assert_envelope(
                &format!("{command} 106 %"),
                &call(command, truing(command, 59.0, 106.0)),
                path,
                &format!(
                    "{field} = 106: relative humidity must be {}",
                    percent_band()
                ),
            );
        }
    }

    /// A failure of `true.plan`'s whole candidate batch (here, a zero the model cannot reach) is
    /// the request's error, with its cause first, never every candidate "unreachable": that
    /// would point a caller at range selection for a problem no range choice can fix.
    #[test]
    fn true_plan_reports_a_whole_batch_failure_as_the_requests_with_its_cause() {
        let mut request = truing("true.plan", 59.0, 50.0);
        request["model"]["muzzle_velocity_fps"] = json!(1100.0);
        request["model"]["ballistic_coefficient"] = json!(0.1);
        request["model"]["zero_distance_yd"] = json!(3000.0);
        let out = call("true.plan", request);
        assert_eq!(out["ok"], false, "{out}");
        let details = &out["error"]["details"];
        assert_eq!(details["reason"], "invalid_request", "{out}");
        assert_eq!(details["rejected_candidates"], json!([]), "{out}");
        let message = out["error"]["message"].as_str().unwrap_or_default();
        assert!(
            message.starts_with("true.plan failed: Cannot find zero angle")
                && message.ends_with("; no candidate range could be evaluated"),
            "{message}"
        );
    }

    fn stored_profile(temperature_c: f64, humidity: f64) -> Value {
        json!({
            "name": "boil", "velocity": 823.0, "bc": 0.475, "mass": 10.9, "diameter": 7.82,
            "drag_model": "G1", "units": "metric", "zero_distance": 100.0,
            "temperature": temperature_c, "pressure": 1013.25, "humidity": humidity,
            "temperature_supplied": true, "pressure_supplied": true,
            "bullet_length": 31.0, "twist_rate": 254.0, "sight_height": 50.0
        })
    }

    #[test]
    fn profile_validate_applies_the_impossible_air_rule_a_flight_does() {
        let out = call("profile.validate", stored_profile(100.0, 100.0));
        assert_eq!(out["ok"], true, "{out}");
        assert_eq!(out["result"]["valid"], false, "{out}");
        let warnings = out["result"]["warnings"].to_string();
        assert!(
            warnings.contains(
                "profile 'boil' field humidity = 100 at 100 °C and 1013.25 hPa: no such air exists"
            ),
            "{warnings}"
        );
        let out = call("profile.validate", stored_profile(100.0, 50.0));
        assert_eq!(out["result"]["valid"], true, "100 C at 50 % exists: {out}");
        let out = call("profile.validate", stored_profile(15.0, 150.0));
        assert!(
            out["result"]["warnings"]
                .to_string()
                .contains(&percent_band()),
            "{out}"
        );
    }

    #[cfg(feature = "profile-export")]
    #[test]
    fn profile_export_a7p_refuses_air_a_flight_refuses() {
        let details = assert_envelope(
            "profile.export_a7p 100 C / 100 %",
            &call("profile.export_a7p", stored_profile(100.0, 100.0)),
            "$.humidity",
            "profile 'boil' field humidity = 100 at 100 °C and 1013.25 hPa: no such air exists",
        );
        assert_eq!(details["reason"], "field", "{details}");
        assert_eq!(details["field"], "humidity", "{details}");
        assert_envelope(
            "profile.export_a7p 106 %",
            &call("profile.export_a7p", stored_profile(15.0, 106.0)),
            "$.humidity",
            &format!(
                "profile 'boil' field humidity = 106: relative humidity must be {}",
                percent_band()
            ),
        );
    }

    #[cfg(feature = "profile-import")]
    #[test]
    fn profile_import_a7p_refuses_air_a_flight_refuses() {
        let a7p = |temperature_c: i64, humidity: i64| {
            let bytes = super::a7p::file(temperature_c, humidity);
            json!({"a7p_base64": super::a7p::base64(&bytes)})
        };
        assert_envelope(
            "profile.import_a7p 100 C / 100 %",
            &call("profile.import_a7p", a7p(100, 100)),
            "a7p:c_zero_air_humidity",
            "c_zero_air_humidity = 100 at 100 °C and 1013.2 hPa: no such air exists",
        );
        assert_envelope(
            "profile.import_a7p 106 %",
            &call("profile.import_a7p", a7p(15, 106)),
            "a7p:c_zero_air_humidity",
            &format!(
                "c_zero_air_humidity = 106: relative humidity must be {}",
                percent_band()
            ),
        );
        let out = call("profile.import_a7p", a7p(100, 50));
        assert_eq!(out["ok"], true, "100 C at 50 % exists: {out}");
    }
}

// --- the profile store --------------------------------------------------------------------------

#[test]
fn profile_show_displays_a_profile_it_would_refuse_to_fly_with_the_refusal_as_a_warning() {
    let home = private_dir("show");
    write_profile(
        home.path(),
        "wet",
        json!({"temperature": 15.0, "pressure": 1013.25, "humidity": 150.0}),
    );
    let output = run_in(home.path(), &["profile", "show", "wet"]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stdout.contains("150"),
        "the stored value is shown: {stdout}"
    );
    assert!(
        stderr.contains(&format!(
            "warning: profile 'wet' field humidity = 150: relative humidity must be {}",
            percent_band()
        )),
        "the refusal is a warning beside it: {stderr}"
    );
    assert!(debug_forms(&stderr).is_empty(), "{stderr}");
    let output = run_in(home.path(), &["profile", "zero-set", "list", "wet"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn profile_save_refuses_air_that_cannot_exist() {
    let home = private_dir("save");
    let output = run_in(
        home.path(),
        &[
            "-u",
            "metric",
            "profile",
            "save",
            "boil",
            "--velocity",
            "823",
            "--bc",
            "0.475",
            "--mass",
            "10.9",
            "--diameter",
            "7.82",
            "--temperature",
            "100",
            "--humidity",
            "100",
        ],
    );
    assert_cli_refusal("profile save at 100 C / 100 %", &output, CLI_IMPOSSIBLE);
    assert!(
        !profile_path(home.path(), "boil").exists(),
        "nothing was written"
    );
}

#[cfg(feature = "profile-import")]
#[test]
fn profile_import_refuses_an_a7p_whose_air_cannot_exist() {
    let home = private_dir("import");
    let file = home.path().join("boil.a7p");
    std::fs::write(&file, a7p::file(100, 100)).unwrap();
    let output = run_in(
        home.path(),
        &[
            "profile",
            "import",
            file.to_str().unwrap(),
            "--name",
            "boil",
        ],
    );
    assert_cli_refusal(
        "profile import boil.a7p",
        &output,
        "c_zero_air_humidity = 100 at 100 °C and 1013.2 hPa: no such air exists",
    );
    assert!(!profile_path(home.path(), "boil").exists());
}

#[cfg(feature = "profile-import")]
mod a7p {
    use ballistics_engine::profile_import::wrap_payload;

    fn varint(mut v: u64, out: &mut Vec<u8>) {
        loop {
            let byte = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(byte);
                return;
            }
            out.push(byte | 0x80);
        }
    }
    fn int(number: u32, value: i64, out: &mut Vec<u8>) {
        varint(u64::from(number) << 3, out);
        varint(value as u64, out);
    }
    fn text(number: u32, s: &str, out: &mut Vec<u8>) {
        varint((u64::from(number) << 3) | 2, out);
        varint(s.len() as u64, out);
        out.extend_from_slice(s.as_bytes());
    }
    fn bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
        varint((u64::from(number) << 3) | 2, out);
        varint(payload.len() as u64, out);
        out.extend_from_slice(payload);
    }

    /// A G7 .a7p whose zero air is `temperature_c` (field 15), 1013.2 hPa (field 16, tenths)
    /// and `humidity` percent (field 17).
    pub(super) fn file(temperature_c: i64, humidity: i64) -> Vec<u8> {
        let mut p = Vec::new();
        text(1, "BOIL 6.5CM", &mut p);
        int(9, 50, &mut p);
        int(10, 800, &mut p);
        int(11, 8230, &mut p);
        int(15, temperature_c, &mut p);
        int(16, 10132, &mut p);
        int(17, humidity, &mut p);
        int(20, 264, &mut p);
        int(21, 1400, &mut p);
        int(22, 1360, &mut p);
        int(24, 1, &mut p);
        let mut packed = Vec::new();
        varint(10_000, &mut packed);
        bytes(26, &packed, &mut p);
        let mut row = Vec::new();
        int(1, 3260, &mut row);
        int(2, 8230, &mut row);
        bytes(27, &row, &mut p);
        let mut payload = Vec::new();
        bytes(1, &p, &mut payload);
        wrap_payload(&payload)
    }

    #[cfg(feature = "bridge")]
    pub(super) fn base64(input: &[u8]) -> String {
        const ALPHABET: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for chunk in input.chunks(3) {
            let b = [
                chunk[0],
                *chunk.get(1).unwrap_or(&0),
                *chunk.get(2).unwrap_or(&0),
            ];
            let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
            for i in 0..4 {
                if i <= chunk.len() {
                    out.push(ALPHABET[((n >> (18 - 6 * i)) & 63) as usize] as char);
                } else {
                    out.push('=');
                }
            }
        }
        out
    }
}

// --- the library's own words ------------------------------------------------------------------

#[test]
fn the_library_states_the_band_in_every_refusal() {
    use ballistics_engine::RelativeHumidity;
    let percent = RelativeHumidity::try_percent(106.0)
        .unwrap_err()
        .to_string();
    assert_eq!(
        percent,
        format!(
            "humidity_percent = 106: relative humidity must be {}",
            percent_band()
        )
    );
    let fraction = RelativeHumidity::try_fraction(1.06)
        .unwrap_err()
        .to_string();
    assert_eq!(
        fraction,
        format!(
            "relative_humidity = 1.06: relative humidity must be {}",
            fraction_band()
        )
    );
    let nan = RelativeHumidity::try_percent(f64::NAN)
        .unwrap_err()
        .to_string();
    assert_eq!(
        nan,
        format!(
            "humidity_percent = NaN: relative humidity must be a finite {}",
            &percent_band()[2..]
        )
    );
}
