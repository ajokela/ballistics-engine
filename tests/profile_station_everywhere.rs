//! Every command that takes a saved profile flies the profile's stored weather (0.50.0).
//!
//! `profile save` stores a temperature, a pressure and its type (absolute or QNH), a humidity,
//! an altitude and optionally a density altitude. Before 0.50.0 only `trajectory
//! --saved-profile`, `plan-truing --profile` and `dsf --saved-profile` flew them; `come-ups`,
//! `range-table`, `wind-card`, `mpbr`, `lead`, `compare`, `hold-corridor`, `mark-to-range`,
//! `bdc-match`, `optimal-zero`, `adaptive-card`, `reticle hold` and `stability` loaded the
//! profile's bullet and zero and flew the flags' defaults, so one profile gave different dope
//! on different commands. Now each station field follows one rule on every one of them: a
//! typed flag wins, then the stored value, then the command's documented default.
//!
//! What these tests hold, for every consumer:
//! 1. a profile with stored weather flies it: the run is byte-identical to the same profile
//!    flown with that weather typed;
//! 2. a typed flag replaces only its own field, and a typed value equal to the default still
//!    wins;
//! 3. a profile without stored weather flies exactly what no profile flies;
//! 4. the zero is solved in the same air (the zero day follows the same rule);
//!
//! and that a stored value is refused only where it is flown: a stored humidity a command
//! never reads (because a flag replaced it, or because the command flies no humidity) no
//! longer stops the run.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every run points `HOME` at a
//! private temp dir.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

use serde_json::Value;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

/// A private `HOME`, removed when dropped.
struct Home(PathBuf);

impl Drop for Home {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

impl Home {
    fn path(&self) -> &Path {
        &self.0
    }

    fn profile_path(&self, name: &str) -> PathBuf {
        self.0
            .join(".ballistics")
            .join("profiles")
            .join(format!("{name}.json"))
    }
}

fn private_home() -> Home {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-profile-station-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    Home(dir)
}

fn strings(parts: &[&[&str]]) -> Vec<String> {
    parts
        .iter()
        .flat_map(|part| part.iter().map(|s| s.to_string()))
        .collect()
}

fn run(home: &Path, args: &[String]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn stdout(home: &Path, args: &[String]) -> String {
    let output = run(home, args);
    assert!(
        output.status.success(),
        "{args:?} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// The bullet and rifle every profile here stores, as `profile save` flags.
const BULLET: &[&str] = &[
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--drag-model",
    "g1",
];
const RIFLE: &[&str] = &[
    "--zero-distance",
    "100",
    "--sight-height",
    "2",
    "--twist-rate",
    "10",
    "--bullet-length",
    "1.24",
];

/// Weather stored by `profile save` (and typed, to fly the same air without it).
const WX: [(&str, &str); 4] = [
    ("--temperature", "95"),
    ("--pressure", "25.5"),
    ("--humidity", "90"),
    ("--altitude", "5000"),
];
/// A stored QNH reading: its type is stored with it.
const QNH: [(&str, &str); 5] = [
    ("--temperature", "40"),
    ("--pressure", "30.10"),
    ("--pressure-type", "qnh"),
    ("--humidity", "20"),
    ("--altitude", "3000"),
];
/// A stored density altitude, with the temperature and humidity it was saved with.
const DA: [(&str, &str); 3] = [
    ("--density-altitude", "7000"),
    ("--temperature", "80"),
    ("--humidity", "30"),
];

fn flags(pairs: &[(&str, &str)]) -> Vec<String> {
    pairs
        .iter()
        .flat_map(|(flag, value)| [flag.to_string(), value.to_string()])
        .collect()
}

/// `pairs` with every field in `replace` replaced by its new value, and the fields in `drop`
/// left out.
fn with_fields(pairs: &[(&str, &str)], replace: &[(&str, &str)], drop: &[&str]) -> Vec<String> {
    let mut out = Vec::new();
    for (flag, value) in pairs {
        if drop.contains(flag) {
            continue;
        }
        let value = replace
            .iter()
            .find(|(f, _)| f == flag)
            .map_or(*value, |(_, v)| v);
        out.push(flag.to_string());
        out.push(value.to_string());
    }
    for (flag, value) in replace {
        if !pairs.iter().any(|(f, _)| f == flag) {
            out.push(flag.to_string());
            out.push(value.to_string());
        }
    }
    out
}

fn save(home: &Home, name: &str, weather: &[String]) {
    let mut args = strings(&[&["profile", "save", name], BULLET, RIFLE]);
    args.extend(weather.iter().cloned());
    let output = run(home.path(), &args);
    assert!(
        output.status.success(),
        "profile save {name} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// A stored profile's JSON, edited by `edit` and written back (a hand-edited file).
fn edit_profile(home: &Home, from: &str, to: &str, edit: impl FnOnce(&mut Value)) {
    let text = std::fs::read_to_string(home.profile_path(from)).unwrap();
    let mut profile: Value = serde_json::from_str(&text).unwrap();
    profile["name"] = Value::String(to.to_string());
    edit(&mut profile);
    std::fs::write(home.profile_path(to), profile.to_string()).unwrap();
}

/// The fixture files `hold-corridor` and `reticle hold` read.
fn fixtures(home: &Home) -> (String, String) {
    let scenarios = home.path().join("scenarios.json");
    std::fs::write(
        &scenarios,
        r#"{"version":1,"nominal":"low","scenarios":[{"name":"low","segments":["4:90:1000"]},{"name":"high","segments":["14:90:1000"]},{"name":"switch","segments":["10:90:400","8:270:1000"]}]}"#,
    )
    .unwrap();
    let reticle = stdout(
        home.path(),
        &strings(&[&[
            "reticle",
            "generate",
            "tree",
            "--rows",
            "8",
            "--row-spacing",
            "1.0",
            "--spread-step",
            "0.5",
            "-o",
            "json",
        ]]),
    );
    let reticle_path = home.path().join("reticle.json");
    std::fs::write(&reticle_path, reticle).unwrap();
    (
        scenarios.to_str().unwrap().to_string(),
        reticle_path.to_str().unwrap().to_string(),
    )
}

/// One command that takes a saved profile.
#[derive(Clone)]
struct Consumer {
    name: &'static str,
    /// The command's own arguments, without the profile flag or any station flag.
    args: Vec<String>,
    /// `--profile` or `--saved-profile`.
    profile_flag: &'static str,
    /// Whether the command flies a humidity (stability's Miller correction is dry).
    flies_humidity: bool,
    /// The same run with the bullet and rifle typed instead of loaded, when the command can run
    /// without a profile.
    bare: Option<Vec<String>>,
}

impl Consumer {
    fn with_profile(&self, profile: &str, extra: &[String]) -> Vec<String> {
        let mut args = self.args.clone();
        args.push(self.profile_flag.to_string());
        args.push(profile.to_string());
        args.extend(extra.iter().cloned());
        args
    }

    fn station_pairs<'a>(&self, pairs: &[(&'a str, &'a str)]) -> Vec<(&'a str, &'a str)> {
        pairs
            .iter()
            .copied()
            .filter(|(flag, _)| self.flies_humidity || *flag != "--humidity")
            .collect()
    }
}

/// Every consumer of a saved profile that flies air, in imperial units.
fn consumers(scenarios: &str, reticle: &str) -> Vec<Consumer> {
    let load = &["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"][..];
    let bare_rifle = &["--zero-distance", "100", "--sight-height", "2"][..];
    let json = &["-o", "json"][..];
    let c = |name, args: Vec<String>, flag, humid, bare: Option<Vec<String>>| Consumer {
        name,
        args,
        profile_flag: flag,
        flies_humidity: humid,
        bare,
    };
    vec![
        c(
            "trajectory",
            strings(&[
                &["trajectory", "--max-range", "600", "--auto-zero", "100"],
                json,
            ]),
            "--saved-profile",
            true,
            None,
        ),
        c(
            "plan-truing",
            strings(&[
                &[
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500,700",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500,700",
                    "--drag-model",
                    "g1",
                ],
                load,
                bare_rifle,
                json,
            ])),
        ),
        c(
            "come-ups",
            strings(&[
                &["come-ups", "--zero-distance", "100", "--end", "600"],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &["come-ups", "--zero-distance", "100", "--end", "600"],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "range-table",
            strings(&[
                &["range-table", "--zero-distance", "100", "--end", "600"],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &["range-table", "--zero-distance", "100", "--end", "600"],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "wind-card",
            strings(&[
                &["wind-card", "--zero-distance", "100", "--end", "600"],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &["wind-card", "--zero-distance", "100", "--end", "600"],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "mpbr",
            strings(&[&["mpbr", "--vital-zone", "8"], json]),
            "--profile",
            true,
            Some(strings(&[
                &["mpbr", "--vital-zone", "8"],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "lead",
            strings(&[&["lead", "--target-speed", "4", "--end", "600"], json]),
            "--profile",
            true,
            Some(strings(&[
                &["lead", "--target-speed", "4", "--end", "600"],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "hold-corridor",
            strings(&[
                &[
                    "hold-corridor",
                    "--scenarios",
                    scenarios,
                    "--ranges",
                    "200,400,600",
                    "--target",
                    "rect:12x18",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "hold-corridor",
                    "--scenarios",
                    scenarios,
                    "--ranges",
                    "200,400,600",
                    "--target",
                    "rect:12x18",
                ],
                load,
                bare_rifle,
                json,
            ])),
        ),
        c(
            "mark-to-range",
            strings(&[&["mark-to-range", "--mark", "2.0", "--mark", "5.0"], json]),
            "--profile",
            true,
            Some(strings(&[
                &["mark-to-range", "--mark", "2.0", "--mark", "5.0"],
                load,
                bare_rifle,
                json,
            ])),
        ),
        c(
            "bdc-match",
            strings(&[
                &[
                    "bdc-match",
                    "--mark-range",
                    "2:300",
                    "--mark-range",
                    "4:500",
                    "--mark-range",
                    "6:700",
                    "--reference-mag",
                    "12",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "bdc-match",
                    "--mark-range",
                    "2:300",
                    "--mark-range",
                    "4:500",
                    "--mark-range",
                    "6:700",
                    "--reference-mag",
                    "12",
                ],
                load,
                bare_rifle,
                json,
            ])),
        ),
        c(
            "optimal-zero",
            strings(&[
                &[
                    "optimal-zero",
                    "--target",
                    "200:10",
                    "--target",
                    "400:12",
                    "--target",
                    "600:18",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "optimal-zero",
                    "--target",
                    "200:10",
                    "--target",
                    "400:12",
                    "--target",
                    "600:18",
                ],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "adaptive-card",
            strings(&[
                &[
                    "adaptive-card",
                    "--zero-distance",
                    "100",
                    "--start",
                    "100",
                    "--end",
                    "600",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "adaptive-card",
                    "--zero-distance",
                    "100",
                    "--start",
                    "100",
                    "--end",
                    "600",
                ],
                load,
                &["--sight-height", "2"],
                json,
            ])),
        ),
        c(
            "reticle-hold",
            strings(&[
                &[
                    "reticle",
                    "hold",
                    "--reticle-json",
                    reticle,
                    "--range",
                    "500",
                ],
                json,
            ]),
            "--profile",
            true,
            Some(strings(&[
                &[
                    "reticle",
                    "hold",
                    "--reticle-json",
                    reticle,
                    "--range",
                    "500",
                ],
                load,
                bare_rifle,
                json,
            ])),
        ),
        c(
            "stability",
            strings(&[
                &["stability", "--twist-rate", "10", "--length", "1.24"],
                json,
            ]),
            "--profile",
            false,
            Some(strings(&[
                &["stability", "--twist-rate", "10", "--length", "1.24"],
                &["-m", "168", "-d", "0.308"],
                json,
            ])),
        ),
    ]
}

/// A home with the shared profiles saved, and the consumers that read them.
fn setup() -> (Home, Vec<Consumer>) {
    let home = private_home();
    let (scenarios, reticle) = fixtures(&home);
    save(&home, "plain", &[]);
    save(&home, "wx", &flags(&WX));
    save(&home, "qnh", &flags(&QNH));
    save(&home, "da", &flags(&DA));
    let list = consumers(&scenarios, &reticle);
    (home, list)
}

/// Run `cases` (a label and two argument lists each) in parallel and return, per case, the
/// label and both stdouts.
fn run_pairs(
    home: &Path,
    cases: Vec<(String, Vec<String>, Vec<String>)>,
) -> Vec<(String, String, String)> {
    std::thread::scope(|scope| {
        let jobs: Vec<_> = cases
            .into_iter()
            .map(|(label, a, b)| {
                scope.spawn(move || {
                    let left = stdout(home, &a);
                    let right = stdout(home, &b);
                    (format!("{label}\n  {a:?}\n  {b:?}"), left, right)
                })
            })
            .collect();
        jobs.into_iter().map(|job| job.join().unwrap()).collect()
    })
}

fn assert_all_equal(results: &[(String, String, String)]) {
    let failures: Vec<&String> = results
        .iter()
        .filter(|(_, a, b)| a != b)
        .map(|(label, _, _)| label)
        .collect();
    assert!(
        failures.is_empty(),
        "{} of {} runs differ where they must be byte-identical:\n{}",
        failures.len(),
        results.len(),
        failures
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
}

fn assert_all_differ(results: &[(String, String, String)]) {
    let failures: Vec<&String> = results
        .iter()
        .filter(|(_, a, b)| a == b)
        .map(|(label, _, _)| label)
        .collect();
    assert!(
        failures.is_empty(),
        "{} of {} runs are identical where the stored weather must move them:\n{}",
        failures.len(),
        results.len(),
        failures
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// (1) A profile's stored weather flies: the run is the same profile's bullet with that weather
/// typed. The stored temperature and pressure were typed when the profile was saved, so they
/// are readings and fly as typed ones do. And the stored weather moves the output.
#[test]
fn every_consumer_flies_the_profiles_stored_weather() {
    let (home, list) = setup();
    let mut same = Vec::new();
    let mut moved = Vec::new();
    for consumer in &list {
        for (profile, pairs) in [("wx", &WX[..]), ("qnh", &QNH[..])] {
            let typed = flags(&consumer.station_pairs(pairs));
            same.push((
                format!(
                    "{}: --profile {profile} flies its stored weather",
                    consumer.name
                ),
                consumer.with_profile(profile, &[]),
                consumer.with_profile("plain", &typed),
            ));
            moved.push((
                format!("{}: --profile {profile} differs from plain", consumer.name),
                consumer.with_profile(profile, &[]),
                consumer.with_profile("plain", &[]),
            ));
        }
    }
    assert_all_equal(&run_pairs(home.path(), same));
    assert_all_differ(&run_pairs(home.path(), moved));
}

/// (2) A typed flag replaces only its own field; every other field still comes from the
/// profile. A typed value equal to the command's default (59 °F, 29.92 inHg, 50 %, 0 ft) still
/// wins over the stored one.
#[test]
fn a_typed_flag_overrides_the_stored_value_field_by_field() {
    let (home, list) = setup();
    let mut same = Vec::new();
    let mut moved = Vec::new();
    for consumer in &list {
        let wx = consumer.station_pairs(&WX);
        for (flag, default) in [
            ("--temperature", "59"),
            ("--pressure", "29.92"),
            ("--humidity", "50"),
            ("--altitude", "0"),
        ] {
            if flag == "--humidity" && !consumer.flies_humidity {
                continue;
            }
            let typed = vec![flag.to_string(), default.to_string()];
            same.push((
                format!("{}: --profile wx {flag} {default}", consumer.name),
                consumer.with_profile("wx", &typed),
                consumer.with_profile("plain", &with_fields(&wx, &[(flag, default)], &[])),
            ));
            if matches!(flag, "--temperature" | "--humidity") {
                moved.push((
                    format!("{}: a typed default {flag} wins", consumer.name),
                    consumer.with_profile("wx", &typed),
                    consumer.with_profile("wx", &[]),
                ));
            }
        }
        // A typed --pressure-type replaces the stored one.
        let qnh = consumer.station_pairs(&QNH);
        let absolute = vec!["--pressure-type".to_string(), "absolute".to_string()];
        same.push((
            format!("{}: --profile qnh --pressure-type absolute", consumer.name),
            consumer.with_profile("qnh", &absolute),
            consumer.with_profile(
                "plain",
                &with_fields(&qnh, &[("--pressure-type", "absolute")], &[]),
            ),
        ));
        moved.push((
            format!("{}: a typed --pressure-type wins", consumer.name),
            consumer.with_profile("qnh", &absolute),
            consumer.with_profile("qnh", &[]),
        ));
        // A typed --pressure is a different reading: the stored type, which describes the stored
        // pressure, does not apply to it.
        let typed_pressure = vec!["--pressure".to_string(), "30.10".to_string()];
        same.push((
            format!(
                "{}: --profile qnh --pressure 30.10 is absolute",
                consumer.name
            ),
            consumer.with_profile("qnh", &typed_pressure),
            consumer.with_profile("plain", &with_fields(&qnh, &[], &["--pressure-type"])),
        ));
    }
    assert_all_equal(&run_pairs(home.path(), same));
    assert_all_differ(&run_pairs(home.path(), moved));
}

/// The station triple a stored density altitude flies, in imperial display units, worked out
/// the way the CLI works it out (its `UnitConverter` factors), printed so it parses back to the
/// same bits.
fn density_altitude_station(da_ft: f64, temperature_f: f64) -> [(String, String); 3] {
    let explicit_c = (temperature_f - 32.0) * 5.0 / 9.0;
    let (altitude_m, temperature_c, pressure_hpa) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(
            da_ft * 0.3048,
            Some(explicit_c),
        );
    [
        (
            "--temperature".to_string(),
            format!("{:?}", temperature_c * 9.0 / 5.0 + 32.0),
        ),
        (
            "--pressure".to_string(),
            format!("{:?}", pressure_hpa / 33.8639),
        ),
        (
            "--altitude".to_string(),
            format!("{:?}", altitude_m / 0.3048),
        ),
    ]
}

/// A stored density altitude flies its derived station (both readings supplied), unless the run
/// types a pressure or an altitude of its own, which it would otherwise supersede. A typed
/// `--pressure-type` alone does not stop it.
#[test]
fn a_stored_density_altitude_flies_unless_the_run_types_its_own_pressure_or_altitude() {
    let (home, list) = setup();
    let station = density_altitude_station(7000.0, 80.0);
    let mut same = Vec::new();
    let mut moved = Vec::new();
    for consumer in &list {
        let derived: Vec<String> = if consumer.name == "trajectory" {
            // trajectory has its own --density-altitude flag.
            flags(&DA)
        } else {
            let mut typed: Vec<String> = station
                .iter()
                .flat_map(|(flag, value)| [flag.clone(), value.clone()])
                .collect();
            if consumer.flies_humidity {
                typed.extend(["--humidity".to_string(), "30".to_string()]);
            }
            typed
        };
        same.push((
            format!("{}: --profile da flies the derived station", consumer.name),
            consumer.with_profile("da", &[]),
            consumer.with_profile("plain", &derived),
        ));
        let stored_air = consumer.station_pairs(&[("--temperature", "80"), ("--humidity", "30")]);
        for blocker in [("--altitude", "0"), ("--pressure", "29.92")] {
            same.push((
                format!(
                    "{}: a typed {} stops the stored DA",
                    consumer.name, blocker.0
                ),
                consumer.with_profile("da", &[blocker.0.to_string(), blocker.1.to_string()]),
                consumer.with_profile("plain", &with_fields(&stored_air, &[blocker], &[])),
            ));
        }
        same.push((
            format!(
                "{}: --pressure-type does not stop the stored DA",
                consumer.name
            ),
            consumer.with_profile("da", &["--pressure-type".into(), "qnh".into()]),
            consumer.with_profile("da", &[]),
        ));
        moved.push((
            format!("{}: the stored DA moves the output", consumer.name),
            consumer.with_profile("da", &[]),
            consumer.with_profile("plain", &flags(&stored_air)),
        ));
    }
    assert_all_equal(&run_pairs(home.path(), same));
    assert_all_differ(&run_pairs(home.path(), moved));
}

/// (3) A profile saved with no weather flags flies exactly what the same bullet flies with no
/// profile, at sea level and at a typed altitude (where the stored, unmarked 59 °F and 29.92 inHg
/// must still fly the ICAO standard there). These runs were also byte-compared with the 0.50.0
/// pre-change build.
#[test]
fn a_profile_without_stored_weather_flies_what_no_profile_flies() {
    let (home, list) = setup();
    let mut same = Vec::new();
    for consumer in &list {
        let Some(bare) = &consumer.bare else { continue };
        for extra in [vec![], vec!["--altitude".to_string(), "5000".to_string()]] {
            let mut bare = bare.clone();
            bare.extend(extra.iter().cloned());
            same.push((
                format!("{}: --profile plain {extra:?}", consumer.name),
                consumer.with_profile("plain", &extra),
                bare,
            ));
        }
    }
    assert_all_equal(&run_pairs(home.path(), same));
}

/// A profile saved in one unit system with no weather flags, flown in the other, flies the run's
/// own default pressure: the 29.92 inHg it stored is the imperial default, not a reading, so a
/// metric run flies 1013.25 hPa, as it does with no profile (it flew 29.92 inHg converted,
/// 1013.2079 hPa, on the commands that honoured a stored pressure before).
#[test]
fn a_weatherless_profile_flies_the_runs_own_default_in_the_other_unit_system() {
    let (home, list) = setup();
    let mut same = Vec::new();
    for consumer in &list {
        let mut metric = vec!["-u".to_string(), "metric".to_string()];
        metric.extend(consumer.with_profile("plain", &[]));
        let mut typed = metric.clone();
        typed.extend(["--pressure".to_string(), "1013.25".to_string()]);
        same.push((
            format!("{}: -u metric --profile plain (imperial)", consumer.name),
            metric,
            typed,
        ));
    }
    assert_all_equal(&run_pairs(home.path(), same));
}

/// The stored weather converts to the run's units like every other stored field.
#[test]
fn the_stored_weather_flies_in_the_runs_units() {
    let (home, list) = setup();
    let mut same = Vec::new();
    let metric_wx = [
        (
            "--temperature",
            format!("{:?}", (95.0_f64 - 32.0) * 5.0 / 9.0),
        ),
        ("--pressure", format!("{:?}", 25.5_f64 * 33.8639)),
        ("--humidity", "90".to_string()),
        ("--altitude", format!("{:?}", 5000.0_f64 * 0.3048)),
    ];
    for consumer in list.iter().filter(|c| {
        matches!(
            c.name,
            "trajectory" | "come-ups" | "mark-to-range" | "hold-corridor" | "stability"
        )
    }) {
        let pairs: Vec<(&str, &str)> = metric_wx
            .iter()
            .map(|(flag, value)| (*flag, value.as_str()))
            .collect();
        let typed = flags(&consumer.station_pairs(&pairs));
        let metric = |args: Vec<String>| {
            let mut full = vec!["-u".to_string(), "metric".to_string()];
            full.extend(args);
            full
        };
        same.push((
            format!("{}: -u metric --profile wx", consumer.name),
            metric(consumer.with_profile("wx", &[])),
            metric(consumer.with_profile("plain", &typed)),
        ));
    }
    assert_all_equal(&run_pairs(home.path(), same));
}

/// (4) The zero day follows the same rule. `trajectory --auto-zero` inherits the shot day's
/// resolved station for every `--zero-*` flag left out, and a typed `--zero-*` replaces only its
/// own field; every other consumer zeroes in the air it flies, which (1) and (2) already hold
/// byte for byte, and `compare` reports each load's zero angle, held below.
#[test]
fn the_zero_day_follows_the_same_rule() {
    let (home, list) = setup();
    let trajectory = list.iter().find(|c| c.name == "trajectory").unwrap();
    let wx = flags(&WX);
    let mut same = Vec::new();
    for zero_day in [
        vec![],
        vec!["--zero-humidity", "10"],
        vec!["--zero-temperature", "59"],
        vec!["--zero-altitude", "0"],
        vec!["--zero-pressure", "29.92"],
    ] {
        let zero_day: Vec<String> = zero_day.iter().map(|s| s.to_string()).collect();
        let mut plain = wx.clone();
        plain.extend(zero_day.iter().cloned());
        same.push((
            format!("trajectory --saved-profile wx {zero_day:?}"),
            trajectory.with_profile("wx", &zero_day),
            trajectory.with_profile("plain", &plain),
        ));
        // A typed shot-day flag the zero day inherits.
        let mut typed_shot = zero_day.clone();
        typed_shot.extend(["--humidity".to_string(), "50".to_string()]);
        same.push((
            format!("trajectory --saved-profile wx --humidity 50 {zero_day:?}"),
            trajectory.with_profile("wx", &typed_shot),
            trajectory.with_profile(
                "plain",
                &[
                    with_fields(&WX, &[("--humidity", "50")], &[]),
                    zero_day.clone(),
                ]
                .concat(),
            ),
        ));
    }
    assert_all_equal(&run_pairs(home.path(), same));
}

fn compare_args(extra: &[&str]) -> Vec<String> {
    strings(&[
        &[
            "compare",
            "--zero-distance",
            "100",
            "--end",
            "600",
            "--load",
            "B:g1:0.500:180:2750:0.308",
            "-o",
            "json",
        ],
        extra,
    ])
}

fn compare_json(home: &Path, extra: &[&str]) -> Value {
    serde_json::from_str(&stdout(home, &compare_args(extra))).unwrap()
}

/// Everything `compare` reports about one load: its entry and its row in every range.
fn compare_load(doc: &Value, name: &str) -> (Value, Vec<Value>) {
    let mut entry = doc["compare"]["loads"]
        .as_array()
        .unwrap()
        .iter()
        .find(|load| load["name"] == name)
        .unwrap_or_else(|| panic!("no load {name}: {doc}"))
        .clone();
    entry["name"] = Value::Null;
    let rows = doc["compare"]["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            let mut load = row["loads"]
                .as_array()
                .unwrap()
                .iter()
                .find(|load| load["name"] == name)
                .unwrap()
                .clone();
            load["name"] = Value::Null;
            for delta in [
                "delta_drop",
                "delta_drift",
                "delta_velocity",
                "delta_energy",
            ] {
                load[delta] = Value::Null;
            }
            load
        })
        .collect();
    (entry, rows)
}

/// `compare` flies each profile load in its own stored weather, zero included, and a `--load`
/// spec in the typed or default air; a load that flies other air than the first load reports it.
#[test]
fn compare_flies_each_profile_load_in_its_stored_weather() {
    let (home, _) = setup();
    let home = home.path();
    let stored = compare_json(home, &["--profile", "wx"]);
    let typed: Vec<&str> = WX.iter().flat_map(|(f, v)| [*f, *v]).collect();
    let typed_run = compare_json(home, &[&["--profile", "plain"][..], &typed].concat());
    let (mut wx_entry, wx_rows) = compare_load(&stored, "wx");
    let (mut plain_entry, plain_rows) = compare_load(&typed_run, "plain");
    // The air block names the air the load flies: both fly WX; only the stored run's differs
    // from the first load's.
    let wx_air = wx_entry["air"].take();
    let plain_air = plain_entry["air"].take();
    assert_eq!(
        wx_entry, plain_entry,
        "the wx load's entry, zero angle included"
    );
    assert_eq!(wx_rows, plain_rows, "the wx load's rows");
    assert!(
        plain_air.is_null(),
        "every load flies the typed air: {typed_run}"
    );
    assert_eq!(
        wx_air["temperature"], 95.0,
        "the wx load reports the air it flies: {stored}"
    );
    assert_eq!(wx_air["humidity"], 90.0, "{stored}");
    assert_eq!(wx_air["altitude"], 5000.0, "{stored}");
    // The --load spec flies the default air either way.
    let default_run = compare_json(home, &["--profile", "plain"]);
    assert_eq!(compare_load(&stored, "B"), compare_load(&default_run, "B"));
    assert_eq!(
        stored["compare"]["conditions"],
        default_run["compare"]["conditions"]
    );
    // And says so.
    let output = run(home, &compare_args(&["--profile", "wx"]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("fly different air"),
        "compare says the loads fly different air:\n{stderr}"
    );
}

/// A profile with no stored weather is the same load as its `--load` spec in `compare`.
#[test]
fn compare_flies_a_weatherless_profile_as_its_load_spec() {
    let (home, _) = setup();
    let home = home.path();
    assert_eq!(
        stdout(home, &compare_args(&["--profile", "plain"])),
        stdout(
            home,
            &compare_args(&["--load", "plain:g1:0.475:168:2700:0.308"])
        )
    );
}

const BAND: &str = "relative humidity must be a percent from 0 to 100 (a sensor reading from -5 \
                    to 105 is read as 0 or 100)";

/// A stored humidity outside the sensor band is refused where it is flown, naming the profile
/// field and the file to edit; where a typed `--humidity` replaces it, or the command flies no
/// humidity at all, it is never read and the run goes ahead.
#[test]
fn a_stored_humidity_is_refused_only_where_it_is_flown() {
    let (home, list) = setup();
    edit_profile(&home, "plain", "h106", |p| p["humidity"] = 106.0.into());
    let path = home.profile_path("h106");
    let refusal = format!(
        "profile 'h106' field humidity = 106: {BAND}; edit the humidity stored in the profile \
         (the profile is stored at {})",
        path.display()
    );
    let mut same = Vec::new();
    for consumer in &list {
        let args = consumer.with_profile("h106", &[]);
        let output = run(home.path(), &args);
        if consumer.flies_humidity {
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(1), "{}: {stderr}", consumer.name);
            assert!(
                stderr
                    .lines()
                    .any(|line| line == format!("Error: {refusal}")),
                "{}: the refusal names the field, the band and the file:\n{stderr}",
                consumer.name
            );
            let typed = vec!["--humidity".to_string(), "50".to_string()];
            same.push((
                format!("{}: --profile h106 --humidity 50", consumer.name),
                consumer.with_profile("h106", &typed),
                consumer.with_profile("plain", &typed),
            ));
        } else {
            same.push((
                format!("{}: --profile h106 never reads the humidity", consumer.name),
                args,
                consumer.with_profile("plain", &[]),
            ));
        }
    }
    // Commands that read a profile and fly no air at all. `dial-plan` reads the optic.
    let mut optic = flags(&[]);
    optic.extend(["--elevation-click".to_string(), "0.1mil".to_string()]);
    save(&home, "optic", &optic);
    edit_profile(&home, "optic", "h106o", |p| p["humidity"] = 106.0.into());
    same.push((
        "dial-plan --profile h106o".to_string(),
        strings(&[&[
            "dial-plan",
            "--profile",
            "h106o",
            "--elevation",
            "5mil",
            "--windage",
            "1mil",
            "--range",
            "600",
        ]]),
        strings(&[&[
            "dial-plan",
            "--profile",
            "optic",
            "--elevation",
            "5mil",
            "--windage",
            "1mil",
            "--range",
            "600",
        ]]),
    ));
    let reticle = home.path().join("reticle.json");
    let reticle = reticle.to_str().unwrap();
    same.push((
        "reticle hold --profile h106 --drop-mil 2".to_string(),
        strings(&[&[
            "reticle",
            "hold",
            "--reticle-json",
            reticle,
            "--drop-mil",
            "2",
            "--profile",
            "h106",
        ]]),
        strings(&[&[
            "reticle",
            "hold",
            "--reticle-json",
            reticle,
            "--drop-mil",
            "2",
            "--profile",
            "plain",
        ]]),
    ));
    assert_all_equal(&run_pairs(home.path(), same));
}

/// Stored air that cannot exist is refused by every consumer that flies it, blamed on the
/// profile field that put it there, before any work; a typed reading that replaces the field
/// lets the run go ahead.
#[test]
fn stored_air_that_cannot_exist_is_refused_naming_the_profile_field() {
    let (home, list) = setup();
    edit_profile(&home, "plain", "boil", |p| {
        p["temperature"] = 212.0.into();
        p["temperature_supplied"] = true.into();
        p["humidity"] = 100.0.into();
    });
    edit_profile(&home, "plain", "hot", |p| {
        p["temperature"] = 250.0.into();
        p["temperature_supplied"] = true.into();
    });
    let mut same = Vec::new();
    for consumer in &list {
        let cases: &[(&str, &str, &[&str])] = if consumer.flies_humidity {
            &[
                (
                    "boil",
                    "profile 'boil' field humidity = 100 at 212 °F and 29.92 inHg: no such air \
                     exists",
                    &["--humidity", "50"],
                ),
                (
                    "hot",
                    "profile 'hot' field humidity = 50 at 250 °F and 29.92 inHg: no such air \
                     exists",
                    &["--temperature", "95"],
                ),
            ]
        } else {
            &[]
        };
        for (profile, cause, fix) in cases {
            let output = run(home.path(), &consumer.with_profile(profile, &[]));
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{} --profile {profile}: {stderr}",
                consumer.name
            );
            let first = stderr
                .lines()
                .find(|line| !line.starts_with("Loaded "))
                .unwrap_or_default();
            assert!(
                first.starts_with(&format!("Error: {cause}")),
                "{} --profile {profile}: the refusal leads with the stored field\n  expected \
                 prefix: Error: {cause}\n  got: {first}",
                consumer.name
            );
            let fix: Vec<String> = fix.iter().map(|s| s.to_string()).collect();
            let mut replaced = fix.clone();
            if *profile == "boil" {
                replaced.extend(["--temperature".to_string(), "212".to_string()]);
            }
            same.push((
                format!("{} --profile {profile} {fix:?}", consumer.name),
                consumer.with_profile(profile, &fix),
                consumer.with_profile("plain", &replaced),
            ));
        }
    }
    assert_all_equal(&run_pairs(home.path(), same));
}

/// Every consumer's `--help` says where an omitted station flag comes from.
#[test]
fn every_consumers_help_names_the_profile_as_the_fallback() {
    let (home, list) = setup();
    let mut missing = Vec::new();
    for consumer in &list {
        let mut args: Vec<String> = consumer
            .args
            .iter()
            .take_while(|arg| !arg.starts_with('-'))
            .cloned()
            .collect();
        args.push("--help".to_string());
        let help = stdout(home.path(), &args);
        let flattened = help.split_whitespace().collect::<Vec<_>>().join(" ");
        let fields: &[&str] = if consumer.flies_humidity {
            &["temperature", "pressure", "humidity", "altitude"]
        } else {
            &["temperature", "pressure", "altitude"]
        };
        for field in fields {
            let phrase = format!("the {}'s stored {field}", consumer.profile_flag);
            if !flattened.contains(&phrase) {
                missing.push(format!("{}: {phrase}", consumer.name));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "help text missing:\n{}",
        missing.join("\n")
    );
}
