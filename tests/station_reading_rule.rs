//! One rule for each station reading on every surface: pressure, altitude and density altitude.
//!
//! - **Pressure.** A station pressure flies when it is finite and above 0 and the air it
//!   describes can exist (the impossible-air rule `StationAir::new` applies). Nothing narrower:
//!   no CLI flag, `--location` column, saved profile, card, truing request, solve or
//!   `atmosphere.density_altitude` request refuses a pressure another surface flies.
//! - **Altitude.** A station altitude flies when it is finite. solve-json has no range of its
//!   own for it.
//! - **Density altitude.** A density altitude flies when it is finite and the station the FAA
//!   rule's inverse gives for it can exist. Beyond that it is refused, naming the reading that
//!   carried it (the flag, the zero day's flag, the `--location` column, the profile field).
//! - **`atmosphere.density_altitude` answers wherever the reading lands.** Its pressure
//!   altitude is the NWS formula inside -5 km to 11 km, bit for bit, and the standard
//!   atmosphere's own altitude (its ICAO layers) outside it.
//! - **Every refusal names the reading where it came from, in that surface's units**, for a
//!   humidity, a temperature and a pressure alike.

#![cfg(all(feature = "bridge", feature = "cli"))]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

use ballistics_engine::atmosphere::{
    nws_pressure_altitude_ft, resolve_atmosphere_for_density_altitude,
};
use ballistics_engine::solve_json::decode_solve_request_v1;
use ballistics_engine::{AtmosphericConditions, RelativeHumidity, StationAtmospherePresence};
use serde_json::{json, Value};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

// --- the readings every surface must treat alike ----------------------------------------------

/// One station: temperature °C, pressure hPa (`None`: omitted, the standard pressure at the
/// altitude), altitude m, humidity percent.
#[derive(Debug, Clone, Copy)]
struct Station {
    temperature_c: f64,
    pressure_hpa: Option<f64>,
    altitude_m: f64,
    humidity: f64,
}

impl Station {
    fn label(&self) -> String {
        format!(
            "{} C, {:?} hPa, {} m, {} %",
            self.temperature_c, self.pressure_hpa, self.altitude_m, self.humidity
        )
    }

    /// The rule: the station air a solve of these readings flies, built by the library.
    fn flies(&self) -> bool {
        AtmosphericConditions {
            temperature: self.temperature_c,
            pressure: self.pressure_hpa.unwrap_or(1013.25),
            humidity: self.humidity,
            altitude: self.altitude_m,
        }
        .station_air(StationAtmospherePresence::new(
            true,
            self.pressure_hpa.is_some(),
        ))
        .is_ok()
    }

    /// The solve-json v1 `atmosphere` object of these readings.
    fn atmosphere(&self) -> Value {
        let mut atmosphere = json!({
            "altitude_m": self.altitude_m,
            "temperature_k": self.temperature_c + 273.15,
            "relative_humidity": self.humidity / 100.0,
        });
        if let Some(hpa) = self.pressure_hpa {
            atmosphere["pressure_pa"] = json!(hpa * 100.0);
        }
        atmosphere
    }
}

/// The pressures the owner's rule must fly or refuse alike: a near vacuum (1 Pa), the
/// stratosphere (5000 Pa), Everest's summit and below it (27000 Pa), a deep mine (110000 Pa) and
/// twice sea level (200000 Pa), dry and at 50 %; and the altitudes -1000 m, 20000 m and 90000 m,
/// with the pressure omitted and with a sea-level pressure given. Then two pressures no air
/// has, 0 and -100 Pa.
fn stations() -> Vec<Station> {
    let mut stations = Vec::new();
    for humidity in [0.0, 50.0] {
        for pressure_pa in [1.0, 5000.0, 27000.0, 110_000.0, 200_000.0] {
            stations.push(Station {
                temperature_c: 15.0,
                pressure_hpa: Some(pressure_pa / 100.0),
                altitude_m: 0.0,
                humidity,
            });
        }
        for altitude_m in [-1000.0, 20_000.0, 90_000.0] {
            for pressure_hpa in [None, Some(1013.25)] {
                stations.push(Station {
                    temperature_c: 15.0,
                    pressure_hpa,
                    altitude_m,
                    humidity,
                });
            }
        }
    }
    for pressure_pa in [0.0, -100.0] {
        stations.push(Station {
            temperature_c: 15.0,
            pressure_hpa: Some(pressure_pa / 100.0),
            altitude_m: 0.0,
            humidity: 0.0,
        });
    }
    stations
}

// --- the wire surfaces ---------------------------------------------------------------------------

fn solve_request(atmosphere: Value) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01089, "diameter_m": 0.00782, "drag_model": "G1",
                       "ballistic_coefficient": 0.475},
        "rifle": {"muzzle_velocity_mps": 823.0},
        "shot": {"max_range_m": 300.0, "zero_distance_m": 100.0},
        "atmosphere": atmosphere,
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": "rk45", "time_step_s": 0.001},
        "effects": {"magnus": false, "coriolis": false, "enhanced_spin_drift": false},
        "sampling": {"interval_m": 100.0}
    })
}

/// A solve of `atmosphere`, as the wire reports it: `Ok(summary)` or `Err(error)`.
fn solve(atmosphere: Value) -> Result<Value, Value> {
    let decoded = decode_solve_request_v1(&solve_request(atmosphere).to_string())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())?;
    ballistics_engine::solve_v1(decoded)
        .map(|success| serde_json::to_value(success).unwrap()["summary"].clone())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())
}

/// `atmosphere.density_altitude` of `atmosphere`, as the wire reports it.
fn density_altitude(atmosphere: &Value) -> Result<Value, Value> {
    use ballistics_engine::atmosphere_service::{
        decode_density_altitude_request_v1, density_altitude_v1,
    };
    decode_density_altitude_request_v1(&json!({ "atmosphere": atmosphere }))
        .and_then(|request| density_altitude_v1(&request))
        .map(|response| serde_json::to_value(response).unwrap())
        .map_err(|envelope| serde_json::to_value(envelope).unwrap()["error"].clone())
}

fn bridge(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap()
}

/// A `card.range_table` request of `station` (°C, hPa, metres), the load in the card's
/// imperial preset.
fn card_request(station: &Station) -> Value {
    let mut request = json!({
        "muzzle_velocity": 2700.0, "ballistic_coefficient": 0.475, "mass": 168.0,
        "diameter": 0.308, "zero_distance": 100.0, "start": 100.0, "end": 300.0, "step": 100.0,
        "temperature": station.temperature_c, "temperature_unit": "celsius",
        "pressure_unit": "hpa", "altitude": station.altitude_m, "altitude_unit": "meters",
        "humidity": station.humidity
    });
    if let Some(hpa) = station.pressure_hpa {
        request["pressure"] = json!(hpa);
    }
    request
}

/// A truing model (`true.*`) of `station`, in the model's imperial units.
fn truing_model(station: &Station) -> Value {
    let mut model = json!({
        "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
        "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
        "zero_distance_yd": 100.0, "sight_height_in": 2.0,
        "temperature_f": station.temperature_c * 9.0 / 5.0 + 32.0,
        "humidity_pct": station.humidity, "altitude_ft": station.altitude_m / 0.3048
    });
    if let Some(hpa) = station.pressure_hpa {
        model["pressure_inhg"] = json!(hpa / 33.8639);
    }
    model
}

fn wind_truing_request(station: &Station) -> Value {
    let mut request = truing_model(station);
    let object = request.as_object_mut().unwrap();
    let bc = object.remove("ballistic_coefficient").unwrap();
    object.insert("bc".into(), bc);
    object.insert(
        "observations".into(),
        json!([{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}]),
    );
    object.insert("twist".into(), json!({"rate_in": 11.0, "right_hand": true}));
    object.insert("earth".into(), Value::Null);
    object.insert("called_crosswind_mph".into(), Value::Null);
    request
}

// --- the CLI ------------------------------------------------------------------------------------

/// A private `HOME`, removed when dropped, so no saved profile of the machine's user can reach a
/// run.
struct Home(PathBuf);

impl Drop for Home {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

impl Home {
    fn new() -> Self {
        static N: AtomicU32 = AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "bx-station-reading-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self) -> &Path {
        &self.0
    }

    fn file(&self, name: &str, contents: &str) -> String {
        let path = self.0.join(name);
        std::fs::write(&path, contents).unwrap();
        path.to_string_lossy().into_owned()
    }

    fn profile_path(&self, name: &str) -> PathBuf {
        self.0
            .join(".ballistics")
            .join("profiles")
            .join(format!("{name}.json"))
    }

    /// Save a profile with `profile save`, then overwrite its stored fields with `edits`, as a
    /// hand-edited file or another tool would leave it.
    fn profile(&self, name: &str, units: &str, flags: &[&str], edits: &[(&str, Value)]) {
        let mut args = vec!["-u", units, "profile", "save", name];
        args.extend_from_slice(match units {
            "metric" => &[
                "-v",
                "823",
                "-b",
                "0.475",
                "-m",
                "10.9",
                "-d",
                "7.82",
                "--zero-distance",
                "91.44",
                "--sight-height",
                "50",
                "--twist-rate",
                "10",
                "--bullet-length",
                "31.5",
            ],
            _ => &[
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
                "--sight-height",
                "2",
                "--twist-rate",
                "10",
                "--bullet-length",
                "1.24",
            ],
        });
        args.extend_from_slice(flags);
        let output = run(self.path(), &args);
        assert!(
            output.status.success(),
            "profile save {flags:?}: {}",
            stderr(&output)
        );
        if !edits.is_empty() {
            let path = self.profile_path(name);
            let mut stored: Value =
                serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
            for (field, value) in edits {
                stored[*field] = value.clone();
            }
            std::fs::write(&path, serde_json::to_string_pretty(&stored).unwrap()).unwrap();
        }
    }
}

fn run(home: &Path, args: &[&str]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The first `Error: ` line of a run.
fn error_line(output: &Output) -> String {
    stderr(output)
        .lines()
        .find(|line| line.starts_with("Error: "))
        .unwrap_or_default()
        .to_string()
}

/// The station flags of `station` for a metric run: every reading typed, the pressure only when
/// given.
fn metric_flags(
    station: &Station,
    temperature: &str,
    pressure: &str,
    altitude: &str,
) -> Vec<String> {
    let mut flags = vec![
        format!("{temperature}={}", station.temperature_c),
        format!("{altitude}={}", station.altitude_m),
    ];
    if let Some(hpa) = station.pressure_hpa {
        flags.push(format!("{pressure}={hpa}"));
    }
    flags
}

/// Run `cases` in parallel: each a label and the argument list, every run in its own `HOME`.
fn run_all(cases: Vec<(String, Vec<String>)>) -> Vec<(String, Output)> {
    std::thread::scope(|scope| {
        let jobs: Vec<_> = cases
            .into_iter()
            .map(|(label, args)| {
                scope.spawn(move || {
                    let home = Home::new();
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    (label, run(home.path(), &args))
                })
            })
            .collect();
        jobs.into_iter().map(|job| job.join().unwrap()).collect()
    })
}

const METRIC_LOAD: [&str; 10] = [
    "-u", "metric", "-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82",
];

/// The CLI runs of `station` that must fly or refuse it as the rule does: `trajectory`, `zero`,
/// a calculator command (`come-ups`), `estimate-bc`, `true-velocity`, `trajectory --auto-zero`'s
/// zero day and a `--location` CSV row. Each is `(name, args)`; the CSV row is written into
/// `home`.
fn cli_runs(station: &Station, home: &Home, index: usize) -> Vec<(String, Vec<String>)> {
    let shot = metric_flags(station, "--temperature", "--pressure", "--altitude");
    let zero_day = metric_flags(
        station,
        "--zero-temperature",
        "--zero-pressure",
        "--zero-altitude",
    );
    let load: Vec<String> = METRIC_LOAD.iter().map(|s| s.to_string()).collect();
    let humidity = format!("--humidity={}", station.humidity);
    let zero_humidity = format!("--zero-humidity={}", station.humidity);
    let with = |head: &[&str], tail: &[String]| -> Vec<String> {
        let mut args: Vec<String> = Vec::new();
        args.extend(load.iter().take(2).cloned());
        args.extend(head.iter().map(|s| s.to_string()));
        args.extend(load.iter().skip(2).cloned());
        args.extend(tail.iter().cloned());
        args
    };
    let mut shot_air = shot.clone();
    shot_air.push(humidity.clone());
    let mut zero_air = zero_day.clone();
    zero_air.push(zero_humidity);
    let csv = home.file(
        &format!("loc-{index}.csv"),
        &format!(
            "NAME,TEMP,PRESSURE,ALT,HUMIDITY\nSITE,{},{},{},{}\n",
            station.temperature_c,
            station
                .pressure_hpa
                .map(|hpa| hpa.to_string())
                .unwrap_or_default(),
            station.altitude_m,
            station.humidity
        ),
    );
    vec![
        (
            "trajectory".into(),
            with(&["trajectory", "--max-range", "300"], &shot_air),
        ),
        (
            "zero".into(),
            with(&["zero", "--target-distance", "200"], &shot_air),
        ),
        (
            "come-ups".into(),
            with(
                &[
                    "come-ups",
                    "--zero-distance",
                    "100",
                    "--end",
                    "300",
                    "--sight-height",
                    "50",
                ],
                &shot_air,
            ),
        ),
        (
            "estimate-bc".into(),
            with(
                &["estimate-bc", "--data", "100,0.02;200,0.15;300,0.45"],
                &shot_air,
            )
            .into_iter()
            .filter(|arg| arg != "-b" && arg != "0.475")
            .collect(),
        ),
        (
            "trajectory --auto-zero (zero day)".into(),
            with(
                &["trajectory", "--max-range", "300", "--auto-zero", "100"],
                &zero_air,
            ),
        ),
        (
            "trajectory --location".into(),
            with(
                &[
                    "trajectory",
                    "--max-range",
                    "300",
                    "--location",
                    &csv,
                    "--site",
                    "SITE",
                ],
                &[],
            ),
        ),
    ]
}

// --- one verdict on every surface ---------------------------------------------------------------

/// Each reading is flown, or refused, alike by the library's station, `solve`,
/// `atmosphere.density_altitude`, `card.range_table`, both truing models and every CLI path. In
/// 0.49.0 the CLI's `--pressure` refused anything outside 8..33 inHg (250..1100 hPa) that a
/// `--location` cell, `--zero-pressure` or `estimate-bc` flew; `atmosphere.density_altitude`
/// refused a station pressure whose pressure altitude was outside -5 km to 11 km; and solve-json
/// refused an altitude outside -5000..84000 m. A surface that refuses names the reading.
#[test]
fn every_surface_flies_or_refuses_each_station_reading_alike() {
    let home = Home::new();
    let stations = stations();
    let mut cli_cases = Vec::new();
    for (index, station) in stations.iter().enumerate() {
        for (name, args) in cli_runs(station, &home, index) {
            cli_cases.push((format!("{index} {name}"), args));
        }
    }
    let cli = run_all(cli_cases);
    let mut flown = 0;
    let mut failures = Vec::new();
    for (index, station) in stations.iter().enumerate() {
        let label = station.label();
        let rule = station.flies();
        flown += usize::from(rule);
        let atmosphere = station.atmosphere();
        let solved = solve(atmosphere.clone());
        if solved.is_ok() != rule {
            failures.push(format!("{label}: solve {solved:?}"));
        }
        let answered = density_altitude(&atmosphere);
        match (&answered, rule) {
            (Ok(answer), true) => {
                for figure in [
                    &answer["pressure_altitude"],
                    &answer["density_altitude"]["faa_rule"],
                    &answer["density_altitude"]["density_matched"],
                ] {
                    for unit in ["m", "ft"] {
                        if !figure[unit].as_f64().is_some_and(f64::is_finite) {
                            failures.push(format!("{label}: density altitude {answer}"));
                        }
                    }
                }
            }
            (Err(_), false) => {}
            _ => failures.push(format!("{label}: atmosphere.density_altitude {answered:?}")),
        }
        let carded = bridge("card.range_table", card_request(station));
        if (carded["ok"] == true) != rule {
            failures.push(format!("{label}: card.range_table {carded}"));
        }
        let model: Result<ballistics_engine::truing::TruingModelInputsV1, _> =
            serde_json::from_value(truing_model(station));
        let checked = model.map_err(|e| e.to_string()).and_then(|m| m.validate());
        if checked.is_ok() != rule {
            failures.push(format!("{label}: true.* model {checked:?}"));
        }
        let wind: Result<ballistics_engine::truing_wind::WindTruingRequest, _> =
            serde_json::from_value(wind_truing_request(station));
        let checked = wind.map_err(|e| e.to_string()).and_then(|w| w.validate());
        if checked.is_ok() != rule {
            failures.push(format!("{label}: true.wind {checked:?}"));
        }
        for (name, output) in cli
            .iter()
            .filter(|(name, _)| name.split(' ').next() == Some(&index.to_string()))
        {
            if output.status.success() != rule {
                failures.push(format!(
                    "{label}: {name}: exit {:?}: {}",
                    output.status.code(),
                    error_line(output)
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} disagreements with the rule:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(
        flown > 10 && flown < stations.len() - 2,
        "{flown} of {}",
        stations.len()
    );
}

/// The C ABI flies, or refuses, the same readings: a non-null trajectory exactly when the rule
/// flies the station.
#[cfg(feature = "ffi")]
#[test]
fn the_c_abi_flies_or_refuses_each_station_reading_alike() {
    use ballistics_engine::ffi::{
        ballistics_calculate_trajectory_with_station_presence, ballistics_free_trajectory_result,
        FFIAtmosphericConditions, FFIBallisticInputs, FFIWindConditions,
        BALLISTICS_STATION_PRESSURE_SUPPLIED, BALLISTICS_STATION_TEMPERATURE_SUPPLIED,
    };
    let inputs = FFIBallisticInputs {
        muzzle_velocity: 823.0,
        muzzle_angle: 0.0,
        bc_value: 0.475,
        bullet_mass: 0.01089,
        bullet_diameter: 0.00782,
        bc_type: 0,
        sight_height: 0.05,
        target_distance: 300.0,
        temperature: 15.0,
        twist_rate: 10.0,
        is_twist_right: 1,
        shooting_angle: 0.0,
        altitude: 0.0,
        latitude: f64::NAN,
        azimuth_angle: 0.0,
        use_rk4: 1,
        use_adaptive_rk45: 0,
        enable_wind_shear: 0,
        enable_trajectory_sampling: 0,
        sample_interval: 0.0,
        enable_pitch_damping: 0,
        enable_precession_nutation: 0,
        enable_spin_drift: 0,
        enable_magnus: 0,
        enable_coriolis: 0,
        shot_azimuth: 0.0,
        cant_angle: 0.0,
        zero_poi_vertical: 0.0,
        zero_poi_horizontal: 0.0,
        sight_offset_lateral: 0.0,
    };
    let wind = FFIWindConditions {
        speed: 0.0,
        direction: 0.0,
        vertical_speed: 0.0,
    };
    for station in stations() {
        let atmosphere = FFIAtmosphericConditions {
            temperature: station.temperature_c,
            pressure: station.pressure_hpa.unwrap_or(1013.25),
            humidity: station.humidity,
            altitude: station.altitude_m,
        };
        let bits = BALLISTICS_STATION_TEMPERATURE_SUPPLIED
            | if station.pressure_hpa.is_some() {
                BALLISTICS_STATION_PRESSURE_SUPPLIED
            } else {
                0
            };
        let inputs = FFIBallisticInputs {
            altitude: station.altitude_m,
            ..inputs
        };
        let result = unsafe {
            ballistics_calculate_trajectory_with_station_presence(
                &inputs,
                &wind,
                &atmosphere,
                300.0,
                1.0,
                std::ptr::null(),
                std::ptr::null(),
                0,
                1.0,
                bits,
            )
        };
        let flew = !result.is_null();
        unsafe { ballistics_free_trajectory_result(result) };
        assert_eq!(flew, station.flies(), "{}", station.label());
    }
}

// --- density altitude -----------------------------------------------------------------------------

/// The density altitudes the rule must fly or refuse alike, m, each without a temperature (the
/// FAA rule's standard temperature) and with 40 °C: far below sea level, the stratosphere, and
/// either side of where the rule's inverse leaves no air (its standard temperature reaches 0 K
/// at a pressure altitude of 145286 ft, and the NWS pressure formula 0 hPa at 145366 ft).
fn density_altitudes() -> Vec<(f64, Option<f64>)> {
    let mut cases = Vec::new();
    for temperature_c in [None, Some(40.0)] {
        for density_altitude_m in [-10_000.0, 20_000.0, 40_000.0, 44_000.0, 44_500.0, 50_000.0] {
            cases.push((density_altitude_m, temperature_c));
        }
    }
    cases.push((60_000.0, Some(40.0)));
    cases
}

/// The rule for a density altitude: the station the FAA rule's inverse gives for it exists,
/// dry.
fn density_altitude_flies(density_altitude_m: f64, temperature_c: Option<f64>) -> bool {
    let (altitude_m, temperature_c, pressure_hpa) =
        resolve_atmosphere_for_density_altitude(density_altitude_m, temperature_c);
    ballistics_engine::StationAir::new(
        altitude_m,
        temperature_c,
        pressure_hpa,
        RelativeHumidity::DRY,
    )
    .is_ok()
}

/// Every surface that takes a density altitude flies it exactly when its station exists, and
/// refuses it naming the reading that carried it: `--density-altitude`, `--zero-density-altitude`,
/// the `--location` column `DA`, a saved profile's `density_altitude`, and `zero`'s flag. In 0.49.0
/// the flags refused anything outside -5000..120000 (feet or metres), which a CSV cell or a
/// profile flew, and a density altitude past the inverse's reach failed as "atmosphere.pressure
/// must be finite", naming a field nobody typed.
#[test]
fn every_surface_flies_or_refuses_each_density_altitude_alike() {
    let home = Home::new();
    let mut cases = Vec::new();
    for (index, &(density_altitude_m, temperature_c)) in density_altitudes().iter().enumerate() {
        let da = density_altitude_m.to_string();
        let temperature: Vec<String> = temperature_c
            .map(|t| vec![format!("--temperature={t}")])
            .unwrap_or_default();
        let zero_temperature: Vec<String> = temperature_c
            .map(|t| vec![format!("--zero-temperature={t}")])
            .unwrap_or_default();
        let head = |command: &[&str]| -> Vec<String> {
            let mut args: Vec<String> = vec!["-u".into(), "metric".into()];
            args.extend(command.iter().map(|s| s.to_string()));
            args.extend(METRIC_LOAD.iter().skip(2).map(|s| s.to_string()));
            args.push("--humidity=0".into());
            args
        };
        let csv = home.file(
            &format!("da-{index}.csv"),
            &format!(
                "NAME,TEMP,DA,HUMIDITY\nSITE,{},{da},0\n",
                temperature_c.map(|t| t.to_string()).unwrap_or_default()
            ),
        );
        let profile = format!("da{index}");
        home.profile(
            &profile,
            "metric",
            &[],
            &[
                ("density_altitude", json!(density_altitude_m)),
                ("humidity", json!(0.0)),
                ("temperature", json!(temperature_c.unwrap_or(15.0))),
                ("temperature_supplied", json!(temperature_c.is_some())),
            ],
        );
        let label_da = |label: &str| {
            format!("{label} = {density_altitude_m} m: no air has this density altitude")
        };
        let mut trajectory = head(&["trajectory", "--max-range", "300"]);
        trajectory.extend(temperature.iter().cloned());
        trajectory.push(format!("--density-altitude={da}"));
        let mut zero_day = head(&["trajectory", "--max-range", "300", "--auto-zero", "100"]);
        zero_day.extend(zero_temperature.iter().cloned());
        zero_day.push(format!("--zero-density-altitude={da}"));
        let mut zero = head(&["zero", "--target-distance", "200"]);
        zero.extend(temperature.iter().cloned());
        zero.push(format!("--density-altitude={da}"));
        let mut located = head(&["trajectory", "--max-range", "300"]);
        located.extend([
            "--location".into(),
            csv.clone(),
            "--site".into(),
            "SITE".into(),
        ]);
        let mut stored: Vec<String> = vec!["-u".into(), "metric".into(), "come-ups".into()];
        stored.extend([
            "--profile".into(),
            profile.clone(),
            "--zero-distance".into(),
            "100".into(),
            "--end".into(),
            "300".into(),
        ]);
        for (name, args, label) in [
            ("trajectory", trajectory, label_da("--density-altitude")),
            ("zero day", zero_day, label_da("--zero-density-altitude")),
            ("zero", zero, label_da("--density-altitude")),
            ("--location", located, label_da("--location CSV column DA")),
            (
                "profile",
                stored,
                label_da(&format!("profile '{profile}' field density_altitude")),
            ),
        ] {
            cases.push((index, density_altitude_m, temperature_c, name, args, label));
        }
    }
    let outputs = std::thread::scope(|scope| {
        let jobs: Vec<_> = cases
            .iter()
            .map(|(_, _, _, _, args, _)| {
                let home = home.path();
                scope.spawn(move || {
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    run(home, &args)
                })
            })
            .collect();
        jobs.into_iter()
            .map(|job| job.join().unwrap())
            .collect::<Vec<_>>()
    });
    let mut failures = Vec::new();
    let mut refused = 0;
    for ((_, density_altitude_m, temperature_c, name, _, label), output) in
        cases.iter().zip(&outputs)
    {
        let rule = density_altitude_flies(*density_altitude_m, *temperature_c);
        let error = error_line(output);
        if output.status.success() != rule {
            failures.push(format!(
                "{density_altitude_m} m, {temperature_c:?} C, {name}: exit {:?}: {error}",
                output.status.code()
            ));
        } else if !rule {
            refused += 1;
            if !error.starts_with(&format!("Error: {label}")) {
                failures.push(format!(
                    "{density_altitude_m} m, {temperature_c:?} C, {name}: expected `{label}`, \
                     got `{error}`"
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} disagreements:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(refused >= 10, "{refused} refusals");
}

// --- the pressure altitude beyond the troposphere --------------------------------------------------

/// The ICAO standard atmosphere at a geopotential height, m, computed here independently of the
/// engine from its published layers: `(T K, P Pa)`.
fn standard_air(height_m: f64) -> (f64, f64) {
    const LAYERS: [(f64, f64); 7] = [
        (0.0, -0.0065),
        (11_000.0, 0.0),
        (20_000.0, 0.001),
        (32_000.0, 0.0028),
        (47_000.0, 0.0),
        (51_000.0, -0.0028),
        (71_000.0, -0.002),
    ];
    let (g, r) = (9.80665, 287.0531);
    let (mut temperature, mut pressure) = (288.15, 101_325.0);
    for (index, &(base, lapse)) in LAYERS.iter().enumerate() {
        let top = LAYERS.get(index + 1).map_or(f64::INFINITY, |next| next.0);
        let rise = height_m.min(top) - base;
        let next = temperature + lapse * rise;
        pressure *= if lapse == 0.0 {
            (-g * rise / (r * temperature)).exp()
        } else {
            (next / temperature).powf(-g / (lapse * r))
        };
        temperature = next;
        if height_m <= top {
            break;
        }
    }
    (temperature, pressure)
}

fn pressure_altitude_m(pressure_pa: f64) -> f64 {
    let answer = density_altitude(&json!({"temperature_k": 216.65, "pressure_pa": pressure_pa,
                                          "relative_humidity": 0.0}))
    .unwrap_or_else(|error| panic!("{pressure_pa} Pa is answered: {error}"));
    answer["pressure_altitude"]["m"].as_f64().unwrap()
}

/// `atmosphere.density_altitude`'s pressure altitude is the standard atmosphere's own altitude
/// of the station pressure wherever it lands: standard pressure at every height up to 84 km
/// reads that height, within the few metres the NWS formula and the ICAO troposphere differ by
/// inside the band, and within a millimetre outside it. 0.49.0 refused every pressure there.
#[test]
fn the_pressure_altitude_reads_standard_pressure_at_its_own_height_in_every_layer() {
    for (height_m, tolerance_m) in [
        (-4000.0, 7.0),
        (0.0, 1e-6),
        (5000.0, 7.0),
        (10_990.0, 7.0),
        (11_010.0, 7.0),
        (15_000.0, 1e-3),
        (20_000.0, 1e-3),
        (25_000.0, 1e-3),
        (32_000.0, 1e-3),
        (40_000.0, 1e-3),
        (47_000.0, 1e-3),
        (49_000.0, 1e-3),
        (51_000.0, 1e-3),
        (60_000.0, 1e-3),
        (71_000.0, 1e-3),
        (80_000.0, 1e-3),
        (84_000.0, 1e-3),
    ] {
        let (_, pressure_pa) = standard_air(height_m);
        let altitude = pressure_altitude_m(pressure_pa);
        assert!(
            (altitude - height_m).abs() < tolerance_m,
            "standard pressure at {height_m} m ({pressure_pa} Pa) reads {altitude} m"
        );
    }
}

/// Inside -5 km to 11 km it is the NWS formula, to the bit, as every 0.49.0 answer was; and it
/// rises steadily as the pressure falls, across that band's edges and every layer boundary, down
/// to 1 Pa.
#[test]
fn the_pressure_altitude_is_the_nws_formula_inside_the_band_and_rises_steadily_beyond_it() {
    let mut inside = 0;
    let mut last = f64::NEG_INFINITY;
    for step in 0..=800 {
        let pressure_pa = 300_000.0 * 10f64.powf(-f64::from(step) / 145.0);
        let answer = density_altitude(&json!({"temperature_k": 250.0, "pressure_pa": pressure_pa,
                                              "relative_humidity": 0.0}))
        .unwrap_or_else(|error| panic!("{pressure_pa} Pa is answered: {error}"));
        let feet = answer["pressure_altitude"]["ft"].as_f64().unwrap();
        let metres = answer["pressure_altitude"]["m"].as_f64().unwrap();
        let nws = nws_pressure_altitude_ft(pressure_pa / 100.0);
        if (-5000.0..=11_000.0).contains(&(nws * 0.3048)) {
            inside += 1;
            assert_eq!(feet.to_bits(), nws.to_bits(), "{pressure_pa} Pa");
        }
        assert!(
            metres.is_finite() && metres >= last,
            "{pressure_pa} Pa reads {metres} m, below {last} m"
        );
        last = metres;
    }
    assert!(inside > 100, "{inside}");
    assert!(last > 70_000.0, "{last}");
}

// --- every refusal names its reading where it came from ---------------------------------------------

/// A typed pressure is held to the one rule on every command that takes one, in the flag's own
/// words and unit, and every reading the rule flies is flown: 5 and 40 inHg used to be refused as
/// "out of range (expected ~8..33 inHg for imperial units)".
#[test]
fn every_pressure_flag_applies_the_one_pressure_rule() {
    let home = Home::new();
    home.profile("p", "imperial", &[], &[]);
    let imperial_load = ["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];
    let commands: Vec<(&str, Vec<&str>)> = vec![
        ("--pressure", vec!["trajectory", "--max-range", "300"]),
        (
            "--zero-pressure",
            vec!["trajectory", "--max-range", "300", "--auto-zero", "100"],
        ),
        ("--pressure", vec!["zero", "--target-distance", "200"]),
        (
            "--pressure",
            vec!["estimate-bc", "--data", "100,0.5;200,2.0;300,5.0"],
        ),
        (
            "--pressure",
            vec!["come-ups", "--zero-distance", "100", "--end", "300"],
        ),
        (
            "--pressure",
            vec!["range-table", "--zero-distance", "100", "--end", "300"],
        ),
        ("--pressure", vec!["mpbr", "--vital-zone", "8"]),
        (
            "--pressure",
            vec!["stability", "--twist-rate", "10", "--length", "1.24"],
        ),
        ("--pressure", vec!["profile", "save", "q"]),
    ];
    for (flag, command) in &commands {
        let mut args: Vec<&str> = command.clone();
        if !command.contains(&"estimate-bc") && !command.contains(&"stability") {
            args.extend_from_slice(&imperial_load);
        } else if command.contains(&"estimate-bc") {
            args.extend_from_slice(&["-v", "2700", "-m", "168", "-d", "0.308"]);
        } else {
            args.extend_from_slice(&["-m", "168", "-d", "0.308"]);
        }
        for good in ["5", "40"] {
            let flag_value = format!("{flag}={good}");
            let mut flown = args.clone();
            flown.extend_from_slice(&[&flag_value, "--humidity", "0"]);
            if command.contains(&"stability") || command.contains(&"profile") {
                flown.retain(|arg| *arg != "--humidity" && *arg != "0");
            }
            let output = run(home.path(), &flown);
            assert!(
                output.status.success(),
                "{command:?} {flag} {good}: {}",
                stderr(&output)
            );
        }
        for (bad, refusal) in [
            ("-1", format!("Error: {flag} -1 inHg is at or below zero")),
            ("0", format!("Error: {flag} 0 inHg is at or below zero")),
            (
                "NaN",
                format!("Error: {flag} NaN inHg is not a finite number"),
            ),
        ] {
            let flag_value = format!("{flag}={bad}");
            let mut refused = args.clone();
            refused.push(&flag_value);
            let output = run(home.path(), &refused);
            assert_eq!(output.status.code(), Some(1), "{command:?} {flag} {bad}");
            assert_eq!(
                error_line(&output),
                refusal,
                "{command:?}:\n{}",
                stderr(&output)
            );
        }
    }
}

/// A typed altitude or density altitude is refused on its own only when it is not a finite
/// number, in the flag's own words and unit, on every command that takes one: in 0.49.0
/// `--altitude NaN` reached the solver as `atmosphere.altitude must be finite`, and
/// `--density-altitude NaN` as clap's range message.
#[test]
fn every_altitude_flag_refuses_only_a_number_that_is_not_finite() {
    let home = Home::new();
    let load = ["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];
    let commands: Vec<(&str, Vec<&str>)> = vec![
        ("--altitude", vec!["trajectory", "--max-range", "300"]),
        (
            "--zero-altitude",
            vec!["trajectory", "--max-range", "300", "--auto-zero", "100"],
        ),
        (
            "--density-altitude",
            vec!["trajectory", "--max-range", "300"],
        ),
        (
            "--zero-density-altitude",
            vec!["trajectory", "--max-range", "300", "--auto-zero", "100"],
        ),
        ("--altitude", vec!["zero", "--target-distance", "200"]),
        (
            "--density-altitude",
            vec!["zero", "--target-distance", "200"],
        ),
        (
            "--altitude",
            vec!["come-ups", "--zero-distance", "100", "--end", "300"],
        ),
        ("--altitude", vec!["profile", "save", "q"]),
        ("--density-altitude", vec!["profile", "save", "q"]),
    ];
    for (flag, command) in &commands {
        let mut args = command.clone();
        args.extend_from_slice(&load);
        for bad in ["NaN", "inf"] {
            let flag_value = format!("{flag}={bad}");
            let mut refused = args.clone();
            refused.push(&flag_value);
            let output = run(home.path(), &refused);
            assert_eq!(output.status.code(), Some(1), "{command:?} {flag} {bad}");
            let shown = if bad == "NaN" { "NaN" } else { "inf" };
            assert_eq!(
                error_line(&output),
                format!("Error: {flag} {shown} ft is not a finite number"),
                "{command:?}:\n{}",
                stderr(&output)
            );
        }
    }
}

/// A `--location` cell is named by its column, in the run's units, for a temperature and a
/// pressure as for a humidity: in 0.49.0 these reached the solver and came back as
/// `atmosphere.pressure = -33.86 hPa: pressure must be greater than 0`.
#[test]
fn a_location_cell_is_refused_naming_its_column() {
    let home = Home::new();
    let csv = home.file(
        "loc.csv",
        "NAME,TEMP,PRESSURE,ALT,HUMIDITY\nP,59,-1,0,50\nT,-500,29.92,0,50\nZ,59,0,0,50\n",
    );
    for (site, refusal) in [
        (
            "P",
            "Error: --location CSV column PRESSURE = -1 inHg: pressure must be greater than 0",
        ),
        (
            "Z",
            "Error: --location CSV column PRESSURE = 0 inHg: pressure must be greater than 0",
        ),
        (
            "T",
            "Error: --location CSV column TEMP = -500 °F: temperature must be above absolute \
             zero (-459.67 °F)",
        ),
    ] {
        let output = run(
            home.path(),
            &[
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
                "--location",
                &csv,
                "--site",
                site,
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{site}");
        assert_eq!(error_line(&output), refusal, "{site}:\n{}", stderr(&output));
    }
}

/// Every command that takes a saved profile names a stored reading it refuses by the profile's
/// field, in the profile's units, for a pressure, a temperature and a density altitude alike. In
/// 0.49.0 only a humidity or air that cannot exist was named; a stored pressure of -1 came back
/// from the solver as `atmosphere.pressure = -33.86 hPa: ...`, and `stability` printed a NaN.
#[test]
fn every_command_names_a_refused_stored_reading_by_its_profile_field() {
    let home = Home::new();
    home.profile(
        "p",
        "imperial",
        &[],
        &[
            ("pressure", json!(-1.0)),
            ("pressure_supplied", json!(true)),
        ],
    );
    home.profile(
        "t",
        "imperial",
        &[],
        &[
            ("temperature", json!(-500.0)),
            ("temperature_supplied", json!(true)),
        ],
    );
    home.profile(
        "d",
        "imperial",
        &[],
        &[("density_altitude", json!(200_000.0))],
    );
    let consumers: Vec<(&str, Vec<&str>)> = vec![
        ("--saved-profile", vec!["trajectory", "--max-range", "300"]),
        (
            "--profile",
            vec!["come-ups", "--zero-distance", "100", "--end", "300"],
        ),
        (
            "--profile",
            vec!["range-table", "--zero-distance", "100", "--end", "300"],
        ),
        (
            "--profile",
            vec!["wind-card", "--zero-distance", "100", "--end", "300"],
        ),
        ("--profile", vec!["mpbr", "--vital-zone", "8"]),
        (
            "--profile",
            vec!["lead", "--target-speed", "4", "--end", "300"],
        ),
        ("--profile", vec!["mark-to-range", "--mark", "2.0"]),
        (
            "--profile",
            vec!["optimal-zero", "--target", "200:10", "--target", "400:12"],
        ),
        (
            "--profile",
            vec![
                "adaptive-card",
                "--zero-distance",
                "100",
                "--start",
                "100",
                "--end",
                "300",
            ],
        ),
        (
            "--profile",
            vec![
                "plan-truing",
                "--measurement-resolution",
                "0.25",
                "--candidate-ranges",
                "300,500",
            ],
        ),
        (
            "--profile",
            vec!["stability", "--twist-rate", "10", "--length", "1.24"],
        ),
    ];
    for (flag, command) in &consumers {
        for (profile, refusal) in [
            (
                "p",
                "Error: profile 'p' field pressure = -1 inHg: pressure must be greater than 0",
            ),
            (
                "t",
                "Error: profile 't' field temperature = -500 °F: temperature must be above \
                 absolute zero (-459.67 °F)",
            ),
            (
                "d",
                "Error: profile 'd' field density_altitude = 200000 ft: no air has this density \
                 altitude",
            ),
        ] {
            let mut args = command.clone();
            args.extend_from_slice(&[flag, profile]);
            let output = run(home.path(), &args);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{command:?} {profile}:\n{}",
                stderr(&output)
            );
            let error = error_line(&output);
            assert!(
                error.starts_with(refusal),
                "{command:?} {profile}: expected `{refusal}`, got `{error}`"
            );
        }
    }
}

/// The zero day names its own readings by its own flags, a density altitude among them, and
/// one it inherits by the shot day's name.
#[test]
fn the_zero_day_names_each_refused_reading_by_its_own_flag() {
    let home = Home::new();
    let base = [
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
    ];
    for (extra, refusal) in [
        (
            vec!["--zero-pressure=-1"],
            "Error: --zero-pressure -1 inHg is at or below zero".to_string(),
        ),
        (
            vec!["--zero-density-altitude=200000"],
            "Error: --zero-density-altitude = 200000 ft: no air has this density altitude"
                .to_string(),
        ),
        (
            vec!["--density-altitude=200000"],
            "Error: --density-altitude = 200000 ft: no air has this density altitude".to_string(),
        ),
    ] {
        let mut args = base.to_vec();
        args.extend(extra.iter().copied());
        let output = run(home.path(), &args);
        assert_eq!(output.status.code(), Some(1), "{extra:?}");
        let error = error_line(&output);
        assert!(
            error.starts_with(&refusal),
            "{extra:?}: expected `{refusal}`, got `{error}`"
        );
    }
}

/// The bridge's `card.*` and `true.*` commands, and the C ABI that carries them, name a refused
/// temperature or pressure by the request's own member, in the request's own units, at its path
/// in `error.details`, as `solve` locates its own. In 0.49.0 they reached the solver and came
/// back as `atmosphere.pressure = -33.86 hPa: ...` with no path.
#[test]
fn the_bridge_names_a_refused_reading_by_its_member_in_its_units() {
    let card = |field: &str, value: f64| {
        let mut request = json!({
            "muzzle_velocity": 2700.0, "ballistic_coefficient": 0.475, "mass": 168.0,
            "diameter": 0.308, "zero_distance": 100.0, "start": 100.0, "end": 300.0,
            "step": 100.0, "humidity": 50.0
        });
        request[field] = json!(value);
        request
    };
    let model = |field: &str, value: f64| {
        let mut model = json!({
            "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
            "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
            "zero_distance_yd": 100.0, "sight_height_in": 2.0,
            "temperature_f": 59.0, "pressure_inhg": 29.92,
            "humidity_pct": 50.0, "altitude_ft": 0.0
        });
        model[field] = json!(value);
        model
    };
    let fit = |field: &str, value: f64| {
        json!({
            "model": model(field, value), "drop_unit": "mil",
            "observations": [{"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1}],
            "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                       "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
            "predictions": []
        })
    };
    let wind = |field: &str, value: f64| {
        let mut request = json!({
            "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
            "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
            "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
            "sight_height_in": 2.0, "temperature_f": 59.0, "pressure_inhg": 29.92,
            "humidity_pct": 50.0, "altitude_ft": 0.0,
            "twist": {"rate_in": 11.0, "right_hand": true},
            "earth": null, "called_crosswind_mph": null
        });
        request[field] = json!(value);
        request
    };
    let absolute_zero = "temperature must be above absolute zero (-459.67 °F)";
    let cases = [
        (
            "card.range_table",
            card("pressure", -1.0),
            "$.pressure",
            "pressure = -1 inHg: pressure must be greater than 0".to_string(),
        ),
        (
            "card.come_ups",
            card("temperature", -500.0),
            "$.temperature",
            format!("temperature = -500 °F: {absolute_zero}"),
        ),
        (
            "true.fit",
            fit("pressure_inhg", 0.0),
            "$.model.pressure_inhg",
            "pressure_inhg = 0 inHg: pressure must be greater than 0".to_string(),
        ),
        (
            "true.fit",
            fit("temperature_f", -500.0),
            "$.model.temperature_f",
            format!("temperature_f = -500 °F: {absolute_zero}"),
        ),
        (
            "true.wind",
            wind("pressure_inhg", -1.0),
            "$.pressure_inhg",
            "pressure_inhg = -1 inHg: pressure must be greater than 0".to_string(),
        ),
        (
            "true.wind",
            wind("temperature_f", -500.0),
            "$.temperature_f",
            format!("temperature_f = -500 °F: {absolute_zero}"),
        ),
    ];
    for (command, request, path, message) in cases {
        let out = bridge(command, request.clone());
        assert_eq!(out["ok"], false, "{command}: {out}");
        let details = &out["error"]["details"]["error"];
        assert_eq!(details["code"], "invalid_value", "{command}: {out}");
        assert_eq!(details["path"], path, "{command}: {out}");
        assert_eq!(details["message"], message, "{command}: {out}");
        assert_eq!(
            out["error"]["message"],
            format!("{command} failed: {message}"),
            "{command}"
        );
        #[cfg(feature = "ffi")]
        {
            let envelope = json!({"api_version": 1, "command": command, "request": request});
            let text = std::ffi::CString::new(envelope.to_string()).unwrap();
            let raw =
                unsafe { ballistics_engine::bridge::ffi::ballistics_bridge_call(text.as_ptr()) };
            let through_ffi: Value =
                serde_json::from_str(unsafe { std::ffi::CStr::from_ptr(raw) }.to_str().unwrap())
                    .unwrap();
            unsafe { ballistics_engine::bridge::ffi::ballistics_bridge_free(raw) };
            assert_eq!(through_ffi, out, "{command}: the C ABI is the bridge");
        }
    }
}

/// Air that cannot exist at a density altitude's station above the troposphere is blamed on the
/// density altitude, with the unit to check: never on a `--pressure` nobody typed.
#[test]
fn impossible_air_at_a_density_altitudes_station_blames_the_density_altitude() {
    let home = Home::new();
    let output = run(
        home.path(),
        &[
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
            "--density-altitude",
            "146000",
            "--temperature",
            "150",
        ],
    );
    assert_eq!(output.status.code(), Some(1));
    let error = error_line(&output);
    assert!(
        error.starts_with("Error: --humidity = 50 at 150 °F and ")
            && error.contains("no such air exists")
            && error.contains("; --density-altitude gives a station pressure of ")
            && error.ends_with(": check the density altitude's unit"),
        "{error}"
    );
}

// --- the figures the docs quote, from the engine -------------------------------------------------

fn read(relative: &str) -> String {
    std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(relative))
        .unwrap_or_else(|error| panic!("{relative}: {error}"))
}

/// `text` with every run of whitespace one space, so a quote may wrap.
fn normalized(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn assert_quotes(document: &str, quotes: &[String]) {
    // A source file's docs and comments, read without their `///`, `//!` and `//` markers, so a
    // quote may wrap there too.
    let text = read(document);
    let text = if document.ends_with(".rs") {
        text.lines()
            .map(|line| {
                let line = line.trim_start();
                ["///", "//!", "//"]
                    .iter()
                    .find_map(|marker| line.strip_prefix(marker))
                    .unwrap_or(line)
            })
            .collect::<Vec<_>>()
            .join("\n")
    } else {
        text
    };
    let text = normalized(&text);
    for quote in quotes {
        println!("{document}: engine measured: {quote}");
        assert!(
            text.contains(&normalized(quote)),
            "{document} does not quote the engine: `{quote}`"
        );
    }
}

/// Published examples must exercise every removed station bound, and distinguish the two
/// pressure-altitude formulas outside the troposphere. Print the engine output for reproduction.
#[test]
fn station_rule_documentation_examples_are_measured() {
    let home = Home::new();
    let mut cli_quotes = Vec::new();
    for extra in [
        vec!["--pressure", "5"],
        vec!["--density-altitude=-60000", "--humidity", "0"],
        vec![
            "--auto-zero",
            "100",
            "--zero-density-altitude=-60000",
            "--zero-humidity",
            "0",
        ],
    ] {
        let mut args = vec![
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
        ];
        args.extend(extra);
        let output = run(home.path(), &args);
        println!(
            "$ ballistics {}\n{}\n{}",
            args.join(" "),
            String::from_utf8_lossy(&output.stdout),
            stderr(&output)
        );
        assert!(output.status.success());
        let stdout = String::from_utf8_lossy(&output.stdout);
        cli_quotes.push(
            stdout
                .lines()
                .find(|line| line.starts_with("║ Impact Velocity:"))
                .unwrap()
                .to_string(),
        );
    }

    let mut solve_rows = Vec::new();
    for altitude in [-6000.0, 90000.0] {
        let atmosphere = json!({"altitude_m": altitude, "temperature_k": 288.15,
            "pressure_pa": 101325.0, "relative_humidity": 0.0});
        let response = bridge("solve", solve_request(atmosphere));
        println!(
            "solve request: {}\nresponse: {}",
            solve_request(json!({"altitude_m": altitude,
            "temperature_k": 288.15, "pressure_pa": 101325.0, "relative_humidity": 0.0})),
            response
        );
        assert_eq!(response["ok"], true);
        let summary = &response["result"]["summary"];
        solve_rows.push(format!(
            "| {altitude:.0} | {:.3} | {:.3} |",
            summary["actual_range_m"].as_f64().unwrap(),
            summary["terminal_speed_mps"].as_f64().unwrap()
        ));
    }

    let mut pressure_rows = Vec::new();
    for (altitude, temperature, pressure) in [
        (20_000.0, standard_air(20_000.0).0, standard_air(20_000.0).1),
        (0.0, 288.15, 1013.25),
        (0.0, 288.15, 202650.0),
    ] {
        let atmosphere = json!({"altitude_m": altitude, "temperature_k": temperature,
            "pressure_pa": pressure, "relative_humidity": 0.0});
        let response = density_altitude(&atmosphere).unwrap();
        println!("atmosphere.density_altitude request: {atmosphere}\nresponse: {response}");
        pressure_rows.push(format!(
            "| {temperature:.2} | {pressure:.5} | {:.0} | {:.0} | {:.0} | {:.0} |",
            response["pressure_altitude"]["m"].as_f64().unwrap(),
            nws_pressure_altitude_ft(pressure / 100.0) * 0.3048,
            response["density_altitude"]["faa_rule"]["m"]
                .as_f64()
                .unwrap(),
            response["density_altitude"]["density_matched"]["m"]
                .as_f64()
                .unwrap()
        ));
    }
    println!(
        "CLI impact lines:\n{}\nSolve rows:\n{}\nPressure rows:\n{}",
        cli_quotes.join("\n"),
        solve_rows.join("\n"),
        pressure_rows.join("\n")
    );
    assert_quotes("CLI_USAGE.md", &cli_quotes);
    assert_quotes("CHANGELOG.md", &cli_quotes);
    for doc in ["CHANGELOG.md", "docs/SOLVE_JSON_V1.md"] {
        assert_quotes(doc, &solve_rows);
    }
    for doc in ["CHANGELOG.md", "docs/ATMOSPHERE_DENSITY_ALTITUDE.md"] {
        assert_quotes(doc, &pressure_rows);
    }
}

/// Overview docs must carry the station rule and format exception, and must not label the
/// layered pressure-altitude result as the NWS continuation.
#[test]
fn station_rule_documentation_names_the_rule_and_format_exception() {
    for document in ["README.md", "src/lib.rs", "include/ballistics_bridge.h"] {
        let text = normalized(&read(document));
        for phrase in [
            "station pressure",
            "finite",
            "density altitude",
            "profile.export_a7p",
            "-100..100",
            "300..1500",
        ] {
            assert!(text.contains(phrase), "{document}: missing {phrase}");
        }
    }
    let doc = read("docs/ATMOSPHERE_DENSITY_ALTITUDE.md");
    assert!(
        !doc.contains("| `pressure_altitude` | NWS pressure altitude of the station pressure. |")
    );
    assert!(doc.contains("`nws_pressure_altitude_ft`"));
}

#[cfg(feature = "profile-export")]
#[test]
fn station_rule_documentation_keeps_the_a7p_format_bounds() {
    let home = Home::new();
    home.profile("format-weather", "metric", &["--humidity", "0"], &[]);
    let profile: Value = serde_json::from_str(
        &std::fs::read_to_string(home.profile_path("format-weather")).unwrap(),
    )
    .unwrap();
    let mut quotes = Vec::new();
    for (field, value, accepted) in [
        ("temperature", -100.0, true),
        ("temperature", 100.0, true),
        ("temperature", -101.0, false),
        ("temperature", 101.0, false),
        ("pressure", 300.0, true),
        ("pressure", 1500.0, true),
        ("pressure", 299.0, false),
        ("pressure", 1501.0, false),
    ] {
        let mut input = profile.clone();
        input[field] = json!(value);
        input[format!("{field}_supplied")] = json!(true);
        let validated = bridge("profile.validate", input.clone());
        assert_eq!(
            validated["result"]["valid"], true,
            "the air exists: {validated}"
        );
        let exported = bridge("profile.export_a7p", input);
        assert_eq!(exported["ok"], accepted, "{field} = {value}: {exported}");
        if !accepted {
            let message = exported["error"]["message"].as_str().unwrap();
            println!("profile.export_a7p {field} = {value}: {message}");
            quotes.push(message.to_string());
        } else {
            println!("profile.export_a7p {field} = {value}: ok");
        }
    }
    assert_quotes("CLI_USAGE.md", &quotes);
}

/// Every figure the CHANGELOG, CLI_USAGE.md, the two wire documents and the `atmosphere` docs
/// quote for the one rule for each station pressure, altitude and density altitude, recomputed
/// from the engine, and every example refusal, run.
#[test]
fn the_docs_quote_the_engine() {
    use ballistics_engine::atmosphere::{pressure_altitude_m, station_for_density_altitude};
    let metres = |metres: f64| format!("{metres:.0}");
    // Standard pressure at 20 km, by the NWS formula continued and by the standard atmosphere.
    let (_, standard_20_km_pa) = standard_air(20_000.0);
    let nws_20_km = nws_pressure_altitude_ft(standard_20_km_pa / 100.0) * 0.3048;
    let walked_20_km = pressure_altitude_m(standard_20_km_pa / 100.0);
    let hpa_as_pa = pressure_altitude_m(10.1325);
    let twice_sea_level = pressure_altitude_m(2026.5);
    // The band's edges: the pressure the NWS formula reads as 11000 m and as -5000 m, and the
    // standard atmosphere's altitude of each.
    let edge = |target_m: f64| {
        let (mut low, mut high) = (1.0_f64, 3000.0_f64);
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if nws_pressure_altitude_ft(middle) * 0.3048 > target_m {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    };
    let (top, bottom) = (edge(11_000.0), edge(-5_000.0));
    let top_step = pressure_altitude_m(top * (1.0 - 1e-12)) - 11_000.0;
    let bottom_step = -5_000.0 - pressure_altitude_m(bottom * (1.0 + 1e-12));
    // The pressure altitude by the FAA rule (the station altitude the inverse returns) of the
    // last density altitude that has air, without and with a temperature.
    let reach = |temperature_c: Option<f64>| {
        let (mut low, mut high) = (40_000.0_f64, 70_000.0_f64);
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if station_for_density_altitude(middle, temperature_c).is_ok() {
                low = middle;
            } else {
                high = middle;
            }
        }
        station_for_density_altitude(low, temperature_c).unwrap().0
    };
    let (no_temperature, with_temperature) = (reach(None), reach(Some(15.0)));
    let feet = |metres: f64| format!("{:.0}", metres / 0.3048);
    // The pressure beyond which the density model gives no density, at 15 °C and dry.
    let mut too_dense_hpa = 1e100_f64;
    while ballistics_engine::atmosphere::calculate_air_density_cipm(15.0, too_dense_hpa, 0.0) > 0.0
    {
        too_dense_hpa *= 1.01;
    }
    let too_dense = format!("{:.0e}", too_dense_hpa * 100.0);
    assert_quotes(
        "CHANGELOG.md",
        &[
            format!(
                "standard pressure at 20 km reads {} m, where the NWS formula continued reads {} m",
                metres(walked_20_km),
                metres(nws_20_km)
            ),
            format!(
                "hPa sent as `pressure_pa` (1013.25 Pa, dry) reads {} m",
                metres(hpa_as_pa)
            ),
            format!(
                "up to a pressure altitude of {} ft ({} m), where the FAA rule's standard \
                 temperature reaches 0 K; with one, up to {} ft ({} m)",
                feet(no_temperature),
                metres(no_temperature),
                feet(with_temperature),
                metres(with_temperature)
            ),
            format!("(above about {too_dense} Pa)"),
        ],
    );
    assert_quotes(
        "CLI_USAGE.md",
        &[format!(
            "of {} ft ({} m), where the FAA rule's standard temperature reaches 0 K; with one, up \
             to {} ft ({} m)",
            feet(no_temperature),
            metres(no_temperature),
            feet(with_temperature),
            metres(with_temperature)
        )],
    );
    assert_quotes(
        "docs/SOLVE_JSON_V1.md",
        &[format!("(above about {too_dense} Pa)")],
    );
    assert_quotes(
        "docs/ATMOSPHERE_DENSITY_ALTITUDE.md",
        &[
            format!(
                "standard pressure at 20 km reads {} m, where the NWS formula continued reads {} m",
                metres(walked_20_km),
                metres(nws_20_km)
            ),
            format!(
                "({} m at 11 km, {} m at −5 km)",
                metres(top_step),
                metres(bottom_step)
            ),
        ],
    );
    assert_quotes(
        "src/atmosphere.rs",
        &[
            format!(
                "(standard pressure at 20 km, 5474.89 Pa, reads {} m)",
                metres(nws_pressure_altitude_ft(54.7489) * 0.3048)
            ),
            format!(
                "the NWS formula reads {:.1} Pa as 11000 m and the standard atmosphere reads it as \
                 {} m, so the figure steps up by {} m there",
                top * 100.0,
                metres(11_000.0 + top_step),
                metres(top_step)
            ),
            format!(
                "at the bottom edge ({:.0} Pa) it steps down by {} m",
                bottom * 100.0,
                metres(bottom_step)
            ),
            format!("(above about {too_dense} Pa"),
        ],
    );
    assert_quotes(
        "src/atmosphere_service.rs",
        &[format!(
            "(1013.25 Pa) reads {} m, twice sea level {} m",
            metres(hpa_as_pa),
            metres(twice_sea_level)
        )],
    );

    // The example runs CLI_USAGE.md and the CHANGELOG print, run.
    let home = Home::new();
    let base = [
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
    ];
    let line = |extra: &[&str], prefix: &str| {
        let mut args = base.to_vec();
        args.extend_from_slice(extra);
        let output = run(home.path(), &args);
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        text.lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("{extra:?}: no `{prefix}` line in:\n{text}"))
            .to_string()
    };
    let flown = line(&["--pressure", "5"], "║ Impact Velocity:");
    let refused = line(&["--pressure=-1"], "Error: ");
    let density_altitude = line(&["--density-altitude", "200000"], "Error: ");
    assert_quotes(
        "CLI_USAGE.md",
        &[flown, refused.clone(), density_altitude.clone()],
    );
    let refused = refused.strip_prefix("Error: ").unwrap().to_string();
    let density_altitude = density_altitude
        .strip_prefix("Error: ")
        .unwrap()
        .to_string();
    assert_quotes("CHANGELOG.md", &[refused, density_altitude]);
}

// --- the C ABI's density-altitude conversions -----------------------------------------------------

/// The C ABI's three density-altitude conversions apply the one density-altitude rule, as every
/// other surface does: where `station_for_density_altitude` gives a station they return it, bit
/// for bit, and where it refuses one (no air has the density altitude, or the explicit
/// temperature is not above absolute zero) all three return NaN, the signal they already gave
/// for a density altitude that is not a number. They used to hand back a station no air has,
/// finite: 44300 m gave -273.26 °C at 1.75e-17 hPa, and 50000 m a temperature of -310.35 °C
/// beside a NaN pressure, where the CLI, the browser terminal and saved profiles refuse.
#[cfg(feature = "ffi")]
#[test]
fn the_c_abis_density_altitude_conversions_apply_the_one_density_altitude_rule() {
    use ballistics_engine::atmosphere::station_for_density_altitude;
    use ballistics_engine::ffi::{
        ballistics_density_altitude_altitude_m, ballistics_density_altitude_pressure_hpa,
        ballistics_density_altitude_temperature_c, FFI_NO_EXPLICIT_TEMPERATURE,
    };
    let mut cases = density_altitudes();
    for density_altitude_m in [0.0, 44_290.0, 44_300.0, 60_000.0] {
        cases.push((density_altitude_m, None));
    }
    for temperature_c in [-273.15, -300.0, f64::INFINITY, f64::NEG_INFINITY] {
        cases.push((1000.0, Some(temperature_c)));
    }
    for density_altitude_m in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        cases.push((density_altitude_m, None));
    }
    let mut refused = 0;
    let mut failures = Vec::new();
    for (density_altitude_m, temperature_c) in cases {
        let explicit = temperature_c.unwrap_or(FFI_NO_EXPLICIT_TEMPERATURE);
        let through_c = [
            ballistics_density_altitude_altitude_m(density_altitude_m, explicit),
            ballistics_density_altitude_temperature_c(density_altitude_m, explicit),
            ballistics_density_altitude_pressure_hpa(density_altitude_m, explicit),
        ];
        let case =
            format!("{density_altitude_m} m, {temperature_c:?} C: (A, T, P) = {through_c:?}");
        match station_for_density_altitude(density_altitude_m, temperature_c) {
            Ok(station) => {
                let station = [station.0, station.1, station.2];
                if through_c.map(f64::to_bits) != station.map(f64::to_bits) {
                    failures.push(format!("{case}: the rule gives {station:?}"));
                }
            }
            Err(error) => {
                refused += 1;
                if !through_c.iter().all(|value| value.is_nan()) {
                    failures.push(format!("{case}: the rule refuses it ({error}), not NaN"));
                }
            }
        }
        // The station every terminal flies or refuses for it, dry.
        let flown = through_c.iter().all(|value| !value.is_nan());
        if flown != density_altitude_flies(density_altitude_m, temperature_c) {
            failures.push(format!("{case}: the terminals disagree"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} disagreements:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert!(refused >= 12, "{refused} refusals");
}

// --- a stored reading's warning -------------------------------------------------------------------

/// A stored reading a flight refuses is warned about in one sentence that says what to edit, by
/// `profile show` (beside the file it is stored in) and by the bridge's `profile.validate`, word
/// for word alike.
#[test]
fn a_refused_stored_reading_is_warned_about_in_one_sentence() {
    let home = Home::new();
    home.profile("p", "imperial", &[], &[("pressure", json!(-1.0))]);
    let warning =
        "profile 'p' field pressure = -1 inHg: pressure must be greater than 0; edit the \
                   humidity, temperature, pressure, altitude or density altitude stored in the \
                   profile";
    let output = run(home.path(), &["-u", "imperial", "profile", "show", "p"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let expected = format!(
        "warning: {warning} (the profile is stored at {})",
        home.profile_path("p").display()
    );
    let text = stderr(&output);
    let warnings: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("warning: "))
        .collect();
    assert_eq!(warnings, [expected.as_str()]);
    let stored: Value =
        serde_json::from_str(&std::fs::read_to_string(home.profile_path("p")).unwrap()).unwrap();
    let validate = bridge("profile.validate", stored);
    assert_eq!(validate["ok"], true, "{validate}");
    assert_eq!(
        validate["result"]["warnings"],
        json!([warning]),
        "{validate}"
    );
}

// --- absurd readings ------------------------------------------------------------------------------

/// `message` writes no number out in more than 12 digits, and quotes each of `quotes`.
fn assert_compact(case: &str, message: &str, quotes: &[&str]) -> Option<String> {
    let longest = message
        .split(|c: char| !c.is_ascii_digit())
        .map(str::len)
        .max()
        .unwrap_or(0);
    if longest > 12 {
        return Some(format!(
            "{case}: a number {longest} digits long in `{message}`"
        ));
    }
    quotes
        .iter()
        .find(|quote| !message.contains(*quote))
        .map(|quote| format!("{case}: `{quote}` is not in `{message}`"))
}

/// An absurd reading is refused in a few characters on every surface, as the SI wire's echo
/// already was (`1e300 K`): a reading, and every altitude its refusal quotes, beyond 10^7 in the
/// surface's unit is written in exponent form to six significant digits, and a typed value
/// beyond 10^7 or below 0.001 is echoed in exponent form, exactly. A pressure is never stated
/// as `0` when it is not 0. Removing the pressure band let `--pressure 1e160` through to the
/// impossible-air rule, whose refusal printed the pressure twice in 161 digits and its pressure
/// altitude in 35.
#[test]
fn absurd_readings_are_refused_in_a_few_characters_on_every_surface() {
    let home = Home::new();
    home.profile(
        "big",
        "metric",
        &[],
        &[("pressure", json!(1e160)), ("humidity", json!(0.0))],
    );
    let metric = |command: &[&str], flags: &[&str]| -> Vec<String> {
        let mut args: Vec<String> = METRIC_LOAD.iter().take(2).map(|s| s.to_string()).collect();
        args.extend(command.iter().map(|s| s.to_string()));
        args.extend(METRIC_LOAD.iter().skip(2).map(|s| s.to_string()));
        args.extend(flags.iter().map(|s| s.to_string()));
        args
    };
    let imperial = |flags: &[&str]| -> Vec<String> {
        let mut args: Vec<String> = [
            "-u",
            "imperial",
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
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        args.extend(flags.iter().map(|s| s.to_string()));
        args
    };
    let trajectory = ["trajectory", "--max-range", "300"];
    let cli: Vec<(&str, Vec<String>, Vec<&str>)> = vec![
        (
            "--pressure 1e160 hPa",
            metric(&trajectory, &["--pressure=1e160", "--humidity=0"]),
            vec![
                "at 15 °C and 1e160 hPa",
                "--pressure = 1e160 hPa is a pressure altitude of -",
            ],
        ),
        (
            "--pressure 1e160 inHg",
            imperial(&["--pressure=1e160", "--humidity=0"]),
            vec![
                "at 59 °F and 1e160 inHg",
                "--pressure = 1e160 inHg is a pressure altitude of -",
            ],
        ),
        (
            "--pressure 1e-300 hPa",
            metric(&trajectory, &["--pressure=1e-300", "--humidity=50"]),
            vec![
                "at 15 °C and 1e-300 hPa",
                "--pressure = 1e-300 hPa is a pressure altitude of ",
            ],
        ),
        (
            "--pressure -1e300",
            metric(&trajectory, &["--pressure=-1e300"]),
            vec!["Error: --pressure -1e300 hPa is at or below zero"],
        ),
        (
            "--temperature -1e300",
            metric(&trajectory, &["--temperature=-1e300"]),
            vec!["Error: --temperature -1e300 C is at or below absolute zero (-273.15 C)"],
        ),
        (
            "--temperature 1e300",
            metric(&trajectory, &["--temperature=1e300", "--humidity=50"]),
            vec!["Error: --humidity = 50 at 1e300 °C and 1013.25 hPa: no such air exists"],
        ),
        (
            "--density-altitude 1e300 m",
            metric(&trajectory, &["--density-altitude=1e300"]),
            vec![
                "Error: --density-altitude = 1e300 m: no air has this density altitude: at its \
                 pressure altitude by the FAA rule, 1e300 m,",
            ],
        ),
        (
            "--density-altitude 1e300 ft",
            imperial(&["--density-altitude=1e300"]),
            vec![
                "Error: --density-altitude = 1e300 ft: no air has this density altitude: at its \
                 pressure altitude by the FAA rule, 1e300 ft,",
            ],
        ),
        (
            "--altitude 1e300 with a QNH",
            metric(
                &trajectory,
                &["--altitude=1e300", "--pressure=1013", "--pressure-type=qnh"],
            ),
            vec!["reduced to the station altitude of 1e300 m it leaves no station pressure"],
        ),
        (
            "a stored pressure of 1e160 hPa",
            metric(
                &[
                    "come-ups",
                    "--profile",
                    "big",
                    "--zero-distance",
                    "100",
                    "--end",
                    "300",
                ],
                &[],
            ),
            vec![
                "at 15 °C and 1e160 hPa",
                "profile 'big' field pressure = 1e160 hPa is a pressure altitude of -",
            ],
        ),
    ];
    let outputs = std::thread::scope(|scope| {
        let jobs: Vec<_> = cli
            .iter()
            .map(|(_, args, _)| {
                let home = home.path();
                scope.spawn(move || {
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    run(home, &args)
                })
            })
            .collect();
        jobs.into_iter()
            .map(|job| job.join().unwrap())
            .collect::<Vec<_>>()
    });
    let mut failures = Vec::new();
    for ((case, _, quotes), output) in cli.iter().zip(&outputs) {
        let error = error_line(output);
        if output.status.code() != Some(1) || error.is_empty() {
            failures.push(format!("{case}: not refused: {}", stderr(output)));
        } else if let Some(failure) = assert_compact(case, &error, quotes) {
            failures.push(failure);
        }
    }
    // The wire surfaces, in their own units.
    let atmosphere = |pressure_pa: f64, humidity: f64| {
        json!({"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": pressure_pa,
               "relative_humidity": humidity})
    };
    for (case, atmosphere, quotes) in [
        (
            "solve, 1e300 Pa",
            atmosphere(1e300, 0.0),
            vec![
                "at 288.15 K and 1e300 Pa",
                "pressure_pa = 1e300 Pa is a pressure altitude of -",
            ],
        ),
        (
            "solve, 1e-300 Pa",
            atmosphere(1e-300, 0.5),
            vec!["at 288.15 K and 1e-300 Pa"],
        ),
        (
            "solve, a QNH at 1e300 m",
            json!({"altitude_m": 1e300, "temperature_k": 288.15, "pressure_pa": 101325.0,
                   "relative_humidity": 0.5, "pressure_reference": "qnh"}),
            vec!["at 1e300 m"],
        ),
    ] {
        for (surface, refusal) in [
            ("solve", solve(atmosphere.clone()).err()),
            (
                "atmosphere.density_altitude",
                density_altitude(&atmosphere).err(),
            ),
        ] {
            let case = format!("{case}, {surface}");
            match refusal {
                None => failures.push(format!("{case}: not refused")),
                Some(error) => {
                    let message = error["message"].as_str().unwrap_or_default();
                    if let Some(failure) = assert_compact(&case, message, &quotes) {
                        failures.push(failure);
                    }
                }
            }
        }
    }
    let card = json!({
        "muzzle_velocity": 2700.0, "ballistic_coefficient": 0.475, "mass": 168.0,
        "diameter": 0.308, "zero_distance": 100.0, "start": 100.0, "end": 300.0,
        "step": 100.0, "humidity": 0.0, "pressure": 1e160
    });
    for command in ["card.range_table", "card.come_ups"] {
        let out = bridge(command, card.clone());
        let message = out["error"]["details"]["error"]["message"]
            .as_str()
            .unwrap_or_default();
        if out["ok"] != false {
            failures.push(format!("{command}: not refused: {out}"));
        } else if let Some(failure) = assert_compact(
            command,
            message,
            &[
                "at 59 °F and 1e160 inHg",
                "pressure = 1e160 inHg is a pressure altitude of -",
            ],
        ) {
            failures.push(failure);
        }
    }
    assert!(
        failures.is_empty(),
        "{} refusals not compact:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn long_decimal_station_readings_have_compact_refusals() {
    let home = Home::new();
    for (field, raw, quote) in [
        ("pressure", format!("1{}", "0".repeat(300)), "1e300 hPa"),
        ("pressure", format!("0.{}1", "0".repeat(299)), "1e-300 hPa"),
        ("temperature", format!("1{}", "0".repeat(300)), "1e300 °C"),
        (
            "humidity",
            format!("105.0000001{}", "0".repeat(300)),
            "105.0000001",
        ),
        (
            "density-altitude",
            format!("1{}", "0".repeat(300)),
            "1e300 m",
        ),
    ] {
        let flag = format!("--{field}={raw}");
        let mut args = METRIC_LOAD.to_vec();
        args.splice(2..2, ["trajectory", "--max-range", "300"]);
        args.push(&flag);
        let output = run(home.path(), &args);
        let message = stderr(&output);
        assert!(!output.status.success(), "{field}: {message}");
        assert!(message.contains(&format!("--{field}")), "{message}");
        assert_eq!(assert_compact(field, &message, &[quote]), None);

        let column = field.replace('-', "_").to_uppercase();
        let csv = home.file("decimal.csv", &format!("NAME,{column}\nDecimal,{raw}\n"));
        args.pop();
        args.extend(["--location", &csv, "--site", "Decimal"]);
        let output = run(home.path(), &args);
        let message = error_line(&output);
        assert!(!output.status.success(), "{field}: {}", stderr(&output));
        assert_eq!(
            assert_compact(field, &message, &[quote, &format!("CSV column {column}")]),
            None
        );

        let stored_field = field.replace('-', "_");
        home.profile(
            "decimal",
            "metric",
            &[],
            &[(stored_field.as_str(), json!(raw.parse::<f64>().unwrap()))],
        );
        let stored_path = home.profile_path("decimal");
        let stored_text = std::fs::read_to_string(&stored_path).unwrap();
        let mut stored: Value = serde_json::from_str(&stored_text).unwrap();
        stored[&stored_field] = json!("decimal-reading");
        let stored_text = stored.to_string().replace("\"decimal-reading\"", &raw);
        std::fs::write(stored_path, &stored_text).unwrap();
        let stored: Value = serde_json::from_str(&stored_text).unwrap();
        let out = bridge("profile.validate", stored);
        let message = out["result"]["warnings"][0].as_str().unwrap();
        assert_eq!(
            assert_compact(field, message, &[quote, &format!("field {stored_field}")]),
            None
        );
    }
}

#[test]
fn absurd_bridge_readings_keep_the_request_field_and_unit() {
    let station = Station {
        temperature_c: 15.0,
        pressure_hpa: Some(1013.25),
        altitude_m: 0.0,
        humidity: 50.0,
    };
    for (card_field, model_field, value, unit) in [
        ("pressure", "pressure_inhg", 1e300, "inHg"),
        ("pressure", "pressure_inhg", 1e-300, "inHg"),
        ("temperature", "temperature_f", -1e300, "°F"),
        ("humidity", "humidity_pct", 1e300, ""),
    ] {
        let mut card = card_request(&station);
        card["temperature_unit"] = json!("fahrenheit");
        card["temperature"] = json!(59.0);
        card["pressure_unit"] = json!("inhg");
        card["pressure"] = json!(29.92);
        card[card_field] = json!(value);
        let mut model = truing_model(&station);
        model[model_field] = json!(value);
        let fit = json!({
            "model": model, "drop_unit": "mil",
            "observations": [{"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1}],
            "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                       "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
            "predictions": []
        });
        let mut wind = wind_truing_request(&station);
        wind[model_field] = json!(value);
        for (command, request, field, path) in [
            (
                "card.range_table",
                card.clone(),
                card_field,
                format!("$.{card_field}"),
            ),
            ("card.come_ups", card, card_field, format!("$.{card_field}")),
            (
                "true.fit",
                fit,
                model_field,
                format!("$.model.{model_field}"),
            ),
            ("true.wind", wind, model_field, format!("$.{model_field}")),
        ] {
            let out = bridge(command, request.clone());
            assert_eq!(out["ok"], false, "{command}: {out}");
            let error = &out["error"]["details"]["error"];
            assert_eq!(error["path"], path, "{command}: {out}");
            let message = error["message"].as_str().unwrap();
            assert_eq!(
                assert_compact(command, message, &[field, &format!("{value:e}"), unit]),
                None
            );
            #[cfg(feature = "ffi")]
            {
                let envelope = json!({"api_version": 1, "command": command, "request": request});
                let text = std::ffi::CString::new(envelope.to_string()).unwrap();
                let raw = unsafe {
                    ballistics_engine::bridge::ffi::ballistics_bridge_call(text.as_ptr())
                };
                let through_ffi: Value = serde_json::from_str(
                    unsafe { std::ffi::CStr::from_ptr(raw) }.to_str().unwrap(),
                )
                .unwrap();
                unsafe { ballistics_engine::bridge::ffi::ballistics_bridge_free(raw) };
                assert_eq!(through_ffi, out, "{command}: the C ABI is the bridge");
            }
        }
    }
}

// --- a refusal found in flight ----------------------------------------------------------------------

/// Air that cannot exist aloft, which a solve finds only in flight, is refused naming the
/// temperature where it came from on every command, as the up-front refusals do: a saved
/// profile's field, a `--location` column, or the flag. At 0.0065 K the standard lapse reaches
/// 0 K 1 m above the station, which a 600 m zero's flight climbs past. Only `trajectory` and
/// `zero` named the source; `come-ups`, `range-table`, `wind-card` and every other command named
/// `--temperature`, a flag nobody typed.
#[test]
fn a_refusal_found_in_flight_names_the_temperature_where_it_came_from() {
    let home = Home::new();
    let cold = [
        "--temperature=-273.1435",
        "--pressure=0.0187",
        "--humidity=0",
    ];
    home.profile("cold", "metric", &cold, &[]);
    let csv = home.file(
        "cold.csv",
        "NAME,TEMP,PRESSURE,HUMIDITY\nSITE,-273.1435,0.0187,0\n",
    );
    let aloft = |label: &str| {
        format!(
            "Error: {label} = -273.14 °C: the air aloft would fall below absolute zero: from this \
             temperature the standard lapse reaches 0 K at 1 m above mean sea level, which this \
             flight reaches (check the temperature and its unit)"
        )
    };
    let with = |head: &[&str], tail: &[&str]| -> Vec<String> {
        let mut args: Vec<String> = METRIC_LOAD.iter().take(2).map(|s| s.to_string()).collect();
        args.extend(head.iter().map(|s| s.to_string()));
        args.extend(METRIC_LOAD.iter().skip(2).map(|s| s.to_string()));
        args.extend(tail.iter().map(|s| s.to_string()));
        args
    };
    let stored = ["--profile", "cold"];
    let mut cases: Vec<(String, Vec<String>, String)> = Vec::new();
    for command in ["come-ups", "range-table", "wind-card"] {
        let head = [command, "--zero-distance", "600", "--end", "800"];
        cases.push((
            format!("{command} --profile"),
            {
                let mut args: Vec<String> =
                    ["-u", "metric"].iter().map(|s| s.to_string()).collect();
                args.extend(head.iter().chain(&stored).map(|s| s.to_string()));
                args
            },
            aloft("profile 'cold' field temperature"),
        ));
        cases.push((
            format!("{command} with flags"),
            with(&head, &cold),
            aloft("--temperature"),
        ));
    }
    cases.push((
        "trajectory --location".into(),
        with(
            &["trajectory", "--max-range", "800", "--auto-zero", "600"],
            &["--location", &csv, "--site", "SITE"],
        ),
        aloft("--location CSV column TEMP"),
    ));
    let scenarios = home.file(
        "scenarios.json",
        r#"{"version":1,"scenarios":[{"name":"steady","segments":["4:90:1000"]}]}"#,
    );
    let reticle = home.file("reticle.json", r#"{"name":"fixture","focal_plane":"ffp","reference_magnification":1.0,"marks":[{"down_mil":0.0,"right_mil":0.0,"kind":"center"}]}"#);
    for head in [
        vec![
            "hold-corridor",
            "--zero-distance",
            "600",
            "--ranges",
            "600,800",
            "--scenarios",
            &scenarios,
        ],
        vec![
            "reticle",
            "hold",
            "--zero-distance",
            "600",
            "--range",
            "800",
            "--reticle-json",
            &reticle,
        ],
        vec![
            "bdc-match",
            "--zero-distance",
            "600",
            "--mark-range",
            "1:600",
            "--mark-range",
            "2:800",
        ],
        vec!["mpbr", "--vital-zone", "200"],
        vec![
            "compare",
            "--zero-distance",
            "600",
            "--end",
            "800",
            "--load",
            "B:G1:0.475:168:2600",
        ],
        vec!["mark-to-range", "--zero-distance", "600", "--mark", "2"],
        vec!["optimal-zero", "--target", "600:200", "--target", "800:200"],
        vec![
            "adaptive-card",
            "--zero-distance",
            "600",
            "--start",
            "600",
            "--end",
            "800",
        ],
        vec![
            "plan-truing",
            "--zero-distance",
            "600",
            "--measurement-resolution",
            "0.25",
            "--candidate-ranges",
            "600,800",
        ],
    ] {
        cases.push((
            head[0].into(),
            [vec!["-u", "metric"], head, stored.to_vec()]
                .concat()
                .iter()
                .map(|s| s.to_string())
                .collect(),
            aloft("profile 'cold' field temperature"),
        ));
    }
    cases.push((
        "zero-day flag".into(),
        with(
            &["trajectory", "--max-range", "800", "--auto-zero", "600"],
            &[
                "--zero-temperature=-273.1435",
                "--zero-pressure=0.0187",
                "--zero-humidity=0",
            ],
        ),
        aloft("--zero-temperature"),
    ));
    let mut failures = Vec::new();
    for (case, args, expected) in &cases {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = run(home.path(), &args);
        let error = error_line(&output);
        let mut expected = expected.clone();
        if case == "compare" {
            expected = expected.replacen("Error: ", "Error: load 'cold': zeroing failed: ", 1);
        } else if case == "optimal-zero" {
            expected.push_str("; no zero distance from 90.0 to 840.0 m could be solved");
        }
        if error != expected {
            failures.push(format!("{case}: expected `{expected}`, got `{error}`"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} misnamed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn truing_and_estimation_keep_the_temperature_flag_when_air_fails_in_flight() {
    let home = Home::new();
    let mut failures = Vec::new();
    for command in [
        vec![
            "true-velocity",
            "--measured-drop",
            "30",
            "--range",
            "900",
            "-b",
            "0.243",
            "-m",
            "175",
            "-d",
            "0.308",
            "--offline",
            "--zero-distance",
            "700",
        ],
        vec![
            "estimate-bc",
            "-v",
            "2700",
            "-m",
            "168",
            "-d",
            "0.308",
            "--data",
            "800,2.5;900,10.2",
            "--zero-range",
            "700",
        ],
    ] {
        let mut args = command.clone();
        args.extend([
            "--temperature=-459.6583",
            "--pressure",
            "0.000552211",
            "--humidity",
            "0",
        ]);
        let output = run(home.path(), &args);
        let error = error_line(&output);
        if output.status.success()
            || !error.contains(
                "--temperature = -459.66 °F: the air aloft would fall below absolute zero",
            )
            || error.contains("temperature_c =")
        {
            failures.push(format!("{}: {}", command[0], stderr(&output)));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The figures the CHANGELOG, CLI_USAGE.md and the C ABI's docs quote for its density-altitude
/// conversions, the refusal found in flight and absurd readings, recomputed from the engine: the
/// station the conversions used to hand back is the unchecked inverse's
/// (`resolve_atmosphere_for_density_altitude`), and every refusal is run.
#[test]
fn the_docs_quote_the_engine_on_the_c_abi_the_flight_and_absurd_readings() {
    let (_, cold_c, cold_hpa) = resolve_atmosphere_for_density_altitude(44_300.0, None);
    let (_, colder_c, colder_hpa) = resolve_atmosphere_for_density_altitude(50_000.0, None);
    assert!(colder_hpa.is_nan());
    assert_quotes(
        "CHANGELOG.md",
        &[format!(
            "DA 44300 m gave {cold_c:.2} °C at {cold_hpa:.1e} hPa, and DA 50000 m {colder_c:.2} \
             °C beside a NaN pressure"
        )],
    );
    assert_quotes(
        "src/ffi.rs",
        &[format!(
            "DA 44300 m returned a station below absolute zero ({cold_c:.2} °C at {cold_hpa:.1e} \
             hPa), and DA 50000 m a temperature of {colder_c:.2} °C beside a NaN pressure"
        )],
    );
    let home = Home::new();
    home.profile(
        "p",
        "metric",
        &[
            "--temperature=-273.1435",
            "--pressure=0.0187",
            "--humidity=0",
        ],
        &[],
    );
    let error = |args: &[&str]| {
        let line = error_line(&run(home.path(), args));
        line.strip_prefix("Error: ")
            .unwrap_or_else(|| panic!("{args:?}: no refusal"))
            .to_string()
    };
    let flight = error(&[
        "-u",
        "metric",
        "come-ups",
        "--profile",
        "p",
        "--zero-distance",
        "600",
        "--end",
        "800",
    ]);
    let flight = flight
        .split_once(": ...")
        .map_or(flight.clone(), |(head, _)| head.to_string());
    let head = flight
        .split_once(" zero:")
        .map(|(head, _)| format!("{head} zero:"))
        .unwrap();
    let with = |flags: &[&str]| {
        let mut args: Vec<&str> = METRIC_LOAD[..2].to_vec();
        args.extend(["trajectory", "--max-range", "300"]);
        args.extend_from_slice(&METRIC_LOAD[2..]);
        args.extend_from_slice(flags);
        error(&args)
    };
    let negative = with(&["--pressure=-1e300"]);
    let tiny = with(&["--pressure=1e-300", "--humidity=50"]);
    let tiny = tiny
        .split_once(": no such air exists")
        .map(|(head, _)| head.split_once(" at ").unwrap().1.to_string())
        .unwrap();
    assert_eq!(tiny, "15 °C and 1e-300 hPa");
    assert_quotes("CHANGELOG.md", &[head, negative, format!("at {tiny}")]);
}

// --- a QNH as it was sent, and a temperature where it came from ---------------------------------

/// The metric CLI runs `qnh_and_temperature` checks a QNH on: every command that takes one, each
/// with the flags it needs to run, the station flags appended.
fn metric_qnh_commands() -> Vec<(&'static str, Vec<&'static str>)> {
    let load = ["-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82"];
    let with = |head: &[&'static str], tail: &[&'static str]| {
        let mut args = vec!["-u", "metric"];
        args.extend_from_slice(head);
        args.extend_from_slice(tail);
        args
    };
    vec![
        (
            "trajectory",
            with(&["trajectory", "--max-range", "300"], &load),
        ),
        ("zero", with(&["zero", "--target-distance", "200"], &load)),
        (
            "come-ups",
            with(
                &["come-ups", "--zero-distance", "100", "--end", "300"],
                &load,
            ),
        ),
        (
            "lead",
            with(&["lead", "--target-speed", "5", "--end", "300"], &load),
        ),
        ("mpbr", with(&["mpbr", "--vital-zone", "20"], &load)),
        (
            "plan-truing",
            with(
                &[
                    "plan-truing",
                    "--measurement-resolution",
                    "0.25",
                    "--candidate-ranges",
                    "300,500",
                ],
                &load,
            ),
        ),
        (
            "estimate-bc",
            with(
                &["estimate-bc", "--data", "100,0.5;200,2.0;300,5.0"],
                &["-v", "823", "-m", "10.9", "-d", "7.82"],
            ),
        ),
        (
            "true-velocity",
            with(
                &[
                    "true-velocity",
                    "--offline",
                    "--measured-drop",
                    "20",
                    "--range",
                    "600",
                ],
                &["--bc", "0.475", "--mass", "10.9", "--diameter", "7.82"],
            ),
        ),
        (
            "true-wind",
            with(
                &["true-wind", "--miss", "500:3", "--twist-rate", "254"],
                &[
                    "--velocity",
                    "823",
                    "--bc",
                    "0.475",
                    "--mass",
                    "10.9",
                    "--diameter",
                    "7.82",
                ],
            ),
        ),
    ]
}

/// A refused QNH is quoted as it was sent, with the station pressure it was reduced to, on every
/// surface: never the reduced station pressure as if it were the reading sent (`pressure_pa =
/// 649.071696485373 Pa` for a QNH of 1013.25 Pa, `--pressure = 2.99 hPa` for a typed 5, or
/// `--pressure = NaN hPa: must be finite` for a QNH with no station pressure). The pressure is
/// blamed when the QNH is itself no sea-level pressure (hPa sent as pascals), the altitude when
/// it is one, on the calculator commands too, and a QNH the altimeter-setting relation leaves no
/// station pressure for is refused the same way.
#[test]
fn a_qnh_refusal_quotes_the_qnh_as_sent_on_every_surface() {
    let mut failures = Vec::new();
    // The SI wire: `solve` and `atmosphere.density_altitude` alike, word for word.
    let wire = |pressure_pa: f64, altitude_m: f64| {
        json!({"altitude_m": altitude_m, "temperature_k": 288.15, "pressure_pa": pressure_pa,
               "pressure_reference": "qnh", "relative_humidity": 0.5})
    };
    for (atmosphere, path, quote, unit) in [
        (
            wire(1013.25, 1500.0),
            "$.atmosphere.pressure_pa",
            "; pressure_pa is a QNH of 1013.25 Pa, and reduced to the station altitude of 1500 m \
             it is a station pressure of ",
            "pressure",
        ),
        (
            wire(101_325.0, 36_000.0),
            "$.atmosphere.altitude_m",
            "; pressure_pa is a QNH of 101325 Pa, and reduced to the station altitude of 36000 m \
             it is a station pressure of ",
            "altitude",
        ),
    ] {
        let solved = solve(atmosphere.clone()).expect_err("solve refuses it");
        let answered = density_altitude(&atmosphere).expect_err("so does the command");
        for (command, error) in [("solve", &solved), ("density_altitude", &answered)] {
            let message = error["message"].as_str().unwrap_or_default();
            if error["path"] != path
                || !message.contains(quote)
                || !message.ends_with(&format!(": check the {unit}'s unit"))
                || message.contains("pressure_pa = ")
            {
                failures.push(format!("{command} {atmosphere}: {error}"));
            }
        }
        if solved != answered {
            failures.push(format!("{atmosphere}: {solved} != {answered}"));
        }
        let mut dry = atmosphere;
        dry["relative_humidity"] = json!(0.0);
        assert!(
            solve(dry.clone()).is_ok(),
            "the dry QNH station flies: {dry}"
        );
        assert!(
            density_altitude(&dry).is_ok(),
            "the dry QNH station is answered: {dry}"
        );
    }

    // The CLI, by flag, zero-day flag, saved profile field and `--location` column.
    let station = |pressure: &'static str, altitude: &'static str, humidity: &'static str| {
        [
            "--pressure",
            pressure,
            "--pressure-type",
            "qnh",
            "--altitude",
            altitude,
            "--temperature",
            "20",
            "--humidity",
            humidity,
        ]
    };
    let troposphere = |label: &str, qnh: &str, altitude: &str, station: &str, unit: &str| {
        (
            format!(
                "; {label} is a QNH of {qnh}, and reduced to the station altitude of {altitude} it \
                 is a station pressure of {station}, a pressure altitude of "
            ),
            format!(": check the {unit}'s unit"),
        )
    };
    let no_station = |label: &str, qnh: &str, altitude: &str, unit: &str| {
        (
            format!(
                "Error: {label} is a QNH of {qnh}, and reduced to the station altitude of \
                 {altitude} it leaves no station pressure: the altimeter-setting relation has none \
                 that high (check the {unit}'s unit)"
            ),
            String::new(),
        )
    };
    let mut cases: Vec<(String, Vec<String>, (String, String))> = Vec::new();
    for (command, args) in metric_qnh_commands() {
        for (flags, expected) in [
            (
                station("5", "1500", "100"),
                troposphere("--pressure", "5 hPa", "1500 m", "2.99 hPa", "pressure"),
            ),
            (
                station("1013.25", "36000", "100"),
                troposphere(
                    "--pressure",
                    "1013.25 hPa",
                    "36000 m",
                    "0.18 hPa",
                    "altitude",
                ),
            ),
            (
                station("5", "20000", "0"),
                no_station("--pressure", "5 hPa", "20000 m", "pressure"),
            ),
            (
                station("1013.25", "50000", "0"),
                no_station("--pressure", "1013.25 hPa", "50000 m", "altitude"),
            ),
        ] {
            let mut run_args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
            run_args.extend(flags.iter().map(|s| s.to_string()));
            cases.push((format!("{command} {flags:?}"), run_args, expected));
        }
    }
    let zero_day: Vec<String> = [
        "-u",
        "metric",
        "trajectory",
        "--max-range",
        "300",
        "-v",
        "823",
        "-b",
        "0.475",
        "-m",
        "10.9",
        "-d",
        "7.82",
        "--auto-zero",
        "100",
        "--zero-pressure",
        "5",
        "--zero-pressure-type",
        "qnh",
        "--zero-altitude",
        "1500",
        "--zero-temperature",
        "20",
        "--zero-humidity",
        "100",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    cases.push((
        "the zero day".into(),
        zero_day,
        troposphere("--zero-pressure", "5 hPa", "1500 m", "2.99 hPa", "pressure"),
    ));
    let home = Home::new();
    home.profile("q", "metric", &station("5", "1500", "0"), &[]);
    let stored = troposphere(
        "profile 'q' field pressure",
        "5 hPa",
        "1500 m",
        "2.99 hPa",
        "pressure",
    );
    for args in [
        vec![
            "-u",
            "metric",
            "come-ups",
            "--profile",
            "q",
            "--zero-distance",
            "100",
            "--end",
            "300",
            "--humidity",
            "100",
        ],
        vec![
            "-u",
            "metric",
            "lead",
            "--profile",
            "q",
            "--target-speed",
            "5",
            "--end",
            "300",
            "--humidity",
            "100",
        ],
        vec![
            "-u",
            "metric",
            "trajectory",
            "--saved-profile",
            "q",
            "--max-range",
            "300",
            "--humidity",
            "100",
        ],
    ] {
        cases.push((
            format!("{args:?}"),
            args.iter().map(|s| s.to_string()).collect(),
            stored.clone(),
        ));
    }
    let csv = home.file(
        "q.csv",
        "NAME,PRESSURE,ALTITUDE,TEMPERATURE,HUMIDITY\nQ,5,1500,20,100\n",
    );
    cases.push((
        "--location".into(),
        [
            "-u",
            "metric",
            "trajectory",
            "--max-range",
            "300",
            "-v",
            "823",
            "-b",
            "0.475",
            "-m",
            "10.9",
            "-d",
            "7.82",
            "--location",
            &csv,
            "--site",
            "Q",
            "--pressure-type",
            "qnh",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect(),
        troposphere(
            "--location CSV column PRESSURE",
            "5 hPa",
            "1500 m",
            "2.99 hPa",
            "pressure",
        ),
    ));
    let outputs = std::thread::scope(|scope| {
        let jobs: Vec<_> = cases
            .iter()
            .map(|(_, args, _)| {
                let home = home.path();
                scope.spawn(move || {
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    run(home, &args)
                })
            })
            .collect();
        jobs.into_iter()
            .map(|job| job.join().unwrap())
            .collect::<Vec<_>>()
    });
    for ((case, _, (quote, tail)), output) in cases.iter().zip(&outputs) {
        let error = error_line(output);
        let quoted = if tail.is_empty() {
            error == *quote
        } else {
            error.contains(quote.as_str()) && error.ends_with(tail.as_str())
        };
        if output.status.code() != Some(1) || !quoted || error.contains("pressure = ") {
            failures.push(format!(
                "{case}: expected `{quote}`…`{tail}`, got `{error}`"
            ));
        }
    }

    // The bridge's `profile.validate`, for the same stored QNH.
    let mut document: Value =
        serde_json::from_str(&std::fs::read_to_string(home.profile_path("q")).unwrap()).unwrap();
    document["humidity"] = json!(100.0);
    let validate = bridge("profile.validate", document);
    let warning = validate["result"]["warnings"][0]
        .as_str()
        .unwrap_or_default();
    if !warning.contains(&stored.0) {
        failures.push(format!("profile.validate: {validate}"));
    }
    assert!(
        failures.is_empty(),
        "{} misquoted:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Air that cannot exist, blamed on the temperature (outside -100 °C to 100 °C, where the
/// humid-air formulas hold), names the temperature where it came from, in the surface's units:
/// its flag, the zero day's flag, a saved profile's field, a `--location` column, the request's
/// member, or the density altitude that gave it, never only the humidity (`--humidity = 100 at
/// 214 °F and 29.92 inHg: no such air exists: ...` named no temperature, so a profile's or a
/// CSV's, or a density altitude, went unnamed, and the bridge blamed `$.temperature` with a
/// message that never mentioned it).
#[test]
fn impossible_air_blamed_on_the_temperature_names_where_it_came_from() {
    let home = Home::new();
    home.profile(
        "hot",
        "imperial",
        &["--temperature", "214", "--humidity", "0"],
        &[],
    );
    home.profile(
        "dap",
        "imperial",
        &["--humidity", "0"],
        &[("density_altitude", json!(145_000.0))],
    );
    let csv = home.file("hot.csv", "NAME,TEMP,HUMIDITY\nH,214,\n");
    let load = ["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];
    let imperial = |label: &str, value: &str| {
        format!(
            "; {label} = {value} is outside -148 °F to 212 °F, where the humid-air formulas hold: \
             check the temperature and its unit"
        )
    };
    let density_altitude = |label: &str| {
        format!(
            "; {label} gives a station temperature of -458.65 °F, outside -148 °F to 212 °F, \
             where the humid-air formulas hold: check the density altitude's unit"
        )
    };
    let with = |head: &[&str], tail: &[&str]| -> Vec<String> {
        head.iter()
            .chain(&load)
            .chain(tail)
            .map(|s| s.to_string())
            .collect()
    };
    let mut cases: Vec<(&str, Vec<String>, String)> = vec![
        (
            "a flag",
            with(
                &["trajectory", "--max-range", "300"],
                &["--temperature", "214", "--humidity", "100"],
            ),
            imperial("--temperature", "214 °F"),
        ),
        (
            "a metric flag",
            [
                "-u",
                "metric",
                "trajectory",
                "--max-range",
                "300",
                "-v",
                "823",
                "-b",
                "0.475",
                "-m",
                "10.9",
                "-d",
                "7.82",
                "--temperature",
                "101",
                "--humidity",
                "100",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            "; --temperature = 101 °C is outside -100 °C to 100 °C, where the humid-air formulas \
             hold: check the temperature and its unit"
                .to_string(),
        ),
        (
            "the zero day's flag",
            with(
                &["trajectory", "--max-range", "300", "--auto-zero", "100"],
                &["--zero-temperature", "214", "--zero-humidity", "100"],
            ),
            imperial("--zero-temperature", "214 °F"),
        ),
        (
            "a come-ups profile field",
            [
                "come-ups",
                "--profile",
                "hot",
                "--zero-distance",
                "100",
                "--end",
                "300",
                "--humidity",
                "100",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            imperial("profile 'hot' field temperature", "214 °F"),
        ),
        (
            "a trajectory profile field",
            [
                "trajectory",
                "--saved-profile",
                "hot",
                "--max-range",
                "300",
                "--humidity",
                "100",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            imperial("profile 'hot' field temperature", "214 °F"),
        ),
        (
            "a --location column",
            with(
                &["trajectory", "--max-range", "300"],
                &["--location", &csv, "--site", "H", "--humidity", "100"],
            ),
            imperial("--location CSV column TEMP", "214 °F"),
        ),
        (
            "a density-altitude flag",
            with(
                &["trajectory", "--max-range", "300"],
                &["--density-altitude=145000"],
            ),
            density_altitude("--density-altitude"),
        ),
        (
            "profile save's density-altitude flag",
            with(&["profile", "save", "dap2"], &["--density-altitude=145000"]),
            density_altitude("--density-altitude"),
        ),
        (
            "a stored density altitude",
            [
                "come-ups",
                "--profile",
                "dap",
                "--zero-distance",
                "100",
                "--end",
                "300",
                "--humidity",
                "50",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            density_altitude("profile 'dap' field density_altitude"),
        ),
    ];
    for command in [
        vec!["range-table", "--zero-distance", "100", "--end", "300"],
        vec!["wind-card", "--zero-distance", "100", "--end", "300"],
        vec!["lead", "--target-speed", "4", "--end", "300"],
        vec!["mpbr", "--vital-zone", "8"],
        vec![
            "compare",
            "--zero-distance",
            "100",
            "--end",
            "300",
            "--load",
            "B:G1:0.475:168:2600",
        ],
        vec!["mark-to-range", "--mark", "2"],
        vec![
            "bdc-match",
            "--mark-range",
            "1:300",
            "--mark-range",
            "2:400",
        ],
        vec!["optimal-zero", "--target", "200:10", "--target", "400:12"],
        vec![
            "adaptive-card",
            "--zero-distance",
            "100",
            "--start",
            "100",
            "--end",
            "300",
        ],
        vec![
            "plan-truing",
            "--measurement-resolution",
            "0.25",
            "--candidate-ranges",
            "300,500",
        ],
    ] {
        cases.push((
            command[0],
            command
                .iter()
                .chain(&["--profile", "hot", "--humidity", "100"])
                .map(|s| s.to_string())
                .collect(),
            imperial("profile 'hot' field temperature", "214 °F"),
        ));
    }
    let mut failures = Vec::new();
    for (case, args, tail) in &cases {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = run(home.path(), &args);
        let error = error_line(&output);
        if output.status.code() != Some(1)
            || !error.contains("no such air exists")
            || !error.ends_with(tail.as_str())
        {
            failures.push(format!("{case}: expected `…{tail}`, got `{error}`"));
        }
    }

    // The wire surfaces name their own member, in their own units.
    let mut card = json!({
        "muzzle_velocity": 2700.0, "ballistic_coefficient": 0.475, "mass": 168.0,
        "diameter": 0.308, "zero_distance": 100.0, "start": 100.0, "end": 300.0,
        "step": 100.0, "temperature": 214.0, "humidity": 100.0
    });
    let model = json!({
        "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
        "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
        "zero_distance_yd": 100.0, "sight_height_in": 2.0,
        "temperature_f": 214.0, "pressure_inhg": 29.92,
        "humidity_pct": 100.0, "altitude_ft": 0.0
    });
    let fit = json!({
        "model": model, "drop_unit": "mil",
        "observations": [{"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1}],
        "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                   "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
        "predictions": []
    });
    for (command, request, path, tail) in [
        (
            "card.range_table",
            card.clone(),
            "$.temperature",
            imperial("temperature", "214 °F"),
        ),
        (
            "true.fit",
            fit,
            "$.model.temperature_f",
            imperial("temperature_f", "214 °F"),
        ),
    ] {
        let out = bridge(command, request);
        let details = &out["error"]["details"]["error"];
        let message = details["message"].as_str().unwrap_or_default();
        if details["path"] != path || !message.ends_with(tail.as_str()) {
            failures.push(format!(
                "{command}: expected `…{tail}` at {path}, got {out}"
            ));
        }
    }
    card["temperature_unit"] = json!("celsius");
    card["temperature"] = json!(101.0);
    let out = bridge("card.come_ups", card);
    let message = out["error"]["details"]["error"]["message"]
        .as_str()
        .unwrap_or_default();
    if !message.ends_with(
        "; temperature = 101 °C is outside -100 °C to 100 °C, where the humid-air formulas hold: \
         check the temperature and its unit",
    ) {
        failures.push(format!("card.come_ups in celsius: {out}"));
    }
    let solved = solve(json!({"altitude_m": 0.0, "temperature_k": 400.0,
                              "relative_humidity": 1.0}))
    .expect_err("solve refuses it");
    if solved["path"] != "$.atmosphere.temperature_k"
        || !solved["message"].as_str().unwrap_or_default().ends_with(
            "; temperature_k = 400 K is outside 173.15 K to 373.15 K, where the humid-air formulas \
             hold: check the temperature and its unit",
        )
    {
        failures.push(format!("solve: {solved}"));
    }
    let mut document: Value =
        serde_json::from_str(&std::fs::read_to_string(home.profile_path("hot")).unwrap()).unwrap();
    document["humidity"] = json!(100.0);
    let validate = bridge("profile.validate", document);
    let warning = validate["result"]["warnings"][0]
        .as_str()
        .unwrap_or_default();
    if !warning.contains(&imperial("profile 'hot' field temperature", "214 °F")) {
        failures.push(format!("profile.validate: {validate}"));
    }
    assert!(
        failures.is_empty(),
        "{} misnamed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The QNH and temperature refusals the CHANGELOG, CLI_USAGE.md and the two wire documents
/// quote, each run.
#[test]
fn the_docs_quote_the_engine_on_a_qnh_and_a_blamed_temperature() {
    let home = Home::new();
    let error = |args: &[&str]| {
        let line = error_line(&run(home.path(), args));
        line.strip_prefix("Error: ")
            .unwrap_or_else(|| panic!("{args:?}: no refusal"))
            .to_string()
    };
    // What a refusal adds to "no such air exists: ...": its last clause.
    let blame = |message: &str| {
        message
            .rsplit_once("; ")
            .map(|(_, blame)| blame.to_string())
            .unwrap_or_else(|| panic!("no reading blamed: {message}"))
    };
    let metric = |flags: &[&str]| {
        let mut args = vec!["-u", "metric", "trajectory", "--max-range", "300"];
        args.extend_from_slice(&METRIC_LOAD[2..]);
        args.extend_from_slice(flags);
        error(&args)
    };
    let qnh = blame(&metric(&[
        "--pressure",
        "5",
        "--pressure-type",
        "qnh",
        "--altitude",
        "1500",
        "--temperature",
        "20",
        "--humidity",
        "100",
    ]));
    let no_station = metric(&[
        "--pressure",
        "5",
        "--pressure-type",
        "qnh",
        "--altitude",
        "20000",
        "--humidity",
        "0",
    ]);
    let hot = metric(&["--temperature", "101", "--humidity", "100"]);
    let density_altitude = blame(&error(&[
        "trajectory",
        "--max-range",
        "300",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--density-altitude",
        "145000",
    ]));
    home.profile(
        "hot",
        "imperial",
        &["--temperature", "214", "--humidity", "0"],
        &[],
    );
    let stored = blame(&error(&[
        "come-ups",
        "--profile",
        "hot",
        "--zero-distance",
        "100",
        "--end",
        "300",
        "--humidity",
        "100",
    ]));
    assert_quotes(
        "CHANGELOG.md",
        &[
            qnh.clone(),
            no_station.clone(),
            density_altitude.clone(),
            stored,
        ],
    );
    assert_quotes(
        "CLI_USAGE.md",
        &[qnh, no_station, density_altitude, format!("Error: {hot}")],
    );

    let temperature = blame(
        solve(json!({"altitude_m": 0.0, "temperature_k": 400.0, "relative_humidity": 1.0}))
            .expect_err("solve refuses it")["message"]
            .as_str()
            .unwrap(),
    );
    let wire_qnh = solve(json!({"altitude_m": 1500.0, "temperature_k": 288.15,
                                "pressure_pa": 1013.25, "pressure_reference": "qnh",
                                "relative_humidity": 0.5}))
    .expect_err("solve refuses it")["message"]
        .as_str()
        .unwrap()
        .to_string();
    let quoted = "pressure_pa is a QNH of 1013.25 Pa, and reduced to the station altitude of \
                  1500 m it is a station pressure of";
    assert!(wire_qnh.contains(quoted), "{wire_qnh}");
    for document in [
        "CHANGELOG.md",
        "docs/SOLVE_JSON_V1.md",
        "docs/ATMOSPHERE_DENSITY_ALTITUDE.md",
    ] {
        assert_quotes(document, &[temperature.clone(), quoted.to_string()]);
    }
}

#[test]
fn bridge_and_json_ffi_name_air_refused_in_flight_in_request_units() {
    let station = Station {
        temperature_c: -273.1435,
        pressure_hpa: Some(0.0187),
        humidity: 0.0,
        altitude_m: 0.0,
    };
    let mut card = card_request(&station);
    card["temperature"] = json!(station.temperature_c * 9.0 / 5.0 + 32.0);
    card["temperature_unit"] = json!("fahrenheit");
    card["zero_distance"] = json!(700.0);
    card["end"] = json!(900.0);
    card["wind_speeds"] = json!([5.0]);
    let mut model = truing_model(&station);
    model["zero_distance_yd"] = json!(700.0);
    let mut wind = wind_truing_request(&station);
    wind["zero_distance_yd"] = json!(700.0);
    let mut cases = vec![
        (
            "card.come_ups",
            card.clone(),
            "$.temperature",
            "temperature",
            "°F",
        ),
        (
            "card.range_table",
            card.clone(),
            "$.temperature",
            "temperature",
            "°F",
        ),
        (
            "card.wind",
            card.clone(),
            "$.temperature",
            "temperature",
            "°F",
        ),
        (
            "true.fit",
            json!({"model": model, "drop_unit": "mil",
            "observations": [{"range_yd": 800.0, "drop": 2.0, "sigma": 0.1}, {"range_yd": 900.0, "drop": 3.0, "sigma": 0.1}],
            "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                       "ballistic_coefficient": {"mean": 0.475, "sigma": 0.01}}, "predictions": []}),
            "$.model.temperature_f",
            "temperature_f",
            "°F",
        ),
        (
            "true.plan",
            json!({"model": model, "candidate_ranges_yd": [700.0, 900.0],
            "observation_count": 2, "minimum_separation_yd": 0.0, "measurement_sigma_1sd": 0.1,
            "drop_unit": "mil"}),
            "$.model.temperature_f",
            "temperature_f",
            "°F",
        ),
        (
            "true.dsf",
            json!({"model": model, "range_yd": 800.0, "observed_drop": 2.0, "drop_unit": "mil"}),
            "$.model.temperature_f",
            "temperature_f",
            "°F",
        ),
        ("true.wind", wind, "$.temperature_f", "temperature_f", "°F"),
    ];
    #[cfg(feature = "pdf")]
    {
        let mut pdf = card.clone();
        pdf["wind_speeds"] = json!([]);
        cases.push(("card.pdf", pdf, "$.temperature", "temperature", "°F"));
    }
    card["temperature"] = json!(station.temperature_c);
    card["temperature_unit"] = json!("celsius");
    card["altitude_unit"] = json!("meters");
    cases.push((
        "card.range_table",
        card,
        "$.temperature",
        "temperature",
        "°C",
    ));
    let mut failures = Vec::new();
    for (command, request, path, field, unit) in cases {
        let out = bridge(command, request.clone());
        let detail = &out["error"]["details"]["error"];
        let message = detail["message"].as_str().unwrap_or_default();
        if detail["path"] != path
            || detail["code"] != "invalid_value"
            || !message.starts_with(&format!("{field} = "))
            || !message.contains(unit)
            || !message.contains("air aloft would fall below absolute zero")
        {
            failures.push(format!("{command}: {out}"));
        }
        #[cfg(feature = "ffi")]
        {
            let envelope = json!({"api_version": 1, "command": command, "request": request});
            let text = std::ffi::CString::new(envelope.to_string()).unwrap();
            let raw =
                unsafe { ballistics_engine::bridge::ffi::ballistics_bridge_call(text.as_ptr()) };
            let ffi: Value =
                serde_json::from_str(unsafe { std::ffi::CStr::from_ptr(raw) }.to_str().unwrap())
                    .unwrap();
            unsafe { ballistics_engine::bridge::ffi::ballistics_bridge_free(raw) };
            assert_eq!(ffi, out, "{command}");
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn dsf_names_the_stored_temperature_when_the_flight_refuses_it() {
    let home = Home::new();
    home.profile(
        "cold",
        "metric",
        &[
            "--temperature=-273.1435",
            "--pressure=0.0187",
            "--humidity=0",
        ],
        &[("zero_distance", json!(600.0))],
    );
    let output = run(
        home.path(),
        &[
            "-u",
            "metric",
            "dsf",
            "--saved-profile",
            "cold",
            "--range",
            "800",
            "--observed-drop",
            "2mil",
        ],
    );
    let error = stderr(&output);
    assert!(!output.status.success(), "{error}");
    assert!(error.contains("profile 'cold' field temperature = -273.14 °C: the air aloft would fall below absolute zero"), "{error}");
}

#[test]
fn density_altitude_sources_name_the_derived_temperature_in_flight() {
    let home = Home::new();
    home.profile(
        "derived",
        "metric",
        &["--density-altitude", "44196", "--humidity", "0"],
        &[("zero_distance", Value::Null), ("auto_zero", Value::Null)],
    );
    let csv = home.file("derived.csv", "NAME,DA,HUMIDITY\nSITE,44196,0\n");
    let load = ["-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82"];
    let mut failures = Vec::new();
    for (extra, label) in [
        (
            vec!["--density-altitude", "44196", "--humidity", "0"],
            "--density-altitude",
        ),
        (
            vec!["--location", &csv, "--site", "SITE"],
            "--location CSV column DA",
        ),
        (
            vec!["--saved-profile", "derived"],
            "profile 'derived' field density_altitude",
        ),
    ] {
        let args = [
            vec![
                "-u",
                "metric",
                "trajectory",
                "--angle",
                "40",
                "--max-range",
                "20000",
            ],
            load.to_vec(),
            extra,
        ]
        .concat();
        let output = run(home.path(), &args);
        let error = error_line(&output);
        if output.status.success()
            || !error.contains(&format!("{label} gives a station temperature of "))
            || !error.contains("air aloft would fall below absolute zero")
            || error.contains("--temperature =")
        {
            failures.push(format!("{label}: {error}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn overriding_zero_pressure_preserves_an_inherited_density_altitude_temperature_source() {
    let home = Home::new();
    let output = run(
        home.path(),
        &[
            "-u",
            "metric",
            "trajectory",
            "-v",
            "823",
            "-b",
            "0.475",
            "-m",
            "10.9",
            "-d",
            "7.82",
            "--density-altitude",
            "44196",
            "--humidity",
            "0",
            "--auto-zero",
            "100",
            "--zero-pressure",
            "1e160",
        ],
    );
    let error = error_line(&output);
    assert!(!output.status.success(), "{error}");
    assert!(
        error.contains("--density-altitude gives a station temperature of "),
        "{error}"
    );
    assert!(error.contains("no such air exists"), "{error}");
}

#[test]
fn overriding_zero_altitude_preserves_an_inherited_density_altitude_pressure_source() {
    let home = Home::new();
    let output = run(
        home.path(),
        &[
            "-u",
            "metric",
            "trajectory",
            "-v",
            "823",
            "-b",
            "0.475",
            "-m",
            "10.9",
            "-d",
            "7.82",
            "--density-altitude",
            "25000",
            "--temperature",
            "40",
            "--humidity",
            "0",
            "--auto-zero",
            "100",
            "--zero-altitude",
            "0",
            "--zero-humidity",
            "100",
        ],
    );
    let error = error_line(&output);
    assert!(!output.status.success(), "{error}");
    assert!(
        error.contains("--density-altitude gives a station pressure of "),
        "{error}"
    );
    assert!(
        !error.contains("--zero-altitude gives a station pressure"),
        "{error}"
    );
    assert!(error.contains("no such air exists"), "{error}");
}
