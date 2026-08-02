//! MBA-1368: earth-fixed compass wind bearings re-referenced against the shot azimuth.
//!
//! Pins: `--wind-ref shooter` (and omission) is byte-identical to pre-1368 output; in
//! compass mode every wind direction — the single flag and each segment — is an
//! absolute bearing converted ONCE at the input boundary (`relative = bearing - shot
//! azimuth`, wind-FROM both sides); the ticket-mandated regression (wind FROM north +
//! shot due north = pure headwind, per the 0.19.0 sign convention); clock positions
//! rejected in compass mode (shooter-relative by definition); a missing
//! `--shot-direction` is a hard error naming the flag; Monte Carlo converts BEFORE
//! sampling (seeded WEZ sweep equality); and the solve-json v1 `wind_reference` field
//! (decode, conversion into the resolved echo, compass-without-azimuth error; omitted
//! and explicit "shooter"/"compass" solve identically and differ ONLY in whether
//! `resolved_request.wind.wind_reference` itself is echoed — the 0.33.0 Task 1
//! resolved-request-completeness echo).

use std::io::Write as _;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn private_home() -> std::path::PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-compass-wind-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(home: &std::path::Path, args: &[&str]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn run_ok(home: &std::path::Path, args: &[&str]) -> String {
    let out = run(home, args);
    assert!(
        out.status.success(),
        "command {:?} failed: {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

const LOAD: &[&str] = &["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];

fn trajectory_json(home: &std::path::Path, extra: &[&str]) -> String {
    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&["--max-range", "600", "--wind-speed", "10", "--full"]);
    args.extend_from_slice(extra);
    args.extend_from_slice(&["-o", "json"]);
    run_ok(home, &args)
}

/// Explicit `--wind-ref shooter` is byte-identical to omitting the flag, with and
/// without a shot direction on the line.
#[test]
fn explicit_shooter_mode_is_byte_identical_to_the_default() {
    let home = private_home();
    assert_eq!(
        trajectory_json(&home, &["--wind-direction", "90"]),
        trajectory_json(&home, &["--wind-direction", "90", "--wind-ref", "shooter"]),
    );
    assert_eq!(
        trajectory_json(&home, &["--wind-direction", "90", "--shot-direction", "220"]),
        trajectory_json(
            &home,
            &["--wind-direction", "90", "--shot-direction", "220", "--wind-ref", "shooter"],
        ),
    );
}

/// The ticket-mandated regression pin: wind FROM north (bearing 0) with the shot fired
/// due north (azimuth 0) is a PURE HEADWIND (relative 0, the 0.19.0 convention) — and
/// more generally bearing == azimuth collapses to a headwind at any azimuth.
#[test]
fn bearing_equal_to_azimuth_is_a_pure_headwind() {
    let home = private_home();

    let north_north = trajectory_json(
        &home,
        &["--wind-ref", "compass", "--shot-direction", "0", "--wind-direction", "0"],
    );
    let headwind = trajectory_json(&home, &["--wind-direction", "0", "--shot-direction", "0"]);
    assert_eq!(north_north, headwind);

    let sw_sw = trajectory_json(
        &home,
        &["--wind-ref", "compass", "--shot-direction", "220", "--wind-direction", "220"],
    );
    let headwind_220 =
        trajectory_json(&home, &["--wind-direction", "0", "--shot-direction", "220"]);
    assert_eq!(sw_sw, headwind_220);
}

/// Bearing 90 (wind from east) on a northbound shot = wind from the shooter's RIGHT
/// (relative 90): byte-identical to the shooter-relative 90 solve, and the drift sign
/// is pinned negative (positive x = right; wind FROM the right drifts x negative).
#[test]
fn bearing_90_on_a_northbound_shot_is_wind_from_the_right() {
    let home = private_home();

    let compass = trajectory_json(
        &home,
        &["--wind-ref", "compass", "--shot-direction", "0", "--wind-direction", "90"],
    );
    let shooter = trajectory_json(&home, &["--wind-direction", "90", "--shot-direction", "0"]);
    assert_eq!(compass, shooter);

    let v: serde_json::Value = serde_json::from_str(&compass).unwrap();
    let last = v["trajectory"].as_array().unwrap().last().unwrap().clone();
    let x = last["x"].as_f64().unwrap();
    assert!(x < 0.0, "wind from the shooter's right must drift x negative, got {x}");
}

/// Segments in compass mode are each re-referenced: a bearing-annotated segmented run
/// equals the same run with pre-converted shooter-relative angles.
#[test]
fn compass_segments_are_each_re_referenced() {
    let home = private_home();

    let compass = trajectory_json(
        &home,
        &[
            "--wind-ref", "compass", "--shot-direction", "90",
            "--wind-segment", "10:180:300", "--wind-segment", "5:45:600",
        ],
    );
    // bearing 180 - azimuth 90 = relative 90; bearing 45 - 90 = -45 -> 315.
    let shooter = trajectory_json(
        &home,
        &[
            "--shot-direction", "90",
            "--wind-segment", "10:90:300", "--wind-segment", "5:315:600",
        ],
    );
    assert_eq!(compass, shooter);
}

/// Clock positions are shooter-relative by definition: compass mode rejects them on
/// the flag AND inside segments, with the documented message.
#[test]
fn clock_positions_are_rejected_in_compass_mode() {
    let home = private_home();

    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&[
        "--max-range", "400", "--wind-speed", "10",
        "--wind-ref", "compass", "--shot-direction", "0", "--wind-direction", "3oc",
    ]);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("clock positions are shooter-relative by definition"),
        "{err}"
    );

    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&[
        "--max-range", "400",
        "--wind-ref", "compass", "--shot-direction", "0", "--wind-segment", "10:3oc:400",
    ]);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("clock positions are shooter-relative by definition"),
        "{err}"
    );
}

/// Compass mode without a shot azimuth is a hard error naming --shot-direction, on
/// trajectory and monte-carlo both (never a silent treat-as-shooter-relative).
#[test]
fn compass_without_shot_direction_is_a_hard_error()  {
    let home = private_home();

    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&["--max-range", "400", "--wind-ref", "compass", "--wind-direction", "90"]);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("--shot-direction"), "{err}");

    let mut args = vec!["monte-carlo"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&[
        "-n", "50", "--wind-speed", "10", "--wind-direction", "90", "--wind-ref", "compass",
    ]);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("--shot-direction"), "{err}");
}

/// Monte Carlo converts the base direction BEFORE sampling: the internally seeded WEZ
/// sweep is byte-identical between compass mode and the pre-converted shooter run —
/// impossible if the conversion happened after (or inside) the dispersion draws.
#[test]
fn monte_carlo_converts_before_sampling_seeded_wez_equality() {
    let home = private_home();

    let wez = |extra: &[&str]| -> String {
        let mut args = vec!["monte-carlo"];
        args.extend_from_slice(LOAD);
        args.extend_from_slice(&[
            "-n", "300", "--wind-speed", "10", "--wind-std", "2",
            "--wind-direction-std", "15", "--velocity-std", "10",
            "--wez", "--target-size", "2",
            "--wez-start", "200", "--wez-end", "600", "--wez-step", "200", "-o", "full",
        ]);
        args.extend_from_slice(extra);
        run_ok(&home, &args)
    };

    let compass = wez(&["--wind-ref", "compass", "--shot-direction", "30", "--wind-direction", "120"]);
    let shooter = wez(&["--wind-direction", "90"]);
    assert_eq!(compass, shooter);
}

// ---------------------------------------------------------------------------
// solve-json v1 wire surface
// ---------------------------------------------------------------------------

fn solve_json(request: &serde_json::Value) -> serde_json::Value {
    let mut child = Command::new(BIN)
        .args(["solve-json"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn solve-json");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(request.to_string().as_bytes())
        .unwrap();
    let out = child.wait_with_output().expect("solve-json output");
    serde_json::from_slice(&out.stdout).expect("solve-json emits one JSON envelope")
}

/// Strip the Task 1 (0.33.0) `resolved_request.wind.wind_reference` echo so two
/// responses that differ only in whether the caller explicitly named a wind reference
/// mode can still be compared byte-for-byte on everything else.
fn without_wind_reference_echo(mut response: serde_json::Value) -> serde_json::Value {
    if let Some(wind) = response["resolved_request"]["wind"].as_object_mut() {
        wind.remove("wind_reference");
    }
    response
}

fn base_request() -> serde_json::Value {
    serde_json::json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.031,
            "drag_model": "G7", "ballistic_coefficient": 0.243
        },
        "rifle": {
            "muzzle_velocity_mps": 823.0, "sight_height_m": 0.05, "muzzle_height_m": 1.2,
            "twist_rate_m_per_turn": 0.2032, "twist_direction": "right"
        },
        "shot": { "max_range_m": 600.0 },
        "atmosphere": {
            "altitude_m": 0.0, "temperature_k": 293.15,
            "pressure_pa": 100000.0, "relative_humidity": 0.4
        },
        "wind": {
            "speed_mps": 4.47,
            "direction_from_rad": std::f64::consts::FRAC_PI_2,
            "vertical_speed_mps": 0.0
        },
        "solver": { "method": "rk45", "time_step_s": 0.001 },
        "effects": { "magnus": false, "coriolis": false, "enhanced_spin_drift": false },
        "sampling": { "interval_m": 100.0 }
    })
}

/// wind_reference decodes (not unknown_field); compass converts the constant direction
/// into the resolved echo AND the physics (byte-identity with the pre-converted
/// shooter request); explicit "shooter"/"compass" and omitted solve identically and
/// differ ONLY in whether `resolved_request.wind.wind_reference` itself is echoed (the
/// 0.33.0 Task 1 resolved-request-completeness echo).
#[test]
fn solve_json_wind_reference_converts_at_the_wire() {
    // Omitted == explicit "shooter" apart from the Task 1 echo: strip
    // `resolved_request.wind.wind_reference` and the two envelopes are byte-for-byte
    // identical, because "shooter" is the behavioral default.
    let omitted = solve_json(&base_request());
    let mut shooter = base_request();
    shooter["wind"]["wind_reference"] = serde_json::json!("shooter");
    let shooter_response = solve_json(&shooter);
    assert!(omitted["resolved_request"]["wind"]
        .get("wind_reference")
        .is_none());
    assert_eq!(
        shooter_response["resolved_request"]["wind"]["wind_reference"],
        "shooter"
    );
    assert_eq!(omitted, without_wind_reference_echo(shooter_response));
    assert_eq!(omitted["status"], "ok");

    // Compass: bearing PI (from south) on a shot fired at azimuth PI/2 (east) =
    // relative PI/2 (from the shooter's right) — identical to the shooter request
    // with direction PI/2 and the same shot azimuth, apart from the same echo
    // difference as above.
    let mut compass = base_request();
    compass["shot"]["shot_azimuth_rad"] = serde_json::json!(std::f64::consts::FRAC_PI_2);
    compass["wind"]["wind_reference"] = serde_json::json!("compass");
    compass["wind"]["direction_from_rad"] = serde_json::json!(std::f64::consts::PI);
    let compass_response = solve_json(&compass);

    let mut reference = base_request();
    reference["shot"]["shot_azimuth_rad"] = serde_json::json!(std::f64::consts::FRAC_PI_2);
    let reference_response = solve_json(&reference);

    assert_eq!(compass_response["status"], "ok", "{compass_response}");
    assert!(reference_response["resolved_request"]["wind"]
        .get("wind_reference")
        .is_none());
    assert_eq!(
        compass_response["resolved_request"]["wind"]["wind_reference"],
        "compass"
    );
    assert_eq!(
        without_wind_reference_echo(compass_response.clone()),
        reference_response
    );

    // The resolved echo carries the CONVERTED shooter-relative direction.
    let resolved_dir = compass_response["resolved_request"]["wind"]["direction_from_rad"]
        .as_f64()
        .expect("resolved constant wind direction");
    assert!((resolved_dir - std::f64::consts::FRAC_PI_2).abs() < 1e-12, "{resolved_dir}");

    // Segments convert per segment too.
    let mut segmented = base_request();
    segmented["shot"]["shot_azimuth_rad"] = serde_json::json!(std::f64::consts::FRAC_PI_2);
    segmented["wind"] = serde_json::json!({
        "wind_reference": "compass",
        "segments": [
            { "until_distance_m": 300.0, "speed_mps": 4.47, "direction_from_rad": std::f64::consts::PI },
            { "until_distance_m": 700.0, "speed_mps": 2.0, "direction_from_rad": std::f64::consts::FRAC_PI_2 }
        ]
    });
    let seg_response = solve_json(&segmented);
    assert_eq!(seg_response["status"], "ok", "{seg_response}");
    let segs = seg_response["resolved_request"]["wind"]["segments"].as_array().unwrap();
    assert!((segs[0]["direction_from_rad"].as_f64().unwrap() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    // PI/2 - PI/2 = 0 (a pure headwind segment).
    assert!(segs[1]["direction_from_rad"].as_f64().unwrap().abs() < 1e-12);
}

/// Compass without shot.shot_azimuth_rad is a structured error at
/// $.wind.wind_reference; an unknown reference value is rejected naming the field.
#[test]
fn solve_json_compass_requires_shot_azimuth_and_validates_values() {
    let mut compass = base_request();
    compass["wind"]["wind_reference"] = serde_json::json!("compass");
    let response = solve_json(&compass);
    assert_eq!(response["status"], "error", "{response}");
    assert_eq!(response["error"]["path"], "$.wind.wind_reference", "{response}");
    assert!(
        response["error"]["message"]
            .as_str()
            .unwrap()
            .contains("shot.shot_azimuth_rad"),
        "{response}"
    );

    // An explicitly PROVIDED azimuth of 0.0 is a real bearing (due north), not "unset".
    let mut north = base_request();
    north["shot"]["shot_azimuth_rad"] = serde_json::json!(0.0);
    north["wind"]["wind_reference"] = serde_json::json!("compass");
    north["wind"]["direction_from_rad"] = serde_json::json!(0.0);
    let response = solve_json(&north);
    assert_eq!(response["status"], "ok", "{response}");

    let mut bogus = base_request();
    bogus["wind"]["wind_reference"] = serde_json::json!("magnetic");
    let response = solve_json(&bogus);
    assert_eq!(response["status"], "error");
    assert_eq!(response["error"]["path"], "$.wind.wind_reference", "{response}");
}
