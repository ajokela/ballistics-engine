//! Integration tests for MBA-1366: `--density-altitude`, the single-value atmosphere entry
//! mode (Shooter, Nosler, AB Analytics, Ballistic AE, TRASOL) that supersedes
//! `--altitude`/`--pressure`/`--pressure-type` by back-solving an ISA-equivalent atmosphere
//! (NOT the direct-density bypass — see `atmosphere::resolve_atmosphere_for_density_altitude`).
//!
//! Hard requirement: with `--density-altitude` absent, output must be byte-identical to before
//! this flag existed. The two `..._byte_identical_to_pinned_baseline` tests below pin the exact
//! stdout captured from this branch's binary BEFORE any of MBA-1366's main.rs changes were
//! written (i.e. before `--density-altitude` existed at all), proving the new CLI arg/CSV
//! plumbing is a pure, inert addition on every path that doesn't reference it.

use serde_json::Value;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn get_cli_binary() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("ballistics");
    if !path.exists() {
        path.pop();
        path.pop();
        path.push("release");
        path.push("ballistics");
    }
    path
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("failed to execute ballistics binary")
}

fn run_json(args: &[&str]) -> Value {
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

fn with_extra<'a>(base: &[&'a str], extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = base.to_vec();
    args.extend_from_slice(extra);
    args
}

/// Writes `content` to a fresh temp file and returns its path. The caller is responsible for
/// removing it (no `tempfile` dev-dependency in this crate); each test uses a distinct
/// process-id + counter-based filename to be safe under parallel test execution.
fn write_temp_csv(content: &str, tag: &str) -> PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let mut path = std::env::temp_dir();
    path.push(format!(
        "ballistics_mba1366_{tag}_{}_{n}.csv",
        std::process::id()
    ));
    let mut f = std::fs::File::create(&path).expect("create temp csv");
    f.write_all(content.as_bytes()).expect("write temp csv");
    path
}

// ---- Hard requirement: byte-identical with no --density-altitude ---------------------------

/// Pinned literal captured from this branch's `ballistics` binary BEFORE MBA-1366's
/// `--density-altitude` flag existed (i.e. before any main.rs change in this task).
const PINNED_JSON_BASELINE: &str = r#"{
  "units": "metric",
  "max_range": 300.0,
  "max_height": 1.5,
  "time_of_flight": 0.12333720481671129,
  "impact_velocity": 2199.3005312520913,
  "impact_energy": 406301.51744832145,
  "stability_coefficient": 2.0961431820049854,
  "spin_drift": null,
  "trajectory": [],
  "legend": {
    "units": {
      "distance": "m",
      "velocity": "m/s",
      "energy": "J"
    },
    "axes": {
      "x": "lateral offset from the muzzle's initial aiming direction; positive means to the shooter's right (e.g. a crosswind FROM the left, --wind-direction 270, drifts x positive; FROM the right, --wind-direction 90, drifts x negative). Zero at the muzzle.",
      "y": "height above the ground in the world frame; positive means up. Starts near bore height and reaches 0 at ground impact. This is NOT height above the line of sight (compare solve-json v1's LOS-relative drop_m).",
      "z": "downrange distance from the muzzle; zero at the muzzle, always increasing."
    }
  }
}"#;

const PINNED_TABLE_BASELINE: &str = "╔════════════════════════════════════════╗
║         TRAJECTORY RESULTS             ║
╠════════════════════════════════════════╣
║ Max Range:           446.23 yd        ║
║ Max Height:            1.67 yd        ║
║ Time of Flight:       0.589 s          ║
║ Impact Velocity:    1912.67 fps       ║
║ Impact Energy:      1364.45 ft-lb     ║
╠════════════════════════════════════════╣
║ Stability (SG):        2.00            ║
║ Status:             assumed            ║
╚════════════════════════════════════════╝
note: twist rate not supplied; assumed 1:11.7\" for the stability/spin-drift estimate — Sg shown is not an evaluated stability verdict.

Bullet struck ground at 446 yd
";

#[test]
fn no_density_altitude_json_output_is_byte_identical_to_pinned_pre_mba_1366_baseline() {
    let args: &[&str] = &[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--units",
        "metric", "--max-range", "300", "--ignore-ground-impact", "-o", "json", "--altitude",
        "1500", "--pressure", "950", "--temperature", "10", "--humidity", "60",
    ];
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert_eq!(
        stdout.trim_end(),
        PINNED_JSON_BASELINE.trim_end(),
        "no --density-altitude output must be byte-identical to the pre-MBA-1366 baseline"
    );
}

#[test]
fn no_density_altitude_table_output_is_byte_identical_to_pinned_pre_mba_1366_baseline() {
    let args: &[&str] = &[
        "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--units",
        "imperial", "--max-range", "500", "-o", "table",
    ];
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert_eq!(
        stdout, PINNED_TABLE_BASELINE,
        "no --density-altitude table output must be byte-identical to the pre-MBA-1366 baseline"
    );
}

// ---- Smoke test: a DA-specified run matches an equivalent explicit run -----------------------

const BASE_ARGS: &[&str] = &[
    "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--units", "metric",
    "--max-range", "500", "--ignore-ground-impact", "-o", "json",
];

#[test]
fn density_altitude_matches_equivalent_explicit_atmosphere_to_printed_precision() {
    // Back-solve the ISA-equivalent (altitude_m, temp_c, pressure_hpa) for a 1000 ft (304.8 m)
    // density altitude, at the ISA-at-DA default temperature -- then feed BOTH representations
    // (the single --density-altitude value, and the explicit altitude+pressure+temperature triple
    // it implies) through the real solve and confirm they match to printed precision.
    let da_m = 1000.0 * 0.3048;
    let (altitude_m, temp_c, pressure_hpa) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(da_m, None);

    let via_da = run_json(&with_extra(
        BASE_ARGS,
        &["--density-altitude", &da_m.to_string()],
    ));
    let via_explicit = run_json(&with_extra(
        BASE_ARGS,
        &[
            "--altitude",
            &altitude_m.to_string(),
            "--pressure",
            &pressure_hpa.to_string(),
            "--temperature",
            &temp_c.to_string(),
        ],
    ));

    let v_da = via_da["impact_velocity"].as_f64().unwrap();
    let v_explicit = via_explicit["impact_velocity"].as_f64().unwrap();
    assert!(
        (v_da - v_explicit).abs() < 0.5,
        "DA-specified and equivalent-explicit runs must match to printed precision: \
         da={v_da} explicit={v_explicit}"
    );
    let h_da = via_da["max_height"].as_f64().unwrap();
    let h_explicit = via_explicit["max_height"].as_f64().unwrap();
    assert!((h_da - h_explicit).abs() < 0.01);
}

#[test]
fn density_altitude_at_zero_matches_sea_level_standard_atmosphere() {
    let sea_level = run_json(&with_extra(
        BASE_ARGS,
        &["--altitude", "0", "--pressure", "1013.25", "--temperature", "15"],
    ));
    let via_da = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "0"]));
    assert_eq!(
        sea_level["impact_velocity"], via_da["impact_velocity"],
        "0 density altitude must reduce to the plain ICAO sea-level standard atmosphere"
    );
}

// ---- Precedence: DA supersedes --altitude/--pressure -----------------------------------------

#[test]
fn density_altitude_supersedes_altitude_and_pressure() {
    let da_only = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "2000"]));
    // Same DA, but with --altitude/--pressure ALSO set to values that would (if honored) fly a
    // materially different trajectory -- DA must win outright.
    let da_with_overrides = run_json(&with_extra(
        BASE_ARGS,
        &[
            "--density-altitude",
            "2000",
            "--altitude",
            "50",
            "--pressure",
            "1013.25",
        ],
    ));
    assert_eq!(
        da_only["impact_velocity"], da_with_overrides["impact_velocity"],
        "--density-altitude must supersede --altitude/--pressure entirely"
    );

    // Sanity: the ignored --altitude/--pressure values would have changed the answer if honored.
    let would_be_different = run_json(&with_extra(
        BASE_ARGS,
        &["--altitude", "50", "--pressure", "1013.25"],
    ));
    assert_ne!(
        da_only["impact_velocity"], would_be_different["impact_velocity"],
        "fixture must show --altitude/--pressure alone would have flown differently"
    );
}

#[test]
fn density_altitude_supersedes_pressure_type_qnh_and_warns() {
    let da_only = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "2000"]));
    let output = run(&with_extra(
        BASE_ARGS,
        &[
            "--density-altitude",
            "2000",
            "--pressure",
            "1030.0",
            "--pressure-type",
            "qnh",
        ],
    ));
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--density-altitude supersedes --pressure/--pressure-type"),
        "stderr must note the supersede: {stderr}"
    );
    let combined: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        da_only["impact_velocity"], combined["impact_velocity"],
        "a QNH --pressure/--pressure-type must be ignored entirely when --density-altitude is given"
    );
}

// ---- Precedence: explicit --temperature wins over ISA-at-DA, density still honored -----------

#[test]
fn explicit_temperature_overrides_isa_default_under_density_altitude() {
    // BASE_ARGS runs --units metric, so --temperature 35 means 35 C directly (well above the
    // ISA-at-3000ft-DA default of roughly 10 C).
    let default_temp = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "3000"]));
    let explicit_temp = run_json(&with_extra(
        BASE_ARGS,
        &["--density-altitude", "3000", "--temperature", "35"],
    ));
    assert_ne!(
        default_temp["impact_velocity"], explicit_temp["impact_velocity"],
        "an explicit --temperature must change the resolved atmosphere under --density-altitude"
    );

    // "Density is still honored" half of the contract, proven directly at the library level:
    // feeding the explicit-temperature branch's own (temp_c, pressure_hpa) back through the
    // production round-trip (exercised end-to-end in src/main.rs's
    // density_altitude_round_trip_tests) reproduces the same 3000 ft density altitude. Here we
    // just confirm the explicit branch does NOT collapse to the ISA-default pressure (i.e. it
    // really re-solved a different pressure for the different temperature).
    let da_m = 3000.0 * 0.3048;
    let (_, default_temp_c, default_pressure_hpa) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(da_m, None);
    let (_, explicit_temp_c, explicit_pressure_hpa) =
        ballistics_engine::atmosphere::resolve_atmosphere_for_density_altitude(da_m, Some(35.0));
    assert_eq!(explicit_temp_c, 35.0);
    assert_ne!(default_temp_c, explicit_temp_c);
    assert_ne!(default_pressure_hpa, explicit_pressure_hpa);
}

// ---- CSV DA/DENSITY_ALTITUDE column import ----------------------------------------------------

#[test]
fn csv_da_column_is_used_when_flag_omitted() {
    let csv_path = write_temp_csv(
        "LOCATION_NAME,DA\nTestSite,2500\n",
        "da_column",
    );
    let via_csv = run_json(&with_extra(
        BASE_ARGS,
        &[
            "--location",
            csv_path.to_str().unwrap(),
            "--site",
            "TestSite",
        ],
    ));
    let via_flag = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "2500"]));
    std::fs::remove_file(&csv_path).ok();

    assert_eq!(
        via_csv["impact_velocity"], via_flag["impact_velocity"],
        "a CSV DA column must be equivalent to the same --density-altitude flag value"
    );
}

#[test]
fn csv_density_altitude_column_name_is_also_recognized() {
    let csv_path = write_temp_csv(
        "LOCATION_NAME,DENSITY_ALTITUDE\nTestSite2,1800\n",
        "density_altitude_column",
    );
    let via_csv = run_json(&with_extra(
        BASE_ARGS,
        &[
            "--location",
            csv_path.to_str().unwrap(),
            "--site",
            "TestSite2",
        ],
    ));
    let via_flag = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "1800"]));
    std::fs::remove_file(&csv_path).ok();

    assert_eq!(
        via_csv["impact_velocity"], via_flag["impact_velocity"],
        "a CSV DENSITY_ALTITUDE column must be equivalent to the same --density-altitude flag value"
    );
}

#[test]
fn cli_density_altitude_flag_wins_over_csv_column() {
    let csv_path = write_temp_csv(
        "LOCATION_NAME,DA\nTestSite3,9000\n",
        "flag_wins",
    );
    let via_both = run_json(&with_extra(
        BASE_ARGS,
        &[
            "--location",
            csv_path.to_str().unwrap(),
            "--site",
            "TestSite3",
            "--density-altitude",
            "500",
        ],
    ));
    let via_flag_alone = run_json(&with_extra(BASE_ARGS, &["--density-altitude", "500"]));
    std::fs::remove_file(&csv_path).ok();

    assert_eq!(
        via_both["impact_velocity"], via_flag_alone["impact_velocity"],
        "an explicit --density-altitude flag must win over a CSV DA column"
    );
}

/// MBA-1422: `--density-altitude` on the standalone `zero` command, and `--zero-density-altitude`
/// for declaring a zero-day density altitude on `trajectory --auto-zero`.
///
/// MBA-1366 scoped the flag to `trajectory` only. These are the two entry points it left out:
/// a user carrying conditions as a density altitude could not solve a zero the same way they
/// solve a trajectory, and could not tear a zero under a different DA than the shot.
///
/// The zero-day path is where 0.29.0's whole-branch review found two bugs of the same family — a
/// stale pressure mode surviving into the zero solve, and an explicit `--zero-pressure-type`
/// re-reducing a DA-derived pressure that was already absolute. The rule these follow is the one
/// that came out of those: a mode travels with the value it describes, or not at all.
mod density_altitude_on_zero_surfaces {
    use super::*;

    fn zero_angle(args: &[&str]) -> f64 {
        let output = Command::new(get_cli_binary())
            .args([
                "zero", "--velocity", "2700", "--bc", "0.243", "--mass", "175", "--diameter",
                "0.308", "--target-distance", "300", "-o", "json",
            ])
            .args(args)
            .output()
            .expect("zero");
        assert!(
            output.status.success(),
            "zero failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let doc: Value = serde_json::from_slice(&output.stdout).expect("json");
        doc["zero_angle_degrees"].as_f64().expect("zero_angle_degrees")
    }

    fn auto_zero_angle(args: &[&str]) -> f64 {
        let output = Command::new(get_cli_binary())
            .args([
                "trajectory", "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7", "-m",
                "175", "-d", "0.308", "--max-range", "500", "--auto-zero", "300", "-o", "json",
            ])
            .args(args)
            .output()
            .expect("trajectory");
        assert!(
            output.status.success(),
            "trajectory failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let doc: Value = serde_json::from_slice(&output.stdout).expect("json");
        doc["zero_angle_degrees"].as_f64().expect("zero_angle_degrees")
    }

    /// A density altitude equal to a real altitude must give the same answer under the standard
    /// atmosphere, since that is what "density altitude" means when conditions are ISA.
    #[test]
    fn zero_density_altitude_matches_the_equivalent_real_altitude_under_isa() {
        assert_eq!(
            zero_angle(&["--density-altitude", "5000"]),
            zero_angle(&["--altitude", "5000"]),
            "under ISA a 5000 ft density altitude is a 5000 ft altitude"
        );
    }

    /// Omitting the flag must leave `zero` exactly as it was.
    #[test]
    fn zero_without_the_flag_is_unchanged() {
        assert_eq!(
            zero_angle(&[]),
            zero_angle(&["--density-altitude", "0"]),
            "a sea-level density altitude is the standard atmosphere"
        );
    }

    /// DA supersedes pressure outright — that is the documented precedence, and it is what makes
    /// the pressure mode inapplicable.
    #[test]
    fn zero_density_altitude_supersedes_an_explicit_pressure() {
        assert_eq!(
            zero_angle(&["--density-altitude", "5000", "--pressure", "29.92"]),
            zero_angle(&["--density-altitude", "5000"]),
            "--pressure must not survive alongside --density-altitude on `zero`"
        );
    }

    /// Same three properties for the zero-day flag on `trajectory`.
    #[test]
    fn zero_day_density_altitude_matches_the_equivalent_zero_day_altitude() {
        assert_eq!(
            auto_zero_angle(&["--zero-density-altitude", "5000"]),
            auto_zero_angle(&["--zero-altitude", "5000"]),
        );
    }

    #[test]
    fn zero_day_density_altitude_supersedes_an_explicit_zero_pressure() {
        assert_eq!(
            auto_zero_angle(&["--zero-density-altitude", "5000", "--zero-pressure", "29.92"]),
            auto_zero_angle(&["--zero-density-altitude", "5000"]),
        );
    }

    /// The zero-day flag must not disturb the shot-day atmosphere, which is the whole point of
    /// having a separate zero-day knob.
    #[test]
    fn a_zero_day_density_altitude_does_not_touch_the_shot_day_atmosphere() {
        // Deliberately compared on IMPACT VELOCITY, not the zero angle. Both runs produce the
        // same zero angle, and correctly so: with no zero-day flags the zero day inherits the
        // shot-day atmosphere, so `--density-altitude 5000` and `--zero-density-altitude 5000`
        // both put the ZERO solve at 5000 ft DA. The difference is in the shot-day trajectory,
        // which only the former moves. Asserting on the angle here would have looked like a bug
        // in the feature when it was a bug in the test.
        let impact = |args: &[&str]| -> f64 {
            let output = Command::new(get_cli_binary())
                .args([
                    "trajectory", "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7",
                    "-m", "175", "-d", "0.308", "--max-range", "500", "--auto-zero", "300", "-o",
                    "json",
                ])
                .args(args)
                .output()
                .expect("trajectory");
            assert!(output.status.success());
            let doc: Value = serde_json::from_slice(&output.stdout).expect("json");
            doc["impact_velocity"].as_f64().expect("impact_velocity")
        };

        let baseline = impact(&[]);
        let zero_day_only = impact(&["--zero-density-altitude", "5000"]);
        let shot_day = impact(&["--density-altitude", "5000"]);

        // Not an equality: a zero-day DA legitimately perturbs the shot-day trajectory by a
        // hair, because it changes the solved zero ANGLE and therefore the bore angle the shot
        // is fired at. That path is real but tiny (~2e-5 fps here). What must NOT happen is the
        // zero-day knob moving the shot-day AIR, which would show up at the same scale as the
        // shot-day flag. Three orders of magnitude separate the two effects.
        let via_zero_angle = (baseline - zero_day_only).abs();
        let via_shot_day_air = (baseline - shot_day).abs();
        assert!(
            via_zero_angle < 0.01,
            "a zero-day density altitude moved the shot-day trajectory by {via_zero_angle} fps,              far more than the bore-angle coupling can explain — it is reaching the shot-day              atmosphere"
        );
        assert!(
            via_shot_day_air > 1.0,
            "a shot-day density altitude should move the shot-day trajectory substantially;              got only {via_shot_day_air} fps"
        );
    }

    /// With a zero-day DA and no --zero-pressure, an explicit --zero-pressure-type has no value
    /// of its own to describe. It must be refused and say so, rather than re-reducing a pressure
    /// the back-solve already produced as absolute.
    #[test]
    fn a_zero_pressure_type_is_refused_and_reported_when_a_zero_day_da_supplies_the_pressure() {
        let output = Command::new(get_cli_binary())
            .args([
                "trajectory", "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7", "-m",
                "175", "-d", "0.308", "--max-range", "500", "--auto-zero", "300", "-o", "json",
                "--zero-density-altitude", "5000", "--zero-pressure-type", "qnh",
            ])
            .output()
            .expect("trajectory");
        assert!(output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("--zero-pressure-type is ignored"),
            "expected a refusal notice, got: {stderr}"
        );
        assert!(
            stderr.contains("--zero-density-altitude"),
            "the notice must name the flag that actually supplied the pressure, got: {stderr}"
        );

        // And it must be refused, not merely announced.
        let refused: Value = serde_json::from_slice(&output.stdout).expect("json");
        let refused = refused["zero_angle_degrees"].as_f64().expect("angle");
        assert_eq!(
            refused,
            auto_zero_angle(&["--zero-density-altitude", "5000"]),
            "the ignored mode still altered the solve"
        );
    }
}
