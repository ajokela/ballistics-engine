//! MBA-1349: the seven ticket-mandated acceptance criteria for robust hold corridors,
//! plus the `hold-corridor` CLI surface.
//!
//! Criteria 2-7 are proved at library level in `src/wind_scenarios.rs`'s own test module
//! (they are properties of the corridor math and need no process boundary). This file
//! carries the one that genuinely needs an INDEPENDENT construction — criterion 1, "each
//! scenario row matches a direct trajectory solve with the same segmented wind" — because
//! proving it with the same code path that produced the row would be circular. Here the
//! reference trajectory is built from the public `TrajectorySolver` API by hand.
//!
//! It then pins the CLI: both output shapes, the versioned JSON, order insensitivity
//! through the file, and every rejection.

use std::process::Command;

use ballistics_engine::cli_api::{
    calculate_zero_angle_with_conditions, AtmosphericConditions, BallisticInputs, TrajectorySolver,
};
use ballistics_engine::wind::parse_wind_segment_str;
use ballistics_engine::wind_scenarios::{
    solve_robust_hold, CorridorLoad, NamedWindScenario, RobustHoldReportV1, RobustHoldRequest,
    TargetSpec, WindScenarioSetV1,
};
use ballistics_engine::{DragModel, WindConditions};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn load() -> CorridorLoad {
    CorridorLoad {
        muzzle_velocity_mps: 823.0,
        bc: 0.475,
        drag_model: DragModel::G1,
        mass_kg: 0.010_886,
        diameter_m: 0.007_823,
        bullet_length_m: 0.031,
        zero_distance_m: 91.44,
        sight_height_m: 0.0508,
        temperature_c: 15.0,
        pressure_hpa: 1013.25,
        humidity_pct: 50.0,
        altitude_m: 0.0,
    }
}

fn scenario(name: &str, tokens: &[&str]) -> NamedWindScenario {
    NamedWindScenario {
        name: name.to_string(),
        segments: tokens
            .iter()
            .map(|token| parse_wind_segment_str(token, true).unwrap())
            .collect(),
    }
}

/// A reference hold, built from the public solver API rather than from the corridor code.
///
/// Deliberately a second implementation: the zero is solved in calm air, the muzzle angle
/// is applied, the segmented wind is set, one solve is sampled, and the hold is read by
/// interpolation. If `solve_robust_hold` ever stopped doing exactly that — re-zeroing per
/// scenario, say, or dropping the segments — this would diverge.
fn reference_hold(
    load: &CorridorLoad,
    segments: &[ballistics_engine::wind::WindSegment],
    range_m: f64,
    max_range_m: f64,
) -> (f64, f64) {
    let atmosphere = AtmosphericConditions {
        temperature: load.temperature_c,
        pressure: load.pressure_hpa,
        humidity: load.humidity_pct,
        altitude: load.altitude_m,
    };
    let inputs = BallisticInputs {
        bc_value: load.bc,
        bc_type: load.drag_model,
        bullet_mass: load.mass_kg,
        muzzle_velocity: load.muzzle_velocity_mps,
        bullet_diameter: load.diameter_m,
        bullet_length: load.bullet_length_m,
        sight_height: load.sight_height_m,
        use_rk4: true,
        ..Default::default()
    };
    let zero_angle = calculate_zero_angle_with_conditions(
        inputs.clone(),
        load.zero_distance_m,
        load.sight_height_m,
        WindConditions::default(),
        atmosphere.clone(),
    )
    .expect("zero");

    let mut inputs = inputs;
    inputs.muzzle_angle = zero_angle;
    inputs.enable_trajectory_sampling = true;
    inputs.sample_interval = 0.9144;

    let mut solver = TrajectorySolver::new(inputs, WindConditions::default(), atmosphere);
    solver.set_max_range(max_range_m);
    solver.set_time_step(0.001);
    solver.set_wind_segments(segments.to_vec());
    let result = solver.solve().expect("solve");
    let samples = result.sampled_points.expect("sampled points");

    let index = samples
        .partition_point(|s| s.distance_m < range_m)
        .min(samples.len() - 1);
    let (lo, hi) = (&samples[index - 1], &samples[index]);
    let t = (range_m - lo.distance_m) / (hi.distance_m - lo.distance_m);
    let lerp = |a: f64, b: f64| a + (b - a) * t;
    (
        lerp(lo.drop_m, hi.drop_m) / range_m * 1000.0,
        lerp(lo.wind_drift_m, hi.wind_drift_m) / range_m * 1000.0,
    )
}

/// **Acceptance criterion 1**: every scenario row equals a direct trajectory solve with the
/// same segmented wind.
#[test]
fn every_scenario_row_matches_a_direct_segmented_wind_solve() {
    let ranges_m = vec![182.88, 365.76, 548.64];
    let scenarios = vec![
        scenario("low", &["4:90:1000"]),
        scenario("high", &["14:90:1000"]),
        scenario("switch", &["10:90:400", "8:270:1000"]),
    ];
    let request = RobustHoldRequest {
        scenarios: WindScenarioSetV1 {
            version: 1,
            scenarios: scenarios.clone(),
            nominal: None,
        },
        ranges_m: ranges_m.clone(),
        target: None,
        load: load(),
    };
    let report = solve_robust_hold(&request).unwrap();
    let max_range_m = ranges_m.iter().fold(0.0_f64, |a, &b| a.max(b)) * 1.02;

    for row in &report.rows {
        for hold in &row.scenarios {
            let source = scenarios
                .iter()
                .find(|s| s.name == hold.name)
                .expect("scenario by name");
            let (elevation_mil, windage_mil) =
                reference_hold(&load(), &source.segments, row.range_m, max_range_m);
            assert_eq!(
                hold.elevation_mil, elevation_mil,
                "'{}' elevation at {} m must equal a direct solve",
                hold.name, row.range_m
            );
            assert_eq!(
                hold.windage_mil, windage_mil,
                "'{}' windage at {} m must equal a direct solve",
                hold.name, row.range_m
            );
        }
    }

    // Not a vacuous check: the scenarios really do produce different windage.
    let windages: Vec<f64> = report.rows[2]
        .scenarios
        .iter()
        .map(|h| h.windage_mil)
        .collect();
    assert!(
        windages
            .iter()
            .any(|w| (w - windages[0]).abs() > 0.3),
        "{windages:?}"
    );
}

fn write_scenarios(dir: &std::path::Path, body: &str) -> std::path::PathBuf {
    let path = dir.join("scenarios.json");
    std::fs::write(&path, body).unwrap();
    path
}

fn tempfile_dir() -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "hold-corridor-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(BIN).args(args).output().expect("hold-corridor");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

const THREE_SCENARIOS: &str = r#"{"version":1,"nominal":"low","scenarios":[
  {"name":"low","segments":["4:90:1000"]},
  {"name":"high","segments":["14:90:1000"]},
  {"name":"switch","segments":["10:90:400","8:270:1000"]}]}"#;

/// The CLI renders both shapes, and `-o json` really is the versioned report schema.
#[test]
fn cli_emits_the_table_and_the_versioned_json() {
    let dir = tempfile_dir();
    let path = write_scenarios(&dir, THREE_SCENARIOS);
    let p = path.to_str().unwrap();

    let (table, stderr, ok) = run(&[
        "hold-corridor", "--scenarios", p, "--ranges", "200,400,600", "--target", "rect:12x18",
        "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--zero-distance", "100",
    ]);
    assert!(ok, "{stderr}");
    assert!(table.contains("Robust Hold Corridor"), "{table}");
    assert!(table.contains("Minimax hold"), "{table}");
    assert!(table.contains("Holding the nominal"), "{table}");
    assert!(table.contains("Fits target"), "{table}");
    assert!(
        table.contains("NOT a probability interval"),
        "the non-goal must be stated in the output: {table}"
    );

    let (json, stderr, ok) = run(&[
        "hold-corridor", "--scenarios", p, "--ranges", "200,400,600", "--target", "circle:10",
        "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--zero-distance", "100",
        "-o", "json",
    ]);
    assert!(ok, "{stderr}");
    let report: RobustHoldReportV1 = serde_json::from_str(&json).unwrap();
    assert_eq!(report.version, 1);
    assert_eq!(report.metric, "circular");
    assert_eq!(report.scenario_names, vec!["high", "low", "switch"]);
    assert_eq!(report.nominal.as_deref(), Some("low"));
    assert_eq!(report.rows.len(), 3);
    assert!(matches!(report.target, Some(TargetSpec::Circle { .. })));
    for row in &report.rows {
        assert_eq!(row.scenarios.len(), 3);
        assert!(row.fits_target.is_some());
        // The corridor really contains every scenario, through the CLI too.
        for hold in &row.scenarios {
            assert!(hold.elevation_mil >= row.elevation_min_mil);
            assert!(hold.elevation_mil <= row.elevation_max_mil);
            assert!(hold.windage_mil >= row.windage_min_mil);
            assert!(hold.windage_mil <= row.windage_max_mil);
        }
        assert!(row.worst_case_miss_mil <= row.nominal_worst_case_miss_mil.unwrap() + 1e-12);
    }

    // No target: no fit verdict, and the metric falls back to per-axis.
    let (json, _, ok) = run(&[
        "hold-corridor", "--scenarios", p, "--ranges", "400",
        "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--zero-distance", "100",
        "-o", "json",
    ]);
    assert!(ok);
    let report: RobustHoldReportV1 = serde_json::from_str(&json).unwrap();
    assert_eq!(report.metric, "rectangular");
    assert!(report.rows[0].fits_target.is_none());
    assert!(report.target.is_none());
}

/// **Acceptance criterion 7, through the file**: reordering the scenarios in the document
/// changes nothing at all — not one byte of output.
#[test]
fn reordering_the_scenario_file_changes_no_output_byte() {
    let dir = tempfile_dir();
    let forward = write_scenarios(&dir, THREE_SCENARIOS);
    let reversed_dir = tempfile_dir();
    let reversed = write_scenarios(
        &reversed_dir,
        r#"{"version":1,"nominal":"low","scenarios":[
  {"name":"switch","segments":["10:90:400","8:270:1000"]},
  {"name":"high","segments":["14:90:1000"]},
  {"name":"low","segments":["4:90:1000"]}]}"#,
    );

    let output_for = |path: &std::path::Path| {
        let (stdout, stderr, ok) = run(&[
            "hold-corridor", "--scenarios", path.to_str().unwrap(), "--ranges", "200,400,600",
            "--target", "rect:12x18", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308",
            "--zero-distance", "100", "-o", "json",
        ]);
        assert!(ok, "{stderr}");
        stdout
    };
    assert_eq!(output_for(&forward), output_for(&reversed));
}

/// **Acceptance criterion 6, through the CLI**: caps, an unknown version, malformed
/// segments and bad targets are all structured errors — and they arrive before any solving.
#[test]
fn cli_rejects_caps_and_malformed_input() {
    let dir = tempfile_dir();
    let base = |extra: Vec<&str>, path: &str| {
        let mut args: Vec<&str> = vec!["hold-corridor", "--scenarios", path];
        args.extend(extra);
        args.extend_from_slice(&[
            "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--zero-distance", "100",
        ]);
        run(&args)
    };

    // Nine scenarios: rejected on count.
    let nine: Vec<String> = (0..9)
        .map(|i| format!(r#"{{"name":"s{i}","segments":["5:90:1000"]}}"#))
        .collect();
    let path = write_scenarios(
        &dir,
        &format!(r#"{{"version":1,"scenarios":[{}]}}"#, nine.join(",")),
    );
    let (_, stderr, ok) = base(vec!["--ranges", "400"], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("exceeds the limit of 8"), "{stderr}");

    // 65 ranges.
    let path = write_scenarios(&tempfile_dir(), THREE_SCENARIOS);
    let ranges: Vec<String> = (1..=65).map(|i| (i * 10).to_string()).collect();
    let joined = ranges.join(",");
    let (_, stderr, ok) = base(vec!["--ranges", &joined], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("exceeds the limit of 64"), "{stderr}");

    // A future version, even with otherwise-valid content.
    let path = write_scenarios(
        &tempfile_dir(),
        r#"{"version":2,"scenarios":[{"name":"a","segments":["5:90:1000"]}]}"#,
    );
    let (_, stderr, ok) = base(vec!["--ranges", "400"], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("unsupported wind scenario set version 2"), "{stderr}");

    // A malformed segment token names the scenario and the segment index.
    let path = write_scenarios(
        &tempfile_dir(),
        r#"{"version":1,"scenarios":[{"name":"bad","segments":["5:90"]}]}"#,
    );
    let (_, stderr, ok) = base(vec!["--ranges", "400"], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("scenario 'bad' segment #0"), "{stderr}");

    // A nominal that is not in the set lists what is.
    let path = write_scenarios(
        &tempfile_dir(),
        r#"{"version":1,"nominal":"nope","scenarios":[{"name":"a","segments":["5:90:1000"]}]}"#,
    );
    let (_, stderr, ok) = base(vec!["--ranges", "400"], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("is not in the set"), "{stderr}");

    // Bad targets.
    let path = write_scenarios(&tempfile_dir(), THREE_SCENARIOS);
    for spec in ["rect:12", "blob:3", "circle:-2", "rect:0x5"] {
        let (_, stderr, ok) = base(
            vec!["--ranges", "400", "--target", spec],
            path.to_str().unwrap(),
        );
        assert!(!ok, "target '{spec}' was accepted");
        assert!(stderr.contains("invalid target"), "{spec}: {stderr}");
    }

    // Duplicate ranges.
    let (_, stderr, ok) = base(vec!["--ranges", "400,600,400"], path.to_str().unwrap());
    assert!(!ok);
    assert!(stderr.contains("more than once"), "{stderr}");

    // A missing file names the path.
    let (_, stderr, ok) = base(vec!["--ranges", "400"], "/nonexistent/scenarios.json");
    assert!(!ok);
    assert!(stderr.contains("could not read scenario file"), "{stderr}");

    // CSV and PDF have no form.
    for format in ["csv", "pdf"] {
        let (_, stderr, ok) = base(
            vec!["--ranges", "400", "-o", format],
            path.to_str().unwrap(),
        );
        assert!(!ok, "{format} was accepted");
        assert!(stderr.contains("no "), "{format}: {stderr}");
    }
}

/// A cap violation is rejected on a document whose scenarios would each take a real solve,
/// and it returns fast — the ticket's "before unbounded work begins".
#[test]
fn cap_rejection_does_not_solve_anything() {
    let nine: Vec<String> = (0..9)
        .map(|i| format!(r#"{{"name":"s{i}","segments":["5:90:2000"]}}"#))
        .collect();
    let path = write_scenarios(
        &tempfile_dir(),
        &format!(r#"{{"version":1,"scenarios":[{}]}}"#, nine.join(",")),
    );
    // Sixty-four ranges out to 1500 yd would be nine long solves if any solving happened.
    let ranges: Vec<String> = (1..=64).map(|i| (i * 20).to_string()).collect();
    let joined = ranges.join(",");
    let start = std::time::Instant::now();
    let (_, stderr, ok) = run(&[
        "hold-corridor", "--scenarios", path.to_str().unwrap(), "--ranges", &joined,
        "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--zero-distance", "100",
    ]);
    let elapsed = start.elapsed();
    assert!(!ok);
    assert!(stderr.contains("exceeds the limit of 8"), "{stderr}");
    // Generous: one scenario's solve alone is far slower than this on any machine that can
    // run the rest of the suite.
    assert!(
        elapsed.as_secs() < 5,
        "the cap must be enforced before solving; took {elapsed:?}"
    );
}

/// Metric units are honoured end to end: segment speeds in m/s, ranges in meters, target
/// dimensions in centimeters.
#[test]
fn metric_units_are_applied_to_the_whole_request() {
    let path = write_scenarios(
        &tempfile_dir(),
        r#"{"version":1,"scenarios":[
  {"name":"calm","segments":["2:90:1000"]},
  {"name":"gusty","segments":["6:90:1000"]}]}"#,
    );
    let (json, stderr, ok) = run(&[
        "--units", "metric", "hold-corridor", "--scenarios", path.to_str().unwrap(),
        "--ranges", "200,400,600", "--target", "circle:30",
        "-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82", "--zero-distance", "100",
        "-o", "json",
    ]);
    assert!(ok, "{stderr}");
    let report: RobustHoldReportV1 = serde_json::from_str(&json).unwrap();
    // Ranges are already meters in metric mode, so they pass through unscaled.
    assert_eq!(report.ranges_m, vec![200.0, 400.0, 600.0]);
    assert_eq!(report.target, Some(TargetSpec::Circle { diameter_m: 0.30 }));
    // A 4 m/s spread at 600 m is a real corridor.
    assert!(report.rows[2].windage_max_mil - report.rows[2].windage_min_mil > 0.2);
}
