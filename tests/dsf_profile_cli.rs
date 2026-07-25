//! End-to-end test of the `dsf` CLI verb (MBA-1357): profile persistence
//! (`dsf_points`), `trajectory --saved-profile` auto-apply (drop bends, velocity/TOF
//! stay byte-identical), `profile show` rendering, `--clear-dsf`, the supersonic
//! rejection gate, and old-profile-without-the-key forward compatibility.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation
//! points `HOME` at a private temp dir (same pattern as `tests/a7p_import_cli.rs` /
//! `tests/compare_cli.rs`).

use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn cli() -> Command {
    Command::new(BIN)
}

fn tempfile_dir() -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "dsf-cli-test-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn profile_path(home: &std::path::Path, name: &str) -> std::path::PathBuf {
    home.join(".ballistics").join("profiles").join(format!("{name}.json"))
}

fn profile_json(home: &std::path::Path, name: &str) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(profile_path(home, name)).unwrap()).unwrap()
}

/// A muzzle velocity chosen just above Mach 1.2 (standard-atmosphere Mach 1 is ~1116
/// fps) so the shot crosses into the DSF-eligible Mach <= 1.2 band within the first
/// couple hundred yards, without needing a long, ground-impact-prone flight to get
/// there (a real .308 168gr load's BC and a modest zero distance keep the physics
/// unremarkable otherwise).
fn save_test_profile(home: &std::path::Path, name: &str) -> Output {
    cli()
        .env("HOME", home)
        .args([
            "profile",
            "save",
            name,
            "-v",
            "1300",
            "-b",
            "0.4",
            "-m",
            "168",
            "-d",
            "0.308",
            "--drag-model",
            "g1",
            "--zero-distance",
            "300",
            "--sight-height",
            "2.0",
        ])
        .output()
        .expect("spawn profile save")
}

fn full_trajectory_json(home: &std::path::Path, name: &str, max_range: f64, output: &str) -> Output {
    cli()
        .env("HOME", home)
        .args(["trajectory", "--saved-profile", name])
        .arg("--max-range")
        .arg(max_range.to_string())
        .args(["-o", output, "--full"])
        .output()
        .expect("spawn trajectory")
}

/// `trajectory --saved-profile NAME --sample-trajectory --full -o csv`: the PDF dope
/// card's data source and the Table's "Sampled Trajectory" section both read
/// `TrajectoryResult.sampled_points`, the SAME field this CSV path renders — a separate
/// array from `.points` (used by `full_trajectory_json` above), and the one Critical #2
/// found `apply_dsf` was silently skipping.
fn full_trajectory_csv_sampled(home: &std::path::Path, name: &str, max_range: f64) -> Output {
    cli()
        .env("HOME", home)
        .args(["trajectory", "--saved-profile", name])
        .arg("--max-range")
        .arg(max_range.to_string())
        .args(["--sample-trajectory", "-o", "csv", "--full"])
        .output()
        .expect("spawn trajectory sampled csv")
}

fn parse_json(out: &Output) -> serde_json::Value {
    serde_json::from_slice(&out.stdout)
        .unwrap_or_else(|e| panic!("invalid JSON ({e}): {}", String::from_utf8_lossy(&out.stdout)))
}

/// Find the first point in a full-trajectory JSON's `trajectory` array whose velocity
/// is comfortably subsonic (well under the Mach 1.2 gate — standard-atmosphere Mach 1
/// is ~1116 fps, so 1000 fps is safely under both the gate and the "silent accept"
/// 0.9-1.2 band, giving predictable test behavior).
fn pick_subsonic_point(traj: &serde_json::Value) -> (f64, f64) {
    let points = traj["trajectory"].as_array().expect("trajectory array (need --full)");
    for p in points {
        let v = p["velocity"].as_f64().unwrap();
        if v < 1000.0 {
            return (p["z"].as_f64().unwrap(), p["y"].as_f64().unwrap());
        }
    }
    panic!(
        "no subsonic (< 1000 fps) point found in {} points: {traj}",
        points.len()
    );
}

#[test]
fn dsf_verb_derives_saves_and_trajectory_auto_applies() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-roundtrip");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    // Baseline (pre-DSF) full trajectory.
    let baseline_out = full_trajectory_json(&home, "dsf-roundtrip", 1500.0, "json");
    assert!(baseline_out.status.success(), "{}", String::from_utf8_lossy(&baseline_out.stderr));
    let baseline = parse_json(&baseline_out);
    let (range_yd, y_yd) = pick_subsonic_point(&baseline);

    // LOS height (world frame, yd) = default bore height (60 in, no --bore-height flag
    // on `profile save`) + the --sight-height given above (2.0 in) — see
    // solve_profile_for_dsf's doc comment for why the bore-height default applies.
    let los_yd = (60.0 + 2.0) * 0.0254 / 0.9144;
    let predicted_drop_in = (los_yd - y_yd) * 36.0;
    assert!(
        predicted_drop_in > 0.0,
        "expected a real drop by the chosen subsonic point: {predicted_drop_in} in at {range_yd} yd"
    );

    // Observed drop deliberately 30% larger than predicted -> dsf ~= 1.3 (well inside
    // the engine's (0.5, 2.0) validity band).
    let observed_drop_in = predicted_drop_in * 1.3;
    let observed_token = format!("{observed_drop_in:.4}in");

    let dsf_out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", "dsf-roundtrip"])
        .arg("--range")
        .arg(format!("{range_yd:.2}"))
        .arg("--observed-drop")
        .arg(&observed_token)
        .output()
        .expect("spawn dsf");
    assert!(dsf_out.status.success(), "{}", String::from_utf8_lossy(&dsf_out.stderr));
    let dsf_stdout = String::from_utf8_lossy(&dsf_out.stdout);
    assert!(dsf_stdout.contains("DSF table for 'dsf-roundtrip'"), "{dsf_stdout}");
    // MBA-1357 Task 2 review, Important #1: the add/supersede announcement belongs on
    // stdout (the plan: "say so on stdout"), not stderr.
    assert!(dsf_stdout.contains("Added DSF point"), "{dsf_stdout}");

    // Saved profile now carries dsf_points with the expected ratio.
    let saved = profile_json(&home, "dsf-roundtrip");
    let points = saved["dsf_points"].as_array().expect("dsf_points present");
    assert_eq!(points.len(), 1);
    let dsf_ratio = points[0]["dsf"].as_f64().unwrap();
    assert!((dsf_ratio - 1.3).abs() < 0.05, "dsf ratio should be ~1.3: {dsf_ratio}");

    // `trajectory --saved-profile` now auto-applies: drop bends, velocity/TOF/max_range
    // stay byte-identical (the drop-only invariant, end-to-end through the real binary).
    let corrected_out = full_trajectory_json(&home, "dsf-roundtrip", 1500.0, "json");
    assert!(corrected_out.status.success());
    let corrected = parse_json(&corrected_out);
    assert_eq!(
        baseline["impact_velocity"], corrected["impact_velocity"],
        "impact velocity must stay identical"
    );
    assert_eq!(
        baseline["time_of_flight"], corrected["time_of_flight"],
        "time of flight must stay identical"
    );
    assert_eq!(
        baseline["max_range"], corrected["max_range"],
        "max_range must stay identical"
    );

    let base_points = baseline["trajectory"].as_array().unwrap();
    let corr_points = corrected["trajectory"].as_array().unwrap();
    assert_eq!(base_points.len(), corr_points.len());
    let mut any_drop_changed = false;
    for (bp, cp) in base_points.iter().zip(corr_points.iter()) {
        let bv = bp["velocity"].as_f64().unwrap();
        let cv = cp["velocity"].as_f64().unwrap();
        assert!((bv - cv).abs() < 1e-9, "per-point velocity must stay identical: {bv} vs {cv}");
        let bz = bp["z"].as_f64().unwrap();
        let cz = cp["z"].as_f64().unwrap();
        assert!((bz - cz).abs() < 1e-9, "per-point downrange must stay identical: {bz} vs {cz}");
        let by = bp["y"].as_f64().unwrap();
        let cy = cp["y"].as_f64().unwrap();
        if (by - cy).abs() > 1e-6 {
            any_drop_changed = true;
        }
    }
    assert!(any_drop_changed, "at least one point's drop must change after DSF auto-apply");

    // JSON gets no DSF text anywhere (purity rule) — check the raw stdout, not just
    // the parsed value, so an accidental extra top-level field would also be caught.
    let corrected_raw = String::from_utf8_lossy(&corrected_out.stdout);
    assert!(!corrected_raw.to_lowercase().contains("dsf"), "{corrected_raw}");

    // CSV: same purity rule.
    let csv_out = full_trajectory_json(&home, "dsf-roundtrip", 1500.0, "csv");
    assert!(csv_out.status.success());
    let csv_stdout = String::from_utf8_lossy(&csv_out.stdout);
    assert!(!csv_stdout.to_lowercase().contains("dsf"), "{csv_stdout}");

    // Table output: the note IS present (table-output-only).
    let table_out = cli()
        .env("HOME", &home)
        .args(["trajectory", "--saved-profile", "dsf-roundtrip", "--max-range", "1500"])
        .output()
        .expect("spawn trajectory table");
    assert!(table_out.status.success());
    let table_stdout = String::from_utf8_lossy(&table_out.stdout);
    assert!(table_stdout.contains("DSF table active (1 points"), "{table_stdout}");

    // `profile show` renders the table.
    let show_out = cli()
        .env("HOME", &home)
        .args(["profile", "show", "dsf-roundtrip"])
        .output()
        .expect("spawn profile show");
    assert!(show_out.status.success());
    let show_stdout = String::from_utf8_lossy(&show_out.stdout);
    assert!(show_stdout.contains("DSF table (1 points"), "{show_stdout}");

    // `profile save --clear-dsf` removes it.
    let clear_out = cli()
        .env("HOME", &home)
        .args([
            "profile", "save", "dsf-roundtrip", "-v", "2600", "-b", "0.3", "-m", "168", "-d",
            "0.308", "--drag-model", "g1", "--zero-distance", "300", "--sight-height", "2.0",
            "--clear-dsf",
        ])
        .output()
        .expect("spawn profile save --clear-dsf");
    assert!(clear_out.status.success(), "{}", String::from_utf8_lossy(&clear_out.stderr));
    let cleared = profile_json(&home, "dsf-roundtrip");
    assert!(cleared.get("dsf_points").is_none(), "dsf_points key must be gone: {cleared}");

    // And the auto-apply note is gone too.
    let table_out2 = cli()
        .env("HOME", &home)
        .args(["trajectory", "--saved-profile", "dsf-roundtrip", "--max-range", "1500"])
        .output()
        .expect("spawn trajectory table post-clear");
    assert!(table_out2.status.success());
    assert!(!String::from_utf8_lossy(&table_out2.stdout).contains("DSF table active"));
}

/// MBA-1357 Task 2 review, Critical #2: `apply_dsf` only rescaled `TrajectoryResult.points`,
/// never `sampled_points` — the array `--sample-trajectory` (Table's "Sampled Trajectory"
/// section, CSV `--full`, and the PDF dope card, which *always* requires sampling) reads.
/// Confirmed live pre-fix: with an active DSF table, sampled CSV rows were byte-identical
/// to the no-DSF baseline. This proves the sampled/CSV path now bends too.
#[test]
fn dsf_auto_apply_bends_sampled_trajectory_csv_rows() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-sampled");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let baseline_out = full_trajectory_csv_sampled(&home, "dsf-sampled", 1500.0);
    assert!(baseline_out.status.success(), "{}", String::from_utf8_lossy(&baseline_out.stderr));
    let baseline_csv = String::from_utf8_lossy(&baseline_out.stdout).into_owned();
    let mut lines = baseline_csv.lines();
    let header = lines.next().expect("csv header");
    assert!(header.starts_with("distance_yd,drop_in,drift_in,velocity_fps"), "{header}");

    // Pick a subsonic (< 1000 fps), already-below-LOS (positive drop) sampled row to
    // derive a DSF point from, same approach as pick_subsonic_point for the JSON path.
    let (range_yd, predicted_drop_in) = lines
        .clone()
        .find_map(|line| {
            let cols: Vec<&str> = line.split(',').collect();
            let vel: f64 = cols.get(3)?.parse().ok()?;
            let drop: f64 = cols.get(1)?.parse().ok()?;
            if vel < 1000.0 && drop > 0.0 {
                Some((cols[0].parse::<f64>().ok()?, drop))
            } else {
                None
            }
        })
        .expect("no subsonic, below-LOS sampled row found");

    let observed_drop_in = predicted_drop_in * 1.3;
    let dsf_out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", "dsf-sampled"])
        .arg("--range")
        .arg(format!("{range_yd:.2}"))
        .arg("--observed-drop")
        .arg(format!("{observed_drop_in:.4}in"))
        .output()
        .expect("spawn dsf");
    assert!(dsf_out.status.success(), "{}", String::from_utf8_lossy(&dsf_out.stderr));

    let corrected_out = full_trajectory_csv_sampled(&home, "dsf-sampled", 1500.0);
    assert!(corrected_out.status.success());
    let corrected_csv = String::from_utf8_lossy(&corrected_out.stdout).into_owned();
    let mut corrected_lines = corrected_csv.lines();
    assert_eq!(corrected_lines.next(), Some(header), "csv header must be unchanged");

    let mut any_drop_changed = false;
    for (base_line, corr_line) in baseline_csv.lines().skip(1).zip(corrected_lines) {
        let base_cols: Vec<&str> = base_line.split(',').collect();
        let corr_cols: Vec<&str> = corr_line.split(',').collect();
        // distance, velocity, energy, time (indices 0, 3, 4, 5) are the drop-only
        // invariant's untouched fields; drift (2) is untouched too.
        assert_eq!(base_cols[0], corr_cols[0], "distance must be byte-identical");
        assert_eq!(base_cols[2], corr_cols[2], "drift must be byte-identical");
        assert_eq!(base_cols[3], corr_cols[3], "velocity must be byte-identical");
        assert_eq!(base_cols[4], corr_cols[4], "energy must be byte-identical");
        assert_eq!(base_cols[5], corr_cols[5], "time must be byte-identical");

        let base_drop: f64 = base_cols[1].parse().unwrap();
        let corr_drop: f64 = corr_cols[1].parse().unwrap();
        if (base_drop - corr_drop).abs() > 1e-6 {
            any_drop_changed = true;
        }
    }
    assert!(
        any_drop_changed,
        "at least one sampled row's drop must change after DSF auto-apply \
         (pre-fix this CSV was byte-identical to the no-DSF baseline): \n{corrected_csv}"
    );
}

#[test]
fn profile_save_without_clear_dsf_carries_the_table_forward() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-carry-forward");
    assert!(save.status.success());

    let baseline_out = full_trajectory_json(&home, "dsf-carry-forward", 1500.0, "json");
    let baseline = parse_json(&baseline_out);
    let (range_yd, y_yd) = pick_subsonic_point(&baseline);
    let los_yd = (60.0 + 2.0) * 0.0254 / 0.9144;
    let predicted_drop_in = (los_yd - y_yd) * 36.0;
    let observed_token = format!("{:.4}in", predicted_drop_in * 1.2);

    let dsf_out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", "dsf-carry-forward"])
        .arg("--range")
        .arg(format!("{range_yd:.2}"))
        .arg("--observed-drop")
        .arg(&observed_token)
        .output()
        .expect("spawn dsf");
    assert!(dsf_out.status.success(), "{}", String::from_utf8_lossy(&dsf_out.stderr));

    // Re-save the SAME profile (e.g. tweaking an unrelated field) without --clear-dsf.
    let resave = cli()
        .env("HOME", &home)
        .args([
            "profile", "save", "dsf-carry-forward", "-v", "2600", "-b", "0.3", "-m", "168", "-d",
            "0.308", "--drag-model", "g1", "--zero-distance", "300", "--sight-height", "2.5",
        ])
        .output()
        .expect("spawn profile save (re-save)");
    assert!(resave.status.success(), "{}", String::from_utf8_lossy(&resave.stderr));

    let after_resave = profile_json(&home, "dsf-carry-forward");
    let points = after_resave["dsf_points"]
        .as_array()
        .expect("dsf_points must survive an ordinary re-save");
    assert_eq!(points.len(), 1);
}

#[test]
fn old_profile_without_dsf_points_key_loads_clean() {
    let home = tempfile_dir();
    let profiles_dir = home.join(".ballistics").join("profiles");
    std::fs::create_dir_all(&profiles_dir).unwrap();
    // Hand-write a pre-MBA-1357 profile with no `dsf_points` key at all, simulating a
    // profile saved before this field existed.
    let old_profile = serde_json::json!({
        "name": "legacy",
        "velocity": 2700.0,
        "bc": 0.475,
        "mass": 168.0,
        "diameter": 0.308,
        "drag_model": "G1",
        "units": "imperial",
        "temperature": 59.0,
        "pressure": 29.92,
        "humidity": 50.0,
        "altitude": 0.0
    });
    std::fs::write(
        profiles_dir.join("legacy.json"),
        serde_json::to_string_pretty(&old_profile).unwrap(),
    )
    .unwrap();

    let show = cli()
        .env("HOME", &home)
        .args(["profile", "show", "legacy"])
        .output()
        .expect("spawn profile show");
    assert!(show.status.success(), "{}", String::from_utf8_lossy(&show.stderr));
    assert!(!String::from_utf8_lossy(&show.stdout).contains("DSF table"));

    let traj = cli()
        .env("HOME", &home)
        .args(["trajectory", "--saved-profile", "legacy", "--max-range", "500"])
        .output()
        .expect("spawn trajectory");
    assert!(traj.status.success(), "{}", String::from_utf8_lossy(&traj.stderr));
    assert!(!String::from_utf8_lossy(&traj.stdout).contains("DSF table active"));
}

#[test]
fn supersonic_observation_is_rejected_with_exact_error_and_not_saved() {
    let home = tempfile_dir();
    // A normal, full-velocity .308 168gr load (unlike save_test_profile's deliberately
    // near-transonic muzzle velocity, used elsewhere to reach the subsonic band
    // quickly) — at 100 yd this is still comfortably supersonic (Mach ~2.2).
    let save = cli()
        .env("HOME", &home)
        .args([
            "profile", "save", "dsf-supersonic", "-v", "2700", "-b", "0.475", "-m", "168", "-d",
            "0.308", "--drag-model", "g1", "--zero-distance", "300", "--sight-height", "2.0",
        ])
        .output()
        .expect("spawn profile save");
    assert!(save.status.success());

    let out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", "dsf-supersonic", "--range", "100", "--observed-drop", "0.5in"])
        .output()
        .expect("spawn dsf");
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(
            "error: observation is supersonic (Mach"
        ) && stderr.contains(
            "calibrate muzzle velocity first (true-velocity), then collect DSF points at Mach <= 1.2"
        ),
        "{stderr}"
    );

    // The rejected observation must not have been saved (profile still has no dsf_points).
    let saved = profile_json(&home, "dsf-supersonic");
    assert!(saved.get("dsf_points").is_none(), "{saved}");
}

/// MBA-1357 Task 2 review, Critical #1: the 90%-of-solved-range warning must require
/// BOTH (a) the observation being beyond the trajectory's Mach-1.0 crossing and (b)
/// beyond 90% of the trajectory's OWN solved max range — not merely "close to however
/// far we chose to solve". `save_test_profile`'s load (v=1300 fps, BC 0.4, .308/168gr,
/// 300 yd zero) solves (with no `--max-range` of its own — `dsf` has none — the solve
/// envelope defaults to the profile's own configured max range, 1000 yd, same as
/// `trajectory --saved-profile`) to a natural ground-impact stop around 403 yd, and
/// crosses Mach 1.0 around 165 yd. Before the fix, the solve envelope was sized to
/// `--range * 1.05`, so `range / solved_max` was always in [0.952, 1.0] — this warning
/// fired unconditionally, including at the 340 yd mid-range point below.
#[test]
fn beyond_90pct_warning_requires_both_gates_not_just_range_ratio() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-90pct-gate");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    // Mid-range: subsonic (well past the ~165 yd Mach-1 crossing) but well under 90% of
    // the profile's own solved max range (~403 yd -> 90% ~= 363 yd). Must NOT warn.
    let mid = cli()
        .env("HOME", &home)
        .args([
            "dsf", "--saved-profile", "dsf-90pct-gate", "--range", "340", "--observed-drop",
            "20.7in",
        ])
        .output()
        .expect("spawn dsf mid-range");
    assert!(mid.status.success(), "{}", String::from_utf8_lossy(&mid.stderr));
    let mid_stderr = String::from_utf8_lossy(&mid.stderr);
    assert!(
        !mid_stderr.contains("beyond 90%"),
        "mid-range subsonic observation on a long-range profile must NOT warn: {mid_stderr}"
    );

    // Beyond 90% of the same profile's solved max range (and still subsonic, past the
    // crossing). Must warn.
    let far = cli()
        .env("HOME", &home)
        .args([
            "dsf", "--saved-profile", "dsf-90pct-gate", "--range", "380", "--observed-drop",
            "50.6in",
        ])
        .output()
        .expect("spawn dsf near-envelope");
    assert!(far.status.success(), "{}", String::from_utf8_lossy(&far.stderr));
    let far_stderr = String::from_utf8_lossy(&far.stderr);
    assert!(
        far_stderr.contains("beyond 90%"),
        "observation beyond 90% of the solved range must warn: {far_stderr}"
    );
}

#[test]
fn come_ups_auto_applies_the_dsf_table_with_a_table_only_note() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-come-ups");
    assert!(save.status.success());

    let come_ups = |output: &str| -> Output {
        cli()
            .env("HOME", &home)
            .args([
                "come-ups", "--profile", "dsf-come-ups", "--zero-distance", "300", "--start",
                "100", "--end", "400", "--step", "100", "-o", output,
            ])
            .output()
            .expect("spawn come-ups")
    };

    let baseline = come_ups("json");
    assert!(baseline.status.success(), "{}", String::from_utf8_lossy(&baseline.stderr));
    let baseline_json = parse_json(&baseline);

    // dsf mil (angular) at 400 yd; using --observed-drop in mil directly avoids
    // needing to hand-derive the linear predicted drop like the trajectory tests above.
    let dsf_out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", "dsf-come-ups", "--range", "400", "--observed-drop", "5mil"])
        .output()
        .expect("spawn dsf");
    assert!(dsf_out.status.success(), "{}", String::from_utf8_lossy(&dsf_out.stderr));

    let corrected = come_ups("json");
    assert!(corrected.status.success());
    let corrected_json = parse_json(&corrected);

    let base_rows = baseline_json["data"].as_array().unwrap();
    let corr_rows = corrected_json["data"].as_array().unwrap();
    assert_eq!(base_rows.len(), corr_rows.len());
    let mut any_drop_changed = false;
    for (b, c) in base_rows.iter().zip(corr_rows.iter()) {
        let bv = b["velocity"].as_f64().unwrap();
        let cv = c["velocity"].as_f64().unwrap();
        assert!((bv - cv).abs() < 1e-6, "velocity must stay identical: {bv} vs {cv}");
        let bt = b["time"].as_f64().unwrap();
        let ct = c["time"].as_f64().unwrap();
        assert!((bt - ct).abs() < 1e-6, "time must stay identical: {bt} vs {ct}");
        let bd = b["drop"].as_f64().unwrap();
        let cd = c["drop"].as_f64().unwrap();
        if (bd - cd).abs() > 1e-6 {
            any_drop_changed = true;
        }
    }
    assert!(any_drop_changed, "at least one row's drop must change after DSF auto-apply");

    // JSON gets no DSF text (purity rule).
    let corrected_raw = String::from_utf8_lossy(&corrected.stdout);
    assert!(!corrected_raw.to_lowercase().contains("dsf"), "{corrected_raw}");

    // Table output: the note IS present.
    let table = come_ups("table");
    assert!(table.status.success());
    let table_stdout = String::from_utf8_lossy(&table.stdout);
    assert!(table_stdout.contains("DSF table active (1 points"), "{table_stdout}");

    // CSV: purity rule again.
    let csv = come_ups("csv");
    assert!(csv.status.success());
    assert!(!String::from_utf8_lossy(&csv.stdout).to_lowercase().contains("dsf"));
}

// -----------------------------------------------------------------------------
// MBA-1405 Task 2: the `dsf` verb's validity-window note (stdout, always — `dsf`
// has no `--output` flag, so it is always the "table path"). Uses the solve the
// verb already performs (`solve_profile_for_dsf`'s `TrajectoryResult`) — no extra
// trajectory solve.
// -----------------------------------------------------------------------------

/// `save_test_profile`'s load (v=1300 fps, BC 0.4, .308/168gr, 300 yd zero) goes
/// supersonic-to-subsonic well inside its own solved range, so
/// `TrajectoryResult::mach_0_9_distance_m` is `Some` and the verb must print the
/// "at or beyond" window note. The 307.6 yd figure is a deterministic pin (90% of
/// this exact profile's downward Mach-0.9 crossing) — a real physics regression, not
/// a hand-picked constant.
#[test]
fn dsf_window_note_reports_the_start_of_the_subsonic_window() {
    let home = tempfile_dir();
    let save = save_test_profile(&home, "dsf-window-some");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let out = cli()
        .env("HOME", &home)
        .args([
            "dsf", "--saved-profile", "dsf-window-some", "--range", "340", "--observed-drop",
            "20.7in",
        ])
        .output()
        .expect("spawn dsf");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("note: DSF window: at or beyond 307.6 yd (90% of the Mach 0.9 distance)"),
        "{stdout}"
    );
}

/// A profile launched already subsonic (900 fps; standard-atmosphere Mach 1 is ~1116
/// fps, so this starts at Mach ~0.81) never records a downward Mach-0.9 CROSSING —
/// `MachTransitionTracker::record_downward_crossings` only fires on a transition from
/// at-or-above the threshold to below it, and this load never starts above 0.9. So
/// `mach_0_9_distance_m` stays `None` and the verb must print the "no subsonic
/// window" note instead.
#[test]
fn dsf_window_note_reports_no_subsonic_window_for_an_always_subsonic_load() {
    let home = tempfile_dir();
    let save = cli()
        .env("HOME", &home)
        .args([
            "profile", "save", "dsf-window-none", "-v", "900", "-b", "0.4", "-m", "168", "-d",
            "0.308", "--drag-model", "g1", "--zero-distance", "100", "--sight-height", "2.0",
        ])
        .output()
        .expect("spawn profile save");
    assert!(save.status.success(), "{}", String::from_utf8_lossy(&save.stderr));

    let out = cli()
        .env("HOME", &home)
        .args([
            "dsf", "--saved-profile", "dsf-window-none", "--range", "150", "--observed-drop",
            "1.5mil",
        ])
        .output()
        .expect("spawn dsf");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("note: no subsonic window inside the solved range"),
        "{stdout}"
    );
    assert!(
        !stdout.contains("note: DSF window: at or beyond"),
        "must not print the Some-branch note: {stdout}"
    );
}
