//! MBA-1361: the `reticle` command family end to end — generators, the hold point, the
//! saved-profile attachment, and the solve-json wire block.
//!
//! The library-level focal-plane / nearest-mark / off-reticle math is pinned in
//! `src/reticle.rs`'s own tests; this file pins the SURFACES: that the native binary emits
//! the shared formatter verbatim (the `drag_curve_cli` contract, for the same MBA-1418
//! reason), that `generate -o json` round-trips into `hold --reticle-json` and
//! `profile save --reticle-json`, and that a request without a `reticle` block leaves the
//! solve-json response byte-identical.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every profile invocation
//! points `HOME` at a private temp dir (same pattern as `tests/dsf_profile_cli.rs`).

use std::process::Command;

use ballistics_engine::reticle::{
    format_reticle_description, format_reticle_hold, hold_point_in_reticle, FocalPlane, MarkKind,
    ReticleDescription, ReticleFormat, ReticleMark,
};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn tempfile_dir(tag: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "reticle-cli-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(BIN).args(args).output().expect("reticle");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

fn sfp_fixture() -> ReticleDescription {
    ReticleDescription {
        name: "fixture".to_string(),
        focal_plane: FocalPlane::Second,
        reference_magnification: 10.0,
        marks: vec![
            ReticleMark::new(0.0, 0.0, MarkKind::Center),
            ReticleMark::new(2.0, 0.0, MarkKind::Hash),
            ReticleMark::new(4.0, 0.0, MarkKind::Hash),
            ReticleMark::new(2.0, 1.0, MarkKind::Dot),
        ],
    }
}

fn write_fixture(dir: &std::path::Path, reticle: &ReticleDescription) -> std::path::PathBuf {
    let path = dir.join("reticle.json");
    std::fs::write(&path, serde_json::to_string(reticle).unwrap()).unwrap();
    path
}

/// The binary's stdout IS the shared formatter's string, byte for byte, for both the hold
/// and the description renderings. Equality rather than fragment checks, because a
/// fragment check is exactly what let the `recoil` CSV header divergence survive.
#[test]
fn native_output_is_the_shared_formatter_verbatim() {
    let dir = tempfile_dir("formatter");
    let reticle = sfp_fixture();
    let path = write_fixture(&dir, &reticle);

    for (flag, format) in [("table", ReticleFormat::Table), ("json", ReticleFormat::Json)] {
        let (stdout, stderr, ok) = run(&[
            "reticle",
            "hold",
            "--reticle-json",
            path.to_str().unwrap(),
            "--drop-mil",
            "3.6",
            "--wind-mil",
            "1.4",
            "--mag",
            "5",
            "-o",
            flag,
        ]);
        assert!(ok, "hold -o {flag} failed: {stderr}");
        let expected = format_reticle_hold(
            &hold_point_in_reticle(3.6, 1.4, 5.0, &reticle).unwrap(),
            &reticle,
            5.0,
            format,
        );
        assert_eq!(stdout, expected, "hold -o {flag} must be the shared formatter");
    }

    for (flag, format) in [("table", ReticleFormat::Table), ("json", ReticleFormat::Json)] {
        let (stdout, stderr, ok) = run(&[
            "reticle",
            "generate",
            "mil-grid",
            "--spacing",
            "1",
            "--extent",
            "5",
            "-o",
            flag,
        ]);
        assert!(ok, "generate -o {flag} failed: {stderr}");
        let expected =
            format_reticle_description(&ReticleDescription::mil_grid(1.0, 5.0).unwrap(), format);
        assert_eq!(
            stdout, expected,
            "generate -o {flag} must be the shared formatter"
        );
    }
}

/// `generate -o json` emits exactly what `hold --reticle-json` consumes — the round trip
/// the two commands are advertised to form.
#[test]
fn generate_json_round_trips_into_hold() {
    let dir = tempfile_dir("roundtrip");
    let (stdout, stderr, ok) = run(&[
        "reticle",
        "generate",
        "tree",
        "--rows",
        "4",
        "--row-spacing",
        "1.0",
        "--spread-step",
        "0.5",
        "--name",
        "generated",
        "-o",
        "json",
    ]);
    assert!(ok, "generate failed: {stderr}");
    let parsed: ReticleDescription = serde_json::from_str(&stdout).unwrap();
    assert_eq!(parsed.name, "generated");

    let path = dir.join("tree.json");
    std::fs::write(&path, &stdout).unwrap();
    let (hold_out, stderr, ok) = run(&[
        "reticle",
        "hold",
        "--reticle-json",
        path.to_str().unwrap(),
        "--drop-mil",
        "3.0",
        "--wind-mil",
        "1.5",
        "--mag",
        "10",
    ]);
    assert!(ok, "hold failed: {stderr}");
    // Row 3 (down 3.0) carries dots at +/-0.5, 1.0, 1.5 — the hold lands exactly on one.
    assert!(hold_out.contains("distance from hold: 0.000 mil"), "{hold_out}");
}

/// SFP subtension rescaling is visible through the CLI, not just the library: the same
/// hold reads a different mark at a different magnification.
#[test]
fn sfp_rescaling_changes_the_named_mark_through_the_cli() {
    let dir = tempfile_dir("sfp");
    let path = write_fixture(&dir, &sfp_fixture());
    let hold_at = |mag: &str| {
        let (stdout, stderr, ok) = run(&[
            "reticle",
            "hold",
            "--reticle-json",
            path.to_str().unwrap(),
            "--drop-mil",
            "4.0",
            "--mag",
            mag,
            "-o",
            "json",
        ]);
        assert!(ok, "hold at {mag}x failed: {stderr}");
        serde_json::from_str::<serde_json::Value>(&stdout).unwrap()
    };

    // At the reference magnification, 4 mil is the 4 mil mark (index 2).
    let at_reference = hold_at("10");
    assert_eq!(at_reference["nearest_mark"]["index"], 2);
    assert_eq!(at_reference["mark_scale"], 1.0);

    // At half of it the etched 2 mil mark covers 4 mil, so the same hold reads mark 1.
    let at_half = hold_at("5");
    assert_eq!(at_half["nearest_mark"]["index"], 1);
    assert_eq!(at_half["mark_scale"], 2.0);
    assert_eq!(at_half["nearest_mark"]["true_down_mil"], 4.0);
    // The hold point itself never moves — only the marks are rescaled.
    assert_eq!(at_half["hold"]["down_mil"], 4.0);
    assert_eq!(at_reference["hold"]["down_mil"], 4.0);
}

/// An FFP reticle reads the same mark at every magnification.
#[test]
fn ffp_hold_is_magnification_invariant_through_the_cli() {
    let dir = tempfile_dir("ffp");
    let mut reticle = sfp_fixture();
    reticle.focal_plane = FocalPlane::First;
    let path = write_fixture(&dir, &reticle);

    let mut seen = Vec::new();
    for mag in ["3", "10", "25"] {
        let (stdout, stderr, ok) = run(&[
            "reticle",
            "hold",
            "--reticle-json",
            path.to_str().unwrap(),
            "--drop-mil",
            "3.6",
            "--mag",
            mag,
            "-o",
            "json",
        ]);
        assert!(ok, "hold at {mag}x failed: {stderr}");
        let mut value: serde_json::Value = serde_json::from_str(&stdout).unwrap();
        // Everything except the echoed magnification must be identical.
        value["magnification"] = serde_json::Value::Null;
        seen.push(value);
    }
    assert_eq!(seen[0], seen[1]);
    assert_eq!(seen[1], seen[2]);
}

/// `--range` solves the trajectory and reads the hold off it, and the answer is the drop
/// `come-ups` reports for the same load and zero (both go through one hold curve).
#[test]
fn range_solve_matches_the_come_ups_drop() {
    let dir = tempfile_dir("range");
    let path = write_fixture(&dir, &ReticleDescription::mil_grid(1.0, 12.0).unwrap());
    let (stdout, stderr, ok) = run(&[
        "reticle",
        "hold",
        "--reticle-json",
        path.to_str().unwrap(),
        "--range",
        "600",
        "--mag",
        "10",
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
        "-o",
        "json",
    ]);
    assert!(ok, "range solve failed: {stderr}");
    let hold: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let down = hold["hold"]["down_mil"].as_f64().unwrap();

    let (come_ups, stderr, ok) = run(&[
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
        "600",
        "--end",
        "600",
        "--step",
        "100",
        "-o",
        "json",
    ]);
    assert!(ok, "come-ups failed: {stderr}");
    let rows: serde_json::Value = serde_json::from_str(&come_ups).unwrap();
    let expected = rows["data"][0]["drop"]
        .as_f64()
        .unwrap_or_else(|| panic!("come-ups shape changed: {come_ups}"));
    assert!(
        (down - expected).abs() < 0.02,
        "reticle --range drop {down} mil should match the come-up {expected} mil"
    );
}

/// Contradictory and incomplete invocations are rejected rather than silently resolved.
#[test]
fn rejects_ambiguous_and_incomplete_invocations() {
    let dir = tempfile_dir("reject");
    let path = write_fixture(&dir, &sfp_fixture());
    let p = path.to_str().unwrap();

    // Both hold sources at once.
    let (_, stderr, ok) = run(&[
        "reticle", "hold", "--reticle-json", p, "--drop-mil", "1", "--range", "500",
    ]);
    assert!(!ok);
    assert!(stderr.contains("mutually exclusive"), "{stderr}");

    // Neither hold source.
    let (_, stderr, ok) = run(&["reticle", "hold", "--reticle-json", p]);
    assert!(!ok);
    assert!(stderr.contains("firing solution is required"), "{stderr}");

    // No reticle at all.
    let (_, stderr, ok) = run(&["reticle", "hold", "--drop-mil", "1"]);
    assert!(!ok);
    assert!(stderr.contains("reticle is required"), "{stderr}");

    // --wind-mil is meaningless with --range (the wind comes from the wind flags there).
    let (_, stderr, ok) = run(&[
        "reticle", "hold", "--reticle-json", p, "--range", "500", "--wind-mil", "1",
    ]);
    assert!(!ok);
    assert!(stderr.contains("--wind-mil applies to --drop-mil only"), "{stderr}");

    // --range without a load.
    let (_, stderr, ok) = run(&["reticle", "hold", "--reticle-json", p, "--range", "500"]);
    assert!(!ok);
    assert!(stderr.contains("--range requires"), "{stderr}");

    // Non-physical magnification (rejected on both planes).
    let (_, stderr, ok) = run(&[
        "reticle", "hold", "--reticle-json", p, "--drop-mil", "1", "--mag", "0",
    ]);
    assert!(!ok, "{stderr}");

    // CSV/PDF have no reticle form.
    for format in ["csv", "pdf"] {
        let (_, stderr, ok) = run(&[
            "reticle", "hold", "--reticle-json", p, "--drop-mil", "1", "-o", format,
        ]);
        assert!(!ok);
        assert!(stderr.contains("no "), "{stderr}");
    }

    // A malformed description file is rejected at the flag that named it.
    let bad = dir.join("bad.json");
    std::fs::write(&bad, "{\"name\":\"x\"}").unwrap();
    let (_, stderr, ok) = run(&[
        "reticle",
        "hold",
        "--reticle-json",
        bad.to_str().unwrap(),
        "--drop-mil",
        "1",
    ]);
    assert!(!ok);
    assert!(stderr.contains("invalid reticle JSON"), "{stderr}");
}

/// A reticle attaches to a saved profile, survives an unrelated re-save, and clears on
/// request — the same three-way treatment the DSF table gets.
#[test]
fn profile_attachment_round_trips_and_carries_forward() {
    let home = tempfile_dir("profile");
    let reticle = sfp_fixture();
    let path = write_fixture(&home, &reticle);

    let save = |extra: &[&str]| {
        let mut args: Vec<&str> = vec![
            "profile", "save", "scoped", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308",
        ];
        args.extend_from_slice(extra);
        let output = Command::new(BIN)
            .env("HOME", &home)
            .args(&args)
            .output()
            .expect("profile save");
        assert!(
            output.status.success(),
            "profile save failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    };
    let stored = || -> serde_json::Value {
        serde_json::from_str(
            &std::fs::read_to_string(
                home.join(".ballistics").join("profiles").join("scoped.json"),
            )
            .unwrap(),
        )
        .unwrap()
    };

    // A profile saved without the flag has no `reticle` key at all.
    save(&[]);
    assert!(stored().get("reticle").is_none(), "{}", stored());

    save(&["--reticle-json", path.to_str().unwrap()]);
    let attached: ReticleDescription =
        serde_json::from_value(stored()["reticle"].clone()).unwrap();
    assert_eq!(attached, reticle);

    // An unrelated re-save carries it forward.
    save(&["--sight-height", "2.2"]);
    assert_eq!(
        serde_json::from_value::<ReticleDescription>(stored()["reticle"].clone()).unwrap(),
        reticle
    );

    // `reticle hold --profile` finds it.
    let output = Command::new(BIN)
        .env("HOME", &home)
        .args(["reticle", "hold", "--profile", "scoped", "--drop-mil", "4", "--mag", "5", "-o", "json"])
        .output()
        .expect("hold --profile");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value =
        serde_json::from_str(&String::from_utf8_lossy(&output.stdout)).unwrap();
    assert_eq!(value["reticle"], "fixture");
    assert_eq!(value["nearest_mark"]["index"], 1);

    // --clear-reticle drops the key entirely.
    save(&["--clear-reticle"]);
    assert!(stored().get("reticle").is_none(), "{}", stored());

    // A profile with no reticle names the fix rather than solving with nothing.
    let output = Command::new(BIN)
        .env("HOME", &home)
        .args(["reticle", "hold", "--profile", "scoped", "--drop-mil", "4"])
        .output()
        .expect("hold --profile");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("no attached reticle"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    // The two attachment flags contradict each other.
    let output = Command::new(BIN)
        .env("HOME", &home)
        .args([
            "profile", "save", "scoped", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308",
            "--reticle-json", path.to_str().unwrap(), "--clear-reticle",
        ])
        .output()
        .expect("profile save");
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("mutually exclusive"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

const SOLVE_REQUEST: &str = r#"{"schema_version":1,
 "projectile":{"mass_kg":0.01134,"diameter_m":0.00782,"drag_model":"G7","ballistic_coefficient":0.243},
 "rifle":{"muzzle_velocity_mps":823.0,"sight_height_m":0.05,"muzzle_height_m":1.2},
 "shot":{"max_range_m":700.0},
 "atmosphere":{"altitude_m":0.0,"temperature_k":293.15,"pressure_pa":100000.0,"relative_humidity":0.4},
 "wind":{"speed_mps":4.47,"direction_from_rad":1.5707963267948966},
 "solver":{"method":"rk45","time_step_s":0.001},
 "effects":{"magnus":false,"coriolis":false,"enhanced_spin_drift":false},
 "sampling":{"interval_m":100.0}}"#;

fn solve_json(request: &str) -> (serde_json::Value, bool) {
    use std::io::Write;
    use std::process::Stdio;
    let mut child = Command::new(BIN)
        .arg("solve-json")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("solve-json");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(request.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    (
        serde_json::from_slice(&output.stdout).unwrap_or_else(|_| {
            panic!(
                "solve-json emitted non-JSON: {}",
                String::from_utf8_lossy(&output.stdout)
            )
        }),
        output.status.success(),
    )
}

/// The `reticle` block is versioned-OPTIONAL: a request without it produces a response
/// with no `reticle_hold` key at all, and adding one changes nothing else in the response.
#[test]
fn solve_json_reticle_block_is_purely_additive() {
    let (plain, ok) = solve_json(SOLVE_REQUEST);
    assert!(ok, "{plain}");
    assert!(
        plain.get("reticle_hold").is_none(),
        "absent request block must leave the key absent, not null: {plain}"
    );

    let with_reticle = SOLVE_REQUEST.replace(
        "\"sampling\":{\"interval_m\":100.0}",
        "\"sampling\":{\"interval_m\":100.0},\
         \"reticle\":{\"range_m\":600.0,\"magnification\":5.0,\
         \"description\":{\"name\":\"fixture\",\"focal_plane\":\"sfp\",\
         \"reference_magnification\":10.0,\"marks\":[\
         {\"down_mil\":0.0,\"right_mil\":0.0,\"kind\":\"center\"},\
         {\"down_mil\":2.0,\"right_mil\":0.0,\"kind\":\"hash\",\"label\":\"600\"},\
         {\"down_mil\":4.0,\"right_mil\":0.0,\"kind\":\"hash\"}]}}",
    );
    let (mut held, ok) = solve_json(&with_reticle);
    assert!(ok, "{held}");
    let hold = held["reticle_hold"].clone();
    assert!(hold.is_object(), "{held}");
    assert_eq!(hold["range_m"], 600.0);
    assert_eq!(hold["magnification"], 5.0);
    assert_eq!(hold["mark_scale"], 2.0);
    // The hold must agree with the response's own samples at that range.
    let samples = held["samples"].as_array().unwrap();
    let at_600 = samples
        .iter()
        .find(|s| (s["distance_m"].as_f64().unwrap() - 600.0).abs() < 1e-9)
        .expect("a 600 m sample");
    let expected_down = at_600["drop_m"].as_f64().unwrap() / 600.0 * 1000.0;
    assert!(
        (hold["down_mil"].as_f64().unwrap() - expected_down).abs() < 1e-9,
        "{hold} vs sample {at_600}"
    );

    // Everything else is unchanged.
    held["reticle_hold"] = serde_json::Value::Null;
    let mut plain_padded = plain.clone();
    plain_padded["reticle_hold"] = serde_json::Value::Null;
    assert_eq!(held, plain_padded, "the reticle block must not move the solve");
}

/// A hold requested beyond the sampled trajectory is a structured error, never an
/// extrapolation.
#[test]
fn solve_json_rejects_a_reticle_range_outside_the_trajectory() {
    let request = SOLVE_REQUEST.replace(
        "\"sampling\":{\"interval_m\":100.0}",
        "\"sampling\":{\"interval_m\":100.0},\
         \"reticle\":{\"range_m\":5000.0,\"magnification\":5.0,\
         \"description\":{\"name\":\"x\",\"focal_plane\":\"ffp\",\
         \"reference_magnification\":1.0,\"marks\":[\
         {\"down_mil\":2.0,\"right_mil\":0.0,\"kind\":\"hash\"}]}}",
    );
    let (response, ok) = solve_json(&request);
    assert!(!ok, "{response}");
    assert_eq!(response["status"], "error");
    assert_eq!(response["error"]["code"], "invalid_value");
    assert_eq!(response["error"]["path"], "$.reticle.range_m");

    // An unknown key inside the strict envelope is rejected too.
    let request = SOLVE_REQUEST.replace(
        "\"sampling\":{\"interval_m\":100.0}",
        "\"sampling\":{\"interval_m\":100.0},\
         \"reticle\":{\"range_m\":600.0,\"magnification\":5.0,\"mag\":3.0,\
         \"description\":{\"name\":\"x\",\"focal_plane\":\"ffp\",\
         \"reference_magnification\":1.0,\"marks\":[\
         {\"down_mil\":2.0,\"right_mil\":0.0,\"kind\":\"hash\"}]}}",
    );
    let (response, ok) = solve_json(&request);
    assert!(!ok, "{response}");
    assert_eq!(response["status"], "error");
}
