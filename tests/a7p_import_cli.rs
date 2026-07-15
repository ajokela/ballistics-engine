//! End-to-end test of `ballistics profile import` against a synthetic .a7p.
#![cfg(feature = "profile-import")]

use std::process::Command;

fn cli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ballistics"))
}

// Spec-derived encoder (same shape as the unit-test helpers).
fn enc_varint(mut v: u64, out: &mut Vec<u8>) {
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
fn enc_i32(number: u32, value: i64, out: &mut Vec<u8>) {
    enc_varint(u64::from(number) << 3, out);
    enc_varint(value as u64, out);
}
fn enc_str(number: u32, s: &str, out: &mut Vec<u8>) {
    enc_varint((u64::from(number) << 3) | 2, out);
    enc_varint(s.len() as u64, out);
    out.extend_from_slice(s.as_bytes());
}
fn enc_bytes(number: u32, payload: &[u8], out: &mut Vec<u8>) {
    enc_varint((u64::from(number) << 3) | 2, out);
    enc_varint(payload.len() as u64, out);
    out.extend_from_slice(payload);
}

fn synthetic_a7p() -> Vec<u8> {
    let mut p = Vec::new();
    enc_str(1, "CLI 6.5CM", &mut p);
    enc_str(3, "140 ELD-M", &mut p);
    enc_i32(9, 50, &mut p); // 50 mm sight height
    enc_i32(10, 800, &mut p); // 8.00 in twist
    enc_i32(11, 8230, &mut p); // 823.0 m/s
    enc_i32(15, 15, &mut p);
    enc_i32(16, 10132, &mut p); // 1013.2 hPa
    enc_i32(17, 50, &mut p);
    enc_i32(20, 264, &mut p); // 0.264 in
    enc_i32(21, 1400, &mut p); // 140 gr
    enc_i32(22, 1360, &mut p); // 1.36 in
    enc_i32(24, 1, &mut p); // G7
    let mut packed = Vec::new();
    enc_varint(10_000, &mut packed);
    enc_bytes(26, &packed, &mut p);
    let mut row = Vec::new();
    enc_i32(1, 3260, &mut row); // G7 BC 0.326
    enc_i32(2, 8230, &mut row);
    enc_bytes(27, &row, &mut p);
    let mut payload = Vec::new();
    enc_bytes(1, &p, &mut payload);
    ballistics_engine::profile_import::wrap_payload(&payload)
}

#[test]
fn import_dry_run_saves_nothing_and_prints_report() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .arg("--dry-run")
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Imported fields:"), "{stdout}");
    assert!(stdout.contains("823.0 m/s"), "{stdout}");
    assert!(stdout.contains("Dry run"), "{stdout}");
    assert!(
        !home.join(".ballistics").join("profiles").join("CLI 6.5CM.json").exists(),
        "dry run must not write"
    );
}

#[test]
fn import_then_show_roundtrip() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "imported-cm"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Saved profile"), "{stdout}");

    let show = cli()
        .env("HOME", &home)
        .args(["profile", "show", "imported-cm"])
        .output()
        .unwrap();
    assert!(show.status.success(), "{show:?}");
    let shown = String::from_utf8_lossy(&show.stdout);
    assert!(shown.contains("imported-cm"), "{shown}");
    assert!(shown.contains("G7"), "{shown}");
}

#[test]
fn corrupted_checksum_warns_by_default_and_fails_with_strict() {
    let home = tempfile_dir();
    let mut bytes = synthetic_a7p();
    bytes[0] = if bytes[0] == b'0' { b'1' } else { b'0' };
    let a7p_path = home.join("bad.a7p");
    std::fs::write(&a7p_path, &bytes).unwrap();

    let lenient = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .arg("--dry-run")
        .output()
        .unwrap();
    assert!(lenient.status.success(), "{lenient:?}");
    let stdout = String::from_utf8_lossy(&lenient.stdout);
    assert!(stdout.contains("checksum mismatch"), "{stdout}");

    let strict = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--dry-run", "--strict"])
        .output()
        .unwrap();
    assert!(!strict.status.success(), "strict must fail on checksum mismatch");
}

#[test]
fn garbage_file_is_a_clean_error() {
    let home = tempfile_dir();
    let path = home.join("junk.a7p");
    std::fs::write(&path, b"not a profile").unwrap();
    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("a7p"), "{stderr}");
}

#[test]
fn import_with_slash_name_is_sanitized_with_note() {
    let home = tempfile_dir();
    let a7p_path = home.join("test.a7p");
    std::fs::write(&a7p_path, synthetic_a7p()).unwrap();

    let output = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "match/338"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("note: profile name sanitized to 'match_338' for the profile store"),
        "{stdout}"
    );

    let show = cli()
        .env("HOME", &home)
        .args(["profile", "show", "match_338"])
        .output()
        .unwrap();
    assert!(show.status.success(), "{show:?}");
}

// ============================================================================
// MBA-1323 Phase 2: multi-row BC and CUSTOM drag-curve consumption end-to-end
// ============================================================================

/// A synthetic G1 profile with TWO widely-separated BC rows (0.500 @ 900 m/s "muzzle" and
/// 0.200 @ 300 m/s "slow") so a run that honors bc_segments diverges visibly from a run using
/// only the scalar (fastest-row) BC.
fn multi_row_g1_a7p() -> Vec<u8> {
    let mut p = Vec::new();
    enc_str(1, "MULTI BC", &mut p);
    enc_str(3, "TEST BULLET", &mut p);
    enc_i32(9, 50, &mut p); // 50 mm sight height
    enc_i32(10, 800, &mut p); // 8.00 in twist
    enc_i32(11, 9000, &mut p); // 900.0 m/s muzzle velocity
    enc_i32(15, 15, &mut p); // 15 C
    enc_i32(16, 10132, &mut p); // 1013.2 hPa
    enc_i32(17, 50, &mut p); // 50%
    enc_i32(20, 264, &mut p); // 0.264 in
    enc_i32(21, 1400, &mut p); // 140.0 gr
    enc_i32(22, 1360, &mut p); // 1.360 in
    enc_i32(24, 0, &mut p); // G1
    let mut row1 = Vec::new();
    enc_i32(1, 5000, &mut row1); // BC 0.5000
    enc_i32(2, 9000, &mut row1); // at 900.0 m/s
    enc_bytes(27, &row1, &mut p);
    let mut row2 = Vec::new();
    enc_i32(1, 2000, &mut row2); // BC 0.2000
    enc_i32(2, 3000, &mut row2); // at 300.0 m/s
    enc_bytes(27, &row2, &mut p);
    let mut payload = Vec::new();
    enc_bytes(1, &p, &mut payload);
    ballistics_engine::profile_import::wrap_payload(&payload)
}

/// A synthetic CUSTOM (bc_type=2) profile: three (Cd, Mach) rows decoding to
/// (0.23, 0.5), (0.45, 1.2), (0.28, 3.0) — deliberately out of Mach order in the file to
/// exercise the importer's sort-before-validate step.
fn custom_curve_a7p() -> Vec<u8> {
    let mut p = Vec::new();
    enc_str(1, "CUSTOM CURVE", &mut p);
    enc_str(3, "TEST BULLET", &mut p);
    enc_i32(9, 50, &mut p);
    enc_i32(10, 800, &mut p);
    enc_i32(11, 8000, &mut p); // 800.0 m/s muzzle velocity
    enc_i32(15, 15, &mut p);
    enc_i32(16, 10132, &mut p);
    enc_i32(17, 50, &mut p);
    enc_i32(20, 264, &mut p);
    enc_i32(21, 1400, &mut p);
    enc_i32(22, 1360, &mut p);
    enc_i32(24, 2, &mut p); // CUSTOM
    for (cd, mach) in [(2800i64, 30000i64), (2300, 5000), (4500, 12000)] {
        let mut row = Vec::new();
        enc_i32(1, cd, &mut row);
        enc_i32(2, mach, &mut row);
        enc_bytes(27, &row, &mut p);
    }
    let mut payload = Vec::new();
    enc_bytes(1, &p, &mut payload);
    ballistics_engine::profile_import::wrap_payload(&payload)
}

/// (impact_velocity, time_of_flight, max_height) parsed from a `trajectory -o json` stdout.
fn key_metrics(stdout: &[u8]) -> (f64, f64, f64) {
    let json: serde_json::Value = serde_json::from_slice(stdout).expect("valid trajectory JSON");
    (
        json["impact_velocity"].as_f64().expect("impact_velocity"),
        json["time_of_flight"].as_f64().expect("time_of_flight"),
        json["max_height"].as_f64().expect("max_height"),
    )
}

#[test]
fn multi_row_import_makes_bc_segments_live_in_trajectory() {
    let home = tempfile_dir();
    let a7p_path = home.join("multi.a7p");
    std::fs::write(&a7p_path, multi_row_g1_a7p()).unwrap();

    let import = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "multi-bc"])
        .output()
        .unwrap();
    assert!(import.status.success(), "{import:?}");

    let profile_path = home
        .join(".ballistics")
        .join("profiles")
        .join("multi-bc.json");
    let saved = std::fs::read_to_string(&profile_path).unwrap();
    let mut profile_json: serde_json::Value = serde_json::from_str(&saved).unwrap();
    assert!(
        profile_json.get("bc_segments").is_some(),
        "profile JSON must contain bc_segments: {saved}"
    );
    let segs = profile_json["bc_segments"].as_array().unwrap();
    assert_eq!(segs.len(), 2, "both coef rows must be mapped: {saved}");

    let with_segments = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "multi-bc",
            "--units",
            "metric",
            "--max-range",
            "700",
            "-o",
            "json",
        ])
        .output()
        .unwrap();
    assert!(with_segments.status.success(), "{with_segments:?}");

    // Single-BC equivalent: same profile, minus bc_segments, so the solver falls back to
    // the scalar `bc` (the fastest row, retained for back-compat) alone.
    profile_json.as_object_mut().unwrap().remove("bc_segments");
    std::fs::write(
        &profile_path,
        serde_json::to_string_pretty(&profile_json).unwrap(),
    )
    .unwrap();

    let single_bc = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "multi-bc",
            "--units",
            "metric",
            "--max-range",
            "700",
            "-o",
            "json",
        ])
        .output()
        .unwrap();
    assert!(single_bc.status.success(), "{single_bc:?}");

    let (v_with, tof_with, _) = key_metrics(&with_segments.stdout);
    let (v_without, tof_without, _) = key_metrics(&single_bc.stdout);
    assert!(
        (v_with - v_without).abs() > 1.0 || (tof_with - tof_without).abs() > 0.01,
        "bc_segments must change the trajectory: with(v={v_with}, tof={tof_with}) \
         without(v={v_without}, tof={tof_without})"
    );
}

#[test]
fn custom_curve_import_matches_native_drag_table() {
    let home = tempfile_dir();
    let a7p_path = home.join("custom.a7p");
    std::fs::write(&a7p_path, custom_curve_a7p()).unwrap();

    let import = cli()
        .env("HOME", &home)
        .args(["profile", "import"])
        .arg(&a7p_path)
        .args(["--name", "custom-curve"])
        .output()
        .unwrap();
    assert!(import.status.success(), "{import:?}");

    let profile_path = home
        .join(".ballistics")
        .join("profiles")
        .join("custom-curve.json");
    let saved = std::fs::read_to_string(&profile_path).unwrap();
    let profile_json: serde_json::Value = serde_json::from_str(&saved).unwrap();
    assert_eq!(profile_json["drag_model"], "CUSTOM");
    assert_eq!(profile_json["bc"].as_f64(), Some(0.0)); // documented inert sentinel
    let curve = profile_json["drag_curve"].as_array().unwrap();
    assert_eq!(curve.len(), 3);

    let via_profile = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "custom-curve",
            "--units",
            "metric",
            "--max-range",
            "500",
            "-o",
            "json",
        ])
        .output()
        .unwrap();
    assert!(via_profile.status.success(), "{via_profile:?}");

    // The SAME curve (Mach 0.5/1.2/3.0 -> Cd 0.23/0.45/0.28), supplied via the native
    // --drag-table CSV path instead of the profile's embedded drag_curve.
    let csv_path = home.join("curve.csv");
    std::fs::write(&csv_path, "mach,cd\n0.5,0.23\n1.2,0.45\n3.0,0.28\n").unwrap();

    let via_drag_table = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "custom-curve",
            "--units",
            "metric",
            "--max-range",
            "500",
            "-o",
            "json",
            "--drag-table",
        ])
        .arg(&csv_path)
        .output()
        .unwrap();
    assert!(via_drag_table.status.success(), "{via_drag_table:?}");

    let (v1, tof1, h1) = key_metrics(&via_profile.stdout);
    let (v2, tof2, h2) = key_metrics(&via_drag_table.stdout);
    assert!(
        (v1 - v2).abs() < 1e-6 && (tof1 - tof2).abs() < 1e-9 && (h1 - h2).abs() < 1e-9,
        "profile-sourced and --drag-table-sourced curves must produce the identical \
         trajectory: profile(v={v1}, tof={tof1}, h={h1}) drag-table(v={v2}, tof={tof2}, h={h2})"
    );
}

#[test]
fn phase1_profile_without_v2_fields_solves_identically_to_explicit_flags() {
    let home = tempfile_dir();
    let profiles_dir = home.join(".ballistics").join("profiles");
    std::fs::create_dir_all(&profiles_dir).unwrap();
    // A Phase-1-era profile JSON: no bc_segments/drag_curve keys at all. temperature/
    // pressure/humidity/altitude are set to the metric standard atmosphere on purpose: the
    // `trajectory` command's --saved-profile support only ever merges velocity/bc/mass/
    // diameter/drag_model from a saved JSON profile (a PRE-EXISTING gap, not touched by
    // MBA-1323 Phase 2 — environmental fields there only come from --profile/--location CSV
    // rows or explicit --temperature/-pressure/-humidity/-altitude flags), so both this
    // fixture and the `via_flags` run below fall back to the same standard-atmosphere
    // defaults for those four fields regardless of what is written here.
    let phase1_json = r#"{
        "name": "phase1-legacy",
        "velocity": 823.0,
        "bc": 0.326,
        "mass": 9.0718474,
        "diameter": 6.7056,
        "drag_model": "G7",
        "units": "metric",
        "temperature": 15.0,
        "pressure": 1013.25,
        "humidity": 50.0,
        "altitude": 0.0
    }"#;
    std::fs::write(profiles_dir.join("phase1-legacy.json"), phase1_json).unwrap();

    let via_profile = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--saved-profile",
            "phase1-legacy",
            "--units",
            "metric",
            "--max-range",
            "500",
            "-o",
            "json",
        ])
        .output()
        .unwrap();
    assert!(via_profile.status.success(), "{via_profile:?}");

    // Deliberately omit --temperature/-pressure/-humidity/-altitude (see comment above) so
    // both runs resolve the same standard-atmosphere defaults; only velocity/bc/mass/
    // diameter/drag_model — the fields --saved-profile actually merges — are passed
    // explicitly here, matching the profile's stored values.
    let via_flags = cli()
        .env("HOME", &home)
        .args([
            "trajectory",
            "--units",
            "metric",
            "--max-range",
            "500",
            "-o",
            "json",
            "--velocity",
            "823.0",
            "--bc",
            "0.326",
            "--mass",
            "9.0718474",
            "--diameter",
            "6.7056",
            "--drag-model",
            "g7",
        ])
        .output()
        .unwrap();
    assert!(via_flags.status.success(), "{via_flags:?}");

    let (v1, tof1, h1) = key_metrics(&via_profile.stdout);
    let (v2, tof2, h2) = key_metrics(&via_flags.stdout);
    assert!(
        (v1 - v2).abs() < 1e-6 && (tof1 - tof2).abs() < 1e-9 && (h1 - h2).abs() < 1e-9,
        "a Phase-1-era profile (no v2 fields) must solve identically to the same values \
         passed explicitly: profile(v={v1}, tof={tof1}, h={h1}) flags(v={v2}, tof={tof2}, h={h2})"
    );
}

/// Unique temp dir per test without external crates.
fn tempfile_dir() -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "a7p-import-test-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}
