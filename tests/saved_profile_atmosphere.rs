//! MBA-1417: `trajectory --saved-profile` loads a profile's atmosphere as a SET.
//!
//! `profile save` has always stored temperature, pressure, humidity and altitude, and since
//! 0.29.0 also `pressure_reference` and `density_altitude` — but the trajectory path loaded none
//! of them, so a profile carried its ballistics and silently dropped its conditions.
//!
//! The delicate part is not loading the values, it is the QUALIFIERS. A pressure MODE describes
//! one specific pressure VALUE and may never be inherited apart from it: MBA-1397 measured that
//! exact mistake at 77 fps of silent error at 300 yd, when a profile saved with
//! `--pressure-type qnh` applied its mode to a CLI-supplied absolute reading. Density altitude is
//! the same shape with a larger blast radius, because it supersedes altitude and pressure
//! outright rather than merely describing them.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation points `HOME`
//! at a private temp dir (same pattern as `tests/dsf_profile_cli.rs`).

use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn private_home() -> std::path::PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-profile-atmo-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Save a profile carrying a full non-standard atmosphere at 5000 ft.
fn save_profile(home: &std::path::Path, name: &str, extra: &[&str]) {
    let output = Command::new(BIN)
        .env("HOME", home)
        .args([
            "profile", "save", name, "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7",
            "--mass", "175", "--diameter", "0.308",
        ])
        .args(extra)
        .output()
        .expect("profile save");
    assert!(
        output.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_stdout(home: &std::path::Path, args: &[&str]) -> String {
    let output = Command::new(BIN)
        .env("HOME", home)
        .args(["trajectory", "--max-range", "500", "-o", "json"])
        .args(args)
        .output()
        .expect("trajectory");
    assert!(
        output.status.success(),
        "trajectory failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn impact_velocity(home: &std::path::Path, args: &[&str]) -> f64 {
    let doc: serde_json::Value =
        serde_json::from_str(&run_stdout(home, args)).expect("json output");
    doc["impact_velocity"].as_f64().expect("impact_velocity")
}

const ATMO: [&str; 8] = [
    "--temperature",
    "59",
    "--pressure",
    "24.90",
    "--altitude",
    "5000",
    "--humidity",
    "50",
];

/// The stored conditions must actually reach the solve — the whole point of the ticket.
#[test]
fn a_profiles_stored_atmosphere_reaches_the_trajectory() {
    let home = private_home();
    save_profile(&home, "atmo", &ATMO);

    let with_profile = impact_velocity(&home, &["--saved-profile", "atmo"]);
    let bare = impact_velocity(
        &home,
        &[
            "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7", "-m", "175", "-d",
            "0.308",
        ],
    );

    assert!(
        (with_profile - bare).abs() > 1.0,
        "the profile's 5000 ft atmosphere changed nothing ({with_profile} vs {bare} fps)"
    );
}

/// A stored pressure mode applies when — and only when — the profile's pressure value is the one
/// in use.
#[test]
fn a_stored_pressure_mode_applies_to_its_own_stored_pressure() {
    let home = private_home();
    let mut args = ATMO.to_vec();
    args.extend_from_slice(&["--pressure-type", "qnh"]);
    save_profile(&home, "qnh", &args);
    save_profile(&home, "abs", &ATMO);

    let qnh = impact_velocity(&home, &["--saved-profile", "qnh"]);
    let absolute = impact_velocity(&home, &["--saved-profile", "abs"]);

    assert!(
        qnh > absolute,
        "a stored qnh mode should reduce the stored pressure to a thinner station pressure and \
         retain more velocity; got qnh {qnh} vs absolute {absolute}"
    );
}

/// The regression that motivated the whole set-wise design. A CLI pressure means the stored mode
/// describes some OTHER reading, so it must not be applied to it. Measured at 77 fps when this
/// was wrong.
#[test]
fn a_stored_mode_is_not_applied_to_a_cli_supplied_pressure() {
    let home = private_home();
    let mut args = ATMO.to_vec();
    args.extend_from_slice(&["--pressure-type", "qnh"]);
    save_profile(&home, "qnh", &args);

    let inherited = impact_velocity(&home, &["--saved-profile", "qnh", "--pressure", "24.90"]);
    let explicit_absolute = impact_velocity(&home, &[
        "--saved-profile",
        "qnh",
        "--pressure",
        "24.90",
        "--pressure-type",
        "absolute",
    ]);

    assert_eq!(
        inherited, explicit_absolute,
        "a CLI --pressure must be treated as absolute unless the CLI also declares a mode; the \
         profile's stored qnh leaked onto a value it does not describe"
    );
}

/// An explicit CLI mode still wins over the stored one.
#[test]
fn an_explicit_pressure_type_overrides_the_stored_one() {
    let home = private_home();
    let mut args = ATMO.to_vec();
    args.extend_from_slice(&["--pressure-type", "qnh"]);
    save_profile(&home, "qnh", &args);

    let stored = impact_velocity(&home, &["--saved-profile", "qnh"]);
    let overridden = impact_velocity(&home, &[
        "--saved-profile",
        "qnh",
        "--pressure-type",
        "absolute",
    ]);

    assert_ne!(
        stored, overridden,
        "--pressure-type absolute must override the profile's stored qnh"
    );
}

/// Profiles that never set an atmosphere take the standard defaults, so wiring the fields must
/// leave those runs byte-identical rather than nudging every existing profile.
#[test]
fn a_profile_without_a_stored_atmosphere_is_byte_identical_to_no_profile() {
    let home = private_home();
    save_profile(&home, "plain", &[]);

    let via_profile = run_stdout(&home, &["--saved-profile", "plain"]);
    let bare = run_stdout(
        &home,
        &[
            "--velocity", "2700", "--bc", "0.243", "--drag-model", "g7", "-m", "175", "-d",
            "0.308",
        ],
    );

    assert_eq!(
        via_profile, bare,
        "a profile storing only the standard atmosphere must not change the output"
    );
}

/// Density altitude supersedes altitude and pressure outright, so a stored DA must not override
/// a pressure the user typed in this run.
#[test]
fn a_stored_density_altitude_yields_to_a_cli_supplied_pressure() {
    let home = private_home();
    save_profile(&home, "da", &["--density-altitude", "6000"]);

    let stored_da = impact_velocity(&home, &["--saved-profile", "da"]);
    let with_cli_pressure =
        impact_velocity(&home, &["--saved-profile", "da", "--pressure", "29.92"]);

    assert_ne!(
        stored_da, with_cli_pressure,
        "a stored density altitude must yield to a pressure supplied in this run, or it would \
         silently supersede an input the user just typed"
    );
}
