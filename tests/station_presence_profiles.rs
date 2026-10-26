//! 0.49.0: saved profiles record whether their temperature and pressure were typed.
//!
//! A profile saved before 0.49.0 cannot tell a typed 59 °F / 29.92 inHg from the default it
//! stored for an omitted flag, so those files keep the legacy reading. `profile save` now writes
//! `temperature_supplied` / `pressure_supplied` for a typed value (and nothing at all for an
//! omitted one, so such a file is byte-identical to before), and `.a7p` import marks both,
//! because the file carries real zero-condition readings and has no concept of omission.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation points
//! `HOME` at a private temp dir.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn private_home() -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-station-presence-profile-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn save(home: &Path, name: &str, extra: &[&str]) -> String {
    let output = Command::new(BIN)
        .env("HOME", home)
        .args([
            "profile",
            "save",
            name,
            "--velocity",
            "2700",
            "--bc",
            "0.243",
            "--drag-model",
            "g7",
            "--mass",
            "175",
            "--diameter",
            "0.308",
            "--altitude",
            "5000",
        ])
        .args(extra)
        .output()
        .expect("profile save");
    assert!(
        output.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(
        home.join(".ballistics/profiles")
            .join(format!("{name}.json")),
    )
    .expect("saved profile file")
}

/// Every save stamps `created` with the time; replace it so two saves compare.
fn normalize_created(json: &str) -> String {
    let start = json.find("\"created\": \"").expect("created field") + "\"created\": \"".len();
    let end = start + json[start..].find('"').unwrap();
    format!("{}TIMESTAMP{}", &json[..start], &json[end..])
}

/// What ab4eb0c's `profile save` wrote for this command, `created` normalized.
const AB4EB0C_PLAIN: &str = r#"{
  "name": "plain",
  "velocity": 2700.0,
  "bc": 0.243,
  "mass": 175.0,
  "diameter": 0.308,
  "drag_model": "G7",
  "twist_rate": null,
  "sight_height": null,
  "zero_distance": null,
  "units": "imperial",
  "temperature": 59.0,
  "pressure": 29.92,
  "humidity": 50.0,
  "altitude": 5000.0,
  "bullet_name": null,
  "created": "TIMESTAMP"
}"#;

#[test]
fn a_profile_saved_without_typed_air_is_byte_identical_to_before() {
    let home = private_home();
    let json = save(&home, "plain", &[]);
    assert_eq!(normalize_created(&json), AB4EB0C_PLAIN);
}

#[test]
fn a_typed_pressure_or_temperature_is_marked() {
    let home = private_home();
    let pressure_only = save(&home, "p", &["--pressure", "29.92"]);
    assert!(
        pressure_only.contains("\"pressure_supplied\": true"),
        "{pressure_only}"
    );
    assert!(
        !pressure_only.contains("temperature_supplied"),
        "{pressure_only}"
    );

    let temperature_only = save(&home, "t", &["--temperature", "59"]);
    assert!(
        temperature_only.contains("\"temperature_supplied\": true"),
        "{temperature_only}"
    );
    assert!(
        !temperature_only.contains("pressure_supplied"),
        "{temperature_only}"
    );

    let both = save(&home, "b", &["--temperature", "59", "--pressure", "29.92"]);
    assert!(both.contains("\"temperature_supplied\": true"), "{both}");
    assert!(both.contains("\"pressure_supplied\": true"), "{both}");
}

#[test]
fn a_legacy_file_without_markers_still_loads() {
    let home = private_home();
    let dir = home.join(".ballistics/profiles");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("legacy.json"),
        AB4EB0C_PLAIN.replace("plain", "legacy"),
    )
    .unwrap();
    let output = Command::new(BIN)
        .env("HOME", &home)
        .args(["profile", "show", "legacy"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
}

#[test]
fn the_markers_round_trip_and_are_omitted_when_false() {
    use ballistics_engine::profile::ProfileData;
    let legacy: ProfileData = serde_json::from_str(AB4EB0C_PLAIN).unwrap();
    assert!(!legacy.temperature_supplied);
    assert!(!legacy.pressure_supplied);
    let written = serde_json::to_string(&legacy).unwrap();
    assert!(!written.contains("_supplied"), "{written}");

    let marked = ProfileData {
        pressure_supplied: true,
        ..legacy
    };
    let written = serde_json::to_string(&marked).unwrap();
    assert!(written.contains("\"pressure_supplied\":true"), "{written}");
    assert!(!written.contains("temperature_supplied"), "{written}");
    let back: ProfileData = serde_json::from_str(&written).unwrap();
    assert!(back.pressure_supplied && !back.temperature_supplied);
}
