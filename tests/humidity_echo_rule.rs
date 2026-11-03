//! One echo rule for a humidity in the sensor band.
//!
//! A stored humidity of 100-105 % (or -5-0 %) is a sensor reading of saturated (or dry) air,
//! and every flight flies it as exactly 100 % (or 0 %). Two kinds of surface show it:
//!
//! - one that reports the air FLOWN echoes the end it flies as: `trajectory`'s stderr
//!   `Effective values` line, whatever the humidity's source (a flag, a `--location` CSV, a
//!   saved profile);
//! - one that displays STORED data shows the stored value, and when that is outside 0-100 says
//!   what it flies as: `profile show`, and the bridge's `profile.validate` and
//!   `profile.normalize`, which echo the stored document unchanged and add `flown_humidity`.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};

struct PrivateDir(PathBuf);

impl Drop for PrivateDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[allow(dead_code)]
fn private_dir() -> PrivateDir {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-echo-rule-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    PrivateDir(dir)
}

/// A stored profile document at `humidity`, as `profile save` writes one (imperial).
#[allow(dead_code)]
fn stored_profile(humidity: f64) -> serde_json::Value {
    serde_json::json!({
        "name": "band", "velocity": 2700.0, "bc": 0.243, "mass": 168.0,
        "diameter": 0.308, "drag_model": "G7", "units": "imperial",
        "temperature": 59.0, "pressure": 29.92, "humidity": humidity, "altitude": 0.0
    })
}

#[cfg(feature = "cli")]
mod cli {
    use super::*;
    use std::path::Path;
    use std::process::{Command, Output};

    const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

    fn home_with(humidity: f64) -> PrivateDir {
        let home = private_dir();
        let dir = home.0.join(".ballistics").join("profiles");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("band.json"), stored_profile(humidity).to_string()).unwrap();
        home
    }

    fn run_at(home: &Path, args: &[&str]) -> Output {
        Command::new(BIN)
            .env("HOME", home)
            .args(args)
            .output()
            .expect("spawn ballistics")
    }

    /// `(stored, how profile show prints it, what it flies as when that differs)`.
    const STORED: [(f64, &str, Option<&str>); 6] = [
        (102.0, "102.0%", Some("100.0%")),
        (105.0, "105.0%", Some("100.0%")),
        (-5.0, "-5.0%", Some("0.0%")),
        (-0.0, "0.0%", None),
        (50.0, "50.0%", None),
        (100.0, "100.0%", None),
    ];

    #[test]
    fn profile_show_displays_the_stored_humidity_and_what_a_band_reading_flies_as() {
        for (stored, shown, flies) in STORED {
            let home = home_with(stored);
            let output = run_at(&home.0, &["profile", "show", "band"]);
            assert!(output.status.success(), "{stored}");
            let stdout = String::from_utf8_lossy(&output.stdout);
            let line = stdout
                .lines()
                .find(|line| line.contains("Humidity:"))
                .unwrap_or_else(|| panic!("no humidity line: {stdout}"));
            let value = line
                .split("Humidity:")
                .nth(1)
                .unwrap()
                .trim_end_matches('║')
                .trim();
            assert_eq!(value, shown, "stored {stored}: {line}");
            let flies_line = stdout.lines().find(|line| line.contains("flies as:"));
            match flies {
                Some(end) => {
                    let flies_line = flies_line
                        .unwrap_or_else(|| panic!("stored {stored}: no flies-as line: {stdout}"));
                    assert!(flies_line.contains(end), "stored {stored}: {flies_line}");
                }
                None => assert!(flies_line.is_none(), "stored {stored}: {stdout}"),
            }
        }
    }

    #[test]
    fn a_flight_of_a_stored_band_reading_echoes_the_end_it_flies_as() {
        for (stored, flown) in [
            (102.0, "100.0"),
            (-5.0, "0.0"),
            (-0.0, "0.0"),
            (65.0, "65.0"),
        ] {
            let home = home_with(stored);
            let location = home.0.join("loc.csv");
            std::fs::write(&location, "NAME,TEMPERATURE\nS1,95\n").unwrap();
            let output = run_at(
                &home.0,
                &[
                    "trajectory",
                    "--saved-profile",
                    "band",
                    "--max-range",
                    "300",
                    "--location",
                    location.to_str().unwrap(),
                    "--site",
                    "S1",
                    "-o",
                    "json",
                ],
            );
            assert!(output.status.success(), "{stored}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            let line = stderr
                .lines()
                .find(|line| line.contains("altitude=") && line.contains("temp="))
                .unwrap_or_else(|| panic!("no location line: {stderr}"));
            assert!(
                line.ends_with(&format!("humidity={flown}%")),
                "stored {stored} flies {flown} %: {line}"
            );
        }
    }
}

#[cfg(feature = "bridge")]
mod bridge {
    use super::*;
    use ballistics_engine::bridge::{bridge_call, BRIDGE_API_VERSION};
    use serde_json::{json, Value};

    fn call(command: &str, humidity: f64) -> Value {
        let raw = bridge_call(
            &json!({
                "api_version": BRIDGE_API_VERSION,
                "command": command,
                "request": stored_profile(humidity),
            })
            .to_string(),
        );
        let out: Value = serde_json::from_str(&raw).unwrap();
        assert_eq!(out["ok"], true, "{command} {humidity}: {out}");
        out["result"].clone()
    }

    #[test]
    fn the_profile_documents_keep_the_stored_humidity_and_say_what_a_band_reading_flies_as() {
        for (stored, flown) in [
            (102.0, Some(100.0)),
            (-5.0, Some(0.0)),
            (50.0, None),
            (100.0, None),
        ] {
            let validate = call("profile.validate", stored);
            let normalize = call("profile.normalize", stored);
            assert_eq!(validate["valid"], true, "{stored}: {validate}");
            assert_eq!(
                validate["normalized"]["humidity"], stored,
                "stored, unchanged"
            );
            assert_eq!(
                normalize["profile"]["humidity"], stored,
                "stored, unchanged"
            );
            for result in [&validate, &normalize] {
                match flown {
                    Some(end) => assert_eq!(result["flown_humidity"], end, "{stored}: {result}"),
                    None => assert!(
                        result.get("flown_humidity").is_none(),
                        "{stored}: a humidity in 0-100 needs no note: {result}"
                    ),
                }
            }
        }
        // A humidity beyond the band flies as nothing: no note, and validate says why.
        let validate = call("profile.validate", 150.0);
        assert_eq!(validate["valid"], false);
        assert!(validate.get("flown_humidity").is_none(), "{validate}");
        assert!(call("profile.normalize", 150.0)
            .get("flown_humidity")
            .is_none());
    }
}
