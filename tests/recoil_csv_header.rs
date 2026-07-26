//! MBA-1418: the `recoil` CSV header, pinned.
//!
//! The browser terminal emitted `velocity_m/s` in metric where native emits `velocity_mps`.
//! Machine output that differs by surface breaks anyone parsing both, and the divergence
//! survived because the imperial spelling (`fps`) happened to agree on both sides — so any test
//! that only exercised imperial would have passed throughout.
//!
//! This pins the native header, which is the side the WASM layer was corrected to match. The
//! browser surface is compile-verified only (its test module is browser-gated), so this file is
//! the machine-checkable half of the parity.

use std::path::PathBuf;
use std::process::Command;

fn cli() -> PathBuf {
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

fn header(units: &str) -> String {
    let output = Command::new(cli())
        .args([
            "--units", units, "recoil", "-b", "150", "-c", "45", "-v", "2700", "-f", "8", "-o",
            "csv",
        ])
        .output()
        .expect("recoil");
    assert!(
        output.status.success(),
        "recoil failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .next()
        .expect("a header line")
        .to_string()
}

#[test]
fn the_metric_csv_header_uses_ascii_unit_spellings() {
    let header = header("metric");
    assert_eq!(
        header,
        "bullet_weight_g,charge_weight_g,velocity_mps,firearm_weight_kg,gas_velocity_mps,\
         recoil_velocity_mps,recoil_energy_j,impulse_ns"
    );
}

#[test]
fn the_imperial_csv_header_uses_ascii_unit_spellings() {
    let header = header("imperial");
    assert_eq!(
        header,
        "bullet_weight_gr,charge_weight_gr,velocity_fps,firearm_weight_lb,gas_velocity_fps,\
         recoil_velocity_fps,recoil_energy_ftlb,impulse_lbs"
    );
}

/// The specific shape of the bug: a `/` in a machine header. It is not a comma, so it never broke
/// a parse outright — it just silently disagreed with the other surface.
#[test]
fn no_csv_header_field_contains_a_slash() {
    for units in ["metric", "imperial"] {
        let header = header(units);
        assert!(
            !header.contains('/'),
            "{units} header contains a '/': {header}"
        );
    }
}
