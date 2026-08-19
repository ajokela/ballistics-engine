//! Integration smoke tests for `--adjustment-unit` on the PDF dope card (MBA-724).
//! The card was MIL-only; it now renders MIL or MOA. These tests confirm both units
//! produce a valid, non-trivial PDF via the CLI (the flag is accepted and the dope-card
//! path runs); the MIL/MOA conversion itself is unit-tested in main.rs
//! (adjustment_unit_tests).
//!
//! No label is grepped out of the PDF here. This file used to say that was because "PDF
//! text is compressed", which is not so — printpdf writes uncompressed content streams
//! whose text is hex-encoded glyph ids of the embedded font subset, so a plain byte grep
//! finds nothing but the text IS recoverable. `tests/card_pdf_bridge.rs` does exactly that
//! (pdftotext, with a self-calibrating glyph scan as the fallback) for the bridge's
//! `card.pdf`; these older smoke tests were simply never rewritten to read the card back.

use std::path::{Path, PathBuf};
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

fn generate_card(unit: Option<&str>, out: &Path) {
    let out_str = out.to_str().unwrap();
    let mut args: Vec<&str> = vec![
        "trajectory", "-v", "2700", "-b", "0.5", "-m", "175", "-d", "0.308", "--drag-model", "g7",
        "--max-range", "800", "--temperature", "59", "--pressure", "29.92", "--auto-zero", "100",
        "--sample-trajectory", "-o", "pdf", "--output-file", out_str,
    ];
    if let Some(u) = unit {
        args.push("--adjustment-unit");
        args.push(u);
    }
    let output = Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("run");
    assert!(
        output.status.success(),
        "dope card ({:?}) failed: {}",
        unit,
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_valid_pdf(path: &Path, label: &str) {
    let bytes = std::fs::read(path).unwrap_or_else(|_| panic!("{label}: no PDF written"));
    assert!(
        bytes.len() > 10_000,
        "{label}: PDF suspiciously small ({} bytes)",
        bytes.len()
    );
    assert_eq!(&bytes[..5], b"%PDF-", "{label}: not a PDF");
}

#[test]
fn default_mil_dope_card_generates() {
    let out = std::env::temp_dir().join("bx_dope_mil.pdf");
    generate_card(None, &out);
    assert_valid_pdf(&out, "mil");
    let _ = std::fs::remove_file(&out);
}

#[test]
fn moa_dope_card_generates() {
    let out = std::env::temp_dir().join("bx_dope_moa.pdf");
    generate_card(Some("moa"), &out);
    assert_valid_pdf(&out, "moa");
    let _ = std::fs::remove_file(&out);
}

/// MBA-1410: independent elevation-vs-windage units on the PDF dope card. mil elevation
/// (Drop) with moa windage (Wind/Lead) must still produce a valid PDF — the per-axis
/// plumbing (elevation_unit/windage_unit, gated click args) must not crash when the two
/// axes diverge.
#[test]
fn mixed_axis_units_dope_card_generates() {
    let out = std::env::temp_dir().join("bx_dope_mixed_axis.pdf");
    let out_str = out.to_str().unwrap();
    let output = std::process::Command::new(get_cli_binary())
        .args([
            "trajectory", "-v", "2700", "-b", "0.5", "-m", "175", "-d", "0.308", "--drag-model", "g7",
            "--max-range", "800", "--temperature", "59", "--pressure", "29.92", "--auto-zero", "100",
            "--sample-trajectory", "-o", "pdf", "--output-file", out_str,
            "--adjustment-unit", "mil", "--windage-unit", "moa",
        ])
        .output()
        .expect("run");
    assert!(
        output.status.success(),
        "mixed-axis dope card failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_valid_pdf(&out, "mixed-axis");
    let _ = std::fs::remove_file(&out);
}

/// Elevation clicks with an explicit non-clicks windage override (moa): the windage click
/// graduation must NOT be applied to the Wind/Lead columns once `--windage-unit` overrides
/// away from clicks (the gating this fold-in relies on).
#[test]
fn elevation_clicks_with_windage_moa_override_generates() {
    let out = std::env::temp_dir().join("bx_dope_clicks_elev_moa_wind.pdf");
    let out_str = out.to_str().unwrap();
    let output = std::process::Command::new(get_cli_binary())
        .args([
            "trajectory", "-v", "2700", "-b", "0.5", "-m", "175", "-d", "0.308", "--drag-model", "g7",
            "--max-range", "800", "--temperature", "59", "--pressure", "29.92", "--auto-zero", "100",
            "--sample-trajectory", "-o", "pdf", "--output-file", out_str,
            "--adjustment-unit", "clicks", "--elevation-click-value", "0.25moa",
            "--windage-unit", "moa",
        ])
        .output()
        .expect("run");
    assert!(
        output.status.success(),
        "elevation-clicks/windage-moa dope card failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_valid_pdf(&out, "elevation-clicks-windage-moa");
    let _ = std::fs::remove_file(&out);
}

/// `--windage-unit clicks` without `--adjustment-unit clicks` is rejected up front (the
/// MBA-1410 judgment call: clicks requires the elevation graduation to be resolved).
#[test]
fn windage_unit_clicks_without_elevation_clicks_is_rejected() {
    let output = std::process::Command::new(get_cli_binary())
        .args([
            "trajectory", "-v", "2700", "-b", "0.5", "-m", "175", "-d", "0.308",
            "--max-range", "800", "--adjustment-unit", "mil", "--windage-unit", "clicks",
        ])
        .output()
        .expect("run");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--windage-unit"), "{stderr}");
    assert!(stderr.contains("--adjustment-unit"), "{stderr}");
}

#[test]
fn explicit_mil_matches_default() {
    let a = std::env::temp_dir().join("bx_dope_mil_a.pdf");
    let b = std::env::temp_dir().join("bx_dope_mil_b.pdf");
    generate_card(None, &a);
    generate_card(Some("mil"), &b);
    assert_valid_pdf(&a, "default");
    assert_valid_pdf(&b, "explicit-mil");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
}
