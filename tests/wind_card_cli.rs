// MBA-727: oblique wind-card angles. The default (no angle flags) must stay
// byte-identical to the legacy full-value card; oblique angles are real solves.
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");
const BASE: &[&str] = &[
    "wind-card", "-v", "2700", "-m", "168", "-d", "0.308", "-b", "0.5",
    "--zero-distance", "100", "--end", "600",
];

fn run(extra: &[&str]) -> Output {
    Command::new(BIN).args(BASE).args(extra).output().expect("spawn")
}

/// Largest drift magnitude parsed from a card's table body (MIL numbers).
///
/// Each row is `│ Range │ drift @ ws1 │ drift @ ws2 │ ... │`; skip the leading
/// empty segment before the first `│` and the Range column itself (up to
/// `--end`, e.g. 600) so it can't be mistaken for a drift value.
fn max_drift(stdout: &str) -> f64 {
    stdout
        .lines()
        .filter(|l| l.contains('│'))
        .flat_map(|l| l.split('│').skip(2))
        .filter_map(|c| c.trim().parse::<f64>().ok())
        .fold(0.0_f64, |a, b| a.max(b.abs()))
}

#[test]
fn default_card_is_byte_identical_to_wind_angle_90_numbers() {
    // Same numbers; only labels may differ between legacy and angle-aware forms.
    let legacy = run(&[]);
    let ninety = run(&["--wind-angle", "90"]);
    assert!(legacy.status.success() && ninety.status.success());
    let (l, n) = (
        String::from_utf8_lossy(&legacy.stdout).to_string(),
        String::from_utf8_lossy(&ninety.stdout).to_string(),
    );
    assert!(l.contains("full-value"), "legacy default must keep legacy labeling");
    let (dl, dn) = (max_drift(&l), max_drift(&n));
    assert!(dl > 0.0, "legacy card must contain drift numbers");
    assert!((dl - dn).abs() < 1e-9, "90° numbers must equal legacy: {dl} vs {dn}");
}

#[test]
fn oblique_45_is_near_sin45_of_full_value() {
    let full = max_drift(&String::from_utf8_lossy(&run(&[]).stdout));
    let out45 = run(&["--wind-angle", "45"]);
    assert!(out45.status.success());
    let d45 = max_drift(&String::from_utf8_lossy(&out45.stdout));
    let expected = full * 45f64.to_radians().sin();
    assert!(
        (d45 - expected).abs() <= 0.05 * expected,
        "45° drift {d45} should be within 5% of sin45×full {expected}"
    );
}

#[test]
fn head_and_tail_winds_have_near_zero_drift() {
    for a in ["0", "180"] {
        let out = run(&["--wind-angle", a]);
        assert!(out.status.success());
        let d = max_drift(&String::from_utf8_lossy(&out.stdout));
        assert!(d < 0.05, "angle {a} drift should be ~0, got {d}");
    }
}

#[test]
fn multi_angle_emits_one_card_per_angle() {
    let out = run(&["--wind-angles", "30,60,90"]);
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    for a in ["30", "60", "90"] {
        assert!(
            s.contains(&format!("{a}°")) || s.contains(&format!("angle {a}")),
            "output must label the {a}° card"
        );
    }
}

#[test]
fn oblique_json_is_valid_and_angle_tagged() {
    let out = run(&["--wind-angle", "45", "-o", "json"]);
    assert!(out.status.success());
    let v: serde_json::Value =
        serde_json::from_slice(&out.stdout).expect("wind-card JSON must parse");
    assert_eq!(v["wind_angle"].as_f64(), Some(45.0));
}

#[test]
fn conflicting_angle_flags_are_rejected() {
    let out = run(&["--wind-angle", "45", "--wind-angles", "30,60"]);
    assert!(!out.status.success());
}
