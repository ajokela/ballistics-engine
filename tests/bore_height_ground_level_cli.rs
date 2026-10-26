//! A bore at ground level (`--bore-height 0`) on the native CLI `trajectory` command.
//!
//! The ground plane sits at height 0 and the flight stops when it reaches it, so a bore at or
//! below 0 has no flight to record unless `--ignore-ground-impact` removes the ground. Two
//! defects made that combination fail:
//!
//! - the `--auto-zero` solve hardcoded a ground threshold of 0 and ignored
//!   `--ignore-ground-impact`, so `--bore-height 0 --auto-zero N --ignore-ground-impact` died
//!   with "No trajectory points calculated" even though the flight itself honored the flag;
//! - without the flag, the same input surfaced only that opaque solver error ("...generated"
//!   under `--use-euler`) instead of saying what was wrong and how to fix it.
//!
//! These tests run the real binary on the load from the report (.308, 175 gr, 2600 fps,
//! G1 0.505) and its metric twin.
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn run(args: &[&str]) -> Output {
    Command::new(BIN)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

const IMPERIAL: &[&str] = &[
    "trajectory",
    "-v",
    "2600",
    "-b",
    "0.505",
    "-m",
    "175",
    "-d",
    "0.308",
];

const METRIC: &[&str] = &[
    "--units",
    "metric",
    "trajectory",
    "-v",
    "792",
    "-b",
    "0.505",
    "-m",
    "11.34",
    "-d",
    "7.82",
];

/// RK45 (default), fixed-step RK4 and Euler.
const INTEGRATORS: &[&[&str]] = &[&[], &["--use-rk4-fixed"], &["--use-euler"]];

/// Run an `--auto-zero ... -o json` command and return the solved zero angle in degrees.
fn zero_angle_deg(args: &[&str]) -> f64 {
    let o = run(args);
    assert!(
        o.status.success(),
        "{args:?} failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let json: serde_json::Value =
        serde_json::from_slice(&o.stdout).expect("trajectory -o json emits JSON");
    let angle = json["zero_angle_degrees"]
        .as_f64()
        .expect("an --auto-zero run reports zero_angle_degrees");
    assert!(angle.is_finite(), "{args:?}: zero angle {angle}");
    angle
}

/// The reported command: `--bore-height 0 --auto-zero 1500 --ignore-ground-impact` must solve,
/// in both unit systems and under every integrator. A bore at ground level is the limit of a
/// tiny positive bore, so the solved zero angle must match a bore a thousandth of an inch (or a
/// hundredth of a millimetre) up, which solved before the fix.
#[test]
fn auto_zero_from_ground_level_with_ignore_ground_impact_succeeds() {
    for (base, zero, tiny_bore) in [(IMPERIAL, "1500", "0.001"), (METRIC, "1372", "0.01")] {
        for integrator in INTEGRATORS {
            let flags = |bore: &'static str| {
                let mut v = base.to_vec();
                v.extend_from_slice(&[
                    "--bore-height",
                    bore,
                    "--auto-zero",
                    zero,
                    "--ignore-ground-impact",
                    "--max-range",
                    zero,
                    "-o",
                    "json",
                ]);
                v.extend_from_slice(integrator);
                v
            };
            let ground = zero_angle_deg(&flags("0"));
            let tiny = zero_angle_deg(&flags(tiny_bore));
            assert!(
                (ground - tiny).abs() < 1e-6,
                "{:?}: ground-level zero {ground} deg vs tiny-bore zero {tiny} deg",
                flags("0")
            );
        }
    }
}

/// Without `--ignore-ground-impact`, a bore at or below the ground plane is rejected up front
/// with an error that names the flag and the way out, not the solver's opaque
/// "No trajectory points calculated/generated".
#[test]
fn bore_on_or_below_ground_without_ignore_is_rejected_with_actionable_error() {
    let extras: &[&[&str]] = &[
        &[],
        &["--auto-zero", "1500"],
        &["--use-rk4-fixed"],
        &["--use-euler"],
    ];
    for base in [IMPERIAL, METRIC] {
        // `=` form: clap would read a separate "-1" as a flag.
        for bore in ["--bore-height=0", "--bore-height=-1", "--muzzle-height=0"] {
            for extra in extras {
                let mut args = base.to_vec();
                args.push(bore);
                args.extend_from_slice(extra);
                let o = run(&args);
                let stderr = String::from_utf8_lossy(&o.stderr);
                assert!(!o.status.success(), "{args:?} should be rejected");
                assert!(
                    stderr.contains("--bore-height/--muzzle-height must be above 0")
                        && stderr.contains("--ignore-ground-impact"),
                    "{args:?}: unexpected error: {stderr}"
                );
                assert!(
                    !stderr.contains("No trajectory points"),
                    "{args:?}: still the opaque solver error: {stderr}"
                );
            }
        }
    }
}

/// The new check does not get in the way of the path that already worked: a ground-level bore
/// with `--ignore-ground-impact` and no `--auto-zero` still flies, from height 0.
#[test]
fn ground_level_with_ignore_ground_impact_still_flies() {
    let mut args = IMPERIAL.to_vec();
    args.extend_from_slice(&[
        "--bore-height",
        "0",
        "--ignore-ground-impact",
        "--full",
        "-o",
        "json",
    ]);
    let o = run(&args);
    assert!(
        o.status.success(),
        "{args:?} failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&o.stdout).expect("JSON output");
    let first_y = json["trajectory"][0]["y"]
        .as_f64()
        .expect("--full emits trajectory points with a y height");
    assert_eq!(first_y, 0.0, "the flight starts on the ground plane");
}
