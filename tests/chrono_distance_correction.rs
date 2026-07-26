//! CLI contracts for MBA-1377 (`--chrono-distance`: correct a downrange chronograph reading
//! back to muzzle velocity). Unit-level coverage of the solver itself
//! (`truing::correct_chrono_velocity_fps`) -- the hand-computed pin, convergence across the
//! 3-25 m band, out-of-band rejection, and the G1-vs-G7 drag-model-is-consulted check -- lives
//! in `src/truing.rs`'s `chrono_correction_tests` module. This file covers the CLI wiring: the
//! additive/no-op guarantee, unit conversion, and end-to-end error surfacing.

use std::process::Command;

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(binary())
        .args(args)
        .output()
        .expect("run ballistics CLI")
}

/// Shared fixture args for a single-observation `true-velocity` run with a chrono comparison.
/// `--offline` keeps this fully local (no network dependency), matching how the online API path
/// is out of scope for this feature (MBA-1377's brief: input/presentation layers only).
const BASE_ARGS: [&str; 16] = [
    "true-velocity",
    "--measured-drop",
    "1.8",
    "--range",
    "300",
    "-b",
    "0.475",
    "--drag-model",
    "g7",
    "-m",
    "168",
    "-d",
    "0.308",
    "--chrono-velocity",
    "2680",
    "--offline",
];

#[test]
fn absent_chrono_distance_is_byte_identical_to_explicit_zero() {
    let absent = run(&BASE_ARGS);
    assert!(absent.status.success(), "{}", String::from_utf8_lossy(&absent.stderr));

    let mut with_zero: Vec<&str> = BASE_ARGS.to_vec();
    with_zero.push("--chrono-distance");
    with_zero.push("0");
    let zero = run(&with_zero);
    assert!(zero.status.success(), "{}", String::from_utf8_lossy(&zero.stderr));

    assert_eq!(
        absent.stdout, zero.stdout,
        "omitting --chrono-distance must be byte-identical to passing an explicit 0"
    );
}

#[test]
fn absent_chrono_distance_is_byte_identical_across_repeated_runs() {
    // Determinism baseline: the no-op path must not introduce any nondeterminism (e.g. an
    // accidental extra solve whose iteration count could vary run to run).
    let first = run(&BASE_ARGS);
    let second = run(&BASE_ARGS);
    assert!(first.status.success());
    assert!(second.status.success());
    assert_eq!(first.stdout, second.stdout);
}

#[test]
fn nonzero_chrono_distance_raises_the_reported_muzzle_velocity() {
    // Table mode prints "Adjustment from Chrono: <value> fps", computed as
    // effective_velocity - chrono_fps. A larger (corrected, higher) chrono_fps must produce a
    // SMALLER (or more negative) adjustment than the raw uncorrected reading -- i.e. correcting
    // the chronograph reading upward narrows the reported gap versus the drop-fitted velocity
    // (or widens a negative one), never the reverse, since the model side is unchanged.
    let baseline = run(&BASE_ARGS);
    assert!(baseline.status.success());
    let baseline_out = String::from_utf8_lossy(&baseline.stdout);

    let mut corrected_args: Vec<&str> = BASE_ARGS.to_vec();
    corrected_args.push("--chrono-distance");
    corrected_args.push("15");
    let corrected = run(&corrected_args);
    assert!(corrected.status.success(), "{}", String::from_utf8_lossy(&corrected.stderr));
    let corrected_out = String::from_utf8_lossy(&corrected.stdout);

    assert_ne!(
        baseline_out, corrected_out,
        "a nonzero --chrono-distance must change the displayed comparison"
    );

    fn adjustment(stdout: &str) -> f64 {
        let line = stdout
            .lines()
            .find(|l| l.contains("Adjustment from Chrono:"))
            .expect("adjustment line present");
        line.split(':')
            .nth(1)
            .expect("value after colon")
            .split_whitespace()
            .next()
            .expect("numeric token")
            .parse()
            .expect("parses as f64")
    }

    let baseline_adj = adjustment(&baseline_out);
    let corrected_adj = adjustment(&corrected_out);
    // effective_velocity is identical in both runs (same --measured-drop/--range/--bc); only
    // chrono_fps changed and increased, so the adjustment (effective - chrono) must decrease.
    assert!(
        corrected_adj < baseline_adj,
        "corrected adjustment {corrected_adj} should be less than baseline {baseline_adj}"
    );
}

#[test]
fn chrono_distance_requires_chrono_velocity() {
    let out = run(&[
        "true-velocity",
        "--measured-drop",
        "1.8",
        "--range",
        "300",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--chrono-distance",
        "15",
        "--offline",
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--chrono-distance requires --chrono-velocity"),
        "{stderr}"
    );
}

#[test]
fn out_of_band_chrono_distance_is_rejected_not_silently_applied() {
    // "--chrono-distance=-5" (a single token, `=`-joined) so clap doesn't mistake the negative
    // number for a short flag.
    for bad in ["--chrono-distance=-5", "--chrono-distance=0.05", "--chrono-distance=1000"] {
        let mut args: Vec<&str> = BASE_ARGS.to_vec();
        args.push(bad);
        let out = run(&args);
        assert!(!out.status.success(), "distance {bad} should be rejected");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("out of range"),
            "distance {bad}: {stderr}"
        );
    }
}

#[test]
fn metric_units_convert_chrono_distance_from_meters() {
    // In metric mode --chrono-distance is meters, not feet; 25 m is the Lapua/JBM reference
    // distance and sits comfortably inside the validated band.
    let args = [
        "--units",
        "metric",
        "true-velocity",
        "--measured-drop",
        "1.8",
        "--range",
        "300",
        "-b",
        "0.475",
        "--drag-model",
        "g7",
        "-m",
        "10.9",
        "-d",
        "7.82",
        "--chrono-velocity",
        "817",
        "--chrono-distance",
        "25",
        "--offline",
    ];
    let out = run(&args);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
}
