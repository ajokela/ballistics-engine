//! MBA-1339: the native CLI `--bore-height` flag was unified onto inches (imperial) / mm
//! (metric) — matching `--sight-height` and the WASM `--muzzle-height` flag (previously feet
//! / metres) — and gained a `--muzzle-height` alias. These tests run the real binary and pin
//! the units, the preserved default, and the alias. The >1000 m bore-height sanity warning is
//! used as a clean unit discriminator.
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn out(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().expect("spawn ballistics")
}

fn ok_stdout(o: &Output) -> String {
    assert!(
        o.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    String::from_utf8_lossy(&o.stdout).into_owned()
}

const IMPERIAL: &[&str] = &[
    "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308", "--max-range", "300",
];

/// Imperial `--bore-height` is inches: 10000 in = 254 m (below the >1000 m warning). Under the
/// old feet units it would be 3048 m and trip the warning; 50000 in = 1270 m still warns.
#[test]
fn imperial_bore_height_is_inches() {
    let small = out(&[IMPERIAL, &["--bore-height", "10000"]].concat());
    assert!(small.status.success());
    assert!(
        !String::from_utf8_lossy(&small.stderr).contains("above ground"),
        "10000 in = 254 m must not trip the >1000 m warning (units are inches, not feet)"
    );
    let big = out(&[IMPERIAL, &["--bore-height", "50000"]].concat());
    assert!(
        String::from_utf8_lossy(&big.stderr).contains("above ground"),
        "50000 in = 1270 m should trip the >1000 m bore-height warning"
    );
}

/// Metric `--bore-height` is millimetres: 1500 mm = 1.5 m (no warning). Under the old metres
/// units it would be 1500 m and warn.
#[test]
fn metric_bore_height_is_mm() {
    const METRIC: &[&str] = &[
        "--units", "metric", "trajectory", "-v", "823", "-b", "0.475", "-m", "10.9", "-d", "7.82",
        "--max-range", "300",
    ];
    let o = out(&[METRIC, &["--bore-height", "1500"]].concat());
    assert!(o.status.success());
    assert!(
        !String::from_utf8_lossy(&o.stderr).contains("above ground"),
        "1500 mm = 1.5 m must not warn (units are mm, not metres)"
    );
}

/// The default bore height is unchanged: a run with no flag equals an explicit 60-inch bore
/// (the old 5 ft default = 1.524 m).
#[test]
fn imperial_default_is_preserved() {
    let default = ok_stdout(&out(IMPERIAL));
    let sixty = ok_stdout(&out(&[IMPERIAL, &["--bore-height", "60"]].concat()));
    assert_eq!(
        default, sixty,
        "default bore height must equal a 60-inch bore (unchanged behaviour)"
    );
}

/// `--muzzle-height` is a native alias for `--bore-height`, on identical units.
#[test]
fn muzzle_height_is_alias_for_bore_height() {
    let bore = ok_stdout(&out(&[IMPERIAL, &["--bore-height", "40"]].concat()));
    let muzzle = ok_stdout(&out(&[IMPERIAL, &["--muzzle-height", "40"]].concat()));
    assert_eq!(bore, muzzle, "--muzzle-height must alias --bore-height");
}
