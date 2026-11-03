//! The relative-humidity axis of `error-budget --sigma` and
//! `tolerance --domain` says its unit, refuses a value no fraction can be, and can be searched
//! from a saturated or a dry nominal.
//!
//! Both commands take a solve-json v1 request FILE, whose `relative_humidity` is a fraction
//! (0.5 = 50 %), and their axes are that file's own axes, in its own units. So the humidity
//! sigma and domain are fractions too, unlike every `--humidity` flag. A percent typed there
//! (`--sigma relative-humidity=5`, `--domain relative-humidity=20:80`) used to be accepted or
//! refused for an unrelated reason; it is now refused with the unit named. There is no percent
//! alias.
//!
//! A sensor reading of 100-105 % (or -5-0 %) flies as saturated (or dry) air, so its request
//! resolves to a humidity of exactly 1 (or 0), the edge of the axis. `tolerance` searches such
//! a nominal in the one direction there is, and reports the other as the physical boundary,
//! the way `error-budget` already takes a one-sided difference there. A calm wind (0 m/s) is
//! the same case for `wind-speed`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

struct PrivateDir(PathBuf);

impl Drop for PrivateDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn private_dir() -> PrivateDir {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-perturbation-unit-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    PrivateDir(dir)
}

/// A small request at `relative_humidity` (a fraction) and `wind_mps` of crosswind.
fn request(dir: &PrivateDir, relative_humidity: f64, wind_mps: f64) -> PathBuf {
    let path = dir.0.join(format!("r-{relative_humidity}-{wind_mps}.json"));
    let value = serde_json::json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                       "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
        "shot": {"max_range_m": 700.0},
        "atmosphere": {"temperature_k": 308.15, "relative_humidity": relative_humidity},
        "wind": {"speed_mps": wind_mps, "direction_from_rad": std::f64::consts::FRAC_PI_2},
        "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
    });
    std::fs::write(&path, value.to_string()).unwrap();
    path
}

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(BIN).args(args).output().expect("run");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

fn error_budget(request: &Path, sigma: &str) -> (String, String, bool) {
    run(&[
        "error-budget",
        "--request",
        request.to_str().unwrap(),
        "--ranges",
        "600",
        "--sigma",
        sigma,
        "--sigma",
        "wind-speed=0.5",
    ])
}

fn tolerance(
    request: &Path,
    axis: &str,
    domain: Option<&str>,
    json: bool,
) -> (String, String, bool) {
    let mut args = vec![
        "tolerance",
        "--request",
        request.to_str().unwrap(),
        "--range",
        "600",
        "--target",
        "rect:8x8",
        "--axis",
        axis,
    ];
    if let Some(domain) = domain {
        args.extend(["--domain", domain]);
    }
    if json {
        args.extend(["-o", "json"]);
    }
    run(&args)
}

#[test]
fn a_humidity_sigma_of_one_or_more_is_refused_as_a_fraction() {
    let dir = private_dir();
    let r = request(&dir, 0.5, 4.0);
    for sigma in [
        "relative-humidity=5",
        "relative-humidity=1",
        "relative-humidity=50",
    ] {
        let (_, stderr, ok) = error_budget(&r, sigma);
        assert!(!ok, "--sigma {sigma} must be refused");
        assert!(
            stderr.contains("fraction") && stderr.contains("0.05 = 5 percentage points"),
            "--sigma {sigma}: the refusal names the unit: {stderr}"
        );
    }
    let (_, stderr, ok) = error_budget(&r, "relative-humidity=0.05");
    assert!(ok, "a sigma of 0.05 (5 points) runs: {stderr}");
}

#[test]
fn a_humidity_domain_outside_zero_to_one_is_refused_as_a_fraction_before_anything_else() {
    let dir = private_dir();
    let r = request(&dir, 0.5, 4.0);
    for domain in [
        "relative-humidity=20:80",
        "relative-humidity=0:1.0000001",
        "relative-humidity=-0.1:0.9",
    ] {
        let (_, stderr, ok) = tolerance(&r, "relative-humidity", Some(domain), false);
        assert!(!ok, "--domain {domain} must be refused");
        assert!(
            stderr.contains("a fraction from 0 to 1 (0.5 = 50 %)"),
            "--domain {domain}: the refusal names the unit and the range: {stderr}"
        );
        assert!(
            !stderr.contains("strictly inside"),
            "--domain {domain}: refused for its unit, not for where the nominal sits: {stderr}"
        );
    }
    let (_, stderr, ok) = tolerance(&r, "wind-speed", Some("wind-speed=-5:20"), false);
    assert!(!ok, "a negative wind-speed domain must be refused");
    assert!(stderr.contains("m/s"), "names the unit: {stderr}");
}

#[test]
fn tolerance_names_the_humidity_unit_in_its_hint_and_its_table() {
    let dir = private_dir();
    let r = request(&dir, 0.5, 4.0);
    let (_, stderr, ok) = tolerance(&r, "relative-humidity", None, false);
    assert!(!ok);
    assert!(
        stderr.contains("--domain relative-humidity=LO:HI (as a fraction, 0.5 = 50 %)"),
        "{stderr}"
    );
    let (stdout, stderr, ok) = tolerance(
        &r,
        "relative-humidity",
        Some("relative-humidity=0.2:0.8"),
        false,
    );
    assert!(ok, "{stderr}");
    assert!(
        stdout.contains("axis relative-humidity (nominal 0.500000 fraction,"),
        "{stdout}"
    );
}

#[test]
fn tolerance_searches_a_saturated_or_dry_nominal_toward_the_side_there_is() {
    let dir = private_dir();
    // (request humidity, domain, the side at the physical boundary). 1.02 and -0.02 are sensor
    // readings that resolve to exactly 1 and 0.
    for (sent, domain, closed) in [
        (1.0, "relative-humidity=0.5:1", "far"),
        (1.02, "relative-humidity=0.5:1", "far"),
        (0.0, "relative-humidity=0:0.5", "near"),
        (-0.02, "relative-humidity=0:0.5", "near"),
    ] {
        let r = request(&dir, sent, 4.0);
        let (stdout, stderr, ok) = tolerance(&r, "relative-humidity", Some(domain), true);
        assert!(ok, "RH {sent} over {domain}: {stderr}");
        let report: serde_json::Value = serde_json::from_str(&stdout).unwrap();
        let axis = &report["axes"][0];
        let open = if closed == "far" { "near" } else { "far" };
        assert_eq!(
            axis[format!("{closed}_at_physical_boundary")],
            true,
            "{axis}"
        );
        assert_eq!(
            axis[format!("{open}_at_physical_boundary")],
            false,
            "{axis}"
        );
        assert!(axis[format!("{closed}_bound")].is_null(), "{axis}");
        assert_eq!(axis[format!("{closed}_has_no_effect")], false, "{axis}");

        let (table, _, ok) = tolerance(&r, "relative-humidity", Some(domain), false);
        assert!(ok);
        assert!(
            table.contains("at the axis's physical boundary"),
            "RH {sent}: the table says why that side was not searched: {table}"
        );
    }
    // A calm wind is the same case for wind-speed.
    let calm = request(&dir, 0.5, 0.0);
    let (stdout, stderr, ok) = tolerance(&calm, "wind-speed", Some("wind-speed=0:20"), true);
    assert!(ok, "a calm nominal: {stderr}");
    let report: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(report["axes"][0]["near_at_physical_boundary"], true);
    assert!(
        report["axes"][0]["far_bound"].is_number(),
        "the far side is searched and found: {}",
        report["axes"][0]
    );
}

#[test]
fn the_sigma_and_domain_help_give_the_humidity_unit() {
    for command in ["error-budget", "tolerance"] {
        let (stdout, _, ok) = run(&[command, "--help"]);
        assert!(ok);
        assert!(
            stdout.contains("relative-humidity") && stdout.contains("a fraction"),
            "{command} --help: {stdout}"
        );
    }
}
