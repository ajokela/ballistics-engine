//! Integration tests for MBA-1397: `--pressure-type <absolute|qnh>` (and the paired
//! `--zero-pressure-type` auto-zero override), the sea-level-corrected altimeter-setting
//! (QNH) pressure mode.
//!
//! Hard requirement: with the toggle absent or explicitly `absolute`, every surface must be
//! byte-identical to pre-MBA-1397 output. The `*_byte_identical_to_omitted` tests below
//! compare raw stdout bytes (not just parsed field subsets) to prove that.

use serde_json::Value;
use std::path::PathBuf;
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

fn run(args: &[&str]) -> std::process::Output {
    Command::new(get_cli_binary())
        .args(args)
        .output()
        .expect("failed to execute ballistics binary")
}

fn run_json(args: &[&str]) -> Value {
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "stdout is not JSON: {error}; stdout={}",
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

const BASE_TRAJECTORY_ARGS: &[&str] = &[
    "trajectory",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--units",
    "metric",
    "--max-range",
    "300",
    "--ignore-ground-impact",
    "-o",
    "json",
    "--altitude",
    "1500",
    "--pressure",
    "1030.0",
];

const BASE_ZERO_ARGS: &[&str] = &[
    "zero",
    "-v",
    "2700",
    "-b",
    "0.475",
    "-m",
    "168",
    "-d",
    "0.308",
    "--units",
    "metric",
    "--target-distance",
    "300",
    "-o",
    "json",
    "--altitude",
    "1500",
    "--pressure",
    "1030.0",
];

fn with_extra<'a>(base: &[&'a str], extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = base.to_vec();
    args.extend_from_slice(extra);
    args
}

// ---- Hard requirement: byte-identical when absent/absolute ---------------------------------

#[test]
fn trajectory_omitted_pressure_type_is_byte_identical_to_explicit_absolute() {
    let omitted = run(&with_extra(BASE_TRAJECTORY_ARGS, &["--full"]));
    let explicit = run(&with_extra(
        BASE_TRAJECTORY_ARGS,
        &["--full", "--pressure-type", "absolute"],
    ));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --pressure-type must be byte-identical to --pressure-type absolute"
    );
}

#[test]
fn zero_omitted_pressure_type_is_byte_identical_to_explicit_absolute() {
    let omitted = run(BASE_ZERO_ARGS);
    let explicit = run(&with_extra(BASE_ZERO_ARGS, &["--pressure-type", "absolute"]));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --pressure-type must be byte-identical to --pressure-type absolute on zero"
    );
}

#[test]
fn trajectory_auto_zero_omitted_zero_pressure_type_is_byte_identical_to_explicit_absolute() {
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "1500",
        "--pressure",
        "1013.25",
        "--auto-zero",
        "300",
        "--zero-pressure",
        "1030.0",
    ];
    let omitted = run(args);
    let explicit = run(&with_extra(args, &["--zero-pressure-type", "absolute"]));
    assert!(omitted.status.success());
    assert!(explicit.status.success());
    assert_eq!(
        omitted.stdout, explicit.stdout,
        "omitted --zero-pressure-type must be byte-identical to --zero-pressure-type absolute"
    );
}

// ---- QNH must actually change the physics (reduced, not raw, pressure reaches the solve) ---

/// Hand-computed pin (see also src/atmosphere.rs's
/// `qnh_reduction_matches_hand_computed_value_at_nonzero_altitude`): reducing a QNH of
/// 1030.0 hPa at 1500 m gives 859.575312392644... hPa station pressure -- strictly, and
/// substantially, lower than the raw 1030.0 hPa. Lower pressure means less dense air, less
/// drag, and a HIGHER retained impact velocity at a fixed range.
#[test]
fn trajectory_qnh_mode_yields_a_higher_impact_velocity_than_absolute_at_the_same_reading() {
    let absolute = run_json(BASE_TRAJECTORY_ARGS);
    let qnh = run_json(&with_extra(BASE_TRAJECTORY_ARGS, &["--pressure-type", "qnh"]));

    let v_absolute = absolute["impact_velocity"].as_f64().unwrap();
    let v_qnh = qnh["impact_velocity"].as_f64().unwrap();
    assert!(
        v_qnh > v_absolute + 10.0,
        "QNH-reduced (lower) pressure must retain velocity better than treating the same \
         reading as absolute station pressure: absolute={v_absolute} qnh={v_qnh}"
    );
}

#[test]
fn zero_qnh_mode_yields_a_smaller_zero_angle_than_absolute_at_the_same_reading() {
    let absolute = run_json(BASE_ZERO_ARGS);
    let qnh = run_json(&with_extra(BASE_ZERO_ARGS, &["--pressure-type", "qnh"]));

    let angle_absolute = absolute["zero_angle_degrees"].as_f64().unwrap();
    let angle_qnh = qnh["zero_angle_degrees"].as_f64().unwrap();
    // Less dense air (QNH correctly reduced) means less drop over the same zero distance,
    // so a SMALLER compensating zero angle is needed.
    assert!(
        angle_qnh < angle_absolute,
        "QNH-reduced pressure must need a smaller zero angle: absolute={angle_absolute} \
         qnh={angle_qnh}"
    );
    assert!((angle_absolute - angle_qnh).abs() > 1e-4);
}

/// `--zero-pressure-type` is independent of the shot-day `--pressure-type`: this test uses a
/// shot day at the plain ICAO default (so the flown trajectory is unaffected) and only
/// changes how the auto-zero SOLVE interprets `--zero-pressure`. The two runs must therefore
/// differ (proving the zero-day reduction actually ran) even though the difference is a
/// second-order effect on the muzzle angle rather than a first-order density change.
#[test]
fn zero_pressure_type_qnh_changes_the_auto_zero_result_vs_absolute() {
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "1500",
        "--pressure",
        "1013.25",
        "--auto-zero",
        "300",
        "--zero-pressure",
        "1030.0",
    ];
    let absolute = run_json(&with_extra(args, &["--zero-pressure-type", "absolute"]));
    let qnh = run_json(&with_extra(args, &["--zero-pressure-type", "qnh"]));
    let v_absolute = absolute["impact_velocity"].as_f64().unwrap();
    let v_qnh = qnh["impact_velocity"].as_f64().unwrap();
    assert_ne!(
        v_absolute, v_qnh,
        "zero-day QNH reduction must change the solved zero angle, and therefore the flown \
         trajectory, versus treating the same zero-day reading as absolute"
    );
}

#[test]
fn qnh_pressure_at_sea_level_matches_absolute_exactly() {
    // At sea level (altitude 0) the QNH reduction is the identity: QNH IS station pressure.
    let args: &[&str] = &[
        "trajectory",
        "-v",
        "2700",
        "-b",
        "0.475",
        "-m",
        "168",
        "-d",
        "0.308",
        "--units",
        "metric",
        "--max-range",
        "300",
        "--ignore-ground-impact",
        "-o",
        "json",
        "--altitude",
        "0",
        "--pressure",
        "1005.0",
    ];
    let absolute = run(args);
    let qnh = run(&with_extra(args, &["--pressure-type", "qnh"]));
    assert!(absolute.status.success());
    assert!(qnh.status.success());
    assert_eq!(
        absolute.stdout, qnh.stdout,
        "at sea level, QNH and absolute must be byte-identical (QNH reduces to itself)"
    );
}

/// MBA-1421: what `--pressure-type` means for a pressure that came from a `--location` CSV
/// rather than from `--pressure`.
///
/// This was recorded as unspecified after the 0.29.0 review, because the pairing has the same
/// shape as two bugs that release fixed (an interpretive MODE reaching a VALUE it does not
/// describe). It turns out the engine already does the defensible thing: the CLI mode applies
/// to the CSV value. Unlike a profile-STORED mode — which describes a value the profile also
/// stored, and which 0.29.0 correctly stopped applying to a CLI `--pressure` — a mode typed on
/// the command line is the user's present-tense declaration about whatever pressure this run
/// uses, and there is exactly one pressure in the run for it to describe.
///
/// These tests exist so that stops being an accident of implementation.
mod csv_supplied_pressure_mode {
    use super::*;
    use std::io::Write;

    /// Writes a location CSV and returns its path. Named per test, matching the temp_dir idiom
    /// the PDF tests already use, so parallel tests cannot collide on one fixture file.
    fn location_csv(test_name: &str, pressure_inhg: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("bx_loc_{test_name}.csv"));
        let mut file = std::fs::File::create(&path).expect("create location csv");
        writeln!(file, "NAME,LATITUDE,ALTITUDE,TEMPERATURE,PRESSURE,HUMIDITY")
            .expect("write header");
        writeln!(file, "TestSite,45.0,5000,59,{pressure_inhg},50").expect("write row");
        path
    }

    fn impact_velocity(args: &[&str]) -> f64 {
        run_json(args)["impact_velocity"]
            .as_f64()
            .expect("impact_velocity")
    }

    fn base() -> Vec<&'static str> {
        vec![
            "trajectory",
            "--bc",
            "0.243",
            "--drag-model",
            "g7",
            "--velocity",
            "2700",
            "-m",
            "175",
            "-d",
            "0.308",
            "--max-range",
            "500",
            "-o",
            "json",
        ]
    }

    /// CONTROL, and the reason this module exists at all. An earlier investigation of this
    /// ticket used `--location-name` (a PDF header display override) instead of `--site` (the
    /// row selector), so the CSV was never loaded, both runs were plain no-CSV runs, and they
    /// matched — which looked exactly like "the mode is inert". Any test that compares two
    /// location runs must first prove the location file changes the answer at all.
    #[test]
    fn the_location_csv_actually_reaches_the_atmosphere() {
        let csv = location_csv("reaches_atmosphere", "24.90");
        let csv = csv.to_str().expect("utf-8 path");

        let mut with_csv = base();
        with_csv.extend_from_slice(&["--location", csv, "--site", "TestSite"]);

        let without = impact_velocity(&base());
        let with = impact_velocity(&with_csv);

        assert!(
            (without - with).abs() > 1.0,
            "the location CSV changed nothing ({without} vs {with} fps) — the fixture is not \
             actually loading, and every comparison below would be vacuous"
        );
    }

    /// The declared mode reaches the CSV-supplied value, and in the correct direction: reading
    /// 24.90 inHg at 5000 ft as a sea-level altimeter setting reduces it to a THINNER station
    /// pressure than taking it as station pressure directly, so drag falls and the projectile
    /// retains more velocity.
    #[test]
    fn a_cli_pressure_mode_applies_to_a_csv_supplied_pressure() {
        let csv = location_csv("mode_applies", "24.90");
        let csv = csv.to_str().expect("utf-8 path");

        let mut absolute = base();
        absolute.extend_from_slice(&[
            "--location",
            csv,
            "--site",
            "TestSite",
            "--pressure-type",
            "absolute",
        ]);
        let mut qnh = base();
        qnh.extend_from_slice(&[
            "--location",
            csv,
            "--site",
            "TestSite",
            "--pressure-type",
            "qnh",
        ]);

        let absolute = impact_velocity(&absolute);
        let qnh = impact_velocity(&qnh);

        assert!(
            qnh > absolute,
            "QNH should reduce the CSV pressure to a thinner station pressure and retain MORE \
             velocity; got qnh {qnh} vs absolute {absolute}"
        );
        assert!(
            (qnh - absolute).abs() > 50.0,
            "the two modes differ by only {:.2} fps on a 5000 ft fixture — too small to \
             distinguish an applied mode from an ignored one",
            (qnh - absolute).abs()
        );
    }

    /// Omitting the mode must equal declaring `absolute`, so a CSV run without the flag is
    /// unchanged from before `--pressure-type` existed.
    #[test]
    fn omitting_the_mode_matches_absolute_on_a_csv_pressure() {
        let csv = location_csv("omitted_matches_absolute", "24.90");
        let csv = csv.to_str().expect("utf-8 path");

        let mut omitted = base();
        omitted.extend_from_slice(&["--location", csv, "--site", "TestSite"]);
        let mut absolute = base();
        absolute.extend_from_slice(&[
            "--location",
            csv,
            "--site",
            "TestSite",
            "--pressure-type",
            "absolute",
        ]);

        assert_eq!(
            impact_velocity(&omitted),
            impact_velocity(&absolute),
            "an omitted --pressure-type must be exactly `absolute` on a CSV pressure too"
        );
    }
}

/// MBA-1416: `--pressure-type` on the ten calculator subcommands that previously accepted
/// `--pressure` and silently treated it as absolute station pressure.
///
/// The risk this closes is inconsistency rather than a wrong number in isolation: someone who
/// learns `--pressure-type` on `trajectory` reasonably assumes `come-ups` honours it too, and a
/// weather-report barometer value entered at elevation is wrong by the ISA reduction on every
/// command that lacks it.
///
/// Two properties per command, and both matter. That the flag REACHES the physics catches the
/// failure where an argument parses but never reaches the air — the exact class this ticket is
/// about. That omitting it is byte-identical to `absolute` protects every existing invocation.
mod pressure_mode_on_calculator_subcommands {
    use super::*;

    /// (name, args). Each set is a minimal run for that subcommand at 5000 ft with a 24.90 inHg
    /// reading, which is a large enough reduction that an applied mode cannot be confused with
    /// numerical noise.
    fn cases() -> Vec<(&'static str, Vec<&'static str>)> {
        let atmo = ["--pressure", "24.90", "--altitude", "5000"];
        let bullet = [
            "--bc", "0.243", "--drag-model", "g7", "--velocity", "2700", "-m", "175", "-d",
            "0.308",
        ];
        let mut out: Vec<(&'static str, Vec<&'static str>)> = Vec::new();
        let mut push = |name: &'static str, extra: Vec<&'static str>, with_bullet: bool| {
            let mut args = vec![name];
            if with_bullet {
                args.extend_from_slice(&bullet);
            }
            args.extend_from_slice(&atmo);
            args.extend(extra);
            out.push((name, args));
        };
        push("mpbr", vec!["--vital-zone", "8"], true);
        push("come-ups", vec!["--zero-distance", "100"], true);
        push("lead", vec!["--target-speed", "10"], true);
        push("wind-card", vec!["--zero-distance", "100"], true);
        push("range-table", vec!["--zero-distance", "100"], true);
        push(
            "stability",
            vec![
                "--twist-rate", "10", "--length", "1.24", "--mass", "175", "--diameter", "0.308",
                "--velocity", "2700",
            ],
            false,
        );
        push(
            "true-velocity",
            vec![
                "--measured-drop", "30", "--range", "500", "--bc", "0.243", "--mass", "175",
                "--diameter", "0.308",
                // --offline: without it this command calls a live API by DEFAULT, which made
                // this suite network-dependent — it went green four CI runs in a row and then
                // failed on an API timeout. A hermetic test may depend on neither $HOME (the
                // ToS lesson) nor the network, and the pressure-mode behaviour under test is
                // entirely in the local solve path.
                "--offline",
            ],
            false,
        );
        push(
            "plan-truing",
            vec![
                "--measurement-resolution", "0.25", "--candidate-ranges", "300,500,700",
                "--velocity", "2700", "--bc", "0.243", "--mass", "175", "--diameter", "0.308",
            ],
            false,
        );
        push(
            "compare",
            vec![
                "--zero-distance", "100", "--load", "A:g7:0.243:2700:175:0.308", "--load",
                "B:g7:0.250:2750:180:0.308",
            ],
            false,
        );
        push(
            "estimate-bc",
            vec![
                "--velocity", "2700", "--mass", "175", "--diameter", "0.308", "--data",
                "300,10;500,30",
            ],
            false,
        );
        out
    }

    /// A private HOME with the terms of service pre-accepted.
    ///
    /// `true-velocity` prompts for ToS acceptance on first use and stores the answer under
    /// `$HOME/.ballistics`. A developer machine has already answered it, so this suite passed
    /// locally and failed on a clean CI runner — the test was reading state from the machine it
    /// ran on. Same fixture as `tests/integration_cli.rs::accepted_tos_home`.
    fn accepted_tos_home() -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU32, Ordering};
        static N: AtomicU32 = AtomicU32::new(0);
        let home = std::env::temp_dir().join(format!(
            "bx-qnh-tos-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let config_dir = home.join(".ballistics");
        std::fs::create_dir_all(&config_dir).expect("create config dir");
        std::fs::write(
            config_dir.join("tos.json"),
            r#"{
  "accepted": true,
  "accepted_version": "2026-01-26",
  "accepted_at": "test",
  "terms_hash": "test"
}"#,
        )
        .expect("write tos.json");
        home
    }

    fn stdout_of(args: &[&str]) -> String {
        let output = Command::new(get_cli_binary())
            .env("HOME", accepted_tos_home())
            .args(args)
            .output()
            .expect("failed to execute ballistics binary");
        assert!(
            output.status.success(),
            "`{}` failed: {}",
            args.join(" "),
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8_lossy(&output.stdout).into_owned()
    }

    /// The declared mode must change the answer. An argument that parses but never reaches the
    /// atmosphere is precisely the defect this ticket exists to close.
    #[test]
    fn qnh_reaches_the_physics_on_every_subcommand() {
        for (name, args) in cases() {
            let mut absolute = args.clone();
            absolute.extend_from_slice(&["--pressure-type", "absolute"]);
            let mut qnh = args.clone();
            qnh.extend_from_slice(&["--pressure-type", "qnh"]);

            assert_ne!(
                stdout_of(&absolute),
                stdout_of(&qnh),
                "`{name} --pressure-type qnh` produced identical output to `absolute`; the flag \
                 parses but is not reaching the atmosphere"
            );
        }
    }

    /// Every pre-existing invocation must be untouched, byte for byte.
    #[test]
    fn omitting_the_mode_is_byte_identical_to_absolute_on_every_subcommand() {
        for (name, args) in cases() {
            let mut absolute = args.clone();
            absolute.extend_from_slice(&["--pressure-type", "absolute"]);

            assert_eq!(
                stdout_of(&args),
                stdout_of(&absolute),
                "`{name}` changed when --pressure-type was omitted; the default must remain \
                 exactly `absolute`"
            );
        }
    }
}
