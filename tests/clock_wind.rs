//! MBA-1367: clock-position wind direction entry with explicit marker syntax.
//!
//! CLI-level pins: marked clock forms (`3oc`, `10h30`, standalone-only `10:30`) map to
//! the same solves as their degree equivalents; bare numbers stay degrees everywhere;
//! `--wind-segment 10:30:400` keeps its historical SPEED:ANGLE:DIST meaning while
//! `10:3oc:400` parses the marked form; invalid clock tokens fail loudly with helpful
//! messages; and the `!= 0.0` sentinel is gone — an explicit `--wind-direction 0` (or
//! `12oc`) now wins over a location CSV's WIND_DIR instead of being silently dropped
//! (the enumerated, documented behavior fix).

use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn private_home() -> std::path::PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-clock-wind-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn run(home: &std::path::Path, args: &[&str]) -> Output {
    Command::new(BIN)
        .env("HOME", home)
        .args(args)
        .output()
        .expect("spawn ballistics")
}

fn run_ok(home: &std::path::Path, args: &[&str]) -> String {
    let out = run(home, args);
    assert!(
        out.status.success(),
        "command {:?} failed: {}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

const LOAD: &[&str] = &["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];

fn trajectory_json(home: &std::path::Path, extra: &[&str]) -> String {
    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&["--max-range", "600", "--wind-speed", "10"]);
    args.extend_from_slice(extra);
    args.extend_from_slice(&["-o", "json"]);
    run_ok(home, &args)
}

/// Every marked clock form solves byte-identically to its degree equivalent, on both
/// the trajectory and come-ups surfaces; bare degrees are untouched.
#[test]
fn clock_forms_match_their_degree_equivalents() {
    let home = private_home();

    for (clock, degrees) in [
        ("3oc", "90"),
        ("6oc", "180"),
        ("9oc", "270"),
        ("10h30", "315"),
        ("10:30", "315"),
    ] {
        assert_eq!(
            trajectory_json(&home, &["--wind-direction", clock]),
            trajectory_json(&home, &["--wind-direction", degrees]),
            "trajectory {clock} != {degrees} deg"
        );
    }

    // Bare numbers stay degrees (identity with themselves formatted differently).
    assert_eq!(
        trajectory_json(&home, &["--wind-direction", "90"]),
        trajectory_json(&home, &["--wind-direction", "90.0"]),
    );

    // come-ups goes through the same clap parser.
    let come_ups = |dir: &str| -> String {
        let mut args = vec!["come-ups"];
        args.extend_from_slice(LOAD);
        args.extend_from_slice(&[
            "--zero-distance", "100", "--start", "200", "--end", "400", "--step", "100",
            "--wind-speed", "10", "--wind-direction", dir, "-o", "json",
        ]);
        run_ok(&home, &args)
    };
    assert_eq!(come_ups("9oc"), come_ups("270"));
}

/// Inside --wind-segment the ANGLE token takes colon-FREE marked forms only, and the
/// historical three-number form keeps its meaning: `10:30:400` is SPEED=10 ANGLE=30.
#[test]
fn wind_segment_angle_clock_forms_and_legacy_numeric_meaning() {
    let home = private_home();

    assert_eq!(
        trajectory_json(&home, &["--wind-segment", "10:3oc:400"]),
        trajectory_json(&home, &["--wind-segment", "10:90:400"]),
    );
    assert_eq!(
        trajectory_json(&home, &["--wind-segment", "10:10h30:400"]),
        trajectory_json(&home, &["--wind-segment", "10:315:400"]),
    );
    // 10:30:400 is still SPEED:ANGLE:DIST — equal to the explicit 30-degree segment
    // and NOT to the 10h30 (315 deg) clock segment.
    assert_eq!(
        trajectory_json(&home, &["--wind-segment", "10:30:400"]),
        trajectory_json(&home, &["--wind-segment", "10:30.0:400"]),
    );
    assert_ne!(
        trajectory_json(&home, &["--wind-segment", "10:30:400"]),
        trajectory_json(&home, &["--wind-segment", "10:10h30:400"]),
    );

    // Clock-form range errors surface through the segment parser's String boundary.
    let mut args = vec!["trajectory"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&["--max-range", "600", "--wind-segment", "10:13oc:400"]);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("hour must be 1-12"), "{err}");
}

/// Invalid clock tokens on the standalone flag are rejected at parse time with the
/// helpful typed messages.
#[test]
fn invalid_clock_tokens_fail_loudly() {
    let home = private_home();
    for (bad, needle) in [
        ("13oc", "hour must be 1-12"),
        ("0oc", "hour must be 1-12"),
        ("3h60", "minutes must be 00-59"),
        ("3:60", "minutes must be 00-59"),
        ("bogus", "expected degrees"),
    ] {
        let mut args = vec!["trajectory"];
        args.extend_from_slice(LOAD);
        args.extend_from_slice(&["--max-range", "400", "--wind-direction", bad]);
        let out = run(&home, &args);
        assert!(!out.status.success(), "{bad} unexpectedly accepted");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(err.contains(needle), "{bad}: {err}");
    }
}

/// The sentinel fix (MBA-1367): presence, not `!= 0.0`, decides whether the CLI wind
/// direction beats a location CSV's WIND_DIR.
///
/// Unchanged: omitting the flag still inherits the CSV value. FIXED (documented
/// behavior change): an explicit `--wind-direction 0` — and `12oc`, which maps to 0 —
/// now WINS over the CSV instead of being silently replaced by it.
#[test]
fn explicit_zero_direction_beats_the_location_csv() {
    let home = private_home();
    let csv = home.join("loc.csv");
    // No atmosphere columns: only WIND_DIR differs from a CSV-less run.
    std::fs::write(&csv, "SITE,WIND_DIR\nS1,270\n").unwrap();
    let csv = csv.to_str().unwrap();

    let with_csv =
        |extra: &[&str]| trajectory_json(&home, &[&["--location", csv, "--site", "S1"], extra].concat());

    // Omitted flag: CSV wind direction still inherited (byte-identical to explicit 270).
    assert_eq!(with_csv(&[]), trajectory_json(&home, &["--wind-direction", "270"]));

    // Explicit 0 now wins over the CSV: identical to a CSV-less headwind run...
    let explicit_zero = with_csv(&["--wind-direction", "0"]);
    assert_eq!(explicit_zero, trajectory_json(&home, &["--wind-direction", "0"]));
    // ...and demonstrably NOT the CSV's 270.
    assert_ne!(explicit_zero, trajectory_json(&home, &["--wind-direction", "270"]));

    // 12 o'clock == 0 degrees survives the same path (the sentinel would have
    // silently dropped it).
    assert_eq!(with_csv(&["--wind-direction", "12oc"]), explicit_zero);
}

/// `profile save --wind-direction 6oc` stores plain degrees in the profile.
#[test]
fn profile_save_stores_clock_entry_as_degrees() {
    let home = private_home();
    let mut args = vec!["profile", "save", "clocked"];
    args.extend_from_slice(LOAD);
    args.extend_from_slice(&["--wind-speed", "10", "--wind-direction", "6oc"]);
    let out = run(&home, &args);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    let stored =
        std::fs::read_to_string(home.join(".ballistics/profiles/clocked.json")).unwrap();
    let v: serde_json::Value = serde_json::from_str(&stored).unwrap();
    assert_eq!(v["wind_direction"].as_f64(), Some(180.0));
}
