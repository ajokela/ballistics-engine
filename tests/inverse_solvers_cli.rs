//! MBA-1362: `mark-to-range`, `bdc-match`, and `optimal-zero`.
//!
//! The three inverse solvers share ONE drop-vs-range helper (`HoldCurve` in main.rs:
//! forward interpolation plus a monotone bisection inverse), so the tests here are written
//! to catch the failure that sharing is meant to prevent — a solver disagreeing with the
//! forward model it is supposedly inverting. `mark-to-range` is checked against `come-ups`'
//! own drop table, `bdc-match` against a magnification constructed from that same table,
//! and `optimal-zero` against explicitly-zeroed `come-ups` runs at the bracket endpoints.
//!
//! All three are CLI-only this train (WASM parity is a tracked follow-up per the ticket),
//! so every test drives the binary.

use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

const LOAD: [&str; 8] = ["-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308"];

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(BIN).args(args).output().expect("solver");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

fn json(args: &[&str]) -> serde_json::Value {
    let (stdout, stderr, ok) = run(args);
    assert!(ok, "{args:?} failed: {stderr}");
    serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("{args:?} emitted non-JSON ({e}): {stdout}"))
}

/// Angular drops (mil) at `ranges`, straight from `come-ups`, for a given zero.
///
/// A fixed 100-unit step, NOT the spacing of `ranges`: `come-ups` samples the trajectory at
/// its `--step` interval and snaps each requested range to the nearest sample, so a coarse
/// step silently reports a neighbouring range's drop. Every range asked for here is a
/// multiple of 100 from the start, so each lands on a sample exactly.
fn come_up_drops(zero: &str, ranges: &[u32]) -> Vec<f64> {
    let start = ranges.iter().min().unwrap().to_string();
    let end = ranges.iter().max().unwrap().to_string();
    let step = "100".to_string();
    let mut args: Vec<&str> = vec!["come-ups"];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&[
        "--zero-distance",
        zero,
        "--start",
        &start,
        "--end",
        &end,
        "--step",
        &step,
        "-o",
        "json",
    ]);
    let table = json(&args);
    ranges
        .iter()
        .map(|range| {
            table["data"]
                .as_array()
                .unwrap()
                .iter()
                .find(|row| (row["range"].as_f64().unwrap() - f64::from(*range)).abs() < 1e-9)
                .unwrap_or_else(|| panic!("come-ups produced no row for {range}: {table}"))["drop"]
                .as_f64()
                .unwrap()
        })
        .collect()
}

/// `mark-to-range` inverts the drop table `come-ups` produces: feeding it the come-up at a
/// known range must give that range back.
#[test]
fn mark_to_range_inverts_a_known_drop_table() {
    let ranges = [300u32, 500, 700];
    let drops = come_up_drops("100", &ranges);

    let mark_args: Vec<String> = drops.iter().map(|d| format!("{d}")).collect();
    let mut args: Vec<&str> = vec!["mark-to-range"];
    for mark in &mark_args {
        args.push("--mark");
        args.push(mark);
    }
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100", "-o", "json"]);
    let report = json(&args);

    let marks = report["marks"].as_array().unwrap();
    assert_eq!(marks.len(), ranges.len());
    for (mark, expected) in marks.iter().zip(ranges) {
        assert_eq!(mark["status"], "reached", "{mark}");
        let solved = mark["range"].as_f64().unwrap();
        assert!(
            (solved - f64::from(expected)).abs() < 0.5,
            "mark {} mil should invert to {expected} yd, got {solved}",
            mark["true_mil"]
        );
        // The per-mark ballistics come along and are self-consistent.
        assert!(mark["time"].as_f64().unwrap() > 0.0);
        assert!(mark["velocity"].as_f64().unwrap() > 0.0);
        assert!(mark["energy"].as_f64().unwrap() > 0.0);
    }
}

/// An SFP mark's TRUE subtension is `nominal * reference_mag / mag`, so the same etched
/// mark inverts to a longer range at lower magnification. The conversion is exact.
#[test]
fn mark_to_range_applies_sfp_scaling_to_the_subtension() {
    let at = |mag: &str| {
        let mut args: Vec<&str> = vec![
            "mark-to-range",
            "--mark",
            "2.0",
            "--focal-plane",
            "sfp",
            "--reference-mag",
            "10",
            "--mag",
            mag,
        ];
        args.extend_from_slice(&LOAD);
        args.extend_from_slice(&["--zero-distance", "100", "-o", "json"]);
        json(&args)
    };
    let reference = at("10");
    let half = at("5");
    assert_eq!(reference["marks"][0]["true_mil"], 2.0);
    assert_eq!(half["marks"][0]["true_mil"], 4.0);
    assert_eq!(half["mark_scale"], 2.0);
    // Twice the subtension is a longer range.
    assert!(
        half["marks"][0]["range"].as_f64().unwrap()
            > reference["marks"][0]["range"].as_f64().unwrap()
    );

    // The equivalent FFP mark is the reference case exactly.
    let mut args: Vec<&str> = vec!["mark-to-range", "--mark", "2.0", "--mag", "5"];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100", "-o", "json"]);
    let ffp = json(&args);
    assert_eq!(
        ffp["marks"][0]["range"],
        reference["marks"][0]["range"],
        "FFP subtensions do not depend on magnification"
    );
}

/// Marks the load cannot reach are REPORTED, not dropped: the table keeps one row per
/// supplied mark, in order, whatever the outcome.
#[test]
fn mark_to_range_reports_unreachable_marks_instead_of_dropping_them() {
    let mut args: Vec<&str> = vec![
        "mark-to-range",
        "--mark",
        "-1.0",
        "--mark",
        "2.0",
        "--mark",
        "200.0",
        "--max-range",
        "800",
    ];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100", "-o", "json"]);
    let report = json(&args);
    let marks = report["marks"].as_array().unwrap();
    assert_eq!(marks.len(), 3, "no mark may be silently dropped: {report}");
    assert_eq!(marks[0]["status"], "not_a_holdover");
    assert!(marks[0]["far_zero_range"].as_f64().unwrap() > 0.0);
    assert_eq!(marks[1]["status"], "reached");
    assert_eq!(marks[2]["status"], "beyond_max_range");
    assert!(marks[2]["max_drop_mil"].as_f64().unwrap() > 0.0);
    assert!(marks[2]["max_range"].as_f64().unwrap() >= 800.0);
}

/// `bdc-match` recovers a constructed magnification EXACTLY.
///
/// Build the marks from the load's own solved drops as `s_i = t_i * M / reference_mag`.
/// Then the least-squares slope through the origin is `reference_mag / M` by construction,
/// so the fitted magnification is `M` — with zero residual at every mark.
#[test]
fn bdc_match_recovers_a_constructed_magnification_exactly() {
    let ranges = [300u32, 500, 700];
    let drops = come_up_drops("100", &ranges);
    let reference_mag = 12.0_f64;
    let constructed_mag = 7.5_f64;

    let pairs: Vec<String> = drops
        .iter()
        .zip(ranges)
        .map(|(drop, range)| format!("{}:{range}", drop * constructed_mag / reference_mag))
        .collect();
    let mut args: Vec<&str> = vec!["bdc-match"];
    for pair in &pairs {
        args.push("--mark-range");
        args.push(pair);
    }
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--reference-mag", "12", "--zero-distance", "100", "-o", "json"]);
    let report = json(&args);

    let fitted = report["fitted_magnification"].as_f64().unwrap();
    assert!(
        (fitted - constructed_mag).abs() < 1e-6,
        "fitted {fitted} should recover the constructed {constructed_mag}"
    );
    assert!(
        report["worst_residual_mil"].as_f64().unwrap() < 1e-9,
        "a constructed reticle fits perfectly: {report}"
    );
    assert_eq!(report["fit_warning"], false);
    for mark in report["marks"].as_array().unwrap() {
        assert!(mark["residual_mil"].as_f64().unwrap().abs() < 1e-9, "{mark}");
    }
}

/// A reticle that does not fit the load at any magnification says so rather than
/// presenting the least-bad answer as a calibration.
#[test]
fn bdc_match_warns_when_no_magnification_fits() {
    let mut args: Vec<&str> = vec![
        "bdc-match",
        "--mark-range",
        "2:300",
        "--mark-range",
        "4:500",
        "--mark-range",
        "6:700",
    ];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--reference-mag", "12", "--zero-distance", "100", "-o", "json"]);
    let report = json(&args);
    // Evenly spaced marks cannot follow a curved drop line: some residual must remain.
    assert!(report["worst_residual_mil"].as_f64().unwrap() > 0.2, "{report}");
    assert_eq!(report["fit_warning"], true);

    let (table, _, ok) = run(&args[..args.len() - 2]);
    assert!(ok);
    assert!(table.contains("WARNING"), "{table}");
    assert!(table.contains("least-bad compromise"), "{table}");
}

/// `optimal-zero` beats BOTH bracket endpoints on a three-target case, which is the whole
/// claim the command makes.
#[test]
fn optimal_zero_beats_both_endpoint_zeroes() {
    let ranges = [200u32, 400, 600];
    let mut args: Vec<&str> = vec![
        "optimal-zero",
        "--target",
        "200:10",
        "--target",
        "400:12",
        "--target",
        "600:18",
    ];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["-o", "json"]);
    let report = json(&args);

    let optimal_max = report["max_hold_mil"].as_f64().unwrap();
    let optimal_zero = report["optimal_zero"].as_f64().unwrap();
    assert!(
        optimal_zero > 200.0 && optimal_zero < 600.0,
        "the optimal zero should land between the extremes, got {optimal_zero}"
    );

    for endpoint in ["200", "600"] {
        let worst = come_up_drops(endpoint, &ranges)
            .iter()
            .fold(0.0_f64, |acc, d| acc.max(d.abs()));
        assert!(
            optimal_max < worst,
            "zeroed at {endpoint} the worst hold is {worst} mil; the optimal zero \
             ({optimal_zero}) must beat that, got {optimal_max}"
        );
    }

    // The min-max solution balances the extremes: the two outer targets end up needing
    // (very nearly) the same magnitude of hold, in opposite directions.
    let targets = report["targets"].as_array().unwrap();
    let near = targets[0]["hold_mil"].as_f64().unwrap();
    let far = targets[2]["hold_mil"].as_f64().unwrap();
    assert!(near < 0.0 && far > 0.0, "{report}");
    assert!(
        (near.abs() - far.abs()).abs() < 0.01,
        "a Chebyshev-centred zero equalizes the extremes: {near} vs {far}"
    );
}

/// The fit / vital-zone verdict is real: a generous vital zone fits, a tight one does not,
/// and the answer is per-target.
#[test]
fn optimal_zero_reports_the_vital_zone_verdict() {
    let solve = |vital: &str| {
        let mut args: Vec<&str> = vec![
            "optimal-zero",
            "--target",
            "200",
            "--target",
            "300",
            "--target",
            "400",
            "--vital",
            vital,
        ];
        args.extend_from_slice(&LOAD);
        args.extend_from_slice(&["-o", "json"]);
        json(&args)
    };
    let generous = solve("40");
    assert_eq!(generous["all_targets_fit"], true, "{generous}");
    assert!(generous["targets"]
        .as_array()
        .unwrap()
        .iter()
        .all(|t| t["fits"] == true));

    let tight = solve("1");
    assert_eq!(tight["all_targets_fit"], false, "{tight}");
}

/// Same inputs, same output — byte for byte. The golden-section search uses fixed
/// constants precisely so this holds.
#[test]
fn all_three_solvers_are_deterministic() {
    let cases: Vec<Vec<&str>> = vec![
        {
            let mut a: Vec<&str> = vec!["mark-to-range", "--mark", "2", "--mark", "4", "--mark", "6"];
            a.extend_from_slice(&LOAD);
            a.extend_from_slice(&["--zero-distance", "100", "-o", "json"]);
            a
        },
        {
            let mut a: Vec<&str> = vec![
                "bdc-match",
                "--mark-range",
                "2:300",
                "--mark-range",
                "4:500",
                "--mark-range",
                "6:700",
            ];
            a.extend_from_slice(&LOAD);
            a.extend_from_slice(&["--reference-mag", "12", "--zero-distance", "100", "-o", "json"]);
            a
        },
        {
            let mut a: Vec<&str> = vec![
                "optimal-zero",
                "--target",
                "200:10",
                "--target",
                "450:12",
                "--target",
                "700:18",
            ];
            a.extend_from_slice(&LOAD);
            a.extend_from_slice(&["-o", "json"]);
            a
        },
    ];
    for case in cases {
        let (first, _, ok) = run(&case);
        assert!(ok);
        let (second, _, ok) = run(&case);
        assert!(ok);
        assert_eq!(first, second, "{case:?} is not deterministic");
    }
}

/// `optimal-zero` does not depend on the ORDER the targets were typed in.
#[test]
fn optimal_zero_is_insensitive_to_target_order() {
    let solve = |targets: [&str; 3]| {
        let mut args: Vec<&str> = vec!["optimal-zero"];
        for target in targets {
            args.push("--target");
            args.push(target);
        }
        args.extend_from_slice(&LOAD);
        args.extend_from_slice(&["-o", "json"]);
        json(&args)
    };
    assert_eq!(
        solve(["200:10", "450:12", "700:18"]),
        solve(["700:18", "200:10", "450:12"])
    );
}

/// Every rejection path, with an explanation rather than a silent fallback.
#[test]
fn rejection_cases() {
    // bdc-match on an FFP reticle is meaningless, and says why.
    let mut args: Vec<&str> = vec!["bdc-match", "--mark-range", "2:300", "--mark-range", "4:500", "--focal-plane", "ffp"];
    args.extend_from_slice(&LOAD);
    args.push("--zero-distance");
    args.push("100");
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("SECOND focal plane"), "{stderr}");
    assert!(stderr.contains("mark-to-range"), "the error should name the alternative: {stderr}");

    // One pair cannot constrain a fit.
    let mut args: Vec<&str> = vec!["bdc-match", "--mark-range", "2:300"];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100"]);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("at least two"), "{stderr}");

    // Malformed pair grammar.
    let mut args: Vec<&str> = vec!["bdc-match", "--mark-range", "2", "--mark-range", "4:500"];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100"]);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("MIL:RANGE"), "{stderr}");

    // mark-to-range: the two subtension sources are exclusive, and one is required.
    let mut args: Vec<&str> = vec!["mark-to-range"];
    args.extend_from_slice(&LOAD);
    args.extend_from_slice(&["--zero-distance", "100"]);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("--mark"), "{stderr}");

    // mark-to-range without a zero has nothing to invert against.
    let mut args: Vec<&str> = vec!["mark-to-range", "--mark", "2"];
    args.extend_from_slice(&LOAD);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("--zero-distance"), "{stderr}");

    // optimal-zero target-count bounds.
    let mut args: Vec<&str> = vec!["optimal-zero", "--target", "300:10"];
    args.extend_from_slice(&LOAD);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("2 to 16"), "{stderr}");

    let many: Vec<String> = (1..=17).map(|i| format!("{}:10", i * 50)).collect();
    let mut args: Vec<&str> = vec!["optimal-zero"];
    for target in &many {
        args.push("--target");
        args.push(target);
    }
    args.extend_from_slice(&LOAD);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("2 to 16"), "{stderr}");

    // A target with no height and no --vital cannot be judged.
    let mut args: Vec<&str> = vec!["optimal-zero", "--target", "200", "--target", "400"];
    args.extend_from_slice(&LOAD);
    let (_, stderr, ok) = run(&args);
    assert!(!ok);
    assert!(stderr.contains("--vital"), "{stderr}");

    // CSV and PDF have no form on any of the three.
    for command in ["mark-to-range", "bdc-match", "optimal-zero"] {
        let mut args: Vec<&str> = match command {
            "mark-to-range" => vec!["mark-to-range", "--mark", "2", "--zero-distance", "100"],
            "bdc-match" => vec![
                "bdc-match",
                "--mark-range",
                "2:300",
                "--mark-range",
                "4:500",
                "--zero-distance",
                "100",
            ],
            _ => vec!["optimal-zero", "--target", "200:10", "--target", "400:12"],
        };
        args.extend_from_slice(&LOAD);
        args.extend_from_slice(&["-o", "csv"]);
        let (_, stderr, ok) = run(&args);
        assert!(!ok, "{command} accepted CSV");
        assert!(stderr.contains("no CSV form"), "{command}: {stderr}");
    }
}

/// WASM parity is a tracked follow-up, so none of the three exist in the browser terminal
/// yet. Pinning it keeps the CHANGELOG's claim honest and makes the eventual addition a
/// deliberate edit rather than an accident.
#[test]
fn the_three_solvers_are_native_only_this_train() {
    // Nothing to drive the WASM dispatch from a host test, so assert the native side of
    // the contract instead: they ARE real native subcommands.
    for command in ["mark-to-range", "bdc-match", "optimal-zero"] {
        let (stdout, _, ok) = run(&[command, "--help"]);
        assert!(ok, "{command} --help should succeed");
        assert!(stdout.contains("MBA-1362"), "{command} help: {stdout}");
    }
}
