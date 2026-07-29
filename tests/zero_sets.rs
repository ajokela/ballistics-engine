//! MBA-1360: multiple named zeroes and per-load POI (dial-correction) offsets in saved
//! profiles.
//!
//! Covers the CLI surface end-to-end: `profile zero-set add|remove|list` round-trip,
//! `--zero-set` selection changing dial output by exactly the stored mils (and its
//! composition with the MBA-1358 tracking CF: bias added to the true angular need
//! BEFORE the CF division), the set's zero_distance feeding the auto-zero, unknown-name
//! hard errors listing the available sets, byte-identity when no set is selected, and
//! the `--profile` CSV row's V_OFFSET_MIL/H_OFFSET_MIL columns forming a selectable
//! ephemeral set.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation
//! points `HOME` at a private temp dir (the tests/dsf_profile_cli.rs pattern).

use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn private_home() -> std::path::PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-zero-sets-{}-{}",
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

fn save_base_profile(home: &std::path::Path, name: &str) {
    let out = run(
        home,
        &[
            "profile", "save", name, "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308",
            "--zero-distance", "100", "--auto-zero", "100",
        ],
    );
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
}

fn come_up_rows(json: &str) -> Vec<(f64, f64)> {
    let v: serde_json::Value = serde_json::from_str(json).expect("come-ups JSON");
    v["data"]
        .as_array()
        .expect("data array")
        .iter()
        .map(|row| {
            (
                row["range"].as_f64().unwrap(),
                row["drop"].as_f64().unwrap(),
            )
        })
        .collect()
}

const COME_UPS_ARGS: &[&str] = &[
    "--zero-distance",
    "100",
    "--start",
    "200",
    "--end",
    "500",
    "--step",
    "100",
    "-o",
    "json",
];

/// `profile zero-set add|list|remove` round-trips, upserts by name, and errors on a
/// missing removal target.
#[test]
fn zero_set_management_round_trips() {
    let home = private_home();
    save_base_profile(&home, "mgmt");

    run_ok(
        &home,
        &[
            "profile", "zero-set", "add", "mgmt", "--name", "suppressed", "--zero-distance",
            "200", "--poi-up", "-0.3", "--poi-right", "0.1", "--notes", "suppressed load",
        ],
    );
    run_ok(
        &home,
        &["profile", "zero-set", "add", "mgmt", "--name", "match", "--poi-up", "0.25"],
    );

    let list = run_ok(&home, &["profile", "zero-set", "list", "mgmt"]);
    assert!(list.contains("suppressed"), "{list}");
    assert!(list.contains("zero 200 yd"), "{list}");
    assert!(list.contains("up -0.30 mil"), "{list}");
    assert!(list.contains("right +0.10 mil"), "{list}");
    assert!(list.contains("suppressed load"), "{list}");
    assert!(list.contains("match"), "{list}");
    assert!(list.contains("up +0.25 mil"), "{list}");

    // Upsert by name: replacing announces itself and does not duplicate.
    let out = run(
        &home,
        &["profile", "zero-set", "add", "mgmt", "--name", "match", "--poi-up", "0.5"],
    );
    assert!(out.status.success());
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("Replaced zero set 'match'"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let list = run_ok(&home, &["profile", "zero-set", "list", "mgmt"]);
    assert_eq!(list.matches("match").count(), 1, "{list}");
    assert!(list.contains("up +0.50 mil"), "{list}");

    // Remove; removing again is a hard error naming the remaining sets.
    run_ok(&home, &["profile", "zero-set", "remove", "mgmt", "--name", "match"]);
    let out = run(&home, &["profile", "zero-set", "remove", "mgmt", "--name", "match"]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("available zero sets: suppressed"), "{err}");

    // Removing the last set drops the key from the stored JSON entirely.
    run_ok(&home, &["profile", "zero-set", "remove", "mgmt", "--name", "suppressed"]);
    let stored = std::fs::read_to_string(home.join(".ballistics/profiles/mgmt.json")).unwrap();
    assert!(!stored.contains("zero_sets"), "{stored}");
}

/// Selecting a set changes the come-ups dial column by EXACTLY the stored mils, and a
/// profile that merely CARRIES sets (none selected) is byte-identical to one without.
#[test]
fn zero_set_selection_shifts_dials_by_exactly_the_stored_mils() {
    let home = private_home();
    save_base_profile(&home, "sel");
    save_base_profile(&home, "plain");
    run_ok(
        &home,
        &["profile", "zero-set", "add", "sel", "--name", "match", "--poi-up", "0.25"],
    );

    let mut args = vec!["come-ups", "--profile", "sel"];
    args.extend_from_slice(COME_UPS_ARGS);
    let base = run_ok(&home, &args);

    // Byte-identity: carrying zero sets changes nothing while none is selected.
    let mut plain_args = vec!["come-ups", "--profile", "plain"];
    plain_args.extend_from_slice(COME_UPS_ARGS);
    assert_eq!(base, run_ok(&home, &plain_args));

    let mut sel_args = vec!["come-ups", "--profile", "sel", "--zero-set", "match"];
    sel_args.extend_from_slice(COME_UPS_ARGS);
    let selected = run_ok(&home, &sel_args);

    let base_rows = come_up_rows(&base);
    let sel_rows = come_up_rows(&selected);
    assert_eq!(base_rows.len(), sel_rows.len());
    for ((r0, d0), (r1, d1)) in base_rows.iter().zip(&sel_rows) {
        assert_eq!(r0, r1);
        assert!(
            ((d1 - d0) - 0.25).abs() < 1e-9,
            "at {r0} yd: base {d0}, selected {d1} — expected exactly +0.25 mil"
        );
    }
}

/// The MBA-1358 composition pin at the CLI level: the bias is added to the TRUE angular
/// need BEFORE the CF division, so selected output == (base + 0.25) / 0.95 exactly.
#[test]
fn zero_set_bias_is_added_before_the_tracking_cf_division() {
    let home = private_home();
    save_base_profile(&home, "cfsel");
    run_ok(
        &home,
        &["profile", "zero-set", "add", "cfsel", "--name", "match", "--poi-up", "0.25"],
    );

    let mut base_args = vec!["come-ups", "--profile", "cfsel"];
    base_args.extend_from_slice(COME_UPS_ARGS);
    let base = come_up_rows(&run_ok(&home, &base_args));

    let mut cf_args = vec![
        "come-ups", "--profile", "cfsel", "--zero-set", "match", "--elevation-cf", "0.95",
    ];
    cf_args.extend_from_slice(COME_UPS_ARGS);
    let corrected = come_up_rows(&run_ok(&home, &cf_args));

    for ((r0, d0), (_, d1)) in base.iter().zip(&corrected) {
        let expected = (d0 + 0.25) / 0.95;
        assert!(
            (d1 - expected).abs() < 1e-9,
            "at {r0} yd: got {d1}, want (true + bias)/cf = {expected}"
        );
        // Pin the ORDER, not just a formula: the wrong order differs measurably.
        let wrong = d0 / 0.95 + 0.25;
        assert!((d1 - wrong).abs() > 1e-4, "at {r0} yd: bias must precede the CF division");
    }
}

/// A set's zero_distance feeds the trajectory auto-zero exactly as the legacy field
/// would; an explicit --auto-zero still wins.
#[test]
fn zero_set_distance_feeds_the_auto_zero() {
    let home = private_home();
    save_base_profile(&home, "dist");
    run_ok(
        &home,
        &["profile", "zero-set", "add", "dist", "--name", "far", "--zero-distance", "300"],
    );

    let zero_angle = |extra: &[&str]| -> f64 {
        let mut args = vec!["trajectory", "--saved-profile", "dist", "--max-range", "600"];
        args.extend_from_slice(extra);
        args.extend_from_slice(&["-o", "json"]);
        let v: serde_json::Value = serde_json::from_str(&run_ok(&home, &args)).unwrap();
        v["zero_angle_degrees"].as_f64().expect("zero angle echo")
    };

    let master = zero_angle(&[]);
    let far = zero_angle(&["--zero-set", "far"]);
    let explicit300 = zero_angle(&["--auto-zero", "300"]);
    let flag_wins = zero_angle(&["--zero-set", "far", "--auto-zero", "100"]);

    assert!(far > master, "300 yd zero must dial above the 100 yd master zero");
    assert!((far - explicit300).abs() < 1e-12, "set zero == explicit --auto-zero 300");
    assert!((flag_wins - master).abs() < 1e-12, "explicit --auto-zero must beat the set");
}

/// Unknown set names are a hard error listing the available names; --zero-set with no
/// profile at all is a hard error naming the requirement (MBA-1425: never silent).
#[test]
fn unknown_zero_set_is_a_hard_error_listing_names() {
    let home = private_home();
    save_base_profile(&home, "errs");
    run_ok(
        &home,
        &["profile", "zero-set", "add", "errs", "--name", "a", "--poi-up", "0.1"],
    );
    run_ok(
        &home,
        &["profile", "zero-set", "add", "errs", "--name", "b", "--poi-right", "0.2"],
    );

    let mut args = vec!["come-ups", "--profile", "errs", "--zero-set", "nope"];
    args.extend_from_slice(COME_UPS_ARGS);
    let out = run(&home, &args);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("zero set 'nope' not found"), "{err}");
    assert!(err.contains("available zero sets: a, b"), "{err}");

    let out = run(
        &home,
        &[
            "trajectory", "-v", "2700", "-b", "0.475", "-m", "168", "-d", "0.308",
            "--zero-set", "x", "--max-range", "400",
        ],
    );
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("requires a profile"), "{err}");
}

/// The `--profile` CSV row's V_OFFSET_MIL/H_OFFSET_MIL columns (allowlisted since
/// MBA-614, unused until now) form an ephemeral zero set named after the row: it is
/// selectable, listed in the unknown-name error, and absent columns produce no set.
#[test]
fn profile_csv_offsets_form_a_selectable_set() {
    let home = private_home();
    let csv = home.join("gun.csv");
    std::fs::write(
        &csv,
        "BULLET,VELOCITY,BC,ZERO_RANGE,V_OFFSET_MIL,H_OFFSET_MIL\n\
         R1,2650,0.450,200,0.3,-0.1\n\
         R2,2600,0.430,150,,\n",
    )
    .unwrap();
    let csv = csv.to_str().unwrap();

    // Selecting the row-named set succeeds (the biases ride the PDF dial columns; here
    // we pin selection + resolution, which must not error).
    let out = run(
        &home,
        &[
            "trajectory", "--profile", csv, "--profile-row", "R1", "-m", "175", "-d", "0.308",
            "--zero-set", "R1", "--max-range", "400", "-o", "json",
        ],
    );
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    // Unknown name: the CSV-derived set shows up in the available list.
    let out = run(
        &home,
        &[
            "trajectory", "--profile", csv, "--profile-row", "R1", "-m", "175", "-d", "0.308",
            "--zero-set", "bogus", "--max-range", "400",
        ],
    );
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("available zero sets: R1"), "{err}");

    // A row without the offset columns has no set to select.
    let out = run(
        &home,
        &[
            "trajectory", "--profile", csv, "--profile-row", "R2", "-m", "175", "-d", "0.308",
            "--zero-set", "R2", "--max-range", "400",
        ],
    );
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("has no zero sets") || err.contains("requires a profile"), "{err}");
}

/// wind-card and range-table inherit the windage/elevation corrections through the same
/// shared boundary: selecting a windage-only set shifts every wind-card drift dial by
/// exactly the stored mils.
#[test]
fn wind_card_drift_shifts_by_the_windage_correction() {
    let home = private_home();
    save_base_profile(&home, "wc");
    run_ok(
        &home,
        &["profile", "zero-set", "add", "wc", "--name", "load2", "--poi-right", "0.2"],
    );

    let card = |extra: &[&str]| -> String {
        let mut args = vec!["wind-card", "--profile", "wc"];
        args.extend_from_slice(extra);
        args.extend_from_slice(&[
            "--zero-distance", "100", "--start", "300", "--end", "300", "--step", "100",
            "--wind-speeds", "10", "-o", "json",
        ]);
        run_ok(&home, &args)
    };

    let base: serde_json::Value = serde_json::from_str(&card(&[])).unwrap();
    let sel: serde_json::Value =
        serde_json::from_str(&card(&["--zero-set", "load2"])).unwrap();
    let get = |v: &serde_json::Value| -> f64 {
        // wind-card JSON rows carry one `wind_<speed>` drift value per speed column.
        v["data"].as_array().unwrap()[0]["wind_10"].as_f64().unwrap()
    };
    let d0 = get(&base);
    let d1 = get(&sel);
    assert!(
        ((d1 - d0) - 0.2).abs() < 1e-9,
        "wind-card drift must shift by exactly +0.2 mil (base {d0}, selected {d1})"
    );
}
