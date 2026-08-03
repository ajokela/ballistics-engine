//! MBA-1348 review fix I1: the twelve turret/hold `profile save` fields (clicks-per-rev,
//! zero-stop, elevation/windage travel, dialed turret state, reticle hold bounds) — and,
//! as a necessary consequence, `--elevation-click`/`--windage-click` themselves — carry
//! forward across an unrelated re-save, exactly like DSF/CF/reticle/zero-sets already do.
//! `tests/reticle_cli.rs`'s `profile_attachment_round_trips_and_carries_forward` is the
//! precedent this mirrors: an explicit flag always wins; an omitted one carries the
//! existing profile's stored value forward; `--clear-turret` clears all twelve turret/hold
//! fields at once (but never the click graduation, which has no `--clear` of its own).
//!
//! This carry-forward logic lives in the `profile save` CLI handler itself (`main.rs`'s
//! `ProfileAction::Save` match arm), not in any unit-testable pure function — `optic_profile`
//! and the `require_angular_pair`/`require_hold_bounds` helpers it calls are private to the
//! binary crate and never see `--clear-turret` or an existing on-disk profile at all, only
//! the already-merged `ProfileData`. So these tests run the actual compiled binary end to
//! end (the `reticle_cli.rs`/`dsf_profile_cli.rs` pattern), never a `ProfileData` literal.
//!
//! Hermetic: profiles live under `$HOME/.ballistics/profiles`, so every invocation points
//! `HOME` at a private per-test temp dir.

use std::path::Path;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_ballistics");

fn tempfile_dir(tag: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "turret-carry-forward-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The profile-field JSON keys `--clear-turret` clears — exactly the twelve `ProfileData`
/// fields declared alongside `reticle` in `main.rs` (NOT `elevation_click`/`windage_click`,
/// which are a separate pair with their own, flag-less carry-forward contract).
const TWELVE_KEYS: [&str; 12] = [
    "clicks_per_revolution",
    "zero_stop",
    "elevation_travel_up_mil",
    "elevation_travel_down_mil",
    "windage_travel_left_mil",
    "windage_travel_right_mil",
    "turret_elevation_dialed_mil",
    "turret_windage_dialed_mil",
    "hold_bound_up_mil",
    "hold_bound_down_mil",
    "hold_bound_left_mil",
    "hold_bound_right_mil",
];

/// One consistent, fully-valid baseline covering all twelve fields plus the click
/// graduation they depend on — distinct, nonzero values on every axis (so a swapped-axis
/// carry-forward bug would show up as a wrong number, not a coincidentally-matching one).
fn baseline_turret_flags() -> Vec<&'static str> {
    vec![
        "--elevation-click",
        "0.1mil",
        "--clicks-per-rev",
        "12",
        "--zero-stop",
        "true",
        "--travel-up",
        "28mil",
        "--travel-down",
        "5mil",
        "--windage-travel-left",
        "10mil",
        "--windage-travel-right",
        "11mil",
        "--turret-elev",
        "3.5mil",
        "--turret-wind",
        "-0.2mil",
        "--hold-up",
        "4mil",
        "--hold-down",
        "12mil",
        "--hold-left",
        "6mil",
        "--hold-right",
        "6.5mil",
    ]
}

/// `velocity` is a separate parameter (not baked into a fixed base-args vector) so a
/// "re-save that changes only --velocity" test can vary it without passing `-v` twice —
/// clap rejects a required option repeated in one invocation.
fn save(home: &Path, name: &str, velocity: &str, extra: &[&str]) {
    let mut args: Vec<&str> = vec![
        "profile", "save", name, "-v", velocity, "-b", "0.475", "-m", "175", "-d", "0.308",
    ];
    args.extend_from_slice(extra);
    let output = Command::new(BIN)
        .env("HOME", home)
        .args(&args)
        .output()
        .expect("run profile save");
    assert!(
        output.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Like `save`, but asserts the invocation FAILS and returns its stderr.
fn save_expect_failure(home: &Path, name: &str, extra: &[&str]) -> String {
    let mut args: Vec<&str> = vec![
        "profile", "save", name, "-v", "2700", "-b", "0.475", "-m", "175", "-d", "0.308",
    ];
    args.extend_from_slice(extra);
    let output = Command::new(BIN)
        .env("HOME", home)
        .args(&args)
        .output()
        .expect("run profile save");
    assert!(
        !output.status.success(),
        "expected `profile save {extra:?}` to fail, but it exited 0: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn stored(home: &Path, name: &str) -> serde_json::Value {
    serde_json::from_str(
        &std::fs::read_to_string(
            home.join(".ballistics").join("profiles").join(format!("{name}.json")),
        )
        .unwrap(),
    )
    .unwrap()
}

/// (a) A re-save touching only `--velocity` must carry all twelve turret/hold fields (and
/// the click graduation they depend on) forward unchanged, AND must itself succeed — which
/// is the only externally-observable proof that `optic_profile()` + `validate()` still
/// assembled cleanly from the carried-forward data, since `profile save` runs that check
/// before ever writing the file (`optic_profile`/`OpticProfile` are private to the binary
/// crate, unreachable from an external integration test any other way).
#[test]
fn turret_fields_carry_forward_on_an_unrelated_resave() {
    let home = tempfile_dir("a");
    save(&home, "rifle", "2700", &baseline_turret_flags());
    let baseline = stored(&home, "rifle");

    save(&home, "rifle", "2705", &[]);
    let after = stored(&home, "rifle");

    assert_eq!(after["velocity"], 2705.0);
    assert_eq!(
        after["elevation_click"], baseline["elevation_click"],
        "click graduation must carry forward too -- the twelve fields depend on it"
    );
    for key in TWELVE_KEYS {
        assert_eq!(
            after[key], baseline[key],
            "field {key} did not carry forward: before={baseline} after={after}"
        );
    }
}

/// (b) One turret flag given on a re-save updates that field; the other eleven (a
/// DIFFERENT field on each axis than the one being changed) carry forward unchanged.
#[test]
fn one_turret_flag_updates_while_the_rest_carry_forward() {
    let home = tempfile_dir("b");
    save(&home, "rifle", "2700", &baseline_turret_flags());
    let baseline = stored(&home, "rifle");

    save(&home, "rifle", "2700", &["--travel-up", "30mil"]);
    let after = stored(&home, "rifle");

    assert_eq!(after["elevation_travel_up_mil"], 30.0, "the given flag must win");
    for key in TWELVE_KEYS {
        if key == "elevation_travel_up_mil" {
            continue;
        }
        assert_eq!(
            after[key], baseline[key],
            "field {key} unexpectedly changed on an unrelated flag: {after}"
        );
    }
    assert_eq!(after["elevation_click"], baseline["elevation_click"]);
}

/// (c) `--clear-turret` removes all twelve fields at once but leaves the click
/// graduation alone (it isn't one of the twelve), the profile still saves and loads, and
/// clearing does NOT trip the "eleven fields need a click" gate -- nothing is left set
/// for that gate to fire on.
#[test]
fn clear_turret_removes_all_twelve_but_keeps_the_click() {
    let home = tempfile_dir("c");
    save(&home, "rifle", "2700", &baseline_turret_flags());

    save(&home, "rifle", "2700", &["--clear-turret"]);
    let after = stored(&home, "rifle");

    for key in TWELVE_KEYS {
        assert!(
            after.get(key).is_none(),
            "field {key} should be absent after --clear-turret: {after}"
        );
    }
    assert_eq!(after["elevation_click"], "0.1mil", "click graduation is not one of the twelve");

    let output = Command::new(BIN)
        .env("HOME", &home)
        .args(["profile", "show", "rifle"])
        .output()
        .expect("run profile show");
    assert!(
        output.status.success(),
        "profile show failed after --clear-turret: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// (d) `--clear-turret` combined with any of the twelve turret flags in the same
/// invocation is a usage error (ambiguous intent), naming both flags -- the same
/// treatment `--reticle-json`+`--clear-reticle` already gets.
#[test]
fn clear_turret_conflicts_with_any_turret_flag() {
    let home = tempfile_dir("d");
    save(&home, "rifle", "2700", &baseline_turret_flags());

    let err = save_expect_failure(&home, "rifle", &["--clear-turret", "--travel-up", "10mil"]);
    assert!(err.contains("--clear-turret"), "{err}");
    assert!(err.contains("--travel-up"), "{err}");

    // Multiple conflicting flags at once are all named, not just the first found.
    let err = save_expect_failure(
        &home,
        "rifle",
        &["--clear-turret", "--hold-up", "1mil", "--zero-stop", "true"],
    );
    assert!(err.contains("--clear-turret"), "{err}");
    assert!(err.contains("--hold-up"), "{err}");
    assert!(err.contains("--zero-stop"), "{err}");
}
