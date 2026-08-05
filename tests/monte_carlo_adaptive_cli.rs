//! MBA-1352 (0.33.0 decision-support Plan C Task 5): the `monte-carlo` CLI surface over Task
//! 4's adaptive driver (`run_monte_carlo_adaptive_seeded`) and Task 2/3's Wilson/confidence-
//! sequence statistics (`mc_stats.rs`).
//!
//! Fixture: `-v 800 -b 0.5 -m 10 -d 7.62` under `--units metric` is
//! `BallisticInputs::default()`'s own ballistic parameters (velocity/BC/mass/diameter), the
//! same bullet `cli_api.rs`'s own `monte_carlo_seeded_tests::seeded_test_fixture` uses --
//! `monte-carlo` has no `--drag-model` flag (G1 only, unlike `adaptive-card`'s borrowed .308
//! G7 fixture), so this reuses the library's own default bullet instead of inventing one.
//! `--target-distance 300` matches that module's `loose_params()` target plane. Every adaptive
//! run below pins a small, explicit `--min-samples`/`--mc-batch-size`/`--max-samples` (never
//! the library's 1000-floor `McConvergence::default()` -- Task 4's own report clocked that at
//! 1000+ unoptimised trajectory solves per run) and a fixed `--seed`, so the whole file stays
//! fast: every number asserted below (samples/attempts/arrivals/stop_reason) was read back from
//! the real binary before being written here, not guessed.

use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

/// The shared load fixture, as CLI args: `BallisticInputs::default()`'s own bullet at 800 m/s.
const LOAD_ARGS: &[&str] = &["-v", "800", "-b", "0.5", "-m", "10", "-d", "7.62"];

fn run(tail: &[&str]) -> (String, String, bool) {
    let mut args: Vec<&str> = vec!["--units", "metric", "monte-carlo"];
    args.extend_from_slice(LOAD_ARGS);
    args.extend_from_slice(tail);
    let output = Command::new(bin()).args(&args).output().expect("run ballistics");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

// ---- (a) default (fixed-count) run gains the additive Wilson CI line, and --seed makes the
// whole run -- new line included -- reproducible byte-for-byte. `-n 50` (far below the 1000
// default) keeps this fast; the CI line's presence does not depend on sample count. ----

#[test]
fn default_run_gains_the_ci_line_and_is_reproducible_with_a_seed() {
    let tail = ["-n", "50", "--seed", "1352", "--target-distance", "300"];

    let (out1, err1, ok1) = run(&tail);
    assert!(ok1, "stderr: {err1}");
    assert!(out1.contains("Hit Probability:"), "point-estimate line missing: {out1}");
    assert!(
        out1.contains("Hit probability 95% CI: ["),
        "missing the additive Wilson CI line: {out1}"
    );
    assert!(out1.contains("(Wilson, n=50)"), "{out1}");
    // The CI line reads UNDER the point-estimate line, not interleaved with the box.
    let hp_idx = out1.find("Hit Probability:").unwrap();
    let ci_idx = out1.find("Hit probability 95% CI:").unwrap();
    assert!(ci_idx > hp_idx, "CI line must come after the point estimate: {out1}");

    // Same --seed, same everything else -> identical stdout, including the new line: the
    // numeric estimates are exactly as reproducible as they were before this line existed.
    let (out2, err2, ok2) = run(&tail);
    assert!(ok2, "stderr: {err2}");
    assert_eq!(out1, out2, "same --seed must reproduce stdout byte-for-byte");
}

/// A run with no `--target-distance` has no hit-probability line at all today, and must still
/// have none -- the CI line is conditioned on the exact same thing the point estimate is.
#[test]
fn no_target_distance_means_no_hit_probability_and_no_ci_line() {
    let (out, err, ok) = run(&["-n", "20", "--seed", "1"]);
    assert!(ok, "stderr: {err}");
    assert!(!out.contains("Hit Probability"), "{out}");
    assert!(!out.contains("CI:"), "{out}");
}

/// The fixed-count JSON's additive `hit_probability_ci` key: present (as explicit `null`,
/// matching `hit_probability`'s own null-not-absent convention) when there is no target
/// distance, and a fully-populated, internally-consistent object when there is.
#[test]
fn fixed_count_json_carries_the_additive_hit_probability_ci_object() {
    let (stdout, stderr, ok) = run(&["-n", "30", "--seed", "3", "-o", "json"]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert!(v["hit_probability"].is_null());
    assert!(
        v.as_object().unwrap().contains_key("hit_probability_ci"),
        "hit_probability_ci must be present (even if null), not absent: {stdout}"
    );
    assert!(v["hit_probability_ci"].is_null());

    let (stdout, stderr, ok) =
        run(&["-n", "30", "--seed", "3", "--target-distance", "300", "-o", "json"]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    let p_hat = v["hit_probability"].as_f64().expect("hit_probability");
    let ci = &v["hit_probability_ci"];
    assert_eq!(ci["method"], "wilson_fixed_n");
    assert_eq!(ci["confidence_percent"], 95);
    assert_eq!(ci["samples"], 30);
    let lo = ci["low"].as_f64().expect("low");
    let hi = ci["high"].as_f64().expect("high");
    assert!(lo <= p_hat && p_hat <= hi, "p_hat {p_hat} not inside [{lo}, {hi}]");
    assert!((0.0..=1.0).contains(&lo) && (0.0..=1.0).contains(&hi));
}

// ---- (b) --adaptive -o json: the full report, every field covered (not just the brief's
// minimum four) -- attempts/arrivals included, since the wire grew those two fields after the
// brief was written (Task 4 fix round, immediately after `samples`). Numbers below
// (samples/attempts/arrivals/stop_reason) were read back from a real run at this exact seed
// before being pinned; see the module doc comment. ----

#[test]
fn adaptive_json_report_parses_and_covers_every_field() {
    let (stdout, stderr, ok) = run(&[
        "--target-distance", "300",
        "--adaptive",
        "--seed", "99",
        "--min-samples", "20",
        "--mc-batch-size", "20",
        "--max-samples", "500",
        "--target-ci-half-width", "0.3",
        "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");

    // schema / method / assumptions -- method and assumptions[0]/[1] are quoted verbatim in
    // CLI_USAGE.md, so pinning them here also keeps the docs honest against a future wording
    // change (a diff would fail this test before it could go stale in the docs).
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["method"], "anytime_beta_binomial_mixture_cs_v1");
    let assumptions = v["assumptions"].as_array().expect("assumptions array");
    assert_eq!(assumptions.len(), 4);
    assert_eq!(
        assumptions[0],
        "Sampling uncertainty only: intervals cover Monte Carlo sampling error, not model error in the trajectory solver or its inputs."
    );
    assert_eq!(
        assumptions[1],
        "Anytime-valid stopping: the beta-binomial mixture confidence sequence keeps its coverage guarantee despite stopping the moment the target half-width is met."
    );
    assert_eq!(
        assumptions[2],
        "Input dispersions are the independent normal distributions declared in MonteCarloParams; correlations between inputs are not modeled."
    );
    assert_eq!(
        assumptions[3],
        "Continuous statistics are streaming Welford moments over trials that reached the target plane, reported with sample (n-1) standard deviations; hit probability's denominator includes all trials."
    );

    // confidence / estimate / interval
    assert_eq!(v["confidence_percent"], 95);
    let p_hat = v["hit_probability"].as_f64().expect("hit_probability");
    let ci_low = v["ci_low"].as_f64().expect("ci_low");
    let ci_high = v["ci_high"].as_f64().expect("ci_high");
    assert!(ci_low <= p_hat && p_hat <= ci_high, "{ci_low} <= {p_hat} <= {ci_high}");
    assert!((0.0..=1.0).contains(&ci_low) && (0.0..=1.0).contains(&ci_high));

    // THE WIRE GREW: samples/attempts/arrivals, in that declaration order, immediately
    // following each other. `attempts >= samples >= arrivals` always; this fixture's stopping
    // controls force attempts == samples == 20 deterministically (checked once, at the first
    // batch boundary), and a real target-plane shortfall at this seed makes arrivals strictly
    // less than samples -- exercising the exact gap MBA-1352's accumulated decisions call out
    // (attempts - samples = drop count, samples - arrivals = shortfall count).
    assert_eq!(v["samples"], 20);
    assert_eq!(v["attempts"], 20);
    let arrivals = v["arrivals"].as_u64().expect("arrivals");
    assert!(arrivals > 0, "arrivals must be positive for this fixture");
    assert!(arrivals < 20, "this fixture must have a real target-plane shortfall (arrivals < samples)");
    let (attempts, samples) = (v["attempts"].as_u64().unwrap(), v["samples"].as_u64().unwrap());
    assert!(attempts >= samples && samples >= arrivals, "ordering invariant violated");

    assert_eq!(v["stop_reason"], "target_half_width_met");
    assert_eq!(v["hit_radius_m"], 0.3);
    assert_eq!(v["target_distance_m"], 300.0);

    // Six Welford stats: present, finite, and the four "at target" fields show real dispersion
    // (arrivals > 1 above guarantees a defined Bessel-corrected sd).
    for key in [
        "mean_impact_velocity_mps",
        "std_impact_velocity_mps",
        "mean_drop_at_target_m",
        "std_drop_at_target_m",
        "mean_wind_drift_at_target_m",
        "std_wind_drift_at_target_m",
    ] {
        assert!(v[key].as_f64().expect(key).is_finite(), "{key} must be finite: {stdout}");
    }
    assert!(v["mean_impact_velocity_mps"].as_f64().unwrap() > 0.0);
    assert!(v["std_impact_velocity_mps"].as_f64().unwrap() > 0.0);
    assert!(v["std_drop_at_target_m"].as_f64().unwrap() > 0.0);
    assert!(v["std_wind_drift_at_target_m"].as_f64().unwrap() > 0.0);

    // No key beyond this exact set (a stray field would be a silent wire change).
    let keys: std::collections::BTreeSet<&str> =
        v.as_object().unwrap().keys().map(String::as_str).collect();
    let expected: std::collections::BTreeSet<&str> = [
        "schema_version", "method", "assumptions", "confidence_percent", "hit_probability",
        "ci_low", "ci_high", "samples", "attempts", "arrivals", "stop_reason", "hit_radius_m",
        "target_distance_m", "mean_impact_velocity_mps", "std_impact_velocity_mps",
        "mean_drop_at_target_m", "std_drop_at_target_m", "mean_wind_drift_at_target_m",
        "std_wind_drift_at_target_m",
    ]
    .into_iter()
    .collect();
    assert_eq!(keys, expected, "AdaptiveMcReportV1's JSON key set changed underneath this pin");
}

// ---- (c) a very loose --target-ci-half-width is met immediately, at min-samples. ----

#[test]
fn adaptive_loose_half_width_stops_at_min_samples() {
    let (stdout, stderr, ok) = run(&[
        "--target-distance", "300",
        "--adaptive",
        "--seed", "5",
        "--min-samples", "20",
        "--mc-batch-size", "20",
        "--max-samples", "5000",
        "--target-ci-half-width", "0.5",
        "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["stop_reason"], "target_half_width_met");
    assert_eq!(v["samples"], 20, "must stop at exactly min_samples on this easy a target");
    assert_eq!(v["attempts"], 20);
}

// ---- (d) --confidence 85 is a usage error naming --confidence and listing the three choices.

#[test]
fn bad_confidence_value_is_a_usage_error_naming_the_flag_and_choices() {
    let (_stdout, stderr, ok) = run(&["--confidence", "85"]);
    assert!(!ok, "85 is not one of the three permitted confidence levels");
    assert!(stderr.contains("--confidence"), "{stderr}");
    assert!(stderr.contains("90"), "{stderr}");
    assert!(stderr.contains("95"), "{stderr}");
    assert!(stderr.contains("99"), "{stderr}");
}

/// `--target-ci-half-width`'s own positivity check must likewise name the flag (the
/// dial-plan/adaptive-card budget-flags precedent), including for a negative value, which
/// requires `allow_hyphen_values` to even reach the validator instead of clap's unrelated
/// "unexpected argument" message.
#[test]
fn non_positive_half_width_is_a_usage_error_naming_the_flag() {
    for bad in ["0", "-1", "-0.5"] {
        let (_stdout, stderr, ok) = run(&["--adaptive", "--target-ci-half-width", bad]);
        assert!(!ok, "{bad} must be rejected");
        assert!(stderr.contains("--target-ci-half-width"), "value {bad}: {stderr}");
    }
}

// ---- (e) an unreachable half-width honestly caps at --max-samples; exit 0 (Plan B's
// budget-not-met precedent: an honest cap is a successful run, not an error). 1500 is the
// brief's own number for this case; at ~1ms/trial in this fixture that is roughly a second. ----

#[test]
fn adaptive_impossible_half_width_caps_at_max_samples_and_still_exits_zero() {
    let (stdout, stderr, ok) = run(&[
        "--target-distance", "300",
        "--adaptive",
        "--seed", "7",
        "--min-samples", "20",
        "--mc-batch-size", "100",
        "--max-samples", "1500",
        "--target-ci-half-width", "0.000001",
        "-o", "json",
    ]);
    assert!(ok, "an honest max-samples cap is still exit 0; stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["stop_reason"], "max_samples_reached");
    assert_eq!(v["samples"], 1500);
    assert_eq!(v["attempts"], 1500);
}

// ---- supporting coverage: the --adaptive-only knobs require --adaptive, naming it. ----

#[test]
fn adaptive_only_flags_require_adaptive_and_name_it() {
    for flag_value in [
        ["--min-samples", "50"],
        ["--max-samples", "50000"],
        ["--mc-batch-size", "250"],
    ] {
        let (_stdout, stderr, ok) = run(&flag_value);
        assert!(!ok, "{flag_value:?} without --adaptive must be a usage error");
        assert!(stderr.contains("--adaptive"), "{flag_value:?}: {stderr}");
    }
    // --target-ci-half-width too, checked separately since it takes a value already used above.
    let (_stdout, stderr, ok) = run(&["--target-ci-half-width", "0.1"]);
    assert!(!ok);
    assert!(stderr.contains("--adaptive"), "{stderr}");
}

// ---- supporting coverage: --adaptive and --wez are mutually exclusive (the WEZ sweep stays
// fixed-count, spec 14). ----

#[test]
fn adaptive_and_wez_conflict() {
    let (_stdout, stderr, ok) = run(&["--adaptive", "--wez", "--target-size", "18x30"]);
    assert!(!ok, "--adaptive and --wez must conflict");
    assert!(stderr.contains("--adaptive"), "{stderr}");
    assert!(stderr.contains("--wez"), "{stderr}");
}

// ---- supporting coverage: --adaptive has no per-trial data to tabulate, so -o statistics is a
// clean, named usage error rather than a silent fallback or a fabricated CSV row. ----

#[test]
fn adaptive_statistics_output_is_a_clean_named_error() {
    let (_stdout, stderr, ok) = run(&[
        "--target-distance", "300",
        "--adaptive",
        "--seed", "1",
        "--min-samples", "20",
        "--mc-batch-size", "20",
        "--max-samples", "40",
        "-o", "statistics",
    ]);
    assert!(!ok, "--output statistics must be rejected under --adaptive");
    assert!(stderr.contains("--output statistics"), "{stderr}");
    assert!(stderr.contains("--adaptive"), "{stderr}");
}

// ---- supporting coverage: --adaptive's table output is a plain summary block, not the
// fixed-count box, and states the same facts -o json does (method, level, p_hit + CI, the
// three cardinalities, stop reason, geometry, and the six Welford stats). ----

#[test]
fn adaptive_table_output_is_a_summary_block_with_every_reported_fact() {
    let (stdout, stderr, ok) = run(&[
        "--target-distance", "300",
        "--adaptive",
        "--seed", "5",
        "--min-samples", "20",
        "--mc-batch-size", "20",
        "--max-samples", "5000",
        "--target-ci-half-width", "0.5",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(stdout.contains("Adaptive Monte Carlo"), "{stdout}");
    assert!(stdout.contains("anytime_beta_binomial_mixture_cs_v1"), "{stdout}");
    assert!(stdout.contains("95%"), "{stdout}");
    assert!(stdout.contains("hit probability:"), "{stdout}");
    assert!(stdout.contains("CI: ["), "{stdout}");
    assert!(stdout.contains("samples:"), "{stdout}");
    assert!(stdout.contains("attempts:"), "{stdout}");
    assert!(stdout.contains("arrivals:"), "{stdout}");
    assert!(stdout.contains("stop reason:"), "{stdout}");
    assert!(stdout.contains("target_half_width_met"), "{stdout}");
    assert!(stdout.contains("hit radius:"), "{stdout}");
    assert!(stdout.contains("target distance:"), "{stdout}");
    assert!(stdout.contains("impact velocity:"), "{stdout}");
    assert!(stdout.contains("drop at target:"), "{stdout}");
    assert!(stdout.contains("drift at target:"), "{stdout}");
    // Not the fixed-count ASCII box.
    assert!(!stdout.contains('║'), "{stdout}");
}

// ---- supporting coverage: `-o json` and `-o full` are the same value under --adaptive. ----

#[test]
fn adaptive_o_json_and_o_full_are_the_same_report() {
    let common = [
        "--target-distance", "300",
        "--adaptive",
        "--seed", "5",
        "--min-samples", "20",
        "--mc-batch-size", "20",
        "--max-samples", "5000",
        "--target-ci-half-width", "0.5",
    ];

    let (json_out, stderr, ok) = {
        let mut t = common.to_vec();
        t.extend(["-o", "json"]);
        run(&t)
    };
    assert!(ok, "stderr: {stderr}");

    let (full_out, stderr, ok) = {
        let mut t = common.to_vec();
        t.extend(["-o", "full"]);
        run(&t)
    };
    assert!(ok, "stderr: {stderr}");

    assert_eq!(json_out, full_out, "-o json must be an alias for -o full, not a distinct shape");
}

// ---- review fix round (task-5-review.md), I1: `json` must be discoverable in --help's own
// possible-values list, and the variant's doc comment must read as user-facing help, not a
// maintainer note naming private functions/structs. ----

#[test]
fn output_help_lists_json_and_carries_no_maintainer_prose() {
    let output = Command::new(bin())
        .args(["monte-carlo", "--help"])
        .output()
        .expect("run ballistics --help");
    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);

    assert!(help.contains("- json:"), "json must be a listed possible value: {help}");
    // `full` still works (checked functionally by adaptive_o_json_and_o_full_are_the_same_report
    // above) but is a backward-compatible alias now, not the primary/listed spelling.
    assert!(!help.contains("- full:"), "{help}");

    // Scope the maintainer-note-leak check to just the `-o, --output` flag's own block: several
    // OTHER flags on this subcommand (e.g. --adaptive) legitimately cite MBA-1352 in their own
    // user-facing help, so checking the whole page would false-positive on those.
    let output_block_start = help.find("-o, --output").expect("--output flag block");
    let output_block = &help[output_block_start..];
    let output_block_end = output_block.find("-h, --help").unwrap_or(output_block.len());
    let output_block = &output_block[..output_block_end];

    // The old doc comment named these private items and a ticket ID as CLI help text, right on
    // the `full`/`json` value's own line; none of that is user-facing information and must not
    // reappear in this flag's help.
    for leaked in [
        "run_monte_carlo_adaptive",
        "sample_one_trial",
        "MonteCarloTrialSampler",
        "dispatch match",
        "MBA-1352",
    ] {
        assert!(!output_block.contains(leaked), "maintainer-note leak {leaked:?}: {output_block}");
    }
}

// ---- review fix round, I2: --confidence/--seed are no-ops under --wez (spec 14 scopes WEZ out
// of sampling statistics); both must say so on stderr rather than silently doing nothing, and
// the run must still succeed. Fast fixture: a single wez step, tiny -n. ----

#[test]
fn wez_with_confidence_or_seed_notes_the_no_op_and_still_succeeds() {
    let wez_tail = [
        "-n", "5", "--wez", "--target-size", "30x30",
        "--wez-start", "100", "--wez-end", "100", "--wez-step", "100",
    ];

    let mut with_both = wez_tail.to_vec();
    with_both.extend(["--confidence", "90", "--seed", "7"]);
    let (stdout, stderr, ok) = run(&with_both);
    assert!(ok, "a disclosed no-op must still succeed; stderr: {stderr}");
    assert!(
        stderr.contains("--confidence") && stderr.contains("--wez") && stderr.contains("no effect"),
        "{stderr}"
    );
    assert!(
        stderr.contains("--seed") && stderr.contains("no effect"),
        "{stderr}"
    );
    assert!(stdout.contains("WEZ sweep:"), "the sweep itself must still run: {stdout}");

    // Neither flag given: no note at all.
    let (_stdout, stderr, ok) = run(&wez_tail);
    assert!(ok, "stderr: {stderr}");
    assert!(!stderr.contains("no effect"), "{stderr}");
}

// ---- review fix round, I3: the max_samples < min_samples cross-check must name the CLI flags
// (--max-samples/--min-samples), never the library's McConvergence field names, and must say
// when one side of the comparison is an untouched default rather than something the user typed.

#[test]
fn cross_field_sample_bounds_error_names_the_flags_and_states_defaults() {
    // Both explicit: plain flag-named values, no struct name anywhere.
    let (_stdout, stderr, ok) = run(&["--adaptive", "--min-samples", "100", "--max-samples", "50"]);
    assert!(!ok, "50 < 100 must be rejected");
    assert!(stderr.contains("--max-samples (50)"), "{stderr}");
    assert!(stderr.contains("--min-samples (100)"), "{stderr}");
    assert!(!stderr.contains("McConvergence"), "the library struct name must not leak: {stderr}");

    // Only --min-samples given: the untouched --max-samples default must say so.
    let (_stdout, stderr, ok) = run(&["--adaptive", "--min-samples", "200000"]);
    assert!(!ok);
    assert!(stderr.contains("--max-samples (default 100000)"), "{stderr}");
    assert!(stderr.contains("--min-samples (200000)"), "{stderr}");

    // Only --max-samples given: the untouched --min-samples default must say so.
    let (_stdout, stderr, ok) = run(&["--adaptive", "--max-samples", "500"]);
    assert!(!ok);
    assert!(stderr.contains("--max-samples (500)"), "{stderr}");
    assert!(stderr.contains("--min-samples (default 1000)"), "{stderr}");

    // One side left at its default (here --max-samples's 100000, well above a tiny explicit
    // --min-samples) must NOT error -- a sanity check that resolving a default doesn't itself
    // trip the new pre-check. `--min-samples 20`/loose half-width keeps this fast (~20 trials)
    // rather than running out the full 1000-sample default floor this fix must not require.
    let (_stdout, stderr, ok) = run(&[
        "--target-distance", "300", "--adaptive", "--seed", "1",
        "--min-samples", "20", "--mc-batch-size", "20", "--target-ci-half-width", "0.9",
    ]);
    assert!(ok, "one defaulted side must remain valid; stderr: {stderr}");
}
