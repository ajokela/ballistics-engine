//! Integration tests for the bridge's `true.*` truing commands: `true.fit`, `true.wind`,
//! `true.tall_target`, `true.dsf`, `true.plan`.
//!
//! The `_matches_cli` tests at the bottom are golden cross-checks (same pattern as
//! `tests/card_bridge_golden.rs`): they run the real `ballistics` binary and the bridge
//! over the same inputs and assert the numbers agree. The CLI is the reference
//! implementation; without these, a future change to either `main.rs` or
//! `truing_service.rs` could silently diverge from the other with nothing to catch it —
//! exactly how this branch's three CLI-vs-library field-loss regressions were only ever
//! caught by manual worktree diffing. Requires both the `bridge` and `cli` features (the
//! default set) because they execute the actual binary via `CARGO_BIN_EXE`.

#![cfg(all(feature = "bridge", feature = "cli"))]

use ballistics_engine::bridge::bridge_call;
use serde_json::{json, Value};
use std::process::Command;

fn call(command: &str, request: Value) -> Value {
    let envelope = json!({"api_version": 1, "command": command, "request": request});
    serde_json::from_str(&bridge_call(&envelope.to_string())).unwrap()
}

/// Spawn the real `ballistics` binary.
fn cli() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ballistics"))
}

const PARITY_TOL: f64 = 1e-9;

/// Same shape as `card_bridge_golden.rs`'s `assert_close`, but over two bare `f64`s
/// (one of the two truing CLI outputs being compared is parsed from plain-text stdout,
/// not JSON, so there is no pair of `serde_json::Value`s to hand it).
fn assert_close_f64(label: &str, cli_value: f64, bridge_value: f64) {
    let scale = cli_value.abs().max(bridge_value.abs()).max(1.0);
    assert!(
        (cli_value - bridge_value).abs() <= PARITY_TOL * scale,
        "{label}: CLI {cli_value} vs bridge {bridge_value}"
    );
}

fn model() -> Value {
    json!({
        "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
        "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
        "zero_distance_yd": 100.0, "sight_height_in": 2.0,
        "temperature_f": 59.0, "pressure_inhg": 29.92,
        "humidity_pct": 50.0, "altitude_ft": 0.0
    })
}

/// Same bullet as `model()`, but a 900 yd zero rather than 100 yd — DELIBERATELY NOT
/// `model()` itself (`true.fit`'s tests depend on that one staying at a 100 yd zero).
///
/// The `dsf` solve path runs with ground-impact detection unconditionally on, terminating
/// the trajectory once it has dropped ~60 in below the muzzle line (a fixed historical
/// default — saved profiles predate the CLI's unified `--bore-height` flag and never
/// stored one). With `model()`'s 100 yd zero, the bore is barely angled up, so the
/// trajectory sinks through that ground plane at ~516 yd — while still Mach 1.6+, well
/// short of the Mach <= 1.2 transonic band `dsf` requires. Swept every reachable range for
/// that fixture and confirmed there is no way to land in-band with it (see the task-10
/// report). A long zero angles the bore up substantially, keeping the projectile above the
/// ground plane far enough downrange to actually go transonic before impact; 900 yd was
/// chosen empirically (swept 700/800/900/1000/1100 yd zeros) as the one that lands
/// squarely mid-band rather than hugging an edge. Do not "simplify" this back to `model()`.
fn dsf_model() -> Value {
    json!({
        "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
        "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
        "zero_distance_yd": 900.0, "sight_height_in": 2.0,
        "temperature_f": 59.0, "pressure_inhg": 29.92,
        "humidity_pct": 50.0, "altitude_ft": 0.0
    })
}

#[test]
fn true_fit_is_advertised_in_capabilities() {
    let v = call("meta.capabilities", json!(null));
    let names: Vec<&str> = v["result"]["commands"]
        .as_array().unwrap().iter().map(|c| c.as_str().unwrap()).collect();
    assert!(names.contains(&"true.fit"), "capabilities were {names:?}");
}

#[test]
fn true_fit_returns_estimates_with_their_uncertainty() {
    // Observations derived from the engine's own solve on the model below:
    // .308 caliber, 168 gr, 0.308 in dia, G7 BC 0.243, 2700 fps, 100 yd zero, ICAO conditions.
    // Drops (mil) at range: 400 yd -> 2.6882, 600 yd -> 4.8626, 800 yd -> 7.6598
    let v = call("true.fit", json!({
        "model": model(),
        "drop_unit": "mil",
        "observations": [
            {"range_yd": 400.0, "drop": 2.6882, "sigma": 0.1},
            {"range_yd": 600.0, "drop": 4.8626, "sigma": 0.1},
            {"range_yd": 800.0, "drop": 7.6598, "sigma": 0.1}
        ],
        "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 20.0},
                   "ballistic_coefficient": {"mean": 0.243, "sigma": 0.01}},
        "predictions": []
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    let r = &v["result"];
    assert!(r["map_muzzle_velocity_fps"].as_f64().unwrap() > 0.0);
    assert!(r["map_ballistic_coefficient"].as_f64().unwrap() > 0.0);

    // The honesty invariant: a point estimate never travels without its uncertainty.
    // Assert that the approximation is available and carries valid intervals for both parameters.
    assert_eq!(r["approximation"]["status"], "available",
        "approximation must be available, got: {}", r["approximation"]["status"]);
    let details = &r["approximation"]["details"];

    // Verify that both parameters have finite intervals.
    let mv_lower = details["muzzle_velocity_interval_95"]["lower"].as_f64().unwrap();
    let mv_upper = details["muzzle_velocity_interval_95"]["upper"].as_f64().unwrap();
    assert!(mv_lower.is_finite() && mv_upper.is_finite(),
        "muzzle velocity interval bounds must be finite");
    assert!(mv_lower < mv_upper, "muzzle velocity lower bound must be less than upper");

    let bc_lower = details["ballistic_coefficient_interval_95"]["lower"].as_f64().unwrap();
    let bc_upper = details["ballistic_coefficient_interval_95"]["upper"].as_f64().unwrap();
    assert!(bc_lower.is_finite() && bc_upper.is_finite(),
        "ballistic coefficient interval bounds must be finite");
    assert!(bc_lower < bc_upper, "ballistic coefficient lower bound must be less than upper");
}

#[test]
fn true_fit_rejects_a_null_payload_and_unknown_fields() {
    let v = call("true.fit", json!(null));
    assert_eq!(v["error"]["code"], "invalid_request");

    let v = call("true.fit", json!({
        "model": model(), "drop_unit": "mil",
        "observations": [{"range_yd": 600.0, "drop": 4.3, "sigma": 0.1, "typo": 1}],
        "priors": {"muzzle_velocity_fps": null, "ballistic_coefficient": null},
        "predictions": []
    }));
    assert_eq!(v["error"]["code"], "invalid_request");
}

#[test]
fn true_wind_fits_a_crosswind_from_an_observed_miss() {
    let v = call("true.wind", json!({
        "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
        "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
        "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
        "sight_height_in": 2.0, "temperature_f": 59.0, "pressure_inhg": 29.92,
        "humidity_pct": 50.0, "altitude_ft": 0.0,
        "twist": {"rate_in": 11.0, "right_hand": true},
        "earth": null, "called_crosswind_mph": null
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    assert_eq!(v["result"]["solutions"].as_array().unwrap().len(), 1);

    // Physics band, not a bare "is a number": sensitivity_m_per_mph ~= 0.05263 m/mph,
    // no_wind_lateral_m (spin drift, 1:11 twist) ~= 0.0533 m, observed miss = 0.315 m,
    // so wind-attributable miss = 0.315 - 0.0533 = 0.2617 m, / 0.05263 ~= 4.97 mph.
    // A wide band (not an exact float) so this isn't brittle against future solver tuning,
    // but tight enough that a sign flip, unit error, or 10x scale error fails it.
    let mean = v["result"]["mean_crosswind_mph"].as_f64().unwrap();
    assert!(mean > 4.5 && mean < 5.5, "mean_crosswind_mph out of band: {mean}");
    // Sign check: a rightward miss must solve to a positive crosswind (from the left,
    // pushing right) -- catches a sign flip, which a bare is_number() cannot.
    assert!(mean > 0.0, "expected a positive (from-the-left) crosswind, got {mean}");

    // Convergence check: the fit must actually reproduce the observation, within the
    // solver's own WIND_SOLVE_TOLERANCE_M (1.0e-5 m), loosened slightly for headroom.
    let solution = &v["result"]["solutions"][0];
    let modeled = solution["modeled_miss_right_m"].as_f64().unwrap();
    assert!(
        (modeled - 0.315).abs() < 1.0e-4,
        "modeled_miss_right_m did not reproduce the observed 0.315 m miss: {modeled}"
    );
    let residual = solution["residual_m"].as_f64().unwrap();
    assert!(residual.abs() < 1.0e-4, "residual_m too large: {residual}");
}

#[test]
fn true_wind_is_advertised_in_capabilities() {
    let v = call("meta.capabilities", json!(null));
    let names: Vec<&str> = v["result"]["commands"]
        .as_array().unwrap().iter().map(|c| c.as_str().unwrap()).collect();
    assert!(names.contains(&"true.wind"), "capabilities were {names:?}");
}

#[test]
fn true_tall_target_computes_a_correction_factor() {
    let v = call("true.tall_target", json!({
        "dialed": 10.0, "measured": 30.0, "range": 100.0,
        "unit": "mil", "metric": false
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    // 30 in at 100 yd is 30/36 yd, i.e. 8.333333333333334 mil of ACTUAL travel against 10 mil
    // dialed, so the correction factor is 0.8333333333333334. This is exact rational
    // arithmetic (30/36 * 10), so assert it tightly rather than with slack.
    let actual = v["result"]["actual"].as_f64().unwrap();
    let cf = v["result"]["correction_factor"].as_f64().unwrap();
    assert!((actual - 8.333333333333334).abs() < 1e-9, "actual travel was {actual}");
    assert!((cf - 0.8333333333333334).abs() < 1e-9, "cf was {cf}");
    assert_eq!(v["result"]["within_accepted_band"], true);
}

#[test]
fn true_tall_target_rejects_clicks() {
    let v = call("true.tall_target", json!({
        "dialed": 10.0, "measured": 30.0, "range": 100.0,
        "unit": "clicks", "metric": false
    }));
    assert_eq!(v["ok"], false);
    assert_eq!(v["error"]["code"], "command_failed");
    assert!(
        v["error"]["message"].as_str().unwrap().contains("clicks"),
        "expected the clicks guard, got: {}", v["error"]["message"]
    );
}

#[test]
fn true_dsf_derives_a_point_without_touching_a_profile() {
    // 900 yd zero, 950 yd observation: predicted drop is 0.9468 mil there (read off a
    // first run, per this fixture's own doc comment on `dsf_model`), so 1.0 mil observed
    // sits close to it without being an exact, untestable match — dsf lands at ~1.056,
    // comfortably inside the sane band below rather than hugging either edge of it.
    let v = call("true.dsf", json!({
        "model": dsf_model(), "range_yd": 950.0,
        "observed_drop": 1.0, "drop_unit": "mil"
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    let mach = v["result"]["mach"].as_f64().unwrap();
    let dsf = v["result"]["dsf"].as_f64().unwrap();
    // The observation must land in the transonic band DSF exists for: at or below the
    // DSF_MACH_CEILING of 1.2, and still moving. `> 0.0` would pass on a supersonic or
    // nonsense Mach, which is the case the service is supposed to REFUSE.
    assert!(mach > 0.5 && mach <= 1.2, "mach was {mach}, outside the DSF band");
    // dsf is observed/predicted drop. A sane correction is near unity; 0.5..2.0 still
    // catches a unit error or an inverted ratio, which `> 0.0` would not.
    assert!(dsf > 0.5 && dsf < 2.0, "dsf was {dsf}");
    assert!(v["result"]["warnings"].is_array());
}

#[test]
fn all_five_true_commands_are_advertised() {
    let v = call("meta.capabilities", json!(null));
    let names: Vec<&str> = v["result"]["commands"]
        .as_array().unwrap().iter().map(|c| c.as_str().unwrap()).collect();
    for c in ["true.fit", "true.wind", "true.tall_target", "true.dsf", "true.plan"] {
        assert!(names.contains(&c), "{c} missing from {names:?}");
    }
}

// ============================================================================
// true.plan
// ============================================================================

/// Same bullet as `model()`. Five candidate ranges, 3-station design, 0.1 mil measurement
/// resolution: the 600/800/1000 yd trio is the informative one for this load (farther
/// stations separate MV and BC better than the two short ones), so a regression that
/// mangles the sensitivity math or the station selection changes which three ranges come
/// back, not just the numbers attached to them. Values below were read off a real run of
/// this exact fixture (see `dsf_model`'s neighboring doc comments for why this crate reads
/// real solver output rather than hand-deriving it).
#[test]
fn true_plan_selects_the_farthest_stations_for_a_joint_mv_bc_design() {
    let v = call("true.plan", json!({
        "model": model(),
        "candidate_ranges_yd": [200.0, 400.0, 600.0, 800.0, 1000.0],
        "observation_count": 3,
        "minimum_separation_yd": 100.0,
        "measurement_sigma_1sd": 0.1,
        "drop_unit": "mil"
    }));
    assert_eq!(v["ok"], true, "response was {v}");
    let r = &v["result"];

    assert_eq!(r["mode"], "joint_mv_bc", "{r}");
    assert_eq!(r["search_strategy"], "exhaustive", "{r}");
    // All 5 candidates are unique and reachable at this load -- none rejected.
    assert_eq!(r["eligible_candidate_count"], 5, "{r}");
    assert!(r["rejected_candidates"].as_array().unwrap().is_empty(), "{r}");

    // The optimizer picks the three FARTHEST candidates, not the three closest or a
    // scattered set -- this is the substantive claim a broken sensitivity calculation
    // would falsify first.
    let selected_ranges: Vec<f64> = r["selected_stations"]
        .as_array().unwrap().iter()
        .map(|s| s["range_yd"].as_f64().unwrap())
        .collect();
    assert_eq!(selected_ranges, vec![600.0, 800.0, 1000.0], "{r}");
    let unselected: Vec<f64> = r["unselected_candidate_ranges_yd"]
        .as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    assert_eq!(unselected, vec![200.0, 400.0], "{r}");

    // BC sensitivity ratio and the weak-axis singular value are real physics-derived
    // numbers, not merely "some positive float" -- a wide-but-real band that would catch a
    // sign flip or an order-of-magnitude regression while surviving ordinary solver tuning.
    let sensitivity_ratio = r["information"]["sensitivity_ratio"].as_f64().unwrap();
    assert!(
        sensitivity_ratio > 0.25 && sensitivity_ratio < 0.4,
        "sensitivity_ratio out of band: {sensitivity_ratio}"
    );
    let min_singular = r["information"]["minimum_singular_value"].as_f64().unwrap();
    assert!(
        min_singular > 10.0 && min_singular < 25.0,
        "minimum_singular_value out of band: {min_singular}"
    );
}

#[test]
fn true_plan_rejects_unknown_fields_on_an_otherwise_valid_request() {
    let mut req = json!({
        "model": model(),
        "candidate_ranges_yd": [200.0, 400.0, 600.0],
        "observation_count": 2,
        "minimum_separation_yd": 100.0,
        "measurement_sigma_1sd": 0.1,
        "drop_unit": "mil"
    });
    // Confirm the request is valid BEFORE poisoning it, so rejection below can only be
    // attributed to the unknown field.
    assert_eq!(call("true.plan", req.clone())["ok"], true);
    req["typo_field"] = json!(1.0);
    let v = call("true.plan", req);
    assert_eq!(v["error"]["code"], "invalid_request", "{v}");
}

#[test]
fn true_plan_reports_insufficient_candidates_with_rejection_diagnostics() {
    // 200 yd duplicated (rejected as duplicate_range) and -50 yd invalid (rejected as
    // invalid_range) leave only ONE unique reachable candidate for a 3-station request.
    let v = call("true.plan", json!({
        "model": model(),
        "candidate_ranges_yd": [200.0, 200.0, -50.0],
        "observation_count": 3,
        "minimum_separation_yd": 100.0,
        "measurement_sigma_1sd": 0.1,
        "drop_unit": "mil"
    }));
    assert_eq!(v["ok"], false, "{v}");
    assert_eq!(v["error"]["code"], "command_failed", "{v}");
    assert_eq!(
        v["error"]["details"]["code"], "insufficient_reachable_candidates",
        "{v}"
    );
    let rejected = v["error"]["details"]["rejected_candidates"].as_array().unwrap();
    let by_index = |index: i64| {
        rejected.iter().find(|c| c["input_index"] == index).unwrap_or_else(|| {
            panic!("no rejected candidate at input_index {index} in {rejected:?}")
        })
    };
    assert_eq!(by_index(1)["reason"], "duplicate_range", "{rejected:?}");
    assert_eq!(by_index(2)["reason"], "invalid_range", "{rejected:?}");
    assert_eq!(by_index(2)["range_yd"], -50.0, "{rejected:?}");
}

// ============================================================================
// CLI-vs-bridge parity (whole-branch review, Finding 2)
// ============================================================================

/// `tall-target` has no `-o json`; the CLI's only output is the plain-text report
/// (`Commands::TallTarget`'s handler, `src/main.rs`). Pull the "Actual travel" and
/// "Correction factor" numbers back out of it so they can be compared against the
/// bridge's JSON result.
fn parse_tall_target_stdout(stdout: &str) -> (f64, f64) {
    let mut actual = None;
    let mut correction_factor = None;
    for line in stdout.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("Actual travel:") {
            actual = rest.split_whitespace().next().and_then(|s| s.parse::<f64>().ok());
        } else if let Some(rest) = line.strip_prefix("Correction factor (actual / dialed):") {
            correction_factor = rest.split_whitespace().next().and_then(|s| s.parse::<f64>().ok());
        }
    }
    (
        actual.expect("CLI printed an 'Actual travel' line"),
        correction_factor.expect("CLI printed a 'Correction factor' line"),
    )
}

#[test]
fn true_tall_target_matches_cli() {
    // 45 in at 100 yd against a 10 mil dial: actual = (45/36/100)*1000 = 12.5 mil exactly,
    // cf = 1.25 exactly (both terminate at <= 2 decimal places), so the CLI's rounded
    // 2/4-decimal display loses no information a tight tolerance would need to absorb —
    // this really is the exact-arithmetic check the design spec calls for, not merely a
    // close one.
    let out = cli()
        .args([
            "tall-target", "--dialed", "10", "--measured", "45", "--range", "100",
            "--unit", "mil",
        ])
        .output()
        .expect("run tall-target");
    assert!(
        out.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (cli_actual, cli_cf) = parse_tall_target_stdout(&String::from_utf8_lossy(&out.stdout));

    let v = call(
        "true.tall_target",
        json!({"dialed": 10.0, "measured": 45.0, "range": 100.0, "unit": "mil", "metric": false}),
    );
    assert_eq!(v["ok"], true, "bridge true.tall_target failed: {v}");
    let bridge_actual = v["result"]["actual"].as_f64().unwrap();
    let bridge_cf = v["result"]["correction_factor"].as_f64().unwrap();

    assert_close_f64("actual", cli_actual, bridge_actual);
    assert_close_f64("correction_factor", cli_cf, bridge_cf);
}

/// Unique-per-call `$HOME` so `profile save`/`dsf` (which persist to
/// `$HOME/.ballistics/profiles`) never collide with another test running concurrently —
/// same isolation pattern `tests/dsf_workflow.rs` uses for the same reason.
fn temp_home_dir(tag: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU32, Ordering};
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "truing-bridge-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn true_dsf_matches_cli() {
    // `dsf --saved-profile` has no bare-model equivalent in the bridge (its whole point is
    // the profile-shaped extras `DsfSolveInputs` carries — see `truing_service.rs`'s module
    // doc comment), so the only way to drive it through the CLI is a real saved profile
    // under an isolated $HOME, same as `tests/dsf_workflow.rs`. This profile mirrors
    // `dsf_model()` above field-for-field (same 900 yd zero, for the same reason
    // `dsf_model()`'s own doc comment gives: it is what lands the observation in-band).
    let home = temp_home_dir("dsf");
    let name = "parity";

    let save_out = cli()
        .env("HOME", &home)
        .args([
            "profile", "save", name,
            "-v", "2700", "-b", "0.243", "-m", "168", "-d", "0.308",
            "--drag-model", "g7",
            "--zero-distance", "900",
            "--sight-height", "2.0",
            "--temperature", "59",
            "--pressure", "29.92",
            "--humidity", "50",
            "--altitude", "0",
        ])
        .output()
        .expect("spawn profile save");
    assert!(
        save_out.status.success(),
        "profile save failed: {}",
        String::from_utf8_lossy(&save_out.stderr)
    );

    let dsf_out = cli()
        .env("HOME", &home)
        .args(["dsf", "--saved-profile", name, "--range", "950"])
        .arg("--observed-drop=1.0mil")
        .output()
        .expect("spawn dsf");
    assert!(
        dsf_out.status.success(),
        "dsf failed: {}",
        String::from_utf8_lossy(&dsf_out.stderr)
    );

    let profile_path = home.join(".ballistics").join("profiles").join(format!("{name}.json"));
    let profile: Value =
        serde_json::from_str(&std::fs::read_to_string(&profile_path).expect("read saved profile"))
            .expect("saved profile is valid JSON");
    let points = profile["dsf_points"].as_array().expect("saved profile has dsf_points");
    let point = points.last().expect("dsf wrote at least one point");
    let cli_mach = point["mach"].as_f64().unwrap();
    let cli_dsf = point["dsf"].as_f64().unwrap();

    let v = call(
        "true.dsf",
        json!({"model": dsf_model(), "range_yd": 950.0, "observed_drop": 1.0, "drop_unit": "mil"}),
    );
    assert_eq!(v["ok"], true, "bridge true.dsf failed: {v}");
    let bridge_mach = v["result"]["mach"].as_f64().unwrap();
    let bridge_dsf = v["result"]["dsf"].as_f64().unwrap();

    assert_close_f64("mach", cli_mach, bridge_mach);
    assert_close_f64("dsf", cli_dsf, bridge_dsf);

    let _ = std::fs::remove_dir_all(&home);
}

#[test]
fn true_plan_matches_cli() {
    // Same load and design parameters as `true_plan_selects_the_farthest_stations...`.
    // `plan-truing -o json` restructures the payload (unit-suffixed keys, a nested
    // `information`/`optimizer` split -- see `display_truing_experiment_plan` in
    // `main.rs`) rather than serializing `TruingExperimentPlanV1` verbatim, so this
    // compares the underlying VALUES field-by-field instead of the two JSON documents
    // wholesale.
    let out = cli()
        .args([
            "plan-truing",
            "-v", "2700", "-b", "0.243", "--drag-model", "g7", "-m", "168", "-d", "0.308",
            "--candidate-ranges", "200,400,600,800,1000",
            "--observation-count", "3",
            "--minimum-separation", "100",
            "--measurement-resolution", "0.1",
            "--drop-unit", "mil",
            "--zero-distance", "100",
            "--sight-height", "2.0",
            "--temperature", "59",
            "--pressure", "29.92",
            "--humidity", "50",
            "--altitude", "0",
            "-o", "json",
        ])
        .output()
        .expect("run plan-truing");
    assert!(
        out.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let cli_json: Value =
        serde_json::from_slice(&out.stdout).expect("plan-truing -o json printed valid JSON");

    let v = call("true.plan", json!({
        "model": model(),
        "candidate_ranges_yd": [200.0, 400.0, 600.0, 800.0, 1000.0],
        "observation_count": 3,
        "minimum_separation_yd": 100.0,
        "measurement_sigma_1sd": 0.1,
        "drop_unit": "mil"
    }));
    assert_eq!(v["ok"], true, "bridge true.plan failed: {v}");
    let r = &v["result"];

    assert_eq!(cli_json["mode"], r["mode"], "mode: cli {cli_json} vs bridge {r}");
    assert_eq!(
        cli_json["requested_observation_count"], r["requested_observation_count"],
        "requested_observation_count"
    );
    assert_eq!(
        cli_json["optimizer"]["eligible_candidate_count"], r["eligible_candidate_count"],
        "eligible_candidate_count"
    );
    assert_eq!(
        cli_json["optimizer"]["search_strategy"], r["search_strategy"],
        "search_strategy"
    );
    assert_close_f64(
        "sensitivity_ratio",
        cli_json["information"]["bc_sensitivity_ratio"].as_f64().unwrap(),
        r["information"]["sensitivity_ratio"].as_f64().unwrap(),
    );
    assert_close_f64(
        "minimum_singular_value",
        cli_json["information"]["minimum_singular_value"].as_f64().unwrap(),
        r["information"]["minimum_singular_value"].as_f64().unwrap(),
    );

    let cli_stations = cli_json["selected_stations"].as_array().unwrap();
    let bridge_stations = r["selected_stations"].as_array().unwrap();
    assert_eq!(cli_stations.len(), bridge_stations.len(), "selected station count");
    for (cli_station, bridge_station) in cli_stations.iter().zip(bridge_stations) {
        assert_eq!(
            cli_station["input_index"], bridge_station["input_index"],
            "input_index: {cli_station} vs {bridge_station}"
        );
        assert_close_f64(
            "range_yd",
            cli_station["range_yd"].as_f64().unwrap(),
            bridge_station["range_yd"].as_f64().unwrap(),
        );
        assert_close_f64(
            "predicted_drop",
            cli_station["predicted_drop_mil"].as_f64().unwrap(),
            bridge_station["predicted_drop"].as_f64().unwrap(),
        );
        assert_close_f64(
            "scaled_mv_sensitivity",
            cli_station["scaled_mv_sensitivity"].as_f64().unwrap(),
            bridge_station["scaled_mv_sensitivity"].as_f64().unwrap(),
        );
    }
}
