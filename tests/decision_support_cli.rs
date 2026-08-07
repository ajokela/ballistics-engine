//! 0.33.0 decision-support Tasks 14-15: the `explain`, `tolerance` and `error-budget` CLI
//! surfaces over `explain_difference` (MBA-1345), `tolerance_envelope` (MBA-1350) and
//! `error_budget_with_target` (MBA-1347).
//!
//! Fixtures are deliberately small (short `max_range_m`, one or two ranges, one to two sources):
//! `explain` solves on the order of 70 times per call, `tolerance` bisects once per axis, and
//! `error-budget` costs one `central_difference` call per DECLARED source (2 real solves for an
//! ordinary axis, ~35-40 for a `requires_rezero` axis on a request that also carries
//! `zero_distance_m` -- see `src/error_budget.rs`'s own module doc, "Cost" section, for the
//! measured numbers), so this file's own runtime is dominated by fixture SIZE, not test COUNT.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU32, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_ballistics")
}

/// A fresh, uniquely-named scratch directory per call -- safe under `cargo test`'s default
/// parallel test execution (mirrors `tests/hold_corridor.rs`'s own helper).
fn tempfile_dir() -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "bx-decision-support-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_json(dir: &Path, name: &str, value: &serde_json::Value) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, value.to_string()).unwrap();
    path
}

fn run(args: &[&str]) -> (String, String, bool) {
    let output = Command::new(bin()).args(args).output().expect("run ballistics");
    (
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

/// A small solve-json v1 request. NOTE: the decoder requires the exact-case built-in spellings
/// (`"G1"`/`"G2"`/`"G5"`/`"G6"`/`"G7"`/`"G8"`/`"GI"`/`"GS"`/`"RA4"`) -- the task-14 brief's
/// own literal used lowercase `"g7"`, which `decode_solve_request_v1` rejects with
/// `InvalidValue` (the same brief-fixture mistake Tasks 7 and 12 already hit and fixed; see
/// their reports and `src/tolerance.rs`'s own test module doc comment).
fn small_request(mv: f64) -> serde_json::Value {
    serde_json::json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                       "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": mv, "sight_height_m": 0.05},
        "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
        "atmosphere": {}, "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
        "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
    })
}

/// A vacuum (drag-free), flat, unzeroed shot -- `drop(v) = 0.5 * g * (x / v)^2` exactly, and no
/// wind at all so `windage_m` is exactly `0.0` everywhere (mirrors `src/tolerance.rs`'s own
/// `vacuum_resolved` test fixture, which this file reuses to get a GUARANTEED, analytically
/// known bound instead of guessing at real-drag-model physics).
fn vacuum_request() -> serde_json::Value {
    serde_json::json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.01, "diameter_m": 0.0077, "drag_model": "G1",
                       "ballistic_coefficient": 100.0},
        "rifle": {"muzzle_velocity_mps": 800.0, "sight_height_m": 0.0},
        "shot": {"max_range_m": 500.0, "muzzle_angle_rad": 0.0},
        "atmosphere": {}, "wind": {}, "solver": {}, "effects": {},
        "sampling": {"interval_m": 5.0}
    })
}

/// A QNH-referenced atmosphere fixture, copied verbatim (field for field) from
/// `src/error_budget.rs`'s own already-verified `qnh_resolved` test fixture: `Altitude` is
/// refused there with `KernelError::AxisUnsupportedForRequest` (a STATIC guard on the resolved
/// atmosphere's own reference mode -- it costs zero real solves, never even reaching
/// `TrajectorySolver`), so declaring it alongside a real source is the cheapest possible way to
/// exercise the "unavailable, not silently dropped" path end to end through this CLI.
fn qnh_request() -> serde_json::Value {
    serde_json::json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                       "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
        "shot": {"max_range_m": 900.0},
        "atmosphere": {"altitude_m": 500.0, "temperature_k": 288.0, "pressure_pa": 101325.0,
                       "pressure_reference": "qnh"},
        "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
    })
}

// ---- Step 1 tests from the task-14 brief (a floor -- see the two adaptations noted inline) ----

/// Brief's own test, byte-identical in spirit. Adapted only for: (a) the `"g7"` -> `"G7"` decode
/// fix noted on `small_request` above, and (b) a unique temp dir per `tempfile_dir()`.
#[test]
fn explain_emits_a_versioned_report_with_a_remainder() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let b = write_json(&dir, "b.json", &small_request(870.0));
    let (stdout, stderr, ok) = run(&[
        "explain", "--a", a.to_str().unwrap(), "--b", b.to_str().unwrap(),
        "--ranges", "600", "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["method"], "symmetric_group_counterfactual");
    assert!(v["rows"][0]["interaction_remainder"].is_object());
    assert!(!v["assumptions"].as_array().unwrap().is_empty());
}

/// Brief's own test, adapted for the mandatory `--domain` this task's real library signature
/// requires (`tolerance_envelope`'s `domains: &[(InputAxis, (f64, f64))]` -- see the task-14
/// report's "divergences from the brief" section): the brief's literal invocation, run
/// unmodified, now fails fast with a clear "--axis wind-speed has no matching --domain" error
/// rather than solving anything, so `--domain wind-speed=0:20` is added here. Also extended
/// with the schema_version/assumptions checks the task instructions ask for on BOTH
/// subcommands (the brief's own tolerance test checked neither).
#[test]
fn tolerance_emits_a_versioned_report() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (stdout, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "wind-speed", "--domain", "wind-speed=0:20",
        "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["method"], "one_variable_deterministic_bisection");
    assert!(v["axes"][0]["axis"].is_string());
    assert!(!v["assumptions"].as_array().unwrap().is_empty());
}

// ---- Beyond the brief: the properties named in the task-14 instructions ----

/// Delta #1: every `--axis` needs a matching `--domain`; the CLI must say so by name (not just
/// forward the library's axis-less `KernelError::InvalidDomain`, which cannot mention the
/// `--domain` flag at all).
#[test]
fn tolerance_missing_domain_names_the_axis_and_the_flag() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "wind-speed",
    ]);
    assert!(!ok, "a missing --domain must be rejected");
    assert!(stderr.contains("wind-speed"), "{stderr}");
    assert!(stderr.contains("--domain"), "{stderr}");
}

/// Delta #4: an unrecognized `--axis` is a hard error listing the accepted names, generated
/// from `InputAxis::ALL` (the full round-trip over all 32 is pinned at the unit level --
/// `axis_kebab_names_round_trip_every_member_of_all`, `src/main.rs` -- this just confirms the
/// CLI surface actually reaches that machinery).
#[test]
fn tolerance_unknown_axis_lists_valid_names() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "not-a-real-axis",
    ]);
    assert!(!ok);
    assert!(stderr.contains("unknown axis 'not-a-real-axis'"), "{stderr}");
    assert!(stderr.contains("wind-speed"), "{stderr}");
    assert!(stderr.contains("muzzle-velocity-mps"), "{stderr}");
}

/// Review minor (b): a repeated `--axis` would otherwise double a real bisection for no reason
/// (the cost note's whole concern) -- rejected instead of silently deduplicated.
#[test]
fn tolerance_rejects_a_repeated_axis() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75",
        "--axis", "wind-speed", "--axis", "wind-speed", "--domain", "wind-speed=0:20",
    ]);
    assert!(!ok, "a repeated --axis must be rejected");
    assert!(stderr.contains("wind-speed"), "{stderr}");
    assert!(stderr.contains("more than once"), "{stderr}");
}

/// Review minor (b): a repeated `--domain` for the same axis is rejected rather than silently
/// last-wins -- a caller who fat-fingered a second `--domain` should not get an unexplained,
/// possibly different, search domain.
#[test]
fn tolerance_rejects_a_repeated_domain_for_one_axis() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "wind-speed",
        "--domain", "wind-speed=0:20", "--domain", "wind-speed=1:19",
    ]);
    assert!(!ok, "a repeated --domain for one axis must be rejected");
    assert!(stderr.contains("wind-speed"), "{stderr}");
    assert!(stderr.contains("more than once"), "{stderr}");
}

/// Review minor (b): a `--domain` for an axis never requested via `--axis` is rejected rather
/// than silently ignored -- most likely a typo'd or forgotten `--axis` on the caller's part.
#[test]
fn tolerance_rejects_a_domain_for_an_axis_not_requested() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "wind-speed",
        "--domain", "wind-speed=0:20", "--domain", "temperature=250:310",
    ]);
    assert!(!ok, "a --domain for an axis outside --axis must be rejected");
    assert!(stderr.contains("temperature"), "{stderr}");
    assert!(stderr.contains("not requested"), "{stderr}");
}

/// Delta #2: a skipped axis (here, `latitude_rad` supplied on `a` only) must show up in the
/// TABLE output, not just the JSON -- and the interaction remainder must be unmistakable for an
/// eighth group contribution. `a`/`b` are otherwise byte-identical, so this is cleanly isolated:
/// nothing else could put a "skipped" line in this table.
#[test]
fn explain_table_labels_the_remainder_and_renders_skipped_axes() {
    let dir = tempfile_dir();
    let mut a = small_request(823.0);
    a["atmosphere"] = serde_json::json!({"latitude_rad": 0.5});
    let a_path = write_json(&dir, "a.json", &a);
    let b_path = write_json(&dir, "b.json", &small_request(823.0));

    let (table, stderr, ok) = run(&[
        "explain", "--a", a_path.to_str().unwrap(), "--b", b_path.to_str().unwrap(),
        "--ranges", "600",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(
        table.contains("interaction remainder") && table.contains("unexplained by any single group"),
        "the interaction remainder must be on its own labelled line: {table}"
    );
    assert!(
        table.contains("skipped axes"),
        "a skipped axis must be rendered in the table, not dropped: {table}"
    );
    assert!(table.contains("latitude"), "{table}");

    // The same fact holds structurally in the JSON form, not just as table prose.
    let (json, stderr, ok) = run(&[
        "explain", "--a", a_path.to_str().unwrap(), "--b", b_path.to_str().unwrap(),
        "--ranges", "600", "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&json).unwrap();
    let skipped = v["skipped_axes"].as_array().unwrap();
    assert!(!skipped.is_empty(), "{json}");
    assert!(skipped.iter().any(|s| s["axis"] == "latitude"), "{json}");
}

/// Review I2: `explain` reads TWO files, so a decode/solve error must name WHICH one failed --
/// `describe_solve_error` alone has no path, so without this the error for a broken `--a` is
/// indistinguishable from the same error for a broken `--b`. Checked in both positions so the
/// fix is not just "always blames --a" (or always blames whichever is read first).
#[test]
fn explain_names_which_of_the_two_files_is_broken() {
    let dir = tempfile_dir();
    let mut broken = small_request(823.0);
    // The same "g7" lowercase defect this file's own `small_request` doc warns about --
    // reused here deliberately, since it is a real, guaranteed-to-fail `InvalidValue`.
    broken["projectile"]["drag_model"] = serde_json::json!("g7");
    let broken_path = write_json(&dir, "broken.json", &broken);
    let good_path = write_json(&dir, "good.json", &small_request(870.0));
    let (broken_str, good_str) = (broken_path.to_str().unwrap(), good_path.to_str().unwrap());

    let (_, stderr, ok) = run(&["explain", "--a", broken_str, "--b", good_str, "--ranges", "600"]);
    assert!(!ok, "a broken --a must be rejected");
    assert!(
        stderr.contains(broken_str),
        "the error must name the broken file's own path: {stderr}"
    );
    assert!(
        !stderr.contains(good_str),
        "the error must not also name the OTHER, valid file: {stderr}"
    );

    let (_, stderr, ok) = run(&["explain", "--a", good_str, "--b", broken_str, "--ranges", "600"]);
    assert!(!ok, "a broken --b must be rejected");
    assert!(
        stderr.contains(broken_str),
        "the error must name the broken file's own path regardless of which flag carried it: \
         {stderr}"
    );
    assert!(!stderr.contains(good_str), "{stderr}");
}

/// Delta #3, situations 1 and 3 together (a found bound vs. no measurable effect at all): the
/// vacuum fixture gives an ANALYTICALLY known crossing for `muzzle-velocity-mps` (drop-only,
/// since there is no wind), paired with `target-distance` in the SAME call, which -- being
/// `requires_rezero: false` -- can never move the impact observed at a fixed `--range` at all.
/// Both directions of both axes share the identical `None`/`None` `_bound` shape in the library
/// return; the property under test is that the table renders them with DIFFERENT text, never
/// the generic "unbounded" phrasing for the axis that has no effect at all.
#[test]
fn tolerance_table_distinguishes_a_found_bound_from_no_measurable_effect() {
    let dir = tempfile_dir();
    let req = write_json(&dir, "vacuum.json", &vacuum_request());
    let req_str = req.to_str().unwrap();

    let args = [
        "--units", "metric", "tolerance", "--request", req_str,
        "--range", "400", "--target", "rect:20x20",
        "--axis", "muzzle-velocity-mps", "--domain", "muzzle-velocity-mps=700:1100",
        "--axis", "target-distance", "--domain", "target-distance=450:600",
    ];

    let (table, stderr, ok) = run(&args);
    assert!(ok, "stderr: {stderr}");
    assert!(table.contains("crosses the Bottom edge"), "{table}");
    assert!(table.contains("crosses the Top edge"), "{table}");
    assert!(table.contains("no measurable effect on the impact"), "{table}");
    assert!(
        !table.contains("no bound within the configured domain"),
        "neither axis in this fixture is merely unbounded-with-effect, so that phrase must not \
         appear at all here (it would make the two situations indistinguishable): {table}"
    );
    // Review I1: target-distance is ALSO `unbounded_in_domain: true` (both bounds `None`), but
    // with NO effect in either direction -- it must not additionally get the "stays inside the
    // target" summary line, which would misread as reassurance about a change that provably
    // does nothing at all. muzzle-velocity-mps has real bounds, so it was never eligible for
    // that line either; this fixture's whole table should therefore contain it zero times.
    assert!(
        !table.contains("unbounded in domain:"),
        "no axis in this fixture should get the unbounded-in-domain summary line -- \
         muzzle-velocity-mps has real bounds, and target-distance has no effect at all: {table}"
    );

    let mut json_args: Vec<&str> = args.to_vec();
    json_args.extend(["-o", "json"]);
    let (stdout, stderr, ok) = run(&json_args);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let axes = v["axes"].as_array().unwrap();

    let mv = axes
        .iter()
        .find(|a| a["axis"] == "muzzle_velocity_mps")
        .expect("muzzle-velocity-mps axis in report");
    assert!(mv["near_bound"].is_number() && mv["far_bound"].is_number(), "{mv}");
    assert_eq!(mv["near_limiting_boundary"], "bottom");
    assert_eq!(mv["far_limiting_boundary"], "top");
    assert_eq!(mv["unbounded_in_domain"], false);
    assert_eq!(mv["near_has_no_effect"], false);
    assert_eq!(mv["far_has_no_effect"], false);

    let td = axes
        .iter()
        .find(|a| a["axis"] == "target_distance")
        .expect("target-distance axis in report");
    assert!(td["near_bound"].is_null() && td["far_bound"].is_null(), "{td}");
    assert_eq!(td["unbounded_in_domain"], true);
    assert_eq!(td["near_has_no_effect"], true);
    assert_eq!(td["far_has_no_effect"], true);
}

/// Delta #3, situations 2 and 5 together (unbounded-but-has-effect vs. structurally
/// unavailable): a tiny wind-speed domain around a calm-ish 3 m/s nominal never leaves a huge
/// 50x50 m target (so it genuinely has an effect, it just never leaves this target), paired
/// with a categorical axis (`magnus-enabled`), which `tolerance_envelope` checks BEFORE ever
/// consulting `--domain` at all and reports as unavailable, not as a zero-width bound. A
/// `--domain` is still required and supplied for it by this CLI's own uniform contract (its
/// bounds simply go unused).
#[test]
fn tolerance_table_distinguishes_unbounded_from_unavailable() {
    let dir = tempfile_dir();
    let req = write_json(&dir, "req.json", &small_request(823.0));
    let req_str = req.to_str().unwrap();

    let args = [
        "--units", "metric", "tolerance", "--request", req_str,
        "--range", "600", "--target", "rect:5000x5000",
        "--axis", "wind-speed", "--domain", "wind-speed=2.9:3.1",
        "--axis", "magnus-enabled", "--domain", "magnus-enabled=0:1",
    ];

    let (table, stderr, ok) = run(&args);
    assert!(ok, "stderr: {stderr}");
    assert!(table.contains("no bound within the configured domain"), "{table}");
    assert!(
        !table.contains("no measurable effect on the impact"),
        "wind-speed does move the impact within this domain, it just never leaves this huge \
         target: {table}"
    );
    // Review I1's positive counterpart: wind-speed IS a genuine "unbounded but has effect" axis
    // (unlike target-distance in the sibling test above), so it must still get the summary
    // line -- the fix must not over-suppress this, only the both-no-effect case.
    assert!(
        table.contains("unbounded in domain:"),
        "wind-speed has a real effect in both directions and never leaves this huge target, so \
         it must still get the summary line: {table}"
    );
    assert!(table.contains("unavailable axes"), "{table}");
    assert!(table.contains("could not be searched at all"), "{table}");
    assert!(table.contains("magnus-enabled"), "{table}");

    let mut json_args: Vec<&str> = args.to_vec();
    json_args.extend(["-o", "json"]);
    let (stdout, stderr, ok) = run(&json_args);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();

    let ws = v["axes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["axis"] == "wind_speed")
        .expect("wind-speed axis in report");
    assert_eq!(ws["unbounded_in_domain"], true);
    assert_eq!(ws["near_has_no_effect"], false);
    assert_eq!(ws["far_has_no_effect"], false);

    let unavailable = v["unavailable_axes"].as_array().unwrap();
    assert!(
        unavailable
            .iter()
            .any(|u| u["axis"] == "magnus_enabled" && u["code"] == "categorical_axis"),
        "{unavailable:?}"
    );
    // An unavailable axis never also appears in `axes` -- the two lists are disjoint.
    assert!(
        v["axes"].as_array().unwrap().iter().all(|a| a["axis"] != "magnus_enabled"),
        "{v}"
    );
}

// ---- Task 15: error-budget (`error_budget_with_target`, MBA-1347) ----

/// Task-15 brief's own test (Step 1), corrected for the REAL canonical axis spelling:
/// `muzzle-velocity-mps` (from `InputAxis::MuzzleVelocityMps` via `axis_kebab_name`), not the
/// brief's stale `muzzle-velocity` -- see the task-14 report's own note flagging this for Task
/// 15. Extended with the schema_version/method checks the task instructions ask for on every
/// subcommand, plus (since this call supplies `--target`) the p_hit/gain presence checks.
#[test]
fn error_budget_ranks_sources_and_states_assumptions() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (stdout, stderr, ok) = run(&[
        "error-budget", "--request", a.to_str().unwrap(), "--ranges", "600",
        "--sigma", "muzzle-velocity-mps=5.0", "--sigma", "wind-speed=1.5",
        "--target", "rect:0.5x0.75", "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).expect("json");
    assert_eq!(v["schema_version"], 1);
    assert_eq!(
        v["method"],
        "central_difference_first_order_propagation_gl20_panelled_pm6sigma"
    );
    let sources = v["rows"][0]["sources"].as_array().unwrap();
    assert_eq!(sources.len(), 2, "both declared sources must be preserved");
    // Ranked descending by variance share.
    let first = sources[0]["variance_share"].as_f64().unwrap();
    let second = sources[1]["variance_share"].as_f64().unwrap();
    assert!(first >= second);
    assert!(
        v["rows"][0]["p_hit"].is_number(),
        "a target was supplied: p_hit must be reported, not null: {v}"
    );
    assert!(
        sources[0]["p_hit_gain_if_perfect"].is_number(),
        "a target was supplied: the per-source gain must be reported, not null: {v}"
    );
    assert!(v["rows"][0]["priority_statement"].as_str().unwrap().len() > 10);
    assert!(v["assumptions"].as_array().unwrap().iter()
        .any(|s| s.as_str().unwrap().contains("independent")));
}

/// Correction #3 (task-15 brief): the table must show `p_hit` and each source's
/// hit-probability gain ONLY when a target was supplied -- never a fabricated number when one
/// was not, and never silently missing when one WAS. A single cheap (`WindSpeed`, not
/// `requires_rezero`) source is enough to exercise this; ranking multiple sources is already
/// covered by the test above.
#[test]
fn error_budget_table_shows_p_hit_and_gain_only_with_a_target() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let a_str = a.to_str().unwrap();

    let (with_target, stderr, ok) = run(&[
        "error-budget", "--request", a_str, "--ranges", "600",
        "--sigma", "wind-speed=1.5", "--target", "rect:0.5x0.75",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(with_target.contains("p_hit:"), "{with_target}");
    assert!(with_target.contains("gain if perfected"), "{with_target}");

    let (without_target, stderr, ok) = run(&[
        "error-budget", "--request", a_str, "--ranges", "600", "--sigma", "wind-speed=1.5",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(!without_target.contains("p_hit:"), "{without_target}");
    assert!(!without_target.contains("gain if perfected"), "{without_target}");
}

/// The unavailable-source path, end to end: `Altitude` is declared alongside a real source
/// (`MuzzleVelocityMps`) on the QNH fixture. The table must name the unavailable axis with its
/// reason, and it must NEVER appear in the ranked list -- a silently dropped source would read
/// as "contributes no uncertainty," the exact wrong answer `error_budget.rs`'s own module doc
/// says this feature exists to prevent. Checked structurally in the JSON too, not just by
/// scraping table text.
#[test]
fn error_budget_names_an_unavailable_source_and_excludes_it_from_the_ranked_list() {
    let dir = tempfile_dir();
    let req = write_json(&dir, "qnh.json", &qnh_request());
    let req_str = req.to_str().unwrap();

    let (table, stderr, ok) = run(&[
        "error-budget", "--request", req_str, "--ranges", "300",
        "--sigma", "altitude=50.0", "--sigma", "muzzle-velocity-mps=5.0",
    ]);
    assert!(ok, "stderr: {stderr}");
    assert!(table.contains("unavailable sources"), "{table}");
    assert!(table.contains("altitude"), "{table}");
    assert!(table.contains("QNH"), "the reason should name the mechanism: {table}");
    let ranked_section = table.split("unavailable sources").next().unwrap();
    assert!(
        !ranked_section.contains("altitude"),
        "an unavailable source must never appear in the ranked list: {table}"
    );
    assert!(ranked_section.contains("muzzle-velocity-mps"), "{table}");

    let (stdout, stderr, ok) = run(&[
        "error-budget", "--request", req_str, "--ranges", "300",
        "--sigma", "altitude=50.0", "--sigma", "muzzle-velocity-mps=5.0", "-o", "json",
    ]);
    assert!(ok, "stderr: {stderr}");
    let v: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    let unavailable = v["unavailable_sources"].as_array().unwrap();
    assert_eq!(unavailable.len(), 1, "{v}");
    assert_eq!(unavailable[0]["axis"], "altitude");
    assert_eq!(unavailable[0]["code"], "axis_unsupported_for_request");

    let sources = v["rows"][0]["sources"].as_array().unwrap();
    assert_eq!(sources.len(), 1, "only the one evaluable source should be ranked: {v}");
    assert_eq!(sources[0]["axis"], "muzzle_velocity_mps");
}

/// A duplicate `--sigma` axis would double-count that axis's variance and corrupt its own
/// leave-one-out counterfactual (`error_budget_with_target`'s own contract) -- rejected by the
/// CLI before any solve, naming the axis and the flag, same as `tolerance --axis`/`--domain`.
#[test]
fn error_budget_rejects_a_duplicate_sigma_axis() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "error-budget", "--request", a.to_str().unwrap(), "--ranges", "600",
        "--sigma", "wind-speed=1.5", "--sigma", "wind-speed=2.0",
    ]);
    assert!(!ok, "a duplicate --sigma axis must be rejected");
    assert!(stderr.contains("wind-speed"), "{stderr}");
    assert!(stderr.contains("more than once"), "{stderr}");
}

/// A negative sigma is physically meaningless (a one-sigma uncertainty is a magnitude) and is
/// rejected by name -- not forwarded to `KernelError::NonFinite`'s generic Display text, which
/// describes an unrelated failure (a computed derivative that came out non-finite). See
/// `parse_sigma_arg`'s own doc comment in `src/main.rs`.
#[test]
fn error_budget_rejects_a_negative_sigma() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let (_, stderr, ok) = run(&[
        "error-budget", "--request", a.to_str().unwrap(), "--ranges", "600",
        "--sigma", "wind-speed=-1.0",
    ]);
    assert!(!ok, "a negative sigma must be rejected");
    assert!(stderr.contains("wind-speed"), "{stderr}");
}

/// CSV and PDF have no form for any of the three subcommands (none of the three reports is
/// tabular in a CSV-friendly sense, matching `hold-corridor`'s identical rejection).
#[test]
fn csv_and_pdf_are_rejected_for_all_three_subcommands() {
    let dir = tempfile_dir();
    let a = write_json(&dir, "a.json", &small_request(823.0));
    let b = write_json(&dir, "b.json", &small_request(870.0));

    for format in ["csv", "pdf"] {
        let (_, stderr, ok) = run(&[
            "explain", "--a", a.to_str().unwrap(), "--b", b.to_str().unwrap(),
            "--ranges", "600", "-o", format,
        ]);
        assert!(!ok, "explain -o {format} was accepted");
        assert!(stderr.contains("no "), "{format}: {stderr}");

        let (_, stderr, ok) = run(&[
            "tolerance", "--request", a.to_str().unwrap(), "--range", "600",
            "--target", "rect:0.5x0.75", "--axis", "wind-speed", "--domain", "wind-speed=0:20",
            "-o", format,
        ]);
        assert!(!ok, "tolerance -o {format} was accepted");
        assert!(stderr.contains("no "), "{format}: {stderr}");

        let (_, stderr, ok) = run(&[
            "error-budget", "--request", a.to_str().unwrap(), "--ranges", "600",
            "--sigma", "wind-speed=1.5", "-o", format,
        ]);
        assert!(!ok, "error-budget -o {format} was accepted");
        assert!(stderr.contains("no "), "{format}: {stderr}");
    }
}

/// A missing request file names the path, for all three subcommands, matching
/// `hold-corridor`'s own "could not read" convention.
#[test]
fn a_missing_request_file_names_the_path() {
    let (_, stderr, ok) = run(&[
        "explain", "--a", "/nonexistent/a.json", "--b", "/nonexistent/b.json", "--ranges", "600",
    ]);
    assert!(!ok);
    assert!(stderr.contains("/nonexistent/a.json"), "{stderr}");

    let (_, stderr, ok) = run(&[
        "tolerance", "--request", "/nonexistent/req.json", "--range", "600",
        "--target", "rect:0.5x0.75", "--axis", "wind-speed", "--domain", "wind-speed=0:20",
    ]);
    assert!(!ok);
    assert!(stderr.contains("/nonexistent/req.json"), "{stderr}");

    let (_, stderr, ok) = run(&[
        "error-budget", "--request", "/nonexistent/req.json", "--ranges", "600",
        "--sigma", "wind-speed=1.5",
    ]);
    assert!(!ok);
    assert!(stderr.contains("/nonexistent/req.json"), "{stderr}");
}
