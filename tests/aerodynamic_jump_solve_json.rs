//! `effects.aerodynamic_jump` on the solve-json v1 wire (MBA-959).
//!
//! The engine has applied Litz's crosswind aerodynamic jump since MBA-959, and the native CLI
//! (`--enable-aerodynamic-jump`) reaches it. The JSON solve path — the mobile embedding
//! surface — hard-coded `enable_aerodynamic_jump: false`, and `EffectsV1` is
//! `deny_unknown_fields`, so an app built on the bridge could not ask for it at all: the field
//! was a parse error, not a no-op. Alfredo Mendiola Loyola reported exactly that from an
//! Android app. This pins the wire contract that closes it:
//!
//!   (a) omitting the field is byte-identical to every earlier response, echo included — the
//!       key appears nowhere in the envelope;
//!   (b) an explicit `false` solves identically to omitting it, but is still echoed, so a
//!       round-tripped request says the same thing as the one that produced it;
//!   (c) `true` is echoed and raises the `experimental_effect` warning its two experimental
//!       siblings raise — the Litz form is a regression fitted near Sg ~ 1.75, not a
//!       derivation, and the bridge says so rather than shipping it silently;
//!   (d) the jump actually MOVES the numbers, vertically and only vertically, and the
//!       response reports how much in `summary.aerodynamic_jump_moa` — a field that exists
//!       because the jump is otherwise unobservable, being folded into every drop rather
//!       than appearing as its own column;
//!   (e) it is computed from the MUZZLE crosswind: a headwind produces a present, exactly
//!       zero jump, which is a different fact from the effect never running, and the response
//!       distinguishes them (`Some(0.0)` vs absent);
//!   (f) enabling it without `rifle.twist_rate_m_per_turn` or `projectile.length_m` does NOT
//!       disable or zero it — it computes one from an assumed 1:12 barrel — so that case
//!       raises `aerodynamic_jump_assumed_geometry`. This is the failure shape of MBA-1484
//!       arriving on a second effect, and the warning is the only thing standing between a
//!       caller and a confident number for a rifle they do not own;
//!   (g) `null` is rejected at the exact path, like every other optional boolean.

use ballistics_engine::solve_json::{decode_solve_request_v1, SolveSuccessV1};
use ballistics_engine::solve_v1;

/// A .308 175 gr at 800 m in a 10 mph pure crosswind from the right, 1:10 twist stated.
///
/// Everything the jump reads is supplied explicitly — twist, twist direction, bullet length —
/// so these assertions describe the estimator rather than the defaults. The
/// assumed-geometry test below is the one that drops them, deliberately.
fn request_json(effects: &str, rifle_extra: &str, projectile_extra: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782{projectile_extra},
                            "drag_model": "G7", "ballistic_coefficient": 0.243}},
            "rifle": {{"muzzle_velocity_mps": 792.48, "sight_height_m": 0.0508{rifle_extra}}},
            "shot": {{"max_range_m": 800.0, "zero_distance_m": 100.0}},
            "atmosphere": {{"altitude_m": 0.0, "temperature_k": 288.15, "pressure_pa": 101325.0,
                            "relative_humidity": 0.0}},
            "wind": {{"speed_mps": 4.4704, "direction_from_rad": 1.5707963267948966}},
            "solver": {{"method": "rk45"}},
            "effects": {{{effects}}},
            "sampling": {{"interval_m": 100.0}}
        }}"#
    )
}

/// Twist and length stated — the configuration every test here uses but the last two.
const STATED_TWIST: &str = r#", "twist_rate_m_per_turn": 0.254, "twist_direction": "right""#;
const STATED_LENGTH: &str = r#", "length_m": 0.0315"#;

fn solve(json: &str) -> SolveSuccessV1 {
    let request = decode_solve_request_v1(json).expect("request must decode");
    solve_v1(request).expect("request must solve")
}

fn stated(effects: &str) -> SolveSuccessV1 {
    solve(&request_json(effects, STATED_TWIST, STATED_LENGTH))
}

fn drop_at_terminal(response: &SolveSuccessV1) -> f64 {
    response
        .samples
        .last()
        .expect("a solved trajectory has samples")
        .drop_m
}

fn windage_at_terminal(response: &SolveSuccessV1) -> f64 {
    response
        .samples
        .last()
        .expect("a solved trajectory has samples")
        .windage_m
}

// (a)
#[test]
fn an_omitted_field_leaves_no_trace_in_the_envelope() {
    let omitted = stated("");

    assert_eq!(omitted.resolved_request.effects.aerodynamic_jump, None);
    assert_eq!(omitted.summary.aerodynamic_jump_moa, None);

    // Not just the echo: the key must not appear ANYWHERE, or a caller diffing responses
    // across an engine upgrade sees a change they did not ask for.
    let wire = serde_json::to_string(&omitted).expect("response must serialize");
    assert!(
        !wire.contains("aerodynamic_jump"),
        "an omitted flag must leave the response byte-identical to one from before the field \
         existed, but the envelope mentions it: {wire}"
    );
}

// (b)
#[test]
fn an_explicit_false_solves_like_omission_but_is_still_echoed() {
    let omitted = stated("");
    let explicit = stated(r#""aerodynamic_jump": false"#);

    assert_eq!(
        explicit.resolved_request.effects.aerodynamic_jump,
        Some(false),
        "a supplied false must survive into the echo -- this is why the resolved field is \
         Option<bool> rather than the bare bool its three older siblings use"
    );
    assert_eq!(explicit.summary.aerodynamic_jump_moa, None);
    assert_eq!(
        drop_at_terminal(&explicit),
        drop_at_terminal(&omitted),
        "an explicit false must solve exactly as an omitted field does"
    );
}

// (c)
#[test]
fn enabling_it_echoes_true_and_warns_that_the_model_is_experimental() {
    let enabled = stated(r#""aerodynamic_jump": true"#);

    assert_eq!(
        enabled.resolved_request.effects.aerodynamic_jump,
        Some(true)
    );

    let experimental: Vec<_> = enabled
        .warnings
        .iter()
        .filter(|notice| notice.code == "experimental_effect")
        .filter(|notice| notice.path.as_deref() == Some("$.effects.aerodynamic_jump"))
        .collect();
    assert_eq!(
        experimental.len(),
        1,
        "exactly one experimental_effect warning, at the flag's own path: {:?}",
        enabled.warnings
    );
}

// (d)
#[test]
fn the_jump_moves_the_drop_and_only_the_drop() {
    let off = stated("");
    let on = stated(r#""aerodynamic_jump": true"#);

    let delta_drop = drop_at_terminal(&on) - drop_at_terminal(&off);
    let delta_windage = windage_at_terminal(&on) - windage_at_terminal(&off);

    // A right-hand twist in a wind from the right jumps the bullet UP, so it drops LESS.
    // Measured at ~0.108 m for this load; the threshold is loose enough to survive an
    // integrator tweak and tight enough that a silently-inert flag fails.
    assert!(
        delta_drop < -0.05,
        "aerodynamic jump must raise the impact by a visible margin, but drop moved by \
         {delta_drop} m"
    );
    assert!(
        delta_windage.abs() < 1e-3,
        "aerodynamic jump is a vertical effect; windage moved by {delta_windage} m"
    );

    // And the response says how much, rather than leaving the caller to difference two solves.
    let reported = on
        .summary
        .aerodynamic_jump_moa
        .expect("an applied jump must be reported");
    assert!(
        reported > 0.3 && reported < 0.6,
        "Litz gives ~0.46 MOA for this load and wind; got {reported}"
    );
    assert_eq!(off.summary.aerodynamic_jump_moa, None);
}

// (e)
#[test]
fn a_headwind_gives_a_present_zero_rather_than_an_absent_jump() {
    let headwind = request_json(r#""aerodynamic_jump": true"#, STATED_TWIST, STATED_LENGTH)
        .replace(
            "\"direction_from_rad\": 1.5707963267948966",
            "\"direction_from_rad\": 0.0",
        );
    let response = solve(&headwind);

    let reported = response
        .summary
        .aerodynamic_jump_moa
        .expect("the effect ran; it simply had no crosswind to act on");
    assert_eq!(
        reported, 0.0,
        "no crosswind at the muzzle means no jump, exactly"
    );

    // The distinction this field exists for: ran-and-was-zero is not the same as never-ran,
    // and a caller debugging \"I enabled it and nothing changed\" needs to tell them apart.
    assert_ne!(response.summary.aerodynamic_jump_moa, None);
}

// (f)
#[test]
fn enabling_it_without_twist_warns_that_the_barrel_was_assumed() {
    for (rifle_extra, projectile_extra, what) in [
        ("", STATED_LENGTH, "twist omitted"),
        (STATED_TWIST, "", "length omitted"),
        ("", "", "both omitted"),
    ] {
        let response = solve(&request_json(
            r#""aerodynamic_jump": true"#,
            rifle_extra,
            projectile_extra,
        ));

        let assumed: Vec<_> = response
            .warnings
            .iter()
            .filter(|n| n.code == "aerodynamic_jump_assumed_geometry")
            .collect();
        assert_eq!(
            assumed.len(),
            1,
            "{what}: the caller must be told the jump was computed from assumed geometry, \
             got {:?}",
            response.warnings
        );
        assert_eq!(
            assumed[0].path.as_deref(),
            Some("$.effects.aerodynamic_jump")
        );

        // The point of the warning: the effect is NOT disabled by the omission. It produces a
        // confident number for a barrel the request never described.
        let reported = response
            .summary
            .aerodynamic_jump_moa
            .expect("an assumed barrel still produces a jump -- that is the hazard");
        assert!(reported > 0.0, "{what}: expected a non-zero jump");
    }
}

#[test]
fn a_fully_specified_request_raises_no_assumed_geometry_warning() {
    let response = stated(r#""aerodynamic_jump": true"#);
    assert!(
        !response
            .warnings
            .iter()
            .any(|n| n.code == "aerodynamic_jump_assumed_geometry"),
        "twist and length were both stated: {:?}",
        response.warnings
    );
}

/// Stating the twist has to reach the estimator, or the warning above would be the only thing
/// distinguishing a real barrel from the assumed one.
#[test]
fn the_stated_twist_changes_the_jump() {
    let fast = solve(&request_json(
        r#""aerodynamic_jump": true"#,
        r#", "twist_rate_m_per_turn": 0.2032, "twist_direction": "right""#,
        STATED_LENGTH,
    ));
    let slow = solve(&request_json(
        r#""aerodynamic_jump": true"#,
        r#", "twist_rate_m_per_turn": 0.3048, "twist_direction": "right""#,
        STATED_LENGTH,
    ));

    let (fast_moa, slow_moa) = (
        fast.summary.aerodynamic_jump_moa.expect("1:8 must report"),
        slow.summary.aerodynamic_jump_moa.expect("1:12 must report"),
    );
    assert!(
        fast_moa > slow_moa,
        "a faster twist raises Sg and so raises the Litz jump: 1:8 gave {fast_moa}, 1:12 gave \
         {slow_moa}"
    );
}

// (g)
#[test]
fn an_explicit_null_is_rejected_at_the_field_path() {
    let json = request_json(r#""aerodynamic_jump": null"#, STATED_TWIST, STATED_LENGTH);
    let error = decode_solve_request_v1(&json).expect_err("null must not decode");
    assert_eq!(
        error.error.path(),
        Some("$.effects.aerodynamic_jump"),
        "a null must fail at its own path, not somewhere generic: {error:?}"
    );
}

#[test]
fn a_non_boolean_is_rejected_at_the_field_path() {
    let json = request_json(r#""aerodynamic_jump": "yes""#, STATED_TWIST, STATED_LENGTH);
    let error = decode_solve_request_v1(&json).expect_err("a string must not decode");
    assert_eq!(error.error.path(), Some("$.effects.aerodynamic_jump"));
}
