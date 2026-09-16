//! An omitted `rifle.twist_rate_m_per_turn` under a spin-driven effect (MBA-1484).
//!
//! `twist_rate_m_per_turn` is optional, which reads like "leave the twist out of it". It is
//! not: `resolve_rifle` substitutes 0.3048 m — exactly 1:12 inches — and every model that
//! reads the twist then runs against that barrel. The early return in
//! `TrajectorySolver::apply_spin_drift` for a non-positive twist never fires on this path,
//! because the default is applied long before the solver sees the field.
//!
//! Before these warnings existed the only trace of that was a `default_applied` assumption
//! notice, which says a default was applied without saying that anything now depends on it.
//! A caller enabling `enhanced_spin_drift` on a 1:8 barrel they never described got a
//! confident drift number for a 1:12 one.
//!
//! MBA-1484 adds two codes. `spin_effect_assumed_twist_rate` covers the two opt-in spin
//! effects. `stability_factor_assumed_twist_rate` covers `summary.stability_factor`, which is
//! not opt-in at all — Sg is a function of the twist and is computed on EVERY solve, so an
//! omitted field publishes the Sg of a 1:12 barrel no matter what `effects` says. Together
//! with the pre-existing `aerodynamic_jump_assumed_geometry` that makes three assumed-barrel
//! codes, all keyed to the CONSUMER rather than to the missing field, so more than one can
//! fire from a single omission. This pins the fix:
//!
//!   (a) the hazard itself — an omitted twist is BIT-IDENTICAL to a stated 1:12 and is not
//!       the same answer as a stated 1:8, so the response cannot be told apart from a
//!       deliberate 1:12 request by its numbers alone;
//!   (b) `enhanced_spin_drift` without a twist raises `spin_effect_assumed_twist_rate` at
//!       `$.effects.enhanced_spin_drift`;
//!   (c) `magnus` without a twist raises it at `$.effects.magnus` — smaller in absolute
//!       terms on a flat-fire shot, equally twist-bound;
//!   (d) stating the twist raises neither code, for either flag;
//!   (e) omitting the twist with NO effect at all raises no spin-effect warning, because the
//!       trajectory really is identical for every twist rate — but it does raise the Sg
//!       warning, because Sg is not, and that is the case this file originally got wrong;
//!   (f) the codes are distinct but NOT mutually exclusive: `aerodynamic_jump` keeps its own
//!       `aerodynamic_jump_assumed_geometry`, and a request enabling jump and `magnus`
//!       together with the twist omitted gets all three at once, each at its own path;
//!   (g) the warnings are the whole fix — the request is still solved, and the 1:12 default
//!       is still materialized in `resolved_request`, byte-for-byte as before.

use ballistics_engine::solve_json::{decode_solve_request_v1, SolveSuccessV1};
use ballistics_engine::solve_v1;

/// A .308 175 gr at 800 m in a 10 mph pure crosswind from the right.
///
/// `rifle_extra` is the only thing that varies: the twist is either stated or left out, which
/// is the entire subject of this file.
fn request_json(effects: &str, rifle_extra: &str) -> String {
    format!(
        r#"{{
            "schema_version": 1,
            "projectile": {{"mass_kg": 0.01134, "diameter_m": 0.00782, "length_m": 0.0315,
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

/// 1:12 inches — the value `resolve_rifle` substitutes when the field is absent, stated out
/// loud so the two requests differ only in whether the caller said it.
const STATED_1_IN_12: &str = r#", "twist_rate_m_per_turn": 0.3048"#;
/// 1:8 inches: a real barrel a .308 shooter might own, and not the default.
const STATED_1_IN_8: &str = r#", "twist_rate_m_per_turn": 0.2032"#;

const CODE: &str = "spin_effect_assumed_twist_rate";
const SG_CODE: &str = "stability_factor_assumed_twist_rate";
const JUMP_CODE: &str = "aerodynamic_jump_assumed_geometry";

fn solve(json: &str) -> SolveSuccessV1 {
    let request = decode_solve_request_v1(json).expect("request must decode");
    solve_v1(request).expect("request must solve")
}

fn warnings_with_code<'a>(
    response: &'a SolveSuccessV1,
    code: &str,
) -> Vec<&'a ballistics_engine::solve_json::SolveNoticeV1> {
    response
        .warnings
        .iter()
        .filter(|n| n.code == code)
        .collect()
}

// (a)
#[test]
fn an_omitted_twist_solves_as_a_stated_1_in_12_and_not_as_the_real_barrel() {
    let omitted = solve(&request_json(r#""enhanced_spin_drift": true"#, ""));
    let stated_default = solve(&request_json(
        r#""enhanced_spin_drift": true"#,
        STATED_1_IN_12,
    ));
    let stated_fast = solve(&request_json(
        r#""enhanced_spin_drift": true"#,
        STATED_1_IN_8,
    ));

    let (omitted_drift, default_drift, fast_drift) = (
        omitted
            .summary
            .spin_drift_m
            .expect("drift must be reported"),
        stated_default
            .summary
            .spin_drift_m
            .expect("drift must be reported"),
        stated_fast
            .summary
            .spin_drift_m
            .expect("drift must be reported"),
    );

    // Bit-identical, not merely close: there is no numeric residue of the omission for a
    // caller to notice. This is why the warning has to carry the information.
    assert_eq!(
        omitted_drift, default_drift,
        "an omitted twist must solve exactly as a stated 1:12 does -- if these ever diverge, \
         the default is no longer 0.3048 m and this file's premise needs revisiting"
    );
    assert_ne!(
        omitted_drift, fast_drift,
        "a stated 1:8 must reach the drift, or the twist would not be load-bearing at all"
    );
    assert!(
        fast_drift > omitted_drift * 1.5,
        "the barrel alone should move the drift substantially: omitted {omitted_drift} m vs \
         1:8 {fast_drift} m"
    );
}

// (b)
#[test]
fn enhanced_spin_drift_without_a_twist_warns_that_the_barrel_was_assumed() {
    let response = solve(&request_json(r#""enhanced_spin_drift": true"#, ""));

    let assumed = warnings_with_code(&response, CODE);
    assert_eq!(
        assumed.len(),
        1,
        "the caller must be told the drift was computed from an assumed twist, got {:?}",
        response.warnings
    );
    assert_eq!(
        assumed[0].path.as_deref(),
        Some("$.effects.enhanced_spin_drift"),
        "the warning belongs at the flag that made the twist load-bearing"
    );

    // The point of the warning: the omission did not disable the effect.
    assert!(
        response.summary.spin_drift_m.is_some_and(|d| d != 0.0),
        "an assumed barrel still produces drift -- that is the hazard"
    );
}

// (c)
#[test]
fn magnus_without_a_twist_warns_that_the_barrel_was_assumed() {
    let response = solve(&request_json(r#""magnus": true"#, ""));

    let assumed = warnings_with_code(&response, CODE);
    assert_eq!(
        assumed.len(),
        1,
        "Magnus is driven by the spin the rifling imparts, so it depends on the twist as much \
         as the drift does, got {:?}",
        response.warnings
    );
    assert_eq!(assumed[0].path.as_deref(), Some("$.effects.magnus"));
}

// (d)
#[test]
fn a_stated_twist_raises_no_assumed_barrel_code_for_either_flag() {
    for effects in [r#""enhanced_spin_drift": true"#, r#""magnus": true"#] {
        for twist in [STATED_1_IN_12, STATED_1_IN_8] {
            let response = solve(&request_json(effects, twist));
            // Both codes, including the always-on Sg one: stating the twist is precisely what
            // makes every assumed-barrel warning inapplicable. Stating the DEFAULT value must
            // silence them too -- the warning is about what the caller said, not about which
            // number came out.
            for code in [CODE, SG_CODE] {
                assert!(
                    warnings_with_code(&response, code).is_empty(),
                    "the twist was stated ({effects}, {twist}), so {code} must not fire: {:?}",
                    response.warnings
                );
            }
        }
    }
}

// (e)
#[test]
fn an_omitted_twist_alone_spares_the_trajectory_but_not_the_reported_sg() {
    // With no spin effect enabled the TRAJECTORY is genuinely twist-independent, so the
    // spin-effect code would be describing a hazard that is not there.
    for effects in [
        "",
        r#""magnus": false, "enhanced_spin_drift": false"#,
        r#""coriolis": false"#,
    ] {
        let response = solve(&request_json(effects, ""));
        assert!(
            warnings_with_code(&response, CODE).is_empty(),
            "no spin effect ran in this request ({effects}): {:?}",
            response.warnings
        );
    }

    let omitted = solve(&request_json("", ""));
    let fast = solve(&request_json("", STATED_1_IN_8));
    let (omitted_last, fast_last) = (
        omitted.samples.last().expect("samples"),
        fast.samples.last().expect("samples"),
    );
    assert_eq!(omitted_last.drop_m, fast_last.drop_m);
    assert_eq!(omitted_last.windage_m, fast_last.windage_m);

    // But the SUMMARY is not twist-independent, which is the claim this file used to make and
    // which was false. Sg is computed on every solve from the resolved twist, so an omitted
    // field publishes the Sg of a barrel the caller never described -- and it is a far bigger
    // discrepancy than the Magnus case that does get a warning.
    let (omitted_sg, fast_sg) = (
        omitted.summary.stability_factor.expect("Sg is reported"),
        fast.summary.stability_factor.expect("Sg is reported"),
    );
    assert_ne!(
        omitted_sg, fast_sg,
        "if these ever match, Sg has stopped reading the twist and this warning is obsolete"
    );
    assert!(
        fast_sg > omitted_sg * 2.0,
        "the barrel alone moves the reported Sg across the marginal/comfortable line: \
         omitted {omitted_sg} vs 1:8 {fast_sg}"
    );

    // So it must be warned, on exactly the requests that publish such an Sg -- including this
    // one, which enables no effect at all.
    let sg_warnings = warnings_with_code(&omitted, SG_CODE);
    assert_eq!(
        sg_warnings.len(),
        1,
        "a solve that reports an Sg for an assumed barrel must say so even with no effects \
         enabled: {:?}",
        omitted.warnings
    );
    assert_eq!(
        sg_warnings[0].path.as_deref(),
        Some("$.rifle.twist_rate_m_per_turn"),
        "there is no flag to hang it on, so it names the omitted field itself"
    );

    // And it must NOT fire once the caller states the twist, or it is noise rather than a
    // warning.
    assert!(
        warnings_with_code(&fast, SG_CODE).is_empty(),
        "the twist was stated: {:?}",
        fast.warnings
    );
}

// (e), continued: the Sg code is independent of the effects flags, not a side effect of them.
#[test]
fn the_stability_factor_warning_tracks_the_twist_not_the_effects() {
    for effects in [
        "",
        r#""coriolis": false"#,
        r#""magnus": true"#,
        r#""enhanced_spin_drift": true"#,
        r#""aerodynamic_jump": true"#,
    ] {
        let omitted = solve(&request_json(effects, ""));
        assert_eq!(
            warnings_with_code(&omitted, SG_CODE).len(),
            1,
            "Sg was reported from the assumed barrel ({effects}): {:?}",
            omitted.warnings
        );

        let stated = solve(&request_json(effects, STATED_1_IN_8));
        assert!(
            warnings_with_code(&stated, SG_CODE).is_empty(),
            "the twist was stated ({effects}): {:?}",
            stated.warnings
        );
    }
}

// (f)
#[test]
fn aerodynamic_jump_alone_raises_its_own_code_and_not_the_spin_effect_one() {
    let response = solve(&request_json(r#""aerodynamic_jump": true"#, ""));

    assert_eq!(
        warnings_with_code(&response, JUMP_CODE).len(),
        1,
        "the jump's own warning still covers it: {:?}",
        response.warnings
    );
    // Not a claim that the codes are mutually exclusive -- see the test below, which shows
    // they are not. Neither spin effect is enabled in THIS request, so the spin-effect code
    // has no consumer to describe and would be reporting a model that did not run.
    assert!(
        warnings_with_code(&response, CODE).is_empty(),
        "no spin effect ran, so nothing should claim one was computed from an assumed \
         barrel: {:?}",
        response.warnings
    );
}

// (f), continued: the property the codes actually have. This is the case the previous
// revision of this file asserted was impossible, in a failure message, while the binary
// produced it.
#[test]
fn one_omitted_twist_raises_every_consumers_code_at_once() {
    let response = solve(&request_json(
        r#""magnus": true, "aerodynamic_jump": true"#,
        "",
    ));

    // Three consumers of one missing field -- Magnus, the jump, and the always-on Sg -- so
    // three warnings, each naming the consumer it belongs to rather than the field.
    for (code, path) in [
        (CODE, "$.effects.magnus"),
        (JUMP_CODE, "$.effects.aerodynamic_jump"),
        (SG_CODE, "$.rifle.twist_rate_m_per_turn"),
    ] {
        let found = warnings_with_code(&response, code);
        assert_eq!(
            found.len(),
            1,
            "{code} must be present exactly once -- the codes are keyed to the consumer, not \
             to the missing field, so they are NOT mutually exclusive and a caller must match \
             on the one it cares about: {:?}",
            response.warnings
        );
        assert_eq!(
            found[0].path.as_deref(),
            Some(path),
            "{code} belongs at {path}, so the caller can tell the consumers apart"
        );
    }

    // Every one of them names the same single omitted field, which is exactly why the codes
    // rather than the messages have to carry the distinction.
    for code in [CODE, JUMP_CODE, SG_CODE] {
        assert!(
            warnings_with_code(&response, code)[0]
                .message
                .contains("rifle.twist_rate_m_per_turn"),
            "{code} should name the field that was left out: {:?}",
            response.warnings
        );
    }
}

// (g)
#[test]
fn the_warning_is_the_whole_fix_and_changes_nothing_else() {
    let response = solve(&request_json(r#""enhanced_spin_drift": true"#, ""));

    // Not a rejection: solve-json v1 is a published contract, and requests written before
    // this warning existed must keep solving.
    assert_eq!(
        response.resolved_request.rifle.twist_rate_m_per_turn, 0.3048,
        "the default itself is untouched -- the silence was the defect, not the default"
    );
    assert!(
        response
            .assumptions
            .iter()
            .any(|n| n.code == "default_applied"
                && n.path.as_deref() == Some("$.rifle.twist_rate_m_per_turn")),
        "the pre-existing assumption notice stays: {:?}",
        response.assumptions
    );

    // The numbers are the same ones the request would have produced before the warning
    // existed -- checked against the explicit 1:12 request, which raises no new warning and
    // is therefore unaffected by this change.
    let stated = solve(&request_json(
        r#""enhanced_spin_drift": true"#,
        STATED_1_IN_12,
    ));
    assert_eq!(response.summary.spin_drift_m, stated.summary.spin_drift_m);
    assert_eq!(
        response.samples.last().expect("samples").drop_m,
        stated.samples.last().expect("samples").drop_m
    );
}
