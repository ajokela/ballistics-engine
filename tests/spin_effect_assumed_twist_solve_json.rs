//! An omitted `rifle.twist_rate_m_per_turn` under a spin-driven effect (MBA-1484).
//!
//! `twist_rate_m_per_turn` is optional, which reads like "leave the twist out of it". It is
//! not: `resolve_rifle` substitutes 0.3048 m — exactly 1:12 inches — and every model that
//! reads the twist then runs against that barrel. The early return in
//! `TrajectorySolver::apply_spin_drift` for a non-positive twist never fires on this path,
//! because the default is applied long before the solver sees the field.
//!
//! Before this warning existed the only trace of that was a `default_applied` assumption
//! notice, which says a default was applied without saying that anything now depends on it.
//! A caller enabling `enhanced_spin_drift` on a 1:8 barrel they never described got a
//! confident drift number for a 1:12 one. This pins the fix:
//!
//!   (a) the hazard itself — an omitted twist is BIT-IDENTICAL to a stated 1:12 and is not
//!       the same answer as a stated 1:8, so the response cannot be told apart from a
//!       deliberate 1:12 request by its numbers alone;
//!   (b) `enhanced_spin_drift` without a twist raises `spin_effect_assumed_twist_rate` at
//!       `$.effects.enhanced_spin_drift`;
//!   (c) `magnus` without a twist raises it at `$.effects.magnus` — smaller in absolute
//!       terms on a flat-fire shot, equally twist-bound;
//!   (d) stating the twist raises nothing, for either flag;
//!   (e) omitting the twist with no twist-reading effect raises nothing either: the
//!       trajectory is identical for every twist rate when nothing consumes it, and every
//!       request written before this warning existed must keep its warning list;
//!   (f) `aerodynamic_jump` keeps its own distinct code and does not also raise this one, so
//!       the code identifies which model was computed against the assumed barrel;
//!   (g) the warning is the whole fix — the request is still solved, and the 1:12 default is
//!       still materialized in `resolved_request`, byte-for-byte as before.

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
        "Magnus is driven by the spin the rifling imparts, so it depends on the twist just as \
         the drift does, got {:?}",
        response.warnings
    );
    assert_eq!(assumed[0].path.as_deref(), Some("$.effects.magnus"));
}

// (d)
#[test]
fn a_stated_twist_raises_nothing_for_either_flag() {
    for effects in [r#""enhanced_spin_drift": true"#, r#""magnus": true"#] {
        for twist in [STATED_1_IN_12, STATED_1_IN_8] {
            let response = solve(&request_json(effects, twist));
            assert!(
                warnings_with_code(&response, CODE).is_empty(),
                "the twist was stated ({effects}, {twist}): {:?}",
                response.warnings
            );
        }
    }
}

// (e)
#[test]
fn an_omitted_twist_alone_raises_nothing() {
    // Not an oversight: with no twist-reading effect enabled the twist reaches nothing, and
    // the trajectory is identical for every value of it. Warning here would fire on every
    // request written before this code existed and would describe no hazard.
    for effects in [
        "",
        r#""magnus": false, "enhanced_spin_drift": false"#,
        r#""coriolis": false"#,
    ] {
        let response = solve(&request_json(effects, ""));
        assert!(
            warnings_with_code(&response, CODE).is_empty(),
            "nothing in this request reads the twist ({effects}): {:?}",
            response.warnings
        );
    }

    // And the claim behind that: with no spin effect, the twist genuinely does not move the
    // numbers, so its absence is not information the caller is missing.
    let omitted = solve(&request_json("", ""));
    let fast = solve(&request_json("", STATED_1_IN_8));
    let (omitted_last, fast_last) = (
        omitted.samples.last().expect("samples"),
        fast.samples.last().expect("samples"),
    );
    assert_eq!(omitted_last.drop_m, fast_last.drop_m);
    assert_eq!(omitted_last.windage_m, fast_last.windage_m);
}

// (f)
#[test]
fn aerodynamic_jump_keeps_its_own_code_and_does_not_also_raise_this_one() {
    let response = solve(&request_json(r#""aerodynamic_jump": true"#, ""));

    assert_eq!(
        warnings_with_code(&response, "aerodynamic_jump_assumed_geometry").len(),
        1,
        "the jump's own warning still covers it: {:?}",
        response.warnings
    );
    assert!(
        warnings_with_code(&response, CODE).is_empty(),
        "one omitted field must not produce two warnings saying the same thing -- the codes \
         are distinct so a caller can tell WHICH model used the assumed barrel: {:?}",
        response.warnings
    );
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
