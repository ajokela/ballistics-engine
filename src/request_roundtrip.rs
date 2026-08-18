//! Reverse conversion from the canonical resolved request back into a solvable request
//! (Phase 0 of the 0.33.0 decision-support train).
//!
//! The resolved request was otherwise output-only. The perturbation kernel needs to take a
//! resolved request, change one input, and re-solve, which is impossible without this
//! direction. Every resolved value is carried across explicitly: a silently lossy
//! conversion would misattribute the dropped field's effect to whatever the caller happened
//! to be perturbing.
//!
//! Two fields are deliberately NOT carried straight across, because the resolved sibling
//! they'd ride along with is already the post-transform value, and re-supplying the
//! original input mode would apply that transform a second time:
//!
//! - `atmosphere.pressure_reference`: when the original request declared `"qnh"`,
//!   [`ResolvedAtmosphereV1::pressure_pa`] is already the REDUCED absolute station pressure
//!   (see `resolve_atmosphere`'s QNH branch) -- echoing `"qnh"` back alongside that
//!   already-reduced value would reduce it a second time. The rebuilt request always states
//!   `pressure_pa` as absolute (the omitted-field default), which is what the resolved
//!   value already is.
//! - `wind.wind_reference`: when the original request declared `"compass"`, every resolved
//!   wind direction is already converted to shooter-relative (see `resolve_wind`'s
//!   `to_relative`) -- echoing `"compass"` back alongside an already-relative direction
//!   would re-reference it against the shot azimuth a second time. The rebuilt request
//!   always states directions as shooter-relative (the omitted-field default), which is
//!   what the resolved values already are.
//!
//! `ResolvedShotV1` carries both `zero_distance_m` (caller intent) and `muzzle_angle_rad`
//! (the effective angle after zeroing). Both are carried onto the rebuilt request: an
//! explicit `muzzle_angle_rad` always takes priority over `zero_distance_m` at resolve time
//! (see `resolve_shot` and `solve_v1`'s zero-search gate), so supplying both reproduces the
//! exact original angle -- bit-identical, not just numerically re-converged -- while still
//! preserving the original zeroing intent as metadata rather than dropping it.

use crate::solve_json::*;

impl From<&ResolvedSolveRequestV1> for SolveRequestV1 {
    fn from(r: &ResolvedSolveRequestV1) -> Self {
        SolveRequestV1 {
            schema_version: SchemaVersionV1,
            projectile: ProjectileV1 {
                mass_kg: r.projectile.mass_kg,
                diameter_m: r.projectile.diameter_m,
                length_m: r.projectile.length_m,
                drag_model: r.projectile.drag_model,
                ballistic_coefficient: r.projectile.ballistic_coefficient,
            },
            rifle: RifleV1 {
                muzzle_velocity_mps: r.rifle.muzzle_velocity_mps,
                sight_height_m: Some(r.rifle.sight_height_m),
                muzzle_height_m: Some(r.rifle.muzzle_height_m),
                twist_rate_m_per_turn: Some(r.rifle.twist_rate_m_per_turn),
                twist_direction: Some(r.rifle.twist_direction),
                sight_offset_lateral_m: r.rifle.sight_offset_lateral_m,
            },
            shot: ShotV1 {
                max_range_m: r.shot.max_range_m,
                zero_distance_m: r.shot.zero_distance_m,
                // Both are carried: zero_distance_m is caller intent, muzzle_angle_rad is
                // the angle actually integrated after zeroing. An explicit muzzle_angle_rad
                // always wins at resolve time (resolve_shot / solve_v1), so this reproduces
                // the exact original angle with no re-zero, whether or not zero_distance_m
                // is also present.
                muzzle_angle_rad: Some(r.shot.muzzle_angle_rad),
                aim_azimuth_rad: Some(r.shot.aim_azimuth_rad),
                shot_azimuth_rad: Some(r.shot.shot_azimuth_rad),
                shooting_angle_rad: Some(r.shot.shooting_angle_rad),
                cant_angle_rad: Some(r.shot.cant_angle_rad),
                target_height_m: Some(r.shot.target_height_m),
                ground_threshold_m: Some(r.shot.ground_threshold_m),
                zero_poi_up_m: r.shot.zero_poi_up_m,
                zero_poi_right_m: r.shot.zero_poi_right_m,
                drops_reference: r.shot.drops_reference,
            },
            atmosphere: AtmosphereV1 {
                altitude_m: Some(r.atmosphere.altitude_m),
                temperature_k: Some(r.atmosphere.temperature_k),
                pressure_pa: Some(r.atmosphere.pressure_pa),
                // See the module doc: pressure_pa above is already absolute station
                // pressure; echoing a "qnh" reference back would reduce it a second time.
                pressure_reference: None,
                relative_humidity: Some(r.atmosphere.relative_humidity),
                latitude_rad: r.atmosphere.latitude_rad,
            },
            wind: wind_from_resolved(&r.wind),
            solver: SolverV1 {
                method: Some(r.solver.method),
                time_step_s: Some(r.solver.time_step_s),
            },
            effects: EffectsV1 {
                magnus: Some(r.effects.magnus),
                coriolis: Some(r.effects.coriolis),
                enhanced_spin_drift: Some(r.effects.enhanced_spin_drift),
            },
            sampling: SamplingV1 {
                interval_m: Some(r.sampling.interval_m),
            },
            reticle: r.reticle.clone(),
            // Carried straight across: segments are always regenerated from the PUBLISHED
            // ballistic_coefficient (see `apply_bc5d_correction`), so re-applying the
            // table on a re-solve reproduces — never compounds — the correction. Dropping
            // it instead would misattribute the correction's whole effect to whatever a
            // perturbation caller happened to be perturbing.
            corrections: r.corrections.clone(),
        }
    }
}

fn wind_from_resolved(w: &ResolvedWindV1) -> WindV1 {
    match w {
        ResolvedWindV1::Constant(c) => WindV1 {
            speed_mps: Some(c.speed_mps),
            direction_from_rad: Some(c.direction_from_rad),
            vertical_speed_mps: Some(c.vertical_speed_mps),
            segments: None,
            // See the module doc: direction_from_rad above is already shooter-relative;
            // echoing a "compass" reference back would re-reference it a second time.
            wind_reference: None,
        },
        ResolvedWindV1::Segmented(s) => WindV1 {
            speed_mps: None,
            direction_from_rad: None,
            vertical_speed_mps: None,
            segments: Some(
                s.segments
                    .iter()
                    .map(|g| WindSegmentV1 {
                        until_distance_m: g.until_distance_m,
                        speed_mps: g.speed_mps,
                        direction_from_rad: g.direction_from_rad,
                        vertical_speed_mps: Some(g.vertical_speed_mps),
                    })
                    .collect(),
            ),
            // See the module doc and the constant-wind arm above.
            wind_reference: None,
        },
    }
}

#[cfg(test)]
mod tests {
    use crate::solve_json::decode_solve_request_v1;
    use crate::solve_json::{PressureReferenceV1, ResolvedWindV1, SolveRequestV1, WindReferenceV1};
    use crate::solve_v1::solve_v1;

    fn sample_json() -> String {
        serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"temperature_k": 288.0, "pressure_pa": 101325.0},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string()
    }

    /// Resolution is idempotent through a round-trip: re-solving a request
    /// rebuilt from a resolved request must resolve to exactly the same values.
    /// This is the acceptance gate for Phase 0.
    ///
    /// Compares the whole success envelope, not just `resolved_request`: some effects (the
    /// windage-zero convergence bias, see `roundtrip_preserves_the_windage_zero_bias` below)
    /// never appear in `resolved_request` at all -- on the very first solve, not only after a
    /// round-trip -- so `resolved_request` equality alone cannot catch every regression.
    #[test]
    fn resolution_is_idempotent_through_roundtrip() {
        let first = solve_v1(decode_solve_request_v1(&sample_json()).unwrap()).unwrap();
        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        let second = solve_v1(rebuilt).unwrap();
        assert_eq!(
            serde_json::to_value(&first.resolved_request).unwrap(),
            serde_json::to_value(&second.resolved_request).unwrap(),
            "resolved request changed after a round-trip"
        );
        assert_eq!(
            first.summary, second.summary,
            "summary changed after a round-trip"
        );
        assert_eq!(
            first.samples, second.samples,
            "samples changed after a round-trip"
        );
    }

    /// The windage-zero convergence bias (`sight_offset_lateral_m` / `zero_poi_right_m`,
    /// applied via `BallisticInputs::windage_zero_bias_rad`) is a term
    /// `calculate_and_set_zero_angle` adds to azimuth ALONGSIDE the elevation search -- it is
    /// not carried by `muzzle_angle_rad`, and (unlike the elevation) it never appears in
    /// `resolved_request` at all, on the first solve or any later one. Skipping the elevation
    /// search on a round-tripped request must not also skip this separate term, or an
    /// offset-mounted sight / deliberate horizontal zero bias would silently stop converging
    /// the moment a resolved request round-trips. `resolved_request` alone cannot see this
    /// (compare `first.resolved_request` above with `second.resolved_request` below: they are
    /// byte-identical even when this regresses), so this compares the solved trajectory too.
    #[test]
    fn roundtrip_preserves_the_windage_zero_bias() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                      "sight_offset_lateral_m": 0.03},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0, "zero_poi_right_m": 0.02},
            "atmosphere": {"temperature_k": 288.0, "pressure_pa": 101325.0},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let first = solve_v1(decode_solve_request_v1(&json).unwrap()).unwrap();
        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        let second = solve_v1(rebuilt).unwrap();

        assert_eq!(
            serde_json::to_value(&first.resolved_request).unwrap(),
            serde_json::to_value(&second.resolved_request).unwrap(),
            "resolved request changed after a round-trip"
        );
        assert_eq!(
            first.summary, second.summary,
            "summary changed after a round-trip -- the windage-zero bias may have been dropped"
        );
        assert_eq!(
            first.samples, second.samples,
            "samples changed after a round-trip -- the windage-zero bias may have been dropped"
        );
        // Sanity check: the bias is actually nonzero in this fixture, so the comparisons
        // above are exercising the real thing rather than two agreeing zeros.
        let windage_m = first
            .samples
            .last()
            .expect("at least one sample")
            .windage_m;
        assert!(
            windage_m.abs() > 0.1,
            "fixture must produce a non-negligible windage-zero bias to be a meaningful test, \
             got {windage_m} m"
        );
    }

    /// A zeroed solve must not silently re-zero on the way back.
    #[test]
    fn roundtrip_preserves_the_effective_muzzle_angle() {
        let first = solve_v1(decode_solve_request_v1(&sample_json()).unwrap()).unwrap();
        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        assert_eq!(
            rebuilt.shot.muzzle_angle_rad,
            Some(first.resolved_request.shot.muzzle_angle_rad)
        );
        assert_eq!(
            rebuilt.shot.zero_distance_m,
            first.resolved_request.shot.zero_distance_m
        );
    }

    /// A QNH-declared pressure is reduced to absolute station pressure exactly once.
    /// `ResolvedAtmosphereV1::pressure_pa` is already that reduced value; the rebuilt
    /// request deliberately does not echo `pressure_reference: "qnh"` back alongside it
    /// (see the module doc), so the round-tripped resolved echo differs in exactly that one
    /// field. What must NOT differ is the reduced `pressure_pa` value itself: if the
    /// rebuilt request echoed the mode back too, the second resolve would reduce it a
    /// second time and silently corrupt it.
    #[test]
    fn roundtrip_does_not_double_reduce_a_qnh_pressure() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {"altitude_m": 500.0, "temperature_k": 288.0, "pressure_pa": 101325.0,
                           "pressure_reference": "qnh"},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let first = solve_v1(decode_solve_request_v1(&json).unwrap()).unwrap();
        assert_eq!(
            first.resolved_request.atmosphere.pressure_reference,
            Some(PressureReferenceV1::Qnh)
        );

        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        let second = solve_v1(rebuilt).unwrap();

        // The one deliberate difference: the reference mode is not echoed back (its
        // transform is already baked into pressure_pa, so re-supplying it would be a
        // second application, not a no-op).
        assert_eq!(second.resolved_request.atmosphere.pressure_reference, None);
        // What actually matters -- the physical quantity -- is unchanged.
        assert_eq!(
            second.resolved_request.atmosphere.pressure_pa,
            first.resolved_request.atmosphere.pressure_pa,
            "a round-tripped QNH pressure must not be reduced a second time"
        );
        assert_eq!(
            second.resolved_request.atmosphere.altitude_m,
            first.resolved_request.atmosphere.altitude_m
        );
        assert_eq!(
            second.resolved_request.atmosphere.temperature_k,
            first.resolved_request.atmosphere.temperature_k
        );
    }

    /// A compass-declared wind direction is converted to shooter-relative exactly once.
    /// The resolved direction is already that converted value; the rebuilt request
    /// deliberately does not echo `wind_reference: "compass"` back alongside it (see the
    /// module doc), so the round-tripped resolved echo differs in exactly that one field.
    /// What must NOT differ is the converted `direction_from_rad` value itself: if the
    /// rebuilt request echoed the mode back too, the second resolve would re-reference it
    /// against the shot azimuth a second time and silently corrupt it.
    #[test]
    fn roundtrip_does_not_double_reference_a_compass_wind() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "shot_azimuth_rad": 0.3},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": 1.0, "wind_reference": "compass"},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let first = solve_v1(decode_solve_request_v1(&json).unwrap()).unwrap();
        let ResolvedWindV1::Constant(first_wind) = &first.resolved_request.wind else {
            panic!("constant wind expected");
        };
        assert_eq!(first_wind.wind_reference, Some(WindReferenceV1::Compass));
        // Sanity check: compass mode actually converted the direction (it is not simply
        // echoing the 1.0 rad bearing the request supplied).
        assert_ne!(first_wind.direction_from_rad, 1.0);
        let first_direction_from_rad = first_wind.direction_from_rad;

        let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
        let second = solve_v1(rebuilt).unwrap();
        let ResolvedWindV1::Constant(second_wind) = &second.resolved_request.wind else {
            panic!("constant wind expected");
        };

        // The one deliberate difference: the reference mode is not echoed back.
        assert_eq!(second_wind.wind_reference, None);
        // What actually matters -- the physical quantity -- is unchanged.
        assert_eq!(
            second_wind.direction_from_rad, first_direction_from_rad,
            "a round-tripped compass wind must not be re-referenced a second time"
        );
    }
}
