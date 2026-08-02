//! Request-level perturbation kernel (Phase 1 of the 0.33.0 decision-support train).
pub mod taxonomy;
pub use taxonomy::{axes_in_group, axis_meta, AxisKind, AxisMeta, InputAxis, InputGroup};

// 0.33.0 decision-support Task 5: read/write access to a single taxonomy axis on a canonical
// request. No feature gate: must compile for wasm32 (depends only on solve_json/solve_v1,
// both unconditional).
pub mod access;
pub use access::{read_axis, with_axis, AxisValue, KernelError};

// 0.33.0 decision-support Task 6: the OUTPUT side of the kernel -- solve a request once and
// read checked observations (drop, windage, time, velocity) at caller-selected ranges. No
// feature gate: must compile for wasm32 (depends only on solve_json/solve_v1/
// trajectory_observation/cli_api, all unconditional).
use crate::cli_api::ZeroTargetFrame;
use crate::solve_json::SolveRequestV1;
use crate::TrajectorySolver;

/// One solved observation at one requested range.
///
/// SI throughout. `drop_m` is positive BELOW the line of sight and `windage_m` is positive to
/// the RIGHT of the line of sight from the shooter's perspective -- the same convention as
/// [`crate::TrajectoryResult::observation_at_range_checked`]'s
/// [`crate::trajectory_observation::TrajectoryObservation`], whose `drop_m`/`windage_m` are
/// copied here unchanged, and as every other consumer of that type in this crate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    pub range_m: f64,
    pub drop_m: f64,
    pub windage_m: f64,
    pub time_s: f64,
    pub velocity_mps: f64,
}

/// Solve `req` once and read a checked observation at each of `ranges_m`, in the order given.
///
/// # Which `TrajectoryResult` path this uses, and why
///
/// This task's brief assumed a `SolveSuccessV1::trajectory_result_for_observation()` accessor;
/// no such method exists. [`crate::solve_json::SolveSuccessV1`] does not expose a raw
/// [`crate::TrajectoryResult`] at all -- only the wire-shaped, already-downsampled
/// `samples: Vec<TrajectorySampleV1>` (the brief's option (a), ruled out for that reason).
///
/// Instead this follows the brief's fallback (option (b)):
/// [`crate::solve_v1::prepare_request`] (widened from module-private to `pub(crate)` for this
/// task) resolves `req` into the same `BallisticInputs`/`WindConditions`/
/// `AtmosphericConditions` triple that [`crate::solve_v1::solve_v1`] itself builds a solver
/// from, and this function drives that [`TrajectorySolver`] itself: one solve, not two, and
/// not a second, independent entry point that could silently drift from what `solve_v1` would
/// report for the same request. The zero-handling `match` below is a deliberate line-for-line
/// mirror of the one in `solve_v1` (`src/solve_v1.rs`); both sides carry a cross-reference
/// comment so a future edit to one is not silently missed in the other.
///
/// # Errors
///
/// - [`KernelError::Solve`] if resolving or solving `req` fails -- the same failure modes
///   `solve_v1` itself reports (invalid or conflicting fields, a zero search that does not
///   converge, a non-finite effective muzzle angle).
/// - [`KernelError::Observation`] if any requested range cannot be read off the resulting
///   trajectory -- most notably a range outside `[0, actual_range]` -- via
///   [`crate::TrajectoryResult::observation_at_range_checked`], used specifically because it
///   ERRORS on an out-of-range query instead of clamping to the nearest sample (the pattern the
///   CLI card commands use, `src/main.rs:17734`).
pub fn evaluate(req: &SolveRequestV1, ranges_m: &[f64]) -> Result<Vec<Observation>, KernelError> {
    let prepared = crate::solve_v1::prepare_request(req)
        .map_err(|e| KernelError::Solve(format!("{:?}: {}", e.error.code, e.error.message)))?;

    let max_range_m = prepared.resolved_request.shot.max_range_m;
    let time_step_s = prepared.resolved_request.solver.time_step_s;
    let zero_distance_m = prepared.resolved_request.shot.zero_distance_m;
    let target_height_m = prepared.resolved_request.shot.target_height_m;

    let mut solver = TrajectorySolver::new_with_resolved_station_atmosphere(
        prepared.inputs,
        prepared.wind,
        prepared.atmosphere,
    );
    solver.set_max_range(max_range_m);
    solver.set_time_step(time_step_s);
    if !prepared.wind_segments.is_empty() {
        solver.set_wind_segments(prepared.wind_segments);
    }

    // Mirrors solve_v1's zero-handling match exactly (src/solve_v1.rs -- see the
    // cross-reference comment there): an explicit muzzle_angle_rad on the ORIGINAL
    // (unresolved) request always wins over zero_distance_m for elevation, but a present
    // zero_distance_m always contributes its windage-zero bias regardless of whether the
    // elevation search ran. If that match ever changes, this one must change with it.
    match (zero_distance_m, req.shot.muzzle_angle_rad) {
        (Some(distance_m), None) => {
            let effective_angle = solver
                .calculate_and_set_zero_angle(
                    distance_m,
                    target_height_m,
                    ZeroTargetFrame::WorldVertical,
                )
                .map_err(|e| KernelError::Solve(e.to_string()))?;
            if !effective_angle.is_finite() {
                return Err(KernelError::Solve(
                    "zero-angle calculation returned a non-finite effective muzzle angle"
                        .to_string(),
                ));
            }
        }
        (Some(distance_m), Some(_)) => solver.apply_windage_zero_bias(distance_m),
        (None, _) => {}
    }

    let result = solver
        .solve()
        .map_err(|e| KernelError::Solve(e.to_string()))?;

    let mut out = Vec::with_capacity(ranges_m.len());
    for &range_m in ranges_m {
        let o = result
            .observation_at_range_checked(range_m)
            .map_err(|e| KernelError::Observation(e.to_string()))?;
        out.push(Observation {
            range_m,
            drop_m: o.drop_m,
            windage_m: o.windage_m,
            time_s: o.time_s,
            velocity_mps: o.speed_mps,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod eval_tests {
    use super::*;

    /// Shared fixture: a G7 .308-class projectile, zeroed at `zero_distance_m`, with a
    /// 90-degree (full) crosswind so wind-drift-dependent assertions have something nonzero to
    /// check.
    fn base_request_json(max_range_m: f64, zero_distance_m: f64, interval_m: f64) -> String {
        serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": max_range_m, "zero_distance_m": zero_distance_m},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": interval_m}
        })
        .to_string()
    }

    #[test]
    fn evaluate_returns_one_observation_per_requested_range() {
        let req =
            crate::solve_json::decode_solve_request_v1(&base_request_json(900.0, 100.0, 25.0))
                .unwrap();
        let obs = evaluate(&req, &[300.0, 600.0, 800.0]).expect("evaluate");
        assert_eq!(obs.len(), 3);
        assert!((obs[0].range_m - 300.0).abs() < 1e-6);
        assert!((obs[1].range_m - 600.0).abs() < 1e-6);
        assert!((obs[2].range_m - 800.0).abs() < 1e-6);
        // Drop grows with range (positive = below line of sight) once past the zero.
        assert!(obs[2].drop_m > obs[0].drop_m);
        // Velocity decays.
        assert!(obs[2].velocity_mps < obs[0].velocity_mps);
        assert!(obs
            .iter()
            .all(|o| o.time_s.is_finite() && o.windage_m.is_finite()));
    }

    #[test]
    fn a_range_beyond_the_trajectory_is_an_error_not_a_clamp() {
        let req =
            crate::solve_json::decode_solve_request_v1(&base_request_json(300.0, 100.0, 25.0))
                .unwrap();
        assert!(matches!(
            evaluate(&req, &[5000.0]),
            Err(KernelError::Observation(_))
        ));
    }

    /// Independent oracle: `solve_v1`'s OWN wire samples, computed via a completely different
    /// call path (`TrajectoryResult::sample_observations` on a regular grid, not a targeted
    /// per-range query), must describe the identical physical point that `evaluate` reports for
    /// the same range on the same request. `solve_v1`'s mapping from `TrajectoryObservation` to
    /// `TrajectorySampleV1` is unchanged, existing, and load-bearing for the whole solve-json v1
    /// wire format, so it is a trustworthy check for `drop_m`/`windage_m` being swapped, and for
    /// this module's mirrored zero-handling silently diverging from solve_v1's original.
    #[test]
    fn evaluate_matches_solve_v1s_own_reported_sample_at_the_same_distance() {
        let json = base_request_json(900.0, 100.0, 300.0);
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let via_solve_v1 = crate::solve_v1::solve_v1(req.clone()).expect("solve_v1");
        // 300.0 is an exact regular-grid point for interval_m=300 starting at the muzzle
        // (x=0.0): grid index 1, strictly before the ~900 m terminal, so it is read straight
        // off the trajectory rather than needing end-of-grid special-casing.
        let sample_at_300 = via_solve_v1
            .samples
            .iter()
            .find(|s| s.distance_m == 300.0)
            .expect("300 m must be an exact grid point for this fixture");
        // Sanity check that this fixture actually distinguishes the two fields: with a
        // 90-degree crosswind and a zero well short of 300 m, drop and windage here must both
        // be unambiguously nonzero (and different in magnitude), or a drop/windage swap could
        // hide behind two coincidentally-similar numbers.
        assert!(sample_at_300.drop_m.abs() > 0.1);
        assert!(sample_at_300.windage_m.abs() > 0.01);

        let obs = evaluate(&req, &[300.0]).expect("evaluate");
        assert_eq!(obs.len(), 1);
        assert_eq!(
            obs[0].drop_m, sample_at_300.drop_m,
            "drop_m diverged from solve_v1's own reported sample"
        );
        assert_eq!(
            obs[0].windage_m, sample_at_300.windage_m,
            "windage_m diverged from solve_v1's own reported sample (possible field swap)"
        );
        assert_eq!(obs[0].time_s, sample_at_300.time_s);
        assert_eq!(obs[0].velocity_mps, sample_at_300.speed_mps);
    }

    /// `evaluate` must return observations in the CALLER's requested order, not sorted, and
    /// must not deduplicate a repeated range. Each output element is tied to its TRUE range via
    /// an independent physical fact (velocity strictly decreases with range for a supersonic
    /// flat-fire rifle round over this bracket) rather than only checking that `range_m` echoes
    /// the query, which a bug that transposes a value together with its label would still pass.
    #[test]
    fn evaluate_preserves_the_caller_supplied_range_order_and_repeats() {
        let req =
            crate::solve_json::decode_solve_request_v1(&base_request_json(900.0, 100.0, 25.0))
                .unwrap();
        let requested = [800.0, 300.0, 800.0, 500.0];
        let obs = evaluate(&req, &requested).expect("evaluate");
        assert_eq!(obs.len(), requested.len());
        for (o, &want) in obs.iter().zip(requested.iter()) {
            assert_eq!(o.range_m, want);
        }
        assert!(
            obs[0].velocity_mps < obs[1].velocity_mps,
            "800 m must be slower than 300 m"
        );
        assert_eq!(
            obs[0].velocity_mps, obs[2].velocity_mps,
            "two queries at the same 800 m range must agree exactly"
        );
        assert!(
            obs[1].velocity_mps > obs[3].velocity_mps,
            "300 m must be faster than 500 m"
        );
    }

    /// `evaluate`'s mirrored zero-handling must actually RUN the zero search when
    /// `zero_distance_m` is present and no explicit `muzzle_angle_rad` was supplied -- not skip
    /// it, and not converge it against the wrong distance.
    ///
    /// Independent physical check: the zero search converges the bullet's absolute
    /// world-vertical height at `zero_distance_m` to `shot.target_height_m` (default `0`; see
    /// `docs/SOLVE_JSON_V1.md`'s `target_height_m` row -- it is a height above the ground
    /// datum, NOT automatically the sight line). `drop_m` is `line_of_sight_height_m -
    /// vertical_m`, and `line_of_sight_height_m` is `muzzle_height_m + sight_height_m`. This
    /// fixture sets `sight_height_m: 0.0` (no bore/sight offset) specifically so
    /// `line_of_sight_height_m == 0 == target_height_m`: with the offset eliminated, a
    /// converged zero and a near-zero `drop_m` at that range are the same fact, so this
    /// isolates "did the elevation search run and converge at the right distance" without
    /// dragging in the separate sight-height-vs-target-height relationship (a nonzero
    /// `sight_height_m` alone, with `target_height_m` still `0`, leaves a `drop_m` residual of
    /// approximately `sight_height_m` at the zero range by design -- that is not this test).
    /// `find_zero_angle` (`src/cli_api.rs`) converges once the height error is under 1e-4 m.
    #[test]
    fn evaluate_actually_zeroes_drop_is_near_zero_at_the_zero_distance() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.0},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 137.0},
            "atmosphere": {}, "wind": {},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let obs = evaluate(&req, &[137.0]).expect("evaluate");
        assert_eq!(obs.len(), 1);
        assert!(
            obs[0].drop_m.abs() < 1e-3,
            "a rifle zeroed at 137 m must show ~0 drop there, got {}",
            obs[0].drop_m
        );
    }

    /// A request that fails at the RESOLVE stage (not the observation-lookup stage) must
    /// surface as `KernelError::Solve`, exercising the `prepare_request` error path rather than
    /// only `observation_at_range_checked`'s. `decode_solve_request_v1` does not itself bound
    /// `relative_humidity` (only mass/diameter/length/BC are checked at decode time -- see
    /// `validate_request_ranges`), so this fixture reaches `evaluate` and fails only once
    /// `prepare_request` -> `resolve_atmosphere` validates the physical range.
    #[test]
    fn an_invalid_request_is_a_solve_error_not_a_panic_or_silent_success() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 300.0},
            "atmosphere": {"relative_humidity": 1.5},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        assert!(matches!(
            evaluate(&req, &[100.0]),
            Err(KernelError::Solve(_))
        ));
    }
}
