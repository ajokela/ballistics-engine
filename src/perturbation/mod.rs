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
// trajectory_observation, all unconditional; solve_v1 itself depends on cli_api, likewise
// unconditional).
use crate::solve_json::{SolveErrorEnvelopeV1, SolveRequestV1};

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
/// Instead this follows the brief's fallback (option (b)), sharing `solve_v1`'s own building
/// blocks rather than re-implementing them: [`crate::solve_v1::prepare_request`] resolves `req`
/// exactly as `solve_v1` does, and [`crate::solve_v1::build_zeroed_solver`] builds and zeroes a
/// `TrajectorySolver` exactly as `solve_v1` does (both widened from module-private to
/// `pub(crate)` for this task). This function then calls `TrajectorySolver::solve` itself and
/// reads observations straight off the result: one solve per `evaluate` call, not two, and not
/// `solve_v1` internally.
///
/// (Revision note: an earlier version of this function hand-duplicated `solve_v1`'s
/// zero-handling `match` here instead of sharing `build_zeroed_solver`, kept honest only by an
/// outboard cross-check test. Review judged that duplication risk too costly to keep -- see
/// "Drop reference plane" below for a related place this function deliberately does NOT mirror
/// `solve_v1` -- so it was replaced with this structural share.)
///
/// # Drop reference plane
///
/// `evaluate` always reports `drop_m` perpendicular to the line of sight, regardless of the
/// request's `shot.drops_reference`. `solve_v1` applies the MBA-1403 target-plane transform
/// (`drop_m /= shooting_angle_rad.cos()`) as a wire-only rescaling at its own boundary
/// (`src/solve_v1.rs`, applied once to each sample after `TrajectoryResult::sample_observations`
/// runs); `evaluate` deliberately does not reproduce it, so the kernel always speaks one
/// geometry no matter which reference plane the original request asked for on output. This is a
/// decision, not an oversight -- see `evaluate_ignores_drops_reference_and_always_reports_los_drop`
/// in this module's tests, which pins it down. A caller that wants the target-plane number can
/// still derive it from `drop_m / shooting_angle_rad.cos()`, using the SAME `shooting_angle_rad`
/// this request resolved to (`ResolvedShotV1::shooting_angle_rad`) -- and should do that
/// explicitly rather than have it baked silently into `drop_m` here, because a later task
/// differentiates observations with respect to taxonomy axes including `ShootingAngle`
/// (`InputAxis::ShootingAngle`, `taxonomy.rs`): folding the target-plane transform into `drop_m`
/// in this function would silently add an extra `d/dtheta[sec(theta)]` term to that derivative,
/// and only for target-referenced requests.
///
/// # Errors
///
/// - [`KernelError::Solve`] if resolving or solving `req` fails -- the same failure modes
///   `solve_v1` itself reports (invalid or conflicting fields, a zero search that does not
///   converge, a non-finite effective muzzle angle), uniformly converted from
///   [`SolveErrorEnvelopeV1`] regardless of which stage (`prepare_request`,
///   `build_zeroed_solver`, or `TrajectorySolver::solve`) produced it.
/// - [`KernelError::Observation`] if any requested range cannot be read off the resulting
///   trajectory -- most notably a range outside `[0, actual_range]` -- via
///   [`crate::TrajectoryResult::observation_at_range_checked`], used specifically because it
///   ERRORS on an out-of-range query instead of clamping to the nearest sample (the pattern the
///   CLI card commands use, `src/main.rs:17734`).
pub fn evaluate(req: &SolveRequestV1, ranges_m: &[f64]) -> Result<Vec<Observation>, KernelError> {
    let prepared = crate::solve_v1::prepare_request(req).map_err(kernel_solve_error)?;

    let max_range_m = prepared.resolved_request.shot.max_range_m;
    let time_step_s = prepared.resolved_request.solver.time_step_s;
    let zero_distance_m = prepared.resolved_request.shot.zero_distance_m;
    let target_height_m = prepared.resolved_request.shot.target_height_m;

    let (solver, _effective_angle) = crate::solve_v1::build_zeroed_solver(
        prepared.inputs,
        prepared.wind,
        prepared.atmosphere,
        prepared.wind_segments,
        max_range_m,
        time_step_s,
        zero_distance_m,
        target_height_m,
        req.shot.muzzle_angle_rad,
    )
    .map_err(kernel_solve_error)?;

    let result = solver
        .solve()
        .map_err(crate::solve_v1::solve_failed)
        .map_err(kernel_solve_error)?;

    let mut out = Vec::with_capacity(ranges_m.len());
    for &range_m in ranges_m {
        let o = result
            .observation_at_range_checked(range_m)
            .map_err(KernelError::Observation)?;
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

/// Convert a v1 solve-error envelope into a `KernelError::Solve`, uniformly across every stage
/// `evaluate` can fail at (M3 review fix): `prepare_request`, `build_zeroed_solver`, and
/// `TrajectorySolver::solve` (via `solve_v1::solve_failed`) all go through this one conversion
/// instead of three independently hand-formatted strings. Carries `e.error.code` through
/// alongside the message (review fix I4(a)) rather than only formatting it into the string, so
/// `KernelError::is_domain_rejection` can tell an `InvalidValue` domain rejection apart from a
/// genuine solve failure without parsing text.
fn kernel_solve_error(e: SolveErrorEnvelopeV1) -> KernelError {
    KernelError::Solve { code: e.error.code, message: e.error.message }
}

// 0.33.0 decision-support Task 7: derived numerics over the kernel -- central-difference
// derivatives (for an error budget) and monotone bisection (for tolerance envelopes), built
// only on evaluate/read_axis/with_axis from Tasks 5-6. No feature gate: must compile for wasm32
// (same unconditional dependency chain as Task 6).
pub mod derive;
pub use derive::{
    bisect_axis, central_difference, DifferenceScheme, Derivative, BISECTION_MAX_ITERATIONS,
};

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
        // M1 review fix: `KernelError::Observation` can carry any `TrajectoryObservationError`
        // variant, so `matches!(.., Err(Observation(_)))` alone cannot distinguish an
        // out-of-range query from, say, `NonMonotonicTrajectory`. Match the specific `OutOfRange`
        // variant directly (I4(a) review fix: `Observation` now wraps the structured error
        // itself rather than a pre-formatted string, so this checks the actual variant instead
        // of substring-matching its rendered message).
        match evaluate(&req, &[5000.0]) {
            Err(KernelError::Observation(
                crate::trajectory_observation::TrajectoryObservationError::OutOfRange {
                    requested_m,
                    ..
                },
            )) => {
                assert_eq!(requested_m, 5000.0);
            }
            other => panic!(
                "expected Err(KernelError::Observation(OutOfRange {{ .. }})) naming an \
                 out-of-range query, got {other:?}"
            ),
        }
    }

    /// Independent oracle, `(Some(zero_distance_m), None)` arm (elevation search runs):
    /// `solve_v1`'s OWN wire samples, computed via a completely different call path
    /// (`TrajectoryResult::sample_observations` on a regular grid, not a targeted per-range
    /// query), must describe the identical physical point that `evaluate` reports for the same
    /// range on the same request. `solve_v1`'s mapping from `TrajectoryObservation` to
    /// `TrajectorySampleV1` is unchanged, existing, and load-bearing for the whole solve-json v1
    /// wire format, so it is a trustworthy check for `drop_m`/`windage_m` being swapped. See
    /// `evaluate_matches_solve_v1_when_an_explicit_angle_and_zero_distance_are_both_present` and
    /// `evaluate_matches_solve_v1_when_no_zero_distance_is_present_at_all` below for the other
    /// two zero-handling arms (I2 review fix -- this file previously tested only this one).
    #[test]
    fn evaluate_matches_solve_v1_when_zero_distance_alone_searches_the_elevation() {
        let json = base_request_json(900.0, 100.0, 300.0);
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        assert_eq!(
            req.shot.muzzle_angle_rad, None,
            "fixture assumption: no explicit angle, exercising the (Some, None) arm"
        );
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
        // be unambiguously nonzero, AND (M2 review fix) far enough apart from EACH OTHER that a
        // drop/windage swap could not hide behind two coincidentally-similar numbers -- the
        // original two independent-threshold checks alone could not catch that.
        assert!(sample_at_300.drop_m.abs() > 0.1);
        assert!(sample_at_300.windage_m.abs() > 0.01);
        assert!(
            (sample_at_300.drop_m - sample_at_300.windage_m).abs() > 0.05,
            "drop_m ({}) and windage_m ({}) must differ meaningfully, or a swap between them \
             would be invisible to the two independent-threshold checks above",
            sample_at_300.drop_m,
            sample_at_300.windage_m
        );

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

    /// Independent oracle, `(Some(zero_distance_m), Some(explicit_angle))` arm (I2 review fix):
    /// the original test suite only ever exercised the `(Some, None)` arm above --
    /// `base_request_json` always supplies `zero_distance_m`, and the one fixture omitting it
    /// fails inside `prepare_request` before the zero-handling match ever runs, so `(Some,
    /// Some)` and `(None, _)` were both silently untested despite the report claiming otherwise.
    /// This arm matters because it is exactly the shape `From<&ResolvedSolveRequestV1> for
    /// SolveRequestV1` (`request_roundtrip.rs`) produces on every zeroed solve: the elevation
    /// search does NOT run (the supplied angle is used directly); only the windage-zero bias
    /// (from `sight_offset_lateral_m`/`zero_poi_right_m`) applies. The fixture makes that bias
    /// non-negligible so a dropped `apply_windage_zero_bias` call would show up as a
    /// `windage_m` mismatch against `solve_v1`, not a coincidental match against zero.
    #[test]
    fn evaluate_matches_solve_v1_when_an_explicit_angle_and_zero_distance_are_both_present() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                      "sight_offset_lateral_m": 0.03},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0, "muzzle_angle_rad": 0.01,
                     "zero_poi_right_m": 0.02},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 300.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        assert!(
            req.shot.zero_distance_m.is_some() && req.shot.muzzle_angle_rad.is_some(),
            "fixture assumption: both fields present, exercising the (Some, Some) arm"
        );

        let via_solve_v1 = crate::solve_v1::solve_v1(req.clone()).expect("solve_v1");
        let sample_at_300 = via_solve_v1
            .samples
            .iter()
            .find(|s| s.distance_m == 300.0)
            .expect("300 m must be an exact grid point for this fixture");
        assert!(
            sample_at_300.windage_m.abs() > 0.005,
            "fixture must produce a non-negligible windage-zero bias, got {}",
            sample_at_300.windage_m
        );

        let obs = evaluate(&req, &[300.0]).expect("evaluate");
        assert_eq!(obs.len(), 1);
        assert_eq!(obs[0].drop_m, sample_at_300.drop_m);
        assert_eq!(
            obs[0].windage_m, sample_at_300.windage_m,
            "windage-zero bias diverged from solve_v1 -- apply_windage_zero_bias may not have \
             run for this arm"
        );
        assert_eq!(obs[0].time_s, sample_at_300.time_s);
        assert_eq!(obs[0].velocity_mps, sample_at_300.speed_mps);
    }

    /// Independent oracle, `(None, _)` arm (I2 review fix, continued): no `zero_distance_m` at
    /// all, with an explicit `muzzle_angle_rad` -- the elevation search does not run, and there
    /// is no `zero_distance_m` to gate a windage-zero bias on either, so the match falls through
    /// to its no-op arm entirely.
    #[test]
    fn evaluate_matches_solve_v1_when_no_zero_distance_is_present_at_all() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "muzzle_angle_rad": 0.01},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 300.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        assert_eq!(
            req.shot.zero_distance_m, None,
            "fixture assumption: no zero distance at all, exercising the (None, _) arm"
        );

        let via_solve_v1 = crate::solve_v1::solve_v1(req.clone()).expect("solve_v1");
        let sample_at_300 = via_solve_v1
            .samples
            .iter()
            .find(|s| s.distance_m == 300.0)
            .expect("300 m must be an exact grid point for this fixture");

        let obs = evaluate(&req, &[300.0]).expect("evaluate");
        assert_eq!(obs.len(), 1);
        assert_eq!(obs[0].drop_m, sample_at_300.drop_m);
        assert_eq!(obs[0].windage_m, sample_at_300.windage_m);
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

    /// `evaluate`'s zero-handling must actually RUN the zero search when `zero_distance_m` is
    /// present and no explicit `muzzle_angle_rad` was supplied -- not skip it, and not converge
    /// it against the wrong distance.
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
        // I4(a) review fix: `Solve` now carries the structured `SolveErrorCodeV1` alongside the
        // message; an out-of-range relative humidity is exactly the `InvalidValue` domain
        // rejection `require_range` produces (`solve_v1.rs`), so pin that down too rather than
        // only the outer variant.
        match evaluate(&req, &[100.0]) {
            Err(KernelError::Solve { code, .. }) => {
                assert_eq!(code, crate::solve_json::SolveErrorCodeV1::InvalidValue);
            }
            other => panic!("expected Err(KernelError::Solve {{ .. }}), got {other:?}"),
        }
    }

    /// I1 review fix: `evaluate` must ignore `shot.drops_reference` and always report the
    /// LOS-perpendicular `drop_m` (see `evaluate`'s doc comment, "Drop reference plane"), rather
    /// than silently reproducing (or half-reproducing) `solve_v1`'s MBA-1403 target-plane
    /// rescaling. Uses a 30-degree shooting angle so `cos(shooting_angle_rad)` (~0.866) is far
    /// enough from `1.0` that a target-plane rescaling would be an obvious ~15.5% disagreement,
    /// not noise.
    #[test]
    fn evaluate_ignores_drops_reference_and_always_reports_los_drop() {
        let shooting_angle_rad: f64 = 30.0_f64.to_radians();
        // An explicit muzzle_angle_rad rather than zero_distance_m: this test is about the
        // drops_reference transform, not zeroing, and a 30-degree incline combined with a
        // WorldVertical zero search at only 100 m is not guaranteed to bracket a convergent
        // angle (the trial trajectory is integrated AT the incline -- see
        // build_zeroed_solver/zero_trial_height_at -- so it is a materially different search
        // than a level zero). Supplying the angle directly sidesteps that entirely.
        let build = |drops_reference: Option<&str>| -> String {
            let mut shot = serde_json::json!({
                "max_range_m": 900.0,
                "muzzle_angle_rad": 0.02,
                "shooting_angle_rad": shooting_angle_rad,
            });
            if let Some(dr) = drops_reference {
                shot["drops_reference"] = serde_json::json!(dr);
            }
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": shot,
                "atmosphere": {}, "wind": {},
                "solver": {}, "effects": {}, "sampling": {"interval_m": 300.0}
            })
            .to_string()
        };

        let req_los = crate::solve_json::decode_solve_request_v1(&build(None)).unwrap();
        let req_target =
            crate::solve_json::decode_solve_request_v1(&build(Some("target"))).unwrap();

        let obs_los = evaluate(&req_los, &[300.0]).expect("evaluate los");
        let obs_target = evaluate(&req_target, &[300.0]).expect("evaluate target");
        assert_eq!(
            obs_los[0].drop_m, obs_target[0].drop_m,
            "evaluate must report the same LOS-perpendicular drop_m regardless of \
             shot.drops_reference"
        );

        // Cross-check against solve_v1's OWN target-plane transform: its wire sample under
        // "target" must equal evaluate's LOS drop_m divided by cos(shooting_angle_rad) --
        // confirming evaluate's number really is the untransformed LOS quantity, not an
        // accidental match.
        let via_solve_v1_target =
            crate::solve_v1::solve_v1(req_target.clone()).expect("solve_v1 target");
        let resolved_angle = via_solve_v1_target.resolved_request.shot.shooting_angle_rad;
        let sample = via_solve_v1_target
            .samples
            .iter()
            .find(|s| s.distance_m == 300.0)
            .expect("300 m must be an exact grid point for this fixture");
        let expected_target_drop = obs_los[0].drop_m / resolved_angle.cos();
        assert!(
            (sample.drop_m - expected_target_drop).abs() < 1e-9,
            "solve_v1's target-plane drop_m ({}) should equal evaluate's LOS drop_m ({}) \
             divided by cos(shooting_angle_rad) ({})",
            sample.drop_m,
            obs_los[0].drop_m,
            resolved_angle.cos()
        );
        // Sanity: the two planes actually differ measurably at this angle, or the test would
        // pass trivially even with a broken transform on either side.
        assert!(
            (sample.drop_m - obs_los[0].drop_m).abs() > 0.05,
            "LOS and target-plane drop must differ meaningfully at a 30 degree shooting angle"
        );
    }
}
