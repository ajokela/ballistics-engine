//! Derived numerics over the perturbation kernel: central-difference derivatives (feeding an
//! uncertainty/error budget) and monotone bisection (feeding tolerance envelopes).
//!
//! Both operations are built ONLY on the Task 5/6 primitives -- [`read_axis`], [`with_axis`],
//! and [`evaluate`] -- and inherit their semantics unchanged rather than reimplementing or
//! "correcting" anything:
//!
//! - `drop_m` stays LOS-perpendicular (see `evaluate`'s "Drop reference plane" doc comment in
//!   `mod.rs`): neither function in this file looks at `shot.drops_reference` at all.
//! - The specialized [`KernelError`] variants `with_axis` uses to refuse a physically-invalid
//!   axis/request combination -- `AxisUnsupportedForRequest` (`Altitude` under QNH pressure,
//!   `ShotAzimuth` under compass wind) and `AxisAbsent` (the three wind axes under segmented
//!   wind) -- propagate out of [`central_difference`] and [`bisect_axis`] unchanged via `?`.
//!   Neither function catches, retries, or maps them onto a derivative/bisection result: a
//!   caller must be able to tell "this axis cannot be perturbed on this request" apart from
//!   "the effect is zero" or "there is no crossing".
//!
//! # Step convention
//!
//! [`central_difference`] follows the one step-size heuristic the taxonomy already defines
//! (`axis_meta(axis).kind`'s `default_rel_step`/`min_abs_step`, `src/perturbation/taxonomy.rs`):
//! `h = (|x| * default_rel_step).max(min_abs_step)`, overridable by the caller's explicit
//! `step`. No second heuristic is introduced here.
//!
//! # Cost
//!
//! A central difference costs exactly two solves PER AXIS, not two solves per range: both
//! `evaluate` calls below take the whole `ranges_m` slice at once, so N ranges still cost 2
//! solves total (`evaluate` itself runs one `TrajectorySolver::solve` per call, however many
//! ranges are then read off the single result). For a `requires_rezero` axis (`taxonomy.rs`),
//! each of those two solves is itself preceded by up to 60 trial solves inside the elevation
//! search (`find_zero_angle`, `src/cli_api.rs`) -- unavoidable here, not something this task
//! changes. [`bisect_axis`] pays that same per-solve cost once per bisection iteration (up to
//! 80, capped the same way as the existing inverse-solver search,
//! `HoldCurve::range_for_angular_drop_mil` in `src/main.rs`), always at the single `range_m`
//! the caller asked to bisect at.
//!
//! # Bisection contract
//!
//! [`bisect_axis`] assumes `predicate` changes truth value at most once across `domain` (hence
//! "monotone" in this module's summary) and returns the crossing to within `tolerance`. When
//! `predicate` holds the SAME value at both ends of `domain`, there is nothing to bisect:
//! [`bisect_axis`] returns `Ok(None)` rather than inventing a bound, so callers can report
//! "stays inside throughout the configured domain" honestly instead of being handed a
//! fabricated one.

use crate::perturbation::access::{read_axis, with_axis, AxisValue, KernelError};
use crate::perturbation::taxonomy::{axis_meta, AxisKind, InputAxis};
use crate::perturbation::{evaluate, Observation};
use crate::solve_json::ResolvedSolveRequestV1;

/// First derivative of impact with respect to one axis, at one range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Derivative {
    pub axis: InputAxis,
    pub range_m: f64,
    pub d_drop_d_x: f64,
    pub d_windage_d_x: f64,
    pub step_used: f64,
}

/// Central difference: `(f(x+h) - f(x-h)) / 2h`, covering every range in `ranges_m` with
/// exactly two solves total (one per side) -- see the module doc's "Cost" section. The step
/// follows the crate convention `h = (|x| * default_rel_step).max(min_abs_step)` (`axis_meta`,
/// `taxonomy.rs`) unless the caller supplies an explicit `step`.
///
/// # Errors
///
/// - [`KernelError::CategoricalAxis`] if `axis_meta(axis).kind` is `AxisKind::Categorical`
///   (spec D7: categorical axes are never differentiated).
/// - [`KernelError::AxisAbsent`] if `axis` has no value on `base` -- `read_axis` returns `None`
///   here directly, without ever reaching `with_axis` -- most notably the three wind axes under
///   segmented wind.
/// - [`KernelError::TypeMismatch`] if `read_axis` ever returns a non-`Scalar` value for a
///   continuous axis. Not reachable with the current taxonomy (every `Continuous` axis reads
///   back as `Scalar` or `None`; `Flag`/`DragModel`/`TwistDirection` only ever come from
///   `Categorical` axes, already rejected above), kept as a defensive catch-all so a future
///   axis added under the wrong `AxisKind` fails loudly here instead of miscomputing silently.
/// - [`KernelError::NonFinite`] if the computed or caller-supplied step is not a finite
///   positive number, or if either resulting derivative is not finite.
/// - [`KernelError::AxisUnsupportedForRequest`], [`KernelError::Solve`], or
///   [`KernelError::Observation`], propagated unchanged from [`with_axis`]/[`evaluate`] on
///   either perturbed solve -- see the module doc.
pub fn central_difference(
    base: &ResolvedSolveRequestV1,
    axis: InputAxis,
    ranges_m: &[f64],
    step: Option<f64>,
) -> Result<Vec<Derivative>, KernelError> {
    let (rel, min_abs) = match axis_meta(axis).kind {
        AxisKind::Continuous { default_rel_step, min_abs_step, .. } => (default_rel_step, min_abs_step),
        AxisKind::Categorical => return Err(KernelError::CategoricalAxis(axis)),
    };
    let x = match read_axis(base, axis) {
        Some(AxisValue::Scalar(x)) => x,
        // Defensive only -- see the "not reachable" note in the doc comment above: every
        // Categorical axis (the only source of Flag/DragModel/TwistDirection) already returned
        // above.
        Some(_) => return Err(KernelError::TypeMismatch(axis)),
        None => return Err(KernelError::AxisAbsent(axis)),
    };
    let h = step.unwrap_or_else(|| (x.abs() * rel).max(min_abs));
    if !(h.is_finite() && h > 0.0) {
        return Err(KernelError::NonFinite(axis));
    }

    let plus = evaluate(&with_axis(base, axis, AxisValue::Scalar(x + h))?, ranges_m)?;
    let minus = evaluate(&with_axis(base, axis, AxisValue::Scalar(x - h))?, ranges_m)?;

    let mut out = Vec::with_capacity(ranges_m.len());
    for (p, m) in plus.iter().zip(minus.iter()) {
        let d_drop = (p.drop_m - m.drop_m) / (2.0 * h);
        let d_wind = (p.windage_m - m.windage_m) / (2.0 * h);
        if !d_drop.is_finite() || !d_wind.is_finite() {
            return Err(KernelError::NonFinite(axis));
        }
        out.push(Derivative {
            axis,
            range_m: p.range_m,
            d_drop_d_x: d_drop,
            d_windage_d_x: d_wind,
            step_used: h,
        });
    }
    Ok(out)
}

/// Bisect `axis` over `domain` at a single `range_m` until `predicate` (evaluated on the
/// observation at that range) changes truth value, to within `tolerance`. See the module doc's
/// "Bisection contract" for what `None` means and what `predicate` must satisfy.
///
/// # Errors
///
/// - [`KernelError::CategoricalAxis`] if `axis` is categorical.
/// - Any error [`with_axis`] or [`evaluate`] produce while probing `domain` -- including
///   [`KernelError::AxisUnsupportedForRequest`] and [`KernelError::AxisAbsent`] -- propagated
///   unchanged from the first probe that hits them (`domain.0`, checked before any bisection
///   step runs), so a request/axis combination `with_axis` refuses is reported as that specific
///   refusal, never silently folded into "no crossing" (`Ok(None)`).
pub fn bisect_axis(
    base: &ResolvedSolveRequestV1,
    axis: InputAxis,
    range_m: f64,
    domain: (f64, f64),
    predicate: &dyn Fn(&Observation) -> bool,
    tolerance: f64,
) -> Result<Option<f64>, KernelError> {
    if matches!(axis_meta(axis).kind, AxisKind::Categorical) {
        return Err(KernelError::CategoricalAxis(axis));
    }
    let at = |v: f64| -> Result<bool, KernelError> {
        let obs = evaluate(&with_axis(base, axis, AxisValue::Scalar(v))?, &[range_m])?;
        Ok(predicate(&obs[0]))
    };
    let (mut lo, mut hi) = domain;
    let lo_state = at(lo)?;
    if lo_state == at(hi)? {
        return Ok(None);
    }
    // 80 iterations matches the existing inverse-solver cap
    // (HoldCurve::range_for_angular_drop_mil's INVERSE_MAX_ITERATIONS, src/main.rs).
    for _ in 0..80 {
        if (hi - lo).abs() <= tolerance {
            break;
        }
        let mid = 0.5 * (lo + hi);
        if at(mid)? == lo_state {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(Some(0.5 * (lo + hi)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbation::InputAxis;

    fn resolved(mv: f64) -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": mv, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    /// More muzzle velocity means less drop at a fixed range: the derivative is negative.
    #[test]
    fn drop_derivative_wrt_muzzle_velocity_is_negative() {
        let r = resolved(823.0);
        let d = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], None).unwrap();
        assert_eq!(d.len(), 1);
        assert!(d[0].d_drop_d_x < 0.0, "expected negative, got {}", d[0].d_drop_d_x);
        assert!(d[0].step_used > 0.0);
    }

    /// Categorical axes must be refused, never silently differentiated (spec D7).
    #[test]
    fn categorical_axes_cannot_be_differentiated() {
        let r = resolved(823.0);
        let e = central_difference(&r, InputAxis::CoriolisEnabled, &[600.0], None);
        assert!(matches!(e, Err(KernelError::CategoricalAxis(_))));
    }

    /// Bisection finds the muzzle velocity at which drop crosses a chosen threshold.
    #[test]
    fn bisect_finds_the_crossing() {
        let r = resolved(823.0);
        let base = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], None).unwrap();
        let _ = base;
        let target_drop = 2.0_f64;
        let found = bisect_axis(&r, InputAxis::MuzzleVelocityMps, 600.0, (600.0, 1100.0),
                                &|o: &Observation| o.drop_m < target_drop, 0.05).unwrap();
        let mv = found.expect("a crossing exists in this domain");
        assert!(mv > 600.0 && mv < 1100.0);
    }

    /// Against the vacuum oracle the derivative has a closed form:
    /// drop = 0.5*g*(x/v)^2  =>  d(drop)/dv = -g*x^2/v^3.
    #[test]
    fn central_difference_matches_the_vacuum_analytic_derivative() {
        let json = serde_json::json!({
            "schema_version": 1,
            // bc_value huge => retardation negligible, same trick as the analytic_vacuum fuzzer
            "projectile": {"mass_kg": 0.01, "diameter_m": 0.0077, "drag_model": "G1",
                           "ballistic_coefficient": 100.0},
            "rifle": {"muzzle_velocity_mps": 800.0, "sight_height_m": 0.0},
            "shot": {"max_range_m": 500.0, "muzzle_angle_rad": 0.0},
            "atmosphere": {}, "wind": {}, "solver": {}, "effects": {},
            "sampling": {"interval_m": 5.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        let x = 400.0_f64;
        let v = 800.0_f64;
        let expected = -9.80665 * x * x / (v * v * v);
        let d = central_difference(&r, InputAxis::MuzzleVelocityMps, &[x], None).unwrap();
        let rel = ((d[0].d_drop_d_x - expected) / expected).abs();
        assert!(rel < 0.02, "expected ~{expected}, got {} (rel {rel})", d[0].d_drop_d_x);
    }

    /// None of the tests above ever pass more than one range to `central_difference`, even
    /// though the whole point of taking a slice (see the module doc's "Cost" section) is to
    /// cover every range from a SINGLE pair of solves. Extend the vacuum oracle across two
    /// ranges at once: both must match their own closed form, tagged with the CALLER's range
    /// (not the loop index or the other side's), and the longer range's derivative must be the
    /// larger one in magnitude (sensitivity to muzzle velocity grows with range^2) -- this would
    /// catch a bug that zipped `plus`/`minus` against the wrong range, or that reused `h` from
    /// one range's scale for another's.
    #[test]
    fn central_difference_matches_the_vacuum_oracle_at_every_requested_range() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.01, "diameter_m": 0.0077, "drag_model": "G1",
                           "ballistic_coefficient": 100.0},
            "rifle": {"muzzle_velocity_mps": 800.0, "sight_height_m": 0.0},
            "shot": {"max_range_m": 500.0, "muzzle_angle_rad": 0.0},
            "atmosphere": {}, "wind": {}, "solver": {}, "effects": {},
            "sampling": {"interval_m": 5.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        let v = 800.0_f64;
        let ranges = [200.0_f64, 400.0_f64];
        let d = central_difference(&r, InputAxis::MuzzleVelocityMps, &ranges, None).unwrap();
        assert_eq!(d.len(), 2);
        for (i, &x) in ranges.iter().enumerate() {
            assert_eq!(
                d[i].range_m, x,
                "derivative {i} must be tagged with the range it was requested at"
            );
            let expected = -9.80665 * x * x / (v * v * v);
            let rel = ((d[i].d_drop_d_x - expected) / expected).abs();
            assert!(
                rel < 0.02,
                "range {x}: expected ~{expected}, got {} (rel {rel})",
                d[i].d_drop_d_x
            );
        }
        assert!(
            d[1].d_drop_d_x.abs() > d[0].d_drop_d_x.abs() * 3.0,
            "sensitivity to muzzle velocity should grow with range^2: at {} m got {}, at {} m got {}",
            ranges[0], d[0].d_drop_d_x, ranges[1], d[1].d_drop_d_x
        );
    }

    /// The tests above only ever check `d_drop_d_x`; `d_windage_d_x` is never exercised, so a
    /// copy-paste bug that computed it from `drop_m` instead of `windage_m` (or one that always
    /// left it at 0) would pass every other test in this file. Wind speed pins it down: the
    /// `resolved` fixture's wind is a 90-degree (full) crosswind (see `mod.rs`'s
    /// `base_request_json` doc for the same convention), so more wind speed must push windage
    /// further in the SAME direction as the baseline windage, and that push must dominate
    /// whatever tiny secondary effect wind speed has on drop.
    #[test]
    fn windage_derivative_wrt_wind_speed_dominates_and_matches_the_baseline_sign() {
        let r = resolved(823.0);
        let baseline_req: crate::solve_json::SolveRequestV1 = (&r).into();
        let baseline = evaluate(&baseline_req, &[600.0]).expect("baseline evaluate");
        assert!(
            baseline[0].windage_m.abs() > 0.01,
            "fixture must have non-negligible baseline windage, got {}",
            baseline[0].windage_m
        );

        let d = central_difference(&r, InputAxis::WindSpeed, &[600.0], None).unwrap();
        assert_eq!(d.len(), 1);
        assert_eq!(
            d[0].d_windage_d_x.signum(),
            baseline[0].windage_m.signum(),
            "more crosswind should push windage further the SAME way: derivative {}, baseline {}",
            d[0].d_windage_d_x,
            baseline[0].windage_m
        );
        assert!(
            d[0].d_windage_d_x.abs() > d[0].d_drop_d_x.abs() * 5.0,
            "a pure crosswind should move windage far more than it moves drop: \
             d_windage_d_x={}, d_drop_d_x={}",
            d[0].d_windage_d_x,
            d[0].d_drop_d_x
        );
    }

    /// `step_used` is part of the public contract of a `Derivative` -- a caller building an
    /// error budget reads it directly to scale an input uncertainty. Pin the formula down
    /// independently of any particular axis's physical response, so a bug that reports a
    /// different `h` than the one actually used in the `2h` division cannot hide behind the
    /// vacuum oracle's tolerance.
    #[test]
    fn step_used_follows_the_crate_convention() {
        let r = resolved(823.0);
        let x = match read_axis(&r, InputAxis::WindSpeed).unwrap() {
            AxisValue::Scalar(x) => x,
            other => panic!("WindSpeed must read back as a scalar, got {other:?}"),
        };
        let (rel, min_abs) = match axis_meta(InputAxis::WindSpeed).kind {
            AxisKind::Continuous { default_rel_step, min_abs_step, .. } => {
                (default_rel_step, min_abs_step)
            }
            AxisKind::Categorical => panic!("WindSpeed must be continuous"),
        };
        let expected_h = (x.abs() * rel).max(min_abs);
        let d = central_difference(&r, InputAxis::WindSpeed, &[600.0], None).unwrap();
        assert_eq!(d[0].step_used, expected_h);
    }

    /// If `predicate` never flips across `domain`, `bisect_axis` must say so with `None` -- not
    /// fabricate a bound by returning some point inside the domain anyway (the exact failure
    /// mode this function's doc comment calls out). A target so large that drop can never reach
    /// it over a 600 m shot holds `true` at both ends of the domain, so there is nothing to
    /// bisect.
    #[test]
    fn bisect_axis_returns_none_when_the_predicate_never_flips() {
        let r = resolved(823.0);
        let found = bisect_axis(
            &r,
            InputAxis::MuzzleVelocityMps,
            600.0,
            (600.0, 1100.0),
            &|o: &Observation| o.drop_m < 1000.0,
            0.05,
        )
        .unwrap();
        assert!(
            found.is_none(),
            "predicate holds everywhere on this domain; bisect_axis must report None, not a \
             fabricated crossing, got {found:?}"
        );
    }

    /// Independent closed-form check for `bisect_axis` itself (not just `central_difference`):
    /// reuse the near-vacuum trick (huge BC => drag negligible) so `drop = 0.5*g*(x/v)^2` is
    /// invertible in closed form: `v = x*sqrt(g/(2*drop))`. The domain `(500, 2000)` is
    /// deliberately asymmetric around that root (~886 m/s, well off the domain's own midpoint
    /// of 1250), so a bug that returned the domain midpoint instead of converging would be off
    /// by roughly 41% -- nowhere close to the 2% tolerance this file uses elsewhere for the
    /// same physical approximation. No `zero_distance_m` is present (only an explicit
    /// `muzzle_angle_rad`), so -- like the derivative oracle above -- no re-zero search runs and
    /// `MuzzleVelocityMps` bisects cheaply here despite being a `requires_rezero` axis in
    /// general.
    #[test]
    fn bisect_axis_converges_to_the_vacuum_analytic_root_not_the_domain_midpoint() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.01, "diameter_m": 0.0077, "drag_model": "G1",
                           "ballistic_coefficient": 100.0},
            "rifle": {"muzzle_velocity_mps": 800.0, "sight_height_m": 0.0},
            "shot": {"max_range_m": 500.0, "muzzle_angle_rad": 0.0},
            "atmosphere": {}, "wind": {}, "solver": {}, "effects": {},
            "sampling": {"interval_m": 5.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let x = 400.0_f64;
        let target_drop = 1.0_f64;
        let domain = (500.0_f64, 2000.0_f64);
        let expected_v = x * (9.80665_f64 / (2.0 * target_drop)).sqrt();
        // Sanity check on the fixture itself: the analytic root must not be suspiciously close
        // to the domain's midpoint, or this test would not actually distinguish "converged" from
        // "returned the midpoint".
        let midpoint = 0.5 * (domain.0 + domain.1);
        assert!(
            (expected_v - midpoint).abs() / expected_v > 0.2,
            "fixture must keep the analytic root well away from the domain midpoint: \
             root {expected_v}, midpoint {midpoint}"
        );

        let found = bisect_axis(
            &r,
            InputAxis::MuzzleVelocityMps,
            x,
            domain,
            &|o: &Observation| o.drop_m < target_drop,
            0.05,
        )
        .unwrap()
        .expect("a crossing exists in this domain");

        let rel = ((found - expected_v) / expected_v).abs();
        assert!(
            rel < 0.02,
            "expected the bisection to converge near the analytic root ~{expected_v}, got \
             {found} (rel {rel})"
        );
    }

    /// `with_axis` refuses `Altitude` on a QNH-referenced request (`access.rs`) because the
    /// rebuilt request cannot re-derive the original altimeter setting. `central_difference`
    /// must PROPAGATE that refusal, not swallow it into a derivative of zero -- a caller needs
    /// to tell "this axis cannot be perturbed here" apart from "the effect is zero".
    #[test]
    fn central_difference_propagates_axis_unsupported_for_request() {
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
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let e = central_difference(&r, InputAxis::Altitude, &[300.0], None);
        match e {
            Err(KernelError::AxisUnsupportedForRequest { axis: InputAxis::Altitude, reason }) => {
                assert!(reason.to_lowercase().contains("qnh"), "reason should name QNH: {reason}");
            }
            other => panic!("expected AxisUnsupportedForRequest, got {other:?}"),
        }
    }

    /// Same guarantee via `bisect_axis`: the refusal must surface on the very first domain
    /// probe, not be absorbed into `Ok(None)`.
    #[test]
    fn bisect_axis_propagates_axis_unsupported_for_request() {
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
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let e = bisect_axis(
            &r,
            InputAxis::Altitude,
            300.0,
            (400.0, 600.0),
            &|o: &Observation| o.drop_m < 1.0,
            0.5,
        );
        match e {
            Err(KernelError::AxisUnsupportedForRequest { axis: InputAxis::Altitude, reason }) => {
                assert!(reason.to_lowercase().contains("qnh"), "reason should name QNH: {reason}");
            }
            other => panic!("expected AxisUnsupportedForRequest, got {other:?}"),
        }
    }

    /// `read_axis` returns `None` for the three wind axes under segmented wind (no single
    /// scalar to perturb, taxonomy.rs Known Limitation (c)). `central_difference` must turn
    /// that into `KernelError::AxisAbsent` via its OWN `None` branch -- this is a different code
    /// path than `with_axis`'s segmented-wind guard, so it needs its own test.
    #[test]
    fn central_difference_reports_axis_absent_for_wind_axes_under_segmented_wind() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {},
            "wind": {"segments": [{"until_distance_m": 900.0, "speed_mps": 3.0,
                                    "direction_from_rad": 1.0}]},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let e = central_difference(&r, InputAxis::WindSpeed, &[300.0], None);
        assert!(matches!(e, Err(KernelError::AxisAbsent(InputAxis::WindSpeed))));
    }

    /// Same guarantee via `bisect_axis`, which never calls `read_axis` at all -- its
    /// `AxisAbsent` must come through `with_axis`'s own segmented-wind guard instead, exercised
    /// on the very first probe.
    #[test]
    fn bisect_axis_reports_axis_absent_for_wind_axes_under_segmented_wind() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {},
            "wind": {"segments": [{"until_distance_m": 900.0, "speed_mps": 3.0,
                                    "direction_from_rad": 1.0}]},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let e = bisect_axis(
            &r,
            InputAxis::WindSpeed,
            300.0,
            (0.0, 10.0),
            &|o: &Observation| o.windage_m < 0.0,
            0.5,
        );
        assert!(matches!(e, Err(KernelError::AxisAbsent(InputAxis::WindSpeed))));
    }
}
