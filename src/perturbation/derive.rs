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
//! A central difference costs exactly two solves PER AXIS in the common case, not two solves
//! per range: both `evaluate` calls below take the whole `ranges_m` slice at once, so N ranges
//! still cost 2 solves total (`evaluate` itself runs one `TrajectorySolver::solve` per call,
//! however many ranges are then read off the single result). When the one-sided fallback below
//! fires, a THIRD solve (of the unperturbed request, at `x` itself) is needed -- still
//! independent of how many ranges are requested. For a `requires_rezero` axis (`taxonomy.rs`),
//! each solve is itself preceded by up to 60 trial solves inside the elevation search
//! (`find_zero_angle`, `src/cli_api.rs`) -- unavoidable here, not something this task changes.
//! [`bisect_axis`] pays that same per-solve cost once per bisection iteration (up to
//! [`BISECTION_MAX_ITERATIONS`], capped the same way as the existing inverse-solver search,
//! `HoldCurve::range_for_angular_drop_mil` in `src/main.rs`), always at the single `range_m`
//! the caller asked to bisect at.
//!
//! # One-sided fallback
//!
//! Several continuous axes have a hard physical domain narrower than all reals: `WindSpeed` is
//! a non-negative magnitude, `RelativeHumidity` is confined to `[0, 1]`, and `TargetDistance`
//! (== `shot.max_range_m`) cannot shrink below a range the caller is asking about. A central
//! difference at or near such a boundary needs a perturbed value on the wrong side of it --
//! still air (`speed_mps: 0.0`, the default absent any `wind` block at all) is the most
//! ordinary example, not an exotic one: `WindSpeed`'s `min_abs_step` (0.05 m/s) makes the minus
//! side `-0.05`, which `resolve_wind`'s `require_non_negative("$.wind.speed_mps")`
//! (`src/solve_v1.rs`) rejects outright.
//!
//! When exactly one of the two perturbed solves fails to evaluate and the other succeeds,
//! [`central_difference`] falls back to a one-sided difference using the side that worked plus
//! the UNPERTURBED value at `x` itself: `(f(x+h) - f(x)) / h`
//! ([`DifferenceScheme::ForwardOneSided`]) if the minus side failed, or
//! `(f(x) - f(x-h)) / h` ([`DifferenceScheme::BackwardOneSided`]) if the plus side failed. This
//! is not merely an accommodation for the solver rejecting an unphysical input: at a hard
//! domain edge (wind speed pinned at exactly zero) a symmetric difference is not merely
//! blocked, it would be answering the wrong question, since windage as a function of SIGNED
//! wind speed is not smooth across that boundary (it is V-shaped, not linear) in the first
//! place.
//!
//! Which scheme actually ran is part of [`Derivative`]'s public contract (its `scheme` field):
//! a one-sided difference has different (generally larger, `O(h)` rather than `O(h^2)`)
//! truncation error than a central one, and a caller building an error budget or reporting a
//! method must be able to tell them apart rather than silently trusting every [`Derivative`] as
//! if it were central.
//!
//! Only VALUE-dependent failures trigger this fallback -- specifically, a failure from
//! [`evaluate`] on an already-successfully-built perturbed request (a physically invalid
//! quantity, or a query range that fell outside a shrunk `max_range_m`). Failures from
//! [`with_axis`] itself (`AxisUnsupportedForRequest`, `AxisAbsent`, `TypeMismatch`) do not
//! depend on the perturbed value at all -- they would occur identically for `x+h` and `x-h`,
//! since they are checked from `axis` and `base`'s OTHER fields before the value is even
//! considered -- so they continue to propagate immediately with no fallback attempted, exactly
//! as before this fix.
//!
//! If BOTH perturbed sides fail to evaluate, there is no data left to build even a one-sided
//! difference from: [`central_difference`] returns [`KernelError::StepOutOfDomain`] rather than
//! propagating whichever of the two solve errors happened to be checked first (arbitrary, and
//! neither one alone explains the real problem, which is the step, not either individual
//! solve).
//!
//! # Bisection contract
//!
//! [`bisect_axis`] assumes `predicate` changes truth value at most once across `domain` (hence
//! "monotone" in this module's summary) and returns the crossing to within `tolerance`.
//!
//! `Ok(None)` means ONLY that `predicate` did not change truth value across `domain` -- nothing
//! more. It does NOT tell the caller which of two opposite facts that is: `predicate` could be
//! true at both ends (e.g. "stays inside a tolerance band throughout this domain") or FALSE at
//! both ends (e.g. "stays outside it throughout"), and `Ok(None)` looks identical either way. A
//! caller that needs to know which one happened must check `predicate` at an endpoint itself
//! (or already know the end state some other way) -- `bisect_axis` deliberately does not
//! resolve that ambiguity. This matters most for a naturally two-sided predicate like
//! `|drop - nominal| <= tolerance`, which is true in the middle and false at BOTH ends:
//! widening `domain` to find the edge of such a band and getting `Ok(None)` back means only "no
//! edge in this domain," never "true throughout" -- treating it as the latter would be exactly
//! the fabricated-bound failure this function exists to avoid, just relocated into the caller's
//! interpretation of a correct result. Contrast `HoldCurve::range_for_angular_drop_mil`
//! (`src/main.rs`), which instead reports each out-of-domain case as its own distinct outcome;
//! `bisect_axis` does not do that here, so its caller must.

use crate::perturbation::access::{read_axis, with_axis, AxisValue, KernelError};
use crate::perturbation::taxonomy::{axis_meta, AxisKind, InputAxis};
use crate::perturbation::{evaluate, Observation};
use crate::solve_json::ResolvedSolveRequestV1;

/// Which finite-difference formula actually produced a [`Derivative`] (review fix I4).
///
/// A one-sided scheme has different, generally larger (`O(h)` vs. `O(h^2)`) truncation error
/// than a central one -- see the module doc's "One-sided fallback" section for when and why
/// each one is chosen. This is part of the public contract precisely so a caller building an
/// error budget (or reporting a method, MBA-1347) does not have to silently assume every
/// [`Derivative`] is central.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DifferenceScheme {
    /// `(f(x+h) - f(x-h)) / 2h` -- both perturbed sides evaluated successfully.
    Central,
    /// `(f(x+h) - f(x)) / h` -- the backward side (`x-h`) failed to evaluate (most often
    /// because it left the axis's physical domain), so the forward side and the unperturbed
    /// value at `x` were used instead.
    ForwardOneSided,
    /// `(f(x) - f(x-h)) / h` -- the forward side (`x+h`) failed to evaluate (most often because
    /// it left the axis's physical domain), so the backward side and the unperturbed value at
    /// `x` were used instead.
    BackwardOneSided,
}

/// First derivative of impact with respect to one axis, at one range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Derivative {
    pub axis: InputAxis,
    pub range_m: f64,
    pub d_drop_d_x: f64,
    pub d_windage_d_x: f64,
    pub step_used: f64,
    pub scheme: DifferenceScheme,
}

/// Central difference: `(f(x+h) - f(x-h)) / 2h`, covering every range in `ranges_m` with
/// exactly two solves total (one per side) in the common case -- see the module doc's "Cost"
/// section. The step follows the crate convention
/// `h = (|x| * default_rel_step).max(min_abs_step)` (`axis_meta`, `taxonomy.rs`) unless the
/// caller supplies an explicit `step`. When one perturbed side leaves the axis's physical
/// domain, falls back to a one-sided difference -- see the module doc's "One-sided fallback"
/// section and [`Derivative::scheme`](Derivative#structfield.scheme).
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
/// - [`KernelError::AxisUnsupportedForRequest`] or [`KernelError::AxisAbsent`], propagated
///   unchanged from [`with_axis`] -- these depend only on `axis`/`base`, never on the perturbed
///   value, so they surface immediately with no fallback attempted (see the module doc).
/// - [`KernelError::StepOutOfDomain`] if BOTH perturbed sides fail to evaluate -- there is no
///   data left to build even a one-sided difference from.
/// - [`KernelError::Solve`] or [`KernelError::Observation`], propagated unchanged from
///   [`evaluate`], if the UNPERTURBED value at `x` itself also fails to evaluate during a
///   one-sided fallback (see the module doc) -- a degenerate case distinct from
///   `StepOutOfDomain` in that even the base value is unusable, not just the step.
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

    // with_axis's own failure modes here (AxisUnsupportedForRequest, AxisAbsent, TypeMismatch)
    // depend only on `axis` and `base`'s OTHER fields, never on the specific value written --
    // x+h and x-h would fail identically -- so they propagate immediately via `?`, from
    // whichever side is built first, exactly as before this function grew a one-sided fallback.
    let plus_req = with_axis(base, axis, AxisValue::Scalar(x + h))?;
    let minus_req = with_axis(base, axis, AxisValue::Scalar(x - h))?;

    // Unlike with_axis's structural errors above, a failure HERE is about the specific
    // perturbed VALUE (e.g. a negative wind speed, or a query range that fell outside a shrunk
    // max_range_m) -- see the module doc's "One-sided fallback".
    let plus_result = evaluate(&plus_req, ranges_m);
    let minus_result = evaluate(&minus_req, ranges_m);

    // `hi`/`lo` always mean "the sample at the higher x" / "the sample at the lower x", so the
    // derivative below is always (hi - lo) / denom regardless of which scheme produced them:
    //   Central:  hi = f(x+h), lo = f(x-h), denom = 2h
    //   Forward:  hi = f(x+h), lo = f(x),   denom = h   (minus side left the domain)
    //   Backward: hi = f(x),   lo = f(x-h), denom = h   (plus side left the domain)
    let (hi, lo, denom, scheme) = match (plus_result, minus_result) {
        (Ok(p), Ok(m)) => (p, m, 2.0 * h, DifferenceScheme::Central),
        (Ok(p), Err(_)) => {
            let base_req: crate::solve_json::SolveRequestV1 = base.into();
            let f_x = evaluate(&base_req, ranges_m)?;
            (p, f_x, h, DifferenceScheme::ForwardOneSided)
        }
        (Err(_), Ok(m)) => {
            let base_req: crate::solve_json::SolveRequestV1 = base.into();
            let f_x = evaluate(&base_req, ranges_m)?;
            (f_x, m, h, DifferenceScheme::BackwardOneSided)
        }
        (Err(_), Err(_)) => {
            return Err(KernelError::StepOutOfDomain { axis, attempted: h });
        }
    };

    debug_assert_eq!(hi.len(), lo.len());
    let mut out = Vec::with_capacity(ranges_m.len());
    for (a, b) in hi.iter().zip(lo.iter()) {
        let d_drop = (a.drop_m - b.drop_m) / denom;
        let d_wind = (a.windage_m - b.windage_m) / denom;
        if !d_drop.is_finite() || !d_wind.is_finite() {
            return Err(KernelError::NonFinite(axis));
        }
        out.push(Derivative {
            axis,
            range_m: a.range_m,
            d_drop_d_x: d_drop,
            d_windage_d_x: d_wind,
            step_used: h,
            scheme,
        });
    }
    Ok(out)
}

/// Iteration cap for [`bisect_axis`]'s search -- matches the existing inverse-solver cap
/// (`HoldCurve::range_for_angular_drop_mil`'s `INVERSE_MAX_ITERATIONS`, `src/main.rs`).
pub const BISECTION_MAX_ITERATIONS: u32 = 80;

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
    for _ in 0..BISECTION_MAX_ITERATIONS {
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
        // Review fix (post-I4): the common both-sides-valid path must still report Central, not
        // default or drift to a one-sided scheme when nothing forced one.
        assert_eq!(d[0].scheme, DifferenceScheme::Central);
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
        // Read back off the RESOLVED request rather than hardcoding 800.0 a second time here:
        // otherwise a future edit to the fixture's muzzle_velocity_mps would silently compare
        // the real derivative against a closed form for the WRONG v (review minor fix).
        let v = r.rifle.muzzle_velocity_mps;
        let expected = -9.80665 * x * x / (v * v * v);
        let d = central_difference(&r, InputAxis::MuzzleVelocityMps, &[x], None).unwrap();
        let rel = ((d[0].d_drop_d_x - expected) / expected).abs();
        assert!(rel < 0.02, "expected ~{expected}, got {} (rel {rel})", d[0].d_drop_d_x);
        assert_eq!(d[0].scheme, DifferenceScheme::Central);
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
        // Read back off the RESOLVED request rather than hardcoding 800.0 a second time here
        // (review minor fix -- same rationale as the single-range oracle test above).
        let v = r.rifle.muzzle_velocity_mps;
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
        assert!(d.iter().all(|d| d.scheme == DifferenceScheme::Central));
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
    /// error budget reads it directly to scale an input uncertainty. This test only pins the
    /// REPORTED value against an independently recomputed formula, decoupled from any
    /// particular axis's physical response; it does NOT by itself prove the same `h` was the
    /// one actually used inside the `2h` division (a bug that reports one `h` but divides by
    /// another would still pass this test unchanged) -- that is what the vacuum oracle test
    /// above protects, by checking the resulting NUMBER rather than this metadata field
    /// (review fix: this comment previously overclaimed what this test covers).
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

    /// I1 review fix: `Ok(None)` means only "no flip across this domain" -- NOT "predicate
    /// holds throughout". This is the mirror of the test above: a target so far below any
    /// reachable drop that the predicate is FALSE at both ends must ALSO come back as `None`,
    /// indistinguishable at the type level from the "true at both ends" case above. That
    /// asymmetry (identical `Ok(None)` for two opposite facts) is exactly what the module doc's
    /// "Bisection contract" section now warns callers about -- this test only pins down that
    /// `bisect_axis` itself does not silently pick one interpretation (e.g. by fabricating a
    /// crossing when the predicate happens to be false everywhere).
    #[test]
    fn bisect_axis_returns_none_when_the_predicate_is_false_at_both_ends() {
        let r = resolved(823.0);
        let found = bisect_axis(
            &r,
            InputAxis::MuzzleVelocityMps,
            600.0,
            (600.0, 1100.0),
            &|o: &Observation| o.drop_m < -1000.0, // never true: drop cannot be this negative
            0.05,
        )
        .unwrap();
        assert!(
            found.is_none(),
            "predicate is false everywhere on this domain; bisect_axis must report None, got \
             {found:?}"
        );
    }

    /// I3 review fix: deleting `bisect_axis`'s categorical guard fails nothing today -- with it
    /// gone, a categorical axis would fall through to `with_axis`, which reports `TypeMismatch`
    /// (since `AxisValue::Scalar` never matches a categorical axis's expected representation),
    /// not the `CategoricalAxis` this function's own `# Errors` section promises. Pin it down
    /// directly and independently of `central_difference`'s own (separately tested) categorical
    /// guard.
    #[test]
    fn bisect_axis_refuses_categorical_axes() {
        let r = resolved(823.0);
        let e = bisect_axis(
            &r,
            InputAxis::CoriolisEnabled,
            600.0,
            (0.0, 1.0),
            &|o: &Observation| o.drop_m < 1.0,
            0.1,
        );
        assert!(matches!(e, Err(KernelError::CategoricalAxis(InputAxis::CoriolisEnabled))));
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

    /// I2 review fix: no test anywhere in this file ever passes `Some(...)` for the explicit
    /// `step` parameter, so a mutation that silently discards the caller's step and falls back
    /// to the default formula would pass every other test in this file unchanged. MBA-1347 is
    /// the obvious consumer of an explicit step ("perturb by 1 sigma"). The chosen step (2.0) is
    /// deliberately different from what the default formula would compute here (0.8, per
    /// `central_difference_matches_the_vacuum_analytic_derivative`), so a bug that ignores
    /// `step` cannot coincidentally pass by landing on the same number anyway.
    #[test]
    fn explicit_step_overrides_the_default_and_still_matches_the_vacuum_oracle() {
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
        let v = r.rifle.muzzle_velocity_mps;
        let x = 400.0_f64;
        let custom_step = 2.0_f64;
        let d = central_difference(&r, InputAxis::MuzzleVelocityMps, &[x], Some(custom_step))
            .unwrap();
        assert_eq!(
            d[0].step_used, custom_step,
            "an explicit step must override the default formula, not just be ignored"
        );
        let expected = -9.80665 * x * x / (v * v * v);
        let rel = ((d[0].d_drop_d_x - expected) / expected).abs();
        assert!(rel < 0.02, "expected ~{expected}, got {} (rel {rel})", d[0].d_drop_d_x);
        assert_eq!(d[0].scheme, DifferenceScheme::Central);
    }

    /// I2 review fix, continued: a non-finite or non-positive explicit step must be rejected,
    /// not silently used -- dividing by zero/NaN, or by a negative number (which would silently
    /// flip the derivative's sign without changing anything else about the result's shape).
    #[test]
    fn non_finite_or_non_positive_explicit_step_is_rejected() {
        let r = resolved(823.0);
        let nan = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], Some(f64::NAN));
        assert!(matches!(nan, Err(KernelError::NonFinite(InputAxis::MuzzleVelocityMps))));
        let zero = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], Some(0.0));
        assert!(matches!(zero, Err(KernelError::NonFinite(InputAxis::MuzzleVelocityMps))));
        let negative = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], Some(-1.0));
        assert!(matches!(negative, Err(KernelError::NonFinite(InputAxis::MuzzleVelocityMps))));
    }

    /// I4 review fix: `WindSpeed` at 0.0 (still air, the default absent any `wind` block at
    /// all, and the SAME magnitude used by BOTH of this file's vacuum-oracle fixtures, so this
    /// is not an exotic edge case) cannot be centrally differentiated. The minus side needs
    /// `-min_abs_step` (-0.05 m/s), which `resolve_wind`'s `require_non_negative
    /// ("$.wind.speed_mps")` (`src/solve_v1.rs`) rejects, since wind speed is a non-negative
    /// magnitude. Falling back to a ONE-SIDED forward difference is not just an accommodation
    /// for the solver rejecting an unphysical input -- it is the numerically correct thing to
    /// do here, since windage as a function of SIGNED wind speed is V-shaped (not smooth) at
    /// zero, so a central difference straddling it would differentiate through a kink.
    ///
    /// Speed and direction must be supplied TOGETHER (`resolve_wind` rejects one without the
    /// other), so still air is expressed here as an EXPLICIT `speed_mps: 0.0` paired with a
    /// 90-degree crosswind direction, rather than omitting the `wind` block entirely: the
    /// omitted-block default direction is a pure head/tailwind (no crosswind component at all),
    /// which would make `d_windage_d_x` genuinely, correctly zero regardless of whether the
    /// fallback logic is right -- that would test nothing.
    #[test]
    fn wind_speed_falls_back_to_one_sided_in_still_air() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {},
            "wind": {"speed_mps": 0.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(
            match &r.wind {
                crate::solve_json::ResolvedWindV1::Constant(c) => c.speed_mps,
                crate::solve_json::ResolvedWindV1::Segmented(_) => panic!("constant wind expected"),
            },
            0.0,
            "fixture assumption: still air"
        );

        let d = central_difference(&r, InputAxis::WindSpeed, &[600.0], None);
        let d = d.unwrap_or_else(|e| panic!("still air must fall back, not error, got {e:?}"));
        assert_eq!(d.len(), 1);
        assert_eq!(
            d[0].scheme,
            DifferenceScheme::ForwardOneSided,
            "the minus side (-0.05 m/s) leaves WindSpeed's domain, so this must be a forward \
             one-sided fallback, not Central and not Backward"
        );
        assert!(d[0].d_windage_d_x.is_finite() && d[0].d_windage_d_x != 0.0);
    }

    /// I4 review fix, continued: `RelativeHumidity` hits the same class of domain violation at
    /// its OTHER boundary. `resolve_atmosphere` validates it to `require_range(..., 0.0, 1.0)`
    /// (`src/solve_v1.rs`); at the dry-air boundary (0.0, itself a common, unremarkable request
    /// -- not every caller supplies humidity, and dry air is a normal thing to model
    /// explicitly), the minus side (`-min_abs_step` = -0.001) falls outside that range.
    #[test]
    fn relative_humidity_falls_back_to_one_sided_at_the_dry_air_boundary() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"relative_humidity": 0.0},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(
            r.atmosphere.relative_humidity, 0.0,
            "fixture assumption: an explicit 0.0 must resolve to exactly 0.0, not the 0.5 \
             literal default"
        );

        let d = central_difference(&r, InputAxis::RelativeHumidity, &[600.0], None);
        let d = d.unwrap_or_else(|e| panic!("the dry-air boundary must fall back, not error, got {e:?}"));
        assert_eq!(d.len(), 1);
        assert_eq!(
            d[0].scheme,
            DifferenceScheme::ForwardOneSided,
            "the minus side (-0.001) leaves RelativeHumidity's [0, 1] domain"
        );
        assert!(d[0].d_drop_d_x.is_finite());
    }

    /// I4 review fix, continued: `TargetDistance` (== `shot.max_range_m`) hits the same class
    /// of domain violation through a DIFFERENT `KernelError` variant than WindSpeed/
    /// RelativeHumidity above -- `Observation` (a query range now outside the computed
    /// trajectory), raised AFTER a successful solve, rather than `Solve` (a validation failure
    /// before one). Differentiating drop's sensitivity to the target distance AT that same
    /// target distance is the natural use of this axis (e.g. "how much does drop at 900 m
    /// change if the target were a little closer or farther"); querying a range close to `x`
    /// means the minus side (a slightly SMALLER max_range_m) can no longer see it.
    ///
    /// x = 900, h = max(900*1e-3, 0.5) = 0.9, so minus = 899.1 and plus = 900.9. A query at
    /// 899.5 m is comfortably inside plus's and the unperturbed x's trajectories but strictly
    /// beyond minus's -- deliberately not exactly 900.0, to avoid depending on whether an exact
    /// upper-endpoint query is treated as inclusive.
    #[test]
    fn target_distance_falls_back_to_one_sided_when_queried_near_its_own_max_range() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {}, "wind": {}, "solver": {}, "effects": {},
            "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(r.shot.max_range_m, 900.0);

        let d = central_difference(&r, InputAxis::TargetDistance, &[899.5], None);
        let d = d.unwrap_or_else(|e| {
            panic!("a range just inside the unperturbed max_range_m must fall back, not error, got {e:?}")
        });
        assert_eq!(d.len(), 1);
        assert_eq!(
            d[0].scheme,
            DifferenceScheme::ForwardOneSided,
            "the minus side (max_range_m = 899.1) no longer covers the 899.5 m query, so this \
             must be a forward one-sided fallback"
        );
        assert!(d[0].d_drop_d_x.is_finite());
    }

    /// I4 review fix, continued: when BOTH perturbed sides fail, there is no data left to build
    /// even a one-sided difference from. An oversized EXPLICIT step pushes `RelativeHumidity`
    /// (at the 0.0 boundary) out of its `[0, 1]` domain in BOTH directions at once (x+h = 2.0,
    /// x-h = -2.0), so this cannot fall back either way. `central_difference` must report the
    /// distinct `StepOutOfDomain` -- not propagate whichever of the two solve errors happened to
    /// be checked first (arbitrary, and neither alone names the real problem), and not silently
    /// claim a zero-effect derivative.
    #[test]
    fn both_sides_out_of_domain_reports_step_out_of_domain() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"relative_humidity": 0.0},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(r.atmosphere.relative_humidity, 0.0, "fixture assumption");

        let e = central_difference(&r, InputAxis::RelativeHumidity, &[600.0], Some(2.0));
        match e {
            Err(KernelError::StepOutOfDomain { axis: InputAxis::RelativeHumidity, attempted }) => {
                assert_eq!(attempted, 2.0);
            }
            other => panic!("expected StepOutOfDomain, got {other:?}"),
        }
    }
}
