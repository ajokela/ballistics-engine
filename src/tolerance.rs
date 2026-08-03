//! MBA-1350: how wrong may one input be before the shot leaves the target?
//!
//! This is the deterministic inverse of a WEZ sweep. For a nominal (already dialed-in, centred)
//! solution and an explicit target, each requested axis is bisected outward -- toward its
//! configured domain's lower bound (`near_bound`) and, independently, toward its upper bound
//! (`far_bound`) -- until the impact crosses the target boundary. Bounds are strictly ONE
//! VARIABLE AT A TIME: two inputs each at their own individual limit will generally miss even
//! though neither alone would, and no probability is attached to any bound reported here. Both
//! of those are stated in [`ToleranceReportV1::assumptions`] itself, not only in this comment --
//! see `the_report_refuses_to_imply_probability` in this module's tests.
//!
//! Built entirely on the existing kernel: [`bisect_axis`] for the search,
//! [`evaluate`]/[`read_axis`]/[`with_axis`] for reading and rebuilding requests, and
//! [`TargetGeometryV1`] (Task 11, `crate::error_budget`) for the target shape -- reused verbatim,
//! not re-defined, including its "always centred on the nominal impact point" semantics.
//!
//! # The central hazard: `Ok(None)` means two different things
//!
//! [`bisect_axis`]'s own contract (`crate::perturbation::derive`'s module doc, "Bisection
//! contract") is explicit that `Ok(None)` means ONLY "the predicate did not change truth value
//! across this domain" -- nothing more. In general, for a two-sided "stays inside a region"
//! predicate, that single fact is consistent with TWO opposite situations that look identical at
//! the type level: the region is never exited, or the search never started inside it at all.
//!
//! **In this module specifically, the second situation has exactly one door.** `target`
//! ([`TargetGeometryV1`]) has no offset field at all -- every variant is, by that type's own
//! doc, ALWAYS defined centred on the nominal impact point. That makes the nominal itself, read
//! against its own centred target, inside by construction for any target with positive area
//! (`dy == dz == 0.0` trivially satisfies `dy.abs() <= height_m / 2.0` etc. whenever `height_m`/
//! `width_m`/`radius_m` is positive). The ONLY way the nominal can fail that check is a
//! DEGENERATE target -- non-positive width, height, or radius, which contains no point at all,
//! not even its own centre. A target genuinely offset from the nominal, or a nominal that is
//! "not really centred" for some other reason, is not an expressible input to this API at all --
//! `TargetGeometryV1` has no field that could carry it -- so that reading of "outside throughout"
//! cannot be the live case here, and this doc previously implied otherwise.
//!
//! **This module still checks it explicitly rather than assuming it**, via a pre-flight step
//! that evaluates the shared anchor point itself -- the IDENTICAL `with_axis`/`evaluate`
//! reconstruction [`bisect_axis`] uses for `domain.0` in both the near and far searches (see
//! [`tolerance_envelope`]'s "Pre-flight" step below) -- BEFORE trusting any `Ok(None)` it
//! returns. This closes the only door that DOES exist (a degenerate target, or a reconstruction
//! bug in `with_axis`/`evaluate` themselves) rather than guarding against a routinely-reachable
//! "off-centre nominal" scenario the type system cannot currently express; it is still worth
//! keeping for that narrower reason, and is still what makes the following distinction sound
//! rather than assumed. If the anchor check reads as "outside" --
//! [`ToleranceAxisV1::nominal_inside_target`] is `false` -- the search is never even run: see
//! correctness requirement 1 below. Only once the anchor is CONFIRMED "inside" does an `Ok(None)`
//! from a search direction unambiguously mean "stays inside throughout that direction," recorded
//! as [`ToleranceAxisV1::unbounded_in_domain`]. The two outcomes share the same `_bound: None`
//! shape but are otherwise completely distinct fields -- see
//! `a_degenerate_target_is_flagged_nominal_outside_not_confused_with_unbounded` and
//! `an_axis_that_never_exits_is_flagged_not_bounded` in this module's tests, which produce
//! IDENTICAL `near_bound`/`far_bound` (`None`/`None`) from two DIFFERENT root causes and assert
//! that `nominal_inside_target`/`unbounded_in_domain` tell them apart.
//!
//! One limitation this module inherits rather than solves: [`bisect_axis`] itself assumes the
//! predicate changes truth value AT MOST ONCE across a search direction. If an axis's effect on
//! impact were non-monotonic over the ENTIRE configured domain (physically unusual, but not
//! impossible for a very wide domain), a single bisection could in principle miss an excursion
//! that both re-enters "inside" before reaching the domain edge -- this is a pre-existing,
//! documented property of [`bisect_axis`] itself (see its module doc), not something specific to
//! this module's use of it, and not something a bounded amount of extra work here can fully
//! close without a much more expensive full-domain scan. Choosing a domain scoped to where the
//! axis is plausibly monotonic is the caller's responsibility, exactly as it already is for every
//! other consumer of [`bisect_axis`].
//!
//! # The four correctness requirements
//!
//! 1. **The nominal must actually read as inside before a bound means anything.** See "The
//!    central hazard" above -- enforced by the pre-flight check, not assumed.
//! 2. **Bounds are one axis at a time and may not be assumed simultaneously.** Stated in
//!    [`ToleranceReportV1::assumptions`].
//! 3. **No probability is implied.** Also stated in `assumptions`.
//! 4. **A bound is never extrapolated beyond the caller's configured domain.** [`bisect_axis`]
//!    itself never looks outside `domain`; this module additionally never INVENTS a domain the
//!    caller did not supply -- see "Domains are validated up front, per axis" below and
//!    [`KernelError::InvalidDomain`].
//!
//! # Domains are validated up front, per axis
//!
//! Unlike an earlier draft of this feature, there is no implicit fallback domain (e.g.
//! `nominal * 0.5 ..= nominal * 1.5`): that specific formula is not just unspecified but actively
//! wrong whenever an axis's nominal value is exactly `0.0` (routine for `WindDirection`, `Cant`,
//! `ShootingAngle`, and several others), where it collapses to a zero-width `(0.0, 0.0)` domain
//! that can never be searched. Every axis in `axes` must have a corresponding entry in `domains`
//! with both bounds finite, the lower bound strictly less than the upper, and the axis's own
//! current value strictly between them (not merely `<=`/`>=`: a nominal value sitting exactly AT
//! one edge would make that direction's search a zero-width probe, i.e. `bisect_axis`'s two
//! endpoints would be the same point and it would trivially -- and misleadingly -- report
//! `Ok(None)`). A missing or invalid domain returns [`KernelError::InvalidDomain`] immediately,
//! before any solve.
//!
//! # Unavailable axes
//!
//! [`bisect_axis`]/[`with_axis`] can legitimately refuse to search a declared axis:
//! [`KernelError::CategoricalAxis`] (an effect toggle -- no numeric domain to bisect),
//! [`KernelError::AxisAbsent`] (a wind axis under segmented wind -- no single scalar value to
//! hold at a nominal or perturb), or [`KernelError::AxisUnsupportedForRequest`] (`Altitude` under
//! a QNH-referenced atmosphere, `ShotAzimuth` under compass-referenced wind). Every one of these
//! is recorded in [`ToleranceReportV1::unavailable_axes`] (axis, machine-readable
//! [`UnavailableReasonCodeV1`], human-readable reason) and the rest of the report is still
//! produced from whatever axes DID evaluate -- silently dropping an unavailable axis would look
//! identical to "this axis was searched and found to have no bound," which is a different fact.
//!
//! This reuses `crate::error_budget`'s existing four-way classification verbatim (widened to
//! `pub(crate)` for this module) rather than defining a second, independently maintained copy of
//! the same split -- see that function's own doc comment for the one caveat this creates (two of
//! its four reason strings read as written for differentiation/uncertainty, and
//! `StepOutOfDomain` is not actually reachable through `bisect_axis` today, only through
//! `central_difference`; both are accepted trade-offs of sharing one classifier).
//!
//! Any OTHER error (`Solve`, `Observation`, the defensive `TypeMismatch`/`NonFinite`,
//! `DuplicateAxis`, or this module's own [`KernelError::InvalidDomain`]) is a genuine failure --
//! most notably a domain whose lower bound, for the `TargetDistance` axis, dips below the
//! caller's own `range_m` (an internally inconsistent request: "how close could the target
//! plausibly be" bisected past "the range I am currently asking about") -- and propagates
//! immediately, aborting the whole report rather than mislabelling it as a per-axis refusal. See
//! `a_genuine_observation_error_propagates_not_recorded_as_unavailable`.
//!
//! # `TargetDistance` cannot answer "how wrong may my range estimate be"
//!
//! The ticket's own motivating example is a range-estimate question -- "my rangefinder read
//! 600 m; how far off could that reading be before I miss?" The axis that superficially matches,
//! `InputAxis::TargetDistance` (`shot.max_range_m`), does NOT answer it. `TargetDistance` is
//! `requires_rezero: false` (`crate::perturbation::axis_meta`): perturbing it changes only how
//! far the trajectory is COMPUTED, never the muzzle angle or any sight correction, so it has NO
//! effect at all on the impact observed at a fixed `range_m` as long as the perturbed
//! `max_range_m` stays at or above `range_m`. Bisecting it reports either a hard
//! [`KernelError::Observation`] (once the search dips below `range_m`, see above) or
//! `unbounded_in_domain: true` with [`ToleranceAxisV1::near_has_no_effect`]/
//! [`ToleranceAxisV1::far_has_no_effect`] both `true` -- literally correct about this specific
//! axis, but not an answer to the range-estimate question, and a caller who reads
//! `unbounded_in_domain: true` alone (without also checking the `_has_no_effect` flags) could
//! easily mistake "this axis has no effect" for "your range estimate can be off by any amount
//! and still hit," which is not what it means here.
//!
//! The axis that actually answers a range-estimate question is `InputAxis::ZeroDistance`
//! (`shot.zero_distance_m`, `requires_rezero: true`): perturbing it re-runs the elevation search
//! for a DIFFERENT assumed zero distance, producing a different effective muzzle angle, and only
//! THEN observes the impact at the caller's own fixed, unchanged `range_m` -- i.e. "if I had
//! dialed for a different distance than the true one, how far off would I land at the true
//! range," which is what a rangefinder error actually does to a shot. See
//! `target_distance_axis_shows_no_measurable_effect_not_a_generic_unbounded_claim` in this
//! module's tests, which pins both the `TargetDistance` "no effect" finding and the contrast
//! against an axis that genuinely does move the impact but merely stays inside a large target.
//!
//! # Cost
//!
//! Per axis: one pre-flight solve, then up to [`crate::perturbation::derive::BISECTION_MAX_ITERATIONS`]
//! solves for EACH of the two search directions (far fewer in practice -- an `Ok(None)` result
//! costs only the two domain-endpoint evaluations, and a found crossing converges to this
//! module's own domain-relative tolerance in on the order of twenty iterations for a typical
//! domain width, not the full cap), plus one more solve for each direction that DID find a bound
//! (to classify which edge of the target it crosses). A `requires_rezero` axis
//! (`crate::perturbation::axis_meta`) multiplies each of those solves by the elevation search's
//! own cost, exactly as it does for [`crate::perturbation::derive::central_difference`] and
//! [`crate::error_budget::error_budget`] -- unavoidable here, not something this module changes.

use serde::{Deserialize, Serialize};

use crate::error_budget::{unavailable_reason, TargetGeometryV1, UnavailableReasonCodeV1};
use crate::perturbation::access::{read_axis, with_axis, AxisValue, KernelError};
use crate::perturbation::derive::bisect_axis;
use crate::perturbation::taxonomy::{axis_meta, AxisKind, InputAxis};
use crate::perturbation::{evaluate, Observation};
use crate::solve_json::ResolvedSolveRequestV1;
use crate::trajectory_observation::TrajectoryObservationError;

/// Schema version for [`ToleranceReportV1`].
pub const TOLERANCE_SCHEMA_VERSION_V1: u32 = 1;

/// Which edge of a [`TargetGeometryV1::Rect`] a found bound crosses, from the shooter's own
/// point of view: `Top`/`Bottom` along the drop axis, `Left`/`Right` along the windage axis.
/// `drop_m` is positive BELOW the line of sight and `windage_m` is positive to the shooter's
/// RIGHT (see [`crate::perturbation::Observation`]), so more drop than nominal crosses the
/// `Bottom` edge and more rightward windage than nominal crosses the `Right` edge.
/// [`TargetGeometryV1::Circle`] has no distinct edges and always reports `Radial`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LimitingBoundaryV1 {
    Top,
    Bottom,
    Left,
    Right,
    Radial,
}

/// One axis's tolerance envelope: how far it may move from its own current value, in EACH
/// direction independently, before the impact crosses `target`'s boundary.
///
/// # Reading these fields together
///
/// - **A real bound was found** in a direction: that direction's `_bound` field is `Some`
///   (verified, in this module's tests, by re-solving AT the bound via a path independent of
///   `with_axis`/`bisect_axis` and confirming it sits on the boundary).
/// - **Confirmed to stay inside throughout the configured domain**: `nominal_inside_target` is
///   `true` and the corresponding `_bound` field(s) are `None`. `unbounded_in_domain` is `true`
///   exactly when this holds in BOTH directions.
/// - **The nominal itself does not read as inside `target`**: `nominal_inside_target` is `false`.
///   Neither direction is even searched -- both `_bound` fields are `None` and
///   `unbounded_in_domain` is `false`, NOT `true`. This is the module's central distinction; see
///   the module doc's "The central hazard" section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ToleranceAxisV1 {
    pub axis: InputAxis,
    /// This axis's own current (unperturbed) value, in its own physical unit
    /// (`crate::perturbation::axis_meta(axis).kind`).
    pub nominal: f64,
    /// Whether the impact, with this axis pinned AT `nominal` via the identical
    /// `with_axis`/`evaluate` reconstruction the bisection itself uses for the shared anchor
    /// point, reads as inside `target`. Correctness requirement 1: bounds are only searched for
    /// when this is `true`. Virtually always `true` in ordinary use -- perturbing an axis back
    /// to its own current value must reproduce the request's own nominal impact -- `false` is
    /// reserved for a target with no positive area (a [`TargetGeometryV1`] that cannot contain
    /// any point, not even its own centre); see this module's tests for how that is exercised
    /// without relying on self-consistency.
    pub nominal_inside_target: bool,
    /// The bound found searching FROM `nominal` TOWARD this axis's configured domain's lower
    /// bound, or `None` -- see the struct doc for what `None` means here.
    pub near_bound: Option<f64>,
    /// The bound found searching FROM `nominal` TOWARD this axis's configured domain's upper
    /// bound, or `None` -- see the struct doc for what `None` means here.
    pub far_bound: Option<f64>,
    /// Which edge of `target` `near_bound` crosses, or `None` exactly when `near_bound` is
    /// `None`.
    pub near_limiting_boundary: Option<LimitingBoundaryV1>,
    /// Which edge of `target` `far_bound` crosses, or `None` exactly when `far_bound` is `None`.
    pub far_limiting_boundary: Option<LimitingBoundaryV1>,
    /// `true` exactly when `nominal_inside_target` is `true` AND both `near_bound` and
    /// `far_bound` are `None`: the impact is confirmed to stay inside `target` across the WHOLE
    /// configured domain (subject to `bisect_axis`'s own one-crossing assumption -- see the
    /// module doc). Deliberately NOT simply `near_bound.is_none() && far_bound.is_none()`: that
    /// would also be `true` when `nominal_inside_target` is `false`, which is the exact
    /// conflation this ticket exists to prevent.
    pub unbounded_in_domain: bool,
    /// `true` exactly when `near_bound` is `None` AND the impact at the domain's own lower edge
    /// is indistinguishable (within `1e-9` m on both drop and windage) from the nominal impact:
    /// this axis has NO MEASURABLE EFFECT on the observed impact in this direction at all, as
    /// distinct from an axis that DOES move the impact but never far enough to leave `target`.
    /// See "`TargetDistance` cannot answer..." in the module doc for the motivating example
    /// (`TargetDistance`'s `requires_rezero: false` means it never affects the impact at a fixed
    /// `range_m` at all, so `unbounded_in_domain: true` for it means "this axis is provably
    /// irrelevant here," not "any error in it is safe"). Always `false` when `near_bound` is
    /// `Some` (a bound exists, so the axis clearly has SOME effect) or `nominal_inside_target` is
    /// `false` (the direction was never searched).
    pub near_has_no_effect: bool,
    /// As `near_has_no_effect`, for the far direction.
    pub far_has_no_effect: bool,
    /// The smallest linear half-extent of `target` (`min(height_m, width_m) / 2.0` for a
    /// [`TargetGeometryV1::Rect`], `radius_m` for a [`TargetGeometryV1::Circle`]) -- a
    /// convenience summary of how tight the target is, identical across every axis in a report
    /// since it depends only on `target`.
    pub margin_linear_m: f64,
}

/// One requested axis [`tolerance_envelope`] could not search at all, and why -- distinct from an
/// axis that WAS searched and found to have no bound (see [`ToleranceAxisV1`]). Never silently
/// dropped: see this module's "Unavailable axes" doc section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnavailableAxisV1 {
    pub axis: InputAxis,
    pub code: UnavailableReasonCodeV1,
    pub reason: String,
}

/// One-variable tolerance envelope report (MBA-1350): how far each requested input may drift
/// from its own current value, one at a time, before the impact leaves an explicit target.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ToleranceReportV1 {
    pub schema_version: u32,
    pub method: String,
    /// Always states, in the payload itself (not only in prose documentation): bounds are one
    /// axis at a time and may not be assumed to hold simultaneously; no probability is implied;
    /// a bound is never extrapolated past the caller's configured domain; and the nominal is
    /// confirmed inside the target before any bound is searched for. See
    /// `assumptions_cover_all_four_correctness_requirements` in this module's tests.
    pub assumptions: Vec<String>,
    pub range_m: f64,
    /// Axes [`tolerance_envelope`] could not search. See "Unavailable axes" in the module doc.
    pub unavailable_axes: Vec<UnavailableAxisV1>,
    /// One entry per requested axis that WAS searched, in the order requested.
    pub axes: Vec<ToleranceAxisV1>,
}

/// Bisection tolerance for one axis's search: proportional to the CALLER's own configured domain
/// width rather than a single fixed absolute number, so it neither wastes iterations converging
/// an oversized domain (e.g. a pressure domain spanning tens of kPa) to an absurdly tight
/// absolute tolerance, nor undershoots a narrow one (e.g. a `+/-0.05` rad cant domain). Mirrors
/// the crate's existing step-size convention (`axis_meta`'s `h = (|x| * rel).max(min_abs)`,
/// `crate::perturbation::taxonomy`) applied to a domain WIDTH instead of a value. Comfortably
/// achievable within `BISECTION_MAX_ITERATIONS` for any domain width a caller would plausibly
/// configure: halving a domain 80 times shrinks it by a factor of `2^80`, far more than the
/// `1e6` this needs.
fn bisection_tolerance(lo: f64, hi: f64) -> f64 {
    ((hi - lo).abs() * 1e-6).max(1e-9)
}

/// The smallest linear half-extent of `target` -- see [`ToleranceAxisV1::margin_linear_m`].
fn target_margin_linear_m(target: TargetGeometryV1) -> f64 {
    match target {
        TargetGeometryV1::Rect { width_m, height_m } => {
            (height_m.max(0.0) / 2.0).min(width_m.max(0.0) / 2.0)
        }
        TargetGeometryV1::Circle { radius_m } => radius_m.max(0.0),
    }
}

/// Whether `o`'s (drop, windage) deviation from `nominal` falls inside `target`, which is always
/// centred on `nominal` -- see [`TargetGeometryV1`]'s own doc. A target with no positive area
/// (non-positive width/height/radius) can never contain any point, not even `nominal` itself
/// (`dy == dz == 0`): this is the mechanism behind [`ToleranceAxisV1::nominal_inside_target`]
/// being `false` for a degenerate target, mirroring `crate::error_budget::p_hit_bivariate`'s
/// identical guard for the identical type.
fn inside(o: &Observation, nominal: &Observation, target: TargetGeometryV1) -> bool {
    let dy = o.drop_m - nominal.drop_m;
    let dz = o.windage_m - nominal.windage_m;
    match target {
        TargetGeometryV1::Rect { width_m, height_m } => {
            width_m > 0.0
                && height_m > 0.0
                && dy.abs() <= height_m / 2.0
                && dz.abs() <= width_m / 2.0
        }
        TargetGeometryV1::Circle { radius_m } => {
            radius_m > 0.0 && (dy * dy + dz * dz).sqrt() <= radius_m
        }
    }
}

/// Whether `o`'s (drop, windage) is indistinguishable from `nominal`'s to within `1e-9` m on
/// each -- see [`ToleranceAxisV1::near_has_no_effect`]/[`ToleranceAxisV1::far_has_no_effect`].
/// `1e-9` m is far tighter than any genuine physical sensitivity this crate's taxonomy produces
/// (see `crate::perturbation::taxonomy`'s `min_abs_step` values, none smaller than `1e-7`) and
/// far looser than floating-point noise from an identical, deterministic re-solve (no
/// axis/request combination this module reaches involves randomness), so this only fires for an
/// axis that is truly, not merely negligibly, disconnected from the observed impact.
fn observation_matches_nominal(o: &Observation, nominal: &Observation) -> bool {
    (o.drop_m - nominal.drop_m).abs() < 1e-9 && (o.windage_m - nominal.windage_m).abs() < 1e-9
}

/// Which edge of `target` `o` (already known to be outside, or exactly on the boundary of,
/// `target` relative to `nominal`) crosses. For a [`TargetGeometryV1::Rect`], compares
/// `dy.abs() / (height_m / 2)` against `dz.abs() / (width_m / 2)` (whichever ratio is larger is
/// the edge actually being exceeded), rewritten as `dy.abs() * width_m >= dz.abs() * height_m`
/// to avoid dividing by a half-extent that could be zero.
fn limiting_boundary(
    o: &Observation,
    nominal: &Observation,
    target: TargetGeometryV1,
) -> LimitingBoundaryV1 {
    let dy = o.drop_m - nominal.drop_m;
    let dz = o.windage_m - nominal.windage_m;
    match target {
        TargetGeometryV1::Circle { .. } => LimitingBoundaryV1::Radial,
        TargetGeometryV1::Rect { width_m, height_m } => {
            if dy.abs() * width_m.max(0.0) >= dz.abs() * height_m.max(0.0) {
                if dy > 0.0 {
                    LimitingBoundaryV1::Bottom
                } else {
                    LimitingBoundaryV1::Top
                }
            } else if dz > 0.0 {
                LimitingBoundaryV1::Right
            } else {
                LimitingBoundaryV1::Left
            }
        }
    }
}

/// How far may each of `axes` drift from its own current value, one at a time, before the impact
/// leaves `target` -- see the module doc for the full contract, and especially "The central
/// hazard" for what `unbounded_in_domain`/`nominal_inside_target` mean together.
///
/// `domains` supplies, for each axis in `axes`, the `(lower, upper)` bound to search within; see
/// "Domains are validated up front, per axis" in the module doc.
///
/// # Errors
///
/// - [`KernelError::Observation`] immediately, before any solve, if `range_m` is not finite or
///   falls outside `[0, base.shot.max_range_m]` -- the same check
///   `crate::error_budget::error_budget` applies to its own `ranges_m`.
/// - [`KernelError::InvalidDomain`] immediately, before any solve for that axis, if `domains` has
///   no entry for one of `axes`, either bound is not finite, the lower bound is not strictly less
///   than the upper, or the axis's own current value does not sit strictly between them.
/// - Otherwise, propagates any [`KernelError`] that is not one of the four structural refusals
///   recorded in [`ToleranceReportV1::unavailable_axes`] instead -- see "Unavailable axes" in the
///   module doc.
pub fn tolerance_envelope(
    base: &ResolvedSolveRequestV1,
    axes: &[InputAxis],
    range_m: f64,
    target: TargetGeometryV1,
    domains: &[(InputAxis, (f64, f64))],
) -> Result<ToleranceReportV1, KernelError> {
    // Validate range_m up front, before any solve -- identical rationale and construction to
    // error_budget's own check on its (plural) ranges_m.
    if !range_m.is_finite() {
        return Err(KernelError::Observation(TrajectoryObservationError::NonFiniteQuery {
            distance_m: range_m,
        }));
    }
    if range_m < 0.0 || range_m > base.shot.max_range_m {
        return Err(KernelError::Observation(TrajectoryObservationError::OutOfRange {
            requested_m: range_m,
            minimum_m: 0.0,
            maximum_m: base.shot.max_range_m,
        }));
    }

    // The fixed reference point every axis's "inside" check is measured against: base's own
    // resolved inputs, solved directly (not through with_axis), exactly once.
    let nominal = evaluate(&base.into(), &[range_m])?[0];
    let margin_linear_m = target_margin_linear_m(target);

    let mut out = Vec::with_capacity(axes.len());
    let mut unavailable = Vec::new();

    for &axis in axes {
        // 1. Categorical axes have no numeric domain to bisect (mirrors bisect_axis's own
        //    guard); record, do not abort the rest of the report.
        if matches!(axis_meta(axis).kind, AxisKind::Categorical) {
            let (code, reason) = unavailable_reason(&KernelError::CategoricalAxis(axis))
                .expect("CategoricalAxis is always classified as unavailable");
            unavailable.push(UnavailableAxisV1 { axis, code, reason });
            continue;
        }

        // 2. Read the current value. `None` means the axis is structurally absent on this
        //    request (the three wind axes under segmented wind) -- record, do not abort.
        let nominal_value = match read_axis(base, axis) {
            Some(AxisValue::Scalar(x)) => x,
            // Unreachable with the current taxonomy (every Continuous axis reads back as Scalar
            // or None; Flag/DragModel/TwistDirection only ever come from Categorical axes,
            // already handled above) -- kept as a defensive catch-all, mirroring
            // central_difference's identical fallback (crate::perturbation::derive).
            Some(_) => return Err(KernelError::TypeMismatch(axis)),
            None => {
                let (code, reason) = unavailable_reason(&KernelError::AxisAbsent(axis))
                    .expect("AxisAbsent is always classified as unavailable");
                unavailable.push(UnavailableAxisV1 { axis, code, reason });
                continue;
            }
        };

        // 3. The domain must be explicitly configured and well-formed -- see "Domains are
        //    validated up front, per axis" in the module doc. Never an invented default.
        let (lo, hi) = match domains.iter().find(|(a, _)| *a == axis).map(|(_, d)| *d) {
            Some(d) => d,
            None => {
                return Err(KernelError::InvalidDomain {
                    axis,
                    reason: "no domain was configured for this axis",
                })
            }
        };
        if !(lo.is_finite() && hi.is_finite()) {
            return Err(KernelError::InvalidDomain {
                axis,
                reason: "domain bounds must both be finite",
            });
        }
        if lo >= hi {
            return Err(KernelError::InvalidDomain {
                axis,
                reason: "the domain's lower bound must be strictly less than its upper bound",
            });
        }
        if !(nominal_value > lo && nominal_value < hi) {
            return Err(KernelError::InvalidDomain {
                axis,
                reason: "the axis's own current value must sit strictly inside the configured \
                         domain, or one search direction would be a zero-width probe",
            });
        }

        // 4. Pre-flight: evaluate AT the current value via the IDENTICAL with_axis/evaluate path
        //    bisect_axis uses for the shared anchor point (domain.0 in both directions below).
        //    Surfaces AxisUnsupportedForRequest/AxisAbsent for this axis/request combination
        //    before any Ok(None) is trusted -- both depend only on axis/base, never on the
        //    specific value written (see with_axis's own doc), so a successful pre-flight here
        //    guarantees neither can recur below for ANY value bisect_axis might try -- AND
        //    supplies the reconstructed "at nominal" observation correctness requirement 1
        //    checks against.
        let pre = with_axis(base, axis, AxisValue::Scalar(nominal_value))
            .and_then(|req| evaluate(&req, &[range_m]));
        let at_nominal = match pre {
            Ok(obs) => obs[0],
            Err(e) => match unavailable_reason(&e) {
                Some((code, reason)) => {
                    unavailable.push(UnavailableAxisV1 { axis, code, reason });
                    continue;
                }
                None => return Err(e),
            },
        };

        // Correctness requirement 1: the nominal must actually read as inside before any bound
        // means anything -- checked via the SAME reconstruction path as domain.0, not merely
        // assumed. See the module doc's "The central hazard" section.
        let nominal_inside_target = inside(&at_nominal, &nominal, target);
        if !nominal_inside_target {
            out.push(ToleranceAxisV1 {
                axis,
                nominal: nominal_value,
                nominal_inside_target: false,
                near_bound: None,
                far_bound: None,
                near_limiting_boundary: None,
                far_limiting_boundary: None,
                unbounded_in_domain: false,
                near_has_no_effect: false,
                far_has_no_effect: false,
                margin_linear_m,
            });
            continue;
        }

        // With the anchor confirmed inside, any Ok(None) bisect_axis returns below is now
        // unambiguous: "stays inside throughout that direction," not "outside throughout" (the
        // conflation this module exists to prevent -- see the module doc). Any error from these
        // two calls, after a successful pre-flight, can only be a genuine Solve/Observation
        // failure (AxisUnsupportedForRequest/AxisAbsent/CategoricalAxis are all already ruled
        // out above), so it propagates directly via `?` -- see
        // `a_genuine_observation_error_propagates_not_recorded_as_unavailable`.
        let tol = bisection_tolerance(lo, hi);
        let pred = |o: &Observation| inside(o, &nominal, target);
        let near = bisect_axis(base, axis, range_m, (nominal_value, lo), &pred, tol)?;
        let far = bisect_axis(base, axis, range_m, (nominal_value, hi), &pred, tol)?;

        // Beyond finding a bound (or confirming its absence), each None direction is checked
        // for whether the axis has ANY measurable effect on the impact at all -- see
        // `near_has_no_effect`/`far_has_no_effect`'s own doc and the module doc's "TargetDistance
        // cannot answer..." section for why this distinction matters. Both endpoint observations
        // (`lo`/`hi`) are recomputed here rather than threaded out of `bisect_axis` (which does
        // evaluate them internally to decide `Ok(None)`, but does not expose them) -- at most two
        // extra solves per axis, only in the already-cheaper `None` case.
        let (near_limiting_boundary, near_has_no_effect) = match near {
            Some(v) => {
                let o = evaluate(&with_axis(base, axis, AxisValue::Scalar(v))?, &[range_m])?[0];
                (Some(limiting_boundary(&o, &nominal, target)), false)
            }
            None => {
                let o = evaluate(&with_axis(base, axis, AxisValue::Scalar(lo))?, &[range_m])?[0];
                (None, observation_matches_nominal(&o, &nominal))
            }
        };
        let (far_limiting_boundary, far_has_no_effect) = match far {
            Some(v) => {
                let o = evaluate(&with_axis(base, axis, AxisValue::Scalar(v))?, &[range_m])?[0];
                (Some(limiting_boundary(&o, &nominal, target)), false)
            }
            None => {
                let o = evaluate(&with_axis(base, axis, AxisValue::Scalar(hi))?, &[range_m])?[0];
                (None, observation_matches_nominal(&o, &nominal))
            }
        };

        out.push(ToleranceAxisV1 {
            axis,
            nominal: nominal_value,
            nominal_inside_target: true,
            near_bound: near,
            far_bound: far,
            near_limiting_boundary,
            far_limiting_boundary,
            unbounded_in_domain: near.is_none() && far.is_none(),
            near_has_no_effect,
            far_has_no_effect,
            margin_linear_m,
        });
    }

    Ok(ToleranceReportV1 {
        schema_version: TOLERANCE_SCHEMA_VERSION_V1,
        method: "one_variable_deterministic_bisection".to_string(),
        assumptions: vec![
            "Each bound holds ONE input at its limit while every other input stays at its \
             nominal value. Bounds from different axes may NOT be assumed to hold \
             simultaneously: two inputs each at their own individual limit will generally miss \
             even though neither alone would."
                .to_string(),
            "No probability distribution is assumed or implied by any bound here. These are \
             deterministic limits of a one-variable search, not confidence intervals or a \
             probability of hit."
                .to_string(),
            "A bound is reported only when found strictly within the axis's own configured \
             search domain. It is never extrapolated beyond that domain: 'no bound within the \
             domain' (unbounded_in_domain) is reported as exactly that fact, never as a guessed \
             number."
                .to_string(),
            "Before any bound is searched for, the nominal solution's own impact is confirmed to \
             read as inside the target (nominal_inside_target); an axis for which that check \
             fails is reported as such instead of a fabricated or misleading bound."
                .to_string(),
        ],
        range_m,
        unavailable_axes: unavailable,
        axes: out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_budget::TargetGeometryV1;
    use crate::perturbation::InputAxis;

    // The brief's own fixture literal uses a lowercase "g7"; the decoder requires exact-case
    // "G1"/"G6"/"G7"/"G8" (confirmed against every other fixture in this crate's
    // perturbation/error_budget test suites, which all use uppercase, and against
    // `docs/SOLVE_JSON_V1.md`) and rejects lowercase with `InvalidValue` -- fixed to uppercase
    // here rather than reproducing a decode failure the brief did not anticipate (the same fix
    // Task 7's report records making to its own brief-supplied fixtures).
    fn resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 600.0},
            "atmosphere": {}, "wind": {"speed_mps": 3.0,
                                       "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 10.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    /// A vacuum (drag-free), flat, unzeroed shot: `drop(v) = 0.5 * g * (x / v)^2` exactly (same
    /// trick `central_difference`'s own analytic oracle test uses,
    /// `src/perturbation/derive.rs`), and no wind at all so `windage_m` is exactly `0.0`
    /// everywhere. No `zero_distance_m` at all, so `with_axis`'s rezero-clearing never fires for
    /// `MuzzleVelocityMps` (a `requires_rezero` axis) -- the trajectory really is the closed
    /// form, not an approximation of one perturbed by re-zero search noise.
    fn vacuum_resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
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
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    // ---- Step 1 tests, verbatim from the brief ----

    /// Acceptance criterion: a bigger target can never shrink a one-variable bound.
    #[test]
    fn a_larger_target_never_shrinks_a_bound() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let small = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.2, height_m: 0.3 }, &domains).unwrap();
        let big = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.6, height_m: 0.9 }, &domains).unwrap();
        let s = &small.axes[0];
        let b = &big.axes[0];
        if let (Some(sf), Some(bf)) = (s.far_bound, b.far_bound) {
            assert!(bf >= sf - 1e-6, "larger target shrank the bound: {sf} -> {bf}");
        }
    }

    /// A bound that does not exist in the domain is reported as such, not invented.
    ///
    /// Review fix: also confirms `near_has_no_effect`/`far_has_no_effect` are BOTH `false` here
    /// -- WindSpeed genuinely does move the impact over this domain (see the probed sweep in
    /// `windage_dominant_axis_reports_left_or_right_not_top_or_bottom`'s doc), it just never
    /// moves it far enough to leave this deliberately huge target. Contrast
    /// `target_distance_axis_shows_no_measurable_effect_not_a_generic_unbounded_claim`, where the
    /// SAME `unbounded_in_domain: true` shape comes from an axis that has genuinely zero effect
    /// -- these two tests together are the discriminating pair for that new field.
    #[test]
    fn an_axis_that_never_exits_is_flagged_not_bounded() {
        let r = resolved();
        // A huge target cannot be missed by a small wind change.
        let domains = [(InputAxis::WindSpeed, (2.9_f64, 3.1_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 50.0, height_m: 50.0 }, &domains).unwrap();
        assert!(rep.axes[0].unbounded_in_domain);
        assert!(rep.axes[0].near_bound.is_none() && rep.axes[0].far_bound.is_none());
        assert!(
            !rep.axes[0].near_has_no_effect && !rep.axes[0].far_has_no_effect,
            "WindSpeed does move the impact within this domain -- it just never leaves this \
             huge target -- so has_no_effect must be false in both directions"
        );
    }

    /// Review fix #3: the ticket's own flagship example is a range-estimate question ("how far
    /// off could my rangefinder reading be"), and the axis that superficially matches --
    /// `TargetDistance` (`shot.max_range_m`) -- does not answer it: it is `requires_rezero:
    /// false`, so perturbing it never changes the muzzle angle or any sight correction, and
    /// therefore never changes the impact observed at a FIXED `range_m` at all (as long as the
    /// perturbed value stays >= `range_m`). Both directions here report `unbounded_in_domain:
    /// true`, and this test's whole point is that `near_has_no_effect`/`far_has_no_effect` being
    /// ALSO `true` is what tells a caller this is "provably irrelevant," not "any range error is
    /// safe." `ZeroDistance` (`shot.zero_distance_m`, `requires_rezero: true`) is the axis that
    /// actually re-zeroes for a different assumed distance and DOES move the impact at the true,
    /// fixed `range_m` -- checked here for contrast, with real, non-`None` bounds.
    #[test]
    fn target_distance_axis_shows_no_measurable_effect_not_a_generic_unbounded_claim() {
        let r = resolved(); // max_range_m = 900.0, zero_distance_m = 600.0
        let target = TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.5 };

        let td_domains = [(InputAxis::TargetDistance, (500.0_f64, 1100.0_f64))];
        let td = tolerance_envelope(&r, &[InputAxis::TargetDistance], 400.0, target, &td_domains)
            .unwrap();
        let a = &td.axes[0];
        assert!(a.nominal_inside_target);
        assert!(a.near_bound.is_none() && a.far_bound.is_none());
        assert!(a.unbounded_in_domain);
        assert!(
            a.near_has_no_effect && a.far_has_no_effect,
            "TargetDistance must show NO measurable effect in either direction, not merely an \
             unbounded one -- it cannot change the impact observed at a fixed range_m at all"
        );

        // Contrast: ZeroDistance, over the SAME true observation range, DOES move the impact,
        // and finds real bounds -- proving `has_no_effect` genuinely discriminates rather than
        // being true for every `requires_rezero: false`-adjacent axis or every wide domain.
        let zd_domains = [(InputAxis::ZeroDistance, (400.0_f64, 800.0_f64))];
        let zd = tolerance_envelope(&r, &[InputAxis::ZeroDistance], 400.0, target, &zd_domains)
            .unwrap();
        let b = &zd.axes[0];
        assert!(b.nominal_inside_target);
        assert!(
            b.near_bound.is_some() && b.far_bound.is_some(),
            "ZeroDistance must produce real bounds: re-zeroing for a different assumed distance \
             DOES move the impact observed at the true, fixed range_m"
        );
        assert!(!b.unbounded_in_domain);
    }

    #[test]
    fn the_report_refuses_to_imply_probability() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.25 }, &domains).unwrap();
        assert_eq!(rep.method, "one_variable_deterministic_bisection");
        assert!(rep.assumptions.iter().any(|s| s.contains("probability")));
        assert!(rep.assumptions.iter().any(|s| s.contains("simultaneously")));
    }

    // ---- Beyond the brief: the properties named in the task instructions ----

    /// All four correctness requirements must be stated in the payload itself, not only in
    /// prose documentation -- extends the brief's own probability/simultaneity check (which this
    /// duplicates for a different fixture) to the two the brief's test does not cover: never
    /// extrapolating past the configured domain, and confirming the nominal reads as inside
    /// before any bound is searched for.
    ///
    /// Review fix: the fourth check previously asserted only `contains("nominal")`, which
    /// `assumptions[0]` ("...stays at its **nominal** value...") already satisfies on its own --
    /// deleting `assumptions[3]` entirely left this test green. Pinned per-index instead, on a
    /// substring unique to each sentence (`"nominal_inside_target"`, the literal field name,
    /// appears ONLY in `assumptions[3]`), plus an exact length so a deleted or reordered
    /// assumption is caught even if its specific substring happened to still appear elsewhere.
    #[test]
    fn assumptions_cover_all_four_correctness_requirements() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.2, height_m: 0.3 }, &domains).unwrap();
        assert_eq!(rep.assumptions.len(), 4, "exactly four assumption sentences are expected");
        assert!(
            rep.assumptions[0].contains("simultaneously"),
            "one-variable-at-a-time must be stated"
        );
        assert!(rep.assumptions[1].contains("probability"), "no probability must be stated");
        assert!(
            rep.assumptions[2].to_lowercase().contains("domain"),
            "never-extrapolate-beyond-domain must be stated"
        );
        assert!(
            rep.assumptions[3].contains("nominal_inside_target"),
            "the nominal-inside precondition must be stated, discriminated from assumptions[0]'s \
             own unrelated use of the word \"nominal\""
        );
    }

    /// THE central distinction this ticket exists to enforce: "stays inside throughout" and
    /// "the nominal itself is not inside" produce the IDENTICAL `near_bound`/`far_bound` shape
    /// (`None`/`None`) from two completely different root causes. A degenerate (zero-area)
    /// target can never contain any point, not even the exact nominal impact
    /// (`dy == dz == 0.0`) -- contrast directly against `an_axis_that_never_exits_is_flagged_not_bounded`
    /// above, which reaches the same `None`/`None` shape via the OPPOSITE fact (a huge target
    /// that is never exited).
    #[test]
    fn a_degenerate_target_is_flagged_nominal_outside_not_confused_with_unbounded() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.0, height_m: 0.3 }, &domains).unwrap();
        let a = &rep.axes[0];
        assert!(
            !a.nominal_inside_target,
            "a zero-area target can never contain even the nominal impact"
        );
        assert!(
            a.near_bound.is_none() && a.far_bound.is_none(),
            "bounds must not be fabricated when the nominal itself is not inside"
        );
        assert!(
            !a.unbounded_in_domain,
            "this is NOT the same fact as 'stays inside throughout' -- conflating the two is \
             the exact failure this ticket exists to prevent"
        );
    }

    /// Independent oracle #1 (analytic): against the vacuum closed form, BOTH the near and far
    /// bound have an exact root, checked to a tight relative tolerance, with the domain's own
    /// midpoint (900.0) chosen far from either root so a "return the midpoint" bug -- which
    /// would still satisfy a bare `lo < x < hi` check -- is caught (mirrors
    /// `bisect_axis_converges_to_the_vacuum_analytic_root_not_the_domain_midpoint`,
    /// `src/perturbation/derive.rs`).
    #[test]
    fn bisection_bounds_match_the_vacuum_analytic_crossing() {
        let r = vacuum_resolved();
        let v0 = r.rifle.muzzle_velocity_mps;
        let x = 400.0_f64;
        const G: f64 = 9.80665;
        let drop = |v: f64| 0.5 * G * (x / v) * (x / v);
        let drop0 = drop(v0);

        let half_height = 0.1_f64;
        let domains = [(InputAxis::MuzzleVelocityMps, (700.0_f64, 1100.0_f64))];
        let rep = tolerance_envelope(
            &r,
            &[InputAxis::MuzzleVelocityMps],
            x,
            TargetGeometryV1::Rect { width_m: 1.0e6, height_m: 2.0 * half_height },
            &domains,
        )
        .unwrap();
        let a = &rep.axes[0];
        assert!(a.nominal_inside_target);
        // Review minor #5: every other test in this module happens to use range_m = 600.0, so
        // `range_m` echoing the caller's input was never checked against any OTHER value --
        // this fixture's own 400.0 closes that gap.
        assert_eq!(rep.range_m, x);

        // far (v > v0): drop DECREASES, crossing where drop(v) = drop0 - half_height.
        let expected_far = x / (2.0 * (drop0 - half_height) / G).sqrt();
        // near (v < v0): drop INCREASES, crossing where drop(v) = drop0 + half_height.
        let expected_near = x / (2.0 * (drop0 + half_height) / G).sqrt();

        let far = a.far_bound.expect("a crossing exists well within (700, 1100)");
        let near = a.near_bound.expect("a crossing exists well within (700, 1100)");
        let rel_far = ((far - expected_far) / expected_far).abs();
        let rel_near = ((near - expected_near) / expected_near).abs();
        assert!(rel_far < 0.02, "far: expected ~{expected_far}, got {far} (rel {rel_far})");
        assert!(rel_near < 0.02, "near: expected ~{expected_near}, got {near} (rel {rel_near})");
        assert!((far - 900.0).abs() > 30.0, "far bound must not be the domain midpoint");
        assert!((near - 900.0).abs() > 30.0, "near bound must not be the domain midpoint");

        assert_eq!(a.far_limiting_boundary, Some(LimitingBoundaryV1::Top));
        assert_eq!(a.near_limiting_boundary, Some(LimitingBoundaryV1::Bottom));
    }

    /// Review fix #7 (part 1): pins `bisection_tolerance`'s formula directly, independent of any
    /// downstream bisection -- catches a hardcoded constant or a different scaling factor.
    #[test]
    fn bisection_tolerance_scales_with_domain_width_not_a_flat_constant() {
        assert_eq!(bisection_tolerance(0.0, 20.0), 20.0 * 1e-6);
        assert_eq!(bisection_tolerance(700.0, 1100.0), 400.0 * 1e-6);
        // Floor: a domain so narrow that width * 1e-6 would underflow below what 80 bisection
        // iterations could usefully resolve is clamped to 1e-9, not left arbitrarily small.
        assert_eq!(bisection_tolerance(800.0, 800.0 + 1e-5), 1e-9);
    }

    /// Review fix #7 (part 2): the domain-proportional formula is undefended end-to-end by every
    /// OTHER test in this module, because every domain they use is wide enough that even a flat
    /// `1e-4` tolerance would still run several real bisection iterations and land close to the
    /// true crossing -- the two formulas are indistinguishable from the outside on a wide domain.
    /// This test uses a domain SO narrow (`5e-5` per direction) that a flat `1e-4` tolerance
    /// would satisfy `bisect_axis`'s very FIRST check (`(hi - lo).abs() <= tolerance`) with ZERO
    /// bisection steps, returning the sub-domain's own exact algebraic midpoint verbatim --
    /// verified by literally making that mutation (`bisection_tolerance` hardcoded to `1e-4`),
    /// which reproducibly returns `far_bound` exactly equal to `800.000025` (`diff = 0e0`, not
    /// even one ULP off) -- then reverting. Under the shipped domain-proportional formula
    /// (`tol = 0.0001 * 1e-6`, floored to `1e-9`), real bisection runs and the result differs
    /// from that same exact midpoint by `~8.6e-6` -- more than four orders of magnitude past any
    /// floating-point noise floor, so `> 1e-7` below is a comfortable, non-flaky margin.
    #[test]
    fn a_narrow_domain_gets_real_bisection_not_an_immediate_midpoint_return() {
        let r = vacuum_resolved();
        let x = 400.0_f64;
        let half_height = 5e-8_f64;
        let domains = [(InputAxis::MuzzleVelocityMps, (799.99995_f64, 800.00005_f64))];
        let rep = tolerance_envelope(
            &r,
            &[InputAxis::MuzzleVelocityMps],
            x,
            TargetGeometryV1::Rect { width_m: 1.0e6, height_m: 2.0 * half_height },
            &domains,
        )
        .unwrap();
        let a = &rep.axes[0];
        let far_midpoint = 0.5 * (800.0_f64 + 800.00005_f64);
        let far = a.far_bound.expect("a crossing must exist this close to nominal");
        assert!(
            (far - far_midpoint).abs() > 1e-7,
            "far_bound ({far}) must differ meaningfully from the domain's own exact midpoint \
             ({far_midpoint}) -- an unchanged difference of ~0 would mean bisect_axis returned \
             the midpoint verbatim without ever refining it, exactly what a flat 1e-4 tolerance \
             does on a domain this narrow"
        );
    }

    /// Independent oracle #2 (re-solve): a found bound is verified by re-solving at it, and at
    /// points just inside and just outside it, through a path that does NOT go through
    /// `with_axis`/`evaluate`/`bisect_axis` at all -- a hand-built request decoded and solved via
    /// `solve_v1` directly (mirroring `mod.rs`'s own
    /// `evaluate_matches_solve_v1_when_zero_distance_alone_searches_the_elevation` cross-check
    /// pattern), reading the sample off `solve_v1`'s OWN wire output at an exact grid point
    /// (600 m, an exact multiple of this fixture's 10 m sampling interval).
    #[test]
    fn a_found_bound_sits_on_the_target_boundary_verified_independently_of_with_axis() {
        let r = resolved();
        let target = TargetGeometryV1::Rect { width_m: 0.2, height_m: 0.3 };
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0, target, &domains).unwrap();
        let a = &rep.axes[0];
        assert!(a.nominal_inside_target);
        let far = a.far_bound.expect("this target/domain combination must produce a far bound");

        fn independent_deviation(wind_speed: f64) -> (f64, f64) {
            let json = serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0, "zero_distance_m": 600.0},
                "atmosphere": {},
                "wind": {"speed_mps": wind_speed,
                         "direction_from_rad": std::f64::consts::FRAC_PI_2},
                "solver": {}, "effects": {}, "sampling": {"interval_m": 10.0}
            })
            .to_string();
            let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
            let solved = crate::solve_v1::solve_v1(req).unwrap();
            let sample = solved
                .samples
                .iter()
                .find(|s| s.distance_m == 600.0)
                .expect("600 m must be an exact grid point for this fixture");
            (sample.drop_m, sample.windage_m)
        }

        let half_width = 0.2 / 2.0;
        let half_height = 0.3 / 2.0;
        let (nominal_drop, nominal_windage) = independent_deviation(3.0);
        let ratio_at = |wind_speed: f64| -> f64 {
            let (d, w) = independent_deviation(wind_speed);
            let ry = (d - nominal_drop).abs() / half_height;
            let rz = (w - nominal_windage).abs() / half_width;
            ry.max(rz)
        };

        let ratio_far = ratio_at(far);
        assert!(
            (ratio_far - 1.0).abs() < 1e-3,
            "found bound does not sit on the target boundary independently: ratio {ratio_far}"
        );

        let just_inside = 3.0 + (far - 3.0) * 0.99;
        let just_outside = 3.0 + (far - 3.0) * 1.01;
        assert!(
            ratio_at(just_inside) < 1.0,
            "a point just inside the found bound must independently read as inside the target"
        );
        assert!(
            ratio_at(just_outside) > 1.0,
            "a point just past the found bound must independently read as outside the target"
        );
    }

    /// The brief's own acceptance criterion, generalized from two points to a genuine sweep:
    /// growing the target can never shrink a one-variable bound, checked across many sizes,
    /// consecutively -- two isolated points cannot catch a bound that stops being monotone only
    /// partway through a range. Also asserts that once a direction becomes unbounded (no
    /// crossing found within the domain) at some size, it stays unbounded at every LARGER size:
    /// the literal meaning of "the crossing moved even farther away, past the domain edge."
    /// WindSpeed's domain `(0, 20)` around a nominal of `3.0` is asymmetric (17 units of room
    /// above, only 3 below), so the near direction saturates (becomes unbounded) at a much
    /// smaller scale than the far direction -- the geometric sweep below is wide enough to
    /// observe both transitions, confirmed by the closing assertion that both actually occurred.
    #[test]
    fn a_larger_target_never_shrinks_a_bound_across_a_size_sweep() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let scales: Vec<f64> = (0..12).map(|i| 0.5 * 1.6_f64.powi(i)).collect();

        let mut prev_far: Option<f64> = None;
        let mut prev_near: Option<f64> = None;
        let mut far_became_unbounded = false;
        let mut near_became_unbounded = false;

        for &scale in &scales {
            let target =
                TargetGeometryV1::Rect { width_m: 0.2 * scale, height_m: 0.3 * scale };
            let rep =
                tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0, target, &domains).unwrap();
            let a = &rep.axes[0];
            assert!(
                a.nominal_inside_target,
                "scale {scale}: nominal must read as inside for any positive-size target"
            );

            if far_became_unbounded {
                assert!(
                    a.far_bound.is_none(),
                    "scale {scale}: far bound reappeared after becoming unbounded at a smaller \
                     scale"
                );
            } else if let (Some(pf), Some(f)) = (prev_far, a.far_bound) {
                assert!(f >= pf - 1e-6, "far bound shrank at scale {scale}: {pf} -> {f}");
            }
            if a.far_bound.is_none() {
                far_became_unbounded = true;
            }
            prev_far = a.far_bound;

            if near_became_unbounded {
                assert!(
                    a.near_bound.is_none(),
                    "scale {scale}: near bound reappeared after becoming unbounded at a smaller \
                     scale"
                );
            } else if let (Some(pn), Some(n)) = (prev_near, a.near_bound) {
                assert!(
                    n <= pn + 1e-6,
                    "near bound shrank (moved toward nominal) at scale {scale}: {pn} -> {n}"
                );
            }
            if a.near_bound.is_none() {
                near_became_unbounded = true;
            }
            prev_near = a.near_bound;
        }

        assert!(
            far_became_unbounded,
            "sweep never reached an unbounded far regime -- the monotonicity check on that \
             transition never actually ran"
        );
        assert!(
            near_became_unbounded,
            "sweep never reached an unbounded near regime -- the monotonicity check on that \
             transition never actually ran"
        );
        // Sanity: the sweep must also have exercised at least one genuinely bounded pair on each
        // side, or the "shrank" assertions above would never fire either.
        assert!(scales[0] < 1.0, "sweep must start comfortably inside the bounded regime");
    }

    /// Same acceptance criterion, briefly, for the OTHER `TargetGeometryV1` shape -- a circle's
    /// single `radius_m` grows strictly with the same scale factor, so this is a smaller sweep
    /// than the rectangle's (which independently varies two dimensions), not a second full
    /// property test.
    #[test]
    fn a_larger_circle_never_shrinks_a_bound_across_a_size_sweep() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let scales = [0.5_f64, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
        let mut prev_far: Option<f64> = None;
        let mut far_became_unbounded = false;
        for &scale in &scales {
            let target = TargetGeometryV1::Circle { radius_m: 0.15 * scale };
            let rep =
                tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0, target, &domains).unwrap();
            let a = &rep.axes[0];
            if far_became_unbounded {
                assert!(a.far_bound.is_none(), "scale {scale}: far bound reappeared");
            } else if let (Some(pf), Some(f)) = (prev_far, a.far_bound) {
                assert!(f >= pf - 1e-6, "far bound shrank at scale {scale}: {pf} -> {f}");
            }
            if a.far_bound.is_none() {
                far_became_unbounded = true;
            }
            prev_far = a.far_bound;
        }
        assert!(far_became_unbounded, "circle sweep never reached an unbounded far regime");
    }

    /// `axis`/`nominal` must each be tied to the axis they describe, not swapped or copied from
    /// a sibling entry -- WindSpeed (nominal 3.0) and MuzzleVelocityMps (nominal 823.0) are far
    /// enough apart that a transposition is unmistakable. `margin_linear_m` is independently
    /// pinned as a TARGET-derived constant identical across both axes (`min(0.6, 0.4) / 2 =
    /// 0.2`), not a copy of `nominal` or any other axis-specific field.
    #[test]
    fn multiple_axes_are_each_independently_computed_and_correctly_tagged() {
        let r = resolved();
        let domains = [
            (InputAxis::WindSpeed, (0.0_f64, 20.0_f64)),
            (InputAxis::MuzzleVelocityMps, (400.0_f64, 1200.0_f64)),
        ];
        let rep = tolerance_envelope(
            &r,
            &[InputAxis::WindSpeed, InputAxis::MuzzleVelocityMps],
            600.0,
            TargetGeometryV1::Rect { width_m: 0.4, height_m: 0.6 },
            &domains,
        )
        .unwrap();
        assert_eq!(rep.axes.len(), 2);
        assert_eq!(rep.axes[0].axis, InputAxis::WindSpeed);
        assert_eq!(rep.axes[1].axis, InputAxis::MuzzleVelocityMps);

        let expected_ws = match read_axis(&r, InputAxis::WindSpeed).unwrap() {
            AxisValue::Scalar(x) => x,
            other => panic!("WindSpeed must read back as a scalar, got {other:?}"),
        };
        let expected_mv = match read_axis(&r, InputAxis::MuzzleVelocityMps).unwrap() {
            AxisValue::Scalar(x) => x,
            other => panic!("MuzzleVelocityMps must read back as a scalar, got {other:?}"),
        };
        assert_eq!(rep.axes[0].nominal, expected_ws);
        assert_eq!(rep.axes[1].nominal, expected_mv);
        assert_ne!(
            rep.axes[0].nominal, rep.axes[1].nominal,
            "fixture sanity: the two axes must have DIFFERENT nominal values or a \
             transposition between them would be invisible"
        );
        assert!(rep.axes[0].nominal_inside_target);
        assert!(rep.axes[1].nominal_inside_target);
        assert!((rep.axes[0].margin_linear_m - 0.2).abs() < 1e-9);
        assert!((rep.axes[1].margin_linear_m - 0.2).abs() < 1e-9);
        assert_eq!(rep.range_m, 600.0);
    }

    /// `near_limiting_boundary`/`far_limiting_boundary` must independently discriminate Left/
    /// Right from Top/Bottom, not default to one or the other. A pure crosswind's effect on
    /// windage dominates its (tiny, secondary) effect on drop by many orders of magnitude (see
    /// `windage_derivative_wrt_wind_speed_dominates_and_matches_the_baseline_sign`,
    /// `src/perturbation/derive.rs`), so any WindSpeed crossing here must be windage-driven: more
    /// wind pushes windage more negative (see this module's own probed sweep), so the far
    /// direction (more wind) exits `Left` and the near direction (less wind) exits `Right`.
    #[test]
    fn windage_dominant_axis_reports_left_or_right_not_top_or_bottom() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.2, height_m: 0.3 }, &domains).unwrap();
        let a = &rep.axes[0];
        assert!(a.far_bound.is_some() && a.near_bound.is_some());
        assert_eq!(a.far_limiting_boundary, Some(LimitingBoundaryV1::Left));
        assert_eq!(a.near_limiting_boundary, Some(LimitingBoundaryV1::Right));
    }

    /// The complementary case for the SAME field: `bisection_bounds_match_the_vacuum_analytic_crossing`
    /// above already asserts `Top`/`Bottom` for a wind-free vacuum fixture where `windage_m` is
    /// exactly `0.0` throughout, so a crossing driven by muzzle velocity alone (which barely
    /// moves windage at all -- see `d_windage_d_x` in the same vacuum context) cannot be
    /// misclassified as `Left`/`Right`. Restated here as its own assertion, independent of the
    /// closed-form numeric check, so a mutation that hardcodes `Radial`/`Left` for every `Rect`
    /// would be caught even if the numeric oracle above were skipped.
    #[test]
    fn drop_dominant_axis_reports_top_or_bottom_not_left_or_right() {
        let r = vacuum_resolved();
        let domains = [(InputAxis::MuzzleVelocityMps, (700.0_f64, 1100.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::MuzzleVelocityMps], 400.0,
            TargetGeometryV1::Rect { width_m: 1.0e6, height_m: 0.2 }, &domains).unwrap();
        let a = &rep.axes[0];
        assert!(a.far_bound.is_some() && a.near_bound.is_some());
        assert_eq!(a.far_limiting_boundary, Some(LimitingBoundaryV1::Top));
        assert_eq!(a.near_limiting_boundary, Some(LimitingBoundaryV1::Bottom));
    }

    /// A [`TargetGeometryV1::Circle`] has no distinct edges and must always report `Radial`,
    /// regardless of which direction (windage- or drop-dominant) the crossing came from. Also
    /// pins `target_margin_linear_m`'s `Circle` arm (review minor #6): previously only the
    /// `Rect` arm was checked numerically (`multiple_axes_are_each_independently_computed_...`),
    /// so a bug specific to the `Circle` branch (e.g. returning `radius_m` un-halved, or a
    /// hardcoded `0.0`) could have passed every existing test.
    #[test]
    fn a_circle_target_always_reports_radial() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.15 }, &domains).unwrap();
        let a = &rep.axes[0];
        assert!(a.far_bound.is_some() && a.near_bound.is_some());
        assert_eq!(a.far_limiting_boundary, Some(LimitingBoundaryV1::Radial));
        assert_eq!(a.near_limiting_boundary, Some(LimitingBoundaryV1::Radial));
        assert!(
            (a.margin_linear_m - 0.15).abs() < 1e-9,
            "Circle's margin_linear_m must be exactly its radius_m, got {}",
            a.margin_linear_m
        );
    }

    /// Review fix #2: every existing fixture that reaches `limiting_boundary` has one deviation
    /// at (or indistinguishable from) zero -- the windage-dominant fixture has `dy ~ 0`, both
    /// drop-dominant fixtures use a wind-free vacuum with `dz` exactly `0`. Under any of those,
    /// swapping the cross-multiplication's `width_m`/`height_m` (`dy.abs() * height_m >=
    /// dz.abs() * width_m` instead of the correct `dy.abs() * width_m >= dz.abs() * height_m`)
    /// still passes, because one side of the comparison collapses to (approximately) zero either
    /// way. This calls `limiting_boundary` directly (it is private, and this test lives in the
    /// same module) with both deviations simultaneously non-trivial and an aspect-ratio-skewed
    /// target, so the two formulas disagree: `dy = 0.1` against `height_m = 0.2` (ratio `1.0`)
    /// and `dz = 0.5` against `width_m = 2.0` (ratio `0.5`) -- correct: `ratio_y > ratio_z` =>
    /// `Bottom`; swapped: `dy.abs()*height_m = 0.02 < dz.abs()*width_m = 1.0` => `Left`/`Right`
    /// branch entirely, i.e. a DIFFERENT edge family, not just a different edge.
    #[test]
    fn limiting_boundary_uses_the_correct_axis_for_each_deviation_not_swapped() {
        let nominal = Observation {
            range_m: 600.0,
            drop_m: 0.0,
            windage_m: 0.0,
            time_s: 1.0,
            velocity_mps: 500.0,
        };
        let o = Observation {
            range_m: 600.0,
            drop_m: 0.1,
            windage_m: 0.5,
            time_s: 1.0,
            velocity_mps: 500.0,
        };
        let target = TargetGeometryV1::Rect { width_m: 2.0, height_m: 0.2 };
        assert_eq!(
            limiting_boundary(&o, &nominal, target),
            LimitingBoundaryV1::Bottom,
            "dy/height_m ratio (1.0) exceeds dz/width_m ratio (0.5): must be a drop-family edge"
        );
        // Mirror on the drop side (dy < 0) to also confirm Top is reachable from this same
        // aspect-ratio-skewed target, not just Bottom.
        let o_top = Observation { drop_m: -0.1, ..o };
        assert_eq!(limiting_boundary(&o_top, &nominal, target), LimitingBoundaryV1::Top);
    }

    /// A categorical axis has no numeric domain to bisect -- must be recorded, not silently
    /// dropped from the report and not a hard failure of the whole call.
    #[test]
    fn a_categorical_axis_is_recorded_unavailable_not_silently_dropped_or_hard_failed() {
        let r = resolved();
        let rep = tolerance_envelope(&r, &[InputAxis::CoriolisEnabled], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.25 }, &[]).unwrap();
        assert!(rep.axes.is_empty(), "a categorical axis must never appear in `axes`");
        assert_eq!(rep.unavailable_axes.len(), 1);
        assert_eq!(rep.unavailable_axes[0].axis, InputAxis::CoriolisEnabled);
        assert_eq!(rep.unavailable_axes[0].code, UnavailableReasonCodeV1::CategoricalAxis);
        assert!(!rep.unavailable_axes[0].reason.is_empty());
    }

    /// The three wind axes have no single scalar value under segmented wind -- `read_axis`
    /// returns `None` directly (never even reaching `bisect_axis`), and that must be recorded
    /// too, with no `domains` entry required for an axis that never gets that far.
    #[test]
    fn a_wind_axis_under_segmented_wind_is_recorded_unavailable() {
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
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 300.0,
            TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.5 }, &[]).unwrap();
        assert!(rep.axes.is_empty());
        assert_eq!(rep.unavailable_axes.len(), 1);
        assert_eq!(rep.unavailable_axes[0].axis, InputAxis::WindSpeed);
        assert_eq!(rep.unavailable_axes[0].code, UnavailableReasonCodeV1::AxisAbsent);
    }

    /// `Altitude` under a QNH-referenced atmosphere is refused by `with_axis` itself
    /// (`AxisUnsupportedForRequest`) -- reached through this module's pre-flight step, not a
    /// hard failure of the whole report.
    #[test]
    fn altitude_under_qnh_is_recorded_unavailable_not_hard_failed() {
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
        let domains = [(InputAxis::Altitude, (0.0_f64, 1000.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::Altitude], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap();
        assert!(rep.axes.is_empty());
        assert_eq!(rep.unavailable_axes.len(), 1);
        assert_eq!(rep.unavailable_axes[0].axis, InputAxis::Altitude);
        assert_eq!(
            rep.unavailable_axes[0].code,
            UnavailableReasonCodeV1::AxisUnsupportedForRequest
        );
        assert!(rep.unavailable_axes[0].reason.to_lowercase().contains("qnh"));
    }

    /// A domain whose lower bound, for `TargetDistance`, dips below the caller's OWN `range_m`
    /// is an internally inconsistent request (bisecting "how close could the target be" past
    /// "the range I am asking about"): `bisect_axis` hits a genuine `Observation::OutOfRange`
    /// partway through the near-direction search, which must propagate as a hard failure of the
    /// whole call, never be swallowed into `unavailable_axes` (it is not one of the four
    /// structural refusals).
    #[test]
    fn a_genuine_observation_error_propagates_not_recorded_as_unavailable() {
        let r = resolved(); // max_range_m = 900.0, zero_distance_m = 600.0
        let domains = [(InputAxis::TargetDistance, (500.0_f64, 1100.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::TargetDistance], 600.0,
            TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.5 }, &domains).unwrap_err();
        match err {
            KernelError::Observation(TrajectoryObservationError::OutOfRange { requested_m, .. }) => {
                assert_eq!(requested_m, 600.0);
            }
            other => panic!("expected Observation(OutOfRange {{ .. }}), got {other:?}"),
        }
    }

    /// `domains` must supply an entry for every requested continuous, present axis -- no
    /// implicit fallback (see the module doc's "Domains are validated up front" section).
    #[test]
    fn a_missing_domain_is_reported_as_invalid_not_defaulted() {
        let r = resolved();
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &[]).unwrap_err();
        match err {
            KernelError::InvalidDomain { axis, .. } => assert_eq!(axis, InputAxis::WindSpeed),
            other => panic!("expected InvalidDomain, got {other:?}"),
        }
    }

    /// The axis's own current value (`3.0`) must sit strictly inside the configured domain, or
    /// one search direction degenerates to a zero-width probe.
    #[test]
    fn nominal_outside_the_configured_domain_is_rejected() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (5.0_f64, 20.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap_err();
        assert!(matches!(err, KernelError::InvalidDomain { axis: InputAxis::WindSpeed, .. }));
    }

    /// An inverted or zero-width domain (`lo >= hi`) is rejected outright.
    #[test]
    fn an_inverted_domain_is_rejected() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (20.0_f64, 0.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap_err();
        assert!(matches!(err, KernelError::InvalidDomain { axis: InputAxis::WindSpeed, .. }));
    }

    /// Review fix (most important finding): `nominal_outside_the_configured_domain_is_rejected`
    /// puts the nominal (`3.0`) OUTSIDE its domain `(5.0, 20.0)` entirely, and
    /// `an_inverted_domain_is_rejected` is caught by the EARLIER `lo >= hi` check before the
    /// nominal-strictness check ever runs -- neither test can tell a strict `>`/`<` comparison
    /// apart from a relaxed `>=`/`<=` one. Relaxing `tolerance.rs`'s validation to `>=`/`<=`
    /// leaves ALL other tests in this module green: with `domains = [(WindSpeed, (3.0, 20.0))]`
    /// and nominal `3.0`, `bisect_axis` would be called with `(nominal_value, lo) = (3.0, 3.0)`
    /// -- a zero-width probe whose two endpoints trivially agree, so it would report `Ok(None)`
    /// for a direction that was never actually searched even one step into, and (if the far
    /// direction also found nothing) `unbounded_in_domain: true` -- "your wind call is robust in
    /// every direction," fabricated from a search that never moved. This test, and its mirror
    /// below for the upper edge, place the nominal EXACTLY at one edge and require
    /// `InvalidDomain` -- these fail under the `>=`/`<=` relaxation described above (verified by
    /// making that exact edit, confirming exactly these two tests fail with all others still
    /// green, then reverting to byte-identical -- see the task report).
    #[test]
    fn nominal_at_the_lower_domain_edge_is_rejected() {
        let r = resolved(); // WindSpeed nominal = 3.0
        let domains = [(InputAxis::WindSpeed, (3.0_f64, 20.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap_err();
        assert!(matches!(err, KernelError::InvalidDomain { axis: InputAxis::WindSpeed, .. }));
    }

    /// Mirror of the above at the UPPER edge -- see that test's doc for the full rationale.
    #[test]
    fn nominal_at_the_upper_domain_edge_is_rejected() {
        let r = resolved(); // WindSpeed nominal = 3.0
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 3.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap_err();
        assert!(matches!(err, KernelError::InvalidDomain { axis: InputAxis::WindSpeed, .. }));
    }

    /// `range_m` itself must be a queryable point on the base trajectory, exactly as
    /// `error_budget` requires of its own `ranges_m` -- checked directly, before any per-axis
    /// work, rather than surfacing as a mysterious per-axis refusal.
    #[test]
    fn an_out_of_range_query_is_rejected_directly() {
        let r = resolved(); // max_range_m = 900.0
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let err = tolerance_envelope(&r, &[InputAxis::WindSpeed], 5000.0,
            TargetGeometryV1::Circle { radius_m: 0.3 }, &domains).unwrap_err();
        match err {
            KernelError::Observation(TrajectoryObservationError::OutOfRange { requested_m, .. }) => {
                assert_eq!(requested_m, 5000.0);
            }
            other => panic!("expected Observation(OutOfRange {{ .. }}), got {other:?}"),
        }
    }

    /// Every public field on [`ToleranceReportV1`]/[`ToleranceAxisV1`]/[`UnavailableAxisV1`]
    /// round-trips through JSON, and the externally-tagged `snake_case` axis name is exactly
    /// what a consumer of this wire type would expect -- pins the `#[serde(rename_all =
    /// "snake_case")]` convention this crate uses throughout, not just that serialization
    /// succeeds at all.
    ///
    /// Compares fields individually with a `1e-9` tolerance on the bisection-derived `f64`s
    /// rather than a whole-struct `assert_eq!`: `serde_json` 1.0.149's float parser round-trips
    /// SOME specific `f64` bit patterns (arbitrary bisection results among them) one ULP off --
    /// serializing produces the shortest correctly-round-tripping decimal string, but re-parsing
    /// that exact string yields a different bit pattern one ULP away. Confirmed with a
    /// standalone reproduction outside this crate on the exact bit pattern this test's own
    /// `near_bound` produced (`serde_json::to_string`/`from_str` on a bare `f64`), matching the
    /// identical upstream characteristic `error_budget.rs`'s
    /// `p_hit_and_gain_round_trip_through_json_when_present` already documents -- not a defect
    /// in this module's `Serialize`/`Deserialize` derives.
    #[test]
    fn the_report_round_trips_through_json() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Rect { width_m: 0.2, height_m: 0.3 }, &domains).unwrap();
        let json = serde_json::to_string(&rep).unwrap();
        assert!(json.contains("\"axis\":\"wind_speed\""));
        let back: ToleranceReportV1 = serde_json::from_str(&json).unwrap();

        assert_eq!(back.schema_version, rep.schema_version);
        assert_eq!(back.method, rep.method);
        assert_eq!(back.assumptions, rep.assumptions);
        assert_eq!(back.range_m, rep.range_m);
        assert_eq!(back.unavailable_axes, rep.unavailable_axes);
        assert_eq!(back.axes.len(), rep.axes.len());
        let (ba, ra) = (&back.axes[0], &rep.axes[0]);
        assert_eq!(ba.axis, ra.axis);
        assert_eq!(ba.nominal_inside_target, ra.nominal_inside_target);
        assert_eq!(ba.unbounded_in_domain, ra.unbounded_in_domain);
        assert_eq!(ba.near_limiting_boundary, ra.near_limiting_boundary);
        assert_eq!(ba.far_limiting_boundary, ra.far_limiting_boundary);
        assert!((ba.nominal - ra.nominal).abs() < 1e-9);
        assert!((ba.margin_linear_m - ra.margin_linear_m).abs() < 1e-9);
        assert!((ba.near_bound.unwrap() - ra.near_bound.unwrap()).abs() < 1e-9);
        assert!((ba.far_bound.unwrap() - ra.far_bound.unwrap()).abs() < 1e-9);
    }

    /// `schema_version` must be the crate's declared constant, not a stray literal that could
    /// silently drift from it.
    #[test]
    fn schema_version_matches_the_declared_constant() {
        let r = resolved();
        let domains = [(InputAxis::WindSpeed, (0.0_f64, 20.0_f64))];
        let rep = tolerance_envelope(&r, &[InputAxis::WindSpeed], 600.0,
            TargetGeometryV1::Circle { radius_m: 0.25 }, &domains).unwrap();
        assert_eq!(rep.schema_version, TOLERANCE_SCHEMA_VERSION_V1);
    }
}
