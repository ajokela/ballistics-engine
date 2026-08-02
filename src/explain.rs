//! MBA-1345: explain why two fully resolved solutions differ.
//!
//! For each input group the difference is attributed by a SYMMETRIC counterfactual: swap the
//! group from B into A, and independently from A into B, and average. That makes the answer
//! independent of replacement order, which a one-directional swap is not. Whatever the group
//! contributions fail to explain is reported as an explicit interaction remainder and is never
//! distributed across groups: for correlated inputs there is no unique causal attribution, and
//! pretending otherwise is the failure this design exists to avoid.
//!
//! # The derived-value exclusion rule
//!
//! A resolved request echoes some axes as plain, independent inputs, but a few are actually
//! COMPUTED from a different axis's value. A group swap must exclude any axis whose resolved
//! value is computed from another axis whenever that source axis is:
//!
//! (a) in a DIFFERENT group, or
//! (b) itself excluded (from the SAME group).
//!
//! Same-group derivation where the source is ALSO swapped normally is self-consistent and
//! needs no special handling -- both values move together, so the pairing the derivation
//! depends on is preserved. That is why, for instance, `Temperature`'s ICAO-standard default
//! from `Altitude` is ordinarily harmless: swapping the whole `Atmosphere` group moves both
//! (when `Altitude` is not itself excluded -- see the QNH case below, which is exactly the
//! exception).
//!
//! Three axes are known to violate this, found one at a time by successive reviews, each
//! documented below as its own subsection: `MuzzleAngle` (case (a): derived from
//! `MuzzleVelocity` and `Atmosphere`, both different groups), `WindDirection` under compass
//! wind (case (a): derived from `ShotGeometry`'s `ShotAzimuth`), and `Pressure` under a QNH
//! atmosphere (case (b): derived from `Atmosphere`'s own `Altitude`, which is excluded rather
//! than swapped). There is no taxonomy metadata recording "derived from" relationships -- each
//! is handled as an explicit, documented special case in `plan_exclusions` rather than a
//! generic mechanism, since three known instances do not justify inventing a
//! derivation-tracking data model the taxonomy does not otherwise have. A future instance
//! should be diagnosed against this same rule and added the same way.
//!
//! ## Instance 1: `MuzzleAngle` is derived, not independent, once a zero distance is present
//!
//! `ZeroSightGeometry` lists `MuzzleAngle` -- but on a RESOLVED request, `muzzle_angle_rad` is
//! not an independent input whenever `zero_distance_m` is also present: it is the
//! ALREADY-SEARCHED elevation, a function of muzzle velocity and atmosphere as much as of any
//! sight/zero knob in this group (review C1; rule case (a) -- `MuzzleVelocity` and `Atmosphere`
//! are different groups). Naively applying it as an ordinary axis would import muzzle-velocity's
//! and atmosphere's own differences into `ZeroSightGeometry`'s bucket -- concretely, swapping a
//! group in which every visible sight/zero setting is byte-identical between the two requests
//! would still report a non-zero contribution, because the group also carried the OTHER
//! request's fully-baked elevation. A shooter reading "5 cm of this is your zero/sight geometry"
//! would go check a scope that never moved.
//!
//! `taxonomy.rs`'s `axes_in_group(ZeroSightGeometry)` lists `MuzzleAngle` FIRST specifically to
//! fix this: `SightHeight` (always present, `requires_rezero`) is applied second, and its own
//! re-resolve clears and re-derives the elevation for the DESTINATION's own muzzle velocity and
//! atmosphere, overwriting whatever `MuzzleAngle`'s write introduced -- see that match arm's own
//! comment for the full mechanism, and `every_groups_contribution_matches_an_independent_recomputation`
//! below for the acceptance test (identical sight/zero inputs must report an EXACTLY zero
//! `ZeroSightGeometry`, not merely a small one). For an angle-only request (no `zero_distance_m`
//! on the destination), nothing later in the group clears the angle, so `MuzzleAngle`'s own
//! write correctly stands: there it genuinely IS the independent input, and it must still be
//! swapped in that case -- see `muzzle_angle_is_still_swapped_for_an_angle_only_request` (review
//! round 3, N3): the fix above is only correct BECAUSE the later clear-gate is conditioned on
//! `zero_distance_m` being present, not unconditionally on axis order, and there was no test
//! pinning that down until this one.
//!
//! ## Instance 2: `WindDirection` is derived from `ShotAzimuth` under compass-referenced wind
//!
//! Under `wind_reference: "compass"`, the resolved (shooter-relative) `direction_from_rad` is
//! `compass_bearing_to_shooter_relative_rad(bearing, shot_azimuth_rad)` =
//! `(bearing - shot_azimuth_rad).rem_euclid(2*pi)` (`solve_v1.rs`'s `resolve_wind`, `wind.rs`) --
//! a function of `ShotAzimuth`, which lives in `ShotGeometry`, a DIFFERENT group (review round
//! 3, N1; rule case (a)). [`with_axis`] already refuses to swap `ShotAzimuth` itself under
//! compass wind, protecting `ShotGeometry`'s own number, but nothing stopped `WindDirection`
//! from moving on ITS OWN: two compass-referenced requests with an IDENTICAL raw wind bearing
//! but different `shot_azimuth_rad` resolve to DIFFERENT `direction_from_rad` values purely
//! because of the azimuth difference, and swapping `WindDirection` would report that as a `Wind`
//! contribution -- the exact `MuzzleAngle` shape, relabelled onto a different pair of groups.
//! `plan_exclusions` now excludes `WindDirection` (both legs, both directions) whenever EITHER
//! `a` or `b` is compass-referenced, regardless of whether the raw bearing or the azimuth
//! actually differ -- the same unconditional-on-values pattern the `Altitude`/`ShotAzimuth`
//! guards already use. See `wind_direction_is_excluded_under_compass_wind_even_with_an_identical_raw_bearing`.
//!
//! ## Instance 3: `Pressure` is derived from `Altitude` under a QNH atmosphere
//!
//! Under `pressure_reference: "qnh"`, the resolved `pressure_pa` is
//! `reduce_qnh_to_station_pressure(qnh_value, altitude_m)` (`solve_v1.rs`'s `resolve_atmosphere`)
//! -- a function of `Altitude`, which lives in the SAME group, `Atmosphere` (review round 3, N2;
//! rule case (b)). Same-group derivation is ordinarily fine (see the rule above), but `Altitude`
//! itself is excluded whenever `pressure_reference` is QNH -- so the pairing breaks: `Pressure`
//! would still move while its source, `Altitude`, stays fixed on both legs, leaking a value that
//! is only physically valid at the OTHER request's (unswapped) altitude into `Atmosphere`'s
//! contribution, in violation of the third `assumptions` entry's promise that an excluded axis's
//! effect lands in the remainder, not partially in its group. `plan_exclusions` now excludes
//! `Pressure` too, using the exact same QNH condition that already excludes `Altitude`. See
//! `pressure_exclusion_does_not_leak_a_partial_effect_into_atmosphere` (a QNH-vs-QNH fixture,
//! added alongside the existing QNH-vs-absolute `altitude_exclusion_does_not_leak_a_partial_effect_into_atmosphere`,
//! which cannot exercise this specific hazard because a plain absolute pressure has no altitude
//! dependency to leak in the first place).
//!
//! **Known limitation, not fixed:** `Temperature` has the identical hazard when it is OMITTED
//! from the raw request under a QNH atmosphere -- its ICAO-standard default
//! (`calculate_icao_standard_atmosphere(altitude_m)`) is ALSO a function of `Altitude`. Unlike
//! `Pressure`'s dependency, this one is not detectable from a resolved request alone: the
//! resolved `temperature_k` echo carries only the concrete numeric value, with no record of
//! whether it came from an explicit input or from this default, and the assumption notice that
//! WOULD say so is not part of `ResolvedSolveRequestV1` at all. Excluding `Temperature`
//! unconditionally whenever `Altitude` is QNH-excluded would trade this rare, real leak for a
//! much more common one (silently discarding a genuinely independent temperature difference on
//! every ordinary QNH request that happens to also specify a temperature), which is a worse
//! trade. Left undetected and undocumented anywhere else but here.
//!
//! # Symmetric exclusion: an axis unusable on either side must be excluded on BOTH legs
//!
//! [`with_axis`] can refuse an axis for a request-specific structural reason:
//! `KernelError::AxisUnsupportedForRequest` for `Altitude` under a QNH-referenced atmosphere or
//! `ShotAzimuth` under compass-referenced wind, and `KernelError::AxisAbsent` for the three wind
//! axes when the request being written TO uses segmented wind (`access.rs`'s module doc).
//! [`read_axis`] separately returns `None` for the same three wind axes when the request being
//! read FROM is segmented, and also for a handful of ordinarily-optional scalar fields
//! (`length_m`, `latitude_rad`, the two zero-POI offsets, the lateral sight offset, and -- true
//! of the code though not named in `read_axis`'s own doc comment -- `zero_distance_m` on a
//! request zeroed purely by an explicit angle) that just were not supplied on a request.
//!
//! An earlier revision of this module decided this PER LEG: forward excluded an axis only if
//! `A` itself refused it (or lacked a value to read from `B`), and backward only if `B` did.
//! That is wrong (review I1): if only ONE leg excludes an axis, the two legs stop measuring the
//! same counterfactual. A QNH-referenced `A` compared against a non-QNH `B`, for instance, has
//! forward keep `A`'s own altitude (refused) while backward hands `B` the whole of `A`'s
//! altitude (not refused) -- so `contribution(Atmosphere)` ends up as the full
//! temperature/pressure/humidity effect PLUS HALF of a completely separate altitude difference,
//! while the report's own `assumptions` claim that effect went to the interaction remainder.
//!
//! `plan_exclusions` fixes this by deciding exclusions for a WHOLE group ONCE, from the two
//! ORIGINAL requests together, before either swap direction runs: an axis is excluded from BOTH
//! legs if EITHER `a` or `b` would refuse it as a destination, or if it is present on exactly
//! one of `a`/`b` (review I2 -- the same splitting hazard, via `read_axis` instead of
//! `with_axis`: an ordinarily-optional field supplied on one saved profile and omitted on the
//! other, such as `latitude_rad`, would otherwise have forward keep `A`'s own value while
//! backward overwrites `B` with it, again half-attributing a real difference with no record of
//! it at all). Both cases are recorded in [`SolutionDiffReportV1::skipped_axes`] -- ALWAYS as a
//! pair, one entry per direction, so the report always explains both legs' identical treatment
//! of that axis -- with a `reason` reworded for this whole-group context rather than reused
//! verbatim from `with_axis`'s own per-axis-perturbation wording (review I1's second point: a
//! reason like "perturbing altitude would move air density without moving pressure" is
//! misleading here, where `Pressure` is a SEPARATE axis in the same group and IS copied
//! normally). An axis absent on BOTH `a` and `b` is not reported at all: there is nothing to
//! attribute either way, the same as if neither request had ever mentioned it.
//!
//! Either way, an exclusion never aborts the comparison -- only that one axis is left out of
//! that one group, on both legs. A GENUINE failure (most notably a `solve_v1` re-resolve failing
//! partway through applying a group -- for instance the pre-existing magnus/enhanced-spin-drift
//! conflict noted in `taxonomy.rs`'s Known Limitation (d), which a group swap can walk straight
//! into by turning one flag on while the request still carries the other from before) is NOT one
//! of these cases and propagates as an error, uncaught, exactly as `central_difference` and
//! `bisect_axis` already do for their own non-refusal failures.
//!
//! # Why exclusions are decided from the original requests, not the accumulated one
//!
//! `plan_exclusions` probes [`with_axis`] against `a` and `b` exactly as passed into
//! [`explain_difference`] -- never a request `swap_group` has already partly rewritten. This
//! matters because a group with more than one axis applies its writes one at a time,
//! re-resolving between them (`swap_group`, via a real `solve_v1` call) so each subsequent
//! [`with_axis`] call sees the previous write. That re-resolve goes through the reverse
//! conversion (`request_roundtrip.rs`), which ALWAYS clears `pressure_reference` and
//! `wind_reference` back to the omitted-field default (that module's own doc: the resolved
//! values already have the transform baked in, so echoing the mode back would apply it a second
//! time). That is correct for solving, but it means the very fact `with_axis`'s guards key on --
//! "was this request originally QNH- or compass-referenced" -- would be erased by the FIRST
//! re-resolve within a group, not preserved across it, if it were checked against that
//! progressively-modified request instead.
//!
//! `ShotGeometry` lists `TargetDistance`, `ShootingAngle` and `Cant` before `ShotAzimuth`: three
//! unrelated axes, each triggering a re-resolve, would be applied before it. Checking the
//! compass guard against a progressively-modified request would have it silently pass by the
//! time `ShotAzimuth` is reached, building exactly the physically-inverted counterfactual the
//! guard exists to prevent -- with no error and no skip recorded, because from `with_axis`'s
//! point of view at that moment the request no longer looks compass-referenced at all. Deciding
//! every axis's exclusion up front, from `a`/`b` as originally given, makes this impossible by
//! construction: there is no accumulated, partly-laundered request for the check to see in the
//! first place. `swap_group` itself never re-derives a refusal at all -- it trusts the
//! `excluded` list `plan_exclusions` already computed, safely, because a refusal can only become
//! MORE permissive as a group's later axes are applied (round-tripping only ever clears a
//! reference-mode echo, it never introduces one): an axis `plan_exclusions` already cleared for
//! both `a` and `b` is guaranteed to still be writable on `swap_group`'s running request,
//! however many earlier axes in the same group have already re-resolved it. The test
//! `shot_azimuth_is_refused_even_when_other_shot_geometry_axes_are_applied_first` pins this down
//! by exercising exactly that four-axis ordering.

use serde::{Deserialize, Serialize};

use crate::perturbation::access::{read_axis, with_axis, KernelError};
use crate::perturbation::taxonomy::{axes_in_group, InputAxis, InputGroup};
use crate::perturbation::{evaluate, kernel_solve_error, Observation};
use crate::solve_json::{PressureReferenceV1, ResolvedSolveRequestV1, SolveRequestV1, WindReferenceV1};

/// The difference between two [`Observation`]s at one range, `b - a`. SI throughout, same sign
/// convention as [`Observation`] itself.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct DeltaV1 {
    pub drop_m: f64,
    pub windage_m: f64,
    pub time_s: f64,
    pub velocity_mps: f64,
}

impl DeltaV1 {
    fn between(a: &Observation, b: &Observation) -> Self {
        DeltaV1 {
            drop_m: b.drop_m - a.drop_m,
            windage_m: b.windage_m - a.windage_m,
            time_s: b.time_s - a.time_s,
            velocity_mps: b.velocity_mps - a.velocity_mps,
        }
    }
    fn mean(x: Self, y: Self) -> Self {
        DeltaV1 {
            drop_m: 0.5 * (x.drop_m + y.drop_m),
            windage_m: 0.5 * (x.windage_m + y.windage_m),
            time_s: 0.5 * (x.time_s + y.time_s),
            velocity_mps: 0.5 * (x.velocity_mps + y.velocity_mps),
        }
    }
    fn neg(self) -> Self {
        DeltaV1 {
            drop_m: -self.drop_m,
            windage_m: -self.windage_m,
            time_s: -self.time_s,
            velocity_mps: -self.velocity_mps,
        }
    }
    fn add(self, o: Self) -> Self {
        DeltaV1 {
            drop_m: self.drop_m + o.drop_m,
            windage_m: self.windage_m + o.windage_m,
            time_s: self.time_s + o.time_s,
            velocity_mps: self.velocity_mps + o.velocity_mps,
        }
    }
    fn sub(self, o: Self) -> Self {
        DeltaV1 {
            drop_m: self.drop_m - o.drop_m,
            windage_m: self.windage_m - o.windage_m,
            time_s: self.time_s - o.time_s,
            velocity_mps: self.velocity_mps - o.velocity_mps,
        }
    }
}

/// Which leg of a group's symmetric counterfactual swap a [`SkippedAxisV1`] occurred on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SwapDirectionV1 {
    /// Building A with the group's axes replaced by B's values (measured against A).
    Forward,
    /// Building B with the group's axes replaced by A's values (measured against B, negated).
    Backward,
}

/// One taxonomy axis that could not be carried across during a group's counterfactual swap, and
/// why -- see the module doc's "Symmetric exclusion" section. Always recorded in pairs, one per
/// [`SwapDirectionV1`], even when only one direction independently hit the refusal or absence:
/// the OTHER direction is excluded too, to keep both legs measuring the same counterfactual.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SkippedAxisV1 {
    pub group: InputGroup,
    pub axis: InputAxis,
    pub direction: SwapDirectionV1,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupContributionV1 {
    pub group: InputGroup,
    pub delta: DeltaV1,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolutionDiffRowV1 {
    pub range_m: f64,
    pub total: DeltaV1,
    pub contributions: Vec<GroupContributionV1>,
    pub interaction_remainder: DeltaV1,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolutionDiffReportV1 {
    pub schema_version: u32,
    pub method: String,
    pub assumptions: Vec<String>,
    /// Axes excluded from a group's swap because one or both requests could not carry them --
    /// see the module doc's "Symmetric exclusion". Always excluded from BOTH swap directions,
    /// and reported as a pair (one `SkippedAxisV1` per direction), even when only one direction
    /// independently hit the refusal or absence. Never distributed or hidden: whatever effect a
    /// skipped axis would have had is left inside `SolutionDiffRowV1::interaction_remainder`,
    /// not silently absorbed into that axis's own group.
    pub skipped_axes: Vec<SkippedAxisV1>,
    pub rows: Vec<SolutionDiffRowV1>,
}

/// Word an axis's exclusion reason for the whole-group context, rather than reusing
/// `with_axis`'s own per-axis-perturbation wording verbatim (review I1): `with_axis`'s own
/// `Altitude` reason talks about perturbing altitude "without moving" pressure, which is
/// misleading here, where `Pressure` is a SEPARATE axis in the same group that DOES get copied
/// normally. `own_refusal` is this DIRECTION's own refusal (if any); `other_refusal` is the
/// OTHER direction's, used to explain a "sympathetic" exclusion -- one leg that would not have
/// refused on its own, excluded anyway to keep both legs measuring the same counterfactual.
fn describe_refusal(axis: InputAxis, own_refusal: &Option<String>, other_refusal: &Option<String>) -> String {
    let context = match axis {
        InputAxis::Altitude => {
            "altitude is tied to pressure under a QNH-referenced atmosphere; this group swap \
             cannot re-derive that relationship for a different altitude. \
             Temperature/RelativeHumidity/Latitude in the same group ARE copied normally; \
             Pressure is ALSO excluded whenever the two requests' altitudes actually differ \
             (see Pressure's own skipped_axes entry when that applies), but is copied normally \
             alongside Altitude's exclusion when the altitudes happen to match"
        }
        InputAxis::ShotAzimuth => {
            "the shot azimuth is tied to an earth-fixed wind bearing under compass-referenced \
             wind; this group swap cannot re-derive that bearing for a different azimuth, even \
             though TargetDistance/ShootingAngle/Cant/AimAzimuth/TargetHeight in the same group \
             ARE copied normally"
        }
        _ => "this axis cannot be swapped for one of the two requests being compared",
    };
    match own_refusal {
        Some(reason) => format!("{context} ({reason})"),
        None => format!(
            "excluded to match the other swap direction, which cannot swap this axis: {context} \
             ({})",
            other_refusal.as_deref().unwrap_or("no further detail available")
        ),
    }
}

/// Word a wind axis's exclusion reason (taxonomy.rs Known Limitation (c)): `own_segmented` is
/// whether THIS direction's destination is segmented (blocking the write); `other_segmented` is
/// whether the OTHER direction's destination is segmented (which, as this direction's source,
/// blocks reading a value to copy).
fn describe_wind_absence(own_segmented: bool, other_segmented: bool) -> String {
    match (own_segmented, other_segmented) {
        (true, true) => "both requests being compared use segmented wind; neither has a \
                         single scalar value for this axis"
            .to_string(),
        (true, false) => "this request's wind is segmented, so there is no single scalar \
                          field to write this axis into"
            .to_string(),
        (false, true) => "the other request's wind is segmented, so there is no single \
                          scalar value to read for this axis"
            .to_string(),
        (false, false) => {
            unreachable!("describe_wind_absence called when neither side is segmented")
        }
    }
}

/// True when `r`'s wind was entered compass-referenced. Shares `access.rs`'s own
/// `wind_reference_of` (review round 4 -- widened to `pub(crate)` for exactly this reuse)
/// rather than re-matching the two wind shapes a second time -- the same fact [`with_axis`]
/// checks to refuse `ShotAzimuth` there, because `WindDirection`'s OWN resolved value depends on
/// it too (module doc, "The derived-value exclusion rule", Instance 2).
fn is_compass_referenced(r: &ResolvedSolveRequestV1) -> bool {
    crate::perturbation::access::wind_reference_of(&r.wind) == Some(WindReferenceV1::Compass)
}

/// True when `r`'s atmosphere was entered as a QNH altimeter setting. Read directly from the
/// public `ResolvedAtmosphereV1::pressure_reference` echo -- the same fact [`with_axis`] checks
/// to refuse `Altitude` (`access.rs`) -- because `Pressure`'s OWN resolved value depends on it
/// too (module doc, "The derived-value exclusion rule", Instance 3).
fn is_qnh_referenced(r: &ResolvedSolveRequestV1) -> bool {
    r.atmosphere.pressure_reference == Some(PressureReferenceV1::Qnh)
}

/// Decide, for one group, which axes must be excluded from BOTH legs of the symmetric swap
/// between `a` and `b`, and produce the [`SkippedAxisV1`] entries explaining why -- see the
/// module doc's "Symmetric exclusion" section for why a per-leg decision is not good enough,
/// and "The derived-value exclusion rule" for the two checks that are not about refusal or
/// absence at all.
///
/// Five structurally different reasons exclude an axis, each checked directly against `a` and
/// `b` as originally given (see the module doc's "Why exclusions are decided from the original
/// requests"):
///
/// - One of the three wind axes, where either request's wind is segmented (taxonomy.rs Known
///   Limitation (c)) -- detected via [`read_axis`] returning `None`.
/// - `WindDirection`, derived from `ShotGeometry`'s `ShotAzimuth` whenever either request is
///   compass-referenced (review round 3, N1) -- the derived-value rule's case (a).
/// - `Pressure`, derived from this SAME group's `Altitude` whenever either request is
///   QNH-referenced, and `Altitude` is itself excluded in that case (review round 3, N2) -- the
///   derived-value rule's case (b).
/// - An ordinarily-optional axis present on exactly one of `a`/`b` (review I2).
/// - [`with_axis`] refusing the axis as a destination for `a` XOR `b` (review I1) -- a
///   QNH-referenced `Altitude` or a compass-referenced `ShotAzimuth`. Probed with each
///   request's OWN current value for that axis (guaranteed well-typed and finite), since
///   `with_axis`'s guards depend only on `axis` and the request's own structural fields, never
///   on the value being written.
///
/// An axis absent on BOTH `a` and `b` triggers none of these and is silently left out of both
/// the exclusion list and the report: there is nothing to attribute either way.
///
/// # Errors
///
/// Only an UNEXPECTED `with_axis` failure (not `AxisUnsupportedForRequest`) propagates -- for
/// instance a `TypeMismatch`/`NonFinite`, which should not be reachable here at all, since every
/// probed value came from `read_axis` for the identical axis on the identical request.
fn plan_exclusions(
    a: &ResolvedSolveRequestV1,
    b: &ResolvedSolveRequestV1,
    group: InputGroup,
) -> Result<(Vec<InputAxis>, Vec<SkippedAxisV1>), KernelError> {
    let mut excluded = Vec::new();
    let mut skipped = Vec::new();

    for &axis in axes_in_group(group) {
        let a_value = read_axis(a, axis);
        let b_value = read_axis(b, axis);

        if matches!(
            axis,
            InputAxis::WindSpeed | InputAxis::WindDirection | InputAxis::WindVertical
        ) {
            if a_value.is_none() || b_value.is_none() {
                excluded.push(axis);
                skipped.push(SkippedAxisV1 {
                    group,
                    axis,
                    direction: SwapDirectionV1::Forward,
                    reason: describe_wind_absence(a_value.is_none(), b_value.is_none()),
                });
                skipped.push(SkippedAxisV1 {
                    group,
                    axis,
                    direction: SwapDirectionV1::Backward,
                    reason: describe_wind_absence(b_value.is_none(), a_value.is_none()),
                });
                continue;
            }
            // N1 (review round 3, narrowed by F1): WindDirection's resolved value is ALSO
            // derived from ShotAzimuth -- ShotGeometry, a DIFFERENT group -- whenever either
            // side is compass-referenced (module doc, derived-value rule, Instance 2). But the
            // derivation only actually CONTAMINATES the swap when the azimuths differ: if
            // a.shot.shot_azimuth_rad == b.shot.shot_azimuth_rad, the compass-to-shooter-relative
            // shift is the SAME on both sides, so the difference between the two resolved
            // directions is exactly the difference between the two raw bearings -- a genuine
            // Wind-group effect, not azimuth contamination (F1: excluding unconditionally on the
            // mode flag over-fires, under-attributing Wind on the most natural same-shot-
            // direction compass comparison). `with_axis`'s OWN Altitude/ShotAzimuth guards do
            // not need this refinement -- they only ever see ONE request, so they cannot compare
            // -- but `plan_exclusions` has both and must.
            if axis == InputAxis::WindDirection
                && (is_compass_referenced(a) || is_compass_referenced(b))
                && a.shot.shot_azimuth_rad != b.shot.shot_azimuth_rad
            {
                excluded.push(axis);
                let reason = "the resolved wind direction is derived from shot_azimuth_rad (a \
                              ShotGeometry axis) under compass-referenced wind; swapping it \
                              here would attribute part of a shot-azimuth difference to Wind \
                              instead of leaving it in the interaction remainder"
                    .to_string();
                skipped.push(SkippedAxisV1 {
                    group,
                    axis,
                    direction: SwapDirectionV1::Forward,
                    reason: reason.clone(),
                });
                skipped.push(SkippedAxisV1 {
                    group,
                    axis,
                    direction: SwapDirectionV1::Backward,
                    reason,
                });
            }
            continue;
        }

        // N2 (review round 3, narrowed by F1): Pressure's resolved value is derived from
        // Altitude -- the SAME group, but Altitude is itself excluded whenever either side is
        // QNH-referenced, so the pairing that would ordinarily make same-group derivation
        // harmless is broken (module doc, derived-value rule, Instance 3). But this only
        // actually breaks the pairing when the altitudes differ: if
        // a.atmosphere.altitude_m == b.atmosphere.altitude_m, the QNH-reduced pressures differ
        // only by the raw altimeter setting, and the swapped value stays physically valid at
        // the (unmoved) destination altitude -- excluding unconditionally on the mode flag
        // over-fires, under-attributing Atmosphere on the most natural same-location QNH
        // comparison (F1). See the WindDirection guard above for why `with_axis`'s own
        // unconditional guards do not need (and cannot use) this same refinement.
        if axis == InputAxis::Pressure
            && (is_qnh_referenced(a) || is_qnh_referenced(b))
            && a.atmosphere.altitude_m != b.atmosphere.altitude_m
        {
            excluded.push(axis);
            let reason = "the resolved station pressure is derived from altitude_m under a \
                          QNH-referenced atmosphere, and Altitude is excluded from this same \
                          comparison; swapping Pressure alone would move a value that is only \
                          physically valid at the OTHER request's (unswapped) altitude"
                .to_string();
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Forward,
                reason: reason.clone(),
            });
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Backward,
                reason,
            });
            continue;
        }

        if a_value.is_some() != b_value.is_some() {
            excluded.push(axis);
            let reason = "this axis is present on only one of the two requests being compared \
                          (absent on the other), so there is no symmetric value to swap in \
                          either direction"
                .to_string();
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Forward,
                reason: reason.clone(),
            });
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Backward,
                reason,
            });
            continue;
        }
        let (Some(a_value), Some(b_value)) = (a_value, b_value) else {
            continue; // absent on BOTH sides -- nothing to attribute, not reported.
        };

        let a_refusal = match with_axis(a, axis, a_value) {
            Ok(_) => None,
            Err(KernelError::AxisUnsupportedForRequest { reason, .. }) => Some(reason.to_string()),
            Err(other) => return Err(other),
        };
        let b_refusal = match with_axis(b, axis, b_value) {
            Ok(_) => None,
            Err(KernelError::AxisUnsupportedForRequest { reason, .. }) => Some(reason.to_string()),
            Err(other) => return Err(other),
        };
        if a_refusal.is_some() || b_refusal.is_some() {
            excluded.push(axis);
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Forward,
                reason: describe_refusal(axis, &a_refusal, &b_refusal),
            });
            skipped.push(SkippedAxisV1 {
                group,
                axis,
                direction: SwapDirectionV1::Backward,
                reason: describe_refusal(axis, &b_refusal, &a_refusal),
            });
        }
    }
    Ok((excluded, skipped))
}

/// Build the [`SolveRequestV1`] representing `dst` with every NON-EXCLUDED axis of `group`
/// copied from `src`, applying each axis one at a time and re-resolving between writes (via a
/// real `solve_v1` call) so a later axis in the same group sees the effect of an earlier one --
/// most importantly, so `ZeroSightGeometry`'s `SightHeight` re-derives the elevation for the
/// destination's own physics after `MuzzleAngle` (listed first specifically for this) has
/// written the source's (see the module doc's "`MuzzleAngle` is derived" section). Does not
/// itself evaluate a trajectory -- the caller does that once, over the whole `ranges_m` slice,
/// exactly as [`central_difference`](crate::perturbation::central_difference) already does for
/// its own pair of counterfactual requests.
///
/// `excluded` is computed ONCE per group by `plan_exclusions`, from `a`/`b` exactly as passed
/// into [`explain_difference`], before either swap direction runs -- never re-derived here. See
/// the module doc's "Why exclusions are decided from the original requests" for why that
/// ordering is required for correctness, not just convenience: every axis reaching this loop is
/// guaranteed still writable on the running request, however many earlier axes in this same
/// group have already re-resolved it.
///
/// # Errors
///
/// Any [`with_axis`] or `solve_v1` failure propagates via `?` unchanged: by the time an axis
/// reaches this loop it is not supposed to be refusable, so a failure here is either a genuine
/// bug (a `read_axis`/`with_axis` mismatch) or an unrelated solve failure (such as the
/// magnus/enhanced-spin-drift conflict noted in `taxonomy.rs`'s Known Limitation (d)) -- never a
/// structurally unrepresentable axis, which `plan_exclusions` has already filtered out.
fn swap_group(
    dst: &ResolvedSolveRequestV1,
    src: &ResolvedSolveRequestV1,
    group: InputGroup,
    excluded: &[InputAxis],
) -> Result<SolveRequestV1, KernelError> {
    let mut current = dst.clone();
    for &axis in axes_in_group(group) {
        if excluded.contains(&axis) {
            continue;
        }
        let Some(v) = read_axis(src, axis) else {
            // plan_exclusions already excludes any axis present on exactly one side; a `None`
            // reaching here means it is absent on BOTH sides -- a plain no-op, nothing to copy.
            continue;
        };
        let req = with_axis(&current, axis, v)?;
        current = crate::solve_v1::solve_v1(req)
            .map_err(kernel_solve_error)?
            .resolved_request;
    }
    Ok((&current).into())
}

/// Attribute the difference between two fully resolved solve results to the seven input groups
/// (MBA-1345).
///
/// For each [`InputGroup`] `g`, the forward leg replaces `g`'s axes on `a` with `b`'s values and
/// measures the change against `a`; the backward leg replaces `g`'s axes on `b` with `a`'s
/// values and measures the change against `b`, negated. The group's reported contribution is
/// the MEAN of the two, which is what makes the decomposition independent of which request is
/// treated as the "before" and which as the "after". Whatever the seven groups do not explain
/// -- genuine nonlinear interaction between them -- is reported once per range, as
/// `SolutionDiffRowV1::interaction_remainder`, and is never distributed across the groups: for
/// correlated inputs there is no unique causal attribution, and pretending otherwise is the
/// failure this design exists to avoid.
///
/// `ranges_m` is passed straight through to two calls to [`evaluate`] per group (one per swap
/// direction), plus two more for `a`/`b` themselves -- 16 solves for the fixed seven-group
/// taxonomy, each covering every range in `ranges_m` at once rather than once per range (see
/// [`evaluate`]'s own cost note). That 16 is NOT the total, though: building each group's
/// counterfactual (`swap_group`) re-resolves via a full `solve_v1` call after EVERY applied
/// axis, regardless of whether that specific axis is `requires_rezero` -- not just the ones that
/// invalidate a zero -- because [`with_axis`] always needs a freshly resolved request to apply
/// its NEXT write onto. For a request with N present axes (up to the full 32; 5-6 are commonly
/// absent, e.g. `length_m`/`latitude_rad`/the POI offsets, so N is often closer to 27), that is
/// up to N solves per direction, 2N total, ON TOP of the 16 -- so a typical call costs on the
/// order of 16 + 2*27 = 70 solves, not 16. `plan_exclusions` itself adds no solves (it only
/// calls [`read_axis`]/[`with_axis`], never `solve_v1`), so excluding an axis only ever reduces
/// this total, never increases it.
///
/// # Axes that cannot be swapped
///
/// See the module doc's "Symmetric exclusion" section for the full explanation. In short:
/// [`with_axis`] can refuse an axis for a structural reason specific to `a` or `b` (a
/// QNH-referenced altitude, a compass-referenced shot azimuth, one of the three wind axes under
/// segmented wind), and an ordinarily-optional axis can be present on only one of `a`/`b`. An
/// axis affected by either is recorded in the returned report's `skipped_axes` -- one entry per
/// swap direction, always both together, even though only one direction may have hit the
/// refusal or absence directly -- and excluded from BOTH legs of that group's swap, so the two
/// legs keep measuring the same counterfactual; the rest of that axis's group is still
/// attributed normally, and a refusal never aborts the comparison. An axis absent from BOTH `a`
/// and `b` (no `length_m` supplied on either, for instance) is skipped without being recorded:
/// there is nothing there to attribute either way.
///
/// # Errors
///
/// Returns whatever [`evaluate`], [`with_axis`], or `plan_exclusions` error on `a`, `b`, or
/// either swapped counterfactual, MINUS the refusal/absence cases above
/// (`KernelError::AxisUnsupportedForRequest`, and the presence-asymmetry case, which never even
/// reaches a `KernelError`), which are caught and recorded rather than returned. Every other
/// error -- most notably a `solve_v1` re-resolve failing validation partway through building a
/// group's counterfactual, such as the pre-existing magnus/enhanced-spin-drift conflict noted in
/// `taxonomy.rs`'s Known Limitation (d) -- propagates unchanged: it is a genuine failure, not a
/// structurally unrepresentable axis, and must not be silently absorbed into a skip or a zero
/// contribution.
///
/// `a` and `b` are read-only throughout: every counterfactual request is built from a clone
/// (`swap_group`), never a mutation of either input.
pub fn explain_difference(
    a: &ResolvedSolveRequestV1,
    b: &ResolvedSolveRequestV1,
    ranges_m: &[f64],
) -> Result<SolutionDiffReportV1, KernelError> {
    let obs_a = evaluate(&a.into(), ranges_m)?;
    let obs_b = evaluate(&b.into(), ranges_m)?;

    let mut skipped_axes = Vec::new();
    let mut per_group: Vec<(InputGroup, Vec<DeltaV1>)> = Vec::with_capacity(InputGroup::ALL.len());
    for &group in InputGroup::ALL {
        let (excluded, mut group_skips) = plan_exclusions(a, b, group)?;
        skipped_axes.append(&mut group_skips);

        // Each swapped request is evaluated exactly ONCE over the whole `ranges_m` slice below
        // -- never once per range -- so N ranges do not multiply the cost of a
        // `requires_rezero` axis's zero search inside `swap_group`.
        let fwd_req = swap_group(a, b, group, &excluded)?;
        let bwd_req = swap_group(b, a, group, &excluded)?;
        let fwd_obs = evaluate(&fwd_req, ranges_m)?;
        let bwd_obs = evaluate(&bwd_req, ranges_m)?;

        let mut deltas = Vec::with_capacity(ranges_m.len());
        for i in 0..ranges_m.len() {
            let forward = DeltaV1::between(&obs_a[i], &fwd_obs[i]); // A with g<-B, vs A
            let backward = DeltaV1::between(&obs_b[i], &bwd_obs[i]).neg(); // B with g<-A, vs B, negated
            deltas.push(DeltaV1::mean(forward, backward));
        }
        per_group.push((group, deltas));
    }

    let mut rows = Vec::with_capacity(ranges_m.len());
    for (i, &range_m) in ranges_m.iter().enumerate() {
        let total = DeltaV1::between(&obs_a[i], &obs_b[i]);
        let contributions: Vec<GroupContributionV1> = per_group
            .iter()
            .map(|(g, d)| GroupContributionV1 {
                group: *g,
                delta: d[i],
            })
            .collect();
        let summed = contributions
            .iter()
            .fold(DeltaV1::default(), |acc, c| acc.add(c.delta));
        rows.push(SolutionDiffRowV1 {
            range_m,
            total,
            contributions,
            interaction_remainder: total.sub(summed),
        });
    }

    Ok(SolutionDiffReportV1 {
        schema_version: 1,
        method: "symmetric_group_counterfactual".to_string(),
        assumptions: vec![
            "Group contributions are symmetric counterfactuals: the mean of swapping the group \
             in each direction, so the result does not depend on replacement order."
                .to_string(),
            "Nonlinear interaction between groups is reported as an explicit interaction \
             remainder and is NOT distributed across groups. For correlated inputs no unique \
             causal attribution exists."
                .to_string(),
            "An axis that could not be swapped for one or both requests (see skipped_axes) is \
             excluded from BOTH directions of that axis's group, so the two directions keep \
             measuring the same counterfactual; any real effect the axis would have had is \
             folded into the interaction remainder, not silently dropped and not partially \
             attributed to its group."
                .to_string(),
        ],
        skipped_axes,
        rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn resolved(mv: f64, temp: f64) -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": mv, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"temperature_k": temp}, "wind": {"speed_mps": 3.0,
                           "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    /// Varies ONLY `ballistic_coefficient` (ProjectileDrag) and wind speed (Wind), holding
    /// muzzle velocity, temperature, and every sight/zero knob identical to `resolved()`'s own
    /// baseline -- see `every_groups_contribution_matches_an_independent_recomputation_bc_and_wind_fixture`.
    fn resolved_with_bc_and_wind_speed(
        ballistic_coefficient: f64,
        wind_speed_mps: f64,
    ) -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": ballistic_coefficient},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"temperature_k": 288.0},
            "wind": {"speed_mps": wind_speed_mps,
                     "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    fn resolved_with_effects(
        mv: f64,
        magnus: bool,
        enhanced_spin_drift: bool,
    ) -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": mv, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {"temperature_k": 288.0},
            "wind": {},
            "solver": {},
            "effects": {"magnus": magnus, "enhanced_spin_drift": enhanced_spin_drift},
            "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    /// Acceptance criterion: identical requests produce zero deltas and zero contributions, and
    /// nothing is reported as skipped (both requests are ordinary constant-wind, absolute-
    /// pressure, shooter-relative requests, so no refusal or structural absence should fire).
    #[test]
    fn identical_requests_produce_zero_everything() {
        let a = resolved(823.0, 288.0);
        let rep = explain_difference(&a, &a, &[300.0, 600.0]).unwrap();
        assert!(
            rep.skipped_axes.is_empty(),
            "expected no skipped axes for two ordinary, identical requests, got {:?}",
            rep.skipped_axes
        );
        for row in &rep.rows {
            assert!(row.total.drop_m.abs() < 1e-9, "total drop {}", row.total.drop_m);
            assert!(row.interaction_remainder.drop_m.abs() < 1e-9);
            for c in &row.contributions {
                assert!(
                    c.delta.drop_m.abs() < 1e-9,
                    "{:?} contributed {}",
                    c.group,
                    c.delta.drop_m
                );
            }
        }
    }

    /// Swapping the two requests negates every reported quantity.
    #[test]
    fn the_decomposition_is_antisymmetric() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let ab = explain_difference(&a, &b, &[600.0]).unwrap();
        let ba = explain_difference(&b, &a, &[600.0]).unwrap();
        assert!((ab.rows[0].total.drop_m + ba.rows[0].total.drop_m).abs() < 1e-6);
        for (x, y) in ab.rows[0]
            .contributions
            .iter()
            .zip(ba.rows[0].contributions.iter())
        {
            assert_eq!(x.group, y.group);
            assert!(
                (x.delta.drop_m + y.delta.drop_m).abs() < 1e-6,
                "{:?} not antisymmetric",
                x.group
            );
        }
    }

    /// Contributions plus the remainder must reconstruct the total exactly.
    #[test]
    fn contributions_plus_remainder_equal_the_total() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let rep = explain_difference(&a, &b, &[600.0]).unwrap();
        let row = &rep.rows[0];
        let sum: f64 = row.contributions.iter().map(|c| c.delta.drop_m).sum();
        assert!((sum + row.interaction_remainder.drop_m - row.total.drop_m).abs() < 1e-9);
    }

    #[test]
    fn the_report_states_its_method_and_assumptions() {
        let a = resolved(823.0, 288.0);
        let rep = explain_difference(&a, &a, &[300.0]).unwrap();
        assert_eq!(rep.method, "symmetric_group_counterfactual");
        assert!(rep.assumptions.iter().any(|s| s.contains("interaction")));
    }

    /// The two source requests must never be mutated by the comparison -- `explain_difference`
    /// takes `&ResolvedSolveRequestV1`, and every counterfactual is built from a clone
    /// (`swap_group`), never a mutation of `a`/`b` themselves.
    #[test]
    fn source_requests_are_not_mutated_by_the_comparison() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let a_before = a.clone();
        let b_before = b.clone();
        let _ = explain_difference(&a, &b, &[300.0, 600.0]).unwrap();
        assert_eq!(a, a_before, "the first source request changed");
        assert_eq!(b, b_before, "the second source request changed");
    }

    /// A forward-only (or backward-only) implementation would report either endpoint instead of
    /// their mean. Recompute the MuzzleVelocity group's forward and backward legs directly,
    /// using the same building blocks `explain_difference` uses internally, and confirm (a) they
    /// genuinely differ here -- or averaging them would be indistinguishable from reporting
    /// either alone -- and (b) the reported contribution is exactly their mean, not either one.
    ///
    /// The sanity gate below requires `|forward - backward| > 2e-6`, not `1e-6`: since
    /// `mean - forward == (backward - forward) / 2` (and symmetrically for `mean - backward`),
    /// a gate of exactly `1e-6` on `|forward - backward|` would only guarantee
    /// `|mean - forward| > 5e-7`, half of the `1e-6` margin the two endpoint assertions below
    /// actually need -- a value that just cleared the old, looser gate could still trip a
    /// false failure on those. `2e-6` is the smallest gate that genuinely implies both.
    #[test]
    fn contribution_is_the_mean_of_a_genuinely_different_forward_and_backward() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let ranges = [600.0];

        let obs_a = evaluate(&(&a).into(), &ranges).unwrap();
        let obs_b = evaluate(&(&b).into(), &ranges).unwrap();
        let (excluded, _skipped) =
            plan_exclusions(&a, &b, InputGroup::MuzzleVelocity).unwrap();
        let fwd_req = swap_group(&a, &b, InputGroup::MuzzleVelocity, &excluded).unwrap();
        let bwd_req = swap_group(&b, &a, InputGroup::MuzzleVelocity, &excluded).unwrap();
        let fwd_obs = evaluate(&fwd_req, &ranges).unwrap();
        let bwd_obs = evaluate(&bwd_req, &ranges).unwrap();
        let forward = fwd_obs[0].drop_m - obs_a[0].drop_m;
        let backward = obs_b[0].drop_m - bwd_obs[0].drop_m;

        assert!(
            (forward - backward).abs() > 2e-6,
            "forward ({forward}) and backward ({backward}) must genuinely differ here, or this \
             test cannot distinguish 'the mean of the two' from 'either endpoint alone'"
        );

        let rep = explain_difference(&a, &b, &ranges).unwrap();
        let mv = rep.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::MuzzleVelocity)
            .unwrap();
        let expected_mean = 0.5 * (forward + backward);
        assert!(
            (mv.delta.drop_m - expected_mean).abs() < 1e-9,
            "reported {} expected mean {}",
            mv.delta.drop_m,
            expected_mean
        );
        assert!(
            (mv.delta.drop_m - forward).abs() > 1e-6,
            "reported value must not equal the forward leg alone"
        );
        assert!(
            (mv.delta.drop_m - backward).abs() > 1e-6,
            "reported value must not equal the backward leg alone"
        );
    }

    /// Cross-checks EVERY group's reported contribution against an independent recomputation
    /// (same primitives -- `plan_exclusions` + `swap_group` + `evaluate` -- but the mean is
    /// taken here, outside `explain_difference`'s own code path) and the remainder against the
    /// same total-minus-sum formula computed independently. This is the check that would fail
    /// under: a forward-only computation (the independent mean would disagree with a
    /// forward-only report), a remainder silently zeroed or redistributed into the groups (the
    /// independent remainder would disagree with a tampered report), or two groups'
    /// contributions swapped (checked BY GROUP IDENTITY). Returns the report so callers can
    /// layer their own fixture-specific "this group must be zero/non-negligible" assertions on
    /// top, since which groups should be non-negligible depends on which axes the caller's own
    /// fixture actually varies.
    fn assert_matches_independent_recomputation(
        a: &crate::solve_json::ResolvedSolveRequestV1,
        b: &crate::solve_json::ResolvedSolveRequestV1,
        range_m: f64,
    ) -> SolutionDiffReportV1 {
        let ranges = [range_m];
        let rep = explain_difference(a, b, &ranges).unwrap();
        let row = &rep.rows[0];

        let obs_a = evaluate(&a.into(), &ranges).unwrap();
        let obs_b = evaluate(&b.into(), &ranges).unwrap();

        let mut independent_sum = DeltaV1::default();
        for &group in InputGroup::ALL {
            let (excluded, _skipped) = plan_exclusions(a, b, group).unwrap();
            let fwd_req = swap_group(a, b, group, &excluded).unwrap();
            let bwd_req = swap_group(b, a, group, &excluded).unwrap();
            let fwd_obs = evaluate(&fwd_req, &ranges).unwrap();
            let bwd_obs = evaluate(&bwd_req, &ranges).unwrap();
            let forward = DeltaV1::between(&obs_a[0], &fwd_obs[0]);
            let backward = DeltaV1::between(&obs_b[0], &bwd_obs[0]).neg();
            let expected = DeltaV1::mean(forward, backward);

            let reported = row
                .contributions
                .iter()
                .find(|c| c.group == group)
                .unwrap_or_else(|| panic!("{group:?} missing from the report"));
            assert!(
                (reported.delta.drop_m - expected.drop_m).abs() < 1e-9,
                "{group:?}: reported drop {} independent {}",
                reported.delta.drop_m,
                expected.drop_m
            );
            assert!(
                (reported.delta.windage_m - expected.windage_m).abs() < 1e-9,
                "{group:?}: windage mismatch"
            );
            assert!(
                (reported.delta.time_s - expected.time_s).abs() < 1e-9,
                "{group:?}: time mismatch"
            );
            assert!(
                (reported.delta.velocity_mps - expected.velocity_mps).abs() < 1e-9,
                "{group:?}: velocity mismatch"
            );

            independent_sum = independent_sum.add(expected);
        }

        let total = DeltaV1::between(&obs_a[0], &obs_b[0]);
        let expected_remainder = total.sub(independent_sum);
        assert!(
            (row.interaction_remainder.drop_m - expected_remainder.drop_m).abs() < 1e-9,
            "reported remainder {} independent remainder {}",
            row.interaction_remainder.drop_m,
            expected_remainder.drop_m
        );

        rep
    }

    /// `resolved()` only ever varies `mv` and `temp`, so a group whose axes are entirely
    /// independent of both -- ProjectileDrag, Wind, ShotGeometry, Effects, AND (review C1)
    /// ZeroSightGeometry -- must come back EXACTLY zero: nothing to swap, since every one of its
    /// axis VALUES is identical on `a` and `b`. Before the C1 fix, `ZeroSightGeometry` was NOT
    /// zero here despite that: it also carried `MuzzleAngle`, the RESOLVED elevation, which is a
    /// function of muzzle velocity and air density -- so it silently imported the OTHER
    /// request's baked-in elevation even though every explicit sight/zero knob was identical.
    /// `axes_in_group(ZeroSightGeometry)` now lists `MuzzleAngle` first specifically so a
    /// following `requires_rezero` axis re-derives the elevation for the destination's own
    /// physics instead (see the module doc). MuzzleVelocity and Atmosphere differ directly and
    /// must be non-negligible and pairwise distinguishable, so a label swap involving either
    /// would be an unmissable mismatch against the per-group check inside
    /// `assert_matches_independent_recomputation`.
    #[test]
    fn every_groups_contribution_matches_an_independent_recomputation() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let rep = assert_matches_independent_recomputation(&a, &b, 600.0);
        let row = &rep.rows[0];

        for g in [
            InputGroup::ProjectileDrag,
            InputGroup::ZeroSightGeometry,
            InputGroup::Wind,
            InputGroup::ShotGeometry,
            InputGroup::Effects,
        ] {
            let c = row.contributions.iter().find(|x| x.group == g).unwrap();
            assert!(
                c.delta.drop_m.abs() < 1e-9,
                "{g:?} does not differ between a and b in this fixture at all, expected \
                 exactly 0, got {}",
                c.delta.drop_m
            );
        }
        let mv = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::MuzzleVelocity)
            .unwrap();
        let atmosphere = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        for (name, c) in [("MuzzleVelocity", mv), ("Atmosphere", atmosphere)] {
            assert!(
                c.delta.drop_m.abs() > 0.01,
                "{name} should have a substantial, non-negligible contribution here, got {}",
                c.delta.drop_m
            );
        }
        assert!(
            (mv.delta.drop_m - atmosphere.delta.drop_m).abs() > 0.01,
            "MuzzleVelocity ({}) and Atmosphere ({}) should be distinguishable, not just both \
             'non-negligible'",
            mv.delta.drop_m,
            atmosphere.delta.drop_m
        );

        // The remainder itself must also be genuinely nonzero here (real second-order
        // interaction between muzzle velocity and temperature/elevation) -- not identically
        // zero, and not silently redistributed into the groups above (which would make THIS
        // remainder read back as zero while the groups' sum still matched `total`). Measured at
        // ~9.2e-4 after the C1 fix (it was ~4.9e-2 before: MuzzleVelocity's own re-zero and
        // Atmosphere's own "not re-zeroed" treatment now correctly explain almost all of the
        // elevation response themselves, instead of a chunk of it being force-fed through a
        // ZeroSightGeometry that should have been zero) -- 1e-4 stays comfortably below that
        // measured value while remaining far above floating-point noise (~1e-9-1e-12).
        assert!(
            row.interaction_remainder.drop_m.abs() > 1e-4,
            "expected a non-negligible interaction remainder in this fixture, got {}",
            row.interaction_remainder.drop_m
        );
    }

    /// Review minor: `resolved()`'s fixture never varies `ballistic_coefficient` or wind speed,
    /// so a mislabeling among ProjectileDrag/Wind (or between either of them and any of the
    /// other always-zero groups) would be INVISIBLE to the test above -- every one of those
    /// groups reports zero there either way. A second fixture, varying only BC and wind speed
    /// (holding mv/temp/sight/zero geometry identical), exercises ProjectileDrag and Wind with a
    /// substantial, distinguishable value instead, closing that gap.
    #[test]
    fn every_groups_contribution_matches_an_independent_recomputation_bc_and_wind_fixture() {
        let a = resolved_with_bc_and_wind_speed(0.243, 3.0);
        let b = resolved_with_bc_and_wind_speed(0.300, 6.0);
        let rep = assert_matches_independent_recomputation(&a, &b, 600.0);
        let row = &rep.rows[0];

        for g in [
            InputGroup::MuzzleVelocity,
            InputGroup::ZeroSightGeometry,
            InputGroup::Atmosphere,
            InputGroup::ShotGeometry,
            InputGroup::Effects,
        ] {
            let c = row.contributions.iter().find(|x| x.group == g).unwrap();
            assert!(
                c.delta.drop_m.abs() < 1e-9,
                "{g:?} does not differ between a and b in this fixture at all, expected \
                 exactly 0, got {}",
                c.delta.drop_m
            );
        }
        let drag = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::ProjectileDrag)
            .unwrap();
        let wind = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Wind)
            .unwrap();
        // ProjectileDrag (a higher BC) is measured on drop, its dominant effect; Wind (a
        // stronger crosswind) is measured on windage, ITS dominant effect -- checking each
        // group in the quantity it actually moves, rather than requiring both to be large in
        // the SAME quantity, so this does not depend on the two effects happening to be
        // comparable in magnitude.
        assert!(
            drag.delta.drop_m.abs() > 0.01,
            "ProjectileDrag should have a substantial drop contribution here, got {}",
            drag.delta.drop_m
        );
        assert!(
            wind.delta.windage_m.abs() > 0.01,
            "Wind should have a substantial windage contribution here, got {}",
            wind.delta.windage_m
        );
        // And each group must not be silently carrying the OTHER's signature: ProjectileDrag's
        // own windage move and Wind's own drop move should both be comparatively small, which
        // is what a labeling swap between these two specific groups would violate.
        assert!(
            drag.delta.windage_m.abs() < wind.delta.windage_m.abs(),
            "ProjectileDrag's windage move ({}) should be smaller than Wind's own ({})",
            drag.delta.windage_m,
            wind.delta.windage_m
        );
        assert!(
            wind.delta.drop_m.abs() < drag.delta.drop_m.abs(),
            "Wind's drop move ({}) should be smaller than ProjectileDrag's own ({})",
            wind.delta.drop_m,
            drag.delta.drop_m
        );
    }

    /// Review round 3, N3: moving `MuzzleAngle` to the front of `ZeroSightGeometry` (C1) is
    /// only safe because a LATER `requires_rezero` axis's clear-gate fires and re-derives the
    /// angle -- and that gate is conditioned on `zero_distance_m` being present. For an
    /// angle-only request (`zero_distance_m == None`), the gate never fires, so `MuzzleAngle`
    /// genuinely IS the independent input and must still be swapped. Two angle-only requests
    /// with distinct `muzzle_angle_rad` values, both with no `zero_distance_m` at all, must
    /// still show a substantial `ZeroSightGeometry` contribution -- confirming the C1 reorder
    /// did not silently turn `MuzzleAngle` into a no-op for the one case it must still cover.
    #[test]
    fn muzzle_angle_is_still_swapped_for_an_angle_only_request() {
        let build = |muzzle_angle_rad: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0, "muzzle_angle_rad": muzzle_angle_rad},
                "atmosphere": {"temperature_k": 288.0},
                "wind": {},
                "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(0.010)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(0.030)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert_eq!(
            a.shot.zero_distance_m, None,
            "fixture assumption: angle-only, no zero distance at all"
        );
        assert_eq!(
            b.shot.zero_distance_m, None,
            "fixture assumption: angle-only, no zero distance at all"
        );

        // Uses the shared helper (same primitives as explain_difference, mean taken outside its
        // own code path) so this ALSO confirms the swapped value matches an independent
        // recomputation, not just that it is non-zero.
        let rep = assert_matches_independent_recomputation(&a, &b, 600.0);
        let row = &rep.rows[0];
        let zero_sight_geometry = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::ZeroSightGeometry)
            .unwrap();
        assert!(
            zero_sight_geometry.delta.drop_m.abs() > 0.01,
            "MuzzleAngle must still be swapped for an angle-only request -- expected a \
             substantial ZeroSightGeometry contribution, got {}",
            zero_sight_geometry.delta.drop_m
        );
    }

    /// Requirement (1): `with_axis` refuses `Altitude` under a QNH-referenced atmosphere. The
    /// refusal must be recorded, not silently dropped, and must not abort the comparison. Per
    /// review I1, exclusion must apply to BOTH swap directions (`b` is not itself QNH-referenced,
    /// so only the Forward leg would refuse Altitude on its own -- Backward must be excluded
    /// anyway, to keep the two legs measuring the same counterfactual).
    #[test]
    fn altitude_is_skipped_and_recorded_under_qnh_pressure() {
        let qnh_json = serde_json::json!({
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
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = resolved(870.0, 300.0);
        assert_eq!(
            b.atmosphere.pressure_reference, None,
            "fixture assumption: b is NOT QNH-referenced, so only a's QNH-ness drives this test"
        );

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("a refused axis must not abort the whole comparison");
        for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
            let hit = rep
                .skipped_axes
                .iter()
                .find(|s| {
                    s.group == InputGroup::Atmosphere
                        && s.axis == InputAxis::Altitude
                        && s.direction == direction
                })
                .unwrap_or_else(|| {
                    panic!(
                        "expected a {direction:?}-direction Altitude skip under QNH pressure \
                         (review I1: both directions must be excluded, not just the one that \
                         independently refused), got {:?}",
                        rep.skipped_axes
                    )
                });
            assert!(
                hit.reason.to_lowercase().contains("qnh"),
                "reason should name QNH: {}",
                hit.reason
            );
        }
    }

    /// Review I1's explicit ask: not just that Altitude is recorded as skipped, but that the
    /// resulting Atmosphere CONTRIBUTION is actually correct -- reflecting ONLY the other axes,
    /// with no partial leakage from the excluded Altitude difference. Proven by showing the
    /// reported contribution is UNCHANGED when b's altitude is moved from one value to ANOTHER:
    /// if Altitude's difference were still (even partially) leaking into the contribution,
    /// moving it would change the reported number; since Altitude is excluded from BOTH legs
    /// either way, it must not.
    ///
    /// Deliberately compares TWO altitudes that BOTH differ from `a`'s (1200 and 900), neither
    /// equal to `a`'s own 500 -- review round 4, F1 excludes `Pressure` only when the two
    /// requests' altitudes actually differ, so comparing against an altitude that happens to
    /// MATCH `a`'s would stop excluding `Pressure` in that one leg and confound this test with a
    /// second, unrelated effect (see `pressure_exclusion_does_not_leak_a_partial_effect_into_atmosphere`'s
    /// own doc comment for that failure mode in detail). `b`'s pressure here is a fixed,
    /// altitude-INDEPENDENT absolute 98000 Pa (not QNH), so unlike that other test, changing
    /// `b`'s altitude between 1200 and 900 changes nothing about `b`'s own actual physics either
    /// (`calculate_atmosphere` ignores `altitude_m` whenever temperature and pressure are both
    /// given directly, which they always are here) -- so this one tolerates NO residual at all,
    /// not even the small second-order kind the Pressure test documents.
    #[test]
    fn altitude_exclusion_does_not_leak_a_partial_effect_into_atmosphere() {
        let qnh_json = serde_json::json!({
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
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json).unwrap(),
        )
        .unwrap()
        .resolved_request;

        let b_json = |altitude_m: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0},
                "atmosphere": {"altitude_m": altitude_m, "temperature_k": 300.0,
                               "pressure_pa": 98000.0},
                "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
            })
            .to_string()
        };
        let b_altitude_1200 = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&b_json(1200.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b_altitude_900 = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&b_json(900.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        // Fixture assumption: NEITHER b variant's altitude may equal a's 500 -- see the doc
        // comment above for why that specific coincidence would confound this test.
        assert_ne!(b_altitude_1200.atmosphere.altitude_m, a.atmosphere.altitude_m);
        assert_ne!(b_altitude_900.atmosphere.altitude_m, a.atmosphere.altitude_m);

        let rep_1200 = explain_difference(&a, &b_altitude_1200, &[300.0]).unwrap();
        let rep_900 = explain_difference(&a, &b_altitude_900, &[300.0]).unwrap();

        // Sanity: Altitude must be excluded on BOTH directions in BOTH comparisons, so this is
        // a fair, apples-to-apples check of the same code path (not, say, one comparison
        // exercising the exclusion and the other happening to skip it because the values
        // coincided).
        for rep in [&rep_1200, &rep_900] {
            for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
                assert!(
                    rep.skipped_axes.iter().any(|s| s.group == InputGroup::Atmosphere
                        && s.axis == InputAxis::Altitude
                        && s.direction == direction),
                    "expected Altitude excluded on {direction:?} in both comparisons"
                );
            }
        }

        let atmosphere_1200 = rep_1200.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        let atmosphere_900 = rep_900.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        assert!(
            (atmosphere_1200.delta.drop_m - atmosphere_900.delta.drop_m).abs() < 1e-9,
            "Atmosphere's contribution changed when b's altitude changed (1200 -> 900 m) even \
             though Altitude is excluded from both comparisons -- {} vs {} -- a partial \
             altitude effect must be leaking through",
            atmosphere_1200.delta.drop_m,
            atmosphere_900.delta.drop_m
        );
        // And it must be genuinely non-zero (temperature alone is a real effect here), or the
        // equality above would be trivially true for the wrong reason (both sides zero).
        assert!(
            atmosphere_900.delta.drop_m.abs() > 0.001,
            "expected a non-negligible Atmosphere contribution from temperature alone, got {}",
            atmosphere_900.delta.drop_m
        );
    }

    /// Review round 3, N2 (derived-value exclusion rule, Instance 3): the Altitude test above
    /// only exercises a QNH-vs-ABSOLUTE comparison, where `b`'s plain absolute pressure has no
    /// altitude dependency to leak in the first place -- it cannot catch Pressure itself still
    /// being swapped. This uses a QNH-vs-QNH comparison instead: both `b` variants share the
    /// SAME altitude (1200 m, always different from `a`'s 500 m -- review round 4, F1, the same
    /// reasoning as the Altitude test's own doc comment: comparing against an altitude that
    /// happens to match `a`'s would stop excluding `Pressure` in that leg and confound the
    /// comparison) but differ in their RAW QNH value, so their resolved `pressure_pa` differs
    /// too. If Pressure were still being copied (the bug this fix addresses), that resolved
    /// difference would move Atmosphere's contribution; since Pressure is excluded in BOTH
    /// comparisons here (both `b` variants sit at 1200 m, never 500 m), it must not.
    ///
    /// Unlike the Altitude test, this ONE still tolerates a small residual, not exact equality:
    /// `Altitude`/`Pressure` are excluded from being SWAPPED, but each `b` variant's baseline
    /// trajectory is still solved at ITS OWN real, unswapped (and, here, differing-by-QNH-value)
    /// pressure, and the swapped `Temperature` axis's effect on drop is itself (very slightly) a
    /// function of the ambient density it is evaluated in -- a genuine, small, second-order
    /// interaction between a held-fixed axis and a swapped one, not a leak (unlike the Altitude
    /// test, where `b`'s pressure is a FIXED absolute value that does not move with altitude at
    /// all, so there is no equivalent second-order channel there). I measured both sides with a
    /// temporary, uncommitted probe before choosing this test's tolerance: with the fix
    /// disabled, this fixture's two comparisons differed by `~9.7e-3` (a first-order effect, the
    /// shape a real leak takes). With the fix applied, they differ by `~1.8e-4` -- over 50x
    /// smaller. `2e-3` sits with comfortable margin above the fixed residual and comfortably
    /// below the broken one.
    ///
    /// See `pressure_and_altitude_exclusion_give_an_exactly_zero_atmosphere_contribution_when_nothing_else_differs`
    /// below for a complementary, bit-EXACT version of this same property (the reviewer's own
    /// suggestion): when the two requests differ ONLY in altitude, every remaining Atmosphere
    /// axis is identical, so there is no second-order channel left and the contribution must be
    /// precisely `0.0`.
    #[test]
    fn pressure_exclusion_does_not_leak_a_partial_effect_into_atmosphere() {
        let qnh_json = |altitude_m: f64, temperature_k: f64, raw_qnh_pa: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0},
                "atmosphere": {"altitude_m": altitude_m, "temperature_k": temperature_k,
                               "pressure_pa": raw_qnh_pa, "pressure_reference": "qnh"},
                "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json(500.0, 288.0, 101_325.0))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b_qnh_101325 = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json(1200.0, 300.0, 101_325.0))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b_qnh_105000 = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json(1200.0, 300.0, 105_000.0))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;

        // Fixture assumptions: both b variants sit at the SAME altitude (1200 m), which must
        // NOT equal a's (500 m) -- see the doc comment above -- and their resolved pressure_pa
        // must differ, purely from the raw QNH value, confirming this fixture actually
        // exercises "does the specific excluded Pressure value leak."
        assert_eq!(b_qnh_101325.atmosphere.altitude_m, b_qnh_105000.atmosphere.altitude_m);
        assert_ne!(b_qnh_101325.atmosphere.altitude_m, a.atmosphere.altitude_m);
        assert_ne!(
            b_qnh_101325.atmosphere.pressure_pa, b_qnh_105000.atmosphere.pressure_pa,
            "fixture assumption: different raw QNH at the same altitude must resolve to a \
             different station pressure"
        );

        let rep_101325 = explain_difference(&a, &b_qnh_101325, &[300.0]).unwrap();
        let rep_105000 = explain_difference(&a, &b_qnh_105000, &[300.0]).unwrap();

        for rep in [&rep_101325, &rep_105000] {
            for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
                assert!(
                    rep.skipped_axes.iter().any(|s| s.group == InputGroup::Atmosphere
                        && s.axis == InputAxis::Pressure
                        && s.direction == direction),
                    "expected Pressure excluded on {direction:?} under QNH-vs-QNH; got {:?}",
                    rep.skipped_axes
                );
            }
        }

        let atmosphere_101325 = rep_101325.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        let atmosphere_105000 = rep_105000.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        // See the doc comment above for exactly why this is 2e-3, not exact equality: measured
        // empirically at ~1.8e-4 with the fix applied vs. ~9.7e-3 with it disabled, on this
        // exact fixture.
        assert!(
            (atmosphere_101325.delta.drop_m - atmosphere_105000.delta.drop_m).abs() < 2e-3,
            "Atmosphere's contribution changed by more than the expected small second-order \
             residual when b's raw QNH value changed (101325 -> 105000 Pa) at the SAME altitude, \
             even though Pressure and Altitude are both excluded -- {} vs {} -- a first-order \
             pressure effect looks like it is leaking through",
            atmosphere_101325.delta.drop_m,
            atmosphere_105000.delta.drop_m
        );
        assert!(
            atmosphere_105000.delta.drop_m.abs() > 0.001,
            "expected a non-negligible Atmosphere contribution from temperature alone, got {}",
            atmosphere_105000.delta.drop_m
        );
    }

    /// The reviewer's own suggestion (review round 4): a bit-EXACT complement to the tolerance
    /// test above. Two QNH-referenced requests differing ONLY in altitude (same explicit
    /// temperature, same everything else) exclude both `Altitude` and `Pressure`, leaving every
    /// OTHER `Atmosphere` axis byte-identical between them -- so `swap_group` reproduces the
    /// destination exactly, both legs are exact no-ops, and the contribution must be precisely
    /// `0.0`, not merely small. This is the same shape as
    /// `wind_direction_is_excluded_under_compass_wind_even_with_an_identical_raw_bearing` and the
    /// post-C1 `ZeroSightGeometry` assertion in `every_groups_contribution_matches_an_independent_recomputation`:
    /// a real leak would produce a large, first-order non-zero here, not a rounding-sized one.
    #[test]
    fn pressure_and_altitude_exclusion_give_an_exactly_zero_atmosphere_contribution_when_nothing_else_differs(
    ) {
        let qnh_json = |altitude_m: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0},
                "atmosphere": {"altitude_m": altitude_m, "temperature_k": 288.0,
                               "pressure_pa": 101_325.0, "pressure_reference": "qnh"},
                "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json(500.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&qnh_json(1200.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert_ne!(a.atmosphere.altitude_m, b.atmosphere.altitude_m, "fixture assumption");
        assert_eq!(
            a.atmosphere.temperature_k, b.atmosphere.temperature_k,
            "fixture assumption: temperature must be identical, or Atmosphere would have a \
             real, non-excluded axis to report and the contribution would not be zero"
        );

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("Altitude/Pressure exclusion must not abort the comparison");
        let atmosphere = rep.rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        assert_eq!(
            atmosphere.delta.drop_m, 0.0,
            "Atmosphere's drop contribution must be exactly zero when altitude is the ONLY \
             thing that differs and both Altitude and Pressure are excluded, got {}",
            atmosphere.delta.drop_m
        );
        assert_eq!(
            atmosphere.delta.windage_m, 0.0,
            "Atmosphere's windage contribution must be exactly zero for the same reason, got {}",
            atmosphere.delta.windage_m
        );

        // Non-triviality: unlike the other exact-zero tests in this file, a nonzero TOTAL isn't
        // the right guard here -- altitude alone, once excluded, is expected to affect nothing
        // at all, so a zero total is normal, not suspicious. What WOULD be suspicious is the
        // exclusion mechanism never actually engaging (e.g. a refactor that silently turned
        // plan_exclusions into a no-op, which would ALSO report an exact zero here, for the
        // wrong reason). Guard against that directly: both axes must actually appear in
        // skipped_axes.
        for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
            for axis in [InputAxis::Altitude, InputAxis::Pressure] {
                assert!(
                    rep.skipped_axes.iter().any(|s| s.group == InputGroup::Atmosphere
                        && s.axis == axis
                        && s.direction == direction),
                    "expected {axis:?} excluded on {direction:?}; got {:?}",
                    rep.skipped_axes
                );
            }
        }
    }

    /// Review round 4, F1: the Pressure exclusion must NOT fire when the two requests'
    /// altitudes MATCH -- only a differing altitude actually breaks the QNH-to-station-pressure
    /// pairing. Two QNH-referenced requests at the SAME altitude but DIFFERENT raw altimeter
    /// settings must attribute that difference to `Atmosphere` normally: excluding `Pressure`
    /// unconditionally on the QNH mode flag alone (the bug this test guards against) would
    /// UNDER-attribute a genuine altimeter-setting difference into the interaction remainder on
    /// the most natural same-location QNH comparison anyone would run.
    #[test]
    fn pressure_is_still_swapped_when_the_altitude_matches() {
        let build = |raw_qnh_pa: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0},
                "atmosphere": {"altitude_m": 500.0, "temperature_k": 288.0,
                               "pressure_pa": raw_qnh_pa, "pressure_reference": "qnh"},
                "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(101_325.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(98_000.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert_eq!(
            a.atmosphere.altitude_m, b.atmosphere.altitude_m,
            "fixture assumption: identical altitude on both sides"
        );
        assert_ne!(
            a.atmosphere.pressure_pa, b.atmosphere.pressure_pa,
            "fixture assumption: different raw QNH must still resolve to different station \
             pressure at the same altitude"
        );

        let rep = explain_difference(&a, &b, &[300.0]).unwrap();

        assert!(
            !rep.skipped_axes.iter().any(|s| s.group == InputGroup::Atmosphere
                && s.axis == InputAxis::Pressure),
            "Pressure must not be excluded when altitude matches on both sides; got {:?}",
            rep.skipped_axes
        );

        let atmosphere = rep
            .rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Atmosphere)
            .unwrap();
        assert!(
            atmosphere.delta.drop_m.abs() > 0.001,
            "Atmosphere's contribution must be substantial and non-zero here -- the raw QNH \
             genuinely differs and the shared altitude means Pressure is safe to swap -- got {}",
            atmosphere.delta.drop_m
        );
    }

    /// Review I2: an ordinarily-optional axis present on exactly one of the two requests splits
    /// the two legs exactly like a `with_axis` refusal does -- forward would keep A's own value
    /// (nothing to copy from B) while backward would overwrite B with A's, silently
    /// half-attributing whatever that value's effect is, with NO record at all. `latitude_rad`
    /// is a convenient example: present here on `a`, absent on `b` -- an ordinary difference
    /// between two saved profiles, nothing to do with Coriolis specifically (Coriolis is left
    /// disabled on both, so this is not entangled with `resolve_effects`'s separate requirement
    /// that Coriolis needs a latitude on any request that enables it).
    #[test]
    fn an_axis_present_on_only_one_side_is_skipped_and_recorded_symmetrically() {
        let with_latitude_json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {"temperature_k": 288.0, "latitude_rad": 0.7},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&with_latitude_json).unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert_eq!(
            a.atmosphere.latitude_rad,
            Some(0.7),
            "fixture assumption: a supplies latitude_rad"
        );

        let b = resolved(870.0, 300.0);
        assert_eq!(
            b.atmosphere.latitude_rad, None,
            "fixture assumption: b omits latitude_rad entirely"
        );

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("a presence-only-on-one-side axis must not abort the whole comparison");

        for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
            assert!(
                rep.skipped_axes.iter().any(|s| s.group == InputGroup::Atmosphere
                    && s.axis == InputAxis::Latitude
                    && s.direction == direction),
                "expected a {direction:?} Latitude skip when it is present on only one side; \
                 got {:?}",
                rep.skipped_axes
            );
        }
    }

    /// This fixture's shape comes from the ordering hazard documented in this module's doc
    /// comment: `ShotGeometry` lists `TargetDistance`, `ShootingAngle` and `Cant` before
    /// `ShotAzimuth`. Deciding exclusions from `a`/`b` directly (`plan_exclusions`, never a
    /// request `swap_group` has already partly rewritten) makes that hazard structurally
    /// impossible now -- there is no accumulated, partly-laundered request for a check to see,
    /// regardless of axis order -- so this specific fixture can no longer directly exercise an
    /// order-dependent failure the way an earlier, per-swap_group-call revision of this module
    /// could. What it still usefully pins down is review I1's symmetric-exclusion fix: `b`'s
    /// wind is shooter-relative, so ONLY the Forward leg independently refuses `ShotAzimuth` (a
    /// compass-referenced `a`), yet Backward must be excluded too, with a reason that still
    /// names the compass refusal that drove it.
    #[test]
    fn shot_azimuth_is_refused_even_when_other_shot_geometry_axes_are_applied_first() {
        let compass_json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "shooting_angle_rad": 0.05, "cant_angle_rad": 0.01,
                     "shot_azimuth_rad": 0.3},
            "atmosphere": {"temperature_k": 288.0},
            "wind": {"speed_mps": 3.0, "direction_from_rad": 1.0, "wind_reference": "compass"},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&compass_json).unwrap(),
        )
        .unwrap()
        .resolved_request;
        match &a.wind {
            crate::solve_json::ResolvedWindV1::Constant(c) => assert_eq!(
                c.wind_reference,
                Some(crate::solve_json::WindReferenceV1::Compass),
                "fixture assumption: a's wind must be compass-referenced"
            ),
            crate::solve_json::ResolvedWindV1::Segmented(_) => panic!("constant wind expected"),
        }

        let shooter_json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 950.0, "shooting_angle_rad": 0.08, "cant_angle_rad": 0.02,
                     "shot_azimuth_rad": 0.9},
            "atmosphere": {"temperature_k": 288.0},
            "wind": {"speed_mps": 3.0, "direction_from_rad": 1.0},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&shooter_json).unwrap(),
        )
        .unwrap()
        .resolved_request;

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("a refused axis must not abort the whole comparison");

        let hit = rep
            .skipped_axes
            .iter()
            .find(|s| {
                s.group == InputGroup::ShotGeometry
                    && s.axis == InputAxis::ShotAzimuth
                    && s.direction == SwapDirectionV1::Forward
            })
            .unwrap_or_else(|| {
                panic!(
                    "expected a Forward-direction ShotAzimuth skip under compass wind, even \
                     though TargetDistance/ShootingAngle/Cant are applied (and each \
                     re-resolved) first within the same ShotGeometry group; got skipped_axes = \
                     {:?}",
                    rep.skipped_axes
                )
            });
        assert!(
            hit.reason.to_lowercase().contains("compass"),
            "reason should name compass wind: {}",
            hit.reason
        );
        // Review I1: even though `b` is shooter-relative and would not refuse ShotAzimuth on
        // its own, Backward must ALSO be excluded -- otherwise forward keeps A's azimuth
        // (refused) while backward hands B the whole of A's azimuth (not refused), splitting
        // the two legs into different counterfactuals.
        let backward_hit = rep
            .skipped_axes
            .iter()
            .find(|s| s.group == InputGroup::ShotGeometry
                && s.axis == InputAxis::ShotAzimuth
                && s.direction == SwapDirectionV1::Backward)
            .unwrap_or_else(|| {
                panic!(
                    "expected a Backward-direction ShotAzimuth skip too (review I1: symmetric \
                     exclusion), even though b's wind is shooter-relative; got skipped_axes = \
                     {:?}",
                    rep.skipped_axes
                )
            });
        assert!(
            backward_hit.reason.to_lowercase().contains("compass"),
            "the sympathetic Backward exclusion's reason should still explain the compass \
             refusal that drove it: {}",
            backward_hit.reason
        );
    }

    /// Review round 3, N1 (derived-value exclusion rule, Instance 2): `with_axis` already
    /// refuses to swap `ShotAzimuth` under compass wind, protecting `ShotGeometry`'s own
    /// number, but nothing protected `WindDirection`, whose OWN resolved (shooter-relative)
    /// value is `(bearing - shot_azimuth_rad).rem_euclid(2*pi)` -- also a function of
    /// `shot_azimuth_rad`. Two compass-referenced requests with the SAME raw wind bearing but
    /// DIFFERENT `shot_azimuth_rad` resolve to DIFFERENT `direction_from_rad`, purely from the
    /// azimuth difference; `Wind`'s contribution must still be exactly zero, not a re-labelled
    /// azimuth effect.
    #[test]
    fn wind_direction_is_excluded_under_compass_wind_even_with_an_identical_raw_bearing() {
        let build = |shot_azimuth_rad: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0, "shot_azimuth_rad": shot_azimuth_rad},
                "atmosphere": {"temperature_k": 288.0},
                "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2,
                         "wind_reference": "compass"},
                "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(0.0)).unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(std::f64::consts::FRAC_PI_4))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;

        // Sanity: the SAME raw bearing (both pi/2) must resolve to DIFFERENT shooter-relative
        // directions, purely from the azimuth difference -- confirming this fixture actually
        // exercises the hazard, not a fixture that happens to agree for an unrelated reason.
        let (a_dir, b_dir) = match (&a.wind, &b.wind) {
            (
                crate::solve_json::ResolvedWindV1::Constant(ca),
                crate::solve_json::ResolvedWindV1::Constant(cb),
            ) => (ca.direction_from_rad, cb.direction_from_rad),
            _ => panic!("constant wind expected on both sides"),
        };
        assert!(
            (a_dir - b_dir).abs() > 0.1,
            "fixture must produce genuinely different resolved wind directions from the SAME \
             raw bearing (a: {a_dir}, b: {b_dir}), or this test proves nothing"
        );

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("a derived-value exclusion must not abort the whole comparison");

        // Non-triviality (small consistency note, review round 4): the total must be genuinely
        // non-zero, or a fixture that collapsed to two identical trajectories overall would also
        // report Wind == 0.0 for the wrong reason (matching the pattern the Pressure tests
        // already use). This fixture's total IS expected to differ: a's and b's OWN effective
        // (shooter-relative) wind angle differ, since the SAME raw bearing is transformed by
        // DIFFERENT shot azimuths -- the sanity check above already pins that resolved-direction
        // difference down; this confirms it actually shows up in the un-swapped baseline too.
        assert!(
            rep.rows[0].total.windage_m.abs() > 0.01,
            "expected a's and b's own (un-swapped) trajectories to differ meaningfully in \
             windage, since their effective wind angles differ -- got total windage {}",
            rep.rows[0].total.windage_m
        );

        let wind = rep
            .rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Wind)
            .unwrap();
        assert_eq!(
            wind.delta.drop_m, 0.0,
            "Wind's drop contribution must be exactly zero (WindDirection excluded), got {}",
            wind.delta.drop_m
        );
        assert_eq!(
            wind.delta.windage_m, 0.0,
            "Wind's windage contribution must be exactly zero (WindDirection excluded), got {}",
            wind.delta.windage_m
        );

        for direction in [SwapDirectionV1::Forward, SwapDirectionV1::Backward] {
            let hit = rep
                .skipped_axes
                .iter()
                .find(|s| {
                    s.group == InputGroup::Wind
                        && s.axis == InputAxis::WindDirection
                        && s.direction == direction
                })
                .unwrap_or_else(|| {
                    panic!(
                        "expected a {direction:?} WindDirection skip under compass wind; got \
                         {:?}",
                        rep.skipped_axes
                    )
                });
            assert!(
                hit.reason.to_lowercase().contains("shot_azimuth")
                    || hit.reason.to_lowercase().contains("azimuth"),
                "reason should name the shot-azimuth dependency: {}",
                hit.reason
            );
        }
    }

    /// Review round 4, F1: the WindDirection exclusion must NOT fire when the two requests'
    /// shot azimuths MATCH -- only a differing azimuth actually entangles the compass-to-
    /// shooter-relative transform. Two compass-referenced requests with the SAME
    /// `shot_azimuth_rad` but DIFFERENT raw wind bearings must attribute that difference to
    /// `Wind` normally: excluding `WindDirection` unconditionally on the compass mode flag
    /// alone (the bug this test guards against) would UNDER-attribute a genuine wind-bearing
    /// difference into the interaction remainder on the most natural same-shot-direction
    /// compass comparison anyone would run.
    #[test]
    fn wind_direction_is_still_swapped_when_the_shot_azimuth_matches() {
        let build = |bearing_rad: f64| {
            serde_json::json!({
                "schema_version": 1,
                "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                               "ballistic_coefficient": 0.243},
                "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
                "shot": {"max_range_m": 900.0, "shot_azimuth_rad": 0.3},
                "atmosphere": {"temperature_k": 288.0},
                "wind": {"speed_mps": 3.0, "direction_from_rad": bearing_rad,
                         "wind_reference": "compass"},
                "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
            })
            .to_string()
        };
        let a = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(std::f64::consts::FRAC_PI_2))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&build(std::f64::consts::FRAC_PI_4))
                .unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert_eq!(
            a.shot.shot_azimuth_rad, b.shot.shot_azimuth_rad,
            "fixture assumption: identical shot azimuth on both sides"
        );

        let rep = explain_difference(&a, &b, &[300.0]).unwrap();

        assert!(
            !rep.skipped_axes.iter().any(|s| s.group == InputGroup::Wind
                && s.axis == InputAxis::WindDirection),
            "WindDirection must not be excluded when the shot azimuth matches on both sides; \
             got {:?}",
            rep.skipped_axes
        );

        let wind = rep
            .rows[0]
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::Wind)
            .unwrap();
        assert!(
            wind.delta.windage_m.abs() > 0.01,
            "Wind's windage contribution must be substantial and non-zero here -- the raw \
             bearing genuinely differs and the shared azimuth means WindDirection is safe to \
             swap -- got {}",
            wind.delta.windage_m
        );
    }

    /// Requirement (1)'s third example, both detection paths at once: `a` is constant wind and
    /// `b` is segmented, so the Forward leg (reading FROM b) hits `read_axis` returning `None`,
    /// and the Backward leg (writing TO b) hits `with_axis`'s own `AxisAbsent`.
    #[test]
    fn wind_axes_are_skipped_and_recorded_under_segmented_wind_on_either_side() {
        let a = resolved(823.0, 288.0);
        let segmented_json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 850.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {"temperature_k": 295.0},
            "wind": {"segments": [{"until_distance_m": 900.0, "speed_mps": 4.0,
                                    "direction_from_rad": 1.2}]},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let b = crate::solve_v1::solve_v1(
            crate::solve_json::decode_solve_request_v1(&segmented_json).unwrap(),
        )
        .unwrap()
        .resolved_request;
        assert!(matches!(
            b.wind,
            crate::solve_json::ResolvedWindV1::Segmented(_)
        ));

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("segmented wind must not abort the whole comparison");

        for axis in [
            InputAxis::WindSpeed,
            InputAxis::WindDirection,
            InputAxis::WindVertical,
        ] {
            assert!(
                rep.skipped_axes.iter().any(|s| s.group == InputGroup::Wind
                    && s.axis == axis
                    && s.direction == SwapDirectionV1::Forward),
                "{axis:?}: expected a Forward skip (the source, b, is segmented); got {:?}",
                rep.skipped_axes
            );
            assert!(
                rep.skipped_axes.iter().any(|s| s.group == InputGroup::Wind
                    && s.axis == axis
                    && s.direction == SwapDirectionV1::Backward),
                "{axis:?}: expected a Backward skip (the destination, b, is segmented); got {:?}",
                rep.skipped_axes
            );
        }
    }

    /// A genuine, non-refusal solve failure caused by the counterfactual construction itself
    /// must propagate, never be silently caught like a `with_axis` refusal. `Effects` lists
    /// `MagnusEnabled` before `EnhancedSpinDriftEnabled`; swapping Magnus on (from B into A)
    /// while A still carries its own `enhanced_spin_drift = true` from before produces an
    /// intermediate request with both flags true, which `solve_v1` rejects as `ConflictingFields`
    /// (`taxonomy.rs`'s Known Limitation (d)) -- even though the FINAL state after the whole
    /// group is applied (magnus on, enhanced spin drift off) would have been perfectly valid.
    #[test]
    fn a_genuine_solve_conflict_from_group_swap_propagates_not_silently_skipped() {
        let a = resolved_with_effects(823.0, false, true);
        let b = resolved_with_effects(823.0, true, false);
        let e = explain_difference(&a, &b, &[300.0]);
        match e {
            Err(KernelError::Solve { code, .. }) => {
                assert_eq!(code, crate::solve_json::SolveErrorCodeV1::ConflictingFields);
            }
            other => panic!(
                "expected a genuine ConflictingFields solve error to propagate, got {other:?}"
            ),
        }
    }
}
