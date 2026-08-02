//! MBA-1345: explain why two fully resolved solutions differ.
//!
//! For each input group the difference is attributed by a SYMMETRIC counterfactual: swap the
//! group from B into A, and independently from A into B, and average. That makes the answer
//! independent of replacement order, which a one-directional swap is not. Whatever the group
//! contributions fail to explain is reported as an explicit interaction remainder and is never
//! distributed across groups: for correlated inputs there is no unique causal attribution, and
//! pretending otherwise is the failure this design exists to avoid.
//!
//! # Axes `with_axis` refuses, and axes that are simply absent
//!
//! Swapping a whole group means writing every one of its axes in turn (`axes_in_group`), and
//! [`with_axis`] can legitimately refuse an individual axis for a specific request:
//! `KernelError::AxisUnsupportedForRequest` for `Altitude` under a QNH-referenced atmosphere or
//! `ShotAzimuth` under compass-referenced wind (`access.rs`'s module doc), and
//! `KernelError::AxisAbsent` for the three wind axes when the request being written TO uses
//! segmented wind. [`read_axis`] separately returns `None` for the same three wind axes when
//! the request being read FROM is segmented, and also for a handful of ordinarily-optional
//! scalar fields (`length_m`, `latitude_rad`, the two zero-POI offsets, the lateral sight
//! offset, and -- true of the code though not named in `read_axis`'s own doc comment --
//! `zero_distance_m` on a request zeroed purely by an explicit angle) that just were not
//! supplied on this particular request.
//!
//! This module treats those two kinds of "nothing to copy" differently, deliberately:
//!
//! - An axis `with_axis` refuses, or one of the three wind axes with no scalar to read or write
//!   because either side uses segmented wind, is recorded in
//!   [`SolutionDiffReportV1::skipped_axes`] with the reason and which swap direction hit it,
//!   then left out of that one swap -- the rest of the group's axes are still attempted.
//!   Recording it matters: a silently skipped axis would let part of the real difference get
//!   folded into the interaction remainder with no indication that anything was left out, which
//!   is exactly the kind of unmarked gap this feature exists to prevent.
//! - An ordinarily-optional axis simply absent from the source request is skipped WITHOUT a
//!   report. There is nothing to attribute either way: the field does not exist on that
//!   request, the same as if the caller had never asked about it.
//!
//! Either way, a refusal or an absence never aborts the comparison -- only that one axis is
//! left out of that one group, for that one swap direction. A GENUINE failure (most notably a
//! `solve_v1` re-resolve failing partway through applying a group -- for instance the
//! pre-existing magnus/enhanced-spin-drift conflict noted in `taxonomy.rs`'s Known Limitation
//! (d), which a group swap can walk straight into by turning one flag on while the request
//! still carries the other from before) is NOT one of these two cases and propagates as an
//! error, uncaught, exactly as `central_difference` and `bisect_axis` already do for their own
//! non-refusal failures.
//!
//! # Why the refusal check reads the request as originally given, not the accumulated one
//!
//! A group with more than one axis applies its writes one at a time, re-resolving between them
//! (`swap_group` below, via a real `solve_v1` call) so each subsequent [`with_axis`] call sees
//! the previous write -- most importantly so a `requires_rezero` axis's cleared
//! `muzzle_angle_rad` is replaced by an ACTUALLY re-zeroed angle, not a placeholder, before the
//! next axis reads it back. That re-resolve goes through the reverse conversion
//! (`request_roundtrip.rs`), which ALWAYS clears `pressure_reference` and `wind_reference` back
//! to the omitted-field default (see that module's doc: the resolved values already have the
//! transform baked in, so echoing the mode back would apply it a second time). That is correct
//! for solving, but it means the very fact [`with_axis`]'s guards key on -- "was this request
//! originally QNH- or compass-referenced" -- is erased by the FIRST re-resolve within a group,
//! not preserved across it.
//!
//! `ShotGeometry` lists `ShotAzimuth` fourth (`TargetDistance`, `ShootingAngle`, `Cant`,
//! `ShotAzimuth`, `AimAzimuth`, `TargetHeight`): three unrelated axes, each triggering a
//! re-resolve, are applied before it. Checking the compass guard against the
//! progressively-modified request would have it silently pass by the time `ShotAzimuth` is
//! reached, building exactly the physically-inverted counterfactual the guard exists to
//! prevent -- with no error and no skip recorded, because from `with_axis`'s point of view at
//! that moment the request no longer looks compass-referenced at all. `swap_group` avoids this
//! by probing every axis's refusal against the group's `dst` argument exactly as passed in --
//! captured once, never itself round-tripped -- before applying the write to the accumulated
//! request. A refusal can only become MORE permissive as a group's axes are applied
//! (round-tripping only ever clears a reference-mode echo, it never introduces one), so an axis
//! that passes the probe is guaranteed to still be writable on the accumulated request. The test
//! `shot_azimuth_is_refused_even_when_other_shot_geometry_axes_are_applied_first` pins this down
//! by exercising exactly the four-axis ordering above.

use serde::{Deserialize, Serialize};

use crate::perturbation::access::{read_axis, with_axis, KernelError};
use crate::perturbation::taxonomy::{axes_in_group, InputAxis, InputGroup};
use crate::perturbation::{evaluate, kernel_solve_error, Observation};
use crate::solve_json::{ResolvedSolveRequestV1, SolveRequestV1};

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

/// One taxonomy axis that could not be carried across during a group's counterfactual swap,
/// and why -- see the module doc's "Axes `with_axis` refuses, and axes that are simply absent".
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
    /// Axes that could not be carried across for one or both swap directions of some group --
    /// see the module doc. Never distributed or hidden: whatever effect a skipped axis would
    /// have had is left inside `SolutionDiffRowV1::interaction_remainder`, not silently
    /// absorbed into that axis's own group.
    pub skipped_axes: Vec<SkippedAxisV1>,
    pub rows: Vec<SolutionDiffRowV1>,
}

/// Build the [`SolveRequestV1`] representing `dst` with every axis of `group` copied from
/// `src`, applying each axis one at a time and re-resolving between writes (via a real
/// `solve_v1` call) so a later axis in the same group sees the effect of an earlier one. Does
/// not itself evaluate a trajectory -- the caller does that once, over the whole `ranges_m`
/// slice, exactly as [`central_difference`](crate::perturbation::central_difference) already
/// does for its own pair of counterfactual requests.
///
/// Every axis's refusal is probed against `dst` exactly as passed in here, never the
/// accumulated running request -- see the module doc's "Why the refusal check reads the request
/// as originally given" section for why checking the accumulated request instead would silently
/// miss a refusal that should fire once an earlier axis in the SAME group has triggered a
/// re-resolve.
///
/// A refused or structurally-absent axis (`with_axis`'s `AxisUnsupportedForRequest`/
/// `AxisAbsent`, or `read_axis` returning `None` for one of the three wind axes because the
/// SOURCE request uses segmented wind) is pushed onto `skipped`, tagged with `group` and
/// `direction`, and left out of this swap -- the running request simply keeps `dst`'s own value
/// for that one axis. Every other axis in `group` is still attempted. Any OTHER error
/// propagates via `?`: it is a genuine failure (a `read_axis`/`with_axis` mismatch that should
/// not be reachable here, or a `solve_v1` re-resolve failing for a reason unrelated to this
/// taxonomy, such as the magnus/enhanced-spin-drift conflict noted in `taxonomy.rs`'s Known
/// Limitation (d)), never a structurally unrepresentable axis, and must be reported, not
/// swallowed.
fn swap_group(
    dst: &ResolvedSolveRequestV1,
    src: &ResolvedSolveRequestV1,
    group: InputGroup,
    direction: SwapDirectionV1,
    skipped: &mut Vec<SkippedAxisV1>,
) -> Result<SolveRequestV1, KernelError> {
    let mut current = dst.clone();
    for &axis in axes_in_group(group) {
        let Some(v) = read_axis(src, axis) else {
            // The three wind axes read back `None` under segmented wind (taxonomy.rs Known
            // Limitation (c)) -- structurally absent, not merely unsupplied -- and that is the
            // ONLY reason read_axis ever returns None for one of these three specifically (see
            // access.rs: `wind_speed`/`wind_dir`/`wind_vert` are `Some` for every constant-wind
            // request, never `None` for any other reason). Worth recording, unlike every other
            // axis read_axis skips here, which is a genuinely-absent optional field with no
            // bearing on wind shape at all (module doc).
            if matches!(
                axis,
                InputAxis::WindSpeed | InputAxis::WindDirection | InputAxis::WindVertical
            ) {
                skipped.push(SkippedAxisV1 {
                    group,
                    axis,
                    direction,
                    reason: "the source request's wind is segmented; there is no single \
                             scalar value for this axis to copy"
                        .to_string(),
                });
            }
            continue;
        };
        // Probe against `dst` EXACTLY as first passed in, never the running `current` -- see
        // the module doc's "Why the refusal check reads the request as originally given".
        if let Err(e) = with_axis(dst, axis, v) {
            match e {
                KernelError::AxisUnsupportedForRequest { reason, .. } => {
                    skipped.push(SkippedAxisV1 {
                        group,
                        axis,
                        direction,
                        reason: reason.to_string(),
                    });
                    continue;
                }
                KernelError::AxisAbsent(_) => {
                    skipped.push(SkippedAxisV1 {
                        group,
                        axis,
                        direction,
                        reason: "the destination request's wind is segmented; there is no \
                                 single scalar field for this axis to overwrite"
                            .to_string(),
                    });
                    continue;
                }
                other => return Err(other),
            }
        }
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
/// direction), plus two more for `a`/`b` themselves -- 16 solves total for the fixed
/// seven-group taxonomy, each covering every range in `ranges_m` at once rather than once per
/// range (see [`evaluate`]'s own cost note). A `requires_rezero` axis inside a multi-axis group
/// additionally costs one re-resolve per axis while that group's own counterfactual request is
/// being built (`swap_group`) -- unavoidable here, since only a real solve produces the
/// correctly re-zeroed angle the next axis in the same group needs to start from.
///
/// # Axes that cannot be swapped
///
/// See the module doc for the full explanation. In short: [`with_axis`] can refuse an axis for
/// a structural reason specific to `a` or `b` (a QNH-referenced altitude, a compass-referenced
/// shot azimuth, one of the three wind axes under segmented wind). A refused or structurally
/// absent axis is recorded in the returned report's `skipped_axes`, and the rest of that axis's
/// group is still attributed normally -- a refusal never aborts the comparison. An
/// ordinarily-optional axis that is simply absent from `a` or `b` (no `length_m` supplied, for
/// instance) is skipped without being recorded: there is nothing there to attribute either way.
///
/// # Errors
///
/// Returns whatever [`evaluate`] or [`with_axis`] error on `a`, `b`, or either swapped
/// counterfactual, MINUS the two refusal cases above (`KernelError::AxisUnsupportedForRequest`,
/// `KernelError::AxisAbsent`), which are caught and recorded rather than returned. Every other
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
        // Each swapped request is evaluated exactly ONCE over the whole `ranges_m` slice below
        // -- never once per range -- so N ranges do not multiply the cost of a
        // `requires_rezero` axis's zero search inside `swap_group`.
        let fwd_req = swap_group(a, b, group, SwapDirectionV1::Forward, &mut skipped_axes)?;
        let bwd_req = swap_group(b, a, group, SwapDirectionV1::Backward, &mut skipped_axes)?;
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
             left out of its group's contribution for that direction; any real effect it would \
             have had is folded into the interaction remainder, not silently dropped."
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
    #[test]
    fn contribution_is_the_mean_of_a_genuinely_different_forward_and_backward() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let ranges = [600.0];

        let obs_a = evaluate(&(&a).into(), &ranges).unwrap();
        let obs_b = evaluate(&(&b).into(), &ranges).unwrap();
        let mut skipped = Vec::new();
        let fwd_req = swap_group(
            &a,
            &b,
            InputGroup::MuzzleVelocity,
            SwapDirectionV1::Forward,
            &mut skipped,
        )
        .unwrap();
        let bwd_req = swap_group(
            &b,
            &a,
            InputGroup::MuzzleVelocity,
            SwapDirectionV1::Backward,
            &mut skipped,
        )
        .unwrap();
        let fwd_obs = evaluate(&fwd_req, &ranges).unwrap();
        let bwd_obs = evaluate(&bwd_req, &ranges).unwrap();
        let forward = fwd_obs[0].drop_m - obs_a[0].drop_m;
        let backward = obs_b[0].drop_m - bwd_obs[0].drop_m;

        assert!(
            (forward - backward).abs() > 1e-6,
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
    /// (same primitives -- `swap_group` + `evaluate` -- but the mean is taken here, outside
    /// `explain_difference`'s own code path) and the remainder against the same total-minus-sum
    /// formula computed independently. This would fail under: a forward-only computation (the
    /// independent mean would disagree with a forward-only report), a remainder silently zeroed
    /// or redistributed into the groups (the independent remainder would disagree with a
    /// tampered report), or two groups' contributions swapped (checked BY GROUP IDENTITY, and
    /// MuzzleVelocity/Atmosphere are the only two groups with a non-negligible true contribution
    /// in this fixture, so a swap involving either is a large, unmissable mismatch).
    #[test]
    fn every_groups_contribution_matches_an_independent_recomputation() {
        let a = resolved(823.0, 288.0);
        let b = resolved(870.0, 300.0);
        let ranges = [600.0];
        let rep = explain_difference(&a, &b, &ranges).unwrap();
        let row = &rep.rows[0];

        let obs_a = evaluate(&(&a).into(), &ranges).unwrap();
        let obs_b = evaluate(&(&b).into(), &ranges).unwrap();

        let mut independent_sum = DeltaV1::default();
        for &group in InputGroup::ALL {
            let mut skipped = Vec::new();
            let fwd_req =
                swap_group(&a, &b, group, SwapDirectionV1::Forward, &mut skipped).unwrap();
            let bwd_req =
                swap_group(&b, &a, group, SwapDirectionV1::Backward, &mut skipped).unwrap();
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

        // `resolved()` only ever varies `mv` and `temp`, so a group whose axes are entirely
        // independent of both -- ProjectileDrag, Wind, ShotGeometry, Effects -- must come back
        // EXACTLY zero (nothing to swap, since every one of its axis values is identical on a
        // and b). MuzzleVelocity and Atmosphere differ directly; ZeroSightGeometry is NOT
        // independent despite `sight_height_m`/`zero_distance_m` etc. being identical, because
        // it also carries `MuzzleAngle` -- the RESOLVED, already-zeroed elevation, which is a
        // function of muzzle velocity and air density and therefore genuinely differs between a
        // and b even though every explicit zero/sight knob does not. All three -- MuzzleVelocity,
        // Atmosphere, ZeroSightGeometry -- must be non-negligible AND pairwise distinguishable in
        // magnitude, so a label swap among any pair of the seven groups would be an unmissable
        // mismatch against the per-group check above, not something these sanity bounds could
        // mask by coincidence.
        for g in [
            InputGroup::ProjectileDrag,
            InputGroup::Wind,
            InputGroup::ShotGeometry,
            InputGroup::Effects,
        ] {
            let c = row.contributions.iter().find(|x| x.group == g).unwrap();
            assert!(
                c.delta.drop_m.abs() < 1e-9,
                "{g:?} does not differ between a and b in this fixture at all, expected exactly \
                 0, got {}",
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
        let zero_sight_geometry = row
            .contributions
            .iter()
            .find(|c| c.group == InputGroup::ZeroSightGeometry)
            .unwrap();
        for (name, c) in [
            ("MuzzleVelocity", mv),
            ("Atmosphere", atmosphere),
            ("ZeroSightGeometry", zero_sight_geometry),
        ] {
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
        assert!(
            (mv.delta.drop_m - zero_sight_geometry.delta.drop_m).abs() > 0.01,
            "MuzzleVelocity ({}) and ZeroSightGeometry ({}) should be distinguishable",
            mv.delta.drop_m,
            zero_sight_geometry.delta.drop_m
        );

        // The remainder itself must also be genuinely nonzero here (real second-order
        // interaction between muzzle velocity and temperature/elevation) -- not identically
        // zero, and not silently redistributed into the groups above (which would make THIS
        // remainder read back as zero while the groups' sum still matched `total`).
        assert!(
            row.interaction_remainder.drop_m.abs() > 0.01,
            "expected a non-negligible interaction remainder in this fixture, got {}",
            row.interaction_remainder.drop_m
        );
    }

    /// Requirement (1): `with_axis` refuses `Altitude` under a QNH-referenced atmosphere. The
    /// refusal must be recorded, not silently dropped, and must not abort the comparison.
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

        let rep = explain_difference(&a, &b, &[300.0])
            .expect("a refused axis must not abort the whole comparison");
        let hit = rep
            .skipped_axes
            .iter()
            .find(|s| {
                s.group == InputGroup::Atmosphere
                    && s.axis == InputAxis::Altitude
                    && s.direction == SwapDirectionV1::Forward
            })
            .unwrap_or_else(|| {
                panic!(
                    "expected a Forward-direction Altitude skip under QNH pressure, got {:?}",
                    rep.skipped_axes
                )
            });
        assert!(
            hit.reason.to_lowercase().contains("qnh"),
            "reason should name QNH: {}",
            hit.reason
        );
    }

    /// Requirement (1) plus the ordering hazard documented in this module's doc comment:
    /// `ShotGeometry` lists `TargetDistance`, `ShootingAngle` and `Cant` before `ShotAzimuth`,
    /// each triggering a re-resolve that launders away the compass-wind echo `with_axis`'s
    /// guard depends on. If the refusal were probed against the progressively-modified request
    /// instead of `dst` as originally given, this specific ordering would silently defeat the
    /// guard -- ShotAzimuth would swap normally, building a physically-inverted counterfactual
    /// with no error and no skip. This test would fail under exactly that regression.
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
        // The backward leg's destination (b) is shooter-relative, so ShotAzimuth must swap
        // normally there -- no Backward-direction skip for this axis.
        assert!(
            !rep.skipped_axes.iter().any(|s| s.axis == InputAxis::ShotAzimuth
                && s.direction == SwapDirectionV1::Backward),
            "did not expect a Backward-direction ShotAzimuth skip: b is shooter-relative"
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
