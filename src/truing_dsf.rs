//! Mach-keyed drop-scale-factor (DSF) truing table (MBA-1357).
//!
//! Applied Ballistics' published two-stage truing workflow calibrates muzzle velocity
//! (MV) first — that fixes the supersonic drag curve against a chronograph/observed-drop
//! comparison at Mach >= 1.2. Below that, as the bullet moves through the transonic
//! region and into the subsonic regime, no single MV correction can fix drop
//! discrepancies that grow with range: the residual is a slowly-varying function of
//! Mach, not a constant offset. AB's second stage records a handful of *observed drop /
//! predicted drop* ratios at specific (subsonic-or-transonic) Mach numbers and uses them
//! to scale predicted drop at nearby Mach numbers on later solves.
//!
//! This module is a **cleanroom reimplementation** of that workflow's *shape*, not a
//! bit-for-bit copy of AB's unpublished interpolation — Kestrel/AB do not publish their
//! exact curve. The design decision unique to this implementation is the **anchor**:
//! the table's Mach domain is `(0, 1.2)` (points at/above Mach 1.2 belong to MV truing,
//! not here — [`DsfTable::from_points`] rejects them), and every table implicitly
//! continues with a DSF of `1.0` at Mach 1.2 — the exact boundary where MV calibration
//! takes over. That implicit anchor point `(1.2, 1.0)` is never stored in
//! [`DsfTable::points`]; it exists only inside [`DsfTable::factor_at`]'s interpolation so
//! the transition from "supersonic, MV-trued, unscaled" to "transonic/subsonic,
//! DSF-scaled" is continuous — a shot solved at Mach 1.1999 and one solved at Mach 1.2001
//! get (to floating-point precision) the same drop. This is a functional-equivalence
//! choice made for this engine, not a replication of AB's internal method.
//!
//! Below the lowest recorded point, [`DsfTable::factor_at`] flat-clamps to that point's
//! DSF — there is no data past it, and AB's guidance is that further subsonic drop
//! continues to track the last-calibrated regime rather than drift back toward identity.
//!
//! [`apply_dsf`] is a **drop-only** post-processing step over an already-solved
//! [`crate::TrajectoryResult`]: it rescales each point's vertical position relative to
//! the line of sight by the DSF at that point's Mach, and touches nothing else —
//! velocity, kinetic energy, time, and downrange/windage position are byte-identical
//! before and after. Per-point Mach is computed the same way the solver's own
//! diagnostics compute it (see [`apply_dsf`]'s doc comment for the exact fields), NOT
//! from a re-derived per-altitude local speed of sound the engine does not store per
//! point.
//!
//! No feature gate: this module must compile for `wasm32-unknown-unknown`. It is
//! fs-free (profile persistence of a table's points is the caller's job, e.g.
//! `main.rs`'s saved-profile handling in a later task).

use serde::{Deserialize, Serialize};

use crate::cli_api::TrajectoryResult;

/// Upper bound (exclusive) of the Mach domain a [`DsfPoint`] may describe. Observations at
/// or above this Mach belong to muzzle-velocity truing, not the DSF table; it doubles as
/// the implicit anchor's Mach coordinate (`(DSF_MACH_CEILING, 1.0)`) in
/// [`DsfTable::factor_at`].
pub const DSF_MACH_CEILING: f64 = 1.2;

/// DSF value of the implicit anchor at [`DSF_MACH_CEILING`] — identity, matching the
/// MV-trued supersonic regime this table hands off from.
pub const DSF_ANCHOR_VALUE: f64 = 1.0;

/// Exclusive lower bound a point's `dsf` must clear.
pub const DSF_MIN: f64 = 0.5;

/// Exclusive upper bound a point's `dsf` must clear.
pub const DSF_MAX: f64 = 2.0;

/// Maximum number of distinct points a [`DsfTable`] may hold.
pub const DSF_MAX_POINTS: usize = 6;

/// A new point within this many Mach units of an existing one supersedes it in
/// [`DsfTable::upsert`] instead of being appended.
pub const DSF_SUPERSEDE_TOLERANCE_MACH: f64 = 0.05;

/// One observed drop-scale-factor keyed to the Mach number it was recorded at.
///
/// `mach` must satisfy `0 < mach < 1.2`; `dsf` must be finite and satisfy
/// `0.5 < dsf < 2.0`. Both bounds are enforced by [`DsfTable::from_points`] and
/// [`DsfTable::upsert`] — this struct itself carries no invariant beyond the field types
/// (serde needs to deserialize arbitrary saved-profile content before it can be
/// validated).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DsfPoint {
    pub mach: f64,
    pub dsf: f64,
}

/// What [`DsfTable::upsert`] did with the incoming point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UpsertOutcome {
    /// No existing point was within [`DSF_SUPERSEDE_TOLERANCE_MACH`] Mach; the point was
    /// added as a new entry.
    Appended,
    /// An existing point was within [`DSF_SUPERSEDE_TOLERANCE_MACH`] Mach and was replaced.
    /// `old` is the point that was overwritten.
    Replaced { old: DsfPoint },
}

/// A validated, Mach-sorted table of up to [`DSF_MAX_POINTS`] [`DsfPoint`]s.
///
/// Construct via [`DsfTable::from_points`] (bulk, e.g. loading a saved profile) or by
/// starting from an empty table (`DsfTable::from_points(Vec::new())`, infallible) and
/// growing it with [`DsfTable::upsert`] (the `dsf` CLI verb's one-observation-at-a-time
/// path in a later task).
#[derive(Debug, Clone, PartialEq)]
pub struct DsfTable {
    /// Sorted ascending by `mach`. Never exceeds [`DSF_MAX_POINTS`] entries.
    points: Vec<DsfPoint>,
}

fn validate_point(point: &DsfPoint) -> Result<(), String> {
    if !point.mach.is_finite() || point.mach <= 0.0 || point.mach >= DSF_MACH_CEILING {
        return Err(format!(
            "DSF point Mach {} is out of range: must be finite and satisfy 0 < mach < {DSF_MACH_CEILING} \
             (observations at/above Mach {DSF_MACH_CEILING} belong to muzzle-velocity truing, not the DSF table)",
            point.mach
        ));
    }
    if !point.dsf.is_finite() || point.dsf <= DSF_MIN || point.dsf >= DSF_MAX {
        return Err(format!(
            "DSF value {} is out of range: must be finite and satisfy {DSF_MIN} < dsf < {DSF_MAX}",
            point.dsf
        ));
    }
    Ok(())
}

fn sort_by_mach(points: &mut [DsfPoint]) {
    points.sort_by(|a, b| {
        a.mach
            .partial_cmp(&b.mach)
            .expect("DsfPoint.mach is validated finite before insertion")
    });
}

/// Linear interpolation of `y` at `x`, between `(x0, y0)` and `(x1, y1)`.
fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if x1 == x0 {
        return y0;
    }
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

impl DsfTable {
    /// Validate, cap-check, and sort a bulk set of points (e.g. deserialized from a saved
    /// profile). Rejects any point failing validation (see [`DsfPoint`] bounds) and rejects more than
    /// [`DSF_MAX_POINTS`] points outright. Does NOT dedupe near-Mach points against each
    /// other — that supersede behavior is [`DsfTable::upsert`]'s job; a caller assembling
    /// points one at a time should use `upsert`, not construct a `Vec` with near-duplicates
    /// and pass it here.
    pub fn from_points(points: Vec<DsfPoint>) -> Result<DsfTable, String> {
        if points.len() > DSF_MAX_POINTS {
            return Err(format!(
                "DSF table supports at most {DSF_MAX_POINTS} points; got {} (remove one first, e.g. --clear-dsf)",
                points.len()
            ));
        }
        for point in &points {
            validate_point(point)?;
        }
        let mut sorted = points;
        sort_by_mach(&mut sorted);
        Ok(DsfTable { points: sorted })
    }

    /// Insert or supersede a point. A new point within [`DSF_SUPERSEDE_TOLERANCE_MACH`]
    /// Mach of an existing one replaces it (returning [`UpsertOutcome::Replaced`] with the
    /// overwritten point); otherwise it is appended (returning
    /// [`UpsertOutcome::Appended`]), unless the table is already at [`DSF_MAX_POINTS`]
    /// distinct points, in which case this errors naming the cap.
    pub fn upsert(&mut self, point: DsfPoint) -> Result<UpsertOutcome, String> {
        validate_point(&point)?;

        if let Some(existing) = self
            .points
            .iter_mut()
            .find(|p| (p.mach - point.mach).abs() <= DSF_SUPERSEDE_TOLERANCE_MACH)
        {
            let old = *existing;
            *existing = point;
            sort_by_mach(&mut self.points);
            return Ok(UpsertOutcome::Replaced { old });
        }

        if self.points.len() >= DSF_MAX_POINTS {
            return Err(format!(
                "DSF table already holds the maximum {DSF_MAX_POINTS} points; remove one first \
                 (e.g. --clear-dsf) before adding another"
            ));
        }

        self.points.push(point);
        sort_by_mach(&mut self.points);
        Ok(UpsertOutcome::Appended)
    }

    /// The drop-scale-factor at `mach`.
    ///
    /// - Identity (`1.0`) at or above [`DSF_MACH_CEILING`], and for an empty table at any
    ///   Mach — there is nothing to scale by.
    /// - Flat-clamped to the lowest point's `dsf` at or below the lowest point's Mach.
    /// - Piecewise-linear between successive points.
    /// - Piecewise-linear between the highest point and the implicit anchor
    ///   `(DSF_MACH_CEILING, DSF_ANCHOR_VALUE)` for Mach between the highest point and the
    ///   ceiling (this includes single-point tables, where "highest" and "lowest" are the
    ///   same point).
    pub fn factor_at(&self, mach: f64) -> f64 {
        if !mach.is_finite() || mach >= DSF_MACH_CEILING || self.points.is_empty() {
            return DSF_ANCHOR_VALUE;
        }

        let lowest = self.points[0];
        if mach <= lowest.mach {
            return lowest.dsf;
        }

        for pair in self.points.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            if mach <= hi.mach {
                return lerp(lo.mach, lo.dsf, hi.mach, hi.dsf, mach);
            }
        }

        // mach is above every explicit key (but still below the ceiling, checked above):
        // interpolate against the implicit anchor.
        let highest = *self.points.last().expect("checked non-empty above");
        lerp(
            highest.mach,
            highest.dsf,
            DSF_MACH_CEILING,
            DSF_ANCHOR_VALUE,
            mach,
        )
    }

    /// The table's points, sorted ascending by Mach.
    pub fn points(&self) -> &[DsfPoint] {
        &self.points
    }
}

/// Apply a DSF table to an already-solved trajectory, IN PLACE, scaling only each
/// point's drop below the line of sight — in BOTH `result.points` and, when present,
/// `result.sampled_points`.
///
/// For each `point` in `result.points`:
/// 1. Per-point Mach is `point.velocity_magnitude / result.station_speed_of_sound_mps` —
///    the same frozen "station" speed of sound the solver itself divides into
///    `velocity_magnitude` for its own per-point Mach diagnostics (Mach-transition
///    tracking, pitch-damping, precession/nutation; see the "MBA-1136 (rank 30)" comments
///    next to `resolved_atmosphere()` in `cli_api.rs`'s integration loops). The engine
///    does NOT store a re-derived per-altitude local speed of sound on `TrajectoryPoint`
///    itself, so this is "the way the solver does it", not a sea-level constant and not a
///    new per-point atmosphere recompute. It is also the same divisor the truing
///    observation path uses to derive an observation's Mach (`trajectory_observation.rs`),
///    so a DSF point keyed at derivation time lands back on the identical Mach at
///    application time — a per-point local recompute here would skew the two apart.
/// 2. `drop = result.line_of_sight_height_m - point.position.y` (drop below the
///    horizontal line of sight, in the solver's ground-referenced frame — the same
///    `drop_offset - y` convention `cli_api::fit_value_at` uses for BC-fit drop curves).
/// 3. `point.position.y` is rewritten so that the (possibly rescaled) drop is
///    `drop * table.factor_at(mach)`.
///
/// `result.sampled_points` (populated when `--sample-trajectory` is requested; read by
/// the Table's "Sampled Trajectory" section, CSV `--full`, and the PDF dope card — the
/// PDF dope card *always* requires sampling) is a SEPARATE `Vec<TrajectorySample>` from
/// `points` and was, until MBA-1357 Task 2's review (Critical #2), left untouched by this
/// function — those outputs silently rendered untrued drops even with an active DSF
/// table. Each `TrajectorySample` already stores its drop directly as `drop_m` (the same
/// `LOS - actual` sign convention derived from `points` above — see
/// `trajectory_sampling.rs`'s `sample_trajectory` doc comment), so no position
/// reconstruction is needed: `sample.drop_m` is simply multiplied by
/// `table.factor_at(mach)`, with `mach` computed from `sample.velocity_mps` via the
/// identical frozen `station_speed_of_sound_mps` divisor used for `points` above. This
/// mirrors `run_sampled_trajectory`'s (come-ups' own sampled-trajectory path, `main.rs`)
/// hand-rolled version of the same transform, so both paths now agree. `None` stays
/// `None` — nothing to scale.
///
/// Nothing else is touched: `position.x` (downrange), `position.z` (windage/lateral),
/// `velocity_magnitude`, `kinetic_energy`, and `time` are byte-identical to their
/// pre-call values on `points`; `distance_m`, `wind_drift_m`, `velocity_mps`,
/// `energy_j`, `time_s`, and `flags` are byte-identical on `sampled_points`; every
/// top-level scalar on `result` itself (`time_of_flight`, `impact_velocity`,
/// `impact_energy`, `max_range`, `max_height`, ...) is untouched too.
pub fn apply_dsf(result: &mut TrajectoryResult, table: &DsfTable) {
    let line_of_sight_height_m = result.line_of_sight_height_m;
    let station_speed_of_sound_mps = result.station_speed_of_sound_mps;

    for point in result.points.iter_mut() {
        let mach = if station_speed_of_sound_mps > 0.0 {
            point.velocity_magnitude / station_speed_of_sound_mps
        } else {
            0.0
        };
        let factor = table.factor_at(mach);
        let drop = line_of_sight_height_m - point.position.y;
        point.position.y = line_of_sight_height_m - drop * factor;
    }

    if let Some(samples) = result.sampled_points.as_mut() {
        for sample in samples.iter_mut() {
            let mach = if station_speed_of_sound_mps > 0.0 {
                sample.velocity_mps / station_speed_of_sound_mps
            } else {
                0.0
            };
            sample.drop_m *= table.factor_at(mach);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli_api::{TrajectoryPoint};
    use crate::trajectory_observation::TrajectoryTermination;
    use crate::trajectory_sampling::{TrajectoryFlag, TrajectorySample};
    use nalgebra::Vector3;

    fn pt(mach: f64, dsf: f64) -> DsfPoint {
        DsfPoint { mach, dsf }
    }

    // ---- factor_at semantics ----

    #[test]
    fn factor_at_identity_at_and_above_ceiling() {
        let table = DsfTable::from_points(vec![pt(0.9, 1.2)]).unwrap();
        assert_eq!(table.factor_at(1.2), 1.0);
        assert_eq!(table.factor_at(1.5), 1.0);
        assert_eq!(table.factor_at(3.0), 1.0);
    }

    #[test]
    fn factor_at_empty_table_is_always_identity() {
        let table = DsfTable::from_points(vec![]).unwrap();
        assert_eq!(table.factor_at(0.5), 1.0);
        assert_eq!(table.factor_at(1.0), 1.0);
        assert_eq!(table.factor_at(1.2), 1.0);
    }

    #[test]
    fn factor_at_single_point_interpolates_to_the_implicit_anchor() {
        // (0.9, 1.15) -> anchor (1.2, 1.0). Halfway (mach 1.05) is halfway between the two
        // DSF values.
        let table = DsfTable::from_points(vec![pt(0.9, 1.15)]).unwrap();
        let expected_half = 1.15 + (1.0 - 1.15) * 0.5;
        assert!((table.factor_at(1.05) - expected_half).abs() < 1e-12);
        // At the point itself: its own dsf.
        assert_eq!(table.factor_at(0.9), 1.15);
        // Continuity at the ceiling boundary: interpolating right up to 1.2 approaches 1.0.
        let near_ceiling = table.factor_at(1.2 - 1e-9);
        assert!((near_ceiling - 1.0).abs() < 1e-6);
    }

    #[test]
    fn factor_at_linear_between_two_keys() {
        // (0.8, 1.2) and (1.0, 1.05); at mach 0.9 (halfway) expect halfway between the DSFs.
        let table = DsfTable::from_points(vec![pt(0.8, 1.2), pt(1.0, 1.05)]).unwrap();
        let expected = 1.2 + (1.05 - 1.2) * 0.5;
        assert!((table.factor_at(0.9) - expected).abs() < 1e-12);
        // Exactly at a key: that key's own dsf.
        assert_eq!(table.factor_at(0.8), 1.2);
        assert_eq!(table.factor_at(1.0), 1.05);
    }

    #[test]
    fn factor_at_flat_clamp_below_lowest() {
        let table = DsfTable::from_points(vec![pt(0.8, 1.2), pt(1.0, 1.05)]).unwrap();
        assert_eq!(table.factor_at(0.5), 1.2);
        assert_eq!(table.factor_at(0.0001), 1.2);
    }

    #[test]
    fn factor_at_interpolates_between_highest_key_and_anchor() {
        // Highest key (1.0, 1.05) -> anchor (1.2, 1.0). At mach 1.1 (halfway).
        let table = DsfTable::from_points(vec![pt(0.8, 1.2), pt(1.0, 1.05)]).unwrap();
        let expected = 1.05 + (1.0 - 1.05) * 0.5;
        assert!((table.factor_at(1.1) - expected).abs() < 1e-12);
    }

    // ---- validation rejections ----

    #[test]
    fn from_points_rejects_mach_at_or_above_ceiling() {
        assert!(DsfTable::from_points(vec![pt(1.2, 1.1)]).is_err());
        assert!(DsfTable::from_points(vec![pt(1.3, 1.1)]).is_err());
    }

    #[test]
    fn from_points_rejects_non_positive_mach() {
        assert!(DsfTable::from_points(vec![pt(0.0, 1.1)]).is_err());
        assert!(DsfTable::from_points(vec![pt(-0.5, 1.1)]).is_err());
    }

    #[test]
    fn from_points_rejects_dsf_out_of_range() {
        assert!(DsfTable::from_points(vec![pt(0.9, 0.0)]).is_err());
        assert!(DsfTable::from_points(vec![pt(0.9, -1.0)]).is_err());
        assert!(DsfTable::from_points(vec![pt(0.9, 0.5)]).is_err()); // exclusive bound
        assert!(DsfTable::from_points(vec![pt(0.9, 2.0)]).is_err()); // exclusive bound
        assert!(DsfTable::from_points(vec![pt(0.9, 2.5)]).is_err());
        assert!(DsfTable::from_points(vec![pt(0.9, f64::NAN)]).is_err());
    }

    #[test]
    fn from_points_rejects_more_than_six_points() {
        let points: Vec<DsfPoint> = (0..7).map(|i| pt(0.1 + i as f64 * 0.1, 1.1)).collect();
        let err = DsfTable::from_points(points).unwrap_err();
        assert!(
            err.contains('6'),
            "error should name the 6-point cap: {err}"
        );
    }

    #[test]
    fn from_points_sorts_ascending_by_mach() {
        let table = DsfTable::from_points(vec![pt(0.9, 1.1), pt(0.3, 1.3), pt(0.6, 1.2)]).unwrap();
        let machs: Vec<f64> = table.points().iter().map(|p| p.mach).collect();
        assert_eq!(machs, vec![0.3, 0.6, 0.9]);
    }

    // ---- upsert ----

    #[test]
    fn upsert_appends_when_no_existing_point_is_within_tolerance() {
        let mut table = DsfTable::from_points(vec![pt(0.5, 1.1)]).unwrap();
        let outcome = table.upsert(pt(0.8, 1.2)).unwrap();
        assert_eq!(outcome, UpsertOutcome::Appended);
        assert_eq!(table.points().len(), 2);
    }

    #[test]
    fn upsert_replaces_within_tolerance() {
        let mut table = DsfTable::from_points(vec![pt(0.5, 1.1)]).unwrap();
        let new_point = pt(0.53, 1.25); // within 0.05 of 0.5
        let outcome = table.upsert(new_point).unwrap();
        match outcome {
            UpsertOutcome::Replaced { old } => assert_eq!(old, pt(0.5, 1.1)),
            other => panic!("expected Replaced, got {other:?}"),
        }
        assert_eq!(table.points().len(), 1);
        assert_eq!(table.points()[0], new_point);
    }

    #[test]
    fn upsert_boundary_just_outside_tolerance_appends() {
        let mut table = DsfTable::from_points(vec![pt(0.5, 1.1)]).unwrap();
        let outcome = table.upsert(pt(0.551, 1.2)).unwrap(); // 0.051 away: outside tolerance
        assert_eq!(outcome, UpsertOutcome::Appended);
        assert_eq!(table.points().len(), 2);
    }

    #[test]
    fn upsert_errors_at_seventh_distinct_point_naming_the_cap() {
        let mut table = DsfTable::from_points(
            (0..6).map(|i| pt(0.1 + i as f64 * 0.15, 1.1)).collect(),
        )
        .unwrap();
        assert_eq!(table.points().len(), 6);
        // Far from every existing point (nearest is 0.85, 0.15 away — outside the 0.05 tolerance).
        let err = table.upsert(pt(1.0, 1.3)).unwrap_err();
        assert!(
            err.contains('6'),
            "error should name the 6-point cap: {err}"
        );
        assert_eq!(table.points().len(), 6, "rejected point must not be added");
    }

    #[test]
    fn upsert_rejects_invalid_point_without_mutating_table() {
        let mut table = DsfTable::from_points(vec![pt(0.5, 1.1)]).unwrap();
        assert!(table.upsert(pt(1.2, 1.1)).is_err());
        assert!(table.upsert(pt(0.6, 3.0)).is_err());
        assert_eq!(table.points().len(), 1, "invalid upsert must not mutate the table");
    }

    // ---- apply_dsf drop-only invariant ----

    fn trajectory_point(time: f64, x: f64, y: f64, z: f64, velocity_magnitude: f64) -> TrajectoryPoint {
        TrajectoryPoint {
            time,
            position: Vector3::new(x, y, z),
            velocity_magnitude,
            kinetic_energy: 0.5 * 0.01 * velocity_magnitude * velocity_magnitude,
            drag_coefficient: None,
        }
    }

    fn trajectory_sample(
        distance_m: f64,
        drop_m: f64,
        wind_drift_m: f64,
        velocity_mps: f64,
        time_s: f64,
        flags: Vec<TrajectoryFlag>,
    ) -> TrajectorySample {
        TrajectorySample {
            distance_m,
            drop_m,
            wind_drift_m,
            velocity_mps,
            energy_j: 0.5 * 0.01 * velocity_mps * velocity_mps,
            time_s,
            flags,
        }
    }

    fn fixture_result(points: Vec<TrajectoryPoint>) -> TrajectoryResult {
        TrajectoryResult {
            max_range: 500.0,
            max_height: 2.0,
            time_of_flight: 1.234,
            impact_velocity: 300.0,
            impact_energy: 1800.0,
            projectile_mass_kg: 0.01,
            line_of_sight_height_m: 0.05,
            station_speed_of_sound_mps: 340.0,
            termination: TrajectoryTermination::MaxRange,
            points,
            sampled_points: None,
            min_pitch_damping: None,
            transonic_mach: None,
            angular_state: None,
            max_yaw_angle: None,
            max_precession_angle: None,
            aerodynamic_jump: None,
            mach_1_2_distance_m: None,
            mach_1_0_distance_m: None,
            mach_0_9_distance_m: None,
        }
    }

    #[test]
    fn apply_dsf_scales_only_drop_leaving_everything_else_byte_identical() {
        let sos = 340.0;
        let points = vec![
            // mach 1.3: supersonic, above the ceiling -> factor 1.0 (untouched drop too).
            trajectory_point(0.0, 0.0, 0.05, 0.0, 1.3 * sos),
            // mach 0.9: within the table -> interpolated factor.
            trajectory_point(0.5, 250.0, 0.02, 1.0, 0.9 * sos),
            // mach 0.5: below the lowest key -> flat-clamped factor.
            trajectory_point(1.0, 500.0, -1.0, 2.0, 0.5 * sos),
        ];
        let mut original = fixture_result(points);
        // MBA-1357 Task 2 review, Critical #2: sampled_points is a SEPARATE array from
        // `points` and must be scaled too (same Mach coverage: 1.3/0.9/0.5, so the same
        // expected_factors below apply to both).
        original.sampled_points = Some(vec![
            // drop_m values chosen to match the drop_before values the points loop below
            // derives (los - y): 0.0, 0.03, 1.05 — so the same expected_factors apply.
            trajectory_sample(0.0, 0.0, 0.0, 1.3 * sos, 0.0, vec![]),
            trajectory_sample(250.0, 0.03, 1.0, 0.9 * sos, 0.5, vec![TrajectoryFlag::MachTransition]),
            trajectory_sample(500.0, 1.05, 2.0, 0.5 * sos, 1.0, vec![TrajectoryFlag::Apex]),
        ]);
        let table = DsfTable::from_points(vec![pt(0.8, 1.2), pt(1.0, 1.05)]).unwrap();

        let mut scaled = original.clone();
        apply_dsf(&mut scaled, &table);

        for (orig, new) in original.points.iter().zip(scaled.points.iter()) {
            assert_eq!(orig.time, new.time, "time must be byte-identical");
            assert_eq!(
                orig.velocity_magnitude, new.velocity_magnitude,
                "velocity must be byte-identical"
            );
            assert_eq!(
                orig.kinetic_energy, new.kinetic_energy,
                "energy must be byte-identical"
            );
            assert_eq!(orig.position.x, new.position.x, "downrange must be byte-identical");
            assert_eq!(orig.position.z, new.position.z, "windage must be byte-identical");
        }
        // Top-level fields are untouched too — every one of them.
        assert_eq!(original.max_range, scaled.max_range);
        assert_eq!(original.max_height, scaled.max_height);
        assert_eq!(original.time_of_flight, scaled.time_of_flight);
        assert_eq!(original.impact_velocity, scaled.impact_velocity);
        assert_eq!(original.impact_energy, scaled.impact_energy);
        assert_eq!(original.projectile_mass_kg, scaled.projectile_mass_kg);
        assert_eq!(original.line_of_sight_height_m, scaled.line_of_sight_height_m);
        assert_eq!(
            original.station_speed_of_sound_mps,
            scaled.station_speed_of_sound_mps
        );
        assert_eq!(original.termination, scaled.termination);
        assert_eq!(original.min_pitch_damping, scaled.min_pitch_damping);
        assert_eq!(original.transonic_mach, scaled.transonic_mach);
        assert_eq!(original.max_yaw_angle, scaled.max_yaw_angle);
        assert_eq!(original.max_precession_angle, scaled.max_precession_angle);
        // AerodynamicJumpComponents has no PartialEq; the fixture carries None and
        // apply_dsf must leave it that way.
        assert!(original.aerodynamic_jump.is_none() && scaled.aerodynamic_jump.is_none());

        let los = original.line_of_sight_height_m;
        let mach_09_factor = 1.2 + (1.05 - 1.2) * 0.5; // mach 0.9: halfway between the two keys
        let expected_factors = [1.0, mach_09_factor, 1.2 /* flat clamp below lowest key */];
        for (i, (orig, new)) in original.points.iter().zip(scaled.points.iter()).enumerate() {
            let drop_before = los - orig.position.y;
            let drop_after = los - new.position.y;
            let expected_drop = drop_before * expected_factors[i];
            assert!(
                (drop_after - expected_drop).abs() < 1e-9,
                "point {i}: expected scaled drop {expected_drop}, got {drop_after}"
            );
        }
        // The untouched (mach >= 1.2) point's position.y must be exactly unchanged.
        assert_eq!(original.points[0].position.y, scaled.points[0].position.y);

        // Critical #2: sampled_points scales the SAME way, and every other field on
        // each sample is byte-identical.
        let orig_samples = original.sampled_points.as_ref().unwrap();
        let scaled_samples = scaled.sampled_points.as_ref().unwrap();
        assert_eq!(orig_samples.len(), scaled_samples.len());
        for (i, (orig, new)) in orig_samples.iter().zip(scaled_samples.iter()).enumerate() {
            assert_eq!(orig.distance_m, new.distance_m, "sample {i}: distance_m must be byte-identical");
            assert_eq!(
                orig.wind_drift_m, new.wind_drift_m,
                "sample {i}: wind_drift_m must be byte-identical"
            );
            assert_eq!(
                orig.velocity_mps, new.velocity_mps,
                "sample {i}: velocity_mps must be byte-identical"
            );
            assert_eq!(orig.energy_j, new.energy_j, "sample {i}: energy_j must be byte-identical");
            assert_eq!(orig.time_s, new.time_s, "sample {i}: time_s must be byte-identical");
            assert_eq!(orig.flags, new.flags, "sample {i}: flags must be byte-identical");

            let expected_drop = orig.drop_m * expected_factors[i];
            assert!(
                (new.drop_m - expected_drop).abs() < 1e-9,
                "sample {i}: expected scaled drop_m {expected_drop}, got {}",
                new.drop_m
            );
        }
        // The untouched (mach >= 1.2) sample's drop_m must be exactly unchanged.
        assert_eq!(orig_samples[0].drop_m, scaled_samples[0].drop_m);
    }

    #[test]
    fn apply_dsf_leaves_sampled_points_none_when_absent() {
        // Same non-empty table as the invariant test above, but sampled_points is None
        // (e.g. a solve without --sample-trajectory) — apply_dsf must not panic or
        // conjure a Some, it must stay None.
        let points = vec![trajectory_point(0.5, 250.0, 0.02, 1.0, 0.9 * 340.0)];
        let original = fixture_result(points);
        assert!(original.sampled_points.is_none());
        let table = DsfTable::from_points(vec![pt(0.8, 1.2), pt(1.0, 1.05)]).unwrap();

        let mut scaled = original.clone();
        apply_dsf(&mut scaled, &table);

        assert!(scaled.sampled_points.is_none(), "None must stay None");
    }

    #[test]
    fn apply_dsf_with_empty_table_leaves_drop_unchanged() {
        let points = vec![trajectory_point(0.5, 250.0, 0.02, 1.0, 0.9 * 340.0)];
        let original = fixture_result(points);
        let table = DsfTable::from_points(vec![]).unwrap();

        let mut scaled = original.clone();
        apply_dsf(&mut scaled, &table);

        assert_eq!(original.points[0].position.y, scaled.points[0].position.y);
    }
}
