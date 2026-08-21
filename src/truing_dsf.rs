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

use crate::cli_api::{
    calculate_zero_angle_with_conditions, AtmosphericConditions, BallisticInputs,
    BcReferenceStandard, DropsReference, TrajectoryPoint, TrajectoryResult, TrajectorySolver,
    WindConditions,
};
use crate::truing::{fallback_bullet_length_m, DragModelArg, TruingModelInputsV1};
use crate::DragModel;

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

/// Linearly interpolate `(position.y, velocity_magnitude)` at horizontal distance
/// `target_dist_m` from a solved trajectory's points (`position.x` = downrange). Mirrors
/// `cli_api::fit_value_at`'s interpolation (private to that module), but resolves both
/// quantities from the same bracketing pair in one pass since the `dsf` verb needs drop
/// AND Mach at the identical range. `None` if the trajectory never reaches `target_dist_m`.
pub fn interpolate_position_and_velocity(
    points: &[TrajectoryPoint],
    target_dist_m: f64,
) -> Option<(f64, f64)> {
    for i in 0..points.len() {
        if points[i].position.x >= target_dist_m {
            if i == 0 {
                return Some((points[0].position.y, points[0].velocity_magnitude));
            }
            let p1 = &points[i - 1];
            let p2 = &points[i];
            let dx = p2.position.x - p1.position.x;
            if dx.abs() < 1e-9 {
                return Some((p2.position.y, p2.velocity_magnitude));
            }
            let t = (target_dist_m - p1.position.x) / dx;
            let y = p1.position.y + t * (p2.position.y - p1.position.y);
            let v = p1.velocity_magnitude + t * (p2.velocity_magnitude - p1.velocity_magnitude);
            return Some((y, v));
        }
    }
    None
}

/// Whether an observation range is beyond 90% of the trajectory's solved max range —
/// past this point the solution's reliability degrades (short-range extrapolation of a
/// trajectory that terminated, e.g., at ground impact just past the observation).
///
/// Moved out of `main.rs` alongside [`dsf_observation_warrants_90pct_warning`] (MBA-1357
/// Task 8): it isn't one of that task's four named helpers, but
/// `dsf_observation_warrants_90pct_warning` calls it, and a library function cannot call a
/// private binary-crate function, so it had to move too. Verbatim, unchanged.
pub fn dsf_observation_beyond_90pct(range_m: f64, solved_max_range_m: f64) -> bool {
    solved_max_range_m > 0.0 && range_m > 0.9 * solved_max_range_m
}

/// The downrange distance (meters) where the trajectory's station Mach first drops
/// below 1.0 (the "crossed_subsonic" transition), linearly interpolated between the
/// bracketing solved points.
///
/// Mirrors `MachTransitionTracker::record_downward_crossings`'s `crossed_subsonic` event
/// (`cli_api.rs` ~1744), reimplemented here because that tracker is a private
/// `cli_api.rs` type, and the crossing distances it collects (`transonic_distances`)
/// aren't retained on `TrajectoryResult` itself — they're only consumed to flag
/// `TrajectorySample`s (`trajectory_sampling::add_trajectory_flags`), which requires
/// trajectory sampling to be enabled. `solve_profile_for_dsf` solves with sampling off
/// (it only needs `result.points`), so this recomputes the same crossing from the plain
/// points array using the identical Mach divisor `apply_dsf` and the observation path
/// use (`velocity_magnitude / station_speed_of_sound_mps`).
///
/// Returns `None` if the trajectory never goes subsonic within the solved points (still
/// supersonic/transonic at the last point, an empty/degenerate solve, or a non-finite
/// station speed of sound).
pub fn mach_1_crossing_range_m(result: &TrajectoryResult) -> Option<f64> {
    let sos = result.station_speed_of_sound_mps;
    if sos <= 0.0 || !sos.is_finite() {
        return None;
    }

    let mut previous: Option<(f64, f64)> = None; // (downrange_m, mach)
    for point in &result.points {
        let mach = point.velocity_magnitude / sos;
        if let Some((prev_x, prev_mach)) = previous {
            if prev_mach >= 1.0 && mach < 1.0 {
                let denom = prev_mach - mach;
                if denom.abs() < f64::EPSILON {
                    return Some(point.position.x);
                }
                let t = (prev_mach - 1.0) / denom;
                return Some(prev_x + t * (point.position.x - prev_x));
            }
        }
        previous = Some((point.position.x, mach));
    }
    None
}

/// Whether the `dsf` verb's "solution reliability degrades" warning should fire.
///
/// Two independent gates, BOTH required (MBA-1357 Task 2 review, Critical #1): a solve
/// envelope sized to comfortably exceed a typical observation made the 90%-of-solved-
/// range check alone fire unconditionally, since the envelope tracked the observation
/// itself rather than the profile's own real reach. Requiring the observation to also be
/// beyond the trajectory's Mach-1.0 crossing ties the warning to an actual downrange-
/// position judgment (deep into the subsonic regime the DSF table's low end targets),
/// not merely wherever the caller happened to stop solving.
pub fn dsf_observation_warrants_90pct_warning(
    range_m: f64,
    mach_1_crossing_range_m: Option<f64>,
    solved_max_range_m: f64,
) -> bool {
    let beyond_mach_1_crossing = mach_1_crossing_range_m
        .map(|crossing_m| range_m > crossing_m)
        .unwrap_or(false);
    beyond_mach_1_crossing && dsf_observation_beyond_90pct(range_m, solved_max_range_m)
}

/// Full input set for [`solve_for_dsf`] — the scalar-BC model (mirroring
/// [`TruingModelInputsV1`]'s fields) plus every profile field the CLI's historical
/// `solve_profile_for_dsf` fed into the physics that `TruingModelInputsV1` alone has no
/// slot for (MBA-1357 Task 8 review, Finding 1). `None` on any `Option` field means
/// exactly what it meant to the historical code when a profile didn't carry that field —
/// the same physically neutral default, documented per field below — so a profile that
/// sets none of them solves byte-identically to a bare converted `TruingModelInputsV1`,
/// and one that does gets ALL of it honored, not silently dropped.
///
/// Three fields the old profile-driven solve also read are intentionally NOT carried
/// here, confirmed empirically rather than just by reading the physics: a stored
/// twist-rate override, twist direction, and a stored bullet-length override. Every place
/// any of the three reaches the trajectory (`enable_magnus`'s Magnus force,
/// `enable_aerodynamic_jump`'s Litz jump estimator, `enable_precession_nutation`'s
/// spin-rate term) is gated behind one of those flags, all permanently off in this solve
/// (see [`solve_for_dsf`]); the only other consumer is the CLI's own separately-computed
/// "Stability (SG)" display line, which `dsf` never prints. Verified by diffing
/// `trajectory --saved-profile` output between two otherwise-identical saved profiles
/// differing only in `--twist-rate` (7 vs 20) and, separately, only in `--bullet-length`
/// (1.0in vs 2.5in): `Max Range`, `Max Height`, `Zero Angle`, `Time of Flight`, `Impact
/// Velocity`, `Impact Energy`, and the ground-impact range were byte-identical in both
/// pairs; only the SG display line (not part of `TrajectoryResult`) differed.
#[derive(Debug, Clone)]
pub struct DsfSolveInputs {
    /// Muzzle velocity, feet/second.
    pub muzzle_velocity_fps: f64,
    /// Nominal scalar ballistic coefficient — superseded by `bc_segments`/
    /// `custom_drag_table` when either is `Some`, same precedence the solver applies.
    pub ballistic_coefficient: f64,
    pub drag_model: DragModelArg,
    /// Bullet mass, grains.
    pub mass_gr: f64,
    /// Bullet diameter, inches.
    pub diameter_in: f64,
    /// Sight height over bore, inches.
    pub sight_height_in: f64,
    /// Ambient temperature, degrees Fahrenheit.
    pub temperature_f: f64,
    /// Station pressure, inches of mercury.
    pub pressure_inhg: f64,
    /// Relative humidity, percent (0 through 100).
    pub humidity_pct: f64,
    /// Altitude, feet.
    pub altitude_ft: f64,

    /// Wind speed, meters/second. `None` (the default) is calm — byte-identical to a bare
    /// `TruingModelInputsV1` solve.
    pub wind_speed_mps: Option<f64>,
    /// Wind direction, radians, wind-FROM convention (0 = headwind). `None` is 0.0.
    pub wind_direction_rad: Option<f64>,
    /// Uphill (positive) / downhill (negative) shooting angle, radians. `None` is level
    /// (0.0). Materially changes predicted drop when set — this is the field most likely
    /// to matter of everything in this struct.
    pub shooting_angle_rad: Option<f64>,
    /// Deliberate vertical point-of-impact offset AT THE ZERO RANGE, meters (MBA-1359
    /// semantics — see [`crate::cli_api::BallisticInputs::zero_poi_vertical_m`]). `None`
    /// is 0.0 (no bias). Directly shifts predicted drop, unlike the two lateral-only
    /// fields below.
    pub zero_poi_vertical_m: Option<f64>,
    /// Deliberate horizontal point-of-impact offset at the zero range, meters. `None` is
    /// 0.0. Carried for full-fidelity `TrajectoryResult` output (lateral position); inert
    /// for the vertical drop/Mach the `dsf` command itself reads off the result.
    pub zero_poi_horizontal_m: Option<f64>,
    /// Lateral sight-to-bore mount offset, meters (MBA-1396). `None` is 0.0. Same
    /// "carried for fidelity, lateral-only" note as `zero_poi_horizontal_m`.
    pub sight_offset_lateral_m: Option<f64>,
    /// Which standard atmosphere `ballistic_coefficient`/`bc_segments` are referenced to.
    /// `None` is ICAO (matches every profile saved before this field existed). A real,
    /// non-cosmetic ~1.8% retardation difference when Army Standard Metro.
    pub bc_reference_standard: Option<BcReferenceStandard>,
    /// Velocity-banded BC schedule (MBA-1323 Phase 2). `None` uses the scalar
    /// `ballistic_coefficient` alone. When `Some` and non-empty, this REPLACES the scalar
    /// BC for the solve — the solver's own segments-then-scalar precedence, same as the
    /// profile path always used. No separate `use_bc_segments` bool: an empty/absent
    /// schedule is a no-op regardless of such a flag (the solver gates on
    /// `use_bc_segments && !segments.is_empty()`), and a populated schedule is honored
    /// unconditionally once present, so the historical `profile.use_bc_segments` flag
    /// could never actually suppress a populated one — carrying it here would be
    /// redundant, not lossy.
    pub bc_segments: Option<Vec<crate::BCSegmentData>>,
    /// Full Mach/Cd drag curve (MBA-1323 Phase 2, `.a7p` CUSTOM import). `None` uses
    /// `ballistic_coefficient`/`bc_segments` instead. When `Some`, replaces the BC model
    /// entirely, same as the profile path.
    pub custom_drag_table: Option<crate::drag::DragTable>,
    /// Resolved zero distance, YARDS (the same imperial convention every other field here
    /// that has one uses). `None` means no zero is configured — the solve stays flat
    /// (`muzzle_angle` 0.0), byte-identical to the CLI's historical behaviour for a
    /// profile with neither `auto_zero` nor `zero_distance` set and no `--zero-set`
    /// selected (MBA-1357 Task 8 review, Finding 2).
    pub zero_distance_yd: Option<f64>,
}

impl From<&TruingModelInputsV1> for DsfSolveInputs {
    /// Widen a bare scalar-BC model into the full [`DsfSolveInputs`] shape a caller with
    /// no profile-shaped extras reaches (Task 9's `true.dsf` bridge command). EXPLICIT
    /// bridge defaults — spelled out here because an app author needs to know what this
    /// assumes: **no wind, a level shot (0.0 shooting angle), no zero-POI bias, no
    /// sight-mount offset, ICAO BC reference, no BC-segment schedule, no custom drag
    /// curve.** `zero_distance_yd` carries `inputs.zero_distance_yd` (mandatory on
    /// `TruingModelInputsV1`, so always `Some` here — never the "flat, no zero" case).
    fn from(inputs: &TruingModelInputsV1) -> Self {
        DsfSolveInputs {
            muzzle_velocity_fps: inputs.muzzle_velocity_fps,
            ballistic_coefficient: inputs.ballistic_coefficient,
            drag_model: inputs.drag_model,
            mass_gr: inputs.mass_gr,
            diameter_in: inputs.diameter_in,
            sight_height_in: inputs.sight_height_in,
            temperature_f: inputs.temperature_f,
            pressure_inhg: inputs.pressure_inhg,
            humidity_pct: inputs.humidity_pct,
            altitude_ft: inputs.altitude_ft,
            wind_speed_mps: None,
            wind_direction_rad: None,
            shooting_angle_rad: None,
            zero_poi_vertical_m: None,
            zero_poi_horizontal_m: None,
            sight_offset_lateral_m: None,
            bc_reference_standard: None,
            bc_segments: None,
            custom_drag_table: None,
            zero_distance_yd: Some(inputs.zero_distance_yd),
        }
    }
}

/// Solve a [`DsfSolveInputs`]'s own trajectory for the `dsf` command's derivation step
/// (MBA-1357 Task 8), given plain values directly rather than a saved `Profile` — the
/// JSON bridge cannot construct a `Profile`, and must not read one from disk.
///
/// This is the library half of what was `main.rs`'s private `solve_profile_for_dsf`. The
/// CLI keeps its saved-profile path by converting a loaded profile into `DsfSolveInputs`
/// (ALL of the profile fields the historical solve honored — see that struct's own doc
/// comment for exactly which, and which three are deliberately absent because they're
/// provably inert here) and delegating. A bridge caller with only a `TruingModelInputsV1`
/// can go through `DsfSolveInputs::from` instead — see that impl's doc comment for the
/// explicit defaults it assumes.
///
/// `max_range_m` is the solve envelope; `inputs.zero_distance_yd` (`None` = flat, no zero
/// applied — see that field's doc comment) supplies the zero.
///
/// Advanced physics toggles this struct has no field for (Magnus/Coriolis/spin-drift/
/// aerodynamic-jump/wind-shear/pitch-damping/precession/powder-sensitivity/`cd_scale`,
/// plus twist rate/direction and bullet length — see `DsfSolveInputs`'s doc comment) stay
/// off/neutral, matching `solve_profile_for_dsf`'s own historical behaviour for every
/// profile-only command (`come-ups`, `lead`, `mpbr`).
pub fn solve_for_dsf(
    inputs: &DsfSolveInputs,
    max_range_m: f64,
) -> Result<TrajectoryResult, String> {
    let velocity_m = inputs.muzzle_velocity_fps * 0.3048;
    let mass_kg = inputs.mass_gr * crate::constants::GRAINS_TO_KG;
    let diameter_m = inputs.diameter_in * 0.0254;
    let sight_height_m = inputs.sight_height_in * 0.0254;
    let bullet_length_m = fallback_bullet_length_m(diameter_m, mass_kg);

    // Saved profiles predate the CLI's unified --bore-height flag and never stored one;
    // solve_profile_for_dsf always fell back to this same constant (60 in) regardless of
    // profile contents, so reusing it here is not a narrowing — it was already
    // unconditional.
    let bore_height_m = 60.0 * 0.0254;

    let temperature_c = (inputs.temperature_f - 32.0) * 5.0 / 9.0;
    let pressure_hpa = inputs.pressure_inhg * 33.8639;
    let altitude_m = inputs.altitude_ft * 0.3048;

    let drag_model = match inputs.drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let wind_speed_m = inputs.wind_speed_mps.unwrap_or(0.0);
    let wind_direction_rad = inputs.wind_direction_rad.unwrap_or(0.0);
    let shooting_angle_rad = inputs.shooting_angle_rad.unwrap_or(0.0);
    let zero_poi_vertical_m = inputs.zero_poi_vertical_m.unwrap_or(0.0);
    let zero_poi_horizontal_m = inputs.zero_poi_horizontal_m.unwrap_or(0.0);
    let sight_offset_lateral_m = inputs.sight_offset_lateral_m.unwrap_or(0.0);
    let bc_reference_standard = inputs
        .bc_reference_standard
        .unwrap_or(BcReferenceStandard::Icao);
    let bc_segments_data = inputs.bc_segments.clone();
    let use_bc_segments = bc_segments_data.is_some();

    let wind = WindConditions {
        speed: wind_speed_m,
        direction: wind_direction_rad,
        vertical_speed: 0.0,
    };
    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        // NOTE: matches solve_profile_for_dsf's (and run_trajectory's) own convention —
        // the same raw (0-100) value feeds both AtmosphericConditions.humidity (percent,
        // what the solve actually reads) and BallisticInputs.humidity below (nominally a
        // 0-1 fraction per its own doc comment). Replicated rather than "fixed" — see
        // that helper's own note; BallisticInputs.humidity is otherwise inert for a plain
        // (non-Monte-Carlo) solve.
        humidity: inputs.humidity_pct,
        altitude: altitude_m,
    };

    let mut ballistic_inputs = BallisticInputs {
        bc_value: inputs.ballistic_coefficient,
        bc_type: drag_model,
        bc_reference_standard,
        bullet_mass: mass_kg,
        muzzle_velocity: velocity_m,
        bullet_diameter: diameter_m,
        bullet_length: bullet_length_m,

        muzzle_angle: 0.0,
        target_distance: max_range_m,
        azimuth_angle: 0.0,
        shot_azimuth: 0.0,
        shooting_angle: shooting_angle_rad,
        cant_angle: 0.0,
        sight_height: sight_height_m,
        sight_offset_lateral_m,
        muzzle_height: bore_height_m,
        target_height: 0.0,
        zero_poi_vertical_m,
        zero_poi_horizontal_m,
        // Ground-impact detection ON (0.0), matching solve_profile_for_dsf's override of
        // the engine's own default (-100.0, effectively disabled) — the `dsf` command has
        // no way to disable it.
        ground_threshold: 0.0,

        altitude: altitude_m,
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity: inputs.humidity_pct,
        latitude: None,

        wind_speed: wind_speed_m,
        wind_angle: wind_direction_rad,

        // Confirmed inert for this solve's own drop/Mach output (see DsfSolveInputs's doc
        // comment) — kept at the same historical default this helper always fell back to
        // when a profile carried no override (twist_right stays the historical `true`
        // default too, for the same reason).
        twist_rate: crate::stability::default_twist_inches(diameter_m, mass_kg, velocity_m),
        is_twist_right: true,
        caliber_inches: diameter_m / 0.0254,
        weight_grains: mass_kg / crate::constants::GRAINS_TO_KG,
        manufacturer: None,
        bullet_model: None,
        bullet_id: None,
        bullet_cluster: None,

        use_rk4: true,
        use_adaptive_rk45: true,

        enable_advanced_effects: false,
        enable_magnus: false,
        enable_coriolis: false,
        use_powder_sensitivity: false,
        powder_temp_sensitivity: 0.0,
        powder_temp: 0.0,
        powder_temp_curve: None,
        powder_curve_temp_c: None,
        tipoff_yaw: 0.0,
        cd_delta2: 7.5,
        tipoff_decay_distance: 50.0,

        use_bc_segments,
        bc_segments: None,
        bc_segments_data,
        use_enhanced_spin_drift: false,
        use_form_factor: false,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        enable_trajectory_sampling: false,
        sample_interval: 0.0,
        drops_reference: DropsReference::Los,
        enable_pitch_damping: false,
        enable_precession_nutation: false,
        enable_aerodynamic_jump: false,
        use_cluster_bc: false,

        custom_drag_table: inputs.custom_drag_table.clone(),
        cd_scale: 1.0,

        bc_type_str: None,
    };

    // Zero angle: flat (0.0) when no zero is configured, same as a bare
    // `trajectory --saved-profile NAME` with no --auto-zero (MBA-1357 Task 8 review,
    // Finding 2 — restores the historical Option<f64> semantics the first version of this
    // function dropped).
    if let Some(zero_distance_yd) = inputs.zero_distance_yd {
        let zero_distance_m = zero_distance_yd * 0.9144;
        ballistic_inputs.muzzle_angle = calculate_zero_angle_with_conditions(
            ballistic_inputs.clone(),
            zero_distance_m,
            bore_height_m + sight_height_m,
            wind.clone(),
            atmosphere.clone(),
        )
        .map_err(|e| e.to_string())?;
        // MBA-1359: the returned angle carries the vertical POI bias; the horizontal bias
        // is an azimuth term the flight inputs must apply themselves.
        ballistic_inputs.azimuth_angle += ballistic_inputs.windage_zero_bias_rad(zero_distance_m);
    }

    let mut solver = TrajectorySolver::new(ballistic_inputs, wind, atmosphere);
    solver.set_max_range(max_range_m);
    solver.set_time_step(0.001);
    solver.solve().map_err(|e| e.to_string())
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
