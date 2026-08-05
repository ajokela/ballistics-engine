//! `HoldCurve`: a solved, sampled drop-vs-range curve read by interpolation, and the
//! sampled-trajectory helpers it is built on (MBA-1361/MBA-1362).
//!
//! Promoted out of the CLI binary (0.33.0 decision-support Task 8) so a library consumer can
//! reuse it without the `cli` feature -- starting with the constant-drop range-card work
//! [`HoldCurve`]'s own doc comment already nominates as its next consumer. No feature gate:
//! must compile for wasm32 (pure math over already-resolved inputs; no fs, no clap -- CLI
//! argument resolution, `InverseSolverLoadArgs::resolve`, stays in `main.rs` and constructs
//! [`HoldCurveLoad`] there).

use crate::constants::GRAINS_TO_KG;
use crate::trajectory_observation::{bracket_param, Bracket};
use crate::truing::fallback_bullet_length_m;
use crate::truing_dsf::DsfTable;
use crate::{
    trajectory_sampling, AtmosphericConditions, BCSegmentData, BallisticInputs, DragModel,
    TrajectorySolver, WindConditions,
};
use std::error::Error;

/// Shared: build BallisticInputs + atmosphere + wind from common parameters (all in metric)
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the shared CLI trajectory compatibility helper"
)]
pub fn build_trajectory_components(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModel,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    wind_speed: f64,
    wind_direction: f64,
    max_range: f64,
    sample_interval: f64,
    // MBA-1323 Phase 2: saved-profile velocity-BC segments / Mach-Cd drag curve, already
    // resolved to engine shapes by the caller (bc_segments_from_profile /
    // drag_table_from_profile). `None` for every caller that does not (yet) source these from
    // a saved profile — see handle_come_ups/handle_lead for the callers that do.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<crate::drag::DragTable>,
    // MBA-1359: deliberate POI offset at the zero range (meters), carried on the inputs so
    // the zero solves that reuse these components inherit it. 0.0 = no offset.
    zero_poi_vertical_m: f64,
    zero_poi_horizontal_m: f64,
    // MBA-1396: lateral sight-to-bore mount offset (meters, positive = sight right of
    // bore). 0.0 = sight directly above the bore.
    sight_offset_lateral_m: f64,
) -> (BallisticInputs, WindConditions, AtmosphericConditions) {
    let drag_model_enum = drag_model;
    let wind_direction_rad = wind_direction.to_radians();
    let use_bc_segments = bc_segments_data.is_some();

    let inputs = BallisticInputs {
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass,
        muzzle_velocity: velocity,
        bullet_diameter: diameter,
        bullet_length: fallback_bullet_length_m(diameter, mass),
        muzzle_angle: 0.0,
        target_distance: max_range,
        sight_height,
        sight_offset_lateral_m,
        zero_poi_vertical_m,
        zero_poi_horizontal_m,
        altitude,
        temperature,
        pressure,
        humidity,
        wind_speed,
        wind_angle: wind_direction_rad,
        use_rk4: true,           // Required for non-Euler solver
        use_adaptive_rk45: true, // Use RK45 adaptive (default solver)
        enable_trajectory_sampling: true,
        sample_interval,
        caliber_inches: diameter / 0.0254,
        weight_grains: mass / GRAINS_TO_KG,
        twist_rate: 12.0,
        is_twist_right: true,
        use_bc_segments,
        bc_segments_data,
        custom_drag_table,
        ..Default::default()
    };

    // wind_direction enters from the CLI in degrees; both engine structures use radians.
    let wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction_rad,
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature,
        pressure,
        humidity,
        altitude,
    };

    (inputs, wind, atmosphere)
}

/// Run a trajectory and return sampled points at the given zero angle
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the shared sampled-trajectory compatibility helper"
)]
pub fn run_sampled_trajectory(
    velocity: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    drag_model: DragModel,
    sight_height: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    altitude: f64,
    wind_speed: f64,
    wind_direction: f64,
    max_range: f64,
    sample_interval: f64,
    zero_angle_rad: f64,
    // MBA-1323 Phase 2: see build_trajectory_components's doc comment on these two.
    bc_segments_data: Option<Vec<BCSegmentData>>,
    custom_drag_table: Option<crate::drag::DragTable>,
    // MBA-1357: saved profile's DSF table, already validated by the caller. `None` for
    // every call site except come-ups (the only sampled-trajectory command with a DSF
    // auto-apply story so far).
    dsf_table: Option<&DsfTable>,
    // MBA-1359: deliberate POI offset at the zero range (meters) plus the distance the
    // caller's zero solve used (Some iff a zero was solved). The vertical offset already
    // rides inside `zero_angle_rad` (the zero solve returns it biased); the horizontal
    // offset becomes an azimuth bias here, exactly as calculate_and_set_zero_angle applies
    // it to its own solver. (0.0, 0.0, _) and (_, _, None) are exact no-ops.
    zero_poi_vertical_m: f64,
    zero_poi_horizontal_m: f64,
    // MBA-1396: lateral sight-mount offset (meters). Physically displaces the initial
    // lateral position, and (via windage_zero_bias_rad) joins the azimuth convergence
    // below when a zero was solved.
    sight_offset_lateral_m: f64,
    zero_solve_distance_m: Option<f64>,
) -> Result<Vec<trajectory_sampling::TrajectorySample>, Box<dyn Error>> {
    let (mut inputs, wind, atmosphere) = build_trajectory_components(
        velocity,
        bc,
        mass,
        diameter,
        drag_model,
        sight_height,
        temperature,
        pressure,
        humidity,
        altitude,
        wind_speed,
        wind_direction,
        max_range,
        sample_interval,
        bc_segments_data,
        custom_drag_table,
        zero_poi_vertical_m,
        zero_poi_horizontal_m,
        sight_offset_lateral_m,
    );
    inputs.muzzle_angle = zero_angle_rad;
    if let Some(zero_distance_m) = zero_solve_distance_m {
        inputs.azimuth_angle += inputs.windage_zero_bias_rad(zero_distance_m);
    }

    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(max_range);
    solver.set_time_step(0.001);
    let result = solver.solve()?;

    let mut samples = result.sampled_points.unwrap_or_default();
    // MBA-1357: drop-only correction, same convention as apply_dsf — TrajectorySample's
    // `drop_m` uses the identical `LOS - actual` sign (see sample_trajectory's own doc
    // comment in trajectory_sampling.rs), so scaling it directly by the DSF factor is
    // the exact per-point transform apply_dsf performs on a full TrajectoryResult.
    // Per-point Mach uses the SAME frozen station speed of sound divisor apply_dsf uses
    // (velocity_mps / station_speed_of_sound_mps), not a per-altitude recompute.
    if let Some(table) = dsf_table {
        let station_sos = result.station_speed_of_sound_mps;
        for s in samples.iter_mut() {
            let mach = if station_sos > 0.0 {
                s.velocity_mps / station_sos
            } else {
                0.0
            };
            s.drop_m *= table.factor_at(mach);
        }
    }

    Ok(samples)
}

/// Everything one sampled hold curve needs, already in METRIC (MBA-1361/MBA-1362).
///
/// Flat and small on purpose: it is built once at a CLI boundary and handed to
/// [`HoldCurve::solve`], which owns the zero solve so every consumer zeroes identically.
#[derive(Debug, Clone)]
pub struct HoldCurveLoad {
    pub velocity_mps: f64,
    pub bc: f64,
    pub mass_kg: f64,
    pub diameter_m: f64,
    pub drag_model: DragModel,
    pub sight_height_m: f64,
    pub zero_distance_m: f64,
    pub temperature_c: f64,
    pub pressure_hpa: f64,
    pub humidity: f64,
    pub altitude_m: f64,
    pub wind_speed_mps: f64,
    pub wind_direction_deg: f64,
}

/// One point read off a [`HoldCurve`], in the ANGULAR units a hold is expressed in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoldPoint {
    pub range_m: f64,
    /// Milliradians below the line of sight (positive = the shooter holds over).
    pub drop_mil: f64,
    /// Milliradians of lateral deflection, positive = the bullet goes RIGHT.
    pub wind_mil: f64,
    pub velocity_mps: f64,
    pub energy_j: f64,
    pub time_s: f64,
}

/// A solved, sampled drop-vs-range curve expressed in ANGULAR units (MBA-1361/MBA-1362).
///
/// One solve, sampled finely, then read by interpolation — the same shape `come-ups` uses,
/// but keyed on angle rather than linear drop because that is what a reticle mark is.
///
/// **This is THE drop-vs-range helper**, forward ([`Self::at_range`]) and inverse
/// ([`Self::range_for_angular_drop_mil`]), and all four consumers go through it:
/// `reticle hold --range` (MBA-1361) plus `mark-to-range`, `bdc-match` and `optimal-zero`
/// (MBA-1362). Sharing it is the point — three separately-written root finds would be
/// three chances for "the drop at 500 yards" to mean three different things. The future
/// constant-drop range-card ticket is mark-to-range with uniformly spaced marks and should
/// reuse it too rather than growing a fourth.
pub struct HoldCurve {
    samples: Vec<trajectory_sampling::TrajectorySample>,
}

/// Result of inverting the hold curve for one mark subtension (MBA-1362).
///
/// The two failure arms are REPORTED, never silently dropped: a mark the load cannot reach
/// is information a shooter needs (it is the reticle's usable range limit), and dropping it
/// would quietly shorten the answer table.
#[derive(Debug, Clone, PartialEq)]
pub enum MarkToRangeOutcome {
    /// The mark's subtension is matched at this range.
    Reached(Box<HoldPoint>),
    /// The subtension corresponds to no range past the far zero. Angular drop is exactly
    /// zero AT the far zero and only grows past it, so a mark at or ABOVE the optical
    /// center (a non-positive subtension) is never matched — those are under-holds for
    /// ranges inside the zero, where drop-vs-range is not monotone and "the range for this
    /// mark" is not well defined.
    InsideZero { far_zero_range_m: f64 },
    /// Angular drop never grows to the subtension within the searched trajectory.
    BeyondSearch {
        max_range_m: f64,
        max_drop_mil: f64,
    },
}

impl HoldCurve {
    /// Sample interval used by every hold curve, meters (~1 yard).
    ///
    /// Fine enough that linear interpolation between neighbours is well below the
    /// resolution any reticle can be read to, and coarse enough that a 1500 m curve is a
    /// few thousand points.
    pub const SAMPLE_INTERVAL_M: f64 = 0.9144;

    /// Solve once and sample out to `max_range_m`.
    pub fn solve(load: &HoldCurveLoad, max_range_m: f64) -> Result<Self, Box<dyn Error>> {
        if !max_range_m.is_finite() || max_range_m <= 0.0 {
            return Err("hold curve max range must be finite and greater than zero".into());
        }
        let zero_inputs = BallisticInputs {
            bc_value: load.bc,
            bc_type: load.drag_model,
            bullet_mass: load.mass_kg,
            muzzle_velocity: load.velocity_mps,
            bullet_diameter: load.diameter_m,
            bullet_length: fallback_bullet_length_m(load.diameter_m, load.mass_kg),
            sight_height: load.sight_height_m,
            use_rk4: true,
            ..Default::default()
        };
        let atmosphere = AtmosphericConditions {
            temperature: load.temperature_c,
            pressure: load.pressure_hpa,
            humidity: load.humidity,
            altitude: load.altitude_m,
        };
        let zero_angle = crate::calculate_zero_angle_with_conditions(
            zero_inputs,
            load.zero_distance_m,
            load.sight_height_m,
            WindConditions::default(),
            atmosphere,
        )?;

        let samples = run_sampled_trajectory(
            load.velocity_mps,
            load.bc,
            load.mass_kg,
            load.diameter_m,
            load.drag_model,
            load.sight_height_m,
            load.temperature_c,
            load.pressure_hpa,
            load.humidity,
            load.altitude_m,
            load.wind_speed_mps,
            load.wind_direction_deg,
            max_range_m,
            Self::SAMPLE_INTERVAL_M,
            zero_angle,
            None,
            None,
            None,
            0.0,
            0.0,
            0.0,
            Some(load.zero_distance_m),
        )?;
        if samples.len() < 2 {
            return Err(
                "the trajectory produced too few sampled points to read a hold from".into(),
            );
        }
        Ok(Self { samples })
    }

    /// The furthest range this curve reaches, meters.
    pub fn max_sampled_range_m(&self) -> f64 {
        self.samples.last().map_or(0.0, |s| s.distance_m)
    }

    /// This curve's own sample ranges, in order, meters.
    ///
    /// Exact multiples of [`Self::SAMPLE_INTERVAL_M`] (`i as f64 * SAMPLE_INTERVAL_M` for
    /// `i = 0..N`), the same arithmetic sequence a caller would otherwise have to reproduce by
    /// hand to reason about where this curve was actually verified -- an additive accessor so
    /// no consumer needs read access to the private `samples` field just to answer "which
    /// ranges did this curve solve at."
    pub fn sample_ranges_m(&self) -> Vec<f64> {
        self.samples.iter().map(|s| s.distance_m).collect()
    }

    /// Linearly interpolate the angular hold at `range_m`.
    ///
    /// `None` when the range is outside the sampled span or non-positive (an angular drop
    /// is undefined at the muzzle — it divides by the range).
    pub fn at_range(&self, range_m: f64) -> Option<HoldPoint> {
        if !range_m.is_finite() || range_m <= 0.0 {
            return None;
        }
        let samples = &self.samples;
        let Bracket::Inside { lo, t } =
            bracket_param(samples.len(), |i| samples[i].distance_m, range_m)
        else {
            return None;
        };
        let hi = lo + 1;
        let lerp = |a: f64, b: f64| a + (b - a) * t;
        let drop_m = lerp(samples[lo].drop_m, samples[hi].drop_m);
        let drift_m = lerp(samples[lo].wind_drift_m, samples[hi].wind_drift_m);
        Some(HoldPoint {
            range_m,
            // Milliradian small-angle definition: 1 mil subtends 1/1000 of the range.
            drop_mil: drop_m / range_m * 1000.0,
            wind_mil: drift_m / range_m * 1000.0,
            velocity_mps: lerp(samples[lo].velocity_mps, samples[hi].velocity_mps),
            energy_j: lerp(samples[lo].energy_j, samples[hi].energy_j),
            time_s: lerp(samples[lo].time_s, samples[hi].time_s),
        })
    }

    /// The downrange distance of the FAR zero crossing, meters — the point past which
    /// angular drop grows monotonically with range.
    ///
    /// Angular drop is not monotone over the whole flight: it starts large and positive at
    /// the muzzle (the bullet is a sight height below the line of sight, divided by a tiny
    /// range), falls through zero at the near zero, goes negative while the bullet rides
    /// above the line of sight, and returns through zero at the far zero. Only past that
    /// second crossing is the inverse below single-valued, so the search domain starts
    /// there rather than at the muzzle.
    pub fn far_zero_range_m(&self) -> f64 {
        let mut far = self.samples.first().map_or(0.0, |s| s.distance_m);
        for sample in &self.samples {
            if sample.distance_m > 0.0 && sample.drop_m <= 0.0 {
                far = sample.distance_m;
            }
        }
        far
    }

    /// Bisection cap for the inverse below. The curve is monotone on the searched interval,
    /// so bisection halves the bracket every step; 80 iterations takes any realistic
    /// bracket far below the tolerance and exists only as a runaway guard.
    const INVERSE_MAX_ITERATIONS: u32 = 80;

    /// Bracket width (meters) at which the inverse stops regardless of residual — a
    /// hundredth of a millimeter, four orders below anything a range card resolves.
    const INVERSE_TOLERANCE_M: f64 = 1.0e-5;

    /// Invert the curve: the range at which the angular drop equals `target_mil`.
    ///
    /// Bisection over the interpolated curve on `[far zero, furthest sample]`, where drop
    /// is monotone increasing in range. Both out-of-domain cases come back as their own
    /// outcome rather than as a clamped range.
    pub fn range_for_angular_drop_mil(&self, target_mil: f64) -> MarkToRangeOutcome {
        let far_zero_range_m = self.far_zero_range_m();
        let max_range_m = self.max_sampled_range_m();
        let drop_at = |range_m: f64| self.at_range(range_m).map(|p| p.drop_mil);
        let max_drop_mil = drop_at(max_range_m).unwrap_or(f64::NEG_INFINITY);

        if !target_mil.is_finite() || target_mil <= 0.0 {
            return MarkToRangeOutcome::InsideZero { far_zero_range_m };
        }
        if target_mil > max_drop_mil {
            return MarkToRangeOutcome::BeyondSearch {
                max_range_m,
                max_drop_mil,
            };
        }
        // Start the bracket just past the far zero: exactly at it the drop is 0, and any
        // positive target is therefore bracketed.
        let mut lo = far_zero_range_m.max(Self::SAMPLE_INTERVAL_M);
        if drop_at(lo).is_none_or(|drop| drop >= target_mil) {
            return MarkToRangeOutcome::InsideZero { far_zero_range_m };
        }
        let mut hi = max_range_m;
        for _ in 0..Self::INVERSE_MAX_ITERATIONS {
            if hi - lo <= Self::INVERSE_TOLERANCE_M {
                break;
            }
            let mid = 0.5 * (lo + hi);
            match drop_at(mid) {
                Some(drop) if drop < target_mil => lo = mid,
                Some(_) => hi = mid,
                None => break,
            }
        }
        match self.at_range(0.5 * (lo + hi)) {
            Some(point) => MarkToRangeOutcome::Reached(Box::new(point)),
            None => MarkToRangeOutcome::BeyondSearch {
                max_range_m,
                max_drop_mil,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A representative metric load (a .308-class 168 gr match bullet at a 100 m zero),
    /// solved with a nonzero crosswind so both `drop_mil` and `wind_mil` are exercised.
    fn oracle_test_load() -> HoldCurveLoad {
        HoldCurveLoad {
            velocity_mps: 800.0,
            bc: 0.223,
            mass_kg: 0.0109,
            diameter_m: 0.00782,
            drag_model: DragModel::G7,
            sight_height_m: 0.045,
            zero_distance_m: 100.0,
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            humidity: 50.0,
            altitude_m: 0.0,
            wind_speed_mps: 3.0,
            wind_direction_deg: 90.0,
        }
    }

    /// Reproduces exactly what [`HoldCurve::solve`] does internally, as a SEPARATE, from-
    /// scratch call -- it never reads any `HoldCurve`'s private `samples` field. This gives
    /// `at_range` an independent direct-solve oracle to be checked against, rather than a
    /// self-consistency check of `HoldCurve` against its own stored state.
    fn independent_direct_solve(
        load: &HoldCurveLoad,
        max_range_m: f64,
    ) -> Vec<trajectory_sampling::TrajectorySample> {
        let zero_inputs = BallisticInputs {
            bc_value: load.bc,
            bc_type: load.drag_model,
            bullet_mass: load.mass_kg,
            muzzle_velocity: load.velocity_mps,
            bullet_diameter: load.diameter_m,
            bullet_length: fallback_bullet_length_m(load.diameter_m, load.mass_kg),
            sight_height: load.sight_height_m,
            use_rk4: true,
            ..Default::default()
        };
        let atmosphere = AtmosphericConditions {
            temperature: load.temperature_c,
            pressure: load.pressure_hpa,
            humidity: load.humidity,
            altitude: load.altitude_m,
        };
        let zero_angle = crate::calculate_zero_angle_with_conditions(
            zero_inputs,
            load.zero_distance_m,
            load.sight_height_m,
            WindConditions::default(),
            atmosphere,
        )
        .expect("zero solve should succeed for a realistic load");

        run_sampled_trajectory(
            load.velocity_mps,
            load.bc,
            load.mass_kg,
            load.diameter_m,
            load.drag_model,
            load.sight_height_m,
            load.temperature_c,
            load.pressure_hpa,
            load.humidity,
            load.altitude_m,
            load.wind_speed_mps,
            load.wind_direction_deg,
            max_range_m,
            HoldCurve::SAMPLE_INTERVAL_M,
            zero_angle,
            None,
            None,
            None,
            0.0,
            0.0,
            0.0,
            Some(load.zero_distance_m),
        )
        .expect("sampled trajectory should succeed for a realistic load")
    }

    #[test]
    fn at_range_matches_an_independent_direct_solves_sampled_observation() {
        let load = oracle_test_load();
        let max_range_m = 1500.0;

        let hold_curve = HoldCurve::solve(&load, max_range_m).expect("hold curve should solve");
        let oracle_samples = independent_direct_solve(&load, max_range_m);

        // A probe range read off the oracle's own grid rather than a hand-picked constant --
        // it only needs to land exactly on a sample node so `at_range`'s bracket search
        // resolves with t == 1.0 (the exact-endpoint-hit convention Task 7 disclosed: an exact
        // key hit lands on the upper endpoint of its bracket, not the lower one) and its own
        // linear interpolation is not itself a source of disagreement with the oracle.
        let probe = &oracle_samples[oracle_samples.len() / 2];
        let probe_range_m = probe.distance_m;
        assert!(probe_range_m > 0.0 && probe_range_m.is_finite());

        let point = hold_curve
            .at_range(probe_range_m)
            .expect("probe range must be inside the sampled span");

        let expected_drop_mil = probe.drop_m / probe_range_m * 1000.0;
        let expected_wind_mil = probe.wind_drift_m / probe_range_m * 1000.0;

        assert!((point.range_m - probe_range_m).abs() < 1e-9);
        assert!((point.drop_mil - expected_drop_mil).abs() < 1e-9);
        assert!((point.wind_mil - expected_wind_mil).abs() < 1e-9);
        assert!((point.velocity_mps - probe.velocity_mps).abs() < 1e-9);
        assert!((point.energy_j - probe.energy_j).abs() < 1e-9);
        assert!((point.time_s - probe.time_s).abs() < 1e-9);
    }

    #[test]
    fn at_range_outside_the_sampled_span_returns_none() {
        let load = oracle_test_load();
        let max_range_m = 1500.0;
        let hold_curve = HoldCurve::solve(&load, max_range_m).expect("hold curve should solve");

        let beyond = hold_curve.max_sampled_range_m() + 10_000.0;
        assert_eq!(hold_curve.at_range(beyond), None);
    }

    /// `sample_ranges_m` must reproduce the exact arithmetic sequence `HoldCurve::solve` built
    /// internally: index 0 at the muzzle (`0.0`), every later index an exact multiple of
    /// `SAMPLE_INTERVAL_M`, and the last entry equal to `max_sampled_range_m()`.
    #[test]
    fn sample_ranges_m_is_the_exact_multiples_of_the_sample_interval() {
        let load = oracle_test_load();
        let hold_curve = HoldCurve::solve(&load, 1500.0).expect("hold curve should solve");

        let ranges = hold_curve.sample_ranges_m();
        assert!(!ranges.is_empty());
        assert_eq!(ranges[0], 0.0);
        assert_eq!(*ranges.last().unwrap(), hold_curve.max_sampled_range_m());
        for (i, &r) in ranges.iter().enumerate() {
            let expected = i as f64 * HoldCurve::SAMPLE_INTERVAL_M;
            assert!(
                (r - expected).abs() < 1e-9,
                "index {i}: got {r}, expected {expected} ({} * SAMPLE_INTERVAL_M)",
                i
            );
        }
    }
}
