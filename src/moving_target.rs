//! Moving-target lead calculation (MBA-1287).
//!
//! Computes the hold ("lead") required to hit a target moving at constant speed along a
//! straight ground track, using the engine's own wind-aware time of flight.
//!
//! # Conventions (locked — see CLI_USAGE.md "Moving-Target Lead")
//! - `angle_deg` is the direction of target TRAVEL relative to the line of sight:
//!   `0` = directly away (outbound), `90` = full crossing left→right,
//!   `180` = directly toward (inbound), `270` = crossing right→left.
//! - Positive lead = hold in the target's direction of travel (right for `90°`).
//! - `lead_mil = lead/range·1000`; `lead_moa = lead/range·3438.0` — the CLI dial
//!   convention shared with the dope-card/come-up tables (MBA-724), so MOA/MIL is
//!   exactly 3.438, not the exact-angle 3437.7467.
//! - The lead is PURE target motion, additive to the wind-corrected hold: time of
//!   flight comes from the wind-aware solve, but wind drift itself stays in the
//!   wind column of your dope, not in the lead.
//! - Non-perpendicular motion changes the intercept distance; `calculate_lead`
//!   iterates `R = R₀ + v_radial·TOF(R)` until the correction is below 0.1 m and
//!   reports TOF/angular lead at the corrected (intercept) range.

use crate::cli_api::{
    AtmosphericConditions, BallisticInputs, BallisticsError, TrajectoryPoint, TrajectorySolver,
    WindConditions,
};
use std::fmt;

/// Milliradian hold per unit offset/range ratio.
pub const MIL_PER_UNIT_RATIO: f64 = 1000.0;
/// MOA hold per unit offset/range ratio — the engine's CLI dial convention (MBA-724),
/// deliberately 3438.0 (not 3437.7467) to match every printed table.
pub const MOA_PER_UNIT_RATIO: f64 = 3438.0;

/// Convergence tolerance for the intercept-range iteration, meters.
const RANGE_TOLERANCE_M: f64 = 0.1;
/// Iteration cap; the ticket's acceptance demands convergence in <10 for 45° outbound.
const MAX_ITERATIONS: u32 = 10;

/// Pure angular/linear lead from an already-known time of flight (no solve).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeadComponents {
    /// Signed linear lead at the target, meters (positive = direction of travel).
    pub lead_m: f64,
    /// Signed angular lead, milliradians.
    pub lead_mil: f64,
    /// Signed angular lead, MOA (CLI dial convention: ratio × 3438.0).
    pub lead_moa: f64,
}

/// Full lead solution from `calculate_lead`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeadSolution {
    /// Time of flight to the corrected (intercept) range, seconds — wind-aware.
    pub time_of_flight_s: f64,
    /// Signed linear lead at the target, meters (positive = direction of travel).
    pub lead_m: f64,
    /// Signed angular lead, milliradians.
    pub lead_mil: f64,
    /// Signed angular lead, MOA (ratio × 3438.0 — see module docs).
    pub lead_moa: f64,
    /// Intercept distance after radial-motion correction, meters.
    pub corrected_range_m: f64,
    /// Fixed-point iterations used (0 when the target has no radial component).
    pub iterations: u32,
}

/// Typed failure modes for `calculate_lead`.
#[derive(Debug)]
pub enum LeadError {
    /// Non-finite or out-of-domain input (negative speed, range ≤ 0, non-finite angle).
    InvalidInput(String),
    /// An inbound target closes to (or past) the shooter before intercept: the
    /// corrected range collapsed to ≤ 0 m.
    TargetOvertakesShooter { corrected_range_m: f64 },
    /// The intercept-range iteration failed to reach `RANGE_TOLERANCE_M` within
    /// `MAX_ITERATIONS` (physically requires a radial speed comparable to the
    /// projectile's remaining velocity).
    Convergence { iterations: u32, residual_m: f64 },
    /// The underlying trajectory solve failed (invalid ballistic inputs, etc.).
    Solver(BallisticsError),
}

impl fmt::Display for LeadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LeadError::InvalidInput(msg) => write!(f, "invalid lead input: {msg}"),
            LeadError::TargetOvertakesShooter { corrected_range_m } => write!(
                f,
                "target reaches the shooter before intercept (corrected range {corrected_range_m:.1} m)"
            ),
            LeadError::Convergence { iterations, residual_m } => write!(
                f,
                "intercept range failed to converge after {iterations} iterations (residual {residual_m:.2} m)"
            ),
            LeadError::Solver(e) => write!(f, "trajectory solve failed: {e}"),
        }
    }
}

impl std::error::Error for LeadError {}

impl From<BallisticsError> for LeadError {
    fn from(e: BallisticsError) -> Self {
        LeadError::Solver(e)
    }
}

/// Signed lateral/radial split of the target velocity, m/s.
/// Radial positive = receding; lateral positive = moving right.
fn velocity_components(speed_mps: f64, angle_deg: f64) -> (f64, f64) {
    let a = angle_deg.to_radians();
    (speed_mps * a.cos(), speed_mps * a.sin())
}

/// Pure lead math from a known time of flight — shared by `calculate_lead` and the
/// PDF dope-card path (which already has per-row TOF from its sampled trajectory).
/// Uses the small-angle linear ratio `lead/range`, matching every other hold table.
pub fn lead_from_tof(speed_mps: f64, angle_deg: f64, tof_s: f64, range_m: f64) -> LeadComponents {
    let (_, v_lateral) = velocity_components(speed_mps, angle_deg);
    let lead_m = v_lateral * tof_s;
    let ratio = if range_m < 1.0 { 0.0 } else { lead_m / range_m };
    LeadComponents {
        lead_m,
        lead_mil: ratio * MIL_PER_UNIT_RATIO,
        lead_moa: ratio * MOA_PER_UNIT_RATIO,
    }
}

/// Linear-interpolated time of flight at downrange distance `x_m`, from a solved
/// trajectory's points. Returns `None` when `x_m` lies beyond the solved span.
/// (MBA-1218 guarantees the final point sits exactly at the solve's max range.)
fn tof_at(points: &[TrajectoryPoint], x_m: f64) -> Option<f64> {
    if points.is_empty() || x_m < 0.0 {
        return None;
    }
    if x_m <= points[0].position.x {
        return Some(points[0].time);
    }
    for w in points.windows(2) {
        let (p1, p2) = (&w[0], &w[1]);
        if p2.position.x >= x_m {
            let dx = p2.position.x - p1.position.x;
            let t = if dx.abs() < 1e-12 { 0.0 } else { (x_m - p1.position.x) / dx };
            return Some(p1.time + t * (p2.time - p1.time));
        }
    }
    None
}

/// Compute the lead solution for a target moving at `speed_mps` along `angle_deg`
/// (see module docs for the datum) at initial range `range_m`.
///
/// Runs ONE wind-aware solve past the range (with outbound headroom), interpolates
/// time of flight from its points, and fixed-point iterates the intercept range for
/// radial motion. Errors are typed (`LeadError`).
pub fn calculate_lead(
    inputs: BallisticInputs,
    wind: WindConditions,
    atmo: AtmosphericConditions,
    speed_mps: f64,
    angle_deg: f64,
    range_m: f64,
) -> Result<LeadSolution, LeadError> {
    if !speed_mps.is_finite() || speed_mps < 0.0 {
        return Err(LeadError::InvalidInput(format!(
            "target speed must be finite and non-negative, got {speed_mps}"
        )));
    }
    if !angle_deg.is_finite() {
        return Err(LeadError::InvalidInput(format!(
            "target angle must be finite, got {angle_deg}"
        )));
    }
    if !range_m.is_finite() || range_m <= 0.0 {
        return Err(LeadError::InvalidInput(format!(
            "range must be finite and positive, got {range_m}"
        )));
    }

    let (v_radial, _) = velocity_components(speed_mps, angle_deg);

    // One solve with headroom for outbound intercepts: 6 s of target motion covers
    // any realistic small-arms TOF at leadable ranges.
    let max_range = range_m + (v_radial.max(0.0) * 6.0).max(20.0);
    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(max_range);
    let result = solver.solve()?;
    let points = &result.points;

    let beyond_solved = |r: f64| LeadError::Convergence {
        iterations: MAX_ITERATIONS,
        residual_m: r - max_range,
    };

    // Fixed-point iteration on the intercept range (no-op for pure crossing motion).
    // The radial guard uses an epsilon, NOT != 0.0: cos(90°.to_radians()) is ~6.1e-17
    // in floating point, and a bare != 0.0 would make every perpendicular call "iterate"
    // once (reporting iterations = 1 and a ~1e-14 m range shift) — the acceptance test
    // requires iterations == 0 for pure crossing. 1e-9 m/s of radial speed is far below
    // physical relevance.
    let mut corrected = range_m;
    let mut iterations = 0u32;
    if v_radial.abs() > 1e-9 && speed_mps > 0.0 {
        loop {
            let tof = tof_at(points, corrected).ok_or_else(|| beyond_solved(corrected))?;
            let next = range_m + v_radial * tof;
            if next <= 0.0 {
                return Err(LeadError::TargetOvertakesShooter { corrected_range_m: next });
            }
            let residual = (next - corrected).abs();
            iterations += 1;
            corrected = next;
            if residual < RANGE_TOLERANCE_M {
                break;
            }
            if iterations >= MAX_ITERATIONS {
                return Err(LeadError::Convergence { iterations, residual_m: residual });
            }
        }
    }

    let tof = tof_at(points, corrected).ok_or_else(|| beyond_solved(corrected))?;
    let components = lead_from_tof(speed_mps, angle_deg, tof, corrected);
    Ok(LeadSolution {
        time_of_flight_s: tof,
        lead_m: components.lead_m,
        lead_mil: components.lead_mil,
        lead_moa: components.lead_moa,
        corrected_range_m: corrected,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtmosphericConditions, BallisticInputs, DragModel, WindConditions};

    fn base_inputs() -> BallisticInputs {
        let mut i = BallisticInputs::default();
        i.muzzle_velocity = 800.0;
        i.bc_value = 0.5;
        i.bc_type = DragModel::G7;
        i.bullet_mass = 0.0109;
        i.bullet_diameter = 0.00782;
        i.bullet_length = 0.0309;
        i.sight_height = 0.05;
        i.use_rk4 = true;
        i
    }

    #[test]
    fn lead_from_tof_perpendicular_is_speed_times_tof() {
        let c = lead_from_tof(3.0, 90.0, 0.5, 300.0);
        assert!((c.lead_m - 1.5).abs() < 1e-12);
        assert!((c.lead_mil - 1.5 / 300.0 * 1000.0).abs() < 1e-12);
        assert!((c.lead_moa - 1.5 / 300.0 * 3438.0).abs() < 1e-12);
    }

    #[test]
    fn moa_mil_ratio_is_exactly_3_438() {
        let c = lead_from_tof(5.0, 90.0, 0.7, 400.0);
        assert!((c.lead_moa / c.lead_mil - 3.438).abs() < 1e-12);
    }

    #[test]
    fn pure_radial_motion_has_zero_lead() {
        for angle in [0.0, 180.0] {
            let c = lead_from_tof(10.0, angle, 0.5, 300.0);
            assert!(c.lead_m.abs() < 1e-9, "angle {angle} must give zero lead");
        }
    }

    #[test]
    fn right_to_left_lead_is_negative() {
        let c = lead_from_tof(3.0, 270.0, 0.5, 300.0);
        assert!(c.lead_m < 0.0, "270 deg = target moving left => negative (hold left)");
    }

    #[test]
    fn zero_speed_gives_zero_lead_solution() {
        let s = calculate_lead(
            base_inputs(),
            WindConditions::default(),
            AtmosphericConditions::default(),
            0.0,
            90.0,
            300.0,
        )
        .expect("zero speed must solve");
        assert_eq!(s.lead_m, 0.0);
        assert_eq!(s.lead_mil, 0.0);
        assert!((s.corrected_range_m - 300.0).abs() < 1e-9);
    }

    #[test]
    fn outbound_45_converges_fast_and_grows_range() {
        let s = calculate_lead(
            base_inputs(),
            WindConditions::default(),
            AtmosphericConditions::default(),
            15.0,
            45.0,
            600.0,
        )
        .expect("outbound must converge");
        assert!(s.iterations < 10, "iterations {}", s.iterations);
        assert!(s.corrected_range_m > 600.0, "outbound target => longer intercept");
        // residual check: one more application of the map moves R by < 0.1 m
    }

    #[test]
    fn inbound_gives_shorter_corrected_range() {
        let s = calculate_lead(
            base_inputs(),
            WindConditions::default(),
            AtmosphericConditions::default(),
            15.0,
            180.0,
            600.0,
        )
        .expect("inbound must converge");
        assert!(s.corrected_range_m < 600.0);
        assert!(s.lead_m.abs() < 1e-9, "pure inbound has no lateral lead");
    }

    #[test]
    fn over_closing_target_is_a_typed_error() {
        // Radial closing speed so large the intercept range collapses to <= 0.
        let r = calculate_lead(
            base_inputs(),
            WindConditions::default(),
            AtmosphericConditions::default(),
            2000.0,
            180.0,
            600.0,
        );
        match r {
            Err(LeadError::TargetOvertakesShooter { .. }) | Err(LeadError::Convergence { .. }) => {}
            other => panic!("expected typed overtake/convergence error, got {other:?}"),
        }
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let bad = |speed: f64, angle: f64, range: f64| {
            calculate_lead(
                base_inputs(),
                WindConditions::default(),
                AtmosphericConditions::default(),
                speed,
                angle,
                range,
            )
        };
        assert!(matches!(bad(f64::NAN, 90.0, 300.0), Err(LeadError::InvalidInput(_))));
        assert!(matches!(bad(-1.0, 90.0, 300.0), Err(LeadError::InvalidInput(_))));
        assert!(matches!(bad(3.0, f64::INFINITY, 300.0), Err(LeadError::InvalidInput(_))));
        assert!(matches!(bad(3.0, 90.0, 0.0), Err(LeadError::InvalidInput(_))));
    }
}
