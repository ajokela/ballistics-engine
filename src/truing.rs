//! Velocity / BC truing core (MBA-1343 Phase A).
//!
//! The multi-observation joint MV+BC calibration (MBA-1316) extracted from the
//! CLI binary so non-CLI front ends (e.g. the WASM terminal) can reuse the
//! exact compute path. All rendering (table / JSON / CSV) stays with the front
//! ends; this module goes as far as building a [`MultiTruingReport`]. The
//! module is silent: it never writes to stdout/stderr, so progress reporting
//! is the caller's job (see [`validate_truing_observations`] for how the CLI
//! sequences its progress line without double-reporting validation errors).

use std::error::Error;

use serde::Serialize;

use crate::cli_api::UnitSystem;
use crate::constants::GRAINS_TO_KG;
use crate::{
    AtmosphericConditions, BCSegmentData, BallisticInputs, DragModel, TrajectorySolver,
    WindConditions,
};

/// Drag-model selector for the truing commands, in the shape the CLI exposes
/// (a plain G1/G7 choice, lowercase-parsed by clap/the WASM terminal). It maps
/// onto the engine's [`DragModel`] inside the solvers; unlike [`DragModel`] it
/// deliberately offers only the two standard models the truing forward model
/// supports.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum DragModelArg {
    G1,
    G7,
}

/// Unit in which observed drops are supplied to (and residuals reported from)
/// `true-velocity`'s multi-observation joint calibration (MBA-1316). The
/// historical single-observation path is always MIL, so `mil` is the default and
/// leaves that path byte-identical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum DropUnit {
    /// Milliradians (angular)
    Mil,
    /// Minutes of angle (angular, true 1/60°)
    Moa,
    /// Inches below line of sight (linear drop at the target)
    In,
}

impl DropUnit {
    /// Short label used in tables/JSON/CSV.
    pub fn label(self) -> &'static str {
        match self {
            DropUnit::Mil => "mil",
            DropUnit::Moa => "moa",
            DropUnit::In => "in",
        }
    }

    /// Parse a drop-unit token (`"mil"` | `"moa"` | `"in"`, case-insensitive)
    /// without going through clap, for non-CLI front ends (e.g. the WASM
    /// terminal parser).
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "mil" => Ok(DropUnit::Mil),
            "moa" => Ok(DropUnit::Moa),
            "in" => Ok(DropUnit::In),
            _ => Err(format!("invalid drop unit '{s}': expected mil, moa, or in")),
        }
    }

    /// Express a linear vertical drop below the line of sight (`drop_m`, meters)
    /// at horizontal distance `z_m` (meters) in this unit. Angular units use the
    /// same small-angle convention the engine's MIL output has always used
    /// (`drop/range`), so `mil` matches the legacy single-observation contract
    /// exactly; `moa` is true minutes of angle and `in` is the linear drop itself.
    pub fn express_drop_m(self, drop_m: f64, z_m: f64) -> f64 {
        match self {
            // (drop_m / z_m) radians * 1000 mrad/rad
            DropUnit::Mil => (drop_m / z_m) * 1000.0,
            // (drop_m / z_m) radians * (180/pi) deg/rad * 60 min/deg
            DropUnit::Moa => (drop_m / z_m) * (180.0 / std::f64::consts::PI) * 60.0,
            // linear drop, meters -> inches
            DropUnit::In => drop_m / 0.0254,
        }
    }
}

/// Resolve a bullet length (meters) for a bullet whose length the shooter did not supply (MBA-1135).
///
/// Prefers the mass-based physical estimate ([`crate::stability::estimate_bullet_length_m`]),
/// falling back to the historical 4.5-caliber heuristic only when mass is unavailable so the caller
/// always gets a positive length. `diameter_m` and `mass_kg` are SI.
pub fn fallback_bullet_length_m(diameter_m: f64, mass_kg: f64) -> f64 {
    let est = crate::stability::estimate_bullet_length_m(diameter_m, mass_kg);
    if est > 0.0 {
        est
    } else {
        diameter_m * 4.5
    }
}

/// Shared solver assembly for the truing forward model (MBA-1316; extracted from
/// [`solve_trajectory_drop`] for MBA-1405 so a second caller can read off the full
/// [`crate::TrajectoryResult`] — e.g. its labeled Mach crossings — instead of just
/// the drop/range pair. SI conversion, base [`BallisticInputs`], zero-angle solve,
/// and `max_range`/`time_step` are all identical to what `solve_trajectory_drop`
/// has always used, so both callers stay physically identical.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing velocity-truing compatibility helper"
)]
fn solve_truing_trajectory(
    velocity_fps: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    range_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<crate::TrajectoryResult, Box<dyn Error>> {
    // Convert to SI units
    let velocity_ms = velocity_fps * 0.3048;
    let mass_kg = mass_gr * GRAINS_TO_KG;
    let diameter_m = diameter_in * 0.0254;
    let zero_m = zero_distance_yd * 0.9144;
    let range_m = range_yd * 0.9144;
    let sight_height_m = sight_height_in * 0.0254;
    let altitude_m = altitude_ft * 0.3048;
    let temperature_c = (temperature_f - 32.0) * 5.0 / 9.0;
    let pressure_hpa = pressure_inhg * 33.8639; // Convert inHg to hPa

    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    // Create base inputs - match defaults used by trajectory command
    let mut inputs = BallisticInputs {
        muzzle_velocity: velocity_ms,
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        bullet_diameter: diameter_m,
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg), // MBA-1135 mass-based estimate
        sight_height: sight_height_m,
        target_distance: range_m + 100.0, // Overshoot to ensure we have data
        use_bc_segments: bc_segments.is_some(),
        bc_segments_data: bc_segments.clone(),
        use_rk4: true,
        muzzle_angle: 0.0,    // Will be set by zero angle calculation
        ..Default::default()  // Uses muzzle_height: 0.0 by default
    };

    // Set up atmospheric conditions
    // AtmosphericConditions expects: temperature in Celsius, pressure in hPa, humidity 0-100, altitude in meters
    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity, // Already 0-100 from input
        altitude: altitude_m,
    };

    let wind = WindConditions::default();

    // Calculate zero angle for the zero distance
    // Target height is sight_height because the bullet must cross the LOS at zero distance
    // The LOS is at y = sight_height (sight is above bore by sight_height)
    // So the bullet (starting at y = 0 = bore level) must rise to y = sight_height at zero distance
    let zero_angle = crate::calculate_zero_angle_with_conditions(
        inputs.clone(),
        zero_m,
        sight_height_m, // target height at zero distance (LOS height)
        wind.clone(),
        atmosphere.clone(),
    )?;

    // Set the calculated zero angle
    inputs.muzzle_angle = zero_angle;

    // Create solver and solve
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(range_m + 100.0);
    solver.set_time_step(0.0001);

    Ok(solver.solve()?)
}

/// Shared forward-model core for velocity/BC truing (MBA-1316).
///
/// Solves the real trajectory for the given `(velocity_fps, bc)` candidate under
/// the supplied load/atmosphere and returns `(drop_m, z_m)` where `drop_m` is the
/// linear vertical distance below the line of sight (positive = below LOS) at the
/// target range and `z_m` is the horizontal distance actually reached (~range_m).
/// Callers convert to MIL / MOA / inches as needed. This is the exact assembly
/// the single-observation binary search has always used, factored out so the
/// multi-observation joint fit reuses identical physics.
///
/// When `interpolate` is false the returned point is the first trajectory sample
/// at or past the target range (the historical behaviour). When true, the drop
/// and horizontal distance are linearly interpolated to land exactly on the
/// target range, removing the ~v*dt spatial quantization that otherwise
/// stair-steps the cost surface — essential for the multi-observation joint fit,
/// whose finite-difference Jacobian would otherwise be dominated by that noise.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing velocity-truing compatibility helper"
)]
pub(crate) fn solve_trajectory_drop(
    velocity_fps: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    range_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
    interpolate: bool,
) -> Result<(f64, f64), Box<dyn Error>> {
    let range_m = range_yd * 0.9144;
    let sight_height_m = sight_height_in * 0.0254;

    let result = solve_truing_trajectory(
        velocity_fps,
        bc,
        drag_model,
        mass_gr,
        diameter_in,
        zero_distance_yd,
        range_yd,
        sight_height_in,
        temperature_f,
        pressure_inhg,
        humidity,
        altitude_ft,
        bc_segments,
    )?;

    // Find the first point at or past the target range.
    let idx = result
        .points
        .iter()
        .position(|p| p.position.x >= range_m)
        .ok_or("Trajectory didn't reach target range")?;

    // Determine the horizontal distance z and bullet height at the target range.
    let (z, bullet_y) = if interpolate && idx > 0 {
        // Linearly interpolate between the bracketing samples to land exactly on
        // range_m (removes ~v*dt spatial quantization).
        let p0 = &result.points[idx - 1];
        let p1 = &result.points[idx];
        let x0 = p0.position.x;
        let x1 = p1.position.x;
        let denom = x1 - x0;
        if denom.abs() < f64::EPSILON {
            (p1.position.x, p1.position.y)
        } else {
            let frac = (range_m - x0) / denom;
            let y = p0.position.y + frac * (p1.position.y - p0.position.y);
            (range_m, y)
        }
    } else {
        let p = &result.points[idx];
        (p.position.x, p.position.y)
    };

    // Calculate drop relative to the line of sight. The trajectory was already launched at the
    // solved zero angle, so the LOS for this zero solve is horizontal at sight_height_m.
    // Drop = LOS height - bullet position (positive = below LOS)
    let drop_m = sight_height_m - bullet_y;

    Ok((drop_m, z))
}

/// Calculate drop at a given muzzle velocity using trajectory solver
/// Returns drop in MILs at the target range
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments preserve the existing velocity-truing compatibility helper"
)]
pub(crate) fn calculate_drop_at_velocity(
    velocity_fps: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    range_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<f64, Box<dyn Error>> {
    // Preserve the historical MIL contract by delegating to the shared linear-drop
    // core (MBA-1316). `interpolate = false` reproduces the original "first point
    // at or past the range" sampling, so this path stays byte-identical.
    let (drop_m, z) = solve_trajectory_drop(
        velocity_fps,
        bc,
        drag_model,
        mass_gr,
        diameter_in,
        zero_distance_yd,
        range_yd,
        sight_height_in,
        temperature_f,
        pressure_inhg,
        humidity,
        altitude_ft,
        bc_segments,
        false,
    )?;

    // Convert to MILs: mil = (drop_inches / 36 / range_yards) * 1000
    // Or equivalently: mil = (drop_m / range_m) * 1000
    let drop_mil = (drop_m / z) * 1000.0;

    Ok(drop_mil)
}

/// Minimum sane chronograph screen distance from the muzzle, in meters (MBA-1377). Below this
/// the correction is negligible and a nonzero value more likely reflects a data-entry mistake
/// than a real chronograph placement, so [`correct_chrono_velocity_fps`] rejects it outright
/// rather than silently applying a near-zero (and therefore noise-dominated) adjustment.
pub const MIN_CHRONO_DISTANCE_M: f64 = 0.3; // ~1 ft

/// Maximum sane chronograph screen distance from the muzzle, in meters (MBA-1377).
/// Comfortably past Lapua/JBM's 25 m reference distance -- past this a "chronograph" reading
/// is no longer describing a screen/radar setup this correction is meant to undo.
pub const MAX_CHRONO_DISTANCE_M: f64 = 30.0; // ~98 ft

/// Result of back-solving a muzzle velocity from a chronograph reading taken downrange
/// (MBA-1377; see [`correct_chrono_velocity_fps`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChronoCorrection {
    /// The back-solved muzzle velocity, in feet per second.
    pub muzzle_velocity_fps: f64,
    /// Number of secant iterations the solve actually took.
    pub iterations: u32,
}

/// Forward-model core for the screen-distance correction (MBA-1377): launches
/// `muzzle_velocity_fps` on a flat (`muzzle_angle = 0`) shot through the SAME BC / drag model /
/// atmosphere the rest of the command is using, and returns the velocity magnitude at
/// `screen_distance_m` (linearly interpolated between the bracketing trajectory samples).
///
/// Deliberately independent of [`solve_truing_trajectory`]: that helper's zero-angle search and
/// sight-height bookkeeping are irrelevant here. Chronograph screen distances are short enough
/// (3-25 m; validated range [`MIN_CHRONO_DISTANCE_M`]..[`MAX_CHRONO_DISTANCE_M`]) that gravity
/// drop is negligible for the velocity this returns, so firing along the bore line is the
/// correct simplification for this calculation, not an approximation of convenience.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the rest of this module's forward-model helpers"
)]
fn velocity_at_distance_fps(
    muzzle_velocity_fps: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
    screen_distance_m: f64,
) -> Result<f64, Box<dyn Error>> {
    let velocity_ms = muzzle_velocity_fps * 0.3048;
    let mass_kg = mass_gr * GRAINS_TO_KG;
    let diameter_m = diameter_in * 0.0254;
    let altitude_m = altitude_ft * 0.3048;
    let temperature_c = (temperature_f - 32.0) * 5.0 / 9.0;
    let pressure_hpa = pressure_inhg * 33.8639; // inHg to hPa

    let drag_model_enum = match drag_model {
        DragModelArg::G1 => DragModel::G1,
        DragModelArg::G7 => DragModel::G7,
    };

    let inputs = BallisticInputs {
        muzzle_velocity: velocity_ms,
        bc_value: bc,
        bc_type: drag_model_enum,
        bullet_mass: mass_kg,
        bullet_diameter: diameter_m,
        bullet_length: fallback_bullet_length_m(diameter_m, mass_kg), // MBA-1135 mass-based estimate
        target_distance: screen_distance_m + 2.0, // overshoot so the sample bracket exists
        use_bc_segments: bc_segments.is_some(),
        bc_segments_data: bc_segments.clone(),
        use_rk4: true,
        muzzle_angle: 0.0, // bore-line: drop over 3-25 m is negligible for velocity
        ..Default::default()
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity,
        altitude: altitude_m,
    };

    let wind = WindConditions::default();

    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(screen_distance_m + 2.0);
    solver.set_time_step(0.0001);
    let result = solver.solve()?;

    let idx = result
        .points
        .iter()
        .position(|p| p.position.x >= screen_distance_m)
        .ok_or("trajectory did not reach the chronograph screen distance")?;

    if idx == 0 {
        return Ok(result.points[0].velocity_magnitude / 0.3048);
    }

    // Linearly interpolate to land exactly on screen_distance_m (same spatial-quantization fix
    // as solve_trajectory_drop's `interpolate = true` path) so the secant solver in
    // correct_chrono_velocity_fps sees a smooth function of the candidate muzzle velocity.
    let p0 = &result.points[idx - 1];
    let p1 = &result.points[idx];
    let x0 = p0.position.x;
    let x1 = p1.position.x;
    let v_ms = if (x1 - x0).abs() < f64::EPSILON {
        p1.velocity_magnitude
    } else {
        let frac = (screen_distance_m - x0) / (x1 - x0);
        p0.velocity_magnitude + frac * (p1.velocity_magnitude - p0.velocity_magnitude)
    };
    Ok(v_ms / 0.3048)
}

/// Back-solve a true muzzle velocity from a chronograph reading taken `screen_distance_m` (SI
/// meters) downrange (MBA-1377).
///
/// Most chronographs read 10-15 ft (or 25 m) downrange rather than at the muzzle, so the raw
/// reading is slightly LOW; JBM ("Distance to Chronograph"), Lapua (V25m to V0), Ballistic AE,
/// ColdBore, Remington Shoot!, and Hornady all correct for this. Flat-fire point-mass
/// deceleration is monotone in muzzle velocity over the validated screen-distance band
/// ([`MIN_CHRONO_DISTANCE_M`]..[`MAX_CHRONO_DISTANCE_M`]), so the secant method converges in a
/// handful of iterations (McCoy, *Modern Exterior Ballistics*; JBM's published calculator).
///
/// This is a pure input-side transform: `velocity_at_distance_fps` runs the EXISTING forward
/// trajectory model (the same BC / `drag_model` / atmosphere the rest of the command is
/// configured with) to predict the screen-distance velocity for a candidate muzzle velocity,
/// and the candidate is iterated until that prediction matches the measured reading. Callers
/// apply this once, before the corrected value is used for display/comparison --
/// `solve_truing_trajectory`'s drop-based solves never receive a screen distance.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the rest of this module's forward-model helpers"
)]
pub fn correct_chrono_velocity_fps(
    measured_velocity_fps: f64,
    screen_distance_m: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<ChronoCorrection, Box<dyn Error>> {
    if !measured_velocity_fps.is_finite() || measured_velocity_fps <= 0.0 {
        return Err("--chrono-velocity must be positive and finite".into());
    }
    if !(MIN_CHRONO_DISTANCE_M..=MAX_CHRONO_DISTANCE_M).contains(&screen_distance_m) {
        return Err(format!(
            "--chrono-distance is out of range: {screen_distance_m:.3} m downrange (expected \
             {MIN_CHRONO_DISTANCE_M}..{MAX_CHRONO_DISTANCE_M} m; most chronographs read 10-15 \
             ft / ~3-5 m, and Lapua/JBM's reference distance is 25 m)"
        )
        .into());
    }

    // Residual: forward-model screen velocity for candidate muzzle velocity `v0`, minus the
    // measured reading. `correct_chrono_velocity_fps` iterates on `v0` until this is ~0. (Not
    // JS/Python `eval` -- a plain Rust closure over a local candidate value.)
    let residual = |v0: f64| -> Result<f64, Box<dyn Error>> {
        let v_screen = velocity_at_distance_fps(
            v0,
            bc,
            drag_model,
            mass_gr,
            diameter_in,
            temperature_f,
            pressure_inhg,
            humidity,
            altitude_ft,
            bc_segments,
            screen_distance_m,
        )?;
        Ok(v_screen - measured_velocity_fps)
    };

    // Secant method. Drag only decelerates, so the true muzzle velocity always exceeds the
    // downrange reading -- seed the second guess 1% above the first, comfortably inside the
    // near-linear region for the sub-percent corrections this band produces.
    let mut v_prev = measured_velocity_fps;
    let mut f_prev = residual(v_prev)?;
    let mut v_curr = measured_velocity_fps * 1.01;
    let mut f_curr = residual(v_curr)?;

    const TOL_FPS: f64 = 1e-4;
    const MAX_ITERATIONS: u32 = 20;
    let mut iterations = 0u32;

    while f_curr.abs() > TOL_FPS && iterations < MAX_ITERATIONS {
        let denom = f_curr - f_prev;
        if denom.abs() < f64::EPSILON {
            return Err("chronograph correction failed to converge (stalled secant step)".into());
        }
        let v_next = v_curr - f_curr * (v_curr - v_prev) / denom;
        v_prev = v_curr;
        f_prev = f_curr;
        v_curr = v_next;
        f_curr = residual(v_curr)?;
        iterations += 1;
    }

    if f_curr.abs() > TOL_FPS {
        return Err(format!(
            "chronograph correction did not converge after {iterations} iterations"
        )
        .into());
    }

    Ok(ChronoCorrection {
        muzzle_velocity_fps: v_curr,
        iterations,
    })
}

/// Result of the classic single-observation velocity truing
/// ([`calculate_true_velocity_local`]).
#[derive(Debug, Clone)]
pub struct TrueVelocityLocalResult {
    /// The fitted effective muzzle velocity, in feet per second.
    pub effective_velocity_fps: f64,
    /// Number of binary-search iterations actually run.
    pub iterations: i32,
    /// Signed residual at the returned velocity, in MIL:
    /// `calculated_drop_mil - measured_drop_mil`. Positive means the model still
    /// predicts MORE drop than was measured at the returned velocity; negative
    /// means less.
    pub final_error_mil: f64,
    /// The model-predicted drop (MIL) at the returned velocity.
    pub calculated_drop_mil: f64,
    /// Convergence quality: `"high"` (converged, |error| < 0.005 mil),
    /// `"medium"` (converged within the 0.01 mil tolerance, or stopped early
    /// with |error| < 0.1 mil), or `"low"` (did not converge; |error| >= 0.1 mil).
    pub confidence: String,
}

/// Calculate the effective muzzle velocity that reproduces `measured_drop_mil`
/// at `range_yd`, via binary search over the real trajectory solver
/// (the classic single-observation `true-velocity` path).
///
/// Returns a [`TrueVelocityLocalResult`] carrying the fitted velocity (fps),
/// the iteration count, the signed final residual in MIL (positive = the model
/// still predicts more drop than measured), the model-predicted drop at the
/// returned velocity, and a `"high"`/`"medium"`/`"low"` confidence label — see
/// the field docs on [`TrueVelocityLocalResult`] for the exact banding.
///
/// Errors on degenerate inputs (`range_yd` non-positive or non-finite,
/// `measured_drop_mil` non-finite) and on any trajectory-solver failure.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable true-velocity CLI command shape"
)]
pub fn calculate_true_velocity_local(
    measured_drop_mil: f64,
    range_yd: f64,
    bc: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_distance_yd: f64,
    sight_height_in: f64,
    temperature_f: f64,
    pressure_inhg: f64,
    humidity: f64,
    altitude_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<TrueVelocityLocalResult, Box<dyn Error>> {
    // Reject degenerate inputs instead of letting NaN/inf flow through the
    // solver and come back as Ok(NaN). The native CLI's clap range validators
    // cannot produce these; direct library / WASM callers can.
    if !range_yd.is_finite() || range_yd <= 0.0 {
        return Err("range must be positive and finite".into());
    }
    if !measured_drop_mil.is_finite() {
        return Err("measured drop must be finite".into());
    }

    // Binary search between velocity bounds
    let mut velocity_low = 1500.0;
    let mut velocity_high = 4500.0;
    let tolerance_mil = 0.01; // 0.01 MIL tolerance
    let max_iterations = 50;

    let mut iterations = 0;
    let mut last_error = 0.0;
    let mut last_calculated_drop = 0.0;

    for i in 0..max_iterations {
        iterations = i + 1;
        let test_velocity = (velocity_low + velocity_high) / 2.0;

        // Run trajectory at test velocity
        let calculated_drop_mil = calculate_drop_at_velocity(
            test_velocity,
            bc,
            drag_model,
            mass_gr,
            diameter_in,
            zero_distance_yd,
            range_yd,
            sight_height_in,
            temperature_f,
            pressure_inhg,
            humidity,
            altitude_ft,
            bc_segments,
        )?;

        last_calculated_drop = calculated_drop_mil;
        let error = calculated_drop_mil - measured_drop_mil;
        last_error = error;

        if error.abs() < tolerance_mil {
            // Converged
            let confidence = if error.abs() < 0.005 {
                "high"
            } else {
                "medium"
            };

            return Ok(TrueVelocityLocalResult {
                effective_velocity_fps: test_velocity,
                iterations,
                final_error_mil: error,
                calculated_drop_mil,
                confidence: confidence.to_string(),
            });
        }

        // Higher calculated drop = bullet is slower = need higher velocity
        // Lower calculated drop = bullet is faster = need lower velocity
        if calculated_drop_mil > measured_drop_mil {
            // Bullet dropping more than observed = slower than actual
            // Need higher velocity
            velocity_low = test_velocity;
        } else {
            // Bullet dropping less = faster than actual
            // Need lower velocity
            velocity_high = test_velocity;
        }

        // Check for convergence issues
        if (velocity_high - velocity_low).abs() < 0.5 {
            break;
        }
    }

    // Did not converge within tolerance, return best estimate
    let final_velocity = (velocity_low + velocity_high) / 2.0;
    let confidence = if last_error.abs() < 0.1 {
        "medium"
    } else {
        "low"
    };

    Ok(TrueVelocityLocalResult {
        effective_velocity_fps: final_velocity,
        iterations,
        final_error_mil: last_error,
        calculated_drop_mil: last_calculated_drop,
        confidence: confidence.to_string(),
    })
}

// ============================================================================
// MBA-1316: multi-observation joint MV + BC calibration (truing v2)
// ============================================================================
//
// The classic `true-velocity` path fits muzzle velocity from a single observed
// drop while BC is held fixed. Mid-range (fully supersonic) drops are dominated
// by time of flight and therefore constrain muzzle velocity; long-range /
// transonic drops are where BC bites. With observations spread across those
// regimes we can separate the two. When the spread is too narrow to tell them
// apart we refuse the joint fit and fit muzzle velocity only, saying so.

// Solver / fit bounds and identifiability gates. See each constant's doc; the
// empirical basis for the gate values (calibrated against real .30-cal data,
// MBA-1316):
//   * 300/600/900 -> sens 0.29, cond 112 -> BC recovered to ~2%   (JOINT)
//   * 300/600/900/1000 -> sens 0.33, cond 214 -> BC to <0.1%      (JOINT)
//   * 300/400/500 -> sens 0.16, cond 282 -> BC error ~3%          (MV-ONLY)
//   * 200/250/300 -> sens 0.12, cond 2337 -> BC error ~12%        (MV-ONLY)
// (.308 168gr @ 2700/0.475, ~0.03 mil observation noise.) The gates sit between
// the "recovers to ~2%" band and the "several-percent garbage" band. They are
// heuristics: a set that just clears them still carries more BC uncertainty
// than a set that clears them comfortably.

/// Lower muzzle-velocity fit bound (fps): brackets subsonic pistol loads.
pub(crate) const TRUING_MV_MIN_FPS: f64 = 1000.0;
/// Upper muzzle-velocity fit bound (fps): brackets hyper-velocity varmint loads.
pub(crate) const TRUING_MV_MAX_FPS: f64 = 5000.0;
/// Lower BC fit bound; matches the CLI's `--bc` validator minimum.
pub(crate) const TRUING_BC_MIN: f64 = 0.05;
/// Upper BC fit bound; matches the CLI's `--bc` validator maximum.
pub(crate) const TRUING_BC_MAX: f64 = 2.0;
/// Iteration cap shared by the MV-only and joint fitters.
pub(crate) const TRUING_MAX_ITERS: usize = 40;

/// Identifiability gate: minimum
/// `sensitivity_ratio = ||bc*d(drop)/d(bc)|| / ||mv*d(drop)/d(mv)||` — how much
/// a fractional BC change moves the predicted drops relative to a fractional MV
/// change. Below this, the observations barely constrain BC and the joint fit
/// is refused (MV-only fallback).
pub(crate) const TRUING_MIN_BC_SENSITIVITY_RATIO: f64 = 0.20;
/// Identifiability gate: maximum `condition_number = (1+|c|)/(1-|c|)` (|c| =
/// correlation of the raw MV/BC Jacobian columns). Above this the two columns
/// are collinear — a BC change can be undone by an MV change — and the joint
/// fit is refused (MV-only fallback). Both this gate and
/// [`TRUING_MIN_BC_SENSITIVITY_RATIO`] must pass to attempt the joint fit.
pub(crate) const TRUING_MAX_CONDITION_NUMBER: f64 = 1.0e3;

/// A single observed impact used for truing: range (internal yards) and the
/// measured drop below line of sight expressed in the caller's drop unit.
#[derive(Debug, Clone, Copy)]
pub struct TruingObservation {
    pub range_yd: f64,
    pub drop: f64,
}

/// Complete scalar-BC forward-model input used by the truing design and
/// uncertainty APIs added in MBA-1346/MBA-1353.
///
/// All values use the truing core's historical imperial units.  A front end is
/// expected to convert its display units once at the boundary.  V1 deliberately
/// models one scalar G1/G7 coefficient: a velocity-banded BC schedule or custom
/// Mach/Cd deck does not have a single BC parameter to identify and must not be
/// silently reduced to one.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct TruingModelInputsV1 {
    /// Nominal muzzle velocity about which a plan is designed, in feet/second.
    pub muzzle_velocity_fps: f64,
    /// Nominal scalar ballistic coefficient about which a plan is designed.
    pub ballistic_coefficient: f64,
    pub drag_model: DragModelArg,
    /// Bullet mass in grains.
    pub mass_gr: f64,
    /// Bullet diameter in inches.
    pub diameter_in: f64,
    /// Zero distance in yards.
    pub zero_distance_yd: f64,
    /// Sight height over bore in inches.
    pub sight_height_in: f64,
    /// Ambient temperature in degrees Fahrenheit.
    pub temperature_f: f64,
    /// Station pressure in inches of mercury.
    pub pressure_inhg: f64,
    /// Relative humidity in percent (0 through 100).
    pub humidity_pct: f64,
    /// Altitude in feet.
    pub altitude_ft: f64,
}

impl TruingModelInputsV1 {
    /// Validate the scalar-BC model before an expensive design or uncertainty
    /// calculation.  These bounds mirror the existing truing solver where
    /// applicable, while rejecting NaN/infinity at the library boundary.
    pub fn validate(&self) -> Result<(), String> {
        if !self.muzzle_velocity_fps.is_finite()
            || !(TRUING_MV_MIN_FPS..=TRUING_MV_MAX_FPS).contains(&self.muzzle_velocity_fps)
        {
            return Err(format!(
                "nominal muzzle velocity must be finite and within {TRUING_MV_MIN_FPS:.0}..={TRUING_MV_MAX_FPS:.0} fps"
            ));
        }
        if !self.ballistic_coefficient.is_finite()
            || !(TRUING_BC_MIN..=TRUING_BC_MAX).contains(&self.ballistic_coefficient)
        {
            return Err(format!(
                "nominal ballistic coefficient must be finite and within {TRUING_BC_MIN:.2}..={TRUING_BC_MAX:.1}"
            ));
        }
        for (name, value) in [
            ("bullet mass", self.mass_gr),
            ("bullet diameter", self.diameter_in),
            ("zero distance", self.zero_distance_yd),
            ("sight height", self.sight_height_in),
            ("pressure", self.pressure_inhg),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(format!("{name} must be positive and finite"));
            }
        }
        if !self.temperature_f.is_finite() {
            return Err("temperature must be finite".to_string());
        }
        if !self.humidity_pct.is_finite() || !(0.0..=100.0).contains(&self.humidity_pct) {
            return Err("humidity must be finite and within 0..=100 percent".to_string());
        }
        if !self.altitude_ft.is_finite() {
            return Err("altitude must be finite".to_string());
        }
        Ok(())
    }

    /// Borrow this owned public shape as the existing internal forward model.
    /// The local `None` is intentional: V1 estimates one scalar BC and therefore
    /// cannot accept a velocity-banded BC schedule.
    pub(crate) fn with_forward_model<T>(
        &self,
        drop_unit: DropUnit,
        f: impl FnOnce(&TruingForwardModel<'_>) -> T,
    ) -> T {
        let no_bc_segments = None;
        let model = TruingForwardModel {
            drag_model: self.drag_model,
            mass_gr: self.mass_gr,
            diameter_in: self.diameter_in,
            zero_yd: self.zero_distance_yd,
            sight_in: self.sight_height_in,
            temp_f: self.temperature_f,
            press_inhg: self.pressure_inhg,
            humidity: self.humidity_pct,
            alt_ft: self.altitude_ft,
            bc_segments: &no_bc_segments,
            drop_unit,
        };
        f(&model)
    }
}

/// Parse an `--observed RANGE:DROP` token. RANGE is in the caller's distance
/// units (yards imperial / meters metric) and is normalized to internal yards;
/// DROP stays in the caller's drop unit. Returns a user-facing error string on
/// malformed input so the CLI can report cleanly instead of panicking.
pub fn parse_truing_observation(s: &str, units: UnitSystem) -> Result<TruingObservation, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return Err(format!(
            "invalid --observed '{s}': expected RANGE:DROP (e.g. 600:5.1)"
        ));
    }
    let range: f64 = parts[0]
        .trim()
        .parse()
        .map_err(|_| format!("invalid --observed range '{}' in '{s}'", parts[0]))?;
    let drop: f64 = parts[1]
        .trim()
        .parse()
        .map_err(|_| format!("invalid --observed drop '{}' in '{s}'", parts[1]))?;
    if !range.is_finite() || !drop.is_finite() {
        return Err(format!("invalid --observed '{s}': values must be finite"));
    }
    let range_yd = match units {
        UnitSystem::Imperial => range,
        UnitSystem::Metric => range / 0.9144,
    };
    Ok(TruingObservation { range_yd, drop })
}

/// Validate a truing observation set: every range finite and positive, every
/// drop finite and non-zero, and no two observations at (numerically) the same
/// range. [`run_multi_observation_truing_core`] runs this itself, so callers
/// don't have to — it is public so a front end can pre-validate when it wants
/// to sequence its own progress output strictly after validation (the native
/// CLI prints its "Fitting N observations..." progress line only for sets that
/// will actually be fitted). Error strings are the stable user-facing ones.
pub fn validate_truing_observations(observations: &[TruingObservation]) -> Result<(), String> {
    for o in observations {
        if !o.range_yd.is_finite() || o.range_yd <= 0.0 {
            return Err(format!(
                "observation range must be a positive finite distance (got {})",
                o.range_yd
            ));
        }
        if !o.drop.is_finite() || o.drop == 0.0 {
            return Err(
                "observation drop must be non-zero (a zero drop carries no truing information)"
                    .to_string(),
            );
        }
    }
    for i in 0..observations.len() {
        for j in (i + 1)..observations.len() {
            if (observations[i].range_yd - observations[j].range_yd).abs() < 1e-6 {
                return Err(format!(
                    "duplicate observation range ({:.3} yd internal): each observation must be at a distinct range",
                    observations[i].range_yd
                ));
            }
        }
    }
    Ok(())
}

/// Fixed load / atmosphere for a truing fit. The two free parameters (muzzle
/// velocity and BC) are supplied per prediction so the fitter can vary them.
pub(crate) struct TruingForwardModel<'a> {
    pub drag_model: DragModelArg,
    pub mass_gr: f64,
    pub diameter_in: f64,
    pub zero_yd: f64,
    pub sight_in: f64,
    pub temp_f: f64,
    pub press_inhg: f64,
    pub humidity: f64,
    pub alt_ft: f64,
    pub bc_segments: &'a Option<Vec<BCSegmentData>>,
    pub drop_unit: DropUnit,
}

impl TruingForwardModel<'_> {
    /// Predicted drop (in the configured drop unit) at `range_yd` for the given
    /// muzzle velocity and BC, using the real trajectory solver.
    pub fn predict(&self, mv_fps: f64, bc: f64, range_yd: f64) -> Result<f64, Box<dyn Error>> {
        self.predict_in_unit(mv_fps, bc, range_yd, self.drop_unit)
    }

    /// `predict` expressed in an explicit unit, independent of the configured
    /// `--drop-unit`. The identifiability diagnostics use this to stay in mil space
    /// (MBA-1337 t1): the linear `in` unit weights each Jacobian row by its range,
    /// which shifts the column correlation. NOTE this makes the gate DIAGNOSTICS
    /// unit-invariant; the MV-only operating point and the least-squares cost are
    /// still minimized in the display unit (deliberately out of scope — changing
    /// them changes fitted numbers), so extreme near-boundary sets can still differ.
    pub fn predict_in_unit(
        &self,
        mv_fps: f64,
        bc: f64,
        range_yd: f64,
        unit: DropUnit,
    ) -> Result<f64, Box<dyn Error>> {
        let (drop_m, z_m) = solve_trajectory_drop(
            mv_fps,
            bc,
            self.drag_model,
            self.mass_gr,
            self.diameter_in,
            self.zero_yd,
            range_yd,
            self.sight_in,
            self.temp_f,
            self.press_inhg,
            self.humidity,
            self.alt_ft,
            self.bc_segments,
            true, // interpolate: smooth forward model for the fitter
        )?;
        Ok(unit.express_drop_m(drop_m, z_m))
    }

    /// Predict several ranges from one zero solve and one trajectory integration.
    /// `None` marks a candidate the trajectory did not physically reach (or an
    /// invalid non-positive/non-finite candidate); other solver failures remain
    /// errors because they invalidate the entire nominal/perturbed model.
    ///
    /// This is the high-throughput path for experiment design and predictive
    /// bands.  It deliberately reproduces [`solve_trajectory_drop`]'s assembly
    /// and interpolation, but integrates through the farthest candidate once.
    pub(crate) fn predict_many_in_unit(
        &self,
        mv_fps: f64,
        bc: f64,
        ranges_yd: &[f64],
        unit: DropUnit,
    ) -> Result<Vec<Option<f64>>, Box<dyn Error>> {
        if ranges_yd.is_empty() {
            return Ok(Vec::new());
        }
        let max_range_yd = ranges_yd
            .iter()
            .copied()
            .filter(|r| r.is_finite() && *r > 0.0)
            .fold(0.0_f64, f64::max);
        if max_range_yd <= 0.0 {
            return Ok(vec![None; ranges_yd.len()]);
        }

        let velocity_ms = mv_fps * 0.3048;
        let mass_kg = self.mass_gr * GRAINS_TO_KG;
        let diameter_m = self.diameter_in * 0.0254;
        let zero_m = self.zero_yd * 0.9144;
        let max_range_m = max_range_yd * 0.9144;
        let sight_height_m = self.sight_in * 0.0254;
        let altitude_m = self.alt_ft * 0.3048;
        let temperature_c = (self.temp_f - 32.0) * 5.0 / 9.0;
        let pressure_hpa = self.press_inhg * 33.8639;
        let drag_model = match self.drag_model {
            DragModelArg::G1 => DragModel::G1,
            DragModelArg::G7 => DragModel::G7,
        };

        let mut inputs = BallisticInputs {
            muzzle_velocity: velocity_ms,
            bc_value: bc,
            bc_type: drag_model,
            bullet_mass: mass_kg,
            bullet_diameter: diameter_m,
            bullet_length: fallback_bullet_length_m(diameter_m, mass_kg),
            sight_height: sight_height_m,
            target_distance: max_range_m + 100.0,
            use_bc_segments: self.bc_segments.is_some(),
            bc_segments_data: self.bc_segments.clone(),
            use_rk4: true,
            muzzle_angle: 0.0,
            ..Default::default()
        };
        let atmosphere = AtmosphericConditions {
            temperature: temperature_c,
            pressure: pressure_hpa,
            humidity: self.humidity,
            altitude: altitude_m,
        };
        let wind = WindConditions::default();
        inputs.muzzle_angle = crate::calculate_zero_angle_with_conditions(
            inputs.clone(),
            zero_m,
            sight_height_m,
            wind.clone(),
            atmosphere.clone(),
        )?;

        let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
        solver.set_max_range(max_range_m + 100.0);
        solver.set_time_step(0.0001);
        let result = solver.solve()?;

        let predictions = ranges_yd
            .iter()
            .map(|range_yd| {
                if !range_yd.is_finite() || *range_yd <= 0.0 {
                    return None;
                }
                let range_m = *range_yd * 0.9144;
                let idx = result.points.iter().position(|p| p.position.x >= range_m)?;
                let (z_m, bullet_y) = if idx > 0 {
                    let p0 = &result.points[idx - 1];
                    let p1 = &result.points[idx];
                    let span = p1.position.x - p0.position.x;
                    if span.abs() < f64::EPSILON {
                        (p1.position.x, p1.position.y)
                    } else {
                        let fraction = (range_m - p0.position.x) / span;
                        (
                            range_m,
                            p0.position.y + fraction * (p1.position.y - p0.position.y),
                        )
                    }
                } else {
                    let p = &result.points[idx];
                    (p.position.x, p.position.y)
                };
                let drop_m = sight_height_m - bullet_y;
                Some(unit.express_drop_m(drop_m, z_m))
            })
            .collect();
        Ok(predictions)
    }

    /// Sum of squared residuals (predicted - observed) over all observations.
    pub fn cost(&self, mv: f64, bc: f64, obs: &[TruingObservation]) -> Result<f64, Box<dyn Error>> {
        let mut c = 0.0;
        for o in obs {
            let r = self.predict(mv, bc, o.range_yd)? - o.drop;
            c += r * r;
        }
        Ok(c)
    }
}

/// One finite-difference row of the truing forward model.  The planner and
/// uncertainty engine deliberately share this helper with the point fitter so
/// that experiment-design claims and posterior diagnostics describe the exact
/// same physics and perturbation sizes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TruingJacobianRow {
    pub predicted_drop: f64,
    pub d_drop_d_mv: f64,
    pub d_drop_d_bc: f64,
}

pub(crate) fn truing_jacobian_row(
    model: &TruingForwardModel<'_>,
    mv_fps: f64,
    bc: f64,
    range_yd: f64,
    unit: DropUnit,
) -> Result<TruingJacobianRow, Box<dyn Error>> {
    let hmv = (mv_fps * 1e-3).max(0.5);
    let hbc = (bc * 1e-3).max(1e-4);
    let predicted_drop = model.predict_in_unit(mv_fps, bc, range_yd, unit)?;
    let d_drop_d_mv = (model.predict_in_unit(mv_fps + hmv, bc, range_yd, unit)?
        - model.predict_in_unit(mv_fps - hmv, bc, range_yd, unit)?)
        / (2.0 * hmv);
    let d_drop_d_bc = (model.predict_in_unit(mv_fps, bc + hbc, range_yd, unit)?
        - model.predict_in_unit(mv_fps, bc - hbc, range_yd, unit)?)
        / (2.0 * hbc);
    Ok(TruingJacobianRow {
        predicted_drop,
        d_drop_d_mv,
        d_drop_d_bc,
    })
}

/// Batched form of [`truing_jacobian_row`].  Five integrations (nominal,
/// MV+/-, BC+/-) cover every requested range.  A row is `None` if any of those
/// trajectories cannot reach that candidate, preventing one-sided or silently
/// fabricated sensitivity estimates.
pub(crate) fn truing_jacobian_rows(
    model: &TruingForwardModel<'_>,
    mv_fps: f64,
    bc: f64,
    ranges_yd: &[f64],
    unit: DropUnit,
) -> Result<Vec<Option<TruingJacobianRow>>, Box<dyn Error>> {
    let hmv = (mv_fps * 1e-3).max(0.5);
    let hbc = (bc * 1e-3).max(1e-4);
    let nominal = model.predict_many_in_unit(mv_fps, bc, ranges_yd, unit)?;
    let mv_plus = model.predict_many_in_unit(mv_fps + hmv, bc, ranges_yd, unit)?;
    let mv_minus = model.predict_many_in_unit(mv_fps - hmv, bc, ranges_yd, unit)?;
    let bc_plus = model.predict_many_in_unit(mv_fps, bc + hbc, ranges_yd, unit)?;
    let bc_minus = model.predict_many_in_unit(mv_fps, bc - hbc, ranges_yd, unit)?;

    Ok((0..ranges_yd.len())
        .map(|i| {
            Some(TruingJacobianRow {
                predicted_drop: nominal[i]?,
                d_drop_d_mv: (mv_plus[i]? - mv_minus[i]?) / (2.0 * hmv),
                d_drop_d_bc: (bc_plus[i]? - bc_minus[i]?) / (2.0 * hbc),
            })
        })
        .collect())
}

/// One-parameter (muzzle velocity) least-squares fit with BC held fixed. Damped
/// Gauss-Newton with a central finite-difference derivative; robust because drop
/// is monotonic in muzzle velocity. Returns `(mv_fps, iterations, converged)`.
pub(crate) fn fit_truing_mv_only(
    model: &TruingForwardModel<'_>,
    obs: &[TruingObservation],
    bc: f64,
    mv_init: f64,
) -> Result<(f64, usize, bool), Box<dyn Error>> {
    let mut mv = mv_init.clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
    let mut converged = false;
    let mut iters = 0;
    for i in 0..TRUING_MAX_ITERS {
        iters = i + 1;
        let h = (mv * 1e-3).max(0.5);
        let mut num = 0.0;
        let mut den = 0.0;
        for o in obs {
            let r = model.predict(mv, bc, o.range_yd)? - o.drop;
            let dp = model.predict(mv + h, bc, o.range_yd)?;
            let dm = model.predict(mv - h, bc, o.range_yd)?;
            let j = (dp - dm) / (2.0 * h);
            num += j * r;
            den += j * j;
        }
        if den < 1e-12 {
            break;
        }
        // Gauss-Newton step, limited to keep the solver in a sane regime.
        let step = (-num / den).clamp(-300.0, 300.0);
        let new_mv = (mv + step).clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
        if (new_mv - mv).abs() < 0.05 {
            mv = new_mv;
            converged = true;
            break;
        }
        mv = new_mv;
    }
    Ok((mv, iters, converged))
}

/// Identifiability diagnostics for BC at the operating point `(mv, bc)`.
///
/// Returns `(sensitivity_ratio, condition_number)`:
///
/// * `sensitivity_ratio = ||bc*d(drop)/d(bc)|| / ||mv*d(drop)/d(mv)||` — the
///   relative influence of a fractional BC change vs a fractional MV change on
///   the predicted drops. Small => BC barely moves the drops (weak-signal mode).
///
/// * `condition_number = (1+|c|)/(1-|c|)` where `c` is the correlation between
///   the raw MV and BC Jacobian columns. Large => the two columns point the same
///   way, so a BC change can be undone by an MV change and the pair is not
///   separable (collinearity mode). This is magnitude-independent, unlike the
///   scaled-normal-matrix condition, so it actually tracks observation spread.
pub(crate) fn truing_identifiability(
    model: &TruingForwardModel<'_>,
    obs: &[TruingObservation],
    mv: f64,
    bc: f64,
) -> Result<(f64, f64), Box<dyn Error>> {
    let (mut n_mv, mut n_bc, mut cross) = (0.0, 0.0, 0.0);
    // Differentiate in mil space regardless of --drop-unit (MBA-1337 t1) so these
    // gate diagnostics do not shift with the display unit. (The operating point mv
    // comes from a display-unit fit, so full invariance is not guaranteed — see
    // predict_in_unit's note.)
    let unit = DropUnit::Mil;
    for o in obs {
        let row = truing_jacobian_row(model, mv, bc, o.range_yd, unit)?;
        n_mv += row.d_drop_d_mv * row.d_drop_d_mv;
        n_bc += row.d_drop_d_bc * row.d_drop_d_bc;
        cross += row.d_drop_d_mv * row.d_drop_d_bc;
    }
    let norm_mv = n_mv.sqrt();
    let norm_bc = n_bc.sqrt();
    let sensitivity_ratio = if mv * norm_mv > 0.0 {
        (bc * norm_bc) / (mv * norm_mv)
    } else {
        0.0
    };
    // Correlation between the raw Jacobian columns.
    let condition_number = if norm_mv > 0.0 && norm_bc > 0.0 {
        let c = (cross / (norm_mv * norm_bc)).clamp(-1.0, 1.0).abs();
        if (1.0 - c) > 1e-15 {
            (1.0 + c) / (1.0 - c)
        } else {
            f64::INFINITY
        }
    } else {
        // One column is numerically zero: the missing parameter is unconstrained.
        f64::INFINITY
    };
    Ok((sensitivity_ratio, condition_number))
}

/// Two-parameter (muzzle velocity, BC) joint fit via Levenberg-Marquardt
/// (damped Gauss-Newton) with central finite-difference Jacobian. Only the
/// diagonal of the normal matrix is damped (classic Marquardt scaling). Returns
/// `(mv_fps, bc, iterations, converged)`.
pub(crate) fn fit_truing_joint(
    model: &TruingForwardModel<'_>,
    obs: &[TruingObservation],
    mv_init: f64,
    bc_init: f64,
) -> Result<(f64, f64, usize, bool), Box<dyn Error>> {
    let mut mv = mv_init.clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
    let mut bc = bc_init.clamp(TRUING_BC_MIN, TRUING_BC_MAX);
    let mut lambda = 1e-3;
    let mut cur_cost = model.cost(mv, bc, obs)?;
    let mut converged = false;
    let mut iters = 0;

    for it in 0..TRUING_MAX_ITERS {
        iters = it + 1;
        // Build J^T J (a00,a01,a11) and J^T r (g0,g1).
        let (mut a00, mut a01, mut a11) = (0.0, 0.0, 0.0);
        let (mut g0, mut g1) = (0.0, 0.0);
        for o in obs {
            let row = truing_jacobian_row(model, mv, bc, o.range_yd, model.drop_unit)?;
            let r = row.predicted_drop - o.drop;
            let jmv = row.d_drop_d_mv;
            let jbc = row.d_drop_d_bc;
            a00 += jmv * jmv;
            a01 += jmv * jbc;
            a11 += jbc * jbc;
            g0 += jmv * r;
            g1 += jbc * r;
        }

        // Inner loop: grow lambda until a step reduces the cost.
        let mut accepted = false;
        for _ in 0..30 {
            let m00 = a00 + lambda * a00.max(1e-12);
            let m11 = a11 + lambda * a11.max(1e-12);
            let det = m00 * m11 - a01 * a01;
            if det.abs() < 1e-20 {
                lambda *= 10.0;
                continue;
            }
            // delta = -(JtJ + lambda*diag)^-1 * Jtr
            let dmv = -(m11 * g0 - a01 * g1) / det;
            let dbc = -(-a01 * g0 + m00 * g1) / det;
            let nmv = (mv + dmv).clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
            let nbc = (bc + dbc).clamp(TRUING_BC_MIN, TRUING_BC_MAX);
            let nc = model.cost(nmv, nbc, obs)?;
            if nc < cur_cost {
                let rel_change =
                    (nmv - mv).abs() / mv.max(1.0) + (nbc - bc).abs() / bc.max(1e-3);
                mv = nmv;
                bc = nbc;
                cur_cost = nc;
                lambda = (lambda * 0.5).max(1e-9);
                accepted = true;
                if rel_change < 1e-6 {
                    converged = true;
                }
                break;
            }
            lambda *= 4.0;
            if lambda > 1e12 {
                break;
            }
        }
        if !accepted {
            // No downhill step exists within damping range: we are at (or
            // numerically indistinguishable from) a local optimum.
            converged = true;
            break;
        }
        if converged {
            break;
        }
    }
    Ok((mv, bc, iters, converged))
}

/// Final report produced by a multi-observation truing fit.
#[derive(Debug, Clone)]
pub struct MultiTruingReport {
    pub fitted_mv_fps: f64,
    pub fitted_bc: f64,
    pub bc_input: f64,
    pub bc_fitted: bool,
    pub observations: Vec<TruingObservation>,
    pub predicted: Vec<f64>,
    pub residuals: Vec<f64>,
    pub rms: f64,
    pub iterations: usize,
    pub converged: bool,
    pub sensitivity_ratio: f64,
    pub condition_number: f64,
    pub quality: String,
    pub reason: String,
    /// Downward Mach-1.2 crossing distance (meters) for the finally fitted
    /// (`fitted_mv_fps`, `fitted_bc`) load, re-solved out to a fixed 3000 yd
    /// envelope independent of the observation set (the window describes the
    /// load, not wherever the caller happened to shoot). Feeds
    /// [`mv_calibration_window`] for the MV-calibration window report line
    /// (MBA-1405 Task 2). `None` when the trajectory never crosses within that
    /// envelope.
    pub mach_1_2_distance_m: Option<f64>,
    /// The downrange distance (meters) the fixed-envelope re-solve actually
    /// reached — names the range checked in the "no MV window" report note
    /// when `mach_1_2_distance_m` is `None` (MBA-1405 Task 2).
    pub window_solved_range_m: f64,
    /// Mach number at the first trajectory point of the fixed-envelope re-solve
    /// (`points[0].velocity_magnitude / station_speed_of_sound_mps`, the same
    /// convention `apply_dsf`/`truing_dsf` use for per-point Mach). Distinguishes,
    /// when `mach_1_2_distance_m` is `None`, a load that is still supersonic at
    /// the end of the window envelope (muzzle Mach >= 1.2) from one that launches
    /// below Mach 1.2 and never crosses downward at all (muzzle Mach < 1.2) — the
    /// two have opposite report text (MBA-1405 Task 2 fix pass).
    pub muzzle_mach: f64,
}

/// Orchestrate the multi-observation joint MV+BC calibration and build the
/// final [`MultiTruingReport`] (MBA-1343: the compute half only — rendering
/// stays with the caller).
///
/// `observations` is the already-parsed observation set, primary
/// (`--range`/`--measured-drop`) first, then each `--observed` impact in order
/// (use [`parse_truing_observation`] to build them from `RANGE:DROP` tokens);
/// every drop is already in `drop_unit` and every range in internal yards. The
/// set is validated with [`validate_truing_observations`] before any fitting.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable true-velocity CLI command shape"
)]
pub fn run_multi_observation_truing_core(
    observations: &[TruingObservation],
    drop_unit: DropUnit,
    bc_input: f64,
    drag_model: DragModelArg,
    mass_gr: f64,
    diameter_in: f64,
    zero_yd: f64,
    sight_in: f64,
    temp_f: f64,
    press_inhg: f64,
    humidity: f64,
    alt_ft: f64,
    bc_segments: &Option<Vec<BCSegmentData>>,
) -> Result<MultiTruingReport, Box<dyn Error>> {
    // Validate: finite, positive range, non-zero drop, no duplicate ranges.
    validate_truing_observations(observations)?;
    let observations: Vec<TruingObservation> = observations.to_vec();

    let model = TruingForwardModel {
        drag_model,
        mass_gr,
        diameter_in,
        zero_yd,
        sight_in,
        temp_f,
        press_inhg,
        humidity,
        alt_ft,
        bc_segments,
        drop_unit,
    };

    // Step 1: MV-only fit holding BC at the supplied value. This is always
    // well-posed and gives a good operating point for the identifiability check
    // and (if BC is identifiable) the joint fit.
    let mv_init = (TRUING_MV_MIN_FPS + TRUING_MV_MAX_FPS) / 2.0;
    let (mv0, mv_iters, mv_conv) = fit_truing_mv_only(&model, &observations, bc_input, mv_init)?;
    let rms_mv_only = rms_at(&model, &observations, mv0, bc_input)?;

    // Step 2: is BC identifiable from this observation set?
    let (sensitivity_ratio, condition_number) =
        truing_identifiability(&model, &observations, mv0, bc_input)?;
    let bc_identifiable = sensitivity_ratio >= TRUING_MIN_BC_SENSITIVITY_RATIO
        && condition_number <= TRUING_MAX_CONDITION_NUMBER
        && condition_number.is_finite();

    // Step 3: joint fit when identifiable, with a guard against a worse or
    // out-of-bounds result (never report a garbage joint fit).
    let mut fitted_mv = mv0;
    let mut fitted_bc = bc_input;
    let mut bc_fitted = false;
    let mut iterations = mv_iters;
    // Start from the MV-only fitter's own flag (MBA-1337 t2): both MV-only outcomes
    // (gate-refused joint, or joint rejected as worse) report THIS fit's convergence.
    // The accepted-joint branch overwrites it with the joint fitter's flag.
    let mut converged = mv_conv;
    let mut reason = String::new();

    if bc_identifiable {
        let (mv_j, bc_j, iters_j, conv_j) =
            fit_truing_joint(&model, &observations, mv0, bc_input)?;
        let rms_joint = rms_at(&model, &observations, mv_j, bc_j)?;
        let bc_at_bound = bc_j <= TRUING_BC_MIN * 1.001 || bc_j >= TRUING_BC_MAX * 0.999;
        if !bc_at_bound && rms_joint <= rms_mv_only + 1e-9 {
            fitted_mv = mv_j;
            fitted_bc = bc_j;
            bc_fitted = true;
            iterations = iters_j;
            converged = conv_j;
        } else {
            // Joint fit did not help (or ran to a bound): keep the honest
            // MV-only answer rather than a false-precision BC.
            reason = if bc_at_bound {
                format!(
                    "joint fit drove BC to a bound ({bc_j:.3}); BC held at input {bc_input:.3}"
                )
            } else {
                format!(
                    "joint fit did not improve on the MV-only solution; BC held at input {bc_input:.3}"
                )
            };
        }
    } else {
        reason = if !condition_number.is_finite() || condition_number > TRUING_MAX_CONDITION_NUMBER
        {
            format!(
                "observation ranges are too similar to separate MV from BC (condition {condition_number:.3e} > {TRUING_MAX_CONDITION_NUMBER:.0e}); BC held at input {bc_input:.3}"
            )
        } else {
            format!(
                "observations do not constrain BC (BC sensitivity ratio {sensitivity_ratio:.4} < {TRUING_MIN_BC_SENSITIVITY_RATIO:.2} threshold); BC held at input {bc_input:.3}. Add a longer-range / transonic observation to fit BC."
            )
        };
    }

    // Final residuals at the reported parameters. Alongside the display-unit RMS,
    // accumulate a mil-equivalent RMS: the quality bands were calibrated in mil
    // (~0.03 mil observation noise), so banding must not shift with --drop-unit
    // (MBA-1337 t1). Reported numbers stay in the user's unit.
    let mut predicted = Vec::with_capacity(observations.len());
    let mut residuals = Vec::with_capacity(observations.len());
    let mut sse = 0.0;
    let mut sse_mil = 0.0;
    for o in &observations {
        let p = model.predict(fitted_mv, fitted_bc, o.range_yd)?;
        let r = p - o.drop;
        let r_mil = match drop_unit {
            DropUnit::Mil => r,
            // moa -> mil: divide by (180/pi)*60/1000 moa-per-mil.
            DropUnit::Moa => r / ((180.0 / std::f64::consts::PI) * 60.0 / 1000.0),
            // inches -> meters -> small-angle mil at this observation's range.
            DropUnit::In => r * 0.0254 / (o.range_yd * 0.9144) * 1000.0,
        };
        predicted.push(p);
        residuals.push(r);
        sse += r * r;
        sse_mil += r_mil * r_mil;
    }
    let rms = (sse / observations.len() as f64).sqrt();
    let rms_mil = (sse_mil / observations.len() as f64).sqrt();

    let quality = truing_quality_line(
        bc_fitted,
        rms,
        rms_mil,
        drop_unit,
        condition_number,
        converged,
        observations.len(),
    );

    // MBA-1405 Task 2: one extra solve at the finally fitted (mv, bc), out to a
    // fixed generous envelope independent of the observation set, purely to read
    // off the downward Mach-1.2 crossing for the MV-calibration window report line.
    let (mach_1_2_distance_m, window_solved_range_m, muzzle_mach) =
        truing_mv_window_mach_1_2_crossing(&model, fitted_mv, fitted_bc)?;

    let report = MultiTruingReport {
        fitted_mv_fps: fitted_mv,
        fitted_bc,
        bc_input,
        bc_fitted,
        observations,
        predicted,
        residuals,
        rms,
        iterations,
        converged,
        sensitivity_ratio,
        condition_number,
        quality,
        reason,
        mach_1_2_distance_m,
        window_solved_range_m,
        muzzle_mach,
    };

    Ok(report)
}

/// Downrange envelope (yards) for the extra solve [`truing_mv_window_mach_1_2_crossing`]
/// runs to locate the downward Mach-1.2 crossing for the MV-calibration window report
/// line (MBA-1405 Task 2): generous enough to cover essentially every practical
/// small-arms transonic transition, independent of wherever the caller's own
/// observations happen to sit — the window describes the load, not the observation set.
const MV_WINDOW_SOLVE_MAX_RANGE_YD: f64 = 3000.0;

/// Re-solve the finally fitted `(mv, bc)` load out to [`MV_WINDOW_SOLVE_MAX_RANGE_YD`]
/// to read off the downward Mach-1.2 crossing distance (meters) for the MV-calibration
/// window (MBA-1405 Task 2). Returns `(mach_1_2_distance_m, solved_max_range_m,
/// muzzle_mach)`; the second element lets the report's "no MV window" note name the
/// range actually checked rather than the nominal envelope (relevant if the solve
/// terminates early, e.g. ground impact). The third element is the Mach number at the
/// first trajectory point (`points[0].velocity_magnitude / station_speed_of_sound_mps`,
/// the same convention `truing_dsf::apply_dsf` uses) — it lets the caller distinguish,
/// when `mach_1_2_distance_m` is `None`, "still supersonic through the whole envelope"
/// (muzzle Mach >= 1.2) from "never reaches Mach 1.2 at all" (muzzle Mach < 1.2), which
/// need opposite report text (MBA-1405 Task 2 fix pass).
fn truing_mv_window_mach_1_2_crossing(
    model: &TruingForwardModel<'_>,
    mv_fps: f64,
    bc: f64,
) -> Result<(Option<f64>, f64, f64), Box<dyn Error>> {
    let result = solve_truing_trajectory(
        mv_fps,
        bc,
        model.drag_model,
        model.mass_gr,
        model.diameter_in,
        model.zero_yd,
        MV_WINDOW_SOLVE_MAX_RANGE_YD,
        model.sight_in,
        model.temp_f,
        model.press_inhg,
        model.humidity,
        model.alt_ft,
        model.bc_segments,
    )?;
    let muzzle_mach = result.points.first().map_or(0.0, |p| {
        if result.station_speed_of_sound_mps > 0.0 {
            p.velocity_magnitude / result.station_speed_of_sound_mps
        } else {
            0.0
        }
    });
    Ok((result.mach_1_2_distance_m, result.max_range, muzzle_mach))
}

/// RMS of residuals at a candidate `(mv, bc)`.
pub(crate) fn rms_at(
    model: &TruingForwardModel<'_>,
    obs: &[TruingObservation],
    mv: f64,
    bc: f64,
) -> Result<f64, Box<dyn Error>> {
    let mut sse = 0.0;
    for o in obs {
        let r = model.predict(mv, bc, o.range_yd)? - o.drop;
        sse += r * r;
    }
    Ok((sse / obs.len() as f64).sqrt())
}

/// Plain-language quality assessment for the fit. `rms` is in the user's drop unit
/// (displayed); `rms_mil` is the mil-equivalent the bands were calibrated against
/// (~0.03 mil observation noise), so the quality word is unit-invariant (MBA-1337 t1).
pub(crate) fn truing_quality_line(
    bc_fitted: bool,
    rms: f64,
    rms_mil: f64,
    drop_unit: DropUnit,
    condition_number: f64,
    converged: bool,
    n_obs: usize,
) -> String {
    let unit = drop_unit.label();
    let n_params = if bc_fitted { 2 } else { 1 };
    // Exactly-determined fit (MBA-1337 t3): zero degrees of freedom drive the
    // residuals to ~0 by construction — an "excellent" RMS validates nothing.
    if n_obs == n_params {
        return format!(
            "{} fit is exactly determined ({n_obs} observations, {n_params} fitted \
             parameters): residuals are zero by construction and do not validate the \
             fit; add an observation to assess quality",
            if bc_fitted { "Joint MV+BC" } else { "MV-only" }
        );
    }
    let quality = if rms_mil < 0.05 {
        "excellent"
    } else if rms_mil < 0.15 {
        "good"
    } else if rms_mil < 0.4 {
        "fair"
    } else {
        "poor (observations may be inconsistent)"
    };
    let nonconv = if converged { "" } else { " (did not fully converge)" };
    if bc_fitted {
        let cond = if condition_number.is_finite() {
            format!("{condition_number:.0}")
        } else {
            "inf".to_string()
        };
        format!(
            "Joint MV+BC fit, {quality}: RMS residual {rms:.3} {unit}, conditioning {cond}{nonconv}"
        )
    } else {
        format!("MV-only fit, {quality}: RMS residual {rms:.3} {unit} (BC held fixed){nonconv}")
    }
}

/// Fraction of the downward Mach 1.2 crossing distance that bounds the near edge of the MV
/// (muzzle-velocity) calibration window (MBA-1405): observations closer than this fraction of
/// the crossing are considered too far from the transonic region to usefully calibrate MV.
const MV_CALIBRATION_WINDOW_START_FRACTION: f64 = 0.90;

/// Fraction of the downward Mach 0.9 crossing distance at which the DSF (drag-scale-factor)
/// window is considered to start (MBA-1405): DSF observations nearer the muzzle than this are
/// still comfortably supersonic/transonic and not representative of the DSF-correction regime.
const DSF_WINDOW_START_FRACTION: f64 = 0.90;

/// The MV (muzzle-velocity) calibration validity window: `(start_m, end_m)`, where `end_m` is
/// the downward Mach 1.2 crossing distance and `start_m` is 90% of it. `None` when the
/// trajectory never crosses Mach 1.2 (nothing to bound the window with) — mirrors
/// `TrajectoryResult::mach_1_2_distance_m`'s own `None` case, MBA-1405. `pub` (not
/// `pub(crate)`): consumed by the native CLI (`main.rs`) and the WASM terminal
/// (`wasm.rs`), both separate crates from this lib (MBA-1405 Task 2).
pub fn mv_calibration_window(mach_1_2_distance_m: Option<f64>) -> Option<(f64, f64)> {
    let d = mach_1_2_distance_m?;
    Some((MV_CALIBRATION_WINDOW_START_FRACTION * d, d))
}

/// The downrange distance (m) at which the DSF (drag-scale-factor) window starts: 90% of the
/// downward Mach 0.9 crossing distance. `None` when the trajectory never crosses Mach 0.9 —
/// mirrors `TrajectoryResult::mach_0_9_distance_m`'s own `None` case, MBA-1405. `pub`: consumed
/// by the native CLI's `dsf` verb (MBA-1405 Task 2).
pub fn dsf_window_start(mach_0_9_distance_m: Option<f64>) -> Option<f64> {
    Some(DSF_WINDOW_START_FRACTION * mach_0_9_distance_m?)
}

#[cfg(test)]
mod window_helper_tests {
    use super::*;

    #[test]
    fn mv_calibration_window_is_90_to_100_percent_of_the_1_2_crossing() {
        assert_eq!(mv_calibration_window(Some(1000.0)), Some((900.0, 1000.0)));
        assert_eq!(mv_calibration_window(Some(671.7257336844475)), Some((0.9 * 671.7257336844475, 671.7257336844475)));
    }

    #[test]
    fn mv_calibration_window_passes_through_none() {
        assert_eq!(mv_calibration_window(None), None);
    }

    #[test]
    fn dsf_window_start_is_90_percent_of_the_0_9_crossing() {
        assert_eq!(dsf_window_start(Some(1000.0)), Some(900.0));
        assert_eq!(dsf_window_start(Some(806.5709746782849)), Some(0.9 * 806.5709746782849));
    }

    #[test]
    fn dsf_window_start_passes_through_none() {
        assert_eq!(dsf_window_start(None), None);
    }

    #[test]
    fn mv_calibration_window_exact_arithmetic_at_a_representative_distance() {
        // Exact f64 arithmetic check (not just structural equality): 0.90 * 500.0 = 450.0
        // has an exact binary representation, so this pins the literal formula, not just
        // its shape.
        let (start, end) = mv_calibration_window(Some(500.0)).expect("Some");
        assert_eq!(start.to_bits(), 450.0_f64.to_bits());
        assert_eq!(end.to_bits(), 500.0_f64.to_bits());
    }

    #[test]
    fn dsf_window_start_exact_arithmetic_at_a_representative_distance() {
        let start = dsf_window_start(Some(500.0)).expect("Some");
        assert_eq!(start.to_bits(), 450.0_f64.to_bits());
    }
}

/// MBA-1377: `correct_chrono_velocity_fps` (screen-distance chronograph correction).
#[cfg(test)]
mod chrono_correction_tests {
    use super::*;

    /// Hand-computed pin, matching a real chronograph shot: 168 gr, 0.308", BC 0.475 G7, a
    /// reading taken 15 ft downrange, dry standard atmosphere (59 F / 29.92124356615747 inHg
    /// == 1013.25 hPa exactly / 0% humidity / sea level).
    ///
    /// Arithmetic (all from published/documented physical constants, not this function).
    /// Speed of sound (dry air, ICAO): `a = sqrt(gamma * R_air * T_K)`, gamma = 1.4, R_air =
    /// 287.0531 J/(kg K) (see `atmosphere.rs`), T_K = 288.15 (59 F == 15 C exactly), giving
    /// `a = sqrt(1.4 * 287.0531 * 288.15) = 340.294124 m/s = 1116.450575 fps`.
    ///
    /// The measured (screen) velocity is chosen so its Mach number lands exactly on a G7
    /// table grid point -- Mach 2.40, Cd = 0.2752 (data/g7.csv row "2.40,0.2752") -- which
    /// sidesteps interpolation uncertainty entirely (cubic Hermite interpolation is exact at
    /// its own knots): `v_screen = 2.40 * 1116.450575 = 2679.481380 fps`.
    ///
    /// Air density ratio vs the 1.225 kg/m^3 CD_TO_RETARD reference, via the crate's own
    /// published CIPM-2007 formula (`calculate_air_density_cimp`, dry air, sea level, 1013.25
    /// hPa, 15 C): `rho = 1.225521 kg/m^3`, `ratio = rho / 1.225 = 1.00042559`.
    ///
    /// Deceleration (McCoy retardation form, see the `constants::CD_TO_RETARD` doc comment):
    /// `a_ft_s2 = Cd * v_fps^2 * (rho/1.225) * CD_TO_RETARD / BC`, CD_TO_RETARD = 2.08551e-4.
    /// Over a flight this short, Cd and the density ratio are effectively constant, so
    /// `dv/dx = a/v = -(CD_TO_RETARD * Cd * ratio / BC) * v` (linear in v), giving the
    /// classic short-range exponential-decay solution `v(x) = v0 * exp(-C*x)`, x in FEET:
    /// `C = CD_TO_RETARD * Cd * ratio / BC = 2.08551e-4 * 0.2752 * 1.00042559 / 0.475
    /// = 1.2087929e-4` (1/ft), so `v0 = v_screen * exp(C * 15)
    /// = 2679.481380 * exp(1.2087929e-4 * 15) = 2684.344194 fps`.
    ///
    /// The real solver differs only by the (second-order, and here tiny) fact that its Cd is
    /// evaluated continuously at the LOCAL velocity along the RK45 path rather than held fixed
    /// at v_screen's value; the measured residual is ~0.0024 fps, so 0.05 fps is a safe,
    /// still-tight pin tolerance (0.001% of muzzle velocity).
    #[test]
    fn hand_computed_pin_168gr_308_bc475_g7_15ft() {
        let measured_screen_velocity_fps = 2679.481379882625;
        let correction = correct_chrono_velocity_fps(
            measured_screen_velocity_fps,
            15.0 * 0.3048, // 15 ft screen distance, in meters
            0.475,         // BC
            DragModelArg::G7,
            168.0,               // mass, grains
            0.308,               // diameter, inches
            59.0,                // temperature, F (== 15 C exactly)
            1013.25 / 33.8639,   // pressure, inHg (== 1013.25 hPa exactly)
            0.0,                 // humidity, % (dry air keeps the arithmetic exact)
            0.0,                 // altitude, ft
            &None,
        )
        .expect("hand-computed case must converge");

        let expected_v0 = 2684.344194108996;
        assert!(
            (correction.muzzle_velocity_fps - expected_v0).abs() < 0.05,
            "got {}, expected {} +/- 0.05",
            correction.muzzle_velocity_fps,
            expected_v0
        );
        // The corrected muzzle velocity must exceed the raw downrange reading -- drag only
        // decelerates, so a chronograph reading downrange always understates true MV.
        assert!(correction.muzzle_velocity_fps > measured_screen_velocity_fps);
        // McCoy/JBM: this band converges in a handful of secant iterations, not dozens.
        assert!(
            correction.iterations <= 6,
            "expected fast convergence, took {} iterations",
            correction.iterations
        );
    }

    /// Realistic smoke case from the task brief: 168 gr .308, BC 0.475 G7, 2680 fps measured
    /// at 15 ft, default atmosphere. The corrected muzzle velocity must be higher than the
    /// measured value by a small, physically plausible margin (a few fps, not tens).
    #[test]
    fn smoke_2680fps_at_15ft_corrects_upward_by_a_plausible_margin() {
        let correction = correct_chrono_velocity_fps(
            2680.0,
            15.0 * 0.3048,
            0.475,
            DragModelArg::G7,
            168.0,
            0.308,
            59.0,
            29.92,
            50.0,
            0.0,
            &None,
        )
        .expect("realistic case must converge");

        let margin = correction.muzzle_velocity_fps - 2680.0;
        assert!(
            (0.5..15.0).contains(&margin),
            "expected a plausible few-fps correction, got {margin} fps (V0={})",
            correction.muzzle_velocity_fps
        );
    }

    /// Convergence across the full stated screen-distance band (3-25 m): every distance in
    /// range solves without error, in a small number of secant iterations, and the correction
    /// grows monotonically with distance (more travel -> more velocity lost -> a bigger
    /// muzzle-velocity correction to explain the same downrange reading).
    #[test]
    fn convergence_across_the_3_to_25_m_band() {
        let mut last_correction = 0.0;
        for distance_m in [3.0, 5.0, 10.0, 15.0, 20.0, 25.0] {
            let correction = correct_chrono_velocity_fps(
                2680.0,
                distance_m,
                0.475,
                DragModelArg::G7,
                168.0,
                0.308,
                59.0,
                29.92,
                50.0,
                0.0,
                &None,
            )
            .unwrap_or_else(|e| panic!("distance {distance_m} m must converge: {e}"));
            assert!(
                correction.iterations <= 6,
                "distance {distance_m} m took {} iterations",
                correction.iterations
            );
            let delta = correction.muzzle_velocity_fps - 2680.0;
            assert!(
                delta > last_correction,
                "correction must grow with distance: {delta} <= {last_correction} at {distance_m} m"
            );
            last_correction = delta;
        }
    }

    /// Absent distance is handled entirely by the CLI/WASM callers (they special-case
    /// `None`/`0.0` as a byte-identical no-op and never call into this solver at all); at this
    /// layer, an explicit zero screen distance is simply out of the validated band and must be
    /// rejected rather than silently treated as "at the muzzle" by a solver that has no
    /// special-case for it.
    #[test]
    fn zero_distance_is_rejected_by_the_low_level_solver() {
        let err = correct_chrono_velocity_fps(
            2680.0, 0.0, 0.475, DragModelArg::G7, 168.0, 0.308, 59.0, 29.92, 50.0, 0.0, &None,
        )
        .unwrap_err();
        assert!(err.to_string().contains("out of range"), "{err}");
    }

    /// Out-of-band screen distances (negative, too-short, too-long) are rejected outright
    /// rather than silently clamped or extrapolated into a garbage muzzle velocity.
    #[test]
    fn out_of_band_distances_are_rejected() {
        for bad_distance_m in [-5.0, -0.001, 0.05, 30.001, 100.0, f64::NAN, f64::INFINITY] {
            let err = correct_chrono_velocity_fps(
                2680.0,
                bad_distance_m,
                0.475,
                DragModelArg::G7,
                168.0,
                0.308,
                59.0,
                29.92,
                50.0,
                0.0,
                &None,
            )
            .unwrap_err();
            assert!(
                err.to_string().contains("out of range"),
                "distance {bad_distance_m}: {err}"
            );
        }
    }

    /// A non-positive or non-finite measured velocity is rejected outright.
    #[test]
    fn invalid_measured_velocity_is_rejected() {
        for bad_v in [0.0, -100.0, f64::NAN, f64::INFINITY] {
            let err = correct_chrono_velocity_fps(
                bad_v, 4.572, 0.475, DragModelArg::G7, 168.0, 0.308, 59.0, 29.92, 50.0, 0.0, &None,
            )
            .unwrap_err();
            assert!(err.to_string().contains("must be positive and finite"), "{err}");
        }
    }

    /// The correction must actually consult the configured drag model, not a hardcoded G1: the
    /// SAME (measured velocity, distance, BC) triple fed through G1 and G7 must back-solve to
    /// DIFFERENT muzzle velocities, since the two standard drag curves are numerically distinct
    /// at this Mach range.
    #[test]
    fn drag_model_is_actually_consulted_g1_vs_g7_diverge() {
        let g1 = correct_chrono_velocity_fps(
            2680.0,
            15.0 * 0.3048,
            0.475,
            DragModelArg::G1,
            168.0,
            0.308,
            59.0,
            29.92,
            50.0,
            0.0,
            &None,
        )
        .expect("G1 must converge");
        let g7 = correct_chrono_velocity_fps(
            2680.0,
            15.0 * 0.3048,
            0.475,
            DragModelArg::G7,
            168.0,
            0.308,
            59.0,
            29.92,
            50.0,
            0.0,
            &None,
        )
        .expect("G7 must converge");

        assert!(
            (g1.muzzle_velocity_fps - g7.muzzle_velocity_fps).abs() > 0.05,
            "G1 ({}) and G7 ({}) corrections should differ measurably for the same inputs",
            g1.muzzle_velocity_fps,
            g7.muzzle_velocity_fps
        );
    }
}
