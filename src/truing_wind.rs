//! MBA-1392: back-solve the effective crosswind from an observed horizontal miss.
//!
//! The rest of the truing family fits muzzle velocity and BC from VERTICAL drop. This
//! module fits the other axis: given where a group actually landed left/right of the aim
//! point, it reports the constant crosswind that reproduces that miss through the real
//! forward model — the number a shooter compares against the wind they called, so their
//! wind calls can be calibrated instead of guessed at.
//!
//! # Sign conventions (the whole set, in one place)
//!
//! * **Observed miss** ([`WindObservation::miss_right_m`], CLI `--miss`): signed, POSITIVE
//!   = the group landed RIGHT of the aim point.
//! * **Solved crosswind** ([`WindTruingSolution::solved_crosswind_mph`]): signed, and the
//!   sign follows the deflection it produces — POSITIVE = a wind FROM the shooter's LEFT
//!   (9 o'clock) that pushes impacts RIGHT; NEGATIVE = a wind FROM the shooter's RIGHT
//!   (3 o'clock) pushing impacts LEFT. It is a full-value crosswind (a 90-degree wind),
//!   not a half-value component of some other bearing.
//! * That maps onto the engine's wind-FROM convention ([`crate::wind::wind_vector`],
//!   `0 = headwind`, `PI/2 = from the right`, `3*PI/2 = from the left`) — the convention
//!   established by the 0.19.0 wind-direction sign fix, which flipped `0` from tailwind
//!   to headwind. A positive solved crosswind is therefore direction `3*PI/2`.
//! * **Twist** ([`crate::truing::TruingTwist::right_hand`]): a right-hand twist drifts
//!   RIGHT (positive lateral), a left-hand twist drifts LEFT.
//! * **Shot azimuth** ([`crate::truing::TruingEarthFrame::shot_azimuth_deg`]): compass
//!   bearing fired ALONG, 0 = North, 90 = East.
//!
//! # What the solved wind actually contains
//!
//! A horizontal miss is not purely wind. Spin drift is always modelled here (a twist rate
//! is required, precisely so it can be), and Coriolis is modelled when a latitude and shot
//! azimuth are supplied. Anything the model was not given data for stays ABSORBED in the
//! solved crosswind, and the report says so ([`WindTruingReport::unsubtracted_effects`])
//! rather than quietly presenting a contaminated number as pure wind.
//!
//! # No scope-tracking correction
//!
//! `--miss` values are LINEAR measurements off the target (inches), not dial readings, so
//! the MBA-1358 tracking correction factor does NOT apply to them — only DIALED
//! observations are CF-converted. This module deliberately has no CF input.

use std::error::Error;

use serde::{Deserialize, Serialize};

use crate::cli_api::UnitSystem;
use crate::truing::{
    DragModelArg, TruingEarthFrame, TruingEnvironment, TruingTwist, TRUING_BC_MAX, TRUING_BC_MIN,
    TRUING_MV_MAX_FPS, TRUING_MV_MIN_FPS,
};
use crate::{BCSegmentData, WindConditions};

/// Miles per hour to meters per second (exact, by definition of the international mile).
/// `pub` so the front ends convert `--called-wind` with the same factor this module
/// renders with, rather than each re-typing the literal.
pub const MPH_TO_MPS: f64 = 0.44704;

/// Widest crosswind the solver will bracket, in mph (signed, so the bracket is
/// `-100..=+100`). Comfortably past any wind a rifle shooter reports; a miss that needs
/// more than this is a data-entry or sign error, not a wind call, and is rejected with a
/// diagnostic rather than solved into a fantasy number.
pub const MAX_SOLVABLE_CROSSWIND_MPH: f64 = 100.0;

/// Convergence tolerance on the modelled-minus-observed lateral miss, in meters
/// (0.01 mm — four orders of magnitude finer than anyone can measure a group centre).
/// `pub` because it is the published meaning of a converged
/// [`WindTruingSolution::residual_m`].
pub const WIND_SOLVE_TOLERANCE_M: f64 = 1.0e-5;

/// Bracket width (mph) at which the root find stops regardless of residual: the two ends
/// are numerically the same wind, so further bisection cannot improve the answer.
const WIND_SOLVE_MIN_BRACKET_MPH: f64 = 1.0e-9;

/// Iteration cap for the per-observation root find. Lateral deflection is monotone and
/// near-linear in crosswind speed (the Didion lag-time relation), so the bracketed
/// false-position iteration below normally converges in well under ten evaluations; this
/// is the runaway guard, not the expected count.
const WIND_SOLVE_MAX_ITERATIONS: u32 = 60;

/// Central-difference step (mph) used to measure how strongly the observation constrains
/// the wind. Large enough that the trajectory integrator's own sampling noise does not
/// dominate the difference, small enough to stay in the locally-linear regime.
const WIND_SENSITIVITY_STEP_MPH: f64 = 0.5;

/// Guide value for the wind-truing validity note: below this many inches of lateral
/// movement per mph of crosswind, the observation barely constrains the wind and the
/// number solved from it is weakly identified. Same spirit as the MV-calibration window
/// (MBA-1405) — a stated band the report checks each observation against, not a hard gate.
pub const MIN_WIND_SENSITIVITY_IN_PER_MPH: f64 = 0.25;

/// A single observed horizontal miss used to back-solve effective wind (MBA-1392).
///
/// SI throughout; front ends convert their display units once at the boundary (see
/// [`parse_wind_observation`]).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindObservation {
    /// Range at which the group centre was measured, meters.
    pub range_m: f64,
    /// Signed horizontal miss of the group centre, meters. POSITIVE = RIGHT of aim.
    pub miss_right_m: f64,
    /// Optional one-standard-deviation measurement error of `miss_right_m`, meters.
    /// Supply it on every observation or on none — see [`WindTruingRequest::validate`].
    pub sigma_m: Option<f64>,
}

/// Parse a `--miss RANGE:RIGHT[:SIGMA]` token.
///
/// `RANGE` is in the caller's distance units (yards imperial / meters metric); `RIGHT` and
/// `SIGMA` are LINEAR INCHES in both unit systems, matching the existing `--drop-unit in`
/// contract (a tape measurement off the target is reported in inches whatever the range
/// unit is). Returns a user-facing error string on malformed input.
pub fn parse_wind_observation(s: &str, units: UnitSystem) -> Result<WindObservation, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 && parts.len() != 3 {
        return Err(format!(
            "invalid --miss '{s}': expected RANGE:RIGHT_IN[:SIGMA] (e.g. 600:8.5 or 600:8.5:0.75)"
        ));
    }
    let range: f64 = parts[0]
        .trim()
        .parse()
        .map_err(|_| format!("invalid --miss range '{}' in '{s}'", parts[0]))?;
    let miss_in: f64 = parts[1]
        .trim()
        .parse()
        .map_err(|_| format!("invalid --miss offset '{}' in '{s}'", parts[1]))?;
    let sigma_in: Option<f64> = match parts.get(2) {
        Some(token) => Some(
            token
                .trim()
                .parse()
                .map_err(|_| format!("invalid --miss sigma '{token}' in '{s}'"))?,
        ),
        None => None,
    };
    if !range.is_finite() || !miss_in.is_finite() || sigma_in.is_some_and(|v| !v.is_finite()) {
        return Err(format!("invalid --miss '{s}': values must be finite"));
    }
    let range_m = match units {
        UnitSystem::Imperial => range * 0.9144,
        UnitSystem::Metric => range,
    };
    Ok(WindObservation {
        range_m,
        miss_right_m: miss_in * 0.0254,
        sigma_m: sigma_in.map(|v| v * 0.0254),
    })
}

/// Everything the wind fit needs: the observed misses plus the load, rifle, atmosphere and
/// the opt-in earth frame (MBA-1392).
///
/// Imperial units for the load/atmosphere fields, matching the truing core's historical
/// internal convention; the observations themselves are SI. Unlike the drop-based truing
/// commands this carries a KNOWN muzzle velocity — wind is the unknown being fitted, so
/// velocity is an input, not an output.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindTruingRequest {
    /// One or more observed horizontal misses. Ranges must be distinct.
    pub observations: Vec<WindObservation>,
    /// Known muzzle velocity, feet/second (true up first with `true-velocity` if unsure).
    pub muzzle_velocity_fps: f64,
    /// Scalar ballistic coefficient for `drag_model`.
    pub bc: f64,
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
    /// Barrel twist. REQUIRED: spin drift is a lateral effect of the same order as a
    /// light wind at long range, so without it the fit would silently report spin drift
    /// as wind.
    pub twist: TruingTwist,
    /// Latitude + shot azimuth. `None` leaves Coriolis unmodelled and absorbed into the
    /// solved wind (the report names it as unsubtracted).
    pub earth: Option<TruingEarthFrame>,
    /// The wind the shooter CALLED, mph, in the same signed convention as the solved
    /// value. `Some` adds a wind-call correction factor (solved / called) to the report.
    pub called_crosswind_mph: Option<f64>,
}

impl WindTruingRequest {
    /// Validate the whole request before any (expensive) trajectory work begins.
    ///
    /// Mirrors [`crate::truing::TruingModelInputsV1::validate`] plus the observation-set
    /// rules from [`crate::truing::validate_truing_observations`]: finite positive ranges,
    /// finite misses, no duplicate ranges. Sigmas are all-or-none — a half-weighted set
    /// would silently mix an inverse-variance mean with unit weights.
    pub fn validate(&self) -> Result<(), String> {
        if self.observations.is_empty() {
            return Err("at least one observed horizontal miss is required".to_string());
        }
        if !self.muzzle_velocity_fps.is_finite()
            || !(TRUING_MV_MIN_FPS..=TRUING_MV_MAX_FPS).contains(&self.muzzle_velocity_fps)
        {
            return Err(format!(
                "muzzle velocity must be finite and within {TRUING_MV_MIN_FPS:.0}..={TRUING_MV_MAX_FPS:.0} fps"
            ));
        }
        if !self.bc.is_finite() || !(TRUING_BC_MIN..=TRUING_BC_MAX).contains(&self.bc) {
            return Err(format!(
                "ballistic coefficient must be finite and within {TRUING_BC_MIN:.2}..={TRUING_BC_MAX:.1}"
            ));
        }
        for (name, value) in [
            ("bullet mass", self.mass_gr),
            ("bullet diameter", self.diameter_in),
            ("zero distance", self.zero_distance_yd),
            ("sight height", self.sight_height_in),
            ("pressure", self.pressure_inhg),
            ("twist rate", self.twist.rate_in),
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
        if let Some(earth) = self.earth {
            if !earth.latitude_deg.is_finite() || !(-90.0..=90.0).contains(&earth.latitude_deg) {
                return Err("latitude must be finite and within -90..=90 degrees".to_string());
            }
            if !earth.shot_azimuth_deg.is_finite() {
                return Err("shot azimuth must be finite".to_string());
            }
        }
        if let Some(called) = self.called_crosswind_mph {
            if !called.is_finite() || called == 0.0 {
                return Err(
                    "the called wind must be finite and non-zero (a zero call has no \
                     correction factor)"
                        .to_string(),
                );
            }
        }
        for observation in &self.observations {
            if !observation.range_m.is_finite() || observation.range_m <= 0.0 {
                return Err(format!(
                    "observation range must be a positive finite distance (got {})",
                    observation.range_m
                ));
            }
            if !observation.miss_right_m.is_finite() {
                return Err("observed horizontal miss must be finite".to_string());
            }
            if observation
                .sigma_m
                .is_some_and(|sigma| !sigma.is_finite() || sigma <= 0.0)
            {
                return Err("an observed-miss sigma must be positive and finite".to_string());
            }
        }
        for i in 0..self.observations.len() {
            for j in (i + 1)..self.observations.len() {
                if (self.observations[i].range_m - self.observations[j].range_m).abs() < 1e-6 {
                    return Err(format!(
                        "duplicate observation range ({:.3} m): each observed miss must be at a \
                         distinct range",
                        self.observations[i].range_m
                    ));
                }
            }
        }
        let with_sigma = self
            .observations
            .iter()
            .filter(|o| o.sigma_m.is_some())
            .count();
        if with_sigma != 0 && with_sigma != self.observations.len() {
            return Err(
                "supply a sigma on every observed miss or on none: mixing weighted and \
                 unweighted observations would silently combine inverse-variance weights \
                 with unit weights"
                    .to_string(),
            );
        }
        Ok(())
    }
}

/// The wind fitted from ONE observed miss (MBA-1392).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct WindTruingSolution {
    /// Observation range, meters.
    pub range_m: f64,
    /// The observed miss that was fitted, meters (positive = right).
    pub observed_miss_right_m: f64,
    /// The supplied measurement sigma, meters, if any.
    pub sigma_m: Option<f64>,
    /// The fitted constant crosswind, mph, signed (positive = from the left, pushing right).
    pub solved_crosswind_mph: f64,
    /// The model's lateral miss at `solved_crosswind_mph`, meters — equals the observed
    /// miss to within [`WIND_SOLVE_TOLERANCE_M`] on a converged solve.
    pub modeled_miss_right_m: f64,
    /// `modeled_miss_right_m - observed_miss_right_m`, meters.
    pub residual_m: f64,
    /// The model's lateral miss with ZERO wind, meters: the part of the observation
    /// attributed to spin drift (always) and Coriolis (when an earth frame was supplied)
    /// rather than to wind.
    pub no_wind_lateral_m: f64,
    /// How far the impact moves per mph of crosswind at the solution, meters/mph — the
    /// identifiability measure behind the report's weak-signal note.
    pub sensitivity_m_per_mph: f64,
    /// The observation sigma propagated into wind units, mph
    /// (`sigma_m / |sensitivity_m_per_mph|`); `None` when no sigma was supplied.
    pub solved_sigma_mph: Option<f64>,
    /// Root-find iterations actually run.
    pub iterations: u32,
    /// Whether the root find hit its tolerance (`false` = the reported value is the best
    /// estimate at the iteration cap).
    pub converged: bool,
}

/// The complete wind-truing result: one fit per observation plus the combined answer.
#[derive(Debug, Clone, Serialize)]
pub struct WindTruingReport {
    /// Per-observation fits, in the order the observations were supplied.
    pub solutions: Vec<WindTruingSolution>,
    /// Combined effective crosswind, mph, signed.
    pub mean_crosswind_mph: f64,
    /// One-sigma uncertainty of `mean_crosswind_mph`, mph; `Some` only when every
    /// observation carried a sigma.
    pub mean_sigma_mph: Option<f64>,
    /// `true` when the mean is inverse-variance weighted (all sigmas supplied), `false`
    /// when it is the plain arithmetic mean.
    pub inverse_variance_weighted: bool,
    /// The wind the shooter called, mph, if supplied.
    pub called_crosswind_mph: Option<f64>,
    /// `mean_crosswind_mph / called_crosswind_mph`: >1 means the shooter under-called the
    /// wind, <1 means they over-called it, negative means they called the wrong side.
    pub wind_call_factor: Option<f64>,
    /// Lateral effects the model actually accounted for, so they are NOT in the solved wind.
    pub subtracted_effects: Vec<String>,
    /// Lateral effects the model had no data for, which therefore ARE absorbed into the
    /// solved wind. Empty means everything this model knows about was subtracted.
    pub unsubtracted_effects: Vec<String>,
}

/// Wind conditions describing a signed full-value crosswind (MBA-1392).
///
/// Uses the engine's wind-FROM convention: `PI/2` is a wind FROM the shooter's RIGHT
/// (pushing impacts LEFT) and `3*PI/2` is FROM the LEFT (pushing impacts RIGHT). The
/// signed input follows the DEFLECTION, so positive maps to `3*PI/2`. The pair is
/// continuous through zero — both directions give the zero vector at zero speed — so the
/// root find sees a smooth function across the sign change.
fn crosswind_conditions(signed_mph: f64) -> WindConditions {
    WindConditions {
        speed: signed_mph.abs() * MPH_TO_MPS,
        direction: if signed_mph < 0.0 {
            std::f64::consts::FRAC_PI_2
        } else {
            3.0 * std::f64::consts::FRAC_PI_2
        },
        vertical_speed: 0.0,
    }
}

/// The FORWARD direction of [`solve_wind_truing`]: the lateral miss (meters, positive =
/// right of the line of sight) this request's load / atmosphere / twist / earth frame
/// predicts at `range_m` under a signed `crosswind_mph` (MBA-1392).
///
/// `request.observations` and `request.called_crosswind_mph` are ignored — only the model
/// half is used — and no validation is run, so a caller that has not validated its request
/// will simply see the solver's own error. Public because "what miss does N mph produce?"
/// is the question a wind-call drill actually asks, and because it lets any caller check
/// that the inversion round-trips on its own model rather than trusting it to.
pub fn modeled_miss_right_m(
    request: &WindTruingRequest,
    crosswind_mph: f64,
    range_m: f64,
) -> Result<f64, Box<dyn Error>> {
    let no_bc_segments: Option<Vec<BCSegmentData>> = None;
    let env = TruingEnvironment {
        wind: crosswind_conditions(crosswind_mph),
        twist: Some(request.twist),
        earth: request.earth,
    };
    let sample = crate::truing::solve_trajectory_sample(
        request.muzzle_velocity_fps,
        request.bc,
        request.drag_model,
        request.mass_gr,
        request.diameter_in,
        request.zero_distance_yd,
        range_m / 0.9144,
        request.sight_height_in,
        request.temperature_f,
        request.pressure_inhg,
        request.humidity_pct,
        request.altitude_ft,
        &no_bc_segments,
        &env,
        true, // interpolate: land exactly on the requested range
    )?;
    Ok(sample.lateral_m)
}

/// Back-solve the effective crosswind from every observed miss in `request` (MBA-1392).
///
/// Each observation is fitted independently with a bracketed root find on the constant
/// crosswind speed, against the real forward trajectory model (the truing core's own
/// solver assembly, zero-angle solve and atmosphere — the same one the drop-based truing
/// commands use, now sampled on the lateral axis and carrying spin drift, and optionally
/// Coriolis). See [`modeled_miss_right_m`] for that forward direction on its own.
///
/// Errors on an invalid request, on any trajectory-solver failure, and when an observed
/// miss cannot be produced by any crosswind inside `+/-`[`MAX_SOLVABLE_CROSSWIND_MPH`] —
/// which in practice means the miss was entered with the wrong sign, or is not a wind
/// effect at all.
pub fn solve_wind_truing(request: &WindTruingRequest) -> Result<WindTruingReport, Box<dyn Error>> {
    request.validate()?;

    // No BC5D velocity-banded segments on this path: the wind fit exposes one scalar BC,
    // matching the scalar-BC truing model (a banded schedule has no single BC to pair
    // with the fitted wind), and keeps the command entirely offline.
    let no_bc_segments: Option<Vec<BCSegmentData>> = None;

    let environment = |crosswind_mph: f64| TruingEnvironment {
        wind: crosswind_conditions(crosswind_mph),
        twist: Some(request.twist),
        earth: request.earth,
    };

    // Lateral (McCoy z, positive = right) at `range_yd` for a candidate crosswind.
    let lateral_at = |crosswind_mph: f64, range_yd: f64| -> Result<f64, Box<dyn Error>> {
        let sample = crate::truing::solve_trajectory_sample(
            request.muzzle_velocity_fps,
            request.bc,
            request.drag_model,
            request.mass_gr,
            request.diameter_in,
            request.zero_distance_yd,
            range_yd,
            request.sight_height_in,
            request.temperature_f,
            request.pressure_inhg,
            request.humidity_pct,
            request.altitude_ft,
            &no_bc_segments,
            &environment(crosswind_mph),
            true, // interpolate: land exactly on the observation range
        )?;
        Ok(sample.lateral_m)
    };

    let mut solutions = Vec::with_capacity(request.observations.len());
    for observation in &request.observations {
        let range_yd = observation.range_m / 0.9144;
        solutions.push(solve_one_observation(observation, range_yd, &lateral_at)?);
    }

    // Combine. Inverse-variance weighting only when every observation carried a sigma
    // (all-or-none, enforced by `validate`), otherwise the plain arithmetic mean. A
    // supplied sigma that could not be expressed in wind units is an error, NOT a quiet
    // demotion to unit weights: the caller asked for a weighted answer.
    let inverse_variance_weighted = solutions.iter().all(|s| s.sigma_m.is_some());
    if inverse_variance_weighted && solutions.iter().any(|s| s.solved_sigma_mph.is_none()) {
        return Err(
            "an observation does not move with crosswind at all, so its measurement sigma \
             cannot be expressed in wind units — drop that observation or its sigma"
                .into(),
        );
    }
    let (mean_crosswind_mph, mean_sigma_mph) = if inverse_variance_weighted {
        let mut weight_sum = 0.0;
        let mut weighted = 0.0;
        for solution in &solutions {
            let sigma = solution
                .solved_sigma_mph
                .expect("checked by inverse_variance_weighted");
            let weight = 1.0 / (sigma * sigma);
            weight_sum += weight;
            weighted += weight * solution.solved_crosswind_mph;
        }
        if weight_sum > 0.0 && weight_sum.is_finite() {
            (weighted / weight_sum, Some((1.0 / weight_sum).sqrt()))
        } else {
            return Err(
                "observed-miss sigmas produced a degenerate weighting (check that every \
                 sigma is positive and that the observations move with wind at all)"
                    .into(),
            );
        }
    } else {
        let sum: f64 = solutions.iter().map(|s| s.solved_crosswind_mph).sum();
        (sum / solutions.len() as f64, None)
    };

    let wind_call_factor = request
        .called_crosswind_mph
        .map(|called| mean_crosswind_mph / called);

    // Spin drift is always modelled (the twist rate is required for exactly that reason);
    // Coriolis only when an earth frame was supplied. Anything not listed as subtracted
    // is, by construction, still inside the solved wind.
    let mut subtracted_effects = vec!["spin drift".to_string()];
    let mut unsubtracted_effects = Vec::new();
    if request.earth.is_some() {
        subtracted_effects.push("Coriolis".to_string());
    } else {
        unsubtracted_effects
            .push("Coriolis (supply --latitude and --shot-direction to subtract it)".to_string());
    }

    Ok(WindTruingReport {
        solutions,
        mean_crosswind_mph,
        mean_sigma_mph,
        inverse_variance_weighted,
        called_crosswind_mph: request.called_crosswind_mph,
        wind_call_factor,
        subtracted_effects,
        unsubtracted_effects,
    })
}

/// Fit one observation. Bracketed false position (Illinois), mirroring the bracketing and
/// best-estimate-at-the-cap discipline of
/// [`crate::truing::calculate_true_velocity_local`]: a verified sign change over the full
/// solvable band, then an iteration that can never leave that bracket.
fn solve_one_observation(
    observation: &WindObservation,
    range_yd: f64,
    lateral_at: &impl Fn(f64, f64) -> Result<f64, Box<dyn Error>>,
) -> Result<WindTruingSolution, Box<dyn Error>> {
    let target = observation.miss_right_m;
    let residual = |crosswind_mph: f64| -> Result<f64, Box<dyn Error>> {
        Ok(lateral_at(crosswind_mph, range_yd)? - target)
    };

    let mut low = -MAX_SOLVABLE_CROSSWIND_MPH;
    let mut high = MAX_SOLVABLE_CROSSWIND_MPH;
    let mut f_low = residual(low)?;
    let mut f_high = residual(high)?;
    if f_low > 0.0 || f_high < 0.0 {
        return Err(format!(
            "no crosswind within +/-{MAX_SOLVABLE_CROSSWIND_MPH:.0} mph reproduces a {:.2} in \
             miss at {range_yd:.0} yd (that band spans {:.2} to {:.2} in of deflection) — check \
             the sign of --miss (positive = impact RIGHT of aim), the twist hand, and the load",
            target / 0.0254,
            (f_low + target) / 0.0254,
            (f_high + target) / 0.0254,
        )
        .into());
    }

    let mut solved = 0.0;
    let mut f_solved = 0.0;
    let mut iterations = 0u32;
    let mut converged = false;
    while iterations < WIND_SOLVE_MAX_ITERATIONS {
        iterations += 1;
        let denom = f_high - f_low;
        let mut candidate = if denom.abs() > f64::MIN_POSITIVE {
            high - f_high * (high - low) / denom
        } else {
            0.5 * (low + high)
        };
        // Never step outside the bracket (that is the whole point of keeping one).
        if !candidate.is_finite() || candidate <= low || candidate >= high {
            candidate = 0.5 * (low + high);
        }
        let f = residual(candidate)?;
        solved = candidate;
        f_solved = f;
        if f.abs() <= WIND_SOLVE_TOLERANCE_M || (high - low) <= WIND_SOLVE_MIN_BRACKET_MPH {
            converged = true;
            break;
        }
        // Illinois: halve the retained endpoint's function value so one side cannot stall.
        if f < 0.0 {
            low = candidate;
            f_low = f;
            f_high *= 0.5;
        } else {
            high = candidate;
            f_high = f;
            f_low *= 0.5;
        }
    }

    // How hard the observation actually pushes back on the wind, measured at the solution.
    let plus = lateral_at(solved + WIND_SENSITIVITY_STEP_MPH, range_yd)?;
    let minus = lateral_at(solved - WIND_SENSITIVITY_STEP_MPH, range_yd)?;
    let sensitivity_m_per_mph = (plus - minus) / (2.0 * WIND_SENSITIVITY_STEP_MPH);

    // The lateral the model produces with NO wind: spin drift, plus Coriolis when the
    // earth frame was supplied. This is what the fit accounted for instead of wind.
    let no_wind_lateral_m = lateral_at(0.0, range_yd)?;

    let solved_sigma_mph = observation.sigma_m.and_then(|sigma| {
        let slope = sensitivity_m_per_mph.abs();
        (slope > 0.0).then_some(sigma / slope)
    });

    Ok(WindTruingSolution {
        range_m: observation.range_m,
        observed_miss_right_m: target,
        sigma_m: observation.sigma_m,
        solved_crosswind_mph: solved,
        modeled_miss_right_m: target + f_solved,
        residual_m: f_solved,
        no_wind_lateral_m,
        sensitivity_m_per_mph,
        solved_sigma_mph,
        iterations,
        converged,
    })
}

/// Which rendering [`format_wind_truing_report`] should produce. Front-end-agnostic so the
/// native CLI's `OutputFormat` and the WASM terminal's `--output` string map onto ONE
/// formatter and cannot drift apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindTruingOutput {
    Table,
    Json,
    Csv,
}

/// Display units for a wind-truing rendering. Ranges follow the unit system; the observed
/// miss is always inches (a linear tape measurement, exactly like `--drop-unit in`).
struct WindTruingUnits {
    range_label: &'static str,
    speed_label: &'static str,
    range_scale: f64,
    speed_scale: f64,
}

impl WindTruingUnits {
    fn for_system(units: UnitSystem) -> Self {
        match units {
            UnitSystem::Imperial => Self {
                range_label: "yd",
                speed_label: "mph",
                range_scale: 1.0 / 0.9144,
                speed_scale: 1.0,
            },
            UnitSystem::Metric => Self {
                range_label: "m",
                speed_label: "m/s",
                range_scale: 1.0,
                speed_scale: MPH_TO_MPS,
            },
        }
    }

    fn range(&self, range_m: f64) -> f64 {
        range_m * self.range_scale
    }

    fn speed(&self, mph: f64) -> f64 {
        mph * self.speed_scale
    }
}

/// Meters to inches, for the linear miss columns.
fn inches(meters: f64) -> f64 {
    meters / 0.0254
}

/// The wind-truing report as a JSON document (MBA-1392).
///
/// Split out from [`format_wind_truing_report`] so the conditional fields (the optional
/// sigma, called wind and correction factor, which are `null` when absent) are testable on
/// the host — `wasm.rs` is wasm32-gated, so a formatter that only existed inside it could
/// never be asserted against natively. Same reason `drag_coefficient_json_value` exists.
pub fn wind_truing_json_value(report: &WindTruingReport, units: UnitSystem) -> serde_json::Value {
    let u = WindTruingUnits::for_system(units);
    let observations: Vec<serde_json::Value> = report
        .solutions
        .iter()
        .map(|s| {
            serde_json::json!({
                format!("range_{}", u.range_label): u.range(s.range_m),
                "miss_right_in": inches(s.observed_miss_right_m),
                "miss_sigma_in": s.sigma_m.map(inches),
                "no_wind_lateral_in": inches(s.no_wind_lateral_m),
                "solved_crosswind": u.speed(s.solved_crosswind_mph),
                "solved_crosswind_sigma": s.solved_sigma_mph.map(|v| u.speed(v)),
                "sensitivity_in_per_mph": inches(s.sensitivity_m_per_mph),
                "residual_in": inches(s.residual_m),
                "iterations": s.iterations,
                "converged": s.converged,
            })
        })
        .collect();

    serde_json::json!({
        "effective_crosswind": u.speed(report.mean_crosswind_mph),
        "effective_crosswind_sigma": report.mean_sigma_mph.map(|v| u.speed(v)),
        "inverse_variance_weighted": report.inverse_variance_weighted,
        "called_crosswind": report.called_crosswind_mph.map(|v| u.speed(v)),
        "wind_call_factor": report.wind_call_factor,
        "observations": observations,
        "effects_subtracted": report.subtracted_effects,
        "effects_not_subtracted": report.unsubtracted_effects,
        "legend": {
            "units": {
                "range": u.range_label,
                "miss": "in",
                "wind_speed": u.speed_label,
            },
            "signs": "--miss positive = impact right of aim; solved crosswind positive = \
                      wind from the shooter's left (9 o'clock) pushing impacts right",
        },
    })
}

/// Render a [`WindTruingReport`] (MBA-1392).
///
/// ONE formatter for both front ends: the native CLI prints the returned string and the
/// WASM terminal returns it, so the two surfaces are byte-identical by construction rather
/// than by a replicated printer that has to be kept in sync.
pub fn format_wind_truing_report(
    report: &WindTruingReport,
    units: UnitSystem,
    output: WindTruingOutput,
) -> String {
    let u = WindTruingUnits::for_system(units);
    match output {
        WindTruingOutput::Json => {
            match serde_json::to_string_pretty(&wind_truing_json_value(report, units)) {
                Ok(s) => format!("{s}\n"),
                Err(e) => format!("Error serializing JSON: {e}\n"),
            }
        }
        WindTruingOutput::Csv => {
            let mut out = String::new();
            out.push_str(&format!(
                "range_{},miss_right_in,miss_sigma_in,no_wind_lateral_in,solved_crosswind_{},\
                 sensitivity_in_per_mph,residual_in,iterations,converged\n",
                u.range_label, u.speed_label
            ));
            for s in &report.solutions {
                out.push_str(&format!(
                    "{:.1},{:+.3},{},{:+.3},{:+.3},{:.4},{:+.4},{},{}\n",
                    u.range(s.range_m),
                    inches(s.observed_miss_right_m),
                    match s.sigma_m {
                        Some(sigma) => format!("{:.3}", inches(sigma)),
                        None => String::new(),
                    },
                    inches(s.no_wind_lateral_m),
                    u.speed(s.solved_crosswind_mph),
                    inches(s.sensitivity_m_per_mph),
                    inches(s.residual_m),
                    s.iterations,
                    s.converged,
                ));
            }
            out.push('\n');
            out.push_str(&format!(
                "effective_crosswind_{},effective_crosswind_sigma_{},inverse_variance_weighted,\
                 called_crosswind_{},wind_call_factor\n",
                u.speed_label, u.speed_label, u.speed_label
            ));
            out.push_str(&format!(
                "{:+.3},{},{},{},{}\n",
                u.speed(report.mean_crosswind_mph),
                match report.mean_sigma_mph {
                    Some(sigma) => format!("{:.3}", u.speed(sigma)),
                    None => String::new(),
                },
                report.inverse_variance_weighted,
                match report.called_crosswind_mph {
                    Some(called) => format!("{:+.3}", u.speed(called)),
                    None => String::new(),
                },
                match report.wind_call_factor {
                    Some(factor) => format!("{factor:.4}"),
                    None => String::new(),
                },
            ));
            out
        }
        WindTruingOutput::Table => {
            let mut out = String::new();
            out.push('\n');
            out.push_str("=== EFFECTIVE WIND TRUING (from observed horizontal miss) ===\n");
            out.push('\n');
            out.push_str(&format!(
                "  {:>10}  {:>12}  {:>14}  {:>16}  {:>10}\n",
                format!("Range ({})", u.range_label),
                "Miss (in)",
                "Spin/Cor (in)",
                format!("Wind ({})", u.speed_label),
                "Resid (in)",
            ));
            out.push_str(&format!("  {}\n", "-".repeat(70)));
            for s in &report.solutions {
                out.push_str(&format!(
                    "  {:>10.1}  {:>+12.2}  {:>+14.2}  {:>+16.2}  {:>+10.3}\n",
                    u.range(s.range_m),
                    inches(s.observed_miss_right_m),
                    inches(s.no_wind_lateral_m),
                    u.speed(s.solved_crosswind_mph),
                    inches(s.residual_m),
                ));
            }
            out.push_str(&format!("  {}\n", "-".repeat(70)));
            out.push('\n');
            let n = report.solutions.len();
            out.push_str(&format!(
                "  Effective crosswind: {:>+8.2} {}{}\n",
                u.speed(report.mean_crosswind_mph),
                u.speed_label,
                match report.mean_sigma_mph {
                    Some(sigma) => format!(
                        "  +/- {:.2} {} (inverse-variance weighted over {n} observations)",
                        u.speed(sigma),
                        u.speed_label
                    ),
                    None if n > 1 => format!("  (mean of {n} observations)"),
                    None => String::new(),
                }
            ));
            if let (Some(called), Some(factor)) =
                (report.called_crosswind_mph, report.wind_call_factor)
            {
                out.push_str(&format!(
                    "  Called wind:         {:>+8.2} {}   ->  wind-call correction factor {:.2}\n",
                    u.speed(called),
                    u.speed_label,
                    factor
                ));
                out.push_str(&format!(
                    "    (multiply your wind calls by {factor:.2} to match what actually hit)\n"
                ));
            }
            if !report.subtracted_effects.is_empty() {
                out.push_str(&format!(
                    "  Effects subtracted:  {}\n",
                    report.subtracted_effects.join(", ")
                ));
            }
            if !report.unsubtracted_effects.is_empty() {
                out.push_str(&format!(
                    "  NOT subtracted (absorbed into the solved wind): {}\n",
                    report.unsubtracted_effects.join(", ")
                ));
            }
            for s in &report.solutions {
                if inches(s.sensitivity_m_per_mph).abs() < MIN_WIND_SENSITIVITY_IN_PER_MPH {
                    out.push_str(&format!(
                        "  note: the observation at {:.1} {} moves only {:.2} in per mph of \
                         crosswind (guide: {MIN_WIND_SENSITIVITY_IN_PER_MPH:.2} in/mph); the \
                         wind fitted from it is weakly identified\n",
                        u.range(s.range_m),
                        u.range_label,
                        inches(s.sensitivity_m_per_mph).abs(),
                    ));
                }
                if !s.converged {
                    out.push_str(&format!(
                        "  note: the fit at {:.1} {} did not fully converge after {} iterations; \
                         the value shown is the best estimate\n",
                        u.range(s.range_m),
                        u.range_label,
                        s.iterations,
                    ));
                }
            }
            out.push('\n');
            out.push_str(
                "  Signs: --miss positive = impact RIGHT of aim. Solved wind positive = wind\n\
                 \x20        FROM the shooter's LEFT (9 o'clock) pushing impacts right; negative\n\
                 \x20        = FROM the right pushing left. Wind-FROM convention throughout\n\
                 \x20        (0 = headwind, as of the 0.19.0 wind-direction sign fix).\n",
            );
            out.push('\n');
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_request(observations: Vec<WindObservation>) -> WindTruingRequest {
        WindTruingRequest {
            observations,
            muzzle_velocity_fps: 2700.0,
            bc: 0.475,
            drag_model: DragModelArg::G7,
            mass_gr: 168.0,
            diameter_in: 0.308,
            zero_distance_yd: 100.0,
            sight_height_in: 2.0,
            temperature_f: 59.0,
            pressure_inhg: 29.92,
            humidity_pct: 50.0,
            altitude_ft: 0.0,
            twist: TruingTwist {
                rate_in: 11.0,
                right_hand: true,
            },
            earth: None,
            called_crosswind_mph: None,
        }
    }

    /// The forward model the fit inverts — used to manufacture "observed" misses for the
    /// round-trip tests from a KNOWN wind.
    fn modeled_lateral_m(request: &WindTruingRequest, crosswind_mph: f64, range_m: f64) -> f64 {
        modeled_miss_right_m(request, crosswind_mph, range_m).expect("forward model must solve")
    }

    /// Round trip: take a KNOWN crosswind, read the lateral miss it produces at three
    /// ranges out of the forward model, feed those misses back in as observations, and
    /// recover the wind. Tolerance is 0.02 mph — three orders of magnitude tighter than
    /// anyone's wind call, so this pins the inversion, not just its ballpark.
    #[test]
    fn round_trip_recovers_a_known_crosswind_at_three_ranges() {
        let known_mph = 7.5;
        let ranges_m = [274.32, 457.2, 640.08]; // 300 / 500 / 700 yd
        let template = base_request(Vec::new());
        let observations = ranges_m
            .iter()
            .map(|range_m| WindObservation {
                range_m: *range_m,
                miss_right_m: modeled_lateral_m(&template, known_mph, *range_m),
                sigma_m: None,
            })
            .collect();

        let report = solve_wind_truing(&base_request(observations)).expect("wind fit must solve");
        assert_eq!(report.solutions.len(), 3);
        for solution in &report.solutions {
            assert!(solution.converged, "{solution:?}");
            assert!(
                (solution.solved_crosswind_mph - known_mph).abs() < 0.02,
                "recovered {} mph at {} m, expected {known_mph}",
                solution.solved_crosswind_mph,
                solution.range_m
            );
        }
        assert!((report.mean_crosswind_mph - known_mph).abs() < 0.02);
        assert!(!report.inverse_variance_weighted);
        assert!(report.mean_sigma_mph.is_none());
    }

    /// Sign pins, both directions. A miss to the RIGHT must solve to a POSITIVE wind
    /// (from the shooter's left, pushing right); a miss to the LEFT must solve NEGATIVE.
    /// These two assertions are the contract the help text documents.
    #[test]
    fn miss_right_solves_positive_and_miss_left_solves_negative() {
        let template = base_request(Vec::new());
        let range_m = 457.2; // 500 yd
        let right_miss = modeled_lateral_m(&template, 8.0, range_m);
        let left_miss = modeled_lateral_m(&template, -8.0, range_m);
        assert!(right_miss > 0.0, "a left-hand wind must push impacts right");
        assert!(left_miss < 0.0, "a right-hand wind must push impacts left");

        let right = solve_wind_truing(&base_request(vec![WindObservation {
            range_m,
            miss_right_m: right_miss,
            sigma_m: None,
        }]))
        .expect("right-miss fit must solve");
        assert!(
            right.mean_crosswind_mph > 0.0,
            "right miss must solve to a positive (left-hand, right-pushing) wind, got {}",
            right.mean_crosswind_mph
        );
        assert!((right.mean_crosswind_mph - 8.0).abs() < 0.02);

        let left = solve_wind_truing(&base_request(vec![WindObservation {
            range_m,
            miss_right_m: left_miss,
            sigma_m: None,
        }]))
        .expect("left-miss fit must solve");
        assert!(
            left.mean_crosswind_mph < 0.0,
            "left miss must solve to a negative (right-hand, left-pushing) wind, got {}",
            left.mean_crosswind_mph
        );
        assert!((left.mean_crosswind_mph + 8.0).abs() < 0.02);
    }

    /// Spin-drift subtraction: a ZERO-wind trajectory still lands right of aim (right-hand
    /// twist). Feeding that pure spin drift in as the observed miss must solve to ~zero
    /// wind — the drift is attributed to spin, not to a phantom wind.
    #[test]
    fn pure_spin_drift_solves_to_zero_wind() {
        let template = base_request(Vec::new());
        let range_m = 640.08; // 700 yd
        let spin_only = modeled_lateral_m(&template, 0.0, range_m);
        assert!(
            spin_only > 0.05,
            "a 1:11 right-hand twist must drift measurably right at 700 yd, got {spin_only} m"
        );

        let report = solve_wind_truing(&base_request(vec![WindObservation {
            range_m,
            miss_right_m: spin_only,
            sigma_m: None,
        }]))
        .expect("spin-only fit must solve");
        assert!(
            report.mean_crosswind_mph.abs() < 0.02,
            "pure spin drift must solve to ~0 wind, got {}",
            report.mean_crosswind_mph
        );
        // ... and the report must say the drift was accounted for, not silently eaten.
        assert!(report
            .solutions
            .iter()
            .all(|s| (s.no_wind_lateral_m - spin_only).abs() < 1e-9));
        assert!(report
            .subtracted_effects
            .iter()
            .any(|e| e.contains("spin drift")));
    }

    /// A left-hand twist drifts the other way, so the SAME right-of-aim miss must solve to
    /// a stronger wind than a right-hand twist needs. Pins that the twist hand actually
    /// reaches the forward model rather than being cosmetic.
    #[test]
    fn twist_hand_changes_the_solved_wind() {
        let range_m = 640.08;
        let observation = WindObservation {
            range_m,
            miss_right_m: 0.25,
            sigma_m: None,
        };
        let right_hand = solve_wind_truing(&base_request(vec![observation])).expect("solve");
        let mut left = base_request(vec![observation]);
        left.twist.right_hand = false;
        let left_hand = solve_wind_truing(&left).expect("solve");
        assert!(
            left_hand.mean_crosswind_mph > right_hand.mean_crosswind_mph + 0.1,
            "left-hand twist ({}) must need more right-pushing wind than right-hand ({})",
            left_hand.mean_crosswind_mph,
            right_hand.mean_crosswind_mph
        );
    }

    /// Coriolis subtraction: supplying a latitude and shot azimuth changes the zero-wind
    /// lateral (Coriolis is now modelled), so the same observed miss solves to a different
    /// wind, and the report moves Coriolis from "not subtracted" to "subtracted".
    #[test]
    fn coriolis_is_subtracted_when_latitude_and_azimuth_are_supplied() {
        let range_m = 914.4; // 1000 yd, where Coriolis is actually measurable
        let observation = WindObservation {
            range_m,
            miss_right_m: 0.30,
            sigma_m: None,
        };
        let without = solve_wind_truing(&base_request(vec![observation])).expect("solve");
        let mut with_earth = base_request(vec![observation]);
        with_earth.earth = Some(TruingEarthFrame {
            latitude_deg: 45.0,
            shot_azimuth_deg: 90.0, // due East, where the Coriolis lateral is largest
        });
        let with = solve_wind_truing(&with_earth).expect("solve");

        assert!(without
            .unsubtracted_effects
            .iter()
            .any(|e| e.contains("Coriolis")));
        assert!(without.subtracted_effects.iter().all(|e| e != "Coriolis"));
        assert!(with.unsubtracted_effects.is_empty());
        assert!(with.subtracted_effects.iter().any(|e| e == "Coriolis"));
        assert!(
            (with.solutions[0].no_wind_lateral_m - without.solutions[0].no_wind_lateral_m).abs()
                > 1e-4,
            "modelling Coriolis must change the zero-wind lateral"
        );
        assert!(
            (with.mean_crosswind_mph - without.mean_crosswind_mph).abs() > 1e-3,
            "modelling Coriolis must change the solved wind"
        );
    }

    /// The wind-call correction factor is solved / called, and its sign survives: calling
    /// the wrong SIDE of the wind gives a negative factor rather than a plausible-looking
    /// positive one.
    #[test]
    fn wind_call_factor_is_solved_over_called_and_keeps_its_sign() {
        let template = base_request(Vec::new());
        let range_m = 457.2;
        let miss = modeled_lateral_m(&template, 9.0, range_m);
        let observation = WindObservation {
            range_m,
            miss_right_m: miss,
            sigma_m: None,
        };

        let mut under_called = base_request(vec![observation]);
        under_called.called_crosswind_mph = Some(6.0);
        let report = solve_wind_truing(&under_called).expect("solve");
        let factor = report.wind_call_factor.expect("factor");
        assert!(
            (factor - 9.0 / 6.0).abs() < 0.01,
            "expected ~1.5, got {factor}"
        );

        let mut wrong_side = base_request(vec![observation]);
        wrong_side.called_crosswind_mph = Some(-6.0);
        let flipped = solve_wind_truing(&wrong_side)
            .expect("solve")
            .wind_call_factor
            .expect("factor");
        assert!(flipped < 0.0, "a wrong-side call must read negative: {flipped}");
    }

    /// Sigmas: all-or-none. Every observation weighted -> inverse-variance mean with a
    /// reported sigma; none weighted -> plain mean; a mix is a hard error rather than a
    /// silent blend of weighting schemes.
    #[test]
    fn sigmas_are_all_or_none_and_drive_inverse_variance_weighting() {
        let template = base_request(Vec::new());
        let near = 274.32;
        let far = 640.08;
        let weighted = vec![
            WindObservation {
                range_m: near,
                miss_right_m: modeled_lateral_m(&template, 6.0, near),
                sigma_m: Some(0.25 * 0.0254),
            },
            WindObservation {
                range_m: far,
                miss_right_m: modeled_lateral_m(&template, 6.0, far),
                sigma_m: Some(0.25 * 0.0254),
            },
        ];
        let report = solve_wind_truing(&base_request(weighted)).expect("solve");
        assert!(report.inverse_variance_weighted);
        let sigma = report.mean_sigma_mph.expect("weighted mean sigma");
        assert!(sigma > 0.0 && sigma.is_finite());
        // The long-range observation is far more sensitive to wind, so its propagated
        // sigma must be the smaller of the two (it carries the most weight).
        let near_sigma = report.solutions[0].solved_sigma_mph.expect("sigma");
        let far_sigma = report.solutions[1].solved_sigma_mph.expect("sigma");
        assert!(far_sigma < near_sigma, "{far_sigma} !< {near_sigma}");
        assert!(sigma <= far_sigma + 1e-12);

        let mixed = base_request(vec![
            WindObservation {
                range_m: near,
                miss_right_m: 0.1,
                sigma_m: Some(0.006),
            },
            WindObservation {
                range_m: far,
                miss_right_m: 0.2,
                sigma_m: None,
            },
        ]);
        let error = mixed.validate().unwrap_err();
        assert!(error.contains("every observed miss or on none"), "{error}");
    }

    /// An unreachable miss (wrong sign, or simply not a wind effect) is rejected with a
    /// diagnostic naming the solvable band, instead of being clamped into a fake answer.
    #[test]
    fn an_unreachable_miss_is_rejected_with_the_solvable_band() {
        let error = solve_wind_truing(&base_request(vec![WindObservation {
            range_m: 274.32,
            miss_right_m: 25.0, // 25 metres right at 300 yd: no wind does that
            sigma_m: None,
        }]))
        .unwrap_err()
        .to_string();
        assert!(error.contains("no crosswind within"), "{error}");
        assert!(error.contains("check the sign of --miss"), "{error}");
    }

    /// MBA-1358 / design R5-DIRECTION: `--miss` values are LINEAR inches off the target,
    /// not dial readings, so a scope tracking correction factor must NOT touch them. The
    /// structural guarantee is that no CF can reach this solver at all — there is no field
    /// for one — and this test pins the consequence: a windage CF applied the way the
    /// DIALED truing path applies it (observation x CF) would move the answer, so it must
    /// never be applied here.
    #[test]
    fn windage_cf_does_not_alter_the_wind_solve() {
        let template = base_request(Vec::new());
        let range_m = 457.2;
        let miss = modeled_lateral_m(&template, 7.0, range_m);
        let observation = WindObservation {
            range_m,
            miss_right_m: miss,
            sigma_m: None,
        };
        let solved = solve_wind_truing(&base_request(vec![observation]))
            .expect("solve")
            .mean_crosswind_mph;
        assert!((solved - 7.0).abs() < 0.02);

        // What the dialed path would have done to a 0.95 CF observation. It changes the
        // answer materially, which is exactly why linear inputs must be left alone.
        let windage_cf = 0.95;
        let cf_applied = solve_wind_truing(&base_request(vec![WindObservation {
            range_m,
            miss_right_m: miss * windage_cf,
            sigma_m: None,
        }]))
        .expect("solve")
        .mean_crosswind_mph;
        assert!(
            (cf_applied - solved).abs() > 0.1,
            "a CF-scaled observation must NOT be equivalent to the linear one \
             ({cf_applied} vs {solved}); --miss therefore takes no CF"
        );
    }

    /// The JSON emit helper: conditional fields are explicit `null`s (never dropped keys),
    /// so consumers can tell "not supplied" from "absent field".
    #[test]
    fn json_value_nulls_absent_optional_fields() {
        let template = base_request(Vec::new());
        let range_m = 457.2;
        let report = solve_wind_truing(&base_request(vec![WindObservation {
            range_m,
            miss_right_m: modeled_lateral_m(&template, 5.0, range_m),
            sigma_m: None,
        }]))
        .expect("solve");
        let value = wind_truing_json_value(&report, UnitSystem::Imperial);
        assert!(value["called_crosswind"].is_null());
        assert!(value["wind_call_factor"].is_null());
        assert!(value["effective_crosswind_sigma"].is_null());
        assert!(value["observations"][0]["miss_sigma_in"].is_null());
        assert!(value["observations"][0]["solved_crosswind_sigma"].is_null());
        assert_eq!(value["legend"]["units"]["wind_speed"], "mph");
        assert_eq!(value["legend"]["units"]["miss"], "in");
        assert_eq!(
            value["effective_crosswind"].as_f64().expect("f64").round(),
            5.0
        );

        // Metric renders the same wind in m/s and the same ranges in meters, while the
        // linear miss stays in inches (the --drop-unit in precedent).
        let metric = wind_truing_json_value(&report, UnitSystem::Metric);
        assert_eq!(metric["legend"]["units"]["wind_speed"], "m/s");
        assert_eq!(metric["legend"]["units"]["range"], "m");
        assert_eq!(metric["legend"]["units"]["miss"], "in");
        let mps = metric["effective_crosswind"].as_f64().expect("f64");
        let mph = value["effective_crosswind"].as_f64().expect("f64");
        assert!((mps - mph * MPH_TO_MPS).abs() < 1e-12);
    }

    /// Parser contract: RANGE follows the unit system, the offset and sigma are inches in
    /// both, and malformed tokens are rejected with a usable message.
    #[test]
    fn parse_wind_observation_units_and_errors() {
        let imperial = parse_wind_observation("600:8.5", UnitSystem::Imperial).expect("parse");
        assert!((imperial.range_m - 600.0 * 0.9144).abs() < 1e-12);
        assert!((imperial.miss_right_m - 8.5 * 0.0254).abs() < 1e-12);
        assert!(imperial.sigma_m.is_none());

        let metric = parse_wind_observation("550:-8.5:0.75", UnitSystem::Metric).expect("parse");
        assert!((metric.range_m - 550.0).abs() < 1e-12);
        assert!((metric.miss_right_m + 8.5 * 0.0254).abs() < 1e-12);
        assert!((metric.sigma_m.expect("sigma") - 0.75 * 0.0254).abs() < 1e-12);

        for bad in ["600", "600:8.5:0.1:2", "600:right", "abc:8.5", "600:nan"] {
            assert!(
                parse_wind_observation(bad, UnitSystem::Imperial).is_err(),
                "'{bad}' should not parse"
            );
        }
    }

    /// Validation rejects degenerate requests before spending a single trajectory solve.
    #[test]
    fn validation_rejects_degenerate_requests() {
        assert!(base_request(Vec::new())
            .validate()
            .unwrap_err()
            .contains("at least one"));

        let duplicate = base_request(vec![
            WindObservation {
                range_m: 457.2,
                miss_right_m: 0.2,
                sigma_m: None,
            },
            WindObservation {
                range_m: 457.2,
                miss_right_m: 0.3,
                sigma_m: None,
            },
        ]);
        assert!(duplicate
            .validate()
            .unwrap_err()
            .contains("duplicate observation range"));

        let mut bad_twist = base_request(vec![WindObservation {
            range_m: 457.2,
            miss_right_m: 0.2,
            sigma_m: None,
        }]);
        bad_twist.twist.rate_in = 0.0;
        assert!(bad_twist
            .validate()
            .unwrap_err()
            .contains("twist rate must be positive"));

        let mut zero_call = base_request(vec![WindObservation {
            range_m: 457.2,
            miss_right_m: 0.2,
            sigma_m: None,
        }]);
        zero_call.called_crosswind_mph = Some(0.0);
        assert!(zero_call.validate().unwrap_err().contains("non-zero"));
    }

    /// Bridge contract: a `WindTruingRequest` deserializes from JSON with the field names
    /// documented on the struct, and a solved `WindTruingReport` serializes back out.
    #[test]
    fn request_deserializes_and_report_serializes() {
        let json = serde_json::json!({
            "observations": [{"range_m": 457.2, "miss_right_m": 0.315, "sigma_m": null}],
            "muzzle_velocity_fps": 2700.0, "bc": 0.243, "drag_model": "g7",
            "mass_gr": 168.0, "diameter_in": 0.308, "zero_distance_yd": 100.0,
            "sight_height_in": 2.0, "temperature_f": 59.0, "pressure_inhg": 29.92,
            "humidity_pct": 50.0, "altitude_ft": 0.0,
            "twist": {"rate_in": 11.0, "right_hand": true},
            "earth": null, "called_crosswind_mph": null
        });
        let req: WindTruingRequest =
            serde_json::from_value(json).expect("request deserializes");
        let report = solve_wind_truing(&req).expect("solves");
        let out = serde_json::to_value(&report).expect("report serializes");
        assert!(out["mean_crosswind_mph"].is_number());
        assert!(out["solutions"].as_array().unwrap().len() == 1);
    }
}
