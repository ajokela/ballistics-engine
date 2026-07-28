//! WEZ (Weapon Employment Zone) sweep core -- MBA-1317, extracted MBA-1343 Phase B.
//!
//! `monte-carlo --wez` reports hit probability vs range for a fixed target size, treating the
//! shooter's wind-CALL error (how well they estimate the current wind) as a source of dispersion
//! distinct from the ballistic --wind-std (gust-to-gust physical variability). See the "WEZ" doc
//! section in CLI_USAGE.md for a worked example.
//!
//! Extracted from the CLI binary so non-CLI front ends (e.g. the WASM terminal) can reuse the
//! exact compute path. All rendering (summary table / statistics CSV / full JSON) stays with the
//! front ends; this module goes as far as building a [`WezResult`].

use std::error::Error;

use nalgebra::Vector3;
use serde::Serialize;

use crate::cli_api::UnitSystem;
use crate::drag::DragTable;
use crate::{
    AtmosphericConditions, BallisticInputs, BallisticsError, DragModel, MonteCarloParams,
    MonteCarloResults, TrajectorySolver, WindConditions,
};

/// A parsed `--target-size` value, still in the CLI's chosen unit (inches imperial / cm
/// metric) -- call [`TargetSize::to_metric`] before using it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TargetSize {
    /// Full width (lateral) x height (vertical), e.g. an 18"x30" plate.
    Rect { width: f64, height: f64 },
    /// A circular radius, matching `--target-radius`'s existing hit semantics but expressed in
    /// target-size units instead of range units.
    Radius(f64),
}

/// Parse a `--target-size` argument: `WIDTHxHEIGHT` (e.g. `18x30`) for a rectangle, or a bare
/// number (e.g. `12`) for a circular radius fallback. Case-insensitive on the `x` separator.
pub fn parse_target_size(spec: &str) -> Result<TargetSize, String> {
    let trimmed = spec.trim();
    if trimmed.is_empty() {
        return Err("expected a size like \"18x30\" or a single radius like \"12\"".to_string());
    }

    let x_positions: Vec<usize> = trimmed
        .char_indices()
        .filter(|(_, c)| *c == 'x' || *c == 'X')
        .map(|(i, _)| i)
        .collect();

    match x_positions.len() {
        0 => {
            let radius: f64 = trimmed
                .parse()
                .map_err(|_| format!("\"{trimmed}\" is not a number or a WIDTHxHEIGHT pair"))?;
            if !(radius.is_finite() && radius > 0.0) {
                return Err(format!(
                    "radius must be a positive, finite number, got {radius}"
                ));
            }
            Ok(TargetSize::Radius(radius))
        }
        1 => {
            let idx = x_positions[0];
            let width_str = &trimmed[..idx];
            let height_str = &trimmed[idx + 1..];
            let width: f64 = width_str
                .trim()
                .parse()
                .map_err(|_| format!("\"{}\" is not a valid width", width_str.trim()))?;
            let height: f64 = height_str
                .trim()
                .parse()
                .map_err(|_| format!("\"{}\" is not a valid height", height_str.trim()))?;
            if !(width.is_finite() && width > 0.0 && height.is_finite() && height > 0.0) {
                return Err(format!(
                    "width and height must be positive, finite numbers, got {width}x{height}"
                ));
            }
            Ok(TargetSize::Rect { width, height })
        }
        _ => Err(format!(
            "\"{trimmed}\" has more than one 'x' separator; expected WIDTHxHEIGHT or a single radius"
        )),
    }
}

/// A [`TargetSize`] converted to meters, ready for [`MonteCarloResults`]'s
/// hit-probability methods.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TargetSizeMetric {
    Rect { width_m: f64, height_m: f64 },
    Radius { radius_m: f64 },
}

/// WEZ target-size length (MBA-1317): inches under imperial, CENTIMETERS (not mm -- target
/// sizes like an 18"x30" plate are naturally cm-scale under metric) under metric.
fn target_size_to_metric(val: f64, units: UnitSystem) -> f64 {
    match units {
        UnitSystem::Metric => val * 0.01,     // cm to meters
        UnitSystem::Imperial => val * 0.0254, // inches to meters
    }
}

impl TargetSize {
    pub fn to_metric(self, units: UnitSystem) -> TargetSizeMetric {
        match self {
            TargetSize::Rect { width, height } => TargetSizeMetric::Rect {
                width_m: target_size_to_metric(width, units),
                height_m: target_size_to_metric(height, units),
            },
            TargetSize::Radius(radius) => TargetSizeMetric::Radius {
                radius_m: target_size_to_metric(radius, units),
            },
        }
    }
}

/// Which WEZ variance-attribution bucket a miss-variance source belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum WezErrorBucket {
    /// The shooter's wind-call error (`--wind-call-error`).
    WindCall,
    /// Muzzle-velocity standard deviation (`--velocity-std`).
    MvSd,
    /// Everything else: mechanical/ammo group dispersion (angle, azimuth, BC) plus the
    /// *ballistic* (non-call) share of wind uncertainty (`--wind-std`, `--wind-direction-std`).
    Other,
}

impl std::fmt::Display for WezErrorBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // MBA-1337 w2: one spelling everywhere. These must match the serde
        // rename_all = "snake_case" values the -o full JSON contract shipped with
        // (0.25.0), so summary/CSV/JSON all agree on the same strings.
        let label = match self {
            WezErrorBucket::WindCall => "wind_call",
            WezErrorBucket::MvSd => "mv_sd",
            WezErrorBucket::Other => "other",
        };
        write!(f, "{label}")
    }
}

/// Per-range WEZ miss-variance attribution shares. `wind_call + mv_sd + other` sums to ~1.0
/// whenever at least one modeled source has nonzero uncertainty; all fields are exactly 0.0 in
/// the fully deterministic (zero-uncertainty) case, where there is no dominant source.
#[derive(Debug, Clone, Copy, Default, Serialize)]
struct WezVarianceShares {
    wind_call: f64,
    mv_sd: f64,
    other: f64,
}

impl WezVarianceShares {
    /// The largest nonzero share, or `None` if every share is zero (nothing to attribute).
    fn dominant(&self) -> Option<WezErrorBucket> {
        [
            (WezErrorBucket::WindCall, self.wind_call),
            (WezErrorBucket::MvSd, self.mv_sd),
            (WezErrorBucket::Other, self.other),
        ]
        .into_iter()
        .filter(|(_, share)| *share > 0.0)
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(bucket, _)| bucket)
    }
}

/// Solve a single deterministic trajectory and return its target-plane impact position.
///
/// The caller is responsible for setting `inputs.muzzle_velocity` to whatever value it wants
/// solved (e.g. baseline + one sigma for the MV-SD sensitivity). The real Monte Carlo sampler
/// instead applies a sampled velocity *delta* after `TrajectorySolver::new` resolves any
/// powder-temperature curve (MBA-1176), because `TrajectorySolver` doesn't expose that resolved
/// value to callers outside `cli_api`. That distinction is a no-op here: `monte-carlo` (and so
/// `--wez`) never sets `powder_temp_curve` on its `BallisticInputs`, so there is no curve for
/// `TrajectorySolver::new` to resolve and a plain pre-construction velocity assignment is
/// equivalent to the sampler's post-construction delta.
///
/// Propagates `solve()`'s error instead of unwrapping it: neither the baseline solve (whose
/// inputs are the CLI's raw, not-yet-validated `-v`/`-a`/etc. values -- clap's own range checks
/// permit e.g. `-v 0`, which `solve()` rejects) nor a perturbed one-sigma solve (MV/BC/angle/wind
/// nudged off a valid baseline) is guaranteed to stay within `solve()`'s validity gate (MBA-1317).
fn wez_solve_target_plane(
    inputs: BallisticInputs,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
    solver_max_range: f64,
    target_distance_m: f64,
) -> Result<Vector3<f64>, BallisticsError> {
    let mut solver = TrajectorySolver::new(inputs, wind, atmosphere);
    solver.set_max_range(solver_max_range);
    let result = solver.solve()?;
    // A successful `solve()` always produces a non-empty trajectory (every solve path errors
    // out on an empty point list before returning `Ok`), so `position_at_range` -- which only
    // returns `None` for an empty trajectory, clamping to the last point otherwise -- cannot
    // fail here. See cli_api::TrajectoryResult::position_at_range.
    Ok(result
        .position_at_range(target_distance_m)
        .expect("WEZ attribution solve: non-empty trajectory always has a last point"))
}

/// One source's target-plane variance contribution: `sigma`'s one-standard-deviation
/// displacement from `baseline`, or `0.0` when `sigma` is non-positive (that source is disabled).
/// `inputs` and `wind` must already have that one sigma applied by the caller.
///
/// Propagates a failed perturbed solve rather than unwrapping it (MBA-1317) -- a one-sigma nudge
/// off a valid baseline is not itself guaranteed to stay within `solve()`'s validity gate.
fn wez_source_variance(
    sigma: f64,
    inputs: BallisticInputs,
    wind: WindConditions,
    atmosphere: &AtmosphericConditions,
    solver_max_range: f64,
    target_distance_m: f64,
    baseline: &Vector3<f64>,
) -> Result<f64, BallisticsError> {
    if sigma.is_nan() || sigma <= 0.0 {
        return Ok(0.0);
    }
    let perturbed = wez_solve_target_plane(
        inputs,
        wind,
        atmosphere.clone(),
        solver_max_range,
        target_distance_m,
    )?;
    let dy = perturbed.y - baseline.y;
    let dz = perturbed.z - baseline.z;
    Ok(dy * dy + dz * dz)
}

/// Linearized (one-sigma finite-difference) WEZ miss-variance attribution at a single range.
///
/// A full "decomposed re-run" attribution -- zero out each source in turn and re-run the whole
/// Monte Carlo sample set for it -- would multiply the sweep's cost by the number of buckets,
/// which does not fit a "finishes in seconds" default sweep. Instead this holds every input at
/// its baseline value except one, nudges that one input by exactly its own standard deviation,
/// and measures the resulting shift in the target-plane impact point. The Monte Carlo sampler's
/// inputs are independent Gaussians and the ballistic response is smooth; at the magnitude of a
/// single sigma it is locally close to linear, so each source's one-sigma displacement `(dy, dz)`
/// is a standard first-order estimate of its variance contribution, and treating the sources as
/// independent gives `Var(total) ~= sum_i (dy_i^2 + dz_i^2)`. This costs a handful of
/// deterministic trajectory solves per range instead of a full extra Monte Carlo sub-run per
/// bucket.
///
/// Takes the (already solved) undispersed `baseline` position at `target_distance_m` from the
/// caller rather than solving it again -- the caller needs that same baseline to compute `p_hit`
/// (see `compute_wez`) and to decide whether attribution is meaningful at all when the
/// baseline doesn't reach this range.
///
/// Errors if any one-sigma perturbed solve fails validation (MBA-1317); the caller only invokes
/// this once the *un*perturbed baseline solve has already succeeded, but a sigma nudge can still
/// push an input (e.g. muzzle velocity, BC) out of `solve()`'s valid range.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the Monte Carlo sampler's own parameter set (MBA-1317)"
)]
fn wez_variance_shares(
    base_inputs: &BallisticInputs,
    base_wind: &WindConditions,
    atmosphere: &AtmosphericConditions,
    solver_max_range: f64,
    target_distance_m: f64,
    baseline: &Vector3<f64>,
    velocity_std_dev: f64,
    angle_std_dev_rad: f64,
    bc_std_dev: f64,
    azimuth_std_dev_rad: f64,
    wind_speed_std_dev: f64,
    wind_call_error_std_dev: f64,
    wind_direction_std_dev_rad: f64,
) -> Result<WezVarianceShares, BallisticsError> {
    // MV SD bucket: muzzle-velocity dispersion.
    let mv_sd_var = {
        let mut inputs = base_inputs.clone();
        inputs.muzzle_velocity = (inputs.muzzle_velocity + velocity_std_dev).max(0.0);
        wez_source_variance(
            velocity_std_dev,
            inputs,
            base_wind.clone(),
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?
    };

    // Other/group bucket: elevation, azimuth, and BC dispersion (mechanical/ammo "group"), plus
    // the ballistic (non-call) share of wind uncertainty.
    let mut other_var = 0.0;
    {
        let mut inputs = base_inputs.clone();
        inputs.muzzle_angle += angle_std_dev_rad;
        other_var += wez_source_variance(
            angle_std_dev_rad,
            inputs,
            base_wind.clone(),
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?;
    }
    {
        let mut inputs = base_inputs.clone();
        inputs.bc_value = (inputs.bc_value + bc_std_dev).max(0.01);
        other_var += wez_source_variance(
            bc_std_dev,
            inputs,
            base_wind.clone(),
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?;
    }
    {
        let mut inputs = base_inputs.clone();
        inputs.azimuth_angle += azimuth_std_dev_rad;
        other_var += wez_source_variance(
            azimuth_std_dev_rad,
            inputs,
            base_wind.clone(),
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?;
    }
    {
        let mut wind = base_wind.clone();
        wind.direction += wind_direction_std_dev_rad;
        other_var += wez_source_variance(
            wind_direction_std_dev_rad,
            base_inputs.clone(),
            wind,
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?;
    }
    {
        let mut wind = base_wind.clone();
        wind.speed += wind_speed_std_dev;
        other_var += wez_source_variance(
            wind_speed_std_dev,
            base_inputs.clone(),
            wind,
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?;
    }

    // Wind-call bucket: the shooter's own wind-speed estimation error, kept separate from the
    // ballistic wind-speed uncertainty above even though both perturb the same physical channel.
    let wind_call_var = {
        let mut wind = base_wind.clone();
        wind.speed += wind_call_error_std_dev;
        wez_source_variance(
            wind_call_error_std_dev,
            base_inputs.clone(),
            wind,
            atmosphere,
            solver_max_range,
            target_distance_m,
            baseline,
        )?
    };

    let total = wind_call_var + mv_sd_var + other_var;
    if total.is_nan() || total <= 0.0 {
        return Ok(WezVarianceShares::default());
    }
    Ok(WezVarianceShares {
        wind_call: wind_call_var / total,
        mv_sd: mv_sd_var / total,
        other: other_var / total,
    })
}

/// WEZ hit probability: the fraction of `results`' samples whose ABSOLUTE target-plane position
/// -- reconstructed as `baseline + (that sample's deviation from baseline)`, since
/// [`MonteCarloResults::impact_positions`] stores only the deviation -- falls
/// within `target_size`, centered on the fixed line of sight (`line_of_sight_height_m` vertically,
/// `z = 0` laterally).
///
/// This is deliberately NOT [`MonteCarloResults::hit_probability`] /
/// `rect_hit_probability`, which measure the miss distance from that SAME range's own baseline
/// (i.e. assume the shooter re-dials elevation perfectly for every range). A WEZ sweep instead
/// answers "how far can I hit this target size with ONE hold", so it must also count the
/// systematic ballistic drop below the fixed line of sight as a source of misses, not just random
/// dispersion -- see the module doc comment above `compute_wez`.
///
/// A sample that never reached the target plane keeps `MonteCarloResults`'s sentinel deviation
/// (`TARGET_NOT_REACHED_SENTINEL_M`, roughly -1e9 m). Added to any finite baseline that stays a
/// miss by a vast margin, so it is correctly excluded here without a separate check.
fn wez_p_hit(
    results: &MonteCarloResults,
    baseline: &Vector3<f64>,
    line_of_sight_height_m: f64,
    target_size: TargetSizeMetric,
) -> f64 {
    if results.impact_positions.is_empty() {
        return 0.0;
    }
    let hits = results
        .impact_positions
        .iter()
        .filter(|deviation| {
            let absolute_y = baseline.y + deviation.y;
            let absolute_z = baseline.z + deviation.z;
            let drop_from_los = absolute_y - line_of_sight_height_m;
            match target_size {
                TargetSizeMetric::Rect { width_m, height_m } => {
                    drop_from_los.abs() <= height_m / 2.0 && absolute_z.abs() <= width_m / 2.0
                }
                TargetSizeMetric::Radius { radius_m } => {
                    (drop_from_los * drop_from_los + absolute_z * absolute_z).sqrt() <= radius_m
                }
            }
        })
        .count();
    hits as f64 / results.impact_positions.len() as f64
}

/// One range step of a WEZ sweep.
#[derive(Debug, Clone, Serialize)]
pub struct WezRow {
    pub range_m: f64,
    pub p_hit: f64,
    pub dominant_error_source: Option<WezErrorBucket>,
    pub wind_call_share: f64,
    pub mv_sd_share: f64,
    pub other_share: f64,
    /// The undispersed baseline trajectory did not reach this range, so `*_share` and
    /// `dominant_error_source` above are not meaningful (left at their zero/`None` default).
    /// `p_hit` is unaffected -- it comes from the fully-dispersed Monte Carlo run directly.
    pub attribution_unavailable: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct WezTargetSizeJson {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub width_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub height_m: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radius_m: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct WezResult {
    pub target_size: WezTargetSizeJson,
    pub wind_speed_std_mps: f64,
    pub wind_call_error_mps: f64,
    /// `sqrt(wind_speed_std_mps^2 + wind_call_error_mps^2)`: the effective wind-speed standard
    /// deviation actually fed to the underlying Monte Carlo sampler at each range step.
    pub combined_wind_speed_std_mps: f64,
    pub num_sims_per_step: usize,
    pub rows: Vec<WezRow>,
}

/// Run a WEZ sweep and return its per-range rows plus the sweep-level parameters as a
/// [`WezResult`], leaving all rendering (summary table / statistics CSV / full JSON) to the
/// caller. All inputs are metric (the CLI converts from user units before calling).
///
/// # Parameters
///
/// NOTE the angle-like parameters (`angle`, `cant`, `wind_direction`, `angle_std`,
/// `wind_direction_std`) are in DEGREES — this function converts to radians itself,
/// unlike [`BallisticInputs`]/[`WindConditions`] elsewhere in the crate, which carry
/// radians. This mirrors the CLI flag set the sweep was extracted from (MBA-1317).
///
/// * `velocity` — muzzle velocity, m/s.
/// * `angle` — launch (elevation) angle held for every sweep step, DEGREES.
/// * `bc` — ballistic coefficient (dimensionless; referenced to `drag_model`).
/// * `mass` — bullet mass, kg.
/// * `diameter` — bullet diameter, m.
/// * `num_sims` — Monte Carlo samples per range step.
/// * `velocity_std` — muzzle-velocity standard deviation, m/s.
/// * `angle_std` — elevation-angle standard deviation, DEGREES (the derived
///   azimuth dispersion is half of it, matching the base Monte Carlo command).
/// * `bc_std` — BC standard deviation (dimensionless).
/// * `wind_std` — ballistic (gust-to-gust) wind-speed standard deviation, m/s.
/// * `wind_direction_std` — wind-direction standard deviation, DEGREES.
/// * `wind_speed` — base wind speed, m/s.
/// * `wind_direction` — base wind direction, DEGREES (wind-FROM: 0 = headwind,
///   90 = from the right).
/// * `wind_vertical` — base vertical wind, m/s, positive = updraft.
/// * `wind_call_error` — the shooter's wind-CALL error, m/s; composed with
///   `wind_std` in quadrature (see [`WezResult::combined_wind_speed_std_mps`]).
/// * `target_size` — the target box/radius, already in meters
///   ([`TargetSize::to_metric`]).
/// * `wez_start` / `wez_end` / `wez_step` — sweep bounds and step, meters
///   (`wez_end` inclusive).
/// * `drag_model` — the G-model `bc` is referenced to ([`DragModel::G1`] /
///   [`DragModel::G7`]); ignored for drag whenever `custom_drag_table` is set.
/// * `custom_drag_table` — optional Mach-keyed Cd deck replacing the G-model +
///   BC drag entirely.
/// * `cd_scale` — whole-curve multiplier on `custom_drag_table`'s interpolated Cd (MBA-1356);
///   `1.0` = neutral. Inert when `custom_drag_table` is `None`.
/// * `cant` — rifle cant, DEGREES, positive = clockwise from the shooter.
/// * `sight_offset_lateral_m` — lateral sight-to-bore mount offset, METERS, positive =
///   sight right of bore (MBA-1396); displaces every sample's initial lateral position
///   exactly like the trajectory command. `0.0` = neutral (byte-identical).
///
/// A fixed, distinct seed per range step (`0x57_45_5A_00 ^ step_index`) keeps a sweep
/// reproducible run-to-run while still drawing independent samples at each range.
#[allow(
    clippy::too_many_arguments,
    reason = "flat arguments mirror the stable Monte Carlo CLI command shape (MBA-1317)"
)]
pub fn compute_wez(
    velocity: f64,
    angle: f64,
    bc: f64,
    mass: f64,
    diameter: f64,
    num_sims: usize,
    velocity_std: f64,
    angle_std: f64,
    bc_std: f64,
    wind_std: f64,
    wind_direction_std: f64,
    wind_speed: f64,
    wind_direction: f64,
    wind_vertical: f64,
    wind_call_error: f64,
    target_size: TargetSizeMetric,
    wez_start: f64,
    wez_end: f64,
    wez_step: f64,
    drag_model: DragModel,
    custom_drag_table: Option<DragTable>,
    cd_scale: f64,
    cant: f64,
    sight_offset_lateral_m: f64,
) -> Result<WezResult, Box<dyn Error>> {
    if !(wez_step > 0.0 && wez_step.is_finite()) {
        return Err("--wez-step must be a positive, finite distance".into());
    }
    if !wez_start.is_finite() || !wez_end.is_finite() || wez_end < wez_start {
        return Err("--wez-end must be finite and >= --wez-start".into());
    }

    // Same bore-height/ground convention as the base `monte-carlo` command (MBA-967).
    let bore_height_metric = 1.5_f64;
    let base_inputs = BallisticInputs {
        muzzle_velocity: velocity,
        muzzle_angle: angle.to_radians(),
        bc_value: bc,
        bc_type: drag_model,
        bullet_mass: mass,
        bullet_diameter: diameter,
        muzzle_height: bore_height_metric,
        ground_threshold: 0.0,
        custom_drag_table,
        cd_scale,
        cant_angle: cant.to_radians(),
        sight_offset_lateral_m,
        ..Default::default()
    };
    let base_wind = WindConditions {
        speed: wind_speed,
        direction: wind_direction.to_radians(),
        vertical_speed: wind_vertical,
    };

    // The shooter's wind-call error is a dispersion source distinct from the ballistic
    // (gust-to-gust) wind-speed uncertainty --wind-std already models, but both perturb the same
    // physical channel (wind speed fed to the solve). As independent random errors they compose
    // in quadrature -- not by simple addition -- into the effective standard deviation the
    // underlying Monte Carlo sampler uses.
    let combined_wind_speed_std = wind_std.hypot(wind_call_error);

    // Matches run_monte_carlo's own convention: horizontal (azimuth) aim dispersion defaults to
    // half of the vertical (elevation) dispersion.
    let angle_std_rad = angle_std.to_radians();
    let azimuth_std_dev = angle_std_rad * 0.5;
    let wind_direction_std_rad = wind_direction_std.to_radians();

    // The fixed reference the WEZ target box is centered on: the horizontal line of sight,
    // extended straight (not the curved bullet path). This is what makes the zero-uncertainty
    // case a genuine step function -- ballistic drop below this fixed line, not just random
    // dispersion, can carry the bullet outside the box as range grows. It does NOT change per
    // range step: a WEZ sweep answers "how far can I engage this target size with ONE hold",
    // the classic point-blank-range question, not "assuming I re-dial for every range".
    let atmosphere = AtmosphericConditions {
        temperature: base_inputs.temperature,
        pressure: base_inputs.pressure,
        humidity: base_inputs.humidity_percent(),
        altitude: base_inputs.altitude,
    };
    let line_of_sight_height_m = base_inputs.muzzle_height + base_inputs.sight_height;

    let mut ranges_m = Vec::new();
    let mut next = wez_start;
    // Guard against an unbounded loop from a step so small that floating-point addition never
    // advances `next` past `wez_end`.
    for _ in 0..100_000 {
        if next > wez_end + wez_step * 1e-9 {
            break;
        }
        ranges_m.push(next);
        next += wez_step;
    }

    let mut rows = Vec::with_capacity(ranges_m.len());
    for (step_index, &range_m) in ranges_m.iter().enumerate() {
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let baseline = wez_solve_target_plane(
            base_inputs.clone(),
            base_wind.clone(),
            atmosphere.clone(),
            solver_max_range,
            range_m,
        )?;
        let baseline_reached = baseline.x >= range_m - 1e-6;

        let mc_params = MonteCarloParams {
            num_simulations: num_sims,
            velocity_std_dev: velocity_std,
            angle_std_dev: angle_std_rad,
            bc_std_dev: bc_std,
            wind_speed_std_dev: combined_wind_speed_std,
            target_distance: Some(range_m),
            base_wind_speed: wind_speed,
            base_wind_direction: wind_direction.to_radians(),
            azimuth_std_dev,
        };

        // A fixed, distinct seed per range step keeps a sweep reproducible run-to-run while
        // still drawing independent samples at each range.
        let seed = 0x57_45_5A_00_u64 ^ (step_index as u64);
        let p_hit = match crate::run_monte_carlo_with_wind_and_direction_std_dev_seeded(
            base_inputs.clone(),
            base_wind.clone(),
            mc_params,
            wind_direction_std_rad,
            seed,
        ) {
            Ok(results) => {
                wez_p_hit(&results, &baseline, line_of_sight_height_m, target_size)
            }
            // The baseline never reached this range plane at all -> every sample is a definite
            // miss for it.
            Err(_) => 0.0,
        };

        let (shares, attribution_unavailable) = if baseline_reached {
            (
                wez_variance_shares(
                    &base_inputs,
                    &base_wind,
                    &atmosphere,
                    solver_max_range,
                    range_m,
                    &baseline,
                    velocity_std,
                    angle_std_rad,
                    bc_std,
                    azimuth_std_dev,
                    wind_std,
                    wind_call_error,
                    wind_direction_std_rad,
                )?,
                false,
            )
        } else {
            (WezVarianceShares::default(), true)
        };

        rows.push(WezRow {
            range_m,
            p_hit,
            dominant_error_source: shares.dominant(),
            wind_call_share: shares.wind_call,
            mv_sd_share: shares.mv_sd,
            other_share: shares.other,
            attribution_unavailable,
        });
    }

    Ok(WezResult {
        target_size: match target_size {
            TargetSizeMetric::Rect { width_m, height_m } => WezTargetSizeJson {
                width_m: Some(width_m),
                height_m: Some(height_m),
                radius_m: None,
            },
            TargetSizeMetric::Radius { radius_m } => WezTargetSizeJson {
                width_m: None,
                height_m: None,
                radius_m: Some(radius_m),
            },
        },
        wind_speed_std_mps: wind_std,
        wind_call_error_mps: wind_call_error,
        combined_wind_speed_std_mps: combined_wind_speed_std,
        num_sims_per_step: num_sims,
        rows,
    })
}

#[cfg(test)]
mod wez_tests {
    use super::*;

    // A modest .308/168gr load, zeroed at 300 m with a shallow elevation that keeps the
    // trajectory well above ground for the whole 50-600 m range these tests sweep -- chosen with
    // `ballistics zero` (see CLI_USAGE.md's WEZ worked example for the imperial equivalent).
    fn test_base_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 823.0, // ~2700 fps
            muzzle_angle: 0.001274, // ~0.073 degrees: a 300 m zero for this load
            bc_value: 0.475,
            bullet_mass: 0.010_886, // 168 gr
            bullet_diameter: 0.007_82, // .308 in
            muzzle_height: 1.5,
            ground_threshold: 0.0,
            ..Default::default()
        }
    }

    fn test_atmosphere(inputs: &BallisticInputs) -> AtmosphericConditions {
        AtmosphericConditions {
            temperature: inputs.temperature,
            pressure: inputs.pressure,
            humidity: inputs.humidity_percent(),
            altitude: inputs.altitude,
        }
    }

    // ---- parse_target_size --------------------------------------------------------------

    #[test]
    fn parse_target_size_accepts_a_wxh_rectangle() {
        assert_eq!(
            parse_target_size("18x30").unwrap(),
            TargetSize::Rect {
                width: 18.0,
                height: 30.0
            }
        );
        // Case-insensitive separator and surrounding whitespace.
        assert_eq!(
            parse_target_size(" 18.5X30.25 ").unwrap(),
            TargetSize::Rect {
                width: 18.5,
                height: 30.25
            }
        );
    }

    #[test]
    fn parse_target_size_accepts_a_single_radius() {
        assert_eq!(parse_target_size("12").unwrap(), TargetSize::Radius(12.0));
        assert_eq!(parse_target_size(" 0.5 ").unwrap(), TargetSize::Radius(0.5));
    }

    #[test]
    fn parse_target_size_rejects_garbage() {
        for bad in [
            "",
            "   ",
            "abc",
            "18xthirty",
            "eighteenx30",
            "18x30x40",
            "0",
            "-5",
            "18x-5",
            "18x0",
            "NaN",
        ] {
            assert!(
                parse_target_size(bad).is_err(),
                "expected an error for {bad:?}"
            );
        }
    }

    // ---- WEZ hit-probability step function -----------------------------------------------

    #[test]
    fn zero_uncertainty_is_a_step_function_in_range() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let atmosphere = test_atmosphere(&inputs);
        // 18x30 box: 0.4572 m x 0.762 m.
        let target = TargetSizeMetric::Rect {
            width_m: 0.4572,
            height_m: 0.762,
        };
        let los_height_m = inputs.muzzle_height + inputs.sight_height;

        let mc_params = MonteCarloParams {
            num_simulations: 20,
            velocity_std_dev: 0.0,
            angle_std_dev: 0.0,
            bc_std_dev: 0.0,
            wind_speed_std_dev: 0.0,
            target_distance: None,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0,
        };

        let mut p_hits = Vec::new();
        for &range_m in &[50.0_f64, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0] {
            let solver_max_range = range_m.max(1000.0) * 2.0;
            let baseline = wez_solve_target_plane(
                inputs.clone(),
                wind.clone(),
                atmosphere.clone(),
                solver_max_range,
                range_m,
            )
            .expect("valid test baseline solve");
            let mut params = mc_params.clone();
            params.target_distance = Some(range_m);
            let results = crate::run_monte_carlo_with_wind_and_direction_std_dev_seeded(
                inputs.clone(),
                wind.clone(),
                params,
                0.0,
                0xA11CE,
            )
            .expect("zero-uncertainty solve");
            let p_hit = wez_p_hit(&results, &baseline, los_height_m, target);
            // Every one of the (identical, undispersed) samples must agree: exactly a hit or
            // exactly a miss, never a fractional probability.
            assert!(
                p_hit == 0.0 || p_hit == 1.0,
                "range {range_m} m: expected a step (0.0 or 1.0), got {p_hit}"
            );
            p_hits.push((range_m, p_hit));
        }

        assert!(
            p_hits.iter().any(|&(_, p)| p == 1.0),
            "expected at least one in-box range close to the muzzle: {p_hits:?}"
        );
        assert!(
            p_hits.iter().any(|&(_, p)| p == 0.0),
            "expected at least one out-of-box range far downrange: {p_hits:?}"
        );
        // Once it steps down to a miss, a plain (unheld) trajectory that has already passed its
        // zero does not come back into a fixed-size box further downrange.
        let first_miss = p_hits.iter().position(|&(_, p)| p == 0.0);
        if let Some(idx) = first_miss {
            assert!(
                p_hits[idx..].iter().all(|&(_, p)| p == 0.0),
                "expected the box exit to be permanent for the rest of the sweep: {p_hits:?}"
            );
        }
    }

    // ---- P(hit) monotonicity with real dispersion -----------------------------------------

    #[test]
    fn p_hit_is_monotone_non_increasing_with_range() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let atmosphere = test_atmosphere(&inputs);
        let target = TargetSizeMetric::Rect {
            width_m: 0.4572,
            height_m: 0.762,
        };
        let los_height_m = inputs.muzzle_height + inputs.sight_height;
        let wind_call_error = 1.5_f64; // m/s
        let wind_std = 0.5_f64; // m/s
        let combined_wind_std = wind_std.hypot(wind_call_error);

        let mc_params = MonteCarloParams {
            num_simulations: 500, // a fixed seed keeps this run-to-run deterministic
            velocity_std_dev: 1.0,
            angle_std_dev: 0.001,
            bc_std_dev: 0.01,
            wind_speed_std_dev: combined_wind_std,
            target_distance: None,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.0005,
        };

        let ranges_m = [100.0_f64, 200.0, 300.0, 400.0, 500.0, 600.0];
        let mut p_hits = Vec::new();
        for (step_index, &range_m) in ranges_m.iter().enumerate() {
            let solver_max_range = range_m.max(1000.0) * 2.0;
            let baseline = wez_solve_target_plane(
                inputs.clone(),
                wind.clone(),
                atmosphere.clone(),
                solver_max_range,
                range_m,
            )
            .expect("valid test baseline solve");
            let mut params = mc_params.clone();
            params.target_distance = Some(range_m);
            let seed = 0x57_45_5A_00_u64 ^ (step_index as u64);
            let results = crate::run_monte_carlo_with_wind_and_direction_std_dev_seeded(
                inputs.clone(),
                wind.clone(),
                params,
                0.0,
                seed,
            )
            .expect("dispersed solve");
            p_hits.push(wez_p_hit(&results, &baseline, los_height_m, target));
        }

        // A large sample count plus a fixed seed makes this close to the noiseless limit, but a
        // finite Monte Carlo estimate can still tick up by a hair at the boundary between two
        // adjacent steps -- allow a small generous tolerance rather than asserting exact
        // non-increase (MBA-1317 test spec).
        let tolerance = 0.03;
        for pair in p_hits.windows(2) {
            assert!(
                pair[1] <= pair[0] + tolerance,
                "P(hit) rose more than the allowed jitter: {p_hits:?}"
            );
        }
        // The overall trend across the full sweep must be a clear decline.
        assert!(
            p_hits.first().unwrap() - p_hits.last().unwrap() > 0.2,
            "expected a clear overall decline across the sweep: {p_hits:?}"
        );
    }

    // ---- Variance-attribution shares --------------------------------------------------------

    #[test]
    fn variance_shares_sum_to_one_when_multiple_sources_are_active() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let atmosphere = test_atmosphere(&inputs);
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let baseline = wez_solve_target_plane(
            inputs.clone(),
            wind.clone(),
            atmosphere.clone(),
            solver_max_range,
            range_m,
        )
        .expect("valid test baseline solve");

        let shares = wez_variance_shares(
            &inputs,
            &wind,
            &atmosphere,
            solver_max_range,
            range_m,
            &baseline,
            /* velocity_std_dev */ 1.0,
            /* angle_std_dev_rad */ 0.001,
            /* bc_std_dev */ 0.01,
            /* azimuth_std_dev_rad */ 0.0005,
            /* wind_speed_std_dev */ 0.4,
            /* wind_call_error_std_dev */ 1.2,
            /* wind_direction_std_dev_rad */ 0.02,
        )
        .expect("valid test attribution solve");

        let sum = shares.wind_call + shares.mv_sd + shares.other;
        assert!(
            (sum - 1.0).abs() < 1e-9,
            "shares should sum to ~1.0, got {sum} ({shares:?})"
        );
        for share in [shares.wind_call, shares.mv_sd, shares.other] {
            assert!((0.0..=1.0).contains(&share), "share out of range: {share}");
        }
        assert!(shares.dominant().is_some());
    }

    #[test]
    fn variance_shares_are_all_zero_with_no_dispersion_sources() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let atmosphere = test_atmosphere(&inputs);
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let baseline = wez_solve_target_plane(
            inputs.clone(),
            wind.clone(),
            atmosphere.clone(),
            solver_max_range,
            range_m,
        )
        .expect("valid test baseline solve");

        let shares = wez_variance_shares(
            &inputs,
            &wind,
            &atmosphere,
            solver_max_range,
            range_m,
            &baseline,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
        .expect("valid test attribution solve");

        assert_eq!(shares.wind_call, 0.0);
        assert_eq!(shares.mv_sd, 0.0);
        assert_eq!(shares.other, 0.0);
        assert!(shares.dominant().is_none());
    }

    #[test]
    fn wind_call_bucket_dominates_when_it_is_the_only_active_source() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let atmosphere = test_atmosphere(&inputs);
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let baseline = wez_solve_target_plane(
            inputs.clone(),
            wind.clone(),
            atmosphere.clone(),
            solver_max_range,
            range_m,
        )
        .expect("valid test baseline solve");

        let shares = wez_variance_shares(
            &inputs,
            &wind,
            &atmosphere,
            solver_max_range,
            range_m,
            &baseline,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            /* wind_call_error_std_dev */ 3.0,
            0.0,
        )
        .expect("valid test attribution solve");

        assert!((shares.wind_call - 1.0).abs() < 1e-9);
        assert_eq!(shares.mv_sd, 0.0);
        assert_eq!(shares.other, 0.0);
        assert_eq!(shares.dominant(), Some(WezErrorBucket::WindCall));
    }
}
