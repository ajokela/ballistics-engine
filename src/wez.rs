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
use crate::perturbation::{central_difference, InputAxis, KernelError};
use crate::solve_json::{
    DragModelV1, ResolvedAtmosphereV1, ResolvedConstantWindV1, ResolvedEffectsV1,
    ResolvedProjectileV1, ResolvedRifleV1, ResolvedSamplingV1, ResolvedShotV1,
    ResolvedSolveRequestV1, ResolvedSolverV1, ResolvedWindV1, SchemaVersionV1, SolverMethodV1,
    TwistDirectionV1,
};
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

/// Map the engine's own `DragModel` to solve-json v1's `DragModelV1`, when possible.
///
/// `None` when `model` has no v1 counterpart: v1 only has `G1`/`G6`/`G7`/`G8`
/// ([`crate::solve_json::DragModelV1`]), while the engine's own enum also has `G2`, `G5`, `GI`,
/// `GS`, and `RA4`. An exhaustive match with no wildcard arm, so a future engine `DragModel`
/// variant fails to compile here until it is explicitly placed on one side or the other, rather
/// than silently falling through to `None`.
fn wez_kernel_drag_model(model: DragModel) -> Option<DragModelV1> {
    match model {
        DragModel::G1 => Some(DragModelV1::G1),
        DragModel::G6 => Some(DragModelV1::G6),
        DragModel::G7 => Some(DragModelV1::G7),
        DragModel::G8 => Some(DragModelV1::G8),
        DragModel::G2 | DragModel::G5 | DragModel::GI | DragModel::GS | DragModel::RA4 => None,
    }
}

/// Build the [`ResolvedSolveRequestV1`] the shared perturbation kernel needs to attribute
/// variance for one WEZ range step (0.33.0 decision-support D2), mirroring `base_inputs`/
/// `base_wind`'s already-fully-resolved values field for field -- same units, same defaults --
/// rather than re-deriving them through solve-json's own defaulting a second time. Every field
/// `compute_wez` does not expose as a sweep parameter is left at exactly the value
/// `BallisticInputs::default()`/`WindConditions::default()` already carry, which is also the
/// value `solve_v1::prepare_request`'s own hardcoded mappings for engine-only fields (powder
/// sensitivity, tipoff, wind shear, pitch damping, ...) always resolve to regardless of the
/// request -- so the two paths agree everywhere the kernel has no way to differ, and this
/// function only needs to carry the fields WEZ actually threads through.
///
/// `solver_max_range` is threaded through separately rather than read off `base_inputs` because
/// it is not a `BallisticInputs` field at all: `compute_wez` derives it per range step
/// (`range_m.max(1000.0) * 2.0`) and applies it via `TrajectorySolver::set_max_range`, so the
/// resolved request's `shot.max_range_m` must match that same per-step value for the kernel's
/// perturbed solves to see the identical trajectory truncation the row's own baseline solve
/// used.
///
/// Returns `None` when this configuration cannot be represented on the solve-json v1 wire
/// contract at all: a loaded custom drag table (`base_inputs.custom_drag_table`, which replaces
/// G-model+BC drag entirely -- v1 has no field for a custom deck at all), or an engine
/// `DragModel` v1 has no variant for ([`wez_kernel_drag_model`]). Differentiating under the
/// wrong drag physics would silently misattribute variance rather than fail loudly, so this is
/// treated by the caller the same as a kernel structural refusal
/// (`WezRow::attribution_unavailable`), not answered with a plausible wrong number. Both are
/// real, CLI/WASM-reachable configurations (`monte-carlo --wez --drag-table ...` /
/// the WASM terminal's `--drag-model g2`/`gi`/`gs`/`g5`/`ra4` on `--wez`), not hypothetical.
fn wez_resolved_request(
    base_inputs: &BallisticInputs,
    base_wind: &WindConditions,
    solver_max_range: f64,
) -> Option<ResolvedSolveRequestV1> {
    if base_inputs.custom_drag_table.is_some() {
        return None;
    }
    let drag_model = wez_kernel_drag_model(base_inputs.bc_type)?;

    Some(ResolvedSolveRequestV1 {
        schema_version: SchemaVersionV1,
        projectile: ResolvedProjectileV1 {
            mass_kg: base_inputs.bullet_mass,
            diameter_m: base_inputs.bullet_diameter,
            length_m: Some(base_inputs.bullet_length),
            drag_model,
            ballistic_coefficient: base_inputs.bc_value,
        },
        rifle: ResolvedRifleV1 {
            muzzle_velocity_mps: base_inputs.muzzle_velocity,
            sight_height_m: base_inputs.sight_height,
            muzzle_height_m: base_inputs.muzzle_height,
            // BallisticInputs stores twist rate in inches/turn; solve-json v1 is SI.
            twist_rate_m_per_turn: base_inputs.twist_rate * 0.0254,
            twist_direction: if base_inputs.is_twist_right {
                TwistDirectionV1::Right
            } else {
                TwistDirectionV1::Left
            },
            sight_offset_lateral_m: Some(base_inputs.sight_offset_lateral_m),
        },
        shot: ResolvedShotV1 {
            max_range_m: solver_max_range,
            // compute_wez always supplies an explicit muzzle angle directly and never zeroes
            // (see its own doc comment) -- no re-zero search may run on any perturbed solve.
            zero_distance_m: None,
            muzzle_angle_rad: base_inputs.muzzle_angle,
            aim_azimuth_rad: base_inputs.azimuth_angle,
            shot_azimuth_rad: base_inputs.shot_azimuth,
            shooting_angle_rad: base_inputs.shooting_angle,
            cant_angle_rad: base_inputs.cant_angle,
            target_height_m: base_inputs.target_height,
            ground_threshold_m: base_inputs.ground_threshold,
            zero_poi_up_m: Some(base_inputs.zero_poi_vertical_m),
            zero_poi_right_m: Some(base_inputs.zero_poi_horizontal_m),
            drops_reference: None,
        },
        atmosphere: ResolvedAtmosphereV1 {
            altitude_m: base_inputs.altitude,
            temperature_k: base_inputs.temperature + 273.15,
            pressure_pa: base_inputs.pressure * 100.0,
            relative_humidity: base_inputs.humidity,
            latitude_rad: base_inputs.latitude.map(f64::to_radians),
            pressure_reference: None,
        },
        wind: ResolvedWindV1::Constant(ResolvedConstantWindV1 {
            speed_mps: base_wind.speed,
            direction_from_rad: base_wind.direction,
            vertical_speed_mps: base_wind.vertical_speed,
            wind_reference: None,
        }),
        solver: ResolvedSolverV1 {
            // BallisticInputs::default() runs adaptive RK45 (use_rk4 && use_adaptive_rk45,
            // both true, and compute_wez never overrides either), so the kernel's perturbed
            // solves must use the same method to match the row's own baseline solve.
            method: SolverMethodV1::Rk45,
            time_step_s: 0.001, // ignored by RK45; any positive finite value is inert here.
        },
        effects: ResolvedEffectsV1 {
            magnus: false,
            coriolis: false,
            enhanced_spin_drift: false,
        },
        // Unused by `perturbation::evaluate` (it queries a specific range off the solved
        // trajectory directly, never the regular sampling grid); any positive value is inert.
        sampling: ResolvedSamplingV1 { interval_m: 10.0 },
        reticle: None,
    })
}

/// Central-difference WEZ miss-variance attribution at a single range (0.33.0 decision-support
/// D2), replacing the retired `wez_source_variance`'s one-sided, one-sigma solves with the
/// shared perturbation kernel's central differences -- see the module doc's opening paragraph
/// for why the two numerical methods needed to converge onto one.
///
/// For each source with standard deviation `sigma` and kernel axis `axis`, the squared
/// one-sigma displacement is `(d.d_drop_d_x * sigma)^2 + (d.d_windage_d_x * sigma)^2`, where `d`
/// is [`central_difference`]`(base_resolved, axis, &[target_distance_m], None)`'s single
/// [`crate::perturbation::Derivative`] -- the same delta-method approximation
/// `crate::error_budget` uses, applied here to WEZ's three collapsed buckets instead of
/// `error_budget`'s per-source ranking. Treating the sources as independent gives
/// `Var(total) ~= sum_i displacement_i^2`, exactly as the retired one-sided version did; only
/// how each `displacement_i` is measured has changed. A source whose `sigma` is non-positive or
/// NaN short-circuits to a `0.0` contribution WITHOUT calling the kernel at all (mirroring
/// `wez_source_variance`'s own guard), so a genuinely-disabled source costs nothing rather than
/// one wasted pair of solves.
///
/// Bucket assignment is unchanged from the retired one-sided version: `MvSd` <- muzzle velocity
/// ([`InputAxis::MuzzleVelocityMps`]); `WindCall` <- the wind-call error channel
/// ([`InputAxis::WindSpeed`], `wind_call_error_std_dev`); `Other` <- the accumulated sum of
/// muzzle angle ([`InputAxis::MuzzleAngle`]), ballistic coefficient
/// ([`InputAxis::BallisticCoefficient`]), aim azimuth ([`InputAxis::AimAzimuth`] -- WEZ's
/// `azimuth_angle` is the small horizontal AIMING offset, not the compass-referenced
/// `InputAxis::ShotAzimuth` Coriolis uses), wind direction ([`InputAxis::WindDirection`]), and
/// the ballistic (non-call) share of wind speed ([`InputAxis::WindSpeed`] again,
/// `wind_speed_std_dev`).
///
/// `WindSpeed` is the ONE axis two different sources share (`wind_speed_std_dev` into `Other`,
/// `wind_call_error_std_dev` into `WindCall`), and unlike every other pair of sources here they
/// are not independently measured: the SAME [`Derivative`](crate::perturbation::Derivative) is
/// computed at most once (see the code below) and both displacements are scaled off it. A real,
/// worth-stating consequence (review finding, not part of the original D2 design intent): the
/// wind-call share and the ballistic-wind PORTION of `Other` are now in an EXACT
/// `(wind_call_error_std_dev / wind_speed_std_dev)^2` ratio whenever both are positive --
/// `displacement_call^2 / displacement_wind^2 = (d * sigma_call)^2 / (d * sigma_wind)^2 =
/// (sigma_call / sigma_wind)^2`, since the shared `d` cancels out of the ratio entirely, leaving
/// pure sigma-squared algebra with no ballistic response left in it. The retired one-sided
/// version did NOT have this property: it perturbed `wind.speed` twice, independently, at two
/// different absolute deltas (`+sigma_call` and `+sigma_wind`), so its two displacements were
/// two genuinely different finite differences, only APPROXIMATELY proportional to
/// `sigma^2` in the locally-linear regime -- see
/// `wind_call_and_ballistic_wind_shares_are_exactly_proportional_to_sigma_squared` in this
/// module's tests.
///
/// Returns `Ok(None)` when ANY of the seven sources hits a structural kernel refusal --
/// [`KernelError::AxisUnsupportedForRequest`], [`KernelError::AxisAbsent`],
/// [`KernelError::CategoricalAxis`], or [`KernelError::StepOutOfDomain`], classified via
/// `crate::error_budget::unavailable_reason` (the same four-way split `error_budget` and
/// `crate::tolerance::tolerance_envelope` already share) -- rather than silently normalizing
/// over whichever sources DID evaluate, which would misattribute the missing source's variance
/// to the other two buckets with no visible sign anything was skipped. `compute_wez` maps this
/// `None` onto [`WezRow::attribution_unavailable`], its existing "nothing to attribute here"
/// contract -- the same flag already used when the baseline does not reach this range. In
/// practice none of the seven axes above carries an `AxisUnsupportedForRequest` guard (those
/// are `Altitude` under QNH pressure and `ShotAzimuth` under compass wind; [`wez_resolved_request`]
/// never sets either reference mode) and the resolved wind is always `Constant` (never
/// `AxisAbsent`'s segmented-wind trigger), so this path is defensive rather than routinely hit --
/// see the task report for the specific fixtures that DO reach it. Short-circuits on the FIRST
/// such refusal rather than evaluating the remaining sources, since the row's attribution is
/// already going to be reported unavailable either way.
///
/// # Errors
/// Propagates any other [`KernelError`] (`Solve`, `Observation`, `TypeMismatch`, `NonFinite`) --
/// a genuine solver or trajectory failure, not a normal "this input cannot be perturbed here"
/// fact -- exactly as a failed perturbed solve propagated out of the retired one-sided
/// `wez_source_variance` (MBA-1317).
#[allow(
    clippy::too_many_arguments,
    reason = "flat sigma list mirrors the Monte Carlo sampler's own parameter set (MBA-1317)"
)]
fn wez_variance_shares(
    base_resolved: &ResolvedSolveRequestV1,
    target_distance_m: f64,
    velocity_std_dev: f64,
    angle_std_dev_rad: f64,
    bc_std_dev: f64,
    azimuth_std_dev_rad: f64,
    wind_speed_std_dev: f64,
    wind_call_error_std_dev: f64,
    wind_direction_std_dev_rad: f64,
) -> Result<Option<WezVarianceShares>, KernelError> {
    // One source's squared one-sigma displacement from an already-computed derivative, or
    // `0.0` when `sigma` is non-positive or NaN (that source is disabled) -- mirrors
    // `wez_source_variance`'s own guard, without needing the kernel call that produced `d` to
    // know anything about `sigma`.
    let displacement_sq = |d: &crate::perturbation::Derivative, sigma: f64| -> f64 {
        if sigma.is_nan() || sigma <= 0.0 {
            0.0
        } else {
            (d.d_drop_d_x * sigma).powi(2) + (d.d_windage_d_x * sigma).powi(2)
        }
    };

    // One axis's central-difference derivative, independent of any particular source's sigma.
    // `Ok(None)` signals a structural refusal the caller must treat as attribution-unavailable
    // rather than a genuine error (see the doc comment above).
    let derivative_for =
        |axis: InputAxis| -> Result<Option<crate::perturbation::Derivative>, KernelError> {
            match central_difference(base_resolved, axis, &[target_distance_m], None) {
                Ok(d) => Ok(Some(d[0])),
                Err(e) if crate::error_budget::unavailable_reason(&e).is_some() => Ok(None),
                Err(e) => Err(e),
            }
        };

    // A single source with sigma `sigma` and kernel axis `axis`: `Ok(Some(0.0))` without
    // calling the kernel at all when `sigma` is non-positive or NaN (that source is disabled),
    // so a genuinely-disabled source costs nothing rather than one wasted pair of solves.
    let contribution = |axis: InputAxis, sigma: f64| -> Result<Option<f64>, KernelError> {
        if sigma.is_nan() || sigma <= 0.0 {
            return Ok(Some(0.0));
        }
        Ok(derivative_for(axis)?.map(|d| displacement_sq(&d, sigma)))
    };

    // MV SD bucket: muzzle-velocity dispersion.
    let Some(mv_sd_var) = contribution(InputAxis::MuzzleVelocityMps, velocity_std_dev)? else {
        return Ok(None);
    };

    // Other/group bucket: elevation, aim azimuth, and BC dispersion (mechanical/ammo "group").
    // The ballistic (non-call) share of wind uncertainty is folded in just below, alongside the
    // WindCall bucket, since both need the SAME WindSpeed derivative (review finding M1).
    let mut other_var = 0.0;
    for (axis, sigma) in [
        (InputAxis::MuzzleAngle, angle_std_dev_rad),
        (InputAxis::BallisticCoefficient, bc_std_dev),
        (InputAxis::AimAzimuth, azimuth_std_dev_rad),
        (InputAxis::WindDirection, wind_direction_std_dev_rad),
    ] {
        let Some(v) = contribution(axis, sigma)? else {
            return Ok(None);
        };
        other_var += v;
    }

    // Wind speed and wind-call bucket: `wind_speed_std_dev` (ballistic, folded into `other_var`
    // above) and `wind_call_error_std_dev` (the WindCall bucket) perturb the SAME
    // `InputAxis::WindSpeed` channel -- kept as separate attribution buckets (the shooter's own
    // wind-call estimation error vs. physical gust-to-gust variability) but sharing one
    // derivative computed AT MOST ONCE (review finding M1: this used to call the kernel twice
    // with identical arguments), only when at least one of the two sigmas is actually active.
    // See the doc comment above for the exact-proportionality consequence this has.
    let wind_call_var = if (wind_speed_std_dev.is_nan() || wind_speed_std_dev <= 0.0)
        && (wind_call_error_std_dev.is_nan() || wind_call_error_std_dev <= 0.0)
    {
        0.0
    } else {
        let Some(d) = derivative_for(InputAxis::WindSpeed)? else {
            return Ok(None);
        };
        other_var += displacement_sq(&d, wind_speed_std_dev);
        displacement_sq(&d, wind_call_error_std_dev)
    };

    let total = wind_call_var + mv_sd_var + other_var;
    if total.is_nan() || total <= 0.0 {
        return Ok(Some(WezVarianceShares::default()));
    }
    Ok(Some(WezVarianceShares {
        wind_call: wind_call_var / total,
        mv_sd: mv_sd_var / total,
        other: other_var / total,
    }))
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
    /// `*_share` and `dominant_error_source` above are not meaningful (left at their
    /// zero/`None` default) for any of THREE reasons, and this flag does not distinguish them:
    /// (1) the undispersed baseline trajectory did not reach this range; (2) central-difference
    /// attribution (0.33.0 decision-support D2) hit a structural kernel refusal on one of its
    /// seven sources (`crate::perturbation::KernelError::AxisUnsupportedForRequest`/
    /// `AxisAbsent`/`CategoricalAxis`/`StepOutOfDomain`); or (3) this configuration cannot be
    /// represented on the shared kernel's solve-json v1 wire contract at all -- a loaded custom
    /// drag table, or an engine drag model solve-json v1 has no variant for (`G2`/`G5`/`GI`/
    /// `GS`/`RA4`; see `wez_resolved_request`). Reason (3) is the one most likely to surprise a
    /// caller: every row of a `--drag-table` or unsupported-`--drag-model` sweep reads `n/a`
    /// here, which is NOT a claim that the bullet fails to reach that range. `p_hit` is
    /// unaffected in every case -- it comes from the fully-dispersed Monte Carlo run directly,
    /// never from the kernel.
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
            match wez_resolved_request(&base_inputs, &base_wind, solver_max_range) {
                Some(resolved) => match wez_variance_shares(
                    &resolved,
                    range_m,
                    velocity_std,
                    angle_std_rad,
                    bc_std,
                    azimuth_std_dev,
                    wind_std,
                    wind_call_error,
                    wind_direction_std_rad,
                )? {
                    Some(shares) => (shares, false),
                    // A structural kernel refusal (see wez_variance_shares's doc) -- treated
                    // the same as the baseline-not-reached case just below: nothing to
                    // attribute, but p_hit is unaffected.
                    None => (WezVarianceShares::default(), true),
                },
                // This configuration cannot be represented on the solve-json v1 wire contract
                // at all (a custom drag table, or a DragModel v1 has no variant for) -- see
                // wez_resolved_request's doc.
                None => (WezVarianceShares::default(), true),
            }
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
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let resolved = wez_resolved_request(&inputs, &wind, solver_max_range)
            .expect("test_base_inputs's default G1 drag model is always kernel-representable");

        let shares = wez_variance_shares(
            &resolved,
            range_m,
            /* velocity_std_dev */ 1.0,
            /* angle_std_dev_rad */ 0.001,
            /* bc_std_dev */ 0.01,
            /* azimuth_std_dev_rad */ 0.0005,
            /* wind_speed_std_dev */ 0.4,
            /* wind_call_error_std_dev */ 1.2,
            /* wind_direction_std_dev_rad */ 0.02,
        )
        .expect("valid test attribution solve")
        .expect("attribution must be available for this fixture");

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
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let resolved = wez_resolved_request(&inputs, &wind, solver_max_range)
            .expect("test_base_inputs's default G1 drag model is always kernel-representable");

        let shares = wez_variance_shares(
            &resolved, range_m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        )
        .expect("valid test attribution solve")
        .expect("attribution must be available for this fixture");

        assert_eq!(shares.wind_call, 0.0);
        assert_eq!(shares.mv_sd, 0.0);
        assert_eq!(shares.other, 0.0);
        assert!(shares.dominant().is_none());
    }

    #[test]
    fn wind_call_bucket_dominates_when_it_is_the_only_active_source() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let resolved = wez_resolved_request(&inputs, &wind, solver_max_range)
            .expect("test_base_inputs's default G1 drag model is always kernel-representable");

        let shares = wez_variance_shares(
            &resolved,
            range_m,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            /* wind_call_error_std_dev */ 3.0,
            0.0,
        )
        .expect("valid test attribution solve")
        .expect("attribution must be available for this fixture");

        assert!((shares.wind_call - 1.0).abs() < 1e-9);
        assert_eq!(shares.mv_sd, 0.0);
        assert_eq!(shares.other, 0.0);
        assert_eq!(shares.dominant(), Some(WezErrorBucket::WindCall));
    }

    /// Review finding M1: `WindSpeed`'s derivative is now computed AT MOST ONCE and reused for
    /// both `wind_speed_std_dev` (folded into `Other`) and `wind_call_error_std_dev`
    /// (`WindCall`). With every OTHER source's sigma at zero, `other_share` is ENTIRELY the
    /// ballistic wind-speed contribution, so the ratio between the two buckets isolates exactly
    /// `(wind_call_error_std_dev / wind_speed_std_dev)^2` -- the shared derivative cancels out
    /// of the ratio entirely, leaving pure sigma-squared algebra with no ballistic response left
    /// in it. This is a real, stated behavior change from the retired one-sided version, which
    /// perturbed `wind.speed` independently at two different absolute deltas (`+sigma_call` and
    /// `+sigma_wind`) and so was only approximately proportional to `sigma^2` in the
    /// locally-linear regime, not exactly.
    #[test]
    fn wind_call_and_ballistic_wind_shares_are_exactly_proportional_to_sigma_squared() {
        let inputs = test_base_inputs();
        let wind = WindConditions::default();
        let range_m: f64 = 300.0;
        let solver_max_range = range_m.max(1000.0) * 2.0;
        let resolved = wez_resolved_request(&inputs, &wind, solver_max_range)
            .expect("test_base_inputs's default G1 drag model is always kernel-representable");

        let wind_speed_std_dev = 0.4_f64;
        let wind_call_error_std_dev = 1.2_f64;
        let shares = wez_variance_shares(
            &resolved,
            range_m,
            0.0,
            0.0,
            0.0,
            0.0,
            wind_speed_std_dev,
            wind_call_error_std_dev,
            0.0,
        )
        .expect("valid test attribution solve")
        .expect("attribution must be available for this fixture");

        assert!(shares.other > 0.0, "fixture must produce a nonzero ballistic-wind share");
        let expected_ratio = (wind_call_error_std_dev / wind_speed_std_dev).powi(2);
        let actual_ratio = shares.wind_call / shares.other;
        assert!(
            (actual_ratio - expected_ratio).abs() < 1e-9,
            "expected wind_call/other == (sigma_call/sigma_wind)^2 == {expected_ratio}, got \
             {actual_ratio}"
        );
    }

    /// `wez_resolved_request`'s output must describe the SAME physical trajectory
    /// `wez_solve_target_plane` computes directly from `BallisticInputs` -- otherwise
    /// central-difference attribution would measure sensitivity against a subtly different
    /// configuration than the one WEZ actually simulates for `p_hit`. Independent oracle:
    /// `TrajectoryObservation::drop_m` is LOS-perpendicular and positive BELOW the line of
    /// sight (`line_of_sight_height_m - position.y`, `trajectory_observation.rs`);
    /// `wez_solve_target_plane`'s raw `position.y` is world-frame height, so
    /// `line_of_sight_height_m - baseline.y` is the independently-derived equivalent.
    /// `windage_m` is `position.z` directly (no sign flip, no scale change), so it compares to
    /// `baseline.z` with no transform at all.
    ///
    /// M3 review fix: `cant_angle`, `sight_offset_lateral_m`, and `azimuth_angle` are
    /// deliberately nonzero here (the first two are directly user-settable on `monte-carlo
    /// --wez`, `--cant`/`--sight-offset`; `azimuth_angle` is not exposed as a baseline there but
    /// `wez_resolved_request` must still map it correctly for any `BallisticInputs`). With all
    /// three left at `test_base_inputs()`'s shared 0.0 default, a mis-mapping -- e.g.
    /// `cant_angle_rad` accidentally written to `shooting_angle_rad`, a DIFFERENT `ResolvedShotV1`
    /// field -- would leave every assertion in this file green anyway, since both fields agree
    /// at the neutral value; nonzero values are the only way this oracle can actually catch that
    /// class of bug.
    #[test]
    fn resolved_request_matches_wez_solve_target_plane_baseline() {
        let mut inputs = test_base_inputs();
        inputs.cant_angle = 5.0_f64.to_radians();
        inputs.sight_offset_lateral_m = 0.02;
        inputs.azimuth_angle = 0.001;
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

        let resolved = wez_resolved_request(&inputs, &wind, solver_max_range)
            .expect("test_base_inputs's default G1 drag model is always kernel-representable");
        let req: crate::solve_json::SolveRequestV1 = (&resolved).into();
        let obs = crate::perturbation::evaluate(&req, &[range_m]).expect("kernel evaluate");
        assert_eq!(obs.len(), 1);

        let line_of_sight_height_m = inputs.muzzle_height + inputs.sight_height;
        let expected_drop_m = line_of_sight_height_m - baseline.y;
        assert!(
            (obs[0].drop_m - expected_drop_m).abs() < 1e-6,
            "kernel drop_m {} disagrees with wez_solve_target_plane-derived {} -- \
             wez_resolved_request likely mismaps a field",
            obs[0].drop_m,
            expected_drop_m
        );
        assert!(
            (obs[0].windage_m - baseline.z).abs() < 1e-6,
            "kernel windage_m {} disagrees with wez_solve_target_plane baseline.z {} -- \
             wez_resolved_request likely mismaps a field",
            obs[0].windage_m,
            baseline.z
        );
    }

    /// `DragModel::G2`/`G5`/`GI`/`GS`/`RA4` have no solve-json v1 counterpart (v1 only has
    /// `G1`/`G6`/`G7`/`G8`) -- differentiating under the wrong drag physics would silently
    /// misattribute variance, so this is treated the same as a kernel structural refusal.
    /// Reachable via the WASM terminal's `--drag-model` on `--wez` (the native CLI has no such
    /// flag and always sweeps G1). `p_hit` still comes from the real Monte Carlo run, so it
    /// must be unaffected by attribution being unavailable.
    #[test]
    fn attribution_unavailable_when_drag_model_has_no_kernel_representation() {
        let result = compute_wez(
            823.0,
            0.0,
            0.243,
            0.0113,
            0.00782,
            20,
            5.0,
            0.0001,
            0.005,
            1.0,
            0.05,
            3.0,
            90.0,
            0.0,
            1.5,
            TargetSizeMetric::Rect { width_m: 0.5, height_m: 0.75 },
            300.0,
            300.0,
            100.0,
            DragModel::G2,
            None,
            1.0,
            0.0,
            0.0,
        )
        .expect("compute_wez");
        let row = result.rows.first().expect("one row");
        assert!(
            row.attribution_unavailable,
            "G2 has no solve-json v1 drag-model counterpart; attribution cannot run"
        );
        assert!(row.p_hit.is_finite(), "p_hit comes from the real Monte Carlo run, unaffected");
    }

    /// A loaded custom drag table replaces G-model+BC drag entirely, and solve-json v1 has no
    /// field for a custom deck at all -- the kernel cannot represent this configuration, so
    /// attribution is unavailable rather than silently computed under the wrong (G-model+BC)
    /// physics. Reachable via `--drag-table`/`--cd-scale` on `monte-carlo --wez`.
    #[test]
    fn attribution_unavailable_with_a_custom_drag_table() {
        let table = crate::drag::DragTable::new(
            vec![0.0, 0.5, 1.0, 1.5, 2.0, 3.0],
            vec![0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        );
        let result = compute_wez(
            823.0,
            0.0,
            0.243,
            0.0113,
            0.00782,
            20,
            5.0,
            0.0001,
            0.005,
            1.0,
            0.05,
            3.0,
            90.0,
            0.0,
            1.5,
            TargetSizeMetric::Rect { width_m: 0.5, height_m: 0.75 },
            300.0,
            300.0,
            100.0,
            DragModel::G7,
            Some(table),
            1.0,
            0.0,
            0.0,
        )
        .expect("compute_wez");
        let row = result.rows.first().expect("one row");
        assert!(
            row.attribution_unavailable,
            "a custom drag table has no solve-json v1 representation; attribution cannot run"
        );
        assert!(row.p_hit.is_finite(), "p_hit comes from the real Monte Carlo run, unaffected");
    }

    /// Characterization: pins the attribution shares for one fixed configuration.
    /// D2 changes these values (one-sided -> central differences). When this test
    /// fails during the upgrade, record the before/after in the commit message and
    /// update the expected constants -- do not silently re-baseline.
    ///
    /// Deliberately written against the public compute_wez so the refactor of the
    /// private attributor cannot change what this test measures.
    ///
    /// Range deviation from the task-13 brief's literal worked example: the brief's own
    /// numbers (angle 0.0, a single row at 600 m) are not reachable -- a perfectly level
    /// (angle = 0) shot from a 1.5 m bore height with the default ground_threshold hits the
    /// ground at ~365 m (verified via `trajectory -v 823 -a 0 -b 0.243 -m 11.3 -d 7.82
    /// --units metric --bore-height 1500 --max-range 700 -o json` => `"max_range":
    /// 364.756..."`), so `attribution_unavailable` is `true` at 600 m for TODAY's code too,
    /// before any D2 change. 300 m keeps every other parameter from the brief unchanged and
    /// sits comfortably inside that ~365 m ceiling.
    #[test]
    fn attribution_shares_for_a_fixed_configuration() {
        let result = compute_wez(
            823.0,      // velocity
            0.0,        // angle
            0.243,      // bc
            0.0113,     // mass
            0.00782,    // diameter
            20,         // num_sims
            5.0,        // velocity_std
            0.0001,     // angle_std
            0.005,      // bc_std
            1.0,        // wind_std
            0.05,       // wind_direction_std
            3.0,        // wind_speed
            90.0,       // wind_direction
            0.0,        // wind_vertical
            1.5,        // wind_call_error
            TargetSizeMetric::Rect { width_m: 0.5, height_m: 0.75 },
            300.0,      // wez_start (brief specifies 600.0; unreachable at angle=0, see above)
            300.0,      // wez_end
            100.0,      // wez_step
            DragModel::G7,
            None,       // custom_drag_table
            1.0,        // cd_scale
            0.0,        // cant
            0.0,        // sight_offset_lateral_m
        ).expect("compute_wez");

        let row = result.rows.first().expect("one row");
        assert!(!row.attribution_unavailable, "attribution must be available here");
        eprintln!("CHARACTERIZATION wind_call={} mv_sd={} other={}",
                  row.wind_call_share, row.mv_sd_share, row.other_share);
        assert!((row.wind_call_share + row.mv_sd_share + row.other_share - 1.0).abs() < 1e-9);

        // AFTER (central differences, D2). BEFORE (one-sided, pre-D2) was:
        //   wind_call=0.6820206535206624 mv_sd=0.012457193963072112 other=0.30552215251626547
        // Measured gap between the two methods on this fixture: wind_call ~4.91e-4,
        // mv_sd ~5.12e-4, other ~2.06e-5 (the tightest of the three). The 1e-6 tolerance below
        // is ~20x tighter than that smallest gap, so reverting the central-difference calls in
        // wez_variance_shares back to one-sided solves fails this test on EVERY bucket, not
        // just the two with a larger gap.
        assert!((row.wind_call_share - 0.681_529_624_943_378).abs() < 1e-6,
                 "update the characterization constant");
        assert!((row.mv_sd_share - 0.012_968_843_319_743_649).abs() < 1e-6,
                 "update the characterization constant");
        assert!((row.other_share - 0.305_501_531_736_878_2).abs() < 1e-6,
                 "update the characterization constant");
    }
}
