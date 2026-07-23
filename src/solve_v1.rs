//! Transport-free implementation of the solve-json v1 trajectory service.
//!
//! This module owns the semantic boundary between the stable, explicit-SI v1 data-transfer
//! objects and the engine's legacy internal input types.  It deliberately performs no JSON I/O,
//! filesystem access, network access, profile lookup, or terminal output.

use crate::solve_json::{
    AtmosphereV1, DragModelV1, EffectsV1, ProjectileV1, ResolvedAtmosphereV1,
    ResolvedConstantWindV1, ResolvedEffectsV1, ResolvedProjectileV1, ResolvedRifleV1,
    ResolvedSamplingV1, ResolvedSegmentedWindV1, ResolvedShotV1, ResolvedSolveRequestV1,
    ResolvedSolverV1, ResolvedWindSegmentV1, ResolvedWindV1, RifleV1, SampleFlagV1, SamplingV1,
    SchemaVersionV1, ShotV1, SolveErrorCodeV1, SolveErrorEnvelopeV1, SolveErrorV1, SolveNoticeV1,
    SolveRequestV1, SolveSuccessV1, SolveSummaryV1, SolverMethodV1, SolverV1, SuccessStatusV1,
    TerminationReasonV1, TrajectorySampleV1, TwistDirectionV1, WindV1, MAX_SOLVE_JSON_SAMPLES_V1,
};
use crate::trajectory_observation::{
    TrajectoryObservation, TrajectoryObservationError, TrajectoryObservationFlag,
    TrajectoryTermination,
};
use crate::wind::WindSegment;
use crate::{
    AtmosphericConditions, BallisticInputs, BallisticsError, DragModel, TrajectorySolver,
    WindConditions,
};

const DEFAULT_SIGHT_HEIGHT_M: f64 = 0.05;
const DEFAULT_MUZZLE_HEIGHT_M: f64 = 0.0;
const DEFAULT_TWIST_RATE_M_PER_TURN: f64 = 0.3048;
const DEFAULT_ALTITUDE_M: f64 = 0.0;
const DEFAULT_RELATIVE_HUMIDITY: f64 = 0.5;
const DEFAULT_TIME_STEP_S: f64 = 0.001;
const DEFAULT_SAMPLE_INTERVAL_M: f64 = 10.0;
const DEFAULT_GROUND_THRESHOLD_M: f64 = -100.0;

const MIN_ICAO_ALTITUDE_M: f64 = -5_000.0;
const MAX_ICAO_ALTITUDE_M: f64 = 84_000.0;
const METERS_PER_INCH: f64 = 0.0254;
const KILOMETRES_PER_METRE: f64 = 0.001;
const SECONDS_PER_HOUR: f64 = 3_600.0;
const PASCALS_PER_HECTOPASCAL: f64 = 100.0;
const KELVIN_OFFSET_C: f64 = 273.15;
const KG_PER_GRAIN: f64 = 0.000_064_798_91;

/// Stable assumption code emitted for a literal protocol default.
pub const ASSUMPTION_DEFAULT_APPLIED: &str = "default_applied";
/// Stable assumption code emitted when station temperature is resolved from ICAO atmosphere.
pub const ASSUMPTION_ICAO_STANDARD_TEMPERATURE: &str = "icao_standard_temperature";
/// Stable assumption code emitted when station pressure is resolved from ICAO atmosphere.
pub const ASSUMPTION_ICAO_STANDARD_PRESSURE: &str = "icao_standard_pressure";
/// Stable assumption code emitted when projectile length is inferred from mass and diameter.
pub const ASSUMPTION_ESTIMATED_PROJECTILE_LENGTH: &str = "estimated_projectile_length";

/// Stable warning code for segmented wind that becomes calm before the requested maximum range.
pub const WARNING_PARTIAL_WIND_COVERAGE: &str = "partial_wind_coverage";
/// Stable warning code emitted for each opt-in experimental physical effect.
pub const WARNING_EXPERIMENTAL_EFFECT: &str = "experimental_effect";
/// Stable warning code for an explicitly supplied fixed step that adaptive RK45 ignores.
pub const WARNING_RK45_TIME_STEP_IGNORED: &str = "rk45_time_step_ignored";

#[derive(Debug)]
struct PreparedSolveV1 {
    resolved_request: ResolvedSolveRequestV1,
    assumptions: Vec<SolveNoticeV1>,
    warnings: Vec<SolveNoticeV1>,
    inputs: BallisticInputs,
    wind: WindConditions,
    wind_segments: Vec<WindSegment>,
    atmosphere: AtmosphericConditions,
}

/// Execute one semantically validated solve-json v1 request.
///
/// The function is transport-free: callers are responsible for decoding and encoding the DTOs.
/// All dimensional request values are SI.  Valid requests produce a fully resolved success DTO;
/// invalid requests and engine failures produce the stable v1 error envelope.
pub fn solve_v1(request: SolveRequestV1) -> Result<SolveSuccessV1, SolveErrorEnvelopeV1> {
    let mut prepared = prepare_request(&request)?;
    let max_range_m = prepared.resolved_request.shot.max_range_m;
    let time_step_s = prepared.resolved_request.solver.time_step_s;
    let sample_interval_m = prepared.resolved_request.sampling.interval_m;
    let zero_distance_m = prepared.resolved_request.shot.zero_distance_m;
    let target_height_m = prepared.resolved_request.shot.target_height_m;
    let enhanced_spin_drift = prepared.resolved_request.effects.enhanced_spin_drift;

    // Retain the exact mapped inputs for summary diagnostics. TrajectorySolver owns its copy and
    // may update the private effective muzzle angle during zeroing.
    let summary_inputs = prepared.inputs.clone();
    let station_temperature_c = prepared.atmosphere.temperature;
    let station_pressure_hpa = prepared.atmosphere.pressure;

    let mut solver = TrajectorySolver::new_with_resolved_station_atmosphere(
        prepared.inputs,
        prepared.wind,
        prepared.atmosphere,
    );
    solver.set_max_range(max_range_m);
    solver.set_time_step(time_step_s);
    if !prepared.wind_segments.is_empty() {
        solver.set_wind_segments(prepared.wind_segments);
    }

    if let Some(distance_m) = zero_distance_m {
        let effective_angle = solver
            .calculate_and_set_zero_angle(distance_m, target_height_m)
            .map_err(solve_failed)?;
        if !effective_angle.is_finite() {
            return Err(solve_failed_message(
                "zero-angle calculation returned a non-finite effective muzzle angle",
            ));
        }
        prepared.resolved_request.shot.muzzle_angle_rad = effective_angle;
    }

    let result = solver.solve().map_err(solve_failed)?;
    let observations = result
        .sample_observations(sample_interval_m, MAX_SOLVE_JSON_SAMPLES_V1)
        .map_err(map_observation_error)?;

    let samples = observations
        .iter()
        .cloned()
        .map(observation_to_wire)
        .collect::<Vec<_>>();
    let terminal = observations.last().ok_or_else(|| {
        internal_error("trajectory observation sampling returned no terminal observation")
    })?;

    let maximum_height_m =
        maximum_world_height_m(&result, prepared.resolved_request.shot.shooting_angle_rad)?;

    let stability = crate::spin_drift::effective_sg_from_inputs(
        &summary_inputs,
        station_temperature_c,
        station_pressure_hpa,
    );
    let stability_factor = if stability.is_finite() && stability >= 0.0 {
        Some(stability)
    } else {
        None
    };
    let spin_drift_m = if enhanced_spin_drift {
        stability_factor
            .map(|sg| {
                crate::spin_drift::litz_drift_meters(
                    sg,
                    terminal.time_s,
                    summary_inputs.is_twist_right,
                )
            })
            .filter(|drift| drift.is_finite())
    } else {
        None
    };

    if enhanced_spin_drift && stability_factor.is_some() && spin_drift_m.is_none() {
        return Err(internal_error(
            "enhanced spin-drift summary calculation produced a non-finite value",
        ));
    }

    let summary = SolveSummaryV1 {
        actual_range_m: terminal.distance_m,
        maximum_height_m,
        time_of_flight_s: terminal.time_s,
        terminal_speed_mps: terminal.speed_mps,
        terminal_energy_j: terminal.energy_j,
        stability_factor,
        spin_drift_m,
        termination: termination_to_wire(result.termination),
    };

    let success = SolveSuccessV1 {
        schema_version: SchemaVersionV1,
        engine_version: env!("CARGO_PKG_VERSION").to_owned(),
        status: SuccessStatusV1::Ok,
        resolved_request: prepared.resolved_request,
        assumptions: prepared.assumptions,
        warnings: prepared.warnings,
        summary,
        samples,
    };
    success.validate_for_serialization()?;
    Ok(success)
}

fn prepare_request(request: &SolveRequestV1) -> Result<PreparedSolveV1, SolveErrorEnvelopeV1> {
    let mut assumptions = Vec::new();
    let mut warnings = Vec::new();

    let (projectile, bullet_length_m) = resolve_projectile(&request.projectile, &mut assumptions)?;
    let rifle = resolve_rifle(&request.rifle, &mut assumptions)?;
    let shot = resolve_shot(&request.shot, rifle.muzzle_height_m, &mut assumptions)?;
    let atmosphere = resolve_atmosphere(&request.atmosphere, &mut assumptions)?;
    let wind_coverage_distance_m = shot.zero_distance_m.unwrap_or(0.0).max(shot.max_range_m);
    let (resolved_wind, wind, wind_segments) = resolve_wind(
        &request.wind,
        wind_coverage_distance_m,
        &mut assumptions,
        &mut warnings,
    )?;
    let solver = resolve_solver(&request.solver, &mut assumptions, &mut warnings)?;
    let effects = resolve_effects(
        &request.effects,
        atmosphere.latitude_rad,
        &mut assumptions,
        &mut warnings,
    )?;
    let sampling = resolve_sampling(&request.sampling, &mut assumptions)?;

    let resolved_request = ResolvedSolveRequestV1 {
        schema_version: SchemaVersionV1,
        projectile,
        rifle,
        shot,
        atmosphere,
        wind: resolved_wind,
        solver,
        effects,
        sampling,
    };

    let temperature_c = checked_conversion(
        resolved_request.atmosphere.temperature_k - KELVIN_OFFSET_C,
        "$.atmosphere.temperature_k",
        "temperature conversion to Celsius overflowed",
    )?;
    let pressure_hpa = checked_conversion(
        resolved_request.atmosphere.pressure_pa / PASCALS_PER_HECTOPASCAL,
        "$.atmosphere.pressure_pa",
        "pressure conversion to hectopascals overflowed",
    )?;
    let humidity_percent = checked_conversion(
        resolved_request.atmosphere.relative_humidity * 100.0,
        "$.atmosphere.relative_humidity",
        "relative-humidity conversion to percent overflowed",
    )?;
    let twist_rate_inches = checked_conversion(
        resolved_request.rifle.twist_rate_m_per_turn / METERS_PER_INCH,
        "$.rifle.twist_rate_m_per_turn",
        "twist-rate conversion to inches per turn overflowed",
    )?;
    let caliber_inches = checked_conversion(
        resolved_request.projectile.diameter_m / METERS_PER_INCH,
        "$.projectile.diameter_m",
        "diameter conversion to inches overflowed",
    )?;
    let weight_grains = checked_conversion(
        resolved_request.projectile.mass_kg / KG_PER_GRAIN,
        "$.projectile.mass_kg",
        "mass conversion to grains overflowed",
    )?;
    let latitude_degrees = resolved_request
        .atmosphere
        .latitude_rad
        .map(f64::to_degrees)
        .map(|value| {
            checked_conversion(
                value,
                "$.atmosphere.latitude_rad",
                "latitude conversion to degrees overflowed",
            )
        })
        .transpose()?;

    let (wind_speed, wind_direction) = match &resolved_request.wind {
        ResolvedWindV1::Constant(constant) => (constant.speed_mps, constant.direction_from_rad),
        ResolvedWindV1::Segmented(_) => (0.0, 0.0),
    };

    let inputs = BallisticInputs {
        bc_value: resolved_request.projectile.ballistic_coefficient,
        bc_type: drag_model_to_engine(resolved_request.projectile.drag_model),
        bullet_mass: resolved_request.projectile.mass_kg,
        muzzle_velocity: resolved_request.rifle.muzzle_velocity_mps,
        bullet_diameter: resolved_request.projectile.diameter_m,
        bullet_length: bullet_length_m,
        muzzle_angle: resolved_request.shot.muzzle_angle_rad,
        target_distance: resolved_request.shot.max_range_m,
        azimuth_angle: resolved_request.shot.aim_azimuth_rad,
        shot_azimuth: resolved_request.shot.shot_azimuth_rad,
        shooting_angle: resolved_request.shot.shooting_angle_rad,
        cant_angle: resolved_request.shot.cant_angle_rad,
        sight_height: resolved_request.rifle.sight_height_m,
        muzzle_height: resolved_request.rifle.muzzle_height_m,
        target_height: resolved_request.shot.target_height_m,
        ground_threshold: resolved_request.shot.ground_threshold_m,
        altitude: resolved_request.atmosphere.altitude_m,
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity: resolved_request.atmosphere.relative_humidity,
        latitude: latitude_degrees,
        wind_speed,
        wind_angle: wind_direction,
        twist_rate: twist_rate_inches,
        is_twist_right: resolved_request.rifle.twist_direction == TwistDirectionV1::Right,
        caliber_inches,
        weight_grains,
        manufacturer: None,
        bullet_model: None,
        bullet_id: None,
        bullet_cluster: None,
        use_rk4: resolved_request.solver.method != SolverMethodV1::Euler,
        use_adaptive_rk45: resolved_request.solver.method == SolverMethodV1::Rk45,
        enable_advanced_effects: resolved_request.effects.magnus
            || resolved_request.effects.coriolis,
        enable_magnus: resolved_request.effects.magnus,
        enable_coriolis: resolved_request.effects.coriolis,
        use_powder_sensitivity: false,
        powder_temp_sensitivity: 0.0,
        powder_temp: temperature_c,
        powder_temp_curve: None,
        powder_curve_temp_c: None,
        tipoff_yaw: 0.0,
        tipoff_decay_distance: 50.0,
        use_bc_segments: false,
        bc_segments: None,
        bc_segments_data: None,
        use_enhanced_spin_drift: resolved_request.effects.enhanced_spin_drift,
        use_form_factor: false,
        enable_wind_shear: false,
        wind_shear_model: "none".to_owned(),
        enable_trajectory_sampling: false,
        sample_interval: resolved_request.sampling.interval_m,
        enable_pitch_damping: false,
        enable_precession_nutation: false,
        enable_aerodynamic_jump: false,
        use_cluster_bc: false,
        custom_drag_table: None,
        // MBA-1356: solve-json v1 has no custom-drag-table field (see the module doc — the
        // public JSON contract deliberately does not expose custom decks), so cd_scale is
        // always inert here; keep it at the neutral default.
        cd_scale: 1.0,
        bc_type_str: None,
    };

    let atmosphere = AtmosphericConditions {
        temperature: temperature_c,
        pressure: pressure_hpa,
        humidity: humidity_percent,
        altitude: resolved_request.atmosphere.altitude_m,
    };

    Ok(PreparedSolveV1 {
        resolved_request,
        assumptions,
        warnings,
        inputs,
        wind,
        wind_segments,
        atmosphere,
    })
}

fn resolve_projectile(
    projectile: &ProjectileV1,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> Result<(ResolvedProjectileV1, f64), SolveErrorEnvelopeV1> {
    require_positive("$.projectile.mass_kg", projectile.mass_kg)?;
    require_positive("$.projectile.diameter_m", projectile.diameter_m)?;
    require_positive(
        "$.projectile.ballistic_coefficient",
        projectile.ballistic_coefficient,
    )?;
    if let Some(length_m) = projectile.length_m {
        require_positive("$.projectile.length_m", length_m)?;
    }

    let bullet_length_m = match projectile.length_m {
        Some(length_m) => length_m,
        None => {
            let estimated = crate::stability::estimate_bullet_length_m(
                projectile.diameter_m,
                projectile.mass_kg,
            );
            if !estimated.is_finite() || estimated <= 0.0 {
                return Err(invalid_value(
                    "$.projectile",
                    "projectile mass and diameter could not produce a finite positive length estimate",
                ));
            }
            assumptions.push(notice(
                ASSUMPTION_ESTIMATED_PROJECTILE_LENGTH,
                format!(
                    "Projectile length was omitted; the engine estimated {estimated:.12} m from mass and diameter."
                ),
                "$.projectile.length_m",
            ));
            estimated
        }
    };

    Ok((
        ResolvedProjectileV1 {
            mass_kg: projectile.mass_kg,
            diameter_m: projectile.diameter_m,
            length_m: projectile.length_m,
            drag_model: projectile.drag_model,
            ballistic_coefficient: projectile.ballistic_coefficient,
        },
        bullet_length_m,
    ))
}

fn resolve_rifle(
    rifle: &RifleV1,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedRifleV1, SolveErrorEnvelopeV1> {
    require_positive("$.rifle.muzzle_velocity_mps", rifle.muzzle_velocity_mps)?;
    if let Some(value) = rifle.sight_height_m {
        require_non_negative("$.rifle.sight_height_m", value)?;
    }
    if let Some(value) = rifle.muzzle_height_m {
        require_finite("$.rifle.muzzle_height_m", value)?;
    }
    if let Some(value) = rifle.twist_rate_m_per_turn {
        require_positive("$.rifle.twist_rate_m_per_turn", value)?;
    }

    let sight_height_m = literal_default(
        rifle.sight_height_m,
        DEFAULT_SIGHT_HEIGHT_M,
        "$.rifle.sight_height_m",
        "Sight height defaulted to 0.05 m.",
        assumptions,
    );
    let muzzle_height_m = literal_default(
        rifle.muzzle_height_m,
        DEFAULT_MUZZLE_HEIGHT_M,
        "$.rifle.muzzle_height_m",
        "Muzzle height defaulted to 0 m.",
        assumptions,
    );
    let twist_rate_m_per_turn = literal_default(
        rifle.twist_rate_m_per_turn,
        DEFAULT_TWIST_RATE_M_PER_TURN,
        "$.rifle.twist_rate_m_per_turn",
        "Twist rate defaulted to 0.3048 m per turn.",
        assumptions,
    );
    let twist_direction = match rifle.twist_direction {
        Some(direction) => direction,
        None => {
            assumptions.push(notice(
                ASSUMPTION_DEFAULT_APPLIED,
                "Twist direction defaulted to right-hand.",
                "$.rifle.twist_direction",
            ));
            TwistDirectionV1::Right
        }
    };

    Ok(ResolvedRifleV1 {
        muzzle_velocity_mps: rifle.muzzle_velocity_mps,
        sight_height_m,
        muzzle_height_m,
        twist_rate_m_per_turn,
        twist_direction,
    })
}

fn resolve_shot(
    shot: &ShotV1,
    muzzle_height_m: f64,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedShotV1, SolveErrorEnvelopeV1> {
    require_positive("$.shot.max_range_m", shot.max_range_m)?;
    if let Some(value) = shot.zero_distance_m {
        require_positive("$.shot.zero_distance_m", value)?;
        if !(value * 2.0).is_finite() {
            return Err(invalid_value(
                "$.shot.zero_distance_m",
                "zero distance is too large to construct a finite trial trajectory range",
            ));
        }
    }
    if let Some(value) = shot.muzzle_angle_rad {
        require_finite("$.shot.muzzle_angle_rad", value)?;
    }
    if shot.zero_distance_m.is_some() && shot.muzzle_angle_rad.is_some() {
        return Err(conflicting_fields(
            "$.shot",
            "zero_distance_m and muzzle_angle_rad cannot both be supplied",
        ));
    }

    for (path, value) in [
        ("$.shot.aim_azimuth_rad", shot.aim_azimuth_rad),
        ("$.shot.shot_azimuth_rad", shot.shot_azimuth_rad),
        ("$.shot.shooting_angle_rad", shot.shooting_angle_rad),
        ("$.shot.cant_angle_rad", shot.cant_angle_rad),
        ("$.shot.target_height_m", shot.target_height_m),
        ("$.shot.ground_threshold_m", shot.ground_threshold_m),
    ] {
        if let Some(value) = value {
            require_finite(path, value)?;
        }
    }

    let muzzle_angle_rad = match (shot.zero_distance_m, shot.muzzle_angle_rad) {
        (_, Some(angle)) => angle,
        (Some(_), None) => 0.0, // Replaced with the calculated effective angle before solving.
        (None, None) => {
            assumptions.push(notice(
                ASSUMPTION_DEFAULT_APPLIED,
                "Muzzle angle defaulted to 0 rad.",
                "$.shot.muzzle_angle_rad",
            ));
            0.0
        }
    };
    let aim_azimuth_rad = literal_default(
        shot.aim_azimuth_rad,
        0.0,
        "$.shot.aim_azimuth_rad",
        "Aim azimuth defaulted to 0 rad.",
        assumptions,
    );
    let shot_azimuth_rad = literal_default(
        shot.shot_azimuth_rad,
        0.0,
        "$.shot.shot_azimuth_rad",
        "Shot azimuth defaulted to 0 rad (north).",
        assumptions,
    );
    let shooting_angle_rad = literal_default(
        shot.shooting_angle_rad,
        0.0,
        "$.shot.shooting_angle_rad",
        "Shooting angle defaulted to 0 rad.",
        assumptions,
    );
    let cant_angle_rad = literal_default(
        shot.cant_angle_rad,
        0.0,
        "$.shot.cant_angle_rad",
        "Cant angle defaulted to 0 rad.",
        assumptions,
    );
    let target_height_m = literal_default(
        shot.target_height_m,
        0.0,
        "$.shot.target_height_m",
        "Target height defaulted to 0 m above the local ground datum.",
        assumptions,
    );
    let ground_threshold_m = literal_default(
        shot.ground_threshold_m,
        DEFAULT_GROUND_THRESHOLD_M,
        "$.shot.ground_threshold_m",
        "Ground threshold defaulted to -100 m.",
        assumptions,
    );
    if ground_threshold_m >= muzzle_height_m {
        return Err(invalid_value(
            "$.shot.ground_threshold_m",
            "ground threshold must be below the resolved muzzle height",
        ));
    }

    Ok(ResolvedShotV1 {
        max_range_m: shot.max_range_m,
        zero_distance_m: shot.zero_distance_m,
        muzzle_angle_rad,
        aim_azimuth_rad,
        shot_azimuth_rad,
        shooting_angle_rad,
        cant_angle_rad,
        target_height_m,
        ground_threshold_m,
    })
}

fn resolve_atmosphere(
    atmosphere: &AtmosphereV1,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedAtmosphereV1, SolveErrorEnvelopeV1> {
    if let Some(value) = atmosphere.altitude_m {
        require_range(
            "$.atmosphere.altitude_m",
            value,
            MIN_ICAO_ALTITUDE_M,
            MAX_ICAO_ALTITUDE_M,
        )?;
    }
    if let Some(value) = atmosphere.temperature_k {
        require_positive("$.atmosphere.temperature_k", value)?;
    }
    if let Some(value) = atmosphere.pressure_pa {
        require_positive("$.atmosphere.pressure_pa", value)?;
    }
    if let Some(value) = atmosphere.relative_humidity {
        require_range("$.atmosphere.relative_humidity", value, 0.0, 1.0)?;
    }
    if let Some(value) = atmosphere.latitude_rad {
        require_range(
            "$.atmosphere.latitude_rad",
            value,
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::FRAC_PI_2,
        )?;
    }

    let altitude_m = match atmosphere.altitude_m {
        Some(value) => value,
        None => {
            assumptions.push(notice(
                ASSUMPTION_DEFAULT_APPLIED,
                "Station altitude defaulted to 0 m.",
                "$.atmosphere.altitude_m",
            ));
            DEFAULT_ALTITUDE_M
        }
    };
    let (standard_temperature_k, standard_pressure_pa) =
        crate::atmosphere::calculate_icao_standard_atmosphere(altitude_m);
    let temperature_k = match atmosphere.temperature_k {
        Some(value) => value,
        None => {
            assumptions.push(notice(
                ASSUMPTION_ICAO_STANDARD_TEMPERATURE,
                format!(
                    "Station temperature was omitted; ICAO standard temperature at {altitude_m:.6} m is {standard_temperature_k:.12} K."
                ),
                "$.atmosphere.temperature_k",
            ));
            standard_temperature_k
        }
    };
    let pressure_pa = match atmosphere.pressure_pa {
        Some(value) => value,
        None => {
            assumptions.push(notice(
                ASSUMPTION_ICAO_STANDARD_PRESSURE,
                format!(
                    "Station pressure was omitted; ICAO standard pressure at {altitude_m:.6} m is {standard_pressure_pa:.12} Pa."
                ),
                "$.atmosphere.pressure_pa",
            ));
            standard_pressure_pa
        }
    };
    let relative_humidity = literal_default(
        atmosphere.relative_humidity,
        DEFAULT_RELATIVE_HUMIDITY,
        "$.atmosphere.relative_humidity",
        "Relative humidity defaulted to 0.5.",
        assumptions,
    );

    Ok(ResolvedAtmosphereV1 {
        altitude_m,
        temperature_k,
        pressure_pa,
        relative_humidity,
        latitude_rad: atmosphere.latitude_rad,
    })
}

fn resolve_wind(
    wind: &WindV1,
    coverage_distance_m: f64,
    assumptions: &mut Vec<SolveNoticeV1>,
    warnings: &mut Vec<SolveNoticeV1>,
) -> Result<(ResolvedWindV1, WindConditions, Vec<WindSegment>), SolveErrorEnvelopeV1> {
    if let Some(segments) = &wind.segments {
        if wind.speed_mps.is_some()
            || wind.direction_from_rad.is_some()
            || wind.vertical_speed_mps.is_some()
        {
            return Err(conflicting_fields(
                "$.wind",
                "segments cannot be combined with constant-wind fields",
            ));
        }
        if segments.is_empty() {
            return Err(invalid_value(
                "$.wind.segments",
                "segmented wind must contain at least one segment",
            ));
        }

        let mut resolved_segments = Vec::with_capacity(segments.len());
        let mut engine_segments = Vec::with_capacity(segments.len());
        let mut previous_until_m = 0.0;
        for (index, segment) in segments.iter().enumerate() {
            let base = format!("$.wind.segments[{index}]");
            require_positive(
                &format!("{base}.until_distance_m"),
                segment.until_distance_m,
            )?;
            require_non_negative(&format!("{base}.speed_mps"), segment.speed_mps)?;
            require_finite(
                &format!("{base}.direction_from_rad"),
                segment.direction_from_rad,
            )?;
            if let Some(value) = segment.vertical_speed_mps {
                require_finite(&format!("{base}.vertical_speed_mps"), value)?;
            }
            if index > 0 && segment.until_distance_m <= previous_until_m {
                return Err(invalid_value(
                    format!("{base}.until_distance_m"),
                    "wind-segment boundaries must be strictly increasing",
                ));
            }
            previous_until_m = segment.until_distance_m;

            let vertical_speed_mps = match segment.vertical_speed_mps {
                Some(value) => value,
                None => {
                    assumptions.push(notice(
                        ASSUMPTION_DEFAULT_APPLIED,
                        "Segment vertical wind speed defaulted to 0 m/s.",
                        format!("{base}.vertical_speed_mps"),
                    ));
                    0.0
                }
            };
            let speed_kmh = checked_conversion(
                segment.speed_mps * SECONDS_PER_HOUR * KILOMETRES_PER_METRE,
                &format!("{base}.speed_mps"),
                "wind speed conversion to kilometres per hour overflowed",
            )?;
            let angle_deg = checked_conversion(
                segment.direction_from_rad.to_degrees(),
                &format!("{base}.direction_from_rad"),
                "wind direction conversion to degrees overflowed",
            )?;

            resolved_segments.push(ResolvedWindSegmentV1 {
                until_distance_m: segment.until_distance_m,
                speed_mps: segment.speed_mps,
                direction_from_rad: segment.direction_from_rad,
                vertical_speed_mps,
            });
            engine_segments.push(WindSegment {
                speed_kmh,
                angle_deg,
                until_m: segment.until_distance_m,
                vertical_mps: vertical_speed_mps,
            });
        }

        let final_until_m = resolved_segments
            .last()
            .expect("non-empty segmented wind was checked")
            .until_distance_m;
        if final_until_m < coverage_distance_m {
            warnings.push(notice(
                WARNING_PARTIAL_WIND_COVERAGE,
                format!(
                    "Segmented wind ends at {final_until_m:.12} m; wind is calm from that boundary through the required {coverage_distance_m:.12} m solve/zero range."
                ),
                "$.wind.segments",
            ));
        }

        Ok((
            ResolvedWindV1::Segmented(ResolvedSegmentedWindV1 {
                segments: resolved_segments,
            }),
            WindConditions::default(),
            engine_segments,
        ))
    } else {
        match (wind.speed_mps, wind.direction_from_rad) {
            (Some(_), None) | (None, Some(_)) => {
                return Err(conflicting_fields(
                    "$.wind",
                    "speed_mps and direction_from_rad must be supplied together",
                ));
            }
            (None, None) if wind.vertical_speed_mps.is_some() => {
                return Err(conflicting_fields(
                    "$.wind",
                    "vertical_speed_mps requires speed_mps and direction_from_rad",
                ));
            }
            _ => {}
        }

        if let Some(value) = wind.speed_mps {
            require_non_negative("$.wind.speed_mps", value)?;
        }
        if let Some(value) = wind.direction_from_rad {
            require_finite("$.wind.direction_from_rad", value)?;
        }
        if let Some(value) = wind.vertical_speed_mps {
            require_finite("$.wind.vertical_speed_mps", value)?;
        }

        let (speed_mps, direction_from_rad, vertical_speed_mps) =
            if let (Some(speed), Some(direction)) = (wind.speed_mps, wind.direction_from_rad) {
                let vertical = match wind.vertical_speed_mps {
                    Some(value) => value,
                    None => {
                        assumptions.push(notice(
                            ASSUMPTION_DEFAULT_APPLIED,
                            "Constant vertical wind speed defaulted to 0 m/s.",
                            "$.wind.vertical_speed_mps",
                        ));
                        0.0
                    }
                };
                (speed, direction, vertical)
            } else {
                assumptions.push(notice(
                    ASSUMPTION_DEFAULT_APPLIED,
                    "Wind defaulted to still air.",
                    "$.wind",
                ));
                (0.0, 0.0, 0.0)
            };

        Ok((
            ResolvedWindV1::Constant(ResolvedConstantWindV1 {
                speed_mps,
                direction_from_rad,
                vertical_speed_mps,
            }),
            WindConditions {
                speed: speed_mps,
                direction: direction_from_rad,
                vertical_speed: vertical_speed_mps,
            },
            Vec::new(),
        ))
    }
}

fn resolve_solver(
    solver: &SolverV1,
    assumptions: &mut Vec<SolveNoticeV1>,
    warnings: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedSolverV1, SolveErrorEnvelopeV1> {
    if let Some(value) = solver.time_step_s {
        require_positive("$.solver.time_step_s", value)?;
    }
    let method = match solver.method {
        Some(method) => method,
        None => {
            assumptions.push(notice(
                ASSUMPTION_DEFAULT_APPLIED,
                "Solver method defaulted to adaptive RK45.",
                "$.solver.method",
            ));
            SolverMethodV1::Rk45
        }
    };
    let time_step_s = literal_default(
        solver.time_step_s,
        DEFAULT_TIME_STEP_S,
        "$.solver.time_step_s",
        "Solver time step defaulted to 0.001 s.",
        assumptions,
    );
    if method == SolverMethodV1::Rk45 && solver.time_step_s.is_some() {
        warnings.push(notice(
            WARNING_RK45_TIME_STEP_IGNORED,
            "Adaptive RK45 owns its time step; the explicitly supplied fixed time step is retained in resolved_request but ignored by integration.",
            "$.solver.time_step_s",
        ));
    }

    Ok(ResolvedSolverV1 {
        method,
        time_step_s,
    })
}

fn resolve_effects(
    effects: &EffectsV1,
    latitude_rad: Option<f64>,
    assumptions: &mut Vec<SolveNoticeV1>,
    warnings: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedEffectsV1, SolveErrorEnvelopeV1> {
    let magnus = bool_default(
        effects.magnus,
        "$.effects.magnus",
        "Magnus effect defaulted to disabled.",
        assumptions,
    );
    let coriolis = bool_default(
        effects.coriolis,
        "$.effects.coriolis",
        "Coriolis effect defaulted to disabled.",
        assumptions,
    );
    let enhanced_spin_drift = bool_default(
        effects.enhanced_spin_drift,
        "$.effects.enhanced_spin_drift",
        "Enhanced spin drift defaulted to disabled.",
        assumptions,
    );

    if magnus && enhanced_spin_drift {
        return Err(conflicting_fields(
            "$.effects",
            "magnus and enhanced_spin_drift cannot both be enabled",
        ));
    }
    if coriolis && latitude_rad.is_none() {
        return Err(invalid_value(
            "$.atmosphere.latitude_rad",
            "latitude_rad is required when the Coriolis effect is enabled",
        ));
    }

    for (enabled, path, name) in [
        (magnus, "$.effects.magnus", "Magnus"),
        (
            enhanced_spin_drift,
            "$.effects.enhanced_spin_drift",
            "enhanced spin drift",
        ),
    ] {
        if enabled {
            warnings.push(notice(
                WARNING_EXPERIMENTAL_EFFECT,
                format!("The {name} effect is an opt-in experimental v1 model."),
                path,
            ));
        }
    }

    Ok(ResolvedEffectsV1 {
        magnus,
        coriolis,
        enhanced_spin_drift,
    })
}

fn resolve_sampling(
    sampling: &SamplingV1,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> Result<ResolvedSamplingV1, SolveErrorEnvelopeV1> {
    if let Some(value) = sampling.interval_m {
        require_positive("$.sampling.interval_m", value)?;
    }
    let interval_m = literal_default(
        sampling.interval_m,
        DEFAULT_SAMPLE_INTERVAL_M,
        "$.sampling.interval_m",
        "Sampling interval defaulted to 10 m.",
        assumptions,
    );
    Ok(ResolvedSamplingV1 { interval_m })
}

fn maximum_world_height_m(
    result: &crate::TrajectoryResult,
    shooting_angle_rad: f64,
) -> Result<f64, SolveErrorEnvelopeV1> {
    let mut maximum = f64::NEG_INFINITY;
    for point in &result.points {
        let height = crate::atmosphere::shot_frame_altitude(
            0.0,
            point.position.x,
            point.position.y,
            shooting_angle_rad,
        );
        if !height.is_finite() {
            return Err(internal_error(
                "world-vertical maximum-height projection produced a non-finite value",
            ));
        }
        maximum = maximum.max(height);
    }
    if maximum.is_finite() {
        Ok(maximum)
    } else {
        Err(internal_error(
            "trajectory contained no finite points for maximum-height calculation",
        ))
    }
}

fn observation_to_wire(observation: TrajectoryObservation) -> TrajectorySampleV1 {
    TrajectorySampleV1 {
        distance_m: observation.distance_m,
        time_s: observation.time_s,
        speed_mps: observation.speed_mps,
        energy_j: observation.energy_j,
        drop_m: observation.drop_m,
        windage_m: observation.windage_m,
        mach: observation.mach,
        flags: observation
            .flags
            .into_iter()
            .map(observation_flag_to_wire)
            .collect(),
    }
}

fn observation_flag_to_wire(flag: TrajectoryObservationFlag) -> SampleFlagV1 {
    match flag {
        TrajectoryObservationFlag::Transonic => SampleFlagV1::Transonic,
        TrajectoryObservationFlag::Subsonic => SampleFlagV1::Subsonic,
        TrajectoryObservationFlag::Terminal => SampleFlagV1::Terminal,
        TrajectoryObservationFlag::GroundThreshold => SampleFlagV1::GroundThreshold,
    }
}

fn termination_to_wire(termination: TrajectoryTermination) -> TerminationReasonV1 {
    match termination {
        TrajectoryTermination::MaxRange => TerminationReasonV1::MaxRange,
        TrajectoryTermination::GroundThreshold => TerminationReasonV1::GroundThreshold,
        TrajectoryTermination::TimeLimit => TerminationReasonV1::TimeLimit,
        TrajectoryTermination::VelocityFloor => TerminationReasonV1::VelocityFloor,
    }
}

fn drag_model_to_engine(model: DragModelV1) -> DragModel {
    match model {
        DragModelV1::G1 => DragModel::G1,
        DragModelV1::G6 => DragModel::G6,
        DragModelV1::G7 => DragModel::G7,
        DragModelV1::G8 => DragModel::G8,
    }
}

fn map_observation_error(error: TrajectoryObservationError) -> SolveErrorEnvelopeV1 {
    let message = error.to_string();
    match error {
        TrajectoryObservationError::SampleLimitExceeded { .. }
        | TrajectoryObservationError::AllocationFailed { .. } => error_at(
            SolveErrorCodeV1::ResourceLimit,
            message,
            "$.sampling.interval_m",
        ),
        TrajectoryObservationError::InvalidInterval { .. }
        | TrajectoryObservationError::UnrepresentableGrid { .. } => error_at(
            SolveErrorCodeV1::InvalidValue,
            message,
            "$.sampling.interval_m",
        ),
        TrajectoryObservationError::EmptyTrajectory
        | TrajectoryObservationError::NonFiniteQuery { .. }
        | TrajectoryObservationError::OutOfRange { .. }
        | TrajectoryObservationError::NonMonotonicTrajectory { .. }
        | TrajectoryObservationError::NonFiniteState { .. }
        | TrajectoryObservationError::InvalidState { .. }
        | TrajectoryObservationError::InvalidMetadata { .. }
        | TrajectoryObservationError::NonFiniteObservation { .. } => internal_error(message),
    }
}

fn solve_failed(error: BallisticsError) -> SolveErrorEnvelopeV1 {
    solve_failed_message(error.to_string())
}

fn solve_failed_message(message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(SolveErrorV1::new(SolveErrorCodeV1::SolveFailed, message))
}

fn internal_error(message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(SolveErrorV1::new(SolveErrorCodeV1::InternalError, message))
}

fn invalid_value(path: impl Into<String>, message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    error_at(SolveErrorCodeV1::InvalidValue, message, path)
}

fn conflicting_fields(path: impl Into<String>, message: impl Into<String>) -> SolveErrorEnvelopeV1 {
    error_at(SolveErrorCodeV1::ConflictingFields, message, path)
}

fn error_at(
    code: SolveErrorCodeV1,
    message: impl Into<String>,
    path: impl Into<String>,
) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(SolveErrorV1::new(code, message).at_path(path))
}

fn notice(
    code: impl Into<String>,
    message: impl Into<String>,
    path: impl Into<String>,
) -> SolveNoticeV1 {
    SolveNoticeV1 {
        code: code.into(),
        message: message.into(),
        path: Some(path.into()),
    }
}

fn literal_default(
    value: Option<f64>,
    default: f64,
    path: &'static str,
    message: &'static str,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> f64 {
    match value {
        Some(value) => value,
        None => {
            assumptions.push(notice(ASSUMPTION_DEFAULT_APPLIED, message, path));
            default
        }
    }
}

fn bool_default(
    value: Option<bool>,
    path: &'static str,
    message: &'static str,
    assumptions: &mut Vec<SolveNoticeV1>,
) -> bool {
    match value {
        Some(value) => value,
        None => {
            assumptions.push(notice(ASSUMPTION_DEFAULT_APPLIED, message, path));
            false
        }
    }
}

fn checked_conversion(value: f64, path: &str, message: &str) -> Result<f64, SolveErrorEnvelopeV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid_value(path, message))
    }
}

fn require_finite(path: &str, value: f64) -> Result<(), SolveErrorEnvelopeV1> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(invalid_value(path, "value must be finite"))
    }
}

fn require_positive(path: &str, value: f64) -> Result<(), SolveErrorEnvelopeV1> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(invalid_value(
            path,
            "value must be finite and greater than zero",
        ))
    }
}

fn require_non_negative(path: &str, value: f64) -> Result<(), SolveErrorEnvelopeV1> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(invalid_value(path, "value must be finite and non-negative"))
    }
}

fn require_range(
    path: &str,
    value: f64,
    minimum: f64,
    maximum: f64,
) -> Result<(), SolveErrorEnvelopeV1> {
    if value.is_finite() && (minimum..=maximum).contains(&value) {
        Ok(())
    } else {
        Err(invalid_value(
            path,
            format!("value must be finite and in the inclusive range [{minimum}, {maximum}]"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_request() -> SolveRequestV1 {
        SolveRequestV1 {
            schema_version: SchemaVersionV1,
            projectile: ProjectileV1 {
                mass_kg: 0.01134,
                diameter_m: 0.00782,
                length_m: None,
                drag_model: DragModelV1::G7,
                ballistic_coefficient: 0.243,
            },
            rifle: RifleV1 {
                muzzle_velocity_mps: 823.0,
                sight_height_m: None,
                muzzle_height_m: None,
                twist_rate_m_per_turn: None,
                twist_direction: None,
            },
            shot: ShotV1 {
                max_range_m: 1_000.0,
                zero_distance_m: None,
                muzzle_angle_rad: None,
                aim_azimuth_rad: None,
                shot_azimuth_rad: None,
                shooting_angle_rad: None,
                cant_angle_rad: None,
                target_height_m: None,
                ground_threshold_m: None,
            },
            atmosphere: AtmosphereV1::default(),
            wind: WindV1::default(),
            solver: SolverV1::default(),
            effects: EffectsV1::default(),
            sampling: SamplingV1::default(),
        }
    }

    #[test]
    fn defaults_are_resolved_in_deterministic_request_field_order() {
        let request = minimal_request();
        let prepared = prepare_request(&request).expect("valid request");
        let code_paths = prepared
            .assumptions
            .iter()
            .map(|notice| (notice.code.as_str(), notice.path.as_deref().unwrap()))
            .collect::<Vec<_>>();

        assert_eq!(
            code_paths,
            vec![
                (
                    ASSUMPTION_ESTIMATED_PROJECTILE_LENGTH,
                    "$.projectile.length_m"
                ),
                (ASSUMPTION_DEFAULT_APPLIED, "$.rifle.sight_height_m"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.rifle.muzzle_height_m"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.rifle.twist_rate_m_per_turn"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.rifle.twist_direction"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.muzzle_angle_rad"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.aim_azimuth_rad"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.shot_azimuth_rad"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.shooting_angle_rad"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.cant_angle_rad"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.target_height_m"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.shot.ground_threshold_m"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.atmosphere.altitude_m"),
                (
                    ASSUMPTION_ICAO_STANDARD_TEMPERATURE,
                    "$.atmosphere.temperature_k"
                ),
                (
                    ASSUMPTION_ICAO_STANDARD_PRESSURE,
                    "$.atmosphere.pressure_pa"
                ),
                (ASSUMPTION_DEFAULT_APPLIED, "$.atmosphere.relative_humidity"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.wind"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.solver.method"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.solver.time_step_s"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.effects.magnus"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.effects.coriolis"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.effects.enhanced_spin_drift"),
                (ASSUMPTION_DEFAULT_APPLIED, "$.sampling.interval_m"),
            ]
        );
        assert_eq!(
            prepared.resolved_request.rifle.twist_rate_m_per_turn,
            0.3048
        );
        assert!((prepared.inputs.twist_rate - 12.0).abs() < 1.0e-12);
        assert_eq!(prepared.atmosphere.humidity, 50.0);
        assert_eq!(prepared.inputs.humidity, 0.5);
    }

    #[test]
    fn explicit_standard_station_values_at_altitude_remain_authoritative() {
        let mut request = minimal_request();
        request.atmosphere.altitude_m = Some(1_500.0);
        request.atmosphere.temperature_k = Some(288.15);
        request.atmosphere.pressure_pa = Some(101_325.0);

        let prepared = prepare_request(&request).expect("valid explicit station conditions");
        assert_eq!(prepared.atmosphere.temperature, 15.0);
        assert_eq!(prepared.atmosphere.pressure, 1013.25);
        assert_eq!(prepared.resolved_request.atmosphere.temperature_k, 288.15);
        assert_eq!(prepared.resolved_request.atmosphere.pressure_pa, 101_325.0);
        assert!(!prepared
            .assumptions
            .iter()
            .any(|notice| notice.code == ASSUMPTION_ICAO_STANDARD_TEMPERATURE
                || notice.code == ASSUMPTION_ICAO_STANDARD_PRESSURE));
    }

    #[test]
    fn partial_segmented_wind_maps_units_and_warns_without_extending_coverage() {
        let mut request = minimal_request();
        request.wind.segments = Some(vec![crate::solve_json::WindSegmentV1 {
            until_distance_m: 300.0,
            speed_mps: 5.0,
            direction_from_rad: std::f64::consts::FRAC_PI_2,
            vertical_speed_mps: None,
        }]);

        let prepared = prepare_request(&request).expect("partial wind is valid");
        assert_eq!(prepared.wind_segments.len(), 1);
        assert_eq!(prepared.wind_segments[0].speed_kmh, 18.0);
        assert_eq!(prepared.wind_segments[0].angle_deg, 90.0);
        assert_eq!(prepared.wind_segments[0].until_m, 300.0);
        assert_eq!(prepared.wind_segments[0].vertical_mps, 0.0);
        assert!(prepared.warnings.iter().any(|notice| {
            notice.code == WARNING_PARTIAL_WIND_COVERAGE
                && notice.path.as_deref() == Some("$.wind.segments")
        }));
        assert!(matches!(
            prepared.resolved_request.wind,
            ResolvedWindV1::Segmented(_)
        ));
    }

    #[test]
    fn partial_wind_coverage_includes_a_zero_beyond_the_output_range() {
        let mut request = minimal_request();
        request.shot.max_range_m = 200.0;
        request.shot.zero_distance_m = Some(500.0);
        request.wind.segments = Some(vec![crate::solve_json::WindSegmentV1 {
            until_distance_m: 300.0,
            speed_mps: 5.0,
            direction_from_rad: 0.0,
            vertical_speed_mps: Some(0.0),
        }]);

        let prepared = prepare_request(&request).expect("partial zero wind is valid");
        assert!(prepared.warnings.iter().any(|notice| {
            notice.code == WARNING_PARTIAL_WIND_COVERAGE
                && notice.path.as_deref() == Some("$.wind.segments")
                && notice.message.contains("500.000000000000 m")
        }));
    }

    #[test]
    fn semantic_conflicts_and_invalid_values_keep_exact_paths() {
        let mut request = minimal_request();
        request.shot.zero_distance_m = Some(100.0);
        request.shot.muzzle_angle_rad = Some(0.01);
        let error = prepare_request(&request).expect_err("zero and angle conflict");
        assert_eq!(error.error.code, SolveErrorCodeV1::ConflictingFields);
        assert_eq!(error.error.path(), Some("$.shot"));

        let mut request = minimal_request();
        request.atmosphere.relative_humidity = Some(1.01);
        let error = prepare_request(&request).expect_err("humidity above one");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.error.path(), Some("$.atmosphere.relative_humidity"));

        let mut request = minimal_request();
        request.effects.coriolis = Some(true);
        let error = prepare_request(&request).expect_err("Coriolis needs latitude");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.error.path(), Some("$.atmosphere.latitude_rad"));

        let mut request = minimal_request();
        request.shot.zero_distance_m = Some(f64::MAX);
        let error = prepare_request(&request).expect_err("zero trial range must remain finite");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.error.path(), Some("$.shot.zero_distance_m"));
    }

    #[test]
    fn enum_mappers_cover_every_engine_metadata_variant() {
        assert_eq!(drag_model_to_engine(DragModelV1::G1), DragModel::G1);
        assert_eq!(drag_model_to_engine(DragModelV1::G6), DragModel::G6);
        assert_eq!(drag_model_to_engine(DragModelV1::G7), DragModel::G7);
        assert_eq!(drag_model_to_engine(DragModelV1::G8), DragModel::G8);

        for (engine, wire) in [
            (
                TrajectoryTermination::MaxRange,
                TerminationReasonV1::MaxRange,
            ),
            (
                TrajectoryTermination::GroundThreshold,
                TerminationReasonV1::GroundThreshold,
            ),
            (
                TrajectoryTermination::TimeLimit,
                TerminationReasonV1::TimeLimit,
            ),
            (
                TrajectoryTermination::VelocityFloor,
                TerminationReasonV1::VelocityFloor,
            ),
        ] {
            assert_eq!(termination_to_wire(engine), wire);
        }
    }
}
