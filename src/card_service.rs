//! Transport-free card-generation services (come-ups, range table, wind card).
//!
//! These replicate the CLI's card commands exactly — same zero solve, same sampled
//! trajectory, same nearest-sample row selection, same adjustment/bias/CF/click
//! ordering via [`crate::adjustment`] — as versioned request/response types the
//! bridge can expose to embedded consumers. The CLI remains the reference
//! implementation; `tests/card_bridge_golden.rs` asserts row-for-row agreement.
//!
//! ## Unit convention
//!
//! Unlike solve-json (explicit SI), card requests are denominated in the declared
//! `units` system, exactly like the CLI flags they mirror: imperial = fps, grains,
//! inches, yards, °F, inHg, mph; metric = m/s, grams, mm, meters, °C, hPa, m/s.
//! A DOPE card is a display artifact; its inputs and outputs share the shooter's
//! unit world, and this keeps the request shape identical to the documented CLI
//! surface.

use serde::{Deserialize, Serialize};

use crate::adjustment::{
    adjustment_display, adjustment_unit_label, parse_click_value, windage_adjustment_display,
    AdjustmentUnit, ClickValue,
};
use crate::hold_curve::run_sampled_trajectory;
use crate::{AtmosphericConditions, BallisticInputs, DragModel, WindConditions};

/// Schema version of the card request/response contract.
pub const CARD_SCHEMA_VERSION_V1: u32 = 1;

const GRAINS_TO_KG: f64 = crate::constants::GRAINS_TO_KG;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CardUnits {
    #[default]
    Imperial,
    Metric,
}

/// Everything the three card surfaces share: load, zero, atmosphere, display axes.
/// Field values are in the `units` system (see module docs).
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CardRequestV1 {
    #[serde(default)]
    pub units: CardUnits,

    // Load
    pub muzzle_velocity: f64,
    pub ballistic_coefficient: f64,
    /// grains (imperial) / grams (metric)
    pub mass: f64,
    /// inches (imperial) / mm (metric)
    pub diameter: f64,
    #[serde(default)]
    pub drag_model: DragModelV1,
    /// inches (imperial) / mm (metric); CLI default 1.5 in
    #[serde(default = "default_sight_height")]
    pub sight_height: f64,

    // Zero
    pub zero_distance: f64,
    /// Deliberate POI offset at the zero range, in the units system's LINEAR drop
    /// unit (inches / mm); 0 = zeroed dead-on.
    #[serde(default)]
    pub zero_poi_vertical: f64,
    #[serde(default)]
    pub zero_poi_horizontal: f64,
    /// Lateral sight-to-bore mount offset (inches / mm).
    #[serde(default)]
    pub sight_offset_lateral: f64,

    // Atmosphere (CLI defaults: 59 °F / 15 °C, 29.92 inHg / 1013.25 hPa, 50 %, 0 alt)
    #[serde(default)]
    pub temperature: Option<f64>,
    #[serde(default)]
    pub pressure: Option<f64>,
    #[serde(default = "default_humidity")]
    pub humidity: f64,
    #[serde(default)]
    pub altitude: f64,

    // Wind for the flight (mph / m/s; wind-FROM degrees). The wind card ignores
    // these and sweeps its own speed list.
    #[serde(default)]
    pub wind_speed: f64,
    #[serde(default)]
    pub wind_direction_deg: f64,

    // Card domain (yards / meters)
    pub start: f64,
    pub end: f64,
    pub step: f64,

    // Display axes
    #[serde(default)]
    pub adjustment_unit: AdjustmentUnit,
    /// Windage axis unit; defaults to the elevation axis.
    #[serde(default)]
    pub windage_unit: Option<AdjustmentUnit>,
    /// Turret graduations as suffixed strings ("0.1mil", "0.25moa"); required when
    /// the corresponding axis unit is `clicks`. Windage falls back to elevation.
    #[serde(default)]
    pub elevation_click_value: Option<String>,
    #[serde(default)]
    pub windage_click_value: Option<String>,
    /// Scope tracking correction factors (MBA-1358); 1.0 = tracks true.
    #[serde(default = "default_cf")]
    pub elevation_cf: f64,
    #[serde(default = "default_cf")]
    pub windage_cf: f64,
    /// Selected zero set's dial corrections, true angular MIL (MBA-1360).
    #[serde(default)]
    pub zero_set_elevation_bias_mil: f64,
    #[serde(default)]
    pub zero_set_windage_bias_mil: f64,

    /// Wind card only: the sweep of wind speeds (mph / m/s), one output column each.
    #[serde(default)]
    pub wind_speeds: Vec<f64>,
    /// Wind card only: wind-FROM angles in degrees; default `[90]` (full-value from
    /// the right), matching the CLI.
    #[serde(default)]
    pub wind_angles_deg: Vec<f64>,
}

fn default_sight_height() -> f64 {
    // The service can't know units at field-default time; resolved in `resolve()`.
    f64::NAN
}
fn default_humidity() -> f64 {
    50.0
}
fn default_cf() -> f64 {
    1.0
}

/// Drag model selector mirroring the CLI's `--drag-model` strings.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DragModelV1 {
    G1,
    #[default]
    G7,
    G2,
    G5,
    G6,
    G8,
    GI,
    GS,
    RA4,
}

impl From<DragModelV1> for DragModel {
    fn from(v: DragModelV1) -> Self {
        match v {
            DragModelV1::G1 => DragModel::G1,
            DragModelV1::G7 => DragModel::G7,
            DragModelV1::G2 => DragModel::G2,
            DragModelV1::G5 => DragModel::G5,
            DragModelV1::G6 => DragModel::G6,
            DragModelV1::G8 => DragModel::G8,
            DragModelV1::GI => DragModel::GI,
            DragModelV1::GS => DragModel::GS,
            DragModelV1::RA4 => DragModel::RA4,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct CardRowV1 {
    /// Range in the request's distance unit (yd / m).
    pub range: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub drop_linear: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub drop_adj: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub come_up: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_linear: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_adj: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub velocity: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub energy: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time: Option<f64>,
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub wind_columns: Vec<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CardUnitsBlockV1 {
    pub distance: &'static str,
    pub velocity: &'static str,
    pub energy: &'static str,
    pub drop: &'static str,
    pub wind_speed: &'static str,
    pub elevation_adjustment: String,
    pub windage_adjustment: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct CardResponseV1 {
    pub schema_version: u32,
    pub kind: &'static str,
    pub zero_distance: f64,
    pub units: CardUnitsBlockV1,
    /// Wind card only: the swept speeds, one per `wind_columns` entry.
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub wind_speeds: Vec<f64>,
    /// Wind card only: one block of rows per wind angle.
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub wind_angles_deg: Vec<f64>,
    pub rows: Vec<CardRowV1>,
    /// Wind card with multiple angles: rows for angles beyond the first.
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub extra_angle_rows: Vec<Vec<CardRowV1>>,
}

#[derive(Debug, thiserror::Error)]
pub enum CardServiceError {
    #[error("invalid request: {0}")]
    InvalidRequest(String),
    #[error("zero solve failed: {0}")]
    ZeroSolve(String),
    #[error("trajectory failed: {0}")]
    Trajectory(String),
}

// --- Unit conversions: byte-identical constants to the CLI's UnitConverter ---

struct Units {
    imperial: bool,
}

impl Units {
    fn velocity_to_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 0.3048 } else { v }
    }
    fn mass_to_kg(&self, v: f64) -> f64 {
        if self.imperial { v * GRAINS_TO_KG } else { v * 0.001 }
    }
    fn distance_to_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 0.9144 } else { v }
    }
    fn small_len_to_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 0.0254 } else { v * 0.001 }
    }
    fn wind_to_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 0.44704 } else { v }
    }
    fn temperature_to_metric(&self, v: f64) -> f64 {
        if self.imperial { (v - 32.0) * 5.0 / 9.0 } else { v }
    }
    fn pressure_to_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 33.8639 } else { v }
    }
    fn velocity_from_metric(&self, v: f64) -> f64 {
        if self.imperial { v / 0.3048 } else { v }
    }
    fn distance_from_metric(&self, v: f64) -> f64 {
        if self.imperial { v / 0.9144 } else { v }
    }
    fn energy_from_metric(&self, v: f64) -> f64 {
        if self.imperial { v * 0.737562 } else { v }
    }
    fn drop_linear_from_metric(&self, v: f64) -> f64 {
        if self.imperial { v / 0.0254 } else { v * 1000.0 }
    }
}

struct Resolved {
    u: Units,
    velocity_m: f64,
    mass_kg: f64,
    diameter_m: f64,
    sight_height_m: f64,
    zero_distance_m: f64,
    zero_poi_vertical_m: f64,
    zero_poi_horizontal_m: f64,
    sight_offset_lateral_m: f64,
    temperature_c: f64,
    pressure_hpa: f64,
    end_m: f64,
    sample_m: f64,
    elevation_click: Option<ClickValue>,
    windage_click: Option<ClickValue>,
    windage_unit: AdjustmentUnit,
    drag_model: DragModel,
}

fn resolve(req: &CardRequestV1) -> Result<Resolved, CardServiceError> {
    let imperial = req.units == CardUnits::Imperial;
    let u = Units { imperial };

    for (name, v) in [
        ("muzzle_velocity", req.muzzle_velocity),
        ("ballistic_coefficient", req.ballistic_coefficient),
        ("mass", req.mass),
        ("diameter", req.diameter),
        ("zero_distance", req.zero_distance),
        ("start", req.start),
        ("end", req.end),
        ("step", req.step),
    ] {
        if !v.is_finite() || v <= 0.0 {
            return Err(CardServiceError::InvalidRequest(format!(
                "{name} must be finite and positive, got {v}"
            )));
        }
    }
    if req.end < req.start {
        return Err(CardServiceError::InvalidRequest(
            "end must be >= start".into(),
        ));
    }
    if !(0.5..=1.5).contains(&req.elevation_cf) || !(0.5..=1.5).contains(&req.windage_cf) {
        return Err(CardServiceError::InvalidRequest(
            "tracking correction factors must lie in (0.5, 1.5)".into(),
        ));
    }

    let sight_height = if req.sight_height.is_nan() {
        if imperial {
            1.5
        } else {
            38.0
        }
    } else {
        req.sight_height
    };
    let temperature = req.temperature.unwrap_or(if imperial { 59.0 } else { 15.0 });
    let pressure = req.pressure.unwrap_or(if imperial { 29.92 } else { 1013.25 });

    let parse_click = |label: &str, s: &Option<String>| -> Result<Option<ClickValue>, CardServiceError> {
        s.as_deref()
            .map(|raw| {
                parse_click_value(raw).map_err(|e| {
                    CardServiceError::InvalidRequest(format!("{label}: {e}"))
                })
            })
            .transpose()
    };
    let elevation_click = parse_click("elevation_click_value", &req.elevation_click_value)?;
    let windage_click =
        parse_click("windage_click_value", &req.windage_click_value)?.or(elevation_click);

    if req.adjustment_unit == AdjustmentUnit::Clicks && elevation_click.is_none() {
        return Err(CardServiceError::InvalidRequest(
            "adjustment_unit 'clicks' requires elevation_click_value".into(),
        ));
    }
    let windage_unit = req.windage_unit.unwrap_or(req.adjustment_unit);
    if windage_unit == AdjustmentUnit::Clicks && windage_click.is_none() {
        return Err(CardServiceError::InvalidRequest(
            "windage_unit 'clicks' requires a click graduation".into(),
        ));
    }

    Ok(Resolved {
        velocity_m: u.velocity_to_metric(req.muzzle_velocity),
        mass_kg: u.mass_to_kg(req.mass),
        diameter_m: u.small_len_to_metric(req.diameter),
        sight_height_m: u.small_len_to_metric(sight_height),
        zero_distance_m: u.distance_to_metric(req.zero_distance),
        zero_poi_vertical_m: u.small_len_to_metric(req.zero_poi_vertical),
        zero_poi_horizontal_m: u.small_len_to_metric(req.zero_poi_horizontal),
        sight_offset_lateral_m: u.small_len_to_metric(req.sight_offset_lateral),
        temperature_c: u.temperature_to_metric(temperature),
        pressure_hpa: u.pressure_to_metric(pressure),
        end_m: u.distance_to_metric(req.end),
        sample_m: u.distance_to_metric(req.step),
        elevation_click,
        windage_click,
        windage_unit,
        drag_model: req.drag_model.into(),
        u,
    })
}

fn solve_zero(req: &CardRequestV1, r: &Resolved) -> Result<f64, CardServiceError> {
    let zero_inputs = BallisticInputs {
        bc_value: req.ballistic_coefficient,
        bc_type: r.drag_model,
        bullet_mass: r.mass_kg,
        muzzle_velocity: r.velocity_m,
        bullet_diameter: r.diameter_m,
        bullet_length: crate::truing::fallback_bullet_length_m(r.diameter_m, r.mass_kg),
        sight_height: r.sight_height_m,
        use_rk4: true,
        zero_poi_vertical_m: r.zero_poi_vertical_m,
        zero_poi_horizontal_m: r.zero_poi_horizontal_m,
        sight_offset_lateral_m: r.sight_offset_lateral_m,
        ..Default::default()
    };
    let atmosphere = AtmosphericConditions {
        temperature: r.temperature_c,
        pressure: r.pressure_hpa,
        humidity: req.humidity,
        altitude: req.altitude,
    };
    crate::calculate_zero_angle_with_conditions(
        zero_inputs,
        r.zero_distance_m,
        r.sight_height_m,
        WindConditions::default(),
        atmosphere,
    )
    .map_err(|e| CardServiceError::ZeroSolve(e.to_string()))
}

#[allow(clippy::too_many_arguments)]
fn sampled(
    req: &CardRequestV1,
    r: &Resolved,
    wind_speed_m: f64,
    wind_direction_deg: f64,
    zero_angle: f64,
) -> Result<Vec<crate::trajectory_sampling::TrajectorySample>, CardServiceError> {
    run_sampled_trajectory(
        r.velocity_m,
        req.ballistic_coefficient,
        r.mass_kg,
        r.diameter_m,
        r.drag_model,
        r.sight_height_m,
        r.temperature_c,
        r.pressure_hpa,
        req.humidity,
        req.altitude,
        wind_speed_m,
        wind_direction_deg,
        r.end_m * 1.1,
        r.sample_m,
        zero_angle,
        None,
        None,
        None,
        r.zero_poi_vertical_m,
        r.zero_poi_horizontal_m,
        r.sight_offset_lateral_m,
        Some(r.zero_distance_m),
    )
    .map_err(|e| CardServiceError::Trajectory(e.to_string()))
}

fn nearest(
    samples: &[crate::trajectory_sampling::TrajectorySample],
    range_m: f64,
) -> Option<&crate::trajectory_sampling::TrajectorySample> {
    samples.iter().min_by(|a, b| {
        (a.distance_m - range_m)
            .abs()
            .partial_cmp(&(b.distance_m - range_m).abs())
            .unwrap()
    })
}

fn units_block(req: &CardRequestV1, windage_unit: AdjustmentUnit) -> CardUnitsBlockV1 {
    let imperial = req.units == CardUnits::Imperial;
    CardUnitsBlockV1 {
        distance: if imperial { "yd" } else { "m" },
        velocity: if imperial { "fps" } else { "m/s" },
        energy: if imperial { "ft-lb" } else { "J" },
        drop: if imperial { "in" } else { "mm" },
        wind_speed: if imperial { "mph" } else { "m/s" },
        elevation_adjustment: adjustment_unit_label(req.adjustment_unit),
        windage_adjustment: adjustment_unit_label(windage_unit),
    }
}

/// Come-ups: elevation dial per range plus the incremental come-up between rows.
/// Mirrors the CLI's `come-ups` command row-for-row.
pub fn come_ups_v1(req: &CardRequestV1) -> Result<CardResponseV1, CardServiceError> {
    let r = resolve(req)?;
    let zero_angle = solve_zero(req, &r)?;
    let samples = sampled(
        req,
        &r,
        r.u.wind_to_metric(req.wind_speed),
        req.wind_direction_deg,
        zero_angle,
    )?;

    let mut rows = Vec::new();
    let mut prev_drop_adj: f64 = 0.0;
    let mut current_range = req.start;
    while current_range <= req.end + 0.1 {
        let range_m = r.u.distance_to_metric(current_range);
        if let Some(sample) = nearest(&samples, range_m) {
            if (sample.distance_m - range_m).abs() < r.sample_m * 1.5 {
                let drop_yd = r.u.distance_from_metric(sample.drop_m);
                let range_display = r.u.distance_from_metric(sample.distance_m);
                let drop_adj = adjustment_display(
                    drop_yd,
                    range_display,
                    req.adjustment_unit,
                    r.elevation_click,
                    req.zero_set_elevation_bias_mil,
                    req.elevation_cf,
                )
                .value;
                rows.push(CardRowV1 {
                    range: current_range,
                    drop_linear: None,
                    drop_adj: Some(drop_adj),
                    come_up: Some(drop_adj - prev_drop_adj),
                    wind_linear: None,
                    wind_adj: None,
                    velocity: Some(r.u.velocity_from_metric(sample.velocity_mps)),
                    energy: Some(r.u.energy_from_metric(sample.energy_j)),
                    time: Some(sample.time_s),
                    wind_columns: Vec::new(),
                });
                prev_drop_adj = drop_adj;
            }
        }
        current_range += req.step;
    }

    Ok(CardResponseV1 {
        schema_version: CARD_SCHEMA_VERSION_V1,
        kind: "come_ups",
        zero_distance: req.zero_distance,
        units: units_block(req, req.adjustment_unit),
        wind_speeds: Vec::new(),
        wind_angles_deg: Vec::new(),
        rows,
        extra_angle_rows: Vec::new(),
    })
}

/// Range table: drop and wind on possibly different display axes.
/// Mirrors the CLI's `range-table` command row-for-row (two solves: with wind for
/// drift, without wind for the pure elevation axis).
pub fn range_table_v1(req: &CardRequestV1) -> Result<CardResponseV1, CardServiceError> {
    let r = resolve(req)?;
    let zero_angle = solve_zero(req, &r)?;
    let wind_samples = sampled(
        req,
        &r,
        r.u.wind_to_metric(req.wind_speed),
        req.wind_direction_deg,
        zero_angle,
    )?;
    let no_wind_samples = sampled(req, &r, 0.0, 0.0, zero_angle)?;

    let mut rows = Vec::new();
    let mut current_range = req.start;
    while current_range <= req.end + 0.1 {
        let range_m = r.u.distance_to_metric(current_range);
        if let (Some(nw), Some(w)) = (nearest(&no_wind_samples, range_m), nearest(&wind_samples, range_m)) {
            if (nw.distance_m - range_m).abs() < r.sample_m * 1.5 {
                let range_display = r.u.distance_from_metric(nw.distance_m);
                let drop_yd = r.u.distance_from_metric(nw.drop_m);
                let drop_adj = adjustment_display(
                    drop_yd,
                    range_display,
                    req.adjustment_unit,
                    r.elevation_click,
                    req.zero_set_elevation_bias_mil,
                    req.elevation_cf,
                )
                .value;
                let drift_yd = r.u.distance_from_metric(w.wind_drift_m);
                let wind_adj = windage_adjustment_display(
                    drift_yd,
                    range_display,
                    r.windage_unit,
                    r.windage_click,
                    req.zero_set_windage_bias_mil,
                    req.windage_cf,
                )
                .value;
                rows.push(CardRowV1 {
                    range: current_range,
                    drop_linear: Some(r.u.drop_linear_from_metric(nw.drop_m)),
                    drop_adj: Some(drop_adj),
                    come_up: None,
                    wind_linear: Some(r.u.drop_linear_from_metric(w.wind_drift_m)),
                    wind_adj: Some(wind_adj),
                    velocity: Some(r.u.velocity_from_metric(nw.velocity_mps)),
                    energy: Some(r.u.energy_from_metric(nw.energy_j)),
                    time: Some(nw.time_s),
                    wind_columns: Vec::new(),
                });
            }
        }
        current_range += req.step;
    }

    Ok(CardResponseV1 {
        schema_version: CARD_SCHEMA_VERSION_V1,
        kind: "range_table",
        zero_distance: req.zero_distance,
        units: units_block(req, r.windage_unit),
        wind_speeds: Vec::new(),
        wind_angles_deg: Vec::new(),
        rows,
        extra_angle_rows: Vec::new(),
    })
}

/// Wind card: drift matrix, ranges x wind speeds, per wind-FROM angle.
/// Mirrors the CLI's `wind-card` command cell-for-cell. The windage axis uses
/// `adjustment_unit` (the CLI's wind-card has a single unit for the matrix).
pub fn wind_card_v1(req: &CardRequestV1) -> Result<CardResponseV1, CardServiceError> {
    if req.wind_speeds.is_empty() {
        return Err(CardServiceError::InvalidRequest(
            "wind_card requires at least one entry in wind_speeds".into(),
        ));
    }
    let r = resolve(req)?;
    let zero_angle = solve_zero(req, &r)?;

    let angles: Vec<f64> = if req.wind_angles_deg.is_empty() {
        vec![90.0]
    } else {
        req.wind_angles_deg.clone()
    };

    let mut ranges: Vec<f64> = Vec::new();
    let mut current = req.start;
    while current <= req.end + 0.1 {
        ranges.push(current);
        current += req.step;
    }

    let mut per_angle_rows: Vec<Vec<CardRowV1>> = Vec::new();
    for &angle_deg in &angles {
        let mut all_drifts: Vec<Vec<f64>> = vec![Vec::new(); ranges.len()];
        for &ws in &req.wind_speeds {
            let samples = sampled(req, &r, r.u.wind_to_metric(ws), angle_deg, zero_angle)?;
            for (ri, &range_display) in ranges.iter().enumerate() {
                let range_m = r.u.distance_to_metric(range_display);
                let drift_adj = match nearest(&samples, range_m) {
                    Some(sample) if (sample.distance_m - range_m).abs() < r.sample_m * 1.5 => {
                        let drift_yd = r.u.distance_from_metric(sample.wind_drift_m);
                        adjustment_display(
                            drift_yd,
                            range_display,
                            req.adjustment_unit,
                            r.windage_click,
                            req.zero_set_windage_bias_mil,
                            req.windage_cf,
                        )
                        .value
                    }
                    _ => 0.0,
                };
                all_drifts[ri].push(drift_adj);
            }
        }
        per_angle_rows.push(
            ranges
                .iter()
                .enumerate()
                .map(|(i, &range)| CardRowV1 {
                    range,
                    drop_linear: None,
                    drop_adj: None,
                    come_up: None,
                    wind_linear: None,
                    wind_adj: None,
                    velocity: None,
                    energy: None,
                    time: None,
                    wind_columns: all_drifts[i].clone(),
                })
                .collect(),
        );
    }

    let first = per_angle_rows.remove(0);
    Ok(CardResponseV1 {
        schema_version: CARD_SCHEMA_VERSION_V1,
        kind: "wind_card",
        zero_distance: req.zero_distance,
        units: units_block(req, req.adjustment_unit),
        wind_speeds: req.wind_speeds.clone(),
        wind_angles_deg: angles,
        rows: first,
        extra_angle_rows: per_angle_rows,
    })
}
