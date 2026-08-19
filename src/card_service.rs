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
use crate::{AtmosphericConditions, BallisticInputs, BCSegmentData, DragModel, WindConditions};

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

    /// Explicit velocity-banded BC schedule, mirroring the CLI's repeatable
    /// `--bc-segment VMIN:VMAX:BC`. Velocities are in the request's `units` velocity
    /// unit (fps imperial / m/s metric), exactly like the CLI flag. When supplied it
    /// wins over `bc5d_table_path`, and — CLI parity — the scalar
    /// `ballistic_coefficient` remains the interior-gap fallback unchanged. The
    /// schedule feeds BOTH the zero solve and every sampled trajectory.
    #[serde(default)]
    pub bc_segments: Option<Vec<CardBcSegmentV1>>,
    /// Filesystem path to a caliber-specific BC5D correction table
    /// (`bc5d_<caliber>.bin`, the exact format the CLI's `--bc-table-dir`
    /// consumes). The table is CRC-verified, a velocity-keyed BC schedule is
    /// generated for this load (same ladder as the CLI), and the scalar BC gets the
    /// table's muzzle correction as the interior-gap fallback. Ignored when
    /// `bc_segments` is supplied. Not available on wasm32 builds (no filesystem);
    /// there it is rejected with an invalid-request error.
    #[serde(default)]
    pub bc5d_table_path: Option<String>,

    /// Presentation-only options for the printable PDF dope card (`card.pdf`). Nothing in
    /// here touches the ballistics: the numbers in the PDF are the rows
    /// [`range_table_v1`] returns, computed by the same call.
    ///
    /// It lives on the shared request — rather than in a `card.pdf`-only wrapper type — so
    /// an app stores ONE `CardRequestV1` per saved card and replays that exact document
    /// against either surface: `card.range_table` for the screen, `card.pdf` for the
    /// printout. The on-screen commands ignore this block.
    ///
    /// Deliberately NOT `#[cfg(feature = "pdf")]`, even though everything that consumes it
    /// is: `CardRequestV1` denies unknown fields, so gating the field would make a
    /// pdf-less build REJECT a stored request that carries it, breaking exactly the
    /// round-trip the field exists to support. A pdf-less build parses it and never reads
    /// it; only the `card.pdf` command itself disappears.
    #[serde(default)]
    pub pdf: Option<PdfCardOptionsV1>,
}

/// Presentation knobs for the PDF dope card: the header/footer labels the printed card
/// carries and the table's font size. The ballistics axes it prints (elevation unit,
/// windage unit, click graduations, tracking CFs, zero-set biases) are NOT here — they are
/// already fields on [`CardRequestV1`], shared with the on-screen card, precisely so the
/// two cannot disagree about what a row means.
///
/// Every field is optional and defaults to the CLI dope card's own default.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PdfCardOptionsV1 {
    /// Leading text of header line 1 — the CLI's rifle/profile name. Default
    /// `"Dope Card"`. Long headers are truncated to fit the page.
    #[serde(default)]
    pub title: Option<String>,
    /// Header `Loc:` label. Default empty.
    #[serde(default)]
    pub location: Option<String>,
    /// Footer `Powder:` label. Default empty.
    #[serde(default)]
    pub powder: Option<String>,
    /// Footer `Bullet:` label. Default empty.
    #[serde(default)]
    pub bullet: Option<String>,
    /// Crossing-target speed for the Lead column, in the request's `units` speed unit
    /// (mph imperial / m/s metric) — the same convention as the CLI's `--target-speed`
    /// (MBA-1325). The lead is `speed * time-of-flight` for a full-value 90° crossing,
    /// held on the windage axis, exactly as `trajectory -o pdf` computes it.
    ///
    /// Absent means "this card carries no lead data": the Lead column renders as em-dashes
    /// rather than a column of confident-looking zeroes. An explicit `0.0` is a different
    /// statement — a stationary target, whose lead genuinely is zero — and prints zeroes,
    /// matching the CLI's `--target-speed 0` default.
    #[serde(default)]
    pub target_speed: Option<f64>,
    /// Data-table font scale. Must be finite and within
    /// [`crate::pdf_dope_card::FONT_SCALE_RANGE`] (0.5..=3.0); out-of-band values are
    /// rejected rather than silently clamped, so what a stored request asks for is what it
    /// gets. Mutually exclusive with `font_preset`. Default 1.0.
    #[serde(default)]
    pub font_scale: Option<f64>,
    /// Named font scale: `"small"` (0.8), `"medium"` (1.0), `"large"` (1.4) — the CLI's
    /// `--font-preset` vocabulary, single-letter aliases included. An unrecognized preset
    /// is an error here (the CLI warns and falls back to medium; a stored request should
    /// not silently render at a size it did not ask for). Mutually exclusive with
    /// `font_scale`.
    #[serde(default)]
    pub font_preset: Option<String>,
    /// Render the data rows in bold. Default false.
    #[serde(default)]
    pub bold_data: bool,
}

/// One velocity band of an explicit velocity-keyed BC schedule
/// (`CardRequestV1::bc_segments`). `velocity_min`/`velocity_max` are in the
/// request's `units` velocity unit and must satisfy `velocity_min < velocity_max`;
/// `bc` is the BC (for the request's `drag_model`) that applies inside the band.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CardBcSegmentV1 {
    pub velocity_min: f64,
    pub velocity_max: f64,
    pub bc: f64,
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
    /// The scalar BC these rows were actually computed with: the request's published
    /// `ballistic_coefficient` unless a `bc5d_table_path` applied its muzzle correction, in
    /// which case it is the corrected value (the same scalar the solve and every sampled
    /// flight ran with).
    ///
    /// A saved card must be able to state the BC its numbers came from — the printed
    /// footer's `BC:` — long after the correction table it used has been replaced on disk.
    /// Without this on the response there is nothing for an app to store, and a reprint
    /// could only quote the BC the request nominated.
    pub bc_for_solve: f64,
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
    /// PDF rendering itself failed (font load/parse, or an empty row set). Gated with the
    /// only code path that can produce it, so a build without `pdf` sees an unchanged enum.
    #[cfg(feature = "pdf")]
    #[error("pdf generation failed: {0}")]
    Pdf(String),
    /// The card is too big to print ([`MAX_PDF_ROWS`] / [`MAX_PDF_PAGES`]). A resource
    /// refusal, not a failure: the request was well formed and the caller can be told
    /// exactly what is true of the card (how many rows, how many pages), so the message is
    /// carried verbatim rather than prefixed.
    #[cfg(feature = "pdf")]
    #[error("{0}")]
    TooLarge(String),
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
    /// Scalar BC actually fed to the zero solve and every sampled run. Equals the
    /// request's `ballistic_coefficient` unless a BC5D table applied its muzzle
    /// correction (CLI parity: the scalar is the schedule's interior-gap fallback).
    bc_for_solve: f64,
    /// Velocity-keyed BC schedule (velocities in fps, the unit the solver compares
    /// against), from explicit `bc_segments` or generated from `bc5d_table_path`.
    bc_segments_fps: Option<Vec<BCSegmentData>>,
}

fn resolve(req: &CardRequestV1) -> Result<Resolved, CardServiceError> {
    resolve_inner(req, true)
}

/// [`resolve`] WITHOUT the BC-schedule step — the request's declared axes, click
/// graduations and unit conversions, and nothing that touches the filesystem.
///
/// This is the resolve `card.pdf` performs when it prints caller-supplied rows: those rows
/// were computed elsewhere (they are the stored `card.range_table` response), so opening
/// `bc5d_table_path` here would do nothing but couple a reprint to whether that file is
/// still on the device — which is exactly the coupling printing stored rows exists to
/// break. `bc_for_solve` is left at the request's published BC and is not used for the
/// footer on that path (the stored card's own `bc_for_solve` is).
#[cfg(feature = "pdf")]
fn resolve_axes_only(req: &CardRequestV1) -> Result<Resolved, CardServiceError> {
    resolve_inner(req, false)
}

fn resolve_inner(req: &CardRequestV1, load_bc_schedule: bool) -> Result<Resolved, CardServiceError> {
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

    let (bc_for_solve, bc_segments_fps) = if load_bc_schedule {
        resolve_bc_schedule(req, imperial)?
    } else {
        (req.ballistic_coefficient, None)
    };

    Ok(Resolved {
        bc_for_solve,
        bc_segments_fps,
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

/// Resolve the velocity-keyed BC schedule for a card request, mirroring the CLI's
/// precedence exactly: explicit `bc_segments` win outright (and leave the scalar BC
/// untouched, like `--bc-segment` restoring the pre-table BC); otherwise a
/// `bc5d_table_path` generates the same segment ladder `--bc-table-dir` does and
/// applies the table's muzzle correction to the scalar fallback BC.
///
/// Returns `(scalar BC for the solve, optional fps-keyed schedule)`.
fn resolve_bc_schedule(
    req: &CardRequestV1,
    imperial: bool,
) -> Result<(f64, Option<Vec<BCSegmentData>>), CardServiceError> {
    if let Some(segments) = req.bc_segments.as_ref().filter(|s| !s.is_empty()) {
        // Display velocity -> fps: the exact factor the CLI's parse_bc_segment uses.
        let to_fps = if imperial { 1.0 } else { 3.280_839_895 };
        let mut converted = Vec::with_capacity(segments.len());
        for (index, segment) in segments.iter().enumerate() {
            if !segment.velocity_min.is_finite()
                || !segment.velocity_max.is_finite()
                || !segment.bc.is_finite()
            {
                return Err(CardServiceError::InvalidRequest(format!(
                    "bc_segments[{index}]: velocity_min, velocity_max, and bc must be finite"
                )));
            }
            if segment.velocity_min >= segment.velocity_max {
                return Err(CardServiceError::InvalidRequest(format!(
                    "bc_segments[{index}]: velocity_min must be < velocity_max"
                )));
            }
            if segment.bc <= 0.0 {
                return Err(CardServiceError::InvalidRequest(format!(
                    "bc_segments[{index}]: bc must be > 0"
                )));
            }
            converted.push(BCSegmentData {
                velocity_min: segment.velocity_min * to_fps,
                velocity_max: segment.velocity_max * to_fps,
                bc_value: segment.bc,
            });
        }
        return Ok((req.ballistic_coefficient, Some(converted)));
    }

    let Some(path) = req.bc5d_table_path.as_deref() else {
        return Ok((req.ballistic_coefficient, None));
    };

    #[cfg(target_arch = "wasm32")]
    {
        let _ = path;
        Err(CardServiceError::InvalidRequest(
            "bc5d_table_path is not supported on this target (no filesystem); supply \
             bc_segments instead"
                .into(),
        ))
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        // The table's CONTENT must be for this shot's caliber, not merely CRC-valid: a
        // wrong-caliber table silently biases every row of the card (see
        // Bc5dTable::ensure_caliber_matches), so a mismatch is refused here rather than
        // applied — and rather than quietly falling back to an uncorrected card, which
        // the caller could not distinguish from a corrected one.
        let diameter_in = if imperial { req.diameter } else { req.diameter / 25.4 };
        let table = crate::bc_table_5d::path_cache::load_verified_for_caliber(
            std::path::Path::new(path),
            diameter_in,
        )
        .map_err(|e| CardServiceError::InvalidRequest(format!("bc5d_table_path: {e}")))?;

        // The BC5D axes are grains + fps; convert from the request's units the same
        // way the CLI does. The v2 tables carry only G1/G7 planes — anything else is
        // typed as G1, matching the CLI/WASM coercion.
        let weight_grains = if imperial { req.mass } else { req.mass * 15.4324 };
        let muzzle_fps = if imperial {
            req.muzzle_velocity
        } else {
            req.muzzle_velocity * 3.280_839_895
        };
        let drag_type = if req.drag_model == DragModelV1::G7 { "G7" } else { "G1" };
        let base_bc = req.ballistic_coefficient;

        match table.generate_segments(base_bc, drag_type, weight_grains, Some(muzzle_fps)) {
            Some(segments) => {
                // CLI parity: the scalar BC becomes the muzzle-corrected fallback for
                // interior coverage gaps (main.rs: trued_bc = base * muzzle_correction).
                let fallback_bc =
                    table.get_effective_bc(weight_grains, base_bc, muzzle_fps, muzzle_fps, drag_type);
                Ok((fallback_bc, Some(segments)))
            }
            // A neutral table (every sampled cell ~= 1.0) carries no correction:
            // keep the published scalar BC, exactly like the CLI.
            None => Ok((base_bc, None)),
        }
    }
}

fn solve_zero(req: &CardRequestV1, r: &Resolved) -> Result<f64, CardServiceError> {
    let zero_inputs = BallisticInputs {
        bc_value: r.bc_for_solve,
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
        // The zero MUST be solved with the same velocity-keyed BC schedule the
        // sampled flights use (the 0.22.11 auto-zero lesson: a schedule that changes
        // early-flight drag would otherwise mis-zero every row).
        use_bc_segments: r.bc_segments_fps.is_some(),
        bc_segments_data: r.bc_segments_fps.clone(),
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
        r.bc_for_solve,
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
        // The same schedule as the zero solve above, on every card surface
        // (come-ups, range table, wind card) consistently.
        r.bc_segments_fps.clone(),
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
        bc_for_solve: r.bc_for_solve,
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
    // Discards the PDF-only extras; see `range_table_rows`.
    range_table_rows(req, None).map(|(card, _lead, _bc)| card)
}

/// The one range-table computation, shared by the on-screen `card.range_table` and the
/// printable `card.pdf`.
///
/// Returns `(the card, the Lead column, the scalar BC the solve ran with)`:
///
/// * The card is exactly what `card.range_table` serializes — `range_table_v1` is a
///   discarding wrapper, so the PDF's Range/Drop/Wind figures ARE the on-screen figures.
///   They are taken from here rather than recomputed because recomputing them is the one
///   way the printed card could ever disagree with the rows the shooter already read.
/// * The Lead column parallels `card.rows` one-to-one. `lead_target_speed` is in the
///   request's `units` speed unit; `None` (no lead requested) yields all `None`, which the
///   PDF renders as em-dashes rather than fake zeroes.
/// * The BC is `Resolved::bc_for_solve` — the published BC unless a BC5D table applied its
///   muzzle correction — which the PDF footer prints. The footer must state the BC the
///   numbers above it came from, not the one the request nominated.
fn range_table_rows(
    req: &CardRequestV1,
    lead_target_speed: Option<f64>,
) -> Result<(CardResponseV1, Vec<Option<f64>>, f64), CardServiceError> {
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
    let lead_speed_mps = lead_target_speed.map(|speed| r.u.wind_to_metric(speed));

    let mut rows = Vec::new();
    let mut lead_adj = Vec::new();
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
                // Lead is the windage-axis hold for a full-value 90° crossing target:
                // `speed * time-of-flight` (MBA-1287's shared moving-target math), on the
                // no-wind sample whose time this row already reports. It is a COMPONENT
                // hold composed on top of the wind dial, which already carries the
                // zero-set bias, so it stays bias-free (bias 0.0) while still being
                // divided by the windage tracking CF — the same treatment
                // `main.rs::dope_card_row_from_sample` gives the CLI's Lead column.
                lead_adj.push(lead_speed_mps.map(|speed_mps| {
                    let lead_display = r.u.distance_from_metric(
                        crate::lead_from_tof(speed_mps, 90.0, nw.time_s, nw.distance_m).lead_m,
                    );
                    windage_adjustment_display(
                        lead_display,
                        range_display,
                        r.windage_unit,
                        r.windage_click,
                        0.0,
                        req.windage_cf,
                    )
                    .value
                }));
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

    debug_assert_eq!(rows.len(), lead_adj.len(), "Lead column must parallel rows");
    Ok((
        CardResponseV1 {
            schema_version: CARD_SCHEMA_VERSION_V1,
            kind: "range_table",
            zero_distance: req.zero_distance,
            bc_for_solve: r.bc_for_solve,
            units: units_block(req, r.windage_unit),
            wind_speeds: Vec::new(),
            wind_angles_deg: Vec::new(),
            rows,
            extra_angle_rows: Vec::new(),
        },
        lead_adj,
        r.bc_for_solve,
    ))
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
        bc_for_solve: r.bc_for_solve,
        units: units_block(req, req.adjustment_unit),
        wind_speeds: req.wind_speeds.clone(),
        wind_angles_deg: angles,
        rows: first,
        extra_angle_rows: per_angle_rows,
    })
}

// ---------------------------------------------------------------------------------------
// Printable PDF dope card (`card.pdf`), feature `pdf`.
//
// Two ways in, one renderer:
//
// * REPRINT — the caller supplies the rows: `StoredCardV1` carries the stored
//   `card.range_table` response verbatim. Nothing is solved, no correction table is opened,
//   and the footer's BC and provenance come from that document. Only under this construction
//   is a reprint of a saved card actually a reprint: it cannot drift when the engine is
//   bumped, or when the BC5D table file at the stored path is overwritten in place by a
//   table-set refresh.
// * SOLVE — no rows supplied: `range_table_rows` runs, which is the same call
//   `card.range_table` makes (`range_table_v1` is a discarding wrapper over it). This is the
//   pre-existing behaviour, unchanged, and what the CLI-shaped caller gets.
//
// Neither path recomputes an adjustment, a unit conversion or a trajectory that the card it
// prints already performed.
// ---------------------------------------------------------------------------------------

/// The one card kind `card.pdf` can print, in the `kind` spelling
/// [`CardResponseV1::kind`] uses.
///
/// A come-ups card's Come-Up column and a wind card's swept drift matrix are not columns a
/// range-table card has, so a request for either is REFUSED rather than answered with a
/// range table — an `ok` response whose defining field was silently ignored is the defect
/// (a stored wind card exported as a range table asserted 0.0 drift on every row while the
/// screen showed up to -0.42 MIL).
#[cfg(feature = "pdf")]
pub const PDF_CARD_KIND: &str = "range_table";

/// Rows a printable card may carry.
///
/// Checked from the row count itself, before a document exists — the byte cap
/// (`bridge::MAX_PDF_BYTES`) can only refuse a 26 MB `Vec<u8>` that has already been built
/// and paginated. At ~0.5 KiB/row over the ~815 KiB font floor this bound keeps every
/// accepted card comfortably inside that cap, which remains the backstop for the other way
/// to make a huge document: header/footer labels, which are drawn verbatim on every page.
#[cfg(feature = "pdf")]
pub const MAX_PDF_ROWS: usize = 5_000;

/// Pages a printable card may run to. Bounds the one row set the row cap above cannot:
/// few enough rows, but at a font scale that fits very few of them per page.
#[cfg(feature = "pdf")]
pub const MAX_PDF_PAGES: usize = 60;

/// Where the printed rows came from — reported so a caller can verify it got a reprint
/// rather than a re-solve, which is otherwise indistinguishable from the document.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PdfRowSource {
    /// Solved by this call, exactly as `card.range_table` would for the same request.
    Solve,
    /// Supplied by the caller: a stored `card.range_table` response, printed as-is.
    StoredRows,
}

#[cfg(feature = "pdf")]
impl PdfRowSource {
    /// Wire spelling for the bridge's `source` field.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Solve => "solve",
            Self::StoredRows => "stored_rows",
        }
    }
}

/// A card whose rows already exist: print THESE numbers instead of solving.
///
/// `card` is the stored `card.range_table` response, pasted verbatim — an app keeps that
/// document for every saved card and shows it on screen, so handing the same bytes back is
/// what makes the paper and the screen the same card by construction rather than by
/// coincidence. The two provenance strings are printed in the footer so a card in a
/// shooter's pocket can be reconciled with a screen afterwards; both are optional, and an
/// absent or empty one prints nothing at all rather than a placeholder.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StoredCardV1 {
    /// The stored response, exactly as this engine emitted it.
    pub card: StoredCardResponseV1,
    /// Engine version that produced the rows (`engine_version` from the bridge envelope of
    /// the call that produced them). Absent/empty prints no engine line.
    #[serde(default)]
    pub engine_version: Option<String>,
    /// Correction-table version the rows were solved against — the published BC5D table-set
    /// version an app records at save time. Absent/empty prints no table line, which is the
    /// honest rendering of "this card used no correction table".
    #[serde(default)]
    pub bc5d_table_version: Option<String>,
}

/// A stored [`CardResponseV1`], as a deserializable mirror.
///
/// [`CardResponseV1`] is `Serialize` only (and its `units` block holds `&'static str`s), so
/// a stored response is read back through this twin. It accepts the response document
/// unchanged, field for field: `kind` and `units` are load-bearing here — they say what the
/// rows ARE and what a row MEANS — while the descriptive scalars are optional so a card
/// stored by an older engine still prints.
///
/// **Deliberately NOT `deny_unknown_fields`** (unlike [`CardRequestV1`], where the caller
/// really is the author of the document). What arrives here is this engine's own output,
/// round-tripped through an app's storage — never hand-written input — so strictness buys
/// no validation and costs a hard cross-version break: the FIRST field ever added to
/// `CardResponseV1` would make every card saved by the newer platform completely
/// unexportable on the older one. The two apps carry independent engine pins and ship
/// through separate stores, so that staggering is the normal case, not a downgrade. An
/// unknown field is therefore ignored, exactly as the apps' own decoders ignore it, and the
/// stored cells still reach the paper. The same reasoning applies to
/// [`StoredCardUnitsBlockV1`] and [`StoredCardRowV1`] below.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct StoredCardResponseV1 {
    /// Which card this is. Anything but [`PDF_CARD_KIND`] is refused.
    pub kind: String,
    /// Column labels and the range unit of the stored rows. The printed headings are
    /// these, not the request's, because these are what the rows were computed in.
    pub units: StoredCardUnitsBlockV1,
    /// The rows to print, in order.
    pub rows: Vec<StoredCardRowV1>,
    #[serde(default)]
    pub schema_version: Option<u32>,
    #[serde(default)]
    pub zero_distance: Option<f64>,
    /// The BC the stored rows were computed with; printed in the footer. Absent (a card
    /// saved by an engine before [`CardResponseV1::bc_for_solve`] existed) falls back to the
    /// request's published `ballistic_coefficient`.
    #[serde(default)]
    pub bc_for_solve: Option<f64>,
    /// Present only on a wind card, and therefore a refusal here.
    #[serde(default)]
    pub wind_speeds: Vec<f64>,
    #[serde(default)]
    pub wind_angles_deg: Vec<f64>,
    #[serde(default)]
    pub extra_angle_rows: Vec<Vec<StoredCardRowV1>>,
}

/// The `units` block of a stored response (see [`CardUnitsBlockV1`], which this mirrors).
/// Unknown fields are ignored, for the reason [`StoredCardResponseV1`] gives.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct StoredCardUnitsBlockV1 {
    /// `"yd"` / `"m"` — the unit the stored `range` values are in.
    pub distance: String,
    /// Drop column label ("MIL"/"MOA"/"SMOA"/"IPHY"/"CLICKS").
    pub elevation_adjustment: String,
    /// Wind and Lead column label; may differ from the elevation label (MBA-1410).
    pub windage_adjustment: String,
    #[serde(default)]
    pub velocity: Option<String>,
    #[serde(default)]
    pub energy: Option<String>,
    #[serde(default)]
    pub drop: Option<String>,
    #[serde(default)]
    pub wind_speed: Option<String>,
}

/// One stored row (see [`CardRowV1`], which this mirrors). Only `range` is required; a
/// column the stored card did not carry stays `None` and prints as an em-dash rather than a
/// fabricated zero. Unknown fields are ignored, for the reason [`StoredCardResponseV1`]
/// gives.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone, Default, PartialEq, Deserialize)]
pub struct StoredCardRowV1 {
    pub range: f64,
    #[serde(default)]
    pub drop_linear: Option<f64>,
    #[serde(default)]
    pub drop_adj: Option<f64>,
    #[serde(default)]
    pub come_up: Option<f64>,
    #[serde(default)]
    pub wind_linear: Option<f64>,
    #[serde(default)]
    pub wind_adj: Option<f64>,
    #[serde(default)]
    pub velocity: Option<f64>,
    #[serde(default)]
    pub energy: Option<f64>,
    /// Time of flight, seconds. The Lead column of a reprint is derived from this — pure
    /// arithmetic on a number the stored card already carries, not a new trajectory. A row
    /// without it prints an em-dash for Lead.
    #[serde(default)]
    pub time: Option<f64>,
    #[serde(default)]
    pub wind_columns: Vec<f64>,
}

/// A rendered PDF dope card: the document plus the facts a transport needs to describe it
/// without re-parsing the bytes.
#[cfg(feature = "pdf")]
#[derive(Debug, Clone)]
pub struct PdfCardV1 {
    pub pdf_bytes: Vec<u8>,
    /// Pages in `pdf_bytes`, from [`crate::pdf_dope_card::dope_card_page_count`] — the
    /// same function the generator paginated by, not a second estimate of it.
    pub page_count: usize,
    /// Rows printed. On the solve path this is the on-screen `card.range_table` row count;
    /// on the reprint path it is the stored row count, so a caller can confirm every stored
    /// row reached the paper.
    pub row_count: usize,
    /// Solved here, or printed from caller-supplied rows.
    pub source: PdfRowSource,
    /// Characters of the header title (`pdf.title`) the card's font could not draw —
    /// distinct, in order of first use, and empty when the title printed in full.
    ///
    /// Not an error: the document is complete and every ROW printed. But a title is how a
    /// shooter tells one card from another, and Liberation Sans covers no CJK, Arabic,
    /// Hebrew, Thai or emoji, so a name in any of those used to vanish from the paper with
    /// `ok: true` and nothing said. Each such character prints as
    /// [`crate::pdf_dope_card::UNPRINTABLE_SUBSTITUTE`] and is named here, so a caller that
    /// accepts any non-empty card name (both apps do) can warn at export time instead of
    /// handing over an unidentifiable card.
    pub unprintable_title_chars: String,
}

/// Resolve the effective table font scale from the mutually exclusive `font_scale` /
/// `font_preset` options.
#[cfg(feature = "pdf")]
fn resolve_font_scale(opts: &PdfCardOptionsV1) -> Result<f32, CardServiceError> {
    use crate::pdf_dope_card::{FontSizePreset, FONT_SCALE_RANGE};

    match (opts.font_scale, opts.font_preset.as_deref()) {
        (Some(_), Some(_)) => Err(CardServiceError::InvalidRequest(
            "pdf.font_scale and pdf.font_preset are mutually exclusive; supply one".into(),
        )),
        (Some(scale), None) => {
            let scale = scale as f32;
            if !scale.is_finite() || !FONT_SCALE_RANGE.contains(&scale) {
                return Err(CardServiceError::InvalidRequest(format!(
                    "pdf.font_scale must be finite and within {}..={}, got {}",
                    FONT_SCALE_RANGE.start(),
                    FONT_SCALE_RANGE.end(),
                    opts.font_scale.unwrap_or(f64::NAN)
                )));
            }
            Ok(scale)
        }
        (None, Some(preset)) => FontSizePreset::from_str(preset).map(|p| p.scale()).ok_or_else(|| {
            CardServiceError::InvalidRequest(format!(
                "pdf.font_preset '{preset}' is not one of small, medium, large"
            ))
        }),
        (None, None) => Ok(1.0),
    }
}

/// The rows to print plus everything the header/footer states about where they came from.
/// Filled by exactly one of the two paths in [`pdf_card_v1`], which then share one
/// renderer, so a reprint and a solve cannot be laid out or labelled differently.
#[cfg(feature = "pdf")]
struct RowsToPrint {
    rows: Vec<crate::card::CardRow>,
    source: PdfRowSource,
    /// Footer `BC:` — the BC these rows were computed with.
    bc: f64,
    /// Footer `Engine:` — empty prints nothing.
    engine_version: String,
    /// Footer `Table:` — empty prints nothing.
    table_version: String,
    /// Drop column label.
    elevation_unit_label: String,
    /// Wind and Lead column label.
    windage_unit_label: String,
}

/// Turn a stored response's rows into printable rows, deriving only the Lead column.
///
/// The Drop and Wind cells are the stored values, untouched. Lead is `target_speed x stored
/// ToF` held on the windage axis — the same [`crate::lead_from_tof`] arithmetic
/// `range_table_rows` performs, applied to the time of flight the stored row already
/// carries. No trajectory, no zero solve, no correction table.
#[cfg(feature = "pdf")]
fn stored_rows_to_print(
    stored: &StoredCardResponseV1,
    req: &CardRequestV1,
    r: &Resolved,
    target_speed: Option<f64>,
) -> Vec<crate::card::CardRow> {
    use crate::card::CardRow;

    let lead_speed_mps = target_speed.map(|speed| r.u.wind_to_metric(speed));
    stored
        .rows
        .iter()
        .map(|row| CardRow {
            range: row.range,
            drop_linear: row.drop_linear,
            drop_adj: row.drop_adj,
            come_up: row.come_up,
            wind_linear: row.wind_linear,
            wind_adj: row.wind_adj,
            velocity: row.velocity,
            energy: row.energy,
            time: row.time,
            lead_adj: lead_speed_mps.zip(row.time).map(|(speed_mps, tof_s)| {
                let range_m = r.u.distance_to_metric(row.range);
                let lead_display =
                    r.u.distance_from_metric(crate::lead_from_tof(speed_mps, 90.0, tof_s, range_m).lead_m);
                // Bias-free (it composes on top of the wind dial, which carries the
                // zero-set bias) but still divided by the windage tracking CF — the exact
                // treatment `range_table_rows` and the CLI's Lead column give it.
                windage_adjustment_display(
                    lead_display,
                    row.range,
                    r.windage_unit,
                    r.windage_click,
                    0.0,
                    req.windage_cf,
                )
                .value
            }),
            wind_columns: row.wind_columns.clone(),
        })
        .collect()
}

/// Reject a stored card that is not the range table it must be, cell values included.
#[cfg(feature = "pdf")]
fn validate_stored_card(
    stored: &StoredCardV1,
    req: &CardRequestV1,
    r: &Resolved,
) -> Result<(), CardServiceError> {
    let card = &stored.card;
    if card.kind != PDF_CARD_KIND {
        return Err(CardServiceError::InvalidRequest(format!(
            "card.pdf prints a {PDF_CARD_KIND} card; the stored card's kind is '{}'. Its \
             columns are not a range table's, so it is refused rather than reprinted as one",
            card.kind
        )));
    }
    // The same "this is not a range table" test the REQUEST gets (see the `wind_speeds` /
    // `wind_angles_deg` loop in `pdf_card_v1`), field for field. The two halves used to
    // disagree — `wind_angles_deg` was refused on the request and accepted here — and the
    // stored half is the one that has to hold when a future kind starts emitting it.
    for (field, present) in [
        ("wind_speeds", !card.wind_speeds.is_empty()),
        ("wind_angles_deg", !card.wind_angles_deg.is_empty()),
        ("extra_angle_rows", !card.extra_angle_rows.is_empty()),
    ] {
        if present {
            return Err(CardServiceError::InvalidRequest(format!(
                "the stored card carries {field}, a wind card's defining field, so it is not the \
                 {PDF_CARD_KIND} response it claims to be"
            )));
        }
    }
    if card.rows.is_empty() {
        return Err(CardServiceError::InvalidRequest(
            "the stored card has no rows; an empty dope card is not a document".into(),
        ));
    }

    // A row means something only in a unit. If the stored labels and the request's own axes
    // disagree, the two are not the same card — printing anyway would put one unit's numbers
    // under another unit's heading, which is precisely the failure the shared request shape
    // exists to prevent.
    let elevation = adjustment_unit_label(req.adjustment_unit);
    let windage = adjustment_unit_label(r.windage_unit);
    let distance = if req.units == CardUnits::Imperial { "yd" } else { "m" };
    for (field, stored_label, request_label) in [
        ("units.distance", card.units.distance.as_str(), distance),
        (
            "units.elevation_adjustment",
            card.units.elevation_adjustment.as_str(),
            elevation.as_str(),
        ),
        (
            "units.windage_adjustment",
            card.units.windage_adjustment.as_str(),
            windage.as_str(),
        ),
    ] {
        if stored_label != request_label {
            return Err(CardServiceError::InvalidRequest(format!(
                "the stored card's {field} is '{stored_label}' but this request's own axes say \
                 '{request_label}' — the stored rows and the request are not the same card"
            )));
        }
    }

    // A non-finite cell would print "NaN"/"inf" on a field card.
    for (index, row) in card.rows.iter().enumerate() {
        for (field, value) in [
            ("range", Some(row.range)),
            ("drop_adj", row.drop_adj),
            ("wind_adj", row.wind_adj),
            ("time", row.time),
        ] {
            if value.is_some_and(|v| !v.is_finite()) {
                return Err(CardServiceError::InvalidRequest(format!(
                    "the stored card's rows[{index}].{field} is not a finite number"
                )));
            }
        }
    }
    Ok(())
}

/// Render the printable PDF dope card for a card request.
///
/// `stored` decides where the numbers come from, and it is the whole point of this surface:
///
/// * `Some(card)` — REPRINT. The rows are the caller's: the stored `card.range_table`
///   response for this same request. Nothing is solved, `bc5d_table_path` is never opened,
///   and the footer's BC, engine version and table version are the stored card's. A saved
///   card therefore reprints identically after an engine bump, after the correction table at
///   the stored path is overwritten in place, and even after that file is deleted.
/// * `None` — SOLVE. `range_table_rows` runs, i.e. literally the rows [`range_table_v1`]
///   would return for this request, and the footer states THIS build's version. Unchanged
///   behaviour for a caller that has no stored rows.
///
/// Either way this prints a [`PDF_CARD_KIND`] card and nothing else: a request carrying a
/// wind card's `wind_speeds`/`wind_angles_deg`, or a stored card of another kind, is refused.
/// Both paths map their rows onto the CLI's Range/Drop/Wind/Lead dope card, plus the Lead
/// column that `CardRequestV1::pdf`'s `target_speed` asks for. The Range column is
/// denominated in the stored/requested distance unit (yards imperial / metres metric),
/// unlike `trajectory -o pdf`, whose dope card is always yards.
///
/// The header/footer block is always imperial, matching both CLI PDF call sites: a metric
/// request's velocity/temperature/pressure/altitude/wind/weight are converted for display
/// only. The `Solver:` label reports this build (`online` with the `online` feature,
/// otherwise `offline`) and the timestamp is generation time, so neither is caller-settable
/// — and neither is a number a shooter dials.
#[cfg(feature = "pdf")]
pub fn pdf_card_v1(
    req: &CardRequestV1,
    stored: Option<&StoredCardV1>,
) -> Result<PdfCardV1, CardServiceError> {
    use crate::pdf_dope_card::{
        calculate_density_altitude, dope_card_page_count, generate_dope_card_pdf, DopeCardConfig,
        RangeUnit,
    };

    let opts = req.pdf.clone().unwrap_or_default();
    // Validate presentation options BEFORE spending a zero solve and two trajectories on a
    // request that cannot be rendered anyway.
    let font_scale = resolve_font_scale(&opts)?;
    // A non-finite crossing speed would print "inf"/"NaN" in every Lead cell rather than
    // failing. (A NEGATIVE speed is accepted, matching the CLI: it reads as the mover
    // crossing the other way and simply flips the sign of the hold.)
    if opts.target_speed.is_some_and(|speed| !speed.is_finite()) {
        return Err(CardServiceError::InvalidRequest(
            "pdf.target_speed must be finite".into(),
        ));
    }

    // A wind card's defining fields cannot be honoured by a range-table PDF. Refuse them
    // instead of returning a document whose Wind column contradicts the screen.
    for (field, present) in [
        ("wind_speeds", !req.wind_speeds.is_empty()),
        ("wind_angles_deg", !req.wind_angles_deg.is_empty()),
    ] {
        if present {
            return Err(CardServiceError::InvalidRequest(format!(
                "card.pdf prints a {PDF_CARD_KIND} card; this request carries {field}, a wind \
                 card's defining field, which a range-table PDF cannot show and would silently \
                 ignore. There is no wind_card or come_ups PDF in this build"
            )));
        }
    }

    let to_print = match stored {
        Some(stored) => {
            // Axes and click graduations only: no BC schedule, so no correction table is
            // opened for a card whose numbers are already decided.
            let r = resolve_axes_only(req)?;
            validate_stored_card(stored, req, &r)?;
            RowsToPrint {
                rows: stored_rows_to_print(&stored.card, req, &r, opts.target_speed),
                source: PdfRowSource::StoredRows,
                bc: stored.card.bc_for_solve.unwrap_or(req.ballistic_coefficient),
                engine_version: stored.engine_version.clone().unwrap_or_default(),
                table_version: stored.bc5d_table_version.clone().unwrap_or_default(),
                elevation_unit_label: stored.card.units.elevation_adjustment.clone(),
                windage_unit_label: stored.card.units.windage_adjustment.clone(),
            }
        }
        None => {
            let (card, lead_adj, bc_for_solve) = range_table_rows(req, opts.target_speed)?;
            if card.rows.is_empty() {
                // `generate_dope_card_pdf` refuses an empty row set, and rightly: an empty
                // dope card is not a document. Name the cause the caller can act on instead.
                return Err(CardServiceError::InvalidRequest(format!(
                    "no card rows in {}..={} — the trajectory does not reach this range domain, \
                     or step is coarser than the samples",
                    req.start, req.end
                )));
            }
            RowsToPrint {
                rows: card
                    .rows
                    .iter()
                    .zip(&lead_adj)
                    .map(|(row, lead)| crate::card::CardRow {
                        range: row.range,
                        drop_linear: row.drop_linear,
                        drop_adj: row.drop_adj,
                        come_up: row.come_up,
                        wind_linear: row.wind_linear,
                        wind_adj: row.wind_adj,
                        velocity: row.velocity,
                        energy: row.energy,
                        time: row.time,
                        lead_adj: *lead,
                        wind_columns: row.wind_columns.clone(),
                    })
                    .collect(),
                source: PdfRowSource::Solve,
                bc: bc_for_solve,
                // No table version is knowable from a path, so none is claimed; the engine
                // version is this build's, which is what produced these rows.
                engine_version: env!("CARGO_PKG_VERSION").to_string(),
                table_version: String::new(),
                // The labels `card.range_table`'s `units` block reports for this request:
                // Drop in the elevation unit, Wind AND Lead in the (possibly different,
                // MBA-1410) windage unit.
                elevation_unit_label: card.units.elevation_adjustment.clone(),
                windage_unit_label: card.units.windage_adjustment.clone(),
            }
        }
    };

    // Refuse on the row/page count now that the rows are known — before a document is built.
    // The byte cap downstream can only measure a document it already paid for.
    let page_count = dope_card_page_count(to_print.rows.len(), font_scale);
    if to_print.rows.len() > MAX_PDF_ROWS || page_count > MAX_PDF_PAGES {
        return Err(CardServiceError::TooLarge(format!(
            "this card has too many rows to print: {} rows, {page_count} pages (the limits are \
             {MAX_PDF_ROWS} rows and {MAX_PDF_PAGES} pages)",
            to_print.rows.len()
        )));
    }

    let imperial = req.units == CardUnits::Imperial;
    // The atmosphere defaults `resolve()` applies, restated for the header so the printed
    // conditions are the ones the solve used rather than blanks.
    let temperature = req.temperature.unwrap_or(if imperial { 59.0 } else { 15.0 });
    let pressure = req.pressure.unwrap_or(if imperial { 29.92 } else { 1013.25 });
    // inHg <-> hPa via the CLI dope card's own factor (main.rs uses 33.8639 on both PDF
    // call sites), NOT pdf_dope_card::INHG_TO_HPA — matching the shipped header exactly.
    let pressure_inhg = if imperial { pressure } else { pressure / 33.8639 };
    let pressure_hpa = if imperial { pressure * 33.8639 } else { pressure };
    let temperature_f = if imperial { temperature } else { temperature * 9.0 / 5.0 + 32.0 };
    // NOTE: `CardRequestV1::altitude` is fed to the solve unconverted (see `solve_zero` /
    // `sampled`), i.e. it is METRES in both unit systems — the one field in this request
    // that does not follow the module's units convention. The header reports the altitude
    // the solve actually used, so it converts from metres regardless of `units`.
    let altitude_ft = req.altitude / 0.3048;

    let config = DopeCardConfig {
        rifle_name: opts.title.clone().unwrap_or_else(|| "Dope Card".to_string()),
        location: opts.location.clone().unwrap_or_default(),
        density_altitude_ft: calculate_density_altitude(altitude_ft, pressure_inhg, temperature_f),
        pressure_inhg,
        pressure_hpa,
        temperature_f,
        altitude_ft,
        wind_speed_mph: if imperial { req.wind_speed } else { req.wind_speed / 0.44704 },
        // Absent target speed prints 0 here (there is no lead to state) while the Lead
        // column itself stays em-dashed; an explicit 0.0 prints the same 0 and zeroes.
        target_speed_mph: match opts.target_speed {
            Some(speed) if imperial => speed,
            Some(speed) => speed / 0.44704,
            None => 0.0,
        },
        solver_mode: if cfg!(feature = "online") { "online".to_string() } else { "offline".to_string() },
        powder: opts.powder.clone().unwrap_or_default(),
        bullet: opts.bullet.clone().unwrap_or_default(),
        weight_gr: if imperial { req.mass } else { req.mass * 15.4324 },
        // The BC these rows came from: the stored card's own, or this solve's (which is the
        // muzzle-corrected value when a BC5D table applied one).
        bc: to_print.bc,
        drag_model: DragModel::from(req.drag_model).to_string(),
        velocity_fps: if imperial { req.muzzle_velocity } else { req.muzzle_velocity / 0.3048 },
        font_scale,
        bold_data: opts.bold_data,
        // The card's own axes: Drop in the elevation unit, Wind AND Lead in the (possibly
        // different, MBA-1410) windage unit — from the document that owns the rows.
        elevation_unit_label: to_print.elevation_unit_label.clone(),
        windage_unit_label: to_print.windage_unit_label.clone(),
        // Provenance on the paper, so a printed card and a screen can be reconciled later.
        engine_version: to_print.engine_version.clone(),
        table_version: to_print.table_version.clone(),
    };

    let range_unit = if imperial { RangeUnit::Yards } else { RangeUnit::Meters };
    // Read BEFORE the document is built, off the same string the header will draw, so the
    // report describes this card's own title rather than a truncation of it.
    let unprintable_title_chars = crate::pdf_dope_card::unprintable_chars(&config.rifle_name);
    let pdf_bytes = generate_dope_card_pdf(&config, &to_print.rows, range_unit)
        .map_err(|e| CardServiceError::Pdf(e.to_string()))?;

    Ok(PdfCardV1 {
        page_count,
        row_count: to_print.rows.len(),
        source: to_print.source,
        pdf_bytes,
        unprintable_title_chars,
    })
}
