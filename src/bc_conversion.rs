//! Ballistic-coefficient conversion between standard reference drag families (MBA-1375).
//!
//! A reference-family BC appears in the point-mass retardation as `Cd_ref(Mach) / BC`.
//! Consequently, conversion at one Mach is an exact table-ratio operation, while fitting one
//! scalar BC across velocity bands is a one-parameter least-squares problem in reciprocal BC.
//! This module owns that shared, filesystem-free math so native and WASM front ends cannot drift.

use crate::drag::{get_drag_coefficient, reference_drag_table, DragTable};
use crate::{BCSegmentData, DragModel};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;

/// Schema version carried by every serializable conversion report.
pub const BC_CONVERSION_SCHEMA_VERSION_V1: u32 = 1;

/// A conversion curve must contain more than the legacy six-row placeholder decks.
///
/// Point count is a capability guard, not a quality metric: the normal shape checks below still
/// require aligned, finite, positive, strictly ascending data. Seven is deliberately the smallest
/// value that rejects every retired six-point placeholder while leaving table-density policy to
/// provenance review and golden tests.
pub const MIN_BC_CONVERSION_TABLE_POINTS: usize = 7;

/// The canonical output forms shared by native and WASM front ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcConversionFormat {
    Table,
    Csv,
    Json,
}

/// Exact conversion of one scalar BC at one Mach number.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScalarBcConversionV1 {
    pub schema_version: u32,
    pub source_drag_model: String,
    pub target_drag_model: String,
    /// Conventional small-arms BC in lb/in^2.
    pub source_bc: f64,
    /// Conventional small-arms BC in lb/in^2.
    pub target_bc: f64,
    pub mach: f64,
    /// Present when the velocity convenience API was used; canonical unit is ft/s.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub velocity_fps: Option<f64>,
    pub source_cd: f64,
    pub target_cd: f64,
    /// `target_cd / source_cd`.
    pub conversion_ratio: f64,
}

/// One input [`BCSegmentData`] band converted to a constant BC in another family.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BcSegmentConversionV1 {
    pub velocity_min_fps: f64,
    pub velocity_max_fps: f64,
    pub source_bc: f64,
    pub target_bc: f64,
    /// Width-normalized RMS relative retardation error inside this band.
    pub relative_rms: f64,
}

/// A velocity-banded conversion into one requested target family.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BandedBcConversionV1 {
    pub schema_version: u32,
    pub source_drag_model: String,
    pub target_drag_model: String,
    pub speed_of_sound_fps: f64,
    pub mach_min: f64,
    pub mach_max: f64,
    /// Same order and bounds as the caller's input segments.
    pub segments: Vec<BcSegmentConversionV1>,
    /// Width-normalized RMS relative retardation error across all converted bands.
    pub relative_rms: f64,
    pub integration_evaluations: usize,
}

impl BandedBcConversionV1 {
    /// Return the converted schedule in the engine's established velocity-band type.
    pub fn converted_segments(&self) -> Vec<BCSegmentData> {
        self.segments
            .iter()
            .map(|segment| BCSegmentData {
                velocity_min: segment.velocity_min_fps,
                velocity_max: segment.velocity_max_fps,
                bc_value: segment.target_bc,
            })
            .collect()
    }
}

/// One candidate family's best single scalar BC over a source band schedule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BcFamilyFitV1 {
    pub candidate_drag_model: String,
    pub fitted_bc: f64,
    /// Width-normalized RMS relative retardation error across the covered velocity range.
    pub relative_rms: f64,
    pub integration_evaluations: usize,
}

/// Ranking of candidate families by their one-scalar fit to a source band schedule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BcFamilyRecommendationV1 {
    pub schema_version: u32,
    pub source_drag_model: String,
    pub speed_of_sound_fps: f64,
    pub mach_min: f64,
    pub mach_max: f64,
    /// Ascending by residual, then by the caller's candidate order for exact ties.
    pub fits: Vec<BcFamilyFitV1>,
    pub recommended: BcFamilyFitV1,
}

/// One banded operation's requested conversion and family recommendation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BcSegmentAnalysisV1 {
    pub conversion: BandedBcConversionV1,
    pub recommendation: BcFamilyRecommendationV1,
}

/// Backward-compatible descriptive alias for [`BcSegmentAnalysisV1`].
pub type BcBandedAnalysisV1 = BcSegmentAnalysisV1;

/// Type-erased report accepted by the shared renderer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
pub enum BcConversionReportV1 {
    Scalar { result: ScalarBcConversionV1 },
    Banded { result: BcSegmentAnalysisV1 },
    Recommendation { result: BcFamilyRecommendationV1 },
}

impl From<BcSegmentAnalysisV1> for BcConversionReportV1 {
    fn from(result: BcSegmentAnalysisV1) -> Self {
        Self::Banded { result }
    }
}

/// Why a BC conversion or family fit was rejected.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum BcConversionError {
    #[error("{field} must be finite and greater than zero, got {value}")]
    InvalidPositiveValue { field: &'static str, value: f64 },
    #[error("Mach must be finite and non-negative, got {mach}")]
    InvalidMach { mach: f64 },
    #[error("no BC segments were supplied")]
    NoSegments,
    #[error(
        "BC segment #{index} has invalid bounds {velocity_min_fps}..{velocity_max_fps} ft/s; bounds must be finite, non-negative, and max must exceed min"
    )]
    InvalidSegmentBounds {
        index: usize,
        velocity_min_fps: f64,
        velocity_max_fps: f64,
    },
    #[error("BC segment #{index} must have a finite positive BC, got {bc}")]
    InvalidSegmentBc { index: usize, bc: f64 },
    #[error(
        "BC segments #{first_index} and #{second_index} overlap ({first_max_fps} > {second_min_fps} ft/s)"
    )]
    OverlappingSegments {
        first_index: usize,
        second_index: usize,
        first_max_fps: f64,
        second_min_fps: f64,
    },
    #[error("no candidate drag families were supplied")]
    NoCandidates,
    #[error("candidate drag family {drag_model} was supplied more than once")]
    DuplicateCandidate { drag_model: String },
    #[error("{drag_model} reference drag table is unusable: {reason}")]
    InvalidDragTable { drag_model: String, reason: String },
    #[error(
        "Mach {mach} is outside the {drag_model} reference table domain [{mach_min}, {mach_max}]; conversion does not clamp or extrapolate"
    )]
    MachOutsideTableDomain {
        drag_model: String,
        mach: f64,
        mach_min: f64,
        mach_max: f64,
    },
    #[error("BC conversion produced a non-finite result while computing {context}")]
    NonFiniteResult { context: &'static str },
    #[error("could not serialize BC conversion report: {message}")]
    Serialization { message: String },
}

/// Convert one scalar BC at a stated Mach number.
pub fn convert_bc_at_mach(
    source_bc: f64,
    source: DragModel,
    target: DragModel,
    mach: f64,
) -> Result<ScalarBcConversionV1, BcConversionError> {
    convert_bc_at_mach_impl(source_bc, source, target, mach, None)
}

/// Convert one scalar BC at a stated velocity and speed of sound, both in ft/s.
pub fn convert_bc_at_velocity(
    source_bc: f64,
    source: DragModel,
    target: DragModel,
    velocity_fps: f64,
    speed_of_sound_fps: f64,
) -> Result<ScalarBcConversionV1, BcConversionError> {
    validate_positive("velocity_fps", velocity_fps)?;
    validate_positive("speed_of_sound_fps", speed_of_sound_fps)?;
    convert_bc_at_mach_impl(
        source_bc,
        source,
        target,
        velocity_fps / speed_of_sound_fps,
        Some(velocity_fps),
    )
}

/// Convert every source band to a constant target-family BC over the same velocity bounds.
pub fn convert_bc_segments(
    segments: &[BCSegmentData],
    source: DragModel,
    target: DragModel,
    speed_of_sound_fps: f64,
) -> Result<BandedBcConversionV1, BcConversionError> {
    let validated = validate_segments_for_pair(segments, source, target, speed_of_sound_fps)?;
    let mut converted_by_input = vec![None; segments.len()];
    let mut total_squared_error = CompensatedSum::default();
    let mut total_weight = CompensatedSum::default();
    let mut integration_evaluations = 0usize;

    for segment in &validated.segments {
        let (target_bc, relative_rms, squared_error, weight, evaluations) = if source == target {
            // Preserve both the band BC and its identity residual exactly. The domain and
            // table were still validated above, so this does not bypass strictness.
            (
                segment.bc,
                0.0,
                0.0,
                segment.velocity_max_fps - segment.velocity_min_fps,
                0,
            )
        } else {
            let samples = collect_relative_samples(segment, source, target, speed_of_sound_fps)?;
            let solution = solve_relative_fit(&samples)?;
            (
                solution.fitted_bc,
                solution.relative_rms,
                solution.weighted_squared_error,
                solution.total_weight,
                solution.integration_evaluations,
            )
        };

        total_squared_error.add(squared_error);
        total_weight.add(weight);
        integration_evaluations += evaluations;
        converted_by_input[segment.original_index - 1] = Some(BcSegmentConversionV1 {
            velocity_min_fps: segment.velocity_min_fps,
            velocity_max_fps: segment.velocity_max_fps,
            source_bc: segment.bc,
            target_bc,
            relative_rms,
        });
    }

    let relative_rms = normalized_rms(total_squared_error.value(), total_weight.value())?;
    let converted_segments = converted_by_input
        .into_iter()
        .map(|segment| {
            segment.ok_or(BcConversionError::NonFiniteResult {
                context: "converted BC segment ordering",
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(BandedBcConversionV1 {
        schema_version: BC_CONVERSION_SCHEMA_VERSION_V1,
        source_drag_model: source.to_string(),
        target_drag_model: target.to_string(),
        speed_of_sound_fps,
        mach_min: validated.mach_min,
        mach_max: validated.mach_max,
        segments: converted_segments,
        relative_rms,
        integration_evaluations,
    })
}

/// Fit one scalar BC in `candidate` across all supplied source-family BC bands.
pub fn fit_bc_family(
    segments: &[BCSegmentData],
    source: DragModel,
    candidate: DragModel,
    speed_of_sound_fps: f64,
) -> Result<BcFamilyFitV1, BcConversionError> {
    let validated = validate_segments_for_pair(segments, source, candidate, speed_of_sound_fps)?;
    fit_validated_segments(&validated.segments, source, candidate, speed_of_sound_fps)
}

/// Rank candidate families by one-scalar relative-retardation least-squares residual.
pub fn recommend_bc_family(
    segments: &[BCSegmentData],
    source: DragModel,
    candidates: &[DragModel],
    speed_of_sound_fps: f64,
) -> Result<BcFamilyRecommendationV1, BcConversionError> {
    validate_positive("speed_of_sound_fps", speed_of_sound_fps)?;
    if candidates.is_empty() {
        return Err(BcConversionError::NoCandidates);
    }
    for (index, candidate) in candidates.iter().enumerate() {
        if candidates[..index].iter().any(|seen| seen == candidate) {
            return Err(BcConversionError::DuplicateCandidate {
                drag_model: candidate.to_string(),
            });
        }
    }

    // Validate source bands once up front to establish the report range. Each candidate fit
    // additionally validates its own table and pair domain before sampling.
    let source_validated =
        validate_segments_for_pair(segments, source, source, speed_of_sound_fps)?;
    let mut fits = Vec::with_capacity(candidates.len());
    for &candidate in candidates {
        fits.push(fit_bc_family(
            segments,
            source,
            candidate,
            speed_of_sound_fps,
        )?);
    }
    // `sort_by` is stable: an exact residual tie retains caller candidate order.
    fits.sort_by(|left, right| left.relative_rms.total_cmp(&right.relative_rms));
    let recommended = fits
        .first()
        .cloned()
        .ok_or(BcConversionError::NoCandidates)?;

    Ok(BcFamilyRecommendationV1 {
        schema_version: BC_CONVERSION_SCHEMA_VERSION_V1,
        source_drag_model: source.to_string(),
        speed_of_sound_fps,
        mach_min: source_validated.mach_min,
        mach_max: source_validated.mach_max,
        fits,
        recommended,
    })
}

/// Perform the requested band conversion and candidate-family recommendation in one call.
pub fn analyze_bc_segments(
    segments: &[BCSegmentData],
    source: DragModel,
    target: DragModel,
    candidates: &[DragModel],
    speed_of_sound_fps: f64,
) -> Result<BcSegmentAnalysisV1, BcConversionError> {
    Ok(BcSegmentAnalysisV1 {
        conversion: convert_bc_segments(segments, source, target, speed_of_sound_fps)?,
        recommendation: recommend_bc_family(segments, source, candidates, speed_of_sound_fps)?,
    })
}

/// Render a report in canonical ft/s units, identically for native and WASM callers.
pub fn format_bc_conversion_report(
    report: &BcConversionReportV1,
    format: BcConversionFormat,
) -> Result<String, BcConversionError> {
    match format {
        BcConversionFormat::Json => {
            let mut rendered = serde_json::to_string_pretty(report).map_err(|error| {
                BcConversionError::Serialization {
                    message: error.to_string(),
                }
            })?;
            rendered.push('\n');
            Ok(rendered)
        }
        BcConversionFormat::Table => Ok(format_report_table(report)),
        BcConversionFormat::Csv => Ok(format_report_csv(report)),
    }
}

fn convert_bc_at_mach_impl(
    source_bc: f64,
    source: DragModel,
    target: DragModel,
    mach: f64,
    velocity_fps: Option<f64>,
) -> Result<ScalarBcConversionV1, BcConversionError> {
    validate_positive("source_bc", source_bc)?;
    validate_mach(mach)?;
    validate_mach_for_model(source, mach)?;
    validate_mach_for_model(target, mach)?;

    let source_cd = get_drag_coefficient(mach, &source);
    let target_cd = if source == target {
        source_cd
    } else {
        get_drag_coefficient(mach, &target)
    };
    validate_computed_cd(source_cd, "source drag coefficient")?;
    validate_computed_cd(target_cd, "target drag coefficient")?;

    // Special-case identity so the public invariant is bit-exact rather than merely within
    // floating-point tolerance after multiplying by Cd/Cd.
    let (conversion_ratio, target_bc) = if source == target {
        (1.0, source_bc)
    } else {
        let ratio = target_cd / source_cd;
        let converted = source_bc * ratio;
        if !ratio.is_finite() || !converted.is_finite() || converted <= 0.0 {
            return Err(BcConversionError::NonFiniteResult {
                context: "scalar target BC",
            });
        }
        (ratio, converted)
    };

    Ok(ScalarBcConversionV1 {
        schema_version: BC_CONVERSION_SCHEMA_VERSION_V1,
        source_drag_model: source.to_string(),
        target_drag_model: target.to_string(),
        source_bc,
        target_bc,
        mach,
        velocity_fps,
        source_cd,
        target_cd,
        conversion_ratio,
    })
}

fn validate_positive(field: &'static str, value: f64) -> Result<(), BcConversionError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(BcConversionError::InvalidPositiveValue { field, value })
    }
}

fn validate_mach(mach: f64) -> Result<(), BcConversionError> {
    if mach.is_finite() && mach >= 0.0 {
        Ok(())
    } else {
        Err(BcConversionError::InvalidMach { mach })
    }
}

fn validate_computed_cd(cd: f64, context: &'static str) -> Result<(), BcConversionError> {
    if cd.is_finite() && cd > 0.0 {
        Ok(())
    } else {
        Err(BcConversionError::NonFiniteResult { context })
    }
}

fn validate_mach_for_model(drag_model: DragModel, mach: f64) -> Result<(), BcConversionError> {
    let (mach_min, mach_max) = validate_drag_table(drag_model)?;
    if mach < mach_min || mach > mach_max {
        return Err(BcConversionError::MachOutsideTableDomain {
            drag_model: drag_model.to_string(),
            mach,
            mach_min,
            mach_max,
        });
    }
    Ok(())
}

fn validate_drag_table(drag_model: DragModel) -> Result<(f64, f64), BcConversionError> {
    let table = reference_drag_table(&drag_model);
    validate_drag_table_shape(drag_model, table)?;
    Ok((
        table.mach_values[0],
        table.mach_values[table.mach_values.len() - 1],
    ))
}

fn validate_drag_table_shape(
    drag_model: DragModel,
    table: &DragTable,
) -> Result<(), BcConversionError> {
    let invalid = |reason: String| BcConversionError::InvalidDragTable {
        drag_model: drag_model.to_string(),
        reason,
    };
    if table.mach_values.len() != table.cd_values.len() {
        return Err(invalid(format!(
            "{} Mach values but {} Cd values",
            table.mach_values.len(),
            table.cd_values.len()
        )));
    }
    if table.mach_values.len() < MIN_BC_CONVERSION_TABLE_POINTS {
        return Err(invalid(format!(
            "needs at least {MIN_BC_CONVERSION_TABLE_POINTS} rows for conversion (legacy sparse/placeholder decks are refused), got {}",
            table.mach_values.len()
        )));
    }
    for (index, (&mach, &cd)) in table.mach_values.iter().zip(&table.cd_values).enumerate() {
        if !mach.is_finite() || mach < 0.0 {
            return Err(invalid(format!(
                "row {} Mach must be finite and non-negative, got {mach}",
                index + 1
            )));
        }
        if !cd.is_finite() || cd <= 0.0 {
            return Err(invalid(format!(
                "row {} Cd must be finite and positive, got {cd}",
                index + 1
            )));
        }
        if index > 0 && mach <= table.mach_values[index - 1] {
            return Err(invalid(format!(
                "Mach values must strictly ascend; row {} ({mach}) follows {}",
                index + 1,
                table.mach_values[index - 1]
            )));
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct ValidatedSegment {
    /// One-based index in the caller's slice, used in public errors and output restoration.
    original_index: usize,
    velocity_min_fps: f64,
    velocity_max_fps: f64,
    bc: f64,
}

struct ValidatedSegments {
    /// Ascending by velocity, regardless of caller order.
    segments: Vec<ValidatedSegment>,
    mach_min: f64,
    mach_max: f64,
}

fn validate_segments_for_pair(
    segments: &[BCSegmentData],
    source: DragModel,
    target: DragModel,
    speed_of_sound_fps: f64,
) -> Result<ValidatedSegments, BcConversionError> {
    validate_positive("speed_of_sound_fps", speed_of_sound_fps)?;
    if segments.is_empty() {
        return Err(BcConversionError::NoSegments);
    }

    let mut validated = Vec::with_capacity(segments.len());
    for (zero_based_index, segment) in segments.iter().enumerate() {
        let index = zero_based_index + 1;
        if !segment.velocity_min.is_finite()
            || !segment.velocity_max.is_finite()
            || segment.velocity_min < 0.0
            || segment.velocity_max <= segment.velocity_min
        {
            return Err(BcConversionError::InvalidSegmentBounds {
                index,
                velocity_min_fps: segment.velocity_min,
                velocity_max_fps: segment.velocity_max,
            });
        }
        if !segment.bc_value.is_finite() || segment.bc_value <= 0.0 {
            return Err(BcConversionError::InvalidSegmentBc {
                index,
                bc: segment.bc_value,
            });
        }
        validated.push(ValidatedSegment {
            original_index: index,
            velocity_min_fps: segment.velocity_min,
            velocity_max_fps: segment.velocity_max,
            bc: segment.bc_value,
        });
    }
    validated.sort_by(|left, right| {
        left.velocity_min_fps
            .total_cmp(&right.velocity_min_fps)
            .then_with(|| left.velocity_max_fps.total_cmp(&right.velocity_max_fps))
            .then_with(|| left.original_index.cmp(&right.original_index))
    });
    for pair in validated.windows(2) {
        let (left, right) = (pair[0], pair[1]);
        if right.velocity_min_fps < left.velocity_max_fps {
            return Err(BcConversionError::OverlappingSegments {
                first_index: left.original_index,
                second_index: right.original_index,
                first_max_fps: left.velocity_max_fps,
                second_min_fps: right.velocity_min_fps,
            });
        }
    }

    let source_domain = validate_drag_table(source)?;
    let target_domain = if source == target {
        source_domain
    } else {
        validate_drag_table(target)?
    };
    for segment in &validated {
        let low_mach = segment.velocity_min_fps / speed_of_sound_fps;
        let high_mach = segment.velocity_max_fps / speed_of_sound_fps;
        validate_mach_in_domain(source, low_mach, source_domain)?;
        validate_mach_in_domain(source, high_mach, source_domain)?;
        validate_mach_in_domain(target, low_mach, target_domain)?;
        validate_mach_in_domain(target, high_mach, target_domain)?;
    }

    let mach_min = validated[0].velocity_min_fps / speed_of_sound_fps;
    let mach_max = validated
        .iter()
        .map(|segment| segment.velocity_max_fps)
        .max_by(f64::total_cmp)
        .ok_or(BcConversionError::NoSegments)?
        / speed_of_sound_fps;

    Ok(ValidatedSegments {
        segments: validated,
        mach_min,
        mach_max,
    })
}

fn validate_mach_in_domain(
    drag_model: DragModel,
    mach: f64,
    (mach_min, mach_max): (f64, f64),
) -> Result<(), BcConversionError> {
    if mach < mach_min || mach > mach_max {
        Err(BcConversionError::MachOutsideTableDomain {
            drag_model: drag_model.to_string(),
            mach,
            mach_min,
            mach_max,
        })
    } else {
        Ok(())
    }
}

#[derive(Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        // Neumaier summation is deterministic and remains effective when the next term is
        // larger than the running sum (which can happen across very uneven velocity bands).
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

#[derive(Debug, Clone, Copy)]
struct WeightedRelativeSample {
    weight: f64,
    /// `candidate_cd * source_bc / source_cd`; candidate relative error is `q*z - 1`,
    /// where `q = 1 / candidate_bc` is the sole fitted parameter.
    z: f64,
}

struct RelativeSamples {
    values: Vec<WeightedRelativeSample>,
}

struct RelativeFitSolution {
    fitted_bc: f64,
    relative_rms: f64,
    weighted_squared_error: f64,
    total_weight: f64,
    integration_evaluations: usize,
}

// Eight-point Gauss-Legendre nodes and weights on [-1, 1]. The integration interval is
// additionally split at every knot from both live drag tables, so the sharp transonic shape is
// never averaged across a table discontinuity in interpolation slope.
const GAUSS_NODES_8: [f64; 8] = [
    -0.960_289_856_497_536_3,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_329,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const GAUSS_WEIGHTS_8: [f64; 8] = [
    0.101_228_536_290_376_3,
    0.222_381_034_453_374_5,
    0.313_706_645_877_887_3,
    0.362_683_783_378_362,
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

fn collect_relative_samples(
    segment: &ValidatedSegment,
    source: DragModel,
    candidate: DragModel,
    speed_of_sound_fps: f64,
) -> Result<RelativeSamples, BcConversionError> {
    let mut breakpoints = vec![segment.velocity_min_fps, segment.velocity_max_fps];
    append_table_breakpoints(
        &mut breakpoints,
        reference_drag_table(&source),
        speed_of_sound_fps,
        segment.velocity_min_fps,
        segment.velocity_max_fps,
    );
    if source != candidate {
        append_table_breakpoints(
            &mut breakpoints,
            reference_drag_table(&candidate),
            speed_of_sound_fps,
            segment.velocity_min_fps,
            segment.velocity_max_fps,
        );
    }
    breakpoints.sort_by(f64::total_cmp);
    breakpoints.dedup_by(|left, right| left.to_bits() == right.to_bits());

    let mut values = Vec::with_capacity((breakpoints.len().saturating_sub(1)) * 8);
    for interval in breakpoints.windows(2) {
        let low = interval[0];
        let high = interval[1];
        if high <= low {
            continue;
        }
        let midpoint = 0.5 * (low + high);
        let half_width = 0.5 * (high - low);
        for (&node, &base_weight) in GAUSS_NODES_8.iter().zip(&GAUSS_WEIGHTS_8) {
            let velocity_fps = midpoint + half_width * node;
            let mach = velocity_fps / speed_of_sound_fps;
            let source_cd = get_drag_coefficient(mach, &source);
            validate_computed_cd(source_cd, "band source drag coefficient")?;
            let z = if source == candidate {
                // Avoid Cd/Cd roundoff; this is also what makes a constant schedule's
                // same-family fit resolve to exactly its input BC at the samples.
                segment.bc
            } else {
                let candidate_cd = get_drag_coefficient(mach, &candidate);
                validate_computed_cd(candidate_cd, "band candidate drag coefficient")?;
                candidate_cd * segment.bc / source_cd
            };
            let weight = half_width * base_weight;
            if !z.is_finite() || z <= 0.0 || !weight.is_finite() || weight <= 0.0 {
                return Err(BcConversionError::NonFiniteResult {
                    context: "band integration sample",
                });
            }
            values.push(WeightedRelativeSample { weight, z });
        }
    }
    if values.is_empty() {
        return Err(BcConversionError::NonFiniteResult {
            context: "band integration grid",
        });
    }
    Ok(RelativeSamples { values })
}

fn append_table_breakpoints(
    breakpoints: &mut Vec<f64>,
    table: &DragTable,
    speed_of_sound_fps: f64,
    velocity_min_fps: f64,
    velocity_max_fps: f64,
) {
    for &mach in &table.mach_values {
        let velocity_fps = mach * speed_of_sound_fps;
        if velocity_fps > velocity_min_fps && velocity_fps < velocity_max_fps {
            breakpoints.push(velocity_fps);
        }
    }
}

fn solve_relative_fit(samples: &RelativeSamples) -> Result<RelativeFitSolution, BcConversionError> {
    let mut weighted_z = CompensatedSum::default();
    let mut weighted_z_squared = CompensatedSum::default();
    let mut total_weight = CompensatedSum::default();
    for sample in &samples.values {
        weighted_z.add(sample.weight * sample.z);
        weighted_z_squared.add(sample.weight * sample.z * sample.z);
        total_weight.add(sample.weight);
    }
    let numerator = weighted_z.value();
    let denominator = weighted_z_squared.value();
    let total_weight = total_weight.value();
    if !numerator.is_finite()
        || !denominator.is_finite()
        || !total_weight.is_finite()
        || numerator <= 0.0
        || denominator <= 0.0
        || total_weight <= 0.0
    {
        return Err(BcConversionError::NonFiniteResult {
            context: "least-squares moments",
        });
    }

    let reciprocal_bc = numerator / denominator;
    let fitted_bc = 1.0 / reciprocal_bc;
    if !reciprocal_bc.is_finite()
        || reciprocal_bc <= 0.0
        || !fitted_bc.is_finite()
        || fitted_bc <= 0.0
    {
        return Err(BcConversionError::NonFiniteResult {
            context: "least-squares fitted BC",
        });
    }

    // Compute the residual directly at the same quadrature nodes rather than subtracting two
    // nearly equal moment products. This keeps exact-family fits near machine epsilon instead of
    // magnifying cancellation into a misleading ~1e-8 RMS score.
    let mut squared_error = CompensatedSum::default();
    for sample in &samples.values {
        let relative_error = reciprocal_bc * sample.z - 1.0;
        squared_error.add(sample.weight * relative_error * relative_error);
    }
    let weighted_squared_error = squared_error.value();
    let relative_rms = normalized_rms(weighted_squared_error, total_weight)?;

    Ok(RelativeFitSolution {
        fitted_bc,
        relative_rms,
        weighted_squared_error,
        total_weight,
        integration_evaluations: samples.values.len(),
    })
}

fn normalized_rms(squared_error: f64, weight: f64) -> Result<f64, BcConversionError> {
    if !squared_error.is_finite() || squared_error < 0.0 || !weight.is_finite() || weight <= 0.0 {
        return Err(BcConversionError::NonFiniteResult {
            context: "normalized relative RMS residual",
        });
    }
    let rms = (squared_error / weight).sqrt();
    if rms.is_finite() {
        Ok(rms)
    } else {
        Err(BcConversionError::NonFiniteResult {
            context: "normalized relative RMS residual",
        })
    }
}

fn fit_validated_segments(
    segments: &[ValidatedSegment],
    source: DragModel,
    candidate: DragModel,
    speed_of_sound_fps: f64,
) -> Result<BcFamilyFitV1, BcConversionError> {
    let mut all_samples = RelativeSamples { values: Vec::new() };
    for segment in segments {
        all_samples.values.extend(
            collect_relative_samples(segment, source, candidate, speed_of_sound_fps)?.values,
        );
    }
    let solution = solve_relative_fit(&all_samples)?;
    Ok(BcFamilyFitV1 {
        candidate_drag_model: candidate.to_string(),
        fitted_bc: solution.fitted_bc,
        relative_rms: solution.relative_rms,
        integration_evaluations: solution.integration_evaluations,
    })
}

fn format_report_table(report: &BcConversionReportV1) -> String {
    let mut output = String::new();
    match report {
        BcConversionReportV1::Scalar { result } => {
            let _ = writeln!(output, "BC drag-family conversion");
            let _ = writeln!(
                output,
                "source  target  mach      velocity_fps  source_bc  source_cd  target_cd  ratio      target_bc"
            );
            let velocity = result
                .velocity_fps
                .map(|value| format!("{value:.3}"))
                .unwrap_or_else(|| "-".to_string());
            let _ = writeln!(
                output,
                "{:<6}  {:<6}  {:>8.5}  {:>12}  {:>9.6}  {:>9.6}  {:>9.6}  {:>9.6}  {:>9.6}",
                result.source_drag_model,
                result.target_drag_model,
                result.mach,
                velocity,
                result.source_bc,
                result.source_cd,
                result.target_cd,
                result.conversion_ratio,
                result.target_bc,
            );
        }
        BcConversionReportV1::Banded { result } => {
            let conversion = &result.conversion;
            let recommendation = &result.recommendation;
            let _ = writeln!(
                output,
                "BC band conversion {} -> {} (speed_of_sound_fps={:.3}, Mach {:.5}..{:.5})",
                conversion.source_drag_model,
                conversion.target_drag_model,
                conversion.speed_of_sound_fps,
                conversion.mach_min,
                conversion.mach_max,
            );
            let _ = writeln!(
                output,
                "velocity_min_fps  velocity_max_fps  source_bc  target_bc  relative_rms"
            );
            for segment in &conversion.segments {
                let _ = writeln!(
                    output,
                    "{:>16.3}  {:>16.3}  {:>9.6}  {:>9.6}  {:>12.8}",
                    segment.velocity_min_fps,
                    segment.velocity_max_fps,
                    segment.source_bc,
                    segment.target_bc,
                    segment.relative_rms,
                );
            }
            let _ = writeln!(
                output,
                "combined_relative_rms={:.8}",
                conversion.relative_rms
            );
            append_recommendation_table(&mut output, recommendation);
        }
        BcConversionReportV1::Recommendation { result } => {
            append_recommendation_table(&mut output, result);
        }
    }
    output
}

fn append_recommendation_table(output: &mut String, recommendation: &BcFamilyRecommendationV1) {
    let _ = writeln!(
        output,
        "Family recommendation (source={}, speed_of_sound_fps={:.3}, Mach {:.5}..{:.5})",
        recommendation.source_drag_model,
        recommendation.speed_of_sound_fps,
        recommendation.mach_min,
        recommendation.mach_max,
    );
    let _ = writeln!(output, "candidate  fitted_bc  relative_rms  recommended");
    for fit in &recommendation.fits {
        let selected = fit.candidate_drag_model == recommendation.recommended.candidate_drag_model
            && fit.fitted_bc.to_bits() == recommendation.recommended.fitted_bc.to_bits()
            && fit.relative_rms.to_bits() == recommendation.recommended.relative_rms.to_bits();
        let _ = writeln!(
            output,
            "{:<9}  {:>9.6}  {:>12.8}  {}",
            fit.candidate_drag_model,
            fit.fitted_bc,
            fit.relative_rms,
            if selected { "yes" } else { "no" },
        );
    }
}

fn format_report_csv(report: &BcConversionReportV1) -> String {
    let mut output = String::new();
    match report {
        BcConversionReportV1::Scalar { result } => {
            let _ = writeln!(
                output,
                "schema_version,source_drag_model,target_drag_model,mach,velocity_fps,source_bc,source_cd,target_cd,conversion_ratio,target_bc"
            );
            let velocity = result
                .velocity_fps
                .map(|value| format!("{value:.9}"))
                .unwrap_or_default();
            let _ = writeln!(
                output,
                "{},{},{},{:.9},{},{:.9},{:.9},{:.9},{:.9},{:.9}",
                result.schema_version,
                result.source_drag_model,
                result.target_drag_model,
                result.mach,
                velocity,
                result.source_bc,
                result.source_cd,
                result.target_cd,
                result.conversion_ratio,
                result.target_bc,
            );
        }
        BcConversionReportV1::Banded { result } => {
            let conversion = &result.conversion;
            let recommendation = &result.recommendation;
            let _ = writeln!(
                output,
                "record_type,schema_version,source_drag_model,target_drag_model,speed_of_sound_fps,mach_min,mach_max,velocity_min_fps,velocity_max_fps,source_bc,target_bc,candidate_drag_model,fitted_bc,relative_rms,recommended,integration_evaluations"
            );
            for segment in &conversion.segments {
                let _ = writeln!(
                    output,
                    "segment,{},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},,,{:.12},,",
                    conversion.schema_version,
                    conversion.source_drag_model,
                    conversion.target_drag_model,
                    conversion.speed_of_sound_fps,
                    conversion.mach_min,
                    conversion.mach_max,
                    segment.velocity_min_fps,
                    segment.velocity_max_fps,
                    segment.source_bc,
                    segment.target_bc,
                    segment.relative_rms,
                );
            }
            let summary = [
                "conversion_summary".to_string(),
                conversion.schema_version.to_string(),
                conversion.source_drag_model.clone(),
                conversion.target_drag_model.clone(),
                format!("{:.9}", conversion.speed_of_sound_fps),
                format!("{:.9}", conversion.mach_min),
                format!("{:.9}", conversion.mach_max),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                format!("{:.12}", conversion.relative_rms),
                String::new(),
                conversion.integration_evaluations.to_string(),
            ];
            let _ = writeln!(output, "{}", summary.join(","));
            for fit in &recommendation.fits {
                let selected = fit.candidate_drag_model
                    == recommendation.recommended.candidate_drag_model
                    && fit.fitted_bc.to_bits() == recommendation.recommended.fitted_bc.to_bits()
                    && fit.relative_rms.to_bits()
                        == recommendation.recommended.relative_rms.to_bits();
                let _ = writeln!(
                    output,
                    "family_fit,{},{},,{:.9},{:.9},{:.9},,,,,{},{:.9},{:.12},{},{}",
                    recommendation.schema_version,
                    recommendation.source_drag_model,
                    recommendation.speed_of_sound_fps,
                    recommendation.mach_min,
                    recommendation.mach_max,
                    fit.candidate_drag_model,
                    fit.fitted_bc,
                    fit.relative_rms,
                    selected,
                    fit.integration_evaluations,
                );
            }
        }
        BcConversionReportV1::Recommendation { result } => {
            let _ = writeln!(
                output,
                "schema_version,source_drag_model,speed_of_sound_fps,mach_min,mach_max,candidate_drag_model,fitted_bc,relative_rms,recommended,integration_evaluations"
            );
            for fit in &result.fits {
                let selected = fit.candidate_drag_model == result.recommended.candidate_drag_model
                    && fit.fitted_bc.to_bits() == result.recommended.fitted_bc.to_bits()
                    && fit.relative_rms.to_bits() == result.recommended.relative_rms.to_bits();
                let _ = writeln!(
                    output,
                    "{},{},{:.9},{:.9},{:.9},{},{:.9},{:.12},{},{}",
                    result.schema_version,
                    result.source_drag_model,
                    result.speed_of_sound_fps,
                    result.mach_min,
                    result.mach_max,
                    fit.candidate_drag_model,
                    fit.fitted_bc,
                    fit.relative_rms,
                    selected,
                    fit.integration_evaluations,
                );
            }
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    const SOS_FPS: f64 = 1_100.0;

    fn segment(min: f64, max: f64, bc: f64) -> BCSegmentData {
        BCSegmentData {
            velocity_min: min,
            velocity_max: max,
            bc_value: bc,
        }
    }

    fn constant_schedule(bc: f64) -> Vec<BCSegmentData> {
        vec![
            segment(660.0, 1_100.0, bc),
            segment(1_100.0, 1_650.0, bc),
            segment(1_650.0, 2_200.0, bc),
            segment(2_200.0, 2_750.0, bc),
            segment(2_750.0, 3_300.0, bc),
        ]
    }

    #[test]
    fn scalar_conversion_preserves_reference_retardation() {
        let converted =
            convert_bc_at_mach(0.500, DragModel::G1, DragModel::G7, 2.0).expect("convert");
        let source_retardation = converted.source_cd / converted.source_bc;
        let target_retardation = converted.target_cd / converted.target_bc;
        assert!(
            (source_retardation - target_retardation).abs() <= 2.0 * f64::EPSILON,
            "{source_retardation} != {target_retardation}"
        );
        assert!((converted.target_bc - 0.251_095_382_541).abs() < 1e-12);
    }

    #[test]
    fn jbm_g1_0492_at_3000_fps_converts_to_g7_0242() {
        // Widely documented JBM conversion example under standard-atmosphere speed of sound.
        // Pinning the unrounded engine value catches both a ratio inversion and table drift;
        // the public-facing three-decimal answer is G7 0.242.
        let converted = convert_bc_at_velocity(
            0.492,
            DragModel::G1,
            DragModel::G7,
            3_000.0,
            crate::constants::SPEED_OF_SOUND_MPS / crate::constants::FPS_TO_MPS,
        )
        .expect("JBM example");
        assert!(
            (converted.target_bc - 0.242_195_467_3).abs() < 5e-10,
            "got {} at Mach {}",
            converted.target_bc,
            converted.mach
        );
        assert_eq!((converted.target_bc * 1_000.0).round() / 1_000.0, 0.242);
    }

    #[test]
    fn scalar_round_trip_recovers_input_at_the_same_mach() {
        for mach in [0.5, 0.95, 1.0, 1.25, 2.0, 3.0, 5.0] {
            let forward =
                convert_bc_at_mach(0.537, DragModel::G1, DragModel::G7, mach).expect("G1 to G7");
            let backward =
                convert_bc_at_mach(forward.target_bc, DragModel::G7, DragModel::G1, mach)
                    .expect("G7 to G1");
            assert!(
                (backward.target_bc - 0.537).abs() < 2e-15,
                "Mach {mach}: {}",
                backward.target_bc
            );
        }
    }

    #[test]
    fn conversion_ratio_is_mach_dependent_not_the_cluster_shortcut() {
        let transonic =
            convert_bc_at_mach(1.0, DragModel::G1, DragModel::G7, 1.0).expect("transonic");
        let supersonic =
            convert_bc_at_mach(1.0, DragModel::G1, DragModel::G7, 2.0).expect("supersonic");
        assert!((transonic.conversion_ratio - supersonic.conversion_ratio).abs() > 0.25);
        assert!((transonic.conversion_ratio - 1.0 / 1.98).abs() > 0.25);
    }

    #[test]
    fn scalar_and_banded_identity_preserve_bc_bits() {
        let source_bc = f64::from_bits(0x3f_df_7c_ed_91_68_72_b0);
        let scalar = convert_bc_at_mach(source_bc, DragModel::G7, DragModel::G7, 2.1)
            .expect("scalar identity");
        assert_eq!(scalar.target_bc.to_bits(), source_bc.to_bits());
        assert_eq!(scalar.conversion_ratio.to_bits(), 1.0f64.to_bits());

        let bands = vec![
            segment(2_000.0, 2_400.0, source_bc),
            segment(1_500.0, 2_000.0, 0.211_111_111_111_111_1),
        ];
        let converted = convert_bc_segments(&bands, DragModel::G7, DragModel::G7, SOS_FPS)
            .expect("band identity");
        for (actual, expected) in converted.segments.iter().zip(&bands) {
            assert_eq!(actual.target_bc.to_bits(), expected.bc_value.to_bits());
            assert_eq!(actual.relative_rms.to_bits(), 0.0f64.to_bits());
        }
        assert_eq!(converted.relative_rms.to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn conversion_rejects_endpoint_clamping_and_invalid_values() {
        assert!(matches!(
            convert_bc_at_mach(0.5, DragModel::G1, DragModel::G7, 5.000_001),
            Err(BcConversionError::MachOutsideTableDomain { .. })
        ));
        assert!(matches!(
            convert_bc_at_mach(0.5, DragModel::G1, DragModel::G7, -0.1),
            Err(BcConversionError::InvalidMach { .. })
        ));
        assert!(matches!(
            convert_bc_at_mach(f64::NAN, DragModel::G1, DragModel::G7, 2.0),
            Err(BcConversionError::InvalidPositiveValue {
                field: "source_bc",
                ..
            })
        ));
        assert!(matches!(
            convert_bc_at_velocity(0.5, DragModel::G1, DragModel::G7, 2_500.0, f64::INFINITY),
            Err(BcConversionError::InvalidPositiveValue {
                field: "speed_of_sound_fps",
                ..
            })
        ));
        assert!(matches!(
            convert_bc_segments(
                &[segment(0.0, f64::INFINITY, 0.5)],
                DragModel::G1,
                DragModel::G7,
                SOS_FPS,
            ),
            Err(BcConversionError::InvalidSegmentBounds { .. })
        ));
    }

    #[test]
    fn sparse_and_placeholder_tables_are_refused() {
        for points in [2usize, 6] {
            let table = DragTable::new(
                (0..points).map(|index| index as f64).collect(),
                vec![0.2; points],
            );
            let error = validate_drag_table_shape(DragModel::G1, &table)
                .expect_err("sparse table must be refused");
            assert!(matches!(error, BcConversionError::InvalidDragTable { .. }));
            assert!(error.to_string().contains("placeholder"));
        }

        let minimum_valid = DragTable::new(
            (0..MIN_BC_CONVERSION_TABLE_POINTS)
                .map(|index| index as f64)
                .collect(),
            vec![0.2; MIN_BC_CONVERSION_TABLE_POINTS],
        );
        validate_drag_table_shape(DragModel::G1, &minimum_valid).expect("seven rows are enough");
    }

    #[test]
    fn unsorted_bands_preserve_output_order_and_overlaps_are_rejected() {
        let unsorted = vec![
            segment(2_200.0, 2_800.0, 0.51),
            segment(700.0, 1_200.0, 0.43),
            segment(1_400.0, 2_000.0, 0.47),
        ];
        let converted = convert_bc_segments(&unsorted, DragModel::G1, DragModel::G7, SOS_FPS)
            .expect("unsorted input");
        for (actual, expected) in converted.segments.iter().zip(&unsorted) {
            assert_eq!(
                actual.velocity_min_fps.to_bits(),
                expected.velocity_min.to_bits()
            );
            assert_eq!(
                actual.velocity_max_fps.to_bits(),
                expected.velocity_max.to_bits()
            );
            assert_eq!(actual.source_bc.to_bits(), expected.bc_value.to_bits());
        }

        let overlap = vec![
            segment(1_000.0, 2_000.0, 0.5),
            segment(1_900.0, 2_500.0, 0.5),
        ];
        assert!(matches!(
            convert_bc_segments(&overlap, DragModel::G1, DragModel::G7, SOS_FPS),
            Err(BcConversionError::OverlappingSegments { .. })
        ));
    }

    #[test]
    fn band_split_and_permutation_leave_fit_nearly_unchanged() {
        let one = vec![segment(660.0, 3_300.0, 0.5)];
        let split = vec![
            segment(660.0, 1_430.0, 0.5),
            segment(1_430.0, 2_310.0, 0.5),
            segment(2_310.0, 3_300.0, 0.5),
        ];
        let permuted = vec![split[2].clone(), split[0].clone(), split[1].clone()];
        let fit_one = fit_bc_family(&one, DragModel::G1, DragModel::G7, SOS_FPS).expect("one");
        let fit_split =
            fit_bc_family(&split, DragModel::G1, DragModel::G7, SOS_FPS).expect("split");
        let fit_permuted =
            fit_bc_family(&permuted, DragModel::G1, DragModel::G7, SOS_FPS).expect("permuted");
        assert!((fit_one.fitted_bc - fit_split.fitted_bc).abs() < 1e-11);
        assert!((fit_one.relative_rms - fit_split.relative_rms).abs() < 1e-11);
        assert_eq!(
            fit_split.fitted_bc.to_bits(),
            fit_permuted.fitted_bc.to_bits()
        );
        assert_eq!(
            fit_split.relative_rms.to_bits(),
            fit_permuted.relative_rms.to_bits()
        );
    }

    #[test]
    fn constant_g1_schedule_recommends_g1() {
        let recommendation = recommend_bc_family(
            &constant_schedule(0.5),
            DragModel::G1,
            &[DragModel::G1, DragModel::G7],
            SOS_FPS,
        )
        .expect("recommend");
        assert_eq!(recommendation.recommended.candidate_drag_model, "G1");
        assert!((recommendation.recommended.fitted_bc - 0.5).abs() < 1e-15);
        assert!(recommendation.recommended.relative_rms < 1e-15);
        assert!(recommendation.fits[1].relative_rms > 0.05);
    }

    #[test]
    fn g1_bands_synthesized_from_constant_g7_recommend_g7() {
        let g7_bc = 0.25;
        let mut g1_bands = Vec::new();
        let mut low_mach = 0.6;
        while low_mach < 3.0 {
            let high_mach = low_mach + 0.1;
            let midpoint = 0.5 * (low_mach + high_mach);
            let as_g1 = convert_bc_at_mach(g7_bc, DragModel::G7, DragModel::G1, midpoint)
                .expect("synthetic G1 point");
            g1_bands.push(segment(
                low_mach * SOS_FPS,
                high_mach * SOS_FPS,
                as_g1.target_bc,
            ));
            low_mach = high_mach;
        }
        let recommendation = recommend_bc_family(
            &g1_bands,
            DragModel::G1,
            &[DragModel::G1, DragModel::G7],
            SOS_FPS,
        )
        .expect("recommend");
        assert_eq!(recommendation.recommended.candidate_drag_model, "G7");
        assert!(
            (recommendation.recommended.fitted_bc - g7_bc).abs() < 0.003,
            "fitted {}",
            recommendation.recommended.fitted_bc
        );
        assert!(recommendation.recommended.relative_rms < recommendation.fits[1].relative_rms);
    }

    #[test]
    fn recommendation_is_scale_invariant() {
        let original = constant_schedule(0.5);
        let scaled: Vec<_> = original
            .iter()
            .map(|band| segment(band.velocity_min, band.velocity_max, band.bc_value * 2.5))
            .collect();
        let candidates = [DragModel::G1, DragModel::G7];
        let first =
            recommend_bc_family(&original, DragModel::G1, &candidates, SOS_FPS).expect("original");
        let second =
            recommend_bc_family(&scaled, DragModel::G1, &candidates, SOS_FPS).expect("scaled");
        assert_eq!(
            first.recommended.candidate_drag_model,
            second.recommended.candidate_drag_model
        );
        for (left, right) in first.fits.iter().zip(&second.fits) {
            assert_eq!(left.candidate_drag_model, right.candidate_drag_model);
            assert!((right.fitted_bc / left.fitted_bc - 2.5).abs() < 2e-14);
            assert!((left.relative_rms - right.relative_rms).abs() < 2e-15);
        }
    }

    #[test]
    fn combined_formatters_are_complete_newline_terminated_documents() {
        let analysis = analyze_bc_segments(
            &constant_schedule(0.5),
            DragModel::G1,
            DragModel::G7,
            &[DragModel::G1, DragModel::G7],
            SOS_FPS,
        )
        .expect("analysis");
        let report = BcConversionReportV1::Banded { result: analysis };

        let json = format_bc_conversion_report(&report, BcConversionFormat::Json).expect("json");
        assert!(json.starts_with('{') && json.ends_with("\n"));
        let document: serde_json::Value = serde_json::from_str(&json).expect("pure JSON");
        assert_eq!(document["mode"], "banded");
        assert_eq!(
            document["result"]["conversion"]["schema_version"],
            BC_CONVERSION_SCHEMA_VERSION_V1
        );
        assert_eq!(
            document["result"]["recommendation"]["schema_version"],
            BC_CONVERSION_SCHEMA_VERSION_V1
        );

        let table = format_bc_conversion_report(&report, BcConversionFormat::Table).expect("table");
        assert!(table.ends_with('\n'));
        assert!(table.contains("velocity_min_fps"));
        assert!(table.contains("Family recommendation"));

        let csv = format_bc_conversion_report(&report, BcConversionFormat::Csv).expect("csv");
        assert!(csv.ends_with('\n'));
        assert!(csv.starts_with("record_type,schema_version,source_drag_model"));
        assert!(csv.contains("segment,1,G1,G7"));
        assert!(csv.contains("conversion_summary,1,G1,G7"));
        assert!(csv.contains("family_fit,1,G1"));
        for row in csv.lines() {
            assert_eq!(
                row.split(',').count(),
                16,
                "banded CSV row does not match its header: {row}"
            );
        }
        let summary = csv
            .lines()
            .find(|row| row.starts_with("conversion_summary,"))
            .expect("conversion summary row");
        let summary_fields: Vec<_> = summary.split(',').collect();
        let BcConversionReportV1::Banded { result } = &report else {
            unreachable!("test report is banded")
        };
        assert_eq!(
            summary_fields[13],
            format!("{:.12}", result.conversion.relative_rms)
        );
        assert!(!summary_fields[15].is_empty());
    }

    #[test]
    fn candidate_validation_rejects_empty_and_duplicate_sets() {
        let schedule = constant_schedule(0.5);
        assert!(matches!(
            recommend_bc_family(&schedule, DragModel::G1, &[], SOS_FPS),
            Err(BcConversionError::NoCandidates)
        ));
        assert!(matches!(
            recommend_bc_family(
                &schedule,
                DragModel::G1,
                &[DragModel::G7, DragModel::G7],
                SOS_FPS,
            ),
            Err(BcConversionError::DuplicateCandidate { .. })
        ));
    }
}
