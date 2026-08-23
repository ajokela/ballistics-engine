//! Observation-range experiment design for joint muzzle-velocity/BC truing (MBA-1346).
//!
//! The planner evaluates the same trajectory forward model and finite-difference
//! Jacobian as the truing fitter.  It therefore recommends *where to measure*,
//! without inventing observations or mutating the supplied profile.  Candidate
//! ranges that are malformed, duplicated, or unreachable are retained as
//! diagnostics rather than silently disappearing.

use std::error::Error;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::truing::{
    truing_jacobian_rows, DropUnit, TruingJacobianRow, TruingModelInputsV1,
    TRUING_MAX_CONDITION_NUMBER, TRUING_MIN_BC_SENSITIVITY_RATIO,
};

/// Maximum raw combination count for exact exhaustive design search.
///
/// Larger design spaces use a deterministic greedy construction followed by
/// deterministic one-for-one exchanges.  The trajectory/Jacobian evaluation is
/// cached per candidate, so neither strategy re-runs the forward model while
/// comparing station sets.
pub const TRUING_PLAN_EXHAUSTIVE_COMBINATION_LIMIT_V1: u64 = 100_000;

const DUPLICATE_RANGE_TOLERANCE_YD: f64 = 1.0e-6;
const SEPARATION_TOLERANCE_YD: f64 = 1.0e-9;
const INFORMATION_EIGENVALUE_FLOOR: f64 = 1.0e-12;

/// Input to the v1 range-design planner.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TruingExperimentPlanRequestV1 {
    /// Nominal scalar-BC load and atmospheric model.
    pub model: TruingModelInputsV1,
    /// Discrete ranges available to the shooter, in yards.
    ///
    /// Use [`discretize_truing_range_interval_v1`] when a facility is described
    /// as an interval rather than as discrete target stations.
    pub candidate_ranges_yd: Vec<f64>,
    /// Exact number of stations the returned design must contain.
    pub observation_count: usize,
    /// Smallest permitted distance between any two selected stations, in yards.
    pub minimum_separation_yd: f64,
    /// One-standard-deviation measurement resolution, expressed in `drop_unit`.
    pub measurement_sigma_1sd: f64,
    /// Unit in which drop is measured and `measurement_sigma_1sd` is expressed.
    pub drop_unit: DropUnit,
}

/// Whether the selected station set can locally separate MV and BC at the
/// supplied nominal profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingPlanModeV1 {
    JointMvBc,
    MvOnly,
}

impl fmt::Display for TruingPlanModeV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::JointMvBc => "joint_mv_bc",
            Self::MvOnly => "mv_only",
        })
    }
}

/// Search strategy used for the discrete station-set optimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingPlanSearchStrategyV1 {
    Exhaustive,
    GreedyExchange,
}

/// Why an input candidate could not participate in the experiment design.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingCandidateRejectionReasonV1 {
    InvalidRange,
    DuplicateRange,
    Unreachable,
}

/// A supplied candidate that was excluded before station-set optimization.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RejectedTruingCandidateV1 {
    /// Zero-based position in `candidate_ranges_yd`.
    pub input_index: usize,
    /// Original value supplied by the caller.  It may be NaN or infinite, so
    /// JSON renderers should sanitize it when `reason == invalid_range`.
    pub range_yd: f64,
    pub reason: TruingCandidateRejectionReasonV1,
    pub detail: String,
}

/// Local sensitivity and information contribution for a selected target station.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TruingPlanStationV1 {
    /// Zero-based position in the caller's candidate list.
    pub input_index: usize,
    pub range_yd: f64,
    /// Nominal predicted drop, expressed in the plan's drop unit.
    pub predicted_drop: f64,
    /// Change in predicted drop, measured in observation sigmas, for a unit
    /// fractional change in MV (`MV * d(drop)/d(MV) / sigma`).
    pub scaled_mv_sensitivity: f64,
    /// Change in predicted drop, measured in observation sigmas, for a unit
    /// fractional change in BC (`BC * d(drop)/d(BC) / sigma`).
    pub scaled_bc_sensitivity: f64,
    /// Loss of local Gaussian information when this station is removed, in
    /// nats: `0.5 * (log det(I + F_all) - log det(I + F_without))`.  `I` is an
    /// explicit identity reference information matrix in fractional MV/BC
    /// coordinates; this is a design score, not a fitted posterior claim.
    pub leave_one_out_information_gain_nats: f64,
}

/// Information diagnostics for the selected design.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TruingPlanInformationV1 {
    /// `||BC*d(drop)/d(BC)|| / ||MV*d(drop)/d(MV)||`, after weighting by
    /// the declared observation sigma.
    pub sensitivity_ratio: f64,
    /// Condition number of the independently column-normalized MV/BC normal
    /// matrix. `None` means one column is zero or the design is rank deficient.
    pub condition_number: Option<f64>,
    /// Smallest singular value of the observation-sigma-weighted fractional
    /// Jacobian.  Its unit is inverse observation sigma (dimensionless here).
    pub minimum_singular_value: f64,
    pub maximum_singular_value: f64,
    /// Local 1-sigma uncertainty along the weakest MV/BC fractional-parameter
    /// direction (`1 / minimum_singular_value`).  It scales linearly with the
    /// declared measurement sigma. `None` denotes a rank-deficient design.
    pub weak_axis_fractional_sigma_1sd: Option<f64>,
    /// Natural log of the unregularized 2x2 information determinant. `None`
    /// denotes a rank-deficient design.
    pub log_determinant: Option<f64>,
    /// `0.5 * log det(I + F)`, a finite local-Gaussian information summary in
    /// nats.  `I` is an identity reference information matrix in fractional
    /// MV/BC coordinates (equivalently, unit reference covariance for those
    /// fractional perturbations), not a hidden truing prior.
    pub expected_information_gain_nats: f64,
}

/// Successful v1 experiment design.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TruingExperimentPlanV1 {
    pub mode: TruingPlanModeV1,
    pub selected_stations: Vec<TruingPlanStationV1>,
    /// Reachable, unique candidates not selected by the optimizer.
    pub unselected_candidate_ranges_yd: Vec<f64>,
    /// Invalid, duplicated, and forward-model-unreachable input candidates.
    pub rejected_candidates: Vec<RejectedTruingCandidateV1>,
    pub information: TruingPlanInformationV1,
    pub requested_observation_count: usize,
    pub minimum_separation_yd: f64,
    pub measurement_sigma_1sd: f64,
    pub measurement_drop_unit: String,
    pub eligible_candidate_count: usize,
    pub raw_combination_count: u64,
    pub evaluated_design_count: u64,
    pub search_strategy: TruingPlanSearchStrategyV1,
    /// Explicit limitations or actions needed before attempting a joint fit.
    pub warnings: Vec<String>,
}

/// Stable categories for planner failures that occur before a design can be returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingPlanErrorCodeV1 {
    InvalidRequest,
    InsufficientReachableCandidates,
    NoFeasibleDesign,
}

/// Structured planner error.  Candidate diagnostics are retained when the
/// request was valid but too few stations survived preprocessing.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TruingPlanErrorV1 {
    pub code: TruingPlanErrorCodeV1,
    pub message: String,
    pub rejected_candidates: Vec<RejectedTruingCandidateV1>,
}

impl fmt::Display for TruingPlanErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl Error for TruingPlanErrorV1 {}

impl TruingPlanErrorV1 {
    fn invalid(message: impl Into<String>) -> Self {
        Self {
            code: TruingPlanErrorCodeV1::InvalidRequest,
            message: message.into(),
            rejected_candidates: Vec::new(),
        }
    }

    /// Structured detail payload for a JSON bridge error envelope, mirroring
    /// [`crate::truing_service::DsfServiceErrorV1::failure_details`]: an INHERENT method
    /// (not a trait impl — a blanket trait impl collides with specific ones under E0119,
    /// and this type must keep compiling with the `bridge` feature off). Surfaces the
    /// machine-readable [`TruingPlanErrorCodeV1`] under the same `reason` key its siblings
    /// (`DsfServiceErrorV1::failure_details`, `OpticError::failure_details`) use as their
    /// top-level discriminator — including the `rejected_candidates` entries nested in this
    /// same payload, which already use `reason` per item — so a caller can branch on
    /// `error.details.reason` uniformly across the whole `true.*` family instead of
    /// special-casing `true.plan`. Also carries the rejected-candidate diagnostics
    /// (`InsufficientReachableCandidates`/`NoFeasibleDesign` both retain them) a wizard
    /// would branch on to suggest a different candidate range or observation count.
    pub fn failure_details(&self) -> Option<serde_json::Value> {
        Some(serde_json::json!({
            "reason": self.code,
            "rejected_candidates": self.rejected_candidates,
        }))
    }
}

/// Explicitly discretize an inclusive facility interval into available ranges.
///
/// The end point is included when it lies on the step grid (within floating
/// point tolerance); it is not appended as an off-grid extra station.  This
/// makes interval expansion reproducible and keeps every selected range on a
/// declared station grid.
pub fn discretize_truing_range_interval_v1(
    start_yd: f64,
    end_yd: f64,
    step_yd: f64,
) -> Result<Vec<f64>, TruingPlanErrorV1> {
    if !start_yd.is_finite() || start_yd <= 0.0 {
        return Err(TruingPlanErrorV1::invalid(
            "range interval start must be positive and finite",
        ));
    }
    if !end_yd.is_finite() || end_yd < start_yd {
        return Err(TruingPlanErrorV1::invalid(
            "range interval end must be finite and no smaller than its start",
        ));
    }
    if !step_yd.is_finite() || step_yd <= 0.0 {
        return Err(TruingPlanErrorV1::invalid(
            "range interval step must be positive and finite",
        ));
    }

    let step_count = ((end_yd - start_yd) / step_yd).floor();
    if !step_count.is_finite() || step_count > 100_000.0 {
        return Err(TruingPlanErrorV1::invalid(
            "range interval expands to more than 100001 candidates",
        ));
    }
    let count = step_count as usize + 1;
    let mut ranges = Vec::with_capacity(count);
    for i in 0..count {
        let range = start_yd + i as f64 * step_yd;
        if range <= end_yd + SEPARATION_TOLERANCE_YD {
            ranges.push(range.min(end_yd));
        }
    }
    Ok(ranges)
}

#[derive(Debug, Clone)]
struct CandidateInformation {
    input_index: usize,
    range_yd: f64,
    predicted_drop: f64,
    scaled_mv: f64,
    scaled_bc: f64,
}

#[derive(Debug, Clone, Copy)]
struct InformationMetrics {
    f00: f64,
    f01: f64,
    f11: f64,
    sensitivity_ratio: f64,
    condition_number: Option<f64>,
    min_singular: f64,
    max_singular: f64,
    log_determinant: Option<f64>,
    regularized_log_determinant: f64,
    expected_information_gain_nats: f64,
    joint_identifiable: bool,
}

impl InformationMetrics {
    fn from_indices(candidates: &[CandidateInformation], indices: &[usize]) -> Self {
        let mut f00 = 0.0;
        let mut f01 = 0.0;
        let mut f11 = 0.0;
        for &index in indices {
            let row = &candidates[index];
            f00 += row.scaled_mv * row.scaled_mv;
            f01 += row.scaled_mv * row.scaled_bc;
            f11 += row.scaled_bc * row.scaled_bc;
        }
        Self::from_information_matrix(f00, f01, f11)
    }

    fn from_information_matrix(f00: f64, f01: f64, f11: f64) -> Self {
        let trace = f00 + f11;
        let discriminant = ((f00 - f11) * (f00 - f11) + 4.0 * f01 * f01).sqrt();
        let lambda_max = ((trace + discriminant) * 0.5).max(0.0);
        let determinant = f00 * f11 - f01 * f01;
        // `trace - discriminant` loses nearly every significant bit for the
        // weak axis of an ill-conditioned design.  For a positive-semidefinite
        // 2x2 information matrix, det = lambda_max * lambda_min is the stable
        // way to recover that eigenvalue.  A tiny negative determinant from
        // accumulated roundoff represents the rank-one limit.
        let lambda_min = if lambda_max > 0.0 {
            (determinant / lambda_max).max(0.0)
        } else {
            0.0
        };
        let min_singular = lambda_min.sqrt();
        let max_singular = lambda_max.sqrt();

        let norm_mv = f00.max(0.0).sqrt();
        let norm_bc = f11.max(0.0).sqrt();
        let sensitivity_ratio = if norm_mv > 0.0 {
            norm_bc / norm_mv
        } else {
            0.0
        };
        let condition_number = if norm_mv > 0.0 && norm_bc > 0.0 {
            let correlation = (f01 / (norm_mv * norm_bc)).clamp(-1.0, 1.0).abs();
            if 1.0 - correlation > 1.0e-15 {
                Some((1.0 + correlation) / (1.0 - correlation))
            } else {
                None
            }
        } else {
            None
        };

        let log_determinant =
            (determinant > INFORMATION_EIGENVALUE_FLOOR).then(|| determinant.ln());
        let regularized_determinant = (1.0 + f00) * (1.0 + f11) - f01 * f01;
        let regularized_log_determinant = regularized_determinant.max(1.0).ln();
        let expected_information_gain_nats = 0.5 * regularized_log_determinant;
        let joint_identifiable = sensitivity_ratio >= TRUING_MIN_BC_SENSITIVITY_RATIO
            && condition_number.is_some_and(|condition| condition <= TRUING_MAX_CONDITION_NUMBER)
            && lambda_min > INFORMATION_EIGENVALUE_FLOOR;

        Self {
            f00,
            f01,
            f11,
            sensitivity_ratio,
            condition_number,
            min_singular,
            max_singular,
            log_determinant,
            regularized_log_determinant,
            expected_information_gain_nats,
            joint_identifiable,
        }
    }
}

#[derive(Debug, Clone)]
struct ScoredDesign {
    indices: Vec<usize>,
    information: InformationMetrics,
}

fn design_is_better(
    candidate: &ScoredDesign,
    incumbent: &ScoredDesign,
    candidates: &[CandidateInformation],
) -> bool {
    // An identifiable design always wins over an MV-only design.  Within the
    // same class, use regularized D-optimality, then maximin singular value,
    // then lower normalized condition, then lexicographically smaller ranges.
    if candidate.information.joint_identifiable != incumbent.information.joint_identifiable {
        return candidate.information.joint_identifiable;
    }
    let ordering = candidate
        .information
        .regularized_log_determinant
        .total_cmp(&incumbent.information.regularized_log_determinant);
    if !ordering.is_eq() {
        return ordering.is_gt();
    }
    let ordering = candidate
        .information
        .min_singular
        .total_cmp(&incumbent.information.min_singular);
    if !ordering.is_eq() {
        return ordering.is_gt();
    }
    match (
        candidate.information.condition_number,
        incumbent.information.condition_number,
    ) {
        (Some(a), Some(b)) if a.total_cmp(&b).is_ne() => return a < b,
        (Some(_), None) => return true,
        (None, Some(_)) => return false,
        _ => {}
    }

    let mut candidate_ranges: Vec<f64> = candidate
        .indices
        .iter()
        .map(|&i| candidates[i].range_yd)
        .collect();
    let mut incumbent_ranges: Vec<f64> = incumbent
        .indices
        .iter()
        .map(|&i| candidates[i].range_yd)
        .collect();
    candidate_ranges.sort_by(f64::total_cmp);
    incumbent_ranges.sort_by(f64::total_cmp);
    candidate_ranges
        .iter()
        .zip(&incumbent_ranges)
        .find_map(|(a, b)| {
            let ordering = a.total_cmp(b);
            (!ordering.is_eq()).then_some(ordering.is_lt())
        })
        .unwrap_or(false)
}

fn separated(
    candidates: &[CandidateInformation],
    indices: &[usize],
    minimum_separation_yd: f64,
) -> bool {
    for left in 0..indices.len() {
        for right in (left + 1)..indices.len() {
            let distance =
                (candidates[indices[left]].range_yd - candidates[indices[right]].range_yd).abs();
            if distance + SEPARATION_TOLERANCE_YD < minimum_separation_yd {
                return false;
            }
        }
    }
    true
}

fn maximum_compatible_count_with_fixed(
    candidates: &[CandidateInformation],
    fixed: &[usize],
    minimum_separation_yd: f64,
) -> usize {
    if !separated(candidates, fixed, minimum_separation_yd) {
        return 0;
    }
    let mut available: Vec<usize> = (0..candidates.len())
        .filter(|index| !fixed.contains(index))
        .filter(|index| {
            fixed.iter().all(|fixed_index| {
                (candidates[*index].range_yd - candidates[*fixed_index].range_yd).abs()
                    + SEPARATION_TOLERANCE_YD
                    >= minimum_separation_yd
            })
        })
        .collect();
    available.sort_by(|a, b| {
        candidates[*a]
            .range_yd
            .total_cmp(&candidates[*b].range_yd)
            .then_with(|| candidates[*a].input_index.cmp(&candidates[*b].input_index))
    });

    let mut count = fixed.len();
    let mut last_selected_range: Option<f64> = None;
    for index in available {
        let range = candidates[index].range_yd;
        if last_selected_range
            .is_none_or(|last| range - last + SEPARATION_TOLERANCE_YD >= minimum_separation_yd)
        {
            last_selected_range = Some(range);
            count += 1;
        }
    }
    count
}

fn binomial_capped(n: usize, k: usize, cap: u64) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut value: u128 = 1;
    for i in 1..=k {
        value = value * (n - k + i) as u128 / i as u128;
        if value > cap as u128 {
            return cap.saturating_add(1);
        }
    }
    value as u64
}

fn exhaustive_design(
    candidates: &[CandidateInformation],
    observation_count: usize,
    minimum_separation_yd: f64,
) -> (Option<ScoredDesign>, u64) {
    struct Search<'a> {
        candidates: &'a [CandidateInformation],
        target_count: usize,
        minimum_separation_yd: f64,
        evaluated: u64,
        best: Option<ScoredDesign>,
    }

    impl Search<'_> {
        fn visit(&mut self, start: usize, selected: &mut Vec<usize>) {
            if selected.len() == self.target_count {
                self.evaluated += 1;
                let design = ScoredDesign {
                    information: InformationMetrics::from_indices(self.candidates, selected),
                    indices: selected.clone(),
                };
                if self
                    .best
                    .as_ref()
                    .is_none_or(|best| design_is_better(&design, best, self.candidates))
                {
                    self.best = Some(design);
                }
                return;
            }

            let needed = self.target_count - selected.len();
            if self.candidates.len().saturating_sub(start) < needed {
                return;
            }
            for index in start..=self.candidates.len() - needed {
                if selected.iter().all(|selected_index| {
                    (self.candidates[index].range_yd - self.candidates[*selected_index].range_yd)
                        .abs()
                        + SEPARATION_TOLERANCE_YD
                        >= self.minimum_separation_yd
                }) {
                    selected.push(index);
                    self.visit(index + 1, selected);
                    selected.pop();
                }
            }
        }
    }

    let mut search = Search {
        candidates,
        target_count: observation_count,
        minimum_separation_yd,
        evaluated: 0,
        best: None,
    };
    search.visit(0, &mut Vec::with_capacity(observation_count));
    (search.best, search.evaluated)
}

fn greedy_exchange_design(
    candidates: &[CandidateInformation],
    observation_count: usize,
    minimum_separation_yd: f64,
) -> (Option<ScoredDesign>, u64) {
    let mut selected = Vec::with_capacity(observation_count);
    let mut evaluated = 0_u64;

    while selected.len() < observation_count {
        let mut best_addition: Option<ScoredDesign> = None;
        for index in 0..candidates.len() {
            if selected.contains(&index) {
                continue;
            }
            let mut trial = selected.clone();
            trial.push(index);
            if !separated(candidates, &trial, minimum_separation_yd)
                || maximum_compatible_count_with_fixed(candidates, &trial, minimum_separation_yd)
                    < observation_count
            {
                continue;
            }
            evaluated += 1;
            let design = ScoredDesign {
                information: InformationMetrics::from_indices(candidates, &trial),
                indices: trial,
            };
            if best_addition
                .as_ref()
                .is_none_or(|best| design_is_better(&design, best, candidates))
            {
                best_addition = Some(design);
            }
        }
        let Some(best_addition) = best_addition else {
            return (None, evaluated);
        };
        selected = best_addition.indices;
    }

    let mut current = ScoredDesign {
        information: InformationMetrics::from_indices(candidates, &selected),
        indices: selected,
    };
    // Each accepted exchange strictly improves a finite score tuple, so this
    // terminates.  The cap is a defensive bound against pathological exact-tie
    // cycles caused by future score changes.
    for _ in 0..candidates.len().saturating_mul(observation_count).max(1) {
        let mut best_exchange: Option<ScoredDesign> = None;
        for selected_position in 0..current.indices.len() {
            for replacement in 0..candidates.len() {
                if current.indices.contains(&replacement) {
                    continue;
                }
                let mut trial = current.indices.clone();
                trial[selected_position] = replacement;
                if !separated(candidates, &trial, minimum_separation_yd) {
                    continue;
                }
                evaluated += 1;
                let design = ScoredDesign {
                    information: InformationMetrics::from_indices(candidates, &trial),
                    indices: trial,
                };
                if design_is_better(&design, &current, candidates)
                    && best_exchange
                        .as_ref()
                        .is_none_or(|best| design_is_better(&design, best, candidates))
                {
                    best_exchange = Some(design);
                }
            }
        }
        let Some(next) = best_exchange else {
            break;
        };
        current = next;
    }
    (Some(current), evaluated)
}

fn optimize_design(
    candidates: &[CandidateInformation],
    observation_count: usize,
    minimum_separation_yd: f64,
    raw_combination_count: u64,
) -> (TruingPlanSearchStrategyV1, Option<ScoredDesign>, u64) {
    if raw_combination_count <= TRUING_PLAN_EXHAUSTIVE_COMBINATION_LIMIT_V1 {
        let (design, evaluated) =
            exhaustive_design(candidates, observation_count, minimum_separation_yd);
        (TruingPlanSearchStrategyV1::Exhaustive, design, evaluated)
    } else {
        let (design, evaluated) =
            greedy_exchange_design(candidates, observation_count, minimum_separation_yd);
        (
            TruingPlanSearchStrategyV1::GreedyExchange,
            design,
            evaluated,
        )
    }
}

fn selected_input_indices(
    design: &ScoredDesign,
    candidates: &[CandidateInformation],
) -> Vec<usize> {
    let mut indices: Vec<usize> = design
        .indices
        .iter()
        .map(|&index| candidates[index].input_index)
        .collect();
    indices.sort_unstable();
    indices
}

fn resolution_changes_design(
    candidates: &[CandidateInformation],
    baseline: &ScoredDesign,
    observation_count: usize,
    minimum_separation_yd: f64,
    raw_combination_count: u64,
    sigma_multiplier: f64,
) -> bool {
    // Sensitivities are divided by sigma, so changing sigma by `m` rescales
    // every information row by `1/m`.  Re-running the discrete optimizer is
    // cheap because the expensive trajectory/Jacobian values remain cached.
    let mut rescaled = candidates.to_vec();
    for candidate in &mut rescaled {
        candidate.scaled_mv /= sigma_multiplier;
        candidate.scaled_bc /= sigma_multiplier;
    }
    let (_, alternative, _) = optimize_design(
        &rescaled,
        observation_count,
        minimum_separation_yd,
        raw_combination_count,
    );
    alternative.is_some_and(|alternative| {
        alternative.information.joint_identifiable != baseline.information.joint_identifiable
            || selected_input_indices(&alternative, &rescaled)
                != selected_input_indices(baseline, candidates)
    })
}

fn scaled_candidate(
    input_index: usize,
    range_yd: f64,
    row: TruingJacobianRow,
    model: TruingModelInputsV1,
    measurement_sigma_1sd: f64,
) -> CandidateInformation {
    CandidateInformation {
        input_index,
        range_yd,
        predicted_drop: row.predicted_drop,
        scaled_mv: model.muzzle_velocity_fps * row.d_drop_d_mv / measurement_sigma_1sd,
        scaled_bc: model.ballistic_coefficient * row.d_drop_d_bc / measurement_sigma_1sd,
    }
}

/// Recommend an exact-size set of observation ranges for MV/BC truing.
///
/// This function is deterministic for a fixed request.  Every expensive
/// trajectory/Jacobian evaluation happens once per unique valid candidate.
/// Short-range or otherwise collinear facilities return an explicit MV-only
/// plan, not a falsely precise BC recommendation.
pub fn plan_truing_experiment_v1(
    request: &TruingExperimentPlanRequestV1,
) -> Result<TruingExperimentPlanV1, TruingPlanErrorV1> {
    request
        .model
        .validate()
        .map_err(TruingPlanErrorV1::invalid)?;
    if request.observation_count < 2 {
        return Err(TruingPlanErrorV1::invalid(
            "observation_count must be at least two for an MV/BC experiment design",
        ));
    }
    if request.candidate_ranges_yd.is_empty() {
        return Err(TruingPlanErrorV1::invalid(
            "candidate_ranges_yd must not be empty",
        ));
    }
    if !request.minimum_separation_yd.is_finite() || request.minimum_separation_yd < 0.0 {
        return Err(TruingPlanErrorV1::invalid(
            "minimum_separation_yd must be finite and non-negative",
        ));
    }
    if !request.measurement_sigma_1sd.is_finite() || request.measurement_sigma_1sd <= 0.0 {
        return Err(TruingPlanErrorV1::invalid(
            "measurement_sigma_1sd must be positive and finite",
        ));
    }

    let mut rejected = Vec::new();
    let mut valid_inputs: Vec<(usize, f64)> = Vec::new();
    for (input_index, &range_yd) in request.candidate_ranges_yd.iter().enumerate() {
        if !range_yd.is_finite() || range_yd <= 0.0 {
            rejected.push(RejectedTruingCandidateV1 {
                input_index,
                range_yd,
                reason: TruingCandidateRejectionReasonV1::InvalidRange,
                detail: "candidate range must be positive and finite".to_string(),
            });
            continue;
        }
        if let Some((duplicate_index, duplicate_range)) = valid_inputs
            .iter()
            .find(|(_, prior)| (range_yd - *prior).abs() < DUPLICATE_RANGE_TOLERANCE_YD)
        {
            rejected.push(RejectedTruingCandidateV1 {
                input_index,
                range_yd,
                reason: TruingCandidateRejectionReasonV1::DuplicateRange,
                detail: format!(
                    "duplicates candidate #{duplicate_index} at {duplicate_range:.6} yd"
                ),
            });
            continue;
        }
        valid_inputs.push((input_index, range_yd));
    }
    valid_inputs.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let ranges_yd: Vec<f64> = valid_inputs.iter().map(|(_, range)| *range).collect();
    let batched_rows = request
        .model
        .with_forward_model(request.drop_unit, |model| {
            truing_jacobian_rows(
                model,
                request.model.muzzle_velocity_fps,
                request.model.ballistic_coefficient,
                &ranges_yd,
                request.drop_unit,
            )
        });
    let mut candidates = Vec::with_capacity(valid_inputs.len());
    match batched_rows {
        Ok(rows) => {
            for ((input_index, range_yd), row) in valid_inputs.into_iter().zip(rows) {
                match row {
                    Some(row)
                        if row.predicted_drop.is_finite()
                            && row.d_drop_d_mv.is_finite()
                            && row.d_drop_d_bc.is_finite() =>
                    {
                        candidates.push(scaled_candidate(
                            input_index,
                            range_yd,
                            row,
                            request.model,
                            request.measurement_sigma_1sd,
                        ));
                    }
                    Some(_) => rejected.push(RejectedTruingCandidateV1 {
                        input_index,
                        range_yd,
                        reason: TruingCandidateRejectionReasonV1::Unreachable,
                        detail: "forward model returned a non-finite prediction or sensitivity"
                            .to_string(),
                    }),
                    None => rejected.push(RejectedTruingCandidateV1 {
                        input_index,
                        range_yd,
                        reason: TruingCandidateRejectionReasonV1::Unreachable,
                        detail: "nominal or perturbed trajectory did not reach candidate range"
                            .to_string(),
                    }),
                }
            }
        }
        Err(error) => {
            let detail = error.to_string();
            for (input_index, range_yd) in valid_inputs {
                rejected.push(RejectedTruingCandidateV1 {
                    input_index,
                    range_yd,
                    reason: TruingCandidateRejectionReasonV1::Unreachable,
                    detail: detail.clone(),
                });
            }
        }
    }
    rejected.sort_by_key(|candidate| candidate.input_index);

    if candidates.len() < request.observation_count {
        return Err(TruingPlanErrorV1 {
            code: TruingPlanErrorCodeV1::InsufficientReachableCandidates,
            message: format!(
                "requested {} observation stations but only {} unique reachable candidates remain",
                request.observation_count,
                candidates.len()
            ),
            rejected_candidates: rejected,
        });
    }
    if maximum_compatible_count_with_fixed(&candidates, &[], request.minimum_separation_yd)
        < request.observation_count
    {
        return Err(TruingPlanErrorV1 {
            code: TruingPlanErrorCodeV1::NoFeasibleDesign,
            message: format!(
                "no {}-station design satisfies the {:.3} yd minimum separation",
                request.observation_count, request.minimum_separation_yd
            ),
            rejected_candidates: rejected,
        });
    }

    let raw_combination_count =
        binomial_capped(candidates.len(), request.observation_count, u64::MAX - 1);
    let (search_strategy, design, evaluated_design_count) = optimize_design(
        &candidates,
        request.observation_count,
        request.minimum_separation_yd,
        raw_combination_count,
    );

    let Some(mut design) = design else {
        return Err(TruingPlanErrorV1 {
            code: TruingPlanErrorCodeV1::NoFeasibleDesign,
            message: format!(
                "no {}-station design satisfies the {:.3} yd minimum separation",
                request.observation_count, request.minimum_separation_yd
            ),
            rejected_candidates: rejected,
        });
    };
    design.indices.sort_by(|a, b| {
        candidates[*a]
            .range_yd
            .total_cmp(&candidates[*b].range_yd)
            .then_with(|| candidates[*a].input_index.cmp(&candidates[*b].input_index))
    });

    let selected_input_indices: Vec<usize> = design
        .indices
        .iter()
        .map(|&index| candidates[index].input_index)
        .collect();
    let selected_stations = design
        .indices
        .iter()
        .map(|&index| {
            let station = &candidates[index];
            let without_f00 = design.information.f00 - station.scaled_mv * station.scaled_mv;
            let without_f01 = design.information.f01 - station.scaled_mv * station.scaled_bc;
            let without_f11 = design.information.f11 - station.scaled_bc * station.scaled_bc;
            let without = InformationMetrics::from_information_matrix(
                without_f00.max(0.0),
                without_f01,
                without_f11.max(0.0),
            );
            TruingPlanStationV1 {
                input_index: station.input_index,
                range_yd: station.range_yd,
                predicted_drop: station.predicted_drop,
                scaled_mv_sensitivity: station.scaled_mv,
                scaled_bc_sensitivity: station.scaled_bc,
                leave_one_out_information_gain_nats: (design
                    .information
                    .expected_information_gain_nats
                    - without.expected_information_gain_nats)
                    .max(0.0),
            }
        })
        .collect();
    let unselected_candidate_ranges_yd = candidates
        .iter()
        .filter(|candidate| !selected_input_indices.contains(&candidate.input_index))
        .map(|candidate| candidate.range_yd)
        .collect();

    let mode = if design.information.joint_identifiable {
        TruingPlanModeV1::JointMvBc
    } else {
        TruingPlanModeV1::MvOnly
    };
    let mut warnings = vec![format!(
        "measurement model assumes independent 1-sigma {:.6} {} drop uncertainty at every station",
        request.measurement_sigma_1sd,
        request.drop_unit.label()
    )];
    warnings.push(
        "information-gain values use an identity reference information matrix in fractional MV/BC coordinates; they are experiment-design scores, not posterior intervals or an undeclared fit prior"
            .to_string(),
    );
    for sigma_multiplier in [0.5, 2.0] {
        if resolution_changes_design(
            &candidates,
            &design,
            request.observation_count,
            request.minimum_separation_yd,
            raw_combination_count,
            sigma_multiplier,
        ) {
            warnings.push(format!(
                "recommendation is sensitive to measurement resolution: at {sigma_multiplier:.1}x the declared sigma the optimizer changes the station set or joint/MV-only classification; rerun with that sigma before collecting data"
            ));
        }
    }
    if design.information.min_singular > INFORMATION_EIGENVALUE_FLOOR.sqrt() {
        let weak_axis_sigma = 1.0 / design.information.min_singular;
        warnings.push(format!(
            "local weak-axis fractional 1-sigma is {weak_axis_sigma:.4}; it scales linearly with measurement resolution ({:.4} at 0.5x sigma, {:.4} at 2x sigma)",
            weak_axis_sigma * 0.5,
            weak_axis_sigma * 2.0
        ));
    } else {
        warnings.push(
            "local weak-axis fractional 1-sigma is unbounded at this measurement resolution"
                .to_string(),
        );
    }
    if mode == TruingPlanModeV1::MvOnly {
        let condition = design
            .information
            .condition_number
            .map(|value| format!("{value:.3e}"))
            .unwrap_or_else(|| "unbounded".to_string());
        warnings.push(format!(
            "available ranges do not reliably separate MV from BC at the nominal profile (BC sensitivity ratio {:.4}, condition {condition}); use this as an MV-only design or add a longer-range station",
            design.information.sensitivity_ratio
        ));
    } else if design.information.sensitivity_ratio < TRUING_MIN_BC_SENSITIVITY_RATIO * 1.5 {
        warnings.push(format!(
            "BC sensitivity ratio {:.4} is close to the {:.2} identifiability threshold; small profile or measurement changes may make the design MV-only",
            design.information.sensitivity_ratio, TRUING_MIN_BC_SENSITIVITY_RATIO
        ));
    }

    Ok(TruingExperimentPlanV1 {
        mode,
        selected_stations,
        unselected_candidate_ranges_yd,
        rejected_candidates: rejected,
        information: TruingPlanInformationV1 {
            sensitivity_ratio: design.information.sensitivity_ratio,
            condition_number: design.information.condition_number,
            minimum_singular_value: design.information.min_singular,
            maximum_singular_value: design.information.max_singular,
            weak_axis_fractional_sigma_1sd: (design.information.min_singular
                > INFORMATION_EIGENVALUE_FLOOR.sqrt())
            .then(|| 1.0 / design.information.min_singular),
            log_determinant: design.information.log_determinant,
            expected_information_gain_nats: design.information.expected_information_gain_nats,
        },
        requested_observation_count: request.observation_count,
        minimum_separation_yd: request.minimum_separation_yd,
        measurement_sigma_1sd: request.measurement_sigma_1sd,
        measurement_drop_unit: request.drop_unit.label().to_string(),
        eligible_candidate_count: candidates.len(),
        raw_combination_count,
        evaluated_design_count,
        search_strategy,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthetic_candidates(count: usize) -> Vec<CandidateInformation> {
        (0..count)
            .map(|index| {
                let range = 100.0 + index as f64 * 50.0;
                CandidateInformation {
                    input_index: index,
                    range_yd: range,
                    predicted_drop: range / 100.0,
                    scaled_mv: 1.0 + index as f64 * 0.1,
                    scaled_bc: (index as f64 * 0.35).powi(2) + 0.05,
                }
            })
            .collect()
    }

    #[test]
    fn interval_discretization_stays_on_grid() {
        assert_eq!(
            discretize_truing_range_interval_v1(100.0, 325.0, 100.0).unwrap(),
            vec![100.0, 200.0, 300.0]
        );
        assert_eq!(
            discretize_truing_range_interval_v1(100.0, 300.0, 100.0).unwrap(),
            vec![100.0, 200.0, 300.0]
        );
    }

    #[test]
    fn exhaustive_search_honors_count_and_separation() {
        let candidates = synthetic_candidates(8);
        let (design, evaluated) = exhaustive_design(&candidates, 3, 150.0);
        let design = design.expect("feasible design");
        assert_eq!(design.indices.len(), 3);
        assert!(separated(&candidates, &design.indices, 150.0));
        assert!(evaluated > 0);
    }

    #[test]
    fn large_space_uses_deterministic_feasible_greedy_exchange() {
        let candidates = synthetic_candidates(25);
        assert!(
            binomial_capped(
                candidates.len(),
                6,
                TRUING_PLAN_EXHAUSTIVE_COMBINATION_LIMIT_V1
            ) > TRUING_PLAN_EXHAUSTIVE_COMBINATION_LIMIT_V1
        );
        let (first, _) = greedy_exchange_design(&candidates, 6, 100.0);
        let (second, _) = greedy_exchange_design(&candidates, 6, 100.0);
        let first = first.expect("feasible design");
        let second = second.expect("feasible design");
        assert_eq!(first.indices, second.indices);
        assert_eq!(first.indices.len(), 6);
        assert!(separated(&candidates, &first.indices, 100.0));
    }

    #[test]
    fn information_metrics_distinguish_collinear_rows() {
        let collinear = vec![
            CandidateInformation {
                input_index: 0,
                range_yd: 100.0,
                predicted_drop: 0.0,
                scaled_mv: 1.0,
                scaled_bc: 1.0,
            },
            CandidateInformation {
                input_index: 1,
                range_yd: 200.0,
                predicted_drop: 0.0,
                scaled_mv: 2.0,
                scaled_bc: 2.0,
            },
        ];
        let spread = vec![
            collinear[0].clone(),
            CandidateInformation {
                scaled_bc: 4.0,
                ..collinear[1].clone()
            },
        ];
        let weak = InformationMetrics::from_indices(&collinear, &[0, 1]);
        let strong = InformationMetrics::from_indices(&spread, &[0, 1]);
        assert!(weak.condition_number.is_none());
        assert_eq!(weak.min_singular, 0.0);
        assert!(strong.condition_number.is_some());
        assert!(strong.min_singular > weak.min_singular);
    }
}
