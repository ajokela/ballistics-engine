//! Uncertainty-aware joint muzzle-velocity / ballistic-coefficient truing.
//!
//! This is an opt-in companion to the historical truing API.  It performs a
//! weighted two-parameter MAP fit using absolute, caller-supplied observation
//! standard deviations and optional explicit independent normal priors.  The
//! legacy point fitter and its output remain untouched.
//!
//! V1 uses a local Laplace/Gaussian approximation at the constrained MAP.  It
//! is intentionally conservative about when that approximation is reported:
//! the fitted MAP is still returned if the information matrix is singular or
//! the optimum is on a parameter bound, but intervals and predictive bands are
//! replaced by a structured approximation failure.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::truing::{
    truing_jacobian_rows, DropUnit, TruingForwardModel, TruingModelInputsV1, TRUING_BC_MAX,
    TRUING_BC_MIN, TRUING_MAX_CONDITION_NUMBER, TRUING_MIN_BC_SENSITIVITY_RATIO, TRUING_MV_MAX_FPS,
    TRUING_MV_MIN_FPS,
};

/// Schema version for the uncertainty-aware library result.
pub const TRUING_UNCERTAINTY_SCHEMA_VERSION_V1: u32 = 1;

/// The interval probability reported by V1.
pub const TRUING_UNCERTAINTY_INTERVAL_LEVEL_V1: f64 = 0.95;

/// Iteration cap for the separately opt-in weighted joint MAP optimizer.
///
/// Weakly identifiable short-range likelihoods have long, shallow valleys and
/// need more iterations than the legacy heuristic-gated point fitter, which
/// refuses those joint fits entirely.
pub const TRUING_UNCERTAINTY_MAX_ITERS_V1: usize = 100;

const NORMAL_95_TWO_SIDED_Z: f64 = 1.959_963_984_540_054;
// Fixed coordinate scales make the two columns of the normal equations
// numerically comparable without changing the physical result.
const MV_COORDINATE_SCALE_FPS: f64 = 100.0;
const BC_COORDINATE_SCALE: f64 = 0.1;
const INFORMATION_RELATIVE_EIGEN_TOLERANCE: f64 = 1.0e-12;
const MAP_SCALED_GRADIENT_TOLERANCE: f64 = 1.0e-6;
// The real trajectory surface has a much finer numerical texture than the
// deliberately broad finite-difference steps used to estimate its covariance.
// If LM cannot satisfy the gradient test, a direct-objective pattern poll must
// reduce every scaled step to this radius and find no improvement larger than
// the absolute chi-square tolerance before the MAP is accepted.
const MAP_OBJECTIVE_INITIAL_POLL_RADIUS: f64 = 1.0e-2;
const MAP_OBJECTIVE_MIN_POLL_RADIUS: f64 = 1.0e-7;
const MAP_OBJECTIVE_IMPROVEMENT_TOLERANCE: f64 = 1.0e-8;
const MAP_OBJECTIVE_MAX_POLL_EVALUATIONS: usize = 1_024;

/// One measured drop and its absolute one-standard-deviation uncertainty.
///
/// `drop` and `sigma` are both expressed in the request's [`DropUnit`].  The
/// range is always in internal yards.  Repeated ranges are accepted: repeated
/// shots tighten the same parameter combination and permit chi-square
/// consistency checks, but they do not, by themselves, separate MV from BC.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WeightedTruingObservationV1 {
    pub range_yd: f64,
    pub drop: f64,
    pub sigma: f64,
}

/// An explicit independent normal prior, in the parameter's physical units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NormalPriorV1 {
    pub mean: f64,
    pub sigma: f64,
}

/// Optional independent priors used by the V1 MAP approximation.
///
/// There are no hidden priors.  `None` means that parameter contributes no
/// prior precision or penalty.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TruingPriorsV1 {
    pub muzzle_velocity_fps: Option<NormalPriorV1>,
    pub ballistic_coefficient: Option<NormalPriorV1>,
}

/// A range at which to propagate parameter uncertainty into predicted drop.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TruingPredictionRequestV1 {
    pub range_yd: f64,
    /// Optional absolute measurement sigma in the request's drop unit.  When
    /// present, V1 reports both the latent model band and the wider band for a
    /// future observation.  It is never inferred from residuals.
    pub future_observation_sigma: Option<f64>,
}

/// Complete request for weighted joint MV+BC MAP truing.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UncertaintyTruingRequestV1 {
    pub model: TruingModelInputsV1,
    pub drop_unit: DropUnit,
    pub observations: Vec<WeightedTruingObservationV1>,
    pub priors: TruingPriorsV1,
    pub predictions: Vec<TruingPredictionRequestV1>,
}

/// A two-sided interval from the local Gaussian approximation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct GaussianIntervalV1 {
    pub estimate: f64,
    pub standard_deviation: f64,
    pub lower: f64,
    pub upper: f64,
    pub probability: f64,
}

impl GaussianIntervalV1 {
    fn from_variance(estimate: f64, variance: f64) -> Option<Self> {
        if !estimate.is_finite() || !variance.is_finite() || variance < 0.0 {
            return None;
        }
        let standard_deviation = variance.sqrt();
        let half_width = NORMAL_95_TWO_SIDED_Z * standard_deviation;
        Some(Self {
            estimate,
            standard_deviation,
            lower: estimate - half_width,
            upper: estimate + half_width,
            probability: TRUING_UNCERTAINTY_INTERVAL_LEVEL_V1,
        })
    }
}

/// Physical-coordinate posterior covariance for `(MV fps, BC)`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TruingCovarianceV1 {
    pub mv_variance_fps2: f64,
    pub mv_bc_covariance_fps: f64,
    pub bc_variance: f64,
}

/// Successful local Gaussian approximation around the MAP.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TruingGaussianApproximationV1 {
    pub covariance: TruingCovarianceV1,
    pub muzzle_velocity_interval_95: GaussianIntervalV1,
    pub ballistic_coefficient_interval_95: GaussianIntervalV1,
    /// Posterior correlation between fitted MV and BC.
    pub mv_bc_correlation: f64,
    /// Condition number of the posterior information matrix in the documented
    /// scaled coordinates (`100 fps`, `0.1 BC`).
    pub scaled_information_condition_number: f64,
}

/// Why a local Gaussian approximation was withheld even though a MAP is
/// available.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingApproximationFailureCodeV1 {
    OptimizerDidNotConverge,
    MapAtParameterBound,
    RankDeficientInformation,
    NonFiniteInformation,
}

/// Structured approximation failure.  Consumers can branch on `code` while
/// still presenting the explanatory `message`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct TruingApproximationFailureV1 {
    pub code: TruingApproximationFailureCodeV1,
    pub message: String,
}

/// Availability of the Laplace/Gaussian approximation.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "snake_case", tag = "status", content = "details")]
pub enum TruingApproximationV1 {
    Available(TruingGaussianApproximationV1),
    Unavailable(TruingApproximationFailureV1),
}

/// Input and residual diagnostics at one observation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct WeightedTruingObservationResultV1 {
    pub range_yd: f64,
    pub observed_drop: f64,
    pub sigma: f64,
    pub predicted_drop: f64,
    pub residual: f64,
    pub standardized_residual: f64,
}

/// Parameter-propagated drop uncertainty at one requested range.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TruingPredictiveBandV1 {
    pub range_yd: f64,
    pub predicted_drop: f64,
    /// Band for the latent model prediction.  `None` means the Gaussian
    /// approximation was unavailable.
    pub latent_interval_95: Option<GaussianIntervalV1>,
    /// Band for a future measured impact, including the explicitly supplied
    /// `future_observation_sigma`.  This remains `None` when no future sigma was
    /// requested or the Gaussian approximation was unavailable.
    pub future_observation_interval_95: Option<GaussianIntervalV1>,
}

/// Machine-readable warning categories emitted alongside the fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingUncertaintyWarningCodeV1 {
    OptimizerDidNotConverge,
    ObjectiveMeshConvergence,
    WeakBcSensitivity,
    IllConditionedData,
    MvPriorDominated,
    BcPriorDominated,
    GaussianApproximationUnavailable,
    IntervalCrossesFitBounds,
    LowEffectiveDegreesOfFreedom,
    PredictionOutsideObservedDomain,
}

/// Numerical criterion that verified the reported MAP.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TruingMapConvergenceCriterionV1 {
    /// The scaled Gauss-Newton half-gradient met its tolerance.
    ScaledGradient,
    /// A deterministic direct-objective pattern search found no material
    /// improvement at its minimum mesh radius.
    ObjectiveMesh,
}

/// A warning with stable code and human-readable context.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct TruingUncertaintyWarningV1 {
    pub code: TruingUncertaintyWarningCodeV1,
    pub message: String,
}

/// Fit and approximation diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TruingUncertaintyDiagnosticsV1 {
    /// Data chi-square, `sum((prediction - observation) / sigma)^2`.
    pub chi_square: f64,
    /// Explicit prior contribution to the penalized chi-square.
    pub prior_penalty: f64,
    /// `chi_square + prior_penalty`, minimized by the MAP fitter.
    pub penalized_chi_square: f64,
    /// Effective number of parameters learned from the observations,
    /// `trace(I_data * I_posterior^-1)`.
    pub effective_parameter_count: Option<f64>,
    /// Effective residual degrees of freedom, `n - effective_parameter_count`.
    pub effective_degrees_of_freedom: Option<f64>,
    /// Data chi-square divided by effective residual degrees of freedom.
    pub reduced_chi_square: Option<f64>,
    /// Fractional BC sensitivity relative to fractional MV sensitivity.
    pub bc_sensitivity_ratio: f64,
    /// Collinearity condition diagnostic for the weighted data Jacobian.
    pub data_condition_number: f64,
    /// Infinity norm of the penalized-chi-square half-gradient in scaled
    /// coordinates at the reported MAP.  This is an informative local-model
    /// diagnostic.  It is the primary convergence test, but the deliberately
    /// broad finite-difference stencil can disagree with the trajectory
    /// solver's much finer numerical surface.  In that case
    /// [`TruingMapConvergenceCriterionV1::ObjectiveMesh`] records a direct-
    /// objective verification instead.
    pub map_scaled_gradient_inf_norm: f64,
    /// Criterion that verified the reported MAP, or `None` if the optimizer
    /// did not converge.
    pub map_convergence_criterion: Option<TruingMapConvergenceCriterionV1>,
    /// Final scaled pattern-poll radius for objective-mesh convergence.
    /// `100 fps` and `0.1 BC` are one scaled unit.  `None` when the gradient
    /// criterion was used or the poll did not reach its resolution target.
    pub map_objective_poll_radius: Option<f64>,
    /// Largest direct chi-square improvement found in the final pattern poll.
    /// Objective-mesh convergence requires this to be no greater than
    /// `1e-8` in absolute penalized chi-square units.
    pub map_max_objective_poll_improvement: Option<f64>,
    /// Number of direct objective evaluations used by the fallback poll.
    pub map_objective_poll_evaluations: usize,
}

/// Result of uncertainty-aware joint MV+BC truing.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct UncertaintyTruingReportV1 {
    pub schema_version: u32,
    pub drop_unit: DropUnit,
    pub map_muzzle_velocity_fps: f64,
    pub map_ballistic_coefficient: f64,
    pub iterations: usize,
    pub converged: bool,
    pub priors: TruingPriorsV1,
    pub observations: Vec<WeightedTruingObservationResultV1>,
    pub diagnostics: TruingUncertaintyDiagnosticsV1,
    pub approximation: TruingApproximationV1,
    pub predictive_bands: Vec<TruingPredictiveBandV1>,
    pub warnings: Vec<TruingUncertaintyWarningV1>,
}

/// Input or forward-model failure before a report can be produced.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum UncertaintyTruingErrorV1 {
    #[error("invalid uncertainty-truing request: {0}")]
    InvalidInput(String),
    #[error("truing forward model failed: {0}")]
    ForwardModel(String),
}

/// A symmetric 2x2 matrix -- posterior/information matrices in this module, and (0.33.0
/// decision-support Task 10, MBA-1347) an impact-space covariance in the `error_budget` module.
///
/// `pub(crate)` together with its three fields and `add_assign`: `error_budget` needs to
/// accumulate several declared sources' contributions into one covariance and then read its
/// eigenvalues, reusing this type's existing, reviewed eigenvalue arithmetic (see
/// `largest_smallest_eigenvalues` below) rather than writing a third copy of it
/// (`crate::monte_carlo::calculate_confidence_ellipse` already has its own, sample-based copy).
/// This is the same precedent as `kernel_solve_error`/`wind_reference_of` in
/// `crate::perturbation`, widened from private for the identical reason in earlier tasks on this
/// branch. Widening is behavior-preserving: every existing use of this type is within this same
/// module, where field/method privacy was never enforced against it anyway.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Symmetric2 {
    pub(crate) a00: f64,
    pub(crate) a01: f64,
    pub(crate) a11: f64,
}

impl Symmetric2 {
    fn determinant(self) -> f64 {
        self.a00 * self.a11 - self.a01 * self.a01
    }

    pub(crate) fn add_assign(&mut self, rhs: Self) {
        self.a00 += rhs.a00;
        self.a01 += rhs.a01;
        self.a11 += rhs.a11;
    }

    /// Eigenvalues of this symmetric 2x2 matrix, largest first, each clamped to non-negative.
    /// Never fails.
    ///
    /// `pub(crate)` (0.33.0 decision-support Task 10, MBA-1347): shares `inverse_with_condition`
    /// below's eigenvalue arithmetic (trace, `hypot`-based discriminant, and the
    /// `determinant = largest * smallest` division -- more accurate than subtracting two nearly
    /// equal values in the naive quadratic formula) but with different error semantics.
    /// `inverse_with_condition` treats a singular or non-positive-definite matrix as a hard
    /// failure, because it specifically inverts an information matrix that this module's MAP fit
    /// needs to be well-conditioned. An error-budget impact covariance is routinely and
    /// legitimately singular (exactly one nonzero-sigma source produces an exact rank-1
    /// covariance, since it is a single outer product `sigma^2 * v * v^T`) or exactly zero
    /// (every declared sigma is zero) -- both normal, expected inputs there, not failures. This
    /// method never errors: a negative discriminant or determinant can only arise here from
    /// floating-point rounding on a matrix that is mathematically positive-semidefinite by
    /// construction (a sum of such outer products), so it is clamped to zero rather than
    /// reported as a failure.
    pub(crate) fn largest_smallest_eigenvalues(self) -> (f64, f64) {
        let trace = self.a00 + self.a11;
        let discriminant = (self.a00 - self.a11).hypot(2.0 * self.a01);
        let largest = (0.5 * (trace + discriminant)).max(0.0);
        let determinant = self.determinant();
        let smallest = if largest > 0.0 {
            (determinant / largest).max(0.0)
        } else {
            0.0
        };
        (largest, smallest)
    }

    fn inverse_with_condition(self) -> Result<(Self, f64), TruingApproximationFailureCodeV1> {
        if !self.a00.is_finite() || !self.a01.is_finite() || !self.a11.is_finite() {
            return Err(TruingApproximationFailureCodeV1::NonFiniteInformation);
        }
        let trace = self.a00 + self.a11;
        let discriminant = (self.a00 - self.a11).hypot(2.0 * self.a01);
        let largest = 0.5 * (trace + discriminant);
        let determinant = self.determinant();
        if !largest.is_finite() || !determinant.is_finite() {
            return Err(TruingApproximationFailureCodeV1::NonFiniteInformation);
        }
        if largest <= 0.0 || determinant <= 0.0 {
            return Err(TruingApproximationFailureCodeV1::RankDeficientInformation);
        }
        // det = lambda_max * lambda_min is more accurate than subtracting two
        // nearly equal values from the quadratic formula.
        let smallest = determinant / largest;
        if !smallest.is_finite()
            || smallest <= 0.0
            || smallest / largest <= INFORMATION_RELATIVE_EIGEN_TOLERANCE
        {
            return Err(TruingApproximationFailureCodeV1::RankDeficientInformation);
        }
        let inverse = Self {
            a00: self.a11 / determinant,
            a01: -self.a01 / determinant,
            a11: self.a00 / determinant,
        };
        if !inverse.a00.is_finite() || !inverse.a01.is_finite() || !inverse.a11.is_finite() {
            return Err(TruingApproximationFailureCodeV1::NonFiniteInformation);
        }
        Ok((inverse, largest / smallest))
    }
}

struct Evaluation {
    data_information: Symmetric2,
    posterior_information: Symmetric2,
    gradient: [f64; 2],
    observation_results: Vec<WeightedTruingObservationResultV1>,
    chi_square: f64,
    prior_penalty: f64,
    bc_sensitivity_ratio: f64,
    data_condition_number: f64,
}

#[derive(Debug, Clone, Copy)]
struct MapFitResult {
    mv: f64,
    bc: f64,
    iterations: usize,
    convergence_criterion: Option<TruingMapConvergenceCriterionV1>,
    objective_poll_radius: Option<f64>,
    max_objective_poll_improvement: Option<f64>,
    objective_poll_evaluations: usize,
}

#[derive(Debug, Clone, Copy)]
struct ObjectivePollResult {
    mv: f64,
    bc: f64,
    converged: bool,
    final_radius: f64,
    max_final_improvement: f64,
    evaluations: usize,
}

impl Evaluation {
    fn penalized_chi_square(&self) -> f64 {
        self.chi_square + self.prior_penalty
    }
}

/// Run the opt-in uncertainty-aware joint MV+BC truing path.
///
/// Observation sigmas are absolute and known.  Consequently posterior
/// covariance is the inverse likelihood-plus-prior information matrix and is
/// **not** rescaled by residual RMS.  No identifiability heuristic changes the
/// requested two-parameter fit; weak or collinear data instead produce broad
/// uncertainty, warnings, or a structured approximation failure.
pub fn run_uncertainty_truing_v1(
    request: &UncertaintyTruingRequestV1,
) -> Result<UncertaintyTruingReportV1, UncertaintyTruingErrorV1> {
    validate_request(request)?;

    request
        .model
        .with_forward_model(request.drop_unit, |model| run_with_model(request, model))
}

fn run_with_model(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
) -> Result<UncertaintyTruingReportV1, UncertaintyTruingErrorV1> {
    let fit = fit_map(request, model)?;
    let map_mv = fit.mv;
    let map_bc = fit.bc;
    let iterations = fit.iterations;
    let converged = fit.convergence_criterion.is_some();
    let evaluation = evaluate(request, model, map_mv, map_bc)?;
    let mut warnings = Vec::new();

    if !converged {
        let gradient_norm = evaluation.gradient[0]
            .abs()
            .max(evaluation.gradient[1].abs());
        push_warning(
            &mut warnings,
            TruingUncertaintyWarningCodeV1::OptimizerDidNotConverge,
            format!(
                "joint MV+BC MAP optimizer stopped after {iterations} iterations with scaled gradient {gradient_norm:.3e} above {MAP_SCALED_GRADIENT_TOLERANCE:.1e}"
            ),
        );
    } else if fit.convergence_criterion == Some(TruingMapConvergenceCriterionV1::ObjectiveMesh) {
        let radius = fit.objective_poll_radius.unwrap_or(f64::NAN);
        let improvement = fit.max_objective_poll_improvement.unwrap_or(f64::NAN);
        push_warning(
            &mut warnings,
            TruingUncertaintyWarningCodeV1::ObjectiveMeshConvergence,
            format!(
                "LM's broad-stencil gradient test stalled (final norm {:.3e}); a direct-objective pattern search found no penalized-chi-square improvement above {:.1e} at scaled radius {radius:.1e} (largest observed {improvement:.3e})",
                evaluation.gradient[0]
                    .abs()
                    .max(evaluation.gradient[1].abs()),
                MAP_OBJECTIVE_IMPROVEMENT_TOLERANCE,
            ),
        );
    }
    if evaluation.bc_sensitivity_ratio < TRUING_MIN_BC_SENSITIVITY_RATIO {
        push_warning(
            &mut warnings,
            TruingUncertaintyWarningCodeV1::WeakBcSensitivity,
            format!(
                "observations weakly constrain BC: fractional sensitivity ratio {:.4} is below {:.2}",
                evaluation.bc_sensitivity_ratio, TRUING_MIN_BC_SENSITIVITY_RATIO
            ),
        );
    }
    if !evaluation.data_condition_number.is_finite()
        || evaluation.data_condition_number > TRUING_MAX_CONDITION_NUMBER
    {
        push_warning(
            &mut warnings,
            TruingUncertaintyWarningCodeV1::IllConditionedData,
            format!(
                "weighted observation Jacobian cannot cleanly separate MV from BC (condition {:.3e})",
                evaluation.data_condition_number
            ),
        );
    }

    let at_bound = parameter_at_bound(map_mv, map_bc);
    let approximation_result = if !converged {
        Err(TruingApproximationFailureV1 {
            code: TruingApproximationFailureCodeV1::OptimizerDidNotConverge,
            message: "MAP optimizer did not converge; covariance around an unverified stationary point is withheld".to_string(),
        })
    } else if at_bound {
        Err(TruingApproximationFailureV1 {
            code: TruingApproximationFailureCodeV1::MapAtParameterBound,
            message: "MAP lies on or numerically near a fit bound; an unconstrained Gaussian approximation would be misleading".to_string(),
        })
    } else {
        build_approximation(map_mv, map_bc, evaluation.posterior_information)
    };

    let (approximation, covariance_q) = match approximation_result {
        Ok((gaussian, covariance_q)) => {
            warn_prior_dominance(request, &gaussian, &mut warnings);
            if gaussian.muzzle_velocity_interval_95.lower < TRUING_MV_MIN_FPS
                || gaussian.muzzle_velocity_interval_95.upper > TRUING_MV_MAX_FPS
                || gaussian.ballistic_coefficient_interval_95.lower < TRUING_BC_MIN
                || gaussian.ballistic_coefficient_interval_95.upper > TRUING_BC_MAX
            {
                push_warning(
                    &mut warnings,
                    TruingUncertaintyWarningCodeV1::IntervalCrossesFitBounds,
                    "local Gaussian interval crosses a constrained fit bound; interpret its tails cautiously".to_string(),
                );
            }
            (
                TruingApproximationV1::Available(gaussian),
                Some(covariance_q),
            )
        }
        Err(failure) => {
            push_warning(
                &mut warnings,
                TruingUncertaintyWarningCodeV1::GaussianApproximationUnavailable,
                failure.message.clone(),
            );
            (TruingApproximationV1::Unavailable(failure), None)
        }
    };

    let effective_parameter_count = covariance_q.map(|covariance| {
        (evaluation.data_information.a00 * covariance.a00
            + 2.0 * evaluation.data_information.a01 * covariance.a01
            + evaluation.data_information.a11 * covariance.a11)
            .clamp(0.0, 2.0)
    });
    let effective_degrees_of_freedom =
        effective_parameter_count.map(|count| request.observations.len() as f64 - count);
    let reduced_chi_square = effective_degrees_of_freedom
        .filter(|dof| *dof > f64::EPSILON)
        .map(|dof| evaluation.chi_square / dof);
    if effective_degrees_of_freedom.is_some_and(|dof| dof <= 1.0) {
        push_warning(
            &mut warnings,
            TruingUncertaintyWarningCodeV1::LowEffectiveDegreesOfFreedom,
            "one or fewer effective residual degrees of freedom: residual fit quality is weakly validated".to_string(),
        );
    }

    let predictive_bands =
        build_predictive_bands(request, model, map_mv, map_bc, covariance_q, &mut warnings)?;

    let penalized_chi_square = evaluation.penalized_chi_square();

    Ok(UncertaintyTruingReportV1 {
        schema_version: TRUING_UNCERTAINTY_SCHEMA_VERSION_V1,
        drop_unit: request.drop_unit,
        map_muzzle_velocity_fps: map_mv,
        map_ballistic_coefficient: map_bc,
        iterations,
        converged,
        priors: request.priors,
        observations: evaluation.observation_results,
        diagnostics: TruingUncertaintyDiagnosticsV1 {
            chi_square: evaluation.chi_square,
            prior_penalty: evaluation.prior_penalty,
            penalized_chi_square,
            effective_parameter_count,
            effective_degrees_of_freedom,
            reduced_chi_square,
            bc_sensitivity_ratio: evaluation.bc_sensitivity_ratio,
            data_condition_number: evaluation.data_condition_number,
            map_scaled_gradient_inf_norm: evaluation.gradient[0]
                .abs()
                .max(evaluation.gradient[1].abs()),
            map_convergence_criterion: fit.convergence_criterion,
            map_objective_poll_radius: fit.objective_poll_radius,
            map_max_objective_poll_improvement: fit.max_objective_poll_improvement,
            map_objective_poll_evaluations: fit.objective_poll_evaluations,
        },
        approximation,
        predictive_bands,
        warnings,
    })
}

fn validate_request(request: &UncertaintyTruingRequestV1) -> Result<(), UncertaintyTruingErrorV1> {
    request
        .model
        .validate()
        .map_err(UncertaintyTruingErrorV1::InvalidInput)?;
    if request.observations.len() < 2 {
        return Err(UncertaintyTruingErrorV1::InvalidInput(
            "at least two weighted observations are required for a joint MV+BC fit".to_string(),
        ));
    }
    for (index, observation) in request.observations.iter().enumerate() {
        if !observation.range_yd.is_finite() || observation.range_yd <= 0.0 {
            return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
                "observation {} range must be positive and finite",
                index + 1
            )));
        }
        if !observation.drop.is_finite() {
            return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
                "observation {} drop must be finite",
                index + 1
            )));
        }
        if !observation.sigma.is_finite() || observation.sigma <= 0.0 {
            return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
                "observation {} sigma must be positive and finite",
                index + 1
            )));
        }
    }
    validate_prior(
        "muzzle-velocity",
        request.priors.muzzle_velocity_fps,
        TRUING_MV_MIN_FPS,
        TRUING_MV_MAX_FPS,
    )?;
    validate_prior(
        "ballistic-coefficient",
        request.priors.ballistic_coefficient,
        TRUING_BC_MIN,
        TRUING_BC_MAX,
    )?;
    for (index, prediction) in request.predictions.iter().enumerate() {
        if !prediction.range_yd.is_finite() || prediction.range_yd <= 0.0 {
            return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
                "prediction {} range must be positive and finite",
                index + 1
            )));
        }
        if prediction
            .future_observation_sigma
            .is_some_and(|sigma| !sigma.is_finite() || sigma <= 0.0)
        {
            return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
                "prediction {} future-observation sigma must be positive and finite",
                index + 1
            )));
        }
    }
    Ok(())
}

fn validate_prior(
    name: &str,
    prior: Option<NormalPriorV1>,
    lower: f64,
    upper: f64,
) -> Result<(), UncertaintyTruingErrorV1> {
    let Some(prior) = prior else {
        return Ok(());
    };
    if !prior.mean.is_finite() || !(lower..=upper).contains(&prior.mean) {
        return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
            "{name} prior mean must be finite and within {lower}..={upper}"
        )));
    }
    if !prior.sigma.is_finite() || prior.sigma <= 0.0 {
        return Err(UncertaintyTruingErrorV1::InvalidInput(format!(
            "{name} prior sigma must be positive and finite"
        )));
    }
    Ok(())
}

fn fit_map(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
) -> Result<MapFitResult, UncertaintyTruingErrorV1> {
    let mut mv = request.model.muzzle_velocity_fps;
    let mut bc = request.model.ballistic_coefficient;
    let mut lambda = 1.0e-6;
    let mut current = objective(request, model, mv, bc)?;
    let mut iterations = 0;
    let mut convergence_criterion = None;

    for _ in 0..TRUING_UNCERTAINTY_MAX_ITERS_V1 {
        iterations += 1;
        let evaluation = evaluate(request, model, mv, bc)?;
        if evaluation.gradient[0]
            .abs()
            .max(evaluation.gradient[1].abs())
            <= MAP_SCALED_GRADIENT_TOLERANCE
        {
            convergence_criterion = Some(TruingMapConvergenceCriterionV1::ScaledGradient);
            break;
        }

        let mut accepted = false;
        for _ in 0..30 {
            let information = evaluation.posterior_information;
            let damped = Symmetric2 {
                a00: information.a00 + lambda * information.a00.max(1.0e-12),
                a01: information.a01,
                a11: information.a11 + lambda * information.a11.max(1.0e-12),
            };
            let determinant = damped.determinant();
            if !determinant.is_finite() || determinant.abs() < 1.0e-24 {
                lambda *= 10.0;
                continue;
            }
            let delta_mv_coordinate = -(damped.a11 * evaluation.gradient[0]
                - damped.a01 * evaluation.gradient[1])
                / determinant;
            let delta_bc_coordinate = -(-damped.a01 * evaluation.gradient[0]
                + damped.a00 * evaluation.gradient[1])
                / determinant;
            let next_mv = (mv + MV_COORDINATE_SCALE_FPS * delta_mv_coordinate)
                .clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
            let next_bc = (bc + BC_COORDINATE_SCALE * delta_bc_coordinate)
                .clamp(TRUING_BC_MIN, TRUING_BC_MAX);
            let next = objective(request, model, next_mv, next_bc)?;
            if objective_improvement_is_material(current, next) {
                mv = next_mv;
                bc = next_bc;
                current = next;
                lambda = (lambda * 0.5).max(1.0e-12);
                accepted = true;
                break;
            }
            // If clamping or floating-point resolution leaves the candidate at
            // the current optimum, increasing damping cannot reveal a lower
            // point.  The direct-objective poll below still verifies the full
            // two-dimensional neighborhood before accepting convergence.
            if next_mv == mv
                && next_bc == bc
                && evaluation.gradient[0]
                    .abs()
                    .max(evaluation.gradient[1].abs())
                    <= MAP_SCALED_GRADIENT_TOLERANCE
            {
                convergence_criterion = Some(TruingMapConvergenceCriterionV1::ScaledGradient);
                break;
            }
            lambda *= 4.0;
        }
        if convergence_criterion.is_some() {
            break;
        }
        if !accepted {
            break;
        }
    }
    if convergence_criterion.is_none() {
        let final_evaluation = evaluate(request, model, mv, bc)?;
        if final_evaluation.gradient[0]
            .abs()
            .max(final_evaluation.gradient[1].abs())
            <= MAP_SCALED_GRADIENT_TOLERANCE
        {
            convergence_criterion = Some(TruingMapConvergenceCriterionV1::ScaledGradient);
        }
    }

    let mut objective_poll_radius = None;
    let mut max_objective_poll_improvement = None;
    let mut objective_poll_evaluations = 0;
    if convergence_criterion.is_none() {
        let poll = polish_and_verify_objective_mesh(request, model, mv, bc, current)?;
        mv = poll.mv;
        bc = poll.bc;
        objective_poll_evaluations = poll.evaluations;
        if poll.converged {
            convergence_criterion = Some(TruingMapConvergenceCriterionV1::ObjectiveMesh);
            objective_poll_radius = Some(poll.final_radius);
            max_objective_poll_improvement = Some(poll.max_final_improvement);
        }
    }

    Ok(MapFitResult {
        mv,
        bc,
        iterations,
        convergence_criterion,
        objective_poll_radius,
        max_objective_poll_improvement,
        objective_poll_evaluations,
    })
}

fn objective_improvement_is_material(current: f64, candidate: f64) -> bool {
    current - candidate > MAP_OBJECTIVE_IMPROVEMENT_TOLERANCE
}

/// Polish an LM candidate against the actual penalized chi-square surface and
/// verify local convergence independently of the finite-difference Jacobian.
///
/// The eight polling directions contain a positive spanning coordinate basis;
/// the local information eigenvectors accelerate progress along the strongly
/// correlated MV/BC valley.
/// A point is accepted only after no direction yields a material improvement
/// at the minimum mesh radius.
fn polish_and_verify_objective_mesh(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
    mut mv: f64,
    mut bc: f64,
    mut current: f64,
) -> Result<ObjectivePollResult, UncertaintyTruingErrorV1> {
    let information = evaluate(request, model, mv, bc)?.posterior_information;
    let angle = 0.5 * (2.0 * information.a01).atan2(information.a00 - information.a11);
    let (sin, cos) = angle.sin_cos();
    // (cos, sin) is the large-information direction and (-sin, cos) the weak
    // direction.  Both are unit vectors in scaled optimizer coordinates.
    let directions = [
        (1.0, 0.0),
        (-1.0, 0.0),
        (0.0, 1.0),
        (0.0, -1.0),
        (cos, sin),
        (-cos, -sin),
        (-sin, cos),
        (sin, -cos),
    ];

    let mut radius = MAP_OBJECTIVE_INITIAL_POLL_RADIUS;
    let mut evaluations = 0;
    let mut max_final_improvement = f64::INFINITY;
    while evaluations < MAP_OBJECTIVE_MAX_POLL_EVALUATIONS {
        let mut best = current;
        let mut best_point = (mv, bc);
        for (mv_direction, bc_direction) in directions {
            if evaluations >= MAP_OBJECTIVE_MAX_POLL_EVALUATIONS {
                break;
            }
            let candidate_mv = (mv + radius * mv_direction * MV_COORDINATE_SCALE_FPS)
                .clamp(TRUING_MV_MIN_FPS, TRUING_MV_MAX_FPS);
            let candidate_bc = (bc + radius * bc_direction * BC_COORDINATE_SCALE)
                .clamp(TRUING_BC_MIN, TRUING_BC_MAX);
            if candidate_mv == mv && candidate_bc == bc {
                continue;
            }
            let candidate = objective(request, model, candidate_mv, candidate_bc)?;
            evaluations += 1;
            if candidate < best {
                best = candidate;
                best_point = (candidate_mv, candidate_bc);
            }
        }

        let improvement = (current - best).max(0.0);
        if objective_improvement_is_material(current, best) {
            mv = best_point.0;
            bc = best_point.1;
            current = best;
            continue;
        }
        max_final_improvement = improvement;
        if radius <= MAP_OBJECTIVE_MIN_POLL_RADIUS {
            return Ok(ObjectivePollResult {
                mv,
                bc,
                converged: true,
                final_radius: radius,
                max_final_improvement,
                evaluations,
            });
        }
        radius = (radius * 0.25).max(MAP_OBJECTIVE_MIN_POLL_RADIUS);
    }

    Ok(ObjectivePollResult {
        mv,
        bc,
        converged: false,
        final_radius: radius,
        max_final_improvement,
        evaluations,
    })
}

fn objective(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
    mv: f64,
    bc: f64,
) -> Result<f64, UncertaintyTruingErrorV1> {
    let mut objective = 0.0;
    let ranges: Vec<f64> = request
        .observations
        .iter()
        .map(|observation| observation.range_yd)
        .collect();
    let predictions = model
        .predict_many_in_unit(mv, bc, &ranges, request.drop_unit)
        .map_err(forward_error)?;
    for (observation, prediction) in request.observations.iter().zip(predictions) {
        let prediction = prediction.ok_or_else(|| unreachable_range_error(observation.range_yd))?;
        let standardized = (prediction - observation.drop) / observation.sigma;
        objective += standardized * standardized;
    }
    if let Some(prior) = request.priors.muzzle_velocity_fps {
        let standardized = (mv - prior.mean) / prior.sigma;
        objective += standardized * standardized;
    }
    if let Some(prior) = request.priors.ballistic_coefficient {
        let standardized = (bc - prior.mean) / prior.sigma;
        objective += standardized * standardized;
    }
    if objective.is_finite() {
        Ok(objective)
    } else {
        Err(UncertaintyTruingErrorV1::ForwardModel(
            "non-finite penalized chi-square".to_string(),
        ))
    }
}

fn evaluate(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
    mv: f64,
    bc: f64,
) -> Result<Evaluation, UncertaintyTruingErrorV1> {
    let mut data_information = Symmetric2::default();
    let mut gradient = [0.0, 0.0];
    let mut results = Vec::with_capacity(request.observations.len());
    let mut chi_square = 0.0;

    // Fractional-sensitivity terms retain the physical parameter scaling even
    // though the optimizer itself works in fixed numerical coordinates.
    let (mut fractional_mv_norm2, mut fractional_bc_norm2): (f64, f64) = (0.0, 0.0);
    let ranges: Vec<f64> = request
        .observations
        .iter()
        .map(|observation| observation.range_yd)
        .collect();
    let rows =
        truing_jacobian_rows(model, mv, bc, &ranges, request.drop_unit).map_err(forward_error)?;
    for (observation, row) in request.observations.iter().zip(rows) {
        let row = row.ok_or_else(|| unreachable_range_error(observation.range_yd))?;
        let residual = row.predicted_drop - observation.drop;
        let standardized_residual = residual / observation.sigma;
        let j_mv = row.d_drop_d_mv * MV_COORDINATE_SCALE_FPS / observation.sigma;
        let j_bc = row.d_drop_d_bc * BC_COORDINATE_SCALE / observation.sigma;
        data_information.a00 += j_mv * j_mv;
        data_information.a01 += j_mv * j_bc;
        data_information.a11 += j_bc * j_bc;
        gradient[0] += j_mv * standardized_residual;
        gradient[1] += j_bc * standardized_residual;
        chi_square += standardized_residual * standardized_residual;
        fractional_mv_norm2 += (row.d_drop_d_mv * mv / observation.sigma).powi(2);
        fractional_bc_norm2 += (row.d_drop_d_bc * bc / observation.sigma).powi(2);
        results.push(WeightedTruingObservationResultV1 {
            range_yd: observation.range_yd,
            observed_drop: observation.drop,
            sigma: observation.sigma,
            predicted_drop: row.predicted_drop,
            residual,
            standardized_residual,
        });
    }

    let mut prior_information = Symmetric2::default();
    let mut prior_penalty = 0.0;
    if let Some(prior) = request.priors.muzzle_velocity_fps {
        let j = MV_COORDINATE_SCALE_FPS / prior.sigma;
        let standardized = (mv - prior.mean) / prior.sigma;
        prior_information.a00 += j * j;
        gradient[0] += j * standardized;
        prior_penalty += standardized * standardized;
    }
    if let Some(prior) = request.priors.ballistic_coefficient {
        let j = BC_COORDINATE_SCALE / prior.sigma;
        let standardized = (bc - prior.mean) / prior.sigma;
        prior_information.a11 += j * j;
        gradient[1] += j * standardized;
        prior_penalty += standardized * standardized;
    }
    let mut posterior_information = data_information;
    posterior_information.add_assign(prior_information);

    let bc_sensitivity_ratio = if fractional_mv_norm2 > 0.0 {
        (fractional_bc_norm2 / fractional_mv_norm2).sqrt()
    } else {
        0.0
    };
    let data_condition_number = column_condition(data_information);

    Ok(Evaluation {
        data_information,
        posterior_information,
        gradient,
        observation_results: results,
        chi_square,
        prior_penalty,
        bc_sensitivity_ratio,
        data_condition_number,
    })
}

fn column_condition(information: Symmetric2) -> f64 {
    if information.a00 <= 0.0 || information.a11 <= 0.0 {
        return f64::INFINITY;
    }
    let correlation = (information.a01 / (information.a00 * information.a11).sqrt())
        .clamp(-1.0, 1.0)
        .abs();
    if 1.0 - correlation <= 1.0e-15 {
        f64::INFINITY
    } else {
        (1.0 + correlation) / (1.0 - correlation)
    }
}

fn build_approximation(
    mv: f64,
    bc: f64,
    information: Symmetric2,
) -> Result<(TruingGaussianApproximationV1, Symmetric2), TruingApproximationFailureV1> {
    let (covariance_q, condition) = information
        .inverse_with_condition()
        .map_err(|code| TruingApproximationFailureV1 {
            code,
            message: match code {
                TruingApproximationFailureCodeV1::OptimizerDidNotConverge => {
                    "MAP optimizer did not converge".to_string()
                }
                TruingApproximationFailureCodeV1::MapAtParameterBound => {
                    "MAP is at a constrained parameter bound".to_string()
                }
                TruingApproximationFailureCodeV1::RankDeficientInformation => {
                    "likelihood-plus-prior information is rank deficient or numerically singular; collect more separated ranges or add an explicit prior".to_string()
                }
                TruingApproximationFailureCodeV1::NonFiniteInformation => {
                    "likelihood-plus-prior information or its inverse is non-finite".to_string()
                }
            },
        })?;
    let covariance = TruingCovarianceV1 {
        mv_variance_fps2: covariance_q.a00 * MV_COORDINATE_SCALE_FPS.powi(2),
        mv_bc_covariance_fps: covariance_q.a01 * MV_COORDINATE_SCALE_FPS * BC_COORDINATE_SCALE,
        bc_variance: covariance_q.a11 * BC_COORDINATE_SCALE.powi(2),
    };
    let mv_interval = GaussianIntervalV1::from_variance(mv, covariance.mv_variance_fps2)
        .ok_or_else(|| TruingApproximationFailureV1 {
            code: TruingApproximationFailureCodeV1::NonFiniteInformation,
            message: "MV marginal variance is invalid".to_string(),
        })?;
    let bc_interval =
        GaussianIntervalV1::from_variance(bc, covariance.bc_variance).ok_or_else(|| {
            TruingApproximationFailureV1 {
                code: TruingApproximationFailureCodeV1::NonFiniteInformation,
                message: "BC marginal variance is invalid".to_string(),
            }
        })?;
    let correlation = covariance.mv_bc_covariance_fps
        / (covariance.mv_variance_fps2 * covariance.bc_variance).sqrt();
    if !correlation.is_finite() {
        return Err(TruingApproximationFailureV1 {
            code: TruingApproximationFailureCodeV1::NonFiniteInformation,
            message: "MV/BC posterior correlation is non-finite".to_string(),
        });
    }
    Ok((
        TruingGaussianApproximationV1 {
            covariance,
            muzzle_velocity_interval_95: mv_interval,
            ballistic_coefficient_interval_95: bc_interval,
            mv_bc_correlation: correlation.clamp(-1.0, 1.0),
            scaled_information_condition_number: condition,
        },
        covariance_q,
    ))
}

fn build_predictive_bands(
    request: &UncertaintyTruingRequestV1,
    model: &TruingForwardModel<'_>,
    mv: f64,
    bc: f64,
    covariance_q: Option<Symmetric2>,
    warnings: &mut Vec<TruingUncertaintyWarningV1>,
) -> Result<Vec<TruingPredictiveBandV1>, UncertaintyTruingErrorV1> {
    let observed_min = request
        .observations
        .iter()
        .map(|observation| observation.range_yd)
        .fold(f64::INFINITY, f64::min);
    let observed_max = request
        .observations
        .iter()
        .map(|observation| observation.range_yd)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut warned_extrapolation = false;
    let mut bands = Vec::with_capacity(request.predictions.len());
    let prediction_ranges: Vec<f64> = request
        .predictions
        .iter()
        .map(|prediction| prediction.range_yd)
        .collect();
    let rows = truing_jacobian_rows(model, mv, bc, &prediction_ranges, request.drop_unit)
        .map_err(forward_error)?;
    for (prediction, row) in request.predictions.iter().zip(rows) {
        let row = row.ok_or_else(|| unreachable_range_error(prediction.range_yd))?;
        if !warned_extrapolation
            && (prediction.range_yd < observed_min || prediction.range_yd > observed_max)
        {
            push_warning(
                warnings,
                TruingUncertaintyWarningCodeV1::PredictionOutsideObservedDomain,
                "one or more predictive ranges lie outside the observed range domain; local linear uncertainty may understate nonlinear extrapolation risk".to_string(),
            );
            warned_extrapolation = true;
        }
        let latent_interval_95 = covariance_q.and_then(|covariance| {
            let g_mv = row.d_drop_d_mv * MV_COORDINATE_SCALE_FPS;
            let g_bc = row.d_drop_d_bc * BC_COORDINATE_SCALE;
            let variance = g_mv * g_mv * covariance.a00
                + 2.0 * g_mv * g_bc * covariance.a01
                + g_bc * g_bc * covariance.a11;
            // Tiny negative roundoff can occur when the covariance is highly
            // correlated; a material negative variance is treated as unavailable.
            let tolerance = 1.0e-12
                * (g_mv * g_mv * covariance.a00)
                    .abs()
                    .max((g_bc * g_bc * covariance.a11).abs())
                    .max(1.0);
            let variance = if variance >= 0.0 {
                variance
            } else if variance >= -tolerance {
                0.0
            } else {
                return None;
            };
            GaussianIntervalV1::from_variance(row.predicted_drop, variance)
        });
        let future_observation_interval_95 =
            match (latent_interval_95, prediction.future_observation_sigma) {
                (Some(latent), Some(sigma)) => GaussianIntervalV1::from_variance(
                    row.predicted_drop,
                    latent.standard_deviation.powi(2) + sigma.powi(2),
                ),
                _ => None,
            };
        bands.push(TruingPredictiveBandV1 {
            range_yd: prediction.range_yd,
            predicted_drop: row.predicted_drop,
            latent_interval_95,
            future_observation_interval_95,
        });
    }
    Ok(bands)
}

fn parameter_at_bound(mv: f64, bc: f64) -> bool {
    let mv_tolerance = 1.0e-6 * (TRUING_MV_MAX_FPS - TRUING_MV_MIN_FPS);
    let bc_tolerance = 1.0e-6 * (TRUING_BC_MAX - TRUING_BC_MIN);
    mv - TRUING_MV_MIN_FPS <= mv_tolerance
        || TRUING_MV_MAX_FPS - mv <= mv_tolerance
        || bc - TRUING_BC_MIN <= bc_tolerance
        || TRUING_BC_MAX - bc <= bc_tolerance
}

fn warn_prior_dominance(
    request: &UncertaintyTruingRequestV1,
    approximation: &TruingGaussianApproximationV1,
    warnings: &mut Vec<TruingUncertaintyWarningV1>,
) {
    if let Some(prior) = request.priors.muzzle_velocity_fps {
        if approximation.covariance.mv_variance_fps2 >= 0.8 * prior.sigma.powi(2) {
            push_warning(
                warnings,
                TruingUncertaintyWarningCodeV1::MvPriorDominated,
                "MV posterior width remains close to its explicit prior width; observations add little marginal MV information".to_string(),
            );
        }
    }
    if let Some(prior) = request.priors.ballistic_coefficient {
        if approximation.covariance.bc_variance >= 0.8 * prior.sigma.powi(2) {
            push_warning(
                warnings,
                TruingUncertaintyWarningCodeV1::BcPriorDominated,
                "BC posterior width remains close to its explicit prior width; observations add little marginal BC information".to_string(),
            );
        }
    }
}

fn push_warning(
    warnings: &mut Vec<TruingUncertaintyWarningV1>,
    code: TruingUncertaintyWarningCodeV1,
    message: String,
) {
    if warnings.iter().any(|warning| warning.code == code) {
        return;
    }
    warnings.push(TruingUncertaintyWarningV1 { code, message });
}

fn forward_error(error: Box<dyn std::error::Error>) -> UncertaintyTruingErrorV1 {
    UncertaintyTruingErrorV1::ForwardModel(error.to_string())
}

fn unreachable_range_error(range_yd: f64) -> UncertaintyTruingErrorV1 {
    UncertaintyTruingErrorV1::ForwardModel(format!(
        "trajectory did not reach requested range {range_yd:.3} yd"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn request_deserializes_from_the_wire_and_rejects_unknown_fields() {
        let json = serde_json::json!({
            "model": {
                "muzzle_velocity_fps": 2700.0, "ballistic_coefficient": 0.243,
                "drag_model": "g7", "mass_gr": 168.0, "diameter_in": 0.308,
                "zero_distance_yd": 100.0, "sight_height_in": 2.0,
                "temperature_f": 59.0, "pressure_inhg": 29.92,
                "humidity_pct": 50.0, "altitude_ft": 0.0
            },
            "drop_unit": "mil",
            "observations": [{"range_yd": 600.0, "drop": 4.2, "sigma": 0.1}],
            "priors": {"muzzle_velocity_fps": {"mean": 2700.0, "sigma": 15.0},
                       "ballistic_coefficient": null},
            "predictions": [{"range_yd": 800.0, "future_observation_sigma": null}]
        });
        let req: UncertaintyTruingRequestV1 =
            serde_json::from_value(json.clone()).expect("valid request deserializes");
        assert_eq!(req.observations.len(), 1);
        assert_eq!(req.model.muzzle_velocity_fps, 2700.0);

        // A misspelled field must be rejected, not silently dropped.
        let mut bad = json;
        bad["observations"][0]["sigmaa"] = serde_json::json!(0.2);
        assert!(serde_json::from_value::<UncertaintyTruingRequestV1>(bad).is_err());
    }
}
