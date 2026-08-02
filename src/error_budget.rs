//! MBA-1347: a per-input error budget and measurement-priority report.
//!
//! A shooter usually has time to improve exactly ONE input before a shot: rezero the sight, get
//! a better wind call, chronograph another string, and so on. This module propagates each
//! DECLARED per-input uncertainty (a caller-supplied one-sigma value for one [`InputAxis`]) to
//! impact covariance at the ranges of interest, via the same central-difference kernel the rest
//! of the 0.33.0 decision-support train is built on ([`central_difference`]), then ranks the
//! sources so the report can end with a concrete answer: which single input is worth improving
//! here, and which ones are not.
//!
//! Sources are preserved individually and NEVER collapsed into an "other" bucket -- that is this
//! report's whole reason to exist, as distinct from the existing WEZ (`monte-carlo --wez`)
//! attribution, which lumps everything not explicitly modelled into dispersion the caller cannot
//! attribute to a specific input at all. See [`error_budget`]'s doc comment for the full
//! contract.
//!
//! # Unavailable sources are recorded, never silently dropped
//!
//! [`central_difference`] can legitimately refuse to differentiate a declared axis:
//! [`KernelError::AxisUnsupportedForRequest`] (`Altitude` under a QNH-referenced atmosphere,
//! `ShotAzimuth` under compass-referenced wind), [`KernelError::AxisAbsent`] (a wind axis under
//! segmented wind), [`KernelError::CategoricalAxis`] (an effect toggle), or
//! [`KernelError::StepOutOfDomain`] (both perturbed sides, using the axis's OWN default step --
//! `error_budget` always requests [`central_difference`]'s `None`-step formula, never a custom
//! one, so the declared sigma has no bearing on whether this fires -- left the axis's physical
//! domain). Every one of these is recorded as an
//! [`UnavailableSourceV1`] (axis, declared sigma, a machine-readable
//! [`UnavailableReasonCodeV1`], and a human-readable reason) in
//! [`ErrorBudgetReportV1::unavailable_sources`], and the rest of the report is still produced
//! from whatever sources DID evaluate.
//!
//! **Silently skipping an unavailable axis would report "this input contributes no
//! uncertainty," the one wrong answer this ticket exists to prevent** -- a source that could not
//! be measured must never look identical to a source that WAS measured and found to contribute
//! exactly zero. [`SourceContributionV1`] is never constructed for an axis that failed; it only
//! ever describes an axis [`central_difference`] actually evaluated.
//!
//! Any OTHER error `central_difference` reports ([`KernelError::Solve`],
//! [`KernelError::Observation`], or the two defensive variants
//! [`KernelError::TypeMismatch`]/[`KernelError::NonFinite`]) is a genuine solver or trajectory
//! failure, not a normal "this input cannot be perturbed here" fact, and [`error_budget`]
//! propagates it unchanged rather than folding it into `unavailable_sources`. The classification
//! (see the private `unavailable_reason` below) is an exhaustive match with no wildcard arm, so
//! a future [`KernelError`] variant fails to compile here until it is explicitly placed in one
//! bucket or the other, rather than silently defaulting into whichever this match's last arm
//! happens to be. [`KernelError::DuplicateAxis`] is classified there too even though
//! `error_budget` itself constructs and returns it before that match ever runs (see "Sources are
//! validated up front" below) -- the exhaustiveness is over the whole `KernelError` type, not
//! only the subset `central_difference` can produce.
//!
//! One caller mistake is deliberately NOT laundered through this mechanism: a range in
//! `ranges_m` beyond `base.shot.max_range_m` would otherwise fail BOTH perturbed sides of EVERY
//! axis identically (an [`KernelError::Observation`] domain rejection on each side becomes
//! [`KernelError::StepOutOfDomain`]), recording every declared source as unavailable with a
//! misleading reason that blames the axis's own step when the real problem is that the caller
//! queried past the trajectory. See "Sources and ranges are validated up front" below --
//! `error_budget` rejects that case directly instead.
//!
//! # Sources and ranges are validated up front
//!
//! Before any [`central_difference`] call: every `range_m` in `ranges_m` must be finite and in
//! `[0, base.shot.max_range_m]`, or [`error_budget`] returns [`KernelError::Observation`]
//! immediately (an honest "this range cannot be observed on this trajectory," not a per-axis
//! "unavailable" that hides the real cause). Every declared `sigma` must be finite and
//! non-negative, and the same [`InputAxis`] must not appear twice in `sources` (two entries
//! would double-count that axis's variance and make its own leave-one-out counterfactual
//! ambiguous), or [`error_budget`] returns [`KernelError::NonFinite`] /
//! [`KernelError::DuplicateAxis`] respectively.
//!
//! # Ranking is deterministic
//!
//! Sources are ranked by [`SourceContributionV1::variance_share`], descending. Ties (equal
//! shares) break on a fixed, declaration-order-independent key (the axis's own `Debug` name), so
//! `error_budget(base, &[(A, sa), (B, sb)], ranges)` and
//! `error_budget(base, &[(B, sb), (A, sa)], ranges)` produce IDENTICAL orderings even when two
//! sources happen to contribute exactly the same share. A real-physics fixture essentially never
//! produces an exact tie in floating point, so a declaration-order test alone (the brief's own
//! `ranking_is_invariant_to_declaration_order`) would still pass with NO tie-break at all, as
//! long as Rust's sort remains stable and the two shares genuinely differ -- see
//! `tied_variance_shares_break_deterministically_regardless_of_input_order` in this module's
//! tests, which constructs a genuine tie directly against the sort function itself, and would
//! fail without the tie-break.
//!
//! # Cost
//!
//! One [`central_difference`] call per DECLARED source, covering every requested range at once
//! -- never once per source per range. `ranges_m` is passed through to the kernel unchanged and
//! the resulting `Vec<Derivative>` is indexed by range when building each row, so the cost of a
//! report is independent of how many ranges are requested and scales only with the number of
//! DECLARED sources.
//!
//! Per source, that is 2 real trajectory solves in the common (central-difference) case, or 3 if
//! one side fell outside the axis's physical domain and the kernel fell back to a one-sided
//! difference (see [`DifferenceScheme`]). If the axis `requires_rezero`
//! (`crate::perturbation::axis_meta(axis).requires_rezero`) and the request
//! carries a `zero_distance_m`, each of those solves is itself preceded by a fresh elevation
//! search of up to 60 trial solves (`find_zero_angle`, `src/cli_api.rs`) -- unavoidable, and not
//! something this module changes; see `crate::perturbation::derive`'s own module doc for where
//! that number comes from.
//!
//! Measured, not guessed (a previous task's cost doc on this branch understated its own number
//! by 5x before being corrected by direct measurement, so this number was obtained the same
//! way): a temporary instrumented low-level solve counter was added to `TrajectorySolver::solve`,
//! run once against this module's own three-source test fixture
//! (`every_declared_source_appears_individually`: `MuzzleVelocityMps` and
//! `BallisticCoefficient`, both `requires_rezero`, plus `WindSpeed`, which is not), then removed
//! (the working tree was diffed against the pre-instrumentation state afterward to confirm a
//! byte-for-byte revert of `src/cli_api.rs`). Measured results:
//!
//! - All three sources together at a single range: **79 low-level solves.**
//! - The IDENTICAL three sources requested over FOUR ranges instead of one: **79 again** --
//!   confirming the per-range independence above empirically, not just structurally.
//! - Decomposed by declaring each source alone (also at one range): `WindSpeed` (not
//!   `requires_rezero`) cost 2 (exactly the "2 solves, common case" above); `MuzzleVelocityMps`
//!   cost 37; `BallisticCoefficient` cost 40. `37 + 40 + 2 = 79`, matching the three-source
//!   total exactly -- confirming each declared source's cost is independent of, and additive
//!   with, every other declared source's, as the "one central_difference call per source" design
//!   above requires. The two `requires_rezero` axes each cost roughly 17-19 trial solves per
//!   perturbed side beyond the one real solve (`(37 - 2) / 2 = 17.5` average,
//!   `(40 - 2) / 2 = 19` average) -- well under the 60-iteration cap, not close to it, for this
//!   fixture's zero geometry.
//!
//! # Why the differencing step ignores the declared sigma
//!
//! `error_budget` always calls [`central_difference`] with `step: None` -- the axis's own small
//! default, never `Some(sigma)`. The delta method's Jacobian is the local SLOPE of impact at the
//! nominal point; the declared sigma only enters afterward, scaling that slope via
//! `J * Sigma * J^T`. The consequence: a large declared sigma is linearly extrapolated from a
//! slope measured over a much SMALLER window (`WindSpeed`'s own default step is 0.05 m/s,
//! regardless of whether the caller declared a 1 m/s or a 10 m/s wind-call sigma) -- this is
//! exactly the local-linearity limitation the `assumptions` payload already discloses, not a
//! separate concern. Using `Some(sigma)` instead was considered and rejected: it would push
//! `WindSpeed`/`RelativeHumidity`-style sigmas straight out of their physical domain on an
//! ordinary still-air or dry-air request (see the "One-sided fallback" section in
//! `crate::perturbation::derive`'s own module doc), making sources vanish into one-sided
//! fallbacks or `StepOutOfDomain` far more often -- the opposite of what a report whose whole
//! purpose is surfacing every declared source should do.

use serde::{Deserialize, Serialize};

use crate::perturbation::access::KernelError;
use crate::perturbation::derive::{central_difference, DifferenceScheme, Derivative};
use crate::perturbation::taxonomy::InputAxis;
use crate::solve_json::ResolvedSolveRequestV1;
use crate::special::normal_cdf;
use crate::trajectory_observation::TrajectoryObservationError;
use crate::truing_uncertainty::Symmetric2;

/// Schema version for [`ErrorBudgetReportV1`].
pub const ERROR_BUDGET_SCHEMA_VERSION_V1: u32 = 1;

/// The chi-square critical value for a 95% confidence region with 2 degrees of freedom -- the
/// same constant `crate::monte_carlo::calculate_confidence_ellipse` uses for its own,
/// sample-based ellipse.
const CHI2_95_2DOF: f64 = 5.991;

/// 95% confidence ellipse for a 2-dof (drop, windage) impact covariance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Ellipse95V1 {
    pub semi_major_m: f64,
    pub semi_minor_m: f64,
    /// Radians, measured from the drop axis toward the windage axis.
    pub rotation_rad: f64,
    pub area_m2: f64,
}

/// A target shape for [`p_hit_bivariate`] / [`error_budget_with_target`], always centred on the
/// nominal (zero-mean) impact point -- there is no separate "offset from point of aim" field, so
/// the reported hit probability implicitly assumes a well-zeroed rifle aimed at the target's own
/// centre. `width_m`/`height_m`/`radius_m` are clamped to non-negative internally by
/// [`p_hit_bivariate`]; a negative value is treated as zero rather than producing an inverted or
/// NaN result.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TargetGeometryV1 {
    /// Drop extent `height_m`, windage extent `width_m` -- matching this module's
    /// (drop, windage) axis order, not (x, y) or (width, height) screen convention.
    Rect { width_m: f64, height_m: f64 },
    Circle { radius_m: f64 },
}

/// One declared source's contribution to impact variance at one range.
///
/// Constructed only for a source [`central_difference`] actually evaluated -- an axis it refused
/// is recorded in [`UnavailableSourceV1`] instead, never here with a fabricated zero.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceContributionV1 {
    pub axis: InputAxis,
    /// The caller-declared one-sigma uncertainty for this axis, in the axis's own physical unit
    /// (`crate::perturbation::axis_meta(axis).kind`).
    pub sigma: f64,
    pub d_drop_d_x: f64,
    pub d_windage_d_x: f64,
    /// Which finite-difference scheme produced the two derivatives above. A one-sided scheme
    /// (`ForwardOneSided`/`BackwardOneSided`) has larger truncation error than `Central` -- see
    /// this module's "Unavailable sources" doc section and [`DifferenceScheme`]'s own doc for
    /// when the kernel falls back to one. Still air (`wind.speed_mps: 0.0`) is the ordinary,
    /// ROUTINE case for `WindSpeed`, not an exotic one -- wind-call uncertainty is this report's
    /// flagship use case, so a non-`Central` scheme here is expected, not a warning sign by
    /// itself.
    pub scheme: DifferenceScheme,
    /// Fraction of total impact variance (`sigma_drop_m^2 + sigma_windage_m^2`) attributable to
    /// this source. Sums to 1.0 across a row's `sources` whenever at least one source has a
    /// nonzero contribution.
    pub variance_share: f64,
    /// Reduction in the row's 95% ellipse area if this source alone were measured perfectly
    /// (its sigma set to zero, every other source unchanged). Always `>= 0.0`.
    ///
    /// **Degenerate with two or fewer sources.** The remaining covariance after removing any ONE
    /// of exactly two declared sources is an outer product (rank <= 1, area exactly 0 -- see
    /// `Symmetric2::largest_smallest_eigenvalues`'s doc for why that is normal, not a
    /// bug), so with exactly two sources BOTH report a reduction equal to the FULL ellipse area,
    /// even when one dominates the variance share overwhelmingly (e.g. a 99%/1% split). That is
    /// literally true (perfecting either one alone does leave a zero-area ellipse) but reads, next
    /// to a ranking, as a tie where none exists on `variance_share`. A declaration of two sources
    /// (e.g. muzzle velocity and wind call, the most likely pair this ticket sees) always has this
    /// property; prefer `variance_share` to distinguish sources when `sources.len() <= 2`. With
    /// three or more sources the reduction is generically discriminating.
    pub ellipse_area_reduction_m2: f64,
    /// The hit-probability gain over the row's target if THIS source alone were measured
    /// perfectly (its sigma set to zero, every other source unchanged), i.e.
    /// `p_hit(without this source) - p_hit(with every declared source)`. `Some` exactly when
    /// [`error_budget_with_target`] was given a target; `None` (never a fabricated `0.0`) when
    /// no target was supplied, matching [`ErrorBudgetRowV1::p_hit`].
    ///
    /// **Always `>= 0.0`** (see [`p_hit_bivariate`]'s doc: shrinking a target-centred impact
    /// covariance in the Loewner order cannot reduce the mass of a symmetric normal over a
    /// symmetric convex target -- Anderson's theorem). The unclamped value is checked against a
    /// small negative tolerance before being clamped to zero (a `debug_assert!` in
    /// [`error_budget_with_target`]), so a bug that made it SYSTEMATICALLY negative -- as
    /// opposed to a few ULP of quadrature noise -- fails a debug/test build loudly rather than
    /// being silently laundered into `0.0`.
    ///
    /// Shares the same two/fewer-sources degeneracy as `ellipse_area_reduction_m2` above (with
    /// exactly two declared sources, removing either one leaves the same rank-1 covariance
    /// shape), though the two fields are not numerically equal to each other -- one is an area in
    /// m^2, the other a probability in `[0, 1]`.
    pub p_hit_gain_if_perfect: Option<f64>,
}

/// Which structural refusal made a source unavailable -- the machine-readable counterpart to
/// [`UnavailableSourceV1::reason`]'s prose. Named identically to the [`KernelError`] variant it
/// comes from. Added at the same time as the rest of this `V1` type (not held back for a later
/// schema revision) specifically because adding a field to an already-shipped `V1` wire type
/// would be a breaking change; this crate had not shipped `error_budget` on any released version
/// when this field was added, so there is no such constraint yet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UnavailableReasonCodeV1 {
    AxisUnsupportedForRequest,
    AxisAbsent,
    CategoricalAxis,
    StepOutOfDomain,
}

/// A declared source [`error_budget`] could not evaluate for this request, and why.
///
/// See this module's "Unavailable sources" doc section. An axis appearing here never also
/// appears in any row's `sources` -- the two are disjoint by construction. This list is the same
/// for every row for the three refusals that depend only on `axis` and `base`'s OTHER fields
/// (`code` other than `StepOutOfDomain`): those never depend on which ranges were requested.
/// `StepOutOfDomain` specifically COULD in principle depend on range (a query near a perturbed
/// request's own shrunk domain -- see `crate::perturbation::derive`'s
/// `target_distance_falls_back_to_one_sided_when_queried_near_its_own_max_range`), but
/// `error_budget` rejects the common, obvious version of that (a range beyond the BASE request's
/// own `max_range_m`) up front instead of letting it reach this list at all -- see
/// `error_budget`'s "Sources and ranges are validated up front" doc section.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnavailableSourceV1 {
    pub axis: InputAxis,
    pub sigma: f64,
    pub code: UnavailableReasonCodeV1,
    pub reason: String,
}

/// The impact covariance and ranked sources at one requested range.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ErrorBudgetRowV1 {
    pub range_m: f64,
    pub sigma_drop_m: f64,
    pub sigma_windage_m: f64,
    pub covariance_m2: f64,
    pub ellipse_95: Ellipse95V1,
    /// Probability the impact (drop, windage) falls inside the row's target, computed by
    /// [`p_hit_bivariate`] from this row's own `sigma_drop_m`/`sigma_windage_m`/`covariance_m2`.
    /// `Some` exactly when [`error_budget_with_target`] was given a target; `None` (never a
    /// fabricated number) when no target was supplied -- see that function's "Hit probability"
    /// doc section.
    pub p_hit: Option<f64>,
    /// Ranked by [`SourceContributionV1::variance_share`], descending, most-informative first --
    /// see this module's "Ranking is deterministic" doc section.
    pub sources: Vec<SourceContributionV1>,
    /// A plain-language statement of which single input is most worth improving at this range,
    /// or that none of the declared sources contributes (or that none could be evaluated at
    /// all -- worded differently; see [`error_budget`]).
    pub priority_statement: String,
}

/// Per-input uncertainty propagation and measurement-priority report (MBA-1347).
///
/// Carries [`method`](ErrorBudgetReportV1#structfield.method) and
/// [`assumptions`](ErrorBudgetReportV1#structfield.assumptions) in the payload itself, not only
/// in prose documentation -- see [`error_budget`]'s "Honesty" doc section and the
/// `the_report_declares_independence_and_linearity` test.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ErrorBudgetReportV1 {
    pub schema_version: u32,
    pub method: String,
    pub assumptions: Vec<String>,
    /// Declared sources this report could not evaluate. See this module's "Unavailable sources"
    /// doc section -- never silently omitted.
    pub unavailable_sources: Vec<UnavailableSourceV1>,
    pub rows: Vec<ErrorBudgetRowV1>,
}

/// One declared source's central-difference Jacobian across every requested range, computed
/// exactly once -- see this module's "Cost" doc section.
struct AxisJacobian {
    axis: InputAxis,
    sigma: f64,
    /// One entry per `ranges_m`, same order -- `central_difference`'s own contract.
    derivatives: Vec<Derivative>,
}

/// One axis's contribution at a single range, gathered from an `AxisJacobian` for row-building.
#[derive(Debug, Clone, Copy)]
struct Entry {
    axis: InputAxis,
    sigma: f64,
    scheme: DifferenceScheme,
    d_drop_d_x: f64,
    d_windage_d_x: f64,
}

/// Classify a `central_difference` failure as either "this source is unavailable, but the rest
/// of the report should still be produced" (`Some((code, reason))`) or "this is a genuine
/// failure the caller must see" (`None`).
///
/// `Some` for the four structural refusals this module's doc names explicitly
/// (`AxisUnsupportedForRequest`, `AxisAbsent`, `CategoricalAxis`, `StepOutOfDomain`); `None` for
/// everything else (`Solve`, `Observation`, `TypeMismatch`, `NonFinite`, `DuplicateAxis`), which
/// the caller must propagate. Exhaustive over every `KernelError` variant with no wildcard arm,
/// so a future variant fails to compile here until it is explicitly placed in one bucket or the
/// other -- this is what forced `DuplicateAxis` to be classified here even though
/// `error_budget`'s own up-front validation constructs and returns it directly, never through
/// `central_difference` at all.
fn unavailable_reason(e: &KernelError) -> Option<(UnavailableReasonCodeV1, String)> {
    match e {
        KernelError::AxisUnsupportedForRequest { reason, .. } => {
            Some((UnavailableReasonCodeV1::AxisUnsupportedForRequest, reason.to_string()))
        }
        KernelError::AxisAbsent(_) => Some((
            UnavailableReasonCodeV1::AxisAbsent,
            "this axis has no single scalar value on this request (for the three wind axes, \
             this means the wind is declared as a segmented profile rather than a constant \
             speed/direction)"
                .to_string(),
        )),
        KernelError::CategoricalAxis(_) => Some((
            UnavailableReasonCodeV1::CategoricalAxis,
            "this axis is categorical (a toggle or an enumerated choice), not a continuous \
             quantity, and cannot be assigned a one-sigma uncertainty or differentiated"
                .to_string(),
        )),
        KernelError::StepOutOfDomain { attempted, .. } => Some((
            UnavailableReasonCodeV1::StepOutOfDomain,
            format!(
                "central differencing needs to perturb this axis by {attempted:.6} in both \
                 directions from its nominal value, using the axis's own default step (this does \
                 not depend on the declared sigma), and both directions left its physical \
                 domain; its sensitivity could not be measured at this operating point"
            ),
        )),
        KernelError::Solve { .. }
        | KernelError::Observation(_)
        | KernelError::TypeMismatch(_)
        | KernelError::NonFinite(_)
        | KernelError::DuplicateAxis(_) => None,
    }
}

/// Sum of `sigma^2 * [d_drop, d_windage] * [d_drop, d_windage]^T` over `entries`, as a
/// [`Symmetric2`] (`a00`/`a11` are drop/windage variance, `a01` is their covariance) --
/// optionally skipping one axis entirely (equivalent to, but cheaper than, substituting a zero
/// sigma for it), used to compute the "if this source were measured perfectly" ellipse.
fn accumulate(entries: &[Entry], exclude: Option<InputAxis>) -> Symmetric2 {
    let mut total = Symmetric2::default();
    for e in entries {
        if Some(e.axis) == exclude {
            continue;
        }
        let s2 = e.sigma * e.sigma;
        total.add_assign(Symmetric2 {
            a00: e.d_drop_d_x * e.d_drop_d_x * s2,
            a01: e.d_drop_d_x * e.d_windage_d_x * s2,
            a11: e.d_windage_d_x * e.d_windage_d_x * s2,
        });
    }
    total
}

/// The 95% ellipse for an impact covariance `cov` (`a00`/`a11` = drop/windage variance, `a01` =
/// their covariance). Never fails -- see `Symmetric2::largest_smallest_eigenvalues`.
fn ellipse_95(cov: Symmetric2) -> Ellipse95V1 {
    let (largest, smallest) = cov.largest_smallest_eigenvalues();
    let semi_major_m = (CHI2_95_2DOF * largest).sqrt();
    let semi_minor_m = (CHI2_95_2DOF * smallest).sqrt();
    // Standard closed form for a symmetric 2x2 matrix's eigenvector angle -- equivalent to (and
    // more robust at cov.a01 == 0 than) solving tan(theta) = (largest - a00) / a01 directly.
    let rotation_rad = 0.5 * (2.0 * cov.a01).atan2(cov.a00 - cov.a11);
    Ellipse95V1 {
        semi_major_m,
        semi_minor_m,
        rotation_rad,
        area_m2: std::f64::consts::PI * semi_major_m * semi_minor_m,
    }
}

/// Rank by [`SourceContributionV1::variance_share`], descending, breaking ties on the axis's own
/// `Debug` name so the order never depends on how the caller declared `sources` -- see this
/// module's "Ranking is deterministic" doc section.
///
/// Uses `f64::total_cmp` rather than `partial_cmp` -- `variance_share` is finite and
/// non-negative for any ordinarily-sized declared sigma, but an astronomically large (still
/// finite) sigma can overflow `sigma * sigma` to infinity and produce a NaN share
/// (`inf / inf`); this crate ships a fuzz suite, and a library function must not panic on an
/// extreme-but-technically-valid input. `total_cmp` gives NaN a well-defined (if not physically
/// meaningful) place in the order instead.
fn sort_by_variance_share_desc(sources: &mut [SourceContributionV1]) {
    sources.sort_by(|x, y| {
        y.variance_share
            .total_cmp(&x.variance_share)
            .then_with(|| format!("{:?}", x.axis).cmp(&format!("{:?}", y.axis)))
    });
}

/// Compose the row's plain-language measurement-priority statement.
fn build_priority_statement(
    sources: &[SourceContributionV1],
    unavailable: &[UnavailableSourceV1],
    range_m: f64,
) -> String {
    let unavailable_note = || -> String {
        if unavailable.is_empty() {
            String::new()
        } else {
            format!(
                " {} declared source{} could not be evaluated at all -- see \
                 unavailable_sources.",
                unavailable.len(),
                if unavailable.len() == 1 { "" } else { "s" },
            )
        }
    };
    match sources.first() {
        Some(top) if top.variance_share > 0.0 => {
            let caveat = if top.scheme != DifferenceScheme::Central {
                " (from a one-sided approximation, not a central difference -- see assumptions)"
            } else {
                ""
            };
            format!(
                "{:?} dominates at {range_m:.0} m ({:.1}% of impact variance){caveat}. \
                 Measuring it better is the highest-value single improvement here.{}",
                top.axis,
                top.variance_share * 100.0,
                unavailable_note(),
            )
        }
        Some(_) => format!(
            "No declared source contributes uncertainty at {range_m:.0} m (every evaluated \
             sigma is zero).{}",
            unavailable_note(),
        ),
        None if !unavailable.is_empty() => format!(
            "None of the declared sources could be evaluated at {range_m:.0} m -- see \
             unavailable_sources for why."
        ),
        None => format!("No sources were declared for this report at {range_m:.0} m."),
    }
}

// ---------------------------------------------------------------------------------------------
// Hit probability (MBA-1347, Task 11): the mass of a bivariate normal impact distribution over a
// target, and the value-of-information question built on it -- see `p_hit_bivariate`'s own doc
// comment for the math and `error_budget_with_target`'s "Hit probability" doc section for how a
// report uses it.
// ---------------------------------------------------------------------------------------------

/// 20-node Gauss-Legendre quadrature on `[-1, 1]`: the positive half of the standard abscissas
/// and weights (the negative half is generated by symmetry in [`gauss_legendre_20`]). Fixed
/// order, so [`p_hit_bivariate`]'s result is deterministic and reproducible across platforms --
/// the order (and the panel-splitting built on it below) is named in
/// [`ErrorBudgetReportV1::method`] whenever a target is supplied (the `"_gl20_pm6sigma"` suffix).
const GL20_X: [f64; 10] = [
    0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271,
    0.636_053_680_726_515, 0.7463319064601508, 0.8391169718222188, 0.912_234_428_251_326,
    0.9639719272779138, 0.9931285991850949,
];

/// Weights matching [`GL20_X`], same order.
const GL20_W: [f64; 10] = [
    0.1527533871307258, 0.1491729864726037, 0.142_096_109_318_382, 0.1316886384491766,
    0.1181945319615184, 0.1019301198172404, 0.0832767415767048, 0.0626720483341091,
    0.0406014298003869, 0.0176140071391521,
];

/// Integrate `f` over `[lo, hi]` with the fixed 20-node Gauss-Legendre rule above.
///
/// [`p_hit_bivariate`] calls this once per SMOOTH panel of its integration domain, never once
/// over the whole domain in a single call -- see that function's doc comment for why a single
/// panel is not accurate enough here. Returns `0.0` for an empty or inverted interval (`hi <=
/// lo`), which the caller legitimately produces whenever two panel boundaries coincide (e.g. a
/// correlation-crossing point that lands exactly on the target's own edge).
fn gauss_legendre_20(lo: f64, hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    if hi <= lo {
        return 0.0;
    }
    let mid = 0.5 * (hi + lo);
    let half = 0.5 * (hi - lo);
    let mut acc = 0.0;
    for k in 0..10 {
        for sign in [-1.0f64, 1.0f64] {
            let u = mid + sign * half * GL20_X[k];
            acc += GL20_W[k] * half * f(u);
        }
    }
    acc
}

/// The windage interval `[a, b]` `target` admits at drop-offset `u`: constant for a rectangle,
/// the chord of the circle at that height for a circle. Total and well-defined (never NaN) for
/// any `u`, including `|u|` beyond a circle's radius (returns `(0.0, 0.0)`, an empty interval,
/// rather than requiring the caller to pre-filter).
fn windage_bounds_at(u: f64, target: TargetGeometryV1) -> (f64, f64) {
    match target {
        TargetGeometryV1::Rect { width_m, .. } => {
            let half_w = width_m.max(0.0) / 2.0;
            (-half_w, half_w)
        }
        TargetGeometryV1::Circle { radius_m } => {
            let r = radius_m.max(0.0);
            let x = (r * r - u * u).max(0.0).sqrt();
            (-x, x)
        }
    }
}

/// P(impact falls inside `target`) for a bivariate normal impact distribution centred at the
/// origin -- the nominal (zero-mean) trajectory solution -- with drop variance `var_drop`,
/// windage variance `var_wind`, and drop/windage covariance `cov`. `target` is always centred on
/// that same origin; see [`TargetGeometryV1`]'s doc for why there is no separate aim-point
/// offset.
///
/// # The math (pinned, MBA-1347 spec section 6.2)
///
/// A correlated covariance does NOT let the rectangle probability separate into a product of two
/// [`normal_cdf`] differences (`uncorrelated_rectangle_matches_the_separable_closed_form` in this
/// module's tests exists specifically to demonstrate the separable form is only valid at zero
/// correlation, and `a_strongly_correlated_case_differs_materially_from_the_wrong_separable_approximation`
/// shows how far a real, moderately-correlated case departs from it). Instead this integrates
/// over the drop axis and, at each drop value `u`, applies the CONDITIONAL normal distribution of
/// windage given that drop:
///
/// `P = integral of phi(u) * [Phi(beta(u)) - Phi(alpha(u))] du`,
///
/// where `phi` is the drop marginal's density, `Phi` is [`normal_cdf`], and `alpha(u)`/`beta(u)`
/// are the target's windage bounds at drop `u` (constant for a rectangle; the circle's chord for
/// a circle), expressed in units of the conditional windage standard deviation and offset by the
/// conditional mean `rho * (sigma_windage / sigma_drop) * u`.
///
/// # Two degenerate covariances, handled explicitly (not merely "does not panic")
///
/// - **Zero total variance** (`var_drop <= 0.0 && var_wind <= 0.0`): a deterministic impact
///   exactly at the origin, which is always the target's own centre here -- inside by
///   definition, so this returns `1.0` outright without touching the quadrature below.
/// - **Zero drop variance alone** (`var_drop <= 0.0`, `var_wind > 0.0`): drop is deterministic at
///   0 but windage is not. The quadrature above integrates OVER the drop axis, which cannot
///   represent a Dirac delta; instead of letting `sd -> 0` silently zero out every quadrature
///   node's density (which the naive translation of the formula above does, and which would
///   wrongly report `0.0` regardless of target size -- caught by
///   `drop_deterministic_windage_random_matches_closed_form_not_hardcoded_zero` in this module's
///   tests), this evaluates the windage marginal directly at drop `= 0`.
///
/// A single declared source (one nonzero-sigma axis) produces a RANK-1 covariance -- `cov`
/// exactly `+-sigma_drop * sigma_windage` before the correlation clamp below -- which is the
/// routine case [`error_budget_with_target`] hits every time it prices "if this source alone
/// were perfected" against a row with exactly one OTHER remaining source (see
/// `SourceContributionV1::p_hit_gain_if_perfect`'s doc). `rho` is clamped to
/// `[-0.999_999, 0.999_999]` so the conditional variance below is never exactly zero, avoiding a
/// division by zero without needing a third special case for perfect correlation.
///
/// # Why the quadrature is PANELLED, not a single 20-node call over the whole domain
///
/// The spec pins "fixed-order Gauss-Legendre quadrature over a truncated +/-6 sigma domain," but
/// a single 20-node rule spread across the WHOLE `[-6 sigma_drop, 6 sigma_drop]` interval is not
/// merely imprecise, it is badly wrong for realistic inputs, for two DIFFERENT reasons, both
/// found by comparing that naive translation against an independent fine (4000-point-per-smooth-
/// piece composite Simpson) reference over a broad sweep of target sizes and correlations, not
/// guessed:
///
/// **First: the target boundary is a genuine discontinuity in the naive formulation.** Written
/// as "integrate over the whole +/-6 sigma domain, contributing zero outside the target's own
/// drop extent," the integrand jumps from a generic nonzero value to zero exactly at the
/// target's edge. Gauss-Legendre quadrature assumes smoothness across its whole panel; a hidden
/// jump degrades it to first-order accuracy. Measured on this module's own pinned
/// zero-correlation rectangle test (drop/windage sigma 0.10 m/0.20 m, target 0.30 m x 0.40 m):
/// the naive single-panel translation is wrong by **8.6e-3** against the closed form (the test
/// requires `< 1e-6`) -- for a SMALL circular target relative to sigma (radius one-fifth of
/// sigma: radius 0.02 with sigma 0.1 in both axes) the naive version places every one of its 20
/// nodes outside the target entirely and returns exactly **0.0** against a true value near
/// 0.020.
///
/// Fix: restrict the integration domain to the target's own drop extent intersected with the
/// +/-6 sigma truncation (`[lo, hi]` below) -- outside that range the contribution is EXACTLY
/// zero, not approximately zero, so there is nothing to lose by not integrating there at all.
/// This alone brings the zero-correlation rectangle case to ~1e-16 (machine precision, since the
/// integrand reduces to a constant windage factor times a plain Gaussian bump, which 20-point
/// Gauss-Legendre integrates essentially exactly).
///
/// **Second: near-perfect correlation creates a separate, INTERNAL sharp transition the target
/// boundary fix does not touch.** As `|rho| -> 1`, the conditional windage standard deviation
/// `sigma_w * sqrt(1 - rho^2) -> 0`, so the bracketed `[Phi(beta(u)) - Phi(alpha(u))]` factor
/// above becomes an increasingly steep (though still, short of exactly `rho = +-1`, smooth)
/// sigmoid in `u`, centred wherever the conditional mean crosses the window bound -- a location
/// that can fall anywhere inside the domain, not just at its edges. This is not a rare input:
/// EVERY "if this source alone were perfected" comparison in [`error_budget_with_target`]
/// evaluates a covariance with exactly one remaining source, which is exactly rank-1 (`rho` at
/// the +-0.999_999 clamp). Measured on that exact shape (a real single-source covariance,
/// `sigma_drop` = 18.5, `sigma_windage` = 6.0, `rho` clamped to -0.999_999, swept over target
/// heights/widths from 0.05x to 20x each sigma): the boundary-restricted-but-still-single-panel
/// quadrature is wrong by up to **0.28** against the fine reference at the swept extremes
/// (height 20x sigma_drop, width 1x sigma_windage), and by up to **7.4e-2** even restricted to
/// height/width within 0.5x-4x of the natural (sigma_drop, sigma_windage) scale -- nowhere near a
/// contrived corner.
///
/// Fix: when `|rho|` is non-negligible, ALSO split the domain at the (closed-form) drop value(s)
/// where the conditional mean crosses the window's bound -- `+-half_width * sigma_drop / (rho *
/// sigma_windage)` for a rectangle (the window bound is constant), or `+-radius / sqrt(1 + (rho *
/// sigma_windage / sigma_drop)^2)` for a circle (from solving the circle's own chord equation) --
/// giving the sigmoid its own smooth sub-panel instead of sharing one with the flat shoulder on
/// either side. This brings the worst case measured over a broad synthetic stress sweep (several
/// (sigma_drop, sigma_windage) magnitudes including the realistic 18.5/6.0 pair above, target
/// heights/widths from 0.05x to 20x sigma, `|rho|` up to 0.999_999, both shapes) down under
/// **1e-3** (6.1e-4 for rectangles, 2.3e-4 for circles) -- and realistic (roughly comparable
/// width/height, moderate correlation) target shapes measured one to two further orders of
/// magnitude better than that worst case.
///
/// This is a more careful IMPLEMENTATION of the pinned formula, not a different one: the number
/// of panel BOUNDARIES is bounded at compile time (at most 4: the two domain edges -- the
/// target's own edge only ever contributes 0 extra boundaries since it already bounds `[lo,
/// hi]` -- plus up to two correlation-crossing points), so there are at most 3 panels, each
/// integrated by the exact same 20-node rule named in the spec. The result therefore stays
/// deterministic and its cost stays bounded (at most 3 panels * 20 nodes * 2 [`normal_cdf`]
/// calls per node = 120 evaluations of [`normal_cdf`], negligible next to the real trajectory
/// solves [`error_budget_with_target`] needs to build the covariance in the first place).
///
/// # Bounded and monotone
///
/// Always clamped to `[0.0, 1.0]` before returning. Monotone non-decreasing in target size for a
/// FIXED covariance (a bigger rectangle or circle strictly contains a smaller one centred at the
/// same origin, so the region of integration only grows) -- see
/// `p_hit_is_bounded_and_grows_with_target_size` (circle, verbatim from the spec) and
/// `p_hit_grows_with_target_size_for_a_rectangle_too` in this module's tests.
pub fn p_hit_bivariate(var_drop: f64, var_wind: f64, cov: f64, target: TargetGeometryV1) -> f64 {
    let sd = var_drop.max(0.0).sqrt();
    let sw = var_wind.max(0.0).sqrt();

    // Zero total variance: see doc comment above.
    if sd <= 0.0 && sw <= 0.0 {
        return 1.0;
    }

    // Zero drop variance alone: see doc comment above -- the windage marginal evaluated directly
    // at drop = 0, not routed through a quadrature that cannot represent a Dirac delta.
    if sd <= 0.0 {
        let (a, b) = windage_bounds_at(0.0, target);
        return (normal_cdf(b / sw) - normal_cdf(a / sw)).clamp(0.0, 1.0);
    }

    let rho = if sw > 0.0 { (cov / (sd * sw)).clamp(-0.999_999, 0.999_999) } else { 0.0 };
    let cond_sw = sw * (1.0 - rho * rho).max(0.0).sqrt();

    // Restrict the domain to the target's own drop extent intersected with the +/-6 sigma
    // truncation -- outside it the contribution is exactly zero, not approximately zero (see
    // doc comment's "First: the target boundary..." section above).
    let edge = match target {
        TargetGeometryV1::Rect { height_m, .. } => height_m.max(0.0) / 2.0,
        TargetGeometryV1::Circle { radius_m } => radius_m.max(0.0),
    };
    let hi = (6.0 * sd).min(edge);
    if hi <= 0.0 {
        return 0.0;
    }
    let lo = -hi;

    // Extra panel boundaries at the conditional-mean/window-bound crossing point(s), needed only
    // when |rho| is non-negligible -- see doc comment's "Second: near-perfect correlation..."
    // section above.
    let mut panel_bounds = vec![lo, hi];
    if rho.abs() > 1e-9 {
        let k = rho * (sw / sd); // conditional mean of windage given drop is k * u
        let c = match target {
            TargetGeometryV1::Rect { width_m, .. } => (width_m.max(0.0) / 2.0) / k,
            TargetGeometryV1::Circle { radius_m } => radius_m.max(0.0) / (1.0 + k * k).sqrt(),
        };
        for candidate in [c, -c] {
            if candidate.is_finite() && candidate > lo && candidate < hi {
                panel_bounds.push(candidate);
            }
        }
    }
    panel_bounds.sort_by(f64::total_cmp);
    panel_bounds.dedup();

    let sqrt_2pi = (2.0 * std::f64::consts::PI).sqrt();
    let mut acc = 0.0;
    for w in panel_bounds.windows(2) {
        acc += gauss_legendre_20(w[0], w[1], |u| {
            let density = (-0.5 * (u / sd) * (u / sd)).exp() / (sd * sqrt_2pi);
            let (a, b) = windage_bounds_at(u, target);
            let mean = rho * (sw / sd) * u;
            let p = if cond_sw > 0.0 {
                normal_cdf((b - mean) / cond_sw) - normal_cdf((a - mean) / cond_sw)
            } else if mean >= a && mean <= b {
                1.0
            } else {
                0.0
            };
            density * p
        });
    }
    acc.clamp(0.0, 1.0)
}

/// Propagate each declared per-input uncertainty in `sources` to impact covariance at every
/// range in `ranges_m`, via central differences through the real solver
/// ([`central_difference`]), and rank the sources by their share of impact variance.
///
/// `sources` is `(axis, sigma)` pairs: `sigma` is the caller's one-sigma uncertainty for that
/// axis, in the axis's own physical unit (`crate::perturbation::axis_meta(axis).kind`). Every
/// `sigma` must be finite and non-negative and every axis must appear at most once -- both are
/// validated up front, before any solve; see "Sources and ranges are validated up front" and
/// `# Errors` below. (Declaring the same axis twice would double-count its variance and make its
/// own leave-one-out counterfactual ambiguous -- which of the two entries would "removing this
/// source" mean? -- so it is rejected rather than given an arbitrary answer.)
///
/// # Honesty
///
/// [`ErrorBudgetReportV1::method`] and [`ErrorBudgetReportV1::assumptions`] state, in the
/// payload itself (not only in prose documentation), that: sources are treated as INDEPENDENT
/// (no correlation between them is modelled); propagation is FIRST-ORDER/local-linear about the
/// nominal solution, using the axis's own small default differencing step regardless of the
/// declared sigma (not exact for large or non-Gaussian input uncertainty -- see "Why the
/// differencing step ignores the declared sigma" above); the 95% ellipse assumes an
/// approximately Gaussian impact distribution; a source's derivative may be one-sided rather
/// than central; and an unavailable source is not the same fact as a zero-contribution one. See
/// `the_report_declares_independence_and_linearity` in this module's tests.
///
/// # Unavailable sources
///
/// See this module's top-level "Unavailable sources" doc section.
///
/// # Errors
///
/// - [`KernelError::Observation`] immediately, before any solve, if any `range_m` in `ranges_m`
///   is not finite or falls outside `[0, base.shot.max_range_m]` -- see "Sources and ranges are
///   validated up front" above.
/// - [`KernelError::NonFinite`] immediately if any declared `sigma` is not finite or is
///   negative.
/// - [`KernelError::DuplicateAxis`] immediately if the same axis appears more than once in
///   `sources`.
/// - Otherwise, propagates any [`KernelError`] from [`central_difference`] that is not one of
///   the four structural refusals recorded in [`ErrorBudgetReportV1::unavailable_sources`]
///   instead -- a genuine solver or trajectory failure, not a normal "this input cannot be
///   perturbed here" fact.
///
/// A thin wrapper over [`error_budget_with_target`] passing `None` -- every row's `p_hit` and
/// every source's `p_hit_gain_if_perfect` come back `None` (never a fabricated number), and
/// `method`/`assumptions` say nothing about hit probability. Call
/// [`error_budget_with_target`] directly to also get those.
pub fn error_budget(
    base: &ResolvedSolveRequestV1,
    sources: &[(InputAxis, f64)],
    ranges_m: &[f64],
) -> Result<ErrorBudgetReportV1, KernelError> {
    error_budget_with_target(base, sources, ranges_m, None)
}

/// As [`error_budget`], but additionally reports hit probability over `target` when it is
/// `Some`: each row's [`ErrorBudgetRowV1::p_hit`], and each of its sources'
/// [`SourceContributionV1::p_hit_gain_if_perfect`] -- the hit-probability gain if that source
/// alone were measured perfectly, the value-of-information number this ticket exists to answer.
/// [`error_budget`] is a thin wrapper passing `None`. Validation, ranking, unavailable-source
/// handling, and cost are otherwise IDENTICAL to [`error_budget`] and documented on it and this
/// module's top-level doc comment; this doc comment covers only what `target` adds.
///
/// # Hit probability
///
/// When `target` is `Some`, each row's `p_hit` is [`p_hit_bivariate`] evaluated at that row's own
/// impact covariance (`sigma_drop_m`, `sigma_windage_m`, `covariance_m2` -- the same numbers the
/// row already reports, not a separately recomputed covariance). Each source's
/// `p_hit_gain_if_perfect` is the SAME row's `p_hit` if that one source's sigma were zero instead
/// (every other declared source unchanged) minus the row's actual `p_hit`, clamped to `>= 0.0`
/// after a `debug_assert!` that the unclamped value is not more than a small tolerance below
/// zero -- see [`SourceContributionV1::p_hit_gain_if_perfect`]'s own doc for why a systematically
/// negative raw value must fail loudly rather than be silently clamped away.
/// [`ErrorBudgetReportV1::method`] gains a `"_gl20_pm6sigma"` suffix and
/// [`ErrorBudgetReportV1::assumptions`] gains one more entry naming the quadrature and the
/// target-centred-on-aim-point assumption -- see
/// `the_report_names_the_quadrature_and_the_aim_point_assumption_only_when_a_target_is_supplied`
/// in this module's tests.
pub fn error_budget_with_target(
    base: &ResolvedSolveRequestV1,
    sources: &[(InputAxis, f64)],
    ranges_m: &[f64],
    target: Option<TargetGeometryV1>,
) -> Result<ErrorBudgetReportV1, KernelError> {
    // Validate ranges_m up front (review I3): an out-of-range query would otherwise fail BOTH
    // perturbed sides of EVERY axis identically (a genuine Observation domain rejection on each
    // side collapses to StepOutOfDomain), recording every declared source as unavailable with a
    // reason that blames the axis's own step -- laundering a caller mistake (a range beyond the
    // trajectory) into a plausible-looking per-axis explanation. Reject it directly instead, the
    // same way `evaluate`'s own `observation_at_range_checked` would once a solve actually ran.
    for &range_m in ranges_m {
        if !range_m.is_finite() {
            return Err(KernelError::Observation(TrajectoryObservationError::NonFiniteQuery {
                distance_m: range_m,
            }));
        }
        if range_m < 0.0 || range_m > base.shot.max_range_m {
            return Err(KernelError::Observation(TrajectoryObservationError::OutOfRange {
                requested_m: range_m,
                minimum_m: 0.0,
                maximum_m: base.shot.max_range_m,
            }));
        }
    }

    // Validate sources up front (review I4 + pre-existing sigma check): every sigma finite and
    // non-negative, and no axis declared twice. Checked together, before any solve.
    for (i, &(axis, sigma)) in sources.iter().enumerate() {
        if !(sigma.is_finite() && sigma >= 0.0) {
            return Err(KernelError::NonFinite(axis));
        }
        if sources[..i].iter().any(|&(earlier_axis, _)| earlier_axis == axis) {
            return Err(KernelError::DuplicateAxis(axis));
        }
    }

    let mut jac: Vec<AxisJacobian> = Vec::with_capacity(sources.len());
    let mut unavailable: Vec<UnavailableSourceV1> = Vec::new();

    for &(axis, sigma) in sources {
        // One central_difference call per DECLARED source, covering every range in ranges_m at
        // once -- never re-derived per range. See this module's "Cost" doc section.
        match central_difference(base, axis, ranges_m, None) {
            Ok(derivatives) => jac.push(AxisJacobian { axis, sigma, derivatives }),
            Err(e) => match unavailable_reason(&e) {
                Some((code, reason)) => {
                    unavailable.push(UnavailableSourceV1 { axis, sigma, code, reason })
                }
                None => return Err(e),
            },
        }
    }

    let mut rows = Vec::with_capacity(ranges_m.len());
    for (i, &range_m) in ranges_m.iter().enumerate() {
        let entries: Vec<Entry> = jac
            .iter()
            .map(|j| {
                let d = &j.derivatives[i];
                debug_assert_eq!(
                    d.range_m, range_m,
                    "central_difference's Nth derivative must be tagged with ranges_m's Nth range"
                );
                Entry {
                    axis: j.axis,
                    sigma: j.sigma,
                    scheme: d.scheme,
                    d_drop_d_x: d.d_drop_d_x,
                    d_windage_d_x: d.d_windage_d_x,
                }
            })
            .collect();

        let total_cov = accumulate(&entries, None);
        let total_var = total_cov.a00 + total_cov.a11;
        let full_ellipse = ellipse_95(total_cov);
        // `Some(p_hit)` iff a target was supplied -- see "Hit probability" doc section above.
        let p_hit = target.map(|t| p_hit_bivariate(total_cov.a00, total_cov.a11, total_cov.a01, t));

        let mut sources_out: Vec<SourceContributionV1> = entries
            .iter()
            .map(|e| {
                let s2 = e.sigma * e.sigma;
                let this_var = e.d_drop_d_x * e.d_drop_d_x * s2 + e.d_windage_d_x * e.d_windage_d_x * s2;
                let reduced = accumulate(&entries, Some(e.axis));
                let reduced_ellipse = ellipse_95(reduced);
                // `target.zip(p_hit)` is `Some` exactly when `target` is (`p_hit` is computed
                // from `target` a few lines up), so this never diverges from `p_hit`'s own
                // Some-ness.
                let p_hit_gain_if_perfect = target.zip(p_hit).map(|(t, base_p_hit)| {
                    let raw = p_hit_bivariate(reduced.a00, reduced.a11, reduced.a01, t) - base_p_hit;
                    debug_assert!(
                        raw > -1e-2,
                        "perfecting {:?} at range {range_m} produced a meaningfully negative raw \
                         p_hit gain ({raw}) before clamping -- perfecting a source shrinks the \
                         impact covariance in the Loewner order, and shrinking a target-centred \
                         covariance that way cannot reduce the mass of a symmetric normal over a \
                         symmetric convex target (Anderson's theorem), so a value this far below \
                         zero means the quadrature or the excluded-source covariance is wrong, \
                         not ordinary numerical noise (measured worst-case quadrature error is \
                         under 1e-3 across a broad stress sweep, and orders of magnitude better \
                         for realistic target shapes -- see p_hit_bivariate's doc comment)",
                        e.axis
                    );
                    raw.max(0.0)
                });
                SourceContributionV1 {
                    axis: e.axis,
                    sigma: e.sigma,
                    d_drop_d_x: e.d_drop_d_x,
                    d_windage_d_x: e.d_windage_d_x,
                    scheme: e.scheme,
                    variance_share: if total_var > 0.0 { this_var / total_var } else { 0.0 },
                    ellipse_area_reduction_m2: (full_ellipse.area_m2 - reduced_ellipse.area_m2)
                        .max(0.0),
                    p_hit_gain_if_perfect,
                }
            })
            .collect();

        sort_by_variance_share_desc(&mut sources_out);
        let priority_statement = build_priority_statement(&sources_out, &unavailable, range_m);

        rows.push(ErrorBudgetRowV1 {
            range_m,
            sigma_drop_m: total_cov.a00.sqrt(),
            sigma_windage_m: total_cov.a11.sqrt(),
            covariance_m2: total_cov.a01,
            ellipse_95: full_ellipse,
            p_hit,
            sources: sources_out,
            priority_statement,
        });
    }

    let mut method = "central_difference_first_order_propagation".to_string();
    let mut assumptions = vec![
        "Declared sources are treated as independent; correlations between them are not \
         modelled."
            .to_string(),
        "Propagation is first-order (local linear) about the nominal solution, evaluated by \
         central differences through the real solver using each axis's own small default \
         step -- independent of the declared sigma, never a step scaled to it. A large \
         declared sigma is therefore a linear extrapolation from a slope measured over a \
         much smaller window, which is not exact for large or non-Gaussian input \
         uncertainty."
            .to_string(),
        "The 95% ellipse uses the chi-square 2-dof critical value 5.991 and assumes an \
         approximately Gaussian impact distribution."
            .to_string(),
        "A source's derivative may come from a one-sided (forward- or backward-only) \
         difference rather than a central one when its nominal value sits at a physical \
         domain boundary (for example, still air for wind speed); see that source's scheme \
         field. A one-sided difference has larger truncation error than a central one."
            .to_string(),
        "A source listed in unavailable_sources could not be evaluated for this request and \
         is excluded from every row's variance and ranking. That is not the same fact as a \
         source contributing zero -- it means this report cannot currently measure that \
         source's effect at all."
            .to_string(),
    ];
    if target.is_some() {
        method.push_str("_gl20_pm6sigma");
        assumptions.push(
            "Hit probability is the bivariate-normal mass over the target, computed by 20-point \
             Gauss-Legendre quadrature per smooth sub-interval, truncated at +/-6 sigma and \
             split at the target's own edge and (when two sources are strongly correlated) at \
             the drop values where the conditional windage window is crossed, so the fixed-order \
             rule is never applied across a hidden discontinuity or an unresolved sharp \
             transition. It reflects the declared input uncertainty only -- not model error -- \
             and assumes the target is centred on the aim point (the nominal trajectory's own \
             impact point), not offset from it."
                .to_string(),
        );
    }

    Ok(ErrorBudgetReportV1 {
        schema_version: ERROR_BUDGET_SCHEMA_VERSION_V1,
        method,
        assumptions,
        unavailable_sources: unavailable,
        rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbation::InputAxis;

    fn resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {}, "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    /// A QNH-referenced atmosphere fixture: `Altitude` is refused by `with_axis`
    /// (`AxisUnsupportedForRequest`) and must show up as unavailable, not silently vanish.
    fn qnh_resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {"altitude_m": 500.0, "temperature_k": 288.0, "pressure_pa": 101325.0,
                           "pressure_reference": "qnh"},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    // ---- Step 1 tests, verbatim from the brief (fixture literal fixed: FRAC_PI_2 instead of a
    // hardcoded pi/2 float, which clippy's approx_constant lint rejects -- see the task report).

    /// Acceptance criterion: a zero-uncertainty source contributes exactly zero.
    #[test]
    fn a_zero_sigma_source_contributes_exactly_zero() {
        let r = resolved();
        let rep = error_budget(&r, &[(InputAxis::MuzzleVelocityMps, 0.0),
                                     (InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        let mv = rep.rows[0].sources.iter()
            .find(|s| s.axis == InputAxis::MuzzleVelocityMps).unwrap();
        assert_eq!(mv.variance_share, 0.0);
        assert_eq!(mv.sigma, 0.0);
        // Beyond the brief: a zero-sigma source must not fabricate a zero derivative either, and
        // its ellipse-area reduction (a second, independently-computed quantity derived from the
        // SAME sigma=0) must also be exactly zero, not merely small.
        assert!(mv.d_drop_d_x != 0.0, "the real derivative must still be reported");
        assert_eq!(mv.ellipse_area_reduction_m2, 0.0);
    }

    /// Sources are preserved individually -- never collapsed into an "other" bucket.
    #[test]
    fn every_declared_source_appears_individually() {
        let r = resolved();
        let declared = [(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::WindSpeed, 1.0),
                        (InputAxis::BallisticCoefficient, 0.005)];
        let rep = error_budget(&r, &declared, &[600.0]).unwrap();
        assert_eq!(rep.rows[0].sources.len(), declared.len());
        for (axis, _) in declared {
            assert!(rep.rows[0].sources.iter().any(|s| s.axis == axis), "{axis:?} missing");
        }
    }

    /// Ranking must not depend on the order sources were declared in.
    #[test]
    fn ranking_is_invariant_to_declaration_order() {
        let r = resolved();
        let a = error_budget(&r, &[(InputAxis::MuzzleVelocityMps, 5.0),
                                   (InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        let b = error_budget(&r, &[(InputAxis::WindSpeed, 1.0),
                                   (InputAxis::MuzzleVelocityMps, 5.0)], &[600.0]).unwrap();
        let order_a: Vec<_> = a.rows[0].sources.iter().map(|s| s.axis).collect();
        let order_b: Vec<_> = b.rows[0].sources.iter().map(|s| s.axis).collect();
        assert_eq!(order_a, order_b);
    }

    /// Variance shares are normalised.
    #[test]
    fn variance_shares_sum_to_one() {
        let r = resolved();
        let rep = error_budget(&r, &[(InputAxis::MuzzleVelocityMps, 5.0),
                                     (InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        let sum: f64 = rep.rows[0].sources.iter().map(|s| s.variance_share).sum();
        assert!((sum - 1.0).abs() < 1e-9, "shares summed to {sum}");
    }

    #[test]
    fn the_report_declares_independence_and_linearity() {
        let r = resolved();
        let rep = error_budget(&r, &[(InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        assert_eq!(rep.method, "central_difference_first_order_propagation");
        assert!(rep.assumptions.iter().any(|s| s.contains("independent")));
        assert!(rep.assumptions.iter().any(|s| s.to_lowercase().contains("linear")));
    }

    /// The honesty requirement names three specific claims the payload must carry: independence,
    /// first-order/local-linearity, AND the Gaussian-ellipse assumption. The test above (verbatim
    /// from the brief) only pins down the first two; this pins down the third explicitly, plus
    /// the two extra assumptions this implementation adds (non-central schemes,
    /// unavailable-is-not-zero) -- so all five payload strings are independently checked for their
    /// OWN specific content, not just "the list is non-empty."
    #[test]
    fn the_report_declares_the_gaussian_ellipse_assumption_and_the_two_added_caveats() {
        let r = resolved();
        let rep = error_budget(&r, &[(InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        assert!(
            rep.assumptions.iter().any(|s| s.contains("5.991") && s.to_lowercase().contains("gaussian")),
            "no assumption states the chi-square constant and the Gaussian-ellipse assumption: {:#?}",
            rep.assumptions
        );
        assert!(
            rep.assumptions.iter().any(|s| s.to_lowercase().contains("one-sided")),
            "no assumption warns that a source's derivative may be one-sided: {:#?}",
            rep.assumptions
        );
        assert!(
            rep.assumptions.iter().any(|s| s.contains("unavailable_sources")
                && s.to_lowercase().contains("not the same fact as")),
            "no assumption distinguishes an unavailable source from a zero-contribution one: {:#?}",
            rep.assumptions
        );
    }

    // ---- Beyond the brief: the four blind spots the task explicitly calls out, plus the
    // determinism/cost/honesty requirements it names as "things the brief does not know."

    /// (1) A one-sided scheme must be visible in the report, not silently indistinguishable from
    /// a central one. Still air is this report's flagship scenario (a wind-call sigma with
    /// `speed_mps: 0.0`), and it deterministically forces `ForwardOneSided` (see
    /// `crate::perturbation::derive`'s own `wind_speed_falls_back_to_one_sided_in_still_air`).
    #[test]
    fn a_non_central_scheme_is_surfaced_in_the_source_and_the_priority_statement() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0},
            "atmosphere": {},
            "wind": {"speed_mps": 0.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let rep = error_budget(&r, &[(InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        let ws = &rep.rows[0].sources[0];
        assert_eq!(ws.axis, InputAxis::WindSpeed);
        assert_eq!(ws.scheme, DifferenceScheme::ForwardOneSided);
        assert!(
            rep.rows[0].priority_statement.contains("one-sided"),
            "priority_statement should flag a one-sided dominant source: {}",
            rep.rows[0].priority_statement
        );
        assert!(rep.assumptions.iter().any(|s| s.to_lowercase().contains("one-sided")));
    }

    /// (2) An unavailable source must be RECORDED, not silently dropped -- and must be
    /// distinguishable from a source that evaluated and contributed zero. Also confirms the rest
    /// of the report (a source that DID evaluate) is unaffected by a sibling's unavailability.
    #[test]
    fn an_unavailable_source_is_recorded_not_silently_dropped() {
        let r = qnh_resolved();
        let rep = error_budget(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::Altitude, 50.0)],
            &[300.0],
        )
        .unwrap();

        // Altitude must NEVER appear as an evaluated source...
        assert!(rep.rows[0].sources.iter().all(|s| s.axis != InputAxis::Altitude));
        // ...but MUST appear, explicitly, as unavailable, with a reason naming the mechanism.
        let skipped = rep
            .unavailable_sources
            .iter()
            .find(|u| u.axis == InputAxis::Altitude)
            .expect("Altitude must be recorded as unavailable, not dropped");
        assert_eq!(skipped.sigma, 50.0);
        assert_eq!(skipped.code, UnavailableReasonCodeV1::AxisUnsupportedForRequest);
        assert!(skipped.reason.to_lowercase().contains("qnh"), "{}", skipped.reason);
        // The unavailable list is not merely non-empty by accident -- it must be EXACTLY the one
        // axis that was actually refused, not every declared axis.
        assert_eq!(rep.unavailable_sources.len(), 1);

        // The sibling source that DID evaluate must be reported normally and dominate (it is the
        // only evaluated source), not be suppressed or zeroed by Altitude's unavailability.
        assert_eq!(rep.rows[0].sources.len(), 1);
        let mv = &rep.rows[0].sources[0];
        assert_eq!(mv.axis, InputAxis::MuzzleVelocityMps);
        assert!(mv.variance_share > 0.0);
        assert!(
            rep.rows[0].priority_statement.contains("could not be evaluated"),
            "priority_statement should mention the unavailable source too: {}",
            rep.rows[0].priority_statement
        );
    }

    /// (2), continued: when EVERY declared source is unavailable, the report must still succeed
    /// (not error out) with an empty `sources` list, a priority statement that says so plainly
    /// (never confused with "every sigma is zero"), and a well-defined (zero) ellipse.
    #[test]
    fn every_source_unavailable_still_produces_a_well_formed_report() {
        let r = qnh_resolved();
        let rep = error_budget(&r, &[(InputAxis::Altitude, 50.0)], &[300.0]).unwrap();
        assert!(rep.rows[0].sources.is_empty());
        assert_eq!(rep.unavailable_sources.len(), 1);
        assert_eq!(rep.unavailable_sources[0].axis, InputAxis::Altitude);
        assert_eq!(
            rep.unavailable_sources[0].code,
            UnavailableReasonCodeV1::AxisUnsupportedForRequest
        );
        assert_eq!(rep.rows[0].ellipse_95.area_m2, 0.0);
        assert!(rep.rows[0].priority_statement.contains("None of the declared sources"));
    }

    /// (2), continued once more: the classification itself, tested directly and exhaustively
    /// against every `KernelError` variant, independent of any physics fixture. This is the
    /// authoritative proof that a genuine (non-structural) failure is classified as "propagate,"
    /// not "unavailable" -- constructing a REAL end-to-end fixture that makes
    /// `central_difference` return a genuine `Solve`/`Observation` failure through
    /// `error_budget`'s fixed (default-step) call is only reachable via a knife-edge convergence
    /// boundary that would make the test fragile against unrelated solver changes; the
    /// exhaustive match in production code (no wildcard arm) is the stronger guarantee here,
    /// since it also protects a FUTURE `KernelError` variant at compile time. See the task
    /// report for the full reasoning.
    ///
    /// (Review I4): `DuplicateAxis` is in the GENUINE bucket, not the structural one --
    /// `error_budget` constructs and returns it directly from its own up-front validation
    /// (never through `central_difference`), so by the time anything reaches this
    /// classification it must propagate, exactly like a malformed request.
    ///
    /// Also verifies the SPECIFIC `UnavailableReasonCodeV1` returned for each structural
    /// refusal, not just that one was returned -- `is_some()` alone would not catch two codes
    /// swapped between variants (e.g. `AxisAbsent` mistakenly tagged `CategoricalAxis`).
    #[test]
    fn unavailable_reason_classifies_every_kernel_error_variant() {
        use crate::solve_json::SolveErrorCodeV1;
        use crate::trajectory_observation::TrajectoryObservationError;

        let structural = [
            (
                KernelError::AxisUnsupportedForRequest { axis: InputAxis::Altitude, reason: "x" },
                UnavailableReasonCodeV1::AxisUnsupportedForRequest,
            ),
            (KernelError::AxisAbsent(InputAxis::WindSpeed), UnavailableReasonCodeV1::AxisAbsent),
            (
                KernelError::CategoricalAxis(InputAxis::CoriolisEnabled),
                UnavailableReasonCodeV1::CategoricalAxis,
            ),
            (
                KernelError::StepOutOfDomain { axis: InputAxis::RelativeHumidity, attempted: 2.0 },
                UnavailableReasonCodeV1::StepOutOfDomain,
            ),
        ];
        for (e, expected_code) in &structural {
            match unavailable_reason(e) {
                Some((code, _)) => assert_eq!(
                    code, *expected_code,
                    "{e:?} classified with the wrong UnavailableReasonCodeV1"
                ),
                None => panic!("{e:?} must be classified as unavailable (recorded), not propagated"),
            }
        }

        let genuine = [
            KernelError::Solve { code: SolveErrorCodeV1::SolveFailed, message: "x".into() },
            KernelError::Observation(TrajectoryObservationError::NonMonotonicTrajectory {
                index: 3,
                previous_distance_m: 10.0,
                distance_m: 9.0,
            }),
            KernelError::TypeMismatch(InputAxis::Mass),
            KernelError::NonFinite(InputAxis::Mass),
            KernelError::DuplicateAxis(InputAxis::Mass),
        ];
        for e in &genuine {
            assert!(
                unavailable_reason(e).is_none(),
                "{e:?} must be classified as a genuine failure (propagated), not recorded"
            );
        }
    }

    /// (3) Ranking must be deterministic even when two sources GENUINELY tie -- a real physics
    /// fixture cannot reliably produce a bit-exact tie, so this constructs one directly against
    /// the sort function itself. Without the tie-break (comparing only `variance_share`), Rust's
    /// stable sort would preserve each input's own relative order for the tied pair, so the two
    /// differently-ordered inputs below would disagree -- this test fails under that
    /// implementation (verified while developing it) and passes only because the tie is broken
    /// on a fixed key.
    #[test]
    fn tied_variance_shares_break_deterministically_regardless_of_input_order() {
        fn stub(axis: InputAxis, variance_share: f64) -> SourceContributionV1 {
            SourceContributionV1 {
                axis,
                sigma: 1.0,
                d_drop_d_x: 0.0,
                d_windage_d_x: 0.0,
                scheme: DifferenceScheme::Central,
                variance_share,
                ellipse_area_reduction_m2: 0.0,
                p_hit_gain_if_perfect: None,
            }
        }
        let mut a = vec![stub(InputAxis::WindSpeed, 0.5), stub(InputAxis::MuzzleVelocityMps, 0.5)];
        let mut b = vec![stub(InputAxis::MuzzleVelocityMps, 0.5), stub(InputAxis::WindSpeed, 0.5)];
        sort_by_variance_share_desc(&mut a);
        sort_by_variance_share_desc(&mut b);
        let order_a: Vec<_> = a.iter().map(|s| s.axis).collect();
        let order_b: Vec<_> = b.iter().map(|s| s.axis).collect();
        assert_eq!(order_a, order_b, "a tie must break the same way regardless of input order");
        // Pin down WHICH order, not just "some order both agree on": "MuzzleVelocityMps" sorts
        // before "WindSpeed" as a Debug string.
        assert_eq!(order_a, vec![InputAxis::MuzzleVelocityMps, InputAxis::WindSpeed]);
    }

    /// (3), continued: a three-way tie (beyond a simple pairwise swap) stays fully deterministic
    /// across every rotation of the input order.
    #[test]
    fn a_three_way_tie_is_fully_deterministic_across_every_rotation() {
        fn stub(axis: InputAxis) -> SourceContributionV1 {
            SourceContributionV1 {
                axis,
                sigma: 1.0,
                d_drop_d_x: 0.0,
                d_windage_d_x: 0.0,
                scheme: DifferenceScheme::Central,
                variance_share: 1.0 / 3.0,
                ellipse_area_reduction_m2: 0.0,
                p_hit_gain_if_perfect: None,
            }
        }
        let axes = [InputAxis::WindSpeed, InputAxis::MuzzleVelocityMps, InputAxis::Mass];
        let mut orders = Vec::new();
        for rotation in 0..axes.len() {
            let mut rotated: Vec<SourceContributionV1> =
                (0..axes.len()).map(|k| stub(axes[(k + rotation) % axes.len()])).collect();
            sort_by_variance_share_desc(&mut rotated);
            orders.push(rotated.iter().map(|s| s.axis).collect::<Vec<_>>());
        }
        for w in orders.windows(2) {
            assert_eq!(w[0], w[1], "every rotation of a full tie must sort identically");
        }
    }

    /// Would a test notice two sources' contributions transposed? Verified against
    /// INDEPENDENTLY computed quantities (direct `central_difference` calls made in this test,
    /// not `error_budget`'s own internal numbers) for two axes with genuinely different
    /// sensitivities, rather than only checking self-consistency.
    #[test]
    fn two_sources_contributions_are_not_transposed() {
        let r = resolved();
        let mv_sigma = 5.0_f64;
        let ws_sigma = 1.0_f64;
        let rep = error_budget(
            &r,
            &[(InputAxis::MuzzleVelocityMps, mv_sigma), (InputAxis::WindSpeed, ws_sigma)],
            &[600.0],
        )
        .unwrap();

        let mv_deriv = central_difference(&r, InputAxis::MuzzleVelocityMps, &[600.0], None)
            .unwrap()[0];
        let ws_deriv = central_difference(&r, InputAxis::WindSpeed, &[600.0], None).unwrap()[0];

        let mv_var = (mv_deriv.d_drop_d_x * mv_sigma).powi(2)
            + (mv_deriv.d_windage_d_x * mv_sigma).powi(2);
        let ws_var = (ws_deriv.d_drop_d_x * ws_sigma).powi(2)
            + (ws_deriv.d_windage_d_x * ws_sigma).powi(2);
        let independent_total = mv_var + ws_var;

        let mv_row = rep.rows[0].sources.iter().find(|s| s.axis == InputAxis::MuzzleVelocityMps)
            .unwrap();
        let ws_row = rep.rows[0].sources.iter().find(|s| s.axis == InputAxis::WindSpeed).unwrap();

        // Raw derivatives and sigma must match the independently-computed kernel call exactly --
        // this is what a transposition (assigning WindSpeed's numbers to the MuzzleVelocityMps
        // row or vice versa) would break.
        assert_eq!(mv_row.d_drop_d_x, mv_deriv.d_drop_d_x);
        assert_eq!(mv_row.d_windage_d_x, mv_deriv.d_windage_d_x);
        assert_eq!(mv_row.sigma, mv_sigma);
        assert_eq!(ws_row.d_drop_d_x, ws_deriv.d_drop_d_x);
        assert_eq!(ws_row.d_windage_d_x, ws_deriv.d_windage_d_x);
        assert_eq!(ws_row.sigma, ws_sigma);

        // variance_share compared against a total computed OUTSIDE error_budget entirely.
        assert!((mv_row.variance_share - mv_var / independent_total).abs() < 1e-9);
        assert!((ws_row.variance_share - ws_var / independent_total).abs() < 1e-9);

        // Sanity: the two shares are not (nearly) equal, or a transposition would be invisible
        // to the assertions above.
        assert!(
            (mv_row.variance_share - ws_row.variance_share).abs() > 0.05,
            "fixture must give the two sources distinguishably different shares: mv={}, ws={}",
            mv_row.variance_share,
            ws_row.variance_share
        );
    }

    /// A second independent check that shares sum to one, using a total computed OUTSIDE
    /// `error_budget`'s own arithmetic (three sources this time, not two, and reusing the
    /// independent per-source variances rather than re-deriving `error_budget`'s own `total_var`
    /// field).
    #[test]
    fn variance_shares_sum_to_one_against_an_independently_recomputed_total() {
        let r = resolved();
        let declared = [
            (InputAxis::MuzzleVelocityMps, 5.0),
            (InputAxis::WindSpeed, 1.0),
            (InputAxis::BallisticCoefficient, 0.005),
        ];
        let rep = error_budget(&r, &declared, &[600.0]).unwrap();

        let mut independent_total = 0.0;
        let mut independent_var = std::collections::HashMap::new();
        for &(axis, sigma) in &declared {
            let d = central_difference(&r, axis, &[600.0], None).unwrap()[0];
            let v = (d.d_drop_d_x * sigma).powi(2) + (d.d_windage_d_x * sigma).powi(2);
            independent_total += v;
            independent_var.insert(axis, v);
        }

        let mut share_sum = 0.0;
        for s in &rep.rows[0].sources {
            let expected_share = independent_var[&s.axis] / independent_total;
            assert!(
                (s.variance_share - expected_share).abs() < 1e-9,
                "{:?}: report said {}, independently expected {}",
                s.axis,
                s.variance_share,
                expected_share
            );
            share_sum += s.variance_share;
        }
        assert!((share_sum - 1.0).abs() < 1e-9, "shares summed to {share_sum}");
    }

    /// (I6, review round) FOUR previously-unasserted public payload fields --
    /// `sigma_drop_m`, `sigma_windage_m`, `covariance_m2`, and `ellipse_95.rotation_rad` -- each
    /// checked against an INDEPENDENTLY computed value, never against `error_budget`'s own
    /// internal `Symmetric2`/`accumulate` bookkeeping. Each of these mutations passes all of this
    /// module's OTHER tests: swapping `sigma_drop_m`/`sigma_windage_m`; forcing `rotation_rad` to
    /// `0.0` unconditionally; forcing `covariance_m2` to `0.0` unconditionally. (The fifth
    /// previously-unasserted field, `ellipse_area_reduction_m2`, is covered by the next test,
    /// which needs 3+ sources to be discriminating at all -- see that field's own doc comment.)
    ///
    /// `Cant` is declared at a NONZERO baseline (`cant_angle_rad: 0.5`, ~29 degrees) specifically
    /// so its derivative has comparable, clearly nonzero components in BOTH drop and windage: at
    /// a baseline cant of exactly zero, canting the rifle by an infinitesimal `dtheta` rotates a
    /// purely-vertical drop vector into windage to FIRST order while changing its own magnitude
    /// only to SECOND order (`d(drop)/d(cant) ~ 0`, `d(windage)/d(cant) ~ drop`), which would
    /// make the covariance term a floating-point-noise-sized artifact rather than a robustly
    /// nonzero cross term -- at a nonzero baseline, both `sin(cant)` and `cos(cant)` are
    /// appreciable, giving `Cant` a genuinely mixed drop/windage sensitivity.
    #[test]
    fn sigma_covariance_and_rotation_are_verified_independently_not_against_themselves() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0, "cant_angle_rad": 0.5},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 25.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let declared = [
            (InputAxis::MuzzleVelocityMps, 5.0),
            (InputAxis::WindSpeed, 1.0),
            (InputAxis::Cant, 0.01),
        ];
        let rep = error_budget(&r, &declared, &[600.0]).unwrap();

        // Independent oracle: direct central_difference calls, summed by hand right here -- not
        // error_budget's own accumulate()/Symmetric2 code path.
        let mut var_drop = 0.0_f64;
        let mut var_wind = 0.0_f64;
        let mut cov = 0.0_f64;
        for &(axis, sigma) in &declared {
            let d = central_difference(&r, axis, &[600.0], None).unwrap()[0];
            let s2 = sigma * sigma;
            var_drop += d.d_drop_d_x * d.d_drop_d_x * s2;
            var_wind += d.d_windage_d_x * d.d_windage_d_x * s2;
            cov += d.d_drop_d_x * d.d_windage_d_x * s2;
        }

        let row = &rep.rows[0];

        // sigma_drop_m / sigma_windage_m: independently computed AND distinguishably different
        // (measured ~0.045 vs ~0.217 for this fixture), so a swap between them is not invisible.
        assert!(
            (row.sigma_drop_m - var_drop.sqrt()).abs() < 1e-9,
            "sigma_drop_m: got {}, expected {}",
            row.sigma_drop_m,
            var_drop.sqrt()
        );
        assert!(
            (row.sigma_windage_m - var_wind.sqrt()).abs() < 1e-9,
            "sigma_windage_m: got {}, expected {}",
            row.sigma_windage_m,
            var_wind.sqrt()
        );
        assert!(
            (row.sigma_drop_m - row.sigma_windage_m).abs()
                > 0.1 * row.sigma_drop_m.max(row.sigma_windage_m),
            "fixture must give sigma_drop_m and sigma_windage_m distinguishably different \
             values, or a swap between them would be invisible here: drop={}, wind={}",
            row.sigma_drop_m,
            row.sigma_windage_m
        );

        // covariance_m2: independently computed, and clearly nonzero (measured ~ -1.5e-4, about
        // 7.7% of var_drop -- not a floating-point-noise-sized artifact).
        assert!(
            (row.covariance_m2 - cov).abs() < 1e-9 * cov.abs().max(1.0),
            "covariance_m2: got {}, expected {}",
            row.covariance_m2,
            cov
        );
        assert!(cov.abs() > 1e-5, "fixture must give a clearly nonzero covariance: {cov}");

        // rotation_rad: verified via the DEFINITION of the major-axis angle (a Rayleigh-quotient
        // check), not by re-deriving the same atan2 formula error_budget itself uses to compute
        // it -- the variance PROJECTED along the direction rotation_rad points must equal the
        // ellipse's own major-axis eigenvalue (semi_major_m^2 / CHI2_95_2DOF). Forcing
        // rotation_rad to 0.0 would project onto the drop axis instead, giving var_drop
        // (~0.002) rather than the true lambda_max (~0.047 for this fixture) -- clearly caught.
        let (s, c) = row.ellipse_95.rotation_rad.sin_cos();
        let var_along_major = c * c * var_drop + 2.0 * c * s * cov + s * s * var_wind;
        let lambda_max = row.ellipse_95.semi_major_m.powi(2) / CHI2_95_2DOF;
        assert!(
            (var_along_major - lambda_max).abs() < 1e-9 * lambda_max.max(1.0),
            "rotation_rad ({}) does not point along the major axis: variance projected along it \
             is {}, but the major-axis eigenvalue is {}",
            row.ellipse_95.rotation_rad,
            var_along_major,
            lambda_max
        );
        assert_ne!(
            row.ellipse_95.rotation_rad, 0.0,
            "fixture must give a genuinely nonzero rotation"
        );
    }

    /// (I7, review round) `ellipse_area_reduction_m2` is DEGENERATE with two or fewer declared
    /// sources (documented on the field: removing either of exactly two leaves a rank-1
    /// covariance, area exactly zero, so BOTH report the full ellipse area as their own
    /// reduction even at a lopsided variance split). With three or more sources it is
    /// generically discriminating; this proves that, and checks every source's value against an
    /// INDEPENDENT oracle that does not call `Symmetric2`/`ellipse_95` at all (a hand-rolled
    /// trace/determinant formula -- the same shape `monte_carlo::calculate_confidence_ellipse`
    /// uses, deliberately a DIFFERENT numerical path than this module's own).
    #[test]
    fn ellipse_area_reduction_is_nonzero_and_discriminating_with_three_or_more_sources() {
        fn independent_ellipse_area(var_drop: f64, var_wind: f64, cov: f64) -> f64 {
            let trace = var_drop + var_wind;
            let det = (var_drop * var_wind - cov * cov).max(0.0);
            let disc = ((trace * trace / 4.0) - det).max(0.0).sqrt();
            let l1 = (trace / 2.0 + disc).max(0.0);
            let l2 = (trace / 2.0 - disc).max(0.0);
            std::f64::consts::PI * (CHI2_95_2DOF * l1).sqrt() * (CHI2_95_2DOF * l2).sqrt()
        }

        let r = resolved();
        // One dominant source (MV) plus two much smaller ones -- a realistic "mostly one input
        // matters" declaration, exactly the shape I7 warns reads as a false tie under the
        // two-source degeneracy above.
        let declared = [
            (InputAxis::MuzzleVelocityMps, 5.0),
            (InputAxis::WindSpeed, 0.3),
            (InputAxis::BallisticCoefficient, 0.001),
        ];
        let rep = error_budget(&r, &declared, &[600.0]).unwrap();
        assert_eq!(rep.rows[0].sources.len(), 3);

        let mut derivs = std::collections::HashMap::new();
        for &(axis, sigma) in &declared {
            let d = central_difference(&r, axis, &[600.0], None).unwrap()[0];
            derivs.insert(axis, (sigma, d.d_drop_d_x, d.d_windage_d_x));
        }
        let variance_excluding = |exclude: Option<InputAxis>| -> (f64, f64, f64) {
            let mut vd = 0.0;
            let mut vw = 0.0;
            let mut cv = 0.0;
            for (&axis, &(sigma, dd, dw)) in &derivs {
                if Some(axis) == exclude {
                    continue;
                }
                let s2 = sigma * sigma;
                vd += dd * dd * s2;
                vw += dw * dw * s2;
                cv += dd * dw * s2;
            }
            (vd, vw, cv)
        };
        let (fvd, fvw, fcv) = variance_excluding(None);
        let full_area = independent_ellipse_area(fvd, fvw, fcv);

        let mut reductions = Vec::new();
        for s in &rep.rows[0].sources {
            let (vd, vw, cv) = variance_excluding(Some(s.axis));
            let reduced_area = independent_ellipse_area(vd, vw, cv);
            let expected_reduction = (full_area - reduced_area).max(0.0);
            assert!(
                (s.ellipse_area_reduction_m2 - expected_reduction).abs()
                    < 1e-9 * expected_reduction.max(1.0),
                "{:?}: got {}, independently expected {}",
                s.axis,
                s.ellipse_area_reduction_m2,
                expected_reduction
            );
            reductions.push((s.axis, s.ellipse_area_reduction_m2));
        }

        // Discriminating, not a "0.0 unconditionally" mutation (which would fail the exact
        // check above already) nor a same-for-everyone tie (the I7 degeneracy, which no longer
        // applies with 3+ sources).
        assert!(
            reductions.iter().all(|&(_, red)| red > 0.0),
            "every reduction should be positive with 3+ sources: {reductions:?}"
        );
        let first = reductions[0].1;
        assert!(
            reductions.iter().any(|&(_, red)| (red - first).abs() > 1e-6 * first.max(1.0)),
            "reductions must discriminate between sources with 3+ declared, not all be equal: \
             {reductions:?}"
        );
    }

    /// (4) Jacobian reuse across ranges: each row must carry its OWN correctly-indexed
    /// derivative (checked against an independent `central_difference` call at that same single
    /// range), and different ranges must give genuinely different numbers -- catching a bug that
    /// indexed every row from row 0's derivative instead of `derivatives[i]`.
    #[test]
    fn error_budget_computes_a_correctly_indexed_row_per_requested_range() {
        let r = resolved();
        let ranges = [300.0_f64, 600.0_f64, 850.0_f64];
        let rep = error_budget(&r, &[(InputAxis::MuzzleVelocityMps, 5.0)], &ranges).unwrap();
        assert_eq!(rep.rows.len(), ranges.len());

        for (i, &range_m) in ranges.iter().enumerate() {
            assert_eq!(rep.rows[i].range_m, range_m);
            let expected = central_difference(&r, InputAxis::MuzzleVelocityMps, &[range_m], None)
                .unwrap()[0];
            let got = &rep.rows[i].sources[0];
            assert_eq!(got.d_drop_d_x, expected.d_drop_d_x, "range {range_m}");
            assert_eq!(got.d_windage_d_x, expected.d_windage_d_x, "range {range_m}");
        }
        // Sensitivity to muzzle velocity grows with range: row 2 must differ meaningfully from
        // row 0, or a bug that copied row 0 into every row would pass the per-row checks above
        // (they'd coincidentally still match central_difference at range 300 for EVERY row only
        // if 300 were queried each time -- it is not, so this also guards indexing directly).
        assert!(
            rep.rows[2].sources[0].d_drop_d_x.abs() > rep.rows[0].sources[0].d_drop_d_x.abs() * 2.0
        );
    }

    /// Declaring a non-finite or negative sigma is rejected up front, before any solve -- it is
    /// not a real one-sigma uncertainty regardless of whether anything downstream would panic on
    /// it. (An earlier revision of this comment cited `variance_share`'s sort comparator
    /// panicking on NaN as the reason; that specific mechanism no longer applies now that the
    /// sort uses `f64::total_cmp` (review M1), but the input is still nonsensical and still
    /// rejected -- see `an_astronomically_large_sigma_does_not_panic_the_sort` below for the
    /// remaining, more extreme case `total_cmp` exists to handle instead of panicking on.)
    #[test]
    fn a_non_finite_or_negative_sigma_is_rejected() {
        let r = resolved();
        let nan = error_budget(&r, &[(InputAxis::WindSpeed, f64::NAN)], &[600.0]);
        assert!(matches!(nan, Err(KernelError::NonFinite(InputAxis::WindSpeed))));
        let neg = error_budget(&r, &[(InputAxis::WindSpeed, -1.0)], &[600.0]);
        assert!(matches!(neg, Err(KernelError::NonFinite(InputAxis::WindSpeed))));
        let inf = error_budget(&r, &[(InputAxis::WindSpeed, f64::INFINITY)], &[600.0]);
        assert!(matches!(inf, Err(KernelError::NonFinite(InputAxis::WindSpeed))));
    }

    /// (M1, review round) An astronomically large (but finite, non-negative -- so it PASSES the
    /// validation above) declared sigma can overflow `sigma * sigma` to infinity, making
    /// `variance_share` a NaN (`inf / inf`). `sort_by_variance_share_desc` must not panic on
    /// that: this crate ships a fuzz suite, and a library function panicking on an
    /// extreme-but-technically-valid input is a real failure mode.
    #[test]
    fn an_astronomically_large_sigma_does_not_panic_the_sort() {
        let r = resolved();
        let huge = 1e200_f64;
        assert!(
            huge.is_finite() && (huge * huge).is_infinite(),
            "fixture assumption: sigma^2 must overflow to infinity"
        );
        let rep = error_budget(
            &r,
            &[(InputAxis::MuzzleVelocityMps, huge), (InputAxis::WindSpeed, 1.0)],
            &[600.0],
        )
        .unwrap();
        // Does not panic; the NaN share's exact placement by total_cmp is not physically
        // meaningful for this pathological input and is not asserted, only that the call
        // returns normally with every declared source still present.
        assert_eq!(rep.rows[0].sources.len(), 2);
    }

    /// (I4, review round) The same axis declared twice must be rejected up front, not silently
    /// double-counted -- with duplicates, `accumulate(&entries, Some(axis))` (keyed on axis)
    /// would exclude BOTH entries when pricing "if this source were measured perfectly" for
    /// EITHER one, computing a counterfactual against a baseline where the caller's OTHER
    /// declaration of the same axis was ALSO perfected -- not what either individual declaration
    /// means.
    #[test]
    fn a_duplicate_axis_declaration_is_rejected() {
        let r = resolved();
        let e = error_budget(
            &r,
            &[(InputAxis::WindSpeed, 1.0), (InputAxis::WindSpeed, 2.0)],
            &[600.0],
        );
        assert!(matches!(e, Err(KernelError::DuplicateAxis(InputAxis::WindSpeed))));
    }

    /// (I4, continued) The duplicate check finds a repeat anywhere in the list, not just
    /// adjacent entries, and names the axis that was actually repeated.
    #[test]
    fn a_duplicate_axis_is_found_even_when_not_adjacent() {
        let r = resolved();
        let e = error_budget(
            &r,
            &[
                (InputAxis::MuzzleVelocityMps, 5.0),
                (InputAxis::BallisticCoefficient, 0.005),
                (InputAxis::MuzzleVelocityMps, 6.0),
            ],
            &[600.0],
        );
        assert!(matches!(e, Err(KernelError::DuplicateAxis(InputAxis::MuzzleVelocityMps))));
    }

    /// (I3, review round) A range beyond the BASE request's own `max_range_m` must be rejected
    /// directly as `Observation(OutOfRange)`, not silently converted into "every declared source
    /// is unavailable" -- without this check, `central_difference` would fail BOTH perturbed
    /// sides of EVERY axis identically on the same out-of-range observation
    /// (`StepOutOfDomain`), recording each one as unavailable with a reason blaming that axis's
    /// own step, when the real cause is the query itself and has nothing to do with any axis.
    #[test]
    fn a_range_beyond_max_range_m_is_rejected_directly_not_laundered_per_axis() {
        let r = resolved(); // max_range_m: 900.0
        let e = error_budget(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::WindSpeed, 1.0)],
            &[600.0, 5000.0],
        );
        match e {
            Err(KernelError::Observation(TrajectoryObservationError::OutOfRange {
                requested_m,
                maximum_m,
                ..
            })) => {
                assert_eq!(requested_m, 5000.0);
                assert_eq!(maximum_m, 900.0);
            }
            other => panic!(
                "expected Err(KernelError::Observation(OutOfRange {{ .. }})) naming the \
                 out-of-range query, got {other:?}"
            ),
        }
    }

    /// (I3, continued) A negative range is likewise rejected directly, matching
    /// `observation_at_range_checked`'s own `distance_m < minimum_m` check.
    #[test]
    fn a_negative_range_is_rejected_directly() {
        let r = resolved();
        let e = error_budget(&r, &[(InputAxis::WindSpeed, 1.0)], &[-1.0]);
        assert!(matches!(
            e,
            Err(KernelError::Observation(TrajectoryObservationError::OutOfRange {
                requested_m,
                ..
            })) if requested_m == -1.0
        ));
    }

    /// (I1, review round) The `StepOutOfDomain` reason names the REAL trigger (the axis's own
    /// default differencing step) and does not blame the declared sigma, which has zero
    /// influence on whether this fires -- `error_budget` always requests `central_difference`'s
    /// default step, never a custom one.
    #[test]
    fn step_out_of_domain_reason_blames_the_default_step_not_the_sigma() {
        let e = KernelError::StepOutOfDomain { axis: InputAxis::RelativeHumidity, attempted: 2.0 };
        let (code, reason) = unavailable_reason(&e).unwrap();
        assert_eq!(code, UnavailableReasonCodeV1::StepOutOfDomain);
        // Names the real trigger (the axis's own default step)...
        assert!(reason.to_lowercase().contains("default step"), "{reason}");
        // ...and explicitly disclaims sigma as the cause, rather than implying sigma's
        // MAGNITUDE is what pushed the step out of domain (the false claim I1 found in this
        // module's top-of-file doc comment, now corrected there too). The word "sigma" appearing
        // at all is fine and expected here -- the message says "does not depend on the declared
        // sigma" precisely to rule that out -- so this checks for the wrong-causation PHRASING,
        // not for the word's mere presence.
        assert!(
            !reason.to_lowercase().contains("sigma larger"),
            "reason should not claim the declared sigma's SIZE caused this: {reason}"
        );
        assert!(
            reason.to_lowercase().contains("does not depend on the declared sigma"),
            "{reason}"
        );
    }

    /// Declaring zero sources is well-defined, not a panic or a spurious error.
    #[test]
    fn declaring_no_sources_is_well_defined() {
        let r = resolved();
        let rep = error_budget(&r, &[], &[600.0]).unwrap();
        assert!(rep.rows[0].sources.is_empty());
        assert!(rep.unavailable_sources.is_empty());
        assert_eq!(rep.rows[0].ellipse_95.area_m2, 0.0);
        assert!(rep.rows[0].priority_statement.contains("No sources were declared"));
    }

    /// The 95% ellipse for a single dominant source matches a hand-computable closed form: with
    /// only one source, the impact covariance is an exact rank-1 outer product
    /// `sigma^2 * [dd, dw] * [dd, dw]^T`, whose only nonzero eigenvalue is
    /// `sigma^2 * (dd^2 + dw^2)` -- an independent check on `Symmetric2::largest_smallest_eigenvalues`
    /// and `ellipse_95` together, not just on `error_budget`'s own bookkeeping.
    ///
    /// (Review I5): the minor axis and area are zero in EXACT arithmetic for a rank-1 covariance,
    /// but NOT bit-exactly in IEEE-754: `determinant = fl(dd^2 s^2 * dw^2 s^2) -
    /// fl((dd * dw * s^2)^2)` is two differently-rounded computations of the same mathematical
    /// product and can land a few ULP to either side of zero depending on `dd`/`dw`'s exact bit
    /// pattern (a catastrophic-cancellation difference of two near-equal drops, hence
    /// libm/compiler/platform sensitive) -- this file's own tests run on Darwin, CI's
    /// `ci.yml` runs on `ubuntu-latest`/glibc. A prior version of this test asserted bit-exact
    /// `== 0.0`, which is a genuine, unverified CI coin-flip, not a property this code actually
    /// guarantees (`largest_smallest_eigenvalues` clamps a NEGATIVE determinant to exactly zero,
    /// but a tiny POSITIVE one survives as a tiny nonzero `smallest`). Bound both RELATIVE to the
    /// major axis instead.
    #[test]
    fn single_source_ellipse_matches_the_closed_form_rank_one_case() {
        let r = resolved();
        let sigma = 1.0_f64;
        let rep = error_budget(&r, &[(InputAxis::WindSpeed, sigma)], &[600.0]).unwrap();
        let d = central_difference(&r, InputAxis::WindSpeed, &[600.0], None).unwrap()[0];

        let expected_largest = sigma * sigma * (d.d_drop_d_x.powi(2) + d.d_windage_d_x.powi(2));
        let expected_major = (CHI2_95_2DOF * expected_largest).sqrt();
        let ellipse = rep.rows[0].ellipse_95;
        assert!(
            ellipse.semi_minor_m <= 1e-9 * expected_major,
            "semi_minor_m should be negligible relative to the major axis for a rank-1 \
             covariance, got minor={}, major={}",
            ellipse.semi_minor_m,
            expected_major
        );
        assert!(
            ellipse.area_m2 <= 1e-9 * std::f64::consts::PI * expected_major * expected_major,
            "area_m2 should be negligible relative to a major-axis-radius circle for a rank-1 \
             covariance, got area={}, major={}",
            ellipse.area_m2,
            expected_major
        );
        assert!(
            (ellipse.semi_major_m - expected_major).abs() < 1e-9,
            "got {}, expected {}",
            ellipse.semi_major_m,
            expected_major
        );
    }

    /// These are wire-shaped `V1` report types (matching the crate's `solve_json` convention):
    /// the honesty requirement's payload claims are only real if they actually survive a JSON
    /// round trip, not merely present on the in-memory Rust struct. Also pins the wire spelling
    /// of `DifferenceScheme` (snake_case, matching `InputAxis`/`InputGroup`'s existing
    /// convention) and of an unavailable source, since both are newly serializable additions
    /// this task made to the kernel's public types.
    #[test]
    fn the_report_round_trips_through_json_including_unavailable_sources_and_scheme() {
        let r = qnh_resolved();
        let rep = error_budget(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::Altitude, 50.0)],
            &[300.0],
        )
        .unwrap();

        let json = serde_json::to_string(&rep).expect("report must serialize");
        assert!(json.contains("\"unavailable_sources\""));
        assert!(
            json.contains("\"scheme\":\"central\"") || json.contains("\"scheme\": \"central\""),
            "DifferenceScheme must serialize snake_case, got: {json}"
        );

        let round_tripped: ErrorBudgetReportV1 =
            serde_json::from_str(&json).expect("report must deserialize back");
        assert_eq!(round_tripped, rep);
    }

    // ---- Task 11 (MBA-1347) Step 1 tests, verbatim from the brief.

    /// With zero correlation the rectangle probability separates, so we can check the
    /// quadrature against the closed-form product of two normal CDFs.
    #[test]
    fn uncorrelated_rectangle_matches_the_separable_closed_form() {
        use crate::special::normal_cdf;
        let (sd, sw) = (0.10_f64, 0.20_f64);
        let (w, h) = (0.30_f64, 0.40_f64);
        let want = (normal_cdf(h / 2.0 / sd) - normal_cdf(-h / 2.0 / sd))
            * (normal_cdf(w / 2.0 / sw) - normal_cdf(-w / 2.0 / sw));
        let got = p_hit_bivariate(
            sd * sd,
            sw * sw,
            0.0,
            TargetGeometryV1::Rect { width_m: w, height_m: h },
        );
        assert!((got - want).abs() < 1e-6, "got {got} want {want}");
    }

    #[test]
    fn p_hit_is_bounded_and_grows_with_target_size() {
        let small = p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Circle { radius_m: 0.1 });
        let big = p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Circle { radius_m: 0.5 });
        assert!((0.0..=1.0).contains(&small) && (0.0..=1.0).contains(&big));
        assert!(big > small);
    }

    /// Perfecting a source can never reduce hit probability.
    #[test]
    fn perfecting_a_source_never_lowers_p_hit() {
        let r = resolved();
        let rep = error_budget_with_target(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::WindSpeed, 1.5)],
            &[600.0],
            Some(TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.75 }),
        )
        .unwrap();
        for s in &rep.rows[0].sources {
            let gain = s.p_hit_gain_if_perfect.expect("target supplied");
            assert!(gain >= -1e-9, "{:?} reported a negative gain {gain}", s.axis);
        }
    }

    // ---- Task 11: beyond the brief. The brief pins the math and gives exactly the three tests
    // above; the task additionally names four correctness requirements (bounded+monotone,
    // gain-never-negative, degenerate-covariance-no-NaN, unavailable-distinguishable-from-zero)
    // and requires the quadrature itself to be checked against an independent oracle at zero
    // correlation AND shown to differ materially from the (wrong) separable approximation once
    // correlated. Every new public field also gets at least one assertion that would fail if
    // that field were replaced by a constant or a copy of a neighbouring field -- see the task
    // report for the full field-to-test mapping.

    /// Requirement: the quadrature must actually be doing something, not silently degenerating
    /// to the (wrong) separable approximation once the two axes are correlated. `rho = 0.9` on
    /// the SAME (sd, sw, w, h) as the zero-correlation test above -- measured difference from the
    /// wrong separable product is ~0.025, comfortably over the 0.01 threshold asserted here.
    #[test]
    fn a_strongly_correlated_case_differs_materially_from_the_wrong_separable_approximation() {
        let (sd, sw) = (0.10_f64, 0.20_f64);
        let (w, h) = (0.30_f64, 0.40_f64);
        let rho = 0.9_f64;
        let cov = rho * sd * sw;
        let got = p_hit_bivariate(
            sd * sd,
            sw * sw,
            cov,
            TargetGeometryV1::Rect { width_m: w, height_m: h },
        );
        let wrong_separable = (normal_cdf(h / 2.0 / sd) - normal_cdf(-h / 2.0 / sd))
            * (normal_cdf(w / 2.0 / sw) - normal_cdf(-w / 2.0 / sw));
        assert!((0.0..=1.0).contains(&got), "got={got} out of bounds");
        assert!(
            (got - wrong_separable).abs() > 0.01,
            "a correlated case (rho={rho}) should differ materially from the separable \
             approximation, proving the conditional decomposition changes the answer: got={got} \
             wrong_separable={wrong_separable} diff={}",
            (got - wrong_separable).abs()
        );
    }

    /// Requirement 1 (bounded + monotone), rectangle variant -- the brief's own pinned test only
    /// covers the circle branch.
    #[test]
    fn p_hit_grows_with_target_size_for_a_rectangle_too() {
        let small =
            p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Rect { width_m: 0.1, height_m: 0.1 });
        let big =
            p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.5 });
        assert!((0.0..=1.0).contains(&small) && (0.0..=1.0).contains(&big));
        assert!(big > small, "small={small} big={big}");
        // Also pins width_m/height_m against a transposition: growing ONLY the width, or ONLY
        // the height, must each independently grow p_hit.
        let wider =
            p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.1 });
        let taller =
            p_hit_bivariate(0.01, 0.01, 0.0, TargetGeometryV1::Rect { width_m: 0.1, height_m: 0.5 });
        assert!(wider > small, "wider={wider} small={small}");
        assert!(taller > small, "taller={taller} small={small}");
    }

    /// Requirement 3 (degenerate covariance, no NaN), the case beyond what the brief's pinned
    /// `p_hit_bivariate` code actually gets right: `var_drop == 0.0` with `var_wind > 0.0` is a
    /// deterministic drop (always exactly at the nominal point) with real windage uncertainty.
    /// The general quadrature integrates OVER the drop axis and cannot represent a Dirac delta,
    /// so a naive translation that just lets `sd -> 0` zero out the density at every quadrature
    /// node returns a hardcoded `0.0` regardless of target size -- this pins the actual (correct)
    /// closed-form behaviour instead, for both target shapes, and confirms it is NOT stuck at a
    /// constant.
    #[test]
    fn drop_deterministic_windage_random_matches_closed_form_not_hardcoded_zero() {
        let sw = 0.2_f64;
        // Tall enough that the deterministic drop = 0 is always inside the rectangle.
        let target = TargetGeometryV1::Rect { width_m: 0.3, height_m: 10.0 };
        let got = p_hit_bivariate(0.0, sw * sw, 0.0, target);
        let want = normal_cdf(0.15 / sw) - normal_cdf(-0.15 / sw);
        assert!((got - want).abs() < 1e-9, "got={got} want={want}");
        assert!(got > 0.0 && got < 1.0, "fixture must give a non-degenerate probability: {got}");

        let wider = p_hit_bivariate(
            0.0,
            sw * sw,
            0.0,
            TargetGeometryV1::Rect { width_m: 1.0, height_m: 10.0 },
        );
        assert!(wider > got, "a hardcoded-zero bug would make wider == got == 0.0: {wider} {got}");

        let r = 0.15_f64;
        let got_circle = p_hit_bivariate(0.0, sw * sw, 0.0, TargetGeometryV1::Circle { radius_m: r });
        let want_circle = normal_cdf(r / sw) - normal_cdf(-r / sw);
        assert!((got_circle - want_circle).abs() < 1e-9, "got={got_circle} want={want_circle}");
    }

    /// The symmetric degenerate case (`var_wind == 0.0`, `var_drop > 0.0`): this one is already
    /// handled correctly by the general quadrature without a special branch (see
    /// `p_hit_bivariate`'s doc comment), but is pinned here directly for completeness rather than
    /// only indirectly through `error_budget_with_target`.
    #[test]
    fn windage_deterministic_drop_random_matches_closed_form() {
        let sd = 0.15_f64;
        let target = TargetGeometryV1::Rect { width_m: 10.0, height_m: 0.4 };
        let got = p_hit_bivariate(sd * sd, 0.0, 0.0, target);
        let want = normal_cdf(0.2 / sd) - normal_cdf(-0.2 / sd);
        assert!((got - want).abs() < 1e-9, "got={got} want={want}");

        let r = 0.2_f64;
        let got_circle = p_hit_bivariate(sd * sd, 0.0, 0.0, TargetGeometryV1::Circle { radius_m: r });
        let want_circle = normal_cdf(r / sd) - normal_cdf(-r / sd);
        assert!((got_circle - want_circle).abs() < 1e-9, "got={got_circle} want={want_circle}");
        assert!(!got.is_nan() && !got_circle.is_nan());
    }

    /// Requirement 3, continued: `var_drop == var_wind == 0.0` (zero total variance) must give
    /// exactly `1.0` regardless of target size or shape (the nominal point is always the
    /// target's own centre) -- and a nonsensical nonzero `cov` supplied alongside two zero
    /// variances (an inconsistent, invalid covariance a caller should never construct) must not
    /// leak through into a NaN.
    #[test]
    fn zero_total_variance_means_a_deterministic_impact_at_the_nominal_point() {
        for target in [
            TargetGeometryV1::Rect { width_m: 0.001, height_m: 0.001 },
            TargetGeometryV1::Rect { width_m: 5.0, height_m: 5.0 },
            TargetGeometryV1::Circle { radius_m: 0.001 },
            TargetGeometryV1::Circle { radius_m: 5.0 },
        ] {
            let got = p_hit_bivariate(0.0, 0.0, 0.0, target);
            assert_eq!(
                got, 1.0,
                "{target:?}: a deterministic impact at the nominal point (always the target's \
                 own centre here) must be inside with probability exactly 1.0, got {got}"
            );
        }
        let got = p_hit_bivariate(0.0, 0.0, 999.0, TargetGeometryV1::Circle { radius_m: 1.0 });
        assert_eq!(got, 1.0);
        assert!(!got.is_nan());
    }

    /// Requirement 3, continued once more: a SINGLE declared source through the real
    /// `error_budget_with_target` path gives an exact rank-1 covariance (this is the routine
    /// case, not an exception -- see `p_hit_bivariate`'s doc comment). Must not NaN or panic, and
    /// the resulting gain (perfecting the only evaluated source collapses the covariance to
    /// exactly zero) must still be non-negative.
    #[test]
    fn a_single_declared_source_gives_a_well_formed_rank_one_p_hit_not_nan_or_negative() {
        let r = resolved();
        let target = TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.75 };
        let rep =
            error_budget_with_target(&r, &[(InputAxis::WindSpeed, 1.5)], &[600.0], Some(target))
                .unwrap();
        let row = &rep.rows[0];
        let p_hit = row.p_hit.expect("target supplied");
        assert!(!p_hit.is_nan(), "single-source rank-1 covariance must not produce NaN");
        assert!((0.0..=1.0).contains(&p_hit), "p_hit out of bounds: {p_hit}");
        let gain = row.sources[0].p_hit_gain_if_perfect.expect("target supplied");
        assert!(!gain.is_nan());
        assert!(gain >= 0.0, "gain={gain}");
    }

    /// Requirement 4: an unavailable source must remain distinguishable from an evaluated source
    /// that happens to have zero gain -- it has NO `p_hit_gain_if_perfect` at all (a different
    /// TYPE, `UnavailableSourceV1`, with no such field), never a fabricated `Some(0.0)`. The
    /// sibling that DID evaluate must still report a real, independently-checked gain, unaffected
    /// by the other source's unavailability -- extending Task 10's
    /// `an_unavailable_source_is_recorded_not_silently_dropped` to the new target/p_hit
    /// machinery specifically.
    #[test]
    fn an_unavailable_source_has_no_gain_field_while_the_evaluated_sibling_does_when_a_target_is_supplied()
    {
        let r = qnh_resolved();
        let target = TargetGeometryV1::Circle { radius_m: 0.5 };
        let rep = error_budget_with_target(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::Altitude, 50.0)],
            &[300.0],
            Some(target),
        )
        .unwrap();

        assert_eq!(rep.unavailable_sources.len(), 1);
        assert_eq!(rep.unavailable_sources[0].axis, InputAxis::Altitude);
        // `UnavailableSourceV1` has no `p_hit_gain_if_perfect` field at all -- this line would
        // simply fail to COMPILE if one were mistakenly added there, which is a stronger
        // guarantee than any runtime assertion could give.

        assert_eq!(rep.rows[0].sources.len(), 1);
        let mv = &rep.rows[0].sources[0];
        assert_eq!(mv.axis, InputAxis::MuzzleVelocityMps);
        let gain = mv.p_hit_gain_if_perfect.expect("target supplied, and MV evaluated");
        // With exactly one EVALUATED source, excluding it leaves a zero covariance (Altitude
        // never contributes any variance at all, evaluated or not), so the independent oracle
        // for "perfecting the only evaluated source" is exactly `1.0 - the row's own p_hit`.
        let base_p_hit = rep.rows[0].p_hit.expect("target supplied");
        let expected = (1.0 - base_p_hit).max(0.0);
        assert!((gain - expected).abs() < 1e-9, "gain={gain} expected={expected}");
        assert!(
            gain > 0.0,
            "a real muzzle-velocity uncertainty against a finite target should show a strictly \
             positive gain, not merely >= 0: {gain}"
        );
    }

    /// `p_hit`/`p_hit_gain_if_perfect` field pin (would fail if either were replaced by a
    /// constant): `None` on both the row and its sources when no target is supplied, mirroring
    /// `error_budget`'s documented "thin wrapper passing None" contract.
    #[test]
    fn p_hit_fields_are_none_when_no_target_is_supplied() {
        let r = resolved();
        let rep = error_budget(&r, &[(InputAxis::WindSpeed, 1.0)], &[600.0]).unwrap();
        assert_eq!(rep.rows[0].p_hit, None);
        assert_eq!(rep.rows[0].sources[0].p_hit_gain_if_perfect, None);
    }

    /// Field pin for `ErrorBudgetRowV1::p_hit`: checked against an INDEPENDENT
    /// `p_hit_bivariate` call built from this row's own `sigma_drop_m`/`sigma_windage_m`/
    /// `covariance_m2` (not `error_budget_with_target`'s own internal bookkeeping), with the
    /// target deliberately sized relative to this fixture's own dispersion so `p_hit` lands
    /// strictly between 0 and 1 -- neither saturated value would distinguish a real computation
    /// from a `0.0`/`1.0` constant.
    #[test]
    fn p_hit_is_computed_from_the_rows_own_covariance_not_a_constant() {
        let r = resolved();
        let declared = [(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::WindSpeed, 1.0)];
        let probe = error_budget(&r, &declared, &[600.0]).unwrap();
        let row0 = &probe.rows[0];
        let target = TargetGeometryV1::Rect {
            width_m: 2.0 * row0.sigma_windage_m,
            height_m: 2.0 * row0.sigma_drop_m,
        };

        let rep = error_budget_with_target(&r, &declared, &[600.0], Some(target)).unwrap();
        let row = &rep.rows[0];
        let got = row.p_hit.expect("target supplied");

        assert!((0.01..0.99).contains(&got), "fixture should give a non-degenerate p_hit: {got}");
        let oracle = p_hit_bivariate(
            row.sigma_drop_m * row.sigma_drop_m,
            row.sigma_windage_m * row.sigma_windage_m,
            row.covariance_m2,
            target,
        );
        assert!(
            (got - oracle).abs() < 1e-9,
            "p_hit ({got}) does not match an independent p_hit_bivariate call on this row's own \
             covariance ({oracle})"
        );
    }

    /// Field pin for `SourceContributionV1::p_hit_gain_if_perfect`: three declared sources (so
    /// the field is generically discriminating, not degenerate the way it would be with only
    /// two -- same reasoning as `ellipse_area_reduction_is_nonzero_and_discriminating_with_three_or_more_sources`
    /// above), each checked against an independent oracle built from direct `central_difference`
    /// calls summed by hand, never `error_budget_with_target`'s own `accumulate`/`Symmetric2`
    /// path.
    #[test]
    fn p_hit_gain_if_perfect_matches_an_independent_oracle_and_discriminates_with_three_sources() {
        let r = resolved();
        let declared = [
            (InputAxis::MuzzleVelocityMps, 5.0),
            (InputAxis::WindSpeed, 0.3),
            (InputAxis::BallisticCoefficient, 0.001),
        ];
        let probe = error_budget(&r, &declared, &[600.0]).unwrap();
        let row0 = &probe.rows[0];
        let target = TargetGeometryV1::Rect {
            width_m: 2.0 * row0.sigma_windage_m,
            height_m: 2.0 * row0.sigma_drop_m,
        };

        let rep = error_budget_with_target(&r, &declared, &[600.0], Some(target)).unwrap();
        let row = &rep.rows[0];
        assert_eq!(row.sources.len(), 3);
        let base_p_hit = row.p_hit.expect("target supplied");

        let mut derivs = std::collections::HashMap::new();
        for &(axis, sigma) in &declared {
            let d = central_difference(&r, axis, &[600.0], None).unwrap()[0];
            derivs.insert(axis, (sigma, d.d_drop_d_x, d.d_windage_d_x));
        }
        let variance_excluding = |exclude: Option<InputAxis>| -> (f64, f64, f64) {
            let mut vd = 0.0;
            let mut vw = 0.0;
            let mut cv = 0.0;
            for (&axis, &(sigma, dd, dw)) in &derivs {
                if Some(axis) == exclude {
                    continue;
                }
                let s2 = sigma * sigma;
                vd += dd * dd * s2;
                vw += dw * dw * s2;
                cv += dd * dw * s2;
            }
            (vd, vw, cv)
        };
        let (fvd, fvw, fcv) = variance_excluding(None);
        let oracle_base_p_hit = p_hit_bivariate(fvd, fvw, fcv, target);
        assert!((base_p_hit - oracle_base_p_hit).abs() < 1e-9);

        let mut gains = Vec::new();
        for s in &row.sources {
            let (vd, vw, cv) = variance_excluding(Some(s.axis));
            let expected = (p_hit_bivariate(vd, vw, cv, target) - oracle_base_p_hit).max(0.0);
            let got = s.p_hit_gain_if_perfect.expect("target supplied");
            assert!(
                (got - expected).abs() < 1e-9,
                "{:?}: got {got}, independently expected {expected}",
                s.axis
            );
            gains.push((s.axis, got));
        }
        assert!(
            gains.iter().any(|&(_, g)| g > 0.0),
            "at least one source should show a positive gain: {gains:?}"
        );
        let first = gains[0].1;
        assert!(
            gains.iter().any(|&(_, g)| (g - first).abs() > 1e-6 * first.max(1.0)),
            "gains must discriminate between sources, not all be equal: {gains:?}"
        );
    }

    /// `method`/`assumptions` field pin: the `"_gl20_pm6sigma"` suffix and the new hit-probability
    /// assumption must appear ONLY when a target is supplied, never when it is not (which would
    /// also break `the_report_declares_independence_and_linearity`'s exact-equality check on
    /// `method`).
    #[test]
    fn the_report_names_the_quadrature_and_the_aim_point_assumption_only_when_a_target_is_supplied()
    {
        let r = resolved();
        let sources = [(InputAxis::WindSpeed, 1.0)];

        let without = error_budget_with_target(&r, &sources, &[600.0], None).unwrap();
        assert_eq!(without.method, "central_difference_first_order_propagation");
        assert!(!without.assumptions.iter().any(|s| s.to_lowercase().contains("gauss-legendre")));

        let with_target = error_budget_with_target(
            &r,
            &sources,
            &[600.0],
            Some(TargetGeometryV1::Circle { radius_m: 0.3 }),
        )
        .unwrap();
        assert_eq!(
            with_target.method,
            "central_difference_first_order_propagation_gl20_pm6sigma"
        );
        assert!(
            with_target
                .assumptions
                .iter()
                .any(|s| s.to_lowercase().contains("gauss-legendre") && s.contains("aim point")),
            "{:#?}",
            with_target.assumptions
        );
    }

    /// Wire-shape pin for the new `TargetGeometryV1` type: externally-tagged, snake_case variant
    /// names, matching this module's existing `DifferenceScheme`/`UnavailableReasonCodeV1`
    /// convention.
    #[test]
    fn target_geometry_serializes_snake_case_externally_tagged() {
        let rect = TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.75 };
        let json = serde_json::to_string(&rect).unwrap();
        assert!(json.contains("\"rect\""), "{json}");
        assert!(json.contains("\"width_m\":0.5") || json.contains("\"width_m\": 0.5"), "{json}");
        assert!(
            json.contains("\"height_m\":0.75") || json.contains("\"height_m\": 0.75"),
            "{json}"
        );
        let back: TargetGeometryV1 = serde_json::from_str(&json).expect("must deserialize");
        assert_eq!(back, rect);

        let circle = TargetGeometryV1::Circle { radius_m: 0.3 };
        let json2 = serde_json::to_string(&circle).unwrap();
        assert!(json2.contains("\"circle\""), "{json2}");
        assert!(
            json2.contains("\"radius_m\":0.3") || json2.contains("\"radius_m\": 0.3"),
            "{json2}"
        );
        let back2: TargetGeometryV1 = serde_json::from_str(&json2).expect("must deserialize");
        assert_eq!(back2, circle);
    }

    /// The honesty requirement's payload claims for hit probability are only real if `p_hit`/
    /// `p_hit_gain_if_perfect` survive a JSON round trip as actual present (`Some`) values, not
    /// only the already-covered `None` case in
    /// `the_report_round_trips_through_json_including_unavailable_sources_and_scheme` above.
    ///
    /// Deliberately does NOT assert full-struct `assert_eq!(round_tripped, rep)` the way that
    /// other test does: this fixture's `WindSpeed` source has a `d_drop_d_x` of
    /// `3.7609987435516246e-6` (a near-zero value -- wind primarily affects windage, not drop, so
    /// a small residual drop sensitivity is ordinary, not a red flag), and serde_json 1.0.149's
    /// float parser round-trips that SPECIFIC value one ULP off (`3.7609987435516246e-6`
    /// serializes correctly via `ryu`,
    /// but re-parsing that exact string back yields bit pattern `...ece00001` instead of the
    /// original `...ece00000`, confirmed with a standalone minimal reproduction outside this
    /// crate) -- a pre-existing float-formatting edge case in a dependency, unrelated to this
    /// task and already latent in `d_drop_d_x` before this task added anything. A bit-exact
    /// struct comparison here would make THIS test flaky on a fact that has nothing to do with
    /// what it is actually checking, so it compares only the fields it cares about, with a
    /// tolerance many orders of magnitude looser than one ULP.
    #[test]
    fn p_hit_and_gain_round_trip_through_json_when_present() {
        let r = resolved();
        let target = TargetGeometryV1::Rect { width_m: 0.5, height_m: 0.75 };
        let rep = error_budget_with_target(
            &r,
            &[(InputAxis::MuzzleVelocityMps, 5.0), (InputAxis::WindSpeed, 1.5)],
            &[600.0],
            Some(target),
        )
        .unwrap();

        let json = serde_json::to_string(&rep).expect("report must serialize");
        assert!(json.contains("\"p_hit\":"), "{json}");
        assert!(json.contains("\"p_hit_gain_if_perfect\":"), "{json}");

        let round_tripped: ErrorBudgetReportV1 =
            serde_json::from_str(&json).expect("report must deserialize back");

        let want_p_hit = rep.rows[0].p_hit.expect("target supplied");
        let got_p_hit = round_tripped.rows[0].p_hit.expect("must round-trip as Some");
        assert!((got_p_hit - want_p_hit).abs() < 1e-9, "got={got_p_hit} want={want_p_hit}");

        assert_eq!(round_tripped.rows[0].sources.len(), rep.rows[0].sources.len());
        for (got_s, want_s) in
            round_tripped.rows[0].sources.iter().zip(rep.rows[0].sources.iter())
        {
            assert_eq!(got_s.axis, want_s.axis);
            let want_gain = want_s.p_hit_gain_if_perfect.expect("target supplied");
            let got_gain = got_s.p_hit_gain_if_perfect.expect("must round-trip as Some");
            assert!(
                (got_gain - want_gain).abs() < 1e-9,
                "{:?}: got={got_gain} want={want_gain}",
                got_s.axis
            );
        }
    }
}
