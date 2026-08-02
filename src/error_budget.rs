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
//! [`KernelError::StepOutOfDomain`] (both perturbed sides left the axis's physical domain, e.g.
//! a declared sigma larger than the distance from the nominal value to a hard boundary). Every
//! one of these is recorded as an
//! [`UnavailableSourceV1`] (axis, declared sigma, and a human-readable reason) in
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
//! happens to be.
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

use serde::{Deserialize, Serialize};

use crate::perturbation::access::KernelError;
use crate::perturbation::derive::{central_difference, DifferenceScheme, Derivative};
use crate::perturbation::taxonomy::InputAxis;
use crate::solve_json::ResolvedSolveRequestV1;
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
    pub ellipse_area_reduction_m2: f64,
}

/// A declared source [`error_budget`] could not evaluate for this request, and why.
///
/// See this module's "Unavailable sources" doc section. An axis appearing here never also
/// appears in any row's `sources` -- the two are disjoint by construction. This list is the same
/// for every row: whether `central_difference` can differentiate a given axis on a given request
/// does not depend on which ranges were requested (see this module's "Cost" section).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnavailableSourceV1 {
    pub axis: InputAxis,
    pub sigma: f64,
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
/// of the report should still be produced" (`Some(reason)`) or "this is a genuine failure the
/// caller must see" (`None`).
///
/// `Some` for the four structural refusals this module's doc names explicitly
/// (`AxisUnsupportedForRequest`, `AxisAbsent`, `CategoricalAxis`, `StepOutOfDomain`); `None` for
/// everything else (`Solve`, `Observation`, `TypeMismatch`, `NonFinite`), which the caller must
/// propagate. Exhaustive over every `KernelError` variant with no wildcard arm, so a future
/// variant fails to compile here until it is explicitly placed in one bucket or the other.
fn unavailable_reason(e: &KernelError) -> Option<String> {
    match e {
        KernelError::AxisUnsupportedForRequest { reason, .. } => Some(reason.to_string()),
        KernelError::AxisAbsent(_) => Some(
            "this axis has no single scalar value on this request (for the three wind axes, \
             this means the wind is declared as a segmented profile rather than a constant \
             speed/direction)"
                .to_string(),
        ),
        KernelError::CategoricalAxis(_) => Some(
            "this axis is categorical (a toggle or an enumerated choice), not a continuous \
             quantity, and cannot be assigned a one-sigma uncertainty or differentiated"
                .to_string(),
        ),
        KernelError::StepOutOfDomain { attempted, .. } => Some(format!(
            "central differencing needs to perturb this axis by {attempted:.6} in both \
             directions from its nominal value, and both directions left its physical domain; \
             its sensitivity could not be measured at this operating point"
        )),
        KernelError::Solve { .. }
        | KernelError::Observation(_)
        | KernelError::TypeMismatch(_)
        | KernelError::NonFinite(_) => None,
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
fn sort_by_variance_share_desc(sources: &mut [SourceContributionV1]) {
    sources.sort_by(|x, y| {
        y.variance_share
            .partial_cmp(&x.variance_share)
            .expect("variance_share is always finite and non-negative by construction")
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
                "{:?} dominates at {range_m:.0} m ({:.0}% of impact variance){caveat}. \
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

/// Propagate each declared per-input uncertainty in `sources` to impact covariance at every
/// range in `ranges_m`, via central differences through the real solver
/// ([`central_difference`]), and rank the sources by their share of impact variance.
///
/// `sources` is `(axis, sigma)` pairs: `sigma` is the caller's one-sigma uncertainty for that
/// axis, in the axis's own physical unit (`crate::perturbation::axis_meta(axis).kind`).
/// Declaring the same axis twice is neither
/// detected nor rejected -- the two entries are treated as fully independent sources, so their
/// variance would be double-counted, which is not physically meaningful. Callers must not do
/// this.
///
/// # Honesty
///
/// [`ErrorBudgetReportV1::method`] and [`ErrorBudgetReportV1::assumptions`] state, in the
/// payload itself (not only in prose documentation), that: sources are treated as INDEPENDENT
/// (no correlation between them is modelled); propagation is FIRST-ORDER/local-linear about the
/// nominal solution (not exact for large or non-Gaussian input uncertainty); the 95% ellipse
/// assumes an approximately Gaussian impact distribution; a source's derivative may be one-sided
/// rather than central; and an unavailable source is not the same fact as a zero-contribution
/// one. See `the_report_declares_independence_and_linearity` in this module's tests.
///
/// # Unavailable sources
///
/// See this module's top-level "Unavailable sources" doc section.
///
/// # Errors
///
/// Propagates any [`KernelError`] from [`central_difference`] that is not one of the four
/// structural refusals recorded in [`ErrorBudgetReportV1::unavailable_sources`] instead -- a
/// genuine solver or trajectory failure, not a normal "this input cannot be perturbed here"
/// fact. Also returns [`KernelError::NonFinite`] immediately, before any solve is attempted, if
/// any declared `sigma` is not finite or is negative.
pub fn error_budget(
    base: &ResolvedSolveRequestV1,
    sources: &[(InputAxis, f64)],
    ranges_m: &[f64],
) -> Result<ErrorBudgetReportV1, KernelError> {
    let mut jac: Vec<AxisJacobian> = Vec::with_capacity(sources.len());
    let mut unavailable: Vec<UnavailableSourceV1> = Vec::new();

    for &(axis, sigma) in sources {
        if !(sigma.is_finite() && sigma >= 0.0) {
            return Err(KernelError::NonFinite(axis));
        }
        // One central_difference call per DECLARED source, covering every range in ranges_m at
        // once -- never re-derived per range. See this module's "Cost" doc section.
        match central_difference(base, axis, ranges_m, None) {
            Ok(derivatives) => jac.push(AxisJacobian { axis, sigma, derivatives }),
            Err(e) => match unavailable_reason(&e) {
                Some(reason) => unavailable.push(UnavailableSourceV1 { axis, sigma, reason }),
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

        let mut sources_out: Vec<SourceContributionV1> = entries
            .iter()
            .map(|e| {
                let s2 = e.sigma * e.sigma;
                let this_var = e.d_drop_d_x * e.d_drop_d_x * s2 + e.d_windage_d_x * e.d_windage_d_x * s2;
                let reduced_ellipse = ellipse_95(accumulate(&entries, Some(e.axis)));
                SourceContributionV1 {
                    axis: e.axis,
                    sigma: e.sigma,
                    d_drop_d_x: e.d_drop_d_x,
                    d_windage_d_x: e.d_windage_d_x,
                    scheme: e.scheme,
                    variance_share: if total_var > 0.0 { this_var / total_var } else { 0.0 },
                    ellipse_area_reduction_m2: (full_ellipse.area_m2 - reduced_ellipse.area_m2)
                        .max(0.0),
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
            sources: sources_out,
            priority_statement,
        });
    }

    Ok(ErrorBudgetReportV1 {
        schema_version: ERROR_BUDGET_SCHEMA_VERSION_V1,
        method: "central_difference_first_order_propagation".to_string(),
        assumptions: vec![
            "Declared sources are treated as independent; correlations between them are not \
             modelled."
                .to_string(),
            "Propagation is first-order (local linear) about the nominal solution, evaluated by \
             central differences through the real solver. It is not exact for large or \
             non-Gaussian input uncertainty."
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
        ],
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
    #[test]
    fn unavailable_reason_classifies_every_kernel_error_variant() {
        use crate::solve_json::SolveErrorCodeV1;
        use crate::trajectory_observation::TrajectoryObservationError;

        let structural = [
            KernelError::AxisUnsupportedForRequest { axis: InputAxis::Altitude, reason: "x" },
            KernelError::AxisAbsent(InputAxis::WindSpeed),
            KernelError::CategoricalAxis(InputAxis::CoriolisEnabled),
            KernelError::StepOutOfDomain { axis: InputAxis::RelativeHumidity, attempted: 2.0 },
        ];
        for e in &structural {
            assert!(
                unavailable_reason(e).is_some(),
                "{e:?} must be classified as unavailable (recorded), not propagated"
            );
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

    /// Declaring a non-finite or negative sigma is rejected up front, before any solve --
    /// otherwise a NaN sigma would silently poison `variance_share`'s comparator in the sort
    /// (`partial_cmp` returns `None` for NaN, and the implementation unwraps it).
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
    #[test]
    fn single_source_ellipse_matches_the_closed_form_rank_one_case() {
        let r = resolved();
        let sigma = 1.0_f64;
        let rep = error_budget(&r, &[(InputAxis::WindSpeed, sigma)], &[600.0]).unwrap();
        let d = central_difference(&r, InputAxis::WindSpeed, &[600.0], None).unwrap()[0];

        let expected_largest = sigma * sigma * (d.d_drop_d_x.powi(2) + d.d_windage_d_x.powi(2));
        // The minor axis is exactly zero for a rank-1 covariance, so the area (semi_major *
        // semi_minor * pi) is exactly zero too, regardless of how large the major axis is.
        assert_eq!(rep.rows[0].ellipse_95.semi_minor_m, 0.0);
        assert_eq!(rep.rows[0].ellipse_95.area_m2, 0.0);
        let expected_major = (CHI2_95_2DOF * expected_largest).sqrt();
        assert!(
            (rep.rows[0].ellipse_95.semi_major_m - expected_major).abs() < 1e-9,
            "got {}, expected {}",
            rep.rows[0].ellipse_95.semi_major_m,
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
}
