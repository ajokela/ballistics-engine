//! Shared display-ready row type behind the CLI's card-shaped table/CSV/JSON surfaces
//! (0.33.0 decision-support Plan B Task 9).
//!
//! `come-ups`, `range-table`, `wind-card` and `compare` each grew their own
//! function-local row struct (`ComeUpRow`/`RangeRow`/`WindRow`/`LoadRow`) that said the
//! same handful of things -- range, drop, wind, velocity, energy, time -- a different way,
//! which blocked any shared card machinery between them. `CardRow` is the display-ready
//! superset: every surface populates only the fields it has ever had and leaves the rest
//! `None` / empty, so each surface's existing table/CSV/JSON writer keeps reading the
//! identical numbers it always did (pinned byte-identical by `tests/card_golden_cli.rs`).
//!
//! No feature gate: this module must compile for `wasm32-unknown-unknown` with no default
//! features (pure data, no `fs`, no `clap`, no `pdf`). Task 10 rewrites the PDF dope card on
//! `&[CardRow]`; Task 11 grows an adaptive-card engine in this module.

/// One card row. Display-ready values in the surface's chosen units (exactly what the
/// legacy per-surface structs stored), so rendering is unchanged; range is f64 metres-
/// or-display per the surface's existing convention — DO NOT re-convert anything.
#[derive(Debug, Clone)]
pub struct CardRow {
    pub range: f64,
    pub drop_linear: Option<f64>,
    pub drop_adj: Option<f64>,
    pub come_up: Option<f64>,
    pub wind_linear: Option<f64>,
    pub wind_adj: Option<f64>,
    pub velocity: Option<f64>,
    pub energy: Option<f64>,
    pub time: Option<f64>,
    pub lead_adj: Option<f64>,
    /// wind-card's per-speed drift columns; empty elsewhere.
    pub wind_columns: Vec<f64>,
}

// ---------------------------------------------------------------------------
// Adaptive range-card engine (0.33.0 decision-support Plan B Task 11, MBA-1351).
//
// `use` items sit here rather than at the top of the file so this task's addition is a
// pure append -- Task 9's `CardRow` block above is untouched, and a parallel lane editing
// it merges cleanly. Rust does not care where module-level items appear.
// ---------------------------------------------------------------------------

use crate::adjustment::{click_size_mil, quantize_angle, ClickBase, ClickValue};
use crate::hold_curve::HoldCurve;
use std::fmt;

/// Schema version of [`AdaptiveCardReportV1`]. Bump only for a breaking shape change.
pub const ADAPTIVE_CARD_SCHEMA_VERSION_V1: u32 = 1;

/// The reconstruction error a shooter is willing to accept, per axis, **in the card's
/// printed adjustment unit** (mil for [`CardAdjustmentUnit::Mil`], the locked-3438 MOA for
/// [`CardAdjustmentUnit::Moa`]). The caller converts; this engine never guesses a unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveBudget {
    pub elevation: f64,
    pub windage: f64,
}

/// Everything one adaptive card needs beyond the solved curve itself.
///
/// `elevation_cf` / `windage_cf` are the scope's tracking correction factors (MBA-1358).
/// Both are VALIDATED here against [`crate::adjustment::tracking_cf_in_range`]'s locked
/// `(0.5, 1.5)` band and rejected with [`CardError::InvalidTrackingCf`] -- this is a public
/// library API that language bindings call without the CLI's own validation, and an
/// out-of-band CF fails silently rather than loudly. `bias_mil` is the selected zero set's elevation dial
/// correction (MBA-1360) in true angular mil; it applies to the ELEVATION axis only, which
/// is what "zero-set bias as a drop-equivalent" means.
#[derive(Debug, Clone)]
pub struct AdaptiveRequest<'a> {
    /// `(start, end)` in meters. Both must be finite, positive, and increasing.
    pub domain_m: (f64, f64),
    /// Ranges that must appear as rows whatever the error says (a known dope point, a
    /// target distance). Validated into the domain, never silently dropped.
    pub anchors_m: Vec<f64>,
    pub budget: AdaptiveBudget,
    /// Upper bound on printed rows. The mandatory seed (both domain ends plus every anchor)
    /// is never truncated to honour it -- see [`AdaptiveCardReportV1::rows_capped`].
    pub max_rows: usize,
    /// `(elevation, windage)` turret graduations. `Some` snaps every printed row onto the
    /// detent lattice, which is what a shooter can actually dial; `None` prints the
    /// unrounded angle.
    pub click: Option<(&'a ClickValue, &'a ClickValue)>,
    pub elevation_cf: f64,
    pub windage_cf: f64,
    pub bias_mil: f64,
}

/// The printed adjustment unit of a card's dial columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CardAdjustmentUnit {
    Mil,
    Moa,
}

impl CardAdjustmentUnit {
    /// Multiplier from the curve's true angular mil into this printed unit.
    ///
    /// `Moa` uses this crate's LOCKED printed-table constant `3438` (MBA-724), deliberately
    /// not the exact-angle 3437.7467 -- every printed MOA column in the crate is drawn on
    /// that ratio, and an adaptive card that measured its error on a different one would be
    /// auditing numbers nobody prints.
    pub fn from_mil_factor(&self) -> f64 {
        match self {
            Self::Mil => 1.0,
            Self::Moa => 3438.0 / 1000.0,
        }
    }

    /// The [`ClickBase`] this printed unit quantizes on, so the synthetic printed-space
    /// graduation handed to [`quantize_angle`] is self-consistent rather than mislabelled.
    fn click_base(&self) -> ClickBase {
        match self {
            Self::Mil => ClickBase::Mil,
            Self::Moa => ClickBase::Moa,
        }
    }
}

/// Every way an [`adaptive_card`] request can be rejected. All five are structured: a range
/// a shooter asked for and did not get back is information, never a silent drop.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CardError {
    /// The domain is not a forward, positive, finite interval. (Angular drop divides by the
    /// range, so a start at or below zero has no angular value at all.)
    EmptyOrInvertedDomain { start_m: f64, end_m: f64 },
    AnchorOutsideDomain {
        anchor_m: f64,
        start_m: f64,
        end_m: f64,
    },
    NonPositiveBudget { axis: &'static str, value: f64 },
    ZeroMaxRows,
    /// The domain runs past the last sampled point of the curve, where there is no ground
    /// truth to verify against.
    DomainOutsideCurve { requested_m: f64, curve_max_m: f64 },
    /// A tracking correction factor outside [`crate::adjustment::tracking_cf_in_range`]'s
    /// locked `(0.5, 1.5)` band (MBA-1358).
    ///
    /// Checked rather than assumed because the failure is SILENT AND CONFIDENT, not loud: a
    /// finite but out-of-band CF -- the realistic slip of passing the percentage `95` where
    /// the ratio `0.95` belongs, or an `INFINITY` -- divides every printed value to ~0, so
    /// every measured error is ~0 and the engine would hand back a two-row card of near-zero
    /// dial values reporting `budget_met: true`. That is a wrong answer on the one field
    /// this whole module exists to make trustworthy, and it is worse than the NaN a zero CF
    /// produces, which at least reports `budget_met: false`.
    InvalidTrackingCf { axis: &'static str, value: f64 },
}

impl fmt::Display for CardError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyOrInvertedDomain { start_m, end_m } => write!(
                f,
                "card domain {start_m} m to {end_m} m must be a forward interval with a positive start"
            ),
            Self::AnchorOutsideDomain {
                anchor_m,
                start_m,
                end_m,
            } => write!(
                f,
                "anchor {anchor_m} m lies outside the card domain {start_m} m to {end_m} m"
            ),
            Self::NonPositiveBudget { axis, value } => {
                write!(f, "{axis} budget {value} must be positive and finite")
            }
            Self::ZeroMaxRows => write!(f, "a card needs room for at least one row"),
            Self::DomainOutsideCurve {
                requested_m,
                curve_max_m,
            } => write!(
                f,
                "range {requested_m} m is past the curve's last sampled point at {curve_max_m} m"
            ),
            Self::InvalidTrackingCf { axis, value } => write!(
                f,
                "{axis} tracking correction factor {value} must be finite and between 0.5 and 1.5 \
                 (it is a ratio such as 0.95, not a percentage)"
            ),
        }
    }
}

impl std::error::Error for CardError {}

/// What an adaptive card is, plus what its own numbers were measured to be worth.
///
/// The error fields are MEASURED by a dense verification pass over the declared grid after
/// the rows are final -- they are not the loop's running estimate, and `budget_met` is
/// recomputed from that pass rather than inherited from the insertion loop.
#[derive(Debug, Clone)]
pub struct AdaptiveCardReportV1 {
    pub schema_version: u32,
    /// Stable identifier for how these rows were chosen.
    pub method: String,
    /// Exactly five entries; see [`adaptive_card`]. Index-stable: a consumer may quote
    /// `assumptions[3]` and mean the half-click floor.
    pub assumptions: Vec<String>,
    /// Printed rows, ascending by range (meters). Dial values are quantized when clicks
    /// were supplied.
    pub rows: Vec<CardRow>,
    pub budget_met: bool,
    /// The row cap stopped refinement while violations remained.
    pub rows_capped: bool,
    /// Worst measured elevation error over the audited points, in the card's printed unit.
    pub worst_elevation_error: f64,
    /// Worst measured windage error over the audited points, in the card's printed unit.
    pub worst_windage_error: f64,
    /// Range of the single worst point by budget-normalized excess (ties to the lowest
    /// range). The two per-axis maxima above may each occur elsewhere, which is exactly why
    /// they are reported as their own scalars.
    pub worst_error_range_m: f64,
    /// Spacing of the declared verification grid, meters. The honesty claim extends to this
    /// grid and no further.
    pub verification_grid_step_m: f64,
}

/// One axis's printed-value pipeline, in the LOCKED composition order the CLI's
/// `adjustment_display` boundary already uses (MBA-1360 x MBA-1358): the zero-set bias joins
/// the TRUE angular need first, the tracking correction divides second, click quantization
/// happens last on the corrected value. The order is load-bearing and must not move.
#[derive(Debug, Clone, Copy)]
struct PrintedAxis {
    unit_factor: f64,
    bias_mil: f64,
    cf: f64,
    /// The detent graduation expressed in the PRINTED unit, so quantization and the error
    /// metric live in the same space. `None` prints the unrounded angle.
    click: Option<ClickValue>,
}

impl PrintedAxis {
    fn new(unit: CardAdjustmentUnit, bias_mil: f64, cf: f64, click: Option<&ClickValue>) -> Self {
        let unit_factor = unit.from_mil_factor();
        Self {
            unit_factor,
            bias_mil,
            cf,
            // click_size_mil is the crate's one click -> mil converter (locked 3438 for MOA);
            // scaling its result by the printed unit factor lands the graduation in printed
            // space without a second, divergent conversion table.
            click: click.map(|c| ClickValue {
                size: click_size_mil(c) * unit_factor,
                base: unit.click_base(),
            }),
        }
    }

    /// The unquantized printed value -- what the card WOULD say with infinite dial
    /// resolution. This is the ground truth every error below is measured against.
    fn exact(&self, true_mil: f64) -> f64 {
        let printed = true_mil * self.unit_factor;
        // Skipping a zero bias is bit-exact on purpose: an unconditional `+ 0.0` flips a
        // -0.0 to +0.0 and changes rendered bytes (the same rule `adjustment_display` keeps).
        let biased = if self.bias_mil != 0.0 {
            printed + self.bias_mil * self.unit_factor
        } else {
            printed
        };
        biased / self.cf
    }

    /// The value actually PRINTED on the row: [`Self::exact`] snapped to the detent lattice
    /// when the optic's clicks are known.
    fn printed(&self, true_mil: f64) -> f64 {
        let exact = self.exact(true_mil);
        match &self.click {
            Some(c) => quantize_angle(exact, c).clicks as f64 * c.size,
            None => exact,
        }
    }

    /// Half the detent spacing: the error a quantized row carries AT ITS OWN RANGE, which no
    /// number of extra rows can reduce.
    fn half_click_floor(&self) -> f64 {
        self.click.map_or(0.0, |c| c.size / 2.0)
    }
}

/// One audited range: the exact printed values, the values the card prints there, and the
/// display extras a [`CardRow`] wants. Precomputed once -- the curve is never re-queried
/// inside the insertion loop, which keeps each pass O(points) and makes the loop's result
/// independent of how many passes it takes.
#[derive(Debug, Clone, Copy)]
struct AuditPoint {
    range_m: f64,
    exact_elevation: f64,
    exact_windage: f64,
    printed_elevation: f64,
    printed_windage: f64,
    drop_linear_m: f64,
    wind_linear_m: f64,
    velocity_mps: f64,
    energy_j: f64,
    time_s: f64,
}

/// Per-axis absolute reconstruction error at one audited point.
type AxisErrors = (f64, f64);

/// Sort ascending and drop exact duplicates. `total_cmp` gives a total order with no
/// comparator-contract hazard; NaN cannot reach here (validation rejects it).
fn sorted_dedup(mut values: Vec<f64>) -> Vec<f64> {
    values.sort_by(f64::total_cmp);
    values.dedup();
    values
}

/// The hold curve's own sample points inside `[start_m, end_m]`.
///
/// [`HoldCurve`] samples at exact multiples of [`HoldCurve::SAMPLE_INTERVAL_M`] (the sampler
/// builds its distance list as `i as f64 * step`), so recomputing that arithmetic sequence
/// reproduces the native grid bit-for-bit without needing access to the curve's private
/// sample vector -- pinned by `verification_grid_lands_on_the_curves_native_samples`.
fn native_grid_m(start_m: f64, end_m: f64) -> Vec<f64> {
    let step = HoldCurve::SAMPLE_INTERVAL_M;
    // Index 0 is the muzzle, where an angular hold is undefined; start at 1.
    let first = ((start_m / step).ceil() as i64).max(1);
    let last = (end_m / step).floor() as i64;
    let mut grid = Vec::new();
    for i in first..=last {
        let g = i as f64 * step;
        // Re-check the bounds rather than trusting ceil/floor at the endpoints.
        if g >= start_m && g <= end_m {
            grid.push(g);
        }
    }
    grid
}

/// Reconstruct the printed card at every audited point and measure the per-axis error.
///
/// At a row the reconstruction IS the row's printed value, so the error there is the
/// quantization residual -- not zero. That is deliberate: hiding it would hide the
/// half-click floor, the one error extra rows cannot fix.
fn sweep(audit: &[AuditPoint], rows: &[usize]) -> Vec<AxisErrors> {
    let mut errors = vec![(0.0, 0.0); audit.len()];
    for &r in rows {
        let p = &audit[r];
        errors[r] = (
            (p.printed_elevation - p.exact_elevation).abs(),
            (p.printed_windage - p.exact_windage).abs(),
        );
    }
    for pair in rows.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let (a, b) = (&audit[lo], &audit[hi]);
        let span = b.range_m - a.range_m;
        for (k, point) in audit.iter().enumerate().take(hi).skip(lo + 1) {
            let t = if span > 0.0 {
                (point.range_m - a.range_m) / span
            } else {
                0.0
            };
            let elevation = a.printed_elevation + (b.printed_elevation - a.printed_elevation) * t;
            let windage = a.printed_windage + (b.printed_windage - a.printed_windage) * t;
            errors[k] = (
                (elevation - point.exact_elevation).abs(),
                (windage - point.exact_windage).abs(),
            );
        }
    }
    errors
}

/// Loop instrumentation, for the termination tests only. Not part of the public report:
/// how many passes the search took is an implementation detail, but it is the ONE
/// observable that distinguishes "stopped because the error is irreducible" from "spun
/// until the runaway backstop caught it".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct LoopTrace {
    iterations: usize,
    iteration_cap: usize,
}

/// Build an adaptive range card and measure what it is worth.
///
/// Greedy worst-point insertion: start from the domain ends plus every anchor, then
/// repeatedly add the single audited point whose linearly-interpolated printed value is
/// furthest outside budget, until nothing violates, the row cap binds, or the remaining
/// error is irreducible. A SEPARATE dense pass over the declared grid then measures the
/// finished card, and `budget_met` comes from that measurement.
///
/// Everything is measured in printed-value space -- the same zero-set bias, tracking
/// correction and click quantization the rows carry -- so the reported error is the error a
/// shooter interpolating the printed card actually makes.
///
/// # What this does and does not buy you
///
/// It buys a MEASURED error bound, guaranteed anchors, and no step to guess. It does not
/// reliably buy a shorter card than a well-chosen fixed step: a single insertion can at
/// best quarter an interval's error (a bisection), so on a trajectory whose curvature
/// barely varies the whole card doubles at once while a uniform card may pick any row
/// count. See `fixed_step_comparison_is_measured_not_assumed` for the measurements.
///
/// # Errors
///
/// Returns [`CardError`] for an inverted or non-positive domain, an anchor outside it, a
/// non-positive budget, a zero row cap, a domain running past the curve's last sample, or a
/// tracking correction factor outside the locked `(0.5, 1.5)` band. Every one of these is
/// checked before any work is done, on every build profile -- a caller cannot reach the
/// solver with a request that would produce a confidently wrong card.
pub fn adaptive_card(
    curve: &HoldCurve,
    req: &AdaptiveRequest,
    unit: CardAdjustmentUnit,
) -> Result<AdaptiveCardReportV1, CardError> {
    adaptive_card_traced(curve, req, unit).map(|(report, _)| report)
}

/// [`adaptive_card`] plus its loop trace. Private: see [`LoopTrace`].
fn adaptive_card_traced(
    curve: &HoldCurve,
    req: &AdaptiveRequest,
    unit: CardAdjustmentUnit,
) -> Result<(AdaptiveCardReportV1, LoopTrace), CardError> {
    let (start_m, end_m) = req.domain_m;

    if req.max_rows == 0 {
        return Err(CardError::ZeroMaxRows);
    }
    for (axis, value) in [
        ("elevation", req.budget.elevation),
        ("windage", req.budget.windage),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(CardError::NonPositiveBudget { axis, value });
        }
    }
    // Enforced, not assumed: an out-of-band CF silently produces a confident, wrong card
    // (see `CardError::InvalidTrackingCf`), and a `debug_assert` is compiled out of exactly
    // the builds that ship. `tracking_cf_in_range` is the crate's ONE locked band (MBA-1358),
    // shared with the CLI and the WASM terminal -- reused here, never restated as a literal.
    for (axis, value) in [
        ("elevation", req.elevation_cf),
        ("windage", req.windage_cf),
    ] {
        if !crate::adjustment::tracking_cf_in_range(value) {
            return Err(CardError::InvalidTrackingCf { axis, value });
        }
    }
    if !start_m.is_finite() || !end_m.is_finite() || start_m <= 0.0 || end_m <= start_m {
        return Err(CardError::EmptyOrInvertedDomain { start_m, end_m });
    }
    let curve_max_m = curve.max_sampled_range_m();
    if end_m > curve_max_m {
        return Err(CardError::DomainOutsideCurve {
            requested_m: end_m,
            curve_max_m,
        });
    }
    for &anchor_m in &req.anchors_m {
        if !anchor_m.is_finite() || anchor_m < start_m || anchor_m > end_m {
            return Err(CardError::AnchorOutsideDomain {
                anchor_m,
                start_m,
                end_m,
            });
        }
    }
    debug_assert!(req.bias_mil.is_finite(), "zero-set bias must be finite");

    let (elevation_click, windage_click) = match req.click {
        Some((e, w)) => (Some(e), Some(w)),
        None => (None, None),
    };
    let elevation = PrintedAxis::new(unit, req.bias_mil, req.elevation_cf, elevation_click);
    // The zero-set bias is an elevation-only dial correction; windage carries its own
    // tracking correction but no bias on this interface.
    let windage = PrintedAxis::new(unit, 0.0, req.windage_cf, windage_click);

    // Audited points = the curve's native grid inside the domain, UNION the mandatory rows.
    // The union matters: a domain end or an anchor need not land on the grid, and a
    // quantized row carries error at its own range, so leaving one unaudited would under-
    // report the very floor `assumptions[3]` warns about. Every point the loop may insert
    // therefore already lives in this fixed, precomputed set.
    let mut seeds = vec![start_m, end_m];
    seeds.extend_from_slice(&req.anchors_m);
    let mut ranges = native_grid_m(start_m, end_m);
    ranges.extend_from_slice(&seeds);
    let ranges = sorted_dedup(ranges);

    let mut audit = Vec::with_capacity(ranges.len());
    for range_m in ranges {
        let point = curve
            .at_range(range_m)
            .ok_or(CardError::DomainOutsideCurve {
                requested_m: range_m,
                curve_max_m,
            })?;
        audit.push(AuditPoint {
            range_m,
            exact_elevation: elevation.exact(point.drop_mil),
            exact_windage: windage.exact(point.wind_mil),
            printed_elevation: elevation.printed(point.drop_mil),
            printed_windage: windage.printed(point.wind_mil),
            drop_linear_m: point.drop_mil / 1000.0 * range_m,
            wind_linear_m: point.wind_mil / 1000.0 * range_m,
            velocity_mps: point.velocity_mps,
            energy_j: point.energy_j,
            time_s: point.time_s,
        });
    }

    // Seed rows: both domain ends plus every anchor. `audit` is sorted and contains each of
    // them, so a partition point gives the seed indices in ascending order.
    let mut rows: Vec<usize> = sorted_dedup(seeds)
        .iter()
        .map(|target| {
            audit
                .partition_point(|p| p.range_m < *target)
                .min(audit.len() - 1)
        })
        .collect();
    rows.dedup();
    let mut is_row = vec![false; audit.len()];
    for &r in &rows {
        is_row[r] = true;
    }
    debug_assert_eq!(rows.first(), Some(&0), "the domain start must seed row 0");
    debug_assert_eq!(
        rows.last(),
        Some(&(audit.len() - 1)),
        "the domain end must seed the last row"
    );

    // Runaway backstop, NOT the termination argument. Every iteration either stops or
    // inserts an audited point that was not already a row, so at most `audit.len()`
    // insertions are possible and the irreducible-error stop below must fire first; this
    // cap only bounds the damage if that reasoning is ever broken by a later edit.
    let iteration_cap = 2 * audit.len() + 8;
    let mut iterations = 0usize;
    let mut rows_capped = false;

    while iterations < iteration_cap {
        iterations += 1;
        let errors = sweep(&audit, &rows);

        let mut any_violation = false;
        let mut worst: Option<(f64, usize)> = None;
        for (k, &(elevation_error, windage_error)) in errors.iter().enumerate() {
            if elevation_error <= req.budget.elevation && windage_error <= req.budget.windage {
                continue;
            }
            any_violation = true;
            if is_row[k] {
                continue;
            }
            // Two axes with different budgets reduce to one comparable number by
            // budget-normalized excess. `audit` is ascending and the comparison is strict,
            // so ties keep the lowest range.
            let excess =
                (elevation_error / req.budget.elevation).max(windage_error / req.budget.windage);
            if worst.is_none_or(|(best, _)| excess > best) {
                worst = Some((excess, k));
            }
        }

        if !any_violation {
            break;
        }
        // IRREDUCIBLE-ERROR STOP. Violations remain, but every one of them is AT a row, so
        // there is nothing left to insert -- the residue is the quantization floor, not a
        // shortage of rows. Without this arm the loop re-measures the same state forever;
        // deleting it makes `quantization_floor_is_honest_and_terminates` run to the runaway
        // backstop (50 iterations for 21 audited points) and fail, which is how that test
        // earns its keep.
        let Some((_, insert_at)) = worst else {
            break;
        };
        if rows.len() >= req.max_rows {
            rows_capped = true;
            break;
        }
        rows.insert(rows.partition_point(|&r| r < insert_at), insert_at);
        is_row[insert_at] = true;
    }

    // Dense verification pass -- deliberately SEPARATE from the loop above. It re-measures
    // the finished card from scratch; nothing about `budget_met` is inherited from the
    // search, which is what stops a loop that ended early (capped or irreducible) from
    // claiming a tolerance it never reached.
    let final_errors = sweep(&audit, &rows);
    let mut worst_elevation_error = 0.0_f64;
    let mut worst_windage_error = 0.0_f64;
    let mut worst_excess = f64::NEG_INFINITY;
    let mut worst_error_range_m = audit[0].range_m;
    for (k, &(elevation_error, windage_error)) in final_errors.iter().enumerate() {
        worst_elevation_error = worst_elevation_error.max(elevation_error);
        worst_windage_error = worst_windage_error.max(windage_error);
        let excess =
            (elevation_error / req.budget.elevation).max(windage_error / req.budget.windage);
        if excess > worst_excess {
            worst_excess = excess;
            worst_error_range_m = audit[k].range_m;
        }
    }
    let budget_met =
        worst_elevation_error <= req.budget.elevation && worst_windage_error <= req.budget.windage;

    // The floor is a property of the axis, not of any one row: assert it here so the
    // invariant `assumptions[3]` states is checked on every debug build.
    debug_assert!(
        elevation.half_click_floor() >= 0.0 && windage.half_click_floor() >= 0.0,
        "a click graduation is positive by construction"
    );

    let card_rows = rows
        .iter()
        .map(|&r| {
            let p = &audit[r];
            CardRow {
                range: p.range_m,
                drop_linear: Some(p.drop_linear_m),
                drop_adj: Some(p.printed_elevation),
                come_up: None,
                wind_linear: Some(p.wind_linear_m),
                wind_adj: Some(p.printed_windage),
                velocity: Some(p.velocity_mps),
                energy: Some(p.energy_j),
                time: Some(p.time_s),
                lead_adj: None,
                wind_columns: Vec::new(),
            }
        })
        .collect();

    let report = AdaptiveCardReportV1 {
        schema_version: ADAPTIVE_CARD_SCHEMA_VERSION_V1,
        method: "greedy_worst_point_insertion_on_holdcurve_grid_v1".to_string(),
        assumptions: adaptive_card_assumptions(),
        rows: card_rows,
        budget_met,
        rows_capped,
        worst_elevation_error,
        worst_windage_error,
        worst_error_range_m,
        verification_grid_step_m: HoldCurve::SAMPLE_INTERVAL_M,
    };
    Ok((
        report,
        LoopTrace {
            iterations,
            iteration_cap,
        },
    ))
}

/// The five index-stable claims [`AdaptiveCardReportV1`] ships with. Written once, here, so
/// the report and its pinning test cannot drift apart by editing only one of them.
fn adaptive_card_assumptions() -> Vec<String> {
    [
        "Verification is limited to the hold curve's declared sample grid (verification_grid_step_m) together with the card's own rows; no claim is made about ranges between those audited points.",
        "The reader of this card interpolates linearly between adjacent rows.",
        "Errors are measured in printed-value space -- the same zero-set bias, tracking-correction division and click quantization the printed rows carry -- so a constant zero-set bias cancels out of the interpolation error and the tracking correction factor is already inside the numbers being compared.",
        "Rows quantized to an optic's clicks carry an irreducible error of up to half a click at the rows themselves, which no number of added rows can remove.",
        "A budget tighter than that half-click floor is reported as budget_met: false; the requested tolerance is never silently relaxed.",
    ]
    .iter()
    .map(|s| (*s).to_string())
    .collect()
}

#[cfg(test)]
mod adaptive_card_tests {
    use super::*;
    use crate::hold_curve::HoldCurveLoad;
    use crate::DragModel;

    /// The same representative .308-class load `hold_curve`'s own tests use, so a curve
    /// difference between the two modules cannot masquerade as a card-engine difference.
    fn test_load() -> HoldCurveLoad {
        HoldCurveLoad {
            velocity_mps: 800.0,
            bc: 0.223,
            mass_kg: 0.0109,
            diameter_m: 0.00782,
            drag_model: DragModel::G7,
            sight_height_m: 0.045,
            zero_distance_m: 100.0,
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            humidity: 50.0,
            altitude_m: 0.0,
            wind_speed_mps: 3.0,
            wind_direction_deg: 90.0,
        }
    }

    fn test_curve(max_range_m: f64) -> HoldCurve {
        HoldCurve::solve(&test_load(), max_range_m).expect("hold curve should solve")
    }

    /// An unbiased, uncorrected, unquantized request: in `Mil` the printed value is then
    /// exactly the curve's `drop_mil` / `wind_mil`, which keeps the independent checks below
    /// free of any conversion the engine could also get wrong.
    fn plain_request(domain_m: (f64, f64), budget: f64, max_rows: usize) -> AdaptiveRequest<'static> {
        AdaptiveRequest {
            domain_m,
            anchors_m: Vec::new(),
            budget: AdaptiveBudget {
                elevation: budget,
                windage: budget,
            },
            max_rows,
            click: None,
            elevation_cf: 1.0,
            windage_cf: 1.0,
            bias_mil: 0.0,
        }
    }

    /// The TEST's own audited set -- `(range_m, drop_mil, wind_mil)` at the curve's native
    /// grid inside the domain plus the two domain ends. Written out longhand rather than
    /// calling `native_grid_m`, so a bug in the engine's grid construction cannot hide
    /// behind the check that uses it.
    fn independent_audited(curve: &HoldCurve, start_m: f64, end_m: f64) -> Vec<(f64, f64, f64)> {
        let step = HoldCurve::SAMPLE_INTERVAL_M;
        let mut ranges = vec![start_m];
        let mut i = 1_i64;
        loop {
            let g = i as f64 * step;
            if g > end_m {
                break;
            }
            if g > start_m {
                ranges.push(g);
            }
            i += 1;
        }
        if ranges.last().is_none_or(|&last| last < end_m) {
            ranges.push(end_m);
        }
        ranges
            .into_iter()
            .map(|g| {
                let p = curve.at_range(g).expect("audited range must be on the curve");
                (g, p.drop_mil, p.wind_mil)
            })
            .collect()
    }

    /// The TEST's own reconstruction check: its own bracket search and hand-written lerp,
    /// deliberately NOT the engine's `sweep`. Returns the worst (elevation, windage)
    /// absolute error over `audited`.
    fn independent_worst_error(
        audited: &[(f64, f64, f64)],
        row_ranges: &[f64],
        row_elevation: &[f64],
        row_windage: &[f64],
    ) -> (f64, f64) {
        let mut worst = (0.0_f64, 0.0_f64);
        for &(g, drop_mil, wind_mil) in audited {
            // The last row at or below `g`, clamped so the final row still brackets.
            let lo = row_ranges
                .partition_point(|&r| r <= g)
                .saturating_sub(1)
                .min(row_ranges.len() - 2);
            let hi = lo + 1;
            let span = row_ranges[hi] - row_ranges[lo];
            let t = if span > 0.0 {
                (g - row_ranges[lo]) / span
            } else {
                0.0
            };
            let elevation = row_elevation[lo] + (row_elevation[hi] - row_elevation[lo]) * t;
            let windage = row_windage[lo] + (row_windage[hi] - row_windage[lo]) * t;
            worst.0 = worst.0.max((elevation - drop_mil).abs());
            worst.1 = worst.1.max((windage - wind_mil).abs());
        }
        worst
    }

    /// Smallest evenly-spaced card meeting `budget` on the same audited points, found by
    /// trial. `None` if no step up to `max_n` rows manages it.
    fn smallest_uniform_card(
        curve: &HoldCurve,
        audited: &[(f64, f64, f64)],
        start_m: f64,
        end_m: f64,
        budget: f64,
        max_n: usize,
    ) -> Option<usize> {
        for n in 2..=max_n {
            let ranges: Vec<f64> = (0..n)
                .map(|i| start_m + (end_m - start_m) * i as f64 / (n - 1) as f64)
                .collect();
            let points: Vec<_> = ranges
                .iter()
                .map(|&r| curve.at_range(r).expect("uniform row on the curve"))
                .collect();
            let elevation: Vec<f64> = points.iter().map(|p| p.drop_mil).collect();
            let windage: Vec<f64> = points.iter().map(|p| p.wind_mil).collect();
            let (worst_elevation, worst_windage) =
                independent_worst_error(audited, &ranges, &elevation, &windage);
            if worst_elevation <= budget && worst_windage <= budget {
                return Some(n);
            }
        }
        None
    }

    fn row_columns(report: &AdaptiveCardReportV1) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            report.rows.iter().map(|r| r.range).collect(),
            report
                .rows
                .iter()
                .map(|r| r.drop_adj.expect("adaptive rows always carry a dial value"))
                .collect(),
            report
                .rows
                .iter()
                .map(|r| r.wind_adj.expect("adaptive rows always carry a dial value"))
                .collect(),
        )
    }

    /// Spec 8.2 acceptance: a card that claims a met budget must survive an audit written
    /// by someone other than the engine. Every audited point is re-measured here with the
    /// test's own grid, bracket search and lerp.
    #[test]
    fn verification_pass_confirms_every_audited_point_within_bounds() {
        let curve = test_curve(900.0);
        let budget = 0.1;
        let req = plain_request((200.0, 800.0), budget, 500);
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");

        assert!(report.budget_met, "0.1 mil over 200-800 m should be reachable");
        assert!(!report.rows_capped);
        assert_eq!(report.schema_version, ADAPTIVE_CARD_SCHEMA_VERSION_V1);
        assert_eq!(report.verification_grid_step_m, HoldCurve::SAMPLE_INTERVAL_M);

        let (ranges, elevation, windage) = row_columns(&report);
        let audited = independent_audited(&curve, 200.0, 800.0);
        let (worst_elevation, worst_windage) =
            independent_worst_error(&audited, &ranges, &elevation, &windage);

        assert!(
            worst_elevation <= budget,
            "independent audit found {worst_elevation} mil of elevation error, budget {budget}"
        );
        assert!(
            worst_windage <= budget,
            "independent audit found {worst_windage} mil of windage error, budget {budget}"
        );
        // The engine's own measurement must not be optimistic relative to the independent one.
        assert!(report.worst_elevation_error >= worst_elevation - 1e-12);
        assert!(report.worst_windage_error >= worst_windage - 1e-12);
    }

    /// A tighter tolerance can never buy a shorter card. Swept across four budgets rather
    /// than compared at two points, so a non-monotone middle cannot slip through.
    #[test]
    fn tightening_the_budget_never_decreases_row_count() {
        let curve = test_curve(900.0);
        let counts: Vec<usize> = [0.4, 0.2, 0.1, 0.05]
            .iter()
            .map(|&budget| {
                let req = plain_request((200.0, 800.0), budget, 800);
                let report =
                    adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");
                assert!(report.budget_met, "{budget} mil should be reachable unquantized");
                report.rows.len()
            })
            .collect();

        for pair in counts.windows(2) {
            assert!(
                pair[1] >= pair[0],
                "row counts must not decrease as the budget tightens: {counts:?}"
            );
        }
        assert!(
            counts[3] > counts[0],
            "a 8x tighter budget should cost rows: {counts:?}"
        );
    }

    /// The search really does adapt: rows bunch up where the curve bends. Measured as the
    /// mean row spacing over the far half of the card against the near half, which is the
    /// property that survives whatever the total row count turns out to be.
    #[test]
    fn adaptive_rows_concentrate_where_the_curve_bends() {
        let curve = test_curve(900.0);
        for budget in [0.1, 0.05, 0.02] {
            let req = plain_request((200.0, 800.0), budget, 800);
            let report =
                adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");
            assert!(report.budget_met, "{budget} mil should be reachable");
            assert!(report.rows.len() >= 4, "need enough rows to halve");

            let gaps: Vec<f64> = report.rows.windows(2).map(|p| p[1].range - p[0].range).collect();
            let half = gaps.len() / 2;
            let near: f64 = gaps[..half].iter().sum::<f64>() / half as f64;
            let far: f64 = gaps[gaps.len() - half..].iter().sum::<f64>() / half as f64;
            assert!(
                far < near,
                "budget {budget}: far-half spacing {far:.1} m is not tighter than near-half \
                 {near:.1} m -- the card is not adapting, gaps {gaps:?}"
            );
        }
    }

    /// The fixed-step comparison, MEASURED rather than assumed -- and it does not come out
    /// the way the task brief predicted.
    ///
    /// The brief specified this as `smooth_trajectory_beats_fixed_step`: adaptive at 0.1 mil
    /// on 200-800 m should need fewer rows than the smallest uniform card meeting the same
    /// budget. It does not. Measured over five domains x four budgets, greedy worst-point
    /// insertion lost 10, tied 5 and won 5, and its row counts cluster on 5 / 9 / 17.
    ///
    /// That is a property of the pinned algorithm, not a defect in it. One insertion splits
    /// an interval into parts of length `l*h` and `(1-l)*h` carrying `l^2` and `(1-l)^2` of
    /// the old error, so the best any single insertion can do to an interval's error is
    /// divide it by four -- a bisection. When a trajectory's curvature barely varies (over
    /// 200-800 m this load's does so by under 2x, which is exactly what "smooth" means),
    /// every interval needs refining at once and the card doubles, while a uniform card is
    /// free to pick any row count at all. Adaptive placement wins where curvature varies
    /// sharply; on a smooth mid-range trajectory its value is the MEASURED error bound, the
    /// anchors and not having to guess a step -- not a shorter card.
    ///
    /// Two assertions, bounding the finding from BOTH sides, because a one-sided bound is
    /// not a pin. The upper bound is a regression backstop (bisection granularity must never
    /// cost more than a doubling). The directional one is the finding itself: adaptive does
    /// not meaningfully beat uniform here. Neither pins exact arithmetic a physics change
    /// would move -- the `+ 1` slack keeps a ULP-level shift at the current 5-vs-5 tie from
    /// firing spuriously -- but if adaptive ever genuinely wins, the directional assertion
    /// fails, and that is the signal that this comment, the public `adaptive_card` docs and
    /// the product guidance built on them have all gone stale.
    #[test]
    fn fixed_step_comparison_is_measured_not_assumed() {
        let curve = test_curve(900.0);
        let (start_m, end_m, budget) = (200.0, 800.0, 0.1);

        let req = plain_request((start_m, end_m), budget, 800);
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");
        assert!(report.budget_met);
        let adaptive_rows = report.rows.len();

        // Smallest uniform card meeting the same budget on the same audited points, found
        // by trial in the test -- the most generous possible fixed-step baseline.
        let audited = independent_audited(&curve, start_m, end_m);
        let uniform_rows = smallest_uniform_card(&curve, &audited, start_m, end_m, budget, 400)
            .expect("some uniform step must meet the budget");

        assert!(
            adaptive_rows <= 2 * uniform_rows,
            "adaptive used {adaptive_rows} rows against a {uniform_rows}-row uniform card; \
             bisection granularity should never cost more than a doubling"
        );
        assert!(
            adaptive_rows + 1 >= uniform_rows,
            "adaptive ({adaptive_rows} rows) now beats uniform ({uniform_rows} rows) at the \
             brief's own parameters -- the insertion rule has been improved, so the finding in \
             this test's doc comment, the \"not a shorter card\" disclosure on `adaptive_card`, \
             and the product guidance built on it are ALL stale and must be revisited"
        );
    }

    /// Anchors are promises, not suggestions -- and the same request must produce the same
    /// card every time (no RNG, no hash-order dependence anywhere in the search).
    #[test]
    fn anchors_always_present_and_determinism() {
        let curve = test_curve(900.0);
        let anchors = vec![300.0, 512.5, 777.0];
        let mut req = plain_request((200.0, 800.0), 0.15, 500);
        req.anchors_m.clone_from(&anchors);

        let first = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");
        let second = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");

        for anchor in &anchors {
            assert!(
                first.rows.iter().any(|row| row.range == *anchor),
                "anchor {anchor} m is missing from the card"
            );
        }
        // Both domain ends are rows too.
        assert_eq!(first.rows.first().map(|r| r.range), Some(200.0));
        assert_eq!(first.rows.last().map(|r| r.range), Some(800.0));

        assert_eq!(first.method, second.method);
        assert_eq!(first.assumptions, second.assumptions);
        assert_eq!(first.budget_met, second.budget_met);
        assert_eq!(first.rows_capped, second.rows_capped);
        assert_eq!(first.rows.len(), second.rows.len());
        // Bit-for-bit, not approximately: determinism means identical, not close.
        assert_eq!(
            first.worst_elevation_error.to_bits(),
            second.worst_elevation_error.to_bits()
        );
        assert_eq!(
            first.worst_windage_error.to_bits(),
            second.worst_windage_error.to_bits()
        );
        assert_eq!(
            first.worst_error_range_m.to_bits(),
            second.worst_error_range_m.to_bits()
        );
        for (a, b) in first.rows.iter().zip(second.rows.iter()) {
            assert_eq!(a.range.to_bits(), b.range.to_bits());
            assert_eq!(
                a.drop_adj.map(f64::to_bits),
                b.drop_adj.map(f64::to_bits)
            );
            assert_eq!(
                a.wind_adj.map(f64::to_bits),
                b.wind_adj.map(f64::to_bits)
            );
        }
    }

    /// A budget below the half-click floor cannot be met by adding rows, and the search
    /// must SAY so rather than grind. The iteration assertion is the fault-injection probe:
    /// with the irreducible-error stop removed the loop stops making progress and runs to
    /// the runaway backstop, blowing this bound.
    #[test]
    fn quantization_floor_is_honest_and_terminates() {
        let curve = test_curve(900.0);
        let step = HoldCurve::SAMPLE_INTERVAL_M;
        // Domain ends chosen ON the native grid so the audited set is exactly 21 points,
        // whichever way the ceil/floor rounds at the endpoints.
        let (start_m, end_m) = (330.0 * step, 350.0 * step);
        let audited_points = 21usize;

        let click = ClickValue {
            size: 0.1,
            base: ClickBase::Mil,
        };
        let half_click = 0.05;
        let budget = 0.001;
        let req = AdaptiveRequest {
            domain_m: (start_m, end_m),
            anchors_m: Vec::new(),
            budget: AdaptiveBudget {
                elevation: budget,
                windage: budget,
            },
            max_rows: 500, // deliberately not binding: the stop must come from the floor
            click: Some((&click, &click)),
            elevation_cf: 1.0,
            windage_cf: 1.0,
            bias_mil: 0.0,
        };

        let (report, trace) = adaptive_card_traced(&curve, &req, CardAdjustmentUnit::Mil)
            .expect("card should build");

        // Termination, stated as a bound rather than demonstrated by hanging: every
        // iteration inserts a distinct audited point or stops, so one pass per point plus
        // the final deciding pass is the most the search can legitimately take.
        assert!(
            trace.iterations <= audited_points + 1,
            "search took {} iterations for {audited_points} audited points (cap {}) -- \
             the irreducible-error stop is not firing",
            trace.iterations,
            trace.iteration_cap
        );
        assert!(
            trace.iterations < trace.iteration_cap,
            "the runaway backstop, not the irreducible-error stop, ended the search"
        );

        assert!(!report.budget_met, "0.001 mil is under the 0.05 mil floor");
        assert!(
            !report.rows_capped,
            "the row cap was not binding; the stop must be attributed to the floor"
        );
        assert!(report.rows.len() <= audited_points);

        // The reported worst error IS the measured half-click residue: recomputed here from
        // the printed rows alone, with no help from the engine.
        let mut row_worst = (0.0_f64, 0.0_f64);
        for row in &report.rows {
            let point = curve.at_range(row.range).expect("row on the curve");
            row_worst.0 = row_worst
                .0
                .max((row.drop_adj.expect("dial") - point.drop_mil).abs());
            row_worst.1 = row_worst
                .1
                .max((row.wind_adj.expect("dial") - point.wind_mil).abs());
        }
        assert!(
            (report.worst_elevation_error - row_worst.0).abs() < 1e-12,
            "worst elevation {} vs independently measured row residue {}",
            report.worst_elevation_error,
            row_worst.0
        );
        assert!(report.worst_elevation_error > budget);
        assert!(
            report.worst_elevation_error <= half_click + 1e-12,
            "residue {} exceeded the half-click floor",
            report.worst_elevation_error
        );
        assert!(report.worst_windage_error <= half_click + 1e-12);
    }

    /// The row cap is honoured and admitted to, never papered over with a met budget.
    #[test]
    fn max_rows_caps_with_capped_flag() {
        let curve = test_curve(900.0);
        let req = plain_request((200.0, 800.0), 0.001, 5);
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");

        assert!(report.rows_capped);
        assert!(!report.budget_met);
        assert_eq!(report.rows.len(), 5);
        assert!(report.worst_elevation_error > 0.001);
        assert!(report.worst_error_range_m >= 200.0 && report.worst_error_range_m <= 800.0);
    }

    /// The method string and all five assumptions, pinned by length AND by exact content at
    /// every index -- a consumer that quotes `assumptions[3]` must keep getting the
    /// half-click floor and not whatever a later edit shuffled into that slot.
    #[test]
    fn report_carries_method_and_all_five_assumptions() {
        let curve = test_curve(900.0);
        let req = plain_request((200.0, 400.0), 0.2, 50);
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");

        assert_eq!(
            report.method,
            "greedy_worst_point_insertion_on_holdcurve_grid_v1"
        );
        assert_eq!(report.assumptions.len(), 5);
        assert_eq!(
            report.assumptions[0],
            "Verification is limited to the hold curve's declared sample grid (verification_grid_step_m) together with the card's own rows; no claim is made about ranges between those audited points."
        );
        assert_eq!(
            report.assumptions[1],
            "The reader of this card interpolates linearly between adjacent rows."
        );
        assert_eq!(
            report.assumptions[2],
            "Errors are measured in printed-value space -- the same zero-set bias, tracking-correction division and click quantization the printed rows carry -- so a constant zero-set bias cancels out of the interpolation error and the tracking correction factor is already inside the numbers being compared."
        );
        assert_eq!(
            report.assumptions[3],
            "Rows quantized to an optic's clicks carry an irreducible error of up to half a click at the rows themselves, which no number of added rows can remove."
        );
        assert_eq!(
            report.assumptions[4],
            "A budget tighter than that half-click floor is reported as budget_met: false; the requested tolerance is never silently relaxed."
        );
    }

    /// The verification grid's one load-bearing assumption: the curve's native samples sit
    /// at exact multiples of `SAMPLE_INTERVAL_M`, so `native_grid_m` reproduces them without
    /// reaching into the curve's private sample vector.
    ///
    /// Discriminating, not merely consistent: the interpolated drop is piecewise linear with
    /// its kinks AT the sample nodes, so three probes straddling a claimed node are
    /// measurably non-collinear while three probes inside one claimed interval are collinear
    /// to floating-point noise. A grid that was offset from the real nodes would swap those
    /// two outcomes.
    #[test]
    fn verification_grid_lands_on_the_curves_native_samples() {
        let curve = test_curve(900.0);
        let step = HoldCurve::SAMPLE_INTERVAL_M;
        let max_m = curve.max_sampled_range_m();

        // The last sample is an exact multiple of the step, bit-for-bit.
        let index = (max_m / step).round();
        assert_eq!((index * step).to_bits(), max_m.to_bits());
        assert!(curve.at_range(max_m).is_some());
        assert!(curve.at_range(max_m + step).is_none());

        // Linear drop at a range, reconstructed from the angular reading.
        let drop_m_at = |range_m: f64| {
            let p = curve.at_range(range_m).expect("probe on the curve");
            p.drop_mil / 1000.0 * range_m
        };
        let node = 800.0 * step; // a claimed node, out where the curve bends hardest
        let delta = step / 4.0;

        let bend_at = |centre: f64| {
            let mid = drop_m_at(centre);
            let avg = 0.5 * (drop_m_at(centre - delta) + drop_m_at(centre + delta));
            (mid - avg).abs()
        };
        let at_node = bend_at(node);
        let inside_interval = bend_at(node + step / 2.0);

        assert!(
            inside_interval < 1e-12,
            "probes inside one claimed sample interval were not collinear ({inside_interval} m) \
             -- the reconstructed grid is offset from the curve's real nodes"
        );
        assert!(
            at_node > 1e-9 && at_node > 100.0 * inside_interval.max(f64::MIN_POSITIVE),
            "no interpolation kink at the claimed node ({at_node} m) -- \
             the reconstructed grid is offset from the curve's real nodes"
        );

        // And the reconstructed grid is inside the domain, ascending, and step-spaced.
        let grid = native_grid_m(300.0, 300.0 + 10.0 * step);
        assert!(grid.len() >= 10);
        for pair in grid.windows(2) {
            assert!((pair[1] - pair[0] - step).abs() < 1e-12);
        }
    }

    /// The locked composition order (MBA-1360 x MBA-1358 x MBA-724): bias joins the TRUE
    /// angular need first, the tracking correction divides second, quantization is last.
    /// The wrong-order value is computed too, so the test provably fails if the pipeline
    /// is ever reordered.
    #[test]
    fn printed_pipeline_keeps_the_locked_bias_then_cf_then_quantize_order() {
        let curve = test_curve(900.0);
        let click = ClickValue {
            size: 0.1,
            base: ClickBase::Mil,
        };
        // Magnitudes chosen so the two composition orders land on DIFFERENT detents: they
        // differ by `bias * (1/cf - 1)`, which must clear a half click (0.05 mil) or
        // quantization would erase the very thing this test is trying to observe.
        let (bias_mil, cf) = (2.0, 0.9);
        let req = AdaptiveRequest {
            domain_m: (300.0, 600.0),
            anchors_m: Vec::new(),
            budget: AdaptiveBudget {
                elevation: 0.2,
                windage: 0.2,
            },
            max_rows: 200,
            click: Some((&click, &click)),
            elevation_cf: cf,
            windage_cf: 1.0,
            bias_mil,
        };
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).expect("card should build");

        let row = &report.rows[1];
        let true_mil = curve.at_range(row.range).expect("row on the curve").drop_mil;

        let right_order = ((true_mil + bias_mil) / cf / 0.1).round() * 0.1;
        let wrong_order = (((true_mil / cf) + bias_mil) / 0.1).round() * 0.1;
        let dialed = row.drop_adj.expect("dial");

        assert!(
            (dialed - right_order).abs() < 1e-12,
            "printed {dialed} is not (true + bias) / cf quantized ({right_order})"
        );
        assert!(
            (right_order - wrong_order).abs() > 1e-9,
            "the chosen bias/CF make both orders agree; this test would not catch a swap"
        );

        // Windage carries its own CF and no bias, and is quantized on the same lattice.
        let true_wind = curve.at_range(row.range).expect("row on the curve").wind_mil;
        let expected_wind = (true_wind / 0.1).round() * 0.1;
        assert!((row.wind_adj.expect("dial") - expected_wind).abs() < 1e-12);
    }

    /// MOA cards are drawn on this crate's locked printed-table ratio (MBA-724), never on
    /// the exact-angle 3437.7467.
    #[test]
    fn moa_cards_use_the_locked_3438_ratio() {
        assert_eq!(CardAdjustmentUnit::Mil.from_mil_factor(), 1.0);
        assert_eq!(CardAdjustmentUnit::Moa.from_mil_factor(), 3438.0 / 1000.0);
        assert_ne!(
            CardAdjustmentUnit::Moa.from_mil_factor(),
            3437.7467 / 1000.0
        );

        let curve = test_curve(900.0);
        let req = plain_request((300.0, 600.0), 0.5, 200);
        let report = adaptive_card(&curve, &req, CardAdjustmentUnit::Moa).expect("card should build");
        let row = &report.rows[0];
        let true_mil = curve.at_range(row.range).expect("row on the curve").drop_mil;
        assert_eq!(
            row.drop_adj.expect("dial").to_bits(),
            (true_mil * (3438.0 / 1000.0)).to_bits()
        );
    }

    /// Every rejection is structured and specific -- a range a shooter asked for and did
    /// not get back must come back as a reason, never as a silently shortened card.
    #[test]
    fn request_validation_reports_every_structured_error() {
        let curve = test_curve(900.0);
        let curve_max_m = curve.max_sampled_range_m();

        let mut req = plain_request((200.0, 800.0), 0.1, 0);
        assert_eq!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).unwrap_err(),
            CardError::ZeroMaxRows
        );

        req = plain_request((200.0, 800.0), 0.1, 50);
        req.budget.elevation = 0.0;
        assert_eq!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).unwrap_err(),
            CardError::NonPositiveBudget {
                axis: "elevation",
                value: 0.0
            }
        );
        req = plain_request((200.0, 800.0), 0.1, 50);
        req.budget.windage = f64::NAN;
        assert!(matches!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil),
            Err(CardError::NonPositiveBudget { axis: "windage", .. })
        ));

        for domain in [(800.0, 200.0), (0.0, 500.0), (-10.0, 500.0), (300.0, 300.0)] {
            let req = plain_request(domain, 0.1, 50);
            assert!(
                matches!(
                    adaptive_card(&curve, &req, CardAdjustmentUnit::Mil),
                    Err(CardError::EmptyOrInvertedDomain { .. })
                ),
                "domain {domain:?} must be rejected"
            );
        }

        let req = plain_request((200.0, curve_max_m + 1.0), 0.1, 50);
        assert_eq!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).unwrap_err(),
            CardError::DomainOutsideCurve {
                requested_m: curve_max_m + 1.0,
                curve_max_m
            }
        );

        let mut req = plain_request((200.0, 800.0), 0.1, 50);
        req.anchors_m = vec![900.0];
        assert_eq!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).unwrap_err(),
            CardError::AnchorOutsideDomain {
                anchor_m: 900.0,
                start_m: 200.0,
                end_m: 800.0
            }
        );

        // Every variant renders as a sentence, so a CLI can print the reason verbatim.
        assert!(CardError::ZeroMaxRows.to_string().contains("at least one row"));
    }

    /// An out-of-band ELEVATION tracking CF is rejected with the exact variant and payload.
    ///
    /// `95.0` is the specific realistic slip this guards: a percentage typed where the ratio
    /// `0.95` belongs. Unvalidated it does not blow up -- it divides every printed value to
    /// ~0, so every measured error is ~0 and the card comes back `budget_met: true` while
    /// being entirely wrong. The assertion below therefore also pins that the request is
    /// refused rather than answered.
    #[test]
    fn out_of_band_elevation_tracking_cf_is_rejected() {
        let curve = test_curve(900.0);
        for bad in [95.0, 0.0, 0.5, 1.5, 2.0, f64::INFINITY, f64::NAN] {
            let mut req = plain_request((200.0, 800.0), 0.1, 50);
            req.elevation_cf = bad;
            let err = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil)
                .expect_err("an out-of-band elevation CF must be refused, not answered");
            match err {
                CardError::InvalidTrackingCf { axis, value } => {
                    assert_eq!(axis, "elevation");
                    // NaN never equals itself; compare bit patterns so the payload is pinned
                    // for every case including the non-finite ones.
                    assert_eq!(value.to_bits(), bad.to_bits(), "payload must echo the input");
                }
                other => panic!("expected InvalidTrackingCf for {bad}, got {other:?}"),
            }
        }
        // The band's interior is accepted, so the guard is not simply refusing everything.
        let mut req = plain_request((200.0, 800.0), 0.1, 50);
        req.elevation_cf = 0.95;
        assert!(adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).is_ok());
    }

    /// Same for the WINDAGE axis -- a per-axis check, because one shared guard covering only
    /// the elevation field would pass an elevation-only test and still ship the bug.
    #[test]
    fn out_of_band_windage_tracking_cf_is_rejected() {
        let curve = test_curve(900.0);
        for bad in [95.0, 0.0, 0.5, 1.5, 2.0, f64::INFINITY, f64::NAN] {
            let mut req = plain_request((200.0, 800.0), 0.1, 50);
            req.windage_cf = bad;
            let err = adaptive_card(&curve, &req, CardAdjustmentUnit::Mil)
                .expect_err("an out-of-band windage CF must be refused, not answered");
            match err {
                CardError::InvalidTrackingCf { axis, value } => {
                    assert_eq!(axis, "windage");
                    assert_eq!(value.to_bits(), bad.to_bits(), "payload must echo the input");
                }
                other => panic!("expected InvalidTrackingCf for {bad}, got {other:?}"),
            }
        }
        let mut req = plain_request((200.0, 800.0), 0.1, 50);
        req.windage_cf = 1.05;
        assert!(adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).is_ok());

        // Elevation is reported first when both axes are bad, so the message names one
        // concrete axis rather than a vague "a tracking factor".
        let mut req = plain_request((200.0, 800.0), 0.1, 50);
        req.elevation_cf = 95.0;
        req.windage_cf = 95.0;
        assert_eq!(
            adaptive_card(&curve, &req, CardAdjustmentUnit::Mil).unwrap_err(),
            CardError::InvalidTrackingCf {
                axis: "elevation",
                value: 95.0
            }
        );
        // And it renders as a sentence that names the ratio-vs-percentage trap.
        let text = CardError::InvalidTrackingCf {
            axis: "elevation",
            value: 95.0,
        }
        .to_string();
        assert!(text.contains("elevation") && text.contains("0.5") && text.contains("1.5"), "{text}");
    }
}
