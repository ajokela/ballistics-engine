//! MBA-1349: robust hold corridors across a bounded set of NAMED segmented-wind scenarios.
//!
//! Field users usually have several concrete plausible wind calls — a low call and a high
//! call, or two different downrange patterns — rather than a distribution they actually
//! believe. A nominal trajectory hides that ambiguity; Monte Carlo demands a probability
//! model they do not have. This module takes the middle path: solve a handful of named
//! scenarios exactly, report every one of them, and report the corridor they span.
//!
//! # No probabilities, anywhere
//!
//! Nothing here assigns a likelihood to a scenario, interpolates between scenarios, or
//! folds the finite set into a standard deviation. The corridor is the range spanned by
//! the hypotheses the user supplied, at the ranges they asked about — no more, and
//! deliberately no more. A three-scenario corridor is not a confidence interval, and this
//! module never presents it as one. Statistical dispersion remains `monte-carlo --wez`'s
//! job.
//!
//! # What a "hold" is here
//!
//! The same convention the rest of the 0.31.0 reticle work uses (see [`crate::reticle`]):
//! angular milliradians the shooter must apply, with **elevation positive = hold UP** (the
//! bullet fell below the line of sight) and **windage positive = hold RIGHT** (the wind
//! pushed the bullet right).
//!
//! # Determinism
//!
//! Scenarios are sorted by name before anything is solved, so no result — corridor,
//! minimax hold, or which scenario is named as the worst case — can depend on the order
//! they appeared in the file. Every tie-break downstream is therefore over a fixed order.
//!
//! # The zero is solved ONCE, in calm air
//!
//! A rifle has one zero; the scenarios are hypotheses about the shot, not about the
//! zeroing session. Re-zeroing per scenario would fold each hypothesis into its own datum
//! and shrink the corridor toward nothing, which is exactly the ambiguity this command
//! exists to show. So the muzzle angle is solved once against
//! [`WindConditions::default()`] and reused verbatim by every scenario.

use std::error::Error;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::cli_api::{AtmosphericConditions, BallisticInputs, TrajectorySolver, UnitSystem};
use crate::drag_model::DragModel;
use crate::wind::{validate_wind_segments, WindSegment};
use crate::WindConditions;

/// Schema version carried by [`WindScenarioSetV1`] and enforced on load.
pub const WIND_SCENARIO_SET_VERSION: u32 = 1;

/// Schema version carried by [`RobustHoldReportV1`].
pub const ROBUST_HOLD_REPORT_VERSION: u32 = 1;

/// Most scenarios one run will accept.
///
/// Each scenario is a full trajectory solve, and the point of the feature is a handful of
/// concrete calls a person can actually reason about — past this, the honest tool is a
/// statistical one. Enforced BEFORE any solving so an oversized set costs nothing.
pub const MAX_WIND_SCENARIOS: usize = 8;

/// Most ranges one run will sample. Same rationale: a bounded, checkable table, not a
/// swept surface.
pub const MAX_CORRIDOR_RANGES: usize = 64;

/// Slack (meters) allowed when deciding whether a miss is inside a target.
///
/// Boundary contact COUNTS AS A FIT — a hold that puts the worst-case impact exactly on the
/// edge of the target is, by the geometry, still on the target. This epsilon exists only so
/// that a case constructed to sit exactly on the boundary is not thrown off it by the last
/// bit of a floating-point division; it is a nanometer, far below any physical meaning.
pub const TARGET_FIT_EPSILON_M: f64 = 1.0e-9;

/// Trajectory sample spacing used for every scenario, meters (~1 yard).
///
/// The requested ranges are read by linear interpolation between neighbouring samples;
/// at this spacing that is far below the resolution any hold is read to.
const SAMPLE_INTERVAL_M: f64 = 0.9144;

/// One named wind hypothesis: a label plus the downrange segments that define it.
#[derive(Debug, Clone, PartialEq)]
pub struct NamedWindScenario {
    /// Human label, carried through the whole report so provenance is never lost.
    pub name: String,
    /// Downrange wind segments, in the same representation
    /// [`crate::wind::parse_wind_segment_str`] produces.
    pub segments: Vec<WindSegment>,
}

/// A bounded, named set of wind hypotheses (MBA-1349).
#[derive(Debug, Clone, PartialEq)]
pub struct WindScenarioSetV1 {
    /// Must equal [`WIND_SCENARIO_SET_VERSION`]; anything else is rejected rather than
    /// best-effort parsed.
    pub version: u32,
    pub scenarios: Vec<NamedWindScenario>,
    /// Which scenario is the user's own best guess, if they named one. Used only as the
    /// comparison baseline for the minimax hold — it is never weighted, and never treated
    /// as more likely.
    pub nominal: Option<String>,
}

/// The target a corridor is judged against.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "snake_case")]
pub enum TargetSpec {
    /// A rectangle centered on the point of aim.
    Rect { width_m: f64, height_m: f64 },
    /// A circle centered on the point of aim.
    Circle { diameter_m: f64 },
}

/// Which norm the minimax hold minimizes.
///
/// This is not cosmetic: the two target shapes have genuinely different optimal holds, and
/// picking one silently would give a rectangle the circle's answer (or vice versa).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorridorMetric {
    /// Per-axis (L-infinity). The two axes are independent for a rectangle, so the optimum
    /// is the per-axis MIDPOINT OF THE EXTREMES on each axis separately, and the objective
    /// it minimizes is `max over scenarios of max(|Δelevation|, |Δwindage|)`.
    ///
    /// This is also the default when no target was given: with no shape to judge against,
    /// treating the axes independently is the assumption that adds the least.
    Rectangular,
    /// Euclidean (L2). For a circle the optimum is the center of the MINIMUM ENCLOSING
    /// CIRCLE of the scenario holds, and the objective is the largest Euclidean distance
    /// from the hold to any scenario — i.e. that circle's radius.
    Circular,
}

impl TargetSpec {
    /// The metric this shape implies.
    pub fn metric(self) -> CorridorMetric {
        match self {
            TargetSpec::Rect { .. } => CorridorMetric::Rectangular,
            TargetSpec::Circle { .. } => CorridorMetric::Circular,
        }
    }
}

/// The rifle, load and atmosphere every scenario is solved with. SI throughout.
#[derive(Debug, Clone, PartialEq)]
pub struct CorridorLoad {
    pub muzzle_velocity_mps: f64,
    pub bc: f64,
    pub drag_model: DragModel,
    pub mass_kg: f64,
    pub diameter_m: f64,
    pub bullet_length_m: f64,
    pub zero_distance_m: f64,
    pub sight_height_m: f64,
    pub temperature_c: f64,
    pub pressure_hpa: f64,
    /// Relative humidity, percent (0 through 100).
    pub humidity_pct: f64,
    pub altitude_m: f64,
}

/// One robust-hold run.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustHoldRequest {
    pub scenarios: WindScenarioSetV1,
    /// Ranges to report, meters. Order is preserved in the report; duplicates are rejected.
    pub ranges_m: Vec<f64>,
    pub target: Option<TargetSpec>,
    pub load: CorridorLoad,
}

/// Why a robust-hold run was rejected. Typed, and every variant is raised BEFORE any
/// trajectory is solved except [`Self::SolveFailed`].
#[derive(Debug, Clone, PartialEq)]
pub enum WindScenarioError {
    /// The set declared a version this build does not implement.
    UnsupportedVersion { version: u32, expected: u32 },
    /// No scenarios at all.
    NoScenarios,
    /// More scenarios than [`MAX_WIND_SCENARIOS`].
    TooManyScenarios { count: usize, max: usize },
    /// A scenario carried an empty or whitespace-only name.
    EmptyScenarioName { index: usize },
    /// Two scenarios shared a name, so the report could not keep their provenance apart.
    DuplicateScenarioName { name: String },
    /// A scenario carried no wind segments.
    NoSegments { name: String },
    /// A segment token could not be parsed.
    MalformedSegment {
        scenario: String,
        index: usize,
        message: String,
    },
    /// A segment parsed but violated [`crate::wind::validate_wind_segments`].
    InvalidSegment { scenario: String, message: String },
    /// `nominal` named a scenario that is not in the set.
    UnknownNominal { name: String, available: Vec<String> },
    /// No ranges were requested.
    NoRanges,
    /// More ranges than [`MAX_CORRIDOR_RANGES`].
    TooManyRanges { count: usize, max: usize },
    /// A requested range was not finite and positive.
    InvalidRange { value: f64 },
    /// The same range was requested twice.
    DuplicateRange { value: f64 },
    /// A target dimension was not finite and positive.
    InvalidTarget { message: String },
    /// A load or atmosphere value was not usable.
    InvalidLoad { field: &'static str, value: f64 },
    /// The scenario set could not be decoded at all.
    MalformedDocument { message: String },
    /// A scenario's trajectory could not be solved, or did not reach a requested range.
    SolveFailed { scenario: String, message: String },
}

impl fmt::Display for WindScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WindScenarioError::UnsupportedVersion { version, expected } => write!(
                f,
                "unsupported wind scenario set version {version}; this build implements {expected}"
            ),
            WindScenarioError::NoScenarios => {
                write!(f, "the scenario set contains no scenarios")
            }
            WindScenarioError::TooManyScenarios { count, max } => write!(
                f,
                "{count} scenarios exceeds the limit of {max}: a hold corridor is meant to \
                 span a handful of concrete wind calls, not a swept distribution"
            ),
            WindScenarioError::EmptyScenarioName { index } => {
                write!(f, "scenario #{index} has an empty name")
            }
            WindScenarioError::DuplicateScenarioName { name } => write!(
                f,
                "two scenarios are both named '{name}'; names carry provenance through the \
                 whole report and must be unique"
            ),
            WindScenarioError::NoSegments { name } => {
                write!(f, "scenario '{name}' has no wind segments")
            }
            WindScenarioError::MalformedSegment {
                scenario,
                index,
                message,
            } => write!(f, "scenario '{scenario}' segment #{index}: {message}"),
            WindScenarioError::InvalidSegment { scenario, message } => {
                write!(f, "scenario '{scenario}': {message}")
            }
            WindScenarioError::UnknownNominal { name, available } => write!(
                f,
                "nominal scenario '{name}' is not in the set (available: {})",
                available.join(", ")
            ),
            WindScenarioError::NoRanges => write!(f, "no ranges were requested"),
            WindScenarioError::TooManyRanges { count, max } => {
                write!(f, "{count} ranges exceeds the limit of {max}")
            }
            WindScenarioError::InvalidRange { value } => write!(
                f,
                "range {value} must be finite and greater than zero"
            ),
            WindScenarioError::DuplicateRange { value } => {
                write!(f, "range {value} was requested more than once")
            }
            WindScenarioError::InvalidTarget { message } => write!(f, "invalid target: {message}"),
            WindScenarioError::InvalidLoad { field, value } => {
                write!(f, "{field} ({value}) must be finite and greater than zero")
            }
            WindScenarioError::MalformedDocument { message } => {
                write!(f, "malformed wind scenario set: {message}")
            }
            WindScenarioError::SolveFailed { scenario, message } => {
                write!(f, "scenario '{scenario}' could not be solved: {message}")
            }
        }
    }
}

impl Error for WindScenarioError {}

/// One scenario's hold at one range.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScenarioHold {
    pub name: String,
    /// Milliradians the shooter holds UP.
    pub elevation_mil: f64,
    /// Milliradians the shooter holds RIGHT.
    pub windage_mil: f64,
}

/// The corridor and chosen hold at one range.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CorridorRow {
    pub range_m: f64,
    /// Every scenario, in sorted-by-name order.
    pub scenarios: Vec<ScenarioHold>,
    pub elevation_min_mil: f64,
    pub elevation_max_mil: f64,
    pub windage_min_mil: f64,
    pub windage_max_mil: f64,
    /// The hold that minimizes the worst case under the active metric.
    pub minimax_elevation_mil: f64,
    pub minimax_windage_mil: f64,
    /// Worst-case deviation from the minimax hold to any scenario, under the active
    /// metric: the larger per-axis half-span for [`CorridorMetric::Rectangular`], the
    /// minimum-enclosing-circle radius for [`CorridorMetric::Circular`].
    pub worst_case_miss_mil: f64,
    /// Per-axis worst cases from the same hold, always reported: these are what the
    /// rectangular fit check uses, and they stay meaningful under either metric.
    pub worst_case_elevation_miss_mil: f64,
    pub worst_case_windage_miss_mil: f64,
    /// Which scenario is furthest from the minimax hold under the active metric. Ties
    /// resolve to the first in sorted-by-name order.
    pub worst_case_scenario: String,
    /// The same worst-case measure computed from the DESIGNATED NOMINAL scenario's hold,
    /// when one was named. Present so the report can show what the minimax buys.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nominal_worst_case_miss_mil: Option<f64>,
    /// Whether the whole corridor is inside the target when held at the minimax hold.
    /// `None` when no target was supplied. Boundary contact counts as a fit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fits_target: Option<bool>,
}

/// A complete robust-hold report (MBA-1349).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RobustHoldReportV1 {
    /// Always [`ROBUST_HOLD_REPORT_VERSION`].
    pub version: u32,
    /// The requested ranges, meters, in the order they were asked for.
    pub ranges_m: Vec<f64>,
    /// Scenario names in the internal sorted order — the order every row's
    /// [`CorridorRow::scenarios`] uses.
    pub scenario_names: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nominal: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target: Option<TargetSpec>,
    /// `"rectangular"` or `"circular"` — which norm the minimax holds minimize.
    pub metric: String,
    pub rows: Vec<CorridorRow>,
}

/// Decode a wind scenario set from JSON TEXT (MBA-1349).
///
/// Text, not a path: this module stays fs-free so it compiles for wasm32, exactly like
/// [`crate::drag_file`] and [`crate::truing_dsf`]. File reading belongs to the front end.
///
/// Segments are written as the SAME `SPEED:ANGLE:UNTIL[:VERTICAL]` tokens `--wind-segment`
/// takes, so there is one wind grammar to learn and one parser to trust:
///
/// ```json
/// {"version": 1, "nominal": "low",
///  "scenarios": [
///    {"name": "low",  "segments": ["5:90:400", "3:90:800"]},
///    {"name": "high", "segments": ["12:90:400", "9:270:800"]}]}
/// ```
///
/// `units` selects how SPEED and UNTIL are read (mph/yards imperial, m/s and meters
/// metric), matching `--wind-segment` exactly.
pub fn parse_wind_scenario_set(
    text: &str,
    units: UnitSystem,
) -> Result<WindScenarioSetV1, WindScenarioError> {
    #[derive(Deserialize)]
    struct WireScenario {
        name: String,
        segments: Vec<String>,
    }
    #[derive(Deserialize)]
    struct WireSet {
        version: u32,
        scenarios: Vec<WireScenario>,
        #[serde(default)]
        nominal: Option<String>,
    }

    let wire: WireSet =
        serde_json::from_str(text).map_err(|e| WindScenarioError::MalformedDocument {
            message: e.to_string(),
        })?;
    // Version FIRST: a future document may legitimately carry shapes this build does not
    // know, and must get the stable version error rather than a confusing field error.
    if wire.version != WIND_SCENARIO_SET_VERSION {
        return Err(WindScenarioError::UnsupportedVersion {
            version: wire.version,
            expected: WIND_SCENARIO_SET_VERSION,
        });
    }
    // The cap is checked before any segment is parsed, so an oversized document costs one
    // JSON decode and nothing more.
    if wire.scenarios.len() > MAX_WIND_SCENARIOS {
        return Err(WindScenarioError::TooManyScenarios {
            count: wire.scenarios.len(),
            max: MAX_WIND_SCENARIOS,
        });
    }

    let imperial = matches!(units, UnitSystem::Imperial);
    let mut scenarios = Vec::with_capacity(wire.scenarios.len());
    for (index, scenario) in wire.scenarios.into_iter().enumerate() {
        if scenario.name.trim().is_empty() {
            return Err(WindScenarioError::EmptyScenarioName { index });
        }
        let mut segments = Vec::with_capacity(scenario.segments.len());
        for (segment_index, token) in scenario.segments.iter().enumerate() {
            let segment = crate::wind::parse_wind_segment_str(token, imperial).map_err(|e| {
                WindScenarioError::MalformedSegment {
                    scenario: scenario.name.clone(),
                    index: segment_index,
                    message: e,
                }
            })?;
            segments.push(segment);
        }
        scenarios.push(NamedWindScenario {
            name: scenario.name,
            segments,
        });
    }

    Ok(WindScenarioSetV1 {
        version: WIND_SCENARIO_SET_VERSION,
        scenarios,
        nominal: wire.nominal,
    })
}

/// Parse a `rect:WxH` or `circle:D` target specification.
///
/// Dimensions are inches (imperial) or centimeters (metric), matching `mpbr --vital-zone`.
pub fn parse_target_spec(spec: &str, units: UnitSystem) -> Result<TargetSpec, WindScenarioError> {
    let invalid = |message: String| WindScenarioError::InvalidTarget { message };
    let to_m = |value: f64| match units {
        UnitSystem::Imperial => value * 0.0254,
        UnitSystem::Metric => value * 0.01,
    };
    let (shape, rest) = spec.split_once(':').ok_or_else(|| {
        invalid(format!(
            "'{spec}': expected rect:WIDTHxHEIGHT or circle:DIAMETER"
        ))
    })?;
    let positive = |token: &str, what: &str| -> Result<f64, WindScenarioError> {
        let value: f64 = token
            .trim()
            .parse()
            .map_err(|_| invalid(format!("'{spec}': {what} '{token}' is not a number")))?;
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid(format!(
                "'{spec}': {what} must be finite and greater than zero"
            )));
        }
        Ok(value)
    };
    match shape.trim().to_lowercase().as_str() {
        "rect" | "rectangle" => {
            let (w, h) = rest
                .split_once(['x', 'X'])
                .ok_or_else(|| invalid(format!("'{spec}': expected rect:WIDTHxHEIGHT")))?;
            Ok(TargetSpec::Rect {
                width_m: to_m(positive(w, "width")?),
                height_m: to_m(positive(h, "height")?),
            })
        }
        "circle" | "circ" => Ok(TargetSpec::Circle {
            diameter_m: to_m(positive(rest, "diameter")?),
        }),
        other => Err(invalid(format!(
            "'{spec}': unknown shape '{other}' (expected rect or circle)"
        ))),
    }
}

impl RobustHoldRequest {
    /// Validate everything cheap BEFORE any trajectory work begins.
    ///
    /// That ordering is the point of the ticket's cap criterion: an oversized or malformed
    /// request must cost a validation pass, never eight full solves.
    pub fn validate(&self) -> Result<(), WindScenarioError> {
        if self.scenarios.version != WIND_SCENARIO_SET_VERSION {
            return Err(WindScenarioError::UnsupportedVersion {
                version: self.scenarios.version,
                expected: WIND_SCENARIO_SET_VERSION,
            });
        }
        if self.scenarios.scenarios.is_empty() {
            return Err(WindScenarioError::NoScenarios);
        }
        if self.scenarios.scenarios.len() > MAX_WIND_SCENARIOS {
            return Err(WindScenarioError::TooManyScenarios {
                count: self.scenarios.scenarios.len(),
                max: MAX_WIND_SCENARIOS,
            });
        }
        let mut seen: Vec<&str> = Vec::with_capacity(self.scenarios.scenarios.len());
        for (index, scenario) in self.scenarios.scenarios.iter().enumerate() {
            if scenario.name.trim().is_empty() {
                return Err(WindScenarioError::EmptyScenarioName { index });
            }
            if seen.contains(&scenario.name.as_str()) {
                return Err(WindScenarioError::DuplicateScenarioName {
                    name: scenario.name.clone(),
                });
            }
            seen.push(&scenario.name);
            if scenario.segments.is_empty() {
                return Err(WindScenarioError::NoSegments {
                    name: scenario.name.clone(),
                });
            }
            validate_wind_segments(&scenario.segments).map_err(|e| {
                WindScenarioError::InvalidSegment {
                    scenario: scenario.name.clone(),
                    message: format!(
                        "segment #{} is invalid ({:?} violates {:?})",
                        e.index, e.field, e.rule
                    ),
                }
            })?;
        }
        if let Some(nominal) = self.scenarios.nominal.as_deref() {
            if !self
                .scenarios
                .scenarios
                .iter()
                .any(|scenario| scenario.name == nominal)
            {
                return Err(WindScenarioError::UnknownNominal {
                    name: nominal.to_string(),
                    available: self
                        .scenarios
                        .scenarios
                        .iter()
                        .map(|scenario| scenario.name.clone())
                        .collect(),
                });
            }
        }

        if self.ranges_m.is_empty() {
            return Err(WindScenarioError::NoRanges);
        }
        if self.ranges_m.len() > MAX_CORRIDOR_RANGES {
            return Err(WindScenarioError::TooManyRanges {
                count: self.ranges_m.len(),
                max: MAX_CORRIDOR_RANGES,
            });
        }
        let mut seen_ranges: Vec<f64> = Vec::with_capacity(self.ranges_m.len());
        for &range in &self.ranges_m {
            if !range.is_finite() || range <= 0.0 {
                return Err(WindScenarioError::InvalidRange { value: range });
            }
            if seen_ranges.contains(&range) {
                return Err(WindScenarioError::DuplicateRange { value: range });
            }
            seen_ranges.push(range);
        }

        if let Some(target) = self.target {
            let bad = |message: &str| WindScenarioError::InvalidTarget {
                message: message.to_string(),
            };
            match target {
                TargetSpec::Rect { width_m, height_m } => {
                    if !width_m.is_finite() || width_m <= 0.0 {
                        return Err(bad("width must be finite and greater than zero"));
                    }
                    if !height_m.is_finite() || height_m <= 0.0 {
                        return Err(bad("height must be finite and greater than zero"));
                    }
                }
                TargetSpec::Circle { diameter_m } => {
                    if !diameter_m.is_finite() || diameter_m <= 0.0 {
                        return Err(bad("diameter must be finite and greater than zero"));
                    }
                }
            }
        }

        let load = &self.load;
        for (field, value) in [
            ("muzzle velocity", load.muzzle_velocity_mps),
            ("ballistic coefficient", load.bc),
            ("bullet mass", load.mass_kg),
            ("bullet diameter", load.diameter_m),
            ("bullet length", load.bullet_length_m),
            ("zero distance", load.zero_distance_m),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(WindScenarioError::InvalidLoad { field, value });
            }
        }
        for (field, value) in [
            ("sight height", load.sight_height_m),
            ("temperature", load.temperature_c),
            ("pressure", load.pressure_hpa),
            ("humidity", load.humidity_pct),
            ("altitude", load.altitude_m),
        ] {
            if !value.is_finite() {
                return Err(WindScenarioError::InvalidLoad { field, value });
            }
        }
        Ok(())
    }
}

/// Solve one robust-hold run (MBA-1349).
///
/// Every scenario gets exactly ONE trajectory solve through the existing segmented-wind
/// machinery, sampled at roughly one-yard intervals and read at each requested range by
/// linear interpolation. The muzzle angle is solved once in calm air (see the module header) and
/// reused, so the corridor reflects only the wind hypotheses.
pub fn solve_robust_hold(
    request: &RobustHoldRequest,
) -> Result<RobustHoldReportV1, WindScenarioError> {
    request.validate()?;

    // Determinism: sort by name once, up front. Everything downstream — row order, the
    // corridor extremes, the worst-case scenario tie-break — is then over a fixed order,
    // so the caller's ordering cannot influence any output.
    let mut scenarios = request.scenarios.scenarios.clone();
    scenarios.sort_by(|a, b| a.name.cmp(&b.name));
    let scenario_names: Vec<String> = scenarios
        .iter()
        .map(|scenario| scenario.name.clone())
        .collect();

    let load = &request.load;
    let atmosphere = AtmosphericConditions {
        temperature: load.temperature_c,
        pressure: load.pressure_hpa,
        humidity: load.humidity_pct,
        altitude: load.altitude_m,
    };
    let base_inputs = BallisticInputs {
        bc_value: load.bc,
        bc_type: load.drag_model,
        bullet_mass: load.mass_kg,
        muzzle_velocity: load.muzzle_velocity_mps,
        bullet_diameter: load.diameter_m,
        bullet_length: load.bullet_length_m,
        sight_height: load.sight_height_m,
        use_rk4: true,
        ..Default::default()
    };
    let zero_angle = crate::cli_api::calculate_zero_angle_with_conditions(
        base_inputs.clone(),
        load.zero_distance_m,
        load.sight_height_m,
        WindConditions::default(),
        atmosphere.clone(),
    )
    .map_err(|e| WindScenarioError::SolveFailed {
        scenario: "<zero>".to_string(),
        message: e.to_string(),
    })?;

    let furthest_m = request
        .ranges_m
        .iter()
        .fold(0.0_f64, |acc, &range| acc.max(range));
    let max_range_m = furthest_m * 1.02;

    // Per scenario: one solve, then the hold at each requested range.
    let mut per_scenario: Vec<Vec<(f64, f64)>> = Vec::with_capacity(scenarios.len());
    for scenario in &scenarios {
        let mut inputs = base_inputs.clone();
        inputs.muzzle_angle = zero_angle;
        inputs.enable_trajectory_sampling = true;
        inputs.sample_interval = SAMPLE_INTERVAL_M;

        let mut solver =
            TrajectorySolver::new(inputs, WindConditions::default(), atmosphere.clone());
        solver.set_max_range(max_range_m);
        solver.set_time_step(0.001);
        solver.set_wind_segments(scenario.segments.clone());
        let result = solver.solve().map_err(|e| WindScenarioError::SolveFailed {
            scenario: scenario.name.clone(),
            message: e.to_string(),
        })?;
        let samples = result.sampled_points.unwrap_or_default();
        if samples.len() < 2 {
            return Err(WindScenarioError::SolveFailed {
                scenario: scenario.name.clone(),
                message: "the trajectory produced too few sampled points".to_string(),
            });
        }
        let mut holds = Vec::with_capacity(request.ranges_m.len());
        for &range_m in &request.ranges_m {
            let hold = interpolate_hold_mil(&samples, range_m).ok_or_else(|| {
                WindScenarioError::SolveFailed {
                    scenario: scenario.name.clone(),
                    message: format!(
                        "the trajectory does not reach {range_m:.1} m (it was sampled to {:.1} m)",
                        samples.last().map_or(0.0, |s| s.distance_m)
                    ),
                }
            })?;
            holds.push(hold);
        }
        per_scenario.push(holds);
    }

    let metric = request
        .target
        .map_or(CorridorMetric::Rectangular, TargetSpec::metric);
    let nominal_index = request
        .scenarios
        .nominal
        .as_deref()
        .and_then(|name| scenarios.iter().position(|s| s.name == name));

    let mut rows = Vec::with_capacity(request.ranges_m.len());
    for (range_index, &range_m) in request.ranges_m.iter().enumerate() {
        let points: Vec<(f64, f64)> = per_scenario
            .iter()
            .map(|holds| holds[range_index])
            .collect();

        let (elevation_min_mil, elevation_max_mil) = span(points.iter().map(|p| p.0));
        let (windage_min_mil, windage_max_mil) = span(points.iter().map(|p| p.1));

        let (minimax_elevation_mil, minimax_windage_mil) = match metric {
            CorridorMetric::Rectangular => (
                0.5 * (elevation_min_mil + elevation_max_mil),
                0.5 * (windage_min_mil + windage_max_mil),
            ),
            CorridorMetric::Circular => minimum_enclosing_circle_center(&points),
        };
        let hold = (minimax_elevation_mil, minimax_windage_mil);

        let worst_case_miss_mil = worst_case(&points, hold, metric);
        let worst_case_elevation_miss_mil = points
            .iter()
            .fold(0.0_f64, |acc, p| acc.max((p.0 - hold.0).abs()));
        let worst_case_windage_miss_mil = points
            .iter()
            .fold(0.0_f64, |acc, p| acc.max((p.1 - hold.1).abs()));
        // First in sorted-by-name order among ties.
        let worst_index = points
            .iter()
            .enumerate()
            .fold((0usize, f64::NEG_INFINITY), |(best, best_d), (i, p)| {
                let d = deviation(*p, hold, metric);
                if d > best_d {
                    (i, d)
                } else {
                    (best, best_d)
                }
            })
            .0;

        let nominal_worst_case_miss_mil =
            nominal_index.map(|index| worst_case(&points, points[index], metric));

        let fits_target = request.target.map(|target| {
            let to_linear_m = |mil: f64| mil / 1000.0 * range_m;
            match target {
                TargetSpec::Rect { width_m, height_m } => {
                    to_linear_m(worst_case_elevation_miss_mil)
                        <= height_m / 2.0 + TARGET_FIT_EPSILON_M
                        && to_linear_m(worst_case_windage_miss_mil)
                            <= width_m / 2.0 + TARGET_FIT_EPSILON_M
                }
                TargetSpec::Circle { diameter_m } => {
                    to_linear_m(worst_case_miss_mil) <= diameter_m / 2.0 + TARGET_FIT_EPSILON_M
                }
            }
        });

        rows.push(CorridorRow {
            range_m,
            scenarios: scenario_names
                .iter()
                .zip(&points)
                .map(|(name, point)| ScenarioHold {
                    name: name.clone(),
                    elevation_mil: point.0,
                    windage_mil: point.1,
                })
                .collect(),
            elevation_min_mil,
            elevation_max_mil,
            windage_min_mil,
            windage_max_mil,
            minimax_elevation_mil,
            minimax_windage_mil,
            worst_case_miss_mil,
            worst_case_elevation_miss_mil,
            worst_case_windage_miss_mil,
            worst_case_scenario: scenario_names[worst_index].clone(),
            nominal_worst_case_miss_mil,
            fits_target,
        });
    }

    Ok(RobustHoldReportV1 {
        version: ROBUST_HOLD_REPORT_VERSION,
        ranges_m: request.ranges_m.clone(),
        scenario_names,
        nominal: request.scenarios.nominal.clone(),
        target: request.target,
        metric: match metric {
            CorridorMetric::Rectangular => "rectangular".to_string(),
            CorridorMetric::Circular => "circular".to_string(),
        },
        rows,
    })
}

/// Linearly interpolate `(elevation_mil, windage_mil)` at `range_m`.
///
/// `TrajectorySample::drop_m` is positive BELOW the line of sight and `wind_drift_m`
/// positive to the RIGHT, so both map straight onto the hold convention. Milliradians use
/// the small-angle definition (1 mil subtends 1/1000 of the range) — the same conversion
/// every other 0.31.0 surface uses.
fn interpolate_hold_mil(
    samples: &[crate::trajectory_sampling::TrajectorySample],
    range_m: f64,
) -> Option<(f64, f64)> {
    if !range_m.is_finite() || range_m <= 0.0 {
        return None;
    }
    let first = samples.first()?;
    let last = samples.last()?;
    if range_m < first.distance_m || range_m > last.distance_m {
        return None;
    }
    let index = samples
        .partition_point(|s| s.distance_m < range_m)
        .min(samples.len() - 1);
    let (lo, hi) = if index == 0 {
        (&samples[0], &samples[0])
    } else {
        (&samples[index - 1], &samples[index])
    };
    let width = hi.distance_m - lo.distance_m;
    let t = if width.abs() < f64::EPSILON {
        0.0
    } else {
        (range_m - lo.distance_m) / width
    };
    let lerp = |a: f64, b: f64| a + (b - a) * t;
    Some((
        lerp(lo.drop_m, hi.drop_m) / range_m * 1000.0,
        lerp(lo.wind_drift_m, hi.wind_drift_m) / range_m * 1000.0,
    ))
}

fn span(values: impl Iterator<Item = f64>) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for value in values {
        if value < lo {
            lo = value;
        }
        if value > hi {
            hi = value;
        }
    }
    (lo, hi)
}

/// Distance from `hold` to `point` under `metric`.
fn deviation(point: (f64, f64), hold: (f64, f64), metric: CorridorMetric) -> f64 {
    let de = (point.0 - hold.0).abs();
    let dw = (point.1 - hold.1).abs();
    match metric {
        CorridorMetric::Rectangular => de.max(dw),
        CorridorMetric::Circular => (de * de + dw * dw).sqrt(),
    }
}

/// Largest deviation from `hold` to any point, under `metric`.
fn worst_case(points: &[(f64, f64)], hold: (f64, f64), metric: CorridorMetric) -> f64 {
    points
        .iter()
        .fold(0.0_f64, |acc, &p| acc.max(deviation(p, hold, metric)))
}

/// Center of the minimum enclosing circle of a small point set.
///
/// Exhaustive over the candidate circles a minimum enclosing circle can be — one point,
/// the diameter of a pair, or the circumcircle of a triple — and the smallest candidate
/// that contains every point wins. With at most [`MAX_WIND_SCENARIOS`] points that is
/// under a hundred candidates, so the exact brute force is both cheaper to reason about
/// and *deterministic*, unlike the usual randomized Welzl construction. Determinism
/// matters here: the ticket requires that reordering scenarios cannot change the chosen
/// hold.
fn minimum_enclosing_circle_center(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if n == 1 {
        return points[0];
    }
    // Slack when testing containment, in the same milliradian units as the holds: a
    // circumcircle's own defining points must count as inside despite rounding.
    const CONTAINS_EPSILON: f64 = 1.0e-9;
    let mut best: Option<((f64, f64), f64)> = None;
    let mut consider = |center: (f64, f64), radius: f64| {
        if !center.0.is_finite() || !center.1.is_finite() || !radius.is_finite() {
            return;
        }
        if points
            .iter()
            .any(|p| deviation(*p, center, CorridorMetric::Circular) > radius + CONTAINS_EPSILON)
        {
            return;
        }
        match best {
            Some((_, best_radius)) if best_radius <= radius => {}
            _ => best = Some((center, radius)),
        }
    };

    for i in 0..n {
        for j in (i + 1)..n {
            let center = (
                0.5 * (points[i].0 + points[j].0),
                0.5 * (points[i].1 + points[j].1),
            );
            let radius = deviation(points[i], center, CorridorMetric::Circular);
            consider(center, radius);
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                if let Some(center) = circumcenter(points[i], points[j], points[k]) {
                    let radius = deviation(points[i], center, CorridorMetric::Circular);
                    consider(center, radius);
                }
            }
        }
    }
    // Degenerate fallback: every point coincident, or all candidates non-finite.
    best.map_or(points[0], |(center, _)| center)
}

/// Circumcenter of three points, or `None` when they are collinear.
fn circumcenter(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> Option<(f64, f64)> {
    let d = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
    if d.abs() < 1.0e-15 {
        return None;
    }
    let a2 = a.0 * a.0 + a.1 * a.1;
    let b2 = b.0 * b.0 + b.1 * b.1;
    let c2 = c.0 * c.0 + c.1 * c.1;
    Some((
        (a2 * (b.1 - c.1) + b2 * (c.1 - a.1) + c2 * (a.1 - b.1)) / d,
        (a2 * (c.0 - b.0) + b2 * (a.0 - c.0) + c2 * (b.0 - a.0)) / d,
    ))
}

/// Rendering shape for the robust-hold report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustHoldFormat {
    Table,
    Json,
}

/// Render a robust-hold report, identically on every surface (MBA-1349).
///
/// THE formatter: the native CLI prints what this returns. The command is native-only this
/// train, but the sharing discipline is set up now so the WASM follow-up is a dispatch
/// entry rather than a second copy of the format strings — the divergence the `recoil` CSV
/// header suffered (MBA-1418) is exactly what that costs.
///
/// Returned strings are newline-terminated.
pub fn format_robust_hold_report(
    report: &RobustHoldReportV1,
    format: RobustHoldFormat,
    units: UnitSystem,
) -> String {
    let (dist_unit, size_unit) = match units {
        UnitSystem::Imperial => ("yd", "in"),
        UnitSystem::Metric => ("m", "cm"),
    };
    let range_display = |range_m: f64| match units {
        UnitSystem::Imperial => range_m / 0.9144,
        UnitSystem::Metric => range_m,
    };
    let size_display = |value_m: f64| match units {
        UnitSystem::Imperial => value_m / 0.0254,
        UnitSystem::Metric => value_m * 100.0,
    };

    match format {
        RobustHoldFormat::Json => {
            // The versioned wire shape IS the report struct: one definition, so the table
            // below and any machine consumer cannot describe different numbers.
            format!(
                "{}\n",
                serde_json::to_string_pretty(report)
                    .unwrap_or_else(|_| "{\"error\":\"serialization failed\"}".to_string())
            )
        }
        RobustHoldFormat::Table => {
            let mut out = String::new();
            out.push_str("Robust Hold Corridor\n");
            out.push_str("====================\n\n");
            out.push_str(&format!(
                "Scenarios ({}): {}\n",
                report.scenario_names.len(),
                report.scenario_names.join(", ")
            ));
            if let Some(nominal) = &report.nominal {
                out.push_str(&format!("Nominal:        {nominal}\n"));
            }
            match report.target {
                Some(TargetSpec::Rect { width_m, height_m }) => out.push_str(&format!(
                    "Target:         rectangle {:.1} x {:.1} {} (per-axis / L-inf metric)\n",
                    size_display(width_m),
                    size_display(height_m),
                    size_unit
                )),
                Some(TargetSpec::Circle { diameter_m }) => out.push_str(&format!(
                    "Target:         circle {:.1} {} across (Euclidean / L2 metric)\n",
                    size_display(diameter_m),
                    size_unit
                )),
                None => out.push_str(
                    "Target:         none — the minimax hold uses the per-axis metric\n",
                ),
            }
            out.push_str(
                "\nHolds are milliradians: elevation positive = hold UP, windage positive = \
                 hold RIGHT.\nThe corridor is the span of the scenarios you supplied. It is \
                 NOT a probability interval.\n\n",
            );

            for row in &report.rows {
                out.push_str(&format!(
                    "--- {:.0} {} ---\n",
                    range_display(row.range_m),
                    dist_unit
                ));
                out.push_str("  Scenario              Elev(mil)  Wind(mil)\n");
                for hold in &row.scenarios {
                    out.push_str(&format!(
                        "  {:<20}  {:>9.3}  {:>9.3}\n",
                        hold.name, hold.elevation_mil, hold.windage_mil
                    ));
                }
                out.push_str(&format!(
                    "  Corridor              {:>9.3}  {:>9.3}   (min)\n",
                    row.elevation_min_mil, row.windage_min_mil
                ));
                out.push_str(&format!(
                    "                        {:>9.3}  {:>9.3}   (max)\n",
                    row.elevation_max_mil, row.windage_max_mil
                ));
                out.push_str(&format!(
                    "  Minimax hold          {:>9.3}  {:>9.3}\n",
                    row.minimax_elevation_mil, row.minimax_windage_mil
                ));
                out.push_str(&format!(
                    "  Worst case from it    {:>9.3} mil ({:.2} {}), scenario '{}'\n",
                    row.worst_case_miss_mil,
                    size_display(row.worst_case_miss_mil / 1000.0 * row.range_m),
                    size_unit,
                    row.worst_case_scenario
                ));
                if let Some(nominal_worst) = row.nominal_worst_case_miss_mil {
                    out.push_str(&format!(
                        "  Holding the nominal   {:>9.3} mil ({:.2} {})\n",
                        nominal_worst,
                        size_display(nominal_worst / 1000.0 * row.range_m),
                        size_unit
                    ));
                }
                if let Some(fits) = row.fits_target {
                    out.push_str(&format!(
                        "  Fits target           {}\n",
                        if fits {
                            "yes — one hold covers every scenario"
                        } else {
                            "NO — no single hold covers every scenario here"
                        }
                    ));
                }
                out.push('\n');
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segments(tokens: &[&str]) -> Vec<WindSegment> {
        tokens
            .iter()
            .map(|token| crate::wind::parse_wind_segment_str(token, true).unwrap())
            .collect()
    }

    fn scenario(name: &str, tokens: &[&str]) -> NamedWindScenario {
        NamedWindScenario {
            name: name.to_string(),
            segments: segments(tokens),
        }
    }

    pub(super) fn test_load() -> CorridorLoad {
        CorridorLoad {
            muzzle_velocity_mps: 823.0,
            bc: 0.475,
            drag_model: DragModel::G1,
            mass_kg: 0.010_886,
            diameter_m: 0.007_823,
            bullet_length_m: 0.031,
            zero_distance_m: 91.44,
            sight_height_m: 0.0508,
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            humidity_pct: 50.0,
            altitude_m: 0.0,
        }
    }

    fn request(scenarios: Vec<NamedWindScenario>, nominal: Option<&str>) -> RobustHoldRequest {
        RobustHoldRequest {
            scenarios: WindScenarioSetV1 {
                version: 1,
                scenarios,
                nominal: nominal.map(str::to_string),
            },
            ranges_m: vec![182.88, 365.76, 548.64],
            target: None,
            load: test_load(),
        }
    }

    // ---- Acceptance criterion 2: the corridor contains every scenario ----
    #[test]
    fn corridor_contains_every_scenario_at_every_range() {
        let report = solve_robust_hold(&request(
            vec![
                scenario("low", &["4:90:1000"]),
                scenario("high", &["14:90:1000"]),
                scenario("switch", &["10:90:400", "8:270:1000"]),
            ],
            None,
        ))
        .unwrap();
        assert_eq!(report.version, ROBUST_HOLD_REPORT_VERSION);
        assert_eq!(report.rows.len(), 3);
        for row in &report.rows {
            assert_eq!(row.scenarios.len(), 3);
            for hold in &row.scenarios {
                assert!(
                    hold.elevation_mil >= row.elevation_min_mil
                        && hold.elevation_mil <= row.elevation_max_mil,
                    "elevation {} outside [{}, {}]",
                    hold.elevation_mil,
                    row.elevation_min_mil,
                    row.elevation_max_mil
                );
                assert!(
                    hold.windage_mil >= row.windage_min_mil
                        && hold.windage_mil <= row.windage_max_mil,
                    "windage {} outside [{}, {}]",
                    hold.windage_mil,
                    row.windage_min_mil,
                    row.windage_max_mil
                );
            }
            // The scenarios genuinely disagree, so this is not a vacuous containment.
            assert!(row.windage_max_mil - row.windage_min_mil > 0.1, "{row:?}");
        }
    }

    // ---- Acceptance criterion 3: one scenario collapses the corridor ----
    #[test]
    fn a_single_scenario_gives_a_zero_width_corridor() {
        let report =
            solve_robust_hold(&request(vec![scenario("only", &["9:90:1000"])], None)).unwrap();
        for row in &report.rows {
            assert_eq!(row.elevation_min_mil, row.elevation_max_mil);
            assert_eq!(row.windage_min_mil, row.windage_max_mil);
            assert_eq!(row.minimax_elevation_mil, row.scenarios[0].elevation_mil);
            assert_eq!(row.minimax_windage_mil, row.scenarios[0].windage_mil);
            assert_eq!(row.worst_case_miss_mil, 0.0);
            assert_eq!(row.worst_case_elevation_miss_mil, 0.0);
            assert_eq!(row.worst_case_windage_miss_mil, 0.0);
        }
    }

    // ---- Acceptance criterion 4: minimax never loses to the nominal hold ----
    #[test]
    fn minimax_worst_case_never_exceeds_the_nominal_holds() {
        for (target, label) in [
            (None, "no target"),
            (
                Some(TargetSpec::Rect {
                    width_m: 0.3,
                    height_m: 0.5,
                }),
                "rect",
            ),
            (Some(TargetSpec::Circle { diameter_m: 0.3 }), "circle"),
        ] {
            for nominal in ["low", "high", "switch"] {
                let mut req = request(
                    vec![
                        scenario("low", &["4:90:1000"]),
                        scenario("high", &["14:90:1000"]),
                        scenario("switch", &["10:90:400", "8:270:1000"]),
                    ],
                    Some(nominal),
                );
                req.target = target;
                let report = solve_robust_hold(&req).unwrap();
                for row in &report.rows {
                    let nominal_worst = row.nominal_worst_case_miss_mil.unwrap();
                    assert!(
                        row.worst_case_miss_mil <= nominal_worst + 1e-12,
                        "{label}/{nominal} at {} m: minimax {} > nominal {}",
                        row.range_m,
                        row.worst_case_miss_mil,
                        nominal_worst
                    );
                }
            }
        }
    }

    // ---- Acceptance criterion 7: reordering changes nothing ----
    #[test]
    fn reordering_scenarios_changes_nothing() {
        let a = solve_robust_hold(&request(
            vec![
                scenario("low", &["4:90:1000"]),
                scenario("high", &["14:90:1000"]),
                scenario("switch", &["10:90:400", "8:270:1000"]),
            ],
            Some("high"),
        ))
        .unwrap();
        let b = solve_robust_hold(&request(
            vec![
                scenario("switch", &["10:90:400", "8:270:1000"]),
                scenario("high", &["14:90:1000"]),
                scenario("low", &["4:90:1000"]),
            ],
            Some("high"),
        ))
        .unwrap();
        assert_eq!(a, b, "scenario order must not affect any output");
        // And the internal order is alphabetical, not the input order.
        assert_eq!(a.scenario_names, vec!["high", "low", "switch"]);
    }

    // ---- Acceptance criterion 5: target fit, both shapes, including the boundary ----
    #[test]
    fn target_fit_is_correct_for_both_shapes_including_boundary_contact() {
        // A GENUINELY TWO-AXIS spread: two crosswinds spread the hold in windage, and a
        // third scenario carrying a strong vertical wind spreads it in elevation. Both
        // axes must move, or the circle (L2) and rectangle (L-infinity) metrics coincide
        // and the strict L2 > L-infinity assertion at the end proves nothing — an earlier
        // pure-crosswind fixture made exactly that mistake.
        let base = || {
            request(
                vec![
                    scenario("cross-light", &["4:90:1000"]),
                    scenario("cross-heavy", &["14:90:1000"]),
                    scenario("updraft", &["4:90:1000:10"]),
                ],
                None,
            )
        };
        let mut probe = base();
        probe.ranges_m = vec![548.64];
        let report = solve_robust_hold(&probe).unwrap();
        let row = &report.rows[0];
        let half_span_wind = row.worst_case_windage_miss_mil;
        let half_span_elev = row.worst_case_elevation_miss_mil;
        // Precondition: both axes really spread. If the vertical wind stopped moving the
        // elevation hold this fails LOUDLY rather than letting the metric test degenerate.
        assert!(
            half_span_wind > 0.05 && half_span_elev > 0.05,
            "fixture must spread BOTH axes: windage {half_span_wind}, elevation {half_span_elev}"
        );
        let range_m = row.range_m;
        let wind_linear = half_span_wind / 1000.0 * range_m;
        let elev_linear = half_span_elev / 1000.0 * range_m;

        // RECTANGLE, exactly on the boundary: width = 2 * the windage half-span.
        let mut exact = base();
        exact.ranges_m = vec![548.64];
        exact.target = Some(TargetSpec::Rect {
            width_m: 2.0 * wind_linear,
            height_m: (2.0 * elev_linear).max(0.01),
        });
        assert_eq!(
            solve_robust_hold(&exact).unwrap().rows[0].fits_target,
            Some(true),
            "boundary contact counts as a fit"
        );

        // A hair narrower does not fit.
        let mut narrow = exact.clone();
        narrow.target = Some(TargetSpec::Rect {
            width_m: 2.0 * wind_linear * 0.99,
            height_m: (2.0 * elev_linear).max(0.01),
        });
        assert_eq!(
            solve_robust_hold(&narrow).unwrap().rows[0].fits_target,
            Some(false)
        );

        // A hair wider does.
        let mut wide = exact.clone();
        wide.target = Some(TargetSpec::Rect {
            width_m: 2.0 * wind_linear * 1.01,
            height_m: (2.0 * elev_linear).max(0.01) * 1.01,
        });
        assert_eq!(
            solve_robust_hold(&wide).unwrap().rows[0].fits_target,
            Some(true)
        );

        // CIRCLE: the metric changes, so re-read the radius under it.
        let mut circle_probe = base();
        circle_probe.ranges_m = vec![548.64];
        circle_probe.target = Some(TargetSpec::Circle { diameter_m: 1.0 });
        let circle_row = solve_robust_hold(&circle_probe).unwrap().rows[0].clone();
        let radius_linear = circle_row.worst_case_miss_mil / 1000.0 * range_m;

        let mut exact_circle = circle_probe.clone();
        exact_circle.target = Some(TargetSpec::Circle {
            diameter_m: 2.0 * radius_linear,
        });
        assert_eq!(
            solve_robust_hold(&exact_circle).unwrap().rows[0].fits_target,
            Some(true),
            "boundary contact counts as a fit for a circle too"
        );

        let mut small_circle = circle_probe.clone();
        small_circle.target = Some(TargetSpec::Circle {
            diameter_m: 2.0 * radius_linear * 0.99,
        });
        assert_eq!(
            solve_robust_hold(&small_circle).unwrap().rows[0].fits_target,
            Some(false)
        );

        // And the two metrics really are different: the circular objective is the
        // Euclidean radius, which for a genuine two-axis spread STRICTLY exceeds the
        // larger half-span (the L-infinity answer). With both axes confirmed nonzero
        // above, a non-strict `>=` here would also pass if the code computed L-infinity,
        // so the strict `>` is what actually distinguishes the two metrics.
        let l_inf = half_span_wind.max(half_span_elev);
        assert!(
            circle_row.worst_case_miss_mil > l_inf + 1e-6,
            "circular metric must be a true L2 radius strictly above the L-inf half-span: \
             L2 {} vs L-inf {}",
            circle_row.worst_case_miss_mil,
            l_inf
        );
    }

    // ---- Acceptance criterion 6: caps and malformed segments, before any work ----
    #[test]
    fn caps_and_malformed_input_are_structured_errors_before_any_solve() {
        let nine: Vec<NamedWindScenario> = (0..9)
            .map(|i| scenario(&format!("s{i}"), &["5:90:1000"]))
            .collect();
        assert_eq!(
            solve_robust_hold(&request(nine, None)).unwrap_err(),
            WindScenarioError::TooManyScenarios { count: 9, max: 8 }
        );

        let mut too_many_ranges = request(vec![scenario("a", &["5:90:1000"])], None);
        too_many_ranges.ranges_m = (1..=65).map(f64::from).collect();
        assert_eq!(
            solve_robust_hold(&too_many_ranges).unwrap_err(),
            WindScenarioError::TooManyRanges { count: 65, max: 64 }
        );

        let mut no_scenarios = request(vec![], None);
        no_scenarios.ranges_m = vec![100.0];
        assert_eq!(
            solve_robust_hold(&no_scenarios).unwrap_err(),
            WindScenarioError::NoScenarios
        );

        let mut bad_segment = request(vec![scenario("a", &["5:90:1000"])], None);
        bad_segment.scenarios.scenarios[0].segments[0].until_m = -1.0;
        assert!(matches!(
            solve_robust_hold(&bad_segment).unwrap_err(),
            WindScenarioError::InvalidSegment { .. }
        ));

        let mut empty_segments = request(vec![scenario("a", &["5:90:1000"])], None);
        empty_segments.scenarios.scenarios[0].segments.clear();
        assert!(matches!(
            solve_robust_hold(&empty_segments).unwrap_err(),
            WindScenarioError::NoSegments { .. }
        ));

        let mut duplicate = request(
            vec![scenario("a", &["5:90:1000"]), scenario("a", &["9:90:1000"])],
            None,
        );
        duplicate.ranges_m = vec![100.0];
        assert!(matches!(
            solve_robust_hold(&duplicate).unwrap_err(),
            WindScenarioError::DuplicateScenarioName { .. }
        ));

        let unknown_nominal = request(vec![scenario("a", &["5:90:1000"])], Some("nope"));
        assert!(matches!(
            solve_robust_hold(&unknown_nominal).unwrap_err(),
            WindScenarioError::UnknownNominal { .. }
        ));

        let mut bad_version = request(vec![scenario("a", &["5:90:1000"])], None);
        bad_version.scenarios.version = 2;
        assert_eq!(
            solve_robust_hold(&bad_version).unwrap_err(),
            WindScenarioError::UnsupportedVersion {
                version: 2,
                expected: 1
            }
        );

        let mut duplicate_range = request(vec![scenario("a", &["5:90:1000"])], None);
        duplicate_range.ranges_m = vec![100.0, 200.0, 100.0];
        assert_eq!(
            solve_robust_hold(&duplicate_range).unwrap_err(),
            WindScenarioError::DuplicateRange { value: 100.0 }
        );

        let mut bad_load = request(vec![scenario("a", &["5:90:1000"])], None);
        bad_load.load.muzzle_velocity_mps = 0.0;
        assert!(matches!(
            solve_robust_hold(&bad_load).unwrap_err(),
            WindScenarioError::InvalidLoad { .. }
        ));
    }

    #[test]
    fn parsing_enforces_the_version_and_the_scenario_cap_before_segments() {
        // A future version is rejected even though its scenario shapes are unknown.
        let err = parse_wind_scenario_set(
            r#"{"version":2,"scenarios":[{"name":"a","segments":["x"]}]}"#,
            UnitSystem::Imperial,
        )
        .unwrap_err();
        assert_eq!(
            err,
            WindScenarioError::UnsupportedVersion {
                version: 2,
                expected: 1
            },
            "the version check must precede segment parsing"
        );

        // Nine scenarios are rejected on count, not on their (valid) contents.
        let scenarios: Vec<String> = (0..9)
            .map(|i| format!(r#"{{"name":"s{i}","segments":["5:90:1000"]}}"#))
            .collect();
        let doc = format!(
            r#"{{"version":1,"scenarios":[{}]}}"#,
            scenarios.join(",")
        );
        assert_eq!(
            parse_wind_scenario_set(&doc, UnitSystem::Imperial).unwrap_err(),
            WindScenarioError::TooManyScenarios { count: 9, max: 8 }
        );

        let err = parse_wind_scenario_set(
            r#"{"version":1,"scenarios":[{"name":"a","segments":["nope"]}]}"#,
            UnitSystem::Imperial,
        )
        .unwrap_err();
        assert!(matches!(err, WindScenarioError::MalformedSegment { .. }), "{err}");

        let err = parse_wind_scenario_set(
            r#"{"version":1,"scenarios":[{"name":"  ","segments":["5:90:1000"]}]}"#,
            UnitSystem::Imperial,
        )
        .unwrap_err();
        assert_eq!(err, WindScenarioError::EmptyScenarioName { index: 0 });

        assert!(matches!(
            parse_wind_scenario_set("not json", UnitSystem::Imperial).unwrap_err(),
            WindScenarioError::MalformedDocument { .. }
        ));

        // A valid document, with the units applied the way --wind-segment applies them.
        let set = parse_wind_scenario_set(
            r#"{"version":1,"nominal":"low",
                "scenarios":[{"name":"low","segments":["10:90:400"]}]}"#,
            UnitSystem::Imperial,
        )
        .unwrap();
        assert_eq!(set.nominal.as_deref(), Some("low"));
        assert_eq!(set.scenarios.len(), 1);
        // 10 mph -> km/h, 400 yd -> m.
        assert!((set.scenarios[0].segments[0].speed_kmh - 16.09344).abs() < 1e-9);
        assert!((set.scenarios[0].segments[0].until_m - 365.76).abs() < 1e-9);
    }

    #[test]
    fn target_spec_parsing() {
        assert_eq!(
            parse_target_spec("rect:12x18", UnitSystem::Imperial).unwrap(),
            TargetSpec::Rect {
                width_m: 12.0 * 0.0254,
                height_m: 18.0 * 0.0254
            }
        );
        assert_eq!(
            parse_target_spec("circle:10", UnitSystem::Metric).unwrap(),
            TargetSpec::Circle { diameter_m: 0.1 }
        );
        for bad in ["rect", "rect:12", "circle:-1", "blob:3", "rect:0x5", "circle:abc"] {
            assert!(
                parse_target_spec(bad, UnitSystem::Imperial).is_err(),
                "'{bad}' should be rejected"
            );
        }
    }

    #[test]
    fn minimum_enclosing_circle_is_exact_and_order_independent() {
        // Three points on a unit circle: the enclosing circle is that circle.
        let points = vec![(1.0, 0.0), (-0.5, 0.866_025_403_784_438_6), (-0.5, -0.866_025_403_784_438_6)];
        let center = minimum_enclosing_circle_center(&points);
        assert!(center.0.abs() < 1e-9 && center.1.abs() < 1e-9, "{center:?}");

        // Two points: the circle on their diameter, not the bounding box's center.
        let pair = vec![(0.0, 0.0), (2.0, 4.0)];
        assert_eq!(minimum_enclosing_circle_center(&pair), (1.0, 2.0));

        // A point inside three others does not move the answer, in any order. Compared
        // with a tolerance, not bit-for-bit: the candidate enumeration order changes with
        // the input order, so the winning circle can differ in its last ULP. The
        // ORDER-INDEPENDENCE the ticket requires is guaranteed one level up, by
        // `solve_robust_hold` sorting scenarios by name before this is ever called (see
        // `reordering_scenarios_changes_nothing`, which asserts exact equality of the
        // whole report).
        let mut with_interior = points.clone();
        with_interior.push((0.1, -0.05));
        let a = minimum_enclosing_circle_center(&with_interior);
        with_interior.reverse();
        let b = minimum_enclosing_circle_center(&with_interior);
        assert!(
            (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12,
            "{a:?} vs {b:?}"
        );

        assert_eq!(minimum_enclosing_circle_center(&[(3.0, 4.0)]), (3.0, 4.0));
        assert_eq!(minimum_enclosing_circle_center(&[]), (0.0, 0.0));
    }

    #[test]
    fn formatter_renders_both_shapes_and_json_is_the_wire_schema() {
        let mut req = request(
            vec![
                scenario("low", &["4:90:1000"]),
                scenario("high", &["14:90:1000"]),
            ],
            Some("low"),
        );
        req.target = Some(TargetSpec::Rect {
            width_m: 0.3,
            height_m: 0.5,
        });
        let report = solve_robust_hold(&req).unwrap();

        // `-o json` IS the versioned wire schema: it deserializes straight back into the
        // report type, field for field. Values are compared with a tolerance rather than
        // bit-for-bit because serde_json's default float parser is accurate to within an
        // ULP, not exactly round-tripping (that needs its `float_roundtrip` feature) —
        // which is a property of the JSON reader, not of this schema.
        let json = format_robust_hold_report(&report, RobustHoldFormat::Json, UnitSystem::Imperial);
        let back: RobustHoldReportV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(back.version, report.version);
        assert_eq!(back.scenario_names, report.scenario_names);
        assert_eq!(back.nominal, report.nominal);
        assert_eq!(back.metric, report.metric);
        assert_eq!(back.rows.len(), report.rows.len());
        for (got, want) in back.rows.iter().zip(&report.rows) {
            assert_eq!(got.worst_case_scenario, want.worst_case_scenario);
            assert_eq!(got.fits_target, want.fits_target);
            assert_eq!(got.scenarios.len(), want.scenarios.len());
            for (a, b) in got.scenarios.iter().zip(&want.scenarios) {
                assert_eq!(a.name, b.name);
                assert!((a.elevation_mil - b.elevation_mil).abs() < 1e-12);
                assert!((a.windage_mil - b.windage_mil).abs() < 1e-12);
            }
            for (a, b) in [
                (got.range_m, want.range_m),
                (got.elevation_min_mil, want.elevation_min_mil),
                (got.elevation_max_mil, want.elevation_max_mil),
                (got.windage_min_mil, want.windage_min_mil),
                (got.windage_max_mil, want.windage_max_mil),
                (got.minimax_elevation_mil, want.minimax_elevation_mil),
                (got.minimax_windage_mil, want.minimax_windage_mil),
                (got.worst_case_miss_mil, want.worst_case_miss_mil),
            ] {
                assert!((a - b).abs() < 1e-12, "{a} vs {b}");
            }
        }

        let table =
            format_robust_hold_report(&report, RobustHoldFormat::Table, UnitSystem::Imperial);
        assert!(table.contains("Robust Hold Corridor"));
        assert!(table.contains("Minimax hold"));
        assert!(table.contains("Holding the nominal"));
        assert!(table.contains("Fits target"));
        assert!(
            table.contains("NOT a probability interval"),
            "the non-goal must be stated where it is read: {table}"
        );
        assert!(table.ends_with('\n'));
    }
}
