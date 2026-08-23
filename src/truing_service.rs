//! Transport-free tall-target and `dsf` derivation services.
//!
//! ## Tall-target
//!
//! Replicates the CLI's `tall-target` arithmetic exactly — same imperial/metric divisor
//! branch, same four validation guards, same tracking-CF band check via
//! [`crate::profile::validate_tracking_cf`] — as a versioned request/result pair the bridge
//! can expose to embedded consumers. The CLI's `Commands::TallTarget` arm keeps its
//! `println!` formatting but calls [`tall_target_v1`] for the arithmetic, so the two callers
//! cannot drift apart.
//!
//! ## `dsf` derivation (MBA-1357 Task 9)
//!
//! [`derive_dsf_point_v1`] is the model-driven entry point a JSON bridge caller with only a
//! bare [`crate::truing::TruingModelInputsV1`] reaches: it widens the model into
//! [`crate::truing_dsf::DsfSolveInputs`] via the existing `From` impl (no wind, level shot,
//! no offsets, ICAO reference, no BC segments — see that impl's own doc comment), solves,
//! and derives one Mach-keyed DSF point.
//!
//! The CLI's saved-profile `dsf` command has a richer input than any bare
//! `TruingModelInputsV1` can express — wind, shooting angle, zero-POI offsets, sight-mount
//! offset, BC reference, BC segments, a custom drag curve (MBA-1357 Task 8's
//! `DsfSolveInputs`, and the two regressions its review caught when that fidelity was first
//! silently narrowed). So the CLI keeps building its own `DsfSolveInputs` from the loaded
//! profile and solving directly via [`crate::truing_dsf::solve_for_dsf`] — unchanged from
//! Task 8 — rather than going through [`derive_dsf_point_v1`]'s narrower model. What it DOES
//! share with [`derive_dsf_point_v1`], so the two callers' arithmetic cannot drift apart, is
//! [`derive_dsf_point_from_solve_v1`]: the post-solve derivation (interpolation, Mach
//! computation, both gates, the DSF ratio) factored out of [`derive_dsf_point_v1`] so it can
//! run against a `TrajectoryResult` from EITHER solve path. `derive_dsf_point_v1` is a thin
//! "build inputs, solve, then call the shared derivation" wrapper around it.
//!
//! **Derive only.** Neither function ever writes a DSF table, opens a profile, touches the
//! filesystem, or calls `std::process::exit` — persistence and printing are the CLI's job.
//!
//! No feature gate: must compile for wasm32-unknown-unknown. The tall-target half depends
//! only on `crate::adjustment` and `crate::profile`; the `dsf` half depends on
//! `crate::truing`, `crate::truing_dsf`, and `crate::cli_api::TrajectoryResult`, all
//! unconditional and filesystem-free.

use serde::{Deserialize, Serialize};

use crate::adjustment::{drop_to_adjustment, AdjustmentUnit};
use crate::cli_api::TrajectoryResult;
use crate::optic::{plan_corrections, AngularCorrection, DialPlanReportV1, OpticError, OpticProfile, Preferences};
use crate::truing::{DropUnit, TruingModelInputsV1};
use crate::truing_dsf::{
    dsf_observation_warrants_90pct_warning, interpolate_position_and_velocity,
    mach_1_crossing_range_m, solve_for_dsf, DsfSolveInputs, DSF_MACH_CEILING,
};

/// Request for a tall-target correction-factor computation.
///
/// `metric` selects the measured-travel unit exactly like the CLI's `--units` flag:
/// `false` measures `measured` in inches at yards (divisor 36.0), `true` measures it in
/// centimeters at meters (divisor 100.0).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TallTargetRequestV1 {
    pub dialed: f64,
    pub measured: f64,
    pub range: f64,
    pub unit: AdjustmentUnit,
    pub metric: bool,
}

/// Result of a tall-target correction-factor computation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct TallTargetResultV1 {
    pub dialed: f64,
    pub actual: f64,
    pub correction_factor: f64,
    pub within_accepted_band: bool,
}

/// Errors returned by [`tall_target_v1`]. One variant per guard (not a single stringly-typed
/// variant) so a caller — in particular a future JSON bridge — can match exhaustively instead
/// of sniffing message text; adding a fifth guard here forces every `match` on this type to be
/// updated at compile time.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum TallTargetErrorV1 {
    #[error("clicks is not an angular unit; enter the dialed travel in mil, moa, smoa, or iphy")]
    ClicksNotAngular,
    #[error("dialed must be a positive angular travel")]
    InvalidDialed,
    #[error("measured must be a positive measured travel")]
    InvalidMeasured,
    #[error("range must be at least 1 yard/meter")]
    InvalidRange,
}

/// Computes the tall-target correction factor: same arithmetic and validation guards as the
/// CLI's `tall-target` command, minus the printing.
pub fn tall_target_v1(req: &TallTargetRequestV1) -> Result<TallTargetResultV1, TallTargetErrorV1> {
    if req.unit == AdjustmentUnit::Clicks {
        return Err(TallTargetErrorV1::ClicksNotAngular);
    }
    if !req.dialed.is_finite() || req.dialed <= 0.0 {
        return Err(TallTargetErrorV1::InvalidDialed);
    }
    if !req.measured.is_finite() || req.measured <= 0.0 {
        return Err(TallTargetErrorV1::InvalidMeasured);
    }
    if !req.range.is_finite() || req.range < 1.0 {
        return Err(TallTargetErrorV1::InvalidRange);
    }
    let drop_len = if req.metric {
        req.measured / 100.0
    } else {
        req.measured / 36.0
    };
    let actual = drop_to_adjustment(drop_len, req.range, req.unit);
    let correction_factor = actual / req.dialed;
    Ok(TallTargetResultV1 {
        dialed: req.dialed,
        actual,
        correction_factor,
        within_accepted_band: crate::profile::validate_tracking_cf(
            correction_factor,
            "the computed factor",
        )
        .is_ok(),
    })
}

// ============================================================================
// `dsf` derivation (MBA-1357 Task 9)
// ============================================================================

/// Mirrors the CLI's own `DSF_DEFAULT_SOLVE_ENVELOPE` (1000 yd) — the solve envelope
/// `trajectory --saved-profile NAME`/`come-ups --profile NAME` use by default, absent an
/// explicit `--max-range`. A bare model-driven request has no saved profile's own
/// configured max range to reuse, so this is the model-side equivalent: solve out to 1000
/// yd, or 5% past the requested observation range when that range itself exceeds 1000 yd,
/// so the solved points always bracket the requested range.
const DSF_DEFAULT_SOLVE_ENVELOPE_YD: f64 = 1000.0;

fn dsf_solve_envelope_m(range_yd: f64) -> f64 {
    let default_envelope_m = DSF_DEFAULT_SOLVE_ENVELOPE_YD * 0.9144;
    let range_m = range_yd * 0.9144;
    if range_m > default_envelope_m {
        range_m * 1.05
    } else {
        default_envelope_m
    }
}

/// Request for a `dsf` point derivation from a bare scalar-BC model (MBA-1357 Task 9).
///
/// `range_yd` is the observation's horizontal range, yards. `observed_drop` is the
/// physically observed drop at that range, in `drop_unit` — the same
/// numeric-value-plus-unit-suffix pair the CLI's `--observed-drop` flag parses (e.g.
/// `5.1mil`, `17.4moa`, `42.0in`).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DsfDeriveRequestV1 {
    pub model: TruingModelInputsV1,
    pub range_yd: f64,
    pub observed_drop: f64,
    pub drop_unit: DropUnit,
}

/// One derived Mach-keyed DSF point, plus any non-fatal warnings the derivation raised.
///
/// `predicted_drop` and `observed_drop` are both expressed in `drop_unit` (mirroring the
/// request), so a caller can compare them directly.
#[derive(Debug, Clone, Serialize)]
pub struct DsfPointResultV1 {
    pub mach: f64,
    pub dsf: f64,
    pub predicted_drop: f64,
    pub observed_drop: f64,
    pub drop_unit: DropUnit,
    pub warnings: Vec<String>,
}

/// Errors returned by [`derive_dsf_point_v1`] and [`derive_dsf_point_from_solve_v1`]. One
/// variant per failure mode (not a single stringly-typed variant) so a caller — in
/// particular the JSON bridge — can match exhaustively instead of sniffing message text.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum DsfServiceErrorV1 {
    #[error("invalid dsf request: {0}")]
    InvalidInput(String),
    #[error(
        "observation at Mach {mach:.2} is supersonic (ceiling {ceiling:.2}); \
         use true.fit to correct muzzle velocity instead"
    )]
    Supersonic { mach: f64, ceiling: f64 },
    #[error("trajectory does not reach {range_yd:.0} yd (solved to {solved_yd:.0} yd)")]
    OutOfRange { range_yd: f64, solved_yd: f64 },
    #[error("predicted drop at {range_yd:.0} yd is zero or non-finite; no DSF ratio exists")]
    DegenerateDrop { range_yd: f64 },
    #[error("dsf forward model failed: {0}")]
    ForwardModel(String),
}

impl DsfServiceErrorV1 {
    /// Structured detail payload for a JSON bridge error envelope (MBA-1357 Task 2's
    /// ruling: an INHERENT method, not a trait impl, so this stays in `truing_service` and
    /// keeps compiling with the `bridge` feature off).
    pub fn failure_details(&self) -> Option<serde_json::Value> {
        Some(match self {
            Self::InvalidInput(_) => serde_json::json!({"reason": "invalid_input"}),
            Self::Supersonic { mach, ceiling } => {
                serde_json::json!({"reason": "supersonic", "mach": mach, "ceiling": ceiling})
            }
            Self::OutOfRange {
                range_yd,
                solved_yd,
            } => serde_json::json!({
                "reason": "out_of_range",
                "range_yd": range_yd,
                "solved_yd": solved_yd
            }),
            Self::DegenerateDrop { range_yd } => {
                serde_json::json!({"reason": "degenerate_drop", "range_yd": range_yd})
            }
            Self::ForwardModel(_) => serde_json::json!({"reason": "forward_model"}),
        })
    }
}

/// Derive a `dsf` point from a bare scalar-BC model (MBA-1357 Task 9) — the JSON bridge's
/// entry point, with no saved profile behind it.
///
/// Builds [`DsfSolveInputs`] from `req.model` via the existing `From` impl (bridge
/// defaults: no wind, level shot, no zero-POI/sight-mount offsets, ICAO BC reference, no
/// BC segments, no custom drag curve — see that impl's own doc comment), sizes a solve
/// envelope the same way the CLI's own default does, solves, and hands the result to
/// [`derive_dsf_point_from_solve_v1`] for the actual derivation.
///
/// **Derive only** — never writes a DSF table, opens a profile, touches the filesystem, or
/// calls `std::process::exit`.
pub fn derive_dsf_point_v1(req: &DsfDeriveRequestV1) -> Result<DsfPointResultV1, DsfServiceErrorV1> {
    req.model
        .validate()
        .map_err(DsfServiceErrorV1::InvalidInput)?;
    if !req.range_yd.is_finite() || req.range_yd <= 0.0 {
        return Err(DsfServiceErrorV1::InvalidInput(
            "range_yd must be a positive, finite distance".to_string(),
        ));
    }
    if !req.observed_drop.is_finite() {
        return Err(DsfServiceErrorV1::InvalidInput(
            "observed_drop must be finite".to_string(),
        ));
    }

    let solve_inputs: DsfSolveInputs = (&req.model).into();
    let max_range_m = dsf_solve_envelope_m(req.range_yd);
    let result =
        solve_for_dsf(&solve_inputs, max_range_m).map_err(DsfServiceErrorV1::ForwardModel)?;

    let range_m = req.range_yd * 0.9144;
    derive_dsf_point_from_solve_v1(&result, range_m, req.observed_drop, req.drop_unit)
}

/// The post-solve half of a `dsf` point derivation: interpolation, Mach computation, both
/// gates, and the DSF ratio, run against an ALREADY-SOLVED `result` (MBA-1357 Task 9).
///
/// Shared by [`derive_dsf_point_v1`] (which solves from a bare `TruingModelInputsV1`) and
/// the CLI's saved-profile `dsf` command (which solves via the profile-rich
/// [`DsfSolveInputs`] instead, preserving every field Task 8 fixed — see this module's own
/// doc comment) so the two callers' arithmetic cannot drift apart.
///
/// `range_m` is the observation's horizontal range in meters. The two gates carried over
/// from the CLI exactly:
/// - A supersonic observation (`mach > DSF_MACH_CEILING`) is an error — that's
///   `true-velocity`'s job, not DSF's.
/// - The "beyond Mach 1 AND beyond 90% of the solved range" case (a COMPOUND condition,
///   not the 90% ratio alone) is a warning appended to the result's `warnings`, not a
///   failure.
/// - A zero or non-finite predicted drop is an error; no DSF ratio exists.
pub fn derive_dsf_point_from_solve_v1(
    result: &TrajectoryResult,
    range_m: f64,
    observed_drop: f64,
    drop_unit: DropUnit,
) -> Result<DsfPointResultV1, DsfServiceErrorV1> {
    let range_yd = range_m / 0.9144;

    let (position_y, velocity_mag) = interpolate_position_and_velocity(&result.points, range_m)
        .ok_or(DsfServiceErrorV1::OutOfRange {
            range_yd,
            solved_yd: result.max_range / 0.9144,
        })?;

    let predicted_drop_m = result.line_of_sight_height_m - position_y;
    let mach = if result.station_speed_of_sound_mps > 0.0 {
        velocity_mag / result.station_speed_of_sound_mps
    } else {
        0.0
    };

    // Gate 1: reject a supersonic observation outright — that's true-velocity's job.
    if mach > DSF_MACH_CEILING {
        return Err(DsfServiceErrorV1::Supersonic {
            mach,
            ceiling: DSF_MACH_CEILING,
        });
    }

    // Gate 2: warn when the observation is BOTH beyond the trajectory's Mach-1.0 crossing
    // AND beyond 90% of the solved max range (a compound condition, not the 90% ratio
    // alone — see dsf_observation_warrants_90pct_warning's own doc comment).
    let mut warnings = Vec::new();
    let crossing_m = mach_1_crossing_range_m(result);
    if dsf_observation_warrants_90pct_warning(range_m, crossing_m, result.max_range) {
        warnings.push(format!(
            "observation at {range_yd:.0} yd is beyond 90% of the solved range; solution \
             reliability degrades past this point"
        ));
    }

    let predicted_value = drop_unit.express_drop_m(predicted_drop_m, range_m);
    if !predicted_value.is_finite() || predicted_value == 0.0 {
        return Err(DsfServiceErrorV1::DegenerateDrop { range_yd });
    }
    let dsf = observed_drop / predicted_value;

    Ok(DsfPointResultV1 {
        mach,
        dsf,
        predicted_drop: predicted_value,
        observed_drop,
        drop_unit,
        warnings,
    })
}

// ============================================================================
// `dial-plan` (MBA-1348's `plan_corrections`, bridge wrapper)
// ============================================================================

/// Yards-to-metres, matching `dsf_solve_envelope_m`'s own constant above: `plan_corrections`
/// takes `range_m` in metres, but every other `true.*` command except `true.wind` is
/// imperial on the wire, so the request stays imperial and this is the one conversion the
/// service performs.
const YARDS_TO_METRES: f64 = 0.9144;

/// Serde default for `DialPlanRequestV1::elevation_cf`/`windage_cf`: 1.0, the neutral
/// tracking-correction-factor value. Most callers have no tall-target CF measured yet, so
/// requiring the field would be a trap — see the struct's own doc comment.
fn unity_cf() -> f64 {
    1.0
}

/// Request for [`dial_plan_v1`]: the bridge's wrapper around [`plan_corrections`]'s six
/// flat parameters (MBA-1348 Task 6's CLI wraps the same six from `--elevation`/`--windage`,
/// `--profile`/inline optic flags, `--range`, and `--prefer-hold`/`--max-hold`).
///
/// `range_yd` is in YARDS: `plan_corrections` itself takes metres, but every other `true.*`
/// command except `true.wind` is imperial on the wire, so this one stays imperial too and
/// [`dial_plan_v1`] converts at the boundary. `elevation_cf`/`windage_cf` default to `1.0`
/// (the neutral value) rather than being required — the CLI's own `--profile`-less inline
/// mode has no tracking-CF flags at all and always passes 1.0, and most callers have no
/// tall-target correction factor measured yet, so requiring the field here would be a trap.
///
/// The optic is supplied INLINE (`optic: OpticProfile`), never loaded from a named profile:
/// the CLI's `--profile` mode reads a saved profile from disk (`load_profile`), and that
/// filesystem dependency must not enter this module — see the module's own doc comment.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DialPlanRequestV1 {
    pub correction: AngularCorrection,
    pub optic: OpticProfile,
    /// Range in YARDS for the linear-miss-at-range residual. See the struct doc comment.
    pub range_yd: f64,
    #[serde(default = "unity_cf")]
    pub elevation_cf: f64,
    #[serde(default = "unity_cf")]
    pub windage_cf: f64,
    #[serde(default)]
    pub preferences: Preferences,
}

/// Turn a TRUE angular [`AngularCorrection`] into ranked dial/hold/hybrid plans for an
/// inline [`OpticProfile`] — the JSON bridge's entry point to [`plan_corrections`].
///
/// Contains NO arithmetic beyond the yards-to-metres conversion: everything else is
/// [`plan_corrections`] itself, so there is no second implementation of the CF rule or the
/// ranking logic to drift from the CLI's `dial-plan` command or from `plan_corrections`'s
/// own tests.
pub fn dial_plan_v1(req: &DialPlanRequestV1) -> Result<DialPlanReportV1, OpticError> {
    plan_corrections(
        req.correction,
        &req.optic,
        req.range_yd * YARDS_TO_METRES,
        req.elevation_cf,
        req.windage_cf,
        &req.preferences,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tall_target_matches_the_cli_arithmetic_and_rejects_clicks() {
        // 30 inches measured at 100 yd against 10 mil dialed.
        let req = TallTargetRequestV1 {
            dialed: 10.0,
            measured: 30.0,
            range: 100.0,
            unit: AdjustmentUnit::Mil,
            metric: false,
        };
        let r = tall_target_v1(&req).expect("computes");
        let expected_actual =
            crate::adjustment::drop_to_adjustment(30.0 / 36.0, 100.0, AdjustmentUnit::Mil);
        assert!((r.actual - expected_actual).abs() < 1e-12);
        assert!((r.correction_factor - expected_actual / 10.0).abs() < 1e-12);

        // Clicks is not an angular unit.
        let bad = TallTargetRequestV1 {
            unit: AdjustmentUnit::Clicks,
            ..req
        };
        assert_eq!(tall_target_v1(&bad), Err(TallTargetErrorV1::ClicksNotAngular));

        // Non-positive and sub-1 range are rejected, matching the CLI guards, each with its
        // own variant (not a shared stringly-typed error).
        assert_eq!(
            tall_target_v1(&TallTargetRequestV1 {
                dialed: 0.0,
                ..req
            }),
            Err(TallTargetErrorV1::InvalidDialed)
        );
        assert_eq!(
            tall_target_v1(&TallTargetRequestV1 {
                measured: -1.0,
                ..req
            }),
            Err(TallTargetErrorV1::InvalidMeasured)
        );
        assert_eq!(
            tall_target_v1(&TallTargetRequestV1 {
                range: 0.5,
                ..req
            }),
            Err(TallTargetErrorV1::InvalidRange)
        );
    }

    #[test]
    fn service_reproduces_the_cli_imperial_and_metric_ratios() {
        // Imperial: inches at yards, 36 in/yd. Metric: cm at metres, 100 cm/m.
        for (metric, measured, divisor) in [(false, 30.0, 36.0), (true, 76.2, 100.0)] {
            let r = tall_target_v1(&TallTargetRequestV1 {
                dialed: 10.0,
                measured,
                range: 100.0,
                unit: AdjustmentUnit::Mil,
                metric,
            })
            .expect("computes");
            let expected = crate::adjustment::drop_to_adjustment(
                measured / divisor,
                100.0,
                AdjustmentUnit::Mil,
            );
            assert!((r.actual - expected).abs() < 1e-12, "metric={metric}");
        }
    }


    // ---- dsf derivation (MBA-1357 Task 9) ----
    //
    // DsfSolveInputs (like the historical profile-driven solve it replaced) hardcodes a
    // 60-inch bore/muzzle height with ground-impact detection at height 0 (see
    // solve_for_dsf's own doc comment) — a typical supersonic rifle load hits that
    // "ground" well before it decelerates into the DSF table's Mach < 1.2 domain, so a
    // muzzle velocity picked to be transonic-quickly (~Mach 1.4 at the muzzle) is needed
    // to reach BOTH a clearly-supersonic short-range point and a transonic longer-range
    // point within the solved envelope. Confirmed empirically (not asserted blind): at
    // this load, Mach crosses 1.2 around 186 yd and the trajectory reaches "ground" around
    // 327 yd, so 50 yd (Mach ~1.37) and 200 yd (Mach ~1.18) below exercise exactly the two
    // gates without ever hitting `OutOfRange`.
    use crate::truing::DragModelArg;

    fn test_model() -> TruingModelInputsV1 {
        TruingModelInputsV1 {
            muzzle_velocity_fps: 1600.0,
            ballistic_coefficient: 0.243,
            drag_model: DragModelArg::G7,
            mass_gr: 168.0,
            diameter_in: 0.308,
            zero_distance_yd: 100.0,
            sight_height_in: 2.0,
            temperature_f: 59.0,
            pressure_inhg: 29.92,
            humidity_pct: 50.0,
            altitude_ft: 0.0,
        }
    }

    #[test]
    fn dsf_derives_a_point_and_refuses_a_supersonic_observation() {
        let model = test_model();
        // Transonic-band observation: derives a point.
        let ok = derive_dsf_point_v1(&DsfDeriveRequestV1 {
            model,
            range_yd: 200.0,
            observed_drop: 12.0,
            drop_unit: DropUnit::Mil,
        })
        .expect("derives");
        assert!(ok.mach <= crate::truing_dsf::DSF_MACH_CEILING);
        assert!(ok.dsf > 0.0);

        // Supersonic observation is true-velocity's job, not DSF's.
        let err = derive_dsf_point_v1(&DsfDeriveRequestV1 {
            model,
            range_yd: 50.0,
            observed_drop: 0.5,
            drop_unit: DropUnit::Mil,
        });
        assert!(err.is_err(), "a supersonic observation must be refused");
        assert!(matches!(err, Err(DsfServiceErrorV1::Supersonic { .. })));
    }

    #[test]
    fn dsf_rejects_invalid_range_and_observed_drop() {
        let model = test_model();
        assert!(matches!(
            derive_dsf_point_v1(&DsfDeriveRequestV1 {
                model,
                range_yd: 0.0,
                observed_drop: 1.0,
                drop_unit: DropUnit::Mil,
            }),
            Err(DsfServiceErrorV1::InvalidInput(_))
        ));
        assert!(matches!(
            derive_dsf_point_v1(&DsfDeriveRequestV1 {
                model,
                range_yd: 200.0,
                observed_drop: f64::NAN,
                drop_unit: DropUnit::Mil,
            }),
            Err(DsfServiceErrorV1::InvalidInput(_))
        ));
        // An invalid model (e.g. non-positive mass) is rejected via TruingModelInputsV1's
        // own validate(), not silently solved.
        let mut bad_model = model;
        bad_model.mass_gr = -1.0;
        assert!(matches!(
            derive_dsf_point_v1(&DsfDeriveRequestV1 {
                model: bad_model,
                range_yd: 200.0,
                observed_drop: 1.0,
                drop_unit: DropUnit::Mil,
            }),
            Err(DsfServiceErrorV1::InvalidInput(_))
        ));
    }

    #[test]
    fn dsf_derive_from_solve_matches_derive_dsf_point_v1() {
        // The CLI's profile-driven path calls derive_dsf_point_from_solve_v1 directly
        // against its own DsfSolveInputs solve; a bare-model DsfSolveInputs conversion
        // (via From<&TruingModelInputsV1>) must reproduce derive_dsf_point_v1's result
        // exactly, since that's what derive_dsf_point_v1 does internally.
        let model = test_model();
        let range_yd = 200.0;
        let observed_drop = 12.0;
        let drop_unit = DropUnit::Mil;

        let via_v1 = derive_dsf_point_v1(&DsfDeriveRequestV1 {
            model,
            range_yd,
            observed_drop,
            drop_unit,
        })
        .expect("derives");

        let solve_inputs: DsfSolveInputs = (&model).into();
        let max_range_m = dsf_solve_envelope_m(range_yd);
        let result = solve_for_dsf(&solve_inputs, max_range_m).expect("solves");
        let via_solve =
            derive_dsf_point_from_solve_v1(&result, range_yd * 0.9144, observed_drop, drop_unit)
                .expect("derives");

        assert_eq!(via_v1.mach, via_solve.mach);
        assert_eq!(via_v1.dsf, via_solve.dsf);
        assert_eq!(via_v1.predicted_drop, via_solve.predicted_drop);
        assert_eq!(via_v1.warnings, via_solve.warnings);
    }

    #[test]
    fn dsf_out_of_range_and_degenerate_drop_errors() {
        let model = test_model();
        // Far beyond anything this load reaches (max range ~327 yd).
        let err = derive_dsf_point_v1(&DsfDeriveRequestV1 {
            model,
            range_yd: 20_000.0,
            observed_drop: 12.0,
            drop_unit: DropUnit::Mil,
        });
        assert!(matches!(err, Err(DsfServiceErrorV1::OutOfRange { .. })));

        // Zero observed drop is a perfectly valid observation; predicted drop is never
        // zero at 200 yd (well past the 100 yd zero crossing) for a real trajectory, so
        // DegenerateDrop is a genuine edge case (only reachable right at the zero-crossing
        // range, where floating-point exactness makes it impractical to force
        // deterministically) rather than something this test can trigger. Confirm zero
        // observed_drop alone still derives a (zero) DSF rather than being rejected as
        // invalid input.
        let ok = derive_dsf_point_v1(&DsfDeriveRequestV1 {
            model,
            range_yd: 200.0,
            observed_drop: 0.0,
            drop_unit: DropUnit::Mil,
        })
        .expect("zero observed drop still derives");
        assert_eq!(ok.dsf, 0.0);
    }
}
