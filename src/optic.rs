//! Turret and reticle geometry: click detents, revolutions, physical travel, current dial
//! state, and reticle hold limits (MBA-1348).
//!
//! Before this module the crate's entire optic model was `crate::adjustment::ClickValue`
//! plus two tracking correction factors (elevation/windage CF), applied ad hoc wherever a
//! dialed number needed to become a true angular one. Real turrets have more structure
//! than that: click detents grouped into revolutions, an optional zero stop, finite
//! mechanical travel, and a *current* dialed offset from zero. Real reticles have a finite
//! usable hold extent too. `OpticProfile` gives all of that a home so a later dial/hold/
//! hybrid engagement planner has somewhere to read it from. This module is pure data and
//! validation — it does not itself decide how to engage a target.
//!
//! # MIL, always
//!
//! Every angular field in this module — travel, turret state, hold bounds — is in
//! milliradians. There is no unit-selection knob here; presenting other units (MOA, SMOA,
//! whole clicks) is a front-end concern, the same way `crate::adjustment` already handles
//! it for the CLI and the WASM terminal.
//!
//! # DIAL-space vs TRUE-angular — the one distinction this module exists to encode
//!
//! Two different kinds of "mil" appear in this crate and they are NOT interchangeable:
//!
//! - **Turret quantities** — `clicks_per_revolution`, `TravelLimits`, `TurretState`, and
//!   clicks generally — live in **DIAL-space**: the number engraved on the turret, the
//!   number a shooter actually dials. A scope with a tracking correction factor (CF)
//!   below 1.0 delivers LESS true angular motion per dialed unit than the engraving
//!   claims, so it takes MORE dial to reach a given true angle. This crate's established
//!   convention for that relationship — see `crate::adjustment::zero_banner_dial_values`
//!   and `crate::truing::scale_report_dial_values` — is that dial-space OUTPUTS are
//!   obtained by DIVIDING a true angular value by the CF. This module's turret types
//!   inherit that convention unchanged: they do not carry or apply a CF themselves, they
//!   simply live in the same DIAL-space that convention produces.
//! - **Reticle holds** elsewhere in the crate (`crate::reticle::ReticleHold`) are **TRUE
//!   angular** and are never CF-scaled — a mil subtension etched in glass does not care
//!   how the turret happens to be calibrated.
//!
//! A hold and a dialed correction for the identical point of impact are therefore
//! numerically different whenever CF != 1.0. Code that mixes the two spaces without an
//! explicit conversion will silently compute nonsense.
//!
//! # Travel and turret state are measured from the CURRENT ZERO, not the mechanical stop
//!
//! `TravelLimits` and `TurretState` are offsets from wherever the shooter has zeroed the
//! rifle, not from the turret's mechanical bottom. That is the fact a shooter actually
//! knows and can act on in the field — "I have 28 mil of up travel left from my zero" —
//! not a fact about the turret's total mechanical range, which also depends on where in
//! that range the zero happens to sit and is not tracked here.
//!
//! # `reticle_hold_bounds` is an explicit input, never derived
//!
//! See `HoldBounds` for why: `crate::reticle`'s `off_reticle` bounding box is a property
//! of which marks someone chose to author, grown by a fixed margin, not a property of the
//! scope's actual usable optical extent. It must never be reused as a stand-in for a real
//! spec.
//!
//! # `plan_corrections` and the CF rule
//!
//! [`plan_corrections`] turns a TRUE angular [`AngularCorrection`] into ranked, executable
//! [`DialPlan`]s: dial the whole correction in whole clicks, hold the whole correction on
//! the reticle, or split it (dial what the turret can reach, hold the TRUE angular
//! remainder). Getting the DIAL-space/TRUE-angular distinction above right in every arm is
//! the entire reason this function exists, so the rule is restated here in full and
//! followed literally everywhere below:
//!
//! - Turret travel, turret state, and click counts all live in DIAL space, as established
//!   above.
//! - A tracking correction factor (CF) maps one space to the other. To EXECUTE a TRUE
//!   angular need of `corr_true` mil on a turret whose CF is `cf`, the DIAL target is
//!   `corr_dial = corr_true / cf` — the CF DIVIDES going from a true need to a dial target,
//!   the same direction as `crate::adjustment::zero_banner_dial_values` and
//!   `crate::truing::scale_report_dial_values`. That target is quantized onto whole clicks
//!   IN DIAL SPACE, via `quantize_angle(corr_dial, &ClickValue { size:
//!   click_size_mil(click), base: ClickBase::Mil })` — a SYNTHETIC click value in mil so
//!   the reconstruction identity (`clicks as f64 * size + residual == corr_dial`) holds in
//!   mil regardless of whether the real click graduation is mil, MOA, or SMOA. Whatever
//!   whole click count results, the dial then EXECUTES `clicks as f64 *
//!   click_size_mil(click) * cf` TRUE angular mil — the CF MULTIPLIES going from dial
//!   clicks back to true angular, the OPPOSITE direction from how it was applied going in.
//! - Reticle holds are TRUE angular and are NEVER CF-scaled, anywhere: a hold component is
//!   computed and bounds-checked entirely in true mil, never quantized onto a click.
//!
//! `cf_dial_space_worked_example` (this module's test suite) pins the hand-derived numbers
//! for a CF of 0.98: a 5.0 true-mil correction on 0.1-mil clicks needs a dial target of
//! `5.0 / 0.98 = 5.10204...` mil, which quantizes to 51 clicks; those 51 clicks EXECUTE
//! `51 * 0.1 * 0.98 = 4.998` true mil, leaving a 0.002 true-mil hybrid hold. Swapping the
//! multiply/divide directions, or holding the DIAL-space remainder (`corr_dial - clicks *
//! size`) instead of the TRUE one (`corr_true - dial_mil_true`), silently produces a
//! different, wrong number here — this module exists specifically to make that mistake
//! impossible to make quietly.
//!
//! # Honesty: nothing is silently clamped
//!
//! A plan whose dial component cannot reach its target given the optic's declared travel
//! is still returned — clamped to the travel it actually has — but its `feasible` field is
//! `false` and the clamp is recorded in `limits_hit` as a [`LimitViolation`]. The same is
//! true of a hold that exceeds `reticle_hold_bounds` or [`Preferences::max_hold_mil`]. A
//! plan's `residual_mil` is always computed against what it actually executes, never
//! against the original request, so a clamped plan's residual honestly reflects the
//! resulting miss (see `infeasible_is_reported_never_silently_clamped`).
//!
//! Declared travel/hold data is trusted exactly as given; its ABSENCE is a different,
//! weaker claim than a KNOWN, EXCEEDED limit, and gets its own [`LimitKind`]
//! (`NoTravelData`/`NoHoldBoundData`) rather than being silently treated as "unlimited."
//! Whether a missing or exceeded limit turns off a strategy's OWN `feasible` flag depends
//! on what that strategy actually promises: [`Strategy::DialAll`] promises to deliver the
//! WHOLE correction by dialing alone, so it cannot affirm that promise without knowing
//! (and fitting inside) its travel. [`Strategy::HoldAll`] and [`Strategy::Hybrid`] promise
//! (respectively) the whole correction, or the dial shortfall's exact remainder, via the
//! reticle, so THEY cannot affirm their promise without knowing (and fitting inside) hold
//! bounds — but a `Hybrid` plan's dial component absorbing less than the full correction
//! because travel ran out is not a broken promise (the hold is DEFINED to absorb exactly
//! that shortfall, by construction), so a missing or exceeded TRAVEL limit is recorded on
//! a `Hybrid` plan purely for disclosure and never by itself makes it infeasible — only
//! its own hold not fitting does.
//!
//! `OpticProfile::zero_stop` is never read by `plan_corrections` at all, at travel-known or
//! travel-`None` alike: it is descriptive turret metadata (see its own doc comment), and
//! this module's ONLY source of truth for what the turret can physically reach is
//! `elevation_travel`/`windage_travel` being `Some` or `None` — a zero-stopped turret with
//! no declared `elevation_travel` is treated exactly like a non-zero-stopped one with no
//! declared travel: unverifiable, not "safely bounded at zero because it says zero_stop."
//!
//! # Correction-space is not reticle-space (the hold-bound mapping)
//!
//! `hold_mil` above is this module's own correction-space: `+` means "the point of impact
//! needs to move up/right," matching `AngularCorrection` and a dial's own convention (a
//! positive dial adjustment moves point of impact up/right). `HoldBounds`
//! (`up_mil`/`down_mil`/`left_mil`/`right_mil`) is `crate::reticle`'s RETICLE-space
//! instead, pinned at `src/reticle.rs:30-42`: `down_mil` is positive BELOW the optical
//! center, `right_mil` is positive to the shooter's RIGHT of it, and a holdover mark
//! (compensating a bullet that falls) sits at POSITIVE `down_mil` — below center. These
//! two spaces are OPPOSITELY signed on BOTH axes, not just one: a positive (UP)
//! correction's hold mark sits BELOW center (consumes `down_mil`), exactly like every BDC
//! reticle's long-range holdover marks, which sit below center, never above it; a positive
//! (RIGHT) correction's hold mark sits to the LEFT of center (consumes `left_mil`) by the
//! identical argument, mirrored. `AxisComputation` names its two fields
//! `bound_for_positive_correction` / `bound_for_negative_correction` (rather than e.g.
//! "`hold_up`") specifically so this mapping is stated at the type level and cannot
//! silently re-invert; see the
//! `hold_bound_mapping_matches_reticle_space_not_correction_space` test, which is built to
//! fail if it does.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::adjustment::{click_size_mil, quantize_angle, ClickBase, ClickValue};

/// A shooter's complete turret and reticle geometry (MBA-1348).
///
/// Every angular field anywhere in this type is in MIL. Turret-related fields
/// (`clicks_per_revolution`, `elevation_travel`, `windage_travel`, `turret_state`) are in
/// DIAL-space; `reticle_hold_bounds` is TRUE angular. See the module docs for what that
/// distinction means and why it matters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OpticProfile {
    /// Elevation turret's click graduation (its engraved size, e.g. 0.1 mil or 1/4 MOA).
    /// Serializes as `crate::adjustment::ClickValue`'s own suffixed string (`"0.1mil"`),
    /// not a structural object. `validate` requires a positive size.
    pub elevation_click: ClickValue,
    /// Windage turret's click graduation. Same serialized form and `validate` rule as
    /// `elevation_click`.
    pub windage_click: ClickValue,
    /// Click detents per full turret revolution, for turrets whose cap marks revolutions
    /// at all (many hunting turrets do not, hence `Option`). Must be at least 1 when
    /// present — `Some(0)` is rejected by `validate`, since a revolution with zero clicks
    /// in it is not a revolution.
    pub clicks_per_revolution: Option<u32>,
    /// Whether the elevation turret hard-stops at `elevation_travel.down_mil`, so the
    /// shooter cannot dial below zero at all. Descriptive metadata about the mechanism:
    /// `validate` does not require `elevation_travel` to be `Some` just because
    /// `zero_stop` is `true`, or vice versa.
    pub zero_stop: bool,
    /// Elevation travel remaining in each direction **from the current zero setting**
    /// (not from the mechanical bottom of the turret) — DIAL-space, mil.
    pub elevation_travel: Option<TravelLimits>,
    /// Windage travel from the current zero — DIAL-space, mil. `down_mil` is LEFT travel,
    /// `up_mil` is RIGHT travel (see `TravelLimits`).
    pub windage_travel: Option<TravelLimits>,
    /// The turret's current dialed offset from zero — DIAL-space, mil.
    pub turret_state: Option<TurretState>,
    /// The reticle's usable hold extent — TRUE angular, mil. See `HoldBounds`: this is an
    /// EXPLICIT input and is never derived from a `crate::reticle::ReticleDescription`.
    pub reticle_hold_bounds: Option<HoldBounds>,
}

/// Remaining mechanical turret travel in each direction from the current zero setting,
/// DIAL-space mil. Both fields are magnitudes — never negative — measured outward from
/// zero: `down_mil` is travel available going down (elevation) or left (windage);
/// `up_mil` is travel available going up (elevation) or right (windage). `validate`
/// rejects negative and non-finite values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TravelLimits {
    pub down_mil: f64,
    pub up_mil: f64,
}

/// The turret's current dialed position, as signed offsets from zero — DIAL-space, mil.
/// Positive `elevation_mil` is dialed up from zero; positive `windage_mil` is dialed right
/// from zero. Unlike `TravelLimits` these ARE signed (a turret can be dialed to either
/// side of zero), so `validate` does not reject a negative value on its own — it only
/// rejects a state whose magnitude on an axis exceeds that axis's declared
/// `TravelLimits`, when both are present.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TurretState {
    pub elevation_mil: f64,
    pub windage_mil: f64,
}

/// A reticle's usable hold extent in each direction, TRUE angular mil, relative to the
/// reticle's own center. All four fields are magnitudes — never negative.
///
/// This is an EXPLICIT input describing the scope's actual usable optical extent —
/// typically read off the manufacturer's spec sheet or a bench measurement — and must
/// NEVER be derived from `crate::reticle::ReticleDescription`. That type's `off_reticle`
/// bounding box is grown from wherever someone chose to author hash marks, padded by a
/// fixed margin fraction; it describes an authoring choice, not the glass. A reticle can
/// be drawn with marks stopping well short of where the image circle actually vignettes,
/// or drawn past it. Treating one as a proxy for the other would silently misstate how
/// far a shooter can actually hold.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HoldBounds {
    pub up_mil: f64,
    pub down_mil: f64,
    pub left_mil: f64,
    pub right_mil: f64,
}

/// Why an `OpticProfile` failed `validate`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum OpticError {
    /// A field that must be finite was NaN or infinite.
    #[error("{field} must be finite")]
    NonFinite { field: &'static str },
    /// A travel or hold-bound magnitude was negative.
    #[error("{field} must not be negative (got {value})")]
    NegativeLimit { field: &'static str, value: f64 },
    /// A click graduation's `size` was zero or negative. `ClickValue`'s fields are `pub`,
    /// so this is reachable by direct construction even though `parse_click_value`
    /// (`crate::adjustment`) already excludes it on every string-parsed or deserialized
    /// click — see that function and `ClickValue`'s `Deserialize` impl.
    #[error("{field} must be a positive click size (got {size})")]
    NonPositiveClickSize { field: &'static str, size: f64 },
    /// `clicks_per_revolution` was `Some(0)`. Zero clicks is not a revolution.
    #[error("clicks_per_revolution must be at least 1 when present, got 0")]
    ZeroClicksPerRevolution,
    /// `turret_state`'s value on `axis` fell outside that axis's declared
    /// `TravelLimits` (both measured from the same zero).
    #[error(
        "{axis} turret state at {dialed_mil} mil from zero is outside its travel of \
         -{down_mil}..={up_mil} mil"
    )]
    StateOutsideTravel {
        axis: &'static str,
        dialed_mil: f64,
        down_mil: f64,
        up_mil: f64,
    },
    /// A tracking correction factor (`elevation_cf`/`windage_cf` in [`plan_corrections`])
    /// was zero, negative, or non-finite (MBA-1348 review fix). `plan_corrections` divides
    /// a TRUE angular correction by the CF to get a DIAL target (see the module's "CF
    /// rule" doc section); a zero or negative CF turns that into infinity, NaN, or a
    /// direction-inverting garbage value that would silently propagate into a "plan" that
    /// looks like ordinary output. Rejected outright rather than merely warned about --
    /// `crate::adjustment::tracking_cf_in_range`'s tighter `(0.5, 1.5)` plausibility band
    /// is deliberately NOT re-enforced here, matching this crate's existing precedent
    /// (`crate::truing::scale_report_dial_values` also takes a bare, unchecked `dial_cf`)
    /// of treating that band as an advisory, caller/CLI-level concern, not a hard
    /// library-level bound.
    #[error("{field} must be a positive, finite tracking factor (got {value})")]
    NonPositiveTrackingFactor { field: &'static str, value: f64 },
}

impl OpticProfile {
    /// Checks internal consistency: every angular field is finite, both click sizes are
    /// positive, every travel/hold-bound magnitude is non-negative, `clicks_per_revolution`
    /// is a sane count, and — when both are present — `turret_state` lies within
    /// `elevation_travel` / `windage_travel`.
    ///
    /// This cannot and does not check that the profile matches any real, physical scope;
    /// it only rejects shapes that are self-contradictory or would corrupt downstream
    /// arithmetic (a zero-length revolution, NaN propagation, a dialed state the turret
    /// could not physically reach given its own declared travel, or a zero/negative click
    /// size that would turn a later `quantize_angle` division into infinity or NaN).
    pub fn validate(&self) -> Result<(), OpticError> {
        require_finite("elevation_click.size", self.elevation_click.size)?;
        require_positive_click_size("elevation_click.size", self.elevation_click.size)?;
        require_finite("windage_click.size", self.windage_click.size)?;
        require_positive_click_size("windage_click.size", self.windage_click.size)?;

        if self.clicks_per_revolution == Some(0) {
            return Err(OpticError::ZeroClicksPerRevolution);
        }

        if let Some(travel) = &self.elevation_travel {
            validate_travel("elevation_travel.down_mil", "elevation_travel.up_mil", travel)?;
        }
        if let Some(travel) = &self.windage_travel {
            validate_travel("windage_travel.down_mil", "windage_travel.up_mil", travel)?;
        }
        if let Some(state) = &self.turret_state {
            require_finite("turret_state.elevation_mil", state.elevation_mil)?;
            require_finite("turret_state.windage_mil", state.windage_mil)?;
        }
        if let Some(bounds) = &self.reticle_hold_bounds {
            for (field, value) in [
                ("reticle_hold_bounds.up_mil", bounds.up_mil),
                ("reticle_hold_bounds.down_mil", bounds.down_mil),
                ("reticle_hold_bounds.left_mil", bounds.left_mil),
                ("reticle_hold_bounds.right_mil", bounds.right_mil),
            ] {
                require_finite(field, value)?;
                require_non_negative(field, value)?;
            }
        }

        // Cross-field: the dialed state must be reachable given the declared travel, on
        // each axis independently. Runs last, after every field involved has already been
        // proven finite and non-negative above, so no NaN/negative can slip into the
        // comparison.
        if let (Some(state), Some(travel)) = (&self.turret_state, &self.elevation_travel) {
            check_state_within_travel("elevation", state.elevation_mil, travel)?;
        }
        if let (Some(state), Some(travel)) = (&self.turret_state, &self.windage_travel) {
            check_state_within_travel("windage", state.windage_mil, travel)?;
        }

        Ok(())
    }
}

/// `value` must be finite, or `field` (its dotted name, e.g. `"elevation_travel.up_mil"`)
/// is reported non-finite.
fn require_finite(field: &'static str, value: f64) -> Result<(), OpticError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(OpticError::NonFinite { field })
    }
}

/// `value` must be a non-negative magnitude, or `field` is reported as a negative limit.
/// Callers must check finiteness first: `NaN >= 0.0` is `false`, so this alone would
/// misreport a non-finite value as merely negative.
fn require_non_negative(field: &'static str, value: f64) -> Result<(), OpticError> {
    if value >= 0.0 {
        Ok(())
    } else {
        Err(OpticError::NegativeLimit { field, value })
    }
}

/// `size` must be strictly positive (a zero or negative click graduation is meaningless
/// and would make `quantize_angle`'s `angle / click.size` divide by zero or flip sign).
/// Callers must check finiteness first, for the same reason as `require_non_negative`:
/// `NaN > 0.0` is `false`, so this alone would misreport a non-finite size as merely
/// non-positive.
fn require_positive_click_size(field: &'static str, size: f64) -> Result<(), OpticError> {
    if size > 0.0 {
        Ok(())
    } else {
        Err(OpticError::NonPositiveClickSize { field, size })
    }
}

/// `value` must be a strictly positive, finite tracking correction factor, or `field` is
/// reported via [`OpticError::NonPositiveTrackingFactor`] (MBA-1348 review fix). A single
/// check for both finiteness and positivity, since a CF is consumed by DIVISION
/// (`plan_corrections`' "the CF rule"): zero divides to infinity, negative flips every
/// direction, and either would otherwise reach `quantize_angle` as silent garbage.
fn require_positive_tracking_factor(field: &'static str, value: f64) -> Result<(), OpticError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(OpticError::NonPositiveTrackingFactor { field, value })
    }
}

/// Finiteness and non-negativity for one `TravelLimits`, under the caller's own dotted
/// field names for `down_mil` / `up_mil` (so the same helper serves both the elevation
/// and windage axes without hand-rolling the check twice).
fn validate_travel(
    down_field: &'static str,
    up_field: &'static str,
    travel: &TravelLimits,
) -> Result<(), OpticError> {
    require_finite(down_field, travel.down_mil)?;
    require_finite(up_field, travel.up_mil)?;
    require_non_negative(down_field, travel.down_mil)?;
    require_non_negative(up_field, travel.up_mil)?;
    Ok(())
}

/// `dialed_mil` (already known finite) must lie within `[-travel.down_mil,
/// travel.up_mil]` (already known finite and non-negative), the closed interval the
/// turret can physically reach from zero on this axis.
fn check_state_within_travel(
    axis: &'static str,
    dialed_mil: f64,
    travel: &TravelLimits,
) -> Result<(), OpticError> {
    if dialed_mil > travel.up_mil || dialed_mil < -travel.down_mil {
        Err(OpticError::StateOutsideTravel {
            axis,
            dialed_mil,
            down_mil: travel.down_mil,
            up_mil: travel.up_mil,
        })
    } else {
        Ok(())
    }
}

/// Splits a DIAL-space click count measured from zero into `(revolutions,
/// clicks_within_the_revolution)`, for a turret whose cap marks revolutions every
/// `clicks_per_revolution` clicks — e.g. rendering "2 revolutions + 7 clicks" instead of a
/// bare click count the shooter has to do the division on in their head under time
/// pressure.
///
/// `clicks_from_zero` is DIAL-space, measured from the CURRENT ZERO (see the module
/// docs) — like every other turret quantity here, NOT from the turret's mechanical
/// bottom.
///
/// Returns `None` for:
/// - a negative `clicks_from_zero`. Turret revolution indicators count UP from zero, so a
///   setting below zero has no revolution to report: "-7 clicks" is unambiguous, while
///   any revolutions-plus-clicks decomposition of a negative count invents a direction of
///   revolution the turret does not have. Callers with a below-zero setting should render
///   the plain (negative) click count directly instead of calling this function.
/// - `clicks_per_revolution == 0`, which has no meaningful revolution length (and would
///   be a division by zero). `OpticProfile::validate` rejects `Some(0)` for the analogous
///   reason, but this free function takes a bare `u32`, not an `OpticProfile`, and so has
///   no other way to reject an invalid revolution length.
///
/// Otherwise returns `Some((revolutions, clicks_in_revolution))` with
/// `clicks_in_revolution < clicks_per_revolution` and the reconstruction identity
/// `revolutions * clicks_per_revolution + clicks_in_revolution == clicks_from_zero`
/// (ordinary non-negative integer division and remainder).
pub fn revolution_annotation(
    clicks_from_zero: i64,
    clicks_per_revolution: u32,
) -> Option<(u32, u32)> {
    if clicks_from_zero < 0 || clicks_per_revolution == 0 {
        return None;
    }
    let cpr = i64::from(clicks_per_revolution);
    let revolutions = clicks_from_zero / cpr;
    let clicks_in_revolution = clicks_from_zero % cpr;
    // clicks_from_zero >= 0 and cpr > 0 (checked above), so both quotient and remainder
    // are non-negative; a value that does not fit u32 is not a real turret setting.
    Some((revolutions as u32, clicks_in_revolution as u32))
}

// ============================================================================================
// MBA-1348 Task 4: `plan_corrections` -- dial / hold / hybrid execution planning
// ============================================================================================

/// Schema version for [`DialPlanReportV1`]'s payload shape (MBA-1348).
pub const DIAL_PLAN_SCHEMA_VERSION_V1: u32 = 1;

/// A TRUE angular correction to deliver on target -- the output of a solve, not yet
/// realized as any turret or reticle instruction. Positive `elevation_mil` is UP; positive
/// `windage_mil` is RIGHT. Never CF-scaled -- see the module's "CF rule" doc section.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AngularCorrection {
    pub elevation_mil: f64,
    pub windage_mil: f64,
}

/// Caller preferences steering [`plan_corrections`]'s ranking -- never its underlying
/// arithmetic, which is identical regardless of these values. See "Ranking is
/// deterministic" in [`plan_corrections`]'s own doc comment.
///
/// `#[derive(Default)]` (`prefer_hold: false, max_hold_mil: None` -- prefer dial, no extra
/// hold cap) rather than a hand-written impl: `bool`'s and `Option`'s own defaults already
/// produce exactly this (MBA-1348 review fix -- confirmed equivalent, not just
/// clippy-driven).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct Preferences {
    /// When ranking ties on `residual_linear_at_range_m`, prefer holding over dialing.
    /// `false` (the [`Default`] impl's value) prefers dialing.
    pub prefer_hold: bool,
    /// An optional blanket cap on any single hold component's magnitude, TRUE angular mil,
    /// applied together with (the tighter of it and) `OpticProfile::reticle_hold_bounds`.
    /// `None` (the default) applies no cap beyond `reticle_hold_bounds` itself.
    pub max_hold_mil: Option<f64>,
}

/// Which of a shot's two independent angular axes an [`AxisInstruction`] or
/// [`LimitViolation`] belongs to (MBA-1348). Unrelated to `crate::reticle_import`'s
/// private, importer-internal `Axis` of the same name in a different module.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    Elevation,
    Windage,
}

/// How a [`DialPlan`] proposes to deliver a correction (MBA-1348).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Strategy {
    /// Dial the whole correction in whole clicks; hold nothing.
    DialAll,
    /// Dial nothing; hold the whole correction on the reticle.
    HoldAll,
    /// Dial the whole clicks the turret can reach (travel-clamped if necessary); hold the
    /// TRUE angular remainder on the reticle. Always reconstructs the correction exactly --
    /// see the module's "CF rule" doc section.
    Hybrid,
}

/// Why a [`DialPlan`] component was clamped, or could not be verified feasible (MBA-1348).
/// Paired with the offending [`Axis`] in a [`LimitViolation`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LimitKind {
    /// The dial target exceeded the axis's declared `TravelLimits` and was clamped to it.
    TravelExceeded,
    /// The hold component exceeded the axis's declared `HoldBounds` and/or
    /// `Preferences::max_hold_mil`. Holds are never clamped -- see `Strategy::HoldAll`.
    HoldBoundExceeded,
    /// A nonzero dial target was needed, but the optic profile declares no `TravelLimits`
    /// for this axis, so reachability could not be verified. NOT the same claim as
    /// `TravelExceeded` -- see the module's "Honesty" doc section.
    NoTravelData,
    /// A nonzero hold was needed, but neither `HoldBounds` nor `Preferences::max_hold_mil`
    /// is declared for this axis, so fitting the reticle could not be verified.
    NoHoldBoundData,
}

/// One recorded reason a [`DialPlan`] was clamped or could not be verified feasible
/// (MBA-1348). Always DISCLOSED in `DialPlan::limits_hit`, never silently absorbed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LimitViolation {
    pub axis: Axis,
    pub kind: LimitKind,
    /// The signed value that triggered this record, BEFORE any clamp: DIAL-space mil for
    /// `TravelExceeded`/`NoTravelData` (comparable to `available_mil`, also dial-space);
    /// TRUE angular mil for `HoldBoundExceeded`/`NoHoldBoundData`.
    pub needed_mil: f64,
    /// The declared limit's magnitude on the exceeded side, when known. `None` for
    /// `NoTravelData`/`NoHoldBoundData` -- that absence is the entire point of those kinds.
    pub available_mil: Option<f64>,
}

/// `AxisInstruction::direction`: which way the DELTA actually dialed turns the turret
/// (MBA-1348). Wire form is lowercase (`"up"`/`"down"`/`"left"`/`"right"`) via
/// `#[serde(rename_all)]` -- unchanged from this type's original plain-`&str` form (MBA-1348
/// review fix M1), but now a real enum so `AxisInstruction`/`DialPlan`/`DialPlanReportV1`
/// can derive `Deserialize` again (a `&'static str` field cannot: serde's blanket impl for
/// `&'a str` ties `'a` to the deserializer's own input lifetime, which is essentially never
/// `'static`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Direction {
    Up,
    Down,
    Left,
    Right,
}

/// One axis's executable instruction within a [`DialPlan`] (MBA-1348).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AxisInstruction {
    pub axis: Axis,
    /// For the DELTA actually dialed (`delta_clicks`'s sign) -- not for the correction's
    /// own sign.
    pub direction: Direction,
    /// Clicks to turn from `OpticProfile::turret_state` (or zero, if absent) to
    /// `target_clicks_from_zero`. Positive is up/right, matching `direction`.
    pub delta_clicks: i64,
    /// The DIAL setting this instruction ends at, measured from zero -- NOT from
    /// `turret_state`. See the `turret_state_shifts_delta_but_not_target` test.
    pub target_clicks_from_zero: i64,
    /// `revolution_annotation(target_clicks_from_zero, clicks_per_revolution)`: `None`
    /// when `clicks_per_revolution` is undeclared, or `target_clicks_from_zero` is
    /// negative (see that function's own contract).
    pub end_revolution: Option<(u32, u32)>,
    /// What the dial EXECUTES in TRUE angular mil: `target_clicks_from_zero as f64 *
    /// click_size_mil(click) * cf`. Zero for `Strategy::HoldAll`.
    pub dial_mil_true: f64,
    /// The TRUE angular reticle hold component, mil (+ = over/right). Zero for
    /// `Strategy::DialAll`.
    pub hold_mil: f64,
    /// `correction − dial_mil_true − hold_mil`: the honest, post-execution miss on this
    /// axis. Exactly `0.0` for `Strategy::Hybrid` (an identity by construction -- `hold_mil`
    /// is DEFINED as whatever makes this exact -- not an approximation that merely happens
    /// to come out small).
    pub residual_mil: f64,
}

/// One ranked, executable plan for delivering a two-axis [`AngularCorrection`] (MBA-1348).
/// `instructions` is always exactly `[elevation, windage]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DialPlan {
    pub strategy: Strategy,
    pub instructions: [AxisInstruction; 2],
    /// RSS of both axes' `residual_mil`, mil-to-linear at `DialPlanReportV1::range_m` by
    /// the small-angle approximation (assumption 0). The ranking key AMONG EQUALLY-FEASIBLE
    /// plans -- `feasible` itself is checked first (MBA-1348 review fix I4) -- see
    /// [`plan_corrections`]'s doc comment.
    pub residual_linear_at_range_m: f64,
    /// `false` whenever a [`LimitViolation`] this strategy's OWN promise depends on is
    /// present in `limits_hit` -- see the module's "Honesty" doc section for what each
    /// strategy promises, and therefore what gates it.
    pub feasible: bool,
    pub limits_hit: Vec<LimitViolation>,
}

/// [`plan_corrections`]'s versioned report (MBA-1348). Carries `method` and `assumptions`
/// in the payload itself -- this train's cross-cutting honesty principle (spec §2) -- not
/// only in prose documentation.
///
/// Derives `Deserialize` as well as `Serialize` (MBA-1348 review fix M1): the original
/// `AxisInstruction::direction: &'static str` field made this impossible (see
/// [`Direction`]'s own doc), which is why that field is now a proper enum. Nothing in this
/// train currently reads a `DialPlanReportV1` back from JSON -- Task 6's CLI only ever
/// writes one out -- but there is no longer a structural reason it couldn't.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DialPlanReportV1 {
    pub schema_version: u32,
    pub method: String,
    pub assumptions: Vec<String>,
    pub range_m: f64,
    /// Ranked best-first -- see "Ranking is deterministic" in [`plan_corrections`]'s doc
    /// comment.
    pub plans: Vec<DialPlan>,
}

/// The largest whole-click count reachable within `boundary_mil` of DIAL-space travel
/// without exceeding it (MBA-1348).
///
/// Ordinarily `(boundary_mil / click_mil).floor()`, EXCEPT when that ratio is within one
/// billionth (`1e-9`) of a click of a whole number, where it rounds to that whole number
/// instead -- so floating-point noise in the division (e.g. `0.4 / 0.1` landing a hair
/// under `4.0`) cannot silently discard a click of real, declared travel (MBA-1348 review
/// fix M2: corrected from an earlier, wrong "click-thousandth" description of this same
/// `1e-9` tolerance). A genuinely fractional boundary (e.g. 0.45 mil of travel on a 0.1 mil
/// click) still floors: 4, never rounding UP past a hard mechanical limit to 5.
fn max_clicks_within(boundary_mil: f64, click_mil: f64) -> i64 {
    let raw = boundary_mil / click_mil;
    let nearest = raw.round();
    if (nearest * click_mil - boundary_mil).abs() <= click_mil * 1e-9 {
        nearest as i64
    } else {
        raw.floor() as i64
    }
}

/// Quantizes `corr_true / cf` onto whole DIAL-space clicks (the CF rule's "way in": see the
/// module doc), clamping to `travel` when declared and the unclamped target would exceed
/// it, and recording exactly one [`LimitViolation`] when either the clamp happened or a
/// nonzero target could not be checked because `travel` is `None` (MBA-1348). Returns
/// `(clamped_target_clicks, corr_dial, violation)` -- `corr_dial` (the pre-quantization
/// continuous dial-space value) is reused by callers as `LimitViolation::needed_mil`.
fn quantize_and_clamp(
    axis: Axis,
    corr_true: f64,
    click_mil: f64,
    cf: f64,
    travel: Option<&TravelLimits>,
) -> (i64, f64, Option<LimitViolation>) {
    let corr_dial = corr_true / cf;
    let synthetic = ClickValue { size: click_mil, base: ClickBase::Mil };
    let target = quantize_angle(corr_dial, &synthetic).clicks;
    let violation = match travel {
        None if target == 0 => None,
        None => Some(LimitViolation {
            axis,
            kind: LimitKind::NoTravelData,
            needed_mil: corr_dial,
            available_mil: None,
        }),
        Some(t) => {
            let max_up = max_clicks_within(t.up_mil, click_mil);
            let max_down = max_clicks_within(t.down_mil, click_mil);
            if target > max_up {
                return (
                    max_up,
                    corr_dial,
                    Some(LimitViolation {
                        axis,
                        kind: LimitKind::TravelExceeded,
                        needed_mil: corr_dial,
                        available_mil: Some(t.up_mil),
                    }),
                );
            } else if target < -max_down {
                return (
                    -max_down,
                    corr_dial,
                    Some(LimitViolation {
                        axis,
                        kind: LimitKind::TravelExceeded,
                        needed_mil: corr_dial,
                        available_mil: Some(t.down_mil),
                    }),
                );
            }
            None
        }
    };
    (target, corr_dial, violation)
}

/// Bounds-checks a TRUE angular hold value against this axis's directional `HoldBounds`
/// magnitude and `Preferences::max_hold_mil` (the tighter of the two, when both are
/// declared), recording exactly one [`LimitViolation`] when the hold exceeds it, or (for a
/// nonzero hold) when NEITHER is declared at all (MBA-1348). Holds are never clamped -- see
/// `Strategy::HoldAll` -- so this never alters `hold_true_mil`, only reports on it.
fn check_hold_bounds(
    axis: Axis,
    hold_true_mil: f64,
    bound_positive: Option<f64>,
    bound_negative: Option<f64>,
    max_hold_mil: Option<f64>,
) -> Option<LimitViolation> {
    if hold_true_mil == 0.0 {
        return None;
    }
    let directional_bound = if hold_true_mil > 0.0 { bound_positive } else { bound_negative };
    let effective_bound = match (directional_bound, max_hold_mil) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    match effective_bound {
        None => Some(LimitViolation {
            axis,
            kind: LimitKind::NoHoldBoundData,
            needed_mil: hold_true_mil,
            available_mil: None,
        }),
        Some(bound) if hold_true_mil.abs() > bound => Some(LimitViolation {
            axis,
            kind: LimitKind::HoldBoundExceeded,
            needed_mil: hold_true_mil,
            available_mil: Some(bound),
        }),
        Some(_) => None,
    }
}

fn end_revolution_for(target_clicks: i64, cpr: Option<u32>) -> Option<(u32, u32)> {
    cpr.and_then(|c| revolution_annotation(target_clicks, c))
}

/// `Direction::Up`/`Down` for elevation, `Left`/`Right` for windage, from the sign of
/// `delta_clicks`. Zero is reported as the positive variant (`Up`/`Right`) -- an
/// instruction to turn zero clicks has no natural direction of its own.
fn direction_for(axis: Axis, delta_clicks: i64) -> Direction {
    match (axis, delta_clicks < 0) {
        (Axis::Elevation, true) => Direction::Down,
        (Axis::Elevation, false) => Direction::Up,
        (Axis::Windage, true) => Direction::Left,
        (Axis::Windage, false) => Direction::Right,
    }
}

/// Per-axis DIAL-space quantization plus every input a [`Strategy`] needs for that axis,
/// computed exactly once per axis per [`plan_corrections`] call and shared across all three
/// strategies (MBA-1348) -- `target_clicks`/`travel_violation` are IDENTICAL for
/// `Strategy::DialAll` and `Strategy::Hybrid` by construction, never recomputed twice.
struct AxisComputation {
    axis: Axis,
    corr_true: f64,
    click_mil: f64,
    cf: f64,
    /// The whole-click target from zero, after travel-clamping (if travel is declared).
    target_clicks: i64,
    travel_violation: Option<LimitViolation>,
    /// `OpticProfile::turret_state`'s click-quantized position on this axis, or 0.
    state_clicks: i64,
    cpr: Option<u32>,
    /// The `HoldBounds` magnitude that bounds a POSITIVE `hold_mil` on this axis --
    /// `down_mil` for elevation, `left_mil` for windage (MBA-1348 review fix C1: NOT
    /// `up_mil`/`right_mil` -- correction-space and reticle-space are opposite-signed
    /// conventions, see the module's "Honesty" doc section, which cites `src/reticle.rs`).
    bound_for_positive_correction: Option<f64>,
    /// The `HoldBounds` magnitude that bounds a NEGATIVE `hold_mil` on this axis --
    /// `up_mil` for elevation, `right_mil` for windage. See `bound_for_positive_correction`.
    bound_for_negative_correction: Option<f64>,
}

// Private per-axis helper: every parameter is a distinct, independently-varying input this
// module's own tests exercise in isolation (click graduation, CF, travel, turret state,
// hold bounds, revolution count) -- collapsing them into a params struct would just move
// the same nine fields into a second, purely-local type with no independent meaning.
#[allow(clippy::too_many_arguments)]
fn compute_axis(
    axis: Axis,
    corr_true: f64,
    click: &ClickValue,
    cf: f64,
    travel: Option<&TravelLimits>,
    state_mil: Option<f64>,
    bound_for_positive_correction: Option<f64>,
    bound_for_negative_correction: Option<f64>,
    cpr: Option<u32>,
) -> AxisComputation {
    let click_mil = click_size_mil(click);
    let (target_clicks, _corr_dial, travel_violation) =
        quantize_and_clamp(axis, corr_true, click_mil, cf, travel);
    let state_clicks = state_mil.map_or(0, |mil| {
        quantize_angle(mil, &ClickValue { size: click_mil, base: ClickBase::Mil }).clicks
    });
    AxisComputation {
        axis,
        corr_true,
        click_mil,
        cf,
        target_clicks,
        travel_violation,
        state_clicks,
        cpr,
        bound_for_positive_correction,
        bound_for_negative_correction,
    }
}

/// Builds one axis's [`AxisInstruction`] for `strategy` from a shared [`AxisComputation`],
/// plus any [`LimitViolation`]s to fold into the plan's `limits_hit` and whether THIS axis,
/// under THIS strategy, satisfies the strategy's own feasibility promise (MBA-1348) -- see
/// the module's "Honesty" doc section for why `Strategy::Hybrid`'s feasibility ignores its
/// own (purely disclosed) travel violation.
fn build_axis_instruction(
    strategy: Strategy,
    ac: &AxisComputation,
    prefs: &Preferences,
) -> (AxisInstruction, Vec<LimitViolation>, bool) {
    match strategy {
        Strategy::DialAll => {
            let clicks = ac.target_clicks;
            let dial_true = clicks as f64 * ac.click_mil * ac.cf;
            let residual = ac.corr_true - dial_true;
            let delta_clicks = clicks - ac.state_clicks;
            let feasible = ac.travel_violation.is_none();
            let violations = ac.travel_violation.into_iter().collect();
            let instr = AxisInstruction {
                axis: ac.axis,
                direction: direction_for(ac.axis, delta_clicks),
                delta_clicks,
                target_clicks_from_zero: clicks,
                end_revolution: end_revolution_for(clicks, ac.cpr),
                dial_mil_true: dial_true,
                hold_mil: 0.0,
                residual_mil: residual,
            };
            (instr, violations, feasible)
        }
        Strategy::HoldAll => {
            let hold = ac.corr_true;
            let hold_violation = check_hold_bounds(
                ac.axis,
                hold,
                ac.bound_for_positive_correction,
                ac.bound_for_negative_correction,
                prefs.max_hold_mil,
            );
            let feasible = hold_violation.is_none();
            let target_clicks = 0_i64;
            let delta_clicks = target_clicks - ac.state_clicks;
            let violations = hold_violation.into_iter().collect();
            let instr = AxisInstruction {
                axis: ac.axis,
                direction: direction_for(ac.axis, delta_clicks),
                delta_clicks,
                target_clicks_from_zero: target_clicks,
                end_revolution: end_revolution_for(target_clicks, ac.cpr),
                dial_mil_true: 0.0,
                hold_mil: hold,
                residual_mil: 0.0,
            };
            (instr, violations, feasible)
        }
        Strategy::Hybrid => {
            let clicks = ac.target_clicks;
            let dial_true = clicks as f64 * ac.click_mil * ac.cf;
            let hold = ac.corr_true - dial_true;
            let residual = ac.corr_true - dial_true - hold;
            let hold_violation = check_hold_bounds(
                ac.axis,
                hold,
                ac.bound_for_positive_correction,
                ac.bound_for_negative_correction,
                prefs.max_hold_mil,
            );
            // Hybrid's OWN promise is that dial+hold reconstruct the correction exactly
            // (always true by construction -- `residual` above), so a travel violation on
            // its dial component is disclosed but does not gate ITS feasibility; only its
            // hold not fitting does. See the module's "Honesty" doc section.
            let feasible = hold_violation.is_none();
            let delta_clicks = clicks - ac.state_clicks;
            let mut violations: Vec<LimitViolation> = ac.travel_violation.into_iter().collect();
            violations.extend(hold_violation);
            let instr = AxisInstruction {
                axis: ac.axis,
                direction: direction_for(ac.axis, delta_clicks),
                delta_clicks,
                target_clicks_from_zero: clicks,
                end_revolution: end_revolution_for(clicks, ac.cpr),
                dial_mil_true: dial_true,
                hold_mil: hold,
                residual_mil: residual,
            };
            (instr, violations, feasible)
        }
    }
}

/// mil→linear at `range_m` by the small-angle approximation, RSS'd across both axes:
/// `sqrt((e/1000*range)^2 + (w/1000*range)^2)` (assumption 0).
fn residual_linear_at_range(e_res_mil: f64, w_res_mil: f64, range_m: f64) -> f64 {
    let e = e_res_mil / 1000.0 * range_m;
    let w = w_res_mil / 1000.0 * range_m;
    (e * e + w * w).sqrt()
}

fn build_plan(
    strategy: Strategy,
    elevation: &AxisComputation,
    windage: &AxisComputation,
    prefs: &Preferences,
    range_m: f64,
) -> DialPlan {
    let (e_instr, e_viol, e_feasible) = build_axis_instruction(strategy, elevation, prefs);
    let (w_instr, w_viol, w_feasible) = build_axis_instruction(strategy, windage, prefs);
    let mut limits_hit = e_viol;
    limits_hit.extend(w_viol);
    DialPlan {
        strategy,
        residual_linear_at_range_m: residual_linear_at_range(
            e_instr.residual_mil,
            w_instr.residual_mil,
            range_m,
        ),
        instructions: [e_instr, w_instr],
        feasible: e_feasible && w_feasible,
        limits_hit,
    }
}

/// Ascending sort key for `prefs.prefer_hold`: `false` (prefer dial) ranks
/// `DialAll < Hybrid < HoldAll`; `true` reverses that. See
/// `ranking_is_deterministic_and_preference_respected`.
fn preference_rank(strategy: Strategy, prefer_hold: bool) -> u8 {
    let dial_first = match strategy {
        Strategy::DialAll => 0,
        Strategy::Hybrid => 1,
        Strategy::HoldAll => 2,
    };
    if prefer_hold { 2 - dial_first } else { dial_first }
}

/// `Strategy`'s own declaration order (`DialAll, HoldAll, Hybrid`), the final ranking
/// tiebreak once residual and preference are both exhausted.
fn declaration_rank(strategy: Strategy) -> u8 {
    match strategy {
        Strategy::DialAll => 0,
        Strategy::HoldAll => 1,
        Strategy::Hybrid => 2,
    }
}

/// Turns a TRUE angular [`AngularCorrection`] into ranked, executable dial/hold/hybrid
/// plans for a real optic (MBA-1348) -- see the module's "`plan_corrections` and the CF
/// rule" and "Honesty" doc sections for the conventions this function follows in every arm.
///
/// Returns exactly three plans, one per [`Strategy`], in `DialPlanReportV1::plans`, ranked
/// best-first: `feasible: true` plans always sort before `feasible: false` ones (MBA-1348
/// review fix I4 -- `HoldAll`/`Hybrid` carry a residual of exactly `0.0` even when
/// infeasible, so residual alone cannot be the primary key without risking an unexecutable
/// `plans[0]`); among equally-feasible plans, ascending
/// [`DialPlan::residual_linear_at_range_m`] (via `f64::total_cmp`, so ranking never panics
/// regardless of input); ties are broken by `prefs.prefer_hold` (`false` ranks
/// `DialAll < Hybrid < HoldAll`, `true` reverses that), and any remaining tie by
/// `Strategy`'s own declaration order (`DialAll, HoldAll, Hybrid`) -- see
/// `ranking_is_deterministic_and_preference_respected` and
/// `infeasible_plans_never_outrank_a_feasible_one`.
///
/// `Err(OpticError)` when `optic.validate()` fails; when `correction`, `range_m`, or
/// `prefs.max_hold_mil` is not finite; when `range_m` / `max_hold_mil` is negative; or when
/// `elevation_cf` / `windage_cf` is zero, negative, or non-finite
/// ([`OpticError::NonPositiveTrackingFactor`], MBA-1348 review fix I5 -- the CF rule
/// divides by it, so a non-positive value is not merely implausible, it is a hard
/// arithmetic failure). Otherwise `elevation_cf`/`windage_cf` are trusted exactly as given,
/// like `crate::truing::scale_report_dial_values` -- this function does not re-apply
/// `crate::adjustment::tracking_cf_in_range`'s tighter `(0.5, 1.5)` plausibility band;
/// enforcing that band (if at all) remains a caller/CLI concern, matching that existing
/// precedent.
pub fn plan_corrections(
    correction: AngularCorrection,
    optic: &OpticProfile,
    range_m: f64,
    elevation_cf: f64,
    windage_cf: f64,
    prefs: &Preferences,
) -> Result<DialPlanReportV1, OpticError> {
    optic.validate()?;
    require_finite("correction.elevation_mil", correction.elevation_mil)?;
    require_finite("correction.windage_mil", correction.windage_mil)?;
    require_finite("range_m", range_m)?;
    require_non_negative("range_m", range_m)?;
    require_positive_tracking_factor("elevation_cf", elevation_cf)?;
    require_positive_tracking_factor("windage_cf", windage_cf)?;
    if let Some(max_hold) = prefs.max_hold_mil {
        require_finite("max_hold_mil", max_hold)?;
        require_non_negative("max_hold_mil", max_hold)?;
    }

    let (hold_up, hold_down, hold_left, hold_right) = match &optic.reticle_hold_bounds {
        Some(b) => (Some(b.up_mil), Some(b.down_mil), Some(b.left_mil), Some(b.right_mil)),
        None => (None, None, None, None),
    };

    // MBA-1348 review fix (C1, CRITICAL): a POSITIVE (up/right) correction-space
    // `hold_mil` is realized by a reticle mark on the OPPOSITE-named side -- a holdover
    // for an UP correction sits BELOW center (`down_mil`); a hold for a RIGHT correction
    // sits to the LEFT of center (`left_mil`). See the module's "Honesty" doc section
    // (which cites `src/reticle.rs:30-42`) for the full derivation. `hold_up`/`hold_down`/
    // `hold_left`/`hold_right` below are therefore passed CROSSED, not straight.
    let elevation = compute_axis(
        Axis::Elevation,
        correction.elevation_mil,
        &optic.elevation_click,
        elevation_cf,
        optic.elevation_travel.as_ref(),
        optic.turret_state.as_ref().map(|s| s.elevation_mil),
        hold_down,
        hold_up,
        optic.clicks_per_revolution,
    );
    let windage = compute_axis(
        Axis::Windage,
        correction.windage_mil,
        &optic.windage_click,
        windage_cf,
        optic.windage_travel.as_ref(),
        optic.turret_state.as_ref().map(|s| s.windage_mil),
        hold_left,
        hold_right,
        optic.clicks_per_revolution,
    );

    let mut plans = vec![
        build_plan(Strategy::DialAll, &elevation, &windage, prefs, range_m),
        build_plan(Strategy::HoldAll, &elevation, &windage, prefs, range_m),
        build_plan(Strategy::Hybrid, &elevation, &windage, prefs, range_m),
    ];
    // MBA-1348 review fix (I4): `feasible` must be the PRIMARY key. `HoldAll`/`Hybrid`
    // carry a residual of exactly 0.0 even when INFEASIBLE (residual is about
    // reconstruction, not achievability), so without this, an infeasible plan could rank
    // ahead of the one plan that actually works -- and Task 6's CLI surfaces `plans[0]` as
    // THE recommendation. `false < true`, so `!feasible` sorts feasible plans first.
    plans.sort_by(|a, b| {
        (!a.feasible)
            .cmp(&!b.feasible)
            .then_with(|| a.residual_linear_at_range_m.total_cmp(&b.residual_linear_at_range_m))
            .then_with(|| {
                preference_rank(a.strategy, prefs.prefer_hold)
                    .cmp(&preference_rank(b.strategy, prefs.prefer_hold))
            })
            .then_with(|| declaration_rank(a.strategy).cmp(&declaration_rank(b.strategy)))
    });

    Ok(DialPlanReportV1 {
        schema_version: DIAL_PLAN_SCHEMA_VERSION_V1,
        method: "dial_space_quantization_v1".to_string(),
        assumptions: vec![
            "Linear miss at range uses the small-angle approximation (mil / 1000 * range); \
             it is not exact at extreme angles."
                .to_string(),
            "Elevation and windage are planned independently; no cant-induced coupling \
             between axes is modeled."
                .to_string(),
            "Reticle holds are assumed continuous and unquantized, unlike turret clicks."
                .to_string(),
            "Travel limits and turret state are trusted exactly as declared in the optic \
             profile, not sensed or independently verified."
                .to_string(),
            "MOA-graduated clicks convert to milliradians using the locked printed-table \
             constant 3438, not the exact geometric 3437.7467."
                .to_string(),
        ],
        range_m,
        plans,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adjustment::ClickBase;

    /// A full, realistic profile: 0.1 mil clicks on both turrets, 10 clicks/revolution, a
    /// zero stop, 0.4 mil of down travel and 28.0 of up (typical of a zero-stopped
    /// elevation turret with most of its range above zero), +-6 mil of windage travel,
    /// dialed to zero on both axes, and modest hold bounds beyond the turrets' own range.
    fn baseline_profile() -> OpticProfile {
        OpticProfile {
            elevation_click: ClickValue { size: 0.1, base: ClickBase::Mil },
            windage_click: ClickValue { size: 0.1, base: ClickBase::Mil },
            clicks_per_revolution: Some(10),
            zero_stop: true,
            elevation_travel: Some(TravelLimits { down_mil: 0.4, up_mil: 28.0 }),
            windage_travel: Some(TravelLimits { down_mil: 6.0, up_mil: 6.0 }),
            turret_state: Some(TurretState { elevation_mil: 0.0, windage_mil: 0.0 }),
            reticle_hold_bounds: Some(HoldBounds {
                up_mil: 5.0,
                down_mil: 10.0,
                left_mil: 6.0,
                right_mil: 6.0,
            }),
        }
    }

    #[test]
    fn validate_accepts_a_full_realistic_profile() {
        assert_eq!(baseline_profile().validate(), Ok(()));
    }

    #[test]
    fn validate_rejects_negative_travel() {
        let mut down_negative = baseline_profile();
        down_negative.elevation_travel = Some(TravelLimits { down_mil: -0.4, up_mil: 28.0 });
        assert!(
            matches!(
                down_negative.validate(),
                Err(OpticError::NegativeLimit { field: "elevation_travel.down_mil", value })
                    if value == -0.4
            ),
            "{:?}",
            down_negative.validate()
        );

        let mut up_negative = baseline_profile();
        up_negative.windage_travel = Some(TravelLimits { down_mil: 6.0, up_mil: -6.0 });
        assert!(
            matches!(
                up_negative.validate(),
                Err(OpticError::NegativeLimit { field: "windage_travel.up_mil", value })
                    if value == -6.0
            ),
            "{:?}",
            up_negative.validate()
        );

        let mut hold_negative = baseline_profile();
        hold_negative.reticle_hold_bounds = Some(HoldBounds {
            up_mil: -5.0,
            down_mil: 10.0,
            left_mil: 6.0,
            right_mil: 6.0,
        });
        assert!(
            matches!(
                hold_negative.validate(),
                Err(OpticError::NegativeLimit { field: "reticle_hold_bounds.up_mil", .. })
            ),
            "{:?}",
            hold_negative.validate()
        );
    }

    #[test]
    fn validate_rejects_zero_clicks_per_revolution() {
        let mut zero_cpr = baseline_profile();
        zero_cpr.clicks_per_revolution = Some(0);
        assert!(matches!(
            zero_cpr.validate(),
            Err(OpticError::ZeroClicksPerRevolution)
        ));

        // Fencepost: 1 is the smallest ACCEPTED value.
        let mut one_cpr = baseline_profile();
        one_cpr.clicks_per_revolution = Some(1);
        assert_eq!(one_cpr.validate(), Ok(()));

        // A turret with no revolution marks at all is also fine.
        let mut no_cpr = baseline_profile();
        no_cpr.clicks_per_revolution = None;
        assert_eq!(no_cpr.validate(), Ok(()));
    }

    #[test]
    fn validate_rejects_non_positive_click_size() {
        // `ClickValue`'s fields are `pub`, so a zero or negative size is reachable by
        // direct construction even though `parse_click_value` already excludes it.
        for bad_size in [0.0, -0.1] {
            let mut elevation_bad = baseline_profile();
            elevation_bad.elevation_click.size = bad_size;
            assert!(
                matches!(
                    elevation_bad.validate(),
                    Err(OpticError::NonPositiveClickSize { field: "elevation_click.size", size })
                        if size == bad_size
                ),
                "size {bad_size}: {:?}",
                elevation_bad.validate()
            );

            let mut windage_bad = baseline_profile();
            windage_bad.windage_click.size = bad_size;
            assert!(
                matches!(
                    windage_bad.validate(),
                    Err(OpticError::NonPositiveClickSize { field: "windage_click.size", size })
                        if size == bad_size
                ),
                "size {bad_size}: {:?}",
                windage_bad.validate()
            );
        }

        // Fencepost: a very small but positive size is accepted.
        let mut tiny = baseline_profile();
        tiny.elevation_click.size = f64::MIN_POSITIVE;
        assert_eq!(tiny.validate(), Ok(()));
    }

    #[test]
    fn validate_rejects_non_finite_fields() {
        type Mutator = fn(&mut OpticProfile);
        let cases: &[(&str, Mutator)] = &[
            ("elevation_click.size", |p| p.elevation_click.size = f64::NAN),
            ("windage_click.size", |p| p.windage_click.size = f64::INFINITY),
            ("elevation_travel.down_mil", |p| {
                p.elevation_travel.as_mut().unwrap().down_mil = f64::NAN
            }),
            ("elevation_travel.up_mil", |p| {
                p.elevation_travel.as_mut().unwrap().up_mil = f64::INFINITY
            }),
            ("windage_travel.down_mil", |p| {
                p.windage_travel.as_mut().unwrap().down_mil = f64::NEG_INFINITY
            }),
            ("windage_travel.up_mil", |p| {
                p.windage_travel.as_mut().unwrap().up_mil = f64::NAN
            }),
            ("turret_state.elevation_mil", |p| {
                p.turret_state.as_mut().unwrap().elevation_mil = f64::NAN
            }),
            ("turret_state.windage_mil", |p| {
                p.turret_state.as_mut().unwrap().windage_mil = f64::INFINITY
            }),
            ("reticle_hold_bounds.up_mil", |p| {
                p.reticle_hold_bounds.as_mut().unwrap().up_mil = f64::NAN
            }),
            ("reticle_hold_bounds.down_mil", |p| {
                p.reticle_hold_bounds.as_mut().unwrap().down_mil = f64::NAN
            }),
            ("reticle_hold_bounds.left_mil", |p| {
                p.reticle_hold_bounds.as_mut().unwrap().left_mil = f64::NAN
            }),
            ("reticle_hold_bounds.right_mil", |p| {
                p.reticle_hold_bounds.as_mut().unwrap().right_mil = f64::NAN
            }),
        ];

        for (field, mutate) in cases {
            let mut profile = baseline_profile();
            mutate(&mut profile);
            let result = profile.validate();
            assert!(
                matches!(&result, Err(OpticError::NonFinite { field: f }) if f == field),
                "field {field}: expected Err(NonFinite {{ field: {field:?} }}), got {result:?}"
            );
        }
    }

    #[test]
    fn validate_rejects_turret_state_outside_travel() {
        let mut above_up = baseline_profile();
        above_up.turret_state = Some(TurretState { elevation_mil: 28.1, windage_mil: 0.0 });
        assert!(
            matches!(
                above_up.validate(),
                Err(OpticError::StateOutsideTravel {
                    axis: "elevation",
                    dialed_mil,
                    down_mil,
                    up_mil
                }) if dialed_mil == 28.1 && down_mil == 0.4 && up_mil == 28.0
            ),
            "{:?}",
            above_up.validate()
        );

        let mut below_down = baseline_profile();
        below_down.turret_state = Some(TurretState { elevation_mil: -0.5, windage_mil: 0.0 });
        assert!(matches!(
            below_down.validate(),
            Err(OpticError::StateOutsideTravel { axis: "elevation", .. })
        ));

        let mut windage_out = baseline_profile();
        windage_out.turret_state = Some(TurretState { elevation_mil: 0.0, windage_mil: 6.1 });
        assert!(matches!(
            windage_out.validate(),
            Err(OpticError::StateOutsideTravel { axis: "windage", .. })
        ));

        // Exactly at either boundary is accepted -- a closed interval.
        let mut at_boundary = baseline_profile();
        at_boundary.turret_state = Some(TurretState { elevation_mil: 28.0, windage_mil: -6.0 });
        assert_eq!(at_boundary.validate(), Ok(()));

        // No declared travel on an axis means nothing to be "outside" on that axis.
        let mut untethered = baseline_profile();
        untethered.elevation_travel = None;
        untethered.turret_state = Some(TurretState { elevation_mil: 999.0, windage_mil: 0.0 });
        assert_eq!(untethered.validate(), Ok(()));
    }

    #[test]
    fn revolution_annotation_matches_known_points_at_cpr_ten() {
        assert_eq!(revolution_annotation(0, 10), Some((0, 0)));
        assert_eq!(revolution_annotation(9, 10), Some((0, 9)));
        assert_eq!(revolution_annotation(10, 10), Some((1, 0)));
        assert_eq!(revolution_annotation(27, 10), Some((2, 7)));
    }

    #[test]
    fn revolution_annotation_reconstructs_the_input_over_a_range() {
        for cpr in [1_u32, 3, 10, 60] {
            for input in 0_i64..100 {
                let (revolutions, clicks_in_rev) = revolution_annotation(input, cpr)
                    .unwrap_or_else(|| panic!("expected Some for input {input}, cpr {cpr}"));
                assert!(
                    clicks_in_rev < cpr,
                    "clicks_in_rev {clicks_in_rev} must be < cpr {cpr}"
                );
                assert_eq!(
                    i64::from(revolutions) * i64::from(cpr) + i64::from(clicks_in_rev),
                    input,
                    "reconstruction broke at input {input}, cpr {cpr}"
                );
            }
        }
    }

    #[test]
    fn revolution_annotation_negative_input_is_none() {
        for clicks in [-1_i64, -7, -10, -27, i64::MIN] {
            assert_eq!(
                revolution_annotation(clicks, 10),
                None,
                "clicks_from_zero {clicks} must be None, not a misleading revolution count"
            );
        }
    }

    #[test]
    fn revolution_annotation_zero_clicks_per_revolution_is_none_not_a_panic() {
        assert_eq!(revolution_annotation(0, 0), None);
        assert_eq!(revolution_annotation(5, 0), None);
        assert_eq!(revolution_annotation(-5, 0), None);
    }

    // MBA-1348: OpticProfile's JSON wire form -- Task 5 stores this in profiles, so its
    // shape (including how absent turret data is represented) is pinned here, not left
    // to whatever serde's defaults happen to produce.
    #[test]
    fn optic_profile_round_trips_through_json() {
        let profile = baseline_profile();
        let json = serde_json::to_string(&profile).unwrap();
        let parsed: OpticProfile = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed, profile);
    }

    #[test]
    fn optic_profile_json_pins_click_fields_to_the_suffixed_string() {
        let profile = baseline_profile();
        let json = serde_json::to_value(&profile).unwrap();
        assert_eq!(json["elevation_click"], serde_json::json!("0.1mil"));
        assert_eq!(json["windage_click"], serde_json::json!("0.1mil"));
    }

    #[test]
    fn optic_profile_json_pins_option_field_presence_and_absence() {
        // Every optional field is `Some` in the baseline profile: present as its real
        // value, not omitted.
        let present = serde_json::to_value(baseline_profile()).unwrap();
        for field in [
            "clicks_per_revolution",
            "elevation_travel",
            "windage_travel",
            "turret_state",
            "reticle_hold_bounds",
        ] {
            let value = &present[field];
            assert!(
                !value.is_null(),
                "{field} should be present as a real value in the baseline profile, got {value:?}"
            );
        }

        // `None` serializes as an explicit JSON `null` -- there is no
        // `#[serde(skip_serializing_if = "Option::is_none")]` on this type, so a `None`
        // field key is still PRESENT in the object, just null-valued, not omitted.
        let mut profile = baseline_profile();
        profile.clicks_per_revolution = None;
        profile.elevation_travel = None;
        profile.windage_travel = None;
        profile.turret_state = None;
        profile.reticle_hold_bounds = None;
        let absent = serde_json::to_value(&profile).unwrap();
        for field in [
            "clicks_per_revolution",
            "elevation_travel",
            "windage_travel",
            "turret_state",
            "reticle_hold_bounds",
        ] {
            assert_eq!(
                absent.get(field),
                Some(&serde_json::Value::Null),
                "{field} should be a present-but-null key, not an omitted one"
            );
        }

        // And the None-shaped profile still round-trips.
        let json = serde_json::to_string(&profile).unwrap();
        let parsed: OpticProfile = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed, profile);
    }
}

// MBA-1348 Task 4: `plan_corrections` and its supporting types. Kept in its own test module,
// separate from `mod tests` above (Task 3's turret-model tests), matching this crate's
// existing practice of splitting a large file's tests by concern (e.g.
// `error_budget::window_helper_tests`).
#[cfg(test)]
mod plan_corrections_tests {
    use super::*;
    use crate::adjustment::ClickBase;

    /// Same numbers as `tests::baseline_profile` (Task 3) -- 0.1 mil clicks on both
    /// turrets, 10 clicks/revolution, a zero stop, 0.4 mil of down travel and 28.0 of up,
    /// +-6 mil of windage travel, dialed to zero on both axes, and hold bounds of 5 up /
    /// 10 down / 6 left / 6 right -- redefined locally (rather than reused across modules)
    /// so this module doesn't depend on `tests`' private helper visibility.
    fn baseline_profile() -> OpticProfile {
        OpticProfile {
            elevation_click: ClickValue { size: 0.1, base: ClickBase::Mil },
            windage_click: ClickValue { size: 0.1, base: ClickBase::Mil },
            clicks_per_revolution: Some(10),
            zero_stop: true,
            elevation_travel: Some(TravelLimits { down_mil: 0.4, up_mil: 28.0 }),
            windage_travel: Some(TravelLimits { down_mil: 6.0, up_mil: 6.0 }),
            turret_state: Some(TurretState { elevation_mil: 0.0, windage_mil: 0.0 }),
            reticle_hold_bounds: Some(HoldBounds {
                up_mil: 5.0,
                down_mil: 10.0,
                left_mil: 6.0,
                right_mil: 6.0,
            }),
        }
    }

    fn plan_for(strategy: Strategy, report: &DialPlanReportV1) -> DialPlan {
        report
            .plans
            .iter()
            .find(|p| p.strategy == strategy)
            .unwrap_or_else(|| panic!("no {strategy:?} plan in report"))
            .clone()
    }

    // Spec §7 acceptance criterion: "exact-click dope produces zero residual."
    //
    // `elevation_mil` is constructed as `23.0 * 0.1`, NOT the literal `2.3`, so it is
    // BIT-IDENTICAL to `clicks as f64 * click.size` after quantization -- matching
    // `crate::adjustment::tests::quantize_exact_multiples_have_residual_exactly_zero`'s own
    // established pattern for the same reason: `2.3_f64 - 23.0 * 0.1 == -4.440892098500626e-16`,
    // NOT zero, so the decimal literal and the click-multiple product are numerically equal
    // but not bit-identical. This is a hand-verified floating-point fact (checked directly
    // in Python's IEEE-754 binary64 floats, identical semantics to Rust `f64`), not a
    // stylistic choice -- the brief's own "residual_mil == 0.0 bit-exact" requirement is
    // only achievable this way.
    #[test]
    fn exact_click_dope_has_residual_exactly_zero() {
        let corr = AngularCorrection { elevation_mil: 23.0 * 0.1, windage_mil: 0.0 };
        let optic = baseline_profile();
        let report =
            plan_corrections(corr, &optic, 600.0, 1.0, 1.0, &Preferences::default()).unwrap();
        let dial_all = plan_for(Strategy::DialAll, &report);
        assert_eq!(dial_all.instructions[0].target_clicks_from_zero, 23);
        assert_eq!(dial_all.instructions[0].residual_mil, 0.0, "must be bit-exact zero");
        assert_eq!(dial_all.instructions[0].residual_mil.to_bits(), 0.0_f64.to_bits());
        assert!(dial_all.feasible, "{dial_all:?}");
        assert!(dial_all.limits_hit.is_empty(), "{dial_all:?}");
    }

    // Spec §7: "fractional-click dope reports the chosen rounding and remaining
    // angular/linear error at range."
    //
    // Hand-derived (and Python-verified against IEEE-754 binary64 arithmetic): 2.34 / 0.1 =
    // 23.4, rounds to 23 clicks; residual = 2.34 - 23*0.1 = 0.03999999999999959 (~0.04,
    // matching `2.34 - 2.3` to within 1e-12 as the brief specifies); linear error at 600 m
    // = |residual| / 1000 * 600 = 0.023999999999999754 (~0.024, within 1e-9).
    #[test]
    fn fractional_click_reports_rounding_and_linear_error() {
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        let optic = baseline_profile();
        let report =
            plan_corrections(corr, &optic, 600.0, 1.0, 1.0, &Preferences::default()).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        assert_eq!(dial_all.instructions[0].target_clicks_from_zero, 23);
        let expected_residual = 2.34 - 2.3;
        assert!(
            (dial_all.instructions[0].residual_mil - expected_residual).abs() < 1e-12,
            "residual {} vs hand-derived {}",
            dial_all.instructions[0].residual_mil,
            expected_residual
        );
        assert!(
            (dial_all.residual_linear_at_range_m - 0.024).abs() < 1e-9,
            "residual_linear_at_range_m = {}",
            dial_all.residual_linear_at_range_m
        );

        let hybrid = plan_for(Strategy::Hybrid, &report);
        assert_eq!(
            hybrid.instructions[0].residual_mil, 0.0,
            "Hybrid's residual is an identity by construction, always bit-exact zero"
        );
        assert!(
            (hybrid.instructions[0].hold_mil - 0.04).abs() < 1e-9,
            "hold_mil = {}",
            hybrid.instructions[0].hold_mil
        );
    }

    // Spec §7: "MIL and MOA produce physically equivalent results."
    #[test]
    fn mil_and_moa_optics_are_physically_equivalent() {
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        let prefs = Preferences::default();

        let mil_optic = baseline_profile();
        let mil_report = plan_corrections(corr, &mil_optic, 100.0, 1.0, 1.0, &prefs).unwrap();

        let mut moa_optic = baseline_profile();
        moa_optic.elevation_click = ClickValue { size: 0.25, base: ClickBase::Moa };
        let moa_report = plan_corrections(corr, &moa_optic, 100.0, 1.0, 1.0, &prefs).unwrap();

        for (label, report, click) in [
            ("mil", &mil_report, &mil_optic.elevation_click),
            ("moa", &moa_report, &moa_optic.elevation_click),
        ] {
            let hybrid = plan_for(Strategy::Hybrid, report);
            let e = &hybrid.instructions[0];
            assert!(
                (e.dial_mil_true + e.hold_mil - corr.elevation_mil).abs() < 1e-12,
                "{label}: dial_true {} + hold {} should reconstruct corr {} to 1e-12",
                e.dial_mil_true,
                e.hold_mil,
                corr.elevation_mil
            );
            assert_eq!(e.residual_mil, 0.0, "{label}: Hybrid residual must be bit-exact zero");

            let dial_all = plan_for(Strategy::DialAll, report);
            let click_mil = click_size_mil(click);
            assert!(
                dial_all.instructions[0].residual_mil.abs() <= click_mil / 2.0,
                "{label}: DialAll residual {} must be within half a click ({})",
                dial_all.instructions[0].residual_mil,
                click_mil / 2.0
            );
        }
    }

    // The brief's pinned worked example -- hand-derived numbers, verified against Rust's
    // exact IEEE-754 f64 semantics before writing this test:
    //   corr_dial  = 5.0 / 0.98               = 5.1020408163265305
    //   clicks     = round(corr_dial / 0.1)   = 51
    //   dial_true  = 51 * 0.1 * 0.98          = 4.998
    //   hybrid hold = 5.0 - 4.998             = 0.0019999999999997797  (~0.002, 1e-12)
    // If the CF multiply/divide directions were swapped, or Hybrid held the DIAL-space
    // remainder instead of the TRUE one, these numbers would come out different and wrong
    // -- see this task's required fault-injection verification.
    #[test]
    fn cf_dial_space_worked_example() {
        let corr = AngularCorrection { elevation_mil: 5.0, windage_mil: 0.0 };
        let optic = baseline_profile();
        let report =
            plan_corrections(corr, &optic, 100.0, 0.98, 1.0, &Preferences::default()).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        let e = &dial_all.instructions[0];
        assert_eq!(e.target_clicks_from_zero, 51, "51 clicks: 5.0/0.98 = 5.10204... -> round");
        assert!(
            (e.dial_mil_true - 4.998).abs() < 1e-12,
            "dial_mil_true = {} (expected 51*0.1*0.98 = 4.998)",
            e.dial_mil_true
        );

        let hybrid = plan_for(Strategy::Hybrid, &report);
        let eh = &hybrid.instructions[0];
        assert_eq!(eh.target_clicks_from_zero, 51);
        assert!(
            (eh.hold_mil - 0.002).abs() < 1e-12,
            "hybrid hold_mil = {} (expected 5.0 - 4.998 = 0.002)",
            eh.hold_mil
        );
        assert_eq!(eh.residual_mil, 0.0, "Hybrid residual must be bit-exact zero");
        assert!(hybrid.feasible, "{hybrid:?}");
    }

    // Spec §7: "revolution and click counts are correct across zero-stop and travel
    // limits" + the brief's specific cpr-10 / 0.4-down / 28-up worked scenario.
    #[test]
    fn revolutions_and_zero_stop() {
        let optic = baseline_profile(); // cpr 10, up 28.0, down 0.4, click 0.1 mil
        let prefs = Preferences::default();

        // --- 27 clicks: end_revolution == Some((2, 7)) ---
        // corr constructed as 27.0 * 0.1 for the same bit-exactness reason as the
        // exact-click test above.
        let up_corr = AngularCorrection { elevation_mil: 27.0 * 0.1, windage_mil: 0.0 };
        let up_report = plan_corrections(up_corr, &optic, 100.0, 1.0, 1.0, &prefs).unwrap();
        let dial_all = plan_for(Strategy::DialAll, &up_report);
        assert_eq!(dial_all.instructions[0].target_clicks_from_zero, 27);
        assert_eq!(dial_all.instructions[0].end_revolution, Some((2, 7)));
        assert!(dial_all.feasible, "{dial_all:?}");

        // --- down 1.0 mil needs -10 clicks; only 0.4 mil (4 clicks) of down travel ---
        let down_corr = AngularCorrection { elevation_mil: -1.0, windage_mil: 0.0 };
        let down_report = plan_corrections(down_corr, &optic, 100.0, 1.0, 1.0, &prefs).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &down_report);
        let e = &dial_all.instructions[0];
        assert_eq!(e.target_clicks_from_zero, -4, "clamped to the 0.4 mil / 0.1 mil = 4 clicks available");
        assert!(!dial_all.feasible, "DialAll must be infeasible when travel-clamped");
        assert_eq!(dial_all.limits_hit.len(), 1);
        assert!(matches!(
            dial_all.limits_hit[0],
            LimitViolation {
                axis: Axis::Elevation,
                kind: LimitKind::TravelExceeded,
                needed_mil,
                available_mil: Some(available_mil),
            } if (needed_mil - (-1.0)).abs() < 1e-12 && available_mil == 0.4
        ), "{:?}", dial_all.limits_hit[0]);

        let hybrid = plan_for(Strategy::Hybrid, &down_report);
        let eh = &hybrid.instructions[0];
        assert_eq!(eh.target_clicks_from_zero, -4, "Hybrid dials the same clamped -4 clicks");
        assert!(
            eh.hold_mil.is_sign_negative() && (eh.hold_mil - (-0.6)).abs() < 1e-9,
            "Hybrid holds the rest: hold_mil = {} (expected -0.6)",
            eh.hold_mil
        );
        assert_eq!(eh.residual_mil, 0.0, "Hybrid residual must be bit-exact zero even when clamped");
        // The travel violation is still recorded on Hybrid (disclosure)...
        assert!(
            hybrid.limits_hit.iter().any(|v| matches!(v.kind, LimitKind::TravelExceeded)),
            "{:?}",
            hybrid.limits_hit
        );
        // ...but Hybrid's own feasibility depends ONLY on whether the hold fits bounds, NOT
        // on the travel clamp -- this is the crux of the brief's "feasible iff hold fits
        // bounds". The hold is -0.6 (a DOWN hold), which consumes baseline's up_mil bound
        // (5.0, the mark ABOVE center a downward correction's hold sits at -- see the
        // module's "Honesty" doc section) -- |-0.6| fits comfortably within 5.0.
        assert!(
            hybrid.feasible,
            "Hybrid must remain feasible: the hold (-0.6) fits within the 5.0 mil up_mil \
             bound, even though its dial component was travel-clamped: {hybrid:?}"
        );
    }

    // Spec §7: "an infeasible request is reported as infeasible rather than silently
    // clamped" -- the brief's own acceptance test name.
    #[test]
    fn infeasible_is_reported_never_silently_clamped() {
        let optic = baseline_profile(); // down travel 0.4 mil
        let corr = AngularCorrection { elevation_mil: -1.0, windage_mil: 0.0 };
        let report =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default()).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        assert!(!dial_all.feasible);
        let e = &dial_all.instructions[0];
        // The clamped dial only executes -4 clicks * 0.1 mil * 1.0 cf = -0.4 true mil, NOT
        // the requested -1.0. The reported residual must equal the REAL shortfall against
        // what was actually executed, never something smaller that pretends the clamp
        // didn't happen.
        // Review fix M4: tightened from a 1e-12 tolerance to bit-exact -- hand-verified
        // (Python, identical IEEE-754 semantics) that `-4.0 * 0.1 == -0.4` and
        // `-1.0 - (-0.4) == -0.6` are both exact for these specific numbers, so an
        // approximate assert here was strictly weaker than what the numbers actually give.
        let clamped_dial_true = e.target_clicks_from_zero as f64 * 0.1 * 1.0;
        assert_eq!(clamped_dial_true, -0.4);
        let honest_residual = corr.elevation_mil - clamped_dial_true;
        assert_eq!(
            e.residual_mil, honest_residual,
            "residual_mil must equal corr_true - clamped_dial_true exactly, not a smaller, \
             optimistic number"
        );
        assert!(
            e.residual_mil.abs() > 0.5,
            "a silently-optimistic implementation might under-report this; the real miss \
             is large (-0.6 mil): residual_mil = {}",
            e.residual_mil
        );

        // No plan anywhere in the ranked list is feasible: true while ALSO carrying a
        // recorded violation that its own strategy's feasibility depends on (DialAll's own
        // travel violation, HoldAll's own hold violation). This is a general property of
        // the whole report, not just DialAll.
        for plan in &report.plans {
            match plan.strategy {
                Strategy::DialAll | Strategy::HoldAll => {
                    assert_eq!(
                        plan.feasible,
                        plan.limits_hit.is_empty(),
                        "{:?}: DialAll/HoldAll feasibility must exactly track limits_hit",
                        plan.strategy
                    );
                }
                Strategy::Hybrid => {} // see revolutions_and_zero_stop: gated on hold only
            }
        }
    }

    // Spec §7 + the brief's own ranking rule: ascending residual_linear_at_range_m, ties
    // broken by `prefer_hold` (false => DialAll < Hybrid < HoldAll, true => reverse), final
    // tie by Strategy's declaration order. No float sort_by panic regardless (total_cmp).
    #[test]
    fn ranking_is_deterministic_and_preference_respected() {
        // --- A full 3-way tie: an exact-click correction makes ALL THREE strategies'
        // residual (and therefore residual_linear_at_range_m) exactly 0.0 -- DialAll lands
        // on an exact click (no quantization error), and HoldAll/Hybrid are ALWAYS exactly
        // 0.0 by construction regardless of the correction. This isolates the preference
        // tiebreak completely: nothing here is decided by the primary residual key.
        let optic = baseline_profile();
        let corr = AngularCorrection { elevation_mil: 23.0 * 0.1, windage_mil: 0.0 };

        let prefer_dial = Preferences { prefer_hold: false, max_hold_mil: None };
        let report_dial =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &prefer_dial).unwrap();
        let strategies: Vec<Strategy> = report_dial.plans.iter().map(|p| p.strategy).collect();
        assert_eq!(
            strategies,
            vec![Strategy::DialAll, Strategy::Hybrid, Strategy::HoldAll],
            "prefer_hold=false must rank DialAll < Hybrid < HoldAll when fully tied"
        );

        let prefer_hold = Preferences { prefer_hold: true, max_hold_mil: None };
        let report_hold =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &prefer_hold).unwrap();
        let strategies: Vec<Strategy> = report_hold.plans.iter().map(|p| p.strategy).collect();
        assert_eq!(
            strategies,
            vec![Strategy::HoldAll, Strategy::Hybrid, Strategy::DialAll],
            "prefer_hold=true must reverse the order to HoldAll < Hybrid < DialAll"
        );

        // --- A NON-tied scenario: DialAll's residual (~0.04) is strictly worse than
        // HoldAll/Hybrid's (always exactly 0.0), so the PRIMARY ascending-residual key must
        // place DialAll last regardless of preference -- preference only ever breaks ties,
        // it never overrides a genuine residual difference.
        let frac_corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        for prefs in [prefer_dial, prefer_hold] {
            let report =
                plan_corrections(frac_corr, &optic, 100.0, 1.0, 1.0, &prefs).unwrap();
            assert_eq!(
                report.plans.last().unwrap().strategy,
                Strategy::DialAll,
                "DialAll's nonzero residual must always rank it last, prefer_hold={}",
                prefs.prefer_hold
            );
            // Ranked ascending by residual_linear_at_range_m, non-decreasing throughout.
            for pair in report.plans.windows(2) {
                assert!(
                    pair[0].residual_linear_at_range_m <= pair[1].residual_linear_at_range_m,
                    "{:?}",
                    report.plans
                );
            }
        }
    }

    // Spec §2 (every report carries method + assumptions) + the Plan-A lesson: a
    // `contains("word")` check alone is satisfiable by a DIFFERENT sentence, so this pins
    // length AND full-string equality at every index (mirrors card.rs's own
    // `report_carries_method_and_all_five_assumptions`) -- a consumer that quotes
    // `assumptions[2]` must keep getting the unquantized-holds sentence, not whatever a
    // later edit shuffled into that slot.
    #[test]
    fn report_carries_method_and_all_five_assumptions() {
        let optic = baseline_profile();
        let corr = AngularCorrection { elevation_mil: 2.3, windage_mil: 0.0 };
        let report =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default()).unwrap();

        assert_eq!(report.schema_version, DIAL_PLAN_SCHEMA_VERSION_V1);
        assert_eq!(report.method, "dial_space_quantization_v1");
        assert_eq!(report.assumptions.len(), 5, "{:?}", report.assumptions);
        assert_eq!(
            report.assumptions[0],
            "Linear miss at range uses the small-angle approximation (mil / 1000 * range); it is not exact at extreme angles."
        );
        assert_eq!(
            report.assumptions[1],
            "Elevation and windage are planned independently; no cant-induced coupling between axes is modeled."
        );
        assert_eq!(
            report.assumptions[2],
            "Reticle holds are assumed continuous and unquantized, unlike turret clicks."
        );
        assert_eq!(
            report.assumptions[3],
            "Travel limits and turret state are trusted exactly as declared in the optic profile, not sensed or independently verified."
        );
        assert_eq!(
            report.assumptions[4],
            "MOA-graduated clicks convert to milliradians using the locked printed-table constant 3438, not the exact geometric 3437.7467."
        );
    }

    // The brief's own acceptance test: turret_state shifts delta_clicks but never
    // target_clicks_from_zero.
    #[test]
    fn turret_state_shifts_delta_but_not_target() {
        let corr = AngularCorrection { elevation_mil: 23.0 * 0.1, windage_mil: 0.0 };
        let prefs = Preferences::default();

        let zeroed = baseline_profile(); // turret_state elevation_mil: 0.0
        let report_zeroed = plan_corrections(corr, &zeroed, 100.0, 1.0, 1.0, &prefs).unwrap();
        let dial_zeroed = plan_for(Strategy::DialAll, &report_zeroed);

        let mut dialed = baseline_profile();
        dialed.turret_state = Some(TurretState { elevation_mil: 1.0, windage_mil: 0.0 });
        let report_dialed = plan_corrections(corr, &dialed, 100.0, 1.0, 1.0, &prefs).unwrap();
        let dial_dialed = plan_for(Strategy::DialAll, &report_dialed);

        assert_eq!(dial_zeroed.instructions[0].target_clicks_from_zero, 23);
        assert_eq!(dial_dialed.instructions[0].target_clicks_from_zero, 23,
            "target_clicks_from_zero must NOT move just because turret_state changed");
        assert_eq!(dial_zeroed.instructions[0].delta_clicks, 23);
        assert_eq!(
            dial_dialed.instructions[0].delta_clicks, 13,
            "1.0 mil dialed = 10 clicks of state; delta_clicks must drop by exactly 10 (23 -> 13)"
        );

        // Review fix I3: `end_revolution` must be computed from `target_clicks_from_zero`
        // (23, unchanged by turret_state -- cpr 10 -> revolution_annotation(23, 10) ==
        // (2, 3)), never from `delta_clicks` (23 for the zeroed case but 13 for the dialed
        // one, which at cpr 10 would give the WRONG (1, 3)). A silent switch to delta would
        // misreport a real, dialed-scope shooter's revolution count.
        assert_eq!(dial_zeroed.instructions[0].end_revolution, Some((2, 3)));
        assert_eq!(
            dial_dialed.instructions[0].end_revolution,
            Some((2, 3)),
            "end_revolution must stay (2, 3) with turret_state dialed, not shift to \
             delta_clicks=13's (1, 3)"
        );
    }

    // ---- Beyond the brief: NoTravelData / NoHoldBoundData disclosure ----
    //
    // Neither LimitKind is exercised by the brief's nine named tests, but both appear in
    // the `LimitKind` interface the brief itself specifies, so their behavior is a real
    // part of this task's surface, not an incidental implementation detail. Absence of
    // declared data is a WEAKER claim than a known, exceeded limit (see the module's
    // "Honesty" doc section): it is disclosed in `limits_hit`, and it gates exactly the
    // strategies whose OWN feasibility promise depends on that data.

    #[test]
    fn missing_travel_data_is_disclosed_and_gates_dial_all_but_not_hybrid() {
        let mut optic = baseline_profile();
        optic.elevation_travel = None; // no travel data at all on this axis
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        let report =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default()).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        assert!(
            dial_all.limits_hit.iter().any(|v| matches!(
                v,
                LimitViolation { axis: Axis::Elevation, kind: LimitKind::NoTravelData, available_mil: None, .. }
            )),
            "{:?}",
            dial_all.limits_hit
        );
        assert!(!dial_all.feasible, "DialAll cannot affirm feasibility without travel data");

        let hybrid = plan_for(Strategy::Hybrid, &report);
        assert!(
            hybrid.limits_hit.iter().any(|v| matches!(v.kind, LimitKind::NoTravelData)),
            "Hybrid still discloses the missing data: {:?}",
            hybrid.limits_hit
        );
        assert!(
            hybrid.feasible,
            "Hybrid's feasibility does not depend on travel data at all, only its hold \
             fitting bounds (it does, here): {hybrid:?}"
        );

        // A trivial (zero) correction needs no travel data at all, so nothing is disclosed.
        let zero_corr = AngularCorrection { elevation_mil: 0.0, windage_mil: 0.0 };
        let zero_report =
            plan_corrections(zero_corr, &optic, 100.0, 1.0, 1.0, &Preferences::default())
                .unwrap();
        let zero_dial_all = plan_for(Strategy::DialAll, &zero_report);
        assert!(zero_dial_all.limits_hit.is_empty(), "{:?}", zero_dial_all.limits_hit);
        assert!(zero_dial_all.feasible);
    }

    #[test]
    fn missing_hold_bound_data_is_disclosed_and_gates_hold_all_and_hybrid() {
        let mut optic = baseline_profile();
        optic.reticle_hold_bounds = None; // no hold bound data at all
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        let prefs = Preferences::default(); // max_hold_mil also None
        let report = plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &prefs).unwrap();

        let hold_all = plan_for(Strategy::HoldAll, &report);
        assert!(
            hold_all.limits_hit.iter().any(|v| matches!(
                v,
                LimitViolation { axis: Axis::Elevation, kind: LimitKind::NoHoldBoundData, available_mil: None, .. }
            )),
            "{:?}",
            hold_all.limits_hit
        );
        assert!(!hold_all.feasible);

        let hybrid = plan_for(Strategy::Hybrid, &report);
        assert!(
            hybrid.limits_hit.iter().any(|v| matches!(v.kind, LimitKind::NoHoldBoundData)),
            "{:?}",
            hybrid.limits_hit
        );
        assert!(!hybrid.feasible, "Hybrid's hold ALSO cannot be verified here: {hybrid:?}");

        // But `max_hold_mil` alone is enough to make the hold checkable again -- 3.0 is
        // comfortably above the 2.34 mil hold actually needed, so this must fit, not merely
        // become checkable-and-still-exceeded.
        let capped = Preferences { prefer_hold: false, max_hold_mil: Some(3.0) };
        let capped_report = plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &capped).unwrap();
        let capped_hold_all = plan_for(Strategy::HoldAll, &capped_report);
        assert!(
            capped_hold_all.limits_hit.is_empty(),
            "2.34 fits within the 3.0 cap, so nothing should be recorded at all: {:?}",
            capped_hold_all.limits_hit
        );
        assert!(capped_hold_all.feasible, "{:?}", capped_hold_all);
    }

    #[test]
    fn max_clicks_within_handles_exact_and_fractional_boundaries() {
        // Exact multiples: floating-point noise in the division must not discard a click.
        assert_eq!(max_clicks_within(0.4, 0.1), 4);
        assert_eq!(max_clicks_within(28.0, 0.1), 280);
        assert_eq!(max_clicks_within(6.0, 0.1), 60);
        // A genuinely fractional boundary floors -- never rounds up past a hard limit.
        assert_eq!(max_clicks_within(0.45, 0.1), 4);
        assert_eq!(max_clicks_within(0.05, 0.1), 0);
    }

    #[test]
    fn preferences_default_prefers_dial_with_no_hold_cap() {
        let p = Preferences::default();
        assert!(!p.prefer_hold);
        assert_eq!(p.max_hold_mil, None);
    }

    #[test]
    fn plan_corrections_rejects_an_invalid_profile() {
        let mut optic = baseline_profile();
        optic.elevation_click.size = 0.0; // rejected by OpticProfile::validate
        let corr = AngularCorrection { elevation_mil: 1.0, windage_mil: 0.0 };
        let result = plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default());
        assert!(matches!(result, Err(OpticError::NonPositiveClickSize { .. })), "{result:?}");
    }

    // ==== 2026-08 review fixes ====
    //
    // A reviewer hand-verified the CF rule and found the arithmetic core correct, but
    // found the hold-bound mapping inverted on BOTH axes (C1, critical) and several gaps
    // this task's original 14 tests could not have caught because nothing exercised an
    // asymmetric hold bound, a nonzero windage correction, HoldAll's own hold value, an
    // infeasible-but-zero-residual plan, a non-positive CF, or a non-unit-CF travel
    // violation. Every test below is a discriminating test for one of those gaps: each is
    // built so the specific bug it targets makes it fail, not just a generic re-assertion
    // of already-covered behavior.

    // C1 (CRITICAL): `hold_mil` is correction-space (+ = up/right, matching a dial's own
    // convention); `HoldBounds` is `crate::reticle`'s RETICLE-space (`src/reticle.rs:30-42`:
    // `down_mil` positive BELOW center, `right_mil` positive to the shooter's RIGHT). These
    // are OPPOSITELY signed on both axes: a positive (up) correction's hold mark sits BELOW
    // center (consumes `down_mil`, not `up_mil`) -- exactly how a BDC reticle's long-range
    // holdover marks sit below center, never above. A positive (right) correction's hold
    // mark sits to the LEFT of center (consumes `left_mil`, not `right_mil`) by the mirrored
    // argument. See the module's "Honesty" doc section for the full derivation.
    //
    // Baseline's elevation hold bounds are deliberately asymmetric (5.0 up / 10.0 down), so
    // this test actually discriminates the two possible mappings: a +7.0 mil correction
    // fits (7.0 <= 10.0 available below center) while a -7.0 mil correction does not
    // (7.0 > 5.0 available above center) -- the OLD, inverted mapping gets BOTH of these
    // backwards (it would have accepted -7.0 and rejected +7.0). Windage needs its own
    // asymmetric profile since baseline's windage bounds (6.0/6.0) are symmetric and
    // therefore cannot discriminate the mapping on their own.
    #[test]
    fn hold_bound_mapping_matches_reticle_space_not_correction_space() {
        let optic = baseline_profile(); // elevation hold bounds: up 5.0 / down 10.0
        let prefs = Preferences::default();

        let up_report = plan_corrections(
            AngularCorrection { elevation_mil: 7.0, windage_mil: 0.0 },
            &optic,
            100.0,
            1.0,
            1.0,
            &prefs,
        )
        .unwrap();
        let up_hold_all = plan_for(Strategy::HoldAll, &up_report);
        assert!(
            up_hold_all.feasible,
            "+7.0 mil (up) consumes down_mil=10.0 (available) and must fit: {up_hold_all:?}"
        );
        assert!(up_hold_all.limits_hit.is_empty(), "{:?}", up_hold_all.limits_hit);

        let down_report = plan_corrections(
            AngularCorrection { elevation_mil: -7.0, windage_mil: 0.0 },
            &optic,
            100.0,
            1.0,
            1.0,
            &prefs,
        )
        .unwrap();
        let down_hold_all = plan_for(Strategy::HoldAll, &down_report);
        assert!(
            !down_hold_all.feasible,
            "-7.0 mil (down) consumes up_mil=5.0 (available) and must NOT fit: {down_hold_all:?}"
        );
        let down_violation = down_hold_all
            .limits_hit
            .iter()
            .find(|v| v.axis == Axis::Elevation && matches!(v.kind, LimitKind::HoldBoundExceeded))
            .unwrap_or_else(|| panic!("{:?}", down_hold_all.limits_hit));
        assert_eq!(down_violation.available_mil, Some(5.0));

        // Windage: same asymmetry, mirrored -- a +right correction consumes left_mil, a
        // -left correction consumes right_mil.
        let mut windage_optic = baseline_profile();
        {
            let bounds = windage_optic.reticle_hold_bounds.as_mut().unwrap();
            bounds.left_mil = 8.0;
            bounds.right_mil = 3.0;
        }

        let right_report = plan_corrections(
            AngularCorrection { elevation_mil: 0.0, windage_mil: 7.0 },
            &windage_optic,
            100.0,
            1.0,
            1.0,
            &prefs,
        )
        .unwrap();
        let right_hold_all = plan_for(Strategy::HoldAll, &right_report);
        assert!(
            right_hold_all.feasible,
            "+7.0 mil (right) consumes left_mil=8.0 (available) and must fit: {right_hold_all:?}"
        );

        let left_report = plan_corrections(
            AngularCorrection { elevation_mil: 0.0, windage_mil: -7.0 },
            &windage_optic,
            100.0,
            1.0,
            1.0,
            &prefs,
        )
        .unwrap();
        let left_hold_all = plan_for(Strategy::HoldAll, &left_report);
        assert!(
            !left_hold_all.feasible,
            "-7.0 mil (left) consumes right_mil=3.0 (available) and must NOT fit: {left_hold_all:?}"
        );
        let left_violation = left_hold_all
            .limits_hit
            .iter()
            .find(|v| v.axis == Axis::Windage && matches!(v.kind, LimitKind::HoldBoundExceeded))
            .unwrap_or_else(|| panic!("{:?}", left_hold_all.limits_hit));
        assert_eq!(left_violation.available_mil, Some(3.0));
    }

    // The "danger case" the review specifically flagged: a travel-clamped Hybrid whose
    // hold the glass genuinely cannot support must be `feasible: false` -- the inverted
    // mapping could certify exactly that as feasible (whichever bound is larger always
    // "fits", regardless of which side the hold is actually on).
    #[test]
    fn travel_clamped_hybrid_with_an_unsupportable_hold_is_infeasible() {
        let optic = baseline_profile(); // down travel 0.4 mil; hold bounds up 5.0 / down 10.0
        // -7.4 mil needs -74 clicks; only 4 clicks (0.4 mil) of down travel exist, so
        // Hybrid clamps to -4 clicks (-0.4 true mil) and must hold the remaining exactly
        // -7.0 (hand-verified: -4.0*0.1 == -0.4 exactly, -7.4-(-0.4) == -7.0 exactly).
        // This magnitude is chosen to land BETWEEN the two candidate bounds and so
        // genuinely discriminate them: the CORRECT mapping checks a down-hold against
        // up_mil=5.0 (7.0 > 5.0 -> infeasible, correctly, since the glass cannot show 7.0
        // mil of upward hold), while the OLD, inverted mapping checked it against
        // down_mil=10.0 (7.0 <= 10.0 -> would have wrongly certified this feasible). See
        // the fault-injection transcript in the task report: re-inverting the mapping
        // makes this test fail, exactly as this comment predicts.
        let corr = AngularCorrection { elevation_mil: -7.4, windage_mil: 0.0 };
        let report =
            plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default()).unwrap();
        let hybrid = plan_for(Strategy::Hybrid, &report);
        let e = &hybrid.instructions[0];
        assert_eq!(e.target_clicks_from_zero, -4);
        assert_eq!(e.hold_mil, -7.0, "hold_mil = {}", e.hold_mil);
        assert!(
            !hybrid.feasible,
            "a -7.0 mil hold exceeds the correct up_mil=5.0 bound and must be infeasible, \
             even though it would fit the WRONG down_mil=10.0 bound an inverted mapping \
             would have checked it against: {hybrid:?}"
        );
    }

    // I1: the original 14 tests all used `windage_mil: 0.0`, so `instructions[1]` was never
    // read -- a bug that transposed elevation's click/CF into the windage arm (or dropped
    // the windage term from the RSS) would have passed every one of them. Deliberately
    // different per-axis click gradutions and CFs so a transposition changes the numbers.
    //
    // Hand-derived (Python-verified against IEEE-754 binary64 arithmetic):
    //   Elevation: 2.34 true mil / cf 0.98 -> corr_dial = 2.387755102040816
    //              / 0.1 mil click -> round to 24 clicks
    //              dial_true = 24*0.1*0.98 = 2.3520000000000003
    //              residual  = 2.34 - 2.352 = -0.012000000000000455 (elevation_cf != 1.0
    //              here, unlike most other tests, so this residual differs from the
    //              cf=1.0, 23-click result used elsewhere in this file)
    //   Windage:  -1.3 true mil / cf 1.05 -> corr_dial = -1.2380952380952381
    //              0.25 MOA click -> click_mil = 0.25*1000/3438 = 0.07271669575334497
    //              corr_dial / click_mil = -17.026... -> round to -17 clicks
    //              dial_true = -17*0.07271669575334497*1.05 = -1.2979930191972078
    //              residual  = -1.3 - (-1.2979930191972078) = -0.0020069808027922686
    //   residual_linear_at_range_m @ 400 m:
    //              e_lin = 0.012000000000000455/1000*400 = 0.0048000000000000004
    //              w_lin = 0.0020069808027922686/1000*400 = 0.0008027923211169075
    //              rss   = sqrt(e_lin^2 + w_lin^2) = 0.004866669858419206
    #[test]
    fn windage_arm_uses_its_own_click_cf_and_contributes_to_the_rss() {
        let mut optic = baseline_profile();
        optic.windage_click = ClickValue { size: 0.25, base: ClickBase::Moa };
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: -1.3 };
        let report = plan_corrections(corr, &optic, 400.0, 0.98, 1.05, &Preferences::default())
            .unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        let e = &dial_all.instructions[0];
        assert_eq!(e.axis, Axis::Elevation);
        assert_eq!(e.target_clicks_from_zero, 24);
        assert_eq!(e.direction, Direction::Up);

        let w = &dial_all.instructions[1];
        assert_eq!(w.axis, Axis::Windage);
        assert_eq!(w.target_clicks_from_zero, -17);
        assert!(
            (w.dial_mil_true - (-1.2979930191972078)).abs() < 1e-12,
            "windage dial_mil_true = {}",
            w.dial_mil_true
        );
        assert_eq!(w.hold_mil, 0.0, "DialAll never holds");
        assert_eq!(w.direction, Direction::Left, "-17 clicks must read as Left, not Up/Down");

        assert!(
            (dial_all.residual_linear_at_range_m - 0.004866669858419206).abs() < 1e-9,
            "residual_linear_at_range_m = {} -- both axes' residuals are nonzero here, so \
             this must be a real two-term RSS, not a single-axis pass-through",
            dial_all.residual_linear_at_range_m
        );

        // Hybrid's windage hold is nonzero and independently exercises the SAME per-axis
        // parameters through `build_axis_instruction`'s separate Hybrid arm.
        let hybrid = plan_for(Strategy::Hybrid, &report);
        let wh = &hybrid.instructions[1];
        assert!(
            (wh.hold_mil - (-0.0020069808027922686)).abs() < 1e-9,
            "hybrid windage hold_mil = {}",
            wh.hold_mil
        );
        assert_eq!(wh.residual_mil, 0.0);
    }

    // I2: HoldAll's hold was never asserted in the original 14 tests, so CF-scaling it
    // (`corr_true / cf` instead of `corr_true`) would have passed all of them -- the exact
    // bug the brief names ("Reticle holds are TRUE angular and are NEVER CF-scaled").
    // Checked at cf = 0.98 specifically: at cf == 1.0 a scaling bug is invisible (dividing
    // by 1.0 is a no-op).
    #[test]
    fn hold_all_hold_is_never_cf_scaled() {
        let optic = baseline_profile();
        let corr = AngularCorrection { elevation_mil: 5.0, windage_mil: 0.0 };
        let report = plan_corrections(corr, &optic, 100.0, 0.98, 1.0, &Preferences::default())
            .unwrap();
        let hold_all = plan_for(Strategy::HoldAll, &report);
        assert_eq!(
            hold_all.instructions[0].hold_mil, 5.0,
            "must equal corr_true exactly (bit-exact copy, no arithmetic at all) -- NOT \
             corr_true / cf = {}",
            5.0_f64 / 0.98
        );
    }

    // I4: `HoldAll`/`Hybrid` carry `residual_mil == 0.0` (and therefore
    // `residual_linear_at_range_m == 0.0`) even when INFEASIBLE -- residual is about exact
    // reconstruction, not achievability. Without `feasible` as the PRIMARY ranking key, an
    // infeasible zero-residual plan could rank ahead of the only plan that actually works,
    // and Task 6's CLI surfaces `plans[0]` as THE recommendation.
    #[test]
    fn infeasible_plans_never_outrank_a_feasible_one() {
        let optic = baseline_profile();
        // 2.34 (not an exact click multiple) so Hybrid's hold is nonzero, and DialAll's
        // residual is nonzero -- both meaningfully exercised, not vacuously all-zero.
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: 0.0 };
        // max_hold_mil: Some(0.0) makes ANY nonzero hold infeasible on HoldAll and Hybrid.
        // DialAll never touches hold at all, so its feasibility depends only on travel,
        // which comfortably fits (23 clicks * 0.1 mil = 2.3 <= 28.0 mil up travel).
        let prefs = Preferences { prefer_hold: false, max_hold_mil: Some(0.0) };
        let report = plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &prefs).unwrap();

        let dial_all = plan_for(Strategy::DialAll, &report);
        let hold_all = plan_for(Strategy::HoldAll, &report);
        let hybrid = plan_for(Strategy::Hybrid, &report);
        assert!(dial_all.feasible, "{dial_all:?}");
        assert!(!hold_all.feasible, "{hold_all:?}");
        assert!(!hybrid.feasible, "{hybrid:?}");
        assert_eq!(hold_all.residual_linear_at_range_m, 0.0);
        assert_eq!(hybrid.residual_linear_at_range_m, 0.0);
        assert!(dial_all.residual_linear_at_range_m > 0.0);

        assert_eq!(
            report.plans[0].strategy,
            Strategy::DialAll,
            "the only feasible plan must rank first, never an infeasible zero-residual one \
             (a residual-only ranking would put Hybrid or HoldAll here instead): {:?}",
            report.plans
        );
        assert!(report.plans[0].feasible);
    }

    // I5: a zero, negative, or non-finite CF divides `corr_true` into infinity, NaN, or a
    // direction-inverting result -- a hard arithmetic failure, not merely an implausible
    // tracking factor, so it is rejected outright rather than only warned about.
    #[test]
    fn non_positive_or_non_finite_tracking_factor_is_rejected() {
        let optic = baseline_profile();
        let corr = AngularCorrection { elevation_mil: 1.0, windage_mil: 0.0 };
        for bad_cf in [0.0, -1.0, -0.5, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let elevation_result =
                plan_corrections(corr, &optic, 100.0, bad_cf, 1.0, &Preferences::default());
            assert!(
                matches!(
                    elevation_result,
                    Err(OpticError::NonPositiveTrackingFactor { field: "elevation_cf", .. })
                ),
                "cf={bad_cf}: {elevation_result:?}"
            );
            let windage_result =
                plan_corrections(corr, &optic, 100.0, 1.0, bad_cf, &Preferences::default());
            assert!(
                matches!(
                    windage_result,
                    Err(OpticError::NonPositiveTrackingFactor { field: "windage_cf", .. })
                ),
                "cf={bad_cf}: {windage_result:?}"
            );
        }
        // A positive CF far outside `tracking_cf_in_range`'s advisory (0.5, 1.5) band is
        // NOT rejected here -- only <= 0 or non-finite is a hard library-level error; the
        // plausibility band stays an advisory, caller/CLI-level concern (see the module doc
        // and `OpticError::NonPositiveTrackingFactor`'s own doc comment).
        let lenient = plan_corrections(corr, &optic, 100.0, 0.001, 1.0, &Preferences::default());
        assert!(lenient.is_ok(), "{lenient:?}");
    }

    // M3: `LimitViolation::needed_mil` must be DIAL-space (`corr_true / cf`), never
    // TRUE-space (`corr_true`) -- indistinguishable at cf == 1.0 (every other test in this
    // file happens to use cf == 1.0 for its travel-violation checks), so this must run at
    // cf != 1.0 to actually discriminate the two.
    #[test]
    fn travel_violation_needed_mil_is_dial_space_not_true_space_at_nonunit_cf() {
        let mut optic = baseline_profile();
        optic.elevation_travel = Some(TravelLimits { down_mil: 0.1, up_mil: 0.1 });
        let corr = AngularCorrection { elevation_mil: 1.0, windage_mil: 0.0 }; // TRUE mil
        let cf = 0.5;
        let report =
            plan_corrections(corr, &optic, 100.0, cf, 1.0, &Preferences::default()).unwrap();
        let dial_all = plan_for(Strategy::DialAll, &report);
        assert!(!dial_all.feasible, "{dial_all:?}");
        let violation = dial_all
            .limits_hit
            .iter()
            .find(|v| matches!(v.kind, LimitKind::TravelExceeded))
            .unwrap_or_else(|| panic!("{:?}", dial_all.limits_hit));
        let corr_dial = corr.elevation_mil / cf; // 2.0
        assert_eq!(
            violation.needed_mil, corr_dial,
            "needed_mil must be DIAL-space (corr_true/cf = {corr_dial}), not TRUE-space \
             ({})",
            corr.elevation_mil
        );
    }

    // M1: `direction` moved from `&'static str` to a real `Direction` enum specifically so
    // `AxisInstruction`/`DialPlan`/`DialPlanReportV1` could derive `Deserialize` again.
    // Pins that the wire form is UNCHANGED (still lowercase `"up"`/`"down"`/`"left"`/
    // `"right"`) and that a full report now genuinely round-trips through JSON.
    #[test]
    fn direction_wire_form_is_unchanged_and_the_full_report_round_trips() {
        for (direction, expected_json) in [
            (Direction::Up, "\"up\""),
            (Direction::Down, "\"down\""),
            (Direction::Left, "\"left\""),
            (Direction::Right, "\"right\""),
        ] {
            assert_eq!(serde_json::to_string(&direction).unwrap(), expected_json);
        }

        let optic = baseline_profile();
        let corr = AngularCorrection { elevation_mil: 2.34, windage_mil: -1.3 };
        let report = plan_corrections(corr, &optic, 100.0, 1.0, 1.0, &Preferences::default())
            .unwrap();
        let json = serde_json::to_string(&report).unwrap();
        let parsed: DialPlanReportV1 = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed, report);
    }
}
