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

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::adjustment::ClickValue;

/// A shooter's complete turret and reticle geometry (MBA-1348).
///
/// Every angular field anywhere in this type is in MIL. Turret-related fields
/// (`clicks_per_revolution`, `elevation_travel`, `windage_travel`, `turret_state`) are in
/// DIAL-space; `reticle_hold_bounds` is TRUE angular. See the module docs for what that
/// distinction means and why it matters.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct OpticProfile {
    /// Elevation turret's click graduation (its engraved size, e.g. 0.1 mil or 1/4 MOA).
    pub elevation_click: ClickValue,
    /// Windage turret's click graduation.
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
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
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
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
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
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
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
}

impl OpticProfile {
    /// Checks internal consistency: every angular field is finite, every travel/hold-bound
    /// magnitude is non-negative, `clicks_per_revolution` is a sane count, and — when both
    /// are present — `turret_state` lies within `elevation_travel` / `windage_travel`.
    ///
    /// This cannot and does not check that the profile matches any real, physical scope;
    /// it only rejects shapes that are self-contradictory or would corrupt downstream
    /// arithmetic (a zero-length revolution, NaN propagation, a dialed state the turret
    /// could not physically reach given its own declared travel).
    pub fn validate(&self) -> Result<(), OpticError> {
        require_finite("elevation_click.size", self.elevation_click.size)?;
        require_finite("windage_click.size", self.windage_click.size)?;

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
}
