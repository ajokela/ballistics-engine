//! Reading and writing one taxonomy axis on a canonical request.
//!
//! `with_axis` produces a full `SolveRequestV1` via the Phase 0 reverse conversion
//! (`request_roundtrip.rs`) and then overwrites exactly one field, so every other input is
//! carried across unchanged.
//!
//! Two axes are enum-valued rather than scalar or boolean (`DragModel`: G1/G6/G7/G8;
//! `TwistDirection`: left/right), so [`AxisValue`] carries a dedicated variant for each
//! instead of round-tripping through a string -- a string round-trip would turn a typo into a
//! runtime [`KernelError::TypeMismatch`] instead of a compile error.
//!
//! `with_axis` also refuses two axis/request combinations that would otherwise silently build
//! a physically WRONG counterfactual, because [the reverse conversion](crate::solve_json)
//! (see `request_roundtrip.rs`'s module doc) always emits absolute station pressure and
//! shooter-relative wind regardless of how the original request entered them:
//!
//! - `Altitude` when the original request declared a QNH pressure: perturbing altitude would
//!   change air density-by-altitude without moving the QNH-referenced station pressure the
//!   way the caller means -- the opposite of the intended counterfactual.
//! - `ShotAzimuth` when the original request declared compass-referenced wind: perturbing the
//!   shot azimuth would rotate the wind WITH the rifle instead of keeping it earth-fixed --
//!   physically inverted.
//!
//! Both are detected from the ORIGINAL resolved request's echoed reference mode (added in
//! Task 1 precisely so this is detectable) and rejected with
//! [`KernelError::AxisUnsupportedForRequest`] rather than producing a silently-wrong request.
//! See `taxonomy.rs`'s "KNOWN LIMITATIONS" comment, items (a) and (b).

use crate::perturbation::taxonomy::{axis_meta, InputAxis};
use crate::solve_json::{
    DragModelV1, PressureReferenceV1, ResolvedSolveRequestV1, ResolvedWindV1, SolveRequestV1,
    TwistDirectionV1, WindReferenceV1,
};

/// One taxonomy axis's value, read from or written to a request.
///
/// `Scalar`/`Flag` cover every continuous and boolean axis. `DragModel`/`TwistDirection`
/// cover the two enum-valued axes losslessly -- see the module doc for why a `String` variant
/// would be the wrong choice.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AxisValue {
    Scalar(f64),
    Flag(bool),
    DragModel(DragModelV1),
    TwistDirection(TwistDirectionV1),
}

#[derive(Debug, Clone, PartialEq)]
pub enum KernelError {
    CategoricalAxis(InputAxis),
    AxisAbsent(InputAxis),
    TypeMismatch(InputAxis),
    /// The axis is well-formed and present, but this particular request's OTHER inputs make
    /// perturbing it physically wrong rather than merely unrepresentable -- see the guards
    /// documented on `with_axis`.
    AxisUnsupportedForRequest {
        axis: InputAxis,
        reason: &'static str,
    },
    Solve(String),
    Observation(String),
    NonFinite(InputAxis),
}

impl std::fmt::Display for KernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KernelError::CategoricalAxis(a) => {
                write!(f, "axis {a:?} is categorical and cannot be differentiated")
            }
            KernelError::AxisAbsent(a) => write!(f, "axis {a:?} is not present in this request"),
            KernelError::TypeMismatch(a) => write!(f, "value type does not match axis {a:?}"),
            KernelError::AxisUnsupportedForRequest { axis, reason } => {
                write!(f, "axis {axis:?} is not supported for this request: {reason}")
            }
            KernelError::Solve(m) => write!(f, "solve failed: {m}"),
            KernelError::Observation(m) => write!(f, "observation failed: {m}"),
            KernelError::NonFinite(a) => write!(f, "axis {a:?} produced a non-finite result"),
        }
    }
}
impl std::error::Error for KernelError {}

/// The request's wind-reference echo, whichever wind shape it resolved to.
fn wind_reference_of(w: &ResolvedWindV1) -> Option<WindReferenceV1> {
    match w {
        ResolvedWindV1::Constant(c) => c.wind_reference,
        ResolvedWindV1::Segmented(s) => s.wind_reference,
    }
}

pub fn read_axis(r: &ResolvedSolveRequestV1, axis: InputAxis) -> Option<AxisValue> {
    use InputAxis::*;
    let wind_speed = match &r.wind {
        ResolvedWindV1::Constant(c) => Some(c.speed_mps),
        ResolvedWindV1::Segmented(_) => None,
    };
    let wind_dir = match &r.wind {
        ResolvedWindV1::Constant(c) => Some(c.direction_from_rad),
        ResolvedWindV1::Segmented(_) => None,
    };
    let wind_vert = match &r.wind {
        ResolvedWindV1::Constant(c) => Some(c.vertical_speed_mps),
        ResolvedWindV1::Segmented(_) => None,
    };
    Some(match axis {
        Mass => AxisValue::Scalar(r.projectile.mass_kg),
        Diameter => AxisValue::Scalar(r.projectile.diameter_m),
        Length => AxisValue::Scalar(r.projectile.length_m?),
        BallisticCoefficient => AxisValue::Scalar(r.projectile.ballistic_coefficient),
        TwistRate => AxisValue::Scalar(r.rifle.twist_rate_m_per_turn),
        TwistDirection => AxisValue::TwistDirection(r.rifle.twist_direction),
        DragModel => AxisValue::DragModel(r.projectile.drag_model),
        MuzzleVelocityMps => AxisValue::Scalar(r.rifle.muzzle_velocity_mps),
        SightHeight => AxisValue::Scalar(r.rifle.sight_height_m),
        ZeroDistance => AxisValue::Scalar(r.shot.zero_distance_m?),
        ZeroPoiUp => AxisValue::Scalar(r.shot.zero_poi_up_m?),
        ZeroPoiRight => AxisValue::Scalar(r.shot.zero_poi_right_m?),
        SightOffsetLateral => AxisValue::Scalar(r.rifle.sight_offset_lateral_m?),
        MuzzleHeight => AxisValue::Scalar(r.rifle.muzzle_height_m),
        MuzzleAngle => AxisValue::Scalar(r.shot.muzzle_angle_rad),
        Altitude => AxisValue::Scalar(r.atmosphere.altitude_m),
        Temperature => AxisValue::Scalar(r.atmosphere.temperature_k),
        Pressure => AxisValue::Scalar(r.atmosphere.pressure_pa),
        RelativeHumidity => AxisValue::Scalar(r.atmosphere.relative_humidity),
        Latitude => AxisValue::Scalar(r.atmosphere.latitude_rad?),
        WindSpeed => AxisValue::Scalar(wind_speed?),
        WindDirection => AxisValue::Scalar(wind_dir?),
        WindVertical => AxisValue::Scalar(wind_vert?),
        TargetDistance => AxisValue::Scalar(r.shot.max_range_m),
        ShootingAngle => AxisValue::Scalar(r.shot.shooting_angle_rad),
        Cant => AxisValue::Scalar(r.shot.cant_angle_rad),
        ShotAzimuth => AxisValue::Scalar(r.shot.shot_azimuth_rad),
        AimAzimuth => AxisValue::Scalar(r.shot.aim_azimuth_rad),
        TargetHeight => AxisValue::Scalar(r.shot.target_height_m),
        MagnusEnabled => AxisValue::Flag(r.effects.magnus),
        CoriolisEnabled => AxisValue::Flag(r.effects.coriolis),
        EnhancedSpinDriftEnabled => AxisValue::Flag(r.effects.enhanced_spin_drift),
    })
}

pub fn with_axis(
    r: &ResolvedSolveRequestV1,
    axis: InputAxis,
    v: AxisValue,
) -> Result<SolveRequestV1, KernelError> {
    // Physics guards (see the module doc): detected from the ORIGINAL resolved request's
    // echoed reference mode, before anything else, so a dangerous combination never reaches
    // the reverse conversion below.
    if axis == InputAxis::Altitude
        && r.atmosphere.pressure_reference == Some(PressureReferenceV1::Qnh)
    {
        return Err(KernelError::AxisUnsupportedForRequest {
            axis,
            reason: "the original request declared a QNH pressure_reference; the rebuilt \
                     request always carries absolute station pressure (request_roundtrip.rs \
                     cannot re-derive the original altimeter setting), so perturbing altitude \
                     would change air density-by-altitude without moving the QNH-referenced \
                     station pressure the way the caller means",
        });
    }
    if axis == InputAxis::ShotAzimuth
        && wind_reference_of(&r.wind) == Some(WindReferenceV1::Compass)
    {
        return Err(KernelError::AxisUnsupportedForRequest {
            axis,
            reason: "the original request declared compass-referenced wind; the rebuilt \
                     request always carries shooter-relative wind (request_roundtrip.rs \
                     cannot re-derive the original earth-fixed bearing), so perturbing the \
                     shot azimuth would rotate the wind WITH the rifle instead of keeping it \
                     earth-fixed",
        });
    }

    let mut req: SolveRequestV1 = r.into();
    let meta = axis_meta(axis);
    let scalar = |v: AxisValue| -> Result<f64, KernelError> {
        match v {
            AxisValue::Scalar(x) if x.is_finite() => Ok(x),
            AxisValue::Scalar(_) => Err(KernelError::NonFinite(axis)),
            _ => Err(KernelError::TypeMismatch(axis)),
        }
    };
    let flag = |v: AxisValue| -> Result<bool, KernelError> {
        match v {
            AxisValue::Flag(b) => Ok(b),
            _ => Err(KernelError::TypeMismatch(axis)),
        }
    };
    let drag_model = |v: AxisValue| -> Result<DragModelV1, KernelError> {
        match v {
            AxisValue::DragModel(m) => Ok(m),
            _ => Err(KernelError::TypeMismatch(axis)),
        }
    };
    let twist_direction = |v: AxisValue| -> Result<TwistDirectionV1, KernelError> {
        match v {
            AxisValue::TwistDirection(d) => Ok(d),
            _ => Err(KernelError::TypeMismatch(axis)),
        }
    };
    // Changing a zero-affecting axis invalidates the stored effective angle: drop it so
    // the service re-zeroes from zero_distance_m rather than reusing a stale angle.
    if meta.requires_rezero && req.shot.zero_distance_m.is_some() {
        req.shot.muzzle_angle_rad = None;
    }
    use InputAxis::*;
    match axis {
        Mass => req.projectile.mass_kg = scalar(v)?,
        Diameter => req.projectile.diameter_m = scalar(v)?,
        Length => req.projectile.length_m = Some(scalar(v)?),
        BallisticCoefficient => req.projectile.ballistic_coefficient = scalar(v)?,
        TwistRate => req.rifle.twist_rate_m_per_turn = Some(scalar(v)?),
        TwistDirection => req.rifle.twist_direction = Some(twist_direction(v)?),
        DragModel => req.projectile.drag_model = drag_model(v)?,
        MuzzleVelocityMps => req.rifle.muzzle_velocity_mps = scalar(v)?,
        SightHeight => req.rifle.sight_height_m = Some(scalar(v)?),
        ZeroDistance => req.shot.zero_distance_m = Some(scalar(v)?),
        ZeroPoiUp => req.shot.zero_poi_up_m = Some(scalar(v)?),
        ZeroPoiRight => req.shot.zero_poi_right_m = Some(scalar(v)?),
        SightOffsetLateral => req.rifle.sight_offset_lateral_m = Some(scalar(v)?),
        MuzzleHeight => req.rifle.muzzle_height_m = Some(scalar(v)?),
        MuzzleAngle => req.shot.muzzle_angle_rad = Some(scalar(v)?),
        Altitude => req.atmosphere.altitude_m = Some(scalar(v)?),
        Temperature => req.atmosphere.temperature_k = Some(scalar(v)?),
        Pressure => req.atmosphere.pressure_pa = Some(scalar(v)?),
        RelativeHumidity => req.atmosphere.relative_humidity = Some(scalar(v)?),
        Latitude => req.atmosphere.latitude_rad = Some(scalar(v)?),
        WindSpeed => req.wind.speed_mps = Some(scalar(v)?),
        WindDirection => req.wind.direction_from_rad = Some(scalar(v)?),
        WindVertical => req.wind.vertical_speed_mps = Some(scalar(v)?),
        TargetDistance => req.shot.max_range_m = scalar(v)?,
        ShootingAngle => req.shot.shooting_angle_rad = Some(scalar(v)?),
        Cant => req.shot.cant_angle_rad = Some(scalar(v)?),
        ShotAzimuth => req.shot.shot_azimuth_rad = Some(scalar(v)?),
        AimAzimuth => req.shot.aim_azimuth_rad = Some(scalar(v)?),
        TargetHeight => req.shot.target_height_m = Some(scalar(v)?),
        MagnusEnabled => req.effects.magnus = Some(flag(v)?),
        CoriolisEnabled => req.effects.coriolis = Some(flag(v)?),
        EnhancedSpinDriftEnabled => req.effects.enhanced_spin_drift = Some(flag(v)?),
    }
    Ok(req)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbation::InputAxis;

    /// Non-default values for every field this file's tests touch, including the five axes
    /// added to the taxonomy after this task's brief was written (TwistRate, TwistDirection,
    /// MuzzleHeight, DragModel, TargetHeight), so a mixed-up field mapping would show up as a
    /// wrong value rather than hiding behind a coincidental default.
    fn resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                      "muzzle_height_m": 0.02, "twist_rate_m_per_turn": 0.2794,
                      "twist_direction": "left"},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0, "target_height_m": 0.3},
            "atmosphere": {}, "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        }).to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        crate::solve_v1::solve_v1(req).unwrap().resolved_request
    }

    #[test]
    fn read_then_write_is_identity() {
        let r = resolved();
        let v = read_axis(&r, InputAxis::MuzzleVelocityMps).expect("axis present");
        let rebuilt = with_axis(&r, InputAxis::MuzzleVelocityMps, v).unwrap();
        assert_eq!(rebuilt.rifle.muzzle_velocity_mps, r.rifle.muzzle_velocity_mps);
    }

    #[test]
    fn writing_an_axis_changes_only_that_axis() {
        let r = resolved();
        let changed = with_axis(&r, InputAxis::MuzzleVelocityMps, AxisValue::Scalar(900.0)).unwrap();
        let baseline: crate::solve_json::SolveRequestV1 = (&r).into();
        assert_eq!(changed.rifle.muzzle_velocity_mps, 900.0);
        assert_eq!(changed.atmosphere.pressure_pa, baseline.atmosphere.pressure_pa);
        assert_eq!(changed.wind.speed_mps, baseline.wind.speed_mps);
        assert_eq!(changed.shot.max_range_m, baseline.shot.max_range_m);
        // The five axes added after the brief was written must also be carried unchanged.
        assert_eq!(changed.projectile.drag_model, baseline.projectile.drag_model);
        assert_eq!(changed.rifle.twist_rate_m_per_turn, baseline.rifle.twist_rate_m_per_turn);
        assert_eq!(changed.rifle.twist_direction, baseline.rifle.twist_direction);
        assert_eq!(changed.rifle.muzzle_height_m, baseline.rifle.muzzle_height_m);
        assert_eq!(changed.shot.target_height_m, baseline.shot.target_height_m);
    }

    #[test]
    fn writing_a_flag_into_a_scalar_axis_is_a_type_error() {
        let r = resolved();
        let e = with_axis(&r, InputAxis::MuzzleVelocityMps, AxisValue::Flag(true));
        assert!(matches!(e, Err(KernelError::TypeMismatch(_))));
    }

    #[test]
    fn effect_flags_round_trip() {
        let r = resolved();
        let changed = with_axis(&r, InputAxis::CoriolisEnabled, AxisValue::Flag(true)).unwrap();
        assert_eq!(changed.effects.coriolis, Some(true));
    }

    /// DragModel is a four-way enum (G1/G6/G7/G8), not a bool or float: AxisValue needs a
    /// dedicated variant to carry it losslessly (see the module doc).
    #[test]
    fn drag_model_axis_reads_and_writes_the_enum_value() {
        let r = resolved();
        let v = read_axis(&r, InputAxis::DragModel).expect("axis present");
        assert_eq!(v, AxisValue::DragModel(DragModelV1::G7));

        let changed =
            with_axis(&r, InputAxis::DragModel, AxisValue::DragModel(DragModelV1::G1)).unwrap();
        assert_eq!(changed.projectile.drag_model, DragModelV1::G1);
        // Nothing else moved.
        assert_eq!(changed.projectile.mass_kg, r.projectile.mass_kg);
        assert_eq!(
            changed.projectile.ballistic_coefficient,
            r.projectile.ballistic_coefficient
        );
    }

    /// TwistDirection is a two-way enum (left/right); same rationale as DragModel above.
    #[test]
    fn twist_direction_axis_reads_and_writes_the_enum_value() {
        let r = resolved();
        let v = read_axis(&r, InputAxis::TwistDirection).expect("axis present");
        assert_eq!(v, AxisValue::TwistDirection(TwistDirectionV1::Left));

        let changed = with_axis(
            &r,
            InputAxis::TwistDirection,
            AxisValue::TwistDirection(TwistDirectionV1::Right),
        )
        .unwrap();
        assert_eq!(changed.rifle.twist_direction, Some(TwistDirectionV1::Right));
    }

    /// Extending AxisValue with enum variants must not weaken the existing type-checking:
    /// a scalar is still rejected for an enum-valued axis.
    #[test]
    fn writing_a_scalar_into_the_drag_model_axis_is_a_type_error() {
        let r = resolved();
        let e = with_axis(&r, InputAxis::DragModel, AxisValue::Scalar(1.0));
        assert!(matches!(
            e,
            Err(KernelError::TypeMismatch(InputAxis::DragModel))
        ));
    }

    /// ...and the reverse: an enum value is rejected for a scalar axis.
    #[test]
    fn writing_a_drag_model_into_a_scalar_axis_is_a_type_error() {
        let r = resolved();
        let e = with_axis(&r, InputAxis::Mass, AxisValue::DragModel(DragModelV1::G1));
        assert!(matches!(e, Err(KernelError::TypeMismatch(InputAxis::Mass))));
    }

    /// Read/write identity for three of the five axes added after the brief was written,
    /// exercising each field mapping directly rather than only checking non-interference.
    #[test]
    fn twist_rate_and_muzzle_height_and_target_height_round_trip() {
        let r = resolved();
        for axis in [InputAxis::TwistRate, InputAxis::MuzzleHeight, InputAxis::TargetHeight] {
            let v = read_axis(&r, axis).unwrap_or_else(|| panic!("{axis:?} should be present"));
            let rebuilt = with_axis(&r, axis, v).unwrap();
            match (axis, v) {
                (InputAxis::TwistRate, AxisValue::Scalar(x)) => {
                    assert_eq!(rebuilt.rifle.twist_rate_m_per_turn, Some(x))
                }
                (InputAxis::MuzzleHeight, AxisValue::Scalar(x)) => {
                    assert_eq!(rebuilt.rifle.muzzle_height_m, Some(x))
                }
                (InputAxis::TargetHeight, AxisValue::Scalar(x)) => {
                    assert_eq!(rebuilt.shot.target_height_m, Some(x))
                }
                _ => panic!("expected a scalar for {axis:?}"),
            }
        }
    }

    /// read_axis returns None for the three wind axes under segmented wind: there is no
    /// single scalar to perturb (taxonomy.rs Known Limitation (c)). This is intended -- a
    /// caller must treat the None as "axis absent," never invent a uniform per-segment
    /// perturbation.
    #[test]
    fn wind_axes_are_absent_under_segmented_wind() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {},
            "wind": {"segments": [{"until_distance_m": 900.0, "speed_mps": 3.0,
                                    "direction_from_rad": 1.0}]},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert!(matches!(r.wind, ResolvedWindV1::Segmented(_)));

        assert_eq!(read_axis(&r, InputAxis::WindSpeed), None);
        assert_eq!(read_axis(&r, InputAxis::WindDirection), None);
        assert_eq!(read_axis(&r, InputAxis::WindVertical), None);
    }

    /// Physics guard (a): perturbing altitude on a QNH-pressure request would silently build
    /// a wrong counterfactual (the rebuilt request always carries absolute station pressure --
    /// see request_roundtrip.rs -- so density-by-altitude would move but the QNH-referenced
    /// station pressure would not, backwards from what a QNH-entering caller means). Detected
    /// from the ORIGINAL resolved request's pressure_reference echo and rejected outright.
    #[test]
    fn altitude_axis_is_unsupported_when_the_original_pressure_was_qnh() {
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
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(r.atmosphere.pressure_reference, Some(PressureReferenceV1::Qnh));

        let e = with_axis(&r, InputAxis::Altitude, AxisValue::Scalar(600.0));
        match e {
            Err(KernelError::AxisUnsupportedForRequest { axis: InputAxis::Altitude, reason }) => {
                assert!(
                    reason.to_lowercase().contains("qnh"),
                    "reason should name QNH: {reason}"
                );
            }
            other => panic!("expected AxisUnsupportedForRequest, got {other:?}"),
        }
    }

    /// The QNH guard must be specific to QNH mode: an ordinary absolute-pressure request (the
    /// default, and every request that predates MBA-1397) must still allow perturbing
    /// altitude normally.
    #[test]
    fn altitude_axis_is_perturbable_when_pressure_is_absolute() {
        let r = resolved();
        assert_eq!(r.atmosphere.pressure_reference, None);
        let changed = with_axis(&r, InputAxis::Altitude, AxisValue::Scalar(1200.0)).unwrap();
        assert_eq!(changed.atmosphere.altitude_m, Some(1200.0));
    }

    /// Physics guard (b): perturbing shot azimuth on a compass-referenced-wind request would
    /// silently build a wrong counterfactual (the rebuilt request always carries
    /// shooter-relative wind -- see request_roundtrip.rs -- so the wind would rotate WITH the
    /// rifle instead of staying earth-fixed). Detected from the ORIGINAL resolved request's
    /// wind_reference echo and rejected outright.
    #[test]
    fn shot_azimuth_axis_is_unsupported_when_the_original_wind_was_compass_referenced() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "shot_azimuth_rad": 0.3},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": 1.0, "wind_reference": "compass"},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;

        let e = with_axis(&r, InputAxis::ShotAzimuth, AxisValue::Scalar(0.5));
        match e {
            Err(KernelError::AxisUnsupportedForRequest {
                axis: InputAxis::ShotAzimuth,
                reason,
            }) => {
                assert!(
                    reason.to_lowercase().contains("compass"),
                    "reason should name compass wind: {reason}"
                );
            }
            other => panic!("expected AxisUnsupportedForRequest, got {other:?}"),
        }
    }

    /// The compass guard must be specific to compass mode: ordinary shooter-relative wind
    /// (the default) must still allow perturbing shot azimuth normally.
    #[test]
    fn shot_azimuth_axis_is_perturbable_when_wind_is_shooter_relative() {
        let r = resolved();
        let changed = with_axis(&r, InputAxis::ShotAzimuth, AxisValue::Scalar(0.5)).unwrap();
        assert_eq!(changed.shot.shot_azimuth_rad, Some(0.5));
    }
}
