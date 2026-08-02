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
//!
//! A third check keeps `with_axis` consistent with `read_axis` rather than physics-aware: the
//! three wind axes have no single scalar to write when the resolved wind is
//! [`ResolvedWindV1::Segmented`] (taxonomy.rs Known Limitation (c)). `read_axis` already
//! returns `None` for that combination; `with_axis` now returns [`KernelError::AxisAbsent`]
//! for the same combination instead of silently writing a constant-wind field (`speed_mps`
//! etc.) alongside the still-present `wind.segments`, which `solve_v1`'s `resolve_wind` would
//! otherwise reject downstream as an opaque `$.wind` "segments cannot be combined with
//! constant-wind fields" error.
//!
//! The zero-affecting-axis rezero clear (see `with_axis`'s body) is evaluated from
//! `req.shot.zero_distance_m` **after** the axis has been written, not before: for the
//! `ZeroDistance` axis itself, the value that matters is the NEW distance being written, not
//! whatever `zero_distance_m` was on the original request (which may have been absent).
//! Gating on the pre-write value would let a request originally specified purely by
//! `muzzle_angle_rad` silently keep that stale angle after a `ZeroDistance` write, so
//! `solve_v1` would skip the elevation search entirely (it takes an explicit angle over a
//! zero distance whenever both are present) and the new zero distance would have no effect on
//! elevation at all.

use crate::perturbation::taxonomy::{axis_meta, InputAxis};
use crate::solve_json::{
    DragModelV1, PressureReferenceV1, ResolvedSolveRequestV1, ResolvedWindV1, SolveErrorCodeV1,
    SolveRequestV1, TwistDirectionV1, WindReferenceV1,
};
use crate::trajectory_observation::TrajectoryObservationError;

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
    /// Reserved for a later (differentiation) task: an axis whose `axis_meta(axis).kind` is
    /// `AxisKind::Categorical` was asked to be treated as continuous. Not constructed by
    /// `read_axis`/`with_axis` themselves, which accept categorical axes just like any other.
    CategoricalAxis(InputAxis),
    /// The axis is structurally unavailable on THIS request's resolved shape -- currently only
    /// the three wind axes when the resolved wind is `ResolvedWindV1::Segmented`, mirroring
    /// `read_axis` returning `None` for the identical condition (see the module doc).
    AxisAbsent(InputAxis),
    /// The `AxisValue` variant supplied to `with_axis` does not match what `axis` expects
    /// (e.g. a `Scalar` for `DragModel`, or a `Flag` for any continuous axis).
    TypeMismatch(InputAxis),
    /// The axis is well-formed and present, but this particular request's OTHER inputs make
    /// perturbing it physically wrong rather than merely unrepresentable -- see the guards
    /// documented on `with_axis`.
    AxisUnsupportedForRequest {
        axis: InputAxis,
        reason: &'static str,
    },
    /// The rebuilt request failed to resolve or solve. Not constructed here -- `with_axis` only
    /// builds the request, it does not solve it; [`crate::perturbation::evaluate`] (0.33.0
    /// decision-support Task 6) is the constructor, uniformly for every stage that can fail
    /// (`solve_v1::prepare_request`, `solve_v1::build_zeroed_solver`, and
    /// `TrajectorySolver::solve`).
    ///
    /// `code` is carried alongside `message` (review fix I4(a), `derive.rs`) specifically so a
    /// caller -- most notably `central_difference`'s one-sided fallback -- can tell a domain
    /// rejection (`SolveErrorCodeV1::InvalidValue`, produced by `require_range`/
    /// `require_non_negative`/`require_positive` in `solve_v1.rs`) apart from a genuine solver
    /// failure (`SolveFailed`, `ResourceLimit`, `InternalError`) or a malformed-request bug
    /// (`ConflictingFields`, `MissingField`, `UnknownField`, ...) that merely happened to be
    /// triggered by one particular perturbed value. Collapsing this to a `String` (the previous
    /// shape) made that distinction unrecoverable except by parsing the message.
    Solve {
        code: SolveErrorCodeV1,
        message: String,
    },
    /// Post-solve observation extraction failed -- most notably a requested range outside the
    /// computed trajectory. Not constructed here for the same reason as `Solve`;
    /// [`crate::perturbation::evaluate`] constructs it directly from a
    /// [`crate::trajectory_observation::TrajectoryObservationError`], preserved WHOLE (not
    /// collapsed to its `Display` string) for the same reason as `Solve`'s `code` above: only
    /// [`TrajectoryObservationError::OutOfRange`] is a domain rejection eligible for
    /// `central_difference`'s one-sided fallback -- `NonMonotonicTrajectory`, `NonFiniteState`,
    /// `SampleLimitExceeded`, `AllocationFailed`, and the rest are genuine bugs that must
    /// propagate, not be silently reinterpreted as "this side left the domain."
    Observation(TrajectoryObservationError),
    /// A `Scalar` value supplied to `with_axis` was not finite (NaN or infinite).
    NonFinite(InputAxis),
    /// Both perturbed sides of a central difference failed to produce a usable observation --
    /// neither `x + attempted` nor `x - attempted` evaluates -- so not even a one-sided
    /// difference can be built (`central_difference`'s domain fallback,
    /// `src/perturbation/derive.rs`). Distinct from a derivative of zero: this means
    /// "undifferentiable with this step," never "no effect."
    StepOutOfDomain { axis: InputAxis, attempted: f64 },
    /// The same axis was declared more than once in a caller-supplied source list. Not
    /// constructed by `read_axis`/`with_axis`/`central_difference` (none of which see more than
    /// one axis at a time) -- `crate::error_budget::error_budget` (0.33.0 decision-support Task
    /// 10, MBA-1347 review) constructs this from its own up-front validation, because two
    /// entries for the same axis would double-count that axis's variance and corrupt its
    /// leave-one-out counterfactual (removing "the" declaration for that axis is ambiguous when
    /// there are two). `KernelError` had not shipped on any released version when this variant
    /// was added, so there is no compatibility concern in extending it.
    DuplicateAxis(InputAxis),
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
            KernelError::Solve { code, message } => write!(f, "solve failed ({code:?}): {message}"),
            KernelError::Observation(e) => write!(f, "observation failed: {e}"),
            KernelError::NonFinite(a) => write!(f, "axis {a:?} produced a non-finite result"),
            KernelError::StepOutOfDomain { axis, attempted } => write!(
                f,
                "axis {axis:?} could not be differentiated with step {attempted}: both the \
                 forward and backward perturbed values failed to evaluate"
            ),
            KernelError::DuplicateAxis(a) => {
                write!(f, "axis {a:?} was declared more than once")
            }
        }
    }
}
impl std::error::Error for KernelError {}

impl KernelError {
    /// True when this error means the specific perturbed VALUE fell outside the axis's physical
    /// domain -- an `InvalidValue` solve rejection (from `require_range`/`require_non_negative`/
    /// `require_positive` in `solve_v1.rs`), or an `OutOfRange` observation query (a requested
    /// range that fell outside a trajectory shrunk by the perturbation) -- as opposed to a
    /// genuine solver or trajectory failure that merely happened to be triggered by one
    /// particular perturbed value.
    ///
    /// This is the ONLY thing [`crate::perturbation::central_difference`]'s one-sided fallback
    /// (`derive.rs`) may gate on. Review fix I4(a): an earlier revision gated on `Err(_)`
    /// (any failure at all on one side), which silently reinterpreted genuine bugs -- a
    /// non-convergent zero search, a non-finite trajectory state, a sample-limit overrun -- as
    /// if they were domain boundaries, answering with a fabricated derivative instead of
    /// reporting the real failure.
    pub fn is_domain_rejection(&self) -> bool {
        matches!(
            self,
            KernelError::Solve { code: SolveErrorCodeV1::InvalidValue, .. }
                | KernelError::Observation(TrajectoryObservationError::OutOfRange { .. })
        )
    }
}

/// The request's wind-reference echo, whichever wind shape it resolved to.
///
/// `pub(crate)` (0.33.0 decision-support Task 9, review round 4): `explain.rs`'s
/// `is_compass_referenced` needs the identical lookup to detect the `WindDirection` derived-
/// value hazard under compass wind (a different axis than the `ShotAzimuth` guard just below,
/// which is what this function was written for), so it shares this one instead of duplicating
/// the two-arm match a second time.
pub(crate) fn wind_reference_of(w: &ResolvedWindV1) -> Option<WindReferenceV1> {
    match w {
        ResolvedWindV1::Constant(c) => c.wind_reference,
        ResolvedWindV1::Segmented(s) => s.wind_reference,
    }
}

/// Read the current value of one taxonomy axis off a resolved request.
///
/// Returns `None` in two DIFFERENT situations a caller should not conflate:
/// - the axis is an optional request field that was never supplied (`Length`, `Latitude`,
///   `ZeroPoiUp`/`ZeroPoiRight`, `SightOffsetLateral`) -- the axis exists, it just has no
///   value on this particular request; or
/// - the axis is structurally unavailable given how a sibling field resolved -- currently only
///   `WindSpeed`/`WindDirection`/`WindVertical` under `ResolvedWindV1::Segmented`, where there
///   is no single scalar wind value to read at all (taxonomy.rs Known Limitation (c)).
///
/// Either way, `None` means "there is nothing to perturb here": a caller sweeping
/// `InputAxis::ALL` should skip the axis, never invent a value for it.
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

/// Rebuild `r` as a solvable [`SolveRequestV1`] with exactly one axis overwritten to `v`.
///
/// Every other input is carried across unchanged via the reverse conversion
/// (`request_roundtrip.rs`). Writing a `requires_rezero` axis (see
/// [`crate::perturbation::axis_meta`]) while a `zero_distance_m` is present on the REBUILT
/// request clears the carried effective `muzzle_angle_rad`, so the next solve re-zeroes at
/// the (possibly just-changed) distance instead of reusing a stale angle -- see the module
/// doc for why this is checked after the axis is written, not before.
///
/// # Errors
///
/// - [`KernelError::AxisUnsupportedForRequest`] if `axis` is well-formed and present, but this
///   request's OTHER inputs make perturbing it physically wrong (the two guards in the module
///   doc: `Altitude` under QNH pressure, `ShotAzimuth` under compass wind).
/// - [`KernelError::AxisAbsent`] if `axis` is structurally unavailable on this request (the
///   three wind axes under segmented wind), mirroring `read_axis` returning `None` for the
///   same condition.
/// - [`KernelError::TypeMismatch`] if `v`'s kind does not match `axis` (e.g. a `Scalar` for
///   `DragModel`).
/// - [`KernelError::NonFinite`] if `v` is a non-finite `Scalar`.
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
    // Third check (see the module doc): keeps with_axis consistent with read_axis, which
    // already returns None for these three axes under segmented wind. Without this, writing
    // e.g. WindSpeed here would set req.wind.speed_mps alongside the still-present
    // req.wind.segments, which solve_v1's resolve_wind rejects downstream as an opaque
    // "segments cannot be combined with constant-wind fields" error instead of this specific,
    // named one.
    if matches!(
        axis,
        InputAxis::WindSpeed | InputAxis::WindDirection | InputAxis::WindVertical
    ) && matches!(r.wind, ResolvedWindV1::Segmented(_))
    {
        return Err(KernelError::AxisAbsent(axis));
    }

    let mut req: SolveRequestV1 = r.into();
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
    // Changing a zero-affecting axis invalidates the stored effective angle: drop it so the
    // service re-zeroes from zero_distance_m rather than reusing a stale angle. Checked AFTER
    // the match (not before, which was a bug -- see the module doc): for the ZeroDistance axis
    // itself, this must react to the NEW distance just written, not whatever zero_distance_m
    // was on the original request.
    if axis_meta(axis).requires_rezero && req.shot.zero_distance_m.is_some() {
        req.shot.muzzle_angle_rad = None;
    }
    Ok(req)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::perturbation::InputAxis;

    /// Non-default values for every axis this file's tests touch, INCLUDING the five
    /// previously-omitted optional fields (`length_m`, `zero_poi_up_m`, `zero_poi_right_m`,
    /// `sight_offset_lateral_m`, `latitude_rad`) so that all 32 taxonomy axes -- not just 27
    /// of them -- read back `Some(..)` here. Without these five, `read_axis` returns `None`
    /// for `Length`/`ZeroPoiUp`/`ZeroPoiRight`/`SightOffsetLateral`/`Latitude` on this fixture
    /// and a naive per-axis loop (see `every_present_axis_round_trips_without_disturbing_any_other_field`
    /// below) would silently skip exactly the pair -- `ZeroPoiUp`/`ZeroPoiRight` -- most likely
    /// to be transposed by a future edit.
    fn resolved() -> crate::solve_json::ResolvedSolveRequestV1 {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243, "length_m": 0.032},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                      "muzzle_height_m": 0.02, "twist_rate_m_per_turn": 0.2794,
                      "twist_direction": "left", "sight_offset_lateral_m": 0.03},
            "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0, "target_height_m": 0.3,
                     "zero_poi_up_m": 0.01, "zero_poi_right_m": 0.02},
            "atmosphere": {"latitude_rad": 0.6},
            "wind": {"speed_mps": 3.0, "direction_from_rad": std::f64::consts::FRAC_PI_2},
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

    /// ...and the two enum-valued axes are not interchangeable with EACH OTHER either: both
    /// are "enum-shaped" AxisValue variants now, so it's worth confirming DragModel's closure
    /// rejects a TwistDirection value (and not just a Scalar/Flag) as a TypeMismatch.
    #[test]
    fn writing_a_twist_direction_into_the_drag_model_axis_is_a_type_error() {
        let r = resolved();
        let e = with_axis(
            &r,
            InputAxis::DragModel,
            AxisValue::TwistDirection(TwistDirectionV1::Left),
        );
        assert!(matches!(
            e,
            Err(KernelError::TypeMismatch(InputAxis::DragModel))
        ));
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
    /// perturbation. with_axis must AGREE with read_axis about this (I4): it must not
    /// silently write a constant-wind field alongside the still-present wind.segments, which
    /// solve_v1 would later reject downstream as an opaque, differently-worded error.
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

        for axis in [
            InputAxis::WindSpeed,
            InputAxis::WindDirection,
            InputAxis::WindVertical,
        ] {
            let e = with_axis(&r, axis, AxisValue::Scalar(1.0));
            assert!(
                matches!(e, Err(KernelError::AxisAbsent(a)) if a == axis),
                "{axis:?}: expected AxisAbsent, got {e:?}"
            );
        }
    }

    /// Regression test (I1): the rezero-clear used to be gated on the PRE-write
    /// `zero_distance_m`, evaluated before the match instead of after. For a request
    /// originally specified purely by `muzzle_angle_rad` (no `zero_distance_m` at all),
    /// writing `ZeroDistance` would then see the ORIGINAL (absent) zero distance, decide there
    /// was nothing to gate on, and leave the carried angle in place. The rebuilt request would
    /// carry both a brand new `zero_distance_m` AND the stale `muzzle_angle_rad`, and
    /// `solve_v1` always prefers an explicit angle over a zero distance when both are present
    /// (see the module doc), so the elevation search would never run: the new zero distance
    /// would have zero effect on elevation, and a sensitivity sweep over `ZeroDistance` would
    /// silently report "no effect" instead of the truth.
    #[test]
    fn writing_zero_distance_onto_an_angle_only_request_clears_the_carried_angle() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "muzzle_angle_rad": 0.01},
            "atmosphere": {}, "wind": {},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        // Sanity check on the fixture: this request has NO zero_distance_m at all, only a
        // directly-supplied angle, so the bug's precondition actually holds.
        assert_eq!(r.shot.zero_distance_m, None);
        assert_eq!(r.shot.muzzle_angle_rad, 0.01);

        let rebuilt = with_axis(&r, InputAxis::ZeroDistance, AxisValue::Scalar(100.0)).unwrap();
        assert_eq!(rebuilt.shot.zero_distance_m, Some(100.0));
        assert_eq!(
            rebuilt.shot.muzzle_angle_rad, None,
            "the carried angle must be cleared so solve_v1 actually re-zeroes at the new \
             distance instead of skipping the elevation search because an explicit angle was \
             still present"
        );
    }

    /// I2, assertion 1 of 3 ("clears"): a `requires_rezero` axis clears the carried
    /// `muzzle_angle_rad` when a `zero_distance_m` is present.
    #[test]
    fn rezero_axis_clears_the_carried_muzzle_angle_when_a_zero_distance_is_present() {
        let r = resolved(); // zero_distance_m = Some(100.0) in this fixture.
        assert!(r.shot.zero_distance_m.is_some());
        assert!(
            axis_meta(InputAxis::Mass).requires_rezero,
            "fixture assumption: Mass must be a requires_rezero axis for this test to mean \
             anything"
        );

        let changed = with_axis(&r, InputAxis::Mass, AxisValue::Scalar(0.02)).unwrap();
        assert_eq!(changed.shot.muzzle_angle_rad, None);
    }

    /// I2, assertion 2 of 3 ("preserves"): a non-`requires_rezero` axis leaves the carried
    /// `muzzle_angle_rad` untouched, even when a `zero_distance_m` is present. (Perturbing
    /// `MuzzleVelocityMps` in `writing_an_axis_changes_only_that_axis` does NOT demonstrate
    /// this -- that axis IS a rezero axis, so its angle is expected to change; that test never
    /// asserts on `muzzle_angle_rad` at all.)
    #[test]
    fn non_rezero_axis_preserves_the_carried_muzzle_angle() {
        let r = resolved(); // zero_distance_m = Some(100.0) in this fixture.
        assert!(r.shot.zero_distance_m.is_some());
        assert!(
            !axis_meta(InputAxis::WindSpeed).requires_rezero,
            "fixture assumption: WindSpeed must NOT be a requires_rezero axis for this test to \
             mean anything"
        );

        let changed = with_axis(&r, InputAxis::WindSpeed, AxisValue::Scalar(5.0)).unwrap();
        assert_eq!(changed.shot.muzzle_angle_rad, Some(r.shot.muzzle_angle_rad));
    }

    /// I2, assertion 3 of 3 ("the gate"): a `requires_rezero` axis does NOT clear the carried
    /// angle when there is no `zero_distance_m` to begin with -- the clearing is gated on
    /// `zero_distance_m` being present, not on `requires_rezero` alone.
    #[test]
    fn rezero_axis_does_not_clear_the_angle_when_no_zero_distance_is_present() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "muzzle_angle_rad": 0.01},
            "atmosphere": {}, "wind": {},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(r.shot.zero_distance_m, None);
        assert!(axis_meta(InputAxis::Mass).requires_rezero);

        let changed = with_axis(&r, InputAxis::Mass, AxisValue::Scalar(0.02)).unwrap();
        assert_eq!(changed.shot.muzzle_angle_rad, Some(r.shot.muzzle_angle_rad));
    }

    /// Guards the "changes only that axis" property: for each axis present on the fixture,
    /// read its current value and write that SAME value back. Since nothing numerically
    /// changes, the rebuilt request must be identical to `SolveRequestV1::from(&r)` in EVERY
    /// field except `shot.muzzle_angle_rad`, which a `requires_rezero` axis deliberately
    /// clears whenever `zero_distance_m` is present (see the module doc).
    ///
    /// CORRECTNESS CAVEAT (found by review, do not re-widen this claim): this loop is
    /// structurally BLIND to a *self-consistent* transposition, where `read_axis` and
    /// `with_axis` are both wrong about the same axis in the same way -- e.g. both mapping
    /// `ZeroPoiUp` to `zero_poi_right_m`. In that case `read_axis(ZeroPoiUp)` returns
    /// `zero_poi_right_m`'s current value, `with_axis` writes that SAME value back onto that
    /// SAME (wrong) field, and the result is a byte-identical no-op: there is nothing for a
    /// whole-struct diff to see. The same blindness covers a two-axis alias (e.g. `AimAzimuth`
    /// reading AND writing `shot_azimuth_rad`). This loop therefore does NOT prove the mapping
    /// is correct, only that whichever field each axis actually touches, it touches only that
    /// one. `every_axis_writes_to_its_own_named_destination_field` below is what actually
    /// pins each axis to its correct field, via a hardcoded per-axis destination that is
    /// independent of (and so cannot silently agree with) a bug in `read_axis`.
    #[test]
    fn every_present_axis_round_trips_without_disturbing_any_other_field() {
        let r = resolved();
        let baseline: crate::solve_json::SolveRequestV1 = (&r).into();
        // In THIS fixture zero_distance_m is present, so "will the angle be cleared" reduces
        // to "is this a requires_rezero axis" -- computed once, outside the loop, from the
        // same (pre-write) state with_axis itself will see for every axis except
        // ZeroDistance, whose own write cannot change whether zero_distance_m.is_some() (it
        // was already Some, and it is written back Some).
        assert!(r.shot.zero_distance_m.is_some());

        let mut exercised = 0usize;
        for &axis in InputAxis::ALL {
            let Some(v) = read_axis(&r, axis) else {
                continue;
            };
            exercised += 1;
            let rebuilt = with_axis(&r, axis, v).unwrap();

            let mut expected = baseline.clone();
            if axis_meta(axis).requires_rezero {
                expected.shot.muzzle_angle_rad = None;
            }
            assert_eq!(
                rebuilt, expected,
                "writing {axis:?} back onto its own current value changed a field other than \
                 itself (and, for a rezero axis, the carried angle)"
            );
        }
        // Prerequisite the fixture must satisfy for the loop above to mean anything: every one
        // of the 32 taxonomy axes must actually be present (Some) on this fixture, or the loop
        // would silently skip whichever axes read back None -- exactly how a naive version of
        // this fixture would have skipped the ZeroPoiUp/ZeroPoiRight pair before it was
        // enriched with the five previously-omitted optional fields (see `resolved`'s doc).
        assert_eq!(
            exercised,
            InputAxis::ALL.len(),
            "fixture does not make every axis readable -- this test is silently under-covering"
        );
    }

    /// The actual mapping check (I3, corrected): for every one of the 32 axes, write a
    /// SENTINEL value distinct from whatever that axis's own field currently holds, then
    /// assert -- via a hardcoded, explicitly-named `match axis { .. }`, the same pattern
    /// `twist_rate_and_muzzle_height_and_target_height_round_trip` already uses for three
    /// axes -- that the sentinel landed in that axis's OWN specific destination field.
    ///
    /// This is a second, independent statement of the axis -> field mapping, written directly
    /// against the real `SolveRequestV1` field paths rather than derived from `read_axis`. That
    /// independence is the whole point: `every_present_axis_round_trips_without_disturbing_any_other_field`
    /// above reads a value FROM `read_axis` and writes it back, so if `read_axis` and
    /// `with_axis` are both wrong about the same axis in the same way, that loop's read and
    /// write agree with each other and the bug is invisible. This test never calls `read_axis`
    /// at all: the expected destination for every axis is spelled out here by hand, so it
    /// cannot silently agree with a mistake made in `access.rs`'s implementation.
    ///
    /// Distinctness matters: a sentinel MUST differ from the field's pre-write value, or a
    /// transposition that leaves the correct field untouched (because it wrote the sentinel to
    /// the WRONG field instead) would go unnoticed -- the untouched field would coincidentally
    /// already equal the "expected" value. `MuzzleAngle`/`Temperature`/`Pressure` add a fixed
    /// delta to whatever the fixture's resolved value happens to be (avoiding a dependency on
    /// the exact auto-zeroed angle or ICAO default) instead of a hardcoded literal.
    #[test]
    fn every_axis_writes_to_its_own_named_destination_field() {
        let r = resolved();
        use InputAxis::*;

        for &axis in InputAxis::ALL {
            match axis {
                Mass => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.0199)).unwrap();
                    assert_eq!(rebuilt.projectile.mass_kg, 0.0199, "{axis:?}");
                }
                Diameter => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.0090)).unwrap();
                    assert_eq!(rebuilt.projectile.diameter_m, 0.0090, "{axis:?}");
                }
                Length => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.040)).unwrap();
                    assert_eq!(rebuilt.projectile.length_m, Some(0.040), "{axis:?}");
                }
                BallisticCoefficient => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.300)).unwrap();
                    assert_eq!(rebuilt.projectile.ballistic_coefficient, 0.300, "{axis:?}");
                }
                TwistRate => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.3556)).unwrap();
                    assert_eq!(rebuilt.rifle.twist_rate_m_per_turn, Some(0.3556), "{axis:?}");
                }
                TwistDirection => {
                    let rebuilt = with_axis(
                        &r,
                        axis,
                        AxisValue::TwistDirection(TwistDirectionV1::Right),
                    )
                    .unwrap();
                    assert_eq!(
                        rebuilt.rifle.twist_direction,
                        Some(TwistDirectionV1::Right),
                        "{axis:?}"
                    );
                }
                DragModel => {
                    let rebuilt =
                        with_axis(&r, axis, AxisValue::DragModel(DragModelV1::G1)).unwrap();
                    assert_eq!(rebuilt.projectile.drag_model, DragModelV1::G1, "{axis:?}");
                }
                MuzzleVelocityMps => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(900.0)).unwrap();
                    assert_eq!(rebuilt.rifle.muzzle_velocity_mps, 900.0, "{axis:?}");
                }
                SightHeight => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.06)).unwrap();
                    assert_eq!(rebuilt.rifle.sight_height_m, Some(0.06), "{axis:?}");
                }
                ZeroDistance => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(200.0)).unwrap();
                    assert_eq!(rebuilt.shot.zero_distance_m, Some(200.0), "{axis:?}");
                }
                ZeroPoiUp => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.05)).unwrap();
                    assert_eq!(rebuilt.shot.zero_poi_up_m, Some(0.05), "{axis:?}");
                }
                ZeroPoiRight => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.06)).unwrap();
                    assert_eq!(rebuilt.shot.zero_poi_right_m, Some(0.06), "{axis:?}");
                }
                SightOffsetLateral => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.10)).unwrap();
                    assert_eq!(rebuilt.rifle.sight_offset_lateral_m, Some(0.10), "{axis:?}");
                }
                MuzzleHeight => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.05)).unwrap();
                    assert_eq!(rebuilt.rifle.muzzle_height_m, Some(0.05), "{axis:?}");
                }
                MuzzleAngle => {
                    let sentinel = r.shot.muzzle_angle_rad + 0.01;
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(sentinel)).unwrap();
                    assert_eq!(rebuilt.shot.muzzle_angle_rad, Some(sentinel), "{axis:?}");
                }
                Altitude => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(1500.0)).unwrap();
                    assert_eq!(rebuilt.atmosphere.altitude_m, Some(1500.0), "{axis:?}");
                }
                Temperature => {
                    let sentinel = r.atmosphere.temperature_k + 5.0;
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(sentinel)).unwrap();
                    assert_eq!(rebuilt.atmosphere.temperature_k, Some(sentinel), "{axis:?}");
                }
                Pressure => {
                    let sentinel = r.atmosphere.pressure_pa + 500.0;
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(sentinel)).unwrap();
                    assert_eq!(rebuilt.atmosphere.pressure_pa, Some(sentinel), "{axis:?}");
                }
                RelativeHumidity => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.7)).unwrap();
                    assert_eq!(rebuilt.atmosphere.relative_humidity, Some(0.7), "{axis:?}");
                }
                Latitude => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.3)).unwrap();
                    assert_eq!(rebuilt.atmosphere.latitude_rad, Some(0.3), "{axis:?}");
                }
                WindSpeed => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(7.0)).unwrap();
                    assert_eq!(rebuilt.wind.speed_mps, Some(7.0), "{axis:?}");
                }
                WindDirection => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.4)).unwrap();
                    assert_eq!(rebuilt.wind.direction_from_rad, Some(0.4), "{axis:?}");
                }
                WindVertical => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(1.5)).unwrap();
                    assert_eq!(rebuilt.wind.vertical_speed_mps, Some(1.5), "{axis:?}");
                }
                TargetDistance => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(1000.0)).unwrap();
                    assert_eq!(rebuilt.shot.max_range_m, 1000.0, "{axis:?}");
                }
                ShootingAngle => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.2)).unwrap();
                    assert_eq!(rebuilt.shot.shooting_angle_rad, Some(0.2), "{axis:?}");
                }
                Cant => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.15)).unwrap();
                    assert_eq!(rebuilt.shot.cant_angle_rad, Some(0.15), "{axis:?}");
                }
                ShotAzimuth => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.25)).unwrap();
                    assert_eq!(rebuilt.shot.shot_azimuth_rad, Some(0.25), "{axis:?}");
                }
                AimAzimuth => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.35)).unwrap();
                    assert_eq!(rebuilt.shot.aim_azimuth_rad, Some(0.35), "{axis:?}");
                }
                TargetHeight => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Scalar(0.5)).unwrap();
                    assert_eq!(rebuilt.shot.target_height_m, Some(0.5), "{axis:?}");
                }
                MagnusEnabled => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Flag(true)).unwrap();
                    assert_eq!(rebuilt.effects.magnus, Some(true), "{axis:?}");
                }
                CoriolisEnabled => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Flag(true)).unwrap();
                    assert_eq!(rebuilt.effects.coriolis, Some(true), "{axis:?}");
                }
                EnhancedSpinDriftEnabled => {
                    let rebuilt = with_axis(&r, axis, AxisValue::Flag(true)).unwrap();
                    assert_eq!(rebuilt.effects.enhanced_spin_drift, Some(true), "{axis:?}");
                }
            }
        }
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

    /// The QNH guard must be specific to QNH mode: an OMITTED pressure_reference (the default,
    /// and every request that predates MBA-1397) must still allow perturbing altitude
    /// normally. This only proves the guard doesn't fire on absence; see
    /// `altitude_axis_is_perturbable_when_pressure_is_explicitly_absolute` below for the
    /// stronger claim that it doesn't fire on an explicit non-QNH value either.
    #[test]
    fn altitude_axis_is_perturbable_when_pressure_is_absolute() {
        let r = resolved();
        assert_eq!(r.atmosphere.pressure_reference, None);
        let changed = with_axis(&r, InputAxis::Altitude, AxisValue::Scalar(1200.0)).unwrap();
        assert_eq!(changed.atmosphere.altitude_m, Some(1200.0));
    }

    /// Stronger version of the test above: an EXPLICIT `pressure_reference: "absolute"` (not
    /// just an omitted field defaulting to it) must still allow perturbing altitude normally.
    /// Guards against a hypothetical `!= Some(Qnh)` -> `== None` typo that happens to pass the
    /// omitted-field test above but would reject this explicit-but-equivalent case.
    #[test]
    fn altitude_axis_is_perturbable_when_pressure_is_explicitly_absolute() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0},
            "atmosphere": {"altitude_m": 500.0, "temperature_k": 288.0, "pressure_pa": 101325.0,
                           "pressure_reference": "absolute"},
            "wind": {}, "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        assert_eq!(
            r.atmosphere.pressure_reference,
            Some(PressureReferenceV1::Absolute)
        );

        let changed = with_axis(&r, InputAxis::Altitude, AxisValue::Scalar(600.0)).unwrap();
        assert_eq!(changed.atmosphere.altitude_m, Some(600.0));
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

    /// The compass guard must be specific to compass mode: an OMITTED wind_reference (the
    /// default) must still allow perturbing shot azimuth normally. This only proves the guard
    /// doesn't fire on absence; see
    /// `shot_azimuth_axis_is_perturbable_when_wind_is_explicitly_shooter_relative` below for
    /// the stronger claim that it doesn't fire on an explicit non-compass value either.
    #[test]
    fn shot_azimuth_axis_is_perturbable_when_wind_is_shooter_relative() {
        let r = resolved();
        assert_eq!(
            match &r.wind {
                ResolvedWindV1::Constant(c) => c.wind_reference,
                ResolvedWindV1::Segmented(_) => panic!("constant wind expected"),
            },
            None
        );
        let changed = with_axis(&r, InputAxis::ShotAzimuth, AxisValue::Scalar(0.5)).unwrap();
        assert_eq!(changed.shot.shot_azimuth_rad, Some(0.5));
    }

    /// Stronger version of the test above: an EXPLICIT `wind_reference: "shooter"` (not just
    /// an omitted field defaulting to it) must still allow perturbing shot azimuth normally.
    /// Guards against a hypothetical `!= Some(Compass)` -> `== None` typo that happens to pass
    /// the omitted-field test above but would reject this explicit-but-equivalent case.
    #[test]
    fn shot_azimuth_axis_is_perturbable_when_wind_is_explicitly_shooter_relative() {
        let json = serde_json::json!({
            "schema_version": 1,
            "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                           "ballistic_coefficient": 0.243},
            "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05},
            "shot": {"max_range_m": 900.0, "shot_azimuth_rad": 0.3},
            "atmosphere": {},
            "wind": {"speed_mps": 3.0, "direction_from_rad": 1.0, "wind_reference": "shooter"},
            "solver": {}, "effects": {}, "sampling": {"interval_m": 50.0}
        })
        .to_string();
        let req = crate::solve_json::decode_solve_request_v1(&json).unwrap();
        let r = crate::solve_v1::solve_v1(req).unwrap().resolved_request;
        let ResolvedWindV1::Constant(wind) = &r.wind else {
            panic!("constant wind expected");
        };
        assert_eq!(wind.wind_reference, Some(WindReferenceV1::Shooter));

        let changed = with_axis(&r, InputAxis::ShotAzimuth, AxisValue::Scalar(0.5)).unwrap();
        assert_eq!(changed.shot.shot_azimuth_rad, Some(0.5));
    }
}
