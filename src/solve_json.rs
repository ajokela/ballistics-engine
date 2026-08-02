//! Versioned, binding-neutral JSON data-transfer objects for trajectory solves.
//!
//! These types deliberately do not expose [`crate::BallisticInputs`]. The public JSON contract
//! has its own names, units, defaults, and compatibility policy so engine internals can evolve
//! without silently changing requests produced by language bindings.

use serde::{
    de::{self, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};
use serde_json::{Map, Value};
use std::{fmt, num::NonZeroUsize};

/// The only solve-json schema version understood by this module.
pub const SOLVE_JSON_SCHEMA_VERSION_V1: u32 = 1;

/// Maximum number of trajectory observations in one solve-json v1 success response.
///
/// Service implementations must reject a response above this limit with
/// [`SolveSuccessV1::validate_for_serialization`] instead of truncating it.
pub const MAX_SOLVE_JSON_SAMPLES_V1: usize = 10_000;

/// Deserialize a request member that may be omitted but may not be JSON `null`.
///
/// Serde supplies [`Option::None`] from the field's `default` only when the member is absent.
/// When the member is present, deserialize `T` directly so `null` is rejected for scalar,
/// enum, object, and array values alike.
fn deserialize_present<'de, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(deserializer).map(Some)
}

/// An invariant solve-json v1 schema discriminator.
///
/// This unit type has no invalid public state: it always serializes as the JSON integer `1` and
/// deserializes only from that integer.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct SchemaVersionV1;

impl SchemaVersionV1 {
    /// Return the integer represented on the wire.
    pub const fn get(self) -> u32 {
        SOLVE_JSON_SCHEMA_VERSION_V1
    }
}

impl Serialize for SchemaVersionV1 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_u32(SOLVE_JSON_SCHEMA_VERSION_V1)
    }
}

impl<'de> Deserialize<'de> for SchemaVersionV1 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SchemaVersionVisitor;

        impl<'de> Visitor<'de> for SchemaVersionVisitor {
            type Value = SchemaVersionV1;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("the integer 1")
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                if value == u64::from(SOLVE_JSON_SCHEMA_VERSION_V1) {
                    Ok(SchemaVersionV1)
                } else {
                    Err(E::custom(format!(
                        "unsupported schema_version {value}; expected {SOLVE_JSON_SCHEMA_VERSION_V1}"
                    )))
                }
            }

            fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                if value == i64::from(SOLVE_JSON_SCHEMA_VERSION_V1) {
                    Ok(SchemaVersionV1)
                } else {
                    Err(E::custom(format!(
                        "unsupported schema_version {value}; expected {SOLVE_JSON_SCHEMA_VERSION_V1}"
                    )))
                }
            }
        }

        deserializer.deserialize_any(SchemaVersionVisitor)
    }
}

/// A complete v1 trajectory-solve request.
///
/// All dimensional values use SI units named in the field. Every section is mandatory, even
/// when all fields in a section have defaults; this keeps the top-level request shape explicit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolveRequestV1 {
    /// Must be exactly [`SOLVE_JSON_SCHEMA_VERSION_V1`].
    pub schema_version: SchemaVersionV1,
    pub projectile: ProjectileV1,
    pub rifle: RifleV1,
    pub shot: ShotV1,
    pub atmosphere: AtmosphereV1,
    pub wind: WindV1,
    pub solver: SolverV1,
    pub effects: EffectsV1,
    pub sampling: SamplingV1,
    /// Optional reticle hold-point request (MBA-1361). Absent (the historical shape) leaves
    /// both the solve and the response byte-identical; present adds
    /// [`SolveSuccessV1::reticle_hold`] and echoes itself at
    /// [`ResolvedSolveRequestV1::reticle`] (0.33.0 decision-support Task 1), and nothing
    /// else. It is a pure post-processing read of the solved samples — it cannot change
    /// the trajectory.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub reticle: Option<ReticleRequestV1>,
}

/// Ask for a reticle hold point alongside the trajectory (MBA-1361).
///
/// `description` is the shared [`crate::reticle::ReticleDescription`] schema — the very same
/// JSON `ballistics reticle generate -o json` emits — so a service and a CLI user exchange
/// one representation. It is deliberately permissive about extra keys inside the
/// description (front ends carry render metadata there); the envelope around it stays
/// strict.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReticleRequestV1 {
    /// Range at which the hold is evaluated, meters. Must be inside the sampled trajectory.
    pub range_m: f64,
    /// The optic's magnification in use. Must be finite and greater than zero on both
    /// focal planes.
    pub magnification: f64,
    pub description: crate::reticle::ReticleDescription,
}

/// The reticle hold point for a solved trajectory (MBA-1361).
///
/// Angles are milliradians from the optical center: `down_mil` positive BELOW center,
/// `right_mil` positive to the shooter's RIGHT.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReticleHoldV1 {
    /// Echo of the requested evaluation range, meters.
    pub range_m: f64,
    /// Echo of the requested magnification.
    pub magnification: f64,
    pub down_mil: f64,
    pub right_mil: f64,
    /// The subtension scale applied to the marks (`reference_magnification / magnification`
    /// for SFP, exactly `1.0` for FFP).
    pub mark_scale: f64,
    /// Index into the request's `description.marks` of the nearest mark, in TRUE angular
    /// space.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_mark_index: Option<usize>,
    /// That mark's label, when it carries one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_mark_label: Option<String>,
    pub nearest_mark_distance_mil: f64,
    /// True when the hold has run outside the marked area (see
    /// [`crate::reticle::ReticleHold::off_reticle`] for the exact rule).
    pub off_reticle: bool,
}

/// Projectile inputs supported by solve-json v1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectileV1 {
    pub mass_kg: f64,
    pub diameter_m: f64,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub length_m: Option<f64>,
    pub drag_model: DragModelV1,
    pub ballistic_coefficient: f64,
}

/// Built-in reference-projectile drag models supported by solve-json v1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DragModelV1 {
    #[serde(rename = "G1")]
    G1,
    #[serde(rename = "G6")]
    G6,
    #[serde(rename = "G7")]
    G7,
    #[serde(rename = "G8")]
    G8,
}

/// Rifle and sight geometry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RifleV1 {
    pub muzzle_velocity_mps: f64,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub sight_height_m: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub muzzle_height_m: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub twist_rate_m_per_turn: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub twist_direction: Option<TwistDirectionV1>,
    /// Lateral offset between the sight axis and the bore axis, meters (MBA-1396,
    /// offset-mounted optics): positive = the sight sits RIGHT of the bore. The bore
    /// starts that far left of the line of sight, and a solved zero adds the windage
    /// convergence (`offset / zero_distance`) so the trajectory crosses the LOS laterally
    /// at the zero range. Omitted (the default) is byte-identical to pre-MBA-1396
    /// behavior. Echoed at [`ResolvedRifleV1::sight_offset_lateral_m`] when supplied
    /// (0.33.0 decision-support Task 1).
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub sight_offset_lateral_m: Option<f64>,
}

/// Direction of rifling twist as viewed from the breech toward the muzzle.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TwistDirectionV1 {
    Left,
    #[default]
    Right,
}

/// Shot geometry and termination range.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShotV1 {
    pub max_range_m: f64,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub zero_distance_m: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub muzzle_angle_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub aim_azimuth_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub shot_azimuth_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub shooting_angle_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub cant_angle_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub target_height_m: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub ground_threshold_m: Option<f64>,
    /// Deliberate vertical POI offset AT the zero range, meters (MBA-1359, Kestrel "zero
    /// height"): positive = the rifle is deliberately zeroed to impact HIGH by this much at
    /// `zero_distance_m`. Meaningful only when `zero_distance_m` is supplied. Omitted (the
    /// default) is byte-identical to pre-MBA-1359 behavior. Echoed at
    /// [`ResolvedShotV1::zero_poi_up_m`] when supplied (0.33.0 decision-support Task 1).
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub zero_poi_up_m: Option<f64>,
    /// Deliberate horizontal POI offset AT the zero range, meters (MBA-1359, Kestrel "zero
    /// offset"): positive = impacts RIGHT by this much at `zero_distance_m`. Same contract
    /// as `zero_poi_up_m`.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub zero_poi_right_m: Option<f64>,
    /// Which plane sampled `drop_m` values are referenced to (MBA-1403). `"los"` (the
    /// default when omitted) keeps the historical LOS-perpendicular drop byte-identically;
    /// `"target"` reports drop as vertical in the target plane — the LOS-perpendicular
    /// drop scaled by `1 / cos(shooting_angle_rad)` (JBM's "target plane" reference).
    /// Output-mode toggle only: the solved trajectory, `windage_m`, the summary block and
    /// zeroing semantics are unchanged. Echoed at [`ResolvedShotV1::drops_reference`]
    /// when supplied (0.33.0 decision-support Task 1).
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub drops_reference: Option<DropsReferenceV1>,
}

/// Wire values for [`ShotV1::drops_reference`] (MBA-1403).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DropsReferenceV1 {
    /// Drop measured perpendicular to the line of sight (the historical default).
    #[default]
    Los,
    /// Drop measured vertically in the target plane: LOS-perpendicular drop scaled by
    /// `1 / cos(shooting_angle_rad)`.
    Target,
}

/// Atmospheric station conditions.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtmosphereV1 {
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub altitude_m: Option<f64>,
    /// Authoritative station temperature, or `None` to resolve ICAO standard temperature at
    /// `altitude_m`.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub temperature_k: Option<f64>,
    /// Authoritative station pressure, or `None` to resolve ICAO standard pressure at
    /// `altitude_m`. Interpreted per `pressure_reference` (MBA-1397): `absolute` (the
    /// default) means this value already IS station pressure; `qnh` means it is a
    /// sea-level-corrected altimeter setting that must be reduced to station pressure at
    /// `altitude_m` before use.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub pressure_pa: Option<f64>,
    /// Whether `pressure_pa` is absolute station pressure or a sea-level-corrected altimeter
    /// setting (QNH, MBA-1397). `None` (the omitted-field default, and every request from
    /// before this field existed) means [`PressureReferenceV1::Absolute`] — byte-identical to
    /// pre-MBA-1397 behavior. Has no effect when `pressure_pa` is omitted: an omitted pressure
    /// resolves to the ICAO standard station pressure either way.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub pressure_reference: Option<PressureReferenceV1>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub relative_humidity: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub latitude_rad: Option<f64>,
}

/// Whether an `AtmosphereV1.pressure_pa` value is absolute station pressure or a sea-level-
/// corrected altimeter setting (QNH) that must be reduced to station pressure before use
/// (MBA-1397). Mirrors [`crate::atmosphere::PressureReferenceMode`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PressureReferenceV1 {
    #[default]
    Absolute,
    Qnh,
}

/// Constant or downrange-segmented wind.
///
/// Omit all fields for still air. A constant wind supplies both `speed_mps` and
/// `direction_from_rad`. Segmented wind supplies `segments` and no constant-wind fields.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindV1 {
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub speed_mps: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub direction_from_rad: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub vertical_speed_mps: Option<f64>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub segments: Option<Vec<WindSegmentV1>>,
    /// Which frame every wind direction in this request is entered in (MBA-1368).
    /// Omitted (or `"shooter"`) = shooter-relative wind-FROM radians, byte-identical to
    /// pre-1368 behavior. `"compass"` = earth-fixed bearings (0 = north) — the constant
    /// `direction_from_rad` AND every segment's — derived shooter-relative at resolve
    /// time as `bearing - shot.shot_azimuth_rad` (normalized to [0, 2π)); the RESOLVED
    /// wind echo therefore reports the converted shooter-relative direction (the QNH
    /// fold-into-resolved-value precedent). Compass mode requires an explicit
    /// `shot.shot_azimuth_rad` (a hard error otherwise, never a silent
    /// treat-as-shooter-relative). The mode itself is echoed at
    /// [`ResolvedConstantWindV1::wind_reference`] /
    /// [`ResolvedSegmentedWindV1::wind_reference`] when supplied (0.33.0
    /// decision-support Task 1).
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub wind_reference: Option<WindReferenceV1>,
}

/// Wire values for [`WindV1::wind_reference`] (MBA-1368).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WindReferenceV1 {
    /// Shooter-relative wind-FROM directions (the historical default).
    #[default]
    Shooter,
    /// Earth-fixed compass bearings, re-referenced against `shot.shot_azimuth_rad`.
    Compass,
}

/// One wind segment, active through `until_distance_m`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindSegmentV1 {
    pub until_distance_m: f64,
    pub speed_mps: f64,
    pub direction_from_rad: f64,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub vertical_speed_mps: Option<f64>,
}

/// Numerical integration configuration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolverV1 {
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub method: Option<SolverMethodV1>,
    /// Fixed step for RK4 and Euler. RK45 owns its adaptive step size.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub time_step_s: Option<f64>,
}

/// Integration algorithms exposed by solve-json v1.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolverMethodV1 {
    Euler,
    Rk4,
    #[default]
    Rk45,
}

/// Optional physical effects supported by solve-json v1.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectsV1 {
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub magnus: Option<bool>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub coriolis: Option<bool>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub enhanced_spin_drift: Option<bool>,
}

/// Regular downrange output sampling configuration.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SamplingV1 {
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub interval_m: Option<f64>,
}

/// A solve request after the service has applied every documented v1 default.
///
/// This is intentionally distinct from [`SolveRequestV1`]. Request DTOs retain whether a
/// defaulted field was omitted, while resolved DTOs require a concrete value for every default
/// that affected the solve. Semantically optional inputs remain optional.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedSolveRequestV1 {
    pub schema_version: SchemaVersionV1,
    pub projectile: ResolvedProjectileV1,
    pub rifle: ResolvedRifleV1,
    pub shot: ResolvedShotV1,
    pub atmosphere: ResolvedAtmosphereV1,
    pub wind: ResolvedWindV1,
    pub solver: ResolvedSolverV1,
    pub effects: ResolvedEffectsV1,
    pub sampling: ResolvedSamplingV1,
    /// Echo of the reticle hold-point request (MBA-1361), when one was supplied. Present
    /// only when the raw request supplied it — completes the resolved request as a
    /// full description of the solve.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reticle: Option<ReticleRequestV1>,
}

/// Resolved projectile inputs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedProjectileV1 {
    pub mass_kg: f64,
    pub diameter_m: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub length_m: Option<f64>,
    pub drag_model: DragModelV1,
    pub ballistic_coefficient: f64,
}

/// Resolved rifle and sight geometry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedRifleV1 {
    pub muzzle_velocity_mps: f64,
    pub sight_height_m: f64,
    pub muzzle_height_m: f64,
    pub twist_rate_m_per_turn: f64,
    pub twist_direction: TwistDirectionV1,
    /// Lateral sight offset for offset-mounted optics. Echoed so the resolved request
    /// is a complete description of the solve (required for counterfactual re-solve).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sight_offset_lateral_m: Option<f64>,
}

/// Resolved shot geometry.
///
/// `zero_distance_m` records caller zeroing intent, while `muzzle_angle_rad` is always the
/// effective angle used by the engine. A zero-distance solve therefore populates both fields.
/// A raw request may also supply both together (0.33.0 decision-support Task 2's
/// `From<&ResolvedSolveRequestV1> for SolveRequestV1` always does, to round-trip a resolved
/// request): `muzzle_angle_rad` then wins outright and no zero search runs, so
/// `zero_distance_m` carries no effect beyond recording the original intent.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedShotV1 {
    pub max_range_m: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_distance_m: Option<f64>,
    pub muzzle_angle_rad: f64,
    pub aim_azimuth_rad: f64,
    pub shot_azimuth_rad: f64,
    pub shooting_angle_rad: f64,
    pub cant_angle_rad: f64,
    pub target_height_m: f64,
    pub ground_threshold_m: f64,
    /// Echo of the requested deliberate vertical POI offset at the zero range, meters
    /// (MBA-1359). Present only when the raw request supplied it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_poi_up_m: Option<f64>,
    /// Echo of the requested deliberate horizontal POI offset at the zero range, meters
    /// (MBA-1359). Present only when the raw request supplied it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_poi_right_m: Option<f64>,
    /// Echo of the requested drops-reference plane (MBA-1403). Present only when the raw
    /// request supplied it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub drops_reference: Option<DropsReferenceV1>,
}

/// Resolved station conditions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedAtmosphereV1 {
    pub altitude_m: f64,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub relative_humidity: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latitude_rad: Option<f64>,
    /// Echo of the requested pressure reference mode (MBA-1397). Present only when the raw
    /// request supplied it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_reference: Option<PressureReferenceV1>,
}

/// Resolved constant or segmented wind.
///
/// The untagged representation retains the request's object shape while making the selected wind
/// variant explicit in Rust. Still air is a constant wind with three zero values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ResolvedWindV1 {
    Constant(ResolvedConstantWindV1),
    Segmented(ResolvedSegmentedWindV1),
}

/// Resolved still-air or constant-wind values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedConstantWindV1 {
    pub speed_mps: f64,
    pub direction_from_rad: f64,
    pub vertical_speed_mps: f64,
    /// Echo of the requested wind-direction reference frame (MBA-1368). Present only when
    /// the raw request supplied it. `direction_from_rad` above is always already converted
    /// to shooter-relative, regardless of this value.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_reference: Option<WindReferenceV1>,
}

/// Resolved downrange wind segments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedSegmentedWindV1 {
    pub segments: Vec<ResolvedWindSegmentV1>,
    /// Echo of the requested wind-direction reference frame (MBA-1368). Present only when
    /// the raw request supplied it. Each segment's `direction_from_rad` above is always
    /// already converted to shooter-relative, regardless of this value.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wind_reference: Option<WindReferenceV1>,
}

/// One resolved wind segment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedWindSegmentV1 {
    pub until_distance_m: f64,
    pub speed_mps: f64,
    pub direction_from_rad: f64,
    pub vertical_speed_mps: f64,
}

/// Resolved numerical integration configuration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedSolverV1 {
    pub method: SolverMethodV1,
    pub time_step_s: f64,
}

/// Resolved physical-effect switches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedEffectsV1 {
    pub magnus: bool,
    pub coriolis: bool,
    pub enhanced_spin_drift: bool,
}

/// Resolved result-sampling configuration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedSamplingV1 {
    pub interval_m: f64,
}

/// A successful solve-json v1 response.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolveSuccessV1 {
    pub schema_version: SchemaVersionV1,
    pub engine_version: String,
    pub status: SuccessStatusV1,
    pub resolved_request: ResolvedSolveRequestV1,
    #[serde(default)]
    pub assumptions: Vec<SolveNoticeV1>,
    #[serde(default)]
    pub warnings: Vec<SolveNoticeV1>,
    pub summary: SolveSummaryV1,
    #[serde(default, serialize_with = "serialize_solve_samples_v1")]
    pub samples: Vec<TrajectorySampleV1>,
    /// The reticle hold point (MBA-1361), present only when the request carried a
    /// `reticle` block. Every response that predates the field, and every request without
    /// one, is byte-identical.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reticle_hold: Option<ReticleHoldV1>,
}

fn serialize_solve_samples_v1<S>(
    samples: &[TrajectorySampleV1],
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    if samples.len() > MAX_SOLVE_JSON_SAMPLES_V1 {
        return Err(serde::ser::Error::custom(format_args!(
            "solve-json v1 response sample limit of {MAX_SOLVE_JSON_SAMPLES_V1} exceeded: response has {} samples",
            samples.len()
        )));
    }
    samples.serialize(serializer)
}

impl SolveSuccessV1 {
    /// Validate service-level response limits immediately before serialization.
    ///
    /// The exact limit is accepted. A response with more samples returns a structured
    /// [`SolveErrorCodeV1::ResourceLimit`] error at the sampling interval that requested the
    /// oversized result. The response is never silently truncated.
    pub fn validate_for_serialization(&self) -> Result<(), SolveErrorEnvelopeV1> {
        if self.samples.len() <= MAX_SOLVE_JSON_SAMPLES_V1 {
            return Ok(());
        }

        Err(SolveErrorEnvelopeV1::new(
            SolveErrorV1::new(
                SolveErrorCodeV1::ResourceLimit,
                format!(
                    "solve-json v1 response sample limit of {MAX_SOLVE_JSON_SAMPLES_V1} exceeded: response has {} samples",
                    self.samples.len()
                ),
            )
            .at_path("$.sampling.interval_m"),
        ))
    }
}

/// The success discriminator serialized in a response envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SuccessStatusV1 {
    Ok,
}

/// A machine-readable assumption or warning.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolveNoticeV1 {
    pub code: String,
    pub message: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
}

/// Aggregate observations for a completed solve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolveSummaryV1 {
    pub actual_range_m: f64,
    /// Greatest world-vertical projectile height above the request's local ground/reference
    /// datum. This is not shot-frame Y and is not height above the line of sight.
    pub maximum_height_m: f64,
    pub time_of_flight_s: f64,
    pub terminal_speed_mps: f64,
    pub terminal_energy_j: f64,
    /// Muzzle gyroscopic stability factor Sg evaluated with the resolved projectile, muzzle
    /// velocity, twist, and station atmosphere; absent only when it cannot be calculated.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stability_factor: Option<f64>,
    /// Signed gyroscopic spin-drift contribution at the terminal sample, positive to the
    /// shooter's right. This excludes wind drift and is absent when the effect is disabled or
    /// cannot be calculated.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spin_drift_m: Option<f64>,
    /// Equivalent horizontal range for an inclined shot (MBA-1395): the flat-fire range
    /// whose angular elevation correction against the same zero matches the inclined
    /// solution's at the terminal range — the BDC "shoot-to" range (SIG AMR / Leica EHR
    /// style, angular-match inversion, not the rifleman's-rule cosine). Present only when
    /// `shot.shooting_angle_rad != 0`, a `zero_distance_m` was solved, and the inverse is
    /// well-defined (terminal past the zero range, positive correction); absent otherwise
    /// — responses that predate the field, and every flat solve, are byte-identical.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equivalent_horizontal_range_m: Option<f64>,
    pub termination: TerminationReasonV1,
}

/// Why the trajectory ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TerminationReasonV1 {
    MaxRange,
    GroundThreshold,
    TimeLimit,
    VelocityFloor,
}

/// One regularly sampled trajectory observation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrajectorySampleV1 {
    pub distance_m: f64,
    pub time_s: f64,
    pub speed_mps: f64,
    pub energy_j: f64,
    /// Positive means below the line of sight.
    pub drop_m: f64,
    /// Positive means right of the line of sight from the shooter's perspective.
    pub windage_m: f64,
    pub mach: f64,
    #[serde(default)]
    pub flags: Vec<SampleFlagV1>,
}

/// Stable annotations for trajectory samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SampleFlagV1 {
    Transonic,
    Subsonic,
    Terminal,
    GroundThreshold,
}

/// A failed solve-json v1 response.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolveErrorEnvelopeV1 {
    pub schema_version: SchemaVersionV1,
    pub status: ErrorStatusV1,
    pub error: SolveErrorV1,
}

impl SolveErrorEnvelopeV1 {
    /// Wrap a protocol error in a v1 envelope.
    pub fn new(error: SolveErrorV1) -> Self {
        Self {
            schema_version: SchemaVersionV1,
            status: ErrorStatusV1::Error,
            error,
        }
    }
}

/// The error discriminator serialized in a response envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ErrorStatusV1 {
    Error,
}

/// A stable, machine-readable protocol error.
///
/// Location state is private so callers cannot construct a path together with a source location,
/// a partial line/column pair, or a zero source coordinate. The JSON representation retains the
/// flat `path`, `line`, and `column` fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolveErrorV1 {
    pub code: SolveErrorCodeV1,
    pub message: String,
    location: SolveErrorLocationV1,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SolveErrorLocationV1 {
    None,
    Path(String),
    Source {
        line: NonZeroUsize,
        column: NonZeroUsize,
    },
}

/// Why a source location could not be attached to a protocol error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveErrorLocationErrorV1 {
    ZeroLine,
    ZeroColumn,
}

impl fmt::Display for SolveErrorLocationErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroLine => formatter.write_str("error line must be one-based"),
            Self::ZeroColumn => formatter.write_str("error column must be one-based"),
        }
    }
}

impl std::error::Error for SolveErrorLocationErrorV1 {}

impl SolveErrorV1 {
    /// Construct an error without input-location information.
    pub fn new(code: SolveErrorCodeV1, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            location: SolveErrorLocationV1::None,
        }
    }

    /// Attach a JSONPath-like input location, replacing any source location.
    pub fn at_path(mut self, path: impl Into<String>) -> Self {
        self.location = SolveErrorLocationV1::Path(path.into());
        self
    }

    /// Attach a one-based malformed-JSON source location, replacing any JSON path.
    pub fn at_location(
        mut self,
        line: usize,
        column: usize,
    ) -> Result<Self, SolveErrorLocationErrorV1> {
        let line = NonZeroUsize::new(line).ok_or(SolveErrorLocationErrorV1::ZeroLine)?;
        let column = NonZeroUsize::new(column).ok_or(SolveErrorLocationErrorV1::ZeroColumn)?;
        self.location = SolveErrorLocationV1::Source { line, column };
        Ok(self)
    }

    /// Attach a parser location while keeping the public wire coordinates one-based.
    ///
    /// `serde_json` reports column zero for some end-of-file errors. Normalize either zero
    /// coordinate to one instead of allowing malformed input to trigger a panic.
    fn at_parser_location(mut self, line: usize, column: usize) -> Self {
        self.location = SolveErrorLocationV1::Source {
            line: NonZeroUsize::new(line).unwrap_or(NonZeroUsize::MIN),
            column: NonZeroUsize::new(column).unwrap_or(NonZeroUsize::MIN),
        };
        self
    }

    /// Return the JSONPath-like input location, when present.
    pub fn path(&self) -> Option<&str> {
        match &self.location {
            SolveErrorLocationV1::Path(path) => Some(path),
            SolveErrorLocationV1::None | SolveErrorLocationV1::Source { .. } => None,
        }
    }

    /// Return the one-based malformed-JSON source line, when present.
    pub fn line(&self) -> Option<usize> {
        match &self.location {
            SolveErrorLocationV1::Source { line, .. } => Some(line.get()),
            SolveErrorLocationV1::None | SolveErrorLocationV1::Path(_) => None,
        }
    }

    /// Return the one-based malformed-JSON source column, when present.
    pub fn column(&self) -> Option<usize> {
        match &self.location {
            SolveErrorLocationV1::Source { column, .. } => Some(column.get()),
            SolveErrorLocationV1::None | SolveErrorLocationV1::Path(_) => None,
        }
    }
}

impl Serialize for SolveErrorV1 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("SolveErrorV1", 5)?;
        state.serialize_field("code", &self.code)?;
        state.serialize_field("message", &self.message)?;
        state.serialize_field("path", &self.path())?;
        state.serialize_field("line", &self.line())?;
        state.serialize_field("column", &self.column())?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for SolveErrorV1 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct SolveErrorWireV1 {
            code: SolveErrorCodeV1,
            message: String,
            path: Option<String>,
            line: Option<usize>,
            column: Option<usize>,
        }

        let wire = SolveErrorWireV1::deserialize(deserializer)?;
        let location = match (wire.path, wire.line, wire.column) {
            (None, None, None) => SolveErrorLocationV1::None,
            (Some(path), None, None) => SolveErrorLocationV1::Path(path),
            (None, Some(line), Some(column)) => {
                let line = NonZeroUsize::new(line)
                    .ok_or_else(|| de::Error::custom("error line must be one-based"))?;
                let column = NonZeroUsize::new(column)
                    .ok_or_else(|| de::Error::custom("error column must be one-based"))?;
                SolveErrorLocationV1::Source { line, column }
            }
            _ => {
                return Err(de::Error::custom(
                    "error location must contain either path or both line and column",
                ));
            }
        };

        Ok(Self {
            code: wire.code,
            message: wire.message,
            location,
        })
    }
}

/// Stable solve-json v1 error categories.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolveErrorCodeV1 {
    InvalidJson,
    UnsupportedSchemaVersion,
    UnknownField,
    MissingField,
    InvalidValue,
    ConflictingFields,
    ResourceLimit,
    SolveFailed,
    IoError,
    InternalError,
}

/// Decode and structurally validate one solve-json v1 request.
///
/// In addition to Serde's `deny_unknown_fields` enforcement, this entry point attaches exact
/// JSONPath-like locations to unknown fields, missing required fields, invalid enums, and an
/// unsupported schema version. Transport implementations should use this function rather than
/// deserializing [`SolveRequestV1`] directly when they need a protocol error envelope.
pub fn decode_solve_request_v1(input: &str) -> Result<SolveRequestV1, SolveErrorEnvelopeV1> {
    let value: Value = serde_json::from_str(input).map_err(|error| {
        let error = SolveErrorV1::new(SolveErrorCodeV1::InvalidJson, error.to_string())
            .at_parser_location(error.line(), error.column());
        envelope(error)
    })?;

    validate_request_shape(&value)?;

    let request: SolveRequestV1 = serde_json::from_value(value).map_err(|error| {
        envelope(SolveErrorV1::new(SolveErrorCodeV1::InvalidValue, error.to_string()).at_path("$"))
    })?;

    validate_request_ranges(&request)?;

    Ok(request)
}

/// Physically meaningful bounds for solve-json v1's numeric inputs (MBA-1413).
///
/// Deliberately enormous — the point is to exclude values that are not projectiles at all, not
/// to police calibers. Every bound below admits everything from an airgun pellet to naval
/// artillery, so a request rejected here was never going to produce a meaningful answer.
///
/// Found by fuzzing, which reached a request declaring a `mass_kg` of 1.1e-66 — about
/// 10^40 times lighter than a proton. The engine accepted it and solved. The visible symptom was
/// narrower: at that magnitude `serde_json`'s float parser is a bit off re-reading its own
/// output, so the request did not survive a JSON round trip and broke the protocol's stability
/// invariant. Tightening the float handling would have addressed the symptom; the actual defect
/// is that a number no projectile could have was accepted as a projectile.
mod limits {
    /// 1 mg to 100 kg. A .17 cal pellet is ~500 mg; a 16-inch naval shell is ~1200 kg, so the
    /// upper bound is the one a caller might conceivably reach — it is set past small artillery
    /// on purpose and can be raised without ceremony.
    pub const MASS_KG: (f64, f64) = (1.0e-6, 100.0);
    /// 0.1 mm to 1 m.
    pub const DIAMETER_M: (f64, f64) = (1.0e-4, 1.0);
    /// 0.1 mm to 10 m. Only checked when supplied.
    pub const LENGTH_M: (f64, f64) = (1.0e-4, 10.0);
    /// Ballistic coefficient, lb/in². Real values run ~0.1–1.5; the bound is far past both ends.
    pub const BALLISTIC_COEFFICIENT: (f64, f64) = (1.0e-4, 100.0);

    // Deliberately NOT bounded here: muzzle_velocity_mps and max_range_m. The MCP server
    // documents — and tests — a split where a structurally valid request the engine cannot
    // solve returns a tool error rather than a protocol error, and an absurd muzzle velocity is
    // its worked example of that case. Bounding those two fields here would reclassify that
    // example as invalid params and change a contract MCP clients may rely on, which is a
    // separate decision from rejecting values that cannot describe a projectile.
}

fn require_range(
    value: f64,
    (min, max): (f64, f64),
    path: &str,
) -> Result<(), SolveErrorEnvelopeV1> {
    if !value.is_finite() {
        return Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            format!("{path} must be a finite number"),
            path,
        ));
    }
    if value < min || value > max {
        return Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            format!("{path} must be between {min} and {max}, got {value}"),
            path,
        ));
    }
    Ok(())
}

/// Reject requests whose numbers cannot describe a projectile, before any solve runs.
///
/// Runs after deserialization so the fields are typed and each error can carry the exact
/// JSONPath the caller used, matching how the shape errors above report.
fn validate_request_ranges(request: &SolveRequestV1) -> Result<(), SolveErrorEnvelopeV1> {
    require_range(
        request.projectile.mass_kg,
        limits::MASS_KG,
        "$.projectile.mass_kg",
    )?;
    require_range(
        request.projectile.diameter_m,
        limits::DIAMETER_M,
        "$.projectile.diameter_m",
    )?;
    require_range(
        request.projectile.ballistic_coefficient,
        limits::BALLISTIC_COEFFICIENT,
        "$.projectile.ballistic_coefficient",
    )?;
    if let Some(length_m) = request.projectile.length_m {
        require_range(length_m, limits::LENGTH_M, "$.projectile.length_m")?;
    }
    Ok(())
}

fn envelope(error: SolveErrorV1) -> SolveErrorEnvelopeV1 {
    SolveErrorEnvelopeV1::new(error)
}

fn protocol_error(
    code: SolveErrorCodeV1,
    message: impl Into<String>,
    path: impl Into<String>,
) -> SolveErrorEnvelopeV1 {
    envelope(SolveErrorV1::new(code, message).at_path(path))
}

fn validate_request_shape(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let root = require_object(value, "$")?;

    // Dispatch on the version before applying any v1-only shape rules. A future-version
    // request may legitimately contain fields unknown to v1 and must still receive the stable
    // unsupported-version error rather than an order-dependent unknown-field error.
    validate_schema_version(root)?;
    validate_members(
        root,
        "$",
        &[
            "schema_version",
            "projectile",
            "rifle",
            "shot",
            "atmosphere",
            "wind",
            "solver",
            "effects",
            "sampling",
            // MBA-1361: optional, additive. Omitting it is the historical shape.
            "reticle",
        ],
        &[
            "schema_version",
            "projectile",
            "rifle",
            "shot",
            "atmosphere",
            "wind",
            "solver",
            "effects",
            "sampling",
        ],
    )?;

    validate_projectile(required_value(root, "projectile", "$")?)?;
    validate_rifle(required_value(root, "rifle", "$")?)?;
    validate_shot(required_value(root, "shot", "$")?)?;
    validate_atmosphere(required_value(root, "atmosphere", "$")?)?;
    validate_wind(required_value(root, "wind", "$")?)?;
    validate_solver(required_value(root, "solver", "$")?)?;
    validate_effects(required_value(root, "effects", "$")?)?;
    validate_sampling(required_value(root, "sampling", "$")?)?;
    if let Some(reticle) = root.get("reticle") {
        validate_reticle(reticle)?;
    }
    Ok(())
}

/// Shape-validate an optional `reticle` block (MBA-1361).
///
/// The ENVELOPE is strict (exactly `range_m`, `magnification`, `description`, all
/// required); the description itself is handed to the shared reticle schema, whose own
/// `validate()` runs at solve time. Splitting it that way keeps the wire contract tight
/// without making the engine the arbiter of a front end's render metadata.
fn validate_reticle(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.reticle";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &["range_m", "magnification", "description"],
        &["range_m", "magnification", "description"],
    )?;
    validate_required_numbers(object, path, &["range_m", "magnification"])?;
    require_object(required_value(object, "description", path)?, "$.reticle.description")?;
    Ok(())
}

fn validate_schema_version(root: &Map<String, Value>) -> Result<(), SolveErrorEnvelopeV1> {
    let value = required_value(root, "schema_version", "$")?;
    let version = if let Some(version) = value.as_i64() {
        i128::from(version)
    } else if let Some(version) = value.as_u64() {
        i128::from(version)
    } else {
        return Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            "schema_version must be the integer 1",
            "$.schema_version",
        ));
    };
    if version != i128::from(SOLVE_JSON_SCHEMA_VERSION_V1) {
        return Err(protocol_error(
            SolveErrorCodeV1::UnsupportedSchemaVersion,
            format!(
                "unsupported schema_version {version}; expected {SOLVE_JSON_SCHEMA_VERSION_V1}"
            ),
            "$.schema_version",
        ));
    }
    Ok(())
}

fn validate_projectile(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.projectile";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &[
            "mass_kg",
            "diameter_m",
            "length_m",
            "drag_model",
            "ballistic_coefficient",
        ],
        &[
            "mass_kg",
            "diameter_m",
            "drag_model",
            "ballistic_coefficient",
        ],
    )?;
    validate_required_numbers(
        object,
        path,
        &["mass_kg", "diameter_m", "ballistic_coefficient"],
    )?;
    validate_optional_number(object, path, "length_m")?;
    validate_string_enum(
        required_value(object, "drag_model", path)?,
        "$.projectile.drag_model",
        &["G1", "G6", "G7", "G8"],
    )
}

fn validate_rifle(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.rifle";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &[
            "muzzle_velocity_mps",
            "sight_height_m",
            "muzzle_height_m",
            "twist_rate_m_per_turn",
            "twist_direction",
            "sight_offset_lateral_m",
        ],
        &["muzzle_velocity_mps"],
    )?;
    validate_required_numbers(object, path, &["muzzle_velocity_mps"])?;
    validate_optional_numbers(
        object,
        path,
        &[
            "sight_height_m",
            "muzzle_height_m",
            "twist_rate_m_per_turn",
            "sight_offset_lateral_m",
        ],
    )?;
    if let Some(direction) = object.get("twist_direction") {
        validate_string_enum(direction, "$.rifle.twist_direction", &["left", "right"])?;
    }
    Ok(())
}

fn validate_shot(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.shot";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &[
            "max_range_m",
            "zero_distance_m",
            "muzzle_angle_rad",
            "aim_azimuth_rad",
            "shot_azimuth_rad",
            "shooting_angle_rad",
            "cant_angle_rad",
            "target_height_m",
            "ground_threshold_m",
            "zero_poi_up_m",
            "zero_poi_right_m",
            "drops_reference",
        ],
        &["max_range_m"],
    )?;
    validate_required_numbers(object, path, &["max_range_m"])?;
    validate_optional_number(object, path, "zero_distance_m")?;
    validate_optional_number(object, path, "muzzle_angle_rad")?;
    validate_optional_numbers(
        object,
        path,
        &[
            "aim_azimuth_rad",
            "shot_azimuth_rad",
            "shooting_angle_rad",
            "cant_angle_rad",
            "target_height_m",
            "ground_threshold_m",
            "zero_poi_up_m",
            "zero_poi_right_m",
        ],
    )?;
    // MBA-1403: string enum, not a number — same shape-validation pattern as
    // $.rifle.twist_direction.
    if let Some(reference) = object.get("drops_reference") {
        validate_string_enum(reference, "$.shot.drops_reference", &["los", "target"])?;
    }
    Ok(())
}

fn validate_atmosphere(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.atmosphere";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &[
            "altitude_m",
            "temperature_k",
            "pressure_pa",
            "pressure_reference",
            "relative_humidity",
            "latitude_rad",
        ],
        &[],
    )?;
    validate_optional_numbers(
        object,
        path,
        &[
            "altitude_m",
            "temperature_k",
            "pressure_pa",
            "relative_humidity",
            "latitude_rad",
        ],
    )?;
    // MBA-1397: string enum, not a number -- same shape-validation pattern as
    // $.shot.drops_reference / $.wind.wind_reference. This field has existed on
    // AtmosphereV1 and been consumed by resolve_atmosphere since MBA-1397, but was never
    // added to this hand-maintained allowlist, so decode_solve_request_v1 rejected it as
    // an unknown field even though direct SolveRequestV1 construction always accepted it.
    if let Some(reference) = object.get("pressure_reference") {
        validate_string_enum(
            reference,
            "$.atmosphere.pressure_reference",
            &["absolute", "qnh"],
        )?;
    }
    Ok(())
}

fn validate_wind(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.wind";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &[
            "speed_mps",
            "direction_from_rad",
            "vertical_speed_mps",
            "segments",
            "wind_reference",
        ],
        &[],
    )?;
    for field in ["speed_mps", "direction_from_rad", "vertical_speed_mps"] {
        validate_optional_number(object, path, field)?;
    }
    // MBA-1368: wind_reference is a closed string enum (shooter | compass).
    if let Some(reference) = object.get("wind_reference") {
        validate_string_enum(reference, "$.wind.wind_reference", &["shooter", "compass"])?;
    }

    if let Some(segments) = object.get("segments") {
        let Some(segments) = segments.as_array() else {
            return Err(protocol_error(
                SolveErrorCodeV1::InvalidValue,
                "segments must be an array",
                "$.wind.segments",
            ));
        };
        for (index, segment) in segments.iter().enumerate() {
            let segment_path = format!("$.wind.segments[{index}]");
            let segment = require_object(segment, &segment_path)?;
            validate_members(
                segment,
                &segment_path,
                &[
                    "until_distance_m",
                    "speed_mps",
                    "direction_from_rad",
                    "vertical_speed_mps",
                ],
                &["until_distance_m", "speed_mps", "direction_from_rad"],
            )?;
            validate_required_numbers(
                segment,
                &segment_path,
                &["until_distance_m", "speed_mps", "direction_from_rad"],
            )?;
            validate_optional_number(segment, &segment_path, "vertical_speed_mps")?;
        }
    }
    Ok(())
}

fn validate_solver(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.solver";
    let object = require_object(value, path)?;
    validate_members(object, path, &["method", "time_step_s"], &[])?;
    validate_optional_number(object, path, "time_step_s")?;
    if let Some(method) = object.get("method") {
        validate_string_enum(method, "$.solver.method", &["euler", "rk4", "rk45"])?;
    }
    Ok(())
}

fn validate_effects(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.effects";
    let object = require_object(value, path)?;
    validate_members(
        object,
        path,
        &["magnus", "coriolis", "enhanced_spin_drift"],
        &[],
    )?;
    validate_optional_booleans(object, path, &["magnus", "coriolis", "enhanced_spin_drift"])?;

    if object.get("magnus").and_then(Value::as_bool) == Some(true)
        && object.get("enhanced_spin_drift").and_then(Value::as_bool) == Some(true)
    {
        return Err(protocol_error(
            SolveErrorCodeV1::ConflictingFields,
            "magnus and enhanced_spin_drift cannot both be enabled",
            "$.effects",
        ));
    }

    Ok(())
}

fn validate_sampling(value: &Value) -> Result<(), SolveErrorEnvelopeV1> {
    let path = "$.sampling";
    let object = require_object(value, path)?;
    validate_members(object, path, &["interval_m"], &[])?;
    validate_optional_number(object, path, "interval_m")
}

fn require_object<'a>(
    value: &'a Value,
    path: &str,
) -> Result<&'a Map<String, Value>, SolveErrorEnvelopeV1> {
    value
        .as_object()
        .ok_or_else(|| protocol_error(SolveErrorCodeV1::InvalidValue, "expected an object", path))
}

fn required_value<'a>(
    object: &'a Map<String, Value>,
    field: &str,
    parent_path: &str,
) -> Result<&'a Value, SolveErrorEnvelopeV1> {
    object.get(field).ok_or_else(|| {
        protocol_error(
            SolveErrorCodeV1::MissingField,
            format!("missing required field `{field}`"),
            child_path(parent_path, field),
        )
    })
}

fn validate_members(
    object: &Map<String, Value>,
    path: &str,
    allowed: &[&str],
    required: &[&str],
) -> Result<(), SolveErrorEnvelopeV1> {
    if let Some(field) = object
        .keys()
        .find(|field| !allowed.contains(&field.as_str()))
    {
        return Err(protocol_error(
            SolveErrorCodeV1::UnknownField,
            format!("unknown field `{field}`"),
            child_path(path, field),
        ));
    }

    if let Some(field) = required.iter().find(|field| !object.contains_key(**field)) {
        return Err(protocol_error(
            SolveErrorCodeV1::MissingField,
            format!("missing required field `{field}`"),
            child_path(path, field),
        ));
    }
    Ok(())
}

fn validate_string_enum(
    value: &Value,
    path: &str,
    allowed: &[&str],
) -> Result<(), SolveErrorEnvelopeV1> {
    let Some(value) = value.as_str() else {
        return Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            "expected a string enum value",
            path,
        ));
    };
    if !allowed.contains(&value) {
        return Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            format!(
                "invalid value `{value}`; expected one of {}",
                allowed.join(", ")
            ),
            path,
        ));
    }
    Ok(())
}

fn validate_required_numbers(
    object: &Map<String, Value>,
    parent_path: &str,
    fields: &[&str],
) -> Result<(), SolveErrorEnvelopeV1> {
    for field in fields {
        let value = required_value(object, field, parent_path)?;
        validate_number(value, &child_path(parent_path, field))?;
    }
    Ok(())
}

fn validate_optional_numbers(
    object: &Map<String, Value>,
    parent_path: &str,
    fields: &[&str],
) -> Result<(), SolveErrorEnvelopeV1> {
    for field in fields {
        validate_optional_number(object, parent_path, field)?;
    }
    Ok(())
}

fn validate_optional_number(
    object: &Map<String, Value>,
    parent_path: &str,
    field: &str,
) -> Result<(), SolveErrorEnvelopeV1> {
    if let Some(value) = object.get(field) {
        validate_number(value, &child_path(parent_path, field))?;
    }
    Ok(())
}

fn validate_number(value: &Value, path: &str) -> Result<(), SolveErrorEnvelopeV1> {
    if value.is_number() {
        Ok(())
    } else {
        Err(protocol_error(
            SolveErrorCodeV1::InvalidValue,
            "expected a number",
            path,
        ))
    }
}

fn validate_optional_booleans(
    object: &Map<String, Value>,
    parent_path: &str,
    fields: &[&str],
) -> Result<(), SolveErrorEnvelopeV1> {
    for field in fields {
        if let Some(value) = object.get(*field) {
            if !value.is_boolean() {
                return Err(protocol_error(
                    SolveErrorCodeV1::InvalidValue,
                    "expected a boolean",
                    child_path(parent_path, field),
                ));
            }
        }
    }
    Ok(())
}

fn child_path(parent: &str, field: &str) -> String {
    format!("{parent}.{field}")
}

/// MBA-1413: physical bounds on the projectile fields.
#[cfg(test)]
mod request_range_tests {
    use super::*;

    /// A complete, ordinary request. Every range test below starts from this and perturbs one
    /// field, so a bound that accidentally rejects real data fails loudly here first.
    fn valid_request_json() -> String {
        r#"{"schema_version":1,
            "projectile":{"mass_kg":0.01134,"diameter_m":0.00782,"length_m":0.031,
                          "drag_model":"G7","ballistic_coefficient":0.243},
            "rifle":{"muzzle_velocity_mps":823.0},
            "shot":{"max_range_m":1000.0},
            "atmosphere":{},"wind":{},"solver":{},"effects":{},"sampling":{}}"#
            .to_string()
    }

    fn decode_err(json: &str) -> SolveErrorEnvelopeV1 {
        decode_solve_request_v1(json).expect_err("request should have been rejected")
    }

    #[test]
    fn an_ordinary_request_still_decodes() {
        decode_solve_request_v1(&valid_request_json()).expect("a real load must not be rejected");
    }

    /// The exact input cargo-fuzz found (MBA-1413). It declared a projectile mass of about
    /// 1.1e-66 kg — roughly 10^40 times lighter than a proton — and the engine accepted and
    /// solved it. The visible symptom was that the request did not survive a JSON round trip,
    /// because serde_json's float parser is a bit off re-reading its own output at that
    /// magnitude; the actual defect was accepting the number at all.
    #[test]
    fn the_fuzz_reproducer_is_now_a_clean_typed_rejection() {
        let reproducer = r#"{"schema_version":1,
            "projectile":{"mass_kg":0.011366666667e-64,"diameter_m":0.00782,
                          "drag_model":"G7","ballistic_coefficient":1.2e2},
            "rifle":{"muzzle_velocity_mps":823.0},
            "shot":{"max_range_m":100.0},
            "atmosphere":{},"wind":{},"solver":{},"effects":{},"sampling":{}}"#;

        let envelope = decode_err(reproducer);
        assert_eq!(envelope.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(envelope.error.path(), Some("$.projectile.mass_kg"));
    }

    /// An error envelope must itself round-trip, since that is the invariant the fuzz target
    /// asserts on the rejection branch.
    #[test]
    fn the_rejection_envelope_round_trips() {
        let envelope = decode_err(
            r#"{"schema_version":1,
                "projectile":{"mass_kg":1.0e-66,"diameter_m":0.00782,
                              "drag_model":"G7","ballistic_coefficient":0.243},
                "rifle":{"muzzle_velocity_mps":823.0},
                "shot":{"max_range_m":100.0},
                "atmosphere":{},"wind":{},"solver":{},"effects":{},"sampling":{}}"#,
        );
        let encoded = serde_json::to_string(&envelope).expect("serialize");
        let decoded: SolveErrorEnvelopeV1 = serde_json::from_str(&encoded).expect("deserialize");
        assert_eq!(decoded, envelope);
    }

    #[test]
    fn each_bounded_field_reports_its_own_path() {
        for (field, bad_value, path) in [
            ("mass_kg", "1.0e-66", "$.projectile.mass_kg"),
            ("diameter_m", "1.0e-9", "$.projectile.diameter_m"),
            (
                "ballistic_coefficient",
                "1.0e6",
                "$.projectile.ballistic_coefficient",
            ),
            ("length_m", "1.0e-9", "$.projectile.length_m"),
        ] {
            let json = valid_request_json().replace(
                &format!("\"{field}\":{}", default_for(field)),
                &format!("\"{field}\":{bad_value}"),
            );
            let envelope = decode_err(&json);
            assert_eq!(
                envelope.error.path(),
                Some(path),
                "wrong path for {field}"
            );
            assert_eq!(envelope.error.code, SolveErrorCodeV1::InvalidValue);
        }
    }

    fn default_for(field: &str) -> &'static str {
        match field {
            "mass_kg" => "0.01134",
            "diameter_m" => "0.00782",
            "ballistic_coefficient" => "0.243",
            "length_m" => "0.031",
            other => panic!("no default recorded for {other}"),
        }
    }

    /// Muzzle velocity is deliberately unbounded: the MCP server documents and tests a split
    /// where a structurally valid request the engine cannot solve returns a tool error rather
    /// than a protocol error, using an absurd muzzle velocity as its example.
    #[test]
    fn muzzle_velocity_is_deliberately_left_unbounded() {
        let json = valid_request_json().replace("823.0", "1.0e308");
        decode_solve_request_v1(&json)
            .expect("muzzle velocity must stay a solve-time concern, not a protocol one");
    }
}
