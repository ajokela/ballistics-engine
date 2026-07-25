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
}

/// Resolved shot geometry.
///
/// `zero_distance_m` records caller zeroing intent, while `muzzle_angle_rad` is always the
/// effective angle used by the engine. A zero-distance solve therefore populates both fields.
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
}

/// Resolved downrange wind segments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedSegmentedWindV1 {
    pub segments: Vec<ResolvedWindSegmentV1>,
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

    serde_json::from_value(value).map_err(|error| {
        envelope(SolveErrorV1::new(SolveErrorCodeV1::InvalidValue, error.to_string()).at_path("$"))
    })
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
        ],
        &["muzzle_velocity_mps"],
    )?;
    validate_required_numbers(object, path, &["muzzle_velocity_mps"])?;
    validate_optional_numbers(
        object,
        path,
        &["sight_height_m", "muzzle_height_m", "twist_rate_m_per_turn"],
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
        ],
    )
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
    )
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
        ],
        &[],
    )?;
    for field in ["speed_mps", "direction_from_rad", "vertical_speed_mps"] {
        validate_optional_number(object, path, field)?;
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
