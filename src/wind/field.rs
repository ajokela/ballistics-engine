//! Deterministic wind fields: the [`WindField`] trait and its implementations.
//!
//! Every field here is stateless: the wind depends only on the time and position passed in, so
//! integrator substeps, snapshots and seeded reproduction give identical results. Nothing reads
//! a clock, the network or a random source.
//!
//! Conventions match the rest of [`crate::wind`]: world-frame McCoy coordinates (x downrange,
//! y up, z right), speeds in m/s, and wind-FROM directions in degrees with 0° a headwind and 90°
//! wind from the right. A field built from a direction gives the same vector as
//! [`wind_vector`] for that speed and direction.

use super::{wind_vector, WindSock};
use core::f64::consts::PI;
use nalgebra::Vector3;
use thiserror::Error;

/// Earth's rotation rate, rad/s (the value the crate's Coriolis paths use).
const EARTH_ROTATION_RAD_PER_S: f64 = 7.2921159e-5;

/// A deterministic, stateless wind field: same time and position always
/// give the same wind, so integrator substeps, snapshots and seeded
/// reproduction are exact.
pub trait WindField {
    /// World-frame air velocity (m/s) at a time and position.
    fn wind(&self, time_s: f64, position_m: Vector3<f64>) -> Vector3<f64>;
}

impl<F: WindField + ?Sized> WindField for &F {
    fn wind(&self, time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        (**self).wind(time_s, position_m)
    }
}

impl<F: WindField + ?Sized> WindField for Box<F> {
    fn wind(&self, time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        (**self).wind(time_s, position_m)
    }
}

/// The existing segmented wind as a field: the wind at downrange position `x`, exactly
/// [`WindSock::vector_for_range_stateless`]. Time, height and cross-range position are ignored.
impl WindField for WindSock {
    fn wind(&self, _time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        self.vector_for_range_stateless(position_m.x)
    }
}

/// Uniform wind everywhere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantField(pub Vector3<f64>);

impl WindField for ConstantField {
    fn wind(&self, _time_s: f64, _position_m: Vector3<f64>) -> Vector3<f64> {
        self.0
    }
}

/// Reads a field at a shifted height, so a field and its caller can use different height
/// datums: `wind(t, p)` is the inner field's wind at `p.y + offset_m`, with time and the
/// horizontal position unchanged.
///
/// `offset_m` is the height of the caller's y = 0 in the field's datum. Convert in one place
/// only: wrapping twice, or also shifting a profile's levels, moves every wind by the offset
/// again. The usual cases:
/// - y is height above sea level, as in rockets-engine: a [`LogProfile`] or [`EkmanProfile`],
///   which measure height from the ground, takes `-` the ground elevation. A sounding
///   [`ProfileField`] in metres above sea level needs no wrapper there.
/// - y = 0 is the ground at elevation `E`: a sea-level sounding [`ProfileField`] takes `+E`.
/// - this crate's [`TrajectorySolver`](crate::TrajectorySolver), whose y = 0 is the level
///   `muzzle_height` is measured from (the bullet starts at y = `muzzle_height`): a
///   [`LogProfile`] or [`EkmanProfile`] takes the bore's real height above the ground minus
///   `muzzle_height`, which is 0 when `muzzle_height` is that real height and the real height
///   when `muzzle_height` is left at 0, as `BallisticInputs::default()` leaves it. A sea-level
///   sounding there takes `+` the site's `altitude` (plus that same difference).
///
/// Elevations are those of the place the positions are, not of the station that measured a
/// sounding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeightOffset<F> {
    field: F,
    offset_m: f64,
}

impl<F> HeightOffset<F> {
    /// Wrap `field` so it is read `offset_m` metres higher. The offset must be finite.
    pub fn new(field: F, offset_m: f64) -> Result<Self, WindFieldError> {
        if !offset_m.is_finite() {
            return Err(WindFieldError::InvalidParameter {
                param: WindFieldParam::OffsetM,
                rule: WindFieldRule::Finite,
            });
        }
        Ok(Self { field, offset_m })
    }

    /// The wrapped field.
    pub fn field(&self) -> &F {
        &self.field
    }

    /// The height added to every query, metres.
    pub fn offset_m(&self) -> f64 {
        self.offset_m
    }

    /// Unwrap the field.
    pub fn into_inner(self) -> F {
        self.field
    }
}

impl<F: WindField> WindField for HeightOffset<F> {
    fn wind(&self, time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        self.field.wind(
            time_s,
            Vector3::new(position_m.x, position_m.y + self.offset_m, position_m.z),
        )
    }
}

/// One level of a vertical profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileLevel {
    /// Height of the level, metres (the y coordinate of the field's positions).
    pub altitude_m: f64,
    /// Horizontal wind speed, m/s (nonnegative).
    pub speed_mps: f64,
    /// Wind-FROM direction, degrees: 0° = headwind, 90° = from the right, matching
    /// [`WindSegment`](super::WindSegment).
    pub direction_from_deg: f64,
}

/// Piecewise-linear vertical profile: speed interpolates linearly in
/// altitude, direction along the shortest circular arc; nearest level holds
/// outside the range.
///
/// The field reads only the height (`position_m.y`), on the same datum as the levels'
/// `altitude_m`. A query exactly at a level gives that level's own vector, and exactly opposite
/// neighbouring directions turn clockwise (increasing FROM angle). Vertical wind is zero.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileField {
    /// Levels sorted by strictly ascending altitude.
    levels: Vec<ProfileLevel>,
}

impl ProfileField {
    /// Build a profile from levels in any order.
    ///
    /// Every value must be finite and every speed `>= 0`; these checks run before sorting, so a
    /// [`WindFieldError::InvalidLevel`] index refers to the caller's own ordering. After
    /// sorting, two levels at the same altitude are rejected with
    /// [`WindFieldError::DuplicateAltitude`].
    pub fn new(levels: Vec<ProfileLevel>) -> Result<Self, WindFieldError> {
        if levels.is_empty() {
            return Err(WindFieldError::EmptyProfile);
        }
        for (index, level) in levels.iter().enumerate() {
            let invalid = |param, rule| Err(WindFieldError::InvalidLevel { index, param, rule });
            if !level.altitude_m.is_finite() {
                return invalid(WindFieldParam::AltitudeM, WindFieldRule::Finite);
            }
            if !level.speed_mps.is_finite() || level.speed_mps < 0.0 {
                return invalid(
                    WindFieldParam::SpeedMps,
                    WindFieldRule::FiniteAndNonNegative,
                );
            }
            if !level.direction_from_deg.is_finite() {
                return invalid(WindFieldParam::DirectionFromDeg, WindFieldRule::Finite);
            }
        }

        let mut indexed: Vec<(usize, ProfileLevel)> = levels.into_iter().enumerate().collect();
        indexed.sort_by(|a, b| a.1.altitude_m.total_cmp(&b.1.altitude_m));
        for pair in indexed.windows(2) {
            let ((ia, a), (ib, b)) = (pair[0], pair[1]);
            // `==`, not total ordering: -0.0 and 0.0 are the same altitude.
            if a.altitude_m == b.altitude_m {
                return Err(WindFieldError::DuplicateAltitude {
                    first_index: ia.min(ib),
                    second_index: ia.max(ib),
                    altitude_m: a.altitude_m,
                });
            }
        }

        Ok(Self {
            levels: indexed.into_iter().map(|(_, level)| level).collect(),
        })
    }

    /// The levels, sorted by ascending altitude.
    pub fn levels(&self) -> &[ProfileLevel] {
        &self.levels
    }

    /// Wind (m/s) at an altitude. A NaN altitude gives zero wind, as
    /// [`WindSock::vector_for_range_stateless`] does for a NaN range.
    pub fn wind_at(&self, altitude_m: f64) -> Vector3<f64> {
        if altitude_m.is_nan() {
            return Vector3::zeros();
        }
        let levels = &self.levels;
        let first = levels[0];
        if altitude_m <= first.altitude_m {
            return horizontal(first.speed_mps, first.direction_from_deg);
        }
        for pair in levels.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            // Strictly below the upper level: a query exactly at a level is answered by that
            // level itself (the next segment's lower end, or the top level below), so every
            // level reproduces its own speed and direction exactly.
            if altitude_m < hi.altitude_m {
                if altitude_m == lo.altitude_m {
                    return horizontal(lo.speed_mps, lo.direction_from_deg);
                }
                // Halving keeps the span finite for levels near ±f64::MAX and is exact otherwise.
                let f = (altitude_m / 2.0 - lo.altitude_m / 2.0)
                    / (hi.altitude_m / 2.0 - lo.altitude_m / 2.0);
                let speed = lo.speed_mps + f * (hi.speed_mps - lo.speed_mps);
                // Interpolate between directions reduced to [0, 360), so a huge angle cannot
                // absorb the turn. Exactly opposite directions turn clockwise (the FROM angle
                // increases), however the angles are written.
                let lo_dir = lo.direction_from_deg.rem_euclid(360.0);
                let mut delta = hi.direction_from_deg.rem_euclid(360.0) - lo_dir;
                if delta > 180.0 {
                    delta -= 360.0;
                } else if delta <= -180.0 {
                    delta += 360.0;
                }
                return horizontal(speed, lo_dir + f * delta);
            }
        }
        let last = levels[levels.len() - 1];
        horizontal(last.speed_mps, last.direction_from_deg)
    }
}

impl WindField for ProfileField {
    fn wind(&self, _time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        self.wind_at(position_m.y)
    }
}

/// Log-law boundary layer: u(z) = u_ref · ln(z/z0) / ln(z_ref/z0), z > z0;
/// zero below z0. Constant direction with height.
///
/// The field reads only the height (`position_m.y`), which must be the height above the ground;
/// [`HeightOffset`] converts another frame's heights. Vertical wind is zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogProfile {
    reference_speed_mps: f64,
    reference_height_m: f64,
    roughness_m: f64,
    direction_from_deg: f64,
}

impl LogProfile {
    /// Build a log-law profile.
    ///
    /// `reference_speed_mps` must be finite and `>= 0`, `reference_height_m` finite and `> 0`,
    /// `roughness_m` (z0) finite, `> 0` and below the reference height, and
    /// `direction_from_deg` finite.
    pub fn new(
        reference_speed_mps: f64,
        reference_height_m: f64,
        roughness_m: f64,
        direction_from_deg: f64,
    ) -> Result<Self, WindFieldError> {
        let invalid = |param, rule| Err(WindFieldError::InvalidParameter { param, rule });
        if !reference_speed_mps.is_finite() || reference_speed_mps < 0.0 {
            return invalid(
                WindFieldParam::ReferenceSpeedMps,
                WindFieldRule::FiniteAndNonNegative,
            );
        }
        if !reference_height_m.is_finite() || reference_height_m <= 0.0 {
            return invalid(
                WindFieldParam::ReferenceHeightM,
                WindFieldRule::FiniteAndPositive,
            );
        }
        // The log difference is the profile's denominator; it must be a usable positive number,
        // which also rejects a roughness indistinguishable from the reference height.
        let log_span = reference_height_m.ln() - roughness_m.ln();
        if !roughness_m.is_finite()
            || roughness_m <= 0.0
            || roughness_m >= reference_height_m
            || !(log_span.is_finite() && log_span > 0.0)
        {
            return invalid(
                WindFieldParam::RoughnessM,
                WindFieldRule::FiniteAndPositiveBelowReferenceHeight,
            );
        }
        if !direction_from_deg.is_finite() {
            return invalid(WindFieldParam::DirectionFromDeg, WindFieldRule::Finite);
        }
        Ok(Self {
            reference_speed_mps,
            reference_height_m,
            roughness_m,
            direction_from_deg,
        })
    }

    /// Wind speed at the reference height, m/s.
    pub fn reference_speed_mps(&self) -> f64 {
        self.reference_speed_mps
    }

    /// Reference height, metres.
    pub fn reference_height_m(&self) -> f64 {
        self.reference_height_m
    }

    /// Aerodynamic roughness length z0, metres.
    pub fn roughness_m(&self) -> f64 {
        self.roughness_m
    }

    /// Wind-FROM direction, degrees (constant with height).
    pub fn direction_from_deg(&self) -> f64 {
        self.direction_from_deg
    }

    /// Wind speed (m/s) at an altitude: zero at or below the roughness length, and for a NaN or
    /// infinite altitude (the log law has no finite value there).
    pub fn speed_at(&self, altitude_m: f64) -> f64 {
        if !altitude_m.is_finite() || altitude_m <= self.roughness_m {
            return 0.0;
        }
        // Differences of logs rather than logs of ratios: z / z0 overflows for a tiny z0.
        let z0_ln = self.roughness_m.ln();
        self.reference_speed_mps * (altitude_m.ln() - z0_ln)
            / (self.reference_height_m.ln() - z0_ln)
    }

    /// Wind (m/s) at an altitude.
    pub fn wind_at(&self, altitude_m: f64) -> Vector3<f64> {
        horizontal(self.speed_at(altitude_m), self.direction_from_deg)
    }
}

impl WindField for LogProfile {
    fn wind(&self, _time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        self.wind_at(position_m.y)
    }
}

/// Ekman spiral (constant eddy viscosity K; Holton's classic solution):
/// u = G(1 − e^(−γz)·cos(γz)), v = sign(f)·G·e^(−γz)·sin(γz),
/// γ = sqrt(|f|/2K), f = 2Ω·sin(latitude). Zero at the surface, approaches
/// the geostrophic wind aloft, rotated to `direction_from_deg`.
///
/// `u` runs along the geostrophic flow and `v` to its left, so in the northern hemisphere the
/// wind near the surface is backed up to 45° from the geostrophic direction and veers with
/// height; the southern hemisphere mirrors this. This constant-K solution is linear near the
/// ground, not logarithmic; [`LogProfile`] is the surface-layer model.
///
/// The field reads only the height (`position_m.y`), which must be the height above the ground;
/// [`HeightOffset`] converts another frame's heights. An infinite height is the geostrophic
/// limit. Vertical wind is zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EkmanProfile {
    geostrophic_mps: f64,
    latitude_deg: f64,
    eddy_viscosity_m2s: f64,
    direction_from_deg: f64,
}

impl EkmanProfile {
    /// Build an Ekman profile.
    ///
    /// `geostrophic_mps` must be finite and `>= 0`, `latitude_deg` finite, nonzero and strictly
    /// between -90 and 90 (the solution has no boundary layer at the equator), and
    /// `eddy_viscosity_m2s` finite and `> 0`. `direction_from_deg`, the wind-FROM direction of
    /// the geostrophic wind aloft, must be finite. A latitude so close to zero that the Coriolis
    /// parameter underflows counts as the equator, and a viscosity so small that the decay rate
    /// overflows is rejected.
    pub fn new(
        geostrophic_mps: f64,
        latitude_deg: f64,
        eddy_viscosity_m2s: f64,
        direction_from_deg: f64,
    ) -> Result<Self, WindFieldError> {
        let invalid = |param, rule| Err(WindFieldError::InvalidParameter { param, rule });
        if !geostrophic_mps.is_finite() || geostrophic_mps < 0.0 {
            return invalid(
                WindFieldParam::GeostrophicMps,
                WindFieldRule::FiniteAndNonNegative,
            );
        }
        if !latitude_deg.is_finite() || latitude_deg == 0.0 || latitude_deg.abs() >= 90.0 {
            return invalid(
                WindFieldParam::LatitudeDeg,
                WindFieldRule::FiniteNonzeroAndOffThePoles,
            );
        }
        if coriolis_parameter(latitude_deg) == 0.0 {
            return invalid(
                WindFieldParam::LatitudeDeg,
                WindFieldRule::FiniteNonzeroAndOffThePoles,
            );
        }
        if !eddy_viscosity_m2s.is_finite()
            || eddy_viscosity_m2s <= 0.0
            || !decay_rate(latitude_deg, eddy_viscosity_m2s).is_finite()
        {
            return invalid(
                WindFieldParam::EddyViscosityM2s,
                WindFieldRule::FiniteAndPositive,
            );
        }
        if !direction_from_deg.is_finite() {
            return invalid(WindFieldParam::DirectionFromDeg, WindFieldRule::Finite);
        }
        Ok(Self {
            geostrophic_mps,
            latitude_deg,
            eddy_viscosity_m2s,
            direction_from_deg,
        })
    }

    /// Geostrophic wind speed aloft, m/s.
    pub fn geostrophic_mps(&self) -> f64 {
        self.geostrophic_mps
    }

    /// Latitude, degrees (its sign sets the turning sense).
    pub fn latitude_deg(&self) -> f64 {
        self.latitude_deg
    }

    /// Constant eddy viscosity K, m²/s.
    pub fn eddy_viscosity_m2s(&self) -> f64 {
        self.eddy_viscosity_m2s
    }

    /// Wind-FROM direction of the geostrophic wind aloft, degrees.
    pub fn direction_from_deg(&self) -> f64 {
        self.direction_from_deg
    }

    /// The along-flow and left-of-flow components `(u, v)`, m/s.
    fn components_at(&self, altitude_m: f64) -> (f64, f64) {
        if altitude_m.is_nan() || altitude_m <= 0.0 {
            return (0.0, 0.0);
        }
        let g = self.geostrophic_mps;
        let x = decay_rate(self.latitude_deg, self.eddy_viscosity_m2s) * altitude_m;
        // Past this, e^-x is exactly zero: the geostrophic limit, also for an infinite height
        // (where cos(x) would be NaN).
        if x > 746.0 {
            return (g, 0.0);
        }
        let decay = (-x).exp();
        let (sin_x, cos_x) = x.sin_cos();
        // 1 - e^-x cos x. Near the surface write it as 2 sin^2(x/2) - cos x * expm1(-x), which
        // has no cancellation, so the 45° surface turn survives at sub-millimetre heights.
        let along = if x < 1.0 {
            let half = (0.5 * x).sin();
            2.0 * half * half - cos_x * (-x).exp_m1()
        } else {
            1.0 - decay * cos_x
        };
        let sign = coriolis_parameter(self.latitude_deg).signum();
        (g * along, sign * g * decay * sin_x)
    }

    /// Wind (m/s) at an altitude: zero at or below the surface and for a NaN altitude, and the
    /// geostrophic wind for an infinite one.
    pub fn wind_at(&self, altitude_m: f64) -> Vector3<f64> {
        let (u, v) = self.components_at(altitude_m);
        let from_rad = direction_to_rad(self.direction_from_deg);
        let (sin, cos) = from_rad.sin_cos();
        // Along-flow unit vector: wind_vector(1, from_rad, 0) = (-cos, 0, -sin). Its left
        // (counterclockwise seen from above, with z to the right) is (-sin, 0, cos).
        Vector3::new(-u * cos - v * sin, 0.0, -u * sin + v * cos)
    }
}

impl WindField for EkmanProfile {
    fn wind(&self, _time_s: f64, position_m: Vector3<f64>) -> Vector3<f64> {
        self.wind_at(position_m.y)
    }
}

/// Which constructor input failed validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindFieldParam {
    /// [`ProfileLevel::altitude_m`].
    AltitudeM,
    /// [`ProfileLevel::speed_mps`].
    SpeedMps,
    /// A wind-FROM direction ([`ProfileLevel`], [`LogProfile`] or [`EkmanProfile`]).
    DirectionFromDeg,
    /// [`LogProfile`] reference speed.
    ReferenceSpeedMps,
    /// [`LogProfile`] reference height.
    ReferenceHeightM,
    /// [`LogProfile`] roughness length.
    RoughnessM,
    /// [`EkmanProfile`] geostrophic speed.
    GeostrophicMps,
    /// [`EkmanProfile`] latitude.
    LatitudeDeg,
    /// [`EkmanProfile`] eddy viscosity.
    EddyViscosityM2s,
    /// [`HeightOffset`] offset.
    OffsetM,
}

impl WindFieldParam {
    /// The parameter name, as it appears in error messages.
    pub fn name(self) -> &'static str {
        match self {
            WindFieldParam::AltitudeM => "altitude_m",
            WindFieldParam::SpeedMps => "speed_mps",
            WindFieldParam::DirectionFromDeg => "direction_from_deg",
            WindFieldParam::ReferenceSpeedMps => "reference_speed_mps",
            WindFieldParam::ReferenceHeightM => "reference_height_m",
            WindFieldParam::RoughnessM => "roughness_m",
            WindFieldParam::GeostrophicMps => "geostrophic_mps",
            WindFieldParam::LatitudeDeg => "latitude_deg",
            WindFieldParam::EddyViscosityM2s => "eddy_viscosity_m2s",
            WindFieldParam::OffsetM => "offset_m",
        }
    }
}

impl core::fmt::Display for WindFieldParam {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.name())
    }
}

/// The rule a constructor input violated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindFieldRule {
    /// The value must be finite.
    Finite,
    /// The value must be finite and `>= 0`.
    FiniteAndNonNegative,
    /// The value must be finite and `> 0`.
    FiniteAndPositive,
    /// The roughness length must be finite, `> 0` and below the reference height.
    FiniteAndPositiveBelowReferenceHeight,
    /// The latitude must be finite, nonzero and strictly between -90 and 90 degrees.
    FiniteNonzeroAndOffThePoles,
}

impl core::fmt::Display for WindFieldRule {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            WindFieldRule::Finite => "finite",
            WindFieldRule::FiniteAndNonNegative => "finite and non-negative",
            WindFieldRule::FiniteAndPositive => "finite and greater than zero",
            WindFieldRule::FiniteAndPositiveBelowReferenceHeight => {
                "finite, greater than zero and below the reference height"
            }
            WindFieldRule::FiniteNonzeroAndOffThePoles => {
                "finite, nonzero and strictly between -90 and 90 degrees"
            }
        })
    }
}

/// A wind-field constructor input that failed validation.
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum WindFieldError {
    /// A [`ProfileField`] needs at least one level.
    #[error("wind profile needs at least one level")]
    EmptyProfile,
    /// A [`ProfileLevel`] value is invalid. `index` is the level's position in the vector as
    /// supplied, before sorting.
    #[error("wind profile levels[{index}].{param} must be {rule}")]
    InvalidLevel {
        index: usize,
        param: WindFieldParam,
        rule: WindFieldRule,
    },
    /// Two [`ProfileLevel`]s share an altitude; altitudes must be strictly ascending once
    /// sorted. The indices are positions in the vector as supplied.
    #[error(
        "wind profile levels[{first_index}] and levels[{second_index}] share altitude_m {altitude_m}; altitudes must be strictly ascending"
    )]
    DuplicateAltitude {
        first_index: usize,
        second_index: usize,
        altitude_m: f64,
    },
    /// A [`LogProfile`], [`EkmanProfile`] or [`HeightOffset`] parameter is invalid.
    #[error("{param} must be {rule}")]
    InvalidParameter {
        param: WindFieldParam,
        rule: WindFieldRule,
    },
}

/// The Coriolis parameter f = 2Ω·sin(latitude), 1/s.
fn coriolis_parameter(latitude_deg: f64) -> f64 {
    2.0 * EARTH_ROTATION_RAD_PER_S * latitude_deg.to_radians().sin()
}

/// The Ekman decay rate γ = sqrt(|f| / 2K), 1/m.
fn decay_rate(latitude_deg: f64, eddy_viscosity_m2s: f64) -> f64 {
    (coriolis_parameter(latitude_deg).abs() / (2.0 * eddy_viscosity_m2s)).sqrt()
}

/// Degrees to radians with the same arithmetic as [`WindSegment`](super::WindSegment)'s
/// conversion, so equal angles give bit-identical vectors.
fn direction_to_rad(direction_deg: f64) -> f64 {
    if direction_deg.abs() <= f64::MAX / PI {
        direction_deg * PI / 180.0
    } else {
        direction_deg / 180.0 * PI
    }
}

/// Horizontal wind from a speed and a wind-FROM direction in degrees.
fn horizontal(speed_mps: f64, direction_from_deg: f64) -> Vector3<f64> {
    wind_vector(speed_mps, direction_to_rad(direction_from_deg), 0.0)
}

#[cfg(test)]
mod tests {
    use super::super::WindSegment;
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    fn level(altitude_m: f64, speed_mps: f64, direction_from_deg: f64) -> ProfileLevel {
        ProfileLevel {
            altitude_m,
            speed_mps,
            direction_from_deg,
        }
    }

    fn at_height(y: f64) -> Vector3<f64> {
        Vector3::new(0.0, y, 0.0)
    }

    /// Horizontal speed of a wind vector.
    fn speed(w: Vector3<f64>) -> f64 {
        w.x.hypot(w.z)
    }

    /// Wind-FROM direction of a wind vector, degrees in [0, 360).
    fn direction_from(w: Vector3<f64>) -> f64 {
        (w.z.atan2(w.x).to_degrees() + 180.0).rem_euclid(360.0)
    }

    /// Signed difference a - b, wrapped to (-180, 180].
    fn angle_diff(a: f64, b: f64) -> f64 {
        let d = (a - b).rem_euclid(360.0);
        if d > 180.0 {
            d - 360.0
        } else {
            d
        }
    }

    // ---- Ported from winds-engine tests/profile.rs ----

    #[test]
    fn wind_direction_round_trips() {
        // FROM 0° blows toward −X (a headwind); FROM 90° toward −Z.
        let headwind = horizontal(10.0, 0.0);
        assert!((headwind.x + 10.0).abs() < 1e-12);
        assert!(headwind.z.abs() < 1e-12);
        assert!(direction_from(headwind).abs() < 1e-9);
        let cross = horizontal(10.0, 90.0);
        assert!(cross.x.abs() < 1e-12);
        assert!((cross.z + 10.0).abs() < 1e-12);
        for deg in [1.0, 90.0, 179.9, 181.0, 270.0, 359.9] {
            let w = horizontal(7.0, deg);
            assert!((speed(w) - 7.0).abs() < 1e-9);
            let back = angle_diff(direction_from(w), deg);
            assert!(back.abs() < 1e-9, "deg {deg} back {back}");
        }
    }

    #[test]
    fn constant_field_is_constant() {
        let field = ConstantField(horizontal(5.0, 45.0));
        assert_eq!(
            field.wind(0.0, Vector3::zeros()),
            field.wind(99.0, Vector3::new(1.0, 2.0, 3.0))
        );
    }

    #[test]
    fn profile_interpolates_speed_linearly() {
        let profile =
            ProfileField::new(vec![level(0.0, 2.0, 0.0), level(1000.0, 12.0, 0.0)]).unwrap();
        assert!((speed(profile.wind_at(500.0)) - 7.0).abs() < 1e-12);
        // Below and above the levels, the nearest holds.
        assert!((speed(profile.wind_at(-50.0)) - 2.0).abs() < 1e-12);
        assert!((speed(profile.wind_at(5000.0)) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn profile_interpolates_direction_along_the_short_arc() {
        // 350° to 10°: through 0°, never through 180°.
        let profile =
            ProfileField::new(vec![level(0.0, 5.0, 350.0), level(1000.0, 5.0, 10.0)]).unwrap();
        let mid = direction_from(profile.wind_at(500.0));
        assert!(!(10.0..=350.0).contains(&mid), "mid {mid}");
        let quarter = direction_from(profile.wind_at(250.0));
        assert!(quarter > 350.0, "quarter {quarter}");
        // And the reverse arc, 10° down to 350°.
        let reverse =
            ProfileField::new(vec![level(0.0, 5.0, 10.0), level(1000.0, 5.0, 350.0)]).unwrap();
        assert!(angle_diff(direction_from(reverse.wind_at(500.0)), 0.0).abs() < 1e-9);
    }

    #[test]
    fn profile_validates_its_levels() {
        assert!(ProfileField::new(vec![]).is_err());
        assert!(ProfileField::new(vec![level(0.0, -1.0, 0.0)]).is_err());
        assert!(ProfileField::new(vec![level(0.0, 1.0, 0.0), level(0.0, 2.0, 0.0)]).is_err());
        // Unsorted input is accepted and sorted.
        let profile =
            ProfileField::new(vec![level(1000.0, 12.0, 0.0), level(0.0, 2.0, 0.0)]).unwrap();
        assert!((speed(profile.wind_at(500.0)) - 7.0).abs() < 1e-12);
        assert_eq!(profile.levels()[0].altitude_m, 0.0);
        assert_eq!(profile.levels()[1].altitude_m, 1000.0);
    }

    #[test]
    fn log_profile_hits_its_reference_and_roughness() {
        let log = LogProfile::new(8.0, 10.0, 0.03, 270.0).unwrap();
        assert!(
            (log.speed_at(10.0) - 8.0).abs() < 1e-12,
            "{}",
            log.speed_at(10.0)
        );
        assert_eq!(log.speed_at(0.03), 0.0);
        assert_eq!(log.speed_at(0.01), 0.0);
        // Doubling height adds a constant increment in a log law.
        let increment = log.speed_at(20.0) - log.speed_at(10.0);
        assert!((increment - (log.speed_at(40.0) - log.speed_at(20.0))).abs() < 1e-12);
        assert!(LogProfile::new(8.0, 10.0, 20.0, 0.0).is_err());
        assert!(LogProfile::new(-1.0, 10.0, 0.03, 0.0).is_err());
    }

    #[test]
    fn ekman_profile_has_the_right_endpoints() {
        let ekman = EkmanProfile::new(12.0, 45.0, 5.0, 0.0).unwrap();
        // Calm at the surface, geostrophic aloft.
        assert!(speed(ekman.wind_at(0.0)) < 1e-12);
        let aloft = ekman.wind_at(5_000.0);
        assert!(
            (speed(aloft) - 12.0).abs() / 12.0 < 1e-3,
            "{}",
            speed(aloft)
        );
        // Aloft it comes FROM the requested geostrophic direction.
        let dir_err = angle_diff(direction_from(aloft), 0.0);
        assert!(dir_err.abs() < 1.0, "dir_err {dir_err}");
        // Near the ground the speed is between zero and geostrophic.
        let low = ekman.wind_at(100.0);
        assert!(speed(low) > 0.0 && speed(low) < 12.0);
        assert!(EkmanProfile::new(12.0, 90.0, 5.0, 0.0).is_err());
        assert!(EkmanProfile::new(12.0, 45.0, 0.0, 0.0).is_err());
    }

    // ---- Conventions ----

    #[test]
    fn directions_match_wind_vector_exactly() {
        for speed_mps in [0.0, 3.0, 12.5] {
            let at0 = ProfileField::new(vec![level(0.0, speed_mps, 0.0)]).unwrap();
            assert_eq!(at0.wind_at(10.0), wind_vector(speed_mps, 0.0, 0.0));
            let at90 = ProfileField::new(vec![level(0.0, speed_mps, 90.0)]).unwrap();
            assert_eq!(at90.wind_at(10.0), wind_vector(speed_mps, FRAC_PI_2, 0.0));
        }
        // A level and a WindSegment with the same speed and angle give the same vector.
        let segment_speed_kmh = 36.0;
        let sock = WindSock::try_new(vec![WindSegment::new(segment_speed_kmh, 33.0, 1e6)]).unwrap();
        let from_sock = sock.vector_for_range_stateless(1.0);
        let profile = ProfileField::new(vec![level(
            0.0,
            segment_speed_kmh * (1000.0 / 3600.0),
            33.0,
        )])
        .unwrap();
        assert_eq!(profile.wind_at(0.0), from_sock);
    }

    #[test]
    fn analytic_profiles_use_the_same_direction_convention() {
        // Log: constant direction, so a 90° profile blows toward -z at every height.
        let log = LogProfile::new(8.0, 10.0, 0.03, 90.0).unwrap();
        let w = log.wind_at(50.0);
        assert!(w.x.abs() < 1e-12 && w.z < 0.0 && w.y == 0.0);
        assert!((direction_from(w) - 90.0).abs() < 1e-9);
        // Ekman aloft: FROM the requested direction, for several directions.
        for dir in [0.0, 90.0, 180.0, 270.0, 33.0] {
            let ekman = EkmanProfile::new(10.0, 45.0, 5.0, dir).unwrap();
            let aloft = ekman.wind_at(20_000.0);
            assert!(angle_diff(direction_from(aloft), dir).abs() < 1e-6, "{dir}");
            assert_eq!(aloft.y, 0.0);
        }
    }

    #[test]
    fn ekman_backs_near_the_surface_in_the_north_and_veers_in_the_south() {
        // Holton: in the northern hemisphere the surface wind is turned toward low pressure,
        // 45° to the LEFT of (backed from) the geostrophic wind, and the wind veers with
        // height. "Left" here is counterclockwise seen from above with z to the right, which
        // lowers the wind-FROM angle. The southern hemisphere mirrors it.
        let geo_dir = 200.0;
        let north = EkmanProfile::new(10.0, 45.0, 5.0, geo_dir).unwrap();
        let south = EkmanProfile::new(10.0, -45.0, 5.0, geo_dir).unwrap();
        let near = 0.01; // metres: gamma * z << 1, so u ≈ v
        assert!((angle_diff(direction_from(north.wind_at(near)), geo_dir) + 45.0).abs() < 0.01);
        assert!((angle_diff(direction_from(south.wind_at(near)), geo_dir) - 45.0).abs() < 0.01);
        // Veering with height in the north: the FROM angle rises from about geo - 45°.
        let mut previous = angle_diff(direction_from(north.wind_at(1.0)), geo_dir);
        for z in [50.0, 100.0, 200.0, 400.0] {
            let now = angle_diff(direction_from(north.wind_at(z)), geo_dir);
            assert!(now > previous, "z {z}: {now} <= {previous}");
            previous = now;
        }
    }

    // ---- Adapter ----

    #[test]
    fn wind_sock_adapter_matches_the_stateless_lookup() {
        let sock = WindSock::try_new(vec![
            WindSegment::new(10.0, 90.0, 100.0),
            WindSegment {
                speed_kmh: 20.0,
                angle_deg: 0.0,
                until_m: 300.0,
                vertical_mps: 1.5,
            },
            WindSegment::new(5.0, 270.0, 600.0),
        ])
        .unwrap();
        for x in [
            -50.0,
            0.0,
            99.999,
            100.0,
            100.001,
            250.0,
            299.999,
            300.0,
            450.0,
            599.0,
            600.0,
            1e6,
            f64::NAN,
        ] {
            for (t, y, z) in [(0.0, 0.0, 0.0), (12.5, 40.0, -3.0), (-1.0, -5.0, 7.0)] {
                let got = sock.wind(t, Vector3::new(x, y, z));
                let want = sock.vector_for_range_stateless(x);
                assert_eq!(got, want, "x {x} t {t} y {y} z {z}");
            }
        }
        // Through a trait object too.
        let field: &dyn WindField = &sock;
        assert_eq!(
            field.wind(0.0, Vector3::new(150.0, 0.0, 0.0)),
            sock.vector_for_range_stateless(150.0)
        );
    }

    // ---- Validation ----

    #[test]
    fn profile_validation_reports_index_param_and_rule() {
        let ok = level(0.0, 1.0, 0.0);
        let cases = [
            (
                level(f64::NAN, 1.0, 0.0),
                WindFieldParam::AltitudeM,
                WindFieldRule::Finite,
            ),
            (
                level(f64::INFINITY, 1.0, 0.0),
                WindFieldParam::AltitudeM,
                WindFieldRule::Finite,
            ),
            (
                level(0.0, -0.5, 0.0),
                WindFieldParam::SpeedMps,
                WindFieldRule::FiniteAndNonNegative,
            ),
            (
                level(0.0, f64::INFINITY, 0.0),
                WindFieldParam::SpeedMps,
                WindFieldRule::FiniteAndNonNegative,
            ),
            (
                level(0.0, f64::NAN, 0.0),
                WindFieldParam::SpeedMps,
                WindFieldRule::FiniteAndNonNegative,
            ),
            (
                level(0.0, 1.0, f64::NAN),
                WindFieldParam::DirectionFromDeg,
                WindFieldRule::Finite,
            ),
            (
                level(0.0, 1.0, f64::NEG_INFINITY),
                WindFieldParam::DirectionFromDeg,
                WindFieldRule::Finite,
            ),
        ];
        for (bad, param, rule) in cases {
            // The bad level is second in the caller's order but would sort first.
            let mut good = ok;
            good.altitude_m = 500.0;
            let err = ProfileField::new(vec![good, bad]).unwrap_err();
            assert_eq!(
                err,
                WindFieldError::InvalidLevel {
                    index: 1,
                    param,
                    rule
                }
            );
        }
        assert_eq!(
            ProfileField::new(vec![]).unwrap_err(),
            WindFieldError::EmptyProfile
        );
        // Unsorted duplicates are found after sorting and reported in caller order.
        let err = ProfileField::new(vec![
            level(300.0, 1.0, 0.0),
            level(100.0, 1.0, 0.0),
            level(200.0, 1.0, 0.0),
            level(100.0, 2.0, 0.0),
        ])
        .unwrap_err();
        assert_eq!(
            err,
            WindFieldError::DuplicateAltitude {
                first_index: 1,
                second_index: 3,
                altitude_m: 100.0
            }
        );
        // -0.0 and 0.0 are the same altitude.
        assert!(matches!(
            ProfileField::new(vec![level(0.0, 1.0, 0.0), level(-0.0, 1.0, 0.0)]),
            Err(WindFieldError::DuplicateAltitude { .. })
        ));
        // Messages follow the module's style.
        assert_eq!(
            WindFieldError::InvalidLevel {
                index: 2,
                param: WindFieldParam::SpeedMps,
                rule: WindFieldRule::FiniteAndNonNegative
            }
            .to_string(),
            "wind profile levels[2].speed_mps must be finite and non-negative"
        );
    }

    #[test]
    fn analytic_validation_reports_param_and_rule() {
        use WindFieldParam as P;
        use WindFieldRule as R;
        let log = |s, h, z0, d| LogProfile::new(s, h, z0, d).unwrap_err();
        let invalid = |param, rule| WindFieldError::InvalidParameter { param, rule };
        assert_eq!(
            log(-1.0, 10.0, 0.03, 0.0),
            invalid(P::ReferenceSpeedMps, R::FiniteAndNonNegative)
        );
        assert_eq!(
            log(f64::NAN, 10.0, 0.03, 0.0),
            invalid(P::ReferenceSpeedMps, R::FiniteAndNonNegative)
        );
        assert_eq!(
            log(8.0, 0.0, 0.03, 0.0),
            invalid(P::ReferenceHeightM, R::FiniteAndPositive)
        );
        assert_eq!(
            log(8.0, f64::INFINITY, 0.03, 0.0),
            invalid(P::ReferenceHeightM, R::FiniteAndPositive)
        );
        for z0 in [0.0, -0.1, 10.0, 11.0, f64::NAN] {
            assert_eq!(
                log(8.0, 10.0, z0, 0.0),
                invalid(P::RoughnessM, R::FiniteAndPositiveBelowReferenceHeight),
                "z0 {z0}"
            );
        }
        assert_eq!(
            log(8.0, 10.0, 0.03, f64::NAN),
            invalid(P::DirectionFromDeg, R::Finite)
        );

        let ekman = |g, lat, k, d| EkmanProfile::new(g, lat, k, d).unwrap_err();
        assert_eq!(
            ekman(-1.0, 45.0, 5.0, 0.0),
            invalid(P::GeostrophicMps, R::FiniteAndNonNegative)
        );
        for lat in [0.0, -0.0, 90.0, -90.0, 120.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                ekman(12.0, lat, 5.0, 0.0),
                invalid(P::LatitudeDeg, R::FiniteNonzeroAndOffThePoles),
                "lat {lat}"
            );
        }
        for k in [0.0, -5.0, f64::NAN] {
            assert_eq!(
                ekman(12.0, 45.0, k, 0.0),
                invalid(P::EddyViscosityM2s, R::FiniteAndPositive)
            );
        }
        assert_eq!(
            ekman(12.0, 45.0, 5.0, f64::INFINITY),
            invalid(P::DirectionFromDeg, R::Finite)
        );
        // Zero speeds are valid and give calm air.
        assert_eq!(
            LogProfile::new(0.0, 10.0, 0.03, 0.0).unwrap().wind_at(5.0),
            Vector3::zeros()
        );
        assert_eq!(
            speed(
                EkmanProfile::new(0.0, 45.0, 5.0, 0.0)
                    .unwrap()
                    .wind_at(100.0)
            ),
            0.0
        );
    }

    // ---- Determinism and edge positions ----

    #[test]
    fn fields_ignore_time_and_horizontal_position() {
        let profile =
            ProfileField::new(vec![level(0.0, 2.0, 30.0), level(800.0, 9.0, 80.0)]).unwrap();
        let log = LogProfile::new(8.0, 10.0, 0.03, 270.0).unwrap();
        let ekman = EkmanProfile::new(12.0, 45.0, 5.0, 120.0).unwrap();
        let fields: [&dyn WindField; 3] = [&profile, &log, &ekman];
        for field in fields {
            let reference = field.wind(0.0, at_height(321.0));
            for (t, x, z) in [(1.0, 10.0, -4.0), (1e6, -300.0, 99.0), (-7.0, 0.0, 0.0)] {
                assert_eq!(field.wind(t, Vector3::new(x, 321.0, z)), reference);
            }
            // Repeated calls are bit-identical.
            assert_eq!(field.wind(3.0, at_height(321.0)), reference);
        }
    }

    #[test]
    fn nan_height_gives_calm_air() {
        let profile = ProfileField::new(vec![level(0.0, 2.0, 30.0)]).unwrap();
        let log = LogProfile::new(8.0, 10.0, 0.03, 270.0).unwrap();
        let ekman = EkmanProfile::new(12.0, 45.0, 5.0, 120.0).unwrap();
        assert_eq!(profile.wind_at(f64::NAN), Vector3::zeros());
        assert_eq!(log.wind_at(f64::NAN), Vector3::zeros());
        assert_eq!(log.speed_at(f64::NAN), 0.0);
        assert_eq!(ekman.wind_at(f64::NAN), Vector3::zeros());
        // Below the surface the analytic profiles are calm too.
        assert_eq!(ekman.wind_at(-10.0), Vector3::zeros());
        assert_eq!(log.wind_at(-10.0), Vector3::zeros());
    }

    #[test]
    fn single_level_profile_is_uniform_with_height() {
        let profile = ProfileField::new(vec![level(150.0, 4.0, 45.0)]).unwrap();
        let w = profile.wind_at(150.0);
        for y in [-1e3, 0.0, 149.0, 151.0, 1e5] {
            assert_eq!(profile.wind_at(y), w);
        }
    }

    // ---- Review fixes: exact levels, extremes, numerics ----

    #[test]
    fn every_level_reproduces_its_own_vector_exactly() {
        // Bottom, interior (after a 350 -> 10 wrap) and top levels, plus a raw angle >= 360.
        let levels = vec![
            level(0.0, 3.0, 350.0),
            level(100.0, 5.0, 10.0),
            level(250.0, 8.0, 100.0),
            level(330.0, 6.5, 370.0),
            level(400.0, 2.0, 359.0),
        ];
        let profile = ProfileField::new(levels.clone()).unwrap();
        for l in levels {
            assert_eq!(
                profile.wind_at(l.altitude_m),
                wind_vector(l.speed_mps, direction_to_rad(l.direction_from_deg), 0.0),
                "level at {} m",
                l.altitude_m
            );
        }
    }

    #[test]
    fn opposite_directions_turn_clockwise_however_written() {
        // Exactly opposite neighbours: the FROM angle increases by 90° at the midpoint.
        for (lo, hi) in [
            (0.0, 180.0),
            (0.0, -180.0),
            (180.0, 0.0),
            (180.0, 360.0),
            (-180.0, 0.0),
            (190.0, 10.0),
            (190.0, 370.0),
        ] {
            let profile =
                ProfileField::new(vec![level(0.0, 5.0, lo), level(100.0, 5.0, hi)]).unwrap();
            let expected = (lo + 90.0_f64).rem_euclid(360.0);
            let mid = direction_from(profile.wind_at(50.0));
            assert!(
                angle_diff(mid, expected).abs() < 1e-9,
                "{lo} -> {hi}: {mid}"
            );
        }
    }

    #[test]
    fn huge_level_angles_still_interpolate() {
        let lo_dir = 1e300_f64;
        let profile =
            ProfileField::new(vec![level(0.0, 5.0, lo_dir), level(100.0, 5.0, 10.0)]).unwrap();
        let lo_n = lo_dir.rem_euclid(360.0);
        let mut delta = 10.0 - lo_n;
        if delta > 180.0 {
            delta -= 360.0;
        } else if delta <= -180.0 {
            delta += 360.0;
        }
        for frac in [0.25, 0.5, 0.75] {
            let want = (lo_n + frac * delta).rem_euclid(360.0);
            let got = direction_from(profile.wind_at(100.0 * frac));
            assert!(
                angle_diff(got, want).abs() < 1e-6,
                "frac {frac}: {got} vs {want}"
            );
        }
    }

    #[test]
    fn extreme_altitude_spans_stay_finite() {
        let profile =
            ProfileField::new(vec![level(-1e308, 2.0, 0.0), level(1e308, 6.0, 0.0)]).unwrap();
        assert!((speed(profile.wind_at(0.0)) - 4.0).abs() < 1e-12);
        assert_eq!(profile.wind_at(1e308), horizontal(6.0, 0.0));
        assert_eq!(profile.wind_at(f64::INFINITY), horizontal(6.0, 0.0));
        assert_eq!(profile.wind_at(f64::NEG_INFINITY), horizontal(2.0, 0.0));
    }

    #[test]
    fn log_profile_is_finite_for_every_accepted_roughness() {
        for z0 in [1e-310, 1e-300, 1e-12, 0.03, 9.999] {
            let log = LogProfile::new(8.0, 10.0, z0, 30.0).unwrap();
            assert!((log.speed_at(10.0) - 8.0).abs() < 1e-12, "z0 {z0}");
            for z in [z0 * 2.0, 1.0, 1e9, f64::MAX] {
                let w = log.wind_at(z);
                assert!(w.iter().all(|c| c.is_finite()), "z0 {z0} z {z}: {w:?}");
            }
        }
        // A roughness indistinguishable from the reference height on a log scale is refused;
        // anything accepted must still hit its reference speed.
        let z0 = f64::from_bits(10.0_f64.to_bits() - 1);
        match LogProfile::new(8.0, 10.0, z0, 0.0) {
            Err(e) => assert_eq!(
                e,
                WindFieldError::InvalidParameter {
                    param: WindFieldParam::RoughnessM,
                    rule: WindFieldRule::FiniteAndPositiveBelowReferenceHeight
                }
            ),
            Ok(log) => assert!((log.speed_at(10.0) - 8.0).abs() < 1e-9),
        }
        // No finite value at an infinite height: calm, like NaN.
        let log = LogProfile::new(8.0, 10.0, 0.03, 30.0).unwrap();
        assert_eq!(log.wind_at(f64::INFINITY), Vector3::zeros());
        assert_eq!(log.wind_at(f64::NEG_INFINITY), Vector3::zeros());
    }

    #[test]
    fn ekman_rejects_underflowing_latitude_and_overflowing_decay() {
        let invalid = |param, rule| WindFieldError::InvalidParameter { param, rule };
        assert_eq!(
            EkmanProfile::new(12.0, 1e-322, 5.0, 0.0).unwrap_err(),
            invalid(
                WindFieldParam::LatitudeDeg,
                WindFieldRule::FiniteNonzeroAndOffThePoles
            )
        );
        assert_eq!(
            EkmanProfile::new(12.0, 45.0, 1e-320, 0.0).unwrap_err(),
            invalid(
                WindFieldParam::EddyViscosityM2s,
                WindFieldRule::FiniteAndPositive
            )
        );
        // Small but representable values are still accepted.
        assert!(EkmanProfile::new(12.0, 1e-6, 5.0, 0.0).is_ok());
        assert!(EkmanProfile::new(12.0, 45.0, 1e-300, 0.0).is_ok());
    }

    #[test]
    fn ekman_is_finite_at_every_height_and_geostrophic_at_infinity() {
        for k in [1e-300, 5.0, 1e6] {
            let ekman = EkmanProfile::new(12.0, -30.0, k, 75.0).unwrap();
            for z in [1e-300, 1e-15, 1.0, 1e4, 1e200, f64::MAX, f64::INFINITY] {
                let w = ekman.wind_at(z);
                assert!(w.iter().all(|c| c.is_finite()), "k {k} z {z}: {w:?}");
            }
            let top = ekman.wind_at(f64::INFINITY);
            assert!((speed(top) - 12.0).abs() < 1e-12);
            assert!(angle_diff(direction_from(top), 75.0).abs() < 1e-9);
        }
    }

    #[test]
    fn ekman_surface_turn_holds_at_tiny_heights() {
        let north = EkmanProfile::new(10.0, 45.0, 5.0, 0.0).unwrap();
        let south = EkmanProfile::new(10.0, -45.0, 5.0, 0.0).unwrap();
        for z in [1e-15, 1e-13, 1e-9, 1e-5] {
            let n = angle_diff(direction_from(north.wind_at(z)), 0.0);
            let s = angle_diff(direction_from(south.wind_at(z)), 0.0);
            assert!((n + 45.0).abs() < 1e-3, "north z {z}: {n}");
            assert!((s - 45.0).abs() < 1e-3, "south z {z}: {s}");
        }
    }

    // ---- Height datum ----

    #[test]
    fn height_offset_reads_the_inner_field_higher() {
        let profile =
            ProfileField::new(vec![level(0.0, 2.0, 30.0), level(800.0, 9.0, 80.0)]).unwrap();
        let shifted = HeightOffset::new(&profile, 250.0).unwrap();
        for (t, x, y, z) in [
            (0.0, 0.0, 0.0, 0.0),
            (3.5, 120.0, 40.0, -2.0),
            (1.0, -5.0, -300.0, 9.0),
        ] {
            assert_eq!(
                shifted.wind(t, Vector3::new(x, y, z)),
                profile.wind(t, Vector3::new(x, y + 250.0, z))
            );
        }
        assert_eq!(shifted.offset_m(), 250.0);
        assert_eq!(*shifted.field(), &profile);
        // A zero offset is the field itself.
        let same = HeightOffset::new(&profile, 0.0).unwrap();
        assert_eq!(same.wind(0.0, at_height(123.0)), profile.wind_at(123.0));
    }

    #[test]
    fn a_sea_level_sounding_matches_a_ground_frame_through_the_site_elevation() {
        // Levels in metres above sea level, the first at the station's 287 m surface. A launch
        // frame with y = 0 on the ground at that elevation reads the surface level at y = 0 and
        // the 610 m level 323 m up.
        let msl = ProfileField::new(vec![
            level(287.0, 0.5, 120.0),
            level(305.0, 3.0867, 125.0),
            level(610.0, 3.6011, 145.0),
        ])
        .unwrap();
        let site = HeightOffset::new(msl.clone(), 287.0).unwrap();
        assert_eq!(site.wind(0.0, Vector3::zeros()), horizontal(0.5, 120.0));
        assert_eq!(site.wind(0.0, at_height(323.0)), horizontal(3.6011, 145.0));
        // Without the offset the ground frame would read the surface value up to 287 m.
        assert_eq!(
            msl.wind(0.0, at_height(287.0 - 1.0)),
            horizontal(0.5, 120.0)
        );
        assert_ne!(msl.wind(0.0, at_height(323.0)), horizontal(3.6011, 145.0));
    }

    #[test]
    fn a_frame_with_y_zero_at_the_bore_reads_a_log_profile_at_the_bore_height() {
        let log = LogProfile::new(8.0, 10.0, 0.03, 90.0).unwrap();
        let from_bore = HeightOffset::new(log, 1.5).unwrap();
        assert_eq!(from_bore.wind(0.0, Vector3::zeros()), log.wind_at(1.5));
        assert!(speed(from_bore.wind(0.0, Vector3::zeros())) > 0.0);
        assert_eq!(from_bore.into_inner(), log);
    }

    #[test]
    fn a_sea_level_frame_reads_a_log_profile_through_minus_the_ground_elevation() {
        // y is height above sea level and the ground sits 1609 m up: the log profile, measured
        // from the ground, takes -1609; a sea-level sounding profile is read as it is.
        let ground_msl = 1609.0;
        let log = LogProfile::new(8.0, 10.0, 0.03, 270.0).unwrap();
        let in_msl_frame = HeightOffset::new(log, -ground_msl).unwrap();
        assert_eq!(
            in_msl_frame.wind(0.0, at_height(ground_msl + 10.0)),
            log.wind_at(10.0)
        );
        assert_eq!(
            in_msl_frame.wind(0.0, at_height(ground_msl)),
            Vector3::zeros()
        );
        assert!((speed(in_msl_frame.wind(0.0, at_height(ground_msl + 10.0))) - 8.0).abs() < 1e-12);

        let sounding_msl =
            ProfileField::new(vec![level(1609.0, 1.0, 200.0), level(2609.0, 9.0, 250.0)]).unwrap();
        assert_eq!(
            sounding_msl.wind(0.0, at_height(2109.0)),
            sounding_msl.wind_at(2109.0)
        );
        // The double-shift trap: wrapping a profile that is already in the frame's datum moves
        // every wind by the offset.
        let double = HeightOffset::new(&sounding_msl, ground_msl).unwrap();
        assert_ne!(
            double.wind(0.0, at_height(2109.0)),
            sounding_msl.wind_at(2109.0)
        );
    }

    #[test]
    fn height_offset_must_be_finite() {
        for offset in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                HeightOffset::new(ConstantField(Vector3::zeros()), offset).unwrap_err(),
                WindFieldError::InvalidParameter {
                    param: WindFieldParam::OffsetM,
                    rule: WindFieldRule::Finite
                }
            );
        }
        assert_eq!(
            WindFieldError::InvalidParameter {
                param: WindFieldParam::OffsetM,
                rule: WindFieldRule::Finite
            }
            .to_string(),
            "offset_m must be finite"
        );
    }

    #[test]
    fn boxed_and_borrowed_fields_are_fields() {
        let boxed: Vec<Box<dyn WindField>> = vec![
            Box::new(ConstantField(Vector3::new(1.0, 2.0, 3.0))),
            Box::new(LogProfile::new(8.0, 10.0, 0.03, 0.0).unwrap()),
        ];
        let p = at_height(25.0);
        assert_eq!(boxed[0].wind(0.0, p), Vector3::new(1.0, 2.0, 3.0));
        fn evaluate<F: WindField>(field: F, p: Vector3<f64>) -> Vector3<f64> {
            field.wind(0.0, p)
        }
        assert_eq!(evaluate(&boxed[1], p), boxed[1].wind(0.0, p));
        assert_eq!(
            evaluate(boxed.into_iter().next().unwrap(), p),
            Vector3::new(1.0, 2.0, 3.0)
        );
    }
}
