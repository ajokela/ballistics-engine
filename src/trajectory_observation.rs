//! Checked, full-state observations over completed trajectories.
//!
//! The legacy [`crate::TrajectoryResult::position_at_range`] API intentionally clamps queries
//! beyond the computed trajectory.  The APIs in this module are for protocol and laboratory
//! consumers that need explicit range errors, finite values, and an exact terminal sample.

use crate::cli_api::{TrajectoryPoint, TrajectoryResult};
use crate::trajectory_sampling::MAX_TRAJECTORY_SAMPLES;
use thiserror::Error;

/// Why trajectory integration stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TrajectoryTermination {
    MaxRange,
    GroundThreshold,
    TimeLimit,
    VelocityFloor,
}

/// Stable annotations attached to a checked trajectory observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TrajectoryObservationFlag {
    /// Mach is in the conventional inclusive transonic band, `0.8..=1.2`.
    Transonic,
    /// Mach is below `1.0`; this intentionally overlaps `Transonic` from `0.8` to `1.0`.
    Subsonic,
    /// This is the exact final observation returned by the solver.
    Terminal,
    /// The terminal observation is an impact at the configured ground threshold.
    GroundThreshold,
}

/// A finite, sight-referenced observation of the complete trajectory state.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryObservation {
    pub distance_m: f64,
    pub time_s: f64,
    pub speed_mps: f64,
    pub energy_j: f64,
    /// Positive values are below the line of sight.
    pub drop_m: f64,
    /// Positive values are right of the line of sight from the shooter's perspective.
    pub windage_m: f64,
    pub mach: f64,
    pub flags: Vec<TrajectoryObservationFlag>,
}

/// Failure to produce a checked trajectory observation.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum TrajectoryObservationError {
    #[error("trajectory contains no points")]
    EmptyTrajectory,

    #[error("observation range must be finite (got {distance_m})")]
    NonFiniteQuery { distance_m: f64 },

    #[error("sample interval must be finite and greater than zero (got {interval_m})")]
    InvalidInterval { interval_m: f64 },

    #[error(
        "requested range {requested_m} m is outside the computed trajectory [{minimum_m}, {maximum_m}] m"
    )]
    OutOfRange {
        requested_m: f64,
        minimum_m: f64,
        maximum_m: f64,
    },

    #[error(
        "trajectory distance is not strictly increasing at point {index}: {previous_distance_m} m then {distance_m} m"
    )]
    NonMonotonicTrajectory {
        index: usize,
        previous_distance_m: f64,
        distance_m: f64,
    },

    #[error("trajectory point {index} contains a non-finite {field}")]
    NonFiniteState { index: usize, field: &'static str },

    #[error("trajectory point {index} contains an invalid {field} value ({value})")]
    InvalidState {
        index: usize,
        field: &'static str,
        value: f64,
    },

    #[error("trajectory metadata field {field} has an invalid value ({value})")]
    InvalidMetadata { field: &'static str, value: f64 },

    #[error("computed observation field {field} is not finite")]
    NonFiniteObservation { field: &'static str },

    #[error("trajectory observation limit of {limit} exceeded (would produce {requested})")]
    SampleLimitExceeded { requested: usize, limit: usize },

    #[error("could not reserve storage for {requested} trajectory observations")]
    AllocationFailed { requested: usize },

    #[error(
        "sample interval {interval_m} m cannot produce a strictly increasing distance at grid index {index} ({previous_distance_m} m then {distance_m} m)"
    )]
    UnrepresentableGrid {
        interval_m: f64,
        index: usize,
        previous_distance_m: f64,
        distance_m: f64,
    },
}

impl TrajectoryResult {
    /// Return a checked full-state observation at an in-range downrange distance.
    ///
    /// Unlike [`Self::position_at_range`], this method never clamps an out-of-range request to
    /// the final point.
    pub fn observation_at_range_checked(
        &self,
        distance_m: f64,
    ) -> Result<TrajectoryObservation, TrajectoryObservationError> {
        validate_trajectory(self)?;
        observation_at_range_validated(self, distance_m)
    }

    /// Sample checked observations on a regular distance grid and include the exact terminal
    /// point once.
    ///
    /// Regular grid points are strictly before the actual reached range.  If the terminal is
    /// off-grid it is appended; if it is on-grid it appears only as the final terminal sample.
    /// `max_samples` lets protocol layers enforce a limit smaller than the engine-wide hard cap.
    pub fn sample_observations(
        &self,
        interval_m: f64,
        max_samples: usize,
    ) -> Result<Vec<TrajectoryObservation>, TrajectoryObservationError> {
        validate_trajectory(self)?;

        let effective_limit = max_samples.min(MAX_TRAJECTORY_SAMPLES);
        let count = projected_observation_count(self, interval_m, effective_limit)?;
        let mut observations = Vec::new();
        observations
            .try_reserve_exact(count)
            .map_err(|_| TrajectoryObservationError::AllocationFailed { requested: count })?;

        let first_distance = self.points[0].position.x;
        let terminal_distance = self.points[self.points.len() - 1].position.x;
        let regular_count = count.saturating_sub(1);

        let mut previous_distance_m = None;
        for index in 0..regular_count {
            let distance_m = first_distance + index as f64 * interval_m;
            if let Some(previous_distance_m) = previous_distance_m {
                if distance_m <= previous_distance_m {
                    return Err(TrajectoryObservationError::UnrepresentableGrid {
                        interval_m,
                        index,
                        previous_distance_m,
                        distance_m,
                    });
                }
            }
            // The count projection removes a rounded terminal grid point. Refuse any remaining
            // collision instead of silently thinning the caller's requested grid.
            if distance_m >= terminal_distance {
                return Err(TrajectoryObservationError::UnrepresentableGrid {
                    interval_m,
                    index,
                    previous_distance_m: previous_distance_m.unwrap_or(first_distance),
                    distance_m,
                });
            }
            observations.push(observation_at_range_validated(self, distance_m)?);
            previous_distance_m = Some(distance_m);
        }

        // The exact stored endpoint is authoritative.  A rounded regular-grid calculation is
        // never allowed to replace it or create a repeated terminal observation.
        if observations
            .last()
            .is_some_and(|observation| observation.distance_m == terminal_distance)
        {
            if let Some(last) = observations.last_mut() {
                *last = observation_at_range_validated(self, terminal_distance)?;
            }
        } else {
            observations.push(observation_at_range_validated(self, terminal_distance)?);
        }

        Ok(observations)
    }
}

fn validate_trajectory(result: &TrajectoryResult) -> Result<(), TrajectoryObservationError> {
    if result.points.is_empty() {
        return Err(TrajectoryObservationError::EmptyTrajectory);
    }

    validate_metadata("projectile_mass_kg", result.projectile_mass_kg, |value| {
        value > 0.0
    })?;
    validate_metadata(
        "line_of_sight_height_m",
        result.line_of_sight_height_m,
        |_| true,
    )?;
    validate_metadata(
        "station_speed_of_sound_mps",
        result.station_speed_of_sound_mps,
        |value| value > 0.0,
    )?;

    for (index, point) in result.points.iter().enumerate() {
        validate_point(index, point)?;
        if index > 0 {
            let previous_distance_m = result.points[index - 1].position.x;
            if point.position.x <= previous_distance_m {
                return Err(TrajectoryObservationError::NonMonotonicTrajectory {
                    index,
                    previous_distance_m,
                    distance_m: point.position.x,
                });
            }
        }
    }

    Ok(())
}

fn validate_metadata(
    field: &'static str,
    value: f64,
    predicate: impl FnOnce(f64) -> bool,
) -> Result<(), TrajectoryObservationError> {
    if value.is_finite() && predicate(value) {
        Ok(())
    } else {
        Err(TrajectoryObservationError::InvalidMetadata { field, value })
    }
}

fn validate_point(index: usize, point: &TrajectoryPoint) -> Result<(), TrajectoryObservationError> {
    for (field, value) in [
        ("time", point.time),
        ("position.x", point.position.x),
        ("position.y", point.position.y),
        ("position.z", point.position.z),
        ("velocity_magnitude", point.velocity_magnitude),
        ("kinetic_energy", point.kinetic_energy),
    ] {
        if !value.is_finite() {
            return Err(TrajectoryObservationError::NonFiniteState { index, field });
        }
    }

    for (field, value) in [
        ("time", point.time),
        ("velocity_magnitude", point.velocity_magnitude),
        ("kinetic_energy", point.kinetic_energy),
    ] {
        if value < 0.0 {
            return Err(TrajectoryObservationError::InvalidState {
                index,
                field,
                value,
            });
        }
    }

    Ok(())
}

fn observation_at_range_validated(
    result: &TrajectoryResult,
    distance_m: f64,
) -> Result<TrajectoryObservation, TrajectoryObservationError> {
    if !distance_m.is_finite() {
        return Err(TrajectoryObservationError::NonFiniteQuery { distance_m });
    }

    let minimum_m = result.points[0].position.x;
    let maximum_m = result.points[result.points.len() - 1].position.x;
    if distance_m < minimum_m || distance_m > maximum_m {
        return Err(TrajectoryObservationError::OutOfRange {
            requested_m: distance_m,
            minimum_m,
            maximum_m,
        });
    }

    let upper_index = result
        .points
        .partition_point(|point| point.position.x < distance_m);

    let (time_s, vertical_m, windage_m, speed_mps) = if upper_index < result.points.len()
        && result.points[upper_index].position.x == distance_m
    {
        let point = &result.points[upper_index];
        (
            point.time,
            point.position.y,
            point.position.z,
            point.velocity_magnitude,
        )
    } else {
        // In-range non-exact observations necessarily have a point on each side.
        let lower = &result.points[upper_index - 1];
        let upper = &result.points[upper_index];
        let bracket_span_m = upper.position.x - lower.position.x;
        require_finite_observation("interpolation_span_m", bracket_span_m)?;
        let bracket_offset_m = distance_m - lower.position.x;
        require_finite_observation("interpolation_offset_m", bracket_offset_m)?;
        let alpha = bracket_offset_m / bracket_span_m;
        require_finite_observation("interpolation_fraction", alpha)?;
        (
            checked_lerp(lower.time, upper.time, alpha, "time_s")?,
            checked_lerp(lower.position.y, upper.position.y, alpha, "drop_m")?,
            checked_lerp(lower.position.z, upper.position.z, alpha, "windage_m")?,
            checked_lerp(
                lower.velocity_magnitude,
                upper.velocity_magnitude,
                alpha,
                "speed_mps",
            )?,
        )
    };

    let energy_j = 0.5 * result.projectile_mass_kg * speed_mps * speed_mps;
    require_finite_observation("energy_j", energy_j)?;
    let drop_m = result.line_of_sight_height_m - vertical_m;
    require_finite_observation("drop_m", drop_m)?;
    let mach = speed_mps / result.station_speed_of_sound_mps;
    require_finite_observation("mach", mach)?;

    for (field, value) in [
        ("distance_m", distance_m),
        ("time_s", time_s),
        ("speed_mps", speed_mps),
        ("windage_m", windage_m),
    ] {
        require_finite_observation(field, value)?;
    }

    let terminal = distance_m == maximum_m;
    let mut flags = Vec::with_capacity(4);
    if (0.8..=1.2).contains(&mach) {
        flags.push(TrajectoryObservationFlag::Transonic);
    }
    if mach < 1.0 {
        flags.push(TrajectoryObservationFlag::Subsonic);
    }
    if terminal {
        flags.push(TrajectoryObservationFlag::Terminal);
        if result.termination == TrajectoryTermination::GroundThreshold {
            flags.push(TrajectoryObservationFlag::GroundThreshold);
        }
    }

    Ok(TrajectoryObservation {
        distance_m,
        time_s,
        speed_mps,
        energy_j,
        drop_m,
        windage_m,
        mach,
        flags,
    })
}

fn checked_lerp(
    lower: f64,
    upper: f64,
    alpha: f64,
    field: &'static str,
) -> Result<f64, TrajectoryObservationError> {
    let value = lower + alpha * (upper - lower);
    require_finite_observation(field, value)?;
    Ok(value)
}

fn require_finite_observation(
    field: &'static str,
    value: f64,
) -> Result<(), TrajectoryObservationError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(TrajectoryObservationError::NonFiniteObservation { field })
    }
}

fn projected_observation_count(
    result: &TrajectoryResult,
    interval_m: f64,
    limit: usize,
) -> Result<usize, TrajectoryObservationError> {
    if !interval_m.is_finite() || interval_m <= 0.0 {
        return Err(TrajectoryObservationError::InvalidInterval { interval_m });
    }

    let first_distance = result.points[0].position.x;
    let terminal_distance = result.points[result.points.len() - 1].position.x;
    let span = terminal_distance - first_distance;
    if !span.is_finite() {
        return Err(TrajectoryObservationError::NonFiniteObservation {
            field: "trajectory_span_m",
        });
    }

    let regular_count_f64 = if span > 0.0 {
        (span / interval_m).ceil().max(1.0)
    } else {
        0.0
    };
    if !regular_count_f64.is_finite() || regular_count_f64 > usize::MAX as f64 {
        return Err(TrajectoryObservationError::SampleLimitExceeded {
            requested: usize::MAX,
            limit,
        });
    }

    // Bound work before correcting a possible one-ULP quotient overshoot. In particular, never
    // walk a huge projected count down one element at a time: a finite sub-ULP interval at a
    // large nonzero starting distance can otherwise turn this check into an unbounded loop.
    if regular_count_f64 > limit as f64 {
        let requested = if regular_count_f64 >= usize::MAX as f64 {
            usize::MAX
        } else {
            (regular_count_f64 as usize).saturating_add(1)
        };
        return Err(TrajectoryObservationError::SampleLimitExceeded { requested, limit });
    }

    let mut regular_count = regular_count_f64 as usize;
    // `ceil(span / interval)` can round upward when a mathematically on-grid terminal is
    // represented just above an integer quotient.  Count only generated points that are
    // strictly before the authoritative endpoint.
    if regular_count > 0 {
        let last_regular = first_distance + (regular_count - 1) as f64 * interval_m;
        if last_regular >= terminal_distance {
            regular_count -= 1;
        }
    }
    // Division can also round a quotient just above an integer down to that integer. Probe the
    // next grid point once so that a representable point strictly before the terminal is not
    // omitted. With the count already bounded by `limit` (and therefore exactly representable as
    // f64), one correction in either direction covers the quotient's possible rounding error.
    let next_regular = first_distance + regular_count as f64 * interval_m;
    if next_regular < terminal_distance {
        regular_count = regular_count.saturating_add(1);
    }

    let requested = regular_count.saturating_add(1);
    if requested > limit {
        Err(TrajectoryObservationError::SampleLimitExceeded { requested, limit })
    } else {
        Ok(requested)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    fn point(time: f64, x: f64, y: f64, z: f64, speed: f64) -> TrajectoryPoint {
        TrajectoryPoint {
            time,
            position: Vector3::new(x, y, z),
            velocity_magnitude: speed,
            kinetic_energy: 0.5 * 0.02 * speed * speed,
            drag_coefficient: None,
        }
    }

    fn result(
        points: Vec<TrajectoryPoint>,
        termination: TrajectoryTermination,
    ) -> TrajectoryResult {
        let terminal = points.last().expect("test trajectory must not be empty");
        TrajectoryResult {
            max_range: terminal.position.x,
            max_height: points
                .iter()
                .map(|point| point.position.y)
                .fold(f64::NEG_INFINITY, f64::max),
            time_of_flight: terminal.time,
            impact_velocity: terminal.velocity_magnitude,
            impact_energy: terminal.kinetic_energy,
            projectile_mass_kg: 0.02,
            line_of_sight_height_m: 1.0,
            station_speed_of_sound_mps: 340.0,
            termination,
            points,
            sampled_points: None,
            min_pitch_damping: None,
            transonic_mach: None,
            angular_state: None,
            max_yaw_angle: None,
            max_precession_angle: None,
            aerodynamic_jump: None,
            mach_1_2_distance_m: None,
            mach_1_0_distance_m: None,
            mach_0_9_distance_m: None,
        }
    }

    fn synthetic_result() -> TrajectoryResult {
        result(
            vec![
                point(0.0, 0.0, 0.5, -0.4, 680.0),
                point(2.0, 100.0, 1.5, 0.4, 340.0),
            ],
            TrajectoryTermination::MaxRange,
        )
    }

    #[test]
    fn interpolates_full_state_with_documented_drop_and_windage_signs() {
        let trajectory = synthetic_result();

        let first = trajectory
            .observation_at_range_checked(25.0)
            .expect("in-range observation");
        assert_relative_eq!(first.time_s, 0.5);
        assert_relative_eq!(first.speed_mps, 595.0);
        assert_relative_eq!(first.energy_j, 3540.25);
        assert_relative_eq!(first.drop_m, 0.25);
        assert_relative_eq!(first.windage_m, -0.2);
        assert_relative_eq!(first.mach, 1.75);

        let second = trajectory
            .observation_at_range_checked(75.0)
            .expect("in-range observation");
        assert_relative_eq!(second.drop_m, -0.25);
        assert_relative_eq!(second.windage_m, 0.2);
    }

    #[test]
    fn preserves_exact_endpoints_and_marks_only_the_terminal_endpoint() {
        let trajectory = synthetic_result();

        let muzzle = trajectory
            .observation_at_range_checked(0.0)
            .expect("muzzle endpoint");
        assert_eq!(muzzle.time_s, 0.0);
        assert_eq!(muzzle.speed_mps, 680.0);
        assert!(!muzzle.flags.contains(&TrajectoryObservationFlag::Terminal));

        let terminal = trajectory
            .observation_at_range_checked(100.0)
            .expect("terminal endpoint");
        assert_eq!(terminal.time_s, 2.0);
        assert_eq!(terminal.speed_mps, 340.0);
        assert_eq!(terminal.drop_m, -0.5);
        assert_eq!(terminal.windage_m, 0.4);
        assert!(terminal
            .flags
            .contains(&TrajectoryObservationFlag::Transonic));
        assert!(terminal
            .flags
            .contains(&TrajectoryObservationFlag::Terminal));
    }

    #[test]
    fn rejects_out_of_range_and_non_finite_queries_instead_of_clamping() {
        let trajectory = synthetic_result();

        for distance_m in [-0.001, 100.001] {
            assert!(matches!(
                trajectory.observation_at_range_checked(distance_m),
                Err(TrajectoryObservationError::OutOfRange { .. })
            ));
        }
        for distance_m in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                trajectory.observation_at_range_checked(distance_m),
                Err(TrajectoryObservationError::NonFiniteQuery { .. })
            ));
        }

        // The established convenience API retains its historical clamping behavior.
        assert_eq!(
            trajectory.position_at_range(101.0),
            Some(Vector3::new(100.0, 1.5, 0.4))
        );
    }

    #[test]
    fn regular_grid_appends_an_off_grid_terminal_and_deduplicates_an_on_grid_terminal() {
        let off_grid = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, 95.0, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::GroundThreshold,
        );
        let samples = off_grid
            .sample_observations(30.0, 10)
            .expect("off-grid samples");
        assert_eq!(
            samples
                .iter()
                .map(|sample| sample.distance_m)
                .collect::<Vec<_>>(),
            vec![0.0, 30.0, 60.0, 90.0, 95.0]
        );
        let terminal = samples.last().expect("terminal sample");
        assert!(terminal
            .flags
            .contains(&TrajectoryObservationFlag::Terminal));
        assert!(terminal
            .flags
            .contains(&TrajectoryObservationFlag::GroundThreshold));

        let on_grid = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, 90.0, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        let samples = on_grid
            .sample_observations(30.0, 10)
            .expect("on-grid samples");
        assert_eq!(
            samples
                .iter()
                .map(|sample| sample.distance_m)
                .collect::<Vec<_>>(),
            vec![0.0, 30.0, 60.0, 90.0]
        );

        let fractional = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, 0.15, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert_eq!(
            fractional
                .sample_observations(0.1, 10)
                .expect("fractional samples")
                .iter()
                .map(|sample| sample.distance_m)
                .collect::<Vec<_>>(),
            vec![0.0, 0.1, 0.15]
        );

        let rounded_grid_point = 3.0_f64 * 0.3;
        let just_past_grid = f64::from_bits(rounded_grid_point.to_bits() + 1);
        let rounded_quotient = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, just_past_grid, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert_eq!(
            rounded_quotient
                .sample_observations(0.3, 5)
                .expect("a rounded quotient retains the last regular grid point")
                .iter()
                .map(|sample| sample.distance_m.to_bits())
                .collect::<Vec<_>>(),
            [0.0, 0.3, 0.6, rounded_grid_point, just_past_grid]
                .map(f64::to_bits)
                .to_vec()
        );
    }

    #[test]
    fn rejects_non_finite_state_metadata_and_derived_values() {
        let mut non_finite_state = synthetic_result();
        non_finite_state.points[1].position.y = f64::NAN;
        assert!(matches!(
            non_finite_state.observation_at_range_checked(50.0),
            Err(TrajectoryObservationError::NonFiniteState {
                index: 1,
                field: "position.y"
            })
        ));

        let mut invalid_metadata = synthetic_result();
        invalid_metadata.station_speed_of_sound_mps = 0.0;
        assert!(matches!(
            invalid_metadata.observation_at_range_checked(50.0),
            Err(TrajectoryObservationError::InvalidMetadata {
                field: "station_speed_of_sound_mps",
                ..
            })
        ));

        let mut overflowing_energy = synthetic_result();
        overflowing_energy.points[0].velocity_magnitude = f64::MAX;
        overflowing_energy.points[0].kinetic_energy = f64::MAX;
        assert!(matches!(
            overflowing_energy.observation_at_range_checked(0.0),
            Err(TrajectoryObservationError::NonFiniteObservation { field: "energy_j" })
        ));

        let overflowing_bracket = result(
            vec![
                point(0.0, -f64::MAX, 1.0, 0.0, 500.0),
                point(1.0, f64::MAX, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert!(matches!(
            overflowing_bracket.observation_at_range_checked(0.0),
            Err(TrajectoryObservationError::NonFiniteObservation {
                field: "interpolation_span_m"
            })
        ));
    }

    #[test]
    fn enforces_caller_and_engine_sample_caps_before_allocation() {
        let small = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, 4.0, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert_eq!(
            small
                .sample_observations(1.0, 5)
                .expect("exact caller limit")
                .len(),
            5
        );
        assert!(matches!(
            small.sample_observations(1.0, 4),
            Err(TrajectoryObservationError::SampleLimitExceeded {
                requested: 5,
                limit: 4
            })
        ));

        let at_engine_limit = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, (MAX_TRAJECTORY_SAMPLES - 1) as f64, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert_eq!(
            projected_observation_count(&at_engine_limit, 1.0, MAX_TRAJECTORY_SAMPLES),
            Ok(MAX_TRAJECTORY_SAMPLES)
        );

        let above_engine_limit = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, MAX_TRAJECTORY_SAMPLES as f64, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );
        assert!(matches!(
            projected_observation_count(
                &above_engine_limit,
                1.0,
                MAX_TRAJECTORY_SAMPLES
            ),
            Err(TrajectoryObservationError::SampleLimitExceeded {
                requested,
                limit: MAX_TRAJECTORY_SAMPLES
            }) if requested == MAX_TRAJECTORY_SAMPLES + 1
        ));
    }

    #[test]
    fn rejects_huge_or_unrepresentable_grids_in_bounded_work() {
        let first_distance = 1.0e300_f64;
        let terminal_distance = f64::from_bits(first_distance.to_bits() + 1);
        let span = terminal_distance - first_distance;
        let trajectory = result(
            vec![
                point(0.0, first_distance, 1.0, 0.0, 500.0),
                point(1.0, terminal_distance, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );

        assert!(matches!(
            trajectory.sample_observations(span / 1.0e15, MAX_TRAJECTORY_SAMPLES),
            Err(TrajectoryObservationError::SampleLimitExceeded { .. })
        ));
        assert!(matches!(
            trajectory.sample_observations(span / 10.0, 20),
            Err(TrajectoryObservationError::UnrepresentableGrid { index: 1, .. })
        ));
    }

    #[test]
    fn interval_ratio_underflow_still_includes_both_endpoints() {
        let terminal_distance = f64::from_bits(1);
        let trajectory = result(
            vec![
                point(0.0, 0.0, 1.0, 0.0, 500.0),
                point(1.0, terminal_distance, 0.5, 0.0, 400.0),
            ],
            TrajectoryTermination::MaxRange,
        );

        let observations = trajectory
            .sample_observations(f64::MAX, 2)
            .expect("a positive span retains its first and terminal observations");
        assert_eq!(observations.len(), 2);
        assert_eq!(observations[0].distance_m.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            observations[1].distance_m.to_bits(),
            terminal_distance.to_bits()
        );
    }

    #[test]
    fn rejects_duplicate_or_reversing_distances() {
        for terminal_x in [0.0, -1.0] {
            let trajectory = result(
                vec![
                    point(0.0, 0.0, 1.0, 0.0, 500.0),
                    point(1.0, terminal_x, 0.5, 0.0, 400.0),
                ],
                TrajectoryTermination::MaxRange,
            );
            assert!(matches!(
                trajectory.observation_at_range_checked(0.0),
                Err(TrajectoryObservationError::NonMonotonicTrajectory { .. })
            ));
        }
    }
}
