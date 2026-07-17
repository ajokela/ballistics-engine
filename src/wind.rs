use nalgebra::Vector3;
use std::cmp::Ordering;
use std::f64::consts::PI;
use thiserror::Error;

/// Conversion constant from KMH to MPS
const KMH_TO_MPS: f64 = 1000.0 / 3600.0;

/// THE wind-vector builder (McCoy frame: x downrange, y up, z right).
/// Horizontal wind uses the wind-FROM convention (0 = headwind, PI/2 = from
/// the right); `vertical_mps` is positive-updraft and lands on y UNSCALED —
/// boundary-layer shear models scale horizontal flow only (MBA-728 decision).
pub fn wind_vector(speed_mps: f64, direction_rad: f64, vertical_mps: f64) -> Vector3<f64> {
    Vector3::new(
        -speed_mps * direction_rad.cos(),
        vertical_mps,
        -speed_mps * direction_rad.sin(),
    )
}

/// One downrange wind segment. `vertical_mps` (m/s, positive = updraft) feeds
/// straight into the segment's wind vector via [`wind_vector`] (MBA-728);
/// boundary-layer shear scales horizontal wind only, so vertical passes
/// through unscaled wherever shear is applied on top of a segment.
///
/// This matches the Python WindSock interface.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct WindSegment {
    pub speed_kmh: f64,
    pub angle_deg: f64,
    pub until_m: f64,
    pub vertical_mps: f64,
}

impl WindSegment {
    /// The historical 3-field constructor (vertical 0.0) — used by every
    /// pre-existing call site.
    pub fn new(speed_kmh: f64, angle_deg: f64, until_m: f64) -> Self {
        Self {
            speed_kmh,
            angle_deg,
            until_m,
            vertical_mps: 0.0,
        }
    }
}

/// Sort wind segments by their `until_distance_m` threshold.
///
/// Shared by [`WindSock`] and the low-level trajectory integrator so every segmented-wind path
/// applies the same interval ordering.
pub(crate) fn sort_wind_segments_by_distance(segments: &mut [WindSegment]) {
    segments.sort_by(|a, b| match (a.until_m.is_nan(), b.until_m.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => {
            a.until_m
                .partial_cmp(&b.until_m)
                .expect("non-NaN distances are ordered")
        }
    });
}

/// Which [`WindSegment`] field failed validation (MBA-1338).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindSegmentField {
    SpeedKmh,
    AngleDeg,
    UntilM,
    VerticalMps,
}

impl WindSegmentField {
    /// The struct field name, as it appears in error messages and the public API.
    pub fn name(self) -> &'static str {
        match self {
            WindSegmentField::SpeedKmh => "speed_kmh",
            WindSegmentField::AngleDeg => "angle_deg",
            WindSegmentField::UntilM => "until_m",
            WindSegmentField::VerticalMps => "vertical_mps",
        }
    }
}

impl std::fmt::Display for WindSegmentField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// The validation rule a [`WindSegment`] field violated (MBA-1338).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindSegmentRule {
    /// The value must be finite (`angle_deg`, `vertical_mps`).
    Finite,
    /// The value must be finite and `>= 0` (`speed_kmh`).
    FiniteAndNonNegative,
    /// The value must be finite and `> 0` (`until_m`).
    FiniteAndPositive,
}

impl std::fmt::Display for WindSegmentRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Phrasing matches the historical String errors exactly so consumers that surface
        // Display output (FastSolution, solver boundaries, bindings) see unchanged messages.
        f.write_str(match self {
            WindSegmentRule::Finite => "finite",
            WindSegmentRule::FiniteAndNonNegative => "finite and non-negative",
            WindSegmentRule::FiniteAndPositive => "finite and greater than zero",
        })
    }
}

/// A malformed [`WindSegment`] rejected at a checked construction/solve boundary (MBA-1338).
///
/// `index` is the segment's position in the vector **as supplied by the caller** —
/// validation runs before any sorting or normalization, so the index always refers to
/// the caller's own ordering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[error("wind.segments[{index}].{field} must be {rule}")]
pub struct WindSegmentError {
    pub index: usize,
    pub field: WindSegmentField,
    pub rule: WindSegmentRule,
}

/// Validate every segment, reporting the first violation as a typed error carrying the
/// caller's segment index, the offending field, and the violated rule (MBA-1338).
///
/// Checked APIs ([`WindSock::try_new`], `try_integrate_trajectory`,
/// `try_solve_trajectory_rust`) call this before sorting, vector precomputation, or
/// producing any trajectory points.
pub fn validate_wind_segments(segments: &[WindSegment]) -> Result<(), WindSegmentError> {
    for (index, segment) in segments.iter().enumerate() {
        let err = |field, rule| Err(WindSegmentError { index, field, rule });

        if !segment.speed_kmh.is_finite() || segment.speed_kmh < 0.0 {
            return err(
                WindSegmentField::SpeedKmh,
                WindSegmentRule::FiniteAndNonNegative,
            );
        }
        if !segment.angle_deg.is_finite() {
            return err(WindSegmentField::AngleDeg, WindSegmentRule::Finite);
        }
        if !segment.until_m.is_finite() || segment.until_m <= 0.0 {
            return err(WindSegmentField::UntilM, WindSegmentRule::FiniteAndPositive);
        }
        if !segment.vertical_mps.is_finite() {
            return err(WindSegmentField::VerticalMps, WindSegmentRule::Finite);
        }
    }

    Ok(())
}

/// Wind condition handler for trajectory calculations
#[derive(Debug, Clone)]
pub struct WindSock {
    /// Sorted wind segments by distance
    winds: Vec<WindSegment>,
    /// Precomputed wind vector for each segment (parallel to `winds`). The Monte-Carlo RK4
    /// kernel queries wind 4x per step, so caching avoids recomputing sin/cos every call.
    wind_vecs: Vec<Vector3<f64>>,
    /// Current segment index
    current: usize,
    /// Distance where next segment starts
    next_range: f64,
    /// Current wind vector
    current_vec: Vector3<f64>,
    /// Validation result captured before sorting, preserving the caller's segment index.
    validation_error: Option<WindSegmentError>,
}

impl WindSock {
    /// Create a new WindSock from wind segments.
    ///
    /// **Legacy/unchecked entry point** (MBA-1338): this constructor is infallible for API
    /// compatibility — a malformed segment is captured as a deferred error that the solver
    /// consumes at its pre-integration validation boundary, but a direct library consumer
    /// gets a constructed sock with no structured signal. Prefer [`WindSock::try_new`]
    /// (or `TryFrom<Vec<WindSegment>>`), which rejects malformed segments up front.
    ///
    /// Args:
    ///     segments: List of (speed_kmh, angle_deg, until_distance_m) tuples
    pub fn new(segments: Vec<WindSegment>) -> Self {
        // Validate before sorting so an error can identify the segment index supplied by
        // the caller; defer it rather than fail (see doc above).
        let validation_error = validate_wind_segments(&segments).err();
        Self::build_unchecked(segments, validation_error)
    }

    /// Create a new WindSock, rejecting malformed segments with a typed error (MBA-1338).
    ///
    /// Validation runs **before** sorting and wind-vector precomputation, and the error's
    /// `index` refers to the caller's own segment ordering. Rejects non-finite or negative
    /// `speed_kmh`, non-finite `angle_deg`, non-finite or non-positive `until_m`, and
    /// non-finite `vertical_mps`.
    pub fn try_new(segments: Vec<WindSegment>) -> Result<Self, WindSegmentError> {
        validate_wind_segments(&segments)?;
        Ok(Self::build_unchecked(segments, None))
    }

    /// Shared constructor body: sort, precompute vectors, seed the cursor.
    fn build_unchecked(
        mut segments: Vec<WindSegment>,
        validation_error: Option<WindSegmentError>,
    ) -> Self {
        // Sort segments by distance, handling NaN safely by treating it as greater than any value
        sort_wind_segments_by_distance(&mut segments);

        // Precompute each segment's wind vector once (depends only on its speed/angle).
        let wind_vecs: Vec<Vector3<f64>> = segments.iter().map(Self::calc_vec).collect();

        let (current, next_range, current_vec) = if segments.is_empty() {
            (0, f64::INFINITY, Vector3::zeros())
        } else {
            (0, segments[0].until_m, wind_vecs[0])
        };

        WindSock {
            winds: segments,
            wind_vecs,
            current,
            next_range,
            current_vec,
            validation_error,
        }
    }

    /// Return any malformed-segment error captured before normalization.
    pub(crate) fn validate_segments(&self) -> Result<(), String> {
        self.validation_error
            .map_or(Ok(()), |err| Err(err.to_string()))
    }

    /// Calculate wind vector from wind segment
    fn calc_vec(seg: &WindSegment) -> Vector3<f64> {
        // Convert kmh to m/s
        let speed_mps = seg.speed_kmh * KMH_TO_MPS;
        // Preserve the historical multiply-then-divide result for ordinary angles, while
        // avoiding intermediate overflow for otherwise-valid finite values near f64::MAX.
        let angle_rad = if seg.angle_deg.abs() <= f64::MAX / PI {
            seg.angle_deg * PI / 180.0
        } else {
            seg.angle_deg / 180.0 * PI
        };

        // Wind convention (matching trajectory coordinates):
        // 0° = headwind (from front, affects -x downrange)
        // 90° = wind from right (affects -z lateral)
        // 180° = tailwind (from back, affects +x downrange)
        // 270° = wind from left (affects +z lateral)
        //
        // McCoy convention: x=downrange, y=vertical, z=lateral. Vertical (MBA-728) passes
        // straight through per-segment; it is not derived from speed_kmh/angle_deg.
        wind_vector(speed_mps, angle_rad, seg.vertical_mps)
    }

    /// Upper bound (m/s) on the wind speed any segment can contribute, horizontal
    /// plus vertical. Feeds the solver's integration divergence guard (MBA-1293).
    pub fn max_speed_mps(&self) -> f64 {
        self.winds
            .iter()
            .map(|seg| seg.speed_kmh.abs() * KMH_TO_MPS + seg.vertical_mps.abs())
            .fold(0.0, f64::max)
    }

    /// Crosswind at the muzzle in the aerodynamic-jump convention: positive means wind from
    /// the right. `None` means there are no segmented winds, while a real zero-crosswind
    /// segment returns `Some(0.0)` so callers do not incorrectly fall back to scalar wind.
    pub(crate) fn muzzle_crosswind_from_right_mps(&self) -> Option<f64> {
        (!self.winds.is_empty()).then(|| -self.vector_for_range_stateless(0.0)[2])
    }

    /// Get wind vector for a given range
    ///
    /// Note: This modifies internal state and expects monotonically increasing ranges
    /// For trajectory integration, we need a stateless version
    pub fn vector_for_range(&mut self, range_m: f64) -> Vector3<f64> {
        // Handle NaN
        if range_m.is_nan() {
            return Vector3::zeros();
        }

        // Advance the cursor across however many segments the query skipped (a single `if`
        // returned a stale vector when a monotonic query jumped past a whole short segment).
        while range_m >= self.next_range && self.current < self.winds.len() {
            self.current += 1;
            if self.current >= self.winds.len() {
                self.current_vec = Vector3::zeros();
                self.next_range = f64::INFINITY;
            } else {
                self.current_vec = self.wind_vecs[self.current];
                self.next_range = self.winds[self.current].until_m;
            }
        }

        self.current_vec
    }

    /// Get wind vector for a given range (stateless version)
    ///
    /// This version doesn't modify internal state and is safe for numerical integration
    /// where the same range might be queried multiple times or out of order
    pub fn vector_for_range_stateless(&self, range_m: f64) -> Vector3<f64> {
        // Handle NaN
        if range_m.is_nan() {
            return Vector3::zeros();
        }

        // Find the appropriate segment (precomputed vector — no per-call trig).
        for (i, segment) in self.winds.iter().enumerate() {
            if range_m < segment.until_m {
                return self.wind_vecs[i];
            }
        }

        // Beyond all segments
        Vector3::zeros()
    }
}

/// Checked conversion mirroring [`WindSock::try_new`] (MBA-1338).
impl TryFrom<Vec<WindSegment>> for WindSock {
    type Error = WindSegmentError;

    fn try_from(segments: Vec<WindSegment>) -> Result<Self, Self::Error> {
        WindSock::try_new(segments)
    }
}

/// Parse a `"SPEED:ANGLE:UNTIL_DISTANCE[:VERTICAL]"` string into a [`WindSegment`]
/// `(speed_kmh, angle_deg, until_distance_m, vertical_mps)`.
///
/// `imperial`: when true, SPEED is mph and UNTIL_DISTANCE is yards; otherwise
/// SPEED is m/s and UNTIL_DISTANCE is meters. ANGLE is always degrees in the
/// wind-FROM convention (0 = headwind, 90 = from the right). The optional 4th
/// field, VERTICAL, is ALWAYS m/s (positive = updraft, raises POI) regardless of
/// `imperial` — it does not follow --units, matching how [`WindSegment::vertical_mps`]
/// stores it. This speed-in-display-units-but-vertical-always-m/s asymmetry is
/// unit-honest (it mirrors the struct field name) even though it reads oddly next
/// to SPEED. Omitting the 4th field keeps the historical 3-field behavior
/// (vertical wind 0.0). Shared by the CLI (`--wind-segment`) and the WASM
/// front-ends so they parse identically.
pub fn parse_wind_segment_str(s: &str, imperial: bool) -> Result<WindSegment, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 && parts.len() != 4 {
        return Err(format!(
            "invalid wind segment '{s}': expected SPEED:ANGLE:UNTIL_DISTANCE[:VERTICAL] \
             (three or four colon-separated numbers; the optional 4th field VERTICAL is always \
             m/s, positive = updraft, regardless of --units)"
        ));
    }
    let num = |i: usize, name: &str| -> Result<f64, String> {
        parts[i].trim().parse::<f64>().map_err(|_| {
            format!("invalid wind segment '{s}': {name} '{}' is not a number", parts[i])
        })
    };
    let speed = num(0, "speed")?;
    let angle = num(1, "angle")?;
    let until = num(2, "until-distance")?;
    let vertical = if parts.len() == 4 {
        num(3, "vertical")?
    } else {
        0.0
    };
    if !speed.is_finite() || !angle.is_finite() || !until.is_finite() || !vertical.is_finite() {
        return Err(format!(
            "invalid wind segment '{s}': speed, angle, until-distance, and vertical (m/s, \
             positive = updraft) must be finite numbers"
        ));
    }
    if speed < 0.0 {
        return Err(format!("invalid wind segment '{s}': speed must be >= 0"));
    }
    if until <= 0.0 {
        return Err(format!("invalid wind segment '{s}': until-distance must be > 0"));
    }
    let (speed_kmh, until_m) = if imperial {
        (speed * 1.609344, until * 0.9144) // mph -> km/h, yards -> meters
    } else {
        (speed * 3.6, until) // m/s -> km/h, meters -> meters
    };
    let mut segment = WindSegment::new(speed_kmh, angle, until_m);
    segment.vertical_mps = vertical;
    Ok(segment)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn segment_sort_is_stable_and_places_nan_endpoints_last() {
        let mut segments = vec![
            WindSegment::new(10.0, 0.0, f64::NAN),
            WindSegment::new(20.0, 0.0, 100.0),
            WindSegment::new(30.0, 0.0, 100.0),
            WindSegment::new(40.0, 0.0, f64::INFINITY),
            WindSegment::new(50.0, 0.0, f64::NEG_INFINITY),
            WindSegment::new(60.0, 0.0, f64::NAN),
        ];

        sort_wind_segments_by_distance(&mut segments);

        assert_eq!(segments[0].speed_kmh, 50.0); // -inf first
        assert_eq!(segments[1].speed_kmh, 20.0); // equal endpoints retain input order
        assert_eq!(segments[2].speed_kmh, 30.0);
        assert_eq!(segments[3].speed_kmh, 40.0); // +inf after finite endpoints
        assert_eq!(segments[4].speed_kmh, 10.0); // NaNs last and stable
        assert_eq!(segments[5].speed_kmh, 60.0);
        assert!(segments[4].until_m.is_nan() && segments[5].until_m.is_nan());
    }

    #[test]
    fn test_wind_sock_empty() {
        let sock = WindSock::new(vec![]);
        assert_eq!(sock.vector_for_range_stateless(50.0), Vector3::zeros());
        assert_eq!(sock.muzzle_crosswind_from_right_mps(), None);
    }

    #[test]
    fn try_new_accepts_valid_and_empty_segments() {
        assert!(WindSock::try_new(vec![]).is_ok());
        let sock = WindSock::try_new(vec![WindSegment::new(16.0934, 90.0, 100.0)]).unwrap();
        assert!(sock.vector_for_range_stateless(50.0).norm() > 0.0);
        // TryFrom mirrors try_new.
        assert!(WindSock::try_from(vec![WindSegment::new(10.0, 0.0, 50.0)]).is_ok());
        assert!(WindSock::try_from(vec![WindSegment::new(f64::NAN, 0.0, 50.0)]).is_err());
    }

    #[test]
    fn try_new_rejects_each_field_with_typed_index_field_rule() {
        // MBA-1338 acceptance criteria: non-finite or negative speed_kmh, non-finite
        // angle_deg, non-finite or non-positive until_m, non-finite vertical_mps.
        let ok = WindSegment::new(10.0, 0.0, 100.0);
        let cases: Vec<(WindSegment, WindSegmentField, WindSegmentRule)> = vec![
            (
                WindSegment::new(f64::NAN, 0.0, 100.0),
                WindSegmentField::SpeedKmh,
                WindSegmentRule::FiniteAndNonNegative,
            ),
            (
                WindSegment::new(-1.0, 0.0, 100.0),
                WindSegmentField::SpeedKmh,
                WindSegmentRule::FiniteAndNonNegative,
            ),
            (
                WindSegment::new(10.0, f64::INFINITY, 100.0),
                WindSegmentField::AngleDeg,
                WindSegmentRule::Finite,
            ),
            (
                WindSegment::new(10.0, 0.0, 0.0),
                WindSegmentField::UntilM,
                WindSegmentRule::FiniteAndPositive,
            ),
            (
                WindSegment::new(10.0, 0.0, f64::NEG_INFINITY),
                WindSegmentField::UntilM,
                WindSegmentRule::FiniteAndPositive,
            ),
            (
                WindSegment {
                    vertical_mps: f64::NAN,
                    ..ok
                },
                WindSegmentField::VerticalMps,
                WindSegmentRule::Finite,
            ),
        ];
        for (bad, field, rule) in cases {
            // Put the bad segment at index 1 behind a valid one so the reported index is
            // meaningful (and would be perturbed if validation ran after sorting).
            let err = WindSock::try_new(vec![ok, bad]).unwrap_err();
            assert_eq!(err.index, 1, "index for {field}");
            assert_eq!(err.field, field);
            assert_eq!(err.rule, rule);
        }
    }

    #[test]
    fn try_new_error_index_is_the_callers_presort_index() {
        // The bad segment sorts FIRST by until_m; the reported index must still be the
        // caller's (2), proving validation runs before sorting.
        let err = WindSock::try_new(vec![
            WindSegment::new(10.0, 0.0, 500.0),
            WindSegment::new(20.0, 0.0, 300.0),
            WindSegment::new(f64::NAN, 0.0, 1.0),
        ])
        .unwrap_err();
        assert_eq!(err.index, 2);
        assert_eq!(err.field, WindSegmentField::SpeedKmh);
    }

    #[test]
    fn wind_segment_error_display_matches_legacy_strings() {
        // The deferred WindSock::new path and FastSolution surface this Display output;
        // it must stay byte-identical to the historical String errors.
        let err = WindSock::try_new(vec![WindSegment::new(f64::NAN, 0.0, 100.0)]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "wind.segments[0].speed_kmh must be finite and non-negative"
        );
        let err = WindSock::try_new(vec![WindSegment::new(10.0, 0.0, -5.0)]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "wind.segments[0].until_m must be finite and greater than zero"
        );
        let err = WindSock::try_new(vec![WindSegment {
            vertical_mps: f64::INFINITY,
            ..WindSegment::new(10.0, 0.0, 100.0)
        }])
        .unwrap_err();
        assert_eq!(
            err.to_string(),
            "wind.segments[0].vertical_mps must be finite"
        );
    }

    #[test]
    fn infallible_new_still_defers_the_same_error() {
        // Legacy behavior preserved: WindSock::new constructs, and the solver-boundary
        // validate_segments() reports the identical message try_new would have raised.
        let sock = WindSock::new(vec![WindSegment::new(10.0, f64::NAN, 100.0)]);
        assert_eq!(
            sock.validate_segments().unwrap_err(),
            "wind.segments[0].angle_deg must be finite"
        );
    }

    #[test]
    fn muzzle_crosswind_distinguishes_an_explicit_zero_segment() {
        let sock = WindSock::new(vec![WindSegment::new(0.0, 90.0, 100.0)]);
        assert_eq!(sock.muzzle_crosswind_from_right_mps(), Some(0.0));
    }

    #[test]
    fn muzzle_crosswind_uses_the_sorted_muzzle_segment_and_wind_from_sign() {
        let from_right = WindSock::new(vec![
            WindSegment::new(32.18688, 270.0, 5000.0),
            WindSegment::new(16.09344, 90.0, 100.0),
        ]);
        let right_mps = from_right.muzzle_crosswind_from_right_mps().unwrap();
        assert!((right_mps - 4.4704).abs() < 1e-12);

        let from_left = WindSock::new(vec![WindSegment::new(16.09344, 270.0, 100.0)]);
        let left_mps = from_left.muzzle_crosswind_from_right_mps().unwrap();
        assert!((left_mps + 4.4704).abs() < 1e-12);
    }

    #[test]
    fn test_wind_sock_single_segment() {
        // 16.0934 kmh (10 mph) @ 90° until 100m
        let sock = WindSock::new(vec![WindSegment::new(16.0934, 90.0, 100.0)]);

        // Should have wind before 100m
        let vec_50 = sock.vector_for_range_stateless(50.0);
        println!("vec_50 = [{}, {}, {}]", vec_50[0], vec_50[1], vec_50[2]);
        assert!(vec_50.norm() > 0.0);
        // 90° wind from right (crosswind, McCoy): negative Z (lateral), zero Y, near-zero X (downrange)
        assert!(
            vec_50[2] < 0.0,
            "Z (lateral) should be negative for 90° wind, got {}",
            vec_50[2]
        );
        assert_eq!(vec_50[1], 0.0); // Zero Y component (vertical_mps defaults to 0.0, MBA-728)
        assert!(
            vec_50[0].abs() < 0.01,
            "X (downrange) should be nearly zero for 90° wind, got {}",
            vec_50[0]
        );

        // No wind after 100m
        let vec_150 = sock.vector_for_range_stateless(150.0);
        assert_eq!(vec_150, Vector3::zeros());
    }

    #[test]
    fn test_wind_sock_multiple_segments() {
        // Multiple wind segments (in kmh)
        let sock = WindSock::new(vec![
            WindSegment::new(16.0934, 90.0, 50.0),  // 10 mph @ 90° until 50m
            WindSegment::new(24.1401, 45.0, 100.0), // 15 mph @ 45° until 100m
            WindSegment::new(8.0467, 180.0, 200.0), // 5 mph @ 180° until 200m
        ]);

        // Test each segment
        let vec_25 = sock.vector_for_range_stateless(25.0);
        println!("vec_25 = [{}, {}, {}]", vec_25[0], vec_25[1], vec_25[2]);
        assert!(vec_25.norm() > 0.0);
        assert!(vec_25[2] < 0.0, "90° wind should have negative Z (lateral)"); // 90° wind from right

        let vec_75 = sock.vector_for_range_stateless(75.0);
        println!("vec_75 = [{}, {}, {}]", vec_75[0], vec_75[1], vec_75[2]);
        assert!(vec_75.norm() > vec_25.norm()); // 15 mph > 10 mph
        assert!(vec_75[0] < 0.0); // 45° wind has negative X component
        assert!(vec_75[2] < 0.0); // 45° wind has negative Z component

        let vec_150 = sock.vector_for_range_stateless(150.0);
        println!("vec_150 = [{}, {}, {}]", vec_150[0], vec_150[1], vec_150[2]);
        assert!(vec_150.norm() < vec_75.norm()); // 5 mph < 15 mph
        assert!(
            vec_150[2].abs() < 0.01,
            "180° wind should have near-zero Z (lateral), got {}",
            vec_150[2]
        ); // 180° wind (from behind)
        assert!(
            vec_150[0] > 0.0,
            "180° wind should have positive X (tailwind, downrange), got {}",
            vec_150[0]
        ); // Tailwind

        let vec_250 = sock.vector_for_range_stateless(250.0);
        assert_eq!(vec_250, Vector3::zeros()); // Beyond all segments
    }

    #[test]
    fn test_wind_conversion() {
        // Test conversion: 16.0934 km/h = 4.47 m/s
        let sock = WindSock::new(vec![WindSegment::new(16.0934, 0.0, 100.0)]);
        let vec = sock.vector_for_range_stateless(50.0);

        let expected_speed = 16.0934 * KMH_TO_MPS;
        assert!((vec.norm() - expected_speed).abs() < 0.01);
    }

    #[test]
    fn test_wind_sock_boundary_is_upper_exclusive() {
        // A segment's `until_distance_m` is exclusive: a query exactly at the
        // boundary rolls to the next segment.
        let sock = WindSock::new(vec![
            WindSegment::new(16.0934, 90.0, 100.0),
            WindSegment::new(32.1868, 270.0, 200.0),
        ]);
        // Just below 100 m -> first segment (90deg, negative Z).
        assert!(sock.vector_for_range_stateless(99.999)[2] < 0.0);
        // Exactly 100 m -> second segment (270deg, positive Z).
        assert!(sock.vector_for_range_stateless(100.0)[2] > 0.0);
        // Beyond the last boundary -> zero.
        assert_eq!(sock.vector_for_range_stateless(200.0), Vector3::zeros());
    }

    #[test]
    fn calc_vec_passes_through_segment_vertical_unscaled() {
        // MBA-728: a segment's vertical_mps must land unchanged on the wind vector's Y
        // component, independent of speed/angle.
        let seg = WindSegment {
            speed_kmh: 16.0934,
            angle_deg: 90.0,
            until_m: 100.0,
            vertical_mps: 3.0,
        };
        let vec = WindSock::calc_vec(&seg);
        assert_eq!(vec[1], 3.0);
    }

    #[test]
    fn calc_vec_zero_vertical_segment_keeps_zero_y() {
        // Upgrades (does not replace) the zero-Y assertions above: the historical
        // 3-field constructor still yields vertical_mps == 0.0 -> Y == 0.0.
        let seg = WindSegment::new(16.0934, 90.0, 100.0);
        assert_eq!(seg.vertical_mps, 0.0);
        let vec = WindSock::calc_vec(&seg);
        assert_eq!(vec[1], 0.0);
    }

    #[test]
    fn test_parse_wind_segment_str_units() {
        // Imperial: 10 mph -> 16.0934 km/h, 100 yd -> 91.44 m.
        let seg = parse_wind_segment_str("10:90:100", true).unwrap();
        assert!((seg.speed_kmh - 16.09344).abs() < 1e-4);
        assert_eq!(seg.angle_deg, 90.0);
        assert!((seg.until_m - 91.44).abs() < 1e-4);

        // Metric: 5 m/s -> 18 km/h, 200 m stays 200 m.
        let seg = parse_wind_segment_str("5:270:200", false).unwrap();
        assert!((seg.speed_kmh - 18.0).abs() < 1e-9);
        assert_eq!(seg.angle_deg, 270.0);
        assert!((seg.until_m - 200.0).abs() < 1e-9);

        // Malformed inputs are rejected.
        assert!(parse_wind_segment_str("10:90", true).is_err()); // too few fields
        assert!(parse_wind_segment_str("10:bad:100", true).is_err()); // non-numeric
        assert!(parse_wind_segment_str("10:90:0", true).is_err()); // zero until-distance
        assert!(parse_wind_segment_str("-3:90:100", true).is_err()); // negative speed
        // Non-finite values must be rejected (NaN comparisons would slip past < / <=).
        assert!(parse_wind_segment_str("10:nan:5000", true).is_err());
        assert!(parse_wind_segment_str("10:90:nan", true).is_err());
        assert!(parse_wind_segment_str("inf:90:100", true).is_err());
    }

    #[test]
    fn test_parse_wind_segment_str_vertical_field() {
        // MBA-728: 3-field input is unchanged (backward compat) -> vertical_mps == 0.0.
        let seg = parse_wind_segment_str("10:90:100", true).unwrap();
        assert_eq!(seg.vertical_mps, 0.0);
        let seg = parse_wind_segment_str("5:270:200", false).unwrap();
        assert_eq!(seg.vertical_mps, 0.0);

        // 4-field input parses the vertical component, m/s, regardless of --units.
        let seg = parse_wind_segment_str("10:90:100:5", true).unwrap();
        assert_eq!(seg.vertical_mps, 5.0);
        // The rest of the fields still go through the imperial conversion.
        assert!((seg.speed_kmh - 16.09344).abs() < 1e-4);
        assert!((seg.until_m - 91.44).abs() < 1e-4);

        // Vertical is NOT converted by --units (always m/s): metric and imperial parses of the
        // same "...:5" 4th field land on the identical vertical_mps.
        let seg_metric = parse_wind_segment_str("5:270:200:5", false).unwrap();
        assert_eq!(seg_metric.vertical_mps, 5.0);

        // Negative vertical (downdraft) is valid.
        let seg_neg = parse_wind_segment_str("10:90:100:-3.5", true).unwrap();
        assert_eq!(seg_neg.vertical_mps, -3.5);

        // A non-numeric or non-finite 4th field is rejected.
        assert!(parse_wind_segment_str("10:90:100:bad", true).is_err());
        assert!(parse_wind_segment_str("10:90:100:nan", true).is_err());
        assert!(parse_wind_segment_str("10:90:100:inf", true).is_err());

        // Too many fields is still rejected.
        assert!(parse_wind_segment_str("10:90:100:5:1", true).is_err());
    }
}
