use nalgebra::Vector3;
use std::cmp::Ordering;
use std::f64::consts::PI;

/// Conversion constant from KMH to MPS
const KMH_TO_MPS: f64 = 1000.0 / 3600.0;

/// Wind segment: (speed_kmh, angle_deg, until_distance_m)
/// This matches the Python WindSock interface
pub type WindSegment = (f64, f64, f64);

/// Sort wind segments by their `until_distance_m` threshold.
///
/// Shared by [`WindSock`] and the low-level trajectory integrator so every segmented-wind path
/// applies the same interval ordering.
pub(crate) fn sort_wind_segments_by_distance(segments: &mut [WindSegment]) {
    segments.sort_by(|a, b| match (a.2.is_nan(), b.2.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => {
            a.2.partial_cmp(&b.2)
                .expect("non-NaN distances are ordered")
        }
    });
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
}

impl WindSock {
    /// Create a new WindSock from wind segments
    ///
    /// Args:
    ///     segments: List of (speed_kmh, angle_deg, until_distance_m) tuples
    pub fn new(mut segments: Vec<WindSegment>) -> Self {
        // Sort segments by distance, handling NaN safely by treating it as greater than any value
        sort_wind_segments_by_distance(&mut segments);

        // Precompute each segment's wind vector once (depends only on its speed/angle).
        let wind_vecs: Vec<Vector3<f64>> = segments.iter().map(Self::calc_vec).collect();

        let (current, next_range, current_vec) = if segments.is_empty() {
            (0, f64::INFINITY, Vector3::zeros())
        } else {
            (0, segments[0].2, wind_vecs[0])
        };

        WindSock {
            winds: segments,
            wind_vecs,
            current,
            next_range,
            current_vec,
        }
    }

    /// Calculate wind vector from wind segment
    fn calc_vec(seg: &WindSegment) -> Vector3<f64> {
        let (speed_kmh, angle_deg, _) = *seg;

        // Convert kmh to m/s
        let speed_mps = speed_kmh * KMH_TO_MPS;
        let angle_rad = angle_deg * PI / 180.0;

        // Wind convention (matching trajectory coordinates):
        // 0° = headwind (from front, affects -x downrange)
        // 90° = wind from right (affects -z lateral)
        // 180° = tailwind (from back, affects +x downrange)
        // 270° = wind from left (affects +z lateral)
        //
        // McCoy convention: x=downrange, y=vertical, z=lateral
        Vector3::new(
            -speed_mps * angle_rad.cos(), // x (downrange - head/tail component)
            0.0,                          // y (vertical)
            -speed_mps * angle_rad.sin(), // z (lateral - crosswind component)
        )
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
                self.next_range = self.winds[self.current].2;
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
            if range_m < segment.2 {
                return self.wind_vecs[i];
            }
        }

        // Beyond all segments
        Vector3::zeros()
    }
}

/// Parse a `"SPEED:ANGLE:UNTIL_DISTANCE"` string into a [`WindSegment`]
/// `(speed_kmh, angle_deg, until_distance_m)`.
///
/// `imperial`: when true, SPEED is mph and UNTIL_DISTANCE is yards; otherwise
/// SPEED is m/s and UNTIL_DISTANCE is meters. ANGLE is always degrees in the
/// wind-FROM convention (0 = headwind, 90 = from the right). Shared by the CLI
/// (`--wind-segment`) and the WASM front-ends so they parse identically.
pub fn parse_wind_segment_str(s: &str, imperial: bool) -> Result<WindSegment, String> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 3 {
        return Err(format!(
            "invalid wind segment '{s}': expected SPEED:ANGLE:UNTIL_DISTANCE (three colon-separated numbers)"
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
    if !speed.is_finite() || !angle.is_finite() || !until.is_finite() {
        return Err(format!(
            "invalid wind segment '{s}': speed, angle, and until-distance must be finite numbers"
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
    Ok((speed_kmh, angle, until_m))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn segment_sort_is_stable_and_places_nan_endpoints_last() {
        let mut segments = vec![
            (10.0, 0.0, f64::NAN),
            (20.0, 0.0, 100.0),
            (30.0, 0.0, 100.0),
            (40.0, 0.0, f64::INFINITY),
            (50.0, 0.0, f64::NEG_INFINITY),
            (60.0, 0.0, f64::NAN),
        ];

        sort_wind_segments_by_distance(&mut segments);

        assert_eq!(segments[0].0, 50.0); // -inf first
        assert_eq!(segments[1].0, 20.0); // equal endpoints retain input order
        assert_eq!(segments[2].0, 30.0);
        assert_eq!(segments[3].0, 40.0); // +inf after finite endpoints
        assert_eq!(segments[4].0, 10.0); // NaNs last and stable
        assert_eq!(segments[5].0, 60.0);
        assert!(segments[4].2.is_nan() && segments[5].2.is_nan());
    }

    #[test]
    fn test_wind_sock_empty() {
        let sock = WindSock::new(vec![]);
        assert_eq!(sock.vector_for_range_stateless(50.0), Vector3::zeros());
    }

    #[test]
    fn test_wind_sock_single_segment() {
        // 16.0934 kmh (10 mph) @ 90° until 100m
        let sock = WindSock::new(vec![(16.0934, 90.0, 100.0)]);

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
        assert_eq!(vec_50[1], 0.0); // Zero Y component
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
            (16.0934, 90.0, 50.0),  // 10 mph @ 90° until 50m
            (24.1401, 45.0, 100.0), // 15 mph @ 45° until 100m
            (8.0467, 180.0, 200.0), // 5 mph @ 180° until 200m
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
        let sock = WindSock::new(vec![(16.0934, 0.0, 100.0)]);
        let vec = sock.vector_for_range_stateless(50.0);

        let expected_speed = 16.0934 * KMH_TO_MPS;
        assert!((vec.norm() - expected_speed).abs() < 0.01);
    }

    #[test]
    fn test_wind_sock_boundary_is_upper_exclusive() {
        // A segment's `until_distance_m` is exclusive: a query exactly at the
        // boundary rolls to the next segment.
        let sock = WindSock::new(vec![(16.0934, 90.0, 100.0), (32.1868, 270.0, 200.0)]);
        // Just below 100 m -> first segment (90deg, negative Z).
        assert!(sock.vector_for_range_stateless(99.999)[2] < 0.0);
        // Exactly 100 m -> second segment (270deg, positive Z).
        assert!(sock.vector_for_range_stateless(100.0)[2] > 0.0);
        // Beyond the last boundary -> zero.
        assert_eq!(sock.vector_for_range_stateless(200.0), Vector3::zeros());
    }

    #[test]
    fn test_parse_wind_segment_str_units() {
        // Imperial: 10 mph -> 16.0934 km/h, 100 yd -> 91.44 m.
        let (kmh, ang, until) = parse_wind_segment_str("10:90:100", true).unwrap();
        assert!((kmh - 16.09344).abs() < 1e-4);
        assert_eq!(ang, 90.0);
        assert!((until - 91.44).abs() < 1e-4);

        // Metric: 5 m/s -> 18 km/h, 200 m stays 200 m.
        let (kmh, ang, until) = parse_wind_segment_str("5:270:200", false).unwrap();
        assert!((kmh - 18.0).abs() < 1e-9);
        assert_eq!(ang, 270.0);
        assert!((until - 200.0).abs() < 1e-9);

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
}
