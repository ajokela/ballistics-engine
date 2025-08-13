use nalgebra::Vector3;
use std::f64::consts::PI;

/// Conversion constant from KMH to MPS
const KMH_TO_MPS: f64 = 1000.0 / 3600.0;

/// Wind segment: (speed_kmh, angle_deg, until_distance_m)
/// This matches the Python WindSock interface
pub type WindSegment = (f64, f64, f64);

/// Wind condition handler for trajectory calculations
#[derive(Debug, Clone)]
pub struct WindSock {
    /// Sorted wind segments by distance
    winds: Vec<WindSegment>,
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
        // Sort segments by distance
        segments.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        
        let (current, next_range, current_vec) = if segments.is_empty() {
            (0, f64::INFINITY, Vector3::zeros())
        } else {
            let vec = Self::calc_vec(&segments[0]);
            (0, segments[0].2, vec)
        };
        
        WindSock {
            winds: segments,
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
        
        // Wind vector points in the direction the wind is blowing TO
        // So a 90° wind blows from right to left (negative Z)
        Vector3::new(
            -speed_mps * angle_rad.cos(),
            0.0,
            -speed_mps * angle_rad.sin(),
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
        
        // Check if we need to advance to next segment
        if range_m >= self.next_range {
            self.current += 1;
            if self.current >= self.winds.len() {
                self.current_vec = Vector3::zeros();
                self.next_range = f64::INFINITY;
            } else {
                let seg = &self.winds[self.current];
                self.current_vec = Self::calc_vec(seg);
                self.next_range = seg.2;
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
        
        // Find the appropriate segment
        for segment in &self.winds {
            if range_m < segment.2 {
                return Self::calc_vec(segment);
            }
        }
        
        // Beyond all segments
        Vector3::zeros()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
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
        assert!(vec_50.norm() > 0.0);
        assert!(vec_50[0].abs() < 0.01); // Nearly zero X component
        assert_eq!(vec_50[1], 0.0); // Zero Y component
        assert!(vec_50[2] < 0.0); // Negative Z (wind from right)
        
        // No wind after 100m
        let vec_150 = sock.vector_for_range_stateless(150.0);
        assert_eq!(vec_150, Vector3::zeros());
    }
    
    #[test]
    fn test_wind_sock_multiple_segments() {
        // Multiple wind segments (in kmh)
        let sock = WindSock::new(vec![
            (16.0934, 90.0, 50.0),   // 10 mph @ 90° until 50m
            (24.1401, 45.0, 100.0),  // 15 mph @ 45° until 100m
            (8.0467, 180.0, 200.0),  // 5 mph @ 180° until 200m
        ]);
        
        // Test each segment
        let vec_25 = sock.vector_for_range_stateless(25.0);
        assert!(vec_25.norm() > 0.0);
        assert!(vec_25[2] < 0.0); // 90° wind
        
        let vec_75 = sock.vector_for_range_stateless(75.0);
        assert!(vec_75.norm() > vec_25.norm()); // 15 mph > 10 mph
        assert!(vec_75[0] < 0.0); // 45° wind has X component
        assert!(vec_75[2] < 0.0); // 45° wind has Z component
        
        let vec_150 = sock.vector_for_range_stateless(150.0);
        assert!(vec_150.norm() < vec_75.norm()); // 5 mph < 15 mph
        assert!(vec_150[0] > 0.0); // 180° wind (from behind)
        
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
}