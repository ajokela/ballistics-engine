use nalgebra::Vector3;
use std::collections::HashSet;

// Constants for unit conversions
const YARDS_TO_METERS: f64 = 0.9144;

/// Trajectory flags for notable events
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum TrajectoryFlag {
    ZeroCrossing,
    MachTransition,
    Apex,
}

impl TrajectoryFlag {
    pub fn to_string(&self) -> String {
        match self {
            TrajectoryFlag::ZeroCrossing => "zero_crossing".to_string(),
            TrajectoryFlag::MachTransition => "mach_transition".to_string(),
            TrajectoryFlag::Apex => "apex".to_string(),
        }
    }
}

/// Single trajectory sample point
#[derive(Debug, Clone)]
pub struct TrajectorySample {
    pub distance_m: f64,
    pub drop_m: f64,
    pub wind_drift_m: f64,
    pub velocity_mps: f64,
    pub energy_j: f64,
    pub time_s: f64,
    pub flags: Vec<TrajectoryFlag>,
}

/// Trajectory solution data for sampling
#[derive(Debug, Clone)]
pub struct TrajectoryData {
    pub times: Vec<f64>,
    pub positions: Vec<Vector3<f64>>,  // [x, y, z] positions
    pub velocities: Vec<Vector3<f64>>, // [vx, vy, vz] velocities
    pub transonic_distances: Vec<f64>, // Distances where mach transitions occur
}

/// Output data for trajectory sampling
#[derive(Debug, Clone)]
pub struct TrajectoryOutputs {
    pub target_distance_horiz_m: f64,
    pub target_vertical_height_m: f64,
    pub time_of_flight_s: f64,
    pub max_ord_dist_horiz_m: f64,
}

/// Sample trajectory at regular distance intervals with vectorized operations
pub fn sample_trajectory(
    trajectory_data: &TrajectoryData,
    outputs: &TrajectoryOutputs,
    step_m: f64,
    mass_kg: f64,
) -> Vec<TrajectorySample> {
    let step_size = if step_m <= 0.0 {
        return Vec::new();
    } else if step_m < 0.1 {
        0.1
    } else {
        step_m
    };
    
    let max_dist = outputs.target_distance_horiz_m;
    if max_dist < 1e-9 {
        return Vec::new();
    }
    
    // Extract trajectory arrays for vectorized operations
    let x_vals: Vec<f64> = trajectory_data.positions.iter().map(|p| p.x).collect();
    let y_vals: Vec<f64> = trajectory_data.positions.iter().map(|p| p.y).collect();
    let z_vals: Vec<f64> = trajectory_data.positions.iter().map(|p| p.z).collect();
    
    // Calculate speeds and energies
    let speeds: Vec<f64> = trajectory_data.velocities.iter()
        .map(|v| v.norm())
        .collect();
    let energies: Vec<f64> = speeds.iter()
        .map(|&speed| 0.5 * mass_kg * speed * speed)
        .collect();
    
    // Generate sampling distances
    let num_steps = ((max_dist / step_size) + 0.5) as usize + 1;
    let distances: Vec<f64> = (0..num_steps)
        .map(|i| i as f64 * step_size)
        .take_while(|&d| d <= max_dist + step_size * 0.5)
        .collect();
    
    // Vectorized interpolation for all trajectory data
    let mut samples = Vec::with_capacity(distances.len());
    
    for &distance in &distances {
        let y_interp = interpolate(&x_vals, &y_vals, distance);
        let wind_drift = interpolate(&x_vals, &z_vals, distance);
        let velocity = interpolate(&x_vals, &speeds, distance);
        let time = interpolate(&x_vals, &trajectory_data.times, distance);
        let energy = interpolate(&x_vals, &energies, distance);
        
        // Calculate line-of-sight y-coordinate and drop
        let los_y = outputs.target_vertical_height_m * distance / max_dist;
        let drop = los_y - y_interp;
        
        samples.push(TrajectorySample {
            distance_m: distance,
            drop_m: drop,
            wind_drift_m: wind_drift,
            velocity_mps: velocity,
            energy_j: energy,
            time_s: time,
            flags: Vec::new(), // Flags will be added later
        });
    }
    
    // Add flags using vectorized detection
    add_trajectory_flags(&mut samples, &trajectory_data.transonic_distances, outputs.max_ord_dist_horiz_m);
    
    samples
}

/// Linear interpolation function optimized for trajectory data
fn interpolate(x_vals: &[f64], y_vals: &[f64], x: f64) -> f64 {
    if x_vals.is_empty() || y_vals.is_empty() {
        return 0.0;
    }
    
    if x_vals.len() != y_vals.len() {
        return 0.0;
    }
    
    if x <= x_vals[0] {
        return y_vals[0];
    }
    
    if x >= x_vals[x_vals.len() - 1] {
        return y_vals[y_vals.len() - 1];
    }
    
    // Binary search for the correct interval
    let mut left = 0;
    let mut right = x_vals.len() - 1;
    
    while right - left > 1 {
        let mid = (left + right) / 2;
        if x_vals[mid] <= x {
            left = mid;
        } else {
            right = mid;
        }
    }
    
    // Linear interpolation
    let x1 = x_vals[left];
    let x2 = x_vals[right];
    let y1 = y_vals[left];
    let y2 = y_vals[right];
    
    if (x2 - x1).abs() < f64::EPSILON {
        return y1;
    }
    
    y1 + (y2 - y1) * (x - x1) / (x2 - x1)
}

/// Add trajectory flags using vectorized detection algorithms
fn add_trajectory_flags(
    samples: &mut [TrajectorySample],
    transonic_distances: &[f64],
    apex_distance: f64,
) {
    let tolerance = 1e-6;
    
    // 1. Zero crossings - vectorized detection
    detect_zero_crossings(samples, tolerance);
    
    // 2. Mach transitions
    for &transonic_dist in transonic_distances {
        if let Some(idx) = find_closest_sample_index(samples, transonic_dist) {
            samples[idx].flags.push(TrajectoryFlag::MachTransition);
        }
    }
    
    // 3. Apex
    if apex_distance > 0.0 {
        if let Some(idx) = find_closest_sample_index(samples, apex_distance) {
            samples[idx].flags.push(TrajectoryFlag::Apex);
        }
    }
}

/// Detect zero crossings in trajectory drop values using vectorized operations
fn detect_zero_crossings(samples: &mut [TrajectorySample], tolerance: f64) {
    if samples.len() < 2 {
        return;
    }
    
    let drops: Vec<f64> = samples.iter().map(|s| s.drop_m).collect();
    
    // Find crossing indices where drop changes sign
    for i in 0..(drops.len() - 1) {
        let current = drops[i];
        let next = drops[i + 1];
        
        // Check for sign change crossings
        let crosses_zero = (current < -tolerance && next >= -tolerance) ||
                          (current > tolerance && next <= tolerance);
        
        if crosses_zero {
            samples[i + 1].flags.push(TrajectoryFlag::ZeroCrossing);
        }
    }
    
    // Find points very close to zero
    for (i, &drop) in drops.iter().enumerate() {
        if drop.abs() <= tolerance {
            samples[i].flags.push(TrajectoryFlag::ZeroCrossing);
        }
    }
    
    // Remove duplicate zero crossing flags
    for sample in samples.iter_mut() {
        let mut unique_flags = Vec::new();
        let mut seen = HashSet::new();
        
        for flag in &sample.flags {
            if seen.insert(flag.clone()) {
                unique_flags.push(flag.clone());
            }
        }
        sample.flags = unique_flags;
    }
}

/// Find the closest sample index to a given distance
fn find_closest_sample_index(samples: &[TrajectorySample], target_distance: f64) -> Option<usize> {
    if samples.is_empty() {
        return None;
    }
    
    // Binary search for the closest distance
    let distances: Vec<f64> = samples.iter().map(|s| s.distance_m).collect();
    
    let mut left = 0;
    let mut right = distances.len();
    
    while left < right {
        let mid = (left + right) / 2;
        if distances[mid] < target_distance {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    
    // Find the closest point (could be left-1 or left)
    let mut best_idx = left.min(distances.len() - 1);
    
    if left > 0 {
        let left_dist = (distances[left - 1] - target_distance).abs();
        let right_dist = (distances[best_idx] - target_distance).abs();
        
        if left_dist < right_dist {
            best_idx = left - 1;
        }
    }
    
    Some(best_idx)
}

/// Convert trajectory samples to Python-compatible format
pub fn trajectory_samples_to_dicts(samples: &[TrajectorySample]) -> Vec<TrajectoryDict> {
    samples.iter().map(|sample| {
        TrajectoryDict {
            distance_m: sample.distance_m,
            drop_m: sample.drop_m,
            wind_drift_m: sample.wind_drift_m,
            velocity_mps: sample.velocity_mps,
            energy_j: sample.energy_j,
            time_s: sample.time_s,
            flags: sample.flags.iter().map(|f| f.to_string()).collect(),
        }
    }).collect()
}

/// Python-compatible trajectory sample structure
#[derive(Debug, Clone)]
pub struct TrajectoryDict {
    pub distance_m: f64,
    pub drop_m: f64,
    pub wind_drift_m: f64,
    pub velocity_mps: f64,
    pub energy_j: f64,
    pub time_s: f64,
    pub flags: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_interpolate() {
        let x_vals = vec![0.0, 1.0, 2.0, 3.0];
        let y_vals = vec![0.0, 10.0, 20.0, 30.0];
        
        assert_eq!(interpolate(&x_vals, &y_vals, 0.5), 5.0);
        assert_eq!(interpolate(&x_vals, &y_vals, 1.5), 15.0);
        assert_eq!(interpolate(&x_vals, &y_vals, 2.5), 25.0);
        
        // Test boundary conditions
        assert_eq!(interpolate(&x_vals, &y_vals, -1.0), 0.0);  // Below range
        assert_eq!(interpolate(&x_vals, &y_vals, 4.0), 30.0);  // Above range
    }
    
    #[test]
    fn test_find_closest_sample_index() {
        let samples = vec![
            TrajectorySample {
                distance_m: 0.0,
                drop_m: 0.0,
                wind_drift_m: 0.0,
                velocity_mps: 100.0,
                energy_j: 1000.0,
                time_s: 0.0,
                flags: Vec::new(),
            },
            TrajectorySample {
                distance_m: 10.0,
                drop_m: -1.0,
                wind_drift_m: 0.1,
                velocity_mps: 95.0,
                energy_j: 950.0,
                time_s: 0.1,
                flags: Vec::new(),
            },
            TrajectorySample {
                distance_m: 20.0,
                drop_m: -4.0,
                wind_drift_m: 0.2,
                velocity_mps: 90.0,
                energy_j: 900.0,
                time_s: 0.2,
                flags: Vec::new(),
            },
        ];
        
        assert_eq!(find_closest_sample_index(&samples, 5.0), Some(0));
        assert_eq!(find_closest_sample_index(&samples, 12.0), Some(1));
        assert_eq!(find_closest_sample_index(&samples, 18.0), Some(2));
    }
    
    #[test]
    fn test_detect_zero_crossings() {
        let mut samples = vec![
            TrajectorySample {
                distance_m: 0.0,
                drop_m: 1.0,  // Positive
                wind_drift_m: 0.0,
                velocity_mps: 100.0,
                energy_j: 1000.0,
                time_s: 0.0,
                flags: Vec::new(),
            },
            TrajectorySample {
                distance_m: 10.0,
                drop_m: -0.5,  // Negative - crossing here
                wind_drift_m: 0.1,
                velocity_mps: 95.0,
                energy_j: 950.0,
                time_s: 0.1,
                flags: Vec::new(),
            },
            TrajectorySample {
                distance_m: 20.0,
                drop_m: -2.0,  // Still negative
                wind_drift_m: 0.2,
                velocity_mps: 90.0,
                energy_j: 900.0,
                time_s: 0.2,
                flags: Vec::new(),
            },
        ];
        
        detect_zero_crossings(&mut samples, 1e-6);
        
        // Should have a zero crossing flag at index 1
        assert!(!samples[0].flags.contains(&TrajectoryFlag::ZeroCrossing));
        assert!(samples[1].flags.contains(&TrajectoryFlag::ZeroCrossing));
        assert!(!samples[2].flags.contains(&TrajectoryFlag::ZeroCrossing));
    }
    
    #[test]
    fn test_sample_trajectory_basic() {
        // Create simple test trajectory data
        let trajectory_data = TrajectoryData {
            times: vec![0.0, 1.0, 2.0],
            positions: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(100.0, 10.0, 1.0),
                Vector3::new(200.0, 5.0, 2.0),
            ],
            velocities: vec![
                Vector3::new(100.0, 10.0, 1.0),
                Vector3::new(95.0, 5.0, 1.0),
                Vector3::new(90.0, 0.0, 1.0),
            ],
            transonic_distances: vec![150.0],
        };
        
        let outputs = TrajectoryOutputs {
            target_distance_horiz_m: 200.0,
            target_vertical_height_m: 0.0,
            time_of_flight_s: 2.0,
            max_ord_dist_horiz_m: 100.0,
        };
        
        let samples = sample_trajectory(&trajectory_data, &outputs, 50.0, 0.1);
        
        // Should have samples at 0, 50, 100, 150, 200 meters
        assert_eq!(samples.len(), 5);
        assert_eq!(samples[0].distance_m, 0.0);
        assert_eq!(samples[1].distance_m, 50.0);
        assert_eq!(samples[2].distance_m, 100.0);
        assert_eq!(samples[3].distance_m, 150.0);
        assert_eq!(samples[4].distance_m, 200.0);
        
        // Check that interpolation is working
        assert!(samples[1].velocity_mps > 90.0 && samples[1].velocity_mps < 100.0);
        
        // Check flags
        assert!(samples[2].flags.contains(&TrajectoryFlag::Apex)); // At apex distance
        assert!(samples[3].flags.contains(&TrajectoryFlag::MachTransition)); // At transonic distance
    }
}