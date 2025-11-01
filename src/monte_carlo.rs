//! Monte Carlo simulation support with statistical analysis
//!
//! This module provides core statistical functions for Monte Carlo trajectory analysis,
//! including CEP (Circular Error Probable), confidence ellipses, and trajectory evaluation.
//!
//! MBA-157: Upstreamed from ballistics_rust for shared use across the ecosystem

use crate::BallisticInputs;
use crate::fast_trajectory::{fast_integrate, FastIntegrationParams};
use crate::wind::WindSock;
use crate::atmosphere::calculate_atmosphere;
use crate::constants::{FPS_TO_MPS, GRAINS_TO_KG};
use nalgebra::Vector3;

const YARDS_TO_METERS: f64 = 0.9144;

/// Simple trajectory output for Monte Carlo analysis
#[derive(Debug, Clone)]
pub struct TrajectoryOutput {
    pub drop: f64,           // meters
    pub wind_drift: f64,     // meters
    pub time: f64,           // seconds
    pub velocity: f64,       // m/s
    pub energy: f64,         // joules
    pub mach: f64,           // mach number
    pub spin_drift: f64,     // meters
    pub distance: f64,       // meters
}

/// Solve trajectory for Monte Carlo run
///
/// This function evaluates a single trajectory with the given inputs and returns
/// simplified output suitable for statistical analysis.
pub fn solve_trajectory_for_monte_carlo(inputs: &BallisticInputs) -> Result<TrajectoryOutput, String> {
    // Convert inputs to metric
    let target_distance_m = inputs.target_distance * YARDS_TO_METERS;
    let muzzle_velocity_mps = inputs.muzzle_velocity * FPS_TO_MPS;
    let mass_kg = inputs.bullet_mass * GRAINS_TO_KG;

    // Calculate atmosphere at altitude
    let (air_density, speed_of_sound) = calculate_atmosphere(
        inputs.altitude * 0.3048,  // feet to meters
        Some(inputs.temperature),
        Some(inputs.pressure),
        inputs.humidity,
    );

    // Create wind segments
    let wind_segments = vec![(
        0.0,  // range_m
        inputs.wind_speed * 0.44704,  // wind_speed (mph to m/s)
        inputs.wind_angle,  // wind_angle
    )];
    let wind_sock = WindSock::new(wind_segments);

    // Set up initial state
    let muzzle_angle_rad = inputs.muzzle_angle;
    let initial_velocity = Vector3::new(
        0.0,
        muzzle_velocity_mps * muzzle_angle_rad.sin(),
        muzzle_velocity_mps * muzzle_angle_rad.cos(),
    );

    let initial_position = Vector3::new(0.0, inputs.sight_height * 0.0254, 0.0);
    let mut initial_state_array = [0.0; 6];
    initial_state_array[0..3].copy_from_slice(&[initial_position.x, initial_position.y, initial_position.z]);
    initial_state_array[3..6].copy_from_slice(&[initial_velocity.x, initial_velocity.y, initial_velocity.z]);

    // Get atmospheric parameters (temperature, pressure, density, sound_speed)
    let temp_c = inputs.temperature;
    let pressure_hpa = inputs.pressure;

    // Create integration params
    let params = FastIntegrationParams {
        initial_state: initial_state_array,
        t_span: (0.0, 30.0),
        horiz: target_distance_m,
        vert: 0.0,  // Target at ground level
        atmo_params: (temp_c, pressure_hpa, air_density, speed_of_sound),
    };

    // Solve trajectory
    let solution = fast_integrate(inputs, &wind_sock, params);

    if solution.t.is_empty() {
        return Err("Empty trajectory solution".to_string());
    }

    // Get final state
    // FastSolution.y is Vec<Vec<f64>> where y[i] is the ith state variable across all time points
    let final_idx = solution.t.len() - 1;

    let final_x = solution.y[0][final_idx];  // lateral drift
    let final_y = solution.y[1][final_idx];  // vertical
    let final_z = solution.y[2][final_idx];  // downrange

    let final_vx = solution.y[3][final_idx];
    let final_vy = solution.y[4][final_idx];
    let final_vz = solution.y[5][final_idx];

    let final_speed = (final_vx * final_vx + final_vy * final_vy + final_vz * final_vz).sqrt();
    let final_mach = final_speed / speed_of_sound;
    let final_energy = 0.5 * mass_kg * final_speed * final_speed;

    // Calculate line-of-sight drop
    let sight_height_m = inputs.sight_height * 0.0254;
    let los_y = sight_height_m + (0.0 - sight_height_m) * (final_z / target_distance_m);
    let drop = los_y - final_y;

    Ok(TrajectoryOutput {
        drop,
        wind_drift: final_x,
        time: solution.t[final_idx],
        velocity: final_speed,
        energy: final_energy,
        mach: final_mach,
        spin_drift: final_x,  // Approximation for now
        distance: final_z,
    })
}

/// Calculate CEP (Circular Error Probable) from impact points
///
/// CEP is the radius of a circle centered at the mean point of impact,
/// within which 50% of the shots fall. It's a standard measure of precision.
pub fn calculate_cep(wind_drift_values: &[f64], drop_values: &[f64]) -> f64 {
    if wind_drift_values.len() != drop_values.len() || wind_drift_values.is_empty() {
        return 0.0;
    }

    // Calculate mean point of impact
    let mean_x = wind_drift_values.iter().sum::<f64>() / wind_drift_values.len() as f64;
    let mean_y = drop_values.iter().sum::<f64>() / drop_values.len() as f64;

    // Calculate distance from each point to mean
    let mut distances: Vec<f64> = wind_drift_values.iter()
        .zip(drop_values.iter())
        .map(|(x, y)| {
            let dx = x - mean_x;
            let dy = y - mean_y;
            (dx * dx + dy * dy).sqrt()
        })
        .collect();

    // Sort distances to find median (50th percentile)
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // CEP is the median distance from center
    percentile(&distances, 0.50)
}

/// Calculate 95% confidence ellipse parameters using covariance matrix
///
/// Returns (center_x, center_y, semi_major_axis, semi_minor_axis, rotation_degrees)
pub fn calculate_confidence_ellipse(wind_drift_values: &[f64], drop_values: &[f64]) -> (f64, f64, f64, f64, f64) {
    if wind_drift_values.len() != drop_values.len() || wind_drift_values.len() < 2 {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }

    let n = wind_drift_values.len() as f64;

    // Calculate means
    let mean_x = wind_drift_values.iter().sum::<f64>() / n;
    let mean_y = drop_values.iter().sum::<f64>() / n;

    // Calculate covariance matrix elements
    let mut cov_xx = 0.0;
    let mut cov_yy = 0.0;
    let mut cov_xy = 0.0;

    for (x, y) in wind_drift_values.iter().zip(drop_values.iter()) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        cov_xx += dx * dx;
        cov_yy += dy * dy;
        cov_xy += dx * dy;
    }

    cov_xx /= n - 1.0;
    cov_yy /= n - 1.0;
    cov_xy /= n - 1.0;

    // Calculate eigenvalues of covariance matrix
    // For 2x2 matrix: [[cov_xx, cov_xy], [cov_xy, cov_yy]]
    let trace = cov_xx + cov_yy;
    let det = cov_xx * cov_yy - cov_xy * cov_xy;
    let discriminant = (trace * trace / 4.0 - det).max(0.0).sqrt();

    let lambda1 = trace / 2.0 + discriminant;  // Larger eigenvalue
    let lambda2 = trace / 2.0 - discriminant;  // Smaller eigenvalue

    // 95% confidence interval chi-square value for 2 DOF is 5.991
    let scale_factor = 5.991_f64.sqrt();
    let semi_major = lambda1.max(0.0).sqrt() * scale_factor;
    let semi_minor = lambda2.max(0.0).sqrt() * scale_factor;

    // Calculate rotation angle (angle of major axis)
    let rotation_rad = if cov_xy.abs() < 1e-10 {
        if cov_xx >= cov_yy { 0.0 } else { std::f64::consts::PI / 2.0 }
    } else {
        ((lambda1 - cov_xx) / cov_xy).atan()
    };

    let rotation_deg = rotation_rad.to_degrees();

    (mean_x, mean_y, semi_major, semi_minor, rotation_deg)
}

/// Sample points for visualization (limit to avoid huge payloads)
pub fn sample_points_for_visualization(
    wind_drift_values: &[f64],
    drop_values: &[f64],
    max_points: usize
) -> Vec<(f64, f64)> {
    let n = wind_drift_values.len();
    if n == 0 {
        return Vec::new();
    }

    if n <= max_points {
        // Return all points
        wind_drift_values.iter()
            .zip(drop_values.iter())
            .map(|(x, y)| (*x, *y))
            .collect()
    } else {
        // Sample evenly spaced points
        let step = n as f64 / max_points as f64;
        (0..max_points)
            .map(|i| {
                let idx = (i as f64 * step) as usize;
                (wind_drift_values[idx], drop_values[idx])
            })
            .collect()
    }
}

/// Calculate percentile from sorted values
pub fn percentile(sorted_values: &[f64], p: f64) -> f64 {
    if sorted_values.is_empty() {
        return 0.0;
    }

    if sorted_values.len() == 1 {
        return sorted_values[0];
    }

    let rank = p * (sorted_values.len() - 1) as f64;
    let lower_idx = rank.floor() as usize;
    let upper_idx = rank.ceil() as usize;
    let fraction = rank - lower_idx as f64;

    if lower_idx == upper_idx {
        sorted_values[lower_idx]
    } else {
        sorted_values[lower_idx] * (1.0 - fraction) + sorted_values[upper_idx] * fraction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_cep() {
        let wind_drift = vec![0.0, 1.0, -1.0, 0.5, -0.5];
        let drop = vec![0.0, 0.5, -0.5, 1.0, -1.0];

        let cep = calculate_cep(&wind_drift, &drop);
        assert!(cep > 0.0);
        assert!(cep < 2.0); // Reasonable range
    }

    #[test]
    fn test_calculate_confidence_ellipse() {
        let wind_drift = vec![0.0, 1.0, -1.0, 0.5, -0.5];
        let drop = vec![0.0, 0.5, -0.5, 1.0, -1.0];

        let (cx, cy, major, minor, _rotation) = calculate_confidence_ellipse(&wind_drift, &drop);

        // Center should be near origin
        assert!(cx.abs() < 0.5);
        assert!(cy.abs() < 0.5);

        // Axes should be positive
        assert!(major > 0.0);
        assert!(minor > 0.0);
        assert!(major >= minor); // Major axis should be >= minor axis
    }

    #[test]
    fn test_sample_points() {
        let wind_drift = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let drop = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5];

        let sampled = sample_points_for_visualization(&wind_drift, &drop, 3);
        assert_eq!(sampled.len(), 3);
    }

    #[test]
    fn test_percentile() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        assert_eq!(percentile(&values, 0.0), 1.0);
        assert_eq!(percentile(&values, 0.5), 3.0);
        assert_eq!(percentile(&values, 1.0), 5.0);
    }
}
