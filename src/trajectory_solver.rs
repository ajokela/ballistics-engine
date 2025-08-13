use crate::BallisticInputs;
use crate::constants::{FPS_TO_MPS, MPS_TO_FPS, GRAINS_TO_KG};
use nalgebra::Vector3;

// Constants for unit conversions
const YARDS_TO_METERS: f64 = 0.9144;
const JOULES_TO_FTLBS: f64 = 0.737562149;
const METERS_TO_INCHES: f64 = 39.3701;

/// Initial conditions for trajectory solving
#[derive(Debug, Clone)]
pub struct InitialConditions {
    pub mass_kg: f64,
    pub muzzle_velocity_mps: f64,
    pub target_distance_los_m: f64,
    pub muzzle_angle_rad: f64,
    pub muzzle_energy_j: f64,
    pub muzzle_energy_ftlbs: f64,
    pub target_horizontal_dist_m: f64,
    pub target_vertical_height_m: f64,
    pub initial_state: [f64; 6],  // [x, y, z, vx, vy, vz]
    pub t_span: (f64, f64),
    pub omega_vector: Option<Vector3<f64>>,
    pub stability_coefficient: f64,
    pub atmo_params: (f64, f64, f64, f64),  // altitude, temp_c, pressure_hpa, density_ratio
    pub air_density: f64,
    pub speed_of_sound: f64,
}

/// Result of trajectory post-processing
#[derive(Debug, Clone)]
pub struct TrajectoryResult {
    pub muzzle_energy_j: f64,
    pub muzzle_energy_ftlbs: f64,
    pub target_distance_los_m: f64,
    pub target_distance_horiz_m: f64,
    pub target_vertical_height_m: f64,
    pub time_of_flight_s: f64,
    pub drop_m: f64,
    pub drop_in: f64,
    pub wind_drift_m: f64,
    pub wind_drift_in: f64,
    pub max_ord_m: f64,
    pub max_ord_in: f64,
    pub max_ord_dist_horiz_m: f64,
    pub final_vel_mps: f64,
    pub final_vel_fps: f64,
    pub final_energy_j: f64,
    pub final_energy_ftlbs: f64,
    pub air_density_kg_m3: f64,
    pub speed_of_sound_mps: f64,
    pub barrel_angle_rad: f64,
}

/// Prepare initial conditions for trajectory solving
pub fn prepare_initial_conditions(
    inputs: &BallisticInputs,
    zero_angle_rad: f64,
    atmo_params: (f64, f64, f64, f64),
    air_density: f64,
    speed_of_sound: f64,
    stability_coefficient: f64,
) -> InitialConditions {
    let mass_kg = inputs.bullet_mass * GRAINS_TO_KG;
    
    // Adjust muzzle velocity (basic implementation - could be enhanced)
    let mv_mps = inputs.muzzle_velocity * FPS_TO_MPS;
    
    // Calculate target coordinates
    let target_dist_m_los = inputs.target_distance * YARDS_TO_METERS;
    let target_horizontal_dist_m = target_dist_m_los;
    let target_vertical_height_m = 0.0; // Simplified for now
    
    // Calculate muzzle angle
    let muzzle_angle_rad = zero_angle_rad + inputs.muzzle_angle.to_radians();
    
    // Calculate energies
    let muzzle_energy_j = 0.5 * mass_kg * mv_mps * mv_mps;
    let muzzle_energy_ftlbs = muzzle_energy_j * JOULES_TO_FTLBS;
    
    // Set up initial velocity vector
    let initial_vel = Vector3::new(
        mv_mps * muzzle_angle_rad.cos(),
        mv_mps * muzzle_angle_rad.sin(),
        0.0,
    );
    
    // Initial state: [x, y, z, vx, vy, vz]
    let initial_state = [
        0.0, 0.0, 0.0,
        initial_vel.x, initial_vel.y, initial_vel.z,
    ];
    
    // Estimate maximum time
    let initial_vx = initial_vel.x;
    let max_time = if initial_vx > 1e-6 && target_horizontal_dist_m > 0.0 {
        let est_min = target_horizontal_dist_m / initial_vx;
        (est_min * 3.0).max(10.0)
    } else {
        10.0
    };
    let t_span = (0.0, max_time);
    
    // Omega vector for Coriolis (simplified - would need more complex calculation)
    let omega_vector = if inputs.enable_advanced_effects {
        // Simplified Coriolis vector calculation
        let latitude_rad = inputs.latitude.unwrap_or(0.0).to_radians();
        let earth_rotation_rate = 7.2921159e-5; // rad/s
        Some(Vector3::new(
            0.0,
            earth_rotation_rate * latitude_rad.cos(),
            earth_rotation_rate * latitude_rad.sin(),
        ))
    } else {
        None
    };
    
    InitialConditions {
        mass_kg,
        muzzle_velocity_mps: mv_mps,
        target_distance_los_m: target_dist_m_los,
        muzzle_angle_rad,
        muzzle_energy_j,
        muzzle_energy_ftlbs,
        target_horizontal_dist_m,
        target_vertical_height_m,
        initial_state,
        t_span,
        omega_vector,
        stability_coefficient,
        atmo_params,
        air_density,
        speed_of_sound,
    }
}

/// Find trajectory apex using Brent's method root finding
pub fn find_trajectory_apex(
    trajectory_points: &[(f64, [f64; 6])], // (time, state) pairs
    target_horizontal_dist_m: f64,
) -> (f64, f64) { // (max_ordinate_m, max_ordinate_x_m)
    let mut max_ord_m = 0.0;
    let mut max_ord_x_m = 0.0;
    
    // Find the highest point that occurs before the target
    for &(_, state) in trajectory_points {
        let x = state[0];
        let y = state[1];
        
        if x <= target_horizontal_dist_m + 1e-6 && y > max_ord_m {
            max_ord_m = y;
            max_ord_x_m = x;
        }
    }
    
    (max_ord_m, max_ord_x_m)
}

/// Brent's method for root finding (simplified implementation)
pub fn brent_root_find<F>(
    f: F,
    mut a: f64,
    mut b: f64,
    tolerance: f64,
    max_iterations: usize,
) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let mut fa = f(a);
    let mut fb = f(b);
    
    // Ensure the root is bracketed
    if fa * fb > 0.0 {
        return Err("Root not bracketed".to_string());
    }
    
    // Ensure |f(a)| >= |f(b)|
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;
    
    for _ in 0..max_iterations {
        if fb.abs() < tolerance {
            return Ok(b);
        }
        
        if fa.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        
        let tolerance_scaled = 2.0 * f64::EPSILON * b.abs() + 0.5 * tolerance;
        let m = 0.5 * (c - b);
        
        if m.abs() <= tolerance_scaled {
            return Ok(b);
        }
        
        if e.abs() >= tolerance_scaled && fc.abs() > fb.abs() {
            let s = fb / fc;
            let mut p;
            let mut q;
            
            if (a - c).abs() < f64::EPSILON {
                // Linear interpolation
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation
                q = fc / fa;
                let r = fb / fa;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            
            let s = e;
            e = d;
            
            if 2.0 * p < 3.0 * m * q - (tolerance_scaled * q).abs() 
                && p < (0.5 * s * q).abs() {
                d = p / q;
            } else {
                d = m;
                e = d;
            }
        } else {
            d = m;
            e = d;
        }
        
        a = b;
        fa = fb;
        
        if d.abs() > tolerance_scaled {
            b += d;
        } else if m > 0.0 {
            b += tolerance_scaled;
        } else {
            b -= tolerance_scaled;
        }
        
        fb = f(b);
        
        if (fc * fb) > 0.0 {
            c = a;
            fc = fa;
            e = b - a;
            d = e;
        }
    }
    
    Err("Maximum iterations exceeded".to_string())
}

/// Post-process trajectory solution to create final results
pub fn post_process_trajectory(
    trajectory_points: &[(f64, [f64; 6])], // (time, state) pairs
    initial_conditions: &InitialConditions,
    inputs: &BallisticInputs,
    target_hit_time: Option<f64>,
    ground_hit_time: Option<f64>,
) -> Result<TrajectoryResult, String> {
    // Determine final time and state
    let (final_time, final_state) = if let Some(hit_time) = target_hit_time {
        // Interpolate state at target hit time
        let final_state = interpolate_trajectory_state(trajectory_points, hit_time)?;
        (hit_time, final_state)
    } else if ground_hit_time.is_some() {
        return Err("Projectile impacted ground before reaching target".to_string());
    } else if let Some((time, state)) = trajectory_points.last() {
        (*time, *state)
    } else {
        return Err("No trajectory data available".to_string());
    };
    
    // Check if target was reached
    let distance_err = final_state[0] - initial_conditions.target_horizontal_dist_m;
    if distance_err.abs() >= 1e-3 && final_state[0] < initial_conditions.target_horizontal_dist_m - 1e-3 {
        return Err(format!(
            "Target horizontal distance ({:.2}m) not reached. Max distance: {:.2}m",
            initial_conditions.target_horizontal_dist_m, final_state[0]
        ));
    }
    
    // Find trajectory apex
    let (max_ord_m, max_ord_x_m) = find_trajectory_apex(
        trajectory_points,
        initial_conditions.target_horizontal_dist_m,
    );
    
    // Calculate final results
    let final_y_m = final_state[1];
    let final_z_m = final_state[2];
    let drop_m = initial_conditions.target_vertical_height_m - final_y_m;
    
    // Calculate wind drift including spin drift
    let mut wind_drift_m = final_z_m;
    if inputs.enable_advanced_effects {
        // Add spin drift using existing function
        use crate::stability::compute_spin_drift;
        wind_drift_m += compute_spin_drift(
            final_time,
            initial_conditions.stability_coefficient,
            inputs.twist_rate,
            inputs.is_twist_right,
        );
    }
    
    // Calculate final velocity and energy
    let final_vel = Vector3::new(final_state[3], final_state[4], final_state[5]);
    let final_vel_mag = final_vel.norm();
    let final_energy_j = 0.5 * initial_conditions.mass_kg * final_vel_mag * final_vel_mag;
    
    Ok(TrajectoryResult {
        muzzle_energy_j: initial_conditions.muzzle_energy_j,
        muzzle_energy_ftlbs: initial_conditions.muzzle_energy_ftlbs,
        target_distance_los_m: initial_conditions.target_distance_los_m,
        target_distance_horiz_m: initial_conditions.target_horizontal_dist_m,
        target_vertical_height_m: initial_conditions.target_vertical_height_m,
        time_of_flight_s: final_time,
        drop_m,
        drop_in: drop_m * METERS_TO_INCHES,
        wind_drift_m,
        wind_drift_in: wind_drift_m * METERS_TO_INCHES,
        max_ord_m,
        max_ord_in: max_ord_m * METERS_TO_INCHES,
        max_ord_dist_horiz_m: max_ord_x_m,
        final_vel_mps: final_vel_mag,
        final_vel_fps: final_vel_mag * MPS_TO_FPS,
        final_energy_j,
        final_energy_ftlbs: final_energy_j * JOULES_TO_FTLBS,
        air_density_kg_m3: initial_conditions.air_density,
        speed_of_sound_mps: initial_conditions.speed_of_sound,
        barrel_angle_rad: initial_conditions.muzzle_angle_rad,
    })
}

/// Interpolate trajectory state at a specific time
fn interpolate_trajectory_state(
    trajectory_points: &[(f64, [f64; 6])],
    target_time: f64,
) -> Result<[f64; 6], String> {
    if trajectory_points.is_empty() {
        return Err("No trajectory points available".to_string());
    }
    
    // Find bounding points
    let mut lower_idx = 0;
    let mut upper_idx = trajectory_points.len() - 1;
    
    for (i, &(time, _)) in trajectory_points.iter().enumerate() {
        if time <= target_time {
            lower_idx = i;
        }
        if time >= target_time && upper_idx == trajectory_points.len() - 1 {
            upper_idx = i;
            break;
        }
    }
    
    if lower_idx == upper_idx {
        return Ok(trajectory_points[lower_idx].1);
    }
    
    // Linear interpolation
    let (t1, state1) = trajectory_points[lower_idx];
    let (t2, state2) = trajectory_points[upper_idx];
    
    if (t2 - t1).abs() < f64::EPSILON {
        return Ok(state1);
    }
    
    let alpha = (target_time - t1) / (t2 - t1);
    let mut interpolated_state = [0.0; 6];
    
    for i in 0..6 {
        interpolated_state[i] = state1[i] + alpha * (state2[i] - state1[i]);
    }
    
    Ok(interpolated_state)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_brent_root_find() {
        // Test with simple quadratic: x^2 - 4 = 0, root at x = 2
        let f = |x: f64| x * x - 4.0;
        let root = brent_root_find(f, 1.0, 3.0, 1e-6, 100).unwrap();
        assert!((root - 2.0).abs() < 1e-6);
    }
    
    #[test]
    fn test_interpolate_trajectory_state() {
        let points = vec![
            (0.0, [0.0, 0.0, 0.0, 100.0, 50.0, 0.0]),
            (1.0, [100.0, 45.0, 0.0, 99.0, 40.0, 0.0]),
            (2.0, [200.0, 80.0, 0.0, 98.0, 30.0, 0.0]),
        ];
        
        let result = interpolate_trajectory_state(&points, 1.5).unwrap();
        
        // Should be halfway between points at t=1.0 and t=2.0
        assert!((result[0] - 150.0).abs() < 1e-10); // x position
        assert!((result[1] - 62.5).abs() < 1e-10);  // y position
        assert!((result[3] - 98.5).abs() < 1e-10);  // vx velocity
    }
    
    #[test]
    fn test_find_trajectory_apex() {
        let points = vec![
            (0.0, [0.0, 0.0, 0.0, 100.0, 50.0, 0.0]),
            (1.0, [100.0, 45.0, 0.0, 99.0, 40.0, 0.0]),
            (2.0, [200.0, 80.0, 0.0, 98.0, 30.0, 0.0]), // Peak here
            (3.0, [300.0, 75.0, 0.0, 97.0, 20.0, 0.0]),
            (4.0, [400.0, 60.0, 0.0, 96.0, 10.0, 0.0]),
        ];
        
        let (max_ord, max_ord_x) = find_trajectory_apex(&points, 500.0);
        
        assert!((max_ord - 80.0).abs() < 1e-10);
        assert!((max_ord_x - 200.0).abs() < 1e-10);
    }
}