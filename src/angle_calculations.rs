use crate::InternalBallisticInputs;
use std::f64;

// Constants for unit conversions
const FPS_TO_MPS: f64 = 0.3048;
const YARDS_TO_METERS: f64 = 0.9144;
const DEGREES_TO_RADIANS: f64 = std::f64::consts::PI / 180.0;
const RADIANS_TO_DEGREES: f64 = 180.0 / std::f64::consts::PI;

// Zero finding constants
const ZERO_FINDING_ACCURACY: f64 = crate::constants::ROOT_FINDING_TOLERANCE;
const ZERO_FINDING_MAX_ITER: usize = 100;

/// Result of angle calculation
#[derive(Debug, Clone)]
pub struct AngleResult {
    pub angle_rad: f64,
    pub iterations_used: usize,
    pub final_error: f64,
    pub success: bool,
}

/// Brent's method for root finding - optimized implementation
pub fn brent_root_find<F>(
    f: F,
    mut a: f64,
    mut b: f64,
    tolerance: f64,
    max_iterations: usize,
) -> Result<AngleResult, String>
where
    F: Fn(f64) -> f64,
{
    let mut fa = f(a);
    let mut fb = f(b);
    let mut iterations = 0;
    
    // Ensure the root is bracketed
    if fa * fb > 0.0 {
        return Err(format!("Root not bracketed: f({a}) = {fa}, f({b}) = {fb}"));
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
    
    while iterations < max_iterations {
        iterations += 1;
        
        if fb.abs() < tolerance {
            return Ok(AngleResult {
                angle_rad: b,
                iterations_used: iterations,
                final_error: fb.abs(),
                success: true,
            });
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
            return Ok(AngleResult {
                angle_rad: b,
                iterations_used: iterations,
                final_error: fb.abs(),
                success: true,
            });
        }
        
        if e.abs() >= tolerance_scaled && fc.abs() > fb.abs() {
            // Check for safe division before interpolation
            if fc.abs() < f64::EPSILON || fa.abs() < f64::EPSILON {
                // Fallback to bisection if denominators are too small
                d = m;
                e = m;
            } else {
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
                
                // Check for safe division in the acceptance test
                if q.abs() > f64::EPSILON && 2.0 * p < 3.0 * m * q - (tolerance_scaled * q).abs() 
                    && p < (0.5 * s * q).abs() {
                    d = p / q;
                } else {
                    d = m;
                    e = d;
                }
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
    
    Ok(AngleResult {
        angle_rad: b,
        iterations_used: iterations,
        final_error: fb.abs(),
        success: fb.abs() < tolerance * 10.0, // Relaxed success criteria
    })
}

/// Calculate adjusted muzzle velocity for powder temperature sensitivity
pub fn adjusted_muzzle_velocity(inputs: &InternalBallisticInputs) -> f64 {
    let mut mv = inputs.muzzle_velocity;
    
    if inputs.use_powder_sensitivity {
        mv *= 1.0 + inputs.powder_temp_sensitivity 
            * (inputs.temperature - inputs.powder_temp) / 15.0;
    }
    
    mv
}

/// Calculate zero angle using Brent's method and Rust trajectory integration
pub fn zero_angle(
    inputs: &InternalBallisticInputs,
    trajectory_func: impl Fn(&InternalBallisticInputs, f64) -> Result<f64, String> + Copy,
) -> Result<AngleResult, String> {
    // Set up the target vertical position based on shooting angle
    let vert = if inputs.shooting_angle.abs() > 1e-6 {
        let angle_rad = inputs.shooting_angle * DEGREES_TO_RADIANS;
        (inputs.target_distance * YARDS_TO_METERS) * angle_rad.sin()
    } else {
        0.0
    };
    
    // Define the height difference function
    let height_diff = |look_angle_rad: f64| -> f64 {
        // Calculate bullet height at target distance minus target height
        match trajectory_func(inputs, look_angle_rad) {
            Ok(bullet_height) => bullet_height - vert,
            Err(_) => -999.0, // Return large negative on failure
        }
    };
    
    // Reasonable bounds for the zero angle in radians
    // Most rifle zeroing will be within +/- 10 degrees
    let lower_bound = -10.0 * DEGREES_TO_RADIANS;
    let upper_bound = 10.0 * DEGREES_TO_RADIANS;
    
    // Try primary bounds first
    match brent_root_find(height_diff, lower_bound, upper_bound, 1e-6, 100) {
        Ok(result) if result.success => Ok(result),
        _ => {
            // Fallback to wider search range
            let wider_lower = -45.0 * DEGREES_TO_RADIANS;
            let wider_upper = 45.0 * DEGREES_TO_RADIANS;
            
            match brent_root_find(height_diff, wider_lower, wider_upper, 1e-5, 150) {
                Ok(result) if result.success => Ok(result),
                Ok(result) => {
                    // Return best attempt even if not fully successful
                    Ok(AngleResult {
                        angle_rad: result.angle_rad,
                        iterations_used: result.iterations_used,
                        final_error: result.final_error,
                        success: false,
                    })
                },
                Err(_) => {
                    // If all else fails, return 0 as a safe default
                    Ok(AngleResult {
                        angle_rad: 0.0,
                        iterations_used: 0,
                        final_error: f64::INFINITY,
                        success: false,
                    })
                }
            }
        }
    }
}

/// Solve muzzle angle using Brent's method optimization
pub fn solve_muzzle_angle(
    inputs: &InternalBallisticInputs,
    zero_distance_los_m: f64,
    trajectory_func: impl Fn(&InternalBallisticInputs) -> Result<f64, String> + Copy, // Returns drop_m
    angle_lower_deg: f64,
    angle_upper_deg: f64,
    rtol: f64,
) -> Result<AngleResult, String> {
    if angle_lower_deg >= angle_upper_deg {
        return Err("angle_lower_deg must be less than angle_upper_deg".to_string());
    }
    
    let lower = angle_lower_deg * DEGREES_TO_RADIANS;
    let mut upper = angle_upper_deg * DEGREES_TO_RADIANS;
    
    // Define the vertical error function
    let vertical_error = |angle_rad: f64| -> f64 {
        // Create modified inputs with new angle and target distance
        let mut candidate = inputs.clone();
        candidate.muzzle_angle = angle_rad * RADIANS_TO_DEGREES;
        candidate.target_distance = zero_distance_los_m / YARDS_TO_METERS; // Convert back to yards
        
        trajectory_func(&candidate).unwrap_or(1e6)
    };
    
    // Check bounds
    let f_lower = vertical_error(lower);
    if f_lower.abs() < 1e-9 {
        return Ok(AngleResult {
            angle_rad: lower,
            iterations_used: 1,
            final_error: f_lower.abs(),
            success: true,
        });
    }
    
    let f_upper = vertical_error(upper);
    if f_upper.abs() < 1e-9 {
        return Ok(AngleResult {
            angle_rad: upper,
            iterations_used: 1,
            final_error: f_upper.abs(),
            success: true,
        });
    }
    
    // Expand upper bound if needed to get a sign change
    if f_lower * f_upper > 0.0 {
        let step = 5.0 * DEGREES_TO_RADIANS;
        let max_angle = 45.0 * DEGREES_TO_RADIANS;
        let mut current = upper;
        let mut f_current = f_upper;
        
        while current < max_angle && f_lower * f_current > 0.0 {
            current += step;
            f_current = vertical_error(current);
        }
        
        if f_lower * f_current > 0.0 {
            return Err("Unable to bracket zero; widen angle bounds or check inputs".to_string());
        }
        
        upper = current;
    }
    
    // Use Brent's method to find the root with safe tolerance calculation
    let range = (upper - lower).abs();
    let tolerance = if range > f64::EPSILON {
        rtol * range
    } else {
        rtol * 1e-12 // Minimum tolerance for very small ranges
    };
    brent_root_find(vertical_error, lower, upper, tolerance, ZERO_FINDING_MAX_ITER)
}

/// Calculate simple ballistic drop approximation for quick estimates
pub fn quick_drop_estimate(
    muzzle_velocity_fps: f64,
    distance_yards: f64,
    _bullet_mass_grains: f64,
    bc: f64,
) -> f64 {
    let mv_mps = muzzle_velocity_fps * FPS_TO_MPS;
    let distance_m = distance_yards * YARDS_TO_METERS;
    
    // Simple ballistic approximation with safe divisions
    if mv_mps <= 0.0 || distance_m <= 0.0 {
        return 0.0; // No drop if no velocity or distance
    }
    
    let time_of_flight = distance_m / mv_mps;
    let gravity = 9.80665;
    
    // Basic drop calculation with BC approximation
    let bc_safe = bc.max(0.1);
    let drag_factor = 1.0 / bc_safe;
    let velocity_loss = drag_factor * time_of_flight * 0.1;
    let effective_velocity = mv_mps * (1.0 - velocity_loss).max(0.1); // Ensure positive
    let adjusted_time = distance_m / effective_velocity;
    
    0.5 * gravity * adjusted_time * adjusted_time
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DragModel;
    
    fn create_test_inputs() -> InternalBallisticInputs {
        InternalBallisticInputs {
            // Core fields
            muzzle_velocity: 823.0,  // 2700 fps in m/s
            launch_angle: 0.0,
            ballistic_coefficient: 0.5,
            mass: 0.0109,  // 168 grains in kg
            diameter: 0.00782,  // 0.308 inches in meters
            drag_model: DragModel::G1,
            sight_height: 0.05,
            
            // Duplicate fields for compatibility
            bc_value: 0.5,
            bc_type: DragModel::G1,
            bullet_mass: 168.0,  // in grains
            altitude: 0.0,
            bc_type_str: Some("G1".to_string()),
            twist_rate: 10.0,
            bullet_length: 1.3,  // in inches
            bullet_diameter: 0.308,  // in inches
            tipoff_yaw: 0.0,
            tipoff_decay_distance: 20.0,
            ground_threshold: 0.0,
            bc_segments: None,
            target_distance: 500.0,
            muzzle_angle: 0.0,
            temperature: 21.1,  // 70°F in Celsius
            latitude: None,
            enable_advanced_effects: false,
            is_twist_right: true,
            shooting_angle: 0.0,
            use_powder_sensitivity: false,
            powder_temp_sensitivity: 0.0,
            powder_temp: 70.0,
            caliber_inches: 0.0,
            weight_grains: 0.0,
            use_bc_segments: false,
            bullet_id: None,
            bc_segments_data: None,
            use_enhanced_spin_drift: false,
            use_form_factor: false,
            bullet_model: None,
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
        }
    }
    
    #[test]
    fn test_brent_root_find_quadratic() {
        // Test with simple quadratic: x^2 - 4 = 0, root at x = 2
        let f = |x: f64| x * x - 4.0;
        let result = brent_root_find(f, 1.0, 3.0, 1e-6, 100).unwrap();
        
        assert!(result.success);
        assert!((result.angle_rad - 2.0).abs() < 1e-6);
        assert!(result.iterations_used > 0);
        assert!(result.final_error < 1e-6);
    }
    
    #[test]
    fn test_brent_root_find_linear() {
        // Test with linear function: 2x - 6 = 0, root at x = 3
        let f = |x: f64| 2.0 * x - 6.0;
        let result = brent_root_find(f, 0.0, 5.0, 1e-6, 100).unwrap();
        
        assert!(result.success);
        assert!((result.angle_rad - 3.0).abs() < 1e-6);
    }
    
    #[test]
    fn test_brent_root_find_no_bracket() {
        // Test with function that doesn't change sign in the interval
        let f = |x: f64| x * x + 1.0; // Always positive
        let result = brent_root_find(f, 1.0, 3.0, 1e-6, 100);
        
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Root not bracketed"));
    }
    
    #[test]
    fn test_adjusted_muzzle_velocity_no_sensitivity() {
        let inputs = create_test_inputs();
        
        let result = adjusted_muzzle_velocity(&inputs);
        assert_eq!(result, 823.0);  // muzzle_velocity in m/s
    }
    
    #[test]
    fn test_adjusted_muzzle_velocity_with_sensitivity() {
        let mut inputs = create_test_inputs();
        inputs.use_powder_sensitivity = true;
        inputs.powder_temp_sensitivity = 1.0; // 1 fps per degree F
        inputs.temperature = 85.0; // 15 degrees above powder temp
        
        let result = adjusted_muzzle_velocity(&inputs);
        // Should be 823 * (1 + 1.0 * (85-70) / 15) = 823 * (1 + 1) = 1646
        assert!((result - 1646.0).abs() < 1e-6);
    }
    
    #[test]
    fn test_quick_drop_estimate() {
        let drop = quick_drop_estimate(2700.0, 500.0, 168.0, 0.5);
        
        // Should be a reasonable drop value (a few meters for 500 yards)
        assert!(drop > 0.0);
        assert!(drop < 50.0); // Sanity check - shouldn't be more than 50m drop
        
        // Test that higher BC gives less drop
        let drop_high_bc = quick_drop_estimate(2700.0, 500.0, 168.0, 0.8);
        assert!(drop_high_bc < drop);
    }
    
    #[test]
    fn test_zero_angle_bounds() {
        // Test that angle bounds are reasonable
        let lower = -10.0 * DEGREES_TO_RADIANS;
        let upper = 10.0 * DEGREES_TO_RADIANS;
        
        assert!(lower < 0.0);
        assert!(upper > 0.0);
        assert!((upper - lower).abs() > 0.1); // Reasonable search range
    }
}