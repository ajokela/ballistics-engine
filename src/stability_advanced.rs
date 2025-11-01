// Advanced stability calculations using refined Miller formula and modern corrections
// Based on:
// - Don Miller's refined stability formula (2005)
// - Bryan Litz's stability refinements
// - Geoffrey Kolbe's Bowman-Howell stability improvements
//
// NOTE: Some advanced stability functions are experimental and kept for future use.
#![allow(dead_code)]

/// Advanced stability parameters for different bullet types
#[derive(Debug, Clone)]
pub struct StabilityParameters {
    /// Shape factor for nose profile (1.0 for tangent, 0.9 for secant)
    pub nose_shape_factor: f64,
    /// Boat tail effectiveness factor
    pub boat_tail_factor: f64,
    /// Plastic tip correction
    pub plastic_tip_factor: f64,
    /// Center of pressure adjustment
    pub cop_adjustment: f64,
}

impl StabilityParameters {
    pub fn for_bullet_type(bullet_type: &str, has_boat_tail: bool, has_plastic_tip: bool) -> Self {
        let mut params = match bullet_type.to_lowercase().as_str() {
            "match" | "bthp" => Self {
                nose_shape_factor: 0.95,
                boat_tail_factor: if has_boat_tail { 0.94 } else { 1.0 },
                plastic_tip_factor: 1.0,
                cop_adjustment: 0.98,
            },
            "vld" | "very_low_drag" => Self {
                nose_shape_factor: 0.88,
                boat_tail_factor: if has_boat_tail { 0.92 } else { 1.0 },
                plastic_tip_factor: 1.0,
                cop_adjustment: 0.96,
            },
            "hybrid" => Self {
                nose_shape_factor: 0.91,
                boat_tail_factor: if has_boat_tail { 0.93 } else { 1.0 },
                plastic_tip_factor: 1.0,
                cop_adjustment: 0.97,
            },
            "hunting" => Self {
                nose_shape_factor: 0.98,
                boat_tail_factor: if has_boat_tail { 0.95 } else { 1.0 },
                plastic_tip_factor: if has_plastic_tip { 0.92 } else { 1.0 },
                cop_adjustment: 0.99,
            },
            _ => Self::default(),
        };
        
        if has_plastic_tip && params.plastic_tip_factor == 1.0 {
            params.plastic_tip_factor = 0.95;
        }
        
        params
    }
    
    pub fn default() -> Self {
        Self {
            nose_shape_factor: 1.0,
            boat_tail_factor: 1.0,
            plastic_tip_factor: 1.0,
            cop_adjustment: 1.0,
        }
    }
}

/// Calculate advanced Miller stability with modern corrections
pub fn calculate_advanced_stability(
    mass_grains: f64,
    velocity_fps: f64,
    twist_rate_inches: f64,
    caliber_inches: f64,
    length_inches: f64,
    air_density_kg_m3: f64,
    temperature_k: f64,
    bullet_type: &str,
    has_boat_tail: bool,
    has_plastic_tip: bool,
) -> f64 {
    if twist_rate_inches == 0.0 || caliber_inches == 0.0 || length_inches == 0.0 {
        return 0.0;
    }
    
    let params = StabilityParameters::for_bullet_type(bullet_type, has_boat_tail, has_plastic_tip);
    
    // Calculate base Miller stability
    let sg_base = calculate_miller_refined(
        mass_grains,
        twist_rate_inches,
        caliber_inches,
        length_inches,
        params.nose_shape_factor,
    );
    
    // Apply velocity correction (Miller's refined formula)
    let sg_velocity_corrected = apply_velocity_correction(sg_base, velocity_fps);
    
    // Apply atmospheric corrections
    let sg_atmosphere_corrected = apply_atmospheric_correction(
        sg_velocity_corrected,
        air_density_kg_m3,
        temperature_k,
    );
    
    // Apply boat tail correction if applicable
    let sg_boat_tail = sg_atmosphere_corrected * params.boat_tail_factor;
    
    // Apply plastic tip correction if applicable
    let sg_plastic_tip = sg_boat_tail * params.plastic_tip_factor;
    
    // Apply center of pressure adjustment
    let sg_final = sg_plastic_tip * params.cop_adjustment;
    
    // Apply Bowman-Howell dynamic stability correction for very high velocities
    if velocity_fps > 3000.0 {
        apply_bowman_howell_correction(sg_final, velocity_fps, caliber_inches)
    } else {
        sg_final
    }
}

/// Miller's refined stability formula (2005 version)
fn calculate_miller_refined(
    mass_grains: f64,
    twist_rate_inches: f64,
    caliber_inches: f64,
    length_inches: f64,
    nose_shape_factor: f64,
) -> f64 {
    // Convert to calibers
    let twist_calibers = twist_rate_inches / caliber_inches;
    let length_calibers = length_inches / caliber_inches;
    
    // Miller's constant (refined from original 30)
    const MILLER_CONSTANT: f64 = 30.0;
    
    // Calculate moment of inertia factor
    // For modern bullets: (1 + L²) where L is length in calibers
    let inertia_factor = 1.0 + length_calibers.powi(2);
    
    // Base Miller formula with nose shape correction
    let numerator = MILLER_CONSTANT * mass_grains * nose_shape_factor;
    let denominator = twist_calibers.powi(2) * caliber_inches.powi(3) * 
                     length_calibers * inertia_factor;
    
    if denominator == 0.0 {
        return 0.0;
    }
    
    numerator / denominator
}

/// Velocity correction using Miller's refined approach
fn apply_velocity_correction(sg_base: f64, velocity_fps: f64) -> f64 {
    // Miller's refined velocity correction
    // Uses 2800 fps as reference with cube root relationship
    const VELOCITY_REFERENCE: f64 = 2800.0;
    
    // For velocities below 1400 fps, use modified correction
    if velocity_fps < 1400.0 {
        let velocity_factor = (velocity_fps / 1400.0).powf(0.5);
        sg_base * velocity_factor * 0.5 // Additional reduction for very low velocities
    } else {
        // Standard Miller velocity correction
        let velocity_factor = (velocity_fps / VELOCITY_REFERENCE).powf(1.0 / 3.0);
        sg_base * velocity_factor
    }
}

/// Atmospheric correction for non-standard conditions
fn apply_atmospheric_correction(
    sg: f64,
    air_density_kg_m3: f64,
    temperature_k: f64,
) -> f64 {
    // Standard atmosphere at sea level
    const STD_DENSITY: f64 = 1.225; // kg/m³
    const STD_TEMP: f64 = 288.15;   // K (15°C)
    
    // Density altitude correction
    let density_ratio = STD_DENSITY / air_density_kg_m3;
    let density_correction = density_ratio.sqrt();
    
    // Temperature correction (affects speed of sound and viscosity)
    let temp_ratio = temperature_k / STD_TEMP;
    let temp_correction = temp_ratio.powf(0.17); // Empirical exponent
    
    sg * density_correction * temp_correction
}

/// Bowman-Howell correction for hypervelocity projectiles
fn apply_bowman_howell_correction(sg: f64, velocity_fps: f64, caliber_inches: f64) -> f64 {
    // For velocities above 3000 fps, additional dynamic effects occur
    if velocity_fps <= 3000.0 {
        return sg;
    }
    
    // Calculate Mach number
    let mach = velocity_fps / 1125.0; // Approximate speed of sound in fps
    
    // Hypervelocity correction factor
    let excess_velocity = (velocity_fps - 3000.0) / 1000.0;
    let mach_correction = 1.0 - 0.05 * excess_velocity.min(2.0);
    
    // Small caliber bullets are more affected
    let caliber_factor = if caliber_inches < 0.264 {
        0.95
    } else if caliber_inches < 0.308 {
        0.97
    } else {
        1.0
    };
    
    sg * mach_correction * caliber_factor
}

/// Calculate dynamic stability including yaw effects
pub fn calculate_dynamic_stability(
    static_stability: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    yaw_angle_rad: f64,
    caliber_m: f64,
    mass_kg: f64,
) -> f64 {
    // Dynamic stability accounts for yaw and precession
    
    // Calculate spin parameter
    let spin_param = if velocity_mps > 0.0 {
        spin_rate_rad_s * caliber_m / (2.0 * velocity_mps)
    } else {
        0.0
    };
    
    // Yaw effect on stability
    let yaw_factor = 1.0 - 0.1 * yaw_angle_rad.abs().min(0.1);
    
    // Precession damping factor
    let precession_factor = 1.0 + 0.05 * spin_param.min(0.5);
    
    static_stability * yaw_factor * precession_factor
}

/// Predict stability over trajectory with velocity decay
pub fn predict_stability_at_distance(
    initial_stability: f64,
    initial_velocity_fps: f64,
    current_velocity_fps: f64,
    spin_decay_factor: f64, // Typically 0.95-0.98
) -> f64 {
    if initial_velocity_fps == 0.0 || current_velocity_fps == 0.0 {
        return initial_stability;
    }
    
    // Velocity ratio
    let velocity_ratio = current_velocity_fps / initial_velocity_fps;
    
    // Spin decays slower than velocity
    let spin_ratio = velocity_ratio * spin_decay_factor;
    
    // Stability changes with velocity and spin
    // SG ∝ (spin²/velocity)
    let stability_ratio = spin_ratio.powi(2) / velocity_ratio;
    
    initial_stability * stability_ratio
}

/// Check if bullet will remain stable throughout trajectory
pub fn check_trajectory_stability(
    muzzle_stability: f64,
    muzzle_velocity_fps: f64,
    terminal_velocity_fps: f64,
    spin_decay_factor: f64,
) -> (bool, f64, String) {
    let terminal_stability = predict_stability_at_distance(
        muzzle_stability,
        muzzle_velocity_fps,
        terminal_velocity_fps,
        spin_decay_factor,
    );
    
    let is_stable = terminal_stability >= 1.3; // Minimum for adequate stability
    
    let status = if terminal_stability < 1.0 {
        "UNSTABLE - Bullet will tumble".to_string()
    } else if terminal_stability < 1.3 {
        "MARGINAL - May experience accuracy issues".to_string()
    } else if terminal_stability < 1.5 {
        "ADEQUATE - Acceptable for most conditions".to_string()
    } else if terminal_stability < 2.5 {
        "GOOD - Optimal stability".to_string()
    } else {
        "OVER-STABILIZED - May reduce BC slightly".to_string()
    };
    
    (is_stable, terminal_stability, status)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_advanced_stability() {
        // Test with .308 168gr Match bullet
        let stability = calculate_advanced_stability(
            168.0,      // mass in grains
            2700.0,     // velocity in fps
            10.0,       // twist rate in inches
            0.308,      // caliber in inches
            1.24,       // length in inches
            1.225,      // air density
            288.15,     // temperature in K
            "match",    // bullet type
            true,       // has boat tail
            false,      // no plastic tip
        );

        println!("Calculated stability: {}", stability);

        // Should give stability around 1.4-1.8 for typical .308 Match
        assert!(stability > 1.3);
        assert!(stability < 2.5, "Stability {} exceeds upper bound", stability);
    }

    #[test]
    fn test_stability_prediction() {
        // Use higher initial stability to maintain adequate stability through velocity drop
        let (is_stable, terminal_sg, status) = check_trajectory_stability(
            2.2,        // muzzle stability (well above marginal)
            2700.0,     // muzzle velocity
            1900.0,     // terminal velocity (moderate drop)
            0.98,       // spin decay factor (good spin retention)
        );

        println!("is_stable: {}, terminal_sg: {}, status: {}", is_stable, terminal_sg, status);

        assert!(is_stable, "Expected stable trajectory but got: is_stable={}, terminal_sg={}, status={}", is_stable, terminal_sg, status);
        assert!(terminal_sg > 1.0, "Terminal SG {} too low", terminal_sg);
        assert!(status.contains("ADEQUATE") || status.contains("GOOD") || status.contains("MARGINAL"));
    }
}