use std::f64::consts::PI;
use crate::spin_decay::{update_spin_rate, SpinDecayParameters};
use crate::pitch_damping::{
    calculate_damped_yaw_of_repose, PitchDampingCoefficients,
    calculate_pitch_damping_moment,
};

/// Components of enhanced spin drift calculation
#[derive(Debug, Clone)]
pub struct SpinDriftComponents {
    pub spin_rate_rps: f64,        // Revolutions per second
    pub spin_rate_rad_s: f64,      // Radians per second
    pub stability_factor: f64,     // Gyroscopic stability (Sg)
    pub yaw_of_repose_rad: f64,    // Equilibrium yaw angle
    pub drift_rate_mps: f64,       // Lateral drift rate (m/s)
    pub total_drift_m: f64,        // Total drift at current time
    pub magnus_component_m: f64,   // Magnus effect contribution
    pub gyroscopic_component_m: f64, // Pure gyroscopic drift
    pub pitch_damping_moment: f64, // Pitch damping moment (N⋅m)
    pub yaw_convergence_rate: f64, // Convergence rate to equilibrium (rad/s)
    pub pitch_rate_rad_s: f64,     // Current pitch/yaw rate (rad/s)
}

/// Calculate bullet spin rate from velocity and twist rate
pub fn calculate_spin_rate(velocity_mps: f64, twist_rate_inches: f64) -> (f64, f64) {
    if twist_rate_inches <= 0.0 {
        return (0.0, 0.0);
    }
    
    // Convert velocity to inches/second
    let velocity_ips = velocity_mps * 39.3701;
    
    // Calculate revolutions per second
    let spin_rate_rps = velocity_ips / twist_rate_inches;
    
    // Convert to radians per second
    let spin_rate_rad_s = spin_rate_rps * 2.0 * PI;
    
    (spin_rate_rps, spin_rate_rad_s)
}

/// Calculate dynamic gyroscopic stability factor using Miller formula
pub fn calculate_dynamic_stability(
    bullet_mass_grains: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    caliber_inches: f64,
    length_inches: f64,
    air_density_kg_m3: f64,
) -> f64 {
    if spin_rate_rad_s == 0.0 || velocity_mps == 0.0 {
        return 0.0;
    }
    
    // Convert velocity to fps for Miller formula
    let velocity_fps = velocity_mps * 3.28084;
    
    // Calculate twist rate in calibers
    if caliber_inches > 0.0 {
        // Back-calculate twist rate from spin rate
        let spin_rps = spin_rate_rad_s / (2.0 * PI);
        let velocity_ips = velocity_fps * 12.0; // inches per second
        let twist_inches = if spin_rps > 0.0 { velocity_ips / spin_rps } else { 0.0 };
        let twist_calibers = if twist_inches > 0.0 { twist_inches / caliber_inches } else { 0.0 };
        
        // Length to diameter ratio
        let length_calibers = if caliber_inches > 0.0 { length_inches / caliber_inches } else { 0.0 };
        
        // Miller stability formula (simplified)
        // Sg = 30 * m / (t^2 * d^3 * l * (1 + l^2))
        // Where: m = mass in grains, t = twist in calibers, d = diameter in inches
        //        l = length in calibers
        
        if twist_calibers == 0.0 || length_calibers == 0.0 {
            return 0.0;
        }
        
        let numerator = 30.0 * bullet_mass_grains;
        let denominator = twist_calibers.powi(2) * caliber_inches.powi(3) * 
                         length_calibers * (1.0 + length_calibers.powi(2));
        
        if denominator == 0.0 {
            return 0.0;
        }
        
        // Base stability
        let sg_base = numerator / denominator;
        
        // Velocity correction (compared to standard 2800 fps)
        let velocity_factor = (velocity_fps / 2800.0).powf(1.0 / 3.0);
        
        // Atmospheric correction
        // Standard conditions: 59°F, 29.92 inHg = 1.225 kg/m³
        let density_factor = (1.225 / air_density_kg_m3).sqrt();
        
        // Final stability
        sg_base * velocity_factor * density_factor
    } else {
        0.0
    }
}

/// Calculate the yaw of repose (equilibrium yaw angle)
pub fn calculate_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    wind_velocity_mps: f64,
    pitch_rate_rad_s: f64,
    air_density_kg_m3: f64,
    caliber_inches: f64,
    length_inches: f64,
    mass_grains: f64,
    mach: f64,
    bullet_type: &str,
    use_pitch_damping: bool,
) -> (f64, f64) {
    if stability_factor <= 1.0 || spin_rate_rad_s == 0.0 {
        return (0.0, 0.0);
    }
    
    // Use enhanced calculation with pitch damping if requested
    if use_pitch_damping && mach > 0.0 {
        // Map bullet types for pitch damping
        let damping_type = match bullet_type.to_lowercase().as_str() {
            "match" => "match_boat_tail",
            "hunting" => "hunting",
            "fmj" => "fmj",
            "vld" => "vld",
            _ => "match_boat_tail",
        };
        
        return calculate_damped_yaw_of_repose(
            stability_factor, velocity_mps, spin_rate_rad_s,
            wind_velocity_mps, pitch_rate_rad_s, air_density_kg_m3,
            caliber_inches, length_inches, mass_grains, mach, damping_type
        );
    }
    
    // Original calculation (backward compatibility)
    // Crosswind component creates yaw
    let yaw_rad = if wind_velocity_mps == 0.0 {
        // No wind - use typical value for spin drift
        // Yaw develops due to nose following curved trajectory
        0.002 // ~0.1 degrees typical
    } else {
        // Wind-induced yaw
        if velocity_mps > 0.0 {
            (wind_velocity_mps / velocity_mps).atan()
        } else {
            0.0
        }
    };
    
    // Damping factor based on stability with safe division
    let stability_term = (stability_factor - 1.0).max(0.0).sqrt();
    let damping = 1.0 / (1.0 + stability_term);
    
    (yaw_rad * damping, 0.0)  // No convergence rate in simple model
}

/// Calculate Magnus effect contribution to drift
pub fn calculate_magnus_drift_component(
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    yaw_rad: f64,
    air_density_kg_m3: f64,
    caliber_inches: f64,
    time_s: f64,
    mass_grains: f64,
) -> f64 {
    let diameter_m = caliber_inches * 0.0254;
    let mass_kg = mass_grains * 0.00006479891;  // Convert grains to kg
    
    // Magnus force coefficient (empirical)
    // Varies with Mach number
    let mach = velocity_mps / 343.0; // Approximate speed of sound
    
    let cmag = if mach < 0.8 {
        0.25
    } else if mach < 1.2 {
        // Transonic reduction
        0.15
    } else {
        // Supersonic
        0.10 + 0.05 * ((mach - 1.2) / 2.0).min(1.0)
    };
    
    // Spin ratio
    let spin_ratio = (spin_rate_rad_s * diameter_m / 2.0) / velocity_mps;
    
    // Magnus force
    let magnus_force = if velocity_mps > 0.0 {
        cmag * spin_ratio * yaw_rad * 
        0.5 * air_density_kg_m3 * velocity_mps.powi(2) * 
        PI * (diameter_m / 2.0).powi(2)
    } else {
        0.0
    };
    
    // Convert force to acceleration by dividing by mass
    let magnus_accel = magnus_force / mass_kg;
    
    // Drift over time (simplified - should integrate)
    
    
    0.5 * magnus_accel * time_s.powi(2)
}

/// Calculate pure gyroscopic drift (Poisson effect)
pub fn calculate_gyroscopic_drift(
    stability_factor: f64,
    _yaw_rad: f64,
    velocity_mps: f64,
    time_s: f64,
    is_right_twist: bool,
) -> f64 {
    if stability_factor <= 1.0 || time_s <= 0.0 {
        return 0.0;
    }
    
    // Litz formula is not reliable for subsonic flight. Disable it.
    let velocity_fps = velocity_mps * 3.28084;
    if velocity_fps < 1125.0 {
        return 0.0;
    }
    
    // Direction based on twist
    let sign = if is_right_twist { 1.0 } else { -1.0 };
    
    // Bryan Litz's empirical formula for spin drift
    let base_coefficient = 1.25 * (stability_factor + 1.2);
    let time_factor = time_s.powf(1.83);
    let drift_in = sign * base_coefficient * time_factor;
    
    // Convert to meters
    
    
    drift_in * 0.0254
}

/// Calculate enhanced spin drift with all components
pub fn calculate_enhanced_spin_drift(
    bullet_mass: f64,
    velocity_mps: f64,
    twist_rate: f64,
    bullet_diameter: f64,
    bullet_length: f64,
    is_twist_right: bool,
    time_s: f64,
    air_density: f64,
    crosswind_mps: f64,
    pitch_rate_rad_s: f64,
    use_pitch_damping: bool,
) -> SpinDriftComponents {
    // Calculate initial spin rate (at muzzle)
    let muzzle_velocity = velocity_mps; // Assuming we're passed muzzle velocity
    let (_initial_spin_rps, initial_spin_rad_s) = calculate_spin_rate(muzzle_velocity, twist_rate);
    
    // Apply spin decay based on time of flight
    let decay_params = SpinDecayParameters::from_bullet_type("match"); // Default to match for now
    let current_spin_rad_s = update_spin_rate(
        initial_spin_rad_s,
        time_s,
        velocity_mps,
        air_density,
        bullet_mass * 15.432358, // Convert to grains
        bullet_diameter,
        bullet_length,
        Some(&decay_params),
    );
    
    let spin_rps = current_spin_rad_s / (2.0 * PI);
    let spin_rad_s = current_spin_rad_s;
    
    // Calculate dynamic stability
    let stability = calculate_dynamic_stability(
        bullet_mass,
        velocity_mps,
        spin_rad_s,
        bullet_diameter,
        bullet_length,
        air_density,
    );
    
    // Calculate Mach number for pitch damping
    let mach = velocity_mps / 343.0; // Approximate speed of sound
    
    // Determine bullet type (default to match for now)
    let bullet_type = "match";
    
    // Calculate yaw of repose with pitch damping
    let (yaw_rad, convergence_rate) = calculate_yaw_of_repose(
        stability,
        velocity_mps,
        spin_rad_s,
        crosswind_mps,
        pitch_rate_rad_s,
        air_density,
        bullet_diameter,
        bullet_length,
        bullet_mass,
        mach,
        bullet_type,
        use_pitch_damping,
    );
    
    // Calculate Magnus component
    let magnus_drift = calculate_magnus_drift_component(
        velocity_mps,
        spin_rad_s,
        yaw_rad,
        air_density,
        bullet_diameter,
        time_s,
        bullet_mass,
    );
    
    // Calculate gyroscopic component
    let gyro_drift = calculate_gyroscopic_drift(
        stability,
        yaw_rad,
        velocity_mps,
        time_s,
        is_twist_right,
    );
    
    // Total drift
    let total_drift = magnus_drift + gyro_drift;
    
    // Drift rate (derivative)
    let drift_rate = if time_s > 0.0 {
        total_drift / time_s
    } else {
        0.0
    };
    
    // Calculate pitch damping moment if using enhanced model
    let pitch_damping_moment = if use_pitch_damping && mach > 0.0 {
        let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);
        calculate_pitch_damping_moment(
            pitch_rate_rad_s,
            velocity_mps,
            air_density,
            bullet_diameter * 0.0254, // Convert to meters
            bullet_length * 0.0254,    // Convert to meters
            mach,
            &coeffs,
        )
    } else {
        0.0
    };
    
    SpinDriftComponents {
        spin_rate_rps: spin_rps,
        spin_rate_rad_s: spin_rad_s,
        stability_factor: stability,
        yaw_of_repose_rad: yaw_rad,
        drift_rate_mps: drift_rate,
        total_drift_m: total_drift,
        magnus_component_m: magnus_drift,
        gyroscopic_component_m: gyro_drift,
        pitch_damping_moment,
        yaw_convergence_rate: convergence_rate,
        pitch_rate_rad_s,
    }
}

/// Apply enhanced spin drift acceleration to derivatives
pub fn apply_enhanced_spin_drift(
    derivatives: &mut [f64; 6],
    spin_components: &SpinDriftComponents,
    time_s: f64,
    is_right_twist: bool,
) {
    if time_s > 0.1 {
        // Calculate acceleration from drift
        // Using second derivative of position
        let spin_accel_z = 2.0 * spin_components.drift_rate_mps / time_s;
        
        // Apply based on twist direction
        let sign = if is_right_twist { 1.0 } else { -1.0 };
        derivatives[5] += sign * spin_accel_z;
    }
}

/// Simplified interface for compatibility with existing code
pub fn compute_enhanced_spin_drift_simple(
    time_s: f64,
    stability: f64,
    velocity_mps: f64,
    twist_rate: f64,
    is_twist_right: bool,
    _caliber: f64,
) -> f64 {
    if twist_rate <= 0.0 {
        return 0.0;
    }
    
    // Calculate initial spin rate
    let (_, initial_spin_rad_s) = calculate_spin_rate(velocity_mps, twist_rate);
    
    // Apply simple spin decay (assume 175gr bullet)
    let decay_params = SpinDecayParameters::from_bullet_type("match");
    let spin_rad_s = update_spin_rate(
        initial_spin_rad_s,
        time_s,
        velocity_mps,
        1.225, // Standard air density
        175.0, // Standard bullet weight
        _caliber,
        1.3, // Standard bullet length
        Some(&decay_params),
    );
    
    // Estimate yaw of repose (use simple model for compatibility)
    let (yaw_rad, _) = calculate_yaw_of_repose(
        stability, velocity_mps, spin_rad_s, 0.0,
        0.0, 1.225, _caliber, 1.3, 175.0,
        velocity_mps / 343.0, "match", false
    );
    
    // Calculate gyroscopic drift (primary component)
    
    
    calculate_gyroscopic_drift(
        stability,
        yaw_rad,
        velocity_mps,
        time_s,
        is_twist_right,
    )
}