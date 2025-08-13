//! Spin Decay Physics for Ballistics Calculations
//!
//! This module implements realistic spin decay modeling based on:
//! - Aerodynamic torque opposing spin
//! - Skin friction effects
//! - Velocity-dependent decay rates
//! - Surface roughness effects

use std::f64::consts::PI;

/// Parameters affecting spin decay rate
#[derive(Debug, Clone, Copy)]
pub struct SpinDecayParameters {
    /// Surface roughness in meters (typical: 0.1mm)
    pub surface_roughness: f64,
    /// Very low for streamlined spinning bullets
    pub skin_friction_coefficient: f64,
    /// Shape factor for spin damping
    pub form_factor: f64,
}

impl SpinDecayParameters {
    /// Create default parameters
    pub fn new() -> Self {
        Self {
            surface_roughness: 0.0001,
            skin_friction_coefficient: 0.00001,
            form_factor: 1.0,
        }
    }

    /// Get typical parameters for different bullet types
    pub fn from_bullet_type(bullet_type: &str) -> Self {
        match bullet_type.to_lowercase().as_str() {
            "match" => Self {
                surface_roughness: 0.00005,
                skin_friction_coefficient: 0.000008,
                form_factor: 0.9,
            },
            "hunting" => Self {
                surface_roughness: 0.0001,
                skin_friction_coefficient: 0.00001,
                form_factor: 1.0,
            },
            "fmj" => Self {
                surface_roughness: 0.00015,
                skin_friction_coefficient: 0.000012,
                form_factor: 1.1,
            },
            "cast" => Self {
                surface_roughness: 0.0002,
                skin_friction_coefficient: 0.000015,
                form_factor: 1.2,
            },
            _ => Self::new(),
        }
    }
}

impl Default for SpinDecayParameters {
    fn default() -> Self {
        Self::new()
    }
}

/// Calculate the aerodynamic moment opposing spin
///
/// Based on:
/// - Skin friction on the bullet surface
/// - Pressure drag effects
/// - Magnus moment damping
pub fn calculate_spin_damping_moment(
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    caliber_m: f64,
    length_m: f64,
    decay_params: &SpinDecayParameters,
) -> f64 {
    if spin_rate_rad_s == 0.0 || velocity_mps == 0.0 {
        return 0.0;
    }

    // Reynolds number based on spin
    let radius = caliber_m / 2.0;
    let tangential_velocity = spin_rate_rad_s * radius;
    let _re_spin = air_density_kg_m3 * tangential_velocity * caliber_m / 1.81e-5; // Air viscosity

    // Skin friction coefficient (modified for spinning cylinder)
    let cf = decay_params.skin_friction_coefficient;

    // Surface area
    let surface_area = PI * caliber_m * length_m;

    // Tangential force due to skin friction
    let f_tangential = 0.5 * air_density_kg_m3 * cf * surface_area * tangential_velocity.powi(2);

    // Moment arm is the radius
    let moment_skin = f_tangential * radius * decay_params.form_factor;

    // Additional damping from Magnus effect
    let spin_ratio = tangential_velocity / velocity_mps.max(1.0);
    let magnus_damping_factor = 0.01 * spin_ratio; // Reduced empirical factor
    let moment_magnus = magnus_damping_factor * moment_skin;

    // Total damping moment (always opposes spin)
    moment_skin + moment_magnus
}

/// Calculate moment of inertia about the longitudinal axis
pub fn calculate_moment_of_inertia(
    mass_kg: f64,
    caliber_m: f64,
    _length_m: f64,
    shape: &str,
) -> f64 {
    let radius = caliber_m / 2.0;

    match shape {
        "cylinder" => {
            // Simple cylinder: I = (1/2) * m * r²
            0.5 * mass_kg * radius.powi(2)
        }
        "ogive" => {
            // Ogive shape has less mass at the edges
            0.4 * mass_kg * radius.powi(2)
        }
        "boat_tail" => {
            // Boat tail has even less mass at the rear
            0.35 * mass_kg * radius.powi(2)
        }
        _ => {
            // Default to cylinder
            0.5 * mass_kg * radius.powi(2)
        }
    }
}

/// Calculate the rate of spin decay in rad/s²
pub fn calculate_spin_decay_rate(
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    mass_grains: f64,
    caliber_inches: f64,
    length_inches: f64,
    decay_params: &SpinDecayParameters,
    bullet_shape: &str,
) -> f64 {
    // Convert units
    let mass_kg = mass_grains * 0.00006479891; // grains to kg
    let caliber_m = caliber_inches * 0.0254;
    let length_m = length_inches * 0.0254;

    // Calculate damping moment
    let damping_moment = calculate_spin_damping_moment(
        spin_rate_rad_s,
        velocity_mps,
        air_density_kg_m3,
        caliber_m,
        length_m,
        decay_params,
    );

    // Calculate moment of inertia
    let moment_of_inertia = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, bullet_shape);

    // Angular deceleration = -M / I
    if moment_of_inertia > 0.0 {
        -damping_moment / moment_of_inertia
    } else {
        0.0
    }
}

/// Calculate the spin rate after accounting for decay
///
/// Uses an empirical model based on published ballistics data.
/// Real bullets typically lose 5-15% of spin over a 3-second flight.
pub fn update_spin_rate(
    initial_spin_rad_s: f64,
    time_elapsed_s: f64,
    velocity_mps: f64,
    _air_density_kg_m3: f64,
    mass_grains: f64,
    _caliber_inches: f64,
    _length_inches: f64,
    decay_params: Option<&SpinDecayParameters>,
) -> f64 {
    if time_elapsed_s <= 0.0 {
        return initial_spin_rad_s;
    }

    // Mass factor (heavier bullets retain spin better)
    let mass_factor = (175.0 / mass_grains).sqrt(); // Normalized to 175gr

    // Velocity factor (higher velocity means more decay)
    let velocity_factor = velocity_mps / 850.0; // Normalized to 850 m/s

    // Base decay rate per second (empirical)
    let base_decay_rate = if let Some(params) = decay_params {
        if params.form_factor < 1.0 {
            // Match bullet
            0.025 // 2.5% per second
        } else {
            // Hunting/FMJ bullet
            0.04 // 4% per second
        }
    } else {
        0.04 // Default to hunting bullet
    };

    // Adjusted decay rate
    let decay_rate_per_second = base_decay_rate * mass_factor * velocity_factor;

    // Apply exponential decay
    let decay_factor = (-decay_rate_per_second * time_elapsed_s).exp();

    // Ensure reasonable bounds (minimum 50% retention over any flight)
    initial_spin_rad_s * decay_factor.clamp(0.5, 1.0)
}

/// Calculate a simple correction factor for spin-dependent effects
///
/// This returns a value between 0 and 1 that represents the fraction
/// of initial spin remaining.
pub fn calculate_spin_decay_correction_factor(
    time_elapsed_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    mass_grains: f64,
    caliber_inches: f64,
    length_inches: f64,
    decay_params: Option<&SpinDecayParameters>,
) -> f64 {
    if time_elapsed_s <= 0.0 {
        return 1.0;
    }

    // Initial spin doesn't matter for the ratio calculation
    let initial_spin = 1000.0; // rad/s (arbitrary reference)

    let current_spin = update_spin_rate(
        initial_spin,
        time_elapsed_s,
        velocity_mps,
        air_density_kg_m3,
        mass_grains,
        caliber_inches,
        length_inches,
        decay_params,
    );

    current_spin / initial_spin
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_decay_parameters() {
        let match_params = SpinDecayParameters::from_bullet_type("match");
        assert_eq!(match_params.form_factor, 0.9);
        assert_eq!(match_params.surface_roughness, 0.00005);

        let hunting_params = SpinDecayParameters::from_bullet_type("hunting");
        assert_eq!(hunting_params.form_factor, 1.0);
    }

    #[test]
    fn test_moment_of_inertia() {
        let mass_kg = 0.01134; // 175 grains
        let caliber_m = 0.00782; // .308 inches

        let i_cylinder = calculate_moment_of_inertia(mass_kg, caliber_m, 0.033, "cylinder");
        let i_ogive = calculate_moment_of_inertia(mass_kg, caliber_m, 0.033, "ogive");

        assert!(i_cylinder > i_ogive); // Cylinder has more inertia than ogive
    }

    #[test]
    fn test_spin_decay_realistic() {
        // Test realistic spin decay for a .308 match bullet
        let initial_spin = 2800.0 * 2.0 * PI; // 2800 rev/s to rad/s
        let params = SpinDecayParameters::from_bullet_type("match");

        // After 3 seconds of flight
        let spin_after_3s = update_spin_rate(
            initial_spin,
            3.0,        // 3 seconds
            750.0,      // 750 m/s average velocity
            1.2,        // air density
            175.0,      // 175 grains
            0.308,      // caliber
            1.3,        // length
            Some(&params),
        );

        let decay_percent = (1.0 - spin_after_3s / initial_spin) * 100.0;

        // Should lose between 2% and 15% of spin
        assert!(decay_percent > 2.0 && decay_percent < 15.0);
    }

    #[test]
    fn test_spin_decay_bounds() {
        let initial_spin = 1000.0;
        let params = SpinDecayParameters::new();

        // Test extreme time - should never lose more than 50%
        let spin_long_time = update_spin_rate(
            initial_spin,
            100.0, // Very long flight time
            500.0,
            1.225,
            150.0,
            0.308,
            1.2,
            Some(&params),
        );

        assert!(spin_long_time >= initial_spin * 0.5);
    }
}