//! Pitch Damping Moment Physics for Ballistics Calculations
//!
//! This module implements pitch damping moments that affect:
//! - Dynamic stability during flight
//! - Precession and nutation damping  
//! - Yaw of repose convergence
//! - Transonic stability transitions

use std::f64::consts::PI;

/// Aerodynamic damping coefficients for different flight regimes
#[derive(Debug, Clone, Copy)]
pub struct PitchDampingCoefficients {
    pub subsonic: f64,       // Cmq + Cmα̇ for M < 0.8
    pub transonic_low: f64,  // For 0.8 <= M < 1.0
    pub transonic_high: f64, // For 1.0 <= M < 1.2 (can be destabilizing)
    pub supersonic: f64,     // For M >= 1.2
}

impl Default for PitchDampingCoefficients {
    fn default() -> Self {
        Self {
            subsonic: -0.8,
            transonic_low: -0.3,
            transonic_high: 0.2,
            supersonic: -0.5,
        }
    }
}

impl PitchDampingCoefficients {
    /// Get typical coefficients for different bullet types
    pub fn from_bullet_type(bullet_type: &str) -> Self {
        match bullet_type.to_lowercase().as_str() {
            "match_boat_tail" => Self {
                subsonic: -0.9,
                transonic_low: -0.4,
                transonic_high: 0.1,
                supersonic: -0.6,
            },
            "match_flat_base" => Self {
                subsonic: -0.7,
                transonic_low: -0.2,
                transonic_high: 0.3,
                supersonic: -0.4,
            },
            "vld" => Self {
                // Very Low Drag - more stable
                subsonic: -1.0,
                transonic_low: -0.5,
                transonic_high: -0.1,
                supersonic: -0.7,
            },
            "hunting" => Self {
                subsonic: -0.6,
                transonic_low: -0.1,
                transonic_high: 0.4,
                supersonic: -0.3,
            },
            "fmj" => Self {
                subsonic: -0.7,
                transonic_low: -0.2,
                transonic_high: 0.2,
                supersonic: -0.5,
            },
            _ => Self::default(),
        }
    }
}

/// Calculate pitch damping coefficient based on Mach number
pub fn calculate_pitch_damping_coefficient(mach: f64, coeffs: &PitchDampingCoefficients) -> f64 {
    if mach < 0.8 {
        // Subsonic - stable damping
        coeffs.subsonic
    } else if mach < 1.0 {
        // Lower transonic - decreasing stability
        // Linear interpolation
        let t = (mach - 0.8) / 0.2;
        coeffs.subsonic * (1.0 - t) + coeffs.transonic_low * t
    } else if mach < 1.2 {
        // Upper transonic - potentially destabilizing
        // This is where "transonic jump" occurs
        let t = (mach - 1.0) / 0.2;
        coeffs.transonic_low * (1.0 - t) + coeffs.transonic_high * t
    } else {
        // Supersonic - returns to stable
        // Asymptotic approach to supersonic value
        let t = ((mach - 1.2) / 0.8).min(1.0);
        coeffs.transonic_high * (1.0 - t) + coeffs.supersonic * t
    }
}

/// Calculate the aerodynamic moment opposing pitch motion
pub fn calculate_pitch_damping_moment(
    pitch_rate_rad_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    caliber_m: f64,
    length_m: f64,
    mach: f64,
    coeffs: &PitchDampingCoefficients,
) -> f64 {
    if velocity_mps == 0.0 || pitch_rate_rad_s == 0.0 {
        return 0.0;
    }

    // Get damping coefficient for current Mach
    let cmq = calculate_pitch_damping_coefficient(mach, coeffs);

    // Dynamic pressure
    let q = 0.5 * air_density_kg_m3 * velocity_mps.powi(2);

    // Reference area (cross-sectional)
    let s = PI * (caliber_m / 2.0).powi(2);

    // Reference length (use diameter for missiles/projectiles)
    let d = caliber_m;

    // Non-dimensional pitch rate
    let q_nondim = pitch_rate_rad_s * d / velocity_mps;

    // Pitch damping moment
    // Negative because it opposes motion
    q * s * d * cmq * q_nondim
}

/// Calculate moment of inertia about transverse axis (pitch/yaw)
pub fn calculate_transverse_moment_of_inertia(
    mass_kg: f64,
    caliber_m: f64,
    length_m: f64,
    shape: &str,
) -> f64 {
    let radius = caliber_m / 2.0;

    match shape {
        "cylinder" => {
            // I_transverse = m * (3*r² + L²) / 12
            mass_kg * (3.0 * radius.powi(2) + length_m.powi(2)) / 12.0
        }
        "ogive" => {
            // Ogive has more mass toward the front
            // Approximate as 85% of cylinder value
            let cylinder_i = mass_kg * (3.0 * radius.powi(2) + length_m.powi(2)) / 12.0;
            0.85 * cylinder_i
        }
        "boat_tail" => {
            // Boat tail has less mass at rear
            // Approximate as 80% of cylinder value
            let cylinder_i = mass_kg * (3.0 * radius.powi(2) + length_m.powi(2)) / 12.0;
            0.80 * cylinder_i
        }
        _ => {
            // Default to cylinder
            mass_kg * (3.0 * radius.powi(2) + length_m.powi(2)) / 12.0
        }
    }
}

/// Calculate angular acceleration from moment and inertia
pub fn calculate_angular_acceleration(moment: f64, moment_of_inertia: f64) -> f64 {
    if moment_of_inertia > 0.0 {
        moment / moment_of_inertia
    } else {
        0.0
    }
}

/// Calculate yaw of repose with pitch damping effects
pub fn calculate_damped_yaw_of_repose(
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
) -> (f64, f64) {
    if stability_factor <= 1.0 || spin_rate_rad_s == 0.0 {
        return (0.0, 0.0);
    }

    // Convert units
    let caliber_m = caliber_inches * 0.0254;
    let length_m = length_inches * 0.0254;
    let mass_kg = mass_grains * 0.00006479891;

    // Base yaw from crosswind or trajectory curvature
    let base_yaw_rad = if wind_velocity_mps != 0.0 && velocity_mps > 0.0 {
        (wind_velocity_mps / velocity_mps).atan()
    } else {
        // Natural yaw from curved trajectory
        0.002 // ~0.1 degrees typical
    };

    // Get damping coefficients
    let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);

    // Calculate pitch damping moment
    let damping_moment = calculate_pitch_damping_moment(
        pitch_rate_rad_s,
        velocity_mps,
        air_density_kg_m3,
        caliber_m,
        length_m,
        mach,
        &coeffs,
    );

    // Calculate transverse moment of inertia
    let i_transverse = calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");

    // Angular acceleration from damping
    let angular_accel = calculate_angular_acceleration(damping_moment, i_transverse);

    // Time constant for convergence
    let time_constant = if angular_accel != 0.0 && pitch_rate_rad_s != 0.0 {
        (pitch_rate_rad_s / angular_accel).abs()
    } else {
        10.0 // Default large value
    };

    // Convergence rate (how fast it approaches equilibrium)
    let convergence_rate = 1.0 / time_constant;

    // Modify equilibrium yaw based on stability and damping
    let stability_factor_clamped = stability_factor.min(10.0);
    let mut damping_factor = 1.0 / (1.0 + (stability_factor_clamped - 1.0).sqrt());

    // In transonic region, damping may be positive (destabilizing)
    if (0.8..=1.2).contains(&mach) {
        let cmq = calculate_pitch_damping_coefficient(mach, &coeffs);
        if cmq > 0.0 {
            // Destabilizing - increase effective yaw
            damping_factor *= 1.0 + cmq.abs();
        }
    }

    let equilibrium_yaw_rad = base_yaw_rad * damping_factor;

    (equilibrium_yaw_rad, convergence_rate)
}

/// Calculate precession rate including damping effects
pub fn calculate_precession_with_damping(
    yaw_angle_rad: f64,
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    pitch_damping_moment: f64,
    transverse_inertia: f64,
    spin_inertia: f64,
) -> f64 {
    if spin_rate_rad_s == 0.0 || velocity_mps == 0.0 {
        return 0.0;
    }

    // Basic precession rate (gyroscopic)
    let basic_precession = (spin_inertia * spin_rate_rad_s * yaw_angle_rad.sin()) 
        / (transverse_inertia * velocity_mps);

    // Damping modification
    let damping_factor = if transverse_inertia > 0.0 {
        let factor = 1.0 - (pitch_damping_moment.abs() / (transverse_inertia * velocity_mps));
        factor.max(0.1) // Minimum 10% precession
    } else {
        1.0
    };

    basic_precession * damping_factor
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pitch_damping_coefficient() {
        let coeffs = PitchDampingCoefficients::default();

        // Subsonic
        assert_eq!(calculate_pitch_damping_coefficient(0.5, &coeffs), -0.8);

        // Transonic
        let transonic = calculate_pitch_damping_coefficient(0.9, &coeffs);
        assert!(transonic > -0.8 && transonic < -0.3);

        // Supersonic
        let supersonic = calculate_pitch_damping_coefficient(2.0, &coeffs);
        assert_eq!(supersonic, -0.5);
    }

    #[test]
    fn test_pitch_damping_moment() {
        let coeffs = PitchDampingCoefficients::default();
        let moment = calculate_pitch_damping_moment(
            0.1,     // pitch rate
            300.0,   // velocity
            1.225,   // air density
            0.00782, // caliber (7.82mm)
            0.033,   // length (33mm)
            0.87,    // Mach
            &coeffs,
        );

        // Should be negative (opposing motion)
        assert!(moment < 0.0);
    }
}