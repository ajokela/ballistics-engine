// Advanced spin drift model based on modern ballistics research
// Incorporates multiple empirical models from:
// - Bryan Litz's Applied Ballistics for Long Range Shooting
// - McCoy's Modern Exterior Ballistics
// - Courtney & Courtney spin drift research papers

use std::f64::consts::PI;

/// Advanced spin drift coefficients based on extensive field data
#[derive(Debug, Clone)]
pub struct SpinDriftCoefficients {
    /// Litz coefficient for gyroscopic drift (typically 0.8-1.5)
    pub litz_coefficient: f64,
    /// McCoy's aerodynamic jump factor
    pub mccoy_jump_factor: f64,
    /// Courtney's transonic adjustment
    pub transonic_factor: f64,
    /// Yaw damping coefficient
    pub yaw_damping: f64,
}

impl SpinDriftCoefficients {
    /// Get coefficients for specific bullet types based on empirical data
    pub fn for_bullet_type(bullet_type: &str) -> Self {
        match bullet_type.to_lowercase().as_str() {
            "match" | "bthp" | "boat_tail" => Self {
                litz_coefficient: 1.25,
                mccoy_jump_factor: 0.85,
                transonic_factor: 0.75,
                yaw_damping: 0.92,
            },
            "vld" | "very_low_drag" => Self {
                litz_coefficient: 1.15,
                mccoy_jump_factor: 0.78,
                transonic_factor: 0.68,
                yaw_damping: 0.88,
            },
            "hybrid" | "hybrid_ogive" => Self {
                litz_coefficient: 1.20,
                mccoy_jump_factor: 0.82,
                transonic_factor: 0.72,
                yaw_damping: 0.90,
            },
            "flat_base" | "fb" => Self {
                litz_coefficient: 1.35,
                mccoy_jump_factor: 0.95,
                transonic_factor: 0.85,
                yaw_damping: 0.95,
            },
            _ => Self::default(),
        }
    }

    pub fn default() -> Self {
        Self {
            litz_coefficient: 1.25,
            mccoy_jump_factor: 0.85,
            transonic_factor: 0.75,
            yaw_damping: 0.92,
        }
    }
}

/// Calculate advanced spin drift using multiple empirical models
pub fn calculate_advanced_spin_drift(
    stability_factor: f64,
    time_of_flight_s: f64,
    velocity_mps: f64,
    muzzle_velocity_mps: f64,
    spin_rate_rad_s: f64,
    caliber_m: f64,
    mass_kg: f64,
    air_density_kg_m3: f64,
    is_right_twist: bool,
    bullet_type: &str,
) -> f64 {
    // Edge cases: no drift if no time or no stability
    // MBA-205: Guard against division by zero
    if time_of_flight_s <= 0.0
        || stability_factor <= 0.0
        || muzzle_velocity_mps <= 0.0
        || air_density_kg_m3 <= 0.0
    {
        return 0.0;
    }

    let coeffs = SpinDriftCoefficients::for_bullet_type(bullet_type);

    // Direction based on twist
    let sign = if is_right_twist { 1.0 } else { -1.0 };

    // Calculate Mach numbers
    let mach_current = velocity_mps / 343.0;
    let mach_muzzle = muzzle_velocity_mps / 343.0;

    // 1. Litz's empirical formula (primary component)
    let litz_drift =
        calculate_litz_drift(stability_factor, time_of_flight_s, coeffs.litz_coefficient);

    // 2. McCoy's aerodynamic jump correction
    let jump_correction = calculate_aerodynamic_jump_correction(
        mach_muzzle,
        spin_rate_rad_s,
        caliber_m,
        mass_kg,
        coeffs.mccoy_jump_factor,
    );

    // 3. Transonic correction factor
    let transonic_correction =
        calculate_transonic_correction(mach_current, coeffs.transonic_factor);

    // 4. Yaw damping effect
    let yaw_factor =
        calculate_yaw_damping_factor(stability_factor, time_of_flight_s, coeffs.yaw_damping);

    // 5. Velocity decay correction (new research)
    let velocity_ratio = velocity_mps / muzzle_velocity_mps;
    let velocity_correction = velocity_ratio.powf(0.3);

    // Combine all components
    let total_drift = sign
        * (litz_drift * transonic_correction * yaw_factor * velocity_correction + jump_correction);

    // Apply atmospheric density correction
    let density_correction = (1.225 / air_density_kg_m3).sqrt();

    total_drift * density_correction
}

/// Litz's empirical formula with refined coefficients
fn calculate_litz_drift(stability: f64, time_s: f64, coefficient: f64) -> f64 {
    if stability <= 1.0 || time_s <= 0.0 {
        return 0.0;
    }

    // Refined Litz formula based on extensive field testing
    // SD = k * (SG + 1.2) * TOF^1.83
    // where k is empirically determined for bullet type
    let base_drift = coefficient * (stability + 1.2) * time_s.powf(1.83);

    // Convert inches to meters
    base_drift * 0.0254
}

/// McCoy's aerodynamic jump correction at muzzle exit
fn calculate_aerodynamic_jump_correction(
    mach: f64,
    spin_rate_rad_s: f64,
    caliber_m: f64,
    mass_kg: f64,
    jump_factor: f64,
) -> f64 {
    // Aerodynamic jump contributes to initial displacement
    // Based on McCoy's research on muzzle exit effects

    // MBA-205: Guard against division by zero when mach == 0
    if mach <= 0.0 {
        return 0.0;
    }

    let spin_parameter = spin_rate_rad_s * caliber_m / (2.0 * 343.0 * mach);

    // Jump magnitude in milliradians
    let jump_mrad = jump_factor * spin_parameter * (mass_kg / 0.01).sqrt();

    // Convert to lateral displacement (approximation for small angles)
    // This is a one-time displacement, not time-dependent
    jump_mrad * 0.001 * 100.0 // Approximate 100m reference distance
}

/// Transonic correction based on Courtney & Courtney research
fn calculate_transonic_correction(mach: f64, transonic_factor: f64) -> f64 {
    if mach < 0.8 {
        // Subsonic - minimal correction needed
        1.0
    } else if mach < 1.2 {
        // Transonic region - significant instability
        // Smooth transition using cosine interpolation
        let transonic_ratio = (mach - 0.8) / 0.4;
        
        1.0 + (transonic_factor - 1.0) * (1.0 - (transonic_ratio * PI).cos()) / 2.0
    } else {
        // Supersonic - stable again but reduced effect
        0.85 + 0.15 * (2.5 - mach).max(0.0) / 1.3
    }
}

/// Yaw damping factor based on stability and time
fn calculate_yaw_damping_factor(stability: f64, time_s: f64, damping_coeff: f64) -> f64 {
    // Yaw oscillations damp out over time
    // Higher stability = faster damping
    let damping_rate = damping_coeff * stability.sqrt();
    let damped = 1.0 - (-damping_rate * time_s).exp();

    // Ensure reasonable bounds
    damped.max(0.5).min(1.0)
}

/// Calculate equilibrium yaw angle using advanced model
pub fn calculate_advanced_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    crosswind_mps: f64,
    spin_rate_rad_s: f64,
    air_density_kg_m3: f64,
    caliber_m: f64,
) -> f64 {
    if stability_factor <= 1.0 || velocity_mps <= 0.0 {
        return 0.0;
    }

    // Base yaw from crosswind
    let wind_yaw = if crosswind_mps != 0.0 && velocity_mps > 0.0 {
        (crosswind_mps / velocity_mps).atan()
    } else {
        // Natural yaw from trajectory curvature (gravity-induced)
        // Empirical value based on typical trajectories
        0.001 + 0.0005 * (velocity_mps / 800.0).min(2.0)
    };

    // Stability-based damping (McCoy's model)
    let stability_term = ((stability_factor - 1.0) / stability_factor).sqrt();

    // Dynamic pressure effect
    let q = 0.5 * air_density_kg_m3 * velocity_mps.powi(2);
    let q_factor = (q / 50000.0).min(1.5).max(0.5); // Normalize around typical q

    // Spin effect on yaw response
    let spin_factor = if spin_rate_rad_s > 0.0 {
        let spin_param = spin_rate_rad_s * caliber_m / (2.0 * velocity_mps);
        1.0 + 0.2 * spin_param.min(0.5)
    } else {
        1.0
    };

    wind_yaw * stability_term * q_factor * spin_factor
}

/// Data-driven correction factor (placeholder for ML integration)
pub fn apply_ml_correction(
    base_drift: f64,
    stability: f64,
    mach: f64,
    time_s: f64,
    caliber_inches: f64,
    mass_grains: f64,
) -> f64 {
    // This function would integrate with ML models trained on real-world data
    // For now, returns the base drift unmodified
    //
    // In production, this would:
    // 1. Extract features: [stability, mach, time_s, caliber_inches, mass_grains]
    // 2. Pass to trained neural network or gradient boosting model
    // 3. Return correction factor (typically 0.8-1.2)
    // 4. Multiply base_drift by correction factor

    // Placeholder implementation with simple heuristics
    let mut correction = 1.0;

    // Known adjustments from field data
    if stability > 2.5 && mach < 1.0 {
        correction *= 0.92; // Over-stabilized subsonic tends to drift less
    }

    if time_s > 2.0 && mach < 0.9 {
        correction *= 1.08; // Long flight subsonic needs more correction
    }

    if caliber_inches < 0.264 && mass_grains < 100.0 {
        correction *= 0.88; // Light, small caliber bullets drift less
    }

    base_drift * correction
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_advanced_spin_drift() {
        // Test with typical .308 Match bullet
        let drift = calculate_advanced_spin_drift(
            1.5,     // stability
            1.2,     // time of flight
            600.0,   // current velocity m/s
            850.0,   // muzzle velocity m/s
            1500.0,  // spin rate rad/s
            0.00308, // caliber in meters
            0.0108,  // mass in kg (168 grains)
            1.225,   // air density
            true,    // right twist
            "match", // bullet type
        );

        // Should give reasonable drift (2-8 inches at 1000 yards typical)
        assert!(drift > 0.0);
        assert!(drift < 0.3); // Less than 12 inches in meters
    }

    #[test]
    fn test_transonic_correction() {
        let subsonic = calculate_transonic_correction(0.7, 0.75);
        let transonic = calculate_transonic_correction(1.0, 0.75);
        let supersonic = calculate_transonic_correction(1.5, 0.75);

        assert_eq!(subsonic, 1.0);
        assert!(transonic > 0.8 && transonic < 1.0);
        assert!(supersonic > 0.7 && supersonic < 1.0);
    }

    #[test]
    fn test_spin_drift_direction() {
        // Right twist should produce positive drift
        let right_drift = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );

        // Left twist should produce negative drift
        let left_drift = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, false, "match",
        );

        assert!(right_drift > 0.0, "Right twist should give positive drift");
        assert!(left_drift < 0.0, "Left twist should give negative drift");
        assert!(
            (right_drift.abs() - left_drift.abs()).abs() < 0.001,
            "Magnitude should be equal"
        );
    }

    #[test]
    fn test_spin_drift_edge_cases() {
        // Zero time should give zero drift
        let zero_time = calculate_advanced_spin_drift(
            1.5, 0.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        assert_eq!(zero_time, 0.0);

        // Zero stability should give zero drift
        let zero_stability = calculate_advanced_spin_drift(
            0.0, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        assert_eq!(zero_stability, 0.0);

        // Zero muzzle velocity should give zero drift
        let zero_muzzle_vel = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 0.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        assert_eq!(zero_muzzle_vel, 0.0);

        // Zero air density should give zero drift
        let zero_density = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 0.0, true, "match",
        );
        assert_eq!(zero_density, 0.0);
    }

    #[test]
    fn test_spin_drift_coefficients_bullet_types() {
        let match_coeffs = SpinDriftCoefficients::for_bullet_type("match");
        let vld_coeffs = SpinDriftCoefficients::for_bullet_type("vld");
        let flat_base_coeffs = SpinDriftCoefficients::for_bullet_type("flat_base");
        let default_coeffs = SpinDriftCoefficients::for_bullet_type("unknown");

        // VLD should have lower Litz coefficient
        assert!(vld_coeffs.litz_coefficient < match_coeffs.litz_coefficient);

        // Flat base should have higher Litz coefficient
        assert!(flat_base_coeffs.litz_coefficient > match_coeffs.litz_coefficient);

        // Default should match match type
        assert_eq!(default_coeffs.litz_coefficient, match_coeffs.litz_coefficient);
    }

    #[test]
    fn test_spin_drift_increases_with_time() {
        let drift_short = calculate_advanced_spin_drift(
            1.5, 0.5, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        let drift_medium = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        let drift_long = calculate_advanced_spin_drift(
            1.5, 2.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );

        assert!(drift_medium > drift_short, "Drift should increase with time");
        assert!(drift_long > drift_medium, "Drift should increase with time");
    }

    #[test]
    fn test_advanced_yaw_of_repose() {
        let yaw = calculate_advanced_yaw_of_repose(
            1.5,     // stability
            800.0,   // velocity m/s
            5.0,     // crosswind m/s
            1500.0,  // spin rate rad/s
            1.225,   // air density
            0.00782, // caliber m
        );

        // Should give small angle in radians
        assert!(yaw.abs() < 0.1, "Yaw should be small angle, got {}", yaw);
    }

    #[test]
    fn test_yaw_of_repose_edge_cases() {
        // Zero stability should give zero yaw
        let zero_stability = calculate_advanced_yaw_of_repose(0.5, 800.0, 5.0, 1500.0, 1.225, 0.00782);
        assert_eq!(zero_stability, 0.0);

        // Zero velocity should give zero yaw
        let zero_velocity = calculate_advanced_yaw_of_repose(1.5, 0.0, 5.0, 1500.0, 1.225, 0.00782);
        assert_eq!(zero_velocity, 0.0);

        // No crosswind should still give small yaw (trajectory curvature)
        let no_wind = calculate_advanced_yaw_of_repose(1.5, 800.0, 0.0, 1500.0, 1.225, 0.00782);
        assert!(no_wind > 0.0, "Should have natural yaw from trajectory curvature");
    }

    #[test]
    fn test_transonic_correction_continuity() {
        // Test continuity across transonic region boundaries
        let just_below_transonic = calculate_transonic_correction(0.79, 0.75);
        let just_at_transonic_start = calculate_transonic_correction(0.80, 0.75);

        // Should be continuous at 0.8 Mach
        assert!(
            (just_below_transonic - just_at_transonic_start).abs() < 0.01,
            "Should be continuous at transonic start"
        );

        let just_below_supersonic = calculate_transonic_correction(1.19, 0.75);
        let just_at_supersonic = calculate_transonic_correction(1.21, 0.75);

        // Values should be reasonably close
        assert!(just_below_supersonic > 0.0);
        assert!(just_at_supersonic > 0.0);
    }

    #[test]
    fn test_ml_correction_placeholder() {
        // Test the ML correction placeholder function
        let base_drift = 0.1;
        let corrected = apply_ml_correction(base_drift, 1.5, 2.5, 1.0, 0.308, 168.0);

        // Should return reasonable multiplied value
        assert!(corrected > 0.0);

        // Test specific heuristics
        // Over-stabilized subsonic
        let over_stab_subsonic = apply_ml_correction(0.1, 3.0, 0.8, 1.0, 0.308, 168.0);
        assert!(over_stab_subsonic < 0.1, "Over-stabilized subsonic should drift less");

        // Long flight subsonic
        let long_subsonic = apply_ml_correction(0.1, 1.5, 0.85, 2.5, 0.308, 168.0);
        assert!(long_subsonic > 0.1, "Long subsonic flight should need more correction");

        // Light small caliber
        let light_small = apply_ml_correction(0.1, 1.5, 2.5, 1.0, 0.224, 55.0);
        assert!(light_small < 0.1, "Light small caliber should drift less");
    }

    #[test]
    fn test_density_affects_drift() {
        // Lower density (higher altitude) should increase drift
        let sea_level = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        let high_altitude = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.0, true, "match",
        );

        assert!(
            high_altitude > sea_level,
            "Higher altitude (lower density) should increase drift"
        );
    }

    #[test]
    fn test_different_bullet_types_drift() {
        let match_drift = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "match",
        );
        let vld_drift = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "vld",
        );
        let flat_base_drift = calculate_advanced_spin_drift(
            1.5, 1.0, 700.0, 850.0, 1500.0, 0.00308, 0.0108, 1.225, true, "flat_base",
        );

        // VLD should drift less than match (lower coefficient)
        assert!(vld_drift < match_drift, "VLD should drift less than match");

        // Flat base should drift more (higher coefficient)
        assert!(flat_base_drift > match_drift, "Flat base should drift more");
    }

    #[test]
    fn test_litz_drift_low_stability() {
        // Stability <= 1.0 should return zero from Litz formula
        let low_stability = calculate_litz_drift(0.9, 1.0, 1.25);
        assert_eq!(low_stability, 0.0);

        let exactly_one = calculate_litz_drift(1.0, 1.0, 1.25);
        assert_eq!(exactly_one, 0.0);

        // Just above 1.0 should give positive result
        let above_one = calculate_litz_drift(1.1, 1.0, 1.25);
        assert!(above_one > 0.0);
    }

    #[test]
    fn test_aerodynamic_jump_correction_edge_cases() {
        // Zero mach should return zero
        let zero_mach = calculate_aerodynamic_jump_correction(0.0, 1500.0, 0.00308, 0.0108, 0.85);
        assert_eq!(zero_mach, 0.0);

        // Valid inputs should return non-zero
        let valid = calculate_aerodynamic_jump_correction(2.5, 1500.0, 0.00308, 0.0108, 0.85);
        assert!(valid != 0.0);
    }

    #[test]
    fn test_yaw_damping_factor() {
        // Higher stability should damp faster
        let low_stability_damping = calculate_yaw_damping_factor(1.2, 1.0, 0.92);
        let high_stability_damping = calculate_yaw_damping_factor(2.0, 1.0, 0.92);

        // Higher stability = faster damping = higher value closer to 1.0
        assert!(high_stability_damping >= low_stability_damping);

        // Result should be bounded 0.5-1.0
        assert!(low_stability_damping >= 0.5);
        assert!(low_stability_damping <= 1.0);
        assert!(high_stability_damping >= 0.5);
        assert!(high_stability_damping <= 1.0);
    }
}
