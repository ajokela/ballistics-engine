use crate::InternalBallisticInputs as BallisticInputs;

/// Calculate the gyroscopic stability coefficient (SG) for the bullet.
///
/// This function uses the Miller stability formula. An SG value greater than 1.5
/// is generally considered to indicate adequate stability.
///
/// # Arguments
/// * `inputs` - Ballistic input parameters
/// * `atmo_params` - Atmospheric parameters (altitude, temp_c, pressure_hpa, density_ratio)
///
/// # Returns
/// * Stability coefficient (dimensionless)
pub fn compute_stability_coefficient(
    inputs: &BallisticInputs,
    atmo_params: (f64, f64, f64, f64),
) -> f64 {
    // Check for required parameters
    if inputs.twist_rate == 0.0 || inputs.bullet_length == 0.0 || inputs.bullet_diameter == 0.0 {
        return 0.0;
    }

    // Pre-calculated constants for efficiency
    const MILLER_CONST: f64 = 30.0;
    const VEL_REF_FPS: f64 = 2800.0;
    const TEMP_REF_K: f64 = 288.15; // 15°C
    const PRESS_REF_HPA: f64 = 1013.25;

    // Calculate intermediate values
    // Convert twist rate from inches to meters for consistency
    let twist_rate_m = inputs.twist_rate.abs() * 0.0254; // inches to meters
    let twist_calibers = twist_rate_m / inputs.bullet_diameter;
    let length_calibers = inputs.bullet_length / inputs.bullet_diameter;

    // Convert units for Miller formula
    let mass_grains = inputs.bullet_mass / 0.00006479891; // kg to grains
    let diameter_inches = inputs.bullet_diameter / 0.0254; // meters to inches
    let velocity_fps = inputs.muzzle_velocity * 3.28084; // m/s to fps

    // Miller stability formula components
    let mass_term = MILLER_CONST * mass_grains;
    let geom_term = twist_calibers.powi(2)
        * diameter_inches.powi(3)
        * length_calibers
        * (1.0 + length_calibers.powi(2));

    if geom_term == 0.0 {
        return 0.0;
    }

    // Extract atmospheric parameters
    let (_, temp_c, current_press_hpa, _) = atmo_params;
    let temp_k = temp_c + 273.15;

    // Atmospheric density correction factor
    // Ratio of reference density to current density
    let density_correction = (temp_k / TEMP_REF_K) * (PRESS_REF_HPA / current_press_hpa);

    // Velocity correction factor
    let velocity_correction = (velocity_fps / VEL_REF_FPS).powf(1.0 / 3.0);

    // Final stability calculation

    (mass_term / geom_term) * velocity_correction * density_correction
}

/// Calculate spin drift in meters using Litz approximation.
///
/// # Arguments
/// * `time_s` - Time of flight in seconds
/// * `stability` - Stability coefficient
/// * `twist_rate` - Twist rate in inches (calibers per turn)
/// * `is_twist_right` - True for right-hand twist, false for left-hand
///
/// # Returns
/// * Spin drift in meters
pub fn compute_spin_drift(
    time_s: f64,
    stability: f64,
    twist_rate: f64,
    is_twist_right: bool,
) -> f64 {
    compute_spin_drift_with_decay(time_s, stability, twist_rate, is_twist_right, None)
}

/// Calculate spin drift with optional spin decay modeling
pub fn compute_spin_drift_with_decay(
    time_s: f64,
    stability: f64,
    twist_rate: f64,
    is_twist_right: bool,
    decay_factor: Option<f64>, // Optional spin decay factor (0-1)
) -> f64 {
    if stability == 0.0 || time_s <= 0.0 || twist_rate == 0.0 {
        return 0.0;
    }

    let sign = if is_twist_right { 1.0 } else { -1.0 };

    // Litz empirical spin drift: inches = 1.25 * (SG + 1.2) * TOF^1.83.
    // Keep this summary/API path consistent with cli_api::apply_spin_drift.
    let scaling_factor = 1.25;
    let base_drift = sign * scaling_factor * (stability + 1.2) * time_s.powf(1.83);

    // Apply spin decay if provided
    let effective_drift = if let Some(decay) = decay_factor {
        base_drift * decay.max(0.0).min(1.0)
    } else {
        base_drift
    };

    // Convert inches to meters
    effective_drift * 0.0254
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_inputs() -> BallisticInputs {
        BallisticInputs {
            muzzle_velocity: 823.0, // 2700 fps in m/s
            bc_value: 0.5,
            bullet_mass: 0.0109,      // 168 grains in kg
            bullet_diameter: 0.00782, // 0.308 inches in meters
            bullet_length: 0.033,     // in meters (1.3 inches)
            twist_rate: 10.0,
            ..Default::default()
        }
    }

    #[test]
    fn test_compute_stability_coefficient() {
        let inputs = create_test_inputs();
        let atmo_params = (0.0, 15.0, 1013.25, 1.0); // Standard conditions

        let stability = compute_stability_coefficient(&inputs, atmo_params);

        // Debug output
        println!("Computed stability: {}", stability);

        // Should be a reasonable stability value
        assert!(stability > 0.0);
        assert!(stability < 10.0); // Sanity check

        // Test with typical values should give stability around 1.5-2.5
        assert!(stability > 1.0);
        assert!(stability < 3.0);
    }

    #[test]
    fn test_compute_stability_coefficient_zero_cases() {
        let mut inputs = create_test_inputs();
        let atmo_params = (0.0, 15.0, 1013.25, 1.0);

        // Test with zero twist rate
        inputs.twist_rate = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);

        // Test with zero bullet length
        inputs = create_test_inputs();
        inputs.bullet_length = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);

        // Test with zero bullet diameter
        inputs = create_test_inputs();
        inputs.bullet_diameter = 0.0;
        assert_eq!(compute_stability_coefficient(&inputs, atmo_params), 0.0);
    }

    #[test]
    fn test_compute_stability_coefficient_atmospheric_effects() {
        let inputs = create_test_inputs();

        // Standard conditions
        let standard_atmo = (0.0, 15.0, 1013.25, 1.0);
        let standard_stability = compute_stability_coefficient(&inputs, standard_atmo);

        // High altitude (lower pressure, lower temperature)
        let high_alt_atmo = (3000.0, 5.0, 900.0, 1.0);
        let high_alt_stability = compute_stability_coefficient(&inputs, high_alt_atmo);

        // High altitude should have higher stability due to lower air density
        assert!(high_alt_stability > standard_stability);

        // Hot conditions (higher temperature)
        let hot_atmo = (0.0, 35.0, 1013.25, 1.0);
        let hot_stability = compute_stability_coefficient(&inputs, hot_atmo);

        // Hot conditions should have higher stability due to lower air density
        assert!(hot_stability > standard_stability);
    }

    #[test]
    fn test_compute_spin_drift() {
        let time_s = 1.5;
        let stability = 2.0;
        let twist_rate = 10.0;

        // Test right-hand twist
        let drift_right = compute_spin_drift(time_s, stability, twist_rate, true);
        assert!(drift_right > 0.0); // Should drift to the right (positive)

        // Test left-hand twist
        let drift_left = compute_spin_drift(time_s, stability, twist_rate, false);
        assert!(drift_left < 0.0); // Should drift to the left (negative)
        assert!((drift_left + drift_right).abs() < 1e-10); // Should be equal magnitude

        // Litz spin drift remains a small correction for this flight time.
        assert!(drift_right.abs() < 0.25); // Less than 25cm for 1.5s flight
    }

    #[test]
    fn test_compute_spin_drift_zero_cases() {
        // Test with zero stability
        assert_eq!(compute_spin_drift(1.5, 0.0, 10.0, true), 0.0);

        // Test with zero time
        assert_eq!(compute_spin_drift(0.0, 2.0, 10.0, true), 0.0);

        // Test with negative time
        assert_eq!(compute_spin_drift(-1.0, 2.0, 10.0, true), 0.0);

        // Test with zero twist rate
        assert_eq!(compute_spin_drift(1.5, 2.0, 0.0, true), 0.0);
    }

    #[test]
    fn test_compute_spin_drift_scaling() {
        let stability = 2.0;
        let twist_rate = 10.0;

        // Test time scaling
        let drift_1s = compute_spin_drift(1.0, stability, twist_rate, true);
        let drift_2s = compute_spin_drift(2.0, stability, twist_rate, true);

        // Drift should increase with time (non-linearly due to 1.83 exponent)
        assert!(drift_2s > drift_1s);

        // Test stability scaling
        let drift_low_stability = compute_spin_drift(1.5, 1.0, twist_rate, true);
        let drift_high_stability = compute_spin_drift(1.5, 3.0, twist_rate, true);

        // Higher stability should produce more drift
        assert!(drift_high_stability > drift_low_stability);
    }
}
