//! Spin Decay Physics for Ballistics Calculations
//!
//! This module implements realistic spin decay modeling based on:
//! - Aerodynamic torque opposing spin
//! - Calibrated roll-damping coefficients
//! - Velocity-dependent decay rates
//! - Projectile shape and surface-finish presets

use std::f64::consts::PI;

const MATCH_REFERENCE_DECAY_RATE_PER_SECOND: f64 = 0.025;
const GENERAL_REFERENCE_DECAY_RATE_PER_SECOND: f64 = 0.04;

/// Parameters affecting spin decay rate
#[derive(Debug, Clone, Copy)]
pub struct SpinDecayParameters {
    /// Surface roughness in meters (typical: 0.1mm)
    pub surface_roughness: f64,
    /// Effective dimensionless roll-damping coefficient magnitude (`|C_lp|`) before applying
    /// `form_factor`. The field name is retained for API compatibility.
    pub skin_friction_coefficient: f64,
    /// Shape factor for spin damping
    pub form_factor: f64,
}

impl SpinDecayParameters {
    /// Create default parameters
    pub fn new() -> Self {
        // The effective coefficient (coefficient * form factor) is calibrated against the
        // empirical 4%/s reference decay used by update_spin_rate for a 175gr .308 ogive.
        Self {
            surface_roughness: 0.0001,
            skin_friction_coefficient: 0.00363,
            form_factor: 1.0,
        }
    }

    /// Get typical parameters for different bullet types
    pub fn from_bullet_type(bullet_type: &str) -> Self {
        // Match bullets calibrate to the empirical 2.5%/s reference; the other presets calibrate
        // to 4%/s after their shape form factor is applied.
        match bullet_type.to_lowercase().as_str() {
            "match" => Self {
                surface_roughness: 0.00005,
                skin_friction_coefficient: 0.00252,
                form_factor: 0.9,
            },
            "hunting" => Self {
                surface_roughness: 0.0001,
                skin_friction_coefficient: 0.00363,
                form_factor: 1.0,
            },
            "fmj" => Self {
                surface_roughness: 0.00015,
                skin_friction_coefficient: 0.00330,
                form_factor: 1.1,
            },
            "cast" => Self {
                surface_roughness: 0.0002,
                skin_friction_coefficient: 0.00303,
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

/// Calculate the magnitude of the aerodynamic moment opposing spin.
///
/// Uses the conventional roll-damping relation
/// `M = 1/4 * rho * V * S * d^2 * |C_lp| * |p|`, where `S` is projectile reference area.
/// Here `C_lp` is the derivative with respect to reduced spin `p*d/(2*V)`; conventions using
/// `p*d/V` report a coefficient smaller by a factor of two.
/// The returned value is nonnegative; [`calculate_spin_decay_rate`] applies the direction that
/// opposes the signed spin rate.
pub fn calculate_spin_damping_moment(
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    caliber_m: f64,
    length_m: f64,
    decay_params: &SpinDecayParameters,
) -> f64 {
    if !spin_rate_rad_s.is_finite()
        || spin_rate_rad_s == 0.0
        || !velocity_mps.is_finite()
        || velocity_mps <= 0.0
        || !air_density_kg_m3.is_finite()
        || air_density_kg_m3 <= 0.0
        || !caliber_m.is_finite()
        || caliber_m <= 0.0
        || !length_m.is_finite()
        || length_m <= 0.0
        || !decay_params.skin_friction_coefficient.is_finite()
        || decay_params.skin_friction_coefficient <= 0.0
        || !decay_params.form_factor.is_finite()
        || decay_params.form_factor <= 0.0
    {
        return 0.0;
    }

    let reference_area = PI * (caliber_m / 2.0).powi(2);
    let roll_damping_coefficient =
        decay_params.skin_friction_coefficient * decay_params.form_factor;

    0.25 * air_density_kg_m3
        * velocity_mps
        * reference_area
        * caliber_m.powi(2)
        * roll_damping_coefficient
        * spin_rate_rad_s.abs()
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
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
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

    // Apply the damping-moment magnitude opposite to either signed spin direction.
    if moment_of_inertia > 0.0 && spin_rate_rad_s.is_finite() {
        -spin_rate_rad_s.signum() * damping_moment / moment_of_inertia
    } else {
        0.0
    }
}

/// Calculate the spin rate after accounting for decay
///
/// Uses an empirical model based on published ballistics data.
/// Real bullets typically lose 5-15% of spin over a 3-second flight.
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
pub fn update_spin_rate(
    initial_spin_rad_s: f64,
    time_elapsed_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    mass_grains: f64,
    caliber_inches: f64,
    length_inches: f64,
    decay_params: Option<&SpinDecayParameters>,
) -> f64 {
    if time_elapsed_s <= 0.0 {
        return initial_spin_rad_s;
    }

    // Mass factor (heavier bullets retain spin better)
    let mass_factor = (175.0 / mass_grains).sqrt(); // Normalized to 175gr

    // Velocity factor (higher velocity means more decay)
    let velocity_factor = velocity_mps / 850.0; // Normalized to 850 m/s

    // Air density scaling: higher density = more aerodynamic damping
    // Normalized to standard sea-level density (1.225 kg/m^3)
    let density_factor = if air_density_kg_m3 > 0.0 {
        air_density_kg_m3 / 1.225
    } else {
        1.0 // Standard conditions fallback
    };

    // Surface area factor from caliber and length
    // Larger surface area relative to a reference .308 bullet = more skin-friction damping.
    // Use sqrt of the ratio: skin-friction torque scales with surface area but rotational
    // inertia also grows with size, so the net effect on decay rate is sub-linear.
    // Reference: .308 caliber (0.308") x 1.3" length
    let ref_surface = PI * 0.308 * 1.3; // reference lateral surface area (inches^2)
    let surface_factor = if caliber_inches > 0.0 && length_inches > 0.0 {
        let bullet_surface = PI * caliber_inches * length_inches;
        (bullet_surface / ref_surface).sqrt()
    } else {
        1.0 // Standard conditions fallback
    };

    // Base decay rate per second (empirical)
    let base_decay_rate = if let Some(params) = decay_params {
        if params.form_factor < 1.0 {
            // Match bullet
            MATCH_REFERENCE_DECAY_RATE_PER_SECOND
        } else {
            // Hunting/FMJ bullet
            GENERAL_REFERENCE_DECAY_RATE_PER_SECOND
        }
    } else {
        GENERAL_REFERENCE_DECAY_RATE_PER_SECOND
    };

    // Adjusted decay rate blending empirical model with physical parameters.
    // At standard conditions (sea-level density, .308 reference bullet) the
    // density_factor and surface_factor are both 1.0, preserving legacy behavior.
    let decay_rate_per_second =
        base_decay_rate * mass_factor * velocity_factor * density_factor * surface_factor;

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
            3.0,   // 3 seconds
            750.0, // 750 m/s average velocity
            1.2,   // air density
            175.0, // 175 grains
            0.308, // caliber
            1.3,   // length
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

    #[test]
    fn test_spin_damping_moment() {
        let params = SpinDecayParameters::from_bullet_type("match");

        // Test with typical values
        let moment = calculate_spin_damping_moment(
            1000.0,  // spin rate rad/s
            800.0,   // velocity m/s
            1.225,   // air density
            0.00782, // caliber in meters (.308")
            0.033,   // length in meters
            &params,
        );

        // Moment should be positive (opposes spin)
        assert!(moment > 0.0);
        assert!(moment < 1.0); // Should be small for a bullet

        // Test zero spin
        let zero_moment = calculate_spin_damping_moment(0.0, 800.0, 1.225, 0.00782, 0.033, &params);
        assert_eq!(zero_moment, 0.0);

        // Test zero velocity
        let zero_vel_moment =
            calculate_spin_damping_moment(1000.0, 0.0, 1.225, 0.00782, 0.033, &params);
        assert_eq!(zero_vel_moment, 0.0);
    }

    #[test]
    fn test_spin_decay_rate() {
        let params = SpinDecayParameters::from_bullet_type("fmj");

        let decay_rate = calculate_spin_decay_rate(
            1000.0, // spin rate rad/s
            800.0,  // velocity m/s
            1.225,  // air density
            168.0,  // mass grains
            0.308,  // caliber inches
            1.2,    // length inches
            &params,
            "boat_tail",
        );

        // Decay rate should be negative (spin decreases)
        assert!(decay_rate < 0.0);
        assert!(decay_rate > -1000.0); // Should be reasonable magnitude
    }

    #[test]
    fn physical_spin_decay_matches_empirical_reference_rate() {
        let spin_rate = 17_000.0;
        let velocity = 800.0;
        let cases = [
            SpinDecayParameters::from_bullet_type("match"),
            SpinDecayParameters::from_bullet_type("hunting"),
        ];

        for params in cases {
            let actual_rate = calculate_spin_decay_rate(
                spin_rate, velocity, 1.225, 175.0, 0.308, 1.3, &params, "ogive",
            );
            let empirical_spin_after_one_second = update_spin_rate(
                spin_rate,
                1.0,
                velocity,
                1.225,
                175.0,
                0.308,
                1.3,
                Some(&params),
            );
            let expected_rate = (empirical_spin_after_one_second / spin_rate).ln() * spin_rate;

            assert!(
                (actual_rate - expected_rate).abs() <= expected_rate.abs() * 0.05,
                "physical rate {actual_rate} did not match empirical reference {expected_rate}"
            );
        }
    }

    #[test]
    fn roll_damping_uses_canonical_reduced_spin_moment() {
        let params = SpinDecayParameters {
            surface_roughness: 0.0,
            skin_friction_coefficient: 0.01,
            form_factor: 1.0,
        };
        let moment = calculate_spin_damping_moment(
            17_000.0,
            800.0,
            1.225,
            0.308 * 0.0254,
            1.3 * 0.0254,
            &params,
        );
        let expected = 1.225_300_524_995_314e-4;

        assert!((moment - expected).abs() <= expected * 1e-12);
    }

    #[test]
    fn roll_damping_moment_has_physical_scaling() {
        let params = SpinDecayParameters::from_bullet_type("match");
        let caliber = 0.308 * 0.0254;
        let length = 1.3 * 0.0254;
        let moment = |spin, velocity, density, diameter, projectile_length| {
            calculate_spin_damping_moment(
                spin,
                velocity,
                density,
                diameter,
                projectile_length,
                &params,
            )
        };
        let base = moment(17_000.0, 800.0, 1.225, caliber, length);
        let doubled_coefficient_params = SpinDecayParameters {
            skin_friction_coefficient: 2.0 * params.skin_friction_coefficient,
            ..params
        };
        let doubled_form_factor_params = SpinDecayParameters {
            form_factor: 2.0 * params.form_factor,
            ..params
        };
        let cases = [
            ("spin", moment(34_000.0, 800.0, 1.225, caliber, length), 2.0),
            (
                "velocity",
                moment(17_000.0, 1_600.0, 1.225, caliber, length),
                2.0,
            ),
            (
                "density",
                moment(17_000.0, 800.0, 2.45, caliber, length),
                2.0,
            ),
            (
                "caliber",
                moment(17_000.0, 800.0, 1.225, 2.0 * caliber, length),
                16.0,
            ),
            (
                "length",
                moment(17_000.0, 800.0, 1.225, caliber, 2.0 * length),
                1.0,
            ),
            (
                "coefficient",
                calculate_spin_damping_moment(
                    17_000.0,
                    800.0,
                    1.225,
                    caliber,
                    length,
                    &doubled_coefficient_params,
                ),
                2.0,
            ),
            (
                "form factor",
                calculate_spin_damping_moment(
                    17_000.0,
                    800.0,
                    1.225,
                    caliber,
                    length,
                    &doubled_form_factor_params,
                ),
                2.0,
            ),
        ];

        for (name, moment, expected_ratio) in cases {
            let actual_ratio = moment / base;
            assert!(
                (actual_ratio - expected_ratio).abs() <= expected_ratio * 1e-12,
                "{name} scaling was {actual_ratio}, expected {expected_ratio}"
            );
        }
    }

    #[test]
    fn spin_decay_always_opposes_spin_direction() {
        let params = SpinDecayParameters::from_bullet_type("match");
        let positive =
            calculate_spin_decay_rate(17_000.0, 800.0, 1.225, 175.0, 0.308, 1.3, &params, "ogive");
        let negative =
            calculate_spin_decay_rate(-17_000.0, 800.0, 1.225, 175.0, 0.308, 1.3, &params, "ogive");

        assert!(positive < 0.0);
        assert!(negative > 0.0);
        assert!((positive + negative).abs() <= positive.abs() * 1e-12);
    }

    #[test]
    fn roll_damping_rejects_nonphysical_inputs() {
        let params = SpinDecayParameters::from_bullet_type("match");
        let invalid_states = [
            (0.0, 800.0, 1.225, 0.00782, 0.033),
            (17_000.0, 0.0, 1.225, 0.00782, 0.033),
            (17_000.0, -800.0, 1.225, 0.00782, 0.033),
            (17_000.0, 800.0, 0.0, 0.00782, 0.033),
            (17_000.0, 800.0, -1.225, 0.00782, 0.033),
            (17_000.0, 800.0, 1.225, 0.0, 0.033),
            (17_000.0, 800.0, 1.225, -0.00782, 0.033),
            (17_000.0, 800.0, 1.225, 0.00782, 0.0),
            (17_000.0, 800.0, 1.225, 0.00782, -0.033),
        ];

        for (spin, velocity, density, caliber, length) in invalid_states {
            let moment =
                calculate_spin_damping_moment(spin, velocity, density, caliber, length, &params);
            assert_eq!(
                moment, 0.0,
                "nonphysical state produced damping moment {moment}"
            );
        }
    }

    #[test]
    fn test_different_bullet_types() {
        // Test all bullet type parameters
        let types = ["match", "hunting", "fmj", "cast", "unknown"];

        for bullet_type in &types {
            let params = SpinDecayParameters::from_bullet_type(bullet_type);
            assert!(params.surface_roughness > 0.0);
            assert!(params.skin_friction_coefficient > 0.0);
            assert!(params.form_factor > 0.0);
        }
    }

    #[test]
    fn test_moment_of_inertia_shapes() {
        let mass_kg = 0.01;
        let caliber_m = 0.008;
        let length_m = 0.03;

        let i_cylinder = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "cylinder");
        let i_ogive = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");
        let i_boat_tail = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "boat_tail");
        let i_default = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "unknown");

        // Check relative magnitudes
        assert!(i_cylinder > i_ogive);
        assert!(i_ogive > i_boat_tail);
        assert_eq!(i_cylinder, i_default); // Unknown defaults to cylinder

        // Check absolute values are reasonable
        assert!(i_cylinder > 0.0);
        assert!(i_boat_tail > 0.0);
    }

    #[test]
    fn test_spin_decay_correction_factor() {
        let params = SpinDecayParameters::from_bullet_type("match");

        // At time zero, factor should be 1.0
        let factor_t0 = calculate_spin_decay_correction_factor(
            0.0,
            800.0,
            1.225,
            175.0,
            0.308,
            1.3,
            Some(&params),
        );
        assert_eq!(factor_t0, 1.0);

        // After some time, factor should be less than 1.0 but greater than 0.5
        let factor_t3 = calculate_spin_decay_correction_factor(
            3.0,
            800.0,
            1.225,
            175.0,
            0.308,
            1.3,
            Some(&params),
        );
        assert!(factor_t3 < 1.0);
        assert!(factor_t3 > 0.5);

        // Factor should decrease with time
        let factor_t1 = calculate_spin_decay_correction_factor(
            1.0,
            800.0,
            1.225,
            175.0,
            0.308,
            1.3,
            Some(&params),
        );
        let factor_t2 = calculate_spin_decay_correction_factor(
            2.0,
            800.0,
            1.225,
            175.0,
            0.308,
            1.3,
            Some(&params),
        );
        assert!(factor_t1 > factor_t2);
        assert!(factor_t2 > factor_t3);
    }

    #[test]
    fn test_default_impl() {
        let params1 = SpinDecayParameters::new();
        let params2 = SpinDecayParameters::default();

        assert_eq!(params1.surface_roughness, params2.surface_roughness);
        assert_eq!(
            params1.skin_friction_coefficient,
            params2.skin_friction_coefficient
        );
        assert_eq!(params1.form_factor, params2.form_factor);
    }

    #[test]
    fn test_mass_factor_effects() {
        let params = SpinDecayParameters::from_bullet_type("match");

        // Light bullet (55gr)
        let spin_light =
            update_spin_rate(1000.0, 2.0, 800.0, 1.225, 55.0, 0.224, 0.9, Some(&params));

        // Heavy bullet (300gr)
        let spin_heavy =
            update_spin_rate(1000.0, 2.0, 800.0, 1.225, 300.0, 0.338, 1.8, Some(&params));

        // Heavy bullet should retain more spin (decay less)
        assert!(spin_heavy > spin_light);
    }

    #[test]
    fn test_velocity_factor_effects() {
        let params = SpinDecayParameters::from_bullet_type("hunting");

        // Low velocity
        let spin_low_vel =
            update_spin_rate(1000.0, 2.0, 400.0, 1.225, 175.0, 0.308, 1.3, Some(&params));

        // High velocity
        let spin_high_vel =
            update_spin_rate(1000.0, 2.0, 1200.0, 1.225, 175.0, 0.308, 1.3, Some(&params));

        // Higher velocity should cause more decay (less spin remaining)
        assert!(spin_low_vel > spin_high_vel);
    }
}
