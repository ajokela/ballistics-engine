// Advanced spin drift model based on modern ballistics research
// Incorporates multiple empirical models from:
// - Bryan Litz's Applied Ballistics for Long Range Shooting
// - McCoy's Modern Exterior Ballistics
// - Courtney & Courtney spin drift research papers

/// Legacy experimental coefficients retained for API compatibility.
///
/// [`calculate_advanced_spin_drift`] no longer applies these values because doing so would
/// reweight the already calibrated Litz total-drift fit.
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

    /// Return the legacy default coefficient set.
    ///
    /// This inherent constructor is retained in addition to [`Default`] for source compatibility
    /// with callers that invoke it through a fully qualified inherent-method path.
    #[allow(clippy::should_implement_trait)] // The trait is implemented below; this preserves API.
    pub fn default() -> Self {
        <Self as Default>::default()
    }
}

impl Default for SpinDriftCoefficients {
    fn default() -> Self {
        Self {
            litz_coefficient: 1.25,
            mccoy_jump_factor: 0.85,
            transonic_factor: 0.75,
            yaw_damping: 0.92,
        }
    }
}

/// Calculate signed spin drift using the calibrated Litz total-drift fit.
///
/// The Litz `1.25 * (Sg + 1.2) * TOF^1.83` relation already fits total observed drift, including
/// the effects of velocity loss and yaw damping. The additional state arguments are retained for
/// API compatibility, but they must not be stacked onto the fitted total as independent
/// multipliers. Atmospheric effects belong in the supplied `stability_factor` and time of flight.
/// Muzzle velocity and density must remain finite and positive solely to preserve the legacy
/// invalid-input contract; within that valid domain they do not rescale the result.
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
pub fn calculate_advanced_spin_drift(
    stability_factor: f64,
    time_of_flight_s: f64,
    _velocity_mps: f64,
    muzzle_velocity_mps: f64,
    _spin_rate_rad_s: f64,
    _caliber_m: f64,
    _mass_kg: f64,
    air_density_kg_m3: f64,
    is_right_twist: bool,
    _bullet_type: &str,
) -> f64 {
    if !stability_factor.is_finite()
        || stability_factor <= 1.0
        || !time_of_flight_s.is_finite()
        || time_of_flight_s <= 0.0
        || !muzzle_velocity_mps.is_finite()
        || muzzle_velocity_mps <= 0.0
        || !air_density_kg_m3.is_finite()
        || air_density_kg_m3 <= 0.0
    {
        return 0.0;
    }

    crate::spin_drift::litz_drift_meters(stability_factor, time_of_flight_s, is_right_twist)
}

/// Estimate flat-fire yaw of repose for the representative projectile in
/// [`crate::precession_nutation::PrecessionNutationParams::default`].
///
/// Retained for source compatibility. Crosswind is a transient and is ignored; density must
/// already be represented in the supplied local `stability_factor`; and caliber alone cannot
/// determine the missing inertia ratio. Use
/// [`crate::precession_nutation::calculate_limit_cycle_yaw_with_inertias`] for a
/// projectile-specific result.
#[deprecated(
    since = "0.22.18",
    note = "use precession_nutation::calculate_limit_cycle_yaw_with_inertias"
)]
pub fn calculate_advanced_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    _crosswind_mps: f64,
    spin_rate_rad_s: f64,
    _air_density_kg_m3: f64,
    _caliber_m: f64,
) -> f64 {
    let reference = crate::precession_nutation::PrecessionNutationParams::default();
    crate::precession_nutation::calculate_limit_cycle_yaw_with_inertias(
        velocity_mps,
        spin_rate_rad_s,
        stability_factor,
        reference.spin_inertia,
        reference.transverse_inertia,
    )
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
    fn advanced_spin_drift_matches_canonical_litz_total() {
        let cases = [
            (1.50, 0.25, 850.0, 850.0),
            (1.74, 1.75, 440.0, 800.0),
            (1.90, 1.40, 343.0, 850.0),
        ];

        for (stability, time, velocity, muzzle_velocity) in cases {
            for is_right_twist in [false, true] {
                let expected =
                    crate::spin_drift::litz_drift_meters(stability, time, is_right_twist);
                let actual = calculate_advanced_spin_drift(
                    stability,
                    time,
                    velocity,
                    muzzle_velocity,
                    0.0,
                    0.308 * 0.0254,
                    175.0 * crate::constants::GRAINS_TO_KG,
                    1.225,
                    is_right_twist,
                    "match",
                );

                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "advanced API changed calibrated Litz total: actual={actual} expected={expected}"
                );
            }
        }
    }

    #[test]
    fn litz_total_ignores_legacy_refinement_arguments() {
        let stability = 1.74;
        let time = 1.75;
        let expected = crate::spin_drift::litz_drift_meters(stability, time, true);
        let legacy_states = [
            (440.0, 800.0, 0.0, 0.00782, 0.01134, 1.225, "match"),
            (343.0, 900.0, 19_000.0, 0.00556, 0.00500, 1.000, "vld"),
            (
                100.0,
                700.0,
                -25_000.0,
                0.01270,
                0.04860,
                0.900,
                "flat_base",
            ),
            (
                1_000.0, 1_200.0, 2_000.0, 0.02000, 0.10000, 1.400, "unknown",
            ),
        ];

        for (velocity, muzzle_velocity, spin, caliber, mass, density, bullet_type) in legacy_states
        {
            let actual = calculate_advanced_spin_drift(
                stability,
                time,
                velocity,
                muzzle_velocity,
                spin,
                caliber,
                mass,
                density,
                true,
                bullet_type,
            );
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
    }

    #[test]
    fn advanced_spin_drift_has_no_fixed_jump_intercept() {
        let time = 1e-9;
        let expected = crate::spin_drift::litz_drift_meters(1.74, time, true);
        let actual = calculate_advanced_spin_drift(
            1.74,
            time,
            800.0,
            800.0,
            17_000.0,
            0.308 * 0.0254,
            175.0 * crate::constants::GRAINS_TO_KG,
            1.225,
            true,
            "match",
        );

        assert_eq!(actual.to_bits(), expected.to_bits());
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
        assert_eq!(
            default_coeffs.litz_coefficient,
            match_coeffs.litz_coefficient
        );
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

        assert!(
            drift_medium > drift_short,
            "Drift should increase with time"
        );
        assert!(drift_long > drift_medium, "Drift should increase with time");
    }

    #[test]
    #[allow(deprecated)]
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
    #[allow(deprecated)]
    fn advanced_yaw_of_repose_matches_gravity_gyroscopic_scaling() {
        let calculate = |stability, velocity, crosswind, spin, density, caliber| {
            calculate_advanced_yaw_of_repose(stability, velocity, crosswind, spin, density, caliber)
        };

        let reference = crate::precession_nutation::PrecessionNutationParams::default();
        let expected = crate::precession_nutation::calculate_limit_cycle_yaw_with_inertias(
            850.0,
            17_522.0,
            2.5,
            reference.spin_inertia,
            reference.transverse_inertia,
        );
        let actual = calculate(2.5, 850.0, 0.0, 17_522.0, 1.225, 0.00782);
        assert_eq!(actual.to_bits(), expected.to_bits());

        for crosswind in [-10.0, 10.0] {
            assert_eq!(
                calculate(2.5, 850.0, crosswind, 17_522.0, 1.225, 0.00782).to_bits(),
                actual.to_bits()
            );
        }

        let fixed_sg_fast = calculate(2.5, 800.0, 0.0, 17_522.0, 1.225, 0.00782);
        let fixed_sg_slow = calculate(2.5, 400.0, 0.0, 17_522.0, 1.225, 0.00782);
        assert!((fixed_sg_slow / fixed_sg_fast - 2.0).abs() < 2e-12);

        let physical_fast = calculate(1.5, 800.0, 0.0, 17_522.0, 1.225, 0.00782);
        let physical_slow = calculate(6.0, 400.0, 0.0, 17_522.0, 1.225, 0.00782);
        assert!((physical_slow / physical_fast - 8.0).abs() < 8e-12);

        let dense = calculate(1.5, 300.0, 0.0, 17_522.0, 1.225, 0.00782);
        let thin_same_sg = calculate(1.5, 300.0, 0.0, 17_522.0, 0.6125, 0.00782);
        let thin_coupled_sg = calculate(3.0, 300.0, 0.0, 17_522.0, 0.6125, 0.00782);
        assert_eq!(thin_same_sg.to_bits(), dense.to_bits());
        assert!((thin_coupled_sg / dense - 2.0).abs() < 2e-12);
    }

    #[test]
    #[allow(deprecated)]
    fn test_yaw_of_repose_edge_cases() {
        // Zero stability should give zero yaw
        let zero_stability =
            calculate_advanced_yaw_of_repose(0.5, 800.0, 5.0, 1500.0, 1.225, 0.00782);
        assert_eq!(zero_stability, 0.0);

        // Zero velocity should give zero yaw
        let zero_velocity = calculate_advanced_yaw_of_repose(1.5, 0.0, 5.0, 1500.0, 1.225, 0.00782);
        assert_eq!(zero_velocity, 0.0);

        let zero_spin = calculate_advanced_yaw_of_repose(1.5, 800.0, 5.0, 0.0, 1.225, 0.00782);
        assert_eq!(zero_spin, 0.0);

        let positive_spin =
            calculate_advanced_yaw_of_repose(1.5, 800.0, 5.0, 1500.0, 1.225, 0.00782);
        let negative_spin =
            calculate_advanced_yaw_of_repose(1.5, 800.0, 5.0, -1500.0, 1.225, 0.00782);
        assert_eq!(negative_spin.to_bits(), positive_spin.to_bits());

        let at_boundary = calculate_advanced_yaw_of_repose(1.0, 800.0, 0.0, 1500.0, 1.225, 0.00782);
        let below_boundary = calculate_advanced_yaw_of_repose(
            1.0_f64.next_down(),
            800.0,
            0.0,
            1500.0,
            1.225,
            0.00782,
        );
        assert!(at_boundary > 0.0);
        assert_eq!(below_boundary, 0.0);

        // No crosswind should still give small yaw (trajectory curvature)
        let no_wind = calculate_advanced_yaw_of_repose(1.5, 800.0, 0.0, 1500.0, 1.225, 0.00782);
        assert!(
            no_wind > 0.0,
            "Should have natural yaw from trajectory curvature"
        );
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
        assert!(
            over_stab_subsonic < 0.1,
            "Over-stabilized subsonic should drift less"
        );

        // Long flight subsonic
        let long_subsonic = apply_ml_correction(0.1, 1.5, 0.85, 2.5, 0.308, 168.0);
        assert!(
            long_subsonic > 0.1,
            "Long subsonic flight should need more correction"
        );

        // Light small caliber
        let light_small = apply_ml_correction(0.1, 1.5, 2.5, 1.0, 0.224, 55.0);
        assert!(light_small < 0.1, "Light small caliber should drift less");
    }
}
