//! Pitch Damping Moment Physics for Ballistics Calculations
//!
//! This module implements pitch damping moments that affect:
//! - Dynamic stability during flight
//! - Precession and nutation damping  
//! - Yaw of repose convergence
//! - Transonic stability transitions

use crate::constants::G_ACCEL_MPS2;
use crate::spin_decay::calculate_moment_of_inertia;
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
        // These values use the q * S * d * (pitch_rate * d / V) convention in
        // `calculate_pitch_damping_moment`. Published small-arms coefficients are order 1-10;
        // a representative 7.62 mm pitch-damping coefficient is about -4.7.
        Self {
            subsonic: -8.0,
            transonic_low: -3.0,
            transonic_high: 2.0,
            supersonic: -5.0,
        }
    }
}

impl PitchDampingCoefficients {
    /// Get typical coefficients for different bullet types
    pub fn from_bullet_type(bullet_type: &str) -> Self {
        match bullet_type.to_lowercase().as_str() {
            "match_boat_tail" => Self {
                subsonic: -9.0,
                transonic_low: -4.0,
                transonic_high: 1.0,
                supersonic: -6.0,
            },
            "match_flat_base" => Self {
                subsonic: -7.0,
                transonic_low: -2.0,
                transonic_high: 3.0,
                supersonic: -4.0,
            },
            "vld" => Self {
                // Very Low Drag - more stable
                subsonic: -10.0,
                transonic_low: -5.0,
                transonic_high: -1.0,
                supersonic: -7.0,
            },
            "hunting" => Self {
                subsonic: -6.0,
                transonic_low: -1.0,
                transonic_high: 4.0,
                supersonic: -3.0,
            },
            "fmj" => Self {
                subsonic: -7.0,
                transonic_low: -2.0,
                transonic_high: 2.0,
                supersonic: -5.0,
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
    _length_m: f64,
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

/// First-order yaw-of-repose magnitude for a flat-fire trajectory.
///
/// Following AMCP 706-238 section 4-11.2, the classical relations
/// `Sg = (Ix * p)^2 / (4 * Iy * M)` and
/// `yaw = Ix * p * g / (M * V)` eliminate the unavailable static-moment slope `M`, leaving
/// `yaw = 4 * Iy * Sg * g / (Ix * |p| * V)`. The API has no flight-path angle, so this uses
/// the flat-fire approximation `cos(theta) = 1`. Twist direction is applied downstream.
pub(crate) fn calculate_gravity_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    mass_kg: f64,
    caliber_m: f64,
    length_m: f64,
) -> f64 {
    if !stability_factor.is_finite()
        || stability_factor <= 1.0
        || !velocity_mps.is_finite()
        || velocity_mps <= 0.0
        || !spin_rate_rad_s.is_finite()
        || spin_rate_rad_s == 0.0
        || !mass_kg.is_finite()
        || mass_kg <= 0.0
        || !caliber_m.is_finite()
        || caliber_m <= 0.0
        || !length_m.is_finite()
        || length_m <= 0.0
    {
        return 0.0;
    }

    let axial_inertia = calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");
    let transverse_inertia =
        calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");
    if axial_inertia <= 0.0 || transverse_inertia <= 0.0 {
        return 0.0;
    }

    4.0 * transverse_inertia * stability_factor * G_ACCEL_MPS2
        / (axial_inertia * spin_rate_rad_s.abs() * velocity_mps)
}

/// Calculate yaw of repose with pitch damping effects.
///
/// Returns the equilibrium yaw and a signed convergence rate in `s^-1`: positive values
/// converge toward equilibrium, while negative values identify a divergent pitch mode.
#[allow(clippy::too_many_arguments)] // Public compatibility API; grouping would be breaking.
pub fn calculate_damped_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    _wind_velocity_mps: f64,
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
    let mass_kg = mass_grains * crate::constants::GRAINS_TO_KG;

    // Crosswind creates an initial transient handled by aerodynamic-jump physics; it is not part
    // of the persistent equilibrium yaw. Use the gravity/gyroscopic balance for repose instead.
    let equilibrium_yaw_rad = calculate_gravity_yaw_of_repose(
        stability_factor,
        velocity_mps,
        spin_rate_rad_s,
        mass_kg,
        caliber_m,
        length_m,
    );

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
    let i_transverse =
        calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");

    // Angular acceleration from damping
    let angular_accel = calculate_angular_acceleration(damping_moment, i_transverse);

    // For q_dot = lambda*q, the convergence rate is -lambda: positive for damping and
    // negative for a destabilizing moment. Preserve the legacy zero-signal fallback.
    let convergence_rate = if angular_accel != 0.0 && pitch_rate_rad_s != 0.0 {
        -angular_accel / pitch_rate_rad_s
    } else {
        0.1
    };

    (equilibrium_yaw_rad, convergence_rate)
}

/// Legacy alias for the slow-mode precession angular frequency in radians per second.
///
/// Pitch damping changes modal amplitude/convergence, not the slow-mode phase frequency. The
/// removed yaw, velocity, and damping-moment arguments could not form a dimensionally valid
/// correction. Use [`crate::precession_nutation::calculate_precession_frequency`] directly.
#[deprecated(
    since = "0.22.18",
    note = "use precession_nutation::calculate_precession_frequency; damping changes modal amplitude, not phase frequency"
)]
pub fn calculate_precession_with_damping(
    spin_rate_rad_s: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
    stability_factor: f64,
) -> f64 {
    crate::precession_nutation::calculate_precession_frequency(
        spin_rate_rad_s,
        spin_inertia,
        transverse_inertia,
        stability_factor,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pitch_damping_coefficient() {
        let coeffs = PitchDampingCoefficients::default();

        // Subsonic
        assert_eq!(calculate_pitch_damping_coefficient(0.5, &coeffs), -8.0);

        // Transonic
        let transonic = calculate_pitch_damping_coefficient(0.9, &coeffs);
        assert!(transonic > -8.0 && transonic < -3.0);

        // Supersonic
        let supersonic = calculate_pitch_damping_coefficient(2.0, &coeffs);
        assert_eq!(supersonic, -5.0);
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

    #[test]
    fn default_308_pitch_rate_damps_on_small_arms_timescale() {
        let pitch_rate = 0.1;
        let velocity = 850.0;
        let density = 1.225;
        let caliber = 0.308 * 0.0254;
        let length = 1.3 * 0.0254;
        let mass = 175.0 * crate::constants::GRAINS_TO_KG;
        let mach = velocity / 343.0;
        let coeffs = PitchDampingCoefficients::default();

        let moment = calculate_pitch_damping_moment(
            pitch_rate, velocity, density, caliber, length, mach, &coeffs,
        );
        let inertia = calculate_transverse_moment_of_inertia(mass, caliber, length, "ogive");
        let angular_accel = calculate_angular_acceleration(moment, inertia);
        let time_constant = (pitch_rate / angular_accel).abs();

        assert!(
            (0.05..=0.25).contains(&time_constant),
            ".308 pitch damping should settle on a small-arms timescale, got tau={time_constant}s"
        );
    }

    #[test]
    fn test_bullet_type_coefficients() {
        let types = [
            "match_boat_tail",
            "match_flat_base",
            "vld",
            "hunting",
            "fmj",
            "unknown",
        ];

        for bullet_type in &types {
            let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);

            // Check that subsonic is always stabilizing (negative)
            assert!(coeffs.subsonic < 0.0);

            // Check that supersonic eventually stabilizes
            assert!(coeffs.supersonic < 0.0);

            // Stable-regime presets must stay on the published small-arms coefficient scale.
            assert!(coeffs.subsonic.abs() >= 3.0);
            assert!(coeffs.supersonic.abs() >= 3.0);

            // Check that VLD is most stable
            if *bullet_type == "vld" {
                let default_coeffs = PitchDampingCoefficients::default();
                assert!(coeffs.subsonic < default_coeffs.subsonic);
            }
        }
    }

    #[test]
    fn test_transonic_instability() {
        let coeffs = PitchDampingCoefficients::from_bullet_type("hunting");

        // Check that transonic high can be destabilizing (positive)
        assert!(coeffs.transonic_high > 0.0);

        // Check coefficient through transonic region
        let mach_1_1 = calculate_pitch_damping_coefficient(1.1, &coeffs);

        // Should be transitioning toward destabilizing
        assert!(mach_1_1 > coeffs.transonic_low);
    }

    #[test]
    fn test_transverse_moment_of_inertia() {
        let mass_kg = 0.01134; // 175 grains
        let caliber_m = 0.00782; // .308"
        let length_m = 0.033; // 1.3"

        let i_cylinder =
            calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "cylinder");
        let i_ogive = calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");
        let i_boat_tail =
            calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "boat_tail");
        let i_unknown =
            calculate_transverse_moment_of_inertia(mass_kg, caliber_m, length_m, "unknown");

        // Check relative magnitudes
        assert!(i_cylinder > i_ogive);
        assert!(i_ogive > i_boat_tail);
        assert_eq!(i_cylinder, i_unknown);

        // Check absolute values are reasonable
        assert!(i_cylinder > 0.0);
        assert!(i_cylinder < 1.0); // Should be small for a bullet
    }

    #[test]
    fn test_angular_acceleration() {
        let moment = -0.001; // Small damping moment
        let inertia = 0.0001; // Small inertia

        let accel = calculate_angular_acceleration(moment, inertia);
        assert_eq!(accel, moment / inertia);

        // Test zero inertia
        let accel_zero = calculate_angular_acceleration(moment, 0.0);
        assert_eq!(accel_zero, 0.0);
    }

    #[test]
    fn test_damped_yaw_of_repose() {
        let (yaw, convergence) = calculate_damped_yaw_of_repose(
            2.5,     // stability factor
            800.0,   // velocity m/s
            19000.0, // spin rate rad/s
            10.0,    // wind velocity m/s
            0.01,    // pitch rate rad/s
            1.225,   // air density
            0.308,   // caliber inches
            1.3,     // length inches
            175.0,   // mass grains
            0.9,     // Mach
            "match_boat_tail",
        );

        // Should have non-zero yaw and convergence
        assert!(yaw > 0.0);
        assert!(yaw < 0.1); // Should be small angle
        assert!(convergence > 0.0);

        // Test with no stability (Sg <= 1)
        let (yaw_unstable, conv_unstable) = calculate_damped_yaw_of_repose(
            0.9,
            800.0,
            19000.0,
            10.0,
            0.01,
            1.225,
            0.308,
            1.3,
            175.0,
            0.9,
            "match_boat_tail",
        );
        assert_eq!(yaw_unstable, 0.0);
        assert_eq!(conv_unstable, 0.0);
    }

    #[test]
    fn damped_yaw_convergence_rate_preserves_stability_sign() {
        let rate = |mach, pitch_rate_rad_s| {
            calculate_damped_yaw_of_repose(
                2.5,
                800.0,
                19_000.0,
                0.0,
                pitch_rate_rad_s,
                1.225,
                0.308,
                1.3,
                175.0,
                mach,
                "fmj",
            )
            .1
        };

        // FMJ uses equal-and-opposite Cmq values at these regime boundaries.
        let damped = rate(1.0, 0.01);
        let divergent = rate(1.2, 0.01);

        assert!(damped > 0.0);
        assert!(divergent < 0.0);
        assert_eq!(damped.to_bits(), (-divergent).to_bits());
    }

    #[test]
    fn crosswind_is_not_persistent_equilibrium_yaw() {
        let calculate = |wind_velocity_mps| {
            calculate_damped_yaw_of_repose(
                2.5,
                300.0,
                19_000.0,
                wind_velocity_mps,
                0.01,
                1.225,
                0.308,
                1.3,
                175.0,
                0.875,
                "match_boat_tail",
            )
        };

        let (calm_yaw, calm_rate) = calculate(0.0);
        let (windy_yaw, windy_rate) = calculate(10.0);

        assert!(
            (windy_yaw - calm_yaw).abs() < 1e-12,
            "crosswind became persistent equilibrium yaw: calm={calm_yaw} windy={windy_yaw}"
        );
        assert_eq!(windy_rate.to_bits(), calm_rate.to_bits());
        assert!(windy_yaw.abs() < 0.003);
    }

    #[test]
    fn gravity_yaw_of_repose_matches_classical_stability_reduction() {
        let stability_factor = 2.5;
        let velocity_mps = 300.0;
        let spin_rate_rad_s = 19_000.0;
        let mass_kg = 175.0 * crate::constants::GRAINS_TO_KG;
        let caliber_m = 0.308 * 0.0254;
        let length_m = 1.3 * 0.0254;

        let actual = calculate_gravity_yaw_of_repose(
            stability_factor,
            velocity_mps,
            spin_rate_rad_s,
            mass_kg,
            caliber_m,
            length_m,
        );

        // Independent expansion of the inertia approximations and the classical flat-fire
        // reduction: yaw = 4 * Iy * Sg * g / (Ix * |p| * V).
        let radius_m = caliber_m / 2.0;
        let axial_inertia = 0.4 * mass_kg * radius_m.powi(2);
        let transverse_inertia =
            0.85 * mass_kg * (3.0 * radius_m.powi(2) + length_m.powi(2)) / 12.0;
        let expected = 4.0 * transverse_inertia * stability_factor * 9.80665
            / (axial_inertia * spin_rate_rad_s * velocity_mps);

        assert!((actual - expected).abs() < 1e-15);
        assert!((actual - 0.000226244442784).abs() < 1e-15);
    }

    #[test]
    #[allow(deprecated)]
    fn test_precession_with_damping() {
        let precession = calculate_precession_with_damping(
            19000.0, // spin rate rad/s
            0.00005, // spin inertia
            0.0001,  // transverse inertia
            2.5,     // stability factor
        );

        assert!(precession > 0.0);

        // Test zero spin
        let precession_zero = calculate_precession_with_damping(0.0, 0.00005, 0.0001, 2.5);
        assert_eq!(precession_zero, 0.0);

        // Test a non-gyroscopically-stable projectile
        let precession_unstable = calculate_precession_with_damping(19000.0, 0.00005, 0.0001, 1.0);
        assert_eq!(precession_unstable, 0.0);
    }

    #[test]
    #[allow(deprecated)]
    fn precession_uses_slow_epicyclic_frequency() {
        let spin_rate_rad_s = 17_522.0;
        let spin_inertia = 6.94e-8;
        let transverse_inertia = 9.13e-7;
        let stability_factor = 2.0;
        let expected = crate::precession_nutation::calculate_precession_frequency(
            spin_rate_rad_s,
            spin_inertia,
            transverse_inertia,
            stability_factor,
        );

        let actual = calculate_precession_with_damping(
            spin_rate_rad_s,
            spin_inertia,
            transverse_inertia,
            stability_factor,
        );

        assert!(
            (actual - expected).abs() <= expected * 1e-12,
            "precession did not use the slow epicyclic rate: actual={actual} expected={expected}"
        );
        assert!((150.0..250.0).contains(&actual));
    }

    #[test]
    fn test_mach_interpolation() {
        let coeffs = PitchDampingCoefficients::default();

        // Test continuity at every piecewise interpolation boundary without assuming a scale.
        let epsilon = 1e-9;
        for boundary in [0.8, 1.0, 1.2, 2.0] {
            let at_boundary = calculate_pitch_damping_coefficient(boundary, &coeffs);
            let below = calculate_pitch_damping_coefficient(boundary - epsilon, &coeffs);
            let above = calculate_pitch_damping_coefficient(boundary + epsilon, &coeffs);

            assert!((at_boundary - below).abs() < 1e-6);
            assert!((at_boundary - above).abs() < 1e-6);
        }
    }

    #[test]
    fn test_pitch_damping_edge_cases() {
        let coeffs = PitchDampingCoefficients::default();

        // Test zero pitch rate
        let moment_zero_pitch =
            calculate_pitch_damping_moment(0.0, 300.0, 1.225, 0.00782, 0.033, 0.87, &coeffs);
        assert_eq!(moment_zero_pitch, 0.0);

        // Test zero velocity
        let moment_zero_vel =
            calculate_pitch_damping_moment(0.1, 0.0, 1.225, 0.00782, 0.033, 0.87, &coeffs);
        assert_eq!(moment_zero_vel, 0.0);
    }

    #[test]
    fn test_default_implementation() {
        let coeffs1 = PitchDampingCoefficients::default();
        let coeffs2 = PitchDampingCoefficients::from_bullet_type("unknown");

        assert_eq!(coeffs1.subsonic, coeffs2.subsonic);
        assert_eq!(coeffs1.transonic_low, coeffs2.transonic_low);
        assert_eq!(coeffs1.transonic_high, coeffs2.transonic_high);
        assert_eq!(coeffs1.supersonic, coeffs2.supersonic);
    }

    #[test]
    fn test_transonic_jump() {
        let _coeffs = PitchDampingCoefficients::from_bullet_type("hunting");

        // In transonic region, check for potential instability
        let (yaw_subsonic, _) = calculate_damped_yaw_of_repose(
            2.5, 250.0, 19000.0, 10.0, 0.01, 1.225, 0.308, 1.3, 175.0, 0.7, "hunting",
        );

        let (yaw_transonic, _) = calculate_damped_yaw_of_repose(
            2.5, 343.0, 19000.0, 10.0, 0.01, 1.225, 0.308, 1.3, 175.0, 1.0, "hunting",
        );

        // Both should be valid but potentially different
        assert!(yaw_subsonic > 0.0);
        assert!(yaw_transonic > 0.0);
    }
}
