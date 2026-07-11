use crate::constants::{AIR_DENSITY_SEA_LEVEL, SPEED_OF_SOUND_MPS};

/// Components of aerodynamic jump calculation
#[derive(Debug, Clone, Copy)]
pub struct AerodynamicJumpComponents {
    pub vertical_jump_moa: f64,    // Vertical displacement in MOA at 100 yards
    pub horizontal_jump_moa: f64,  // Horizontal displacement in MOA at 100 yards
    pub jump_angle_rad: f64,       // Total angular displacement in radians
    pub magnus_component_moa: f64, // Magnus effect contribution
    pub yaw_component_moa: f64,    // Initial yaw contribution
    pub stabilization_factor: f64, // How quickly projectile stabilizes (0-1)
}

/// Bryan Litz's crosswind aerodynamic-jump estimator ("Applied Ballistics").
///
/// Linear regression for the VERTICAL jump a crosswind imparts to a
/// spin-stabilized bullet, in MOA per mph of crosswind:
///
/// ```text
///   Y = 0.01*Sg - 0.0024*L + 0.032        [MOA per mph]
/// ```
///
/// where `sg` is the gyroscopic (Miller) stability factor and `length_calibers`
/// is the bullet length in calibers. The returned value is the SIGNED vertical
/// jump in MOA for the supplied crosswind and twist hand: per Litz a right-twist
/// bullet jumps UP for a crosswind from the right and DOWN for one from the left,
/// so `crosswind_from_right_mph` is positive for a wind coming from the right.
///
/// This is a regression valid mainly near `Sg ~ 1.75`; accuracy degrades for
/// bullets well outside the fitted data set. See MBA-959.
pub fn litz_crosswind_jump_moa(
    sg: f64,
    length_calibers: f64,
    crosswind_from_right_mph: f64,
    is_right_twist: bool,
) -> f64 {
    let y_per_mph = 0.01 * sg - 0.0024 * length_calibers + 0.032;
    let hand = if is_right_twist { 1.0 } else { -1.0 };
    hand * y_per_mph * crosswind_from_right_mph
}

/// Legacy heuristic aerodynamic-jump model. NOTE: the trajectory solver uses the
/// validated [`litz_crosswind_jump_moa`] estimator instead; this self-calibrated
/// model (with a hand-tuned `magnus_enhancement` factor and a known horizontal-sign
/// quirk) is retained only for backward compatibility with external callers.
///
/// Calculate aerodynamic jump for a spinning projectile.
///
/// Aerodynamic jump is the displacement of the projectile's trajectory
/// as it transitions from constrained motion in the barrel to free flight.
pub fn calculate_aerodynamic_jump(
    muzzle_velocity_mps: f64,
    spin_rate_rad_s: f64,
    crosswind_mps: f64,
    caliber_m: f64,
    mass_kg: f64,
    barrel_length_m: f64,
    twist_rate_calibers: f64,
    is_right_twist: bool,
    initial_yaw_rad: f64,
    air_density_kg_m3: f64,
) -> AerodynamicJumpComponents {
    if muzzle_velocity_mps <= 0.0
        || caliber_m <= 0.0
        || mass_kg <= 0.0
        || twist_rate_calibers <= 0.0
    {
        return AerodynamicJumpComponents {
            vertical_jump_moa: 0.0,
            horizontal_jump_moa: 0.0,
            jump_angle_rad: 0.0,
            magnus_component_moa: 0.0,
            yaw_component_moa: 0.0,
            stabilization_factor: 0.0,
        };
    }

    // Calculate Magnus force coefficient
    let mach = muzzle_velocity_mps / SPEED_OF_SOUND_MPS;
    let magnus_coeff = if mach < 0.8 {
        0.25
    } else if mach < 1.2 {
        0.15 // Reduced in transonic
    } else {
        0.20
    };

    // Spin parameter (non-dimensional)
    let spin_param = (spin_rate_rad_s * caliber_m / 2.0) / muzzle_velocity_mps;

    // Effective yaw angle during muzzle exit
    let crosswind_yaw = if crosswind_mps != 0.0 {
        (crosswind_mps / muzzle_velocity_mps).atan()
    } else {
        0.0
    };

    let total_yaw_rad = crosswind_yaw + initial_yaw_rad;

    // Magnus force during barrel exit
    let area = std::f64::consts::PI * (caliber_m / 2.0).powi(2);
    let magnus_force = 0.5
        * air_density_kg_m3
        * muzzle_velocity_mps.powi(2)
        * area
        * magnus_coeff
        * spin_param
        * total_yaw_rad.sin();

    // Time for projectile to clear muzzle
    let exit_time = 2.0 * barrel_length_m / muzzle_velocity_mps;

    // Stabilization distance
    let stabilization_calibers = 20.0 / (twist_rate_calibers / 10.0).sqrt();
    let stabilization_distance = stabilization_calibers * caliber_m;
    let stabilization_time = stabilization_distance / muzzle_velocity_mps;

    // Total effective time
    let effective_time = exit_time + stabilization_time;

    // Calculate jump displacement. Direction comes from the crosswind, falling back to the yaw
    // direction when there is no crosswind — signum(0.0) == +1.0 would otherwise impose a phantom
    // positive direction for pure-yaw (no-wind) inputs.
    let dir_sign = if crosswind_mps != 0.0 {
        crosswind_mps.signum()
    } else {
        total_yaw_rad.signum()
    };
    let vertical_sign = if is_right_twist { dir_sign } else { -dir_sign };

    // Magnus acceleration
    let magnus_accel = magnus_force / mass_kg;

    // Enhanced calculation accounting for barrel exit dynamics
    let lever_factor = (barrel_length_m / caliber_m) * 0.1;
    let magnus_enhancement = 50.0; // Calibrated to match empirical data

    // Vertical displacement
    let mut vertical_jump_m = magnus_enhancement
        * lever_factor
        * vertical_sign
        * magnus_accel.abs()
        * effective_time.powi(2);

    // Add yaw-induced component
    if total_yaw_rad != 0.0 {
        let yaw_contribution = total_yaw_rad.abs() * barrel_length_m * 0.5;
        vertical_jump_m += vertical_sign * yaw_contribution;
    }

    // Horizontal component (smaller effect)
    let horizontal_jump_m = 0.25 * vertical_jump_m * (2.0 * total_yaw_rad).sin();

    // Convert to MOA at 100 yards
    const YARDS_TO_M: f64 = 0.9144;
    const MOA_PER_RADIAN: f64 = 3437.7467707849; // 1 / 0.0002908882

    let range_100y = 100.0 * YARDS_TO_M;
    let vertical_angle_rad = vertical_jump_m / range_100y;
    let horizontal_angle_rad = horizontal_jump_m / range_100y;

    let vertical_jump_moa = vertical_angle_rad * MOA_PER_RADIAN;
    let horizontal_jump_moa = horizontal_angle_rad * MOA_PER_RADIAN;

    // Total jump angle
    let total_jump_rad = (vertical_angle_rad.powi(2) + horizontal_angle_rad.powi(2)).sqrt();

    // Component breakdown
    let magnus_component_moa = vertical_jump_moa.abs() * 0.8;
    let yaw_component_moa = vertical_jump_moa.abs() * 0.2;

    // Stabilization factor
    let caliber_in = caliber_m / 0.0254;
    // mass_kg -> grains uses the kilograms->grains factor (15432.358), NOT the grams->grains
    // factor (15.432); the latter made sg_approx ~1000x too small so stabilization_factor
    // collapsed to ~0 for every real projectile. (Simplified Sg: length terms omitted.)
    let sg_approx =
        30.0 * mass_kg * 15432.358 / (twist_rate_calibers.powi(2) * caliber_in.powi(3));
    let stabilization_factor = (sg_approx / 1.5).min(1.0);

    AerodynamicJumpComponents {
        vertical_jump_moa,
        horizontal_jump_moa,
        jump_angle_rad: total_jump_rad,
        magnus_component_moa,
        yaw_component_moa,
        stabilization_factor,
    }
}

/// Calculate sight corrections needed to compensate for aerodynamic jump
///
/// Aerodynamic jump is an angular muzzle departure, so its equal-and-opposite
/// sight correction is independent of zero range and sight height. Those
/// parameters remain in the public signature for backward compatibility.
pub fn calculate_sight_correction_for_jump(
    jump_components: &AerodynamicJumpComponents,
    zero_range_m: f64,
    _sight_height_m: f64,
) -> (f64, f64) {
    // Preserve the public helper's established invalid-range behavior.
    if !zero_range_m.is_finite() || zero_range_m <= 0.0 {
        return (0.0, 0.0);
    }

    (
        -jump_components.vertical_jump_moa,
        -jump_components.horizontal_jump_moa,
    )
}

/// Calculate sensitivity to crosswind for aerodynamic jump (MOA per mph)
pub fn calculate_crosswind_jump_sensitivity(
    muzzle_velocity_mps: f64,
    spin_rate_rad_s: f64,
    caliber_m: f64,
    mass_kg: f64,
    twist_rate_calibers: f64,
    is_right_twist: bool,
) -> f64 {
    const MPH_TO_MPS: f64 = 0.44704;
    let crosswind_1mph = MPH_TO_MPS;

    let jump = calculate_aerodynamic_jump(
        muzzle_velocity_mps,
        spin_rate_rad_s,
        crosswind_1mph,
        caliber_m,
        mass_kg,
        0.6, // Typical 24" barrel
        twist_rate_calibers,
        is_right_twist,
        0.0, // No initial yaw
        AIR_DENSITY_SEA_LEVEL,
    );

    jump.vertical_jump_moa.abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aerodynamic_jump_zero_conditions() {
        // Test with no crosswind
        let jump = calculate_aerodynamic_jump(
            800.0,   // velocity
            1000.0,  // spin rate
            0.0,     // no crosswind
            0.00762, // .30 cal
            0.01134, // 175gr
            0.6,     // barrel length
            32.47,   // twist rate in calibers
            true,    // right twist
            0.0,     // no initial yaw
            1.225,   // air density
        );

        assert_eq!(jump.vertical_jump_moa, 0.0);
        assert!(jump.horizontal_jump_moa.abs() < 0.001);
    }

    #[test]
    fn test_aerodynamic_jump_with_crosswind() {
        // Test with 10 mph right crosswind
        let jump = calculate_aerodynamic_jump(
            800.0,   // velocity
            17593.0, // spin rate for 1:10 twist
            4.4704,  // 10 mph crosswind
            0.00782, // .308 cal
            0.01134, // 175gr
            0.6096,  // 24" barrel
            32.47,   // twist rate in calibers
            true,    // right twist
            0.0,     // no initial yaw
            1.225,   // air density
        );

        // Right twist + right wind should give positive (upward) jump
        assert!(jump.vertical_jump_moa > 0.0);
        // Just check that we have some stabilization
        assert!(jump.stabilization_factor > 0.0);
    }

    #[test]
    fn test_opposite_twist_direction() {
        let crosswind = 4.4704; // 10 mph

        // Right twist
        let jump_right = calculate_aerodynamic_jump(
            800.0, 17593.0, crosswind, 0.00782, 0.01134, 0.6096, 32.47, true, 0.0, 1.225,
        );

        // Left twist
        let jump_left = calculate_aerodynamic_jump(
            800.0, 17593.0, crosswind, 0.00782, 0.01134, 0.6096, 32.47, false, 0.0, 1.225,
        );

        // Opposite twist should give opposite vertical jump
        assert!((jump_right.vertical_jump_moa + jump_left.vertical_jump_moa).abs() < 0.001);
    }

    #[test]
    fn sight_correction_is_the_equal_and_opposite_jump_angle() {
        let jump = AerodynamicJumpComponents {
            vertical_jump_moa: 0.5,
            horizontal_jump_moa: 0.1,
            jump_angle_rad: 0.0001,
            magnus_component_moa: 0.4,
            yaw_component_moa: 0.1,
            stabilization_factor: 0.9,
        };

        for (zero_range_m, sight_height_m) in [
            (22.86, 0.0),  // 25 yards
            (91.44, 0.05), // 100 yards, 2-inch sight height
            (274.32, 0.1), // 300 yards, tall sight
        ] {
            let (vertical, horizontal) =
                calculate_sight_correction_for_jump(&jump, zero_range_m, sight_height_m);

            assert!(
                (vertical + 0.5).abs() < 1e-12,
                "vertical correction changed with range/height: range={zero_range_m}, height={sight_height_m}, correction={vertical}"
            );
            assert!(
                (horizontal + 0.1).abs() < 1e-12,
                "horizontal correction changed with range/height: range={zero_range_m}, height={sight_height_m}, correction={horizontal}"
            );
        }
    }

    #[test]
    fn sight_correction_rejects_invalid_zero_ranges() {
        let jump = AerodynamicJumpComponents {
            vertical_jump_moa: 0.5,
            horizontal_jump_moa: 0.1,
            jump_angle_rad: 0.0001,
            magnus_component_moa: 0.4,
            yaw_component_moa: 0.1,
            stabilization_factor: 0.9,
        };

        for zero_range_m in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                calculate_sight_correction_for_jump(&jump, zero_range_m, 0.05),
                (0.0, 0.0),
                "invalid zero range must be rejected: {zero_range_m}"
            );
        }
    }

    #[test]
    fn test_crosswind_sensitivity() {
        let sensitivity = calculate_crosswind_jump_sensitivity(
            800.0,   // velocity
            17593.0, // spin rate
            0.00782, // caliber
            0.01134, // mass
            32.47,   // twist rate
            true,    // right twist
        );

        // Should be positive and reasonable (typically 0.01-0.1 MOA/mph)
        assert!(sensitivity > 0.0);
        assert!(sensitivity < 0.5);
    }

    // ---- Litz crosswind aerodynamic-jump estimator (the canonical solver model) ----

    #[test]
    fn litz_matches_the_published_formula() {
        // Y = 0.01*Sg - 0.0024*L + 0.032 [MOA/mph], scaled by crosswind and twist hand.
        // Sg = 1.75, L = 4.0 -> 0.0175 - 0.0096 + 0.032 = 0.0399 MOA/mph.
        let per_mph = 0.01 * 1.75 - 0.0024 * 4.0 + 0.032;
        let got = litz_crosswind_jump_moa(1.75, 4.0, 10.0, true);
        assert!(
            (got - per_mph * 10.0).abs() < 1e-12,
            "got {got}, expected {}",
            per_mph * 10.0
        );
        // Sanity: a few tenths of an MOA at 10 mph.
        assert!((got - 0.399).abs() < 1e-3);
    }

    #[test]
    fn litz_is_linear_in_crosswind() {
        let one = litz_crosswind_jump_moa(1.8, 3.5, 1.0, true);
        let ten = litz_crosswind_jump_moa(1.8, 3.5, 10.0, true);
        assert!((ten - 10.0 * one).abs() < 1e-12);
        assert_eq!(litz_crosswind_jump_moa(1.8, 3.5, 0.0, true), 0.0);
    }

    #[test]
    fn litz_sign_flips_with_wind_side_and_twist() {
        // Wind from the right + right twist -> up (positive).
        let base = litz_crosswind_jump_moa(1.9, 4.0, 12.0, true);
        assert!(base > 0.0);
        // Reversing the wind side flips the sign, same magnitude.
        assert!((litz_crosswind_jump_moa(1.9, 4.0, -12.0, true) + base).abs() < 1e-12);
        // Flipping the twist hand flips the sign.
        assert!((litz_crosswind_jump_moa(1.9, 4.0, 12.0, false) + base).abs() < 1e-12);
    }

    #[test]
    fn litz_regression_can_go_negative_outside_its_fitted_range() {
        // The estimator is a faithful linear fit (not clamped): a very long, marginally
        // stable bullet drives 0.01*Sg - 0.0024*L + 0.032 below zero, reversing the jump.
        // This is the extrapolation regime — see MBA-959.
        let per_mph = 0.01 * 1.0 - 0.0024 * 20.0 + 0.032; // = -0.006
        assert!(per_mph < 0.0);
        let got = litz_crosswind_jump_moa(1.0, 20.0, 10.0, true);
        assert!((got - per_mph * 10.0).abs() < 1e-12);
        assert!(got < 0.0);
    }
}
