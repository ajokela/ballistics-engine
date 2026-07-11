//! Reynolds-number utilities and a legacy correction for genuinely low-Re flow.
//!
//! Standard projectile drag tables already include the Reynolds-number range in which they were
//! measured, so this module leaves ordinary ballistic inputs unchanged. The bounded correction is
//! retained only below `Re=10,000`, outside that standard-table regime.

/// Flow regime classification based on Reynolds number
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
enum FlowRegime {
    Laminar,      // Re < 2000
    Transitional, // 2000 < Re < 5e5
    Turbulent,    // Re > 5e5
}

/// Calculate Reynolds number for a projectile
///
/// Re = ρ × V × L / μ
///
/// # Arguments
/// * `velocity_mps` - Velocity in meters per second
/// * `diameter_m` - Projectile diameter in meters
/// * `air_density_kg_m3` - Air density in kg/m³
/// * `temperature_k` - Temperature in Kelvin
fn calculate_reynolds_number(
    velocity_mps: f64,
    diameter_m: f64,
    air_density_kg_m3: f64,
    temperature_k: f64,
) -> f64 {
    let mu = calculate_air_viscosity(temperature_k);
    air_density_kg_m3 * velocity_mps * diameter_m / mu
}

/// Calculate dynamic viscosity of air using Sutherland's formula
///
/// # Arguments
/// * `temperature_k` - Temperature in Kelvin
///
/// # Returns
/// Dynamic viscosity in Pa·s (kg/m·s)
fn calculate_air_viscosity(temperature_k: f64) -> f64 {
    // Reference values
    const T0: f64 = 273.15; // Reference temperature (K)
    const MU0: f64 = 1.716e-5; // Reference viscosity at T0 (Pa·s)
    const S: f64 = 110.4; // Sutherland's constant for air (K)

    // Sutherland's formula
    MU0 * (T0 + S) / (temperature_k + S) * (temperature_k / T0).powf(1.5)
}

/// Determine flow regime based on Reynolds number
#[allow(dead_code)]
fn get_flow_regime(reynolds_number: f64) -> FlowRegime {
    if reynolds_number < 2000.0 {
        FlowRegime::Laminar
    } else if reynolds_number < 5e5 {
        FlowRegime::Transitional
    } else {
        FlowRegime::Turbulent
    }
}

/// Calculate drag coefficient correction factor based on Reynolds number
///
/// Standard G-model drag is returned unchanged at `Re >= 10,000`. Below that threshold, the
/// legacy low-Re factor is blended back to `1.0` as either Reynolds number or Mach approaches the
/// standard regime.
///
/// # Arguments
/// * `reynolds_number` - Reynolds number
/// * `mach` - Mach number
/// * `base_cd` - Base drag coefficient from standard model
///
/// # Returns
/// Correction factor (multiply by base Cd)
fn reynolds_drag_correction(reynolds_number: f64, mach: f64, _base_cd: f64) -> f64 {
    // Standard projectile drag tables already cover the ordinary ballistic Reynolds-number
    // range. Restrict this opt-in correction to genuinely low-Re flow instead of treating
    // Re=1e6 as the single condition where the tables are valid.
    const FULL_LOW_RE_CORRECTION: f64 = 1e3;
    const STANDARD_TABLE_RE: f64 = 1e4;

    if reynolds_number >= STANDARD_TABLE_RE {
        return 1.0;
    }

    // Only apply corrections for subsonic flow.
    if mach >= 1.0 {
        return 1.0;
    }

    // Preserve the legacy low-Re power-law shape below the standard-table threshold. This is a
    // normalization scale for that residual only, not a claim that G tables require Re=1e6.
    const LOW_RE_POWER_SCALE: f64 = 1e6;

    // Low Reynolds number correction
    let correction = if reynolds_number < 1000.0 {
        // Very low Re - Stokes flow region; Cd increases dramatically
        24.0 / reynolds_number + 1.0
    } else {
        // Low Re - significant viscous effects (power law)
        (reynolds_number / LOW_RE_POWER_SCALE).powf(-0.4)
    };
    // Retain the existing sanity cap for the specialized low-Re path.
    let correction = correction.min(5.0);

    // Fade the specialized low-Re correction into the standard-table regime, avoiding a new
    // discontinuity at Re=1e4. Also fade it through the upper subsonic band so the factor reaches
    // exactly 1.0 at Mach 1 instead of introducing a sonic Cd step.
    const SONIC_FADE_START: f64 = 0.8;
    let reynolds_weight = ((STANDARD_TABLE_RE - reynolds_number)
        / (STANDARD_TABLE_RE - FULL_LOW_RE_CORRECTION))
        .clamp(0.0, 1.0);
    let sonic_weight = ((1.0 - mach) / (1.0 - SONIC_FADE_START)).clamp(0.0, 1.0);

    1.0 + (correction - 1.0) * reynolds_weight * sonic_weight
}

/// Calculate drag coefficient with Reynolds number correction
///
/// # Arguments
/// * `velocity_mps` - Velocity in meters per second
/// * `diameter_m` - Projectile diameter in meters
/// * `air_density_kg_m3` - Air density in kg/m³
/// * `temperature_k` - Temperature in Kelvin
/// * `mach` - Mach number
/// * `base_cd` - Base drag coefficient from standard model
///
/// # Returns
/// Tuple of (corrected_cd, reynolds_number)
fn calculate_corrected_drag(
    velocity_mps: f64,
    diameter_m: f64,
    air_density_kg_m3: f64,
    temperature_k: f64,
    mach: f64,
    base_cd: f64,
) -> (f64, f64) {
    let re = calculate_reynolds_number(velocity_mps, diameter_m, air_density_kg_m3, temperature_k);
    let correction = reynolds_drag_correction(re, mach, base_cd);
    let corrected_cd = base_cd * correction;

    (corrected_cd, re)
}

/// Apply the low-Re correction to a drag coefficient (convenience function).
///
/// Coefficients from standard projectile drag tables are unchanged for ordinary ballistic
/// Reynolds numbers (`Re >= 10,000`) because those empirical tables already include their
/// measured Reynolds-number dependence.
///
/// # Arguments
/// * `base_cd` - Base drag coefficient from G1/G7 model
/// * `velocity_mps` - Velocity in meters per second
/// * `diameter_inches` - Bullet diameter in inches
/// * `air_density_kg_m3` - Air density in kg/m³
/// * `temperature_c` - Temperature in Celsius
/// * `mach` - Mach number
///
/// # Returns
/// Corrected drag coefficient
pub fn apply_reynolds_correction(
    base_cd: f64,
    velocity_mps: f64,
    diameter_inches: f64,
    air_density_kg_m3: f64,
    temperature_c: f64,
    mach: f64,
) -> f64 {
    // Convert units
    let diameter_m = diameter_inches * 0.0254; // inches to meters
    let temperature_k = temperature_c + 273.15; // Celsius to Kelvin

    // Skip correction for very high velocities
    if velocity_mps > 1000.0 || mach > 3.0 {
        return base_cd;
    }

    // Calculate and apply correction
    let (corrected_cd, _) = calculate_corrected_drag(
        velocity_mps,
        diameter_m,
        air_density_kg_m3,
        temperature_k,
        mach,
        base_cd,
    );

    corrected_cd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_air_viscosity() {
        // Test at standard temperature (15°C = 288.15K)
        let mu = calculate_air_viscosity(288.15);
        assert!((mu - 1.789e-5).abs() < 1e-7);
    }

    #[test]
    fn test_reynolds_number() {
        let re = calculate_reynolds_number(100.0, 0.00782, 1.225, 288.15);
        // Should be around 5.4e4
        assert!(re > 5e4 && re < 6e4);
    }

    #[test]
    fn test_flow_regime() {
        assert_eq!(get_flow_regime(1000.0), FlowRegime::Laminar);
        assert_eq!(get_flow_regime(1e5), FlowRegime::Transitional);
        assert_eq!(get_flow_regime(1e6), FlowRegime::Turbulent);
    }

    #[test]
    fn test_reynolds_correction() {
        // Test no correction for supersonic
        let correction = reynolds_drag_correction(1e5, 1.5, 0.5);
        assert_eq!(correction, 1.0);

        // Test correction for low Re
        let correction = reynolds_drag_correction(5e3, 0.5, 0.5);
        assert!(correction > 1.0);

        // Test extreme low Re
        let correction = reynolds_drag_correction(100.0, 0.1, 0.5);
        // Just check that there's significant correction
        assert!(correction > 1.0);
        assert!(correction <= 5.0); // Should be capped
    }

    #[test]
    fn ordinary_rifle_reynolds_range_uses_standard_table_unchanged() {
        let cases = [
            (100.0, 100.0 / 340.0, 53_559.675_375),
            (300.0, 300.0 / 340.0, 160_679.026_126),
        ];

        for (velocity_mps, mach, expected_re) in cases {
            let re = calculate_reynolds_number(velocity_mps, 0.308 * 0.0254, 1.225, 288.15);
            assert!((re - expected_re).abs() < 0.1);

            let base_cd = 0.5;
            let corrected_cd =
                apply_reynolds_correction(base_cd, velocity_mps, 0.308, 1.225, 15.0, mach);

            assert_eq!(
                corrected_cd.to_bits(),
                base_cd.to_bits(),
                "ordinary rifle Re changed drag at Re={re}"
            );
        }
    }

    #[test]
    fn reynolds_correction_is_continuous_at_mach_one() {
        let base_cd = 0.5;
        let below = apply_reynolds_correction(base_cd, 340.0, 0.308, 1.225, 15.0, 1.0 - 1e-6);
        let above = apply_reynolds_correction(base_cd, 340.0, 0.308, 1.225, 15.0, 1.0 + 1e-6);

        assert!(
            (below - above).abs() <= base_cd * 1e-5,
            "sonic Reynolds correction jumped: below={below}, above={above}"
        );
    }

    #[test]
    fn low_re_residual_fades_to_identity_at_mach_one() {
        let below = reynolds_drag_correction(5e3, 1.0 - 1e-8, 0.5);
        let above = reynolds_drag_correction(5e3, 1.0 + 1e-8, 0.5);

        assert!(
            (below - above).abs() <= 1e-6,
            "low-Re correction jumped at Mach 1: below={below}, above={above}"
        );
    }

    #[test]
    fn reynolds_correction_is_continuous_at_standard_table_threshold() {
        let below = reynolds_drag_correction(1e4 - 1e-6, 0.5, 0.5);
        let above = reynolds_drag_correction(1e4 + 1e-6, 0.5, 0.5);

        assert!(
            (below - above).abs() <= 1e-9,
            "Reynolds correction jumped at standard-table threshold: below={below}, above={above}"
        );
    }
}

// Removed Python-specific function
