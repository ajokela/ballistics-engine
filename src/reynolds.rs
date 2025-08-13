//! Reynolds number corrections for low-velocity drag.
//!
//! At low velocities, viscous effects become more significant and the drag coefficient
//! increases. This module provides corrections based on Reynolds number to improve
//! accuracy for subsonic projectiles and end-of-trajectory calculations.

use pyo3::prelude::*;

/// Flow regime classification based on Reynolds number
#[derive(Debug, Clone, Copy, PartialEq)]
enum FlowRegime {
    Laminar,        // Re < 2000
    Transitional,   // 2000 < Re < 5e5
    Turbulent,      // Re > 5e5
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
    const T0: f64 = 273.15;      // Reference temperature (K)
    const MU0: f64 = 1.716e-5;   // Reference viscosity at T0 (Pa·s)
    const S: f64 = 110.4;        // Sutherland's constant for air (K)
    
    // Sutherland's formula
    MU0 * (T0 + S) / (temperature_k + S) * (temperature_k / T0).powf(1.5)
}

/// Determine flow regime based on Reynolds number
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
/// This correction is most significant at low velocities where viscous effects
/// dominate. The correction factor is multiplied by the base drag coefficient.
/// 
/// # Arguments
/// * `reynolds_number` - Reynolds number
/// * `mach` - Mach number
/// * `base_cd` - Base drag coefficient from standard model
/// 
/// # Returns
/// Correction factor (multiply by base Cd)
fn reynolds_drag_correction(
    reynolds_number: f64,
    mach: f64,
    _base_cd: f64,
) -> f64 {
    // Only apply corrections for subsonic flow
    if mach >= 1.0 {
        return 1.0;
    }
    
    // Reference Reynolds number where standard drag curves are valid
    const RE_REF: f64 = 1e6;
    
    if reynolds_number >= RE_REF {
        return 1.0;
    }
    
    // Low Reynolds number correction
    if reynolds_number < 1000.0 {
        // Very low Re - Stokes flow region
        // Cd increases dramatically
        let correction = 24.0 / reynolds_number + 1.0;
        // Limit correction to reasonable values
        correction.min(5.0)
    } else if reynolds_number < 1e4 {
        // Low Re - significant viscous effects
        // Use power law interpolation
        let n = -0.4;
        (reynolds_number / RE_REF).powf(n)
    } else if reynolds_number < 1e5 {
        // Moderate Re - transitional region
        let n = -0.2;
        (reynolds_number / RE_REF).powf(n)
    } else {
        // High Re but below reference
        let n = -0.1;
        (reynolds_number / RE_REF).powf(n)
    }
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

/// Apply Reynolds number correction to drag coefficient (convenience function)
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
    let diameter_m = diameter_inches * 0.0254;  // inches to meters
    let temperature_k = temperature_c + 273.15;  // Celsius to Kelvin
    
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
        let correction = reynolds_drag_correction(1e4, 0.5, 0.5);
        assert!(correction > 1.0);
        
        // Test extreme low Re
        let correction = reynolds_drag_correction(100.0, 0.1, 0.5);
        assert!(correction > 2.0);
        assert!(correction <= 5.0); // Should be capped
    }
}

/// Python-exposed function for applying Reynolds correction
#[pyfunction]
#[pyo3(name = "apply_reynolds_correction_rust")]
pub fn apply_reynolds_correction_py(
    base_cd: f64,
    velocity_mps: f64,
    diameter_inches: f64,
    air_density_kg_m3: f64,
    temperature_c: f64,
    mach: f64,
) -> PyResult<f64> {
    Ok(apply_reynolds_correction(
        base_cd,
        velocity_mps,
        diameter_inches,
        air_density_kg_m3,
        temperature_c,
        mach,
    ))
}