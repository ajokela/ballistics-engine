/// Transonic drag modeling with shock wave effects
/// 
/// This module implements physics-based corrections for drag in the transonic regime
/// (Mach 0.8-1.2) where shock waves significantly affect the drag coefficient.

use pyo3::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectileShape {
    Spitzer,    // Sharp pointed (most common)
    RoundNose,  // Blunt/round nose
    FlatBase,   // Flat base (wadcutter)
    BoatTail,   // Boat tail design
}

impl ProjectileShape {
    /// Parse from string representation
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "spitzer" => Self::Spitzer,
            "round_nose" => Self::RoundNose,
            "flat_base" => Self::FlatBase,
            "boat_tail" => Self::BoatTail,
            _ => Self::Spitzer, // Default
        }
    }
}

/// Calculate Prandtl-Glauert correction factor for compressibility effects
/// 
/// This factor accounts for air compressibility effects in subsonic flow approaching Mach 1.
/// Formula: 1/sqrt(1-M²) where M is Mach number
/// 
/// Physical basis: As flow approaches sonic speeds, local acceleration over the projectile
/// surface causes compressibility effects that increase the effective angle of attack and drag.
/// 
/// Note: This correction is only valid for subsonic flow (Mach < 1.0)
fn prandtl_glauert_correction(mach: f64) -> f64 {
    if mach >= 0.99 {
        // Near Mach 1, the theoretical correction approaches infinity (1/sqrt(1-1²) = 1/0)
        // Cap at 10.0 to prevent numerical divergence while maintaining physical meaning
        // This represents a practical limit where the linear theory breaks down
        return 10.0;
    }
    
    // Classic Prandtl-Glauert compressibility correction factor
    // Accounts for increased pressure coefficients due to compressibility
    let beta = (1.0 - mach * mach).sqrt();
    1.0 / beta
}

/// Get the critical Mach number where transonic drag rise begins
/// 
/// The critical Mach number is where local flow acceleration first reaches Mach 1,
/// causing shock wave formation even though freestream velocity is still subsonic.
/// These values are based on wind tunnel data and computational fluid dynamics.
fn critical_mach_number(shape: ProjectileShape) -> f64 {
    match shape {
        // Sharp, pointed nose delays shock formation due to gradual flow acceleration
        // Range: 0.83-0.87 depending on nose sharpness angle
        // Source: Aberdeen Proving Ground wind tunnel data
        ProjectileShape::Spitzer => 0.85,
        
        // Blunt nose causes earlier shock formation due to rapid flow deceleration
        // Range: 0.72-0.78 depending on nose radius
        // Source: Classical aerodynamics literature (Hoerner, 1965)
        ProjectileShape::RoundNose => 0.75,
        
        // Flat base creates additional pressure drag and early transonic effects
        // Range: 0.68-0.72, includes base drag contributions
        // Source: Experimental ballistics data
        ProjectileShape::FlatBase => 0.70,
        
        // Tapered base (boat tail) reduces pressure drag, delays transonic rise
        // Range: 0.86-0.90 depending on boat tail angle (typically 7-9 degrees)
        // Source: Modern CFD validation studies
        ProjectileShape::BoatTail => 0.88,
    }
}

/// Calculate the transonic drag rise factor
/// 
/// This models the sharp increase in drag as shock waves form and strengthen
/// in the transonic regime. Based on empirical correlations and theoretical
/// models from Anderson's "Modern Compressible Flow" and McCoy's work.
pub fn transonic_drag_rise(mach: f64, shape: ProjectileShape) -> f64 {
    let m_crit = critical_mach_number(shape);
    
    if mach < m_crit {
        // Below critical Mach, no significant drag rise
        return 1.0;
    }
    
    if mach < 1.0 {
        // Subsonic drag rise (shock waves forming locally)
        // Use a more physically accurate model based on empirical data
        
        // Progress through the drag rise region with safe division
        let denominator = 1.0 - m_crit;
        if denominator.abs() < f64::EPSILON {
            return 1.0; // No drag rise if critical Mach is 1.0
        }
        let progress = (mach - m_crit) / denominator;
        
        if progress < 0.0 {
            return 1.0;
        }
        
        // Different rise rates for different shapes
        let rise_factor = match shape {
            ProjectileShape::BoatTail => {
                // Boat tail has gentler drag rise
                1.0 + 1.2 * progress.powi(2)
            },
            ProjectileShape::RoundNose => {
                // Round nose has steeper drag rise
                1.0 + 2.0 * progress.powf(1.5)
            },
            _ => {
                // Spitzer is in between
                1.0 + 1.5 * progress.powf(1.8)
            }
        };
        
        // Add compressibility effects near Mach 1
        let rise_factor = if mach > 0.92 {
            // Smoother transition
            let comp_progress = (mach - 0.92) / 0.08;
            let compressibility = 1.0 + 0.5 * comp_progress.powi(3);
            rise_factor * compressibility
        } else {
            rise_factor
        };
        
        rise_factor.min(2.5)  // More realistic cap
    } else if mach < 1.2 {
        // Transonic/early supersonic (bow shock forming)
        // Peak drag occurs around Mach 1.0-1.1
        
        // Shape-dependent peak location
        let peak_mach = match shape {
            ProjectileShape::Spitzer => 1.05,
            _ => 1.02,
        };
        
        if mach <= peak_mach {
            // Rising to peak
            let base_rise = 1.8; // More realistic peak around 1.8x
            let shape_factor = match shape {
                ProjectileShape::RoundNose => 1.2,
                _ => 1.0,
            };
            base_rise * shape_factor * (1.0 + 0.3 * (mach - 1.0))
        } else {
            // Descending from peak
            let peak_drag = match shape {
                ProjectileShape::RoundNose => 2.2,
                _ => 1.8,
            };
            let decline_rate = 3.0; // How fast drag drops after peak
            peak_drag * (-(decline_rate * (mach - peak_mach))).exp()
        }
    } else {
        // Supersonic (Mach > 1.2)
        // Drag decreases with Mach number as shock angle decreases
        // This is handled by the base drag tables
        1.0
    }
}

/// Calculate the wave drag coefficient component
/// 
/// Wave drag is the additional drag caused by shock waves in transonic/supersonic flow.
/// This is additive to the base drag coefficient.
fn wave_drag_coefficient(mach: f64, shape: ProjectileShape) -> f64 {
    if mach < 0.8 {
        return 0.0;
    }
    
    if mach < 1.0 {
        // Subsonic wave drag (local shocks only)
        let m_crit = critical_mach_number(shape);
        if mach < m_crit {
            return 0.0;
        }
        
        // Gradual onset of wave drag with safe division
        let denominator = 1.0 - m_crit;
        if denominator.abs() < f64::EPSILON {
            return 0.0; // No wave drag if critical Mach is 1.0
        }
        let progress = (mach - m_crit) / denominator;
        let max_subsonic_wave = 0.1;
        max_subsonic_wave * progress.powi(2)
    } else {
        // Supersonic wave drag
        // Based on modified Whitcomb area rule
        let fineness_ratio = 3.5; // Typical for bullets (length/diameter)
        
        // Wave drag coefficient for cone-cylinder
        let cd_wave_base = 0.15 / fineness_ratio;
        
        // Mach number correction (wave drag decreases with Mach)
        // Avoid division by zero at Mach 1.0
        let mach_factor = 1.0 / (mach * mach - 1.0).max(crate::constants::MIN_MACH_THRESHOLD).sqrt();
        
        // Shape correction
        let shape_factor = match shape {
            ProjectileShape::Spitzer => 0.8,     // Good for wave drag
            ProjectileShape::RoundNose => 1.2,   // Poor for wave drag
            ProjectileShape::FlatBase => 1.5,    // Very poor
            ProjectileShape::BoatTail => 0.7,    // Best for wave drag
        };
        
        cd_wave_base * mach_factor * shape_factor
    }
}

/// Apply transonic corrections to a base drag coefficient
/// 
/// This is the main function to use for correcting drag coefficients
/// in the transonic regime.
pub fn transonic_correction(
    mach: f64,
    base_cd: f64,
    shape: ProjectileShape,
    include_wave_drag: bool,
) -> f64 {
    // Get the drag rise factor
    let rise_factor = transonic_drag_rise(mach, shape);
    
    // Apply to base drag
    let mut corrected_cd = base_cd * rise_factor;
    
    // Add wave drag if requested
    if include_wave_drag && mach > 0.8 {
        let wave_cd = wave_drag_coefficient(mach, shape);
        corrected_cd += wave_cd;
    }
    
    corrected_cd
}

/// Estimate projectile shape from physical parameters
/// 
/// This is a simple heuristic based on typical bullet designs.
pub fn get_projectile_shape(caliber: f64, weight_grains: f64, g_model: &str) -> ProjectileShape {
    // G7 is typically used for boat tail bullets
    if g_model == "G7" {
        return ProjectileShape::BoatTail;
    }
    
    // Heavy for caliber often means longer, boat tail design
    let weight_per_caliber = weight_grains / caliber;
    if weight_per_caliber > 500.0 {
        // e.g., 175gr .308
        return ProjectileShape::BoatTail;
    }
    
    // Default to spitzer for most rifle bullets
    if caliber < 0.35 {
        // Rifle calibers
        ProjectileShape::Spitzer
    } else {
        // Larger calibers often have round nose
        ProjectileShape::RoundNose
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prandtl_glauert() {
        // Test subsonic values
        assert!((prandtl_glauert_correction(0.5) - 1.1547).abs() < 0.001);
        assert!((prandtl_glauert_correction(0.8) - 1.6667).abs() < 0.001);
        assert!((prandtl_glauert_correction(0.95) - 3.2026).abs() < 0.001);
        
        // Test near Mach 1 capping
        assert_eq!(prandtl_glauert_correction(0.99), 10.0);
    }

    #[test]
    fn test_critical_mach() {
        assert_eq!(critical_mach_number(ProjectileShape::Spitzer), 0.85);
        assert_eq!(critical_mach_number(ProjectileShape::BoatTail), 0.88);
        assert_eq!(critical_mach_number(ProjectileShape::FlatBase), 0.70);
    }

    #[test]
    fn test_transonic_drag_rise() {
        let shape = ProjectileShape::Spitzer;
        
        // Below critical Mach
        assert_eq!(transonic_drag_rise(0.8, shape), 1.0);
        
        // In transonic rise
        let rise_0_9 = transonic_drag_rise(0.9, shape);
        assert!(rise_0_9 > 1.0 && rise_0_9 < 2.0);
        
        // Near Mach 1
        let rise_0_98 = transonic_drag_rise(0.98, shape);
        assert!(rise_0_98 > 2.0);
        
        // Past peak
        let rise_1_1 = transonic_drag_rise(1.1, shape);
        assert!(rise_1_1 > 1.5 && rise_1_1 < 2.5);
    }

    #[test]
    fn test_projectile_shape_estimation() {
        // G7 should always return boat tail
        assert_eq!(get_projectile_shape(0.308, 175.0, "G7"), ProjectileShape::BoatTail);
        
        // Heavy for caliber
        assert_eq!(get_projectile_shape(0.308, 200.0, "G1"), ProjectileShape::BoatTail);
        
        // Standard rifle bullet
        assert_eq!(get_projectile_shape(0.224, 55.0, "G1"), ProjectileShape::Spitzer);
        
        // Large caliber
        assert_eq!(get_projectile_shape(0.50, 300.0, "G1"), ProjectileShape::RoundNose);
    }
}

/// Python-exposed function for transonic drag correction
#[pyfunction]
#[pyo3(name = "transonic_correction_rust", signature = (mach, base_cd, shape_str=None, include_wave_drag=None))]
pub fn transonic_correction_py(
    mach: f64,
    base_cd: f64,
    shape_str: Option<&str>,
    include_wave_drag: Option<bool>,
) -> PyResult<f64> {
    let shape = shape_str
        .map(ProjectileShape::from_str)
        .unwrap_or(ProjectileShape::Spitzer);
    let include_wave = include_wave_drag.unwrap_or(true);
    
    Ok(transonic_correction(mach, base_cd, shape, include_wave))
}

/// Python-exposed function for projectile shape estimation
#[pyfunction]
#[pyo3(name = "get_projectile_shape_rust")]
pub fn get_projectile_shape_py(
    caliber: f64,
    weight_grains: f64,
    g_model: &str,
) -> PyResult<String> {
    let shape = get_projectile_shape(caliber, weight_grains, g_model);
    let shape_str = match shape {
        ProjectileShape::Spitzer => "spitzer",
        ProjectileShape::RoundNose => "round_nose",
        ProjectileShape::FlatBase => "flat_base",
        ProjectileShape::BoatTail => "boat_tail",
    };
    Ok(shape_str.to_string())
}