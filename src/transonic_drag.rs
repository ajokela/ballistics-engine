/// Transonic drag modeling with shock wave effects
///
/// This module implements physics-based corrections for drag in the transonic regime
/// (Mach 0.8-1.2) where shock waves significantly affect the drag coefficient.

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectileShape {
    Spitzer,   // Sharp pointed (most common)
    RoundNose, // Blunt/round nose
    FlatBase,  // Flat base (wadcutter)
    BoatTail,  // Boat tail design
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
#[allow(dead_code)]
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

fn sonic_drag_rise(shape: ProjectileShape) -> f64 {
    let shape_factor = match shape {
        ProjectileShape::RoundNose => 1.2,
        _ => 1.0,
    };
    1.8 * shape_factor
}

fn peak_drag_mach(shape: ProjectileShape) -> f64 {
    match shape {
        ProjectileShape::Spitzer => 1.05,
        _ => 1.02,
    }
}

/// Calculate the transonic drag rise factor
///
/// This models the sharp increase in drag as shock waves form and strengthen
/// in the transonic regime. Based on empirical correlations and theoretical
/// models from Anderson's "Modern Compressible Flow" and McCoy's work.
pub fn transonic_drag_rise(mach: f64, shape: ProjectileShape) -> f64 {
    let m_crit = critical_mach_number(shape);
    let sonic_rise = sonic_drag_rise(shape);

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
            }
            ProjectileShape::RoundNose => {
                // Round nose has steeper drag rise
                1.0 + 2.0 * progress.powf(1.5)
            }
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

        // Join the transonic branch continuously at Mach 1. Shape-specific caps avoid the old
        // inverted 2.5 -> 1.8/2.16 sonic step (MBA-1162).
        rise_factor.min(sonic_rise)
    } else if mach < 1.2 {
        // Transonic/early supersonic (bow shock forming)
        // Peak drag occurs around Mach 1.0-1.1

        // Shape-dependent peak location
        let peak_mach = peak_drag_mach(shape);

        if mach <= peak_mach {
            // Rising to peak
            sonic_rise * (1.0 + 0.3 * (mach - 1.0))
        } else {
            // Descend from the actual rising-branch endpoint so the peak join is continuous.
            let peak_drag = sonic_rise * (1.0 + 0.3 * (peak_mach - 1.0));
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
        let mach_factor = 1.0
            / (mach * mach - 1.0)
                .max(crate::constants::MIN_MACH_THRESHOLD)
                .sqrt();

        // Shape correction
        let shape_factor = match shape {
            ProjectileShape::Spitzer => 0.8,   // Good for wave drag
            ProjectileShape::RoundNose => 1.2, // Poor for wave drag
            ProjectileShape::FlatBase => 1.5,  // Very poor
            ProjectileShape::BoatTail => 0.7,  // Best for wave drag
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
    // Standard G1/G7/Gx tables already include their transonic drag rise. Only apply this
    // empirical rise when the caller is explicitly adding wave/transonic effects to a table that
    // does not already contain them.
    let mut corrected_cd = if include_wave_drag {
        base_cd * transonic_drag_rise(mach, shape)
    } else {
        base_cd
    };

    // Add wave drag if requested
    if include_wave_drag && mach > 0.8 {
        let wave_cd = wave_drag_coefficient(mach, shape);
        corrected_cd += wave_cd;
    }

    corrected_cd
}

/// Resolve the projectile shape for transonic corrections (MBA-949). A user-supplied bullet_model
/// name (boat/bt, round/rn, flat/fb) takes precedence; otherwise fall back to the caliber/weight/
/// drag-model heuristic below. This name-parsing previously lived in cli_api and fast_trajectory
/// but NOT derivatives, so named shapes were silently ignored on the integrate_trajectory path.
/// Sharing one helper makes all three solver families resolve the same shape for the same inputs.
pub fn resolve_projectile_shape(
    bullet_model: Option<&str>,
    caliber: f64,
    weight_grains: f64,
    g_model: &str,
) -> ProjectileShape {
    if let Some(model) = bullet_model {
        let m = model.to_lowercase();
        if m.contains("boat") || m.contains("bt") {
            return ProjectileShape::BoatTail;
        } else if m.contains("round") || m.contains("rn") {
            return ProjectileShape::RoundNose;
        } else if m.contains("flat") || m.contains("fb") {
            return ProjectileShape::FlatBase;
        }
    }
    get_projectile_shape(caliber, weight_grains, g_model)
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
    fn test_mba949_resolve_projectile_shape() {
        // Named bullet_model takes precedence (case-insensitive, substring match).
        let c = 0.308;
        let w = 168.0;
        assert!(matches!(
            resolve_projectile_shape(Some("Sierra MatchKing Boat Tail"), c, w, "G1"),
            ProjectileShape::BoatTail
        ));
        assert!(matches!(
            resolve_projectile_shape(Some("300gr RN"), c, w, "G1"),
            ProjectileShape::RoundNose
        ));
        assert!(matches!(
            resolve_projectile_shape(Some("flat base"), c, w, "G1"),
            ProjectileShape::FlatBase
        ));
        // No name -> heuristic. None and an unrecognized name both defer to get_projectile_shape,
        // so they must agree with it for the same caliber/weight/drag model.
        let heuristic = get_projectile_shape(c, w, "G7");
        assert_eq!(resolve_projectile_shape(None, c, w, "G7"), heuristic);
        assert_eq!(
            resolve_projectile_shape(Some("unknown"), c, w, "G7"),
            heuristic
        );
    }

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
        let rise_1_0 = transonic_drag_rise(1.0, shape);
        assert!(rise_0_98 > 1.0 && rise_0_98 <= rise_1_0);

        // Past peak
        let rise_1_1 = transonic_drag_rise(1.1, shape);
        assert!(rise_1_1 > 1.5 && rise_1_1 < 2.5);
    }

    #[test]
    fn transonic_drag_rise_is_continuous_at_sonic_and_peak() {
        const EPSILON: f64 = 1e-9;

        for (shape, peak_mach, expected_sonic, expected_peak) in [
            (ProjectileShape::Spitzer, 1.05, 1.8, 1.827),
            (ProjectileShape::RoundNose, 1.02, 2.16, 2.17296),
            (ProjectileShape::FlatBase, 1.02, 1.8, 1.8108),
            (ProjectileShape::BoatTail, 1.02, 1.8, 1.8108),
        ] {
            let sonic = transonic_drag_rise(1.0, shape);
            let just_below_sonic = transonic_drag_rise(1.0 - EPSILON, shape);
            assert!((sonic - expected_sonic).abs() < 1e-12);
            assert!(
                (just_below_sonic - sonic).abs() < 1e-6,
                "{shape:?} discontinuity at Mach 1: below={just_below_sonic}, at={sonic}"
            );

            let peak = transonic_drag_rise(peak_mach, shape);
            let just_above_peak = transonic_drag_rise(peak_mach + EPSILON, shape);
            assert!((peak - expected_peak).abs() < 1e-12);
            assert!(
                (peak - just_above_peak).abs() < 1e-6,
                "{shape:?} discontinuity at peak Mach {peak_mach}: at={peak}, above={just_above_peak}"
            );
            assert!(peak > sonic, "{shape:?} peak must occur above Mach 1");
            assert!(
                transonic_drag_rise(peak_mach + 0.01, shape) < peak,
                "{shape:?} drag rise must descend after its peak"
            );
        }
    }

    #[test]
    fn test_projectile_shape_estimation() {
        // G7 should always return boat tail
        let shape = get_projectile_shape(0.308, 175.0, "G7");
        assert!(matches!(shape, ProjectileShape::BoatTail));

        // Test that we get valid shapes for various inputs
        let shape1 = get_projectile_shape(0.308, 200.0, "G1");
        assert!(matches!(
            shape1,
            ProjectileShape::Spitzer | ProjectileShape::BoatTail | ProjectileShape::RoundNose
        ));

        let shape2 = get_projectile_shape(0.224, 55.0, "G1");
        assert!(matches!(
            shape2,
            ProjectileShape::Spitzer | ProjectileShape::BoatTail | ProjectileShape::RoundNose
        ));

        let shape3 = get_projectile_shape(0.50, 300.0, "G1");
        assert!(matches!(
            shape3,
            ProjectileShape::Spitzer | ProjectileShape::BoatTail | ProjectileShape::RoundNose
        ));
    }
}

// Removed Python-specific function

// Removed Python-specific function
