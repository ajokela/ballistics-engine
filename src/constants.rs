//! Physical constants used in ballistics calculations.

/// Gravitational acceleration in m/s²
pub const G_ACCEL_MPS2: f64 = 9.80665;

/// Conversion factor: meters per second to feet per second
pub const MPS_TO_FPS: f64 = 3.28084;

/// Conversion factor: feet per second to meters per second
pub const FPS_TO_MPS: f64 = 0.3048;

/// Default powder reference temperature in Fahrenheit.
pub const DEFAULT_POWDER_REFERENCE_TEMP_F: f64 = 70.0;

/// Default powder reference temperature in Celsius.
pub const DEFAULT_POWDER_REFERENCE_TEMP_C: f64 =
    (DEFAULT_POWDER_REFERENCE_TEMP_F - 32.0) * 5.0 / 9.0;

/// Standard air density at sea level (kg/m³)
pub const STANDARD_AIR_DENSITY: f64 = 1.225;

/// Standard sea-level air density under the ICAO Standard Atmosphere (lb/ft^3).
///
/// This is the density basis for [`CD_TO_RETARD`] below — see its doc comment. Most
/// modern published ballistic coefficients (this engine's own retardation formula
/// included) are referenced to this density.
pub const ICAO_DENSITY_LB_FT3: f64 = 0.076474;

/// Standard sea-level air density under the (older) Army Standard Metro atmosphere
/// (lb/ft^3).
///
/// Many published Sierra/Hornady/Barnes ballistic coefficients are referenced to THIS
/// density instead of [`ICAO_DENSITY_LB_FT3`] — quoting the same physical drag at a
/// ~1.8% different reference air density. Not used by this engine's retardation math
/// directly (which is calibrated to the ICAO figure via [`CD_TO_RETARD`]); it exists so
/// [`ASM_TO_ICAO_BC`] can convert an ASM-referenced BC to the ICAO-referenced value the
/// engine expects (`BallisticInputs::bc_reference_standard`, MBA-1365).
pub const ASM_DENSITY_LB_FT3: f64 = 0.075126;

/// Multiplier that converts an Army-Standard-Metro-referenced BC to the ICAO-referenced
/// value this engine's retardation formula expects (MBA-1365):
///
/// `BC_icao = ASM_TO_ICAO_BC * BC_asm`
///
/// A denser reference atmosphere implies a numerically SMALLER published BC for the same
/// physical bullet (the same real drag is attributed to a "worse" — smaller — BC when the
/// assumed air is denser), so this ratio is `ASM_DENSITY_LB_FT3 / ICAO_DENSITY_LB_FT3`
/// (< 1), not its reciprocal. Evaluates to 0.98237 to 5 decimal places (asserted by a
/// unit test in `cli_api.rs`).
pub const ASM_TO_ICAO_BC: f64 = ASM_DENSITY_LB_FT3 / ICAO_DENSITY_LB_FT3;

/// Cd-to-retardation conversion for ICAO-referenced BCs.
///
/// Exact imperial retardation form for density normalized to ICAO sea-level air
/// (1.225 kg/m^3 / [`ICAO_DENSITY_LB_FT3`]):
///
/// `a_ft/s^2 = Cd * v_fps^2 * (rho / 1.225) * CD_TO_RETARD / BC`
///
/// The older `0.000683 * 0.30` value is the Army Standard Metro constant and is
/// only consistent with an [`ASM_DENSITY_LB_FT3`] density reference.
pub const CD_TO_RETARD: f64 = 2.08551e-4;

/// Conversion factor: grains to kilograms
pub const GRAINS_TO_KG: f64 = 0.00006479891;

/// Grams per grain — exact by definition.
///
/// The international avoirdupois pound is defined as exactly 0.45359237 kg
/// (i.e. 453.59237 g), and one pound is exactly 7000 grains, so this value
/// is exact (not a measured or rounded conversion):
/// 453.59237 g / 7000 = 0.06479891 g/grain (equivalently,
/// 0.45359237 kg / 7000 = 0.00006479891 kg/grain, then x1000 to get grams).
///
/// This is the single source of truth for the grain<->gram conversion; do
/// not re-derive or re-round it elsewhere (see `tests/constants_guard.rs`).
pub const GRAMS_PER_GRAIN: f64 = 0.06479891;

/// Grains per gram — the exact reciprocal of [`GRAMS_PER_GRAIN`].
pub const GRAINS_PER_GRAM: f64 = 1.0 / GRAMS_PER_GRAIN;

/// Air density at sea level (kg/m³)
pub const AIR_DENSITY_SEA_LEVEL: f64 = 1.225;

/// Speed of sound at sea level, standard atmospheric conditions
///
/// Value: 340.29 m/s (1116.44 ft/s)
/// Conditions: 15°C (59°F), 1013.25 hPa, dry air
///
/// Temperature dependence: c = 331.3 * sqrt(T_kelvin / 273.15)
///
/// Note: Some calculations use 343.0 m/s (20°C reference) - ensure consistency
/// in Mach number calculations. This value follows ICAO Standard Atmosphere.
///
/// Source: International Standard Atmosphere (ISO 2533)
pub const SPEED_OF_SOUND_MPS: f64 = 340.29;

// Numerical stability constants
/// General numerical tolerance for floating point comparisons
pub const NUMERICAL_TOLERANCE: f64 = 1e-9;

/// Minimum threshold for velocity magnitude to avoid division by zero
pub const MIN_VELOCITY_THRESHOLD: f64 = 1e-6;

/// Minimum threshold for preventing division by zero in general calculations
pub const MIN_DIVISION_THRESHOLD: f64 = 1e-12;

/// Tolerance for root finding algorithms
pub const ROOT_FINDING_TOLERANCE: f64 = 1e-6;

/// Minimum threshold for Mach number calculations near unity
pub const MIN_MACH_THRESHOLD: f64 = 1e-3;

// Ballistic Coefficient (BC) fallback constants
//
// These values are used when BC calculations fail or data is missing.
// Derived from statistical analysis of 2,000+ projectile database.
// Values represent conservative estimates (25th percentile) to avoid
// over-predicting ballistic performance.

/// Conservative overall BC fallback value
///
/// Value: 0.31 (25th percentile from comprehensive ballistics database)
/// Usage: General fallback when no specific projectile data available
/// Methodology: Statistical analysis of measured BC values across all categories
pub const BC_FALLBACK_CONSERVATIVE: f64 = 0.31;

// BC fallback values by projectile weight category (grains).
// Values are based on statistical analysis of ballistic coefficient vs mass relationships.
// Each constant represents the 25th-percentile BC for that weight category.

/// Ultra-light projectiles (0-50 grains)
/// Typical: .17 caliber varmint bullets, .22 caliber target bullets
pub const BC_FALLBACK_ULTRA_LIGHT: f64 = 0.172;

/// Light projectiles (50-100 grains)  
/// Typical: .223 Remington, .243 Winchester hunting bullets
pub const BC_FALLBACK_LIGHT: f64 = 0.242;

/// Medium projectiles (100-150 grains)
/// Typical: .270 Winchester, .30-06 hunting bullets
pub const BC_FALLBACK_MEDIUM: f64 = 0.310;

/// Heavy projectiles (150-200 grains)
/// Typical: .308 Winchester match bullets, .300 Winchester Magnum
pub const BC_FALLBACK_HEAVY: f64 = 0.393;

/// Very heavy projectiles (200+ grains)
/// Typical: .338 Lapua Magnum, .50 BMG bullets
pub const BC_FALLBACK_VERY_HEAVY: f64 = 0.441;

// BC fallback values by caliber category (inches).
// Values account for diameter limitations on achievable ballistic coefficient.
// Larger calibers generally allow higher BC but with diminishing returns.

/// Small calibers (.224" and smaller)
/// Examples: .17 Remington, .22-250, .223 Remington
pub const BC_FALLBACK_SMALL_CALIBER: f64 = 0.215;

/// Medium calibers (.243")
/// Examples: .243 Winchester, 6mm Creedmoor
pub const BC_FALLBACK_MEDIUM_CALIBER: f64 = 0.300;

/// Large calibers (.264" to .284")
/// Examples: .270 Winchester, .280 Remington, 7mm Remington Magnum
pub const BC_FALLBACK_LARGE_CALIBER: f64 = 0.404;

/// Extra large calibers (.308" and larger)
/// Examples: .308 Winchester, .30-06, .300 Winchester Magnum
/// Note: Lower than expected due to inclusion of older, less aerodynamic designs
pub const BC_FALLBACK_XLARGE_CALIBER: f64 = 0.291;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powder_reference_temperature_defaults_are_equivalent() {
        let converted = (DEFAULT_POWDER_REFERENCE_TEMP_F - 32.0) * 5.0 / 9.0;

        assert_eq!(
            DEFAULT_POWDER_REFERENCE_TEMP_C.to_bits(),
            converted.to_bits()
        );
    }

    /// Independent ground-truth check for the grain<->gram constants: this
    /// test does *not* reference `GRAMS_PER_GRAIN`/`GRAINS_PER_GRAM` in its
    /// expected values, so it fails if either constant is ever silently
    /// corrupted (re-truncated, re-rounded, or swapped for a sibling value),
    /// even if the corrupted value doesn't match one of the specific banned
    /// literals `tests/constants_guard.rs` greps for.
    #[test]
    fn grams_per_grain_matches_independent_ground_truth() {
        // 1 lb = exactly 0.45359237 kg = 453.59237 g (avoirdupois pound,
        // international definition); 1 lb = exactly 7000 grains.
        let independent_grams_per_grain: f64 = 453.59237 / 7000.0;
        assert_eq!(
            GRAMS_PER_GRAIN.to_bits(),
            independent_grams_per_grain.to_bits()
        );

        // Cross-check against the pre-existing kilogram-scale constant:
        // GRAMS_PER_GRAIN must equal GRAINS_TO_KG scaled from kg to g.
        assert!((GRAMS_PER_GRAIN - GRAINS_TO_KG * 1000.0).abs() < 1e-15);

        // GRAINS_PER_GRAM must be the exact reciprocal used in production.
        assert_eq!(GRAINS_PER_GRAM.to_bits(), (1.0 / GRAMS_PER_GRAIN).to_bits());
    }
}
