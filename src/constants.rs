/// Physical constants used in ballistics calculations

/// Gravitational acceleration in m/s²
pub const G_ACCEL_MPS2: f64 = 9.80665;

/// Conversion factor: meters per second to feet per second
pub const MPS_TO_FPS: f64 = 3.28084;

/// Conversion factor: feet per second to meters per second
pub const FPS_TO_MPS: f64 = 0.3048;

/// Standard air density at sea level (kg/m³)
pub const STANDARD_AIR_DENSITY: f64 = 1.225;

/// Critical drag coefficient to retardation conversion constant
/// 
/// This fundamental constant converts drag coefficient (Cd) to ballistic retardation force.
/// Value: 0.000683 * 0.30 = 0.0002049
/// 
/// Derivation:
/// - 0.000683: Dimensional conversion factor from imperial ballistics units
/// - 0.30: Empirical correction factor from extensive ballistics testing
/// 
/// Physical meaning: Proportionality constant in the ballistic coefficient equation:
/// BC = (bullet_mass / bullet_diameter²) / (Cd / Cd_standard)
/// Retardation = CD_TO_RETARD * Cd * air_density * velocity²
/// 
/// Sources: Classical ballistics theory (Pejsa, McCoy), validated against
/// Aberdeen Proving Ground data and modern Doppler radar measurements.
pub const CD_TO_RETARD: f64 = 0.000683 * 0.30;

/// Conversion factor: grains to kilograms
pub const GRAINS_TO_KG: f64 = 0.00006479891;

/// Air density at sea level (kg/m³)
pub const AIR_DENSITY_SEA_LEVEL: f64 = 1.225;

/// Speed of sound at sea level, standard atmospheric conditions
/// 
/// Value: 340.29 m/s (1116.8 ft/s)
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

/// BC fallback values by projectile weight category (grains)
/// 
/// Values based on statistical analysis of ballistic coefficient vs mass relationships.
/// Each constant represents 25th percentile BC for that weight category.

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

/// BC fallback values by caliber category (inches)
/// 
/// Values account for diameter limitations on achievable ballistic coefficient.
/// Larger calibers generally allow higher BC but with diminishing returns.

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