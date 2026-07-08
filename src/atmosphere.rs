//! Enhanced atmospheric calculations for ballistics.
//!
//! This module provides Rust-accelerated implementations of atmospheric calculations
//! with full ICAO Standard Atmosphere support for improved accuracy at all altitudes.

/// ICAO Standard Atmosphere layer definitions
#[derive(Debug, Clone)]
struct AtmosphereLayer {
    /// Base altitude of this layer (m)
    base_altitude: f64,
    /// Base temperature at layer start (K)
    base_temperature: f64,
    /// Base pressure at layer start (Pa)
    base_pressure: f64,
    /// Temperature lapse rate (K/m)
    lapse_rate: f64,
}

/// ICAO Standard Atmosphere constants
const G_ACCEL_MPS2: f64 = 9.80665;
const R_AIR: f64 = 287.0531; // Specific gas constant for dry air (J/(kg·K))
const GAMMA: f64 = 1.4; // Heat capacity ratio for air
const R_DRY: f64 = 287.05; // Gas constant for dry air
const R_VAPOR: f64 = 461.495; // Gas constant for water vapor

/// CIPM constants for precise air density calculation
const R: f64 = 8.314472; // Universal gas constant
const M_A: f64 = 28.96546e-3; // Molar mass of dry air (kg/mol)
const M_V: f64 = 18.01528e-3; // Molar mass of water vapor (kg/mol)

/// ICAO Standard Atmosphere layer data up to 84 km
/// Pressures calculated using barometric formula between layers
const ICAO_LAYERS: &[AtmosphereLayer] = &[
    // Troposphere (0 - 11 km)
    AtmosphereLayer {
        base_altitude: 0.0,
        base_temperature: 288.15, // 15°C
        base_pressure: 101325.0,  // 1013.25 hPa
        lapse_rate: -0.0065,      // -6.5 K/km
    },
    // Tropopause (11 - 20 km)
    AtmosphereLayer {
        base_altitude: 11000.0,
        base_temperature: 216.65, // -56.5°C
        base_pressure: 22632.1,   // 226.32 hPa
        lapse_rate: 0.0,          // Isothermal
    },
    // Stratosphere 1 (20 - 32 km)
    AtmosphereLayer {
        base_altitude: 20000.0,
        base_temperature: 216.65, // -56.5°C
        base_pressure: 5474.89,   // 54.75 hPa
        lapse_rate: 0.001,        // +1 K/km
    },
    // Stratosphere 2 (32 - 47 km)
    AtmosphereLayer {
        base_altitude: 32000.0,
        base_temperature: 228.65, // -44.5°C
        base_pressure: 868.02,    // 8.68 hPa
        lapse_rate: 0.0028,       // +2.8 K/km
    },
    // Stratopause (47 - 51 km)
    AtmosphereLayer {
        base_altitude: 47000.0,
        base_temperature: 270.65, // -2.5°C
        base_pressure: 110.91,    // 1.11 hPa
        lapse_rate: 0.0,          // Isothermal
    },
    // Mesosphere 1 (51 - 71 km)
    AtmosphereLayer {
        base_altitude: 51000.0,
        base_temperature: 270.65, // -2.5°C
        base_pressure: 66.94,     // 0.67 hPa
        lapse_rate: -0.0028,      // -2.8 K/km
    },
    // Mesosphere 2 (71 - 84 km)
    AtmosphereLayer {
        base_altitude: 71000.0,
        base_temperature: 214.65, // -58.5°C
        base_pressure: 3.96,      // 0.04 hPa
        lapse_rate: -0.002,       // -2.0 K/km
    },
];

/// Calculate ICAO Standard Atmosphere conditions at any altitude.
///
/// This function implements the full ICAO Standard Atmosphere model with all
/// atmospheric layers up to 84 km altitude.
///
/// # Arguments
/// * `altitude_m` - Altitude in meters (0 to 84000)
///
/// # Returns
/// Tuple of (temperature_k, pressure_pa)
fn calculate_icao_standard_atmosphere(altitude_m: f64) -> (f64, f64) {
    // Clamp altitude to valid range
    let altitude = altitude_m.clamp(0.0, 84000.0);

    // Find the appropriate atmospheric layer
    let layer = ICAO_LAYERS
        .iter()
        .rev()
        .find(|layer| altitude >= layer.base_altitude)
        .unwrap_or(&ICAO_LAYERS[0]);

    let height_diff = altitude - layer.base_altitude;
    let temperature = layer.base_temperature + layer.lapse_rate * height_diff;

    let pressure = if layer.lapse_rate.abs() < 1e-10 {
        // Isothermal layer
        layer.base_pressure * (-G_ACCEL_MPS2 * height_diff / (R_AIR * layer.base_temperature)).exp()
    } else {
        // Non-isothermal layer
        let temp_ratio = temperature / layer.base_temperature;
        layer.base_pressure * temp_ratio.powf(-G_ACCEL_MPS2 / (layer.lapse_rate * R_AIR))
    };

    (temperature, pressure)
}

/// Resolve the station-pressure override for an air-density calculation.
///
/// Altitude and pressure are redundant inputs for density. The rule:
/// * An explicitly-supplied pressure is the authoritative STATION pressure (already
///   altitude-reduced); it is returned as `Some` and used directly, so altitude is NOT
///   double-counted.
/// * When pressure is left at the sea-level standard default (≈1013.25 hPa) while a real
///   altitude is given, the caller meant "standard atmosphere at this altitude": return
///   `None` so [`calculate_atmosphere`] derives the station pressure from altitude (ICAO
///   standard) instead of silently using sea-level density.
///
/// Without this, `--altitude` with the default pressure produced sea-level density (altitude
/// had no effect on drag). The ±0.5 hPa tolerance covers the `29.92 inHg ≈ 1013.21 hPa`
/// conversion, and `>1 m` avoids triggering at sea level. (Mirrors the existing
/// `pressure != 29.92` "user override" sentinel used elsewhere in the CLI.)
pub fn resolve_station_pressure(pressure_hpa: f64, altitude_m: f64) -> Option<f64> {
    const SEA_LEVEL_HPA: f64 = 1013.25;
    if (pressure_hpa - SEA_LEVEL_HPA).abs() < 0.5 && altitude_m.abs() > 1.0 {
        None // pressure left at default + real altitude → derive station pressure from altitude
    } else {
        Some(pressure_hpa) // explicit station pressure is authoritative
    }
}

/// Resolve the temperature override for an air-density calculation, mirroring
/// [`resolve_station_pressure`].
///
/// * An explicitly-supplied temperature is authoritative (returned as `Some`).
/// * When temperature is left at the sea-level standard default (15 °C) while a real altitude
///   is given, the caller meant "standard atmosphere at this altitude": return `None` so
///   [`calculate_atmosphere`] applies the ICAO lapse-rate temperature for that altitude
///   (≈ −6.5 °C/km).
///
/// Without this, `--altitude` with the default temperature held the air at 15 °C, which
/// under-estimates density (warm air is thinner) by ~2.4% at 1 km up to ~7% at 3 km versus the
/// standard atmosphere — validated against py_ballisticcalc, which derives both temperature and
/// pressure from altitude. The 0.1 °C tolerance matches the `59 °F = 15.0 °C` default exactly,
/// and `>1 m` avoids triggering at sea level. A shooter at a genuinely non-standard temperature
/// at altitude should pass an explicit temperature (same contract as station pressure).
pub fn resolve_station_temperature(temperature_c: f64, altitude_m: f64) -> Option<f64> {
    const SEA_LEVEL_TEMP_C: f64 = 15.0;
    if (temperature_c - SEA_LEVEL_TEMP_C).abs() < 0.1 && altitude_m.abs() > 1.0 {
        None // temperature left at default + real altitude → derive ICAO lapse temperature
    } else {
        Some(temperature_c) // explicit temperature is authoritative
    }
}

/// Return the station temperature and pressure that [`calculate_atmosphere`] will use after
/// applying the default-at-altitude resolution rules.
pub fn resolve_station_conditions(
    temperature_c: f64,
    pressure_hpa: f64,
    altitude_m: f64,
) -> (f64, f64) {
    let temp_override = resolve_station_temperature(temperature_c, altitude_m);
    let press_override = resolve_station_pressure(pressure_hpa, altitude_m);
    let (std_temp_k, std_pressure_pa) = calculate_icao_standard_atmosphere(altitude_m);
    let temp_c = temp_override.unwrap_or(std_temp_k - 273.15);
    let pressure_hpa = press_override.unwrap_or(std_pressure_pa / 100.0);
    (temp_c, pressure_hpa)
}

/// Enhanced atmospheric calculation with ICAO Standard Atmosphere.
///
/// # Arguments
/// * `altitude_m` - Altitude in meters
/// * `temp_override_c` - Temperature override in Celsius (None for standard)
/// * `press_override_hpa` - Pressure override in hPa (None for standard)
/// * `humidity_percent` - Humidity percentage (0-100)
///
/// # Returns
/// Tuple of (air_density_kg_m3, speed_of_sound_mps)
pub fn calculate_atmosphere(
    altitude_m: f64,
    temp_override_c: Option<f64>,
    press_override_hpa: Option<f64>,
    humidity_percent: f64,
) -> (f64, f64) {
    // Get standard atmosphere conditions or use overrides
    let (temp_k, pressure_pa) = if temp_override_c.is_some() && press_override_hpa.is_some() {
        // Both overrides provided
        (
            temp_override_c.unwrap() + 273.15,
            press_override_hpa.unwrap() * 100.0,
        )
    } else {
        // Get ICAO standard conditions
        let (std_temp_k, std_pressure_pa) = calculate_icao_standard_atmosphere(altitude_m);

        let final_temp_k = if let Some(temp_c) = temp_override_c {
            temp_c + 273.15
        } else {
            std_temp_k
        };

        let final_pressure_pa = if let Some(press_hpa) = press_override_hpa {
            press_hpa * 100.0
        } else {
            std_pressure_pa
        };

        (final_temp_k, final_pressure_pa)
    };

    // Enhanced humidity effects on air density and speed of sound
    let humidity_clamped = humidity_percent.clamp(0.0, 100.0);

    // Calculate saturation vapor pressure (enhanced Magnus formula)
    let temp_c = temp_k - 273.15;
    let es_hpa = if temp_c >= 0.0 {
        // Over water (Arden Buck): es = 6.1121 * exp[(18.678 - T/234.5) * (T/(257.14+T))].
        // The ENTIRE product is the exponent; previously the linear factor sat OUTSIDE exp,
        // over-estimating es by ~7x at 15C and corrupting humidity-dependent air density.
        6.1121 * ((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c))).exp()
    } else {
        // Over ice (Arden Buck): es = 6.1115 * exp[(23.036 - T/333.7) * (T/(279.82+T))].
        6.1115 * ((23.036 - temp_c / 333.7) * (temp_c / (279.82 + temp_c))).exp()
    };

    // Calculate actual vapor pressure
    let vapor_pressure_pa = humidity_clamped / 100.0 * es_hpa * 100.0;
    let dry_pressure_pa = (pressure_pa - vapor_pressure_pa).max(0.0);

    // Calculate air density with humidity effects
    let density = dry_pressure_pa / (R_DRY * temp_k) + vapor_pressure_pa / (R_VAPOR * temp_k);

    // Enhanced speed of sound calculation with humidity effects
    // Speed of sound in moist air (Cramer, 1993)
    // Guard pressure_pa == 0 (a 0 hPa override would otherwise give +Inf -> -Inf gamma -> NaN
    // speed of sound) and cap the mole fraction at the physical maximum of 1.
    let mole_fraction_vapor = (vapor_pressure_pa / pressure_pa.max(f64::MIN_POSITIVE)).min(1.0);
    let temp_c_abs = temp_k;

    // Heat capacity ratio for moist air. The coefficient is for vapor mole fraction.
    let gamma_moist = GAMMA * (1.0 - mole_fraction_vapor * 0.062);

    // Gas constant for moist air. 0.6078 belongs to specific humidity; for mole fraction the
    // dry-air molecular-weight ratio gives approximately 0.378.
    let r_moist = R_AIR * (1.0 + 0.378 * mole_fraction_vapor);

    // Speed of sound in moist air (Cramer, 1993) — physics-based correction only. The
    // gamma_moist / r_moist terms above already yield the correct humid speed of sound; the
    // previous extra empirical humidity_correction factor double-counted the humidity effect
    // (roughly doubling it), over-predicting the public speed_of_sound and the reported MC Mach.
    let speed_of_sound = (gamma_moist * r_moist * temp_c_abs).sqrt();

    (density, speed_of_sound)
}

/// Enhanced air density calculation using CIPM formula with ICAO atmosphere.
///
/// # Arguments
/// * `temp_c` - Temperature in Celsius
/// * `pressure_hpa` - Pressure in hPa
/// * `humidity_percent` - Humidity percentage (0-100)
///
/// # Returns
/// Air density in kg/m³
pub fn calculate_air_density_cimp(temp_c: f64, pressure_hpa: f64, humidity_percent: f64) -> f64 {
    let t_k = temp_c + 273.15;

    // Enhanced saturation vapor pressure calculation
    let p_sv = enhanced_saturation_vapor_pressure(t_k);

    let pressure_pa = pressure_hpa * 100.0;

    // Enhanced enhancement factor with temperature dependence. CIPM constants use Pa.
    let f = enhanced_enhancement_factor(pressure_pa, temp_c);

    // Vapor pressure with clamping. p_sv is in hPa (enhanced_saturation_vapor_pressure
    // returns hPa — its critical-pressure constant is 220640 hPa), so p_v is in hPa too.
    let p_v = humidity_percent.clamp(0.0, 100.0) / 100.0 * f * p_sv;

    // Convert the vapor pressure to Pa BEFORE forming the mole fraction: the divisor below
    // is in Pa. Dividing the hPa p_v by the Pa total made x_v 100x too small, which erased
    // the humidity term and returned essentially dry-air density (e.g. 15 C / 1013.25 hPa /
    // 50% RH gave ~1.2254 instead of the CIPM-2007 moist value ~1.2211 — moist air is
    // LIGHTER than dry air).
    let p_v_pa = p_v * 100.0;

    // Floor the pressure divisor (mirrors calculate_atmosphere): a 0 hPa pressure would
    // otherwise make x_v = +Inf -> NaN density. No-op for all valid (>0) pressures.
    let p_pa = pressure_pa.max(f64::MIN_POSITIVE);

    // Mole fraction of water vapor (capped at the physical maximum of 1)
    let x_v = (p_v_pa / p_pa).min(1.0);

    // Enhanced compressibility factor. CIPM virial constants use Pa.
    let z = enhanced_compressibility_factor(p_pa, t_k, x_v);

    // Calculate density with enhanced precision
    // Note: parentheses are important here for correct operator precedence
    ((p_pa * M_A) / (z * R * t_k)) * (1.0 - x_v * (1.0 - M_V / M_A))
}

/// Enhanced saturation vapor pressure calculation.
/// Uses the IAPWS-IF97 formulation for high precision.
#[inline(always)]
fn enhanced_saturation_vapor_pressure(t_k: f64) -> f64 {
    // IAPWS-IF97 coefficients for better accuracy
    const A: [f64; 6] = [
        -7.85951783,
        1.84408259,
        -11.7866497,
        22.6807411,
        -15.9618719,
        1.80122502,
    ];

    // Ensure temperature is positive and reasonable
    let t_k_safe = t_k.max(173.15); // -100°C minimum

    let tau = 1.0 - t_k_safe / 647.096; // Critical temperature of water
    let ln_p_ratio = (647.096 / t_k_safe)
        * (A[0] * tau
            + A[1] * tau.powf(1.5)
            + A[2] * tau.powf(3.0)
            + A[3] * tau.powf(3.5)
            + A[4] * tau.powf(4.0)
            + A[5] * tau.powf(7.5));

    220640.0 * ln_p_ratio.exp() // Critical pressure in hPa (22.064 MPa)
}

/// Enhanced enhancement factor with altitude and temperature dependence.
#[inline(always)]
fn enhanced_enhancement_factor(p: f64, t: f64) -> f64 {
    const ALPHA: f64 = 1.00062;
    const BETA: f64 = 3.14e-8;
    const GAMMA: f64 = 5.6e-7;
    const DELTA: f64 = 1.2e-10; // Additional temperature term

    ALPHA + BETA * p + GAMMA * t * t + DELTA * p * t
}

/// Enhanced compressibility factor with improved accuracy.
#[inline(always)]
fn enhanced_compressibility_factor(p: f64, t_k: f64, x_v: f64) -> f64 {
    // Enhanced virial coefficients for better accuracy
    const A0: f64 = 1.58123e-6;
    const A1: f64 = -2.9331e-8;
    const A2: f64 = 1.1043e-10;
    const B0: f64 = 5.707e-6;
    const B1: f64 = -2.051e-8;
    const C0: f64 = 1.9898e-4;
    const C1: f64 = -2.376e-6;
    const D: f64 = 1.83e-11;
    const E: f64 = -0.765e-8;

    // Additional third-order terms for enhanced accuracy
    const F0: f64 = 2.1e-12;
    const F1: f64 = -1.1e-14;

    // Ensure temperature is positive
    let t_k_safe = t_k.max(173.15); // -100°C minimum
    let t = t_k_safe - 273.15;
    let p_t = p / t_k_safe;

    let z_second_order =
        1.0 - p_t * (A0 + A1 * t + A2 * t * t + (B0 + B1 * t) * x_v + (C0 + C1 * t) * x_v * x_v);

    let z_third_order = p_t * p_t * (D + E * x_v * x_v);

    // Enhanced fourth-order correction
    let z_fourth_order = p_t * p_t * p_t * (F0 + F1 * x_v * x_v * x_v);

    z_second_order + z_third_order + z_fourth_order
}

/// Enhanced local atmospheric calculation with variable lapse rates.
///
/// # Arguments
/// * `altitude_m` - Altitude in meters
/// * `base_alt` - Base altitude for calculation
/// * `base_temp_c` - Base temperature in Celsius
/// * `base_press_hpa` - Base pressure in hPa
/// * `base_ratio` - Base density ratio
///
/// # Returns
/// Tuple of (air_density_kg_m3, speed_of_sound_mps)
pub fn get_local_atmosphere(
    altitude_m: f64,
    base_alt: f64,
    base_temp_c: f64,
    base_press_hpa: f64,
    base_ratio: f64,
) -> (f64, f64) {
    // Round altitude to the nearest meter for caching in Python
    let altitude_m_rounded = altitude_m.round();
    let height_diff = altitude_m_rounded - base_alt;

    // Determine appropriate lapse rate based on altitude
    let lapse_rate = determine_local_lapse_rate(altitude_m_rounded);

    // Calculate temperature with variable lapse rate
    let temp_c = base_temp_c + lapse_rate * height_diff;
    let temp_k = temp_c + 273.15;
    let base_temp_k = base_temp_c + 273.15;

    // Calculate pressure using barometric formula
    let pressure_hpa = if lapse_rate.abs() < 1e-10 {
        // Isothermal atmosphere
        base_press_hpa * (-G_ACCEL_MPS2 * height_diff / (R_AIR * base_temp_k)).exp()
    } else {
        // Non-isothermal atmosphere
        let temp_ratio = temp_k / base_temp_k;
        base_press_hpa * temp_ratio.powf(-G_ACCEL_MPS2 / (lapse_rate * R_AIR))
    };

    // Enhanced density calculation
    let density_ratio = base_ratio * (base_temp_k * pressure_hpa) / (base_press_hpa * temp_k);
    let density = density_ratio * 1.225;

    // Enhanced speed of sound calculation
    let speed_of_sound = (temp_k * 401.874).sqrt(); // More precise constant

    (density, speed_of_sound)
}

/// Determine local lapse rate based on altitude and atmospheric layer.
#[inline(always)]
fn determine_local_lapse_rate(altitude_m: f64) -> f64 {
    // Find the current atmospheric layer to get appropriate lapse rate
    let layer = ICAO_LAYERS
        .iter()
        .rev()
        .find(|layer| altitude_m >= layer.base_altitude)
        .unwrap_or(&ICAO_LAYERS[0]);

    layer.lapse_rate
}

/// Direct atmosphere calculation for simple cases.
///
/// # Arguments
/// * `density` - Pre-computed air density
/// * `speed_of_sound` - Pre-computed speed of sound
///
/// # Returns
/// Tuple of (air_density, speed_of_sound) - just passes through the values
#[inline(always)]
pub fn get_direct_atmosphere(density: f64, speed_of_sound: f64) -> (f64, f64) {
    (density, speed_of_sound)
}

/// Legacy function name for backwards compatibility
pub fn calculate_air_density_cipm(temp_c: f64, pressure_hpa: f64, humidity_percent: f64) -> f64 {
    calculate_air_density_cimp(temp_c, pressure_hpa, humidity_percent)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_icao_standard_atmosphere() {
        // Test sea level
        let (temp, press) = calculate_icao_standard_atmosphere(0.0);
        assert!((temp - 288.15).abs() < 0.01);
        assert!((press - 101325.0).abs() < 1.0);

        // Test tropopause
        let (temp_11km, press_11km) = calculate_icao_standard_atmosphere(11000.0);
        assert!((temp_11km - 216.65).abs() < 0.01);
        assert!(press_11km < 101325.0);

        // Test stratosphere
        let (temp_25km, _) = calculate_icao_standard_atmosphere(25000.0);
        assert!(temp_25km > 216.65); // Temperature increases in stratosphere
    }

    #[test]
    fn test_enhanced_atmosphere_sea_level() {
        let (density, speed) = calculate_atmosphere(0.0, None, None, 0.0);
        assert!((density - 1.225).abs() < 0.01);
        assert!((speed - 340.0).abs() < 1.0);
    }

    #[test]
    fn test_resolve_station_pressure_contract() {
        // Default sea-level pressure + real altitude => derive from altitude (None).
        assert_eq!(resolve_station_pressure(1013.25, 2000.0), None);
        // 29.92 inHg ≈ 1013.21 hPa is also treated as the default (within tolerance).
        assert_eq!(resolve_station_pressure(1013.21, 2000.0), None);
        // An explicit, non-default station pressure is authoritative (Some, used directly).
        assert_eq!(resolve_station_pressure(850.0, 2000.0), Some(850.0));
        // At/near sea level the default is used directly (no derivation needed).
        assert_eq!(resolve_station_pressure(1013.25, 0.0), Some(1013.25));
    }

    #[test]
    fn test_altitude_affects_density_with_default_pressure() {
        // Regression: with the default pressure, altitude MUST lower density (previously the
        // air-density path ignored altitude whenever pressure was the sea-level default).
        let press = resolve_station_pressure(1013.25, 0.0);
        let (rho_sea, _) = calculate_atmosphere(0.0, Some(15.0), press, 50.0);
        let press_alt = resolve_station_pressure(1013.25, 2000.0);
        let (rho_2km, _) = calculate_atmosphere(2000.0, Some(15.0), press_alt, 50.0);
        assert!(
            rho_2km < rho_sea * 0.9,
            "density at 2000 m ({rho_2km}) should be well below sea level ({rho_sea})"
        );

        // But an explicit station pressure stays authoritative (no altitude double-count):
        // density with an explicit pressure is independent of the altitude field.
        let p = resolve_station_pressure(900.0, 2000.0);
        let (rho_a, _) = calculate_atmosphere(2000.0, Some(15.0), p, 50.0);
        let (rho_b, _) = calculate_atmosphere(0.0, Some(15.0), p, 50.0);
        assert!(
            (rho_a - rho_b).abs() < 1e-9,
            "explicit pressure must ignore altitude"
        );
    }

    #[test]
    fn test_resolve_station_temperature_contract() {
        // Default 15 C + real altitude => derive ICAO lapse temperature (None).
        assert_eq!(resolve_station_temperature(15.0, 2000.0), None);
        // An explicit, non-default temperature is authoritative (Some, used directly).
        assert_eq!(resolve_station_temperature(-5.0, 2000.0), Some(-5.0));
        assert_eq!(resolve_station_temperature(30.0, 2000.0), Some(30.0));
        // At/near sea level the default is used directly (no derivation needed).
        assert_eq!(resolve_station_temperature(15.0, 0.0), Some(15.0));
    }

    #[test]
    fn test_altitude_only_default_matches_full_icao_standard() {
        // Regression: resolving BOTH temperature and pressure for an altitude-only query (defaults
        // left in place) must equal the fully-standard atmosphere at that altitude — i.e. altitude
        // now drives temperature (ICAO lapse) AND pressure, not just pressure. Validated against
        // py_ballisticcalc to ~0.04%. Previously the air held 15 C, leaving density ~7% too thin
        // (warm) at 3 km.
        for alt in [1000.0, 2000.0, 2500.0, 3000.0] {
            let t = resolve_station_temperature(15.0, alt);
            let p = resolve_station_pressure(1013.25, alt);
            let (rho_resolved, _) = calculate_atmosphere(alt, t, p, 0.0);
            let (rho_std, _) = calculate_atmosphere(alt, None, None, 0.0);
            assert!(
                (rho_resolved - rho_std).abs() < 1e-9,
                "alt {alt}: altitude-only default density {rho_resolved} should equal the full \
                 ICAO standard {rho_std}"
            );
            // And it must be denser than the old temperature-held-at-15C behavior (colder = denser).
            let (rho_warm, _) = calculate_atmosphere(alt, Some(15.0), p, 0.0);
            assert!(
                rho_resolved > rho_warm,
                "alt {alt}: lapse-temperature density {rho_resolved} should exceed 15 C-held {rho_warm}"
            );
        }
    }

    #[test]
    fn test_enhanced_atmosphere_with_humidity() {
        let (density_dry, speed_dry) = calculate_atmosphere(0.0, None, None, 0.0);
        let (density_humid, speed_humid) = calculate_atmosphere(0.0, None, None, 80.0);

        // Humid air should be less dense
        assert!(density_humid < density_dry);
        // Humid air should have slightly higher speed of sound
        assert!(speed_humid > speed_dry);
    }

    #[test]
    fn test_enhanced_atmosphere_stratosphere() {
        // Test in stratosphere where temperature increases
        let (density_20km, speed_20km) = calculate_atmosphere(20000.0, None, None, 0.0);
        let (density_30km, speed_30km) = calculate_atmosphere(30000.0, None, None, 0.0);

        // Density should decrease with altitude
        assert!(density_30km < density_20km);
        // Speed of sound should increase due to temperature increase
        assert!(speed_30km > speed_20km);
    }

    #[test]
    fn test_enhanced_cimp_density() {
        let density = calculate_air_density_cimp(15.0, 1013.25, 0.0);
        assert!((density - 1.225).abs() < 0.01);

        // Test with humidity
        let density_humid = calculate_air_density_cimp(15.0, 1013.25, 50.0);
        assert!(density_humid < density);
    }

    #[test]
    fn test_cipm_moist_air_matches_python_reference() {
        // Regression for the hPa/Pa mole-fraction slip: p_v (hPa) was divided by the total
        // pressure in Pa, making x_v 100x too small and erasing the humidity effect entirely.
        // Reference values computed with the validated Python implementation
        // (ballistics.physics.atmosphere_icao.calculate_air_density_cipm_icao), same cases as
        // the Flask suite's tests/test_atmosphere.py::TestCalculateAirDensityCIPM. Tolerance
        // 0.1% (matches that suite's Rust-vs-Python assertion).
        let cases = [
            (15.0, 1013.25, 50.0, 1.221125867723075),
            (30.0, 1000.0, 80.0, 1.1344071877123691),
            (-10.0, 1020.0, 20.0, 1.3500610713710515),
        ];
        for (temp_c, pressure_hpa, humidity_pct, expected) in cases {
            let density = calculate_air_density_cipm(temp_c, pressure_hpa, humidity_pct);
            let rel_err = ((density - expected) / expected).abs();
            assert!(
                rel_err < 1e-3,
                "CIPM density at {temp_c} C / {pressure_hpa} hPa / {humidity_pct}% RH: \
                 got {density}, expected {expected} (rel err {rel_err:.2e} >= 1e-3)"
            );
        }

        // Moist air must be materially lighter than dry air at the same temp/pressure
        // (the broken version returned a difference of only ~4e-5 kg/m^3).
        let dry = calculate_air_density_cipm(15.0, 1013.25, 0.0);
        let moist = calculate_air_density_cipm(15.0, 1013.25, 50.0);
        assert!(
            dry - moist > 3e-3,
            "humidity effect too small: dry {dry} vs 50% RH {moist}"
        );
    }

    #[test]
    fn test_variable_lapse_rates() {
        // Test that lapse rates change appropriately with altitude
        let lapse_tropo = determine_local_lapse_rate(5000.0);
        let lapse_strato = determine_local_lapse_rate(25000.0);

        assert!((lapse_tropo - (-0.0065)).abs() < 0.0001);
        assert!(lapse_strato > 0.0); // Positive lapse rate in stratosphere
    }
}
