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
const R_AIR: f64 = 287.0531;  // Specific gas constant for dry air (J/(kg·K))
const GAMMA: f64 = 1.4;       // Heat capacity ratio for air
const R_DRY: f64 = 287.05;    // Gas constant for dry air
const R_VAPOR: f64 = 461.495; // Gas constant for water vapor

/// CIPM constants for precise air density calculation
const R: f64 = 8.314472;      // Universal gas constant
const M_A: f64 = 28.96546e-3; // Molar mass of dry air (kg/mol)
const M_V: f64 = 18.01528e-3; // Molar mass of water vapor (kg/mol)

/// ICAO Standard Atmosphere layer data up to 84 km
/// Pressures calculated using barometric formula between layers
const ICAO_LAYERS: &[AtmosphereLayer] = &[
    // Troposphere (0 - 11 km)
    AtmosphereLayer {
        base_altitude: 0.0,
        base_temperature: 288.15,   // 15°C
        base_pressure: 101325.0,    // 1013.25 hPa
        lapse_rate: -0.0065,        // -6.5 K/km
    },
    // Tropopause (11 - 20 km)
    AtmosphereLayer {
        base_altitude: 11000.0,
        base_temperature: 216.65,   // -56.5°C
        base_pressure: 22632.1,     // 226.32 hPa
        lapse_rate: 0.0,            // Isothermal
    },
    // Stratosphere 1 (20 - 32 km)
    AtmosphereLayer {
        base_altitude: 20000.0,
        base_temperature: 216.65,   // -56.5°C
        base_pressure: 5474.89,     // 54.75 hPa
        lapse_rate: 0.001,          // +1 K/km
    },
    // Stratosphere 2 (32 - 47 km)
    AtmosphereLayer {
        base_altitude: 32000.0,
        base_temperature: 228.65,   // -44.5°C
        base_pressure: 868.02,      // 8.68 hPa
        lapse_rate: 0.0028,         // +2.8 K/km
    },
    // Stratopause (47 - 51 km)
    AtmosphereLayer {
        base_altitude: 47000.0,
        base_temperature: 270.65,   // -2.5°C
        base_pressure: 110.91,      // 1.11 hPa
        lapse_rate: 0.0,            // Isothermal
    },
    // Mesosphere 1 (51 - 71 km)
    AtmosphereLayer {
        base_altitude: 51000.0,
        base_temperature: 270.65,   // -2.5°C
        base_pressure: 66.94,       // 0.67 hPa
        lapse_rate: -0.0028,        // -2.8 K/km
    },
    // Mesosphere 2 (71 - 84 km)
    AtmosphereLayer {
        base_altitude: 71000.0,
        base_temperature: 214.65,   // -58.5°C
        base_pressure: 3.96,        // 0.04 hPa
        lapse_rate: -0.002,         // -2.0 K/km
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
    let layer = ICAO_LAYERS.iter()
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
        (temp_override_c.unwrap() + 273.15, press_override_hpa.unwrap() * 100.0)
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
        // Over water (Arden Buck equation)
        6.1121 * (18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)).exp()
    } else {
        // Over ice (Arden Buck equation)
        6.1115 * (23.036 - temp_c / 333.7) * (temp_c / (279.82 + temp_c)).exp()
    };
    
    // Calculate actual vapor pressure
    let vapor_pressure_pa = humidity_clamped / 100.0 * es_hpa * 100.0;
    let dry_pressure_pa = (pressure_pa - vapor_pressure_pa).max(0.0);
    
    // Calculate air density with humidity effects
    let density = dry_pressure_pa / (R_DRY * temp_k) + vapor_pressure_pa / (R_VAPOR * temp_k);
    
    // Enhanced speed of sound calculation with humidity effects
    // Speed of sound in moist air (Cramer, 1993)
    let mole_fraction_vapor = vapor_pressure_pa / pressure_pa;
    let temp_c_abs = temp_k;
    
    // Heat capacity ratio for moist air
    let gamma_moist = GAMMA * (1.0 - mole_fraction_vapor * 0.11);
    
    // Gas constant for moist air
    let r_moist = R_AIR * (1.0 + 0.6078 * mole_fraction_vapor);
    
    // Speed of sound with enhanced humidity correction
    let speed_of_sound_base = (gamma_moist * r_moist * temp_c_abs).sqrt();
    
    // Additional humidity correction for molecular effects
    let humidity_correction = 1.0 + 0.0001 * humidity_clamped * (temp_c / 20.0);
    let speed_of_sound = speed_of_sound_base * humidity_correction;
    
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
pub fn calculate_air_density_cimp(
    temp_c: f64,
    pressure_hpa: f64,
    humidity_percent: f64,
) -> f64 {
    let t_k = temp_c + 273.15;
    
    // Enhanced saturation vapor pressure calculation
    let p_sv = enhanced_saturation_vapor_pressure(t_k);
    
    // Enhanced enhancement factor with temperature dependence
    let f = enhanced_enhancement_factor(pressure_hpa, temp_c);
    
    // Vapor pressure with clamping
    let p_v = humidity_percent.clamp(0.0, 100.0) / 100.0 * f * p_sv;
    
    // Mole fraction of water vapor
    let x_v = p_v / pressure_hpa;
    
    // Enhanced compressibility factor
    let z = enhanced_compressibility_factor(pressure_hpa, t_k, x_v);
    
    // Calculate density with enhanced precision
    // Note: parentheses are important here for correct operator precedence
    let density = ((pressure_hpa * M_A) / (z * R * t_k)) * (1.0 - x_v * (1.0 - M_V / M_A));
    
    // Convert from SI units (pressure in Pa) to final density in kg/m³
    // pressure_hpa is in hPa, so multiply by 100 to get Pa
    density * 100.0
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
        1.80122502
    ];
    
    // Ensure temperature is positive and reasonable
    let t_k_safe = t_k.max(173.15); // -100°C minimum
    
    let tau = 1.0 - t_k_safe / 647.096; // Critical temperature of water
    let ln_p_ratio = (647.096 / t_k_safe) * (
        A[0] * tau +
        A[1] * tau.powf(1.5) +
        A[2] * tau.powf(3.0) +
        A[3] * tau.powf(3.5) +
        A[4] * tau.powf(4.0) +
        A[5] * tau.powf(7.5)
    );
    
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
    
    let z_second_order = 1.0 - p_t * (
        A0 + A1 * t + A2 * t * t + 
        (B0 + B1 * t) * x_v + 
        (C0 + C1 * t) * x_v * x_v
    );
    
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
    let layer = ICAO_LAYERS.iter()
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
pub fn calculate_air_density_cipm(
    temp_c: f64,
    pressure_hpa: f64,
    humidity_percent: f64,
) -> f64 {
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
    fn test_variable_lapse_rates() {
        // Test that lapse rates change appropriately with altitude
        let lapse_tropo = determine_local_lapse_rate(5000.0);
        let lapse_strato = determine_local_lapse_rate(25000.0);
        
        assert!((lapse_tropo - (-0.0065)).abs() < 0.0001);
        assert!(lapse_strato > 0.0); // Positive lapse rate in stratosphere
    }
}