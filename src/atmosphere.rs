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

    // Humidity clamp shared by the CIPM density and the moist speed of sound.
    let humidity_clamped = humidity_percent.clamp(0.0, 100.0);
    let temp_c = temp_k - 273.15;

    // Density: CIPM-2007 is the single canonical humid-air density model. Every solver
    // (cli_api / monte_carlo / ffi / fast_trajectory) reaches this one formula through
    // calculate_atmosphere, so there is no second (Arden-Buck ideal-gas) density path to drift.
    let density = calculate_air_density_cimp(temp_c, pressure_pa / 100.0, humidity_clamped);

    // Speed of sound in moist air (Cramer, 1993). Extracted into `moist_speed_of_sound` so the
    // integrators can share it; its vapor pressure comes from the SAME IAPWS saturation formula
    // (`enhanced_saturation_vapor_pressure`) + CIPM enhancement factor that the density above
    // uses, so ONE vapor formula feeds both density and c.
    let speed_of_sound = moist_speed_of_sound(temp_k, pressure_pa, humidity_clamped);

    (density, speed_of_sound)
}

/// Speed of sound in moist air (Cramer, 1993).
///
/// The water-vapor mole fraction is derived from the SAME IAPWS saturation vapor pressure
/// (`enhanced_saturation_vapor_pressure`) and CIPM-2007 enhancement factor used by
/// [`calculate_air_density_cimp`], so a single vapor formula feeds both density and c.
///
/// # Arguments
/// * `temp_k` - Temperature in Kelvin
/// * `pressure_pa` - Total (station) pressure in Pa
/// * `humidity_percent` - Relative humidity percentage (0-100)
///
/// # Returns
/// Speed of sound in m/s
pub fn moist_speed_of_sound(temp_k: f64, pressure_pa: f64, humidity_percent: f64) -> f64 {
    let humidity_clamped = humidity_percent.clamp(0.0, 100.0);
    let temp_c = temp_k - 273.15;

    // Water-vapor partial pressure p_v = RH * f * p_sv, matching CIPM's x_v exactly. p_sv is in
    // hPa (enhanced_saturation_vapor_pressure returns hPa), so convert to Pa before forming the
    // mole fraction against the Pa total pressure.
    let p_sv_hpa = enhanced_saturation_vapor_pressure(temp_k);
    let f = enhanced_enhancement_factor(pressure_pa, temp_c);
    let vapor_pressure_pa = humidity_clamped / 100.0 * f * p_sv_hpa * 100.0;

    // Cap the mole fraction at the physical maximum of 1 and guard pressure_pa == 0 (a 0 hPa
    // override would otherwise give +Inf -> NaN speed of sound).
    let mole_fraction_vapor = (vapor_pressure_pa / pressure_pa.max(f64::MIN_POSITIVE)).min(1.0);

    // Heat-capacity ratio and gas constant for moist air (mole-fraction coefficients). 0.378 is
    // the dry-air molecular-weight ratio (0.6078 would belong to specific humidity, not mole
    // fraction).
    let gamma_moist = GAMMA * (1.0 - mole_fraction_vapor * 0.062);
    let r_moist = R_AIR * (1.0 + 0.378 * mole_fraction_vapor);

    (gamma_moist * r_moist * temp_k).sqrt()
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

/// CIPM-2007 enhancement factor `f = alpha + beta*p + gamma*t^2` (p in Pa, t in Celsius).
#[inline(always)]
fn enhanced_enhancement_factor(p: f64, t: f64) -> f64 {
    const ALPHA: f64 = 1.00062;
    const BETA: f64 = 3.14e-8;
    const GAMMA: f64 = 5.6e-7;

    ALPHA + BETA * p + GAMMA * t * t
}

/// CIPM-2007 compressibility factor `Z` (virial expansion, second order in `p/T`).
#[inline(always)]
fn enhanced_compressibility_factor(p: f64, t_k: f64, x_v: f64) -> f64 {
    // CIPM-2007 molar virial coefficients (p in Pa, t in Celsius).
    const A0: f64 = 1.58123e-6;
    const A1: f64 = -2.9331e-8;
    const A2: f64 = 1.1043e-10;
    const B0: f64 = 5.707e-6;
    const B1: f64 = -2.051e-8;
    const C0: f64 = 1.9898e-4;
    const C1: f64 = -2.376e-6;
    const D: f64 = 1.83e-11;
    const E: f64 = -0.765e-8;

    // Ensure temperature is positive
    let t_k_safe = t_k.max(173.15); // -100°C minimum
    let t = t_k_safe - 273.15;
    let p_t = p / t_k_safe;

    let z_second_order =
        1.0 - p_t * (A0 + A1 * t + A2 * t * t + (B0 + B1 * t) * x_v + (C0 + C1 * t) * x_v * x_v);

    let z_third_order = p_t * p_t * (D + E * x_v * x_v);

    z_second_order + z_third_order
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
    let (temp_k, _pressure_pa, density) =
        local_temp_pressure_density(altitude_m, base_alt, base_temp_c, base_press_hpa, base_ratio);

    // Dry speed of sound. 401.874 ~ gamma * R_air; kept exactly for back-compat with existing
    // callers (get_local_atmosphere_humid uses the precise moist formula instead).
    let speed_of_sound = (temp_k * 401.874).sqrt();

    (density, speed_of_sound)
}

/// Humidity-aware companion to [`get_local_atmosphere`]: identical local density, but the speed
/// of sound is the moist-air value ([`moist_speed_of_sound`]) evaluated at the LOCAL temperature
/// and pressure.
///
/// [`get_local_atmosphere`] is intentionally left unchanged (dry speed of sound) for
/// API/back-compat; call this variant only where a real relative humidity is available.
///
/// # Arguments
/// * `altitude_m` - Query altitude in meters
/// * `base_alt` - Base (station) altitude in meters
/// * `base_temp_c` - Base temperature in Celsius
/// * `base_press_hpa` - Base pressure in hPa
/// * `base_ratio` - Base density ratio (density / 1.225)
/// * `humidity_percent` - Relative humidity percentage (0-100)
///
/// # Returns
/// Tuple of (air_density_kg_m3, moist_speed_of_sound_mps)
pub fn get_local_atmosphere_humid(
    altitude_m: f64,
    base_alt: f64,
    base_temp_c: f64,
    base_press_hpa: f64,
    base_ratio: f64,
    humidity_percent: f64,
) -> (f64, f64) {
    let (temp_k, pressure_pa, density) =
        local_temp_pressure_density(altitude_m, base_alt, base_temp_c, base_press_hpa, base_ratio);
    (density, moist_speed_of_sound(temp_k, pressure_pa, humidity_percent))
}

/// Shared local temperature / pressure / density computation for [`get_local_atmosphere`] and
/// [`get_local_atmosphere_humid`]. Returns `(local_temp_k, local_pressure_pa, density_kg_m3)`.
#[inline]
fn local_temp_pressure_density(
    altitude_m: f64,
    base_alt: f64,
    base_temp_c: f64,
    base_press_hpa: f64,
    base_ratio: f64,
) -> (f64, f64, f64) {
    // Round altitude to the nearest meter for caching in Python
    let altitude_m_rounded = altitude_m.round();
    let base_temp_k = base_temp_c + 273.15;

    // A non-finite endpoint would make the boundary walk fail to advance. The
    // previous single-column formula also produced non-finite outputs here.
    if !altitude_m_rounded.is_finite() || !base_alt.is_finite() {
        return (f64::NAN, f64::NAN, f64::NAN);
    }

    let (temp_k, pressure_hpa) = integrate_local_atmosphere_layers(
        base_alt,
        altitude_m_rounded,
        base_temp_k,
        base_press_hpa,
    );

    // Enhanced density calculation
    let density_ratio = base_ratio * (base_temp_k * pressure_hpa) / (base_press_hpa * temp_k);
    let density = density_ratio * 1.225;

    (temp_k, pressure_hpa * 100.0, density)
}

/// Carry arbitrary station temperature and pressure through each crossed ICAO
/// layer. The layer table supplies lapse rates and boundaries only; station
/// deviations from the standard atmosphere remain anchored at `base_alt`.
fn integrate_local_atmosphere_layers(
    base_alt: f64,
    target_alt: f64,
    mut temp_k: f64,
    mut pressure_hpa: f64,
) -> (f64, f64) {
    let mut current_alt = base_alt;

    if target_alt > current_alt {
        while current_alt < target_alt {
            // At an exact boundary, ascent starts in the higher layer.
            let layer_index = ICAO_LAYERS
                .iter()
                .rposition(|layer| current_alt >= layer.base_altitude)
                .unwrap_or(0);
            let segment_end = ICAO_LAYERS
                .get(layer_index + 1)
                .map_or(target_alt, |next| target_alt.min(next.base_altitude));
            (temp_k, pressure_hpa) = integrate_local_atmosphere_segment(
                temp_k,
                pressure_hpa,
                segment_end - current_alt,
                ICAO_LAYERS[layer_index].lapse_rate,
            );
            current_alt = segment_end;
        }
    } else {
        while current_alt > target_alt {
            // At an exact boundary, descent starts in the lower layer. The
            // strict comparison is what makes a cross-layer round trip reversible.
            let layer_index = ICAO_LAYERS
                .iter()
                .rposition(|layer| current_alt > layer.base_altitude)
                .unwrap_or(0);
            let segment_end = if layer_index == 0 {
                target_alt
            } else {
                target_alt.max(ICAO_LAYERS[layer_index].base_altitude)
            };
            (temp_k, pressure_hpa) = integrate_local_atmosphere_segment(
                temp_k,
                pressure_hpa,
                segment_end - current_alt,
                ICAO_LAYERS[layer_index].lapse_rate,
            );
            current_alt = segment_end;
        }
    }

    (temp_k, pressure_hpa)
}

#[inline]
fn integrate_local_atmosphere_segment(
    base_temp_k: f64,
    base_pressure_hpa: f64,
    height_diff: f64,
    lapse_rate: f64,
) -> (f64, f64) {
    let temp_k = base_temp_k + lapse_rate * height_diff;
    let pressure_hpa = if lapse_rate.abs() < 1e-10 {
        base_pressure_hpa * (-G_ACCEL_MPS2 * height_diff / (R_AIR * base_temp_k)).exp()
    } else {
        let temp_ratio = temp_k / base_temp_k;
        base_pressure_hpa * temp_ratio.powf(-G_ACCEL_MPS2 / (lapse_rate * R_AIR))
    };

    (temp_k, pressure_hpa)
}

/// Determine local lapse rate based on altitude and atmospheric layer.
#[cfg(test)]
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

/// A single downrange-referenced atmosphere zone:
/// `(temp_c, pressure_hpa, humidity_percent, until_distance_m)`.
///
/// The T/P/H are the STATION-REFERENCED conditions (defined at the shooter base altitude) that
/// apply from the previous segment's threshold out to `until_distance_m`. This mirrors
/// [`crate::wind::WindSegment`]'s `(speed, angle, until_distance)` shape so the two segmented
/// models compose the same way (wind by X, atmosphere by X, altitude lapse by Y).
pub type AtmoSegment = (f64, f64, f64, f64);

/// Downrange-segmented atmosphere handler (MBA-1137), the density analogue of
/// [`crate::wind::WindSock`].
///
/// Holds a set of station-referenced atmosphere zones ordered by their `until_distance_m`
/// threshold and answers a stateless downrange lookup ([`AtmoSock::atmo_for_range`]).
///
/// The zone T/P/H are the base (shooter-altitude) conditions for that stretch of range; the
/// solver swaps them into the SAME `get_local_atmosphere` altitude-lapse pipeline that a
/// single-station solve uses, so the downrange (X) zone and the vertical (Y) altitude lapse
/// compose orthogonally without double-counting (the zone sets the base tuple, the lapse
/// multiplies on top of it).
#[derive(Debug, Clone)]
pub struct AtmoSock {
    /// Zones sorted ascending by `until_distance_m` (segment slot 3).
    segments: Vec<AtmoSegment>,
}

impl AtmoSock {
    /// Create a new `AtmoSock` from station-referenced atmosphere zones.
    ///
    /// Each segment is `(temp_c, pressure_hpa, humidity_percent, until_distance_m)`. Segments are
    /// sorted by `until_distance_m` (NaN thresholds are ordered last, matching `WindSock::new`).
    pub fn new(mut segments: Vec<AtmoSegment>) -> Self {
        // Sort by until_distance, treating NaN as greater than any value (mirrors WindSock::new).
        segments.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Greater));
        AtmoSock { segments }
    }

    /// True when this sock carries no zones (a lookup falls back to sea-level ISA).
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    /// Stateless downrange lookup of the active zone's `(temp_c, pressure_hpa, humidity_percent)`.
    ///
    /// Selection matches [`crate::wind::WindSock::vector_for_range_stateless`]: the first segment
    /// whose `until_distance_m` STRICTLY exceeds `downrange_m` wins (thresholds are upper-exclusive).
    /// Unlike wind — which returns zero past the last threshold — the LAST zone is used for any
    /// distance at or beyond the final threshold (there is no "zero atmosphere"). An empty sock
    /// returns the sea-level ISA reference `(15 C, 1013.25 hPa, 0% RH)`.
    ///
    /// This is stateless and safe for numerical integration (the same X may be queried repeatedly
    /// or out of order across RK substeps).
    pub fn atmo_for_range(&self, downrange_m: f64) -> (f64, f64, f64) {
        if self.segments.is_empty() {
            return (15.0, 1013.25, 0.0); // sea-level ISA fallback
        }
        // NaN X can't be ordered; use the first (nearest) zone deterministically.
        if downrange_m.is_nan() {
            let s = self.segments[0];
            return (s.0, s.1, s.2);
        }
        for seg in &self.segments {
            if downrange_m < seg.3 {
                return (seg.0, seg.1, seg.2);
            }
        }
        // Beyond the final threshold: hold the last zone (no zeroing).
        let last = self.segments[self.segments.len() - 1];
        (last.0, last.1, last.2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- MBA-1136: CIPM-2007 as the single canonical humid-air density ----

    /// Gate 1: dry sea-level (15 C, 1013.25 hPa, 0% RH) must stay at the ISA reference — density
    /// 1.225 +- 0.002 kg/m^3 and speed of sound 340.3 +- 0.6 m/s. CIPM-2007 at 0% RH reduces to
    /// dry-air ideal gas to within rounding, so this is essentially unchanged from the pre-CIPM
    /// baseline (baseline was 1.225012 / 340.294; now 1.225521 / 340.294 — the tiny density bump
    /// is CIPM compressibility + exact molar mass, the speed of sound is bit-identical).
    #[test]
    fn test_mba1136_dry_sea_level_reference() {
        let (density, sos) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 0.0);
        assert!(
            (density - 1.225).abs() < 0.002,
            "dry sea-level density {density} not within 1.225 +- 0.002"
        );
        assert!(
            (sos - 340.3).abs() < 0.6,
            "dry sea-level speed of sound {sos} not within 340.3 +- 0.6"
        );
    }

    /// Gate 2: humid air (15 C, 1013.25 hPa, 50% RH) is the CIPM-2007 value (~1.2211 +- 0.002),
    /// STRICTLY lighter than dry air at the same T/P, with a speed of sound slightly ABOVE dry.
    #[test]
    fn test_mba1136_humid_lighter_than_dry() {
        let (dry_rho, dry_sos) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 0.0);
        let (moist_rho, moist_sos) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 50.0);

        assert!(
            (moist_rho - 1.2211).abs() < 0.002,
            "50% RH density {moist_rho} not within CIPM 1.2211 +- 0.002"
        );
        assert!(
            moist_rho < dry_rho,
            "moist air ({moist_rho}) must be lighter than dry ({dry_rho})"
        );
        assert!(
            moist_sos > dry_sos,
            "moist speed of sound ({moist_sos}) must exceed dry ({dry_sos})"
        );
    }

    /// Gate 3: density is monotone-decreasing in humidity (100% RH < 50% RH < 0% RH).
    #[test]
    fn test_mba1136_density_monotone_in_humidity() {
        let (rho_0, _) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 0.0);
        let (rho_50, _) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 50.0);
        let (rho_100, _) = calculate_atmosphere(0.0, Some(15.0), Some(1013.25), 100.0);
        assert!(
            rho_100 < rho_50 && rho_50 < rho_0,
            "humidity monotonicity violated: 100%={rho_100}, 50%={rho_50}, 0%={rho_0}"
        );
    }

    /// rank 28: `calculate_atmosphere`'s density is exactly `calculate_air_density_cimp` — there
    /// is a single canonical humid-air density path (no separate Arden-Buck ideal-gas density).
    #[test]
    fn test_mba1136_atmosphere_density_is_cipm() {
        for (t, p, h) in [
            (15.0, 1013.25, 0.0),
            (15.0, 1013.25, 50.0),
            (30.0, 1000.0, 80.0),
            (-10.0, 1020.0, 20.0),
        ] {
            let (density, _) = calculate_atmosphere(0.0, Some(t), Some(p), h);
            let cipm = calculate_air_density_cimp(t, p, h);
            assert_eq!(
                density, cipm,
                "calculate_atmosphere density must equal CIPM at {t}C/{p}hPa/{h}%"
            );
        }
    }

    /// rank 9: the extracted `moist_speed_of_sound` is exactly what `calculate_atmosphere`
    /// returns (behavior-identical extraction), across dry and humid conditions.
    #[test]
    fn test_mba1136_moist_speed_of_sound_extraction() {
        for (t, p, h) in [
            (15.0, 1013.25, 0.0),
            (15.0, 1013.25, 50.0),
            (25.0, 900.0, 100.0),
        ] {
            let (_, sos) = calculate_atmosphere(0.0, Some(t), Some(p), h);
            let extracted = moist_speed_of_sound(t + 273.15, p * 100.0, h);
            assert_eq!(
                sos, extracted,
                "extracted moist_speed_of_sound must match calculate_atmosphere at {t}C/{p}hPa/{h}%"
            );
        }
    }

    /// rank 9: `get_local_atmosphere` is behavior-locked after the shared-helper refactor.
    /// Reference values captured from the pre-refactor implementation
    /// (base: 500 m, 10 C, 950 hPa, ratio 1.05).
    #[test]
    fn test_mba1136_get_local_atmosphere_unchanged() {
        let (d0, c0) = get_local_atmosphere(500.0, 500.0, 10.0, 950.0, 1.05);
        assert!((d0 - 1.286250000000).abs() < 1e-9, "local density@500m drifted: {d0}");
        assert!((c0 - 337.328657395129).abs() < 1e-9, "local sos@500m drifted: {c0}");
        let (d1, c1) = get_local_atmosphere(1500.0, 500.0, 10.0, 950.0, 1.05);
        assert!((d1 - 1.165201643681).abs() < 1e-9, "local density@1500m drifted: {d1}");
        assert!((c1 - 333.434314520866).abs() < 1e-9, "local sos@1500m drifted: {c1}");
    }

    /// rank 9: `get_local_atmosphere_humid` returns the SAME density as `get_local_atmosphere`,
    /// and at 0% RH its speed of sound reduces to the dry value (within the 401.874-vs-gamma*R
    /// constant rounding). At real humidity the speed of sound exceeds the dry value.
    #[test]
    fn test_mba1136_get_local_atmosphere_humid() {
        let (d_dry, c_dry) = get_local_atmosphere(1500.0, 500.0, 10.0, 950.0, 1.05);
        let (d_h0, c_h0) = get_local_atmosphere_humid(1500.0, 500.0, 10.0, 950.0, 1.05, 0.0);
        assert_eq!(d_dry, d_h0, "humid variant must not change density");
        assert!(
            (c_h0 - c_dry).abs() < 1e-3,
            "0% RH humid sos {c_h0} should match dry sos {c_dry}"
        );
        let (_, c_h80) = get_local_atmosphere_humid(1500.0, 500.0, 10.0, 950.0, 1.05, 80.0);
        assert!(c_h80 > c_dry, "humid sos {c_h80} should exceed dry {c_dry}");
    }

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
    fn local_atmosphere_walks_icao_layers_continuously() {
        for altitude_m in [
            10999.0, 11000.0, 11001.0, 11050.0, 19999.0, 20000.0, 20001.0, 25000.0,
        ] {
            let (local_temp_k, local_pressure_pa, _) =
                local_temp_pressure_density(altitude_m, 0.0, 15.0, 1013.25, 1.0);
            let (standard_temp_k, standard_pressure_pa) =
                calculate_icao_standard_atmosphere(altitude_m);

            assert!(
                (local_temp_k - standard_temp_k).abs() < 1e-9,
                "local temperature diverged from ICAO at {altitude_m} m: local={local_temp_k}, standard={standard_temp_k}"
            );
            assert!(
                ((local_pressure_pa - standard_pressure_pa) / standard_pressure_pa).abs() < 5e-5,
                "local pressure diverged from ICAO at {altitude_m} m: local={local_pressure_pa}, standard={standard_pressure_pa}"
            );
        }

        let (density_below, sound_below) =
            get_local_atmosphere(10999.0, 0.0, 15.0, 1013.25, 1.0);
        let (density_above, sound_above) =
            get_local_atmosphere(11001.0, 0.0, 15.0, 1013.25, 1.0);
        assert!(
            (density_below - density_above).abs() < 0.001,
            "density jumped across 11 km: below={density_below}, above={density_above}"
        );
        assert!(
            (sound_below - sound_above).abs() < 0.1,
            "speed of sound jumped across 11 km: below={sound_below}, above={sound_above}"
        );
    }

    #[test]
    fn local_atmosphere_preserves_nonstandard_station_offset_and_round_trips() {
        let base_alt = 7500.0;
        let base_temp_c = 5.0;
        let base_pressure_hpa = 410.0;
        let base_ratio = 0.72;

        let (high_temp_k, high_pressure_pa, high_density) = local_temp_pressure_density(
            25000.0,
            base_alt,
            base_temp_c,
            base_pressure_hpa,
            base_ratio,
        );
        assert!((high_temp_k - 260.4).abs() < 1e-9);
        assert!((high_pressure_pa / 100.0 - 40.505969).abs() < 1e-6);
        assert!((high_density - 0.093076885).abs() < 1e-9);

        let (back_temp_k, back_pressure_pa, back_density) = local_temp_pressure_density(
            base_alt,
            25000.0,
            high_temp_k - 273.15,
            high_pressure_pa / 100.0,
            high_density / 1.225,
        );
        assert!((back_temp_k - (base_temp_c + 273.15)).abs() < 1e-9);
        assert!((back_pressure_pa / 100.0 - base_pressure_hpa).abs() < 1e-8);
        assert!((back_density - base_ratio * 1.225).abs() < 1e-9);
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

    // ---- MBA-1137: AtmoSock stateless downrange lookup (mirrors the WindSock tests) ----

    #[test]
    fn test_atmo_sock_empty_falls_back_to_isa() {
        let sock = AtmoSock::new(vec![]);
        assert!(sock.is_empty());
        // Empty sock returns the sea-level ISA reference regardless of distance.
        assert_eq!(sock.atmo_for_range(0.0), (15.0, 1013.25, 0.0));
        assert_eq!(sock.atmo_for_range(500.0), (15.0, 1013.25, 0.0));
    }

    #[test]
    fn test_atmo_sock_single_segment_holds_beyond_last() {
        // One zone until 100 m; it must apply BOTH before and beyond the threshold (unlike wind,
        // which zeroes past the last segment — atmosphere holds the last zone).
        let sock = AtmoSock::new(vec![(25.0, 1000.0, 30.0, 100.0)]);
        assert_eq!(sock.atmo_for_range(50.0), (25.0, 1000.0, 30.0));
        assert_eq!(sock.atmo_for_range(100.0), (25.0, 1000.0, 30.0)); // beyond last -> hold
        assert_eq!(sock.atmo_for_range(5000.0), (25.0, 1000.0, 30.0));
    }

    #[test]
    fn test_atmo_sock_boundary_is_upper_exclusive() {
        // A zone's until_distance_m is exclusive: a query exactly at the boundary rolls to the
        // next zone (mirrors WindSock::test_wind_sock_boundary_is_upper_exclusive).
        let sock = AtmoSock::new(vec![
            (30.0, 1010.0, 80.0, 100.0), // hot/humid near zone
            (-5.0, 900.0, 10.0, 200.0),  // cold/thin far zone
        ]);
        // Just below 100 m -> first zone.
        assert_eq!(sock.atmo_for_range(99.999), (30.0, 1010.0, 80.0));
        // Exactly 100 m -> second zone.
        assert_eq!(sock.atmo_for_range(100.0), (-5.0, 900.0, 10.0));
        // Beyond the last boundary -> hold the last zone (NOT zeroed).
        assert_eq!(sock.atmo_for_range(200.0), (-5.0, 900.0, 10.0));
        assert_eq!(sock.atmo_for_range(1e6), (-5.0, 900.0, 10.0));
    }

    #[test]
    fn test_atmo_sock_sorts_unordered_segments() {
        // Segments supplied out of order are sorted by until_distance so the lookup is monotone.
        let sock = AtmoSock::new(vec![
            (-5.0, 900.0, 10.0, 200.0),
            (30.0, 1010.0, 80.0, 100.0),
        ]);
        assert_eq!(sock.atmo_for_range(50.0), (30.0, 1010.0, 80.0));
        assert_eq!(sock.atmo_for_range(150.0), (-5.0, 900.0, 10.0));
    }

    #[test]
    fn test_atmo_sock_nan_uses_first_zone() {
        let sock = AtmoSock::new(vec![
            (30.0, 1010.0, 80.0, 100.0),
            (-5.0, 900.0, 10.0, 200.0),
        ]);
        // NaN can't be ordered; deterministically use the nearest (first) zone rather than panic.
        assert_eq!(sock.atmo_for_range(f64::NAN), (30.0, 1010.0, 80.0));
    }
}
