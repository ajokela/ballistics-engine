//! Altitude-dependent wind shear modeling for ballistics.
//!
//! Wind shear refers to the change in wind speed and/or direction with altitude.
//! This is important for long-range ballistics where projectiles reach significant
//! altitudes and experience different wind conditions at different heights.

// Wind shear modeling - now integrated!

use nalgebra::Vector3;
use std::f64::consts::PI;

/// Wind shear model types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WindShearModel {
    None,
    Logarithmic,
    PowerLaw,
    EkmanSpiral,
    CustomLayers,
}

/// Wind conditions at a specific altitude
#[derive(Debug, Clone, Copy)]
pub struct WindLayer {
    pub altitude_m: f64,
    pub speed_mps: f64,
    pub direction_deg: f64, // Direction wind is coming FROM
}

impl WindLayer {
    /// Convert to wind vector [x, y, z] in m/s
    /// STANDARD BALLISTICS CONVENTION: X=downrange, Y=vertical, Z=lateral
    pub fn to_vector(&self) -> Vector3<f64> {
        let ang = self.direction_deg.to_radians();
        Vector3::new(
            -self.speed_mps * ang.cos(), // X (downrange - head/tail component)
            0.0,                         // Y (vertical)
            -self.speed_mps * ang.sin(), // Z (lateral - crosswind component)
        )
    }
}

/// Complete wind shear profile definition
#[derive(Debug, Clone)]
pub struct WindShearProfile {
    pub model: WindShearModel,
    pub surface_wind: WindLayer,
    pub reference_height: f64, // Standard meteorological measurement height
    pub roughness_length: f64, // Surface roughness (0.03 = short grass)
    pub power_exponent: f64,   // Power law exponent (1/7 for neutral stability)
    pub geostrophic_wind: Option<WindLayer>, // Wind above boundary layer
    pub custom_layers: Vec<WindLayer>,
}

impl Default for WindShearProfile {
    fn default() -> Self {
        Self {
            model: WindShearModel::None,
            surface_wind: WindLayer {
                altitude_m: 0.0,
                speed_mps: 0.0,
                direction_deg: 0.0,
            },
            reference_height: 10.0,
            roughness_length: 0.03,
            power_exponent: 1.0 / 7.0,
            geostrophic_wind: None,
            custom_layers: Vec::new(),
        }
    }
}

impl WindShearProfile {
    /// Get wind vector at specified altitude
    pub fn get_wind_at_altitude(&self, altitude_m: f64) -> Vector3<f64> {
        match self.model {
            WindShearModel::None => self.surface_wind.to_vector(),
            WindShearModel::Logarithmic => self.logarithmic_profile(altitude_m),
            WindShearModel::PowerLaw => self.power_law_profile(altitude_m),
            WindShearModel::EkmanSpiral => self.ekman_spiral_profile(altitude_m),
            WindShearModel::CustomLayers => self.interpolate_layers(altitude_m),
        }
    }

    /// Logarithmic wind profile (boundary layer)
    /// U(z) = U_ref * ln(z/z0) / ln(z_ref/z0)
    fn logarithmic_profile(&self, altitude_m: f64) -> Vector3<f64> {
        // Handle negative altitudes (bullet below sight line)
        // Use absolute altitude, but add small offset only if very close to ground
        let effective_altitude = if altitude_m < 0.0 {
            // For negative altitudes, use a small positive value
            0.001 // 1mm above ground
        } else if altitude_m < 0.001 {
            // Very small positive altitudes
            0.001
        } else {
            altitude_m
        };

        // If very close to roughness length, return near-zero wind
        if effective_altitude <= self.roughness_length {
            return Vector3::zeros();
        }

        // Calculate speed ratio
        let speed_ratio = if effective_altitude > self.roughness_length
            && self.reference_height > self.roughness_length
        {
            (effective_altitude / self.roughness_length).ln()
                / (self.reference_height / self.roughness_length).ln()
        } else {
            1.0
        };

        // Apply to surface wind
        self.surface_wind.to_vector() * speed_ratio.max(0.0)
    }

    /// Power law wind profile
    fn power_law_profile(&self, altitude_m: f64) -> Vector3<f64> {
        if altitude_m <= 0.0 {
            return Vector3::zeros();
        }

        // Calculate speed ratio
        let speed_ratio = (altitude_m / self.reference_height).powf(self.power_exponent);

        // Apply to surface wind
        self.surface_wind.to_vector() * speed_ratio
    }

    /// Ekman spiral - wind direction changes with altitude
    fn ekman_spiral_profile(&self, altitude_m: f64) -> Vector3<f64> {
        // Default geostrophic wind if not specified
        let geo_wind = self.geostrophic_wind.unwrap_or({
            WindLayer {
                altitude_m: 1000.0,
                speed_mps: self.surface_wind.speed_mps * 1.5,
                direction_deg: self.surface_wind.direction_deg - 30.0, // 30° backing
            }
        });

        // Ekman layer depth (simplified)
        let ekman_depth = 1000.0; // meters

        if altitude_m >= ekman_depth {
            return geo_wind.to_vector();
        }

        // Interpolate between surface and geostrophic
        let ratio = altitude_m / ekman_depth;

        // Interpolate speed
        let speed = self.surface_wind.speed_mps * (1.0 - ratio) + geo_wind.speed_mps * ratio;

        // Interpolate direction (accounting for circular interpolation)
        let dir1 = self.surface_wind.direction_deg.to_radians();
        let mut dir2 = geo_wind.direction_deg.to_radians();

        // Handle angle wrapping
        if (dir2 - dir1).abs() > PI {
            if dir2 > dir1 {
                dir2 -= 2.0 * PI;
            } else {
                dir2 += 2.0 * PI;
            }
        }

        let direction_rad = dir1 * (1.0 - ratio) + dir2 * ratio;

        // Convert to vector (X=downrange, Y=vertical, Z=lateral)
        Vector3::new(
            -speed * direction_rad.cos(), // X (downrange - head/tail)
            0.0,
            -speed * direction_rad.sin(), // Z (lateral - crosswind)
        )
    }

    /// Interpolate wind from custom altitude layers
    fn interpolate_layers(&self, altitude_m: f64) -> Vector3<f64> {
        if self.custom_layers.is_empty() {
            return self.surface_wind.to_vector();
        }

        // Find bracketing layers
        let mut low_idx = 0;
        let mut high_idx = 0;

        for (i, layer) in self.custom_layers.iter().enumerate() {
            if layer.altitude_m <= altitude_m {
                low_idx = i;
            }
            if layer.altitude_m >= altitude_m {
                high_idx = i;
                break;
            }
        }

        // Same layer or out of bounds
        if low_idx == high_idx {
            return self.custom_layers[low_idx].to_vector();
        }

        // Linear interpolation
        let low_layer = &self.custom_layers[low_idx];
        let high_layer = &self.custom_layers[high_idx];

        // MBA-246: Guard against division by zero when layers have same altitude
        let altitude_diff = high_layer.altitude_m - low_layer.altitude_m;
        if altitude_diff.abs() < 1e-9 {
            return low_layer.to_vector();
        }

        let ratio = (altitude_m - low_layer.altitude_m) / altitude_diff;

        // Interpolate vectors
        let low_vec = low_layer.to_vector();
        let high_vec = high_layer.to_vector();

        low_vec * (1.0 - ratio) + high_vec * ratio
    }
}

/// Extended wind sock with altitude-dependent shear
#[derive(Debug, Clone)]
pub struct WindShearWindSock {
    pub segments: Vec<(f64, f64, f64)>, // (speed_mps, angle_deg, until_range_m)
    pub shear_profile: Option<WindShearProfile>,
    pub shooter_altitude_m: f64,
}

impl WindShearWindSock {
    pub fn new(segments: Vec<(f64, f64, f64)>, shear_profile: Option<WindShearProfile>) -> Self {
        Self {
            segments,
            shear_profile,
            shooter_altitude_m: 0.0,
        }
    }

    pub fn with_shooter_altitude(
        segments: Vec<(f64, f64, f64)>,
        shear_profile: Option<WindShearProfile>,
        shooter_altitude_m: f64,
    ) -> Self {
        Self {
            segments,
            shear_profile,
            shooter_altitude_m,
        }
    }

    /// Get wind vector for 3D position
    /// Standard ballistics coordinate system: X=downrange, Y=vertical, Z=lateral
    pub fn vector_for_position(&self, position: Vector3<f64>) -> Vector3<f64> {
        let range_m = position.x; // X is downrange (McCoy)
        let altitude_m = position.y; // Relative to shooter

        // Get base wind at this range
        let base_wind = self.get_range_wind(range_m);

        if let Some(profile) = &self.shear_profile {
            if profile.model != WindShearModel::None {
                // Apply altitude-dependent shear
                // Calculate absolute altitude by adding shooter's altitude
                let absolute_altitude_m = altitude_m + self.shooter_altitude_m;
                let altitude_vec = profile.get_wind_at_altitude(absolute_altitude_m);

                // Scale the base wind by altitude profile
                if base_wind.norm() > 0.0 {
                    let scale_factor =
                        altitude_vec.norm() / profile.surface_wind.speed_mps.max(1e-9);
                    return base_wind * scale_factor;
                }

                return altitude_vec;
            }
        }

        base_wind
    }

    /// Get wind based on horizontal range
    /// Returns wind vector in standard ballistics coordinates: X=downrange, Y=vertical, Z=lateral
    fn get_range_wind(&self, range_m: f64) -> Vector3<f64> {
        if range_m.is_nan() || self.segments.is_empty() {
            return Vector3::zeros();
        }

        // Find appropriate wind segment
        for &(speed_mps, angle_deg, until_dist) in &self.segments {
            if range_m <= until_dist {
                let ang = angle_deg.to_radians();
                return Vector3::new(
                    -speed_mps * ang.cos(), // X (downrange - head/tail)
                    0.0,
                    -speed_mps * ang.sin(), // Z (lateral - crosswind)
                );
            }
        }

        // Beyond all segments
        Vector3::zeros()
    }
}

/// High-level API function to get wind at arbitrary position
///
/// This is a convenience wrapper that handles wind segments, shear models,
/// and altitude calculations in a single function call.
///
/// # Arguments
/// * `position` - 3D position vector [x_downrange, y_vertical, z_lateral]
/// * `wind_segments` - Wind segments as (speed_kmh, angle_deg, until_distance_m)
/// * `enable_wind_shear` - Whether to apply wind shear modeling
/// * `wind_shear_model` - Model type: "none", "logarithmic", "power_law", "ekman_spiral"
/// * `shooter_altitude_m` - Shooter's altitude above sea level
///
/// # Returns
/// Wind vector in m/s [x_downrange, y_vertical, z_lateral]
pub fn get_wind_at_position(
    position: &Vector3<f64>,
    wind_segments: &[(f64, f64, f64)], // (speed_kmh, angle_deg, until_distance_m)
    enable_wind_shear: bool,
    wind_shear_model: &str,
    shooter_altitude_m: f64,
) -> Vector3<f64> {
    // X IS DOWNRANGE (McCoy)
    let range_m = position[0];
    let altitude_m = position[1]; // Y is vertical, relative to shooter

    // Find appropriate wind segment based on range
    let base_wind = if wind_segments.is_empty() {
        (0.0, 0.0)
    } else {
        // Find the segment that covers this range
        let mut found_wind = (wind_segments[0].0, wind_segments[0].1);
        for seg in wind_segments {
            if range_m <= seg.2 {
                found_wind = (seg.0, seg.1);
                break;
            }
        }
        found_wind
    };

    // Convert base wind from km/h to m/s
    let base_speed_mps = base_wind.0 * 0.2777778; // km/h to m/s
    let base_direction_deg = base_wind.1;

    if !enable_wind_shear || wind_shear_model == "none" {
        // No shear - return constant wind
        let ang = base_direction_deg.to_radians();
        return Vector3::new(
            -base_speed_mps * ang.cos(), // x (downrange)
            0.0,                         // y (vertical)
            -base_speed_mps * ang.sin(), // z (lateral)
        );
    }

    // Create wind shear profile
    let mut profile = WindShearProfile::default();
    profile.model = match wind_shear_model {
        "logarithmic" => WindShearModel::Logarithmic,
        "power_law" | "powerlaw" => WindShearModel::PowerLaw,
        "ekman_spiral" | "ekman" => WindShearModel::EkmanSpiral,
        "custom_layers" | "custom" => WindShearModel::CustomLayers,
        _ => WindShearModel::None,
    };
    profile.surface_wind = WindLayer {
        altitude_m: 0.0,
        speed_mps: base_speed_mps,
        direction_deg: base_direction_deg,
    };

    // Calculate absolute altitude by adding shooter's altitude
    let absolute_altitude_m = altitude_m + shooter_altitude_m;

    // OPTIMIZATION: Skip complex shear for very small altitude changes
    // This avoids numerical issues near ground level
    if absolute_altitude_m.abs() < 0.1 && altitude_m.abs() < 0.1 {
        // Near ground level - use base wind directly with reduction
        let ang = base_direction_deg.to_radians();
        return Vector3::new(
            -base_speed_mps * ang.cos() * 0.5, // Reduced at ground
            0.0,
            -base_speed_mps * ang.sin() * 0.5,
        );
    }

    // OPTIMIZATION: For long-range shots, use simplified model
    // to avoid numerical instability in RK45 integration
    if range_m > 800.0 {
        // Use simplified linear interpolation for stability
        let altitude_factor = (1.0 + absolute_altitude_m / 100.0).min(2.0).max(0.1);
        let sheared_speed = base_speed_mps * altitude_factor;
        let ang = base_direction_deg.to_radians();
        return Vector3::new(-sheared_speed * ang.cos(), 0.0, -sheared_speed * ang.sin());
    }

    // For normal ranges, use full shear model with clamped altitude
    let clamped_altitude = absolute_altitude_m.max(-10.0).min(1000.0);
    profile.get_wind_at_altitude(clamped_altitude)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wind_layer() {
        // Standard ballistics coordinate system: X=downrange, Y=vertical, Z=lateral
        // Wind direction: 0°=headwind, 90°=from right, 180°=tailwind, 270°=from left

        // Test 0° wind (from north/front - headwind)
        let layer_headwind = WindLayer {
            altitude_m: 100.0,
            speed_mps: 10.0,
            direction_deg: 0.0, // Wind from front (headwind)
        };

        let vec = layer_headwind.to_vector();
        assert!(
            (vec.x - (-10.0)).abs() < 1e-6,
            "0° wind should be headwind (negative X downrange)"
        );
        assert_eq!(vec.y, 0.0);
        assert!(
            (vec.z).abs() < 1e-6,
            "0° wind (headwind) should have zero lateral (Z) component"
        );

        // Test 90° wind (from right)
        let layer_crosswind = WindLayer {
            altitude_m: 100.0,
            speed_mps: 10.0,
            direction_deg: 90.0, // Wind from right
        };

        let vec_cross = layer_crosswind.to_vector();
        assert!(
            (vec_cross.z - (-10.0)).abs() < 1e-6,
            "90° wind should have negative Z lateral (from right)"
        );
        assert_eq!(vec_cross.y, 0.0);
        assert!(
            (vec_cross.x).abs() < 1e-6,
            "90° wind (crosswind) should have zero downrange (X) component"
        );
    }

    #[test]
    fn test_logarithmic_profile() {
        let mut profile = WindShearProfile::default();
        profile.model = WindShearModel::Logarithmic;
        profile.surface_wind = WindLayer {
            altitude_m: 0.0,
            speed_mps: 10.0,
            direction_deg: 0.0,
        };

        // Wind should increase with altitude
        let v10 = profile.get_wind_at_altitude(10.0).norm();
        let v50 = profile.get_wind_at_altitude(50.0).norm();
        let v100 = profile.get_wind_at_altitude(100.0).norm();

        assert!(v10 > 0.0);
        assert!(v50 > v10);
        assert!(v100 > v50);
    }

    #[test]
    fn test_power_law_profile() {
        let mut profile = WindShearProfile::default();
        profile.model = WindShearModel::PowerLaw;
        profile.surface_wind = WindLayer {
            altitude_m: 0.0,
            speed_mps: 10.0,
            direction_deg: 0.0,
        };

        // Check power law relationship
        let v100 = profile.get_wind_at_altitude(100.0).norm();
        let expected = 10.0 * (100.0_f64 / 10.0).powf(1.0 / 7.0);
        assert!((v100 - expected).abs() < 0.01);
    }
}
