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
    pub fn to_vector(&self) -> Vector3<f64> {
        let ang = self.direction_deg.to_radians();
        Vector3::new(
            -self.speed_mps * ang.cos(), // x (downrange)
            0.0,                         // y (vertical)
            -self.speed_mps * ang.sin(), // z (lateral)
        )
    }
}

/// Complete wind shear profile definition
#[derive(Debug, Clone)]
pub struct WindShearProfile {
    pub model: WindShearModel,
    pub surface_wind: WindLayer,
    pub reference_height: f64,    // Standard meteorological measurement height
    pub roughness_length: f64,    // Surface roughness (0.03 = short grass)
    pub power_exponent: f64,      // Power law exponent (1/7 for neutral stability)
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
    fn logarithmic_profile(&self, altitude_m: f64) -> Vector3<f64> {
        if altitude_m <= self.roughness_length {
            return Vector3::zeros();
        }

        // Calculate speed ratio
        let speed_ratio = if altitude_m > self.roughness_length && self.reference_height > self.roughness_length {
            (altitude_m / self.roughness_length).ln() / (self.reference_height / self.roughness_length).ln()
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

        // Convert to vector
        Vector3::new(
            -speed * direction_rad.cos(),
            0.0,
            -speed * direction_rad.sin(),
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

        let ratio = (altitude_m - low_layer.altitude_m) / (high_layer.altitude_m - low_layer.altitude_m);

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
    
    pub fn with_shooter_altitude(segments: Vec<(f64, f64, f64)>, shear_profile: Option<WindShearProfile>, shooter_altitude_m: f64) -> Self {
        Self {
            segments,
            shear_profile,
            shooter_altitude_m,
        }
    }

    /// Get wind vector for 3D position
    pub fn vector_for_position(&self, position: Vector3<f64>) -> Vector3<f64> {
        let range_m = position.x;
        let altitude_m = position.y;  // Relative to shooter

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
                    let scale_factor = altitude_vec.norm() / profile.surface_wind.speed_mps.max(1e-9);
                    return base_wind * scale_factor;
                }

                return altitude_vec;
            }
        }

        base_wind
    }

    /// Get wind based on horizontal range
    fn get_range_wind(&self, range_m: f64) -> Vector3<f64> {
        if range_m.is_nan() || self.segments.is_empty() {
            return Vector3::zeros();
        }

        // Find appropriate wind segment
        for &(speed_mps, angle_deg, until_dist) in &self.segments {
            if range_m <= until_dist {
                let ang = angle_deg.to_radians();
                return Vector3::new(
                    -speed_mps * ang.cos(),
                    0.0,
                    -speed_mps * ang.sin(),
                );
            }
        }

        // Beyond all segments
        Vector3::zeros()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wind_layer() {
        let layer = WindLayer {
            altitude_m: 100.0,
            speed_mps: 10.0,
            direction_deg: 0.0, // North wind
        };

        let vec = layer.to_vector();
        assert!((vec.x - (-10.0)).abs() < 1e-6);
        assert_eq!(vec.y, 0.0);
        assert!((vec.z).abs() < 1e-6);
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