//! HTTP client for Flask API communication
//!
//! This module provides an HTTP client for routing trajectory calculations
//! through the Flask API instead of local computation, giving CLI users
//! access to ML-enhanced predictions.

use serde::{Deserialize, Serialize};
use std::time::Duration;

/// Request structure for trajectory calculation via Flask API
#[derive(Debug, Clone, Serialize)]
pub struct TrajectoryRequest {
    /// Ballistic coefficient value
    pub bc_value: f64,
    /// BC type: "G1" or "G7"
    pub bc_type: String,
    /// Bullet mass in grams
    pub bullet_mass: f64,
    /// Muzzle velocity in m/s
    pub muzzle_velocity: f64,
    /// Target distance in meters
    pub target_distance: f64,
    /// Zero range in meters (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zero_range: Option<f64>,
    /// Wind speed in m/s (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_speed: Option<f64>,
    /// Wind angle in degrees (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_angle: Option<f64>,
    /// Temperature in Celsius (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// Pressure in hPa/mbar (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure: Option<f64>,
    /// Humidity percentage 0-100 (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub humidity: Option<f64>,
    /// Altitude in meters (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub altitude: Option<f64>,
    /// Latitude for Coriolis calculations (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latitude: Option<f64>,
    /// Longitude for weather zones (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub longitude: Option<f64>,
    /// Shot direction/azimuth in degrees (optional, 0=North, 90=East)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shot_direction: Option<f64>,
    /// Shooting angle in degrees (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shooting_angle: Option<f64>,
    /// Barrel twist rate in inches per turn (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub twist_rate: Option<f64>,
    /// Bullet diameter in meters (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bullet_diameter: Option<f64>,
    /// Bullet length in meters (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bullet_length: Option<f64>,
    /// Ground threshold in meters (optional, negative = ignore ground impact)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ground_threshold: Option<f64>,
    /// Enable weather zones (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub enable_weather_zones: Option<bool>,
    /// Enable 3D weather corrections (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub enable_3d_weather: Option<bool>,
    /// Wind shear model (optional: none, logarithmic, power_law, ekman_spiral)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_shear_model: Option<String>,
    /// Weather zone interpolation method (optional: linear, cubic, step)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub weather_zone_interpolation: Option<String>,
    /// Sample/trajectory step interval in meters (optional)
    /// Will be converted to yards for API as trajectory_step
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_interval: Option<f64>,
}

/// Response structure from Flask API trajectory calculation
#[derive(Debug, Clone, Deserialize)]
pub struct TrajectoryResponse {
    /// Array of trajectory points
    pub trajectory: Vec<ApiTrajectoryPoint>,
    /// Zero angle in radians
    pub zero_angle: f64,
    /// Total time of flight in seconds
    pub time_of_flight: f64,
    /// BC confidence score (0-1) if available
    #[serde(default)]
    pub bc_confidence: Option<f64>,
    /// List of ML corrections applied
    #[serde(default)]
    pub ml_corrections_applied: Option<Vec<String>>,
    /// Maximum ordinate (height) in meters
    #[serde(default)]
    pub max_ordinate: Option<f64>,
    /// Impact velocity in m/s
    #[serde(default)]
    pub impact_velocity: Option<f64>,
    /// Impact energy in Joules
    #[serde(default)]
    pub impact_energy: Option<f64>,
}

/// A single point in the trajectory from API response
#[derive(Debug, Clone, Deserialize)]
pub struct ApiTrajectoryPoint {
    /// Range/distance in meters
    pub range: f64,
    /// Drop below line of sight in meters (negative = below)
    pub drop: f64,
    /// Wind drift in meters
    pub drift: f64,
    /// Velocity at this point in m/s
    pub velocity: f64,
    /// Kinetic energy at this point in Joules
    pub energy: f64,
    /// Time of flight to this point in seconds
    pub time: f64,
}

/// Request structure for velocity truing via Flask API
#[derive(Debug, Clone, Serialize)]
pub struct TrueVelocityRequest {
    /// Measured drop in MILs at the target range
    pub measured_drop_mil: f64,
    /// Range at which drop was measured (yards)
    pub range_yd: f64,
    /// Ballistic coefficient
    pub bc: f64,
    /// Drag model ("G1" or "G7")
    pub drag_model: String,
    /// Bullet weight in grains
    pub weight_gr: f64,
    /// Bullet caliber/diameter in inches
    pub caliber: f64,
    /// Zero range in yards (default: 100)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub zero_range_yd: Option<f64>,
    /// Chronograph velocity in fps (optional, for comparison)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chrono_velocity_fps: Option<f64>,
    /// Altitude in feet (default: 0)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub altitude_ft: Option<f64>,
    /// Temperature in Fahrenheit (default: 59)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_f: Option<f64>,
    /// Barometric pressure in inHg (default: 29.92)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure_inhg: Option<f64>,
    /// Humidity percentage (default: 50)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub humidity: Option<f64>,
    /// Sight height in inches (default: 2.0)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sight_height_in: Option<f64>,
    /// Enable BC enhancement (default: true)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub use_bc_enhancement: Option<bool>,
}

/// Response structure from Flask API velocity truing
#[derive(Debug, Clone, Deserialize)]
pub struct TrueVelocityResponse {
    /// Effective muzzle velocity in fps
    pub effective_velocity_fps: f64,
    /// Velocity adjustment from chrono (if chrono provided)
    #[serde(default)]
    pub velocity_adjustment_fps: Option<f64>,
    /// Adjustment as percentage (if chrono provided)
    #[serde(default)]
    pub adjustment_percent: Option<f64>,
    /// Confidence in the result ("high", "medium", "low")
    pub confidence: String,
    /// Number of iterations to converge
    pub iterations: i32,
    /// Final error in MILs after convergence
    pub final_error_mil: f64,
    /// Calculated drop at the effective velocity
    pub calculated_drop_mil: f64,
}

/// Error types for API communication
#[derive(Debug)]
pub enum ApiError {
    /// Network connectivity error
    NetworkError(String),
    /// Request timed out
    Timeout,
    /// Invalid or unparseable response
    InvalidResponse(String),
    /// HTTP error from server (status code, message)
    ServerError(u16, String),
    /// Request building error
    RequestError(String),
}

impl std::fmt::Display for ApiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ApiError::NetworkError(msg) => write!(f, "Network error: {}", msg),
            ApiError::Timeout => write!(f, "Request timed out"),
            ApiError::InvalidResponse(msg) => write!(f, "Invalid response: {}", msg),
            ApiError::ServerError(code, msg) => write!(f, "Server error {}: {}", code, msg),
            ApiError::RequestError(msg) => write!(f, "Request error: {}", msg),
        }
    }
}

impl std::error::Error for ApiError {}

/// HTTP client for Flask API communication
pub struct ApiClient {
    base_url: String,
    timeout: Duration,
}

impl ApiClient {
    /// Create a new API client
    ///
    /// # Arguments
    /// * `base_url` - Base URL of the Flask API (e.g., "https://api.ballistics.7.62x51mm.sh")
    /// * `timeout_secs` - Request timeout in seconds
    pub fn new(base_url: &str, timeout_secs: u64) -> Self {
        // Normalize URL by removing trailing slash
        let base_url = base_url.trim_end_matches('/').to_string();

        Self {
            base_url,
            timeout: Duration::from_secs(timeout_secs),
        }
    }

    /// Calculate trajectory via Flask API
    ///
    /// # Arguments
    /// * `request` - Trajectory calculation request parameters
    ///
    /// # Returns
    /// * `Ok(TrajectoryResponse)` - Successful calculation with trajectory data
    /// * `Err(ApiError)` - Error during API communication
    #[cfg(feature = "online")]
    pub fn calculate_trajectory(
        &self,
        request: &TrajectoryRequest,
    ) -> Result<TrajectoryResponse, ApiError> {
        // Flask API uses GET /v1/calculate with query parameters (imperial units)
        let url = format!("{}/v1/calculate", self.base_url);

        // Convert metric values to imperial for API
        let velocity_fps = request.muzzle_velocity / 0.3048; // m/s to fps
        let mass_grains = request.bullet_mass / 0.0647989; // grams to grains
        let distance_yards = request.target_distance / 0.9144; // meters to yards

        let mut req = ureq::get(&url)
            .set("Accept", "application/json")
            .set("User-Agent", &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")))
            .timeout(self.timeout)
            .query("bc_value", &request.bc_value.to_string())
            .query("bc_type", &request.bc_type)
            .query("bullet_mass", &format!("{:.1}", mass_grains))
            .query("muzzle_velocity", &format!("{:.1}", velocity_fps))
            .query("target_distance", &format!("{:.1}", distance_yards));

        // Add optional parameters
        if let Some(zero_range) = request.zero_range {
            let zero_yards = zero_range / 0.9144;
            req = req.query("zero_distance", &format!("{:.1}", zero_yards));
        }
        if let Some(wind_speed) = request.wind_speed {
            let wind_mph = wind_speed * 2.23694; // m/s to mph
            req = req.query("wind_speed", &format!("{:.1}", wind_mph));
        }
        if let Some(wind_angle) = request.wind_angle {
            req = req.query("wind_angle", &format!("{:.1}", wind_angle));
        }
        if let Some(temp) = request.temperature {
            let temp_f = temp * 9.0 / 5.0 + 32.0; // Celsius to Fahrenheit
            req = req.query("temperature", &format!("{:.1}", temp_f));
        }
        if let Some(pressure) = request.pressure {
            let pressure_inhg = pressure / 33.8639; // hPa to inHg
            req = req.query("pressure", &format!("{:.2}", pressure_inhg));
        }
        if let Some(humidity) = request.humidity {
            req = req.query("humidity", &format!("{:.1}", humidity));
        }
        if let Some(altitude) = request.altitude {
            let altitude_ft = altitude / 0.3048; // meters to feet
            req = req.query("altitude", &format!("{:.1}", altitude_ft));
        }
        if let Some(shooting_angle) = request.shooting_angle {
            req = req.query("shooting_angle", &format!("{:.1}", shooting_angle));
        }
        if let Some(latitude) = request.latitude {
            req = req.query("latitude", &format!("{:.2}", latitude));
        }
        if let Some(longitude) = request.longitude {
            req = req.query("longitude", &format!("{:.2}", longitude));
        }
        if let Some(shot_direction) = request.shot_direction {
            req = req.query("shot_direction", &format!("{:.1}", shot_direction));
        }
        if let Some(twist_rate) = request.twist_rate {
            req = req.query("twist_rate", &format!("{:.1}", twist_rate));
        }
        if let Some(diameter) = request.bullet_diameter {
            let diameter_in = diameter / 0.0254; // meters to inches
            req = req.query("bullet_diameter", &format!("{:.3}", diameter_in));
        }
        if let Some(threshold) = request.ground_threshold {
            // Ground threshold is in meters - API expects meters (will be converted by parse_ballistic_inputs)
            // Use a very large negative value when ignoring ground impact
            req = req.query("ground_threshold", &format!("{:.1}", threshold));
        }
        if let Some(enable) = request.enable_weather_zones {
            req = req.query("enable_weather_zones", if enable { "true" } else { "false" });
        }
        if let Some(enable) = request.enable_3d_weather {
            req = req.query("enable_3d_weather", if enable { "true" } else { "false" });
        }
        if let Some(ref model) = request.wind_shear_model {
            req = req.query("wind_shear_model", model);
        }
        if let Some(ref method) = request.weather_zone_interpolation {
            req = req.query("weather_zone_interpolation", method);
        }
        if let Some(sample_interval) = request.sample_interval {
            // Convert from meters to yards for the Flask API (trajectory_step is in yards)
            let step_yards = sample_interval / 0.9144;
            req = req.query("trajectory_step", &format!("{:.4}", step_yards));
        }

        let response = req.call().map_err(|e| match e {
            ureq::Error::Status(code, response) => {
                let body = response.into_string().unwrap_or_default();
                ApiError::ServerError(code, body)
            }
            ureq::Error::Transport(transport) => {
                // Check for timeout by looking at the error message
                let msg = transport.to_string();
                if msg.contains("timed out") || msg.contains("timeout") {
                    ApiError::Timeout
                } else {
                    ApiError::NetworkError(msg)
                }
            }
        })?;

        let body = response
            .into_string()
            .map_err(|e| ApiError::InvalidResponse(e.to_string()))?;

        // Parse the Flask API response and convert to our format
        let api_response: serde_json::Value = serde_json::from_str(&body)
            .map_err(|e| ApiError::InvalidResponse(format!("JSON parse error: {}", e)))?;

        // Convert Flask API response to our TrajectoryResponse format
        self.convert_api_response(&api_response)
    }

    /// Helper to extract a value from a nested {value: x, unit: y} structure or plain number
    #[cfg(feature = "online")]
    fn extract_value(val: &serde_json::Value) -> Option<f64> {
        // Try nested {value: x} first, then plain number
        val.get("value")
            .and_then(|v| v.as_f64())
            .or_else(|| val.as_f64())
    }

    #[cfg(feature = "online")]
    fn convert_api_response(&self, api_response: &serde_json::Value) -> Result<TrajectoryResponse, ApiError> {
        // Get results object
        let results = api_response.get("results");

        // Extract trajectory points from Flask API response
        // The Flask API returns trajectory in "trajectory" array with nested value objects
        let trajectory_array = api_response.get("trajectory")
            .and_then(|t| t.as_array())
            .ok_or_else(|| ApiError::InvalidResponse("Missing trajectory array".to_string()))?;

        let trajectory: Vec<ApiTrajectoryPoint> = trajectory_array
            .iter()
            .filter_map(|point| {
                // Flask API returns nested {value: x, unit: y} objects in imperial units
                let range_yards = point.get("distance")
                    .and_then(Self::extract_value)?;
                let drop_inches = point.get("drop")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let drift_inches = point.get("wind_drift")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let velocity_fps = point.get("velocity")
                    .and_then(Self::extract_value)?;
                let energy_ftlbs = point.get("energy")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let time = point.get("time")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);

                Some(ApiTrajectoryPoint {
                    range: range_yards * 0.9144,        // yards to meters
                    drop: drop_inches * 0.0254,          // inches to meters
                    drift: drift_inches * 0.0254,        // inches to meters
                    velocity: velocity_fps * 0.3048,     // fps to m/s
                    energy: energy_ftlbs * 1.35582,      // ft-lbs to Joules
                    time,
                })
            })
            .collect();

        // Extract summary values from results object
        let zero_angle = results
            .and_then(|r| r.get("barrel_angle"))
            .and_then(Self::extract_value)
            .unwrap_or(0.0)
            .to_radians();

        let time_of_flight = results
            .and_then(|r| r.get("time_of_flight"))
            .and_then(Self::extract_value)
            .unwrap_or_else(|| trajectory.last().map(|p| p.time).unwrap_or(0.0));

        let bc_confidence = api_response.get("bc_confidence")
            .and_then(|v| v.as_f64());

        let ml_corrections = api_response.get("ml_corrections_applied")
            .or_else(|| api_response.get("corrections_applied"))
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_str().map(String::from))
                    .collect()
            });

        // max_height is directly a number in results
        let max_ordinate = results
            .and_then(|r| r.get("max_height"))
            .and_then(|v| v.as_f64())
            .map(|h| h * 0.0254); // inches to meters

        let impact_velocity = results
            .and_then(|r| r.get("final_velocity"))
            .and_then(Self::extract_value)
            .map(|v| v * 0.3048); // fps to m/s

        let impact_energy = results
            .and_then(|r| r.get("final_energy"))
            .and_then(Self::extract_value)
            .map(|e| e * 1.35582); // ft-lbs to Joules

        Ok(TrajectoryResponse {
            trajectory,
            zero_angle,
            time_of_flight,
            bc_confidence,
            ml_corrections_applied: ml_corrections,
            max_ordinate,
            impact_velocity,
            impact_energy,
        })
    }

    /// Check API health
    #[cfg(feature = "online")]
    pub fn health_check(&self) -> Result<bool, ApiError> {
        let url = format!("{}/health", self.base_url);

        let response = ureq::get(&url)
            .timeout(Duration::from_secs(5))
            .call()
            .map_err(|e| match e {
                ureq::Error::Status(code, response) => {
                    let body = response.into_string().unwrap_or_default();
                    ApiError::ServerError(code, body)
                }
                ureq::Error::Transport(transport) => {
                    let msg = transport.to_string();
                    if msg.contains("timed out") || msg.contains("timeout") {
                        ApiError::Timeout
                    } else {
                        ApiError::NetworkError(msg)
                    }
                }
            })?;

        Ok(response.status() == 200)
    }

    /// Calculate true/effective muzzle velocity via Flask API
    ///
    /// # Arguments
    /// * `request` - Velocity truing request parameters
    ///
    /// # Returns
    /// * `Ok(TrueVelocityResponse)` - Successful calculation with effective velocity
    /// * `Err(ApiError)` - Error during API communication
    #[cfg(feature = "online")]
    pub fn true_velocity(
        &self,
        request: &TrueVelocityRequest,
    ) -> Result<TrueVelocityResponse, ApiError> {
        let url = format!("{}/v1/true-velocity", self.base_url);

        let body = serde_json::to_string(request)
            .map_err(|e| ApiError::RequestError(format!("Failed to serialize request: {}", e)))?;

        let response = ureq::post(&url)
            .set("Content-Type", "application/json")
            .set("Accept", "application/json")
            .set("User-Agent", &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")))
            .timeout(self.timeout)
            .send_string(&body)
            .map_err(|e| match e {
                ureq::Error::Status(code, response) => {
                    let body = response.into_string().unwrap_or_default();
                    ApiError::ServerError(code, body)
                }
                ureq::Error::Transport(transport) => {
                    let msg = transport.to_string();
                    if msg.contains("timed out") || msg.contains("timeout") {
                        ApiError::Timeout
                    } else {
                        ApiError::NetworkError(msg)
                    }
                }
            })?;

        let response_body = response
            .into_string()
            .map_err(|e| ApiError::InvalidResponse(format!("Failed to read response: {}", e)))?;

        serde_json::from_str(&response_body)
            .map_err(|e| ApiError::InvalidResponse(format!("JSON parse error: {}", e)))
    }
}

/// Builder for TrajectoryRequest
#[derive(Default)]
pub struct TrajectoryRequestBuilder {
    bc_value: Option<f64>,
    bc_type: Option<String>,
    bullet_mass: Option<f64>,
    muzzle_velocity: Option<f64>,
    target_distance: Option<f64>,
    zero_range: Option<f64>,
    wind_speed: Option<f64>,
    wind_angle: Option<f64>,
    temperature: Option<f64>,
    pressure: Option<f64>,
    humidity: Option<f64>,
    altitude: Option<f64>,
    latitude: Option<f64>,
    longitude: Option<f64>,
    shot_direction: Option<f64>,
    shooting_angle: Option<f64>,
    twist_rate: Option<f64>,
    bullet_diameter: Option<f64>,
    bullet_length: Option<f64>,
    ground_threshold: Option<f64>,
    enable_weather_zones: Option<bool>,
    enable_3d_weather: Option<bool>,
    wind_shear_model: Option<String>,
    weather_zone_interpolation: Option<String>,
    sample_interval: Option<f64>,
}

impl TrajectoryRequestBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn bc_value(mut self, value: f64) -> Self {
        self.bc_value = Some(value);
        self
    }

    pub fn bc_type(mut self, value: &str) -> Self {
        self.bc_type = Some(value.to_string());
        self
    }

    pub fn bullet_mass(mut self, value: f64) -> Self {
        self.bullet_mass = Some(value);
        self
    }

    pub fn muzzle_velocity(mut self, value: f64) -> Self {
        self.muzzle_velocity = Some(value);
        self
    }

    pub fn target_distance(mut self, value: f64) -> Self {
        self.target_distance = Some(value);
        self
    }

    pub fn zero_range(mut self, value: f64) -> Self {
        self.zero_range = Some(value);
        self
    }

    pub fn wind_speed(mut self, value: f64) -> Self {
        self.wind_speed = Some(value);
        self
    }

    pub fn wind_angle(mut self, value: f64) -> Self {
        self.wind_angle = Some(value);
        self
    }

    pub fn temperature(mut self, value: f64) -> Self {
        self.temperature = Some(value);
        self
    }

    pub fn pressure(mut self, value: f64) -> Self {
        self.pressure = Some(value);
        self
    }

    pub fn humidity(mut self, value: f64) -> Self {
        self.humidity = Some(value);
        self
    }

    pub fn altitude(mut self, value: f64) -> Self {
        self.altitude = Some(value);
        self
    }

    pub fn latitude(mut self, value: f64) -> Self {
        self.latitude = Some(value);
        self
    }

    pub fn longitude(mut self, value: f64) -> Self {
        self.longitude = Some(value);
        self
    }

    pub fn shot_direction(mut self, value: f64) -> Self {
        self.shot_direction = Some(value);
        self
    }

    pub fn shooting_angle(mut self, value: f64) -> Self {
        self.shooting_angle = Some(value);
        self
    }

    pub fn twist_rate(mut self, value: f64) -> Self {
        self.twist_rate = Some(value);
        self
    }

    pub fn bullet_diameter(mut self, value: f64) -> Self {
        self.bullet_diameter = Some(value);
        self
    }

    pub fn bullet_length(mut self, value: f64) -> Self {
        self.bullet_length = Some(value);
        self
    }

    pub fn ground_threshold(mut self, value: f64) -> Self {
        self.ground_threshold = Some(value);
        self
    }

    pub fn enable_weather_zones(mut self, value: bool) -> Self {
        self.enable_weather_zones = Some(value);
        self
    }

    pub fn enable_3d_weather(mut self, value: bool) -> Self {
        self.enable_3d_weather = Some(value);
        self
    }

    pub fn wind_shear_model(mut self, value: &str) -> Self {
        self.wind_shear_model = Some(value.to_string());
        self
    }

    pub fn weather_zone_interpolation(mut self, value: &str) -> Self {
        self.weather_zone_interpolation = Some(value.to_string());
        self
    }

    pub fn sample_interval(mut self, value: f64) -> Self {
        self.sample_interval = Some(value);
        self
    }

    /// Build the TrajectoryRequest
    ///
    /// # Returns
    /// * `Ok(TrajectoryRequest)` - Valid request
    /// * `Err(String)` - Missing required fields
    pub fn build(self) -> Result<TrajectoryRequest, String> {
        let bc_value = self.bc_value.ok_or("bc_value is required")?;
        let bc_type = self.bc_type.ok_or("bc_type is required")?;
        let bullet_mass = self.bullet_mass.ok_or("bullet_mass is required")?;
        let muzzle_velocity = self.muzzle_velocity.ok_or("muzzle_velocity is required")?;
        let target_distance = self.target_distance.ok_or("target_distance is required")?;

        Ok(TrajectoryRequest {
            bc_value,
            bc_type,
            bullet_mass,
            muzzle_velocity,
            target_distance,
            zero_range: self.zero_range,
            wind_speed: self.wind_speed,
            wind_angle: self.wind_angle,
            temperature: self.temperature,
            pressure: self.pressure,
            humidity: self.humidity,
            altitude: self.altitude,
            latitude: self.latitude,
            longitude: self.longitude,
            shot_direction: self.shot_direction,
            shooting_angle: self.shooting_angle,
            twist_rate: self.twist_rate,
            bullet_diameter: self.bullet_diameter,
            bullet_length: self.bullet_length,
            ground_threshold: self.ground_threshold,
            enable_weather_zones: self.enable_weather_zones,
            enable_3d_weather: self.enable_3d_weather,
            wind_shear_model: self.wind_shear_model,
            weather_zone_interpolation: self.weather_zone_interpolation,
            sample_interval: self.sample_interval,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_request_builder_required_fields() {
        let result = TrajectoryRequestBuilder::new()
            .bc_value(0.238)
            .bc_type("G7")
            .bullet_mass(9.07) // 140gr in grams
            .muzzle_velocity(860.0)
            .target_distance(1000.0)
            .build();

        assert!(result.is_ok());
        let request = result.unwrap();
        assert_eq!(request.bc_value, 0.238);
        assert_eq!(request.bc_type, "G7");
    }

    #[test]
    fn test_request_builder_missing_fields() {
        let result = TrajectoryRequestBuilder::new()
            .bc_value(0.238)
            .build();

        assert!(result.is_err());
    }

    #[test]
    fn test_request_builder_all_optional_fields() {
        let result = TrajectoryRequestBuilder::new()
            .bc_value(0.238)
            .bc_type("G7")
            .bullet_mass(9.07)
            .muzzle_velocity(860.0)
            .target_distance(1000.0)
            .zero_range(100.0)
            .wind_speed(5.0)
            .wind_angle(90.0)
            .temperature(15.0)
            .pressure(1013.25)
            .humidity(50.0)
            .altitude(500.0)
            .latitude(45.0)
            .shooting_angle(0.0)
            .twist_rate(10.0)
            .bullet_diameter(0.00671)
            .bullet_length(0.035)
            .build();

        assert!(result.is_ok());
        let request = result.unwrap();
        assert_eq!(request.zero_range, Some(100.0));
        assert_eq!(request.wind_speed, Some(5.0));
        assert_eq!(request.latitude, Some(45.0));
    }

    #[test]
    fn test_api_client_url_normalization() {
        let client1 = ApiClient::new("https://api.example.com/", 10);
        assert_eq!(client1.base_url, "https://api.example.com");

        let client2 = ApiClient::new("https://api.example.com", 10);
        assert_eq!(client2.base_url, "https://api.example.com");
    }

    #[test]
    fn test_api_error_display() {
        assert_eq!(
            format!("{}", ApiError::NetworkError("connection refused".to_string())),
            "Network error: connection refused"
        );
        assert_eq!(format!("{}", ApiError::Timeout), "Request timed out");
        assert_eq!(
            format!("{}", ApiError::ServerError(500, "Internal error".to_string())),
            "Server error 500: Internal error"
        );
    }
}
