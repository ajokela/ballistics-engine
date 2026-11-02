//! HTTP client for Flask API communication
//!
//! This module provides an HTTP client for routing trajectory calculations
//! through the Flask API instead of local computation, giving CLI users
//! access to ML-enhanced predictions.

#[cfg(feature = "online")]
use crate::constants::GRAMS_PER_GRAIN;
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
    /// `Some(true)` when `temperature` is a reading the caller typed (0.49.0), so the service
    /// can fly it as typed even at 15 °C. `None` (never `Some(false)`) for an omitted value,
    /// whose request is then exactly what it was before 0.49.0. Sent as the
    /// `temperature_supplied=true` query parameter.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_supplied: Option<bool>,
    /// `Some(true)` when `pressure` is a station reading the caller typed (0.49.0); see
    /// [`Self::temperature_supplied`]. Sent as `pressure_supplied=true`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure_supplied: Option<bool>,
    /// Relative humidity, PERCENT 0-100 (optional). Sent at full precision; a sensor reading
    /// within 5 points of the range is sent as the end it reads as, and a value beyond that is
    /// refused before sending (`HUMIDITY_SENSOR_TOLERANCE_PERCENT`).
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

/// The humidity a request sends, checked by the engine's whole humidity rule before anything
/// goes over the network (R3, review UX-6): the sensor band
/// ([`crate::RelativeHumidity::try_percent`]) AND the impossible-air test every local solve
/// applies ([`crate::StationAir::new`], through
/// [`crate::AtmosphericConditions::air_refusal`]), on the station the request describes.
///
/// A sensor reading within `HUMIDITY_SENSOR_TOLERANCE_PERCENT` of 0-100 is sent as the end it
/// reads as (101.7 as 100), the value every local solve flies. A value beyond the band, or air
/// that cannot exist, is an [`ApiError::RequestError`] in the request's own words:
/// `humidity = 100 at 214 °F and 29.92 inHg: no such air exists: ...`. The station is resolved
/// as a local solve resolves it: an omitted temperature or pressure is the standard default,
/// which a solve flies as the ICAO value at the altitude (a typed one, `*_supplied`, flies as
/// typed), and an omitted humidity is the default 50 %.
///
/// `air` is the request's `(temperature °C, pressure hPa, altitude m)` converted from its own
/// units, `supplied` its temperature and pressure markers, and `spell` states a resolved
/// temperature (°C) and pressure (hPa) in the request's units.
/// How a request states a resolved reading in its own units: a temperature (°C) or a pressure
/// (hPa) as text with the unit named.
#[cfg(feature = "online")]
type ReadingText = fn(f64) -> String;

#[cfg(feature = "online")]
fn sendable_humidity(
    humidity: Option<f64>,
    air: (Option<f64>, Option<f64>, Option<f64>),
    supplied: (Option<bool>, Option<bool>),
    names: [&'static str; 4],
    spell: (ReadingText, ReadingText),
) -> Result<Option<f64>, ApiError> {
    let percent = humidity
        .map(|percent| {
            crate::RelativeHumidity::try_percent(percent)
                .map(crate::RelativeHumidity::percent)
                .map_err(|error| ApiError::RequestError(error.at_field(names[0]).to_string()))
        })
        .transpose()?;
    let (temperature_c, pressure_hpa, altitude_m) = air;
    let conditions = crate::AtmosphericConditions {
        temperature: temperature_c.unwrap_or(15.0),
        pressure: pressure_hpa.unwrap_or(1013.25),
        humidity: percent.unwrap_or(crate::DEFAULT_HUMIDITY_PERCENT),
        altitude: altitude_m.unwrap_or(0.0),
    };
    let presence = crate::StationAtmospherePresence::new(
        temperature_c.is_some() && supplied.0 == Some(true),
        pressure_hpa.is_some() && supplied.1 == Some(true),
    );
    let refusal = conditions
        .air_refusal(
            presence,
            names.map(crate::AirField::named),
            conditions.humidity,
            spell.0,
            spell.1,
        )
        .filter(crate::atmosphere::AirRefusal::is_humidity_rule);
    match refusal {
        Some(refusal) => Err(ApiError::RequestError(refusal.message().to_string())),
        None => Ok(percent),
    }
}

/// A temperature (°C) in a metric request's words.
#[cfg(feature = "online")]
fn celsius_text(celsius: f64) -> String {
    crate::atmosphere::reading_text(celsius, "°C")
}

/// A pressure (hPa) in a metric request's words.
#[cfg(feature = "online")]
fn hpa_text(hpa: f64) -> String {
    crate::atmosphere::reading_text(hpa, "hPa")
}

/// A temperature (°C) in an imperial request's words.
#[cfg(feature = "online")]
fn fahrenheit_text(celsius: f64) -> String {
    crate::atmosphere::reading_text(celsius * 9.0 / 5.0 + 32.0, "°F")
}

/// A pressure (hPa) in an imperial request's words.
#[cfg(feature = "online")]
fn inhg_text(hpa: f64) -> String {
    crate::atmosphere::reading_text(hpa / 33.8639, "inHg")
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
    /// `Some(true)` when `temperature_f` is a reading the caller typed (0.49.0); `None` for an
    /// omitted value, whose request is then exactly what it was before 0.49.0.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_supplied: Option<bool>,
    /// `Some(true)` when `pressure_inhg` is a station reading the caller typed (0.49.0); `None`
    /// for an omitted value.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pressure_supplied: Option<bool>,
    /// Relative humidity, PERCENT 0-100 (default: 50). A sensor reading within 5 points of the
    /// range is sent as the end it reads as, and a value beyond that is refused before sending.
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

impl ApiError {
    /// A CLI-friendly hint for auth/subscription failures, if applicable.
    pub fn cli_hint(&self) -> Option<&'static str> {
        match self {
            ApiError::ServerError(401, _) => Some(
                "Unauthorized - run `ballistics login` (create a token at https://ballisticsinsight.com/account).",
            ),
            ApiError::ServerError(402, _) => Some(
                "This endpoint needs an active Ballistics Insight subscription: https://ballisticsinsight.com/subscribe",
            ),
            _ => None,
        }
    }
}

impl std::error::Error for ApiError {}

/// HTTP client for Flask API communication
pub struct ApiClient {
    base_url: String,
    timeout: Duration,
    #[cfg_attr(not(feature = "online"), allow(dead_code))]
    token: Option<String>,
}

impl ApiClient {
    /// Create a new API client
    ///
    /// # Arguments
    /// * `base_url` - Base URL of the Flask API (e.g., <https://api.ballistics.7.62x51mm.sh>)
    /// * `timeout_secs` - Request timeout in seconds
    pub fn new(base_url: &str, timeout_secs: u64) -> Self {
        // Normalize URL by removing trailing slash
        let base_url = base_url.trim_end_matches('/').to_string();

        Self {
            base_url,
            timeout: Duration::from_secs(timeout_secs),
            token: None,
        }
    }

    /// Attach a Personal Access Token; sent as `Authorization: Bearer <token>` on requests.
    pub fn with_token(mut self, token: Option<String>) -> Self {
        self.token = token;
        self
    }

    /// The Authorization header pair to attach, if a token is configured.
    #[cfg(feature = "online")]
    fn auth_header(&self) -> Option<(&'static str, String)> {
        self.token
            .as_ref()
            .map(|t| ("Authorization", format!("Bearer {}", t)))
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
        // The engine's whole humidity rule (the band and air that cannot exist), before
        // anything is sent; the request is metric.
        let humidity = sendable_humidity(
            request.humidity,
            (request.temperature, request.pressure, request.altitude),
            (request.temperature_supplied, request.pressure_supplied),
            ["humidity", "temperature", "pressure", "altitude"],
            (celsius_text, hpa_text),
        )?;
        // Flask API uses GET /v1/calculate with query parameters (imperial units)
        let url = format!("{}/v1/calculate", self.base_url);

        // Convert metric values to imperial for API
        let velocity_fps = request.muzzle_velocity / 0.3048; // m/s to fps
        let mass_grains = request.bullet_mass / GRAMS_PER_GRAIN; // grams to grains
        let distance_yards = request.target_distance / 0.9144; // meters to yards

        let mut req = ureq::get(&url)
            .config()
            .timeout_global(Some(self.timeout))
            .build()
            .header("Accept", "application/json")
            .header(
                "User-Agent",
                &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")),
            )
            .query("bc_value", request.bc_value.to_string())
            .query("bc_type", &request.bc_type)
            .query("bullet_mass", format!("{:.1}", mass_grains))
            .query("muzzle_velocity", format!("{:.1}", velocity_fps))
            .query("target_distance", format!("{:.1}", distance_yards));

        if let Some((name, value)) = self.auth_header() {
            req = req.header(name, &value);
        }

        // Add optional parameters
        if let Some(zero_range) = request.zero_range {
            let zero_yards = zero_range / 0.9144;
            req = req.query("zero_distance", format!("{:.1}", zero_yards));
        }
        if let Some(wind_speed) = request.wind_speed {
            let wind_mph = wind_speed * 2.23694; // m/s to mph
            req = req.query("wind_speed", format!("{:.1}", wind_mph));
        }
        if let Some(wind_angle) = request.wind_angle {
            req = req.query("wind_angle", format!("{:.1}", wind_angle));
        }
        if let Some(temp) = request.temperature {
            let temp_f = temp * 9.0 / 5.0 + 32.0; // Celsius to Fahrenheit
            req = req.query("temperature", format!("{:.1}", temp_f));
        }
        if let Some(pressure) = request.pressure {
            let pressure_inhg = pressure / 33.8639; // hPa to inHg
            req = req.query("pressure", format!("{:.2}", pressure_inhg));
        }
        // 0.49.0: opt-in presence markers, only for a value the caller typed. An omitted value
        // sends no marker, so its request is exactly what it was before.
        if request.temperature_supplied == Some(true) {
            req = req.query("temperature_supplied", "true");
        }
        if request.pressure_supplied == Some(true) {
            req = req.query("pressure_supplied", "true");
        }
        // Full precision, the shortest decimal that round-trips (humidity train, plan C4): a
        // relative humidity is a percent on both sides, and `{:.1}` sent 37.375 as 37.4.
        if let Some(humidity) = humidity {
            req = req.query("humidity", humidity.to_string());
        }
        if let Some(altitude) = request.altitude {
            let altitude_ft = altitude / 0.3048; // meters to feet
            req = req.query("altitude", format!("{:.1}", altitude_ft));
        }
        if let Some(shooting_angle) = request.shooting_angle {
            req = req.query("shooting_angle", format!("{:.1}", shooting_angle));
        }
        if let Some(latitude) = request.latitude {
            req = req.query("latitude", format!("{:.2}", latitude));
        }
        if let Some(longitude) = request.longitude {
            req = req.query("longitude", format!("{:.2}", longitude));
        }
        if let Some(shot_direction) = request.shot_direction {
            req = req.query("shot_direction", format!("{:.1}", shot_direction));
        }
        if let Some(twist_rate) = request.twist_rate {
            req = req.query("twist_rate", format!("{:.1}", twist_rate));
        }
        if let Some(diameter) = request.bullet_diameter {
            let diameter_in = diameter / 0.0254; // meters to inches
            req = req.query("bullet_diameter", format!("{:.3}", diameter_in));
        }
        if let Some(threshold) = request.ground_threshold {
            // Ground threshold is in meters. --ignore-ground-impact sets f64::NEG_INFINITY, which
            // format!("{:.1}") would render as the literal "-inf" (rejected by most server-side
            // numeric validators, silently dropping the ignore-ground intent); send a large finite
            // negative sentinel for THAT case only. Finite values pass through. Any other
            // non-finite value (NaN or +Inf) is invalid input, not an ignore-ground request, so
            // omit the parameter rather than silently mapping it to the sentinel.
            if threshold.is_finite() {
                req = req.query("ground_threshold", format!("{:.1}", threshold));
            } else if threshold == f64::NEG_INFINITY {
                req = req.query("ground_threshold", "-1000000000.0");
            }
        }
        if let Some(enable) = request.enable_weather_zones {
            req = req.query(
                "enable_weather_zones",
                if enable { "true" } else { "false" },
            );
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
            req = req.query("trajectory_step", format!("{:.4}", step_yards));
        }

        let mut response = req.call().map_err(map_ureq_error)?;

        // In ureq 3.x, 4xx/5xx responses are no longer returned as Err by default;
        // .call() yields Ok(response) for any status. Inspect the status ourselves and
        // map non-2xx to ServerError(code, body) to preserve the previous behavior.
        let status = response.status();
        if !status.is_success() {
            let body = response.body_mut().read_to_string().unwrap_or_default();
            return Err(ApiError::ServerError(status.as_u16(), body));
        }

        let body = response
            .body_mut()
            .read_to_string()
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
    fn convert_api_response(
        &self,
        api_response: &serde_json::Value,
    ) -> Result<TrajectoryResponse, ApiError> {
        // Get results object
        let results = api_response.get("results");

        // Extract trajectory points from Flask API response
        // The Flask API returns trajectory in "trajectory" array with nested value objects
        let trajectory_array = api_response
            .get("trajectory")
            .and_then(|t| t.as_array())
            .ok_or_else(|| ApiError::InvalidResponse("Missing trajectory array".to_string()))?;

        let trajectory: Vec<ApiTrajectoryPoint> = trajectory_array
            .iter()
            .filter_map(|point| {
                // Flask API returns nested {value: x, unit: y} objects in imperial units
                let range_yards = point.get("distance").and_then(Self::extract_value)?;
                let drop_inches = point
                    .get("drop")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let drift_inches = point
                    .get("wind_drift")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let velocity_fps = point.get("velocity").and_then(Self::extract_value)?;
                let energy_ftlbs = point
                    .get("energy")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);
                let time = point
                    .get("time")
                    .and_then(Self::extract_value)
                    .unwrap_or(0.0);

                Some(ApiTrajectoryPoint {
                    range: range_yards * 0.9144,     // yards to meters
                    drop: drop_inches * 0.0254,      // inches to meters
                    drift: drift_inches * 0.0254,    // inches to meters
                    velocity: velocity_fps * 0.3048, // fps to m/s
                    energy: energy_ftlbs * 1.35582,  // ft-lbs to Joules
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

        let bc_confidence = api_response.get("bc_confidence").and_then(|v| v.as_f64());

        let ml_corrections = api_response
            .get("ml_corrections_applied")
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
            .config()
            .timeout_global(Some(Duration::from_secs(5)))
            .build()
            .call()
            .map_err(map_ureq_error)?;

        Ok(response.status().as_u16() == 200)
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

        // The engine's whole humidity rule (the band and air that cannot exist), before
        // anything is sent; the request is imperial, converted as the truing model converts it.
        let request = TrueVelocityRequest {
            humidity: sendable_humidity(
                request.humidity,
                (
                    request.temperature_f.map(|f| (f - 32.0) * 5.0 / 9.0),
                    request.pressure_inhg.map(|inhg| inhg * 33.8639),
                    request.altitude_ft.map(|ft| ft * 0.3048),
                ),
                (request.temperature_supplied, request.pressure_supplied),
                ["humidity", "temperature_f", "pressure_inhg", "altitude_ft"],
                (fahrenheit_text, inhg_text),
            )?,
            ..request.clone()
        };
        let body = serde_json::to_string(&request)
            .map_err(|e| ApiError::RequestError(format!("Failed to serialize request: {}", e)))?;

        let mut req = ureq::post(&url)
            .config()
            .timeout_global(Some(self.timeout))
            .build()
            .header("Content-Type", "application/json")
            .header("Accept", "application/json")
            .header(
                "User-Agent",
                &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")),
            );
        if let Some((name, value)) = self.auth_header() {
            req = req.header(name, &value);
        }
        let mut response = req.send(&body).map_err(map_ureq_error)?;

        // In ureq 3.x, 4xx/5xx responses are returned as Ok; check status explicitly.
        let status = response.status();
        if !status.is_success() {
            let body = response.body_mut().read_to_string().unwrap_or_default();
            return Err(ApiError::ServerError(status.as_u16(), body));
        }

        let response_body = response
            .body_mut()
            .read_to_string()
            .map_err(|e| ApiError::InvalidResponse(format!("Failed to read response: {}", e)))?;

        serde_json::from_str(&response_body)
            .map_err(|e| ApiError::InvalidResponse(format!("JSON parse error: {}", e)))
    }

    /// Generic authenticated POST of a JSON body to `<base_url><path>`, returning the parsed JSON
    /// response. Used by the online reverse-solver subcommands; sends the Bearer token.
    #[cfg(feature = "online")]
    pub fn post_json(
        &self,
        path: &str,
        body: &serde_json::Value,
    ) -> Result<serde_json::Value, ApiError> {
        let url = format!("{}{}", self.base_url, path);
        let body_str = serde_json::to_string(body)
            .map_err(|e| ApiError::RequestError(format!("Failed to serialize request: {}", e)))?;
        let mut req = ureq::post(&url)
            .config()
            .timeout_global(Some(self.timeout))
            .build()
            .header("Content-Type", "application/json")
            .header("Accept", "application/json")
            .header(
                "User-Agent",
                &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")),
            );
        if let Some((name, value)) = self.auth_header() {
            req = req.header(name, &value);
        }
        let mut response = req.send(&body_str).map_err(map_ureq_error)?;
        let status = response.status();
        let text = response
            .body_mut()
            .read_to_string()
            .map_err(|e| ApiError::InvalidResponse(format!("Failed to read response: {}", e)))?;
        if !status.is_success() {
            return Err(ApiError::ServerError(status.as_u16(), text));
        }
        serde_json::from_str(&text)
            .map_err(|e| ApiError::InvalidResponse(format!("JSON parse error: {}", e)))
    }

    /// Generic authenticated GET with query params to `<base_url><path>`, returning the parsed
    /// JSON response. Sibling of `post_json` for GET reverse-solver endpoints. Sends the Bearer.
    #[cfg(feature = "online")]
    pub fn get_json(
        &self,
        path: &str,
        params: &[(String, String)],
    ) -> Result<serde_json::Value, ApiError> {
        let url = format!("{}{}", self.base_url, path);
        let mut req = ureq::get(&url)
            .config()
            .timeout_global(Some(self.timeout))
            .build()
            .header("Accept", "application/json")
            .header(
                "User-Agent",
                &format!("ballistics-cli/{}", env!("CARGO_PKG_VERSION")),
            );
        for (k, v) in params {
            req = req.query(k.as_str(), v.as_str());
        }
        if let Some((name, value)) = self.auth_header() {
            req = req.header(name, &value);
        }
        let mut response = req.call().map_err(map_ureq_error)?;
        let status = response.status();
        let text = response
            .body_mut()
            .read_to_string()
            .map_err(|e| ApiError::InvalidResponse(format!("Failed to read response: {}", e)))?;
        if !status.is_success() {
            return Err(ApiError::ServerError(status.as_u16(), text));
        }
        serde_json::from_str(&text)
            .map_err(|e| ApiError::InvalidResponse(format!("JSON parse error: {}", e)))
    }
}

/// Map a ureq 3.x transport/protocol error to an `ApiError`.
///
/// Note: in ureq 3.x, HTTP 4xx/5xx responses are NOT returned as errors by default
/// (`.call()`/`.send()` yield `Ok(response)` for any status), so the status is checked
/// on the returned response rather than here. This helper only sees genuine
/// transport-level failures.
#[cfg(feature = "online")]
fn map_ureq_error(e: ureq::Error) -> ApiError {
    match e {
        // Any configured timeout (global/connect/etc.) maps to Timeout.
        ureq::Error::Timeout(_) => ApiError::Timeout,
        // Socket/IO errors: distinguish timeouts surfaced as IO errors, otherwise network.
        ureq::Error::Io(io_err) => {
            if io_err.kind() == std::io::ErrorKind::TimedOut {
                ApiError::Timeout
            } else {
                ApiError::NetworkError(io_err.to_string())
            }
        }
        // Other transport-level failures (DNS, connect, proxy, redirects, TLS, etc.).
        other => ApiError::NetworkError(other.to_string()),
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
    temperature_supplied: Option<bool>,
    pressure_supplied: Option<bool>,
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

    /// Mark the temperature as a typed reading (0.49.0). `false` leaves the marker off, which
    /// is what an omitted temperature sends.
    pub fn temperature_supplied(mut self, supplied: bool) -> Self {
        self.temperature_supplied = supplied.then_some(true);
        self
    }

    /// Mark the pressure as a typed station reading (0.49.0). `false` leaves the marker off.
    pub fn pressure_supplied(mut self, supplied: bool) -> Self {
        self.pressure_supplied = supplied.then_some(true);
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
            temperature_supplied: self.temperature_supplied,
            pressure_supplied: self.pressure_supplied,
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
        let result = TrajectoryRequestBuilder::new().bc_value(0.238).build();

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

    /// 0.49.0: the presence markers are opt-in. An omitted value's request serializes exactly
    /// as it did before 0.49.0 (no marker key at all); a typed one carries `true`.
    #[test]
    fn presence_markers_are_sent_only_for_typed_values() {
        let base = || {
            TrajectoryRequestBuilder::new()
                .bc_value(0.475)
                .bc_type("G1")
                .bullet_mass(10.9)
                .muzzle_velocity(823.0)
                .target_distance(600.0)
                .temperature(15.0)
                .pressure(1013.25)
                .altitude(150.0)
        };
        // Byte for byte what ab4eb0c's TrajectoryRequest serialized for these inputs.
        let omitted = serde_json::to_string(&base().build().unwrap()).unwrap();
        assert_eq!(
            omitted,
            r#"{"bc_value":0.475,"bc_type":"G1","bullet_mass":10.9,"muzzle_velocity":823.0,"target_distance":600.0,"temperature":15.0,"pressure":1013.25,"altitude":150.0}"#
        );
        let omitted_explicitly = base()
            .temperature_supplied(false)
            .pressure_supplied(false)
            .build()
            .unwrap();
        assert_eq!(serde_json::to_string(&omitted_explicitly).unwrap(), omitted);

        let typed = base()
            .temperature_supplied(false)
            .pressure_supplied(true)
            .build()
            .unwrap();
        let json = serde_json::to_string(&typed).unwrap();
        assert!(json.contains(r#""pressure_supplied":true"#), "{json}");
        assert!(!json.contains("temperature_supplied"), "{json}");
        let both = base()
            .temperature_supplied(true)
            .pressure_supplied(true)
            .build()
            .unwrap();
        assert_eq!(both.temperature_supplied, Some(true));
        assert_eq!(both.pressure_supplied, Some(true));

        let true_velocity = |temperature_supplied, pressure_supplied| TrueVelocityRequest {
            measured_drop_mil: 3.6,
            range_yd: 500.0,
            bc: 0.475,
            drag_model: "G1".to_string(),
            weight_gr: 168.0,
            caliber: 0.308,
            zero_range_yd: Some(100.0),
            chrono_velocity_fps: None,
            altitude_ft: Some(5000.0),
            temperature_f: Some(59.0),
            pressure_inhg: Some(29.92),
            temperature_supplied,
            pressure_supplied,
            humidity: Some(50.0),
            sight_height_in: Some(1.5),
            use_bc_enhancement: Some(true),
        };
        assert_eq!(
            serde_json::to_string(&true_velocity(None, None)).unwrap(),
            r#"{"measured_drop_mil":3.6,"range_yd":500.0,"bc":0.475,"drag_model":"G1","weight_gr":168.0,"caliber":0.308,"zero_range_yd":100.0,"altitude_ft":5000.0,"temperature_f":59.0,"pressure_inhg":29.92,"humidity":50.0,"sight_height_in":1.5,"use_bc_enhancement":true}"#
        );
        let json = serde_json::to_string(&true_velocity(Some(true), Some(true))).unwrap();
        assert!(json.contains(r#""temperature_supplied":true"#), "{json}");
        assert!(json.contains(r#""pressure_supplied":true"#), "{json}");
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
            format!(
                "{}",
                ApiError::NetworkError("connection refused".to_string())
            ),
            "Network error: connection refused"
        );
        assert_eq!(format!("{}", ApiError::Timeout), "Request timed out");
        assert_eq!(
            format!(
                "{}",
                ApiError::ServerError(500, "Internal error".to_string())
            ),
            "Server error 500: Internal error"
        );
    }

    #[cfg(feature = "online")]
    #[test]
    fn auth_header_present_when_token_set() {
        let c = ApiClient::new("https://x", 5).with_token(Some("bpat_abc".into()));
        assert_eq!(
            c.auth_header(),
            Some(("Authorization", "Bearer bpat_abc".to_string()))
        );
    }

    #[cfg(feature = "online")]
    #[test]
    fn auth_header_absent_without_token() {
        let c = ApiClient::new("https://x", 5);
        assert!(c.auth_header().is_none());
    }

    #[test]
    fn cli_hint_maps_auth_errors() {
        assert!(ApiError::ServerError(401, "x".into())
            .cli_hint()
            .unwrap()
            .contains("login"));
        assert!(ApiError::ServerError(402, "x".into())
            .cli_hint()
            .unwrap()
            .contains("subscription"));
        assert!(ApiError::ServerError(500, "x".into()).cli_hint().is_none());
    }

    /// R3 (whole-branch review UX-6): the client applies the engine's whole humidity rule, the
    /// sensor band AND the impossible-air test, before any byte goes over the network.
    /// 192.0.2.1 is TEST-NET-1 (RFC 5737) and never routed, so a request that tried to connect
    /// would come back as a timeout or a network error, never as a `RequestError`.
    #[cfg(feature = "online")]
    #[test]
    fn air_that_cannot_exist_is_refused_before_anything_is_sent() {
        let client = ApiClient::new("http://192.0.2.1:9", 2);
        let started = std::time::Instant::now();

        // 214 °F (a typo for 21.4) at 100 % and 29.92 inHg: saturated air at 101 °C needs more
        // vapour pressure than the air's total pressure.
        let true_velocity = TrueVelocityRequest {
            measured_drop_mil: 3.6,
            range_yd: 500.0,
            bc: 0.475,
            drag_model: "G1".to_string(),
            weight_gr: 168.0,
            caliber: 0.308,
            zero_range_yd: Some(100.0),
            chrono_velocity_fps: None,
            altitude_ft: Some(0.0),
            temperature_f: Some(214.0),
            pressure_inhg: Some(29.92),
            temperature_supplied: Some(true),
            pressure_supplied: Some(true),
            humidity: Some(100.0),
            sight_height_in: Some(1.5),
            use_bc_enhancement: Some(true),
        };
        match client.true_velocity(&true_velocity) {
            Err(ApiError::RequestError(message)) => assert!(
                message.starts_with("humidity = 100 at 214 °F and 29.92 inHg: no such air exists"),
                "{message}"
            ),
            other => panic!("true_velocity sent impossible air: {other:?}"),
        }

        let trajectory = TrajectoryRequestBuilder::new()
            .bc_value(0.475)
            .bc_type("G1")
            .bullet_mass(10.9)
            .muzzle_velocity(823.0)
            .target_distance(600.0)
            .temperature((214.0 - 32.0) * 5.0 / 9.0)
            .temperature_supplied(true)
            .pressure(1013.25)
            .pressure_supplied(true)
            .humidity(100.0)
            .build()
            .unwrap();
        match client.calculate_trajectory(&trajectory) {
            Err(ApiError::RequestError(message)) => assert!(
                message.starts_with("humidity = 100 at 101.11 °C and 1013.25 hPa: no such air"),
                "{message}"
            ),
            other => panic!("calculate_trajectory sent impossible air: {other:?}"),
        }
        assert!(
            started.elapsed() < std::time::Duration::from_secs(1),
            "both refusals came back without waiting on the network"
        );
    }

    /// The client check flies what a local solve flies: ordinary air passes, a sensor reading
    /// in the band is sent as its end, an omitted humidity stays omitted, and a humidity beyond
    /// the band is refused naming the field.
    #[cfg(feature = "online")]
    #[test]
    fn the_client_check_sends_what_a_local_solve_flies() {
        let names = ["humidity", "temperature", "pressure", "altitude"];
        let spell: (ReadingText, ReadingText) = (celsius_text, hpa_text);
        let typed = (Some(true), Some(true));
        let sea_level = (Some(15.0), Some(1013.25), Some(0.0));
        assert_eq!(
            sendable_humidity(Some(101.7), sea_level, typed, names, spell).unwrap(),
            Some(100.0)
        );
        assert_eq!(
            sendable_humidity(Some(37.375), sea_level, typed, names, spell).unwrap(),
            Some(37.375)
        );
        assert_eq!(
            sendable_humidity(None, (None, None, None), (None, None), names, spell).unwrap(),
            None
        );
        // Saturated air at 50 °C exists; at 101 °C it does not, typed or not.
        assert!(sendable_humidity(
            Some(100.0),
            (Some(50.0), Some(1013.25), Some(0.0)),
            typed,
            names,
            spell
        )
        .is_ok());
        assert!(matches!(
            sendable_humidity(
                Some(100.0),
                (Some(101.0), Some(1013.25), Some(0.0)),
                (None, None),
                names,
                spell
            ),
            Err(ApiError::RequestError(_))
        ));
        match sendable_humidity(Some(150.0), sea_level, typed, names, spell) {
            Err(ApiError::RequestError(message)) => assert!(
                message.starts_with("humidity = 150: relative humidity must be a percent"),
                "{message}"
            ),
            other => panic!("150 % must be refused: {other:?}"),
        }
    }
}
