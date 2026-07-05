// CLI API module - provides simplified interfaces for command-line tool
use crate::cluster_bc::ClusterBCDegradation;
use crate::pitch_damping::{calculate_pitch_damping_coefficient, PitchDampingCoefficients};
use crate::precession_nutation::{
    calculate_combined_angular_motion, AngularState, PrecessionNutationParams,
};
use crate::trajectory_sampling::{
    sample_trajectory, TrajectoryData, TrajectoryOutputs, TrajectorySample,
};
use crate::wind_shear::WindShearModel;
use crate::DragModel;
use nalgebra::Vector3;
use std::error::Error;
use std::fmt;

// Unit system for input/output
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitSystem {
    Imperial,
    Metric,
}

// Output format for results
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutputFormat {
    Table,
    Json,
    Csv,
}

// Error type for CLI operations
#[derive(Debug)]
pub struct BallisticsError {
    message: String,
}

impl fmt::Display for BallisticsError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for BallisticsError {}

impl From<String> for BallisticsError {
    fn from(msg: String) -> Self {
        BallisticsError { message: msg }
    }
}

impl From<&str> for BallisticsError {
    fn from(msg: &str) -> Self {
        BallisticsError {
            message: msg.to_string(),
        }
    }
}

// Ballistic input parameters - MBA-151 Reconciled Structure
// Unified structure used by both ballistics-engine and ballistics_rust
// Duplicates removed, all necessary fields included
#[derive(Debug, Clone)]
pub struct BallisticInputs {
    // Core ballistics parameters (using intuitive names)
    pub bc_value: f64,        // Ballistic coefficient (G1, G7, etc.)
    pub bc_type: DragModel,   // Drag model (G1, G7, G8, etc.)
    pub bullet_mass: f64,     // kg
    pub muzzle_velocity: f64, // m/s
    pub bullet_diameter: f64, // meters
    pub bullet_length: f64,   // meters

    // Targeting and positioning
    pub muzzle_angle: f64,    // radians (launch angle)
    pub target_distance: f64, // meters
    pub azimuth_angle: f64, // horizontal aiming angle in radians (small aim offset within the shot frame)
    /// Compass bearing the shot is fired ALONG, radians, 0 = North, π/2 = East.
    /// Used only by the Coriolis model (Earth-rotation depends on which way downrange
    /// points relative to true North). Distinct from `azimuth_angle`, which is the
    /// small horizontal *aiming* offset and rotates the launch velocity.
    pub shot_azimuth: f64,
    pub shooting_angle: f64,   // uphill/downhill angle in radians
    pub sight_height: f64,     // meters above bore
    pub muzzle_height: f64,    // meters above ground
    pub target_height: f64,    // meters above ground for zeroing
    pub ground_threshold: f64, // meters below which to stop

    // Environmental conditions
    pub altitude: f64,    // meters
    pub temperature: f64, // Celsius
    pub pressure: f64,    // millibars/hPa
    /// Relative humidity as a FRACTION in `[0, 1]` (e.g. 0.5 = 50%). NOTE the scale
    /// differs from [`AtmosphericConditions::humidity`], which is a PERCENT in `[0, 100]`.
    /// The atmosphere helpers (`calculate_air_density_*`) expect percent, so convert via
    /// [`BallisticInputs::humidity_percent`] before passing this value to them (MBA-722).
    pub humidity: f64,
    pub latitude: Option<f64>, // degrees

    // Wind conditions
    pub wind_speed: f64, // m/s
    pub wind_angle: f64, // radians (0=headwind, 90=from right)

    // Bullet characteristics
    pub twist_rate: f64,               // inches per turn
    pub is_twist_right: bool,          // right-hand twist
    pub caliber_inches: f64,           // diameter in inches
    pub weight_grains: f64,            // mass in grains
    pub manufacturer: Option<String>,  // Bullet manufacturer
    pub bullet_model: Option<String>,  // Bullet model name
    pub bullet_id: Option<String>,     // Unique bullet identifier
    pub bullet_cluster: Option<usize>, // BC cluster ID for cluster_bc module

    // Integration method selection
    pub use_rk4: bool,           // Use RK4 integration instead of Euler
    pub use_adaptive_rk45: bool, // Use RK45 adaptive step size integration

    // Advanced effects flags
    pub enable_advanced_effects: bool,
    pub enable_magnus: bool,   // Magnus side force (independent of Coriolis)
    pub enable_coriolis: bool, // Coriolis deflection (requires latitude)
    pub use_powder_sensitivity: bool,
    pub powder_temp_sensitivity: f64,
    pub powder_temp: f64,           // Celsius
    /// Optional measured powder-temperature -> muzzle-velocity curve, as
    /// (temperature_celsius, muzzle_velocity_m_s) points sorted ascending by
    /// temperature. When present it supersedes the linear `powder_temp_sensitivity`
    /// model: the muzzle velocity is interpolated from this table at the ambient
    /// `temperature` (clamped to the endpoints — no extrapolation beyond measured
    /// data). This is the data-driven, non-linear alternative to the constant slope.
    pub powder_temp_curve: Option<Vec<(f64, f64)>>,
    /// Temperature (Celsius) at which to interpolate `powder_temp_curve` — the POWDER
    /// temperature, which may differ from the ambient `temperature` (air). `None` uses
    /// `temperature`. Decouples the velocity lookup from the air-density temperature.
    pub powder_curve_temp_c: Option<f64>,
    pub tipoff_yaw: f64,            // radians
    pub tipoff_decay_distance: f64, // meters
    pub use_bc_segments: bool,
    pub bc_segments: Option<Vec<(f64, f64)>>, // Mach-BC pairs
    pub bc_segments_data: Option<Vec<crate::BCSegmentData>>, // Velocity-BC segments
    pub use_enhanced_spin_drift: bool,
    pub use_form_factor: bool,
    pub enable_wind_shear: bool,
    pub wind_shear_model: String,
    pub enable_trajectory_sampling: bool,
    pub sample_interval: f64, // meters
    pub enable_pitch_damping: bool,
    pub enable_precession_nutation: bool,
    // MBA-959: apply aerodynamic jump as a muzzle launch-angle perturbation.
    // EXPERIMENTAL — the underlying model is heuristic and not yet validated; default OFF.
    pub enable_aerodynamic_jump: bool,
    pub use_cluster_bc: bool, // Use cluster-based BC degradation

    // Custom drag model support
    pub custom_drag_table: Option<crate::drag::DragTable>,

    // Legacy field for compatibility
    pub bc_type_str: Option<String>,
}

impl BallisticInputs {
    /// `humidity` as a PERCENT in `[0, 100]`, clamped — the scale the atmosphere
    /// density helpers expect. Centralizes the 0–1 → 0–100 conversion so callers don't
    /// re-derive it (and can't accidentally feed the raw 0–1 fraction as a percentage).
    /// See the field doc on [`BallisticInputs::humidity`] (MBA-722).
    pub fn humidity_percent(&self) -> f64 {
        (self.humidity * 100.0).clamp(0.0, 100.0)
    }
}

impl Default for BallisticInputs {
    fn default() -> Self {
        let mass_kg = 0.01;
        let diameter_m = 0.00762;
        let bc = 0.5;
        let muzzle_angle_rad = 0.0;
        let bc_type = DragModel::G1;

        Self {
            // Core ballistics parameters
            bc_value: bc,
            bc_type,
            bullet_mass: mass_kg,
            muzzle_velocity: 800.0,
            bullet_diameter: diameter_m,
            bullet_length: diameter_m * 4.5, // Approximate (match the CLI's 4.5-caliber heuristic)

            // Targeting and positioning
            muzzle_angle: muzzle_angle_rad,
            target_distance: 100.0,
            azimuth_angle: 0.0,
            shot_azimuth: 0.0,
            shooting_angle: 0.0,
            sight_height: 0.05,
            muzzle_height: 0.0,       // Default 0 - height is in sight_height
            target_height: 0.0,       // Target at ground level by default
            ground_threshold: -100.0, // Effectively disable ground detection (allow bullet to drop 100m below start)

            // Environmental conditions
            altitude: 0.0,
            temperature: 15.0,
            pressure: 1013.25, // Standard sea level pressure (millibars)
            humidity: 0.5,     // 50% relative humidity
            latitude: None,

            // Wind conditions
            wind_speed: 0.0,
            wind_angle: 0.0,

            // Bullet characteristics
            twist_rate: 12.0, // 1:12" typical
            is_twist_right: true,
            caliber_inches: diameter_m / 0.0254, // Convert to inches
            weight_grains: mass_kg / 0.00006479891, // Convert to grains
            manufacturer: None,
            bullet_model: None,
            bullet_id: None,
            bullet_cluster: None,

            // Integration method selection
            use_rk4: true,           // Use Runge-Kutta methods by default
            use_adaptive_rk45: true, // Default to RK45 adaptive for best accuracy

            // Advanced effects (disabled by default)
            enable_advanced_effects: false,
            enable_magnus: false,
            enable_coriolis: false,
            use_powder_sensitivity: false,
            powder_temp_sensitivity: 0.0,
            powder_temp: 15.0,
            powder_temp_curve: None,
            powder_curve_temp_c: None,
            tipoff_yaw: 0.0,
            tipoff_decay_distance: 50.0,
            use_bc_segments: false,
            bc_segments: None,
            bc_segments_data: None,
            use_enhanced_spin_drift: false,
            use_form_factor: false,
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
            enable_trajectory_sampling: false,
            sample_interval: 10.0, // Default 10 meter intervals
            enable_pitch_damping: false,
            enable_precession_nutation: false,
            enable_aerodynamic_jump: false,
            use_cluster_bc: false, // Disabled by default for backward compatibility

            // Custom drag model support
            custom_drag_table: None,

            // Legacy field for compatibility
            bc_type_str: None,
        }
    }
}

/// Interpolate a muzzle velocity (m/s) from a measured powder-temperature curve at
/// `temp_c` (Celsius). `curve` is `(temperature_celsius, velocity_m_s)` points; it is
/// sorted ascending by temperature before use. Values below the first point or above
/// the last are CLAMPED to the endpoint velocity (no extrapolation beyond measured
/// data), and segments are linearly interpolated. A single point yields a constant.
pub fn interpolate_powder_temp_curve(curve: &[(f64, f64)], temp_c: f64) -> f64 {
    debug_assert!(!curve.is_empty());
    if curve.is_empty() {
        return 0.0;
    }
    // Defensive: accept unsorted input by sorting a local copy only when needed.
    // Callers (CLI/WASM parsers) already sort, so the common path is a no-op scan.
    let mut sorted;
    let pts: &[(f64, f64)] = if curve.windows(2).all(|w| w[0].0 <= w[1].0) {
        curve
    } else {
        sorted = curve.to_vec();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        &sorted
    };
    let n = pts.len();
    if temp_c <= pts[0].0 {
        return pts[0].1; // clamp below the coldest measured point
    }
    if temp_c >= pts[n - 1].0 {
        return pts[n - 1].1; // clamp above the hottest measured point
    }
    for i in 1..n {
        let (t0, v0) = pts[i - 1];
        let (t1, v1) = pts[i];
        if temp_c <= t1 {
            let span = t1 - t0;
            if span.abs() < f64::EPSILON {
                return v1; // coincident temps: avoid divide-by-zero, take the upper
            }
            let f = (temp_c - t0) / span;
            return v0 + f * (v1 - v0);
        }
    }
    pts[n - 1].1
}

// Wind conditions
#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64, // m/s
    // radians, wind-FROM convention: 0 = headwind, PI/2 = from the right,
    // PI = tailwind, 3*PI/2 = from the left (matches WindSock / the bindings).
    pub direction: f64,
}

impl Default for WindConditions {
    fn default() -> Self {
        Self {
            speed: 0.0,
            direction: 0.0,
        }
    }
}

// Atmospheric conditions
#[derive(Debug, Clone)]
pub struct AtmosphericConditions {
    pub temperature: f64, // Celsius
    pub pressure: f64,    // hPa
    /// Relative humidity as a PERCENT in `[0, 100]`. NOTE: [`BallisticInputs::humidity`]
    /// uses a 0–1 FRACTION instead — convert with `BallisticInputs::humidity_percent` when
    /// crossing between them (MBA-722).
    pub humidity: f64,
    pub altitude: f64, // meters
}

impl Default for AtmosphericConditions {
    fn default() -> Self {
        Self {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 0.0,
        }
    }
}

// Trajectory point data
#[derive(Debug, Clone)]
pub struct TrajectoryPoint {
    pub time: f64,
    pub position: Vector3<f64>,
    pub velocity_magnitude: f64,
    pub kinetic_energy: f64,
}

// Trajectory result
#[derive(Debug, Clone)]
pub struct TrajectoryResult {
    pub max_range: f64,
    pub max_height: f64,
    pub time_of_flight: f64,
    pub impact_velocity: f64,
    pub impact_energy: f64,
    pub points: Vec<TrajectoryPoint>,
    pub sampled_points: Option<Vec<TrajectorySample>>, // Trajectory samples at regular intervals
    pub min_pitch_damping: Option<f64>, // Minimum pitch damping coefficient (for stability warning)
    pub transonic_mach: Option<f64>,    // Mach number when entering transonic regime
    pub angular_state: Option<AngularState>, // Final angular state if precession/nutation enabled
    pub max_yaw_angle: Option<f64>,     // Maximum yaw angle during flight (radians)
    pub max_precession_angle: Option<f64>, // Maximum precession angle (radians)
    // MBA-959: aerodynamic-jump components applied at the muzzle (None unless
    // enable_aerodynamic_jump). EXPERIMENTAL.
    pub aerodynamic_jump: Option<crate::aerodynamic_jump::AerodynamicJumpComponents>,
}

impl TrajectoryResult {
    /// Interpolate position at a given downrange distance (X coordinate, McCoy).
    /// Returns the interpolated (x, y, z) position at that range.
    /// If the target range exceeds the trajectory, returns the last point.
    pub fn position_at_range(&self, target_range: f64) -> Option<Vector3<f64>> {
        if self.points.is_empty() {
            return None;
        }

        // Find the two points that bracket the target range
        for i in 0..self.points.len() - 1 {
            let p1 = &self.points[i];
            let p2 = &self.points[i + 1];

            // Check if target range is between these two points (X is downrange)
            if p1.position.x <= target_range && p2.position.x >= target_range {
                // Linear interpolation factor
                let dx = p2.position.x - p1.position.x;
                if dx.abs() < 1e-10 {
                    return Some(p1.position);
                }
                let t = (target_range - p1.position.x) / dx;

                // Interpolate Y and Z, use exact target_range for X
                return Some(Vector3::new(
                    target_range,
                    p1.position.y + t * (p2.position.y - p1.position.y),
                    p1.position.z + t * (p2.position.z - p1.position.z),
                ));
            }
        }

        // Target range is beyond trajectory - return last point
        self.points.last().map(|p| p.position)
    }
}

// Trajectory solver
pub struct TrajectorySolver {
    inputs: BallisticInputs,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
    max_range: f64,
    time_step: f64,
    cluster_bc: Option<ClusterBCDegradation>,
    /// Optional downrange-segmented wind. When `Some`, the per-step wind vector is
    /// looked up by downrange distance from this `WindSock` and the scalar `wind`
    /// field is ignored. When `None`, the constant `wind` vector is used (default),
    /// so a non-segmented solve is numerically identical to pre-feature behavior.
    wind_sock: Option<crate::wind::WindSock>,
}

impl TrajectorySolver {
    pub fn new(
        mut inputs: BallisticInputs,
        wind: WindConditions,
        atmosphere: AtmosphericConditions,
    ) -> Self {
        // Compute derived fields from base units
        inputs.caliber_inches = inputs.bullet_diameter / 0.0254;
        inputs.weight_grains = inputs.bullet_mass / 0.00006479891;

        // Resolve the muzzle velocity for the ambient temperature before integration.
        // A measured powder-temperature -> velocity curve (data-driven, non-linear)
        // takes precedence when supplied; otherwise fall back to the linear
        // powder-temperature-sensitivity model (MBA-963). Both operate in canonical
        // SI (Celsius, m/s) and are applied here so every solver built from these
        // inputs — the main trajectory AND the zero-angle search — sees the same
        // temperature-resolved velocity. In particular, when a zero solve passes the
        // zero-day temperature, the curve automatically yields the zero-day velocity.
        if let Some(curve) = inputs.powder_temp_curve.as_ref() {
            if !curve.is_empty() {
                // Interpolate at the POWDER temperature, which defaults to the ambient
                // air temperature but can be decoupled (powder warmed/cooled relative to
                // the air) via powder_curve_temp_c. Air temperature still drives density
                // separately; this only sets the velocity. Absolute override (idempotent).
                let lookup_c = inputs.powder_curve_temp_c.unwrap_or(inputs.temperature);
                inputs.muzzle_velocity = interpolate_powder_temp_curve(curve, lookup_c);
            }
        } else if inputs.use_powder_sensitivity {
            let temp_delta_c = inputs.temperature - inputs.powder_temp;
            inputs.muzzle_velocity += inputs.powder_temp_sensitivity * temp_delta_c;
        }

        // Initialize cluster BC if enabled
        let cluster_bc = if inputs.use_cluster_bc {
            Some(ClusterBCDegradation::new())
        } else {
            None
        };

        Self {
            inputs,
            wind,
            atmosphere,
            max_range: 1000.0,
            time_step: 0.001,
            cluster_bc,
            wind_sock: None,
        }
    }

    pub fn set_max_range(&mut self, range: f64) {
        self.max_range = range;
    }

    pub fn set_time_step(&mut self, step: f64) {
        self.time_step = step;
    }

    /// Supply downrange-segmented wind. Each segment is `(speed_kmh, angle_deg,
    /// until_distance_m)`; the wind for a given downrange distance is the first
    /// segment whose `until_distance_m` exceeds it (a step function), and wind is
    /// zero beyond the last segment. An empty list clears segmented wind (reverts
    /// to the scalar `wind`). The angle convention matches `WindConditions`
    /// (0 = headwind, 90 = from the right).
    pub fn set_wind_segments(&mut self, segments: Vec<crate::wind::WindSegment>) {
        self.wind_sock = if segments.is_empty() {
            None
        } else {
            Some(crate::wind::WindSock::new(segments))
        };
    }

    /// Effective initial launch direction `(elevation, azimuth)` in radians, including
    /// the aerodynamic-jump muzzle perturbation when `enable_aerodynamic_jump` is set.
    ///
    /// Aerodynamic jump is the fixed angular departure imparted as the projectile
    /// transitions from the constrained bore to free flight; applying it as an initial
    /// launch-angle offset is the physically correct integration point. Returns the bare
    /// `(muzzle_angle, azimuth_angle)` when the flag is off, so a default solve is
    /// numerically identical to pre-feature behavior. (MBA-959)
    fn launch_angles_from(
        &self,
        aj: Option<&crate::aerodynamic_jump::AerodynamicJumpComponents>,
    ) -> (f64, f64) {
        let elev = self.inputs.muzzle_angle;
        let azim = self.inputs.azimuth_angle;
        match aj {
            Some(c) => {
                // vertical_/horizontal_jump_moa ARE the jump angles expressed in MOA.
                const MOA_PER_RAD: f64 = 3437.7467707849;
                (
                    elev + c.vertical_jump_moa / MOA_PER_RAD,
                    azim + c.horizontal_jump_moa / MOA_PER_RAD,
                )
            }
            None => (elev, azim),
        }
    }

    /// Compute the aerodynamic-jump components for the current inputs, or `None` when the
    /// feature is disabled / inputs are degenerate.
    ///
    /// Uses Bryan Litz's crosswind aerodynamic-jump estimator
    /// (`Y = 0.01*Sg - 0.0024*L + 0.032` MOA/mph) fed by the engine's own Miller Sg.
    /// Aerodynamic jump is a vertical effect, so only the elevation is perturbed.
    /// The estimator is a regression best near Sg ~ 1.75 — see MBA-959.
    fn aerodynamic_jump_components(
        &self,
    ) -> Option<crate::aerodynamic_jump::AerodynamicJumpComponents> {
        if !self.inputs.enable_aerodynamic_jump {
            return None;
        }
        // Reject degenerate/non-finite inputs before they can reach the launch angle.
        // A bare `<= 0.0` test lets NaN through (NaN comparisons are always false), and a
        // NaN/Inf here would poison the muzzle angle and collapse the whole trajectory.
        let diameter_m = self.inputs.bullet_diameter;
        if !(self.inputs.twist_rate.is_finite() && self.inputs.twist_rate != 0.0)
            || !(diameter_m.is_finite() && diameter_m > 0.0)
            || !(self.inputs.bullet_length.is_finite() && self.inputs.bullet_length > 0.0)
            || !self.inputs.muzzle_velocity.is_finite()
        {
            return None;
        }

        // Engine's own gyroscopic (Miller) stability factor — same Sg shown elsewhere.
        let (_, _, temp_c, pressure_hpa) = self.resolved_atmosphere();
        let sg = crate::stability::compute_stability_coefficient(
            &self.inputs,
            (self.atmosphere.altitude, temp_c, pressure_hpa, 0.0),
        );
        if !(sg.is_finite() && sg > 0.0) {
            return None;
        }
        let length_calibers = self.inputs.bullet_length / diameter_m;

        // Crosswind-from-the-right (mph) for Litz's estimator. Wind direction uses the
        // wind-FROM convention (0 = headwind, +90deg = from the right), matching the
        // fast-integrate path (fast_trajectory::aerodynamic_jump_launch_offset_rad) and
        // the lateral windage sign, so a from-the-right wind on a right-twist barrel
        // jumps the impact UP and drifts it left.
        const MS_TO_MPH: f64 = 2.236_936_292_054_4;
        let crosswind_from_right_mph = self.wind.speed * self.wind.direction.sin() * MS_TO_MPH;

        let vertical_jump_moa = crate::aerodynamic_jump::litz_crosswind_jump_moa(
            sg,
            length_calibers,
            crosswind_from_right_mph,
            self.inputs.is_twist_right,
        );
        if !vertical_jump_moa.is_finite() {
            return None;
        }

        const MOA_PER_RAD: f64 = 3437.7467707849;
        Some(crate::aerodynamic_jump::AerodynamicJumpComponents {
            vertical_jump_moa,
            // Aerodynamic jump is a vertical effect; the Litz estimator has no horizontal term.
            horizontal_jump_moa: 0.0,
            jump_angle_rad: vertical_jump_moa.abs() / MOA_PER_RAD,
            magnus_component_moa: 0.0,
            yaw_component_moa: 0.0,
            stabilization_factor: (sg / 1.5).clamp(0.0, 1.0),
        })
    }

    fn resolved_atmosphere(&self) -> (f64, f64, f64, f64) {
        let (temp_c, pressure_hpa) = crate::atmosphere::resolve_station_conditions(
            self.atmosphere.temperature,
            self.atmosphere.pressure,
            self.atmosphere.altitude,
        );
        let (density, speed_of_sound) = crate::atmosphere::calculate_atmosphere(
            self.atmosphere.altitude,
            Some(temp_c),
            Some(pressure_hpa),
            self.atmosphere.humidity,
        );
        (density, speed_of_sound, temp_c, pressure_hpa)
    }

    fn gravity_acceleration(&self) -> Vector3<f64> {
        let theta = self.inputs.shooting_angle;
        Vector3::new(
            -crate::constants::G_ACCEL_MPS2 * theta.sin(),
            -crate::constants::G_ACCEL_MPS2 * theta.cos(),
            0.0,
        )
    }

    fn get_wind_at_altitude(&self, altitude_m: f64) -> Vector3<f64> {
        // Scale the operative surface wind by the boundary-layer multiplier. `altitude_m` is the
        // bullet's height relative to the muzzle (McCoy Y). The multiplier is floored at 1.0, so
        // flat-fire trajectories keep ~full wind and only high-arcing shots see increased wind.
        //
        // We build the vector with THIS solver's non-shear sign convention (X=-cos, Z=-sin; see
        // the `wind_vector` used in solve_rk4/solve_euler, matching WindSock) and scale it, so that
        // "shear on" equals "shear off" * ratio (ratio == 1.0 for flat fire). An earlier revision
        // attenuated the wind near the line of sight and flipped its sign relative to the non-shear
        // path; this keeps them sign-consistent.
        // Map the requested model name to the boundary-layer model (MBA-965).
        // Names match wind_shear::get_wind_at_position. Unknown strings should
        // never reach here (the CLI parses an enum), but default to PowerLaw to
        // preserve the historical "exponential" behaviour for any caller that
        // forwards an unexpected value.
        let model = match self.inputs.wind_shear_model.as_str() {
            "logarithmic" => WindShearModel::Logarithmic,
            "power_law" | "powerlaw" | "exponential" => WindShearModel::PowerLaw,
            "ekman_spiral" | "ekman" => WindShearModel::EkmanSpiral,
            "custom_layers" | "custom" => WindShearModel::CustomLayers,
            _ => WindShearModel::PowerLaw,
        };
        let speed_ratio = crate::wind_shear::boundary_layer_speed_ratio(altitude_m, model);

        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindConditions / WindSock); wind enters drag via velocity - wind.
        Vector3::new(
            -self.wind.speed * self.wind.direction.cos() * speed_ratio, // X: downrange head/tail
            0.0,
            -self.wind.speed * self.wind.direction.sin() * speed_ratio, // Z: lateral crosswind
        )
    }

    pub fn solve(&self) -> Result<TrajectoryResult, BallisticsError> {
        let mut result = if self.inputs.use_rk4 {
            if self.inputs.use_adaptive_rk45 {
                self.solve_rk45()?
            } else {
                self.solve_rk4()?
            }
        } else {
            self.solve_euler()?
        };
        self.apply_spin_drift(&mut result);
        Ok(result)
    }

    /// Gyroscopic spin drift via the empirical Litz model, applied in the engine
    /// (not the WASM formatter) so it covers Euler/RK4/RK45 and all consumers.
    /// Uses the canonical SI fields and converts to grains/inches correctly,
    /// avoiding the kg/m-vs-grains/in unit bug in `calculate_enhanced_spin_drift`.
    /// Frame (McCoy): Z = lateral (windage), so drift adds to `position.z`.
    fn apply_spin_drift(&self, result: &mut TrajectoryResult) {
        if !self.inputs.use_enhanced_spin_drift {
            return;
        }
        let d_in = self.inputs.bullet_diameter / 0.0254; // m -> in
        let m_gr = self.inputs.bullet_mass / 0.00006479891; // kg -> grains
        let twist_in = self.inputs.twist_rate; // inches/turn
        if d_in <= 0.0 || m_gr <= 0.0 || twist_in <= 0.0 {
            return;
        }

        // Real length when available, else 4.5 cal (typical match bullet).
        let length_in = if self.inputs.bullet_length > 0.0 {
            self.inputs.bullet_length / 0.0254
        } else {
            4.5 * d_in
        };
        // MBA-942: apply the canonical Miller atmospheric correction (LINEAR in density ratio,
        // = rho0/rho via ideal gas: (T/T0)*(P0/P)), matching stability.rs and py_ballisticcalc.
        // miller_stability returns the bare geometric Sg with no density dependence, so without
        // this the spin drift under-predicts at altitude (Sg should rise as the air thins). At
        // standard sea level (15 C, 1013.25 hPa) the factor is exactly 1.0 — a no-op there.
        let (_, _, temp_c, press_hpa) = self.resolved_atmosphere();
        let temp_k = temp_c + 273.15; // Celsius -> Kelvin
        let density_correction = if press_hpa > 0.0 && temp_k > 0.0 {
            (temp_k / 288.15) * (1013.25 / press_hpa)
        } else {
            1.0
        };
        let sg = crate::spin_drift::miller_stability(d_in, m_gr, twist_in, length_in)
            * density_correction;
        let sign = if self.inputs.is_twist_right {
            1.0
        } else {
            -1.0
        };

        for p in result.points.iter_mut() {
            if p.time <= 0.0 {
                continue;
            }
            let sd_in = 1.25 * (sg + 1.2) * p.time.powf(1.83); // Litz drift, inches
            p.position.z += sign * sd_in * 0.0254; // in -> m, Z = lateral
        }

        // sampled_points are snapshotted from the PRE-drift trajectory inside each solver, so the
        // sampled wind_drift_m column would omit the spin drift that result.points carry. Apply
        // the same Litz drift to keep the two user-facing outputs consistent.
        if let Some(samples) = result.sampled_points.as_mut() {
            for s in samples.iter_mut() {
                if s.time_s <= 0.0 {
                    continue;
                }
                let sd_in = 1.25 * (sg + 1.2) * s.time_s.powf(1.83);
                s.wind_drift_m += sign * sd_in * 0.0254;
            }
        }
    }

    fn solve_euler(&self) -> Result<TrajectoryResult, BallisticsError> {
        // Simple trajectory integration using Euler method
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        let mut position = Vector3::new(
            0.0,
            self.inputs.muzzle_height, // Bore position above ground (NOT + sight_height)
            0.0,
        );
        // Calculate initial velocity components with both elevation and azimuth
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic
                                       // Downrange distances where the projectile crosses Mach 1.2 (transonic) then Mach 1.0
                                       // (subsonic), so the sampled trajectory output can flag those transitions
                                       // (trajectory_sampling::add_trajectory_flags consumes this).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut crossed_transonic = false;
        let mut crossed_subsonic = false;

        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001, // Small initial disturbance
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        // Calculate air density
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) = self.resolved_atmosphere();

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        let wind_vector = Vector3::new(
            -self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
            0.0,
            -self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
        );

        // Pitch-damping coefficients depend only on the (constant) bullet_model; compute once
        // instead of re-deriving them (with a to_lowercase alloc) every integration step.
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );

        // Main integration loop (X is downrange)
        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < 100.0
        {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy =
                0.5 * self.inputs.bullet_mass * velocity_magnitude * velocity_magnitude;

            points.push(TrajectoryPoint {
                time,
                position: position,
                velocity_magnitude,
                kinetic_energy,
            });

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                if !crossed_transonic && mach_here < 1.2 {
                    crossed_transonic = true;
                    transonic_distances.push(position.x);
                }
                if !crossed_subsonic && mach_here < 1.0 {
                    crossed_subsonic = true;
                    transonic_distances.push(position.x);
                }
            }

            // Debug: log first and every 100th point. Debug builds only — this was ungated and
            // polluted release/WASM stderr on the --use-euler path (the other solvers have none).
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
            #[cfg(debug_assertions)]
            if points.len() == 1 || points.len() % 100 == 0 {
                eprintln!("Trajectory point {}: time={:.3}s, downrange={:.2}m, vertical={:.2}m, lateral={:.2}m, vel={:.1}m/s",
                    points.len(), time, position.x, position.y, position.z, velocity_magnitude);
            }

            // Track max height
            if position.y > max_height {
                max_height = position.y;
            }

            // Calculate pitch damping if enabled
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let mach = velocity_magnitude / speed_of_sound;

                    // Calculate spin rate from twist rate and velocity
                    let spin_rate_rad_s = if self.inputs.twist_rate > 0.0 {
                        let velocity_fps = velocity_magnitude * 3.28084;
                        let twist_rate_ft = self.inputs.twist_rate / 12.0;
                        (velocity_fps / twist_rate_ft) * 2.0 * std::f64::consts::PI
                    } else {
                        0.0
                    };

                    // Create precession/nutation parameters
                    let params = PrecessionNutationParams {
                        mass_kg: self.inputs.bullet_mass,
                        caliber_m: self.inputs.bullet_diameter,
                        length_m: self.inputs.bullet_length,
                        spin_rate_rad_s,
                        spin_inertia: 6.94e-8,       // Typical value
                        transverse_inertia: 9.13e-7, // Typical value
                        velocity_mps: velocity_magnitude,
                        air_density_kg_m3: air_density,
                        mach,
                        pitch_damping_coeff: -0.8,
                        nutation_damping_factor: 0.05,
                    };

                    // Update angular state
                    *state = calculate_combined_angular_motion(
                        &params,
                        state,
                        time,
                        self.time_step,
                        0.001, // Initial disturbance
                    );

                    // Track maximums
                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            // Use the same acceleration kernel as RK4/RK45 so all three solvers share ONE drag
            // model. solve_euler previously used a bespoke frontal-area drag (0.5*rho*Cd*A*v^2/m)
            // that IGNORED the ballistic coefficient entirely (diverging up to ~2.3x from the
            // BC-retardation RK4/RK45 path), and also omitted the Magnus/Coriolis terms.
            // calculate_acceleration applies BC-retardation drag, gravity, Coriolis, Magnus, wind
            // shear, and the zero-relative-velocity gravity-only guard.
            let acceleration =
                self.calculate_acceleration(&position, &velocity, air_density,
                    &wind_vector,
                    (speed_of_sound, resolved_temp_c, resolved_press_hpa),
                );

            // Update state
            velocity += acceleration * self.time_step;
            position += velocity * self.time_step;
            time += self.time_step;
        }

        // Get final values
        let last_point = points.last().ok_or("No trajectory points generated")?;

        // Create trajectory sampling data if enabled
        let sampled_points = if self.inputs.enable_trajectory_sampling {
            let trajectory_data = TrajectoryData {
                times: points.iter().map(|p| p.time).collect(),
                positions: points.iter().map(|p| p.position).collect(),
                velocities: points
                    .iter()
                    .map(|p| {
                        // Reconstruct velocity vectors from magnitude (approximate)
                        Vector3::new(0.0, 0.0, p.velocity_magnitude)
                    })
                    .collect(),
                transonic_distances, // populated above at each Mach-threshold crossing
            };

            // For LOS calculation in ground-referenced coordinates:
            // sight_position_m is the sight's actual y-position above ground
            // (muzzle_height + sight_height, not just sight_height)
            // For flat shots, target is at same height as the sight (horizontal LOS)
            let sight_position_m = self.inputs.muzzle_height + self.inputs.sight_height;
            let outputs = TrajectoryOutputs {
                target_distance_horiz_m: last_point.position.x, // X is downrange
                target_vertical_height_m: sight_position_m,
                time_of_flight_s: last_point.time,
                max_ord_dist_horiz_m: max_height,
                sight_height_m: sight_position_m,
            };

            // Sample at specified intervals
            let samples = sample_trajectory(
                &trajectory_data,
                &outputs,
                self.inputs.sample_interval,
                self.inputs.bullet_mass,
            );
            Some(samples)
        } else {
            None
        };

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
        })
    }

    fn solve_rk4(&self) -> Result<TrajectoryResult, BallisticsError> {
        // RK4 trajectory integration for better accuracy
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        // The sight_height affects the LOS calculation and zero angle, not the starting position
        let mut position = Vector3::new(
            0.0,
            self.inputs.muzzle_height, // Bore position above ground (NOT + sight_height)
            0.0,
        );

        // Calculate initial velocity components with both elevation and azimuth
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic
                                       // Downrange distances where the projectile crosses Mach 1.2 (transonic) then Mach 1.0
                                       // (subsonic), so the sampled trajectory output can flag those transitions
                                       // (trajectory_sampling::add_trajectory_flags consumes this).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut crossed_transonic = false;
        let mut crossed_subsonic = false;

        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001, // Small initial disturbance
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        // Calculate air density
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) = self.resolved_atmosphere();

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        let wind_vector = Vector3::new(
            -self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
            0.0,
            -self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
        );

        // Pitch-damping coefficients depend only on the (constant) bullet_model; compute once
        // instead of re-deriving them (with a to_lowercase alloc) every integration step.
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );

        // Main RK4 integration loop (X is downrange)
        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < 100.0
        {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy =
                0.5 * self.inputs.bullet_mass * velocity_magnitude * velocity_magnitude;

            points.push(TrajectoryPoint {
                time,
                position: position,
                velocity_magnitude,
                kinetic_energy,
            });

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                if !crossed_transonic && mach_here < 1.2 {
                    crossed_transonic = true;
                    transonic_distances.push(position.x);
                }
                if !crossed_subsonic && mach_here < 1.0 {
                    crossed_subsonic = true;
                    transonic_distances.push(position.x);
                }
            }

            if position.y > max_height {
                max_height = position.y;
            }

            // Calculate pitch damping if enabled (RK4 solver)
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled (RK4 solver)
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let mach = velocity_magnitude / speed_of_sound;

                    // Calculate spin rate from twist rate and velocity
                    let spin_rate_rad_s = if self.inputs.twist_rate > 0.0 {
                        let velocity_fps = velocity_magnitude * 3.28084;
                        let twist_rate_ft = self.inputs.twist_rate / 12.0;
                        (velocity_fps / twist_rate_ft) * 2.0 * std::f64::consts::PI
                    } else {
                        0.0
                    };

                    // Create precession/nutation parameters
                    let params = PrecessionNutationParams {
                        mass_kg: self.inputs.bullet_mass,
                        caliber_m: self.inputs.bullet_diameter,
                        length_m: self.inputs.bullet_length,
                        spin_rate_rad_s,
                        spin_inertia: 6.94e-8,       // Typical value
                        transverse_inertia: 9.13e-7, // Typical value
                        velocity_mps: velocity_magnitude,
                        air_density_kg_m3: air_density,
                        mach,
                        pitch_damping_coeff: -0.8,
                        nutation_damping_factor: 0.05,
                    };

                    // Update angular state
                    *state = calculate_combined_angular_motion(
                        &params,
                        state,
                        time,
                        self.time_step,
                        0.001, // Initial disturbance
                    );

                    // Track maximums
                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            // RK4 method
            let dt = self.time_step;

            // k1
            let acc1 = self.calculate_acceleration(&position, &velocity, air_density, &wind_vector, (speed_of_sound, resolved_temp_c, resolved_press_hpa));

            // k2
            let pos2 = position + velocity * (dt * 0.5);
            let vel2 = velocity + acc1 * (dt * 0.5);
            let acc2 = self.calculate_acceleration(&pos2, &vel2, air_density, &wind_vector, (speed_of_sound, resolved_temp_c, resolved_press_hpa));

            // k3
            let pos3 = position + vel2 * (dt * 0.5);
            let vel3 = velocity + acc2 * (dt * 0.5);
            let acc3 = self.calculate_acceleration(&pos3, &vel3, air_density, &wind_vector, (speed_of_sound, resolved_temp_c, resolved_press_hpa));

            // k4
            let pos4 = position + vel3 * dt;
            let vel4 = velocity + acc3 * dt;
            let acc4 = self.calculate_acceleration(&pos4, &vel4, air_density, &wind_vector, (speed_of_sound, resolved_temp_c, resolved_press_hpa));

            // Update position and velocity
            position += (velocity + vel2 * 2.0 + vel3 * 2.0 + vel4) * (dt / 6.0);
            velocity += (acc1 + acc2 * 2.0 + acc3 * 2.0 + acc4) * (dt / 6.0);
            time += dt;
        }

        // Get final values
        let last_point = points.last().ok_or("No trajectory points generated")?;

        // Create trajectory sampling data if enabled
        let sampled_points = if self.inputs.enable_trajectory_sampling {
            let trajectory_data = TrajectoryData {
                times: points.iter().map(|p| p.time).collect(),
                positions: points.iter().map(|p| p.position).collect(),
                velocities: points
                    .iter()
                    .map(|p| {
                        // Reconstruct velocity vectors from magnitude (approximate)
                        Vector3::new(0.0, 0.0, p.velocity_magnitude)
                    })
                    .collect(),
                transonic_distances, // populated above at each Mach-threshold crossing
            };

            // For LOS calculation in ground-referenced coordinates:
            // sight_position_m is the sight's actual y-position above ground
            // (muzzle_height + sight_height, not just sight_height)
            // For flat shots, target is at same height as the sight (horizontal LOS)
            let sight_position_m = self.inputs.muzzle_height + self.inputs.sight_height;
            let outputs = TrajectoryOutputs {
                target_distance_horiz_m: last_point.position.x, // X is downrange
                target_vertical_height_m: sight_position_m,
                time_of_flight_s: last_point.time,
                max_ord_dist_horiz_m: max_height,
                sight_height_m: sight_position_m,
            };

            // Sample at specified intervals
            let samples = sample_trajectory(
                &trajectory_data,
                &outputs,
                self.inputs.sample_interval,
                self.inputs.bullet_mass,
            );
            Some(samples)
        } else {
            None
        };

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
        })
    }

    fn solve_rk45(&self) -> Result<TrajectoryResult, BallisticsError> {
        // RK45 adaptive step size integration (Dormand-Prince method)
        let mut time = 0.0;
        // Bullet starts at the BORE position, which is muzzle_height above ground
        // The sight is sight_height ABOVE the bore, so we don't add sight_height here
        let mut position = Vector3::new(
            0.0,
            self.inputs.muzzle_height, // Bore position above ground (NOT + sight_height)
            0.0,
        );

        // Calculate initial velocity components
        // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral (right)
        // Launch direction includes the aerodynamic-jump muzzle perturbation when enabled
        // (a no-op returning the bare muzzle/azimuth angles otherwise). MBA-959. Computed
        // once here and reused for the result so it isn't evaluated twice per solve.
        let aj_components = self.aerodynamic_jump_components();
        let (launch_elev, launch_azim) = self.launch_angles_from(aj_components.as_ref());
        let horizontal_velocity = self.inputs.muzzle_velocity * launch_elev.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * launch_azim.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * launch_elev.sin(), // Y: vertical component
            horizontal_velocity * launch_azim.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut dt = 0.001; // Initial step size
        let tolerance = 1e-6; // Error tolerance
        let safety_factor = 0.9; // Safety factor for step size adjustment
        let max_dt = 0.01; // Maximum step size
        let min_dt = 1e-6; // Minimum step size

        // Add a point counter to debug
        let mut iteration_count = 0;
        const MAX_ITERATIONS: usize = 100000;

        // Air density and wind are constant for the whole solve (self.atmosphere / self.wind
        // are immutable); compute once instead of every iteration (mirrors solve_rk4).
        let (air_density, speed_of_sound, resolved_temp_c, resolved_press_hpa) = self.resolved_atmosphere();
        // 0deg = headwind, 90deg = from the right (McCoy wind-FROM convention, matching
        // WindSock); wind enters drag via velocity - wind. Used when no segmented wind.
        let wind_vector = Vector3::new(
            -self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
            0.0,
            -self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
        );

        // Mach-transition distances for the sampled-output flags (see solve_euler/solve_rk4).
        let mut transonic_distances: Vec<f64> = Vec::new();
        let mut crossed_transonic = false;
        let mut crossed_subsonic = false;

        // Pitch-damping / precession diagnostics (MBA-966). Previously only the
        // Euler and fixed-RK4 solvers tracked these, so the default adaptive
        // RK45 path always reported null even with --enable-pitch-damping /
        // --enable-precession set. Mirror the RK4 tracking here.
        let mut min_pitch_damping = 1.0;
        let mut transonic_mach: Option<f64> = None;
        let pitch_coeffs = PitchDampingCoefficients::from_bullet_type(
            self.inputs.bullet_model.as_deref().unwrap_or("default"),
        );
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001,
                yaw_angle: 0.001,
                pitch_rate: 0.0,
                yaw_rate: 0.0,
                precession_angle: 0.0,
                nutation_phase: 0.0,
            })
        } else {
            None
        };
        let mut max_yaw_angle = 0.0;
        let mut max_precession_angle = 0.0;

        while position.x < self.max_range
            && position.y > self.inputs.ground_threshold
            && time < 100.0
        {
            // X is downrange
            iteration_count += 1;
            if iteration_count > MAX_ITERATIONS {
                break; // Prevent infinite loop
            }

            // Store current point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy = 0.5 * self.inputs.bullet_mass * velocity_magnitude.powi(2);

            points.push(TrajectoryPoint {
                time,
                position: position,
                velocity_magnitude,
                kinetic_energy,
            });

            // Record Mach-transition distances (constant sea-level speed of sound, matching the
            // transonic_mach tracking). Each threshold is recorded once, in descending order.
            {
                let mach_here = if speed_of_sound > 0.0 {
                    velocity_magnitude / speed_of_sound
                } else {
                    0.0
                };
                if !crossed_transonic && mach_here < 1.2 {
                    crossed_transonic = true;
                    transonic_distances.push(position.x);
                }
                if !crossed_subsonic && mach_here < 1.0 {
                    crossed_subsonic = true;
                    transonic_distances.push(position.x);
                }
            }

            if position.y > max_height {
                max_height = position.y;
            }

            // Pitch damping (RK45 solver) — track the minimum coefficient and the
            // Mach at which the projectile enters the transonic band (MBA-966).
            if self.inputs.enable_pitch_damping {
                let mach = velocity_magnitude / speed_of_sound;
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &pitch_coeffs);
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Precession / nutation (RK45 solver). Uses the step `dt` actually
            // taken for this iteration so the angular integration stays in sync
            // with the variable-step trajectory.
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let mach = velocity_magnitude / speed_of_sound;

                    let spin_rate_rad_s = if self.inputs.twist_rate > 0.0 {
                        let velocity_fps = velocity_magnitude * 3.28084;
                        let twist_rate_ft = self.inputs.twist_rate / 12.0;
                        (velocity_fps / twist_rate_ft) * 2.0 * std::f64::consts::PI
                    } else {
                        0.0
                    };

                    let params = PrecessionNutationParams {
                        mass_kg: self.inputs.bullet_mass,
                        caliber_m: self.inputs.bullet_diameter,
                        length_m: self.inputs.bullet_length,
                        spin_rate_rad_s,
                        spin_inertia: 6.94e-8,
                        transverse_inertia: 9.13e-7,
                        velocity_mps: velocity_magnitude,
                        air_density_kg_m3: air_density,
                        mach,
                        pitch_damping_coeff: -0.8,
                        nutation_damping_factor: 0.05,
                    };

                    *state = calculate_combined_angular_motion(&params, state, time, dt, 0.001);

                    if state.yaw_angle.abs() > max_yaw_angle {
                        max_yaw_angle = state.yaw_angle.abs();
                    }
                    if state.precession_angle.abs() > max_precession_angle {
                        max_precession_angle = state.precession_angle.abs();
                    }
                }
            }

            // RK45 step with adaptive step size (air_density / wind_vector hoisted above)
            let (new_pos, new_vel, new_dt) = self.rk45_step(
                &position,
                &velocity,
                dt,
                air_density,
                &wind_vector,
                tolerance,
                (speed_of_sound, resolved_temp_c, resolved_press_hpa),
            );

            // Advance state and time by the dt actually used for THIS step. (Previously dt
            // was overwritten with the adapted next-step size BEFORE `time += dt`, so every
            // reported time advanced by the NEXT step's dt — desyncing time from state and
            // corrupting time_of_flight and per-point / sampled times.)
            position = new_pos;
            velocity = new_vel;
            time += dt;

            // Adapt the step size for the NEXT iteration.
            dt = (safety_factor * new_dt).clamp(min_dt, max_dt);
        }

        // Ensure we have at least one point
        if points.is_empty() {
            return Err(BallisticsError::from("No trajectory points calculated"));
        }

        // Boundary interpolation to exactly max_range (MBA-968). The adaptive
        // loop stores the point at the TOP of each iteration, so the last stored
        // point sits one (possibly large) step SHORT of max_range while the
        // post-loop `position` has just overshot it. Without this, the default
        // RK45 solver reports ~2% short of --max-range, unlike the fixed-step
        // solvers. When the loop exited by crossing the range (not by hitting the
        // ground / time cap / iteration cap), append a linearly-interpolated
        // point at exactly max_range so the reported range matches the request.
        {
            let prev = points.last().unwrap().clone();
            let overshoot_x = position.x;
            let crossed_range = overshoot_x >= self.max_range && prev.position.x < self.max_range;
            if crossed_range {
                let span = overshoot_x - prev.position.x;
                if span > 1e-9 {
                    let frac = (self.max_range - prev.position.x) / span;
                    let interp_pos = prev.position + (position - prev.position) * frac;
                    let interp_vel_mag = prev.velocity_magnitude
                        + (velocity.magnitude() - prev.velocity_magnitude) * frac;
                    let interp_time = prev.time + (time - prev.time) * frac;
                    let interp_ke = 0.5 * self.inputs.bullet_mass * interp_vel_mag * interp_vel_mag;
                    points.push(TrajectoryPoint {
                        time: interp_time,
                        position: interp_pos,
                        velocity_magnitude: interp_vel_mag,
                        kinetic_energy: interp_ke,
                    });
                    if interp_pos.y > max_height {
                        max_height = interp_pos.y;
                    }
                }
            }
        }

        let last_point = points.last().unwrap();

        // Generate sampled trajectory points if enabled
        let sampled_points = if self.inputs.enable_trajectory_sampling {
            // Build trajectory data for sampling
            let trajectory_data = TrajectoryData {
                times: points.iter().map(|p| p.time).collect(),
                positions: points.iter().map(|p| p.position).collect(),
                velocities: points
                    .iter()
                    .map(|p| {
                        // Approximate velocity direction from position changes
                        Vector3::new(0.0, 0.0, p.velocity_magnitude)
                    })
                    .collect(),
                transonic_distances, // populated at each Mach-threshold crossing
            };

            // For LOS calculation in ground-referenced coordinates:
            // sight_position_m is the sight's actual y-position above ground
            // (muzzle_height + sight_height, not just sight_height)
            // For flat shots, target is at same height as the sight (horizontal LOS)
            let sight_position_m = self.inputs.muzzle_height + self.inputs.sight_height;
            let outputs = TrajectoryOutputs {
                target_distance_horiz_m: last_point.position.x,
                target_vertical_height_m: sight_position_m,
                time_of_flight_s: last_point.time,
                max_ord_dist_horiz_m: max_height,
                sight_height_m: sight_position_m,
            };

            let samples = sample_trajectory(
                &trajectory_data,
                &outputs,
                self.inputs.sample_interval,
                self.inputs.bullet_mass,
            );
            Some(samples)
        } else {
            None
        };

        Ok(TrajectoryResult {
            max_range: last_point.position.x, // X is downrange
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping {
                Some(min_pitch_damping)
            } else {
                None
            },
            transonic_mach,
            angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation {
                Some(max_yaw_angle)
            } else {
                None
            },
            max_precession_angle: if self.inputs.enable_precession_nutation {
                Some(max_precession_angle)
            } else {
                None
            },
            aerodynamic_jump: aj_components,
        })
    }

    fn rk45_step(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        dt: f64,
        air_density: f64,
        wind_vector: &Vector3<f64>,
        tolerance: f64,
        resolved_atmo: (f64, f64, f64), // (speed_of_sound, temp_c, press_hpa)
    ) -> (Vector3<f64>, Vector3<f64>, f64) {
        // Dormand-Prince coefficients
        const A21: f64 = 1.0 / 5.0;
        const A31: f64 = 3.0 / 40.0;
        const A32: f64 = 9.0 / 40.0;
        const A41: f64 = 44.0 / 45.0;
        const A42: f64 = -56.0 / 15.0;
        const A43: f64 = 32.0 / 9.0;
        const A51: f64 = 19372.0 / 6561.0;
        const A52: f64 = -25360.0 / 2187.0;
        const A53: f64 = 64448.0 / 6561.0;
        const A54: f64 = -212.0 / 729.0;
        const A61: f64 = 9017.0 / 3168.0;
        const A62: f64 = -355.0 / 33.0;
        const A63: f64 = 46732.0 / 5247.0;
        const A64: f64 = 49.0 / 176.0;
        const A65: f64 = -5103.0 / 18656.0;
        const A71: f64 = 35.0 / 384.0;
        const A73: f64 = 500.0 / 1113.0;
        const A74: f64 = 125.0 / 192.0;
        const A75: f64 = -2187.0 / 6784.0;
        const A76: f64 = 11.0 / 84.0;

        // 5th order coefficients
        const B1: f64 = 35.0 / 384.0;
        const B3: f64 = 500.0 / 1113.0;
        const B4: f64 = 125.0 / 192.0;
        const B5: f64 = -2187.0 / 6784.0;
        const B6: f64 = 11.0 / 84.0;

        // 4th order coefficients for error estimation
        const B1_ERR: f64 = 5179.0 / 57600.0;
        const B3_ERR: f64 = 7571.0 / 16695.0;
        const B4_ERR: f64 = 393.0 / 640.0;
        const B5_ERR: f64 = -92097.0 / 339200.0;
        const B6_ERR: f64 = 187.0 / 2100.0;
        const B7_ERR: f64 = 1.0 / 40.0;

        // Compute RK45 stages
        let k1_v = self.calculate_acceleration(position, velocity, air_density, wind_vector, resolved_atmo);
        let k1_p = *velocity;

        let p2 = position + dt * A21 * k1_p;
        let v2 = velocity + dt * A21 * k1_v;
        let k2_v = self.calculate_acceleration(&p2, &v2, air_density, wind_vector, resolved_atmo);
        let k2_p = v2;

        let p3 = position + dt * (A31 * k1_p + A32 * k2_p);
        let v3 = velocity + dt * (A31 * k1_v + A32 * k2_v);
        let k3_v = self.calculate_acceleration(&p3, &v3, air_density, wind_vector, resolved_atmo);
        let k3_p = v3;

        let p4 = position + dt * (A41 * k1_p + A42 * k2_p + A43 * k3_p);
        let v4 = velocity + dt * (A41 * k1_v + A42 * k2_v + A43 * k3_v);
        let k4_v = self.calculate_acceleration(&p4, &v4, air_density, wind_vector, resolved_atmo);
        let k4_p = v4;

        let p5 = position + dt * (A51 * k1_p + A52 * k2_p + A53 * k3_p + A54 * k4_p);
        let v5 = velocity + dt * (A51 * k1_v + A52 * k2_v + A53 * k3_v + A54 * k4_v);
        let k5_v = self.calculate_acceleration(&p5, &v5, air_density, wind_vector, resolved_atmo);
        let k5_p = v5;

        let p6 = position + dt * (A61 * k1_p + A62 * k2_p + A63 * k3_p + A64 * k4_p + A65 * k5_p);
        let v6 = velocity + dt * (A61 * k1_v + A62 * k2_v + A63 * k3_v + A64 * k4_v + A65 * k5_v);
        let k6_v = self.calculate_acceleration(&p6, &v6, air_density, wind_vector, resolved_atmo);
        let k6_p = v6;

        let p7 = position + dt * (A71 * k1_p + A73 * k3_p + A74 * k4_p + A75 * k5_p + A76 * k6_p);
        let v7 = velocity + dt * (A71 * k1_v + A73 * k3_v + A74 * k4_v + A75 * k5_v + A76 * k6_v);
        let k7_v = self.calculate_acceleration(&p7, &v7, air_density, wind_vector, resolved_atmo);
        let k7_p = v7;

        // 5th order solution
        let new_pos = position + dt * (B1 * k1_p + B3 * k3_p + B4 * k4_p + B5 * k5_p + B6 * k6_p);
        let new_vel = velocity + dt * (B1 * k1_v + B3 * k3_v + B4 * k4_v + B5 * k5_v + B6 * k6_v);

        // 4th order solution for error estimate
        let pos_err = position
            + dt * (B1_ERR * k1_p
                + B3_ERR * k3_p
                + B4_ERR * k4_p
                + B5_ERR * k5_p
                + B6_ERR * k6_p
                + B7_ERR * k7_p);
        let vel_err = velocity
            + dt * (B1_ERR * k1_v
                + B3_ERR * k3_v
                + B4_ERR * k4_v
                + B5_ERR * k5_v
                + B6_ERR * k6_v
                + B7_ERR * k7_v);

        // Estimate error
        let pos_error = (new_pos - pos_err).magnitude();
        let vel_error = (new_vel - vel_err).magnitude();
        let error = (pos_error + vel_error) / (1.0 + position.magnitude() + velocity.magnitude());

        // Calculate new step size
        let dt_new = if error < tolerance {
            dt * (tolerance / error).powf(0.2).min(2.0)
        } else {
            dt * (tolerance / error).powf(0.25).max(0.1)
        };

        (new_pos, new_vel, dt_new)
    }

    fn calculate_acceleration(
        &self,
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        air_density: f64,
        wind_vector: &Vector3<f64>,
        resolved_atmo: (f64, f64, f64), // (speed_of_sound, temp_c, press_hpa) hoisted per-solve
    ) -> Vector3<f64> {
        // Resolve the wind at this point. Downrange-segmented wind (when supplied)
        // takes precedence and is sampled by downrange distance (position.x) per
        // step; otherwise altitude-dependent shear (if enabled); otherwise the
        // constant `wind_vector`. Segmented wind is not combined with shear (the
        // CLI/WASM front-ends reject that combination), so the order is safe.
        let actual_wind = if let Some(ref sock) = self.wind_sock {
            sock.vector_for_range_stateless(position.x)
        } else if self.inputs.enable_wind_shear {
            self.get_wind_at_altitude(position.y)
        } else {
            *wind_vector
        };

        let relative_velocity = velocity - actual_wind;
        let velocity_magnitude = relative_velocity.magnitude();

        if velocity_magnitude < 0.001 {
            return self.gravity_acceleration();
        }

        // Get drag coefficient from drag model (Mach-indexed from drag tables)
        let cd = self.calculate_drag_coefficient(velocity_magnitude, resolved_atmo.0);

        // Convert velocity to fps for BC lookups
        let velocity_fps = velocity_magnitude * 3.28084;

        // Look up BC from segments if available (highest priority - most accurate)
        let base_bc = if let Some(ref segments) = self.inputs.bc_segments_data {
            // Find matching segment for current velocity
            segments
                .iter()
                .find(|seg| velocity_fps >= seg.velocity_min && velocity_fps < seg.velocity_max)
                .map(|seg| seg.bc_value)
                .unwrap_or(self.inputs.bc_value)
        } else {
            self.inputs.bc_value
        };

        // Apply cluster BC correction if enabled (on top of segment BC)
        let effective_bc = if let Some(ref cluster_bc) = self.cluster_bc {
            cluster_bc.apply_correction(
                base_bc,
                self.inputs.caliber_inches, // predict_cluster normalizes against an inches range
                self.inputs.weight_grains,
                velocity_fps,
            )
        } else {
            base_bc
        };
        // Guard bc_value == 0 (allowed on the FFI/WASM surfaces, which lack the CLI's 0.001
        // lower bound): dividing by effective_bc below would be Inf -> NaN. Inert for valid
        // BCs (>= 0.001).
        let effective_bc = effective_bc.max(1e-6);

        // Use proper ballistics retardation formula
        // This matches the proven formula from fast_trajectory.rs
        // The standard retardation factor converts Cd to drag deceleration
        // Note: velocity_fps already calculated above for BC segment lookup
        let cd_to_retard = crate::constants::CD_TO_RETARD;
        let standard_factor = cd * cd_to_retard;
        let density_scale = air_density / 1.225; // Scale relative to standard air (1.225 kg/m³)

        // Drag acceleration in ft/s² then convert to m/s²
        let a_drag_ft_s2 =
            (velocity_fps * velocity_fps) * standard_factor * density_scale / effective_bc;
        let a_drag_m_s2 = a_drag_ft_s2 * 0.3048; // ft/s² to m/s²

        // Apply drag opposite to velocity direction
        let drag_acceleration = -a_drag_m_s2 * (relative_velocity / velocity_magnitude);

        // Total acceleration = drag + gravity. `shooting_angle` rotates gravity into the shot
        // frame for inclined fire; at 0 deg this is the normal vertical-only gravity vector.
        let mut accel = drag_acceleration + self.gravity_acceleration();

        // Coriolis (Earth rotation). McCoy frame: X=downrange, Y=vertical, Z=lateral,
        // azimuth 0 = North. McCoy frame: X=downrange, Y=vertical, Z=lateral.
        if self.inputs.enable_coriolis {
            if let Some(lat_deg) = self.inputs.latitude {
                let omega_earth = 7.2921159e-5_f64; // rad/s
                let lat = lat_deg.to_radians();
                let az = self.inputs.shot_azimuth; // compass bearing (0=N), NOT the aiming offset
                                                   // Earth's angular velocity in the shot frame (X=downrange, Y=up,
                                                   // Z=lateral). Projecting Omega=(0, Ω cosφ, Ω sinφ) [local E,N,U] onto
                                                   // the azimuth-rotated shot axes gives a NEGATIVE lateral component:
                                                   // lateral = downrange × up points East for a North shot, and
                                                   // Omega·East = -Ω cosφ sin(az). The previous code dropped that sign.
                let omega = Vector3::new(
                    omega_earth * lat.cos() * az.cos(),  // X: downrange
                    omega_earth * lat.sin(),             // Y: vertical
                    -omega_earth * lat.cos() * az.sin(), // Z: lateral (MBA-938: corrected sign)
                );
                // Coriolis acceleration is the physical -2 Ω×v (MBA-938). The old +2 with
                // an "output-preserving relabel" justification produced left-ward drift for
                // a North shot in the Northern hemisphere; first principles (and the +Eötvös
                // lift for East shots) require -2 with the corrected omega above.
                accel += -2.0 * omega.cross(velocity);
            }
        }

        // Magnus side force (spinning projectile). SI units in this solver.
        if self.inputs.enable_magnus
            && self.inputs.bullet_diameter > 0.0
            && self.inputs.twist_rate > 0.0
        {
            let (_, spin_rad_s) =
                crate::spin_drift::calculate_spin_rate(velocity_magnitude, self.inputs.twist_rate);
            let (speed_of_sound, temp_c, press_hpa) = resolved_atmo;
            let temp_k = temp_c + 273.15;
            let mach = velocity_magnitude / speed_of_sound;

            // Imperial conversions for the stability / yaw-of-repose helpers.
            let d_in = self.inputs.bullet_diameter / 0.0254;
            let m_gr = self.inputs.bullet_mass / 0.00006479891;
            let l_in = if self.inputs.bullet_length > 0.0 {
                self.inputs.bullet_length / 0.0254
            } else {
                4.5 * d_in
            };
            // MBA-958: apply the canonical linear Miller density correction (T/T0)*(P0/P) to the
            // Magnus/yaw-of-repose Sg too, matching the spin-drift Sg (MBA-942) and stability.rs.
            // No-op at sea-level standard (15 C, 1013.25 hPa -> factor 1.0).
            let density_correction = if press_hpa > 0.0 && temp_k > 0.0 {
                (temp_k / 288.15) * (1013.25 / press_hpa)
            } else {
                1.0
            };
            let sg = crate::spin_drift::miller_stability(d_in, m_gr, self.inputs.twist_rate, l_in)
                * density_correction;

            // Yaw of repose (radians); zero for unstable bullets (Sg <= 1).
            let (yaw_rad, _) = crate::spin_drift::calculate_yaw_of_repose(
                sg,
                velocity_magnitude,
                spin_rad_s,
                0.0, // crosswind handled elsewhere
                0.0, // pitch rate not tracked
                air_density,
                d_in,
                l_in,
                m_gr,
                mach,
                "match",
                false,
            );

            // Proper McCoy Magnus FORCE: F = q S C_Npa (pd/2V) sin(alpha_R).
            let diameter_m = self.inputs.bullet_diameter; // already meters
            let spin_param = spin_rad_s * diameter_m / (2.0 * velocity_magnitude);
            let c_np = crate::derivatives::calculate_magnus_moment_coefficient(mach);
            let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);
            let magnus_force = 0.5
                * air_density
                * velocity_magnitude.powi(2)
                * area
                * c_np
                * spin_param
                * yaw_rad.sin();

            // Horizontal direction perpendicular to velocity. In McCoy (RH) frame,
            // v_unit × up = +Z (right) for a downrange shot, matching spin-drift sign.
            let velocity_unit = relative_velocity / velocity_magnitude;
            let up = Vector3::new(0.0, 1.0, 0.0);
            let mut dir = velocity_unit.cross(&up);
            let dir_norm = dir.norm();
            if dir_norm > 1e-12 && magnus_force.abs() > 1e-12 {
                dir /= dir_norm;
                if !self.inputs.is_twist_right {
                    dir = -dir;
                }
                accel += (magnus_force / self.inputs.bullet_mass) * dir;
            }
        }

        accel
    }

    fn calculate_drag_coefficient(&self, velocity: f64, speed_of_sound: f64) -> f64 {
        let mach = velocity / speed_of_sound;

        // MBA-940: a user-supplied custom drag table is the final Cd, used as-is — no G-model
        // lookup, no transonic shape correction, no form factor. The supplied curve already
        // encodes the projectile's true drag, so applying those would distort/double-count it.
        if let Some(ref table) = self.inputs.custom_drag_table {
            return table.interpolate(mach);
        }

        // Get drag coefficient from the drag tables (Mach-indexed)
        let base_cd = crate::drag::get_drag_coefficient(mach, &self.inputs.bc_type);

        // MBA-948: honor use_form_factor here too — previously only derivatives.rs applied it,
        // so cli_api and fast_trajectory silently ignored the flag. apply_form_factor_to_drag
        // short-circuits when the flag is false, so this is a no-op for every current consumer
        // (the flag is false on all CLI/FFI/WASM/binding surfaces and defaults false).
        crate::form_factor::apply_form_factor_to_drag(
            base_cd,
            self.inputs.bullet_model.as_deref(),
            &self.inputs.bc_type,
            self.inputs.use_form_factor,
        )
    }
}

// Monte Carlo parameters
#[derive(Debug, Clone)]
pub struct MonteCarloParams {
    pub num_simulations: usize,
    pub velocity_std_dev: f64,
    pub angle_std_dev: f64,
    pub bc_std_dev: f64,
    pub wind_speed_std_dev: f64,
    pub target_distance: Option<f64>,
    pub base_wind_speed: f64,
    pub base_wind_direction: f64,
    pub azimuth_std_dev: f64, // Horizontal aiming variation in radians
}

impl Default for MonteCarloParams {
    fn default() -> Self {
        Self {
            num_simulations: 1000,
            velocity_std_dev: 1.0,
            angle_std_dev: 0.001,
            bc_std_dev: 0.01,
            wind_speed_std_dev: 1.0,
            target_distance: None,
            base_wind_speed: 0.0,
            base_wind_direction: 0.0,
            azimuth_std_dev: 0.001, // Default horizontal spread ~0.057 degrees
        }
    }
}

// Monte Carlo results
#[derive(Debug, Clone)]
pub struct MonteCarloResults {
    pub ranges: Vec<f64>,
    pub impact_velocities: Vec<f64>,
    pub impact_positions: Vec<Vector3<f64>>,
}

/// Default hit-zone radius (meters) around the point of aim at the target plane — a 30 cm
/// circle. Shared by the CLI, FFI, and WASM so "hit probability" means the same thing everywhere.
pub const DEFAULT_HIT_RADIUS_M: f64 = 0.3;

impl MonteCarloResults {
    /// Fraction of simulations whose impact at the target plane lands within `hit_radius_m`
    /// of the point of aim. `impact_positions` are deviations from the baseline at the target
    /// plane (the downrange component is 0), so the vector norm is the radial miss distance.
    /// Samples that fall short of the target are clamped to their ground impact (a large
    /// deviation) and so correctly count as misses. Returns 0.0 when there are no samples.
    ///
    /// Single source of truth for hit probability — previously the CLI used a range-precision
    /// notion and the FFI a position notion with a redundant clause, so they disagreed.
    pub fn hit_probability(&self, hit_radius_m: f64) -> f64 {
        if self.impact_positions.is_empty() {
            return 0.0;
        }
        let hits = self
            .impact_positions
            .iter()
            .filter(|p| p.norm() < hit_radius_m)
            .count();
        hits as f64 / self.impact_positions.len() as f64
    }
}

// Run Monte Carlo simulation (backwards compatibility)
pub fn run_monte_carlo(
    base_inputs: BallisticInputs,
    params: MonteCarloParams,
) -> Result<MonteCarloResults, BallisticsError> {
    let base_wind = WindConditions {
        speed: params.base_wind_speed,
        direction: params.base_wind_direction,
    };
    run_monte_carlo_with_wind(base_inputs, base_wind, params)
}

// Run Monte Carlo simulation with wind
pub fn run_monte_carlo_with_wind(
    base_inputs: BallisticInputs,
    base_wind: WindConditions,
    params: MonteCarloParams,
) -> Result<MonteCarloResults, BallisticsError> {
    use rand_distr::{Distribution, Normal};

    let mut rng = rand::rng();
    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut impact_positions = Vec::new();

    let atmosphere = AtmosphericConditions {
        temperature: base_inputs.temperature,
        pressure: base_inputs.pressure,
        humidity: base_inputs.humidity_percent(),
        altitude: base_inputs.altitude,
    };
    let target_hint = params
        .target_distance
        .unwrap_or(base_inputs.target_distance);
    let solver_max_range = target_hint.max(1000.0) * 2.0;

    // First, calculate baseline trajectory with no variations
    let mut baseline_solver =
        TrajectorySolver::new(base_inputs.clone(), base_wind.clone(), atmosphere.clone());
    baseline_solver.set_max_range(solver_max_range);
    let baseline_result = baseline_solver.solve()?;

    // Determine target distance: use explicit target or baseline max range
    let target_distance = params.target_distance.unwrap_or(baseline_result.max_range);

    // Get baseline position at target distance (interpolated)
    let baseline_at_target = baseline_result
        .position_at_range(target_distance)
        .ok_or("Could not interpolate baseline at target distance")?;

    // Create normal distributions for variations
    let velocity_dist = Normal::new(base_inputs.muzzle_velocity, params.velocity_std_dev)
        .map_err(|e| format!("Invalid velocity distribution: {}", e))?;
    let angle_dist = Normal::new(base_inputs.muzzle_angle, params.angle_std_dev)
        .map_err(|e| format!("Invalid angle distribution: {}", e))?;
    let bc_dist = Normal::new(base_inputs.bc_value, params.bc_std_dev)
        .map_err(|e| format!("Invalid BC distribution: {}", e))?;
    let wind_speed_dist = Normal::new(base_wind.speed, params.wind_speed_std_dev)
        .map_err(|e| format!("Invalid wind speed distribution: {}", e))?;
    // MBA-952: wind-direction spread is APPROXIMATED from the wind-SPEED std dev (×0.1), a unit
    // conflation (m/s scaled as radians) — there is no dedicated wind_direction_std_dev field yet.
    // The dead WASM `--wind-dir-std` setter was removed (it set nothing). A proper fix is an
    // API-breaking wind_direction_std_dev on MonteCarloParams plumbed through WASM/FFI/main.
    let wind_dir_dist = Normal::new(base_wind.direction, params.wind_speed_std_dev * 0.1)
        .map_err(|e| format!("Invalid wind direction distribution: {}", e))?;
    let azimuth_dist = Normal::new(base_inputs.azimuth_angle, params.azimuth_std_dev)
        .map_err(|e| format!("Invalid azimuth distribution: {}", e))?;

    for _ in 0..params.num_simulations {
        // Create varied inputs
        let mut inputs = base_inputs.clone();
        inputs.muzzle_velocity = velocity_dist.sample(&mut rng).max(0.0);
        inputs.muzzle_angle = angle_dist.sample(&mut rng);
        inputs.bc_value = bc_dist.sample(&mut rng).max(0.01);
        inputs.azimuth_angle = azimuth_dist.sample(&mut rng); // Add horizontal variation

        // Create varied wind (now based on base wind conditions)
        let wind = WindConditions {
            speed: wind_speed_dist.sample(&mut rng).abs(),
            direction: wind_dir_dist.sample(&mut rng),
        };

        // Run trajectory
        let mut solver = TrajectorySolver::new(inputs, wind, atmosphere.clone());
        solver.set_max_range(solver_max_range);
        match solver.solve() {
            Ok(result) => {
                // MBA-967: do NOT skip samples that fall short of the target. range/velocity are
                // recorded at GROUND IMPACT for EVERY sample, so "Mean Range" is the ground-impact
                // distribution — independent of target_distance and consistent with `trajectory`.
                // All three result vectors still grow together per sample, so the equal-length FFI
                // ABI (exposed under one count) is preserved.
                let deviation = if result.max_range < target_distance {
                    // This sample never reached the target plane -> definite miss. Keep the
                    // encoded miss finite but far outside any practical target radius.
                    Vector3::new(0.0, -1.0e9, 0.0)
                } else {
                    let pos_at_target = match result.position_at_range(target_distance) {
                        Some(p) => p,
                        None => continue, // defensive: skip the whole sample (keeps vectors aligned)
                    };
                    // Deviation from baseline at the SAME target distance (McCoy): X = downrange
                    // (0 here), Y = vertical (elevation), Z = lateral (windage). Muzzle-angle
                    // sampling already models vertical pointing dispersion, so do not add a
                    // second independent vertical pointing draw here.
                    Vector3::new(
                        0.0,
                        pos_at_target.y - baseline_at_target.y,
                        pos_at_target.z - baseline_at_target.z,
                    )
                };

                ranges.push(result.max_range);
                impact_velocities.push(result.impact_velocity);
                impact_positions.push(deviation);
            }
            Err(_) => {
                // Skip failed simulations
                continue;
            }
        }
    }

    if ranges.is_empty() {
        return Err("No successful simulations".into());
    }

    Ok(MonteCarloResults {
        ranges,
        impact_velocities,
        impact_positions,
    })
}

// Calculate zero angle for a target
pub fn calculate_zero_angle(
    inputs: BallisticInputs,
    target_distance: f64,
    target_height: f64,
) -> Result<f64, BallisticsError> {
    calculate_zero_angle_with_conditions(
        inputs,
        target_distance,
        target_height,
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
}

pub fn calculate_zero_angle_with_conditions(
    inputs: BallisticInputs,
    target_distance: f64,
    target_height: f64,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
) -> Result<f64, BallisticsError> {
    // Helper function to get height at target distance for a given angle
    let get_height_at_angle = |angle: f64| -> Result<Option<f64>, BallisticsError> {
        let mut test_inputs = inputs.clone();
        test_inputs.muzzle_angle = angle;
        // MBA-959: zero on the bare bore. Aerodynamic jump is a constant elevation
        // offset, so leaving it on here would let the zero search silently absorb the
        // vertical jump. Disabling it makes AJ an additive POI shift relative to the
        // no-jump zero, regardless of the conditions the caller zeroes in.
        test_inputs.enable_aerodynamic_jump = false;

        let mut solver = TrajectorySolver::new(test_inputs, wind.clone(), atmosphere.clone());
        solver.set_max_range(target_distance * 2.0);
        solver.set_time_step(0.001);
        let result = solver.solve()?;

        // X is downrange in McCoy coordinates
        for i in 0..result.points.len() {
            if result.points[i].position.x >= target_distance {
                if i > 0 {
                    let p1 = &result.points[i - 1];
                    let p2 = &result.points[i];
                    let t = (target_distance - p1.position.x) / (p2.position.x - p1.position.x);
                    return Ok(Some(p1.position.y + t * (p2.position.y - p1.position.y)));
                } else {
                    return Ok(Some(result.points[i].position.y));
                }
            }
        }
        Ok(None)
    };

    // Binary search for the angle that hits the target
    // Use only positive angles to ensure proper ballistic arc (upward trajectory)
    let mut low_angle = 0.0; // radians (horizontal)
    let mut high_angle = 0.2; // radians (about 11 degrees)
    let tolerance = 0.00001; // radians
    let max_iterations = 50;

    // MBA-194: Validate bracketing before starting binary search
    // Check that the target height is actually between low and high angle trajectories
    let low_height = get_height_at_angle(low_angle)?;
    let high_height = get_height_at_angle(high_angle)?;

    match (low_height, high_height) {
        (Some(lh), Some(hh)) => {
            let low_error = lh - target_height;
            let high_error = hh - target_height;

            // For proper bracketing, low angle should undershoot (negative error)
            // and high angle should overshoot (positive error)
            if low_error > 0.0 && high_error > 0.0 {
                // Both angles overshoot - target is too close or height too low
                // This shouldn't happen for typical zeroing, but handle gracefully
                // Try to find a valid bracket by reducing low_angle (can't go negative)
                // Since we can't go below 0, just proceed and let binary search find best
            } else if low_error < 0.0 && high_error < 0.0 {
                // Both angles undershoot - target is beyond effective range
                // Try expanding high_angle up to 45 degrees (0.785 rad)
                let mut expanded = false;
                for multiplier in [2.0, 3.0, 4.0] {
                    let new_high = (high_angle * multiplier).min(0.785);
                    if let Ok(Some(h)) = get_height_at_angle(new_high) {
                        if h - target_height > 0.0 {
                            high_angle = new_high;
                            expanded = true;
                            break;
                        }
                    }
                    if new_high >= 0.785 {
                        break;
                    }
                }
                if !expanded {
                    return Err("Cannot find zero angle: target beyond effective range even at maximum angle".into());
                }
            }
            // If signs are opposite, we have valid bracketing - proceed
        }
        (None, Some(_hh)) => {
            // Low angle doesn't reach target, high does - this is fine
            // Binary search will increase low_angle until trajectory reaches
        }
        (Some(_lh), None) => {
            // High angle doesn't reach target - shouldn't happen
            return Err(
                "Cannot find zero angle: high angle trajectory doesn't reach target distance"
                    .into(),
            );
        }
        (None, None) => {
            // Neither reaches target - target too far
            return Err(
                "Cannot find zero angle: trajectory cannot reach target distance at any angle"
                    .into(),
            );
        }
    }

    for _iteration in 0..max_iterations {
        let mid_angle = (low_angle + high_angle) / 2.0;

        let mut test_inputs = inputs.clone();
        test_inputs.muzzle_angle = mid_angle;
        // MBA-959: zero on the bare bore so aerodynamic jump is not absorbed (see above).
        test_inputs.enable_aerodynamic_jump = false;

        let mut solver = TrajectorySolver::new(test_inputs, wind.clone(), atmosphere.clone());
        // Make sure we calculate far enough to reach the target
        solver.set_max_range(target_distance * 2.0);
        solver.set_time_step(0.001);
        let result = solver.solve()?;

        // Find the height at target distance (X is downrange)
        let mut height_at_target = None;
        for i in 0..result.points.len() {
            if result.points[i].position.x >= target_distance {
                if i > 0 {
                    // Linear interpolation
                    let p1 = &result.points[i - 1];
                    let p2 = &result.points[i];
                    let t = (target_distance - p1.position.x) / (p2.position.x - p1.position.x);
                    height_at_target = Some(p1.position.y + t * (p2.position.y - p1.position.y));
                } else {
                    height_at_target = Some(result.points[i].position.y);
                }
                break;
            }
        }

        match height_at_target {
            Some(height) => {
                let error = height - target_height;
                // MBA-193: Check height error FIRST (primary convergence criterion)
                // Height accuracy is what matters for zeroing - angle tolerance is secondary
                if error.abs() < 0.001 {
                    return Ok(mid_angle);
                }

                // Only use angle tolerance as convergence criterion if we have
                // exhausted angle precision AND height error is still acceptable
                // (within 10mm which is reasonable for long range)
                if (high_angle - low_angle).abs() < tolerance {
                    if error.abs() < 0.01 {
                        // Height error within 10mm - acceptable for practical use
                        return Ok(mid_angle);
                    }
                    // Angle bracket collapsed but the height error is still too large: the
                    // target is not actually reachable / was never bracketed. Returning
                    // Ok(mid_angle) here reported a NOT-zeroed angle as success (callers use
                    // it directly as muzzle_angle); surface it as an error instead.
                    return Err("Zero angle did not converge: residual height error too large (target not reachable / not bracketed)".into());
                }

                if error > 0.0 {
                    high_angle = mid_angle;
                } else {
                    low_angle = mid_angle;
                }
            }
            None => {
                // Trajectory didn't reach target distance, increase angle
                low_angle = mid_angle;

                // MBA-193: Check angle tolerance for None case too
                if (high_angle - low_angle).abs() < tolerance {
                    return Err("Trajectory cannot reach target distance - angle converged without valid solution".into());
                }
            }
        }
    }

    Err("Failed to find zero angle".into())
}

// Estimate BC from trajectory data
pub fn estimate_bc_from_trajectory(
    velocity: f64,
    mass: f64,
    diameter: f64,
    points: &[(f64, f64)], // (distance, drop) pairs
) -> Result<f64, BallisticsError> {
    // Simple BC estimation using least squares
    let mut best_bc = 0.5;
    let mut best_error = f64::MAX;
    let mut found_valid = false;

    // Try different BC values
    for bc in (100..1000).step_by(10) {
        let bc_value = bc as f64 / 1000.0;

        let inputs = BallisticInputs {
            muzzle_velocity: velocity,
            bc_value,
            bullet_mass: mass,
            bullet_diameter: diameter,
            ..Default::default()
        };

        let mut solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
        // Set max range for BC estimation
        solver.set_max_range(points.last().map(|(d, _)| *d * 1.5).unwrap_or(1000.0));

        let result = match solver.solve() {
            Ok(r) => r,
            Err(_) => continue, // Skip this BC value if solve fails
        };

        // Calculate error
        let mut total_error = 0.0;
        for (target_dist, target_drop) in points {
            // Find drop at this distance
            let mut calculated_drop = None;
            for i in 0..result.points.len() {
                if result.points[i].position.x >= *target_dist {
                    if i > 0 {
                        // Linear interpolation
                        let p1 = &result.points[i - 1];
                        let p2 = &result.points[i];
                        let t = (target_dist - p1.position.x) / (p2.position.x - p1.position.x);
                        calculated_drop =
                            Some(-(p1.position.y + t * (p2.position.y - p1.position.y)));
                    } else {
                        calculated_drop = Some(-result.points[i].position.y);
                    }
                    break;
                }
            }

            if let Some(drop) = calculated_drop {
                let error = (drop - target_drop).abs();
                total_error += error * error;
            }
        }

        if total_error < best_error {
            best_error = total_error;
            best_bc = bc_value;
            found_valid = true;
        }
    }

    if !found_valid {
        return Err(BallisticsError::from("Unable to estimate BC from provided data. Check that drop values are in correct units.".to_string()));
    }

    Ok(best_bc)
}

// Add rand dependencies for Monte Carlo
use rand;
use rand_distr;

#[cfg(test)]
mod ground_termination_tests {
    use super::*;

    // Regression lock for the unified ground termination: solve_euler/solve_rk4/solve_rk45 all
    // loop while `position.y > ground_threshold` (default -100.0), so they agree with RK45. A
    // lofted shot that returns to launch level before reaching max_range must keep descending to
    // the -100 m floor instead of stopping at y = 0 — and RK4-fixed and RK45 must behave the same.
    #[test]
    fn rk4_and_rk45_descend_to_ground_threshold() {
        for adaptive in [false, true] {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_angle = 0.1; // ~5.7 deg: arcs up, then descends past launch level
            inputs.use_rk4 = true;
            inputs.use_adaptive_rk45 = adaptive;
            assert_eq!(
                inputs.ground_threshold, -100.0,
                "default ground_threshold is -100 m"
            );

            let mut solver = TrajectorySolver::new(
                inputs,
                WindConditions::default(),
                AtmosphericConditions::default(),
            );
            // Huge max range: termination must be driven by ground_threshold, not the range cap.
            solver.set_max_range(1.0e7);

            let result = solver.solve().expect("solve should succeed");
            let final_y = result
                .points
                .last()
                .expect("trajectory has points")
                .position
                .y;
            assert!(
                final_y < -1.0,
                "adaptive_rk45={adaptive}: final y = {final_y} m; a lofted shot should descend \
                 past launch level toward the ground_threshold floor, not stop at y = 0"
            );
        }
    }
}

#[cfg(test)]
mod coriolis_direction_tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn transonic_crossing_flags_a_sampled_point() {
        // A supersonic shot that slows past Mach 1 must flag a sampled point as a Mach
        // transition. The underlying transonic_distances were a Vec::new() TODO, so this
        // flag was NEVER set regardless of trajectory — this is the regression guard.
        use crate::trajectory_sampling::TrajectoryFlag;
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 850.0; // ~2790 fps, well supersonic
        inputs.bc_value = 0.2; // low BC -> slows past Mach 1 within range
        inputs.bc_type = DragModel::G7;
        inputs.muzzle_angle = 0.03;
        inputs.enable_trajectory_sampling = true;
        inputs.sample_interval = 50.0;
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(2000.0);
        let r = solver.solve().expect("solve");
        let samples = r
            .sampled_points
            .expect("sampling enabled -> sampled_points present");
        assert!(
            samples
                .iter()
                .any(|s| s.flags.contains(&TrajectoryFlag::MachTransition)),
            "a shot that crosses Mach 1 must flag at least one Mach-transition sample"
        );
    }

    #[test]
    fn humidity_percent_converts_and_clamps() {
        // MBA-722: BallisticInputs.humidity is a 0-1 fraction; the helper yields 0-100 percent.
        let mut i = BallisticInputs::default();
        i.humidity = 0.5;
        assert!((i.humidity_percent() - 50.0).abs() < 1e-9, "0.5 -> 50%");
        i.humidity = 0.0;
        assert_eq!(i.humidity_percent(), 0.0);
        i.humidity = 1.0;
        assert_eq!(i.humidity_percent(), 100.0);
        i.humidity = 1.5; // out of range -> clamped, never > 100
        assert_eq!(i.humidity_percent(), 100.0);
    }

    /// Vertical position (m) at a given downrange `range_m`, for a shot fired along
    /// compass bearing `shot_azimuth` (radians, 0=N) with Coriolis enabled.
    fn vertical_at(shot_azimuth: f64, range_m: f64) -> f64 {
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 800.0;
        inputs.bc_value = 0.5;
        inputs.bc_type = DragModel::G7;
        inputs.muzzle_angle = 0.02; // ~20 mrad so it carries well past range_m
        inputs.enable_coriolis = true;
        inputs.latitude = Some(45.0);
        inputs.shot_azimuth = shot_azimuth;
        inputs.ground_threshold = f64::NEG_INFINITY; // never terminate early
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(range_m + 50.0);
        let r = solver.solve().expect("solve");
        let pts = &r.points;
        for i in 1..pts.len() {
            if pts[i].position.x >= range_m {
                let p1 = &pts[i - 1];
                let p2 = &pts[i];
                let t = (range_m - p1.position.x) / (p2.position.x - p1.position.x);
                return p1.position.y + t * (p2.position.y - p1.position.y);
            }
        }
        panic!("range {range_m} not reached");
    }

    /// Regression for the shot-direction Coriolis bug: the Eötvös vertical term
    /// `a_up = +2Ω cosφ v_east` lifts an EAST shot and depresses a WEST shot, so at a
    /// common range east must sit HIGHER than west, with north in between. Before the
    /// fix, `--shot-direction` never reached the solver and E/W/N were identical.
    #[test]
    fn eotvos_east_higher_than_west() {
        let range = 600.0;
        let east = vertical_at(FRAC_PI_2, range); // 90° E
        let west = vertical_at(3.0 * FRAC_PI_2, range); // 270° W
        let north = vertical_at(0.0, range); // 0° N
        assert!(
            east > west,
            "east ({east:.5}) must be higher than west ({west:.5}) at {range} m (Eötvös)"
        );
        assert!(
            east > north && north > west,
            "north ({north:.5}) must lie between east ({east:.5}) and west ({west:.5})"
        );
        assert!(
            (east - west) > 1e-3,
            "E-W vertical separation ({:.6} m) should be physically meaningful, not FP noise",
            east - west
        );
    }
}
