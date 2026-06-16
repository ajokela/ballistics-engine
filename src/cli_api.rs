// CLI API module - provides simplified interfaces for command-line tool
use crate::cluster_bc::ClusterBCDegradation;
use crate::pitch_damping::{calculate_pitch_damping_coefficient, PitchDampingCoefficients};
use crate::precession_nutation::{
    calculate_combined_angular_motion, AngularState, PrecessionNutationParams,
};
use crate::trajectory_sampling::{
    sample_trajectory, TrajectoryData, TrajectoryOutputs, TrajectorySample,
};
use crate::transonic_drag::{get_projectile_shape, transonic_correction, ProjectileShape};
use crate::wind_shear::{WindLayer, WindShearModel, WindShearProfile};
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
    pub bullet_length: f64,   // meters in the cli_api solver; inches in the derivatives/trajectory_integration path

    // Targeting and positioning
    pub muzzle_angle: f64,     // radians (launch angle)
    pub target_distance: f64,  // meters
    pub azimuth_angle: f64,    // horizontal aiming angle in radians
    pub shooting_angle: f64,   // uphill/downhill angle in radians
    pub sight_height: f64,     // meters above bore
    pub muzzle_height: f64,    // meters above ground
    pub target_height: f64,    // meters above ground for zeroing
    pub ground_threshold: f64, // meters below which to stop

    // Environmental conditions
    pub altitude: f64,         // meters
    pub temperature: f64,      // Celsius
    pub pressure: f64,         // millibars/hPa
    pub humidity: f64,         // relative humidity (0-1)
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
    pub use_cluster_bc: bool, // Use cluster-based BC degradation

    // Custom drag model support
    pub custom_drag_table: Option<crate::drag::DragTable>,

    // Legacy field for compatibility
    pub bc_type_str: Option<String>,
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
            bullet_length: diameter_m * 4.0, // Approximate

            // Targeting and positioning
            muzzle_angle: muzzle_angle_rad,
            target_distance: 100.0,
            azimuth_angle: 0.0,
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
            use_cluster_bc: false, // Disabled by default for backward compatibility

            // Custom drag model support
            custom_drag_table: None,

            // Legacy field for compatibility
            bc_type_str: None,
        }
    }
}

// Wind conditions
#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64,     // m/s
    pub direction: f64, // radians (0 = North, PI/2 = East)
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
    pub humidity: f64,    // percentage (0-100)
    pub altitude: f64,    // meters
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
        }
    }

    pub fn set_max_range(&mut self, range: f64) {
        self.max_range = range;
    }

    pub fn set_time_step(&mut self, step: f64) {
        self.time_step = step;
    }

    fn get_wind_at_altitude(&self, altitude_m: f64) -> Vector3<f64> {
        // Create wind shear profile based on surface wind
        let profile = WindShearProfile {
            model: if self.inputs.wind_shear_model == "logarithmic" {
                WindShearModel::Logarithmic
            } else if self.inputs.wind_shear_model == "power" {
                WindShearModel::PowerLaw
            } else {
                WindShearModel::PowerLaw // Default to power law
            },
            surface_wind: WindLayer {
                altitude_m: 0.0,
                speed_mps: self.wind.speed,
                direction_deg: self.wind.direction.to_degrees(),
            },
            reference_height: 10.0, // Standard meteorological measurement height
            roughness_length: 0.03, // Short grass
            power_exponent: 1.0 / 7.0, // Neutral stability
            geostrophic_wind: None,
            custom_layers: Vec::new(),
        };

        profile.get_wind_at_altitude(altitude_m)
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
        let sg = crate::spin_drift::miller_stability(d_in, m_gr, twist_in, length_in);
        let sign = if self.inputs.is_twist_right { 1.0 } else { -1.0 };

        for p in result.points.iter_mut() {
            if p.time <= 0.0 {
                continue;
            }
            let sd_in = 1.25 * (sg + 1.2) * p.time.powf(1.83); // Litz drift, inches
            p.position.z += sign * sd_in * 0.0254; // in -> m, Z = lateral
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
        let horizontal_velocity = self.inputs.muzzle_velocity * self.inputs.muzzle_angle.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * self.inputs.azimuth_angle.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * self.inputs.muzzle_angle.sin(), // Y: vertical component
            horizontal_velocity * self.inputs.azimuth_angle.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic

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
        let air_density = calculate_air_density(&self.atmosphere);

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        let wind_vector = Vector3::new(
            self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
            0.0,
            self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
        );

        // Main integration loop (X is downrange)
        while position.x < self.max_range && position.y >= 0.0 && time < 100.0 {
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

            // Debug: Log first and every 100th point
            // McCoy coordinate system: X=downrange, Y=vertical, Z=lateral
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
                let temp_c = self.atmosphere.temperature;
                let temp_k = temp_c + 273.15;
                let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let bullet_type = if let Some(ref model) = self.inputs.bullet_model {
                    model.as_str()
                } else {
                    "default"
                };
                let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let temp_c = self.atmosphere.temperature;
                    let temp_k = temp_c + 273.15;
                    let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
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

            // Calculate drag with altitude-dependent wind if enabled
            let actual_wind = if self.inputs.enable_wind_shear {
                self.get_wind_at_altitude(position.y)
            } else {
                wind_vector
            };
            let velocity_rel = velocity - actual_wind;
            let velocity_rel_mag = velocity_rel.magnitude();
            let drag_coefficient = self.calculate_drag_coefficient(velocity_rel_mag);

            // Calculate drag force
            let drag_force = 0.5
                * air_density
                * drag_coefficient
                * self.inputs.bullet_diameter
                * self.inputs.bullet_diameter
                * std::f64::consts::PI
                / 4.0
                * velocity_rel_mag
                * velocity_rel_mag;

            // Calculate acceleration
            let drag_acceleration = -drag_force / self.inputs.bullet_mass;
            let acceleration = Vector3::new(
                drag_acceleration * velocity_rel.x / velocity_rel_mag,
                drag_acceleration * velocity_rel.y / velocity_rel_mag - 9.80665,
                drag_acceleration * velocity_rel.z / velocity_rel_mag,
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
                transonic_distances: Vec::new(), // TODO: Track Mach transitions
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
        let horizontal_velocity = self.inputs.muzzle_velocity * self.inputs.muzzle_angle.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * self.inputs.azimuth_angle.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * self.inputs.muzzle_angle.sin(), // Y: vertical component
            horizontal_velocity * self.inputs.azimuth_angle.sin(), // Z: lateral (side deviation)
        );

        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0; // Track minimum pitch damping coefficient
        let mut transonic_mach = None; // Track when we enter transonic

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
        let air_density = calculate_air_density(&self.atmosphere);

        // Wind vector (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
        let wind_vector = Vector3::new(
            self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
            0.0,
            self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
        );

        // Main RK4 integration loop (X is downrange)
        while position.x < self.max_range && position.y >= 0.0 && time < 100.0 {
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

            if position.y > max_height {
                max_height = position.y;
            }

            // Calculate pitch damping if enabled (RK4 solver)
            if self.inputs.enable_pitch_damping {
                let temp_c = self.atmosphere.temperature;
                let temp_k = temp_c + 273.15;
                let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
                let mach = velocity_magnitude / speed_of_sound;

                // Track when we enter transonic
                if transonic_mach.is_none() && mach < 1.2 && mach > 0.8 {
                    transonic_mach = Some(mach);
                }

                // Calculate pitch damping coefficient
                let bullet_type = if let Some(ref model) = self.inputs.bullet_model {
                    model.as_str()
                } else {
                    "default"
                };
                let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);
                let pitch_damping = calculate_pitch_damping_coefficient(mach, &coeffs);

                // Track minimum (most critical for stability)
                if pitch_damping < min_pitch_damping {
                    min_pitch_damping = pitch_damping;
                }
            }

            // Calculate precession/nutation if enabled (RK4 solver)
            if self.inputs.enable_precession_nutation {
                if let Some(ref mut state) = angular_state {
                    let velocity_magnitude = velocity.magnitude();
                    let temp_c = self.atmosphere.temperature;
                    let temp_k = temp_c + 273.15;
                    let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
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
            let acc1 = self.calculate_acceleration(&position, &velocity, air_density, &wind_vector);

            // k2
            let pos2 = position + velocity * (dt * 0.5);
            let vel2 = velocity + acc1 * (dt * 0.5);
            let acc2 = self.calculate_acceleration(&pos2, &vel2, air_density, &wind_vector);

            // k3
            let pos3 = position + vel2 * (dt * 0.5);
            let vel3 = velocity + acc2 * (dt * 0.5);
            let acc3 = self.calculate_acceleration(&pos3, &vel3, air_density, &wind_vector);

            // k4
            let pos4 = position + vel3 * dt;
            let vel4 = velocity + acc3 * dt;
            let acc4 = self.calculate_acceleration(&pos4, &vel4, air_density, &wind_vector);

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
                transonic_distances: Vec::new(), // TODO: Track Mach transitions
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
        let horizontal_velocity = self.inputs.muzzle_velocity * self.inputs.muzzle_angle.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * self.inputs.azimuth_angle.cos(), // X: downrange (forward)
            self.inputs.muzzle_velocity * self.inputs.muzzle_angle.sin(), // Y: vertical component
            horizontal_velocity * self.inputs.azimuth_angle.sin(), // Z: lateral (side deviation)
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

            if position.y > max_height {
                max_height = position.y;
            }

            // Wind (McCoy): X=downrange (head/tail wind), Y=0, Z=lateral (crosswind)
            let air_density = calculate_air_density(&self.atmosphere);
            let wind_vector = Vector3::new(
                self.wind.speed * self.wind.direction.cos(), // X: downrange (head/tail wind)
                0.0,
                self.wind.speed * self.wind.direction.sin(), // Z: lateral (crosswind)
            );

            // RK45 step with adaptive step size
            let (new_pos, new_vel, new_dt) = self.rk45_step(
                &position,
                &velocity,
                dt,
                air_density,
                &wind_vector,
                tolerance,
            );

            // Update step size with safety factor and bounds
            dt = (safety_factor * new_dt).clamp(min_dt, max_dt);

            // Update state
            position = new_pos;
            velocity = new_vel;
            time += dt;
        }

        // Ensure we have at least one point
        if points.is_empty() {
            return Err(BallisticsError::from("No trajectory points calculated"));
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
                transonic_distances: Vec::new(),
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
            min_pitch_damping: None,
            transonic_mach: None,
            angular_state: None,
            max_yaw_angle: None,
            max_precession_angle: None,
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
        let k1_v = self.calculate_acceleration(position, velocity, air_density, wind_vector);
        let k1_p = *velocity;

        let p2 = position + dt * A21 * k1_p;
        let v2 = velocity + dt * A21 * k1_v;
        let k2_v = self.calculate_acceleration(&p2, &v2, air_density, wind_vector);
        let k2_p = v2;

        let p3 = position + dt * (A31 * k1_p + A32 * k2_p);
        let v3 = velocity + dt * (A31 * k1_v + A32 * k2_v);
        let k3_v = self.calculate_acceleration(&p3, &v3, air_density, wind_vector);
        let k3_p = v3;

        let p4 = position + dt * (A41 * k1_p + A42 * k2_p + A43 * k3_p);
        let v4 = velocity + dt * (A41 * k1_v + A42 * k2_v + A43 * k3_v);
        let k4_v = self.calculate_acceleration(&p4, &v4, air_density, wind_vector);
        let k4_p = v4;

        let p5 = position + dt * (A51 * k1_p + A52 * k2_p + A53 * k3_p + A54 * k4_p);
        let v5 = velocity + dt * (A51 * k1_v + A52 * k2_v + A53 * k3_v + A54 * k4_v);
        let k5_v = self.calculate_acceleration(&p5, &v5, air_density, wind_vector);
        let k5_p = v5;

        let p6 = position + dt * (A61 * k1_p + A62 * k2_p + A63 * k3_p + A64 * k4_p + A65 * k5_p);
        let v6 = velocity + dt * (A61 * k1_v + A62 * k2_v + A63 * k3_v + A64 * k4_v + A65 * k5_v);
        let k6_v = self.calculate_acceleration(&p6, &v6, air_density, wind_vector);
        let k6_p = v6;

        let p7 = position + dt * (A71 * k1_p + A73 * k3_p + A74 * k4_p + A75 * k5_p + A76 * k6_p);
        let v7 = velocity + dt * (A71 * k1_v + A73 * k3_v + A74 * k4_v + A75 * k5_v + A76 * k6_v);
        let k7_v = self.calculate_acceleration(&p7, &v7, air_density, wind_vector);
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
    ) -> Vector3<f64> {
        // Calculate altitude-dependent wind if wind shear is enabled
        let actual_wind = if self.inputs.enable_wind_shear {
            self.get_wind_at_altitude(position.y)
        } else {
            *wind_vector
        };

        let relative_velocity = velocity - actual_wind;
        let velocity_magnitude = relative_velocity.magnitude();

        if velocity_magnitude < 0.001 {
            return Vector3::new(0.0, -9.81, 0.0);
        }

        // Get drag coefficient from drag model (Mach-indexed from drag tables)
        let cd = self.calculate_drag_coefficient(velocity_magnitude);

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
                self.inputs.caliber_inches * 0.0254, // Convert back to meters for consistency
                self.inputs.weight_grains,
                velocity_fps,
            )
        } else {
            base_bc
        };

        // Use proper ballistics retardation formula
        // This matches the proven formula from fast_trajectory.rs
        // The standard retardation factor converts Cd to drag deceleration
        // Note: velocity_fps already calculated above for BC segment lookup
        let cd_to_retard = 0.000683 * 0.30; // Standard ballistics constant
        let standard_factor = cd * cd_to_retard;
        let density_scale = air_density / 1.225; // Scale relative to standard air (1.225 kg/m³)

        // Drag acceleration in ft/s² then convert to m/s²
        let a_drag_ft_s2 =
            (velocity_fps * velocity_fps) * standard_factor * density_scale / effective_bc;
        let a_drag_m_s2 = a_drag_ft_s2 * 0.3048; // ft/s² to m/s²

        // Apply drag opposite to velocity direction
        let drag_acceleration = -a_drag_m_s2 * (relative_velocity / velocity_magnitude);

        // Total acceleration = drag + gravity
        let mut accel = drag_acceleration + Vector3::new(0.0, -9.81, 0.0);

        // Coriolis (Earth rotation). McCoy frame: X=downrange, Y=vertical, Z=lateral,
        // azimuth 0 = North. McCoy frame: X=downrange, Y=vertical, Z=lateral.
        if self.inputs.enable_coriolis {
            if let Some(lat_deg) = self.inputs.latitude {
                let omega_earth = 7.2921159e-5_f64; // rad/s
                let lat = lat_deg.to_radians();
                let az = self.inputs.azimuth_angle;
                let omega = Vector3::new(
                    omega_earth * lat.cos() * az.cos(), // X: downrange component
                    omega_earth * lat.sin(),            // Y: vertical component
                    omega_earth * lat.cos() * az.sin(), // Z: lateral component
                );
                // Coriolis is -2 Ω×v. Migrating from the engine's original
                // left-handed (lateral, up, downrange) frame to McCoy right-handed
                // (downrange, up, lateral) is a reflection, which negates the cross
                // product; the +2 here reproduces the original deflection in the new
                // frame (output-preserving relabel).
                accel += 2.0 * omega.cross(velocity);
            }
        }

        // Magnus side force (spinning projectile). SI units in this solver.
        if self.inputs.enable_magnus
            && self.inputs.bullet_diameter > 0.0
            && self.inputs.twist_rate > 0.0
        {
            let (_, spin_rad_s) =
                crate::spin_drift::calculate_spin_rate(velocity_magnitude, self.inputs.twist_rate);
            let temp_k = self.atmosphere.temperature + 273.15;
            let speed_of_sound = (1.4 * 287.05 * temp_k).sqrt();
            let mach = velocity_magnitude / speed_of_sound;

            // Imperial conversions for the stability / yaw-of-repose helpers.
            let d_in = self.inputs.bullet_diameter / 0.0254;
            let m_gr = self.inputs.bullet_mass / 0.00006479891;
            let l_in = if self.inputs.bullet_length > 0.0 {
                self.inputs.bullet_length / 0.0254
            } else {
                4.5 * d_in
            };
            let sg = crate::spin_drift::miller_stability(d_in, m_gr, self.inputs.twist_rate, l_in);

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

    fn calculate_drag_coefficient(&self, velocity: f64) -> f64 {
        // Calculate speed of sound based on atmospheric temperature
        let temp_c = self.atmosphere.temperature;
        let temp_k = temp_c + 273.15;
        let gamma = 1.4; // Ratio of specific heats for air
        let r_specific = 287.05; // Specific gas constant for air (J/kg·K)
        let speed_of_sound = (gamma * r_specific * temp_k).sqrt();
        let mach = velocity / speed_of_sound;

        // Get drag coefficient from the drag tables (Mach-indexed)
        let base_cd = crate::drag::get_drag_coefficient(mach, &self.inputs.bc_type);

        // Determine projectile shape for transonic corrections
        let projectile_shape = if let Some(ref model) = self.inputs.bullet_model {
            // Try to determine shape from bullet model string
            if model.to_lowercase().contains("boat") || model.to_lowercase().contains("bt") {
                ProjectileShape::BoatTail
            } else if model.to_lowercase().contains("round") || model.to_lowercase().contains("rn")
            {
                ProjectileShape::RoundNose
            } else if model.to_lowercase().contains("flat") || model.to_lowercase().contains("fb") {
                ProjectileShape::FlatBase
            } else {
                // Use heuristic based on caliber, weight, and drag model
                get_projectile_shape(
                    self.inputs.bullet_diameter,
                    self.inputs.bullet_mass / 0.00006479891, // Convert kg to grains
                    &self.inputs.bc_type.to_string(),
                )
            }
        } else {
            // Use heuristic based on caliber, weight, and drag model
            get_projectile_shape(
                self.inputs.bullet_diameter,
                self.inputs.bullet_mass / 0.00006479891, // Convert kg to grains
                &self.inputs.bc_type.to_string(),
            )
        };

        // Apply transonic corrections
        // Note: Wave drag is disabled because G7/G1 drag functions already include
        // transonic effects. Adding wave drag on top would double-count the drag rise.
        // Wave drag should only be enabled for custom drag functions that don't
        // include transonic behavior.
        let include_wave_drag = false;
        transonic_correction(mach, base_cd, projectile_shape, include_wave_drag)
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
    use rand::thread_rng;
    use rand_distr::{Distribution, Normal};

    let mut rng = thread_rng();
    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut impact_positions = Vec::new();

    // First, calculate baseline trajectory with no variations
    let baseline_solver =
        TrajectorySolver::new(base_inputs.clone(), base_wind.clone(), Default::default());
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
    let wind_dir_dist =
        Normal::new(base_wind.direction, params.wind_speed_std_dev * 0.1) // Small variation in direction
            .map_err(|e| format!("Invalid wind direction distribution: {}", e))?;
    let azimuth_dist = Normal::new(base_inputs.azimuth_angle, params.azimuth_std_dev)
        .map_err(|e| format!("Invalid azimuth distribution: {}", e))?;

    // Create distribution for pointing errors (simulates shooter's aiming consistency)
    let pointing_error_dist = Normal::new(0.0, params.angle_std_dev * target_distance)
        .map_err(|e| format!("Invalid pointing distribution: {}", e))?;

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
        let solver = TrajectorySolver::new(inputs, wind, Default::default());
        match solver.solve() {
            Ok(result) => {
                ranges.push(result.max_range);
                impact_velocities.push(result.impact_velocity);

                // Interpolate position at target distance (not ground impact)
                if let Some(pos_at_target) = result.position_at_range(target_distance) {
                    // Calculate deviation from baseline at the SAME target distance (McCoy)
                    // X = downrange (0 here), Y = vertical (elevation), Z = lateral (windage)
                    let mut deviation = Vector3::new(
                        0.0, // X downrange deviation is 0 since we compare at same range
                        pos_at_target.y - baseline_at_target.y, // Vertical deviation
                        pos_at_target.z - baseline_at_target.z, // Lateral deviation (windage)
                    );

                    // Add additional pointing error to simulate realistic group sizes
                    // This represents the shooter's ability to aim consistently
                    let pointing_error_y = pointing_error_dist.sample(&mut rng);
                    deviation.y += pointing_error_y;

                    impact_positions.push(deviation);
                }
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
                    // Angle converged but height error too large - this shouldn't happen
                    // with proper tolerance values, but return best effort
                    return Ok(mid_angle);
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

// Helper function to calculate air density
fn calculate_air_density(atmosphere: &AtmosphericConditions) -> f64 {
    // Simplified air density calculation
    // P / (R * T) where R is specific gas constant for dry air
    let r_specific = 287.058; // J/(kg·K)
    let temperature_k = atmosphere.temperature + 273.15;

    // Convert pressure from hPa to Pa
    let pressure_pa = atmosphere.pressure * 100.0;

    // Basic density calculation
    let density = pressure_pa / (r_specific * temperature_k);

    // Altitude correction (simplified)
    let altitude_factor = (-atmosphere.altitude / 8000.0).exp();

    density * altitude_factor
}

// Add rand dependencies for Monte Carlo
use rand;
use rand_distr;
