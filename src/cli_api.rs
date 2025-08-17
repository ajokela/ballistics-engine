// CLI API module - provides simplified interfaces for command-line tool
use crate::DragModel;
use crate::wind_shear::{WindShearProfile, WindShearModel, WindLayer};
use crate::transonic_drag::{transonic_correction, get_projectile_shape, ProjectileShape};
use crate::trajectory_sampling::{sample_trajectory, TrajectoryData, TrajectoryOutputs, TrajectorySample};
use crate::pitch_damping::{calculate_pitch_damping_coefficient, PitchDampingCoefficients};
use crate::precession_nutation::{AngularState, PrecessionNutationParams, calculate_combined_angular_motion};
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
        BallisticsError { message: msg.to_string() }
    }
}

// Ballistic input parameters
#[derive(Debug, Clone)]
pub struct BallisticInputs {
    // Core ballistics parameters
    pub muzzle_velocity: f64,      // m/s
    pub launch_angle: f64,          // radians (same as muzzle_angle)
    pub ballistic_coefficient: f64,
    pub mass: f64,                  // kg
    pub diameter: f64,              // meters
    pub drag_model: DragModel,
    pub sight_height: f64,          // meters
    
    // Additional fields for compatibility
    pub altitude: f64,              // meters
    pub bc_type: DragModel,         // same as drag_model
    pub bc_value: f64,              // same as ballistic_coefficient
    pub caliber_inches: f64,        // diameter in inches
    pub weight_grains: f64,         // mass in grains
    pub bullet_diameter: f64,       // same as diameter
    pub bullet_mass: f64,           // same as mass
    pub bullet_length: f64,         // meters
    pub muzzle_angle: f64,          // same as launch_angle
    pub target_distance: f64,       // meters
    pub azimuth_angle: f64,         // horizontal aiming angle in radians
    pub use_rk4: bool,              // Use RK4 integration instead of Euler
    pub temperature: f64,           // Celsius
    pub twist_rate: f64,            // inches per turn
    pub is_twist_right: bool,       // right-hand twist
    pub shooting_angle: f64,        // uphill/downhill angle in radians
    pub latitude: Option<f64>,      // degrees
    pub ground_threshold: f64,      // meters below which to stop
    
    // Advanced effects flags
    pub enable_advanced_effects: bool,
    pub use_powder_sensitivity: bool,
    pub powder_temp_sensitivity: f64,
    pub powder_temp: f64,          // Celsius
    pub tipoff_yaw: f64,            // radians
    pub tipoff_decay_distance: f64, // meters
    pub use_bc_segments: bool,
    pub bc_segments: Option<Vec<(f64, f64)>>,  // Mach-BC pairs
    pub bc_segments_data: Option<Vec<crate::BCSegmentData>>,  // Velocity-BC segments
    pub use_enhanced_spin_drift: bool,
    pub use_form_factor: bool,
    pub enable_wind_shear: bool,
    pub wind_shear_model: String,
    pub enable_trajectory_sampling: bool,
    pub sample_interval: f64,  // meters
    pub enable_pitch_damping: bool,
    pub enable_precession_nutation: bool,
    
    // Additional data fields
    pub bc_type_str: Option<String>,
    pub bullet_model: Option<String>,
    pub bullet_id: Option<String>,
}

impl Default for BallisticInputs {
    fn default() -> Self {
        let mass_kg = 0.01;
        let diameter_m = 0.00762;
        let bc = 0.5;
        let launch_angle_rad = 0.0;
        let drag_model = DragModel::G1;
        
        Self {
            // Core parameters
            muzzle_velocity: 800.0,
            launch_angle: launch_angle_rad,
            ballistic_coefficient: bc,
            mass: mass_kg,
            diameter: diameter_m,
            drag_model,
            sight_height: 0.05,
            
            // Duplicate fields for compatibility
            altitude: 0.0,
            bc_type: drag_model,
            bc_value: bc,
            caliber_inches: diameter_m / 0.0254,  // Convert to inches
            weight_grains: mass_kg / 0.00006479891,  // Convert to grains
            bullet_diameter: diameter_m,
            bullet_mass: mass_kg,
            bullet_length: diameter_m * 4.0,  // Approximate
            muzzle_angle: launch_angle_rad,
            target_distance: 100.0,
            azimuth_angle: 0.0,
            use_rk4: true,  // Default to RK4 for better accuracy
            temperature: 15.0,
            twist_rate: 12.0,  // 1:12" typical
            is_twist_right: true,
            shooting_angle: 0.0,
            latitude: None,
            ground_threshold: -10.0,
            
            // Advanced effects (disabled by default)
            enable_advanced_effects: false,
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
            sample_interval: 10.0,  // Default 10 meter intervals
            enable_pitch_damping: false,
            enable_precession_nutation: false,
            
            // Optional data
            bc_type_str: None,
            bullet_model: None,
            bullet_id: None,
        }
    }
}

// Wind conditions
#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64,        // m/s
    pub direction: f64,    // radians (0 = North, PI/2 = East)
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
    pub temperature: f64,  // Celsius
    pub pressure: f64,     // hPa
    pub humidity: f64,     // percentage (0-100)
    pub altitude: f64,     // meters
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
    pub sampled_points: Option<Vec<TrajectorySample>>,  // Trajectory samples at regular intervals
    pub min_pitch_damping: Option<f64>,  // Minimum pitch damping coefficient (for stability warning)
    pub transonic_mach: Option<f64>,      // Mach number when entering transonic regime
    pub angular_state: Option<AngularState>,  // Final angular state if precession/nutation enabled
    pub max_yaw_angle: Option<f64>,           // Maximum yaw angle during flight (radians)
    pub max_precession_angle: Option<f64>,    // Maximum precession angle (radians)
}

// Trajectory solver
pub struct TrajectorySolver {
    inputs: BallisticInputs,
    wind: WindConditions,
    atmosphere: AtmosphericConditions,
    max_range: f64,
    time_step: f64,
}

impl TrajectorySolver {
    pub fn new(mut inputs: BallisticInputs, wind: WindConditions, atmosphere: AtmosphericConditions) -> Self {
        // Ensure duplicate fields are synchronized
        inputs.bc_type = inputs.drag_model;
        inputs.bc_value = inputs.ballistic_coefficient;
        inputs.bullet_diameter = inputs.diameter;
        inputs.bullet_mass = inputs.mass;
        inputs.muzzle_angle = inputs.launch_angle;
        inputs.caliber_inches = inputs.diameter / 0.0254;
        inputs.weight_grains = inputs.mass / 0.00006479891;
        
        Self {
            inputs,
            wind,
            atmosphere,
            max_range: 1000.0,
            time_step: 0.001,
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
                WindShearModel::PowerLaw  // Default to power law
            },
            surface_wind: WindLayer {
                altitude_m: 0.0,
                speed_mps: self.wind.speed,
                direction_deg: self.wind.direction.to_degrees(),
            },
            reference_height: 10.0,  // Standard meteorological measurement height
            roughness_length: 0.03,  // Short grass
            power_exponent: 1.0 / 7.0,  // Neutral stability
            geostrophic_wind: None,
            custom_layers: Vec::new(),
        };
        
        profile.get_wind_at_altitude(altitude_m)
    }
    
    pub fn solve(&self) -> Result<TrajectoryResult, BallisticsError> {
        if self.inputs.use_rk4 {
            self.solve_rk4()
        } else {
            self.solve_euler()
        }
    }
    
    fn solve_euler(&self) -> Result<TrajectoryResult, BallisticsError> {
        // Simple trajectory integration using Euler method
        let mut time = 0.0;
        let mut position = Vector3::new(0.0, self.inputs.sight_height, 0.0);
        // Calculate initial velocity components with both elevation and azimuth
        let horizontal_velocity = self.inputs.muzzle_velocity * self.inputs.launch_angle.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * self.inputs.azimuth_angle.sin(),  // X: side deviation (left/right)
            self.inputs.muzzle_velocity * self.inputs.launch_angle.sin(),  // Y: vertical component
            horizontal_velocity * self.inputs.azimuth_angle.cos(),  // Z: forward component (down-range)
        );
        
        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0;  // Track minimum pitch damping coefficient
        let mut transonic_mach = None;    // Track when we enter transonic
        
        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001,  // Small initial disturbance
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
        
        // Wind vector
        let wind_vector = Vector3::new(
            self.wind.speed * self.wind.direction.sin(),
            0.0,
            self.wind.speed * self.wind.direction.cos(),
        );
        
        // Main integration loop
        while position.z < self.max_range && position.y >= 0.0 && time < 100.0 {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy = 0.5 * self.inputs.mass * velocity_magnitude * velocity_magnitude;
            
            points.push(TrajectoryPoint {
                time,
                position: position.clone(),
                velocity_magnitude,
                kinetic_energy,
            });
            
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
                        mass_kg: self.inputs.mass,
                        caliber_m: self.inputs.diameter,
                        length_m: self.inputs.bullet_length,
                        spin_rate_rad_s,
                        spin_inertia: 6.94e-8,      // Typical value
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
                        0.001,  // Initial disturbance
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
                wind_vector.clone()
            };
            let velocity_rel = velocity - actual_wind;
            let velocity_rel_mag = velocity_rel.magnitude();
            let drag_coefficient = self.calculate_drag_coefficient(velocity_rel_mag);
            
            // Calculate drag force
            let drag_force = 0.5 * air_density * drag_coefficient * 
                            self.inputs.diameter * self.inputs.diameter * 
                            std::f64::consts::PI / 4.0 * velocity_rel_mag * velocity_rel_mag;
            
            // Calculate acceleration
            let drag_acceleration = -drag_force / self.inputs.mass;
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
                positions: points.iter().map(|p| p.position.clone()).collect(),
                velocities: points.iter().map(|p| {
                    // Reconstruct velocity vectors from magnitude (approximate)
                    Vector3::new(0.0, 0.0, p.velocity_magnitude)
                }).collect(),
                transonic_distances: Vec::new(),  // TODO: Track Mach transitions
            };
            
            let outputs = TrajectoryOutputs {
                target_distance_horiz_m: last_point.position.z,
                target_vertical_height_m: last_point.position.y,
                time_of_flight_s: last_point.time,
                max_ord_dist_horiz_m: max_height,
            };
            
            // Sample at specified intervals
            let samples = sample_trajectory(&trajectory_data, &outputs, self.inputs.sample_interval, self.inputs.mass);
            Some(samples)
        } else {
            None
        };
        
        Ok(TrajectoryResult {
            max_range: last_point.position.z,
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping { Some(min_pitch_damping) } else { None },
            transonic_mach: transonic_mach,
            angular_state: angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation { Some(max_yaw_angle) } else { None },
            max_precession_angle: if self.inputs.enable_precession_nutation { Some(max_precession_angle) } else { None },
        })
    }
    
    fn solve_rk4(&self) -> Result<TrajectoryResult, BallisticsError> {
        // RK4 trajectory integration for better accuracy
        let mut time = 0.0;
        let mut position = Vector3::new(0.0, self.inputs.sight_height, 0.0);
        
        // Calculate initial velocity components with both elevation and azimuth
        let horizontal_velocity = self.inputs.muzzle_velocity * self.inputs.launch_angle.cos();
        let mut velocity = Vector3::new(
            horizontal_velocity * self.inputs.azimuth_angle.sin(),
            self.inputs.muzzle_velocity * self.inputs.launch_angle.sin(),
            horizontal_velocity * self.inputs.azimuth_angle.cos(),
        );
        
        let mut points = Vec::new();
        let mut max_height = position.y;
        let mut min_pitch_damping = 1.0;  // Track minimum pitch damping coefficient
        let mut transonic_mach = None;    // Track when we enter transonic
        
        // Initialize angular state for precession/nutation tracking
        let mut angular_state = if self.inputs.enable_precession_nutation {
            Some(AngularState {
                pitch_angle: 0.001,  // Small initial disturbance
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
        
        // Wind vector
        let wind_vector = Vector3::new(
            self.wind.speed * self.wind.direction.sin(),
            0.0,
            self.wind.speed * self.wind.direction.cos(),
        );
        
        // Main RK4 integration loop
        while position.z < self.max_range && position.y >= 0.0 && time < 100.0 {
            // Store trajectory point
            let velocity_magnitude = velocity.magnitude();
            let kinetic_energy = 0.5 * self.inputs.mass * velocity_magnitude * velocity_magnitude;
            
            points.push(TrajectoryPoint {
                time,
                position: position.clone(),
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
                        mass_kg: self.inputs.mass,
                        caliber_m: self.inputs.diameter,
                        length_m: self.inputs.bullet_length,
                        spin_rate_rad_s,
                        spin_inertia: 6.94e-8,      // Typical value
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
                        0.001,  // Initial disturbance
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
                positions: points.iter().map(|p| p.position.clone()).collect(),
                velocities: points.iter().map(|p| {
                    // Reconstruct velocity vectors from magnitude (approximate)
                    Vector3::new(0.0, 0.0, p.velocity_magnitude)
                }).collect(),
                transonic_distances: Vec::new(),  // TODO: Track Mach transitions
            };
            
            let outputs = TrajectoryOutputs {
                target_distance_horiz_m: last_point.position.z,
                target_vertical_height_m: last_point.position.y,
                time_of_flight_s: last_point.time,
                max_ord_dist_horiz_m: max_height,
            };
            
            // Sample at specified intervals
            let samples = sample_trajectory(&trajectory_data, &outputs, self.inputs.sample_interval, self.inputs.mass);
            Some(samples)
        } else {
            None
        };
        
        Ok(TrajectoryResult {
            max_range: last_point.position.z,
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
            sampled_points,
            min_pitch_damping: if self.inputs.enable_pitch_damping { Some(min_pitch_damping) } else { None },
            transonic_mach: transonic_mach,
            angular_state: angular_state,
            max_yaw_angle: if self.inputs.enable_precession_nutation { Some(max_yaw_angle) } else { None },
            max_precession_angle: if self.inputs.enable_precession_nutation { Some(max_precession_angle) } else { None },
        })
    }
    
    fn calculate_acceleration(&self, position: &Vector3<f64>, velocity: &Vector3<f64>, air_density: f64, wind_vector: &Vector3<f64>) -> Vector3<f64> {
        // Calculate altitude-dependent wind if wind shear is enabled
        let actual_wind = if self.inputs.enable_wind_shear {
            self.get_wind_at_altitude(position.y)
        } else {
            wind_vector.clone()
        };
        
        let relative_velocity = velocity - &actual_wind;
        let velocity_magnitude = relative_velocity.magnitude();
        
        if velocity_magnitude < 0.001 {
            return Vector3::new(0.0, -9.81, 0.0);
        }
        
        // Calculate drag force
        let cd = self.calculate_drag_coefficient(velocity_magnitude);
        let reference_area = std::f64::consts::PI * (self.inputs.diameter / 2.0).powi(2);
        let drag_magnitude = 0.5 * air_density * velocity_magnitude.powi(2) * cd * reference_area / self.inputs.ballistic_coefficient;
        
        // Drag acts opposite to velocity
        let drag_force = -relative_velocity.normalize() * drag_magnitude;
        
        // Total acceleration = drag/mass + gravity
        let acceleration = drag_force / self.inputs.mass + Vector3::new(0.0, -9.81, 0.0);
        
        acceleration
    }
    
    fn calculate_drag_coefficient(&self, velocity: f64) -> f64 {
        // Calculate speed of sound based on atmospheric temperature
        // Use standard atmosphere temperature at sea level if not available
        let temp_c = self.atmosphere.temperature;
        let temp_k = temp_c + 273.15;
        let gamma = 1.4;  // Ratio of specific heats for air
        let r_specific = 287.05;  // Specific gas constant for air (J/kg·K)
        let speed_of_sound = (gamma * r_specific * temp_k).sqrt();
        let mach = velocity / speed_of_sound;
        
        // Base drag coefficient from drag model
        let base_cd = match self.inputs.drag_model {
            DragModel::G1 => 0.5,
            DragModel::G7 => 0.4,
            _ => 0.45,
        };
        
        // Determine projectile shape for transonic corrections
        let projectile_shape = if let Some(ref model) = self.inputs.bullet_model {
            // Try to determine shape from bullet model string
            if model.to_lowercase().contains("boat") || model.to_lowercase().contains("bt") {
                ProjectileShape::BoatTail
            } else if model.to_lowercase().contains("round") || model.to_lowercase().contains("rn") {
                ProjectileShape::RoundNose
            } else if model.to_lowercase().contains("flat") || model.to_lowercase().contains("fb") {
                ProjectileShape::FlatBase
            } else {
                // Use heuristic based on caliber, weight, and drag model
                get_projectile_shape(
                    self.inputs.diameter,
                    self.inputs.mass / 0.00006479891,  // Convert kg to grains
                    &self.inputs.drag_model.to_string()
                )
            }
        } else {
            // Use heuristic based on caliber, weight, and drag model
            get_projectile_shape(
                self.inputs.diameter,
                self.inputs.mass / 0.00006479891,  // Convert kg to grains
                &self.inputs.drag_model.to_string()
            )
        };
        
        // Apply transonic corrections
        // Enable wave drag if advanced effects are enabled
        let include_wave_drag = self.inputs.enable_advanced_effects;
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
    pub azimuth_std_dev: f64,  // Horizontal aiming variation in radians
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
            azimuth_std_dev: 0.001,  // Default horizontal spread ~0.057 degrees
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
    use rand::{thread_rng, Rng};
    use rand_distr::{Distribution, Normal};
    
    let mut rng = thread_rng();
    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut impact_positions = Vec::new();
    
    // First, calculate baseline trajectory with no variations
    let baseline_solver = TrajectorySolver::new(base_inputs.clone(), base_wind.clone(), Default::default());
    let baseline_result = baseline_solver.solve()?;
    let baseline_impact = baseline_result.points.last()
        .ok_or("No baseline trajectory points")?
        .position.clone();
    
    // Create normal distributions for variations
    let velocity_dist = Normal::new(base_inputs.muzzle_velocity, params.velocity_std_dev)
        .map_err(|e| format!("Invalid velocity distribution: {}", e))?;
    let angle_dist = Normal::new(base_inputs.launch_angle, params.angle_std_dev)
        .map_err(|e| format!("Invalid angle distribution: {}", e))?;
    let bc_dist = Normal::new(base_inputs.ballistic_coefficient, params.bc_std_dev)
        .map_err(|e| format!("Invalid BC distribution: {}", e))?;
    let wind_speed_dist = Normal::new(base_wind.speed, params.wind_speed_std_dev)
        .map_err(|e| format!("Invalid wind speed distribution: {}", e))?;
    let wind_dir_dist = Normal::new(base_wind.direction, params.wind_speed_std_dev * 0.1)  // Small variation in direction
        .map_err(|e| format!("Invalid wind direction distribution: {}", e))?;
    let azimuth_dist = Normal::new(base_inputs.azimuth_angle, params.azimuth_std_dev)
        .map_err(|e| format!("Invalid azimuth distribution: {}", e))?;
    
    // Create distribution for pointing errors (simulates shooter's aiming consistency)
    let distance = baseline_impact.z;  // Distance to target
    let pointing_error_dist = Normal::new(0.0, params.angle_std_dev * distance)
        .map_err(|e| format!("Invalid pointing distribution: {}", e))?;
    
    for _ in 0..params.num_simulations {
        // Create varied inputs
        let mut inputs = base_inputs.clone();
        inputs.muzzle_velocity = velocity_dist.sample(&mut rng).max(0.0);
        inputs.launch_angle = angle_dist.sample(&mut rng);
        inputs.ballistic_coefficient = bc_dist.sample(&mut rng).max(0.01);
        inputs.azimuth_angle = azimuth_dist.sample(&mut rng);  // Add horizontal variation
        
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
                if let Some(last_point) = result.points.last() {
                    // Calculate physical deviation from baseline impact point
                    let mut deviation = Vector3::new(
                        last_point.position.x - baseline_impact.x,  // Lateral deviation
                        last_point.position.y - baseline_impact.y,  // Vertical deviation
                        last_point.position.z - baseline_impact.z,  // Range deviation
                    );
                    
                    // Add additional pointing error to simulate realistic group sizes
                    // This represents the shooter's ability to aim consistently
                    let pointing_error_y = pointing_error_dist.sample(&mut rng);
                    deviation.y += pointing_error_y;
                    
                    impact_positions.push(deviation);
                }
            },
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
    // Binary search for the angle that hits the target
    let mut low_angle = -0.2; // radians (about -11 degrees)
    let mut high_angle = 0.2;  // radians (about 11 degrees)
    let tolerance = 0.00001;   // radians
    let max_iterations = 50;
    
    for iteration in 0..max_iterations {
        let mid_angle = (low_angle + high_angle) / 2.0;
        
        let mut test_inputs = inputs.clone();
        test_inputs.launch_angle = mid_angle;
        
        let mut solver = TrajectorySolver::new(test_inputs, wind.clone(), atmosphere.clone());
        // Make sure we calculate far enough to reach the target
        solver.set_max_range(target_distance * 2.0);
        solver.set_time_step(0.001);
        let result = solver.solve()?;
        
        eprintln!("  Iteration {}: angle = {:.6} rad ({:.4} deg)", 
                 iteration, mid_angle, mid_angle * 180.0 / std::f64::consts::PI);
        
        // Find the height at target distance
        let mut height_at_target = None;
        for i in 0..result.points.len() {
            if result.points[i].position.z >= target_distance {
                if i > 0 {
                    // Linear interpolation
                    let p1 = &result.points[i - 1];
                    let p2 = &result.points[i];
                    let t = (target_distance - p1.position.z) / (p2.position.z - p1.position.z);
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
                eprintln!("    Height at target: {:.6} m, target: {:.6} m, error: {:.6} m", 
                         height, target_height, error);
                if error.abs() < 0.001 {
                    eprintln!("  Converged!");
                    return Ok(mid_angle);
                }
                
                if error > 0.0 {
                    high_angle = mid_angle;
                } else {
                    low_angle = mid_angle;
                }
            },
            None => {
                // Trajectory didn't reach target distance, increase angle
                low_angle = mid_angle;
            }
        }
        
        if (high_angle - low_angle).abs() < tolerance {
            return Ok(mid_angle);
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
    
    // Try different BC values
    for bc in (100..1000).step_by(10) {
        let bc_value = bc as f64 / 1000.0;
        
        let inputs = BallisticInputs {
            muzzle_velocity: velocity,
            ballistic_coefficient: bc_value,
            mass,
            diameter,
            ..Default::default()
        };
        
        let solver = TrajectorySolver::new(inputs, Default::default(), Default::default());
        let result = solver.solve()?;
        
        // Calculate error
        let mut total_error = 0.0;
        for (target_dist, target_drop) in points {
            // Find drop at this distance
            let mut calculated_drop = None;
            for i in 0..result.points.len() {
                if result.points[i].position.z >= *target_dist {
                    if i > 0 {
                        // Linear interpolation
                        let p1 = &result.points[i - 1];
                        let p2 = &result.points[i];
                        let t = (target_dist - p1.position.z) / (p2.position.z - p1.position.z);
                        calculated_drop = Some(-(p1.position.y + t * (p2.position.y - p1.position.y)));
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
        }
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