// CLI API module - provides simplified interfaces for command-line tool
use crate::DragModel;
use nalgebra::Vector3;
use std::error::Error;
use std::fmt;

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
    
    pub fn solve(&self) -> Result<TrajectoryResult, BallisticsError> {
        // Simple trajectory integration using Euler method
        let mut time = 0.0;
        let mut position = Vector3::new(0.0, self.inputs.sight_height, 0.0);
        let mut velocity = Vector3::new(
            self.inputs.muzzle_velocity * self.inputs.launch_angle.cos(),
            self.inputs.muzzle_velocity * self.inputs.launch_angle.sin(),
            0.0,
        );
        
        let mut points = Vec::new();
        let mut max_height = position.y;
        
        // Calculate air density
        let air_density = calculate_air_density(&self.atmosphere);
        
        // Wind vector
        let wind_vector = Vector3::new(
            self.wind.speed * self.wind.direction.sin(),
            0.0,
            self.wind.speed * self.wind.direction.cos(),
        );
        
        // Main integration loop
        while position.x < self.max_range && position.y >= 0.0 && time < 100.0 {
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
            
            // Calculate drag
            let velocity_rel = velocity - wind_vector;
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
        
        Ok(TrajectoryResult {
            max_range: last_point.position.x,
            max_height,
            time_of_flight: last_point.time,
            impact_velocity: last_point.velocity_magnitude,
            impact_energy: last_point.kinetic_energy,
            points,
        })
    }
    
    fn calculate_drag_coefficient(&self, velocity: f64) -> f64 {
        // Simplified drag calculation
        // In reality, this would use the drag tables based on drag_model
        let mach = velocity / 343.0; // Approximate speed of sound
        
        // Basic drag curve approximation
        let base_cd = match self.inputs.drag_model {
            DragModel::G1 => 0.5,
            DragModel::G7 => 0.4,
            _ => 0.45,
        };
        
        // Transonic effects
        if mach > 0.8 && mach < 1.2 {
            base_cd * 1.5
        } else if mach > 1.2 {
            base_cd * 0.8
        } else {
            base_cd
        }
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

// Run Monte Carlo simulation
pub fn run_monte_carlo(
    base_inputs: BallisticInputs,
    params: MonteCarloParams,
) -> Result<MonteCarloResults, BallisticsError> {
    use rand::{thread_rng, Rng};
    use rand_distr::{Distribution, Normal};
    
    let mut rng = thread_rng();
    let mut ranges = Vec::new();
    let mut impact_velocities = Vec::new();
    let mut impact_positions = Vec::new();
    
    // Create normal distributions for variations
    let velocity_dist = Normal::new(base_inputs.muzzle_velocity, params.velocity_std_dev)
        .map_err(|e| format!("Invalid velocity distribution: {}", e))?;
    let angle_dist = Normal::new(base_inputs.launch_angle, params.angle_std_dev)
        .map_err(|e| format!("Invalid angle distribution: {}", e))?;
    let bc_dist = Normal::new(base_inputs.ballistic_coefficient, params.bc_std_dev)
        .map_err(|e| format!("Invalid BC distribution: {}", e))?;
    let wind_dist = Normal::new(0.0, params.wind_speed_std_dev)
        .map_err(|e| format!("Invalid wind distribution: {}", e))?;
    
    for _ in 0..params.num_simulations {
        // Create varied inputs
        let mut inputs = base_inputs.clone();
        inputs.muzzle_velocity = velocity_dist.sample(&mut rng).max(0.0);
        inputs.launch_angle = angle_dist.sample(&mut rng);
        inputs.ballistic_coefficient = bc_dist.sample(&mut rng).max(0.01);
        
        // Create varied wind
        let wind = WindConditions {
            speed: wind_dist.sample(&mut rng).abs(),
            direction: rng.gen_range(0.0..2.0 * std::f64::consts::PI),
        };
        
        // Run trajectory
        let solver = TrajectorySolver::new(inputs, wind, Default::default());
        match solver.solve() {
            Ok(result) => {
                ranges.push(result.max_range);
                impact_velocities.push(result.impact_velocity);
                if let Some(last_point) = result.points.last() {
                    impact_positions.push(last_point.position.clone());
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
                if result.points[i].position.x >= *target_dist {
                    if i > 0 {
                        // Linear interpolation
                        let p1 = &result.points[i - 1];
                        let p2 = &result.points[i];
                        let t = (target_dist - p1.position.x) / (p2.position.x - p1.position.x);
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