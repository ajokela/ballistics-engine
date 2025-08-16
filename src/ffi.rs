//! FFI bindings for iOS/Swift integration

use crate::{
    BallisticInputs, TrajectorySolver, WindConditions, AtmosphericConditions,
    DragModel, calculate_zero_angle_with_conditions,
    run_monte_carlo, MonteCarloParams,
};
use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::ptr;

// FFI-safe structures with C-compatible layouts

#[repr(C)]
pub struct FFIBallisticInputs {
    pub muzzle_velocity: c_double,        // m/s
    pub launch_angle: c_double,           // radians
    pub ballistic_coefficient: c_double,
    pub mass: c_double,                   // kg
    pub diameter: c_double,               // meters
    pub drag_model: c_int,                // 0=G1, 1=G7, 2=CustomDragTable
    pub sight_height: c_double,           // meters
    pub target_distance: c_double,        // meters
    pub temperature: c_double,            // Celsius
    pub twist_rate: c_double,             // inches per turn
    pub is_twist_right: c_int,            // 0=false, 1=true
    pub shooting_angle: c_double,         // uphill/downhill angle in radians
    pub altitude: c_double,               // meters
    pub latitude: c_double,                // degrees (use NAN if not provided)
    pub azimuth_angle: c_double,          // horizontal aiming angle in radians
    pub use_rk4: c_int,                   // 0=Euler, 1=RK4
    pub enable_wind_shear: c_int,         // 0=false, 1=true
    pub enable_trajectory_sampling: c_int, // 0=false, 1=true
    pub sample_interval: c_double,        // meters
    pub enable_pitch_damping: c_int,      // 0=false, 1=true
    pub enable_precession_nutation: c_int,// 0=false, 1=true
    pub enable_spin_drift: c_int,         // 0=false, 1=true
    pub enable_magnus: c_int,             // 0=false, 1=true
    pub enable_coriolis: c_int,           // 0=false, 1=true
}

#[repr(C)]
pub struct FFIWindConditions {
    pub speed: c_double,                  // m/s
    pub direction: c_double,              // radians (0 = North, PI/2 = East)
}

#[repr(C)]
pub struct FFIAtmosphericConditions {
    pub temperature: c_double,            // Celsius
    pub pressure: c_double,               // hPa
    pub humidity: c_double,               // percentage (0-100)
    pub altitude: c_double,               // meters
}

#[repr(C)]
pub struct FFITrajectorySample {
    pub distance: c_double,               // meters
    pub time: c_double,                   // seconds
    pub velocity_mps: c_double,           // meters per second
    pub energy_joules: c_double,          // joules
    pub drop_meters: c_double,            // meters
    pub windage_meters: c_double,         // meters
    pub mach: c_double,                   // Mach number
    pub spin_rate_rps: c_double,          // revolutions per second
}

#[repr(C)]
pub struct FFITrajectoryPoint {
    pub time: c_double,
    pub position_x: c_double,
    pub position_y: c_double,
    pub position_z: c_double,
    pub velocity_magnitude: c_double,
    pub kinetic_energy: c_double,
}

#[repr(C)]
pub struct FFITrajectoryResult {
    pub max_range: c_double,
    pub max_height: c_double,
    pub time_of_flight: c_double,
    pub impact_velocity: c_double,
    pub impact_energy: c_double,
    pub points: *mut FFITrajectoryPoint,
    pub point_count: c_int,
    pub sampled_points: *mut FFITrajectorySample,
    pub sampled_point_count: c_int,
    pub min_pitch_damping: c_double,      // NAN if not calculated
    pub transonic_mach: c_double,         // NAN if not reached
    pub final_pitch_angle: c_double,      // NAN if not calculated
    pub final_yaw_angle: c_double,        // NAN if not calculated
    pub max_yaw_angle: c_double,          // NAN if not calculated
    pub max_precession_angle: c_double,   // NAN if not calculated
}

// Monte Carlo simulation parameters
#[repr(C)]
pub struct FFIMonteCarloParams {
    pub num_simulations: c_int,
    pub velocity_std_dev: c_double,
    pub angle_std_dev: c_double,
    pub bc_std_dev: c_double,
    pub wind_speed_std_dev: c_double,
    pub target_distance: c_double,    // Use NAN if not specified
    pub base_wind_speed: c_double,    // Base wind speed in m/s
    pub base_wind_direction: c_double, // Base wind direction in radians
    pub azimuth_std_dev: c_double,    // Horizontal aiming variation in radians
}

// Monte Carlo simulation results
#[repr(C)]
pub struct FFIMonteCarloResults {
    pub ranges: *mut c_double,
    pub impact_velocities: *mut c_double,
    pub impact_positions_x: *mut c_double,
    pub impact_positions_y: *mut c_double,
    pub impact_positions_z: *mut c_double,
    pub num_results: c_int,
    pub mean_range: c_double,
    pub std_dev_range: c_double,
    pub mean_impact_velocity: c_double,
    pub std_dev_impact_velocity: c_double,
    pub hit_probability: c_double,     // If target_distance was specified
}

// Helper function to convert FFI inputs to internal types
fn convert_inputs(inputs: &FFIBallisticInputs) -> BallisticInputs {
    let mut ballistic_inputs = BallisticInputs::default();
    
    ballistic_inputs.muzzle_velocity = inputs.muzzle_velocity;
    ballistic_inputs.launch_angle = inputs.launch_angle;
    ballistic_inputs.muzzle_angle = inputs.launch_angle;
    ballistic_inputs.azimuth_angle = inputs.azimuth_angle;
    ballistic_inputs.use_rk4 = inputs.use_rk4 != 0;
    ballistic_inputs.ballistic_coefficient = inputs.ballistic_coefficient;
    ballistic_inputs.bc_value = inputs.ballistic_coefficient;
    ballistic_inputs.mass = inputs.mass;
    ballistic_inputs.bullet_mass = inputs.mass;
    ballistic_inputs.diameter = inputs.diameter;
    ballistic_inputs.bullet_diameter = inputs.diameter;
    ballistic_inputs.drag_model = match inputs.drag_model {
        1 => DragModel::G7,
        2 => DragModel::G2,
        3 => DragModel::G5,
        4 => DragModel::G6,
        5 => DragModel::G8,
        6 => DragModel::GI,
        7 => DragModel::GS,
        _ => DragModel::G1,
    };
    ballistic_inputs.bc_type = ballistic_inputs.drag_model.clone();
    ballistic_inputs.sight_height = inputs.sight_height;
    ballistic_inputs.target_distance = inputs.target_distance;
    ballistic_inputs.temperature = inputs.temperature;
    ballistic_inputs.twist_rate = inputs.twist_rate;
    ballistic_inputs.is_twist_right = inputs.is_twist_right != 0;
    ballistic_inputs.shooting_angle = inputs.shooting_angle;
    ballistic_inputs.altitude = inputs.altitude;
    
    if !inputs.latitude.is_nan() {
        ballistic_inputs.latitude = Some(inputs.latitude);
    }
    
    // Set derived values
    ballistic_inputs.caliber_inches = inputs.diameter / 0.0254;
    ballistic_inputs.weight_grains = inputs.mass / 0.00006479891;
    ballistic_inputs.bullet_length = inputs.diameter * 4.0;
    
    // New advanced physics flags
    ballistic_inputs.enable_wind_shear = inputs.enable_wind_shear != 0;
    ballistic_inputs.enable_trajectory_sampling = inputs.enable_trajectory_sampling != 0;
    ballistic_inputs.sample_interval = inputs.sample_interval;
    ballistic_inputs.enable_pitch_damping = inputs.enable_pitch_damping != 0;
    ballistic_inputs.enable_precession_nutation = inputs.enable_precession_nutation != 0;
    ballistic_inputs.use_enhanced_spin_drift = inputs.enable_spin_drift != 0;
    ballistic_inputs.enable_advanced_effects = inputs.enable_magnus != 0 || inputs.enable_coriolis != 0;
    
    ballistic_inputs
}

// Main trajectory calculation function for FFI
#[no_mangle]
pub extern "C" fn ballistics_calculate_trajectory(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    max_range: c_double,
    step_size: c_double,
) -> *mut FFITrajectoryResult {
    if inputs.is_null() {
        return ptr::null_mut();
    }
    
    let inputs = unsafe { &*inputs };
    let ballistic_inputs = convert_inputs(inputs);
    
    let wind_conditions = if wind.is_null() {
        WindConditions::default()
    } else {
        let wind = unsafe { &*wind };
        WindConditions {
            speed: wind.speed,
            direction: wind.direction,
        }
    };
    
    let atmospheric_conditions = if atmosphere.is_null() {
        AtmosphericConditions::default()
    } else {
        let atmo = unsafe { &*atmosphere };
        AtmosphericConditions {
            temperature: atmo.temperature,
            pressure: atmo.pressure,
            humidity: atmo.humidity,
            altitude: atmo.altitude,
        }
    };
    
    // Create solver and calculate trajectory
    let mut solver = TrajectorySolver::new(
        ballistic_inputs,
        wind_conditions,
        atmospheric_conditions,
    );
    
    // Set max range and time step
    solver.set_max_range(max_range);
    solver.set_time_step(step_size / 1000.0); // Convert to time step
    
    match solver.solve() {
        Ok(result) => {
            // Convert trajectory points to FFI format
            let point_count = result.points.len();
            let points = if point_count > 0 {
                let mut ffi_points = Vec::with_capacity(point_count);
                for point in &result.points {
                    ffi_points.push(FFITrajectoryPoint {
                        time: point.time,
                        position_x: point.position[0],
                        position_y: point.position[1],
                        position_z: point.position[2],
                        velocity_magnitude: point.velocity_magnitude,
                        kinetic_energy: point.kinetic_energy,
                    });
                }
                let points_ptr = ffi_points.as_mut_ptr();
                std::mem::forget(ffi_points); // Prevent deallocation
                points_ptr
            } else {
                ptr::null_mut()
            };
            
            // Convert sampled points if available
            let (sampled_points, sampled_point_count) = if let Some(ref samples) = result.sampled_points {
                let mut ffi_samples = Vec::with_capacity(samples.len());
                for sample in samples {
                    ffi_samples.push(FFITrajectorySample {
                        distance: sample.distance_m,
                        time: sample.time_s,
                        velocity_mps: sample.velocity_mps,
                        energy_joules: sample.energy_j,
                        drop_meters: sample.drop_m,
                        windage_meters: sample.wind_drift_m,
                        mach: 0.0,  // Would need to calculate from velocity
                        spin_rate_rps: 0.0,  // Not available in TrajectorySample
                    });
                }
                let count = ffi_samples.len() as c_int;
                let samples_ptr = ffi_samples.as_mut_ptr();
                std::mem::forget(ffi_samples);
                (samples_ptr, count)
            } else {
                (ptr::null_mut(), 0)
            };
            
            // Extract angular state values if available
            let (final_pitch, final_yaw, max_yaw, max_prec) = if let Some(ref angular) = result.angular_state {
                (
                    angular.pitch_angle,
                    angular.yaw_angle,
                    result.max_yaw_angle.unwrap_or(std::f64::NAN),
                    result.max_precession_angle.unwrap_or(std::f64::NAN),
                )
            } else {
                (std::f64::NAN, std::f64::NAN, std::f64::NAN, std::f64::NAN)
            };
            
            // Create result on heap
            let ffi_result = Box::new(FFITrajectoryResult {
                max_range: result.max_range,
                max_height: result.max_height,
                time_of_flight: result.time_of_flight,
                impact_velocity: result.impact_velocity,
                impact_energy: result.impact_energy,
                points,
                point_count: point_count as c_int,
                sampled_points,
                sampled_point_count,
                min_pitch_damping: result.min_pitch_damping.unwrap_or(std::f64::NAN),
                transonic_mach: result.transonic_mach.unwrap_or(std::f64::NAN),
                final_pitch_angle: final_pitch,
                final_yaw_angle: final_yaw,
                max_yaw_angle: max_yaw,
                max_precession_angle: max_prec,
            });
            
            Box::into_raw(ffi_result)
        }
        Err(_) => ptr::null_mut(),
    }
}

// Free allocated trajectory result
#[no_mangle]
pub extern "C" fn ballistics_free_trajectory_result(result: *mut FFITrajectoryResult) {
    if !result.is_null() {
        unsafe {
            let result = Box::from_raw(result);
            if !result.points.is_null() && result.point_count > 0 {
                let points = Vec::from_raw_parts(
                    result.points,
                    result.point_count as usize,
                    result.point_count as usize,
                );
                drop(points);
            }
            if !result.sampled_points.is_null() && result.sampled_point_count > 0 {
                let samples = Vec::from_raw_parts(
                    result.sampled_points,
                    result.sampled_point_count as usize,
                    result.sampled_point_count as usize,
                );
                drop(samples);
            }
            drop(result);
        }
    }
}

// Calculate zero angle for a given target distance
#[no_mangle]
pub extern "C" fn ballistics_calculate_zero_angle(
    inputs: *const FFIBallisticInputs,
    wind: *const FFIWindConditions,
    atmosphere: *const FFIAtmosphericConditions,
    zero_distance: c_double,
) -> c_double {
    if inputs.is_null() {
        return f64::NAN;
    }
    
    let inputs = unsafe { &*inputs };
    let ballistic_inputs = convert_inputs(inputs);
    
    let wind_conditions = if wind.is_null() {
        WindConditions::default()
    } else {
        let wind = unsafe { &*wind };
        WindConditions {
            speed: wind.speed,
            direction: wind.direction,
        }
    };
    
    let atmospheric_conditions = if atmosphere.is_null() {
        AtmosphericConditions::default()
    } else {
        let atmo = unsafe { &*atmosphere };
        AtmosphericConditions {
            temperature: atmo.temperature,
            pressure: atmo.pressure,
            humidity: atmo.humidity,
            altitude: atmo.altitude,
        }
    };
    
    // For zero angle, we want the bullet to hit at sight height at the zero distance
    // This means the bullet crosses the line of sight at the zero distance
    let target_height = ballistic_inputs.sight_height;
    
    eprintln!("FFI: Calculating zero angle for:");
    eprintln!("  Zero distance: {} m", zero_distance);
    eprintln!("  Target height: {} m", target_height);
    eprintln!("  Sight height: {} m", ballistic_inputs.sight_height);
    eprintln!("  Wind speed: {} m/s", wind_conditions.speed);
    eprintln!("  Temperature: {} C", atmospheric_conditions.temperature);
    
    match calculate_zero_angle_with_conditions(
        ballistic_inputs,
        zero_distance,
        target_height,
        wind_conditions,
        atmospheric_conditions,
    ) {
        Ok(angle) => {
            eprintln!("  Calculated angle: {} rad ({} deg)", angle, angle * 180.0 / std::f64::consts::PI);
            angle
        },
        Err(e) => {
            eprintln!("  Error: {:?}", e);
            f64::NAN
        }
    }
}

// Simple trajectory calculation for quick results
#[no_mangle]
pub extern "C" fn ballistics_quick_trajectory(
    muzzle_velocity: c_double,
    bc: c_double,
    sight_height: c_double,
    zero_distance: c_double,
    target_distance: c_double,
) -> c_double {
    // This provides a simple drop calculation at target distance
    // Using simplified ballistic calculations
    
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = muzzle_velocity;
    inputs.ballistic_coefficient = bc;
    inputs.bc_value = bc;
    inputs.sight_height = sight_height;
    inputs.target_distance = target_distance;
    
    let wind = WindConditions::default();
    let atmo = AtmosphericConditions::default();
    
    // First calculate zero angle
    let zero_angle = match calculate_zero_angle_with_conditions(inputs.clone(), zero_distance, sight_height, wind.clone(), atmo.clone()) {
        Ok(angle) => angle,
        Err(_) => return f64::NAN,
    };
    
    // Now calculate trajectory with that zero angle
    inputs.launch_angle = zero_angle;
    inputs.muzzle_angle = zero_angle;
    
    let mut solver = TrajectorySolver::new(inputs, wind, atmo);
    solver.set_max_range(target_distance * 1.1);
    
    match solver.solve() {
        Ok(result) => {
            // Find the drop at target distance
            for point in result.points {
                if point.position[0] >= target_distance {
                    return -point.position[1]; // Return drop (negative y is drop)
                }
            }
            f64::NAN
        }
        Err(_) => f64::NAN,
    }
}

// Monte Carlo simulation
#[no_mangle]
pub extern "C" fn ballistics_monte_carlo(
    inputs: *const FFIBallisticInputs,
    atmosphere: *const FFIAtmosphericConditions,
    params: *const FFIMonteCarloParams,
) -> *mut FFIMonteCarloResults {
    if inputs.is_null() || params.is_null() {
        return ptr::null_mut();
    }
    
    let inputs = unsafe { &*inputs };
    let params = unsafe { &*params };
    
    // Convert FFI inputs to internal types
    let ballistic_inputs = convert_inputs(inputs);
    
    // Note: Atmospheric conditions are already included in the conversion
    // from FFIBallisticInputs (temperature, altitude, etc.)
    
    // Create Monte Carlo parameters
    let mc_params = MonteCarloParams {
        num_simulations: params.num_simulations as usize,
        velocity_std_dev: params.velocity_std_dev,
        angle_std_dev: params.angle_std_dev,
        bc_std_dev: params.bc_std_dev,
        wind_speed_std_dev: params.wind_speed_std_dev,
        target_distance: if params.target_distance.is_nan() {
            None
        } else {
            Some(params.target_distance)
        },
        base_wind_speed: params.base_wind_speed,
        base_wind_direction: params.base_wind_direction,
        azimuth_std_dev: params.azimuth_std_dev,
    };
    
    // Run Monte Carlo simulation
    match run_monte_carlo(ballistic_inputs, mc_params) {
        Ok(results) => {
            let num_results = results.ranges.len() as c_int;
            
            // Calculate statistics
            let mean_range: f64 = results.ranges.iter().sum::<f64>() / num_results as f64;
            let variance_range: f64 = results.ranges.iter()
                .map(|r| (r - mean_range).powi(2))
                .sum::<f64>() / num_results as f64;
            let std_dev_range = variance_range.sqrt();
            
            let mean_velocity: f64 = results.impact_velocities.iter().sum::<f64>() / num_results as f64;
            let variance_velocity: f64 = results.impact_velocities.iter()
                .map(|v| (v - mean_velocity).powi(2))
                .sum::<f64>() / num_results as f64;
            let std_dev_velocity = variance_velocity.sqrt();
            
            // Calculate hit probability if target distance was specified
            let hit_probability = if params.target_distance.is_nan() {
                0.0
            } else {
                let target = params.target_distance;
                let hit_radius = 0.3; // 30cm radius for hit zone
                let hits = results.impact_positions.iter()
                    .filter(|pos| {
                        let distance = (pos.x.powi(2) + pos.y.powi(2)).sqrt();
                        distance < target && pos.norm() < hit_radius
                    })
                    .count();
                hits as f64 / num_results as f64
            };
            
            // Allocate memory for arrays
            let ranges_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap()
                ) as *mut c_double;
                for (i, &range) in results.ranges.iter().enumerate() {
                    *ptr.add(i) = range;
                }
                ptr
            };
            
            let velocities_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap()
                ) as *mut c_double;
                for (i, &vel) in results.impact_velocities.iter().enumerate() {
                    *ptr.add(i) = vel;
                }
                ptr
            };
            
            let pos_x_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap()
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.x;
                }
                ptr
            };
            
            let pos_y_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap()
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.y;
                }
                ptr
            };
            
            let pos_z_ptr = unsafe {
                let ptr = std::alloc::alloc(
                    std::alloc::Layout::array::<c_double>(num_results as usize).unwrap()
                ) as *mut c_double;
                for (i, pos) in results.impact_positions.iter().enumerate() {
                    *ptr.add(i) = pos.z;
                }
                ptr
            };
            
            // Create result structure
            let result = Box::new(FFIMonteCarloResults {
                ranges: ranges_ptr,
                impact_velocities: velocities_ptr,
                impact_positions_x: pos_x_ptr,
                impact_positions_y: pos_y_ptr,
                impact_positions_z: pos_z_ptr,
                num_results,
                mean_range,
                std_dev_range,
                mean_impact_velocity: mean_velocity,
                std_dev_impact_velocity: std_dev_velocity,
                hit_probability,
            });
            
            Box::into_raw(result)
        }
        Err(_) => ptr::null_mut()
    }
}

// Free Monte Carlo results
#[no_mangle]
pub extern "C" fn ballistics_free_monte_carlo_results(results: *mut FFIMonteCarloResults) {
    if results.is_null() {
        return;
    }
    
    unsafe {
        let results = Box::from_raw(results);
        let num = results.num_results as usize;
        
        // Free arrays
        if !results.ranges.is_null() {
            std::alloc::dealloc(
                results.ranges as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap()
            );
        }
        
        if !results.impact_velocities.is_null() {
            std::alloc::dealloc(
                results.impact_velocities as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap()
            );
        }
        
        if !results.impact_positions_x.is_null() {
            std::alloc::dealloc(
                results.impact_positions_x as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap()
            );
        }
        
        if !results.impact_positions_y.is_null() {
            std::alloc::dealloc(
                results.impact_positions_y as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap()
            );
        }
        
        if !results.impact_positions_z.is_null() {
            std::alloc::dealloc(
                results.impact_positions_z as *mut u8,
                std::alloc::Layout::array::<c_double>(num).unwrap()
            );
        }
        
        // Box automatically deallocates the result structure
    }
}

// Get library version
#[no_mangle]
pub extern "C" fn ballistics_get_version() -> *const c_char {
    let version = CString::new("0.3.0").unwrap();
    let ptr = version.as_ptr();
    std::mem::forget(version);
    ptr
}