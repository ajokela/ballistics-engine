//! FFI bindings for iOS/Swift integration

use crate::{
    BallisticInputs, TrajectorySolver, WindConditions, AtmosphericConditions,
    DragModel, calculate_zero_angle_with_conditions,
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
}

// Helper function to convert FFI inputs to internal types
fn convert_inputs(inputs: &FFIBallisticInputs) -> BallisticInputs {
    let mut ballistic_inputs = BallisticInputs::default();
    
    ballistic_inputs.muzzle_velocity = inputs.muzzle_velocity;
    ballistic_inputs.launch_angle = inputs.launch_angle;
    ballistic_inputs.muzzle_angle = inputs.launch_angle;
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
            
            // Create result on heap
            let ffi_result = Box::new(FFITrajectoryResult {
                max_range: result.max_range,
                max_height: result.max_height,
                time_of_flight: result.time_of_flight,
                impact_velocity: result.impact_velocity,
                impact_energy: result.impact_energy,
                points,
                point_count: point_count as c_int,
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

// Get library version
#[no_mangle]
pub extern "C" fn ballistics_get_version() -> *const c_char {
    let version = CString::new("0.1.0").unwrap();
    let ptr = version.as_ptr();
    std::mem::forget(version);
    ptr
}