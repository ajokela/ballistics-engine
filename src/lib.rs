//! # Ballistics Rust Acceleration Module
//! 
//! This module provides Rust-accelerated implementations of ballistics calculations
//! for improved performance while maintaining full compatibility with the Python implementation.
//! 
//! ## Testing Strategy
//! 
//! Due to PyO3's requirements for Python runtime linking, traditional `cargo test` 
//! cannot be used. Instead, all functionality is tested through:
//! 
//! - **Python Integration Tests**: `tests/test_rust_integration.py` (12 tests)
//! - **Python Drag Tests**: `tests/test_rust_drag.py` (11 tests) 
//! - **Performance Validation**: Benchmarking and accuracy comparisons
//! - **Compilation Validation**: `cargo check` ensures code correctness
//! 
//! This approach provides comprehensive validation while working within PyO3's constraints.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule};
use numpy::{PyArray1, PyReadonlyArray1, ToPyArray, PyArrayMethods};
use nalgebra::Vector3;

mod derivatives;
mod atmosphere;
mod drag;
mod stability;
mod trajectory_solver;
mod trajectory_sampling;
mod angle_calculations;
mod constants;
mod transonic_drag;
mod wind;
mod spin_drift;
mod bc_estimation;
mod spin_decay;
mod fast_trajectory;
mod reynolds;
mod wind_shear;
mod pitch_damping;
mod precession_nutation;
mod aerodynamic_jump;
mod monte_carlo;
mod cluster_bc;

#[cfg(test)]
mod drag_pure_test;
mod form_factor;

use derivatives::compute_derivatives;
use drag::{get_drag_coefficient, interpolated_bc_py};
use atmosphere::{calculate_atmosphere as calc_atmosphere, calculate_air_density_cipm, get_local_atmosphere};
use transonic_drag::{transonic_correction_py, get_projectile_shape_py};
use reynolds::apply_reynolds_correction_py;
use stability::{compute_stability_coefficient, compute_spin_drift};
use trajectory_solver::{prepare_initial_conditions, post_process_trajectory};
use trajectory_sampling::{sample_trajectory, trajectory_samples_to_dicts, TrajectoryData, TrajectoryOutputs};
use angle_calculations::{brent_root_find, zero_angle, solve_muzzle_angle, adjusted_muzzle_velocity};
use wind::{WindSock, WindSegment};
use wind_shear::{WindShearModel, WindShearProfile, WindLayer};
use aerodynamic_jump::calculate_aerodynamic_jump;
use spin_decay::SpinDecayParameters;
use bc_estimation::BCSegmentEstimator;
use monte_carlo::{monte_carlo_parallel, monte_carlo_batch, monte_carlo_streaming};

/// Rust implementation of the ballistics derivatives function
#[pyfunction]
#[pyo3(signature = (_t, state, inputs, wind_segments, atmos_params, bc_used, _target_horizontal_dist_m, _target_vertical_height_m, omega_vector=None))]
fn derivatives_rust(
    py: Python,
    _t: f64,
    state: PyReadonlyArray1<f64>,
    inputs: &Bound<'_, PyDict>,
    wind_segments: Vec<WindSegment>,
    atmos_params: PyReadonlyArray1<f64>,
    bc_used: f64,
    _target_horizontal_dist_m: f64,
    _target_vertical_height_m: f64,
    omega_vector: Option<&Bound<'_, PyAny>>,
) -> PyResult<Py<PyArray1<f64>>> {
    // Convert inputs to Rust types
    let state_array = state.as_array();
    let atmos_array = atmos_params.as_array();
    
    if state_array.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("State array must have 6 elements"));
    }
    
    if atmos_array.len() != 4 {
        return Err(pyo3::exceptions::PyValueError::new_err("Atmospheric parameters must have 4 elements"));
    }
    
    // Extract state vector
    let pos = Vector3::new(state_array[0], state_array[1], state_array[2]);
    let vel = Vector3::new(state_array[3], state_array[4], state_array[5]);
    
    // Extract atmospheric parameters
    let atmos_params_tuple = (atmos_array[0], atmos_array[1], atmos_array[2], atmos_array[3]);
    
    // Extract ballistic inputs
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    
    // Create appropriate wind sock based on wind shear settings
    let wind_vector = if ballistic_inputs.enable_wind_shear && ballistic_inputs.wind_shear_model != "none" {
        // Create WindShearWindSock for altitude-dependent wind
        let mut profile = WindShearProfile::default();
        profile.model = match ballistic_inputs.wind_shear_model.as_str() {
            "logarithmic" => WindShearModel::Logarithmic,
            "power_law" => WindShearModel::PowerLaw,
            "ekman_spiral" => WindShearModel::EkmanSpiral,
            _ => WindShearModel::None,
        };
        
        // Set surface wind from first segment
        // Note: The wind speed should be the measurement at reference_height (typically 10m)
        // The altitude_m field is not used in calculations, only the wind vector
        if !wind_segments.is_empty() {
            profile.surface_wind = WindLayer {
                altitude_m: 0.0,  // Not used in wind shear calculations
                speed_mps: wind_segments[0].0 * 0.2777778, // km/h to m/s at reference_height
                direction_deg: wind_segments[0].1,
            };
        }
        
        let wind_shear_sock = wind_shear::WindShearWindSock::with_shooter_altitude(
            wind_segments, 
            Some(profile), 
            ballistic_inputs.altitude
        );
        wind_shear_sock.vector_for_position(pos)
    } else {
        // Use regular WindSock for constant wind
        let wind_sock = WindSock::new(wind_segments);
        wind_sock.vector_for_range_stateless(pos[0])
    };
    
    // Convert omega vector if provided
    let omega_vec = if let Some(omega) = omega_vector {
        // Check if the Python object is None
        if omega.is_none() {
            None
        } else {
            // Try to extract as array
            match omega.extract::<PyReadonlyArray1<f64>>() {
                Ok(omega_array) => {
                    let omega_slice = omega_array.as_array();
                    if omega_slice.len() != 3 {
                        return Err(pyo3::exceptions::PyValueError::new_err("Omega vector must have 3 elements"));
                    }
                    Some(Vector3::new(omega_slice[0], omega_slice[1], omega_slice[2]))
                },
                Err(_) => None  // If extraction fails, treat as None
            }
        }
    } else {
        None
    };
    
    // Calculate derivatives
    let result = compute_derivatives(
        pos,
        vel,
        &ballistic_inputs,
        wind_vector,
        atmos_params_tuple,
        bc_used,
        omega_vec,
        _t,  // Pass time for enhanced spin drift
    )?;
    
    // Convert result back to Python array
    let result_array = [result[0], result[1], result[2], result[3], result[4], result[5]];
    Ok(result_array.to_pyarray_bound(py).into())
}

/// Extract ballistic inputs from Python dictionary
pub fn extract_ballistic_inputs(inputs: &Bound<'_, PyDict>) -> PyResult<BallisticInputs> {
    let bc_value: f64 = inputs.get_item("bc_value")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing required field: bc_value"))?
        .extract()?;
    let bc_type: String = inputs.get_item("bc_type")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing required field: bc_type"))?
        .extract()?;
    let bullet_mass: f64 = inputs.get_item("bullet_mass")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing required field: bullet_mass"))?
        .extract()?;
    let muzzle_velocity: f64 = inputs.get_item("muzzle_velocity")?.map_or(Ok(0.0), |v| v.extract())?;
    let altitude: f64 = inputs.get_item("altitude")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing required field: altitude"))?
        .extract()?;
    let twist_rate: f64 = inputs.get_item("twist_rate")?.map_or(Ok(0.0), |v| v.extract())?;
    let bullet_length: f64 = inputs.get_item("bullet_length")?.map_or(Ok(0.0), |v| v.extract())?;
    let bullet_diameter: f64 = inputs.get_item("bullet_diameter")?.map_or(Ok(0.0), |v| v.extract())?;
    let target_distance: f64 = inputs.get_item("target_distance")?.map_or(Ok(0.0), |v| v.extract())?;
    let muzzle_angle: f64 = inputs.get_item("muzzle_angle")?.map_or(Ok(0.0), |v| v.extract())?;
    let wind_speed: f64 = inputs.get_item("wind_speed")?.map_or(Ok(0.0), |v| v.extract())?;
    let wind_angle: f64 = inputs.get_item("wind_angle")?.map_or(Ok(0.0), |v| v.extract())?;
    let temperature: f64 = inputs.get_item("temperature")?.map_or(Ok(15.0), |v| v.extract())?;
    let pressure: f64 = inputs.get_item("pressure")?.map_or(Ok(1013.25), |v| v.extract())?;
    let humidity: f64 = inputs.get_item("humidity")?.map_or(Ok(0.0), |v| v.extract())?;
    let latitude: Option<f64> = inputs.get_item("latitude")?.map(|v| v.extract()).transpose()?;
    let enable_advanced_effects: bool = inputs.get_item("enable_advanced_effects")?.map_or(Ok(false), |v| v.extract())?;
    let is_twist_right: bool = inputs.get_item("is_twist_right")?.map_or(Ok(true), |v| v.extract())?;
    let shooting_angle: f64 = inputs.get_item("shooting_angle")?.map_or(Ok(0.0), |v| v.extract())?;
    let use_powder_sensitivity: bool = inputs.get_item("use_powder_sensitivity")?.map_or(Ok(false), |v| v.extract())?;
    let powder_temp_sensitivity: f64 = inputs.get_item("powder_temp_sensitivity")?.map_or(Ok(0.0), |v| v.extract())?;
    let powder_temp: f64 = inputs.get_item("powder_temp")?.map_or(Ok(70.0), |v| v.extract())?;
    let tipoff_yaw: f64 = inputs.get_item("tipoff_yaw")?.map_or(Ok(0.0), |v| v.extract())?;
    let tipoff_decay_distance: f64 = inputs.get_item("tipoff_decay_distance")?.map_or(Ok(20.0), |v| v.extract())?;
    let ground_threshold: f64 = inputs.get_item("ground_threshold")?.map_or(Ok(0.0), |v| v.extract())?;
    
    // Extract wind shear feature flags
    let enable_wind_shear: bool = inputs.get_item("enable_wind_shear")?.map_or(Ok(false), |v| v.extract())?;
    let wind_shear_model: String = inputs.get_item("wind_shear_model")?.map_or(Ok("none".to_string()), |v| v.extract())?;
    
    // Extract cluster BC fields
    let use_cluster_bc: bool = inputs.get_item("use_cluster_bc")?.map_or(Ok(false), |v| v.extract())?;
    let bullet_cluster: Option<usize> = inputs.get_item("bullet_cluster")?.map(|v| v.extract()).transpose()?;
    
    // Parse bc_segments if present
    let bc_segments = if let Some(segments) = inputs.get_item("bc_segments")? {
        if segments.is_none() {
            None
        } else {
            Some(extract_bc_segments(&segments)?)
        }
    } else {
        None
    };
    
    let bc_type_enum = match bc_type.as_str() {
        "G1" | "DragModel.G1" | "DragModel::G1" => DragModel::G1,
        "G7" | "DragModel.G7" | "DragModel::G7" => DragModel::G7,
        s if s.contains("G1") => DragModel::G1,
        s if s.contains("G7") => DragModel::G7,
        _ => return Err(pyo3::exceptions::PyValueError::new_err(format!("Invalid BC type: '{bc_type}'. Expected G1 or G7."))),
    };
    
    // Get caliber and weight for drag calculations
    let caliber_inches = inputs.get_item("caliber_inches")?
        .and_then(|v| v.extract::<f64>().ok())
        .unwrap_or(0.0);
    let weight_grains = inputs.get_item("weight_grains")?
        .and_then(|v| v.extract::<f64>().ok())
        .unwrap_or(0.0);
    
    // Get velocity-based BC segments fields
    let use_bc_segments = inputs.get_item("use_bc_segments")?.map_or(Ok(false), |v| v.extract())?;
    
    let bullet_id = inputs.get_item("bullet_id")?.and_then(|v| {
        if v.is_none() {
            None
        } else {
            Some(v.extract::<String>())
        }
    }).transpose()?;
    let bc_segments_data = if let Some(data) = inputs.get_item("bc_segments_data")? {
        if data.is_none() {
            None
        } else {
            Some(extract_bc_segments_data(&data)?)
        }
    } else {
        None
    };
    
    // Get form factor and BC optimization fields
    let use_form_factor = inputs.get_item("use_form_factor")?.map_or(Ok(true), |v| v.extract())?;
    let manufacturer = inputs.get_item("manufacturer")?.and_then(|v| {
        if v.is_none() {
            None
        } else {
            Some(v.extract::<String>())
        }
    }).transpose()?;
    let bullet_model = inputs.get_item("bullet_model")?.and_then(|v| {
        if v.is_none() {
            None
        } else {
            Some(v.extract::<String>())
        }
    }).transpose()?;
    
    Ok(BallisticInputs {
        bc_value,
        bc_type: bc_type_enum,
        bullet_mass,
        muzzle_velocity,
        altitude,
        twist_rate,
        bullet_length,
        bullet_diameter,
        target_distance,
        muzzle_angle,
        wind_speed,
        wind_angle,
        temperature,
        pressure,
        humidity,
        latitude,
        enable_advanced_effects,
        is_twist_right,
        shooting_angle,
        use_powder_sensitivity,
        powder_temp_sensitivity,
        powder_temp,
        tipoff_yaw,
        tipoff_decay_distance,
        ground_threshold,
        bc_segments,
        caliber_inches,
        weight_grains,
        use_bc_segments,
        bullet_id,
        bc_segments_data,
        use_enhanced_spin_drift: false, // Disable enhanced spin drift - use standard Magnus instead
        use_form_factor,
        manufacturer,
        bullet_model,
        enable_wind_shear,
        wind_shear_model,
        use_cluster_bc,
        bullet_cluster,
    })
}

/// Extract BC segments from Python object
fn extract_bc_segments(segments: &Bound<'_, PyAny>) -> PyResult<Vec<(f64, f64)>> {
    if segments.is_none() {
        return Ok(Vec::new());
    }
    
    // Try to extract as Vec<(f64, f64)> first (most common case)
    if let Ok(segments_tuple_list) = segments.extract::<Vec<(f64, f64)>>() {
        Ok(segments_tuple_list)
    } else {
        // Fallback: try to extract as a list of objects and convert each
        Err(pyo3::exceptions::PyValueError::new_err("Could not extract BC segments - expected Vec<(f64, f64)>"))
    }
}

/// Extract velocity-based BC segments data from Python object
fn extract_bc_segments_data(data: &Bound<'_, PyAny>) -> PyResult<Vec<BCSegmentData>> {
    if data.is_none() {
        return Ok(Vec::new());
    }
    
    // Try to iterate over the list
    let mut result = Vec::new();
    
    if let Ok(list) = data.downcast::<pyo3::types::PyList>() {
        for item in list.iter() {
            if let Ok(dict) = item.downcast::<PyDict>() {
                let velocity_min = dict.get_item("velocity_min")?
                    .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("BC segment missing velocity_min"))?
                    .extract::<f64>()?;
                let velocity_max = dict.get_item("velocity_max")?
                    .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("BC segment missing velocity_max"))?
                    .extract::<f64>()?;
                let bc_value = dict.get_item("bc_value")?
                    .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("BC segment missing bc_value"))?
                    .extract::<f64>()?;
                
                result.push(BCSegmentData {
                    velocity_min,
                    velocity_max,
                    bc_value,
                });
            }
        }
    }
    
    Ok(result)
}

/// Get wind vector from Python wind sock object
fn get_wind_vector(wind_sock: &Bound<'_, PyAny>, range_m: f64) -> PyResult<Vector3<f64>> {
    let wind_vector_method = wind_sock.getattr("vector_for_range")?;
    let wind_result = wind_vector_method.call1((range_m,))?;
    
    // Try multiple extraction methods for maximum compatibility
    let wind_vector = if let Ok(np_array) = wind_result.downcast::<PyArray1<f64>>() {
        // Handle numpy array
        let array = np_array.readonly();
        let array_slice = array.as_array();
        if array_slice.len() != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!("Wind vector must have 3 elements, got {}", array_slice.len())));
        }
        Vector3::new(array_slice[0], array_slice[1], array_slice[2])
    } else if let Ok(list_result) = wind_result.extract::<Vec<f64>>() {
        // Handle list of floats
        if list_result.len() != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!("Wind vector must have 3 elements, got {}", list_result.len())));
        }
        Vector3::new(list_result[0], list_result[1], list_result[2])
    } else if let Ok(tuple_result) = wind_result.extract::<(f64, f64, f64)>() {
        // Handle tuple
        Vector3::new(tuple_result.0, tuple_result.1, tuple_result.2)
    } else {
        // Try to extract as array-like sequence
        return Err(pyo3::exceptions::PyValueError::new_err(format!("Could not convert wind vector result to 3-element array. Got type: {}", wind_result.get_type().name()?)))
    };
    
    Ok(wind_vector)
}

/// Ballistic inputs structure
#[derive(Debug, Clone)]
pub struct BallisticInputs {
    bc_value: f64,
    bc_type: DragModel,
    bullet_mass: f64,
    muzzle_velocity: f64,
    altitude: f64,
    twist_rate: f64,
    bullet_length: f64,
    bullet_diameter: f64,
    target_distance: f64,
    muzzle_angle: f64,
    wind_speed: f64,
    wind_angle: f64,
    temperature: f64,
    pressure: f64,
    humidity: f64,
    latitude: Option<f64>,
    enable_advanced_effects: bool,
    is_twist_right: bool,
    shooting_angle: f64,
    use_powder_sensitivity: bool,
    powder_temp_sensitivity: f64,
    powder_temp: f64,
    tipoff_yaw: f64,
    tipoff_decay_distance: f64,
    ground_threshold: f64,
    bc_segments: Option<Vec<(f64, f64)>>,
    // Drag calculation fields
    caliber_inches: f64,
    weight_grains: f64,
    // Velocity-based BC segments
    use_bc_segments: bool,
    bullet_id: Option<String>,
    bc_segments_data: Option<Vec<BCSegmentData>>,
    use_enhanced_spin_drift: bool,
    // Form factor and BC optimization fields
    use_form_factor: bool,
    manufacturer: Option<String>,
    bullet_model: Option<String>,
    // Wind shear fields
    enable_wind_shear: bool,
    wind_shear_model: String,
    // Cluster BC fields
    use_cluster_bc: bool,
    bullet_cluster: Option<usize>,
}

/// Velocity-based BC segment data
#[derive(Debug, Clone)]
pub struct BCSegmentData {
    pub velocity_min: f64,
    pub velocity_max: f64,
    pub bc_value: f64,
}

impl BallisticInputs {
    /// Get bc_type as string for ML model
    pub fn bc_type_str(&self) -> &str {
        match self.bc_type {
            DragModel::G1 => "G1",
            DragModel::G7 => "G7",
        }
    }
}

/// Drag model enumeration
#[derive(Debug, Clone)]
pub enum DragModel {
    G1,
    G7,
}

impl std::fmt::Display for DragModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DragModel::G1 => write!(f, "G1"),
            DragModel::G7 => write!(f, "G7"),
        }
    }
}

/// Python wrapper for drag coefficient calculation
#[pyfunction]
fn get_drag_coefficient_rust(mach: f64, drag_model: &str) -> PyResult<f64> {
    let model = match drag_model.to_uppercase().as_str() {
        "G1" => DragModel::G1,
        "G7" => DragModel::G7,
        _ => return Err(pyo3::exceptions::PyValueError::new_err(format!("Invalid drag model: '{drag_model}'. Expected G1 or G7."))),
    };
    
    Ok(get_drag_coefficient(mach, &model))
}

/// Python wrapper for drag coefficient calculation with transonic corrections
#[pyfunction]
#[pyo3(signature = (mach, drag_model, apply_transonic_correction=true, caliber=None, weight_grains=None))]
fn get_drag_coefficient_rust_transonic(
    mach: f64, 
    drag_model: &str,
    apply_transonic_correction: bool,
    caliber: Option<f64>,
    weight_grains: Option<f64>,
) -> PyResult<f64> {
    let model = match drag_model.to_uppercase().as_str() {
        "G1" => DragModel::G1,
        "G7" => DragModel::G7,
        _ => return Err(pyo3::exceptions::PyValueError::new_err(format!("Invalid drag model: '{drag_model}'. Expected G1 or G7."))),
    };
    
    Ok(drag::get_drag_coefficient_with_transonic(
        mach, 
        &model,
        apply_transonic_correction,
        None, // Let Rust determine shape
        caliber,
        weight_grains,
    ))
}

/// Python wrapper for atmospheric calculation
#[pyfunction]
#[pyo3(signature = (altitude_m, temp_override_c=None, press_override_hpa=None, humidity_percent=0.0))]
fn calculate_atmosphere_rust(
    _py: Python,
    altitude_m: f64,
    temp_override_c: Option<f64>,
    press_override_hpa: Option<f64>,
    humidity_percent: f64,
) -> PyResult<(f64, f64)> {
    Ok(calc_atmosphere(altitude_m, temp_override_c, press_override_hpa, humidity_percent))
}

/// Python wrapper for CIPM air density calculation
#[pyfunction]
fn calculate_air_density_cipm_rust(
    temp_c: f64,
    pressure_hpa: f64,
    humidity_percent: f64,
) -> PyResult<f64> {
    Ok(calculate_air_density_cipm(temp_c, pressure_hpa, humidity_percent))
}

/// Python wrapper for local atmosphere calculation
#[pyfunction]
fn get_local_atmosphere_rust(
    _py: Python,
    altitude_m: f64,
    atmos_params: PyReadonlyArray1<f64>,
) -> PyResult<(f64, f64)> {
    let params = atmos_params.as_array();
    
    // Handle the two different parameter formats
    if params.len() == 2 {
        // Direct atmosphere values
        Ok((params[0], params[1]))
    } else if params.len() == 4 {
        // Calculate from base parameters
        Ok(get_local_atmosphere(
            altitude_m,
            params[0],  // base_alt
            params[1],  // base_temp_c
            params[2],  // base_press_hpa
            params[3],  // base_ratio
        ))
    } else {
        Err(pyo3::exceptions::PyValueError::new_err(
            format!("Atmosphere params must contain 2 or 4 values, got {}", params.len())
        ))
    }
}

/// Python wrapper for stability coefficient calculation
#[pyfunction]
fn compute_stability_coefficient_rust(
    inputs: &Bound<'_, PyDict>,
    atmos_params: PyReadonlyArray1<f64>,
) -> PyResult<f64> {
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    let params = atmos_params.as_array();
    
    if params.len() != 4 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            format!("Atmosphere params must contain 4 values, got {}", params.len())
        ));
    }
    
    let atmo_params_tuple = (params[0], params[1], params[2], params[3]);
    Ok(compute_stability_coefficient(&ballistic_inputs, atmo_params_tuple))
}

/// Python wrapper for spin drift calculation
#[pyfunction]
fn compute_spin_drift_rust(
    time_s: f64,
    stability: f64,
    twist_rate: f64,
    is_twist_right: bool,
) -> PyResult<f64> {
    Ok(compute_spin_drift(time_s, stability, twist_rate, is_twist_right))
}

/// Python wrapper for trajectory solver initial conditions preparation
#[pyfunction]
fn prepare_initial_conditions_rust(
    inputs: &Bound<'_, PyDict>,
    zero_angle_rad: f64,
    atmos_params: PyReadonlyArray1<f64>,
    air_density: f64,
    speed_of_sound: f64,
    stability_coefficient: f64,
) -> PyResult<PyObject> {
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    let params = atmos_params.as_array();
    
    if params.len() != 4 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            format!("Atmosphere params must contain 4 values, got {}", params.len())
        ));
    }
    
    let atmo_params_tuple = (params[0], params[1], params[2], params[3]);
    let conditions = prepare_initial_conditions(
        &ballistic_inputs,
        zero_angle_rad,
        atmo_params_tuple,
        air_density,
        speed_of_sound,
        stability_coefficient,
    );
    
    // Convert to Python dict
    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new_bound(py);
        dict.set_item("mass_kg", conditions.mass_kg)?;
        dict.set_item("muzzle_velocity_mps", conditions.muzzle_velocity_mps)?;
        dict.set_item("target_distance_los_m", conditions.target_distance_los_m)?;
        dict.set_item("muzzle_angle_rad", conditions.muzzle_angle_rad)?;
        dict.set_item("muzzle_energy_j", conditions.muzzle_energy_j)?;
        dict.set_item("muzzle_energy_ftlbs", conditions.muzzle_energy_ftlbs)?;
        dict.set_item("target_horizontal_dist_m", conditions.target_horizontal_dist_m)?;
        dict.set_item("target_vertical_height_m", conditions.target_vertical_height_m)?;
        dict.set_item("initial_state", conditions.initial_state.to_pyarray_bound(py))?;
        dict.set_item("t_span", (conditions.t_span.0, conditions.t_span.1))?;
        dict.set_item("stability_coefficient", conditions.stability_coefficient)?;
        dict.set_item("air_density", conditions.air_density)?;
        dict.set_item("speed_of_sound", conditions.speed_of_sound)?;
        
        if let Some(omega) = conditions.omega_vector {
            dict.set_item("omega_vector", [omega.x, omega.y, omega.z].to_pyarray_bound(py))?;
        } else {
            dict.set_item("omega_vector", py.None())?;
        }
        
        Ok(dict.into())
    })
}

/// Python wrapper for trajectory post-processing
#[pyfunction]
#[pyo3(signature = (trajectory_points, initial_conditions, inputs, target_hit_time=None, ground_hit_time=None))]
fn post_process_trajectory_rust(
    py: Python,
    trajectory_points: PyReadonlyArray1<f64>, // Flattened array: [t1, x1, y1, z1, vx1, vy1, vz1, t2, ...]
    initial_conditions: &Bound<'_, PyDict>,
    inputs: &Bound<'_, PyDict>,
    target_hit_time: Option<f64>,
    ground_hit_time: Option<f64>,
) -> PyResult<PyObject> {
    let points_array = trajectory_points.as_array();
    let n_points = points_array.len() / 7; // 7 values per point: time + 6 state values
    
    if points_array.len() % 7 != 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Trajectory points array length must be divisible by 7"
        ));
    }
    
    // Convert flattened array to vector of (time, state) tuples
    let mut trajectory_points_vec = Vec::new();
    for i in 0..n_points {
        let base_idx = i * 7;
        let time = points_array[base_idx];
        let state = [
            points_array[base_idx + 1], // x
            points_array[base_idx + 2], // y
            points_array[base_idx + 3], // z
            points_array[base_idx + 4], // vx
            points_array[base_idx + 5], // vy
            points_array[base_idx + 6], // vz
        ];
        trajectory_points_vec.push((time, state));
    }
    
    // Extract required values from initial_conditions dict
    let mass_kg: f64 = initial_conditions.get_item("mass_kg")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing mass_kg in initial_conditions"))?
        .extract()?;
    let target_distance_los_m: f64 = initial_conditions.get_item("target_distance_los_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing target_distance_los_m in initial_conditions"))?
        .extract()?;
    let target_horizontal_dist_m: f64 = initial_conditions.get_item("target_horizontal_dist_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing target_horizontal_dist_m in initial_conditions"))?
        .extract()?;
    let target_vertical_height_m: f64 = initial_conditions.get_item("target_vertical_height_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing target_vertical_height_m in initial_conditions"))?
        .extract()?;
    let muzzle_energy_j: f64 = initial_conditions.get_item("muzzle_energy_j")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing muzzle_energy_j in initial_conditions"))?
        .extract()?;
    let muzzle_energy_ftlbs: f64 = initial_conditions.get_item("muzzle_energy_ftlbs")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing muzzle_energy_ftlbs in initial_conditions"))?
        .extract()?;
    let muzzle_angle_rad: f64 = initial_conditions.get_item("muzzle_angle_rad")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing muzzle_angle_rad in initial_conditions"))?
        .extract()?;
    let stability_coefficient: f64 = initial_conditions.get_item("stability_coefficient")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing stability_coefficient in initial_conditions"))?
        .extract()?;
    let air_density: f64 = initial_conditions.get_item("air_density")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing air_density in initial_conditions"))?
        .extract()?;
    let speed_of_sound: f64 = initial_conditions.get_item("speed_of_sound")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing speed_of_sound in initial_conditions"))?
        .extract()?;
    
    // Create a simplified InitialConditions struct for post-processing
    let conditions = trajectory_solver::InitialConditions {
        mass_kg,
        muzzle_velocity_mps: 0.0, // Not needed for post-processing
        target_distance_los_m,
        muzzle_angle_rad,
        muzzle_energy_j,
        muzzle_energy_ftlbs,
        target_horizontal_dist_m,
        target_vertical_height_m,
        initial_state: [0.0; 6], // Not needed for post-processing
        t_span: (0.0, 0.0), // Not needed for post-processing
        omega_vector: None, // Not needed for post-processing
        stability_coefficient,
        atmo_params: (0.0, 0.0, 0.0, 0.0), // Not needed for post-processing
        air_density,
        speed_of_sound,
    };
    
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    
    let result = post_process_trajectory(
        &trajectory_points_vec,
        &conditions,
        &ballistic_inputs,
        target_hit_time,
        ground_hit_time,
    ).map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    
    // Convert result to Python dict
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("muzzle_energy_j", result.muzzle_energy_j)?;
    dict.set_item("muzzle_energy_ftlbs", result.muzzle_energy_ftlbs)?;
    dict.set_item("target_distance_los_m", result.target_distance_los_m)?;
    dict.set_item("target_distance_horiz_m", result.target_distance_horiz_m)?;
    dict.set_item("target_vertical_height_m", result.target_vertical_height_m)?;
    dict.set_item("time_of_flight_s", result.time_of_flight_s)?;
    dict.set_item("drop_m", result.drop_m)?;
    dict.set_item("drop_in", result.drop_in)?;
    dict.set_item("wind_drift_m", result.wind_drift_m)?;
    dict.set_item("wind_drift_in", result.wind_drift_in)?;
    dict.set_item("max_ord_m", result.max_ord_m)?;
    dict.set_item("max_ord_in", result.max_ord_in)?;
    dict.set_item("max_ord_dist_horiz_m", result.max_ord_dist_horiz_m)?;
    dict.set_item("final_vel_mps", result.final_vel_mps)?;
    dict.set_item("final_vel_fps", result.final_vel_fps)?;
    dict.set_item("final_energy_j", result.final_energy_j)?;
    dict.set_item("final_energy_ftlbs", result.final_energy_ftlbs)?;
    dict.set_item("air_density_kg_m3", result.air_density_kg_m3)?;
    dict.set_item("speed_of_sound_mps", result.speed_of_sound_mps)?;
    dict.set_item("barrel_angle_rad", result.barrel_angle_rad)?;
    
    Ok(dict.into())
}

/// Python wrapper for trajectory sampling
#[pyfunction]
fn sample_trajectory_rust(
    py: Python,
    trajectory_times: PyReadonlyArray1<f64>,
    trajectory_positions: PyReadonlyArray1<f64>, // Flattened [x1, y1, z1, x2, y2, z2, ...]
    trajectory_velocities: PyReadonlyArray1<f64>, // Flattened [vx1, vy1, vz1, vx2, vy2, vz2, ...]
    transonic_distances: PyReadonlyArray1<f64>,
    outputs_dict: &Bound<'_, PyDict>,
    step_m: f64,
    mass_kg: f64,
) -> PyResult<PyObject> {
    let times_array = trajectory_times.as_array();
    let positions_array = trajectory_positions.as_array();
    let velocities_array = trajectory_velocities.as_array();
    let transonic_array = transonic_distances.as_array();
    
    // Validate input dimensions
    let n_points = times_array.len();
    if positions_array.len() != n_points * 3 || velocities_array.len() != n_points * 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "Position and velocity arrays must have 3 * n_points elements"
        ));
    }
    
    // Convert flattened arrays to vector of Vector3
    let mut positions = Vec::with_capacity(n_points);
    let mut velocities = Vec::with_capacity(n_points);
    
    for i in 0..n_points {
        let pos_idx = i * 3;
        let vel_idx = i * 3;
        
        positions.push(Vector3::new(
            positions_array[pos_idx],
            positions_array[pos_idx + 1],
            positions_array[pos_idx + 2],
        ));
        
        velocities.push(Vector3::new(
            velocities_array[vel_idx],
            velocities_array[vel_idx + 1],
            velocities_array[vel_idx + 2],
        ));
    }
    
    // Extract outputs from dictionary
    let target_distance_horiz_m: f64 = outputs_dict.get_item("target_distance_horiz_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing target_distance_horiz_m in outputs"))?
        .extract()?;
    let target_vertical_height_m: f64 = outputs_dict.get_item("target_vertical_height_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing target_vertical_height_m in outputs"))?
        .extract()?;
    let time_of_flight_s: f64 = outputs_dict.get_item("time_of_flight_s")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing time_of_flight_s in outputs"))?
        .extract()?;
    let max_ord_dist_horiz_m: f64 = outputs_dict.get_item("max_ord_dist_horiz_m")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing max_ord_dist_horiz_m in outputs"))?
        .extract()?;
    
    // Create trajectory data structure
    let trajectory_data = TrajectoryData {
        times: times_array.to_vec(),
        positions,
        velocities,
        transonic_distances: transonic_array.to_vec(),
    };
    
    let outputs = TrajectoryOutputs {
        target_distance_horiz_m,
        target_vertical_height_m,
        time_of_flight_s,
        max_ord_dist_horiz_m,
    };
    
    // Perform trajectory sampling
    let samples = sample_trajectory(&trajectory_data, &outputs, step_m, mass_kg);
    
    // Convert to Python-compatible format
    let sample_dicts = trajectory_samples_to_dicts(&samples);
    
    // Convert to Python list of dictionaries
    let py_list = pyo3::types::PyList::empty_bound(py);
    
    for sample_dict in &sample_dicts {
        let py_dict = pyo3::types::PyDict::new_bound(py);
        py_dict.set_item("distance_m", sample_dict.distance_m)?;
        py_dict.set_item("drop_m", sample_dict.drop_m)?;
        py_dict.set_item("wind_drift_m", sample_dict.wind_drift_m)?;
        py_dict.set_item("velocity_mps", sample_dict.velocity_mps)?;
        py_dict.set_item("energy_j", sample_dict.energy_j)?;
        py_dict.set_item("time_s", sample_dict.time_s)?;
        py_dict.set_item("flags", &sample_dict.flags)?;
        
        py_list.append(py_dict)?;
    }
    
    Ok(py_list.into())
}

/// Python wrapper for Brent's method root finding
#[pyfunction]
fn brent_root_find_rust(
    py: Python,
    func_name: String,
    a: f64,
    b: f64,
    tolerance: f64,
    max_iterations: usize,
) -> PyResult<PyObject> {
    // This is a simplified wrapper - in practice, we'd need to pass Python functions
    // For now, we'll provide a few predefined test functions
    let f = match func_name.as_str() {
        "quadratic" => |x: f64| x * x - 4.0,
        "linear" => |x: f64| 2.0 * x - 6.0,
        _ => return Err(pyo3::exceptions::PyValueError::new_err("Unknown function name")),
    };
    
    match brent_root_find(f, a, b, tolerance, max_iterations) {
        Ok(result) => {
            let dict = pyo3::types::PyDict::new_bound(py);
            dict.set_item("angle_rad", result.angle_rad)?;
            dict.set_item("iterations_used", result.iterations_used)?;
            dict.set_item("final_error", result.final_error)?;
            dict.set_item("success", result.success)?;
            Ok(dict.into())
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Python wrapper for adjusted muzzle velocity calculation
#[pyfunction]
fn adjusted_muzzle_velocity_rust(inputs: &Bound<'_, PyDict>) -> PyResult<f64> {
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    Ok(adjusted_muzzle_velocity(&ballistic_inputs))
}

/// Simplified Python wrapper for zero angle calculation (for testing)
#[pyfunction]
fn zero_angle_simplified_rust(
    inputs: &Bound<'_, PyDict>,
    _target_height_m: f64,
) -> PyResult<f64> {
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    
    // Simple trajectory function for testing - returns height at target distance
    let trajectory_func = |_inputs: &BallisticInputs, angle_rad: f64| -> Result<f64, String> {
        // Simplified ballistic calculation for demonstration
        let target_distance_m = _inputs.target_distance * 0.9144; // yards to meters
        let muzzle_velocity_mps = _inputs.muzzle_velocity * 0.3048; // fps to mps
        
        let time_of_flight = target_distance_m / (muzzle_velocity_mps * angle_rad.cos());
        let height = muzzle_velocity_mps * angle_rad.sin() * time_of_flight 
                   - 0.5 * 9.80665 * time_of_flight * time_of_flight;
        
        Ok(height)
    };
    
    match zero_angle(&ballistic_inputs, trajectory_func) {
        Ok(result) => Ok(result.angle_rad),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Python wrapper for solve muzzle angle (simplified for testing)
#[pyfunction]
fn solve_muzzle_angle_rust(
    inputs: &Bound<'_, PyDict>,
    zero_distance_los_m: f64,
    angle_lower_deg: f64,
    angle_upper_deg: f64,
    rtol: f64,
) -> PyResult<PyObject> {
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    
    // Simple trajectory function that returns drop_m
    let trajectory_func = |_inputs: &BallisticInputs| -> Result<f64, String> {
        // Simplified drop calculation for demonstration
        let target_distance_m = _inputs.target_distance * 0.9144; // yards to meters
        let muzzle_velocity_mps = _inputs.muzzle_velocity * 0.3048; // fps to mps
        let angle_rad = _inputs.muzzle_angle * std::f64::consts::PI / 180.0;
        
        let time_of_flight = target_distance_m / (muzzle_velocity_mps * angle_rad.cos());
        let height = muzzle_velocity_mps * angle_rad.sin() * time_of_flight 
                   - 0.5 * 9.80665 * time_of_flight * time_of_flight;
        
        // Return drop (negative height)
        Ok(-height)
    };
    
    match solve_muzzle_angle(
        &ballistic_inputs,
        zero_distance_los_m,
        trajectory_func,
        angle_lower_deg,
        angle_upper_deg,
        rtol,
    ) {
        Ok(result) => {
            Python::with_gil(|py| {
                let dict = pyo3::types::PyDict::new_bound(py);
                dict.set_item("angle_rad", result.angle_rad)?;
                dict.set_item("iterations_used", result.iterations_used)?;
                dict.set_item("final_error", result.final_error)?;
                dict.set_item("success", result.success)?;
                Ok(dict.into())
            })
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

/// Python wrapper for WindSock
#[pyclass]
struct PyWindSock {
    inner: WindSock,
}

#[pymethods]
impl PyWindSock {
    #[new]
    fn new(segments: Vec<(f64, f64, f64)>) -> Self {
        PyWindSock {
            inner: WindSock::new(segments),
        }
    }
    
    /// Get wind vector for a given range (stateless version)
    fn vector_for_range(&self, range_m: f64) -> PyResult<Py<PyArray1<f64>>> {
        Python::with_gil(|py| {
            let vec = self.inner.vector_for_range_stateless(range_m);
            let array = [vec.x, vec.y, vec.z];
            Ok(array.to_pyarray_bound(py).into())
        })
    }
}

/// Python wrapper for WindShearWindSock
#[pyclass]
struct PyWindShearWindSock {
    inner: wind_shear::WindShearWindSock,
}

#[pymethods]
impl PyWindShearWindSock {
    #[new]
    fn new(segments: Vec<(f64, f64, f64)>, enable_shear: bool, shear_model: String, shooter_altitude_m: Option<f64>) -> Self {
        let shear_profile = if enable_shear && shear_model != "none" {
            let mut profile = WindShearProfile::default();
            profile.model = match shear_model.as_str() {
                "logarithmic" => WindShearModel::Logarithmic,
                "power_law" => WindShearModel::PowerLaw,
                "ekman_spiral" => WindShearModel::EkmanSpiral,
                _ => WindShearModel::None,
            };
            // Set default surface wind from the first segment
            // Note: The wind speed should be the measurement at reference_height (typically 10m)
            // The altitude_m field is not used in calculations, only the wind vector
            if !segments.is_empty() {
                profile.surface_wind = WindLayer {
                    altitude_m: 0.0,  // Not used in wind shear calculations
                    speed_mps: segments[0].0 * 0.2777778, // km/h to m/s at reference_height
                    direction_deg: segments[0].1,
                };
            }
            Some(profile)
        } else {
            None
        };
        
        PyWindShearWindSock {
            inner: wind_shear::WindShearWindSock::with_shooter_altitude(
                segments, 
                shear_profile, 
                shooter_altitude_m.unwrap_or(0.0)
            ),
        }
    }
    
    /// Get wind vector for a given 3D position
    fn vector_for_position(&self, py: Python, x: f64, y: f64, z: f64) -> PyResult<Py<PyArray1<f64>>> {
        let pos = Vector3::new(x, y, z);
        let vec = self.inner.vector_for_position(pos);
        let array = [vec.x, vec.y, vec.z];
        Ok(array.to_pyarray_bound(py).into())
    }
}

/// Python wrapper for fast trajectory integration
#[pyfunction]
fn fast_integrate_rust(
    py: Python,
    inputs: &Bound<'_, PyDict>,
    wind_segments: Vec<WindSegment>,
    horiz: f64,
    vert: f64,
    initial_state: PyReadonlyArray1<f64>,
    t_span: (f64, f64),
    atmo_params: PyReadonlyArray1<f64>,
) -> PyResult<PyObject> {
    use fast_trajectory::{fast_integrate, FastIntegrationParams};
    
    // Extract inputs
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    let wind_sock = WindSock::new(wind_segments);
    
    // Convert arrays
    let initial_state_array = initial_state.as_array();
    let atmo_params_array = atmo_params.as_array();
    
    if initial_state_array.len() != 6 {
        return Err(pyo3::exceptions::PyValueError::new_err("Initial state must have 6 elements"));
    }
    
    if atmo_params_array.len() != 4 {
        return Err(pyo3::exceptions::PyValueError::new_err("Atmospheric parameters must have 4 elements"));
    }
    
    let mut initial_state_arr = [0.0; 6];
    for i in 0..6 {
        initial_state_arr[i] = initial_state_array[i];
    }
    
    let params = FastIntegrationParams {
        horiz,
        vert,
        initial_state: initial_state_arr,
        t_span,
        atmo_params: (
            atmo_params_array[0],
            atmo_params_array[1],
            atmo_params_array[2],
            atmo_params_array[3],
        ),
    };
    
    // Run integration
    let solution = fast_integrate(&ballistic_inputs, &wind_sock, params);
    
    // Convert to Python dictionary
    let dict = pyo3::types::PyDict::new_bound(py);
    
    // Convert times
    dict.set_item("t", solution.t.to_pyarray_bound(py))?;
    
    // Convert y matrix (6 x n_points) to numpy array
    let n_points = solution.t.len();
    // Convert the row-major solution.y to a 2D numpy array
    let y_array = numpy::PyArray2::from_vec2_bound(py, &solution.y).unwrap();
    dict.set_item("y", y_array)?;
    
    // Convert t_events
    let t_events_list = pyo3::types::PyList::empty_bound(py);
    for events in &solution.t_events {
        t_events_list.append(events.to_pyarray_bound(py))?;
    }
    dict.set_item("t_events", t_events_list)?;
    
    dict.set_item("success", solution.success)?;
    
    Ok(dict.into())
}

/// Python wrapper for BC segment estimation
#[pyfunction]
#[pyo3(signature = (base_bc, caliber, weight, model, drag_model="G1"))]
fn estimate_bc_segments_rust(
    py: Python,
    base_bc: f64,
    caliber: f64,
    weight: f64,
    model: &str,
    drag_model: &str,
) -> PyResult<PyObject> {
    // Estimate BC segments
    let segments = BCSegmentEstimator::estimate_bc_segments(
        base_bc,
        caliber,
        weight,
        model,
        drag_model,
    );
    
    // Convert to Python list of dicts
    let py_list = pyo3::types::PyList::empty_bound(py);
    
    for segment in segments {
        let py_dict = pyo3::types::PyDict::new_bound(py);
        py_dict.set_item("velocity_min", segment.velocity_min)?;
        py_dict.set_item("velocity_max", segment.velocity_max)?;
        py_dict.set_item("bc_value", segment.bc_value)?;
        py_list.append(py_dict)?;
    }
    
    Ok(py_list.into())
}

/// Python wrapper for spin decay calculation
#[pyfunction]
fn calculate_spin_decay_rust(
    initial_spin_rad_s: f64,
    time_elapsed_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    mass_grains: f64,
    caliber_inches: f64,
    length_inches: f64,
    bullet_type: &str,
) -> PyResult<f64> {
    let params = SpinDecayParameters::from_bullet_type(bullet_type);
    
    let current_spin = spin_decay::update_spin_rate(
        initial_spin_rad_s,
        time_elapsed_s,
        velocity_mps,
        air_density_kg_m3,
        mass_grains,
        caliber_inches,
        length_inches,
        Some(&params),
    );
    
    Ok(current_spin)
}

/// Python wrapper for spin decay correction factor
#[pyfunction]
fn calculate_spin_decay_correction_rust(
    time_elapsed_s: f64,
    velocity_mps: f64,
    air_density_kg_m3: f64,
    mass_grains: f64,
    caliber_inches: f64,
    length_inches: f64,
    bullet_type: &str,
) -> PyResult<f64> {
    let params = SpinDecayParameters::from_bullet_type(bullet_type);
    
    let correction = spin_decay::calculate_spin_decay_correction_factor(
        time_elapsed_s,
        velocity_mps,
        air_density_kg_m3,
        mass_grains,
        caliber_inches,
        length_inches,
        Some(&params),
    );
    
    Ok(correction)
}

/// Python wrapper for epicyclic motion calculation
#[pyfunction]
fn calculate_epicyclic_motion_rust(
    py: Python,
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    stability_factor: f64,
    time_s: f64,
    initial_yaw_rad: f64,
) -> PyResult<PyObject> {
    let (pitch, yaw) = precession_nutation::calculate_epicyclic_motion(
        spin_rate_rad_s,
        velocity_mps,
        stability_factor,
        time_s,
        initial_yaw_rad,
    );
    
    // Return as dictionary
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("pitch_rad", pitch)?;
    dict.set_item("yaw_rad", yaw)?;
    dict.set_item("pitch_moa", pitch * 3437.7467707849)?; // rad to MOA
    dict.set_item("yaw_moa", yaw * 3437.7467707849)?;
    
    Ok(dict.into())
}

/// Python wrapper for limit cycle yaw calculation
#[pyfunction]
fn calculate_limit_cycle_yaw_rust(
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    stability_factor: f64,
    crosswind_mps: f64,
) -> PyResult<f64> {
    let yaw = precession_nutation::calculate_limit_cycle_yaw(
        velocity_mps,
        spin_rate_rad_s,
        stability_factor,
        crosswind_mps,
    );
    
    Ok(yaw)
}

/// Python wrapper for aerodynamic jump calculation
#[pyfunction]
fn calculate_aerodynamic_jump_rust(
    py: Python,
    inputs: &Bound<'_, PyDict>,
    crosswind_mps: f64,
    barrel_length_m: f64,
    air_density_kg_m3: f64,
) -> PyResult<PyObject> {
    // Extract inputs
    let ballistic_inputs = extract_ballistic_inputs(inputs)?;
    
    // Calculate spin rate
    let muzzle_velocity_mps = ballistic_inputs.muzzle_velocity * 0.3048; // fps to mps
    let spin_rate_rad_s = if ballistic_inputs.twist_rate > 0.0 {
        let velocity_fps = ballistic_inputs.muzzle_velocity;
        let twist_rate_ft = ballistic_inputs.twist_rate / 12.0;
        let revolutions_per_second = velocity_fps / twist_rate_ft;
        revolutions_per_second * 2.0 * std::f64::consts::PI
    } else {
        0.0
    };
    
    // Convert to metric
    let caliber_m = ballistic_inputs.bullet_diameter * 0.0254; // inches to meters
    let mass_kg = ballistic_inputs.bullet_mass * 0.00006479891; // grains to kg
    let twist_rate_calibers = ballistic_inputs.twist_rate / ballistic_inputs.bullet_diameter;
    
    // Calculate aerodynamic jump
    let jump = calculate_aerodynamic_jump(
        muzzle_velocity_mps,
        spin_rate_rad_s,
        crosswind_mps,
        caliber_m,
        mass_kg,
        barrel_length_m,
        twist_rate_calibers,
        ballistic_inputs.is_twist_right,
        ballistic_inputs.tipoff_yaw,
        air_density_kg_m3,
    );
    
    // Convert to Python dict
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("vertical_jump_moa", jump.vertical_jump_moa)?;
    dict.set_item("horizontal_jump_moa", jump.horizontal_jump_moa)?;
    dict.set_item("jump_angle_rad", jump.jump_angle_rad)?;
    dict.set_item("magnus_component_moa", jump.magnus_component_moa)?;
    dict.set_item("yaw_component_moa", jump.yaw_component_moa)?;
    dict.set_item("stabilization_factor", jump.stabilization_factor)?;
    
    Ok(dict.into())
}

/// Python module definition
#[pymodule]
fn ballistics_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(derivatives_rust, m)?)?;
    m.add_function(wrap_pyfunction!(get_drag_coefficient_rust, m)?)?;
    m.add_function(wrap_pyfunction!(get_drag_coefficient_rust_transonic, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_atmosphere_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_air_density_cipm_rust, m)?)?;
    m.add_function(wrap_pyfunction!(get_local_atmosphere_rust, m)?)?;
    m.add_function(wrap_pyfunction!(compute_stability_coefficient_rust, m)?)?;
    m.add_function(wrap_pyfunction!(compute_spin_drift_rust, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_initial_conditions_rust, m)?)?;
    m.add_function(wrap_pyfunction!(post_process_trajectory_rust, m)?)?;
    m.add_function(wrap_pyfunction!(sample_trajectory_rust, m)?)?;
    m.add_function(wrap_pyfunction!(adjusted_muzzle_velocity_rust, m)?)?;
    m.add_function(wrap_pyfunction!(zero_angle_simplified_rust, m)?)?;
    m.add_function(wrap_pyfunction!(solve_muzzle_angle_rust, m)?)?;
    m.add_function(wrap_pyfunction!(fast_integrate_rust, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_bc_segments_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_aerodynamic_jump_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_spin_decay_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_spin_decay_correction_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_epicyclic_motion_rust, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_limit_cycle_yaw_rust, m)?)?;
    
    // Monte Carlo functions
    m.add_function(wrap_pyfunction!(monte_carlo_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(monte_carlo_batch, m)?)?;
    m.add_function(wrap_pyfunction!(monte_carlo_streaming, m)?)?;
    
    // Aliases for easier access
    m.add("get_drag_coefficient", m.getattr("get_drag_coefficient_rust")?)?;
    m.add("calculate_atmosphere", m.getattr("calculate_atmosphere_rust")?)?;
    m.add("calculate_air_density_cipm", m.getattr("calculate_air_density_cipm_rust")?)?;
    m.add("get_local_atmosphere", m.getattr("get_local_atmosphere_rust")?)?;
    m.add("compute_stability_coefficient", m.getattr("compute_stability_coefficient_rust")?)?;
    m.add("compute_spin_drift", m.getattr("compute_spin_drift_rust")?)?;
    m.add("prepare_initial_conditions", m.getattr("prepare_initial_conditions_rust")?)?;
    m.add("post_process_trajectory", m.getattr("post_process_trajectory_rust")?)?;
    m.add("sample_trajectory", m.getattr("sample_trajectory_rust")?)?;
    m.add("adjusted_muzzle_velocity", m.getattr("adjusted_muzzle_velocity_rust")?)?;
    m.add("zero_angle_simplified", m.getattr("zero_angle_simplified_rust")?)?;
    m.add("solve_muzzle_angle", m.getattr("solve_muzzle_angle_rust")?)?;
    m.add("fast_integrate", m.getattr("fast_integrate_rust")?)?;
    m.add("estimate_bc_segments", m.getattr("estimate_bc_segments_rust")?)?;
    m.add("calculate_aerodynamic_jump", m.getattr("calculate_aerodynamic_jump_rust")?)?;
    m.add("calculate_spin_decay", m.getattr("calculate_spin_decay_rust")?)?;
    m.add("calculate_spin_decay_correction", m.getattr("calculate_spin_decay_correction_rust")?)?;
    m.add("calculate_epicyclic_motion", m.getattr("calculate_epicyclic_motion_rust")?)?;
    m.add("calculate_limit_cycle_yaw", m.getattr("calculate_limit_cycle_yaw_rust")?)?;
    
    // Add WindSock class
    m.add_class::<PyWindSock>()?;
    // Add WindShearWindSock class
    m.add_class::<PyWindShearWindSock>()?;
    
    // Add transonic drag functions
    m.add_function(wrap_pyfunction!(transonic_correction_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_projectile_shape_py, m)?)?;
    
    // Add Reynolds correction function
    m.add_function(wrap_pyfunction!(apply_reynolds_correction_py, m)?)?;
    
    // Add BC interpolation function
    m.add_function(wrap_pyfunction!(interpolated_bc_py, m)?)?;
    
    Ok(())
}