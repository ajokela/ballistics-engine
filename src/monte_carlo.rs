use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::collections::HashMap;
use numpy::PyReadonlyArray2;

use crate::InternalBallisticInputs as BallisticInputs;
use crate::fast_trajectory::{fast_integrate, FastIntegrationParams};
use crate::wind::WindSock;
use crate::atmosphere::calculate_atmosphere;
use crate::constants::{FPS_TO_MPS, GRAINS_TO_KG};
use nalgebra::Vector3;

// Define constants that might be missing
const YARDS_TO_METERS: f64 = 0.9144;

/// Configure Rayon thread pool with proper error handling
/// 
/// This function centralizes thread pool configuration logic to ensure consistent
/// error handling across all Monte Carlo functions.
/// 
/// # Arguments
/// * `num_threads` - Optional number of threads to configure. If None, uses default.
/// 
/// # Returns
/// * `PyResult<()>` - Returns error if thread count is invalid (0), otherwise succeeds
/// 
/// # Behavior
/// - Validates thread count is greater than 0
/// - Attempts to configure global Rayon thread pool
/// - On failure, logs warning but continues with default threading
/// - This ensures Monte Carlo functions remain functional even if thread configuration fails
fn configure_thread_pool(num_threads: Option<usize>) -> PyResult<()> {
    if let Some(n) = num_threads {
        // Validate thread count before attempting to set it
        if n == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "Thread count must be greater than 0"
            ));
        }
        
        // Attempt to configure the global thread pool
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global() 
        {
            // Log warning but continue execution with default threading
            // This approach ensures the Monte Carlo calculations remain functional
            // even if the specific thread configuration cannot be applied
            eprintln!("Warning: Failed to set {n} threads, using default threading: {e}");
        }
    }
    Ok(())
}

/// Simple trajectory output for Monte Carlo
#[derive(Debug, Clone)]
pub struct TrajectoryOutput {
    pub drop: f64,           // meters
    pub wind_drift: f64,      // meters
    pub time: f64,            // seconds
    pub velocity: f64,        // m/s
    pub energy: f64,          // joules
    pub mach: f64,            // mach number
    pub spin_drift: f64,      // meters
    pub distance: f64,        // meters
}

/// Solve trajectory for Monte Carlo run
fn solve_trajectory_for_monte_carlo(inputs: &BallisticInputs) -> Result<TrajectoryOutput, String> {
    // Convert inputs to metric
    let target_distance_m = inputs.target_distance * YARDS_TO_METERS;
    let muzzle_velocity_mps = inputs.muzzle_velocity * FPS_TO_MPS;
    let mass_kg = inputs.bullet_mass * GRAINS_TO_KG;
    
    // Calculate atmosphere at altitude
    let (air_density, speed_of_sound) = calculate_atmosphere(
        inputs.altitude * 0.3048,  // feet to meters
        Some(inputs.temperature),
        Some(inputs.pressure),
        inputs.humidity,
    );
    
    // Create wind segments as tuples (the WindSock expects tuples)
    let wind_segments = vec![(
        0.0,  // range_m
        inputs.wind_speed * 0.44704,  // wind_speed (mph to m/s)
        inputs.wind_angle,  // wind_angle
    )];
    let wind_sock = WindSock::new(wind_segments);
    
    // Set up initial state
    let angle_rad = inputs.muzzle_angle.to_radians();
    let initial_state = [
        0.0, 0.0, 0.0,  // position
        muzzle_velocity_mps * angle_rad.cos(),  // vx
        muzzle_velocity_mps * angle_rad.sin(),  // vy
        0.0,  // vz
    ];
    
    // Time span (estimate)
    let t_span = (0.0, target_distance_m / (muzzle_velocity_mps * 0.5));
    
    // Atmospheric parameters
    let atmo_params = (
        inputs.altitude * 0.3048,
        inputs.temperature,
        inputs.pressure,
        air_density / 1.225,  // density ratio
    );
    
    // Fast integration parameters
    let params = FastIntegrationParams {
        horiz: target_distance_m,
        vert: 0.0,  // Simplified - no vertical target
        initial_state,
        t_span,
        atmo_params,
    };
    
    // Run integration
    let solution = fast_integrate(inputs, &wind_sock, params);
    
    if !solution.success || solution.t.is_empty() {
        return Err("Trajectory integration failed".to_string());
    }
    
    // Extract final values
    let final_idx = solution.t.len() - 1;
    let final_pos = Vector3::new(
        solution.y[0][final_idx],
        solution.y[1][final_idx],
        solution.y[2][final_idx],
    );
    let final_vel = Vector3::new(
        solution.y[3][final_idx],
        solution.y[4][final_idx],
        solution.y[5][final_idx],
    );
    
    let velocity = final_vel.norm();
    let energy = 0.5 * mass_kg * velocity * velocity;
    let mach = velocity / speed_of_sound;
    
    Ok(TrajectoryOutput {
        drop: -final_pos.y,  // Negative y is drop
        wind_drift: final_pos.z,
        time: solution.t[final_idx],
        velocity,
        energy,
        mach,
        spin_drift: 0.0,  // Simplified
        distance: final_pos.x,
    })
}

/// Statistics for a single output field
#[derive(Debug, Clone)]
pub struct FieldStatistics {
    pub mean: f64,
    pub std: f64,
    pub min: f64,
    pub max: f64,
    pub percentiles: Vec<(f64, f64)>,  // (percentile, value) pairs
}

/// Results from Monte Carlo simulation
#[derive(Debug)]
pub struct MonteCarloResults {
    pub statistics: HashMap<String, FieldStatistics>,
    pub valid_runs: usize,
    pub failed_runs: usize,
}

/// Parallel Monte Carlo trajectory evaluation
/// 
/// This function takes pre-generated parameter samples from Python and evaluates
/// trajectories in parallel using Rayon. It focuses solely on the computationally
/// intensive trajectory calculations, leaving statistical modeling to Python.
#[pyfunction]
#[pyo3(signature = (base_inputs, param_samples, param_names, num_threads=None))]
pub fn monte_carlo_parallel(
    py: Python,
    base_inputs: &Bound<'_, PyDict>,  // Dictionary of inputs
    param_samples: PyReadonlyArray2<f64>,  // 2D array: [n_samples x n_params]
    param_names: Vec<String>,        // Parameter names corresponding to columns
    num_threads: Option<usize>,
) -> PyResult<PyObject> {
    // Configure thread pool with proper error handling
    configure_thread_pool(num_threads)?;
    
    // Extract base inputs from dictionary
    let base_inputs = extract_ballistic_inputs_from_dict(base_inputs)?;
    
    // Convert numpy array to ndarray for easier manipulation
    let samples_array = param_samples.as_array();
    let n_samples = samples_array.shape()[0];
    
    // Pre-allocate results storage
    let trajectory_results: Vec<Option<TrajectoryOutput>> = (0..n_samples)
        .into_par_iter()
        .map(|i| {
            // Clone base inputs for this run
            let mut run_inputs = base_inputs.clone();
            
            // Apply parameter perturbations
            for (j, param_name) in param_names.iter().enumerate() {
                let value = samples_array[[i, j]];
                apply_parameter(&mut run_inputs, param_name, value);
            }
            
            // Solve trajectory
            match solve_trajectory_for_monte_carlo(&run_inputs) {
                Ok(output) => Some(output),
                Err(_) => None,  // Failed run
            }
        })
        .collect();
    
    // Calculate statistics
    let (statistics, valid_runs, failed_runs) = calculate_statistics(&trajectory_results);
    
    // Convert to Python dictionary
    let py_dict = PyDict::new_bound(py);
    
    // Add statistics
    let stats_dict = PyDict::new_bound(py);
    for (field_name, field_stats) in statistics {
        let field_dict = PyDict::new_bound(py);
        field_dict.set_item("mean", field_stats.mean)?;
        field_dict.set_item("std", field_stats.std)?;
        field_dict.set_item("min", field_stats.min)?;
        field_dict.set_item("max", field_stats.max)?;
        
        let percentiles_dict = PyDict::new_bound(py);
        for (p, v) in field_stats.percentiles {
            percentiles_dict.set_item(format!("p{}", (p * 100.0) as i32), v)?;
        }
        field_dict.set_item("percentiles", percentiles_dict)?;
        
        stats_dict.set_item(field_name, field_dict)?;
    }
    py_dict.set_item("statistics", stats_dict)?;
    
    // Add metadata
    let metadata = PyDict::new_bound(py);
    metadata.set_item("valid_runs", valid_runs)?;
    metadata.set_item("failed_runs", failed_runs)?;
    metadata.set_item("total_runs", n_samples)?;
    metadata.set_item("success_rate", valid_runs as f64 / n_samples as f64)?;
    py_dict.set_item("metadata", metadata)?;
    
    Ok(py_dict.into())
}

/// Batch Monte Carlo evaluation with multiple parameter sets
/// 
/// This function allows evaluating multiple Monte Carlo scenarios in parallel,
/// useful for sensitivity analysis across different parameter combinations.
#[pyfunction]
#[pyo3(signature = (base_inputs_list, param_samples_list, param_names_list, _runs_per_scenario, num_threads=None))]
pub fn monte_carlo_batch(
    py: Python,
    base_inputs_list: Vec<Bound<'_, PyDict>>,
    param_samples_list: Vec<PyReadonlyArray2<f64>>,
    param_names_list: Vec<Vec<String>>,
    _runs_per_scenario: Vec<usize>,
    num_threads: Option<usize>,
) -> PyResult<Vec<PyObject>> {
    // Configure thread pool with proper error handling
    configure_thread_pool(num_threads)?;
    
    // Since we can't use into_par_iter on borrowed data, process sequentially
    // (still uses parallelism within each monte_carlo_parallel call)
    let mut results = Vec::new();
    
    for i in 0..base_inputs_list.len() {
        let result = monte_carlo_parallel(
            py, 
            &base_inputs_list[i], 
            param_samples_list[i].clone(), 
            param_names_list[i].clone(), 
            None
        ).unwrap_or_else(|_| py.None());
        results.push(result);
    }
    
    Ok(results)
}

/// Extract ballistic inputs from Python dictionary
fn extract_ballistic_inputs_from_dict(dict: &Bound<'_, PyDict>) -> PyResult<BallisticInputs> {
    // Import the extraction function from lib.rs
    use crate::extract_ballistic_inputs;
    extract_ballistic_inputs(dict)
}

/// Apply a parameter value to the inputs structure
fn apply_parameter(inputs: &mut BallisticInputs, param_name: &str, value: f64) {
    match param_name {
        "bc_value" => inputs.bc_value = value,
        "bullet_mass" => inputs.bullet_mass = value,
        "muzzle_velocity" => inputs.muzzle_velocity = value,
        "wind_speed" => inputs.wind_speed = value,
        "wind_angle" => inputs.wind_angle = value,
        "target_distance" => inputs.target_distance = value,
        "muzzle_angle" => inputs.muzzle_angle = value,
        "altitude" => inputs.altitude = value,
        "temperature" => inputs.temperature = value,
        "pressure" => inputs.pressure = value,
        "humidity" => inputs.humidity = value,
        "latitude" => inputs.latitude = Some(value),
        "twist_rate" => inputs.twist_rate = value,
        "bullet_length" => inputs.bullet_length = value,
        "bullet_diameter" => inputs.bullet_diameter = value,
        _ => {} // Ignore unknown parameters
    }
}

/// Calculate statistics from trajectory results using efficient algorithms
fn calculate_statistics(results: &[Option<TrajectoryOutput>]) -> (HashMap<String, FieldStatistics>, usize, usize) {
    let mut statistics = HashMap::new();
    let valid_results: Vec<&TrajectoryOutput> = results.iter().filter_map(|r| r.as_ref()).collect();
    let valid_runs = valid_results.len();
    let failed_runs = results.len() - valid_runs;
    
    if valid_runs == 0 {
        // Return empty statistics with clear indication of failure
        return (statistics, valid_runs, failed_runs);
    }
    
    // Check for minimum viable sample size
    if valid_runs < 2 {
        // For single samples, we can still provide basic statistics but no variance
        if valid_runs == 1 {
            let sample = &valid_results[0];
            let field_names = vec![
                ("drop_m", sample.drop),
                ("wind_drift_m", sample.wind_drift),
                ("time_of_flight_s", sample.time),
                ("final_vel_fps", sample.velocity * 3.28084),
                ("energy_ft_lbs", sample.energy * 0.737562),
                ("mach", sample.mach),
            ];
            
            for (field_name, value) in field_names {
                statistics.insert(
                    field_name.to_string(),
                    FieldStatistics {
                        mean: value,
                        std: 0.0,
                        min: value,
                        max: value,
                        percentiles: vec![(0.50, value)], // Only median makes sense
                    }
                );
            }
        }
        return (statistics, valid_runs, failed_runs);
    }
    
    // Define field names for statistics
    let field_names = vec![
        "drop_m",
        "wind_drift_m", 
        "time_of_flight_s",
        "final_vel_fps",
        "energy_ft_lbs",
        "mach",
    ];
    
    // Calculate statistics for each field
    for field_name in field_names {
        let mut values: Vec<f64> = valid_results.iter().map(|r| {
            match field_name {
                "drop_m" => r.drop,
                "wind_drift_m" => r.wind_drift,
                "time_of_flight_s" => r.time,
                "final_vel_fps" => r.velocity * 3.28084,
                "energy_ft_lbs" => r.energy * 0.737562,
                "mach" => r.mach,
                _ => 0.0,
            }
        }).collect();
        
        if values.is_empty() {
            continue;
        }
        
        // Sort for percentiles
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        
        // Calculate basic statistics
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let variance = if n > 1.0 {
            values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)
        } else {
            0.0 // Single sample has no variance
        };
        let std = variance.sqrt();
        let min = *values.first().unwrap();
        let max = *values.last().unwrap();
        
        // Calculate percentiles
        let percentiles = vec![
            (0.05, percentile(&values, 0.05)),
            (0.10, percentile(&values, 0.10)),
            (0.25, percentile(&values, 0.25)),
            (0.50, percentile(&values, 0.50)),
            (0.75, percentile(&values, 0.75)),
            (0.90, percentile(&values, 0.90)),
            (0.95, percentile(&values, 0.95)),
        ];
        
        statistics.insert(
            field_name.to_string(),
            FieldStatistics {
                mean,
                std,
                min,
                max,
                percentiles,
            },
        );
    }
    
    (statistics, valid_runs, failed_runs)
}

/// Calculate percentile using linear interpolation
fn percentile(sorted_values: &[f64], p: f64) -> f64 {
    if sorted_values.is_empty() {
        return 0.0;
    }
    
    let n = sorted_values.len();
    if n == 1 {
        return sorted_values[0];
    }
    
    let index = p * (n - 1) as f64;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;
    
    if lower == upper {
        sorted_values[lower]
    } else {
        let weight = index - lower as f64;
        sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight
    }
}

/// Memory-efficient streaming statistics calculator using Welford's algorithm
/// This is useful for very large Monte Carlo runs where storing all results
/// would consume too much memory
#[pyfunction]
#[pyo3(signature = (base_inputs, param_samples, param_names, chunk_size=1000, num_threads=None))]
pub fn monte_carlo_streaming(
    py: Python,
    base_inputs: &Bound<'_, PyDict>,
    param_samples: PyReadonlyArray2<f64>,
    param_names: Vec<String>,
    chunk_size: usize,
    num_threads: Option<usize>,
) -> PyResult<PyObject> {
    // Configure thread pool with proper error handling
    configure_thread_pool(num_threads)?;
    
    // Extract base inputs
    let base_inputs = extract_ballistic_inputs_from_dict(base_inputs)?;
    
    let samples_array = param_samples.as_array();
    let n_samples = samples_array.shape()[0];
    
    // Initialize streaming statistics
    let mut streaming_stats = StreamingStats::new();
    let mut valid_runs = 0;
    let mut failed_runs = 0;
    
    // Process in chunks to limit memory usage
    for chunk_start in (0..n_samples).step_by(chunk_size) {
        let chunk_end = (chunk_start + chunk_size).min(n_samples);
        
        let chunk_results: Vec<Option<TrajectoryOutput>> = (chunk_start..chunk_end)
            .into_par_iter()
            .map(|i| {
                let mut run_inputs = base_inputs.clone();
                
                for (j, param_name) in param_names.iter().enumerate() {
                    let value = samples_array[[i, j]];
                    apply_parameter(&mut run_inputs, param_name, value);
                }
                
                solve_trajectory_for_monte_carlo(&run_inputs).ok()
            })
            .collect();
        
        // Update streaming statistics
        for result in chunk_results {
            if let Some(output) = result {
                streaming_stats.update(&output);
                valid_runs += 1;
            } else {
                failed_runs += 1;
            }
        }
    }
    
    // Extract final statistics
    let statistics = streaming_stats.finalize();
    
    // Convert to Python dictionary
    let py_dict = PyDict::new_bound(py);
    
    let stats_dict = PyDict::new_bound(py);
    for (field_name, (mean, std, min, max)) in statistics {
        let field_dict = PyDict::new_bound(py);
        field_dict.set_item("mean", mean)?;
        field_dict.set_item("std", std)?;
        field_dict.set_item("min", min)?;
        field_dict.set_item("max", max)?;
        stats_dict.set_item(field_name, field_dict)?;
    }
    py_dict.set_item("statistics", stats_dict)?;
    
    let metadata = PyDict::new_bound(py);
    metadata.set_item("valid_runs", valid_runs)?;
    metadata.set_item("failed_runs", failed_runs)?;
    metadata.set_item("total_runs", n_samples)?;
    metadata.set_item("success_rate", valid_runs as f64 / n_samples as f64)?;
    py_dict.set_item("metadata", metadata)?;
    
    Ok(py_dict.into())
}

/// Streaming statistics using Welford's online algorithm
struct StreamingStats {
    n: usize,
    means: HashMap<String, f64>,
    m2s: HashMap<String, f64>,  // Sum of squared differences
    mins: HashMap<String, f64>,
    maxs: HashMap<String, f64>,
}

impl StreamingStats {
    fn new() -> Self {
        let fields = vec![
            "drop_m", "wind_drift_m", "time_of_flight_s",
            "final_vel_fps", "energy_ft_lbs", "mach"
        ];
        
        let mut means = HashMap::new();
        let mut m2s = HashMap::new();
        let mut mins = HashMap::new();
        let mut maxs = HashMap::new();
        
        for field in fields {
            means.insert(field.to_string(), 0.0);
            m2s.insert(field.to_string(), 0.0);
            mins.insert(field.to_string(), f64::INFINITY);
            maxs.insert(field.to_string(), f64::NEG_INFINITY);
        }
        
        Self { n: 0, means, m2s, mins, maxs }
    }
    
    fn update(&mut self, output: &TrajectoryOutput) {
        self.n += 1;
        let n = self.n as f64;
        
        let values = vec![
            ("drop_m", output.drop),
            ("wind_drift_m", output.wind_drift),
            ("time_of_flight_s", output.time),
            ("final_vel_fps", output.velocity * 3.28084),
            ("energy_ft_lbs", output.energy * 0.737562),
            ("mach", output.mach),
        ];
        
        for (field, value) in values {
            let field_str = field.to_string();
            
            // Update min/max
            if let Some(min) = self.mins.get_mut(&field_str) {
                *min = min.min(value);
            }
            if let Some(max) = self.maxs.get_mut(&field_str) {
                *max = max.max(value);
            }
            
            // Welford's algorithm for mean and variance
            if let (Some(mean), Some(m2)) = (self.means.get_mut(&field_str), self.m2s.get_mut(&field_str)) {
                let delta = value - *mean;
                *mean += delta / n;
                let delta2 = value - *mean;
                *m2 += delta * delta2;
            }
        }
    }
    
    fn finalize(&self) -> HashMap<String, (f64, f64, f64, f64)> {
        let mut results = HashMap::new();
        
        for (field, mean) in &self.means {
            let variance = if self.n > 1 {
                self.m2s[field] / (self.n - 1) as f64
            } else {
                0.0
            };
            let std = variance.sqrt();
            
            results.insert(
                field.clone(),
                (*mean, std, self.mins[field], self.maxs[field])
            );
        }
        
        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_percentile_calculation() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(percentile(&values, 0.0), 1.0);
        assert_eq!(percentile(&values, 0.5), 3.0);
        assert_eq!(percentile(&values, 1.0), 5.0);
        assert_eq!(percentile(&values, 0.25), 2.0);
        assert_eq!(percentile(&values, 0.75), 4.0);
    }
    
    #[test]
    fn test_streaming_stats() {
        let mut stats = StreamingStats::new();
        
        // Simulate some outputs
        for i in 1..=100 {
            let output = TrajectoryOutput {
                drop: i as f64,
                wind_drift: i as f64 * 0.5,
                time: i as f64 * 0.1,
                velocity: 100.0,
                energy: 1000.0,
                mach: 1.0,
                spin_drift: 0.0,
                distance: i as f64 * 100.0,
            };
            stats.update(&output);
        }
        
        let results = stats.finalize();
        
        // Check mean of drop (should be 50.5)
        let (mean, _, min, max) = results["drop_m"];
        assert!((mean - 50.5).abs() < 0.01);
        assert_eq!(min, 1.0);
        assert_eq!(max, 100.0);
    }
}