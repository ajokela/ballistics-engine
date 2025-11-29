//! Advanced trajectory integration methods (RK4, RK45)
//!
//! This module provides production-grade numerical integration for ballistic trajectories:
//! - RK4: 4th-order Runge-Kutta (fixed step)
//! - RK45: Dormand-Prince adaptive method (same as scipy.integrate.solve_ivp)
//!
//! MBA-155: Upstreamed from ballistics_rust for shared use

use nalgebra::{Vector3, Vector6};
use std::collections::HashMap;

use crate::derivatives::compute_derivatives;
use crate::wind::WindSegment;
use crate::DragModel;
use crate::BallisticInputs;

/// RK4 integration step
fn rk4_step(
    state: &Vector6<f64>,
    t: f64,
    dt: f64,
    params: &TrajectoryParams,
) -> Vector6<f64> {
    // RK4 integration
    let k1 = compute_derivatives_vec(state, t, params);
    let k2 = compute_derivatives_vec(&(state + dt * 0.5 * k1), t + dt * 0.5, params);
    let k3 = compute_derivatives_vec(&(state + dt * 0.5 * k2), t + dt * 0.5, params);
    let k4 = compute_derivatives_vec(&(state + dt * k3), t + dt, params);

    state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

/// Adaptive RK45 integration step (Dormand-Prince method)
fn rk45_step(
    state: &Vector6<f64>,
    t: f64,
    dt: f64,
    params: &TrajectoryParams,
    tol: f64,
) -> (Vector6<f64>, f64, f64) {
    // Dormand-Prince coefficients (same as scipy.integrate.solve_ivp RK45)
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

    // 4th order coefficients (for error estimation)
    const B1_ERR: f64 = 5179.0 / 57600.0;
    const B3_ERR: f64 = 7571.0 / 16695.0;
    const B4_ERR: f64 = 393.0 / 640.0;
    const B5_ERR: f64 = -92097.0 / 339200.0;
    const B6_ERR: f64 = 187.0 / 2100.0;
    const B7_ERR: f64 = 1.0 / 40.0;

    // Compute stages
    let k1 = compute_derivatives_vec(state, t, params);
    let k2 = compute_derivatives_vec(&(state + dt * A21 * k1), t + dt * 0.2, params);
    let k3 = compute_derivatives_vec(&(state + dt * (A31 * k1 + A32 * k2)), t + dt * 0.3, params);
    let k4 = compute_derivatives_vec(&(state + dt * (A41 * k1 + A42 * k2 + A43 * k3)), t + dt * 0.8, params);
    let k5 = compute_derivatives_vec(&(state + dt * (A51 * k1 + A52 * k2 + A53 * k3 + A54 * k4)), t + dt * 8.0/9.0, params);
    let k6 = compute_derivatives_vec(&(state + dt * (A61 * k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5)), t + dt, params);
    let k7 = compute_derivatives_vec(&(state + dt * (A71 * k1 + A73 * k3 + A74 * k4 + A75 * k5 + A76 * k6)), t + dt, params);

    // 5th order solution
    let y_new = state + dt * (B1 * k1 + B3 * k3 + B4 * k4 + B5 * k5 + B6 * k6);

    // 4th order solution for error estimate
    let y_err = state + dt * (B1_ERR * k1 + B3_ERR * k3 + B4_ERR * k4 + B5_ERR * k5 + B6_ERR * k6 + B7_ERR * k7);

    // Error estimate
    let error = (y_new - y_err).norm() / (1.0 + state.norm());

    // Adaptive step size
    let safety = 0.9;
    let dt_new = if error < tol {
        dt * safety * (tol / error).powf(0.2).min(2.0)
    } else {
        dt * safety * (tol / error).powf(0.25).max(0.1)
    };

    (y_new, dt_new, error)
}

/// Parameters for trajectory computation
pub struct TrajectoryParams {
    pub mass_kg: f64,
    pub bc: f64,
    pub drag_model: DragModel,
    pub wind_segments: Vec<WindSegment>,
    pub atmos_params: (f64, f64, f64, f64),
    pub omega_vector: Option<Vector3<f64>>,
    pub enable_spin_drift: bool,
    pub enable_magnus: bool,
    pub enable_coriolis: bool,
    pub target_distance_m: f64,  // Target horizontal distance in meters
    pub enable_wind_shear: bool,
    pub wind_shear_model: String,
    pub shooter_altitude_m: f64,
    pub is_twist_right: bool,  // True for right-hand twist, false for left-hand
    pub custom_drag_table: Option<crate::drag::DragTable>,  // Custom Drag Model (CDM) data
}

/// Convert state to Vector6 and call compute_derivatives
fn compute_derivatives_vec(
    state: &Vector6<f64>,
    t: f64,
    params: &TrajectoryParams,
) -> Vector6<f64> {

    let pos = Vector3::new(state[0], state[1], state[2]);
    let vel = Vector3::new(state[3], state[4], state[5]);

    // Calculate wind at current position with shear support
    let wind_vector = if !params.wind_segments.is_empty() {
        if params.enable_wind_shear && params.wind_shear_model != "none" {
            crate::wind_shear::get_wind_at_position(
                &pos,
                &params.wind_segments,
                params.enable_wind_shear,
                &params.wind_shear_model,
                params.shooter_altitude_m,
            )
        } else {
            // Simple constant wind (original implementation)
            let seg = &params.wind_segments[0];
            let wind_speed_mps = seg.0 * 0.2777778; // km/h to m/s
            let wind_angle_rad = seg.1.to_radians();
            // Z IS DOWNRANGE: x=lateral, y=vertical, z=downrange
            Vector3::new(
                -wind_speed_mps * wind_angle_rad.sin(),  // x (lateral - crosswind component)
                0.0,                                      // y (vertical)
                -wind_speed_mps * wind_angle_rad.cos(),  // z (downrange - head/tail component)
            )
        }
    } else {
        Vector3::zeros()
    };

    // Create a minimal BallisticInputs struct for the derivatives function
    let inputs = BallisticInputs {
        bc_value: params.bc,
        bc_type: params.drag_model.clone(),
        bullet_mass: params.mass_kg / 0.00006479891, // kg to grains
        muzzle_velocity: vel.norm() * 3.28084, // m/s to fps
        bullet_diameter: 0.308, // default
        bullet_length: 1.24, // default
        twist_rate: 10.0, // default
        is_twist_right: params.is_twist_right,
        enable_advanced_effects: params.enable_spin_drift || params.enable_magnus || params.enable_coriolis,
        altitude: params.atmos_params.0,
        temperature: params.atmos_params.1,
        pressure: params.atmos_params.2,
        humidity: params.atmos_params.3,
        tipoff_yaw: 0.0,
        target_distance: 1000.0, // default
        muzzle_angle: 0.0,
        wind_speed: if !params.wind_segments.is_empty() { params.wind_segments[0].0 } else { 0.0 },
        wind_angle: if !params.wind_segments.is_empty() { params.wind_segments[0].1 } else { 0.0 },
        latitude: None,
        shooting_angle: 0.0,
        azimuth_angle: 0.0,
        use_powder_sensitivity: false,
        powder_temp_sensitivity: 0.0,
        powder_temp: 59.0,
        tipoff_decay_distance: 0.0,
        ground_threshold: -1000.0,
        bc_segments: None,
        caliber_inches: 0.308,
        weight_grains: params.mass_kg / 0.00006479891,
        use_bc_segments: false,
        bullet_id: None,
        bc_segments_data: None,
        use_enhanced_spin_drift: params.enable_spin_drift,
        use_form_factor: false,
        manufacturer: None,
        bullet_model: None,
        enable_wind_shear: false,
        wind_shear_model: "none".to_string(),
        use_cluster_bc: false,
        bullet_cluster: None,

        // Pass through custom drag table (CDM) from trajectory parameters
        custom_drag_table: params.custom_drag_table.clone(),

        bc_type_str: None,
        enable_pitch_damping: false,
        enable_precession_nutation: false,
        use_rk4: true,
        use_adaptive_rk45: false,
        enable_trajectory_sampling: false,
        sample_interval: 10.0,
        sight_height: 0.0,
        muzzle_height: 0.0,
        target_height: 0.0,
    };

    // Call compute_derivatives - returns [f64; 6] directly
    let deriv_result = compute_derivatives(
        pos,
        vel,
        &inputs,
        wind_vector,
        params.atmos_params,
        params.bc,
        params.omega_vector,
        t,
    );

    Vector6::new(
        deriv_result[0], deriv_result[1], deriv_result[2],
        deriv_result[3], deriv_result[4], deriv_result[5],
    )
}

/// Main trajectory integration function
pub fn integrate_trajectory(
    initial_state: [f64; 6],
    t_span: (f64, f64),
    params: TrajectoryParams,
    method: &str,
    tolerance: f64,
    max_step: f64,
) -> Vec<(f64, Vector6<f64>)> {
    let mut state = Vector6::new(
        initial_state[0], initial_state[1], initial_state[2],
        initial_state[3], initial_state[4], initial_state[5],
    );

    let mut t = t_span.0;
    let t_end = t_span.1;
    let mut dt = (t_end - t) / 1000.0; // Initial step size

    let mut trajectory = Vec::with_capacity(10000);
    trajectory.push((t, state.clone()));

    match method {
        "RK4" => {
            // Fixed step RK4 with target detection
            dt = dt.min(max_step).min(0.001); // Use smaller steps for accuracy

            while t < t_end {
                if t + dt > t_end {
                    dt = t_end - t;
                }

                let new_state = rk4_step(&state, t, dt, &params);

                // Check if we're about to pass the target (z is downrange)
                if state[2] < params.target_distance_m && new_state[2] >= params.target_distance_m {
                    // Interpolate to find exact target crossing
                    let alpha = (params.target_distance_m - state[2]) / (new_state[2] - state[2]);
                    let dt_to_target = dt * alpha;

                    // Take a smaller step to reach target exactly
                    let final_state = rk4_step(&state, t, dt_to_target, &params);

                    // Ensure we don't overshoot
                    let mut corrected_state = final_state;
                    if corrected_state[2] > params.target_distance_m {
                        corrected_state[2] = params.target_distance_m;
                    }

                    trajectory.push((t + dt_to_target, corrected_state));
                    break;  // Stop at target
                }

                state = new_state;
                t += dt;
                trajectory.push((t, state.clone()));

                // Check if we've reached or passed the target
                if state[2] >= params.target_distance_m {  // z is downrange
                    // Add final point exactly at target
                    let mut final_state = state;
                    final_state[2] = params.target_distance_m;  // z is downrange
                    trajectory.push((t, final_state));
                    break;
                }

                // Check if bullet hit ground
                if state[1] < -1000.0 {
                    break;
                }
            }
        }
        "RK45" | _ => {
            // Adaptive RK45 with better sampling
            let mut last_save_z = 0.0;  // z is downrange
            let save_interval_m = params.target_distance_m / 50.0; // Save ~50 points minimum

            // OPTIMIZATION: Adjust max step size when wind shear is enabled
            // This improves numerical stability at long ranges
            let effective_max_step = if params.enable_wind_shear && params.wind_shear_model != "none" {
                // Use smaller steps for wind shear, but not TOO small
                if params.target_distance_m > 800.0 {
                    0.01  // Smaller steps for long range with shear (10ms)
                } else {
                    0.02  // Normal steps for medium range with shear (20ms)
                }
            } else {
                max_step  // Use provided max_step when no wind shear
            };

            // Set initial step size - ensure it's reasonable
            dt = dt.min(effective_max_step).max(0.0001);  // At least 0.1ms to avoid infinite loops

            // Safety check: maximum iterations to prevent infinite loops
            let max_iterations = 100000;  // Should be more than enough for any realistic trajectory
            let mut iteration_count = 0;

            while t < t_end && iteration_count < max_iterations {
                iteration_count += 1;

                // Limit time step for better resolution
                if t + dt > t_end {
                    dt = t_end - t;
                }

                let (new_state, dt_new, _error) = rk45_step(&state, t, dt, &params, tolerance);

                // Check if we're about to pass the target (z is downrange)
                if state[2] < params.target_distance_m && new_state[2] >= params.target_distance_m {
                    // Interpolate to find exact target crossing
                    let alpha = (params.target_distance_m - state[2]) / (new_state[2] - state[2]);
                    let dt_to_target = dt * alpha;

                    // Take a smaller step to reach target exactly
                    let (final_state, _, _) = rk45_step(&state, t, dt_to_target, &params, tolerance);

                    // Make sure we don't overshoot
                    let mut corrected_state = final_state;
                    if corrected_state[2] > params.target_distance_m {
                        corrected_state[2] = params.target_distance_m;
                    }

                    trajectory.push((t + dt_to_target, corrected_state));
                    break;  // Stop at target - no more points after this
                }

                // Update state
                state = new_state;
                t += dt;

                // Save trajectory point if we've moved enough distance
                if state[2] - last_save_z >= save_interval_m || state[2] >= params.target_distance_m {  // z is downrange
                    trajectory.push((t, state.clone()));
                    last_save_z = state[2];
                }

                // Limit dt for next step - ensure we get enough resolution
                dt = dt_new.min(effective_max_step).max(0.0001); // Use effective max step, min 0.1ms

                // Stop if we've reached the target
                if state[2] >= params.target_distance_m {  // z is downrange
                    // Add final point at target distance
                    let mut final_state = state;
                    final_state[2] = params.target_distance_m;  // z is downrange
                    trajectory.push((t, final_state));
                    break;
                }

                // Check if bullet hit ground
                if state[1] < -1000.0 {
                    break;
                }
            }

            // Warn if we hit the iteration limit
            if iteration_count >= max_iterations {
                eprintln!("WARNING: Trajectory integration hit maximum iteration limit ({} iterations)", max_iterations);
                eprintln!("  Final time: {}, Target time: {}", t, t_end);
                eprintln!("  Final position: z={}, Target: {}m", state[2], params.target_distance_m);
            }
        }
    }

    trajectory
}

/// Python-exposed function for complete trajectory integration
pub fn solve_trajectory_rust(
    initial_state: [f64; 6],
    t_span: (f64, f64),
    mass_kg: f64,
    bc: f64,
    drag_model: DragModel,
    wind_segments: Vec<WindSegment>,
    atmos_params: (f64, f64, f64, f64),
    omega_vector: Option<Vec<f64>>,
    enable_spin_drift: bool,
    enable_magnus: bool,
    enable_coriolis: bool,
    method: String,
    tolerance: f64,
    max_step: f64,
    target_distance_m: f64,
) -> Vec<HashMap<String, f64>> {
    let omega_vec = omega_vector.map(|v| Vector3::new(v[0], v[1], v[2]));

    let params = TrajectoryParams {
        mass_kg,
        bc,
        drag_model,
        wind_segments,
        atmos_params,
        omega_vector: omega_vec,
        enable_spin_drift,
        enable_magnus,
        enable_coriolis,
        target_distance_m,
        enable_wind_shear: false,  // Default for test function
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: 0.0,
        is_twist_right: true,  // Default for test function
        custom_drag_table: None,  // No CDM for test function
    };

    let trajectory = integrate_trajectory(
        initial_state,
        t_span,
        params,
        &method,
        tolerance,
        max_step,
    );

    // Convert to Python-friendly format
    trajectory.into_iter().map(|(t, state)| {
        let mut point = HashMap::new();
        point.insert("t".to_string(), t);
        point.insert("x".to_string(), state[0]);
        point.insert("y".to_string(), state[1]);
        point.insert("z".to_string(), state[2]);
        point.insert("vx".to_string(), state[3]);
        point.insert("vy".to_string(), state[4]);
        point.insert("vz".to_string(), state[5]);
        point
    }).collect()
}
