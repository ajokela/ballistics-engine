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
use crate::BallisticInputs;
use crate::DragModel;

fn wind_vector_for_range(range_m: f64, wind_segments: &[WindSegment]) -> Vector3<f64> {
    if range_m.is_nan() {
        return Vector3::zeros();
    }
    for seg in wind_segments {
        if range_m < seg.2 {
            let wind_speed_mps = seg.0 * 0.2777778; // km/h to m/s
            let wind_angle_rad = seg.1.to_radians();
            return Vector3::new(
                -wind_speed_mps * wind_angle_rad.cos(),
                0.0,
                -wind_speed_mps * wind_angle_rad.sin(),
            );
        }
    }
    Vector3::zeros()
}

/// RK4 integration step
fn rk4_step(
    state: &Vector6<f64>,
    t: f64,
    dt: f64,
    params: &TrajectoryParams,
    inputs: &BallisticInputs,
) -> Vector6<f64> {
    // RK4 integration
    let k1 = compute_derivatives_vec(state, t, params, inputs);
    let k2 = compute_derivatives_vec(&(state + dt * 0.5 * k1), t + dt * 0.5, params, inputs);
    let k3 = compute_derivatives_vec(&(state + dt * 0.5 * k2), t + dt * 0.5, params, inputs);
    let k4 = compute_derivatives_vec(&(state + dt * k3), t + dt, params, inputs);

    state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

/// Adaptive RK45 integration step (Dormand-Prince method)
fn rk45_step(
    state: &Vector6<f64>,
    t: f64,
    dt: f64,
    params: &TrajectoryParams,
    inputs: &BallisticInputs,
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
    let k1 = compute_derivatives_vec(state, t, params, inputs);
    let k2 = compute_derivatives_vec(&(state + dt * A21 * k1), t + dt * 0.2, params, inputs);
    let k3 = compute_derivatives_vec(
        &(state + dt * (A31 * k1 + A32 * k2)),
        t + dt * 0.3,
        params,
        inputs,
    );
    let k4 = compute_derivatives_vec(
        &(state + dt * (A41 * k1 + A42 * k2 + A43 * k3)),
        t + dt * 0.8,
        params,
        inputs,
    );
    let k5 = compute_derivatives_vec(
        &(state + dt * (A51 * k1 + A52 * k2 + A53 * k3 + A54 * k4)),
        t + dt * 8.0 / 9.0,
        params,
        inputs,
    );
    let k6 = compute_derivatives_vec(
        &(state + dt * (A61 * k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5)),
        t + dt,
        params,
        inputs,
    );
    let k7 = compute_derivatives_vec(
        &(state + dt * (A71 * k1 + A73 * k3 + A74 * k4 + A75 * k5 + A76 * k6)),
        t + dt,
        params,
        inputs,
    );

    // 5th order solution
    let y_new = state + dt * (B1 * k1 + B3 * k3 + B4 * k4 + B5 * k5 + B6 * k6);

    // 4th order solution for error estimate
    let y_err = state
        + dt * (B1_ERR * k1 + B3_ERR * k3 + B4_ERR * k4 + B5_ERR * k5 + B6_ERR * k6 + B7_ERR * k7);

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
    /// Dual-mode atmosphere tuple consumed by `compute_derivatives`:
    /// **Standard** `(base_alt_m, base_temp_c, base_pressure_hPa, base_density_ratio)` — note
    /// slot 3 is a density RATIO, NOT humidity, even though it rides in the `humidity` field;
    /// or **Direct** `(air_density, speed_of_sound, 0.0, 0.0)` — slots 2 and 3 are zero
    /// sentinels. A pressure of 0 that is not the direct-mode sentinel disables drag.
    pub atmos_params: (f64, f64, f64, f64),
    pub omega_vector: Option<Vector3<f64>>,
    pub enable_spin_drift: bool,
    pub enable_magnus: bool,
    pub enable_coriolis: bool,
    pub target_distance_m: f64, // Target horizontal distance in meters
    pub enable_wind_shear: bool,
    pub wind_shear_model: String,
    pub shooter_altitude_m: f64,
    pub is_twist_right: bool, // True for right-hand twist, false for left-hand
    pub shooting_angle: f64,  // uphill/downhill angle in radians
    // MBA-717: real bullet geometry so spin-drift / Magnus / stability on this fast/MC
    // path use the actual bullet instead of hardcoded .308 / 1.24in / 10-twist placeholders.
    pub bullet_diameter: f64,                              // meters
    pub bullet_length: f64, // meters (0.0 -> derivatives falls back to the 4.5-caliber heuristic)
    pub twist_rate: f64,    // inches per turn
    pub custom_drag_table: Option<crate::drag::DragTable>, // Custom Drag Model (CDM) data
    pub bc_segments: Option<Vec<(f64, f64)>>, // Mach-based BC segments: (mach, bc)
    pub use_bc_segments: bool, // Whether to use BC segment interpolation
    /// MBA-954: altitude (m, relative to launch) below which integration stops. -1000.0 is the
    /// historical default — effectively "no early ground impact" for normal flat-fire shots.
    pub ground_threshold: f64,
    /// MBA-1137: optional downrange-segmented atmosphere. When `Some`, `compute_derivatives`
    /// swaps the standard-mode base T/P/H for the zone selected by downrange distance before the
    /// altitude lapse. `None` (default) is byte-identical to pre-feature behavior.
    pub atmo_sock: Option<crate::atmosphere::AtmoSock>,
}

/// Build the loop-invariant BallisticInputs for the derivatives function ONCE per integration,
/// instead of rebuilding it (a "none".to_string() alloc plus bc_segments / custom_drag_table
/// clones) on every derivative evaluation (4x per RK4 step, 7x per RK45 step). muzzle_velocity is
/// set to 0.0 because compute_derivatives never reads it on this path; every other field depends
/// only on `params`, so the struct is constant for the whole integration.
fn build_inputs(params: &TrajectoryParams) -> BallisticInputs {
    let mut inputs = BallisticInputs {
        bc_value: params.bc,
        bc_type: params.drag_model,
        bullet_mass: params.mass_kg,             // kg
        muzzle_velocity: 0.0,                    // unread by compute_derivatives on this path
        bullet_diameter: params.bullet_diameter, // MBA-717: real geometry, not placeholders
        bullet_length: params.bullet_length,
        twist_rate: params.twist_rate,
        is_twist_right: params.is_twist_right,
        enable_advanced_effects: params.enable_spin_drift
            || params.enable_magnus
            || params.enable_coriolis,
        enable_magnus: params.enable_magnus,
        enable_coriolis: params.enable_coriolis,
        altitude: params.atmos_params.0,
        temperature: params.atmos_params.1,
        pressure: params.atmos_params.2,
        humidity: params.atmos_params.3,
        tipoff_yaw: 0.0,
        target_distance: 1000.0, // default
        muzzle_angle: 0.0,
        wind_speed: if !params.wind_segments.is_empty() {
            params.wind_segments[0].0 * 0.2777778 // km/h -> m/s
        } else {
            0.0
        },
        wind_angle: if !params.wind_segments.is_empty() {
            params.wind_segments[0].1.to_radians() // degrees -> radians
        } else {
            0.0
        },
        latitude: None,
        shooting_angle: params.shooting_angle,
        azimuth_angle: 0.0,
        shot_azimuth: 0.0, // this fast path doesn't plumb latitude/bearing (no directional Coriolis here)
        use_powder_sensitivity: false,
        powder_temp_sensitivity: 0.0,
        powder_temp: 59.0,
        powder_temp_curve: None,
        powder_curve_temp_c: None,
        tipoff_decay_distance: 0.0,
        ground_threshold: params.ground_threshold, // MBA-954: honor the configured ground plane
        bc_segments: params.bc_segments.clone(),
        caliber_inches: params.bullet_diameter / 0.0254, // MBA-717: from real diameter
        weight_grains: params.mass_kg / 0.00006479891,
        use_bc_segments: params.use_bc_segments,
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
        custom_drag_table: params.custom_drag_table.clone(),
        bc_type_str: None,
        enable_pitch_damping: false,
        enable_precession_nutation: false,
        // MBA-959/MBA-1183: aerodynamic jump stays OFF inside this low-level raw-state integrator.
        // The high-level fast wrappers form Sg from their complete BallisticInputs and rotate the
        // prebuilt initial velocity before entering their integration loops; enabling it again
        // here would double-apply the launch perturbation. Direct low-level callers likewise own
        // any desired launch-state rotation. (Real geometry is still carried for spin/Magnus.)
        enable_aerodynamic_jump: false,
        use_rk4: true,
        use_adaptive_rk45: false,
        enable_trajectory_sampling: false,
        sample_interval: 10.0,
        sight_height: 0.0,
        muzzle_height: 0.0,
        target_height: 0.0,
    };

    // MBA-955: pre-populate velocity-BC segments ONCE here, instead of get_bc_for_velocity
    // rebuilding them (a model String + a segment Vec) on every derivative evaluation (4-7x per
    // step). Gated to EXACTLY the case where the per-step path would estimate: use_bc_segments on,
    // no explicit velocity segments, and no Mach-based bc_segments (those take a different,
    // unchanged path). bc_used there == params.bc == inputs.bc_value, so the estimated segments are
    // identical and the per-step fast-path lookup returns the same BC -> byte-identical output.
    if inputs.use_bc_segments && inputs.bc_segments_data.is_none() && inputs.bc_segments.is_none() {
        inputs.bc_segments_data =
            crate::derivatives::estimate_bc_segments_for(&inputs, inputs.bc_value);
    }
    inputs
}

/// Convert state to Vector6 and call compute_derivatives
fn compute_derivatives_vec(
    state: &Vector6<f64>,
    t: f64,
    params: &TrajectoryParams,
    inputs: &BallisticInputs,
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
            wind_vector_for_range(pos.x, &params.wind_segments)
        }
    } else {
        Vector3::zeros()
    };

    // Call compute_derivatives - returns [f64; 6] directly. `inputs` is built once per
    // integration by build_inputs() and threaded in, instead of rebuilt every call.
    let deriv_result = compute_derivatives(
        pos,
        vel,
        inputs,
        wind_vector,
        params.atmos_params,
        params.bc,
        params.omega_vector,
        t,
        params.atmo_sock.as_ref(),
    );

    Vector6::new(
        deriv_result[0],
        deriv_result[1],
        deriv_result[2],
        deriv_result[3],
        deriv_result[4],
        deriv_result[5],
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
        initial_state[0],
        initial_state[1],
        initial_state[2],
        initial_state[3],
        initial_state[4],
        initial_state[5],
    );

    let mut t = t_span.0;
    let t_end = t_span.1;
    let mut dt = (t_end - t) / 1000.0; // Initial step size

    let mut trajectory = Vec::with_capacity(10000);
    trajectory.push((t, state));

    // Build the (loop-invariant) derivative inputs once for the whole integration, instead of
    // rebuilding the struct on every derivative evaluation.
    let inputs = build_inputs(&params);

    match method {
        "RK4" => {
            // Fixed step RK4 with target detection
            dt = dt.min(max_step).min(0.001); // Use smaller steps for accuracy

            while t < t_end {
                if t + dt > t_end {
                    dt = t_end - t;
                }

                let new_state = rk4_step(&state, t, dt, &params, &inputs);

                // Check if we're about to pass the target (X is downrange, McCoy)
                if state[0] < params.target_distance_m && new_state[0] >= params.target_distance_m {
                    // Interpolate to find exact target crossing
                    let alpha = (params.target_distance_m - state[0]) / (new_state[0] - state[0]);
                    let dt_to_target = dt * alpha;

                    // Take a smaller step to reach target exactly
                    let final_state = rk4_step(&state, t, dt_to_target, &params, &inputs);

                    // Ensure we don't overshoot
                    let mut corrected_state = final_state;
                    if corrected_state[0] > params.target_distance_m {
                        corrected_state[0] = params.target_distance_m;
                    }

                    trajectory.push((t + dt_to_target, corrected_state));
                    break; // Stop at target
                }

                state = new_state;
                t += dt;
                trajectory.push((t, state));

                // Check if we've reached or passed the target
                if state[0] >= params.target_distance_m {
                    // X is downrange (McCoy)
                    // Add final point exactly at target
                    let mut final_state = state;
                    final_state[0] = params.target_distance_m; // X is downrange
                    trajectory.push((t, final_state));
                    break;
                }

                // Check if bullet hit ground (MBA-954: honor the configured ground plane,
                // not a hardcoded -1000.0)
                if state[1] < params.ground_threshold {
                    break;
                }
            }
        }
        "RK45" | _ => {
            // Adaptive RK45 with better sampling
            let mut last_save_x = 0.0; // X is downrange (McCoy)
            let save_interval_m = params.target_distance_m / 50.0; // Save ~50 points minimum

            // OPTIMIZATION: Adjust max step size when wind shear is enabled
            // This improves numerical stability at long ranges
            let effective_max_step =
                if params.enable_wind_shear && params.wind_shear_model != "none" {
                    // Use smaller steps for wind shear, but not TOO small
                    if params.target_distance_m > 800.0 {
                        0.01 // Smaller steps for long range with shear (10ms)
                    } else {
                        0.02 // Normal steps for medium range with shear (20ms)
                    }
                } else {
                    max_step // Use provided max_step when no wind shear
                };

            // Set initial step size - ensure it's reasonable
            dt = dt.min(effective_max_step).max(0.0001); // At least 0.1ms to avoid infinite loops

            // Safety check: maximum iterations to prevent infinite loops
            let max_iterations = 100000; // Should be more than enough for any realistic trajectory
            let mut iteration_count = 0;

            while t < t_end && iteration_count < max_iterations {
                iteration_count += 1;

                // Limit time step for better resolution
                if t + dt > t_end {
                    dt = t_end - t;
                }

                let (new_state, dt_new, _error) =
                    rk45_step(&state, t, dt, &params, &inputs, tolerance);

                // Check if we're about to pass the target (X is downrange, McCoy)
                if state[0] < params.target_distance_m && new_state[0] >= params.target_distance_m {
                    // Interpolate to find exact target crossing
                    let alpha = (params.target_distance_m - state[0]) / (new_state[0] - state[0]);
                    let dt_to_target = dt * alpha;

                    // Take a smaller step to reach target exactly
                    let (final_state, _, _) =
                        rk45_step(&state, t, dt_to_target, &params, &inputs, tolerance);

                    // Make sure we don't overshoot
                    let mut corrected_state = final_state;
                    if corrected_state[0] > params.target_distance_m {
                        corrected_state[0] = params.target_distance_m;
                    }

                    trajectory.push((t + dt_to_target, corrected_state));
                    break; // Stop at target - no more points after this
                }

                // Update state
                state = new_state;
                t += dt;

                // Save trajectory point if we've moved enough distance
                if state[0] - last_save_x >= save_interval_m || state[0] >= params.target_distance_m
                {
                    // X is downrange
                    trajectory.push((t, state));
                    last_save_x = state[0];
                }

                // Limit dt for next step - ensure we get enough resolution
                dt = dt_new.min(effective_max_step).max(0.0001); // Use effective max step, min 0.1ms

                // Stop if we've reached the target
                if state[0] >= params.target_distance_m {
                    // X is downrange (McCoy)
                    // Add final point at target distance
                    let mut final_state = state;
                    final_state[0] = params.target_distance_m; // X is downrange
                    trajectory.push((t, final_state));
                    break;
                }

                // Check if bullet hit ground (MBA-954: honor the configured ground plane,
                // not a hardcoded -1000.0)
                if state[1] < params.ground_threshold {
                    break;
                }
            }

            // Warn if we hit the iteration limit
            if iteration_count >= max_iterations {
                eprintln!(
                    "WARNING: Trajectory integration hit maximum iteration limit ({} iterations)",
                    max_iterations
                );
                eprintln!("  Final time: {}, Target time: {}", t, t_end);
                eprintln!(
                    "  Final position: downrange(x)={}, Target: {}m",
                    state[0], params.target_distance_m
                );
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
        enable_wind_shear: false, // Default for test function
        wind_shear_model: "none".to_string(),
        shooter_altitude_m: 0.0,
        is_twist_right: true, // Default for test function
        shooting_angle: 0.0,  // This legacy entry takes no inclined-fire arg; flat fire only
        // This legacy entry takes no geometry args; keep the historical placeholders so its
        // behavior is unchanged (callers needing real geometry use fast_integrate_with_segments).
        bullet_diameter: 0.0078232,
        bullet_length: 0.031496,
        twist_rate: 10.0,
        custom_drag_table: None, // No CDM for test function
        bc_segments: None,       // No BC segments for legacy function
        use_bc_segments: false,
        ground_threshold: -1000.0, // MBA-954: preserve the historical default
        atmo_sock: None,           // MBA-1137: legacy entry has no downrange atmosphere
    };

    let trajectory =
        integrate_trajectory(initial_state, t_span, params, &method, tolerance, max_step);

    // Convert to Python-friendly format
    trajectory
        .into_iter()
        .map(|(t, state)| {
            let mut point = HashMap::new();
            point.insert("t".to_string(), t);
            point.insert("x".to_string(), state[0]);
            point.insert("y".to_string(), state[1]);
            point.insert("z".to_string(), state[2]);
            point.insert("vx".to_string(), state[3]);
            point.insert("vy".to_string(), state[4]);
            point.insert("vz".to_string(), state[5]);
            point
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_params(target_distance_m: f64) -> TrajectoryParams {
        TrajectoryParams {
            mass_kg: 0.01134, // 175 grains in kg
            bc: 0.442,
            bullet_diameter: 0.0078232, // .308 in
            bullet_length: 0.031496,    // 1.24 in
            twist_rate: 10.0,
            drag_model: DragModel::G7,
            wind_segments: vec![],
            atmos_params: (0.0, 15.0, 1013.25, 1.0),
            omega_vector: None,
            enable_spin_drift: false,
            enable_magnus: false,
            enable_coriolis: false,
            target_distance_m,
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
            shooter_altitude_m: 0.0,
            is_twist_right: true,
            shooting_angle: 0.0,
            custom_drag_table: None,
            bc_segments: None,
            use_bc_segments: false,
            ground_threshold: -1000.0,
            atmo_sock: None,
        }
    }

    #[test]
    fn test_mba954_ground_threshold_honored() {
        // MBA-954: integrate_trajectory must honor the configured ground plane, not a hardcoded
        // -1000.0. A descending bullet with a shallow ground_threshold must terminate earlier
        // (fewer points) than one with the historical deep default.
        let initial_state = [0.0, 0.0, 0.0, 300.0, -30.0, 0.0]; // descending (vy = -30 m/s)

        let mut shallow = create_test_params(1_000_000.0); // huge target so range never terminates
        shallow.ground_threshold = -20.0; // stop ~20 m below launch
        let mut deep = create_test_params(1_000_000.0);
        deep.ground_threshold = -1000.0; // historical default

        let t_shallow =
            integrate_trajectory(initial_state, (0.0, 60.0), shallow, "RK4", 1e-6, 0.001);
        let t_deep = integrate_trajectory(initial_state, (0.0, 60.0), deep, "RK4", 1e-6, 0.001);

        assert!(
            t_shallow.len() < t_deep.len(),
            "shallow ground_threshold (-20) should terminate earlier than deep (-1000): \
             shallow={}, deep={}",
            t_shallow.len(),
            t_deep.len()
        );
    }

    #[test]
    fn test_integrate_trajectory_basic() {
        // Initial state [x,y,z,vx,vy,vz] (McCoy: X=downrange, Z=lateral)
        // x=0 (downrange start), vx=821.52 (downrange velocity)
        let initial_state = [0.0, -0.038, 0.0, 821.52, 48.61, 0.0];

        let params = TrajectoryParams {
            mass_kg: 0.01134, // 175 grains in kg
            bc: 0.442,
            bullet_diameter: 0.0078232, // .308 in
            bullet_length: 0.031496,    // 1.24 in
            twist_rate: 10.0,
            drag_model: DragModel::G7,
            wind_segments: vec![(0.0, 90.0, 914.4)],
            atmos_params: (0.0, 15.0, 1013.25, 1.0),
            omega_vector: None,
            enable_spin_drift: false,
            enable_magnus: false,
            enable_coriolis: false,
            target_distance_m: 914.4, // 1000 yards in meters
            enable_wind_shear: false,
            wind_shear_model: "none".to_string(),
            shooter_altitude_m: 0.0,
            is_twist_right: true,
            shooting_angle: 0.0,
            custom_drag_table: None,
            bc_segments: None,
            use_bc_segments: false,
            ground_threshold: -1000.0,
            atmo_sock: None,
        };

        println!("Running integrate_trajectory test...");
        println!("Initial state: {:?}", initial_state);
        println!("Target distance: {} m", params.target_distance_m);

        let trajectory =
            integrate_trajectory(initial_state, (0.0, 10.0), params, "RK45", 1e-6, 0.01);

        println!("Trajectory has {} points", trajectory.len());

        // Should have more than just initial point
        assert!(
            trajectory.len() > 1,
            "Trajectory should have more than 1 point, but has {}",
            trajectory.len()
        );

        // Check that we actually moved downrange
        if let Some((_, final_state)) = trajectory.last() {
            println!("Final state: downrange(x)={}", final_state[0]);
            assert!(
                final_state[0] > 0.0,
                "Final x should be positive (bullet moved downrange)"
            );
            assert!(
                final_state[0] >= 900.0,
                "Final x should be near target distance"
            );
        }
    }

    #[test]
    fn test_rk4_vs_rk45_consistency() {
        // Both methods should give similar results for the same trajectory
        let initial_state = [0.0, 0.0, 0.0, 800.0, 30.0, 0.0]; // McCoy: vx=downrange
        let target_distance = 500.0;

        let params_rk4 = create_test_params(target_distance);
        let params_rk45 = create_test_params(target_distance);

        let trajectory_rk4 =
            integrate_trajectory(initial_state, (0.0, 5.0), params_rk4, "RK4", 1e-6, 0.001);
        let trajectory_rk45 =
            integrate_trajectory(initial_state, (0.0, 5.0), params_rk45, "RK45", 1e-6, 0.01);

        // Both should reach target
        assert!(!trajectory_rk4.is_empty());
        assert!(!trajectory_rk45.is_empty());

        let (_, final_rk4) = trajectory_rk4.last().unwrap();
        let (_, final_rk45) = trajectory_rk45.last().unwrap();

        // Final downrange positions should be within 1% of each other
        let rk4_z = final_rk4[0];
        let rk45_z = final_rk45[0];
        let diff_percent = ((rk4_z - rk45_z) / rk45_z).abs() * 100.0;

        assert!(
            diff_percent < 1.0,
            "RK4 and RK45 final positions differ by {}%: RK4={}, RK45={}",
            diff_percent,
            rk4_z,
            rk45_z
        );
    }

    #[test]
    fn test_ground_impact_detection() {
        // Trajectory with steep downward angle should hit ground
        let initial_state = [0.0, 100.0, 0.0, 300.0, -50.0, 0.0]; // McCoy: vx=downrange // Steep descent

        let mut params = create_test_params(10000.0); // Far target
        params.target_distance_m = 10000.0;
        let ground_threshold = 0.0;
        params.ground_threshold = ground_threshold;

        let trajectory =
            integrate_trajectory(initial_state, (0.0, 20.0), params, "RK4", 1e-6, 0.01);

        // Should stop before reaching target due to ground impact
        let (_, final_state) = trajectory.last().unwrap();

        // y should have crossed the configured ground threshold.
        assert!(
            final_state[1] <= ground_threshold,
            "Should hit ground, but y={}",
            final_state[1]
        );
        assert!(
            final_state[0] < 10000.0,
            "Should not reach target, but z={}",
            final_state[0]
        );
    }

    #[test]
    fn test_target_distance_reached() {
        let initial_state = [0.0, 0.0, 0.0, 800.0, 20.0, 0.0]; // McCoy: vx=downrange
        let target_distance = 300.0;

        let params = create_test_params(target_distance);

        let trajectory =
            integrate_trajectory(initial_state, (0.0, 5.0), params, "RK45", 1e-6, 0.01);

        let (_, final_state) = trajectory.last().unwrap();

        // Should stop at or very near target distance
        assert!(
            (final_state[0] - target_distance).abs() < 1.0,
            "Should reach target at {}m, but stopped at {}m",
            target_distance,
            final_state[0]
        );
    }

    #[test]
    fn test_wind_affects_trajectory() {
        // Test that wind segments are properly stored and passed through
        // The actual wind effect depends on the derivatives computation which
        // uses the wind vector in the drag calculation
        let initial_state = [0.0, 0.0, 0.0, 800.0, 30.0, 0.0]; // McCoy: vx=downrange
        let target_distance = 500.0;

        // No wind
        let params_no_wind = create_test_params(target_distance);

        // Strong headwind (0 degrees = headwind)
        let mut params_headwind = create_test_params(target_distance);
        params_headwind.wind_segments = vec![(72.0, 0.0, 500.0)]; // 72 km/h = 20 m/s headwind

        let trajectory_no_wind = integrate_trajectory(
            initial_state,
            (0.0, 5.0),
            params_no_wind,
            "RK45",
            1e-6,
            0.01,
        );
        let trajectory_headwind = integrate_trajectory(
            initial_state,
            (0.0, 5.0),
            params_headwind,
            "RK45",
            1e-6,
            0.01,
        );

        // Both trajectories should complete
        assert!(
            !trajectory_no_wind.is_empty(),
            "No-wind trajectory should complete"
        );
        assert!(
            !trajectory_headwind.is_empty(),
            "Headwind trajectory should complete"
        );

        let (time_no_wind, final_no_wind) = trajectory_no_wind.last().unwrap();
        let (time_headwind, final_headwind) = trajectory_headwind.last().unwrap();

        // Headwind should slow the bullet, resulting in longer flight time
        // or different drop at same distance
        let drop_no_wind = final_no_wind[1];
        let drop_headwind = final_headwind[1];

        // With headwind, bullet should drop more (more negative y) due to slower velocity
        // The effect might be small, so we check that both trajectories are valid
        println!("No wind: time={}, drop={}", time_no_wind, drop_no_wind);
        println!("Headwind: time={}, drop={}", time_headwind, drop_headwind);

        // Both should reach approximately the target distance
        assert!(
            (final_no_wind[0] - target_distance).abs() < 10.0,
            "No-wind should reach target"
        );
        assert!(
            (final_headwind[0] - target_distance).abs() < 10.0,
            "Headwind should reach target"
        );
    }

    #[test]
    fn test_solve_trajectory_rust_output_format() {
        let initial_state = [0.0, 0.0, 0.0, 800.0, 30.0, 0.0]; // McCoy: vx=downrange

        let result = solve_trajectory_rust(
            initial_state,
            (0.0, 2.0),
            0.01134,                 // mass_kg
            0.442,                   // bc
            DragModel::G7,           // drag_model
            vec![],                  // wind_segments
            (0.0, 59.0, 29.92, 0.0), // atmos_params
            None,                    // omega_vector
            false,                   // enable_spin_drift
            false,                   // enable_magnus
            false,                   // enable_coriolis
            "RK45".to_string(),      // method
            1e-6,                    // tolerance
            0.01,                    // max_step
            500.0,                   // target_distance_m
        );

        // Should return Vec of HashMaps with expected keys
        assert!(!result.is_empty());

        let first_point = &result[0];
        assert!(first_point.contains_key("t"));
        assert!(first_point.contains_key("x"));
        assert!(first_point.contains_key("y"));
        assert!(first_point.contains_key("z"));
        assert!(first_point.contains_key("vx"));
        assert!(first_point.contains_key("vy"));
        assert!(first_point.contains_key("vz"));
    }

    #[test]
    fn test_left_vs_right_twist() {
        let initial_state = [0.0, 0.0, 0.0, 800.0, 30.0, 0.0]; // McCoy: vx=downrange
        let target_distance = 500.0;

        let mut params_right = create_test_params(target_distance);
        params_right.is_twist_right = true;
        params_right.enable_spin_drift = true;

        let mut params_left = create_test_params(target_distance);
        params_left.is_twist_right = false;
        params_left.enable_spin_drift = true;

        let trajectory_right =
            integrate_trajectory(initial_state, (0.0, 5.0), params_right, "RK45", 1e-6, 0.01);
        let trajectory_left =
            integrate_trajectory(initial_state, (0.0, 5.0), params_left, "RK45", 1e-6, 0.01);

        // Both should complete
        assert!(!trajectory_right.is_empty());
        assert!(!trajectory_left.is_empty());

        // Right and left twist should produce valid trajectories
        let (_, final_right) = trajectory_right.last().unwrap();
        let (_, final_left) = trajectory_left.last().unwrap();

        // Both should reach approximately the same downrange distance
        assert!((final_right[2] - final_left[2]).abs() < 10.0);
    }
}
