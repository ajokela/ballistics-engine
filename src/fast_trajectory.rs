//! Fast trajectory solver for longer ranges.
//!
//! This is a Rust implementation of the fast fixed-step trajectory solver
//! that provides significant performance improvements for long-range calculations.

use nalgebra::Vector3;
use crate::{
    BallisticInputs, DragModel, BCSegmentData,
    atmosphere::get_local_atmosphere,
    drag::get_drag_coefficient,
    wind::WindSock,
    constants::{MPS_TO_FPS, GRAINS_TO_KG, G_ACCEL_MPS2},
};

/// Fast solution container matching Python implementation
#[derive(Debug, Clone)]
pub struct FastSolution {
    /// Time points
    pub t: Vec<f64>,
    /// State vectors at each time point [6 x n_points]
    pub y: Vec<Vec<f64>>,
    /// Event times [target_hit, max_ord, ground_hit]
    pub t_events: [Vec<f64>; 3],
    /// Whether integration succeeded
    pub success: bool,
}

impl FastSolution {
    /// Interpolate solution at time t
    pub fn sol(&self, t_query: &[f64]) -> Vec<Vec<f64>> {
        let mut result = vec![vec![0.0; t_query.len()]; 6];
        
        for (i, &tq) in t_query.iter().enumerate() {
            // Find the right interval using binary search
            let idx = match self.t.binary_search_by(|&t| t.partial_cmp(&tq).unwrap()) {
                Ok(idx) => idx,
                Err(idx) => idx,
            };
            
            if idx == 0 {
                // Before first point
                for j in 0..6 {
                    result[j][i] = self.y[j][0];
                }
            } else if idx >= self.t.len() {
                // After last point
                for j in 0..6 {
                    result[j][i] = self.y[j][self.t.len() - 1];
                }
            } else {
                // Linear interpolation
                let t0 = self.t[idx - 1];
                let t1 = self.t[idx];
                let frac = (tq - t0) / (t1 - t0);
                
                for j in 0..6 {
                    let y0 = self.y[j][idx - 1];
                    let y1 = self.y[j][idx];
                    result[j][i] = y0 + frac * (y1 - y0);
                }
            }
        }
        
        result
    }
    
    /// Convert from row-major to column-major format for compatibility
    pub fn from_trajectory_data(times: Vec<f64>, states: Vec<[f64; 6]>, t_events: [Vec<f64>; 3]) -> Self {
        let n_points = times.len();
        let mut y = vec![vec![0.0; n_points]; 6];
        
        for (i, state) in states.iter().enumerate() {
            for j in 0..6 {
                y[j][i] = state[j];
            }
        }
        
        FastSolution {
            t: times,
            y,
            t_events,
            success: true,
        }
    }
}

/// Fast trajectory integration parameters
pub struct FastIntegrationParams {
    pub horiz: f64,
    pub vert: f64,
    pub initial_state: [f64; 6],
    pub t_span: (f64, f64),
    pub atmo_params: (f64, f64, f64, f64),
}

/// Fast fixed-step integration for longer trajectories
pub fn fast_integrate(
    inputs: &BallisticInputs,
    wind_sock: &WindSock,
    params: FastIntegrationParams,
) -> FastSolution {
    // Extract parameters
    let mass_kg = inputs.bullet_mass * GRAINS_TO_KG;
    let bc = inputs.bc_value;
    let drag_model = &inputs.bc_type;
    
    // Check for BC segments
    let has_bc_segments = inputs.bc_segments.is_some() && !inputs.bc_segments.as_ref().unwrap().is_empty();
    let has_bc_segments_data = inputs.bc_segments_data.is_some() && !inputs.bc_segments_data.as_ref().unwrap().is_empty();
    
    // Time step - adjust based on distance
    let dt = if params.horiz > 200.0 {
        0.001
    } else if params.horiz > 100.0 {
        0.0005
    } else {
        0.0001
    };
    
    // Maximum time based on estimated flight time
    let v0 = Vector3::new(
        params.initial_state[3],
        params.initial_state[4],
        params.initial_state[5]
    ).norm();
    
    let t_max = if v0 > 1e-6 && params.horiz > 0.0 {
        (2.0 * params.horiz / v0).min(params.t_span.1)
    } else {
        params.t_span.1
    };
    
    // Initialize arrays
    let n_steps = ((t_max / dt) as usize) + 1;
    let mut times = Vec::with_capacity(n_steps);
    let mut states = Vec::with_capacity(n_steps);
    
    // Initial state
    times.push(0.0);
    states.push(params.initial_state);
    
    // Get base atmospheric density
    let (base_density, _) = get_local_atmosphere(
        0.0,
        params.atmo_params.0,
        params.atmo_params.1,
        params.atmo_params.2,
        params.atmo_params.3,
    );
    
    // Integration loop
    let mut hit_target = false;
    let mut hit_ground = false;
    let mut max_ord_time = None;
    let mut max_ord_y = 0.0;
    let ground_threshold = inputs.ground_threshold;
    
    // RK4 integration
    for i in 0..n_steps-1 {
        let t = i as f64 * dt;
        let state = states[i];
        
        let pos = Vector3::new(state[0], state[1], state[2]);
        let vel = Vector3::new(state[3], state[4], state[5]);
        
        // Check termination conditions
        if pos.x >= params.horiz {
            hit_target = true;
            times.push(t);
            states.push(state);
            break;
        }
        
        if pos.y <= ground_threshold {
            hit_ground = true;
            times.push(t);
            states.push(state);
            break;
        }
        
        // Track maximum ordinate
        if pos.y > max_ord_y {
            max_ord_y = pos.y;
            max_ord_time = Some(t);
        }
        
        // RK4 step
        let k1 = compute_derivatives(&state, inputs, wind_sock, base_density, drag_model, bc, has_bc_segments, has_bc_segments_data);
        
        let mut state2 = state;
        for j in 0..6 {
            state2[j] = state[j] + 0.5 * dt * k1[j];
        }
        let k2 = compute_derivatives(&state2, inputs, wind_sock, base_density, drag_model, bc, has_bc_segments, has_bc_segments_data);
        
        let mut state3 = state;
        for j in 0..6 {
            state3[j] = state[j] + 0.5 * dt * k2[j];
        }
        let k3 = compute_derivatives(&state3, inputs, wind_sock, base_density, drag_model, bc, has_bc_segments, has_bc_segments_data);
        
        let mut state4 = state;
        for j in 0..6 {
            state4[j] = state[j] + dt * k3[j];
        }
        let k4 = compute_derivatives(&state4, inputs, wind_sock, base_density, drag_model, bc, has_bc_segments, has_bc_segments_data);
        
        // Update state
        let mut new_state = state;
        for j in 0..6 {
            new_state[j] = state[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
        }
        
        times.push(t + dt);
        states.push(new_state);
    }
    
    // Create event arrays
    let t_events = [
        if hit_target { vec![*times.last().unwrap()] } else { vec![] },
        if let Some(t) = max_ord_time { vec![t] } else { vec![] },
        if hit_ground { vec![*times.last().unwrap()] } else { vec![] },
    ];
    
    FastSolution::from_trajectory_data(times, states, t_events)
}

/// Compute derivatives for the state vector
fn compute_derivatives(
    state: &[f64; 6],
    inputs: &BallisticInputs,
    wind_sock: &WindSock,
    base_density: f64,
    drag_model: &DragModel,
    bc: f64,
    has_bc_segments: bool,
    has_bc_segments_data: bool,
) -> [f64; 6] {
    let pos = Vector3::new(state[0], state[1], state[2]);
    let vel = Vector3::new(state[3], state[4], state[5]);
    
    // Get wind vector
    let wind_vector = wind_sock.vector_for_range_stateless(pos.x);
    
    // Velocity relative to air
    let vel_adjusted = vel - wind_vector;
    let v_mag = vel_adjusted.norm();
    
    // Calculate acceleration
    let accel = if v_mag < 1e-6 {
        Vector3::new(0.0, -G_ACCEL_MPS2, 0.0)
    } else {
        // Calculate drag
        let v_fps = v_mag * MPS_TO_FPS;
        let mach = v_mag / 340.0; // Approximate speed of sound
        
        // Get BC value (potentially from segments)
        let bc_current = if has_bc_segments_data && inputs.bc_segments_data.is_some() {
            get_bc_from_velocity_segments(v_fps, inputs.bc_segments_data.as_ref().unwrap())
        } else if has_bc_segments && inputs.bc_segments.is_some() {
            crate::derivatives::interpolated_bc(mach, inputs.bc_segments.as_ref().unwrap(), Some(inputs))
        } else {
            bc
        };
        
        let drag_factor = get_drag_coefficient(mach, drag_model);
        
        // Calculate drag acceleration using proper ballistics formula
        let cd_to_retard = 0.000683 * 0.30;
        let standard_factor = drag_factor * cd_to_retard;
        let density_scale = base_density / 1.225;
        
        // Drag acceleration in ft/s^2
        let a_drag_ft_s2 = (v_fps * v_fps) * standard_factor * density_scale / bc_current;
        
        // Convert to m/s^2 and apply to velocity vector
        let a_drag_m_s2 = a_drag_ft_s2 * 0.3048; // ft/s^2 to m/s^2
        let accel_drag = -a_drag_m_s2 * (vel_adjusted / v_mag);
        
        // Total acceleration
        accel_drag + Vector3::new(0.0, -G_ACCEL_MPS2, 0.0)
    };
    
    // Return derivatives [vx, vy, vz, ax, ay, az]
    [
        vel.x,
        vel.y,
        vel.z,
        accel.x,
        accel.y,
        accel.z,
    ]
}

/// Get BC from velocity-based segments
fn get_bc_from_velocity_segments(velocity_fps: f64, segments: &[BCSegmentData]) -> f64 {
    for segment in segments {
        if velocity_fps >= segment.velocity_min && velocity_fps <= segment.velocity_max {
            return segment.bc_value;
        }
    }
    
    // If no matching segment, use the BC from the closest segment
    if let Some(first) = segments.first() {
        if velocity_fps < first.velocity_min {
            return first.bc_value;
        }
    }
    
    if let Some(last) = segments.last() {
        if velocity_fps > last.velocity_max {
            return last.bc_value;
        }
    }
    
    // Fallback (shouldn't reach here if segments are properly defined)
    0.5
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_fast_solution_interpolation() {
        let times = vec![0.0, 1.0, 2.0];
        let states = vec![
            [0.0, 0.0, 0.0, 100.0, 50.0, 0.0],
            [100.0, 45.0, 0.0, 99.0, 40.0, 0.0],
            [198.0, 80.0, 0.0, 98.0, 30.0, 0.0],
        ];
        
        let solution = FastSolution::from_trajectory_data(times, states, [vec![], vec![], vec![]]);
        
        // Test interpolation at t=1.5
        let result = solution.sol(&[1.5]);
        
        assert!((result[0][0] - 149.0).abs() < 1e-10); // x position
        assert!((result[1][0] - 62.5).abs() < 1e-10);  // y position
        assert!((result[3][0] - 98.5).abs() < 1e-10);  // vx velocity
    }
    
    #[test]
    fn test_bc_from_velocity_segments() {
        let segments = vec![
            BCSegmentData { velocity_min: 0.0, velocity_max: 1000.0, bc_value: 0.5 },
            BCSegmentData { velocity_min: 1000.0, velocity_max: 2000.0, bc_value: 0.52 },
            BCSegmentData { velocity_min: 2000.0, velocity_max: 3000.0, bc_value: 0.55 },
        ];
        
        assert_eq!(get_bc_from_velocity_segments(500.0, &segments), 0.5);
        assert_eq!(get_bc_from_velocity_segments(1500.0, &segments), 0.52);
        assert_eq!(get_bc_from_velocity_segments(2500.0, &segments), 0.55);
        
        // Test edge cases
        assert_eq!(get_bc_from_velocity_segments(-100.0, &segments), 0.5); // Below min
        assert_eq!(get_bc_from_velocity_segments(3500.0, &segments), 0.55); // Above max
    }
}