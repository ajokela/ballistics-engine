//! Fast trajectory solver for longer ranges.
//!
//! This is a Rust implementation of the fast fixed-step trajectory solver
//! that provides significant performance improvements for long-range calculations.

use crate::{
    atmosphere::get_local_atmosphere,
    constants::{G_ACCEL_MPS2, MPS_TO_FPS},
    drag::get_drag_coefficient,
    wind::WindSock,
    BCSegmentData, DragModel, InternalBallisticInputs as BallisticInputs,
};
use nalgebra::Vector3;

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
            // Use unwrap_or to safely handle NaN values by treating them as greater
            let idx = match self
                .t
                .binary_search_by(|&t| t.partial_cmp(&tq).unwrap_or(std::cmp::Ordering::Greater))
            {
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
    pub fn from_trajectory_data(
        times: Vec<f64>,
        states: Vec<[f64; 6]>,
        t_events: [Vec<f64>; 3],
    ) -> Self {
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
    let _mass_kg = inputs.bullet_mass; // SI (kg)
    let bc = inputs.bc_value;
    let drag_model = &inputs.bc_type;

    // Check for BC segments
    let has_bc_segments =
        inputs.bc_segments.is_some() && !inputs.bc_segments.as_ref().unwrap().is_empty();
    let has_bc_segments_data =
        inputs.bc_segments_data.is_some() && !inputs.bc_segments_data.as_ref().unwrap().is_empty();

    // Time step - adjust based on distance
    let dt = if params.horiz > 200.0 {
        0.001
    } else if params.horiz > 100.0 {
        0.0005
    } else {
        0.0001
    };

    // Maximum integration time. This only bounds the step-array pre-allocation; the
    // hit_target / hit_ground early-breaks below terminate the loop for real shots. Estimate the
    // flight time from the HORIZONTAL velocity with a 4x margin: the previous 2x estimate using
    // the FULL muzzle speed truncated long-range trajectories once drag slowed the bullet (real
    // time of flight to the target far exceeds horiz/v0), so Monte Carlo reported impact metrics
    // at a too-short downrange.
    let vx = params.initial_state[3]; // horizontal (downrange) velocity
    let t_max = if vx > 1e-6 && params.horiz > 0.0 {
        (4.0 * params.horiz / vx).min(params.t_span.1)
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

    // Base drag density = the muzzle (shooter-altitude) density. atmo_params.3 is base_ratio
    // = air_density/1.225 at the shooter altitude (the MC caller computes it via
    // calculate_atmosphere). Previously this called get_local_atmosphere with query alt 0.0
    // while base_alt = shooter altitude, which re-scaled that ratio DOWN to sea level —
    // discarding the correct altitude density and inflating drag for every elevated MC run.
    let base_density = params.atmo_params.3 * 1.225;

    // Integration loop
    let mut hit_target = false;
    let mut hit_ground = false;
    let mut max_ord_time = None;
    let mut max_ord_y = 0.0;
    let ground_threshold = inputs.ground_threshold;

    // RK4 integration
    for i in 0..n_steps - 1 {
        let t = i as f64 * dt;
        let state = states[i];

        let pos = Vector3::new(state[0], state[1], state[2]);
        let _vel = Vector3::new(state[3], state[4], state[5]);

        // Check termination conditions (X is downrange, McCoy)
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
        let k1 = compute_derivatives(
            &state,
            inputs,
            wind_sock,
            base_density,
            drag_model,
            bc,
            has_bc_segments,
            has_bc_segments_data,
        );

        let mut state2 = state;
        for j in 0..6 {
            state2[j] = state[j] + 0.5 * dt * k1[j];
        }
        let k2 = compute_derivatives(
            &state2,
            inputs,
            wind_sock,
            base_density,
            drag_model,
            bc,
            has_bc_segments,
            has_bc_segments_data,
        );

        let mut state3 = state;
        for j in 0..6 {
            state3[j] = state[j] + 0.5 * dt * k2[j];
        }
        let k3 = compute_derivatives(
            &state3,
            inputs,
            wind_sock,
            base_density,
            drag_model,
            bc,
            has_bc_segments,
            has_bc_segments_data,
        );

        let mut state4 = state;
        for j in 0..6 {
            state4[j] = state[j] + dt * k3[j];
        }
        let k4 = compute_derivatives(
            &state4,
            inputs,
            wind_sock,
            base_density,
            drag_model,
            bc,
            has_bc_segments,
            has_bc_segments_data,
        );

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
        if hit_target {
            vec![*times.last().unwrap()]
        } else {
            vec![]
        },
        if let Some(t) = max_ord_time {
            vec![t]
        } else {
            vec![]
        },
        if hit_ground {
            vec![*times.last().unwrap()]
        } else {
            vec![]
        },
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

    // Get wind vector (based on downrange distance, which is X coordinate, McCoy)
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

        // Calculate speed of sound from altitude using standard lapse rate
        // atmo_params: (base_alt, base_temp_c, base_press_hpa, base_ratio)
        let altitude = inputs.altitude + pos.y;
        let (_, speed_of_sound) = get_local_atmosphere(
            altitude,
            inputs.altitude, // base_alt approximation
            inputs.temperature,
            inputs.pressure,
            if inputs.humidity > 0.0 { inputs.humidity } else { 1.0 },
        );
        let mach = v_mag / speed_of_sound;

        // Get BC value (potentially from segments)
        let bc_current = if has_bc_segments_data && inputs.bc_segments_data.is_some() {
            get_bc_from_velocity_segments(v_fps, inputs.bc_segments_data.as_ref().unwrap())
        } else if has_bc_segments && inputs.bc_segments.is_some() {
            crate::derivatives::interpolated_bc(
                mach,
                inputs.bc_segments.as_ref().unwrap(),
                Some(inputs),
            )
        } else {
            bc
        };
        // Guard bc_value == 0 (allowed on FFI/WASM/public MC surfaces): the division below
        // would be Inf -> NaN. Mirrors cli_api's effective_bc.max(1e-6); inert for valid BCs.
        let bc_current = bc_current.max(1e-6);

        // Apply the transonic drag-rise correction once (mirrors derivatives.rs / cli_api) so
        // the Monte Carlo / fast path doesn't under-predict drag near Mach 1. SI fallbacks for
        // caliber/weight (SI-only MC callers may leave the imperial fields 0). wave_drag=false:
        // the G1/G7 tables already embed the rise.
        let base_cd = get_drag_coefficient(mach, drag_model);
        let caliber_in = if inputs.caliber_inches > 0.0 {
            inputs.caliber_inches
        } else {
            inputs.bullet_diameter / 0.0254
        };
        let weight_gr = if inputs.weight_grains > 0.0 {
            inputs.weight_grains
        } else {
            inputs.bullet_mass / 0.00006479891
        };
        let shape =
            crate::transonic_drag::get_projectile_shape(caliber_in, weight_gr, &drag_model.to_string());
        let drag_factor = crate::transonic_drag::transonic_correction(mach, base_cd, shape, false);

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
    [vel.x, vel.y, vel.z, accel.x, accel.y, accel.z]
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

/// Fast integration with explicit wind segments using RK45
/// MBA-155: Upstreamed from ballistics_rust
pub fn fast_integrate_with_segments(
    inputs: &BallisticInputs,
    wind_segments: Vec<crate::wind::WindSegment>,
    params: FastIntegrationParams,
) -> FastSolution {
    // Use the RK45 implementation from trajectory_integration module
    use crate::trajectory_integration::{integrate_trajectory, TrajectoryParams};

    // Extract parameters
    let mass_kg = inputs.bullet_mass; // SI (kg)
    let bc = inputs.bc_value;
    let drag_model = inputs.bc_type;

    // Get omega vector if advanced effects enabled
    let omega_vector = if inputs.enable_advanced_effects {
        // Calculate omega based on latitude and shot azimuth
        // The Earth's rotation vector must be projected into the shooter's
        // local frame which depends on azimuth (shooting direction).
        // azimuth_angle: 0 = North, pi/2 = East
        let omega_earth = 7.2921159e-5; // rad/s
        let lat_rad = inputs.latitude.unwrap_or(0.0).to_radians();
        let azimuth = inputs.azimuth_angle; // already in radians
        Some(Vector3::new(
            omega_earth * lat_rad.cos() * azimuth.cos(), // X: downrange component
            omega_earth * lat_rad.sin(),                 // Y: vertical component
            omega_earth * lat_rad.cos() * azimuth.sin(), // Z: lateral component
        ))
    } else {
        None
    };

    // Set up trajectory parameters
    let traj_params = TrajectoryParams {
        mass_kg,
        bc,
        drag_model,
        wind_segments,
        atmos_params: params.atmo_params,
        omega_vector,
        enable_spin_drift: inputs.enable_advanced_effects,
        enable_magnus: inputs.enable_advanced_effects,
        enable_coriolis: inputs.enable_advanced_effects,
        target_distance_m: params.horiz,
        enable_wind_shear: inputs.enable_wind_shear,
        wind_shear_model: inputs.wind_shear_model.clone(),
        shooter_altitude_m: inputs.altitude,
        is_twist_right: inputs.is_twist_right,
        custom_drag_table: inputs.custom_drag_table.clone(),
        bc_segments: inputs.bc_segments.clone(),
        use_bc_segments: inputs.use_bc_segments,
    };

    // Use RK45 adaptive integration
    let trajectory = integrate_trajectory(
        params.initial_state,
        params.t_span,
        traj_params,
        "RK45", // Use RK45 implementation
        1e-6,   // tolerance
        0.01,   // max_step
    );

    // Convert trajectory to FastSolution format
    let n_points = trajectory.len();
    let mut times = Vec::with_capacity(n_points);
    let mut states = Vec::with_capacity(n_points);

    let mut target_hit_time: Option<f64> = None;
    let mut ground_hit_time: Option<f64> = None;
    let mut max_ord_time = None;
    let mut max_ord_y = 0.0;

    for (t, state_vec) in trajectory {
        // Convert Vector6 to array
        let state = [
            state_vec[0],
            state_vec[1],
            state_vec[2],
            state_vec[3],
            state_vec[4],
            state_vec[5],
        ];

        // Check termination conditions
        // McCoy: state[0]=downrange, state[1]=vertical, state[2]=lateral

        // Record FIRST time target is hit
        if target_hit_time.is_none() && state[0] >= params.horiz {
            target_hit_time = Some(t);
        }

        // Record ground hit
        if ground_hit_time.is_none() && state[1] <= inputs.ground_threshold {
            ground_hit_time = Some(t);
        }

        // Track maximum ordinate
        if state[1] > max_ord_y {
            max_ord_y = state[1];
            max_ord_time = Some(t);
        }

        times.push(t);
        states.push(state);
    }

    // Create event arrays
    let t_events = [
        if let Some(t) = target_hit_time {
            vec![t]
        } else {
            vec![]
        },
        if let Some(t) = max_ord_time {
            vec![t]
        } else {
            vec![]
        },
        if let Some(t) = ground_hit_time {
            vec![t]
        } else {
            vec![]
        },
    ];

    FastSolution::from_trajectory_data(times, states, t_events)
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
        assert!((result[1][0] - 62.5).abs() < 1e-10); // y position
        assert!((result[3][0] - 98.5).abs() < 1e-10); // vx velocity
    }

    #[test]
    fn test_bc_from_velocity_segments() {
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 1000.0,
                bc_value: 0.5,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.52,
            },
            BCSegmentData {
                velocity_min: 2000.0,
                velocity_max: 3000.0,
                bc_value: 0.55,
            },
        ];

        assert_eq!(get_bc_from_velocity_segments(500.0, &segments), 0.5);
        assert_eq!(get_bc_from_velocity_segments(1500.0, &segments), 0.52);
        assert_eq!(get_bc_from_velocity_segments(2500.0, &segments), 0.55);

        // Test edge cases
        assert_eq!(get_bc_from_velocity_segments(-100.0, &segments), 0.5); // Below min
        assert_eq!(get_bc_from_velocity_segments(3500.0, &segments), 0.55); // Above max
    }

    #[test]
    fn test_fast_solution_interpolation_edge_cases() {
        let times = vec![0.0, 1.0, 2.0, 3.0];
        let states = vec![
            [0.0, 0.0, 0.0, 800.0, 50.0, 0.0],
            [800.0, 40.0, 100.0, 750.0, 30.0, 0.0],
            [1550.0, 60.0, 200.0, 700.0, 10.0, 0.0],
            [2250.0, 50.0, 300.0, 650.0, -10.0, 0.0],
        ];

        let solution = FastSolution::from_trajectory_data(times, states, [vec![], vec![], vec![]]);

        // Test interpolation before first point
        let result_before = solution.sol(&[-0.5]);
        assert!((result_before[0][0] - 0.0).abs() < 1e-10); // Should clamp to first

        // Test interpolation after last point
        let result_after = solution.sol(&[5.0]);
        assert!((result_after[0][0] - 2250.0).abs() < 1e-10); // Should clamp to last

        // Test interpolation at exact points
        let result_exact = solution.sol(&[1.0]);
        assert!((result_exact[0][0] - 800.0).abs() < 1e-10);

        // Test multiple query points
        let result_multi = solution.sol(&[0.5, 1.5, 2.5]);
        assert_eq!(result_multi[0].len(), 3);
    }

    #[test]
    fn test_fast_solution_from_trajectory_data() {
        let times = vec![0.0, 0.5, 1.0];
        let states = vec![
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [10.0, 11.0, 12.0, 13.0, 14.0, 15.0],
            [20.0, 21.0, 22.0, 23.0, 24.0, 25.0],
        ];
        let t_events = [vec![1.0], vec![0.5], vec![]];

        let solution = FastSolution::from_trajectory_data(times.clone(), states, t_events);

        // Check that data is stored correctly
        assert_eq!(solution.t, times);
        assert_eq!(solution.y.len(), 6); // 6 state components
        assert_eq!(solution.y[0].len(), 3); // 3 time points
        assert!(solution.success);

        // Verify column-major storage
        assert_eq!(solution.y[0][0], 0.0); // x at t=0
        assert_eq!(solution.y[1][0], 1.0); // y at t=0
        assert_eq!(solution.y[0][2], 20.0); // x at t=1.0
    }

    #[test]
    fn test_bc_segments_boundary_conditions() {
        // Test with single segment
        let single_segment = vec![BCSegmentData {
            velocity_min: 1000.0,
            velocity_max: 2000.0,
            bc_value: 0.5,
        }];

        assert_eq!(get_bc_from_velocity_segments(500.0, &single_segment), 0.5); // Below
        assert_eq!(get_bc_from_velocity_segments(1500.0, &single_segment), 0.5); // In range
        assert_eq!(get_bc_from_velocity_segments(2500.0, &single_segment), 0.5); // Above

        // Test with exact boundary values
        // Note: When velocity matches boundary, first matching segment wins
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 999.0,  // Exclusive upper bound to avoid overlap
                bc_value: 0.45,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.50,
            },
        ];

        assert_eq!(get_bc_from_velocity_segments(1000.0, &segments), 0.50); // At second segment start
        assert_eq!(get_bc_from_velocity_segments(0.0, &segments), 0.45); // At min
        assert_eq!(get_bc_from_velocity_segments(999.0, &segments), 0.45); // At first segment max
    }

    #[test]
    fn test_bc_segments_empty_fallback() {
        let empty_segments: Vec<BCSegmentData> = vec![];

        // With empty segments, should return fallback value
        let result = get_bc_from_velocity_segments(1500.0, &empty_segments);
        assert_eq!(result, 0.5); // Fallback value
    }

    #[test]
    fn test_fast_integration_params() {
        // Verify FastIntegrationParams struct can be constructed
        let params = FastIntegrationParams {
            horiz: 1000.0,
            vert: 0.0,
            initial_state: [0.0, 0.0, 0.0, 800.0, 50.0, 0.0], // McCoy: vx=downrange
            t_span: (0.0, 5.0),
            atmo_params: (0.0, 59.0, 29.92, 0.0),
        };

        assert_eq!(params.horiz, 1000.0);
        assert_eq!(params.t_span.0, 0.0);
        assert_eq!(params.t_span.1, 5.0);
        assert_eq!(params.initial_state[3], 800.0); // vx (downrange, McCoy)
    }

    #[test]
    fn test_fast_solution_event_arrays() {
        let times = vec![0.0, 1.0, 2.0];
        let states = vec![
            [0.0, 0.0, 0.0, 800.0, 50.0, 0.0],
            [800.0, 40.0, 500.0, 750.0, 30.0, 0.0],
            [1500.0, 20.0, 1000.0, 700.0, 10.0, 0.0],
        ];

        // Create solution with events
        let t_events = [
            vec![2.0],  // target_hit at t=2
            vec![0.5],  // max_ord at t=0.5
            vec![],     // no ground_hit
        ];

        let solution = FastSolution::from_trajectory_data(times, states, t_events);

        assert_eq!(solution.t_events[0], vec![2.0]); // Target hit
        assert_eq!(solution.t_events[1], vec![0.5]); // Max ordinate
        assert!(solution.t_events[2].is_empty()); // No ground hit
    }
}
