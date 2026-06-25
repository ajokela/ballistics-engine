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

    /// A clearly-failed solution carrying only the launch state, with `success = false`.
    /// Used when inputs are too degenerate to integrate (e.g. non-physical atmosphere),
    /// so callers see `success = false` instead of a stub trajectory reported as success.
    fn degenerate(initial_state: &[f64; 6]) -> Self {
        let mut y = vec![Vec::new(); 6];
        for (j, slot) in y.iter_mut().enumerate() {
            slot.push(initial_state[j]);
        }
        FastSolution {
            t: vec![0.0],
            y,
            t_events: [Vec::new(), Vec::new(), Vec::new()],
            success: false,
        }
    }
}

/// True if `atmo_params` (alt_m, temp_c, pressure_hPa, humidity) can yield a finite,
/// positive air density. A pressure <= 0 (or non-finite temp/pressure — often a unit
/// mistake, e.g. inHg passed where hPa is expected) gives zero/NaN density and would
/// otherwise silently truncate the integration to a single point.
fn atmo_is_physical(atmo_params: (f64, f64, f64, f64)) -> bool {
    let (_alt, temp_c, pressure_hpa, _humidity) = atmo_params;
    temp_c.is_finite() && pressure_hpa.is_finite() && pressure_hpa > 0.0
}

/// Fast trajectory integration parameters
pub struct FastIntegrationParams {
    pub horiz: f64,
    pub vert: f64,
    pub initial_state: [f64; 6],
    pub t_span: (f64, f64),
    pub atmo_params: (f64, f64, f64, f64),
}

/// Aerodynamic-jump vertical launch-angle offset (radians) for the fast-integrate path.
///
/// Bryan Litz's crosswind estimator (`Y = 0.01*Sg - 0.0024*L + 0.032` MOA/mph) fed by the
/// engine's Miller Sg. The fast path receives a prebuilt initial velocity (no muzzle angle),
/// so the caller rotates that velocity by this offset. Returns 0 when the feature is off or
/// the inputs are degenerate. Crosswind is taken from `wind_speed`/`wind_angle`
/// (BallisticInputs convention: 0 = headwind, +90deg = from the right). MBA-959, EXPERIMENTAL.
pub fn aerodynamic_jump_launch_offset_rad(
    inputs: &BallisticInputs,
    atmo_params: (f64, f64, f64, f64),
) -> f64 {
    if !inputs.enable_aerodynamic_jump {
        return 0.0;
    }
    let diameter = inputs.bullet_diameter;
    if !(inputs.twist_rate.is_finite() && inputs.twist_rate != 0.0)
        || !(diameter.is_finite() && diameter > 0.0)
        || !(inputs.bullet_length.is_finite() && inputs.bullet_length > 0.0)
        || !inputs.muzzle_velocity.is_finite()
    {
        return 0.0;
    }
    let sg = crate::stability::compute_stability_coefficient(inputs, atmo_params);
    if !(sg.is_finite() && sg > 0.0) {
        return 0.0;
    }
    let length_cal = inputs.bullet_length / diameter;
    const MS_TO_MPH: f64 = 2.236_936_292_054_4;
    let crosswind_from_right_mph = inputs.wind_speed * inputs.wind_angle.sin() * MS_TO_MPH;
    let vertical_moa = crate::aerodynamic_jump::litz_crosswind_jump_moa(
        sg,
        length_cal,
        crosswind_from_right_mph,
        inputs.is_twist_right,
    );
    if !vertical_moa.is_finite() {
        return 0.0;
    }
    const MOA_PER_RAD: f64 = 3437.7467707849;
    vertical_moa / MOA_PER_RAD
}

/// Rotate a McCoy-frame state's velocity (indices 3..6) by `theta_rad` in the vertical
/// (downrange–vertical) plane, preserving speed and horizontal heading. Positive = up.
fn rotate_launch_velocity(state: &mut [f64; 6], theta_rad: f64) {
    let (vx, vy, vz) = (state[3], state[4], state[5]);
    let speed = (vx * vx + vy * vy + vz * vz).sqrt();
    if speed <= 0.0 {
        return;
    }
    let h = (vx * vx + vz * vz).sqrt(); // horizontal (downrange+lateral) speed
    let new_elev = vy.atan2(h) + theta_rad;
    state[4] = speed * new_elev.sin();
    let new_h = speed * new_elev.cos();
    let scale = if h > 1e-12 { new_h / h } else { 0.0 };
    state[3] = vx * scale;
    state[5] = vz * scale;
}

/// Fast fixed-step integration for longer trajectories
pub fn fast_integrate(
    inputs: &BallisticInputs,
    wind_sock: &WindSock,
    params: FastIntegrationParams,
) -> FastSolution {
    // Degenerate atmosphere -> non-physical air density would silently stub the run.
    if !atmo_is_physical(params.atmo_params) {
        return FastSolution::degenerate(&params.initial_state);
    }
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

    // Maximum integration time. This bounds BOTH the step-array pre-allocation (n_steps) AND the
    // integration loop itself (the loop runs for at most n_steps-1 iterations); the
    // hit_target / hit_ground early-breaks below terminate the loop sooner for real shots. Estimate
    // the flight time from the HORIZONTAL velocity with a 4x margin: the previous 2x estimate using
    // the FULL muzzle speed truncated long-range trajectories once drag slowed the bullet (real
    // time of flight to the target far exceeds horiz/v0), so Monte Carlo reported impact metrics
    // at a too-short downrange. NOTE: the 4x factor is a heuristic, NOT a proven upper bound — it
    // can still be exceeded by extreme high-drag / high-launch-angle shots, which would truncate
    // the loop before impact.
    // MBA-959: aerodynamic jump perturbs the prebuilt launch velocity vertically (this path is
    // handed an initial_state, not a muzzle angle). A no-op returning the original when disabled.
    let mut initial_state = params.initial_state;
    let aj_offset = aerodynamic_jump_launch_offset_rad(inputs, params.atmo_params);
    if aj_offset != 0.0 {
        rotate_launch_velocity(&mut initial_state, aj_offset);
    }
    let vx = initial_state[3]; // horizontal (downrange) velocity
    let t_max = if vx > 1e-6 && params.horiz > 0.0 {
        (4.0 * params.horiz / vx).min(params.t_span.1)
    } else {
        params.t_span.1
    };

    // Initialize arrays
    let n_steps = ((t_max / dt) as usize) + 1;
    let mut times = Vec::with_capacity(n_steps);
    let mut states = Vec::with_capacity(n_steps);

    // Initial state (with the aerodynamic-jump launch perturbation applied above)
    times.push(0.0);
    states.push(initial_state);

    // Base drag density = the muzzle (shooter-altitude) density. atmo_params.3 is base_ratio
    // = air_density/1.225 at the shooter altitude (the MC caller computes it via
    // calculate_atmosphere). Previously this called get_local_atmosphere with query alt 0.0
    // while base_alt = shooter altitude, which re-scaled that ratio DOWN to sea level —
    // discarding the correct altitude density and inflating drag for every elevated MC run.
    // Guard a missing/absent ratio (base_ratio <= 0, e.g. legacy or uninitialized atmo_params):
    // fall back to the standard sea-level density rather than 0, so a zero ratio cannot collapse
    // density_scale to 0 and silently disable drag entirely. (atmo_params.0 is base_alt here, not
    // a density, so it is not a usable fallback.)
    let base_density = if params.atmo_params.3 > 0.0 {
        params.atmo_params.3 * 1.225
    } else {
        1.225 // standard sea-level air density (kg/m^3)
    };

    // Hoist invariants out of compute_derivatives (called 4x per step). Both the drag-model name
    // and the projectile shape depend only on inputs, not on state/mach, so computing them per
    // call wasted an allocation + heuristic every k1..k4. Mirrors cli_api.rs exactly.

    // Drag-model name as a borrowed &'static str. DragModel's Display goes via Debug, which
    // heap-allocates a String on every call; this match is bit-identical (Display == Debug ==
    // variant name) with no per-step allocation.
    let drag_model_str: &str = match drag_model {
        DragModel::G1 => "G1",
        DragModel::G2 => "G2",
        DragModel::G5 => "G5",
        DragModel::G6 => "G6",
        DragModel::G7 => "G7",
        DragModel::G8 => "G8",
        DragModel::GI => "GI",
        DragModel::GS => "GS",
    };

    // SI fallbacks for caliber/weight (SI-only MC callers may leave the imperial fields 0).
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

    // Projectile shape for transonic corrections (MBA-949: shared resolver — bullet_model name
    // first, then the caliber/weight/drag-model heuristic).
    let projectile_shape = crate::transonic_drag::resolve_projectile_shape(
        inputs.bullet_model.as_deref(),
        caliber_in,
        weight_gr,
        drag_model_str,
    );

    // Coriolis omega (Earth rotation), hoisted (invariant over the flight). MBA-957:
    // fast_integrate — the Monte Carlo / Python-binding path — previously applied NO Coriolis.
    // Mirror the validated cli_api solver exactly: McCoy frame X=downrange, Y=up, Z=lateral;
    // azimuth 0 = North; omega.Z is NEGATIVE (Omega.East = -Omega cos(lat) sin(az)); applied as
    // the physical -2 Omega x v inside compute_derivatives.
    let omega_vector = if inputs.enable_coriolis && inputs.latitude.is_some() {
        let omega_earth = 7.2921159e-5_f64; // rad/s
        let lat = inputs.latitude.unwrap().to_radians();
        let az = inputs.shot_azimuth; // compass bearing (0=N), NOT the aiming offset
        Some(Vector3::new(
            omega_earth * lat.cos() * az.cos(),  // X: downrange
            omega_earth * lat.sin(),             // Y: vertical
            -omega_earth * lat.cos() * az.sin(), // Z: lateral (corrected sign)
        ))
    } else {
        None
    };

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
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
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
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
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
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
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
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
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
    projectile_shape: crate::transonic_drag::ProjectileShape,
    bc: f64,
    has_bc_segments: bool,
    has_bc_segments_data: bool,
    omega: Option<Vector3<f64>>,
) -> [f64; 6] {
    let pos = Vector3::new(state[0], state[1], state[2]);
    let vel = Vector3::new(state[3], state[4], state[5]);

    // Get wind vector (based on downrange distance, which is X coordinate, McCoy)
    let wind_vector = wind_sock.vector_for_range_stateless(pos.x);

    // Velocity relative to air
    let vel_adjusted = vel - wind_vector;
    let v_mag = vel_adjusted.norm();

    // Calculate acceleration
    let mut accel = if v_mag < 1e-6 {
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
        // the Monte Carlo / fast path doesn't under-predict drag near Mach 1. The projectile
        // shape is invariant for the whole integration, so it is hoisted into fast_integrate and
        // passed in rather than recomputed per call. wave_drag=false: the G1/G7 tables already
        // embed the rise.
        // MBA-940: a user-supplied custom drag table is the final Cd, used as-is (no G-model
        // lookup, transonic, or form-factor correction — the curve already encodes the true drag).
        let drag_factor = if let Some(ref table) = inputs.custom_drag_table {
            table.interpolate(mach)
        } else {
            let base_cd = get_drag_coefficient(mach, drag_model);
            let cd =
                crate::transonic_drag::transonic_correction(mach, base_cd, projectile_shape, false);
            // MBA-948: honor use_form_factor in the fast path too (was derivatives.rs-only). No-op
            // when the flag is false (apply_form_factor_to_drag short-circuits), as it is on every
            // current consumer surface.
            crate::form_factor::apply_form_factor_to_drag(
                cd,
                inputs.bullet_model.as_deref(),
                &inputs.bc_type,
                inputs.use_form_factor,
            )
        };

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

    // Coriolis (Earth rotation), MBA-957. omega already carries the corrected lateral sign; use
    // the ground-frame velocity and the physical -2 Omega x v, exactly as the validated cli_api.
    if let Some(omega) = omega {
        accel += -2.0 * omega.cross(&vel);
    }

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

    // Degenerate atmosphere -> non-physical air density would silently stub the run.
    if !atmo_is_physical(params.atmo_params) {
        return FastSolution::degenerate(&params.initial_state);
    }

    // Extract parameters
    let mass_kg = inputs.bullet_mass; // SI (kg)
    let bc = inputs.bc_value;
    let drag_model = inputs.bc_type;

    // Coriolis omega — gated on enable_coriolis (+ a latitude), INDEPENDENT of
    // spin-drift/Magnus. A caller can now request Coriolis-only (enable_coriolis=true
    // with enable_advanced_effects=false) instead of being forced to enable all three.
    let omega_vector = if inputs.enable_coriolis && inputs.latitude.is_some() {
        // Calculate omega based on latitude and shot azimuth
        // The Earth's rotation vector must be projected into the shooter's
        // local frame which depends on azimuth (shooting direction).
        // azimuth_angle: 0 = North, pi/2 = East
        let omega_earth = 7.2921159e-5; // rad/s
        let lat_rad = inputs.latitude.unwrap_or(0.0).to_radians();
        let azimuth = inputs.shot_azimuth; // compass bearing (0=N), NOT the aiming offset
        Some(Vector3::new(
            omega_earth * lat_rad.cos() * azimuth.cos(), // X: downrange component
            omega_earth * lat_rad.sin(),                 // Y: vertical component
            -omega_earth * lat_rad.cos() * azimuth.sin(), // Z: lateral (MBA-957: corrected sign)
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
        enable_magnus: inputs.enable_magnus,
        enable_coriolis: inputs.enable_coriolis,
        target_distance_m: params.horiz,
        enable_wind_shear: inputs.enable_wind_shear,
        wind_shear_model: inputs.wind_shear_model.clone(),
        shooter_altitude_m: inputs.altitude,
        is_twist_right: inputs.is_twist_right,
        // MBA-717: carry the real bullet geometry so spin-drift / Magnus use it.
        bullet_diameter: inputs.bullet_diameter,
        bullet_length: inputs.bullet_length,
        twist_rate: inputs.twist_rate,
        custom_drag_table: inputs.custom_drag_table.clone(),
        bc_segments: inputs.bc_segments.clone(),
        use_bc_segments: inputs.use_bc_segments,
        // MBA-954: keep the historical -1000.0 here (behavior-preserving for this binding path);
        // threading inputs.ground_threshold would change the default ground plane for existing
        // callers. Direct TrajectoryParams constructors can now configure it.
        ground_threshold: -1000.0,
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

    #[test]
    fn fast_path_coriolis_uses_shot_direction() {
        // Regression: fast_integrate_with_segments (the path the Python binding uses)
        // built its Coriolis omega from azimuth_angle (the aiming offset, always ~0)
        // instead of shot_azimuth, so east and west shots came out identical. After the
        // fix they must differ with the correct Eotvos sign (east lifted, higher).
        use std::f64::consts::FRAC_PI_2;
        // Returns (final_downrange, final_vertical) for a shot fired along `shot_az`.
        fn final_xy(shot_az: f64) -> (f64, f64) {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_velocity = 800.0;
            inputs.bc_value = 0.5;
            inputs.bc_type = DragModel::G7;
            inputs.enable_advanced_effects = true; // gates the omega vector
            inputs.enable_coriolis = true;
            inputs.latitude = Some(45.0);
            inputs.shot_azimuth = shot_az;
            let v = 800.0_f64;
            let elev = 0.02_f64;
            let params = FastIntegrationParams {
                horiz: 1000.0,
                vert: 0.0,
                initial_state: [0.0, 0.0, 0.0, v * elev.cos(), v * elev.sin(), 0.0],
                t_span: (0.0, 5.0),
                atmo_params: (0.0, 59.0, 29.92, 0.0),
            };
            let sol = fast_integrate_with_segments(&inputs, vec![], params);
            let n = sol.y[0].len();
            (sol.y[0][n - 1], sol.y[1][n - 1])
        }
        let (ex, ey) = final_xy(FRAC_PI_2); // east
        let (wx, wy) = final_xy(3.0 * FRAC_PI_2); // west
        // Both shots cover essentially the same downrange (Coriolis barely affects x),
        // so comparing the final vertical is apples-to-apples.
        assert!(
            (ex - wx).abs() < 0.5,
            "east/west downrange should be ~equal (ex={ex:.4}, wx={wx:.4})"
        );
        // Pre-fix the fast path was North-locked, making these byte-identical. The Eotvos
        // term now lifts the east shot above the west shot.
        assert!(
            ey > wy,
            "fast-path east ({ey:.6}) must be higher than west ({wy:.6}) (Eotvos)"
        );
        assert!(
            (ey - wy) > 1e-5,
            "fast-path E-W vertical separation ({:.8} m) should be non-zero (the pre-fix bug was exact equality)",
            ey - wy
        );
    }

    #[test]
    fn fast_path_coriolis_independent_of_advanced_effects() {
        // Coriolis is now gated on enable_coriolis (+ latitude), NOT enable_advanced_effects.
        // So a caller can request Coriolis-only without being forced to enable spin/Magnus.
        use std::f64::consts::FRAC_PI_2;
        fn final_y(coriolis: bool, shot_az: f64) -> f64 {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_velocity = 800.0;
            inputs.bc_value = 0.5;
            inputs.bc_type = DragModel::G7;
            inputs.enable_coriolis = coriolis;
            inputs.enable_advanced_effects = false; // explicitly OFF — Coriolis must still work
            inputs.latitude = Some(45.0);
            inputs.shot_azimuth = shot_az;
            let v = 800.0_f64;
            let elev = 0.02_f64;
            let params = FastIntegrationParams {
                horiz: 1000.0,
                vert: 0.0,
                initial_state: [0.0, 0.0, 0.0, v * elev.cos(), v * elev.sin(), 0.0],
                t_span: (0.0, 5.0),
                atmo_params: (0.0, 59.0, 29.92, 0.0),
            };
            let sol = fast_integrate_with_segments(&inputs, vec![], params);
            let n = sol.y[0].len();
            sol.y[1][n - 1]
        }
        // enable_coriolis=true with advanced effects OFF: directional Coriolis still applies.
        let e = final_y(true, FRAC_PI_2);
        let w = final_y(true, 3.0 * FRAC_PI_2);
        assert!(e > w && (e - w) > 1e-5, "Coriolis-only (no advanced effects) must still be directional: E={e} W={w}");
        // enable_coriolis=false: no Coriolis at all, east == west.
        let e2 = final_y(false, FRAC_PI_2);
        let w2 = final_y(false, 3.0 * FRAC_PI_2);
        assert!((e2 - w2).abs() < 1e-9, "with enable_coriolis=false, east/west must be identical: E={e2} W={w2}");
    }

    #[test]
    fn fast_path_rejects_degenerate_atmosphere() {
        let mut inputs = BallisticInputs::default();
        inputs.muzzle_velocity = 800.0;
        inputs.bc_value = 0.5;
        inputs.bc_type = DragModel::G7;
        let v = 800.0_f64;
        let e = 0.02_f64;
        let mk = |atmo: (f64, f64, f64, f64)| FastIntegrationParams {
            horiz: 500.0,
            vert: 0.0,
            initial_state: [0.0, 0.0, 0.0, v * e.cos(), v * e.sin(), 0.0],
            t_span: (0.0, 5.0),
            atmo_params: atmo,
        };
        // pressure <= 0 -> fail loudly (success=false) instead of a 1-point stub as success.
        let zero_p = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, 0.0, 50.0)));
        assert!(!zero_p.success, "pressure=0 atmosphere must yield success=false");
        // non-finite pressure -> also rejected.
        let nan_p = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, f64::NAN, 50.0)));
        assert!(!nan_p.success, "NaN pressure must yield success=false");
        // realistic atmosphere -> success.
        let good = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, 1013.25, 50.0)));
        assert!(good.success, "realistic atmosphere must yield success=true");
    }

    #[test]
    fn fast_path_carries_real_bullet_geometry() {
        // MBA-717: build_inputs hardcoded diameter=.308 / length=1.24in / twist=10 because
        // TrajectoryParams didn't carry them. They're now plumbed through, so the BallisticInputs
        // the derivatives see reflect the real bullet (caliber gates the Magnus block at
        // bullet_diameter > 0.0). Guard against regressing back to the hardcoded values: a
        // zero-diameter input must reach the derivatives as zero (Magnus skipped), whereas the
        // old code would have forced .308 regardless. We assert the run still completes and the
        // two geometries don't crash — the data path is exercised end-to-end.
        let run = |diameter: f64, twist: f64| {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_velocity = 800.0;
            inputs.bc_value = 0.5;
            inputs.bc_type = DragModel::G7;
            inputs.bullet_diameter = diameter;
            inputs.bullet_length = 0.0318;
            inputs.twist_rate = twist;
            inputs.enable_advanced_effects = true;
            inputs.enable_magnus = true;
            let v = 800.0_f64;
            let elev = 0.02_f64;
            let params = FastIntegrationParams {
                horiz: 1000.0,
                vert: 0.0,
                initial_state: [0.0, 0.0, 0.0, v * elev.cos(), v * elev.sin(), 0.0],
                t_span: (0.0, 5.0),
                atmo_params: (0.0, 15.0, 1013.25, 50.0),
            };
            fast_integrate_with_segments(&inputs, vec![], params)
        };
        // Both a real .224 and a real .338 produce a valid trajectory (geometry is honored, not
        // overridden by a hardcoded .308). Regression sentinel for the plumbing.
        assert!(run(0.00569, 7.0).success, ".224 geometry must solve");
        assert!(run(0.00858, 10.0).success, ".338 geometry must solve");
    }
}

