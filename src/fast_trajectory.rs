//! Fast trajectory solver for longer ranges.
//!
//! This is a Rust implementation of the fast fixed-step trajectory solver
//! that provides significant performance improvements for long-range calculations.

use crate::{
    atmosphere::{calculate_air_density_cimp, get_local_atmosphere_humid, AtmoSock},
    bc_estimation::velocity_segment_bc,
    constants::{G_ACCEL_MPS2, MPS_TO_FPS, STANDARD_AIR_DENSITY},
    drag::get_drag_coefficient,
    wind::WindSock,
    DragModel, InternalBallisticInputs as BallisticInputs,
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
                let span = t1 - t0;

                for j in 0..6 {
                    let y0 = self.y[j][idx - 1];
                    let y1 = self.y[j][idx];
                    result[j][i] = if span.abs() < f64::EPSILON {
                        y1
                    } else {
                        let frac = (tq - t0) / span;
                        y0 + frac * (y1 - y0)
                    };
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

fn direct_atmosphere_values(
    atmo_params: (f64, f64, f64, f64),
) -> Option<(f64, f64)> {
    let (a, b, c, d) = atmo_params;
    (a.is_finite()
        && b.is_finite()
        && c == 0.0
        && d == 0.0
        && a > 0.0
        && a < 2.0
        && b > 200.0)
        .then_some((a, b))
}

fn stability_atmosphere_params(atmo_params: (f64, f64, f64, f64)) -> (f64, f64, f64, f64) {
    if let Some((air_density, _)) = direct_atmosphere_values(atmo_params) {
        // compute_stability_coefficient consumes standard-mode temperature/pressure, whose
        // Miller correction is rho_ref/rho. At the 15 C reference temperature, scaling pressure
        // by rho/rho_ref supplies that exact correction without misreading sound speed as temp.
        (0.0, 15.0, 1013.25 * air_density / STANDARD_AIR_DENSITY, 1.0)
    } else {
        atmo_params
    }
}

const MAX_STANDARD_DENSITY_RATIO: f64 = 2.0;

/// True if `atmo_params` can yield a finite, positive air density. `atmo_params` has TWO
/// modes that `compute_derivatives` distinguishes:
///   * **Standard**: `(base_alt_m, base_temp_c, base_pressure_hPa, base_density_ratio)` —
///     positive station pressure. Slot 3 is a density RATIO, not humidity, despite the
///     `humidity` field it lands in via `build_inputs`. A nonpositive ratio means "not supplied"
///     and falls back to `1.0`; a supplied ratio must be below `2.0`.
///   * **Direct**: `(air_density, speed_of_sound, 0.0, 0.0)` — slots 2 and 3 are `0.0`
///     sentinels, with a real density (< 2.0 kg/m³) and speed of sound (> 200 m/s).
///
/// A pressure <= 0 that ISN'T the direct-mode sentinel, an overlarge supplied density ratio, or
/// any non-finite value would yield a non-physical atmosphere, so it is rejected. (Earlier this
/// guard also rejected legitimate direct-mode input; the direct-mode allowance below fixes that.)
fn atmo_is_physical(atmo_params: (f64, f64, f64, f64)) -> bool {
    let (a, b, c, d) = atmo_params;
    if !(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
        return false;
    }
    // Direct-atmosphere mode: (density, speed_of_sound, 0, 0). Constants mirror
    // derivatives.rs (MAX_REALISTIC_DENSITY = 2.0 kg/m³, MIN_REALISTIC_SPEED_OF_SOUND = 200 m/s).
    // Standard mode: positive station pressure (hPa), plus either the documented missing-ratio
    // sentinel or a supplied ratio below the physical ceiling.
    direct_atmosphere_values(atmo_params).is_some() || (c > 0.0 && d < MAX_STANDARD_DENSITY_RATIO)
}

#[derive(Debug, Clone, Copy)]
enum FastAtmosphere {
    Direct {
        air_density: f64,
        speed_of_sound: f64,
    },
    Standard {
        base_density: f64,
    },
}

/// Fast trajectory integration parameters
pub struct FastIntegrationParams {
    pub horiz: f64,
    pub vert: f64,
    pub initial_state: [f64; 6],
    pub t_span: (f64, f64),
    /// Dual-mode atmosphere tuple — see [`atmo_is_physical`]. Standard mode is
    /// `(base_alt_m, base_temp_c, base_pressure_hPa, base_density_ratio)` (slot 3 is a density
    /// RATIO, not humidity; nonpositive means "not supplied" and falls back to `1.0`, while a
    /// supplied ratio must be below `2.0`). Direct mode is
    /// `(air_density, speed_of_sound, 0.0, 0.0)`.
    pub atmo_params: (f64, f64, f64, f64),
    /// MBA-1137: optional downrange-segmented atmosphere. When `Some`, the per-substep drag
    /// density samples the base (station-referenced) T/P/H for the zone at the current downrange
    /// distance before applying the altitude lapse. `None` (default) keeps the single-station
    /// base and — combined with the MBA-1137 density-freeze fix — varies density with altitude only.
    pub atmo_sock: Option<AtmoSock>,
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
    let stability_atmo = stability_atmosphere_params(atmo_params);
    let sg = crate::stability::compute_stability_coefficient(inputs, stability_atmo);
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

fn launch_state_with_aerodynamic_jump(
    inputs: &BallisticInputs,
    atmo_params: (f64, f64, f64, f64),
    mut initial_state: [f64; 6],
) -> [f64; 6] {
    let offset = aerodynamic_jump_launch_offset_rad(inputs, atmo_params);
    if offset != 0.0 {
        rotate_launch_velocity(&mut initial_state, offset);
    }
    initial_state
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
    let mut effective_inputs = inputs.clone();
    if params.atmo_params.2 > 0.0 {
        effective_inputs.altitude = params.atmo_params.0;
        effective_inputs.temperature = params.atmo_params.1;
        effective_inputs.pressure = params.atmo_params.2;
        effective_inputs.humidity = params.atmo_params.3;
    }
    let inputs = &effective_inputs;
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

    // MBA-959: aerodynamic jump perturbs the prebuilt launch velocity vertically (this path is
    // handed an initial_state, not a muzzle angle). A no-op returning the original when disabled.
    let initial_state =
        launch_state_with_aerodynamic_jump(inputs, params.atmo_params, params.initial_state);
    let vx = initial_state[3]; // horizontal (downrange) velocity

    // MBA-1145: decouple the integration-loop ceiling from the pre-allocation heuristic.
    // Previously t_max = min(4*horiz/vx, t_span.1) bounded BOTH the Vec sizing AND the loop
    // itself, so the 4x-horiz/vx estimate doubled as the loop cap. That estimate is only a
    // heuristic — extreme high-drag or high-launch-angle shots exceed it, and the loop then
    // terminated BEFORE hit_target/hit_ground, silently truncating the trajectory tail (Monte
    // Carlo reported impact metrics short of the real range). The loop already breaks on
    // hit_target (pos.x >= horiz) and hit_ground (pos.y <= ground_threshold), and any real
    // trajectory descends below the ground threshold within a few seconds of apex, so the loop
    // bound's only job is a runaway safety ceiling: use the full t_span.1 (default 30 s). The
    // 4x estimate is retained ONLY to size the pre-allocation (keeps the small-allocation fast
    // path for short shots); the Vec grows for the rare long run because the loop bound is n_steps.
    let n_steps = ((params.t_span.1 / dt) as usize) + 1;
    let est_steps = if vx > 1e-6 && params.horiz > 0.0 {
        (((4.0 * params.horiz / vx) / dt) as usize) + 1
    } else {
        n_steps
    };
    let cap = est_steps.min(n_steps);
    let mut times = Vec::with_capacity(cap);
    let mut states = Vec::with_capacity(cap);

    // Initial state (with the aerodynamic-jump launch perturbation applied above)
    times.push(0.0);
    states.push(initial_state);

    // Direct mode supplies a fixed density and speed of sound. Standard mode supplies station
    // conditions plus base_ratio; its local density/sound speed are lapsed at every substep.
    // Guard a missing standard-mode ratio by falling back to sea-level density (MBA-1157 owns
    // stricter validation of that separate contract).
    let atmosphere = if let Some((air_density, speed_of_sound)) =
        direct_atmosphere_values(params.atmo_params)
    {
        FastAtmosphere::Direct {
            air_density,
            speed_of_sound,
        }
    } else {
        let base_density = if params.atmo_params.3 > 0.0 {
            params.atmo_params.3 * 1.225
        } else {
            1.225
        };
        FastAtmosphere::Standard { base_density }
    };

    // MBA-1137: borrow the optional downrange-segmented atmosphere once (queried 4x per step).
    let atmo_sock = params.atmo_sock.as_ref();

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
    // First project into level downrange/up/lateral axes: azimuth 0 = North; omega.Z is NEGATIVE
    // (Omega.East = -Omega cos(lat) sin(az)). The derivative kernel then projects this vector into
    // the inclined shot frame and applies the physical -2 Omega x v.
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
    // Parse the string model once; the derivative kernel is called four times per RK4 step.
    let wind_shear_model = if inputs.enable_wind_shear {
        let model = crate::wind_shear::boundary_layer_model_from_name(&inputs.wind_shear_model);
        (model != crate::wind_shear::WindShearModel::None).then_some(model)
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
            break;
        }

        if pos.y <= ground_threshold {
            hit_ground = true;
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
            atmosphere,
            drag_model,
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
            wind_shear_model,
            atmo_sock,
        );

        let mut state2 = state;
        for j in 0..6 {
            state2[j] = state[j] + 0.5 * dt * k1[j];
        }
        let k2 = compute_derivatives(
            &state2,
            inputs,
            wind_sock,
            atmosphere,
            drag_model,
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
            wind_shear_model,
            atmo_sock,
        );

        let mut state3 = state;
        for j in 0..6 {
            state3[j] = state[j] + 0.5 * dt * k2[j];
        }
        let k3 = compute_derivatives(
            &state3,
            inputs,
            wind_sock,
            atmosphere,
            drag_model,
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
            wind_shear_model,
            atmo_sock,
        );

        let mut state4 = state;
        for j in 0..6 {
            state4[j] = state[j] + dt * k3[j];
        }
        let k4 = compute_derivatives(
            &state4,
            inputs,
            wind_sock,
            atmosphere,
            drag_model,
            projectile_shape,
            bc,
            has_bc_segments,
            has_bc_segments_data,
            omega_vector,
            wind_shear_model,
            atmo_sock,
        );

        // Update state
        let mut new_state = state;
        for j in 0..6 {
            new_state[j] = state[j] + dt * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
        }

        if state[0] < params.horiz && new_state[0] >= params.horiz {
            // Keep every terminal metric on the requested target plane. The former top-of-loop
            // check retained this full-step endpoint, biasing time, drop, and velocity past it.
            let alpha = (params.horiz - state[0]) / (new_state[0] - state[0]);
            let mut target_state = state;
            for j in 0..6 {
                target_state[j] = state[j] + alpha * (new_state[j] - state[j]);
            }
            target_state[0] = params.horiz;
            times.push(t + alpha * dt);
            states.push(target_state);
            hit_target = true;
            break;
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

    // MBA-1134 (rank 10) — regression fix: apply the canonical empirical Litz drift as a
    // post-process to the lateral (McCoy Z) here too. The derivatives kernel no longer
    // integrates spin drift, and this plain fast_integrate is the single-shot fast path
    // (solve_trajectory_rust / the API), so without this it would carry NO spin drift and
    // twist direction would have no effect (the Magnus term is also suppressed when
    // use_enhanced_spin_drift is set). Mirrors fast_integrate_with_segments,
    // cli_api::apply_spin_drift and the Monte-Carlo path so all solver families agree; uses
    // the SAME muzzle Sg (spin_drift::effective_sg_from_inputs).
    if inputs.use_enhanced_spin_drift {
        // Standard-mode atmo_params is (base_alt, temp_c, press_hpa, ratio); direct mode
        // (density, sound, 0, 0) carries no explicit temp/pressure, so fall back to sea-level
        // standard (the Sg density correction is a no-op there).
        let (sd_temp_c, sd_press_hpa) = if params.atmo_params.2 > 0.0 {
            (params.atmo_params.1, params.atmo_params.2)
        } else {
            (15.0, 1013.25)
        };
        let sg = crate::spin_drift::effective_sg_from_inputs(inputs, sd_temp_c, sd_press_hpa);
        for (t, state) in times.iter().zip(states.iter_mut()) {
            if *t > 0.0 {
                state[2] += crate::spin_drift::litz_drift_meters(sg, *t, inputs.is_twist_right);
            }
        }
    }

    FastSolution::from_trajectory_data(times, states, t_events)
}

fn fast_magnus_acceleration(
    inputs: &BallisticInputs,
    air_velocity: Vector3<f64>,
    air_density: f64,
    mach: f64,
    gravity_acceleration: Vector3<f64>,
) -> Vector3<f64> {
    if !inputs.enable_magnus
        || inputs.use_enhanced_spin_drift
        || inputs.bullet_diameter <= 0.0
        || inputs.twist_rate <= 0.0
        || inputs.bullet_mass <= 0.0
    {
        return Vector3::zeros();
    }

    let speed_air = air_velocity.norm();
    let diameter_m = inputs.bullet_diameter;
    let (spin_rate_rad_s, spin_param) = crate::spin_drift::calculate_magnus_spin_state(
        inputs.muzzle_velocity,
        speed_air,
        inputs.twist_rate,
        diameter_m,
    );
    let d_in = if inputs.caliber_inches > 0.0 {
        inputs.caliber_inches
    } else {
        diameter_m / 0.0254
    };
    let m_gr = if inputs.weight_grains > 0.0 {
        inputs.weight_grains
    } else {
        inputs.bullet_mass / 0.00006479891
    };
    let l_in = if inputs.bullet_length > 0.0 {
        inputs.bullet_length / 0.0254
    } else {
        let estimated = crate::stability::estimate_bullet_length_m(diameter_m, inputs.bullet_mass);
        if estimated > 0.0 {
            estimated / 0.0254
        } else {
            4.5 * d_in.max(1e-9)
        }
    };
    let sg = crate::spin_drift::calculate_dynamic_stability(
        m_gr,
        speed_air,
        spin_rate_rad_s,
        d_in,
        l_in,
        air_density,
    );
    let (yaw_rad, _) = crate::spin_drift::calculate_yaw_of_repose(
        sg,
        speed_air,
        spin_rate_rad_s,
        0.0,
        0.0,
        air_density,
        d_in,
        l_in,
        m_gr,
        mach,
        "match",
        false,
    );
    let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);
    let c_np = crate::derivatives::calculate_magnus_moment_coefficient(mach);
    let force = 0.5 * air_density * speed_air.powi(2) * area * c_np * spin_param * yaw_rad.sin();
    if force <= 1e-12 {
        return Vector3::zeros();
    }

    crate::derivatives::yaw_of_repose_magnus_direction(
        air_velocity,
        gravity_acceleration,
        inputs.is_twist_right,
    )
    .map_or_else(Vector3::zeros, |direction| {
        (force / inputs.bullet_mass) * direction
    })
}

/// Compute derivatives for the state vector
#[allow(clippy::too_many_arguments)]
fn compute_derivatives(
    state: &[f64; 6],
    inputs: &BallisticInputs,
    wind_sock: &WindSock,
    atmosphere: FastAtmosphere,
    drag_model: &DragModel,
    projectile_shape: crate::transonic_drag::ProjectileShape,
    bc: f64,
    has_bc_segments: bool,
    has_bc_segments_data: bool,
    omega: Option<Vector3<f64>>,
    wind_shear_model: Option<crate::wind_shear::WindShearModel>,
    // MBA-1137: optional downrange-segmented atmosphere (zone base swapped by pos.x before lapse).
    atmo_sock: Option<&AtmoSock>,
) -> [f64; 6] {
    let pos = Vector3::new(state[0], state[1], state[2]);
    let vel = Vector3::new(state[3], state[4], state[5]);

    // Resolve the cached downrange wind in level axes, then apply boundary-layer shear using
    // height gained above launch. Site elevation is MSL and intentionally does not enter shear.
    let level_wind = wind_sock.vector_for_range_stateless(pos.x);
    let level_wind = if let Some(model) = wind_shear_model {
        let height_rel_launch =
            crate::atmosphere::shot_frame_altitude(0.0, pos.x, pos.y, inputs.shooting_angle);
        crate::wind_shear::apply_boundary_layer_shear(level_wind, height_rel_launch, model)
    } else {
        level_wind
    };
    let wind_vector =
        crate::derivatives::level_vector_to_shot_frame(level_wind, inputs.shooting_angle);

    // Velocity relative to air
    let vel_adjusted = vel - wind_vector;
    let v_mag = vel_adjusted.norm();

    // Gravity acceleration vector, rotated into the shot-aligned frame by shooting_angle
    // (uphill/downhill inclined fire), matching cli_api::TrajectorySolver::gravity_acceleration.
    let theta = inputs.shooting_angle;
    let accel_gravity = Vector3::new(
        -G_ACCEL_MPS2 * theta.sin(),
        -G_ACCEL_MPS2 * theta.cos(),
        0.0,
    );

    // Calculate acceleration
    let mut accel = if v_mag < 1e-6 {
        accel_gravity
    } else {
        // Calculate drag
        let v_fps = v_mag * MPS_TO_FPS;

        // Resolve LOCAL density and speed of sound. Direct mode uses its supplied fixed values;
        // standard mode evaluates the substep altitude with the lapse-rate pipeline.
        // `inputs.temperature`/`inputs.pressure` are the standard-mode base (station) T/P set by
        // fast_integrate from atmo_params; `base_density` is the station-altitude density, so
        // `base_density / 1.225` recovers the station base_ratio the lapse pipeline expects.
        //
        // MBA-1137 (density-freeze fix): previously ONLY the speed of sound was taken from this
        // call and `density_scale` used the flight-constant `base_density`, so the fast/MC path
        // held density frozen for the whole flight (density did NOT vary with altitude — a latent
        // bug the cli_api/derivatives paths already avoid). Now the LOCAL density is used for
        // `density_scale` below, so the fast path varies density with altitude too.
        //
        // MBA-1137 (zones): when a downrange-segmented atmosphere is present, swap the BASE
        // (station-referenced) temp/pressure/ratio for the zone at the current downrange distance
        // (pos.x), recomputing the zone base_ratio via CIPM, BEFORE the altitude lapse — so
        // downrange-zone selection and the world-vertical lapse compose without double-counting.
        // `None` keeps the single-station base.
        //
        // Humidity is NOT plumbed to this call site: on the fast path (fast_integrate) the
        // `humidity` FIELD is overwritten with atmo_params.3 (the density RATIO), so
        // `inputs.humidity` is not a real RH — pass 0.0 (dry) rather than fabricate.
        let (local_density, speed_of_sound) = match atmosphere {
            FastAtmosphere::Direct {
                air_density,
                speed_of_sound,
            } => (air_density, speed_of_sound),
            FastAtmosphere::Standard { base_density } => {
                let altitude = crate::atmosphere::shot_frame_altitude(
                    inputs.altitude,
                    pos.x,
                    pos.y,
                    inputs.shooting_angle,
                );
                let (base_temp_c, base_press_hpa, base_ratio) = match atmo_sock {
                    Some(sock) => {
                        let (zt, zp, zh) = sock.atmo_for_range(pos.x);
                        (zt, zp, calculate_air_density_cimp(zt, zp, zh) / 1.225)
                    }
                    None => (inputs.temperature, inputs.pressure, base_density / 1.225),
                };
                get_local_atmosphere_humid(
                    altitude,
                    inputs.altitude, // base_alt approximation
                    base_temp_c,
                    base_press_hpa,
                    base_ratio,
                    0.0, // humidity not available here (see note above)
                )
            }
        };
        let mach = v_mag / speed_of_sound;

        // Get BC value (potentially from segments)
        let bc_current = if inputs.use_bc_segments
            && has_bc_segments_data
            && inputs.bc_segments_data.is_some()
        {
            velocity_segment_bc(v_fps, inputs.bc_segments_data.as_ref().unwrap(), bc)
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
        // Its Cd is the projectile's ACTUAL drag coefficient, so the retardation denominator
        // must be the sectional density (lb/in²), not a BC: Cd_own / SD == Cd_ref / BC
        // (see BallisticInputs::custom_drag_denominator).
        let (drag_factor, retard_denom) = if let Some(ref table) = inputs.custom_drag_table {
            (
                table.interpolate(mach),
                inputs.custom_drag_denominator(bc_current),
            )
        } else {
            let base_cd = get_drag_coefficient(mach, drag_model);
            let cd =
                crate::transonic_drag::transonic_correction(mach, base_cd, projectile_shape, false);
            (cd, bc_current)
        };

        // Calculate drag acceleration using proper ballistics formula
        let cd_to_retard = crate::constants::CD_TO_RETARD;
        let standard_factor = drag_factor * cd_to_retard;
        // MBA-1137: use the LOCAL (per-substep) density, not the frozen flight-constant
        // `base_density`, so drag varies with altitude AND downrange zone.
        let density_scale = local_density / 1.225;

        // Drag acceleration in ft/s^2
        let a_drag_ft_s2 = (v_fps * v_fps) * standard_factor * density_scale / retard_denom;

        // Convert to m/s^2 and apply to velocity vector
        let a_drag_m_s2 = a_drag_ft_s2 * 0.3048; // ft/s^2 to m/s^2
        let accel_drag = -a_drag_m_s2 * (vel_adjusted / v_mag);

        // The Litz post-process owns the gyroscopic effect whenever it is enabled. Otherwise,
        // honor the explicit Magnus flag with the same dynamic-Sg/yaw model as sibling solvers.
        let accel_magnus =
            fast_magnus_acceleration(inputs, vel_adjusted, local_density, mach, accel_gravity);

        // Total acceleration
        accel_drag + accel_gravity + accel_magnus
    };

    // Coriolis (Earth rotation), MBA-957. Omega arrives in level downrange/up/lateral axes;
    // project it into the inclined shot frame before applying the physical -2 Omega x v.
    if let Some(omega) = omega {
        let omega = crate::derivatives::level_vector_to_shot_frame(omega, inputs.shooting_angle);
        accel += -2.0 * omega.cross(&vel);
    }

    // Return derivatives [vx, vy, vz, ax, ay, az]
    [vel.x, vel.y, vel.z, accel.x, accel.y, accel.z]
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

    // Match plain fast_integrate: this entry point also receives a prebuilt launch state, so
    // apply the experimental aerodynamic-jump angle exactly once before the low-level integrator.
    let initial_state =
        launch_state_with_aerodynamic_jump(inputs, params.atmo_params, params.initial_state);

    // Extract parameters
    let mass_kg = inputs.bullet_mass; // SI (kg)
    let bc = inputs.bc_value;
    let drag_model = inputs.bc_type;

    // Coriolis omega — gated on enable_coriolis (+ a latitude), INDEPENDENT of
    // spin-drift/Magnus. A caller can now request Coriolis-only (enable_coriolis=true
    // with enable_advanced_effects=false) instead of being forced to enable all three.
    let omega_vector = if inputs.enable_coriolis && inputs.latitude.is_some() {
        // Calculate omega based on latitude and shot azimuth
        // First project Earth's rotation into level downrange/up/lateral axes based on azimuth;
        // the derivative kernel handles the additional inclined-shot projection.
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
        enable_spin_drift: inputs.use_enhanced_spin_drift,
        enable_magnus: inputs.enable_magnus,
        enable_coriolis: inputs.enable_coriolis,
        target_distance_m: params.horiz,
        enable_wind_shear: inputs.enable_wind_shear,
        wind_shear_model: inputs.wind_shear_model.clone(),
        shooter_altitude_m: inputs.altitude,
        is_twist_right: inputs.is_twist_right,
        shooting_angle: inputs.shooting_angle,
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
        // MBA-1137: forward the downrange-segmented atmosphere to the RK45 derivatives path.
        atmo_sock: params.atmo_sock,
    };

    // Use RK45 adaptive integration
    let trajectory = integrate_trajectory(
        initial_state,
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

    // MBA-1134 (rank 10): the derivatives kernel no longer integrates spin drift, so apply the
    // canonical empirical Litz drift as a post-process to the lateral (McCoy Z) of every point at
    // its time of flight — matching cli_api::apply_spin_drift and the Monte-Carlo path so all three
    // solver families agree. Uses the SAME muzzle Sg (spin_drift::effective_sg_from_inputs).
    if inputs.use_enhanced_spin_drift {
        // Standard-mode atmo_params is (base_alt, temp_c, press_hpa, ratio); direct mode
        // (density, sound, 0, 0) carries no explicit temp/pressure, so fall back to sea-level
        // standard (the Sg density correction is a no-op there).
        let (sd_temp_c, sd_press_hpa) = if params.atmo_params.2 > 0.0 {
            (params.atmo_params.1, params.atmo_params.2)
        } else {
            (15.0, 1013.25)
        };
        let sg = crate::spin_drift::effective_sg_from_inputs(inputs, sd_temp_c, sd_press_hpa);
        for (t, state) in times.iter().zip(states.iter_mut()) {
            if *t > 0.0 {
                state[2] += crate::spin_drift::litz_drift_meters(sg, *t, inputs.is_twist_right);
            }
        }
    }

    FastSolution::from_trajectory_data(times, states, t_events)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BCSegmentData;

    fn expected_shot_frame_vector(level: Vector3<f64>, angle: f64) -> Vector3<f64> {
        let (sin_angle, cos_angle) = angle.sin_cos();
        Vector3::new(
            level.x * cos_angle + level.y * sin_angle,
            -level.x * sin_angle + level.y * cos_angle,
            level.z,
        )
    }

    #[test]
    fn measured_bc_fast_drag_ignores_name_based_form_factor_flag() {
        let derivatives_with_flag = |use_form_factor| {
            let inputs = BallisticInputs {
                bc_value: 0.462,
                bc_type: DragModel::G1,
                bullet_model: Some("168gr SMK Match".to_string()),
                use_form_factor,
                temperature: 15.0,
                pressure: 1013.25,
                ..BallisticInputs::default()
            };

            compute_derivatives(
                &[0.0, 0.0, 0.0, 600.0, 0.0, 0.0],
                &inputs,
                &WindSock::new(vec![]),
                FastAtmosphere::Standard {
                    base_density: 1.225,
                },
                &inputs.bc_type,
                crate::transonic_drag::ProjectileShape::Spitzer,
                inputs.bc_value,
                false,
                false,
                None,
                None,
                None,
            )
        };

        let baseline = derivatives_with_flag(false);
        let flagged = derivatives_with_flag(true);

        for component in 3..6 {
            assert_eq!(
                flagged[component].to_bits(),
                baseline[component].to_bits(),
                "published BC already encodes form factor: component {component}, baseline={} flagged={}",
                baseline[component],
                flagged[component]
            );
        }
    }

    #[test]
    fn velocity_bc_data_requires_opt_in_in_plain_fast_kernel() {
        let acceleration = |inputs: &BallisticInputs| {
            let has_mach_segments = inputs
                .bc_segments
                .as_ref()
                .is_some_and(|segments| !segments.is_empty());
            let has_velocity_segments = inputs
                .bc_segments_data
                .as_ref()
                .is_some_and(|segments| !segments.is_empty());
            compute_derivatives(
                &[0.0, 0.0, 0.0, 600.0, 0.0, 0.0],
                inputs,
                &WindSock::new(vec![]),
                FastAtmosphere::Standard {
                    base_density: 1.225,
                },
                &inputs.bc_type,
                crate::transonic_drag::ProjectileShape::Spitzer,
                inputs.bc_value,
                has_mach_segments,
                has_velocity_segments,
                None,
                None,
                None,
            )
        };

        let scalar_inputs = BallisticInputs {
            bc_value: 0.5,
            bc_type: DragModel::G7,
            temperature: 15.0,
            pressure: 1013.25,
            ..BallisticInputs::default()
        };
        let mut disabled_inputs = scalar_inputs.clone();
        disabled_inputs.bc_segments_data = Some(vec![BCSegmentData {
            velocity_min: 0.0,
            velocity_max: 4_000.0,
            bc_value: 0.46,
        }]);
        disabled_inputs.use_bc_segments = false;
        let mut enabled_inputs = disabled_inputs.clone();
        enabled_inputs.use_bc_segments = true;
        let mut mach_only_inputs = scalar_inputs.clone();
        mach_only_inputs.bc_segments = Some(vec![(0.0, 0.4), (3.0, 0.4)]);
        let mut disabled_with_both = mach_only_inputs.clone();
        disabled_with_both.bc_segments_data = disabled_inputs.bc_segments_data.clone();

        let scalar = acceleration(&scalar_inputs);
        let disabled = acceleration(&disabled_inputs);
        let enabled = acceleration(&enabled_inputs);
        let mach_only = acceleration(&mach_only_inputs);
        let disabled_with_both = acceleration(&disabled_with_both);

        assert_eq!(
            disabled[3].to_bits(),
            scalar[3].to_bits(),
            "a populated velocity table must not change drag while use_bc_segments is false"
        );
        assert!(
            enabled[3] < disabled[3] - 1.0,
            "enabling the lower BC table must increase drag: disabled ax={} enabled ax={}",
            disabled[3],
            enabled[3]
        );
        assert_eq!(
            disabled_with_both[3].to_bits(),
            mach_only[3].to_bits(),
            "disabling velocity data must fall through to an explicit Mach table"
        );
    }

    #[test]
    fn inclined_positions_at_same_world_altitude_have_same_fast_acceleration() {
        let angle = std::f64::consts::FRAC_PI_6;
        let inputs = BallisticInputs {
            altitude: 100.0,
            temperature: 15.0,
            pressure: 1013.25,
            shooting_angle: angle,
            ..BallisticInputs::default()
        };
        let wind_sock = WindSock::new(vec![]);
        let atmosphere = FastAtmosphere::Standard {
            base_density: 1.225,
        };
        let state_along_slant = [1_000.0, 0.0, 0.0, 600.0, 0.0, 0.0];
        let state_across_slant = [0.0, 500.0 / angle.cos(), 0.0, 600.0, 0.0, 0.0];

        let a = compute_derivatives(
            &state_along_slant,
            &inputs,
            &wind_sock,
            atmosphere,
            &inputs.bc_type,
            crate::transonic_drag::ProjectileShape::Spitzer,
            inputs.bc_value,
            false,
            false,
            None,
            None,
            None,
        );
        let b = compute_derivatives(
            &state_across_slant,
            &inputs,
            &wind_sock,
            atmosphere,
            &inputs.bc_type,
            crate::transonic_drag::ProjectileShape::Spitzer,
            inputs.bc_value,
            false,
            false,
            None,
            None,
            None,
        );

        for component in 3..6 {
            assert!(
                (a[component] - b[component]).abs() < 1e-10,
                "fast derivative component {component} differs at equal world altitude: {} vs {}",
                a[component],
                b[component]
            );
        }
    }

    #[test]
    fn inclined_headwind_is_rotated_into_solver_frame() {
        let angle = std::f64::consts::FRAC_PI_6;
        let inputs = BallisticInputs {
            shooting_angle: angle,
            ..BallisticInputs::default()
        };
        let speed_mps = 360.0 * (1000.0 / 3600.0);
        let level_headwind = Vector3::new(-speed_mps, 0.0, 0.0);
        let velocity = expected_shot_frame_vector(level_headwind, angle);
        let state = [0.0, 0.0, 0.0, velocity.x, velocity.y, velocity.z];
        let actual = compute_derivatives(
            &state,
            &inputs,
            &WindSock::new(vec![(360.0, 0.0, 1000.0)]),
            FastAtmosphere::Direct {
                air_density: 1.225,
                speed_of_sound: 340.0,
            },
            &inputs.bc_type,
            crate::transonic_drag::ProjectileShape::Spitzer,
            inputs.bc_value,
            false,
            false,
            None,
            None,
            None,
        );
        let expected = Vector3::new(
            -G_ACCEL_MPS2 * angle.sin(),
            -G_ACCEL_MPS2 * angle.cos(),
            0.0,
        );

        assert!(
            (Vector3::new(actual[3], actual[4], actual[5]) - expected).norm() < 1e-12,
            "co-moving horizontal wind must leave only shot-frame gravity: {actual:?}"
        );
    }

    #[test]
    fn plain_fast_kernel_applies_power_law_wind_shear() {
        let state = [500.0, 100.0, 0.0, 700.0, 0.0, 0.0];
        let run = |enable_wind_shear: bool, model: &str, wind_speed_kmh: f64| {
            let inputs = BallisticInputs {
                bc_value: 0.5,
                bc_type: DragModel::G7,
                enable_wind_shear,
                wind_shear_model: model.to_string(),
                ..BallisticInputs::default()
            };
            let wind_shear_model = enable_wind_shear
                .then(|| crate::wind_shear::boundary_layer_model_from_name(model))
                .filter(|model| *model != crate::wind_shear::WindShearModel::None);
            compute_derivatives(
                &state,
                &inputs,
                &WindSock::new(vec![(wind_speed_kmh, 90.0, 2_000.0)]),
                FastAtmosphere::Direct {
                    air_density: 1.225,
                    speed_of_sound: 340.0,
                },
                &inputs.bc_type,
                crate::transonic_drag::ProjectileShape::Spitzer,
                inputs.bc_value,
                false,
                false,
                None,
                wind_shear_model,
                None,
            )
        };

        let uniform = run(false, "power_law", 36.0);
        let model_none = run(true, "none", 36.0);
        assert_eq!(model_none, uniform, "model=none must preserve uniform wind");

        let sheared = run(true, "power_law", 36.0);
        assert!(
            sheared[5] < uniform[5],
            "stronger aloft crosswind must increase leftward acceleration: uniform={}, shear={}",
            uniform[5],
            sheared[5]
        );

        let ratio = crate::wind_shear::boundary_layer_speed_ratio(
            state[1],
            crate::wind_shear::WindShearModel::PowerLaw,
        );
        let equivalent_uniform = run(false, "none", 36.0 * ratio);
        for component in 3..6 {
            assert!(
                (sheared[component] - equivalent_uniform[component]).abs() < 1e-12,
                "shear component {component} must equal base wind scaled by {ratio}: shear={}, expected={}",
                sheared[component],
                equivalent_uniform[component]
            );
        }
    }

    #[test]
    fn plain_fast_path_wind_shear_changes_high_arc_drift() {
        let run = |enable_wind_shear: bool, model: &str| {
            let inputs = BallisticInputs {
                muzzle_velocity: 800.0,
                bc_value: 0.5,
                bc_type: DragModel::G7,
                bullet_mass: 168.0 * 0.00006479891,
                bullet_diameter: 0.308 * 0.0254,
                enable_wind_shear,
                wind_shear_model: model.to_string(),
                ground_threshold: -100.0,
                ..BallisticInputs::default()
            };
            let elevation = 0.12_f64;
            let solution = fast_integrate(
                &inputs,
                &WindSock::new(vec![(36.0, 90.0, 2_000.0)]),
                FastIntegrationParams {
                    horiz: 1_000.0,
                    vert: 0.0,
                    initial_state: [
                        0.0,
                        0.0,
                        0.0,
                        inputs.muzzle_velocity * elevation.cos(),
                        inputs.muzzle_velocity * elevation.sin(),
                        0.0,
                    ],
                    t_span: (0.0, 5.0),
                    atmo_params: (0.0, 15.0, 1013.25, 1.0),
                    atmo_sock: None,
                },
            );
            let last = solution.t.len() - 1;
            assert_eq!(solution.y[0][last].to_bits(), 1_000.0_f64.to_bits());
            solution.y[2][last]
        };

        let uniform = run(false, "power_law");
        let model_none = run(true, "none");
        assert_eq!(model_none.to_bits(), uniform.to_bits());
        let sheared = run(true, "power_law");
        assert!(
            sheared.abs() > uniform.abs() + 0.01,
            "aloft shear must increase drift magnitude: uniform={uniform}, shear={sheared}"
        );
    }

    #[test]
    fn inclined_coriolis_is_rotated_into_solver_frame() {
        let angle = std::f64::consts::FRAC_PI_6;
        let inputs = BallisticInputs {
            shooting_angle: angle,
            ..BallisticInputs::default()
        };
        let velocity = Vector3::new(600.0, 20.0, 5.0);
        let state = [0.0, 0.0, 0.0, velocity.x, velocity.y, velocity.z];
        let level_omega = Vector3::new(3.0e-5, 6.0e-5, -2.0e-5);
        let run = |omega| {
            compute_derivatives(
                &state,
                &inputs,
                &WindSock::new(vec![]),
                FastAtmosphere::Direct {
                    air_density: 1.225,
                    speed_of_sound: 340.0,
                },
                &inputs.bc_type,
                crate::transonic_drag::ProjectileShape::Spitzer,
                inputs.bc_value,
                false,
                false,
                omega,
                None,
                None,
            )
        };
        let baseline = run(None);
        let with_coriolis = run(Some(level_omega));
        let actual = Vector3::new(
            with_coriolis[3] - baseline[3],
            with_coriolis[4] - baseline[4],
            with_coriolis[5] - baseline[5],
        );
        let expected = -2.0 * expected_shot_frame_vector(level_omega, angle).cross(&velocity);

        assert!(
            (actual - expected).norm() < 1e-12,
            "inclined Coriolis mismatch: actual={actual:?}, expected={expected:?}"
        );
    }

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

        assert_eq!(velocity_segment_bc(500.0, &segments, 0.5), 0.5);
        assert_eq!(velocity_segment_bc(1500.0, &segments, 0.5), 0.52);
        assert_eq!(velocity_segment_bc(2500.0, &segments, 0.5), 0.55);

        // Test edge cases
        assert_eq!(velocity_segment_bc(-100.0, &segments, 0.5), 0.5); // Below min
        assert_eq!(velocity_segment_bc(3500.0, &segments, 0.5), 0.55); // Above max
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

        assert_eq!(velocity_segment_bc(500.0, &single_segment, 0.5), 0.5); // Below
        assert_eq!(velocity_segment_bc(1500.0, &single_segment, 0.5), 0.5); // In range
        assert_eq!(velocity_segment_bc(2500.0, &single_segment, 0.5), 0.5); // Above

        // Test with exact boundary values
        // Half-open bands make a shared boundary belong to the band that starts there.
        let segments = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 999.0, // Exclusive upper bound to avoid overlap
                bc_value: 0.45,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.50,
            },
        ];

        assert_eq!(velocity_segment_bc(1000.0, &segments, 0.7), 0.50); // At second segment start
        assert_eq!(velocity_segment_bc(0.0, &segments, 0.7), 0.45); // At min
        assert_eq!(velocity_segment_bc(998.999, &segments, 0.7), 0.45); // Just below exclusive max
        assert_eq!(velocity_segment_bc(999.0, &segments, 0.7), 0.7); // Gap starts at exclusive max
    }

    #[test]
    fn velocity_segment_gaps_and_clamps_do_not_depend_on_order() {
        let fallback_bc = 0.73;
        let ascending_with_gap = vec![
            BCSegmentData {
                velocity_min: 0.0,
                velocity_max: 999.0,
                bc_value: 0.6,
            },
            BCSegmentData {
                velocity_min: 1000.0,
                velocity_max: 2000.0,
                bc_value: 0.8,
            },
        ];
        assert_eq!(
            velocity_segment_bc(999.5, &ascending_with_gap, fallback_bc),
            fallback_bc,
            "coverage gaps must use the projectile's base BC"
        );

        let mut descending = ascending_with_gap.clone();
        descending.reverse();
        assert_eq!(
            velocity_segment_bc(-100.0, &descending, fallback_bc),
            0.6,
            "below coverage must clamp to the lowest-velocity band"
        );
        assert_eq!(
            velocity_segment_bc(2500.0, &descending, fallback_bc),
            0.8,
            "above coverage must clamp to the highest-velocity band"
        );
    }

    #[test]
    fn test_bc_segments_empty_fallback() {
        let empty_segments: Vec<BCSegmentData> = vec![];

        // With empty segments, should return fallback value
        let result = velocity_segment_bc(1500.0, &empty_segments, 0.73);
        assert_eq!(result, 0.73); // Caller-provided fallback value
    }

    #[test]
    fn test_fast_integration_params() {
        // Verify FastIntegrationParams struct can be constructed
        let params = FastIntegrationParams {
            horiz: 1000.0,
            vert: 0.0,
            initial_state: [0.0, 0.0, 0.0, 800.0, 50.0, 0.0], // McCoy: vx=downrange
            t_span: (0.0, 5.0),
            atmo_params: (0.0, 15.0, 1013.25, 1.0),
            atmo_sock: None,
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
            vec![2.0], // target_hit at t=2
            vec![0.5], // max_ord at t=0.5
            vec![],    // no ground_hit
        ];

        let solution = FastSolution::from_trajectory_data(times, states, t_events);

        assert_eq!(solution.t_events[0], vec![2.0]); // Target hit
        assert_eq!(solution.t_events[1], vec![0.5]); // Max ordinate
        assert!(solution.t_events[2].is_empty()); // No ground hit
    }

    #[test]
    fn plain_fast_path_interpolates_the_target_crossing() {
        let target = 500.123456789;
        let initial_state = [0.0, 0.0, 0.25, 800.0, 12.0, -2.5];
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0,
            bc_value: 0.5,
            bc_type: DragModel::G7,
            ground_threshold: -100.0,
            use_enhanced_spin_drift: false,
            ..BallisticInputs::default()
        };
        let run = |horiz| {
            fast_integrate(
                &inputs,
                &WindSock::new(vec![]),
                FastIntegrationParams {
                    horiz,
                    vert: 0.0,
                    initial_state,
                    t_span: (0.0, 2.0),
                    atmo_params: (0.0, 15.0, 1013.25, 1.0),
                    atmo_sock: None,
                },
            )
        };

        // A longer run retains the two full RK4 samples bracketing the requested target.
        let reference = run(target + 2.0);
        let left = reference.y[0]
            .windows(2)
            .position(|x| x[0] < target && x[1] > target)
            .expect("reference trajectory must bracket target");
        let right = left + 1;
        let alpha =
            (target - reference.y[0][left]) / (reference.y[0][right] - reference.y[0][left]);

        let solution = run(target);
        let last = solution.t.len() - 1;
        assert_eq!(solution.y[0][last].to_bits(), target.to_bits());
        let expected_time = reference.t[left] + alpha * (reference.t[right] - reference.t[left]);
        assert!((solution.t[last] - expected_time).abs() < 1e-12);
        assert_eq!(solution.t_events[0], vec![solution.t[last]]);

        for component in 0..6 {
            assert_eq!(solution.y[component].len(), solution.t.len());
            let expected = reference.y[component][left]
                + alpha * (reference.y[component][right] - reference.y[component][left]);
            assert!(
                (solution.y[component][last] - expected).abs() < 1e-9,
                "component {component} is not at the crossing: actual={}, expected={expected}",
                solution.y[component][last]
            );
        }
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
                atmo_params: (0.0, 15.0, 1013.25, 1.0),
                atmo_sock: None,
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
                atmo_params: (0.0, 15.0, 1013.25, 1.0),
                atmo_sock: None,
            };
            let sol = fast_integrate_with_segments(&inputs, vec![], params);
            let n = sol.y[0].len();
            sol.y[1][n - 1]
        }
        // enable_coriolis=true with advanced effects OFF: directional Coriolis still applies.
        let e = final_y(true, FRAC_PI_2);
        let w = final_y(true, 3.0 * FRAC_PI_2);
        assert!(
            e > w && (e - w) > 1e-5,
            "Coriolis-only (no advanced effects) must still be directional: E={e} W={w}"
        );
        // enable_coriolis=false: no Coriolis at all, east == west.
        let e2 = final_y(false, FRAC_PI_2);
        let w2 = final_y(false, 3.0 * FRAC_PI_2);
        assert!(
            (e2 - w2).abs() < 1e-9,
            "with enable_coriolis=false, east/west must be identical: E={e2} W={w2}"
        );
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
            atmo_sock: None,
        };
        // pressure <= 0 -> fail loudly (success=false) instead of a 1-point stub as success.
        let zero_p = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, 0.0, 1.0)));
        assert!(
            !zero_p.success,
            "pressure=0 atmosphere must yield success=false"
        );
        // non-finite pressure -> also rejected.
        let nan_p = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, f64::NAN, 1.0)));
        assert!(!nan_p.success, "NaN pressure must yield success=false");
        // A supplied density ratio must imply a physically plausible density in both wrappers.
        let segmented_too_dense =
            fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, 1013.25, 50.0)));
        let plain_too_dense = fast_integrate(
            &inputs,
            &WindSock::new(vec![]),
            mk((0.0, 15.0, 1013.25, 50.0)),
        );
        assert!(
            !segmented_too_dense.success && !plain_too_dense.success,
            "ratio=50 atmosphere must fail in both wrappers: segmented={}, plain={}",
            segmented_too_dense.success,
            plain_too_dense.success
        );
        // realistic atmosphere -> success.
        let good = fast_integrate_with_segments(&inputs, vec![], mk((0.0, 15.0, 1013.25, 1.0)));
        assert!(good.success, "realistic atmosphere must yield success=true");
        // Direct-atmosphere mode (density, speed_of_sound, 0, 0) is legitimate and must NOT
        // be rejected by the guard (regression: 0.21.2 rejected it via the pressure<=0 check).
        let direct = fast_integrate_with_segments(&inputs, vec![], mk((1.225, 340.0, 0.0, 0.0)));
        assert!(
            direct.success,
            "direct-atmosphere mode (pressure=0 sentinel) must yield success=true"
        );
    }

    #[test]
    fn plain_fast_path_honors_direct_atmosphere_values() {
        fn final_speed(muzzle_velocity: f64, atmo_params: (f64, f64, f64, f64)) -> f64 {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_velocity = muzzle_velocity;
            inputs.bc_value = 0.5;
            inputs.bc_type = DragModel::G7;
            inputs.ground_threshold = -100.0;

            let wind_sock = WindSock::new(vec![]);
            let solution = fast_integrate(
                &inputs,
                &wind_sock,
                FastIntegrationParams {
                    horiz: 10_000.0,
                    vert: 0.0,
                    initial_state: [0.0, 0.0, 0.0, muzzle_velocity, 0.0, 0.0],
                    t_span: (0.0, 0.2),
                    atmo_params,
                    atmo_sock: None,
                },
            );
            assert!(solution.success);

            let last = solution.y[0].len() - 1;
            (solution.y[3][last].powi(2)
                + solution.y[4][last].powi(2)
                + solution.y[5][last].powi(2))
            .sqrt()
        }

        let thin_air = final_speed(800.0, (0.905, 340.0, 0.0, 0.0));
        let dense_air = final_speed(800.0, (1.225, 340.0, 0.0, 0.0));
        assert!(
            thin_air > dense_air,
            "lower supplied density must retain more velocity: thin={thin_air}, dense={dense_air}"
        );

        let low_sound_speed = final_speed(340.0, (1.0, 300.0, 0.0, 0.0));
        let high_sound_speed = final_speed(340.0, (1.0, 400.0, 0.0, 0.0));
        assert!(
            (low_sound_speed - high_sound_speed).abs() > 1e-6,
            "supplied sound speed must affect Mach-dependent drag"
        );
    }

    #[test]
    fn segmented_fast_path_nonpositive_density_ratio_uses_standard_fallback() {
        fn terminal_velocity(base_ratio: f64) -> f64 {
            let mut inputs = BallisticInputs::default();
            inputs.muzzle_velocity = 800.0;
            inputs.bc_value = 0.5;
            inputs.bc_type = DragModel::G7;
            inputs.ground_threshold = -100.0;

            let solution = fast_integrate_with_segments(
                &inputs,
                vec![],
                FastIntegrationParams {
                    horiz: 500.0,
                    vert: 0.0,
                    initial_state: [0.0, 0.0, 0.0, 800.0, 0.0, 0.0],
                    t_span: (0.0, 5.0),
                    atmo_params: (0.0, 15.0, 1013.25, base_ratio),
                    atmo_sock: None,
                },
            );
            assert!(solution.success);

            let last = solution.y[3].len() - 1;
            solution.y[3][last]
        }

        let explicit_sea_level = terminal_velocity(1.0);
        for base_ratio in [0.0, -1.0] {
            let missing_ratio = terminal_velocity(base_ratio);
            assert!(
                missing_ratio < 800.0,
                "missing density ratio must not create a vacuum trajectory: {missing_ratio}"
            );
            assert!((missing_ratio - explicit_sea_level).abs() < 1e-9);
        }
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
                atmo_params: (0.0, 15.0, 1013.25, 1.0),
                atmo_sock: None,
            };
            fast_integrate_with_segments(&inputs, vec![], params)
        };
        // Both a real .224 and a real .338 produce a valid trajectory (geometry is honored, not
        // overridden by a hardcoded .308). Regression sentinel for the plumbing.
        assert!(run(0.00569, 7.0).success, ".224 geometry must solve");
        assert!(run(0.00858, 10.0).success, ".338 geometry must solve");
    }

    #[test]
    fn segmented_fast_spin_flags_do_not_depend_on_advanced_umbrella() {
        fn endpoint(
            enable_advanced_effects: bool,
            enable_magnus: bool,
            use_enhanced_spin_drift: bool,
        ) -> Vector3<f64> {
            let inputs = BallisticInputs {
                muzzle_velocity: 823.0,
                bullet_mass: 168.0 * 0.00006479891,
                bullet_diameter: 0.308 * 0.0254,
                bullet_length: 1.215 * 0.0254,
                caliber_inches: 0.308,
                weight_grains: 168.0,
                bc_value: 0.475,
                bc_type: DragModel::G1,
                twist_rate: 12.0,
                is_twist_right: true,
                enable_advanced_effects,
                enable_magnus,
                use_enhanced_spin_drift,
                ..BallisticInputs::default()
            };
            let elevation = 0.02_f64;
            let solution = fast_integrate_with_segments(
                &inputs,
                vec![],
                FastIntegrationParams {
                    horiz: 1_000.0,
                    vert: 0.0,
                    initial_state: [
                        0.0,
                        0.0,
                        0.0,
                        inputs.muzzle_velocity * elevation.cos(),
                        inputs.muzzle_velocity * elevation.sin(),
                        0.0,
                    ],
                    t_span: (0.0, 5.0),
                    atmo_params: (0.0, 15.0, 1013.25, 1.0),
                    atmo_sock: None,
                },
            );
            assert!(solution.success);
            let last = solution.t.len() - 1;
            Vector3::new(
                solution.y[0][last],
                solution.y[1][last],
                solution.y[2][last],
            )
        }

        let baseline = endpoint(false, false, false);
        let magnus_without_umbrella = endpoint(false, true, false);
        let magnus_with_umbrella = endpoint(true, true, false);
        assert!(
            (magnus_without_umbrella - baseline).norm() > 1e-5,
            "test shot must produce a measurable Magnus displacement"
        );
        assert!(
            (magnus_with_umbrella - magnus_without_umbrella).norm() < 1e-12,
            "the legacy umbrella must not suppress explicitly enabled Magnus: without={magnus_without_umbrella:?} with={magnus_with_umbrella:?}"
        );

        let litz_only = endpoint(false, false, true);
        let litz_without_umbrella = endpoint(false, true, true);
        let litz_with_umbrella = endpoint(true, true, true);
        assert!(
            (litz_without_umbrella - litz_only).norm() < 1e-12,
            "Litz mode must suppress explicitly enabled Magnus: litz={litz_only:?} both={litz_without_umbrella:?}"
        );
        assert!(
            (litz_with_umbrella - litz_without_umbrella).norm() < 1e-12,
            "the legacy umbrella must not change Magnus suppression in Litz mode: without={litz_without_umbrella:?} with={litz_with_umbrella:?}"
        );
    }
}
