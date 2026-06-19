use crate::atmosphere::{get_direct_atmosphere, get_local_atmosphere};
use crate::bc_estimation::BCSegmentEstimator;
use crate::constants::*;
use crate::drag::get_drag_coefficient_full;
use crate::form_factor::apply_form_factor_to_drag;
use crate::spin_drift::{apply_enhanced_spin_drift, calculate_enhanced_spin_drift};
use crate::InternalBallisticInputs as BallisticInputs;
use nalgebra::Vector3;

// Physics constants
const INCHES_PER_FOOT: f64 = 12.0;
const STANDARD_AIR_DENSITY_METRIC: f64 = 1.225; // kg/m³ at sea level

// Magnus Effect Constants
//
// The Magnus effect causes spinning projectiles to deflect perpendicular to both
// their velocity vector and spin axis due to asymmetric pressure distribution.
// These constants define the Magnus moment coefficient (C_Lα) for different flight regimes.

/// Magnus coefficient for subsonic flow (M < 0.8)
///
/// Value: 0.030 (dimensionless coefficient)
/// Physical basis: Fully developed boundary layer circulation around spinning projectile
/// Regime: Subsonic flow where boundary layer remains attached
/// Source: McCoy's "Modern Exterior Ballistics", validated against wind tunnel data
const MAGNUS_COEFF_SUBSONIC: f64 = 0.030;

/// Magnus coefficient reduction factor for transonic regime (0.8 < M < 1.2)
///
/// Value: 0.0075 (25% of subsonic value at M=1.2)
/// Physical basis: Shock waves disrupt circulation patterns, reducing Magnus effect
/// Effect: Spin drift significantly reduced in transonic flight
/// Source: Experimental spinning projectile studies
const MAGNUS_COEFF_TRANSONIC_REDUCTION: f64 = 0.0075;

/// Base Magnus coefficient for supersonic flow (M > 1.2)
///
/// Value: 0.015 (dimensionless coefficient)
/// Physical basis: Shock-dominated flow with reduced but persistent circulation
/// Effect: Lower Magnus effect than subsonic, but higher than transonic minimum
const MAGNUS_COEFF_SUPERSONIC_BASE: f64 = 0.015;

/// Magnus coefficient scaling factor for high supersonic speeds
///
/// Value: 0.0044 (additional scaling with Mach number)
/// Formula: Magnus_coeff = BASE + SCALE * (M - 1.2) for M > 1.2
/// Physical basis: Partial recovery of circulation effects at higher Mach numbers
const MAGNUS_COEFF_SUPERSONIC_SCALE: f64 = 0.0044;

/// Transonic regime boundaries for Magnus effect calculations
const MAGNUS_TRANSONIC_LOWER: f64 = 0.8; // Lower bound of transonic regime
const MAGNUS_TRANSONIC_UPPER: f64 = 1.2; // Upper bound of transonic regime
const MAGNUS_TRANSONIC_RANGE: f64 = 0.4; // Range width (1.2 - 0.8)
const MAGNUS_SUPERSONIC_RANGE: f64 = 1.8; // Scaling range for supersonic recovery

// Note: These Magnus coefficients are calibrated against real-world spin drift measurements
// from McCoy's "Modern Exterior Ballistics" and experimental data. The dimensionless
// coefficients represent the Magnus moment per unit angle of attack.

// Atmosphere detection thresholds
const MAX_REALISTIC_DENSITY: f64 = 2.0; // kg/m³
const MIN_REALISTIC_SPEED_OF_SOUND: f64 = 200.0; // m/s

/// Calculate spin rate from twist rate and velocity
fn calculate_spin_rate(twist_rate: f64, velocity_mps: f64) -> f64 {
    if twist_rate <= 0.0 {
        return 0.0;
    }

    // Convert velocity to ft/s and twist rate to ft/turn
    let velocity_fps = velocity_mps * MPS_TO_FPS;
    let twist_rate_ft = twist_rate / INCHES_PER_FOOT;

    // Calculate spin rate: revolutions per second = velocity_fps / twist_rate_ft
    // Convert to rad/s: rad/s = (revolutions/s) * 2π
    let revolutions_per_second = velocity_fps / twist_rate_ft;

    revolutions_per_second * 2.0 * std::f64::consts::PI
}

/// Calculate Magnus moment coefficient C_Lα based on Mach number
/// Based on McCoy's 'Modern Exterior Ballistics' and empirical data
pub(crate) fn calculate_magnus_moment_coefficient(mach: f64) -> f64 {
    // Magnus moment coefficient varies with Mach number
    // Values based on empirical data for spitzer bullets

    if mach < MAGNUS_TRANSONIC_LOWER {
        // Subsonic: relatively constant
        MAGNUS_COEFF_SUBSONIC
    } else if mach < MAGNUS_TRANSONIC_UPPER {
        // Transonic: reduced due to shock formation
        // Linear interpolation through transonic region
        MAGNUS_COEFF_SUBSONIC
            - MAGNUS_COEFF_TRANSONIC_REDUCTION * (mach - MAGNUS_TRANSONIC_LOWER)
                / MAGNUS_TRANSONIC_RANGE
    } else {
        // Supersonic: gradually recovers
        MAGNUS_COEFF_SUPERSONIC_BASE
            + MAGNUS_COEFF_SUPERSONIC_SCALE
                * ((mach - MAGNUS_TRANSONIC_UPPER) / MAGNUS_SUPERSONIC_RANGE).min(1.0)
    }
}

/// Compute ballistic derivatives for trajectory integration
pub fn compute_derivatives(
    pos: Vector3<f64>,
    vel: Vector3<f64>,
    inputs: &BallisticInputs,
    wind_vector: Vector3<f64>,
    atmos_params: (f64, f64, f64, f64),
    bc_used: f64,
    omega_vector: Option<Vector3<f64>>,
    time: f64,
) -> [f64; 6] {
    // Gravity acceleration vector
    let accel_gravity = Vector3::new(0.0, -G_ACCEL_MPS2, 0.0);

    // Wind-adjusted velocity
    let velocity_adjusted = vel - wind_vector;
    let speed_air = velocity_adjusted.norm();

    // Initialize drag acceleration
    let mut accel_drag = Vector3::zeros();
    let mut accel_magnus = Vector3::zeros();

    // Calculate drag if velocity is significant
    if speed_air > crate::constants::MIN_VELOCITY_THRESHOLD {
        let v_rel_fps = speed_air * MPS_TO_FPS;

        // Get atmospheric conditions
        let altitude_at_pos = inputs.altitude + pos[1];

        // Check if we have direct atmosphere values
        // Direct atmosphere is indicated by having only 2 parameters where:
        // params[0] = air density, params[1] = speed of sound
        // params[2] and params[3] would be 0.0
        // BUT: we need to check if params[0] is a reasonable density value (< 2.0 kg/m³)
        let (air_density, speed_of_sound, temperature_c) = if atmos_params.0 < MAX_REALISTIC_DENSITY
            && atmos_params.1 > MIN_REALISTIC_SPEED_OF_SOUND
            && atmos_params.2 == 0.0
            && atmos_params.3 == 0.0
        {
            // Direct atmosphere values: atmos_params.1 is the SPEED OF SOUND here, NOT Celsius,
            // so back-compute temperature from it (c = sqrt(1.4*287.05*T_k)) for the Reynolds
            // correction below — which previously read atmos_params.1 as temperature directly.
            let (rho, sound) = get_direct_atmosphere(atmos_params.0, atmos_params.1);
            (rho, sound, sound * sound / (1.4 * 287.05) - 273.15)
        } else {
            // Calculate from base parameters
            let (rho, sound) = get_local_atmosphere(
                altitude_at_pos,
                atmos_params.0, // base_alt
                atmos_params.1, // base_temp_c
                atmos_params.2, // base_press_hpa
                atmos_params.3, // base_ratio
            );
            // LOCAL temperature at the projectile altitude, back-computed from the LOCAL speed of
            // sound (get_local_atmosphere returns density/sound at altitude_at_pos but not temp;
            // its sound = sqrt(1.4*287.05*T_k)). Using base_temp_c here would feed the Reynolds
            // viscosity the shooter-altitude temperature while density/sound are local.
            (rho, sound, sound * sound / (1.4 * 287.05) - 273.15)
        };

        // Calculate Mach number with safe division
        let mach = if speed_of_sound > 1e-9 {
            speed_air / speed_of_sound
        } else {
            0.0 // No meaningful Mach number at zero speed of sound
        };

        // Get drag coefficient with transonic and Reynolds corrections
        let mut drag_factor = get_drag_coefficient_full(
            mach,
            &inputs.bc_type,
            false, // transonic applied exactly once below (was double-applied here + in block)
            false, // Reynolds applied once below (manual block ~243); was double-applied here + there
            None, // let it determine shape
            if inputs.caliber_inches > 0.0 {
                Some(inputs.caliber_inches)
            } else {
                Some(inputs.bullet_diameter / 0.0254) // meters -> inches
            },
            if inputs.weight_grains > 0.0 {
                Some(inputs.weight_grains)
            } else {
                Some(inputs.bullet_mass / 0.00006479891) // kg -> grains
            },
            Some(speed_air),
            Some(air_density),
            Some(atmos_params.1), // temperature in Celsius
        );

        // Apply form factor if enabled
        if inputs.use_form_factor {
            drag_factor = apply_form_factor_to_drag(
                drag_factor,
                inputs.bullet_model.as_deref(),
                &inputs.bc_type,
                true,
            );
        }

        // Get BC value
        let mut bc_val = bc_used;

        if inputs.use_bc_segments {
            // First try velocity-based segments if available
            if inputs.bc_segments_data.is_some() {
                bc_val = get_bc_for_velocity(v_rel_fps, inputs, bc_used);
            } else if let Some(ref segments) = inputs.bc_segments {
                // Fall back to Mach-based segments when use_bc_segments=true but no velocity data
                bc_val = interpolated_bc(mach, segments, Some(inputs));
            } else {
                // No explicit segments - try BC estimation
                bc_val = get_bc_for_velocity(v_rel_fps, inputs, bc_used);
            }
        } else if let Some(ref segments) = inputs.bc_segments {
            // Explicit Mach-based segments (legacy behavior when use_bc_segments=false)
            bc_val = interpolated_bc(mach, segments, Some(inputs));
        }

        // Guard bc_val == 0 (allowed on the FFI/WASM/library surfaces, which lack the CLI's
        // 0.001 floor, and a user-supplied BC segment can be 0): the drag division below would be
        // Inf -> NaN, poisoning the whole trajectory. Mirrors the guards already in
        // cli_api::calculate_acceleration and fast_trajectory::compute_derivatives. Inert for
        // valid BCs (>= 0.001).
        let bc_val = bc_val.max(1e-6);

        // Calculate yaw effect with safe division
        let yaw_deg = if inputs.tipoff_decay_distance.abs() > 1e-9 {
            inputs.tipoff_yaw * (-pos[0] / inputs.tipoff_decay_distance).exp()
        } else {
            inputs.tipoff_yaw // No decay if distance is zero
        };
        let yaw_rad = yaw_deg.to_radians();
        let yaw_multiplier = 1.0 + yaw_rad.powi(2);

        // Calculate density scaling
        let density_scale = air_density / STANDARD_AIR_DENSITY;

        // Apply the transonic drag-rise correction exactly ONCE. The base Cd above is taken
        // WITHOUT transonic correction (apply_transonic_correction=false), so this is the only
        // application. Previously the correction was applied here AND inside
        // get_drag_coefficient_full, which squared the drag-rise factor and double-counted wave
        // drag across the transonic band (Cd ~3x too high near Mach 1). transonic_correction
        // self-gates via the projectile's critical Mach (returns the input unchanged outside the
        // band), and include_wave_drag=false matches cli_api::calculate_drag_coefficient — the
        // G1/G7 tables already embed the transonic rise, so additive wave drag would double-count.
        // Use the same SI fallbacks as the get_drag_coefficient_full call above (and
        // fast_trajectory): an SI-only caller may leave caliber_inches/weight_grains at 0, so
        // derive them from the SI bullet_diameter/bullet_mass rather than feeding zeros into
        // get_projectile_shape (which would mis-classify the shape via weight/caliber).
        let caliber_in = if inputs.caliber_inches > 0.0 {
            inputs.caliber_inches
        } else {
            inputs.bullet_diameter / 0.0254 // meters -> inches
        };
        let weight_gr = if inputs.weight_grains > 0.0 {
            inputs.weight_grains
        } else {
            inputs.bullet_mass / 0.00006479891 // kg -> grains
        };
        // MBA-949: shared resolver so named bullet_model shapes are honored here too (this path
        // previously used only the caliber/weight heuristic and ignored the name).
        let shape = crate::transonic_drag::resolve_projectile_shape(
            inputs.bullet_model.as_deref(),
            caliber_in,
            weight_gr,
            &inputs.bc_type.to_string(),
        );
        let drag_factor =
            crate::transonic_drag::transonic_correction(mach, drag_factor, shape, false);

        // Apply Reynolds correction for low velocities
        let drag_factor = if mach < 1.0 && speed_air < 200.0 {
            // temperature_c is derived per atmosphere mode above (base_temp_c, or back-computed
            // from the speed of sound in direct-atmosphere mode where atmos_params.1 is NOT Celsius).
            crate::reynolds::apply_reynolds_correction(
                drag_factor,
                speed_air,
                caliber_in, // inches, with SI fallback (shared with the transonic block above);
                // apply_reynolds_correction converts to meters internally. SI-only callers leave
                // caliber_inches at 0, which would otherwise feed 0 into the Reynolds calc.
                air_density,
                temperature_c,
                mach,
            )
        } else {
            drag_factor
        };

        // MBA-940: a user-supplied custom drag table overrides the G-model Cd entirely and is used
        // as-is — the transonic/Reynolds/form-factor corrections above are intentionally NOT
        // applied to it (the curve already encodes the projectile's true drag, so applying them
        // would distort/double-count it).
        let drag_factor = match inputs.custom_drag_table {
            Some(ref table) => table.interpolate(mach),
            None => drag_factor,
        };

        // Calculate drag acceleration
        let standard_factor = drag_factor * CD_TO_RETARD;
        let a_drag_ft_s2 =
            (v_rel_fps.powi(2) * standard_factor * yaw_multiplier * density_scale) / bc_val;
        let a_drag_m_s2 = a_drag_ft_s2 * FPS_TO_MPS;

        // Apply drag in opposite direction of relative velocity
        accel_drag = -a_drag_m_s2 * (velocity_adjusted / speed_air);

        // Magnus Effect calculation. Gated on enable_magnus specifically so it is
        // independent of Coriolis (matches the cli_api solver's decoupled flags).
        if inputs.enable_magnus && inputs.bullet_diameter > 0.0 && inputs.twist_rate > 0.0
        {
            // Calculate spin rate from twist rate and velocity
            let spin_rate_rad_s = calculate_spin_rate(inputs.twist_rate, speed_air);

            let c_np = calculate_magnus_moment_coefficient(mach);

            // bullet_diameter is SI (meters)
            let diameter_m = inputs.bullet_diameter;

            // Calculate spin parameter (dimensionless) with safe division
            let spin_param = if speed_air > 1e-9 {
                spin_rate_rad_s * diameter_m / (2.0 * speed_air)
            } else {
                0.0 // No spin effect at zero speed
            };

            // Calculate reference area
            let area = std::f64::consts::PI * (diameter_m / 2.0).powi(2);

            // Yaw of repose for the proper Magnus force. Stability/yaw helpers are
            // imperial: use the explicit imperial mirror fields, and convert the SI
            // bullet_length to inches at this boundary.
            let d_in = inputs.caliber_inches;
            let m_gr = inputs.weight_grains;
            let l_in = if inputs.bullet_length > 0.0 {
                inputs.bullet_length / 0.0254 // meters -> inches
            } else {
                4.5 * d_in.max(1e-9)
            };
            let sg = crate::spin_drift::miller_stability(d_in, m_gr, inputs.twist_rate, l_in);
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

            // Proper McCoy Magnus FORCE: F = q S C_Npa (pd/2V) sin(alpha_R).
            let magnus_force_magnitude =
                0.5 * air_density * speed_air.powi(2) * area * c_np * spin_param * yaw_rad.sin();

            // Magnus force is perpendicular to both velocity and spin axis
            // For a bullet spinning around its axis of travel, the spin vector is aligned with velocity
            let velocity_unit = velocity_adjusted / speed_air;

            // The Magnus force creates lift perpendicular to velocity
            // For right-hand twist, force is to the right when looking downrange
            // We need a vector perpendicular to velocity in the horizontal plane

            // Simplified approach: Magnus primarily causes horizontal drift
            // The force is perpendicular to both spin axis (velocity) and gravity
            let vertical = Vector3::new(0.0, 1.0, 0.0); // Up direction

            // Magnus force direction: velocity × vertical (for right-hand twist)
            let magnus_direction = velocity_unit.cross(&vertical);
            let magnus_norm = magnus_direction.norm();

            if magnus_norm > 1e-12 && magnus_force_magnitude > 1e-12 {
                let magnus_direction = magnus_direction / magnus_norm;

                // Reverse direction for left-hand twist
                let magnus_direction = if inputs.is_twist_right {
                    magnus_direction
                } else {
                    -magnus_direction
                };

                // Convert bullet mass to kg
                let bullet_mass_kg = inputs.bullet_mass; // already kg (SI)

                // Calculate acceleration
                accel_magnus = (magnus_force_magnitude / bullet_mass_kg) * magnus_direction;
            }
        }
    }

    // Total acceleration
    let mut accel = accel_gravity + accel_drag + accel_magnus;

    // Add Coriolis acceleration if omega vector is provided. The physical Coriolis term is
    // -2 Ω×v (MBA-957: the old +2 "frame-relabel" justification was wrong — it flipped the
    // lateral drift; the caller now builds omega with the corrected lateral sign, matching the
    // validated cli_api solver, so the canonical -2 applies directly).
    if let Some(omega) = omega_vector {
        let accel_coriolis = -2.0 * omega.cross(&vel);
        accel += accel_coriolis;
    }

    // Apply enhanced spin drift if enabled
    let mut derivatives = [vel[0], vel[1], vel[2], accel[0], accel[1], accel[2]];

    if inputs.use_enhanced_spin_drift && inputs.enable_advanced_effects && time > 0.0 {
        // Calculate crosswind component
        let velocity_adjusted = vel - wind_vector;
        let crosswind_speed = if velocity_adjusted.norm() > crate::constants::MIN_VELOCITY_THRESHOLD
        {
            let trajectory_unit = velocity_adjusted / velocity_adjusted.norm();
            let crosswind = wind_vector - wind_vector.dot(&trajectory_unit) * trajectory_unit;
            crosswind.norm()
        } else {
            0.0
        };

        // Get air density (already calculated above)
        let air_density = if speed_air > crate::constants::MIN_VELOCITY_THRESHOLD {
            let altitude_at_pos = inputs.altitude + pos[1];
            let (density, _) = if atmos_params.0 < MAX_REALISTIC_DENSITY
                && atmos_params.1 > MIN_REALISTIC_SPEED_OF_SOUND
                && atmos_params.2 == 0.0
                && atmos_params.3 == 0.0
            {
                get_direct_atmosphere(atmos_params.0, atmos_params.1)
            } else {
                get_local_atmosphere(
                    altitude_at_pos,
                    atmos_params.0,
                    atmos_params.1,
                    atmos_params.2,
                    atmos_params.3,
                )
            };
            density
        } else {
            STANDARD_AIR_DENSITY_METRIC // Standard air density
        };

        // Calculate enhanced spin drift components
        // calculate_enhanced_spin_drift is imperial (grains/inches): convert at boundary.
        let spin_components = calculate_enhanced_spin_drift(
            inputs.weight_grains,
            vel.norm(),
            inputs.twist_rate,
            inputs.caliber_inches,
            inputs.bullet_length / 0.0254, // meters -> inches
            inputs.is_twist_right,
            time,
            air_density,
            crosswind_speed,
            0.0,   // pitch_rate_rad_s - we don't track angular rates yet
            false, // use_pitch_damping - disabled for now
        );

        // Apply enhanced spin drift acceleration
        apply_enhanced_spin_drift(
            &mut derivatives,
            &spin_components,
            time,
            inputs.is_twist_right,
        );
    }

    // Return state derivatives: [velocity, acceleration]
    derivatives
}

/// Calculate appropriate BC fallback based on available bullet parameters
fn calculate_bc_fallback(
    bullet_mass: Option<f64>,     // grains
    bullet_diameter: Option<f64>, // inches
    bc_type: Option<&str>,        // "G1" or "G7"
) -> f64 {
    use crate::constants::*;

    // Weight-based fallback (most reliable predictor)
    if let Some(weight) = bullet_mass {
        let base_bc = if weight < 50.0 {
            BC_FALLBACK_ULTRA_LIGHT
        } else if weight < 100.0 {
            BC_FALLBACK_LIGHT
        } else if weight < 150.0 {
            BC_FALLBACK_MEDIUM
        } else if weight < 200.0 {
            BC_FALLBACK_HEAVY
        } else {
            BC_FALLBACK_VERY_HEAVY
        };

        // G7 vs G1 adjustment
        return if let Some(drag_model) = bc_type {
            if drag_model == "G7" {
                base_bc * 0.85 // G7 BCs are typically lower than G1
            } else {
                base_bc
            }
        } else {
            base_bc
        };
    }

    // Caliber-based fallback (second most reliable)
    if let Some(caliber) = bullet_diameter {
        let base_bc = if caliber <= 0.224 {
            BC_FALLBACK_SMALL_CALIBER
        } else if caliber <= 0.243 {
            BC_FALLBACK_MEDIUM_CALIBER
        } else if caliber <= 0.284 {
            BC_FALLBACK_LARGE_CALIBER
        } else {
            BC_FALLBACK_XLARGE_CALIBER
        };

        // G7 vs G1 adjustment
        return if let Some(drag_model) = bc_type {
            if drag_model == "G7" {
                base_bc * 0.85 // G7 BCs are typically lower than G1
            } else {
                base_bc
            }
        } else {
            base_bc
        };
    }

    // Final fallback - conservative overall
    let base_fallback = BC_FALLBACK_CONSERVATIVE;
    if let Some(drag_model) = bc_type {
        if drag_model == "G7" {
            return base_fallback * 0.85;
        }
    }

    base_fallback
}

/// Interpolate ballistic coefficient from segments with dynamic fallback
pub fn interpolated_bc(
    mach: f64,
    segments: &[(f64, f64)],
    inputs: Option<&BallisticInputs>,
) -> f64 {
    if segments.is_empty() {
        // Use dynamic fallback based on bullet characteristics if available
        if let Some(inputs) = inputs {
            let bc_type_str = match inputs.bc_type {
                crate::DragModel::G1 => "G1",
                crate::DragModel::G7 => "G7",
                _ => "G1", // Default to G1 for other models
            };
            return calculate_bc_fallback(
                Some(inputs.weight_grains),  // grains
                Some(inputs.caliber_inches), // inches
                Some(bc_type_str),
            );
        }
        return crate::constants::BC_FALLBACK_CONSERVATIVE; // Conservative fallback based on database analysis
    }

    if segments.len() == 1 {
        return segments[0].1;
    }

    // Ensure ascending-Mach order for interpolation. Fast path: when the segments are
    // already sorted (the common case — they are normalized once at construction), borrow
    // them and skip the per-call heap alloc + O(n log n) sort on the integration hot path.
    let sorted_segments: std::borrow::Cow<[(f64, f64)]> =
        if segments.windows(2).all(|w| w[0].0 <= w[1].0) {
            std::borrow::Cow::Borrowed(segments)
        } else {
            let mut v = segments.to_vec();
            v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            std::borrow::Cow::Owned(v)
        };

    // Handle out-of-range cases first
    if mach <= sorted_segments[0].0 {
        return sorted_segments[0].1;
    }
    if mach >= sorted_segments[sorted_segments.len() - 1].0 {
        return sorted_segments[sorted_segments.len() - 1].1;
    }

    // Find the appropriate segment using binary search
    let idx = sorted_segments.partition_point(|(m, _)| *m <= mach);
    if idx == 0 || idx >= sorted_segments.len() {
        // Should not happen given the checks above
        return sorted_segments[0].1;
    }

    let (mach1, bc1) = sorted_segments[idx - 1];
    let (mach2, bc2) = sorted_segments[idx];

    // Linear interpolation with safe division
    let denominator = mach2 - mach1;
    if denominator.abs() < crate::constants::MIN_DIVISION_THRESHOLD {
        return bc1; // Return first BC value if Mach values are identical
    }
    let t = (mach - mach1) / denominator;
    bc1 + t * (bc2 - bc1)
}

/// Get BC value for current velocity, supporting velocity-based BC segments
fn get_bc_for_velocity(velocity_fps: f64, inputs: &BallisticInputs, bc_used: f64) -> f64 {
    // Check if velocity-based BC segments are enabled
    if !inputs.use_bc_segments {
        return bc_used;
    }

    // Try direct BC segments data first
    if let Some(ref bc_segments_data) = inputs.bc_segments_data {
        for segment in bc_segments_data {
            if velocity_fps >= segment.velocity_min && velocity_fps <= segment.velocity_max {
                return segment.bc_value;
            }
        }
    }

    // Try BC estimation if we have bullet details but no segments. MBA-955: the estimation is
    // factored into estimate_bc_segments_for so the per-integration setup (build_inputs) can
    // pre-populate bc_segments_data ONCE rather than rebuilding it here every step.
    if let Some(segments) = estimate_bc_segments_for(inputs, bc_used) {
        for segment in &segments {
            if velocity_fps >= segment.velocity_min && velocity_fps <= segment.velocity_max {
                return segment.bc_value;
            }
        }
    }

    // Fallback to constant BC
    bc_used
}

/// Estimate velocity-BC segments from bullet characteristics (MBA-955). Extracted from
/// get_bc_for_velocity's slow path so the per-integration setup can compute the segments ONCE
/// (build_inputs pre-populates bc_segments_data) instead of rebuilding them — allocating a model
/// String and a segment Vec — on every derivative evaluation. Returns None when the bullet
/// details needed for estimation are absent (the caller then falls back to the constant BC). The
/// logic is byte-identical to the previous inline slow path.
pub(crate) fn estimate_bc_segments_for(
    inputs: &BallisticInputs,
    bc_used: f64,
) -> Option<Vec<crate::BCSegmentData>> {
    if !(inputs.bullet_diameter > 0.0 && inputs.bullet_mass > 0.0 && bc_used > 0.0) {
        return None;
    }
    // Model string from bullet_id or a generic weight-based description (unchanged).
    let model = if let Some(ref bullet_id) = inputs.bullet_id {
        bullet_id.clone()
    } else {
        format!("{}gr bullet", inputs.weight_grains as i32)
    };
    let bc_type_str = inputs.bc_type_str.as_deref().unwrap_or("G1");
    Some(BCSegmentEstimator::estimate_bc_segments(
        bc_used,
        inputs.caliber_inches,
        inputs.weight_grains,
        &model,
        bc_type_str,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_inputs() -> BallisticInputs {
        // SI-canonical geometry/mass (kg, meters) — same convention as the struct
        // docs and cli_api — plus the explicit imperial mirror fields
        // (caliber_inches/weight_grains) the stability/Magnus helpers read.
        BallisticInputs {
            muzzle_velocity: 800.0,         // m/s
            bc_value: 0.5,
            bullet_mass: 168.0 * 0.00006479891, // kg (168 gr)
            bullet_diameter: 0.308 * 0.0254,    // meters (.308 in)
            bullet_length: 1.215 * 0.0254,      // meters
            caliber_inches: 0.308,
            weight_grains: 168.0,
            altitude: 1000.0,
            ..Default::default()
        }
    }

    #[test]
    fn test_mba955_bc_segments_prepopulate_byte_identical() {
        // MBA-955: pre-populating bc_segments_data once (in build_inputs) must return
        // BYTE-IDENTICAL BC to the old per-step estimation. Build the slow-path inputs
        // (bc_segments_data = None -> get_bc_for_velocity estimates every call) and the
        // pre-populated inputs (bc_segments_data = estimate_bc_segments_for, the same helper
        // build_inputs now calls), and assert get_bc_for_velocity agrees bit-for-bit across the
        // whole velocity range.
        let mut slow = create_test_inputs();
        slow.use_bc_segments = true;
        slow.bc_segments_data = None;
        slow.bc_segments = None;

        let bc_used = slow.bc_value;
        let mut fast = slow.clone();
        fast.bc_segments_data = estimate_bc_segments_for(&fast, bc_used);
        assert!(
            fast.bc_segments_data.is_some(),
            "estimation should yield segments for a valid bullet"
        );

        for v in (200..=3500).step_by(50) {
            let vf = v as f64;
            let a = get_bc_for_velocity(vf, &slow, bc_used);
            let b = get_bc_for_velocity(vf, &fast, bc_used);
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "BC differs at {vf} fps: slow={a} fast={b}"
            );
        }
    }

    #[test]
    fn test_compute_derivatives_basic() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::zeros();
        // Use direct atmosphere values: (air_density, speed_of_sound, 0.0, 0.0)
        let atmos_params = (1.225, 340.0, 0.0, 0.0); // Standard air density and speed of sound
        let bc_used = 0.5;

        let result = compute_derivatives(
            pos,
            vel,
            &inputs,
            wind_vector,
            atmos_params,
            bc_used,
            None,
            0.0,
        );

        // Check that we get velocity and acceleration components
        assert_eq!(result.len(), 6);

        // Velocity components should match input velocity
        assert!((result[0] - vel[0]).abs() < 1e-10);
        assert!((result[1] - vel[1]).abs() < 1e-10);
        assert!((result[2] - vel[2]).abs() < 1e-10);

        // Should have gravitational acceleration
        assert!(result[4] < 0.0); // Negative y acceleration due to gravity

        // Should have drag acceleration opposing motion
        assert!(result[3] < 0.0); // Negative x acceleration due to drag
    }

    #[test]
    fn test_compute_derivatives_with_wind() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::new(10.0, 0.0, 0.0); // Tailwind
        let atmos_params = (1.225, 340.0, 0.0, 0.0); // Standard air density and speed of sound
        let bc_used = 0.5;

        let result = compute_derivatives(
            pos,
            vel,
            &inputs,
            wind_vector,
            atmos_params,
            bc_used,
            None,
            0.0,
        );

        // With tailwind, effective velocity should be lower, thus less drag
        // Just check that we have some drag (negative acceleration)
        assert!(result[3] < 0.0); // Should have drag
    }

    #[test]
    fn test_compute_derivatives_with_coriolis() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(800.0, 0.0, 0.0);
        let inputs = create_test_inputs();
        let wind_vector = Vector3::zeros();
        let atmos_params = (1.225, 340.0, 0.0, 0.0); // Standard air density and speed of sound
        let bc_used = 0.5;
        let omega = Vector3::new(0.0, 0.0, 7.2921e-5); // Earth's rotation

        let result = compute_derivatives(
            pos,
            vel,
            &inputs,
            wind_vector,
            atmos_params,
            bc_used,
            Some(omega),
            0.0,
        );

        // Should have Coriolis effect
        assert!(result[4].abs() > 1e-3); // Should have some y-component from Coriolis
    }

    #[test]
    fn test_interpolated_bc() {
        let segments = vec![(0.5, 0.4), (1.0, 0.5), (1.5, 0.6), (2.0, 0.5)];

        // Test exact matches
        assert!((interpolated_bc(1.0, &segments, None) - 0.5).abs() < 1e-10);

        // Test interpolation
        let bc_075 = interpolated_bc(0.75, &segments, None);
        assert!(bc_075 > 0.4 && bc_075 < 0.5);

        // Test out of range
        assert!((interpolated_bc(0.1, &segments, None) - 0.4).abs() < 1e-10);
        assert!((interpolated_bc(3.0, &segments, None) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_interpolated_bc_edge_cases() {
        // Empty segments
        assert!(
            (interpolated_bc(1.0, &[], None) - crate::constants::BC_FALLBACK_CONSERVATIVE).abs()
                < 1e-10
        );

        // Single segment
        let single = vec![(1.0, 0.7)];
        assert!((interpolated_bc(1.5, &single, None) - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_magnus_effect() {
        let pos = Vector3::new(0.0, 0.0, 0.0);
        let vel = Vector3::new(822.96, 0.0, 0.0); // 2700 fps
        let mut inputs = create_test_inputs();
        inputs.twist_rate = 10.0; // 1:10 twist
        inputs.is_twist_right = true;
        inputs.enable_magnus = true; // decoupled from enable_advanced_effects

        let wind_vector = Vector3::zeros();
        let atmos_params = (1.225, 340.0, 0.0, 0.0); // Standard air density and speed of sound
        let bc_used = 0.5;

        let result = compute_derivatives(
            pos,
            vel,
            &inputs,
            wind_vector,
            atmos_params,
            bc_used,
            None,
            0.0,
        );

        // Magnus is a small lateral (z) acceleration, positive (right) for RH twist,
        // using the proper yaw-of-repose force model (the old 1.8 fudge factor is gone).
        // For this .308/168gr/1:10 case at ~2700 fps it is ~0.003 m/s² — a fraction of
        // gravity, integrating to a sub-inch drift, consistent with the cli_api solver.
        assert!(
            result[5] > 0.0,
            "Magnus should drift right for RH twist, got {}",
            result[5]
        );
        assert!(
            result[5] < 0.05,
            "Magnus accel should be small/physical, got {}",
            result[5]
        );
    }

    #[test]
    fn test_magnus_moment_coefficient() {
        // Test at various Mach numbers with corrected coefficients
        assert!((calculate_magnus_moment_coefficient(0.5) - 0.030).abs() < 0.001); // Subsonic
        assert!((calculate_magnus_moment_coefficient(0.8) - 0.030).abs() < 0.001); // Start of transonic
        assert!((calculate_magnus_moment_coefficient(1.0) - 0.02625).abs() < 0.001); // Mid transonic
        assert!((calculate_magnus_moment_coefficient(1.2) - 0.015).abs() < 0.001); // End of transonic
        assert!((calculate_magnus_moment_coefficient(2.0) - 0.01653).abs() < 0.001);
        // Supersonic
    }
}
