use crate::pitch_damping::{
    calculate_damped_yaw_of_repose, calculate_gravity_yaw_of_repose,
    calculate_pitch_damping_moment, PitchDampingCoefficients,
};
use crate::spin_decay::{update_spin_rate, SpinDecayParameters};
use crate::BallisticInputs;
use std::f64::consts::PI;

/// Components of enhanced spin drift calculation
#[derive(Debug, Clone)]
pub struct SpinDriftComponents {
    pub spin_rate_rps: f64,          // Revolutions per second
    pub spin_rate_rad_s: f64,        // Radians per second
    pub stability_factor: f64,       // Gyroscopic stability (Sg)
    pub yaw_of_repose_rad: f64,      // Equilibrium yaw angle
    pub drift_rate_mps: f64,         // Lateral drift rate (m/s)
    pub total_drift_m: f64,          // Total drift at current time
    pub magnus_component_m: f64,     // Magnus effect contribution
    pub gyroscopic_component_m: f64, // Pure gyroscopic drift
    pub pitch_damping_moment: f64,   // Pitch damping moment (N⋅m)
    pub yaw_convergence_rate: f64,   // Convergence rate to equilibrium (rad/s)
    pub pitch_rate_rad_s: f64,       // Current pitch/yaw rate (rad/s)
}

/// Base Miller gyroscopic stability factor (no velocity/density correction).
/// All inputs imperial: caliber/length in inches, mass in grains, twist in
/// inches-per-turn. Returns 0.0 for non-positive inputs.
///   Sg = 30 m / (t^2 d^3 l (1 + l^2)),  t,l in calibers, d in inches, m in grains
pub(crate) fn miller_stability(
    caliber_in: f64,
    weight_gr: f64,
    twist_in: f64,
    length_in: f64,
) -> f64 {
    if caliber_in <= 0.0 || weight_gr <= 0.0 || twist_in <= 0.0 || length_in <= 0.0 {
        return 0.0;
    }
    let twist_cal = twist_in / caliber_in;
    let l_cal = length_in / caliber_in;
    let denom = twist_cal * twist_cal * caliber_in.powi(3) * l_cal * (1.0 + l_cal * l_cal);
    if denom == 0.0 {
        return 0.0;
    }
    30.0 * weight_gr / denom
}

/// Empirical Litz spin-drift MAGNITUDE in inches from the muzzle gyroscopic stability `sg`
/// and time of flight `t_s` (seconds):
///   `drift_inches = 1.25 * (Sg + 1.2) * t^1.83`
///
/// This is the single source of the Litz drift coefficient (MBA-1134). It is UNSIGNED — the
/// caller applies the twist-direction sign (see [`litz_drift_meters`]). cli_api::apply_spin_drift
/// and the fast / Monte-Carlo path both go through this so the three solver families stay
/// bit-identical on the drift math.
pub fn litz_drift_inches(sg: f64, t_s: f64) -> f64 {
    1.25 * (sg + 1.2) * t_s.powf(1.83)
}

/// Signed Litz spin drift in METERS along McCoy Z (lateral / windage). A right-hand twist
/// drifts to the right (+Z); a left-hand twist drifts left (-Z). MBA-1134.
pub fn litz_drift_meters(sg: f64, t_s: f64, is_twist_right: bool) -> f64 {
    let sign = if is_twist_right { 1.0 } else { -1.0 };
    sign * litz_drift_inches(sg, t_s) * 0.0254
}

/// Canonical muzzle gyroscopic stability Sg for the empirical Litz spin-drift model, shared by
/// cli_api and the fast / Monte-Carlo path so every solver family uses ONE Sg (MBA-1134, rank 31).
///
/// Delegates to [`crate::stability::compute_stability_coefficient`], the single source of truth
/// for Miller Sg. That INCLUDES the `(v/2800)^(1/3)` muzzle-velocity term (matching the reported
/// SG and the aerodynamic-jump Sg) and the linear Miller density correction `(T/T0)*(P0/P)`.
/// `temp_c` / `press_hpa` are the resolved muzzle atmosphere.
///
/// `compute_stability_coefficient` returns 0.0 when `bullet_length` is unset, so this first
/// substitutes a length estimate. MBA-1135: use the mass-based [`crate::stability::estimate_bullet_length_m`]
/// (falling back to the historical 4.5-caliber literal only when mass is unavailable) instead of a
/// mass-blind 4.5-caliber default, so heavier/longer bullets get a physically consistent Sg.
pub fn effective_sg_from_inputs(inputs: &BallisticInputs, temp_c: f64, press_hpa: f64) -> f64 {
    let mut eff = inputs.clone();
    if eff.bullet_length <= 0.0 && eff.bullet_diameter > 0.0 {
        let est = crate::stability::estimate_bullet_length_m(eff.bullet_diameter, eff.bullet_mass);
        eff.bullet_length = if est > 0.0 {
            est
        } else {
            4.5 * eff.bullet_diameter // 4.5 calibers, mass unavailable
        };
    }
    // atmo_params = (altitude, temp_c, press_hpa, _): compute_stability_coefficient reads only
    // temp_c and press_hpa (altitude is ignored there), so pass the resolved muzzle values.
    crate::stability::compute_stability_coefficient(&eff, (eff.altitude, temp_c, press_hpa, 0.0))
}

/// Calculate bullet spin rate from velocity and twist rate
pub fn calculate_spin_rate(velocity_mps: f64, twist_rate_inches: f64) -> (f64, f64) {
    if twist_rate_inches <= 0.0 {
        return (0.0, 0.0);
    }

    // Convert velocity to inches/second
    let velocity_ips = velocity_mps * 39.3701;

    // Calculate revolutions per second
    let spin_rate_rps = velocity_ips / twist_rate_inches;

    // Convert to radians per second
    let spin_rate_rad_s = spin_rate_rps * 2.0 * PI;

    (spin_rate_rps, spin_rate_rad_s)
}

/// Return muzzle-set spin rate and the dimensionless spin parameter at the current airspeed.
///
/// These acceleration kernels do not integrate roll as a state, so spin is held at its launch
/// value while translational airspeed changes. Modeling spin decay would require carrying roll
/// rate through the integrator rather than re-deriving it from current velocity.
pub(crate) fn calculate_magnus_spin_state(
    muzzle_velocity_mps: f64,
    current_velocity_mps: f64,
    twist_rate_inches: f64,
    caliber_m: f64,
) -> (f64, f64) {
    let (_, spin_rate_rad_s) = calculate_spin_rate(muzzle_velocity_mps, twist_rate_inches);
    let spin_parameter = if current_velocity_mps > 1e-9 {
        spin_rate_rad_s * caliber_m / (2.0 * current_velocity_mps)
    } else {
        0.0
    };
    (spin_rate_rad_s, spin_parameter)
}

/// Calculate dynamic gyroscopic stability factor using Miller formula
pub fn calculate_dynamic_stability(
    bullet_mass_grains: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    caliber_inches: f64,
    length_inches: f64,
    air_density_kg_m3: f64,
) -> f64 {
    if spin_rate_rad_s == 0.0 || velocity_mps == 0.0 {
        return 0.0;
    }

    // Convert velocity to fps for Miller formula
    let velocity_fps = velocity_mps * 3.28084;

    // Calculate twist rate in calibers
    if caliber_inches > 0.0 {
        // Back-calculate twist rate from spin rate
        let spin_rps = spin_rate_rad_s / (2.0 * PI);
        let velocity_ips = velocity_fps * 12.0; // inches per second
        let twist_inches = if spin_rps > 0.0 {
            velocity_ips / spin_rps
        } else {
            0.0
        };
        let twist_calibers = if twist_inches > 0.0 {
            twist_inches / caliber_inches
        } else {
            0.0
        };

        // Length to diameter ratio
        let length_calibers = if caliber_inches > 0.0 {
            length_inches / caliber_inches
        } else {
            0.0
        };

        // Miller stability formula (simplified)
        // Sg = 30 * m / (t^2 * d^3 * l * (1 + l^2))
        // Where: m = mass in grains, t = twist in calibers, d = diameter in inches
        //        l = length in calibers

        if twist_calibers == 0.0 || length_calibers == 0.0 {
            return 0.0;
        }

        let numerator = 30.0 * bullet_mass_grains;
        let denominator = twist_calibers.powi(2)
            * caliber_inches.powi(3)
            * length_calibers
            * (1.0 + length_calibers.powi(2));

        if denominator == 0.0 {
            return 0.0;
        }

        // Base stability
        let sg_base = numerator / denominator;

        // Velocity correction (compared to standard 2800 fps)
        let velocity_factor = (velocity_fps / 2800.0).powf(1.0 / 3.0);

        // Atmospheric correction (MBA-942): canonical Miller is LINEAR in density ratio
        // (FTP = (T/T0)*(P0/P) = rho0/rho), matching stability.rs and py_ballisticcalc. The
        // previous sqrt(1.225/rho) under-corrected Sg by ~14% at altitude. Standard
        // conditions: 59°F, 29.92 inHg = 1.225 kg/m³ (sqrt(1)=1, so sea level is unchanged).
        let density_factor = 1.225 / air_density_kg_m3;

        // Final stability
        sg_base * velocity_factor * density_factor
    } else {
        0.0
    }
}

/// Calculate the yaw of repose (equilibrium yaw angle)
pub fn calculate_yaw_of_repose(
    stability_factor: f64,
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    wind_velocity_mps: f64,
    pitch_rate_rad_s: f64,
    air_density_kg_m3: f64,
    caliber_inches: f64,
    length_inches: f64,
    mass_grains: f64,
    mach: f64,
    bullet_type: &str,
    use_pitch_damping: bool,
) -> (f64, f64) {
    if stability_factor <= 1.0 || spin_rate_rad_s == 0.0 {
        return (0.0, 0.0);
    }

    // Use enhanced calculation with pitch damping if requested
    if use_pitch_damping && mach > 0.0 {
        // Map bullet types for pitch damping
        let damping_type = match bullet_type.to_lowercase().as_str() {
            "match" => "match_boat_tail",
            "hunting" => "hunting",
            "fmj" => "fmj",
            "vld" => "vld",
            _ => "match_boat_tail",
        };

        return calculate_damped_yaw_of_repose(
            stability_factor,
            velocity_mps,
            spin_rate_rad_s,
            wind_velocity_mps,
            pitch_rate_rad_s,
            air_density_kg_m3,
            caliber_inches,
            length_inches,
            mass_grains,
            mach,
            damping_type,
        );
    }

    // Crosswind yaw is a muzzle transient, not a persistent equilibrium. The simple and damped
    // paths share the same gravity/gyroscopic repose angle; only the damped path reports a
    // convergence rate.
    let yaw_rad = calculate_gravity_yaw_of_repose(
        stability_factor,
        velocity_mps,
        spin_rate_rad_s,
        mass_grains * 0.00006479891,
        caliber_inches * 0.0254,
        length_inches * 0.0254,
    );

    (yaw_rad, 0.0)
}

/// Calculate Magnus effect contribution to drift
pub fn calculate_magnus_drift_component(
    velocity_mps: f64,
    spin_rate_rad_s: f64,
    yaw_rad: f64,
    air_density_kg_m3: f64,
    caliber_inches: f64,
    time_s: f64,
    mass_grains: f64,
) -> f64 {
    let diameter_m = caliber_inches * 0.0254;
    let mass_kg = mass_grains * 0.00006479891; // Convert grains to kg

    // Magnus force coefficient (empirical)
    // Varies with Mach number
    let mach = velocity_mps / 343.0; // Approximate speed of sound

    let cmag = if mach < 0.8 {
        0.25
    } else if mach < 1.2 {
        // Transonic reduction
        0.15
    } else {
        // Supersonic
        0.10 + 0.05 * ((mach - 1.2) / 2.0).min(1.0)
    };

    // Spin ratio
    let spin_ratio = (spin_rate_rad_s * diameter_m / 2.0) / velocity_mps;

    // Magnus force
    let magnus_force = if velocity_mps > 0.0 {
        cmag * spin_ratio
            * yaw_rad
            * 0.5
            * air_density_kg_m3
            * velocity_mps.powi(2)
            * PI
            * (diameter_m / 2.0).powi(2)
    } else {
        0.0
    };

    // Convert force to acceleration by dividing by mass
    let magnus_accel = magnus_force / mass_kg;

    // Drift over time (simplified - should integrate)

    0.5 * magnus_accel * time_s.powi(2)
}

/// Calculate pure gyroscopic drift (Poisson effect)
pub fn calculate_gyroscopic_drift(
    stability_factor: f64,
    _yaw_rad: f64,
    velocity_mps: f64,
    time_s: f64,
    is_right_twist: bool,
) -> f64 {
    if stability_factor <= 1.0 || time_s <= 0.0 {
        return 0.0;
    }

    // Litz formula is not reliable for subsonic flight. Disable it.
    let velocity_fps = velocity_mps * 3.28084;
    if velocity_fps < 1125.0 {
        return 0.0;
    }

    // Direction based on twist
    let sign = if is_right_twist { 1.0 } else { -1.0 };

    // Bryan Litz's empirical formula for spin drift
    let base_coefficient = 1.25 * (stability_factor + 1.2);
    let time_factor = time_s.powf(1.83);
    let drift_in = sign * base_coefficient * time_factor;

    // Convert to meters

    drift_in * 0.0254
}

/// Calculate enhanced spin drift with all components.
///
/// DEPRECATED (MBA-1134): this in-integration acceleration model is no longer wired into any
/// solver path — spin drift is now the single canonical empirical Litz post-process
/// ([`litz_drift_meters`], via [`effective_sg_from_inputs`]). Retained for backward compatibility
/// and unit tests only. Do NOT reintroduce it into an integration loop alongside the Litz model,
/// or lateral drift will be double-counted.
pub fn calculate_enhanced_spin_drift(
    bullet_mass: f64,
    velocity_mps: f64,
    twist_rate: f64,
    bullet_diameter: f64,
    bullet_length: f64,
    is_twist_right: bool,
    time_s: f64,
    air_density: f64,
    crosswind_mps: f64,
    pitch_rate_rad_s: f64,
    use_pitch_damping: bool,
) -> SpinDriftComponents {
    // Calculate initial spin rate (at muzzle)
    let muzzle_velocity = velocity_mps; // Assuming we're passed muzzle velocity
    let (_initial_spin_rps, initial_spin_rad_s) = calculate_spin_rate(muzzle_velocity, twist_rate);

    // Apply spin decay based on time of flight
    let decay_params = SpinDecayParameters::from_bullet_type("match"); // Default to match for now
    let current_spin_rad_s = update_spin_rate(
        initial_spin_rad_s,
        time_s,
        velocity_mps,
        air_density,
        bullet_mass, // already grains (update_spin_rate wants mass_grains)
        bullet_diameter,
        bullet_length,
        Some(&decay_params),
    );

    let spin_rps = current_spin_rad_s / (2.0 * PI);
    let spin_rad_s = current_spin_rad_s;

    // Calculate dynamic stability
    let stability = calculate_dynamic_stability(
        bullet_mass,
        velocity_mps,
        spin_rad_s,
        bullet_diameter,
        bullet_length,
        air_density,
    );

    // Calculate Mach number for pitch damping
    let mach = velocity_mps / 343.0; // Approximate speed of sound

    // Determine bullet type (default to match for now)
    let bullet_type = "match";

    // Calculate yaw of repose with pitch damping
    let (yaw_rad, convergence_rate) = calculate_yaw_of_repose(
        stability,
        velocity_mps,
        spin_rad_s,
        crosswind_mps,
        pitch_rate_rad_s,
        air_density,
        bullet_diameter,
        bullet_length,
        bullet_mass,
        mach,
        bullet_type,
        use_pitch_damping,
    );

    // Calculate Magnus component
    let magnus_drift = calculate_magnus_drift_component(
        velocity_mps,
        spin_rad_s,
        yaw_rad,
        air_density,
        bullet_diameter,
        time_s,
        bullet_mass,
    );

    // Calculate gyroscopic component
    let gyro_drift =
        calculate_gyroscopic_drift(stability, yaw_rad, velocity_mps, time_s, is_twist_right);

    // Total drift. gyro_drift already carries the twist-direction sign (from
    // calculate_gyroscopic_drift); sign the Magnus term to the SAME convention so both
    // contributions are consistently directed before being summed. (magnus_component_m is
    // kept unsigned below for backward compatibility.)
    let twist_sign = if is_twist_right { 1.0 } else { -1.0 };
    let total_drift = twist_sign * magnus_drift + gyro_drift;

    // Drift rate (derivative)
    let drift_rate = if time_s > 0.0 {
        total_drift / time_s
    } else {
        0.0
    };

    // Calculate pitch damping moment if using enhanced model
    let pitch_damping_moment = if use_pitch_damping && mach > 0.0 {
        let coeffs = PitchDampingCoefficients::from_bullet_type(bullet_type);
        calculate_pitch_damping_moment(
            pitch_rate_rad_s,
            velocity_mps,
            air_density,
            bullet_diameter * 0.0254, // Convert to meters
            bullet_length * 0.0254,   // Convert to meters
            mach,
            &coeffs,
        )
    } else {
        0.0
    };

    SpinDriftComponents {
        spin_rate_rps: spin_rps,
        spin_rate_rad_s: spin_rad_s,
        stability_factor: stability,
        yaw_of_repose_rad: yaw_rad,
        drift_rate_mps: drift_rate,
        total_drift_m: total_drift,
        magnus_component_m: magnus_drift,
        gyroscopic_component_m: gyro_drift,
        pitch_damping_moment,
        yaw_convergence_rate: convergence_rate,
        pitch_rate_rad_s,
    }
}

/// Apply enhanced spin drift acceleration to derivatives.
///
/// DEPRECATED (MBA-1134): companion to [`calculate_enhanced_spin_drift`]; no longer called by any
/// integration path. Spin drift is now the canonical Litz post-process ([`litz_drift_meters`]).
/// Retained for backward compatibility and unit tests only.
pub fn apply_enhanced_spin_drift(
    derivatives: &mut [f64; 6],
    spin_components: &SpinDriftComponents,
    time_s: f64,
    _is_right_twist: bool,
) {
    if time_s > 0.1 {
        // Calculate acceleration from drift
        // Using second derivative of position
        let spin_accel_z = 2.0 * spin_components.drift_rate_mps / time_s;

        // drift_rate_mps already carries the twist-direction sign (set in
        // calculate_enhanced_spin_drift), so apply it directly. Multiplying by the twist
        // sign again here previously CANCELED the gyroscopic sign (sign^2 = +1), so left-
        // and right-twist barrels pushed spin drift the same way.
        derivatives[5] += spin_accel_z;
    }
}

/// Simplified interface for compatibility with existing code
pub fn compute_enhanced_spin_drift_simple(
    time_s: f64,
    stability: f64,
    velocity_mps: f64,
    twist_rate: f64,
    is_twist_right: bool,
    _caliber: f64,
) -> f64 {
    if twist_rate <= 0.0 {
        return 0.0;
    }

    // Calculate initial spin rate
    let (_, initial_spin_rad_s) = calculate_spin_rate(velocity_mps, twist_rate);

    // Apply simple spin decay (assume 175gr bullet)
    let decay_params = SpinDecayParameters::from_bullet_type("match");
    let spin_rad_s = update_spin_rate(
        initial_spin_rad_s,
        time_s,
        velocity_mps,
        1.225, // Standard air density
        175.0, // Standard bullet weight
        _caliber,
        1.3, // Standard bullet length
        Some(&decay_params),
    );

    // Estimate yaw of repose (use simple model for compatibility)
    let (yaw_rad, _) = calculate_yaw_of_repose(
        stability,
        velocity_mps,
        spin_rad_s,
        0.0,
        0.0,
        1.225,
        _caliber,
        1.3,
        175.0,
        velocity_mps / 343.0,
        "match",
        false,
    );

    // Calculate gyroscopic drift (primary component)

    calculate_gyroscopic_drift(stability, yaw_rad, velocity_mps, time_s, is_twist_right)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_spin_rate() {
        // Test with 1:10 twist at 800 m/s
        let (rps, rad_s) = calculate_spin_rate(800.0, 10.0);

        // 800 m/s = 31496 in/s, divided by 10 = 3149.6 rps
        assert!((rps - 3149.6).abs() < 1.0);
        assert!((rad_s - rps * 2.0 * PI).abs() < 0.1);

        // Test with zero twist rate
        let (rps_zero, rad_s_zero) = calculate_spin_rate(800.0, 0.0);
        assert_eq!(rps_zero, 0.0);
        assert_eq!(rad_s_zero, 0.0);
    }

    #[test]
    fn magnus_spin_is_set_at_muzzle_while_spin_parameter_grows_downrange() {
        let muzzle_velocity = 800.0;
        let (muzzle_spin, muzzle_parameter) =
            calculate_magnus_spin_state(muzzle_velocity, muzzle_velocity, 10.0, 0.00782);
        let (downrange_spin, downrange_parameter) =
            calculate_magnus_spin_state(muzzle_velocity, muzzle_velocity / 2.0, 10.0, 0.00782);

        assert_eq!(downrange_spin.to_bits(), muzzle_spin.to_bits());
        assert!((downrange_parameter / muzzle_parameter - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_calculate_dynamic_stability() {
        let sg = calculate_dynamic_stability(
            168.0,   // grains
            800.0,   // m/s
            19792.0, // rad/s (from 1:10 twist)
            0.308,   // inches
            1.2,     // inches
            1.225,   // kg/m³
        );

        // Should be > 1.0 for stable bullet
        assert!(sg > 1.0);
        assert!(sg < 10.0); // Reasonable upper bound
    }

    #[test]
    fn test_calculate_yaw_of_repose() {
        let (yaw, _) = calculate_yaw_of_repose(
            2.5,     // Sg
            800.0,   // velocity m/s
            19792.0, // spin rate rad/s
            10.0,    // crosswind m/s
            0.0,     // pitch rate
            1.225,   // air density
            0.308,   // caliber
            1.2,     // length
            168.0,   // mass
            2.33,    // mach
            "match", // bullet type
            false,   // use pitch damping
        );

        // Should be small but non-zero
        assert!(yaw.abs() > 0.0);
        assert!(yaw.abs() < 0.1); // Less than ~6 degrees
    }

    #[test]
    fn yaw_of_repose_increases_with_physical_spin_stability() {
        let calculate = |stability_factor, spin_rate_rad_s, use_pitch_damping| {
            calculate_yaw_of_repose(
                stability_factor,
                300.0,
                spin_rate_rad_s,
                0.0,
                0.01,
                1.225,
                0.308,
                1.3,
                175.0,
                300.0 / 343.0,
                "match",
                use_pitch_damping,
            )
            .0
        };
        let low_stability: f64 = 1.1;
        let high_stability: f64 = 4.0;
        let low_spin = 19_000.0;
        let high_spin = low_spin * (high_stability / low_stability).sqrt();
        let expected_ratio = (high_stability / low_stability).sqrt();

        for use_pitch_damping in [false, true] {
            let low = calculate(low_stability, low_spin, use_pitch_damping);
            let high = calculate(high_stability, high_spin, use_pitch_damping);

            assert!(high > low, "yaw decreased as physical spin/Sg increased");
            assert!(
                (high / low - expected_ratio).abs() <= expected_ratio * 1e-12,
                "yaw stability ratio was {}, expected {expected_ratio}",
                high / low
            );
        }
    }

    #[test]
    fn yaw_of_repose_has_inverse_cube_velocity_scaling() {
        let calculate = |stability_factor, velocity_mps, use_pitch_damping| {
            calculate_yaw_of_repose(
                stability_factor,
                velocity_mps,
                19_000.0,
                0.0,
                0.01,
                1.225,
                0.308,
                1.3,
                175.0,
                velocity_mps / 343.0,
                "match",
                use_pitch_damping,
            )
            .0
        };

        // With spin fixed, halving velocity increases physical Sg by four. The classical
        // reduction yaw = 4*Iy*Sg*g/(Ix*p*V) must therefore increase by 4*2 = 8 (V^-3).
        for use_pitch_damping in [false, true] {
            let fast = calculate(1.5, 800.0, use_pitch_damping);
            let slow = calculate(6.0, 400.0, use_pitch_damping);

            assert!((slow / fast - 8.0).abs() <= 8e-12);
        }
    }

    #[test]
    fn crosswind_does_not_change_yaw_of_repose_in_either_model() {
        let calculate = |wind_velocity_mps, use_pitch_damping| {
            calculate_yaw_of_repose(
                2.5,
                300.0,
                19_000.0,
                wind_velocity_mps,
                0.01,
                1.225,
                0.308,
                1.3,
                175.0,
                0.875,
                "match",
                use_pitch_damping,
            )
            .0
        };

        let simple_calm = calculate(0.0, false);
        let simple_windy = calculate(10.0, false);
        let damped_calm = calculate(0.0, true);
        let damped_windy = calculate(10.0, true);

        assert_eq!(simple_windy.to_bits(), simple_calm.to_bits());
        assert_eq!(damped_windy.to_bits(), damped_calm.to_bits());
        assert_eq!(damped_calm.to_bits(), simple_calm.to_bits());
        assert!(simple_calm > 0.0 && simple_calm < 0.003);
    }

    #[test]
    fn test_enhanced_spin_drift_calculation() {
        let components = calculate_enhanced_spin_drift(
            168.0, // mass grains
            800.0, // velocity m/s
            10.0,  // twist rate inches
            0.308, // caliber inches
            1.2,   // length inches
            true,  // right twist
            1.0,   // time s
            1.225, // air density
            10.0,  // crosswind
            0.0,   // pitch rate
            false, // use pitch damping
        );

        // Should produce non-zero drift
        assert!(components.total_drift_m.abs() > 0.0);
        assert!(components.spin_rate_rps > 0.0);
        assert!(components.stability_factor > 0.0);
    }

    #[test]
    fn test_litz_drift_helpers_sign_and_magnitude() {
        // litz_drift_inches is unsigned and matches 1.25*(Sg+1.2)*t^1.83 exactly.
        let sg = 2.0_f64;
        let t = 1.5_f64;
        let expected_in = 1.25 * (sg + 1.2) * t.powf(1.83);
        assert!((litz_drift_inches(sg, t) - expected_in).abs() < 1e-12);
        // litz_drift_meters applies the twist sign and the inch->meter conversion.
        let right = litz_drift_meters(sg, t, true);
        let left = litz_drift_meters(sg, t, false);
        assert!((right - expected_in * 0.0254).abs() < 1e-12);
        assert!((right + left).abs() < 1e-12, "left twist must mirror right");
        assert!(right > 0.0 && left < 0.0);
    }

    #[test]
    fn test_effective_sg_from_inputs_includes_velocity_term_and_length_fallback() {
        // effective_sg_from_inputs must (1) equal compute_stability_coefficient, (2) INCLUDE the
        // (v/2800)^(1/3) muzzle-velocity term (so it differs from the bare geometric miller_stability),
        // and (3) apply the 4.5-caliber length fallback when bullet_length is unset.
        let inputs = BallisticInputs {
            muzzle_velocity: 800.0, // 2624.7 fps -> velocity term < 1.0
            bullet_mass: 175.0 * 0.00006479891,
            bullet_diameter: 0.308 * 0.0254,
            bullet_length: 1.24 * 0.0254,
            twist_rate: 10.0,
            ..Default::default()
        };

        let temp_c = 15.0;
        let press_hpa = 1013.25;
        let sg = effective_sg_from_inputs(&inputs, temp_c, press_hpa);

        // (1) identical to compute_stability_coefficient on the same (length-filled) inputs.
        let direct =
            crate::stability::compute_stability_coefficient(&inputs, (0.0, temp_c, press_hpa, 0.0));
        assert!((sg - direct).abs() < 1e-12, "sg {sg} != direct {direct}");

        // (2) includes the velocity term: at sea-level standard the density factor is 1.0, so
        //     sg == bare_geometric_Sg * (v_fps/2800)^(1/3).
        let d_in = inputs.bullet_diameter / 0.0254;
        let m_gr = inputs.bullet_mass / 0.00006479891;
        let l_in = inputs.bullet_length / 0.0254;
        let bare = miller_stability(d_in, m_gr, inputs.twist_rate, l_in);
        let vel_corr = (inputs.muzzle_velocity * 3.28084 / 2800.0).powf(1.0 / 3.0);
        assert!(vel_corr < 1.0, "muzzle < 2800 fps should shrink Sg");
        assert!(
            (sg - bare * vel_corr).abs() < 1e-6,
            "sg {sg} != bare {bare} * vel_corr {vel_corr}"
        );

        // (3) MBA-1135: the zero-length fallback now uses the mass-based length estimate
        // (crate::stability::estimate_bullet_length_m), NOT the old mass-blind 4.5-caliber
        // length. Zero length must reproduce an explicit estimate-length input.
        let mut no_len = inputs.clone();
        no_len.bullet_length = 0.0;
        let sg_fallback = effective_sg_from_inputs(&no_len, temp_c, press_hpa);
        let mut explicit = inputs.clone();
        explicit.bullet_length =
            crate::stability::estimate_bullet_length_m(inputs.bullet_diameter, inputs.bullet_mass);
        let sg_explicit = effective_sg_from_inputs(&explicit, temp_c, press_hpa);
        assert!(
            (sg_fallback - sg_explicit).abs() < 1e-12,
            "zero-length fallback {sg_fallback} != explicit estimate-length {sg_explicit}"
        );
        assert!(sg_fallback > 0.0);
        // And it must differ from the retired 4.5-caliber default (the whole point of MBA-1135).
        let mut old_default = inputs.clone();
        old_default.bullet_length = 4.5 * old_default.bullet_diameter;
        let sg_old = effective_sg_from_inputs(&old_default, temp_c, press_hpa);
        assert!(
            (sg_fallback - sg_old).abs() > 1e-6,
            "mass-based fallback should differ from the old 4.5-cal default"
        );
    }

    #[test]
    fn test_miller_stability_308_168gr() {
        // .308, 168 gr, 1:12 twist, ~1.215 in length -> base Sg (no velocity/density correction)
        // Formula: Sg = 30*m / (t^2 * d^3 * l * (1+l^2)), t and l in calibers
        // twist_cal = 12/0.308 = 38.96, l_cal = 1.215/0.308 = 3.94 -> Sg ~ 1.74
        let sg = miller_stability(0.308, 168.0, 12.0, 1.215);
        assert!(sg > 1.5 && sg < 2.0, "expected base Sg ~1.74, got {}", sg);
    }

    #[test]
    fn test_miller_stability_invalid_inputs_zero() {
        assert_eq!(miller_stability(0.0, 168.0, 12.0, 1.2), 0.0);
        assert_eq!(miller_stability(0.308, 0.0, 12.0, 1.2), 0.0);
        assert_eq!(miller_stability(0.308, 168.0, 0.0, 1.2), 0.0);
        assert_eq!(miller_stability(0.308, 168.0, 12.0, 0.0), 0.0);
    }

    #[test]
    fn test_opposite_twist_directions() {
        // Right twist
        let right_drift = calculate_enhanced_spin_drift(
            168.0, 800.0, 10.0, 0.308, 1.2, true, 1.0, 1.225, 0.0, 0.0, false,
        );

        // Left twist
        let left_drift = calculate_enhanced_spin_drift(
            168.0, 800.0, 10.0, 0.308, 1.2, false, 1.0, 1.225, 0.0, 0.0, false,
        );

        // Should have opposite signs for gyroscopic component
        assert!(right_drift.gyroscopic_component_m * left_drift.gyroscopic_component_m < 0.0);
        assert!(
            (right_drift.gyroscopic_component_m.abs() - left_drift.gyroscopic_component_m.abs())
                .abs()
                < 0.001
        );
    }

    #[test]
    fn test_applied_spin_drift_flips_with_twist() {
        // Regression: the APPLIED lateral acceleration (derivatives[5]) must reverse
        // direction with the twist hand. apply_enhanced_spin_drift previously multiplied by
        // the twist sign a second time, canceling the gyroscopic sign so left- and right-
        // twist barrels pushed the same way. The existing test only checks the component
        // field, never the applied derivative.
        let time_s = 1.0;
        let right = calculate_enhanced_spin_drift(
            168.0, 800.0, 10.0, 0.308, 1.2, true, time_s, 1.225, 0.0, 0.0, false,
        );
        let left = calculate_enhanced_spin_drift(
            168.0, 800.0, 10.0, 0.308, 1.2, false, time_s, 1.225, 0.0, 0.0, false,
        );

        let mut d_right = [0.0_f64; 6];
        let mut d_left = [0.0_f64; 6];
        apply_enhanced_spin_drift(&mut d_right, &right, time_s, true);
        apply_enhanced_spin_drift(&mut d_left, &left, time_s, false);

        assert!(d_right[5].abs() > 0.0, "expected non-zero spin drift accel");
        assert!(d_left[5].abs() > 0.0, "expected non-zero spin drift accel");
        assert!(
            d_right[5] * d_left[5] < 0.0,
            "expected opposite-sign lateral accel for opposite twist, got {} and {}",
            d_right[5],
            d_left[5]
        );
    }
}
