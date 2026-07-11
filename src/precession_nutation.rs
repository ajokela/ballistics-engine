//! Precession and Nutation Physics for Ballistic Projectiles
//!
//! This module implements the complex angular motion of spinning projectiles:
//! - Precession: Slow coning motion of the projectile axis
//! - Nutation: Fast oscillatory motion superimposed on precession
//! - Angular momentum conservation
//! - Gyroscopic effects

// Precession and nutation modeling - now integrated!

use crate::pitch_damping::{calculate_pitch_damping_moment, PitchDampingCoefficients};

/// Estimate the projectile's longitudinal and transverse moments of inertia from its geometry.
///
/// The engine models a typical pointed projectile as an ogive for both axes. Invalid geometry
/// returns zero inertias so callers naturally take the existing no-motion guard instead of
/// propagating non-finite angular state.
pub(crate) fn projectile_moments_of_inertia(
    mass_kg: f64,
    caliber_m: f64,
    length_m: f64,
) -> (f64, f64) {
    if !mass_kg.is_finite()
        || mass_kg <= 0.0
        || !caliber_m.is_finite()
        || caliber_m <= 0.0
        || !length_m.is_finite()
        || length_m <= 0.0
    {
        return (0.0, 0.0);
    }

    let spin_inertia =
        crate::spin_decay::calculate_moment_of_inertia(mass_kg, caliber_m, length_m, "ogive");
    let transverse_inertia = crate::pitch_damping::calculate_transverse_moment_of_inertia(
        mass_kg, caliber_m, length_m, "ogive",
    );

    if spin_inertia.is_finite()
        && spin_inertia > 0.0
        && transverse_inertia.is_finite()
        && transverse_inertia > 0.0
    {
        (spin_inertia, transverse_inertia)
    } else {
        (0.0, 0.0)
    }
}

/// Complete angular state of the projectile
#[derive(Debug, Clone, Copy)]
pub struct AngularState {
    pub pitch_angle: f64,      // Angle between axis and velocity (rad)
    pub yaw_angle: f64,        // Angle in plane perpendicular to velocity (rad)
    pub pitch_rate: f64,       // Rate of pitch angle change (rad/s)
    pub yaw_rate: f64,         // Rate of yaw angle change (rad/s)
    pub precession_angle: f64, // Cumulative precession angle (rad)
    pub nutation_phase: f64,   // Phase of nutation oscillation (rad)
}

/// Parameters for precession and nutation calculations
#[derive(Debug, Clone)]
pub struct PrecessionNutationParams {
    // Projectile properties
    pub mass_kg: f64,
    pub caliber_m: f64,
    pub length_m: f64,
    pub spin_rate_rad_s: f64,

    // Moments of inertia
    pub spin_inertia: f64,       // About longitudinal axis
    pub transverse_inertia: f64, // About transverse axis

    // Flight conditions
    pub velocity_mps: f64,
    pub air_density_kg_m3: f64,
    pub mach: f64,

    // Damping coefficients
    pub pitch_damping_coeff: f64,
    pub nutation_damping_factor: f64, // Fraction of critical damping
}

impl Default for PrecessionNutationParams {
    fn default() -> Self {
        Self {
            mass_kg: 0.01134, // 175 grains
            caliber_m: 0.00782,
            length_m: 0.033,
            spin_rate_rad_s: 17522.0,
            spin_inertia: 6.94e-8,
            transverse_inertia: 9.13e-7,
            velocity_mps: 850.0,
            air_density_kg_m3: 1.225,
            mach: 2.48,
            pitch_damping_coeff: PitchDampingCoefficients::default().subsonic,
            nutation_damping_factor: 0.05,
        }
    }
}

/// The two epicyclic yaw-arm angular frequencies (rad/s) for a gyroscopically stable projectile:
/// the FAST mode (nutation) and the SLOW mode (precession). Standard linearized aeroballistic
/// result from the spinning-projectile yaw equation:
///   phi_{fast,slow} = (Ix * p / 2 Iy) * [1 ± sqrt(1 - 1/Sg)]
/// where Ix/Iy are the spin/transverse moments of inertia, p the spin rate, Sg the (dimensionless)
/// gyroscopic stability factor. Returns (0, 0) when Sg <= 1 (no real epicyclic motion — the
/// projectile is not gyroscopically stable) or the transverse inertia is zero. (MBA-941: the
/// previous per-frequency formulas were dimensionally inconsistent — rad/m and length — and ad hoc.)
pub fn epicyclic_frequencies(
    spin_inertia: f64,
    transverse_inertia: f64,
    spin_rate_rad_s: f64,
    stability_factor: f64,
) -> (f64, f64) {
    if stability_factor <= 1.0 || transverse_inertia == 0.0 {
        return (0.0, 0.0);
    }
    // Ix * p / 2 Iy  [rad/s] — the mean of the two arm rates.
    let arm = (spin_inertia * spin_rate_rad_s) / (2.0 * transverse_inertia);
    let disc = (1.0 - 1.0 / stability_factor).sqrt();
    (arm * (1.0 + disc), arm * (1.0 - disc)) // (fast = nutation, slow = precession)
}

/// Slow-mode (precession) angular frequency in rad/s — the slow coning of the spin axis:
/// phi_slow = (Ix p / 2 Iy)(1 - sqrt(1 - 1/Sg)).
pub fn calculate_precession_frequency(
    spin_rate_rad_s: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
    stability_factor: f64,
) -> f64 {
    epicyclic_frequencies(spin_inertia, transverse_inertia, spin_rate_rad_s, stability_factor).1
}

/// Fast-mode (nutation) angular frequency in rad/s:
/// phi_fast = (Ix p / 2 Iy)(1 + sqrt(1 - 1/Sg)).
pub fn calculate_nutation_frequency(
    spin_rate_rad_s: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
    stability_factor: f64,
) -> f64 {
    epicyclic_frequencies(spin_inertia, transverse_inertia, spin_rate_rad_s, stability_factor).0
}

/// Calculate nutation amplitude with exponential damping
pub fn calculate_nutation_amplitude(
    initial_disturbance_rad: f64,
    time_s: f64,
    nutation_frequency: f64,
    damping_factor: f64,
    spin_rate_rad_s: f64,
) -> f64 {
    if nutation_frequency == 0.0 || spin_rate_rad_s == 0.0 {
        return 0.0;
    }

    // Damping rate
    let damping_rate = damping_factor * nutation_frequency;

    // Exponential decay
    let amplitude = initial_disturbance_rad * (-damping_rate * time_s).exp();

    // Clamp to reasonable bounds
    amplitude.min(0.1) // Max 0.1 rad (~5.7 degrees)
}

/// Calculate the combined precession and nutation motion
pub fn calculate_combined_angular_motion(
    params: &PrecessionNutationParams,
    angular_state: &AngularState,
    time_s: f64,
    dt: f64,
    initial_disturbance: f64,
) -> AngularState {
    // MBA-198: Guard against division by zero in stability calculation
    if params.transverse_inertia == 0.0
        || params.velocity_mps == 0.0
        || params.length_m == 0.0
        || !params.air_density_kg_m3.is_finite()
        || params.air_density_kg_m3 <= 0.0
    {
        // Return unchanged state if invalid parameters
        return *angular_state;
    }

    // Dimensionless gyroscopic stability factor via the engine's canonical dynamic Miller
    // calculation, including the (V/2800 fps)^(1/3) and rho0/rho corrections. Use spin magnitude
    // for stability (which depends on p^2); the signed rate below still controls phase direction.
    let caliber_in = params.caliber_m / 0.0254;
    let length_in = params.length_m / 0.0254;
    let mass_gr = params.mass_kg / 0.00006479891;
    let stability = crate::spin_drift::calculate_dynamic_stability(
        mass_gr,
        params.velocity_mps,
        params.spin_rate_rad_s.abs(),
        caliber_in,
        length_in,
        params.air_density_kg_m3,
    );

    // Precession (slow) and nutation (fast) angular frequencies, both rad/s.
    let omega_p = calculate_precession_frequency(
        params.spin_rate_rad_s,
        params.spin_inertia,
        params.transverse_inertia,
        stability,
    );
    let omega_n = calculate_nutation_frequency(
        params.spin_rate_rad_s,
        params.spin_inertia,
        params.transverse_inertia,
        stability,
    );

    // Nutation amplitude (decaying)
    let nutation_amp = calculate_nutation_amplitude(
        initial_disturbance,
        time_s,
        omega_n,
        params.nutation_damping_factor,
        params.spin_rate_rad_s,
    );

    // Update precession angle
    let new_precession_angle = angular_state.precession_angle + omega_p * dt;

    // Update nutation phase
    let new_nutation_phase = angular_state.nutation_phase + omega_n * dt;

    // Calculate pitch damping moment
    let pitch_moment = calculate_pitch_damping_moment(
        angular_state.pitch_rate,
        params.velocity_mps,
        params.air_density_kg_m3,
        params.caliber_m,
        params.length_m,
        params.mach,
        &PitchDampingCoefficients {
            subsonic: params.pitch_damping_coeff,
            ..Default::default()
        },
    );

    // Angular acceleration from damping
    let pitch_accel = pitch_moment / params.transverse_inertia;

    // Update angular rates
    let new_pitch_rate = angular_state.pitch_rate + pitch_accel * dt;

    // MBA-941: bounded epicyclic yaw. Previously `total_yaw = yaw_angle + nutation_amp*sin(phase)`
    // re-added the nutation to the carried-forward yaw every step, so the yaw random-walked and the
    // precession rate (yaw_rate) was never actually integrated. The yaw is now a bounded function
    // of the cumulative precession/nutation PHASES — a slow precession arm plus the damped fast
    // nutation arm — and yaw_rate is its true time derivative.
    let coning_amp = initial_disturbance;
    let total_yaw =
        coning_amp * new_precession_angle.cos() + nutation_amp * new_nutation_phase.sin();
    let new_yaw_rate = -coning_amp * omega_p * new_precession_angle.sin()
        + nutation_amp * omega_n * new_nutation_phase.cos();

    // Pitch angle evolves more slowly
    let new_pitch = angular_state.pitch_angle + new_pitch_rate * dt;

    AngularState {
        pitch_angle: new_pitch,
        yaw_angle: total_yaw,
        pitch_rate: new_pitch_rate,
        yaw_rate: new_yaw_rate,
        precession_angle: new_precession_angle,
        nutation_phase: new_nutation_phase,
    }
}

/// Calculate the epicyclic (combined precession + nutation) motion
pub fn calculate_epicyclic_motion(
    spin_inertia: f64,
    transverse_inertia: f64,
    spin_rate_rad_s: f64,
    stability_factor: f64,
    time_s: f64,
    initial_yaw_rad: f64,
) -> (f64, f64) {
    // MBA-198: Guard against division by zero
    if stability_factor <= 1.0 || spin_rate_rad_s == 0.0 {
        // Unstable or no spin - no regular motion
        return (initial_yaw_rad, initial_yaw_rad);
    }

    // Fast (nutation) and slow (precession) angular frequencies, both rad/s (MBA-941).
    let (omega_fast, omega_slow) =
        epicyclic_frequencies(spin_inertia, transverse_inertia, spin_rate_rad_s, stability_factor);

    // Amplitude ratio (fast/slow) — the nutation arm shrinks as stability grows.
    let amplitude_ratio = 1.0 / stability_factor;

    // Damping (exponential decay of fast mode)
    let damping_factor = 0.1; // Typical value
    let fast_amplitude = amplitude_ratio * (-damping_factor * omega_fast * time_s).exp();

    // Combined motion
    let slow_phase = omega_slow * time_s;
    let fast_phase = omega_fast * time_s;

    // Epicyclic coordinates
    let yaw = initial_yaw_rad * (slow_phase.cos() + fast_amplitude * fast_phase.cos());
    let pitch = initial_yaw_rad * (slow_phase.sin() + fast_amplitude * fast_phase.sin());

    (pitch, yaw)
}

/// Calculate the limit cycle yaw angle
pub fn calculate_limit_cycle_yaw(
    velocity_mps: f64,
    _spin_rate_rad_s: f64,
    stability_factor: f64,
    crosswind_mps: f64,
) -> f64 {
    // Base yaw from crosswind
    let wind_yaw = if crosswind_mps != 0.0 && velocity_mps > 0.0 {
        (crosswind_mps / velocity_mps).atan()
    } else {
        0.0
    };

    // Yaw of repose (equilibrium)
    let yaw_of_repose = if stability_factor > 1.0 {
        // Typical value for spin-stabilized projectiles
        let repose_factor = 1.0 / (1.0 + 0.5 * (stability_factor - 1.0));
        0.002 * repose_factor // ~0.1 degrees nominal
    } else {
        0.01 // Larger for marginally stable
    };

    wind_yaw + yaw_of_repose
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projectile_inertias_match_reference_and_geometry_scaling() {
        let mass = 0.01134;
        let caliber = 0.00782;
        let length = 0.033;
        let (spin, transverse) = projectile_moments_of_inertia(mass, caliber, length);

        assert!((spin / 6.94e-8 - 1.0).abs() < 0.01);
        assert!((transverse / 9.13e-7 - 1.0).abs() < 0.01);

        let (double_mass_spin, double_mass_transverse) =
            projectile_moments_of_inertia(2.0 * mass, caliber, length);
        assert!((double_mass_spin / spin - 2.0).abs() < 1e-12);
        assert!((double_mass_transverse / transverse - 2.0).abs() < 1e-12);

        let (double_caliber_spin, double_caliber_transverse) =
            projectile_moments_of_inertia(mass, 2.0 * caliber, length);
        assert!((double_caliber_spin / spin - 4.0).abs() < 1e-12);
        assert!(double_caliber_transverse > transverse);

        let (double_length_spin, double_length_transverse) =
            projectile_moments_of_inertia(mass, caliber, 2.0 * length);
        assert!((double_length_spin / spin - 1.0).abs() < 1e-12);
        assert!(double_length_transverse / transverse > 3.5);
    }

    #[test]
    fn projectile_inertias_reject_invalid_geometry() {
        let invalid = [0.0, -1.0, f64::NAN, f64::INFINITY];

        for value in invalid {
            assert_eq!(
                projectile_moments_of_inertia(value, 0.00782, 0.033),
                (0.0, 0.0)
            );
            assert_eq!(
                projectile_moments_of_inertia(0.01134, value, 0.033),
                (0.0, 0.0)
            );
            assert_eq!(
                projectile_moments_of_inertia(0.01134, 0.00782, value),
                (0.0, 0.0)
            );
        }

        assert_eq!(
            projectile_moments_of_inertia(f64::MAX, f64::MAX, f64::MAX),
            (0.0, 0.0)
        );
    }

    #[test]
    fn test_mba941_epicyclic_relations_and_limits() {
        // Validate the corrected frequencies against the EXACT algebraic relations of the standard
        // epicyclic decomposition (no external reference needed): with arm = Ix p / 2 Iy,
        //   fast + slow = Ix p / Iy = 2*arm     (sum of the arm rates)
        //   fast * slow = arm^2 / Sg            (product)
        let (ix, iy, p) = (6.94e-8_f64, 9.13e-7_f64, 17522.0_f64);
        let arm = ix * p / (2.0 * iy);
        for &sg in &[1.5_f64, 2.5, 5.0, 50.0] {
            let (fast, slow) = epicyclic_frequencies(ix, iy, p, sg);
            assert!(fast > slow && slow > 0.0, "expect fast>slow>0 at Sg={sg}");
            assert!(
                ((fast + slow) - 2.0 * arm).abs() < 1e-6 * arm,
                "sum != Ix p / Iy at Sg={sg}"
            );
            assert!(
                (fast * slow - arm * arm / sg).abs() < 1e-6 * arm * arm,
                "product != arm^2 / Sg at Sg={sg}"
            );
        }
        // Marginal stability (Sg -> 1+): the two modes coalesce at Ix p / 2 Iy.
        let (f1, s1) = epicyclic_frequencies(ix, iy, p, 1.0 + 1e-9);
        assert!((f1 - arm).abs() < 1e-3 * arm && (s1 - arm).abs() < 1e-3 * arm);
        // High stability: slow precession -> 0, fast nutation -> Ix p / Iy = 2*arm.
        let (f2, s2) = epicyclic_frequencies(ix, iy, p, 1.0e6);
        assert!(s2 < 1e-3 * arm, "slow precession should vanish at high Sg");
        assert!((f2 - 2.0 * arm).abs() < 1e-3 * arm, "fast -> Ix p / Iy at high Sg");
        // Not gyroscopically stable -> no epicyclic motion.
        assert_eq!(epicyclic_frequencies(ix, iy, p, 0.9), (0.0, 0.0));
    }

    #[test]
    fn test_precession_frequency() {
        // Slow (precession) mode, rad/s: (Ix p / 2 Iy)(1 - sqrt(1 - 1/Sg)).
        let freq = calculate_precession_frequency(17522.0, 6.94e-8, 9.13e-7, 2.5);
        let nut = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 2.5);
        // Positive, and slower than the nutation (fast) mode.
        assert!(
            freq > 0.0 && freq < nut,
            "precession {freq} should satisfy 0 < freq < nutation {nut}"
        );
        // Unstable -> no precession.
        assert_eq!(
            calculate_precession_frequency(17522.0, 6.94e-8, 9.13e-7, 0.9),
            0.0
        );
    }

    #[test]
    fn test_nutation_frequency() {
        // Fast (nutation) mode, rad/s: (Ix p / 2 Iy)(1 + sqrt(1 - 1/Sg)).
        // arm = Ix p / 2 Iy ~= 666 rad/s; fast = arm*(1 + sqrt(1/3)) ~= 1050 rad/s.
        let freq = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 1.5);
        assert!(
            (900.0..1200.0).contains(&freq),
            "nutation freq {freq} rad/s out of expected band"
        );
    }

    #[test]
    fn test_nutation_damping() {
        let initial = 0.01;
        let freq = 3000.0;

        // Check exponential decay
        let amp_0 = calculate_nutation_amplitude(initial, 0.0, freq, 0.05, 17522.0);
        let amp_1 = calculate_nutation_amplitude(initial, 0.1, freq, 0.05, 17522.0);

        assert_eq!(amp_0, initial);
        assert!(amp_1 < amp_0);
        assert!(amp_1 > 0.0);
    }

    #[test]
    fn test_precession_edge_cases() {
        // Unstable / marginally stable -> no regular precession.
        assert_eq!(
            calculate_precession_frequency(17522.0, 6.94e-8, 9.13e-7, 0.9),
            0.0
        );
        assert_eq!(
            calculate_precession_frequency(17522.0, 6.94e-8, 9.13e-7, 1.0),
            0.0
        );
        // Zero transverse inertia -> guarded to 0.
        assert_eq!(
            calculate_precession_frequency(17522.0, 6.94e-8, 0.0, 2.0),
            0.0
        );
    }

    #[test]
    fn test_nutation_edge_cases() {
        // Test unstable projectile (Sg <= 1)
        let freq_unstable = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 0.9);
        assert_eq!(freq_unstable, 0.0);

        // Test marginally stable (Sg = 1)
        let freq_marginal = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 1.0);
        assert_eq!(freq_marginal, 0.0);

        // Test zero transverse inertia
        let freq_zero_inertia = calculate_nutation_frequency(17522.0, 6.94e-8, 0.0, 2.0);
        assert_eq!(freq_zero_inertia, 0.0);
    }

    #[test]
    fn test_nutation_amplitude_bounds() {
        let initial = 0.5; // Large initial disturbance
        let freq = 3000.0;
        let spin = 17522.0;

        // Even with large initial disturbance, should be clamped
        let amp = calculate_nutation_amplitude(initial, 0.0, freq, 0.05, spin);
        assert!(amp <= 0.1); // Max 0.1 rad

        // Test zero frequency
        let amp_zero_freq = calculate_nutation_amplitude(initial, 1.0, 0.0, 0.05, spin);
        assert_eq!(amp_zero_freq, 0.0);

        // Test zero spin
        let amp_zero_spin = calculate_nutation_amplitude(initial, 1.0, freq, 0.05, 0.0);
        assert_eq!(amp_zero_spin, 0.0);
    }

    #[test]
    fn test_epicyclic_motion() {
        let (pitch, yaw) = calculate_epicyclic_motion(
            6.94e-8, // spin inertia
            9.13e-7, // transverse inertia
            17522.0, // spin rate
            2.5,     // stability factor
            0.1,     // time
            0.01,    // initial yaw
        );

        // Bounded by |initial_yaw| * (1 + 1/Sg) (slow arm amplitude 1 + fast arm amplitude 1/Sg).
        let bound = 0.01 * (1.0 + 1.0 / 2.5) + 1e-9;
        assert!(pitch.abs() <= bound, "pitch {pitch} exceeds bound {bound}");
        assert!(yaw.abs() <= bound, "yaw {yaw} exceeds bound {bound}");

        // Unstable case -> returns the initial yaw unchanged.
        let (pitch_unstable, yaw_unstable) =
            calculate_epicyclic_motion(6.94e-8, 9.13e-7, 17522.0, 0.9, 0.1, 0.01);
        assert_eq!(pitch_unstable, 0.01);
        assert_eq!(yaw_unstable, 0.01);
    }

    #[test]
    fn test_limit_cycle_yaw() {
        // Test with crosswind
        let yaw_wind = calculate_limit_cycle_yaw(
            850.0,   // velocity
            17522.0, // spin rate
            2.5,     // stability
            10.0,    // crosswind
        );

        // Should be small but non-zero
        assert!(yaw_wind > 0.0);
        assert!(yaw_wind < 0.1);

        // Test without crosswind
        let yaw_no_wind = calculate_limit_cycle_yaw(850.0, 17522.0, 2.5, 0.0);
        assert!(yaw_no_wind > 0.0);
        assert!(yaw_no_wind < yaw_wind);

        // Test unstable projectile
        let yaw_unstable = calculate_limit_cycle_yaw(850.0, 17522.0, 0.9, 0.0);
        assert_eq!(yaw_unstable, 0.01); // Fixed value for unstable
    }

    #[test]
    fn test_combined_angular_motion() {
        let params = PrecessionNutationParams::default();
        let initial_state = AngularState {
            pitch_angle: 0.001,
            yaw_angle: 0.002,
            pitch_rate: 0.01,
            yaw_rate: 0.01,
            precession_angle: 0.0,
            nutation_phase: 0.0,
        };

        let new_state = calculate_combined_angular_motion(
            &params,
            &initial_state,
            0.1,   // time
            0.001, // dt
            0.001, // initial disturbance
        );

        // Check that nutation phase evolved (it always should with non-zero frequency)
        // Precession might be very small with small yaw angles
        assert!(
            new_state.nutation_phase != initial_state.nutation_phase
                || new_state.precession_angle != initial_state.precession_angle
        );

        // Check reasonable bounds
        assert!(new_state.pitch_angle.abs() < 1.0);
        assert!(new_state.yaw_angle.abs() < 1.0);
    }

    #[test]
    fn combined_motion_applies_velocity_and_density_to_stability() {
        let velocity_mps = 1_000.0;
        let spin_rate_rad_s = 15_095.0;
        let params = PrecessionNutationParams {
            mass_kg: 0.01134,
            caliber_m: 0.00782,
            length_m: 0.033,
            spin_rate_rad_s,
            spin_inertia: 6.94e-8,
            transverse_inertia: 9.13e-7,
            velocity_mps,
            air_density_kg_m3: 1.0,
            mach: velocity_mps / 343.0,
            pitch_damping_coeff: PitchDampingCoefficients::default().subsonic,
            nutation_damping_factor: 0.05,
        };
        let initial_state = AngularState {
            pitch_angle: 0.0,
            yaw_angle: 0.0,
            pitch_rate: 0.0,
            yaw_rate: 0.0,
            precession_angle: 0.0,
            nutation_phase: 0.0,
        };

        let caliber_in = params.caliber_m / 0.0254;
        let length_in = params.length_m / 0.0254;
        let mass_gr = params.mass_kg / 0.00006479891;
        let spin_rps = spin_rate_rad_s / (2.0 * std::f64::consts::PI);
        let twist_in = velocity_mps * 3.28084 * 12.0 / spin_rps;
        let bare_sg = crate::spin_drift::miller_stability(
            caliber_in,
            mass_gr,
            twist_in,
            length_in,
        );
        let velocity_correction = (velocity_mps * 3.28084 / 2800.0).powf(1.0 / 3.0);
        let density_correction = 1.225 / params.air_density_kg_m3;
        let corrected_sg = crate::spin_drift::calculate_dynamic_stability(
            mass_gr,
            velocity_mps,
            spin_rate_rad_s,
            caliber_in,
            length_in,
            params.air_density_kg_m3,
        );
        assert!(bare_sg < 1.0);
        assert!(bare_sg * velocity_correction < 1.0);
        assert!(bare_sg * density_correction < 1.0);
        assert!(corrected_sg > 1.0);

        let dt = 0.0001;
        let actual = calculate_combined_angular_motion(&params, &initial_state, 0.0, dt, 0.001);
        let expected_precession = calculate_precession_frequency(
            spin_rate_rad_s,
            params.spin_inertia,
            params.transverse_inertia,
            corrected_sg,
        ) * dt;
        let expected_nutation = calculate_nutation_frequency(
            spin_rate_rad_s,
            params.spin_inertia,
            params.transverse_inertia,
            corrected_sg,
        ) * dt;

        assert!(
            (actual.precession_angle - expected_precession).abs() < 1e-12,
            "corrected precession phase mismatch: actual={} expected={expected_precession}",
            actual.precession_angle
        );
        assert!(
            (actual.nutation_phase - expected_nutation).abs() < 1e-12,
            "corrected nutation phase mismatch: actual={} expected={expected_nutation}",
            actual.nutation_phase
        );
    }

    #[test]
    fn test_default_params() {
        let params = PrecessionNutationParams::default();

        // Check reasonable default values
        assert!(params.mass_kg > 0.0);
        assert!(params.caliber_m > 0.0);
        assert!(params.length_m > 0.0);
        assert!(params.spin_rate_rad_s > 0.0);
        assert!(params.spin_inertia > 0.0);
        assert!(params.transverse_inertia > 0.0);
        assert!(params.velocity_mps > 0.0);
        assert!(params.air_density_kg_m3 > 0.0);
        assert!(params.mach > 0.0);
        assert!(params.nutation_damping_factor > 0.0);
        assert!(params.nutation_damping_factor < 1.0); // Should be fraction
    }

    #[test]
    fn test_stability_effects() {
        // Higher stability gives a higher nutation (fast-mode) frequency:
        // phi_fast = (Ix p / 2 Iy)(1 + sqrt(1 - 1/Sg)) increases monotonically with Sg.
        let freq_high_stability = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 5.0);
        let freq_low_stability = calculate_nutation_frequency(17522.0, 6.94e-8, 9.13e-7, 1.5);
        assert!(freq_high_stability > freq_low_stability);
    }

    #[test]
    fn test_damping_time_evolution() {
        let initial = 0.01;
        let freq = 3000.0;
        let spin = 17522.0;
        let damping = 0.1;

        // Sample at different times
        let times = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2];
        let mut last_amp = initial;

        for &t in &times[1..] {
            let amp = calculate_nutation_amplitude(initial, t, freq, damping, spin);

            // Should monotonically decrease
            assert!(amp < last_amp);
            assert!(amp >= 0.0);
            last_amp = amp;
        }
    }

    #[test]
    fn test_angular_state_evolution() {
        let params = PrecessionNutationParams {
            mass_kg: 0.01,
            caliber_m: 0.008,
            length_m: 0.03,
            spin_rate_rad_s: 16000.0, // fast enough that the Miller Sg > 1 (gyroscopically stable)
            spin_inertia: 5e-8,
            transverse_inertia: 8e-7,
            velocity_mps: 800.0,
            air_density_kg_m3: 1.2,
            mach: 2.3,
            pitch_damping_coeff: -5.0,
            nutation_damping_factor: 0.08,
        };

        let mut state = AngularState {
            pitch_angle: 0.0,
            yaw_angle: 0.005,
            pitch_rate: 0.0,
            yaw_rate: 0.0,
            precession_angle: 0.0,
            nutation_phase: 0.0,
        };

        // Store initial state for comparison
        let initial_phase = state.nutation_phase;
        let initial_precession = state.precession_angle;

        // Evolve for several timesteps
        let dt = 0.0001;
        for i in 0..100 {
            let time = i as f64 * dt;
            state = calculate_combined_angular_motion(&params, &state, time, dt, 0.002);
        }

        // Should have evolved - at least one of these should change
        assert!(
            state.precession_angle != initial_precession || state.nutation_phase != initial_phase
        );

        // Should remain bounded
        assert!(state.yaw_angle.abs() < 0.1);
        assert!(state.pitch_angle.abs() < 0.1);
    }
}
