//! Precession and Nutation Physics for Ballistic Projectiles
//!
//! This module implements the complex angular motion of spinning projectiles:
//! - Precession: Slow coning motion of the projectile axis
//! - Nutation: Fast oscillatory motion superimposed on precession
//! - Angular momentum conservation
//! - Gyroscopic effects

// Precession and nutation modeling - now integrated!

use crate::pitch_damping::{
    calculate_pitch_damping_moment,
    PitchDampingCoefficients,
};

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
    pub spin_inertia: f64,      // About longitudinal axis
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
            mass_kg: 0.01134,  // 175 grains
            caliber_m: 0.00782,
            length_m: 0.033,
            spin_rate_rad_s: 17522.0,
            spin_inertia: 6.94e-8,
            transverse_inertia: 9.13e-7,
            velocity_mps: 850.0,
            air_density_kg_m3: 1.225,
            mach: 2.48,
            pitch_damping_coeff: -0.8,
            nutation_damping_factor: 0.05,
        }
    }
}

/// Calculate the natural precession frequency
pub fn calculate_precession_frequency(
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
    yaw_angle_rad: f64,
) -> f64 {
    if velocity_mps == 0.0 || transverse_inertia == 0.0 {
        return 0.0;
    }
    
    // Basic gyroscopic precession
    // ωp = (Is * ωs * sin(α)) / (It * V)
    (spin_inertia * spin_rate_rad_s * yaw_angle_rad.sin()) / 
        (transverse_inertia * velocity_mps)
}

/// Calculate the natural nutation frequency
pub fn calculate_nutation_frequency(
    spin_rate_rad_s: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
    stability_factor: f64,
) -> f64 {
    if stability_factor <= 1.0 || transverse_inertia == 0.0 {
        return 0.0;
    }
    
    // Nutation frequency
    // ωn = ωs * sqrt(Is / It) * sqrt(Sg - 1)
    let inertia_ratio = spin_inertia / transverse_inertia;
    spin_rate_rad_s * inertia_ratio.sqrt() * (stability_factor - 1.0).sqrt()
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
    // Calculate stability factor (simplified)
    let stability = (params.spin_inertia * params.spin_rate_rad_s.powi(2)) / 
        (4.0 * params.transverse_inertia * params.velocity_mps.powi(2) / params.length_m);
    
    // Precession frequency
    let omega_p = calculate_precession_frequency(
        params.spin_rate_rad_s,
        params.velocity_mps,
        params.spin_inertia,
        params.transverse_inertia,
        angular_state.yaw_angle,
    );
    
    // Nutation frequency
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
    let new_yaw_rate = omega_p; // Precession rate
    
    // Combined angle with nutation
    // The total yaw is precession + nutation oscillation
    let total_yaw = angular_state.yaw_angle + nutation_amp * new_nutation_phase.sin();
    
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
    spin_rate_rad_s: f64,
    velocity_mps: f64,
    stability_factor: f64,
    time_s: f64,
    initial_yaw_rad: f64,
) -> (f64, f64) {
    if stability_factor <= 1.0 {
        // Unstable - no regular motion
        return (initial_yaw_rad, initial_yaw_rad);
    }
    
    // Frequencies (simplified model)
    // Slow mode (precession)
    let omega_slow = 2.0 * velocity_mps / (stability_factor * spin_rate_rad_s);
    
    // Fast mode (nutation)
    let omega_fast = spin_rate_rad_s * ((stability_factor - 1.0).sqrt()) / stability_factor;
    
    // Amplitude ratio (fast/slow)
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
    fn test_precession_frequency() {
        let freq = calculate_precession_frequency(
            17522.0,  // spin rate
            850.0,    // velocity
            6.94e-8,  // spin inertia
            9.13e-7,  // transverse inertia
            0.002,    // yaw angle
        );
        
        // Should be small for small yaw angles
        assert!(freq.abs() < 1.0);
    }
    
    #[test]
    fn test_nutation_frequency() {
        let freq = calculate_nutation_frequency(
            17522.0,  // spin rate
            6.94e-8,  // spin inertia
            9.13e-7,  // transverse inertia
            1.5,      // stability
        );
        
        // Should be in the kHz range
        assert!(freq > 1000.0);
        assert!(freq < 10000.0);
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
}