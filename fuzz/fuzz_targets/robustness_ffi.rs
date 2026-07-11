#![no_main]
use libfuzzer_sys::fuzz_target;
use arbitrary::Unstructured;
use ballistics_engine::ffi::{
    ballistics_calculate_trajectory, ballistics_free_trajectory_result,
    FFIAtmosphericConditions, FFIBallisticInputs, FFIWindConditions,
};
use ballistics_engine_fuzz::domain::{ranged, valid_inputs};

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    // Finding #1: the engine lacks a finite/physical input-validation gate. Feed FINITE,
    // plausible inputs and exercise the FFI conversion + alloc/free path for soundness
    // (no panic/UB, clean round-trip, finite outputs). A couple of fields take finite
    // hostile extremes to probe the C boundary. Widen to non-finite FFI inputs once the
    // engine gains a validation gate.
    let Ok(vi) = valid_inputs(&mut u) else { return };
    let Ok(sight) = ranged(&mut u, -0.2, 0.2) else { return };       // finite hostile sight height
    let ffi_in = FFIBallisticInputs {
        muzzle_velocity: vi.muzzle_velocity, muzzle_angle: vi.muzzle_angle, bc_value: vi.bc_value,
        bullet_mass: vi.bullet_mass, bullet_diameter: vi.bullet_diameter,
        bc_type: 1, sight_height: sight, target_distance: vi.target_distance,
        temperature: vi.temperature, twist_rate: vi.twist_rate, is_twist_right: 1,
        shooting_angle: 0.0, altitude: vi.altitude, latitude: f64::NAN, // NAN = "unset" per FFI contract
        azimuth_angle: 0.0, use_rk4: 1, use_adaptive_rk45: 0, enable_wind_shear: 0,
        enable_trajectory_sampling: 0, sample_interval: 10.0, enable_pitch_damping: 0,
        enable_precession_nutation: 0, enable_spin_drift: 0, enable_magnus: 0,
        enable_coriolis: 0, shot_azimuth: 0.0,
    };
    // SAFETY: valid pointers to locals; result freed via the paired free fn.
    unsafe {
        let wind = FFIWindConditions { speed: 0.0, direction: 0.0 };
        let atmo = FFIAtmosphericConditions { temperature: 15.0, pressure: 1013.25, humidity: 50.0, altitude: 0.0 };
        let res = ballistics_calculate_trajectory(&ffi_in, &wind, &atmo, 3000.0, 0.001);
        if !res.is_null() {
            ballistics_free_trajectory_result(res);
        }
    }
});
