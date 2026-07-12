#![no_main]
use arbitrary::Unstructured;
use ballistics_engine::ffi::{
    ballistics_calculate_trajectory, ballistics_free_trajectory_result, FFIAtmosphericConditions,
    FFIBallisticInputs, FFIWindConditions,
};
use ballistics_engine_fuzz::domain::{hostile_inputs, ranged};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    // Exercise hostile numeric values through the C conversion and alloc/free boundary. Invalid
    // inputs must return null cleanly; successful results must remain safely releasable.
    let Ok(vi) = hostile_inputs(&mut u) else {
        return;
    };
    let Ok(sight) = ranged(&mut u, -0.2, 0.2) else {
        return;
    }; // finite hostile sight height
    let ffi_in = FFIBallisticInputs {
        muzzle_velocity: vi.muzzle_velocity,
        muzzle_angle: vi.muzzle_angle,
        bc_value: vi.bc_value,
        bullet_mass: vi.bullet_mass,
        bullet_diameter: vi.bullet_diameter,
        bc_type: 1,
        sight_height: sight,
        target_distance: vi.target_distance,
        temperature: vi.temperature,
        twist_rate: vi.twist_rate,
        is_twist_right: 1,
        shooting_angle: 0.0,
        altitude: vi.altitude,
        latitude: f64::NAN, // NAN = "unset" per FFI contract
        azimuth_angle: 0.0,
        use_rk4: 1,
        use_adaptive_rk45: 0,
        enable_wind_shear: 0,
        enable_trajectory_sampling: 0,
        sample_interval: 10.0,
        enable_pitch_damping: 0,
        enable_precession_nutation: 0,
        enable_spin_drift: 0,
        enable_magnus: 0,
        enable_coriolis: 0,
        shot_azimuth: 0.0,
    };
    // FFI step_size is internally divided by 1000 (dt = step_size/1000 s), so the
    // README/tests use [0.1, 1.0] (dt 1e-4..1e-3 s). Fuzz within that sane band.
    // Finding #2 (logged): a far smaller step_size drives unbounded trajectory-point
    // allocation / OOM through this ABI — tracked separately, not re-triggered here.
    let Ok(step) = ranged(&mut u, 0.1, 1.0) else {
        return;
    };
    // SAFETY: valid pointers to locals; result freed via the paired free fn.
    unsafe {
        let wind = FFIWindConditions {
            speed: 0.0,
            direction: 0.0,
        };
        let atmo = FFIAtmosphericConditions {
            temperature: 15.0,
            pressure: 1013.25,
            humidity: 50.0,
            altitude: 0.0,
        };
        let res = ballistics_calculate_trajectory(&ffi_in, &wind, &atmo, 3000.0, step);
        if !res.is_null() {
            ballistics_free_trajectory_result(res);
        }
    }
});
