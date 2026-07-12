#![no_main]
use arbitrary::Unstructured;
use ballistics_engine::ffi::{
    ballistics_calculate_trajectory, ballistics_free_trajectory_result, FFIAtmosphericConditions,
    FFIBallisticInputs, FFIWindConditions, MIN_FFI_STEP_SIZE_MS,
};
use ballistics_engine::MAX_TRAJECTORY_POINTS;
use ballistics_engine_fuzz::domain::{hostile_inputs, ranged, valid_inputs};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    // Keep half the cases physically valid so hostile step-size values reliably reach the FFI
    // preflight; use hostile physics in the other half to retain conversion/validation coverage.
    let Ok(use_valid_physics) = u.arbitrary::<bool>() else {
        return;
    };
    let vi = if use_valid_physics {
        let Ok(inputs) = valid_inputs(&mut u) else {
            return;
        };
        inputs
    } else {
        let Ok(inputs) = hostile_inputs(&mut u) else {
            return;
        };
        inputs
    };
    let Ok(sight) = ranged(&mut u, -0.2, 0.2) else {
        return;
    }; // finite hostile sight height
    let Ok(mode) = u.int_in_range(0u8..=2) else {
        return;
    };
    let (use_rk4, use_adaptive_rk45) = match mode {
        0 => (0, 0),
        1 => (1, 0),
        _ => (1, 1),
    };
    // Keep accepted minimum-step cases cheap enough for sustained sanitizer fuzzing. The exact
    // 250,000-point exhaustion behavior is covered deterministically in the solver unit tests.
    let max_range = vi.target_distance.min(300.0);
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
        use_rk4,
        use_adaptive_rk45,
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
    // Mix the exact invalid boundary classes with the documented valid band. `step_size` is
    // milliseconds; invalid/sub-minimum values must return null without starting integration.
    let Ok(step_class) = u.int_in_range(0u8..=9) else {
        return;
    };
    let step = match step_class {
        0 => f64::NAN,
        1 => f64::INFINITY,
        2 => f64::NEG_INFINITY,
        3 => 0.0,
        4 => -0.0,
        5 => -1.0,
        6 => MIN_FFI_STEP_SIZE_MS / 10.0,
        7 => MIN_FFI_STEP_SIZE_MS - 0.001,
        8 => MIN_FFI_STEP_SIZE_MS,
        _ => {
            let Ok(step) = ranged(&mut u, MIN_FFI_STEP_SIZE_MS, 1.0) else {
                return;
            };
            step
        }
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
        let res = ballistics_calculate_trajectory(&ffi_in, &wind, &atmo, max_range, step);
        let invalid_step = !step.is_finite() || step < MIN_FFI_STEP_SIZE_MS;
        if invalid_step {
            assert!(res.is_null(), "invalid FFI step_size returned a result");
        } else if !res.is_null() {
            assert!((*res).point_count >= 0, "negative FFI point count");
            assert!(
                (*res).point_count as usize <= MAX_TRAJECTORY_POINTS,
                "FFI result exceeded the trajectory-point budget"
            );
            ballistics_free_trajectory_result(res);
        }
    }
});
