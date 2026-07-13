//! Regressions frozen from fuzzing findings.
//!
//! When a harness surfaces a real bug:
//!   1. minimize it:  `cd fuzz && cargo +nightly fuzz tmin <target> <artifact>`
//!   2. decode the minimized input the same way the harness does (domain::*),
//!   3. add a `#[test]` here that reproduces the *fixed* behavior so it can never
//!      silently return. Reference the MBA ticket in the test name/comment.
//!
//! This file starts with only a smoke test so `cargo test` stays green.

use ballistics_engine::{AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions};

#[test]
fn baseline_solve_is_finite() {
    // A plain, valid shot must always produce a finite trajectory — the invariant
    // every fuzz finding is a specialization of.
    let mut b = BallisticInputs::default();
    b.bc_value = 0.5;
    b.bullet_mass = 0.01;
    b.bullet_diameter = 0.0077;
    b.bullet_length = 0.03;
    b.muzzle_velocity = 800.0;
    b.muzzle_angle = 0.02;
    b.target_distance = 500.0;
    b.twist_rate = 10.0;
    b.use_rk4 = true;

    let r = TrajectorySolver::new(b, WindConditions::default(), AtmosphericConditions::default())
        .solve()
        .expect("baseline solve should succeed");
    assert!(r.time_of_flight.is_finite() && r.time_of_flight > 0.0);
    assert!(r.impact_velocity.is_finite() && r.impact_velocity >= 0.0);
    for p in &r.points {
        assert!(p.position.x.is_finite() && p.position.y.is_finite() && p.position.z.is_finite());
    }
}

// ── frozen findings go below, one #[test] each ──

#[test]
fn mba1282_invalid_core_inputs_are_rejected_in_every_solver() {
    type Setter = fn(&mut BallisticInputs, f64);
    let fields: [(&str, Setter); 4] = [
        ("bc_value", |inputs, value| inputs.bc_value = value),
        ("bullet_mass", |inputs, value| inputs.bullet_mass = value),
        ("bullet_diameter", |inputs, value| {
            inputs.bullet_diameter = value
        }),
        ("muzzle_velocity", |inputs, value| {
            inputs.muzzle_velocity = value
        }),
    ];
    let bad_values = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0];
    let modes = [
        ("Euler", false, false),
        ("RK4", true, false),
        ("RK45", true, true),
    ];

    for (mode, use_rk4, use_adaptive_rk45) in modes {
        for (field, set) in fields {
            for value in bad_values {
                let mut inputs = BallisticInputs {
                    use_rk4,
                    use_adaptive_rk45,
                    ..BallisticInputs::default()
                };
                set(&mut inputs, value);

                let error = TrajectorySolver::new(
                    inputs,
                    WindConditions::default(),
                    AtmosphericConditions::default(),
                )
                .solve()
                .expect_err("invalid core input must fail before integration");
                assert!(
                    error.to_string().contains(field),
                    "{mode} error for {field}={value:?} did not name the field: {error}"
                );
            }
        }
    }
}

#[test]
fn mba1282_non_finite_launch_angle_is_rejected() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let inputs = BallisticInputs {
            muzzle_angle: value,
            ..BallisticInputs::default()
        };
        let error = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        )
        .solve()
        .expect_err("non-finite muzzle angle must fail before integration");
        assert!(error.to_string().contains("muzzle_angle"), "{error}");
    }
}

#[test]
fn mba1282_powder_curve_override_is_validated_after_resolution() {
    let valid_override = BallisticInputs {
        muzzle_velocity: f64::NAN,
        powder_temp_curve: Some(vec![(0.0, 700.0), (30.0, 800.0)]),
        ..BallisticInputs::default()
    };
    TrajectorySolver::new(
        valid_override,
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
    .solve()
    .expect("a valid absolute powder-curve override replaces the raw velocity");

    let invalid_override = BallisticInputs {
        powder_temp_curve: Some(vec![(0.0, f64::INFINITY), (30.0, f64::INFINITY)]),
        ..BallisticInputs::default()
    };
    let error = TrajectorySolver::new(
        invalid_override,
        WindConditions::default(),
        AtmosphericConditions::default(),
    )
    .solve()
    .expect_err("the resolved powder-curve velocity must be validated");
    assert!(error.to_string().contains("muzzle_velocity"), "{error}");
}

#[test]
fn mba1282_finite_overflow_is_returned_as_an_error() {
    let inputs = BallisticInputs {
        bullet_mass: 1.0e304,
        use_rk4: false,
        ..BallisticInputs::default()
    };
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_max_range(0.1);
    solver.set_time_step(0.001);

    let error = solver
        .solve()
        .expect_err("finite inputs that overflow must not produce Ok(Inf)");
    assert!(
        error.to_string().contains("non-finite"),
        "unexpected overflow error: {error}"
    );
}

#[test]
fn mba1282_non_finite_segment_state_is_returned_as_an_error() {
    for (mode, use_rk4, use_adaptive_rk45) in [
        ("Euler", false, false),
        ("RK4", true, false),
        ("RK45", true, true),
    ] {
        let inputs = BallisticInputs {
            use_rk4,
            use_adaptive_rk45,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_wind_segments(vec![ballistics_engine::wind::WindSegment::new(
            f64::NAN, 90.0, 100.0,
        )]);

        let error = solver
            .solve()
            .expect_err("a non-finite segmented wind must not leave a stale finite result");
        assert!(
            error.to_string().contains("non-finite state"),
            "unexpected {mode} segmented-wind error: {error}"
        );
    }
}
