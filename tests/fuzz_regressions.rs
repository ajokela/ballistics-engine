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
    let b = BallisticInputs {
        bc_value: 0.5,
        bullet_mass: 0.01,
        bullet_diameter: 0.0077,
        bullet_length: 0.03,
        muzzle_velocity: 800.0,
        muzzle_angle: 0.02,
        target_distance: 500.0,
        twist_rate: 10.0,
        use_rk4: true,
        ..Default::default()
    };

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
fn mba1293_stiff_drag_explosion_is_an_error_not_a_negative_range() {
    // Decoded from the minimized robustness-inputs artifact (fuzz/tests/crash-mba-1293.bin):
    // every field passes the positivity/finiteness gate, but the drag time constant is
    // astronomically shorter than the minimum integration step. The accepted minimum-dt RK45
    // step multiplied speed 13x AND reversed it, and the solve returned
    // Ok(max_range = -50.588...) — the bullet reported 50 m behind the shooter.
    for (mode, use_rk4, use_adaptive_rk45) in [
        ("Euler", false, false),
        ("RK4", true, false),
        ("RK45", true, true),
    ] {
        let inputs = BallisticInputs {
            bc_value: 6.821210274138984e-2,
            bullet_mass: 3.320973714192708e8,      // kg
            bullet_diameter: 3.333333333333333e8,  // m
            muzzle_velocity: 3.333435058584939e8,  // m/s
            muzzle_angle: 0.0,
            target_distance: 10.0,
            use_rk4,
            use_adaptive_rk45,
            ..BallisticInputs::default()
        };
        let mut solver = TrajectorySolver::new(
            inputs,
            WindConditions::default(),
            AtmosphericConditions::default(),
        );
        solver.set_max_range(10.0);

        match solver.solve() {
            Err(_) => {} // clean rejection is the required outcome
            Ok(result) => panic!(
                "{mode}: stiff-input explosion must be an Err, got Ok with max_range {}",
                result.max_range
            ),
        }
    }
}

#[test]
fn mba1293_sane_solves_stay_within_the_speed_budget() {
    // The divergence guard must never fire on legitimate trajectories, including
    // a strong segmented wind (the budget accounts for the strongest segment).
    let inputs = BallisticInputs {
        bc_value: 0.5,
        bullet_mass: 0.01,
        bullet_diameter: 0.0077,
        bullet_length: 0.03,
        muzzle_velocity: 900.0,
        muzzle_angle: 0.03,
        target_distance: 1000.0,
        twist_rate: 10.0,
        use_rk4: true,
        use_adaptive_rk45: true,
        ..Default::default()
    };
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions::default(),
        AtmosphericConditions::default(),
    );
    solver.set_wind_segments(vec![
        ballistics_engine::wind::WindSegment::new(40.0, 90.0, 400.0),
        ballistics_engine::wind::WindSegment::new(60.0, 270.0, 1000.0),
    ]);
    let result = solver.solve().expect("windy 1000 m solve must succeed");
    assert!(result.max_range > 990.0 && result.time_of_flight > 0.5);
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
