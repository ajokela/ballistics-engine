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
