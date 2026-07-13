//! Integration tests for downrange-segmented wind on the `TrajectorySolver`
//! (`set_wind_segments`). Verifies that segmented wind actually varies the wind
//! along the path, that a single full-range segment matches the equivalent
//! scalar wind, and that the wind-FROM angle convention is correct.

use ballistics_engine::{
    AtmosphericConditions, BallisticInputs, TrajectorySolver, WindConditions,
};
use ballistics_engine::wind::WindSegment;
use std::f64::consts::PI;

fn base_inputs() -> BallisticInputs {
    let mut inputs = BallisticInputs::default();
    inputs.muzzle_velocity = 792.48; // 2600 fps
    inputs.bc_value = 0.243;
    inputs.bc_type = ballistics_engine::DragModel::G7;
    inputs.bullet_mass = 0.01134; // ~175 gr
    inputs.bullet_diameter = 0.007823; // .308"
    inputs.muzzle_angle = 0.005; // slight elevation so it flies a while
    inputs
}

/// Final lateral drift (McCoy Z) at ~700 m for the given wind setup.
fn drift_z(wind: WindConditions, segments: Vec<WindSegment>) -> f64 {
    let mut solver = TrajectorySolver::new(base_inputs(), wind, AtmosphericConditions::default());
    solver.set_max_range(800.0);
    if !segments.is_empty() {
        solver.set_wind_segments(segments);
    }
    let r = solver.solve().unwrap();
    r.position_at_range(700.0).unwrap().z
}

#[test]
fn single_full_range_segment_matches_scalar_wind() {
    // 10 mph (4.4704 m/s = 16.09344 km/h) from the right (+90deg), covering the
    // whole flight, must match the equivalent scalar wind.
    let scalar = drift_z(
        WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
        },
        vec![],
    );
    let segmented = drift_z(
        WindConditions::default(),
        vec![WindSegment::new(16.09344, 90.0, 5000.0)],
    );
    assert!(
        (scalar - segmented).abs() < 1e-6,
        "single full-range segment {segmented} != scalar wind {scalar}"
    );
    // Sanity: a right wind drifts the bullet left (negative Z).
    assert!(scalar < 0.0, "wind from the right should drift left, got {scalar}");
}

#[test]
fn segmented_wind_varies_along_path() {
    // Wind only over the first ~half of the flight produces drift strictly
    // between the no-wind (zero) and full-wind cases.
    let full = drift_z(
        WindConditions::default(),
        vec![WindSegment::new(16.09344, 90.0, 5000.0)],
    );
    let none = drift_z(WindConditions::default(), vec![]);
    let half = drift_z(
        WindConditions::default(),
        vec![WindSegment::new(16.09344, 90.0, 400.0)],
    );

    assert!(none.abs() < 1e-9, "no segments -> no wind, got {none}");
    // full and half are negative (drift left); half is smaller in magnitude.
    assert!(full < 0.0 && half < 0.0);
    assert!(
        half.abs() < full.abs(),
        "partial-coverage drift {half} should be smaller than full-coverage {full}"
    );
}

#[test]
fn segment_angle_convention_matches_scalar() {
    // From the right (+90) drifts left (z<0); from the left (270) drifts right (z>0).
    let from_right = drift_z(
        WindConditions::default(),
        vec![WindSegment::new(16.09344, 90.0, 5000.0)],
    );
    let from_left = drift_z(
        WindConditions::default(),
        vec![WindSegment::new(16.09344, 270.0, 5000.0)],
    );
    assert!(from_right < 0.0, "90deg (from right) -> drift left, got {from_right}");
    assert!(from_left > 0.0, "270deg (from left) -> drift right, got {from_left}");
    assert!(
        (from_right + from_left).abs() < 1e-6,
        "left/right crosswind drift should be symmetric ({from_right}, {from_left})"
    );
}

#[test]
fn empty_segments_revert_to_scalar() {
    // Passing an empty segment list clears segmented wind (uses scalar wind).
    let mut solver = TrajectorySolver::new(
        base_inputs(),
        WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
        },
        AtmosphericConditions::default(),
    );
    solver.set_max_range(800.0);
    solver.set_wind_segments(vec![WindSegment::new(99.0, 90.0, 5000.0)]);
    solver.set_wind_segments(vec![]); // clear
    let z = solver.solve().unwrap().position_at_range(700.0).unwrap().z;
    let scalar = drift_z(
        WindConditions {
            speed: 4.4704,
            direction: PI / 2.0,
        },
        vec![],
    );
    assert!((z - scalar).abs() < 1e-6, "cleared segments {z} != scalar {scalar}");
}
