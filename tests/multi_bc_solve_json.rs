//! Multi-BC support via `projectile.bc_segments` on the solve-json v1 wire.
//!
//! `bridge_call()` accepts and returns JSON through the solve-json v1 schema. Before this
//! feature the only way to reach velocity-segmented BCs over that path was the filesystem-
//! bound `corrections.bc5d_table_path`, which reads a proprietary binary `.bin` file.
//! `projectile.bc_segments` makes the same engine slot reachable from plain JSON data.
//!
//! This file pins the full wire contract:
//!
//!   (a) Absent `bc_segments` — no field at all — changes nothing: the solve is
//!       byte-identical to every prior request.
//!   (b) A two-segment ladder produces a trajectory that differs from the single-BC baseline,
//!       proving `use_bc_segments = true` was actually installed.
//!   (c) The resolved request echoes `bc_segments` verbatim so the perturbation / round-trip
//!       kernels re-apply the same ladder; a rebuilt request re-solves identically.
//!   (d) Mutual exclusion with `projectile.drag_table` → `ConflictingFields` error.
//!   (e) Mutual exclusion with `corrections.bc5d_table_path` → `ConflictingFields` error.
//!   (f) Validation errors name the exact offending path.
//!   (g) Bridge `solve` command accepts `bc_segments` and echoes it in the result.
//!   (h) A Rust caller that builds `ProjectileV1` directly (bypassing the decoder) is still
//!       refused by the service when the segments are invalid.

use ballistics_engine::solve_json::{
    decode_solve_request_v1, BcSegmentV1, ProjectileV1, SolveErrorCodeV1, SolveSuccessV1,
};
use ballistics_engine::solve_v1;
use serde_json::{json, Value};

// ---------------------------------------------------------------------------
// Shared fixture helpers
// ---------------------------------------------------------------------------

/// A .308 175 gr at 792 m/s, zeroed at 100 m, solved to 800 m in still ICAO air.
fn base_request() -> Value {
    json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 0.01134,
            "diameter_m": 0.00782,
            "drag_model": "G7",
            "ballistic_coefficient": 0.243
        },
        "rifle": {
            "muzzle_velocity_mps": 792.48,
            "sight_height_m": 0.0508
        },
        "shot": {
            "max_range_m": 800.0,
            "zero_distance_m": 100.0
        },
        "atmosphere": {
            "altitude_m": 0.0,
            "temperature_k": 288.15,
            "pressure_pa": 101325.0,
            "relative_humidity": 0.0
        },
        "wind": {"speed_mps": 0.0, "direction_from_rad": 0.0},
        "solver": {"method": "rk4", "time_step_s": 0.0005},
        "effects": {},
        "sampling": {"interval_m": 100.0}
    })
}

/// Two-segment ladder: high-velocity band uses BC 0.260, low-velocity band uses 0.220.
/// Both differ from the scalar BC 0.243, so the trajectory MUST differ from the baseline.
fn two_segment_ladder() -> Value {
    json!([
        {"velocity_from_mps": 274.3, "velocity_to_mps": 999.0, "ballistic_coefficient": 0.260},
        {"velocity_from_mps": 0.0,   "velocity_to_mps": 274.3, "ballistic_coefficient": 0.220}
    ])
}

fn with_bc_segments(mut request: Value, segments: Value) -> Value {
    request["projectile"]["bc_segments"] = segments;
    request
}

fn decode(
    request: &Value,
) -> Result<ballistics_engine::SolveRequestV1, ballistics_engine::SolveErrorEnvelopeV1> {
    decode_solve_request_v1(&request.to_string())
}

fn solve(request: &Value) -> SolveSuccessV1 {
    let decoded = decode(request).expect("request should decode without error");
    solve_v1(decoded).expect("solve should succeed")
}

fn expect_error(request: &Value, code: SolveErrorCodeV1, path: &str) {
    let json_str = request.to_string();
    let result = decode_solve_request_v1(&json_str);
    let envelope = match result {
        Err(e) => e,
        Ok(typed) => solve_v1(typed).expect_err("service should have rejected the request"),
    };
    assert_eq!(
        envelope.error.code, code,
        "expected error code {:?} but got {:?}: {}",
        code, envelope.error.code, envelope.error.message
    );
    assert_eq!(
        envelope.error.path().map(|p| p.to_owned()),
        Some(path.to_owned()),
        "expected error path {path:?} but got {:?}",
        envelope.error.path()
    );
}

// ---------------------------------------------------------------------------
// (a) Backward compatibility — absent field is a transparent no-op
// ---------------------------------------------------------------------------

#[test]
fn absent_bc_segments_is_byte_identical_to_baseline() {
    let baseline = solve(&base_request());
    let mut with_none = base_request();
    with_none["projectile"]
        .as_object_mut()
        .unwrap()
        .remove("bc_segments");
    let result = solve(&with_none);

    assert_eq!(
        baseline.samples, result.samples,
        "absent bc_segments must not change the trajectory"
    );
    assert!(
        result.resolved_request.projectile.bc_segments.is_none(),
        "bc_segments must be absent from the resolved request when not supplied"
    );
}

// ---------------------------------------------------------------------------
// (b) Two-segment ladder produces a different trajectory
// ---------------------------------------------------------------------------

#[test]
fn two_segment_ladder_produces_different_drops_than_scalar_bc() {
    let baseline = solve(&base_request());
    let with_segments = solve(&with_bc_segments(base_request(), two_segment_ladder()));

    let drops_differ = baseline
        .samples
        .iter()
        .zip(with_segments.samples.iter())
        .any(|(b, s)| (b.drop_m - s.drop_m).abs() > 1e-6);

    assert!(
        drops_differ,
        "a two-segment BC ladder must produce drops that differ from the scalar-BC baseline;\n\
         baseline last drop: {:.6} m, segmented last drop: {:.6} m",
        baseline.samples.last().map(|s| s.drop_m).unwrap_or(f64::NAN),
        with_segments.samples.last().map(|s| s.drop_m).unwrap_or(f64::NAN),
    );
}

// ---------------------------------------------------------------------------
// (c) Resolved request echoes bc_segments; round-trip re-solves identically
// ---------------------------------------------------------------------------

#[test]
fn resolved_request_echoes_bc_segments() {
    let request = with_bc_segments(base_request(), two_segment_ladder());
    let result = solve(&request);

    let echo = result
        .resolved_request
        .projectile
        .bc_segments
        .as_ref()
        .expect("resolved_request must echo bc_segments when supplied");

    assert_eq!(echo.len(), 2);
    assert!((echo[0].velocity_from_mps - 274.3).abs() < 1e-9);
    assert!((echo[0].velocity_to_mps - 999.0).abs() < 1e-9);
    assert!((echo[0].ballistic_coefficient - 0.260).abs() < 1e-9);
    assert!((echo[1].velocity_from_mps - 0.0).abs() < 1e-9);
    assert!((echo[1].velocity_to_mps - 274.3).abs() < 1e-9);
    assert!((echo[1].ballistic_coefficient - 0.220).abs() < 1e-9);
}

#[test]
fn round_trip_resolves_identically() {
    use ballistics_engine::solve_json::SolveRequestV1;

    let request = with_bc_segments(base_request(), two_segment_ladder());
    let first = solve(&request);

    let rebuilt: SolveRequestV1 = (&first.resolved_request).into();
    assert!(
        rebuilt.projectile.bc_segments.is_some(),
        "round-tripped request must carry bc_segments"
    );

    let second = solve_v1(rebuilt).expect("round-tripped request must solve");

    assert_eq!(
        first.samples, second.samples,
        "a round-tripped bc_segments request must re-solve identically"
    );
}

// ---------------------------------------------------------------------------
// (d) Mutual exclusion with projectile.drag_table
// ---------------------------------------------------------------------------

#[test]
fn bc_segments_conflicts_with_drag_table() {
    let mut request = with_bc_segments(base_request(), two_segment_ladder());
    request["projectile"]["drag_table"] = json!({
        "kind": "reference",
        "points": [
            {"mach": 0.0, "cd": 0.2},
            {"mach": 2.0, "cd": 0.3}
        ]
    });
    expect_error(
        &request,
        SolveErrorCodeV1::ConflictingFields,
        "$.projectile.bc_segments",
    );
}

// ---------------------------------------------------------------------------
// (e) Mutual exclusion with corrections.bc5d_table_path
// ---------------------------------------------------------------------------

#[test]
fn bc_segments_conflicts_with_bc5d_table_path() {
    let mut request = with_bc_segments(base_request(), two_segment_ladder());
    request["corrections"] = json!({"bc5d_table_path": "/some/path/bc5d_308.bin"});
    expect_error(
        &request,
        SolveErrorCodeV1::ConflictingFields,
        "$.projectile.bc_segments",
    );
}

// ---------------------------------------------------------------------------
// (f) Validation errors name the exact path
// ---------------------------------------------------------------------------

#[test]
fn bc_segments_must_be_an_array() {
    let mut request = base_request();
    request["projectile"]["bc_segments"] = json!("not_an_array");
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments",
    );
}

#[test]
fn bc_segments_must_not_be_empty() {
    let request = with_bc_segments(base_request(), json!([]));
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments",
    );
}

#[test]
fn bc_segments_segment_missing_required_field() {
    // Missing velocity_to_mps. The validator reports the exact missing field path.
    let request = with_bc_segments(
        base_request(),
        json!([{"velocity_from_mps": 274.3, "ballistic_coefficient": 0.260}]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::MissingField,
        "$.projectile.bc_segments[0].velocity_to_mps",
    );
}

#[test]
fn bc_segments_segment_unknown_field_rejected() {
    // Unknown field. The validator reports the exact unknown field path.
    let request = with_bc_segments(
        base_request(),
        json!([{
            "velocity_from_mps": 274.3,
            "velocity_to_mps": 999.0,
            "ballistic_coefficient": 0.260,
            "unknown_extra": true
        }]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::UnknownField,
        "$.projectile.bc_segments[0].unknown_extra",
    );
}

#[test]
fn bc_segments_non_positive_bc_rejected() {
    for bad_bc in [0.0_f64, -0.1_f64] {
        let request = with_bc_segments(
            base_request(),
            json!([{
                "velocity_from_mps": 0.0,
                "velocity_to_mps": 999.0,
                "ballistic_coefficient": bad_bc
            }]),
        );
        expect_error(
            &request,
            SolveErrorCodeV1::InvalidValue,
            "$.projectile.bc_segments[0].ballistic_coefficient",
        );
    }
}

#[test]
fn bc_segments_velocity_to_must_exceed_velocity_from() {
    // velocity_to_mps == velocity_from_mps (degenerate band).
    let request = with_bc_segments(
        base_request(),
        json!([{
            "velocity_from_mps": 274.3,
            "velocity_to_mps": 274.3,
            "ballistic_coefficient": 0.260
        }]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments[0].velocity_to_mps",
    );

    // velocity_to_mps < velocity_from_mps.
    let request = with_bc_segments(
        base_request(),
        json!([{
            "velocity_from_mps": 500.0,
            "velocity_to_mps": 274.3,
            "ballistic_coefficient": 0.260
        }]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments[0].velocity_to_mps",
    );
}

#[test]
fn bc_segments_must_be_in_descending_order_of_velocity_from() {
    // Ascending order (lowest velocity first) — wrong.
    let request = with_bc_segments(
        base_request(),
        json!([
            {"velocity_from_mps": 0.0,   "velocity_to_mps": 274.3, "ballistic_coefficient": 0.220},
            {"velocity_from_mps": 274.3, "velocity_to_mps": 999.0, "ballistic_coefficient": 0.260}
        ]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments[1].velocity_from_mps",
    );

    // Duplicate velocity_from_mps — also non-descending.
    let request = with_bc_segments(
        base_request(),
        json!([
            {"velocity_from_mps": 274.3, "velocity_to_mps": 999.0, "ballistic_coefficient": 0.260},
            {"velocity_from_mps": 274.3, "velocity_to_mps": 500.0, "ballistic_coefficient": 0.240}
        ]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments[1].velocity_from_mps",
    );
}

#[test]
fn negative_velocity_from_rejected() {
    let request = with_bc_segments(
        base_request(),
        json!([{
            "velocity_from_mps": -1.0,
            "velocity_to_mps": 999.0,
            "ballistic_coefficient": 0.260
        }]),
    );
    expect_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments[0].velocity_from_mps",
    );
}

// ---------------------------------------------------------------------------
// (h) Service refuses a bad segment even when the decoder was bypassed
// ---------------------------------------------------------------------------

#[test]
fn service_refuses_bad_segment_that_bypassed_decoder() {
    use ballistics_engine::solve_json::{
        AtmosphereV1, DragModelV1, EffectsV1, RifleV1, SamplingV1, SchemaVersionV1, ShotV1,
        SolveRequestV1, SolverV1, WindV1,
    };

    // Build a request in Rust with an invalid segment BC (zero is not positive).
    let typed_request = SolveRequestV1 {
        schema_version: SchemaVersionV1,
        projectile: ProjectileV1 {
            mass_kg: 0.01134,
            diameter_m: 0.00782,
            drag_model: DragModelV1::G7,
            ballistic_coefficient: 0.243,
            length_m: None,
            drag_table: None,
            bc_segments: Some(vec![BcSegmentV1 {
                velocity_from_mps: 274.3,
                velocity_to_mps: 999.0,
                ballistic_coefficient: 0.0, // invalid
            }]),
        },
        rifle: RifleV1 {
            muzzle_velocity_mps: 792.48,
            sight_height_m: Some(0.0508),
            muzzle_height_m: None,
            twist_rate_m_per_turn: None,
            twist_direction: None,
            sight_offset_lateral_m: None,
        },
        shot: ShotV1 {
            max_range_m: 800.0,
            zero_distance_m: Some(100.0),
            muzzle_angle_rad: None,
            aim_azimuth_rad: None,
            shot_azimuth_rad: None,
            shooting_angle_rad: None,
            cant_angle_rad: None,
            target_height_m: None,
            ground_threshold_m: None,
            zero_poi_up_m: None,
            zero_poi_right_m: None,
            drops_reference: None,
        },
        atmosphere: AtmosphereV1::default(),
        wind: WindV1::default(),
        solver: SolverV1::default(),
        effects: EffectsV1::default(),
        sampling: SamplingV1::default(),
        reticle: None,
        corrections: None,
    };

    let error = solve_v1(typed_request).expect_err("service must refuse bc=0.0");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(
        error.error.path(),
        Some("$.projectile.bc_segments[0].ballistic_coefficient")
    );
}

// ---------------------------------------------------------------------------
// (g) Bridge solve command accepts bc_segments (only compiled with `bridge` feature)
// ---------------------------------------------------------------------------

#[cfg(feature = "bridge")]
#[test]
fn bridge_solve_accepts_bc_segments() {
    let envelope = json!({
        "api_version": 1,
        "command": "solve",
        "request": with_bc_segments(base_request(), two_segment_ladder()),
    });
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap();

    assert_eq!(out["ok"], json!(true), "bridge solve must succeed: {out}");

    let echo = &out["result"]["resolved_request"]["projectile"]["bc_segments"];
    assert!(
        echo.is_array(),
        "bridge response must echo bc_segments; got: {echo}"
    );
    assert_eq!(
        echo.as_array().unwrap().len(),
        2,
        "bc_segments echo must have 2 segments"
    );
    assert_eq!(echo[0]["ballistic_coefficient"], json!(0.260));
    assert_eq!(echo[1]["ballistic_coefficient"], json!(0.220));
}
