//! `projectile.bc_segments` on the solve-json v1 wire (MBA-1598): a velocity-banded BC ladder,
//! the wire form of the CLI's `--bc-segment`, so an app on the JSON bridge can fly a published
//! multi-BC schedule from plain request data.
//!
//! This file pins the wire contract:
//!
//!   (a) omitting the field changes nothing, and the key appears nowhere in the response;
//!   (b) the ladder really flies: crossing a band edge mid-flight switches to the slower
//!       band's BC, while the trajectory before the edge is untouched, and a band slow enough
//!       to stop the bullet downrange ends the flight as an ordinary early termination;
//!   (c) the resolved request echoes the ladder verbatim and re-solves identically from it,
//!       through the same JSON decoder a caller uses;
//!   (d) the perturbation kernels behind `explain`, `error-budget` and `tolerance` see the
//!       scalar BC only where it flies: in an interior gap, never under a gapless ladder of
//!       two or more bands. Under a one-band ladder the BC axis moves the band with the scalar;
//!   (e) the ladder cannot be combined with `projectile.drag_table`, and a one-band ladder's BC
//!       must equal the scalar BC;
//!   (f) the decoder itself refuses every bad ladder at the exact entry, and the service
//!       refuses the same ladders when a Rust caller skips the decoder;
//!   (g) the bridge `solve` command accepts and echoes a ladder.
//!
//! The m/s-to-fps conversion is pinned in `src/solve_v1.rs`'s unit tests, against a solve
//! built directly from engine `BCSegmentData` in fps.

use ballistics_engine::solve_json::{
    decode_solve_request_v1, BcSegmentV1, SampleFlagV1, SolveErrorCodeV1, SolveErrorEnvelopeV1,
    SolveRequestV1, SolveSuccessV1, TerminationReasonV1, MAX_BC_SEGMENTS_V1,
    ONE_BAND_BC_RELATIVE_TOLERANCE_V1,
};
use ballistics_engine::solve_v1;
use serde_json::{json, Value};

// ---------------------------------------------------------------------------
// Shared fixture helpers
// ---------------------------------------------------------------------------

/// A .308 175 gr at 792 m/s, zeroed at 100 m, solved to 800 m in still ICAO air. It slows
/// through 600 m/s near 330 m and reaches about 375 m/s at 800 m.
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

fn band(min: f64, max: f64, bc: f64) -> Value {
    json!({"velocity_min_mps": min, "velocity_max_mps": max, "ballistic_coefficient": bc})
}

/// Two touching bands with one edge at 600 m/s, which the fixture crosses mid-flight. Both BCs
/// differ from the scalar 0.243, which never flies under this gapless ladder.
fn two_band_ladder() -> Value {
    json!([band(600.0, 999.0, 0.260), band(0.0, 600.0, 0.220)])
}

fn with_bc_segments(mut request: Value, segments: Value) -> Value {
    request["projectile"]["bc_segments"] = segments;
    request
}

fn decode(request: &Value) -> Result<SolveRequestV1, SolveErrorEnvelopeV1> {
    decode_solve_request_v1(&request.to_string())
}

fn solve(request: &Value) -> SolveSuccessV1 {
    let decoded = decode(request).expect("request should decode without error");
    solve_v1(decoded).expect("solve should succeed")
}

/// The decoder alone must refuse the request, with this code at this path.
fn expect_decode_error(request: &Value, code: SolveErrorCodeV1, path: &str) {
    let error = decode(request).expect_err("the decoder must refuse this request");
    assert_eq!(error.error.code, code, "{}", error.error.message);
    assert_eq!(error.error.path(), Some(path), "{}", error.error.message);
}

fn expect_ladder_refused(segments: Value, path: &str) {
    expect_decode_error(
        &with_bc_segments(base_request(), segments),
        SolveErrorCodeV1::InvalidValue,
        path,
    );
}

/// `count` touching 10 m/s bands, fastest first, from 1000 m/s down.
fn touching_bands(count: usize) -> Value {
    Value::Array(
        (0..count)
            .map(|i| {
                let max = 1000.0 - 10.0 * i as f64;
                band(max - 10.0, max, 0.243)
            })
            .collect(),
    )
}

// ---------------------------------------------------------------------------
// (a) Omitting the field leaves no trace
// ---------------------------------------------------------------------------

#[test]
fn omitting_the_ladder_leaves_no_trace_in_the_response() {
    let response = solve(&base_request());
    let text = serde_json::to_string(&response).unwrap();
    assert!(
        !text.contains("bc_segments"),
        "an omitted bc_segments must not appear anywhere in the response"
    );

    // The resolved request survives its own serialization, and the request rebuilt from it
    // survives the decoder, both without growing a ladder.
    let resolved = &response.resolved_request;
    let resolved_text = serde_json::to_string(resolved).unwrap();
    let reparsed: ballistics_engine::ResolvedSolveRequestV1 =
        serde_json::from_str(&resolved_text).expect("resolved request parses back");
    assert_eq!(&reparsed, resolved);

    let rebuilt = SolveRequestV1::from(resolved);
    let redecoded = decode_solve_request_v1(&serde_json::to_string(&rebuilt).unwrap())
        .expect("a rebuilt request must pass the decoder");
    assert_eq!(redecoded, rebuilt);
    assert!(redecoded.projectile.bc_segments.is_none());
    assert_eq!(solve_v1(redecoded).unwrap().samples, response.samples);
}

// ---------------------------------------------------------------------------
// (b) The ladder flies
// ---------------------------------------------------------------------------

#[test]
fn crossing_a_band_edge_mid_flight_flies_the_slower_band() {
    let ladder = solve(&with_bc_segments(base_request(), two_band_ladder()));
    // The fast band alone: a single band clamps to itself at every speed. A one-band ladder
    // must carry the scalar BC, so the scalar moves with it.
    let mut fast_only = with_bc_segments(base_request(), json!([band(600.0, 999.0, 0.260)]));
    fast_only["projectile"]["ballistic_coefficient"] = json!(0.260);
    let fast_only = solve(&fast_only);
    assert_eq!(ladder.samples.len(), fast_only.samples.len());

    // Above the edge and its blend window (at most +/-7.62 m/s) both fly 0.260 exactly.
    let mut crossed = false;
    for (with_low, without) in ladder.samples.iter().zip(&fast_only.samples) {
        if without.speed_mps > 610.0 {
            assert_eq!(
                with_low, without,
                "before the edge the low band must not fly"
            );
        } else {
            crossed = true;
            assert!(
                with_low.speed_mps < without.speed_mps && with_low.drop_m > without.drop_m,
                "past the edge the lower BC must slow the bullet and drop it more: \
                 {with_low:?} vs {without:?}"
            );
        }
    }
    assert!(crossed, "the fixture must slow through the 600 m/s edge");

    let (ladder_end, fast_end) = (
        ladder.samples.last().unwrap(),
        fast_only.samples.last().unwrap(),
    );
    assert!(
        ladder_end.drop_m - fast_end.drop_m > 0.05,
        "the low band must move the drop at 800 m by centimetres: {} vs {}",
        ladder_end.drop_m,
        fast_end.drop_m
    );
}

/// MBA-1599: a band at the bottom of the BC range sheds the bullet's downrange speed within a
/// few metres of its edge. The bullet then falls almost straight down until its downrange
/// position stops changing in floating point. That is an ordinary early termination at the
/// ground threshold where downrange motion stopped, not an `internal_error`.
#[test]
fn a_band_slow_enough_to_stop_the_bullet_ends_at_the_ground_threshold() {
    let mut request = with_bc_segments(
        base_request(),
        json!([band(300.0, 1200.0, 0.243), band(0.0, 300.0, 0.0001)]),
    );
    request["shot"]["max_range_m"] = json!(2000.0);
    let response = solve(&request);

    assert_eq!(
        response.summary.termination,
        TerminationReasonV1::GroundThreshold
    );
    let range = response.summary.actual_range_m;
    assert!(
        (1100.0..1150.0).contains(&range),
        "the bullet stops moving downrange just past the 300 m/s edge, got {range} m"
    );
    assert_eq!(
        response
            .samples
            .iter()
            .filter(|sample| sample.flags.contains(&SampleFlagV1::Terminal))
            .count(),
        1
    );
    let terminal = response.samples.last().expect("terminal sample");
    assert!(terminal.flags.contains(&SampleFlagV1::GroundThreshold));
    assert_eq!(terminal.distance_m, range);
    assert!(
        terminal.speed_mps < 10.0,
        "the bullet ends in a near-vertical fall, got {} m/s",
        terminal.speed_mps
    );
    assert!(response
        .samples
        .windows(2)
        .all(|pair| pair[1].distance_m > pair[0].distance_m));
}

// ---------------------------------------------------------------------------
// (c) Echo and round trip
// ---------------------------------------------------------------------------

#[test]
fn resolved_request_echoes_the_ladder_verbatim() {
    let request = with_bc_segments(base_request(), two_band_ladder());
    let response = solve(&request);

    let echoed = serde_json::to_value(&response.resolved_request).unwrap();
    assert_eq!(
        echoed["projectile"]["bc_segments"],
        request["projectile"]["bc_segments"]
    );
    assert_eq!(
        response.resolved_request.projectile.bc_segments,
        Some(vec![
            BcSegmentV1 {
                velocity_min_mps: 600.0,
                velocity_max_mps: 999.0,
                ballistic_coefficient: 0.260,
            },
            BcSegmentV1 {
                velocity_min_mps: 0.0,
                velocity_max_mps: 600.0,
                ballistic_coefficient: 0.220,
            },
        ])
    );
}

#[test]
fn a_request_rebuilt_from_the_echo_re_solves_identically_through_the_decoder() {
    let first = solve(&with_bc_segments(base_request(), two_band_ladder()));

    let rebuilt = SolveRequestV1::from(&first.resolved_request);
    let redecoded = decode_solve_request_v1(&serde_json::to_string(&rebuilt).unwrap())
        .expect("a rebuilt request must pass the decoder");
    assert_eq!(redecoded, rebuilt);
    assert!(redecoded.projectile.bc_segments.is_some());

    let second = solve_v1(redecoded).expect("the rebuilt request must solve");
    assert_eq!(first.samples, second.samples);
    assert_eq!(first.summary, second.summary);
}

// ---------------------------------------------------------------------------
// (d) The perturbation kernels
// ---------------------------------------------------------------------------

/// The kernels perturb `ballistic_coefficient` on the resolved request and re-solve it. Under a
/// gapless ladder the scalar never flies, so its slope is exactly zero: a known limitation, the
/// same one a `kind: "projectile"` drag table has. With an interior gap the bullet crosses, the
/// scalar flies there and the slope is real, which also proves the re-solves kept the ladder.
#[test]
fn the_perturbation_kernel_sees_the_scalar_bc_only_in_an_interior_gap() {
    use ballistics_engine::perturbation::{central_difference, InputAxis};

    let bc_slope = |ladder: Value| {
        let response = solve(&with_bc_segments(base_request(), ladder));
        central_difference(
            &response.resolved_request,
            InputAxis::BallisticCoefficient,
            &[600.0],
            None,
        )
        .expect("the BC axis is continuous")[0]
            .d_drop_d_x
    };

    assert_eq!(bc_slope(two_band_ladder()), 0.0);

    // A 500-650 m/s gap, which the fixture crosses between about 230 m and 520 m.
    let gapped = bc_slope(json!([band(650.0, 999.0, 0.260), band(0.0, 500.0, 0.220)]));
    assert!(
        gapped < 0.0,
        "a higher scalar BC flown through the gap must reduce the drop: dDrop/dBC = {gapped}"
    );
}

/// A one-band ladder must carry the scalar BC (section (e)), and that one BC flies at every
/// speed, so the BC axis moves the band with the scalar. Each perturbed solve then flies the
/// same BC a request without the ladder would, and the slope is the plain scalar's exactly.
#[test]
fn the_perturbation_kernel_moves_a_one_band_ladder_with_the_scalar_bc() {
    use ballistics_engine::perturbation::{central_difference, InputAxis};

    let bc_slope = |request: &Value| {
        central_difference(
            &solve(request).resolved_request,
            InputAxis::BallisticCoefficient,
            &[600.0],
            None,
        )
        .expect("the BC axis is continuous")[0]
            .d_drop_d_x
    };

    let plain = bc_slope(&base_request());
    assert!(
        plain < 0.0,
        "a higher BC must reduce the drop: dDrop/dBC = {plain}"
    );
    let one_band = bc_slope(&with_bc_segments(
        base_request(),
        json!([band(0.0, 999.0, 0.243)]),
    ));
    assert_eq!(one_band, plain);
}

// ---------------------------------------------------------------------------
// (e) Conflicts: projectile.drag_table, and a one-band ladder against the scalar BC
// ---------------------------------------------------------------------------

#[test]
fn a_ladder_cannot_be_combined_with_a_drag_table() {
    let mut request = with_bc_segments(base_request(), two_band_ladder());
    request["projectile"]["drag_table"] = json!({
        "kind": "reference",
        "points": [
            {"mach": 0.0, "cd": 0.2},
            {"mach": 2.0, "cd": 0.3}
        ]
    });
    let error = solve_v1(decode(&request).expect("each field is valid on its own"))
        .expect_err("the pair must be refused");
    assert_eq!(error.error.code, SolveErrorCodeV1::ConflictingFields);
    assert_eq!(error.error.path(), Some("$.projectile.bc_segments"));
    assert_eq!(
        error.error.message,
        "projectile.bc_segments cannot be combined with projectile.drag_table: a drag table \
         replaces the drag_model curve that the band BCs are measured against"
    );
}

/// A one-band ladder flies its band's BC at every speed, so the scalar never flies beside it,
/// and two different BCs for one projectile can only be a caller mistake. The decoder refuses
/// the pair at the band's BC, and so does the service for a Rust caller that skips the
/// decoder. Equal values fly exactly as the scalar alone, wherever the band's edges sit.
#[test]
fn a_one_band_ladder_must_carry_the_scalar_bc() {
    let path = "$.projectile.bc_segments[0].ballistic_coefficient";
    let message = "a one-band ladder flies its band's BC at every speed, so \
                   bc_segments[0].ballistic_coefficient (0.26) must match \
                   projectile.ballistic_coefficient (0.243): send the same value in both, or \
                   omit bc_segments";

    let mismatched = with_bc_segments(base_request(), json!([band(0.0, 999.0, 0.260)]));
    let error = decode(&mismatched).expect_err("the decoder must refuse two different BCs");
    assert_eq!(error.error.code, SolveErrorCodeV1::ConflictingFields);
    assert_eq!(error.error.path(), Some(path));
    assert_eq!(error.error.message, message);

    let mut skipped = decode(&base_request()).unwrap();
    skipped.projectile.bc_segments = Some(vec![BcSegmentV1 {
        velocity_min_mps: 0.0,
        velocity_max_mps: 999.0,
        ballistic_coefficient: 0.260,
    }]);
    let error = solve_v1(skipped).expect_err("the service must refuse it too");
    assert_eq!(error.error.code, SolveErrorCodeV1::ConflictingFields);
    assert_eq!(error.error.path(), Some(path));
    assert_eq!(error.error.message, message);

    // Below every flight speed, across the whole flight, and above the muzzle velocity.
    let plain = solve(&base_request());
    for (min, max) in [(100.0, 200.0), (0.0, 999.0), (900.0, 1500.0)] {
        let one_band = solve(&with_bc_segments(
            base_request(),
            json!([band(min, max, 0.243)]),
        ));
        assert_eq!(one_band.samples, plain.samples, "band {min} to {max} m/s");
    }

    // The rule is for one band only: a longer ladder's band BCs need not match the scalar.
    solve(&with_bc_segments(base_request(), two_band_ladder()));
}

/// Rounding nobody chose must not refuse a solve: a BC held as an `f32` and widened back, or a
/// value inside one part in a million of the scalar, is the same BC, and the band flies. Just
/// past the tolerance, in either direction, it is a different BC again.
#[test]
fn a_one_band_ladder_forgives_rounding_noise_but_not_a_different_bc() {
    let scalar = 0.243_f64;
    let widened_f32 = f64::from(0.243_f32);
    assert_ne!(widened_f32, scalar, "the fixture needs real f32 noise");
    let inside = scalar * (1.0 + 0.5 * ONE_BAND_BC_RELATIVE_TOLERANCE_V1);
    let below = scalar * (1.0 - 0.5 * ONE_BAND_BC_RELATIVE_TOLERANCE_V1);
    for band_bc in [widened_f32, inside, below] {
        let response = solve(&with_bc_segments(
            base_request(),
            json!([band(0.0, 999.0, band_bc)]),
        ));
        assert_eq!(
            response.resolved_request.projectile.bc_segments.unwrap()[0].ballistic_coefficient,
            band_bc,
            "the band is echoed, and flown, as sent"
        );
    }

    for band_bc in [
        scalar * (1.0 + 2.0 * ONE_BAND_BC_RELATIVE_TOLERANCE_V1),
        scalar * (1.0 - 2.0 * ONE_BAND_BC_RELATIVE_TOLERANCE_V1),
    ] {
        expect_decode_error(
            &with_bc_segments(base_request(), json!([band(0.0, 999.0, band_bc)])),
            SolveErrorCodeV1::ConflictingFields,
            "$.projectile.bc_segments[0].ballistic_coefficient",
        );
    }
}

// ---------------------------------------------------------------------------
// (f) Validation: shape
// ---------------------------------------------------------------------------

#[test]
fn null_ladder_is_not_omission() {
    let mut request = base_request();
    request["projectile"]["bc_segments"] = Value::Null;
    expect_decode_error(
        &request,
        SolveErrorCodeV1::InvalidValue,
        "$.projectile.bc_segments",
    );
}

#[test]
fn the_ladder_must_be_an_array() {
    expect_ladder_refused(json!("not_an_array"), "$.projectile.bc_segments");
    expect_ladder_refused(band(0.0, 999.0, 0.26), "$.projectile.bc_segments");
}

#[test]
fn an_empty_ladder_is_refused() {
    expect_ladder_refused(json!([]), "$.projectile.bc_segments");
}

#[test]
fn a_band_missing_a_member_is_refused_at_that_member() {
    expect_decode_error(
        &with_bc_segments(
            base_request(),
            json!([{"velocity_min_mps": 274.3, "ballistic_coefficient": 0.260}]),
        ),
        SolveErrorCodeV1::MissingField,
        "$.projectile.bc_segments[0].velocity_max_mps",
    );
}

#[test]
fn unknown_and_old_member_names_are_refused_at_their_path() {
    let mut extra = band(274.3, 999.0, 0.260);
    extra["unknown_extra"] = json!(true);
    expect_decode_error(
        &with_bc_segments(base_request(), json!([extra])),
        SolveErrorCodeV1::UnknownField,
        "$.projectile.bc_segments[0].unknown_extra",
    );

    expect_decode_error(
        &with_bc_segments(
            base_request(),
            json!([{"velocity_from_mps": 274.3, "velocity_max_mps": 999.0,
                    "ballistic_coefficient": 0.260}]),
        ),
        SolveErrorCodeV1::UnknownField,
        "$.projectile.bc_segments[0].velocity_from_mps",
    );
}

#[test]
fn a_band_member_must_be_a_number() {
    expect_ladder_refused(
        json!([{"velocity_min_mps": 0.0, "velocity_max_mps": "999",
                "ballistic_coefficient": 0.26}]),
        "$.projectile.bc_segments[0].velocity_max_mps",
    );
}

// ---------------------------------------------------------------------------
// (f) Validation: values, refused by the decoder itself
// ---------------------------------------------------------------------------

#[test]
fn band_bc_is_held_to_the_scalar_bc_bounds() {
    for bad_bc in [0.0, -0.1, 0.00009, 100.5, 1.0e6] {
        expect_ladder_refused(
            json!([band(0.0, 999.0, bad_bc)]),
            "$.projectile.bc_segments[0].ballistic_coefficient",
        );
    }
    // The bounds themselves are legal, as they are for the scalar.
    // Two bands, because a one-band ladder must also match the scalar 0.243.
    for edge_bc in [1.0e-4, 100.0] {
        decode(&with_bc_segments(
            base_request(),
            json!([band(500.0, 999.0, edge_bc), band(0.0, 500.0, 0.243)]),
        ))
        .expect("a band BC at the scalar's bound is legal");
    }

    // The message reads exactly as the scalar's does.
    let error = decode(&with_bc_segments(
        base_request(),
        json!([band(0.0, 999.0, 100.5)]),
    ))
    .unwrap_err();
    let mut scalar = base_request();
    scalar["projectile"]["ballistic_coefficient"] = json!(100.5);
    let scalar_error = decode(&scalar).unwrap_err();
    assert_eq!(
        error.error.message.replace("bc_segments[0].", ""),
        scalar_error.error.message
    );
}

#[test]
fn velocity_max_must_exceed_velocity_min() {
    expect_ladder_refused(
        json!([band(274.3, 274.3, 0.260)]),
        "$.projectile.bc_segments[0].velocity_max_mps",
    );
    expect_ladder_refused(
        json!([band(500.0, 274.3, 0.260)]),
        "$.projectile.bc_segments[0].velocity_max_mps",
    );
}

#[test]
fn velocity_min_must_not_be_negative() {
    expect_ladder_refused(
        json!([band(-1.0, 999.0, 0.260)]),
        "$.projectile.bc_segments[0].velocity_min_mps",
    );
}

#[test]
fn bands_run_from_the_highest_speed_to_the_lowest() {
    let ascending = with_bc_segments(
        base_request(),
        json!([band(0.0, 600.0, 0.220), band(600.0, 999.0, 0.260)]),
    );
    let error = decode(&ascending).expect_err("an ascending ladder must be refused");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(
        error.error.path(),
        Some("$.projectile.bc_segments[1].velocity_min_mps")
    );
    for hint in ["bc_segments[0]", "highest speed to the lowest", "reverse"] {
        assert!(
            error.error.message.contains(hint),
            "the message must say {hint:?}: {}",
            error.error.message
        );
    }

    // A repeated start is not descending either.
    expect_ladder_refused(
        json!([band(274.3, 999.0, 0.260), band(274.3, 500.0, 0.240)]),
        "$.projectile.bc_segments[1].velocity_min_mps",
    );
}

#[test]
fn bands_must_not_overlap() {
    let overlapping = with_bc_segments(
        base_request(),
        json!([band(600.0, 999.0, 0.260), band(400.0, 650.0, 0.220)]),
    );
    let error = decode(&overlapping).expect_err("an overlap must be refused");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(
        error.error.path(),
        Some("$.projectile.bc_segments[1].velocity_max_mps")
    );
    assert!(
        error.error.message.contains("bc_segments[0]"),
        "{}",
        error.error.message
    );

    // A band nested inside the one before it starts above it, so the order rule refuses it;
    // a third band reaching past the second's start is an overlap at that band.
    expect_ladder_refused(
        json!([band(300.0, 999.0, 0.260), band(400.0, 500.0, 0.220)]),
        "$.projectile.bc_segments[1].velocity_min_mps",
    );
    expect_ladder_refused(
        json!([
            band(700.0, 999.0, 0.260),
            band(500.0, 700.0, 0.240),
            band(100.0, 500.1, 0.220)
        ]),
        "$.projectile.bc_segments[2].velocity_max_mps",
    );
}

#[test]
fn touching_edges_and_interior_gaps_are_both_legal() {
    solve(&with_bc_segments(base_request(), two_band_ladder()));
    solve(&with_bc_segments(
        base_request(),
        json!([band(650.0, 999.0, 0.260), band(0.0, 500.0, 0.220)]),
    ));
    // A gap below the slowest band and a ladder that never reaches the muzzle are legal too.
    solve(&with_bc_segments(
        base_request(),
        json!([band(500.0, 700.0, 0.250), band(300.0, 450.0, 0.230)]),
    ));
}

#[test]
fn band_count_is_capped() {
    expect_ladder_refused(
        touching_bands(MAX_BC_SEGMENTS_V1 + 1),
        "$.projectile.bc_segments",
    );
    solve(&with_bc_segments(
        base_request(),
        touching_bands(MAX_BC_SEGMENTS_V1),
    ));
    solve(&with_bc_segments(base_request(), touching_bands(1)));
}

/// Rust callers can build a request without the decoder; the service still refuses a bad
/// ladder at the same path.
#[test]
fn the_service_refuses_a_bad_ladder_that_skipped_the_decoder() {
    let good = decode(&with_bc_segments(base_request(), two_band_ladder())).unwrap();
    let seg = |min: f64, max: f64, bc: f64| BcSegmentV1 {
        velocity_min_mps: min,
        velocity_max_mps: max,
        ballistic_coefficient: bc,
    };
    let cases = [
        (
            vec![seg(274.3, 999.0, 0.0)],
            "$.projectile.bc_segments[0].ballistic_coefficient",
        ),
        (
            vec![seg(274.3, 999.0, f64::NAN)],
            "$.projectile.bc_segments[0].ballistic_coefficient",
        ),
        (
            vec![seg(0.0, f64::INFINITY, 0.26)],
            "$.projectile.bc_segments[0].velocity_max_mps",
        ),
        (
            vec![seg(0.0, 600.0, 0.22), seg(600.0, 999.0, 0.26)],
            "$.projectile.bc_segments[1].velocity_min_mps",
        ),
        (
            vec![seg(600.0, 999.0, 0.26), seg(0.0, 601.0, 0.22)],
            "$.projectile.bc_segments[1].velocity_max_mps",
        ),
        (Vec::new(), "$.projectile.bc_segments"),
        (
            vec![seg(600.0, 999.0, 0.26); MAX_BC_SEGMENTS_V1 + 1],
            "$.projectile.bc_segments",
        ),
    ];
    for (ladder, path) in cases {
        let mut request = good.clone();
        request.projectile.bc_segments = Some(ladder);
        let error = solve_v1(request).expect_err("the service must refuse a bad ladder");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.error.path(), Some(path), "{}", error.error.message);
    }
}

// ---------------------------------------------------------------------------
// (g) Bridge
// ---------------------------------------------------------------------------

#[cfg(feature = "bridge")]
#[test]
fn the_bridge_solve_command_accepts_a_ladder() {
    let envelope = json!({
        "api_version": 1,
        "command": "solve",
        "request": with_bc_segments(base_request(), two_band_ladder()),
    });
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap();

    assert_eq!(out["ok"], json!(true), "bridge solve must succeed: {out}");
    assert_eq!(
        out["result"]["resolved_request"]["projectile"]["bc_segments"],
        two_band_ladder()
    );
}

/// A bad ladder is refused by the decoder, so the bridge reports a rejected request rather
/// than a failed solve.
#[cfg(feature = "bridge")]
#[test]
fn the_bridge_rejects_a_bad_ladder_as_a_request_error() {
    let envelope = json!({
        "api_version": 1,
        "command": "solve",
        "request": with_bc_segments(
            base_request(),
            json!([band(600.0, 999.0, 0.260), band(0.0, 650.0, 0.220)])
        ),
    });
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap();
    assert_eq!(out["ok"], json!(false), "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["path"],
        json!("$.projectile.bc_segments[1].velocity_max_mps"),
        "{out}"
    );
    assert!(
        out["error"]["message"]
            .as_str()
            .is_some_and(|message| message.contains("rejected")),
        "{out}"
    );
}

/// A one-band ladder that disagrees with the scalar BC is refused by the decoder, so the bridge
/// reports a rejected request, as it does for every other bad ladder.
#[cfg(feature = "bridge")]
#[test]
fn the_bridge_rejects_a_one_band_ladder_that_disagrees_with_the_scalar() {
    let envelope = json!({
        "api_version": 1,
        "command": "solve",
        "request": with_bc_segments(base_request(), json!([band(0.0, 274.32, 0.150)])),
    });
    let out: Value = serde_json::from_str(&ballistics_engine::bridge::bridge_call(
        &envelope.to_string(),
    ))
    .unwrap();
    assert_eq!(out["ok"], json!(false), "{out}");
    assert_eq!(
        out["error"]["details"]["error"]["code"],
        json!("conflicting_fields"),
        "{out}"
    );
    assert_eq!(
        out["error"]["details"]["error"]["path"],
        json!("$.projectile.bc_segments[0].ballistic_coefficient"),
        "{out}"
    );
    assert!(
        out["error"]["message"]
            .as_str()
            .is_some_and(|message| message.contains("rejected")),
        "{out}"
    );
}
