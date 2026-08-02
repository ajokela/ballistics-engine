use approx::assert_relative_eq;
use ballistics_engine::solve_json::{
    decode_solve_request_v1, SampleFlagV1, SolveRequestV1, SolveSummaryV1, SolverMethodV1,
    TerminationReasonV1,
};
use ballistics_engine::wind::WindSegment;
use ballistics_engine::{
    calculate_zero_angle_with_conditions, solve_v1, AtmosphericConditions, BallisticInputs,
    DragModel, TrajectorySolver, WindConditions,
};
use serde_json::{json, Value};

fn request_value(method: &str) -> Value {
    json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 0.01134,
            "diameter_m": 0.00782,
            "length_m": 0.031,
            "drag_model": "G7",
            "ballistic_coefficient": 0.243
        },
        "rifle": {
            "muzzle_velocity_mps": 823.0,
            "sight_height_m": 0.05,
            "muzzle_height_m": 1.2,
            "twist_rate_m_per_turn": 0.2032,
            "twist_direction": "right"
        },
        "shot": {
            "max_range_m": 200.0,
            "muzzle_angle_rad": 0.001,
            "aim_azimuth_rad": 0.0,
            "shot_azimuth_rad": 0.0,
            "shooting_angle_rad": 0.0,
            "cant_angle_rad": 0.0,
            "target_height_m": 1.25,
            "ground_threshold_m": -100.0
        },
        "atmosphere": {
            "altitude_m": 0.0,
            "temperature_k": 293.15,
            "pressure_pa": 100000.0,
            "relative_humidity": 0.4
        },
        "wind": {
            "speed_mps": 4.0,
            "direction_from_rad": std::f64::consts::FRAC_PI_2,
            "vertical_speed_mps": 0.25
        },
        "solver": {
            "method": method,
            "time_step_s": 0.001
        },
        "effects": {
            "magnus": false,
            "coriolis": false,
            "enhanced_spin_drift": false
        },
        "sampling": {
            "interval_m": 25.0
        }
    })
}

fn decode(value: &Value) -> SolveRequestV1 {
    decode_solve_request_v1(&serde_json::to_string(value).expect("serialize request"))
        .expect("valid request")
}

fn direct_inputs(method: SolverMethodV1) -> BallisticInputs {
    let (use_rk4, use_adaptive_rk45) = match method {
        SolverMethodV1::Euler => (false, false),
        SolverMethodV1::Rk4 => (true, false),
        SolverMethodV1::Rk45 => (true, true),
    };
    BallisticInputs {
        bc_value: 0.243,
        bc_type: DragModel::G7,
        bullet_mass: 0.01134,
        muzzle_velocity: 823.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.031,
        muzzle_angle: 0.001,
        target_distance: 200.0,
        azimuth_angle: 0.0,
        shot_azimuth: 0.0,
        shooting_angle: 0.0,
        cant_angle: 0.0,
        sight_height: 0.05,
        muzzle_height: 1.2,
        target_height: 1.25,
        ground_threshold: -100.0,
        altitude: 0.0,
        temperature: 20.0,
        pressure: 1000.0,
        humidity: 0.4,
        wind_speed: 4.0,
        wind_angle: std::f64::consts::FRAC_PI_2,
        twist_rate: 8.0,
        is_twist_right: true,
        use_rk4,
        use_adaptive_rk45,
        ..BallisticInputs::default()
    }
}

fn direct_solver(method: SolverMethodV1) -> TrajectorySolver {
    let mut solver = TrajectorySolver::new(
        direct_inputs(method),
        WindConditions {
            speed: 4.0,
            direction: std::f64::consts::FRAC_PI_2,
            vertical_speed: 0.25,
        },
        AtmosphericConditions {
            temperature: 20.0,
            pressure: 1000.0,
            humidity: 40.0,
            altitude: 0.0,
        },
    );
    solver.set_max_range(200.0);
    solver.set_time_step(0.001);
    solver
}

#[test]
fn every_solver_matches_the_direct_crosswind_engine_path() {
    for (wire_method, method) in [
        ("euler", SolverMethodV1::Euler),
        ("rk4", SolverMethodV1::Rk4),
        ("rk45", SolverMethodV1::Rk45),
    ] {
        let success = solve_v1(decode(&request_value(wire_method))).expect("service solve");
        let direct = direct_solver(method).solve().expect("direct solve");
        let observations = direct
            .sample_observations(25.0, 10_000)
            .expect("direct observations");

        assert_eq!(success.samples.len(), observations.len(), "{wire_method}");
        assert_eq!(success.summary.termination, TerminationReasonV1::MaxRange);
        for (actual, expected) in success.samples.iter().zip(observations) {
            assert_relative_eq!(actual.distance_m, expected.distance_m, epsilon = 1e-12);
            assert_relative_eq!(actual.time_s, expected.time_s, epsilon = 1e-12);
            assert_relative_eq!(actual.speed_mps, expected.speed_mps, epsilon = 1e-10);
            assert_relative_eq!(actual.energy_j, expected.energy_j, epsilon = 1e-9);
            assert_relative_eq!(actual.drop_m, expected.drop_m, epsilon = 1e-12);
            assert_relative_eq!(actual.windage_m, expected.windage_m, epsilon = 1e-12);
            assert_relative_eq!(actual.mach, expected.mach, epsilon = 1e-12);
        }
    }
}

#[test]
fn zero_distance_uses_the_same_scalar_wind_and_hits_absolute_target_height() {
    let mut value = request_value("rk45");
    value["shot"]
        .as_object_mut()
        .expect("shot object")
        .remove("muzzle_angle_rad");
    value["shot"]["zero_distance_m"] = json!(100.0);
    value["sampling"]["interval_m"] = json!(10.0);

    let first = solve_v1(decode(&value)).expect("zeroed solve");
    let second = solve_v1(decode(&value)).expect("deterministic repeat");
    assert_eq!(first, second);
    assert_ne!(first.resolved_request.shot.muzzle_angle_rad, 0.0);

    let zero_sample = first
        .samples
        .iter()
        .find(|sample| sample.distance_m == 100.0)
        .expect("zero-distance sample");
    let projectile_height_m = first.resolved_request.rifle.muzzle_height_m
        + first.resolved_request.rifle.sight_height_m
        - zero_sample.drop_m;
    assert_relative_eq!(projectile_height_m, 1.25, epsilon = 1e-4);
}

/// 0.33.0 decision-support Task 2: an explicit `muzzle_angle_rad` must win over
/// `zero_distance_m` end-to-end through `solve_v1`, not merely in `prepare_request`'s shape
/// validation -- otherwise a request rebuilt from a resolved request
/// (`request_roundtrip.rs`) would have its carried effective angle silently overwritten by
/// a fresh zero search the moment it actually solves.
#[test]
fn explicit_muzzle_angle_is_not_overridden_by_a_present_zero_distance() {
    // Nonzero, so equivalent_horizontal_range_m (gated on zero_distance_m being present,
    // independent of whether it re-derived the elevation) is actually exercised below --
    // pinned at 0.0 this test would pass even if zero_distance_m secretly still re-zeroed
    // the elevation, because a flat shot's equivalent_horizontal_range_m is always absent.
    let shooting_angle_rad = 0.3_f64;
    let zero_distance_m = 100.0_f64;
    // World-vertical target height consistent with an inclined zero at this angle (mirrors
    // solve_v1_service.rs's inclined_zero_targets_world_vertical_height): a flat-fire
    // target_height_m (1.25, request_value's default) does not bracket a convergent zero
    // once shooting_angle_rad is nonzero.
    let line_of_sight_y_m = 1.2 + 0.05; // rifle.muzzle_height_m + rifle.sight_height_m
    let target_height_m =
        zero_distance_m * shooting_angle_rad.sin() + line_of_sight_y_m * shooting_angle_rad.cos();

    // request_value's default max_range_m (200.0) leaves only 100 m past the zero -- not
    // enough room for a slightly over-elevated inclined shot to arc back below the line of
    // sight, which equivalent_horizontal_range_m's "positive correction" requires. Widen it.
    let max_range_m = 500.0_f64;

    let mut zero_only_value = request_value("rk45");
    zero_only_value["shot"]["shooting_angle_rad"] = json!(shooting_angle_rad);
    zero_only_value["shot"]["target_height_m"] = json!(target_height_m);
    zero_only_value["shot"]["max_range_m"] = json!(max_range_m);
    zero_only_value["shot"]
        .as_object_mut()
        .expect("shot object")
        .remove("muzzle_angle_rad");
    zero_only_value["shot"]["zero_distance_m"] = json!(zero_distance_m);
    let natural = solve_v1(decode(&zero_only_value)).expect("zero-only solve");
    let natural_angle = natural.resolved_request.shot.muzzle_angle_rad;

    // An angle measurably away from the naturally-solved one (the zero search converges to
    // ~1e-7 rad tolerance, so this is far outside coincidental-match range), but still small
    // enough that the shot drops back below the line of sight before max_range_m -- required
    // for equivalent_horizontal_range_m's "positive correction" to be well-defined below.
    let distinct_angle = natural_angle + 0.002;

    let mut both_value = request_value("rk45");
    both_value["shot"]["shooting_angle_rad"] = json!(shooting_angle_rad);
    both_value["shot"]["target_height_m"] = json!(target_height_m);
    both_value["shot"]["max_range_m"] = json!(max_range_m);
    both_value["shot"]["zero_distance_m"] = json!(zero_distance_m);
    both_value["shot"]["muzzle_angle_rad"] = json!(distinct_angle);
    let both = solve_v1(decode(&both_value)).expect("both fields together must solve");

    let mut angle_only_value = request_value("rk45");
    angle_only_value["shot"]["shooting_angle_rad"] = json!(shooting_angle_rad);
    angle_only_value["shot"]["target_height_m"] = json!(target_height_m);
    angle_only_value["shot"]["max_range_m"] = json!(max_range_m);
    angle_only_value["shot"]["muzzle_angle_rad"] = json!(distinct_angle);
    let angle_only = solve_v1(decode(&angle_only_value)).expect("angle-only solve");

    assert_eq!(
        both.resolved_request.shot.muzzle_angle_rad, distinct_angle,
        "an explicit muzzle_angle_rad must not be silently replaced by a zero search"
    );
    assert_eq!(
        both.resolved_request.shot.zero_distance_m,
        Some(zero_distance_m)
    );

    // zero_distance_m is not fully inert just because an explicit angle is also present
    // (0.33.0 decision-support Task 2, review finding I2): it still gates
    // equivalent_horizontal_range_m. Confirm that's exactly the one summary field that
    // legitimately differs between "zero_distance_m present" and "absent", then require
    // everything else -- including the actual solved trajectory, where zero_distance_m
    // truly has no remaining effect once the elevation is explicit -- to match exactly.
    assert!(
        both.summary.equivalent_horizontal_range_m.is_some(),
        "zero_distance_m present alongside an inclined shot should still gate \
         equivalent_horizontal_range_m even though the elevation was supplied directly"
    );
    assert_eq!(angle_only.summary.equivalent_horizontal_range_m, None);
    assert_eq!(
        SolveSummaryV1 {
            equivalent_horizontal_range_m: None,
            ..both.summary.clone()
        },
        angle_only.summary,
        "summary must be identical apart from equivalent_horizontal_range_m"
    );
    assert_eq!(
        both.samples, angle_only.samples,
        "zero_distance_m must have no effect on the solved trajectory once an explicit \
         angle is present"
    );
}

#[test]
fn calm_zeroed_service_matches_the_direct_configured_solver() {
    let mut value = request_value("rk45");
    value["shot"]
        .as_object_mut()
        .expect("shot object")
        .remove("muzzle_angle_rad");
    value["shot"]["zero_distance_m"] = json!(100.0);
    value["wind"] = json!({});
    let success = solve_v1(decode(&value)).expect("calm zeroed service solve");

    let mut direct_inputs = direct_inputs(SolverMethodV1::Rk45);
    direct_inputs.wind_speed = 0.0;
    direct_inputs.wind_angle = 0.0;
    let atmosphere = AtmosphericConditions {
        temperature: 20.0,
        pressure: 1000.0,
        humidity: 40.0,
        altitude: 0.0,
    };
    let direct_angle = calculate_zero_angle_with_conditions(
        direct_inputs.clone(),
        100.0,
        1.25,
        WindConditions::default(),
        atmosphere.clone(),
    )
    .expect("direct zero");
    assert_relative_eq!(
        success.resolved_request.shot.muzzle_angle_rad,
        direct_angle,
        epsilon = 1e-14
    );

    direct_inputs.muzzle_angle = direct_angle;
    let mut direct = TrajectorySolver::new(direct_inputs, WindConditions::default(), atmosphere);
    direct.set_max_range(200.0);
    direct.set_time_step(0.001);
    let observations = direct
        .solve()
        .expect("direct zeroed solve")
        .sample_observations(25.0, 10_000)
        .expect("direct zeroed observations");
    assert_eq!(success.samples.len(), observations.len());
    for (actual, expected) in success.samples.iter().zip(observations) {
        assert_relative_eq!(actual.distance_m, expected.distance_m, epsilon = 1e-12);
        assert_relative_eq!(actual.time_s, expected.time_s, epsilon = 1e-12);
        assert_relative_eq!(actual.speed_mps, expected.speed_mps, epsilon = 1e-10);
        assert_relative_eq!(actual.drop_m, expected.drop_m, epsilon = 1e-12);
        assert_relative_eq!(actual.windage_m, expected.windage_m, epsilon = 1e-12);
    }
}

#[test]
fn explicit_standard_station_values_remain_authoritative_during_zeroing() {
    let mut explicit_value = request_value("rk45");
    explicit_value["shot"]
        .as_object_mut()
        .expect("shot object")
        .remove("muzzle_angle_rad");
    explicit_value["shot"]["zero_distance_m"] = json!(100.0);
    explicit_value["atmosphere"]["altitude_m"] = json!(2_000.0);
    explicit_value["atmosphere"]["temperature_k"] = json!(288.15);
    explicit_value["atmosphere"]["pressure_pa"] = json!(101_325.0);

    let mut omitted_value = explicit_value.clone();
    omitted_value["atmosphere"]
        .as_object_mut()
        .expect("atmosphere object")
        .remove("temperature_k");
    omitted_value["atmosphere"]
        .as_object_mut()
        .expect("atmosphere object")
        .remove("pressure_pa");

    let explicit = solve_v1(decode(&explicit_value)).expect("explicit-atmosphere zero");
    let omitted = solve_v1(decode(&omitted_value)).expect("ICAO-atmosphere zero");
    assert_eq!(explicit.resolved_request.atmosphere.temperature_k, 288.15);
    assert_eq!(explicit.resolved_request.atmosphere.pressure_pa, 101_325.0);
    assert_ne!(
        omitted.resolved_request.atmosphere.temperature_k,
        explicit.resolved_request.atmosphere.temperature_k
    );
    assert_ne!(
        omitted.resolved_request.atmosphere.pressure_pa,
        explicit.resolved_request.atmosphere.pressure_pa
    );
    assert!(omitted.assumptions.iter().any(|notice| {
        notice.code == "icao_standard_temperature"
            && notice.path.as_deref() == Some("$.atmosphere.temperature_k")
    }));
    assert!(omitted.assumptions.iter().any(|notice| {
        notice.code == "icao_standard_pressure"
            && notice.path.as_deref() == Some("$.atmosphere.pressure_pa")
    }));
    assert_ne!(
        omitted.resolved_request.shot.muzzle_angle_rad,
        explicit.resolved_request.shot.muzzle_angle_rad,
        "zero trials must use the independently resolved station atmosphere"
    );
}

#[test]
fn inclined_zero_targets_world_vertical_height() {
    let mut value = request_value("rk45");
    value["shot"]
        .as_object_mut()
        .expect("shot object")
        .remove("muzzle_angle_rad");
    let shooting_angle = 10.0_f64.to_radians();
    let zero_distance_m = 100.0;
    let line_of_sight_y_m = 1.2 + 0.05;
    let target_world_height_m =
        zero_distance_m * shooting_angle.sin() + line_of_sight_y_m * shooting_angle.cos();
    value["shot"]["zero_distance_m"] = json!(zero_distance_m);
    value["shot"]["shooting_angle_rad"] = json!(shooting_angle);
    value["shot"]["target_height_m"] = json!(target_world_height_m);

    let success = solve_v1(decode(&value)).expect("inclined zero");
    let zero_sample = success
        .samples
        .iter()
        .find(|sample| sample.distance_m == zero_distance_m)
        .expect("zero-distance sample");
    let shot_y_m = line_of_sight_y_m - zero_sample.drop_m;
    let actual_world_height_m =
        zero_distance_m * shooting_angle.sin() + shot_y_m * shooting_angle.cos();
    assert_relative_eq!(actual_world_height_m, target_world_height_m, epsilon = 1e-4);
}

#[test]
fn segmented_wind_matches_direct_engine_and_partial_coverage_warns() {
    let mut value = request_value("rk4");
    value["wind"] = json!({
        "segments": [
            {
                "until_distance_m": 60.0,
                "speed_mps": 2.0,
                "direction_from_rad": 0.0,
                "vertical_speed_mps": 0.1
            },
            {
                "until_distance_m": 120.0,
                "speed_mps": 6.0,
                "direction_from_rad": std::f64::consts::FRAC_PI_2,
                "vertical_speed_mps": -0.2
            }
        ]
    });
    let success = solve_v1(decode(&value)).expect("segmented solve");
    assert!(success
        .warnings
        .iter()
        .any(|warning| warning.code == "partial_wind_coverage"));

    let mut direct = direct_solver(SolverMethodV1::Rk4);
    direct.set_wind_segments(vec![
        WindSegment {
            speed_kmh: 7.2,
            angle_deg: 0.0,
            until_m: 60.0,
            vertical_mps: 0.1,
        },
        WindSegment {
            speed_kmh: 21.6,
            angle_deg: 90.0,
            until_m: 120.0,
            vertical_mps: -0.2,
        },
    ]);
    let observations = direct
        .solve()
        .expect("direct segmented solve")
        .sample_observations(25.0, 10_000)
        .expect("direct samples");
    for (actual, expected) in success.samples.iter().zip(observations) {
        assert_relative_eq!(actual.distance_m, expected.distance_m, epsilon = 1e-12);
        assert_relative_eq!(actual.drop_m, expected.drop_m, epsilon = 1e-11);
        assert_relative_eq!(actual.windage_m, expected.windage_m, epsilon = 1e-11);
        assert_relative_eq!(actual.speed_mps, expected.speed_mps, epsilon = 1e-9);
    }
}

#[test]
fn early_ground_impact_reports_the_exact_terminal_once() {
    let mut value = request_value("rk45");
    value["shot"]["max_range_m"] = json!(500.0);
    value["shot"]["muzzle_angle_rad"] = json!(-0.05);
    value["shot"]["ground_threshold_m"] = json!(0.0);
    value["sampling"]["interval_m"] = json!(30.0);

    let success = solve_v1(decode(&value)).expect("early-ground solve");
    assert_eq!(
        success.summary.termination,
        TerminationReasonV1::GroundThreshold
    );
    assert!(success.summary.actual_range_m < 500.0);
    assert_eq!(
        success
            .samples
            .iter()
            .filter(|sample| sample.flags.contains(&SampleFlagV1::Terminal))
            .count(),
        1
    );
    let terminal = success.samples.last().expect("terminal sample");
    assert!(terminal.flags.contains(&SampleFlagV1::GroundThreshold));
    assert_relative_eq!(
        terminal.distance_m,
        success.summary.actual_range_m,
        epsilon = 1e-12
    );

    let mut direct_inputs = direct_inputs(SolverMethodV1::Rk45);
    direct_inputs.target_distance = 500.0;
    direct_inputs.muzzle_angle = -0.05;
    direct_inputs.ground_threshold = 0.0;
    let mut direct = TrajectorySolver::new(
        direct_inputs,
        WindConditions {
            speed: 4.0,
            direction: std::f64::consts::FRAC_PI_2,
            vertical_speed: 0.25,
        },
        AtmosphericConditions {
            temperature: 20.0,
            pressure: 1000.0,
            humidity: 40.0,
            altitude: 0.0,
        },
    );
    direct.set_max_range(500.0);
    direct.set_time_step(0.001);
    let observations = direct
        .solve()
        .expect("direct early-ground solve")
        .sample_observations(30.0, 10_000)
        .expect("direct early-ground observations");
    assert_eq!(success.samples.len(), observations.len());
    for (actual, expected) in success.samples.iter().zip(observations) {
        assert_relative_eq!(actual.distance_m, expected.distance_m, epsilon = 1e-12);
        assert_relative_eq!(actual.time_s, expected.time_s, epsilon = 1e-12);
        assert_relative_eq!(actual.speed_mps, expected.speed_mps, epsilon = 1e-10);
        assert_relative_eq!(actual.drop_m, expected.drop_m, epsilon = 1e-12);
        assert_relative_eq!(actual.windage_m, expected.windage_m, epsilon = 1e-12);
    }
}

#[test]
fn sample_ceiling_uses_actual_reached_range_and_keeps_exact_boundary() {
    let mut exact_value = request_value("rk45");
    exact_value["shot"]["max_range_m"] = json!(99.99);
    exact_value["sampling"]["interval_m"] = json!(0.01);
    let exact = solve_v1(decode(&exact_value)).expect("exact 10,000-sample response");
    assert_eq!(exact.samples.len(), 10_000);

    let mut oversized_value = exact_value.clone();
    oversized_value["shot"]["max_range_m"] = json!(100.0);
    let oversized = solve_v1(decode(&oversized_value))
        .expect_err("10,001 actual observations must be rejected");
    assert_eq!(
        oversized.error.code,
        ballistics_engine::solve_json::SolveErrorCodeV1::ResourceLimit
    );
    assert_eq!(oversized.error.path(), Some("$.sampling.interval_m"));

    let mut early_ground_value = request_value("rk45");
    early_ground_value["shot"]["max_range_m"] = json!(500.0);
    early_ground_value["shot"]["muzzle_angle_rad"] = json!(-0.05);
    early_ground_value["shot"]["ground_threshold_m"] = json!(0.0);
    early_ground_value["sampling"]["interval_m"] = json!(0.01);
    let early_ground = solve_v1(decode(&early_ground_value))
        .expect("fine grid is valid when actual early-ground response stays under the cap");
    assert_eq!(
        early_ground.summary.termination,
        TerminationReasonV1::GroundThreshold
    );
    assert!(early_ground.samples.len() < 10_000);
}

/// Every field solve_v1 reads from the RAW request must have a resolved echo.
/// This is the guard that keeps Phase 0 true as the schema grows: if you add a
/// request field and consume it in solve_v1, add it here and to the resolved DTO.
#[test]
fn resolved_request_echoes_every_consumed_raw_field() {
    let json = serde_json::json!({
        "schema_version": 1,
        "projectile": {"mass_kg": 0.0113, "diameter_m": 0.00782, "drag_model": "G7",
                       "ballistic_coefficient": 0.243},
        "rifle": {"muzzle_velocity_mps": 823.0, "sight_height_m": 0.05,
                  "sight_offset_lateral_m": 0.012},
        "shot": {"max_range_m": 900.0, "zero_distance_m": 100.0,
                 "zero_poi_up_m": 0.01, "zero_poi_right_m": -0.005,
                 "drops_reference": "target", "shot_azimuth_rad": 0.0},
        "atmosphere": {"pressure_pa": 90000.0, "pressure_reference": "absolute"},
        "wind": {"speed_mps": 3.0, "direction_from_rad": 1.57, "wind_reference": "compass"},
        "solver": {}, "effects": {}, "sampling": {}
    });
    let req = ballistics_engine::decode_solve_request_v1(&json.to_string()).expect("decode");
    let ok = ballistics_engine::solve_v1(req).expect("solve");
    let r = &ok.resolved_request;

    assert_eq!(r.rifle.sight_offset_lateral_m, Some(0.012));
    assert_eq!(r.shot.zero_poi_up_m, Some(0.01));
    assert_eq!(r.shot.zero_poi_right_m, Some(-0.005));
    assert!(r.shot.drops_reference.is_some(), "drops_reference must be echoed");
    assert!(
        r.atmosphere.pressure_reference.is_some(),
        "pressure_reference must be echoed"
    );
}
