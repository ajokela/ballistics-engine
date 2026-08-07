use ballistics_engine::solve_json::{
    decode_solve_request_v1, DragModelV1, ResolvedAtmosphereV1, ResolvedConstantWindV1,
    ResolvedEffectsV1, ResolvedProjectileV1, ResolvedRifleV1, ResolvedSamplingV1, ResolvedShotV1,
    ResolvedSolveRequestV1, ResolvedSolverV1, ResolvedWindSegmentV1, ResolvedWindV1, SampleFlagV1,
    SchemaVersionV1, SolveErrorCodeV1, SolveErrorEnvelopeV1, SolveErrorLocationErrorV1,
    SolveErrorV1, SolveNoticeV1, SolveRequestV1, SolveSuccessV1, SolveSummaryV1, SuccessStatusV1,
    TerminationReasonV1, TrajectorySampleV1, DRAG_MODEL_WIRE_NAMES_V1,
    MAX_SOLVE_JSON_SAMPLES_V1, SOLVE_JSON_SCHEMA_VERSION_V1,
};
use serde_json::{json, Value};

fn request_value() -> Value {
    json!({
        "schema_version": 1,
        "projectile": {
            "mass_kg": 0.01134,
            "diameter_m": 0.00671,
            "length_m": 0.031,
            "drag_model": "G7",
            "ballistic_coefficient": 0.243
        },
        "rifle": {
            "muzzle_velocity_mps": 823.0,
            "sight_height_m": 0.0381,
            "twist_rate_m_per_turn": 0.2032,
            "twist_direction": "right"
        },
        "shot": {
            "max_range_m": 1000.0,
            "zero_distance_m": 100.0,
            "shooting_angle_rad": 0.0,
            "cant_angle_rad": 0.0
        },
        "atmosphere": {
            "altitude_m": 250.0,
            "temperature_k": 288.15,
            "pressure_pa": 101325.0,
            "relative_humidity": 0.5,
            "latitude_rad": std::f64::consts::FRAC_PI_4
        },
        "wind": {
            "speed_mps": 4.4704,
            "direction_from_rad": std::f64::consts::FRAC_PI_2,
            "vertical_speed_mps": 0.0
        },
        "solver": {
            "method": "rk45",
            "time_step_s": 0.001
        },
        "effects": {
            "magnus": false,
            "coriolis": true,
            "enhanced_spin_drift": true
        },
        "sampling": {
            "interval_m": 10.0
        }
    })
}

fn decode(value: &Value) -> Result<SolveRequestV1, SolveErrorEnvelopeV1> {
    decode_solve_request_v1(&serde_json::to_string(value).expect("serialize fixture"))
}

fn resolved_request_for(request: &SolveRequestV1) -> ResolvedSolveRequestV1 {
    let wind = if let Some(segments) = &request.wind.segments {
        ResolvedWindV1::Segmented(ballistics_engine::solve_json::ResolvedSegmentedWindV1 {
            segments: segments
                .iter()
                .map(|segment| ResolvedWindSegmentV1 {
                    until_distance_m: segment.until_distance_m,
                    speed_mps: segment.speed_mps,
                    direction_from_rad: segment.direction_from_rad,
                    vertical_speed_mps: segment.vertical_speed_mps.unwrap_or(0.0),
                })
                .collect(),
            wind_reference: request.wind.wind_reference,
        })
    } else {
        ResolvedWindV1::Constant(ResolvedConstantWindV1 {
            speed_mps: request.wind.speed_mps.unwrap_or(0.0),
            direction_from_rad: request.wind.direction_from_rad.unwrap_or(0.0),
            vertical_speed_mps: request.wind.vertical_speed_mps.unwrap_or(0.0),
            wind_reference: request.wind.wind_reference,
        })
    };

    ResolvedSolveRequestV1 {
        schema_version: SchemaVersionV1,
        projectile: ResolvedProjectileV1 {
            mass_kg: request.projectile.mass_kg,
            diameter_m: request.projectile.diameter_m,
            length_m: request.projectile.length_m,
            drag_model: request.projectile.drag_model,
            ballistic_coefficient: request.projectile.ballistic_coefficient,
        },
        rifle: ResolvedRifleV1 {
            muzzle_velocity_mps: request.rifle.muzzle_velocity_mps,
            sight_height_m: request.rifle.sight_height_m.unwrap_or(0.05),
            muzzle_height_m: request.rifle.muzzle_height_m.unwrap_or(0.0),
            twist_rate_m_per_turn: request.rifle.twist_rate_m_per_turn.unwrap_or(0.3048),
            twist_direction: request.rifle.twist_direction.unwrap_or_default(),
            sight_offset_lateral_m: request.rifle.sight_offset_lateral_m,
        },
        shot: ResolvedShotV1 {
            max_range_m: request.shot.max_range_m,
            zero_distance_m: request.shot.zero_distance_m,
            muzzle_angle_rad: request.shot.muzzle_angle_rad.unwrap_or(0.0),
            aim_azimuth_rad: request.shot.aim_azimuth_rad.unwrap_or(0.0),
            shot_azimuth_rad: request.shot.shot_azimuth_rad.unwrap_or(0.0),
            shooting_angle_rad: request.shot.shooting_angle_rad.unwrap_or(0.0),
            cant_angle_rad: request.shot.cant_angle_rad.unwrap_or(0.0),
            target_height_m: request.shot.target_height_m.unwrap_or(0.0),
            ground_threshold_m: request.shot.ground_threshold_m.unwrap_or(-100.0),
            zero_poi_up_m: request.shot.zero_poi_up_m,
            zero_poi_right_m: request.shot.zero_poi_right_m,
            drops_reference: request.shot.drops_reference,
        },
        atmosphere: ResolvedAtmosphereV1 {
            altitude_m: request.atmosphere.altitude_m.unwrap_or(0.0),
            temperature_k: request.atmosphere.temperature_k.unwrap_or(288.15),
            pressure_pa: request.atmosphere.pressure_pa.unwrap_or(101_325.0),
            relative_humidity: request.atmosphere.relative_humidity.unwrap_or(0.5),
            latitude_rad: request.atmosphere.latitude_rad,
            pressure_reference: request.atmosphere.pressure_reference,
        },
        wind,
        solver: ResolvedSolverV1 {
            method: request.solver.method.unwrap_or_default(),
            time_step_s: request.solver.time_step_s.unwrap_or(0.001),
        },
        effects: ResolvedEffectsV1 {
            magnus: request.effects.magnus.unwrap_or(false),
            coriolis: request.effects.coriolis.unwrap_or(false),
            enhanced_spin_drift: request.effects.enhanced_spin_drift.unwrap_or(false),
        },
        sampling: ResolvedSamplingV1 {
            interval_m: request.sampling.interval_m.unwrap_or(10.0),
        },
        reticle: request.reticle.clone(),
    }
}

fn success_for(request: SolveRequestV1) -> SolveSuccessV1 {
    let resolved_request = resolved_request_for(&request);
    success_for_resolved(resolved_request)
}

fn success_for_resolved(resolved_request: ResolvedSolveRequestV1) -> SolveSuccessV1 {
    SolveSuccessV1 {
        schema_version: SchemaVersionV1,
        engine_version: "0.24.1".to_owned(),
        status: SuccessStatusV1::Ok,
        resolved_request,
        assumptions: vec![SolveNoticeV1 {
            code: "resolved_inputs".to_owned(),
            message: "All omitted inputs were resolved before the solve.".to_owned(),
            path: None,
        }],
        warnings: Vec::new(),
        summary: SolveSummaryV1 {
            actual_range_m: 1000.0,
            maximum_height_m: 3.2,
            time_of_flight_s: 1.6,
            terminal_speed_mps: 360.0,
            terminal_energy_j: 734.8,
            stability_factor: Some(1.5),
            spin_drift_m: Some(0.12),
            equivalent_horizontal_range_m: None,
            termination: TerminationReasonV1::MaxRange,
        },
        samples: vec![TrajectorySampleV1 {
            distance_m: 1000.0,
            time_s: 1.6,
            speed_mps: 360.0,
            energy_j: 734.8,
            drop_m: 8.1,
            windage_m: 0.43,
            mach: 1.06,
            flags: vec![SampleFlagV1::Transonic, SampleFlagV1::Terminal],
        }],
        // MBA-1361: absent unless the request carried a `reticle` block.
        reticle_hold: None,
    }
}

#[test]
fn omitted_defaults_remain_absent_across_request_round_trip() {
    let mut value = request_value();
    value["rifle"] = json!({ "muzzle_velocity_mps": 823.0 });
    value["shot"] = json!({ "max_range_m": 1000.0 });
    value["atmosphere"] = json!({});
    value["wind"] = json!({});
    value["solver"] = json!({});
    value["effects"] = json!({});
    value["sampling"] = json!({});

    let request = decode(&value).expect("valid request");
    assert_eq!(request.rifle.sight_height_m, None);
    assert_eq!(request.rifle.muzzle_height_m, None);
    assert_eq!(request.rifle.twist_rate_m_per_turn, None);
    assert_eq!(request.rifle.twist_direction, None);
    assert_eq!(request.shot.aim_azimuth_rad, None);
    assert_eq!(request.shot.shot_azimuth_rad, None);
    assert_eq!(request.shot.shooting_angle_rad, None);
    assert_eq!(request.shot.cant_angle_rad, None);
    assert_eq!(request.shot.target_height_m, None);
    assert_eq!(request.shot.ground_threshold_m, None);
    assert_eq!(request.atmosphere.altitude_m, None);
    assert_eq!(request.atmosphere.temperature_k, None);
    assert_eq!(request.atmosphere.pressure_pa, None);
    assert_eq!(request.atmosphere.relative_humidity, None);
    assert_eq!(request.wind.segments, None);
    assert_eq!(request.solver.method, None);
    assert_eq!(request.solver.time_step_s, None);
    assert_eq!(request.effects.magnus, None);
    assert_eq!(request.effects.coriolis, None);
    assert_eq!(request.effects.enhanced_spin_drift, None);
    assert_eq!(request.sampling.interval_m, None);

    let encoded = serde_json::to_value(&request).expect("serialize request");
    for section in ["atmosphere", "wind", "solver", "effects", "sampling"] {
        assert_eq!(encoded[section], json!({}));
    }
    assert_eq!(encoded["rifle"], json!({ "muzzle_velocity_mps": 823.0 }));
    assert_eq!(encoded["shot"], json!({ "max_range_m": 1000.0 }));

    let decoded = decode(&encoded).expect("deserialize request");
    assert_eq!(decoded, request);
}

#[test]
fn explicit_values_equal_to_defaults_remain_present() {
    let mut omitted_value = request_value();
    omitted_value["rifle"] = json!({ "muzzle_velocity_mps": 823.0 });
    omitted_value["shot"] = json!({ "max_range_m": 1000.0 });
    omitted_value["atmosphere"] = json!({});
    omitted_value["solver"] = json!({});
    omitted_value["effects"] = json!({});
    omitted_value["sampling"] = json!({});
    let omitted = decode(&omitted_value).expect("omitted defaults are valid");

    let mut explicit_value = omitted_value;
    explicit_value["rifle"] = json!({
        "muzzle_velocity_mps": 823.0,
        "sight_height_m": 0.05,
        "muzzle_height_m": 0.0,
        "twist_rate_m_per_turn": 0.3048,
        "twist_direction": "right"
    });
    explicit_value["shot"] = json!({
        "max_range_m": 1000.0,
        "aim_azimuth_rad": 0.0,
        "shot_azimuth_rad": 0.0,
        "shooting_angle_rad": 0.0,
        "cant_angle_rad": 0.0,
        "target_height_m": 0.0,
        "ground_threshold_m": -100.0
    });
    explicit_value["atmosphere"] = json!({
        "altitude_m": 0.0,
        "temperature_k": 288.15,
        "pressure_pa": 101325.0,
        "relative_humidity": 0.5
    });
    explicit_value["solver"] = json!({ "method": "rk45", "time_step_s": 0.001 });
    explicit_value["effects"] = json!({
        "magnus": false,
        "coriolis": false,
        "enhanced_spin_drift": false
    });
    explicit_value["sampling"] = json!({ "interval_m": 10.0 });

    let explicit = decode(&explicit_value).expect("explicit defaults are valid");
    assert_ne!(explicit, omitted);
    assert_eq!(explicit.rifle.sight_height_m, Some(0.05));
    assert_eq!(explicit.shot.ground_threshold_m, Some(-100.0));
    assert_eq!(explicit.atmosphere.altitude_m, Some(0.0));
    assert_eq!(explicit.solver.method, Some(Default::default()));
    assert_eq!(explicit.effects.magnus, Some(false));
    assert_eq!(explicit.sampling.interval_m, Some(10.0));

    let wire = serde_json::to_value(&explicit).expect("serialize explicit request");
    assert_eq!(wire["rifle"]["sight_height_m"], json!(0.05));
    assert_eq!(wire["shot"]["ground_threshold_m"], json!(-100.0));
    assert_eq!(wire["atmosphere"]["altitude_m"], json!(0.0));
    assert_eq!(wire["solver"]["method"], json!("rk45"));
    assert_eq!(wire["effects"]["magnus"], json!(false));
    assert_eq!(wire["sampling"]["interval_m"], json!(10.0));
}

#[test]
fn omitted_and_explicit_standard_sea_level_values_remain_distinct_at_altitude() {
    let mut omitted_value = request_value();
    omitted_value["atmosphere"] = json!({
        "altitude_m": 1500.0,
        "relative_humidity": 0.5
    });
    let omitted = decode(&omitted_value).expect("omitted atmosphere values are valid");
    assert_eq!(omitted.atmosphere.temperature_k, None);
    assert_eq!(omitted.atmosphere.pressure_pa, None);
    let omitted_wire = serde_json::to_value(&omitted).expect("serialize omitted request");
    assert!(omitted_wire["atmosphere"].get("temperature_k").is_none());
    assert!(omitted_wire["atmosphere"].get("pressure_pa").is_none());

    let mut explicit_value = omitted_value;
    explicit_value["atmosphere"]["temperature_k"] = json!(288.15);
    explicit_value["atmosphere"]["pressure_pa"] = json!(101_325.0);
    let explicit = decode(&explicit_value).expect("explicit station values are valid");
    assert_eq!(explicit.atmosphere.temperature_k, Some(288.15));
    assert_eq!(explicit.atmosphere.pressure_pa, Some(101_325.0));
    assert_ne!(explicit, omitted);
}

#[test]
fn resolved_request_materializes_atmosphere_values_in_success_output() {
    let mut value = request_value();
    value["atmosphere"] = json!({ "altitude_m": 1500.0 });
    let request = decode(&value).expect("request with omitted station values");
    assert_eq!(request.atmosphere.temperature_k, None);
    assert_eq!(request.atmosphere.pressure_pa, None);

    // MBA-1302 performs the ICAO-at-altitude resolution and must populate both fields before
    // constructing a success envelope. These representative values exercise the wire invariant.
    let mut resolved = resolved_request_for(&request);
    resolved.atmosphere.temperature_k = 278.4;
    resolved.atmosphere.pressure_pa = 84_556.0;
    let wire =
        serde_json::to_value(success_for_resolved(resolved)).expect("serialize success output");
    assert_eq!(
        wire["resolved_request"]["atmosphere"]["temperature_k"],
        json!(278.4)
    );
    assert_eq!(
        wire["resolved_request"]["atmosphere"]["pressure_pa"],
        json!(84_556.0)
    );
}

#[test]
fn success_resolved_request_materializes_every_documented_default() {
    let mut value = request_value();
    value["rifle"] = json!({ "muzzle_velocity_mps": 823.0 });
    value["shot"] = json!({ "max_range_m": 1000.0 });
    value["atmosphere"] = json!({});
    value["wind"] = json!({});
    value["solver"] = json!({});
    value["effects"] = json!({});
    value["sampling"] = json!({});

    let request = decode(&value).expect("request with omitted defaults");
    let wire = serde_json::to_value(success_for(request)).expect("serialize success");
    let resolved = &wire["resolved_request"];

    assert_eq!(
        resolved["rifle"],
        json!({
            "muzzle_velocity_mps": 823.0,
            "sight_height_m": 0.05,
            "muzzle_height_m": 0.0,
            "twist_rate_m_per_turn": 0.3048,
            "twist_direction": "right"
        })
    );
    assert_eq!(
        resolved["shot"],
        json!({
            "max_range_m": 1000.0,
            "muzzle_angle_rad": 0.0,
            "aim_azimuth_rad": 0.0,
            "shot_azimuth_rad": 0.0,
            "shooting_angle_rad": 0.0,
            "cant_angle_rad": 0.0,
            "target_height_m": 0.0,
            "ground_threshold_m": -100.0
        })
    );
    assert_eq!(
        resolved["atmosphere"],
        json!({
            "altitude_m": 0.0,
            "temperature_k": 288.15,
            "pressure_pa": 101325.0,
            "relative_humidity": 0.5
        })
    );
    assert_eq!(
        resolved["wind"],
        json!({
            "speed_mps": 0.0,
            "direction_from_rad": 0.0,
            "vertical_speed_mps": 0.0
        })
    );
    assert_eq!(
        resolved["solver"],
        json!({ "method": "rk45", "time_step_s": 0.001 })
    );
    assert_eq!(
        resolved["effects"],
        json!({
            "magnus": false,
            "coriolis": false,
            "enhanced_spin_drift": false
        })
    );
    assert_eq!(resolved["sampling"], json!({ "interval_m": 10.0 }));
}

#[test]
fn resolved_zeroing_records_intent_and_effective_muzzle_angle() {
    let request = decode(&request_value()).expect("zero-distance request");
    let mut resolved = resolved_request_for(&request);
    resolved.shot.muzzle_angle_rad = 0.0123;

    let wire = serde_json::to_value(success_for_resolved(resolved)).expect("serialize success");
    assert_eq!(
        wire["resolved_request"]["shot"]["zero_distance_m"],
        json!(100.0)
    );
    assert_eq!(
        wire["resolved_request"]["shot"]["muzzle_angle_rad"],
        json!(0.0123)
    );
}

#[test]
fn wind_segment_and_vertical_omissions_survive_request_round_trip() {
    let mut absent_value = request_value();
    absent_value["wind"] = json!({});
    let absent = decode(&absent_value).expect("omitted wind is valid");
    assert_eq!(absent.wind.segments, None);

    let mut empty_value = absent_value;
    empty_value["wind"] = json!({ "segments": [] });
    let empty = decode(&empty_value).expect("explicit segment array is structurally valid");
    assert_eq!(empty.wind.segments, Some(Vec::new()));
    assert_ne!(empty, absent);
    assert_eq!(
        serde_json::to_value(&empty).unwrap()["wind"]["segments"],
        json!([])
    );

    let mut segments_value = request_value();
    segments_value["wind"] = json!({
        "segments": [{
            "until_distance_m": 1000.0,
            "speed_mps": 2.0,
            "direction_from_rad": 0.0
        }]
    });
    let omitted_vertical = decode(&segments_value).expect("omitted segment vertical speed");
    assert_eq!(
        omitted_vertical.wind.segments.as_ref().unwrap()[0].vertical_speed_mps,
        None
    );

    segments_value["wind"]["segments"][0]["vertical_speed_mps"] = json!(0.0);
    let explicit_vertical = decode(&segments_value).expect("explicit segment vertical speed");
    assert_eq!(
        explicit_vertical.wind.segments.as_ref().unwrap()[0].vertical_speed_mps,
        Some(0.0)
    );
    assert_ne!(explicit_vertical, omitted_vertical);

    let resolved_wire =
        serde_json::to_value(success_for(omitted_vertical)).expect("serialize resolved segments");
    assert_eq!(
        resolved_wire["resolved_request"]["wind"]["segments"][0]["vertical_speed_mps"],
        json!(0.0)
    );
}

#[test]
fn success_and_error_dtos_round_trip() {
    let request = decode(&request_value()).expect("valid request");
    let success = success_for(request);

    let success_json = serde_json::to_string(&success).expect("serialize success");
    let success_again: SolveSuccessV1 =
        serde_json::from_str(&success_json).expect("deserialize success");
    assert_eq!(success_again, success);

    let failure = SolveErrorEnvelopeV1::new(
        SolveErrorV1::new(SolveErrorCodeV1::InvalidValue, "mass must be positive")
            .at_path("$.projectile.mass_kg"),
    );
    let failure_json = serde_json::to_string(&failure).expect("serialize failure");
    let failure_again: SolveErrorEnvelopeV1 =
        serde_json::from_str(&failure_json).expect("deserialize failure");
    assert_eq!(failure_again, failure);
}

#[test]
fn success_response_sample_limit_accepts_exactly_ten_thousand() {
    assert_eq!(MAX_SOLVE_JSON_SAMPLES_V1, 10_000);

    let request = decode(&request_value()).expect("valid request");
    let mut success = success_for(request);
    let sample = success.samples[0].clone();
    success.samples.resize(MAX_SOLVE_JSON_SAMPLES_V1, sample);

    success
        .validate_for_serialization()
        .expect("the exact response sample limit must be accepted");
    serde_json::to_vec(&success).expect("an at-limit response remains serializable");
}

#[test]
fn success_response_sample_limit_rejects_ten_thousand_and_one() {
    let request = decode(&request_value()).expect("valid request");
    let mut success = success_for(request);
    let sample = success.samples[0].clone();
    success
        .samples
        .resize(MAX_SOLVE_JSON_SAMPLES_V1 + 1, sample);

    let error = success
        .validate_for_serialization()
        .expect_err("an oversized response must fail before serialization");
    assert_eq!(error.error.code, SolveErrorCodeV1::ResourceLimit);
    assert_eq!(error.error.path(), Some("$.sampling.interval_m"));
    assert!(error.error.message.contains("10000"));
    assert!(error.error.message.contains("10001"));

    let serializer_error = serde_json::to_vec(&success)
        .expect_err("the wire type must fail closed if a caller skips explicit validation");
    assert!(serializer_error.to_string().contains("10000"));
    assert!(serializer_error.to_string().contains("10001"));

    let wire = serde_json::to_value(error).expect("serialize structured limit error");
    assert_eq!(wire["status"], json!("error"));
    assert_eq!(wire["error"]["code"], json!("resource_limit"));
}

#[test]
fn every_nested_request_object_rejects_unknown_fields_with_a_path() {
    for section in [
        "projectile",
        "rifle",
        "shot",
        "atmosphere",
        "wind",
        "solver",
        "effects",
        "sampling",
    ] {
        let mut value = request_value();
        value[section]["typo"] = json!(true);

        let error = decode(&value).expect_err("unknown field must fail");
        assert_eq!(error.error.code, SolveErrorCodeV1::UnknownField);
        let expected_path = format!("$.{section}.typo");
        assert_eq!(error.error.path(), Some(expected_path.as_str()));
    }
}

#[test]
fn nested_wind_segment_rejects_unknown_fields_with_an_indexed_path() {
    let mut value = request_value();
    value["wind"] = json!({
        "segments": [{
            "until_distance_m": 500.0,
            "speed_mps": 2.0,
            "direction_from_rad": 0.0,
            "dirction_from_rad": 1.0
        }]
    });

    let error = decode(&value).expect_err("unknown segment field must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::UnknownField);
    assert_eq!(
        error.error.path(),
        Some("$.wind.segments[0].dirction_from_rad")
    );
}

#[test]
fn missing_required_nested_field_has_a_stable_code_and_path() {
    let mut value = request_value();
    value["projectile"]
        .as_object_mut()
        .expect("projectile object")
        .remove("mass_kg");

    let error = decode(&value).expect_err("missing field must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::MissingField);
    assert_eq!(error.error.path(), Some("$.projectile.mass_kg"));
}

#[test]
fn wrong_scalar_type_has_the_exact_nested_path() {
    let mut value = request_value();
    value["projectile"]["mass_kg"] = json!("eleven grams");

    let error = decode(&value).expect_err("wrong scalar type must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(error.error.path(), Some("$.projectile.mass_kg"));
    assert_eq!(error.error.line(), None);
    assert_eq!(error.error.column(), None);
}

#[test]
fn explicit_null_is_not_treated_as_field_omission() {
    let fields = [
        ("projectile", "length_m"),
        ("rifle", "sight_height_m"),
        ("rifle", "muzzle_height_m"),
        ("rifle", "twist_rate_m_per_turn"),
        ("rifle", "twist_direction"),
        ("shot", "zero_distance_m"),
        ("shot", "muzzle_angle_rad"),
        ("shot", "aim_azimuth_rad"),
        ("shot", "shot_azimuth_rad"),
        ("shot", "shooting_angle_rad"),
        ("shot", "cant_angle_rad"),
        ("shot", "target_height_m"),
        ("shot", "ground_threshold_m"),
        ("atmosphere", "altitude_m"),
        ("atmosphere", "temperature_k"),
        ("atmosphere", "pressure_pa"),
        ("atmosphere", "relative_humidity"),
        ("atmosphere", "latitude_rad"),
        ("wind", "speed_mps"),
        ("wind", "direction_from_rad"),
        ("wind", "vertical_speed_mps"),
        ("wind", "segments"),
        ("solver", "method"),
        ("solver", "time_step_s"),
        ("effects", "magnus"),
        ("effects", "coriolis"),
        ("effects", "enhanced_spin_drift"),
        ("sampling", "interval_m"),
    ];

    for (section, field) in fields {
        let mut value = request_value();
        value[section][field] = Value::Null;

        serde_json::from_value::<SolveRequestV1>(value.clone())
            .expect_err("direct serde must reject null too");
        let error = decode(&value).expect_err("null must not mean omitted");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        let expected_path = format!("$.{section}.{field}");
        assert_eq!(error.error.path(), Some(expected_path.as_str()));
    }

    let mut value = request_value();
    value["wind"] = json!({
        "segments": [{
            "until_distance_m": 1000.0,
            "speed_mps": 2.0,
            "direction_from_rad": 0.0,
            "vertical_speed_mps": null
        }]
    });
    serde_json::from_value::<SolveRequestV1>(value.clone())
        .expect_err("direct serde must reject null segment vertical speed");
    let error = decode(&value).expect_err("null segment vertical speed must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
    assert_eq!(
        error.error.path(),
        Some("$.wind.segments[0].vertical_speed_mps")
    );
}

#[test]
fn invalid_enum_values_have_stable_paths() {
    for (section, field, invalid, expected_path) in [
        ("projectile", "drag_model", "G9", "$.projectile.drag_model"),
        ("projectile", "drag_model", "g2", "$.projectile.drag_model"),
        (
            "rifle",
            "twist_direction",
            "clockwise",
            "$.rifle.twist_direction",
        ),
        ("solver", "method", "verlet", "$.solver.method"),
    ] {
        let mut value = request_value();
        value[section][field] = json!(invalid);

        let error = decode(&value).expect_err("invalid enum must fail");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidValue);
        assert_eq!(error.error.path(), Some(expected_path));
    }
}

#[test]
fn unsupported_schema_version_is_distinct_from_an_invalid_value() {
    let mut value = request_value();
    value["schema_version"] = json!(2);
    value["future_v2_section"] = json!({ "new_field": true });

    let error = decode(&value).expect_err("unsupported version must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::UnsupportedSchemaVersion);
    assert_eq!(error.error.path(), Some("$.schema_version"));

    value
        .as_object_mut()
        .expect("request object")
        .remove("future_v2_section");
    let direct_error = serde_json::from_value::<SolveRequestV1>(value)
        .expect_err("direct serde decoding also rejects unsupported versions");
    assert!(direct_error
        .to_string()
        .contains("unsupported schema_version 2"));
}

#[test]
fn every_other_integer_schema_version_is_unsupported() {
    for version in [-1, 0] {
        let mut value = request_value();
        value["schema_version"] = json!(version);

        let error = decode(&value).expect_err("non-v1 integer must fail");
        assert_eq!(error.error.code, SolveErrorCodeV1::UnsupportedSchemaVersion);
        assert_eq!(error.error.path(), Some("$.schema_version"));
    }

    let direct_error = serde_json::from_value::<SchemaVersionV1>(json!(-1))
        .expect_err("the invariant type rejects negative versions");
    assert!(direct_error
        .to_string()
        .contains("unsupported schema_version -1"));
}

#[test]
fn schema_version_invariant_always_serializes_as_integer_one() {
    assert_eq!(SchemaVersionV1.get(), SOLVE_JSON_SCHEMA_VERSION_V1);
    assert_eq!(serde_json::to_value(SchemaVersionV1).unwrap(), json!(1));

    let request = decode(&request_value()).expect("valid request");
    assert_eq!(
        serde_json::to_value(&request).unwrap()["schema_version"],
        json!(1)
    );
    assert_eq!(
        serde_json::to_value(success_for(request)).unwrap()["schema_version"],
        json!(1)
    );
    assert_eq!(
        serde_json::to_value(SolveErrorEnvelopeV1::new(SolveErrorV1::new(
            SolveErrorCodeV1::InternalError,
            "contained failure",
        )))
        .unwrap()["schema_version"],
        json!(1)
    );
}

#[test]
fn all_reference_drag_models_decode_and_round_trip_exact_wire_names() {
    let cases = [
        ("G1", DragModelV1::G1),
        ("G2", DragModelV1::G2),
        ("G5", DragModelV1::G5),
        ("G6", DragModelV1::G6),
        ("G7", DragModelV1::G7),
        ("G8", DragModelV1::G8),
        ("GI", DragModelV1::GI),
        ("GS", DragModelV1::GS),
        ("RA4", DragModelV1::RA4),
    ];
    assert_eq!(
        cases.map(|(wire_name, _)| wire_name),
        DRAG_MODEL_WIRE_NAMES_V1
    );

    for (wire_name, expected) in cases {
        let mut value = request_value();
        value["projectile"]["drag_model"] = json!(wire_name);

        let request = decode(&value).expect("every built-in reference model must decode");
        assert_eq!(request.projectile.drag_model, expected);

        let encoded = serde_json::to_value(&request).expect("serialize request");
        assert_eq!(encoded["projectile"]["drag_model"], json!(wire_name));

        let direct: DragModelV1 =
            serde_json::from_value(json!(wire_name)).expect("direct enum decode");
        assert_eq!(direct, expected);
        assert_eq!(serde_json::to_value(direct).unwrap(), json!(wire_name));
    }
}

#[test]
fn magnus_and_enhanced_spin_drift_are_a_structural_conflict() {
    let mut value = request_value();
    value["effects"]["magnus"] = json!(true);
    value["effects"]["enhanced_spin_drift"] = json!(true);

    let error = decode(&value).expect_err("mutually suppressing effects must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::ConflictingFields);
    assert_eq!(error.error.path(), Some("$.effects"));
}

#[test]
fn termination_and_ground_sample_names_match_engine_metadata() {
    for (termination, wire_name) in [
        (TerminationReasonV1::MaxRange, "max_range"),
        (TerminationReasonV1::GroundThreshold, "ground_threshold"),
        (TerminationReasonV1::TimeLimit, "time_limit"),
        (TerminationReasonV1::VelocityFloor, "velocity_floor"),
    ] {
        assert_eq!(serde_json::to_value(termination).unwrap(), json!(wire_name));
    }
    assert_eq!(
        serde_json::to_value(SampleFlagV1::GroundThreshold).unwrap(),
        json!("ground_threshold")
    );
    for retired_name in ["ground_impact", "velocity_limit", "solver_limit"] {
        assert!(serde_json::from_value::<TerminationReasonV1>(json!(retired_name)).is_err());
    }
}

#[test]
fn error_location_builders_replace_the_other_location_kind() {
    let source = SolveErrorV1::new(SolveErrorCodeV1::InvalidJson, "bad json")
        .at_path("$")
        .at_location(2, 7)
        .expect("one-based source location");
    assert_eq!(source.path(), None);
    assert_eq!(source.line(), Some(2));
    assert_eq!(source.column(), Some(7));
    let source_wire = serde_json::to_value(&source).expect("serialize source error");
    assert_eq!(source_wire["path"], Value::Null);
    assert_eq!(source_wire["line"], json!(2));
    assert_eq!(source_wire["column"], json!(7));

    let path = source.at_path("$.projectile.mass_kg");
    assert_eq!(path.path(), Some("$.projectile.mass_kg"));
    assert_eq!(path.line(), None);
    assert_eq!(path.column(), None);
    let path_wire = serde_json::to_value(path).expect("serialize path error");
    assert_eq!(path_wire["path"], json!("$.projectile.mass_kg"));
    assert_eq!(path_wire["line"], Value::Null);
    assert_eq!(path_wire["column"], Value::Null);
}

#[test]
fn error_location_builder_rejects_zero_coordinates() {
    let error = SolveErrorV1::new(SolveErrorCodeV1::InvalidJson, "bad json");
    assert_eq!(
        error.clone().at_location(0, 1).unwrap_err(),
        SolveErrorLocationErrorV1::ZeroLine
    );
    assert_eq!(
        error.at_location(1, 0).unwrap_err(),
        SolveErrorLocationErrorV1::ZeroColumn
    );
}

#[test]
fn error_wire_rejects_ambiguous_partial_and_zero_locations() {
    let cases = [
        json!({
            "code": "invalid_value",
            "message": "ambiguous",
            "path": "$.projectile.mass_kg",
            "line": 1,
            "column": 2
        }),
        json!({
            "code": "invalid_json",
            "message": "partial",
            "path": null,
            "line": 1,
            "column": null
        }),
        json!({
            "code": "invalid_json",
            "message": "partial",
            "path": null,
            "line": null,
            "column": 2
        }),
        json!({
            "code": "invalid_json",
            "message": "zero line",
            "path": null,
            "line": 0,
            "column": 1
        }),
        json!({
            "code": "invalid_json",
            "message": "zero column",
            "path": null,
            "line": 1,
            "column": 0
        }),
    ];

    for wire in cases {
        assert!(
            serde_json::from_value::<SolveErrorV1>(wire).is_err(),
            "invalid error location must be rejected"
        );
    }
}

#[test]
fn malformed_json_reports_source_location() {
    let error = decode_solve_request_v1("{\n  \"schema_version\": 1,\n}")
        .expect_err("trailing comma must fail");
    assert_eq!(error.error.code, SolveErrorCodeV1::InvalidJson);
    assert!(error.error.line().is_some_and(|line| line > 0));
    assert!(error.error.column().is_some_and(|column| column > 0));
    assert_eq!(error.error.path(), None);
}

#[test]
fn empty_or_newline_terminated_input_returns_an_error_instead_of_panicking() {
    for input in ["", "\n"] {
        let error = decode_solve_request_v1(input).expect_err("empty input must fail");
        assert_eq!(error.error.code, SolveErrorCodeV1::InvalidJson);
        assert!(error.error.line().is_some_and(|line| line > 0));
        assert!(error.error.column().is_some_and(|column| column > 0));
        assert_eq!(error.error.path(), None);
    }
}
