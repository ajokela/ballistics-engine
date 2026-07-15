use ballistics_engine::solve_json::{
    decode_solve_request_v1, SampleFlagV1, SolveErrorCodeV1, SolveRequestV1, TerminationReasonV1,
};
use ballistics_engine::{
    solve_v1, AtmosphericConditions, BallisticInputs, DragModel, TrajectorySolver, WindConditions,
};
use serde_json::Value;
use std::collections::BTreeSet;

const REPRESENTATIVE_REQUEST: &str =
    include_str!("fixtures/solve_json_v1/representative.request.json");
const REPRESENTATIVE_SUCCESS: &str =
    include_str!("fixtures/solve_json_v1/representative.success.json");
const VALIDATION_REQUEST: &str = include_str!("fixtures/solve_json_v1/validation.request.json");
const VALIDATION_ERROR: &str = include_str!("fixtures/solve_json_v1/validation.error.json");
const SCHEMA_REQUEST: &str = include_str!("fixtures/solve_json_v1/schema.request.json");
const SCHEMA_ERROR: &str = include_str!("fixtures/solve_json_v1/schema.error.json");
const RESOURCE_REQUEST: &str = include_str!("fixtures/solve_json_v1/resource.request.json");
const RESOURCE_ERROR: &str = include_str!("fixtures/solve_json_v1/resource.error.json");
const EARLY_REQUEST: &str = include_str!("fixtures/solve_json_v1/early-termination.request.json");
const EARLY_SUCCESS: &str = include_str!("fixtures/solve_json_v1/early-termination.success.json");

/// Fixture numerics are semantic values rather than serialized byte snapshots.
///
/// Every numeric comparison uses:
/// `abs(a-b) <= ABS_TOLERANCE + REL_TOLERANCE * max(abs(a), abs(b))`.
const ABS_TOLERANCE: f64 = 1.0e-10;
const REL_TOLERANCE: f64 = 1.0e-9;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum JsonShape {
    Null,
    Boolean,
    Number,
    String,
    Array(Vec<JsonShape>),
    Object(Vec<(String, JsonShape)>),
}

fn parse_fixture(document: &str) -> Value {
    serde_json::from_str(document).expect("checked solve-json fixture must be valid JSON")
}

fn evaluate_request(document: &str) -> Value {
    match decode_solve_request_v1(document) {
        Ok(request) => match solve_v1(request) {
            Ok(success) => serde_json::to_value(success).expect("success envelope serialization"),
            Err(error) => serde_json::to_value(error).expect("error envelope serialization"),
        },
        Err(error) => serde_json::to_value(error).expect("decode error envelope serialization"),
    }
}

/// Normalize a JSON document to unordered object keys and unique array element shapes.
///
/// Array cardinality and element order are intentionally excluded. All distinct element shapes
/// are retained, so an empty first sample's `flags` cannot hide the terminal sample's string
/// flags later in the response.
fn normalized_shape(value: &Value) -> JsonShape {
    match value {
        Value::Null => JsonShape::Null,
        Value::Bool(_) => JsonShape::Boolean,
        Value::Number(_) => JsonShape::Number,
        Value::String(_) => JsonShape::String,
        Value::Array(values) => JsonShape::Array(
            values
                .iter()
                .map(normalized_shape)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect(),
        ),
        Value::Object(values) => {
            let mut members = values
                .iter()
                .map(|(name, value)| (name.clone(), normalized_shape(value)))
                .collect::<Vec<_>>();
            members.sort_by(|left, right| left.0.cmp(&right.0));
            JsonShape::Object(members)
        }
    }
}

fn assert_same_shape(actual: &Value, expected: &Value, label: &str) {
    assert_eq!(
        normalized_shape(actual),
        normalized_shape(expected),
        "{label} schema shape changed"
    );
}

fn assert_close(actual: f64, expected: f64, path: &str) {
    let difference = (actual - expected).abs();
    let limit = ABS_TOLERANCE + REL_TOLERANCE * actual.abs().max(expected.abs());
    assert!(
        difference <= limit,
        "{path}: actual {actual:.17e}, expected {expected:.17e}, difference {difference:.17e} exceeds {limit:.17e}"
    );
}

fn assert_fixture_numbers_close(actual: &Value, expected: &Value, path: &str) {
    match (actual, expected) {
        (Value::Number(actual), Value::Number(expected)) => assert_close(
            actual.as_f64().expect("actual JSON number must fit f64"),
            expected.as_f64().expect("fixture JSON number must fit f64"),
            path,
        ),
        (Value::Array(actual), Value::Array(expected)) => {
            assert_eq!(
                actual.len(),
                expected.len(),
                "{path}: fixture array length changed"
            );
            for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
                assert_fixture_numbers_close(actual, expected, &format!("{path}[{index}]"));
            }
        }
        (Value::Object(actual), Value::Object(expected)) => {
            for (name, expected) in expected {
                let actual = actual
                    .get(name)
                    .unwrap_or_else(|| panic!("{path}: missing fixture member {name}"));
                assert_fixture_numbers_close(actual, expected, &format!("{path}.{name}"));
            }
        }
        _ => {}
    }
}

/// Compare stable non-numeric fixture values while leaving descriptive messages and the build's
/// engine version free to evolve. Shape comparison separately protects nullability and types.
fn assert_fixture_semantics(actual: &Value, expected: &Value, path: &str) {
    match (actual, expected) {
        (Value::String(actual), Value::String(expected))
            if path != "$.engine_version" && !path.ends_with(".message") =>
        {
            assert_eq!(actual, expected, "{path}: stable string value changed");
        }
        (Value::Bool(actual), Value::Bool(expected)) => {
            assert_eq!(actual, expected, "{path}: boolean value changed");
        }
        (Value::Null, Value::Null) | (Value::Number(_), Value::Number(_)) => {}
        (Value::Array(actual), Value::Array(expected)) => {
            assert_eq!(
                actual.len(),
                expected.len(),
                "{path}: fixture array length changed"
            );
            for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
                assert_fixture_semantics(actual, expected, &format!("{path}[{index}]"));
            }
        }
        (Value::Object(actual), Value::Object(expected)) => {
            for (name, expected) in expected {
                let actual = actual
                    .get(name)
                    .unwrap_or_else(|| panic!("{path}: missing fixture member {name}"));
                assert_fixture_semantics(actual, expected, &format!("{path}.{name}"));
            }
        }
        (Value::String(_), Value::String(_)) => {}
        _ => panic!("{path}: actual and fixture value kinds differ"),
    }
}

fn decode_request(document: &str) -> SolveRequestV1 {
    decode_solve_request_v1(document).expect("fixture request must decode")
}

fn direct_representative_solver() -> TrajectorySolver {
    let inputs = BallisticInputs {
        bc_value: 0.243,
        bc_type: DragModel::G7,
        bullet_mass: 0.01134,
        muzzle_velocity: 823.0,
        bullet_diameter: 0.00782,
        bullet_length: 0.031,
        muzzle_angle: 0.001,
        target_distance: 50.0,
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
        use_rk4: true,
        use_adaptive_rk45: true,
        ..BallisticInputs::default()
    };
    let mut solver = TrajectorySolver::new(
        inputs,
        WindConditions {
            speed: 4.0,
            direction: std::f64::consts::FRAC_PI_2,
            vertical_speed: 0.0,
        },
        AtmosphericConditions {
            temperature: 20.0,
            pressure: 1000.0,
            humidity: 40.0,
            altitude: 0.0,
        },
    );
    solver.set_max_range(50.0);
    solver.set_time_step(0.001);
    solver
}

#[test]
fn checked_fixtures_have_order_independent_v1_shapes() {
    let raw_request = parse_fixture(REPRESENTATIVE_REQUEST);
    let decoded_request = serde_json::to_value(decode_request(REPRESENTATIVE_REQUEST))
        .expect("request serialization");
    assert_same_shape(&decoded_request, &raw_request, "representative request");

    for (label, request, response) in [
        (
            "representative success",
            REPRESENTATIVE_REQUEST,
            REPRESENTATIVE_SUCCESS,
        ),
        ("validation error", VALIDATION_REQUEST, VALIDATION_ERROR),
        ("schema error", SCHEMA_REQUEST, SCHEMA_ERROR),
        ("resource error", RESOURCE_REQUEST, RESOURCE_ERROR),
        ("early termination", EARLY_REQUEST, EARLY_SUCCESS),
    ] {
        let actual = evaluate_request(request);
        let expected = parse_fixture(response);
        assert_same_shape(&actual, &expected, label);
        assert_fixture_semantics(&actual, &expected, "$");
    }
}

#[test]
fn representative_success_matches_tolerant_numeric_fixture() {
    let actual = evaluate_request(REPRESENTATIVE_REQUEST);
    let expected = parse_fixture(REPRESENTATIVE_SUCCESS);
    assert_same_shape(&actual, &expected, "representative success");
    assert_fixture_numbers_close(&actual, &expected, "$");

    assert_eq!(actual["schema_version"], 1);
    assert_eq!(actual["status"], "ok");
    assert_eq!(actual["summary"]["termination"], "max_range");
    assert_eq!(actual["assumptions"][0]["code"], "default_applied");
    assert_eq!(
        actual["assumptions"][0]["path"],
        "$.wind.vertical_speed_mps"
    );
    assert_eq!(actual["warnings"][0]["code"], "rk45_time_step_ignored");
    assert_eq!(actual["warnings"][0]["path"], "$.solver.time_step_s");
}

#[test]
fn representative_fixture_matches_direct_trajectory_solver() {
    let success = solve_v1(decode_request(REPRESENTATIVE_REQUEST)).expect("service solve");
    let direct = direct_representative_solver()
        .solve()
        .expect("direct trajectory solve");
    let observations = direct
        .sample_observations(10.0, 10_000)
        .expect("direct checked observations");

    assert_eq!(success.summary.termination, TerminationReasonV1::MaxRange);
    assert_eq!(success.samples.len(), observations.len());
    for (index, (actual, expected)) in success.samples.iter().zip(observations).enumerate() {
        let base = format!("$.samples[{index}]");
        assert_close(
            actual.distance_m,
            expected.distance_m,
            &format!("{base}.distance_m"),
        );
        assert_close(actual.time_s, expected.time_s, &format!("{base}.time_s"));
        assert_close(
            actual.speed_mps,
            expected.speed_mps,
            &format!("{base}.speed_mps"),
        );
        assert_close(
            actual.energy_j,
            expected.energy_j,
            &format!("{base}.energy_j"),
        );
        assert_close(actual.drop_m, expected.drop_m, &format!("{base}.drop_m"));
        assert_close(
            actual.windage_m,
            expected.windage_m,
            &format!("{base}.windage_m"),
        );
        assert_close(actual.mach, expected.mach, &format!("{base}.mach"));
        assert_eq!(
            actual.flags.contains(&SampleFlagV1::Terminal),
            expected
                .flags
                .contains(&ballistics_engine::TrajectoryObservationFlag::Terminal),
            "{base}.flags terminal marker differs"
        );
    }
}

#[test]
fn validation_schema_and_resource_error_fixtures_lock_code_and_path() {
    for (label, request, fixture, code, path) in [
        (
            "validation",
            VALIDATION_REQUEST,
            VALIDATION_ERROR,
            "invalid_value",
            "$.projectile.mass_kg",
        ),
        (
            "schema",
            SCHEMA_REQUEST,
            SCHEMA_ERROR,
            "unsupported_schema_version",
            "$.schema_version",
        ),
        (
            "resource",
            RESOURCE_REQUEST,
            RESOURCE_ERROR,
            "resource_limit",
            "$.sampling.interval_m",
        ),
    ] {
        let actual = evaluate_request(request);
        let expected = parse_fixture(fixture);
        assert_same_shape(&actual, &expected, label);
        assert_eq!(actual["schema_version"], 1, "{label}");
        assert_eq!(actual["status"], "error", "{label}");
        assert_eq!(actual["error"]["code"], code, "{label}");
        assert_eq!(actual["error"]["path"], path, "{label}");
        assert!(actual["error"]["line"].is_null(), "{label}");
        assert!(actual["error"]["column"].is_null(), "{label}");
    }

    let resource = solve_v1(decode_request(RESOURCE_REQUEST)).expect_err("resource limit");
    assert_eq!(resource.error.code, SolveErrorCodeV1::ResourceLimit);
    assert_eq!(resource.error.path(), Some("$.sampling.interval_m"));
}

#[test]
fn early_termination_fixture_has_one_exact_terminal_sample() {
    let success = solve_v1(decode_request(EARLY_REQUEST)).expect("early-ground solve");
    let actual = serde_json::to_value(&success).expect("early success serialization");
    let expected = parse_fixture(EARLY_SUCCESS);
    assert_same_shape(&actual, &expected, "early termination");
    assert_fixture_numbers_close(&actual, &expected, "$");

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
    assert!(terminal.flags.contains(&SampleFlagV1::Terminal));
    assert!(terminal.flags.contains(&SampleFlagV1::GroundThreshold));
    assert_close(
        terminal.distance_m,
        success.summary.actual_range_m,
        "$.samples[-1].distance_m",
    );
}
